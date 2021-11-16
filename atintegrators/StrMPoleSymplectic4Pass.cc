#include "atelem.cc"
#include "atlalib.cc"
#include "atphyslib.cc"
#include "quadfringe.cc"	/* QuadFringePassP, QuadFringePassN */

#define DRIFT1  0.6756035959798286638
#define DRIFT2 -0.1756035959798286639
#define KICK1   1.351207191959657328
#define KICK2  -1.702414383919314656

struct elem
{
  double
  Length,
    *PolynomA,
    *PolynomB;
  int
  MaxOrder,
    NumIntSteps,
  /* Optional fields */
    FringeQuadEntrance,
    FringeQuadExit;
  double
  *fringeIntM0,
    *fringeIntP0,
    *R1,
    *R2,
    *T1,
    *T2,
    *RApertures,
    *EApertures,
    *KickAngle;
};


void StrMPoleSymplectic4Pass
(double *ps_n, double le, double *A, double *B, int max_order,
 int num_int_steps, int FringeQuadEntrance, int FringeQuadExit,
 /* 0 (no fringe), 1 (lee-whiting) or 2 (lee-whiting+elegant-like) */
 double *fringeIntM0,
 /* I0m/K1, I1m/K1, I2m/K1, I3m/K1, Lambda2m/K1 */
 double *fringeIntP0,
 /* I0p/K1, I1p/K1, I2p/K1, I3p/K1, Lambda2p/K1 */        
 double *T1, double *T2, double *R1, double *R2, double *RApertures,
 double *EApertures, double *KickAngle, int num_particles)
{
  int
    c;
  double
    SL = le/num_int_steps,
    L1 = SL*DRIFT1,
    L2 = SL*DRIFT2,
    K1 = SL*KICK1,
    K2 = SL*KICK2;
  bool
    useLinFrEleEntrance =
    (fringeIntM0 != NULL && fringeIntP0 != NULL && FringeQuadEntrance==2),
    useLinFrEleExit     =
    (fringeIntM0 != NULL && fringeIntP0 != NULL && FringeQuadExit==2);

  if (KickAngle) { // Convert corrector component to polynomial coefficients
    B[0] -= sin(KickAngle[0])/le; 
    A[0] += sin(KickAngle[1])/le;
  }

#pragma omp parallel for if (num_particles > OMP_PARTICLE_THRESHOLD)        \
  default(none)								    \
  shared(ps_n, num_particles, R1, T1, R2, T2, RApertures, EApertures, A, B, \
	 L1, L2, K1, K2, max_order, num_int_steps, FringeQuadEntrance,      \
	 useLinFrEleEntrance, FringeQuadExit, useLinFrEleExit, fringeIntM0, \
	 fringeIntP0)							    \
  private(c)

  for (c = 0; c < num_particles; c++) { /*Loop over particles  */
    double *r6 = ps_n+c*6;
    if(!atIsNaN(r6[0])) {
      int
	m;
      double
	norm = 1.0/(1.0+r6[4]),
	NormL1 = L1*norm,
	NormL2 = L2*norm;

      /*  misalignment at entrance  */
      if (T1) ATaddvv(r6, T1);
      if (R1) ATmultmv(r6, R1);
#if 0
      /* Check physical apertures at the entrance of the magnet */
      if (RApertures) checkiflostRectangularAp(r6,RApertures);
      if (EApertures) checkiflostEllipticalAp(r6,EApertures);
#endif
      if (FringeQuadEntrance && B[1] != 0) {
	if (useLinFrEleEntrance) /*Linear fringe fields from elegant */
	  linearQuadFringeElegantEntrance(r6, B[1], fringeIntM0, fringeIntP0);
	else
	  QuadFringePassP(r6, B[1]);
      }
      /*  integrator  */
      for (m=0; m < num_int_steps; m++) {  /*  Loop over slices */
	fastdrift(r6, NormL1);
	strthinkick(r6, A, B, K1, max_order);
	fastdrift(r6, NormL2);
	strthinkick(r6, A, B, K2, max_order);
	fastdrift(r6, NormL2);
	strthinkick(r6, A, B, K1, max_order);
	fastdrift(r6, NormL1);
      }
      if (FringeQuadExit && B[1]!=0) {
	if (useLinFrEleExit) /*Linear fringe fields from elegant*/
	  linearQuadFringeElegantExit(r6, B[1], fringeIntM0, fringeIntP0);
	else
	  QuadFringePassN(r6, B[1]);
      }
      /* Check physical apertures at the exit of the magnet */
      if (RApertures) checkiflostRectangularAp(r6,RApertures);
      if (EApertures) checkiflostEllipticalAp(r6,EApertures);
      /* Misalignment at exit */
      if (R2) ATmultmv(r6, R2);
      if (T2) ATaddvv(r6, T2);
    }
  }
  if (KickAngle) { /* Remove corrector component in polynomial coefficients */
    B[0] += sin(KickAngle[0])/le; 
    A[0] -= sin(KickAngle[1])/le;
  }
}


#if defined(MATLAB_MEX_FILE) || defined(PYAT)
ExportMode struct elem*
trackFunction(const atElem *ElemData,struct elem *Elem, double *r_in,
	      int num_particles, struct parameters *Param)
{
  if (!Elem) {
    int
      MaxOrder,
      NumIntSteps,
      FringeQuadEntrance,
      FringeQuadExit;
    double
      Length,
      *PolynomA,
      *PolynomB,
      *R1,
      *R2,
      *T1,
      *T2,
      *EApertures,
      *RApertures,
      *fringeIntM0,
      *fringeIntP0,
      *KickAngle;

    Length = atGetDouble(ElemData, (char*)"Length"); check_error();
    PolynomA = atGetDoubleArray(ElemData, (char*)"PolynomA"); check_error();
    PolynomB = atGetDoubleArray(ElemData, (char*)"PolynomB"); check_error();
    MaxOrder = atGetLong(ElemData, (char*)"MaxOrder"); check_error();
    NumIntSteps = atGetLong(ElemData, (char*)"NumIntSteps"); check_error();
    /*optional fields*/
    FringeQuadEntrance =
      atGetOptionalLong(ElemData, (char*)"FringeQuadEntrance", 0);
    FringeQuadExit = atGetOptionalLong(ElemData, (char*)"FringeQuadExit", 0);
    fringeIntM0 = atGetOptionalDoubleArray(ElemData, (char*)"fringeIntM0");
    check_error();
    fringeIntP0 = atGetOptionalDoubleArray(ElemData, (char*)"fringeIntP0");
    check_error();
    R1 = atGetOptionalDoubleArray(ElemData, (char*)"R1"); check_error();
    R2 = atGetOptionalDoubleArray(ElemData, (char*)"R2"); check_error();
    T1 = atGetOptionalDoubleArray(ElemData, (char*)"T1"); check_error();
    T2 = atGetOptionalDoubleArray(ElemData, (char*)"T2"); check_error();
    EApertures = atGetOptionalDoubleArray(ElemData, (char*)"EApertures");
    check_error();
    RApertures = atGetOptionalDoubleArray(ElemData, (char*)"RApertures");
    check_error();
    KickAngle = atGetOptionalDoubleArray(ElemData, (char*)"KickAngle");
    check_error();
        
    Elem  =  (struct elem*)atMalloc(sizeof(struct elem));
    Elem->Length = Length;
    Elem->PolynomA = PolynomA;
    Elem->PolynomB = PolynomB;
    Elem->MaxOrder = MaxOrder;
    Elem->NumIntSteps = NumIntSteps;
    /*optional fields*/
    Elem->FringeQuadEntrance = FringeQuadEntrance;
    Elem->FringeQuadExit = FringeQuadExit;
    Elem->fringeIntM0 = fringeIntM0;
    Elem->fringeIntP0 = fringeIntP0;
    Elem->R1 = R1;
    Elem->R2 = R2;
    Elem->T1 = T1;
    Elem->T2 = T2;
    Elem->EApertures = EApertures;
    Elem->RApertures = RApertures;
    Elem->KickAngle = KickAngle;
  }
  StrMPoleSymplectic4Pass
    (r_in, Elem->Length, Elem->PolynomA, Elem->PolynomB, 
     Elem->MaxOrder, Elem->NumIntSteps, Elem->FringeQuadEntrance, 
     Elem->FringeQuadExit, Elem->fringeIntM0, Elem->fringeIntP0, 
     Elem->T1, Elem->T2, Elem->R1, Elem->R2, 
     Elem->RApertures, Elem->EApertures, Elem->KickAngle, num_particles);
  return Elem;
}

MODULE_DEF(StrMPoleSymplectic4Pass) /* Dummy module initialisation */

#endif /*defined(MATLAB_MEX_FILE) || defined(PYAT)*/
