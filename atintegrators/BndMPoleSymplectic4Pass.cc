#include "elem.cc"
#include "atlalib.cc"

#define DRIFT1    0.6756035959798286638
#define DRIFT2   -0.1756035959798286639
#define KICK1     1.351207191959657328
#define KICK2    -1.702414383919314656

#define SQR(X) ((X)*(X))


struct elem_type* init_mpole(const atElem *ElemData, struct elem_type *Elem)
{
  int
    MaxOrder, NumIntSteps,  FringeBendEntrance, FringeBendExit,
    FringeQuadEntrance, FringeQuadExit;
  double
    Length, BendingAngle, EntranceAngle, ExitAngle, FullGap, FringeInt1,
    FringeInt2, *PolynomA, *PolynomB, *R1, *R2, *T1, *T2, *EApertures,
    *RApertures, *fringeIntM0, *fringeIntP0, *KickAngle;
  elem_mpole *mpole;

  Elem = (struct elem_type*)malloc(sizeof(struct elem_type));
  Elem->mpole_ptr = (struct elem_mpole*)malloc(sizeof(struct elem_mpole));
  mpole            = Elem->mpole_ptr;

  Length        = atGetDouble(ElemData, "Length"); check_error();
  PolynomA      = atGetDoubleArray(ElemData, (char*)"PolynomA");
  check_error();
  PolynomB      = atGetDoubleArray(ElemData, (char*)"PolynomB");
  check_error();
  MaxOrder      = atGetLong(ElemData, (char*)"MaxOrder");
  check_error();
  NumIntSteps   = atGetLong(ElemData, (char*)"NumIntSteps");
  check_error();

  //----------------------------------------------------------------------------
  BendingAngle  = atGetDouble(ElemData, (char*)"BendingAngle");
  check_error();
  EntranceAngle = atGetDouble(ElemData, (char*)"EntranceAngle");
  check_error();
  ExitAngle = atGetDouble(ElemData, (char*)(char*)"ExitAngle");
  check_error();
  FringeBendEntrance =
    atGetOptionalLong(ElemData, (char*)"FringeBendEntrance", 1);
  check_error();
  FringeBendExit     = atGetOptionalLong(ElemData, (char*)"FringeBendExit", 1);
  check_error();
  FullGap            = atGetOptionalDouble(ElemData, (char*)"FullGap", 0);
  check_error();
  FringeInt1         = atGetOptionalDouble(ElemData, (char*)"FringeInt1", 0);
  check_error();
  FringeInt2         = atGetOptionalDouble(ElemData, (char*)"FringeInt2", 0);
  check_error();
  //----------------------------------------------------------------------------

  FringeQuadEntrance =
    atGetOptionalLong(ElemData, (char*)"FringeQuadEntrance", 0);
  check_error();
  FringeQuadExit     = atGetOptionalLong(ElemData, (char*)"FringeQuadExit", 0);
  check_error();
  fringeIntM0        = atGetOptionalDoubleArray(ElemData, (char*)"fringeIntM0");
  check_error();
  fringeIntP0        = atGetOptionalDoubleArray(ElemData, (char*)"fringeIntP0");
  check_error();
  
  R1 = atGetOptionalDoubleArray(ElemData, (char*)"R1");
  check_error();
  R2 = atGetOptionalDoubleArray(ElemData, (char*)"R2");
  check_error();
  T1 = atGetOptionalDoubleArray(ElemData, (char*)"T1");
  check_error();
  T2 = atGetOptionalDoubleArray(ElemData, (char*)"T2");
  check_error();

  EApertures = atGetOptionalDoubleArray(ElemData, (char*)"EApertures");
  check_error();
  RApertures = atGetOptionalDoubleArray(ElemData, (char*)"RApertures");
  check_error();
  KickAngle = atGetOptionalDoubleArray(ElemData, (char*)"KickAngle");
  check_error();
        
  Elem->Length             = Length;
  Elem->R1                 = R1;
  Elem->R2                 = R2;
  Elem->T1                 = T1;
  Elem->T2                 = T2;
  Elem->EApertures         = EApertures;
  Elem->RApertures         = RApertures;

  mpole->PolynomA          = PolynomA;
  mpole->PolynomB          = PolynomB;
  mpole->MaxOrder          = MaxOrder;
  mpole->NumIntSteps       = NumIntSteps;

  //----------------------------------------------------------------------------
  mpole->BendingAngle       = BendingAngle;
  mpole->EntranceAngle      = EntranceAngle;
  mpole->ExitAngle          = ExitAngle;
  mpole->FringeBendEntrance = FringeBendEntrance;
  mpole->FringeBendExit     = FringeBendExit;
  mpole->FullGap            = FullGap;
  mpole->FringeInt1         = FringeInt1;
  mpole->FringeInt2         = FringeInt2;
  //----------------------------------------------------------------------------

  mpole->FringeQuadEntrance = FringeQuadEntrance;
  mpole->FringeQuadExit     = FringeQuadExit;
  mpole->fringeIntM0        = fringeIntM0;
  mpole->fringeIntP0        = fringeIntP0;
  mpole->KickAngle          = KickAngle;

  return Elem;
}

void BendPass(double *ps, const int num_particles, const struct elem_type *Elem)
{	
  int    k;
  double *ps_vec;

  const elem_mpole *mpole = Elem->mpole_ptr;

  const bool
    useLinFrEleEntrance =
    (mpole->fringeIntM0 != NULL && mpole->fringeIntP0 != NULL
     && mpole->FringeQuadEntrance == 2),
    useLinFrEleExit =
    (mpole->fringeIntM0 != NULL && mpole->fringeIntP0 != NULL
     && mpole->FringeQuadExit == 2);

  const double
    SL = Elem->Length/mpole->NumIntSteps,
    L1 = SL*DRIFT1,
    L2 = SL*DRIFT2,
    K1 = SL*KICK1,
    K2 = SL*KICK2,

    irho = mpole->BendingAngle/Elem->Length;

  if (mpole->KickAngle) {
    /* Convert corrector component to polynomial coefficients */
    mpole->PolynomB[0] -= sin(mpole->KickAngle[0])/Elem->Length; 
    mpole->PolynomA[0] += sin(mpole->KickAngle[1])/Elem->Length;
  }
#pragma omp parallel for if (num_particles > OMP_PARTICLE_THRESHOLD)	\
  default(none)								\
  shared(r, num_particles, Elem)					\
  private(k)

  for(k = 0; k < num_particles; k++) {	/* Loop over particles  */
    ps_vec = ps+k*6;
    if (!atIsNaN(ps_vec[0])) {
      int    m;
      double norm = 1.0/(1.0+ps_vec[4]), NormL1 = L1*norm, NormL2 = L2*norm;

      /*  misalignment at entrance  */
      if (Elem->T1) ATaddvv(ps_vec, Elem->T1);
      if (Elem->R1) ATmultmv(ps_vec, Elem->R1);
      /* Check physical apertures at the entrance of the magnet */
      if (Elem->RApertures) checkiflostRectangularAp(ps_vec, Elem->RApertures);
      if (Elem->EApertures) checkiflostEllipticalAp(ps_vec, Elem->EApertures);

      /* edge focus */
      edge_fringe_entrance
	(ps_vec, irho, mpole->EntranceAngle, mpole->FringeInt1, mpole->FullGap,
	 mpole->FringeBendEntrance);

      /* quadrupole gradient fringe entrance*/
      if (mpole->FringeQuadEntrance && mpole->PolynomB[1]!=0) {
	if (useLinFrEleEntrance) /*Linear fringe fields from elegant*/
	  linearQuadFringeElegantEntrance
	    (ps_vec, mpole->PolynomB[1], mpole->fringeIntM0,
	     mpole->fringeIntP0);
	else
	  QuadFringePassP(ps_vec, mpole->PolynomB[1]);
      }    

      /* integrator */
      for(m=0; m < mpole->NumIntSteps; m++) { /* Loop over slices*/
	fastdrift(ps_vec, NormL1);
	bndthinkick(ps_vec, mpole->PolynomA, mpole->PolynomB, K1, irho,
		    mpole->MaxOrder);
	fastdrift(ps_vec, NormL2);
	bndthinkick(ps_vec, mpole->PolynomA, mpole->PolynomB, K2, irho,
		    mpole->MaxOrder);
	fastdrift(ps_vec, NormL2);
	bndthinkick(ps_vec, mpole->PolynomA, mpole->PolynomB, K1, irho,
		    mpole->MaxOrder);
	fastdrift(ps_vec, NormL1);
      }

      /* quadrupole gradient fringe */
      if (mpole->FringeQuadExit && mpole->PolynomB[1]!=0) {
	if (useLinFrEleExit) /*Linear fringe fields from elegant*/
	  linearQuadFringeElegantExit
	    (ps_vec, mpole->PolynomB[1], mpole->fringeIntM0,
	     mpole->fringeIntP0);
	else
	  QuadFringePassN(ps_vec, mpole->PolynomB[1]);
      }

      /* edge focus */
      edge_fringe_exit
	(ps_vec, irho, mpole->ExitAngle, mpole->FringeInt2, mpole->FullGap,
	 mpole->FringeBendExit);

      /* Check physical apertures at the exit of the magnet */
      if (Elem->RApertures) checkiflostRectangularAp(ps_vec, Elem->RApertures);
      if (Elem->EApertures) checkiflostEllipticalAp(ps_vec, Elem->EApertures);

      /* Misalignment at exit */
      if (Elem->R2) ATmultmv(ps_vec, Elem->R2);
      if (Elem->T2) ATaddvv(ps_vec, Elem->T2);
    }
  }
  if (mpole->KickAngle) {
    /* Remove corrector component in polynomial coefficients */
    mpole->PolynomB[0] += sin(mpole->KickAngle[0])/Elem->Length;
    mpole->PolynomA[0] -= sin(mpole->KickAngle[1])/Elem->Length;
  }
}

struct elem_type*
trackFunction(const atElem *ElemData, struct elem_type *Elem,
	      double *ps, int num_particles, struct parameters *Param)
{
  if (!Elem) Elem = init_mpole(ElemData, Elem);
  BendPass(ps, num_particles, Elem);
  return Elem;
}
