#include "atelem.cc"
#include "atlalib.cc"
#include "atphyslib.cc"

#define DRIFT1    0.6756035959798286638
#define DRIFT2   -0.1756035959798286639
#define KICK1     1.351207191959657328
#define KICK2    -1.702414383919314656


#define SQR(X) ((X)*(X))


void BendPass
(double *r, double le, double irho, double *A, double *B, int max_order,
 int num_int_steps, double entrance_angle, double exit_angle,
 int FringeBendEntrance, int FringeBendExit, double fint1, double fint2,
 double gap, int FringeQuadEntrance, int FringeQuadExit, double *fringeIntM0,
 /* I0m/K1, I1m/K1, I2m/K1, I3m/K1, Lambda2m/K1 */
 double *fringeIntP0,
 /* I0p/K1, I1p/K1, I2p/K1, I3p/K1, Lambda2p/K1 */        
 double *T1, double *T2, double *R1, double *R2, double *RApertures,
 double *EApertures,  double *KickAngle, int num_particles)
{	
  int    c;
  double SL = le/num_int_steps;
  double L1 = SL*DRIFT1;
  double L2 = SL*DRIFT2;
  double K1 = SL*KICK1;
  double K2 = SL*KICK2;
  bool useLinFrEleEntrance =
    (fringeIntM0 != NULL && fringeIntP0 != NULL  && FringeQuadEntrance==2);
  bool useLinFrEleExit =
    (fringeIntM0 != NULL && fringeIntP0 != NULL  && FringeQuadExit==2);

  if (KickAngle) { /* Convert corrector component to polynomial coefficients */
    B[0] -= sin(KickAngle[0])/le; 
    A[0] += sin(KickAngle[1])/le;
  }
#pragma omp parallel for if (num_particles > OMP_PARTICLE_THRESHOLD)    \
  default(none)								\
  shared(r, num_particles, R1, T1, R2, T2, RApertures, EApertures,	\
	 irho, gap, A, B, L1, L2, K1, K2, max_order, num_int_steps,	\
	 FringeBendEntrance, entrance_angle, fint1, FringeBendExit,     \
	 exit_angle, fint2, FringeQuadEntrance, useLinFrEleEntrance,    \
	 FringeQuadExit, useLinFrEleExit, fringeIntM0, fringeIntP0)	\
  private(c)
  for(c = 0;c <num_particles;c++)	/* Loop over particles  */
    {
      double *r6 = r+c*6;
      if (!atIsNaN(r6[0])) {
	int m;
	double p_norm = 1.0/(1.0+r6[4]);
	double NormL1 = L1*p_norm;
	double NormL2 = L2*p_norm;
	/*  misalignment at entrance  */
	if (T1) ATaddvv(r6, T1);
	if (R1) ATmultmv(r6, R1);
	/* Check physical apertures at the entrance of the magnet */
	if (RApertures) checkiflostRectangularAp(r6, RApertures);
	if (EApertures) checkiflostEllipticalAp(r6, EApertures);
	/* edge focus */
	edge_fringe_entrance(r6, irho, entrance_angle, fint1, gap,
			     FringeBendEntrance);
	/* quadrupole gradient fringe entrance*/
	if (FringeQuadEntrance && B[1]!=0) {
	  if (useLinFrEleEntrance) /*Linear fringe fields from elegant*/
	    linearQuadFringeElegantEntrance(r6, B[1], fringeIntM0, fringeIntP0);
	  else
	    QuadFringePassP(r6, B[1]);
	}    
	/* integrator */
	for(m=0; m < num_int_steps; m++) /* Loop over slices*/
	  {
	    fastdrift(r6, NormL1);
	    bndthinkick(r6, A, B, K1, irho, max_order);
	    fastdrift(r6, NormL2);
	    bndthinkick(r6, A, B, K2, irho, max_order);
	    fastdrift(r6, NormL2);
	    bndthinkick(r6, A, B, K1, irho, max_order);
	    fastdrift(r6, NormL1);
	  }
	/* quadrupole gradient fringe */
	if (FringeQuadExit && B[1]!=0) {
	  if (useLinFrEleExit) /*Linear fringe fields from elegant*/
	    linearQuadFringeElegantExit(r6, B[1], fringeIntM0, fringeIntP0);
	  else
	    QuadFringePassN(r6, B[1]);
	}
	/* edge focus */
	edge_fringe_exit(r6, irho, exit_angle, fint2, gap, FringeBendExit);
	/* Check physical apertures at the exit of the magnet */
	if (RApertures) checkiflostRectangularAp(r6, RApertures);
	if (EApertures) checkiflostEllipticalAp(r6, EApertures);
	/* Misalignment at exit */
	if (R2) ATmultmv(r6, R2);
	if (T2) ATaddvv(r6, T2);
      }
    }
  if (KickAngle) {  /* Remove corrector component in polynomial coefficients */
    B[0] += sin(KickAngle[0])/le;
    A[0] -= sin(KickAngle[1])/le;
  }
}

struct elem_type* init_bend(const atElem *ElemData, struct elem_type *Elem)
{
  int
    MaxOrder, NumIntSteps,  FringeBendEntrance, FringeBendExit,
    FringeQuadEntrance, FringeQuadExit;
  double
    Length, BendingAngle, EntranceAngle, ExitAngle, FullGap, FringeInt1,
    FringeInt2, *PolynomA, *PolynomB, *R1, *R2, *T1, *T2, *EApertures,
    *RApertures, *fringeIntM0, *fringeIntP0, *KickAngle;
  elem_bend *bend;

  Elem = (struct elem_type*)malloc(sizeof(struct elem_type));
  Elem->bend_ptr = (struct elem_bend*)malloc(sizeof(struct elem_bend));
  bend               = Elem->bend_ptr;

  Length        = atGetDouble(ElemData, "Length"); check_error();
  PolynomA      = atGetDoubleArray(ElemData, (char*)"PolynomA");
  check_error();
  PolynomB      = atGetDoubleArray(ElemData, (char*)"PolynomB");
  check_error();
  MaxOrder      = atGetLong(ElemData, (char*)"MaxOrder");
  check_error();
  NumIntSteps   = atGetLong(ElemData, (char*)"NumIntSteps");
  check_error();
  BendingAngle  = atGetDouble(ElemData, (char*)"BendingAngle");
  check_error();
  EntranceAngle = atGetDouble(ElemData, (char*)"EntranceAngle");
  check_error();
  ExitAngle = atGetDouble(ElemData, (char*)(char*)"ExitAngle");
  check_error();
  /*optional fields*/
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

  bend->PolynomA           = PolynomA;
  bend->PolynomB           = PolynomB;
  bend->MaxOrder           = MaxOrder;
  bend->NumIntSteps        = NumIntSteps;
  bend->BendingAngle       = BendingAngle;
  bend->EntranceAngle      = EntranceAngle;
  bend->ExitAngle          = ExitAngle;
  bend->FringeBendEntrance = FringeBendEntrance;
  bend->FringeBendExit     = FringeBendExit;
  bend->FullGap            = FullGap;
  bend->FringeInt1         = FringeInt1;
  bend->FringeInt2         = FringeInt2;
  bend->FringeQuadEntrance = FringeQuadEntrance;
  bend->FringeQuadExit     = FringeQuadExit;
  bend->fringeIntM0        = fringeIntM0;
  bend->fringeIntP0        = fringeIntP0;
  bend->KickAngle          = KickAngle;

  return Elem;
}

struct elem_type*
trackFunction(const atElem *ElemData, struct elem_type *Elem,
	      double *r_in, int num_particles, struct parameters *Param)
{
  int       irho;
  elem_bend *bend;

  if (!Elem) Elem = init_bend(ElemData, Elem);
  bend = Elem->bend_ptr;
  irho = bend->BendingAngle/Elem->Length;

  BendPass
    (r_in, Elem->Length, irho, bend->PolynomA, bend->PolynomB, bend->MaxOrder,
     bend->NumIntSteps, bend->EntranceAngle, bend->ExitAngle,
     bend->FringeBendEntrance, bend->FringeBendExit, bend->FringeInt1,
     bend->FringeInt2, bend->FullGap, bend->FringeQuadEntrance,
     bend->FringeQuadExit, bend->fringeIntM0, bend->fringeIntP0, Elem->T1,
     Elem->T2, Elem->R1, Elem->R2, Elem->RApertures, Elem->EApertures,
     bend->KickAngle, num_particles);

  return Elem;
}
