#include "atelem.cc"
#include "atlalib.cc"


void IdentityPass(double ps_n[], const double t1[], const double t2[],
		  const double r1[], const double r2[], const double limits[],
		  const double *axesptr, const int num_particles)
{	
  int       j, k;
  arma::vec ps(PS_DIM);

  for (j = 0; j < num_particles; j++) {	// Loop over particles.
    for (k = 0; k < PS_DIM; k++)
      ps(k) = ps_n[j*PS_DIM+k];
    if (!isnan(ps(0))) {
      //  misalignment at entrance.
      if (t1) ps = ps + arrtovec(t1);
      if (r1) ps = arrtomat(r1)*ps;
      // Check physical apertures.
#if 0
      if (limits) checkiflostRectangularAp(r6,limits);
      if (axesptr) checkiflostEllipticalAp(r6,axesptr);
#endif
      // Misalignment at exit.
      if (r2) ps = arrtomat(r2)*ps;
      if (t2) ps = ps + arrtovec(t2);
    }
    for (k = 0; k < PS_DIM; k++)
      ps_n[j*PS_DIM+k] = ps(k);
  }
}


struct elem {
  double
    *R1,
    *R2,
    *T1,
    *T2,
    *EApertures,
    *RApertures;
};


struct elem*
trackFunction(const atElem *ElemData, struct elem *Elem, double ps_in[],
	      int num_particles, struct parameters *Param)
{
  if (!Elem) {
    double *R1, *R2, *T1, *T2, *EApertures, *RApertures;
    R1 = atGetOptionalDoubleArray(ElemData, (char*)"R1"); check_error();
    R2 = atGetOptionalDoubleArray(ElemData, (char*)"R2"); check_error();
    T1 = atGetOptionalDoubleArray(ElemData, (char*)"T1"); check_error();
    T2 = atGetOptionalDoubleArray(ElemData, (char*)"T2"); check_error();
    EApertures = atGetOptionalDoubleArray(ElemData, (char*)"EApertures");
    check_error();
    RApertures = atGetOptionalDoubleArray(ElemData, (char*)"RApertures");
    check_error();
    Elem = (struct elem*)atMalloc(sizeof(struct elem));
    Elem->R1 = R1;
    Elem->R2 = R2;
    Elem->T1 = T1;
    Elem->T2 = T2;
    Elem->EApertures = EApertures;
    Elem->RApertures = RApertures;
  }
  IdentityPass(ps_in, Elem->T1, Elem->T2, Elem->R1, Elem->R2, 
	       Elem->RApertures, Elem->EApertures, num_particles);
  return Elem;
}
