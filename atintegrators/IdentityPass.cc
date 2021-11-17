#include "atelem.cc"
#include "atlalib.cc"
#include "atphyslib.cc"


void IdentityPass(double *r_in, const double *T1, const double *T2,
		  const double *R1, const double *R2, const double *limits,
		  const double *axesptr, int num_particles)
{	
  double *r6;
  int    c;
    
  for (c = 0; c<num_particles; c++) {	/*Loop over particles  */
    r6 = r_in+c*6;
    if (!atIsNaN(r6[0])) {
      /*  misalignment at entrance  */
      if (T1) ATaddvv(r6, T1);
      if (R1) ATmultmv(r6, R1);
      /* Check physical apertures */
      if (limits) checkiflostRectangularAp(r6,limits);
      if (axesptr) checkiflostEllipticalAp(r6,axesptr);
      /* Misalignment at exit */
      if (R2) ATmultmv(r6, R2);
      if (T2) ATaddvv(r6, T2);
    }
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
