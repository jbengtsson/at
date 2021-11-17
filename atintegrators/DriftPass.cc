#include "atelem.cc"
#include "atlalib.cc"
#include "atphyslib.cc"


struct elem {
  double
  Length,
    *R1,
    *R2,
    *T1,
    *T2,
    *EApertures,
    *RApertures;
};


void DriftPass(double *r_in, double le, const double *T1, const double *T2,
	       const double *R1, const double *R2, double *RApertures,
	       double *EApertures, int num_particles)
/* le - physical length
   r_in - 6-by-N matrix of initial conditions reshaped into 
   1-d array of 6*N elements                                                  */
{
  double *r6;
  int    c;
  
#pragma omp parallel for if (num_particles > OMP_PARTICLE_THRESHOLD*10)    \
  default(shared) shared(r_in,num_particles) private(c,r6)
  for (c = 0; c<num_particles; c++) { /*Loop over particles  */
    r6 = r_in+c*6;
    if(!atIsNaN(r6[0])) {
      /*  misalignment at entrance  */
      if (T1) ATaddvv(r6, T1);
      if (R1) ATmultmv(r6, R1);
      /* Check physical apertures at the entrance of the magnet */
      if (RApertures) checkiflostRectangularAp(r6,RApertures);
      if (EApertures) checkiflostEllipticalAp(r6,EApertures);
      atdrift(r6, le);
      /* Check physical apertures at the exit of the magnet */
      if (RApertures) checkiflostRectangularAp(r6,RApertures);
      if (EApertures) checkiflostEllipticalAp(r6,EApertures);
      /* Misalignment at exit */
      if (R2) ATmultmv(r6, R2);
      if (T2) ATaddvv(r6, T2);
    }
  }
}


struct elem*
trackFunction(const atElem *ElemData,struct elem *Elem, double ps_in[],
	      int num_particles, struct parameters *Param)
{
  /*  if (ElemData) {*/
  if (!Elem) {
    double Length;
    double *R1, *R2, *T1, *T2, *EApertures, *RApertures;
    Length=atGetDouble(ElemData,"Length"); check_error();
    R1=atGetOptionalDoubleArray(ElemData, (char*)"R1"); check_error();
    R2=atGetOptionalDoubleArray(ElemData, (char*)"R2"); check_error();
    T1=atGetOptionalDoubleArray(ElemData, (char*)"T1"); check_error();
    T2=atGetOptionalDoubleArray(ElemData, (char*)"T2"); check_error();
    EApertures=atGetOptionalDoubleArray(ElemData, (char*)"EApertures");
    check_error();
    RApertures=atGetOptionalDoubleArray(ElemData, (char*)"RApertures");
    check_error();
    Elem = (struct elem*)atMalloc(sizeof(struct elem));
    Elem->Length=Length;
    Elem->R1=R1;
    Elem->R2=R2;
    Elem->T1=T1;
    Elem->T2=T2;
    Elem->EApertures=EApertures;
    Elem->RApertures=RApertures;
  }
  DriftPass(ps_in, Elem->Length, Elem->T1, Elem->T2, Elem->R1, Elem->R2,
	    Elem->RApertures, Elem->EApertures, num_particles);
  /*  }
      else {
      atFree(Elem->T1);
      atFree(Elem->T2);
      atFree(Elem->R1);
      atFree(Elem->R2);
      atFree(Elem->EApertures);
      atFree(Elem->RApertures);
      }*/
  return Elem;
}
