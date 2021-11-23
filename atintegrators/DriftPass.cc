#include "elem.cc"
#include "tracy-2.cc"


struct elem_type* init_drift(const atElem *ElemData, struct elem_type *Elem)
{
  double Length, *R1, *R2, *T1, *T2, *EApertures, *RApertures;

  Elem            = (struct elem_type*)malloc(sizeof(struct elem_type));
  Elem->drift_ptr = (struct elem_drift*)malloc(sizeof(struct elem_drift));

  Length     = atGetDouble(ElemData, "Length");
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

  Elem->Length     = Length;
  Elem->R1         = R1;
  Elem->R2         = R2;
  Elem->T1         = T1;
  Elem->T2         = T2;
  Elem->EApertures = EApertures;
  Elem->RApertures = RApertures;

  return Elem;
}

void DriftPass(double *ps, const int num_particles,
	       const struct elem_type *Elem)
/* le - physical length
   ps - 6-by-N matrix of initial conditions reshaped into 
   1-d array of 6*N elements                                                  */
{
  int    k;
  double *ps_vec;
  
#pragma omp parallel for if (num_particles > OMP_PARTICLE_THRESHOLD*10)    \
  default(shared) shared(ps, num_particles) private(k, ps_vec)

  for (k = 0; k < num_particles; k++) { /*Loop over particles  */
    ps_vec = ps+k*6;
    if(!atIsNaN(ps_vec[0])) {
      /*  misalignment at entrance  */
      if (Elem->T1) ATaddvv(ps_vec, Elem->T1);
      if (Elem->R1) ATmultmv(ps_vec, Elem->R1);
      /* Check physical apertures at the entrance of the magnet */
      if (Elem->RApertures) checkiflostRectangularAp(ps_vec, Elem->RApertures);
      if (Elem->EApertures) checkiflostEllipticalAp(ps_vec, Elem->EApertures);
      atdrift(ps_vec, Elem->Length);
      /* Check physical apertures at the exit of the magnet */
      if (Elem->RApertures) checkiflostRectangularAp(ps_vec, Elem->RApertures);
      if (Elem->EApertures) checkiflostEllipticalAp(ps_vec, Elem->EApertures);
      /* Misalignment at exit */
      if (Elem->R2) ATmultmv(ps_vec, Elem->R2);
      if (Elem->T2) ATaddvv(ps_vec, Elem->T2);
    }
  }
}

struct elem_type*
trackFunction(const atElem *ElemData, struct elem_type *Elem, double ps[],
	      const int num_particles, const struct parameters *Param)
{
  if (!Elem) Elem = init_drift(ElemData, Elem);
  DriftPass(ps, num_particles, Elem);
  return Elem;
}
