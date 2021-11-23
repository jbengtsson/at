#include "elem.cc"
#include "atlalib.cc"


struct elem_type* init_id(const atElem *ElemData, struct elem_type *Elem)
{
  double *R1, *R2, *T1, *T2, *EApertures, *RApertures;

  Elem = (struct elem_type*)malloc(sizeof(struct elem_type));
  Elem->id_ptr = (struct elem_id*)malloc(sizeof(struct elem_id));

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

  Elem->Length     = 0e0;
  Elem->R1         = R1;
  Elem->R2         = R2;
  Elem->T1         = T1;
  Elem->T2         = T2;
  Elem->EApertures = EApertures;
  Elem->RApertures = RApertures;

  return Elem;
}

void IdentityPass(double ps[], const int num_particles,
		  const struct elem_type *Elem)
{
  int    k;
  double *ps_vec;
    
  for (k = 0; k < num_particles; k++) {	/*Loop over particles  */
    ps_vec = ps+k*6;
    if (!atIsNaN(ps_vec[0])) {
      /*  misalignment at entrance  */
      if (Elem->T1) ATaddvv(ps_vec, Elem->T1);
      if (Elem->R1) ATmultmv(ps_vec, Elem->R1);
      /* Check physical apertures */
      if (Elem->RApertures) checkiflostRectangularAp(ps_vec, Elem->RApertures);
      if (Elem->EApertures) checkiflostEllipticalAp(ps_vec, Elem->EApertures);
      /* Misalignment at exit */
      if (Elem->R2) ATmultmv(ps_vec, Elem->R2);
      if (Elem->T2) ATaddvv(ps_vec, Elem->T2);
    }
  }
}

struct elem_type*
trackFunction(const atElem *ElemData, struct elem_type *Elem,
	      double ps[], int num_particles, struct parameters *Param)
{
  if (!Elem) Elem = init_id(ElemData, Elem);
  IdentityPass(ps, num_particles, Elem);
  return Elem;
}
