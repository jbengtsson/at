#include "elem.cc"
#include "tracy-2.cc"


struct elem_type* init_cav(const atElem *ElemData, struct elem_type *Elem)
{
  double   Length, Voltage, Energy, Frequency, TimeLag;
  elem_cav *cav;

  Elem          = (struct elem_type*)malloc(sizeof(struct elem_type));
  Elem->cav_ptr = (struct elem_cav*)malloc(sizeof(struct elem_cav));
  cav           = Elem->cav_ptr;

  Length    = atGetDouble(ElemData,"Length"); check_error();
  Voltage   = atGetDouble(ElemData,"Voltage"); check_error();
  Energy    = atGetDouble(ElemData,"Energy"); check_error();
  Frequency = atGetDouble(ElemData,"Frequency"); check_error();
  TimeLag   = atGetOptionalDouble(ElemData,"TimeLag",0); check_error();

  Elem->Length   = Length;
  cav->Voltage   = Voltage;
  cav->Energy    = Energy;
  cav->Frequency = Frequency;
  cav->TimeLag   = TimeLag;

  return Elem;
}

void CavityPass(double ps[], const int num_particles,
		const struct elem_type *Elem)
{
  int    k;
  double *ps_vec;

  const elem_cav *cav = Elem->cav_ptr;

  for(k = 0; k < num_particles; k++) {
    ps_vec = ps+k*6;
    if(!atIsNaN(ps_vec[0]))
      cav_pass(ps_vec,  Elem->Length, cav->Voltage, cav->Frequency,
	       cav->TimeLag);
  }
}

struct elem_type*
trackFunction(const atElem *ElemData, struct elem_type *Elem, double *ps,
	      const int num_particles, const struct parameters *Param)
{
  if (!Elem) Elem = init_cav(ElemData, Elem);
  CavityPass(ps, num_particles, Elem);
  return Elem;
}
