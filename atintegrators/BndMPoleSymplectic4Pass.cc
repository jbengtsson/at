#include "elem.cc"
#include "tracy-2.cc"


struct elem_type*
trackFunction(const PyObject *ElemData, struct elem_type *Elem, double *ps,
	      const int num_particles, const struct parameters *Param)
{
  if (!Elem) Elem = init_mpole(ElemData, Elem, true);
  if (Elem) {
    MpolePass(ps, num_particles, Elem);
    return Elem;
  } else
    return NULL;
}
