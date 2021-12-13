#include "tracy-2.cc"


extern "C" struct elem_type*
trackFunction(const PyObject *ElemData, struct elem_type *Elem, double ps[],
	      const int num_particles, const struct parameters *Param)
{
  if (!Elem) Elem = init_mpole(ElemData, Elem, false, false, false, false);
  if (Elem) {
    MpolePass(ps, num_particles, Elem, false);
    return Elem;
  } else
    return NULL;
}
