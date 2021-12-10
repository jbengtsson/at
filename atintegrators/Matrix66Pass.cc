#include "tracy-2.cc"

struct elem_type*
trackFunction(const PyObject *ElemData,	struct elem_type *Elem, double ps[],
	      const int num_particles, const struct parameters *Param)
{
  if (!Elem) Elem = init_M66(ElemData, Elem);
  if (Elem) {
    Matrix66Pass(ps, num_particles, Elem);
    return Elem;
  } else
    return NULL;
}
