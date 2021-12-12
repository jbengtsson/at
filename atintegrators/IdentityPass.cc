#include "tracy-2.cc"


struct elem_type*
trackFunction(const PyObject *ElemData, struct elem_type *Elem, double ps[],
	       const int num_particles, const struct parameters *Param)
{
  if (!Elem) Elem = init_id(ElemData, Elem);
  if (Elem) {
    IdentityPass(ps, num_particles, Elem);
    printf("\nHej pa dig!\n");
    return Elem;
  } else
    return NULL;
}
