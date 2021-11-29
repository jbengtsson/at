#include "elem.cc"
#include "tracy-2.cc"


/* Straight dipole w/ multipole using Symplectic Integration and rotation at
 * dipole faces.
 * Created by Xiaobiao Huang, 7/31/2018 */


struct elem_type*
trackFunction(const PyObject *ElemData, struct elem_type *Elem, double ps[],
	      const int num_particles, const struct parameters *Param)
{
  if (!Elem) Elem = init_mpole(ElemData, Elem, true, true, false, false);
  if (Elem) {
    CBendPass(ps, num_particles, Elem);
    return Elem;
  } else
    return NULL;
}
