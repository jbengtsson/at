/* Exact integrator for different element types
 
   This method will work for a drift, a quadrupole, a sextupole or a
   bending magnet. It distinguishes between these using the Class field
   on the element.

   The 'ExactHamiltonianPass' method uses the square root hamiltonian in
   cartesian co-ordinates (see other notes for derivation). 
   This is equivalent to setting exact=true in MADX-PTC.
   Multipole fringe fields are also enabled for quadrupoles
   (fringe = true option in MADX-PTC). 

   Note that the PolynomB array in the exact cartesian rectangular bend 
   refers to the normalized straight multipole components of the vector 
   potential, so PolynomB(1) should be set to 1/rho (B_bend / Brho). 
   The conjugate momenta in the curvilinear co-ordinate system are not 
   the same as in the cartesian system so PolynomB(1) must be set back 
   to zero when using a curvilinear symplectic integrator method such 
   as the 'BndMPoleSymplectic4E2Pass'. See Forest p362 for a detailed 
   explanation of the vector potential in curvilinear co-ordinates.           */

#include <string.h>

#include "elem.cc"
#include "tracy-2.cc"


struct elem_type*
trackFunction(const PyObject *ElemData, struct elem_type *Elem, double ps[],
	      const int num_particles, const struct parameters *Param)
{
  if (!Elem) Elem = init_H(ElemData, Elem);
  if (Elem) {
    HamPass(ps, num_particles, Elem);
    return Elem;
  } else
    return NULL;
}
