#include "tracy-2.cc"
#include "element_integrators_macros.h"
#include "ElemPass.h"

/*
 * helpers for automatic initialisation
 *
 * these could be avoided by consistent naming
 */
static struct elem_type *init_aper(const PyObject *ElemData, struct elem_type *Elem)
{
  /*
   * why do I need to pass on a zero point:
   * Santa Claus coming to town?
   */
  return init_ap(ElemData, Elem);
}


/*
 * helpers for automatic initialisation
 *
 * Can these be avoided by consistent naming ?
 */
static struct elem_type * init_mpole_rad(const PyObject *ElemData, struct elem_type *Elem)
{
  return init_mpole(ElemData, Elem, true, false, true, false);
}

static struct elem_type * init_bend(const PyObject *ElemData, struct elem_type *Elem)
{
  return init_mpole(ElemData, Elem, true, false, false, false);
}

static struct elem_type * init_bend_rad(const PyObject *ElemData, struct elem_type *Elem)
{
  return init_mpole(ElemData, Elem, true, false, true, false);
}

static struct elem_type * init_bend_exact(const PyObject *ElemData, struct elem_type *Elem)
{
  return init_mpole(ElemData, Elem, true, false, false, true);
}

static struct elem_type * init_cbend(const PyObject *ElemData, struct elem_type *Elem)
{
  return init_mpole(ElemData, Elem, true, true, false, false);
}

/*
 * Following lines uses polymorphims.
 *
 * Could profit if the called functions could be renamed.
 * Or one had to stuck the functions listed here in a separate name space
 */
static struct elem_type * init_mpole(const PyObject *ElemData, struct elem_type *Elem)
{
    return init_mpole(ElemData, Elem, true, false, false, false);
}
static void MpolePassNoRad(double ps[], const int num_particles, struct elem_type *Elem)
{
  return MpolePass(ps, num_particles, Elem, false);
}

static void MpolePassRad(double ps[], const int num_particles, struct elem_type *Elem)
{
  return MpolePass(ps, num_particles, Elem, true);
}
static void MpolePass(double ps[], const int num_particles, struct elem_type *Elem)
{
  MpolePassNoRad(ps, num_particles, Elem);
}
static struct elem_type * init_wig(const PyObject *ElemData, struct elem_type *Elem)
{
  return init_wig(ElemData, Elem, false);
}
static struct elem_type * init_wig_rad(const PyObject *ElemData, struct elem_type *Elem)
{
  return init_wig(ElemData, Elem, true);
}



static struct elem_type * init_H_exact(const PyObject *ElemData, struct elem_type *Elem)
{
  return init_H(ElemData, Elem);
}

#undef ELEM_PASS
#define ELEM_INIT_FUNC_NAME(name) init_ ## name

#define ELEM_PASS(name, pass_name, api_identifier)   	                     \
__BEGIN_DECLS                                                                \
  FUNCDEFMACRO(name){					                     \
    std::cerr << "Start: " << __FUNCTION__ << " Elem " << Elem << std::endl; \
    if(!Elem){ std::cerr <<"Initalising elem to ";                           \
       Elem = ELEM_INIT_FUNC_NAME(name)(ElemData, Elem);                     \
       std::cerr << Elem << std::endl;}					     \
    if(!Elem){ return NULL;};						     \
      pass_name(ps, num_particles, Elem);                                    \
      std::cerr << "END: " << __FUNCTION__ << " returning "                  \
		<< Elem << std::endl;				             \
      return Elem;						             \
  }                                                                          \
__END_DECLS

/* use it here to define the functions themselves */
#include "element_integrators.h"
