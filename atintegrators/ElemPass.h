#ifndef _PyAT_ELEMPASS_H_
#define _PyAT_ELEMPASS_H_ 1

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

/**
 *  Todo:
 *   include definition of parameters of the FUNC DEC MACRO
 */
/*
 * need to solve header organisation
 */
// #include "element_integrator_macros.h"
#define FUNCDEFMACRO(name)			\
    struct elem_type*  \
    ELEM_PASS_FUNC_NAME(name)(const PyObject *ElemData, struct elem_type *Elem, double ps[], \
              const int num_particles, const struct parameters *Param)
#define HEADERMACRO(name) FUNCDEFMACRO(name) ;


#undef ELEM_PASS
#define ELEM_PASS(name, pass_name, api_identifier) HEADERMACRO(name)
/* use it here to make forward declarations for the functions */
#include "element_integrators.h"
__END_DECLS

#endif /* _PyAT_ELEMPASS_H_ */
