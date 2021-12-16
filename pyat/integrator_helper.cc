#include <map>

/**
 * Todo:
 *     Should be improved with std::mv
 */
#include "integrator-src/element_integrators_macros.h"

#define ELEM_PASS(name, pass_name, api_identifier) \
  lut[api_identifier] = ELEM_PASS_FUNC_NAME(name);

#include "integrator-src/ElemPass.h"

static std::map<std::string, track_function>  create_integrators_lookup_table(void){

  std::map<std::string, track_function> lut;

#include "integrator-src/element_integrators.h"

}
std::map<std::string, track_function>  integrators_lookup_table = create_integrators_lookup_table();
