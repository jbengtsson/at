#include <map>

/**
 * Todo:
 *     Should be improved with std::mv
 */

#include "integrator-src/element_integrators_macros.h"
#include "integrator-src/ElemPass.h"



//    std::cerr << "Adding Integrator " << tmp <<  " using '"  << ((void*)f) << "'" <<std::endl;

#define STRINGIFY(a_var) #a_var
#define TOSTRING(a_str) STRINGIFY(a_str)

#undef ELEM_PASS
#define ELEM_PASS(name, pass_name, api_identifier)                    \
{                                                                     \
    std::string ai = api_identifier;			              \
    std::string tmp  = (ai.size() == 0) ? TOSTRING(pass_name) : ai;   \
    track_function f = ( (track_function) ELEM_PASS_FUNC_NAME(name)); \
    lut[tmp] = f;						      \
    /* lut[tmp] = TOSTRING(name);  */				      \
}

static std::map<std::string, track_function>  create_integrators_lookup_table(void)
{
  std::map<std::string, track_function> lut;
  std::map<std::string, std::string> lut2;
  std::string tmp;
  #include "integrator-src/element_integrators.h"
#if 0
  std::cerr << "Known integrators ";
  for(auto it = lut.begin(); it != lut.end(); ++it){
    std::cerr << " Name : " << it->first << ";";
    ;
  }
#endif
  std::cerr << std::endl;
  std::cerr.flush();

  return lut;
}
std::map<std::string, track_function>  integrators_lookup_table = create_integrators_lookup_table();

/**
 * Return function pointer for given (integrator) function name
 *
 * If not found prints known integrators to stderr
 *
 * Can this be simplified ?
 */
track_function lookup_function(std::string fn_name)
{

  std::string buffer;

  auto cnt = integrators_lookup_table.count(fn_name);
  switch(cnt){
  case 1:
    /* found one name ... fine! */
    break;
  case 0:
    std::cerr << "Known integrators ";
    for(auto it = integrators_lookup_table.begin(); it != integrators_lookup_table.end(); ++it){
      std::cerr << " Name : " << it->first << ";";
      ;
    }
    std::cerr << std::endl;
    std::cerr.flush();
    buffer = "PassMethod '" + fn_name + "': library, module or trackFunction not found";
    PyErr_SetString(PyExc_RuntimeError, buffer.c_str());
    return NULL;
    break;
  default:
    {
      std::stringstream msg;
      msg <<  "Sanity Error: PassMethod " << fn_name << ": found " <<  cnt << " times. Expected only one!";
      PyErr_SetString(PyExc_RuntimeError, msg.str().c_str());
    }
    return NULL;
  }

  return integrators_lookup_table[fn_name];
}
