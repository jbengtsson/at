#ifndef _PyAT_UTILS_H_
#define _PyAT_UTILS_H_ 1

/**
 * Todo:
 *     unify with header in thor_scsi/python
 */
#define STRINGIFY(a_var) #a_var
#define TOSTRING(a_str) STRINGIFY(a_str)


#define DISPLAY_VAR(a_var) #a_var  " = " STRINGIFY(a_var)
#define DISPLAY_AT() __FILE__ ":" TOSTRING(__LINE__)
#define DISPLAY_VAR_AT(a_var)  DISPLAY_AT() " "  DISPLAY_VAR(NO_TPSA)

#endif
