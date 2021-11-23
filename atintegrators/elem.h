#include <math.h>


/*----------------------------------------------------*/
/*            For the integrator code                 */
/*----------------------------------------------------*/

#define atIsFinite isfinite
#define atIsNaN    isnan
#define atGetNaN() (NAN)
#define atGetInf() (INFINITY)
#define atMalloc   malloc
#define atCalloc   calloc
#define atFree     free

#define SQR(X) ((X)*(X))
