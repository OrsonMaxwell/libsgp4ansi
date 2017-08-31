#ifndef LIBTLE2TCP_H_
#define LIBTLE2TCP_H_

#include <time.h>
#include <inttypes.h>

#include "types.h"
#include "const.h"

// ************************************************************************* //
//                                 INTERFACE                                 //
// ************************************************************************* //

extern int
orbit_init(sat*);

extern int
azel(void);

#endif /* LIBTLE2TCP_H_ */
