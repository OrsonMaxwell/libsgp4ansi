#ifndef LIBSGP4ANSI_H_
#define LIBSGP4ANSI_H_

#include <time.h>
#include <inttypes.h>

#include "types.h"
#include "const.h"

// ************************************************************************* //
//                                 INTERFACE                                 //
// ************************************************************************* //

extern int
orbit_init(orbit*);

extern int
azel(void);

#endif /* LIBSGP4ANSI_H_ */
