// -*- mode: C++ -*-
#ifndef FLOW_UTIL_H
#define FLOW_UTIL_H
#include <cmath>
#ifndef M_PI
# define M_PI 3.14159265358979323846264338327
#endif

/** @defgroup u_utils Utilities 
    @brief Group of utility classes and functions */
//__________________________________________________________________
/** Normalize the angle @a ang to the interval @f$ [0,2\pi)@f$ 
    @ingroup u_utils
    @param ang Angle to normalize 
    @return the normalised angle */
inline Double_t 
NormalizeAngle(Double_t ang) 
{ 
  while (ang <  0)      ang += 2 * M_PI;
  while (ang >= 2*M_PI) ang -= 2 * M_PI;
  return ang;
}

#endif
//____________________________________________________________________
//
// EOF
//

