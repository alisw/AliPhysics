// -*- mode: C++ -*-
/* Copyright (C) 2007 Christian Holm Christensen <cholm@nbi.dk>
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1 of
 * the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 * USA
 */
#ifndef ALIFMDFLOWUTIL_H
#define ALIFMDFLOWUTIL_H
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

