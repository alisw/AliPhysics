///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2010
//
//    This file is part of starlight.
//
//    starlight is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    starlight is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with starlight. If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////
//
// File and Version Information:
// $Rev:: 213                         $: revision of last commit
// $Author:: butter                   $: author of last commit
// $Date:: 2015-08-15 23:08:02 +0200 #$: date of last commit
//
// Description:
//
//
//
///////////////////////////////////////////////////////////////////////////


#include <iostream>
#include <fstream>
#include <cmath>

#include "randomgenerator.h"


using namespace std;


//USED IN ROOT under TRANDOM3
// Random number generator class based on
//   M. Matsumoto and T. Nishimura,
//   Mersenne Twistor: A 623-diminsionally equidistributed
//   uniform pseudorandom number generator
//   ACM Transactions on Modeling and Computer Simulation,
//   Vol. 8, No. 1, January 1998, pp 3--30.
//
// For more information see the Mersenne Twistor homepage
//   http://www.math.keio.ac.jp/~matumoto/emt.html
//
// Advantage: large period 2**19937-1
//            relativly fast
//              (only two times slower than TRandom, but
//               two times faster than TRandom2)
// Drawback:  a relative large internal state of 624 integers
//
//
// Aug.99 ROOT implementation based on CLHEP by P.Malzacher
//
// the original code contains the following copyright notice:
/* This library is free software; you can redistribute it and/or   */
/* modify it under the terms of the GNU Library General Public     */
/* License as published by the Free Software Foundation; either    */
/* version 2 of the License, or (at your option) any later         */
/* version.                                                        */
/* This library is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of  */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.            */
/* See the GNU Library General Public License for more details.    */
/* You should have received a copy of the GNU Library General      */
/* Public License along with this library; if not, write to the    */
/* Free Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA   */
/* 02111-1307  USA                                                 */
/* Copyright (C) 1997 Makoto Matsumoto and Takuji Nishimura.       */
/* When you use this, send an email to: matumoto@math.keio.ac.jp   */
/* with an appropriate reference to your work.                     */
/////////////////////////////////////////////////////////////////////


void randomGenerator::SetSeed(unsigned int seed)
{
//  Set the random generator sequence
// if seed is 0 (default value) a TUUID is generated and used to fill
// the first 8 integers of the seed array.
// In this case the seed is guaranteed to be unique in space and time.
// Use upgraded seeding procedure to fix a known problem when seeding with values
// with many zero in the bit pattern (like 2**28).
// see http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MT2002/emt19937ar.html

  
  _count624 = 624;
  int i,j;
  
  _Mt[0] = seed;
  j = 1;
  // use multipliers from  Knuth's "Art of Computer Programming" Vol. 2, 3rd Ed. p.106
  for(i=j; i<624; i++) {
    _Mt[i] = (1812433253 * ( _Mt[i-1]  ^ ( _Mt[i-1] >> 30)) + i);
  }
}


double randomGenerator::Rndom(int)
{

//  Machine independent random number generator.
//  Produces uniformly-distributed floating points in [0,1]
//  Method: Mersenne Twistor

   
   unsigned int y;

   const int  kM = 397;
   const int  kN = 624;
   const unsigned int kTemperingMaskB =  0x9d2c5680;
   const unsigned int kTemperingMaskC =  0xefc60000;
   const unsigned int kUpperMask =       0x80000000;
   const unsigned int kLowerMask =       0x7fffffff;
   const unsigned int kMatrixA =         0x9908b0df;

   if (_count624 >= kN) {
      int i;

      for (i=0; i < kN-kM; i++) {
         y = (_Mt[i] & kUpperMask) | (_Mt[i+1] & kLowerMask);
         _Mt[i] = _Mt[i+kM] ^ (y >> 1) ^ ((y & 0x1) ? kMatrixA : 0x0);
      }

      for (   ; i < kN-1    ; i++) {
         y = (_Mt[i] & kUpperMask) | (_Mt[i+1] & kLowerMask);
         _Mt[i] = _Mt[i+kM-kN] ^ (y >> 1) ^ ((y & 0x1) ? kMatrixA : 0x0);
      }

      y = (_Mt[kN-1] & kUpperMask) | (_Mt[0] & kLowerMask);
      _Mt[kN-1] = _Mt[kM-1] ^ (y >> 1) ^ ((y & 0x1) ? kMatrixA : 0x0);
      _count624 = 0;
   }

   y = _Mt[_count624++];
   y ^=  (y >> 11);
   y ^= ((y << 7 ) & kTemperingMaskB );
   y ^= ((y << 15) & kTemperingMaskC );
   y ^=  (y >> 18);

   if (y) return ( (double) y * 2.3283064365386963e-10); // * Power(2,-32)
   return Rndom();
}
