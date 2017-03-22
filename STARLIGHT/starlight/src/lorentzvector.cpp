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
// $Rev:: 28                          $: revision of last commit
// $Author:: bgrube                   $: author of last commit
// $Date:: 2010-12-10 19:30:01 +0100 #$: date of last commit
//
// Description:
//
//
//
///////////////////////////////////////////////////////////////////////////


#include "lorentzvector.h"


lorentzVector::lorentzVector() :
fSpaceVec()
,fTime(0)
{ }


lorentzVector::lorentzVector(double x, double y, double z, double t) :
fSpaceVec(x, y, z)
,fTime(t)
{ }


lorentzVector::~lorentzVector()
{ }


void lorentzVector::SetXYZT(double x, double y, double z, double t)
{
   fSpaceVec.SetVector(x, y, z);
   fTime  = t;
}
