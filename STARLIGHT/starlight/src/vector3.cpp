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


#include "vector3.h"


vector3::vector3() 
{
   _vec[0] = 0;
   _vec[1] = 0;
   _vec[2] = 0;
}


vector3::vector3(double x, double y, double z)
{
   _vec[0] = x;
   _vec[1] = y;
   _vec[2] = z;
}


vector3::~vector3()
{ }


void vector3::SetVector(double x, double y, double z)
{
   _vec[0] = x;
   _vec[1] = y;
   _vec[2] = z;
}


void vector3::SetVector(double *vec)
{
   _vec[0] = vec[0];
   _vec[1] = vec[1];
   _vec[2] = vec[2];
}
