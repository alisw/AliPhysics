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
// $Rev:: 228                         $: revision of last commit
// $Author:: butter                   $: author of last commit
// $Date:: 2016-01-18 18:10:17 +0100 #$: date of last commit
//
// Description:
//    Bessel functions taken from ROOT
//
//
///////////////////////////////////////////////////////////////////////////


#ifndef BESSEL_H
#define BESSEL_H


namespace bessel
{
		double besI0(double x);
		double dbesk0(double x);
		double dbesk1(double x);
		double besI1(double x);
};


#endif  // BESSEL_H
