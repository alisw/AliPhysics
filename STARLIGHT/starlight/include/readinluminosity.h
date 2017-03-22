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


#ifndef READINLUMINOSITY_H
#define READINLUMINOSITY_H


#include "inputParameters.h"
#include "starlightlimits.h"


class readLuminosity
{
 public:
  readLuminosity(const inputParameters& input);
  ~readLuminosity();
  
  void read();
  double *_Warray;
  double *_Yarray;
  double **_Farray; 
  double **_Farray1;
  double **_Farray2;
  
  double _f_max;
  double _f_max1;
  double _f_max2;

  double _fptarray[500][500];

  double _bwnormsave;

 protected:
  const int _ReadInputNPT;
  const int _ReadInputnumy;
  const int _ReadInputnumw;
  const int _ReadInputgg_or_gP;
  const int _ReadInputinterferencemode;
  const std::string _baseFileName;
};


#endif  // READINLUMINOSITY_H
