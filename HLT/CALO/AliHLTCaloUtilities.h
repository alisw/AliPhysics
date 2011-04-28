//-*- Mode: C++ -*-
// $Id: AliHLTCALOUtilities.h 34264 2009-08-14 18:29:23Z odjuvsla $

#ifndef ALIHLTCALOUTILITIES_H
#define ALIHLTCALOUTILITIES_H

/**************************************************************************
 * This file is property of and copyright by the Experimental Nuclear     *
 * Physics Group, Dep. of Physics                                         *
 * University of Oslo, Norway, 2007                                       *
 *                                                                        *
 * Author: Per Thomas Hille <perthi@fys.uio.no> for the ALICE HLT Project.*
 * Contributors are mentioned in the code where appropriate.              *
 * Please report bugs to perthi@fys.uio.no                                *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include <iostream>

using namespace std;

class  AliHLTCaloUtilities
{
 public:
  AliHLTCaloUtilities();
  virtual ~AliHLTCaloUtilities();
  
  template<typename T> 
    static T  MaxValue(T *array, int N)
    {
      T tmpMax = 0;

      for(int i = 0; i < N; i++)
	{
	  if(array[i] > tmpMax)
	    {
	      tmpMax = array[i];
	    }
	}
      return tmpMax;
    }
  
};

#endif
