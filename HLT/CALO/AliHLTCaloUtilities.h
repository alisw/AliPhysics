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

//      AliHLTCaloUtilities
class  AliHLTCaloUtilities
{
 public:
  AliHLTCaloUtilities();
  virtual ~AliHLTCaloUtilities();
 
  static bool CheckFile(const char *fileName, const char *opt);

  //  static bool ScanSingleIntArgument(int argc, const char** argv, const char *name, int *value)
  static bool ScanSingleIntArgument(int argc, const char** argv, const char *name, int *value = 0 );
  static bool ScanSingleFloatArgument(int argc, const char** argv, const char *name, float *value = 0 );
  static bool ScanSingleNameArgument(int argc, const char** argv, const char *name, char *outname = 0 );
  static bool ScanSingleArgument(int argc, const char** argv, const char *name);

  template<typename T> 
    static void  DumpData(T *array, int N, int nPerLine)
    {
      //   cout <<   "DumpData N=  " << N <<endl;
      for(int i= 0; i< N; i++)
	{
	  if((i%nPerLine == 0)  &&  (i != 0))
	    {
	      //printf("\n");
	      cout << endl;
	    }

	  cout << array[i]<< "\t";
	}
      //     printf("\n");
      cout << endl;
    }

  template<typename T> 
    static void  ResetArray(T *array, int N)
    {
      for(int i= 0; i< N; i++)
	{
	  array[i] = 0;
	}
    }
 
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



  
 private:
  static int DoExistArgument(const int argc, const char** argv, const char *argument);

};

#endif
