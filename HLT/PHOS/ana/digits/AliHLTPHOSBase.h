#ifndef ALIHLTPHOSBASE_H
#define ALIHLTPHOSBASE_H

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
#include <Rtypes.h>
#include "TString.h"
#include "AliHLTDataTypes.h"
#include "AliHLTPHOSConstants.h"
#include "AliHLTPHOSConfig.h"
#include "AliHLTPHOSAltroConfig.h"

using namespace PhosHLTConst;
using namespace std;


class AliHLTPHOSConfig;
class AliHLTPHOSAltroConfig;

class  AliHLTPHOSBase : public AliHLTPHOSConfig, public AliHLTPHOSAltroConfig
{
public:
  AliHLTPHOSBase();
  virtual ~AliHLTPHOSBase();


  template<typename T> 
  void  DumpData(T *array, int N, int nPerLine)
  {
    cout <<   "DumpData N=  " << N <<endl;
    for(int i= 0; i< N; i++)
      {
	if((i%nPerLine == 0)  &&  (i != 0))
	  {
	    printf("\n");
	  }

	cout << array[i]<< "\t";

      }
  }

  template<typename T> 
  void  ResetArray(T *array, int N)
  {
    for(int i= 0; i< N; i++)
      {
	array[i] = 0;
      }
  }
 
  template<typename T> 
  T  MaxValue(T *array, int N)
  {
    //   Double_t tmpMax = 0;
      
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

  
protected:
  AliHLTPHOSConfig *fConfigPtr;
  AliHLTPHOSConfig *fAltroConfigPtr;

private:
  AliHLTPHOSBase(const AliHLTPHOSBase & );
  AliHLTPHOSBase & operator = (const AliHLTPHOSBase &);
};

#endif
