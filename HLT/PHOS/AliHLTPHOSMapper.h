#ifndef ALIHLTPHOSMAPPER_H
#define ALIHLTPHOSMAPPER_H

/**************************************************************************
 * This file is property of and copyright by the Experimental Nuclear     *
 * Physics Group, Dep. of Physics                                         *
 * University of Oslo, Norway, 2006                                       *
 *                                                                        * 
 * Author: Per Thomas Hille perthi@fys.uio.no for the ALICE DCS Project.  *
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

//#include "PhosFeeClient.h"

//#include "stdio.h"
//#include <iostream>
#include <cstdlib>
#include <assert.h>
//#include "AliHLTPHOSCommonDefs.h"
//#include "AliHLTPHOSConstants.h"
#include "AliHLTPHOSBase.h"
    
//            PhosHLTConst
using namespace std;
using namespace PhosHLTConst;

class AliHLTPHOSMapper : public AliHLTPHOSBase
{
 public:
  AliHLTPHOSMapper();
  virtual ~AliHLTPHOSMapper();
  void InitAltroMapping(); 

  struct altromap{ 
    //    int mod;
    int zRow;
    int xCol;
    int gain;
    //    int rcu;
    //    int branch;
    //    int card;
    //    int chip;
    //    int chan;
    //    int csp;
    //    int num;
    //    int hid;
  };


altromap *hw2geomapPtr;

};

#endif
