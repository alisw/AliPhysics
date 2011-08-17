// $Id:   $

/**************************************************************************
 * This file is property of and copyright by the Experimental Nuclear     *
 * Physics Group, Dep. of Physics                                         *
 * University of Oslo, Norway, 2006                                       *
 *                                                                        * 
 * Author: Federico Ronchetti								  *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//
//
//

#include "AliHLTCaloGeometry.h"
#include "unistd.h"
#include <iostream>

ClassImp(AliHLTCaloGeometry);

AliHLTCaloGeometry::AliHLTCaloGeometry(TString det) :  
  AliHLTCaloConstantsHandler(det)
  ,AliHLTLogging()
  ,fIsInitialised(kFALSE)
{  
  //see header file for class documentation
}


AliHLTCaloGeometry::~AliHLTCaloGeometry()
{

}





