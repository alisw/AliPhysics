/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// ALICE Reconstruction parameterization:                                    //
//                                                                           //
//                                                                           //
// Base Class for Detector reconstruction parameters                         //
// Revision: cvetan.cheshkov@cern.ch 12/06/2008                              //
// Its structure has been revised and it is interfaced to AliEventInfo.      //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TObjArray.h"
#include "AliDetectorRecoParam.h"

#include "AliRecoParam.h"


ClassImp(AliRecoParam)

AliRecoParam::AliRecoParam(): 
  TNamed("",""),
  fRecoParamArray(0)
{
  // Default constructor
  // ...
}

AliRecoParam::AliRecoParam(const char *detector): 
  TNamed(detector,detector),
  fRecoParamArray(0)
{
  // Default constructor
  // ...
}

AliRecoParam::~AliRecoParam(){
  // Destructor
  // ...
  // Delete the array with the reco-param objects
  if (fRecoParamArray){
    fRecoParamArray->Delete();
    delete fRecoParamArray;
  }
}

void  AliRecoParam::Print(Option_t *option) const {
  //
  // Print reconstruction setup
  //
  printf("AliRecoParam object for %s\n",GetName()); 
  if (!fRecoParamArray) return;
  Int_t nparam = fRecoParamArray->GetEntriesFast();
  for (Int_t iparam=0; iparam<nparam; iparam++){
    AliDetectorRecoParam * param = (AliDetectorRecoParam *)fRecoParamArray->At(iparam);
    if (!param) continue;
    param->Print(option);
  }
}

void  AliRecoParam::AddRecoParam(AliDetectorRecoParam* param){
  // Add an instance of reco params object into
  // the fRecoParamArray
  //
  if (!fRecoParamArray) fRecoParamArray = new TObjArray;
  fRecoParamArray->AddLast(param);
}
