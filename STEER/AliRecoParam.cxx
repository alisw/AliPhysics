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
// ALICE Reconstruction parameterization:
// 
//
// Retrieving paramaters:
// 0.  Read the parameters from database
// 1.  Using the ConfigRecoParam.C script (example : $ALICE_ROOT/macros)
// 2.  Register additional parametes (AliRecoParam::Instance()->RegisterRecoParam(tpcRecoParam);
//
//
// Using the reconstruction parameters
//    AliRecoParam::Instance()->GetRecoParam(detType, eventType)                                                                       //  
//  detType: 
//  1. Detectors - ITS, TPC, TRD, TOF ...
//  2. Process   - V0, Kink, ESDcuts

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TObjArray.h"
#include "AliDetectorRecoParam.h"

#include "AliRecoParam.h"


ClassImp(AliRecoParam)


AliRecoParam* AliRecoParam::fgInstance = 0x0;


//_____________________________________________________________________________
AliRecoParam* AliRecoParam::Instance()
{
  //
  // returns AliRecoParam instance (singleton)
  //
  if (!fgInstance) {
    fgInstance = new AliRecoParam();
  }
  
  return fgInstance;
}



AliRecoParam::AliRecoParam(): 
  TNamed("ALICE","ALICE"),
  fRecoParamArray(0)
{
  //
  //
  //
}

AliRecoParam::~AliRecoParam(){
  //
  //
  //
  if (fRecoParamArray){
    fRecoParamArray->Delete();
    delete fRecoParamArray;
  }
}

void  AliRecoParam::Print(Option_t *option) const{
  //
  // Print reconstruction setup
  //
  printf("AliRecoParam\n"); 
  if (!fRecoParamArray) return;
  Int_t nparam = fRecoParamArray->GetEntriesFast();
  for (Int_t iparam=0; iparam<nparam; iparam++){
    AliDetectorRecoParam * param = (AliDetectorRecoParam *)fRecoParamArray->At(iparam);
    if (!param) continue;
    param->Print(option);
  }
}

void  AliRecoParam::RegisterRecoParam(AliDetectorRecoParam* param){
  //
  //
  //
  if (!fRecoParamArray) fRecoParamArray = new TObjArray;
  fRecoParamArray->AddLast(param);
}

TObjArray * AliRecoParam::GetRecoParam(const char * detType, Int_t *eventType){
  //
  // Get the list of Reconstruction parameters for given detector
  // and event type
  //
  if (!fRecoParamArray) return 0;
  TObjArray * array = 0;
  Int_t nparam = fRecoParamArray->GetEntriesFast();
  for (Int_t iparam=0;iparam<nparam; iparam++){
    AliDetectorRecoParam * param = (AliDetectorRecoParam *)fRecoParamArray->At(iparam);
    if (!param) continue;
    TString str(param->GetName());
    if (!str.Contains(detType)) continue;
    if (!array) array = new TObjArray;
    array->AddLast(param);    
  }
  return array;
}


