/**************************************************************************
 * Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
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
#include <iostream>
#include <TString.h>

#include "AliESDInputHandler.h"
#include "AliAnalysisTaskDummy.h"

ClassImp(AliAnalysisTaskDummy)

AliAnalysisTaskDummy::AliAnalysisTaskDummy() :
  fDebugLevel(0)
{

}

AliAnalysisTaskDummy::~AliAnalysisTaskDummy() {
}

void AliAnalysisTaskDummy::UserExec(Option_t *){
  TString filename = "";
  AliESDInputHandler *inh = dynamic_cast<AliESDInputHandler *>(fInputHandler);
  if(inh) filename = inh->GetInputFileName();
  if(fDebugLevel){
    std::cout << "Running User exec from task AliAnalysisTaskDummy on ESDtree from file" << filename << std::endl;
  }
}
