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
#include <chrono>
#include <thread>
#include <TString.h>

#include "AliInputEventHandler.h"
#include "AliAnalysisTaskDummy.h"

ClassImp(AliAnalysisTaskDummy)

void AliAnalysisTaskDummy::UserExec(Option_t *){
  TString filename = "";
  if(fInputHandler) filename = fInputHandler->GetInputFileName();
  if(fDebugLevel){
    std::cout << "Running User exec from task AliAnalysisTaskDummy on file" << filename << std::endl;
  }
  std::chrono::milliseconds interval(fWaitTime);
  std::this_thread::sleep_for(interval);

}
