/**************************************************************************
 * Copyright(c) 2013-2014, ALICE Experiment at CERN, All rights reserved. *
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

#ifdef __ECLIPSE_IDE
//  few includes and external declarations just for the IDE
#include "Riostream.h"
#include "TSystem.h"
#include "AliAnalysisTaskFlowVectorCorrections.h"
#include "AliAnalysisManager.h"
AliAnalysisDataContainer* AddTaskFlowQnVectorCorrections();

#include "runAnalysis.H"

#endif // ifdef __ECLIPSE_IDE declaration and includes for the ECLIPSE IDE

using std::cout;
using std::endl;

void AddTaskFlowQnVectorCorrectionsToLegoTrain(const char *configpath = ".") {

  /* strange way of including the header file is for lego train scenarios */
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/EVCHAR/FlowVectorCorrections/QnCorrectionsInterface/macros/runAnalysis.H");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/EVCHAR/FlowVectorCorrections/QnCorrectionsInterface/macros/loadRunOptions.C");
  if (!loadRunOptions(kFALSE, configpath)) {
    cout << "ERROR: configuration options not loaded. ABORTING!!!" << endl;
    return;
  }

  gSystem->AddIncludePath("-I$ALICE_PHYSICS/include");

  gSystem->Load("libPWGPPevcharQn.so");
  gSystem->Load("libPWGPPevcharQnInterface.so");

  AliAnalysisManager *mgr;
  if (!bTrainScope) {
    cout << "This macro shall not be used outside the lego train. Aborting" << endl;
    return;
  }
  else {
    mgr = AliAnalysisManager::GetAnalysisManager();
  }

  gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/EVCHAR/FlowVectorCorrections/QnCorrectionsInterface/macros/AddTaskFlowQnVectorCorrections.C");
  AliAnalysisDataContainer *corrTask = AddTaskFlowQnVectorCorrections();

  if (bRunQnVectorAnalysisTask) {
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/EVCHAR/FlowVectorCorrections/QnCorrectionsInterface/macros/AddTaskQnVectorAnalysis.C");
    AliAnalysisTaskQnVectorAnalysis* taskQn = AddTaskQnVectorAnalysis(bUseMultiplicity, b2015DataSet);
    taskQn->SetExpectedCorrectionPass(szCorrectionPass.Data());
    taskQn->SetAlternativeCorrectionPass(szAltCorrectionPass.Data());

    mgr->AddTask(taskQn);

    //create output container
    AliAnalysisDataContainer *cOutputQnAnaEventQA =
      mgr->CreateContainer("QnAnalysisEventQA",
          TList::Class(),
          AliAnalysisManager::kOutputContainer,
          "QnAnalysisEventQA.root");

    mgr->ConnectInput(taskQn,  0, mgr->GetCommonInputContainer());
    mgr->ConnectInput(taskQn,  1, corrTask);
    mgr->ConnectOutput(taskQn, 1, cOutputQnAnaEventQA );
  }
}

