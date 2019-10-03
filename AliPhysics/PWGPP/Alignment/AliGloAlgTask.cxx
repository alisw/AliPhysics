/*************************************************************************
* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
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


#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TSystem.h>
#include <TROOT.h>
#include "AliAnalysisManager.h"
#include "AliESDInputHandler.h"
#include "AliGloAlgTask.h"
#include "AliLog.h"
#include "AliAlgSteer.h"
#include "AliAlgDOFStat.h"

ClassImp(AliGloAlgTask)

//________________________________________________________________________
AliGloAlgTask::AliGloAlgTask(const char *name) 
  : AliAnalysisTaskSE(name)
  ,fOutput(0)
  ,fTrigSel(0)
  ,fAlgSteer(0)
  ,fIniParFileName()
  ,fConfMacroName()
  ,fStopWatch()
  ,fChunks(0)
  ,fApplyMPSolAlignment(kFALSE)
{
  // Constructor
  DefineOutput(1, TList::Class());
  SetConfMacroName();
  //
}

//________________________________________________________________________
AliGloAlgTask::~AliGloAlgTask()
{
  // Destructor
  if (fOutput && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {  //RRR
    printf("Deleteing output\n");
    delete fOutput;
    fOutput = 0;
  }
  //
  delete fAlgSteer;
  //
}

//________________________________________________________________________
void AliGloAlgTask::UserCreateOutputObjects() 
{
  //
  fOutput = new TList();
  fOutput->SetOwner(); 
  //
  PostData(1, fOutput);
  //
}

//________________________________________________________________________
void AliGloAlgTask::NotifyRun()
{
  // attach tree to alignment object
  AliAnalysisManager* anMan = AliAnalysisManager::GetAnalysisManager();
  AliESDInputHandler *handler = (AliESDInputHandler*)anMan->GetInputEventHandler();
  fAlgSteer->SetESDTree(handler->GetTree());
}

//________________________________________________________________________
Bool_t AliGloAlgTask::Notify()
{
  // notify chunk change
  fStopWatch.Stop();
  AliInfoF("Processing chunk %d",fChunks++);
  fStopWatch.Print();
  fStopWatch.Start(kFALSE);
  return kTRUE;
}

//________________________________________________________________________
void AliGloAlgTask::UserExec(Option_t *) 
{
  // Main loop
  //
  AliAnalysisManager* anMan = AliAnalysisManager::GetAnalysisManager();
  AliESDInputHandler *handler = (AliESDInputHandler*)anMan->GetInputEventHandler();
  AliESDEvent* esdEv = (AliESDEvent*)handler->GetEvent();
  AliESDfriend *esdFr = handler->GetESDfriend(); // get the input friend
  //
  if(!esdEv) {AliInfo("no ESD"); return;} 
  if(!esdFr || esdFr->GetNumberOfTracks()<esdEv->GetNumberOfTracks()) {
    AliDebug(3,"no ESDFriend"); 
    return;
  }
  //  AliInfo(Form("Number of ESD tracks in input = %d ",esdEv->GetNumberOfTracks()));
  //  AliInfo(Form("Number of tracks in input friends = %d ",esdFr->GetNumberOfTracks()));
  //
  fAlgSteer->ProcessEvent(esdEv);
  //
}      

//________________________________________________________________________
void AliGloAlgTask::Terminate(Option_t *) 
{
  Printf("Terminating...");
  fAlgSteer->Terminate();
  //
  TH1* hstat = 0;
  AliAlgDOFStat* dofSt= fAlgSteer->GetDOFStat();
  if (dofSt) fOutput->Add(dofSt);
  fAlgSteer->DetachDOFStat();
  //
  hstat = fAlgSteer->GetHistoStat();
  if (hstat) fOutput->Add(hstat);
  fAlgSteer->DetachHistoStat();
  //
  fStopWatch.Stop();
  AliInfoF("Processed %d chunks",fChunks);
  fStopWatch.Print();
  //
}

//_________________________________________________________________________
void AliGloAlgTask::LocalInit()
{
  // init alignment object by running user provided macro containing a function with
  // the same name as the macro and accepting as an argument a pointer to fAlgSteer
  AliInfo("Processing");
  if (fAlgSteer) AliFatal("Something is wrong, alignment was already initialized");
  //
  if (fConfMacroName.IsNull() || gSystem->AccessPathName(fConfMacroName.Data(), kFileExists)) {
    AliFatalF("Cannot find configuration macro %s",fConfMacroName.Data());
  }
  //
  fAlgSteer = new AliAlgSteer();
  gROOT->ProcessLine(Form(".x %s+g((AliAlgSteer*)%p)",fConfMacroName.Data(),fAlgSteer));
  //
  if (!fIniParFileName.IsNull() && !gSystem->AccessPathName(fIniParFileName.Data(), kFileExists)) {
    AliInfoF("Imposing initial parameters from %s",fIniParFileName.Data());
    fAlgSteer->ReadParameters(fIniParFileName.Data());
    if (fApplyMPSolAlignment) fAlgSteer->ApplyAlignmentFromMPSol();
  }
  fAlgSteer->Print();
}
