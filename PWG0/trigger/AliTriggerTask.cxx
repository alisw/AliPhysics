/* $Id: AliTriggerTask.cxx 35782 2009-10-22 11:54:31Z jgrosseo $ */

#include "AliTriggerTask.h"

#include <TCanvas.h>
#include <TFile.h>
#include <TChain.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>

#include <AliLog.h>
#include <AliESDEvent.h>
#include <AliHeader.h>
#include <AliAnalysisManager.h>
#include <AliESDInputHandler.h>
#include <AliESDHeader.h>

ClassImp(AliTriggerTask)

AliTriggerTask::AliTriggerTask(const char* opt) :
  AliAnalysisTask("AliTriggerTask", ""),
  fESD(0),
  fOutput(0),
  fOption(opt),
  fStartTime(0),
  fEndTime(0),
  fNTriggers(0),
  fTriggerList(0),
  fStats(0)
{
  //
  // Constructor. Initialization of pointers
  //

  // Define input and output slots here
  DefineInput(0, TChain::Class());
  DefineOutput(0, TList::Class());
  
  fNTriggers = 13;
  
  static AliPWG0Helper::Trigger triggerList[] = { AliPWG0Helper::kAcceptAll, AliPWG0Helper::kFPANY, AliPWG0Helper::kMB1, AliPWG0Helper::kMB2, AliPWG0Helper::kMB3, AliPWG0Helper::kSPDGFO, AliPWG0Helper::kV0A, AliPWG0Helper::kV0C, AliPWG0Helper::kZDC, AliPWG0Helper::kZDCA, AliPWG0Helper::kZDCC, AliPWG0Helper::kFMDA, AliPWG0Helper::kFMDC };
  fTriggerList = triggerList;
  
  fStats = new TH1*[fNTriggers];
  
  AliLog::SetClassDebugLevel("AliTriggerTask", AliLog::kWarning);
}

AliTriggerTask::~AliTriggerTask()
{
  //
  // Destructor
  //

  // histograms are in the output list and deleted when the output
  // list is deleted by the TSelector dtor

  if (fOutput) {
    delete fOutput;
    fOutput = 0;
  }
}

//________________________________________________________________________
void AliTriggerTask::ConnectInputData(Option_t *)
{
  // Connect ESD
  // Called once

  Printf("AliTriggerTask::ConnectInputData called");

  AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

  if (!esdH) {
    Printf("ERROR: Could not get ESDInputHandler");
  } else {
    fESD = esdH->GetEvent();
    
    TString branches("AliESDHeader Vertex AliMultiplicity ALIESDVZERO ALIESDZDC FMD");

    // Enable only the needed branches
    esdH->SetActiveBranches(branches);
  }
}

void AliTriggerTask::CreateOutputObjects()
{
  // create result objects and add to output list

  Printf("AliTriggerTask::CreateOutputObjects");

  fOutput = new TList;
  fOutput->SetOwner();
  
  if (fStartTime == fEndTime)
    AliWarning("Start and endtime not set. Automatic binning will be used. This does not work in parallel systems");

  Int_t nBins = 1000;
  if (fEndTime - fStartTime > 0)
    nBins = fEndTime - fStartTime;
  for (Int_t i=0; i<fNTriggers; i++)
  {
    fStats[i] = new TH1F(Form("fStats_%d", i), Form("%s;time;counts", AliPWG0Helper::GetTriggerName(fTriggerList[i])), nBins, 0, fEndTime - fStartTime);
    fOutput->Add(fStats[i]);
  }
}

void AliTriggerTask::Exec(Option_t*)
{
  // process the event

  // post the data already here
  PostData(0, fOutput);

  if (!fESD)
  {
    AliError("ESD branch not available");
    return;
  }
  
  // check event type (should be PHYSICS = 7)
  AliESDHeader* esdHeader = fESD->GetHeader();
  if (!esdHeader)
  {
    Printf("ERROR: esdHeader could not be retrieved");
    return;
  }
  
  UInt_t eventType = esdHeader->GetEventType();
  if (eventType != 7)
  {
    Printf("Skipping event because it is of type %d", eventType);
    return;
  }

  //Printf("Trigger classes: %s:", fESD->GetFiredTriggerClasses().Data());
  
  UInt_t timeStamp = fESD->GetTimeStamp() - fStartTime;
  //Printf("%d", timeStamp);
  
  for (Int_t i = 0; i < fNTriggers; i++)
  {
    Bool_t triggered = AliPWG0Helper::IsEventTriggered(fESD, (AliPWG0Helper::Trigger) (fTriggerList[i] | AliPWG0Helper::kOfflineFlag));
    if (triggered)
      fStats[i]->Fill(timeStamp);
    //Printf("%s: %d", AliPWG0Helper::GetTriggerName(fTriggerList[i]), triggered);
  }
}

void AliTriggerTask::Terminate(Option_t *)
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

  fOutput = dynamic_cast<TList*> (GetOutputData(0));
  if (!fOutput)
    Printf("ERROR: fOutput not available");

  if (fOutput)
  {
    for (Int_t i=0; i<fNTriggers; i++)
      fStats[i] = dynamic_cast<TH1*> (fOutput->FindObject(Form("fStats_%d", i)));
  }

  TFile* fout = new TFile("trigger.root", "RECREATE");

  for (Int_t i=0; i<fNTriggers; i++)
    if (fStats[i])
      fStats[i]->Write();

  fout->Write();
  fout->Close();
  
  Int_t nX = (Int_t) TMath::Sqrt(fNTriggers);
  Int_t nY = nX;
  
  while (nX * nY < fNTriggers)
  {
    if (nX == nY)
      nX++;
    else
      nY++;
  }
  
  TCanvas* c = new TCanvas("c", "c", 800, 800);
  c->Divide(nX, nY);
  
  Printf("+++++++++ TRIGGER STATS:");

  Int_t base = 1;
  if (fStats[0])
    base = (Int_t) fStats[0]->Integral();

  Int_t length = fEndTime - fStartTime;
  
  for (Int_t i=0; i<fNTriggers; i++)
    if (fStats[i])
    {
      c->cd(i+1);
      fStats[i]->Draw();
      Printf("%s: %d triggers | %f %% of all triggered | Rate: %f Hz", AliPWG0Helper::GetTriggerName(fTriggerList[i]), (UInt_t) fStats[i]->Integral(), 100.0 * fStats[i]->Integral() / base, (length > 0) ? (fStats[i]->Integral() / length) : -1);
    }

  Printf("Writting result to trigger.root");
}
