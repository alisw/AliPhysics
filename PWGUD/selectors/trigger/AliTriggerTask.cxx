/* $Id: AliTriggerTask.cxx 35782 2009-10-22 11:54:31Z jgrosseo $ */

#include "AliTriggerTask.h"
BREAKMEVERYBADLY!
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
#include <AliPhysicsSelection.h>

ClassImp(AliTriggerTask)

AliTriggerTask::AliTriggerTask(const char* opt) :
  AliAnalysisTask("AliTriggerTask", ""),
  fESD(0),
  fOutput(0),
  fOption(opt),
  fStartTime(0),
  fEndTime(0),
  fUseOrbits(kFALSE),
  fFirstOrbit(0),
  fLastOrbit(0),
  fNTriggers(0),
  fTriggerList(0),
  fNTriggerClasses(0),
  fTriggerClassesList(0),
  fStats(0),
  fPhysicsSelection(0)
{
  //
  // Constructor. Initialization of pointers
  //

  // Define input and output slots here
  DefineInput(0, TChain::Class());
  DefineOutput(0, TList::Class());
  
  fNTriggers = 14;
  
  static AliTriggerAnalysis::Trigger triggerList[] = { AliTriggerAnalysis::kAcceptAll, AliTriggerAnalysis::kFPANY, AliTriggerAnalysis::kMB1, AliTriggerAnalysis::kMB2, AliTriggerAnalysis::kMB3, AliTriggerAnalysis::kSPDGFO, AliTriggerAnalysis::kSPDGFOBits, AliTriggerAnalysis::kV0A, AliTriggerAnalysis::kV0C, AliTriggerAnalysis::kZDC, AliTriggerAnalysis::kZDCA, AliTriggerAnalysis::kZDCC, AliTriggerAnalysis::kFMDA, AliTriggerAnalysis::kFMDC };
  fTriggerList = triggerList;
  
  fNTriggerClasses = 4;
  static const char* triggerClassesList[] = { "CINT1B-ABCE-NOPF-ALL", "CINT1C-ABCE-NOPF-ALL", "CINT1A-ABCE-NOPF-ALL", "CINT1-E-NOPF-ALL" };
  fTriggerClassesList = triggerClassesList;
  
  fStats = new TH1**[fNTriggerClasses];
  for (Int_t i=0; i<fNTriggerClasses; i++)
    fStats[i] = new TH1*[fNTriggers];
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
    fESD = (AliESDEvent*)esdH->GetEvent();
    
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
    nBins = fEndTime - fStartTime + 1;
  if (nBins > 10000)
    nBins = 10000;
  
  Int_t start = 0;
  Int_t end = fEndTime - fStartTime;
  
  if (fUseOrbits)
  {
    start = fStartTime;
    end = fEndTime;
  }
  
  for (Int_t j=0; j<fNTriggerClasses; j++)
  {
    for (Int_t i=0; i<fNTriggers; i++)
    {
      fStats[j][i] = new TH1F(Form("fStats_%d_%d", j, i), Form("%s %s;%s;counts", fTriggerClassesList[j], AliTriggerAnalysis::GetTriggerName(fTriggerList[i]), (fUseOrbits) ? "orbit number" : "time"), nBins, start - 0.5, end + 0.5);
      fOutput->Add(fStats[j][i]);
    }
  }
    
  fFirstOrbit = new TParameter<Long_t> ("fFirstOrbit", 0);
  fLastOrbit = new TParameter<Long_t> ("fLastOrbit", 0);
  fOutput->Add(fFirstOrbit);
  fOutput->Add(fLastOrbit);
  
  fOutput->Add(fPhysicsSelection);
  
  AliLog::SetClassDebugLevel("AliPhysicsSelection", AliLog::kDebug);
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
  
  // fill histograms
  fPhysicsSelection->IsCollisionCandidate(fESD);
  
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
  
  // TODO select on hardware trigger for histograms...
  Int_t triggerClass = 0;
  while (triggerClass < fNTriggerClasses && !fESD->IsTriggerClassFired(fTriggerClassesList[triggerClass]))
    triggerClass++;
  if (triggerClass == fNTriggerClasses)
  {
    Printf("Unknown trigger class %s. Skipping event", fESD->GetFiredTriggerClasses().Data());
    return;
  }
  
  Long64_t timeStamp = 0;
  if (fUseOrbits)
  {
    timeStamp = fESD->GetBunchCrossNumber();
    timeStamp += (Long64_t) 3564 * (fESD->GetOrbitNumber() + fESD->GetPeriodNumber() * 16777215);
    timeStamp = (Long64_t) (25e-9 * timeStamp);
    timeStamp -= fStartTime;
  }
  else
    timeStamp = fESD->GetTimeStamp() - fStartTime;
    
  //Printf("%d", timeStamp);
  
  //Annalisa Time (s) = 1440*period + 88*10-6 * orbit + 25*10-9 *bc
  
  for (Int_t i = 0; i < fNTriggers; i++)
  {
    Bool_t triggered = fPhysicsSelection->GetTriggerAnalysis()->IsOfflineTriggerFired(fESD, fTriggerList[i]);
    if (triggered)
      fStats[triggerClass][i]->Fill(timeStamp);
    //Printf("%s: %d", AliTriggerAnalysis::GetTriggerName(fTriggerList[i]), triggered);
  }
  
  if (fFirstOrbit->GetVal() == 0)
    fFirstOrbit->SetVal(timeStamp);
  else
    fFirstOrbit->SetVal(TMath::Min(fFirstOrbit->GetVal(), (Long_t) timeStamp));
  
  fLastOrbit->SetVal(TMath::Max(fLastOrbit->GetVal(), (Long_t) timeStamp));
}

void AliTriggerTask::Terminate(Option_t *)
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

  fOutput = dynamic_cast<TList*> (GetOutputData(0));
  if (!fOutput)
    Printf("ERROR: fOutput not available");
    
  //fOutput->Print();

  if (fOutput)
  {
    for (Int_t j=0; j<fNTriggerClasses; j++)
      for (Int_t i=0; i<fNTriggers; i++)
        fStats[j][i] = dynamic_cast<TH1*> (fOutput->FindObject(Form("fStats_%d_%d", j, i)));
    fPhysicsSelection = dynamic_cast<AliPhysicsSelection*> (fOutput->FindObject("AliPhysicsSelection"));
  }

  TFile* fout = new TFile("trigger.root", "RECREATE");

  for (Int_t j=0; j<fNTriggerClasses; j++)
    for (Int_t i=0; i<fNTriggers; i++)
      if (fStats[j][i])
        fStats[j][i]->Write();
        
  if (fPhysicsSelection)
  {
    fPhysicsSelection->SaveHistograms("physics_selection");
    fPhysicsSelection->Print();
  }
    
  if (fFirstOrbit)
    fFirstOrbit->Dump();
  if (fLastOrbit)
    fLastOrbit->Dump();
    
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

  Int_t triggerClass = 0;

  Int_t base = 1;
  if (fStats[triggerClass][0])
    base = (Int_t) fStats[triggerClass][0]->Integral();

  Int_t length = fEndTime - fStartTime;
  
  for (Int_t i=0; i<fNTriggers; i++)
    if (fStats[triggerClass][i])
    {
      c->cd(i+1);
      fStats[triggerClass][i]->Draw();
      Printf("%s: %d triggers | %f %% of all triggered | Rate: %f Hz", AliTriggerAnalysis::GetTriggerName(fTriggerList[i]), (UInt_t) fStats[triggerClass][i]->Integral(), 100.0 * fStats[triggerClass][i]->Integral() / base, (length > 0) ? (fStats[triggerClass][i]->Integral() / length) : -1);
    }

  Printf("Writting result to trigger.root");
}
