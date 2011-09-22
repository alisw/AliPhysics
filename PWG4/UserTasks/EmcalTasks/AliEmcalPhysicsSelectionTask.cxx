/* $Id$ */

#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TROOT.h>
#include "AliAnalysisDataContainer.h"
#include "AliAnalysisDataSlot.h"
#include "AliAnalysisManager.h"
#include "AliEmcalPhysicsSelection.h"
#include "AliEmcalPhysicsSelectionTask.h"
#include "AliESDEvent.h"
#include "AliInputEventHandler.h"
#include "AliLog.h"

ClassImp(AliEmcalPhysicsSelectionTask)

//__________________________________________________________________________________________________
AliEmcalPhysicsSelectionTask::AliEmcalPhysicsSelectionTask() :
  AliPhysicsSelectionTask(),
  fDoWriteHistos(1),
  fNCalled(0),
  fNAccepted(0),
  fHAcc(0),
  fHEvtTypes(0)
{
  // Default constructor.
}

//__________________________________________________________________________________________________
AliEmcalPhysicsSelectionTask::AliEmcalPhysicsSelectionTask(const char* opt) : 
  AliPhysicsSelectionTask(),
  fDoWriteHistos(1),
  fNCalled(0),
  fNAccepted(0),
  fHAcc(0),
  fHEvtTypes(0)
{
  // Constructor.

  fOption = opt;
  fPhysicsSelection = new AliEmcalPhysicsSelection;

  AliInputEventHandler* handler = dynamic_cast<AliInputEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if (handler) {
    handler->SetEventSelection(fPhysicsSelection);
    AliInfo("Physics Event Selection enabled.");
  } else {
    AliError("No input event handler connected to analysis manager. No Physics Event Selection.");
  }
  // Define input and output slots here
  DefineOutput(1, TList::Class());
  fBranchNames = "ESD:AliESDRun.,AliESDHeader.,AliMultiplicity.,AliESDVZERO.,"
                 "AliESDZDC.,SPDVertex.,PrimaryVertex.,TPCVertex.,Tracks,SPDPileupVertices";
  
  AliLog::SetClassDebugLevel("AliEmcalPhysicsSelectionTask", AliLog::kWarning);
}

//__________________________________________________________________________________________________
void AliEmcalPhysicsSelectionTask::UserCreateOutputObjects()
{
  // User create outputs.

  AliPhysicsSelectionTask::UserCreateOutputObjects();
  fHAcc = new TH1D("hEvCount",";0=rej/1=acc;#",2,-0.5,1.5);
  fOutput->Add(fHAcc);
  fHEvtTypes = new TH1D("hEvtTypes",";#",8,-0.5,7.5);
  fHEvtTypes->GetXaxis()->SetBinLabel(1,"All");
  fHEvtTypes->GetXaxis()->SetBinLabel(2,"MB");
  fHEvtTypes->GetXaxis()->SetBinLabel(3,"FO");
  fHEvtTypes->GetXaxis()->SetBinLabel(4,"EMC");
  fHEvtTypes->GetXaxis()->SetBinLabel(5,"Good");
  fHEvtTypes->GetXaxis()->SetBinLabel(6,"HC");
  fHEvtTypes->GetXaxis()->SetBinLabel(7,"HT");
  fHEvtTypes->GetXaxis()->SetBinLabel(8,"LED");
  fOutput->Add(fHEvtTypes);
  if (!fDoWriteHistos) {
    fOutput->Remove(fPhysicsSelection);
  }
}

//__________________________________________________________________________________________________
void AliEmcalPhysicsSelectionTask::UserExec(const Option_t *opt)
{
  // User exec.

  AliPhysicsSelectionTask::UserExec(opt);

  ++fNCalled;

  UInt_t res = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  if (res>0) {
    ++fNAccepted;
    fHAcc->Fill(1);
  } else {
    fHAcc->Fill(0);
  }

  AliEmcalPhysicsSelection *ps=static_cast<AliEmcalPhysicsSelection *>(fPhysicsSelection);
  fHEvtTypes->Fill(0);
  if (res&AliVEvent::kAnyINT)
    fHEvtTypes->Fill(1);
  if (ps->IsFastOnly())
    fHEvtTypes->Fill(2);
  if ((res&AliVEvent::kEMC1) || (res&AliVEvent::kEMC7))
    fHEvtTypes->Fill(3);
  if (ps->IsGoodEvent())
    fHEvtTypes->Fill(4);
  if (res&AliEmcalPhysicsSelection::kEmcalHC)
    fHEvtTypes->Fill(5);
  if (res&AliEmcalPhysicsSelection::kEmcalHT)
    fHEvtTypes->Fill(6);
  if (ps->IsLedEvent())
    fHEvtTypes->Fill(7);
}

//__________________________________________________________________________________________________
void AliEmcalPhysicsSelectionTask::Terminate(Option_t *)
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

  AliInfo(Form("Called %d times, accepted %d events", fNCalled, fNAccepted));

  if (!fDoWriteHistos)
    return;

  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutput) {
    AliError("fOutput not available");
    return;
  }

  AliAnalysisDataSlot *oslot = GetOutputSlot(1);
  if (!oslot)
    return;

  AliAnalysisDataContainer *ocont = oslot->GetContainer();
  if (!ocont)
    return;

  TFile *file = OpenFile(1);
  if (!file)
    return;

  TDirectory::TContext context(file); 
  if (AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
    fPhysicsSelection = dynamic_cast<AliPhysicsSelection*> (fOutput->FindObject("AliPhysicsSelection"));
  }
  if (fPhysicsSelection) {
    //fPhysicsSelection->Print();
    fPhysicsSelection->SaveHistograms(Form("%sHists",ocont->GetName()));
    AliInfo(Form("Writing result to %s",file->GetName()));
  }
  fOutput->Remove(fPhysicsSelection);
}
