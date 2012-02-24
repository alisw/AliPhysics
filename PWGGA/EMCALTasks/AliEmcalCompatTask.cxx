// $Id$
//
// Task to setup emcal related global objects
//
//

#include "AliEmcalCompatTask.h"
#include <TClonesArray.h>
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliEsdTrackExt.h"
#include "AliEventplane.h"
#include "AliCentrality.h"
ClassImp(AliEmcalCompatTask)

//________________________________________________________________________
AliEmcalCompatTask::AliEmcalCompatTask() : 
  AliAnalysisTaskSE(),
  fDoCent(1),
  fDoEp(1)
{
  // Constructor.
}

//________________________________________________________________________
AliEmcalCompatTask::AliEmcalCompatTask(const char *name) : 
  AliAnalysisTaskSE(name),
  fDoCent(1),
  fDoEp(1)
{
  // Constructor.

  fBranchNames = "ESD:AliESDHeader.,AliESDRun.,Tracks";
}

//________________________________________________________________________
AliEmcalCompatTask::~AliEmcalCompatTask()
{
  // Destructor.
}

//________________________________________________________________________
void AliEmcalCompatTask::UserExec(Option_t *) 
{
  // Main loop, called for each event.

  AliESDEvent *esdEv = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!esdEv) {
    AliError("Task works only on ESD events, returning");
    return;
  }

  AliAnalysisManager *am = AliAnalysisManager::GetAnalysisManager();
  if (!am) {
    AliError("Manager zero, returning");
    return;
  }

  am->LoadBranch("AliESDHeader.");
  am->LoadBranch("AliESDRun.");

  AliESDHeader *header = esdEv->GetHeader();
  TString title;
  if (header)
    title = header->GetTitle();
  if (title.Length()==0)
    return;

  if (fDoCent) {
    am->LoadBranch("Centrality.");
    AliCentrality *centin = dynamic_cast<AliCentrality*>(esdEv->FindListObject("Centrality"));
    AliCentrality *centout = esdEv->GetCentrality();
    if (centin&&centout&&centout->GetQuality()==999) {
      centout->SetQuality(centin->GetQuality());
      centout->SetCentralityV0M(centin->GetCentralityPercentileUnchecked("V0M"));
      centout->SetCentralityFMD(centin->GetCentralityPercentileUnchecked("FMD"));
      centout->SetCentralityTRK(centin->GetCentralityPercentileUnchecked("TRK"));
      centout->SetCentralityTKL(centin->GetCentralityPercentileUnchecked("TKL"));
      centout->SetCentralityCL0(centin->GetCentralityPercentileUnchecked("CL0"));
      centout->SetCentralityCL1(centin->GetCentralityPercentileUnchecked("CL1"));
      centout->SetCentralityV0MvsFMD(centin->GetCentralityPercentileUnchecked("V0MvsFMD"));
      centout->SetCentralityTKLvsV0M(centin->GetCentralityPercentileUnchecked("TKLvsV0M"));
      centout->SetCentralityZEMvsZDC(centin->GetCentralityPercentileUnchecked("ZEMvsZDC"));
    }
  }

  if (fDoEp) {
    am->LoadBranch("Eventplane.");
    AliEventplane *epin  = dynamic_cast<AliEventplane*>(esdEv->FindListObject("Eventplane"));
    AliEventplane *epout = esdEv->GetEventplane();
    if (epin&&epout&&(epout->GetQVector()==0)&&(epin->GetQVector()!=0)) {
      epout->SetQVector(new TVector2(*epin->GetQVector()));
      epout->SetEventplaneQ(epin->GetEventplane("Q"));
      epout->SetQsub(new TVector2(*epin->GetQsub1()),new TVector2(*epin->GetQsub2()));
      epout->SetQsubRes(epin->GetQsubRes());
    }
  }

  TTree *tree = am->GetTree();
  if (tree&&tree->GetBranch("PicoTracks"))
    am->LoadBranch("PicoTracks");

  if (tree&&tree->GetBranch("Tracks")) {
    am->LoadBranch("Tracks");
    TClonesArray *ts = dynamic_cast<TClonesArray*>(esdEv->FindListObject("Tracks"));
    if (ts) {
      TString clsname(ts->GetClass()->GetName());
      if (clsname == "AliEsdTrackExt") {
        const Int_t N = ts->GetEntries();
        for (Int_t i=0; i<N; ++i) {
          AliEsdTrackExt *t = static_cast<AliEsdTrackExt*>(ts->At(i));
          if (t) {
            t->SetESDEvent(esdEv);
            t->Setup();
          }
        }
      }
    }
  }
}
