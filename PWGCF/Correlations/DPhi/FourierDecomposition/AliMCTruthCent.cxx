//
// Determine MC truth centrality for kinematics-only productions
// pass it via centrality object
//

#include "AliMCTruthCent.h"

#include <TString.h>
#include <TH1.h>
#include <TH2.h>
#include <TChain.h>
#include <TList.h>

#include "AliAnalysisManager.h"
#include "AliVEventHandler.h"
#include "AliMCEvent.h"
#include "AliVParticle.h"
#include "AliLog.h"
#include "AliESDEvent.h"
#include "AliCentrality.h"

ClassImp(AliMCTruthCent)

//________________________________________________________________________
AliMCTruthCent::AliMCTruthCent() :
AliAnalysisTaskSE("AliMCTruthCent"),
fOutputList(0x0),
fFillHistos(kFALSE),
fHMultV0A(0x0), fHMultV0C(0x0), fHMultV0M(0x0), fHMultV0AvsV0C(0x0),
fHCentV0A(0x0), fHCentV0C(0x0), fHCentV0M(0x0), fHCentV0AvsV0C(0x0),
fV0ALo(0.0), fV0AHi(100.0),
fV0CLo(0.0), fV0CHi(100.0),
fV0MLo(0.0), fV0MHi(200.0)
{
  // Constructor.
}

//________________________________________________________________________
AliMCTruthCent::AliMCTruthCent(const char *name) :
AliAnalysisTaskSE(name),
fOutputList(0x0),
fFillHistos(kFALSE),
fHMultV0A(0x0), fHMultV0C(0x0), fHMultV0M(0x0), fHMultV0AvsV0C(0x0),
fHCentV0A(0x0), fHCentV0C(0x0), fHCentV0M(0x0), fHCentV0AvsV0C(0x0),
fV0ALo(0.0), fV0AHi(100.0),
fV0CLo(0.0), fV0CHi(100.0),
fV0MLo(0.0), fV0MHi(200.0)
{
  // Constructor.
  
  // Define input slot
  DefineInput(0, TChain::Class());

}

//________________________________________________________________________
AliMCTruthCent::~AliMCTruthCent()
{
  // Destructor.
  if (fOutputList)
    delete fOutputList;
}

//________________________________________________________________________
void AliMCTruthCent::UserCreateOutputObjects()
{
  // Create my user objects.
  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    AliFatal("No analysis manager!");
    return;
  }
  
  AliVEventHandler *handler = mgr->GetInputEventHandler();
  if (!handler) {
    AliFatal("No input handler!");
    return;
  }
  
  if (!handler->InheritsFrom("AliESDInputHandler")) {
    AliError("Sorry, this only works for kinematics trees ... return");
    return;
  }
  
  if (fFillHistos) {
    fOutputList = new TList();
    fOutputList->SetOwner();
    
    // histograms
    fHMultV0A      = new TH1D("fHMultV0A","V0A multiplicity",500,0.0,500.0);
    fHMultV0C      = new TH1D("fHMultV0C","V0C multiplicity",500,0.0,500.0);
    fHMultV0M      = new TH1D("fHMultV0M","V0M multiplicity",500,0.0,500.0);
    fHMultV0AvsV0C = new TH2D("fHMultV0AvsV0C","V0A vs. C multiplicity",500,0.0,500.0,500,0.0,500.0);
    
    fOutputList->Add(fHMultV0A);
    fOutputList->Add(fHMultV0C);
    fOutputList->Add(fHMultV0M);
    fOutputList->Add(fHMultV0AvsV0C);
    
    fHCentV0A      = new TH1D("fHCentV0A","V0A centrality",100,0.0,100.0);
    fHCentV0C      = new TH1D("fHCentV0C","V0C centrality",100,0.0,100.0);
    fHCentV0M      = new TH1D("fHCentV0M","V0M centrality",100,0.0,100.0);
    fHCentV0AvsV0C = new TH2D("fHCentV0AvsV0C","V0A vs. C centrality",100,0.0,100.0,100,0.0,100.0);
    
    fOutputList->Add(fHCentV0A);
    fOutputList->Add(fHCentV0C);
    fOutputList->Add(fHCentV0M);
    fOutputList->Add(fHCentV0AvsV0C);
    
    PostData(1, fOutputList);
  }
}

//________________________________________________________________________
void AliMCTruthCent::UserExec(Option_t *)
{
  // Main loop, called for each event.
  if (!InputEvent()) {
    AliError("Could not retrieve event! Returning");
    return;
  }
  
  if (!MCEvent()) {
    AliError("Could not retrieve ESD MC event! Returning");
    return;
  }
  
  // can I get the centrality object from the InputEvent?
  // then use it and overwrite the centrality
  AliCentrality *esdCent = InputEvent()->GetCentrality();
  esdCent->SetQuality(0);
  
  // define experiment-like estimators: multiplicity in range of
  // V0A  2.8 < eta <  5.1
  // V0C -3.7 < eta < -1.7
  // V0M: sum of V0A and V0C
  Double_t etaV0ALo = 2.8;
  Double_t etaV0AHi = 5.1;
  Double_t etaV0CLo = -3.7;
  Double_t etaV0CHi = -1.7;
  
  Int_t nV0A = 0;
  Int_t nV0C = 0;
  
  // loop over MC tracks
  Double_t       etaTrack = 0.0;
  AliVParticle * track    = 0x0;
  const Int_t Ntracks = MCEvent()->GetNumberOfTracks();
  for (Int_t iTracks = 0; iTracks < Ntracks; ++iTracks) {
    track = GetTrack(iTracks);
    if (!track)
      continue;
    
    if (track->Charge()==0)
      continue;

    etaTrack = track->Eta();
    
    if (etaTrack>etaV0ALo && etaTrack<etaV0AHi) {
      ++nV0A;
    } else if (etaTrack>etaV0CLo && etaTrack<etaV0CHi) {
      ++nV0C;
    }
  }
  Int_t nV0M = nV0A+nV0C;

  
  // linear centrality approximation
  Double_t centV0A = (fV0AHi-nV0A) / (fV0AHi-fV0ALo);
  centV0A *= 100.0;
  if (centV0A >= 100.0)
    centV0A = 99.9;
  if (centV0A < 0.0)
    centV0A = 0.0;
  esdCent->SetCentralityV0A(centV0A);

  Double_t centV0C = (fV0CHi-nV0C) / (fV0CHi-fV0CLo);
  centV0C *= 100.0;
  if (centV0C >= 100.0)
    centV0C = 99.9;
  if (centV0C < 0.0)
    centV0C = 0.0;
  esdCent->SetCentralityV0C(centV0C);
  
  Double_t centV0M = (fV0MHi-nV0M) / (fV0MHi-fV0MLo);
  centV0M *= 100.0;
  if (centV0M >= 100.0)
    centV0M = 99.9;
  if (centV0M < 0.0)
    centV0M = 0.0;
  esdCent->SetCentralityV0M(centV0M);

  if (fFillHistos) {
    // write out distributions to histograms
    fHMultV0A->Fill(nV0A);
    fHMultV0C->Fill(nV0C);
    fHMultV0M->Fill(nV0M);
    fHMultV0AvsV0C->Fill(nV0A,nV0C);
    
    fHCentV0A->Fill(centV0A);
    fHCentV0C->Fill(centV0C);
    fHCentV0M->Fill(centV0M);
    fHCentV0AvsV0C->Fill(centV0A,centV0C);
    
    PostData(1, fOutputList);
  }

  // next step: get centrality from histogram...
  
}

//________________________________________________________________________
AliVParticle* AliMCTruthCent::GetTrack(Int_t i)
{
  if (!MCEvent()->IsPhysicalPrimary(i))
    return 0;
  
  return MCEvent()->GetTrack(i);
}

