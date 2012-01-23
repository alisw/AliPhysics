
/* $Id$ */

//===================================================================
//  Class AliAnalysisTaskLUT
//
//  Extract ESD muon tracks information and store in ntuple.
//  Used for offline calculation/adjustment of Look-up-Tables for the
//  trigger chambers.
//
//  Clermont-Ferrand 2008
//===================================================================

#define AliAnalysisTaskLUT_cxx

#include "Riostream.h"
#include "TChain.h"
#include "TTree.h"
#include "TNtuple.h"

#include "AliESDEvent.h"
#include "AliESDVertex.h"
#include "AliESDMuonTrack.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisTaskLUT.h"

ClassImp(AliAnalysisTaskLUT)

//________________________________________________________________________
AliAnalysisTaskLUT::AliAnalysisTaskLUT(const char *name) :
  AliAnalysisTask(name,""), 
  fESDEvent(0), 
  fTracksLUT(0)
{
  // Constructor.
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());

  // Output slot #0 writes into a TNtuple container
  DefineOutput(0, TNtuple::Class());

}

//___________________________________________________________________________
void AliAnalysisTaskLUT::ConnectInputData(Option_t *) 
{
  // Input ESD tree
  printf("   ConnectInputData of task %s\n", GetName());

  TTree *tinput = (TTree*)GetInputData(0);
  fESDEvent = new AliESDEvent();
  fESDEvent->ReadFromTree(tinput);

}

//___________________________________________________________________________
void AliAnalysisTaskLUT::CreateOutputObjects() 
{
  // Output object (TNtuple)
  printf("   CreateOutputObjects of task %s\n", GetName());
  
  fTracksLUT = new TNtuple("ntTracksLUT","ntTracksLUT","TrigMask:VertZ:NHit:Circ:StripX:StripY:Dev:Lpt:Hpt:Chi2:Chi2trg:Z:P:Pt:Phi:The:TrackZ:PU:PtU:PhiU:TheU:TrackZU");

}

//________________________________________________________________________
void AliAnalysisTaskLUT::Exec(Option_t *) 
{
  // Execution
  Float_t ntvar[22];

  Float_t px, py, pz, pt;
  Float_t zVertex;
  Int_t dev;

  if (!fESDEvent) return;
  
  Int_t nTracks = fESDEvent->GetNumberOfMuonTracks();

  AliESDVertex* vertex = (AliESDVertex*) fESDEvent->GetVertex();
  zVertex = 0.0;
  if (vertex) {
    zVertex = vertex->GetZv();
  }
  ULong64_t mask = fESDEvent->GetTriggerMask();

  ntvar[ 0] = mask;
  ntvar[ 1] = zVertex;

  for(Int_t iTracks = 0; iTracks < nTracks; iTracks++) {

    AliESDMuonTrack* muonTrack = fESDEvent->GetMuonTrack(iTracks);
    if (muonTrack == 0x0) continue;

    if (!muonTrack->GetMatchTrigger()) continue;

    dev  = muonTrack->LoDev();
    dev *= (Int_t)(TMath::Sign(1.0,muonTrack->GetInverseBendingMomentum()));
    dev += 15;

    ntvar[ 2] = muonTrack->GetNHit();
    ntvar[ 3] = muonTrack->LoCircuit();
    ntvar[ 4] = muonTrack->LoStripX();
    ntvar[ 5] = muonTrack->LoStripY();
    ntvar[ 6] = dev;
    ntvar[ 7] = muonTrack->LoLpt();
    ntvar[ 8] = muonTrack->LoHpt();
    ntvar[ 9] = muonTrack->GetChi2() / (2.0 * muonTrack->GetNHit() - 5);
    ntvar[10] = muonTrack->GetChi2MatchTrigger();
    ntvar[11] = TMath::Sign(1.,muonTrack->GetInverseBendingMomentum());

    // corrected to vertex

    px = muonTrack->Px();
    py = muonTrack->Py();
    pz = muonTrack->Pz();
    pt = TMath::Sqrt(px*px+py*py);
    ntvar[12] = muonTrack->P();
    ntvar[13] = TMath::Sqrt(px*px+py*py);
    ntvar[14] = TMath::ATan2(py,px)*180./TMath::Pi();
    ntvar[15] = TMath::ATan2(pt,pz)*180./TMath::Pi();
    ntvar[16] = muonTrack->GetZ();

    // uncorrected

    px = muonTrack->PxUncorrected();
    py = muonTrack->PyUncorrected();
    pz = muonTrack->PzUncorrected();
    pt = TMath::Sqrt(px*px+py*py);
    ntvar[17] = muonTrack->PUncorrected();
    ntvar[18] = TMath::Sqrt(px*px+py*py);
    ntvar[19] = TMath::ATan2(py,px)*180./TMath::Pi();
    ntvar[20] = TMath::ATan2(pt,pz)*180./TMath::Pi();
    ntvar[21] = muonTrack->GetZUncorrected();

    fTracksLUT->Fill(ntvar);

  } // end ESD track loop

  // Post final data. It will be written to a file with option "RECREATE"
  PostData(0, fTracksLUT);
  
}      

//________________________________________________________________________
void AliAnalysisTaskLUT::Terminate(Option_t *) {
  // End function, empty
}

