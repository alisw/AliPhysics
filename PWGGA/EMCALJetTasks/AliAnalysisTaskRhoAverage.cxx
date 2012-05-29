// $Id$
//
// Calculation of rho, method: sum of all particle pt / full acceptance area
//
// Authors: Salvatore Aiola

#include <TList.h>
#include <TClonesArray.h>
#include <TMath.h>

#include "AliLog.h"
#include "AliAnalysisManager.h"
#include "AliVEventHandler.h"
#include "AliCentrality.h"
#include "AliEmcalJet.h"
#include "AliVTrack.h"
#include "AliVCluster.h"

#include "AliAnalysisTaskRhoAverage.h"

ClassImp(AliAnalysisTaskRhoAverage)

//________________________________________________________________________
AliAnalysisTaskRhoAverage::AliAnalysisTaskRhoAverage() : 
  AliAnalysisTaskRhoBase(),
  fTracksName("tracks"),
  fClustersName("caloClusters"),
  fJetsName("KtJets"),
  fEtaMin(-0.9),
  fEtaMax(0.9),
  fPhiMin(0),
  fPhiMax(2 * TMath::Pi()),
  fPtMin(0.15)
{
  // Constructor
}

//________________________________________________________________________
AliAnalysisTaskRhoAverage::AliAnalysisTaskRhoAverage(const char *name) :
  AliAnalysisTaskRhoBase(name),
  fTracksName("tracks"),
  fClustersName("caloClusters"),
  fJetsName("KtJets"),
  fEtaMin(-0.9),
  fEtaMax(0.9),
  fPhiMin(0),
  fPhiMax(2 * TMath::Pi()),
  fPtMin(0.15)
{
  // Constructor

}

//________________________________________________________________________
void AliAnalysisTaskRhoAverage::UserExec(Option_t *) 
{
  // Main loop, called for each event.
  
  AliAnalysisTaskRhoBase::UserExec("");

  fRho->SetVal(-1);

  TClonesArray *jets = 0;
  TClonesArray *tracks = 0;
  TClonesArray *clusters = 0;

  tracks = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fTracksName));
  if (!tracks) {
  AliError(Form("Pointer to tracks %s == 0", fTracksName.Data() ));
   return;
  }

  if (fClustersName != "") {
    clusters = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fClustersName));
    if (!clusters) {
      AliError(Form("Pointer to clusters %s == 0", fClustersName.Data() ));
      return;
    }
  }

  if (fJetsName != "") {
    jets = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fJetsName));
    if (!jets) {
      AliError(Form("Pointer to jets %s == 0", fJetsName.Data() ));
      return;
    }
  }

  Double_t rho = 0;
  
  Int_t Ntracks = tracks->GetEntries();

  Int_t Nclusters = 0;
  if (clusters)
    Nclusters = clusters->GetEntries();

  Int_t Njets = 0;
  if (jets)
    Njets = jets->GetEntries();

  Float_t maxJetPt = 0;
  Int_t maxJetId = -1;
  AliEmcalJet *maxJet = 0;
  for (Int_t ij = 0; ij < Njets; ij++) {
      
    AliEmcalJet *jet = static_cast<AliEmcalJet*>(jets->At(ij));
  
    if (!jet) {
      AliError(Form("Could not receive jet %d", ij));
      continue;
    } 
  
    if (jet->Pt() > maxJetPt) {
      maxJetPt = jet->Pt();
      maxJetId = ij;
    }
  }

  if (maxJetId >= 0)
    maxJet = static_cast<AliEmcalJet*>(jets->At(maxJetId));

  for (Int_t it = 0; it < Ntracks; it++) {
      
    AliVTrack *track = static_cast<AliVTrack*>(tracks->At(it));
  
    if (!track) {
      AliError(Form("Could not receive track %d", it));
      continue;
    } 

    if (track->Eta() < fEtaMin || track->Eta() > fEtaMax || track->Phi() < fPhiMin || track->Phi() > fPhiMax)
      continue;

    if (track->Pt() < fPtMin)
      continue;

    if (maxJet && IsJetTrack(maxJet, it))
      continue;

    rho += track->Pt();
  }
  
  Double_t vertex[] = {0, 0, 0};
  InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);

  for (Int_t ic = 0; ic < Nclusters; ic++) {
      
    AliVCluster *cluster = static_cast<AliVCluster*>(clusters->At(ic));
  
    if (!cluster) {
      AliError(Form("Could not receive cluster %d", ic));
      continue;
    } 

    Float_t pos[3];
    cluster->GetPosition(pos);
    TVector3 clusVec(pos);

    if (clusVec.Eta() < fEtaMin || clusVec.Eta() > fEtaMax || clusVec.Phi() < fPhiMin || clusVec.Phi() > fPhiMax)
      continue;

    TLorentzVector nPart;
    cluster->GetMomentum(nPart, const_cast<Double_t*>(vertex));

    if (nPart.Et() < fPtMin)
      continue;

    if (maxJet && IsJetCluster(maxJet, ic))
      continue;

    rho += nPart.Et();
  }
 
  Double_t area = (fEtaMax - fEtaMin) * (fPhiMax - fPhiMin);

  if (maxJet)
    area -= maxJet->Area();

  rho /= area;

  fRho->SetVal(rho);
}      

//________________________________________________________________________
Bool_t AliAnalysisTaskRhoAverage::IsJetTrack(AliEmcalJet* jet, Int_t itrack) const
{
  for (Int_t i = 0; i < jet->GetNumberOfTracks(); i++) {
    Int_t ijettrack = jet->TrackAt(i);
    if (ijettrack == itrack)
      return kTRUE;
  }
  return kFALSE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskRhoAverage::IsJetCluster(AliEmcalJet* jet, Int_t iclus) const
{
  for (Int_t i = 0; i < jet->GetNumberOfClusters(); i++) {
    Int_t ijetclus = jet->ClusterAt(i);
    if (ijetclus == iclus)
      return kTRUE;
  }
  return kFALSE;
}


//________________________________________________________________________
void AliAnalysisTaskRhoAverage::Terminate(Option_t *) 
{

}
