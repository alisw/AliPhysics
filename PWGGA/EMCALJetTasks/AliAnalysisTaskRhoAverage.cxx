// $Id$
//
// Calculation of rho, method: sum of all particle pt / full acceptance area.
//
// Authors: S. Aiola

#include "AliAnalysisTaskRhoAverage.h"

#include <TClonesArray.h>
#include <TList.h>
#include <TMath.h>

#include "AliAnalysisManager.h"
#include "AliCentrality.h"
#include "AliEmcalJet.h"
#include "AliLog.h"
#include "AliRhoParameter.h"
#include "AliVCluster.h"
#include "AliVEventHandler.h"
#include "AliVTrack.h"

ClassImp(AliAnalysisTaskRhoAverage)

//________________________________________________________________________
AliAnalysisTaskRhoAverage::AliAnalysisTaskRhoAverage() : 
  AliAnalysisTaskRhoBase(),
  fTracksName(),
  fClustersName(),
  fJetsName(),
  fEtaMin(0),
  fEtaMax(0),
  fPhiMin(0),
  fPhiMax(0),
  fPtMin(0),
  fClusters(0),
  fJets(0),
  fTracks(0)
{
  // Default constructor.
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
  fPtMin(0.15),
  fClusters(0),
  fJets(0),
  fTracks(0)
{
  // Constructor.
}

//________________________________________________________________________
void AliAnalysisTaskRhoAverage::UserExec(Option_t *) 
{
  // Main loop, called for each event.
  
  if (!fIsInit) {
    ExecOnce();
    fIsInit = 1;
  }

  fRho->SetVal(-1);

  Double_t rho = 0;
  
  Int_t Ntracks = 0;
  if (fTracks) 
    Ntracks = fTracks->GetEntriesFast();

  Int_t Nclusters = 0;
  if (fClusters)
    Nclusters = fClusters->GetEntriesFast();

  Int_t Njets = 0;
  if (fJets)
    Njets = fJets->GetEntriesFast();

  Double_t maxJetPt = 0;
  Int_t maxJetId = -1;
  AliEmcalJet *maxJet = 0;
  for (Int_t ij = 0; ij < Njets; ij++) {
      
    AliEmcalJet *jet = static_cast<AliEmcalJet*>(fJets->At(ij));
  
    if (!jet) {
      AliError(Form("%s: Could not receive jet %d", GetName(), ij));
      continue;
    } 
  
    if (jet->Pt() > maxJetPt) {
      maxJetPt = jet->Pt();
      maxJetId = ij;
    }
  }

  if (maxJetId >= 0)
    maxJet = static_cast<AliEmcalJet*>(fJets->At(maxJetId));

  for (Int_t it = 0; it < Ntracks; ++it) {
      
    AliVTrack *track = static_cast<AliVTrack*>(fTracks->At(it));
  
    if (!track) {
      AliError(Form("%s: Could not receive track %d", GetName(), it));
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

  for (Int_t ic = 0; ic < Nclusters; ++ic) {
      
    AliVCluster *cluster = static_cast<AliVCluster*>(fClusters->At(ic));
  
    if (!cluster) {
      AliError(Form("%s: Could not receive cluster %d", GetName(), ic));
      continue;
    } 

    Float_t pos[3];
    cluster->GetPosition(pos);
    TVector3 clusVec(pos);

    if (clusVec.Eta() < fEtaMin || clusVec.Eta() > fEtaMax || 
        clusVec.Phi() < fPhiMin || clusVec.Phi() > fPhiMax)
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

  if (area>0) {
    rho /= area;
    fRho->SetVal(rho);
  } else {
    AliError(Form("%s: Area negative %f", GetName(), area));
  }
}      

//________________________________________________________________________
void AliAnalysisTaskRhoAverage::ExecOnce() 
{
  // Initialize some settings that need to be determined in UserExec.

  AliAnalysisTaskRhoBase::ExecOnce();

  fClusters = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fClustersName));
  if (!fClusters) {
    AliError(Form("%s: Pointer to jets %s == 0", GetName(), fClustersName.Data() ));
    return;
  }

  fJets = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fJetsName));
  if (!fJets) {
    AliError(Form("%s: Pointer to jets %s == 0", GetName(), fJetsName.Data() ));
    return;
  }

  fTracks = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fTracksName));
  if (!fTracks) {
    AliError(Form("%s: Pointer to tracks %s == 0", GetName(), fTracksName.Data() ));
    return;
  }
}

//________________________________________________________________________
Bool_t AliAnalysisTaskRhoAverage::IsJetTrack(AliEmcalJet* jet, Int_t itrack) const
{
  // Return true if track is in jet.

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
  // Return true if cluster is in jet.

  for (Int_t i = 0; i < jet->GetNumberOfClusters(); i++) {
    Int_t ijetclus = jet->ClusterAt(i);
    if (ijetclus == iclus)
      return kTRUE;
  }
  return kFALSE;
}
