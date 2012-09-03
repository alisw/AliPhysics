// $Id$
//
// Calculation of rho, method: median all particle pt / multiplicity density.
//
// Authors: S. Aiola

#include "AliAnalysisTaskRhoAverage.h"

#include <TClonesArray.h>
#include <TMath.h>

#include "AliLog.h"
#include "AliRhoParameter.h"
#include "AliVCluster.h"
#include "AliVTrack.h"

ClassImp(AliAnalysisTaskRhoAverage)

//________________________________________________________________________
AliAnalysisTaskRhoAverage::AliAnalysisTaskRhoAverage() : 
  AliAnalysisTaskRhoBase("AliAnalysisTaskRhoAverage"),
  fRhoType(0),
  fNExclLeadPart(0),
  fUseMedian(kFALSE)
{
  // Default constructor.
}

//________________________________________________________________________
AliAnalysisTaskRhoAverage::AliAnalysisTaskRhoAverage(const char *name, Bool_t histo) :
  AliAnalysisTaskRhoBase(name, histo),
  fRhoType(0),
  fNExclLeadPart(0),
  fUseMedian(kFALSE)
{
  // Constructor.
}

//________________________________________________________________________
Bool_t AliAnalysisTaskRhoAverage::Run() 
{
  // Run the analysis.

  static Double_t rhovec[9999];
  Int_t NpartAcc = 0;

  Int_t   maxPartIds[] = {0, 0};
  Float_t maxPartPts[] = {0, 0};

  // push all jets within selected acceptance into stack

  if (fNExclLeadPart > 0) {

    if (fTracks && (fRhoType == 0 || fRhoType == 1)) {
      
      const Int_t Ntracks = fTracks->GetEntriesFast();
      
      for (Int_t it = 0; it < Ntracks; ++it) {
	
	AliVTrack *track = static_cast<AliVTrack*>(fTracks->At(it));
	
	if (!track) {
	  AliError(Form("%s: Could not receive track %d", GetName(), it));
	  continue;
	} 
	
	if (!AcceptTrack(track))
	  continue;
	
	if (track->Pt() > maxPartPts[0]) {
	  maxPartPts[1] = maxPartPts[0];
	  maxPartIds[1] = maxPartIds[0];
	  maxPartPts[0] = track->Pt();
	  maxPartIds[0] = it+1;
	} 
	else if (track->Pt() > maxPartPts[1]) {
	  maxPartPts[1] = track->Pt();
	  maxPartIds[1] = it+1;
	}
      }
    }

    if (fCaloClusters && (fRhoType == 0 || fRhoType == 2)) {

      const Int_t Nclusters = fCaloClusters->GetEntriesFast();
      
      for (Int_t ic = 0; ic < Nclusters; ++ic) {
	
	AliVCluster *cluster = static_cast<AliVCluster*>(fCaloClusters->At(ic));
	
	if (!cluster) {
	  AliError(Form("%s: Could not receive cluster %d", GetName(), ic));
	  continue;
	} 
	
	if (!AcceptCluster(cluster))
	  continue;
	
	TLorentzVector nPart;
	cluster->GetMomentum(nPart, fVertex);
	
	if (nPart.Pt() > maxPartPts[0]) {
	  maxPartPts[1] = maxPartPts[0];
	  maxPartIds[1] = maxPartIds[0];
	  maxPartPts[0] = nPart.Pt();
	  maxPartIds[0] = -ic-1;
	} 
	else if (nPart.Pt() > maxPartPts[1]) {
	  maxPartPts[1] = nPart.Pt();
	  maxPartIds[1] = -ic-1;
	}
      }
    }
 
    if (fNExclLeadPart < 2) {
      maxPartIds[1] = 0;
      maxPartPts[1] = 0;
    }
  }
  
  if (fTracks && (fRhoType == 0 || fRhoType == 1)) {

    const Int_t Ntracks = fTracks->GetEntriesFast();

    for (Int_t it = 0; it < Ntracks; ++it) {

      // exlcuding lead particles
      if (it == maxPartIds[0]-1 || it == maxPartIds[1]-1)
	continue;
      
      AliVTrack *track = static_cast<AliVTrack*>(fTracks->At(it));
  
      if (!track) {
	AliError(Form("%s: Could not receive track %d", GetName(), it));
	continue;
      } 

      if (!AcceptTrack(track))
	continue;

      rhovec[NpartAcc] = track->Pt();
      ++NpartAcc;
    }
  }

  if (fCaloClusters && (fRhoType == 0 || fRhoType == 2)) {

    const Int_t Nclusters = fCaloClusters->GetEntriesFast();

    for (Int_t ic = 0; ic < Nclusters; ++ic) {

      // exlcuding lead particles
      if (ic == -maxPartIds[0]-1 || ic == -maxPartIds[1]-1)
	continue;
      
      AliVCluster *cluster = static_cast<AliVCluster*>(fCaloClusters->At(ic));
      
      if (!cluster) {
	AliError(Form("%s: Could not receive cluster %d", GetName(), ic));
	continue;
      } 
      
      if (!AcceptCluster(cluster))
	continue;
      
      TLorentzVector nPart;
      cluster->GetMomentum(nPart, fVertex);

      rhovec[NpartAcc] = nPart.Pt();
      ++NpartAcc;
    }
  }

  Double_t rho = 0;

  if (fUseMedian)
    rho = TMath::Median(NpartAcc, rhovec);
  else
    rho = TMath::Mean(NpartAcc, rhovec);
 
  Float_t maxEta = fTrackMaxEta;
  Float_t minEta = fTrackMinEta;
  Float_t maxPhi = fTrackMaxPhi;
  Float_t minPhi = fTrackMinPhi;

  if (maxPhi > TMath::Pi() * 2) maxPhi = TMath::Pi() * 2;
  if (minPhi < 0) minPhi = 0;

  Double_t area = (maxEta - minEta) * (maxPhi - minPhi);

  if (area > 0) {
    rho *= NpartAcc / area;
    fRho->SetVal(rho);
  } 
  else {
    AliError(Form("%s: Area <= 0 %f", GetName(), area));
    return kFALSE;
  }

  if (fScaleFunction) {
    Double_t rhoScaled = rho * GetScaleFactor(fCent);
    fRhoScaled->SetVal(rhoScaled);
  }

  return kTRUE;
}
