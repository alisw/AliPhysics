//_________________________________________________________________________
//  Utility Class for transverse energy studies
//  Selection class for EMCAL
//
//*-- Authors: Oystein Djuvsland (Bergen)
//_________________________________________________________________________


#include "AliAnalysisEtSelectorEmcal.h"
#include "AliAnalysisEtCuts.h"
#include "AliEMCALTrack.h"
#include "TParticle.h"

AliAnalysisEtSelectorEmcal::AliAnalysisEtSelectorEmcal(AliAnalysisEtCuts* cuts):AliAnalysisEtSelector(cuts)
{

}
AliAnalysisEtSelectorEmcal::AliAnalysisEtSelectorEmcal():AliAnalysisEtSelector()
{
}

AliAnalysisEtSelectorEmcal::~AliAnalysisEtSelectorEmcal()
{

}

void AliAnalysisEtSelectorEmcal::Init()
{
    AliAnalysisEtSelector::Init();
}

Int_t AliAnalysisEtSelectorEmcal::Init(const AliESDEvent* event)
{   // Init

    AliAnalysisEtSelector::Init(event);
    Printf("Initializing selector for run: %d", event->GetRunNumber());
    fInitialized = kTRUE;
    return 0;
}

TRefArray* AliAnalysisEtSelectorEmcal::GetClusters()
{   // Get clusters

    if(!fClusterArray) fClusterArray = new TRefArray;

    if(fClusterArray)
    {
        fEvent->GetEMCALClusters(fClusterArray);
    }
    else
    {
        Printf("Could not initialize cluster array");
    }
    return fClusterArray;
}

Bool_t AliAnalysisEtSelectorEmcal::PassMinEnergyCut(const AliESDCaloCluster& cl) const
{
  Float_t pos[3];
  cl.GetPosition(pos);
  TVector3 cp(pos);
  return TMath::Sin(cp.Theta())*cl.E() > fCuts->GetReconstructedEmcalClusterEnergyCut();
}

Bool_t AliAnalysisEtSelectorEmcal::PassMinEnergyCut(const TParticle& p) const
{
  return TMath::Sin(p.Theta())*p.Energy() > fCuts->GetReconstructedEmcalClusterEnergyCut();
}

Bool_t AliAnalysisEtSelectorEmcal::PassDistanceToBadChannelCut(const AliESDCaloCluster& ) const
{
    return kTRUE;
}

Bool_t AliAnalysisEtSelectorEmcal::PassTrackMatchingCut(const AliESDCaloCluster& cluster) const
{

  Int_t nTracksMatched = cluster.GetNTracksMatched();
  if(nTracksMatched == 0){
    return kTRUE;
  }
  
  Int_t trackMatchedIndex = cluster.GetTrackMatchedIndex();
  if(trackMatchedIndex < 0){
    return kTRUE;
  }
  AliVParticle *track = fEvent->GetTrack(trackMatchedIndex);
  if(track->Pt()<0.5) return kTRUE;//Track matches below about 500 MeV are mostly random.  It takes ~460 MeV to reach the EMCal

  //Float_t recoE = cluster.E();
  //Float_t pos[3];

//     //cluster.GetPosition(pos);
//     Int_t trackMatchIdx = cluster.GetTrackMatchedIndex();
//     //Double_t distance = 9999.0;
//     if(trackMatchIdx>-1)
//     {
//       return kTRUE;
//       //distance = CalcTrackClusterDistance(pos, &trackMatchIdx);
//     }


    return kFALSE;
    //return distance > fCuts->GetEmcalTrackDistanceCut();
}

Bool_t AliAnalysisEtSelectorEmcal::CutGeometricalAcceptance(const TParticle& part) const
{
    return TMath::Abs(part.Eta()) < fCuts->GetGeometryEmcalEtaAccCut()
           && part.Phi() < fCuts->GetGeometryEmcalPhiAccMaxCut()*TMath::Pi()/180.
           && part.Phi() > fCuts->GetGeometryEmcalPhiAccMinCut()*TMath::Pi()/180.;
}

Bool_t AliAnalysisEtSelectorEmcal::CutGeometricalAcceptance(const AliVTrack& part) const
{
    return TMath::Abs(part.Eta()) < fCuts->GetGeometryEmcalEtaAccCut()
           && part.Phi() < fCuts->GetGeometryEmcalPhiAccMaxCut()*TMath::Pi()/180.
           && part.Phi() > fCuts->GetGeometryEmcalPhiAccMinCut()*TMath::Pi()/180.;
}


Double_t
AliAnalysisEtSelectorEmcal::CalcTrackClusterDistance(const Float_t clsPos[3],Int_t *trkMatchId) const
{   // calculate distance between cluster and closest track

    Double_t trkPos[3] = {0,0,0};

    Int_t bestTrkMatchId = -1;
    Double_t distance = 9999; // init to a big number

    Double_t dist = 0;
    Double_t distX = 0, distY = 0, distZ = 0;

    for (Int_t iTrack = 0; iTrack < fEvent->GetNumberOfTracks(); iTrack++) {
        AliESDtrack *track = static_cast<AliESDtrack*>( fEvent->GetTrack(iTrack) );
        if (!track) {
            AliError(Form("ERROR: Could not get track %d", iTrack));
            continue;
        }

        // check for approx. eta and phi range before we propagate..
        // TBD

        AliEMCALTrack emctrack(*track);
        if (!emctrack.PropagateToGlobal(clsPos[0],clsPos[1],clsPos[2],0.,0.) ) {
            continue;
        }
        emctrack.GetXYZ(trkPos);

        distX = clsPos[0]-trkPos[0];
        distY = clsPos[1]-trkPos[1];
        distZ = clsPos[2]-trkPos[2];
        dist = TMath::Sqrt(distX*distX + distY*distY + distZ*distZ);

        if (dist < distance) {
            distance = dist;
            bestTrkMatchId = iTrack;
        }
    } // iTrack

    // printf("CalcTrackClusterDistance: bestTrkMatch %d origTrkMatch %d distance %f\n", bestTrkMatchId, *trkMatchId, distance);
    *trkMatchId = bestTrkMatchId;
    return distance;
}





