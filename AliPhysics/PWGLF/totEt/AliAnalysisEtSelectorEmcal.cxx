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
#include <iostream>
#include "AliStack.h"

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
Bool_t AliAnalysisEtSelectorEmcal::PassMinEnergyCut(Double_t e) const
{
  return e > fCuts->GetReconstructedEmcalClusterEnergyCut();
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

Bool_t AliAnalysisEtSelectorEmcal::CutGeometricalAcceptance(const TParticle& part)
{
  float myphi = part.Phi();
  myphi = AliAnalysisEtSelector::ShiftAngle(myphi);
    return TMath::Abs(part.Eta()) < fCuts->GetGeometryEmcalEtaAccCut()
    && myphi < fCuts->GetGeometryEmcalPhiAccMaxCut()*TMath::Pi()/180.
    && myphi > fCuts->GetGeometryEmcalPhiAccMinCut()*TMath::Pi()/180.;
}

Bool_t AliAnalysisEtSelectorEmcal::CutGeometricalAcceptance(const AliVTrack& part)
{
  float myphi = part.Phi();
  myphi = AliAnalysisEtSelector::ShiftAngle(myphi);
  return TMath::Abs(part.Eta()) < fCuts->GetGeometryEmcalEtaAccCut()
    && myphi < fCuts->GetGeometryEmcalPhiAccMaxCut()*TMath::Pi()/180.
    && myphi > fCuts->GetGeometryEmcalPhiAccMinCut()*TMath::Pi()/180.;
}

Bool_t AliAnalysisEtSelectorEmcal::CutGeometricalAcceptance(const AliESDCaloCluster& cluster)
{
  Float_t pos[3];
  cluster.GetPosition(pos);
  TVector3 cp(pos);
  float myphi = cp.Phi();
  myphi = AliAnalysisEtSelector::ShiftAngle(myphi);
  return TMath::Abs(cp.Eta()) < fCuts->GetGeometryEmcalEtaAccCut()
    && myphi < fCuts->GetGeometryEmcalPhiAccMaxCut()*TMath::Pi()/180.
    && myphi > fCuts->GetGeometryEmcalPhiAccMinCut()*TMath::Pi()/180.;
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

UInt_t AliAnalysisEtSelectorEmcal::GetLabel(const AliESDCaloCluster *cluster, AliStack *stack){
  UInt_t correctLabel = TMath::Abs(cluster->GetLabel());
  correctLabel = GetFirstMotherNotFromDetectorCover(correctLabel,*stack);

  //Now we want to see if this particle is really just something that converted in the cover of the detector and if so, override the label
//   if(!stack){return 0;}
//   if( correctLabel > 0 && stack->IsSecondaryFromMaterial(correctLabel)){//if this is flagged as a secondary then we look to see where it really came from
//     TParticle *hitParticle = stack->Particle(correctLabel);
//     if(hitParticle){
//       Bool_t partVtxSecondary = (TMath::Sqrt(hitParticle->Vx()*hitParticle->Vx() + hitParticle->Vy()*hitParticle->Vy()) >420);
//       if(partVtxSecondary){//at this point we have something which converted near the detector.  Let's find the mother particle
// 	UInt_t mothIdx = stack->Particle(correctLabel)->GetMother(0);
// 	if(mothIdx>0){
// 	  TParticle *mother = stack->Particle(mothIdx);
// 	  if(mother){
// 	    if(AliAnalysisEtSelector::CutGeometricalAcceptance(*mother)){//and the mother is in the acceptance
// 	      if( !(mother->GetPdgCode()==fgPi0Code)){//some of these are decays that just happen this far out
// 		cout<<"I am declaring that "<<hitParticle->GetName()<<" with a vertex of "<< TMath::Sqrt(hitParticle->Vx()*hitParticle->Vx() + hitParticle->Vy()*hitParticle->Vy()) <<" is actually "<<mother->GetName()<<endl;
// 		cout<<"ID check "<<mothIdx<<" vs "<<mother->GetUniqueID()<<endl;
// 		//so now we know that the particle originated near the cover and within the acceptance of the detector
// 		return mothIdx;
// 	      }
// 	    }
// 	  }
// 	}
//       }

//     }
//   }
  return  correctLabel ; // DS: should this line be inside n>0 check, and return another value if n<=0 ? 

}
