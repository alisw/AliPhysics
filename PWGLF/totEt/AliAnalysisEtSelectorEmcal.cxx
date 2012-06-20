//_________________________________________________________________________
//  Utility Class for transverse energy studies
//  Selection class for EMCAL
//
//*-- Authors: Oystein Djuvsland (Bergen)
//_________________________________________________________________________


#include "AliAnalysisEtSelectorEmcal.h"
#include "AliAnalysisEtCuts.h"
#include "TParticle.h"
void AliAnalysisEtSelectorEmcal::SetEvent(const AliESDEvent* event)
{
    fEvent = event;
    if(!fInitialized) Init(event);
}

AliAnalysisEtSelectorEmcal::AliAnalysisEtSelectorEmcal(AliAnalysisEtCuts* cuts)
{

}

AliAnalysisEtSelectorEmcal::~AliAnalysisEtSelector()
{

}

void AliAnalysisEtSelectorEmcal::Init()
{
    AliAnalysisEtSelector::Init();
}

Int_t AliAnalysisEtSelectorEmcal::Init(const AliESDEvent* event)
{
    
    AliAnalysisEtSelector::Init(event);
    Printf("Initializing selector for run: %d", event->GetRunNumber());
    fInitialized = kTRUE;
    return 0;
}

TRefArray* AliAnalysisEtSelectorEmcal::GetClusters()
{
    
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

Bool_t AliAnalysisEtSelectorEmcal::CutMinEnergy(const AliESDCaloCluster& cl) const
{
  return cl.E() > fCuts->GetReconstructedEmcalClusterEnergyCut();
}

Bool_t AliAnalysisEtSelectorEmcal::CutMinEnergy(const TParticle& p) const
{
    return p.Energy() > fCuts->GetReconstructedEmcalClusterEnergyCut();
}

Bool_t AliAnalysisEtSelectorEmcal::CutDistanceToBadChannel(const AliESDCaloCluster& ) const
{
    return AliAnalysisEtSelector::CutDistanceToBadChannel();
}

Bool_t AliAnalysisEtSelectorEmcal::CutTrackMatching(const AliESDCaloCluster& ) const
{
    return AliAnalysisEtSelector::CutTrackMatching();
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






