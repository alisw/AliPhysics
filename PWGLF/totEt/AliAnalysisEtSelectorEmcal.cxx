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

AliAnalysisEtSelectorEmcal::AliAnalysisEtSelectorEmcal(AliAnalysisEtCuts* cuts):AliAnalysisEtSelector(cuts)
,fInitialized(kFALSE)
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
{ // Init
    
    AliAnalysisEtSelector::Init(event);
    Printf("Initializing selector for run: %d", event->GetRunNumber());
    fInitialized = kTRUE;
    return 0;
}

TRefArray* AliAnalysisEtSelectorEmcal::GetClusters()
{ // Get clusters
    
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
  return cl.E() > fCuts->GetReconstructedEmcalClusterEnergyCut();
}

Bool_t AliAnalysisEtSelectorEmcal::PassMinEnergyCut(const TParticle& p) const
{
    return p.Energy() > fCuts->GetReconstructedEmcalClusterEnergyCut();
}

Bool_t AliAnalysisEtSelectorEmcal::PassDistanceToBadChannelCut(const AliESDCaloCluster& ) const
{
    return kTRUE;
}

Bool_t AliAnalysisEtSelectorEmcal::PassTrackMatchingCut(const AliESDCaloCluster& ) const
{
    return kTRUE;
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






