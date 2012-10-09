/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes hereby granted      *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//_________________________________________________________________________
// Class for analysis of particle isolation
// Input is selected particles put in AOD branch (AliAODPWG4ParticleCorrelation)
//
//  Class created from old AliPHOSGammaJet 
//  (see AliRoot versions previous Release 4-09)
//
// -- Author: Gustavo Conesa (LNF-INFN) 

//-Yaxian Mao (add the possibility for different IC method with different pt range, 01/10/2010)
//////////////////////////////////////////////////////////////////////////////


// --- ROOT system --- 
#include <TClonesArray.h>
#include <TList.h>
#include <TObjString.h>
#include <TH2F.h>
#include <TClass.h>

// --- Analysis system --- 
#include "AliAnaParticleIsolation.h" 
#include "AliCaloTrackReader.h"
#include "AliIsolationCut.h"
#include "AliNeutralMesonSelection.h"
#include "AliAODPWG4ParticleCorrelation.h"
#include "AliMCAnalysisUtils.h"
#include "AliVTrack.h"
#include "AliVCluster.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"

ClassImp(AliAnaParticleIsolation)

//______________________________________________________________________________
AliAnaParticleIsolation::AliAnaParticleIsolation() : 
AliAnaCaloTrackCorrBaseClass(),   fCalorimeter(""), 
fReMakeIC(0),                     fMakeSeveralIC(0),               
fFillPileUpHistograms(0),
fFillTMHisto(0),                  fFillSSHisto(0),
// Several IC
fNCones(0),                       fNPtThresFrac(0), 
fConeSizes(),                     fPtThresholds(),                 
fPtFractions(),                   fSumPtThresholds(),
// Histograms
fhEIso(0),                        fhPtIso(0),                       
fhPhiIso(0),                      fhEtaIso(0),                     fhEtaPhiIso(0), 
fhEtaPhiNoIso(0), 
fhPtNoIso(0),                     fhPtDecayIso(0),                 fhPtDecayNoIso(0),
fhEtaPhiDecayIso(0),              fhEtaPhiDecayNoIso(0), 
fhConeSumPt(0),                   fhPtInCone(0),                   
fhPtInConePileUp(0),              fhPtInConeCent(0),
fhFRConeSumPt(0),                 fhPtInFRCone(0),                 fhPhiUEConeSumPt(0),
fhEtaUEConeSumPt(0),              fhEtaBand(0),                    fhPhiBand(0),
fhConeSumPtEtaUESub(0),           fhConeSumPtPhiUESub(0),
// MC histograms
fhPtIsoPrompt(0),                 fhPhiIsoPrompt(0),               fhEtaIsoPrompt(0), 
fhPtThresIsolatedPrompt(),        fhPtFracIsolatedPrompt(),        fhPtSumIsolatedPrompt(),
fhPtIsoFragmentation(0),          fhPhiIsoFragmentation(0),        fhEtaIsoFragmentation(0), 
fhPtThresIsolatedFragmentation(), fhPtFracIsolatedFragmentation(), fhPtSumIsolatedFragmentation(),
fhPtIsoPi0(0),                    fhPhiIsoPi0(0),                  fhEtaIsoPi0(0),
fhPtThresIsolatedPi0(),           fhPtFracIsolatedPi0(),           fhPtSumIsolatedPi0(),
fhPtIsoPi0Decay(0),               fhPhiIsoPi0Decay(0),             fhEtaIsoPi0Decay(0),
fhPtThresIsolatedPi0Decay(),      fhPtFracIsolatedPi0Decay(),      fhPtSumIsolatedPi0Decay(),
fhPtIsoEtaDecay(0),               fhPhiIsoEtaDecay(0),             fhEtaIsoEtaDecay(0),
fhPtThresIsolatedEtaDecay(),      fhPtFracIsolatedEtaDecay(),      fhPtSumIsolatedEtaDecay(),
fhPtIsoOtherDecay(0),             fhPhiIsoOtherDecay(0),           fhEtaIsoOtherDecay(0), 
fhPtThresIsolatedOtherDecay(),    fhPtFracIsolatedOtherDecay(),    fhPtSumIsolatedOtherDecay(),
//fhPtIsoConversion(0),             fhPhiIsoConversion(0),           fhEtaIsoConversion(0), 
//fhPtThresIsolatedConversion(),    fhPtFracIsolatedConversion(),    fhPtSumIsolatedConversion(),
fhPtIsoHadron(0),                 fhPhiIsoHadron(0),               fhEtaIsoHadron(0), 
fhPtThresIsolatedHadron(),        fhPtFracIsolatedHadron(),        fhPtSumIsolatedHadron(),
fhPtNoIsoPi0(0),                  fhPtNoIsoPi0Decay(0),             
fhPtNoIsoEtaDecay(0),             fhPtNoIsoOtherDecay(0),
fhPtNoIsoPrompt(0),               fhPtIsoMCPhoton(0),              fhPtNoIsoMCPhoton(0),
//fhPtNoIsoConversion(0),           
fhPtNoIsoFragmentation(0),        fhPtNoIsoHadron(0),
// Hist several IC
fhSumPtLeadingPt(),               fhPtLeadingPt(), 
fhFRSumPtLeadingPt(),             fhFRPtLeadingPt(),
fhPtThresIsolated(),              fhPtFracIsolated(),              fhPtSumIsolated(),
fhEtaPhiPtThresIso(),             fhEtaPhiPtThresDecayIso(),       fhPtPtThresDecayIso(),
fhEtaPhiPtFracIso(),              fhEtaPhiPtFracDecayIso(),        fhPtPtFracDecayIso(),
fhPtPtSumDecayIso(),              fhEtaPhiSumDensityIso(),         fhEtaPhiSumDensityDecayIso(),
fhPtSumDensityIso(),              fhPtSumDensityDecayIso(), 
fhPtFracPtSumIso(),               fhPtFracPtSumDecayIso(),      
fhEtaPhiFracPtSumIso(),           fhEtaPhiFracPtSumDecayIso(),
// Cluster control histograms
fhTrackMatchedDEta(),             fhTrackMatchedDPhi(),           fhTrackMatchedDEtaDPhi(),
fhdEdx(),                         fhEOverP(),                     fhTrackMatchedMCParticle(),
fhELambda0() ,                    fhELambda1(),                   fhELambda0SSBkg(),
fhELambda0TRD(),                  fhELambda1TRD(),
fhELambda0MCPhoton(),             fhELambda0MCPi0(),              fhELambda0MCPi0Decay(),
fhELambda0MCEtaDecay(),           fhELambda0MCOtherDecay(),       fhELambda0MCHadron(),
// Number of local maxima in cluster
fhNLocMax(),
fhELambda0LocMax1(),              fhELambda1LocMax1(),
fhELambda0LocMax2(),              fhELambda1LocMax2(),
fhELambda0LocMaxN(),              fhELambda1LocMaxN(),
// PileUp
fhTimeENoCut(0),                  fhTimeESPD(0),                  fhTimeESPDMulti(0),
fhTimeNPileUpVertSPD(0),          fhTimeNPileUpVertTrack(0),
fhTimeNPileUpVertContributors(0),
fhTimePileUpMainVertexZDistance(0), fhTimePileUpMainVertexZDiamond(0),
// Histograms settings
fHistoNPtSumBins(0),              fHistoPtSumMax(0.),              fHistoPtSumMin(0.),
fHistoNPtInConeBins(0),           fHistoPtInConeMax(0.),           fHistoPtInConeMin(0.)
{
  //default ctor
  
  //Initialize parameters
  InitParameters();
  
  for(Int_t i = 0; i < 5 ; i++)
  { 
    fConeSizes[i]      = 0 ; 
    
    fhPtSumIsolatedPrompt       [i] = 0 ;  
    fhPtSumIsolatedFragmentation[i] = 0 ;  
    fhPtSumIsolatedPi0Decay     [i] = 0 ;  
    fhPtSumIsolatedPi0          [i] = 0 ;  
    fhPtSumIsolatedEtaDecay     [i] = 0 ;  
    fhPtSumIsolatedOtherDecay   [i] = 0 ;  
//    fhPtSumIsolatedConversion   [i] = 0 ;  
    fhPtSumIsolatedHadron      [i] = 0 ;  
    
    for(Int_t j = 0; j < 5 ; j++)
    { 
      fhPtThresIsolated      [i][j] = 0 ;  
      fhPtFracIsolated       [i][j] = 0 ; 
      fhPtSumIsolated        [i][j] = 0 ;
      
      fhEtaPhiPtThresIso     [i][j] = 0 ;
      fhEtaPhiPtThresDecayIso[i][j] = 0 ;
      fhPtPtThresDecayIso    [i][j] = 0 ;
      
      fhEtaPhiPtFracIso      [i][j] = 0 ;
      fhEtaPhiPtFracDecayIso [i][j] = 0 ;
      fhPtPtFracDecayIso     [i][j] = 0 ;
      fhPtPtSumDecayIso      [i][j] = 0 ;
      fhPtSumDensityIso      [i][j] = 0 ;
      fhPtSumDensityDecayIso [i][j] = 0 ;
      fhEtaPhiSumDensityIso      [i][j] = 0 ;
      fhEtaPhiSumDensityDecayIso [i][j] = 0 ;
      fhPtFracPtSumIso           [i][j] = 0 ;
      fhPtFracPtSumDecayIso      [i][j] = 0 ;
      fhEtaPhiFracPtSumIso       [i][j] = 0 ;
      fhEtaPhiFracPtSumDecayIso  [i][j] = 0 ;
      
      fhPtThresIsolatedPrompt       [i][j] = 0 ;  
      fhPtThresIsolatedFragmentation[i][j] = 0 ; 
      fhPtThresIsolatedPi0Decay     [i][j] = 0 ; 
      fhPtThresIsolatedPi0          [i][j] = 0 ; 
      fhPtThresIsolatedEtaDecay     [i][j] = 0 ; 
      fhPtThresIsolatedOtherDecay   [i][j] = 0 ;  
//      fhPtThresIsolatedConversion   [i][j] = 0 ;  
      fhPtThresIsolatedHadron      [i][j] = 0 ;  
      
      fhPtFracIsolatedPrompt        [i][j] = 0 ;  
      fhPtFracIsolatedFragmentation [i][j] = 0 ;  
      fhPtFracIsolatedPi0           [i][j] = 0 ;  
      fhPtFracIsolatedPi0Decay      [i][j] = 0 ;  
      fhPtFracIsolatedEtaDecay      [i][j] = 0 ;  
      fhPtFracIsolatedOtherDecay    [i][j] = 0 ;  
//      fhPtFracIsolatedConversion    [i][j] = 0 ;
      fhPtFracIsolatedHadron       [i][j] = 0 ;  
      
    }  
  } 
  
  for(Int_t i = 0; i < 5 ; i++)
  { 
    fPtFractions    [i] = 0 ; 
    fPtThresholds   [i] = 0 ;
    fSumPtThresholds[i] = 0 ;
  } 

  
  for(Int_t i = 0; i < 2 ; i++)
  { 
    fhTrackMatchedDEta[i] = 0 ;             fhTrackMatchedDPhi[i] = 0 ;           fhTrackMatchedDEtaDPhi  [i] = 0 ;
    fhdEdx            [i] = 0 ;             fhEOverP          [i] = 0 ;           fhTrackMatchedMCParticle[i] = 0 ;
    fhELambda0        [i] = 0 ;             fhELambda1        [i] = 0 ; 
    fhELambda0TRD     [i] = 0 ;             fhELambda1TRD     [i] = 0 ;
    
    fhELambda0MCPhoton  [i] = 0 ;           fhELambda0MCPi0       [i] = 0 ;       fhELambda0MCPi0Decay[i] = 0 ;
    fhELambda0MCEtaDecay[i] = 0 ;           fhELambda0MCOtherDecay[i] = 0 ;       fhELambda0MCHadron  [i] = 0 ;

    
    // Number of local maxima in cluster
    fhNLocMax        [i] = 0 ;
    fhELambda0LocMax1[i] = 0 ;              fhELambda1LocMax1[i] = 0 ;
    fhELambda0LocMax2[i] = 0 ;              fhELambda1LocMax2[i] = 0 ;
    fhELambda0LocMaxN[i] = 0 ;              fhELambda1LocMaxN[i] = 0 ;
    
  } 
  
}

//_________________________________________________________________
void AliAnaParticleIsolation::FillPileUpHistograms(Int_t clusterID) 
{
  // Fill some histograms to understand pile-up
  if(!fFillPileUpHistograms) return;
  
  if(clusterID < 0 ) 
  {
    printf("AliAnaParticleIsolation::FillPileUpHistograms(), ID of cluster = %d, not possible! ", clusterID);
    return;
  }
  
  Int_t iclus = -1;
  TObjArray* clusters = 0x0;
  if     (fCalorimeter == "EMCAL") clusters = GetEMCALClusters();
  else if(fCalorimeter == "PHOS" ) clusters = GetPHOSClusters();
  
  Float_t energy = 0;
  Float_t time   = -1000;

  if(clusters)
  {
    AliVCluster *cluster = FindCluster(clusters,clusterID,iclus); 
    energy = cluster->E();
    time   = cluster->GetTOF()*1e9;
  } 
  
  //printf("E %f, time %f\n",energy,time);
  AliVEvent * event = GetReader()->GetInputEvent();
  
  fhTimeENoCut->Fill(energy,time);
  if(GetReader()->IsPileUpFromSPD())     fhTimeESPD     ->Fill(energy,time);
  if(event->IsPileupFromSPDInMultBins()) fhTimeESPDMulti->Fill(energy,time);
  
  if(energy < 8) return; // Fill time figures for high energy clusters not too close to trigger threshold
  
  AliESDEvent* esdEv = dynamic_cast<AliESDEvent*> (event);
  AliAODEvent* aodEv = dynamic_cast<AliAODEvent*> (event);
  
  // N pile up vertices
  Int_t nVerticesSPD    = -1;
  Int_t nVerticesTracks = -1;
  
  if      (esdEv)
  {
    nVerticesSPD    = esdEv->GetNumberOfPileupVerticesSPD();
    nVerticesTracks = esdEv->GetNumberOfPileupVerticesTracks();
    
  }//ESD
  else if (aodEv)
  {
    nVerticesSPD    = aodEv->GetNumberOfPileupVerticesSPD();
    nVerticesTracks = aodEv->GetNumberOfPileupVerticesTracks();
  }//AOD
  
  fhTimeNPileUpVertSPD  ->Fill(time,nVerticesSPD);
  fhTimeNPileUpVertTrack->Fill(time,nVerticesTracks);
  
  //printf("Is SPD %d, Is SPD Multi %d, n spd %d, n track %d\n", 
  //       GetReader()->IsPileUpFromSPD(),event->IsPileupFromSPDInMultBins(),nVerticesSPD,nVerticesTracks);
  
  Int_t ncont = -1;
  Float_t z1 = -1, z2 = -1;
  Float_t diamZ = -1;
  for(Int_t iVert=0; iVert<nVerticesSPD;iVert++)
  {
    if      (esdEv)
    {
      const AliESDVertex* pv=esdEv->GetPileupVertexSPD(iVert);
      ncont=pv->GetNContributors();
      z1 = esdEv->GetPrimaryVertexSPD()->GetZ();
      z2 = pv->GetZ();
      diamZ = esdEv->GetDiamondZ();
    }//ESD
    else if (aodEv)
    {
      AliAODVertex *pv=aodEv->GetVertex(iVert);
      if(pv->GetType()!=AliAODVertex::kPileupSPD) continue;
      ncont=pv->GetNContributors();
      z1=aodEv->GetPrimaryVertexSPD()->GetZ();
      z2=pv->GetZ();
      diamZ = aodEv->GetDiamondZ();
    }// AOD
    
    Double_t distZ  = TMath::Abs(z2-z1);
    diamZ  = TMath::Abs(z2-diamZ);
    
    fhTimeNPileUpVertContributors  ->Fill(time,ncont);
    fhTimePileUpMainVertexZDistance->Fill(time,distZ);
    fhTimePileUpMainVertexZDiamond ->Fill(time,diamZ);
    
  }// loop
}

//________________________________________________________________________________________________
void AliAnaParticleIsolation::FillTrackMatchingShowerShapeControlHistograms(const Bool_t isolated,
                                                                            const Int_t  clusterID,
                                                                            const Int_t  nMaxima,
                                                                            const Int_t  mcTag,
                                                                            const TObjArray * plCTS, 
                                                                            const TObjArray * plNe, 
                                                                            AliAODPWG4ParticleCorrelation  *pCandidate,
                                                                            const AliCaloTrackReader * reader, 
                                                                            const AliCaloPID * pid)
{
  // Fill Track matching and Shower Shape control histograms  
  if(!fFillTMHisto &&  !fFillSSHisto) return;
  
  if(clusterID < 0 ) 
  {
    printf("AliAnaParticleIsolation::FillTrackMatchingShowerShapeControlHistograms(), ID of cluster = %d, not possible! ", clusterID);
    return;
  }
  
  Int_t iclus = -1;
  TObjArray* clusters = 0x0;
  if     (fCalorimeter == "EMCAL") clusters = GetEMCALClusters();
  else if(fCalorimeter == "PHOS" ) clusters = GetPHOSClusters();
  
  if(clusters)
  {
    
    AliVCluster *cluster = FindCluster(clusters,clusterID,iclus); 
    Float_t energy = cluster->E();
    
    if(fFillSSHisto)
    {
      fhELambda0[isolated]->Fill(energy, cluster->GetM02() );  
      fhELambda1[isolated]->Fill(energy, cluster->GetM20() );  
      
      if(IsDataMC())
      {
        if     (GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPrompt) ||        
                GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCFragmentation)) fhELambda0MCPhoton    [isolated]->Fill(energy, cluster->GetM02());
        else if(GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPi0))           fhELambda0MCPi0       [isolated]->Fill(energy, cluster->GetM02());
        else if(GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPi0Decay))      fhELambda0MCPi0Decay  [isolated]->Fill(energy, cluster->GetM02());
        else if(GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCEtaDecay))      fhELambda0MCEtaDecay  [isolated]->Fill(energy, cluster->GetM02());
        else if(GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCOtherDecay))    fhELambda0MCOtherDecay[isolated]->Fill(energy, cluster->GetM02());
       
        //        else if(GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCConversion))    fhPtNoIsoConversion   ->Fill(energy, cluster->GetM02());
        else if(!GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCElectron))     fhELambda0MCHadron    [isolated]->Fill(energy, cluster->GetM02());        
      
      }
      
      if(fCalorimeter == "EMCAL" && GetModuleNumber(cluster) > 5) // TO DO: CHANGE FOR 2012
      {
        fhELambda0TRD[isolated]->Fill(energy, cluster->GetM02() );  
        fhELambda1TRD[isolated]->Fill(energy, cluster->GetM20() );  
      }
      
      fhNLocMax[isolated]->Fill(energy,nMaxima);
      if     (nMaxima==1) { fhELambda0LocMax1[isolated]->Fill(energy,cluster->GetM02()); fhELambda1LocMax1[isolated]->Fill(energy,cluster->GetM20()); }
      else if(nMaxima==2) { fhELambda0LocMax2[isolated]->Fill(energy,cluster->GetM02()); fhELambda1LocMax2[isolated]->Fill(energy,cluster->GetM20()); }
      else                { fhELambda0LocMaxN[isolated]->Fill(energy,cluster->GetM02()); fhELambda1LocMaxN[isolated]->Fill(energy,cluster->GetM20()); }

      if(isolated==0)
	{
	  //Analyse non-isolated events
	  Int_t   n = 0; 
	  Int_t nfrac = 0;
	  Bool_t  iso  = kFALSE ;
	  Float_t coneptsum = 0 ;
	  GetIsolationCut()->SetPtThresholdMax(1.);
	  GetIsolationCut()->MakeIsolationCut(plCTS,   plNe, 
					      reader, pid,
					      kFALSE, pCandidate, "", 
					      n,nfrac,coneptsum, iso);
	  if (!iso) fhELambda0SSBkg->Fill(energy, cluster->GetM02());
	  
      
	  if(GetDebug() > 0) printf("AliAnaParticleIsolation::MakeAnalysisFillHistograms() - Energy Sum in Isolation Cone %2.2f\n", coneptsum);    
	}
      GetIsolationCut()->SetPtThresholdMax(10000.);

    } // SS histo fill        
    
    
    if(fFillTMHisto)
    {
      Float_t dZ  = cluster->GetTrackDz();
      Float_t dR  = cluster->GetTrackDx();
      
      if(cluster->IsEMCAL() && GetCaloUtils()->IsRecalculationOfClusterTrackMatchingOn())
      {
        dR = 2000., dZ = 2000.;
        GetCaloUtils()->GetEMCALRecoUtils()->GetMatchedResiduals(cluster->GetID(),dZ,dR);
      }
      
      //printf("ParticleIsolation: dPhi %f, dEta %f\n",dR,dZ);
      if(fhTrackMatchedDEta[isolated] && TMath::Abs(dR) < 999)
      {
        fhTrackMatchedDEta[isolated]->Fill(energy,dZ);
        fhTrackMatchedDPhi[isolated]->Fill(energy,dR);
        if(energy > 0.5) fhTrackMatchedDEtaDPhi[isolated]->Fill(dZ,dR);
      }
      
      // Check dEdx and E/p of matched clusters
      
      if(TMath::Abs(dZ) < 0.05 && TMath::Abs(dR) < 0.05)
      {
        
        AliVTrack *track = GetCaloUtils()->GetMatchedTrack(cluster, GetReader()->GetInputEvent());
        
        if(track) 
        {
          Float_t dEdx = track->GetTPCsignal();
          fhdEdx[isolated]->Fill(cluster->E(), dEdx);
          
          Float_t eOverp = cluster->E()/track->P();
          fhEOverP[isolated]->Fill(cluster->E(),  eOverp);
        }
        //else 
        //  printf("AliAnaParticleIsolation::FillTrackMatchingShowerShapeHistograms() - Residual OK but (dR, dZ)= (%2.4f,%2.4f) no track associated WHAT? \n", dR,dZ);
        
        
        if(IsDataMC())
        {
          if ( !GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCConversion)  )
          {
            if       ( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPi0)      ||
                       GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCEta)       ) fhTrackMatchedMCParticle[isolated]->Fill(energy, 2.5 );
            else if  ( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPhoton)    ) fhTrackMatchedMCParticle[isolated]->Fill(energy, 0.5 );
            else if  ( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCElectron)  ) fhTrackMatchedMCParticle[isolated]->Fill(energy, 1.5 );
            else                                                                                   fhTrackMatchedMCParticle[isolated]->Fill(energy, 3.5 );
            
          }
          else
          {
            if       ( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPi0)      ||
                       GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCEta)       ) fhTrackMatchedMCParticle[isolated]->Fill(energy, 6.5 );
            else if  ( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPhoton)    ) fhTrackMatchedMCParticle[isolated]->Fill(energy, 4.5 );
            else if  ( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCElectron)  ) fhTrackMatchedMCParticle[isolated]->Fill(energy, 5.5 );
            else                                                                                   fhTrackMatchedMCParticle[isolated]->Fill(energy, 7.5 );
          }                    
          
        }  // MC           
        
      } // match window            
      
    }// TM histos fill
    
  } // clusters array available
  
}

//______________________________________________________
TObjString *  AliAnaParticleIsolation::GetAnalysisCuts()
{ 
  //Save parameters used for analysis
  TString parList ; //this will be list of parameters used for this analysis.
  const Int_t buffersize = 255;
  char onePar[buffersize] ;
  
  snprintf(onePar, buffersize,"--- AliAnaParticleIsolation ---\n") ;
  parList+=onePar ;	
  snprintf(onePar, buffersize,"Calorimeter: %s\n",fCalorimeter.Data()) ;
  parList+=onePar ;
  snprintf(onePar, buffersize,"fReMakeIC =%d (Flag for reisolation during histogram filling) \n",fReMakeIC) ;
  parList+=onePar ;
  snprintf(onePar, buffersize,"fMakeSeveralIC=%d (Flag for isolation with several cuts at the same time ) \n",fMakeSeveralIC) ;
  parList+=onePar ;  
  snprintf(onePar, buffersize,"fFillTMHisto=%d (Flag for track matching histograms) \n",fFillTMHisto) ;
  parList+=onePar ;
  snprintf(onePar, buffersize,"fFillSSHisto=%d (Flag for shower shape histograms) \n",fFillSSHisto) ;
  parList+=onePar ;
  
  if(fMakeSeveralIC)
  {
    snprintf(onePar, buffersize,"fNCones =%d (Number of cone sizes) \n",fNCones) ;
    parList+=onePar ;
    snprintf(onePar, buffersize,"fNPtThresFrac=%d (Flag for isolation with several cuts at the same time ) \n",fNPtThresFrac) ;
    parList+=onePar ;
    
    for(Int_t icone = 0; icone < fNCones ; icone++)
    {
      snprintf(onePar, buffersize,"fConeSizes[%d]=%1.2f (isolation cone size) \n",icone, fConeSizes[icone]) ;
      parList+=onePar ;	
    }
    for(Int_t ipt = 0; ipt < fNPtThresFrac ; ipt++)
    {
      snprintf(onePar, buffersize,"fPtThresholds[%d]=%1.2f (isolation pt threshold) \n",ipt, fPtThresholds[ipt]) ;
      parList+=onePar ;	
    }
    for(Int_t ipt = 0; ipt < fNPtThresFrac ; ipt++)
    {
      snprintf(onePar, buffersize,"fPtFractions[%d]=%1.2f (isolation pt fraction threshold) \n",ipt, fPtFractions[ipt]) ;
      parList+=onePar ;	
    }
    for(Int_t ipt = 0; ipt < fNPtThresFrac ; ipt++)
    {
      snprintf(onePar, buffersize,"fSumPtThresholds[%d]=%1.2f (isolation sum pt threshold) \n",ipt, fSumPtThresholds[ipt]) ;
      parList+=onePar ;	
    }		
  }
  
  //Get parameters set in base class.
  parList += GetBaseParametersList() ;
  
  //Get parameters set in IC class.
  if(!fMakeSeveralIC)parList += GetIsolationCut()->GetICParametersList() ;
  
  return new TObjString(parList) ;
  
}

//________________________________________________________
TList *  AliAnaParticleIsolation::GetCreateOutputObjects()
{  
  // Create histograms to be saved in output file and 
  // store them in outputContainer
  TList * outputContainer = new TList() ; 
  outputContainer->SetName("IsolatedParticleHistos") ; 
  
  Int_t   nptbins  = GetHistogramRanges()->GetHistoPtBins();
  Int_t   nphibins = GetHistogramRanges()->GetHistoPhiBins();
  Int_t   netabins = GetHistogramRanges()->GetHistoEtaBins();
  Float_t ptmax    = GetHistogramRanges()->GetHistoPtMax();
  Float_t phimax   = GetHistogramRanges()->GetHistoPhiMax();
  Float_t etamax   = GetHistogramRanges()->GetHistoEtaMax();
  Float_t ptmin    = GetHistogramRanges()->GetHistoPtMin();
  Float_t phimin   = GetHistogramRanges()->GetHistoPhiMin();
  Float_t etamin   = GetHistogramRanges()->GetHistoEtaMin();	
  Int_t   ssbins   = GetHistogramRanges()->GetHistoShowerShapeBins(); 
  Float_t ssmax    = GetHistogramRanges()->GetHistoShowerShapeMax();  
  Float_t ssmin    = GetHistogramRanges()->GetHistoShowerShapeMin();
  Int_t   ntimebins= GetHistogramRanges()->GetHistoTimeBins();         
  Float_t timemax  = GetHistogramRanges()->GetHistoTimeMax();         
  Float_t timemin  = GetHistogramRanges()->GetHistoTimeMin();       

  Int_t   nresetabins = GetHistogramRanges()->GetHistoTrackResidualEtaBins();          
  Float_t resetamax   = GetHistogramRanges()->GetHistoTrackResidualEtaMax();          
  Float_t resetamin   = GetHistogramRanges()->GetHistoTrackResidualEtaMin();
  Int_t   nresphibins = GetHistogramRanges()->GetHistoTrackResidualPhiBins();          
  Float_t resphimax   = GetHistogramRanges()->GetHistoTrackResidualPhiMax();          
  Float_t resphimin   = GetHistogramRanges()->GetHistoTrackResidualPhiMin();  
  
  Int_t   ndedxbins   = GetHistogramRanges()->GetHistodEdxBins();         
  Float_t dedxmax     = GetHistogramRanges()->GetHistodEdxMax();         
  Float_t dedxmin     = GetHistogramRanges()->GetHistodEdxMin();
  Int_t   nPoverEbins = GetHistogramRanges()->GetHistoPOverEBins();       
  Float_t pOverEmax   = GetHistogramRanges()->GetHistoPOverEMax();       
  Float_t pOverEmin   = GetHistogramRanges()->GetHistoPOverEMin();
  
  Int_t   nptsumbins    = fHistoNPtSumBins;
  Float_t ptsummax      = fHistoPtSumMax;
  Float_t ptsummin      = fHistoPtSumMin;	
  Int_t   nptinconebins = fHistoNPtInConeBins;
  Float_t ptinconemax   = fHistoPtInConeMax;
  Float_t ptinconemin   = fHistoPtInConeMin;
  
  Float_t ptthre = GetIsolationCut()->GetPtThreshold();
  Float_t ptfrac = GetIsolationCut()->GetPtFraction();
  Float_t r      = GetIsolationCut()->GetConeSize();
  
  if(!fMakeSeveralIC)
  {
    TString hName [] = {"NoIso",""};
    TString hTitle[] = {"Not isolated"  ,"isolated"};
    if(fFillSSHisto)
    { 
    	fhELambda0SSBkg  = new TH2F
      ("hELambda0SSBkg","Non isolated clusters : E vs #lambda_{0}",nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
      fhELambda0SSBkg->SetYTitle("#lambda_{0}^{2}");
      fhELambda0SSBkg->SetXTitle("E (GeV)");
      outputContainer->Add(fhELambda0SSBkg) ;
    }
    
    for(Int_t iso = 0; iso < 2; iso++)
    {
      if(fFillTMHisto)
      {
        fhTrackMatchedDEta[iso]  = new TH2F
        (Form("hTrackMatchedDEta%s",hName[iso].Data()),
         Form("%s - d#eta of cluster-track vs cluster energy for R = %2.2f, p_{T}^{th} = %2.2f, p_{T}^{fr} = %2.2f",hTitle[iso].Data(),r,ptthre,ptfrac),
         nptbins,ptmin,ptmax,nresetabins,resetamin,resetamax); 
        fhTrackMatchedDEta[iso]->SetYTitle("d#eta");
        fhTrackMatchedDEta[iso]->SetXTitle("E_{cluster} (GeV)");
        
        fhTrackMatchedDPhi[iso]  = new TH2F
        (Form("hTrackMatchedDPhi%s",hName[iso].Data()),
         Form("%s - d#phi of cluster-track vs cluster energy for R = %2.2f, p_{T}^{th} = %2.2f, p_{T}^{fr} = %2.2f",hTitle[iso].Data(),r,ptthre,ptfrac),
         nptbins,ptmin,ptmax,nresphibins,resphimin,resphimax); 
        fhTrackMatchedDPhi[iso]->SetYTitle("d#phi (rad)");
        fhTrackMatchedDPhi[iso]->SetXTitle("E_{cluster} (GeV)");
        
        fhTrackMatchedDEtaDPhi[iso]  = new TH2F
        (Form("hTrackMatchedDEtaDPhi%s",hName[iso].Data()),
         Form("%s - d#eta vs d#phi of cluster-track for R = %2.2f, p_{T}^{th} = %2.2f, p_{T}^{fr} = %2.2f",hTitle[iso].Data(),r,ptthre,ptfrac),       
         nresetabins,resetamin,resetamax,nresphibins,resphimin,resphimax); 
        fhTrackMatchedDEtaDPhi[iso]->SetYTitle("d#phi (rad)");
        fhTrackMatchedDEtaDPhi[iso]->SetXTitle("d#eta");   
        
        outputContainer->Add(fhTrackMatchedDEta[iso]) ; 
        outputContainer->Add(fhTrackMatchedDPhi[iso]) ;
        outputContainer->Add(fhTrackMatchedDEtaDPhi[iso]) ;
        
        fhdEdx[iso]  = new TH2F
        (Form("hdEdx%s",hName[iso].Data()),
         Form("%s - Matched track <dE/dx> vs cluster E for R = %2.2f, p_{T}^{th} = %2.2f, p_{T}^{fr} = %2.2f",hTitle[iso].Data(),r,ptthre,ptfrac), 
         nptbins,ptmin,ptmax,ndedxbins, dedxmin, dedxmax); 
        fhdEdx[iso]->SetXTitle("E (GeV)");
        fhdEdx[iso]->SetYTitle("<dE/dx>");
        outputContainer->Add(fhdEdx[iso]);  
        
        fhEOverP[iso]  = new TH2F
        (Form("hEOverP%s",hName[iso].Data()),
         Form("%s - Matched track E/p vs cluster E for R = %2.2f, p_{T}^{th} = %2.2f, p_{T}^{fr} = %2.2f",hTitle[iso].Data(),r,ptthre,ptfrac), 
         nptbins,ptmin,ptmax,nPoverEbins,pOverEmin,pOverEmax); 
        fhEOverP[iso]->SetXTitle("E (GeV)");
        fhEOverP[iso]->SetYTitle("E/p");
        outputContainer->Add(fhEOverP[iso]);   
        
        if(IsDataMC())
        {
          fhTrackMatchedMCParticle[iso]  = new TH2F
          (Form("hTrackMatchedMCParticle%s",hName[iso].Data()),
           Form("%s - Origin of particle vs energy vs cluster E for R = %2.2f, p_{T}^{th} = %2.2f, p_{T}^{fr} = %2.2f",hTitle[iso].Data(),r,ptthre,ptfrac), 
           nptbins,ptmin,ptmax,8,0,8); 
          fhTrackMatchedMCParticle[iso]->SetXTitle("E (GeV)");   
          //fhTrackMatchedMCParticle[iso]->SetYTitle("Particle type");
          
          fhTrackMatchedMCParticle[iso]->GetYaxis()->SetBinLabel(1 ,"Photon");
          fhTrackMatchedMCParticle[iso]->GetYaxis()->SetBinLabel(2 ,"Electron");
          fhTrackMatchedMCParticle[iso]->GetYaxis()->SetBinLabel(3 ,"Meson Merged");
          fhTrackMatchedMCParticle[iso]->GetYaxis()->SetBinLabel(4 ,"Rest");
          fhTrackMatchedMCParticle[iso]->GetYaxis()->SetBinLabel(5 ,"Conv. Photon");
          fhTrackMatchedMCParticle[iso]->GetYaxis()->SetBinLabel(6 ,"Conv. Electron");
          fhTrackMatchedMCParticle[iso]->GetYaxis()->SetBinLabel(7 ,"Conv. Merged");
          fhTrackMatchedMCParticle[iso]->GetYaxis()->SetBinLabel(8 ,"Conv. Rest");
          
          outputContainer->Add(fhTrackMatchedMCParticle[iso]);         
        }
      }
      
      if(fFillSSHisto)
      {
        fhELambda0[iso]  = new TH2F
        (Form("hELambda0%s",hName[iso].Data()),
         Form("%s cluster : E vs #lambda_{0}",hTitle[iso].Data()),nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
        fhELambda0[iso]->SetYTitle("#lambda_{0}^{2}");
        fhELambda0[iso]->SetXTitle("E (GeV)");
        outputContainer->Add(fhELambda0[iso]) ; 
        
        if(IsDataMC())
        {
          fhELambda0MCPhoton[iso]  = new TH2F
          (Form("hELambda0%s_MCPhoton",hName[iso].Data()),
           Form("%s cluster : E vs #lambda_{0}: Origin is final state photon",hTitle[iso].Data()),nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
          fhELambda0MCPhoton[iso]->SetYTitle("#lambda_{0}^{2}");
          fhELambda0MCPhoton[iso]->SetXTitle("E (GeV)");
          outputContainer->Add(fhELambda0MCPhoton[iso]) ; 
          
          fhELambda0MCPi0[iso]  = new TH2F
          (Form("hELambda0%s_MCPi0",hName[iso].Data()),
           Form("%s cluster : E vs #lambda_{0}: Origin is pi0 (2 #gamma)",hTitle[iso].Data()),nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
          fhELambda0MCPi0[iso]->SetYTitle("#lambda_{0}^{2}");
          fhELambda0MCPi0[iso]->SetXTitle("E (GeV)");
          outputContainer->Add(fhELambda0MCPi0[iso]) ; 
          
          fhELambda0MCPi0Decay[iso]  = new TH2F
          (Form("hELambda0%s_MCPi0Decay",hName[iso].Data()),
           Form("%s cluster : E vs #lambda_{0}: Origin is pi0 decay",hTitle[iso].Data()),nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
          fhELambda0MCPi0Decay[iso]->SetYTitle("#lambda_{0}^{2}");
          fhELambda0MCPi0Decay[iso]->SetXTitle("E (GeV)");
          outputContainer->Add(fhELambda0MCPi0Decay[iso]) ; 
          
          fhELambda0MCEtaDecay[iso]  = new TH2F
          (Form("hELambda0%s_MCEtaDecay",hName[iso].Data()),
           Form("%s cluster : E vs #lambda_{0}: Origin is eta decay",hTitle[iso].Data()),nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
          fhELambda0MCEtaDecay[iso]->SetYTitle("#lambda_{0}^{2}");
          fhELambda0MCEtaDecay[iso]->SetXTitle("E (GeV)");
          outputContainer->Add(fhELambda0MCEtaDecay[iso]) ; 
          
          fhELambda0MCOtherDecay[iso]  = new TH2F
          (Form("hELambda0%s_MCOtherDecay",hName[iso].Data()),
           Form("%s cluster : E vs #lambda_{0}: Origin is other decay",hTitle[iso].Data()),nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
          fhELambda0MCOtherDecay[iso]->SetYTitle("#lambda_{0}^{2}");
          fhELambda0MCOtherDecay[iso]->SetXTitle("E (GeV)");
          outputContainer->Add(fhELambda0MCOtherDecay[iso]) ; 
          
          fhELambda0MCHadron[iso]  = new TH2F
          (Form("hELambda0%s_MCHadron",hName[iso].Data()),
           Form("%s cluster : E vs #lambda_{0}: Origin is hadron",hTitle[iso].Data()),nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
          fhELambda0MCHadron[iso]->SetYTitle("#lambda_{0}^{2}");
          fhELambda0MCHadron[iso]->SetXTitle("E (GeV)");
          outputContainer->Add(fhELambda0MCHadron[iso]) ; 
        } 
        
        fhELambda1[iso]  = new TH2F
        (Form("hELambda1%s",hName[iso].Data()),
         Form("%s cluster: E vs #lambda_{1}",hTitle[iso].Data()),nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
        fhELambda1[iso]->SetYTitle("#lambda_{1}^{2}");
        fhELambda1[iso]->SetXTitle("E (GeV)");
        outputContainer->Add(fhELambda1[iso]) ;
        
        if(fCalorimeter=="EMCAL")
        {
          fhELambda0TRD[iso]  = new TH2F
          (Form("hELambda0TRD%s",hName[iso].Data()),
           Form("%s cluster: E vs #lambda_{0}, SM behind TRD",hTitle[iso].Data()),nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
          fhELambda0TRD[iso]->SetYTitle("#lambda_{0}^{2}");
          fhELambda0TRD[iso]->SetXTitle("E (GeV)");
          outputContainer->Add(fhELambda0TRD[iso]) ; 
          
          fhELambda1TRD[iso]  = new TH2F
          (Form("hELambda1TRD%s",hName[iso].Data()),
           Form("%s cluster: E vs #lambda_{1}, SM behind TRD",hTitle[iso].Data()),nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
          fhELambda1TRD[iso]->SetYTitle("#lambda_{1}^{2}");
          fhELambda1TRD[iso]->SetXTitle("E (GeV)");
          outputContainer->Add(fhELambda1TRD[iso]) ;         
        }
        
        fhNLocMax[iso] = new TH2F
        (Form("hNLocMax%s",hName[iso].Data()),
         Form("%s - Number of local maxima in cluster",hTitle[iso].Data()),
         nptbins,ptmin,ptmax,10,0,10); 
        fhNLocMax[iso]->SetYTitle("N maxima");
        fhNLocMax[iso]->SetXTitle("E (GeV)");
        outputContainer->Add(fhNLocMax[iso]) ;       
        
        fhELambda0LocMax1[iso]  = new TH2F
        (Form("hELambda0LocMax1%s",hName[iso].Data()),
         Form("%s cluster (#eta) pairs: E vs #lambda_{0}, 1 Local maxima",hTitle[iso].Data()),nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
        fhELambda0LocMax1[iso]->SetYTitle("#lambda_{0}^{2}");
        fhELambda0LocMax1[iso]->SetXTitle("E (GeV)");
        outputContainer->Add(fhELambda0LocMax1[iso]) ; 
        
        fhELambda1LocMax1[iso]  = new TH2F
        (Form("hELambda1LocMax1%s",hName[iso].Data()),
         Form("%s cluster (#eta) pairs: E vs #lambda_{1}, 1 Local maxima",hTitle[iso].Data()),nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
        fhELambda1LocMax1[iso]->SetYTitle("#lambda_{1}^{2}");
        fhELambda1LocMax1[iso]->SetXTitle("E (GeV)");
        outputContainer->Add(fhELambda1LocMax1[iso]) ; 
        
        fhELambda0LocMax2[iso]  = new TH2F
        (Form("hELambda0LocMax2%s",hName[iso].Data()),
         Form("%s cluster (#eta) pairs: E vs #lambda_{0}, 2 Local maxima",hTitle[iso].Data()),nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
        fhELambda0LocMax2[iso]->SetYTitle("#lambda_{0}^{2}");
        fhELambda0LocMax2[iso]->SetXTitle("E (GeV)");
        outputContainer->Add(fhELambda0LocMax2[iso]) ; 
        
        fhELambda1LocMax2[iso]  = new TH2F
        (Form("hELambda1LocMax2%s",hName[iso].Data()),
         Form("%s cluster (#eta) pairs: E vs #lambda_{1}, 2 Local maxima",hTitle[iso].Data()),nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
        fhELambda1LocMax2[iso]->SetYTitle("#lambda_{1}^{2}");
        fhELambda1LocMax2[iso]->SetXTitle("E (GeV)");
        outputContainer->Add(fhELambda1LocMax2[iso]) ; 
        
        fhELambda0LocMaxN[iso]  = new TH2F
        ( Form("hELambda0LocMaxN%s",hName[iso].Data()),
         Form("%s cluster (#eta) pairs: E vs #lambda_{0}, N>2 Local maxima",hTitle[iso].Data()),nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
        fhELambda0LocMaxN[iso]->SetYTitle("#lambda_{0}^{2}");
        fhELambda0LocMaxN[iso]->SetXTitle("E (GeV)");
        outputContainer->Add(fhELambda0LocMaxN[iso]) ; 
        
        fhELambda1LocMaxN[iso]  = new TH2F
        (Form("hELambda1LocMaxN%s",hName[iso].Data()),
         Form("%s cluster (#eta) pairs: E vs #lambda_{1}, N>2 Local maxima",hTitle[iso].Data()),nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
        fhELambda1LocMaxN[iso]->SetYTitle("#lambda_{1}^{2}");
        fhELambda1LocMaxN[iso]->SetXTitle("E (GeV)");
        outputContainer->Add(fhELambda1LocMaxN[iso]) ; 
        
      }
    } // control histograms for isolated and non isolated objects
    
    fhConeSumPt  = new TH2F("hConePtSum",
                            Form("#Sigma p_{T} in isolation cone for R = %2.2f",r),
                            nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
    fhConeSumPt->SetYTitle("#Sigma p_{T}");
    fhConeSumPt->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fhConeSumPt) ;
    
    fhPtInCone  = new TH2F("hPtInCone",
                           Form("p_{T} in isolation cone for R = %2.2f",r),
                           nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
    fhPtInCone->SetYTitle("p_{T in cone} (GeV/c)");
    fhPtInCone->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fhPtInCone) ;
    
    if(fFillPileUpHistograms)
    {
      fhPtInConePileUp  = new TH2F("hPtInConePileUp",
                             Form("p_{T} in isolation cone for R = %2.2f, from pile-up (SPD)",r),
                             nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
      fhPtInConePileUp->SetYTitle("p_{T in cone} (GeV/c)");
      fhPtInConePileUp->SetXTitle("p_{T} (GeV/c)");
      outputContainer->Add(fhPtInConePileUp) ;      
    }
    
    fhPtInConeCent  = new TH2F("hPtInConeCent",
                               Form("p_{T} in isolation cone for R = %2.2f",r),
                               100,0,100,nptinconebins,ptinconemin,ptinconemax);
    fhPtInConeCent->SetYTitle("p_{T in cone} (GeV/c)");
    fhPtInConeCent->SetXTitle("centrality");
    outputContainer->Add(fhPtInConeCent) ;
    
    fhFRConeSumPt  = new TH2F("hFRConePtSum",
                              Form("#Sigma p_{T} in the forward region isolation cone for R = %2.2f",r),
                              nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
    fhFRConeSumPt->SetYTitle("#Sigma p_{T}");
    fhFRConeSumPt->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fhFRConeSumPt) ;
    
    fhPtInFRCone  = new TH2F("hPtInFRCone",
                             Form("p_{T} in forward region isolation cone for R = %2.2f",r),
                             nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
    fhPtInFRCone->SetYTitle("p_{T in cone} (GeV/c)");
    fhPtInFRCone->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fhPtInFRCone) ; 
    
    fhPhiUEConeSumPt  = new TH2F("hPhiUEConeSumPt",
                                 Form("p_{T} in phi band around isolation cone for R = %2.2f",r),
                                 nptbins,ptmin,ptmax,4*nptinconebins,ptinconemin,4*ptinconemax);
    fhPhiUEConeSumPt->SetYTitle("p_{T in band} (GeV/c)");
    fhPhiUEConeSumPt->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fhPhiUEConeSumPt) ; 
    
    fhEtaUEConeSumPt  = new TH2F("hEtaUEConeSumPt",
                                 Form("p_{T} in eta band around isolation cone for R = %2.2f",r),
                                 nptbins,ptmin,ptmax,4*nptinconebins,ptinconemin,4*ptinconemax);
    fhEtaUEConeSumPt->SetYTitle("p_{T in band} (GeV/c)");
    fhEtaUEConeSumPt->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fhEtaUEConeSumPt) ; 
    
    fhEtaBand  = new TH2F("fhEtaBand",
                          Form("Eta/Phi of particle in Eta band isolation cone for R = %2.2f",r),
                          netabins,etamin,etamax,nphibins,phimin,phimax); 
    fhEtaBand->SetXTitle("#eta");
    fhEtaBand->SetYTitle("#phi");
    outputContainer->Add(fhEtaBand) ; 
    
    fhPhiBand  = new TH2F("fhPhiBand",
                          Form("Eta/Phi of particle in Phi band isolation cone for R = %2.2f",r),
                          netabins,etamin,etamax,nphibins,phimin,phimax); 
    fhPhiBand->SetXTitle("#eta");
    fhPhiBand->SetYTitle("#phi");
    outputContainer->Add(fhPhiBand) ;
    
    fhConeSumPtEtaUESub  = new TH2F("hConeSumPtEtaUESub",
                                    Form("#Sigma p_{T} after bkg subtraction from eta band in the isolation cone for R = %2.2f",r),
                                    nptbins,ptmin,ptmax,2*nptsumbins,-ptsummax,ptsummax);
    fhConeSumPtEtaUESub->SetYTitle("#Sigma p_{T}");
    fhConeSumPtEtaUESub->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fhConeSumPtEtaUESub) ;
    
    fhConeSumPtPhiUESub  = new TH2F("hConeSumPtPhiUESub",
                                    Form("#Sigma p_{T} after bkg subtraction from phi band in the isolation cone for R = %2.2f",r),
                                    nptbins,ptmin,ptmax,2*nptsumbins,-ptsummax,ptsummax);
    fhConeSumPtPhiUESub->SetYTitle("#Sigma p_{T}");
    fhConeSumPtPhiUESub->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fhConeSumPtPhiUESub) ;
    
    fhEIso   = new TH1F("hE",
                        Form("Number of isolated particles vs E for R = %2.2f, p_{T}^{th} = %2.2f, p_{T}^{fr} = %2.2f",r,ptthre,ptfrac),
                        nptbins,ptmin,ptmax); 
    fhEIso->SetYTitle("dN / dE");
    fhEIso->SetXTitle("E (GeV/c)");
    outputContainer->Add(fhEIso) ; 
    
    fhPtIso  = new TH1F("hPt",
                        Form("Number of isolated particles vs p_{T} for R = %2.2f, p_{T}^{th} = %2.2f, p_{T}^{fr} = %2.2f",r,ptthre,ptfrac),
                        nptbins,ptmin,ptmax); 
    fhPtIso->SetYTitle("dN / p_{T}");
    fhPtIso->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fhPtIso) ; 
    
    fhPhiIso  = new TH2F("hPhi",
                         Form("Number of isolated particles vs #phi for R = %2.2f, p_{T}^{th} = %2.2f, p_{T}^{fr} = %2.2f",r,ptthre,ptfrac),
                         nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
    fhPhiIso->SetYTitle("#phi");
    fhPhiIso->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fhPhiIso) ; 
    
    fhEtaIso  = new TH2F("hEta",
                         Form("Number of isolated particles vs #eta for R = %2.2f, p_{T}^{th} = %2.2f, p_{T}^{fr} = %2.2f",r,ptthre,ptfrac),
                         nptbins,ptmin,ptmax,netabins,etamin,etamax); 
    fhEtaIso->SetYTitle("#eta");
    fhEtaIso->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fhEtaIso) ;
    
    fhEtaPhiIso  = new TH2F("hEtaPhiIso",
                            Form("Number of isolated particles #eta vs #phi for R = %2.2f, p_{T}^{th} = %2.2f, p_{T}^{fr} = %2.2f",r,ptthre,ptfrac),
                            netabins,etamin,etamax,nphibins,phimin,phimax); 
    fhEtaPhiIso->SetXTitle("#eta");
    fhEtaPhiIso->SetYTitle("#phi");
    outputContainer->Add(fhEtaPhiIso) ;
    
    fhPtDecayIso  = new TH1F("hPtDecayIso",
                             Form("Number of isolated #pi^{0} decay particles vs p_{T} for R = %2.2f, p_{T}^{th} = %2.2f, p_{T}^{fr} = %2.2f",r,ptthre,ptfrac),
                             nptbins,ptmin,ptmax); 
    fhPtDecayIso->SetYTitle("N");
    fhPtDecayIso->SetXTitle("p_{T}(GeV/c)");
    outputContainer->Add(fhPtDecayIso) ;
    
    fhEtaPhiDecayIso  = new TH2F("hEtaPhiDecayIso",
                                 Form("Number of isolated Pi0 decay particles #eta vs #phi for R = %2.2f, p_{T}^{th} = %2.2f, p_{T}^{fr} = %2.2f",r,ptthre,ptfrac),
                                 netabins,etamin,etamax,nphibins,phimin,phimax); 
    fhEtaPhiDecayIso->SetXTitle("#eta");
    fhEtaPhiDecayIso->SetYTitle("#phi");
    outputContainer->Add(fhEtaPhiDecayIso) ;
    
    if(IsDataMC())
    {
      fhPtIsoPrompt  = new TH1F("hPtMCPrompt","Number of isolated prompt #gamma",nptbins,ptmin,ptmax); 
      fhPtIsoPrompt->SetYTitle("N");
      fhPtIsoPrompt->SetXTitle("p_{T #gamma}(GeV/c)");
      outputContainer->Add(fhPtIsoPrompt) ; 
      
      fhPhiIsoPrompt  = new TH2F
      ("hPhiMCPrompt","Number of isolated prompt #gamma",nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
      fhPhiIsoPrompt->SetYTitle("#phi");
      fhPhiIsoPrompt->SetXTitle("p_{T #gamma} (GeV/c)");
      outputContainer->Add(fhPhiIsoPrompt) ; 
      
      fhEtaIsoPrompt  = new TH2F
      ("hEtaMCPrompt","Number of isolated prompt #gamma",nptbins,ptmin,ptmax,netabins,etamin,etamax); 
      fhEtaIsoPrompt->SetYTitle("#eta");
      fhEtaIsoPrompt->SetXTitle("p_{T #gamma} (GeV/c)");
      outputContainer->Add(fhEtaIsoPrompt) ;
      
      fhPtIsoFragmentation  = new TH1F("hPtMCFragmentation","Number of isolated #gamma",nptbins,ptmin,ptmax); 
      fhPtIsoFragmentation->SetYTitle("N");
      fhPtIsoFragmentation->SetXTitle("p_{T #gamma}(GeV/c)");
      outputContainer->Add(fhPtIsoFragmentation) ; 
      
      fhPhiIsoFragmentation  = new TH2F
      ("hPhiMCFragmentation","Number of isolated fragmentation #gamma",nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
      fhPhiIsoFragmentation->SetYTitle("#phi");
      fhPhiIsoFragmentation->SetXTitle("p_{T #gamma} (GeV/c)");
      outputContainer->Add(fhPhiIsoFragmentation) ; 
      
      fhEtaIsoFragmentation  = new TH2F
      ("hEtaMCFragmentation","Number of isolated fragmentation #gamma",nptbins,ptmin,ptmax,netabins,etamin,etamax); 
      fhEtaIsoFragmentation->SetYTitle("#eta");
      fhEtaIsoFragmentation->SetXTitle("p_{T #gamma} (GeV/c)");
      outputContainer->Add(fhEtaIsoFragmentation) ;
      
      fhPtIsoPi0  = new TH1F("hPtMCPi0","Number of isolated #gamma from #pi^{0} (2 #gamma)",nptbins,ptmin,ptmax); 
      fhPtIsoPi0->SetYTitle("N");
      fhPtIsoPi0->SetXTitle("p_{T #gamma}(GeV/c)");
      outputContainer->Add(fhPtIsoPi0) ; 
      
      fhPhiIsoPi0  = new TH2F
      ("hPhiMCPi0","Number of isolated #gamma from #pi^{0} (2 #gamma)",nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
      fhPhiIsoPi0->SetYTitle("#phi");
      fhPhiIsoPi0->SetXTitle("p_{T #gamma} (GeV/c)");
      outputContainer->Add(fhPhiIsoPi0) ; 
      
      fhEtaIsoPi0  = new TH2F
      ("hEtaMCPi0","Number of isolated #gamma from #pi^{0} (2 #gamma)",nptbins,ptmin,ptmax,netabins,etamin,etamax); 
      fhEtaIsoPi0->SetYTitle("#eta");
      fhEtaIsoPi0->SetXTitle("p_{T #gamma} (GeV/c)");
      outputContainer->Add(fhEtaIsoPi0) ;
      
      fhPtIsoPi0Decay  = new TH1F("hPtMCPi0Decay","Number of isolated #gamma from #pi^{0} decay",nptbins,ptmin,ptmax); 
      fhPtIsoPi0Decay->SetYTitle("N");
      fhPtIsoPi0Decay->SetXTitle("p_{T #gamma}(GeV/c)");
      outputContainer->Add(fhPtIsoPi0Decay) ; 
      
      fhPhiIsoPi0Decay  = new TH2F
      ("hPhiMCPi0Decay","Number of isolated #gamma from #pi^{0} decay",nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
      fhPhiIsoPi0Decay->SetYTitle("#phi");
      fhPhiIsoPi0Decay->SetXTitle("p_{T #gamma} (GeV/c)");
      outputContainer->Add(fhPhiIsoPi0Decay) ; 
      
      fhEtaIsoPi0Decay  = new TH2F
      ("hEtaMCPi0Decay","Number of isolated #gamma from #pi^{0} decay",nptbins,ptmin,ptmax,netabins,etamin,etamax); 
      fhEtaIsoPi0Decay->SetYTitle("#eta");
      fhEtaIsoPi0Decay->SetXTitle("p_{T #gamma} (GeV/c)");
      outputContainer->Add(fhEtaIsoPi0Decay) ;
      
      fhPtIsoEtaDecay  = new TH1F("hPtMCEtaDecay","Number of isolated #gamma from #eta decay",nptbins,ptmin,ptmax); 
      fhPtIsoEtaDecay->SetYTitle("N");
      fhPtIsoEtaDecay->SetXTitle("p_{T #gamma}(GeV/c)");
      outputContainer->Add(fhPtIsoEtaDecay) ; 
      
      fhPhiIsoEtaDecay  = new TH2F
      ("hPhiMCEtaDecay","Number of isolated #gamma from #eta decay",nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
      fhPhiIsoEtaDecay->SetYTitle("#phi");
      fhPhiIsoEtaDecay->SetXTitle("p_{T #gamma} (GeV/c)");
      outputContainer->Add(fhPhiIsoEtaDecay) ; 
      
      fhEtaIsoEtaDecay  = new TH2F
      ("hEtaMCEtaDecay","Number of isolated #gamma from #eta decay",nptbins,ptmin,ptmax,netabins,etamin,etamax); 
      fhEtaIsoEtaDecay->SetYTitle("#eta");
      fhEtaIsoEtaDecay->SetXTitle("p_{T #gamma} (GeV/c)");
      outputContainer->Add(fhEtaIsoEtaDecay) ;
      
      fhPtIsoOtherDecay  = new TH1F("hPtMCOtherDecay","Number of isolated #gamma from non-#pi^{0} decay",nptbins,ptmin,ptmax); 
      fhPtIsoOtherDecay->SetYTitle("N");
      fhPtIsoOtherDecay->SetXTitle("p_{T #gamma}(GeV/c)");
      outputContainer->Add(fhPtIsoOtherDecay) ; 
      
      fhPhiIsoOtherDecay  = new TH2F
      ("hPhiMCOtherDecay","Number of isolated #gamma from non-#pi^{0} decay",nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
      fhPhiIsoOtherDecay->SetYTitle("#phi");
      fhPhiIsoOtherDecay->SetXTitle("p_{T #gamma} (GeV/c)");
      outputContainer->Add(fhPhiIsoOtherDecay) ; 
      
      fhEtaIsoOtherDecay  = new TH2F
      ("hEtaMCOtherDecay","Number of isolated #gamma non-#pi^{0} decay",nptbins,ptmin,ptmax,netabins,etamin,etamax); 
      fhEtaIsoOtherDecay->SetYTitle("#eta");
      fhEtaIsoOtherDecay->SetXTitle("p_{T #gamma} (GeV/c)");
      outputContainer->Add(fhEtaIsoOtherDecay) ;
      
      //      fhPtIsoConversion  = new TH1F("hPtMCConversion","Number of isolated converted #gamma",nptbins,ptmin,ptmax); 
      //      fhPtIsoConversion->SetYTitle("N");
      //      fhPtIsoConversion->SetXTitle("p_{T #gamma}(GeV/c)");
      //      outputContainer->Add(fhPtIsoConversion) ; 
      //      
      //      fhPhiIsoConversion  = new TH2F
      //      ("hPhiMCConversion","Number of isolated converted #gamma",nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
      //      fhPhiIsoConversion->SetYTitle("#phi");
      //      fhPhiIsoConversion->SetXTitle("p_{T #gamma} (GeV/c)");
      //      outputContainer->Add(fhPhiIsoConversion) ; 
      //      
      //      fhEtaIsoConversion  = new TH2F
      //      ("hEtaMCConversion","Number of isolated converted #gamma",nptbins,ptmin,ptmax,netabins,etamin,etamax); 
      //      fhEtaIsoConversion->SetYTitle("#eta");
      //      fhEtaIsoConversion->SetXTitle("p_{T #gamma} (GeV/c)");
      //      outputContainer->Add(fhEtaIsoConversion) ;
      
      fhPtIsoHadron  = new TH1F("hPtMCHadron","Number of isolated non-#gamma particles",nptbins,ptmin,ptmax); 
      fhPtIsoHadron->SetYTitle("N");
      fhPtIsoHadron->SetXTitle("p_{T}(GeV/c)");
      outputContainer->Add(fhPtIsoHadron) ; 
      
      fhPhiIsoHadron  = new TH2F
      ("hPhiMCHadron","Number of isolated non-#gamma particles",nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
      fhPhiIsoHadron->SetYTitle("#phi");
      fhPhiIsoHadron->SetXTitle("p_{T} (GeV/c)");
      outputContainer->Add(fhPhiIsoHadron) ; 
      
      fhEtaIsoHadron  = new TH2F
      ("hEtaMCHadron","Number of isolated non-#gamma particles",nptbins,ptmin,ptmax,netabins,etamin,etamax); 
      fhEtaIsoHadron->SetYTitle("#eta");
      fhEtaIsoHadron->SetXTitle("p_{T} (GeV/c)");
      outputContainer->Add(fhEtaIsoHadron) ;
      
    }//Histos with MC
    
  }
  
  // Not Isolated histograms, reference histograms
  
  fhPtNoIso  = new TH1F("hPtNoIso",
                        Form("Number of not isolated leading particles vs p_{T} for R = %2.2f, p_{T}^{th} = %2.2f, p_{T}^{fr} = %2.2f",r,ptthre,ptfrac),
                        nptbins,ptmin,ptmax); 
  fhPtNoIso->SetYTitle("N");
  fhPtNoIso->SetXTitle("p_{T}(GeV/c)");
  outputContainer->Add(fhPtNoIso) ;
  
  
  fhEtaPhiNoIso  = new TH2F("hEtaPhiNoIso",
                            Form("Number of not isolated leading particles #eta vs #phi for R = %2.2f, p_{T}^{th} = %2.2f, p_{T}^{fr} = %2.2f",r,ptthre,ptfrac),
                            netabins,etamin,etamax,nphibins,phimin,phimax); 
  fhEtaPhiNoIso->SetXTitle("#eta");
  fhEtaPhiNoIso->SetYTitle("#phi");
  outputContainer->Add(fhEtaPhiNoIso) ;    
  
  fhPtDecayNoIso  = new TH1F("hPtDecayNoIso",
                             Form("Number of not isolated leading pi0 decay particles vs p_{T} for R = %2.2f, p_{T}^{th} = %2.2f, p_{T}^{fr} = %2.2f",r,ptthre,ptfrac),
                             nptbins,ptmin,ptmax); 
  fhPtDecayNoIso->SetYTitle("N");
  fhPtDecayNoIso->SetXTitle("p_{T}(GeV/c)");
  outputContainer->Add(fhPtDecayNoIso) ;
  
  fhEtaPhiDecayNoIso  = new TH2F("hEtaPhiDecayNoIso",
                                 Form("Number of not isolated leading Pi0 decay particles #eta vs #phi for R = %2.2f, p_{T}^{th} = %2.2f, p_{T}^{fr} = %2.2f",r,ptthre,ptfrac),
                                 netabins,etamin,etamax,nphibins,phimin,phimax); 
  fhEtaPhiDecayNoIso->SetXTitle("#eta");
  fhEtaPhiDecayNoIso->SetYTitle("#phi");
  outputContainer->Add(fhEtaPhiDecayNoIso) ;
  
 

  if(IsDataMC())
  {
    fhPtNoIsoPi0  = new TH1F
    ("hPtNoIsoPi0","Number of not isolated leading #gamma from #pi^{0} (2 #gamma)",nptbins,ptmin,ptmax); 
    fhPtNoIsoPi0->SetYTitle("N");
    fhPtNoIsoPi0->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fhPtNoIsoPi0) ;
    
    fhPtNoIsoPi0Decay  = new TH1F
    ("hPtNoIsoPi0Decay","Number of not isolated leading #gamma from #pi^{0} decay",nptbins,ptmin,ptmax); 
    fhPtNoIsoPi0Decay->SetYTitle("N");
    fhPtNoIsoPi0Decay->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fhPtNoIsoPi0Decay) ;
    
    fhPtNoIsoEtaDecay  = new TH1F
    ("hPtNoIsoEtaDecay","Number of not isolated leading #gamma from eta decay",nptbins,ptmin,ptmax); 
    fhPtNoIsoEtaDecay->SetYTitle("N");
    fhPtNoIsoEtaDecay->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fhPtNoIsoEtaDecay) ;
    
    fhPtNoIsoOtherDecay  = new TH1F
    ("hPtNoIsoOtherDecay","Number of not isolated leading #gamma from other decay",nptbins,ptmin,ptmax); 
    fhPtNoIsoOtherDecay->SetYTitle("N");
    fhPtNoIsoOtherDecay->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fhPtNoIsoOtherDecay) ;
    
    fhPtNoIsoPrompt  = new TH1F
    ("hPtNoIsoPrompt","Number of not isolated leading prompt #gamma",nptbins,ptmin,ptmax); 
    fhPtNoIsoPrompt->SetYTitle("N");
    fhPtNoIsoPrompt->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fhPtNoIsoPrompt) ;
    
    fhPtIsoMCPhoton  = new TH1F
    ("hPtIsoMCPhoton","Number of isolated leading  #gamma",nptbins,ptmin,ptmax); 
    fhPtIsoMCPhoton->SetYTitle("N");
    fhPtIsoMCPhoton->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fhPtIsoMCPhoton) ;
    
    fhPtNoIsoMCPhoton  = new TH1F
    ("hPtNoIsoMCPhoton","Number of not isolated leading #gamma",nptbins,ptmin,ptmax); 
    fhPtNoIsoMCPhoton->SetYTitle("N");
    fhPtNoIsoMCPhoton->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fhPtNoIsoMCPhoton) ;
    
//    fhPtNoIsoConversion  = new TH1F
//    ("hPtNoIsoConversion","Number of not isolated leading conversion #gamma",nptbins,ptmin,ptmax); 
//    fhPtNoIsoConversion->SetYTitle("N");
//    fhPtNoIsoConversion->SetXTitle("p_{T} (GeV/c)");
//    outputContainer->Add(fhPtNoIsoConversion) ;
    
    fhPtNoIsoFragmentation  = new TH1F
    ("hPtNoIsoFragmentation","Number of not isolated leading fragmentation #gamma",nptbins,ptmin,ptmax); 
    fhPtNoIsoFragmentation->SetYTitle("N");
    fhPtNoIsoFragmentation->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fhPtNoIsoFragmentation) ;
    
    fhPtNoIsoHadron  = new TH1F
    ("hPtNoIsoHadron","Number of not isolated leading hadrons",nptbins,ptmin,ptmax); 
    fhPtNoIsoHadron->SetYTitle("N");
    fhPtNoIsoHadron->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fhPtNoIsoHadron) ;
    
  }//Histos with MC
  
  
  if(fMakeSeveralIC)
  {
    const Int_t buffersize = 255;
    char name[buffersize];
    char title[buffersize];
    for(Int_t icone = 0; icone<fNCones; icone++)
    {	
      // sum pt in cone vs. pt leading
      snprintf(name, buffersize,"hSumPtLeadingPt_Cone_%d",icone);
      snprintf(title, buffersize,"#Sigma p_{T} in isolation cone for R = %2.2f",fConeSizes[icone]);
      fhSumPtLeadingPt[icone]  = new TH2F(name, title,nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
      fhSumPtLeadingPt[icone] ->SetYTitle("#sum_{cone}p_{T} (GeV/c)");//#Sigma p_{T}
      fhSumPtLeadingPt[icone] ->SetXTitle("p_{T}^{leading} (GeV/c)");
      outputContainer->Add(fhSumPtLeadingPt[icone]) ;
   
      // pt in cone vs. pt leading      
      snprintf(name, buffersize,"hPtLeadingPt_Cone_%d",icone);
      snprintf(title, buffersize,"p_{T} in isolation cone for R = %2.2f",fConeSizes[icone]);
      fhPtLeadingPt[icone]  = new TH2F(name, title,  nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax); 
      fhPtLeadingPt[icone] ->SetYTitle("p_{T}^{cone} (GeV/c)");
      fhPtLeadingPt[icone] ->SetXTitle("p_{T}^{leading} (GeV/c)");
      outputContainer->Add(fhPtLeadingPt[icone]) ;    

       // sum pt in cone vs. pt leading in the forward region (for background subtraction studies)
        snprintf(name, buffersize,"hFRSumPtLeadingPt_Cone_%d",icone);
      snprintf(title, buffersize,"#Sigma p_{T} in isolation cone for R = %2.2f",fConeSizes[icone]);
      fhFRSumPtLeadingPt[icone]  = new TH2F(name, title,nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
      fhFRSumPtLeadingPt[icone] ->SetYTitle("#sum_{cone}p_{T} (GeV/c)");//#Sigma p_{T}
      fhFRSumPtLeadingPt[icone] ->SetXTitle("p_{T}^{leading} (GeV/c)");
      outputContainer->Add(fhFRSumPtLeadingPt[icone]) ;
      
      // pt in cone vs. pt leading in the forward region (for background subtraction studies)    
        snprintf(name, buffersize,"hFRPtLeadingPt_Cone_%d",icone);
      snprintf(title, buffersize,"p_{T} in isolation cone for R = %2.2f",fConeSizes[icone]);
      fhFRPtLeadingPt[icone]  = new TH2F(name, title,  nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax); 
      fhFRPtLeadingPt[icone] ->SetYTitle("p_{T}^{cone} (GeV/c)");
      fhFRPtLeadingPt[icone] ->SetXTitle("p_{T}^{leading} (GeV/c)");
      outputContainer->Add(fhFRPtLeadingPt[icone]) ;    
  

      if(IsDataMC())
      {
        snprintf(name, buffersize,"hPtSumPrompt_Cone_%d",icone);
        snprintf(title, buffersize,"Candidate Prompt cone sum p_{T} for R = %2.2f vs candidate p_{T}",fConeSizes[icone]);
        fhPtSumIsolatedPrompt[icone]  = new TH2F(name, title,nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
        fhPtSumIsolatedPrompt[icone]->SetYTitle("#Sigma p_{T} (GeV/c)");
        fhPtSumIsolatedPrompt[icone]->SetXTitle("p_{T} (GeV/c)");
        outputContainer->Add(fhPtSumIsolatedPrompt[icone]) ; 
        
        snprintf(name, buffersize,"hPtSumFragmentation_Cone_%d",icone);
        snprintf(title, buffersize,"Candidate Fragmentation cone sum p_{T} for R = %2.2fvs candidate p_{T}",fConeSizes[icone]);
        fhPtSumIsolatedFragmentation[icone]  = new TH2F(name, title,nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
        fhPtSumIsolatedFragmentation[icone]->SetYTitle("#Sigma p_{T} (GeV/c)");
        fhPtSumIsolatedFragmentation[icone]->SetXTitle("p_{T} (GeV/c)");
        outputContainer->Add(fhPtSumIsolatedFragmentation[icone]) ; 
        
        snprintf(name, buffersize,"hPtSumPi0_Cone_%d",icone);
        snprintf(title, buffersize,"Candidate Pi0 cone sum p_{T} for R = %2.2f vs candidate p_{T}",fConeSizes[icone]);
        fhPtSumIsolatedPi0[icone]  = new TH2F(name, title,nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
        fhPtSumIsolatedPi0[icone]->SetYTitle("#Sigma p_{T} (GeV/c)");
        fhPtSumIsolatedPi0[icone]->SetXTitle("p_{T} (GeV/c)");
        outputContainer->Add(fhPtSumIsolatedPi0[icone]) ; 
        
        snprintf(name, buffersize,"hPtSumPi0Decay_Cone_%d",icone);
        snprintf(title, buffersize,"Candidate Pi0Decay cone sum p_{T} for R = %2.2f vs candidate p_{T}",fConeSizes[icone]);
        fhPtSumIsolatedPi0Decay[icone]  = new TH2F(name, title,nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
        fhPtSumIsolatedPi0Decay[icone]->SetYTitle("#Sigma p_{T} (GeV/c)");
        fhPtSumIsolatedPi0Decay[icone]->SetXTitle("p_{T} (GeV/c)");
        outputContainer->Add(fhPtSumIsolatedPi0Decay[icone]) ; 
        
        snprintf(name, buffersize,"hPtSumEtaDecay_Cone_%d",icone);
        snprintf(title, buffersize,"Candidate EtaDecay cone sum p_{T} for R = %2.2f vs candidate p_{T}",fConeSizes[icone]);
        fhPtSumIsolatedEtaDecay[icone]  = new TH2F(name, title,nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
        fhPtSumIsolatedEtaDecay[icone]->SetYTitle("#Sigma p_{T} (GeV/c)");
        fhPtSumIsolatedEtaDecay[icone]->SetXTitle("p_{T} (GeV/c)");
        outputContainer->Add(fhPtSumIsolatedEtaDecay[icone]) ;         
        
        snprintf(name, buffersize,"hPtSumOtherDecay_Cone_%d",icone);
        snprintf(title, buffersize,"Candidate OtherDecay cone sum p_{T} for R = %2.2f vs candidate p_{T}",fConeSizes[icone]);
        fhPtSumIsolatedOtherDecay[icone]  = new TH2F(name, title,nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
        fhPtSumIsolatedOtherDecay[icone]->SetYTitle("#Sigma p_{T} (GeV/c)");
        fhPtSumIsolatedOtherDecay[icone]->SetXTitle("p_{T} (GeV/c)");
        outputContainer->Add(fhPtSumIsolatedOtherDecay[icone]) ; 
        
//        snprintf(name, buffersize,"hPtSumConversion_Cone_%d",icone);
//        snprintf(title, buffersize,"Candidate Conversion cone sum p_{T} for R = %2.2f vs candidate p_{T}",fConeSizes[icone]);
//        fhPtSumIsolatedConversion[icone]  = new TH2F(name, title,nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
//        fhPtSumIsolatedConversion[icone]->SetYTitle("#Sigma p_{T} (GeV/c)");
//        fhPtSumIsolatedConversion[icone]->SetXTitle("p_{T} (GeV/c)");
//        outputContainer->Add(fhPtSumIsolatedConversion[icone]) ; 
        
        snprintf(name, buffersize,"hPtSumHadron_Cone_%d",icone);
        snprintf(title, buffersize,"Candidate Hadron cone sum p_{T} for R = %2.2f vs candidate p_{T}",fConeSizes[icone]);
        fhPtSumIsolatedHadron[icone]  = new TH2F(name, title,nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
        fhPtSumIsolatedHadron[icone]->SetYTitle("#Sigma p_{T} (GeV/c)");
        fhPtSumIsolatedHadron[icone]->SetXTitle("p_{T} (GeV/c)");
        outputContainer->Add(fhPtSumIsolatedHadron[icone]) ; 
        
      }//Histos with MC
      
      for(Int_t ipt = 0; ipt<fNPtThresFrac;ipt++)
      {   

        snprintf(name, buffersize,"hPtThres_Cone_%d_Pt%d",icone,ipt);
        snprintf(title, buffersize,"Isolated candidate p_{T} distribution for R = %2.2f and p_{T}^{th} = %2.2f GeV/c",fConeSizes[icone],fPtThresholds[ipt]);
        fhPtThresIsolated[icone][ipt]  = new TH1F(name, title,nptbins,ptmin,ptmax);
        fhPtThresIsolated[icone][ipt]->SetXTitle("p_{T} (GeV/c)");
        outputContainer->Add(fhPtThresIsolated[icone][ipt]) ; 
        
        snprintf(name, buffersize,"hPtFrac_Cone_%d_Pt%d",icone,ipt);
        snprintf(title, buffersize,"Isolated candidate p_{T} distribution for R = %2.2f and p_{T}^{fr} = %2.2f GeV/c",fConeSizes[icone],fPtFractions[ipt]);
        fhPtFracIsolated[icone][ipt]  = new TH1F(name, title,nptbins,ptmin,ptmax);
        fhPtFracIsolated[icone][ipt]->SetXTitle("p_{T} (GeV/c)");
        outputContainer->Add(fhPtFracIsolated[icone][ipt]) ; 
        
        
        snprintf(name, buffersize,"hPtSum_Cone_%d_Pt%d",icone,ipt);
        snprintf(title, buffersize,"Isolated candidate p_{T} distribution for R = %2.2f and p_{T}^{sum} = %2.2f GeV/c",fConeSizes[icone],fSumPtThresholds[ipt]);
        fhPtSumIsolated[icone][ipt]  = new TH1F(name, title,nptbins,ptmin,ptmax);
        // fhPtSumIsolated[icone][ipt]->SetYTitle("#Sigma p_{T} (GeV/c)");
        fhPtSumIsolated[icone][ipt]->SetXTitle("p_{T} (GeV/c)");
        outputContainer->Add(fhPtSumIsolated[icone][ipt]) ;
        
        snprintf(name, buffersize,"hPtSumDensity_Cone_%d_Pt%d",icone,ipt);
        snprintf(title, buffersize,"Isolated candidate p_{T} distribution for density in R = %2.2f and p_{T}^{sum} = %2.2f GeV/c",fConeSizes[icone],fSumPtThresholds[ipt]);
        fhPtSumDensityIso[icone][ipt]  = new TH1F(name, title,nptbins,ptmin,ptmax);//,nptsumbins,ptsummin,ptsummax);
        //fhPtSumIsolated[icone][ipt]->SetYTitle("#Sigma p_{T} (GeV/c)");
        fhPtSumDensityIso[icone][ipt]->SetXTitle("p_{T} (GeV/c)");
        outputContainer->Add(fhPtSumDensityIso[icone][ipt]) ;
        
        snprintf(name, buffersize,"hPtFracPtSum_Cone_%d_Pt%d",icone,ipt);
        snprintf(title, buffersize,"Isolated candidate p_{T} distribution for PtFracPtSum in R = %2.2f and p_{T}^{fr} = %2.2f GeV/c",fConeSizes[icone],fPtFractions[ipt]);
        fhPtFracPtSumIso[icone][ipt]  = new TH1F(name, title,nptbins,ptmin,ptmax);//,nptsumbins,ptsummin,ptsummax);
        //fhPtSumIsolated[icone][ipt]->SetYTitle("#Sigma p_{T} (GeV/c)");
        fhPtFracPtSumIso[icone][ipt]->SetXTitle("p_{T} (GeV/c)");
        outputContainer->Add(fhPtFracPtSumIso[icone][ipt]) ;
        
        // pt decays isolated
        snprintf(name, buffersize,"hPtThres_Decay_Cone_%d_Pt%d",icone,ipt);
        snprintf(title, buffersize,"Isolated decay candidate p_{T} distribution for R = %2.2f and p_{T}^{th} = %2.2f GeV/c",fConeSizes[icone],fPtThresholds[ipt]);
        fhPtPtThresDecayIso[icone][ipt]  = new TH1F(name, title,nptbins,ptmin,ptmax);
        fhPtPtThresDecayIso[icone][ipt]->SetXTitle("p_{T} (GeV/c)");
        outputContainer->Add(fhPtPtThresDecayIso[icone][ipt]) ;
        
        snprintf(name, buffersize,"hPtFrac_Decay_Cone_%d_Pt%d",icone,ipt);
        snprintf(title, buffersize,"Isolated decay candidate p_{T} distribution for R = %2.2f and p_{T}^{fr} = %2.2f GeV/c",fConeSizes[icone],fPtFractions[ipt]);
        fhPtPtFracDecayIso[icone][ipt]  = new TH1F(name, title,nptbins,ptmin,ptmax);
        fhPtPtFracDecayIso[icone][ipt]->SetXTitle("p_{T} (GeV/c)");
        outputContainer->Add(fhPtPtFracDecayIso[icone][ipt]) ;
        
        snprintf(name, buffersize,"hPtSum_Decay_Cone_%d_Pt%d",icone,ipt);
        snprintf(title, buffersize,"Isolated decay candidate p_{T} distribution for R = %2.2f and p_{T}^{sum} = %2.2f GeV/c",fConeSizes[icone],fSumPtThresholds[ipt]);
        fhPtPtSumDecayIso[icone][ipt]  = new TH1F(name, title,nptbins,ptmin,ptmax);//,nptsumbins,ptsummin,ptsummax);
        //  fhPtPtSumDecayIso[icone]->SetYTitle("#Sigma p_{T} (GeV/c)");
        fhPtPtSumDecayIso[icone][ipt]->SetXTitle("p_{T} (GeV/c)");
        outputContainer->Add(fhPtPtSumDecayIso[icone][ipt]) ;
        
        snprintf(name, buffersize,"hPtSumDensity_Decay_Cone_%d_Pt%d",icone,ipt);
        snprintf(title, buffersize,"Isolated decay candidate p_{T} distribution for density in R = %2.2f and p_{T}^{sum} = %2.2f GeV/c",fConeSizes[icone],fSumPtThresholds[ipt]);
        fhPtSumDensityDecayIso[icone][ipt]  = new TH1F(name, title,nptbins,ptmin,ptmax);//,nptsumbins,ptsummin,ptsummax);
        //  fhPtPtSumDecayIso[icone]->SetYTitle("#Sigma p_{T} (GeV/c)");
        fhPtSumDensityDecayIso[icone][ipt]->SetXTitle("p_{T} (GeV/c)");
        outputContainer->Add(fhPtSumDensityDecayIso[icone][ipt]) ;
        
        snprintf(name, buffersize,"hPtFracPtSum_Decay_Cone_%d_Pt%d",icone,ipt);
        snprintf(title, buffersize,"Isolated decay candidate p_{T} distribution for PtFracPtSum in R = %2.2f and p_{T}^{fr} = %2.2f GeV/c",fConeSizes[icone],fPtFractions[ipt]);
        fhPtFracPtSumDecayIso[icone][ipt]  = new TH1F(name, title,nptbins,ptmin,ptmax);//,nptsumbins,ptsummin,ptsummax);
        //  fhPtPtSumDecayIso[icone]->SetYTitle("#Sigma p_{T} (GeV/c)");
        fhPtFracPtSumDecayIso[icone][ipt]->SetXTitle("p_{T} (GeV/c)");
        outputContainer->Add(fhPtFracPtSumDecayIso[icone][ipt]) ;
        
        
        // eta:phi
        snprintf(name, buffersize,"hEtaPhiPtThres_Cone_%d_Pt%d",icone,ipt);
        snprintf(title, buffersize,"Isolated candidate #eta:#phi distribution for R = %2.2f and p_{T}^{th} = %2.2f GeV/c",fConeSizes[icone],fPtThresholds[ipt]);
        fhEtaPhiPtThresIso[icone][ipt]  = new TH2F(name, title,netabins,etamin,etamax,nphibins,phimin,phimax);
        fhEtaPhiPtThresIso[icone][ipt]->SetXTitle("#eta");
        fhEtaPhiPtThresIso[icone][ipt]->SetYTitle("#phi");
        outputContainer->Add(fhEtaPhiPtThresIso[icone][ipt]) ;
        
        snprintf(name, buffersize,"hEtaPhiPtFrac_Cone_%d_Pt%d",icone,ipt);
        snprintf(title, buffersize,"Isolated candidate #eta:#phi distribution for R = %2.2f and p_{T}^{fr} = %2.2f GeV/c",fConeSizes[icone],fPtFractions[ipt]);
        fhEtaPhiPtFracIso[icone][ipt]  = new TH2F(name, title,netabins,etamin,etamax,nphibins,phimin,phimax);
        fhEtaPhiPtFracIso[icone][ipt]->SetXTitle("#eta");
        fhEtaPhiPtFracIso[icone][ipt]->SetYTitle("#phi");
        outputContainer->Add(fhEtaPhiPtFracIso[icone][ipt]) ;
        
        snprintf(name, buffersize,"hEtaPhiPtSum_Cone_%d_Pt%d",icone,ipt);
        snprintf(title, buffersize,"Isolated candidate #eta:#phi distribution for R = %2.2f and p_{T}^{sum} = %2.2f GeV/c",fConeSizes[icone],fSumPtThresholds[ipt]);
        fhEtaPhiPtSumIso[icone][ipt]  = new TH2F(name, title,netabins,etamin,etamax,nphibins,phimin,phimax);
        fhEtaPhiPtSumIso[icone][ipt]->SetXTitle("#eta");
        fhEtaPhiPtSumIso[icone][ipt]->SetYTitle("#phi");
        outputContainer->Add(fhEtaPhiPtSumIso[icone][ipt]) ;
        
        snprintf(name, buffersize,"hEtaPhiSumDensity_Cone_%d_Pt%d",icone,ipt);
        snprintf(title, buffersize,"Isolated candidate #eta:#phi distribution for density R = %2.2f and p_{T}^{sum} = %2.2f GeV/c",fConeSizes[icone],fSumPtThresholds[ipt]);
        fhEtaPhiSumDensityIso[icone][ipt]  = new TH2F(name, title,netabins,etamin,etamax,nphibins,phimin,phimax);
        fhEtaPhiSumDensityIso[icone][ipt]->SetXTitle("#eta");
        fhEtaPhiSumDensityIso[icone][ipt]->SetYTitle("#phi");
        outputContainer->Add(fhEtaPhiSumDensityIso[icone][ipt]) ;
        
        snprintf(name, buffersize,"hEtaPhiFracPtSum_Cone_%d_Pt%d",icone,ipt);
        snprintf(title, buffersize,"Isolated candidate #eta:#phi distribution for FracPtSum R = %2.2f and p_{T}^{fr} = %2.2f GeV/c",fConeSizes[icone],fPtFractions[ipt]);
        fhEtaPhiFracPtSumIso[icone][ipt]  = new TH2F(name, title,netabins,etamin,etamax,nphibins,phimin,phimax);
        fhEtaPhiFracPtSumIso[icone][ipt]->SetXTitle("#eta");
        fhEtaPhiFracPtSumIso[icone][ipt]->SetYTitle("#phi");
        outputContainer->Add(fhEtaPhiFracPtSumIso[icone][ipt]) ;
        
        // eta:phi decays
        snprintf(name, buffersize,"hEtaPhiPtThres_Decay_Cone_%d_Pt%d",icone,ipt);
        snprintf(title, buffersize,"Isolated decay candidate #eta:#phi distribution for R = %2.2f and p_{T}^{th} = %2.2f GeV/c",fConeSizes[icone],fPtThresholds[ipt]);
        fhEtaPhiPtThresDecayIso[icone][ipt]  = new TH2F(name, title,netabins,etamin,etamax,nphibins,phimin,phimax);
        fhEtaPhiPtThresDecayIso[icone][ipt]->SetXTitle("#eta");
        fhEtaPhiPtThresDecayIso[icone][ipt]->SetYTitle("#phi");
        outputContainer->Add(fhEtaPhiPtThresDecayIso[icone][ipt]) ;
        
        snprintf(name, buffersize,"hEtaPhiPtFrac_Decay_Cone_%d_Pt%d",icone,ipt);
        snprintf(title, buffersize,"Isolated decay candidate #eta:#phi distribution for R = %2.2f and p_{T}^{fr} = %2.2f GeV/c",fConeSizes[icone],fPtFractions[ipt]);
        fhEtaPhiPtFracDecayIso[icone][ipt]  = new TH2F(name, title,netabins,etamin,etamax,nphibins,phimin,phimax);
        fhEtaPhiPtFracDecayIso[icone][ipt]->SetXTitle("#eta");
        fhEtaPhiPtFracDecayIso[icone][ipt]->SetYTitle("#phi");
        outputContainer->Add(fhEtaPhiPtFracDecayIso[icone][ipt]) ;
        
        
        snprintf(name, buffersize,"hEtaPhiPtSum_Decay_Cone_%d_Pt%d",icone,ipt);
        snprintf(title, buffersize,"Isolated decay candidate #eta:#phi distribution for R = %2.2f and p_{T}^{sum} = %2.2f GeV/c",fConeSizes[icone],fSumPtThresholds[ipt]);
        fhEtaPhiPtSumDecayIso[icone][ipt]  = new TH2F(name, title,netabins,etamin,etamax,nphibins,phimin,phimax);
        fhEtaPhiPtSumDecayIso[icone][ipt]->SetXTitle("#eta");
        fhEtaPhiPtSumDecayIso[icone][ipt]->SetYTitle("#phi");
        outputContainer->Add(fhEtaPhiPtSumDecayIso[icone][ipt]) ;
        
        snprintf(name, buffersize,"hEtaPhiSumDensity_Decay_Cone_%d_Pt%d",icone,ipt);
        snprintf(title, buffersize,"Isolated decay candidate #eta:#phi distribution for density R = %2.2f and p_{T}^{sum} = %2.2f GeV/c",fConeSizes[icone],fSumPtThresholds[ipt]);
        fhEtaPhiSumDensityDecayIso[icone][ipt]  = new TH2F(name, title,netabins,etamin,etamax,nphibins,phimin,phimax);
        fhEtaPhiSumDensityDecayIso[icone][ipt]->SetXTitle("#eta");
        fhEtaPhiSumDensityDecayIso[icone][ipt]->SetYTitle("#phi");
        outputContainer->Add(fhEtaPhiSumDensityDecayIso[icone][ipt]) ;
        
        snprintf(name, buffersize,"hEtaPhiFracPtSum_Decay_Cone_%d_Pt%d",icone,ipt);
        snprintf(title, buffersize,"Isolated decay candidate #eta:#phi distribution for FracPtSum R = %2.2f and p_{T}^{fr} = %2.2f GeV/c",fConeSizes[icone],fPtFractions[ipt]);
        fhEtaPhiFracPtSumDecayIso[icone][ipt]  = new TH2F(name, title,netabins,etamin,etamax,nphibins,phimin,phimax);
        fhEtaPhiFracPtSumDecayIso[icone][ipt]->SetXTitle("#eta");
        fhEtaPhiFracPtSumDecayIso[icone][ipt]->SetYTitle("#phi");
        outputContainer->Add(fhEtaPhiFracPtSumDecayIso[icone][ipt]) ;
        
        
        if(IsDataMC())
        {
          snprintf(name, buffersize,"hPtThresMCPrompt_Cone_%d_Pt%d",icone,ipt);
          snprintf(title, buffersize,"Isolated candidate Prompt p_{T} distribution for cone size %d and p_{T}^{th} %d",icone,ipt);
          fhPtThresIsolatedPrompt[icone][ipt]  = new TH1F(name, title,nptbins,ptmin,ptmax);
          fhPtThresIsolatedPrompt[icone][ipt]->SetXTitle("p_{T} (GeV/c)");
          outputContainer->Add(fhPtThresIsolatedPrompt[icone][ipt]) ; 
          
          snprintf(name, buffersize,"hPtFracMCPrompt_Cone_%d_Pt%d",icone,ipt);
          snprintf(title, buffersize,"Isolated candidate Prompt p_{T} distribution for cone size %d and p_{T}^{th} %d",icone,ipt);
          fhPtFracIsolatedPrompt[icone][ipt]  = new TH1F(name, title,nptbins,ptmin,ptmax);
          fhPtFracIsolatedPrompt[icone][ipt]->SetXTitle("p_{T} (GeV/c)");
          outputContainer->Add(fhPtFracIsolatedPrompt[icone][ipt]) ; 
          
          snprintf(name, buffersize,"hPtThresMCFragmentation_Cone_%d_Pt%d",icone,ipt);
          snprintf(title, buffersize,"Isolated candidate Fragmentation p_{T} distribution for cone size %d and p_{T}^{th} %d",icone,ipt);
          fhPtThresIsolatedFragmentation[icone][ipt]  = new TH1F(name, title,nptbins,ptmin,ptmax);
          fhPtThresIsolatedFragmentation[icone][ipt]->SetXTitle("p_{T} (GeV/c)");
          outputContainer->Add(fhPtThresIsolatedFragmentation[icone][ipt]) ; 
          
          snprintf(name, buffersize,"hPtFracMCFragmentation_Cone_%d_Pt%d",icone,ipt);
          snprintf(title, buffersize,"Isolated candidate Fragmentation p_{T} distribution for cone size %d and p_{T}^{th} %d",icone,ipt);
          fhPtFracIsolatedFragmentation[icone][ipt]  = new TH1F(name, title,nptbins,ptmin,ptmax);
          fhPtFracIsolatedFragmentation[icone][ipt]->SetXTitle("p_{T} (GeV/c)");
          outputContainer->Add(fhPtFracIsolatedFragmentation[icone][ipt]) ; 
          
          snprintf(name, buffersize,"hPtThresMCPi0_Cone_%d_Pt%d",icone,ipt);
          snprintf(title, buffersize,"Isolated candidate Pi0 p_{T} distribution for cone size %d and p_{T}^{th} %d",icone,ipt);
          fhPtThresIsolatedPi0[icone][ipt]  = new TH1F(name, title,nptbins,ptmin,ptmax);
          fhPtThresIsolatedPi0[icone][ipt]->SetXTitle("p_{T} (GeV/c)");
          outputContainer->Add(fhPtThresIsolatedPi0[icone][ipt]) ; 
          
          snprintf(name, buffersize,"hPtFracMCPi0_Cone_%d_Pt%d",icone,ipt);
          snprintf(title, buffersize,"Isolated candidate Pi0 p_{T} distribution for cone size %d and p_{T}^{th} %d",icone,ipt);
          fhPtFracIsolatedPi0[icone][ipt]  = new TH1F(name, title,nptbins,ptmin,ptmax);
          fhPtFracIsolatedPi0[icone][ipt]->SetXTitle("p_{T} (GeV/c)");
          outputContainer->Add(fhPtFracIsolatedPi0[icone][ipt]) ;
          
          snprintf(name, buffersize,"hPtThresMCPi0Decay_Cone_%d_Pt%d",icone,ipt);
          snprintf(title, buffersize,"Isolated candidate Pi0Decay p_{T} distribution for cone size %d and p_{T}^{th} %d",icone,ipt);
          fhPtThresIsolatedPi0Decay[icone][ipt]  = new TH1F(name, title,nptbins,ptmin,ptmax);
          fhPtThresIsolatedPi0Decay[icone][ipt]->SetXTitle("p_{T} (GeV/c)");
          outputContainer->Add(fhPtThresIsolatedPi0Decay[icone][ipt]) ; 
          
          snprintf(name, buffersize,"hPtFracMCPi0Decay_Cone_%d_Pt%d",icone,ipt);
          snprintf(title, buffersize,"Isolated candidate Pi0Decay p_{T} distribution for cone size %d and p_{T}^{th} %d",icone,ipt);
          fhPtFracIsolatedPi0Decay[icone][ipt]  = new TH1F(name, title,nptbins,ptmin,ptmax);
          fhPtFracIsolatedPi0Decay[icone][ipt]->SetXTitle("p_{T} (GeV/c)");
          outputContainer->Add(fhPtFracIsolatedPi0Decay[icone][ipt]) ; 
          
          snprintf(name, buffersize,"hPtThresMCEtaDecay_Cone_%d_Pt%d",icone,ipt);
          snprintf(title, buffersize,"Isolated candidate EtaDecay p_{T} distribution for cone size %d and p_{T}^{th} %d",icone,ipt);
          fhPtThresIsolatedEtaDecay[icone][ipt]  = new TH1F(name, title,nptbins,ptmin,ptmax);
          fhPtThresIsolatedEtaDecay[icone][ipt]->SetXTitle("p_{T} (GeV/c)");
          outputContainer->Add(fhPtThresIsolatedEtaDecay[icone][ipt]) ; 
          
          snprintf(name, buffersize,"hPtFracMCEtaDecay_Cone_%d_Pt%d",icone,ipt);
          snprintf(title, buffersize,"Isolated candidate EtaDecay p_{T} distribution for cone size %d and p_{T}^{th} %d",icone,ipt);
          fhPtFracIsolatedEtaDecay[icone][ipt]  = new TH1F(name, title,nptbins,ptmin,ptmax);
          fhPtFracIsolatedEtaDecay[icone][ipt]->SetXTitle("p_{T} (GeV/c)");
          outputContainer->Add(fhPtFracIsolatedEtaDecay[icone][ipt]) ; 
          
          
          snprintf(name, buffersize,"hPtThresMCOtherDecay_Cone_%d_Pt%d",icone,ipt);
          snprintf(title, buffersize,"Isolated candidate OtherDecay p_{T} distribution for cone size %d and p_{T}^{th} %d",icone,ipt);
          fhPtThresIsolatedOtherDecay[icone][ipt]  = new TH1F(name, title,nptbins,ptmin,ptmax);
          fhPtThresIsolatedOtherDecay[icone][ipt]->SetXTitle("p_{T} (GeV/c)");
          outputContainer->Add(fhPtThresIsolatedOtherDecay[icone][ipt]) ; 
          
          snprintf(name, buffersize,"hPtFracMCOtherDecay_Cone_%d_Pt%d",icone,ipt);
          snprintf(title, buffersize,"Isolated candidate OtherDecay p_{T} distribution for cone size %d and p_{T}^{th} %d",icone,ipt);
          fhPtFracIsolatedOtherDecay[icone][ipt]  = new TH1F(name, title,nptbins,ptmin,ptmax);
          fhPtFracIsolatedOtherDecay[icone][ipt]->SetXTitle("p_{T} (GeV/c)");
          outputContainer->Add(fhPtFracIsolatedOtherDecay[icone][ipt]) ;
          
//          snprintf(name, buffersize,"hPtThresMCConversion_Cone_%d_Pt%d",icone,ipt);
//          snprintf(title, buffersize,"Isolated candidate Conversion p_{T} distribution for cone size %d and p_{T}^{th} %d",icone,ipt);
//          fhPtThresIsolatedConversion[icone][ipt]  = new TH1F(name, title,nptbins,ptmin,ptmax);
//          fhPtThresIsolatedConversion[icone][ipt]->SetXTitle("p_{T} (GeV/c)");
//          outputContainer->Add(fhPtThresIsolatedConversion[icone][ipt]) ; 
          
//          snprintf(name, buffersize,"hPtFracMCConversion_Cone_%d_Pt%d",icone,ipt);
//          snprintf(title, buffersize,"Isolated candidate Conversion p_{T} distribution for cone size %d and p_{T}^{th} %d",icone,ipt);
//          fhPtFracIsolatedConversion[icone][ipt]  = new TH1F(name, title,nptbins,ptmin,ptmax);
//          fhPtFracIsolatedConversion[icone][ipt]->SetXTitle("p_{T} (GeV/c)");
//          outputContainer->Add(fhPtFracIsolatedConversion[icone][ipt]) ;
          
          snprintf(name, buffersize,"hPtThresMCHadron_Cone_%d_Pt%d",icone,ipt);
          snprintf(title, buffersize,"Isolated candidate Hadron p_{T} distribution for cone size %d and p_{T}^{th} %d",icone,ipt);
          fhPtThresIsolatedHadron[icone][ipt]  = new TH1F(name, title,nptbins,ptmin,ptmax);
          fhPtThresIsolatedHadron[icone][ipt]->SetXTitle("p_{T} (GeV/c)");
          outputContainer->Add(fhPtThresIsolatedHadron[icone][ipt]) ; 
          
          snprintf(name, buffersize,"hPtFracMCHadron_Cone_%d_Pt%d",icone,ipt);
          snprintf(title, buffersize,"Isolated candidate Hadron p_{T} distribution for cone size %d and p_{T}^{th} %d",icone,ipt);
          fhPtFracIsolatedHadron[icone][ipt]  = new TH1F(name, title,nptbins,ptmin,ptmax);
          fhPtFracIsolatedHadron[icone][ipt]->SetXTitle("p_{T} (GeV/c)");
          outputContainer->Add(fhPtFracIsolatedHadron[icone][ipt]) ;  
          
        }//Histos with MC
      }//icone loop
    }//ipt loop
  }
  
  if(fFillPileUpHistograms)
  {
    fhTimeENoCut  = new TH2F ("hTimeE_NoCut","time of cluster vs E of clusters, no cut", nptbins,ptmin,ptmax, ntimebins,timemin,timemax); 
    fhTimeENoCut->SetXTitle("E (GeV)");
    fhTimeENoCut->SetYTitle("time (ns)");
    outputContainer->Add(fhTimeENoCut);  
    
    fhTimeESPD  = new TH2F ("hTimeE_SPD","time of cluster vs E of clusters, SPD cut", nptbins,ptmin,ptmax, ntimebins,timemin,timemax); 
    fhTimeESPD->SetXTitle("E (GeV)");
    fhTimeESPD->SetYTitle("time (ns)");
    outputContainer->Add(fhTimeESPD);  
    
    fhTimeESPDMulti  = new TH2F ("hTimeE_SPDMulti","time of cluster vs E of clusters, SPD multi cut", nptbins,ptmin,ptmax, ntimebins,timemin,timemax); 
    fhTimeESPDMulti->SetXTitle("E (GeV)");
    fhTimeESPDMulti->SetYTitle("time (ns)");
    outputContainer->Add(fhTimeESPDMulti);  
    
    fhTimeNPileUpVertSPD  = new TH2F ("hTime_NPileUpVertSPD","time of cluster vs N pile-up SPD vertex", ntimebins,timemin,timemax,50,0,50); 
    fhTimeNPileUpVertSPD->SetYTitle("# vertex ");
    fhTimeNPileUpVertSPD->SetXTitle("time (ns)");
    outputContainer->Add(fhTimeNPileUpVertSPD);  
    
    fhTimeNPileUpVertTrack  = new TH2F ("hTime_NPileUpVertTracks","time of cluster vs N pile-up Tracks vertex", ntimebins,timemin,timemax, 50,0,50 ); 
    fhTimeNPileUpVertTrack->SetYTitle("# vertex ");
    fhTimeNPileUpVertTrack->SetXTitle("time (ns)");
    outputContainer->Add(fhTimeNPileUpVertTrack);  
    
    fhTimeNPileUpVertContributors  = new TH2F ("hTime_NPileUpVertContributors","time of cluster vs N constributors to pile-up SPD vertex", ntimebins,timemin,timemax,50,0,50); 
    fhTimeNPileUpVertContributors->SetYTitle("# vertex ");
    fhTimeNPileUpVertContributors->SetXTitle("time (ns)");
    outputContainer->Add(fhTimeNPileUpVertContributors);  
    
    fhTimePileUpMainVertexZDistance  = new TH2F ("hTime_PileUpMainVertexZDistance","time of cluster vs distance in Z pile-up SPD vertex - main SPD vertex",ntimebins,timemin,timemax,100,0,50); 
    fhTimePileUpMainVertexZDistance->SetYTitle("distance Z (cm) ");
    fhTimePileUpMainVertexZDistance->SetXTitle("time (ns)");
    outputContainer->Add(fhTimePileUpMainVertexZDistance);  
    
    fhTimePileUpMainVertexZDiamond  = new TH2F ("hTime_PileUpMainVertexZDiamond","time of cluster vs distance in Z pile-up SPD vertex - z diamond",ntimebins,timemin,timemax,100,0,50); 
    fhTimePileUpMainVertexZDiamond->SetYTitle("diamond distance Z (cm) ");
    fhTimePileUpMainVertexZDiamond->SetXTitle("time (ns)");
    outputContainer->Add(fhTimePileUpMainVertexZDiamond);  
    
  }
  
  return outputContainer ;
  
}

//__________________________________
void AliAnaParticleIsolation::Init()
{
  // Do some checks and init stuff
  
  // In case of several cone and thresholds analysis, open the cuts for the filling of the 
  // track and cluster reference arrays in cone when done in the MakeAnalysisFillAOD(). 
  // The different cones, thresholds are tested for this list of tracks, clusters.
  if(fMakeSeveralIC)
  {
    printf("AliAnaParticleIsolation::Init() - Open default isolation cuts for multiple Isolation analysis\n");
    GetIsolationCut()->SetPtThreshold(100);
    GetIsolationCut()->SetPtFraction(100);
    GetIsolationCut()->SetConeSize(1);
  }
}

//____________________________________________
void AliAnaParticleIsolation::InitParameters()
{
  
  //Initialize the parameters of the analysis.
  SetInputAODName("PWG4Particle");
  SetAODObjArrayName("IsolationCone");  
  AddToHistogramsName("AnaIsolation_");
  
  fCalorimeter = "PHOS" ;
  fReMakeIC = kFALSE ;
  fMakeSeveralIC = kFALSE ;
  
  //----------- Several IC-----------------
  fNCones             = 5 ; 
  fNPtThresFrac       = 5 ; 
  fConeSizes      [0] = 0.1;     fConeSizes      [1] = 0.2;   fConeSizes      [2] = 0.3; fConeSizes      [3] = 0.4;  fConeSizes      [4] = 0.5;
  fPtThresholds   [0] = 1.;      fPtThresholds   [1] = 2.;    fPtThresholds   [2] = 3.;  fPtThresholds   [3] = 4.;   fPtThresholds   [4] = 5.; 
  fPtFractions    [0] = 0.05;    fPtFractions    [1] = 0.075; fPtFractions    [2] = 0.1; fPtFractions    [3] = 1.25; fPtFractions    [4] = 1.5; 
  fSumPtThresholds[0] = 1.;      fSumPtThresholds[1] = 2.;    fSumPtThresholds[2] = 3.;  fSumPtThresholds[3] = 4.;   fSumPtThresholds[4] = 5.; 
  
  //------------- Histograms settings -------
  fHistoNPtSumBins = 100 ;
  fHistoPtSumMax   = 50 ;
  fHistoPtSumMin   = 0.  ;
  
  fHistoNPtInConeBins = 100 ;
  fHistoPtInConeMax   = 50 ;
  fHistoPtInConeMin   = 0.  ;
  
}

//__________________________________________________
void  AliAnaParticleIsolation::MakeAnalysisFillAOD() 
{
  //Do analysis and fill aods
  //Search for the isolated photon in fCalorimeter with pt > GetMinPt()
  
  if(!GetInputAODBranch())
  {
    printf("AliAnaParticleIsolation::MakeAnalysisFillAOD() - No input particles in AOD with name branch < %s >, STOP \n",GetInputAODName().Data());
    abort();
  }
  
  if(strcmp(GetInputAODBranch()->GetClass()->GetName(), "AliAODPWG4ParticleCorrelation"))
  {
    printf("AliAnaParticleIsolation::MakeAnalysisFillAOD() - Wrong type of AOD object, change AOD class name in input AOD: It should be <AliAODPWG4ParticleCorrelation> and not <%s> \n",GetInputAODBranch()->GetClass()->GetName());
    abort();
  }
  
  Int_t n = 0, nfrac = 0;
  Bool_t  isolated  = kFALSE ;
  Float_t coneptsum = 0 ;
  TObjArray * pl    = 0x0; ; 
  
  //Select the calorimeter for candidate isolation with neutral particles
  if      (fCalorimeter == "PHOS" )
    pl = GetPHOSClusters();
  else if (fCalorimeter == "EMCAL")
    pl = GetEMCALClusters();
  
  //Loop on AOD branch, filled previously in AliAnaPhoton, find leading particle to do isolation only with it
  Double_t ptLeading = 0. ;
  Int_t    idLeading = -1 ;
  TLorentzVector mom ;
  Int_t naod = GetInputAODBranch()->GetEntriesFast();
  if(GetDebug() > 0) printf("AliAnaParticleIsolation::MakeAnalysisFillAOD() - Input aod branch entries %d\n", naod);
  
  for(Int_t iaod = 0; iaod < naod; iaod++)
  {
    AliAODPWG4ParticleCorrelation * aodinput =  (AliAODPWG4ParticleCorrelation*) (GetInputAODBranch()->At(iaod));
    
    //If too small or too large pt, skip
    if(aodinput->Pt() < GetMinPt() || aodinput->Pt() > GetMaxPt() ) continue ; 
    
    //check if it is low pt trigger particle
    if((aodinput->Pt() < GetIsolationCut()->GetPtThreshold() || 
        aodinput->Pt() < GetIsolationCut()->GetSumPtThreshold()) && 
       !fMakeSeveralIC)
    {
      continue ; //trigger should not come from underlying event
    }
    
    //vertex cut in case of mixing
    Int_t check = CheckMixedEventVertex(aodinput->GetCaloLabel(0), aodinput->GetTrackLabel(0));
    if(check ==  0) continue;
    if(check == -1) return;
    
    //find the leading particles with highest momentum
    if ( aodinput->Pt() > ptLeading ) 
    {
      ptLeading = aodinput->Pt() ;
      idLeading = iaod ;
    }
    
    aodinput->SetLeadingParticle(kFALSE);
    
  }//finish searching for leading trigger particle
  
  // Check isolation of leading particle
  if(idLeading < 0) return;
  
  AliAODPWG4ParticleCorrelation * aodinput =  (AliAODPWG4ParticleCorrelation*) (GetInputAODBranch()->At(idLeading));
  aodinput->SetLeadingParticle(kTRUE);
  
  // Check isolation only of clusters in fiducial region
  if(IsFiducialCutOn())
  {
    Bool_t in = GetFiducialCut()->IsInFiducialCut(*aodinput->Momentum(),aodinput->GetDetector()) ;
    if(! in ) return ;
  }
  
  //After cuts, study isolation
  n=0; nfrac = 0; isolated = kFALSE; coneptsum = 0;
  GetIsolationCut()->MakeIsolationCut(GetCTSTracks(),pl,
                                      GetReader(), GetCaloPID(),
                                      kTRUE, aodinput, GetAODObjArrayName(), 
                                      n,nfrac,coneptsum, isolated);
  
  if(!fMakeSeveralIC) aodinput->SetIsolated(isolated);
    
  if(GetDebug() > 1) 
  {
    if(isolated)printf("AliAnaParticleIsolation::MakeAnalysisFillAOD() : Particle %d IS ISOLATED \n",idLeading);
    printf("AliAnaParticleIsolation::MakeAnalysisFillAOD() - End fill AODs \n");  
  }
  
}

//_________________________________________________________
void  AliAnaParticleIsolation::MakeAnalysisFillHistograms() 
{
  //Do analysis and fill histograms
  
  Int_t   n = 0, nfrac = 0;
  Bool_t  isolated  = kFALSE ;
  Float_t coneptsum = 0 ;
  Float_t etaUEptsum = 0 ;
  Float_t phiUEptsum = 0 ;
  
  //Loop on stored AOD 
  Int_t naod = GetInputAODBranch()->GetEntriesFast();
  if(GetDebug() > 0) printf("AliAnaParticleIsolation::MakeAnalysisFillHistograms() - Histo aod branch entries %d\n", naod);
  
  //Get vertex for photon momentum calculation
  Double_t vertex[] = {0,0,0} ; //vertex ;
  if(GetReader()->GetDataType() != AliCaloTrackReader::kMC) 
  {
    GetReader()->GetVertex(vertex);
  }	
  
  for(Int_t iaod = 0; iaod < naod ; iaod++)
  {
    AliAODPWG4ParticleCorrelation* aod =  (AliAODPWG4ParticleCorrelation*) (GetInputAODBranch()->At(iaod));
    
    if(!aod->IsLeadingParticle()) continue; // Try to isolate only leading cluster or track
    
    // Check isolation only of clusters in fiducial region
    if(IsFiducialCutOn())
    {
      Bool_t in = GetFiducialCut()->IsInFiducialCut(*aod->Momentum(),aod->GetDetector()) ;
      if(! in ) continue ;
    }
    
    Bool_t  isolation  = aod->IsIsolated(); 
    Bool_t  decay      = aod->IsTagged();
    Float_t energy     = aod->E();
    Float_t pt         = aod->Pt();
    Float_t phi        = aod->Phi();
    Float_t eta        = aod->Eta();
    Float_t conesize   = GetIsolationCut()->GetConeSize();
    
    //Recover reference arrays with clusters and tracks
    TObjArray * refclusters = aod->GetObjArray(GetAODObjArrayName()+"Clusters");
    TObjArray * reftracks   = aod->GetObjArray(GetAODObjArrayName()+"Tracks");
    
    //If too small or too large pt, skip
    if(pt < GetMinPt() || pt > GetMaxPt() ) continue ; 
    
    // --- In case of redoing isolation from delta AOD ----
    
    if(fMakeSeveralIC) 
    {
      //Analysis of multiple IC at same time
      MakeSeveralICAnalysis(aod);
      continue;
    }
    else if(fReMakeIC)
    {
      //In case a more strict IC is needed in the produced AOD
      n=0; nfrac = 0; isolated = kFALSE; coneptsum = 0;
      GetIsolationCut()->MakeIsolationCut(reftracks,   refclusters, 
                                          GetReader(), GetCaloPID(),
                                          kFALSE, aod, "", 
                                          n,nfrac,coneptsum, isolated);
      fhConeSumPt->Fill(pt,coneptsum);    
      if(GetDebug() > 0) printf("AliAnaParticleIsolation::MakeAnalysisFillHistograms() - Energy Sum in Isolation Cone %2.2f\n", coneptsum);    
    }
    
    // ---------------------------------------------------
    
    //Fill pt distribution of particles in cone
    //Tracks
    coneptsum = 0;
    Double_t sumptFR = 0. ;
    TObjArray * trackList   = GetCTSTracks() ;
    for(Int_t itrack=0; itrack < trackList->GetEntriesFast(); itrack++)
    {
      AliVTrack* track = (AliVTrack *) trackList->At(itrack);
      
      if(!track)
      {
        printf("AliAnaParticleIsolation::MakeAnalysisFillHistograms() - Track not available?");
        continue;
      }
      
      //fill histogram for UE in phi band
      if(track->Eta() > (eta-conesize) && track->Eta()  < (eta+conesize))
      {
        phiUEptsum+=track->Pt();
        fhPhiBand->Fill(track->Eta(),track->Phi());
      }
      //fill histogram for UE in eta band in EMCal acceptance
      if(track->Phi() > (phi-conesize) && track->Phi() < (phi+conesize) && track->Eta() > -0.6 && track->Eta() < 0.6)
      {
        etaUEptsum+=track->Pt();
        fhEtaBand->Fill(track->Eta(),track->Phi());
      }
      
      //fill the histograms at forward range
      Double_t dPhi = phi - track->Phi() + TMath::PiOver2();
      Double_t dEta = eta - track->Eta();
      Double_t arg  = dPhi*dPhi + dEta*dEta;
      if(TMath::Sqrt(arg) < conesize)
      {
        fhPtInFRCone->Fill(pt,TMath::Sqrt(track->Px()*track->Px()+track->Py()*track->Py()));
        sumptFR+=track->Pt();
      }
      
      dPhi = phi - track->Phi() - TMath::PiOver2();
      arg  = dPhi*dPhi + dEta*dEta;
      if(TMath::Sqrt(arg) < conesize)
      {
        fhPtInFRCone->Fill(pt,TMath::Sqrt(track->Px()*track->Px()+track->Py()*track->Py()));
        sumptFR+=track->Pt();
      }      
    }
    
    fhFRConeSumPt->Fill(pt,sumptFR);
    if(reftracks)
    {  
      for(Int_t itrack=0; itrack < reftracks->GetEntriesFast(); itrack++)
      {
        AliVTrack* track = (AliVTrack *) reftracks->At(itrack);
        Float_t pTtrack = track->Pt();
        fhPtInCone->Fill(pt,pTtrack);
        if(fFillPileUpHistograms && GetReader()->IsPileUpFromSPD())  fhPtInConePileUp->Fill(pt,pTtrack);
        if (GetEventCentrality()) fhPtInConeCent->Fill(GetEventCentrality(),pTtrack);
        coneptsum+=pTtrack;
      }
    }
    
    //CaloClusters
    if(refclusters)
    {    
      TLorentzVector mom ;
      for(Int_t icalo=0; icalo < refclusters->GetEntriesFast(); icalo++)
      {
        AliVCluster* calo = (AliVCluster *) refclusters->At(icalo);
        calo->GetMomentum(mom,vertex) ;//Assume that come from vertex in straight line
        
        //fill histogram for UE in phi band
        if(mom.Eta() > (eta-conesize) && mom.Eta()  < (eta+conesize))
        {
          phiUEptsum+=mom.Pt();
          fhPhiBand->Fill(mom.Eta(),mom.Phi());
        }
        //fill histogram for UE in eta band in EMCal acceptance
        if(mom.Phi() > (phi-conesize) && mom.Phi() < (phi+conesize))
        {
          etaUEptsum+=mom.Pt();
          fhEtaBand->Fill(mom.Eta(),mom.Phi());
        }
        
        fhPtInCone->Fill(pt, mom.Pt());
        if(fFillPileUpHistograms && GetReader()->IsPileUpFromSPD())  fhPtInConePileUp->Fill(pt,mom.Pt());
        if (GetEventCentrality()) fhPtInConeCent->Fill(GetEventCentrality(),mom.Pt());
        coneptsum+=mom.Pt();
      }
    }
    
    //normalize phi/eta band per area unit
    fhPhiUEConeSumPt->Fill(pt, phiUEptsum*(TMath::Pi()*conesize*conesize)/(2*conesize*2*TMath::Pi()));
    fhEtaUEConeSumPt->Fill(pt, etaUEptsum*(TMath::Pi()*conesize*conesize)/(2*conesize*2*0.6));
    
    Double_t sumPhiUESub = coneptsum-(phiUEptsum*(TMath::Pi()*conesize*conesize)/(2*conesize*2*TMath::Pi()));
    Double_t sumEtaUESub = coneptsum-(etaUEptsum*(TMath::Pi()*conesize*conesize)/(2*conesize*2*0.6));
    
    fhConeSumPtPhiUESub->Fill(pt,sumPhiUESub);
    fhConeSumPtEtaUESub->Fill(pt,sumEtaUESub);
    
    
    if(GetDebug() > 1) printf("AliAnaParticleIsolation::MakeAnalysisFillHistograms() - Particle %d Energy Sum in Isolation Cone %2.2f\n", iaod, coneptsum);
    
    if(!fReMakeIC) fhConeSumPt->Fill(pt,coneptsum);
    
    Int_t mcTag = aod->GetTag() ;
    Int_t clID  = aod->GetCaloLabel(0) ;
    
    if(GetDebug() > 0) printf(" AliAnaParticleIsolation::MakeAnalysisFillHistograms() - pt %1.1f, eta %1.1f, phi %1.1f\n",pt, eta, phi);
    
    FillTrackMatchingShowerShapeControlHistograms(isolation, clID,aod->GetFiducialArea(),mcTag,reftracks,refclusters,aod,GetReader(), GetCaloPID());
    
    if(isolation)
    {    
      if(GetDebug() > 1) printf("AliAnaParticleIsolation::MakeAnalysisFillHistograms() - Particle %d ISOLATED: fill histograms\n", iaod);
      
      // Fill histograms to undertand pile-up before other cuts applied
      // Remember to relax time cuts in the reader
      FillPileUpHistograms(clID);     
      
      fhEIso      ->Fill(energy);
      fhPtIso     ->Fill(pt);
      fhPhiIso    ->Fill(pt,phi);
      fhEtaIso    ->Fill(pt,eta);
      fhEtaPhiIso ->Fill(eta,phi);
      
      if(decay) 
      {
        fhPtDecayIso->Fill(pt);
        fhEtaPhiDecayIso->Fill(eta,phi);
      }
      
      if(IsDataMC())
      {
        
        if(GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPhoton))
        {
          fhPtIsoMCPhoton  ->Fill(pt);
        }        
        
        if(GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPrompt))
        {
          fhPtIsoPrompt  ->Fill(pt);
          fhPhiIsoPrompt ->Fill(pt,phi);
          fhEtaIsoPrompt ->Fill(pt,eta);
        }
        else if(GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCFragmentation))
        {
          fhPtIsoFragmentation  ->Fill(pt);
          fhPhiIsoFragmentation ->Fill(pt,phi);
          fhEtaIsoFragmentation ->Fill(pt,eta);
        }
        else if(GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPi0))
        {
          fhPtIsoPi0  ->Fill(pt);
          fhPhiIsoPi0 ->Fill(pt,phi);
          fhEtaIsoPi0 ->Fill(pt,eta);
        }        
        else if(GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPi0Decay))
        {
          fhPtIsoPi0Decay  ->Fill(pt);
          fhPhiIsoPi0Decay ->Fill(pt,phi);
          fhEtaIsoPi0Decay ->Fill(pt,eta);
        }
        else if(GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCEtaDecay))
        {
          fhPtIsoEtaDecay  ->Fill(pt);
          fhPhiIsoEtaDecay ->Fill(pt,phi);
          fhEtaIsoEtaDecay ->Fill(pt,eta);
        }        
        else if(GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCOtherDecay))
        {
          fhPtIsoOtherDecay  ->Fill(pt);
          fhPhiIsoOtherDecay ->Fill(pt,phi);
          fhEtaIsoOtherDecay ->Fill(pt,eta);
        }
        //        else if(GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCConversion))
        //        {
        //          fhPtIsoConversion  ->Fill(pt);
        //          fhPhiIsoConversion ->Fill(pt,phi);
        //          fhEtaIsoConversion ->Fill(pt,eta);
        //        }
        else if(!GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCElectron))// anything else but electrons
        {
          fhPtIsoHadron  ->Fill(pt);
          fhPhiIsoHadron ->Fill(pt,phi);
          fhEtaIsoHadron ->Fill(pt,eta);
        }
      }//Histograms with MC
      
    }//Isolated histograms
    else // NON isolated
    {
      if(GetDebug() > 1) printf("AliAnaParticleIsolation::MakeAnalysisFillHistograms() - Particle %d NOT ISOLATED, fill histograms\n", iaod);
      
      fhPtNoIso    ->Fill(pt);
      fhEtaPhiNoIso->Fill(eta,phi);
      
      if(decay) 
      {
        fhPtDecayNoIso    ->Fill(pt);
        fhEtaPhiDecayNoIso->Fill(eta,phi);
      }
      
      if(IsDataMC())
      {
        if     (GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPhoton))        fhPtNoIsoMCPhoton     ->Fill(pt);
        if     (GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPi0))           fhPtNoIsoPi0          ->Fill(pt);
        else if(GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPi0Decay))      fhPtNoIsoPi0Decay     ->Fill(pt);
        else if(GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCEtaDecay))      fhPtNoIsoEtaDecay     ->Fill(pt);
        else if(GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCOtherDecay))    fhPtNoIsoOtherDecay   ->Fill(pt);
        else if(GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPrompt))        fhPtNoIsoPrompt       ->Fill(pt);
        else if(GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCFragmentation)) fhPtNoIsoFragmentation->Fill(pt);
        //        else if(GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCConversion))    fhPtNoIsoConversion   ->Fill(pt);
        else if(!GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCElectron))     fhPtNoIsoHadron      ->Fill(pt);        
      }
    }
  }// aod loop
  
}


//_____________________________________________________________________________________
void  AliAnaParticleIsolation::MakeSeveralICAnalysis(AliAODPWG4ParticleCorrelation* ph) 
{
  
  //Isolation Cut Analysis for both methods and different pt cuts and cones
  Float_t ptC   = ph->Pt();	
  Float_t etaC  = ph->Eta();
  Float_t phiC  = ph->Phi();
  Int_t   tag   = ph->GetTag(); 
  Bool_t  decay = ph->IsTagged();
  
  if(GetDebug() > 0) printf("AliAnaParticleIsolation::MakeSeveralICAnalysis() - Isolate pT %2.2f\n",ptC);
  
  //Keep original setting used when filling AODs, reset at end of analysis  
  Float_t ptthresorg = GetIsolationCut()->GetPtThreshold();
  Float_t ptfracorg  = GetIsolationCut()->GetPtFraction();
  Float_t rorg       = GetIsolationCut()->GetConeSize();
  
  Float_t coneptsum = 0 ; 
  Int_t   n    [10][10];//[fNCones][fNPtThresFrac];
  Int_t   nfrac[10][10];//[fNCones][fNPtThresFrac];
  Bool_t  isolated  = kFALSE;
  Int_t   nCone     = 0;
  Int_t   nFracCone = 0;
 
  // fill hist with all particles before isolation criteria
  fhPtNoIso    ->Fill(ptC);
  fhEtaPhiNoIso->Fill(etaC,phiC);
  
  if(IsDataMC())
  {
    if     (GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPhoton))        fhPtNoIsoMCPhoton     ->Fill(ptC);
    if     (GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0))           fhPtNoIsoPi0          ->Fill(ptC);
    else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0Decay))      fhPtNoIsoPi0Decay     ->Fill(ptC);
    else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCEtaDecay))      fhPtNoIsoEtaDecay     ->Fill(ptC);
    else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCOtherDecay))    fhPtNoIsoOtherDecay   ->Fill(ptC);
    else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPrompt))        fhPtNoIsoPrompt       ->Fill(ptC);
    else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCFragmentation)) fhPtNoIsoFragmentation->Fill(ptC);
//    else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCConversion))    fhPtNoIsoConversion   ->Fill(ptC);
    else if(!GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCElectron))    fhPtNoIsoHadron      ->Fill(ptC);
  }
  
  if(decay) 
  {
    fhPtDecayNoIso    ->Fill(ptC);
    fhEtaPhiDecayNoIso->Fill(etaC,phiC);
  }
  //Get vertex for photon momentum calculation
  Double_t vertex[] = {0,0,0} ; //vertex ;
  if(GetReader()->GetDataType() != AliCaloTrackReader::kMC) 
    {
      GetReader()->GetVertex(vertex);
    }

  //Loop on cone sizes
  for(Int_t icone = 0; icone<fNCones; icone++)
  {
    //Recover reference arrays with clusters and tracks
    TObjArray * refclusters = ph->GetObjArray(GetAODObjArrayName()+"Clusters");
    TObjArray * reftracks   = ph->GetObjArray(GetAODObjArrayName()+"Tracks");
    
    //If too small or too large pt, skip
    if(ptC < GetMinPt() || ptC > GetMaxPt() ) continue ; 
 
   //In case a more strict IC is needed in the produced AOD

    nCone=0; nFracCone = 0; isolated = kFALSE; coneptsum = 0; 
  
  GetIsolationCut()->SetSumPtThreshold(100);
  GetIsolationCut()->SetPtThreshold(100);
  GetIsolationCut()->SetPtFraction(100);
  GetIsolationCut()->SetConeSize(fConeSizes[icone]);
  GetIsolationCut()->MakeIsolationCut(reftracks,   refclusters, 
                                          GetReader(), GetCaloPID(),
                                          kFALSE, ph, "", 
				      nCone,nFracCone,coneptsum, isolated);

        
    fhSumPtLeadingPt[icone]->Fill(ptC,coneptsum);  
    
    // retreive pt tracks to fill histo vs. pt leading
    //Fill pt distribution of particles in cone
    //fhPtLeadingPt(),fhFRSumPtLeadingPt(),fhFRPtLeadingPt(),
    
    //Tracks
    coneptsum = 0;
    Double_t sumptFR = 0. ;
    TObjArray * trackList   = GetCTSTracks() ;
    for(Int_t itrack=0; itrack < trackList->GetEntriesFast(); itrack++)
    {
      AliVTrack* track = (AliVTrack *) trackList->At(itrack);
      //fill the histograms at forward range
      if(!track)
      {
        printf("AliAnaParticleIsolation::MakeAnalysisFillHistograms() - Track not available?");
        continue;
      }
      
      Double_t dPhi = phiC - track->Phi() + TMath::PiOver2();
      Double_t dEta = etaC - track->Eta();
      Double_t arg  = dPhi*dPhi + dEta*dEta;
      if(TMath::Sqrt(arg) < fConeSizes[icone])
      {
        fhFRPtLeadingPt[icone]->Fill(ptC,TMath::Sqrt(track->Px()*track->Px()+track->Py()*track->Py()));
        sumptFR+=track->Pt();
      }
      
      dPhi = phiC - track->Phi() - TMath::PiOver2();
      arg  = dPhi*dPhi + dEta*dEta;
      if(TMath::Sqrt(arg) < fConeSizes[icone])
      {
        fhFRPtLeadingPt[icone]->Fill(ptC,TMath::Sqrt(track->Px()*track->Px()+track->Py()*track->Py()));
        sumptFR+=track->Pt();
      }      
    }
    fhFRSumPtLeadingPt[icone]->Fill(ptC,sumptFR);
    if(reftracks)
    {  
      for(Int_t itrack=0; itrack < reftracks->GetEntriesFast(); itrack++)
      {
        AliVTrack* track = (AliVTrack *) reftracks->At(itrack);
        fhPtLeadingPt[icone]->Fill(ptC,TMath::Sqrt(track->Px()*track->Px()+track->Py()*track->Py()));
        coneptsum+=track->Pt();
      }
    }    
    //CaloClusters
    if(refclusters)
    {    
      TLorentzVector mom ;
      for(Int_t icalo=0; icalo < refclusters->GetEntriesFast(); icalo++)
      {
        AliVCluster* calo = (AliVCluster *) refclusters->At(icalo);
        calo->GetMomentum(mom,vertex) ;//Assume that come from vertex in straight line
        
        fhPtLeadingPt[icone]->Fill(ptC, mom.Pt());
        coneptsum+=mom.Pt();
      }
    }
    ///////////////////


    //Loop on ptthresholds
    for(Int_t ipt = 0; ipt<fNPtThresFrac ;ipt++)
    {
      n    [icone][ipt]=0;
      nfrac[icone][ipt]=0;
      GetIsolationCut()->SetPtThreshold(fPtThresholds[ipt]);
      GetIsolationCut()->SetPtFraction(fPtFractions[ipt]) ;
      GetIsolationCut()->SetSumPtThreshold(fSumPtThresholds[ipt]);
      
      GetIsolationCut()->MakeIsolationCut(reftracks, refclusters,
                                          GetReader(), GetCaloPID(),
                                          kFALSE, ph, "",
                                          n[icone][ipt],nfrac[icone][ipt],coneptsum, isolated);
      
      if(!isolated) continue;
      //Normal ptThreshold cut
      
      if(GetDebug() > 0) printf(" AliAnaParticleIsolation::MakeSeveralICAnalysis() - cone size %1.1f, ptThres  %1.1f, sumptThresh  %1.1f, n %d, nfrac %d, coneptsum %2.2f, isolated %d\n",
                                fConeSizes[icone],fPtThresholds[ipt],fSumPtThresholds[ipt],n[icone][ipt],nfrac[icone][ipt],coneptsum, isolated);
      if(GetDebug() > 0) printf(" AliAnaParticleIsolation::MakeSeveralICAnalysis() - pt %1.1f, eta %1.1f, phi %1.1f\n",ptC, etaC, phiC);
      
      if(n[icone][ipt] == 0) 
      {
        if(GetDebug()>0) printf(" AliAnaParticleIsolation::MakeSeveralICAnalysis() - filling pt threshold loop\n");
        fhPtThresIsolated[icone][ipt]->Fill(ptC);
        fhEtaPhiPtThresIso[icone][ipt]->Fill(etaC,phiC);
        
        if(decay)
        {
          fhPtPtThresDecayIso[icone][ipt]->Fill(ptC);
          //	  fhEtaPhiPtThresDecayIso[icone][ipt]->Fill(etaC,phiC);
        }
        
        if(IsDataMC())
        {
          if     ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPrompt))        fhPtThresIsolatedPrompt[icone][ipt]       ->Fill(ptC) ;
//          else if( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCConversion))    fhPtThresIsolatedConversion[icone][ipt]   ->Fill(ptC) ;
          else if( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCFragmentation)) fhPtThresIsolatedFragmentation[icone][ipt]->Fill(ptC) ;
          else if( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0))           fhPtThresIsolatedPi0[icone][ipt]           ->Fill(ptC) ;
          else if( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0Decay))      fhPtThresIsolatedPi0Decay[icone][ipt]     ->Fill(ptC) ;
          else if( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCEtaDecay))      fhPtThresIsolatedEtaDecay[icone][ipt]     ->Fill(ptC) ;
          else if( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCOtherDecay))    fhPtThresIsolatedOtherDecay[icone][ipt]   ->Fill(ptC) ;
          else if(!GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCElectron))      fhPtThresIsolatedHadron[icone][ipt]      ->Fill(ptC) ;
        }
      }
      
      // pt in cone fraction
      if(nfrac[icone][ipt] == 0)
      {
        if(GetDebug()>0) printf(" AliAnaParticleIsolation::MakeSeveralICAnalysis() - filling frac loop\n");
        fhPtFracIsolated[icone][ipt]->Fill(ptC);
        fhEtaPhiPtFracIso[icone][ipt]->Fill(etaC,phiC);
        
        if(decay)
        {
          fhPtPtFracDecayIso[icone][ipt]->Fill(ptC);
          fhEtaPhiPtFracDecayIso[icone][ipt]->Fill(etaC,phiC);
        }
        
        if(IsDataMC())
        {
          if     ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPrompt))        fhPtFracIsolatedPrompt[icone][ipt]       ->Fill(ptC) ;
//          else if( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCConversion))    fhPtFracIsolatedConversion[icone][ipt]   ->Fill(ptC) ;
          else if( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCFragmentation)) fhPtFracIsolatedFragmentation[icone][ipt]->Fill(ptC) ;
          else if( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0))           fhPtFracIsolatedPi0[icone][ipt]          ->Fill(ptC) ;
          else if( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0Decay))      fhPtFracIsolatedPi0Decay[icone][ipt]     ->Fill(ptC) ;
          else if( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCEtaDecay))      fhPtFracIsolatedEtaDecay[icone][ipt]     ->Fill(ptC) ;
          else if( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCOtherDecay))    fhPtFracIsolatedOtherDecay[icone][ipt]   ->Fill(ptC) ;
          else if(!GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCElectron))      fhPtFracIsolatedHadron[icone][ipt]->Fill(ptC) ;
        }
      }
      
      if(GetDebug()>0) printf(" AliAnaParticleIsolation::MakeSeveralICAnalysis() - checking IC method : %i\n",GetIsolationCut()->GetICMethod());
      
      //Pt threshold on pt cand/ sum in cone histograms
      if(coneptsum<fSumPtThresholds[ipt])
      {//      if((GetIsolationCut()->GetICMethod())==1){//kSumPtIC){
        if(GetDebug()>0) printf(" AliAnaParticleIsolation::MakeSeveralICAnalysis() - filling sum loop\n");
        fhPtSumIsolated[icone][ipt]->Fill(ptC) ;
        fhEtaPhiPtSumIso[icone][ipt]->Fill(etaC, phiC) ;
        if(decay)
        {
          fhPtPtSumDecayIso[icone][ipt]->Fill(ptC);
          fhEtaPhiPtSumDecayIso[icone][ipt]->Fill(etaC, phiC) ;
        }
      }
      
    // pt sum pt frac method
//    if( ((fPtFractions[ipt]*ptC < fSumPtThresholds[ipt]) && (coneptsum < fSumPtThresholds[ipt])) || ((fPtFractions[ipt]*ptC > fSumPtThresholds[ipt]) && (coneptsum < fPtFractions[ipt]*ptC)) )

      if(coneptsum < fPtFractions[ipt]*ptC)
       {
        if(GetDebug()>0) printf(" AliAnaParticleIsolation::MakeSeveralICAnalysis() - filling PtFrac PtSum loop\n");
        fhPtFracPtSumIso[icone][ipt]->Fill(ptC) ;
        fhEtaPhiFracPtSumIso[icone][ipt]->Fill(etaC,phiC) ;
        
        if(decay)
        {
          fhPtFracPtSumDecayIso[icone][ipt]->Fill(ptC);
          fhEtaPhiFracPtSumDecayIso[icone][ipt]->Fill(etaC,phiC);
        }
      }
      
      // density method
      Float_t cellDensity = GetIsolationCut()->GetCellDensity( ph, GetReader());
      if(coneptsum<fSumPtThresholds[ipt]*cellDensity)
      {//(GetIsolationCut()->GetICMethod())==4){//kSumDensityIC) {
        if(GetDebug()>0) printf(" AliAnaParticleIsolation::MakeSeveralICAnalysis() - filling density loop\n");
        fhPtSumDensityIso[icone][ipt]->Fill(ptC) ;
        fhEtaPhiSumDensityIso[icone][ipt]->Fill(etaC,phiC) ;
        
        if(decay)
        {
          fhPtSumDensityDecayIso[icone][ipt]->Fill(ptC);
          fhEtaPhiSumDensityDecayIso[icone][ipt]->Fill(etaC, phiC);
        }
        
      }
    }//pt thresh loop
    
    if(IsDataMC())
    {
      if     ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPrompt))        fhPtSumIsolatedPrompt[icone]       ->Fill(ptC,coneptsum) ;
//      else if( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCConversion))    fhPtSumIsolatedConversion[icone]   ->Fill(ptC,coneptsum) ;
      else if( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCFragmentation)) fhPtSumIsolatedFragmentation[icone]->Fill(ptC,coneptsum) ;
      else if( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0))           fhPtSumIsolatedPi0[icone]          ->Fill(ptC,coneptsum) ;
      else if( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0Decay))      fhPtSumIsolatedPi0Decay[icone]     ->Fill(ptC,coneptsum) ;
      else if( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCEtaDecay))      fhPtSumIsolatedEtaDecay[icone]     ->Fill(ptC,coneptsum) ;
      else if( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCOtherDecay))    fhPtSumIsolatedOtherDecay[icone]   ->Fill(ptC,coneptsum) ;
      else if(!GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCElectron))      fhPtSumIsolatedHadron[icone]->Fill(ptC,coneptsum) ;
    }
    
  }//cone size loop
  
  //Reset original parameters for AOD analysis
  GetIsolationCut()->SetPtThreshold(ptthresorg);
  GetIsolationCut()->SetPtFraction(ptfracorg);
  GetIsolationCut()->SetConeSize(rorg);
  
}

//_____________________________________________________________
void AliAnaParticleIsolation::Print(const Option_t * opt) const
{
  
  //Print some relevant parameters set for the analysis
  if(! opt)
    return;
  
  printf("**** Print %s %s ****\n", GetName(), GetTitle() ) ;
  AliAnaCaloTrackCorrBaseClass::Print(" ");
  
  printf("ReMake Isolation          = %d \n",  fReMakeIC) ;
  printf("Make Several Isolation    = %d \n",  fMakeSeveralIC) ;
  printf("Calorimeter for isolation = %s \n",  fCalorimeter.Data()) ;
  
  if(fMakeSeveralIC)
  {
    printf("N Cone Sizes       =     %d\n", fNCones) ; 
    printf("Cone Sizes          =    \n") ;
    for(Int_t i = 0; i < fNCones; i++)
      printf("  %1.2f;",  fConeSizes[i]) ;
    printf("    \n") ;
    
    printf("N pT thresholds/fractions = %d\n", fNPtThresFrac) ;
    printf(" pT thresholds         =    \n") ;
    for(Int_t i = 0; i < fNPtThresFrac; i++)
      printf("   %2.2f;",  fPtThresholds[i]) ;
    
    printf("    \n") ;
    
    printf(" pT fractions         =    \n") ;
    for(Int_t i = 0; i < fNPtThresFrac; i++)
      printf("   %2.2f;",  fPtFractions[i]) ;
    
    printf("    \n") ;
    
    printf("sum pT thresholds         =    \n") ;
    for(Int_t i = 0; i < fNPtThresFrac; i++)
      printf("   %2.2f;",  fSumPtThresholds[i]) ;
    
    
  }  
  
  printf("Histograms: %3.1f < pT sum < %3.1f,  Nbin = %d\n",    fHistoPtSumMin,    fHistoPtSumMax,    fHistoNPtSumBins   );
  printf("Histograms: %3.1f < pT in cone < %3.1f, Nbin = %d\n", fHistoPtInConeMin, fHistoPtInConeMax, fHistoNPtInConeBins);
  
  printf("    \n") ;
  
} 

