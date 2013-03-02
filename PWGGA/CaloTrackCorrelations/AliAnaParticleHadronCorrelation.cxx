/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//_________________________________________________________________________
// Class for the analysis of particle - hadron correlations
// Particle (for example direct gamma) must be found in a previous analysis 
//-- Author: Gustavo Conesa (LNF-INFN) 

//  Modified by Yaxian Mao:
// 1. add the UE subtraction for corrlation study
// 2. change the correlation variable
// 3. Only use leading particle(cluster/track) as trigger for correlation (2010/07/02)
// 4. Make decay photon-hadron correlations where decay contribute pi0 mass (2010/09/09)
// 5. fill the pout to extract kt at the end, also to study charge asymmetry(2010/10/06) 
// 6. Add the possibility for event selection analysis based on vertex and multiplicity bins (10/10/2010)
// 7. change the way of delta phi cut for UE study due to memory issue (reduce histograms)
// 8. Add the possibility to request the absolute leading particle at the near side or not, set trigger bins, general clean-up (08/2011)
//////////////////////////////////////////////////////////////////////////////


// --- ROOT system ---
//#include "TClonesArray.h"
#include <TClass.h>
#include <TMath.h>
#include <TH2F.h>
#include <TDatabasePDG.h>

//---- ANALYSIS system ----
#include "AliNeutralMesonSelection.h" 
#include "AliAnaParticleHadronCorrelation.h" 
#include "AliCaloTrackReader.h"
#include "AliAODPWG4ParticleCorrelation.h"
#include "AliFiducialCut.h"
#include "AliVTrack.h"
#include "AliVCluster.h"
#include "AliMCAnalysisUtils.h"
#include "TParticle.h"
#include "AliStack.h"
#include "AliAODMCParticle.h"
#include "AliMixedEvent.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliEventplane.h"

ClassImp(AliAnaParticleHadronCorrelation)


//___________________________________________________________________
  AliAnaParticleHadronCorrelation::AliAnaParticleHadronCorrelation(): 
    AliAnaCaloTrackCorrBaseClass(),
    fMinTriggerPt(0.),   
    fMaxAssocPt(1000.),             fMinAssocPt(0.),   
    fDeltaPhiMaxCut(0.),            fDeltaPhiMinCut(0.),   
    fSelectIsolated(0),             fMakeSeveralUE(0),              
    fUeDeltaPhiMaxCut(0.),          fUeDeltaPhiMinCut(0.), 
    fPi0AODBranchName(""),          fNeutralCorr(0),       
    fPi0Trigger(0),                 fDecayTrigger(0),
    fMakeAbsoluteLeading(0),        fMakeNearSideLeading(0),      
    fLeadingTriggerIndex(-1),       fHMPIDCorrelation(0),  fFillBradHisto(0),
    fNAssocPtBins(0),               fAssocPtBinLimit(),
    fCorrelVzBin(0),
    fListMixTrackEvents(),          fListMixCaloEvents(),
    fUseMixStoredInReader(0),       fFillNeutralEventMixPool(0),
    fM02MaxCut(0),                  fM02MinCut(0),  
    fFillPileUpHistograms(0),       
    //Histograms
    fhPtInput(0),                   fhPtFidCut(0),
    fhPtLeading(0),                 fhPtLeadingVtxBC0(0),
    fhPtLeadingVzBin(0),            fhPtLeadingBin(0),                 
    fhPhiLeading(0),                fhEtaLeading(0),   
    fhPtLeadingMC(),
    fhPtLeadingCentrality(0),       fhPtLeadingEventPlane(0), 
    fhLeadingEventPlaneCentrality(0),
    fhPtLeadingMixed(0),            fhPtLeadingMixedVzBin(0), fhPtLeadingMixedBin(0),              
    fhPhiLeadingMixed(0),           fhEtaLeadingMixed(0), 
    fhDeltaPhiDeltaEtaCharged(0),
    fhPhiCharged(0),                fhEtaCharged(0), 
    fhDeltaPhiCharged(0),           fhDeltaEtaCharged(0), 
    fhDeltaPhiChargedPt(0),         fhDeltaPhiUeChargedPt(0), 
    fhUePart(0),
    fhXECharged(0),                 fhXEUeCharged(0),
    fhXEPosCharged(0),              fhXENegCharged(0),
    fhPtHbpXECharged(0),            fhPtHbpXEUeCharged(0),
    fhZTCharged(0),                 fhZTUeCharged(0),
    fhZTPosCharged(0),              fhZTNegCharged(0),
    fhPtHbpZTCharged(0),            fhPtHbpZTUeCharged(0),
    fhXEChargedMC(),                fhDeltaPhiChargedMC(),  
    fhDeltaPhiDeltaEtaChargedPtA3GeV(0),
    fhDeltaPhiChargedPtA3GeV(0),    fhDeltaEtaChargedPtA3GeV(0),
    //Pile-Up
    fhDeltaPhiChargedPileUp(),      fhDeltaEtaChargedPileUp(),
    fhDeltaPhiChargedPtA3GeVPileUp(), fhDeltaEtaChargedPtA3GeVPileUp(),
    fhXEChargedPileUp(),            fhXEUeChargedPileUp(),
    fhZTChargedPileUp(),            fhZTUeChargedPileUp(), 
    fhPtTrigChargedPileUp(),
    fhDeltaPhiChargedOtherBC(),     fhDeltaPhiChargedPtA3GeVOtherBC(),
    fhXEChargedOtherBC(),           fhXEUeChargedOtherBC(),
    fhZTChargedOtherBC(),           fhZTUeChargedOtherBC(),
    fhPtTrigChargedOtherBC(),
    fhDeltaPhiChargedBC0(),         fhDeltaPhiChargedPtA3GeVBC0(),
    fhXEChargedBC0(),               fhXEUeChargedBC0(),
    fhZTChargedBC0(),               fhZTUeChargedBC0(),
    fhPtTrigChargedBC0(),
    fhDeltaPhiChargedVtxBC0(),      fhDeltaPhiChargedPtA3GeVVtxBC0(),
    fhXEChargedVtxBC0(),            fhXEUeChargedVtxBC0(),
    fhZTChargedVtxBC0(),            fhZTUeChargedVtxBC0(),
    fhPtTrigChargedVtxBC0(),
    fhDeltaPhiUeLeftCharged(0),     fhDeltaPhiUeRightCharged(0),
    fhDeltaPhiUeLeftUpCharged(0),   fhDeltaPhiUeRightUpCharged(0),
    fhDeltaPhiUeLeftDownCharged(0), fhDeltaPhiUeRightDownCharged(0),
    fhXEUeLeftCharged(0),           fhXEUeRightCharged(0),
    fhXEUeLeftUpCharged(0),         fhXEUeRightUpCharged(0),
    fhXEUeLeftDownCharged(0),       fhXEUeRightDownCharged(0),
    fhPtHbpXEUeLeftCharged(0),      fhPtHbpXEUeRightCharged(0), 
    fhZTUeLeftCharged(0),           fhZTUeRightCharged(0),
    fhPtHbpZTUeLeftCharged(0),      fhPtHbpZTUeRightCharged(0), 
    fhPtTrigPout(0),                fhPtTrigCharged(0),
    fhTrigDeltaPhiCharged(0x0),     fhTrigDeltaEtaCharged(0x0),
    fhTrigXECorr(0x0),              fhTrigXEUeCorr(0x0),
    fhTrigZTCorr(0x0),              fhTrigZTUeCorr(0x0),
    fhAssocPtBkg(0),                fhDeltaPhiDeltaEtaAssocPtBin(0),
    fhDeltaPhiAssocPtBin(0),        
    fhDeltaPhiAssocPtBinDEta08(0),  fhDeltaPhiAssocPtBinDEta0(0),
    fhDeltaPhiAssocPtBinHMPID(0),   fhDeltaPhiAssocPtBinHMPIDAcc(0),
    fhDeltaPhiBradAssocPtBin(0),    fhDeltaPhiBrad(0),
    fhXEAssocPtBin(0),              fhZTAssocPtBin(0),             
    fhDeltaPhiDeltaEtaNeutral(0), 
    fhPhiNeutral(0),                fhEtaNeutral(0), 
    fhDeltaPhiNeutral(0),           fhDeltaEtaNeutral(0),
    fhDeltaPhiNeutralPt(0),         fhDeltaPhiUeNeutralPt(0), 
    fhXENeutral(0),                 fhXEUeNeutral(0),
    fhPtHbpXENeutral(0),            fhPtHbpXEUeNeutral(0),
    fhZTNeutral(0),                 fhZTUeNeutral(0),
    fhPtHbpZTNeutral(0),            fhPtHbpZTUeNeutral(0),
    fhDeltaPhiUeLeftNeutral(0),     fhDeltaPhiUeRightNeutral(0),
    fhXEUeLeftNeutral(0),           fhXEUeRightNeutral(0),
    fhPtHbpXEUeLeftNeutral(0),      fhPtHbpXEUeRightNeutral(0),
    fhZTUeLeftNeutral(0),           fhZTUeRightNeutral(0),
    fhPtHbpZTUeLeftNeutral(0),      fhPtHbpZTUeRightNeutral(0),
    fhPtPi0DecayRatio(0),
    fhDeltaPhiDecayCharged(0),      fhXEDecayCharged(0), fhZTDecayCharged(0), 
    fhDeltaPhiDecayNeutral(0),      fhXEDecayNeutral(0), fhZTDecayNeutral(0),
    fhDeltaPhiDecayChargedAssocPtBin(0), 
    fhXEDecayChargedAssocPtBin(0),  fhZTDecayChargedAssocPtBin(0),                
    fh2phiLeadingParticle(0x0),     fhMCPtLeading(0),
    fhMCPhiLeading(0),              fhMCEtaLeading(0),
    fhMCEtaCharged(0),              fhMCPhiCharged(0), 
    fhMCDeltaEtaCharged(0),         fhMCDeltaPhiCharged(0x0),
    fhMCDeltaPhiDeltaEtaCharged(0), fhMCDeltaPhiChargedPt(0),
    fhMCPtXECharged(0),             fhMCPtXEUeCharged(0),
    fhMCPtHbpXECharged(0),          fhMCPtHbpXEUeCharged(0),
    fhMCUePart(0),
    fhMCPtZTCharged(0),             fhMCPtHbpZTCharged(0),
    fhMCPtTrigPout(0),              fhMCPtAssocDeltaPhi(0),
    //Mixing
    fhNEventsTrigger(0),            fhNtracksMB(0),                 fhNclustersMB(0),
    fhMixDeltaPhiCharged(0),        fhMixDeltaPhiDeltaEtaCharged(0),
    fhMixXECharged(0),              fhMixXEUeCharged(0),            fhMixHbpXECharged(0),
    fhMixDeltaPhiChargedAssocPtBin(),
    fhMixDeltaPhiChargedAssocPtBinDEta08(),
    fhMixDeltaPhiChargedAssocPtBinDEta0(),
    fhMixDeltaPhiDeltaEtaChargedAssocPtBin(),
    fhEventBin(0),                  fhEventMixBin(0)
{ 
  //Default Ctor 
    
  //Initialize parameters
  InitParameters();
  
  for(Int_t i = 0; i < 7; i++)
  { 
    fhPtLeadingMC[i] = 0;
    fhXEChargedMC[i] = 0;
    fhDeltaPhiChargedMC[i] = 0;
  }
  
  for(Int_t i = 0; i < 7; i++)
  {
    fhPtLeadingPileUp             [i] = 0 ;
    fhDeltaPhiChargedPileUp       [i] = 0 ; fhDeltaEtaChargedPileUp       [i] = 0 ;
    fhXEChargedPileUp             [i] = 0 ; fhXEUeChargedPileUp           [i] = 0 ;
    fhZTChargedPileUp             [i] = 0 ; fhZTUeChargedPileUp           [i] = 0 ;
    fhPtTrigChargedPileUp         [i] = 0 ;
    fhDeltaPhiChargedPtA3GeVPileUp[i] = 0 ; fhDeltaEtaChargedPtA3GeVPileUp[i] = 0 ;
  }
  
}

//_________________________________________________________________
AliAnaParticleHadronCorrelation::~AliAnaParticleHadronCorrelation()
{
  // Remove event containers
  
  if(DoOwnMix())
  {
    if(fListMixTrackEvents)
    {      
      for(Int_t iz=0; iz < GetNZvertBin(); iz++)
      {
        for(Int_t ic=0; ic < GetNCentrBin(); ic++)
        {
          for(Int_t irp=0; irp<GetNRPBin(); irp++)
          {
            Int_t bin = GetEventMixBin(ic, iz, irp);
            fListMixTrackEvents[bin]->Delete() ;
            delete fListMixTrackEvents[bin] ;
          }
        }
      }
    }

    delete[] fListMixTrackEvents;

    if(fListMixCaloEvents)
    {      
      for(Int_t iz=0; iz < GetNZvertBin(); iz++)
      {
        for(Int_t ic=0; ic < GetNCentrBin(); ic++)
        {
          for(Int_t irp=0; irp<GetNRPBin(); irp++)
          {
            Int_t bin = GetEventMixBin(ic, iz, irp);
            fListMixCaloEvents[bin]->Delete() ;
            delete fListMixCaloEvents[bin] ;
          }
        }
      }
    }
  
    delete[] fListMixCaloEvents;
    
  }
}

//______________________________________________________________________________________________________________________________________________________
void AliAnaParticleHadronCorrelation::FillChargedAngularCorrelationHistograms(Float_t ptAssoc,  Float_t ptTrig,      Int_t   bin,
                                                                              Float_t phiAssoc, Float_t phiTrig,     Float_t &     deltaPhi,
                                                                              Float_t etaAssoc, Float_t etaTrig,  
                                                                              Bool_t  decay,    Float_t hmpidSignal, Int_t  outTOF,
                                                                              Int_t nTracks,       Int_t   mcTag)
{
  // Fill angular correlation related histograms
  
  Float_t deltaEta    = etaTrig-etaAssoc;
          deltaPhi    = phiTrig-phiAssoc;
  Float_t deltaPhiOrg = deltaPhi;
  
  if(deltaPhi <= -TMath::PiOver2()) deltaPhi+=TMath::TwoPi();
  if(deltaPhi > 3*TMath::PiOver2()) deltaPhi-=TMath::TwoPi();
  
  fhEtaCharged       ->Fill(ptAssoc,etaAssoc);
  fhPhiCharged       ->Fill(ptAssoc,phiAssoc);
  fhDeltaEtaCharged  ->Fill(ptTrig ,deltaEta);
  fhDeltaPhiCharged  ->Fill(ptTrig ,deltaPhi);
  fhDeltaPhiChargedPt->Fill(ptAssoc, deltaPhi);
  fhDeltaPhiDeltaEtaCharged->Fill(deltaPhi, deltaEta);
  
  if(ptAssoc > 3 )
  {
    fhDeltaEtaChargedPtA3GeV        ->Fill(ptTrig  ,deltaEta);
    fhDeltaPhiChargedPtA3GeV        ->Fill(ptTrig  ,deltaPhi);
    fhDeltaPhiDeltaEtaChargedPtA3GeV->Fill(deltaPhi, deltaEta);    
  }   
  
  // Pile up studies
  
  if(fFillPileUpHistograms)
  {
    if     (outTOF==1)
    {
      fhDeltaPhiChargedOtherBC->Fill(ptTrig ,deltaPhi) ;
      if(ptAssoc > 3 ) fhDeltaPhiChargedPtA3GeVOtherBC->Fill(ptTrig ,deltaPhi) ;
    }
    else if(outTOF==0)
    {
        fhDeltaPhiChargedBC0->Fill(ptTrig ,deltaPhi) ;
        if(ptAssoc > 3 ) fhDeltaPhiChargedPtA3GeVBC0->Fill(ptTrig ,deltaPhi) ;
    }
    
    Int_t vtxBC = GetReader()->GetVertexBC();
    if(vtxBC == 0 || vtxBC==AliVTrack::kTOFBCNA)
    {
      fhDeltaPhiChargedVtxBC0->Fill(ptTrig ,deltaPhi) ;
      if(ptAssoc > 3 ) fhDeltaPhiChargedPtA3GeVVtxBC0->Fill(ptTrig ,deltaPhi) ;
    }
    
    if(GetReader()->IsPileUpFromSPD())               { fhDeltaEtaChargedPileUp[0]->Fill(ptTrig ,deltaEta) ; fhDeltaPhiChargedPileUp[0]->Fill(ptTrig ,deltaPhi) ; }
    if(GetReader()->IsPileUpFromEMCal())             { fhDeltaEtaChargedPileUp[1]->Fill(ptTrig ,deltaEta) ; fhDeltaPhiChargedPileUp[1]->Fill(ptTrig ,deltaPhi) ; }
    if(GetReader()->IsPileUpFromSPDOrEMCal())        { fhDeltaEtaChargedPileUp[2]->Fill(ptTrig ,deltaEta) ; fhDeltaPhiChargedPileUp[2]->Fill(ptTrig ,deltaPhi) ; }
    if(GetReader()->IsPileUpFromSPDAndEMCal())       { fhDeltaEtaChargedPileUp[3]->Fill(ptTrig ,deltaEta) ; fhDeltaPhiChargedPileUp[3]->Fill(ptTrig ,deltaPhi) ; }
    if(GetReader()->IsPileUpFromSPDAndNotEMCal())    { fhDeltaEtaChargedPileUp[4]->Fill(ptTrig ,deltaEta) ; fhDeltaPhiChargedPileUp[4]->Fill(ptTrig ,deltaPhi) ; }
    if(GetReader()->IsPileUpFromEMCalAndNotSPD())    { fhDeltaEtaChargedPileUp[5]->Fill(ptTrig ,deltaEta) ; fhDeltaPhiChargedPileUp[5]->Fill(ptTrig ,deltaPhi) ; }
    if(GetReader()->IsPileUpFromNotSPDAndNotEMCal()) { fhDeltaEtaChargedPileUp[6]->Fill(ptTrig ,deltaEta) ; fhDeltaPhiChargedPileUp[6]->Fill(ptTrig ,deltaPhi) ; }
    
    if(ptAssoc > 3 )
    {
      if(GetReader()->IsPileUpFromSPD())               { fhDeltaEtaChargedPtA3GeVPileUp[0]->Fill(ptTrig ,deltaEta) ; fhDeltaPhiChargedPtA3GeVPileUp[0]->Fill(ptTrig ,deltaPhi) ; }
      if(GetReader()->IsPileUpFromEMCal())             { fhDeltaEtaChargedPtA3GeVPileUp[1]->Fill(ptTrig ,deltaEta) ; fhDeltaPhiChargedPtA3GeVPileUp[1]->Fill(ptTrig ,deltaPhi) ; }
      if(GetReader()->IsPileUpFromSPDOrEMCal())        { fhDeltaEtaChargedPtA3GeVPileUp[2]->Fill(ptTrig ,deltaEta) ; fhDeltaPhiChargedPtA3GeVPileUp[2]->Fill(ptTrig ,deltaPhi) ; }
      if(GetReader()->IsPileUpFromSPDAndEMCal())       { fhDeltaEtaChargedPtA3GeVPileUp[3]->Fill(ptTrig ,deltaEta) ; fhDeltaPhiChargedPtA3GeVPileUp[3]->Fill(ptTrig ,deltaPhi) ; }
      if(GetReader()->IsPileUpFromSPDAndNotEMCal())    { fhDeltaEtaChargedPtA3GeVPileUp[4]->Fill(ptTrig ,deltaEta) ; fhDeltaPhiChargedPtA3GeVPileUp[4]->Fill(ptTrig ,deltaPhi) ; }
      if(GetReader()->IsPileUpFromEMCalAndNotSPD())    { fhDeltaEtaChargedPtA3GeVPileUp[5]->Fill(ptTrig ,deltaEta) ; fhDeltaPhiChargedPtA3GeVPileUp[5]->Fill(ptTrig ,deltaPhi) ; }
      if(GetReader()->IsPileUpFromNotSPDAndNotEMCal()) { fhDeltaEtaChargedPtA3GeVPileUp[6]->Fill(ptTrig ,deltaEta) ; fhDeltaPhiChargedPtA3GeVPileUp[6]->Fill(ptTrig ,deltaPhi) ; }
    }
  }
  
  if(IsDataMC())
  {
    Int_t mcIndex = GetMCTagHistogramIndex(mcTag);
    fhDeltaPhiChargedMC[mcIndex]->Fill(ptTrig , deltaPhi);
  }  
  
  if(fDecayTrigger && decay) fhDeltaPhiDecayCharged   ->Fill(ptTrig  , deltaPhi);      
  
  Double_t  dphiBrad = -100;
  if(fFillBradHisto)
  {
    dphiBrad = atan2(sin(deltaPhiOrg), cos(deltaPhiOrg))/TMath::Pi();//-1 to 1
    if(TMath::Abs(dphiBrad)>0.325 && TMath::Abs(dphiBrad)<0.475)  //Hardcoded values, BAD, FIXME
    {
      fhAssocPtBkg->Fill(ptTrig, ptAssoc);
    }
    
    if(dphiBrad<-1./3) dphiBrad += 2;
    fhDeltaPhiBrad->Fill(ptTrig, dphiBrad);
  }
  
  // Fill histograms in bins of associated particle pT
  if(bin>=0)
  {
      fhDeltaPhiDeltaEtaAssocPtBin    [bin]->Fill(deltaPhi,deltaEta);

      fhDeltaPhiAssocPtBin            [bin]->Fill(ptTrig, deltaPhi);
    
    if(TMath::Abs(deltaEta)> 0.8) 
      fhDeltaPhiAssocPtBinDEta08      [bin]->Fill(ptTrig, deltaPhi);

    if(TMath::Abs(deltaEta)< 0.01) 
      fhDeltaPhiAssocPtBinDEta0       [bin]->Fill(ptTrig, deltaPhi);
    
    if (fFillBradHisto)
      fhDeltaPhiBradAssocPtBin        [bin]->Fill(ptTrig, dphiBrad);
    
    if(fDecayTrigger && decay)
      fhDeltaPhiDecayChargedAssocPtBin[bin]->Fill(ptTrig, deltaPhi);      
    
    if(fHMPIDCorrelation)
    {
      if( hmpidSignal > 0 )
      {
        //printf("Track pt %f with HMPID signal %f \n",pt,hmpidSignal);
        fhDeltaPhiAssocPtBinHMPID[bin]->Fill(ptTrig, deltaPhi);        
      }
      
      if(phiAssoc > 5*TMath::DegToRad() && phiAssoc < 20*TMath::DegToRad())
      {
        //printf("Track pt %f in HMPID acceptance phi %f \n ",pt,phi*TMath::RadToDeg() );
        fhDeltaPhiAssocPtBinHMPIDAcc[bin]->Fill(ptTrig, deltaPhi);        
      }
    }
  }
  
  //fill different multiplicity histogram
  if(DoEventSelect())
  {
    for(Int_t im = 0; im<GetMultiBin(); im++)
    {
      if(nTracks < ( GetMaxMulti() - GetMinMulti() )/GetMultiBin()*(im+1))
      {
        fhTrigDeltaPhiCharged[im]->Fill(ptTrig,deltaPhi);
        fhTrigDeltaEtaCharged[im]->Fill(ptTrig,deltaEta);
      }
    }
  }  
}

//____________________________________________________________________________________________________________________________________________________
Bool_t AliAnaParticleHadronCorrelation::FillChargedMCCorrelationHistograms(Float_t mcAssocPt,      Float_t mcAssocPhi, Float_t mcAssocEta,
                                                                           Float_t mcTrigPt, Float_t mcTrigPhi,  Float_t mcTrigEta)
{
  // Fill MC histograms independently of AOD or ESD
  
  //Select only hadrons in pt range
  if(mcAssocPt < fMinAssocPt || mcAssocPt > fMaxAssocPt) return kTRUE ; // exclude but continue
  
  if(mcAssocPhi < 0) mcAssocPhi+=TMath::TwoPi();    
  
  //remove trigger itself for correlation when use charged triggers 
  if(TMath::Abs(mcAssocPt -mcTrigPt ) < 1e-6 && 
     TMath::Abs(mcAssocPhi-mcTrigPhi) < 1e-6 && 
     TMath::Abs(mcAssocEta-mcTrigEta) < 1e-6)            return kTRUE ; // exclude but continue       
  
  // Absolute leading?
  if( fMakeAbsoluteLeading && mcAssocPt > mcTrigPt )     return kFALSE; // jump event
  
  //jump out this event if near side associated partile pt larger than trigger
  if( fMakeNearSideLeading && mcAssocPt > mcTrigPt && 
     TMath::Abs(mcAssocPhi-mcTrigPhi)<TMath::PiOver2() ) return kFALSE; // jump event
  
  Float_t mcdeltaPhi= mcTrigPhi-mcAssocPhi; 
  if(mcdeltaPhi <= -TMath::PiOver2()) mcdeltaPhi+=TMath::TwoPi();
  if(mcdeltaPhi > 3*TMath::PiOver2()) mcdeltaPhi-=TMath::TwoPi();            
  
  Float_t mcxE    =-mcAssocPt/mcTrigPt*TMath::Cos(mcdeltaPhi);// -(mcAssocPx*pxprim+mcAssocPy*pyprim)/(mcTrigPt*mcTrigPt);  
  Float_t mchbpXE =-100 ;
  if(mcxE > 0 ) mchbpXE = TMath::Log(1./mcxE);
  
  Float_t mczT    = mcAssocPt/mcTrigPt ;
  Float_t mchbpZT =-100 ;
  if(mczT > 0 ) mchbpZT = TMath::Log(1./mczT);
  
  //Selection within angular range
  if( mcdeltaPhi< -TMath::PiOver2())  mcdeltaPhi+=TMath::TwoPi();
  if( mcdeltaPhi>3*TMath::PiOver2())  mcdeltaPhi-=TMath::TwoPi();              
  
  Double_t mcpout = mcAssocPt*TMath::Sin(mcdeltaPhi) ; 
  
  if(GetDebug() > 0 )  
    printf("AliAnaParticleHadronCorrelation::MakeMCChargedCorrelation() - Charged hadron: track Pt %f, track Phi %f, phi trigger %f. Cuts:  delta phi  %2.2f < %2.2f < %2.2f, pT min %2.2f \n",
           mcAssocPt,mcAssocPhi, mcTrigPhi,fDeltaPhiMinCut, mcdeltaPhi, fDeltaPhiMaxCut, GetMinPt());              
  
  // Fill Histograms
  fhMCEtaCharged     ->Fill(mcAssocPt, mcAssocEta);
  fhMCPhiCharged     ->Fill(mcAssocPt, mcAssocPhi);
  fhMCDeltaEtaCharged->Fill(mcTrigPt,    mcTrigEta-mcAssocEta);
  fhMCDeltaPhiCharged->Fill(mcTrigPt,    mcdeltaPhi);
  fhMCPtAssocDeltaPhi->Fill(mcAssocPt, mcdeltaPhi);
  
  fhMCDeltaPhiDeltaEtaCharged->Fill(mcdeltaPhi,mcTrigEta-mcAssocEta);
  
  //delta phi cut for correlation
  if( (mcdeltaPhi > fDeltaPhiMinCut) && (mcdeltaPhi < fDeltaPhiMaxCut) ) 
  {
    fhMCDeltaPhiChargedPt->Fill(mcAssocPt,mcdeltaPhi);
    fhMCPtXECharged      ->Fill(mcTrigPt,mcxE); 
    fhMCPtHbpXECharged   ->Fill(mcTrigPt,mchbpXE);
    fhMCPtZTCharged      ->Fill(mcTrigPt,mczT); 
    fhMCPtHbpZTCharged   ->Fill(mcTrigPt,mchbpZT);
    fhMCPtTrigPout       ->Fill(mcTrigPt, mcpout) ;
  }

  //underlying event
  if ( (mcdeltaPhi > fUeDeltaPhiMinCut) && (mcdeltaPhi < fUeDeltaPhiMaxCut) )   
    {
       Double_t randomphi = gRandom->Uniform(TMath::Pi()/2,3*TMath::Pi()/2);
       Double_t mcUexE = -(mcAssocPt/mcTrigPt)*TMath::Cos(randomphi);
  
       if(mcUexE < 0.) mcUexE = -mcUexE;
    
       fhMCPtXEUeCharged->Fill(mcTrigPt,mcUexE);
       if(mcUexE > 0) fhMCPtHbpXEUeCharged->Fill(mcTrigPt,TMath::Log(1/mcUexE));

       fhMCUePart->Fill(mcTrigPt);
    }
  
  return kTRUE;
} 

//___________________________________________________________________________________________________________________
void AliAnaParticleHadronCorrelation::FillChargedMomentumImbalanceHistograms(Float_t ptTrig,   Float_t ptAssoc, 
                                                                             Float_t xE,       Float_t hbpXE, 
                                                                             Float_t zT,       Float_t hbpZT, 
                                                                             Float_t pout,
                                                                             Int_t   nTracks,  Int_t   charge,
                                                                             Int_t   bin,      Bool_t  decay,
                                                                             Int_t   outTOF,   Int_t   mcTag)

{
  // Fill mostly momentum imbalance related histograms
  
  fhXECharged         ->Fill(ptTrig , xE); 
  fhPtHbpXECharged    ->Fill(ptTrig , hbpXE);
  fhZTCharged         ->Fill(ptTrig , zT); 
  fhPtHbpZTCharged    ->Fill(ptTrig , hbpZT);
  fhPtTrigPout        ->Fill(ptTrig , pout) ;
  fhPtTrigCharged     ->Fill(ptTrig , ptAssoc) ;
  
  // Pile up studies
  if(fFillPileUpHistograms) 
  {
    if     (outTOF==1)
    {
      fhXEChargedOtherBC    ->Fill(ptTrig,xE);
      fhZTChargedOtherBC    ->Fill(ptTrig,zT);
      fhPtTrigChargedOtherBC->Fill(ptTrig,ptAssoc);
    }
    else if(outTOF==0)
    {
      fhXEChargedBC0    ->Fill(ptTrig,xE);
      fhZTChargedBC0    ->Fill(ptTrig,zT);
      fhPtTrigChargedBC0->Fill(ptTrig,ptAssoc);
    }

    Int_t vtxBC = GetReader()->GetVertexBC();
    if(vtxBC == 0 || vtxBC==AliVTrack::kTOFBCNA)
    {
      fhXEChargedVtxBC0    ->Fill(ptTrig,xE);
      fhZTChargedVtxBC0    ->Fill(ptTrig,zT);
      fhPtTrigChargedVtxBC0->Fill(ptTrig,ptAssoc);
    }
       
    if(GetReader()->IsPileUpFromSPD())                { fhXEChargedPileUp[0]->Fill(ptTrig,xE); fhZTChargedPileUp[0]->Fill(ptTrig,zT); fhPtTrigChargedPileUp[0]->Fill(ptTrig,ptAssoc); }
    if(GetReader()->IsPileUpFromEMCal())              { fhXEChargedPileUp[1]->Fill(ptTrig,xE); fhZTChargedPileUp[1]->Fill(ptTrig,zT); fhPtTrigChargedPileUp[1]->Fill(ptTrig,ptAssoc); }
    if(GetReader()->IsPileUpFromSPDOrEMCal())         { fhXEChargedPileUp[2]->Fill(ptTrig,xE); fhZTChargedPileUp[2]->Fill(ptTrig,zT); fhPtTrigChargedPileUp[2]->Fill(ptTrig,ptAssoc); }
    if(GetReader()->IsPileUpFromSPDAndEMCal())        { fhXEChargedPileUp[3]->Fill(ptTrig,xE); fhZTChargedPileUp[3]->Fill(ptTrig,zT); fhPtTrigChargedPileUp[3]->Fill(ptTrig,ptAssoc); }
    if(GetReader()->IsPileUpFromSPDAndNotEMCal())     { fhXEChargedPileUp[4]->Fill(ptTrig,xE); fhZTChargedPileUp[4]->Fill(ptTrig,zT); fhPtTrigChargedPileUp[4]->Fill(ptTrig,ptAssoc); }
    if(GetReader()->IsPileUpFromEMCalAndNotSPD())     { fhXEChargedPileUp[5]->Fill(ptTrig,xE); fhZTChargedPileUp[5]->Fill(ptTrig,zT); fhPtTrigChargedPileUp[5]->Fill(ptTrig,ptAssoc); }
    if(GetReader()->IsPileUpFromNotSPDAndNotEMCal())  { fhXEChargedPileUp[6]->Fill(ptTrig,xE); fhZTChargedPileUp[6]->Fill(ptTrig,zT); fhPtTrigChargedPileUp[6]->Fill(ptTrig,ptAssoc); }
  }
  
  if(IsDataMC())
  {
    Int_t mcIndex = GetMCTagHistogramIndex(mcTag);
    fhXEChargedMC      [mcIndex]->Fill(ptTrig , xE      ); 
  }  
  
  if(fDecayTrigger && decay)
  {          
    fhXEDecayCharged->Fill(ptTrig,xE); 
    fhZTDecayCharged->Fill(ptTrig,zT);
  } // photon decay pi0/eta trigger        
  
  if(bin >= 0 )//away side 
  {
    fhXEAssocPtBin[bin]->Fill(ptTrig, xE) ;
    fhZTAssocPtBin[bin]->Fill(ptTrig, zT) ;
    
    if(fDecayTrigger && decay)
    {          
      fhXEDecayChargedAssocPtBin[bin]->Fill(ptTrig, xE); 
      fhZTDecayChargedAssocPtBin[bin]->Fill(ptTrig, zT);
    }
  }        
  
  if(charge > 0)
  {
    fhXEPosCharged->Fill(ptTrig,xE) ;
    fhZTPosCharged->Fill(ptTrig,zT) ;
  }
  else  
  {
    fhXENegCharged->Fill(ptTrig,xE) ;
    fhZTNegCharged->Fill(ptTrig,zT) ;
  }
  
  //fill different multiplicity histogram
  if(DoEventSelect())
  {
    for(Int_t im=0; im<GetMultiBin(); im++)
    {
      if(nTracks < ( GetMaxMulti() - GetMinMulti() )/GetMultiBin()*(im+1))
      {
        fhTrigXECorr[im]->Fill(ptTrig,xE);
        fhTrigZTCorr[im]->Fill(ptTrig,zT);
      }
    }
  } //multiplicity events selection
} 

//_______________________________________________________________________________________________________________________
void AliAnaParticleHadronCorrelation::FillChargedUnderlyingEventHistograms(Float_t ptTrig,   Float_t ptAssoc,
                                                                           Float_t deltaPhi, Int_t nTracks, Int_t outTOF)
{
  // Fill underlying event histograms
  
  fhDeltaPhiUeChargedPt->Fill(ptAssoc,deltaPhi);
  
  Double_t randomphi = gRandom->Uniform(fDeltaPhiMinCut,fDeltaPhiMaxCut);
  Double_t uexE = -(ptAssoc/ptTrig)*TMath::Cos(randomphi);
  Double_t uezT =   ptAssoc/ptTrig;
  
  if(uexE < 0.) uexE = -uexE;
    
  fhXEUeCharged->Fill(ptTrig,uexE);
  if(uexE > 0) fhPtHbpXEUeCharged->Fill(ptTrig,TMath::Log(1/uexE));
  
  fhZTUeCharged->Fill(ptTrig,uezT);
  if(uexE > 0) fhPtHbpZTUeCharged->Fill(ptTrig,TMath::Log(1/uezT));
  
  // Pile up studies
  
  if(fFillPileUpHistograms)
  {
    if     (outTOF==1)
    {
      fhXEUeChargedOtherBC->Fill(ptTrig,uexE);
      fhZTUeChargedOtherBC->Fill(ptTrig,uezT);
    }
    else if(outTOF==0)
    {
      fhXEUeChargedBC0->Fill(ptTrig,uexE);
      fhZTUeChargedBC0->Fill(ptTrig,uezT);
    }
    
    Int_t vtxBC = GetReader()->GetVertexBC();
    if(vtxBC == 0 || vtxBC==AliVTrack::kTOFBCNA)
    {
      fhXEUeChargedVtxBC0->Fill(ptTrig,uexE);
      fhZTUeChargedVtxBC0->Fill(ptTrig,uezT);
    }

    if(GetReader()->IsPileUpFromSPD())               { fhXEUeChargedPileUp[0]->Fill(ptTrig,uexE); fhZTUeChargedPileUp[0]->Fill(ptTrig,uezT);}
    if(GetReader()->IsPileUpFromEMCal())             { fhXEUeChargedPileUp[1]->Fill(ptTrig,uexE); fhZTUeChargedPileUp[1]->Fill(ptTrig,uezT);}
    if(GetReader()->IsPileUpFromSPDOrEMCal())        { fhXEUeChargedPileUp[2]->Fill(ptTrig,uexE); fhZTUeChargedPileUp[2]->Fill(ptTrig,uezT);}
    if(GetReader()->IsPileUpFromSPDAndEMCal())       { fhXEUeChargedPileUp[3]->Fill(ptTrig,uexE); fhZTUeChargedPileUp[3]->Fill(ptTrig,uezT);}
    if(GetReader()->IsPileUpFromSPDAndNotEMCal())    { fhXEUeChargedPileUp[4]->Fill(ptTrig,uexE); fhZTUeChargedPileUp[4]->Fill(ptTrig,uezT);}
    if(GetReader()->IsPileUpFromEMCalAndNotSPD())    { fhXEUeChargedPileUp[5]->Fill(ptTrig,uexE); fhZTUeChargedPileUp[5]->Fill(ptTrig,uezT);}
    if(GetReader()->IsPileUpFromNotSPDAndNotEMCal()) { fhXEUeChargedPileUp[6]->Fill(ptTrig,uexE); fhZTUeChargedPileUp[6]->Fill(ptTrig,uezT);}
  }
  
  if(DoEventSelect())
  {
    for(Int_t im=0; im<GetMultiBin(); im++)
    {
      if(nTracks < ( GetMaxMulti() - GetMinMulti() )/GetMultiBin()*(im+1))
      {
        fhTrigXEUeCorr[im]->Fill(ptTrig,uexE); // xE? CHECK
        fhTrigZTUeCorr[im]->Fill(ptTrig,uezT); // zT? CHECK
      }
    }
  } //multiplicity events selection
}

//_____________________________________________________________________________________________________
void AliAnaParticleHadronCorrelation::FillChargedUnderlyingEventSidesHistograms(Float_t ptTrig, 
                                                                                Float_t ptAssoc, 
                                                                                Float_t deltaPhi)
{
 // Fill underlying event histograms to the left and right of trigger
  if((deltaPhi<-fUeDeltaPhiMinCut) || (deltaPhi >2*fUeDeltaPhiMaxCut))
  {  
    fhDeltaPhiUeLeftCharged->Fill(ptAssoc,deltaPhi); 
    Double_t randomphi = gRandom->Uniform(fDeltaPhiMinCut,fDeltaPhiMaxCut);
    Double_t uexE = -(ptAssoc/ptTrig)*TMath::Cos(randomphi);
    Double_t uezT =   ptAssoc/ptTrig;
  
    if(uexE < 0.) uexE = -uexE;
    
    fhXEUeLeftCharged->Fill(ptTrig,uexE);
    if(uexE > 0) fhPtHbpXEUeLeftCharged->Fill(ptTrig,TMath::Log(1/uexE));
  
    fhZTUeLeftCharged->Fill(ptTrig,uezT);
    if(uexE > 0) fhPtHbpZTUeLeftCharged->Fill(ptTrig,TMath::Log(1/uezT));
    fhDeltaPhiUeLeftCharged->Fill(ptAssoc, deltaPhi);
  }
  
  if((deltaPhi > fUeDeltaPhiMinCut) && (deltaPhi < fUeDeltaPhiMaxCut))
  {  
    fhDeltaPhiUeRightCharged->Fill(ptAssoc,deltaPhi); 
    Double_t randomphi = gRandom->Uniform(fDeltaPhiMinCut,fDeltaPhiMaxCut);
    Double_t uexE = -(ptAssoc/ptTrig)*TMath::Cos(randomphi);
    Double_t uezT =   ptAssoc/ptTrig;
  
    if(uexE < 0.) uexE = -uexE;
    
    fhXEUeRightCharged->Fill(ptTrig,uexE);
    if(uexE > 0) fhPtHbpXEUeRightCharged->Fill(ptTrig,TMath::Log(1/uexE));
  
    fhZTUeRightCharged->Fill(ptTrig,uezT);
    if(uexE > 0) fhPtHbpZTUeRightCharged->Fill(ptTrig,TMath::Log(1/uezT));
    fhDeltaPhiUeRightCharged->Fill(ptAssoc, deltaPhi);
  }

  if((deltaPhi<-fUeDeltaPhiMinCut) && (deltaPhi >-TMath::Pi()/2))
  {  
    fhDeltaPhiUeLeftDownCharged->Fill(ptAssoc,deltaPhi); 
    Double_t randomphi = gRandom->Uniform(fDeltaPhiMinCut,fDeltaPhiMaxCut);
    Double_t uexE = -(ptAssoc/ptTrig)*TMath::Cos(randomphi);
    
    if(uexE < 0.) uexE = -uexE;
    
    fhXEUeLeftDownCharged->Fill(ptTrig,uexE);
  }
  
  if((deltaPhi>2*fUeDeltaPhiMaxCut) && (deltaPhi <3*TMath::Pi()/2))
  {  
    fhDeltaPhiUeLeftUpCharged->Fill(ptAssoc,deltaPhi); 
    Double_t randomphi = gRandom->Uniform(fDeltaPhiMinCut,fDeltaPhiMaxCut);
    Double_t uexE = -(ptAssoc/ptTrig)*TMath::Cos(randomphi);
    
    if(uexE < 0.) uexE = -uexE;
    
    fhXEUeLeftUpCharged->Fill(ptTrig,uexE);
  }
  
  if((deltaPhi > TMath::Pi()/2) && (deltaPhi < fUeDeltaPhiMaxCut))
  {  
    fhDeltaPhiUeRightUpCharged->Fill(ptAssoc,deltaPhi); 
    Double_t randomphi = gRandom->Uniform(fDeltaPhiMinCut,fDeltaPhiMaxCut);
    Double_t uexE = -(ptAssoc/ptTrig)*TMath::Cos(randomphi);
    
    if(uexE < 0.) uexE = -uexE;
    
    fhXEUeRightUpCharged->Fill(ptTrig,uexE);
  }
  
  if((deltaPhi > fUeDeltaPhiMinCut) && (deltaPhi < TMath::Pi()/2))
  {  
    fhDeltaPhiUeRightDownCharged->Fill(ptAssoc,deltaPhi); 
    Double_t randomphi = gRandom->Uniform(fDeltaPhiMinCut,fDeltaPhiMaxCut);
    Double_t uexE = -(ptAssoc/ptTrig)*TMath::Cos(randomphi);
    
    if(uexE < 0.) uexE = -uexE;
    
    fhXEUeRightDownCharged->Fill(ptTrig,uexE);
  }  
} 

//______________________________________________________________________________________________________________________________
void AliAnaParticleHadronCorrelation::FillDecayPhotonCorrelationHistograms(Float_t ptAssoc,     Float_t phiAssoc, 
                                                                           TLorentzVector mom1, TLorentzVector mom2,
                                                                           Bool_t bChargedOrNeutral)
{
  // Do correlation with decay photons of triggered pi0 or eta
  
  // Calculate the correlation parameters
  Float_t ptDecay1 = mom1.Pt();
  Float_t ptDecay2 = mom2.Pt();
  
  Float_t zTDecay1 = -100, zTDecay2 = -100;
  if(ptDecay1) zTDecay1 = ptAssoc/ptDecay1 ;
  if(ptDecay2) zTDecay2 = ptAssoc/ptDecay2 ;
  
  Float_t deltaPhiDecay1 = mom1.Phi()-phiAssoc;
  if(deltaPhiDecay1< -TMath::PiOver2()) deltaPhiDecay1+=TMath::TwoPi();
  if(deltaPhiDecay1>3*TMath::PiOver2()) deltaPhiDecay1-=TMath::TwoPi();
  
  Float_t deltaPhiDecay2 = mom2.Phi()-phiAssoc;
  if(deltaPhiDecay2< -TMath::PiOver2()) deltaPhiDecay2+=TMath::TwoPi();
  if(deltaPhiDecay2>3*TMath::PiOver2()) deltaPhiDecay2-=TMath::TwoPi();
  
  Float_t xEDecay1   =-zTDecay1*TMath::Cos(deltaPhiDecay1); // -(px*pxTrig+py*pyTrig)/(ptTrig*ptTrig);
  Float_t xEDecay2   =-zTDecay2*TMath::Cos(deltaPhiDecay2); // -(px*pxTrig+py*pyTrig)/(ptTrig*ptTrig);
  
  if(bChargedOrNeutral) // correlate with charges
  {
    fhDeltaPhiDecayCharged->Fill(ptDecay1, deltaPhiDecay1);
    fhDeltaPhiDecayCharged->Fill(ptDecay2, deltaPhiDecay2);
    
    if(GetDebug() > 1) printf("AliAnaParticleHadronCorrelation::FillDecayPhotonHistograms( Charged corr) - deltaPhoton1 = %f, deltaPhoton2 = %f \n", deltaPhiDecay1, deltaPhiDecay2);
    
    if( (deltaPhiDecay1 > fDeltaPhiMinCut) && ( deltaPhiDecay1 < fDeltaPhiMaxCut) )
    {
      fhZTDecayCharged->Fill(ptDecay1,zTDecay1); 
      fhXEDecayCharged->Fill(ptDecay1,xEDecay1); 
    }
    if( (deltaPhiDecay2 > fDeltaPhiMinCut) && ( deltaPhiDecay2 < fDeltaPhiMaxCut) )
    {
      fhZTDecayCharged->Fill(ptDecay2,zTDecay2);
      fhXEDecayCharged->Fill(ptDecay2,xEDecay2);
    }
  }
  else // correlate with neutrals
  {
    fhDeltaPhiDecayCharged->Fill(ptDecay1, deltaPhiDecay1);
    fhDeltaPhiDecayCharged->Fill(ptDecay2, deltaPhiDecay2);
    
    if(GetDebug() > 1) printf("AliAnaParticleHadronCorrelation::FillDecayPhotonHistograms(Neutral corr) - deltaPhoton1 = %f, deltaPhoton2 = %f \n", deltaPhiDecay1, deltaPhiDecay2);
    
    if( (deltaPhiDecay1 > fDeltaPhiMinCut) && ( deltaPhiDecay1 < fDeltaPhiMaxCut) )
    {
      fhZTDecayCharged->Fill(ptDecay1,zTDecay1); 
      fhXEDecayCharged->Fill(ptDecay1,xEDecay1); 
    }
    if( (deltaPhiDecay2 > fDeltaPhiMinCut) && ( deltaPhiDecay2 < fDeltaPhiMaxCut) )
    {
      fhZTDecayCharged->Fill(ptDecay2,zTDecay2);
      fhXEDecayCharged->Fill(ptDecay2,xEDecay2);
    }    
  }
}

//______________________________________________________________________________________________________________________________________________________
void AliAnaParticleHadronCorrelation::FillNeutralAngularCorrelationHistograms(Float_t ptAssoc,  Float_t ptTrig,  
                                                                              Float_t phiAssoc, Float_t phiTrig,  Float_t &     deltaPhi,
                                                                              Float_t etaAssoc, Float_t etaTrig)
{
  // Fill angular correlation related histograms
  
  Float_t deltaEta    = etaTrig-etaAssoc;
  deltaPhi    = phiTrig-phiAssoc;
  
  if(deltaPhi <= -TMath::PiOver2()) deltaPhi+=TMath::TwoPi();
  if(deltaPhi > 3*TMath::PiOver2()) deltaPhi-=TMath::TwoPi();
  
  fhEtaNeutral     ->Fill(ptAssoc,etaAssoc);
  fhPhiNeutral     ->Fill(ptAssoc,phiAssoc);
  fhDeltaEtaNeutral->Fill(ptTrig ,deltaEta);
  fhDeltaPhiNeutral->Fill(ptTrig ,deltaPhi);
  
  if(ptAssoc > 2 ) fhDeltaPhiDeltaEtaNeutral->Fill(deltaPhi, deltaEta);
  
}

//_____________________________________________________________________________________________________________________________
void AliAnaParticleHadronCorrelation::FillNeutralUnderlyingEventSidesHistograms(Float_t ptTrig,   Float_t ptAssoc, 
                                                                                Float_t xE,       Float_t hbpXE, 
                                                                                Float_t zT,       Float_t hbpZT, 
                                                                                Float_t deltaPhi)
{
  // Fill underlying event histograms to the left and right of trigger
  
  if((deltaPhi<-fUeDeltaPhiMinCut) && (deltaPhi >-fUeDeltaPhiMaxCut))
  {  
    fhDeltaPhiUeLeftNeutral->Fill(ptAssoc, deltaPhi);
    fhXEUeLeftNeutral      ->Fill(ptTrig , xE);
    fhPtHbpXEUeLeftNeutral ->Fill(ptTrig , hbpXE);
    fhZTUeLeftNeutral      ->Fill(ptTrig , zT);
    fhPtHbpZTUeLeftNeutral ->Fill(ptTrig , hbpZT);
  }
  
  if((deltaPhi > fUeDeltaPhiMinCut) && (deltaPhi < fUeDeltaPhiMaxCut))
  {  
    fhDeltaPhiUeRightNeutral->Fill(ptAssoc, deltaPhi);
    fhXEUeRightNeutral      ->Fill(ptTrig , xE);
    fhPtHbpXEUeRightNeutral ->Fill(ptTrig , hbpXE);
    fhZTUeRightNeutral      ->Fill(ptTrig , zT);
    fhPtHbpZTUeRightNeutral ->Fill(ptTrig , hbpZT);
  }
} 

//______________________________________________________
void AliAnaParticleHadronCorrelation::FillEventMixPool()
{
  // Fill the pool with tracks if requested
    
  if(DoOwnMix())
  {
    FillChargedEventMixPool();
    
    if(OnlyIsolated() || fFillNeutralEventMixPool)
      FillNeutralEventMixPool();
  }
  
}

//_____________________________________________________________
void AliAnaParticleHadronCorrelation::FillChargedEventMixPool()
{
  // Mixed event pool filling for tracks
    
  if(fUseMixStoredInReader && GetReader()->GetLastTracksMixedEvent() == GetEventNumber())
  {
    //printf("%s : Pool already filled for this event !!!\n",GetInputAODName().Data());
    return ; // pool filled previously for another trigger
  }
  
  Int_t nTracks = GetCTSTracks()->GetEntriesFast();
    
  AliAnalysisManager   * manager      = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler * inputHandler = dynamic_cast<AliInputEventHandler*>(manager->GetInputEventHandler());
  
  if(!inputHandler) return ;
    
  // Do mixing only with MB event (or the chosen mask), if not skip
  if( !(inputHandler->IsEventSelected( ) & GetReader()->GetMixEventTriggerMask()) ) return ;
  
  fhNtracksMB->Fill(nTracks);
  
  Int_t eventBin = GetEventMixBin();
  
  //Check that the bin exists, if not (bad determination of RP, centrality or vz bin) do nothing
  if(eventBin < 0) return;
  
  TObjArray * mixEventTracks = new TObjArray;
  
  if(fUseMixStoredInReader) 
  {
    fListMixTrackEvents[eventBin] = GetReader()->GetListWithMixedEventsForTracks(eventBin);
  }
  
  if(!fListMixTrackEvents[eventBin]) fListMixTrackEvents[eventBin] = new TList();
  
  //printf("%s ***** Pool Event bin : %d - nTracks %d\n",GetInputAODName().Data(),eventBin, GetCTSTracks()->GetEntriesFast());
  
  TList * pool = fListMixTrackEvents[eventBin];
  
  TVector3 p3;  
  for(Int_t ipr = 0;ipr < GetCTSTracks()->GetEntriesFast() ; ipr ++ )
  {
    AliVTrack * track = (AliVTrack *) (GetCTSTracks()->At(ipr)) ;
    
    Double_t mom[3] = {track->Px(),track->Py(),track->Pz()};
    p3.SetXYZ(mom[0],mom[1],mom[2]);
    Float_t pt   = p3.Pt();
    
    //Select only hadrons in pt range
    if(pt < fMinAssocPt || pt > fMaxAssocPt) continue ;
    
    AliAODPWG4Particle * mixedTrack = new AliAODPWG4Particle(mom[0],mom[1],mom[2],0);
    mixedTrack->SetDetector("CTS");
    mixedTrack->SetChargedBit(track->Charge()>0);
    mixEventTracks->Add(mixedTrack);
  }
  
  //Set the event number where the last event was added, to avoid double pool filling
  GetReader()->SetLastTracksMixedEvent(GetEventNumber());
  
  //printf("Add event to pool with %d tracks \n ",mixEventTracks->GetEntries());
  pool->AddFirst(mixEventTracks);
  mixEventTracks = 0;
  
  //printf("Pool size %d, max %d\n",pool->GetSize(), GetNMaxEvMix());

  if(pool->GetSize() > GetNMaxEvMix())
  {//Remove last event
    TClonesArray * tmp = static_cast<TClonesArray*>(pool->Last()) ;
    pool->RemoveLast() ;
    delete tmp ;
  }
}

//_____________________________________________________________
void AliAnaParticleHadronCorrelation::FillNeutralEventMixPool()
{
  // Mixed event pool filling for neutral clusters
  // Right now only for EMCAL and in isolation case
  
  //printf("FillNeutralEventMixPool for %s\n",GetInputAODName().Data());
  
  TObjArray * pl = GetEMCALClusters();
  //if (GetAODObjArrayName.Contains("PHOS") )pl    = GetPHOSClusters();
  //else                                     pl    = GetEMCALClusters();
  
  Int_t nClusters   = pl->GetEntriesFast();
  
  if(fUseMixStoredInReader && GetReader()->GetLastCaloMixedEvent() == GetEventNumber())
  {
    //printf("%s : Pool already filled for this event !!!\n",GetInputAODName().Data());
    return ; // pool filled previously for another trigger
  }
  
  AliAnalysisManager   * manager      = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler * inputHandler = dynamic_cast<AliInputEventHandler*>(manager->GetInputEventHandler());
  
  if(!inputHandler) return ;
    
  // Do mixing only with MB event (or the chosen mask), if not skip
  if( !(inputHandler->IsEventSelected( ) & GetReader()->GetMixEventTriggerMask()) ) return ;
  
  fhNclustersMB->Fill(nClusters);
  
  Int_t eventBin = GetEventMixBin();
  
  //Check that the bin exists, if not (bad determination of RP, centrality or vz bin) do nothing
  if(eventBin < 0) return;
  
  TObjArray * mixEventCalo = new TObjArray;
  
  if(fUseMixStoredInReader) 
  {
    fListMixCaloEvents[eventBin] = GetReader()->GetListWithMixedEventsForCalo(eventBin);
  }
  
  if(!fListMixCaloEvents[eventBin]) fListMixCaloEvents[eventBin] = new TList();
  
  TList * poolCalo = fListMixCaloEvents[eventBin];
  
  TLorentzVector mom;

  for(Int_t ipr = 0;ipr <  nClusters ; ipr ++ )
  {
    AliVCluster * calo = (AliVCluster *) (pl->At(ipr)) ;
  
    // remove matched clusters
    if( IsTrackMatched( calo, GetReader()->GetInputEvent() ) ) continue ;
    
    //Cluster momentum calculation
    if(GetReader()->GetDataType() != AliCaloTrackReader::kMC)
    {
      calo->GetMomentum(mom,GetVertex(0)) ;
    }//Assume that come from vertex in straight line
    else
    {
      Double_t vertex[]={0,0,0};
      calo->GetMomentum(mom,vertex) ;
    }
    
    Float_t pt = mom.Pt();
    //Select only clusters in pt range
    if(pt < fMinAssocPt || pt > fMaxAssocPt) continue ;
    
    AliAODPWG4Particle * mixedCalo = new AliAODPWG4Particle(mom);
    mixedCalo->SetDetector("EMCAL");
    mixEventCalo->Add(mixedCalo);
  }
  
  //Set the event number where the last event was added, to avoid double pool filling
  GetReader()->SetLastCaloMixedEvent(GetEventNumber());
  
  //printf("Add event to pool with %d clusters \n ",mixEventCalo->GetEntries());
  poolCalo->AddFirst(mixEventCalo);
  mixEventCalo = 0;
  
  //printf("Pool size %d, max %d\n",poolCalo->GetSize(), GetNMaxEvMix());
  
  if(poolCalo->GetSize() > GetNMaxEvMix())
  {//Remove last event
    TClonesArray * tmp = static_cast<TClonesArray*>(poolCalo->Last()) ;
    poolCalo->RemoveLast() ;
    delete tmp ;
  }  
}

//____________________________________________________________
TObjString* AliAnaParticleHadronCorrelation::GetAnalysisCuts()
{
  //Save parameters used for analysis
  TString parList ; //this will be list of parameters used for this analysis.
  const Int_t buffersize = 560;
  char onePar[buffersize] ;
  
  snprintf(onePar,buffersize,"--- AliAnaPaticleHadronCorrelation ---\n") ;
  parList+=onePar ;	
  snprintf(onePar,buffersize," Pt Trigger > %3.2f ",    fMinTriggerPt) ; 
  parList+=onePar ;
  snprintf(onePar,buffersize," %3.2f < Pt associated < %3.2f ", fMinAssocPt,   fMaxAssocPt) ; 
  parList+=onePar ;
  snprintf(onePar,buffersize," %3.2f < Phi trigger particle-Hadron < %3.2f ",    fDeltaPhiMinCut,   fDeltaPhiMaxCut) ; 
  parList+=onePar ;
  snprintf(onePar,buffersize," %3.2f < Phi trigger particle-UeHadron <  %3.2f ", fUeDeltaPhiMinCut, fUeDeltaPhiMaxCut) ; 
  parList+=onePar ;
  snprintf(onePar,buffersize,"Isolated Trigger?  %d\n",    fSelectIsolated) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"Several UE?  %d\n",          fMakeSeveralUE) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"Name of AOD Pi0 Branch %s ", fPi0AODBranchName.Data());
  parList+=onePar ;
  snprintf(onePar,buffersize,"Do Decay-hadron correlation ?  pi0 %d, decay %d", fPi0Trigger, fDecayTrigger) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"Select absolute leading for cluster triggers ? %d or Near Side Leading %d \n", 
           fMakeAbsoluteLeading, fMakeNearSideLeading) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"Associated particle pt bins  %d: ", fNAssocPtBins) ;
  parList+=onePar ;
  for (Int_t ibin = 0; ibin<fNAssocPtBins; ibin++) {
    snprintf(onePar,buffersize,"bin %d = [%2.1f,%2.1f];", ibin, fAssocPtBinLimit[ibin], fAssocPtBinLimit[ibin+1]) ;
  }
  parList+=onePar ;
  
  //Get parameters set in base class.
  parList += GetBaseParametersList() ;
  
  //Get parameters set in FiducialCut class (not available yet)
  //parlist += GetFidCut()->GetFidCutParametersList() 
  
  return new TObjString(parList) ;  
  
} 

//________________________________________________________________
TList *  AliAnaParticleHadronCorrelation::GetCreateOutputObjects()
{  
  
  // Create histograms to be saved in output file and 
  // store them in fOutputContainer
  
  TList * outputContainer = new TList() ; 
  outputContainer->SetName("CorrelationHistos") ; 
  
  Int_t   nptbins = GetHistogramRanges()->GetHistoPtBins(); Int_t  nphibins = GetHistogramRanges()->GetHistoPhiBins(); Int_t   netabins = GetHistogramRanges()->GetHistoEtaBins(); Int_t  ndeltaphibins = GetHistogramRanges()->GetHistoDeltaPhiBins(); Int_t   ndeltaetabins = GetHistogramRanges()->GetHistoDeltaEtaBins();
  Float_t ptmax   = GetHistogramRanges()->GetHistoPtMax();  Float_t phimax  = GetHistogramRanges()->GetHistoPhiMax();  Float_t etamax   = GetHistogramRanges()->GetHistoEtaMax(); Float_t deltaphimax  = GetHistogramRanges()->GetHistoDeltaPhiMax();  Float_t deltaetamax   = GetHistogramRanges()->GetHistoDeltaEtaMax();
  Float_t ptmin   = GetHistogramRanges()->GetHistoPtMin();  Float_t phimin  = GetHistogramRanges()->GetHistoPhiMin();  Float_t etamin   = GetHistogramRanges()->GetHistoEtaMin(); Float_t deltaphimin  = GetHistogramRanges()->GetHistoDeltaPhiMin();  Float_t deltaetamin   = GetHistogramRanges()->GetHistoDeltaEtaMin();	
  
  Int_t nMixBins = GetNCentrBin()*GetNZvertBin()*GetNRPBin();
  
  TString nameMC[]     = {"Photon","Pi0","Pi0Decay","EtaDecay","OtherDecay","Electron","Hadron"};
  TString pileUpName[] = {"SPD","EMCAL","SPDOrEMCAL","SPDAndEMCAL","SPDAndNotEMCAL","EMCALAndNotSPD","NotSPDAndNotEMCAL"} ;

  // For vz dependent histograms, if option ON
  Int_t   nz  = 1  ;
  if(fCorrelVzBin) nz = GetNZvertBin();
  TString sz  = "" ;
  TString tz  = "" ;
  
  fhPtInput  = new TH1F("hPtInput","p_{T} distribution of input trigger particles", nptbins,ptmin,ptmax); 
  fhPtInput->SetXTitle("p_{T}^{trig} (GeV/c)");
  outputContainer->Add(fhPtInput);

  fhPtFidCut  = new TH1F("hPtFidCut","p_{T} distribution of input trigger particles after fiducial cut", nptbins,ptmin,ptmax); 
  fhPtFidCut->SetXTitle("p_{T}^{trig} (GeV/c)");
  outputContainer->Add(fhPtFidCut);

  fhPtLeading  = new TH1F("hPtLeading","p_{T} distribution of leading particles", nptbins,ptmin,ptmax); 
  fhPtLeading->SetXTitle("p_{T}^{trig} (GeV/c)");
  outputContainer->Add(fhPtLeading);

  if(IsDataMC())
  {
    for(Int_t i=0; i < 7; i++)
    {
      fhPtLeadingMC[i]  = new TH1F(Form("hPtLeading_MC%s",nameMC[i].Data()),
                                   Form("p_{T} distribution of leading particles, trigger origin is %s",nameMC[i].Data()), 
                                   nptbins,ptmin,ptmax); 
      fhPtLeadingMC[i]->SetXTitle("p_{T}^{trig} (GeV/c)");
      outputContainer->Add(fhPtLeadingMC[i]);
    }
  }
  
  if(fCorrelVzBin)
  {
    fhPtLeadingVzBin  = new TH2F("hPtLeadingVzBin","p_{T} distribution of leading particles vs vz bin", nptbins,ptmin,ptmax,GetNZvertBin(),0,GetNZvertBin()); 
    fhPtLeadingVzBin->SetXTitle("p_{T}^{trig} (GeV/c)");
    fhPtLeadingVzBin->SetYTitle("v_{z} bin");
    outputContainer->Add(fhPtLeadingVzBin);
  }
  
  fhPtLeadingBin  = new TH2F ("hPtLeadingBin","p_{T} distribution of leading particles", nptbins,ptmin,ptmax,nMixBins,0,nMixBins); 
  fhPtLeadingBin->SetXTitle("p_{T}^{trig} (GeV/c)");
  fhPtLeadingBin->SetYTitle("Bin");
  outputContainer->Add(fhPtLeadingBin);

  fhPhiLeading  = new TH2F ("hPhiLeading","#phi distribution of leading Particles",nptbins,ptmin,ptmax, nphibins,phimin,phimax); 
  fhPhiLeading->SetYTitle("#phi (rad)");
  outputContainer->Add(fhPhiLeading);

  fhEtaLeading  = new TH2F ("hEtaLeading","#eta distribution of leading",nptbins,ptmin,ptmax, netabins,etamin,etamax); 
  fhEtaLeading->SetYTitle("#eta ");  
  outputContainer->Add(fhEtaLeading);
  
  fhPtLeadingCentrality   = new TH2F("hPtLeadingCentrality","Leading particle p_{T} vs centrality",nptbins,ptmin,ptmax,100,0.,100) ;
  fhPtLeadingCentrality->SetXTitle("p_{T}^{trig} (GeV/c)");
  fhPtLeadingCentrality->SetYTitle("Centrality (%)");
  outputContainer->Add(fhPtLeadingCentrality) ;  
  
  fhPtLeadingEventPlane  = new TH2F("hPtLeadingEventPlane","Leading particle p_{T} vs event plane angle",nptbins,ptmin,ptmax, 100,0.,TMath::Pi()) ;
  fhPtLeadingEventPlane->SetXTitle("p_{T}^{trig} (GeV/c)");
  fhPtLeadingEventPlane->SetXTitle("EP angle (rad)");
  outputContainer->Add(fhPtLeadingEventPlane) ;
  
  fhLeadingEventPlaneCentrality  = new TH2F("hLeadingEventPlane","Leading particle centrality vs event plane angle",100,0.,100,100,0.,TMath::Pi()) ;
  fhLeadingEventPlaneCentrality->SetXTitle("Centrality (%)");
  fhLeadingEventPlaneCentrality->SetYTitle("EP angle (rad)");
  outputContainer->Add(fhLeadingEventPlaneCentrality) ;
  
  //Correlation with charged hadrons
  if(GetReader()->IsCTSSwitchedOn()) 
  {
    fhDeltaPhiDeltaEtaCharged  = new TH2F
    ("hDeltaPhiDeltaEtaCharged","#eta_{trigger} - #eta_{h^{#pm}} vs #phi_{trigger} - #phi_{h^{#pm}}",
    ndeltaphibins ,deltaphimin,deltaphimax,ndeltaetabins,deltaetamin,deltaetamax); 
    fhDeltaPhiDeltaEtaCharged->SetXTitle("#Delta #phi (rad)");
    fhDeltaPhiDeltaEtaCharged->SetYTitle("#Delta #eta");    
    
    fhDeltaPhiDeltaEtaChargedPtA3GeV  = new TH2F
    ("hDeltaPhiDeltaEtaChargedPtA3GeV","#eta_{trigger} - #eta_{h^{#pm}} vs #phi_{trigger} - #phi_{h^{#pm}, p_{TA}>3 GeV/c}",
     ndeltaphibins ,deltaphimin,deltaphimax,ndeltaetabins,deltaetamin,deltaetamax); 
    fhDeltaPhiDeltaEtaChargedPtA3GeV->SetXTitle("#Delta #phi (rad)");
    fhDeltaPhiDeltaEtaChargedPtA3GeV->SetYTitle("#Delta #eta");    
        
    fhPhiCharged  = new TH2F
    ("hPhiCharged","#phi_{h^{#pm}}  vs p_{T #pm}",
     nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
    fhPhiCharged->SetYTitle("#phi_{h^{#pm}} (rad)");
    fhPhiCharged->SetXTitle("p_{T #pm} (GeV/c)");
    
    fhEtaCharged  = new TH2F
    ("hEtaCharged","#eta_{h^{#pm}}  vs p_{T #pm}",
     nptbins,ptmin,ptmax,netabins,etamin,etamax); 
    fhEtaCharged->SetYTitle("#eta_{h^{#pm}} (rad)");
    fhEtaCharged->SetXTitle("p_{T #pm} (GeV/c)");
    
    fhDeltaPhiCharged  = new TH2F
    ("hDeltaPhiCharged","#phi_{trigger} - #phi_{h^{#pm}} vs p_{T trigger}",
     nptbins,ptmin,ptmax, ndeltaphibins ,deltaphimin,deltaphimax);
    fhDeltaPhiCharged->SetYTitle("#Delta #phi (rad)");
    fhDeltaPhiCharged->SetXTitle("p_{T trigger} (GeV/c)");
    
    fhDeltaPhiChargedPtA3GeV  = new TH2F
    ("hDeltaPhiChargedPtA3GeV","#phi_{trigger} - #phi_{h^{#pm}} vs p_{T trigger}, p_{TA}>3 GeV/c",
     nptbins,ptmin,ptmax, ndeltaphibins ,deltaphimin,deltaphimax);
    fhDeltaPhiChargedPtA3GeV->SetYTitle("#Delta #phi (rad)");
    fhDeltaPhiChargedPtA3GeV->SetXTitle("p_{T trigger} (GeV/c)");
    

    fhDeltaPhiChargedPt  = new TH2F
    ("hDeltaPhiChargedPt","#phi_{trigger} - #phi_{#h^{#pm}} vs p_{T h^{#pm}}",
     nptbins,ptmin,ptmax, ndeltaphibins ,deltaphimin,deltaphimax);
    fhDeltaPhiChargedPt->SetYTitle("#Delta #phi (rad)");
    fhDeltaPhiChargedPt->SetXTitle("p_{T h^{#pm}} (GeV/c)");
    
    fhDeltaPhiUeChargedPt  = new TH2F
    ("hDeltaPhiUeChargedPt","#phi_{trigger} - #phi_{#Ueh^{#pm}} vs p_{T Ueh^{#pm}}",
     nptbins,ptmin,ptmax, ndeltaphibins ,deltaphimin,deltaphimax);
    fhDeltaPhiUeChargedPt->SetYTitle("#Delta #phi (rad)");
    fhDeltaPhiUeChargedPt->SetXTitle("p_{T h^{#pm}} (GeV/c)");
    
    fhUePart  =  new TH1F("hUePart","UE particles distribution vs pt trig",
             nptbins,ptmin,ptmax); 
    fhUePart->SetYTitle("dNch");
    fhUePart->SetXTitle("p_{T trigger} (GeV/c)");
    
    
    fhDeltaEtaCharged  = new TH2F
    ("hDeltaEtaCharged","#eta_{trigger} - #eta_{h^{#pm}} vs p_{T trigger}",
     nptbins,ptmin,ptmax,ndeltaetabins,deltaetamin,deltaetamax);  
    fhDeltaEtaCharged->SetYTitle("#Delta #eta");
    fhDeltaEtaCharged->SetXTitle("p_{T trigger} (GeV/c)");
    
    fhDeltaEtaChargedPtA3GeV  = new TH2F
    ("hDeltaEtaChargedPtA3GeV","#eta_{trigger} - #eta_{h^{#pm}} vs p_{T trigger}, p_{TA}>3 GeV/c",
     nptbins,ptmin,ptmax,ndeltaetabins,deltaetamin,deltaetamax);  
    fhDeltaEtaChargedPtA3GeV->SetYTitle("#Delta #eta");
    fhDeltaEtaChargedPtA3GeV->SetXTitle("p_{T trigger} (GeV/c)");    
    
    fhXECharged  = 
    new TH2F("hXECharged","x_{E} for charged tracks",
             nptbins,ptmin,ptmax,200,0.,2.); 
    fhXECharged->SetYTitle("x_{E}");
    fhXECharged->SetXTitle("p_{T trigger} (GeV/c)");
    
    fhXEUeCharged  = 
    new TH2F("hXEUeCharged","x_{E} for Underlying Event",
             nptbins,ptmin,ptmax,200,0.,2.); 
    fhXEUeCharged->SetYTitle("x_{E}");
    fhXEUeCharged->SetXTitle("p_{T trigger} (GeV/c)");
    
    fhXEPosCharged  = 
    new TH2F("hXEPositiveCharged","x_{E} for positive charged tracks",
             nptbins,ptmin,ptmax,200,0.,2.); 
    fhXEPosCharged->SetYTitle("x_{E}");
    fhXEPosCharged->SetXTitle("p_{T trigger} (GeV/c)");
    
    fhXENegCharged  = 
    new TH2F("hXENegativeCharged","x_{E} for negative charged tracks",
             nptbins,ptmin,ptmax,200,0.,2.); 
    fhXENegCharged->SetYTitle("x_{E}");
    fhXENegCharged->SetXTitle("p_{T trigger} (GeV/c)");
    
    fhPtHbpXECharged  = 
    new TH2F("hHbpXECharged","#xi = ln(1/x_{E}) with charged hadrons",
             nptbins,ptmin,ptmax,200,0.,10.); 
    fhPtHbpXECharged->SetYTitle("ln(1/x_{E})");
    fhPtHbpXECharged->SetXTitle("p_{T trigger} (GeV/c)");
    
    fhPtHbpXEUeCharged  = 
    new TH2F("hHbpXEUeCharged","#xi = ln(1/x_{E}) with charged hadrons,Underlying Event",
             nptbins,ptmin,ptmax,200,0.,10.); 
    fhPtHbpXEUeCharged->SetYTitle("ln(1/x_{E})");
    fhPtHbpXEUeCharged->SetXTitle("p_{T trigger} (GeV/c)");
    
    fhZTCharged  = 
    new TH2F("hZTCharged","z_{T} for charged tracks",
             nptbins,ptmin,ptmax,200,0.,2.); 
    fhZTCharged->SetYTitle("z_{T}");
    fhZTCharged->SetXTitle("p_{T trigger}");
    
    fhZTUeCharged  = 
    new TH2F("hZTUeCharged","z_{T} for Underlying Event",
             nptbins,ptmin,ptmax,200,0.,2.); 
    fhZTUeCharged->SetYTitle("z_{T}");
    fhZTUeCharged->SetXTitle("p_{T trigger} (GeV/c)");
    
    fhZTPosCharged  = 
    new TH2F("hZTPositiveCharged","z_{T} for positive charged tracks",
             nptbins,ptmin,ptmax,200,0.,2.); 
    fhZTPosCharged->SetYTitle("z_{T}");
    fhZTPosCharged->SetXTitle("p_{T trigger} (GeV/c)");
    
    fhZTNegCharged  = 
    new TH2F("hZTNegativeCharged","z_{T} for negative charged tracks",
             nptbins,ptmin,ptmax,200,0.,2.); 
    fhZTNegCharged->SetYTitle("z_{T}");
    fhZTNegCharged->SetXTitle("p_{T trigger} (GeV/c)");
    
    fhPtHbpZTCharged  = 
    new TH2F("hHbpZTCharged","#xi = ln(1/z_{T}) with charged hadrons",
             nptbins,ptmin,ptmax,200,0.,10.); 
    fhPtHbpZTCharged->SetYTitle("ln(1/z_{T})");
    fhPtHbpZTCharged->SetXTitle("p_{T trigger} (GeV/c)");
    
    fhPtHbpZTUeCharged  = 
    new TH2F("hHbpZTUeCharged","#xi = ln(1/z_{T}) with charged hadrons,Underlying Event",
             nptbins,ptmin,ptmax,200,0.,10.); 
    fhPtHbpZTUeCharged->SetYTitle("ln(1/x_{E})");
    fhPtHbpZTUeCharged->SetXTitle("p_{T trigger} (GeV/c)");
    
    fhPtTrigPout  = 
    new TH2F("hPtTrigPout","Pout with triggers",
             nptbins,ptmin,ptmax,2*nptbins,-ptmax,ptmax); 
    fhPtTrigPout->SetYTitle("p_{out} (GeV/c)");
    fhPtTrigPout->SetXTitle("p_{T trigger} (GeV/c)"); 
    
    fhPtTrigCharged  = 
    new TH2F("hPtTrigCharged","trigger and charged tracks pt distribution",
             nptbins,ptmin,ptmax,nptbins,ptmin,ptmax); 
    fhPtTrigCharged->SetYTitle("p_{T h^{#pm}} (GeV/c)");
    fhPtTrigCharged->SetXTitle("p_{T trigger} (GeV/c)");    
	  
    outputContainer->Add(fhDeltaPhiDeltaEtaCharged);
    outputContainer->Add(fhDeltaPhiDeltaEtaChargedPtA3GeV);
    outputContainer->Add(fhPhiCharged) ;
    outputContainer->Add(fhEtaCharged) ;
    outputContainer->Add(fhDeltaPhiCharged) ; 
    outputContainer->Add(fhDeltaPhiChargedPtA3GeV) ; 
    outputContainer->Add(fhDeltaEtaCharged) ;
    outputContainer->Add(fhDeltaEtaChargedPtA3GeV) ;
    outputContainer->Add(fhDeltaPhiChargedPt) ;
    outputContainer->Add(fhDeltaPhiUeChargedPt) ;
    outputContainer->Add(fhUePart);

    outputContainer->Add(fhXECharged) ;
  
    if(IsDataMC())
    {
      for(Int_t i=0; i < 7; i++)
      {
        
        fhDeltaPhiChargedMC[i]  = new TH2F(Form("hDeltaPhiCharged_MC%s",nameMC[i].Data()),
                                     Form("#Delta #phi for charged tracks, trigger origin is %s",nameMC[i].Data()),
                                     nptbins,ptmin,ptmax,ndeltaphibins ,deltaphimin,deltaphimax); 
        fhDeltaPhiChargedMC[i]->SetYTitle("x_{E}");
        fhDeltaPhiChargedMC[i]->SetXTitle("p_{T trigger} (GeV/c)");
        outputContainer->Add(fhDeltaPhiChargedMC[i]) ;
        
        fhXEChargedMC[i]  = new TH2F(Form("hXECharged_MC%s",nameMC[i].Data()),
                                     Form("x_{E} for charged tracks, trigger origin is %s",nameMC[i].Data()),
         nptbins,ptmin,ptmax,200,0.,2.); 
        fhXEChargedMC[i]->SetYTitle("x_{E}");
        fhXEChargedMC[i]->SetXTitle("p_{T trigger} (GeV/c)");
        outputContainer->Add(fhXEChargedMC[i]) ;
      }
    }
  
    outputContainer->Add(fhXEPosCharged) ;
    outputContainer->Add(fhXENegCharged) ;
    outputContainer->Add(fhXEUeCharged) ;
    outputContainer->Add(fhPtHbpXECharged) ;
    outputContainer->Add(fhPtHbpXEUeCharged) ;

    outputContainer->Add(fhZTCharged) ;
    outputContainer->Add(fhZTPosCharged) ;
    outputContainer->Add(fhZTNegCharged) ;
    outputContainer->Add(fhZTUeCharged) ;
    outputContainer->Add(fhPtHbpZTCharged) ;
    outputContainer->Add(fhPtHbpZTUeCharged) ;
    
    outputContainer->Add(fhPtTrigPout) ;
    outputContainer->Add(fhPtTrigCharged) ;
    
    if(fFillPileUpHistograms)
    {
      fhDeltaPhiChargedOtherBC  = new TH2F
      ("hDeltaPhiChargedOtherBC","#phi_{trigger} - #phi_{h^{#pm}} vs p_{T trigger}, track BC!=0",
       nptbins,ptmin,ptmax, ndeltaphibins ,deltaphimin,deltaphimax);
      fhDeltaPhiChargedOtherBC->SetYTitle("#Delta #phi (rad)");
      fhDeltaPhiChargedOtherBC->SetXTitle("p_{T trigger} (GeV/c)");
      
      fhDeltaPhiChargedPtA3GeVOtherBC  = new TH2F
      ("hDeltaPhiChargedPtA3GeVOtherBC","#phi_{trigger} - #phi_{h^{#pm}} vs p_{T trigger}, p_{TA}>3 GeV/c, track BC!=0",
       nptbins,ptmin,ptmax, ndeltaphibins ,deltaphimin,deltaphimax);
      fhDeltaPhiChargedPtA3GeVOtherBC->SetYTitle("#Delta #phi (rad)");
      fhDeltaPhiChargedPtA3GeVOtherBC->SetXTitle("p_{T trigger} (GeV/c)");
      
      fhPtTrigChargedOtherBC  =
      new TH2F("hPtTrigChargedOtherBC","trigger and charged tracks pt distribution, track BC!=0",
               nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
      fhPtTrigChargedOtherBC->SetYTitle("p_{T h^{#pm}} (GeV/c)");
      fhPtTrigChargedOtherBC->SetXTitle("p_{T trigger} (GeV/c)");
      
      fhXEChargedOtherBC  =
      new TH2F("hXEChargedOtherBC","x_{E} for charged tracks, track BC!=0",
               nptbins,ptmin,ptmax,200,0.,2.);
      fhXEChargedOtherBC->SetYTitle("x_{E}");
      fhXEChargedOtherBC->SetXTitle("p_{T trigger} (GeV/c)");
      
      fhXEUeChargedOtherBC  =
      new TH2F("hXEUeChargedOtherBC","x_{E} for Underlying Event, track BC!=0",
               nptbins,ptmin,ptmax,200,0.,2.);
      fhXEUeChargedOtherBC->SetYTitle("x_{E}");
      fhXEUeChargedOtherBC->SetXTitle("p_{T trigger} (GeV/c)");
      
      fhZTChargedOtherBC  =
      new TH2F("hZTChargedOtherBC","z_{T} for charged tracks, track BC!=0",
               nptbins,ptmin,ptmax,200,0.,2.);
      fhZTChargedOtherBC->SetYTitle("z_{T}");
      fhZTChargedOtherBC->SetXTitle("p_{T trigger}");
      
      fhZTUeChargedOtherBC  =
      new TH2F("hZTUeChargedOtherBC","z_{T} for Underlying Event, track BC!=0",
               nptbins,ptmin,ptmax,200,0.,2.);
      fhZTUeChargedOtherBC->SetYTitle("z_{T}");
      fhZTUeChargedOtherBC->SetXTitle("p_{T trigger} (GeV/c)");
      
      outputContainer->Add(fhDeltaPhiChargedOtherBC) ;
      outputContainer->Add(fhDeltaPhiChargedPtA3GeVOtherBC) ;
      outputContainer->Add(fhXEChargedOtherBC) ;
      outputContainer->Add(fhXEUeChargedOtherBC) ;
      outputContainer->Add(fhZTChargedOtherBC) ;
      outputContainer->Add(fhZTUeChargedOtherBC) ;
      outputContainer->Add(fhPtTrigChargedOtherBC) ;    

      fhDeltaPhiChargedBC0  = new TH2F
      ("hDeltaPhiChargedBC0","#phi_{trigger} - #phi_{h^{#pm}} vs p_{T trigger}, track BC==0",
       nptbins,ptmin,ptmax, ndeltaphibins ,deltaphimin,deltaphimax);
      fhDeltaPhiChargedBC0->SetYTitle("#Delta #phi (rad)");
      fhDeltaPhiChargedBC0->SetXTitle("p_{T trigger} (GeV/c)");
      
      fhDeltaPhiChargedPtA3GeVBC0  = new TH2F
      ("hDeltaPhiChargedPtA3GeVBC0","#phi_{trigger} - #phi_{h^{#pm}} vs p_{T trigger}, p_{TA}>3 GeV/c, track BC==0",
       nptbins,ptmin,ptmax, ndeltaphibins ,deltaphimin,deltaphimax);
      fhDeltaPhiChargedPtA3GeVBC0->SetYTitle("#Delta #phi (rad)");
      fhDeltaPhiChargedPtA3GeVBC0->SetXTitle("p_{T trigger} (GeV/c)");
      
      fhPtTrigChargedBC0  =
      new TH2F("hPtTrigChargedBC0","trigger and charged tracks pt distribution, track BC==0",
               nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
      fhPtTrigChargedBC0->SetYTitle("p_{T h^{#pm}} (GeV/c)");
      fhPtTrigChargedBC0->SetXTitle("p_{T trigger} (GeV/c)");
      
      fhXEChargedBC0  =
      new TH2F("hXEChargedBC0","x_{E} for charged tracks, track BC==0",
               nptbins,ptmin,ptmax,200,0.,2.);
      fhXEChargedBC0->SetYTitle("x_{E}");
      fhXEChargedBC0->SetXTitle("p_{T trigger} (GeV/c)");
      
      fhXEUeChargedBC0  =
      new TH2F("hXEUeChargedBC0","x_{E} for Underlying Event, track BC==0",
               nptbins,ptmin,ptmax,200,0.,2.);
      fhXEUeChargedBC0->SetYTitle("x_{E}");
      fhXEUeChargedBC0->SetXTitle("p_{T trigger} (GeV/c)");
      
      fhZTChargedBC0  =
      new TH2F("hZTChargedBC0","z_{T} for charged tracks, track BC==0",
               nptbins,ptmin,ptmax,200,0.,2.);
      fhZTChargedBC0->SetYTitle("z_{T}");
      fhZTChargedBC0->SetXTitle("p_{T trigger}");
      
      fhZTUeChargedBC0  =
      new TH2F("hZTUeChargedBC0","z_{T} for Underlying Event, track BC==0",
               nptbins,ptmin,ptmax,200,0.,2.);
      fhZTUeChargedBC0->SetYTitle("z_{T}");
      fhZTUeChargedBC0->SetXTitle("p_{T trigger} (GeV/c)");
      
      outputContainer->Add(fhDeltaPhiChargedBC0) ;
      outputContainer->Add(fhDeltaPhiChargedPtA3GeVBC0) ;
      outputContainer->Add(fhXEChargedBC0) ;
      outputContainer->Add(fhXEUeChargedBC0) ;
      outputContainer->Add(fhZTChargedBC0) ;
      outputContainer->Add(fhZTUeChargedBC0) ;
      outputContainer->Add(fhPtTrigChargedBC0) ;

      fhPtLeadingVtxBC0  = new TH1F("hPtLeadingVtxBC0","p_{T} distribution of leading particles", nptbins,ptmin,ptmax);
      fhPtLeadingVtxBC0->SetXTitle("p_{T}^{trig} (GeV/c)");
      
      fhDeltaPhiChargedVtxBC0  = new TH2F
      ("hDeltaPhiChargedVtxBC0","#phi_{trigger} - #phi_{h^{#pm}} vs p_{T trigger}, track BC==0",
       nptbins,ptmin,ptmax, ndeltaphibins ,deltaphimin,deltaphimax);
      fhDeltaPhiChargedVtxBC0->SetYTitle("#Delta #phi (rad)");
      fhDeltaPhiChargedVtxBC0->SetXTitle("p_{T trigger} (GeV/c)");
      
      fhDeltaPhiChargedPtA3GeVVtxBC0  = new TH2F
      ("hDeltaPhiChargedPtA3GeVVtxBC0","#phi_{trigger} - #phi_{h^{#pm}} vs p_{T trigger}, p_{TA}>3 GeV/c, track BC==0",
       nptbins,ptmin,ptmax, ndeltaphibins ,deltaphimin,deltaphimax);
      fhDeltaPhiChargedPtA3GeVVtxBC0->SetYTitle("#Delta #phi (rad)");
      fhDeltaPhiChargedPtA3GeVVtxBC0->SetXTitle("p_{T trigger} (GeV/c)");
      
      fhPtTrigChargedVtxBC0  =
      new TH2F("hPtTrigChargedVtxBC0","trigger and charged tracks pt distribution, track BC==0",
               nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
      fhPtTrigChargedVtxBC0->SetYTitle("p_{T h^{#pm}} (GeV/c)");
      fhPtTrigChargedVtxBC0->SetXTitle("p_{T trigger} (GeV/c)");
      
      fhXEChargedVtxBC0  =
      new TH2F("hXEChargedVtxBC0","x_{E} for charged tracks, track BC==0",
               nptbins,ptmin,ptmax,200,0.,2.);
      fhXEChargedVtxBC0->SetYTitle("x_{E}");
      fhXEChargedVtxBC0->SetXTitle("p_{T trigger} (GeV/c)");
      
      fhXEUeChargedVtxBC0  =
      new TH2F("hXEUeChargedVtxBC0","x_{E} for Underlying Event, track BC==0",
               nptbins,ptmin,ptmax,200,0.,2.);
      fhXEUeChargedVtxBC0->SetYTitle("x_{E}");
      fhXEUeChargedVtxBC0->SetXTitle("p_{T trigger} (GeV/c)");
      
      fhZTChargedVtxBC0  =
      new TH2F("hZTChargedVtxBC0","z_{T} for charged tracks, track BC==0",
               nptbins,ptmin,ptmax,200,0.,2.);
      fhZTChargedVtxBC0->SetYTitle("z_{T}");
      fhZTChargedVtxBC0->SetXTitle("p_{T trigger}");
      
      fhZTUeChargedVtxBC0  =
      new TH2F("hZTUeChargedVtxBC0","z_{T} for Underlying Event, track BC==0",
               nptbins,ptmin,ptmax,200,0.,2.);
      fhZTUeChargedVtxBC0->SetYTitle("z_{T}");
      fhZTUeChargedVtxBC0->SetXTitle("p_{T trigger} (GeV/c)");
      
      outputContainer->Add(fhPtLeadingVtxBC0);
      outputContainer->Add(fhDeltaPhiChargedVtxBC0) ;
      outputContainer->Add(fhDeltaPhiChargedPtA3GeVVtxBC0) ;
      outputContainer->Add(fhXEChargedVtxBC0) ;
      outputContainer->Add(fhXEUeChargedVtxBC0) ;
      outputContainer->Add(fhZTChargedVtxBC0) ;
      outputContainer->Add(fhZTUeChargedVtxBC0) ;
      outputContainer->Add(fhPtTrigChargedVtxBC0) ;
      
      for(Int_t i = 0 ; i < 7 ; i++)
      {
        fhPtLeadingPileUp[i]  = new TH1F(Form("hPtLeadingPileUp%s",pileUpName[i].Data()),
                                         Form("p_{T} distribution of leading particles, %s Pile-Up event",pileUpName[i].Data()), nptbins,ptmin,ptmax);
        fhPtLeadingPileUp[i]->SetXTitle("p_{T}^{trig} (GeV/c)");
        outputContainer->Add(fhPtLeadingPileUp[i]);
        
        fhDeltaPhiChargedPileUp[i]  = new TH2F(Form("hDeltaPhiChargedPileUp%s",pileUpName[i].Data()),
                                                    Form("#phi_{trigger} - #phi_{h^{#pm}} vs p_{T trigger}, %s Pile-Up event",pileUpName[i].Data()),
         nptbins,ptmin,ptmax, ndeltaphibins ,deltaphimin,deltaphimax);
        fhDeltaPhiChargedPileUp[i]->SetYTitle("#Delta #phi (rad)");
        fhDeltaPhiChargedPileUp[i]->SetXTitle("p_{T trigger} (GeV/c)");
        outputContainer->Add(fhDeltaPhiChargedPileUp[i]) ;
        
        fhDeltaPhiChargedPtA3GeVPileUp[i]  = new TH2F(Form("hDeltaPhiChargedPtA3GeVPileUp%s",pileUpName[i].Data()),
                                                           Form("#phi_{trigger} - #phi_{h^{#pm}} vs p_{T trigger}, p_{TA}>3 GeV/c, %s Pile-Up event",pileUpName[i].Data()),
         nptbins,ptmin,ptmax, ndeltaphibins ,deltaphimin,deltaphimax);
        fhDeltaPhiChargedPtA3GeVPileUp[i]->SetYTitle("#Delta #phi (rad)");
        fhDeltaPhiChargedPtA3GeVPileUp[i]->SetXTitle("p_{T trigger} (GeV/c)");
        outputContainer->Add(fhDeltaPhiChargedPtA3GeVPileUp[i]) ;
        
        fhDeltaEtaChargedPileUp[i]  = new TH2F(Form("hDeltaEtaChargedPileUp%s",pileUpName[i].Data()),
                                                    Form("#eta_{trigger} - #eta_{h^{#pm}} vs p_{T trigger}, %s Pile-Up event",pileUpName[i].Data()),
         nptbins,ptmin,ptmax,ndeltaetabins,deltaetamin,deltaetamax);  
        fhDeltaEtaChargedPileUp[i]->SetYTitle("#Delta #eta");
        fhDeltaEtaChargedPileUp[i]->SetXTitle("p_{T trigger} (GeV/c)");
        outputContainer->Add(fhDeltaEtaChargedPileUp[i]) ;
        
        fhDeltaEtaChargedPtA3GeVPileUp[i]  = new TH2F(Form("hDeltaEtaChargedPtA3GeVPileUp%s",pileUpName[i].Data()),
                                                           Form("#eta_{trigger} - #eta_{h^{#pm}} vs p_{T trigger}, p_{TA}>3 GeV/c, %s Pile-Up event",pileUpName[i].Data()),
         nptbins,ptmin,ptmax,ndeltaetabins,deltaetamin,deltaetamax);  
        fhDeltaEtaChargedPtA3GeVPileUp[i]->SetYTitle("#Delta #eta");
        fhDeltaEtaChargedPtA3GeVPileUp[i]->SetXTitle("p_{T trigger} (GeV/c)");    
        outputContainer->Add(fhDeltaEtaChargedPtA3GeVPileUp[i]) ;
        
        fhXEChargedPileUp[i]  = new TH2F(Form("hXEChargedPileUp%s",pileUpName[i].Data()),
                                              Form("x_{E} for charged tracks, %s Pile-Up event",pileUpName[i].Data()),
                 nptbins,ptmin,ptmax,200,0.,2.); 
        fhXEChargedPileUp[i]->SetYTitle("x_{E}");
        fhXEChargedPileUp[i]->SetXTitle("p_{T trigger} (GeV/c)");
        outputContainer->Add(fhXEChargedPileUp[i]) ;
        
        fhXEUeChargedPileUp[i]  = new TH2F(Form("hXEUeChargedPileUp%s",pileUpName[i].Data()),
                                                Form("x_{E} for Underlying Event, %s Pile-Up event",pileUpName[i].Data()),
                 nptbins,ptmin,ptmax,200,0.,2.); 
        fhXEUeChargedPileUp[i]->SetYTitle("x_{E}");
        fhXEUeChargedPileUp[i]->SetXTitle("p_{T trigger} (GeV/c)");
        outputContainer->Add(fhXEUeChargedPileUp[i]) ;
        
        fhZTChargedPileUp[i]  = new TH2F(Form("hZTChargedPileUp%s",pileUpName[i].Data()),
                                              Form("z_{T} for charged tracks, %s Pile-Up event",pileUpName[i].Data()),
                 nptbins,ptmin,ptmax,200,0.,2.); 
        fhZTChargedPileUp[i]->SetYTitle("z_{T}");
        fhZTChargedPileUp[i]->SetXTitle("p_{T trigger} (GeV/c)");
        outputContainer->Add(fhZTChargedPileUp[i]) ;
        
        fhZTUeChargedPileUp[i]  = new TH2F(Form("hZTUeChargedPileUp%s",pileUpName[i].Data()),
                                                Form("z_{T} for Underlying Event, %s Pile-Up event",pileUpName[i].Data()),
                 nptbins,ptmin,ptmax,200,0.,2.); 
        fhZTUeChargedPileUp[i]->SetYTitle("z_{T}");
        fhZTUeChargedPileUp[i]->SetXTitle("p_{T trigger} (GeV/c)");
        outputContainer->Add(fhZTUeChargedPileUp[i]) ;
        
        fhPtTrigChargedPileUp[i]  = new TH2F(Form("hPtTrigChargedPileUp%s",pileUpName[i].Data()),
                                                  Form("trigger and charged tracks pt distribution, %s Pile-Up event",pileUpName[i].Data()),
                 nptbins,ptmin,ptmax,nptbins,ptmin,ptmax); 
        fhPtTrigChargedPileUp[i]->SetYTitle("p_{T h^{#pm}} (GeV/c)");
        fhPtTrigChargedPileUp[i]->SetXTitle("p_{T trigger} (GeV/c)");    
        outputContainer->Add(fhPtTrigChargedPileUp[i]) ;
        
      }
    }
    
    if(DoEventSelect())
    { 
      Int_t nMultiBins = GetMultiBin();
      fhTrigDeltaPhiCharged = new TH2F*[nMultiBins] ;
      fhTrigDeltaEtaCharged = new TH2F*[nMultiBins] ;
      fhTrigXECorr          = new TH2F*[nMultiBins] ;
      fhTrigXEUeCorr        = new TH2F*[nMultiBins] ;
      fhTrigZTCorr          = new TH2F*[nMultiBins] ;
      fhTrigZTUeCorr        = new TH2F*[nMultiBins] ;
      
      for(Int_t im=0; im<nMultiBins; im++)
      {
        fhTrigDeltaPhiCharged[im]  = new TH2F 
        (Form("hTrigDeltaPhiCharged_%d",im),Form("hTrigDeltaPhiCharged_%d",im), nptbins,ptmin,ptmax, ndeltaphibins ,deltaphimin,deltaphimax); 
        fhTrigDeltaPhiCharged[im]->SetXTitle("p_{T trigger} (GeV/c)");
        fhTrigDeltaPhiCharged[im]->SetYTitle("#Delta #phi (rad)");
        
        fhTrigDeltaEtaCharged[im]  = new TH2F 
        (Form("hTrigDeltaEtaCharged_%d",im),Form("hTrigDeltaEtaCharged_%d",im), nptbins,ptmin,ptmax, ndeltaetabins ,deltaetamin,deltaetamax); 
        fhTrigDeltaEtaCharged[im]->SetXTitle("p_{T trigger} (GeV/c)");
        fhTrigDeltaEtaCharged[im]->SetYTitle("#Delta #eta");
        
        fhTrigXECorr[im]  = new TH2F
        (Form("hTrigXEPtCorr_%d",im),Form("hTrigXEPtCorr_%d",im), nptbins,ptmin,ptmax,200,0.,2.); 
        fhTrigXECorr[im]->SetYTitle("x_{E trigger h^{#pm}}");
        fhTrigXECorr[im]->SetXTitle("p_{T trigger} (GeV/c)");
        
        fhTrigXEUeCorr[im]  = new TH2F
        (Form("hTrigXEPtUeCorr_%d",im),Form("hTrigXEPtUeCorr_%d",im), nptbins,ptmin,ptmax,200,0.,2.); 
        fhTrigXEUeCorr[im]->SetYTitle("x_{E trigger h^{#pm}}");
        fhTrigXEUeCorr[im]->SetXTitle("p_{T trigger}(GeV/c)");       
        
        fhTrigZTCorr[im]  = new TH2F
        (Form("hTrigZTPtCorr_%d",im),Form("hTrigZTPtCorr_%d",im), nptbins,ptmin,ptmax,200,0.,2.); 
        fhTrigZTCorr[im]->SetYTitle("z_{trigger h^{#pm}}");
        fhTrigZTCorr[im]->SetXTitle("p_{T trigger} (GeV/c)");
        
        fhTrigZTUeCorr[im]  = new TH2F
        (Form("hTrigZTPtUeCorr_%d",im),Form("hTrigZTPtUeCorr_%d",im), nptbins,ptmin,ptmax,200,0.,2.); 
        fhTrigZTUeCorr[im]->SetYTitle("z_{trigger h^{#pm}}");
        fhTrigZTUeCorr[im]->SetXTitle("p_{T trigger} (GeV/c)");               
        
        outputContainer->Add(fhTrigDeltaPhiCharged[im]) ;
        outputContainer->Add(fhTrigDeltaEtaCharged[im]) ;
        outputContainer->Add(fhTrigXECorr[im]);
        outputContainer->Add(fhTrigXEUeCorr[im]);
        outputContainer->Add(fhTrigZTCorr[im]);
        outputContainer->Add(fhTrigZTUeCorr[im]);
      }
    }
    
    if(fFillBradHisto)
    {
      fhAssocPtBkg        = new TH2F("hAssocPtBkg", " Trigger p_{T} vs associated hadron p_{T} from background",
                                     nptbins, ptmin, ptmax,nptbins,ptmin,ptmax);
      fhAssocPtBkg->SetXTitle("p_{T trigger} (GeV/c)");
      fhAssocPtBkg->SetYTitle("p_{T associated} (GeV/c)");
      outputContainer->Add(fhAssocPtBkg) ;
      
      fhDeltaPhiBrad = new TH2F("hDeltaPhiBrad","atan2(sin(#Delta #phi), cos(#Delta #phi))/#pi vs p_{T trigger} ", 
                                nptbins, ptmin, ptmax,288, -1.0/3.0, 5.0/3.0);
      fhDeltaPhiBrad->SetXTitle("p_{T trigger} (GeV/c)");
      fhDeltaPhiBrad->SetYTitle("atan2(sin(#Delta #phi), cos(#Delta #phi))/#pi");
      outputContainer->Add(fhDeltaPhiBrad) ;
    }

    fhDeltaPhiDeltaEtaAssocPtBin = new TH2F*[fNAssocPtBins*nz];
    fhDeltaPhiAssocPtBin       = new TH2F*[fNAssocPtBins*nz];
    fhDeltaPhiAssocPtBinDEta08 = new TH2F*[fNAssocPtBins*nz];
    fhDeltaPhiAssocPtBinDEta0  = new TH2F*[fNAssocPtBins*nz];
    fhXEAssocPtBin             = new TH2F*[fNAssocPtBins*nz];
    fhZTAssocPtBin             = new TH2F*[fNAssocPtBins*nz];
    
    if(fFillBradHisto)  
      fhDeltaPhiBradAssocPtBin = new TH2F*[fNAssocPtBins*nz];
    
    if(fPi0Trigger || fDecayTrigger)
    {
      fhDeltaPhiAssocPtBin       = new TH2F*[fNAssocPtBins*nz];
      fhDeltaPhiAssocPtBinDEta08 = new TH2F*[fNAssocPtBins*nz];
      fhXEAssocPtBin             = new TH2F*[fNAssocPtBins*nz];
      fhZTAssocPtBin             = new TH2F*[fNAssocPtBins*nz];
      fhXEDecayChargedAssocPtBin = new TH2F*[fNAssocPtBins*nz];
      fhZTDecayChargedAssocPtBin = new TH2F*[fNAssocPtBins*nz];
      fhDeltaPhiDecayChargedAssocPtBin = new TH2F*[fNAssocPtBins*nz];
    }

    if(fHMPIDCorrelation)
    {
      fhDeltaPhiAssocPtBinHMPID   = new TH2F*[fNAssocPtBins*nz];
      fhDeltaPhiAssocPtBinHMPIDAcc= new TH2F*[fNAssocPtBins*nz];
    }
    
    for(Int_t i = 0 ; i < fNAssocPtBins ; i++)
    {
      for(Int_t z = 0 ; z < nz ; z++)
      {
        Int_t bin = i*nz+z;
        
        if(fCorrelVzBin)
        {
          sz = "_vz%d"+z;
          tz = ", v_{z} bin "+z;
        }
        
        //printf("iAssoc %d, Vz %d, bin %d - sz %s, tz %s	\n",i,z,bin,sz.Data(),tz.Data());
        
        fhDeltaPhiDeltaEtaAssocPtBin[bin]  = new TH2F(Form("hDeltaPhiDeltaEtaPtAssocPt%2.1f_%2.1f%s", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1],sz.Data()), 
                                                 Form("#Delta #phi vs #Delta #eta for associated p_{T} bin [%2.1f,%2.1f]%s", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1],tz.Data()), 
                                                 ndeltaphibins ,deltaphimin,deltaphimax,ndeltaetabins,deltaetamin,deltaetamax); 
        fhDeltaPhiDeltaEtaAssocPtBin[bin]->SetXTitle("#Delta #phi (rad)");
        fhDeltaPhiDeltaEtaAssocPtBin[bin]->SetYTitle("#Delta #eta");  
        
        fhDeltaPhiAssocPtBin[bin] = new TH2F(Form("hDeltaPhiPtAssocPt%2.1f_%2.1f%s", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1],sz.Data()), 
                                           Form("#Delta #phi vs p_{T trigger} for associated p_{T} bin [%2.1f,%2.1f]%s", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1],tz.Data()), 
                                           nptbins, ptmin, ptmax, ndeltaphibins ,deltaphimin,deltaphimax);
        fhDeltaPhiAssocPtBin[bin]->SetXTitle("p_{T trigger} (GeV/c)");
        fhDeltaPhiAssocPtBin[bin]->SetYTitle("#Delta #phi (rad)");
        
        fhDeltaPhiAssocPtBinDEta08[bin] = new TH2F(Form("hDeltaPhiDeltaEta0.8PtAssocPt%2.1f_%2.1f%s", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1],sz.Data()), 
                                                 Form("#Delta #phi vs p_{T trigger} for associated p_{T} bin [%2.1f,%2.1f]%s, for #Delta #eta > 0.8", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1],tz.Data()), 
                                                 nptbins, ptmin, ptmax, ndeltaphibins ,deltaphimin,deltaphimax);
        fhDeltaPhiAssocPtBinDEta08[bin]->SetXTitle("p_{T trigger} (GeV/c)");
        fhDeltaPhiAssocPtBinDEta08[bin]->SetYTitle("#Delta #phi (rad)");      

        fhDeltaPhiAssocPtBinDEta0[bin] = new TH2F(Form("hDeltaPhiDeltaEta0PtAssocPt%2.1f_%2.1f%s", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1],tz.Data()), 
                                                   Form("#Delta #phi vs p_{T trigger} for associated p_{T} bin [%2.1f,%2.1f]%s, for #Delta #eta = 0.", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1],tz.Data()), 
                                                   nptbins, ptmin, ptmax, ndeltaphibins ,deltaphimin,deltaphimax);
        fhDeltaPhiAssocPtBinDEta0[bin]->SetXTitle("p_{T trigger} (GeV/c)");
        fhDeltaPhiAssocPtBinDEta0[bin]->SetYTitle("#Delta #phi (rad)");    
        
        fhXEAssocPtBin[bin]       = new TH2F(Form("hXEAssocPtBin%1.f_%1.f%s", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1],sz.Data()), 
                                           Form("x_{E} vs p_{T trigger} for associated p_{T} bin [%2.1f,%2.1f]%s", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1],tz.Data()), 
                                           nptbins, ptmin, ptmax,200, 0.0, 2.0);
        fhXEAssocPtBin[bin]->SetXTitle("p_{T trigger} (GeV/c)");
        fhXEAssocPtBin[bin]->SetYTitle("x_{E}");
        
        fhZTAssocPtBin[bin]       = new TH2F(Form("hZTAssocPtBin%1.f_%1.f%s", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1],sz.Data()), 
                                           Form("z_{T} vs p_{T trigger} for associated p_{T} bin [%2.1f,%2.1f]%s", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1],tz.Data()), 
                                           nptbins, ptmin, ptmax,200, 0.0, 2.0);
        fhZTAssocPtBin[bin]->SetXTitle("p_{T trigger} (GeV/c)");
        fhZTAssocPtBin[bin]->SetYTitle("z_{T}");
        
        outputContainer->Add(fhDeltaPhiDeltaEtaAssocPtBin[bin]) ;
        outputContainer->Add(fhDeltaPhiAssocPtBin[bin]) ;
        outputContainer->Add(fhDeltaPhiAssocPtBinDEta08[bin]) ;
        outputContainer->Add(fhDeltaPhiAssocPtBinDEta0[bin]) ;
        outputContainer->Add(fhXEAssocPtBin[bin]);
        outputContainer->Add(fhZTAssocPtBin[bin]);

        if(fPi0Trigger || fDecayTrigger) 
        {
          fhDeltaPhiDecayChargedAssocPtBin[bin] = new TH2F(Form("hDeltaPhiPtDecayChargedAssocPt%2.1f_%2.1f%s", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1],sz.Data()), 
                                                         Form("#Delta #phi vs p_{T trigger} tagged as decay for associated p_{T} bin [%2.1f,%2.1f]%s", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1],tz.Data()), 
                                                         nptbins, ptmin, ptmax, ndeltaphibins ,deltaphimin,deltaphimax);
          fhDeltaPhiDecayChargedAssocPtBin[bin]->SetXTitle("p_{T trigger} (GeV/c)");
          fhDeltaPhiDecayChargedAssocPtBin[bin]->SetYTitle("#Delta #phi (rad)");
          
          fhXEDecayChargedAssocPtBin[bin]       = new TH2F(Form("hXEDecayChargedAssocPtBin%1.f_%1.f%s", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1],sz.Data()), 
                                                         Form("x_{E} vs p_{T trigger} tagged as decay for associated p_{T} bin [%2.1f,%2.1f]%s", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1],tz.Data()), 
                                                         nptbins, ptmin, ptmax,200, 0.0, 2.0);
          fhXEDecayChargedAssocPtBin[bin]->SetXTitle("p_{T trigger} (GeV/c)");
          fhXEDecayChargedAssocPtBin[bin]->SetYTitle("x_{E}");
          
          fhZTDecayChargedAssocPtBin[bin]       = new TH2F(Form("hZTDecayChargedAssocPtBin%1.f_%1.f%s", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1],sz.Data()), 
                                                         Form("z_{T} vs p_{T trigger} tagged as decay for associated p_{T} bin [%2.1f,%2.1f]%s", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1],tz.Data()), 
                                                         nptbins, ptmin, ptmax,200, 0.0, 2.0);
          fhZTDecayChargedAssocPtBin[bin]->SetXTitle("p_{T trigger} (GeV/c)");
          fhZTDecayChargedAssocPtBin[bin]->SetYTitle("z_{T}");
          
          outputContainer->Add(fhDeltaPhiDecayChargedAssocPtBin[bin]) ;
          outputContainer->Add(fhXEDecayChargedAssocPtBin[bin]);
          outputContainer->Add(fhZTDecayChargedAssocPtBin[bin]);
          
        }
        
        if(fFillBradHisto) 
        {
          fhDeltaPhiBradAssocPtBin[bin] = new TH2F(Form("hDeltaPhiBradPtAssocPt%2.1f_%2.1f%s", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1],sz.Data()), 
                                                 Form("atan2(sin(#Delta #phi), cos(#Delta #phi))/#pi vs p_{T trigger} for associated p_{T} bin [%2.1f,%2.1f]%s", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1],tz.Data()), 
                                                 nptbins, ptmin, ptmax,288, -1.0/3.0, 5.0/3.0);
          fhDeltaPhiBradAssocPtBin[bin]->SetXTitle("p_{T trigger} (GeV/c)");
          fhDeltaPhiBradAssocPtBin[bin]->SetYTitle("atan2(sin(#Delta #phi), cos(#Delta #phi))/#pi");
          outputContainer->Add(fhDeltaPhiBradAssocPtBin[bin]) ;
        }       
        
        if(fHMPIDCorrelation)
        {
          fhDeltaPhiAssocPtBinHMPID[bin] = new TH2F(Form("hDeltaPhiPtAssocPt%2.1f_%2.1f%sHMPID", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1],sz.Data()), 
                                                  Form("#Delta #phi vs p_{T trigger} for associated p_{T} bin [%2.1f,%2.1f]%s, with track having HMPID signal", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1],tz.Data()), 
                                                  nptbins, ptmin, ptmax, ndeltaphibins ,deltaphimin,deltaphimax);
          fhDeltaPhiAssocPtBinHMPID[bin]->SetXTitle("p_{T trigger} (GeV/c)" );
          fhDeltaPhiAssocPtBinHMPID[bin]->SetYTitle("#Delta #phi (rad)");      
          
          fhDeltaPhiAssocPtBinHMPIDAcc[bin] = new TH2F(Form("hDeltaPhiPtAssocPt%2.1f_%2.1f%sHMPIDAcc", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1],sz.Data()), 
                                                     Form("#Delta #phi vs p_{T trigger} for associated p_{T} bin [%2.1f,%2.1f]%s, with track within 5<phi<20 deg", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1],tz.Data()), 
                                                     nptbins, ptmin, ptmax, ndeltaphibins ,deltaphimin,deltaphimax);
          fhDeltaPhiAssocPtBinHMPIDAcc[bin]->SetXTitle("p_{T trigger} (GeV/c)");
          fhDeltaPhiAssocPtBinHMPIDAcc[bin]->SetYTitle("#Delta #phi (rad)"); 
          
          outputContainer->Add(fhDeltaPhiAssocPtBinHMPID[bin]) ;
          outputContainer->Add(fhDeltaPhiAssocPtBinHMPIDAcc[bin]) ;
          
        }      
      }
    }
    
    if(fPi0Trigger || fDecayTrigger)
    {
      if(fPi0Trigger)
      {
        fhPtPi0DecayRatio  = new TH2F
        ("hPtPi0DecayRatio","p_{T} of #pi^{0} and the ratio of pt for two decay", 
         nptbins,ptmin,ptmax, 100,0.,2.); 
        fhPtPi0DecayRatio->SetXTitle("p_{T}^{#pi^{0}} (GeV/c)");
        fhPtPi0DecayRatio->SetYTitle("p_{T}^{Decay}/p_{T}^{#pi^{0}}");
        outputContainer->Add(fhPtPi0DecayRatio) ; 
      }
      
      fhDeltaPhiDecayCharged  = new TH2F
      ("hDeltaPhiDecayCharged","#phi_{Decay} - #phi_{h^{#pm}} vs p_{T Decay}",
       nptbins,ptmin,ptmax, ndeltaphibins ,deltaphimin,deltaphimax); 
      fhDeltaPhiDecayCharged->SetYTitle("#Delta #phi (rad)");
      fhDeltaPhiDecayCharged->SetXTitle("p_{T Decay} (GeV/c)");
      
      fhXEDecayCharged  = 
      new TH2F("hXEDecayCharged","x_{E}  Decay",
               nptbins,ptmin,ptmax,200,0.,2.); 
      fhXEDecayCharged->SetYTitle("x_{E}");
      fhXEDecayCharged->SetXTitle("p_{T decay} (GeV/c)");
      
      fhZTDecayCharged  = 
      new TH2F("hZTDecayCharged","z_{trigger h^{#pm}} = p_{T h^{#pm}} / p_{T Decay}",
               nptbins,ptmin,ptmax,200,0.,2.); 
      fhZTDecayCharged->SetYTitle("z_{decay h^{#pm}}");
      fhZTDecayCharged->SetXTitle("p_{T decay} (GeV/c)");      
      
      outputContainer->Add(fhDeltaPhiDecayCharged) ; 
      outputContainer->Add(fhXEDecayCharged) ;
      outputContainer->Add(fhZTDecayCharged) ;
    }    
    
    if(fMakeSeveralUE)
    { 
      fhDeltaPhiUeLeftCharged  = new TH2F
      ("hDeltaPhiUeLeftChargedPt","#phi_{trigger} - #phi_{#Ueh^{#pm}} vs p_{T Ueh^{#pm}} with UE left side range of trigger particles",
       nptbins,ptmin,ptmax, ndeltaphibins ,deltaphimin,deltaphimax);
      fhDeltaPhiUeLeftCharged->SetYTitle("#Delta #phi (rad)");
      fhDeltaPhiUeLeftCharged->SetXTitle("p_{T h^{#pm}} (GeV/c)");
      outputContainer->Add(fhDeltaPhiUeLeftCharged) ;
      
      fhDeltaPhiUeRightCharged  = new TH2F
      ("hDeltaPhiUeRightChargedPt","#phi_{trigger} - #phi_{#Ueh^{#pm}} vs p_{T Ueh^{#pm}} with UE right side range of trigger particles",
       nptbins,ptmin,ptmax, ndeltaphibins ,deltaphimin,deltaphimax);
      fhDeltaPhiUeRightCharged->SetYTitle("#Delta #phi (rad)");
      fhDeltaPhiUeRightCharged->SetXTitle("p_{T h^{#pm}} (GeV/c)");
      outputContainer->Add(fhDeltaPhiUeRightCharged) ;
      
      fhDeltaPhiUeLeftUpCharged  = new TH2F
      ("hDeltaPhiUeLeftUpChargedPt","#phi_{trigger} - #phi_{#Ueh^{#pm}} vs p_{T Ueh^{#pm}} with UE left Up side range of trigger particles",
       nptbins,ptmin,ptmax, ndeltaphibins ,deltaphimin,deltaphimax);
      fhDeltaPhiUeLeftUpCharged->SetYTitle("#Delta #phi (rad)");
      fhDeltaPhiUeLeftUpCharged->SetXTitle("p_{T h^{#pm}} (GeV/c)");
      outputContainer->Add(fhDeltaPhiUeLeftUpCharged) ;
      
      fhDeltaPhiUeRightUpCharged  = new TH2F
      ("hDeltaPhiUeRightUpChargedPt","#phi_{trigger} - #phi_{#Ueh^{#pm}} vs p_{T Ueh^{#pm}} with UE right Up side range of trigger particles",
       nptbins,ptmin,ptmax, ndeltaphibins ,deltaphimin,deltaphimax);
      fhDeltaPhiUeRightUpCharged->SetYTitle("#Delta #phi (rad)");
      fhDeltaPhiUeRightUpCharged->SetXTitle("p_{T h^{#pm}} (GeV/c)");
      outputContainer->Add(fhDeltaPhiUeRightUpCharged) ;
      
      fhDeltaPhiUeLeftDownCharged  = new TH2F
      ("hDeltaPhiUeLeftDownChargedPt","#phi_{trigger} - #phi_{#Ueh^{#pm}} vs p_{T Ueh^{#pm}} with UE left Down side range of trigger particles",
       nptbins,ptmin,ptmax, ndeltaphibins ,deltaphimin,deltaphimax);
      fhDeltaPhiUeLeftDownCharged->SetYTitle("#Delta #phi (rad)");
      fhDeltaPhiUeLeftDownCharged->SetXTitle("p_{T h^{#pm}} (GeV/c)");
      outputContainer->Add(fhDeltaPhiUeLeftDownCharged) ;
      
      fhDeltaPhiUeRightDownCharged  = new TH2F
      ("hDeltaPhiUeRightDownChargedPt","#phi_{trigger} - #phi_{#Ueh^{#pm}} vs p_{T Ueh^{#pm}} with UE right Down side range of trigger particles",
       nptbins,ptmin,ptmax, ndeltaphibins ,deltaphimin,deltaphimax);
      fhDeltaPhiUeRightDownCharged->SetYTitle("#Delta #phi (rad)");
      fhDeltaPhiUeRightDownCharged->SetXTitle("p_{T h^{#pm}} (GeV/c)");
      outputContainer->Add(fhDeltaPhiUeRightDownCharged) ;
      
      fhXEUeLeftCharged  = 
      new TH2F("hXEUeChargedLeft","x_{E} with UE left side of trigger",
               nptbins,ptmin,ptmax,200,0.,2.); 
      fhXEUeLeftCharged->SetYTitle("x_{E Ueh^{#pm}}");
      fhXEUeLeftCharged->SetXTitle("p_{T trigger} (GeV/c)");
      outputContainer->Add(fhXEUeLeftCharged) ;
      
      fhXEUeRightCharged  = 
      new TH2F("hXEUeChargedRight","x_{E h^{#pm}} with UE right side of trigger",
               nptbins,ptmin,ptmax,200,0.,2.); 
      fhXEUeRightCharged->SetYTitle("z_{trigger Ueh^{#pm}}");
      fhXEUeRightCharged->SetXTitle("p_{T trigger} (GeV/c)");
      outputContainer->Add(fhXEUeRightCharged) ;
      
      fhXEUeLeftUpCharged  = 
      new TH2F("hXEUeChargedLeftUp","x_{E} with UE left Up side of trigger",
               nptbins,ptmin,ptmax,200,0.,2.); 
      fhXEUeLeftUpCharged->SetYTitle("x_{E Ueh^{#pm}}");
      fhXEUeLeftUpCharged->SetXTitle("p_{T trigger} (GeV/c)");
      outputContainer->Add(fhXEUeLeftUpCharged) ;
      
      fhXEUeRightUpCharged  = 
      new TH2F("hXEUeChargedRightUp","x_{E h^{#pm}} with UE right Up side of trigger",
               nptbins,ptmin,ptmax,200,0.,2.); 
      fhXEUeRightUpCharged->SetYTitle("z_{trigger Ueh^{#pm}}");
      fhXEUeRightUpCharged->SetXTitle("p_{T trigger} (GeV/c)");
      outputContainer->Add(fhXEUeRightUpCharged) ;
      
      fhXEUeLeftDownCharged  = 
      new TH2F("hXEUeChargedLeftDown","x_{E} with UE left Down side of trigger",
               nptbins,ptmin,ptmax,200,0.,2.); 
      fhXEUeLeftDownCharged->SetYTitle("x_{E Ueh^{#pm}}");
      fhXEUeLeftDownCharged->SetXTitle("p_{T trigger} (GeV/c)");
      outputContainer->Add(fhXEUeLeftDownCharged) ;
      
      fhXEUeRightDownCharged  = 
      new TH2F("hXEUeChargedRightDown","x_{E h^{#pm}} with UE right Down side of trigger",
               nptbins,ptmin,ptmax,200,0.,2.); 
      fhXEUeRightDownCharged->SetYTitle("z_{trigger Ueh^{#pm}}");
      fhXEUeRightDownCharged->SetXTitle("p_{T trigger} (GeV/c)");
      outputContainer->Add(fhXEUeRightDownCharged) ;
      
      fhPtHbpXEUeLeftCharged  = 
      new TH2F("hHbpXEUeChargedLeft","#xi = ln(1/x_{E}) with charged UE left side of trigger",
               nptbins,ptmin,ptmax,200,0.,10.); 
      fhPtHbpXEUeLeftCharged->SetYTitle("ln(1/x_{E})");
      fhPtHbpXEUeLeftCharged->SetXTitle("p_{T trigger} (GeV/c)");
      outputContainer->Add(fhPtHbpXEUeLeftCharged) ;
      
      fhPtHbpXEUeRightCharged  = 
      new TH2F("hHbpXEUeChargedRight","#xi = ln(1/x_{E}) with charged UE right side of trigger",
               nptbins,ptmin,ptmax,200,0.,10.); 
      fhPtHbpXEUeRightCharged->SetYTitle("ln(1/x_{E})");
      fhPtHbpXEUeRightCharged->SetXTitle("p_{T trigger} (GeV/c)");
      outputContainer->Add(fhPtHbpXEUeRightCharged) ;
      
      fhZTUeLeftCharged  = 
      new TH2F("hZTUeChargedLeft","z_{trigger h^{#pm}} = p_{T Ueh^{#pm}} / p_{T trigger} with UE left side of trigger",
               nptbins,ptmin,ptmax,200,0.,2.); 
      fhZTUeLeftCharged->SetYTitle("z_{trigger Ueh^{#pm}}");
      fhZTUeLeftCharged->SetXTitle("p_{T trigger} (GeV/c)");
      outputContainer->Add(fhZTUeLeftCharged) ;
      
      fhZTUeRightCharged  = 
      new TH2F("hZTUeChargedRight","z_{trigger h^{#pm}} = p_{T Ueh^{#pm}} / p_{T trigger} with UE right side of trigger",
               nptbins,ptmin,ptmax,200,0.,2.); 
      fhZTUeRightCharged->SetYTitle("z_{trigger Ueh^{#pm}}");
      fhZTUeRightCharged->SetXTitle("p_{T trigger} (GeV/c)");
      outputContainer->Add(fhZTUeRightCharged) ;      
      
      fhPtHbpZTUeLeftCharged  = 
      new TH2F("hHbpZTUeChargedLeft","#xi = ln(1/z_{T}) with charged UE left side of trigger",
               nptbins,ptmin,ptmax,200,0.,10.); 
      fhPtHbpZTUeLeftCharged->SetYTitle("ln(1/z_{T})");
      fhPtHbpZTUeLeftCharged->SetXTitle("p_{T trigger} (GeV/c)");
      outputContainer->Add(fhPtHbpZTUeLeftCharged) ;
      
      fhPtHbpZTUeRightCharged  = 
      new TH2F("hHbpZTUeChargedRight","#xi = ln(1/z_{T}) with charged UE right side of trigger",
               nptbins,ptmin,ptmax,200,0.,10.); 
      fhPtHbpZTUeRightCharged->SetYTitle("ln(1/z_{T})");
      fhPtHbpZTUeRightCharged->SetXTitle("p_{T trigger} (GeV/c)");
      outputContainer->Add(fhPtHbpZTUeRightCharged) ;
      
    } 
  }  //Correlation with charged hadrons

  //Correlation with neutral hadrons
  if(fNeutralCorr)
  {
    fhDeltaPhiDeltaEtaNeutral  = new TH2F
    ("hDeltaPhiDeltaEtaNeutral","#phi_{trigger} - #phi_{h^{0}} vs #eta_{trigger} - #eta_{h^{0}}",
     ndeltaphibins ,deltaphimin,deltaphimax, ndeltaetabins ,deltaetamin,deltaetamax); 
    fhDeltaPhiDeltaEtaNeutral->SetXTitle("#Delta #phi (rad)");
    fhDeltaPhiDeltaEtaNeutral->SetYTitle("#Delta #eta");   
	  
    fhPhiNeutral  = new TH2F
    ("hPhiNeutral","#phi_{#pi^{0}}  vs p_{T #pi^{0}}",
     nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
    fhPhiNeutral->SetYTitle("#phi_{#pi^{0}} (rad)");
    fhPhiNeutral->SetXTitle("p_{T #pi^{0}} (GeV/c)");
    
    fhEtaNeutral  = new TH2F
    ("hEtaNeutral","#eta_{#pi^{0}}  vs p_{T #pi^{0}}",
     nptbins,ptmin,ptmax,netabins,etamin,etamax); 
    fhEtaNeutral->SetYTitle("#eta_{#pi^{0}} (rad)");
    fhEtaNeutral->SetXTitle("p_{T #pi^{0}} (GeV/c)");
    
    fhDeltaPhiNeutral  = new TH2F
    ("hDeltaPhiNeutral","#phi_{trigger} - #phi_{#pi^{0}} vs p_{T trigger}",
     nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
    fhDeltaPhiNeutral->SetYTitle("#Delta #phi (rad)");
    fhDeltaPhiNeutral->SetXTitle("p_{T trigger} (GeV/c)");
    
    fhDeltaPhiNeutralPt  = new TH2F
    ("hDeltaPhiNeutralPt","#phi_{trigger} - #phi_{#pi^{0}} vs p_{T #pi^{0}}}",
     nptbins,ptmin,ptmax, ndeltaphibins ,deltaphimin,deltaphimax); 
    fhDeltaPhiNeutralPt->SetYTitle("#Delta #phi (rad)");
    fhDeltaPhiNeutralPt->SetXTitle("p_{T h^{0}} (GeV/c)");
    
    fhDeltaPhiUeNeutralPt  = new TH2F
    ("hDeltaPhiUeNeutralPt","#phi_{trigger} - #phi_{#pi^{0}} vs p_{T #pi^{0}}}",
     nptbins,ptmin,ptmax, ndeltaphibins ,deltaphimin,deltaphimax); 
    fhDeltaPhiUeNeutralPt->SetYTitle("#Delta #phi (rad)");
    fhDeltaPhiUeNeutralPt->SetXTitle("p_{T h^{0}} (GeV/c)");
    
    fhDeltaEtaNeutral  = new TH2F
    ("hDeltaEtaNeutral","#eta_{trigger} - #eta_{#pi^{0}} vs p_{T trigger}",
     nptbins,ptmin,ptmax, ndeltaetabins ,deltaetamin,deltaetamax);  
    fhDeltaEtaNeutral->SetYTitle("#Delta #eta");
    fhDeltaEtaNeutral->SetXTitle("p_{T trigger} (GeV/c)");
    
    fhXENeutral  = 
    new TH2F("hXENeutral","x_{E} for #pi^{0} associated",
             nptbins,ptmin,ptmax,200,0.,2.); 
    fhXENeutral->SetYTitle("x_{E}");
    fhXENeutral->SetXTitle("p_{T trigger} (GeV/c)");
    
    fhXEUeNeutral  = 
    new TH2F("hXEUeNeutral","x_{E} for #pi^{0} associated",
             nptbins,ptmin,ptmax,200,0.,2.); 
    fhXEUeNeutral->SetYTitle("x_{E}");
    fhXEUeNeutral->SetXTitle("p_{T trigger} (GeV/c)");
    
    fhPtHbpXENeutral  = 
    new TH2F("hHbpXENeutral","#xi = ln(1/x_{E})for #pi^{0} associated",
             nptbins,ptmin,ptmax,200,0.,10.); 
    fhPtHbpXENeutral->SetYTitle("ln(1/x_{E})");
    fhPtHbpXENeutral->SetXTitle("p_{T trigger} (GeV/c)");
    
    fhPtHbpXEUeNeutral  = 
    new TH2F("hHbpXEUeNeutral","#xi = ln(1/x_{E}) for #pi^{0} associated",
             nptbins,ptmin,ptmax,200,0.,10.); 
    fhPtHbpXEUeNeutral->SetYTitle("ln(1/x_{E})");
    fhPtHbpXEUeNeutral->SetXTitle("p_{T trigger} (GeV/c)");
    
    fhZTNeutral  = 
    new TH2F("hZTNeutral","z_{trigger #pi} = p_{T #pi^{0}} / p_{T trigger} for #pi^{0} associated",
             nptbins,ptmin,ptmax,200,0.,2.); 
    fhZTNeutral->SetYTitle("z_{trigger #pi^{0}}");
    fhZTNeutral->SetXTitle("p_{T trigger} (GeV/c)");
    
    fhZTUeNeutral  = 
    new TH2F("hZTUeNeutral","z_{trigger #pi} = p_{T #pi^{0}} / p_{T trigger} for #pi^{0} associated",
             nptbins,ptmin,ptmax,200,0.,2.); 
    fhZTUeNeutral->SetYTitle("z_{trigger #pi^{0}}");
    fhZTUeNeutral->SetXTitle("p_{T trigger} (GeV/c)");
    
    fhPtHbpZTNeutral  = 
    new TH2F("hHbpZTNeutral","#xi = ln(1/x_{E}) for #pi^{0} associated",
             nptbins,ptmin,ptmax,200,0.,10.); 
    fhPtHbpZTNeutral->SetYTitle("ln(1/z_{T})");
    fhPtHbpZTNeutral->SetXTitle("p_{T trigger}");
    
    fhPtHbpZTUeNeutral  = 
    new TH2F("hHbpZTUeNeutral","#xi = ln(1/x_{E}) for #pi^{0} associated",
             nptbins,ptmin,ptmax,200,0.,10.); 
    fhPtHbpXEUeNeutral->SetYTitle("ln(1/z_{T})");
    fhPtHbpXEUeNeutral->SetXTitle("p_{T trigger} (GeV/c)");    
    
    outputContainer->Add(fhDeltaPhiDeltaEtaNeutral); 
    outputContainer->Add(fhPhiNeutral) ;  
    outputContainer->Add(fhEtaNeutral) ;   
    outputContainer->Add(fhDeltaPhiNeutral) ; 
    outputContainer->Add(fhDeltaPhiNeutralPt) ; 
    outputContainer->Add(fhDeltaPhiUeNeutralPt) ; 
    outputContainer->Add(fhDeltaEtaNeutral) ; 
    outputContainer->Add(fhXENeutral) ;
    outputContainer->Add(fhXEUeNeutral) ;  
    outputContainer->Add(fhPtHbpXENeutral) ;
    outputContainer->Add(fhPtHbpXEUeNeutral) ;    
    outputContainer->Add(fhZTNeutral) ;
    outputContainer->Add(fhZTUeNeutral) ;  
    outputContainer->Add(fhPtHbpZTNeutral) ;
    outputContainer->Add(fhPtHbpZTUeNeutral) ;    
    
    if(fPi0Trigger || fDecayTrigger)
    {
      fhDeltaPhiDecayNeutral  = new TH2F
      ("hDeltaPhiDecayNeutral","#phi_{Decay} - #phi_{h^{0}} vs p_{T Decay}",
       nptbins,ptmin,ptmax, ndeltaphibins ,deltaphimin,deltaphimax);  
      fhDeltaPhiDecayNeutral->SetYTitle("#Delta #phi (rad)");
      fhDeltaPhiDecayNeutral->SetXTitle("p_{T Decay} (GeV/c)");
      
      fhXEDecayNeutral  = 
      new TH2F("hXEDecayNeutral","x_{E} for decay trigger",
               nptbins,ptmin,ptmax,200,0.,2.); 
      fhXEDecayNeutral->SetYTitle("x_{E}");
      fhXEDecayNeutral->SetXTitle("p_{T decay}");
      
      fhZTDecayNeutral  = 
      new TH2F("hZTDecayNeutral","z_{trigger h^{0}} = p_{T h^{0}} / p_{T Decay}",
               nptbins,ptmin,ptmax,200,0.,2.); 
      fhZTDecayNeutral->SetYTitle("z_{h^{0}}");
      fhZTDecayNeutral->SetXTitle("p_{T decay}");      
      
      outputContainer->Add(fhDeltaPhiDecayNeutral) ; 
      outputContainer->Add(fhXEDecayNeutral) ;      
      outputContainer->Add(fhZTDecayNeutral) ;

    }
    
    if(fMakeSeveralUE)
    { 
      fhDeltaPhiUeLeftNeutral  = new TH2F
      ("hDeltaPhiUeLeftNeutralPt","#phi_{trigger} - #phi_{#Ueh^{0}} vs p_{T h^{0}} with neutral UE left side range of trigger particles",
       nptbins,ptmin,ptmax, ndeltaphibins ,deltaphimin,deltaphimax); 
      fhDeltaPhiUeLeftNeutral->SetYTitle("#Delta #phi (rad)");
      fhDeltaPhiUeLeftNeutral->SetXTitle("p_{T h^{0}} (GeV/c)");
      outputContainer->Add(fhDeltaPhiUeLeftNeutral) ;
      
      fhDeltaPhiUeRightNeutral  = new TH2F
      ("hDeltaPhiUeRightNeutralPt","#phi_{trigger} - #phi_{#Ueh^{0}} vs p_{T Ueh^{0}} with neutral UE right side range of trigger particles",
       nptbins,ptmin,ptmax, ndeltaphibins ,deltaphimin,deltaphimax); 
      fhDeltaPhiUeRightNeutral->SetYTitle("#Delta #phi (rad)");
      fhDeltaPhiUeRightNeutral->SetXTitle("p_{T h^{0}} (GeV/c)");
      outputContainer->Add(fhDeltaPhiUeRightNeutral) ;
      
      fhXEUeLeftNeutral  = 
      new TH2F("hXEUeNeutralLeft","x_{E} = p_{T Ueh^{0}} / p_{T trigger} with neutral UE left side of trigger",
               nptbins,ptmin,ptmax,140,0.,2.); 
      fhXEUeLeftNeutral->SetYTitle("z_{trigger Ueh^{0}}");
      fhXEUeLeftNeutral->SetXTitle("p_{T trigger} (GeV/c)");
      outputContainer->Add(fhXEUeLeftNeutral) ;
      
      fhXEUeRightNeutral  = 
      new TH2F("hXEUeNeutralRight","x_{E} = p_{T Ueh^{0}} / p_{T trigger} with neutral UE right side of trigger",
               nptbins,ptmin,ptmax,200,0.,2.); 
      fhXEUeRightNeutral->SetYTitle("z_{trigger Ueh^{0}}");
      fhXEUeRightNeutral->SetXTitle("p_{T trigger} (GeV/c)");
      outputContainer->Add(fhXEUeRightNeutral) ;
      
      fhPtHbpXEUeLeftNeutral  = 
      new TH2F("hHbpXEUeNeutralLeft","#xi = ln(1/x_{E}) with neutral UE left side of trigger",
               nptbins,ptmin,ptmax,200,0.,10.); 
      fhPtHbpXEUeLeftNeutral->SetYTitle("ln(1/x_{E})");
      fhPtHbpXEUeLeftNeutral->SetXTitle("p_{T trigger} (GeV/c)");
      outputContainer->Add(fhPtHbpXEUeLeftNeutral) ;
      
      fhPtHbpXEUeRightNeutral  = 
      new TH2F("hHbpXEUeNeutralRight","#xi = ln(1/x_{E}) with neutral UE right side of trigger",
               nptbins,ptmin,ptmax,200,0.,10.); 
      fhPtHbpXEUeRightNeutral->SetYTitle("ln(1/x_{E})");
      fhPtHbpXEUeRightNeutral->SetXTitle("p_{T trigger} (GeV/c)");
      outputContainer->Add(fhPtHbpXEUeRightNeutral) ;
      
      fhZTUeLeftNeutral  = 
      new TH2F("hZTUeNeutralLeft","z_{trigger h^{0}} = p_{T Ueh^{0}} / p_{T trigger} with neutral UE left side of trigger",
               nptbins,ptmin,ptmax,140,0.,2.); 
      fhZTUeLeftNeutral->SetYTitle("z_{trigger Ueh^{0}}");
      fhZTUeLeftNeutral->SetXTitle("p_{T trigger} (GeV/c)");
      outputContainer->Add(fhZTUeLeftNeutral) ;
      
      fhZTUeRightNeutral  = 
      new TH2F("hZTUeNeutralRight","z_{trigger h^{0}} = p_{T Ueh^{0}} / p_{T trigger} with neutral UE right side of trigger",
               nptbins,ptmin,ptmax,200,0.,2.); 
      fhZTUeRightNeutral->SetYTitle("z_{trigger Ueh^{0}}");
      fhZTUeRightNeutral->SetXTitle("p_{T trigger} (GeV/c)");
      outputContainer->Add(fhZTUeRightNeutral) ;
      
      fhPtHbpZTUeLeftNeutral  = 
      new TH2F("hHbpZTUeNeutralLeft","#xi = ln(1/z_{T}) with neutral UE left side of trigger",
               nptbins,ptmin,ptmax,200,0.,10.); 
      fhPtHbpZTUeLeftNeutral->SetYTitle("ln(1/z_{T})");
      fhPtHbpZTUeLeftNeutral->SetXTitle("p_{T trigger}");
      outputContainer->Add(fhPtHbpZTUeLeftNeutral) ;
      
      fhPtHbpZTUeRightNeutral  = 
      new TH2F("hHbpZTUeNeutralRight","#xi = ln(1/z_{T}) with neutral UE right side of trigger",
               nptbins,ptmin,ptmax,200,0.,10.); 
      fhPtHbpZTUeRightNeutral->SetYTitle("ln(1/z_{T})");
      fhPtHbpZTUeRightNeutral->SetXTitle("p_{T trigger}");
      outputContainer->Add(fhPtHbpZTUeRightNeutral) ;
      
    }  
        
  }//Correlation with neutral hadrons
  
  //if data is MC, fill more histograms
  if(IsDataMC())
  {
    fh2phiLeadingParticle=new TH2F("h2phiLeadingParticle","#phi resolustion for trigger particles",nptbins,ptmin,ptmax,100,-1,1);
    fh2phiLeadingParticle->GetXaxis()->SetTitle("p_{T gen Leading} (GeV/c)");
    fh2phiLeadingParticle->GetYaxis()->SetTitle("(#phi_{rec}-#phi_{gen})/#phi_{gen}");

    fhMCPtLeading  = new TH1F ("hMCPtLeading","MC : p_{T} distribution of leading particles", nptbins,ptmin,ptmax); 
    fhMCPtLeading->SetXTitle("p_{T}^{trig} (GeV/c)");

    fhMCPhiLeading  = new TH2F ("hMCPhiLeading","MC : #phi distribution of leading Particles",nptbins,ptmin,ptmax, nphibins,phimin,phimax); 
    fhMCPhiLeading->SetYTitle("#phi (rad)");
  
    fhMCEtaLeading  = new TH2F ("hMCEtaLeading","MC : #eta distribution of leading",nptbins,ptmin,ptmax, netabins,etamin,etamax); 
    fhMCEtaLeading->SetYTitle("#eta "); 
    
    
    fhMCEtaCharged  = new TH2F
    ("hMCEtaCharged","MC #eta_{h^{#pm}}  vs p_{T #pm}",
     nptbins,ptmin,ptmax,netabins,etamin,etamax); 
    fhMCEtaCharged->SetYTitle("#eta_{h^{#pm}} (rad)");
    fhMCEtaCharged->SetXTitle("p_{T #pm} (GeV/c)");
    
    fhMCPhiCharged  = new TH2F
    ("hMCPhiCharged","#MC phi_{h^{#pm}}  vs p_{T #pm}",
     200,ptmin,ptmax,nphibins,phimin,phimax); 
    fhMCPhiCharged->SetYTitle("MC #phi_{h^{#pm}} (rad)");
    fhMCPhiCharged->SetXTitle("p_{T #pm} (GeV/c)");
    
    fhMCDeltaPhiDeltaEtaCharged  = new TH2F
    ("hMCDeltaPhiDeltaEtaCharged","#MC phi_{trigger} - #phi_{h^{#pm}} vs #eta_{trigger} - #eta_{h^{#pm}}",
     140,-2.,5.,200,-2,2); 
    fhMCDeltaPhiDeltaEtaCharged->SetXTitle("#Delta #phi (rad)");
    fhMCDeltaPhiDeltaEtaCharged->SetYTitle("#Delta #eta");    
    
    fhMCDeltaEtaCharged  = new TH2F
    ("hMCDeltaEtaCharged","MC #eta_{trigger} - #eta_{h^{#pm}} vs p_{T trigger} and p_{T assoc}",
     nptbins,ptmin,ptmax,200,-2,2); 
    fhMCDeltaEtaCharged->SetYTitle("#Delta #eta");
    fhMCDeltaEtaCharged->SetXTitle("p_{T trigger} (GeV/c)");
    
    fhMCDeltaPhiCharged  = new TH2F
    ("hMCDeltaPhiCharged","#phi_{trigger} - #phi_{h^{#pm}} vs p_{T trigger}",
     nptbins,ptmin,ptmax,ndeltaphibins ,deltaphimin,deltaphimax); 
    fhMCDeltaPhiCharged->SetYTitle("#Delta #phi (rad)");
    fhMCDeltaPhiCharged->SetXTitle("p_{T trigger} (GeV/c)");
    
    fhMCDeltaPhiChargedPt  = new TH2F
    ("hMCDeltaPhiChargedPt","MC #phi_{trigger} - #phi_{#h^{#pm}} vs p_{T h^{#pm}}",
     nptbins,ptmin,ptmax,ndeltaphibins ,deltaphimin,deltaphimax); 
    fhMCDeltaPhiChargedPt->SetYTitle("#Delta #phi (rad)");
    fhMCDeltaPhiChargedPt->SetXTitle("p_{T h^{#pm}} (GeV/c)");
    
    fhMCPtXECharged  = 
    new TH2F("hMCPtXECharged","x_{E}",
             nptbins,ptmin,ptmax,200,0.,2.); 
    fhMCPtXECharged->SetYTitle("x_{E}");
    fhMCPtXECharged->SetXTitle("p_{T trigger} (GeV/c)");

    fhMCPtXEUeCharged  = 
    new TH2F("hMCPtXEUeCharged","x_{E}",
             nptbins,ptmin,ptmax,200,0.,2.); 
    fhMCPtXEUeCharged->SetYTitle("x_{E}");
    fhMCPtXEUeCharged->SetXTitle("p_{T trigger} (GeV/c)");
    
    fhMCPtHbpXECharged  = 
    new TH2F("hMCHbpXECharged","MC #xi = ln(1/x_{E}) with charged hadrons",
             nptbins,ptmin,ptmax,200,0.,10.); 
    fhMCPtHbpXECharged->SetYTitle("ln(1/x_{E})");
    fhMCPtHbpXECharged->SetXTitle("p_{T trigger}");

    fhMCPtHbpXEUeCharged =
    new TH2F("hMCPtHbpXEUeCharged","#xi = ln(1/x_{E}) with charged hadrons,Underlying Event",
             nptbins,ptmin,ptmax,200,0.,10.); 
    fhMCPtHbpXEUeCharged->SetYTitle("ln(1/x_{E})");
    fhMCPtHbpXEUeCharged->SetXTitle("p_{T trigger} (GeV/c)");

    fhMCUePart  = 
    new TH1F("hMCUePart","MC UE particles distribution vs pt trig",
             nptbins,ptmin,ptmax); 
    fhMCUePart->SetYTitle("dNch");
    fhMCUePart->SetXTitle("p_{T trigger} (GeV/c)");
    
    fhMCPtZTCharged  = 
    new TH2F("hMCPtZTCharged","z_{T}",
             nptbins,ptmin,ptmax,200,0.,2.); 
    fhMCPtZTCharged->SetYTitle("z_{T}");
    fhMCPtZTCharged->SetXTitle("p_{T trigger} (GeV/c)"); 
    
    fhMCPtHbpZTCharged  = 
    new TH2F("hMCHbpZTCharged","MC #xi = ln(1/z_{T}) with charged hadrons",
             nptbins,ptmin,ptmax,200,0.,10.); 
    fhMCPtHbpZTCharged->SetYTitle("ln(1/z_{T})");
    fhMCPtHbpZTCharged->SetXTitle("p_{T trigger} (GeV/c)");
    
    fhMCPtTrigPout  = 
    new TH2F("hMCPtTrigPout","AOD MC Pout with triggers",
             nptbins,ptmin,ptmax,2*nptbins,-ptmax,ptmax); 
    fhMCPtTrigPout->SetYTitle("p_{out} (GeV/c)");
    fhMCPtTrigPout->SetXTitle("p_{T trigger} (GeV/c)"); 
    
    fhMCPtAssocDeltaPhi  = 
    new TH2F("hMCPtAssocDeltaPhi","AOD MC delta phi with associated charged hadrons",
             nptbins,ptmin,ptmax,ndeltaphibins ,deltaphimin,deltaphimax); 
    fhMCPtAssocDeltaPhi->SetYTitle("#Delta #phi (rad)");
    fhMCPtAssocDeltaPhi->SetXTitle("p_{T trigger} (GeV/c)"); 
        
    outputContainer->Add(fh2phiLeadingParticle);
    outputContainer->Add(fhMCPtLeading);
    outputContainer->Add(fhMCPhiLeading);
    outputContainer->Add(fhMCEtaLeading);
    outputContainer->Add(fhMCDeltaPhiDeltaEtaCharged);
    outputContainer->Add(fhMCPhiCharged) ;
    outputContainer->Add(fhMCEtaCharged) ;
    outputContainer->Add(fhMCDeltaEtaCharged) ;
    outputContainer->Add(fhMCDeltaPhiCharged) ; 
    
    outputContainer->Add(fhMCDeltaPhiChargedPt) ;
    outputContainer->Add(fhMCPtXECharged) ;
    outputContainer->Add(fhMCPtXEUeCharged) ;
    outputContainer->Add(fhMCPtZTCharged) ;
    outputContainer->Add(fhMCPtHbpXECharged) ;
    outputContainer->Add(fhMCPtHbpXEUeCharged);
    outputContainer->Add(fhMCUePart);
    outputContainer->Add(fhMCPtHbpZTCharged) ;
    outputContainer->Add(fhMCPtTrigPout) ;
    outputContainer->Add(fhMCPtAssocDeltaPhi) ;      
  } //for MC histogram
  
  if(DoOwnMix())
  {
    //create event containers
    
    if(!fUseMixStoredInReader || (fUseMixStoredInReader && !GetReader()->ListWithMixedEventsForTracksExists())) 
    {
      Int_t nvz = GetNZvertBin();
      Int_t nrp = GetNRPBin();
      Int_t nce = GetNCentrBin();
      
      fListMixTrackEvents= new TList*[nvz*nrp*nce] ;
      
      for( Int_t ice = 0 ; ice < nce ; ice++ )
      {
        for( Int_t ivz = 0 ; ivz < nvz ; ivz++ )
        {
          for( Int_t irp = 0 ; irp < nrp ; irp++ )
          {
            Int_t bin = GetEventMixBin(ice,ivz,irp); //ic*nvz*nrp+iz*nrp+irp;
            
            //printf("GetCreateOutputObjects - Bins : cent %d, vz %d, RP %d, event %d/%d\n",
            //       ic,iz, irp, bin);
            
            fListMixTrackEvents[bin] = new TList() ;
            fListMixTrackEvents[bin]->SetOwner(kFALSE);
          }
        }
      }    
    }
    
    fhPtLeadingMixed  = new TH1F ("hPtLeadingMixed","p_{T} distribution of leading particles, used for mixing", nptbins,ptmin,ptmax); 
    fhPtLeadingMixed->SetXTitle("p_{T}^{trig} (GeV/c)");

    if(fCorrelVzBin)
    {
      fhPtLeadingMixedVzBin  = new TH2F ("hPtLeadingMixedVzBin","p_{T} distribution of leading particles, used for mixing", nptbins,ptmin,ptmax,GetNZvertBin(),0,GetNZvertBin()); 
      fhPtLeadingMixedVzBin->SetXTitle("p_{T}^{trig} (GeV/c)");
      fhPtLeadingMixedVzBin->SetYTitle("v_{z} bin");
      outputContainer->Add(fhPtLeadingMixedVzBin);
    }
    
    fhPtLeadingMixedBin  = new TH2F ("hPtLeadingMixedBin","p_{T} distribution of leading particles vs mixing bin", nptbins,ptmin,ptmax,nMixBins,0,nMixBins); 
    fhPtLeadingMixedBin->SetXTitle("p_{T}^{trig} (GeV/c)");
    fhPtLeadingMixedBin->SetYTitle("Bin");

    fhPhiLeadingMixed  = new TH2F ("hPhiLeadingMixed","#phi distribution of leading Particles, used for mixing",nptbins,ptmin,ptmax, nphibins,phimin,phimax); 
    fhPhiLeadingMixed->SetYTitle("#phi (rad)");
    
    fhEtaLeadingMixed  = new TH2F ("hEtaLeadingMixed","#eta distribution of leading, used for mixing",nptbins,ptmin,ptmax, netabins,etamin,etamax); 
    fhEtaLeadingMixed->SetYTitle("#eta ");  
    
    outputContainer->Add(fhPtLeadingMixed);
    outputContainer->Add(fhPtLeadingMixedBin);
    outputContainer->Add(fhPhiLeadingMixed);
    outputContainer->Add(fhEtaLeadingMixed);
    
    // Fill the cluster pool only in isolation analysis or if requested
    if( ( OnlyIsolated()        ||  fFillNeutralEventMixPool) &&
        (!fUseMixStoredInReader || (fUseMixStoredInReader && !GetReader()->ListWithMixedEventsForCaloExists())))
    {
      Int_t nvz = GetNZvertBin();
      Int_t nrp = GetNRPBin();
      Int_t nce = GetNCentrBin();
      
      fListMixCaloEvents= new TList*[nvz*nrp*nce] ;
      
      for( Int_t ice = 0 ; ice < nce ; ice++ )
      {
        for( Int_t ivz = 0 ; ivz < nvz ; ivz++ )
        {
          for( Int_t irp = 0 ; irp < nrp ; irp++ )
          {
            Int_t bin = GetEventMixBin(ice,ivz,irp); //ic*nvz*nrp+iz*nrp+irp;
            
            //printf("GetCreateOutputObjects - Bins : cent %d, vz %d, RP %d, event %d/%d\n",
            //       ic,iz, irp, bin);
            
            fListMixCaloEvents[bin] = new TList() ;
            fListMixCaloEvents[bin]->SetOwner(kFALSE);
          }
        }
      }    
    }
    
    //Init the list in the reader if not done previously
    if(fUseMixStoredInReader)
    {
      if( !GetReader()->ListWithMixedEventsForTracksExists() ) 
        GetReader()->SetListWithMixedEventsForTracks(fListMixTrackEvents);
      
      if( !GetReader()->ListWithMixedEventsForCaloExists()   ) 
        GetReader()->SetListWithMixedEventsForCalo  (fListMixCaloEvents );
    }
    
    fhEventBin=new TH1I("hEventBin","Number of real events per bin(cen,vz,rp)",
                        GetNCentrBin()*GetNZvertBin()*GetNRPBin()+1,0, 
                        GetNCentrBin()*GetNZvertBin()*GetNRPBin()+1) ;
    fhEventBin->SetXTitle("bin");
    outputContainer->Add(fhEventBin) ;
    
    fhEventMixBin=new TH1I("hEventMixBin","Number of events  per bin(cen,vz,rp)",
                           GetNCentrBin()*GetNZvertBin()*GetNRPBin()+1,0,
                           GetNCentrBin()*GetNZvertBin()*GetNRPBin()+1) ;
    fhEventMixBin->SetXTitle("bin");
    outputContainer->Add(fhEventMixBin) ;
    
    fhNtracksMB=new TH1F("hNtracksMBEvent","Number of tracks w/ event trigger kMB",2000,0,2000);
    outputContainer->Add(fhNtracksMB);
    
    if(fFillNeutralEventMixPool || OnlyIsolated())
    {
      fhNclustersMB=new TH1F("hNclustersMBEvent","Number of clusters w/ event trigger kMB",2000,0,2000);
      outputContainer->Add(fhNclustersMB);
    }
    
    fhMixDeltaPhiCharged  = new TH2F
    ("hMixDeltaPhiCharged","Mixed event : #phi_{trigger} - #phi_{h^{#pm}} vs p_{T trigger}",
     nptbins,ptmin,ptmax,ndeltaphibins ,deltaphimin,deltaphimax); 
    fhMixDeltaPhiCharged->SetYTitle("#Delta #phi (rad)");
    fhMixDeltaPhiCharged->SetXTitle("p_{T trigger} (GeV/c)");
    outputContainer->Add(fhMixDeltaPhiCharged);
    
    fhMixDeltaPhiDeltaEtaCharged  = new TH2F
    ("hMixDeltaPhiDeltaEtaCharged","Mixed event : #phi_{trigger} - #phi_{h^{#pm}} vs #eta_{trigger} - #eta_{h^{#pm}}",
     ndeltaphibins ,deltaphimin,deltaphimax,ndeltaetabins ,deltaetamin,deltaetamax); 
    fhMixDeltaPhiDeltaEtaCharged->SetXTitle("#Delta #phi (rad)");
    fhMixDeltaPhiDeltaEtaCharged->SetYTitle("#Delta #eta");
    outputContainer->Add(fhMixDeltaPhiDeltaEtaCharged);
    
    fhMixXECharged  = 
    new TH2F("hMixXECharged","Mixed event : x_{E} for charged tracks",
             nptbins,ptmin,ptmax,200,0.,2.); 
    fhMixXECharged->SetYTitle("x_{E}");
    fhMixXECharged->SetXTitle("p_{T trigger} (GeV/c)");
    outputContainer->Add(fhMixXECharged);
    
    fhMixXEUeCharged  =
    new TH2F("hMixXEUeCharged","Mixed event : x_{E} for charged tracks in Ue region",
             nptbins,ptmin,ptmax,200,0.,2.);
    fhMixXEUeCharged->SetYTitle("x_{E}");
    fhMixXEUeCharged->SetXTitle("p_{T trigger} (GeV/c)");
    outputContainer->Add(fhMixXEUeCharged);

    fhMixHbpXECharged  =
    new TH2F("hMixHbpXECharged","mixed event : #xi = ln(1/x_{E}) with charged hadrons",
             nptbins,ptmin,ptmax,200,0.,10.); 
    fhMixHbpXECharged->SetYTitle("ln(1/x_{E})");
    fhMixHbpXECharged->SetXTitle("p_{T trigger} (GeV/c)");
    outputContainer->Add(fhMixHbpXECharged);

    fhMixDeltaPhiChargedAssocPtBin         = new TH2F*[fNAssocPtBins*nz];
    fhMixDeltaPhiChargedAssocPtBinDEta08   = new TH2F*[fNAssocPtBins*nz];
    fhMixDeltaPhiChargedAssocPtBinDEta0    = new TH2F*[fNAssocPtBins*nz];
    fhMixDeltaPhiDeltaEtaChargedAssocPtBin = new TH2F*[fNAssocPtBins*nz];
    
    for(Int_t i = 0 ; i < fNAssocPtBins ; i++)
    {    
      for(Int_t z = 0 ; z < nz ; z++)
      {
        Int_t bin = i*nz+z;
        
        if(fCorrelVzBin)
        {
          sz = "_vz%d"+z;
          tz = ", v_{z} bin "+z;
        }
        
        //printf("MIX : iAssoc %d, Vz %d, bin %d - sz %s, tz %s	\n",i,z,bin,sz.Data(),tz.Data());
        
        fhMixDeltaPhiChargedAssocPtBin[bin] = new TH2F(Form("hMixDeltaPhiChargedAssocPtBin%2.1f_%2.1f%s", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1],sz.Data()), 
                                                     Form("Mixed event #Delta #phi vs p_{T trigger} for associated p_{T} bin [%2.1f,%2.1f]%s", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1],tz.Data()), 
                                                     nptbins, ptmin, ptmax,  ndeltaphibins ,deltaphimin,deltaphimax);
        fhMixDeltaPhiChargedAssocPtBin[bin]->SetXTitle("p_{T trigger} (GeV/c)");
        fhMixDeltaPhiChargedAssocPtBin[bin]->SetYTitle("#Delta #phi (rad)");
        
        fhMixDeltaPhiChargedAssocPtBinDEta08[bin] = new TH2F(Form("hMixDeltaPhiDeltaEta0.8ChargedAssocPtBin%2.1f_%2.1f%s", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1],sz.Data()), 
                                                           Form("Mixed event #Delta #phi vs p_{T trigger} for associated p_{T} bin [%2.1f,%2.1f]%s, for #Delta #eta > 0.8", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1],tz.Data()), 
                                                           nptbins, ptmin, ptmax,  ndeltaphibins ,deltaphimin,deltaphimax);
        fhMixDeltaPhiChargedAssocPtBinDEta08[bin]->SetXTitle("p_{T trigger} (GeV/c)");
        fhMixDeltaPhiChargedAssocPtBinDEta08[bin]->SetYTitle("#Delta #phi (rad)");      
        
        fhMixDeltaPhiChargedAssocPtBinDEta0[bin] = new TH2F(Form("hMixDeltaPhiDeltaEta0ChargedAssocPtBin%2.1f_%2.1f%s", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1],sz.Data()), 
                                                             Form("Mixed event #Delta #phi vs p_{T trigger} for associated p_{T} bin [%2.1f,%2.1f]%s, for #Delta #eta = 0", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1],tz.Data()), 
                                                             nptbins, ptmin, ptmax,  ndeltaphibins ,deltaphimin,deltaphimax);
        fhMixDeltaPhiChargedAssocPtBinDEta0[bin]->SetXTitle("p_{T trigger} (GeV/c)");
        fhMixDeltaPhiChargedAssocPtBinDEta0[bin]->SetYTitle("#Delta #phi (rad)");
        
        fhMixDeltaPhiDeltaEtaChargedAssocPtBin[bin] = new TH2F(Form("hMixDeltaPhiDeltaEtaChargedAssocPtBin%2.1f_%2.1f%s", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1],sz.Data()), 
                                                             Form("Mixed event #Delta #phi vs p_{T trigger} for associated p_{T} bin [%2.1f,%2.1f]%s", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1],tz.Data()), 
                                                             ndeltaphibins ,deltaphimin,deltaphimax,ndeltaetabins ,deltaetamin,deltaetamax); 
        fhMixDeltaPhiDeltaEtaChargedAssocPtBin[bin]->SetXTitle("#Delta #phi (rad)");
        fhMixDeltaPhiDeltaEtaChargedAssocPtBin[bin]->SetYTitle("#Delta #eta");
        
        outputContainer->Add(fhMixDeltaPhiChargedAssocPtBin[bin]);
        outputContainer->Add(fhMixDeltaPhiChargedAssocPtBinDEta08[bin]);
        outputContainer->Add(fhMixDeltaPhiChargedAssocPtBinDEta0[bin]);
        outputContainer->Add(fhMixDeltaPhiDeltaEtaChargedAssocPtBin[bin]);
      }
    }          
  }
  
  return outputContainer;
  
}

//_________________________________________________________________________________________________
Bool_t AliAnaParticleHadronCorrelation::GetDecayPhotonMomentum(const AliAODPWG4Particle* trigger, 
                                                             TLorentzVector & mom1, 
                                                             TLorentzVector & mom2)
{
  // Get the momentum of the pi0/eta assigned decay photons
  // In case of pi0/eta trigger, we may want to check their decay correlation, 
  // get their decay children
  
  Int_t indexPhoton1 = trigger->GetCaloLabel(0);
  Int_t indexPhoton2 = trigger->GetCaloLabel(1);
  Float_t ptTrig     = trigger->Pt();
  
  if(indexPhoton1!=-1 || indexPhoton2!=-1) return kFALSE;
  
  if(GetDebug() > 1) 
    printf("AliAnaParticleHadronCorrelation::GetDecayPhotonMomentum() - indexPhoton1 = %d, indexPhoton2 = %d \n", indexPhoton1, indexPhoton2);
  
  TObjArray * clusters  = 0x0 ;  
  if(trigger->GetDetector()=="EMCAL") clusters = GetEMCALClusters() ;
  else                                clusters = GetPHOSClusters()  ;
  
  for(Int_t iclus = 0; iclus < clusters->GetEntriesFast(); iclus++)
  {
    AliVCluster * photon =  (AliVCluster*) (clusters->At(iclus));	
    if(photon->GetID()==indexPhoton1) 
    {
      photon->GetMomentum(mom1,GetVertex(0)) ;
      if(ptTrig) fhPtPi0DecayRatio->Fill(ptTrig, mom1.Pt()/ptTrig);
    }
    if(photon->GetID()==indexPhoton2) 
    {
      photon->GetMomentum(mom1,GetVertex(0)) ;
      if(ptTrig > 0) fhPtPi0DecayRatio->Fill(ptTrig, mom2.Pt()/ptTrig);
    } 
    
    if(GetDebug() > 1)printf("AliAnaParticleHadronCorrelation::GetDecayPhotonMomentum() - Photon1 = %f, Photon2 = %f \n", mom1.Pt(), mom2.Pt());
    
  } //cluster loop        
  
  return kTRUE;
  
} 

//_____________________________________________________________
Int_t AliAnaParticleHadronCorrelation::GetMCTagHistogramIndex(Int_t mcTag)
{
  // Index of MC histograms depending on MC origin
  
  if     ( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPrompt) ||        
           GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCFragmentation)) return 0;
  else if( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPi0))           return 1;    
  else if( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPi0Decay))      return 2;
  else if( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCEtaDecay))      return 3;
  else if( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCOtherDecay))    return 4;
  else if(!GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCElectron))      return 5;
  else                                                                                    return 6;
  
}

//____________________________________________________
void AliAnaParticleHadronCorrelation::InitParameters()
{
  
  //Initialize the parameters of the analysis.
  SetInputAODName("Particle");
  SetAODObjArrayName("Hadrons");  
  AddToHistogramsName("AnaHadronCorr_");
  
  SetPtCutRange(0.,300);
  fDeltaPhiMinCut       = 1.5 ;
  fDeltaPhiMaxCut       = 4.5 ;
  fSelectIsolated       = kFALSE;
  fMakeSeveralUE        = kFALSE;
  fUeDeltaPhiMinCut     = 1. ;
  fUeDeltaPhiMaxCut     = 1.5 ;
  
  fNeutralCorr          = kFALSE ;
  fPi0Trigger           = kFALSE ;
  fDecayTrigger         = kFALSE ;
  fHMPIDCorrelation     = kFALSE ;
  
  fMakeAbsoluteLeading  = kTRUE;
  fMakeNearSideLeading  = kFALSE;

  fNAssocPtBins         = 9   ;
  fAssocPtBinLimit[0]   = 0.2 ; 
  fAssocPtBinLimit[1]   = 0.5 ; 
  fAssocPtBinLimit[2]   = 1.0 ; 
  fAssocPtBinLimit[3]   = 2.0 ; 
  fAssocPtBinLimit[4]   = 3.0 ; 
  fAssocPtBinLimit[5]   = 4.0 ; 
  fAssocPtBinLimit[6]   = 5.0 ;
  fAssocPtBinLimit[7]   = 6.0 ;
  fAssocPtBinLimit[8]   = 7.0 ;
  fAssocPtBinLimit[9]   = 8.0 ;
  fAssocPtBinLimit[10]  = 9.0 ; 
  fAssocPtBinLimit[11]  = 10.0 ; 
  fAssocPtBinLimit[12]  = 12.0 ; 
  fAssocPtBinLimit[13]  = 14.0 ; 
  fAssocPtBinLimit[14]  = 16.0 ; 
  fAssocPtBinLimit[15]  = 20.0 ; 
  fAssocPtBinLimit[16]  = 30.0 ;
  fAssocPtBinLimit[17]  = 40.0 ;
  fAssocPtBinLimit[18]  = 50.0 ;
  fAssocPtBinLimit[19]  = 200.0 ;
  
  
  fUseMixStoredInReader = kTRUE;
  
  fM02MinCut   = -1 ;
  fM02MaxCut   = -1 ;
  
}

//__________________________________________________________
void  AliAnaParticleHadronCorrelation::MakeAnalysisFillAOD()  
{  
  //Particle-Hadron Correlation Analysis, fill AODs
  
  if(!GetInputAODBranch())
  {
    printf("AliAnaParticleHadronCorrelation::MakeAnalysisFillAOD() - No input particles in AOD with name branch < %s >, STOP \n",GetInputAODName().Data());
    abort();
  }
	
  if(strcmp(GetInputAODBranch()->GetClass()->GetName(), "AliAODPWG4ParticleCorrelation"))
  {
    printf("AliAnaParticleHadronCorrelation::MakeAnalysisFillAOD() - Wrong type of AOD object, change AOD class name in input AOD: It should be <AliAODPWG4ParticleCorrelation> and not <%s> \n",GetInputAODBranch()->GetClass()->GetName());
    abort();
  }
	
  if(GetDebug() > 1)
  {
    printf("AliAnaParticleHadronCorrelation::MakeAnalysisFillAOD() - Begin hadron correlation analysis, fill AODs \n");
    printf("AliAnaParticleHadronCorrelation::MakeAnalysisFillAOD() - In particle branch aod entries %d\n", GetInputAODBranch()->GetEntriesFast());
    printf("AliAnaParticleHadronCorrelation::MakeAnalysisFillAOD() - In CTS aod entries %d\n",   GetCTSTracks()    ->GetEntriesFast());
    printf("AliAnaParticleHadronCorrelation::MakeAnalysisFillAOD() - In EMCAL aod entries %d\n", GetEMCALClusters()->GetEntriesFast());
    printf("AliAnaParticleHadronCorrelation::MakeAnalysisFillAOD() - In PHOS aod entries %d\n",  GetPHOSClusters() ->GetEntriesFast());
  }
  
  //Get the vertex and check it is not too large in z
  Double_t v[3] = {0,0,0}; //vertex ;
  GetReader()->GetVertex(v);
  if(!GetMixedEvent() && TMath::Abs(v[2]) > GetZvertexCut()) return ;   
    
  //Loop on stored AOD particles, find leading trigger
  Double_t ptTrig      = fMinTriggerPt ;
  fLeadingTriggerIndex = -1 ;
  Int_t    naod        = GetInputAODBranch()->GetEntriesFast() ;
  for(Int_t iaod = 0; iaod < naod ; iaod++)
  {
    AliAODPWG4ParticleCorrelation* particle =  (AliAODPWG4ParticleCorrelation*) (GetInputAODBranch()->At(iaod));
    
    // Vertex cut in case of mixing
    Int_t check = CheckMixedEventVertex(particle->GetCaloLabel(0), particle->GetTrackLabel(0));
    if(check ==  0) continue;
    if(check == -1) return;
        
    // find the leading particles with highest momentum
    if (particle->Pt() > ptTrig) 
    {
      ptTrig               = particle->Pt() ;
      fLeadingTriggerIndex = iaod ;
    }
  }// finish search of leading trigger particle
	
  
  //Do correlation with leading particle
  if(fLeadingTriggerIndex >= 0)
  {
	  
    AliAODPWG4ParticleCorrelation* particle =  (AliAODPWG4ParticleCorrelation*) (GetInputAODBranch()->At(fLeadingTriggerIndex));
    
    //check if the particle is isolated or if we want to take the isolation into account
    if(OnlyIsolated() && !particle->IsIsolated()) return;
    
    //Make correlation with charged hadrons
    Bool_t okcharged = kTRUE;
    Bool_t okneutral = kTRUE;
    if(GetReader()->IsCTSSwitchedOn() )
      okcharged = MakeChargedCorrelation(particle, GetCTSTracks(),kFALSE);
    
    TObjArray * pi0list = (TObjArray*) GetAODBranch(fPi0AODBranchName); //For the future, foresee more possible pi0 lists
    if(fNeutralCorr && pi0list && pi0list->GetEntriesFast() > 0)
      okneutral = MakeNeutralCorrelation(particle, pi0list,kFALSE);
    
  }//Correlate leading
  
  if(GetDebug() > 1) printf("AliAnaParticleHadronCorrelation::MakeAnalysisFillAOD() - End fill AODs \n");
  
}

//_________________________________________________________________
void  AliAnaParticleHadronCorrelation::MakeAnalysisFillHistograms()  
{  
  //Particle-Hadron Correlation Analysis, fill histograms
  
  if(!GetInputAODBranch())
  {
    printf("AliAnaParticleHadronCorrelation::MakeAnalysisFillHistograms() - No input particles in AOD with name branch < %s >, STOP \n",GetInputAODName().Data());
    abort();
  }
  
  if(GetDebug() > 1)
  {
    printf("AliAnaParticleHadronCorrelation::MakeAnalysisFillHistograms() - Begin hadron correlation analysis, fill histograms \n");
    printf("AliAnaParticleHadronCorrelation::MakeAnalysisFillHistograms() - In particle branch aod entries %d\n", GetInputAODBranch()->GetEntriesFast());
  }
    
  //Get the vertex and check it is not too large in z
  Double_t v[3] = {0,0,0}; //vertex ;
  GetReader()->GetVertex(v);
  if(!GetMixedEvent() && TMath::Abs(v[2]) > GetZvertexCut()) return ;  
  
  Float_t cen = GetEventCentrality();
  Float_t ep  = GetEventPlaneAngle();
  
  //Loop on stored AOD particles, find leading
  Double_t ptTrig    = fMinTriggerPt;
  if(fLeadingTriggerIndex < 0)
  {
    //Search leading if not done before
    Int_t    naod      = GetInputAODBranch()->GetEntriesFast() ;
    for(Int_t iaod = 0; iaod < naod ; iaod++)
    {	 //loop on input trigger AOD file 
      AliAODPWG4ParticleCorrelation* particle =  (AliAODPWG4ParticleCorrelation*) (GetInputAODBranch()->At(iaod));

      // Vertex cut in case of mixing
      Int_t check = CheckMixedEventVertex(particle->GetCaloLabel(0), particle->GetTrackLabel(0));
      if(check ==  0) continue;
      if(check == -1) return;
            
      //check if the particle is isolated or if we want to take the isolation into account
      if(OnlyIsolated() && !particle->IsIsolated()) continue;
      
      //find the leading particles with highest momentum
      if (particle->Pt() > ptTrig) 
      {
        ptTrig               = particle->Pt() ;
        fLeadingTriggerIndex = iaod ;
      }
      
    }// Finish search of leading trigger particle
  }// Search leading if not done before
  
  if(fLeadingTriggerIndex >= 0 )
  { //using trigger particle to do correlations
    
    AliAODPWG4ParticleCorrelation* particle =  (AliAODPWG4ParticleCorrelation*) (GetInputAODBranch()->At(fLeadingTriggerIndex));

    // check if it was a calorimeter cluster and if the SS cut was requested, if so, apply it
    Int_t clID1  = particle->GetCaloLabel(0) ;
    Int_t clID2  = particle->GetCaloLabel(1) ; // for photon clusters should not be set.
    //printf("Leading for for %s: id1 %d, id2 %d, min %f, max %f, det %s\n",
    //       GetInputAODName().Data(),clID1,clID2,fM02MinCut,fM02MaxCut,(particle->GetDetector()).Data());

    if(clID1 > 0 && clID2 < 0 && fM02MaxCut > 0 && fM02MinCut > 0)
    {
      Int_t iclus = -1;
      TObjArray* clusters = 0x0;
      if     (particle->GetDetector() == "EMCAL") clusters = GetEMCALClusters();
      else if(particle->GetDetector() == "PHOS" ) clusters = GetPHOSClusters();
      
      if(clusters)
      {
        AliVCluster *cluster = FindCluster(clusters,clID1,iclus); 
        Float_t m02 = cluster->GetM02();
        //printf("\t Check m02 = %2.2f\n",m02);
        if(m02 > fM02MaxCut || m02 < fM02MinCut) 
        {
          //printf("\t \t Not accepted\n");
          return;
        }
      }        
    }
    
    // Check if the particle is isolated or if we want to take the isolation into account
    if(OnlyIsolated() && !particle->IsIsolated()) return;
    
    Float_t pt = particle->Pt();
    fhPtInput->Fill(pt);
    
    // Check if trigger is in fiducial region
    if(IsFiducialCutOn())
    {
      Bool_t in = GetFiducialCut()->IsInFiducialCut(*particle->Momentum(),particle->GetDetector()) ;
      if(! in ) return ;
    }
    
    fhPtFidCut->Fill(pt);
        
    // Make correlation with charged hadrons
    Bool_t okcharged = kTRUE;
    Bool_t okneutral = kTRUE;
    if(GetReader()->IsCTSSwitchedOn() )
    {
      okcharged = MakeChargedCorrelation(particle, GetCTSTracks(),kTRUE);
      if(IsDataMC())
      {      
        MakeMCChargedCorrelation(particle);
      }
    }  
    
    TObjArray * pi0list = (TObjArray*) GetAODBranch(fPi0AODBranchName); //For the future, foresee more possible pi0 lists
    if(fNeutralCorr && pi0list)
    {
      if(pi0list->GetEntriesFast() > 0)
        okneutral = MakeNeutralCorrelation(particle, pi0list,kTRUE);
    }
    
    // Fill leading particle histogram if correlation went well and 
    // no problem was found, like not absolute leading, or bad vertex in mixing.
    if(okcharged && okneutral)
    {
      fhPtLeading->Fill(pt);
      fhPtLeadingBin->Fill(pt,GetEventMixBin());
      if(fCorrelVzBin) fhPtLeadingVzBin->Fill(pt,GetEventVzBin());
 
      if(fFillPileUpHistograms)
      {
        Int_t vtxBC = GetReader()->GetVertexBC();
        if(vtxBC == 0 || vtxBC==AliVTrack::kTOFBCNA)     fhPtLeadingVtxBC0->Fill(pt);

        if(GetReader()->IsPileUpFromSPD())               fhPtLeadingPileUp[0]->Fill(pt);
        if(GetReader()->IsPileUpFromEMCal())             fhPtLeadingPileUp[1]->Fill(pt);
        if(GetReader()->IsPileUpFromSPDOrEMCal())        fhPtLeadingPileUp[2]->Fill(pt);
        if(GetReader()->IsPileUpFromSPDAndEMCal())       fhPtLeadingPileUp[3]->Fill(pt);
        if(GetReader()->IsPileUpFromSPDAndNotEMCal())    fhPtLeadingPileUp[4]->Fill(pt);
        if(GetReader()->IsPileUpFromEMCalAndNotSPD())    fhPtLeadingPileUp[5]->Fill(pt);
        if(GetReader()->IsPileUpFromNotSPDAndNotEMCal()) fhPtLeadingPileUp[6]->Fill(pt);
      }
      
      Float_t phi = particle->Phi();
      if(phi<0)phi+=TMath::TwoPi();
      fhPhiLeading->Fill(pt, phi);
      
      fhEtaLeading->Fill(pt, particle->Eta());
      //printf("AliAnaParticleHadronCorrelation::MakeAnalysisFillHistograms() - Leading particle : pt %f, eta %f, phi %f\n",particle->Pt(),particle->Eta(),phi);
      
      if(IsDataMC())
      {
        Int_t mcIndex = GetMCTagHistogramIndex(particle->GetTag());
        fhPtLeadingMC[mcIndex]->Fill(pt);
      }        
            
      fhPtLeadingCentrality        ->Fill(pt,cen);
      fhPtLeadingEventPlane        ->Fill(pt,ep);
      fhLeadingEventPlaneCentrality->Fill(cen,ep);
      
    }//ok charged && neutral
  }//Aod branch loop
  
  //Reinit for next event
  fLeadingTriggerIndex = -1;
  
  if(GetDebug() > 1) printf("AliAnaParticleHadronCorrelation::MakeAnalysisFillHistograms() - End fill histograms \n");
}

//___________________________________________________________________________________________________________
Bool_t  AliAnaParticleHadronCorrelation::MakeChargedCorrelation(AliAODPWG4ParticleCorrelation *aodParticle, 
                                                                const TObjArray* pl, Bool_t bFillHisto)
{  
  // Charged Hadron Correlation Analysis
  if(GetDebug() > 1) 
    printf("AliAnaParticleHadronCorrelation::MakeChargedCorrelation() - Make trigger particle - charged hadron correlation \n");
    
  Float_t phiTrig = aodParticle->Phi();
  Float_t etaTrig = aodParticle->Eta();
  Float_t ptTrig  = aodParticle->Pt();  
  Bool_t   decay  = aodParticle->IsTagged();
  Int_t    mcTag  = aodParticle->GetTag();
  Double_t bz     = GetReader()->GetInputEvent()->GetMagneticField();

  Float_t pt       = -100. ;
  Float_t zT       = -100. ; 
  Float_t xE       = -100. ; 
  Float_t hbpXE    = -100. ; 
  Float_t hbpZT    = -100. ; 
  Float_t phi      = -100. ;
  Float_t eta      = -100. ;
  Float_t pout     = -100. ;
  Float_t deltaPhi = -100. ;
  
  TVector3 p3;  
  TLorentzVector photonMom ;	
  TObjArray * reftracks = 0x0;
  Int_t nrefs           = 0;
  Int_t nTracks         = GetCTSTracks()->GetEntriesFast() ;
  
  // Mixed event settings
  Int_t evtIndex11   = -1 ; // cluster trigger or pi0 trigger 
  Int_t evtIndex12   = -1 ; // pi0 trigger
  Int_t evtIndex13   = -1 ; // charged trigger
  
  Double_t v[3]      = {0,0,0}; //vertex ;
  GetReader()->GetVertex(v);
  
  if (GetMixedEvent()) 
  {
    evtIndex11 = GetMixedEvent()->EventIndexForCaloCluster(aodParticle->GetCaloLabel(0)) ;
    evtIndex12 = GetMixedEvent()->EventIndexForCaloCluster(aodParticle->GetCaloLabel(1)) ;    
    evtIndex13 = GetMixedEvent()->EventIndex(aodParticle->GetTrackLabel(0)) ;
  }
  
  // In case of pi0/eta trigger, we may want to check their decay correlation, 
  // get their decay children
  TLorentzVector decayMom1;
  TLorentzVector decayMom2;
  Bool_t decayFound = kFALSE;
  if(fPi0Trigger && bFillHisto) decayFound = GetDecayPhotonMomentum(aodParticle,decayMom1, decayMom2);

  //-----------------------------------------------------------------------
  //Track loop, select tracks with good pt, phi and fill AODs or histograms
  //-----------------------------------------------------------------------

  for(Int_t ipr = 0;ipr < pl->GetEntriesFast() ; ipr ++ )
  {
    AliVTrack * track = (AliVTrack *) (pl->At(ipr)) ;
    
    Double_t mom[3] = {track->Px(),track->Py(),track->Pz()};
    p3.SetXYZ(mom[0],mom[1],mom[2]);
    pt   = p3.Pt();
    eta  = p3.Eta();
    phi  = p3.Phi() ;
    if(phi < 0) phi+=TMath::TwoPi();
    
    //Select only hadrons in pt range
    if(pt < fMinAssocPt || pt > fMaxAssocPt) continue ;
    
    //remove trigger itself for correlation when use charged triggers    
    if( track->GetID() == aodParticle->GetTrackLabel(0) || track->GetID() == aodParticle->GetTrackLabel(1) ||
        track->GetID() == aodParticle->GetTrackLabel(2) || track->GetID() == aodParticle->GetTrackLabel(3)   ) 
      continue ;
    
    //jump out this event if near side associated particle pt larger than trigger
    if (fMakeNearSideLeading)
    {
      if(pt > ptTrig && TMath::Abs(phi-phiTrig) < TMath::PiOver2())  return kFALSE;
    }
    //jump out this event if there is any other particle with pt larger than trigger
    else if(fMakeAbsoluteLeading)
    {
      if(pt > ptTrig)  return kFALSE;
    }
    
    //Only for mixed event
    Int_t evtIndex2 = 0 ; 
    if (GetMixedEvent()) 
    {
      evtIndex2 = GetMixedEvent()->EventIndex(track->GetID()) ;
      if (evtIndex11 == evtIndex2 || evtIndex12 == evtIndex2 || evtIndex13 == evtIndex2 ) // photon and track from different events
        continue ; 
      //vertex cut
      if (TMath::Abs(GetVertex(evtIndex2)[2]) > GetZvertexCut()) 
        return kFALSE;
    }    
    
    // Fill Histograms
    if(bFillHisto)
    {      

      if(GetDebug() > 2 ) 
        printf("AliAnaParticleHadronCorrelation::MakeChargedCorrelation() - Selected charge for momentum imbalance: pt %2.2f, phi %2.2f, eta %2.2f \n",pt,phi,eta);
            
      // Set the pt associated bin for the defined bins
      Int_t assocBin   = -1; 
      for(Int_t i = 0 ; i < fNAssocPtBins ; i++)
      {
        if(pt > fAssocPtBinLimit[i] && pt < fAssocPtBinLimit[i+1]) assocBin= i; 
      }      
      
      // Assign to the histogram array a bin corresponding to a combination of pTa and vz bins
      Int_t nz = 1;
      Int_t vz = 0;
      
      if(fCorrelVzBin) 
      {
        nz = GetNZvertBin();
        vz = GetEventVzBin();
      }
      
      Int_t bin = assocBin*nz+vz;
      
      //printf("assoc Bin = %d, vZ bin  = %d, bin = %d \n", assocBin,GetEventVzBin(),bin);
      
      ULong_t status = track->GetStatus();
      Bool_t okTOF = ( (status & AliVTrack::kTOFout) == AliVTrack::kTOFout ) ;
      //Double32_t tof = track->GetTOFsignal()*1e-3;
      Int_t trackBC = track->GetTOFBunchCrossing(bz);

      Int_t outTOF = -1;
      if     (okTOF && trackBC!=0) outTOF = 1;
      else if(okTOF && trackBC==0) outTOF = 0;
      
      // Azimuthal Angle
      // calculate deltaPhi for later, shift when needed
      FillChargedAngularCorrelationHistograms(pt,  ptTrig,  bin, phi, phiTrig,  deltaPhi,
                                              eta, etaTrig, decay, track->GetHMPIDsignal(),outTOF,nTracks,mcTag);
      
      // Imbalance zT/xE/pOut
      zT = pt/ptTrig ;
      if(zT > 0 ) hbpZT = TMath::Log(1./zT);
      else        hbpZT =-100;
      
      xE   =-pt/ptTrig*TMath::Cos(deltaPhi); // -(px*pxTrig+py*pyTrig)/(ptTrig*ptTrig);
      //if(xE <0.)xE =-xE;
      if(xE > 0 ) hbpXE = TMath::Log(1./xE); 
      else        hbpXE =-100;
    
      pout = pt*TMath::Sin(deltaPhi) ;
      
      //delta phi cut for momentum imbalance correlation
      if  ( (deltaPhi > fDeltaPhiMinCut)   && (deltaPhi < fDeltaPhiMaxCut)   )
      {
        
        FillChargedMomentumImbalanceHistograms(ptTrig, pt, xE, hbpXE, zT, hbpZT, pout, 
                                               nTracks, track->Charge(), bin, decay,outTOF,mcTag);
        
      }
      
      if ( (deltaPhi > fUeDeltaPhiMinCut) && (deltaPhi < fUeDeltaPhiMaxCut) )
      { //UE study
        
        FillChargedUnderlyingEventHistograms(ptTrig, pt, deltaPhi, nTracks,outTOF);

        fhUePart->Fill(ptTrig);
        
      }
      
      if(fPi0Trigger && decayFound) 
        FillDecayPhotonCorrelationHistograms(pt, phi, decayMom1,decayMom2, kTRUE) ;
      
      //several UE calculation 
      if(fMakeSeveralUE) FillChargedUnderlyingEventSidesHistograms(ptTrig,pt,deltaPhi);
      
    } //Fill histogram 
    else
    {
      nrefs++;
      if(nrefs==1)
      {
        reftracks = new TObjArray(0);
        TString trackname = Form("%s+Tracks", GetAODObjArrayName().Data());
        reftracks->SetName(trackname.Data());
        reftracks->SetOwner(kFALSE);        
      }
      
      reftracks->Add(track);
      
    }//aod particle loop
  }// track loop
  
  //Fill AOD with reference tracks, if not filling histograms
  if(!bFillHisto && reftracks) 
  {
    aodParticle->AddObjArray(reftracks);
  }

  //Own mixed event, add event and remove previous or fill the mixed histograms
  if(DoOwnMix() && bFillHisto)
  {
      MakeChargedMixCorrelation(aodParticle);
  }
  
  return kTRUE;
  
}  


//_________________________________________________________________________________________________________
void AliAnaParticleHadronCorrelation::MakeChargedMixCorrelation(AliAODPWG4ParticleCorrelation *aodParticle) 
{  
  // Mix current trigger with tracks in another MB event
  
  if(GetDebug() > 1) printf("AliAnaParticleHadronCorrelationNew::MakeChargedMixCorrelation() - Make trigger particle - charged hadron mixed event correlation \n");
  
  if(GetMixedEvent()) return;  // This is not the mixed event from general mixing frame
  
  // Get the event with similar caracteristics
  //printf("MakeChargedMixCorrelation for %s\n",GetInputAODName().Data());

  AliAnalysisManager   * manager      = AliAnalysisManager::GetAnalysisManager();
  
  AliInputEventHandler * inputHandler = dynamic_cast<AliInputEventHandler*>(manager->GetInputEventHandler());
  
  if(!inputHandler) return;
  
  if(!(inputHandler->IsEventSelected( ) & GetReader()->GetEventTriggerMask())) return;
    
  // Get the pool, check if it exits
  Int_t eventBin = GetEventMixBin();

  fhEventBin->Fill(eventBin);
  
  //Check that the bin exists, if not (bad determination of RP, centrality or vz bin) do nothing
  if(eventBin < 0) return;
  
  TList * pool     = 0;
  TList * poolCalo = 0;
  if(fUseMixStoredInReader) 
  {
    pool     = GetReader()->GetListWithMixedEventsForTracks(eventBin);
    if(OnlyIsolated() || fFillNeutralEventMixPool) poolCalo = GetReader()->GetListWithMixedEventsForCalo  (eventBin);
  }
  else
  {
    pool     = fListMixTrackEvents[eventBin];
    if(OnlyIsolated()  || fFillNeutralEventMixPool) poolCalo = fListMixCaloEvents [eventBin];
  }
  
  if(!pool) return ;
    
  if((OnlyIsolated()  || fFillNeutralEventMixPool ) && !poolCalo &&
     (GetIsolationCut()->GetParticleTypeInCone()!=AliIsolationCut::AliIsolationCut::kOnlyCharged)) 
    printf("AliAnaParticleHadronCorrelation::MakeChargedMixCorrelation() - Careful, cluster pool not available\n");
  
  Double_t ptTrig  = aodParticle->Pt();
  Double_t etaTrig = aodParticle->Eta();
  Double_t phiTrig = aodParticle->Phi();
  if(phiTrig < 0.) phiTrig+=TMath::TwoPi();
  
  if(GetDebug() > 1) 
    printf("AliAnaParticleHadronCorrelationNew::MakeChargedMixCorrelation() - Pool bin %d size %d, leading trigger pt=%f, phi=%f, eta=%f\n",
           eventBin,pool->GetSize(), ptTrig,phiTrig,etaTrig);
  
  Double_t ptAssoc  = -999.;
  Double_t phiAssoc = -999.;
  Double_t etaAssoc = -999.;
  Double_t deltaPhi = -999.;
  Double_t deltaEta = -999.;
  Double_t xE = -999.;
  Double_t hbpXE = -999.;
      
  //Start from first event in pool except if in this same event the pool was filled
  Int_t ev0 = 0;
  if(GetReader()->GetLastTracksMixedEvent() == GetEventNumber()) ev0 = 1;

  for(Int_t ev=ev0; ev < pool->GetSize(); ev++)
  {
    TObjArray* bgTracks = static_cast<TObjArray*>(pool->At(ev));
    TObjArray* bgCalo   = 0;

    // Check if the particle is isolated in the mixed event, it not, do not fill the histograms
    if(OnlyIsolated() || fFillNeutralEventMixPool)
    {
      if(pool->GetSize()!=poolCalo->GetSize()) 
        printf("AliAnaParticleHadronCorrelationNew::MakeChargedMixCorrelation() - Different size of calo and track pools\n");
      
      bgCalo = static_cast<TObjArray*>(poolCalo->At(ev));
      
      if(!bgCalo) 
        printf("AliAnaParticleHadronCorrelationNew::MakeChargedMixCorrelation() - Event %d in calo pool not available?\n",ev);
      
      if(OnlyIsolated() && bgCalo)
      {
        Int_t n=0; Int_t nfrac = 0; Bool_t isolated = kFALSE; Float_t coneptsum = 0;
        GetIsolationCut()->MakeIsolationCut(bgTracks,bgCalo,
                                            GetReader(), GetCaloPID(),
                                            kFALSE, aodParticle, "", 
                                            n,nfrac,coneptsum, isolated);
        
        //printf("AliAnaParticleHadronCorrelation::MakeChargedMixCorrelation() - Isolated? %d - cone %f, ptthres %f",
        //       isolated,GetIsolationCut()->GetConeSize(),GetIsolationCut()->GetPtThreshold());
        //if(bgTracks)printf(" - n track %d", bgTracks->GetEntriesFast());
        //printf("\n");
        
        if(!isolated) continue ;
      }
    }
    
    fhEventMixBin->Fill(eventBin);
    
    Int_t nTracks=bgTracks->GetEntriesFast();
    //printf("\t Read Pool event %d, nTracks %d\n",ev,nTracks);

    //Check if it is leading if mixed event
    if(fMakeNearSideLeading || fMakeAbsoluteLeading)
    {
      Bool_t leading = kTRUE;
      for(Int_t jlead = 0;jlead <nTracks; jlead++ )
      {
        AliAODPWG4Particle *track = (AliAODPWG4Particle*) bgTracks->At(jlead) ;
        
        ptAssoc  = track->Pt();
        phiAssoc = track->Phi() ;
        
        if(phiAssoc < 0) phiAssoc+=TMath::TwoPi();
        if (fMakeNearSideLeading)
        {
          if(ptAssoc > ptTrig && TMath::Abs(phiAssoc-phiTrig) < TMath::PiOver2())  
          {
            leading = kFALSE;
            break;
          }
        }
        //jump out this event if there is any other particle with pt larger than trigger
        else if(fMakeAbsoluteLeading)
        {
          if(ptAssoc > ptTrig) 
          { 
            leading = kFALSE;
            break;
          }
        }
      }
      
      if(fFillNeutralEventMixPool && bgCalo)
      {
        Int_t nClusters=bgCalo->GetEntriesFast();
        TLorentzVector mom ;
        for(Int_t jlead = 0;jlead <nClusters; jlead++ )
        {
          AliVCluster *cluster = (AliVCluster*) bgCalo->At(jlead) ;
          
          Double_t vertex[]={0,0,0}; // assume 0 vertex
          cluster->GetMomentum(mom,vertex) ;

          ptAssoc  = mom.Pt();
          phiAssoc = mom.Phi() ;
          
          if(phiAssoc < 0) phiAssoc+=TMath::TwoPi();
          if (fMakeNearSideLeading)
          {
            if(ptAssoc > ptTrig && TMath::Abs(phiAssoc-phiTrig) < TMath::PiOver2())
            {
              leading = kFALSE;
              break;
            }
          }
          //jump out this event if there is any other particle with pt larger than trigger
          else if(fMakeAbsoluteLeading)
          {
            if(ptAssoc > ptTrig)
            {
              leading = kFALSE;
              break;
            }
          }
        }
      }
      
      if(!leading) continue; // not leading, check the next event in pool
    
    }
    
    fhPtLeadingMixed   ->Fill(ptTrig);
    fhPhiLeadingMixed  ->Fill(ptTrig, phiTrig);
    fhEtaLeadingMixed  ->Fill(ptTrig, etaTrig);
    fhPtLeadingMixedBin->Fill(ptTrig,eventBin);
    if(fCorrelVzBin)fhPtLeadingMixedVzBin->Fill(ptTrig, GetEventVzBin());

    for(Int_t j1 = 0;j1 <nTracks; j1++ )
    {
      AliAODPWG4Particle *track = (AliAODPWG4Particle*) bgTracks->At(j1) ;
      
      if(!track) continue;
      
      ptAssoc  = track->Pt();
      etaAssoc = track->Eta();
      phiAssoc = track->Phi() ;
      if(phiAssoc < 0) phiAssoc+=TMath::TwoPi();
            
      if(IsFiducialCutOn())
      {
        Bool_t in = GetFiducialCut()->IsInFiducialCut(*aodParticle->Momentum(),"CTS") ;
        if(!in) continue ;
      }
      
      deltaPhi = phiTrig-phiAssoc;
      if(deltaPhi < -TMath::PiOver2())  deltaPhi+=TMath::TwoPi();
      if(deltaPhi > 3*TMath::PiOver2()) deltaPhi-=TMath::TwoPi();
      deltaEta = etaTrig-etaAssoc;
      
      if(GetDebug()>0)
        printf("AliAnaParticleHadronCorrelationNew::MakeChargedMixCorrelation(): deltaPhi= %f, deltaEta=%f\n",deltaPhi, deltaEta);
      
      // Set the pt associated bin for the defined bins
      Int_t assocBin   = -1; 
      for(Int_t i = 0 ; i < fNAssocPtBins ; i++)
      {
        if(ptAssoc > fAssocPtBinLimit[i] && ptAssoc < fAssocPtBinLimit[i+1]) assocBin= i; 
      }      

      // Assign to the histogram array a bin corresponding to a combination of pTa and vz bins
      Int_t nz = 1;
      Int_t vz = 0;
      
      if(fCorrelVzBin) 
      {
        nz = GetNZvertBin();
        vz = GetEventVzBin();
      }
      
      Int_t bin = assocBin*nz+vz;
      
      fhMixDeltaPhiCharged        ->Fill(ptTrig,  deltaPhi);
      fhMixDeltaPhiDeltaEtaCharged->Fill(deltaPhi, deltaEta);

      fhMixDeltaPhiCharged        ->Fill(ptTrig,  deltaPhi);
      fhMixDeltaPhiDeltaEtaCharged->Fill(deltaPhi, deltaEta);

      xE   =-ptAssoc/ptTrig*TMath::Cos(deltaPhi); // -(px*pxTrig+py*pyTrig)/(ptTrig*ptTrig);
      //if(xE <0.)xE =-xE;
      if(xE > 0 ) hbpXE = TMath::Log(1./xE); 
      else        hbpXE =-100;

      if ( (deltaPhi > fDeltaPhiMinCut)   && (deltaPhi < fDeltaPhiMaxCut)   )
      {
        fhMixXECharged->Fill(ptTrig,xE);
        fhMixHbpXECharged->Fill(ptTrig,hbpXE);
      }
      
      if ( (deltaPhi > fUeDeltaPhiMinCut) && (deltaPhi < fUeDeltaPhiMaxCut) )
      {
        //Underlying event region
        Double_t randomphi = gRandom->Uniform(fDeltaPhiMinCut,fDeltaPhiMaxCut);
        Double_t uexE = -(ptAssoc/ptTrig)*TMath::Cos(randomphi);
        
        if(uexE < 0.) uexE = -uexE;
        
        fhMixXEUeCharged->Fill(ptTrig,uexE);
      }
      
      if(bin < 0) continue ; // this pt bin was not considered
      
      if(TMath::Abs(deltaEta) > 0.8) 
        fhMixDeltaPhiChargedAssocPtBinDEta08  [bin]->Fill(ptTrig,   deltaPhi);
      if(TMath::Abs(deltaEta) < 0.01) 
        fhMixDeltaPhiChargedAssocPtBinDEta0   [bin]->Fill(ptTrig,   deltaPhi);
      
        fhMixDeltaPhiChargedAssocPtBin        [bin]->Fill(ptTrig,   deltaPhi);
        fhMixDeltaPhiDeltaEtaChargedAssocPtBin[bin]->Fill(deltaPhi, deltaEta);
      
    } // track loop
  } // mixed event loop
}
  

//________________________________________________________________________________________________________________
Bool_t  AliAnaParticleHadronCorrelation::MakeNeutralCorrelation(AliAODPWG4ParticleCorrelation * const aodParticle, 
                                                                const TObjArray* pi0list, Bool_t bFillHisto)  
{  
  // Neutral Pion Correlation Analysis
  if(GetDebug() > 1) printf("AliAnaParticleHadronCorrelation::MakeNeutralCorrelation() - Make trigger particle - pi0 correlation, %d pi0's \n",
                            pi0list->GetEntriesFast());
  
  Int_t evtIndex11 = 0 ; 
  Int_t evtIndex12 = 0 ; 
  if (GetMixedEvent()) 
  {
    evtIndex11 = GetMixedEvent()->EventIndexForCaloCluster(aodParticle->GetCaloLabel(0)) ;
    evtIndex12 = GetMixedEvent()->EventIndexForCaloCluster(aodParticle->GetCaloLabel(1)) ;    
  }  
  
  Float_t pt   = -100. ;
  Float_t zT   = -100. ; 
  Float_t phi  = -100. ;
  Float_t eta  = -100. ;
  Float_t xE   = -100. ; 
  Float_t hbpXE= -100. ; 
  Float_t hbpZT= -100. ; 

  Float_t ptTrig  = aodParticle->Pt();
  Float_t phiTrig = aodParticle->Phi();
  Float_t etaTrig = aodParticle->Eta();
  Float_t deltaPhi= -100. ;

  TLorentzVector photonMom ;
	
  // In case of pi0/eta trigger, we may want to check their decay correlation, 
  // get their decay children
  TLorentzVector decayMom1;
  TLorentzVector decayMom2;
  Bool_t decayFound = kFALSE;
  if(fPi0Trigger && bFillHisto) decayFound = GetDecayPhotonMomentum(aodParticle,decayMom1, decayMom2);  
  
  TObjArray * refpi0 = 0x0;
  Int_t nrefs        = 0;
  
  //Loop on stored AOD pi0
  
  Int_t naod = pi0list->GetEntriesFast();
  if(GetDebug() > 0) 
    printf("AliAnaParticleHadronCorrelation::MakeNeutralCorrelationFillHistograms() -  aod branch entries %d\n", naod);
  
  for(Int_t iaod = 0; iaod < naod ; iaod++)
  {
    AliAODPWG4Particle* pi0 =  (AliAODPWG4Particle*) (pi0list->At(iaod));
    
    Int_t evtIndex2 = 0 ; 
    Int_t evtIndex3 = 0 ; 
    if (GetMixedEvent()) 
    {
      evtIndex2 = GetMixedEvent()->EventIndexForCaloCluster(pi0->GetCaloLabel(0)) ;
      evtIndex3 = GetMixedEvent()->EventIndexForCaloCluster(pi0->GetCaloLabel(1)) ;
      
      if (evtIndex11 == evtIndex2 || evtIndex12 == evtIndex2 || 
          evtIndex11 == evtIndex3 || evtIndex12 == evtIndex3) // trigger and pi0 are not from different events
        continue ; 
    }      

    pt  = pi0->Pt();
     
    if(pt < fMinAssocPt || pt > fMaxAssocPt) continue ;
    
    //jump out this event if near side associated particle pt larger than trigger
    if (fMakeNearSideLeading)
    {
      if(pt > ptTrig && TMath::Abs(phi-phiTrig) < TMath::PiOver2())  return kFALSE;
    }
    //jump out this event if there is any other particle with pt larger than trigger
    else if(fMakeAbsoluteLeading)
    {
      if(pt > ptTrig)  return kFALSE;
    }
    
    if(bFillHisto)
    {
      phi = pi0->Phi() ;
      eta = pi0->Eta() ;
      
      FillNeutralAngularCorrelationHistograms(pt, ptTrig, phi, phiTrig, deltaPhi, eta, etaTrig);
      
      zT  = pt/ptTrig ;
      xE  =-pt/ptTrig*TMath::Cos(deltaPhi); // -(px*pxTrig+py*pyTrig)/(ptTrig*ptTrig);
      
      //if(xE <0.)xE =-xE;
      
      hbpXE = -100;
      hbpZT = -100;
      
      if(xE > 0 ) hbpXE = TMath::Log(1./xE); 
      if(zT > 0 ) hbpZT = TMath::Log(1./zT); 
      
      if(fPi0Trigger && decayFound)
        FillDecayPhotonCorrelationHistograms(pt, phi, decayMom1,decayMom2,kFALSE) ;
      
      //delta phi cut for correlation
      if( (deltaPhi > fDeltaPhiMinCut) && ( deltaPhi < fDeltaPhiMaxCut) ) 
      {
        fhDeltaPhiNeutralPt->Fill(pt,deltaPhi);
        fhXENeutral        ->Fill(ptTrig,xE); 
        fhPtHbpXENeutral   ->Fill(ptTrig,hbpXE); 
      }
      else if ( (deltaPhi > fUeDeltaPhiMinCut) && (deltaPhi < fUeDeltaPhiMaxCut) )      
      {
        fhDeltaPhiUeNeutralPt->Fill(pt,deltaPhi);
        fhXEUeNeutral        ->Fill(ptTrig,xE);
        fhPtHbpXEUeNeutral   ->Fill(ptTrig,hbpXE); 
      }
      
      //several UE calculation 
      if(fMakeSeveralUE) FillChargedUnderlyingEventSidesHistograms(ptTrig,pt,deltaPhi);

	  }
    else
    {
      nrefs++;
      if(nrefs==1)
      {
        refpi0 = new TObjArray(0);
        refpi0->SetName(GetAODObjArrayName()+"Pi0s");
        refpi0->SetOwner(kFALSE);
      }
      refpi0->Add(pi0);
    }//put references in trigger AOD 
    
    if(GetDebug() > 2 ) 
      printf("AliAnaParticleHadronCorrelation::MakeNeutralCorrelation() - Selected neutral for momentum imbalance: pt %2.2f, phi %2.2f, eta %2.2f \n",pt,phi,eta);
    
  }//loop
  
  return kTRUE;
}
  
//_________________________________________________________________________________________________________
void  AliAnaParticleHadronCorrelation::MakeMCChargedCorrelation(AliAODPWG4ParticleCorrelation *aodParticle)
{  
  // Charged Hadron Correlation Analysis with MC information
  
  if(GetDebug()>1)
    printf("AliAnaParticleHadronCorrelation::MakeMCChargedCorrelation() - Make trigger particle - charged hadron correlation in AOD MC level\n");
  
  AliStack         * stack        = 0x0 ;
  TParticle        * primary      = 0x0 ;   
  TClonesArray     * mcparticles0 = 0x0 ;
  TClonesArray     * mcparticles  = 0x0 ;
  AliAODMCParticle * aodprimary   = 0x0 ; 
  
  Double_t eprim   = 0 ;
  Double_t ptprim  = 0 ;
  Double_t phiprim = 0 ;
  Double_t etaprim = 0 ;
  Int_t    nTracks = 0 ;  
  Int_t iParticle  = 0 ;
  Double_t charge  = 0.;
  
  if(GetReader()->ReadStack())
  {
    nTracks = GetMCStack()->GetNtrack() ;
  }
  else 
  {
    nTracks = GetReader()->GetAODMCParticles()->GetEntriesFast() ;
  }
  //Int_t trackIndex[nTracks];
  
  Int_t label= aodParticle->GetLabel();
  if(label < 0)
  {
    if(GetDebug() > 0) printf("AliAnaParticleHadronCorrelation::MakeMCChargedCorrelation *** bad label ***:  label %d \n", label);
    return;
  }  
  
  if(GetReader()->ReadStack())
  {
    stack =  GetMCStack() ;
    if(!stack) 
    {
      printf(" AliAnaParticleHadronCorrelation::MakeMCChargedCorrelation- Stack not available, is the MC handler called? STOP\n");
      abort();
    }
    
    nTracks=stack->GetNprimary();
    if(label >=  stack->GetNtrack()) 
    {
      if(GetDebug() > 2)  printf("AliAnaParticleHadronCorrelation::MakeMCChargedCorrelation *** large label ***:  label %d, n tracks %d \n", label, stack->GetNtrack());
      return ;
    }
    
    primary = stack->Particle(label);
    if(!primary)
    {
      printf("AliAnaParticleHadronCorrelation::MakeMCChargedCorrelation *** no primary ***:  label %d \n", label);   
      return;
    }
    
    if(primary)
    {
      eprim    = primary->Energy();
      ptprim   = primary->Pt();
      phiprim  = primary->Phi();
      etaprim  = primary->Eta();
      
      if(ptprim < 0.01 || eprim < 0.01) return ;
      
      for (iParticle = 0 ; iParticle <  nTracks ; iParticle++) 
      {
        TParticle * particle = stack->Particle(iParticle);
        TLorentzVector momentum;
        
        //keep only final state particles
        if(particle->GetStatusCode()!=1) continue ;
        
        Int_t pdg = particle->GetPdgCode();	
        
        charge = TDatabasePDG::Instance()->GetParticle(pdg)->Charge();
        
        particle->Momentum(momentum);	
        
        //---------- Charged particles ----------------------
        if(charge != 0)
        {   
          //Particles in CTS acceptance
          Bool_t inCTS =  GetFiducialCut()->IsInFiducialCut(momentum,"CTS");
          
          if(TMath::Abs(pdg) == 11 && stack->Particle(particle->GetFirstMother())->GetPdgCode()==22) continue ;
          
          if(inCTS)
          {            
            if( label!=iParticle) // avoid trigger particle
            {
              if(!FillChargedMCCorrelationHistograms(particle->Pt(),particle->Phi(),particle->Eta(),ptprim,phiprim,etaprim)) return;
            }
          }// in CTS acceptance
        }// charged
      } //track loop
    } //when the leading particles could trace back to MC
  } //ESD MC
  
  else if(GetReader()->ReadAODMCParticles())
  {
    //Get the list of MC particles
    mcparticles0 = GetReader()->GetAODMCParticles();
    if(!mcparticles0) return;
    
    if(label >=mcparticles0->GetEntriesFast())
    {
      if(GetDebug() > 2)  
        printf("AliAnaParticleHadronCorrelation::MakeMCChargedCorrelation *** large label ***:  label %d, n tracks %d \n", label,mcparticles0->GetEntriesFast());
      return;
    }
    
    //Get the particle
    aodprimary = (AliAODMCParticle*) mcparticles0->At(label);
    if(!aodprimary)  
    {
      printf("AliAnaParticleHadronCorrelation::MakeMCChargedCorrelation *** no AOD primary ***:  label %d \n", label);   
      return;
    }
    
    if(aodprimary)
    {
      ptprim  = aodprimary->Pt();
      phiprim = aodprimary->Phi();
      etaprim = aodprimary->Eta();
      eprim   = aodprimary->E();
      
      Bool_t lead = kFALSE;
      
      if(ptprim < 0.01 || eprim < 0.01) return ;
      
      mcparticles= GetReader()->GetAODMCParticles();
      for (iParticle = 0; iParticle < nTracks; iParticle++) 
      {
        AliAODMCParticle *part = (AliAODMCParticle*) mcparticles->At(iParticle);
        
        if (!part->IsPhysicalPrimary()) continue;        
   
        Int_t pdg = part->GetPdgCode();	
        charge = TDatabasePDG::Instance()->GetParticle(pdg)->Charge();
        TLorentzVector momentum(part->Px(),part->Py(),part->Pz(),part->E());   
        
        if(charge != 0)
        {
          if(part->Pt()> GetReader()->GetCTSPtMin())
          {
            //Particles in CTS acceptance
            Bool_t inCTS =  GetFiducialCut()->IsInFiducialCut(momentum,"CTS");
            Int_t indexmother=part->GetMother();
            
            if(indexmother>-1)
            {
              Int_t mPdg = ((AliAODMCParticle*) mcparticles->At(indexmother)) ->GetPdgCode();
              if(TMath::Abs(pdg) == 11 && mPdg == 22) continue;
            }
            
            if(inCTS)
            {            
              if( label!=iParticle) // avoid trigger particle
              {
                if(!FillChargedMCCorrelationHistograms(part->Pt(),part->Phi(),part->Eta(),ptprim,phiprim,etaprim)) return;
                else lead = kTRUE;
              }
            } // in acceptance
          } // min pt cut
        } //only charged particles
      }  //MC particle loop 
      if (lead) 
      {
        fhMCPtLeading->Fill(ptprim);
        fhMCPhiLeading->Fill(ptprim,phiprim);
        fhMCEtaLeading->Fill(ptprim,etaprim);
      }
    } //when the leading particles could trace back to MC
  }// AOD MC
}

//_____________________________________________________________________
void AliAnaParticleHadronCorrelation::Print(const Option_t * opt) const
{
  
  //Print some relevant parameters set for the analysis
  if(! opt)
    return;
  
  printf("**** Print %s %s ****\n", GetName(), GetTitle() ) ;
  AliAnaCaloTrackCorrBaseClass::Print(" ");
  printf("Pt trigger           >  %3.2f\n", fMinTriggerPt) ;
  printf("Pt associated hadron <  %3.2f\n", fMaxAssocPt) ; 
  printf("Pt associated hadron >  %3.2f\n", fMinAssocPt) ;
  printf("Phi trigger particle-Hadron   <  %3.2f\n", fDeltaPhiMaxCut) ; 
  printf("Phi trigger particle-Hadron   >  %3.2f\n", fDeltaPhiMinCut) ;
  printf("Phi trigger particle-UeHadron <  %3.2f\n", fUeDeltaPhiMaxCut) ; 
  printf("Phi trigger particle-UeHadron >  %3.2f\n", fUeDeltaPhiMinCut) ;
  printf("Isolated Trigger?  %d\n"     , fSelectIsolated) ;
  printf("Several UE?  %d\n"           , fMakeSeveralUE) ;
  printf("Name of AOD Pi0 Branch %s \n", fPi0AODBranchName.Data());
  printf("Do Decay-hadron correlation ?  %d\n", fPi0Trigger) ;
  printf("Select absolute leading for cluster triggers ?  %d or Near Side %d\n", 
         fMakeAbsoluteLeading, fMakeNearSideLeading) ;
  printf("Trigger pt bins  %d\n", fNAssocPtBins) ;
  for (Int_t ibin = 0; ibin<fNAssocPtBins; ibin++) {
    printf("\t bin %d = [%2.1f,%2.1f]\n", ibin, fAssocPtBinLimit[ibin], fAssocPtBinLimit[ibin+1]) ;
  }
  
} 

//____________________________________________________________
void AliAnaParticleHadronCorrelation::SetNAssocPtBins(Int_t n)
{
  // Set number of bins
  
  fNAssocPtBins  = n ; 
  
  
  if(n < 20 && n > 0)
  {
    fNAssocPtBins  = n ; 
  }
  else 
  {
    printf("n = larger than 19 or too small, set to 19 \n");
    fNAssocPtBins = 19;
  }
}

//______________________________________________________________________________
void AliAnaParticleHadronCorrelation::SetAssocPtBinLimit(Int_t ibin, Float_t pt)
{ 
  // Set the list of limits for the trigger pt bins
  
  if(ibin <= fNAssocPtBins || ibin >= 0) 
  {
    fAssocPtBinLimit[ibin] = pt  ;
  }
  else 
  {
    printf("AliAnaParticleHadronCorrelation::SetAssocPtBinLimit() - bin  number too large %d > %d or small, nothing done\n", ibin, fNAssocPtBins) ; 
  }
}

