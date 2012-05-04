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
    fDoOwnMix(0),                   fUseTrackMultBins(0),
    fListMixEvents(),               
    //Histograms
    fhPtLeading(0),                 fhPhiLeading(0),       
    fhEtaLeading(0),                fhDeltaPhiDeltaEtaCharged(0),
    fhPhiCharged(0),                fhEtaCharged(0), 
    fhDeltaPhiCharged(0),           fhDeltaEtaCharged(0), 
    fhDeltaPhiChargedPt(0),         fhDeltaPhiUeChargedPt(0), 
    fhXECharged(0),                 fhXEUeCharged(0),
    fhXEPosCharged(0),              fhXENegCharged(0),
    fhPtHbpXECharged(0),            fhPtHbpXEUeCharged(0),
    fhZTCharged(0),                 fhZTUeCharged(0),
    fhZTPosCharged(0),              fhZTNegCharged(0),
    fhPtHbpZTCharged(0),            fhPtHbpZTUeCharged(0),
    fhDeltaPhiUeLeftCharged(0),     fhDeltaPhiUeRightCharged(0),
    fhXEUeLeftCharged(0),           fhXEUeRightCharged(0),
    fhPtHbpXEUeLeftCharged(0),      fhPtHbpXEUeRightCharged(0), 
    fhZTUeLeftCharged(0),           fhZTUeRightCharged(0),
    fhPtHbpZTUeLeftCharged(0),      fhPtHbpZTUeRightCharged(0), 
    fhPtTrigPout(0),                fhPtTrigCharged(0),
    fhTrigDeltaPhiCharged(0x0),     fhTrigDeltaEtaCharged(0x0),
    fhTrigXECorr(0x0),              fhTrigXEUeCorr(0x0),
    fhTrigZTCorr(0x0),              fhTrigZTUeCorr(0x0),
    fhAssocPtBkg(0),  
    fhDeltaPhiAssocPtBin(0),        fhDeltaPhiAssocPtBinHMPID(0), fhDeltaPhiAssocPtBinHMPIDAcc(0),
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
    fh2phiLeadingParticle(0x0),
    fhMCEtaCharged(0),              fhMCPhiCharged(0), 
    fhMCDeltaEtaCharged(0),         fhMCDeltaPhiCharged(0x0),
    fhMCDeltaPhiDeltaEtaCharged(0), fhMCDeltaPhiChargedPt(0),
    fhMCPtXECharged(0),             fhMCPtHbpXECharged(0),          
    fhMCPtZTCharged(0),             fhMCPtHbpZTCharged(0),
    fhMCPtTrigPout(0),              fhMCPtAssocDeltaPhi(0),
    //Mixing
    fhNEventsTrigger(0),            
    fhNtracksAll(0),                fhNtracksTrigger(0),            
    fhNtracksINT(0),               
    fhMixDeltaPhiCharged(0),        fhMixDeltaPhiDeltaEtaCharged(0),
    fhMixDeltaPhiChargedAssocPtBin(),fhMixDeltaPhiDeltaEtaChargedAssocPtBin(0)
{ 
  //Default Ctor
    
  //Initialize parameters
  InitParameters();
}

//_________________________________________________________________
AliAnaParticleHadronCorrelation::~AliAnaParticleHadronCorrelation() 
{
  // Remove event containers
  
  if(fDoOwnMix && fListMixEvents)
  {      
    for(Int_t iz=0; iz < GetNZvertBin(); iz++)
    {
      for(Int_t ic=0; ic < GetNCentrBin(); ic++)
      {
        fListMixEvents[ic*GetNZvertBin()+iz]->Delete() ;
        delete fListMixEvents[ic*GetNZvertBin()+iz] ;
      }
    }
    
    delete[] fListMixEvents;
    
  }
}

//______________________________________________________________________________________________________________________________________________________
void AliAnaParticleHadronCorrelation::FillChargedAngularCorrelationHistograms(const Float_t ptAssoc,  const Float_t ptTrig,   const Int_t   assocBin,
                                                                              const Float_t phiAssoc, const Float_t phiTrig,  Float_t &     deltaPhi,
                                                                              const Float_t etaAssoc, const Float_t etaTrig,  
                                                                              const Bool_t  decay,    const Float_t hmpidSignal, const Int_t nTracks)
{
  // Fill angular correlation related histograms
  
  Float_t deltaEta    = etaTrig-etaAssoc;
          deltaPhi    = phiTrig-phiAssoc;
  Float_t deltaPhiOrg = deltaPhi;
  
  if(deltaPhi <= -TMath::PiOver2()) deltaPhi+=TMath::TwoPi();
  if(deltaPhi > 3*TMath::PiOver2()) deltaPhi-=TMath::TwoPi();
  
  fhEtaCharged     ->Fill(ptAssoc,etaAssoc);
  fhPhiCharged     ->Fill(ptAssoc,phiAssoc);
  fhDeltaEtaCharged->Fill(ptTrig ,deltaEta);
  fhDeltaPhiCharged->Fill(ptTrig ,deltaPhi);
  
  if(ptAssoc > 2 )           fhDeltaPhiDeltaEtaCharged->Fill(deltaPhi, deltaEta);
  
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
  if(assocBin>=0)
  {
    fhDeltaPhiAssocPtBin[assocBin]->Fill(ptTrig, deltaPhi);
    if (fFillBradHisto)        fhDeltaPhiBradAssocPtBin        [assocBin]->Fill(ptTrig, dphiBrad);
    if(fDecayTrigger && decay) fhDeltaPhiDecayChargedAssocPtBin[assocBin]->Fill(ptTrig, deltaPhi);      
    
    if(fHMPIDCorrelation)
    {
      if( hmpidSignal > 0 )
      {
        //printf("Track pt %f with HMPID signal %f \n",pt,hmpidSignal);
        fhDeltaPhiAssocPtBinHMPID[assocBin]->Fill(ptTrig, deltaPhi);        
      }
      
      if(phiAssoc > 5*TMath::DegToRad() && phiAssoc < 20*TMath::DegToRad())
      {
        //printf("Track pt %f in HMPID acceptance phi %f \n ",pt,phi*TMath::RadToDeg() );
        fhDeltaPhiAssocPtBinHMPIDAcc[assocBin]->Fill(ptTrig, deltaPhi);        
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
Bool_t AliAnaParticleHadronCorrelation::FillChargedMCCorrelationHistograms(const Float_t mcAssocPt,      Float_t mcAssocPhi, const Float_t mcAssocEta,
                                                                           const Float_t mcTrigPt, const Float_t mcTrigPhi,  const Float_t mcTrigEta)
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
  
  return kTRUE;
} 

//___________________________________________________________________________________________________________________
void AliAnaParticleHadronCorrelation::FillChargedMomentumImbalanceHistograms(const Float_t ptTrig,   const Float_t ptAssoc, 
                                                                             const Float_t xE,       const Float_t hbpXE, 
                                                                             const Float_t zT,       const Float_t hbpZT, 
                                                                             const Float_t pout,     const Float_t deltaPhi,
                                                                             const Int_t   nTracks,  const Int_t   charge,
                                                                             const Int_t   assocBin, const Bool_t  decay    )

{
  // Fill mostly momentum imbalance related histograms
  
  fhDeltaPhiChargedPt ->Fill(ptAssoc, deltaPhi);
  fhXECharged         ->Fill(ptTrig , xE); 
  fhPtHbpXECharged    ->Fill(ptTrig , hbpXE);
  fhZTCharged         ->Fill(ptTrig , zT); 
  fhPtHbpZTCharged    ->Fill(ptTrig , hbpZT);
  fhPtTrigPout        ->Fill(ptTrig , pout) ;
  fhPtTrigCharged     ->Fill(ptTrig , ptAssoc) ;
  
  if(fDecayTrigger && decay)
  {          
    fhXEDecayCharged->Fill(ptTrig,xE); 
    fhZTDecayCharged->Fill(ptTrig,zT);
  } // photon decay pi0/eta trigger        
  
  if(assocBin >= 0 )//away side 
  {
    fhXEAssocPtBin[assocBin]->Fill(ptTrig, xE) ;
    fhZTAssocPtBin[assocBin]->Fill(ptTrig, zT) ;
    
    if(fDecayTrigger && decay)
    {          
      fhXEDecayChargedAssocPtBin[assocBin]->Fill(ptTrig, xE); 
      fhZTDecayChargedAssocPtBin[assocBin]->Fill(ptTrig, zT);
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
void AliAnaParticleHadronCorrelation::FillChargedUnderlyingEventHistograms(const Float_t ptTrig,   const Float_t ptAssoc,
                                                                           const Float_t deltaPhi, const Int_t nTracks)
{
  // Fill underlying event histograms
  
  fhDeltaPhiUeChargedPt->Fill(ptAssoc,deltaPhi);
  
  Double_t randomphi = gRandom->Uniform(TMath::Pi()/2,3*TMath::Pi()/2);
  Double_t uexE = -(ptAssoc/ptTrig)*TMath::Cos(randomphi);
  Double_t uezT =   ptAssoc/ptTrig;
  
  if(uexE < 0.) uexE = -uexE;
    
  fhXEUeCharged->Fill(ptTrig,uexE);
  if(uexE > 0) fhPtHbpXEUeCharged->Fill(ptTrig,TMath::Log(1/uexE));
  
  fhZTUeCharged->Fill(ptTrig,uezT);
  if(uexE > 0) fhPtHbpZTUeCharged->Fill(ptTrig,TMath::Log(1/uezT));
  
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

//___________________________________________________________________________________________________________________________
void AliAnaParticleHadronCorrelation::FillChargedUnderlyingEventSidesHistograms(const Float_t ptTrig,   const Float_t ptAssoc, 
                                                                                const Float_t xE,       const Float_t hbpXE, 
                                                                                const Float_t zT,       const Float_t hbpZT, 
                                                                                const Float_t deltaPhi)
{
  // Fill underlying event histograms to the left and right of trigger
  
  if((deltaPhi<-fUeDeltaPhiMinCut) && (deltaPhi >-fUeDeltaPhiMaxCut))
  {  
    fhDeltaPhiUeLeftCharged->Fill(ptAssoc, deltaPhi);
    fhXEUeLeftCharged      ->Fill(ptTrig , xE);
    fhPtHbpXEUeLeftCharged ->Fill(ptTrig , hbpXE);
    fhZTUeLeftCharged      ->Fill(ptTrig , zT);
    fhPtHbpZTUeLeftCharged ->Fill(ptTrig , hbpZT);
  }
  
  if((deltaPhi > fUeDeltaPhiMinCut) && (deltaPhi < fUeDeltaPhiMaxCut))
  {  
    fhDeltaPhiUeRightCharged->Fill(ptAssoc, deltaPhi);
    fhXEUeRightCharged      ->Fill(ptTrig , xE);
    fhPtHbpXEUeRightCharged ->Fill(ptTrig , hbpXE);
    fhZTUeRightCharged      ->Fill(ptTrig , zT);
    fhPtHbpZTUeRightCharged ->Fill(ptTrig , hbpZT);
  }
} 

//______________________________________________________________________________________________________________________________
void AliAnaParticleHadronCorrelation::FillDecayPhotonCorrelationHistograms(const Float_t ptAssoc,     const Float_t phiAssoc, 
                                                                           const TLorentzVector mom1, const TLorentzVector mom2,
                                                                           const Bool_t bChargedOrNeutral)
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
    
    if(GetDebug() > 1)printf("AliAnaParticleHadronCorrelation::FillDecayPhotonHistograms( Charged corr) - deltaPhoton1 = %f, deltaPhoton2 = %f \n", deltaPhiDecay1, deltaPhiDecay2);
    
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
    
    if(GetDebug() > 1)printf("AliAnaParticleHadronCorrelation::FillDecayPhotonHistograms(Neutral corr) - deltaPhoton1 = %f, deltaPhoton2 = %f \n", deltaPhiDecay1, deltaPhiDecay2);
    
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
void AliAnaParticleHadronCorrelation::FillNeutralAngularCorrelationHistograms(const Float_t ptAssoc,  const Float_t ptTrig,  
                                                                              const Float_t phiAssoc, const Float_t phiTrig,  Float_t &     deltaPhi,
                                                                              const Float_t etaAssoc, const Float_t etaTrig)
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
void AliAnaParticleHadronCorrelation::FillNeutralUnderlyingEventSidesHistograms(const Float_t ptTrig,   const Float_t ptAssoc, 
                                                                                const Float_t xE,       const Float_t hbpXE, 
                                                                                const Float_t zT,       const Float_t hbpZT, 
                                                                                const Float_t deltaPhi)
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

//_____________________________________________________________
void AliAnaParticleHadronCorrelation::FillChargedEventMixPool()
{
  // Mixed event init
    
  Int_t nTracks = GetCTSTracks()->GetEntriesFast();
  
  fhNtracksAll->Fill(nTracks);
  
  AliAnalysisManager   * manager      = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler * inputHandler = dynamic_cast<AliInputEventHandler*>(manager->GetInputEventHandler());
  
  if(!inputHandler) return ;
  
  if(inputHandler->IsEventSelected( ) & GetReader()->GetEventTriggerMask())
  {
    fhNtracksTrigger->Fill(nTracks);
  }
  
  if(inputHandler->IsEventSelected( ) & AliVEvent::kAnyINT)
  {
    
    fhNtracksINT->Fill(nTracks);
    
    //Get vertex z bin
    Double_t v[3] = {0,0,0}; //vertex ;
    GetReader()->GetVertex(v);
    
    Int_t curZvertBin = (Int_t)(0.5*GetNZvertBin()*(v[2]+GetZvertexCut())/GetZvertexCut()) ;
    
    // centrality or tracks bin
    Int_t curCentrBin = 0;
    if(fUseTrackMultBins)
    { // Track multiplicity bins
      //curCentrBin = (GetTrackMultiplicity()-1)/5; 
      //if(curCentrBin > GetNCentrBin()-1) curCentrBin=GetNCentrBin()-1;
      Int_t trackMult = GetReader()->GetTrackMultiplicity();
      if(trackMult<=5)
        curCentrBin=8;
      else if(trackMult<=10)
        curCentrBin=7;
      else if(trackMult<=15)
        curCentrBin=6;
      else if(trackMult<=20)
        curCentrBin=5;
      else if(trackMult<=30)
        curCentrBin=4;
      else if(trackMult<=40)
        curCentrBin=3;
      else if(trackMult<=55)
        curCentrBin=2;
      else if(trackMult<=70)
        curCentrBin=1 ;
      else curCentrBin=0 ;        
    }
    else // Set centrality based on centrality task
    {
      curCentrBin = GetEventCentrality() * GetNCentrBin() / GetReader()->GetCentralityOpt(); 
      if(GetDebug() > 0 )printf("AliAnaParticleHadronCorrelation::MakeAnalysisFillHistograms() - curCentrBin %d, centrality %d, n bins %d, max bin from centrality %d\n",
                                curCentrBin, GetEventCentrality(), GetNCentrBin(), GetReader()->GetCentralityOpt());        
    }
    
    TObjArray * mixEventTracks = new TObjArray;
    
    //printf("curCen %d, curZ %d, bin %d\n",curCentrBin,curZvertBin,curCentrBin*GetNZvertBin()+curZvertBin);
    
    if(!fListMixEvents[curCentrBin*GetNZvertBin()+curZvertBin]) 
      fListMixEvents[curCentrBin*GetNZvertBin()+curZvertBin]=new TList();
    
    TList * pool = fListMixEvents[curCentrBin*GetNZvertBin()+curZvertBin];
    
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
    
    pool->AddFirst(mixEventTracks);
    mixEventTracks = 0;
    if(pool->GetSize() >= GetNMaxEvMix())
    {//Remove last event
      TClonesArray * tmp = static_cast<TClonesArray*>(pool->Last()) ;
      pool->RemoveLast() ;
      delete tmp ;
    }
    
  } // MB event
  
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
  
  fhPtLeading  = new TH1F ("hPtLeading","p_T distribution of leading particles", nptbins,ptmin,ptmax); 
  fhPtLeading->SetXTitle("p_{T}^{trig} (GeV/c)");
  
  fhPhiLeading  = new TH2F ("hPhiLeading","#phi distribution of leading Particles",nptbins,ptmin,ptmax, nphibins,phimin,phimax); 
  fhPhiLeading->SetYTitle("#phi (rad)");
  
  fhEtaLeading  = new TH2F ("hEtaLeading","#eta distribution of leading",nptbins,ptmin,ptmax, netabins,etamin,etamax); 
  fhEtaLeading->SetYTitle("#eta ");  
  
  outputContainer->Add(fhPtLeading);
  outputContainer->Add(fhPhiLeading);
  outputContainer->Add(fhEtaLeading);
  
  //Correlation with charged hadrons
  if(GetReader()->IsCTSSwitchedOn()) 
  {
    fhDeltaPhiDeltaEtaCharged  = new TH2F
    ("hDeltaPhiDeltaEtaCharged","#eta_{trigger} - #eta_{h^{#pm}} vs #phi_{trigger} - #phi_{h^{#pm}}",
    ndeltaphibins ,deltaphimin,deltaphimax,ndeltaetabins,deltaetamin,deltaetamax); 
    fhDeltaPhiDeltaEtaCharged->SetXTitle("#Delta #phi");
    fhDeltaPhiDeltaEtaCharged->SetYTitle("#Delta #eta");    
    
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
    fhDeltaPhiCharged->SetYTitle("#Delta #phi");
    fhDeltaPhiCharged->SetXTitle("p_{T trigger} (GeV/c)");
    
    fhDeltaPhiChargedPt  = new TH2F
    ("hDeltaPhiChargedPt","#phi_{trigger} - #phi_{#h^{#pm}} vs p_{T h^{#pm}}",
     nptbins,ptmin,ptmax, ndeltaphibins ,deltaphimin,deltaphimax);
    fhDeltaPhiChargedPt->SetYTitle("#Delta #phi");
    fhDeltaPhiChargedPt->SetXTitle("p_{T h^{#pm}} (GeV/c)");
    
    fhDeltaPhiUeChargedPt  = new TH2F
    ("hDeltaPhiUeChargedPt","#phi_{trigger} - #phi_{#Ueh^{#pm}} vs p_{T Ueh^{#pm}}",
     nptbins,ptmin,ptmax, ndeltaphibins ,deltaphimin,deltaphimax);
    fhDeltaPhiUeChargedPt->SetYTitle("#Delta #phi");
    fhDeltaPhiUeChargedPt->SetXTitle("p_{T h^{#pm}} (GeV/c)");
    
    fhDeltaEtaCharged  = new TH2F
    ("hDeltaEtaCharged","#eta_{trigger} - #eta_{h^{#pm}} vs p_{T trigger}",
     nptbins,ptmin,ptmax,ndeltaetabins,deltaetamin,deltaetamax);  
    fhDeltaEtaCharged->SetYTitle("#Delta #eta");
    fhDeltaEtaCharged->SetXTitle("p_{T trigger} (GeV/c)");
    
    fhXECharged  = 
    new TH2F("hXECharged","x_{E} for charged tracks",
             nptbins,ptmin,ptmax,200,0.,2.); 
    fhXECharged->SetYTitle("x_{E}");
    fhXECharged->SetXTitle("p_{T trigger}");
    
    fhXEUeCharged  = 
    new TH2F("hXEUeCharged","x_{E} for Underlying Event",
             nptbins,ptmin,ptmax,200,0.,2.); 
    fhXEUeCharged->SetYTitle("x_{E}");
    fhXEUeCharged->SetXTitle("p_{T trigger}");
    
    fhXEPosCharged  = 
    new TH2F("hXEPositiveCharged","x_{E} for positive charged tracks",
             nptbins,ptmin,ptmax,200,0.,2.); 
    fhXEPosCharged->SetYTitle("x_{E}");
    fhXEPosCharged->SetXTitle("p_{T trigger}");
    
    fhXENegCharged  = 
    new TH2F("hXENegativeCharged","x_{E} for negative charged tracks",
             nptbins,ptmin,ptmax,200,0.,2.); 
    fhXENegCharged->SetYTitle("x_{E}");
    fhXENegCharged->SetXTitle("p_{T trigger}");
    
    fhPtHbpXECharged  = 
    new TH2F("hHbpXECharged","#xi = ln(1/x_{E}) with charged hadrons",
             nptbins,ptmin,ptmax,200,0.,10.); 
    fhPtHbpXECharged->SetYTitle("ln(1/x_{E})");
    fhPtHbpXECharged->SetXTitle("p_{T trigger}");
    
    fhPtHbpXEUeCharged  = 
    new TH2F("hHbpXEUeCharged","#xi = ln(1/x_{E}) with charged hadrons,Underlying Event",
             nptbins,ptmin,ptmax,200,0.,10.); 
    fhPtHbpXEUeCharged->SetYTitle("ln(1/x_{E})");
    fhPtHbpXEUeCharged->SetXTitle("p_{T trigger}");
    
    fhZTCharged  = 
    new TH2F("hZTCharged","z_{T} for charged tracks",
             nptbins,ptmin,ptmax,200,0.,2.); 
    fhZTCharged->SetYTitle("z_{T}");
    fhZTCharged->SetXTitle("p_{T trigger}");
    
    fhZTUeCharged  = 
    new TH2F("hZTUeCharged","z_{T} for Underlying Event",
             nptbins,ptmin,ptmax,200,0.,2.); 
    fhZTUeCharged->SetYTitle("z_{T}");
    fhZTUeCharged->SetXTitle("p_{T trigger}");
    
    fhZTPosCharged  = 
    new TH2F("hZTPositiveCharged","z_{T} for positive charged tracks",
             nptbins,ptmin,ptmax,200,0.,2.); 
    fhZTPosCharged->SetYTitle("z_{T}");
    fhZTPosCharged->SetXTitle("p_{T trigger}");
    
    fhZTNegCharged  = 
    new TH2F("hZTNegativeCharged","z_{T} for negative charged tracks",
             nptbins,ptmin,ptmax,200,0.,2.); 
    fhZTNegCharged->SetYTitle("z_{T}");
    fhZTNegCharged->SetXTitle("p_{T trigger}");
    
    fhPtHbpZTCharged  = 
    new TH2F("hHbpZTCharged","#xi = ln(1/z_{T}) with charged hadrons",
             nptbins,ptmin,ptmax,200,0.,10.); 
    fhPtHbpZTCharged->SetYTitle("ln(1/z_{T})");
    fhPtHbpZTCharged->SetXTitle("p_{T trigger}");
    
    fhPtHbpZTUeCharged  = 
    new TH2F("hHbpZTUeCharged","#xi = ln(1/z_{T}) with charged hadrons,Underlying Event",
             nptbins,ptmin,ptmax,200,0.,10.); 
    fhPtHbpZTUeCharged->SetYTitle("ln(1/x_{E})");
    fhPtHbpZTUeCharged->SetXTitle("p_{T trigger}");
    
    fhPtTrigPout  = 
    new TH2F("hPtTrigPout","Pout with triggers",
             nptbins,ptmin,ptmax,2*nptbins,-ptmax,ptmax); 
    fhPtTrigPout->SetYTitle("p_{out} (GeV/c)");
    fhPtTrigPout->SetXTitle("p_{T trigger} (GeV/c)"); 
    
    fhPtTrigCharged  = 
    new TH2F("hPtTrigCharged","trgger and charged tracks pt distribution",
             nptbins,ptmin,ptmax,nptbins,ptmin,ptmax); 
    fhPtTrigCharged->SetYTitle("p_{T h^{#pm}} (GeV/c)");
    fhPtTrigCharged->SetXTitle("p_{T trigger} (GeV/c)");    
	  
    outputContainer->Add(fhDeltaPhiDeltaEtaCharged);
    outputContainer->Add(fhPhiCharged) ;
    outputContainer->Add(fhEtaCharged) ;
    outputContainer->Add(fhDeltaPhiCharged) ; 
    outputContainer->Add(fhDeltaEtaCharged) ;
    outputContainer->Add(fhDeltaPhiChargedPt) ;
    outputContainer->Add(fhDeltaPhiUeChargedPt) ;

    outputContainer->Add(fhXECharged) ;
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
        fhTrigDeltaPhiCharged[im]->SetYTitle("#Delta #phi");
        
        fhTrigDeltaEtaCharged[im]  = new TH2F 
        (Form("hTrigDeltaEtaCharged_%d",im),Form("hTrigDeltaEtaCharged_%d",im), nptbins,ptmin,ptmax, ndeltaetabins ,deltaetamin,deltaetamax); 
        fhTrigDeltaEtaCharged[im]->SetXTitle("p_{T trigger} (GeV/c)");
        fhTrigDeltaEtaCharged[im]->SetYTitle("#Delta #eta");
        
        fhTrigXECorr[im]  = new TH2F
        (Form("hTrigXEPtCorr_%d",im),Form("hTrigXEPtCorr_%d",im), nptbins,ptmin,ptmax,200,0.,2.); 
        fhTrigXECorr[im]->SetYTitle("x_{E trigger h^{#pm}}");
        fhTrigXECorr[im]->SetXTitle("p_{T trigger}");
        
        fhTrigXEUeCorr[im]  = new TH2F
        (Form("hTrigXEPtUeCorr_%d",im),Form("hTrigXEPtUeCorr_%d",im), nptbins,ptmin,ptmax,200,0.,2.); 
        fhTrigXEUeCorr[im]->SetYTitle("x_{E trigger h^{#pm}}");
        fhTrigXEUeCorr[im]->SetXTitle("p_{T trigger}");       
        
        fhTrigZTCorr[im]  = new TH2F
        (Form("hTrigZTPtCorr_%d",im),Form("hTrigZTPtCorr_%d",im), nptbins,ptmin,ptmax,200,0.,2.); 
        fhTrigZTCorr[im]->SetYTitle("z_{trigger h^{#pm}}");
        fhTrigZTCorr[im]->SetXTitle("p_{T trigger}");
        
        fhTrigZTUeCorr[im]  = new TH2F
        (Form("hTrigZTPtUeCorr_%d",im),Form("hTrigZTPtUeCorr_%d",im), nptbins,ptmin,ptmax,200,0.,2.); 
        fhTrigZTUeCorr[im]->SetYTitle("z_{trigger h^{#pm}}");
        fhTrigZTUeCorr[im]->SetXTitle("p_{T trigger}");               
        
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
      fhAssocPtBkg->SetXTitle("p_{T trigger}");
      fhAssocPtBkg->SetYTitle("p_{T associated}");
      outputContainer->Add(fhAssocPtBkg) ;
      
      fhDeltaPhiBrad = new TH2F("hDeltaPhiBrad","atan2(sin(#Delta #phi), cos(#Delta #phi))/#pi vs p_{T trigger} ", 
                                nptbins, ptmin, ptmax,288, -1.0/3.0, 5.0/3.0);
      fhDeltaPhiBrad->SetXTitle("p_{T trigger}");
      fhDeltaPhiBrad->SetYTitle("atan2(sin(#Delta #phi), cos(#Delta #phi))/#pi");
      outputContainer->Add(fhDeltaPhiBrad) ;
    }
    
    fhDeltaPhiAssocPtBin     = new TH2F*[fNAssocPtBins] ;
    fhXEAssocPtBin           = new TH2F*[fNAssocPtBins] ;
    fhZTAssocPtBin           = new TH2F*[fNAssocPtBins] ;
    
    if(fFillBradHisto)  
      fhDeltaPhiBradAssocPtBin = new TH2F*[fNAssocPtBins] ;
    
    if(fPi0Trigger || fDecayTrigger)
    {
      fhDeltaPhiAssocPtBin       = new TH2F*[fNAssocPtBins] ;
      fhXEAssocPtBin             = new TH2F*[fNAssocPtBins] ;
      fhZTAssocPtBin             = new TH2F*[fNAssocPtBins] ;
      fhXEDecayChargedAssocPtBin = new TH2F*[fNAssocPtBins] ;
      fhZTDecayChargedAssocPtBin = new TH2F*[fNAssocPtBins] ;
      fhDeltaPhiDecayChargedAssocPtBin = new TH2F*[fNAssocPtBins] ;
    }

    if(fHMPIDCorrelation)
    {
      fhDeltaPhiAssocPtBinHMPID   = new TH2F*[fNAssocPtBins] ;
      fhDeltaPhiAssocPtBinHMPIDAcc= new TH2F*[fNAssocPtBins] ;
    }
    
    for(Int_t i = 0 ; i < fNAssocPtBins ; i++)
    {
      fhDeltaPhiAssocPtBin[i] = new TH2F(Form("hDeltaPhiPtAssocPt%2.1f_%2.1f", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1]), 
                                         Form("#Delta #phi vs p_{T trigger} for associated p_{T} bin [%2.1f,%2.1f]", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1]), 
                                         nptbins, ptmin, ptmax, ndeltaphibins ,deltaphimin,deltaphimax);
      fhDeltaPhiAssocPtBin[i]->SetXTitle("p_{T trigger}");
      fhDeltaPhiAssocPtBin[i]->SetYTitle("#Delta #phi");
       
      fhXEAssocPtBin[i]       = new TH2F(Form("hXEAssocPtBin%1.f_%1.f", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1]), 
                                         Form("x_{E} vs p_{T trigger} for associated p_{T} bin [%2.1f,%2.1f]", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1]), 
                                         nptbins, ptmin, ptmax,200, 0.0, 2.0);
      fhXEAssocPtBin[i]->SetXTitle("p_{T trigger}");
      fhXEAssocPtBin[i]->SetYTitle("x_{E}");
      
      fhZTAssocPtBin[i]       = new TH2F(Form("hZTAssocPtBin%1.f_%1.f", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1]), 
                                         Form("z_{T} vs p_{T trigger} for associated p_{T} bin [%2.1f,%2.1f]", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1]), 
                                         nptbins, ptmin, ptmax,200, 0.0, 2.0);
      fhZTAssocPtBin[i]->SetXTitle("p_{T trigger}");
      fhZTAssocPtBin[i]->SetYTitle("z_{T}");
      
      outputContainer->Add(fhDeltaPhiAssocPtBin[i]) ;
      outputContainer->Add(fhXEAssocPtBin[i]);
      outputContainer->Add(fhZTAssocPtBin[i]);
      
      if(fPi0Trigger || fDecayTrigger) 
      {
        fhDeltaPhiDecayChargedAssocPtBin[i] = new TH2F(Form("hDeltaPhiPtDecayChargedAssocPt%2.1f_%2.1f", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1]), 
                                           Form("#Delta #phi vs p_{T trigger} tagged as decay for associated p_{T} bin [%2.1f,%2.1f]", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1]), 
                                           nptbins, ptmin, ptmax, ndeltaphibins ,deltaphimin,deltaphimax);
        fhDeltaPhiDecayChargedAssocPtBin[i]->SetXTitle("p_{T trigger}");
        fhDeltaPhiDecayChargedAssocPtBin[i]->SetYTitle("#Delta #phi");
        
        fhXEDecayChargedAssocPtBin[i]       = new TH2F(Form("hXEDecayChargedAssocPtBin%1.f_%1.f", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1]), 
                                           Form("x_{E} vs p_{T trigger} tagged as decay for associated p_{T} bin [%2.1f,%2.1f]", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1]), 
                                           nptbins, ptmin, ptmax,200, 0.0, 2.0);
        fhXEDecayChargedAssocPtBin[i]->SetXTitle("p_{T trigger}");
        fhXEDecayChargedAssocPtBin[i]->SetYTitle("x_{E}");
        
        fhZTDecayChargedAssocPtBin[i]       = new TH2F(Form("hZTDecayChargedAssocPtBin%1.f_%1.f", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1]), 
                                           Form("z_{T} vs p_{T trigger} tagged as decay for associated p_{T} bin [%2.1f,%2.1f]", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1]), 
                                           nptbins, ptmin, ptmax,200, 0.0, 2.0);
        fhZTDecayChargedAssocPtBin[i]->SetXTitle("p_{T trigger}");
        fhZTDecayChargedAssocPtBin[i]->SetYTitle("z_{T}");
        
        outputContainer->Add(fhDeltaPhiDecayChargedAssocPtBin[i]) ;
        outputContainer->Add(fhXEDecayChargedAssocPtBin[i]);
        outputContainer->Add(fhZTDecayChargedAssocPtBin[i]);
        
      }
      
      if(fFillBradHisto) 
      {
        fhDeltaPhiBradAssocPtBin[i] = new TH2F(Form("hDeltaPhiBradPtAssocPt%2.1f_%2.1f", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1]), 
                                               Form("atan2(sin(#Delta #phi), cos(#Delta #phi))/#pi vs p_{T trigger} for associated p_{T} bin [%2.1f,%2.1f]", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1]), 
                                               nptbins, ptmin, ptmax,288, -1.0/3.0, 5.0/3.0);
        fhDeltaPhiBradAssocPtBin[i]->SetXTitle("p_{T trigger}");
        fhDeltaPhiBradAssocPtBin[i]->SetYTitle("atan2(sin(#Delta #phi), cos(#Delta #phi))/#pi");
        outputContainer->Add(fhDeltaPhiBradAssocPtBin[i]) ;
      }       
      
      if(fHMPIDCorrelation)
      {
        fhDeltaPhiAssocPtBinHMPID[i] = new TH2F(Form("hDeltaPhiPtAssocPt%2.1f_%2.1fHMPID", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1]), 
                                                Form("#Delta #phi vs p_{T trigger} for associated p_{T} bin [%2.1f,%2.1f], with track having HMPID signal", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1]), 
                                                nptbins, ptmin, ptmax, ndeltaphibins ,deltaphimin,deltaphimax);
        fhDeltaPhiAssocPtBinHMPID[i]->SetXTitle("p_{T trigger}");
        fhDeltaPhiAssocPtBinHMPID[i]->SetYTitle("#Delta #phi");      
        
        fhDeltaPhiAssocPtBinHMPIDAcc[i] = new TH2F(Form("hDeltaPhiPtAssocPt%2.1f_%2.1fHMPIDAcc", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1]), 
                                                   Form("#Delta #phi vs p_{T trigger} for associated p_{T} bin [%2.1f,%2.1f], with track within 5<phi<20 deg", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1]), 
                                                   nptbins, ptmin, ptmax, ndeltaphibins ,deltaphimin,deltaphimax);
        fhDeltaPhiAssocPtBinHMPIDAcc[i]->SetXTitle("p_{T trigger}");
        fhDeltaPhiAssocPtBinHMPIDAcc[i]->SetYTitle("#Delta #phi"); 
        
        outputContainer->Add(fhDeltaPhiAssocPtBinHMPID[i]) ;
        outputContainer->Add(fhDeltaPhiAssocPtBinHMPIDAcc[i]) ;
        
      }      
    }
    
    if(fPi0Trigger || fDecayTrigger)
    {
      if(fPi0Trigger)
      {
        fhPtPi0DecayRatio  = new TH2F
        ("hPtPi0DecayRatio","p_T of #pi^{0} and the ratio of pt for two decay", 
         nptbins,ptmin,ptmax, 100,0.,2.); 
        fhPtPi0DecayRatio->SetXTitle("p_{T}^{#pi^{0}} (GeV/c)");
        fhPtPi0DecayRatio->SetYTitle("p_{T}^{Decay}/p_{T}^{#pi^{0}}");
        outputContainer->Add(fhPtPi0DecayRatio) ; 
      }
      
      fhDeltaPhiDecayCharged  = new TH2F
      ("hDeltaPhiDecayCharged","#phi_{Decay} - #phi_{h^{#pm}} vs p_{T Decay}",
       nptbins,ptmin,ptmax, ndeltaphibins ,deltaphimin,deltaphimax); 
      fhDeltaPhiDecayCharged->SetYTitle("#Delta #phi");
      fhDeltaPhiDecayCharged->SetXTitle("p_{T Decay} (GeV/c)");
      
      fhXEDecayCharged  = 
      new TH2F("hXEDecayCharged","x_{E}  Decay",
               nptbins,ptmin,ptmax,200,0.,2.); 
      fhXEDecayCharged->SetYTitle("x_{E}");
      fhXEDecayCharged->SetXTitle("p_{T decay}");
      
      fhZTDecayCharged  = 
      new TH2F("hZTDecayCharged","z_{trigger h^{#pm}} = p_{T h^{#pm}} / p_{T Decay}",
               nptbins,ptmin,ptmax,200,0.,2.); 
      fhZTDecayCharged->SetYTitle("z_{decay h^{#pm}}");
      fhZTDecayCharged->SetXTitle("p_{T decay}");      
      
      outputContainer->Add(fhDeltaPhiDecayCharged) ; 
      outputContainer->Add(fhXEDecayCharged) ;
      outputContainer->Add(fhZTDecayCharged) ;
    }    
    
    if(fMakeSeveralUE)
    { 
      fhDeltaPhiUeLeftCharged  = new TH2F
      ("hDeltaPhiUeLeftChargedPt","#phi_{trigger} - #phi_{#Ueh^{#pm}} vs p_{T Ueh^{#pm}} with UE left side range of trigger particles",
       nptbins,ptmin,ptmax, ndeltaphibins ,deltaphimin,deltaphimax);
      fhDeltaPhiUeLeftCharged->SetYTitle("#Delta #phi");
      fhDeltaPhiUeLeftCharged->SetXTitle("p_{T h^{#pm}} (GeV/c)");
      outputContainer->Add(fhDeltaPhiUeLeftCharged) ;
      
      fhDeltaPhiUeRightCharged  = new TH2F
      ("hDeltaPhiUeRightChargedPt","#phi_{trigger} - #phi_{#Ueh^{#pm}} vs p_{T Ueh^{#pm}} with UE right side range of trigger particles",
       nptbins,ptmin,ptmax, ndeltaphibins ,deltaphimin,deltaphimax);
      fhDeltaPhiUeRightCharged->SetYTitle("#Delta #phi");
      fhDeltaPhiUeRightCharged->SetXTitle("p_{T h^{#pm}} (GeV/c)");
      outputContainer->Add(fhDeltaPhiUeRightCharged) ;
      
      fhXEUeLeftCharged  = 
      new TH2F("hXEUeChargedLeft","x_{E} with UE left side of trigger",
               nptbins,ptmin,ptmax,200,0.,2.); 
      fhXEUeLeftCharged->SetYTitle("x_{E Ueh^{#pm}}");
      fhXEUeLeftCharged->SetXTitle("p_{T trigger}");
      outputContainer->Add(fhXEUeLeftCharged) ;
      
      fhXEUeRightCharged  = 
      new TH2F("hXEUeChargedRight","x_{E h^{#pm}} with UE right side of trigger",
               nptbins,ptmin,ptmax,200,0.,2.); 
      fhXEUeRightCharged->SetYTitle("z_{trigger Ueh^{#pm}}");
      fhXEUeRightCharged->SetXTitle("p_{T trigger}");
      outputContainer->Add(fhXEUeRightCharged) ;
      
      fhPtHbpXEUeLeftCharged  = 
      new TH2F("hHbpXEUeChargedLeft","#xi = ln(1/x_{E}) with charged UE left side of trigger",
               nptbins,ptmin,ptmax,200,0.,10.); 
      fhPtHbpXEUeLeftCharged->SetYTitle("ln(1/x_{E})");
      fhPtHbpXEUeLeftCharged->SetXTitle("p_{T trigger}");
      outputContainer->Add(fhPtHbpXEUeLeftCharged) ;
      
      fhPtHbpXEUeRightCharged  = 
      new TH2F("hHbpXEUeChargedRight","#xi = ln(1/x_{E}) with charged UE right side of trigger",
               nptbins,ptmin,ptmax,200,0.,10.); 
      fhPtHbpXEUeRightCharged->SetYTitle("ln(1/x_{E})");
      fhPtHbpXEUeRightCharged->SetXTitle("p_{T trigger}");
      outputContainer->Add(fhPtHbpXEUeRightCharged) ;
      
      fhZTUeLeftCharged  = 
      new TH2F("hZTUeChargedLeft","z_{trigger h^{#pm}} = p_{T Ueh^{#pm}} / p_{T trigger} with UE left side of trigger",
               nptbins,ptmin,ptmax,200,0.,2.); 
      fhZTUeLeftCharged->SetYTitle("z_{trigger Ueh^{#pm}}");
      fhZTUeLeftCharged->SetXTitle("p_{T trigger}");
      outputContainer->Add(fhZTUeLeftCharged) ;
      
      fhZTUeRightCharged  = 
      new TH2F("hZTUeChargedRight","z_{trigger h^{#pm}} = p_{T Ueh^{#pm}} / p_{T trigger} with UE right side of trigger",
               nptbins,ptmin,ptmax,200,0.,2.); 
      fhZTUeRightCharged->SetYTitle("z_{trigger Ueh^{#pm}}");
      fhZTUeRightCharged->SetXTitle("p_{T trigger}");
      outputContainer->Add(fhZTUeRightCharged) ;      
      
      fhPtHbpZTUeLeftCharged  = 
      new TH2F("hHbpZTUeChargedLeft","#xi = ln(1/z_{T}) with charged UE left side of trigger",
               nptbins,ptmin,ptmax,200,0.,10.); 
      fhPtHbpZTUeLeftCharged->SetYTitle("ln(1/z_{T})");
      fhPtHbpZTUeLeftCharged->SetXTitle("p_{T trigger}");
      outputContainer->Add(fhPtHbpZTUeLeftCharged) ;
      
      fhPtHbpZTUeRightCharged  = 
      new TH2F("hHbpZTUeChargedRight","#xi = ln(1/z_{T}) with charged UE right side of trigger",
               nptbins,ptmin,ptmax,200,0.,10.); 
      fhPtHbpZTUeRightCharged->SetYTitle("ln(1/z_{T})");
      fhPtHbpZTUeRightCharged->SetXTitle("p_{T trigger}");
      outputContainer->Add(fhPtHbpZTUeRightCharged) ;
      
    }  
  }  //Correlation with charged hadrons

  //Correlation with neutral hadrons
  if(fNeutralCorr)
  {
    fhDeltaPhiDeltaEtaNeutral  = new TH2F
    ("hDeltaPhiDeltaEtaNeutral","#phi_{trigger} - #phi_{h^{0}} vs #eta_{trigger} - #eta_{h^{0}}",
     ndeltaphibins ,deltaphimin,deltaphimax, ndeltaetabins ,deltaetamin,deltaetamax); 
    fhDeltaPhiDeltaEtaNeutral->SetXTitle("#Delta #phi");
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
    fhDeltaPhiNeutral->SetYTitle("#Delta #phi");
    fhDeltaPhiNeutral->SetXTitle("p_{T trigger} (GeV/c)");
    
    fhDeltaPhiNeutralPt  = new TH2F
    ("hDeltaPhiNeutralPt","#phi_{trigger} - #phi_{#pi^{0}} vs p_{T #pi^{0}}}",
     nptbins,ptmin,ptmax, ndeltaphibins ,deltaphimin,deltaphimax); 
    fhDeltaPhiNeutralPt->SetYTitle("#Delta #phi");
    fhDeltaPhiNeutralPt->SetXTitle("p_{T h^{0}} (GeV/c)");
    
    fhDeltaPhiUeNeutralPt  = new TH2F
    ("hDeltaPhiUeNeutralPt","#phi_{trigger} - #phi_{#pi^{0}} vs p_{T #pi^{0}}}",
     nptbins,ptmin,ptmax, ndeltaphibins ,deltaphimin,deltaphimax); 
    fhDeltaPhiUeNeutralPt->SetYTitle("#Delta #phi");
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
    fhXENeutral->SetXTitle("p_{T trigger}");
    
    fhXEUeNeutral  = 
    new TH2F("hXEUeNeutral","x_{E} for #pi^{0} associated",
             nptbins,ptmin,ptmax,200,0.,2.); 
    fhXEUeNeutral->SetYTitle("x_{E}");
    fhXEUeNeutral->SetXTitle("p_{T trigger}");
    
    fhPtHbpXENeutral  = 
    new TH2F("hHbpXENeutral","#xi = ln(1/x_{E})for #pi^{0} associated",
             nptbins,ptmin,ptmax,200,0.,10.); 
    fhPtHbpXENeutral->SetYTitle("ln(1/x_{E})");
    fhPtHbpXENeutral->SetXTitle("p_{T trigger}");
    
    fhPtHbpXEUeNeutral  = 
    new TH2F("hHbpXEUeNeutral","#xi = ln(1/x_{E}) for #pi^{0} associated",
             nptbins,ptmin,ptmax,200,0.,10.); 
    fhPtHbpXEUeNeutral->SetYTitle("ln(1/x_{E})");
    fhPtHbpXEUeNeutral->SetXTitle("p_{T trigger}");
    
    fhZTNeutral  = 
    new TH2F("hZTNeutral","z_{trigger #pi} = p_{T #pi^{0}} / p_{T trigger} for #pi^{0} associated",
             nptbins,ptmin,ptmax,200,0.,2.); 
    fhZTNeutral->SetYTitle("z_{trigger #pi^{0}}");
    fhZTNeutral->SetXTitle("p_{T trigger}");
    
    fhZTUeNeutral  = 
    new TH2F("hZTUeNeutral","z_{trigger #pi} = p_{T #pi^{0}} / p_{T trigger} for #pi^{0} associated",
             nptbins,ptmin,ptmax,200,0.,2.); 
    fhZTUeNeutral->SetYTitle("z_{trigger #pi^{0}}");
    fhZTUeNeutral->SetXTitle("p_{T trigger}");
    
    fhPtHbpZTNeutral  = 
    new TH2F("hHbpZTNeutral","#xi = ln(1/x_{E}) for #pi^{0} associated",
             nptbins,ptmin,ptmax,200,0.,10.); 
    fhPtHbpZTNeutral->SetYTitle("ln(1/z_{T})");
    fhPtHbpZTNeutral->SetXTitle("p_{T trigger}");
    
    fhPtHbpZTUeNeutral  = 
    new TH2F("hHbpZTUeNeutral","#xi = ln(1/x_{E}) for #pi^{0} associated",
             nptbins,ptmin,ptmax,200,0.,10.); 
    fhPtHbpXEUeNeutral->SetYTitle("ln(1/z_{T})");
    fhPtHbpXEUeNeutral->SetXTitle("p_{T trigger}");    
    
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
      fhDeltaPhiDecayNeutral->SetYTitle("#Delta #phi");
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
      fhDeltaPhiUeLeftNeutral->SetYTitle("#Delta #phi");
      fhDeltaPhiUeLeftNeutral->SetXTitle("p_{T h^{0}} (GeV/c)");
      outputContainer->Add(fhDeltaPhiUeLeftNeutral) ;
      
      fhDeltaPhiUeRightNeutral  = new TH2F
      ("hDeltaPhiUeRightNeutralPt","#phi_{trigger} - #phi_{#Ueh^{0}} vs p_{T Ueh^{0}} with neutral UE right side range of trigger particles",
       nptbins,ptmin,ptmax, ndeltaphibins ,deltaphimin,deltaphimax); 
      fhDeltaPhiUeRightNeutral->SetYTitle("#Delta #phi");
      fhDeltaPhiUeRightNeutral->SetXTitle("p_{T h^{0}} (GeV/c)");
      outputContainer->Add(fhDeltaPhiUeRightNeutral) ;
      
      fhXEUeLeftNeutral  = 
      new TH2F("hXEUeNeutralLeft","x_{E} = p_{T Ueh^{0}} / p_{T trigger} with neutral UE left side of trigger",
               nptbins,ptmin,ptmax,140,0.,2.); 
      fhXEUeLeftNeutral->SetYTitle("z_{trigger Ueh^{0}}");
      fhXEUeLeftNeutral->SetXTitle("p_{T trigger}");
      outputContainer->Add(fhXEUeLeftNeutral) ;
      
      fhXEUeRightNeutral  = 
      new TH2F("hXEUeNeutralRight","x_{E} = p_{T Ueh^{0}} / p_{T trigger} with neutral UE right side of trigger",
               nptbins,ptmin,ptmax,200,0.,2.); 
      fhXEUeRightNeutral->SetYTitle("z_{trigger Ueh^{0}}");
      fhXEUeRightNeutral->SetXTitle("p_{T trigger}");
      outputContainer->Add(fhXEUeRightNeutral) ;
      
      fhPtHbpXEUeLeftNeutral  = 
      new TH2F("hHbpXEUeNeutralLeft","#xi = ln(1/x_{E}) with neutral UE left side of trigger",
               nptbins,ptmin,ptmax,200,0.,10.); 
      fhPtHbpXEUeLeftNeutral->SetYTitle("ln(1/x_{E})");
      fhPtHbpXEUeLeftNeutral->SetXTitle("p_{T trigger}");
      outputContainer->Add(fhPtHbpXEUeLeftNeutral) ;
      
      fhPtHbpXEUeRightNeutral  = 
      new TH2F("hHbpXEUeNeutralRight","#xi = ln(1/x_{E}) with neutral UE right side of trigger",
               nptbins,ptmin,ptmax,200,0.,10.); 
      fhPtHbpXEUeRightNeutral->SetYTitle("ln(1/x_{E})");
      fhPtHbpXEUeRightNeutral->SetXTitle("p_{T trigger}");
      outputContainer->Add(fhPtHbpXEUeRightNeutral) ;
      
      fhZTUeLeftNeutral  = 
      new TH2F("hZTUeNeutralLeft","z_{trigger h^{0}} = p_{T Ueh^{0}} / p_{T trigger} with neutral UE left side of trigger",
               nptbins,ptmin,ptmax,140,0.,2.); 
      fhZTUeLeftNeutral->SetYTitle("z_{trigger Ueh^{0}}");
      fhZTUeLeftNeutral->SetXTitle("p_{T trigger}");
      outputContainer->Add(fhZTUeLeftNeutral) ;
      
      fhZTUeRightNeutral  = 
      new TH2F("hZTUeNeutralRight","z_{trigger h^{0}} = p_{T Ueh^{0}} / p_{T trigger} with neutral UE right side of trigger",
               nptbins,ptmin,ptmax,200,0.,2.); 
      fhZTUeRightNeutral->SetYTitle("z_{trigger Ueh^{0}}");
      fhZTUeRightNeutral->SetXTitle("p_{T trigger}");
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
    fhMCDeltaPhiDeltaEtaCharged->SetXTitle("#Delta #phi");
    fhMCDeltaPhiDeltaEtaCharged->SetYTitle("#Delta #eta");    
    
    fhMCDeltaEtaCharged  = new TH2F
    ("hMCDeltaEtaCharged","MC #eta_{trigger} - #eta_{h^{#pm}} vs p_{T trigger} and p_{T assoc}",
     nptbins,ptmin,ptmax,200,-2,2); 
    fhMCDeltaEtaCharged->SetYTitle("#Delta #eta");
    fhMCDeltaEtaCharged->SetXTitle("p_{T trigger} (GeV/c)");
    
    fhMCDeltaPhiCharged  = new TH2F
    ("hMCDeltaPhiCharged","#phi_{trigger} - #phi_{h^{#pm}} vs p_{T trigger}",
     nptbins,ptmin,ptmax,ndeltaphibins ,deltaphimin,deltaphimax); 
    fhMCDeltaPhiCharged->SetYTitle("#Delta #phi");
    fhMCDeltaPhiCharged->SetXTitle("p_{T trigger} (GeV/c)");
    
    fhMCDeltaPhiChargedPt  = new TH2F
    ("hMCDeltaPhiChargedPt","MC #phi_{trigger} - #phi_{#h^{#pm}} vs p_{T h^{#pm}}",
     nptbins,ptmin,ptmax,ndeltaphibins ,deltaphimin,deltaphimax); 
    fhMCDeltaPhiChargedPt->SetYTitle("#Delta #phi");
    fhMCDeltaPhiChargedPt->SetXTitle("p_{T h^{#pm}} (GeV/c)");
    
    fhMCPtXECharged  = 
    new TH2F("hMCPtXECharged","x_{E}",
             nptbins,ptmin,ptmax,200,0.,2.); 
    fhMCPtXECharged->SetYTitle("x_{E}");
    fhMCPtXECharged->SetXTitle("p_{T trigger}");  
    
    fhMCPtHbpXECharged  = 
    new TH2F("hMCHbpXECharged","MC #xi = ln(1/x_{E}) with charged hadrons",
             nptbins,ptmin,ptmax,200,0.,10.); 
    fhMCPtHbpXECharged->SetYTitle("ln(1/x_{E})");
    fhMCPtHbpXECharged->SetXTitle("p_{T trigger}");
    
    fhMCPtZTCharged  = 
    new TH2F("hMCPtZTCharged","z_{T}",
             nptbins,ptmin,ptmax,200,0.,2.); 
    fhMCPtZTCharged->SetYTitle("z_{T}");
    fhMCPtZTCharged->SetXTitle("p_{T trigger}"); 
    
    fhMCPtHbpZTCharged  = 
    new TH2F("hMCHbpZTCharged","MC #xi = ln(1/z_{T}) with charged hadrons",
             nptbins,ptmin,ptmax,200,0.,10.); 
    fhMCPtHbpZTCharged->SetYTitle("ln(1/z_{T})");
    fhMCPtHbpZTCharged->SetXTitle("p_{T trigger}");
    
    fhMCPtTrigPout  = 
    new TH2F("hMCPtTrigPout","AOD MC Pout with triggers",
             nptbins,ptmin,ptmax,2*nptbins,-ptmax,ptmax); 
    fhMCPtTrigPout->SetYTitle("p_{out} (GeV/c)");
    fhMCPtTrigPout->SetXTitle("p_{T trigger} (GeV/c)"); 
    
    fhMCPtAssocDeltaPhi  = 
    new TH2F("hMCPtAssocDeltaPhi","AOD MC delta phi with associated charged hadrons",
             nptbins,ptmin,ptmax,ndeltaphibins ,deltaphimin,deltaphimax); 
    fhMCPtAssocDeltaPhi->SetYTitle("#Delta #phi");
    fhMCPtAssocDeltaPhi->SetXTitle("p_{T trigger} (GeV/c)"); 
        
    outputContainer->Add(fh2phiLeadingParticle);
    outputContainer->Add(fhMCDeltaPhiDeltaEtaCharged);
    outputContainer->Add(fhMCPhiCharged) ;
    outputContainer->Add(fhMCEtaCharged) ;
    outputContainer->Add(fhMCDeltaEtaCharged) ;
    outputContainer->Add(fhMCDeltaPhiCharged) ; 
    
    outputContainer->Add(fhMCDeltaPhiChargedPt) ;
    outputContainer->Add(fhMCPtXECharged) ;
    outputContainer->Add(fhMCPtZTCharged) ;
    outputContainer->Add(fhMCPtHbpXECharged) ;
    outputContainer->Add(fhMCPtHbpZTCharged) ;
    outputContainer->Add(fhMCPtTrigPout) ;
    outputContainer->Add(fhMCPtAssocDeltaPhi) ;      
  } //for MC histogram
  
  if(fDoOwnMix)
  {
    //create event containers
    fListMixEvents= new TList*[GetNCentrBin()*GetNZvertBin()] ;
    
    for(Int_t ic=0; ic<GetNCentrBin(); ic++){
      for(Int_t iz=0; iz<GetNZvertBin(); iz++){
        fListMixEvents[ic*GetNZvertBin()+iz] = new TList() ;
        fListMixEvents[ic*GetNZvertBin()+iz]->SetOwner(kFALSE);
      }
    }    
    
    fhNtracksAll=new TH1F("hNtracksAll","Number of tracks w/o event trigger",2000,0,2000);
    outputContainer->Add(fhNtracksAll);
    
    fhNtracksTrigger=new TH1F("hNtracksTriggerEvent","Number of tracks w/ event trigger",2000,0,2000);
    outputContainer->Add(fhNtracksTrigger);
    
    fhNtracksINT=new TH1F("hNtracksMBEvent","Number of tracks w/ event trigger kAnyINT",2000,0,2000);
    outputContainer->Add(fhNtracksINT);
    
    fhMixDeltaPhiCharged  = new TH2F
    ("hMixDeltaPhiCharged","Mixed event : #phi_{trigger} - #phi_{h^{#pm}} vs p_{T trigger}",
     nptbins,ptmin,ptmax,ndeltaphibins ,deltaphimin,deltaphimax); 
    fhMixDeltaPhiCharged->SetYTitle("#Delta #phi");
    fhMixDeltaPhiCharged->SetXTitle("p_{T trigger} (GeV/c)");
    outputContainer->Add(fhMixDeltaPhiCharged);

    fhMixDeltaPhiDeltaEtaCharged  = new TH2F
    ("hMixDeltaPhiDeltaEtaCharged","Mixed event : #phi_{trigger} - #phi_{h^{#pm}} vs #eta_{trigger} - #eta_{h^{#pm}}",
     ndeltaphibins ,deltaphimin,deltaphimax,ndeltaetabins ,deltaetamin,deltaetamax); 
    fhMixDeltaPhiDeltaEtaCharged->SetXTitle("#Delta #phi");
    fhMixDeltaPhiDeltaEtaCharged->SetYTitle("#Delta #eta");
    outputContainer->Add(fhMixDeltaPhiDeltaEtaCharged);

    fhMixDeltaPhiChargedAssocPtBin         = new TH2F*[fNAssocPtBins] ;
    fhMixDeltaPhiDeltaEtaChargedAssocPtBin = new TH2F*[fNAssocPtBins] ;

    for(Int_t i = 0 ; i < fNAssocPtBins ; i++)
    {
      fhMixDeltaPhiChargedAssocPtBin[i] = new TH2F(Form("hMixDeltaPhiChargedAssocPtBin%2.1f_%2.1f", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1]), 
                                         Form("Mixed event #Delta #phi vs p_{T trigger} for associated p_{T} bin [%2.1f,%2.1f]", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1]), 
                                         nptbins, ptmin, ptmax,  ndeltaphibins ,deltaphimin,deltaphimax);
      fhMixDeltaPhiChargedAssocPtBin[i]->SetXTitle("p_{T trigger}");
      fhMixDeltaPhiChargedAssocPtBin[i]->SetYTitle("#Delta #phi");

      fhMixDeltaPhiDeltaEtaChargedAssocPtBin[i] = new TH2F(Form("hMixDeltaPhiDeltaEtaChargedAssocPtBin%2.1f_%2.1f", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1]), 
                                                   Form("Mixed event #Delta #phi vs p_{T trigger} for associated p_{T} bin [%2.1f,%2.1f]", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1]), 
                                                    ndeltaphibins ,deltaphimin,deltaphimax,ndeltaetabins ,deltaetamin,deltaetamax); 
      fhMixDeltaPhiDeltaEtaChargedAssocPtBin[i]->SetXTitle("#Delta #phi");
      fhMixDeltaPhiDeltaEtaChargedAssocPtBin[i]->SetYTitle("#Delta #eta");
      
      outputContainer->Add(fhMixDeltaPhiChargedAssocPtBin[i]);
      outputContainer->Add(fhMixDeltaPhiDeltaEtaChargedAssocPtBin[i]);
      
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


//____________________________________________________
void AliAnaParticleHadronCorrelation::InitParameters()
{
  
  //Initialize the parameters of the analysis.
  SetInputAODName("PWG4Particle");
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
  fAssocPtBinLimit[1]   = 2.0 ; 
  fAssocPtBinLimit[2]   = 4.0 ; 
  fAssocPtBinLimit[3]   = 6.0 ; 
  fAssocPtBinLimit[4]   = 8.0 ; 
  fAssocPtBinLimit[5]   = 10. ; 
  fAssocPtBinLimit[6]   = 12. ;
  fAssocPtBinLimit[7]   = 15. ;
  fAssocPtBinLimit[8]   = 25. ;
  fAssocPtBinLimit[9]   = 50. ;
  
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
	
  if(GetDebug() > 1){
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
  
  // Fill the pool with tracks if requested
  if(fDoOwnMix) FillChargedEventMixPool();

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

    // Check if trigger is in fiducial region
    if(IsFiducialCutOn())
    {
      Bool_t in = GetFiducialCut()->IsInFiducialCut(*particle->Momentum(),particle->GetDetector()) ;
      if(! in ) return ;
    }    
    
    // Check if the particle is isolated or if we want to take the isolation into account
    if(OnlyIsolated() && !particle->IsIsolated()) return;
    
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
      fhPtLeading->Fill(particle->Pt());
      Float_t phi = particle->Phi();
      if(phi<0)phi+=TMath::TwoPi();
      fhPhiLeading->Fill(particle->Pt(), phi);
      fhEtaLeading->Fill(particle->Pt(), particle->Eta());
    }//ok charged && neutral
  }//Aod branch loop
  
  //Reinit for next event
  fLeadingTriggerIndex = -1;
  
  if(GetDebug() > 1) printf("AliAnaParticleHadronCorrelation::MakeAnalysisFillHistograms() - End fill histograms \n");
}

//___________________________________________________________________________________________________________
Bool_t  AliAnaParticleHadronCorrelation::MakeChargedCorrelation(AliAODPWG4ParticleCorrelation *aodParticle, 
                                                                const TObjArray* pl, const Bool_t bFillHisto)
{  
  // Charged Hadron Correlation Analysis
  if(GetDebug() > 1) printf("AliAnaParticleHadronCorrelation::MakeChargedCorrelation() - Make trigger particle - charged hadron correlation \n");
    
  Float_t phiTrig = aodParticle->Phi();
  Float_t etaTrig = aodParticle->Eta();
  Float_t ptTrig  = aodParticle->Pt();  
  Bool_t   decay  = aodParticle->IsTagged();

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

      if(GetDebug() > 2 ) printf("AliAnaParticleHadronCorrelation::MakeChargedCorrelation() - Selected charge for momentum imbalance: pt %2.2f, phi %2.2f, eta %2.2f \n",pt,phi,eta);
            
      // Set the pt associated bin for the defined bins
      Int_t assocBin   = -1; 
      for(Int_t i = 0 ; i < fNAssocPtBins ; i++)
      {
        if(pt > fAssocPtBinLimit[i] && pt < fAssocPtBinLimit[i+1]) assocBin= i; 
      }      
      
      // Azimuthal Angle
      // calculate deltaPhi for later, shift when needed
      FillChargedAngularCorrelationHistograms(pt,  ptTrig,  assocBin, phi, phiTrig,  deltaPhi,
                                              eta, etaTrig, decay, track->GetHMPIDsignal(),nTracks);
      
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
      if      ( (deltaPhi > fDeltaPhiMinCut)   && (deltaPhi < fDeltaPhiMaxCut)   ) 
      {
        
        FillChargedMomentumImbalanceHistograms(ptTrig, pt, xE, hbpXE, zT, hbpZT, pout, deltaPhi, 
                                               nTracks, track->Charge(), assocBin, decay);
        
      } 
      else if ( (deltaPhi > fUeDeltaPhiMinCut) && (deltaPhi < fUeDeltaPhiMaxCut) ) 
      { //UE study
        
        FillChargedUnderlyingEventHistograms(ptTrig, pt, deltaPhi, nTracks);
        
      }
      
      if(fPi0Trigger && decayFound) 
        FillDecayPhotonCorrelationHistograms(pt, phi, decayMom1,decayMom2, kTRUE) ;
      
      //several UE calculation 
      if(fMakeSeveralUE) FillChargedUnderlyingEventSidesHistograms(ptTrig,pt,xE,hbpXE,zT,hbpZT,deltaPhi);
      
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
  if(fDoOwnMix && bFillHisto)
  {
      MakeChargedMixCorrelation(aodParticle);
  }
  
  return kTRUE;
  
}  


//______________________________________________________________________________________________________________________
void AliAnaParticleHadronCorrelation::MakeChargedMixCorrelation(AliAODPWG4ParticleCorrelation *aodParticle)
{  
  if(GetDebug() > 1) printf("AliAnaParticleHadronCorrelationNew::MakeChargedMixCorrelation() - Make trigger particle - charged hadron mixed event correlation \n");
  
  Double_t v[3] = {0,0,0}; //vertex ;
  GetReader()->GetVertex(v);
  
  if(GetMixedEvent()) return;  // This is not the mixed event from general mixing frame
  
  
  // Get the event with similar caracteristics
  
  //Get vertex z bin
  
  Int_t curZvertBin = (Int_t)(0.5*GetNZvertBin()*(v[2]+GetZvertexCut())/GetZvertexCut()) ;
  
  // centrality or tracks bin
  Int_t curCentrBin = 0;
  if(fUseTrackMultBins)
  { // Track multiplicity bins
    //curCentrBin = (GetTrackMultiplicity()-1)/5; 
    //if(curCentrBin > GetNCentrBin()-1) curCentrBin=GetNCentrBin()-1;
    Int_t trackMult = GetReader()->GetTrackMultiplicity();
    if(trackMult<=5)
      curCentrBin=8;
    else if(trackMult<=10)
      curCentrBin=7;
    else if(trackMult<=15)
      curCentrBin=6;
    else if(trackMult<=20)
      curCentrBin=5;
    else if(trackMult<=30)
      curCentrBin=4;
    else if(trackMult<=40)
      curCentrBin=3;
    else if(trackMult<=55)
      curCentrBin=2;
    else if(trackMult<=70)
      curCentrBin=1 ;
    else curCentrBin=0 ;        
  }
  else // Set centrality based on centrality task
  {
    curCentrBin = GetEventCentrality() * GetNCentrBin() / GetReader()->GetCentralityOpt(); 
    if(GetDebug() > 0 )printf("AliAnaParticleHadronCorrelation::MakeAnalysisFillHistograms() - curCentrBin %d, centrality %d, n bins %d, max bin from centrality %d\n",
                              curCentrBin, GetEventCentrality(), GetNCentrBin(), GetReader()->GetCentralityOpt());        
  }
  
  // Get the pool, check if it exits
  if(!fListMixEvents[curCentrBin*GetNZvertBin()+curZvertBin]) return ;
  
  //printf("curCen %d, curZ %d, bin %d\n",curCentrBin,curZvertBin,curCentrBin*GetNZvertBin()+curZvertBin);
  
  TList * pool = fListMixEvents[curCentrBin*GetNZvertBin()+curZvertBin];
    
  Double_t ptTrig  = aodParticle->Pt();
  Double_t etaTrig = aodParticle->Eta();
  Double_t phiTrig = aodParticle->Phi();
  if(phiTrig < 0.) phiTrig+=TMath::TwoPi();
  
  if(GetDebug() > 1) printf("AliAnaParticleHadronCorrelationNew::MakeChargedMixCorrelation() - leading trigger pt=%f, phi=%f, eta=%f\n",ptTrig,phiTrig,etaTrig);
  
  Double_t ptAssoc  = -999.;
  Double_t phiAssoc = -999. ;
  Double_t etaAssoc = -999.;
  Double_t deltaPhi = -999.;
  Double_t deltaEta = -999.;
  
  for(Int_t ev=0; ev <pool->GetSize(); ev++)
  {
    TObjArray* bgTracks = static_cast<TObjArray*>(pool->At(ev));
    
    Int_t nTracks=bgTracks->GetEntriesFast();
    
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

      fhMixDeltaPhiCharged        ->Fill(ptTrig,  deltaPhi);
      fhMixDeltaPhiDeltaEtaCharged->Fill(deltaPhi, deltaEta);

      if(assocBin < 0) continue ; // this pt bin was not considered
      
      fhMixDeltaPhiChargedAssocPtBin[assocBin]        ->Fill(ptTrig, deltaPhi);
      fhMixDeltaPhiDeltaEtaChargedAssocPtBin[assocBin]->Fill(deltaPhi, deltaEta);
      
    } // track loop
  } // mixed event loop
}
  

//________________________________________________________________________________________________________________
Bool_t  AliAnaParticleHadronCorrelation::MakeNeutralCorrelation(AliAODPWG4ParticleCorrelation * const aodParticle, 
                                                                const TObjArray* pi0list, const Bool_t bFillHisto)  
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
      if(fMakeSeveralUE) FillChargedUnderlyingEventSidesHistograms(ptTrig,pt,xE,hbpXE,zT,hbpZT,deltaPhi);

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
    if(!stack) {
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
    mcparticles0 = GetReader()->GetAODMCParticles(0);
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
      
      if(ptprim < 0.01 || eprim < 0.01) return ;
      
      mcparticles= GetReader()->GetAODMCParticles();
      for (Int_t i = 0; i < nTracks; i++) 
      {
        AliAODMCParticle *part = (AliAODMCParticle*) mcparticles->At(i);
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
              }
            } // in acceptance
          } // min pt cut
        } //only charged particles
      }  //MC particle loop      
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
  
  
  if(n < 10 && n > 0)
  {
    fNAssocPtBins  = n ; 
  }
  else 
  {
    printf("n = larger than 9 or too small, set to 9 \n");
    fNAssocPtBins = 9;
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
  else {
    printf("AliAnaParticleHadronCorrelation::SetAssocPtBinLimit() - bin  number too large %d > %d or small, nothing done\n", ibin, fNAssocPtBins) ; 
    
  }
}

