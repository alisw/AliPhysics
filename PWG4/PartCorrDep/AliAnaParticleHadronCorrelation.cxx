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
/* $Id: $ */

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
//////////////////////////////////////////////////////////////////////////////


// --- ROOT system ---
//#include "TClonesArray.h"
#include "TClass.h"
#include "TMath.h"
#include "TH3D.h"

//---- ANALYSIS system ----
#include "AliNeutralMesonSelection.h" 
#include "AliAnaParticleHadronCorrelation.h" 
#include "AliCaloTrackReader.h"
#include "AliCaloPID.h"
#include "AliAODPWG4ParticleCorrelation.h"
#include "AliFiducialCut.h"
#include "AliVTrack.h"
#include "AliVCluster.h"
#include "AliMCAnalysisUtils.h"
#include "TParticle.h"
#include "AliStack.h"
#include "AliAODMCParticle.h"
#include "AliMixedEvent.h"

ClassImp(AliAnaParticleHadronCorrelation)


//____________________________________________________________________________
  AliAnaParticleHadronCorrelation::AliAnaParticleHadronCorrelation(): 
    AliAnaPartCorrBaseClass(),
    fDeltaPhiMaxCut(0.), fDeltaPhiMinCut(0.), fSelectIsolated(0),
    fMakeSeveralUE(0),  
    fUeDeltaPhiMaxCut(0.), fUeDeltaPhiMinCut(0.), 
    fPi0AODBranchName(""),fNeutralCorr(0), fPi0Trigger(0),
   // fMultiBin(0),fNZvertBin(0),fNrpBin(0),fZvtxCut(0.),
   // fUseSelectEvent(kFALSE), 
   // fhNclustersNtracks(0), //fhVertex(0),
    fhPtLeading(0),fhPhiLeading(0),fhEtaLeading(0), 
    fhDeltaPhiDeltaEtaCharged(0),
    fhPhiCharged(0), fhEtaCharged(0), 
    fhDeltaPhiCharged(0), 
    fhDeltaEtaCharged(0), 
    fhDeltaPhiChargedPt(0), 
    fhDeltaPhiUeChargedPt(0), 
    fhPtImbalanceCharged(0), 
    fhPtImbalanceUeCharged(0),
    fhPtImbalancePosCharged(0),fhPtImbalanceNegCharged(0),
    fhPtHbpCharged(0), fhPtHbpUeCharged(0),
    fhDeltaPhiUeLeftCharged(0),fhDeltaPhiUeRightCharged(0),
    fhPtImbalanceUeLeftCharged(0),fhPtImbalanceUeRightCharged(0),
    fhPtHbpUeLeftCharged(0),fhPtHbpUeRightCharged(0), 
    fhPtTrigPout(0), fhPtAssocDeltaPhi(0),
    fhPtTrigCharged(0),
    fhTrigDeltaPhiCharged(0x0), fhTrigDeltaEtaCharged(0x0),fhTrigCorr(0x0),fhTrigUeCorr(0x0),
    fhDeltaPhiDeltaEtaNeutral(0), 
    fhPhiNeutral(0), fhEtaNeutral(0), 
    fhDeltaPhiNeutral(0), fhDeltaEtaNeutral(0),
    fhDeltaPhiNeutralPt(0),fhDeltaPhiUeNeutralPt(0), 
    fhPtImbalanceNeutral(0),fhPtImbalanceUeNeutral(0),
    fhPtHbpNeutral(0), fhPtHbpUeNeutral(0),
    fhDeltaPhiUeLeftNeutral(0),fhDeltaPhiUeRightNeutral(0),
    fhPtImbalanceUeLeftNeutral(0),fhPtImbalanceUeRightNeutral(0),
    fhPtHbpUeLeftNeutral(0),fhPtHbpUeRightNeutral(0),
    fhPtPi0DecayRatio(0),
    fhDeltaPhiDecayCharged(0),
    fhPtImbalanceDecayCharged(0), 
    fhDeltaPhiDecayNeutral(0),
    fhPtImbalanceDecayNeutral(0)
{
  //Default Ctor
  
  //Initialize parameters
  InitParameters();
}

//________________________________________________________________________
TList *  AliAnaParticleHadronCorrelation::GetCreateOutputObjects()
{  
  
  // Create histograms to be saved in output file and 
  // store them in fOutputContainer
  TList * outputContainer = new TList() ; 
  outputContainer->SetName("CorrelationHistos") ; 
  
  Int_t nptbins  = GetHistoPtBins();
  Int_t nphibins = GetHistoPhiBins();
  Int_t netabins = GetHistoEtaBins();
  Float_t ptmax  = GetHistoPtMax();
  Float_t phimax = GetHistoPhiMax();
  Float_t etamax = GetHistoEtaMax();
  Float_t ptmin  = GetHistoPtMin();
  Float_t phimin = GetHistoPhiMin();
  Float_t etamin = GetHistoEtaMin();	

  
//  fhNclustersNtracks  = new TH2F ("Multiplicity","Neutral cluster and charged track multiplicity",1000, 0., 1000.,1000, 0., 1000.); 
//  fhNclustersNtracks->SetYTitle("# of tracks");
//  fhNclustersNtracks->SetXTitle("# of clusters");

  fhPtLeading  = new TH1F ("hPtLeading","p_T distribution of leading particles", nptbins,ptmin,ptmax); 
  fhPtLeading->SetXTitle("p_{T}^{trig} (GeV/c)");
  
  fhPhiLeading  = new TH2F ("hPhiLeading","#phi distribution of leading Particles",nptbins,ptmin,ptmax, nphibins,phimin,phimax); 
  fhPhiLeading->SetYTitle("#phi (rad)");
  
  fhEtaLeading  = new TH2F ("hEtaLeading","#eta distribution of leading",nptbins,ptmin,ptmax, netabins,etamin,etamax); 
  fhEtaLeading->SetYTitle("#eta ");  
 // outputContainer->Add(fhNclustersNtracks);
  outputContainer->Add(fhPtLeading);
  outputContainer->Add(fhPhiLeading);
  outputContainer->Add(fhEtaLeading);
  
  //Correlation with charged hadrons
  if(GetReader()->IsCTSSwitchedOn()) {
    fhDeltaPhiDeltaEtaCharged  = new TH2F
    ("DeltaPhiDeltaEtaCharged","#phi_{trigger} - #phi_{h^{#pm}} vs #eta_{trigger} - #eta_{h^{#pm}}",
     140,-2.,5.,200,-2,2); 
    fhDeltaPhiDeltaEtaCharged->SetXTitle("#Delta #phi");
    fhDeltaPhiDeltaEtaCharged->SetYTitle("#Delta #eta");    
    
    fhPhiCharged  = new TH2F
    ("PhiCharged","#phi_{h^{#pm}}  vs p_{T #pm}",
     nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
    fhPhiCharged->SetYTitle("#phi_{h^{#pm}} (rad)");
    fhPhiCharged->SetXTitle("p_{T #pm} (GeV/c)");
    
    fhEtaCharged  = new TH2F
    ("EtaCharged","#eta_{h^{#pm}}  vs p_{T #pm}",
     nptbins,ptmin,ptmax,netabins,etamin,etamax); 
    fhEtaCharged->SetYTitle("#eta_{h^{#pm}} (rad)");
    fhEtaCharged->SetXTitle("p_{T #pm} (GeV/c)");
    
    fhDeltaPhiCharged  = new TH2F
    ("DeltaPhiCharged","#phi_{trigger} - #phi_{h^{#pm}} vs p_{T trigger}",
     nptbins,ptmin,ptmax,140,-2.,5.); 
    fhDeltaPhiCharged->SetYTitle("#Delta #phi");
    fhDeltaPhiCharged->SetXTitle("p_{T trigger} (GeV/c)");
    
    fhDeltaPhiChargedPt  = new TH2F
    ("DeltaPhiChargedPt","#phi_{trigger} - #phi_{#h^{#pm}} vs p_{T h^{#pm}}",
     nptbins,ptmin,ptmax,140,-2.,5.);
    fhDeltaPhiChargedPt->SetYTitle("#Delta #phi");
    fhDeltaPhiChargedPt->SetXTitle("p_{T h^{#pm}} (GeV/c)");
    
    fhDeltaPhiUeChargedPt  = new TH2F
    ("DeltaPhiUeChargedPt","#phi_{trigger} - #phi_{#Ueh^{#pm}} vs p_{T Ueh^{#pm}}",
     nptbins,ptmin,ptmax,140,-2.,5.);
    fhDeltaPhiUeChargedPt->SetYTitle("#Delta #phi");
    fhDeltaPhiUeChargedPt->SetXTitle("p_{T h^{#pm}} (GeV/c)");
    
    fhDeltaEtaCharged  = new TH2F
    ("DeltaEtaCharged","#eta_{trigger} - #eta_{h^{#pm}} vs p_{T trigger}",
     nptbins,ptmin,ptmax,200,-2,2); 
    fhDeltaEtaCharged->SetYTitle("#Delta #eta");
    fhDeltaEtaCharged->SetXTitle("p_{T trigger} (GeV/c)");
    
    fhPtImbalanceCharged  = 
    new TH2F("CorrelationCharged","z_{trigger h^{#pm}} = p_{T h^{#pm}} / p_{T trigger}",
             nptbins,ptmin,ptmax,200,0.,2.); 
    fhPtImbalanceCharged->SetYTitle("z_{trigger h^{#pm}}");
    fhPtImbalanceCharged->SetXTitle("p_{T trigger}");
    
    fhPtImbalanceUeCharged  = 
    new TH2F("CorrelationUeCharged","z_{trigger h^{#pm}} = p_{T Ueh^{#pm}} / p_{T trigger}",
             nptbins,ptmin,ptmax,200,0.,2.); 
    fhPtImbalanceUeCharged->SetYTitle("z_{trigger Ueh^{#pm}}");
    fhPtImbalanceUeCharged->SetXTitle("p_{T trigger}");
    
    fhPtImbalancePosCharged  = 
    new TH2F("CorrelationPositiveCharged","z_{trigger h^{+}} = p_{T h^{+}} / p_{T trigger}",
             nptbins,ptmin,ptmax,200,0.,2.); 
    fhPtImbalancePosCharged->SetYTitle("z_{trigger h^{+}}");
    fhPtImbalancePosCharged->SetXTitle("p_{T trigger}");
    
    fhPtImbalanceNegCharged  = 
    new TH2F("CorrelationNegativeCharged","z_{trigger h^{-}} = p_{T h^{-}} / p_{T trigger}",
             nptbins,ptmin,ptmax,200,0.,2.); 
    fhPtImbalanceNegCharged->SetYTitle("z_{trigger h^{-}}");
    fhPtImbalanceNegCharged->SetXTitle("p_{T trigger}");
    
    fhPtHbpCharged  = 
    new TH2F("HbpCharged","#xi = ln(1/x_{E}) with charged hadrons",
             nptbins,ptmin,ptmax,200,0.,10.); 
    fhPtHbpCharged->SetYTitle("ln(1/x_{E})");
    fhPtHbpCharged->SetXTitle("p_{T trigger}");
    
    fhPtHbpUeCharged  = 
    new TH2F("HbpUeCharged","#xi = ln(1/x_{E}) with charged hadrons",
             nptbins,ptmin,ptmax,200,0.,10.); 
    fhPtHbpUeCharged->SetYTitle("ln(1/x_{E})");
    fhPtHbpUeCharged->SetXTitle("p_{T trigger}");
    
    fhPtTrigPout  = 
    new TH2F("PtTrigPout","Pout with triggers",
             nptbins,ptmin,ptmax,2*nptbins,-ptmax,ptmax); 
    fhPtTrigPout->SetYTitle("p_{out} (GeV/c)");
    fhPtTrigPout->SetXTitle("p_{T trigger} (GeV/c)"); 
    
    fhPtAssocDeltaPhi  = 
    new TH2F("fhPtAssocDeltaPhi"," charged hadrons vs. delta phi",
             nptbins,ptmin,ptmax,140,-2.,5.); 
    fhPtAssocDeltaPhi->SetXTitle("p_{T h^{#pm}} (GeV/c)");  
    fhPtAssocDeltaPhi->SetYTitle("#Delta #phi (GeV/c)"); 
    
    fhPtTrigCharged  = 
    new TH2F("PtTrigCharged","trgger and charged tracks pt distribution",
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
    outputContainer->Add(fhPtImbalanceCharged) ;
    outputContainer->Add(fhPtImbalancePosCharged) ;
    outputContainer->Add(fhPtImbalanceNegCharged) ;
    outputContainer->Add(fhPtImbalanceUeCharged) ;
    outputContainer->Add(fhPtHbpCharged) ;
    outputContainer->Add(fhPtHbpUeCharged) ;
    outputContainer->Add(fhPtTrigPout) ;
    outputContainer->Add(fhPtAssocDeltaPhi) ;
    outputContainer->Add(fhPtTrigCharged) ;
    
    if(DoEventSelect()){ 
      Int_t nMultiBins = GetMultiBin();
      fhTrigDeltaPhiCharged = new TH2F*[nMultiBins] ;
      fhTrigDeltaEtaCharged = new TH2F*[nMultiBins] ;
      fhTrigCorr = new TH2F*[nMultiBins];
      fhTrigUeCorr = new TH2F*[nMultiBins];
      for(Int_t im=0; im<nMultiBins; im++){
        fhTrigDeltaPhiCharged[im]  = new TH2F 
        (Form("fhTrigDeltaPhiCharged_%d",im),Form("fhTrigDeltaPhiCharged_%d",im), nptbins,ptmin,ptmax, 140,-2.,5.); 
        fhTrigDeltaPhiCharged[im]->SetXTitle("p_{T trigger} (GeV/c)");
        fhTrigDeltaPhiCharged[im]->SetYTitle("#Delta #phi");
        fhTrigDeltaEtaCharged[im]  = new TH2F 
        (Form("fhTrigDeltaEtaCharged_%d",im),Form("fhTrigDeltaEtaCharged_%d",im), nptbins,ptmin,ptmax, 200,-2,2); 
        fhTrigDeltaEtaCharged[im]->SetXTitle("p_{T trigger} (GeV/c)");
        fhTrigDeltaEtaCharged[im]->SetYTitle("#Delta #eta");
        fhTrigCorr[im]  = new TH2F
        (Form("fhTrigPtCorr_%d",im),Form("fhTrigPtCorr_%d",im), nptbins,ptmin,ptmax,200,0.,2.); 
        fhTrigCorr[im]->SetYTitle("z_{trigger h^{#pm}}");
        fhTrigCorr[im]->SetXTitle("p_{T trigger}");
        fhTrigUeCorr[im]  = new TH2F
        (Form("fhTrigPtUeCorr_%d",im),Form("fhTrigPtUeCorr_%d",im), nptbins,ptmin,ptmax,200,0.,2.); 
        fhTrigUeCorr[im]->SetYTitle("z_{trigger h^{#pm}}");
        fhTrigUeCorr[im]->SetXTitle("p_{T trigger}");       

        outputContainer->Add(fhTrigDeltaPhiCharged[im]) ;
        outputContainer->Add(fhTrigDeltaEtaCharged[im]) ;
        outputContainer->Add(fhTrigCorr[im]);
        outputContainer->Add(fhTrigUeCorr[im]);
        
        }
    }
   
    if(fPi0Trigger){
      fhPtPi0DecayRatio  = new TH2F
      ("hPtPi0DecayRatio","p_T of #pi^{0} and the ratio of pt for two decay", 
       nptbins,ptmin,ptmax, 100,0.,2.); 
      fhPtPi0DecayRatio->SetXTitle("p_{T}^{#pi^{0}} (GeV/c)");
      fhPtPi0DecayRatio->SetYTitle("p_{T}^{Decay}/p_{T}^{#pi^{0}}");
      
      fhDeltaPhiDecayCharged  = new TH2F
      ("DeltaPhiDecayCharged","#phi_{Decay} - #phi_{h^{#pm}} vs p_{T Decay}",
       nptbins,ptmin,ptmax,140,-2.,5.); 
      fhDeltaPhiDecayCharged->SetYTitle("#Delta #phi");
      fhDeltaPhiDecayCharged->SetXTitle("p_{T Decay} (GeV/c)");
      
      fhPtImbalanceDecayCharged  = 
      new TH2F("CorrelationDecayCharged","z_{trigger h^{#pm}} = p_{T h^{#pm}} / p_{T Decay}",
               nptbins,ptmin,ptmax,200,0.,2.); 
      fhPtImbalanceDecayCharged->SetYTitle("z_{decay h^{#pm}}");
      fhPtImbalanceDecayCharged->SetXTitle("p_{T decay}");
      
      outputContainer->Add(fhPtPi0DecayRatio) ; 
      outputContainer->Add(fhDeltaPhiDecayCharged) ; 
      outputContainer->Add(fhPtImbalanceDecayCharged) ;
    }    
    
    
    if(fMakeSeveralUE){ 
      fhDeltaPhiUeLeftCharged  = new TH2F
      ("DeltaPhiUeLeftChargedPt","#phi_{trigger} - #phi_{#Ueh^{#pm}} vs p_{T Ueh^{#pm}} with UE left side range of trigger particles",
       nptbins,ptmin,ptmax,140,-2.,5.);
      fhDeltaPhiUeLeftCharged->SetYTitle("#Delta #phi");
      fhDeltaPhiUeLeftCharged->SetXTitle("p_{T h^{#pm}} (GeV/c)");
      outputContainer->Add(fhDeltaPhiUeLeftCharged) ;
      
      fhDeltaPhiUeRightCharged  = new TH2F
      ("DeltaPhiUeRightChargedPt","#phi_{trigger} - #phi_{#Ueh^{#pm}} vs p_{T Ueh^{#pm}} with UE right side range of trigger particles",
       nptbins,ptmin,ptmax,140,-2.,5.);
      fhDeltaPhiUeRightCharged->SetYTitle("#Delta #phi");
      fhDeltaPhiUeRightCharged->SetXTitle("p_{T h^{#pm}} (GeV/c)");
      outputContainer->Add(fhDeltaPhiUeRightCharged) ;
      
      fhPtImbalanceUeLeftCharged  = 
      new TH2F("CorrelationUeChargedLeft","z_{trigger h^{#pm}} = p_{T Ueh^{#pm}} / p_{T trigger} with UE left side of trigger",
               nptbins,ptmin,ptmax,200,0.,2.); 
      fhPtImbalanceUeLeftCharged->SetYTitle("z_{trigger Ueh^{#pm}}");
      fhPtImbalanceUeLeftCharged->SetXTitle("p_{T trigger}");
      outputContainer->Add(fhPtImbalanceUeLeftCharged) ;
      
      fhPtImbalanceUeRightCharged  = 
      new TH2F("CorrelationUeChargedRight","z_{trigger h^{#pm}} = p_{T Ueh^{#pm}} / p_{T trigger} with UE right side of trigger",
               nptbins,ptmin,ptmax,200,0.,2.); 
      fhPtImbalanceUeRightCharged->SetYTitle("z_{trigger Ueh^{#pm}}");
      fhPtImbalanceUeRightCharged->SetXTitle("p_{T trigger}");
      outputContainer->Add(fhPtImbalanceUeRightCharged) ;
      
      fhPtHbpUeLeftCharged  = 
      new TH2F("HbpUeChargedLeft","#xi = ln(1/x_{E}) with charged UE left side of trigger",
               nptbins,ptmin,ptmax,200,0.,10.); 
      fhPtHbpUeLeftCharged->SetYTitle("ln(1/x_{E})");
      fhPtHbpUeLeftCharged->SetXTitle("p_{T trigger}");
      outputContainer->Add(fhPtHbpUeLeftCharged) ;
      
      fhPtHbpUeRightCharged  = 
      new TH2F("HbpUeChargedRight","#xi = ln(1/x_{E}) with charged UE right side of trigger",
               nptbins,ptmin,ptmax,200,0.,10.); 
      fhPtHbpUeRightCharged->SetYTitle("ln(1/x_{E})");
      fhPtHbpUeRightCharged->SetXTitle("p_{T trigger}");
      outputContainer->Add(fhPtHbpUeRightCharged) ;
      
    }  
  }  //Correlation with charged hadrons
  
  //Correlation with neutral hadrons
  if(fNeutralCorr){
    
    fhDeltaPhiDeltaEtaNeutral  = new TH2F
    ("DeltaPhiDeltaEtaNeutral","#phi_{trigger} - #phi_{h^{0}} vs #eta_{trigger} - #eta_{h^{0}}",
     140,-2.,5.,200,-2,2); 
    fhDeltaPhiDeltaEtaNeutral->SetXTitle("#Delta #phi");
    fhDeltaPhiDeltaEtaNeutral->SetYTitle("#Delta #eta");   
	  
    fhPhiNeutral  = new TH2F
    ("PhiNeutral","#phi_{#pi^{0}}  vs p_{T #pi^{0}}",
     nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
    fhPhiNeutral->SetYTitle("#phi_{#pi^{0}} (rad)");
    fhPhiNeutral->SetXTitle("p_{T #pi^{0}} (GeV/c)");
    
    fhEtaNeutral  = new TH2F
    ("EtaNeutral","#eta_{#pi^{0}}  vs p_{T #pi^{0}}",
     nptbins,ptmin,ptmax,netabins,etamin,etamax); 
    fhEtaNeutral->SetYTitle("#eta_{#pi^{0}} (rad)");
    fhEtaNeutral->SetXTitle("p_{T #pi^{0}} (GeV/c)");
    
    fhDeltaPhiNeutral  = new TH2F
    ("DeltaPhiNeutral","#phi_{trigger} - #phi_{#pi^{0}} vs p_{T trigger}",
     nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
    fhDeltaPhiNeutral->SetYTitle("#Delta #phi");
    fhDeltaPhiNeutral->SetXTitle("p_{T trigger} (GeV/c)");
    
    fhDeltaPhiNeutralPt  = new TH2F
    ("DeltaPhiNeutralPt","#phi_{trigger} - #phi_{#pi^{0}} vs p_{T #pi^{0}}}",
     nptbins,ptmin,ptmax,140,-2.,5.); 
    fhDeltaPhiNeutralPt->SetYTitle("#Delta #phi");
    fhDeltaPhiNeutralPt->SetXTitle("p_{T h^{0}} (GeV/c)");
    
    fhDeltaPhiUeNeutralPt  = new TH2F
    ("DeltaPhiUeNeutralPt","#phi_{trigger} - #phi_{#pi^{0}} vs p_{T #pi^{0}}}",
     nptbins,ptmin,ptmax,140,-2.,5.); 
    fhDeltaPhiUeNeutralPt->SetYTitle("#Delta #phi");
    fhDeltaPhiUeNeutralPt->SetXTitle("p_{T h^{0}} (GeV/c)");
    
    fhDeltaEtaNeutral  = new TH2F
    ("DeltaEtaNeutral","#eta_{trigger} - #eta_{#pi^{0}} vs p_{T trigger}",
     nptbins,ptmin,ptmax,200,-2,2); 
    fhDeltaEtaNeutral->SetYTitle("#Delta #eta");
    fhDeltaEtaNeutral->SetXTitle("p_{T trigger} (GeV/c)");
    
    fhPtImbalanceNeutral  = 
    new TH2F("CorrelationNeutral","z_{trigger #pi} = p_{T #pi^{0}} / p_{T trigger}",
             nptbins,ptmin,ptmax,200,0.,2.); 
    fhPtImbalanceNeutral->SetYTitle("z_{trigger #pi^{0}}");
    fhPtImbalanceNeutral->SetXTitle("p_{T trigger}");
    
    fhPtImbalanceUeNeutral  = 
    new TH2F("CorrelationUeNeutral","z_{trigger #pi} = p_{T #pi^{0}} / p_{T trigger}",
             nptbins,ptmin,ptmax,200,0.,2.); 
    fhPtImbalanceUeNeutral->SetYTitle("z_{trigger #pi^{0}}");
    fhPtImbalanceUeNeutral->SetXTitle("p_{T trigger}");
    
    fhPtHbpNeutral  = 
    new TH2F("HbpNeutral","#xi = ln(1/x_{E}) with neutral particles",
             nptbins,ptmin,ptmax,200,0.,10.); 
    fhPtHbpNeutral->SetYTitle("ln(1/x_{E})");
    fhPtHbpNeutral->SetXTitle("p_{T trigger}");
    
    fhPtHbpUeNeutral  = 
    new TH2F("HbpUeNeutral","#xi = ln(1/x_{E}) with neutral particles",
             nptbins,ptmin,ptmax,200,0.,10.); 
    fhPtHbpUeNeutral->SetYTitle("ln(1/x_{E})");
    fhPtHbpUeNeutral->SetXTitle("p_{T trigger}");
    
    
    outputContainer->Add(fhDeltaPhiDeltaEtaNeutral); 
    outputContainer->Add(fhPhiNeutral) ;  
    outputContainer->Add(fhEtaNeutral) ;   
    outputContainer->Add(fhDeltaPhiNeutral) ; 
    outputContainer->Add(fhDeltaPhiNeutralPt) ; 
    outputContainer->Add(fhDeltaPhiUeNeutralPt) ; 
    outputContainer->Add(fhDeltaEtaNeutral) ; 
    outputContainer->Add(fhPtImbalanceNeutral) ;
    outputContainer->Add(fhPtImbalanceUeNeutral) ;  
    outputContainer->Add(fhPtHbpNeutral) ;
    outputContainer->Add(fhPtHbpUeNeutral) ;    
    
    if(fPi0Trigger){
      fhDeltaPhiDecayNeutral  = new TH2F
      ("DeltaPhiDecayNeutral","#phi_{Decay} - #phi_{h^{0}} vs p_{T Decay}",
       nptbins,ptmin,ptmax,140,-2.,5.); 
      fhDeltaPhiDecayNeutral->SetYTitle("#Delta #phi");
      fhDeltaPhiDecayNeutral->SetXTitle("p_{T Decay} (GeV/c)");
      
      fhPtImbalanceDecayNeutral  = 
      new TH2F("CorrelationDecayNeutral","z_{trigger h^{0}} = p_{T h^{0}} / p_{T Decay}",
               nptbins,ptmin,ptmax,200,0.,2.); 
      fhPtImbalanceDecayNeutral->SetYTitle("z_{decay h^{0}}");
      fhPtImbalanceDecayNeutral->SetXTitle("p_{T decay}");
      
      outputContainer->Add(fhDeltaPhiDecayNeutral) ; 
      outputContainer->Add(fhPtImbalanceDecayNeutral) ;
    }
    
    
    if(fMakeSeveralUE){ 
      fhDeltaPhiUeLeftNeutral  = new TH2F
      ("DeltaPhiUeLeftNeutralPt","#phi_{trigger} - #phi_{#Ueh^{0}} vs p_{T h^{0}} with neutral UE left side range of trigger particles",
       nptbins,ptmin,ptmax,140,-2.,5.);
      fhDeltaPhiUeLeftNeutral->SetYTitle("#Delta #phi");
      fhDeltaPhiUeLeftNeutral->SetXTitle("p_{T h^{0}} (GeV/c)");
      outputContainer->Add(fhDeltaPhiUeLeftNeutral) ;
      
      fhDeltaPhiUeRightNeutral  = new TH2F
      ("DeltaPhiUeRightNeutralPt","#phi_{trigger} - #phi_{#Ueh^{0}} vs p_{T Ueh^{0}} with neutral UE right side range of trigger particles",
       nptbins,ptmin,ptmax,140,-2.,5.);
      fhDeltaPhiUeRightNeutral->SetYTitle("#Delta #phi");
      fhDeltaPhiUeRightNeutral->SetXTitle("p_{T h^{0}} (GeV/c)");
      outputContainer->Add(fhDeltaPhiUeRightNeutral) ;
      
      fhPtImbalanceUeLeftNeutral  = 
      new TH2F("CorrelationUeNeutralLeft","z_{trigger h^{0}} = p_{T Ueh^{0}} / p_{T trigger} with neutral UE left side of trigger",
               nptbins,ptmin,ptmax,140,0.,2.); 
      fhPtImbalanceUeLeftNeutral->SetYTitle("z_{trigger Ueh^{0}}");
      fhPtImbalanceUeLeftNeutral->SetXTitle("p_{T trigger}");
      outputContainer->Add(fhPtImbalanceUeLeftNeutral) ;
      
      fhPtImbalanceUeRightNeutral  = 
      new TH2F("CorrelationUeNeutralRight","z_{trigger h^{0}} = p_{T Ueh^{0}} / p_{T trigger} with neutral UE right side of trigger",
               nptbins,ptmin,ptmax,200,0.,2.); 
      fhPtImbalanceUeRightNeutral->SetYTitle("z_{trigger Ueh^{0}}");
      fhPtImbalanceUeRightNeutral->SetXTitle("p_{T trigger}");
      outputContainer->Add(fhPtImbalanceUeRightNeutral) ;
      
      fhPtHbpUeLeftNeutral  = 
      new TH2F("HbpUeNeutralLeft","#xi = ln(1/x_{E}) with neutral UE left side of trigger",
               nptbins,ptmin,ptmax,200,0.,10.); 
      fhPtHbpUeLeftNeutral->SetYTitle("ln(1/x_{E})");
      fhPtHbpUeLeftNeutral->SetXTitle("p_{T trigger}");
      outputContainer->Add(fhPtHbpUeLeftNeutral) ;
      
      fhPtHbpUeRightNeutral  = 
      new TH2F("HbpUeNeutralRight","#xi = ln(1/x_{E}) with neutral UE right side of trigger",
               nptbins,ptmin,ptmax,200,0.,10.); 
      fhPtHbpUeRightNeutral->SetYTitle("ln(1/x_{E})");
      fhPtHbpUeRightNeutral->SetXTitle("p_{T trigger}");
      outputContainer->Add(fhPtHbpUeRightNeutral) ;
      
    }  
    
    //Keep neutral meson selection histograms if requiered
    //Setting done in AliNeutralMesonSelection
    if(GetNeutralMesonSelection()){
      TList * nmsHistos = GetNeutralMesonSelection()->GetCreateOutputObjects() ;
      if(GetNeutralMesonSelection()->AreNeutralMesonSelectionHistosKept())
        for(Int_t i = 0; i < nmsHistos->GetEntries(); i++) outputContainer->Add(nmsHistos->At(i)) ;
      delete nmsHistos;
    }
    
  }//Correlation with neutral hadrons
  
  return outputContainer;
  
}

//____________________________________________________________________________
void AliAnaParticleHadronCorrelation::InitParameters()
{
  
  //Initialize the parameters of the analysis.
  SetInputAODName("PWG4Particle");
  SetAODObjArrayName("Hadrons");  
  AddToHistogramsName("AnaHadronCorr_");

  SetPtCutRange(0.,300);
  fDeltaPhiMinCut = 1.5 ;
  fDeltaPhiMaxCut = 4.5 ;
  fSelectIsolated = kFALSE;
  fMakeSeveralUE = kFALSE;
  fUeDeltaPhiMinCut = 1. ;
  fUeDeltaPhiMaxCut = 1.5 ;
  fNeutralCorr = kFALSE ;
  fPi0Trigger = kFALSE ;

}

//__________________________________________________________________
void AliAnaParticleHadronCorrelation::Print(const Option_t * opt) const
{

  //Print some relevant parameters set for the analysis
  if(! opt)
    return;
  
  printf("**** Print %s %s ****\n", GetName(), GetTitle() ) ;
  AliAnaPartCorrBaseClass::Print(" ");

  printf("Phi trigger particle-Hadron      <     %3.2f\n", fDeltaPhiMaxCut) ; 
  printf("Phi trigger particle-Hadron      >     %3.2f\n", fDeltaPhiMinCut) ;
  printf("Isolated Trigger?  %d\n", fSelectIsolated) ;
  printf("Phi trigger particle-UeHadron    <    %3.2f\n", fUeDeltaPhiMaxCut) ; 
  printf("Phi trigger particle-UeHadron    >    %3.2f\n", fUeDeltaPhiMinCut) ;
  printf("Several UE?  %d\n", fMakeSeveralUE) ;
  printf("Name of AOD Pi0 Branch %s \n",fPi0AODBranchName.Data());
  printf("Do Decay-hadron correlation ?  %d\n", fPi0Trigger) ;

  
} 

//__________________________________________________________________
TObjString* AliAnaParticleHadronCorrelation::GetAnalysisCuts()
{
  //Save parameters used for analysis
  TString parList ; //this will be list of parameters used for this analysis.
  const Int_t buffersize = 255;
  char onePar[buffersize] ;
 
  snprintf(onePar,buffersize,"--- AliAnaPaticleHadronCorrelation ---\n") ;
  parList+=onePar ;	
  snprintf(onePar,buffersize,"Phi trigger particle-Hadron      <     %3.2f ", fDeltaPhiMaxCut) ; 
  parList+=onePar ;
  snprintf(onePar,buffersize,"Phi trigger particle-Hadron      >     %3.2f ", fDeltaPhiMinCut) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"Isolated Trigger?  %d\n", fSelectIsolated) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"Phi trigger particle-UeHadron    <    %3.2f ", fUeDeltaPhiMaxCut) ; 
  parList+=onePar ;
  snprintf(onePar,buffersize,"Phi trigger particle-UeHadron    >    %3.2f ", fUeDeltaPhiMinCut) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"Several UE?  %d\n", fMakeSeveralUE) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"Name of AOD Pi0 Branch %s ",fPi0AODBranchName.Data());
  parList+=onePar ;
  snprintf(onePar,buffersize,"Do Decay-hadron correlation ?  %d", fPi0Trigger) ;
  parList+=onePar ;

  //Get parameters set in base class.
  parList += GetBaseParametersList() ;
  
  //Get parameters set in PID class.
  //parList += GetCaloPID()->GetPIDParametersList() ;
  
  //Get parameters set in FiducialCut class (not available yet)
  //parlist += GetFidCut()->GetFidCutParametersList() 
  
  return new TObjString(parList) ;  

} 


//____________________________________________________________________________
void  AliAnaParticleHadronCorrelation::MakeAnalysisFillAOD()  
{  
  //Particle-Hadron Correlation Analysis, fill AODs
  
  if(!GetInputAODBranch()){
    printf("AliAnaParticleHadronCorrelation::MakeAnalysisFillAOD() - No input particles in AOD with name branch < %s >, STOP \n",GetInputAODName().Data());
    abort();
  }
	
  if(strcmp(GetInputAODBranch()->GetClass()->GetName(), "AliAODPWG4ParticleCorrelation")){
    printf("AliAnaParticleHadronCorrelation::MakeAnalysisFillAOD() - Wrong type of AOD object, change AOD class name in input AOD: It should be <AliAODPWG4ParticleCorrelation> and not <%s> \n",GetInputAODBranch()->GetClass()->GetName());
    abort();
  }
	
  if(GetDebug() > 1){
    printf("AliAnaParticleHadronCorrelation::MakeAnalysisFillAOD() - Begin hadron correlation analysis, fill AODs \n");
    printf("AliAnaParticleHadronCorrelation::MakeAnalysisFillAOD() - In particle branch aod entries %d\n", GetInputAODBranch()->GetEntriesFast());
    printf("AliAnaParticleHadronCorrelation::MakeAnalysisFillAOD() - In CTS aod entries %d\n", GetCTSTracks()->GetEntriesFast());
    printf("AliAnaParticleHadronCorrelation::MakeAnalysisFillAOD() - In EMCAL aod entries %d\n", GetEMCALClusters()->GetEntriesFast());
    printf("AliAnaParticleHadronCorrelation::MakeAnalysisFillAOD() - In PHOS aod entries %d\n", GetPHOSClusters()->GetEntriesFast());
  }
  
  //Loop on stored AOD particles, trigger
  Double_t ptTrig    = 0.;
  Int_t    trigIndex = -1;
  Int_t naod = GetInputAODBranch()->GetEntriesFast();
  //fhNclustersNtracks->Fill(naod, GetCTSTracks()->GetEntriesFast());
  for(Int_t iaod = 0; iaod < naod ; iaod++){
    AliAODPWG4ParticleCorrelation* particle =  (AliAODPWG4ParticleCorrelation*) (GetInputAODBranch()->At(iaod));
    //find the leading particles with highest momentum
    if (particle->Pt()>ptTrig) {
      ptTrig = particle->Pt() ;
      trigIndex = iaod ;
    }
  }//Aod branch loop
	
  //Do correlation with leading particle
  if(trigIndex!=-1){
	  
    AliAODPWG4ParticleCorrelation* particle =  (AliAODPWG4ParticleCorrelation*) (GetInputAODBranch()->At(trigIndex));
    //Make correlation with charged hadrons
    if(GetReader()->IsCTSSwitchedOn() )
      MakeChargedCorrelation(particle, GetCTSTracks(),kFALSE);
    
    TObjArray * pi0list = (TObjArray*) GetAODBranch(fPi0AODBranchName); //For the future, foresee more possible pi0 lists
    if(fNeutralCorr && pi0list && pi0list->GetEntriesFast() > 0)
      MakeNeutralCorrelation(particle, pi0list,kFALSE);
    
  }//Correlate leading
  
  if(GetDebug() > 1) printf("AliAnaParticleHadronCorrelation::MakeAnalysisFillAOD() - End fill AODs \n");
  
}

//____________________________________________________________________________
void  AliAnaParticleHadronCorrelation::MakeAnalysisFillHistograms()  
{  
  //Particle-Hadron Correlation Analysis, fill histograms
  
  if(!GetInputAODBranch()){
    printf("AliAnaParticleHadronCorrelation::MakeAnalysisFillHistograms() - No input particles in AOD with name branch < %s >, STOP \n",GetInputAODName().Data());
    abort();
  }
  
  if(GetDebug() > 1){
    printf("AliAnaParticleHadronCorrelation::MakeAnalysisFillHistograms() - Begin hadron correlation analysis, fill histograms \n");
    printf("AliAnaParticleHadronCorrelation::MakeAnalysisFillHistograms() - In particle branch aod entries %d\n", GetInputAODBranch()->GetEntriesFast());
  }
  
  //Get the vertex and check it is not too large in z
  Double_t v[3] = {0,0,0}; //vertex ;
  GetReader()->GetVertex(v);
  if(!GetMixedEvent() && TMath::Abs(v[2]) > GetZvertexCut()) return ;  
  
  //Loop on stored AOD particles, find leading
  Int_t naod = GetInputAODBranch()->GetEntriesFast();
 // if(naod!=0)fhVertex->Fill(v[0],v[1],v[2]);
  Double_t ptTrig = 0.;
  Int_t trigIndex = -1 ;
  for(Int_t iaod = 0; iaod < naod ; iaod++){	 //loop on input trigger AOD file 
    AliAODPWG4ParticleCorrelation* particle =  (AliAODPWG4ParticleCorrelation*) (GetInputAODBranch()->At(iaod));
    //vertex cut in case of mixing
    if (GetMixedEvent()) {
      Int_t evt=-1;
      Int_t id =-1;
      if     (particle->GetCaloLabel(0)  >= 0 ){
        id=particle->GetCaloLabel(0); 
        if(id >= 0 )evt= GetMixedEvent()-> EventIndexForCaloCluster(id) ;
      }
      else if(particle->GetTrackLabel(0) >= 0 ){
        id=particle->GetTrackLabel(0);
        if(id >= 0 )evt= GetMixedEvent()->EventIndex(id) ;
      }
      else continue;
      
      if (TMath::Abs(GetVertex(evt)[2]) > GetZvertexCut()) 
        return ;
    }

    //check if the particle is isolated or if we want to take the isolation into account
    if(OnlyIsolated() && !particle->IsIsolated()) continue;
    //check if inside the vertex cut
    //find the leading particles with highest momentum
    if (particle->Pt()>ptTrig) {
      ptTrig = particle->Pt() ;
		  trigIndex = iaod ;
    }
  }//finish searching for leading trigger particle
  if(trigIndex!=-1){ //using trigger partilce to do correlations
    AliAODPWG4ParticleCorrelation* particle =  (AliAODPWG4ParticleCorrelation*) (GetInputAODBranch()->At(trigIndex));
	  
    if (GetMixedEvent()) {
      Int_t evt=-1;
      Int_t id = 0;
      if     (particle->GetCaloLabel(0)  >= 0 ){
        id=particle->GetCaloLabel(0); 
        if(id >= 0 )evt= GetMixedEvent()-> EventIndexForCaloCluster(id) ;
      }
      else if(particle->GetTrackLabel(0) >= 0 ){
        id=particle->GetTrackLabel(0);
        if(id >= 0 )evt= GetMixedEvent()->EventIndex(id) ;
      }
      else return;
      
      if (TMath::Abs(GetVertex(evt)[2]) > GetZvertexCut()) 
        return ;
    }
        
    //Fill leading particle histogram   
    fhPtLeading->Fill(particle->Pt());
    Float_t phi = particle->Phi();
    if(phi<0)phi+=TMath::TwoPi();
    fhPhiLeading->Fill(particle->Pt(), phi);
    fhEtaLeading->Fill(particle->Pt(), particle->Eta());
    
    //Make correlation with charged hadrons
    if(GetReader()->IsCTSSwitchedOn() )
      MakeChargedCorrelation(particle, GetCTSTracks(),kTRUE);
    
    TObjArray * pi0list = (TObjArray*) GetAODBranch(fPi0AODBranchName); //For the future, foresee more possible pi0 lists
    if(fNeutralCorr && pi0list && pi0list->GetEntriesFast() > 0)
      MakeNeutralCorrelation(particle, pi0list,kTRUE);
    
  }//Aod branch loop
  
  if(GetDebug() > 1) printf("AliAnaParticleHadronCorrelation::MakeAnalysisFillHistograms() - End fill histograms \n");
  
}

//____________________________________________________________________________
void  AliAnaParticleHadronCorrelation::MakeChargedCorrelation(AliAODPWG4ParticleCorrelation *aodParticle, TObjArray* const pl, const Bool_t bFillHisto)
{  
  // Charged Hadron Correlation Analysis
  if(GetDebug() > 1)printf("AliAnaParticleHadronCorrelation::MakeChargedCorrelation() - Make trigger particle - charged hadron correlation \n");
  
  Int_t evtIndex11 = -1 ; //cluster trigger or pi0 trigger 
  Int_t evtIndex12 = -1 ; // pi0 trigger
  Int_t evtIndex13 = -1 ; // charged trigger
  Int_t indexPhoton1 = -1 ;
  Int_t indexPhoton2 = -1 ;  
  Double_t v[3] = {0,0,0}; //vertex ;
  GetReader()->GetVertex(v);
  if(!GetMixedEvent() && TMath::Abs(v[2]) > GetZvertexCut()) return ;  

  Int_t nTracks = GetCTSTracks()->GetEntriesFast() ;
  
  if (GetMixedEvent()) {
    evtIndex11 = GetMixedEvent()->EventIndexForCaloCluster(aodParticle->GetCaloLabel(0)) ;
    evtIndex12 = GetMixedEvent()->EventIndexForCaloCluster(aodParticle->GetCaloLabel(1)) ;    
    evtIndex13 = GetMixedEvent()->EventIndex(aodParticle->GetTrackLabel(0)) ;
  }
  
  Double_t ptTrig  = aodParticle->Pt();
  Double_t pxTrig  = aodParticle->Px();
  Double_t pyTrig  = aodParticle->Py();
  
  Double_t phiTrig = aodParticle->Phi();
  //Double_t etaTrig = aodParticle->Eta();
  
  Double_t pt   = -100.;
  Double_t px   = -100.;
  Double_t py   = -100.;
  Double_t rat  = -100.; 
  Double_t xE   = -100.; 
  Double_t cosi = -100.; 
  Double_t phi  = -100. ;
  Double_t eta  = -100. ;
  TVector3 p3;  
	
  TObjArray * reftracks = 0x0;
  Int_t nrefs = 0;
	
  Double_t ptDecay1  = 0. ;
  Double_t pxDecay1  = 0. ;
  Double_t pyDecay1  = 0. ;
  Double_t phiDecay1 = 0. ;
  Double_t ptDecay2  = 0. ;
  Double_t pxDecay2  = 0. ;
  Double_t pyDecay2  = 0. ;
  Double_t phiDecay2 = 0. ;
  
  Double_t ratDecay1  = -100.;  
  Double_t ratDecay2  = -100.; 
  Float_t deltaphi    = -100. ;
  Float_t deltaphiDecay1 = -100. ;
  Float_t deltaphiDecay2 = -100. ;
  TObjArray * clusters   = 0x0 ;  
  TLorentzVector photonMom ;	
  
  if(fPi0Trigger){
    indexPhoton1 = aodParticle->GetCaloLabel (0);
    indexPhoton2 = aodParticle->GetCaloLabel (1);
    if(GetDebug() > 1)printf("indexPhoton1 = %d, indexPhoton2 = %d \n", indexPhoton1, indexPhoton2);
    
    if(indexPhoton1!=-1 && indexPhoton2!=-1){
      if(aodParticle->GetDetector()=="EMCAL") clusters = GetEMCALClusters() ;
      else  clusters = GetPHOSClusters() ;
      for(Int_t iclus = 0; iclus < clusters->GetEntriesFast(); iclus++){
        AliVCluster * photon =  (AliVCluster*) (clusters->At(iclus));	
        photon->GetMomentum(photonMom,GetVertex(0)) ;
        if(photon->GetID()==indexPhoton1) {
          ptDecay1  = photonMom.Pt();
          pxDecay1  = photonMom.Px();
          pyDecay1  = photonMom.Py();
          phiDecay1 = photonMom.Phi();
          if(ptTrig && bFillHisto) fhPtPi0DecayRatio->Fill(ptTrig, ptDecay1/ptTrig);
        }
        if(photon->GetID()==indexPhoton2) {
          ptDecay2  = photonMom.Pt();
          pxDecay2  = photonMom.Px();
          pyDecay2  = photonMom.Py();
          phiDecay2 = photonMom.Phi();
          if(ptTrig && bFillHisto) fhPtPi0DecayRatio->Fill(ptTrig, ptDecay2/ptTrig);
        } 
        if(GetDebug() > 1)printf("Photon1 = %f, Photon2 = %f \n", ptDecay1, ptDecay2);
      } //cluster loop        
    } //index of decay photons found
    
  } //make decay-hadron correlation
  
  //Track loop, select tracks with good pt, phi and fill AODs or histograms
  //Int_t currentIndex = -1 ; 
  for(Int_t ipr = 0;ipr < pl->GetEntriesFast() ; ipr ++ ){
    AliVTrack * track = (AliVTrack *) (pl->At(ipr)) ;

    //check if inside the vertex cut
    //printf("charge = %d\n", track->Charge());
    Int_t evtIndex2 = 0 ; 
    if (GetMixedEvent()) {
      evtIndex2 = GetMixedEvent()->EventIndex(track->GetID()) ;
      if (evtIndex11 == evtIndex2 || evtIndex12 == evtIndex2 || evtIndex13 == evtIndex2 ) // photon and track from different events
        continue ; 
      //vertex cut
      if (TMath::Abs(GetVertex(evtIndex2)[2]) > GetZvertexCut()) 
         return ;
      //      if(currentIndex == evtIndex2) // tracks from different event 
      //        continue ;
      //      currentIndex = evtIndex2 ;
    }
    Double_t mom[3] = {track->Px(),track->Py(),track->Pz()};
    p3.SetXYZ(mom[0],mom[1],mom[2]);
    pt   = p3.Pt();
    px   = p3.Px();
    py   = p3.Py();
    eta  = p3.Eta();
    phi  = p3.Phi() ;
    if(phi < 0) phi+=TMath::TwoPi();

    //Select only hadrons in pt range
    if(pt < GetMinPt() || pt > GetMaxPt()) continue ;
    //remove trigger itself for correlation when use charged triggers    
    if( track->GetID() == aodParticle->GetTrackLabel(0) || track->GetID() == aodParticle->GetTrackLabel(1) ||
        track->GetID() == aodParticle->GetTrackLabel(2) || track->GetID() == aodParticle->GetTrackLabel(3)   ) continue ;
       //&&pt==ptTrig && phi==phiTrig && eta==etaTrig) continue ;

    if(IsFiducialCutOn()){
      Bool_t in = GetFiducialCut()->IsInFiducialCut(mom,"CTS") ;
      if(! in ) continue ;
    }    
    //jumped out this event if near side associated partile pt larger than trigger
    if(pt > ptTrig && TMath::Abs(phi-phiTrig)<TMath::PiOver2())  break ;
    rat   = pt/ptTrig ;
    xE    = -(px*pxTrig+py*pyTrig)/(ptTrig*ptTrig);
    if(xE <0.)xE =-xE;
    cosi = TMath::Log(1/xE);
    // printf("rat = %f, xE = %f, cosi =%f \n", rat, xE, cosi);
    // printf("phi = %f \n", phi);
    
     if(fPi0Trigger){
      if(indexPhoton1!=-1 && indexPhoton2!=-1){
        if(ptDecay1) ratDecay1 = pt/ptDecay1 ;
        if(ptDecay2) ratDecay2 = pt/ptDecay2 ; 
        deltaphiDecay1 = phiDecay1-phi;
        deltaphiDecay2 = phiDecay2-phi;
        if(deltaphiDecay1< -TMath::PiOver2()) deltaphiDecay1+=TMath::TwoPi();
        if(deltaphiDecay1>3*TMath::PiOver2()) deltaphiDecay1-=TMath::TwoPi();
        if(deltaphiDecay2< -TMath::PiOver2()) deltaphiDecay2+=TMath::TwoPi();
        if(deltaphiDecay2>3*TMath::PiOver2()) deltaphiDecay2-=TMath::TwoPi();    
      }
    } //do decay-hadron correlation    
            
    //Selection within angular range
    deltaphi = phiTrig-phi;
    if(deltaphi< -TMath::PiOver2()) deltaphi+=TMath::TwoPi();
    if(deltaphi>3*TMath::PiOver2()) deltaphi-=TMath::TwoPi();

    Double_t pout = pt*TMath::Sin(deltaphi) ;
    
    if(GetDebug() > 2)
      printf("AliAnaParticleHadronCorrelation::MakeChargedCorrelation() - Charged hadron: pt %f, phi %f, phi trigger %f. Cuts:  delta phi  %2.2f < %2.2f < %2.2f, pT min %2.2f \n",
             pt,phi, phiTrig,fDeltaPhiMinCut, deltaphi, fDeltaPhiMaxCut, GetMinPt());
    
    if(bFillHisto){
      // Fill Histograms
      fhEtaCharged->Fill(pt,eta);
      fhPhiCharged->Fill(pt,phi);
      fhDeltaEtaCharged->Fill(ptTrig,aodParticle->Eta()-eta);
      fhDeltaPhiCharged->Fill(ptTrig, deltaphi);
      fhDeltaPhiDeltaEtaCharged->Fill(deltaphi,aodParticle->Eta()-eta);
      fhPtAssocDeltaPhi->Fill(pt, deltaphi);
      
      if(GetDebug() > 2 ) printf("AliAnaParticleHadronCorrelation::MakeChargedCorrelation() - Selected charge for momentum imbalance: pt %2.2f, phi %2.2f, eta %2.2f \n",pt,phi,eta);
      //fill different multiplicity histogram
      if(DoEventSelect()){
        for(Int_t im=0; im<GetMultiBin(); im++){
          if(nTracks < ( GetMaxMulti() - GetMinMulti() )/GetMultiBin()*(im+1)){
            fhTrigDeltaPhiCharged[im]->Fill(ptTrig,deltaphi);
            fhTrigDeltaEtaCharged[im]->Fill(ptTrig,aodParticle->Eta()-eta);
         }
        }
      }
      //delta phi cut for correlation
      if( (deltaphi > fDeltaPhiMinCut) && ( deltaphi < fDeltaPhiMaxCut) ) {
        fhDeltaPhiChargedPt->Fill(pt,deltaphi);
        fhPtImbalanceCharged->Fill(ptTrig,xE); 
        fhPtHbpCharged->Fill(ptTrig,cosi);
        fhPtTrigPout->Fill(ptTrig, pout) ;
        fhPtTrigCharged->Fill(ptTrig, pt) ;
        if(track->Charge()>0) fhPtImbalancePosCharged->Fill(ptTrig,xE) ;
        else fhPtImbalanceNegCharged->Fill(ptTrig,xE) ;
        //fill different multiplicity histogram
        if(DoEventSelect()){
          for(Int_t im=0; im<GetMultiBin(); im++){
            if(nTracks < ( GetMaxMulti() - GetMinMulti() )/GetMultiBin()*(im+1))
              fhTrigCorr[im]->Fill(ptTrig,xE);
          }
        } //multiplicity events selection
      } //delta phi cut for correlation
      else if ((deltaphi > fUeDeltaPhiMinCut) && ( deltaphi < fUeDeltaPhiMaxCut)) { //UE study
        fhDeltaPhiUeChargedPt->Fill(pt,deltaphi);
        //fhUePoutPtTrigPtAssoc->Fill(ptTrig, pt, pout) ;
        fhPtImbalanceUeCharged->Fill(ptTrig,xE);
        fhPtHbpUeCharged->Fill(ptTrig,cosi);
        if(DoEventSelect()){
          for(Int_t im=0; im<GetMultiBin(); im++){
            if(nTracks < ( GetMaxMulti() - GetMinMulti() )/GetMultiBin()*(im+1))
              fhTrigUeCorr[im]->Fill(ptTrig,xE);
          }
        } //multiplicity events selection
        
      } //UE study
      
      if(fPi0Trigger){
        if(indexPhoton1!=-1 && indexPhoton2!=-1){
          fhDeltaPhiDecayCharged->Fill(ptDecay1, deltaphiDecay1);
          fhDeltaPhiDecayCharged->Fill(ptDecay2, deltaphiDecay2);
          if(GetDebug() > 1)printf("deltaPhoton1 = %f, deltaPhoton2 = %f \n", deltaphiDecay1, deltaphiDecay2);
          if( (deltaphiDecay1 > fDeltaPhiMinCut) && ( deltaphiDecay1 < fDeltaPhiMaxCut) )
            fhPtImbalanceDecayCharged->Fill(ptDecay1,ratDecay1); 
          if( (deltaphiDecay2 > fDeltaPhiMinCut) && ( deltaphiDecay2 < fDeltaPhiMaxCut) )
            fhPtImbalanceDecayCharged->Fill(ptDecay2,ratDecay2);
          if(GetDebug() > 1)printf("ratPhoton1 = %f, ratPhoton2 = %f \n", pt/ptDecay1, pt/ptDecay2);
        } //index of decay photons found
      } //make decay-hadron correlation          
      
      //several UE calculation 
      if(fMakeSeveralUE){
        if((deltaphi<-fUeDeltaPhiMinCut) && (deltaphi >-fUeDeltaPhiMaxCut)){  
          fhDeltaPhiUeLeftCharged->Fill(pt,deltaphi);
          fhPtImbalanceUeLeftCharged->Fill(ptTrig,rat);
          fhPtHbpUeLeftCharged->Fill(ptTrig,cosi);
        }
        if((deltaphi>fUeDeltaPhiMinCut) && (deltaphi <fUeDeltaPhiMaxCut)){  
          fhDeltaPhiUeRightCharged->Fill(pt,deltaphi);
          fhPtImbalanceUeRightCharged->Fill(ptTrig,rat);
          fhPtHbpUeRightCharged->Fill(ptTrig,cosi);
          
        }
      } //several UE calculation
      
    } //Fill histogram 
    else{
      nrefs++;
      if(nrefs==1){
        reftracks = new TObjArray(0);
        reftracks->SetName(GetAODObjArrayName()+"Tracks");
        reftracks->SetOwner(kFALSE);
      }
      reftracks->Add(track);
    }//aod particle loop
  }// track loop
  
  //Fill AOD with reference tracks, if not filling histograms
  if(!bFillHisto && reftracks) {
    aodParticle->AddObjArray(reftracks);
  }
  
  //delete reftracks;
  
}  



//____________________________________________________________________________
//void  AliAnaParticleHadronCorrelation::MakeNeutralCorrelationFillAOD(AliAODPWG4ParticleCorrelation* const aodParticle,TObjArray* const pl, TString detector)  
//{  
//    // Neutral Pion Correlation Analysis, find pi0, put them in new output aod, if correlation cuts passed
//  if(GetDebug() > 1) printf("AliAnaParticleHadronCorrelation::MakeNeutralCorrelationFillAOD() - Make trigger particle - neutral hadron correlation \n");
//  
//  if(!NewOutputAOD()){
//    printf("AliAnaParticleHadronCorrelation::MakeNeutralCorrelationFillAOD() - Output aod not created, set AOD class name and branch name in the configuration file, STOP! \n");
//    abort();
//  }
//  
//  Double_t phiTrig = aodParticle->Phi();
//  Int_t	tag = 0;
//  TLorentzVector gammai;
//  TLorentzVector gammaj;
//  
//  //Get vertex for photon momentum calculation
//  
//  if(GetReader()->GetDataType()!= AliCaloTrackReader::kMC) 
//  {
//    for (Int_t iev = 0; iev < GetNMixedEvent(); iev++) {
//      if (!GetMixedEvent()) 
//        GetReader()->GetVertex(GetVertex(iev));
//      else 
//        GetMixedEvent()->GetVertexOfEvent(iev)->GetXYZ(GetVertex(iev)); 
//    }
//  }
//  Double_t vertex2[] = {0.0,0.0,0.0} ; //vertex of second input aod
//  if(GetReader()->GetDataType()!= AliCaloTrackReader::kMC) 
//  {
//     if(GetReader()->GetSecondInputAODTree()) GetReader()->GetSecondInputAODVertex(vertex2);
//  }
//	
//    //Cluster loop, select pairs with good pt, phi and fill AODs or histograms
//    //Int_t iEvent= GetReader()->GetEventNumber() ;
//  Int_t nclus = pl->GetEntriesFast();
//  for(Int_t iclus = 0;iclus < nclus ; iclus ++ ){
//    AliVCluster * calo = (AliVCluster *) (pl->At(iclus)) ;
//    
//    Int_t evtIndex1 = 0 ; 
//    if (GetMixedEvent()) {
//      evtIndex1=GetMixedEvent()->EventIndexForCaloCluster(calo->GetID()) ; 
//    }
//    
//
//      //Input from second AOD?
//    Int_t inputi = 0;
//    if     (aodParticle->GetDetector() == "EMCAL" && GetReader()->GetEMCALClustersNormalInputEntries() <= iclus) 
//      inputi = 1 ;
//    else if(aodParticle->GetDetector() == "PHOS"  && GetReader()->GetPHOSClustersNormalInputEntries()  <= iclus) 
//      inputi = 1;
//    
//      //Cluster selection, not charged, with photon or pi0 id and in fiducial cut
//    //FIXME
//    Int_t pdg=0;
//    //if     (inputi == 0 && !SelectCluster(calo, GetVertex(evtIndex1),  gammai, pdg))  
//      continue ;
//    //MEFIX
//    else if(inputi == 1 && !SelectCluster(calo, vertex2, gammai, pdg))  
//      continue ;
//    
//    if(GetDebug() > 2)
//      printf("AliAnaParticleHadronCorrelation::MakeNeutralCorrelationFillAOD() - Neutral cluster in %s: pt %f, phi %f, phi trigger %f. Cuts:  delta phi min %2.2f,  max %2.2f, pT min %2.2f \n",
//             detector.Data(), gammai.Pt(),gammai.Phi(),phiTrig,fDeltaPhiMinCut, fDeltaPhiMaxCut, GetMinPt());
//    
//      //2 gamma overlapped, found with PID
//    if(pdg == AliCaloPID::kPi0){ 
//      
//        //Select only hadrons in pt range
//      if(gammai.Pt() < GetMinPt() || gammai.Pt() > GetMaxPt()) 
//        continue ;
//      
//        //Selection within angular range
//      Float_t phi = gammai.Phi();
//      if(phi < 0) phi+=TMath::TwoPi();
//        //Float_t deltaphi = TMath::Abs(phiTrig-phi);
//        //if( (deltaphi < fDeltaPhiMinCut) || ( deltaphi > fDeltaPhiMaxCut) ) continue ;
//      
//      AliAODPWG4Particle pi0 = AliAODPWG4Particle(gammai);
//        //pi0.SetLabel(calo->GetLabel());
//      pi0.SetPdg(AliCaloPID::kPi0);
//      pi0.SetDetector(detector);
//      
//      if(IsDataMC()){
//        pi0.SetTag(GetMCAnalysisUtils()->CheckOrigin(calo->GetLabel(),GetReader(),inputi));
//        if(GetDebug() > 0) printf("AliAnaParticleHadronCorrelation::MakeNeutralCorrelationFillAOD() - Origin of candidate %d\n",pi0.GetTag());
//      }//Work with stack also 
//       //Set the indeces of the original caloclusters  
//      pi0.SetCaloLabel(calo->GetID(),-1);
//      AddAODParticle(pi0);
//      
//      if(GetDebug() > 2) 
//        printf("AliAnaParticleHadronCorrelation::MakeNeutralCorrelationFillAOD() - Correlated with selected pi0 (pid): pt %f, phi %f\n",pi0.Pt(),pi0.Phi());
//      
//    }// pdg = 111
//    
//      //Make invariant mass analysis
//    else if(pdg == AliCaloPID::kPhoton){	
//        //Search the photon companion in case it comes from  a Pi0 decay
//        //Apply several cuts to select the good pair;
//      for(Int_t jclus = iclus+1; jclus < pl->GetEntries() ; jclus ++ ){
//        AliVCluster * calo2 = (AliVCluster *) (pl->At(jclus)) ;
//        Int_t evtIndex2 = 0 ; 
//        if (GetMixedEvent()) {
//          evtIndex2=GetMixedEvent()->EventIndexForCaloCluster(calo->GetID()) ; 
//        }
//        if (GetMixedEvent() && (evtIndex1 == evtIndex2))
//          continue ;
//        
//          //Input from second AOD?
//        Int_t inputj = 0;
//        if     (aodParticle->GetDetector() == "EMCAL" && GetReader()->GetEMCALClustersNormalInputEntries() <= jclus) 
//          inputj = 1;
//        else if(aodParticle->GetDetector() == "PHOS"  && GetReader()->GetPHOSClustersNormalInputEntries()  <= jclus) 
//          inputj = 1;
//        
//          //Cluster selection, not charged with photon or pi0 id and in fiducial cut
//        Int_t pdgj=0;
//        //FIXME
//        //if     (inputj == 0 && !SelectCluster(calo2, GetVertex(evtIndex2),  gammaj, pdgj))  
//          continue ;
//        //MEFIX
//
//        else if(inputj == 1 && !SelectCluster(calo2, vertex2, gammaj, pdgj))  
//          continue ;
//        //FIXME
//        //if(!SelectCluster(calo2,GetVertex(evtIndex2), gammaj, pdgj)) 
//        //MEFIX
//          continue ;
//        
//        if(pdgj == AliCaloPID::kPhoton ){
//          
//          if((gammai+gammaj).Pt() < GetMinPt() || (gammai+gammaj).Pt() > GetMaxPt()) 
//            continue ;
//          
//            //Selection within angular range
//          Float_t phi = (gammai+gammaj).Phi();
//          if(phi < 0) phi+=TMath::TwoPi();
//            //Float_t deltaphi = TMath::Abs(phiTrig-phi);
//            //if( (deltaphi < fDeltaPhiMinCut) || ( deltaphi > fDeltaPhiMaxCut) ) continue ;
//          
//            //Select good pair (aperture and invariant mass)
//          if(GetNeutralMesonSelection()->SelectPair(gammai, gammaj)){
//            
//            if(GetDebug() > 2 ) printf("AliAnaParticleHadronCorrelation::MakeNeutralCorrelationFillAOD() - Neutral Hadron Correlation: AOD Selected gamma pair: pt %2.2f, phi %2.2f, eta %2.2f, M %2.3f\n",
//                                       (gammai+gammaj).Pt(),(gammai+gammaj).Phi(),(gammai+gammaj).Eta(), (gammai+gammaj).M());
//            
//            TLorentzVector pi0mom = gammai+gammaj;
//            AliAODPWG4Particle pi0 = AliAODPWG4Particle(pi0mom);
//              //pi0.SetLabel(calo->GetLabel());
//            pi0.SetPdg(AliCaloPID::kPi0);
//            pi0.SetDetector(detector);	
//            if(IsDataMC()){
//                //Check origin of the candidates
//              
//              Int_t label1 = calo->GetLabel();
//              Int_t label2 = calo2->GetLabel();
//              Int_t tag1 = GetMCAnalysisUtils()->CheckOrigin(label1, GetReader(), inputi);
//              Int_t tag2 = GetMCAnalysisUtils()->CheckOrigin(label2, GetReader(), inputj);
//              
//              if(GetDebug() > 0) 
//                printf("AliAnaParticleHadronCorrelation::MakeNeutralCorrelationFillAOD() - Origin of: photon1 %d; photon2 %d \n",tag1, tag2);
//              if(GetMCAnalysisUtils()->CheckTagBit(tag1,AliMCAnalysisUtils::kMCPi0Decay) && GetMCAnalysisUtils()->CheckTagBit(tag2,AliMCAnalysisUtils::kMCPi0Decay)){
//                
//                  //Check if pi0 mother is the same
//                if(GetReader()->ReadStack()){ 
//                  TParticle * mother1 = GetMCStack()->Particle(label1);//photon in kine tree
//                  label1 = mother1->GetFirstMother();
//                    //mother1 = GetMCStack()->Particle(label1);//pi0
//                  
//                  TParticle * mother2 = GetMCStack()->Particle(label2);//photon in kine tree
//                  label2 = mother2->GetFirstMother();
//                    //mother2 = GetMCStack()->Particle(label2);//pi0
//                }
//                else if(GetReader()->ReadAODMCParticles()){
//                  AliAODMCParticle * mother1 = (AliAODMCParticle *) (GetReader()->GetAODMCParticles(inputi))->At(label1);//photon in kine tree
//                  label1 = mother1->GetMother();
//                    //mother1 = GetMCStack()->Particle(label1);//pi0
//                  AliAODMCParticle * mother2 = (AliAODMCParticle *) (GetReader()->GetAODMCParticles(inputj))->At(label2);//photon in kine tree
//                  label2 = mother2->GetMother();
//                    //mother2 = GetMCStack()->Particle(label2);//pi0
//                }
//                
//                  //printf("mother1 %d, mother2 %d\n",label1,label2);
//                if(label1 == label2)
//                  GetMCAnalysisUtils()->SetTagBit(tag,AliMCAnalysisUtils::kMCPi0);
//              }
//            }//Work with mc information also   
//            pi0.SetTag(tag);
//              //Set the indeces of the original caloclusters  
//            pi0.SetCaloLabel(calo->GetID(), calo2->GetID());
//            AddAODParticle(pi0);
//            
//            
//          }//Pair selected
//        }//if pair of gammas
//      }//2nd loop
//    }// if pdg = 22
//  }//1st loop
//  
//  if(GetDebug() > 1) 
//    printf("AliAnaParticleHadronCorrelation::MakeNeutralCorrelationFillAOD() - End, %d pi0's found \n",GetOutputAODBranch()->GetEntriesFast());
//}

//____________________________________________________________________________
void  AliAnaParticleHadronCorrelation::MakeNeutralCorrelation(AliAODPWG4ParticleCorrelation * const aodParticle, TObjArray* pi0list, const Bool_t bFillHisto)  
{  
  // Neutral Pion Correlation Analysis
  if(GetDebug() > 1) printf("AliAnaParticleHadronCorrelation::MakeNeutralCorrelation() - Make trigger particle - pi0 correlation, %d pi0's \n",pi0list->GetEntriesFast());
  
  Int_t evtIndex11 = 0 ; 
  Int_t evtIndex12 = 0 ; 
  if (GetMixedEvent()) {
    evtIndex11 = GetMixedEvent()->EventIndexForCaloCluster(aodParticle->GetCaloLabel(0)) ;
    evtIndex12 = GetMixedEvent()->EventIndexForCaloCluster(aodParticle->GetCaloLabel(1)) ;    
  }  
  
  Double_t pt   = -100.;
  Double_t px   = -100.;
  Double_t py   = -100.;
  Double_t rat = -100.; 
  Double_t phi = -100.;
  Double_t eta = -100.;
  Double_t xE  = -100.; 
  Double_t cosi  = -100.; 
  
  Double_t ptTrig  = aodParticle->Pt();
  Double_t phiTrig = aodParticle->Phi();
  Double_t etaTrig = aodParticle->Eta();
  Double_t pxTrig  = aodParticle->Px();
  Double_t pyTrig  = aodParticle->Py();
  
  
  Int_t indexPhoton1 = -1 ;
  Int_t indexPhoton2 = -1 ;    
  Double_t ptDecay1 = 0. ;
  Double_t pxDecay1  = 0. ;
  Double_t pyDecay1  = 0. ;
  Double_t phiDecay1  = 0. ;
  Double_t ptDecay2 = 0. ;
  Double_t pxDecay2  = 0. ;
  Double_t pyDecay2  = 0. ;
  Double_t phiDecay2  = 0. ;
  
  Double_t ratDecay1  = -100.;  
  Double_t ratDecay2  = -100.; 
  Float_t deltaphi = -100. ;
  Float_t deltaphiDecay1 = -100. ;
  Float_t deltaphiDecay2 = -100. ;
  TObjArray * clusters = 0x0 ;  
  TLorentzVector photonMom ;	
  if(fPi0Trigger){
    indexPhoton1 = aodParticle->GetCaloLabel (0);
    indexPhoton2 = aodParticle->GetCaloLabel (1);
    if(GetDebug() > 1)printf("indexPhoton1 = %d, indexPhoton2 = %d \n", indexPhoton1, indexPhoton2);
    
    if(indexPhoton1!=-1 && indexPhoton2!=-1){
      if(aodParticle->GetDetector()=="EMCAL") clusters = GetEMCALClusters() ;
      else                                    clusters = GetPHOSClusters() ;
      for(Int_t iclus = 0; iclus < clusters->GetEntriesFast(); iclus++){
        AliVCluster * photon =  (AliVCluster*) (clusters->At(iclus));	
        photon->GetMomentum(photonMom,GetVertex(0)) ;
        if(photon->GetID()==indexPhoton1) {
          ptDecay1  = photonMom.Pt();
          pxDecay1  = photonMom.Px();
          pyDecay1  = photonMom.Py();
          phiDecay1 = photonMom.Phi();
        }
        if(photon->GetID()==indexPhoton2) {
          ptDecay2  = photonMom.Pt();
          pxDecay2  = photonMom.Px();
          pyDecay2  = photonMom.Py();
          phiDecay2 = photonMom.Phi();
        } 
        if(GetDebug() > 1)printf("Photon1 = %f, Photon2 = %f \n", ptDecay1, ptDecay2);
      } //photonAOD loop        
    } //index of decay photons found
    if(ptTrig && bFillHisto) fhPtPi0DecayRatio->Fill(ptTrig, ptDecay1/ptTrig, ptDecay2/ptTrig);
  } //make decay-hadron correlation
  
  TObjArray * refpi0    =0x0;
  Int_t nrefs = 0;
  
  //Loop on stored AOD pi0
  Int_t naod = pi0list->GetEntriesFast();
  if(GetDebug() > 0) printf("AliAnaParticleHadronCorrelation::MakeNeutralCorrelationFillHistograms() -  aod branch entries %d\n", naod);
  for(Int_t iaod = 0; iaod < naod ; iaod++){
    AliAODPWG4Particle* pi0 =  (AliAODPWG4Particle*) (pi0list->At(iaod));
    
    Int_t evtIndex2 = 0 ; 
    Int_t evtIndex3 = 0 ; 
    if (GetMixedEvent()) {
      evtIndex2 = GetMixedEvent()->EventIndexForCaloCluster(pi0->GetCaloLabel(0)) ;
      evtIndex3 = GetMixedEvent()->EventIndexForCaloCluster(pi0->GetCaloLabel(1)) ;
      
      if (evtIndex11 == evtIndex2 || evtIndex12 == evtIndex2 || evtIndex11 == evtIndex3 || evtIndex12 == evtIndex3) // trigger and pi0 are not from different events
        continue ; 
    }      
    
    //Int_t pdg = pi0->GetPdg();
    //if(pdg != AliCaloPID::kPi0) continue;  
    
    pt  = pi0->Pt();
    px  = pi0->Px();
    py  = pi0->Py();    
    if(pt < GetMinPt() || pt > GetMaxPt()) continue ;
    //jumped out this event if near side associated partile pt larger than trigger
    if(pt > ptTrig && TMath::Abs(phi-phiTrig)<TMath::PiOver2())  break ;

    //Selection within angular range
    phi = pi0->Phi();
    //Float_t deltaphi = TMath::Abs(phiTrig-phi);
    //if( (deltaphi < fDeltaPhiMinCut) || ( deltaphi > fDeltaPhiMaxCut) ) continue ;
    
    if(bFillHisto){
      
      deltaphi = phiTrig-phi;
      if(deltaphi<-TMath::PiOver2()) deltaphi+=TMath::TwoPi();
      if(deltaphi>3*TMath::PiOver2()) deltaphi-=TMath::TwoPi();
      
      rat = pt/ptTrig ;
      phi = pi0->Phi() ;
      eta = pi0->Eta() ;
      xE   = -(px*pxTrig+py*pyTrig)/(ptTrig*ptTrig);
      if(xE <0.)xE =-xE;
      cosi = TMath::Log(1/xE);
      
      if(fPi0Trigger){
        if(indexPhoton1!=-1 && indexPhoton2!=-1){
          if(ptDecay1) ratDecay1 = pt/ptDecay1 ;
          if(ptDecay2) ratDecay2 = pt/ptDecay2 ; 
          deltaphiDecay1 = phiDecay1-phi;
          deltaphiDecay2 = phiDecay2-phi;
          if(deltaphiDecay1< -TMath::PiOver2()) deltaphiDecay1+=TMath::TwoPi();
          if(deltaphiDecay1>3*TMath::PiOver2()) deltaphiDecay1-=TMath::TwoPi();
          if(deltaphiDecay2< -TMath::PiOver2()) deltaphiDecay2+=TMath::TwoPi();
          if(deltaphiDecay2>3*TMath::PiOver2()) deltaphiDecay2-=TMath::TwoPi();   
          fhDeltaPhiDecayNeutral->Fill(ptDecay1, deltaphiDecay1);
          fhDeltaPhiDecayNeutral->Fill(ptDecay2, deltaphiDecay2);
          if(GetDebug() > 1)printf("deltaPhoton1 = %f, deltaPhoton2 = %f \n", deltaphiDecay1, deltaphiDecay2);
          if( (deltaphiDecay1 > fDeltaPhiMinCut) && ( deltaphiDecay1 < fDeltaPhiMaxCut) )
            fhPtImbalanceDecayNeutral->Fill(ptDecay1,ratDecay1); 
          if( (deltaphiDecay2 > fDeltaPhiMinCut) && ( deltaphiDecay2 < fDeltaPhiMaxCut) )
            fhPtImbalanceDecayNeutral->Fill(ptDecay2,ratDecay2);
          if(GetDebug() > 1)printf("ratPhoton1 = %f, ratPhoton2 = %f \n", pt/ptDecay1, pt/ptDecay2);
        }
      } //do decay-hadron correlation
      
      fhEtaNeutral->Fill(pt,eta);
      fhPhiNeutral->Fill(pt,phi);
      fhDeltaEtaNeutral->Fill(ptTrig,etaTrig-eta);
      fhDeltaPhiNeutral->Fill(ptTrig,deltaphi);
      fhDeltaPhiDeltaEtaNeutral->Fill(deltaphi,etaTrig-eta);
      
      //delta phi cut for correlation
      if( (deltaphi > fDeltaPhiMinCut) && ( deltaphi < fDeltaPhiMaxCut) ) {
        fhDeltaPhiNeutralPt->Fill(pt,deltaphi);
        fhPtImbalanceNeutral->Fill(ptTrig,rat); 
        fhPtHbpNeutral->Fill(ptTrig,cosi); 
      }
      else {
        fhDeltaPhiUeNeutralPt->Fill(pt,deltaphi);
        fhPtImbalanceUeNeutral->Fill(ptTrig,rat);
        fhPtHbpUeNeutral->Fill(ptTrig,cosi); 
      }
      //several UE calculation 
      if(fMakeSeveralUE){
        if((deltaphi<-fUeDeltaPhiMinCut) && (deltaphi >-fUeDeltaPhiMaxCut)){  
          fhDeltaPhiUeLeftNeutral->Fill(pt,deltaphi);
          fhPtImbalanceUeLeftNeutral->Fill(ptTrig,rat);
          fhPtHbpUeLeftNeutral->Fill(ptTrig,cosi);
        }
        if((deltaphi>fUeDeltaPhiMinCut) && (deltaphi <fUeDeltaPhiMaxCut)){  
          fhDeltaPhiUeRightNeutral->Fill(pt,deltaphi);
          fhPtImbalanceUeRightNeutral->Fill(ptTrig,rat);
          fhPtHbpUeRightNeutral->Fill(ptTrig,cosi);
        }
      } //several UE calculation
	  }
    else{
      nrefs++;
      if(nrefs==1){
        refpi0 = new TObjArray(0);
        refpi0->SetName(GetAODObjArrayName()+"Pi0s");
        refpi0->SetOwner(kFALSE);
      }
      refpi0->Add(pi0);
    }//put references in trigger AOD 
      
     //if(GetDebug() > 2 ) printf("AliAnaParticleHadronCorrelation::MakeNeutralCorrelation() - Selected neutral for momentum imbalance: pt %2.2f, phi %2.2f, eta %2.2f \n",pt,phi,eta);
      
    }//loop
}
  

//____________________________________________________________________________
//Bool_t  AliAnaParticleHadronCorrelation::SelectCluster(AliVCluster * calo, Double_t *vertex, TLorentzVector & mom, Int_t & pdg) {
//  //Select cluster depending on its pid and acceptance selections
//  
//  //Skip matched clusters with tracks
//  if(IsTrackMatched(calo)) return kFALSE;
//  
//  TString detector = "";
//  if     (calo->IsPHOS())  detector= "PHOS";
//  else if(calo->IsEMCAL()) detector= "EMCAL";
//		
//  //Check PID
//  calo->GetMomentum(mom,vertex);//Assume that come from vertex in straight line
//  pdg = AliCaloPID::kPhoton;   
//  if(IsCaloPIDOn()){
//    //Get most probable PID, 2 options check PID weights (in MC this option is mandatory)
//    //or redo PID, recommended option for EMCal.
//    
//    if(!IsCaloPIDRecalculationOn() || GetReader()->GetDataType() == AliCaloTrackReader::kMC )
//      pdg = GetCaloPID()->GetPdg(detector,calo->GetPID(),mom.E());//PID with weights
//    else
//      pdg = GetCaloPID()->GetPdg(detector,mom,calo);//PID recalculated
//    
//    if(GetDebug() > 5) printf("AliAnaParticleHadronCorrelation::SelectCluster() - PDG of identified particle %d\n",pdg);
//    
//    //If it does not pass pid, skip
//    if(pdg != AliCaloPID::kPhoton && pdg != AliCaloPID::kPi0) {
//      return kFALSE ;
//    }
//  }//PID on
//  
//  //Check acceptance selection
//  if(IsFiducialCutOn()){
//    Bool_t in = GetFiducialCut()->IsInFiducialCut(mom,detector) ;
//    if(! in ) return kFALSE ;
//  }
//  
//  if(GetDebug() > 5) printf("AliAnaParticleHadronCorrelation::SelectCluster() - Correlation photon selection cuts passed: pT %3.2f, pdg %d\n",mom.Pt(), pdg);
//  
//  return kTRUE;
//  
//}
