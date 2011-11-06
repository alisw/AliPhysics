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
// 6. Add the possibility for event selection analysis based on vertex and multiplicity bins (10/10/2010)
// 7. change the way of delta phi cut for UE study due to memory issue (reduce histograms)
// 8. Add the possibility to request the absolute leading particle at the near side or not, set trigger bins, general clean-up (08/2011)
//////////////////////////////////////////////////////////////////////////////


// --- ROOT system ---
//#include "TClonesArray.h"
#include "TClass.h"
#include "TMath.h"
#include "TH3D.h"
#include "TDatabasePDG.h"

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


//___________________________________________________________________
  AliAnaParticleHadronCorrelation::AliAnaParticleHadronCorrelation(): 
    AliAnaPartCorrBaseClass(),
    fDeltaPhiMaxCut(0.),            fDeltaPhiMinCut(0.),   
    fSelectIsolated(0),             fMakeSeveralUE(0),              
    fUeDeltaPhiMaxCut(0.),          fUeDeltaPhiMinCut(0.), 
    fPi0AODBranchName(""),          fNeutralCorr(0),       
    fPi0Trigger(0),                 fMakeAbsoluteLeading(0),        
    fLeadingTriggerIndex(-1),       
    fNAssocPtBins(0),               fAssocPtBinLimit(),
    //Histograms
    fhPtLeading(0),                 fhPhiLeading(0),       
    fhEtaLeading(0),                fhDeltaPhiDeltaEtaCharged(0),
    fhPhiCharged(0),                fhEtaCharged(0), 
    fhDeltaPhiCharged(0),           fhDeltaEtaCharged(0), 
    fhDeltaPhiChargedPt(0),         fhDeltaPhiUeChargedPt(0), 
    fhPtImbalanceCharged(0),        fhPtImbalanceUeCharged(0),
    fhPtImbalancePosCharged(0),     fhPtImbalanceNegCharged(0),
    fhPtHbpCharged(0),              fhPtHbpUeCharged(0),
    fhDeltaPhiUeLeftCharged(0),     fhDeltaPhiUeRightCharged(0),
    fhPtImbalanceUeLeftCharged(0),  fhPtImbalanceUeRightCharged(0),
    fhPtHbpUeLeftCharged(0),        fhPtHbpUeRightCharged(0), 
    fhPtTrigPout(0),                fhPtTrigCharged(0),
    fhTrigDeltaPhiCharged(0x0),     fhTrigDeltaEtaCharged(0x0),
    fhTrigCorr(0x0),                fhTrigUeCorr(0x0),
    fhAssocPt(0),                   fhAssocPtBkg(0),  
    fhDeltaPhiAssocPtBin(0),        fhDeltaPhiBradAssocPtBin(0),
    fhXEAssocPtBin(0),
    fhDeltaPhiDeltaEtaNeutral(0), 
    fhPhiNeutral(0),                fhEtaNeutral(0), 
    fhDeltaPhiNeutral(0),           fhDeltaEtaNeutral(0),
    fhDeltaPhiNeutralPt(0),         fhDeltaPhiUeNeutralPt(0), 
    fhPtImbalanceNeutral(0),        fhPtImbalanceUeNeutral(0),
    fhPtHbpNeutral(0),              fhPtHbpUeNeutral(0),
    fhDeltaPhiUeLeftNeutral(0),     fhDeltaPhiUeRightNeutral(0),
    fhPtImbalanceUeLeftNeutral(0),  fhPtImbalanceUeRightNeutral(0),
    fhPtHbpUeLeftNeutral(0),        fhPtHbpUeRightNeutral(0),
    fhPtPi0DecayRatio(0),
    fhDeltaPhiDecayCharged(0),      fhPtImbalanceDecayCharged(0), 
    fhDeltaPhiDecayNeutral(0),      fhPtImbalanceDecayNeutral(0),
    fh2phiLeadingParticle(0x0),
    fhMCLeadingCount(0),
    fhMCEtaCharged(0),              fhMCPhiCharged(0), 
    fhMCDeltaEtaCharged(0),         fhMCDeltaPhiCharged(0x0),
    fhMCDeltaPhiDeltaEtaCharged(0), fhMCDeltaPhiChargedPt(0),
    fhMCPtImbalanceCharged(0),
    fhMCPtHbpCharged(0),
    fhMCPtTrigPout(0),
    fhMCPtAssocDeltaPhi(0)
{
  //Default Ctor
  
  //Initialize parameters
  InitParameters();
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
  snprintf(onePar,buffersize,"Select absolute leading for cluster triggers ?  %d\n", fMakeAbsoluteLeading) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"Associated particle pt bins  %d: ", fNAssocPtBins) ;
  parList+=onePar ;
  for (Int_t ibin = 0; ibin<fNAssocPtBins; ibin++) {
    snprintf(onePar,buffersize,"bin %d = [%2.1f,%2.1f];", ibin, fAssocPtBinLimit[ibin], fAssocPtBinLimit[ibin+1]) ;
  }
  parList+=onePar ;
  
  //Get parameters set in base class.
  parList += GetBaseParametersList() ;
  
  //Get parameters set in PID class.
  //parList += GetCaloPID()->GetPIDParametersList() ;
  
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
  
  Int_t   nptbins = GetHistoPtBins(); Int_t  nphibins = GetHistoPhiBins(); Int_t   netabins = GetHistoEtaBins();
  Float_t ptmax   = GetHistoPtMax();  Float_t phimax  = GetHistoPhiMax();  Float_t etamax   = GetHistoEtaMax();
  Float_t ptmin   = GetHistoPtMin();  Float_t phimin  = GetHistoPhiMin();  Float_t etamin   = GetHistoEtaMin();	
  
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
    outputContainer->Add(fhPtTrigCharged) ;
    
    if(DoEventSelect()){ 
      Int_t nMultiBins = GetMultiBin();
      fhTrigDeltaPhiCharged = new TH2F*[nMultiBins] ;
      fhTrigDeltaEtaCharged = new TH2F*[nMultiBins] ;
      fhTrigCorr            = new TH2F*[nMultiBins];
      fhTrigUeCorr          = new TH2F*[nMultiBins];
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
    
    fhAssocPt           = new TH2F("fhAssocPt", " Trigger p_{T} vs associated hadron p_{T}",
                                   nptbins, ptmin, ptmax,nptbins,ptmin,ptmax);
    fhAssocPt->SetXTitle("p_{T trigger}");
    fhAssocPt->SetYTitle("p_{T associated}");
    outputContainer->Add(fhAssocPt) ;
    
    fhAssocPtBkg        = new TH2F("fhAssocPtBkg", " Trigger p_{T} vs associated hadron p_{T} from background",
                                   nptbins, ptmin, ptmax,nptbins,ptmin,ptmax);
    fhAssocPtBkg->SetXTitle("p_{T trigger}");
    fhAssocPtBkg->SetYTitle("p_{T associated}");
    outputContainer->Add(fhAssocPtBkg) ;
    
    fhDeltaPhiAssocPtBin     = new TH2F*[fNAssocPtBins] ;
    fhDeltaPhiBradAssocPtBin = new TH2F*[fNAssocPtBins] ;
    fhXEAssocPtBin           = new TH2F*[fNAssocPtBins] ;
    for(Int_t i = 0 ; i < fNAssocPtBins ; i++){
      fhDeltaPhiAssocPtBin[i] = new TH2F(Form("fhDeltaPhiPtAssocPt%2.1f_%2.1f", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1]), 
                                         Form("#Delta #phi vs p_{T trigger} for associated p_{T} bin [%2.1f,%2.1f]", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1]), 
                                         nptbins, ptmin, ptmax,140,-2.,5.);
      fhDeltaPhiAssocPtBin[i]->SetXTitle("p_{T trigger}");
      fhDeltaPhiAssocPtBin[i]->SetYTitle("#Delta #phi");
      
      fhDeltaPhiBradAssocPtBin[i] = new TH2F(Form("fhDeltaPhiBradPtAssocPt%2.1f_%2.1f", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1]), 
                                             Form("atan2(sin(#Delta #phi), cos(#Delta #phi))/#pi vs p_{T trigger} for associated p_{T} bin [%2.1f,%2.1f]", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1]), 
                                             nptbins, ptmin, ptmax,288, -1.0/3.0, 5.0/3.0);
      fhDeltaPhiBradAssocPtBin[i]->SetXTitle("p_{T trigger}");
      fhDeltaPhiBradAssocPtBin[i]->SetYTitle("atan2(sin(#Delta #phi), cos(#Delta #phi))/#pi");
      
      
      fhXEAssocPtBin[i]       = new TH2F(Form("fhXEAssocPtBin%1.f_%1.f", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1]), 
                                         Form("x_[E] vs p_{T trigger} for associated p_{T} bin [%2.1f,%2.1f]", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1]), 
                                         nptbins, ptmin, ptmax,50, 0.0, 2.0);
      fhXEAssocPtBin[i]->SetXTitle("p_{T trigger}");
      fhXEAssocPtBin[i]->SetYTitle("x_{E}");
      
      outputContainer->Add(fhDeltaPhiAssocPtBin[i]) ;
      outputContainer->Add(fhDeltaPhiBradAssocPtBin[i]) ;
      outputContainer->Add(fhXEAssocPtBin[i]);
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
    
    //if data is MC, fill more histograms
    if(IsDataMC()){
      fh2phiLeadingParticle=new TH2F("fh2phiLeadingParticle","#phi resolustion for trigger particles",nptbins,ptmin,ptmax,100,-1,1);
      fh2phiLeadingParticle->GetXaxis()->SetTitle("p_{T gen Leading} (GeV/c)");
      fh2phiLeadingParticle->GetYaxis()->SetTitle("(#phi_{rec}-#phi_{gen})/#phi_{gen}");
      
      fhMCLeadingCount=new TH1F("MCLeadingTriggerCount","MCLeadingTriggerCount",nptbins,ptmin,ptmax);
      fhMCLeadingCount->SetXTitle("p_{T trig}");
      
      fhMCEtaCharged  = new TH2F
      ("MCEtaCharged","MC #eta_{h^{#pm}}  vs p_{T #pm}",
       nptbins,ptmin,ptmax,netabins,etamin,etamax); 
      fhMCEtaCharged->SetYTitle("#eta_{h^{#pm}} (rad)");
      fhMCEtaCharged->SetXTitle("p_{T #pm} (GeV/c)");
      
      fhMCPhiCharged  = new TH2F
      ("MCPhiCharged","#MC phi_{h^{#pm}}  vs p_{T #pm}",
       200,ptmin,ptmax,nphibins,phimin,phimax); 
      fhMCPhiCharged->SetYTitle("MC #phi_{h^{#pm}} (rad)");
      fhMCPhiCharged->SetXTitle("p_{T #pm} (GeV/c)");
      
      fhMCDeltaPhiDeltaEtaCharged  = new TH2F
      ("MCDeltaPhiDeltaEtaCharged","#MC phi_{trigger} - #phi_{h^{#pm}} vs #eta_{trigger} - #eta_{h^{#pm}}",
       140,-2.,5.,200,-2,2); 
      fhMCDeltaPhiDeltaEtaCharged->SetXTitle("#Delta #phi");
      fhMCDeltaPhiDeltaEtaCharged->SetYTitle("#Delta #eta");    
      
      fhMCDeltaEtaCharged  = new TH2F
      ("MCDeltaEtaCharged","MC #eta_{trigger} - #eta_{h^{#pm}} vs p_{T trigger} and p_{T assoc}",
       nptbins,ptmin,ptmax,200,-2,2); 
      fhMCDeltaEtaCharged->SetYTitle("#Delta #eta");
      fhMCDeltaEtaCharged->SetXTitle("p_{T trigger} (GeV/c)");
      
      fhMCDeltaPhiCharged  = new TH2F
      ("MCDeltaPhiCharged","#phi_{trigger} - #phi_{h^{#pm}} vs p_{T trigger}",
       nptbins,ptmin,ptmax,140,-2.,5.); 
      fhMCDeltaPhiCharged->SetYTitle("#Delta #phi");
      fhMCDeltaPhiCharged->SetXTitle("p_{T trigger} (GeV/c)");
      
      fhMCDeltaPhiChargedPt  = new TH2F
      ("MCDeltaPhiChargedPt","MC #phi_{trigger} - #phi_{#h^{#pm}} vs p_{T h^{#pm}}",
       nptbins,ptmin,ptmax,140,-2.,5.);
      fhMCDeltaPhiChargedPt->SetYTitle("#Delta #phi");
      fhMCDeltaPhiChargedPt->SetXTitle("p_{T h^{#pm}} (GeV/c)");
      
      fhMCPtImbalanceCharged  = 
      new TH2F("MCCorrelationCharged","z_{trigger h^{#pm}} = p_{T h^{#pm}} / p_{T trigger}",
               nptbins,ptmin,ptmax,200,0.,2.); 
      fhMCPtImbalanceCharged->SetYTitle("z_{trigger h^{#pm}}");
      fhMCPtImbalanceCharged->SetXTitle("p_{T trigger}");  
      
      fhMCPtHbpCharged  = 
      new TH2F("MCHbpCharged","MC #xi = ln(1/x_{E}) with charged hadrons",
               nptbins,ptmin,ptmax,200,0.,10.); 
      fhMCPtHbpCharged->SetYTitle("ln(1/x_{E})");
      fhMCPtHbpCharged->SetXTitle("p_{T trigger}");
      
      fhMCPtTrigPout  = 
      new TH2F("MCPtTrigPout","AOD MC Pout with triggers",
               nptbins,ptmin,ptmax,2*nptbins,-ptmax,ptmax); 
      fhMCPtTrigPout->SetYTitle("p_{out} (GeV/c)");
      fhMCPtTrigPout->SetXTitle("p_{T trigger} (GeV/c)"); 
      
      fhMCPtAssocDeltaPhi  = 
      new TH2F("fhMCPtAssocDeltaPhi","AOD MC delta phi with associated charged hadrons",
               nptbins,ptmin,ptmax,140,-2.,5.); 
      fhMCPtAssocDeltaPhi->SetYTitle("#Delta #phi");
      fhMCPtAssocDeltaPhi->SetXTitle("p_{T trigger} (GeV/c)"); 
      
      outputContainer->Add(fh2phiLeadingParticle);
      outputContainer->Add(fhMCLeadingCount);
      outputContainer->Add(fhMCDeltaPhiDeltaEtaCharged);
      outputContainer->Add(fhMCPhiCharged) ;
      outputContainer->Add(fhMCEtaCharged) ;
      outputContainer->Add(fhMCDeltaEtaCharged) ;
      outputContainer->Add(fhMCDeltaPhiCharged) ; 
      
      outputContainer->Add(fhMCDeltaPhiChargedPt) ;
      outputContainer->Add(fhMCPtImbalanceCharged) ;
      outputContainer->Add(fhMCPtHbpCharged) ;
      outputContainer->Add(fhMCPtTrigPout) ;
      outputContainer->Add(fhMCPtAssocDeltaPhi) ;      
    } //for MC histogram
    
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
  fMakeAbsoluteLeading  = kTRUE;
  
  fNAssocPtBins        = 7  ;
  fAssocPtBinLimit[0]  = 2.  ; 
  fAssocPtBinLimit[1]  = 4.  ; 
  fAssocPtBinLimit[2]  = 6.  ; 
  fAssocPtBinLimit[3]  = 8.  ; 
  fAssocPtBinLimit[4]  = 10. ; 
  fAssocPtBinLimit[5]  = 13. ; 
  fAssocPtBinLimit[6]  = 16. ;
  fAssocPtBinLimit[7]  = 20. ;
  fAssocPtBinLimit[8]  = 100.;
  fAssocPtBinLimit[9]  = 200.;
  
}

//__________________________________________________________
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
  
  //Get the vertex and check it is not too large in z
  Double_t v[3] = {0,0,0}; //vertex ;
  GetReader()->GetVertex(v);
  if(!GetMixedEvent() && TMath::Abs(v[2]) > GetZvertexCut()) return ;   
  
  //Loop on stored AOD particles, find leading trigger
  Double_t ptTrig      = GetMinPt() ;
  fLeadingTriggerIndex = -1 ;
  Int_t    naod        = GetInputAODBranch()->GetEntriesFast() ;
  for(Int_t iaod = 0; iaod < naod ; iaod++){
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
    
    // find the leading particles with highest momentum
    if (particle->Pt() > ptTrig) {
      ptTrig               = particle->Pt() ;
      fLeadingTriggerIndex = iaod ;
    }
  }// finish search of leading trigger particle
	
  //Do correlation with leading particle
  if(fLeadingTriggerIndex >= 0){
	  
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
  Double_t ptTrig    = GetMinPt() ;
  if(fLeadingTriggerIndex < 0){//Search leading if not done before
    Int_t    naod      = GetInputAODBranch()->GetEntriesFast() ;
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
      if (particle->Pt() > ptTrig) {
        ptTrig               = particle->Pt() ;
        fLeadingTriggerIndex = iaod ;
      }
    }//finish search of leading trigger particle
  }//Search leading if not done before
  
  if(fLeadingTriggerIndex >= 0 ){ //using trigger particle to do correlations
    AliAODPWG4ParticleCorrelation* particle =  (AliAODPWG4ParticleCorrelation*) (GetInputAODBranch()->At(fLeadingTriggerIndex));
	  
    //check if the particle is isolated or if we want to take the isolation into account
    if(OnlyIsolated() && !particle->IsIsolated()) return;
    
    //Make correlation with charged hadrons
    Bool_t okcharged = kTRUE;
    Bool_t okneutral = kTRUE;
    if(GetReader()->IsCTSSwitchedOn() ){
      okcharged = MakeChargedCorrelation(particle, GetCTSTracks(),kTRUE);
      if(IsDataMC()){      
        MakeMCChargedCorrelation(particle);
      }
    }  
    
    TObjArray * pi0list = (TObjArray*) GetAODBranch(fPi0AODBranchName); //For the future, foresee more possible pi0 lists
    if(fNeutralCorr && pi0list){
      if(pi0list->GetEntriesFast() > 0)
        okneutral = MakeNeutralCorrelation(particle, pi0list,kTRUE);
    }
    
    // Fill leading particle histogram if correlation went well and 
    // no problem was found, like not absolute leading, or bad vertex in mixing.
    if(okcharged && okneutral){
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

//_________________________________________________________________________________________________________
Bool_t  AliAnaParticleHadronCorrelation::MakeChargedCorrelation(AliAODPWG4ParticleCorrelation *aodParticle, 
                                                                TObjArray* pl, const Bool_t bFillHisto)
{  
  // Charged Hadron Correlation Analysis
  if(GetDebug() > 1)printf("AliAnaParticleHadronCorrelation::MakeChargedCorrelation() - Make trigger particle - charged hadron correlation \n");
  
  Int_t evtIndex11   = -1 ; //cluster trigger or pi0 trigger 
  Int_t evtIndex12   = -1 ; // pi0 trigger
  Int_t evtIndex13   = -1 ; // charged trigger
  Int_t indexPhoton1 = -1 ;
  Int_t indexPhoton2 = -1 ;  
  
  Double_t v[3]      = {0,0,0}; //vertex ;
  GetReader()->GetVertex(v);
  
  if (GetMixedEvent()) {
    evtIndex11 = GetMixedEvent()->EventIndexForCaloCluster(aodParticle->GetCaloLabel(0)) ;
    evtIndex12 = GetMixedEvent()->EventIndexForCaloCluster(aodParticle->GetCaloLabel(1)) ;    
    evtIndex13 = GetMixedEvent()->EventIndex(aodParticle->GetTrackLabel(0)) ;
  }
  
  Double_t phiTrig = aodParticle->Phi();
  Double_t etaTrig = aodParticle->Eta(); 
  Double_t ptTrig  = aodParticle->Pt();  
  
  Double_t pt             = -100. ;
  Double_t px             = -100. ;
  Double_t py             = -100. ;
  Double_t rat            = -100. ; 
  Double_t xE             = -100. ; 
  Double_t cosi           = -100. ; 
  Double_t phi            = -100. ;
  Double_t eta            = -100. ;
  Double_t pout           = -100. ;
  
  Double_t ptDecay1       = 0. ;
  Double_t pxDecay1       = 0. ;
  Double_t pyDecay1       = 0. ;
  Double_t phiDecay1      = 0. ;
  Double_t ptDecay2       = 0. ;
  Double_t pxDecay2       = 0. ;
  Double_t pyDecay2       = 0. ;
  Double_t phiDecay2      = 0. ;
  
  Double_t ratDecay1      = -100. ;  
  Double_t ratDecay2      = -100. ; 
  Double_t deltaPhi       = -100. ;
  Double_t deltaPhiOrg    = -100. ;
  Double_t deltaPhiDecay1 = -100. ;
  Double_t deltaPhiDecay2 = -100. ;
  
  TVector3 p3;  
  TLorentzVector photonMom ;	
  TObjArray * clusters  = 0x0 ;  
  TObjArray * reftracks = 0x0;
  Int_t nrefs           = 0;
  Int_t nTracks         = GetCTSTracks()->GetEntriesFast() ;
  
  if(fPi0Trigger){
    indexPhoton1 = aodParticle->GetCaloLabel (0);
    indexPhoton2 = aodParticle->GetCaloLabel (1);
    if(GetDebug() > 1)printf("indexPhoton1 = %d, indexPhoton2 = %d \n", indexPhoton1, indexPhoton2);
    
    if(indexPhoton1!=-1 && indexPhoton2!=-1){
      if(aodParticle->GetDetector()=="EMCAL") clusters = GetEMCALClusters() ;
      else                                    clusters = GetPHOSClusters()  ;
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
  for(Int_t ipr = 0;ipr < pl->GetEntriesFast() ; ipr ++ ){
    AliVTrack * track = (AliVTrack *) (pl->At(ipr)) ;
    
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
       track->GetID() == aodParticle->GetTrackLabel(2) || track->GetID() == aodParticle->GetTrackLabel(3)   ) 
      continue ;
    
    if(IsFiducialCutOn()){
      Bool_t in = GetFiducialCut()->IsInFiducialCut(mom,"CTS") ;
      if(! in ) continue ;
    }    
    
    //jump out this event if near side associated particle pt larger than trigger
    if (fMakeAbsoluteLeading){
      if(pt > ptTrig && TMath::Abs(phi-phiTrig)<TMath::PiOver2())  return kFALSE;
    }
    
    //Only for mixed event
    Int_t evtIndex2 = 0 ; 
    if (GetMixedEvent()) {
      evtIndex2 = GetMixedEvent()->EventIndex(track->GetID()) ;
      if (evtIndex11 == evtIndex2 || evtIndex12 == evtIndex2 || evtIndex13 == evtIndex2 ) // photon and track from different events
        continue ; 
      //vertex cut
      if (TMath::Abs(GetVertex(evtIndex2)[2]) > GetZvertexCut()) 
        return kFALSE;
    }    
    
    if(fPi0Trigger){
      if(indexPhoton1!=-1 && indexPhoton2!=-1){
        if(ptDecay1) ratDecay1 = pt/ptDecay1 ;
        if(ptDecay2) ratDecay2 = pt/ptDecay2 ; 
        deltaPhiDecay1 = phiDecay1-phi;
        deltaPhiDecay2 = phiDecay2-phi;
        if(deltaPhiDecay1< -TMath::PiOver2()) deltaPhiDecay1+=TMath::TwoPi();
        if(deltaPhiDecay1>3*TMath::PiOver2()) deltaPhiDecay1-=TMath::TwoPi();
        if(deltaPhiDecay2< -TMath::PiOver2()) deltaPhiDecay2+=TMath::TwoPi();
        if(deltaPhiDecay2>3*TMath::PiOver2()) deltaPhiDecay2-=TMath::TwoPi();    
      }
    } //do decay-hadron correlation    
    
    //Selection within angular range
    deltaPhi    = phiTrig-phi;
    deltaPhiOrg = deltaPhi;
    if(deltaPhi <= -TMath::PiOver2()) deltaPhi+=TMath::TwoPi();
    if(deltaPhi > 3*TMath::PiOver2()) deltaPhi-=TMath::TwoPi();
    
    pout = pt*TMath::Sin(deltaPhi) ;
    rat  = pt/ptTrig ;
    xE   =-pt/ptTrig*TMath::Cos(deltaPhi);
    cosi =-100;
    if(xE > 0 ) cosi = TMath::Log(1/xE); 
    else {
      if(GetDebug() > 1 ) printf("AliAnaParticleHadronCorrelation::MakeChargedCorrelation() - xE=%f negative or null, check!, pT %f, ptTrig %f\n", xE,pt,ptTrig);
      //return kFALSE; // shbould return or continue?
    }
    
    if(GetDebug() > 2)
      printf("AliAnaParticleHadronCorrelation::MakeChargedCorrelation() - Charged hadron: pt %f, phi %f, phi trigger %f. Cuts:  delta phi  %2.2f < %2.2f < %2.2f, pT min %2.2f \n",
             pt,phi, phiTrig,fDeltaPhiMinCut, deltaPhi, fDeltaPhiMaxCut, GetMinPt());
    
    // Fill Histograms
    if(bFillHisto){
      
      Int_t    assocBin   = -1; 
      for(Int_t i = 0 ; i < fNAssocPtBins ; i++){
        if(ptTrig > fAssocPtBinLimit[i] && ptTrig < fAssocPtBinLimit[i+1]) assocBin= i; 
      }
      
      fhAssocPt->Fill(ptTrig,pt);
      
      if(TMath::Cos(deltaPhi) < 0 && assocBin >= 0 )//away side 
        fhXEAssocPtBin[assocBin]->Fill(ptTrig, xE) ;
      
      //Hardcoded values, BAD, FIXME
      Double_t  dphiBrad = atan2(sin(deltaPhiOrg), cos(deltaPhiOrg))/TMath::Pi();//-1 to 1
      if(TMath::Abs(dphiBrad)>0.325 && TMath::Abs(dphiBrad)<0.475){
        fhAssocPtBkg->Fill(ptTrig, pt);
      }
      
      if(dphiBrad<-1./3) dphiBrad += 2;
      if(assocBin>=0){
        fhDeltaPhiBradAssocPtBin[assocBin]->Fill(ptTrig, dphiBrad);
        fhDeltaPhiAssocPtBin    [assocBin]->Fill(ptTrig, deltaPhi);
      }
      
      fhEtaCharged     ->Fill(pt,eta);
      fhPhiCharged     ->Fill(pt,phi);
      fhDeltaEtaCharged->Fill(ptTrig,aodParticle->Eta()-eta);
      fhDeltaPhiCharged->Fill(ptTrig, deltaPhi);
      fhDeltaPhiDeltaEtaCharged->Fill(deltaPhi,aodParticle->Eta()-eta);
      
      if(GetDebug() > 2 ) printf("AliAnaParticleHadronCorrelation::MakeChargedCorrelation() - Selected charge for momentum imbalance: pt %2.2f, phi %2.2f, eta %2.2f \n",pt,phi,eta);
      //fill different multiplicity histogram
      if(DoEventSelect()){
        for(Int_t im=0; im<GetMultiBin(); im++){
          if(nTracks < ( GetMaxMulti() - GetMinMulti() )/GetMultiBin()*(im+1)){
            fhTrigDeltaPhiCharged[im]->Fill(ptTrig,deltaPhi);
            fhTrigDeltaEtaCharged[im]->Fill(ptTrig,aodParticle->Eta()-eta);
          }
        }
      }
      //delta phi cut for correlation
      if( (deltaPhi > fDeltaPhiMinCut) && ( deltaPhi < fDeltaPhiMaxCut) ) {
        fhDeltaPhiChargedPt->Fill(pt,deltaPhi);
        fhPtImbalanceCharged->Fill(ptTrig,xE); 
        fhPtHbpCharged ->Fill(ptTrig, cosi);
        fhPtTrigPout   ->Fill(ptTrig, pout) ;
        fhPtTrigCharged->Fill(ptTrig, pt) ;
        if(track->Charge() > 0) fhPtImbalancePosCharged->Fill(ptTrig,xE) ;
        else                    fhPtImbalanceNegCharged->Fill(ptTrig,xE) ;
        //fill different multiplicity histogram
        if(DoEventSelect()){
          for(Int_t im=0; im<GetMultiBin(); im++){
            if(nTracks < ( GetMaxMulti() - GetMinMulti() )/GetMultiBin()*(im+1))
              fhTrigCorr[im]->Fill(ptTrig,xE);
          }
        } //multiplicity events selection
      } //delta phi cut for correlation
      else if ((deltaPhi > fUeDeltaPhiMinCut) && ( deltaPhi < fUeDeltaPhiMaxCut)) { //UE study
        fhDeltaPhiUeChargedPt->Fill(pt,deltaPhi);
        Double_t randomphi = gRandom->Uniform(TMath::Pi()/2,3*TMath::Pi()/2);
        Double_t uexE = -(pt/ptTrig)*TMath::Cos(randomphi);
        if(uexE < 0.) uexE = -uexE;
        if(GetDebug() > 1)printf("AliAnaParticleHadronCorrelation::MakeChargedCorrelation() - xe = %f, uexE = %f \n", xE, uexE);
        fhPtImbalanceUeCharged->Fill(ptTrig,uexE);
        if(uexE>0)fhPtHbpUeCharged->Fill(ptTrig,TMath::Log(1/uexE));
        if(DoEventSelect()){
          for(Int_t im=0; im<GetMultiBin(); im++){
            if(nTracks < ( GetMaxMulti() - GetMinMulti() )/GetMultiBin()*(im+1))
              fhTrigUeCorr[im]->Fill(ptTrig,xE);
          }
        } //multiplicity events selection
        
      } //UE study
      
      if(fPi0Trigger){
        if(indexPhoton1!=-1 && indexPhoton2!=-1){
          fhDeltaPhiDecayCharged->Fill(ptDecay1, deltaPhiDecay1);
          fhDeltaPhiDecayCharged->Fill(ptDecay2, deltaPhiDecay2);
          if(GetDebug() > 1)printf("AliAnaParticleHadronCorrelation::MakeChargedCorrelation() - deltaPhoton1 = %f, deltaPhoton2 = %f \n", deltaPhiDecay1, deltaPhiDecay2);
          if( (deltaPhiDecay1 > fDeltaPhiMinCut) && ( deltaPhiDecay1 < fDeltaPhiMaxCut) )
            fhPtImbalanceDecayCharged->Fill(ptDecay1,ratDecay1); 
          if( (deltaPhiDecay2 > fDeltaPhiMinCut) && ( deltaPhiDecay2 < fDeltaPhiMaxCut) )
            fhPtImbalanceDecayCharged->Fill(ptDecay2,ratDecay2);
          if(GetDebug() > 1)printf("AliAnaParticleHadronCorrelation::MakeChargedCorrelation() - ratPhoton1 = %f, ratPhoton2 = %f \n", pt/ptDecay1, pt/ptDecay2);
        } //index of decay photons found
      } //make decay-hadron correlation          
      
      //several UE calculation 
      if(fMakeSeveralUE){
        if((deltaPhi<-fUeDeltaPhiMinCut) && (deltaPhi >-fUeDeltaPhiMaxCut)){  
          fhDeltaPhiUeLeftCharged->Fill(pt,deltaPhi);
          fhPtImbalanceUeLeftCharged->Fill(ptTrig,rat);
          fhPtHbpUeLeftCharged->Fill(ptTrig,cosi);
        }
        if((deltaPhi>fUeDeltaPhiMinCut) && (deltaPhi <fUeDeltaPhiMaxCut)){  
          fhDeltaPhiUeRightCharged->Fill(pt,deltaPhi);
          fhPtImbalanceUeRightCharged->Fill(ptTrig,rat);
          fhPtHbpUeRightCharged->Fill(ptTrig,cosi);
          
        }
      } //several UE calculation
      
      //Fill leading particle histogram   
      fhPtLeading->Fill(ptTrig);
      if(phiTrig<0)phiTrig+=TMath::TwoPi();
      fhPhiLeading->Fill(ptTrig, phiTrig);
      fhEtaLeading->Fill(ptTrig, etaTrig);
      
    } //Fill histogram 
    else{
      nrefs++;
      if(nrefs==1){
        reftracks = new TObjArray(0);
        TString trackname = Form("%s+Tracks", GetAODObjArrayName().Data());
        reftracks->SetName(trackname.Data());
        reftracks->SetOwner(kFALSE);
      }
      reftracks->Add(track);
    }//aod particle loop
  }// track loop
  
  //Fill AOD with reference tracks, if not filling histograms
  if(!bFillHisto && reftracks) {
    aodParticle->AddObjArray(reftracks);
  }
  
  return kTRUE;
  
}  

//________________________________________________________________________________________________________________
Bool_t  AliAnaParticleHadronCorrelation::MakeNeutralCorrelation(AliAODPWG4ParticleCorrelation * const aodParticle, 
                                                                TObjArray* pi0list, const Bool_t bFillHisto)  
{  
  // Neutral Pion Correlation Analysis
  if(GetDebug() > 1) printf("AliAnaParticleHadronCorrelation::MakeNeutralCorrelation() - Make trigger particle - pi0 correlation, %d pi0's \n",pi0list->GetEntriesFast());
  
  Int_t evtIndex11 = 0 ; 
  Int_t evtIndex12 = 0 ; 
  if (GetMixedEvent()) {
    evtIndex11 = GetMixedEvent()->EventIndexForCaloCluster(aodParticle->GetCaloLabel(0)) ;
    evtIndex12 = GetMixedEvent()->EventIndexForCaloCluster(aodParticle->GetCaloLabel(1)) ;    
  }  
  
  Double_t pt   = -100. ;
  Double_t px   = -100. ;
  Double_t py   = -100. ;
  Double_t rat  = -100. ; 
  Double_t phi  = -100. ;
  Double_t eta  = -100. ;
  Double_t xE   = -100. ; 
  Double_t cosi = -100. ; 
  
  Double_t ptTrig  = aodParticle->Pt();
  Double_t phiTrig = aodParticle->Phi();
  Double_t etaTrig = aodParticle->Eta();
  Double_t pxTrig  = aodParticle->Px();
  Double_t pyTrig  = aodParticle->Py();
  
  Int_t indexPhoton1 =-1  ;
  Int_t indexPhoton2 =-1  ;    
  Double_t ptDecay1  = 0. ;
  Double_t pxDecay1  = 0. ;
  Double_t pyDecay1  = 0. ;
  Double_t phiDecay1 = 0. ;
  Double_t ptDecay2  = 0. ;
  Double_t pxDecay2  = 0. ;
  Double_t pyDecay2  = 0. ;
  Double_t phiDecay2 = 0. ;
  
  Double_t ratDecay1      = -100. ;  
  Double_t ratDecay2      = -100. ; 
  Double_t deltaPhi       = -100. ;
  Double_t deltaPhiDecay1 = -100. ;
  Double_t deltaPhiDecay2 = -100. ;
  
  TObjArray * clusters = 0x0 ;  
  TLorentzVector photonMom ;
	
  if(fPi0Trigger){
    indexPhoton1 = aodParticle->GetCaloLabel (0);
    indexPhoton2 = aodParticle->GetCaloLabel (1);
    if(GetDebug() > 1)
      printf("AliAnaParticleHadronCorrelation::MakeNeutralCorrelationFillHistograms() - indexPhoton1 = %d, indexPhoton2 = %d \n", indexPhoton1, indexPhoton2);
    
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
        if(GetDebug() > 1)
          printf("AliAnaParticleHadronCorrelation::MakeNeutralCorrelationFillHistograms() - Photon1 = %f, Photon2 = %f \n", ptDecay1, ptDecay2);
      } //photonAOD loop        
    } //index of decay photons found
    if(ptTrig && bFillHisto) fhPtPi0DecayRatio->Fill(ptTrig, ptDecay1/ptTrig, ptDecay2/ptTrig);
  } //make decay-hadron correlation
  
  TObjArray * refpi0    =0x0;
  Int_t nrefs = 0;
  
  //Loop on stored AOD pi0
  Int_t naod = pi0list->GetEntriesFast();
  if(GetDebug() > 0) 
    printf("AliAnaParticleHadronCorrelation::MakeNeutralCorrelationFillHistograms() -  aod branch entries %d\n", naod);
  for(Int_t iaod = 0; iaod < naod ; iaod++){
    AliAODPWG4Particle* pi0 =  (AliAODPWG4Particle*) (pi0list->At(iaod));
    
    Int_t evtIndex2 = 0 ; 
    Int_t evtIndex3 = 0 ; 
    if (GetMixedEvent()) {
      evtIndex2 = GetMixedEvent()->EventIndexForCaloCluster(pi0->GetCaloLabel(0)) ;
      evtIndex3 = GetMixedEvent()->EventIndexForCaloCluster(pi0->GetCaloLabel(1)) ;
      
      if (evtIndex11 == evtIndex2 || evtIndex12 == evtIndex2 || 
          evtIndex11 == evtIndex3 || evtIndex12 == evtIndex3) // trigger and pi0 are not from different events
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
    //Float_t deltaPhi = TMath::Abs(phiTrig-phi);
    //if( (deltaPhi < fDeltaPhiMinCut) || ( deltaPhi > fDeltaPhiMaxCut) ) continue ;
    
    if(bFillHisto){
      
      deltaPhi = phiTrig-phi;
      if(deltaPhi<-TMath::PiOver2()) deltaPhi+=TMath::TwoPi();
      if(deltaPhi>3*TMath::PiOver2()) deltaPhi-=TMath::TwoPi();
      
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
          deltaPhiDecay1 = phiDecay1-phi;
          deltaPhiDecay2 = phiDecay2-phi;
          if(deltaPhiDecay1< -TMath::PiOver2()) deltaPhiDecay1+=TMath::TwoPi();
          if(deltaPhiDecay1>3*TMath::PiOver2()) deltaPhiDecay1-=TMath::TwoPi();
          if(deltaPhiDecay2< -TMath::PiOver2()) deltaPhiDecay2+=TMath::TwoPi();
          if(deltaPhiDecay2>3*TMath::PiOver2()) deltaPhiDecay2-=TMath::TwoPi();   
          fhDeltaPhiDecayNeutral->Fill(ptDecay1, deltaPhiDecay1);
          fhDeltaPhiDecayNeutral->Fill(ptDecay2, deltaPhiDecay2);
          if(GetDebug() > 1)printf("AliAnaParticleHadronCorrelation::MakeNeutralCorrelation() - deltaPhoton1 = %f, deltaPhoton2 = %f \n", deltaPhiDecay1, deltaPhiDecay2);
          if( (deltaPhiDecay1 > fDeltaPhiMinCut) && ( deltaPhiDecay1 < fDeltaPhiMaxCut) )
            fhPtImbalanceDecayNeutral->Fill(ptDecay1,ratDecay1); 
          if( (deltaPhiDecay2 > fDeltaPhiMinCut) && ( deltaPhiDecay2 < fDeltaPhiMaxCut) )
            fhPtImbalanceDecayNeutral->Fill(ptDecay2,ratDecay2);
          if(GetDebug() > 1)printf("AliAnaParticleHadronCorrelation::MakeNeutralCorrelation() - ratPhoton1 = %f, ratPhoton2 = %f \n", pt/ptDecay1, pt/ptDecay2);
        }
      } //do decay-hadron correlation
      
      fhEtaNeutral->Fill(pt,eta);
      fhPhiNeutral->Fill(pt,phi);
      fhDeltaEtaNeutral->Fill(ptTrig,etaTrig-eta);
      fhDeltaPhiNeutral->Fill(ptTrig,deltaPhi);
      fhDeltaPhiDeltaEtaNeutral->Fill(deltaPhi,etaTrig-eta);
      
      //delta phi cut for correlation
      if( (deltaPhi > fDeltaPhiMinCut) && ( deltaPhi < fDeltaPhiMaxCut) ) {
        fhDeltaPhiNeutralPt->Fill(pt,deltaPhi);
        fhPtImbalanceNeutral->Fill(ptTrig,rat); 
        fhPtHbpNeutral->Fill(ptTrig,cosi); 
      }
      else {
        fhDeltaPhiUeNeutralPt->Fill(pt,deltaPhi);
        fhPtImbalanceUeNeutral->Fill(ptTrig,rat);
        fhPtHbpUeNeutral->Fill(ptTrig,cosi); 
      }
      //several UE calculation 
      if(fMakeSeveralUE){
        if((deltaPhi<-fUeDeltaPhiMinCut) && (deltaPhi >-fUeDeltaPhiMaxCut)){  
          fhDeltaPhiUeLeftNeutral->Fill(pt,deltaPhi);
          fhPtImbalanceUeLeftNeutral->Fill(ptTrig,rat);
          fhPtHbpUeLeftNeutral->Fill(ptTrig,cosi);
        }
        if((deltaPhi>fUeDeltaPhiMinCut) && (deltaPhi <fUeDeltaPhiMaxCut)){  
          fhDeltaPhiUeRightNeutral->Fill(pt,deltaPhi);
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
  Double_t pxprim  = 0 ;
  Double_t pyprim  = 0 ;
  Double_t pzprim  = 0 ;
  Int_t    nTracks = 0 ;  
  Int_t iParticle  = 0 ;
  Double_t charge  = 0.;
  
  Double_t mcrat   =-100 ;
  Double_t mcxE    =-100 ;
  Double_t mccosi  =-100 ;
  
  //Track loop, select tracks with good pt, phi and fill AODs or histograms
  //Int_t currentIndex = -1 ; 
  Double_t mcTrackPt  = 0 ;
  Double_t mcTrackPhi = 0 ;
  Double_t mcTrackEta = 0 ;
  Double_t mcTrackPx  = 0 ;
  Double_t mcTrackPy  = 0 ;
  Double_t mcTrackPz  = 0 ;
  
  if(GetReader()->ReadStack()){
    nTracks = GetMCStack()->GetNtrack() ;
  }
  else nTracks = GetReader()->GetAODMCParticles()->GetEntriesFast() ;
  //Int_t trackIndex[nTracks];
  
  Int_t label= aodParticle->GetLabel();
  if(label<0){
    printf("AliAnaParticleHadronCorrelation::MakeMCChargedCorrelation *** bad label ***:  label %d \n", label);
    return;
  }  
  
  if(GetReader()->ReadStack()){
    stack =  GetMCStack() ;
    if(!stack) {
      printf(" AliAnaParticleHadronCorrelation::MakeMCChargedCorrelation- Stack not available, is the MC handler called? STOP\n");
      abort();
    }
    
    nTracks=stack->GetNprimary();
    if(label >=  stack->GetNtrack()) {
      if(GetDebug() > 2)  printf("AliAnaParticleHadronCorrelation::MakeMCChargedCorrelation *** large label ***:  label %d, n tracks %d \n", label, stack->GetNtrack());
      return ;
    }
    primary = stack->Particle(label);
    if(!primary){
      printf("AliAnaParticleHadronCorrelation::MakeMCChargedCorrelation *** no primary ***:  label %d \n", label);   
      return;
    }
    
    eprim    = primary->Energy();
    ptprim   = primary->Pt();
    phiprim  = primary->Phi();
    etaprim  = primary->Eta();
    pxprim   = primary->Px();
    pyprim   = primary->Py();
    pzprim   = primary->Pz(); 
    
    if(primary){
      
      for (iParticle = 0 ; iParticle <  nTracks ; iParticle++) {
        TParticle * particle = stack->Particle(iParticle);
        TLorentzVector momentum;
        //keep only final state particles
        if(particle->GetStatusCode()!=1) continue ;
        Int_t pdg = particle->GetPdgCode();						
        charge = TDatabasePDG::Instance()->GetParticle(pdg)->Charge();
        particle->Momentum(momentum);
        
        //---------- Charged particles ----------------------
        if(charge != 0){   
          //Particles in CTS acceptance
          Bool_t inCTS =  GetFiducialCut()->IsInFiducialCut(momentum,"CTS");
          if(TMath::Abs(pdg) == 11 && stack->Particle(particle->GetFirstMother())->GetPdgCode()==22) continue ;
          if(inCTS&&momentum.Pt() >GetMinPt())
          {            
            mcTrackPt  = particle->Pt();
            mcTrackPhi = particle->Phi();
            mcTrackEta = particle->Eta();
            mcTrackPx  = particle->Px();
            mcTrackPy  = particle->Py();
            mcTrackPz  = particle->Pz();              
            if(mcTrackPhi < 0) mcTrackPhi+=TMath::TwoPi();              
            //Select only hadrons in pt range
            if(mcTrackPt < GetMinPt() || mcTrackPt > GetMaxPt()) continue ;
            //remove trigger itself for correlation when use charged triggers 
            if(label==iParticle && mcTrackPt==ptprim && mcTrackPhi==phiprim && mcTrackEta==etaprim) 
              continue ;                  
            //jumped out this event if near side associated partile pt larger than trigger
            if( mcTrackPt> ptprim && TMath::Abs(mcTrackPhi-phiprim)<TMath::PiOver2()) 
              return ;
            
            mcrat   = mcTrackPt/ptprim ;
            mcxE    = -(mcTrackPx*pxprim+mcTrackPy*pyprim)/(ptprim*ptprim);
            if(mcxE <0.) mcxE =-mcxE;
            mccosi = TMath::Log(1/mcxE);
            // printf("rat = %f, xE = %f, cosi =%f \n", rat, xE, cosi);
            // printf("phi = %f \n", phi);
            
            //Selection within angular range
            Double_t mcdeltaPhi = phiprim-mcTrackPhi;
            if( mcdeltaPhi< -TMath::PiOver2())  mcdeltaPhi+=TMath::TwoPi();
            if( mcdeltaPhi>3*TMath::PiOver2())  mcdeltaPhi-=TMath::TwoPi();              
            Double_t mcpout = mcTrackPt*TMath::Sin(mcdeltaPhi) ;              
            if(GetDebug()>0 )  
              printf("AliAnaParticleHadronCorrelation::MakeMCChargedCorrelation() - Charged hadron: track Pt %f, track Phi %f, phi trigger %f. Cuts:  delta phi  %2.2f < %2.2f < %2.2f, pT min %2.2f \n",
                     mcTrackPt,mcTrackPhi, phiprim,fDeltaPhiMinCut, mcdeltaPhi, fDeltaPhiMaxCut, GetMinPt());              
            // Fill Histograms
            fhMCEtaCharged->Fill(mcTrackPt,mcTrackEta);
            fhMCPhiCharged->Fill(mcTrackPt,mcTrackPhi);
            fhMCDeltaEtaCharged->Fill(ptprim,etaprim-mcTrackEta);
            fhMCDeltaPhiCharged->Fill(ptprim,mcdeltaPhi);
            fhMCPtAssocDeltaPhi->Fill(mcTrackPt, mcdeltaPhi);
            //  fhDeltaPhiCharged->Fill(ptTrig, deltaPhi);
            fhMCDeltaPhiDeltaEtaCharged->Fill(mcdeltaPhi,etaprim-mcTrackEta);
            
            //delta phi cut for correlation
            if( (mcdeltaPhi > fDeltaPhiMinCut) && ( mcdeltaPhi < fDeltaPhiMaxCut) ) {
              fhMCDeltaPhiChargedPt->Fill(mcTrackPt,mcdeltaPhi);
              fhMCPtImbalanceCharged->Fill(ptprim,mcxE); 
              fhMCPtHbpCharged->Fill(ptprim,mccosi);
              fhMCPtTrigPout->Fill(ptprim, mcpout) ;
            }//delta phi cut for correlation
          } //tracks after cuts
        }//Charged
      } //track loop
    } //when the leading particles could trace back to MC
  } //ESD MC
  else if(GetReader()->ReadAODMCParticles()){
    //Get the list of MC particles
    mcparticles0 = GetReader()->GetAODMCParticles(0);
    if(!mcparticles0) return;
    if(label >=mcparticles0->GetEntriesFast()){
      if(GetDebug() > 2)  
        printf("AliAnaParticleHadronCorrelation::MakeMCChargedCorrelation *** large label ***:  label %d, n tracks %d \n", label,mcparticles0->GetEntriesFast());
      return;
    }
    //Get the particle
    aodprimary = (AliAODMCParticle*) mcparticles0->At(label);
    if(!aodprimary)  {
      printf("AliAnaParticleHadronCorrelation::MakeMCChargedCorrelation *** no AOD primary ***:  label %d \n", label);   
      return;
    }
    
    ptprim  = aodprimary->Pt();
    phiprim  = aodprimary->Phi();
    etaprim  =aodprimary->Eta();
    pxprim =aodprimary->Px();
    pyprim =aodprimary->Py();
    pzprim =aodprimary->Pz();  
    if(aodprimary){
      mcparticles= GetReader()->GetAODMCParticles();
      for (Int_t i=0; i<nTracks;i++) {
        AliAODMCParticle *part = (AliAODMCParticle*) mcparticles->At(i);
        if (!part->IsPhysicalPrimary()) continue;        
        Int_t pdg = part->GetPdgCode();	
        charge = TDatabasePDG::Instance()->GetParticle(pdg)->Charge();
        TLorentzVector momentum(part->Px(),part->Py(),part->Pz(),part->E());        
        if(charge != 0){
          if(part->Pt()> GetReader()->GetCTSPtMin()){
            //Particles in CTS acceptance
            Bool_t inCTS =  GetFiducialCut()->IsInFiducialCut(momentum,"CTS");
            Int_t indexmother=part->GetMother();
            if(indexmother>-1)
            {
              Int_t mPdg = ((AliAODMCParticle*) mcparticles->At(indexmother)) ->GetPdgCode();
              if(TMath::Abs(pdg) == 11 && mPdg == 22) continue;
            }
            if(inCTS&&momentum.Pt() >GetMinPt())
            {            
              mcTrackPt=part->Pt();
              mcTrackPhi=part->Phi();
              mcTrackEta=part->Eta();
              mcTrackPx=part->Px();
              mcTrackPy=part->Py();
              mcTrackPz=part->Pz();              
              if(mcTrackPhi < 0) mcTrackPhi+=TMath::TwoPi();              
              //Select only hadrons in pt range
              if(mcTrackPt < GetMinPt() || mcTrackPt > GetMaxPt()) continue ;
              //remove trigger itself for correlation when use charged triggers 
              if(label==i && mcTrackPt==ptprim && mcTrackPhi==phiprim && mcTrackEta==etaprim) 
                continue ;                  
              //jumped out this event if near side associated partile pt larger than trigger
              if( mcTrackPt> ptprim && TMath::Abs(mcTrackPhi-phiprim)<TMath::PiOver2()) 
                return ;
              
              mcrat   = mcTrackPt/ptprim ;
              mcxE    = -(mcTrackPx*pxprim+mcTrackPy*pyprim)/(ptprim*ptprim);
              if(mcxE <0.)mcxE =-mcxE;
              mccosi = TMath::Log(1/mcxE);
              
              //Selection within angular range
              Double_t mcdeltaPhi = phiprim-mcTrackPhi;
              if( mcdeltaPhi< -TMath::PiOver2())  mcdeltaPhi+=TMath::TwoPi();
              if( mcdeltaPhi>3*TMath::PiOver2())  mcdeltaPhi-=TMath::TwoPi();              
              Double_t mcpout = mcTrackPt*TMath::Sin(mcdeltaPhi) ;              
              if(GetDebug()>0)  
                printf("AliAnaParticleHadronCorrelation::MakeMCChargedCorrelation() - Charged hadron: track Pt %f, track Phi %f, phi trigger %f. Cuts:  delta phi  %2.2f < %2.2f < %2.2f, pT min %2.2f \n",
                       mcTrackPt,mcTrackPhi, phiprim,fDeltaPhiMinCut, mcdeltaPhi, fDeltaPhiMaxCut, GetMinPt());              
              // Fill Histograms
              fhMCEtaCharged->Fill(mcTrackPt,mcTrackEta);
              fhMCPhiCharged->Fill(mcTrackPt,mcTrackPhi);
              fhMCDeltaEtaCharged->Fill(ptprim,etaprim-mcTrackEta);
              fhMCDeltaPhiCharged->Fill(ptprim,mcdeltaPhi);
              fhMCPtAssocDeltaPhi->Fill(mcTrackPt, mcdeltaPhi);
              //  fhDeltaPhiCharged->Fill(ptTrig, deltaPhi);
              fhMCDeltaPhiDeltaEtaCharged->Fill(mcdeltaPhi,etaprim-mcTrackEta);
              
              //delta phi cut for correlation
              if( (mcdeltaPhi > fDeltaPhiMinCut) && ( mcdeltaPhi < fDeltaPhiMaxCut) ) {
                fhMCDeltaPhiChargedPt->Fill(mcTrackPt,mcdeltaPhi);
                fhMCPtImbalanceCharged->Fill(ptprim,mcxE); 
                fhMCPtHbpCharged->Fill(ptprim,mccosi);
                fhMCPtTrigPout->Fill(ptprim, mcpout) ;
              }//delta phi cut for correlation
              
            } //tracks after cuts
            
          } //with minimum pt cut
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
  AliAnaPartCorrBaseClass::Print(" ");
  
  printf("Phi trigger particle-Hadron      <     %3.2f\n", fDeltaPhiMaxCut) ; 
  printf("Phi trigger particle-Hadron      >     %3.2f\n", fDeltaPhiMinCut) ;
  printf("Isolated Trigger?  %d\n", fSelectIsolated) ;
  printf("Phi trigger particle-UeHadron    <    %3.2f\n", fUeDeltaPhiMaxCut) ; 
  printf("Phi trigger particle-UeHadron    >    %3.2f\n", fUeDeltaPhiMinCut) ;
  printf("Several UE?  %d\n", fMakeSeveralUE) ;
  printf("Name of AOD Pi0 Branch %s \n",fPi0AODBranchName.Data());
  printf("Do Decay-hadron correlation ?  %d\n", fPi0Trigger) ;
  printf("Select absolute leading for cluster triggers ?  %d\n", fMakeAbsoluteLeading) ;
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

