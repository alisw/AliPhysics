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
//////////////////////////////////////////////////////////////////////////////


// --- ROOT system ---
#include "TH2F.h"
#include "TClonesArray.h"
#include "TClass.h"

//---- ANALYSIS system ----
#include "AliNeutralMesonSelection.h" 
#include "AliAnaParticleHadronCorrelation.h" 
#include "AliCaloTrackReader.h"
#include "AliCaloPID.h"
#include "AliAODPWG4ParticleCorrelation.h"
#include "AliFidutialCut.h"
#include "AliAODTrack.h"
#include "AliAODCaloCluster.h"
#include "AliMCAnalysisUtils.h"
#include "TParticle.h"
#include "AliStack.h"

ClassImp(AliAnaParticleHadronCorrelation)


//____________________________________________________________________________
  AliAnaParticleHadronCorrelation::AliAnaParticleHadronCorrelation() : 
    AliAnaPartCorrBaseClass(),
    fDeltaPhiMaxCut(0.), fDeltaPhiMinCut(0.), fSelectIsolated(0),  
    fhPhiCharged(0), fhPhiNeutral(0), fhEtaCharged(0), fhEtaNeutral(0), 
    fhDeltaPhiCharged(0), fhDeltaPhiNeutral(0), 
    fhDeltaEtaCharged(0), fhDeltaEtaNeutral(0),
    fhDeltaPhiChargedPt(0), fhDeltaPhiNeutralPt(0), 
    fhPtImbalanceNeutral(0), fhPtImbalanceCharged(0)
{
  //Default Ctor
  
  //Initialize parameters
  InitParameters();
}

//____________________________________________________________________________
AliAnaParticleHadronCorrelation::AliAnaParticleHadronCorrelation(const AliAnaParticleHadronCorrelation & g) :   
  AliAnaPartCorrBaseClass(g),
  fDeltaPhiMaxCut(g.fDeltaPhiMaxCut), fDeltaPhiMinCut(g.fDeltaPhiMinCut), 
  fSelectIsolated(g.fSelectIsolated),
  fhPhiCharged(g.fhPhiCharged), fhPhiNeutral(g.fhPhiNeutral), 
  fhEtaCharged(g.fhEtaCharged), fhEtaNeutral(g.fhEtaNeutral), 
  fhDeltaPhiCharged(g.fhDeltaPhiCharged),  
  fhDeltaPhiNeutral(g.fhDeltaPhiNeutral), 
  fhDeltaEtaCharged(g.fhDeltaEtaCharged), 
  fhDeltaEtaNeutral(g.fhDeltaEtaNeutral), 
  fhDeltaPhiChargedPt(g.fhDeltaPhiChargedPt), 
  fhDeltaPhiNeutralPt(g.fhDeltaPhiNeutralPt), 
  fhPtImbalanceNeutral(g.fhPtImbalanceNeutral), 
  fhPtImbalanceCharged(g.fhPtImbalanceCharged)
{
  // cpy ctor
  
}

//_________________________________________________________________________
AliAnaParticleHadronCorrelation & AliAnaParticleHadronCorrelation::operator = (const AliAnaParticleHadronCorrelation & source)
{
  // assignment operator
  
  if(this == &source)return *this;
  ((AliAnaPartCorrBaseClass *)this)->operator=(source);
  
  fDeltaPhiMaxCut = source.fDeltaPhiMaxCut ; 
  fDeltaPhiMinCut = source.fDeltaPhiMinCut ; 
  fSelectIsolated = source.fSelectIsolated ;
  
  fhPhiCharged = source.fhPhiCharged ; fhPhiNeutral = source.fhPhiNeutral ; 
  fhEtaCharged = source.fhEtaCharged ; fhEtaNeutral = source.fhEtaNeutral ; 
  fhDeltaPhiCharged = source.fhDeltaPhiCharged ;  
  fhDeltaPhiNeutral = source.fhDeltaPhiNeutral ; 
  fhDeltaPhiNeutralPt = source.fhDeltaPhiNeutralPt ; 
  fhDeltaEtaCharged = source.fhDeltaEtaCharged ; 
  fhDeltaEtaNeutral = source.fhDeltaEtaNeutral ; 
  fhDeltaPhiChargedPt = source.fhDeltaPhiChargedPt ;
  
  fhPtImbalanceNeutral = source.fhPtImbalanceNeutral ; 
  fhPtImbalanceCharged = source.fhPtImbalanceCharged ; 

  return *this;

}

//________________________________________________________________________
TList *  AliAnaParticleHadronCorrelation::GetCreateOutputObjects()
{  
  
  // Create histograms to be saved in output file and 
  // store them in fOutputContainer
  TList * outputContainer = new TList() ; 
  outputContainer->SetName("CorrelationHistos") ; 
  
  Int_t nptbins  = GetHistoNPtBins();
  Int_t nphibins = GetHistoNPhiBins();
  Int_t netabins = GetHistoNEtaBins();
  Float_t ptmax  = GetHistoPtMax();
  Float_t phimax = GetHistoPhiMax();
  Float_t etamax = GetHistoEtaMax();
  Float_t ptmin  = GetHistoPtMin();
  Float_t phimin = GetHistoPhiMin();
  Float_t etamin = GetHistoEtaMin();	
  
  //Correlation with charged hadrons
  if(GetReader()->IsCTSSwitchedOn()) {
    fhPhiCharged  = new TH2F
      ("PhiCharged","#phi_{h^{#pm}}  vs p_{T trigger}",
       nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
    fhPhiCharged->SetYTitle("#phi_{h^{#pm}} (rad)");
    fhPhiCharged->SetXTitle("p_{T trigger} (GeV/c)");
    
    fhEtaCharged  = new TH2F
      ("EtaCharged","#eta_{h^{#pm}}  vs p_{T trigger}",
       nptbins,ptmin,ptmax,netabins,etamin,etamax); 
    fhEtaCharged->SetYTitle("#eta_{h^{#pm}} (rad)");
    fhEtaCharged->SetXTitle("p_{T trigger} (GeV/c)");
    
    fhDeltaPhiCharged  = new TH2F
      ("DeltaPhiCharged","#phi_{trigger} - #phi_{h^{#pm}} vs p_{T trigger}",
       nptbins,ptmin,ptmax,200,0,6.4); 
    fhDeltaPhiCharged->SetYTitle("#Delta #phi");
    fhDeltaPhiCharged->SetXTitle("p_{T trigger} (GeV/c)");
    
    fhDeltaPhiChargedPt  = new TH2F
      ("DeltaPhiChargedPt","#phi_{trigger} - #phi_{#h^{#pm}} vs p_{T h^{#pm}}",
       nptbins,ptmin,ptmax,200,0,6.4);
    fhDeltaPhiChargedPt->SetYTitle("#Delta #phi");
    fhDeltaPhiChargedPt->SetXTitle("p_{T h^{#pm}} (GeV/c)");
    
    fhDeltaEtaCharged  = new TH2F
      ("DeltaEtaCharged","#eta_{trigger} - #eta_{h^{#pm}} vs p_{T trigger}",
       nptbins,ptmin,ptmax,200,-2,2); 
    fhDeltaEtaCharged->SetYTitle("#Delta #eta");
    fhDeltaEtaCharged->SetXTitle("p_{T trigger} (GeV/c)");
    
    fhPtImbalanceCharged  = 
      new TH2F("CorrelationCharged","z_{trigger h^{#pm}} = p_{T h^{#pm}} / p_{T trigger}",
	       nptbins,ptmin,ptmax,1000,0.,1.2); 
    fhPtImbalanceCharged->SetYTitle("z_{trigger h^{#pm}}");
    fhPtImbalanceCharged->SetXTitle("p_{T trigger}");
    
    outputContainer->Add(fhPhiCharged) ;
    outputContainer->Add(fhEtaCharged) ;
    outputContainer->Add(fhDeltaPhiCharged) ; 
    outputContainer->Add(fhDeltaEtaCharged) ;
    outputContainer->Add(fhPtImbalanceCharged) ;
    outputContainer->Add(fhDeltaPhiChargedPt) ;
  }  //Correlation with charged hadrons
  
  //Correlation with neutral hadrons
  if(GetReader()->IsEMCALSwitchedOn() || GetReader()->IsPHOSSwitchedOn()){
    
    fhPhiNeutral  = new TH2F
      ("PhiNeutral","#phi_{#pi^{0}}  vs p_{T trigger}",
       nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
    fhPhiNeutral->SetYTitle("#phi_{#pi^{0}} (rad)");
    fhPhiNeutral->SetXTitle("p_{T trigger} (GeV/c)");
    
    fhEtaNeutral  = new TH2F
      ("EtaNeutral","#eta_{#pi^{0}}  vs p_{T trigger}",
       nptbins,ptmin,ptmax,netabins,etamin,etamax); 
    fhEtaNeutral->SetYTitle("#eta_{#pi^{0}} (rad)");
    fhEtaNeutral->SetXTitle("p_{T trigger} (GeV/c)");
    
    fhDeltaPhiNeutral  = new TH2F
      ("DeltaPhiNeutral","#phi_{trigger} - #phi_{#pi^{0}} vs p_{T trigger}",
       nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
    fhDeltaPhiNeutral->SetYTitle("#Delta #phi");
    fhDeltaPhiNeutral->SetXTitle("p_{T trigger} (GeV/c)");
    
    fhDeltaPhiNeutralPt  = new TH2F
      ("DeltaPhiNeutralPt","#phi_{trigger} - #phi_{#pi^{0}} vs p_{T #pi^{0}}}",
       nptbins,ptmin,ptmax,200,0,6.4); 
    fhDeltaPhiNeutralPt->SetYTitle("#Delta #phi");
    fhDeltaPhiNeutralPt->SetXTitle("p_{T trigger} (GeV/c)");
    
    fhDeltaEtaNeutral  = new TH2F
      ("DeltaEtaNeutral","#eta_{trigger} - #eta_{#pi^{0}} vs p_{T trigger}",
       nptbins,ptmin,ptmax,200,-2,2); 
    fhDeltaEtaNeutral->SetYTitle("#Delta #eta");
    fhDeltaEtaNeutral->SetXTitle("p_{T trigger} (GeV/c)");
    
    fhPtImbalanceNeutral  = 
      new TH2F("CorrelationNeutral","z_{trigger #pi} = p_{T #pi^{0}} / p_{T trigger}",
	       nptbins,ptmin,ptmax,1000,0.,1.2); 
    fhPtImbalanceNeutral->SetYTitle("z_{trigger #pi^{0}}");
    fhPtImbalanceNeutral->SetXTitle("p_{T trigger}");
    
    outputContainer->Add(fhPhiNeutral) ;  
    outputContainer->Add(fhEtaNeutral) ;   
    outputContainer->Add(fhDeltaPhiNeutral) ; 
    outputContainer->Add(fhDeltaEtaNeutral) ; 
    outputContainer->Add(fhPtImbalanceNeutral) ;
    
    
    //Keep neutral meson selection histograms if requiered
    //Setting done in AliNeutralMesonSelection
    if(GetNeutralMesonSelection()){
      TList * nmsHistos = GetNeutralMesonSelection()->GetCreateOutputObjects() ;
      if(GetNeutralMesonSelection()->AreNeutralMesonSelectionHistosKept())
	for(Int_t i = 0; i < nmsHistos->GetEntries(); i++) outputContainer->Add(nmsHistos->At(i)) ;
    }
	
  }//Correlation with neutral hadrons
  
  return outputContainer;

}

//____________________________________________________________________________
void AliAnaParticleHadronCorrelation::InitParameters()
{
  
  //Initialize the parameters of the analysis.
  SetInputAODName("PWG4Particle");
  SetAODRefArrayName("Hadrons");  
  AddToHistogramsName("AnaHadronCorr_");

  //Correlation with neutrals
  //SetOutputAODClassName("AliAODPWG4Particle");
  //SetOutputAODName("Pi0Correlated");
	
  SetPtCutRange(2,300);
  fDeltaPhiMinCut = 1.5 ;
  fDeltaPhiMaxCut = 4.5 ;
  fSelectIsolated = kFALSE;
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
} 

//____________________________________________________________________________
void  AliAnaParticleHadronCorrelation::MakeAnalysisFillAOD()  
{  
  //Particle-Hadron Correlation Analysis, fill AODs
  
  if(!GetInputAODBranch()){
    printf("AliAnaParticleHadronCorrelation::MakeAnalysisFillAOD() - No input particles in AOD with name branch < %s >, ABORT \n",GetInputAODName().Data());
    abort();
  }
	
  if(strcmp(GetInputAODBranch()->GetClass()->GetName(), "AliAODPWG4ParticleCorrelation")){
	printf("AliAnaParticleHadronCorrelation::MakeAnalysisFillAOD() - Wrong type of AOD object, change AOD class name in input AOD: It should be <AliAODPWG4ParticleCorrelation> and not <%s> \n",GetInputAODBranch()->GetClass()->GetName());
	abort();
  }
	
  if(GetDebug() > 1){
    printf("AliAnaParticleHadronCorrelation::MakeAnalysisFillAOD() - Begin hadron correlation analysis, fill AODs \n");
    printf("AliAnaParticleHadronCorrelation::MakeAnalysisFillAOD() - In particle branch aod entries %d\n", GetInputAODBranch()->GetEntriesFast());
    printf("AliAnaParticleHadronCorrelation::MakeAnalysisFillAOD() - In CTS aod entries %d\n", GetAODCTS()->GetEntriesFast());
    printf("AliAnaParticleHadronCorrelation::MakeAnalysisFillAOD() - In EMCAL aod entries %d\n", GetAODEMCAL()->GetEntriesFast());
    printf("AliAnaParticleHadronCorrelation::MakeAnalysisFillAOD() - In PHOS aod entries %d\n", GetAODPHOS()->GetEntriesFast());
  }
  
  //Loop on stored AOD particles, trigger
  Int_t naod = GetInputAODBranch()->GetEntriesFast();
  for(Int_t iaod = 0; iaod < naod ; iaod++){
    AliAODPWG4ParticleCorrelation* particle =  (AliAODPWG4ParticleCorrelation*) (GetInputAODBranch()->At(iaod));

	//Make correlation with charged hadrons
    if(GetReader()->IsCTSSwitchedOn() )
      MakeChargedCorrelation(particle, GetAODCTS(),kFALSE);
    
    //Make correlation with neutral pions
    //Trigger particle in PHOS, correlation with EMCAL
    if(particle->GetDetector()=="PHOS" && GetReader()->IsEMCALSwitchedOn() && GetAODEMCAL()->GetEntriesFast() > 0)
      MakeNeutralCorrelationFillAOD(particle, GetAODEMCAL(),"EMCAL");
    //Trigger particle in EMCAL, correlation with PHOS
    else if(particle->GetDetector()=="EMCAL" && GetReader()->IsPHOSSwitchedOn() && GetAODPHOS()->GetEntriesFast() > 0)
      MakeNeutralCorrelationFillAOD(particle, GetAODPHOS(),"PHOS");
    //Trigger particle in CTS, correlation with PHOS, EMCAL and CTS
    else if(particle->GetDetector()=="CTS" ){
      if(GetReader()->IsPHOSSwitchedOn() && GetAODPHOS()->GetEntriesFast() > 0) 
	MakeNeutralCorrelationFillAOD(particle, GetAODPHOS(),"PHOS");
      if(GetReader()->IsEMCALSwitchedOn() && GetAODEMCAL()->GetEntriesFast() > 0) 
	MakeNeutralCorrelationFillAOD(particle, GetAODEMCAL(),"EMCAL");
    }
  
	
  }//Aod branch loop
  
  if(GetDebug() > 1) printf("AliAnaParticleHadronCorrelation::MakeAnalysisFillAOD() - End fill AODs \n");
  
}

//____________________________________________________________________________
void  AliAnaParticleHadronCorrelation::MakeAnalysisFillHistograms()  
{  
  //Particle-Hadron Correlation Analysis, fill histograms
  
  if(!GetInputAODBranch()){
    printf("AliAnaParticleHadronCorrelation::MakeAnalysisFillHistograms() - No input particles in AOD with name branch < %s >, ABORT \n",GetInputAODName().Data());
    abort();
  }
  
  if(GetDebug() > 1){
    printf("AliAnaParticleHadronCorrelation::MakeAnalysisFillHistograms() - Begin hadron correlation analysis, fill histograms \n");
    printf("AliAnaParticleHadronCorrelation::MakeAnalysisFillHistograms() - In particle branch aod entries %d\n", GetInputAODBranch()->GetEntriesFast());
  }
  
  //Loop on stored AOD particles
  Int_t naod = GetInputAODBranch()->GetEntriesFast();
  for(Int_t iaod = 0; iaod < naod ; iaod++){	  
    AliAODPWG4ParticleCorrelation* particle =  (AliAODPWG4ParticleCorrelation*) (GetInputAODBranch()->At(iaod));
    
    //check if the particle is isolated or if we want to take the isolation into account
    if(OnlyIsolated() && !particle->IsIsolated()) continue;
    
    //Make correlation with charged hadrons
    TRefArray * reftracks   = particle->GetRefArray(GetAODRefArrayName()+"Tracks");
    if(reftracks){
      if(GetDebug() > 1) printf("AliAnaParticleHadronCorrelation::MakeAnalysisFillHistograms() - Particle %d, In Track Refs  entries %d\n", iaod, reftracks->GetEntriesFast());
      if(reftracks->GetEntriesFast() > 0) MakeChargedCorrelation(particle, reftracks,kTRUE);
    }
    
    //Make correlation with neutral pions
    if(GetOutputAODBranch() && GetOutputAODBranch()->GetEntriesFast() > 0){
      if(GetDebug() > 1) printf("AliAnaParticleHadronCorrelation::MakeAnalysisFillHistograms() - Particle %d, In Cluster Refs entries %d\n",iaod, GetOutputAODBranch()->GetEntriesFast());      
      MakeNeutralCorrelationFillHistograms(particle);
    }
    
  }//Aod branch loop
  
  if(GetDebug() > 1) printf("AliAnaParticleHadronCorrelation::MakeAnalysisFillHistograms() - End fill histograms \n");
  
}

//____________________________________________________________________________
void  AliAnaParticleHadronCorrelation::MakeChargedCorrelation(AliAODPWG4ParticleCorrelation *aodParticle, TRefArray* pl, const Bool_t bFillHisto)
{  
  // Charged Hadron Correlation Analysis
  if(GetDebug() > 1)printf("AliAnaParticleHadronCorrelation::MakeChargedCorrelation() - Make trigger particle - charged hadron correlation \n");
  
  Double_t ptTrig  = aodParticle->Pt();
  Double_t phiTrig = aodParticle->Phi();
  Double_t pt   = -100.;
  Double_t rat  = -100.; 
  Double_t phi  = -100. ;
  Double_t eta  = -100. ;
  Double_t p[3];
  Bool_t   first=kTRUE;
  
  TRefArray * reftracks    =0x0;
  if(!bFillHisto) 
    reftracks    = new TRefArray;
  

  //Track loop, select tracks with good pt, phi and fill AODs or histograms
  for(Int_t ipr = 0;ipr < pl->GetEntries() ; ipr ++ ){
    AliAODTrack * track = (AliAODTrack *) (pl->At(ipr)) ;
    track->GetPxPyPz(p) ;
    TLorentzVector mom(p[0],p[1],p[2],0);
    pt   = mom.Pt();
    eta  = mom.Eta();
    phi  = mom.Phi() ;
    if(phi < 0) phi+=TMath::TwoPi();
    rat   = pt/ptTrig ;
    
    if(IsFidutialCutOn()){
      Bool_t in = GetFidutialCut()->IsInFidutialCut(mom,"CTS") ;
      if(! in ) continue ;
    }    

    //Select only hadrons in pt range
    if(pt < GetMinPt() || pt > GetMaxPt()) continue ;
    
    //Selection within angular range
    Float_t deltaphi = TMath::Abs(phiTrig-phi);
    if( (deltaphi < fDeltaPhiMinCut) || ( deltaphi > fDeltaPhiMaxCut) ) continue ;    
    
    if(GetDebug() > 2)
      printf("AliAnaParticleHadronCorrelation::MakeChargedCorrelation() - Charged hadron: pt %f, phi %f, phi trigger %f. Cuts:  delta phi  %2.2f < %2.2f < %2.2f, pT min %2.2f \n",
	     pt,phi, phiTrig,fDeltaPhiMinCut, deltaphi, fDeltaPhiMaxCut, GetMinPt());
    
    if(bFillHisto){
      // Fill Histograms
      fhEtaCharged->Fill(ptTrig,eta);
      fhPhiCharged->Fill(ptTrig,phi);
      fhDeltaEtaCharged->Fill(ptTrig,aodParticle->Eta()-eta);
      fhDeltaPhiCharged->Fill(ptTrig, deltaphi);
      fhDeltaPhiChargedPt->Fill(pt,deltaphi);
      if(GetDebug() > 2 ) printf("AliAnaParticleHadronCorrelation::MakeChargedCorrelation() - Selected charge for momentum imbalance: pt %2.2f, phi %2.2f, eta %2.2f \n",pt,phi,eta);
      fhPtImbalanceCharged->Fill(ptTrig,rat);
    }
    else{
      //Fill AODs
      
      if(first) {
	new (reftracks) TRefArray(TProcessID::GetProcessWithUID(track)); 
	first = kFALSE;
      }
      
      reftracks->Add(track);
      
    }//aod particle loop
  }// track loop
  
  //Fill AOD with reference tracks, if not filling histograms
  if(!bFillHisto && reftracks->GetEntriesFast() > 0) {
    reftracks->SetName(GetAODRefArrayName()+"Tracks");
    aodParticle->AddRefArray(reftracks);
  }
  
}  

//____________________________________________________________________________
void  AliAnaParticleHadronCorrelation::MakeNeutralCorrelationFillAOD(AliAODPWG4ParticleCorrelation * aodParticle,TRefArray* pl, TString detector)  
{  
  // Neutral Pion Correlation Analysis, find pi0, put them in new output aod, if correlation cuts passed
  if(GetDebug() > 1) printf("AliAnaParticleHadronCorrelation::MakeNeutralCorrelationFillAOD() - Make trigger particle - neutral hadron correlation \n");
  
  if(!NewOutputAOD()){
    printf("AliAnaParticleHadronCorrelation::MakeNeutralCorrelationFillAOD() - Output aod not created, set AOD class name and branch name in the configuration file, ABORT! \n");
    abort();
  }
  
  Double_t phiTrig = aodParticle->Phi();
  Int_t	tag = -1;
  TLorentzVector gammai;
  TLorentzVector gammaj;
  
  Double_t vertex[] = {0,0,0};
  if(!GetReader()->GetDataType()== AliCaloTrackReader::kMC) GetReader()->GetVertex(vertex);
  
  //Cluster loop, select pairs with good pt, phi and fill AODs or histograms
  //Int_t iEvent= GetReader()->GetEventNumber() ;
  Int_t nclus = pl->GetEntriesFast();
  for(Int_t iclus = 0;iclus < nclus ; iclus ++ ){
    AliAODCaloCluster * calo = (AliAODCaloCluster *) (pl->At(iclus)) ;
    
    //Cluster selection, not charged, with photon or pi0 id and in fidutial cut
    Int_t pdg=0;
    if(!SelectCluster(calo, vertex, gammai, pdg)) continue ;

	if(GetDebug() > 2)
      printf("AliAnaParticleHadronCorrelation::MakeNeutralCorrelationFillAOD() - Neutral cluster in %s: pt %f, phi %f, phi trigger %f. Cuts:  delta phi min %2.2f,  max %2.2f, pT min %2.2f \n",
	     detector.Data(), gammai.Pt(),gammai.Phi(),phiTrig,fDeltaPhiMinCut, fDeltaPhiMaxCut, GetMinPt());
    
    //2 gamma overlapped, found with PID
    if(pdg == AliCaloPID::kPi0){ 
      
      //Select only hadrons in pt range
      if(gammai.Pt() < GetMinPt() || gammai.Pt() > GetMaxPt()) continue ;
      
      //Selection within angular range
      Float_t phi = gammai.Phi();
      if(phi < 0) phi+=TMath::TwoPi();
      Float_t deltaphi = TMath::Abs(phiTrig-phi);
      if( (deltaphi < fDeltaPhiMinCut) || ( deltaphi > fDeltaPhiMaxCut) ) continue ;
      
      AliAODPWG4Particle pi0 = AliAODPWG4Particle(gammai);
      //pi0.SetLabel(calo->GetLabel(0));
      pi0.SetPdg(AliCaloPID::kPi0);
      pi0.SetDetector(detector);
      
      if(IsDataMC()){
	pi0.SetTag(GetMCAnalysisUtils()->CheckOrigin(calo->GetLabel(0),GetMCStack()));
	if(GetDebug() > 0) printf("AliAnaParticleHadronCorrelation::MakeNeutralCorrelationFillAOD() - Origin of candidate %d\n",pi0.GetTag());
      }//Work with stack also 
      //Set the indeces of the original caloclusters  
      pi0.SetCaloLabel(calo->GetID(),-1);
      AddAODParticle(pi0);
      
      if(GetDebug() > 2) 
	printf("AliAnaParticleHadronCorrelation::MakeNeutralCorrelationFillAOD() - Correlated with selected pi0 (pid): pt %f, phi %f\n",pi0.Pt(),pi0.Phi());
      
    }// pdg = 111
    
    //Make invariant mass analysis
    else if(pdg == AliCaloPID::kPhoton){	
      //Search the photon companion in case it comes from  a Pi0 decay
      //Apply several cuts to select the good pair;
      for(Int_t jclus = iclus+1; jclus < pl->GetEntries() ; jclus ++ ){
	AliAODCaloCluster * calo2 = (AliAODCaloCluster *) (pl->At(jclus)) ;
	
	//Cluster selection, not charged with photon or pi0 id and in fidutial cut
	Int_t pdgj=0;
	if(!SelectCluster(calo2,vertex, gammaj, pdgj)) continue ;
	
	if(pdgj == AliCaloPID::kPhoton ){
	  
	  if((gammai+gammaj).Pt() < GetMinPt() || (gammai+gammaj).Pt() > GetMaxPt()) continue ;
	  
	  //Selection within angular range
	  Float_t phi = (gammai+gammaj).Phi();
	  if(phi < 0) phi+=TMath::TwoPi();
	  Float_t deltaphi = TMath::Abs(phiTrig-phi);
	  if( (deltaphi < fDeltaPhiMinCut) || ( deltaphi > fDeltaPhiMaxCut) ) continue ;
	  
	  //Select good pair (aperture and invariant mass)
	  if(GetNeutralMesonSelection()->SelectPair(gammai, gammaj)){
	    
	    if(GetDebug() > 2 ) printf("AliAnaParticleHadronCorrelation::MakeNeutralCorrelationFillAOD() - Neutral Hadron Correlation: AOD Selected gamma pair: pt %2.2f, phi %2.2f, eta %2.2f, M %2.3f\n",
				       (gammai+gammaj).Pt(),(gammai+gammaj).Phi(),(gammai+gammaj).Eta(), (gammai+gammaj).M());
	    
	    TLorentzVector pi0mom = gammai+gammaj;
	    AliAODPWG4Particle pi0 = AliAODPWG4Particle(pi0mom);
	    //pi0.SetLabel(calo->GetLabel(0));
	    pi0.SetPdg(AliCaloPID::kPi0);
	    pi0.SetDetector(detector);	
	    if(IsDataMC()){
	      //Check origin of the candidates
	      Int_t tag1 = GetMCAnalysisUtils()->CheckOrigin(calo->GetLabel(0),  GetMCStack());
	      Int_t tag2 = GetMCAnalysisUtils()->CheckOrigin(calo2->GetLabel(0), GetMCStack());
	      
	      if(GetDebug() > 0) printf("AliAnaParticleHadronCorrelation::MakeNeutralCorrelationFillAOD() - Origin of: photon1 %d; photon2 %d \n",tag1, tag2);
	      if(tag1 == AliMCAnalysisUtils::kMCPi0Decay && tag2 == AliMCAnalysisUtils::kMCPi0Decay){
		
		//Check if pi0 mother is the same
		Int_t label1 = calo->GetLabel(0);
		TParticle * mother1 = GetMCStack()->Particle(label1);//photon in kine tree
		label1 = mother1->GetFirstMother();
		//mother1 = GetMCStack()->Particle(label1);//pi0
		
		Int_t label2 = calo2->GetLabel(0);
		TParticle * mother2 = GetMCStack()->Particle(label2);//photon in kine tree
		label2 = mother2->GetFirstMother();
		//mother2 = GetMCStack()->Particle(label2);//pi0
		
		//printf("mother1 %d, mother2 %d\n",label1,label2);
		if(label1 == label2)
		  tag = AliMCAnalysisUtils::kMCPi0;
	      }
	    }//Work with stack also   
	    pi0.SetTag(tag);
	    //Set the indeces of the original caloclusters  
	    pi0.SetCaloLabel(calo->GetID(), calo2->GetID());
	    AddAODParticle(pi0);
	    
	    
	  }//Pair selected
	}//if pair of gammas
      }//2nd loop
    }// if pdg = 22
  }//1st loop
  
  if(GetDebug() > 1) 
    printf("AliAnaParticleHadronCorrelation::MakeNeutralCorrelationFillAOD() - End, %d pi0's found \n",GetOutputAODBranch()->GetEntriesFast());
}

//____________________________________________________________________________
void  AliAnaParticleHadronCorrelation::MakeNeutralCorrelationFillHistograms(AliAODPWG4ParticleCorrelation * aodParticle)  
{  
  // Neutral Pion Correlation Analysis
  if(GetDebug() > 1) printf("AliAnaParticleHadronCorrelation::MakeNeutralCorrelationFillHistogramS() - Make trigger particle - pi0 correlation, %d pi0's \n",GetOutputAODBranch()->GetEntriesFast());
  
  Double_t pt  = -100.;
  Double_t rat = -100.; 
  Double_t phi = -100.;
  Double_t eta = -100.;
  Double_t ptTrig  = aodParticle->Pt();
  Double_t phiTrig = aodParticle->Phi();
  Double_t etaTrig = aodParticle->Eta();
  
  if(!GetOutputAODBranch()){
    printf("AliAnaParticleHadronCorrelation::MakeNeutralCorrelationFillHistograms() - No output pi0 in AOD branch with name < %s >,STOP \n",GetOutputAODName().Data());
    abort();
  }
  
  //Loop on stored AOD pi0
  Int_t naod = GetOutputAODBranch()->GetEntriesFast();
  if(GetDebug() > 0) printf("AliAnaParticleHadronCorrelation::MakeNeutralCorrelationFillHistograms() -  aod branch entries %d\n", naod);
  for(Int_t iaod = 0; iaod < naod ; iaod++){
    AliAODPWG4Particle* pi0 =  (AliAODPWG4Particle*) (GetOutputAODBranch()->At(iaod));
    Int_t pdg = pi0->GetPdg();
    
    if(pdg != AliCaloPID::kPi0) continue;          	
    pt  = pi0->Pt();
    
    if(pt < GetMinPt() || pt > GetMaxPt()) continue ;
    
    //Selection within angular range
    phi = pi0->Phi();
    Float_t deltaphi = TMath::Abs(phiTrig-phi);
    if( (deltaphi < fDeltaPhiMinCut) || ( deltaphi > fDeltaPhiMaxCut) ) continue ;
    
    rat = pt/ptTrig ;
    phi = pi0->Phi() ;
    eta = pi0->Eta() ;
    
    fhEtaNeutral->Fill(ptTrig,eta);
    fhPhiNeutral->Fill(ptTrig,phi);
    fhDeltaEtaNeutral->Fill(ptTrig,etaTrig-eta);
    fhDeltaPhiNeutral->Fill(ptTrig,deltaphi);
    fhDeltaPhiNeutralPt->Fill(pt,deltaphi);
    fhPtImbalanceNeutral->Fill(ptTrig,rat);
    
    if(GetDebug() > 2 ) printf("AliAnaParticleHadronCorrelation::MakeNeutralCorrelation() - Selected neutral for momentum imbalance: pt %2.2f, phi %2.2f, eta %2.2f \n",pt,phi,eta);
    
  }//loop
}


//____________________________________________________________________________
Bool_t  AliAnaParticleHadronCorrelation::SelectCluster(AliAODCaloCluster * calo, Double_t *vertex, TLorentzVector & mom, Int_t & pdg) const {
  //Select cluster depending on its pid and acceptance selections
  
  //Skip matched clusters with tracks
  if(calo->GetNTracksMatched() > 0) return kFALSE;
  
  TString detector = "";
  if(calo->IsPHOSCluster()) detector= "PHOS";
  else if(calo->IsEMCALCluster()) detector= "EMCAL";
		
  //Check PID
  calo->GetMomentum(mom,vertex);//Assume that come from vertex in straight line
  pdg = AliCaloPID::kPhoton;   
  if(IsCaloPIDOn()){
    //Get most probable PID, 2 options check PID weights (in MC this option is mandatory)
    //or redo PID, recommended option for EMCal.
    
    if(!IsCaloPIDRecalculationOn() || GetReader()->GetDataType() == AliCaloTrackReader::kMC )
      pdg = GetCaloPID()->GetPdg(detector,calo->PID(),mom.E());//PID with weights
    else
      pdg = GetCaloPID()->GetPdg(detector,mom,calo);//PID recalculated
    
    if(GetDebug() > 5) printf("AliAnaParticleHadronCorrelation::SelectCluster() - PDG of identified particle %d\n",pdg);
    
    //If it does not pass pid, skip
    if(pdg != AliCaloPID::kPhoton && pdg != AliCaloPID::kPi0) {
      return kFALSE ;
    }
  }//PID on
  
  //Check acceptance selection
  if(IsFidutialCutOn()){
    Bool_t in = GetFidutialCut()->IsInFidutialCut(mom,detector) ;
    if(! in ) return kFALSE ;
  }
  
  if(GetDebug() > 5) printf("AliAnaParticleHadronCorrelation::SelectCluster() - Correlation photon selection cuts passed: pT %3.2f, pdg %d\n",mom.Pt(), pdg);
  
  return kTRUE;
  
}
