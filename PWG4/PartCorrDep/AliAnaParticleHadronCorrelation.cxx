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

//---- ANALYSIS system ----
#include "AliLog.h"
#include "AliNeutralMesonSelection.h" 
#include "AliAnaParticleHadronCorrelation.h" 
#include "AliCaloTrackReader.h"
#include "AliCaloPID.h"
#include "AliAODPWG4ParticleCorrelation.h"

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
		("DeltaPhiChargedPt","#phi_{trigger} - #phi_{#p^{#pm}i} vs p_{T h^{#pm}}",
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
  SetInputAODName("photons");
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
  
  Info("*** Print *** ", "%s %s", GetName(), GetTitle() ) ;
  printf("pT Hadron       >    %2.2f\n", GetMinPt()) ; 
  printf("pT Hadron       <    %2.2f\n", GetMaxPt()) ; 
  printf("Phi trigger particle-Hadron      <     %3.2f\n", fDeltaPhiMaxCut) ; 
  printf("Phi trigger particle-Hadron      >     %3.2f\n", fDeltaPhiMinCut) ;
  printf("Isolated Trigger?  %d\n", fSelectIsolated) ;
} 

//____________________________________________________________________________
void  AliAnaParticleHadronCorrelation::MakeAnalysisFillAOD()  
{  
	//Particle-Hadron Correlation Analysis, fill AODs
	
	if(!GetInputAODBranch())
		AliFatal(Form("ParticleHadronCorrelation::FillAOD: No input particles in AOD with name branch < %s > \n",GetInputAODName().Data()));
	
	if(GetDebug() > 1){
		printf("Begin hadron correlation analysis, fill AODs \n");
		printf("In particle branch aod entries %d\n", GetInputAODBranch()->GetEntriesFast());
		printf("In CTS aod entries %d\n", GetAODCTS()->GetEntriesFast());
		printf("In EMCAL aod entries %d\n", GetAODEMCAL()->GetEntriesFast());
		printf("In PHOS aod entries %d\n", GetAODPHOS()->GetEntriesFast());
	}
	
	//Loop on stored AOD particles, trigger
	Int_t naod = GetInputAODBranch()->GetEntriesFast();
	for(Int_t iaod = 0; iaod < naod ; iaod++){
		AliAODPWG4ParticleCorrelation* particle =  (AliAODPWG4ParticleCorrelation*) (GetInputAODBranch()->At(iaod));
		//Make correlation with charged hadrons
		if(GetReader()->IsCTSSwitchedOn() )
			MakeChargedCorrelation(particle, (TSeqCollection*)GetAODCTS(),kFALSE);
		
		//Make correlation with neutral pions
		//Trigger particle in PHOS, correlation with EMCAL
		if(particle->GetDetector()=="PHOS" && GetReader()->IsEMCALSwitchedOn() && GetAODEMCAL()->GetEntriesFast() > 0)
			MakeNeutralCorrelation(particle,(TSeqCollection*)GetAODEMCAL(),kFALSE);
		//Trigger particle in EMCAL, correlation with PHOS
		else if(particle->GetDetector()=="EMCAL" && GetReader()->IsPHOSSwitchedOn() && GetAODPHOS()->GetEntriesFast() > 0)
			MakeNeutralCorrelation(particle,(TSeqCollection*)GetAODPHOS(),kFALSE);
		//Trigger particle in CTS, correlation with PHOS, EMCAL and CTS
		else if(particle->GetDetector()=="CTS" ){
			if(GetReader()->IsPHOSSwitchedOn() && GetAODPHOS()->GetEntriesFast() > 0) 
				MakeNeutralCorrelation(particle,(TSeqCollection*)GetAODPHOS(),kFALSE);
			if(GetReader()->IsEMCALSwitchedOn() && GetAODEMCAL()->GetEntriesFast() > 0) 
				MakeNeutralCorrelation(particle,(TSeqCollection*)GetAODEMCAL(),kFALSE);
		}
		
	}//Aod branch loop
	
	if(GetDebug() > 1) printf("End hadron correlation analysis, fill AODs \n");
	
}

//____________________________________________________________________________
void  AliAnaParticleHadronCorrelation::MakeAnalysisFillHistograms()  
{  
	//Particle-Hadron Correlation Analysis, fill histograms
	
	if(!GetInputAODBranch())
		AliFatal(Form("ParticleHadronCorrelation::FillHistos: No input particles in AOD with name branch < %s > \n",GetInputAODName().Data()));
	
	if(GetDebug() > 1){
		printf("Begin hadron correlation analysis, fill histograms \n");
		printf("In particle branch aod entries %d\n", GetInputAODBranch()->GetEntriesFast());
	}
	
	//Loop on stored AOD particles
	Int_t naod = GetInputAODBranch()->GetEntriesFast();
	for(Int_t iaod = 0; iaod < naod ; iaod++){
		AliAODPWG4ParticleCorrelation* particle =  (AliAODPWG4ParticleCorrelation*) (GetInputAODBranch()->At(iaod));
		
		if(GetDebug() > 1){
			printf("Particle %d, In Track Refs  entries %d\n", iaod, (particle->GetRefTracks())->GetEntriesFast());
			printf("Particle %d, In Cluster Refs entries %d\n",iaod, (particle->GetRefClusters())->GetEntriesFast());
		}
		
		if(OnlyIsolated() && !particle->IsIsolated()) continue;
		
		//Make correlation with charged hadrons
		if((particle->GetRefTracks())->GetEntriesFast() > 0)
			MakeChargedCorrelation(particle, (TSeqCollection*) (particle->GetRefTracks()),kTRUE);
		
		//Make correlation with neutral pions
		if((particle->GetRefClusters())->GetEntriesFast() > 0)
			MakeNeutralCorrelation(particle,  (TSeqCollection*) (particle->GetRefClusters()), kTRUE);
		
	}//Aod branch loop
	
	if(GetDebug() > 1) printf("End hadron correlation analysis, fill histograms \n");
	
}

//____________________________________________________________________________
void  AliAnaParticleHadronCorrelation::MakeChargedCorrelation(AliAODPWG4ParticleCorrelation *aodParticle, TSeqCollection* pl, const Bool_t bFillHisto)
{  
  // Charged Hadron Correlation Analysis
  if(GetDebug() > 1)printf("Make trigger particle - charged hadron correlation \n");
  
  Double_t ptTrig  = aodParticle->Pt();
  Double_t phiTrig = aodParticle->Phi();
  Double_t pt   = -100.;
  Double_t rat  = -100.; 
  Double_t phi  = -100. ;
  Double_t eta  = -100. ;
  Double_t p[3];
  
  //Track loop, select tracks with good pt, phi and fill AODs or histograms
  for(Int_t ipr = 0;ipr < pl->GetEntries() ; ipr ++ ){
    AliAODTrack * track = (AliAODTrack *) (pl->At(ipr)) ;
    track->GetPxPyPz(p) ;
    TLorentzVector mom(p[0],p[1],p[2],0);
    pt    = mom.Pt();
    eta  = mom.Eta();
    phi  = mom.Phi() ;
    if(phi<0) phi+=TMath::TwoPi();
    rat   = pt/ptTrig ;

    if(IsFidutialCutOn()){
      Bool_t in = GetFidutialCut()->IsInFidutialCut(mom,"CTS") ;
      if(! in ) continue ;
    }    

    //Select only hadrons in pt range
    if(pt < GetMinPt() || pt > GetMaxPt()) continue ;
    
    if(GetDebug() > 2)
      printf("charged hadron: pt %f, phi %f, phi trigger %f. Cuts:  delta phi min %2.2f,  max%2.2f, pT min %2.2f \n",
	     pt,phi,phiTrig,fDeltaPhiMinCut, fDeltaPhiMaxCut, GetMinPt());

    if(bFillHisto){
      // Fill Histograms
      fhEtaCharged->Fill(ptTrig,eta);
      fhPhiCharged->Fill(ptTrig,phi);
      fhDeltaEtaCharged->Fill(ptTrig,aodParticle->Eta()-eta);
      fhDeltaPhiCharged->Fill(ptTrig,phiTrig-phi);
      fhDeltaPhiChargedPt->Fill(pt,phiTrig-phi);
      //Selection within angular range
      if(((phiTrig-phi)> fDeltaPhiMinCut) && ((phiTrig-phi)<fDeltaPhiMaxCut) ){
	if(GetDebug() > 2 ) printf("Selected charge for momentum imbalance: pt %2.2f, phi %2.2f, eta %2.2f ",pt,phi,eta);
	fhPtImbalanceCharged->Fill(ptTrig,rat);
      } 
    }
    else{
      //Fill AODs
      aodParticle->AddTrack(track);
      
    }//aod particle loop
  }// track loop

}  
//____________________________________________________________________________
void  AliAnaParticleHadronCorrelation::MakeNeutralCorrelation(AliAODPWG4ParticleCorrelation * aodParticle,TSeqCollection* pl, const Bool_t bFillHisto)  
{  
  // Neutral Pion Correlation Analysis
  if(GetDebug() > 1) printf("Make trigger particle - neutral hadron correlation \n");
  
  Double_t pt = -100.;
  Double_t rat = -100.; 
  Double_t phi = -100. ;
  Double_t eta = -100. ;
  Double_t ptTrig  = aodParticle->Pt();
  Double_t phiTrig = aodParticle->Phi();
  Double_t etaTrig = aodParticle->Eta();
  
  TLorentzVector gammai;
  TLorentzVector gammaj;
  
  Double_t vertex[] = {0,0,0};

  if(!GetReader()->GetDataType()== AliCaloTrackReader::kMC) GetReader()->GetVertex(vertex);
  
  //Cluster loop, select pairs with good pt, phi and fill AODs or histograms
  for(Int_t iclus = 0;iclus < pl->GetEntries() ; iclus ++ ){
    AliAODCaloCluster * calo = (AliAODCaloCluster *) (pl->At(iclus)) ;
    if(!bFillHisto){

      //Cluster selection, not charged, with photon or pi0 id and in fidutial cut
      Int_t pdg=0;
      if(!SelectCluster(calo, vertex, gammai, pdg)) continue ;
      
    if(GetDebug() > 2)
      printf("neutral cluster: pt %f, phi %f, phi trigger %f. Cuts:  delta phi min %2.2f,  max%2.2f, pT min %2.2f \n",
	     gammai.Pt(),gammai.Phi(),phiTrig,fDeltaPhiMinCut, fDeltaPhiMaxCut, GetMinPt());
    
      //2 gamma overlapped, found with PID
      if(pdg == AliCaloPID::kPi0){ 

	//Select only hadrons in pt range
	if(gammai.Pt() < GetMinPt() || gammai.Pt() > GetMaxPt()) continue ;

	aodParticle->AddCluster(calo);
	if(GetDebug() > 2) printf("Correlated with selected pi0 (pid): pt %f, phi %f",gammai.Pt(),gammai.Phi());

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
	    
	    //Select good pair (aperture and invariant mass)
	    if(GetNeutralMesonSelection()->SelectPair(gammai, gammaj)){

	      if(GetDebug() > 2 ) printf("Neutral Hadron Correlation: AOD Selected gamma pair: pt %2.2f, phi %2.2f, eta %2.2f, M %2.3f",
					 (gammai+gammaj).Pt(),(gammai+gammaj).Phi(),(gammai+gammaj).Eta(), (gammai+gammaj).M());
	      Int_t labels[]={calo->GetLabel(0),calo2->GetLabel(0)};
	      Float_t pid[]={0,0,0,0,0,0,1,0,0,0,0,0,0};//Pi0 weight 1
	      Float_t pos[]={(gammai+gammaj).X(), (gammai+gammaj).Y(), (gammai+gammaj).Z()};

	      AliAODCaloCluster *caloCluster =  new AliAODCaloCluster(0,2,labels,(gammai+gammaj).E(), pos, pid,calo->GetType(),0);
	      aodParticle->AddCluster(caloCluster);
	    }//Pair selected
	  }//if pair of gammas
	}//2nd loop
      }// if pdg = 22
    }// Fill AODs
    else{ //Fill histograms
      
      calo->GetMomentum(gammai,vertex);//Assume that come from vertex in straight line
      pt  = gammai.Pt();
      
      if(pt < GetMinPt() || pt > GetMaxPt()) continue ;
      
      rat = pt/ptTrig ;
      phi = gammai.Phi() ;
      eta = gammai.Eta() ;

      if(GetDebug() > 2 ) printf("Neutral Hadron Correlation: Histograms selected gamma pair: pt %2.2f, phi %2.2f, eta %2.2f",pt,phi,eta);
      
      fhEtaNeutral->Fill(ptTrig,eta);
      fhPhiNeutral->Fill(ptTrig,phi);
      fhDeltaEtaNeutral->Fill(ptTrig,etaTrig-eta);
      fhDeltaPhiNeutral->Fill(ptTrig,phiTrig-phi);
      fhDeltaPhiNeutralPt->Fill(pt,phiTrig-phi);
      //Selection within angular range
      if(((phiTrig-phi)> fDeltaPhiMinCut) && ((phiTrig-phi)<fDeltaPhiMaxCut) ){
	if(GetDebug() > 2 ) printf("Selected neutral for momentum imbalance: pt %2.2f, phi %2.2f, eta %2.2f ",pt,phi,eta);
	fhPtImbalanceNeutral->Fill(ptTrig,rat);
      }    
    }//Fill histograms
  }//1st loop
  
}

//____________________________________________________________________________
Bool_t  AliAnaParticleHadronCorrelation::SelectCluster(AliAODCaloCluster * calo, Double_t *vertex, TLorentzVector & mom, Int_t & pdg) const {
   //Select cluster depending on its pid and acceptance selections
   
   //Skip matched clusters with tracks
  if(calo->GetNTracksMatched() > 0) {
  return kFALSE;} 
   
   //Check PID
   calo->GetMomentum(mom,vertex);//Assume that come from vertex in straight line
   pdg = AliCaloPID::kPhoton;   
   if(IsCaloPIDOn()){
     //Get most probable PID, 2 options check PID weights (in MC this option is mandatory)
     //or redo PID, recommended option for EMCal.
     TString detector = "";
     if(calo->IsPHOSCluster()) detector= "PHOS";
     else if(calo->IsEMCALCluster()) detector= "EMCAL";

     if(!IsCaloPIDRecalculationOn() || GetReader()->GetDataType() == AliCaloTrackReader::kMC )
       pdg = GetCaloPID()->GetPdg(detector,calo->PID(),mom.E());//PID with weights
     else
       pdg = GetCaloPID()->GetPdg(detector,mom,calo);//PID recalculated
     
     if(GetDebug() > 1) printf("PDG of identified particle %d\n",pdg);
     
     //If it does not pass pid, skip
     if(pdg != AliCaloPID::kPhoton || pdg != AliCaloPID::kPi0) {
       return kFALSE ;
     }
   }
   
   //Check acceptance selection
   if(IsFidutialCutOn()){
     Bool_t in = kFALSE;
     if(calo->IsPHOSCluster())
       in = GetFidutialCut()->IsInFidutialCut(mom,"PHOS") ;
     else if(calo->IsEMCALCluster())
       in = GetFidutialCut()->IsInFidutialCut(mom,"EMCAL") ;
     if(! in ) { return kFALSE ;}
   }
   
   if(GetDebug() > 1) printf("Correlation photon selection cuts passed: pT %3.2f, pdg %d\n",mom.Pt(), pdg);
   
   return kTRUE;
   
 }
