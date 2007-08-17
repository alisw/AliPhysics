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
/* $Id$ */

/* History of cvs commits:
 *
 * $Log$
 * Revision 1.4.4.4  2007/07/26 10:32:09  schutz
 * new analysis classes in the the new analysis framework
 *
 *
 */

//_________________________________________________________________________
// Class for the prompt gamma analysis, isolation cut
//
//  Class created from old AliPHOSGammaJet 
//  (see AliRoot versions previous Release 4-09)
//
//*-- Author: Gustavo Conesa (LNF-INFN) 
//////////////////////////////////////////////////////////////////////////////
  
  
// --- ROOT system --- 
#include <TParticle.h>
#include <TH2.h>
#include <TList.h>
#include "AliAnaGammaDirect.h" 
#include "Riostream.h"
#include "AliLog.h"
  
ClassImp(AliAnaGammaDirect)
  
//____________________________________________________________________________
  AliAnaGammaDirect::AliAnaGammaDirect() : 
    TObject(),
    fMinGammaPt(0.),
    fConeSize(0.),fPtThreshold(0.),fPtSumThreshold(0), 
    fICMethod(0),fhNGamma(0),fhPhiGamma(0),fhEtaGamma(0),  
    //kSeveralIC
    fNCones(0),fNPtThres(0), fConeSizes(),  fPtThresholds(), 
    fhPtThresIsolated(), fhPtSumIsolated()
{
  //default ctor
  
  //Initialize parameters
  InitParameters();

}

//____________________________________________________________________________
AliAnaGammaDirect::AliAnaGammaDirect(const AliAnaGammaDirect & g) : 
  TObject(g),
  fMinGammaPt(g.fMinGammaPt), 
  fConeSize(g.fConeSize),
  fPtThreshold(g.fPtThreshold),
  fPtSumThreshold(g.fPtSumThreshold), 
  fICMethod(g.fICMethod),
  fhNGamma(g.fhNGamma),fhPhiGamma(g.fhPhiGamma),fhEtaGamma(g.fhEtaGamma),  
  //kSeveralIC
  fNCones(g.fNCones),fNPtThres(g.fNPtThres), fConeSizes(),fPtThresholds(), 
  fhPtThresIsolated(), fhPtSumIsolated()
{
  // cpy ctor
  
  //kSeveralIC
  for(Int_t i = 0; i < fNCones ; i++){
    fConeSizes[i] =  g.fConeSizes[i];
    fhPtSumIsolated[i] = g.fhPtSumIsolated[i]; 
    for(Int_t j = 0; j < fNPtThres ; j++)
      fhPtThresIsolated[i][j] = g.fhPtThresIsolated[i][j]; 
  }
  
  for(Int_t i = 0; i < fNPtThres ; i++)
    fPtThresholds[i]=   g.fPtThresholds[i];
}

//_________________________________________________________________________
AliAnaGammaDirect & AliAnaGammaDirect::operator = (const AliAnaGammaDirect & source)
{
  // assignment operator
  
  if(&source == this) return *this;

  fMinGammaPt = source.fMinGammaPt ;   
  fConeSize = source.fConeSize ;
  fPtThreshold = source.fPtThreshold ;
  fPtSumThreshold = source.fPtSumThreshold ; 
  fICMethod = source.fICMethod ;
  fhNGamma = source.fhNGamma ; 
  fhPhiGamma = source.fhPhiGamma ;
  fhEtaGamma = source.fhEtaGamma ;
  
  //kSeveralIC
  fNCones = source.fNCones ;
  fNPtThres = source.fNPtThres ; 
   
  for(Int_t i = 0; i < fNCones ; i++){
    fConeSizes[i] =  source.fConeSizes[i];
    fhPtSumIsolated[i] = source.fhPtSumIsolated[i] ;
    for(Int_t j = 0; j < fNPtThres ; j++)
      fhPtThresIsolated[i][j] = source.fhPtThresIsolated[i][j] ;
  }
  
  for(Int_t i = 0; i < fNPtThres ; i++)
    fPtThresholds[i]=   source.fPtThresholds[i];
  
  return *this;
  
}

//____________________________________________________________________________
AliAnaGammaDirect::~AliAnaGammaDirect() 
{
  // Remove all pointers
  
  delete fhNGamma    ;  
  delete fhPhiGamma  ; 
  delete fhEtaGamma   ;  
  
  //kSeveralIC
  delete [] fhPtThresIsolated ;
  delete [] fhPtSumIsolated ;
  
}

//________________________________________________________________________
TList *  AliAnaGammaDirect::GetCreateOutputObjects()
{  

  // Create histograms to be saved in output file and 
  // store them in outputContainer
  TList * outputContainer = new TList() ; 
  outputContainer->SetName("DirectGammaHistos") ; 
  
  //Histograms of highest gamma identified in Event
  fhNGamma  = new TH1F("NGamma","Number of #gamma over calorimeter",240,0,120); 
  fhNGamma->SetYTitle("N");
  fhNGamma->SetXTitle("p_{T #gamma}(GeV/c)");
  outputContainer->Add(fhNGamma) ; 
  
  fhPhiGamma  = new TH2F
    ("PhiGamma","#phi_{#gamma}",200,0,120,200,0,7); 
  fhPhiGamma->SetYTitle("#phi");
  fhPhiGamma->SetXTitle("p_{T #gamma} (GeV/c)");
  outputContainer->Add(fhPhiGamma) ; 
  
  fhEtaGamma  = new TH2F
    ("EtaGamma","#phi_{#gamma}",200,0,120,200,-0.8,0.8); 
  fhEtaGamma->SetYTitle("#eta");
  fhEtaGamma->SetXTitle("p_{T #gamma} (GeV/c)");
  outputContainer->Add(fhEtaGamma) ;

  if(fICMethod == kSeveralIC){
    char name[128];
    char title[128];
    for(Int_t icone = 0; icone<fNCones; icone++){
      sprintf(name,"PtSumIsolated_Cone_%d",icone);
      sprintf(title,"Candidate cone sum p_{T} for cone size %d vs candidate p_{T}",icone);
      fhPtSumIsolated[icone]  = new TH2F(name, title,240,0,120,120,0,10);
      fhPtSumIsolated[icone]->SetYTitle("#Sigma p_{T} (GeV/c)");
      fhPtSumIsolated[icone]->SetXTitle("p_{T} (GeV/c)");
      outputContainer->Add(fhPtSumIsolated[icone]) ; 
      
      for(Int_t ipt = 0; ipt<fNPtThres;ipt++){ 
	sprintf(name,"PtThresIsol_Cone_%d_Pt%d",icone,ipt);
	sprintf(title,"Isolated candidate p_{T} distribution for cone size %d and p_{T}^{th} %d",icone,ipt);
	fhPtThresIsolated[icone][ipt]  = new TH1F(name, title,240,0,120);
	fhPtThresIsolated[icone][ipt]->SetXTitle("p_{T} (GeV/c)");
	outputContainer->Add(fhPtThresIsolated[icone][ipt]) ; 
      }//icone loop
    }//ipt loop
  }

  return outputContainer ;

}

//____________________________________________________________________________
void AliAnaGammaDirect::GetPromptGamma(TClonesArray * pl, TClonesArray * plCTS, TParticle *pGamma, Bool_t &found) const 
{
  //Search for the prompt photon in Calorimeter with pt > fMinGammaPt

  Double_t pt = 0;
  Int_t index = -1; 
  for(Int_t ipr = 0;ipr < pl->GetEntries() ; ipr ++ ){
    TParticle * particle = dynamic_cast<TParticle *>(pl->At(ipr)) ;

    if((particle->Pt() > fMinGammaPt) && (particle->Pt() > pt)){
      index = ipr ;
      pt = particle->Pt();
      pGamma->SetMomentum(particle->Px(),particle->Py(),particle->Pz(),particle->Energy());
      found  = kTRUE;
    }
  }

  //Do Isolation?
  if( ( fICMethod == kPtIC  ||  fICMethod == kSumPtIC )  && found)
    {
      Float_t coneptsum = 0 ;
      Bool_t  icPtThres   = kFALSE;
      Bool_t  icPtSum     = kFALSE;
      MakeIsolationCut(plCTS,pl, pGamma, index, 
		       icPtThres, icPtSum,coneptsum);
      if(fICMethod == kPtIC) //Pt thres method
	found = icPtThres ;
      if(fICMethod == kSumPtIC) //Pt cone sum method
	found = icPtSum ;
    }
  
  if(found){
    AliDebug(1,Form("Cluster with p_{T} larger than %f found in calorimeter ", fMinGammaPt)) ;
    AliDebug(1,Form("Gamma: pt %f, phi %f, eta %f", pGamma->Pt(),pGamma->Phi(),pGamma->Eta())) ;
    //Fill prompt gamma histograms 
    fhNGamma->Fill(pGamma->Pt());
    fhPhiGamma->Fill( pGamma->Pt(),pGamma->Phi());
    fhEtaGamma->Fill(pGamma->Pt(),pGamma->Eta());
  }
  else
    AliDebug(1,Form("NO Cluster with pT larger than %f found in calorimeter ", fMinGammaPt)) ;
}

  //____________________________________________________________________________
void AliAnaGammaDirect::InitParameters()
{
 
  //Initialize the parameters of the analysis.
  fMinGammaPt  = 5. ;

  //Fill particle lists when PID is ok
  fConeSize             = 0.2 ; 
  fPtThreshold         = 2.0; 
  fPtSumThreshold  = 1.; 

  fICMethod = kNoIC; // 0 don't isolate, 1 pt thresh method, 2 cone pt sum method

 //-----------kSeveralIC-----------------
  fNCones           = 4 ; 
  fNPtThres         = 4 ; 
  fConeSizes[0] = 0.1; fConeSizes[1] = 0.2; fConeSizes[2] = 0.3; fConeSizes[3] = 0.4;
  fPtThresholds[0]=1.; fPtThresholds[1]=2.; fPtThresholds[2]=3.; fPtThresholds[3]=4.;

}

//__________________________________________________________________
void  AliAnaGammaDirect::MakeIsolationCut(TClonesArray * plCTS, 
					   TClonesArray * plNe, 
					   TParticle * pCandidate, 
					   Int_t index, 
					   Bool_t  &icmpt,  Bool_t  &icms, 
					   Float_t &coneptsum) const 
{  
  //Search in cone around a candidate particle if it is isolated 
  Float_t phiC  = pCandidate->Phi() ;
  Float_t etaC = pCandidate->Eta() ;
  Float_t pt     = -100. ;
  Float_t eta   = -100.  ;
  Float_t phi    = -100.  ;
  Float_t rad   = -100 ;
  Int_t    n        = 0 ;
  TParticle * particle  = new TParticle;

  coneptsum = 0; 
  icmpt = kFALSE;
  icms  = kFALSE;

  //Check charged particles in cone.
  for(Int_t ipr = 0;ipr < plCTS->GetEntries() ; ipr ++ ){
    particle = dynamic_cast<TParticle *>(plCTS->At(ipr)) ;
    pt    = particle->Pt();
    eta  = particle->Eta();
    phi  = particle->Phi() ;
    
    //Check if there is any particle inside cone with pt larger than  fPtThreshold
    rad = TMath::Sqrt(TMath::Power((eta-etaC),2) +
		      TMath::Power((phi-phiC),2));
    if(rad<fConeSize){
      AliDebug(3,Form("charged in cone pt %f, phi %f, eta %f, R %f ",pt,phi,eta,rad));
      coneptsum+=pt;
      if(pt > fPtThreshold ) n++;
    }
  }// charged particle loop
  
  //Check neutral particles in cone.
  for(Int_t ipr = 0;ipr < plNe->GetEntries() ; ipr ++ ){
    if(ipr != index){//Do not count the candidate
      particle = dynamic_cast<TParticle *>(plNe->At(ipr)) ;
      pt    = particle->Pt();
      eta  = particle->Eta();
      phi  = particle->Phi() ;
      
      //Check if there is any particle inside cone with pt larger than  fPtThreshold
      rad = TMath::Sqrt(TMath::Power((eta-etaC),2) +
			TMath::Power((phi-phiC),2));
      if(rad<fConeSize){
	AliDebug(3,Form("charged in cone pt %f, phi %f, eta %f, R %f ",pt,phi,eta,rad));
	coneptsum+=pt;
	if(pt > fPtThreshold ) n++;
      }
    }
  }// neutral particle loop
  
  if(n == 0) 
    icmpt =  kTRUE ;
  if(coneptsum < fPtSumThreshold)
    icms  =  kTRUE ;

}

//__________________________________________________________________
void  AliAnaGammaDirect::MakeSeveralICAnalysis(TClonesArray * plCalo, TClonesArray * plCTS) 
{
  //Isolation Cut Analysis for both methods and different pt cuts and cones

  if (fICMethod != kSeveralIC)
    AliFatal("Remember to set in config file: directGamma->SetICMethod(kSeveralIC)");
  
  for(Int_t ipr = 0; ipr < plCalo->GetEntries() ; ipr ++ ){
    TParticle * pCandidate = dynamic_cast<TParticle *>(plCalo->At(ipr)) ;
    
    if(pCandidate->Pt() > fMinGammaPt){
      
      Bool_t  icPtThres   = kFALSE;
      Bool_t  icPtSum     = kFALSE;
      
      Float_t ptC      = pCandidate->Pt() ;
   
      fhNGamma->Fill(ptC);
      fhPhiGamma->Fill( ptC,pCandidate->Phi());
      fhEtaGamma->Fill(ptC,pCandidate->Eta());
    
      for(Int_t icone = 0; icone<fNCones; icone++){
	fConeSize=fConeSizes[icone] ;
	Float_t coneptsum = 0 ;
	for(Int_t ipt = 0; ipt<fNPtThres;ipt++){ 
	  fPtThreshold=fPtThresholds[ipt] ;
	  MakeIsolationCut(plCTS,plCalo, pCandidate, ipr, icPtThres, icPtSum,coneptsum);
	  AliDebug(1,Form("Candidate pt %f, pt in cone %f, Isolated? ICPt %d, ICSum %d",
			  pCandidate->Pt(), coneptsum, icPtThres, icPtSum));

	  fhPtThresIsolated[icone][ipt]->Fill(ptC); 
	}//pt thresh loop
	fhPtSumIsolated[icone]->Fill(ptC,coneptsum) ;
      }//cone size loop
    }//min pt candidate
  }//candidate loop
}

void AliAnaGammaDirect::Print(const Option_t * opt) const
{
  
  //Print some relevant parameters set for the analysis
  if(! opt)
    return;
  
  Info("Print", "%s %s", GetName(), GetTitle() ) ;
  
  printf("Min Gamma pT      =     %f\n",  fMinGammaPt) ;
  printf("IC method               =     %d\n", fICMethod) ; 
  printf("Cone Size               =     %f\n", fConeSize) ; 
   if(fICMethod == kPtIC) printf("pT threshold           =     %f\n", fPtThreshold) ;
   if(fICMethod == kSumPtIC) printf("pT sum threshold   =     %f\n", fPtSumThreshold) ;
   
  if(fICMethod == kSeveralIC){
    printf("N Cone Sizes               =     %d\n", fNCones) ; 
    printf("N pT thresholds           =     %d\n", fNPtThres) ;
    printf("Cone Sizes                  =    \n") ;
    for(Int_t i = 0; i < fNCones; i++)
      printf("   %f;",  fConeSizes[i]) ;
    printf("    \n") ;
    for(Int_t i = 0; i < fNPtThres; i++)
      printf("   %f;",  fPtThresholds[i]) ;
  }

  printf("    \n") ;
  
} 
