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
 * Revision 1.6  2007/08/17 12:40:04  schutz
 * New analysis classes by Gustavo Conesa
 *
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
#include "Riostream.h"
#include "TROOT.h"

// --- AliRoot system --- 
#include "AliAnaGammaDirect.h" 
#include "AliLog.h"

ClassImp(AliAnaGammaDirect)
  
//____________________________________________________________________________
  AliAnaGammaDirect::AliAnaGammaDirect() : 
    TObject(),
    fMinGammaPt(0.),
    fConeSize(0.),fPtThreshold(0.),fPtSumThreshold(0), 
    fICMethod(0), fAnaMC(0), fIsolatePi0(0),
    fhNGamma(0),fhPhiGamma(0),fhEtaGamma(0), fhConeSumPt(0), 
    fntuplePrompt(0),
    //kSeveralIC
    fNCones(0),fNPtThres(0), fConeSizes(),  fPtThresholds(), 
    fhPtThresIsolated(), fhPtSumIsolated(), fntSeveralIC()
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
  fAnaMC( g.fAnaMC), 
  fIsolatePi0(g.fIsolatePi0),
  fhNGamma(g.fhNGamma),fhPhiGamma(g.fhPhiGamma),
  fhEtaGamma(g.fhEtaGamma), fhConeSumPt(g.fhConeSumPt),    
  fntuplePrompt(g.fntuplePrompt),
  //kSeveralIC
  fNCones(g.fNCones),fNPtThres(g.fNPtThres), fConeSizes(),fPtThresholds(), 
  fhPtThresIsolated(), fhPtSumIsolated()
{
  // cpy ctor
  
  //kSeveralIC
  for(Int_t i = 0; i < fNCones ; i++){
    fConeSizes[i] =  g.fConeSizes[i];
    fhPtSumIsolated[i] = g.fhPtSumIsolated[i]; 
    fntSeveralIC[i]= g.fntSeveralIC[i];  
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
  fAnaMC = source.fAnaMC ;
  fIsolatePi0 =  source.fIsolatePi0 ;

  fhNGamma = source.fhNGamma ; 
  fhPhiGamma = source.fhPhiGamma ;
  fhEtaGamma = source.fhEtaGamma ;
  fhConeSumPt = source.fhConeSumPt ;

  fntuplePrompt = source.fntuplePrompt ;

  //kSeveralIC
  fNCones = source.fNCones ;
  fNPtThres = source.fNPtThres ; 
   
  for(Int_t i = 0; i < fNCones ; i++){
    fConeSizes[i] =  source.fConeSizes[i];
    fhPtSumIsolated[i] = source.fhPtSumIsolated[i] ;
    fntSeveralIC[i]= source.fntSeveralIC[i];  
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
  delete fhConeSumPt ;
  delete fntuplePrompt    ;  

  //kSeveralIC
  delete [] fhPtThresIsolated ;
  delete [] fhPtSumIsolated ;
  delete [] fntSeveralIC ;

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

  fhConeSumPt  = new TH2F
    ("ConePtSum","#Sigma p_{T}  in cone ",200,0,120,100,0,100);
  fhConeSumPt->SetYTitle("#Sigma p_{T}");
  fhConeSumPt->SetXTitle("p_{T #gamma} (GeV/c)");
  outputContainer->Add(fhConeSumPt) ;
  
  //NTUPLE
  fntuplePrompt = new TNtuple("ntuplePromptGamma", "Tree of prompt #gamma", "ptcluster:phicluster:etacluster:pdg:status:ptprimary:phiprimary:etaprimary:pdgprimary:statusprimary");
  outputContainer->Add(fntuplePrompt) ;

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
      
      sprintf(name,"nt_Cone_%d",icone);
      sprintf(title,"ntuple for cone size %d",icone);
      fntSeveralIC[icone] = new TNtuple(name, title, "ptcand:phicand:etacand:ptsum:type:ncone0:ncone1:ncone2:ncone3:ncone4:ncone5");
      outputContainer->Add(fntSeveralIC[icone]) ; 

      for(Int_t ipt = 0; ipt<fNPtThres;ipt++){ 
	sprintf(name,"PtThresIsol_Cone_%d_Pt%d",icone,ipt);
	sprintf(title,"Isolated candidate p_{T} distribution for cone size %d and p_{T}^{th} %d",icone,ipt);
	fhPtThresIsolated[icone][ipt]  = new TH1F(name, title,240,0,120);
	fhPtThresIsolated[icone][ipt]->SetXTitle("p_{T} (GeV/c)");
	outputContainer->Add(fhPtThresIsolated[icone][ipt]) ; 
      }//icone loop
    }//ipt loop
  }

  gROOT->cd();

  return outputContainer ;

}

//____________________________________________________________________________
void AliAnaGammaDirect::GetPromptGamma(TClonesArray * plCalo, TClonesArray * plCTS, TClonesArray * plPrimCalo, TParticle *pGamma, Bool_t &found) const 
{
  //Search for the prompt photon in Calorimeter with pt > fMinGammaPt

  Double_t pt = 0;
  Int_t n = 0;
  Int_t index = -1; 
  Float_t coneptsum = 0 ;

  for(Int_t ipr = 0;ipr < plCalo->GetEntries() ; ipr ++ ){
    TParticle * particle = dynamic_cast<TParticle *>(plCalo->At(ipr)) ;

    if((particle->Pt() > fMinGammaPt) && (particle->Pt() > pt) && 
       (particle->GetPdgCode() == 22 || (fIsolatePi0 && particle->GetPdgCode() == 111))){
      index = ipr ;
      pt = particle->Pt();
      pGamma->SetMomentum(particle->Px(),particle->Py(),particle->Pz(),particle->Energy());
      pGamma->SetStatusCode(particle->GetStatusCode());
      pGamma->SetPdgCode(particle->GetPdgCode());
      found  = kTRUE;
    }
  }
 
  //Do Isolation?
  if( ( fICMethod == kPtIC  ||  fICMethod == kSumPtIC )  && found)
    {
      Bool_t  icPtThres   = kFALSE;
      Bool_t  icPtSum     = kFALSE;
      MakeIsolationCut(plCTS,plCalo, pGamma, index,n,
		       icPtThres, icPtSum,coneptsum);
      if(fICMethod == kPtIC) //Pt thres method
	found = icPtThres ;
      if(fICMethod == kSumPtIC) //Pt cone sum method
	found = icPtSum ;
    }
  
  if(found){
    AliDebug(1,Form("Cluster with p_{T} larger than the pt treshold %f found in calorimeter ", fMinGammaPt)) ;
    AliDebug(1,Form("Gamma: pt %f, phi %f, eta %f", pGamma->Pt(),pGamma->Phi(),pGamma->Eta())) ;
    //Fill prompt gamma histograms 
    Float_t ptcluster = pGamma->Pt();
    Float_t phicluster = pGamma->Phi();
    Float_t etacluster = pGamma->Eta();
    Int_t statuscluster = pGamma->GetStatusCode();
    Int_t pdgcluster = pGamma->GetPdgCode();

    fhNGamma   ->Fill(ptcluster);
    fhPhiGamma ->Fill(ptcluster,phicluster);
    fhEtaGamma ->Fill(ptcluster,etacluster);
    fhConeSumPt->Fill(ptcluster,coneptsum);

    Float_t ptprimary = 0 ;
    Float_t phiprimary = 0 ;
    Float_t etaprimary = 0 ;
    Int_t pdgprimary = 0 ;
    Int_t statusprimary = 0 ;

    if(fAnaMC){
      TParticle * primary = dynamic_cast<TParticle *>(plPrimCalo->At(index)) ;
      ptprimary = primary->Pt();
      phiprimary = primary->Phi();
      etaprimary = primary->Eta();
      pdgprimary =  TMath::Abs(primary->GetPdgCode()) ;
      statusprimary = primary->GetStatusCode(); // = 2 means decay gamma!!!
 
      AliDebug(1, Form("Identified prompt Gamma pT %2.2f; Primary, pdg %d, pT %2.2f",ptcluster,pdgprimary,ptprimary));
      //printf("Identified prompt Gamma pT %2.2f; Primary, pdg %d, pT %2.2f \n",ptcluster,pdgprimary,ptprimary);
    }

    //Fill ntuple with cluster / MC data
    gROOT->cd();
    fntuplePrompt->Fill(ptcluster,phicluster,etacluster,pdgcluster,statuscluster,ptprimary,phiprimary, etaprimary,pdgprimary,statusprimary);
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
  fConeSize             = 0.4 ; 
  fPtThreshold         = 1.; 
  fPtSumThreshold  = 1.; 

  fICMethod = kNoIC; // 0 don't isolate, 1 pt thresh method, 2 cone pt sum method
  fAnaMC = kFALSE ;
  fIsolatePi0 = kFALSE ;
 //-----------kSeveralIC-----------------
  fNCones           = 5 ; 
  fNPtThres         = 6 ; 
  fConeSizes[0] = 0.1; fConeSizes[1] = 0.2; fConeSizes[2] = 0.3; fConeSizes[3] = 0.4; fConeSizes[4] = 0.5;
  fPtThresholds[0]=0.; fPtThresholds[1]=1.; fPtThresholds[2]=2.; fPtThresholds[3]=3.; fPtThresholds[4]=4.;fPtThresholds[5]=5.;

}

//__________________________________________________________________
void  AliAnaGammaDirect::MakeIsolationCut(TClonesArray * plCTS, 
					  TClonesArray * plNe, 
					  TParticle * pCandidate, 
					  Int_t index, Int_t & n,
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
  TParticle * particle  = new TParticle;

  n = 0 ;
  coneptsum = 0.; 
  icmpt = kFALSE;
  icms  = kFALSE;

  //Check charged particles in cone.
  for(Int_t ipr = 0;ipr < plCTS->GetEntries() ; ipr ++ ){
    particle = dynamic_cast<TParticle *>(plCTS->At(ipr)) ;
    pt    = particle->Pt();
    eta  = particle->Eta();
    phi  = particle->Phi() ;
    
    //Check if there is any particle inside cone with pt larger than  fPtThreshold
    rad = TMath::Sqrt((eta-etaC)*(eta-etaC)+ (phi-phiC)*(phi-phiC));
    if(rad < fConeSize){
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
      rad = TMath::Sqrt((eta-etaC)*(eta-etaC)+ (phi-phiC)*(phi-phiC));
      if(rad < fConeSize){
	AliDebug(3,Form("charged in cone pt %f, phi %f, eta %f, R %f ",pt,phi,eta,rad));
	coneptsum+=pt;
	if(pt > fPtThreshold ) n++;
      }
    }
  }// neutral particle loop
  
  if(n == 0) 
    icmpt =  kTRUE ;
  if(coneptsum <= fPtSumThreshold)
    icms  =  kTRUE ;

}

//__________________________________________________________________
void  AliAnaGammaDirect::MakeSeveralICAnalysis(TClonesArray * plCalo, TClonesArray * plCTS) 
{
  //Isolation Cut Analysis for both methods and different pt cuts and cones
  if (fICMethod != kSeveralIC)
    AliFatal("Remember to set in config file: directGamma->SetICMethod(kSeveralIC)");
  Int_t type = 0;

  //Search maximum energy photon in the event
  Double_t ptC = 0;
  Int_t index = -1; 
  TParticle * pCandidate = new TParticle();
  Bool_t found = kFALSE;
  
  for(Int_t ipr = 0;ipr < plCalo->GetEntries() ; ipr ++ ){
    TParticle * particle = dynamic_cast<TParticle *>(plCalo->At(ipr)) ;
    
    if((particle->Pt() > fMinGammaPt) && (particle->Pt() > ptC) && 
       (particle->GetPdgCode() == 22 ||  (fIsolatePi0 && particle->GetPdgCode() == 111))){
      index = ipr ;
      ptC = particle->Pt();
      pCandidate = particle ;
      found  = kTRUE;
    }
  }
  
  //If there is a large cluster, larger than threshold, study isolation cut
  if(found){
    
    fhNGamma->Fill(ptC);
    fhPhiGamma->Fill(ptC,pCandidate->Phi());
    fhEtaGamma->Fill(ptC,pCandidate->Eta());
    
    Int_t ncone[10][10];//[fNCones][fNPtThres];
    Bool_t  icPtThres   = kFALSE;
    Bool_t  icPtSum     = kFALSE;
    
    for(Int_t icone = 0; icone<fNCones; icone++){
      fConeSize=fConeSizes[icone] ;
      Float_t coneptsum = 0 ;
      for(Int_t ipt = 0; ipt<fNPtThres;ipt++){
	ncone[icone][ipt]=0;
	fPtThreshold=fPtThresholds[ipt] ;
	MakeIsolationCut(plCTS,plCalo, pCandidate, index,  
			 ncone[icone][ipt], icPtThres, icPtSum,coneptsum);
	AliDebug(1,Form("Candidate pt %f, pt in cone %f, Isolated? ICPt %d, ICSum %d",
			pCandidate->Pt(), coneptsum, icPtThres, icPtSum));
// 	if(ptC >15 && ptC < 25 && (icPtThres || icPtSum) && ipt ==0){
// 	  printf("R %0.1f, ptthres %1.1f, ptsum %1.1f, Candidate pt %2.2f,  eta %2.2f, phi %2.2f, pt in cone %2.2f, Isolated? ICPt %d, ICSum %d\n",
// 		 fConeSize,  fPtThreshold, fPtSumThreshold, ptC, pCandidate->Eta(), pCandidate->Phi()*TMath::RadToDeg(), coneptsum, icPtThres, icPtSum);
// 	  //cout<<"mother label "<<pCandidate->GetMother(0)<<endl;
// 	}
	
	fhPtThresIsolated[icone][ipt]->Fill(ptC); 
      }//pt thresh loop
      fhPtSumIsolated[icone]->Fill(ptC,coneptsum) ;
      gROOT->cd();
      fntSeveralIC[icone]->Fill(ptC,pCandidate->Phi(),pCandidate->Eta(), coneptsum,type,ncone[icone][0],ncone[icone][1],ncone[icone][2],ncone[icone][3],ncone[icone][4],ncone[icone][5]);
    }//cone size loop
  }//found high energy gamma in the event
}

//__________________________________________________________________
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
