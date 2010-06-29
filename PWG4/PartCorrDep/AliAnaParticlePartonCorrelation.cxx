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
/* $Id:  $ */

//_________________________________________________________________________
// Class for the analysis of particle-parton correlations
// Particle (for example direct gamma) must be found in a previous analysis 
// -- Author: Gustavo Conesa (LNF-INFN) 
//////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---
//#include "Riostream.h"
#include "TH2F.h"
#include "TParticle.h"
#include "TClass.h"

//---- ANALYSIS system ----
#include "AliAnaParticlePartonCorrelation.h" 
#include "AliStack.h"  
#include "AliAODPWG4ParticleCorrelation.h"

  ClassImp(AliAnaParticlePartonCorrelation)
  

//____________________________________________________________________________
  AliAnaParticlePartonCorrelation::AliAnaParticlePartonCorrelation() : 
    AliAnaPartCorrBaseClass(),   
    fhDeltaEtaNearParton(0), fhDeltaPhiNearParton(0), 
    fhDeltaPtNearParton(0), fhPtRatNearParton(0),
    fhDeltaEtaAwayParton(0), fhDeltaPhiAwayParton(0), 
    fhDeltaPtAwayParton(0), fhPtRatAwayParton(0)
{
  //Default Ctor

  //Initialize parameters
  InitParameters();
}
/*
//____________________________________________________________________________
AliAnaParticlePartonCorrelation::AliAnaParticlePartonCorrelation(const AliAnaParticlePartonCorrelation & g) :   
  AliAnaPartCorrBaseClass(g),   
  fhDeltaEtaNearParton(g.fhDeltaEtaNearParton), fhDeltaPhiNearParton(g.fhDeltaPhiNearParton), 
  fhDeltaPtNearParton(g.fhDeltaPtNearParton), fhPtRatNearParton(g.fhPtRatNearParton),
  fhDeltaEtaAwayParton(g.fhDeltaEtaAwayParton), fhDeltaPhiAwayParton(g.fhDeltaPhiAwayParton), 
  fhDeltaPtAwayParton(g.fhDeltaPtAwayParton), fhPtRatAwayParton(g.fhPtRatAwayParton)
{
  // cpy ctor

}

//_________________________________________________________________________
AliAnaParticlePartonCorrelation & AliAnaParticlePartonCorrelation::operator = (const AliAnaParticlePartonCorrelation & source)
{
  // assignment operator

  if(this == &source)return *this;
  ((AliAnaPartCorrBaseClass *)this)->operator=(source);

  fhDeltaEtaAwayParton = source.fhDeltaEtaAwayParton;
  fhDeltaPhiAwayParton = source.fhDeltaPhiAwayParton;
  fhDeltaPtAwayParton = source.fhDeltaPtAwayParton;
  fhPtRatAwayParton = source.fhPtRatAwayParton;
  fhDeltaEtaNearParton = source.fhDeltaEtaNearParton;
  fhDeltaPhiNearParton = source.fhDeltaPhiNearParton;
  fhDeltaPtNearParton = source.fhDeltaPtNearParton;
  fhPtRatNearParton = source.fhPtRatNearParton;

  return *this;

}
*/

//________________________________________________________________________
TList *  AliAnaParticlePartonCorrelation::GetCreateOutputObjects()
{  
  // Create histograms to be saved in output file 

  TList * outputContainer = new TList() ; 
  outputContainer->SetName("ParticlePartonHistos") ; 
  
  fhDeltaPhiNearParton  = new TH2F
    ("DeltaPhiNearParton","#phi_{particle} - #phi_{parton} vs p_{T particle}",
     200,0,120,200,0,6.4); 
  fhDeltaPhiNearParton->SetYTitle("#Delta #phi");
  fhDeltaPhiNearParton->SetXTitle("p_{T particle} (GeV/c)");
  outputContainer->Add(fhDeltaPhiNearParton);
  
  fhDeltaEtaNearParton  = new TH2F
    ("DeltaEtaNearParton","#eta_{particle} - #eta_{parton} vs p_{T particle}",
     200,0,120,200,-2,2); 
  fhDeltaEtaNearParton->SetYTitle("#Delta #eta");
  fhDeltaEtaNearParton->SetXTitle("p_{T particle} (GeV/c)");
  outputContainer->Add(fhDeltaEtaNearParton);
  
  fhDeltaPtNearParton  = new TH2F
    ("DeltaPtNearParton","#p_{T particle} - #p_{T parton} vs p_{T particle}",
     200,0,120,100,-10,10); 
  fhDeltaPtNearParton->SetYTitle("#Delta #p_{T}");
  fhDeltaPtNearParton->SetXTitle("p_{T particle} (GeV/c)"); 
  outputContainer->Add(fhDeltaPtNearParton);
  
  fhPtRatNearParton  = new TH2F
    ("PtRatNearParton","#p_{T parton} / #p_{T particle} vs p_{T particle}",
     200,0,120,200,0,5); 
  fhPtRatNearParton->SetYTitle("ratio");
  fhPtRatNearParton->SetXTitle("p_{T particle} (GeV/c)");
  outputContainer->Add(fhPtRatNearParton);
  
  fhDeltaPhiAwayParton  = new TH2F
    ("DeltaPhiAwayParton","#phi_{particle} - #phi_{parton} vs p_{T particle}",
     200,0,120,200,0,6.4); 
  fhDeltaPhiAwayParton->SetYTitle("#Delta #phi");
  fhDeltaPhiAwayParton->SetXTitle("p_{T particle} (GeV/c)");
  outputContainer->Add(fhDeltaPhiAwayParton);
  
  fhDeltaEtaAwayParton  = new TH2F
    ("DeltaEtaAwayParton","#eta_{particle} - #eta_{parton} vs p_{T particle}",
     200,0,120,200,-2,2); 
  fhDeltaEtaAwayParton->SetYTitle("#Delta #eta");
  fhDeltaEtaAwayParton->SetXTitle("p_{T particle} (GeV/c)");
  outputContainer->Add(fhDeltaEtaAwayParton);
  
  fhDeltaPtAwayParton  = new TH2F
    ("DeltaPtAwayParton","#p_{T particle} - #p_{T parton} vs p_{T particle}",
     200,0,120,100,-10,10); 
  fhDeltaPtAwayParton->SetYTitle("#Delta #p_{T}");
  fhDeltaPtAwayParton->SetXTitle("p_{T particle} (GeV/c)"); 
  outputContainer->Add(fhDeltaPtAwayParton);
  
  fhPtRatAwayParton  = new TH2F
    ("PtRatAwayParton","#p_{T parton} / #p_{T particle} vs p_{T particle}",
     200,0,120,200,0,5); 
  fhPtRatAwayParton->SetYTitle("ratio");
  fhPtRatAwayParton->SetXTitle("p_{T particle} (GeV/c)");
  outputContainer->Add(fhPtRatAwayParton);
  
  return outputContainer;

}

//____________________________________________________________________________
void AliAnaParticlePartonCorrelation::InitParameters()
{
  
  //Initialize the parameters of the analysis.
  SetInputAODName("PWG4Particle");
  SetAODObjArrayName("Partons");  
  AddToHistogramsName("AnaPartonCorr_");

}

//__________________________________________________________________
void AliAnaParticlePartonCorrelation::Print(const Option_t * opt) const
{
  
  //Print some relevant parameters set for the analysis
  if(! opt)
    return;
  
  printf("**** Print %s %s ****\n", GetName(), GetTitle() ) ;
  AliAnaPartCorrBaseClass::Print(" ");

} 

//__________________________________________________________________
void  AliAnaParticlePartonCorrelation::MakeAnalysisFillAOD()  
{
  //Particle-Parton Correlation Analysis, create AODs
  //Add partons to the reference list of the trigger particle
  //Partons are considered those in the first eight possitions in the stack
  //being 0, and 1 the 2 protons, and 6 and 7 the outgoing final partons.
  if(!GetInputAODBranch()){
    printf("AliAnaParticlePartonCorrelation::MakeAnalysisFillAOD() - No input particles in AOD with name branch < %s > \n",GetInputAODName().Data());
    abort();	
  }
	
  if(strcmp(GetInputAODBranch()->GetClass()->GetName(), "AliAODPWG4ParticleCorrelation")){
	printf("AliAnaParticlePartonCorrelation::MakeAnalysisFillAOD() - Wrong type of AOD object, change AOD class name in input AOD: It should be <AliAODPWG4ParticleCorrelation> and not <%s> \n",GetInputAODBranch()->GetClass()->GetName());
	abort();
  }	
	
  if(GetDebug() > 1){
    printf("AliAnaParticlePartonCorrelation::MakeAnalysisFillAOD() - Begin fill AODs \n");
    printf("AliAnaParticlePartonCorrelation::MakeAnalysisFillAOD() - In particle branch aod entries %d\n", GetInputAODBranch()->GetEntriesFast());
  }
  
  //Loop on stored AOD particles
  Int_t naod = GetInputAODBranch()->GetEntriesFast();
  for(Int_t iaod = 0; iaod < naod ; iaod++){
    AliAODPWG4ParticleCorrelation* particle =  (AliAODPWG4ParticleCorrelation*) (GetInputAODBranch()->At(iaod));
    
    AliStack * stack =  GetMCStack() ;
    if(!stack){ 
      printf("AliAnaParticlePartonCorrelation::MakeAnalysisFillAOD() - No Stack available, STOP\n");
      abort();
    }
    if(stack->GetNtrack() < 8) {
      printf("AliAnaParticlePartonCorrelation::MakeAnalysisFillAOD() *** small number of particles, not a PYTHIA simulation? ***:  n tracks %d \n", stack->GetNprimary());
      continue ;
    }
    
    //Fill AOD reference only with partons
    TParticle * parton = new TParticle ;

    //Array with reference to partons, initialize
    TObjArray * objarray  = new TObjArray;

    for(Int_t ipr = 0;ipr < 8; ipr ++ ){
      parton = stack->Particle(ipr) ;
	  objarray->Add(parton);
    }//parton loop
	
    objarray->SetName(GetAODObjArrayName());
    if(objarray->GetEntriesFast() > 0) particle->AddObjArray(objarray);

  }//Aod branch loop
  
  if(GetDebug() > 1) printf("AliAnaParticlePartonCorrelation::MakeAnalysisFillAOD() - End fill AODs \n");
}

//__________________________________________________________________
void  AliAnaParticlePartonCorrelation::MakeAnalysisFillHistograms()  
{
  //Particle-Parton Correlation Analysis, fill histograms
  if(!GetInputAODBranch()){
    printf("AliAnaParticlePartonCorrelation::MakeAnalysisFillHistograms() - No input particles in AOD with name branch < %s > \n",GetInputAODName().Data());
    abort();	
  }
  if(GetDebug() > 1){
    printf("AliAnaParticlePartonCorrelation::MakeAnalysisFillHistograms() - Begin parton correlation analysis, fill histograms \n");
    printf("AliAnaParticlePartonCorrelation::MakeAnalysisFillHistograms() - In particle branch aod entries %d\n", GetInputAODBranch()->GetEntriesFast());
  }
  
  AliStack * stack =  GetMCStack() ;
  if(!stack) {
    printf("AliAnaParticlePartonCorrelation::MakeAnalysisFillHistograms() - No Stack available, STOP\n");
    abort();
  }
 
  //Loop on stored AOD particles
  Int_t naod = GetInputAODBranch()->GetEntriesFast();
  TParticle *  mom =new TParticle ;
  
  for(Int_t iaod = 0; iaod < naod ; iaod++){
    AliAODPWG4ParticleCorrelation* particle =  (AliAODPWG4ParticleCorrelation*) (GetInputAODBranch()->At(iaod));
    
    Float_t ptTrigg  = particle->Pt();
    Float_t phiTrigg = particle->Phi();
    Float_t etaTrigg = particle->Eta(); 
    Int_t imom = particle->GetLabel();
    Int_t iparent  = 2000;
    Int_t iawayparent = -1;

    TObjArray * objarray = particle->GetObjArray(GetAODObjArrayName());
    if(!(objarray) || (objarray->GetEntriesFast() < 7) ) {
      printf("AliAnaParticlePartonCorrelation::MakeAnalysisFillHistograms() - Reference list with partons not filled, STOP analysis\n");
      abort();
    }

    //Check and get indeces of mother and parton    
    if(imom < 8 ) iparent = imom ;   //mother is already a parton
    else if (imom <  stack->GetNtrack()) {
      mom =  stack->Particle(imom);
      iparent=mom->GetFirstMother();
      //cout<<" iparent "<<iparent<<endl;
      while(iparent > 7 ){
	mom = stack->Particle(iparent);
	imom = iparent ; //Mother label is of the inmediate parton daughter
	iparent = mom->GetFirstMother();
	//cout<<" while iparent "<<iparent<<endl;
      }   
    }
    
    if(GetDebug() > 1) printf("AliAnaParticlePartonCorrelation::MakeAnalysisFillHistograms() - N reference partons %d; labels:  mother %d, parent %d \n", objarray->GetEntriesFast(), imom, iparent);
    
    
    if(iparent < 0 || iparent > 8) { 
      if(GetDebug() > 0 ) printf("AliAnaParticlePartonCorrelation::MakeAnalysisFillHistograms() - Failed to find appropriate parton, index %d", iparent);
      continue ;
    }

    //Near parton is the parton that fragmented and created the mother    
    TParticle * nearParton = (TParticle*) objarray->At(iparent);
    Float_t  ptNearParton    = nearParton->Pt();
    Float_t  phiNearParton   = nearParton->Phi() ;
    Float_t  etaNearParton   = nearParton->Eta() ;
    
    fhDeltaEtaNearParton->Fill(ptTrigg,etaTrigg-etaNearParton);
    fhDeltaPhiNearParton->Fill(ptTrigg,phiTrigg-phiNearParton);
    fhDeltaPtNearParton->Fill(ptTrigg,ptTrigg-ptNearParton);
    fhPtRatNearParton->Fill(ptTrigg,ptNearParton/ptTrigg);
    
    if(iparent == 7) iawayparent =6;
    else if(iparent == 6) iawayparent =7;
    else{
      printf("AliAnaParticlePartonCorrelation::MakeAnalysisFillHistograms() - Parent parton is not final state, skip \n");
      continue;
    }

    //Away parton is the other final parton.
    TParticle * awayParton = (TParticle*) objarray->At(iawayparent);
    Float_t  ptAwayParton    = awayParton->Pt();
    Float_t  phiAwayParton   = awayParton->Phi() ;
    Float_t  etaAwayParton   = awayParton->Eta() ;
    fhDeltaEtaAwayParton->Fill(ptTrigg,etaTrigg-etaAwayParton);
    fhDeltaPhiAwayParton->Fill(ptTrigg,phiTrigg-phiAwayParton);
    fhDeltaPtAwayParton->Fill(ptTrigg,ptTrigg-ptAwayParton);
    fhPtRatAwayParton->Fill(ptTrigg,ptAwayParton/ptTrigg);
    
  }

  if(GetDebug() > 1) printf("AliAnaParticlePartonCorrelation::MakeAnalysisFillHistograms() - End fill histograms \n");
  
} 
