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
 * Revision 1.5  2007/08/17 12:40:04  schutz
 * New analysis classes by Gustavo Conesa
 *
 * Revision 1.4.4.2  2007/07/26 10:32:09  schutz
 * new analysis classes in the the new analysis framework
 *
 *
 */

//_________________________________________________________________________
// Class for the analysis of gamma - hadron correlations
//*-- Author: Gustavo Conesa (LNF-INFN) 
//////////////////////////////////////////////////////////////////////////////


// --- ROOT system ---
#include "Riostream.h"

//---- AliRoot system ----
#include "AliLog.h"
#include "AliNeutralMesonSelection.h" 
#include "AliAnaGammaHadron.h" 

ClassImp(AliAnaGammaHadron)


//____________________________________________________________________________
  AliAnaGammaHadron::AliAnaGammaHadron() : 
    AliAnaGammaCorrelation(),
    fhPhiCharged(0), fhPhiNeutral(0), fhEtaCharged(0), fhEtaNeutral(0), 
    fhDeltaPhiGammaCharged(0),  fhDeltaPhiGammaNeutral(0), 
    fhDeltaEtaGammaCharged(0), fhDeltaEtaGammaNeutral(0),
    fhDeltaPhiChargedPt(0),  
    fhCorrelationGammaNeutral(0), fhCorrelationGammaCharged(0)
{
  //Default Ctor
  
  SetCorrelationType(kHadron);
  //Initialize parameters
  InitParameters();
}

//____________________________________________________________________________
AliAnaGammaHadron::AliAnaGammaHadron(const AliAnaGammaHadron & g) :   
  AliAnaGammaCorrelation(g),
  fhPhiCharged(g.fhPhiCharged), fhPhiNeutral(g.fhPhiNeutral), 
  fhEtaCharged(g.fhEtaCharged), fhEtaNeutral(g.fhEtaNeutral), 
  fhDeltaPhiGammaCharged(g.fhDeltaPhiGammaCharged),  
  fhDeltaPhiGammaNeutral(g.fhDeltaPhiGammaNeutral), 
  fhDeltaEtaGammaCharged(g.fhDeltaEtaGammaCharged), 
  fhDeltaEtaGammaNeutral(g.fhDeltaEtaGammaNeutral), 
  fhDeltaPhiChargedPt(g.fhDeltaPhiChargedPt), 
  fhCorrelationGammaNeutral(g.fhCorrelationGammaNeutral), 
  fhCorrelationGammaCharged(g.fhCorrelationGammaCharged)
{
  // cpy ctor

}

//_________________________________________________________________________
AliAnaGammaHadron & AliAnaGammaHadron::operator = (const AliAnaGammaHadron & source)
{
  // assignment operator

  if(this == &source)return *this;
  ((AliAnaGammaCorrelation *)this)->operator=(source);
  
  fhPhiCharged = source.fhPhiCharged ; fhPhiNeutral = source.fhPhiNeutral ; 
  fhEtaCharged = source.fhEtaCharged ; fhEtaNeutral = source.fhEtaNeutral ; 
  fhDeltaPhiGammaCharged = source.fhDeltaPhiGammaCharged ;  
  fhDeltaPhiGammaNeutral = source.fhDeltaPhiGammaNeutral ; 
  fhDeltaEtaGammaCharged = source.fhDeltaEtaGammaCharged ; 
  fhDeltaEtaGammaNeutral = source.fhDeltaEtaGammaNeutral ; 
  fhDeltaPhiChargedPt = source.fhDeltaPhiChargedPt ;

  fhCorrelationGammaNeutral = source.fhCorrelationGammaNeutral ; 
  fhCorrelationGammaCharged = source.fhCorrelationGammaCharged ; 

  return *this;

}

//____________________________________________________________________________
AliAnaGammaHadron::~AliAnaGammaHadron() 
{
  
  delete fhPhiCharged  ;  
  delete fhPhiNeutral   ; 
  delete fhEtaCharged  ; 
  delete fhEtaNeutral  ; 
  delete fhDeltaPhiGammaCharged  ;  
  delete fhDeltaPhiGammaNeutral   ; 
  delete fhDeltaEtaGammaCharged  ; 
  delete fhDeltaEtaGammaNeutral  ; 
  delete fhDeltaPhiChargedPt  ;

  delete fhCorrelationGammaNeutral  ; 
  delete fhCorrelationGammaCharged  ;
 
}


//________________________________________________________________________
TList *  AliAnaGammaHadron::GetCreateOutputObjects()
{  

  // Create histograms to be saved in output file and 
  // store them in fOutputContainer
  TList * outputContainer = new TList() ; 
  outputContainer->SetName("GammaCorrelationHistos") ; 
  
  fhPhiCharged  = new TH2F
    ("PhiCharged","#phi_{#pi^{#pm}}  vs p_{T #gamma}",
     120,0,120,120,0,7); 
  fhPhiCharged->SetYTitle("#phi_{#pi^{#pm}} (rad)");
  fhPhiCharged->SetXTitle("p_{T #gamma} (GeV/c)");
  
  fhEtaCharged  = new TH2F
    ("EtaCharged","#eta_{#pi^{#pm}}  vs p_{T #gamma}",
     120,0,120,120,-1,1); 
  fhEtaCharged->SetYTitle("#eta_{#pi^{#pm}} (rad)");
  fhEtaCharged->SetXTitle("p_{T #gamma} (GeV/c)");
  
  fhDeltaPhiGammaCharged  = new TH2F
    ("DeltaPhiGammaCharged","#phi_{#gamma} - #phi_{charged #pi} vs p_{T #gamma}",
     200,0,120,200,0,6.4); 
  fhDeltaPhiGammaCharged->SetYTitle("#Delta #phi");
  fhDeltaPhiGammaCharged->SetXTitle("p_{T #gamma} (GeV/c)");

  fhDeltaPhiChargedPt  = new TH2F
    ("DeltaPhiChargedPt","#phi_{#gamma} - #phi_{charged #pi} vs p_{T #pi}",
     200,0,120,200,0,6.4);
  fhDeltaPhiChargedPt->SetYTitle("#Delta #phi");
  fhDeltaPhiChargedPt->SetXTitle("p_{T #pi} (GeV/c)");

  fhDeltaEtaGammaCharged  = new TH2F
    ("DeltaEtaGammaCharged","#eta_{#gamma} - #eta_{#pi^{#pm}} vs p_{T #gamma}",
     200,0,120,200,-2,2); 
  fhDeltaEtaGammaCharged->SetYTitle("#Delta #eta");
  fhDeltaEtaGammaCharged->SetXTitle("p_{T #gamma} (GeV/c)");

  fhCorrelationGammaCharged  = 
    new TH2F("CorrelationGammaCharged","z_{#gamma #pi} = p_{T #pi^{#pm}} / p_{T #gamma}",
	     240,0.,120.,1000,0.,1.2); 
  fhCorrelationGammaCharged->SetYTitle("z_{#gamma #pi}");
  fhCorrelationGammaCharged->SetXTitle("p_{T #gamma}");
  
  outputContainer->Add(fhPhiCharged) ;
  outputContainer->Add(fhEtaCharged) ;
  outputContainer->Add(fhDeltaPhiGammaCharged) ; 
  outputContainer->Add(fhDeltaEtaGammaCharged) ;
  outputContainer->Add(fhCorrelationGammaCharged) ;
  outputContainer->Add(fhDeltaPhiChargedPt) ;

  if(!AreJetsOnlyInCTS()){
    //---- kHadron and kJetLeadCone ----
    fhPhiNeutral  = new TH2F
      ("PhiNeutral","#phi_{#pi^{0}}  vs p_{T #gamma}",
       120,0,120,120,0,7); 
    fhPhiNeutral->SetYTitle("#phi_{#pi^{0}} (rad)");
    fhPhiNeutral->SetXTitle("p_{T #gamma} (GeV/c)");
    
    fhEtaNeutral  = new TH2F
      ("EtaNeutral","#eta_{#pi^{0}}  vs p_{T #gamma}",
       120,0,120,120,-1,1); 
    fhEtaNeutral->SetYTitle("#eta_{#pi^{0}} (rad)");
    fhEtaNeutral->SetXTitle("p_{T #gamma} (GeV/c)");
    
    fhDeltaPhiGammaNeutral  = new TH2F
      ("DeltaPhiGammaNeutral","#phi_{#gamma} - #phi_{#pi^{0}} vs p_{T #gamma}",
       200,0,120,200,0,6.4); 
    fhDeltaPhiGammaNeutral->SetYTitle("#Delta #phi");
    fhDeltaPhiGammaNeutral->SetXTitle("p_{T #gamma} (GeV/c)");
    
    fhDeltaEtaGammaNeutral  = new TH2F
      ("DeltaEtaGammaNeutral","#eta_{#gamma} - #eta_{#pi^{#pm}} vs p_{T #gamma}",
       200,0,120,200,-2,2); 
    fhDeltaEtaGammaNeutral->SetYTitle("#Delta #eta");
    fhDeltaEtaGammaNeutral->SetXTitle("p_{T #gamma} (GeV/c)");
    
    fhCorrelationGammaNeutral  = 
      new TH2F("CorrelationGammaNeutral","z_{#gamma #pi} = p_{T #pi^{0}} / p_{T #gamma}",
	       240,0.,120.,1000,0.,1.2); 
    fhCorrelationGammaNeutral->SetYTitle("z_{#gamma #pi}");
    fhCorrelationGammaNeutral->SetXTitle("p_{T #gamma}");

    outputContainer->Add(fhPhiNeutral) ;  
    outputContainer->Add(fhEtaNeutral) ;   
    outputContainer->Add(fhDeltaPhiGammaNeutral) ; 
    outputContainer->Add(fhDeltaEtaGammaNeutral) ; 
    outputContainer->Add(fhCorrelationGammaNeutral) ;
  }

  SetOutputContainer(outputContainer);
  
  return outputContainer;
}

 //____________________________________________________________________________
void AliAnaGammaHadron::InitParameters()
{
 
  //Initialize the parameters of the analysis.

  SetMinPtHadron(2.)   ;
  SetDeltaPhiCutRange(1.5,4.5);
  SetJetsOnlyInCTS(kFALSE) ;

}

//__________________________________________________________________
void AliAnaGammaHadron::Print(const Option_t * opt) const
{

  //Print some relevant parameters set for the analysis
  if(! opt)
    return;
  
  Info("Print", "%s %s", GetName(), GetTitle() ) ;
  printf("Correlation analysis           =     %d\n", kHadron) ;
  printf("pT Hadron       >    %f\n", GetMinPtHadron()) ; 
  printf("Phi gamma-Hadron      <     %f\n", GetDeltaPhiMaxCut()) ; 
  printf("Phi gamma-Hadron      >     %f\n", GetDeltaPhiMinCut()) ;
  

} 

//____________________________________________________________________________
void  AliAnaGammaHadron::MakeGammaCorrelation(TParticle * pGamma, TClonesArray * plCTS, TClonesArray * plCalo)  
{  
  //Gamma Hadron Correlation Analysis
  AliDebug(2, "Make gamma-hadron correlation");

  MakeGammaChargedCorrelation(pGamma, plCTS);  
  if(!AreJetsOnlyInCTS())
  MakeGammaNeutralCorrelation(pGamma, plCalo);
    

}

//____________________________________________________________________________
void  AliAnaGammaHadron::MakeGammaChargedCorrelation(TParticle * pGamma, TClonesArray * pl)
{  
  //Gamma Charged Hadron Correlation Analysis
  AliDebug(2,"Make gamma-charged hadron correlation");

  Double_t ptg  = pGamma->Pt();
  Double_t phig = pGamma->Phi();
  Double_t pt    = -100.;
  Double_t rat   = -100.; 
  Double_t phi   = -100. ;

  for(Int_t ipr = 0;ipr < pl->GetEntries() ; ipr ++ ){
    
    TParticle * particle = dynamic_cast<TParticle *>(pl->At(ipr)) ;

    pt    = particle->Pt();
    rat   = pt/ptg ;
    phi   = particle->Phi() ;
    
    AliDebug(3,Form("pt %f, phi %f, phi gamma %f. Cuts:  delta phi min %f,  max%f, pT min %f",pt,phi,phig,GetDeltaPhiMinCut(),GetDeltaPhiMaxCut(),GetMinPtHadron()));
    
    fhEtaCharged->Fill(ptg,particle->Eta());
    fhPhiCharged->Fill(ptg,phi);
    fhDeltaEtaGammaCharged->Fill(ptg,pGamma->Eta()-particle->Eta());
    fhDeltaPhiGammaCharged->Fill(ptg,phig-phi);
    fhDeltaPhiChargedPt->Fill(pt,phig-phi);

    //Selection within angular and energy limits
    if(((phig-phi)> GetDeltaPhiMinCut()) && ((phig-phi)<GetDeltaPhiMaxCut()) && pt > GetMinPtHadron()){
      AliDebug(2,Form("Selected: pt %f, phi %f",pt,phi));
      fhCorrelationGammaCharged->Fill(ptg,rat);
    } 
  }//particle loop
}

//____________________________________________________________________________
void  AliAnaGammaHadron::MakeGammaNeutralCorrelation(TParticle * pGamma, TClonesArray * pl)  
{  
  //Gamma Neutral Hadron Correlation Analysis
  AliDebug(2,"Make gamma-neutral hadron correlation");

  Double_t pt = -100.;
  Double_t rat = -100.; 
  Double_t phi = -100. ;
  Double_t ptg  = pGamma->Pt();
  Double_t phig = pGamma->Phi();
  
  TIter next(pl);
  TParticle * particlei = 0;
  TParticle * particlej = 0;
  TLorentzVector gammai;
  TLorentzVector gammaj;

  Int_t iPrimary = -1;
  Int_t ksPdg = 0;
  Int_t jPrimary=-1;
  
  while ( (particlei = (TParticle*)next()) ) {
    iPrimary++;	  
    ksPdg = particlei->GetPdgCode();
    AliDebug(2, Form("neutral particles opposite to gamma: pt %f, pdg %d", particlei->Pt(),ksPdg));
    //2 gamma overlapped, found with PID
    if(ksPdg == 111){ 
      pt  = particlei->Pt();
      rat = pt/ptg ;
      phi = particlei->Phi() ;
      fhEtaNeutral->Fill(ptg,particlei->Eta());
      fhPhiNeutral->Fill(ptg,phi);
      fhDeltaEtaGammaNeutral->Fill(ptg,pGamma->Eta()-particlei->Eta());
      fhDeltaPhiGammaNeutral->Fill(ptg,phig-phi);
      
      //Selection within angular and energy limits
      if( (phig-phi)>GetDeltaPhiMinCut() && (phig-phi)<GetDeltaPhiMaxCut() && pt > GetMinPtHadron()){
	fhCorrelationGammaNeutral ->Fill(ptg,rat);
	AliDebug(2,Form("Selected pi0: pt %f, phi %f",pt,phi));
      }// cuts
    }// pdg = 111

    //Make invariant mass analysis
    else if(ksPdg == 22){//  gamma i
      
      //Search the photon companion in case it comes from  a Pi0 decay
      //Apply several cuts to select the good pair;
      particlei->Momentum(gammai);
       jPrimary=-1;
      TIter next2(pl);
      while ( (particlej = (TParticle*)next2()) ) {
 	jPrimary++;
 	if(jPrimary>iPrimary){
 	  ksPdg = particlej->GetPdgCode();
	  particlej->Momentum(gammaj);
 	  if(ksPdg == 22 ){

	    phi = (gammai+gammaj).Phi();
	    if(phi < 0)
	      phi+=TMath::TwoPi();
	    rat          = (gammai+gammaj).Pt()/ptg ;
	    	    
	    //Fill histograms
	    fhEtaNeutral->Fill(ptg,(gammai+gammaj).Eta());
	    fhPhiNeutral->Fill(ptg,phi);
	    fhDeltaEtaGammaNeutral->Fill(ptg,pGamma->Eta()-(gammai+gammaj).Eta());
	    fhDeltaPhiGammaNeutral->Fill(ptg,phig-phi);

	    //Select good pair (good phit, pt cuts, aperture and invariant mass)
	    if(GetNeutralMesonSelection()->SelectPair(pGamma, gammai, gammaj)){
	      AliDebug(2,Form("Selected gamma pair: pt %f, phi %f",(gammai+gammaj).Pt(),phi));
	      //correlation histogram
	      fhCorrelationGammaNeutral ->Fill(ptg,rat);
	    }//pair selection
	  }//if pair of gammas
	}//jPrimary>iPrimary
      }//while
    }// if pdg = 22
  }//while
  
}
