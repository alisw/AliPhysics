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
/* $Id: AliNeutralMesonSelection.cxx 27413 2008-07-18 13:28:12Z gconesab $ */

//_________________________________________________________________________
// Class that contains methods to select candidate pairs to neutral meson 
// 2 main selections, invariant mass around pi0 (also any other mass),
// apperture angle to distinguish from combinatorial.
//-- Author: Gustavo Conesa (INFN-LNF)

// --- ROOT system ---
#include <TLorentzVector.h>
#include <TH2.h>
#include <TList.h>
 
//---- AliRoot system ----
#include "AliNeutralMesonSelection.h" 

ClassImp(AliNeutralMesonSelection)
  
  
//____________________________________________________________________________
  AliNeutralMesonSelection::AliNeutralMesonSelection() : 
    TObject(),             fAsymmetryCut(1),                fUseAsymmetryCut(0),
    fM(0),                 fInvMassMaxCut(0.),              fInvMassMinCut(0.),           fInvMassMaxCutParam(),
    fAngleMaxParam(),      fUseAngleCut(0),                 
    fShiftMinAngle(0),     fKeepNeutralMesonHistos(0), 
    fhAnglePairNoCut(0),   fhAnglePairOpeningAngleCut(0),   fhAnglePairAsymmetryCut(0),   fhAnglePairAllCut(0), 
    fhInvMassPairNoCut(0), fhInvMassPairOpeningAngleCut(0), fhInvMassPairAsymmetryCut(0), fhInvMassPairAllCut(0),
    fhAsymmetryNoCut(0),   fhAsymmetryOpeningAngleCut(0),   fhAsymmetryAllCut(0),
    fHistoNEBins(0),       fHistoEMax(0.),                  fHistoEMin(0.),
    fHistoNAngleBins(0),   fHistoAngleMax(0.),              fHistoAngleMin(0.),
    fHistoNIMBins(0),      fHistoIMMax(0.),                 fHistoIMMin(0.)
{
  //Default Ctor
  
  //Initialize parameters
  InitParameters();
}

//________________________________________________________________________
TList *  AliNeutralMesonSelection::GetCreateOutputObjects()
{  
  // Create histograms to be saved in output file and 
  // store them in outputContainer of the analysis class that calls this class.
  
  TList * outputContainer = new TList() ; 
  outputContainer->SetName("MesonDecayHistos") ; 
  
  if(fKeepNeutralMesonHistos){
	  
	  outputContainer->SetOwner(kFALSE);
	  
	  fhAnglePairNoCut  = new TH2F
	  ("AnglePairNoCut",
	   "Angle between all #gamma pair vs E_{#pi^{0}}",fHistoNEBins,fHistoEMin,fHistoEMax,fHistoNAngleBins,fHistoAngleMin,fHistoAngleMax); 
	  fhAnglePairNoCut->SetYTitle("Angle (rad)");
	  fhAnglePairNoCut->SetXTitle("E_{ #pi^{0}} (GeV)");
	      
    fhAsymmetryNoCut  = new TH2F
	  ("AsymmetryNoCut","Asymmetry of all #gamma pair vs E_{#pi^{0}}",
	   fHistoNEBins,fHistoEMin,fHistoEMax,100,0,1); 
	  fhAsymmetryNoCut->SetYTitle("Asymmetry");
	  fhAsymmetryNoCut->SetXTitle("E_{ #pi^{0}} (GeV)");    
    
    fhInvMassPairNoCut  = new TH2F
	  ("InvMassPairNoCut","Invariant Mass of all #gamma pair vs E_{#pi^{0}}",
	   fHistoNEBins,fHistoEMin,fHistoEMax,fHistoNIMBins,fHistoIMMin,fHistoIMMax); 
	  fhInvMassPairNoCut->SetYTitle("Invariant Mass (GeV/c^{2})");
	  fhInvMassPairNoCut->SetXTitle("E_{ #pi^{0}} (GeV)");    
    
    outputContainer->Add(fhAnglePairNoCut) ; 
	  outputContainer->Add(fhAsymmetryNoCut) ; 
    outputContainer->Add(fhInvMassPairNoCut) ; 

    if(fUseAngleCut) {
      fhAnglePairOpeningAngleCut  = new TH2F
      ("AnglePairOpeningAngleCut",
       "Angle between all #gamma pair (opening angle) vs E_{#pi^{0}}"
       ,fHistoNEBins,fHistoEMin,fHistoEMax,fHistoNAngleBins,fHistoAngleMin,fHistoAngleMax); 
      fhAnglePairOpeningAngleCut->SetYTitle("Angle (rad)");
      fhAnglePairOpeningAngleCut->SetXTitle("E_{ #pi^{0}} (GeV)");
      
      fhAsymmetryOpeningAngleCut  = new TH2F
      ("AsymmetryOpeningAngleCut",
       "Asymmetry of #gamma pair (angle cut) vs E_{#pi^{0}}",
       fHistoNEBins,fHistoEMin,fHistoEMax,100,0,1); 
      fhAsymmetryOpeningAngleCut->SetYTitle("Asymmetry");
      fhAsymmetryOpeningAngleCut->SetXTitle(" E_{#pi^{0}}(GeV)");   
      
      fhInvMassPairOpeningAngleCut  = new TH2F
      ("InvMassPairOpeningAngleCut",
       "Invariant Mass of #gamma pair (angle cut) vs E_{#pi^{0}}",
       fHistoNEBins,fHistoEMin,fHistoEMax,fHistoNIMBins,fHistoIMMin,fHistoIMMax); 
      fhInvMassPairOpeningAngleCut->SetYTitle("Invariant Mass (GeV/c^{2})");
      fhInvMassPairOpeningAngleCut->SetXTitle(" E_{#pi^{0}}(GeV)");
      
      outputContainer->Add(fhAnglePairOpeningAngleCut) ;
      outputContainer->Add(fhAsymmetryOpeningAngleCut) ;
      outputContainer->Add(fhInvMassPairOpeningAngleCut) ;
    }
    
	  if(fUseAsymmetryCut) {
      fhAnglePairAsymmetryCut  = new TH2F
      ("AnglePairAsymmetryCut",
       "Angle between all #gamma pair (opening angle + asymetry cut) vs E_{#pi^{0}}"
       ,fHistoNEBins,fHistoEMin,fHistoEMax,fHistoNAngleBins,fHistoAngleMin,fHistoAngleMax); 
      fhAnglePairAsymmetryCut->SetYTitle("Angle (rad)");
      fhAnglePairAsymmetryCut->SetXTitle("E_{ #pi^{0}} (GeV)");
      
      fhInvMassPairAsymmetryCut  = new TH2F
      ("InvMassPairAsymmetryCut",
       "Invariant Mass of #gamma pair (opening angle + asymmetry) vs E_{#pi^{0}}",
       fHistoNEBins,fHistoEMin,fHistoEMax,fHistoNIMBins,fHistoIMMin,fHistoIMMax); 
      fhInvMassPairAsymmetryCut->SetYTitle("Invariant Mass (GeV/c^{2})");
      fhInvMassPairAsymmetryCut->SetXTitle("E_{#pi^{0}}(GeV)");      
      
      outputContainer->Add(fhAnglePairAsymmetryCut) ;
      outputContainer->Add(fhInvMassPairAsymmetryCut) ;
      
    }
    
	  fhAnglePairAllCut  = new TH2F
	  ("AnglePairAllCut",
	   "Angle between all #gamma pair (opening angle + asymmetry + inv mass cut) vs E_{#pi^{0}}"
	   ,fHistoNEBins,fHistoEMin,fHistoEMax,fHistoNAngleBins,fHistoAngleMin,fHistoAngleMax); 
	  fhAnglePairAllCut->SetYTitle("Angle (rad)");
	  fhAnglePairAllCut->SetXTitle("E_{ #pi^{0}} (GeV)");    
        
	  fhInvMassPairAllCut  = new TH2F
	  ("InvMassPairAllCut",
	   "Invariant Mass of #gamma pair (opening angle + asymmetry + invmass cut) vs E_{#pi^{0}}",
	   fHistoNEBins,fHistoEMin,fHistoEMax,fHistoNIMBins,fHistoIMMin,fHistoIMMax); 
	  fhInvMassPairAllCut->SetYTitle("Invariant Mass (GeV/c^{2})");
	  fhInvMassPairAllCut->SetXTitle("E_{#pi^{0}}(GeV)");
	  
	  fhAsymmetryAllCut  = new TH2F
	  ("AsymmetryAllCut",
	   "Asymmetry of #gamma pair (opening angle+invmass cut) vs E_{#pi^{0}}",
	   fHistoNEBins,fHistoEMin,fHistoEMax,100,0,1); 
	  fhAsymmetryAllCut->SetYTitle("Asymmetry");
	  fhAsymmetryAllCut->SetXTitle("E_{#pi^{0}}(GeV)");
    
    outputContainer->Add(fhAnglePairAllCut) ; 
	  outputContainer->Add(fhAsymmetryAllCut) ; 
    outputContainer->Add(fhInvMassPairAllCut) ;     

  }
  
  return outputContainer;

}

//____________________________________________________________________________
void AliNeutralMesonSelection::InitParameters()
{
  
  //Initialize the parameters of the analysis.
  fAngleMaxParam.Set(4) ;
  fAngleMaxParam.Reset(0.);

  SetParticle("Pi0");

  fShiftMinAngle   = 0.03;
  
  //Histogrammes settings
  fHistoNEBins     = 200 ;
  fHistoEMax       = 50  ;
  fHistoEMin       = 0.  ;  

  fHistoNAngleBins = 200 ;
  fHistoAngleMax   = 0.5 ;
  fHistoAngleMin   = 0.  ;

}

//__________________________________________________________________________-
Bool_t AliNeutralMesonSelection::IsAngleInWindow(const Float_t angle,const Float_t e) const 
{
 
  // Check if the opening angle of the candidate pairs is inside 
  // our selection window
  // Attention, only valid for Pi0, if needed for Eta need to revise max angle function or change parameters
  	
  Double_t max =  fAngleMaxParam.At(0)*TMath::Exp(fAngleMaxParam.At(1)*e)
    +fAngleMaxParam.At(2)+fAngleMaxParam.At(3)*e;
  Double_t arg = (e*e-2*fM*fM)/(e*e);
  Double_t min = 100. ;
  if(arg>0.)
    min = TMath::ACos(arg)-fShiftMinAngle;
  
  if((angle<max)&&(angle>=min)) return kTRUE  ;
  else                          return kFALSE ;

}

//____________________________________________________________________________
Bool_t  AliNeutralMesonSelection::SelectPair(TLorentzVector gammai, TLorentzVector gammaj, TString calo)  
{  
  
  //Search for the neutral pion within selection cuts
  
  //  Double_t pt  = (gammai+gammaj).Pt();
  Double_t phi = (gammai+gammaj).Phi();
  if(phi < 0)
    phi+=TMath::TwoPi();
  Double_t invmass = (gammai+gammaj).M();
  Double_t angle   = gammaj.Angle(gammai.Vect());
  Double_t e       = (gammai+gammaj).E();
  Double_t asy     = TMath::Abs((gammai-gammaj).E())/(gammai+gammaj).E();

  //Fill histograms with no cuts applied.
  if(fKeepNeutralMesonHistos){
	  fhAnglePairNoCut  ->Fill(e,angle);
	  fhInvMassPairNoCut->Fill(e,invmass);
	  fhAsymmetryNoCut  ->Fill(e,asy);
  }
  
  //Cut on the aperture of the pair
  if(fUseAngleCut){
    if(IsAngleInWindow(angle,e)){
      if(fKeepNeutralMesonHistos ){
        fhAnglePairOpeningAngleCut  ->Fill(e,angle);
        fhInvMassPairOpeningAngleCut->Fill(e,invmass);
        fhAsymmetryOpeningAngleCut  ->Fill(e,asy);
      }
      //AliDebug(2,Form("Angle cut: pt %f, phi %f",pt,phi));
    } else return kFALSE;
  }
  
  // Asymmetry cut
  if(fUseAsymmetryCut){
    if(fAsymmetryCut > asy){
      if(fKeepNeutralMesonHistos){
        fhInvMassPairAsymmetryCut->Fill(e,invmass);
        fhAnglePairAsymmetryCut  ->Fill(e,angle);
      }
    } else return kFALSE;
  }
  
  
  //Cut on the invariant mass of the pair
  
  Float_t invmassmaxcut = fInvMassMaxCut;
  if(calo=="EMCAL" && e > 6.){ // for EMCAL, pi0s, mass depends strongly with energy for e > 6, loose max cut
  
    invmassmaxcut = (fInvMassMaxCutParam[0]+fInvMassMaxCut)+fInvMassMaxCutParam[1]*e+fInvMassMaxCutParam[2]*e*e;
    //printf("e %f, max cut %f, p00 %f,p0 %f,p1 %f,p2 %f\n",
    //       e,invmassmaxcut,fInvMassMaxCut,fInvMassMaxCutParam[0],fInvMassMaxCutParam[1],fInvMassMaxCutParam[2]);
  }
  
  if((invmass > fInvMassMinCut) && (invmass < invmassmaxcut)){ 
    if(fKeepNeutralMesonHistos){
      fhInvMassPairAllCut->Fill(e,invmass);
      fhAnglePairAllCut  ->Fill(e,angle);
      fhAsymmetryAllCut  ->Fill(e,asy);
    }      
    
    //AliDebug(2,Form("IM cut: pt %f, phi %f",pt,phi));
    return kTRUE;
    
  }//(invmass>0.125) && (invmass<0.145)
  else return kFALSE;
  
}

//____________________________________________________________________________
void  AliNeutralMesonSelection::SetParticle(TString particleName){
  // Set some default parameters for selection of pi0 or eta
  
  if(particleName=="Pi0"){
    
    fM               = 0.135 ; // GeV
    fInvMassMaxCut   = 0.16  ; // GeV
    fInvMassMinCut   = 0.11  ; // GeV

    fInvMassMaxCutParam[0] = 0.0   ;
    fInvMassMaxCutParam[1] =-7.e-5 ;
    fInvMassMaxCutParam[2] = 8.e-5 ;    
         
    fAngleMaxParam.AddAt( 0.40, 0) ;
    fAngleMaxParam.AddAt(-0.25, 1) ;
    fAngleMaxParam.AddAt( 0.025,2) ; //for pi0 shift, for eta maybe 0.09 
    fAngleMaxParam.AddAt(-2.e-4,3) ;
    
    fHistoNIMBins    = 150 ;
    fHistoIMMax      = 0.3 ;
    fHistoIMMin      = 0.  ;  
    
  }else if(particleName=="Eta"){
  
    fM               = 0.547 ; // GeV
    fInvMassMaxCut   = 0.590 ; // GeV
    fInvMassMinCut   = 0.510 ; // GeV
    
    fInvMassMaxCutParam[0] = 0.0 ;
    fInvMassMaxCutParam[1] = 0.0 ;
    fInvMassMaxCutParam[2] = 0.0 ;    
    
    fAngleMaxParam.AddAt( 0.40, 0) ; // Same as pi0
    fAngleMaxParam.AddAt(-0.25, 1) ; // Same as pi0
    fAngleMaxParam.AddAt( 0.1,2) ;   // Shifted with respect to pi0
    fAngleMaxParam.AddAt(-2.e-4,3) ; // Same as pi0
  
    fHistoNIMBins    = 200  ; // GeV
    fHistoIMMax      = 0.75 ; // GeV
    fHistoIMMin      = 0.35 ; // GeV 
    
  }
  else 
    printf("AliAnaNeutralMesonSelection::SetParticle(%s) *** Particle NOT defined (Pi0 or Eta), Pi0 settings by default *** \n",particleName.Data());
  
  
}

//__________________________________________________________________
void AliNeutralMesonSelection::Print(const Option_t * opt) const
{

  //Print some relevant parameters set for the analysis
  if(! opt)
    return;
  
  printf("**** Print %s %s ****\n", GetName(), GetTitle() ) ;

  printf("mass : %f  \n", fM );
  printf("Invariant mass limits : %f < m < %f \n", fInvMassMinCut , fInvMassMinCut );
  
  printf("Use asymmetry cut? : %d ; A < %f \n", fUseAngleCut, fAsymmetryCut );
  
  printf("Use angle cut? : %d  \n", fUseAngleCut );
  if(fUseAngleCut){
    printf("Angle selection param: \n");
    printf("p0 :     %f\n", fAngleMaxParam.At(0));
    printf("p1 :     %f\n", fAngleMaxParam.At(1));
    printf("p2 :     %f\n", fAngleMaxParam.At(2));
    printf("p3 :     %f\n", fAngleMaxParam.At(3));
    printf("Min angle shift : %1.2f\n", fShiftMinAngle);
  }
  
  printf("Keep Neutral Meson Histos = %d\n",fKeepNeutralMesonHistos);
  
  if(fKeepNeutralMesonHistos){
    printf("Histograms: %3.1f < E  < %3.1f,  Nbin = %d\n",   fHistoEMin,     fHistoEMax,     fHistoNEBins);
    printf("Histograms: %3.1f < angle < %3.1f, Nbin = %d\n", fHistoAngleMin, fHistoAngleMax, fHistoNAngleBins);
    printf("Histograms: %3.1f < IM < %3.1f, Nbin = %d\n",    fHistoIMMin,    fHistoIMMax,    fHistoNIMBins);    
  }
  
} 
