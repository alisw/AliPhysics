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
    TObject(), fM(0),
    fInvMassMaxCut(0.), fInvMassMinCut(0.),
    fAngleMaxParam(),  fKeepNeutralMesonHistos(0), 
    fhAnglePairNoCut(0), fhAnglePairOpeningAngleCut(0), 
    fhAnglePairAllCut(0), 
    fhInvMassPairNoCut(0), fhInvMassPairOpeningAngleCut(0), 
    fhInvMassPairAllCut(0),
    fHistoNEBins(0),   fHistoEMax(0.),   fHistoEMin(0.),
    fHistoNPtBins(0),  fHistoPtMax(0.),  fHistoPtMin(0.),
    fHistoNAngleBins(0), fHistoAngleMax(0.), fHistoAngleMin(0.),
    fHistoNIMBins(0), fHistoIMMax(0.), fHistoIMMin(0.)
{
  //Default Ctor
  
  //Initialize parameters
  
  // kGammaHadron and kGammaJet 
  fAngleMaxParam.Set(4) ;
  fAngleMaxParam.Reset(0.);
  
  //Initialize parameters
  InitParameters();
}

//____________________________________________________________________________
AliNeutralMesonSelection::AliNeutralMesonSelection(const AliNeutralMesonSelection & g) :   
  TObject(), fM(g.fM),
  fInvMassMaxCut(g.fInvMassMaxCut), fInvMassMinCut(g.fInvMassMinCut),
  fAngleMaxParam(g.fAngleMaxParam),
  fKeepNeutralMesonHistos(g.fKeepNeutralMesonHistos),
  fhAnglePairNoCut(g. fhAnglePairNoCut), 
  fhAnglePairOpeningAngleCut(g. fhAnglePairOpeningAngleCut), 
  fhAnglePairAllCut(g. fhAnglePairAllCut), 
  fhInvMassPairNoCut(g.fhInvMassPairNoCut),  
  fhInvMassPairOpeningAngleCut(g.fhInvMassPairOpeningAngleCut), 
  fhInvMassPairAllCut(g.fhInvMassPairAllCut),
  fHistoNEBins(g.fHistoNEBins),   fHistoEMax(g.fHistoEMax),   fHistoEMin(g.fHistoEMin),
  fHistoNPtBins(g.fHistoNPtBins),   fHistoPtMax(g.fHistoPtMax),   fHistoPtMin(g.fHistoPtMin),
  fHistoNAngleBins(g.fHistoNAngleBins), fHistoAngleMax(g.fHistoAngleMax), fHistoAngleMin(g.fHistoAngleMin),
  fHistoNIMBins(g.fHistoNIMBins), fHistoIMMax(g.fHistoIMMax), fHistoIMMin(g.fHistoIMMin)
{
  // cpy ctor
}

//_________________________________________________________________________
AliNeutralMesonSelection & AliNeutralMesonSelection::operator = (const AliNeutralMesonSelection & source)
{
  // assignment operator
  
  if(this == &source)return *this;
  ((TObject *)this)->operator=(source);

  fM = source.fM ;
  fInvMassMaxCut = source.fInvMassMaxCut ; 
  fInvMassMinCut = source.fInvMassMinCut ;
  fAngleMaxParam = source.fAngleMaxParam ;
  fKeepNeutralMesonHistos = source.fKeepNeutralMesonHistos;
 
  fhAnglePairNoCut = source. fhAnglePairNoCut ; 
  fhAnglePairOpeningAngleCut = source. fhAnglePairOpeningAngleCut ; 
  fhAnglePairAllCut = source. fhAnglePairAllCut ; 
  fhInvMassPairNoCut = source.fhInvMassPairNoCut ; 
  fhInvMassPairOpeningAngleCut = source.fhInvMassPairOpeningAngleCut ; 
  fhInvMassPairAllCut = source.fhInvMassPairAllCut ; 
  
  fHistoNEBins = source.fHistoNEBins;   fHistoEMax = source.fHistoEMax;   fHistoEMin = source.fHistoEMin;
  fHistoNPtBins = source.fHistoNPtBins;   fHistoPtMax = source.fHistoPtMax;   fHistoPtMin = source.fHistoPtMin;
  fHistoNAngleBins = source.fHistoNAngleBins; fHistoAngleMax = source.fHistoAngleMax; fHistoAngleMin = source.fHistoAngleMin;
  fHistoNIMBins = source.fHistoNIMBins; fHistoIMMax = source.fHistoIMMax; fHistoIMMin = source.fHistoIMMin;
  
  return *this;
  
}

//____________________________________________________________________________
AliNeutralMesonSelection::~AliNeutralMesonSelection() 
{
  //dtor

//  if(!fKeepNeutralMesonHistos){
//    //Histograms initialized and filled but not passed to output container
//    //delete here, I am not sure this is correct
//    
//    if(fhAnglePairNoCut) delete fhAnglePairNoCut;
//    if(fhAnglePairOpeningAngleCut) delete fhAnglePairOpeningAngleCut; 
//    if(fhAnglePairAllCut) delete fhAnglePairAllCut;
//    if(fhInvMassPairNoCut) delete fhInvMassPairNoCut;
//    if(fhInvMassPairOpeningAngleCut) delete fhInvMassPairOpeningAngleCut;
//    if(fhInvMassPairAllCut) delete fhInvMassPairAllCut; 
//
//  }
  
}
//________________________________________________________________________
TList *  AliNeutralMesonSelection::GetCreateOutputObjects()
{  
  // Create histograms to be saved in output file and 
  // store them in outputContainer of the analysis class that calls this class.
  
  TList * outputContainer = new TList() ; 
  outputContainer->SetName("MesonDecayHistos") ; 
	if(fKeepNeutralMesonHistos) outputContainer->SetOwner(kFALSE);
	
  fhAnglePairNoCut  = new TH2F
    ("AnglePairNoCut",
     "Angle between all #gamma pair vs E_{#pi^{0}}",fHistoNEBins,fHistoEMin,fHistoEMax,fHistoNAngleBins,fHistoAngleMin,fHistoAngleMax); 
  fhAnglePairNoCut->SetYTitle("Angle (rad)");
  fhAnglePairNoCut->SetXTitle("E_{ #pi^{0}} (GeV)");
  
  fhAnglePairOpeningAngleCut  = new TH2F
    ("AnglePairOpeningAngleCut",
     "Angle between all #gamma pair (opening angle + azimuth cut) vs E_{#pi^{0}}"
     ,fHistoNEBins,fHistoEMin,fHistoEMax,fHistoNAngleBins,fHistoAngleMin,fHistoAngleMax); 
  fhAnglePairOpeningAngleCut->SetYTitle("Angle (rad)");
  fhAnglePairOpeningAngleCut->SetXTitle("E_{ #pi^{0}} (GeV)");
  
  fhAnglePairAllCut  = new TH2F
    ("AnglePairAllCut",
     "Angle between all #gamma pair (opening angle + inv mass cut+azimuth) vs E_{#pi^{0}}"
     ,fHistoNEBins,fHistoEMin,fHistoEMax,fHistoNAngleBins,fHistoAngleMin,fHistoAngleMax); 
  fhAnglePairAllCut->SetYTitle("Angle (rad)");
  fhAnglePairAllCut->SetXTitle("E_{ #pi^{0}} (GeV)");    
  
  //
  fhInvMassPairNoCut  = new TH2F
    ("InvMassPairNoCut","Invariant Mass of all #gamma pair vs E_{#pi^{0}}",
     fHistoNPtBins,fHistoPtMin,fHistoPtMax,fHistoNIMBins,fHistoIMMin,fHistoIMMax); 
  fhInvMassPairNoCut->SetYTitle("Invariant Mass (GeV/c^{2})");
  fhInvMassPairNoCut->SetXTitle("E_{ #pi^{0}} (GeV)");
  
  fhInvMassPairOpeningAngleCut  = new TH2F
    ("InvMassPairOpeningAngleCut",
     "Invariant Mass of #gamma pair (angle cut) vs E_{#pi^{0}}",
     fHistoNPtBins,fHistoPtMin,fHistoPtMax,fHistoNIMBins,fHistoIMMin,fHistoIMMax); 
  fhInvMassPairOpeningAngleCut->SetYTitle("Invariant Mass (GeV/c^{2})");
  fhInvMassPairOpeningAngleCut->SetXTitle(" E_{#pi^{0}}(GeV)");
  
  fhInvMassPairAllCut  = new TH2F
    ("InvMassPairAllCut",
     "Invariant Mass of #gamma pair (opening angle+invmass cut) vs E_{#pi^{0}}",
     fHistoNPtBins,fHistoPtMin,fHistoPtMax,fHistoNIMBins,fHistoIMMin,fHistoIMMax); 
  fhInvMassPairAllCut->SetYTitle("Invariant Mass (GeV/c^{2})");
  fhInvMassPairAllCut->SetXTitle("E_{#pi^{0}}(GeV)");
  
  outputContainer->Add(fhAnglePairNoCut) ; 
  outputContainer->Add(fhAnglePairOpeningAngleCut) ;
  outputContainer->Add(fhAnglePairAllCut) ; 
  
  outputContainer->Add(fhInvMassPairNoCut) ; 
  outputContainer->Add(fhInvMassPairOpeningAngleCut) ; 
  outputContainer->Add(fhInvMassPairAllCut) ; 
  
  return outputContainer;
}

//____________________________________________________________________________
void AliNeutralMesonSelection::InitParameters()
{
  
  //Initialize the parameters of the analysis.
  fKeepNeutralMesonHistos = kFALSE ;
  
  fAngleMaxParam.Set(4) ;
  fAngleMaxParam.AddAt(0.4,0);//={0.4,-0.25,0.025,-2e-4};
  fAngleMaxParam.AddAt(-0.25,1) ;
  fAngleMaxParam.AddAt(0.025,2) ;
  fAngleMaxParam.AddAt(-2e-4,3) ;

  fInvMassMaxCut  = 0.16 ;
  fInvMassMinCut  = 0.11 ;

  fM = 0.1349766;//neutralMeson mass, pi0
  
 //Histogrammes settings
  fHistoNEBins = 100 ;
  fHistoEMax   = 50 ;
  fHistoEMin   = 0.  ;  
  
  fHistoNPtBins = 240 ;
  fHistoPtMax   = 120 ;
  fHistoPtMin   = 0.  ;

  fHistoNAngleBins = 200 ;
  fHistoAngleMax   = 0.2;
  fHistoAngleMin   = 0.  ;

  fHistoNIMBins = 300 ;
  fHistoIMMax   = 0.5   ;
  fHistoIMMin   = 0.  ;  
}

//__________________________________________________________________________-
Bool_t AliNeutralMesonSelection::IsAngleInWindow(const Float_t angle,const Float_t e) const {
  //Check if the opening angle of the candidate pairs is inside 
  //our selection windowd

  Bool_t result = kFALSE;
  Double_t max =  fAngleMaxParam.At(0)*TMath::Exp(fAngleMaxParam.At(1)*e)
    +fAngleMaxParam.At(2)+fAngleMaxParam.At(3)*e;
  Double_t arg = (e*e-2*fM*fM)/(e*e);
  Double_t min = 100. ;
  if(arg>0.)
    min = TMath::ACos(arg);

  if((angle<max)&&(angle>=min))
    result = kTRUE;
 
  return result;
}

//____________________________________________________________________________
Bool_t  AliNeutralMesonSelection::SelectPair(TLorentzVector gammai, TLorentzVector gammaj)  
{  
  
  //Search for the neutral pion within selection cuts
  Bool_t goodpair = kFALSE ;
  
//  Double_t pt  = (gammai+gammaj).Pt();
  Double_t phi = (gammai+gammaj).Phi();
  if(phi < 0)
    phi+=TMath::TwoPi();
  Double_t invmass = (gammai+gammaj).M();
  Double_t angle   = gammaj.Angle(gammai.Vect());
  Double_t e       = (gammai+gammaj).E();
  
  //Fill histograms with no cuts applied.
  fhAnglePairNoCut->Fill(e,angle);
  fhInvMassPairNoCut->Fill(e,invmass);
    
  //Cut on the aperture of the pair
  if(IsAngleInWindow(angle,e)){
    fhAnglePairOpeningAngleCut     ->Fill(e,angle);
    fhInvMassPairOpeningAngleCut->Fill(e,invmass);
    //AliDebug(2,Form("Angle cut: pt %f, phi %f",pt,phi));
    
    //Cut on the invariant mass of the pair
    if((invmass>fInvMassMinCut) && (invmass<fInvMassMaxCut)){ 
      fhInvMassPairAllCut  ->Fill(e,invmass);
      fhAnglePairAllCut       ->Fill(e,angle);
      goodpair = kTRUE;
      //AliDebug(2,Form("IM cut: pt %f, phi %f",pt,phi));
    }//(invmass>0.125) && (invmass<0.145)
  }//Opening angle cut
  
  return goodpair; 
  
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
  printf("Angle selection param: \n");
  printf("p0 :     %f\n", fAngleMaxParam.At(0));
  printf("p1 :     %f\n", fAngleMaxParam.At(1));
  printf("p2 :     %f\n", fAngleMaxParam.At(2));
  printf("p3 :     %f\n", fAngleMaxParam.At(3));

  printf("Keep Neutral Meson Histos = %d\n",fKeepNeutralMesonHistos);
  
  if(fKeepNeutralMesonHistos){
    printf("Histograms: %3.1f < E  < %3.1f,  Nbin = %d\n", fHistoEMin,  fHistoEMax,  fHistoNEBins);
    printf("Histograms: %3.1f < pT < %3.1f,  Nbin = %d\n", fHistoPtMin,  fHistoPtMax,  fHistoNPtBins);
    printf("Histograms: %3.1f < angle < %3.1f, Nbin = %d\n", fHistoAngleMin, fHistoAngleMax, fHistoNAngleBins);
    printf("Histograms: %3.1f < IM < %3.1f, Nbin = %d\n", fHistoIMMin, fHistoIMMax, fHistoNIMBins);
    
  }
  
} 
