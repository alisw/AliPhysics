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

/* History of cvs commits:
 *
 * $Log$
 *
 *
 */

//_________________________________________________________________________
// Class for the analysis of gamma-jet correlations:
// 1)Take the trigger particle stored in AliAODParticleCorrelation,
// 2) Search for the highest pt leading particle opposite to the photon within a phi, pt window
// 3) Take all particles around leading in a cone R with pt larger than threshold and construct the jet
//
//  Class created from old AliPHOSGammaJet 
//  (see AliRoot versions previous Release 4-09)
//
//*-- Author: Gustavo Conesa (LNF-INFN) 
//////////////////////////////////////////////////////////////////////////////


// --- ROOT system ---

//---- AliRoot system ----
#include "AliNeutralMesonSelection.h"
#include "AliAnaParticleJetLeadingCone.h"  
#include "AliLog.h"

ClassImp(AliAnaParticleJetLeadingCone)


//____________________________________________________________________________
  AliAnaParticleJetLeadingCone::AliAnaParticleJetLeadingCone() : 
    AliAnaBaseClass(),  fPbPb(kFALSE),     
    fSeveralConeAndPtCuts(0),  
    fJetCTSRatioMaxCut(0.), fJetCTSRatioMinCut(0.), fJetRatioMaxCut(0.),
    fJetRatioMinCut(0.), 
    fJetNCone(0),fJetNPt(0), fJetCone(0), 
    fJetPtThreshold(0),fJetPtThresPbPb(0),
    fPtJetSelectionCut(0.0), fSelect(0),
    fhPhiCharged(0), fhPhiNeutral(0), fhEtaCharged(0), fhEtaNeutral(0), 
    fhDeltaPhiGammaCharged(0),  fhDeltaPhiGammaNeutral(0), 
    fhDeltaEtaGammaCharged(0), fhDeltaEtaGammaNeutral(0), 
    fhAnglePairLeading(), fhInvMassPairLeading(), 
    fhChargedRatio(0), fhNeutralRatio (0), 
    fhNBkg (0), fhNLeading(0), fhNJet(0), fhJetRatio(0), fhJetPt (0), 
    fhBkgRatio (0), fhBkgPt(0),  fhJetFragment(0), fhBkgFragment(0), 
    fhJetPtDist(0),  fhBkgPtDist(0) 
{
  //Default Ctor

  //Initialize parameters

  SetCorrelationType(kJetLeadCone);

  for(Int_t i = 0; i<10; i++){
    fJetCones[i]         = 0.0 ;
    fJetNameCones[i]     = ""  ;
    fJetPtThres[i]      = 0.0 ;
    fJetNamePtThres[i]  = ""  ;
    if( i < 6 ){
      fJetXMin1[i]     = 0.0 ;
      fJetXMin2[i]     = 0.0 ;
      fJetXMax1[i]     = 0.0 ;
      fJetXMax2[i]     = 0.0 ;
      fBkgMean[i]      = 0.0 ;
      fBkgRMS[i]       = 0.0 ;
      if( i < 2 ){
	fJetE1[i]        = 0.0 ;
	fJetE2[i]        = 0.0 ;
	fJetSigma1[i]    = 0.0 ;
	fJetSigma2[i]    = 0.0 ;
      }
    }
  }

  InitParameters();
}

//____________________________________________________________________________
AliAnaParticleJetLeadingCone::AliAnaParticleJetLeadingCone(const AliAnaParticleJetLeadingCone & g) :   
  AliAnaBaseClass(g), fPbPb(g.fPbPb), 
  fSeveralConeAndPtCuts(g.fSeveralConeAndPtCuts), 
  fJetCTSRatioMaxCut(g.fJetCTSRatioMaxCut),
  fJetCTSRatioMinCut(g.fJetCTSRatioMinCut), fJetRatioMaxCut(g.fJetRatioMaxCut),
  fJetRatioMinCut(g.fJetRatioMinCut),  fJetNCone(g.fJetNCone),
  fJetNPt(g.fJetNPt), fJetCone(g.fJetCone),
  fJetPtThreshold(g.fJetPtThreshold),fJetPtThresPbPb(g.fJetPtThresPbPb),
  fPtJetSelectionCut(g.fPtJetSelectionCut), fSelect(g.fSelect),  
  fhPhiCharged(g.fhPhiCharged), fhPhiNeutral(g.fhPhiNeutral), 
  fhEtaCharged(g.fhEtaCharged), fhEtaNeutral(g.fhEtaNeutral), 
  fhDeltaPhiGammaCharged(g.fhDeltaPhiGammaCharged),  
  fhDeltaPhiGammaNeutral(g.fhDeltaPhiGammaNeutral), 
  fhDeltaEtaGammaCharged(g.fhDeltaEtaGammaCharged), 
  fhDeltaEtaGammaNeutral(g.fhDeltaEtaGammaNeutral), 
  fhAnglePairLeading(g.fhAnglePairLeading), fhInvMassPairLeading(g.fhInvMassPairLeading), 
  fhChargedRatio(g.fhChargedRatio), fhNeutralRatio(g.fhNeutralRatio), 
  fhNBkg(g. fhNBkg), fhNLeading(g. fhNLeading), fhNJet(g.fhNJet), fhJetRatio(g.fhJetRatio), fhJetPt(g.fhJetPt), 
  fhBkgRatio (g.fhBkgRatio), fhBkgPt(g.fhBkgPt),  fhJetFragment(g.fhJetFragment), fhBkgFragment(g.fhBkgFragment), 
  fhJetPtDist(g.fhJetPtDist),  fhBkgPtDist(g.fhBkgPtDist)   
{
  // cpy ctor

  for(Int_t i = 0; i<10; i++){
    fJetCones[i]        = g.fJetCones[i] ;
    fJetNameCones[i]    = g.fJetNameCones[i] ;
    fJetPtThres[i]      = g.fJetPtThres[i] ;
    fJetNamePtThres[i]  = g.fJetNamePtThres[i] ;
    if( i < 6 ){
      fJetXMin1[i]       = g.fJetXMin1[i] ;
      fJetXMin2[i]       = g.fJetXMin2[i] ;
      fJetXMax1[i]       = g.fJetXMax1[i] ;
      fJetXMax2[i]       = g.fJetXMax2[i] ;
      fBkgMean[i]        = g.fBkgMean[i] ;
      fBkgRMS[i]         = g.fBkgRMS[i] ;
      if( i < 2 ){
	fJetE1[i]        = g.fJetE1[i] ;
	fJetE2[i]        = g.fJetE2[i] ;
	fJetSigma1[i]    = g.fJetSigma1[i] ;
	fJetSigma2[i]    = g.fJetSigma2[i] ;
      }
    }          
  } 

  
}

//_________________________________________________________________________
AliAnaParticleJetLeadingCone & AliAnaParticleJetLeadingCone::operator = (const AliAnaParticleJetLeadingCone & source)
{
  // assignment operator

  if(this == &source)return *this;
  ((AliAnaBaseClass *)this)->operator=(source);

  fSeveralConeAndPtCuts = source.fSeveralConeAndPtCuts ; 
  fPbPb = source.fPbPb ;
  fJetCTSRatioMaxCut = source.fJetCTSRatioMaxCut ;
  fJetCTSRatioMinCut = source.fJetCTSRatioMinCut ; fJetRatioMaxCut = source.fJetRatioMaxCut ;
  fJetRatioMinCut = source.fJetRatioMinCut ;  fJetNCone = source.fJetNCone ;
  fJetNPt = source.fJetNPt ; fJetCone = source.fJetCone ; 
  fJetPtThreshold = source.fJetPtThreshold ;
  fJetPtThresPbPb = source.fJetPtThresPbPb ;
  fPtJetSelectionCut = source.fPtJetSelectionCut ;
  fSelect = source.fSelect ;  fhChargedRatio = source.fhChargedRatio ; fhNeutralRatio = source.fhNeutralRatio ; 

  fhPhiCharged = source.fhPhiCharged ; fhPhiNeutral = source.fhPhiNeutral ; 
  fhEtaCharged = source.fhEtaCharged ; fhEtaNeutral = source.fhEtaNeutral ; 
  fhDeltaPhiGammaCharged = source.fhDeltaPhiGammaCharged ;  
  fhDeltaPhiGammaNeutral = source.fhDeltaPhiGammaNeutral ; 
  fhDeltaEtaGammaCharged = source.fhDeltaEtaGammaCharged ; 
  fhDeltaEtaGammaNeutral = source.fhDeltaEtaGammaNeutral ; 
  
  fhAnglePairLeading = source.fhAnglePairLeading ; 
  fhInvMassPairLeading = source.fhInvMassPairLeading ; 
  fhNBkg = source. fhNBkg ; fhNLeading = source. fhNLeading ; 
  fhNJet = source.fhNJet ; fhJetRatio = source.fhJetRatio ; fhJetPt = source.fhJetPt ; 
  fhBkgRatio  = source.fhBkgRatio ; fhBkgPt = source.fhBkgPt ;  
  fhJetFragment = source.fhJetFragment ; fhBkgFragment = source.fhBkgFragment ; 
  fhJetPtDist = source.fhJetPtDist ;  fhBkgPtDist = source.fhBkgPtDist ;


  for(Int_t i = 0; i<10; i++){
    fJetCones[i]        = source.fJetCones[i] ;
    fJetNameCones[i]    = source.fJetNameCones[i] ;
    fJetPtThres[i]      = source.fJetPtThres[i] ;
    fJetNamePtThres[i]  = source.fJetNamePtThres[i] ;
    if( i < 6 ){
      fJetXMin1[i]       = source.fJetXMin1[i] ;
      fJetXMin2[i]       = source.fJetXMin2[i] ;
      fJetXMax1[i]       = source.fJetXMax1[i] ;
      fJetXMax2[i]       = source.fJetXMax2[i] ;
      fBkgMean[i]        = source.fBkgMean[i] ;
      fBkgRMS[i]         = source.fBkgRMS[i] ;
      if( i < 2 ){
	fJetE1[i]        = source.fJetE1[i] ;
	fJetE2[i]        = source.fJetE2[i] ;
	fJetSigma1[i]    = source.fJetSigma1[i] ;
	fJetSigma2[i]    = source.fJetSigma2[i] ;

      }
    }          
  } 

  return *this;

}

//____________________________________________________________________________
AliAnaParticleJetLeadingCone::~AliAnaParticleJetLeadingCone() 
{
   // Remove all pointers except analysis output pointers.
  delete [] fJetE1;  
  delete [] fJetE2;    
  delete [] fJetSigma1;
  delete [] fJetSigma2;
  delete [] fBkgMean; 
  delete [] fBkgRMS;  
  delete [] fJetXMin1;
  delete [] fJetXMin2;
  delete [] fJetXMax1;
  delete [] fJetXMax2; 
  delete [] fJetCones;         
  delete [] fJetNameCones;   
  delete [] fJetPtThres;       
  delete [] fJetNamePtThres;  
}



//________________________________________________________________________
TList *  AliAnaParticleJetLeadingCone::GetCreateOutputObjects()
{  
  // Create histograms to be saved in output file and 
  // store them in fOutputContainer
 
  AliDebug(1,"Init jet in leading cone histograms");
  
  TList * outputContainer = new TList() ; 
  outputContainer->SetName("GammaJetCorrelationHistos") ; 

  fhChargedRatio  = new TH2F
    ("ChargedRatio","p_{T leading charge} /p_{T #gamma} vs p_{T #gamma}",
     120,0,120,120,0,1); 
  fhChargedRatio->SetYTitle("p_{T lead charge} /p_{T #gamma}");
  fhChargedRatio->SetXTitle("p_{T #gamma} (GeV/c)");
  
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
  
  fhDeltaEtaGammaCharged  = new TH2F
    ("DeltaEtaGammaCharged","#eta_{#gamma} - #eta_{#pi^{#pm}} vs p_{T #gamma}",
     200,0,120,200,-2,2); 
  fhDeltaEtaGammaCharged->SetYTitle("#Delta #eta");
  fhDeltaEtaGammaCharged->SetXTitle("p_{T #gamma} (GeV/c)");
  
  outputContainer->Add(fhPhiCharged) ;
  outputContainer->Add(fhEtaCharged) ;
  outputContainer->Add(fhChargedRatio) ;
  outputContainer->Add(fhDeltaPhiGammaCharged) ; 
  outputContainer->Add(fhDeltaEtaGammaCharged) ; 

  if(!AreJetsOnlyInCTS()){
    
    fhNeutralRatio  = new TH2F
      ("NeutralRatio","p_{T leading  #pi^{0}} /p_{T #gamma} vs p_{T #gamma}",
       120,0,120,120,0,1); 
    fhNeutralRatio->SetYTitle("p_{T lead  #pi^{0}} /p_{T #gamma}");
    fhNeutralRatio->SetXTitle("p_{T #gamma} (GeV/c)");
    
    //
    fhAnglePairLeading  = new TH2F
      ("AnglePairLeading",
       "Angle between all #gamma pair finally selected vs p_{T  #pi^{0}}",
       200,0,50,200,0,0.2); 
    fhAnglePairLeading->SetYTitle("Angle (rad)");
    fhAnglePairLeading->SetXTitle("E_{ #pi^{0}} (GeV/c)");
    
    fhInvMassPairLeading  = new TH2F
      ("InvMassPairLeading",
       "Invariant Mass of #gamma pair selected vs p_{T #gamma}",
       120,0,120,360,0,0.5); 
    fhInvMassPairLeading->SetYTitle("Invariant Mass (GeV/c^{2})");
    fhInvMassPairLeading->SetXTitle("p_{T #gamma} (GeV/c)");

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

    outputContainer->Add(fhPhiNeutral) ;  
    outputContainer->Add(fhEtaNeutral) ;  
    outputContainer->Add(fhNeutralRatio) ; 
    outputContainer->Add(fhDeltaPhiGammaNeutral) ; 
    outputContainer->Add(fhDeltaEtaGammaNeutral) ;
    
    outputContainer->Add(fhInvMassPairLeading) ; 
    outputContainer->Add(fhAnglePairLeading) ; 
  }
  
  if(!fSeveralConeAndPtCuts){// not several cones
    
    //Count
    fhNBkg = new TH1F("NBkg","bkg multiplicity",9000,0,9000); 
    fhNBkg->SetYTitle("counts");
    fhNBkg->SetXTitle("N");
    outputContainer->Add(fhNBkg) ; 
    
    fhNLeading  = new TH2F
      ("NLeading","Accepted Jet Leading", 240,0,120,240,0,120); 
    fhNLeading->SetYTitle("p_{T charge} (GeV/c)");
    fhNLeading->SetXTitle("p_{T #gamma}(GeV/c)");
    outputContainer->Add(fhNLeading) ; 
    
    fhNJet  = new TH1F("NJet","Accepted jets",240,0,120); 
    fhNJet->SetYTitle("N");
    fhNJet->SetXTitle("p_{T #gamma}(GeV/c)");
    outputContainer->Add(fhNJet) ; 
    
    //Ratios and Pt dist of reconstructed (not selected) jets
    //Jet
    fhJetRatio  = new TH2F
      ("JetRatio","p_{T jet lead}/p_{T #gamma} vs p_{T #gamma}",
       240,0,120,200,0,10);
    fhJetRatio->SetYTitle("p_{T jet lead #pi^{0}}/p_{T #gamma}");
    fhJetRatio->SetXTitle("p_{T #gamma} (GeV/c)");
    outputContainer->Add(fhJetRatio) ; 
    
    fhJetPt  = new TH2F
      ("JetPt", "p_{T jet lead} vs p_{T #gamma}",240,0,120,400,0,200);
    fhJetPt->SetYTitle("p_{T jet}");
    fhJetPt->SetXTitle("p_{T #gamma} (GeV/c)");
    outputContainer->Add(fhJetPt) ; 
    
    //Bkg
    
    fhBkgRatio  = new TH2F
      ("BkgRatio","p_{T bkg lead}/p_{T #gamma} vs p_{T #gamma}",
       240,0,120,200,0,10);
    fhBkgRatio->SetYTitle("p_{T bkg lead charge}/p_{T #gamma}");
    fhBkgRatio->SetXTitle("p_{T #gamma} (GeV/c)");
    outputContainer->Add(fhBkgRatio) ;
    
    fhBkgPt  = new TH2F
      ("BkgPt","p_{T jet lead} vs p_{T #gamma}",240,0,120,400,0,200);
    fhBkgPt->SetYTitle("p_{T jet lead charge}/p_{T #gamma}");
    fhBkgPt->SetXTitle("p_{T #gamma} (GeV/c)");
    outputContainer->Add(fhBkgPt) ;
    
    //Jet Distributions
    
    fhJetFragment  = 
      new TH2F("JetFragment","x = p_{T i charged}/p_{T #gamma}",
	       240,0.,120.,1000,0.,1.2); 
    fhJetFragment->SetYTitle("x_{T}");
    fhJetFragment->SetXTitle("p_{T #gamma}");
    outputContainer->Add(fhJetFragment) ;
    
    fhBkgFragment  = new TH2F
      ("BkgFragment","x = p_{T i charged}/p_{T #gamma}",
       240,0.,120.,1000,0.,1.2);
    fhBkgFragment->SetYTitle("x_{T}");
    fhBkgFragment->SetXTitle("p_{T #gamma}");
    outputContainer->Add(fhBkgFragment) ;
    
    fhJetPtDist  = 
      new TH2F("JetPtDist","x = p_{T i charged}",240,0.,120.,400,0.,200.); 
    fhJetPtDist->SetXTitle("p_{T #gamma} (GeV/c)");
    outputContainer->Add(fhJetPtDist) ;
    
    fhBkgPtDist  = new TH2F
      ("BkgPtDist","x = p_{T i charged}",240,0.,120.,400,0.,200.); 
    fhBkgPtDist->SetXTitle("p_{T #gamma} (GeV/c)");
    outputContainer->Add(fhBkgPtDist) ;
    
  }//not several cones
  else{ //If we want to study the jet for different cones and pt
    for(Int_t icone = 0; icone<fJetNCone; icone++){//icone
      for(Int_t ipt = 0; ipt<fJetNPt;ipt++){ //ipt
	
	//Jet
	
	fhJetRatios[icone][ipt]  = new TH2F
	  ("JetRatioCone"+fJetNameCones[icone]+"Pt"+fJetNamePtThres[ipt], 
	   "p_{T jet lead #pi^{0}}/p_{T #gamma} vs p_{T #gamma}, cone ="
	   +fJetNameCones[icone]+", pt>" +fJetNamePtThres[ipt]+" GeV/c",
	   240,0,120,200,0,10);
	fhJetRatios[icone][ipt]->
	  SetYTitle("p_{T jet lead #pi^{0}}/p_{T #gamma}");
	fhJetRatios[icone][ipt]->SetXTitle("p_{T #gamma} (GeV/c)");
	outputContainer->Add(fhJetRatios[icone][ipt]) ; 
	
	
	fhJetPts[icone][ipt]  = new TH2F
	  ("JetPtCone"+fJetNameCones[icone]+"Pt"+fJetNamePtThres[ipt], 
	   "p_{T jet lead #pi^{0}}/p_{T #gamma} vs p_{T #gamma}, cone ="
	   +fJetNameCones[icone]+", pt>" +fJetNamePtThres[ipt]+" GeV/c",
	   240,0,120,400,0,200);
	fhJetPts[icone][ipt]->
	  SetYTitle("p_{T jet lead #pi^{0}}/p_{T #gamma}");
	fhJetPts[icone][ipt]->SetXTitle("p_{T #gamma} (GeV/c)");
	outputContainer->Add(fhJetPts[icone][ipt]) ; 
	
	//Bkg
	fhBkgRatios[icone][ipt]  = new TH2F
	  ("BkgRatioCone"+fJetNameCones[icone]+"Pt"+fJetNamePtThres[ipt], 
	   "p_{T bkg lead #pi^{0}}/p_{T #gamma} vs p_{T #gamma}, cone ="
	   +fJetNameCones[icone]+", pt>" +fJetNamePtThres[ipt]+" GeV/c",
	   240,0,120,200,0,10);
	fhBkgRatios[icone][ipt]->
	  SetYTitle("p_{T bkg lead #pi^{0}}/p_{T #gamma}");
	fhBkgRatios[icone][ipt]->SetXTitle("p_{T #gamma} (GeV/c)");
	outputContainer->Add(fhBkgRatios[icone][ipt]) ; 
	
	fhBkgPts[icone][ipt]  = new TH2F
	  ("BkgPtCone"+fJetNameCones[icone]+"Pt"+fJetNamePtThres[ipt], 
	   "p_{T jet lead #pi^{0}}/p_{T #gamma} vs p_{T #gamma}, cone ="
	   +fJetNameCones[icone]+", pt>" +fJetNamePtThres[ipt]+" GeV/c",
	   240,0,120,400,0,200);
	fhBkgPts[icone][ipt]->
	  SetYTitle("p_{T jet lead #pi^{0}}/p_{T #gamma}");
	fhBkgPts[icone][ipt]->SetXTitle("p_{T #gamma} (GeV/c)");
	outputContainer->Add(fhBkgPts[icone][ipt]) ; 
	
	//Counts
	fhNBkgs[icone][ipt]  = new TH1F
	  ("NBkgCone"+fJetNameCones[icone]+"Pt"+fJetNamePtThres[ipt],
	   "bkg multiplicity cone ="+fJetNameCones[icone]+", pt>" 
	   +fJetNamePtThres[ipt]+" GeV/c",9000,0,9000); 
	fhNBkgs[icone][ipt]->SetYTitle("counts");
	fhNBkgs[icone][ipt]->SetXTitle("N");
	outputContainer->Add(fhNBkgs[icone][ipt]) ; 
	
	fhNLeadings[icone][ipt]  = new TH2F
	  ("NLeadingCone"+fJetNameCones[icone]+"Pt"+fJetNamePtThres[ipt],
	   "p_{T #gamma} vs p_{T #pi^{0}} cone ="+fJetNameCones[icone]+", pt>" 
	   +fJetNamePtThres[ipt]+" GeV/c",120,0,120,120,0,120); 
	fhNLeadings[icone][ipt]->SetYTitle("p_{T #pi^{0}}(GeV/c)");
	fhNLeadings[icone][ipt]->SetXTitle("p_{T #gamma}(GeV/c)");
	outputContainer->Add(fhNLeadings[icone][ipt]) ; 
	
	fhNJets[icone][ipt]  = new TH1F
	  ("NJetCone"+fJetNameCones[icone]+"Pt"+fJetNamePtThres[ipt],
	   "Number of neutral jets, cone ="+fJetNameCones[icone]+", pt>" 
	   +fJetNamePtThres[ipt]+" GeV/c",120,0,120); 
	fhNJets[icone][ipt]->SetYTitle("N");
	fhNJets[icone][ipt]->SetXTitle("p_{T #gamma}(GeV/c)");
	outputContainer->Add(fhNJets[icone][ipt]) ; 
	
	//Fragmentation Function
	fhJetFragments[icone][ipt]  = new TH2F
	  ("JetFragmentCone"+fJetNameCones[icone]+"Pt"+fJetNamePtThres[ipt],
	   "x_{T} = p_{T i}/p_{T #gamma}, cone ="+fJetNameCones[icone]+", pt>" 
	   +fJetNamePtThres[ipt]+" GeV/c",120,0.,120.,240,0.,1.2); 
	fhJetFragments[icone][ipt]->SetYTitle("x_{T}");
	fhJetFragments[icone][ipt]->SetXTitle("p_{T #gamma}");
	outputContainer->Add(fhJetFragments[icone][ipt]) ; 
	
	fhBkgFragments[icone][ipt]  = new TH2F
	  ("BkgFragmentCone"+fJetNameCones[icone]+"Pt"+fJetNamePtThres[ipt],
	   "x_{T} = p_{T i}/p_{T #gamma}, cone ="+fJetNameCones[icone]+", pt>" 
	   +fJetNamePtThres[ipt]+" GeV/c",120,0.,120.,240,0.,1.2); 
	fhBkgFragments[icone][ipt]->SetYTitle("x_{T}");
	fhBkgFragments[icone][ipt]->SetXTitle("p_{T #gamma}");
	outputContainer->Add(fhBkgFragments[icone][ipt]) ; 
	
	//Jet particle distribution
	
	fhJetPtDists[icone][ipt]  = new TH2F
	  ("JetPtDistCone"+fJetNameCones[icone]+"Pt"+fJetNamePtThres[ipt],
	   "p_{T i}, cone ="+fJetNameCones[icone]+", pt>" +fJetNamePtThres[ipt]+
	   " GeV/c",120,0.,120.,120,0.,120.); 
	fhJetPtDists[icone][ipt]->SetXTitle("p_{T #gamma} (GeV/c)");
	outputContainer->Add(fhJetPtDists[icone][ipt]) ; 
	
	fhBkgPtDists[icone][ipt]  = new TH2F
	  ("BkgPtDistCone"+fJetNameCones[icone]+"Pt"+fJetNamePtThres[ipt],
	   "p_{T i}, cone ="+fJetNameCones[icone]+", pt>" +fJetNamePtThres[ipt]+
	   " GeV/c",120,0.,120.,120,0.,120.); 
	fhBkgPtDists[icone][ipt]->SetXTitle("p_{T #gamma} (GeV/c)");
	outputContainer->Add(fhBkgPtDists[icone][ipt]) ; 
	
      }//ipt
    } //icone
  }//If we want to study any cone or pt threshold

  return outputContainer;
}

  //____________________________________________________________________________
  void AliAnaParticleJetLeadingCone::InitParameters()
{
 
  //Initialize the parameters of the analysis.
  SetJetsOnlyInCTS(kFALSE) ;
  fPbPb                = kFALSE ;

  SetDeltaPhiCutRange(2.9,3.4) ; 
  SetRatioCutRange(0.1,1.0) ; 

  //Jet selection parameters
  //Fixed cut   
  fJetRatioMaxCut = 1.2 ; 
  fJetRatioMinCut = 0.3 ; 
  fJetCTSRatioMaxCut = 1.2 ;
  fJetCTSRatioMinCut = 0.3 ;
  fSelect         = 0  ;

  //Cut depending on gamma energy

  fPtJetSelectionCut = 20.; //For Low pt jets+BKG, another limits applied
  //Reconstructed jet energy dependence parameters 
  //e_jet = a1+e_gamma b2. 
  //Index 0-> Pt>2 GeV r = 0.3; Index 1-> Pt>0.5 GeV r = 0.3
  fJetE1[0] = -5.75; fJetE1[1] = -4.1;
  fJetE2[0] = 1.005; fJetE2[1] = 1.05;

  //Reconstructed sigma of jet energy dependence parameters 
  //s_jet = a1+e_gamma b2. 
  //Index 0-> Pt>2 GeV r = 0.3; Index 1-> Pt>0.5 GeV r = 0.3
  fJetSigma1[0] = 2.65;   fJetSigma1[1] = 2.75;
  fJetSigma2[0] = 0.0018; fJetSigma2[1] = 0.033;

  //Background mean energy and RMS
  //Index 0-> No BKG; Index 1-> BKG > 2 GeV; 
  //Index 2-> (low pt jets)BKG > 0.5 GeV;
  //Index > 2, same for CTS conf
  fBkgMean[0] = 0.; fBkgMean[1] = 8.8 ; fBkgMean[2] = 69.5;
  fBkgMean[3] = 0.; fBkgMean[4] = 6.4;  fBkgMean[5] = 48.6;
  fBkgRMS[0]  = 0.; fBkgRMS[1]  = 7.5;  fBkgRMS[2]  = 22.0; 
  fBkgRMS[3]  = 0.; fBkgRMS[4]  = 5.4;  fBkgRMS[5]  = 13.2; 

  //Factor x of min/max = E -+ x * sigma. Obtained after selecting the
  //limits for monoenergetic jets.
  //Index 0-> No BKG; Index 1-> BKG > 2 GeV; 
  //Index 2-> (low pt jets) BKG > 0.5 GeV;
  //Index > 2, same for CTS conf

  fJetXMin1[0] =-0.69 ; fJetXMin1[1] = 0.39 ; fJetXMin1[2] =-0.88 ; 
  fJetXMin1[3] =-2.0  ; fJetXMin1[4] =-0.442 ; fJetXMin1[5] =-1.1  ;
  fJetXMin2[0] = 0.066; fJetXMin2[1] = 0.038; fJetXMin2[2] = 0.034; 
  fJetXMin2[3] = 0.25 ; fJetXMin2[4] = 0.113; fJetXMin2[5] = 0.077 ;
  fJetXMax1[0] =-3.8  ; fJetXMax1[1] =-0.76 ; fJetXMax1[2] =-3.6  ; 
  fJetXMax1[3] =-2.7  ; fJetXMax1[4] =-1.21 ; fJetXMax1[5] =-3.7  ;
  fJetXMax2[0] =-0.012; fJetXMax2[1] =-0.022; fJetXMax2[2] = 0.016; 
  fJetXMax2[3] =-0.024; fJetXMax2[4] =-0.008; fJetXMax2[5] = 0.027;


  //Different cones and pt thresholds to construct the jet

  fJetCone        = 0.3  ;
  fJetPtThreshold = 0.5   ;
  fJetPtThresPbPb = 2.   ;
  fJetNCone       = 4    ;
  fJetNPt         = 4    ;
  fJetCones[0]    = 0.2  ; fJetNameCones[0]   = "02" ;
  fJetCones[1]    = 0.3  ; fJetNameCones[1]   = "03" ;
  fJetCones[2]    = 0.4  ; fJetNameCones[2]   = "04" ;
  fJetCones[2]    = 0.5  ; fJetNameCones[2]   = "05" ;

  fJetPtThres[0]  = 0.0  ; fJetNamePtThres[0] = "00" ;
  fJetPtThres[1]  = 0.5  ; fJetNamePtThres[1] = "05" ;
  fJetPtThres[2]  = 1.0  ; fJetNamePtThres[2] = "10" ;
  fJetPtThres[3]  = 2.0  ; fJetNamePtThres[3] = "20" ;
}

//__________________________________________________________________
void AliAnaParticleJetLeadingCone::Print(const Option_t * opt) const
{

  //Print some relevant parameters set for the analysis
  if(! opt)
    return;
  
  Info("Print", "%s %s", GetName(), GetTitle() ) ;
  printf("Correlation analysis           =     %d\n",kJetLeadCone) ;
  
  
  printf("Phi gamma-Leading      <     %f\n", GetDeltaPhiMaxCut()) ; 
  printf("Phi gamma-Leading      >     %f\n", GetDeltaPhiMinCut()) ;
  printf("pT Leading / pT Gamma             <     %f\n", GetRatioMaxCut()) ; 
  printf("pT Leading / pT Gamma             >     %f\n", GetRatioMinCut()) ;
  
  if(fSelect == 2){
    printf("pT Jet / pT Gamma                     <    %f\n", fJetRatioMaxCut) ; 
    printf("pT Jet / pT Gamma                     >    %f\n", fJetRatioMinCut) ;
    printf("pT Jet (Only CTS)/ pT Gamma   <    %f\n", fJetCTSRatioMaxCut) ; 
    printf("pT Jet (Only CTS)/ pT Gamma   >    %f\n", fJetCTSRatioMinCut) ;
  }
  
} 

//__________________________________________________________________
void  AliAnaParticleJetLeadingCone::MakeAnalysisFillAOD() 
{
  
  //Particle-Hadron Correlation Analysis, fill AODs
  if(GetDebug() > 1){
    printf("Begin jet leading cone  correlation analysis, fill AODs \n");
    printf("In particle branch aod entries %d\n", GetAODBranch()->GetEntries());
    printf("In CTS aod entries %d\n", GetAODCTS()->GetEntries());
    printf("In EMCAL aod entries %d\n", GetAODEMCAL()->GetEntries());
  }
 
  TLorentzVector * pLeading = new TLorentVector; //It will contain the kinematics of the found leading particle
  
  //Loop on stored AOD particles, trigger
  Int_t naod = GetAODBranch()->GetEntriesFast();
  for(Int_t iaod = 0; iaod < naod ; iaod++){
    AliAODParticleCorrelation* particle =  dynamic_cast<AliAODParticleCorrelation*> (GetAODBranch()->At(iaod));
    
    //Search leading particles in CTS and EMCAL 
    if(GetLeadingParticle(particle, pLeading)){
      
      if(GetDebug() > 1) printf("Leading: pt %f, phi %f, eta %f", pLeading->Pt(),pLeading->Phi(),pLeading->Eta())) ;
      
    MakeJet(particle, pLeading,"", kFALSE);

    }//Leading
  }//AOD trigger particle loop
  
  if(GetDebug() >1)printf("End of jet leading cone analysis, fill AODs \n");

} 

//__________________________________________________________________
void  AliAnaParticleJetLeadingCone::MakeAnalysisFillHistograms() 
{
  
  //Particle-Hadron Correlation Analysis, fill AODs
  if(GetDebug() > 1){
    printf("Begin jet leading cone  correlation analysis, fill histograms \n");
    printf("In particle branch aod entries %d\n", GetAODBranch()->GetEntries());
    printf("In CTS aod entries %d\n", GetAODCTS()->GetEntries());
    printf("In EMCAL aod entries %d\n", GetAODEMCAL()->GetEntries());
  }
 
  TLorentzVector * pLeading = new TLorentzVector;
  
  //Loop on stored AOD particles, trigger
  Int_t naod = GetAODBranch()->GetEntriesFast();
  for(Int_t iaod = 0; iaod < naod ; iaod++){
    AliAODParticleCorrelation* particle =  dynamic_cast<AliAODParticleCorrelation*> (GetAODBranch()->At(iaod));
    
    Double_t pt = particle->Pt();
    Double_t phi = particle->Phi();
    Double_t eta = particle->Eta();
    
    //Get leading particle, fill histograms
    pLeading = particle->GetLeading();
    TString det = particle->GetLeadingDetector();      

    if(det!="" && pLeading){
      Double_t ptL = pLeading->Pt(); 
      Double_t phiL = pLeading->Phi(); 
      if(phiL < 0 ) phiL+=TMath::TwoPi();
      Double_t etaL = pLeading->Eta(); 
      
      if(GetDebug() > 1) printf("Leading found in %s, with pt %3.2f, phi %2.2f, eta %2.2f",det.Data(), ptL, phiL, etaL);
      if(det == "CTS"){
	fhChargedRatio->Fill(pt,ptL/particle->Pt());
	fhEtaCharged->Fill(pt,etaL);
	fhPhiCharged->Fill(pt,phiL);
	fhDeltaPhiGammaCharged->Fill(pt,phi-phiL);
	fhDeltaEtaGammaCharged->Fill(pt,eta-etaL);
      }
      else if(det== "EMCAL"){
	fhNeutralRatio->Fill(pt,ptL/particle->Pt());
	fhEtaNeutral->Fill(pt,etaL);
	fhPhiNeutral->Fill(pt,phiL);
	fhDeltaPhiGammaNeutral->Fill(pt,phi-phiL);
	fhDeltaEtaGammaNeutral->Fill(pt,eta-etaL);
      }
    }
      
      //Search Jet
      if(!fSeveralConeAndPtCuts)
	MakeJet(particle, pLeading, "", kTRUE);
      else{
	for(Int_t icone = 0; icone<fJetNCone; icone++) {
	  for(Int_t ipt = 0; ipt<fJetNPt;ipt++) {  
	    TString lastname ="Cone"+ fJetNameCones[icone]+"Pt"+ fJetNamePtThres[ipt];
	    fJetCone=fJetCones[icone];
	    fJetPtThreshold=fJetPtThres[ipt];
	    MakeJet(particle, pLeading, lastname);
	  }//icone
	}//ipt
      }//fSeveralConeAndPtCuts
    }//Leading
  }//AOD trigger particle loop
  
  if(GetDebug() >1)printf("End of jet leading cone analysis, fill histograms \n");

} 

//____________________________________________________________________________
Bool_t  AliAnaParticleJetLeadingCone::GetLeadingParticle(AliAODParticleCorrelation *particle, TLorentzVector * pLeading) 
{
  //Search Charged or Neutral leading particle, select the highest one.
  
  TLorentzVector *pLeadingCh = new TLorentzVector;
  TLorentzVector *pLeadingPi0 = new TLorentzVector;
  
  Double_t pt  =  particle->Pt(); 
  Double_t phi = particle->Phi(); 
  Double_t eta = particle->Eta(); 

  GetLeadingCharge(particle, pLeadingCh) ;
  if(!AreJetsOnlyInCTS()) GetLeadingPi0(particle, pLeadingPi0) ;
      
  Double_t ptch = pLeadingCh->Pt(); 
  Double_t phich = pLeadingCh->Phi(); 
  if(phich < 0 ) phich+=TMath::TwoPi();
  Double_t etach = pLeadingCh->Eta(); 
  Double_t ptpi = pLeadingPi0->Pt(); 
  Double_t phipi = pLeadingPi0->Phi(); 
  if(phipi < 0 ) phipi+=TMath::TwoPi();
  Double_t etapi = pLeadingPi0->Eta(); 
      
  if (ptch > 0 || ptpi > 0){
    if((ptch >= ptpi)){
      if(GetDebug() > 1)printf("Leading found in CTS \n");
      pLeading = pLeadingCh;
      particle->SetLeading(pLeading);
      return kTRUE;
    }
    else{
      if(GetDebug() > 1)printf("Leading found in EMCAL \n");
      pLeading = pLeadingPi0;
      particle->SetLeading(pLeading);	 
      return kTRUE;   
    }
    
    if(GetDebug() > 1)printf ("NO LEADING PARTICLE FOUND \n");
    return kFALSE; 
  
}

//____________________________________________________________________________
void  AliAnaParticleJetLeadingCone::GetLeadingCharge(AliAODParticleCorrelation * particle, TLorentzVector * pLeading) 
{  
  //Search for the charged particle with highest pt and with 
  //Phi=Phi_trigger-Pi and pT=0.1E_gamma 

  if(GetAODCTS()){
    Double_t ptTrig = particle->Pt();
    Double_t phiTrig = particle->Phi();
    Double_t rat  = -100 ;
    Double_t ptl  = -100 ;
    Double_t phil = -100 ;
    Double_t pt  = -100.;
    Double_t phi = -100.;
    TVector3 p3;
    
    for(Int_t ipr = 0;ipr < GetAODCTS()->GetEntries() ; ipr ++ ){
      AliAODTrack* track = dynamic_cast<AliAODTrack *>(GetAODCTS()->At(ipr)) ;
      p3.SetXYZ(track->Px(),track->Py(),track->Pz());
      pt    = p3.Pt();
      phi  = p3.Phi() ;
      if(phi<0) phi+=TMath::TwoPi();
      rat  = ptl/particle->Pt() ;
      
      //Selection within angular and energy limits
      if(((phiTrig-phi) > GetDeltaPhiMinCut()) && ((phiTrig-phi)<GetDeltaPhiMaxCut()) &&
	 (rat > GetRatioMinCut()) && (rat < GetRatioMaxCut())  && (pt  > ptl)) {
	phil = phi ;
	ptl  = pt ;
	pLeading->SetMomentum(p3.Px(),p3.Py(),p3.Pz(),0);
      }
    }// track loop
  }//CTS list exist
  
  if(GetDebug() > 1)printf("Leading in CTS: pt %f eta %f phi %f pt/ptTrig %f \n", ptl, pLeading->Eta(), phil,ptl/ptTrig)) ;
 
}

//____________________________________________________________________________
void  AliAnaParticleJetLeadingCone::GetLeadingPi0(AliAODParticleCorrelation * particle, TLorentzVector * pLeading) 
{  
  //Search for the neutral pion with highest pt and with 
  //Phi=Phi_trigger-Pi and pT=0.1E_gamma 
  if(GetAODEMCAL()){
    Double_t ptTrig = particle->Pt();
    Double_t phiTrig = particle->Phi();
    Double_t rat  = -100 ;
    Double_t ptl  = -100 ;
    Double_t phil = -100 ;
    Double_t pt  = -100.;
    Double_t phi = -100.;
    
    TLorentzVector gammai;
    TLorentzVector gammaj;
    
    Double_t vertex[] = {0,0,0};
    
    if(!GetReader()->GetDataType()== AliCaloTrackReader::kMC) GetReader()->GetVertex(vertex);
    
    //Cluster loop, select pairs with good pt, phi and fill AODs or histograms
    for(Int_t iclus = 0;iclus < GetAODEMCAL()->GetEntries() ; iclus ++ ){
      AliAODCaloCluster * calo = dynamic_cast< AliAODCaloCluster *>(GetAODEMCAL()->At(iclus)) ;
      
      //Cluster selection, not charged, with photon or pi0 id and in fidutial cut
      Int_t pdg=0;
      if(!SelectCluster(calo2,vertex, gammaj, pdgj)) continue ;
      
      if(GetDebug() > 2) printf("neutral cluster: pt %f, phi %f \n", gammai.Pt(),gammai.Phi());
	
      //2 gamma overlapped, found with PID
      if(pdg == AliCaloPID::kPi0){ 
	pt = gammai->Pt();
	rat = pt/ptTrig;
	phi = gammai->Phi();
	if(phi<0) phi+=TMath::TwoPi();
	
	//Selection within angular and energy limits
	if(ptl > pt  && rat > GetRatioMinCut()  && rat < GetRatioMaxCut()  && 
	   (phig-phil) > GetDeltaPhiMinCut() && (phig-phil) < GetDeltaPhiMaxCut() )
	  {
	    phi = phil ;
	    pt  = ptl ;
	    pLeading->SetMomentum(gammai.Px(),gammai.Py(),gammai.Pz(),gammai.E());
	  }// cuts
      }// pdg = 111 
      //Make invariant mass analysis
      else if(pdg == AliCaloPID::kPhoton){	
	//Search the photon companion in case it comes from  a Pi0 decay
	//Apply several cuts to select the good pair;
	for(Int_t jclus = iclus+1; jclus < GetAODEMCAL()->GetEntries() ; jclus ++ ){
	  AliAODCaloCluster * calo2 = dynamic_cast< AliAODCaloCluster *>(GetAODEMCAL()->At(jclus)) ;
	  
	  //Cluster selection, not charged with photon or pi0 id and in fidutial cut
	  Int_t pdgj=0;
	  if(!SelectCluster(calo2,vertex, gammaj, pdgj)) continue ;
	  
	  if(pdgj == AliCaloPID::kPhoton ){
	    
	    pt  = (gammai+gammaj).Pt();
	    phi = (gammai+gammaj).Phi();
	    rat = pt/ptTrig;
	    
	    //Selection within angular and energy limits
	    if(ptl > pt  && rat > GetRatioMinCut()  && rat < GetRatioMaxCut()  && 
	       (phig-phil) > GetDeltaPhiMinCut() && (phig-phil) < GetDeltaPhiMaxCut() ){
	      //Select good pair (aperture and invariant mass)
 	      if(GetNeutralMesonSelection()->SelectPair(gammai, gammaj)){
		phi = phil ;
 		pt  = ptl ;
 		pLeading->SetMomentum(gammai.Px(),gammai.Py(),gammai.Pz(),gammai.E());
 	      }//pi0 selection
	      
	      if(GetDebug() > 3 ) printf("Neutral Hadron Correlation: Selected gamma pair: pt %2.2f, phi %2.2f, eta %2.2f, M %2.3f",
					 (gammai+gammaj).Pt(),(gammai+gammaj).Phi(),(gammai+gammaj).Eta(), (gammai+gammaj).M());
	    }//Pair selected as leading
	  }//if pair of gammas
	}//2nd loop
      }// if pdg = 22
    }// 1st Loop
  }//EMCAL list exists

  if(GetDebug()>2) printf("Leading EMCAL: pt %f eta %f phi %f pt/Eg %f \n",  pLeading->Pt(), pLeading->Eta(),  pLeading->Phi(),  pLeading->Pt()/ptTrig) ;
  
}


//__________________________________________________________________________-
Bool_t AliAnaParticleJetLeadingCone::IsJetSelected(const Double_t ptg, const Double_t ptj){
  //Check if the energy of the reconstructed jet is within an energy window

  Double_t par[6];
  Double_t xmax[2];
  Double_t xmin[2];

  Int_t iCTS = 0;
  if(AreJetsOnlyInCTS())
    iCTS = 3 ;

  if(!fPbPb){
    //Phythia alone, jets with pt_th > 0.2, r = 0.3 
    par[0] = fJetE1[0]; par[1] = fJetE2[0]; 
    //Energy of the jet peak
    //e_jet = fJetE1[0]+fJetE2[0]*e_gamma, simulation fit
    par[2] = fJetSigma1[0]; par[3] = fJetSigma2[0];
    //Sigma  of the jet peak
    //sigma_jet = fJetSigma1[0]+fJetSigma2[0]*e_gamma, simulation fit
    par[4] = fBkgMean[0 + iCTS]; par[5] = fBkgRMS[0 + iCTS];
    //Parameters reserved for PbPb bkg.
    xmax[0] = fJetXMax1[0 + iCTS]; xmax[1] = fJetXMax2[0 + iCTS];
    xmin[0] = fJetXMin1[0 + iCTS]; xmin[1] = fJetXMin2[0 + iCTS];
    //Factor that multiplies sigma to obtain the best limits, 
    //by observation, of mono jet ratios (ptjet/ptg)
    //X_jet = fJetX1[0]+fJetX2[0]*e_gamma
   
  }
  else{
    if(ptg > fPtJetSelectionCut){
      //Phythia +PbPb with  pt_th > 2 GeV/c, r = 0.3 
      par[0] = fJetE1[0]; par[1] = fJetE2[0]; 
      //Energy of the jet peak, same as in pp
      //e_jet = fJetE1[0]+fJetE2[0]*e_gamma, simulation fit
      par[2] = fJetSigma1[0]; par[3] = fJetSigma2[0];
      //Sigma  of the jet peak, same as in pp
      //sigma_jet = fJetSigma1[0]+fJetSigma2[0]*e_gamma, simulation fit
      par[4] = fBkgMean[1 + iCTS]; par[5] = fBkgRMS[1 + iCTS];
      //Mean value and RMS of PbPb Bkg 
      xmax[0] = fJetXMax1[1 + iCTS]; xmax[1] = fJetXMax2[1 + iCTS];
      xmin[0] = fJetXMin1[1 + iCTS]; xmin[1] = fJetXMin2[1 + iCTS];
      //Factor that multiplies sigma to obtain the best limits, 
      //by observation, of mono jet ratios (ptjet/ptg) mixed with PbPb Bkg, 
      //pt_th > 2 GeV, r = 0.3
      //X_jet = fJetX1[0]+fJetX2[0]*e_gamma
     
    }
    else{
      //Phythia + PbPb with  pt_th > 0.5 GeV/c, r = 0.3
      par[0] = fJetE1[1]; par[1] = fJetE2[1]; 
      //Energy of the jet peak, pt_th > 2 GeV/c, r = 0.3 
      //e_jet = fJetE1[0]+fJetE2[0]*e_gamma, simulation fit
      par[2] = fJetSigma1[1]; par[3] = fJetSigma2[1];
      //Sigma  of the jet peak, pt_th > 2 GeV/c, r = 0.3
      //sigma_jet = fJetSigma1[0]+fJetSigma2[0]*e_gamma, simulation fit
      par[4] = fBkgMean[2 + iCTS]; par[5] = fBkgRMS[2 + iCTS];
      //Mean value and RMS of PbPb Bkg in a 0.3 cone, pt > 2 GeV.
      xmax[0] = fJetXMax1[2 + iCTS]; xmax[1] = fJetXMax2[2 + iCTS];
      xmin[0] = fJetXMin1[2 + iCTS]; xmin[1] = fJetXMin2[2 + iCTS];
      //Factor that multiplies sigma to obtain the best limits, 
      //by observation, of mono jet ratios (ptjet/ptg) mixed with PbPb Bkg, 
      //pt_th > 2 GeV, r = 0.3
      //X_jet = fJetX1[0]+fJetX2[0]*e_gamma
     
    }//If low pt jet in bkg
  }//if Bkg

 //Calculate minimum and maximum limits of the jet ratio.
  Double_t min = CalculateJetRatioLimit(ptg, par, xmin);
  Double_t max = CalculateJetRatioLimit(ptg, par, xmax);
  
  AliDebug(3,Form("Jet selection?  : Limits min %f, max %f,  pt_jet %f,  pt_gamma %f, pt_jet / pt_gamma %f",min,max,ptj,ptg,ptj/ptg));

  if(( min < ptj/ptg ) && ( max > ptj/ptg))
    return kTRUE;
  else
    return kFALSE;

}

//____________________________________________________________________________
void AliAnaParticleJetLeadingCone::MakeJet(AliAODParticleCorrelation *particle, TParticle* pLeading,TString lastname)
{
  //Fill the jet with the particles around the leading particle with 
  //R=fJetCone and pt_th = fJetPtThres. Calculate the energy of the jet and 
  //check if we select it. Fill jet histograms
  
  TClonesArray * jetList = new TClonesArray("TParticle",1000);
  TClonesArray * bkgList = new TClonesArray("TParticle",1000);

  TLorentzVector jet   (0,0,0,0);  
  TLorentzVector bkg(0,0,0,0);
  TLorentzVector lv (0,0,0,0);

  Double_t ptjet = 0.0;
  Double_t ptbkg = 0.0;
  Int_t n0 = 0;
  Int_t n1 = 0;  
  Bool_t b1 = kFALSE;
  Bool_t b0 = kFALSE;
  
  Double_t ptg  = pGamma->Pt();
  Double_t phig = pGamma->Phi();
  Double_t ptl  = pLeading->Pt();
  Double_t phil = pLeading->Phi();
  Double_t etal = pLeading->Eta();

  Float_t ptcut = fJetPtThreshold;
  if(fPbPb && !fSeveralConeAndPtCuts && ptg > fPtJetSelectionCut)  ptcut = fJetPtThresPbPb ;
 
  //Add charged particles to jet
  TIter nextch(plCTS) ; 
  TParticle * particle = 0 ; 
  while ( (particle = dynamic_cast<TParticle*>(nextch())) ) {
    
    b0 = kFALSE;
    b1 = kFALSE;

    //Particles in jet 
    SetJet(particle, b0, fJetCone, etal, phil) ;  

    if(b0){
      new((*jetList)[n0++]) TParticle(*particle) ;
      particle->Momentum(lv);
      if(particle->Pt() > ptcut ){
	jet+=lv;
	ptjet+=particle->Pt();
      }
    }

    //Background around (phi_gamma-pi, eta_leading)
    SetJet(particle, b1, fJetCone,etal, phig) ;

    if(b1) { 
      new((*bkgList)[n1++]) TParticle(*particle) ;
      particle->Momentum(lv);
      if(particle->Pt() > ptcut ){
	bkg+=lv;
	ptbkg+=particle->Pt();    
      }  
    }
  }

   //Add neutral particles to jet
  TIter nextne(plNe) ; 
  particle = 0 ; 
  while ( (particle = dynamic_cast<TParticle*>(nextne())) ) {
    
    b0 = kFALSE;
    b1 = kFALSE;

    //Particles in jet 
    SetJet(particle, b0, fJetCone, etal, phil) ;  

    if(b0){
      new((*jetList)[n0++]) TParticle(*particle) ;
      particle->Momentum(lv);
      if(particle->Pt() > ptcut ){
	jet+=lv;
	ptjet+=particle->Pt();
      }
    }

    //Background around (phi_gamma-pi, eta_leading)
    SetJet(particle, b1, fJetCone,etal, phig) ;

    if(b1) { 
      new((*bkgList)[n1++]) TParticle(*particle) ;
      particle->Momentum(lv);
      if(particle->Pt() > ptcut ){
	bkg+=lv;
	ptbkg+=particle->Pt();    
      }  
    }
  }
  
  ptjet = jet.Pt();
  ptbkg = bkg.Pt();

  if(ptjet > 0.) {

    AliDebug(2,Form("Gamma   pt %f, Jet pt %f, Bkg pt %f",ptg,ptjet,ptbkg));
    
    //Fill histograms
    
    Double_t ratjet   = ptjet/ptg ;
    Double_t ratbkg  = ptbkg/ptg ;

    dynamic_cast<TH2F*>
      (GetOutputContainer()->FindObject("JetRatio"+lastname))
      ->Fill(ptg,ratjet);	 
    dynamic_cast<TH2F*>
      (GetOutputContainer()->FindObject("JetPt"+lastname))
      ->Fill(ptg,ptjet);
    
    dynamic_cast<TH2F*>
      (GetOutputContainer()->FindObject("BkgRatio"+lastname))
      ->Fill(ptg,ratbkg);
    
    dynamic_cast<TH2F*>
      (GetOutputContainer()->FindObject("BkgPt"+lastname))
      ->Fill(ptg,ptbkg);

    //Jet selection
    Bool_t kSelect = kFALSE;
    if(fSelect == 0)
      kSelect = kTRUE; //Accept all jets, no restriction
    else if(fSelect == 1){
      //Selection with parametrized cuts
      if(IsJetSelected(ptg,ptjet))   kSelect = kTRUE;
    }
    else if(fSelect == 2){
      //Simple selection
      if(!AreJetsOnlyInCTS()){
	if((ratjet <  fJetRatioMaxCut) && (ratjet > fJetRatioMinCut )) kSelect = kTRUE;
      }
      else{
	if((ratjet <  fJetCTSRatioMaxCut) && (ratjet > fJetCTSRatioMinCut )) kSelect = kTRUE;
      }
    }
    else
      AliError("Jet selection option larger than 2, DON'T SELECT JETS");
    
    
     if(kSelect){
       AliDebug(1,Form("Jet Selected: pt %f ", ptjet)) ;
      
       FillJetHistos(jetList, ptg, ptl,"Jet",lastname);
       FillJetHistos(bkgList, ptg, ptl, "Bkg",lastname);
     }
  } //ptjet > 0
  
  jetList ->Delete();
  bkgList ->Delete();
  
}

//___________________________________________________________________
void AliAnaParticleJetLeadingCone::SetJet(TParticle * part, Bool_t & b, Float_t cone, 
			     Double_t eta, Double_t phi)
{

  //Check if the particle is inside the cone defined by the leading particle
  b = kFALSE;
  
  if(phi > TMath::TwoPi())
    phi-=TMath::TwoPi();
  if(phi < 0.)
    phi+=TMath::TwoPi();
  
  Double_t  rad = 10000 + cone;
  
  if(TMath::Abs(part->Phi()-phi) <= (TMath::TwoPi() - cone))
    rad = TMath::Sqrt(TMath::Power(part->Eta()-eta,2)+
		      TMath::Power(part->Phi()-phi,2));
  else{
    if(part->Phi()-phi > TMath::TwoPi() - cone)
      rad = TMath::Sqrt(TMath::Power(part->Eta()-eta,2)+
			TMath::Power((part->Phi()-TMath::TwoPi())-phi,2));
    if(part->Phi()-phi < -(TMath::TwoPi() - cone))
      rad = TMath::Sqrt(TMath::Power(part->Eta()-eta,2)+
			TMath::Power((part->Phi()+TMath::TwoPi())-phi,2));
  }

  if(rad < cone )
    b = kTRUE;
  
}

//____________________________________________________________________________
void AliAnaParticleJetLeadingCone::FillJetHistos(TClonesArray * pl, Double_t ptg, Double_t ptl, TString type, TString lastname)
{
  //Fill histograms wth jet fragmentation 
  //and number of selected jets and leading particles
  //and the background multiplicity
  TParticle * particle = 0 ;
  Int_t ipr = 0;
  Float_t  charge = 0;
  
  TIter next(pl) ; 
  while ( (particle = dynamic_cast<TParticle*>(next())) ) {
    ipr++ ;
    Double_t pt = particle->Pt();
    
    charge = TDatabasePDG::Instance()
      ->GetParticle(particle->GetPdgCode())->Charge();
    if(charge != 0){//Only jet Charged particles 
      dynamic_cast<TH2F*>
     	(GetOutputContainer()->FindObject(type+"Fragment"+lastname))
     	->Fill(ptg,pt/ptg);
      dynamic_cast<TH2F*>
     	(GetOutputContainer()->FindObject(type+"PtDist"+lastname))
     	->Fill(ptg,pt);
    }//charged
    
  }//while
  
  if(type == "Bkg")
    dynamic_cast<TH1F*>
      (GetOutputContainer()->FindObject("NBkg"+lastname))
      ->Fill(ipr);
  else{
    dynamic_cast<TH1F*>
      (GetOutputContainer()->FindObject("NJet"+lastname))->
      Fill(ptg);
    dynamic_cast<TH2F*>
      (GetOutputContainer()->FindObject("NLeading"+lastname))
      ->Fill(ptg,ptl);
  }
  
}


//____________________________________________________________________________
Bool_t  AliAnaParticleJetLeadingConeCorrelation::SelectCluster(AliAODCaloCluster * calo, Double_t *vertex, TLorentzVector & mom, Int_t & pdg){
   //Select cluster depending on its pid and acceptance selections
   
   //Skip matched clusters with tracks
  if(calo->GetNTracksMatched() > 0) return kFALSE;
   
   //Check PID
   calo->GetMomentum(mom,vertex);//Assume that come from vertex in straight line
   pdg = AliCaloPID::kPhoton;   
   if(IsCaloPIDOn()){
     //Get most probable PID, 2 options check PID weights (in MC this option is mandatory)
     //or redo PID, recommended option for EMCal.
     if(!IsCaloPIDRecalculationOn() || GetReader()->GetDataType() == AliCaloTrackReader::kMC )
       pdg = GetCaloPID()->GetPdg("EMCAL",calo->PID(),mom.E());//PID with weights
     else
       pdg = GetCaloPID()->GetPdg("EMCAL",mom,calo->GetM02(),0,0,0,0);//PID recalculated
     
     if(GetDebug() > 1) printf("PDG of identified particle %d\n",pdg);
     //If it does not pass pid, skip
     if(pdg != AliCaloPID::kPhoton || pdg != AliCaloPID::kPi0) {
       return kFALSE ;
   }
   
   //Check acceptance selection
   if(IsFidutialCutOn()){
     Bool_t in = GetFidutialCut()->IsInFidutialCut(mom,"EMCAL") ;
     if(! in ) return kFALSE ;
   }
   
   if(GetDebug() > 1) printf("cluster selection cuts passed: pT %3.2f, pdg %d\n",mom.Pt(), pdg);
   
   return kTRUE;
   
 }
