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


//_________________________________________________________________________
//  Class for Filling JetFinder Plots
// --
//*-- Author: Mark Horner (LBL/UCT)
// --
// --


#include "TMath.h"
#include "AliEMCALJetFinderPlots.h"

ClassImp(AliEMCALJetFinderPlots)
	
AliEMCALJetFinderPlots::AliEMCALJetFinderPlots() :
  fDebug(0),fConeRadius(0.),fNominalEnergy(0.),fOutput(0),
  fhFragmFcn(0),fhPartonFragmFcn(0),fhPartonJT(0),fhPartonPL(0),
  fhJetJT(0),fhJetPL(0),fhJetEt(0),fhJetEtDiff(0),fhJetEta(0),fhJetPhi(0),
  fhPartonEta(0),fhPartonPhi(0),fhEtaDiff(0),fhPhiDiff(0),fhEtaPhiSpread(0),
  fhNJets(0),fhFragmFcn2(0),fhPartonFragmFcn2(0),fhPartonJT2(0),fhPartonPL2(0),
  fhJetJT2(0),fhJetPL2(0),fhJetEt2(0),fhJetEtDiff2(0),fhJetEta2(0),fhJetPhi2(0),
  fhPartonEta2(0),fhPartonPhi2(0),fhEtaDiff2(0),fhPhiDiff2(0),fhEtaPhiSpread2(0),
  fhNJets2(0),fhJetEtSecond2(0),fhJetEtRatio2(0),fhEtaPhiDist2(0),fhInputOutput(0),
  fhRecoBinPt(0),fhRecoBinPtNoBg(0),fhRecoBinPartonPt(0),fhRecoBinJetEt(0),fhRecoBinInputJetEt(0),
  fhJetPT(0),fhPartonPT(0),fhJetPT2(0),fhPartonPT2(0),fhRecoBinFragmFcn(0),fhRecoBinFragmFcnNoBg(0),
  fhRecoBinPartonFragmFcn(0),fhJetInvE(0),fhJetInvE2(0),fhBackHisto(0),fScaleFactor(0),fInitialised(0)
{
	// Constructor to initialise variables
  fInitialised = kFALSE;
  fNominalEnergy = 0.0;
  fConeRadius = 0.3;
  fDebug = 0;
  fOutput=0;
    fhFragmFcn=0;// = new TH1F("hFragmFcn","Fragmentation Function",100,0,1);
    fhPartonFragmFcn=0;// = new TH1F("hPartonFragmFcn","Fragmentation Function",100,0,1);
   fhPartonJT=0;// = new TH1F("hPartonJT","Track Momentum Perpendicular to Parton Axis",100,0.,10.);
    fhPartonPL=0;// = new TH1F("hPartonPL","Track Momentum Parallel to Parton Axis ",100,0.,100.);
    fhJetJT=0;// = new TH1F("hJetJT","Track Momentum Perpendicular to Jet Axis",100,0.,10.);
    fhJetPL=0;// = new TH1F("hJetPL","Track Momentum Parallel to Jet Axis ",100,0.,100.);
    fhJetEt=0;// = new TH1F("hJetEt","E_{T}^{reco}",250,0.,250.);
    fhJetEta=0;// = new       TH1F("hJetEta","#eta_{jet}^{reco}",180,-0.9,0.9);
    fhJetPhi=0;// = new       TH1F("hJetPhi","#phi_{jet}^{reco}",62,0.,3.1);
    fhPartonEta=0;// = new    TH1F("hPartonEta","#eta_{Parton}",180,-0.9,0.9);
    fhPartonPhi=0;// = new    TH1F("hPartonPhi","#phi_{Parton}",62,0.,3.1);
    fhEtaDiff=0;// = new      TH1F("hEtaDiff","#eta_{jet}^{reco}-#eta_{jet}^{input}",100,-0.5,0.5);
    fhPhiDiff=0;//  = new TH1F("hPhiDiff","#phi_{jet}^{reco}-#phi_{jet}^{input}",100,-0.5,0.5);
    fhNJets=0;// = new        TH1F("hNJets","N Reconstructed jets",11,-0.5,10.5);
  fhEtaPhiSpread=0;	    

fhFragmFcn2=0;	// ("hFragmFcn2","Fragmentation Function",100,0,1);
fhPartonFragmFcn2=0;// ("hFragmFcn2","Parton Fragmentation Function",100,0,1);
fhPartonJT2=0;	// ("hPartonJT2","Track Momentum Perpendicular to Parton Axis",100,0.,10.);
fhPartonPL2=0;	// ("hPartonPL2","Track Momentum Parallel to Parton Axis ",100,0.,100.);
fhJetJT2=0;	// ("hJetJT2","Track Momentum Perpendicular to Jet Axis",100,0.,10.);
fhJetPL2=0;	// ("hJetPL2","Track Momentum Parallel to Jet Axis ",100,0.,100.);
fhJetEt2=0;	// ("hJetEt2","E_{T}^{reco}",250,0.,250.);
fhJetEta2=0;	// ("hJetEta2","#eta_{jet}^{reco}",180,-0.9,0.9);
fhJetPhi2=0;	// ("hJetPhi2","#phi_{jet}^{reco}",62,0.,3.1);
fhPartonEta2=0;	// ("hPartonEta2","#eta_{Parton}",180,-0.9,0.9);
fhPartonPhi2=0;	// ("hPartonPhi2","#phi_{Parton}",62,0.,3.1);
fhEtaDiff2=0;	// ("hEtaDiff2","#eta_{jet}^{reco}-#eta_{jet}^{input}",100,-0.5,0.5);
fhPhiDiff2=0;	// ("hPhiDiff2","#phi_{jet}^{reco}-#phi_{jet}^{input}",100,-0.5,0.5);
fhEtaPhiSpread2=0;	// ("hEtaPhiSpread2","#eta - #phi Distribution 
						//of Reconstructed Jets",192,-0.7,0.7,288,pi/3,pi);
fhNJets2=0;	// ("hNJets2","N Reconstructed jets",11,-0.5,10.5);
fhJetEtSecond2=0; //("hJetEtSecond2","E_{T}^{reco}",250,0.,250.); 
fhJetEtRatio2=0;  //("hJetEtRatio2","Ratio of Second Highest to Highest",100,0,1);
fhEtaPhiDist2=0;  //("hEtaPhiDist2","Angular Distance Between First and Second",100,0,3);
fhInputOutput=0;
//     TH2F                            *fhInputOutput;  //("hJetEtRatio2","Ratio of Second Highest to Highest",100,0,1);
     
fhRecoBinPt=0;	       // ("fhRecoBinPt","Reconstructed Pt Distribution",100,0,1);
fhRecoBinPtNoBg=0;	       // ("fhRecoBinPt","Reconstructed Pt Distribution",100,0,1);
fhRecoBinPartonPt=0;    // ("fhRecoBinPartonPt","Input Pt Distribution",100,0,1);
fhRecoBinJetEt=0;       // ("fhRecoJetEt","E_{T}^{reco}",250,0.,250.);
fhRecoBinInputJetEt=0;  // ("fhRecoInputJetEt","E_{T}^{reco}",250,0.,250.);
fhJetPT =0;// new TH1F("hJetPT","P_{T} Distribution",200,0,200);
fhPartonPT =0;// new TH1F("hPartonPT","Parton P_{T} Distribution",200,0,1);
fhJetPT2 =0;// new TH1F("hJetPT","P_{T} Distribution",200,0,200);
fhPartonPT2 =0;// new TH1F("hPartonPT","Parton P_{T} Distribution",200,0,1);
fhRecoBinFragmFcn =0;//new TH1F("fhRecoBinFragmFcn","Reconstructed Frag. Fcn",100,0,1);
fhRecoBinFragmFcnNoBg =0;//new TH1F("fhRecoBinFragmFcn","Reconstructed Frag. Fcn",100,0,1);
fhRecoBinPartonFragmFcn =0;// new TH1F("fhRecoBinPartonFragmFcn","Input Bin Fragm Fcn Distribution",100,0,1);
fhJetInvE=0;// = new TH1F("fhJetInvE","#frac{1}{E_{R}}",100,0,1);
fhJetInvE2=0;// = new TH1F("fhJetInvE","#frac{1}{E_{R}}",100,0,1);
fScaleFactor = 1.0/0.6731;
fhBackHisto=0;

}

AliEMCALJetFinderPlots::AliEMCALJetFinderPlots(const AliEMCALJetFinderPlots& jfp) 
  : TObject(jfp), fDebug(jfp.fDebug),fConeRadius(jfp.fConeRadius),fNominalEnergy(jfp.fNominalEnergy),
    fOutput(jfp.fOutput),fhFragmFcn(jfp.fhFragmFcn),fhPartonFragmFcn(jfp.fhPartonFragmFcn),
    fhPartonJT(jfp.fhPartonJT),fhPartonPL(jfp.fhPartonPL),fhJetJT(jfp.fhJetJT),fhJetPL(jfp.fhJetPL),
    fhJetEt(jfp.fhJetEt),fhJetEtDiff(jfp.fhJetEtDiff),fhJetEta(jfp.fhJetEta),fhJetPhi(jfp.fhJetPhi),
    fhPartonEta(jfp.fhPartonEta),fhPartonPhi(jfp.fhPartonPhi),fhEtaDiff(jfp.fhEtaDiff),
    fhPhiDiff(jfp.fhPhiDiff),fhEtaPhiSpread(jfp.fhEtaPhiSpread),fhNJets(jfp.fhNJets),
    fhFragmFcn2(jfp.fhFragmFcn2),fhPartonFragmFcn2(jfp.fhPartonFragmFcn2),fhPartonJT2(jfp.fhPartonJT2),
    fhPartonPL2(jfp.fhPartonPL2),fhJetJT2(jfp.fhJetJT2),fhJetPL2(jfp.fhJetPL2),fhJetEt2(jfp.fhJetEt2),
    fhJetEtDiff2(jfp.fhJetEtDiff2),fhJetEta2(jfp.fhJetEta2),fhJetPhi2(jfp.fhJetPhi2), fhPartonEta2(jfp.fhPartonEta2),
    fhPartonPhi2(jfp.fhPartonPhi2),fhEtaDiff2(jfp.fhEtaDiff2),fhPhiDiff2(jfp.fhPhiDiff2),
    fhEtaPhiSpread2(jfp.fhEtaPhiSpread2),fhNJets2(jfp.fhNJets2),fhJetEtSecond2(jfp.fhJetEtSecond2),
    fhJetEtRatio2(jfp.fhJetEtRatio2),fhEtaPhiDist2(jfp.fhEtaPhiDist2),fhInputOutput(jfp.fhInputOutput),
    fhRecoBinPt(jfp.fhRecoBinPt),fhRecoBinPtNoBg(jfp.fhRecoBinPtNoBg),fhRecoBinPartonPt(jfp.fhRecoBinPartonPt),
    fhRecoBinJetEt(jfp.fhRecoBinJetEt),fhRecoBinInputJetEt(jfp.fhRecoBinInputJetEt),fhJetPT(jfp.fhJetPT),
    fhPartonPT(jfp.fhPartonPT),fhJetPT2(jfp.fhJetPT2),fhPartonPT2(jfp.fhPartonPT2),fhRecoBinFragmFcn(jfp.fhRecoBinFragmFcn),
    fhRecoBinFragmFcnNoBg(jfp.fhRecoBinFragmFcnNoBg),fhRecoBinPartonFragmFcn(jfp.fhRecoBinPartonFragmFcn),
    fhJetInvE(jfp.fhJetInvE),fhJetInvE2(jfp.fhJetInvE2),fhBackHisto(jfp.fhBackHisto),fScaleFactor(jfp.fScaleFactor),
    fInitialised(jfp.fInitialised)
{
  //copy ctor
}

void AliEMCALJetFinderPlots::InitPlots()
{
//========================= CASE 1 =======================================	
    fhFragmFcn = new TH1F("hFragmFcn","Fragmentation Function",200,0,2);
    fhFragmFcn->Sumw2();
    fhJetPT = new TH1F("hJetPT","P_{T} Distribution",200,0,200);
    fhJetPT->Sumw2();
    fhPartonPT = new TH1F("hPartonPT","Parton P_{T} Distribution",200,0,200);
    fhPartonPT->Sumw2();
    fhPartonFragmFcn = new TH1F("hPartonFragmFcn","Parton Fragmentation Function",200,0,2);
    fhPartonFragmFcn->Sumw2();
    fhPartonJT = new TH1F("hPartonJT","Track Momentum Perpendicular to Parton Axis",100,0.,10.);
    fhPartonJT->Sumw2();
    fhPartonPL = new TH1F("hPartonPL","Track Momentum Parallel to Parton Axis ",100,0.,100.);
    fhPartonPL->Sumw2();
    fhJetJT = new TH1F("hJetJT","Track Momentum Perpendicular to Jet Axis",100,0.,10.);
    fhJetJT->Sumw2();
    fhJetPL = new TH1F("hJetPL","Track Momentum Parallel to Jet Axis ",100,0.,100.);
    fhJetPL->Sumw2();
    fhJetEt = new TH1F("hJetEt","E_{T}^{reco}",250,0.,250.);
    fhJetEt->Sumw2();
    fhJetEtDiff = new TH1F("hJetEtDiff","E_{T}^{reco}-E_{T}^{Parton}",250,-124.,125.);
    fhJetEtDiff->Sumw2();
    fhJetEta = new	TH1F("hJetEta","#eta_{jet}^{reco}",180,-0.9,0.9);
    fhJetEta->Sumw2();
    fhJetPhi = new 	TH1F("hJetPhi","#phi_{jet}^{reco}",62,0.,3.1);
    fhJetPhi->Sumw2();
    fhPartonEta = new	TH1F("hPartonEta","#eta_{Parton}",180,-0.9,0.9);
    fhPartonEta->Sumw2();
    fhPartonPhi = new 	TH1F("hPartonPhi","#phi_{Parton}",62,0.,3.1);
    fhPartonPhi->Sumw2();
    fhEtaDiff = new	TH1F("hEtaDiff","#eta_{jet}^{reco}-#eta_{jet}^{input}",100,-0.5,0.5);
    fhEtaDiff->Sumw2();
    fhPhiDiff  = new TH1F("hPhiDiff","#phi_{jet}^{reco}-#phi_{jet}^{input}",100,-0.5,0.5);
    fhPhiDiff->Sumw2();
    fhNJets = new	TH1F("hNJets","N Reconstructed jets",11,-0.5,10.5);
    fhNJets->Sumw2();
    fhEtaPhiSpread = new TH2F("hEtaPhiSpread","#eta - #phi Distribution of Reconstructed Jets",100,-0.5,0.5,100,-0.5,0.5);
    fhEtaPhiSpread->Sumw2();    
  fhNJets->SetXTitle("N_{jets}^{reco}/event");
  fhNJets->SetYTitle("N_{events}");

  //Jet properties
  fhJetEt->SetFillColor(16);
  fhJetEt->SetXTitle("E_{T}^{reco}");
  
  fhJetEta->SetFillColor(16);
  fhJetEta->SetXTitle("#eta_{jet}^{reco}");
  
  fhJetPhi->SetFillColor(16);
  fhJetPhi->SetXTitle("#phi_{jet}^{reco}");
  
  fhPartonEta->SetFillColor(16);
  fhPartonEta->SetXTitle("#eta_{parton}");
  
  fhPartonPhi->SetFillColor(16);
  fhPartonPhi->SetXTitle("#phi_{parton}");

  fhPartonPL->SetXTitle("p (GeV/c)");
  fhPartonJT->SetXTitle("p (GeV/c)");

  fhPartonFragmFcn->SetXTitle("Z = p_{T}^{Chg}/E_{T}^{parton}");

  //Jet component properties

  fhJetPL->SetXTitle("p (GeV/c)");
  fhJetJT->SetXTitle("p (GeV/c)");
  fhFragmFcn->SetXTitle("Z = p_{T}^{Chg}/E_{T}^{reco}");
  fhPartonFragmFcn->SetXTitle("Z = p_{T}^{Chg}/E_{T}^{reco}");

  fhEtaDiff->SetXTitle("#eta_{jet}^{reco}-#eta_{jet}^{input}");
  fhPhiDiff->SetXTitle("#phi_{jet}^{reco}-#phi_{jet}^{input}");
  fhEtaPhiSpread->SetXTitle("#eta");
  fhEtaPhiSpread->SetYTitle("#phi");

//======================= CASE 2 ======================================
  

fhFragmFcn2 		= new TH1F("hFragmFcn2","Fragmentation Function",200,0,2);
fhFragmFcn2->Sumw2();
    fhJetPT2 = new TH1F("hJetPT2","P_{T} Distribution",200,0,200);
    fhJetPT2->Sumw2();
    fhPartonPT2 = new TH1F("hPartonPT2","Parton P_{T} Distribution",200,0,1);
    fhPartonPT2->Sumw2();
fhPartonFragmFcn2 	= new TH1F("hPartonFragmFcn2","Parton Fragmentation Function",200,0,2);
fhPartonFragmFcn2->Sumw2();
fhPartonJT2 		= new TH1F("hPartonJT2","Track Momentum Perpendicular to Parton Axis",100,0.,10.);
fhPartonJT2->Sumw2();
fhPartonPL2		= new TH1F("hPartonPL2","Track Momentum Parallel to Parton Axis ",100,0.,100.);
fhPartonPL2->Sumw2();
fhJetJT2			= new TH1F("hJetJT2","Track Momentum Perpendicular to Jet Axis",100,0.,10.);
fhJetJT2->Sumw2();
fhJetPL2			= new TH1F("hJetPL2","Track Momentum Parallel to Jet Axis ",100,0.,100.);
fhJetPL2->Sumw2();
fhJetEt2			= new TH1F("hJetEt2","E_{T}^{reco}",250,0.,250.);
fhJetEt2->Sumw2();
    fhJetEtDiff2 = new TH1F("hJetEtDiff2","E_{T}^{reco}-E_{T}^{Parton}",250,-124.,125.);
    fhJetEtDiff2->Sumw2();
fhJetEta2		= new TH1F("hJetEta2","#eta_{jet}^{reco}",180,-0.9,0.9);
fhJetEta2->Sumw2();
fhJetPhi2		= new TH1F("hJetPhi2","#phi_{jet}^{reco}",62,0.,3.1);
fhJetPhi2->Sumw2();
fhPartonEta2		= new TH1F("hPartonEta2","#eta_{Parton}",180,-0.9,0.9);
fhPartonEta2->Sumw2();
fhPartonPhi2		= new TH1F("hPartonPhi2","#phi_{Parton}",62,0.,3.1);
fhPartonPhi2->Sumw2();
fhEtaDiff2		= new TH1F("hEtaDiff2","#eta_{jet}^{reco}-#eta_{jet}^{input}",100,-0.5,0.5);
fhEtaDiff2->Sumw2();
fhPhiDiff2		= new TH1F("hPhiDiff2","#phi_{jet}^{reco}-#phi_{jet}^{input}",100,-0.5,0.5);
fhPhiDiff2->Sumw2();
fhEtaPhiSpread2 = new TH2F("hEtaPhiSpread2","#eta - #phi Distribution of Reconstructed Jets",100,-0.5,0.5,100,-0.5,0.5);
fhEtaPhiSpread2->Sumw2();
fhNJets2			= new TH1F("hNJets2","N Reconstructed jets",11,-0.5,10.5);
fhNJets2->Sumw2();
fhJetEtSecond2		= new TH1F("hJetEtSecond2","E_{T}^{reco}",250,0.,250.); 
fhJetEtSecond2->Sumw2();
fhJetEtRatio2		= new TH1F("hJetEtRatio2","Ratio of Second Highest to Highest",100,0,1);
fhJetEtRatio2->Sumw2();
fhEtaPhiDist2		= new TH1F("hEtaPhiDist2","Angular Distance Between First and Second",100,0,3);
fhEtaPhiDist2->Sumw2();

fhInputOutput=	new TH2F("hInputOutput","Input and Reconstruction Correlations;Input;Output",200,0,200,200,0,200);  //("hJetEtRatio2","Ratio of Second Highest to Highest",100,0,1);

//============================== Reconstruction Bin Comparison  ============================================

fhRecoBinPt =new TH1F("fhRecoBinPt","Reconstructed Pt Distribution",200,0,200);
fhRecoBinPt->Sumw2();
fhRecoBinPtNoBg =new TH1F("fhRecoBinPtNoBg","Reconstructed Pt Distribution Background Subtracted",200,0,200);
fhRecoBinPtNoBg->Sumw2();
fhRecoBinPartonPt = new TH1F("fhRecoBinPartonPt","Input Pt Distribution",200,0,200);
fhRecoBinPartonPt->Sumw2();
fhRecoBinFragmFcn =new TH1F("fhRecoBinFragmFcn","Reconstructed Frag. Fcn",200,0,2);
fhRecoBinFragmFcn->Sumw2();
fhRecoBinFragmFcnNoBg =new TH1F("fhRecoBinFragmFcnNoBg","Reconstructed Frag. Fcn With Background Removed",200,0,2);
fhRecoBinFragmFcnNoBg->Sumw2();
fhRecoBinPartonFragmFcn = new TH1F("fhRecoBinPartonFragmFcn","Input Bin Fragm Fcn Distribution",200,0,2);
fhRecoBinPartonFragmFcn->Sumw2();
fhRecoBinJetEt = new TH1F("fhRecoJetEt","E_{T}^{reco}",250,0.,250.);
fhRecoBinJetEt->Sumw2();
fhRecoBinInputJetEt = new TH1F("fhRecoInputJetEt","E_{T}^{reco}",250,0.,250.);
fhRecoBinInputJetEt->Sumw2();
        

fhJetInvE = new TH1F("fhJetInvE","#frac{1}{E_{R}}",100,0,1);
fhJetInvE->Sumw2();
fhJetInvE2 = new TH1F("fhJetInvE2","#frac{1}{E_{R}}",100,0,1);
fhJetInvE2->Sumw2();
		


  fInitialised = kTRUE;	
  
}

AliEMCALJetFinderPlots::~AliEMCALJetFinderPlots()
{
	// To ensure that all requested memory is returned
delete    fhFragmFcn;// = new TH1F("hFragmFcn","Fragmentation Function",100,0,1);
delete    fhPartonFragmFcn;// = new TH1F("hFragmFcn","Fragmentation Function",100,0,1);
delete   fhPartonJT;// = new TH1F("hPartonJT","Track Momentum Perpendicular to Parton Axis",100,0.,10.);
delete    fhPartonPL;// = new TH1F("hPartonPL","Track Momentum Parallel to Parton Axis ",100,0.,100.);
delete    fhJetJT;// = new TH1F("hJetJT","Track Momentum Perpendicular to Jet Axis",100,0.,10.);
delete    fhJetPL;// = new TH1F("hJetPL","Track Momentum Parallel to Jet Axis ",100,0.,100.);
delete    fhJetEt;// = new TH1F("hJetEt","E_{T}^{reco}",250,0.,250.);
delete	fhJetEtDiff;	// ("hJetEt2","E_{T}^{reco}",250,0.,250.);
delete    fhJetEta;// = new       TH1F("hJetEta","#eta_{jet}^{reco}",180,-0.9,0.9);
delete    fhJetPhi;// = new       TH1F("hJetPhi","#phi_{jet}^{reco}",62,0.,3.1);
delete    fhPartonEta;// = new    TH1F("hPartonEta","#eta_{Parton}",180,-0.9,0.9);
delete    fhPartonPhi;// = new    TH1F("hPartonPhi","#phi_{Parton}",62,0.,3.1);
delete    fhEtaDiff;// = new      TH1F("hEtaDiff","#eta_{jet}^{reco}-#eta_{jet}^{input}",100,-0.5,0.5);
delete    fhPhiDiff;//  = new TH1F("hPhiDiff","#phi_{jet}^{reco}-#phi_{jet}^{input}",100,-0.5,0.5);
delete    fhNJets;// = new        TH1F("hNJets","N Reconstructed jets",11,-0.5,10.5);
delete 	  fhEtaPhiSpread;	    

	delete				fhFragmFcn2;	// ("hFragmFcn2","Fragmentation Function",100,0,1);
	delete				fhPartonFragmFcn2;// ("hFragmFcn2","Parton Fragmentation Function",100,0,1);
	delete				fhPartonJT2;	// ("hPartonJT2","Track Momentum Perpendicular to Parton Axis",100,0.,10.);
	delete				fhPartonPL2;	// ("hPartonPL2","Track Momentum Parallel to Parton Axis ",100,0.,100.);
	delete				fhJetJT2;	// ("hJetJT2","Track Momentum Perpendicular to Jet Axis",100,0.,10.);
	delete				fhJetPL2;	// ("hJetPL2","Track Momentum Parallel to Jet Axis ",100,0.,100.);
	delete				fhJetEt2;	// ("hJetEt2","E_{T}^{reco}",250,0.,250.);
	delete				fhJetEtDiff2;	// ("hJetEt2","E_{T}^{reco}",250,0.,250.);
	delete				fhJetEta2;	// ("hJetEta2","#eta_{jet}^{reco}",180,-0.9,0.9);
	delete				fhJetPhi2;	// ("hJetPhi2","#phi_{jet}^{reco}",62,0.,3.1);
	delete				fhPartonEta2;	// ("hPartonEta2","#eta_{Parton}",180,-0.9,0.9);
	delete				fhPartonPhi2;	// ("hPartonPhi2","#phi_{Parton}",62,0.,3.1);
	delete				fhEtaDiff2;	// ("hEtaDiff2","#eta_{jet}^{reco}-#eta_{jet}^{input}",100,-0.5,0.5);
	delete				fhPhiDiff2;	// ("hPhiDiff2","#phi_{jet}^{reco}-#phi_{jet}^{input}",100,-0.5,0.5);
	delete				fhEtaPhiSpread2;	// ("hEtaPhiSpread2","#eta - #phi Distribution 
							//of Reconstructed Jets",192,-0.7,0.7,288,pi/3,pi);
	delete				fhNJets2;	// ("hNJets2","N Reconstructed jets",11,-0.5,10.5);
	delete				fhJetEtSecond2; //("hJetEtSecond2","E_{T}^{reco}",250,0.,250.); 
	delete				fhJetEtRatio2;  //("hJetEtRatio2","Ratio of Second Highest to Highest",100,0,1);
	delete 				fhEtaPhiDist2;  //("hEtaPhiDist2","Angular Distance Between First and Second",100,0,3);

	delete fhRecoBinPt;          // ("fhRecoBinPt","Reconstructed Pt Distribution",100,0,1);
	delete fhRecoBinPtNoBg;          // ("fhRecoBinPt","Reconstructed Pt Distribution",100,0,1);
	delete fhRecoBinPartonPt;    // ("fhRecoBinPartonPt","Input Pt Distribution",100,0,1);
	delete fhRecoBinJetEt;       // ("fhRecoJetEt","E_{T}^{reco}",250,0.,250.);
	delete fhRecoBinInputJetEt;  // ("fhRecoInputJetEt","E_{T}^{reco}",250,0.,250.);
	                                
 	delete fhJetPT ;// new TH1F("hJetPT","P_{T} Distribution",200,0,200);
	delete fhPartonPT ;// new TH1F("hPartonPT","Parton P_{T} Distribution",200,0,1);
 	delete fhJetPT2 ;// new TH1F("hJetPT","P_{T} Distribution",200,0,200);
	delete fhPartonPT2 ;// new TH1F("hPartonPT","Parton P_{T} Distribution",200,0,1);
	delete fhRecoBinFragmFcn;//new TH1F("fhRecoBinFragmFcn","Reconstructed Frag. Fcn",100,0,1);
	delete fhRecoBinFragmFcnNoBg;//new TH1F("fhRecoBinFragmFcn","Reconstructed Frag. Fcn",100,0,1);
	delete fhRecoBinPartonFragmFcn;// new TH1F("fhRecoBinPartonFragmFcn","Input Bin Fragm Fcn Distribution",100,0,1);
    delete fhJetInvE;// = new TH1F("fhJetInvE","#frac{1}{E_{R}}",100,0,1);
    delete fhJetInvE2;// = new TH1F("fhJetInvE","#frac{1}{E_{R}}",100,0,1);

}	

void AliEMCALJetFinderPlots::FillFromOutput(AliEMCALJetFinderOutput* output, Float_t weight)
{
	// Fill histograms from an output object
if (!fInitialised) InitPlots();	
  fOutput = output;
  if (!fOutput) return;
// Make some temporary histograms to make sure we subtract 
  // background properly
/*
tempFragmFcnNoBg =new TH1F("tempFragmFcnNoBg","Reconstructed Frag. Fcn With Background Removed",200,0,2);
tempPtNoBg =new TH1F("tempPtNoBg","Reconstructed Frag. Fcn With Background Removed",200,0,200);
tempFragmFcnNoBg->Fill(count/(jethighest->Energy()*fScaleFactor),-fhBackHisto->GetBinContent(count));
tempPtNoBg->AddBinContent(count,-fhBackHisto->GetBinContent(count));
*/
  
  fhNJets->Fill(fOutput->GetNJets());
  Bool_t doesJetMeetBinCriteria = 0;
  AliEMCALJet* jethighest=0; 
  AliEMCALJet* jetsecond=0; 
  //  Find Highest and Second Highest Jet 
  // (NB!!!!!!!)  Pointing into the EMCAL!!!!!!
  
// =========================== All cases ===================================
 

	// I will make a little array of jet indices for which jets are in 
	// the EMCAL then my counter can loop from 0 below - but it will 
	// be the index of the array of applicable jets
    
  Int_t appjet[4];
  Int_t numappjet=0;
  
  for (Int_t appc=0;appc<fOutput->GetNJets();appc++)
  { // Check all jets for applicability
	  Float_t eta = fOutput->GetJet(appc)->Eta();
	  Float_t phi = fOutput->GetJet(appc)->Phi();
	  if (eta > -0.7 && eta < 0.7 && phi > 1./3.*TMath::Pi() && phi < TMath::Pi())
	  { // Then jet is applicable
		  appjet[numappjet]=appc;				
		  numappjet++;
	  }
  }
	
  
Float_t et=0;
if (numappjet >=1)
{
	Float_t theta = 2.0*atan(exp(-fOutput->GetParton(0)->Eta()));
	et    = fOutput->GetParton(0)->Energy() * TMath::Sin(theta);
	if (fOutput->GetNJets()>1)
	{
	      	for (Int_t counter = 0; counter<numappjet;counter++)  
		{	  
			if (counter==0)
			{ 		  
				jethighest = fOutput->GetJet(appjet[0]);
	      			jetsecond  = fOutput->GetJet(appjet[1]);
	      		}
	      		if (counter>0)
	      		{
	      			Float_t energyhighest = jethighest->Energy();
	      			Float_t energysecond = jetsecond->Energy();
	      			
	      			if ((fOutput->GetJet(appjet[counter]))->Energy()>energyhighest) 
	      			{
	      				jetsecond=jethighest;
	      				jethighest=fOutput->GetJet(appjet[counter]);
	      			}else if ((fOutput->GetJet(appjet[counter]))->Energy()>energysecond) 
	      			{
	      				jetsecond=fOutput->GetJet(appjet[counter]);
	      			}
	      
			}
	      	}
	}else
	{
		Float_t eta = fOutput->GetJet(0)->Eta();
		Float_t phi = fOutput->GetJet(0)->Phi();
		if (eta > -0.7 && eta < 0.7 && phi > 1./3.*TMath::Pi() && phi < TMath::Pi())
		{ // Then jet is applicable
			jethighest=fOutput->GetJet(0);
			jetsecond=0;
		}else
		{
			Error("FillFromOutput","There is only one jet and it isn't in the area of applicability");
		}
						
	}
	if ( 95.0 < jethighest->Energy()*fScaleFactor && jethighest->Energy()*fScaleFactor < 105.0  )
	{
		doesJetMeetBinCriteria = 1;
		fhRecoBinJetEt->Fill(jethighest->Energy()*fScaleFactor,weight);
		fhRecoBinInputJetEt->Fill(et,weight);
	}
	fhInputOutput->Fill(et,jethighest->Energy());
		
}

if (numappjet > 1)
{	
//========================= CASE 2 ===========================	
  Int_t nPartons = fOutput->GetNPartons();
  fhNJets2->Fill(fOutput->GetNJets());
  AliEMCALParton* parton;

  // End finding highest and second highest and continue
  fhJetEt2->Fill(jethighest->Energy()*fScaleFactor,weight);
  fhJetEtDiff2->Fill(jethighest->Energy()*fScaleFactor-et,weight);
  fhJetInvE2->Fill(1.0/(jethighest->Energy()*fScaleFactor),weight);
  fhJetEta2->Fill(jethighest->Eta(),weight  );
  fhJetPhi2->Fill(jethighest->Phi(),weight );
  if (nPartons ==0) return;
  parton = fOutput->GetParton(0);
  
  fhPartonEta2->Fill( parton->Eta(),weight );
  fhPartonPhi2->Fill( parton->Phi(),weight );

  //hJetEtDiff->Fill( jet->Energy() - parton->Energy() );
  fhEtaDiff2->Fill( jethighest->Eta() - parton->Eta(),weight );
  fhPhiDiff2->Fill( jethighest->Phi() - parton->Phi(),weight );
  fhEtaPhiSpread2->Fill(jethighest->Eta()-parton->Eta(),jethighest->Phi() - parton->Phi());
  fhJetEtSecond2->Fill(jetsecond->Energy()*fScaleFactor,weight); 
  fhJetEtRatio2->Fill(jetsecond->Energy()/jethighest->Energy(),weight); 
  fhEtaPhiDist2->Fill( TMath::Sqrt((jethighest->Eta() - jetsecond->Eta())*(jethighest->Eta() - jetsecond->Eta())
		      + (jethighest->Phi() - jetsecond->Phi())*(jethighest->Phi() - jetsecond->Phi())	  ),weight);  
  /* 
  Float_t *pt,*phi,*eta;
  Int_t *pdg;
  pt  = new Float_t[parton->GetNTracks()];
  eta = new Float_t[parton->GetNTracks()];
  phi = new Float_t[parton->GetNTracks()];
  pdg = new Int_t[parton->GetNTracks()];*/



 Float_t pt[2000];
 Float_t eta[2000];
 Float_t phi[2000];
 Int_t pdg[2000];
  
  parton->GetTrackList(pt,eta,phi,pdg);
  for(Int_t iT=0; iT< parton->GetNTracks() ; iT++ ) 
  {
      if ( (eta[iT]-parton->Eta())*(eta[iT]-parton->Eta())+
           (phi[iT]-parton->Phi())*(phi[iT]-parton->Phi()) >fConeRadius * fConeRadius ) continue; 
      Double_t tt = 2.0*atan(exp(-eta[iT])); // These names are short to make the equation manageable
      Double_t rt = 2.0*atan(exp(-parton->Eta()));
      Double_t ctt = cos(tt);
      Double_t crt = cos(rt);
      Double_t stt = sin(tt);
      Double_t srt = sin(rt);
      Double_t ctp = cos(phi[iT]);
      Double_t crp = cos(parton->Phi());
      Double_t stp = sin(phi[iT]);
      Double_t srp = sin(parton->Phi());
      //Double_t alpha = acos(crp*ctp*srt*stt+srp*stp*srt*stt+crt*ctt);
	  Double_t alpha;   
	  if (TMath::Abs(crp*ctp*srt*stt+srp*stp*srt*stt+crt*ctt) > 0.9990)
	  {
		  alpha = 0.0;
	  }else
	  {
		  alpha = TMath::ACos(crp*ctp*srt*stt+srp*stp*srt*stt+crt*ctt);   
	  }
      Double_t correctp = pt[iT]/stt; 	
      fhPartonPL2->Fill( correctp*cos(alpha),weight);
      if ( (parton->Eta()-eta[iT])*(parton->Eta()-eta[iT]) +  
      (parton->Phi()-phi[iT])*(parton->Phi()-phi[iT]) < 0.2*0.2 ) 
    	      fhPartonJT2->Fill( correctp*sin(alpha),weight);
      fhPartonPT2->Fill(correctp*sin(tt),weight);
      if (fNominalEnergy == 0.0) {
	      fhPartonFragmFcn2->Fill(  correctp*sin(tt)/parton->Energy(),weight  );
      }else 
      {
              fhPartonFragmFcn2->Fill(correctp*sin(tt)/fNominalEnergy,weight);
      }    	  
      if (doesJetMeetBinCriteria)
      {
	      fhRecoBinPartonPt->Fill(correctp*sin(tt),weight);          // ("fhRecoBinPt","Reconstructed Pt Distribution",100,0,1);
      }
  }// loop over tracks

/*
  pt  = new Float_t[jet->NTracks()];
  eta = new Float_t[jet->NTracks()];
  phi = new Float_t[jet->NTracks()];
  pdg = new Int_t[jet->NTracks()];*/
  jethighest->TrackList(pt,eta,phi,pdg);
  for(Int_t iT=0; iT< jethighest->NTracks() ; iT++ )
  {
      	  Double_t tt = 2.0*atan(exp(-eta[iT])); // These names are short to make the equation manageable
	  Double_t rt = 2.0*atan(exp(-jethighest->Eta()));
	  Double_t ctt = cos(tt);                   
	  Double_t crt = cos(rt);                         
	  Double_t stt = sin(tt);                               
	  Double_t srt = sin(rt);
	  Double_t ctp = cos(phi[iT]);
	  Double_t crp = cos(jethighest->Phi());
	  Double_t stp = sin(phi[iT]);
	  Double_t srp = sin(jethighest->Phi());
	 // Double_t alpha = acos(crp*ctp*srt*stt+srp*stp*srt*stt+crt*ctt);   
	  Double_t alpha;   
	  if (TMath::Abs(crp*ctp*srt*stt+srp*stp*srt*stt+crt*ctt) > 0.9990)
	  {
		  alpha = 0.0;
	  }else
	  {
		  alpha = TMath::ACos(crp*ctp*srt*stt+srp*stp*srt*stt+crt*ctt);   
	  }
	  Double_t correctp = pt[iT]/stt;
	  fhJetPL2->Fill( correctp*cos(alpha),weight);      
	  if ( (jethighest->Eta()-eta[iT])*(jethighest->Eta()-eta[iT]) +  
			  (jethighest->Phi()-phi[iT])*(jethighest->Phi()-phi[iT]) < 0.2*0.2 )   
		  fhJetJT2->Fill( correctp*sin(alpha),weight);
	  fhJetPT2->Fill(correctp*sin(tt),weight);
	  if (fNominalEnergy==0.0){
		  fhFragmFcn2->Fill(  correctp*sin(tt)/(jethighest->Energy()*fScaleFactor),weight  );
	  } else
	  {
                  fhFragmFcn2->Fill(  correctp*sin(tt)/fNominalEnergy,weight );
	  }
    	  if (doesJetMeetBinCriteria)
    	  {
	  	  fhRecoBinPt->Fill(correctp*sin(tt),weight);          // ("fhRecoBinPt","Reconstructed Pt Distribution",100,0,1);
	  	  fhRecoBinPtNoBg->Fill(correctp*sin(tt),weight);          // ("fhRecoBinPt","Reconstructed Pt Distribution",100,0,1);
		  fhRecoBinFragmFcn->Fill(  correctp*sin(tt)/(jethighest->Energy()*fScaleFactor),weight  ); // This is the jet fragmentation function
		  fhRecoBinFragmFcnNoBg->Fill(  correctp*sin(tt)/(jethighest->Energy()*fScaleFactor),weight  ); // This is the jet fragmentation function
    	  }
  }// loop over tracks
  }

  if (numappjet == 1)
  {

//========================= CASE 1 ===========================	
  Int_t nPartons = fOutput->GetNPartons();
  if (fOutput->GetNJets()!=1) return;
  AliEMCALParton* parton;
  AliEMCALJet* jet; 
  jet = jethighest;//fOutput->GetJet(0);
  fhJetEt->Fill(jet->Energy()*fScaleFactor,weight);
  fhJetEtDiff->Fill(jethighest->Energy()*fScaleFactor-et,weight);
  fhJetInvE->Fill(1.0/(jethighest->Energy()*fScaleFactor),weight);
  fhJetEta->Fill(jet->Eta(),weight  );
  fhJetPhi->Fill(jet->Phi(),weight );
  if (nPartons ==0) return;
  parton = fOutput->GetParton(0);
  
  fhPartonEta->Fill( parton->Eta(),weight );
  fhPartonPhi->Fill( parton->Phi(),weight );

  //hJetEtDiff->Fill( jet->Energy() - parton->Energy() );
  fhEtaDiff->Fill( jet->Eta() - parton->Eta(),weight );
  fhPhiDiff->Fill( jet->Phi() - parton->Phi(),weight );
  fhEtaPhiSpread->Fill(jet->Eta()-parton->Eta(),jet->Phi() - parton->Phi());
 /* 
  Float_t *pt,*phi,*eta;
  Int_t *pdg;
  pt  = new Float_t[parton->GetNTracks()];
  eta = new Float_t[parton->GetNTracks()];
  phi = new Float_t[parton->GetNTracks()];
  pdg = new Int_t[parton->GetNTracks()];*/



 Float_t pt[2000];
 Float_t eta[2000];
 Float_t phi[2000];
 Int_t pdg[2000];
  
  parton->GetTrackList(pt,eta,phi,pdg);
  for(Int_t iT=0; iT< parton->GetNTracks() ; iT++ ) 
  {
      if ( (eta[iT]-parton->Eta())*(eta[iT]-parton->Eta())+
           (phi[iT]-parton->Phi())*(phi[iT]-parton->Phi()) >fConeRadius * fConeRadius ) continue; 
      Double_t tt = 2.0*atan(exp(-eta[iT])); // These names are short to make the equation manageable
      Double_t rt = 2.0*atan(exp(-parton->Eta()));
      Double_t ctt = cos(tt);
      Double_t crt = cos(rt);
      Double_t stt = sin(tt);
      Double_t srt = sin(rt);
      Double_t ctp = cos(phi[iT]);
      Double_t crp = cos(parton->Phi());
      Double_t stp = sin(phi[iT]);
      Double_t srp = sin(parton->Phi());
     // Double_t alpha = acos(crp*ctp*srt*stt+srp*stp*srt*stt+crt*ctt);
	  Double_t alpha;   
	  if (TMath::Abs(crp*ctp*srt*stt+srp*stp*srt*stt+crt*ctt) > 0.9990)
	  {
		  alpha = 0.0;
	  }else
	  {
		  alpha = TMath::ACos(crp*ctp*srt*stt+srp*stp*srt*stt+crt*ctt); 	  }
      Double_t correctp = pt[iT]/stt; 	
      fhPartonPL->Fill( correctp*cos(alpha),weight);
      if ( (parton->Eta()-eta[iT])*(parton->Eta()-eta[iT]) +  
      (parton->Phi()-phi[iT])*(parton->Phi()-phi[iT]) < 0.2*0.2 ) 
    	      fhPartonJT->Fill( correctp*sin(alpha),weight);
      if (fNominalEnergy == 0.0) {
	      fhPartonFragmFcn->Fill(  correctp*sin(tt)/parton->Energy(),weight  );
	  fhPartonPT->Fill(correctp*sin(tt),weight);
      }else 
      {
              fhPartonFragmFcn->Fill(correctp*sin(tt)/fNominalEnergy,weight);
      }
      if (doesJetMeetBinCriteria)
      {
	      fhRecoBinPartonPt->Fill(correctp*sin(tt),weight);          // ("fhRecoBinPt","Reconstructed Pt Distribution",100,0,1);
      }
  }// loop over tracks

/*
  pt  = new Float_t[jet->NTracks()];
  eta = new Float_t[jet->NTracks()];
  phi = new Float_t[jet->NTracks()];
  pdg = new Int_t[jet->NTracks()];*/
  jet->TrackList(pt,eta,phi,pdg);
  for(Int_t iT=0; iT< jet->NTracks() ; iT++ )
  {
      	  Double_t tt = 2.0*atan(exp(-eta[iT])); // These names are short to make the equation manageable
	  Double_t rt = 2.0*atan(exp(-jet->Eta()));
	  Double_t ctt = cos(tt);                   
	  Double_t crt = cos(rt);                         
	  Double_t stt = sin(tt);                               
	  Double_t srt = sin(rt);
	  Double_t ctp = cos(phi[iT]);
	  Double_t crp = cos(jet->Phi());
	  Double_t stp = sin(phi[iT]);
	  Double_t srp = sin(jet->Phi());
	  //Info("plots","acos(%1.16f)\nstt=%f\npt=%f",crp*ctp*srt*stt+srp*stp*srt*stt+crt*ctt,stt,pt[iT]);
	  //Info("plots","diff to 1 %f",1.0-crp*ctp*srt*stt+srp*stp*srt*stt+crt*ctt);
	  Double_t alpha;   
	  if (TMath::Abs(crp*ctp*srt*stt+srp*stp*srt*stt+crt*ctt) > 0.9990)
	  {
		  alpha = 0.0;
	  }else
	  {
		  alpha = TMath::ACos(crp*ctp*srt*stt+srp*stp*srt*stt+crt*ctt);   
	  }
	  Double_t correctp = pt[iT]/stt;
	  fhJetPL->Fill( correctp*cos(alpha),weight);      
	  if ( (jet->Eta()-eta[iT])*(jet->Eta()-eta[iT]) +  
	       (jet->Phi()-phi[iT])*(jet->Phi()-phi[iT]) < 0.2*0.2 )   
		  fhJetJT->Fill( correctp*sin(alpha),weight);
	  fhJetPT->Fill(correctp*sin(tt),weight);
	  if (fNominalEnergy==0.0){
		  fhFragmFcn->Fill(  correctp*sin(tt)/(jet->Energy()*fScaleFactor),weight  ); // This is the jet fragmentation function
	  } else
	  {
                  fhFragmFcn->Fill(  correctp*sin(tt)/fNominalEnergy,weight );
	  }
    	  if (doesJetMeetBinCriteria)
    	  {
	  	  fhRecoBinPt->Fill(correctp*sin(tt),weight);          // ("fhRecoBinPt","Reconstructed Pt Distribution",100,0,1);
	  	  fhRecoBinPtNoBg->Fill(correctp*sin(tt),weight);          // ("fhRecoBinPt","Reconstructed Pt Distribution",100,0,1);
		  fhRecoBinFragmFcn->Fill(  correctp*sin(tt)/(jet->Energy()*fScaleFactor),weight  ); // This is the jet fragmentation function
		  fhRecoBinFragmFcnNoBg->Fill(  correctp*sin(tt)/(jet->Energy()*fScaleFactor),weight  ); // This is the jet fragmentation function
    	  }
  }// loop over tracks
  }

if (numappjet>=1 && fhBackHisto != 0 && doesJetMeetBinCriteria)
  {
    for (Int_t count=1;count<=100;count++)
    {
	fhRecoBinFragmFcnNoBg->Fill( ((Float_t)count)/(jethighest->Energy()*fScaleFactor),-fhBackHisto->GetBinContent(count)*weight);
	fhRecoBinPtNoBg->AddBinContent(count,-fhBackHisto->GetBinContent(count)*weight);
    } 
  }


}


