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

//
// Classe to create the coktail for physics with muons for pp at 14 TeV
// The followoing sources: 
// jpsi,psiP, upsilon, upsilonP, upsilonPP are added to Pythia
// The free parameeters are :
//      pp reaction cross-section
//      production cross-sections in pp collisions 
// July 07:added heavy quark production from AliGenCorrHF and heavy quark 
//         production switched off in Pythia 
// Aug. 07: added trigger cut on total momentum
// 2009: added possibility to hide x-sections (B. Vulpescu)
// 2009: added possibility to have the cocktail (fast generator and param.) 
// for pp @ 10 TeV or pp @ 14 TeV (N. Bastid)
//-----------------
// 2009:  added polarization (L. Bianchi)
//------------------
// 11/2009: added chi_c1 & chi_c2 (P.Crochet & N.Bastid). 
// Cross-sections for charmonia are now directly taken from the Yellow Report 
// (hep-ph/0311048) Tab.9, page 19. See below for details w.r.t. beam energy. 
// usage: see example of Config in $ALICE_ROOT/prod/LHC09a10/Config.C

#include <TObjArray.h>
#include <TParticle.h>
#include <TF1.h>
#include <TVirtualMC.h>

#include "AliGenCocktailEventHeader.h"

#include "AliGenCocktailEntry.h"
#include "AliGenMUONCocktailpp.h"
#include "AliGenMUONlib.h"
#include "AliGenParam.h"
#include "AliMC.h"
#include "AliRun.h"
#include "AliStack.h"
#include "AliDecayer.h"
#include "AliLog.h"
#include "AliGenCorrHF.h"
#include "AliDecayerPolarized.h"

ClassImp(AliGenMUONCocktailpp)  
  
//________________________________________________________________________
AliGenMUONCocktailpp::AliGenMUONCocktailpp()
    :AliGenCocktail(),
     fDecayer(0),
     fDecayModeResonance(kAll),
     fDecayModePythia(kAll),
     fMuonMultiplicity(0),
     fMuonPtCut(0.),
     fMuonPCut(0.),
     fMuonThetaMinCut(0.),
     fMuonThetaMaxCut(180.),
     fMuonOriginCut(-999.),
     fNSucceded(0),
     fNGenerated(0),

     fJpsiPol(0), 
     fChic1Pol(0), 
     fChic2Pol(0), 
     fPsiPPol(0), 
     fUpsPol(0), 
     fUpsPPol(0), 
     fUpsPPPol(0),
     fPolFrame(0),

// x-sections for pp @ 7 TeV: charmonia from hep-ph/0311048 Tab.9, page 19,
// bottomnium as for 10 TeV
/*     fCMSEnergy(7),
     fSigmaReaction(0.0695),
     fSigmaJPsi(21.8e-6),
     fSigmaChic1(21.1e-6),
     fSigmaChic2(34.9e-6),
     fSigmaPsiP(4.93e-6),
     fSigmaUpsilon(0.463e-6),
     fSigmaUpsilonP(0.154e-6),
     fSigmaUpsilonPP(0.0886e-6),
     fSigmaCCbar(6.91e-3),
     fSigmaBBbar(0.232e-3),
*/

//x-sections for pp @ 10 TeV: charmonia and bottomonia from 14 TeV numbers
// scaled down according to ccbar and bbar cross-sections
     fCMSEnergy(10),
     fSigmaReaction(0.0695),
     fSigmaJPsi(26.06e-6),
     fSigmaChic1(25.18e-6),
     fSigmaChic2(41.58e-6),
     fSigmaPsiP(5.88e-6),
     fSigmaUpsilon(0.658e-6),
     fSigmaUpsilonP(0.218e-6),
     fSigmaUpsilonPP(0.122e-6),
     fSigmaCCbar(8.9e-3),
     fSigmaBBbar(0.33e-3),


//x-sections for pp @ 14 TeV: charmonia from hep-ph/0311048 Tab.9, page 19,
// bottomonium from hep-ph/0311048 Tab.9, page 19 taken inton account that 
// feed-down from chib is included
/*     fCMSEnergy(14),
     fSigmaReaction(0.070),
     fSigmaJPsi(32.9e-6),
     fSigmaChic1(31.8e-6),
     fSigmaChic2(52.5e-6),
     fSigmaPsiP(7.43e-6),
     fSigmaUpsilon(0.989e-6),
     fSigmaUpsilonP(0.502e-6),
     fSigmaUpsilonPP(0.228e-6),
     fSigmaCCbar(11.2e-3),
     fSigmaBBbar(0.51e-3),
*/
     fSigmaSilent(kFALSE)
{
// Constructor
}

//_________________________________________________________________________
AliGenMUONCocktailpp::~AliGenMUONCocktailpp()
// Destructor
{
    
}

//_________________________________________________________________________
void AliGenMUONCocktailpp::SetResPolarization(Double_t JpsiPol, Double_t PsiPPol, Double_t UpsPol,
				Double_t UpsPPol, Double_t UpsPPPol, char *PolFrame){
   
   if (strcmp(PolFrame,"kColSop")==0){
      fJpsiPol  = (JpsiPol>=-1  && JpsiPol<=1)  ? JpsiPol  : 0;
      fPsiPPol  = (PsiPPol>=-1  && PsiPPol<=1)  ? PsiPPol  : 0;
      fUpsPol   = (UpsPol>=-1   && UpsPol<=1)   ? UpsPol   : 0;
      fUpsPPol  = (UpsPPol>=-1  && UpsPPol<=1)  ? UpsPPol  : 0;
      fUpsPPPol = (UpsPPPol>=-1 && UpsPPPol<=1) ? UpsPPPol : 0;
      fPolFrame = 0;
   } else if (strcmp(PolFrame,"kHelicity")==0){
        fJpsiPol  = (JpsiPol>=-1  && JpsiPol<=1)  ? JpsiPol  : 0;
        fPsiPPol  = (PsiPPol>=-1  && PsiPPol<=1)  ? PsiPPol  : 0;
        fUpsPol   = (UpsPol>=-1   && UpsPol<=1)   ? UpsPol   : 0;
        fUpsPPol  = (UpsPPol>=-1  && UpsPPol<=1)  ? UpsPPol  : 0;
        fUpsPPPol = (UpsPPPol>=-1 && UpsPPPol<=1) ? UpsPPPol : 0;
	fPolFrame = 1;
   
   } else {
       AliInfo(Form("The polarization frame is not valid"));
       AliInfo(Form("No polarization will be set"));
       fJpsiPol=0.;
       fPsiPPol=0.;
       fUpsPol=0.;
       fUpsPPol=0.;
       fUpsPPPol=0.;
   }
}

//_________________________________________________________________________
void AliGenMUONCocktailpp::CreateCocktail()
{
    
// These limits are only used for renormalization of quarkonia cross section
// Pythia events are generated in 4pi  
    Double_t ptMin  = fPtMin;
    Double_t ptMax  = fPtMax;
    Double_t yMin   = fYMin;;
    Double_t yMax   = fYMax;;
    Double_t phiMin = fPhiMin*180./TMath::Pi();
    Double_t phiMax = fPhiMax*180./TMath::Pi();
    AliInfo(Form("Ranges pT:%4.1f : %4.1f GeV/c, y:%4.2f : %4.2f, Phi:%5.1f : %5.1f degres",ptMin,ptMax,yMin,yMax,phiMin,phiMax));
    
// Ratios with respect to the reaction cross-section in the 
// kinematics limit of the MUONCocktail
    Double_t ratiojpsi;
    Double_t ratiochic1;
    Double_t ratiochic2;
    Double_t ratiopsiP;
    Double_t ratioupsilon;
    Double_t ratioupsilonP;
    Double_t ratioupsilonPP;
    Double_t ratioccbar;
    Double_t ratiobbbar;

// Beam energy
    Int_t cmsEnergy = fCMSEnergy; 

// Cross sections in barns (from PPR Vol. II p: 552) pp - 14 TeV and 
// corrected from feed down of higher resonances 

    Double_t sigmajpsi      = fSigmaJPsi;  
    Double_t sigmachic1     = fSigmaChic1;  
    Double_t sigmachic2     = fSigmaChic2;  
    Double_t sigmapsiP      = fSigmaPsiP;  
    Double_t sigmaupsilon   = fSigmaUpsilon;  
    Double_t sigmaupsilonP  = fSigmaUpsilonP;  
    Double_t sigmaupsilonPP = fSigmaUpsilonPP;
    Double_t sigmaccbar     = fSigmaCCbar;
    Double_t sigmabbbar     = fSigmaBBbar;

// Cross sections corrected with the BR in mu+mu-
// (only in case of use of AliDecayerPolarized)

    if(fJpsiPol  != 0)	{sigmajpsi      = fSigmaJPsi*0.0593;}
    if(fChic1Pol != 0)	{sigmachic1     = fSigmaChic1*0.;} // tb consistent
    if(fChic2Pol != 0)	{sigmachic2     = fSigmaChic2*0.;} // tb consistent 
    if(fPsiPPol  != 0)	{sigmapsiP      = fSigmaPsiP*0.0075;}
    if(fUpsPol   != 0)	{sigmaupsilon   = fSigmaUpsilon*0.0248;}
    if(fUpsPPol  != 0)	{sigmaupsilonP  = fSigmaUpsilonP*0.0193;}
    if(fUpsPPPol != 0)	{sigmaupsilonPP = fSigmaUpsilonPP*0.0218;}

// MinBias NN cross section @ pp 14 TeV -PR  Vol II p:473
    Double_t sigmaReaction = fSigmaReaction;

    Int_t eincStart = 10;

    AliInfo(Form("the parametrised resonances uses the decay mode %d",fDecayModeResonance));

// Generation using CDF scaled distribution for pp@14TeV (from D. Stocco)
//----------------------------------------------------------------------
// Jpsi generator
    AliGenParam * genjpsi;
    if(cmsEnergy == eincStart){
	genjpsi = new AliGenParam(1, AliGenMUONlib::kJpsi, "CDF pp 10", "Jpsi");
    } else {
	genjpsi = new AliGenParam(1, AliGenMUONlib::kJpsi, "CDF pp ", "Jpsi");
    }
// first step: generation in 4pi
    genjpsi->SetPtRange(0.,100.);
    genjpsi->SetYRange(-8.,8.);
    genjpsi->SetPhiRange(0.,360.);
    genjpsi->SetForceDecay(fDecayModeResonance);
    //if (!gMC) genjpsi->SetDecayer(fDecayer);

    genjpsi->Init();  // generation in 4pi
    ratiojpsi = sigmajpsi / sigmaReaction * genjpsi->GetRelativeArea(ptMin,ptMax,yMin,yMax,phiMin,phiMax); // get weight
    if (!fSigmaSilent) {
      AliInfo(Form("jpsi prod. cross-section in pp %5.3g b",sigmajpsi));
      AliInfo(Form("jpsi prod. probability per collision in acceptance %5.3g",ratiojpsi));
    }
// second step: generation in selected kinematical range
    genjpsi->SetPtRange(ptMin, ptMax);  
    genjpsi->SetYRange(yMin, yMax);
    genjpsi->SetPhiRange(phiMin, phiMax);
    genjpsi->Init(); // generation in selected kinematic range

    TVirtualMCDecayer *Jpsidec = 0;
    if(fJpsiPol != 0){
    AliInfo(Form("******Setting polarized decayer for J/psi"));
    if(fPolFrame==0) {
         Jpsidec = new AliDecayerPolarized(fJpsiPol,AliDecayerPolarized::kColSop,AliDecayerPolarized::kMuon);
	 AliInfo(Form("******Reference frame: %s, alpha: %f","Collins-Soper",fJpsiPol));
	   }
    if(fPolFrame==1) {
         Jpsidec = new AliDecayerPolarized(fJpsiPol,AliDecayerPolarized::kHelicity,AliDecayerPolarized::kMuon);
         AliInfo(Form("******Reference frame: %s, alpha: %f","Helicity",fJpsiPol));
	   }
    Jpsidec->SetForceDecay(kAll);
    Jpsidec->Init();
    genjpsi->SetDecayer(Jpsidec);
    }

    AddGenerator(genjpsi, "Jpsi", ratiojpsi); // Adding Generator


//----------------------------------------------------------------------
// Chic1 generator
    AliGenParam * genchic1;
    if(cmsEnergy == eincStart){
	genchic1 = new AliGenParam(1, AliGenMUONlib::kChic1, "CDF pp 10", "Chic1");
    } else {
	genchic1 = new AliGenParam(1, AliGenMUONlib::kChic1, "CDF pp ", "Chic1");
    }
// first step: generation in 4pi
    genchic1->SetPtRange(0.,100.);
    genchic1->SetYRange(-8.,8.);
    genchic1->SetPhiRange(0.,360.);
    genchic1->SetForceDecay(fDecayModeResonance);
    //if (!gMC) genchic1->SetDecayer(fDecayer);

    genchic1->Init();  // generation in 4pi
    ratiochic1 = sigmachic1 / sigmaReaction * genchic1->GetRelativeArea(ptMin,ptMax,yMin,yMax,phiMin,phiMax); // get weight
    if (!fSigmaSilent) {
      AliInfo(Form("chic1 prod. cross-section in pp %5.3g b",sigmachic1));
      AliInfo(Form("chic1 prod. probability per collision in acceptance %5.3g",ratiochic1));
    }
// second step: generation in selected kinematical range
    genchic1->SetPtRange(ptMin, ptMax);  
    genchic1->SetYRange(yMin, yMax);
    genchic1->SetPhiRange(phiMin, phiMax);
    genchic1->Init(); // generation in selected kinematic range

    TVirtualMCDecayer *Chic1dec = 0;
    if(fChic1Pol != 0){
    AliInfo(Form("******Setting polarized decayer for Chic1"));
    if(fPolFrame==0) {
         Chic1dec = new AliDecayerPolarized(fChic1Pol,AliDecayerPolarized::kColSop,AliDecayerPolarized::kMuon);
	 AliInfo(Form("******Reference frame: %s, alpha: %f","Collins-Soper",fChic1Pol));
	   }
    if(fPolFrame==1) {
         Chic1dec = new AliDecayerPolarized(fChic1Pol,AliDecayerPolarized::kHelicity,AliDecayerPolarized::kMuon);
         AliInfo(Form("******Reference frame: %s, alpha: %f","Helicity",fChic1Pol));
	   }
    Chic1dec->SetForceDecay(kAll);
    Chic1dec->Init();
    genchic1->SetDecayer(Chic1dec);
    }

    AddGenerator(genchic1, "Chic1", ratiochic1); // Adding Generator

//----------------------------------------------------------------------
// Chic2 generator
    AliGenParam * genchic2;
    if(cmsEnergy == eincStart){
	genchic2 = new AliGenParam(1, AliGenMUONlib::kChic2, "CDF pp 10", "Chic2");
    } else {
	genchic2 = new AliGenParam(1, AliGenMUONlib::kChic2, "CDF pp ", "Chic2");
    }
// first step: generation in 4pi
    genchic2->SetPtRange(0.,100.);
    genchic2->SetYRange(-8.,8.);
    genchic2->SetPhiRange(0.,360.);
    genchic2->SetForceDecay(fDecayModeResonance);
    //if (!gMC) genchic1->SetDecayer(fDecayer);

    genchic2->Init();  // generation in 4pi
    ratiochic2 = sigmachic2 / sigmaReaction * genchic2->GetRelativeArea(ptMin,ptMax,yMin,yMax,phiMin,phiMax); // get weight
    if (!fSigmaSilent) {
      AliInfo(Form("chic2 prod. cross-section in pp %5.3g b",sigmachic2));
      AliInfo(Form("chic2 prod. probability per collision in acceptance %5.3g",ratiochic2));
    }
// second step: generation in selected kinematical range
    genchic2->SetPtRange(ptMin, ptMax);  
    genchic2->SetYRange(yMin, yMax);
    genchic2->SetPhiRange(phiMin, phiMax);
    genchic2->Init(); // generation in selected kinematic range

    TVirtualMCDecayer *Chic2dec = 0;
    if(fChic2Pol != 0){
    AliInfo(Form("******Setting polarized decayer for Chic2"));
    if(fPolFrame==0) {
         Chic2dec = new AliDecayerPolarized(fChic2Pol,AliDecayerPolarized::kColSop,AliDecayerPolarized::kMuon);
	 AliInfo(Form("******Reference frame: %s, alpha: %f","Collins-Soper",fChic2Pol));
	   }
    if(fPolFrame==1) {
         Chic2dec = new AliDecayerPolarized(fChic2Pol,AliDecayerPolarized::kHelicity,AliDecayerPolarized::kMuon);
         AliInfo(Form("******Reference frame: %s, alpha: %f","Helicity",fChic2Pol));
	   }
    Chic2dec->SetForceDecay(kAll);
    Chic2dec->Init();
    genchic2->SetDecayer(Chic2dec);
    }

    AddGenerator(genchic2, "Chic2", ratiochic2); // Adding Generator

//------------------------------------------------------------------
// Psi prime generator
    AliGenParam * genpsiP;
    if(cmsEnergy == eincStart){    
	genpsiP = new AliGenParam(1, AliGenMUONlib::kPsiP, "CDF pp 10", "PsiP");
    } else {
	genpsiP = new AliGenParam(1, AliGenMUONlib::kPsiP, "CDF pp", "PsiP");
    }
// first step: generation in 4pi
    genpsiP->SetPtRange(0.,100.);
    genpsiP->SetYRange(-8.,8.);
    genpsiP->SetPhiRange(0.,360.);
    genpsiP->SetForceDecay(fDecayModeResonance);
    //if (!gMC) genpsiP->SetDecayer(fDecayer);

    genpsiP->Init();  // generation in 4pi
    ratiopsiP = sigmapsiP / sigmaReaction * genpsiP->GetRelativeArea(ptMin,ptMax,yMin,yMax,phiMin,phiMax);
    if (!fSigmaSilent) {
      AliInfo(Form("psiP prod. cross-section in pp %5.3g b",sigmapsiP));
      AliInfo(Form("psiP prod. probability per collision in acceptance %5.3g",ratiopsiP));
    }
// second step: generation in selected kinematical range
    genpsiP->SetPtRange(ptMin, ptMax);  
    genpsiP->SetYRange(yMin, yMax);
    genpsiP->SetPhiRange(phiMin, phiMax);
    genpsiP->Init(); // generation in selected kinematic range

     TVirtualMCDecayer *PsiPdec = 0;
     if(fPsiPPol != 0){
      AliInfo(Form("******Setting polarized decayer for psi'"));
      if(fPolFrame==0) {
        PsiPdec = new AliDecayerPolarized(fPsiPPol,AliDecayerPolarized::kColSop,AliDecayerPolarized::kMuon);
        AliInfo(Form("******Reference frame: %s, alpha: %f","Collins-Soper",fPsiPPol));
	  }
      if(fPolFrame==1) {
        PsiPdec = new AliDecayerPolarized(fPsiPPol,AliDecayerPolarized::kHelicity,AliDecayerPolarized::kMuon);
        AliInfo(Form("******Reference frame: %s, alpha: %f","Helicity",fPsiPPol));
	  }
      PsiPdec->SetForceDecay(kAll);
      PsiPdec->Init();
      genpsiP->SetDecayer(PsiPdec);
     }

    AddGenerator(genpsiP, "PsiP", ratiopsiP); // Adding Generator

//------------------------------------------------------------------
// Upsilon 1S generator
    AliGenParam * genupsilon;
    if(cmsEnergy == eincStart) {
	genupsilon = new AliGenParam(1, AliGenMUONlib::kUpsilon, "CDF pp 10", "Upsilon");
    } else {
	genupsilon = new AliGenParam(1, AliGenMUONlib::kUpsilon, "CDF pp", "Upsilon");
    }
// first step: generation in 4pi
    genupsilon->SetPtRange(0.,100.); 
    genupsilon->SetYRange(-8.,8.);
    genupsilon->SetPhiRange(0.,360.);
    genupsilon->SetForceDecay(fDecayModeResonance);
    //if (!gMC) genupsilon->SetDecayer(fDecayer);
    genupsilon->Init();  // generation in 4pi
    ratioupsilon = sigmaupsilon / sigmaReaction * genupsilon->GetRelativeArea(ptMin,ptMax,yMin,yMax,phiMin,phiMax);
    if (!fSigmaSilent) {
      AliInfo(Form("Upsilon 1S prod. cross-section in pp %5.3g b",sigmaupsilon));
      AliInfo(Form("Upsilon 1S prod. probability per collision in acceptance %5.3g",ratioupsilon));
    }
// second step: generation in selected kinematical range
    genupsilon->SetPtRange(ptMin, ptMax);  
    genupsilon->SetYRange(yMin, yMax);
    genupsilon->SetPhiRange(phiMin, phiMax);
    genupsilon->Init(); // generation in selected kinematical range
     
     TVirtualMCDecayer *Upsdec = 0;
     if(fUpsPol != 0){
      AliInfo(Form("******Setting polarized decayer for Upsilon"));
      if(fPolFrame==0) {
        Upsdec = new AliDecayerPolarized(fUpsPol,AliDecayerPolarized::kColSop,AliDecayerPolarized::kMuon);
        AliInfo(Form("******Reference frame: %s, alpha: %f","Collins-Soper",fUpsPol));
	  }
      if(fPolFrame==1) {
        Upsdec = new AliDecayerPolarized(fUpsPol,AliDecayerPolarized::kHelicity,AliDecayerPolarized::kMuon);
        AliInfo(Form("******Reference frame: %s, alpha: %f","Helicity",fUpsPol));
	  }
      Upsdec->SetForceDecay(kAll);
      Upsdec->Init();
      genupsilon->SetDecayer(Upsdec);
     }

    AddGenerator(genupsilon,"Upsilon", ratioupsilon); // Adding Generator
//------------------------------------------------------------------
// Upsilon 2S generator
    AliGenParam * genupsilonP;
    if(cmsEnergy == eincStart){
	genupsilonP = new AliGenParam(1, AliGenMUONlib::kUpsilonP, "CDF pp 10", "UpsilonP");
    } else {
	genupsilonP = new AliGenParam(1, AliGenMUONlib::kUpsilonP, "CDF pp", "UpsilonP");
    }
	
// first step: generation in 4pi
    genupsilonP->SetPtRange(0.,100.);
    genupsilonP->SetYRange(-8.,8.);
    genupsilonP->SetPhiRange(0.,360.);
    genupsilonP->SetForceDecay(fDecayModeResonance);
    //if (!gMC) genupsilonP->SetDecayer(fDecayer);  
    genupsilonP->Init();  // generation in 4pi
    ratioupsilonP = sigmaupsilonP / sigmaReaction * genupsilonP->GetRelativeArea(ptMin,ptMax,yMin,yMax,phiMin,phiMax);
    if (!fSigmaSilent) {
      AliInfo(Form("Upsilon 2S prod. cross-section in pp %5.3g b",sigmaupsilonP));
      AliInfo(Form("Upsilon 2S prod. probability per collision in acceptance %5.3g",ratioupsilonP));
    }
// second step: generation in the kinematical range
    genupsilonP->SetPtRange(ptMin, ptMax);  
    genupsilonP->SetYRange(yMin, yMax);
    genupsilonP->SetPhiRange(phiMin, phiMax);
    genupsilonP->Init(); // generation in selected kinematical range

     TVirtualMCDecayer *UpsPdec = 0;
     if(fUpsPPol != 0){
      AliInfo(Form("******Setting polarized decayer for Upsilon'"));
      if(fPolFrame==0) {
        UpsPdec = new AliDecayerPolarized(fUpsPPol,AliDecayerPolarized::kColSop,AliDecayerPolarized::kMuon);
        AliInfo(Form("******Reference frame: %s, alpha: %f","Collins-Soper",fUpsPPol));
	  }
      if(fPolFrame==1) {
        UpsPdec = new AliDecayerPolarized(fUpsPPol,AliDecayerPolarized::kHelicity,AliDecayerPolarized::kMuon);
        AliInfo(Form("******Reference frame: %s, alpha: %f","Helicity",fUpsPPol));
	  }
      UpsPdec->SetForceDecay(kAll);
      UpsPdec->Init();
      genupsilonP->SetDecayer(UpsPdec);
     }
    
    AddGenerator(genupsilonP,"UpsilonP", ratioupsilonP); // Adding Generator
    
//------------------------------------------------------------------
// Upsilon 3S generator
    AliGenParam * genupsilonPP;
    if(cmsEnergy == eincStart){
	genupsilonPP = new AliGenParam(1, AliGenMUONlib::kUpsilonPP, "CDF pp 10", "UpsilonPP");
    }
    else {
	genupsilonPP = new AliGenParam(1, AliGenMUONlib::kUpsilonPP, "CDF pp", "UpsilonPP");	
    }

// first step: generation in 4pi
    genupsilonPP->SetPtRange(0.,100.); 
    genupsilonPP->SetYRange(-8.,8.);
    genupsilonPP->SetPhiRange(0.,360.);
    genupsilonPP->SetForceDecay(fDecayModeResonance);
    //if (!gMC) genupsilonPP->SetDecayer(fDecayer);  
    genupsilonPP->Init();  // generation in 4pi
    ratioupsilonPP = sigmaupsilonPP / sigmaReaction * genupsilonPP->GetRelativeArea(ptMin,ptMax,yMin,yMax,phiMin,phiMax);
    if (!fSigmaSilent) {
      AliInfo(Form("Upsilon 3S prod. cross-section in pp %5.3g b",sigmaupsilonPP));
      AliInfo(Form("Upsilon 3S prod. probability per collision in acceptance %5.3g",ratioupsilonPP));
    }
// second step: generation in selected kinematical range
    genupsilonPP->SetPtRange(ptMin, ptMax);  
    genupsilonPP->SetYRange(yMin, yMax);
    genupsilonPP->SetPhiRange(phiMin, phiMax);
    genupsilonPP->Init(); // generation in selected kinematical range

     TVirtualMCDecayer *UpsPPdec = 0;
     if(fUpsPPPol != 0){
      AliInfo(Form("******Setting polarized decayer for Upsilon''"));
      if(fPolFrame==0) {
        UpsPPdec = new AliDecayerPolarized(fUpsPPPol,AliDecayerPolarized::kColSop,AliDecayerPolarized::kMuon);
        AliInfo(Form("******Reference frame: %s, alpha: %f","Collins-Soper",fUpsPPPol));
	  }
      if(fPolFrame==1) {
        UpsPPdec = new AliDecayerPolarized(fUpsPPPol,AliDecayerPolarized::kHelicity,AliDecayerPolarized::kMuon);
        AliInfo(Form("******Reference frame: %s, alpha: %f","Helicity",fUpsPPPol));
	  }
      UpsPPdec->SetForceDecay(kAll);
      UpsPPdec->Init();
      genupsilonPP->SetDecayer(UpsPPdec);
     }
    
    AddGenerator(genupsilonPP,"UpsilonPP", ratioupsilonPP); // Adding Generator

//------------------------------------------------------------------
// Generator of charm
    
    AliGenCorrHF *gencharm = new AliGenCorrHF(1, 4, cmsEnergy);  
    gencharm->SetMomentumRange(0,9999);
    gencharm->SetForceDecay(kAll);
    ratioccbar = sigmaccbar/sigmaReaction;
    //if (!gMC) gencharm->SetDecayer(fDecayer);  
    gencharm->Init();
    if (!fSigmaSilent) {
      AliInfo(Form("c-cbar prod. cross-section in pp %5.3g b",sigmaccbar));
      AliInfo(Form("c-cbar prod. probability per collision in acceptance %5.3g",ratioccbar));
    }
    AddGenerator(gencharm,"CorrHFCharm",ratioccbar);

//------------------------------------------------------------------
// Generator of beauty

    AliGenCorrHF *genbeauty = new AliGenCorrHF(1, 5, cmsEnergy);  
    genbeauty->SetMomentumRange(0,9999);
    genbeauty->SetForceDecay(kAll);
    ratiobbbar = sigmabbbar/sigmaReaction;
    //if (!gMC) genbeauty->SetDecayer(fDecayer);  
    genbeauty->Init();
    if (!fSigmaSilent) {
      AliInfo(Form("b-bbar prod. cross-section in pp  %5.3g b",sigmabbbar));
      AliInfo(Form("b-bbar prod. probability per collision in acceptance %5.3g",ratiobbbar));
    }
    AddGenerator(genbeauty,"CorrHFBeauty",ratiobbbar); 

//-------------------------------------------------------------------
// Pythia generator
//
// This has to go into the Config.C
//
//    AliGenPythia *pythia = new AliGenPythia(1);
//    pythia->SetProcess(kPyMbMSEL1);
//    pythia->SetStrucFunc(kCTEQ5L);
//    pythia->SetEnergyCMS(14000.);
//    AliInfo(Form("\n\npythia uses the decay mode %d", GetDecayModePythia()));
//    Decay_t dt = gener->GetDecayModePythia();
//    pythia->SetForceDecay(dt);
//    pythia->SetPtRange(0.,100.);
//    pythia->SetYRange(-8.,8.);
//    pythia->SetPhiRange(0.,360.);
//    pythia->SetPtHard(2.76,-1.0);
//    pythia->SwitchHFOff();
//    pythia->Init(); 
//    AddGenerator(pythia,"Pythia",1);

}

void AliGenMUONCocktailpp::Init()
{
    //
    // Initialisation
    TIter next(fEntries);
    AliGenCocktailEntry *entry;
    if (fStack) {
	while((entry = (AliGenCocktailEntry*)next())) {
	    entry->Generator()->SetStack(fStack);
	}
    }
}

//_________________________________________________________________________
void AliGenMUONCocktailpp::Generate()
{
// Generate event 
    TIter next(fEntries);
    AliGenCocktailEntry *entry = 0;
    AliGenCocktailEntry *preventry = 0;
    AliGenerator* gen = 0;

    if (fHeader) delete fHeader;
    fHeader = new AliGenCocktailEventHeader("MUON Cocktail Header");

    const TObjArray *partArray = gAlice->GetMCApp()->Particles();
    
// Generate the vertex position used by all generators    
    if(fVertexSmear == kPerEvent) Vertex();

// Loop on primordialTrigger: 
// minimum muon multiplicity above a pt cut in a theta acceptance region

    Bool_t primordialTrigger = kFALSE;
    while(!primordialTrigger) {
	//Reseting stack
	AliRunLoader * runloader = AliRunLoader::Instance();
	if (runloader)
	    if (runloader->Stack())
		runloader->Stack()->Clean();
	// Loop over generators and generate events
	Int_t igen = 0;
	Int_t npart = 0;	
	const char* genName = 0;
	while((entry = (AliGenCocktailEntry*)next())) {
	    gen = entry->Generator();
	    genName = entry->GetName();
	    gen->SetVertex(fVertex.At(0), fVertex.At(1), fVertex.At(2));

	    npart = (strcmp(genName,"Pythia") == 0) ? 1 :
		gRandom->Poisson(entry->Rate());

	    if (npart > 0) {
		igen++;	
		if (igen == 1) entry->SetFirst(0);		
		else  entry->SetFirst((partArray->GetEntriesFast())+1);
		
		gen->SetNumberParticles(npart);		
		gen->Generate();
		entry->SetLast(partArray->GetEntriesFast());
		preventry = entry;
	    }
	}  
	next.Reset();

// Testing primordial trigger: Single muons or dimuons with Pt above a Pt cut 
// in the muon spectrometer acceptance
	Int_t iPart;
	fNGenerated++;
	Int_t numberOfMuons=0;Int_t maxPart = partArray->GetEntriesFast();
	for(iPart=0; iPart<maxPart; iPart++){      
	    
	  TParticle *part = gAlice->GetMCApp()->Particle(iPart);
	  if ( TMath::Abs(part->GetPdgCode()) == 13 ){  
	    if((part->Vz() > fMuonOriginCut) && //take only the muons that decayed before the abs + 1 int. length in C abs
	       (part->Theta()*180./TMath::Pi()>fMuonThetaMinCut) &&
	       (part->Theta()*180./TMath::Pi()<fMuonThetaMaxCut) &&
	       (part->Pt()>fMuonPtCut) &&
	       (part->P()>fMuonPCut)) {
	      numberOfMuons++;
	    }
	  }
	}	
	if (numberOfMuons >= fMuonMultiplicity) {
	  primordialTrigger = kTRUE;
	  fHeader->SetNProduced(maxPart);
	}

    }
    fNSucceded++;

    TArrayF eventVertex;
    eventVertex.Set(3);
    for (Int_t j=0; j < 3; j++) eventVertex[j] = fVertex[j];

    fHeader->SetPrimaryVertex(eventVertex);

    gAlice->SetGenEventHeader(fHeader);

//     AliInfo(Form("Generated Events are %d and Succeeded Events are %d",fNGenerated,fNSucceded));
    AliDebug(5,Form("Generated Events are %d and Succeeded Events are %d",fNGenerated,fNSucceded));
}

    
