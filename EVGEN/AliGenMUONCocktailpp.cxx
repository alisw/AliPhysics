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
// Class to create the coktail for physics with muons for pp collisions
// using the The followoing sources: 
// jpsi, psiP, upsilon, upsilonP, upsilonPP, open charm and open beauty
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
//------------------------
// 04/2010:
// - CMS energy passed via parameter
// i.e. gener->SetCMSEnergy(AliGenMUONCocktailpp::kCMS07TeV) in Config.C 
// - resonances now added to the cocktail via AddReso2Generator 
// - cleaning 
// B.Vulpescu & P.Crochet         
 
#include <TObjArray.h>
#include <TParticle.h>
#include <TF1.h>
#include <TVirtualMC.h>
#include <iostream.h>
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

     fCMSEnergyTeV(0),
     fCMSEnergyTeVArray(),
     fSigmaReaction(0),  
     fSigmaReactionArray(),  
     fSigmaJPsi(0),      
     fSigmaJPsiArray(),      
     fSigmaChic1(0),     
     fSigmaChic1Array(),     
     fSigmaChic2(0),     
     fSigmaChic2Array(),     
     fSigmaPsiP(0),      
     fSigmaPsiPArray(),      
     fSigmaUpsilon(0),   
     fSigmaUpsilonArray(),   
     fSigmaUpsilonP(0),  
     fSigmaUpsilonPArray(),  
     fSigmaUpsilonPP(0), 
     fSigmaUpsilonPPArray(), 
     fSigmaCCbar(0),     
     fSigmaCCbarArray(),     
     fSigmaBBbar(0),
     fSigmaBBbarArray(),

     fSigmaSilent(kFALSE)
{
// Constructor

// x-sections for pp @ 7 TeV: charmonia from hep-ph/0311048 Tab.9, page 19,
// bottomnium as for 10 TeV
    fCMSEnergyTeVArray[0] =   7.00;
    fSigmaReactionArray[0] =  0.0695;
    fSigmaJPsiArray[0] =      21.8e-6;
    fSigmaChic1Array[0] =     21.1e-6;
    fSigmaChic2Array[0] =     34.9e-6;
    fSigmaPsiPArray[0] =      4.93e-6;
    fSigmaUpsilonArray[0] =   0.463e-6;
    fSigmaUpsilonPArray[0] =  0.154e-6;
    fSigmaUpsilonPPArray[0] = 0.0886e-6;
    fSigmaCCbarArray[0] =     6.91e-3;
    fSigmaBBbarArray[0] =     0.232e-3;
    
//x-sections for pp @ 10 TeV: charmonia and bottomonia from 14 TeV numbers
// scaled down according to ccbar and bbar cross-sections
    fCMSEnergyTeVArray[1] =   10.00;
    fSigmaReactionArray[1] =  0.0695;
    fSigmaJPsiArray[1] =      26.06e-6;
    fSigmaChic1Array[1] =     25.18e-6;
    fSigmaChic2Array[1] =     41.58e-6;
    fSigmaPsiPArray[1] =      5.88e-6;
    fSigmaUpsilonArray[1] =   0.658e-6;
    fSigmaUpsilonPArray[1] =  0.218e-6;
    fSigmaUpsilonPPArray[1] = 0.122e-6;
    fSigmaCCbarArray[1] =     8.9e-3;
    fSigmaBBbarArray[1] =     0.33e-3;
    
//x-sections for pp @ 14 TeV: charmonia from hep-ph/0311048 Tab.9, page 19,
// bottomonium from hep-ph/0311048 Tab.9, page 19 taken inton account that 
// feed-down from chib is included
    fCMSEnergyTeVArray[2] =   14.00;
    fSigmaReactionArray[2] =  0.070;
    fSigmaJPsiArray[2] =      32.9e-6;
    fSigmaChic1Array[2] =     31.8e-6;
    fSigmaChic2Array[2] =     52.5e-6;
    fSigmaPsiPArray[2] =      7.43e-6;
    fSigmaUpsilonArray[2] =   0.989e-6;
    fSigmaUpsilonPArray[2] =  0.502e-6;
    fSigmaUpsilonPPArray[2] = 0.228e-6;
    fSigmaCCbarArray[2] =     11.2e-3;
    fSigmaBBbarArray[2] =     0.51e-3;
    
}

//_________________________________________________________________________
AliGenMUONCocktailpp::~AliGenMUONCocktailpp()
{
// Destructor

}

//_________________________________________________________________________
void AliGenMUONCocktailpp::SetCMSEnergy(CMSEnergyCode cmsEnergy) 
{
// setter for CMSEnergy and corresponding cross-sections
  fCMSEnergyTeV   = fCMSEnergyTeVArray[cmsEnergy];
  fSigmaReaction  = fSigmaReactionArray[cmsEnergy];  
  fSigmaJPsi      = fSigmaJPsiArray[cmsEnergy];      
  fSigmaChic1     = fSigmaChic1Array[cmsEnergy];     
  fSigmaChic2     = fSigmaChic2Array[cmsEnergy];     
  fSigmaPsiP      = fSigmaPsiPArray[cmsEnergy];      
  fSigmaUpsilon   = fSigmaUpsilonArray[cmsEnergy];   
  fSigmaUpsilonP  = fSigmaUpsilonPArray[cmsEnergy];  
  fSigmaUpsilonPP = fSigmaUpsilonPPArray[cmsEnergy]; 
  fSigmaCCbar     = fSigmaCCbarArray[cmsEnergy];     
  fSigmaBBbar     = fSigmaBBbarArray[cmsEnergy];
}

//_________________________________________________________________________
void AliGenMUONCocktailpp::SetResPolarization(Double_t JpsiPol, Double_t PsiPPol, Double_t UpsPol,
				Double_t UpsPPol, Double_t UpsPPPol, char *PolFrame){
// setter for resonances polarization   
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
// create and add resonances and open HF to the coctail
    Int_t cmsEnergy = Int_t(fCMSEnergyTeV);

// These limits are only used for renormalization of quarkonia cross section
// Pythia events are generated in 4pi  
    Double_t ptMin  = fPtMin;
    Double_t ptMax  = fPtMax;
    Double_t yMin   = fYMin;;
    Double_t yMax   = fYMax;;
    Double_t phiMin = fPhiMin*180./TMath::Pi();
    Double_t phiMax = fPhiMax*180./TMath::Pi();
    AliInfo(Form("Ranges pT:%4.1f : %4.1f GeV/c, y:%4.2f : %4.2f, Phi:%5.1f : %5.1f degres",ptMin,ptMax,yMin,yMax,phiMin,phiMax));
    
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

    if(TMath::Abs(fJpsiPol) > 1.e-30) {sigmajpsi      = fSigmaJPsi*0.0593;}
    if(TMath::Abs(fChic1Pol) > 1.e-30) {sigmachic1     = fSigmaChic1*0.;} // tb consistent
    if(TMath::Abs(fChic2Pol) > 1.e-30) {sigmachic2     = fSigmaChic2*0.;} // tb consistent 
    if(TMath::Abs(fPsiPPol) > 1.e-30) {sigmapsiP      = fSigmaPsiP*0.0075;}
    if(TMath::Abs(fUpsPol) > 1.e-30) {sigmaupsilon   = fSigmaUpsilon*0.0248;}
    if(TMath::Abs(fUpsPPol) > 1.e-30) {sigmaupsilonP  = fSigmaUpsilonP*0.0193;}
    if(TMath::Abs(fUpsPPPol) > 1.e-30) {sigmaupsilonPP = fSigmaUpsilonPP*0.0218;}

    AliInfo(Form("the parametrised resonances uses the decay mode %d",fDecayModeResonance));

// Create and add resonances to the generator
    AliGenParam * genjpsi=0;
    AliGenParam * genchic1=0;
    AliGenParam * genchic2=0;
    AliGenParam * genpsiP=0;
    AliGenParam * genupsilon=0;
    AliGenParam * genupsilonP=0;
    AliGenParam * genupsilonPP=0;
    
    Char_t nameJpsi[10];    
    Char_t nameChic1[10];    
    Char_t nameChic2[10];    
    Char_t namePsiP[10];
    Char_t nameUps[10];    
    Char_t nameUpsP[10];    
    Char_t nameUpsPP[10];    

    sprintf(nameJpsi,"Jpsi");    
    sprintf(nameChic1,"Chic1");
    sprintf(nameChic2,"Chic2");
    sprintf(namePsiP,"PsiP");
    sprintf(nameUps,"Ups");
    sprintf(nameUpsP,"UpsP");
    sprintf(nameUpsPP,"UpsPP");

    if(cmsEnergy == 10){
	genjpsi = new AliGenParam(1, AliGenMUONlib::kJpsi, "CDF pp 10", "Jpsi");
	genchic1 = new AliGenParam(1, AliGenMUONlib::kChic1, "CDF pp 10", "Chic1");
	genchic2 = new AliGenParam(1, AliGenMUONlib::kChic2, "CDF pp 10", "Chic2");
	genpsiP = new AliGenParam(1, AliGenMUONlib::kPsiP, "CDF pp 10", "PsiP");
	genupsilon = new AliGenParam(1, AliGenMUONlib::kUpsilon, "CDF pp 10", "Upsilon");

	genupsilonP = new AliGenParam(1, AliGenMUONlib::kUpsilonP, "CDF pp 10", "UpsilonP");
	genupsilonPP = new AliGenParam(1, AliGenMUONlib::kUpsilonPP, "CDF pp 10", "UpsilonPP");
    } else if (cmsEnergy == 7){
	genjpsi = new AliGenParam(1, AliGenMUONlib::kJpsi, "CDF pp 7", "Jpsi");
	genchic1 = new AliGenParam(1, AliGenMUONlib::kChic1, "CDF pp 7", "Chic1");
	genchic2 = new AliGenParam(1, AliGenMUONlib::kChic2, "CDF pp 7", "Chic2");
	genpsiP = new AliGenParam(1, AliGenMUONlib::kPsiP, "CDF pp 7", "PsiP");

	genupsilon = new AliGenParam(1, AliGenMUONlib::kUpsilon, "CDF pp 7", "Upsilon");
	genupsilonP = new AliGenParam(1, AliGenMUONlib::kUpsilonP, "CDF pp 7", "UpsilonP");
	genupsilonPP = new AliGenParam(1, AliGenMUONlib::kUpsilonPP, "CDF pp 7", "UpsilonPP");
    } else if (cmsEnergy == 14){
	genjpsi = new AliGenParam(1, AliGenMUONlib::kJpsi, "CDF pp ", "Jpsi");
	genchic1 = new AliGenParam(1, AliGenMUONlib::kChic1, "CDF pp ", "Chic1");
	genchic2 = new AliGenParam(1, AliGenMUONlib::kChic2, "CDF pp ", "Chic2");
	genpsiP = new AliGenParam(1, AliGenMUONlib::kPsiP, "CDF pp", "PsiP");

	genupsilon = new AliGenParam(1, AliGenMUONlib::kUpsilon, "CDF pp", "Upsilon");
	genupsilonP = new AliGenParam(1, AliGenMUONlib::kUpsilonP, "CDF pp", "UpsilonP");

	genupsilonPP = new AliGenParam(1, AliGenMUONlib::kUpsilonPP, "CDF pp", "UpsilonPP");	
    }

    AddReso2Generator(nameJpsi,genjpsi,sigmajpsi,fJpsiPol);
    AddReso2Generator(nameChic1,genchic2,sigmachic1,fChic2Pol);
    AddReso2Generator(nameChic2,genpsiP,sigmapsiP,fPsiPPol);    
    AddReso2Generator(namePsiP,genchic1,sigmachic1,fChic1Pol);    
    AddReso2Generator(nameUps,genupsilon,sigmaupsilon,fUpsPol);    
    AddReso2Generator(nameUpsP,genupsilonP,sigmaupsilonP,fUpsPPol);    
    AddReso2Generator(nameUpsPP,genupsilonPP,sigmaupsilonPP,fUpsPPPol);    

//------------------------------------------------------------------
// Generator of charm
    AliGenCorrHF *gencharm = new AliGenCorrHF(1, 4, cmsEnergy);
    gencharm->SetMomentumRange(0,9999);
    gencharm->SetForceDecay(kAll);
    Double_t ratioccbar = sigmaccbar/fSigmaReaction;
    if (!gMC) gencharm->SetDecayer(fDecayer);  
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
    Double_t ratiobbbar = sigmabbbar/fSigmaReaction;
    if (!gMC) genbeauty->SetDecayer(fDecayer);  
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

//-------------------------------------------------------------------
void AliGenMUONCocktailpp::AddReso2Generator(Char_t* nameReso, 
					     AliGenParam* const genReso,
					     Double_t sigmaReso,
					     Double_t polReso)
{
// add resonances to the cocktail
    Double_t phiMin = fPhiMin*180./TMath::Pi();
    Double_t phiMax = fPhiMax*180./TMath::Pi();

    // first step: generation in 4pi
    genReso->SetPtRange(0.,100.); 
    genReso->SetYRange(-8.,8.);
    genReso->SetPhiRange(0.,360.);
    genReso->SetForceDecay(fDecayModeResonance);
    if (!gMC) genReso->SetDecayer(fDecayer);  
    genReso->Init();  // generation in 4pi
// Ratios with respect to the reaction cross-section in the 
// kinematics limit of the MUONCocktail
    Double_t ratioReso = sigmaReso / fSigmaReaction * genReso->GetRelativeArea(fPtMin,fPtMax,fYMin,fYMax,phiMin,phiMax);
    if (!fSigmaSilent) {
      AliInfo(Form("%s prod. cross-section in pp %5.3g b",nameReso,sigmaReso));
      AliInfo(Form("%s prod. probability per collision in acceptance %5.3g",nameReso,ratioReso));
    }
// second step: generation in selected kinematical range
    genReso->SetPtRange(fPtMin, fPtMax);  
    genReso->SetYRange(fYMin, fYMax);
    genReso->SetPhiRange(phiMin, phiMax);
    genReso->Init(); // generation in selected kinematical range

     TVirtualMCDecayer *decReso = 0;
     if(TMath::Abs(polReso) > 1.e-30){
      AliInfo(Form("******Setting polarized decayer for %s''",nameReso));
      if(fPolFrame==0) {
        decReso = new AliDecayerPolarized(polReso,AliDecayerPolarized::kColSop,AliDecayerPolarized::kMuon);
        AliInfo(Form("******Reference frame: %s, alpha: %f","Collins-Soper",polReso));
	  }
      if(fPolFrame==1) {
        decReso = new AliDecayerPolarized(polReso,AliDecayerPolarized::kHelicity,AliDecayerPolarized::kMuon);
        AliInfo(Form("******Reference frame: %s, alpha: %f","Helicity",polReso));
	  }
      decReso->SetForceDecay(kAll);
      decReso->Init();
      genReso->SetDecayer(decReso);
     }
    
    AddGenerator(genReso,nameReso,ratioReso); // Adding Generator    
}

//-------------------------------------------------------------------
void AliGenMUONCocktailpp::Init()
{
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

    
