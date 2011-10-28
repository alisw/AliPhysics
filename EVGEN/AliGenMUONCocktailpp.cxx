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
//-----------------------         
// 10/2011: 
// - added the cocktail for p-Pb & Pb-p @ 8.8 TeV with 4 centrality bins and
// for Pb-Pb @ 2.76 TeV with 11 centrality bins. Bins should be defined also
// in the Config.C with one AliGenMUONCocktailpp per bin. These generators
// included in a AliGenCocktail together with an event generator (e.g. Hijing)
// providing the underlying event and collision centrality. The bin number n
// passed via AliGenMUONCocktailpp::SetCentralityBin(n).
// See details in my presentation at the PWG3-Muon meeting (05.10.2011):
// https://indico.cern.ch/conferenceDisplay.py?confId=157367
// - simplifications and bug fix in CreateCocktail() 
// S. Grigoryan
 
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
     fCentralityBin(0),

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

// x-sections for pp @ 7 TeV: 
// -charmonia: 4pi integral of fit function for inclusive J/psi dsigma/dy LHC data 
// gives 60 mub; so sigma_prompt = 54 mub, while Ref = R.Vogt_arXiv:1003.3497 (Table 2)
// gives 35 mub. Below we use sigma_direct from the Ref scaled by the factor 54/35.
// -bottomonia: 4pi integral of fit function for inclusive Upsilon1S dsigma/dy LHC data
// gives 0.56 mub, sigmas for 2S & 3S obtained using CMS data for ratios 2S/1S & 3S/1S
// -ccbar & bbbar: NLO pQCD computations - http://www-alice.gsi.de/ana/MNR/results.html
    fCMSEnergyTeVArray[0] =   7.00;
    fSigmaReactionArray[0] =  0.070;
    fSigmaJPsiArray[0] =      33.6e-6;
    fSigmaChic1Array[0] =     32.6e-6;
    fSigmaChic2Array[0] =     53.8e-6;
    fSigmaPsiPArray[0] =       7.6e-6;
    fSigmaUpsilonArray[0] =   0.56e-6;
    fSigmaUpsilonPArray[0] =  0.19e-6;
    fSigmaUpsilonPPArray[0] = 0.09e-6;
    fSigmaCCbarArray[0] =     6.91e-3;
    fSigmaBBbarArray[0] =     0.232e-3;
    
//x-sections for pp @ 10 TeV: charmonia and bottomonia from 14 TeV numbers
// scaled down according to ccbar and bbbar cross-sections
    fCMSEnergyTeVArray[1] =   10.00;
    fSigmaReactionArray[1] =  0.070;
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
// bottomonium from hep-ph/0311048 Tab.9, page 19 taken into account that 
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
    fSigmaBBbarArray[2] =     0.445e-3;

// x-sections for Min. Bias p-Pb & Pb-p @ 8.8 TeV: charmonia and bottomonia 
// from 7 TeV numbers scaled according to pQCD ccbar and bbbar x-sections
// and with Glauber scaling
    fCMSEnergyTeVArray[3] =   9.00;           // for 8.8 TeV
    fSigmaReactionArray[3] =  2.10;
    fSigmaJPsiArray[3] =      8.19e-3;        // 208*1.172*33.6e-6
    fSigmaChic1Array[3] =     7.95e-3;
    fSigmaChic2Array[3] =     13.1e-3;
    fSigmaPsiPArray[3] =      1.85e-3;
    fSigmaUpsilonArray[3] =   0.146e-3;       // 208*1.25*0.56e-6
    fSigmaUpsilonPArray[3] =  0.049e-3;
    fSigmaUpsilonPPArray[3] = 0.023e-3;
    fSigmaCCbarArray[3] =     1.68;           // 208*8.1e-3
    fSigmaBBbarArray[3] =     0.061;          // 208*0.29e-3

    fCMSEnergyTeVArray[4] =  -fCMSEnergyTeVArray[3];
    fSigmaReactionArray[4] =  fSigmaReactionArray[3];
    fSigmaJPsiArray[4] =      fSigmaJPsiArray[3];
    fSigmaChic1Array[4] =     fSigmaChic1Array[3];
    fSigmaChic2Array[4] =     fSigmaChic2Array[3];
    fSigmaPsiPArray[4] =      fSigmaPsiPArray[3];
    fSigmaUpsilonArray[4] =   fSigmaUpsilonArray[3];
    fSigmaUpsilonPArray[4] =  fSigmaUpsilonPArray[3];
    fSigmaUpsilonPPArray[4] = fSigmaUpsilonPPArray[3];
    fSigmaCCbarArray[4] =     fSigmaCCbarArray[3];
    fSigmaBBbarArray[4] =     fSigmaBBbarArray[3];

// x-sections for Min. Bias Pb-Pb @ 2.76 TeV: charmonia and bottomonia 
// from 7 TeV numbers scaled according to pQCD ccbar and bbbar x-sections
// and with Glauber scaling
    fCMSEnergyTeVArray[5] =   3.00;           // for 2.76 TeV
    fSigmaReactionArray[5] =  7.65;
    fSigmaJPsiArray[5] =      0.734;          // 208*208*0.505*33.6e-6
    fSigmaChic1Array[5] =     0.712;
    fSigmaChic2Array[5] =     1.175;
    fSigmaPsiPArray[5] =      0.166;
    fSigmaUpsilonArray[5] =   0.0092;         // 208*208*0.379*0.56e-6
    fSigmaUpsilonPArray[5] =  0.0031;
    fSigmaUpsilonPPArray[5] = 0.0015;
    fSigmaCCbarArray[5] =     151.;           // 208*208*3.49e-3
    fSigmaBBbarArray[5] =     3.8;            // 208*208*0.088e-3
    
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

    snprintf(nameJpsi,10, "Jpsi");    
    snprintf(nameChic1,10, "Chic1");
    snprintf(nameChic2,10, "Chic2");
    snprintf(namePsiP,10, "PsiP");
    snprintf(nameUps,10, "Ups");
    snprintf(nameUpsP,10, "UpsP");
    snprintf(nameUpsPP,10, "UpsPP");

    Char_t tname[40] = "";
    if(cmsEnergy == 10)        {snprintf(tname, 40, "CDF pp 10");
    } else if (cmsEnergy == 14){snprintf(tname, 40, "CDF pp");
    } else if (cmsEnergy == 7) {snprintf(tname, 40, "pp 7");
      //    } else if (cmsEnergy == 2) {snprintf(tname, 40, "pp 2.76");
    } else if (cmsEnergy == 9) {snprintf(tname, 40, "pPb 8.8");
      if (fCentralityBin > 0) snprintf(tname, 40, "pPb 8.8c%d",fCentralityBin); 
    } else if (cmsEnergy == -9){snprintf(tname, 40, "Pbp 8.8");
      if (fCentralityBin > 0) snprintf(tname, 40, "Pbp 8.8c%d",fCentralityBin); 
    } else if (cmsEnergy == 3) {snprintf(tname, 40, "PbPb 2.76");
      if (fCentralityBin > 0) snprintf(tname, 40, "PbPb 2.76c%d",fCentralityBin); 
    } else {
	AliError("Initialisation failed, wrong cmsEnergy");
	return;
    }
    genjpsi = new AliGenParam(1, AliGenMUONlib::kJpsi, tname, "Jpsi");
    genchic1 = new AliGenParam(1, AliGenMUONlib::kChic1, tname, "Chic1");
    genchic2 = new AliGenParam(1, AliGenMUONlib::kChic2,  tname, "Chic2");
    genpsiP   = new AliGenParam(1, AliGenMUONlib::kPsiP,   tname, "PsiP");
    genupsilon = new AliGenParam(1, AliGenMUONlib::kUpsilon, tname, "Upsilon");
    genupsilonP = new AliGenParam(1, AliGenMUONlib::kUpsilonP, tname, "UpsilonP");
    genupsilonPP = new AliGenParam(1, AliGenMUONlib::kUpsilonPP, tname, "UpsilonPP");

// Hard process yield per pA or AA collision for i-th centrality bin is R*r[i]*shad[i]
// where R is the ratio of hard and geometrical x-sections, r[i] is the ratio of these
// x-section fractions for given centrality and shad[i] is the shadowing factor (in 4pi).
// The latter is assumed to be the same for HF-hadrons & quarkonia of the same flavour.
    Int_t i = 0;
    Double_t chard[20] = {0};     // charm & beauty shadowing factors are different
    Double_t bhard[20] = {0};
    chard[0] = 1;                 // 1st element for pp and min. bias (MB) collisions
    bhard[0] = 1;

// 4 centrality bins for p-Pb & Pb-p: 0-20-40-60-100 % 
    if (cmsEnergy == 9 || cmsEnergy == -9) {
      const Int_t n9 = 5;         // 1st element for MB collisions
      Double_t r9[n9] = {1, 1.936, 1.473, 0.914, 0.333};        // ratio of hard-over-geo fractions
      Double_t cshad9[n9] = {0.785, 0.715, 0.775, 0.856, 0.951};// EKS98 shadowing factors
      Double_t bshad9[n9] = {0.889, 0.853, 0.884, 0.926, 0.975};
      for(i=0; i<n9; i++) {
	  chard[i] = cshad9[i]*r9[i];   
	  bhard[i] = bshad9[i]*r9[i];
      }
    }

// 11 centrality bins for Pb-Pb: 0-5-10-20-30-40-50-60-70-80-90-100 % 
    if (cmsEnergy == 3) {
      const Int_t n3 = 12;        // 1st element for MB collisions
      Double_t r3[n3] = {1, 4.661, 3.647, 2.551, 1.544, 0.887, 0.474,
			    0.235, 0.106, 0.044, 0.017, 0.007};        // ratio of hard-over-geo fractions
      Double_t cshad3[n3] = {0.662, 0.622, 0.631, 0.650, 0.681, 0.718, 
			     0.760, 0.805, 0.849, 0.888, 0.918, 0.944};// EKS98 shadowing factors
      Double_t bshad3[n3] = {0.874, 0.856, 0.861, 0.869, 0.883, 0.898, 
			     0.915, 0.932, 0.948, 0.962, 0.972, 0.981};
      for(i=0; i<n3; i++) {
	  chard[i] = cshad3[i]*r3[i];   
	  bhard[i] = bshad3[i]*r3[i];
      }
    }

    AddReso2Generator(nameJpsi,genjpsi,chard[fCentralityBin]*sigmajpsi,fJpsiPol);
    AddReso2Generator(nameChic1,genchic1,chard[fCentralityBin]*sigmachic1,fChic1Pol);
    AddReso2Generator(nameChic2,genchic2,chard[fCentralityBin]*sigmachic2,fChic2Pol);
    AddReso2Generator(namePsiP,genpsiP,chard[fCentralityBin]*sigmapsiP,fPsiPPol);

    AddReso2Generator(nameUps,genupsilon,bhard[fCentralityBin]*sigmaupsilon,fUpsPol);
    AddReso2Generator(nameUpsP,genupsilonP,bhard[fCentralityBin]*sigmaupsilonP,fUpsPPol);
    AddReso2Generator(nameUpsPP,genupsilonPP,bhard[fCentralityBin]*sigmaupsilonPP,fUpsPPPol);

//------------------------------------------------------------------
// Generator of charm
    AliGenCorrHF *gencharm = new AliGenCorrHF(1, 4, cmsEnergy);
    gencharm->SetMomentumRange(0,9999);
    gencharm->SetForceDecay(kAll);
    Double_t ratioccbar = chard[fCentralityBin]*sigmaccbar/fSigmaReaction;
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
    Double_t ratiobbbar = bhard[fCentralityBin]*sigmabbbar/fSigmaReaction;
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
      if (decReso) {
	  decReso->SetForceDecay(kAll);
	  decReso->Init();
	  genReso->SetDecayer(decReso);
      }
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
	    gen->SetTime(fTime);

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
    fHeader->SetInteractionTime(fTime);

    gAlice->SetGenEventHeader(fHeader);

//     AliInfo(Form("Generated Events are %d and Succeeded Events are %d",fNGenerated,fNSucceded));
    AliDebug(5,Form("Generated Events are %d and Succeeded Events are %d",fNGenerated,fNSucceded));
}

    
