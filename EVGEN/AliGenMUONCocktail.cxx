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
// Classe to create the MUON coktail for physics in the Alice muon spectrometer
// The followoing muons sources are included in this cocktail: 
//     jpsi, upsilon, non-correlated open and beauty, and muons from pion and kaons.
// The free parameeters are :
//      pp reaction cross-section
//      production cross-sections in pp collisions and 
//      branching ratios in the muon channel
// Hard probes are supposed to scale with Ncoll and hadronic production with (0.8Ncoll+0.2*Npart)
// There is a primordial trigger wiche requires :
//      a minimum muon multiplicity above a pT cut in a theta acceptance cone
//
// Gines Martinez, jan 2004, Nantes  martinez@in2p3.fr


//

#include <TObjArray.h>
#include <TParticle.h>

#include "AliGenCocktailEntry.h"
#include "AliGenMUONCocktail.h"
#include "AliGenMUONlib.h"
#include "AliGenParam.h"
#include "AliMC.h"
#include "AliRun.h"
#include "AliStack.h"

ClassImp(AliGenMUONCocktail)
  
  
//________________________________________________________________________
AliGenMUONCocktail::AliGenMUONCocktail()
  :AliGenCocktail()
{
// Constructor
  fTotalRate =0;  
  fNSucceded=0; 
  fNGenerated=0; 
  fMuonMultiplicity=1;
  fMuonPtCut= 1.;
  fMuonThetaMinCut=171.; 
  fMuonThetaMaxCut=178.;
  fNumberOfCollisions = 400; // Minimum bias Pb+Pb collisions
    //
}
//_________________________________________________________________________
AliGenMUONCocktail::AliGenMUONCocktail(const AliGenMUONCocktail & cocktail):
    AliGenCocktail(cocktail)
{
// Copy constructor
  fTotalRate =0; 
  fNSucceded=0; 
  fNGenerated=0; 
  fMuonMultiplicity=1;
  fMuonPtCut= 1.;
  fMuonThetaMinCut=171.; 
  fMuonThetaMaxCut=178.;
  fNumberOfCollisions = 400; // Minimum bias Pb+Pb collisions
    //
}
//_________________________________________________________________________
AliGenMUONCocktail::~AliGenMUONCocktail()
{
// Destructor
}

//_________________________________________________________________________
void AliGenMUONCocktail::Init()
{
  // Defining MUON physics cocktail

  // Kinematical limits for particle generation
  Float_t ptMin  = fPtMin;
  Float_t ptMax  = fPtMax;
  Float_t yMin   = fYMin;;
  Float_t yMax   = fYMax;;
  Float_t phiMin = fPhiMin*180./TMath::Pi();
  Float_t phiMax = fPhiMax*180./TMath::Pi();
  printf(">>> Kinematical range pT:%f:%f  y:%f:%f Phi:%f:%f\n",ptMin,ptMax,yMin,yMax,phiMin,phiMax);

  Float_t sigmaReaction =   0.072;   // MINB pp at LHC energies 72 mb

  // Generating J/Psi Physics
  AliGenParam * genjpsi = new AliGenParam(1, AliGenMUONlib::kJpsi, "Vogt", "Jpsi");
  // 4pi Generation 
  genjpsi->SetPtRange(0,100.);
  genjpsi->SetYRange(-8.,8);
  genjpsi->SetPhiRange(0.,360.);
  genjpsi->SetForceDecay(kDiMuon);
  genjpsi->SetTrackingFlag(1);
  // Calculation of the paritcle multiplicity per event in the muonic channel
  Float_t ratiojpsi; // Ratio with respect to the reaction cross-section for the muonic channel in the kinematics limit of the MUONCocktail
  Float_t sigmajpsi     = 31.0e-6 * 0.437;   //   section "6.7 Quarkonia Production" table 6.5 for pp  times shadowing
  Float_t brjpsi        = 0.0588;           // Branching Ratio for JPsi
  genjpsi->Init();  // Generating pT and Y parametrsation for the 4pi
  ratiojpsi = sigmajpsi * brjpsi * fNumberOfCollisions / sigmaReaction * genjpsi->GetRelativeArea(ptMin,ptMax,yMin,yMax,phiMin,phiMax);
  printf(">>> ratio jpsi %g et %g Ncol %g sigma %g\n",ratiojpsi,genjpsi->GetRelativeArea(ptMin,ptMax,yMin,yMax,phiMin,phiMax),fNumberOfCollisions, sigmajpsi );
  // Generation in the kinematical limits of AliGenMUONCocktail
  genjpsi->SetPtRange(ptMin, ptMax);  
  genjpsi->SetYRange(yMin, yMax);
  genjpsi->SetPhiRange(phiMin, phiMax);
  genjpsi->Init(); // Generating pT and Y parametrsation in the desired kinematic range
  // Adding Generator 
  AddGenerator(genjpsi, "Jpsi", ratiojpsi);
  fTotalRate+=ratiojpsi;

 // Generating Psi prime Physics
  AliGenParam * genpsiP = new AliGenParam(1, AliGenMUONlib::kPsiP, "Vogt", "PsiP");
  // 4pi Generation 
  genpsiP->SetPtRange(0,100.);
  genpsiP->SetYRange(-8.,8);
  genpsiP->SetPhiRange(0.,360.);
  genpsiP->SetForceDecay(kDiMuon);
  genpsiP->SetTrackingFlag(1);
  // Calculation of the paritcle multiplicity per event in the muonic channel
  Float_t ratiopsiP; // Ratio with respect to the reaction cross-section for the muonic channel in the kinematics limit of the MUONCocktail
  Float_t sigmapsiP     = 4.68e-6 * 0.437;   //   section "6.7 Quarkonia Production" table 6.5 for pp  times shadowing
  Float_t brpsiP        = 0.0103;           // Branching Ratio for PsiP
  genpsiP->Init();  // Generating pT and Y parametrsation for the 4pi
  ratiopsiP = sigmapsiP * brpsiP * fNumberOfCollisions / sigmaReaction * genpsiP->GetRelativeArea(ptMin,ptMax,yMin,yMax,phiMin,phiMax);
  printf(">>> ratio psiP %g et %g Ncol %g sigma %g\n",ratiopsiP,genpsiP->GetRelativeArea(ptMin,ptMax,yMin,yMax,phiMin,phiMax),fNumberOfCollisions, sigmapsiP );
  // Generation in the kinematical limits of AliGenMUONCocktail
  genpsiP->SetPtRange(ptMin, ptMax);  
  genpsiP->SetYRange(yMin, yMax);
  genpsiP->SetPhiRange(phiMin, phiMax);
  genpsiP->Init(); // Generating pT and Y parametrsation in the desired kinematic range
  // Adding Generator 
  AddGenerator(genpsiP, "PsiP", ratiopsiP);
  fTotalRate+=ratiopsiP;


// Generating Upsilon Physics
  AliGenParam * genupsilon = new AliGenParam(1, AliGenMUONlib::kUpsilon, "Vogt", "Upsilon");  
  genupsilon->SetPtRange(0,100.);  
  genupsilon->SetYRange(-8.,8);
  genupsilon->SetPhiRange(0.,360.);
  genupsilon->SetForceDecay(kDiMuon);
  genupsilon->SetTrackingFlag(1);
  Float_t ratioupsilon; // Ratio with respect to the reaction cross-section for the muonic channel in the kinematics limit of the MUONCocktail
  Float_t sigmaupsilon     = 0.501e-6 * 0.674;   //  section "6.7 Quarkonia Production" table 6.5 for pp  times shadowing
  Float_t brupsilon        = 0.0248;  // Branching Ratio for Upsilon
  genupsilon->Init();  // Generating pT and Y parametrsation for the 4pi
  ratioupsilon = sigmaupsilon * brupsilon * fNumberOfCollisions / sigmaReaction * genupsilon->GetRelativeArea(ptMin,ptMax,yMin,yMax,phiMin,phiMax);
  printf(">>> ratio upsilon %g et %g\n",ratioupsilon, genupsilon->GetRelativeArea(ptMin,ptMax,yMin,yMax,phiMin,phiMax));
  genupsilon->SetPtRange(ptMin, ptMax);  
  genupsilon->SetYRange(yMin, yMax);
  genupsilon->SetPhiRange(phiMin, phiMax);
  genupsilon->Init(); // Generating pT and Y parametrsation in the desired kinematic range
  AddGenerator(genupsilon,"Upsilon", ratioupsilon);
  fTotalRate+=ratioupsilon;

// Generating UpsilonP Physics
  AliGenParam * genupsilonP = new AliGenParam(1, AliGenMUONlib::kUpsilonP, "Vogt", "UpsilonP");  
  genupsilonP->SetPtRange(0,100.);  
  genupsilonP->SetYRange(-8.,8);
  genupsilonP->SetPhiRange(0.,360.);
  genupsilonP->SetForceDecay(kDiMuon);
  genupsilonP->SetTrackingFlag(1);
  Float_t ratioupsilonP; // Ratio with respect to the reaction cross-section for the muonic channel in the kinematics limit of the MUONCocktail
  Float_t sigmaupsilonP     = 0.246e-6 * 0.674;   //  section "6.7 Quarkonia Production" table 6.5 for pp  times shadowing
  Float_t brupsilonP        = 0.0131;  // Branching Ratio for UpsilonP
  genupsilonP->Init();  // Generating pT and Y parametrsation for the 4pi
  ratioupsilonP = sigmaupsilonP * brupsilonP * fNumberOfCollisions / sigmaReaction * genupsilonP->GetRelativeArea(ptMin,ptMax,yMin,yMax,phiMin,phiMax);
  printf(">>> ratio upsilonP %g et %g\n",ratioupsilonP, genupsilonP->GetRelativeArea(ptMin,ptMax,yMin,yMax,phiMin,phiMax));
  genupsilonP->SetPtRange(ptMin, ptMax);  
  genupsilonP->SetYRange(yMin, yMax);
  genupsilonP->SetPhiRange(phiMin, phiMax);
  genupsilonP->Init(); // Generating pT and Y parametrsation in the desired kinematic range
  AddGenerator(genupsilonP,"UpsilonP", ratioupsilonP);
  fTotalRate+=ratioupsilonP;


// Generating UpsilonPP Physics
  AliGenParam * genupsilonPP = new AliGenParam(1, AliGenMUONlib::kUpsilonPP, "Vogt", "UpsilonPP");  
  genupsilonPP->SetPtRange(0,100.);  
  genupsilonPP->SetYRange(-8.,8);
  genupsilonPP->SetPhiRange(0.,360.);
  genupsilonPP->SetForceDecay(kDiMuon);
  genupsilonPP->SetTrackingFlag(1);
  Float_t ratioupsilonPP; // Ratio with respect to the reaction cross-section for the muonic channel in the kinematics limit of the MUONCocktail
  Float_t sigmaupsilonPP     = 0.100e-6 * 0.674;   //  section "6.7 Quarkonia Production" table 6.5 for pp  times shadowing
  Float_t brupsilonPP        = 0.0181;  // Branching Ratio for UpsilonPP
  genupsilonPP->Init();  // Generating pT and Y parametrsation for the 4pi
  ratioupsilonPP = sigmaupsilonPP * brupsilonPP * fNumberOfCollisions / sigmaReaction * genupsilonPP->GetRelativeArea(ptMin,ptMax,yMin,yMax,phiMin,phiMax);
  printf(">>> ratio upsilonPP %g et %g\n",ratioupsilonPP, genupsilonPP->GetRelativeArea(ptMin,ptMax,yMin,yMax,phiMin,phiMax));
  genupsilonPP->SetPtRange(ptMin, ptMax);  
  genupsilonPP->SetYRange(yMin, yMax);
  genupsilonPP->SetPhiRange(phiMin, phiMax);
  genupsilonPP->Init(); // Generating pT and Y parametrsation in the desired kinematic range
  AddGenerator(genupsilonPP,"UpsilonPP", ratioupsilonPP);
  fTotalRate+=ratioupsilonPP;


// Generating Charm Physics 
  AliGenParam * gencharm = new AliGenParam(1, AliGenMUONlib::kCharm, "Vogt", "Charm");  
  gencharm->SetPtRange(0,100.);  
  gencharm->SetYRange(-8.,8);
  gencharm->SetPhiRange(0.,360.);
  gencharm->SetForceDecay(kSemiMuonic);
  gencharm->SetTrackingFlag(1);
  Float_t ratiocharm; // Ratio with respect to the reaction cross-section for the muonic channel in the kinematics limit of the MUONCocktail
  Float_t sigmacharm     = 2. * 6.64e-3 * 0.65 ;   // 
  Float_t brcharm        = 0.12;  // Branching Ratio for Charm
  gencharm->Init();  // Generating pT and Y parametrsation for the 4pi
  ratiocharm = sigmacharm * brcharm * fNumberOfCollisions / sigmaReaction * 
    gencharm->GetRelativeArea(ptMin,ptMax,yMin,yMax,phiMin,phiMax);
  gencharm->SetPtRange(ptMin, ptMax);  
  gencharm->SetYRange(yMin, yMax);
  gencharm->SetPhiRange(phiMin, phiMax);
  gencharm->Init(); // Generating pT and Y parametrsation in the desired kinematic range

  printf(">>> ratio charm %f\n",ratiocharm);
  AddGenerator(gencharm,"Charm", ratiocharm);
  fTotalRate+=ratiocharm;

// Generating Beauty Physics "Correlated Pairs"
  AliGenParam * genbeauty = new AliGenParam(2, AliGenMUONlib::kBeauty, "Vogt", "Beauty");  
  genbeauty->SetPtRange(0,100.);  
  genbeauty->SetYRange(-8.,8);
  genbeauty->SetPhiRange(0.,360.);
  genbeauty->SetForceDecay(kSemiMuonic);
  genbeauty->SetTrackingFlag(1);
  Float_t ratiobeauty; // Ratio with respect to the reaction cross-section for the muonic channel in the kinematics limit of the MUONCocktail
  Float_t sigmabeauty     = 2. * 0.210e-3 * 0.84;   // 
  Float_t brbeauty        = 0.15;  // Branching Ratio for Beauty
  genbeauty->Init();  // Generating pT and Y parametrsation for the 4pi
  ratiobeauty = sigmabeauty * brbeauty * fNumberOfCollisions / sigmaReaction * 
    genbeauty->GetRelativeArea(ptMin,ptMax,yMin,yMax,phiMin,phiMax);
  genbeauty->SetPtRange(ptMin, ptMax);  
  genbeauty->SetYRange(yMin, yMax);
  genbeauty->SetPhiRange(phiMin, phiMax);
  genbeauty->Init(); // Generating pT and Y parametrisation in the desired kinematic range

  printf(">>> ratio beauty %f\n",ratiobeauty);
  AddGenerator(genbeauty,"Beauty", ratiobeauty);
  fTotalRate+=ratiobeauty;

// Generating Pion Physics
  AliGenParam * genpion = new AliGenParam(1, AliGenMUONlib::kPion, "Vogt", "Pion");  
  genpion->SetPtRange(0,100.);  
  genpion->SetYRange(-8.,8);
  genpion->SetPhiRange(0.,360.);
  genpion->SetForceDecay(kPiToMu);
  genpion->SetTrackingFlag(1);
  Float_t ratiopion; // Ratio with respect to the reaction cross-section for the muonic channel in the kinematics limit of the MUONCocktail
  Float_t sigmapion     = 0.93e-2; // Valerie presentation Clermont-16-jan-2004 and Alice-int-2002-06
  Float_t brpion        = 0.9999;  // Branching Ratio for Pion 
  genpion->Init();  // Generating pT and Y parametrsation for the 4pi
  ratiopion = sigmapion * brpion *  (0.80*fNumberOfParticipants+0.2*fNumberOfCollisions) / sigmaReaction * genpion->GetRelativeArea(ptMin,ptMax,yMin,yMax,phiMin,phiMax);
  genpion->SetPtRange(ptMin, ptMax);  
  genpion->SetYRange(yMin, yMax);
  genpion->SetPhiRange(phiMin, phiMax);
  genpion->Init(); // Generating pT and Y parametrsation in the desired kinematic range
  printf(">>> ratio pion %f\n",ratiopion);
  AddGenerator(genpion,"Pion", ratiopion);
  fTotalRate+=ratiopion;

// Generating Kaon Physics
  AliGenParam * genkaon = new AliGenParam(1, AliGenMUONlib::kKaon, "Vogt", "Kaon");  
  genkaon->SetPtRange(0,100.);  
  genkaon->SetYRange(-8.,8);
  genkaon->SetPhiRange(0.,360.);
  genkaon->SetForceDecay(kKaToMu);
  genkaon->SetTrackingFlag(1);
  Float_t ratiokaon; // Ratio with respect to the reaction cross-section for the muonic channel in the kinematics limit of the MUONCocktail
  Float_t sigmakaon     = 1.23e-4;   // Valerie presentation Clermont-16-jan-2004 and Alice-int-2002-06 
  Float_t brkaon        = 0.6351 ;  // Branching Ratio for Kaon 
  genkaon->Init();  // Generating pT and Y parametrsation for the 4pi
  ratiokaon = sigmakaon * brkaon * (0.80*fNumberOfParticipants+0.2*fNumberOfCollisions)/ sigmaReaction * genkaon->GetRelativeArea(ptMin,ptMax,yMin,yMax,phiMin,phiMax);
  genkaon->SetPtRange(ptMin, ptMax);  
  genkaon->SetYRange(yMin, yMax);
  genkaon->SetPhiRange(phiMin, phiMax);
  genkaon->Init(); // Generating pT and Y parametrsation in the desired kinematic range
  printf(">>> ratio kaon %f\n",ratiokaon);
  AddGenerator(genkaon,"Kaon", ratiokaon);
  fTotalRate+=ratiokaon;
}

//_________________________________________________________________________
void AliGenMUONCocktail::Generate()
{
  //
// Generate event 
    TIter next(fEntries);
    AliGenCocktailEntry *entry = 0;
    AliGenCocktailEntry *preventry = 0;
    AliGenerator* gen = 0;

    TObjArray *partArray = gAlice->GetMCApp()->Particles();

//
//  Generate the vertex position used by all generators
//    
    if(fVertexSmear == kPerEvent) Vertex();
    Bool_t primordialTrigger = kFALSE;

    while(!primordialTrigger) {

      //Reseting stack
      AliRunLoader * runloader = gAlice->GetRunLoader();
      if (runloader)
	if (runloader->Stack())
	  runloader->Stack()->Reset();
      //
      // Loop over generators and generate events
      Int_t igen=0;
      Int_t npart =0;
      
      while((entry = (AliGenCocktailEntry*)next())) {
	gen = entry->Generator();
	gen->SetVertex(fVertex.At(0), fVertex.At(1), fVertex.At(2));
	if ( (npart = gRandom->Poisson(entry->Rate())) >0 ) {
	  igen++;	
	  if (igen ==1) entry->SetFirst(0);
	  else  entry->SetFirst((partArray->GetEntriesFast())+1);
	  gen->SetNumberParticles(npart);
	  gen->Generate();
	  entry->SetLast(partArray->GetEntriesFast());
	  preventry = entry;
	}
      }  
      next.Reset();

      // Tesitng primordial trigger : Muon  pair in the MUON spectrometer acceptance 171,178 and pTCut
      Int_t iPart;
      fNGenerated++;
      Int_t numberOfMuons=0;
      //      printf(">>>fNGenerated is %d\n",fNGenerated);
      for(iPart=0; iPart<partArray->GetEntriesFast(); iPart++){      
	//	gAlice->GetMCApp()->Particle(iPart)->Print();
	if ( (TMath::Abs(gAlice->GetMCApp()->Particle(iPart)->GetPdgCode())==13)  &&
	     (gAlice->GetMCApp()->Particle(iPart)->Theta()*180./TMath::Pi()>fMuonThetaMinCut) &&
	     (gAlice->GetMCApp()->Particle(iPart)->Theta()*180./TMath::Pi()<fMuonThetaMaxCut) &&
	     (gAlice->GetMCApp()->Particle(iPart)->Pt()>fMuonPtCut)                             ) { 
	  gAlice->GetMCApp()->Particle(iPart)->SetProductionVertex(fVertex.At(0), fVertex.At(1), fVertex.At(2), 0.);   
	  numberOfMuons++;
	}
      }
      //  printf(">>> Number of Muons is %d \n", numberOfMuons);
      if (numberOfMuons >= fMuonMultiplicity ) primordialTrigger = kTRUE;
    }
    //printf(">>> Trigger Accepted!!! \n");
    fNSucceded++;
    //   Float_t Ratio = (Float_t) fNSucceded/fNGenerated;
    //    printf("Generated %d, Succeded %d and Ratio %f\n",fNGenerated, fNSucceded,Ratio);
}






