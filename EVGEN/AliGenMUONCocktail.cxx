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

// Classe to create the MUON coktail for physics in the Alice muon spectrometer
// Gines Martinez, jan 2004, Nantes martinez@in2p3.fr


//

#include <TList.h>
#include <TObjArray.h>
#include <TF1.h>
#include <TParticle.h>

#include "AliGenParam.h"
#include "AliGenMUONlib.h"
#include "AliGenMUONCocktail.h"
#include "AliGenCocktailEntry.h"
#include "AliCollisionGeometry.h"
#include "AliRun.h"
#include "AliMC.h"
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
    delete fEntries;
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

  Float_t sigma_reaction =   0.072;   // MINB pp at LHC energies 72 mb
  
  // Generating J/Psi Physics
  AliGenParam * gen_jpsi = new AliGenParam(1, AliGenMUONlib::kJpsi, "Vogt", "Jpsi");
  // 4pi Generation 
  gen_jpsi->SetPtRange(0,100.);
  gen_jpsi->SetYRange(-8.,8);
  gen_jpsi->SetPhiRange(0.,360.);
  gen_jpsi->SetForceDecay(kDiMuon);
  gen_jpsi->SetTrackingFlag(1);
  // Calculation of the paritcle multiplicity per event in the muonic channel
  Float_t ratio_jpsi; // Ratio with respect to the reaction cross-section for the muonic channel in the kinematics limit of the MUONCocktail
  Float_t sigma_jpsi     = 31.0e-6 * 0.437;   //   section "6.7 Quarkonia Production" table 6.5 for pp  times shadowing
  Float_t br_jpsi        = 0.0588;           // Branching Ratio for JPsi
  gen_jpsi->Init();  // Generating pT and Y parametrsation for the 4pi
  ratio_jpsi = sigma_jpsi * br_jpsi * fNumberOfCollisions / sigma_reaction * gen_jpsi->GetRelativeArea(ptMin,ptMax,yMin,yMax,phiMin,phiMax);
  printf(">>> ratio jpsi %g et %g Ncol %g sigma %g\n",ratio_jpsi,gen_jpsi->GetRelativeArea(ptMin,ptMax,yMin,yMax,phiMin,phiMax),fNumberOfCollisions, sigma_jpsi );
  // Generation in the kinematical limits of AliGenMUONCocktail
  gen_jpsi->SetPtRange(ptMin, ptMax);  
  gen_jpsi->SetYRange(yMin, yMax);
  gen_jpsi->SetPhiRange(phiMin, phiMax);
  gen_jpsi->Init(); // Generating pT and Y parametrsation in the desired kinematic range
  // Adding Generator 
  AddGenerator(gen_jpsi, "Jpsi", ratio_jpsi);
  fTotalRate+=ratio_jpsi;

// Generating Upsilon Physics
  AliGenParam * gen_upsilon = new AliGenParam(1, AliGenMUONlib::kUpsilon, "Vogt", "Upsilon");  
  gen_upsilon->SetPtRange(0,100.);  
  gen_upsilon->SetYRange(-8.,8);
  gen_upsilon->SetPhiRange(0.,360.);
  gen_upsilon->SetForceDecay(kDiMuon);
  gen_upsilon->SetTrackingFlag(1);
  Float_t ratio_upsilon; // Ratio with respect to the reaction cross-section for the muonic channel in the kinematics limit of the MUONCocktail
  Float_t sigma_upsilon     = 0.501e-6 * 0.674;   //  section "6.7 Quarkonia Production" table 6.5 for pp  times shadowing
  Float_t br_upsilon        = 0.0248;  // Branching Ratio for Upsilon
  gen_upsilon->Init();  // Generating pT and Y parametrsation for the 4pi
  ratio_upsilon = sigma_upsilon * br_upsilon * fNumberOfCollisions / sigma_reaction * gen_upsilon->GetRelativeArea(ptMin,ptMax,yMin,yMax,phiMin,phiMax);
  printf(">>> ratio upsilon %g et %g\n",ratio_upsilon, gen_upsilon->GetRelativeArea(ptMin,ptMax,yMin,yMax,phiMin,phiMax));
  gen_upsilon->SetPtRange(ptMin, ptMax);  
  gen_upsilon->SetYRange(yMin, yMax);
  gen_upsilon->SetPhiRange(phiMin, phiMax);
  gen_upsilon->Init(); // Generating pT and Y parametrsation in the desired kinematic range
  AddGenerator(gen_upsilon,"Upsilon", ratio_upsilon);
  fTotalRate+=ratio_upsilon;

// Generating Charm Physics 
  AliGenParam * gen_charm = new AliGenParam(1, AliGenMUONlib::kCharm, "Vogt", "Charm");  
  gen_charm->SetPtRange(0,100.);  
  gen_charm->SetYRange(-8.,8);
  gen_charm->SetPhiRange(0.,360.);
  gen_charm->SetForceDecay(kSemiMuonic);
  gen_charm->SetTrackingFlag(1);
  Float_t ratio_charm; // Ratio with respect to the reaction cross-section for the muonic channel in the kinematics limit of the MUONCocktail
  Float_t sigma_charm     = 2. * 6.64e-3 * 0.65 ;   // 
  Float_t br_charm        = 0.12;  // Branching Ratio for Charm
  gen_charm->Init();  // Generating pT and Y parametrsation for the 4pi
  ratio_charm = sigma_charm * br_charm * fNumberOfCollisions / sigma_reaction * 
    gen_charm->GetRelativeArea(ptMin,ptMax,yMin,yMax,phiMin,phiMax);
  gen_charm->SetPtRange(ptMin, ptMax);  
  gen_charm->SetYRange(yMin, yMax);
  gen_charm->SetPhiRange(phiMin, phiMax);
  gen_charm->Init(); // Generating pT and Y parametrsation in the desired kinematic range

  printf(">>> ratio charm %f\n",ratio_charm);
  AddGenerator(gen_charm,"Charm", ratio_charm);
  fTotalRate+=ratio_charm;

// Generating Beauty Physics "Correlated Pairs"
  AliGenParam * gen_beauty = new AliGenParam(2, AliGenMUONlib::kBeauty, "Vogt", "Beauty");  
  gen_beauty->SetPtRange(0,100.);  
  gen_beauty->SetYRange(-8.,8);
  gen_beauty->SetPhiRange(0.,360.);
  gen_beauty->SetForceDecay(kSemiMuonic);
  gen_beauty->SetTrackingFlag(1);
  Float_t ratio_beauty; // Ratio with respect to the reaction cross-section for the muonic channel in the kinematics limit of the MUONCocktail
  Float_t sigma_beauty     = 2. * 0.210e-3 * 0.84;   // 
  Float_t br_beauty        = 0.15;  // Branching Ratio for Beauty
  gen_beauty->Init();  // Generating pT and Y parametrsation for the 4pi
  ratio_beauty = sigma_beauty * br_beauty * fNumberOfCollisions / sigma_reaction * 
    gen_beauty->GetRelativeArea(ptMin,ptMax,yMin,yMax,phiMin,phiMax);
  gen_beauty->SetPtRange(ptMin, ptMax);  
  gen_beauty->SetYRange(yMin, yMax);
  gen_beauty->SetPhiRange(phiMin, phiMax);
  gen_beauty->Init(); // Generating pT and Y parametrisation in the desired kinematic range

  printf(">>> ratio beauty %f\n",ratio_beauty);
  AddGenerator(gen_beauty,"Beauty", ratio_beauty);
  fTotalRate+=ratio_beauty;

// Generating Pion Physics
  AliGenParam * gen_pion = new AliGenParam(1, AliGenMUONlib::kPion, "Vogt", "Pion");  
  gen_pion->SetPtRange(0,100.);  
  gen_pion->SetYRange(-8.,8);
  gen_pion->SetPhiRange(0.,360.);
  gen_pion->SetForceDecay(kPiToMu);
  gen_pion->SetTrackingFlag(1);
  Float_t ratio_pion; // Ratio with respect to the reaction cross-section for the muonic channel in the kinematics limit of the MUONCocktail
  Float_t sigma_pion     = 1.6e-2;   // A ojo TO be studied in detail.  
  Float_t br_pion        = 0.9999;  // Branching Ratio for Pion 
  gen_pion->Init();  // Generating pT and Y parametrsation for the 4pi
  ratio_pion = sigma_pion * br_pion * fNumberOfParticipants / sigma_reaction * gen_pion->GetRelativeArea(ptMin,ptMax,yMin,yMax,phiMin,phiMax);
  gen_pion->SetPtRange(ptMin, ptMax);  
  gen_pion->SetYRange(yMin, yMax);
  gen_pion->SetPhiRange(phiMin, phiMax);
  gen_pion->Init(); // Generating pT and Y parametrsation in the desired kinematic range
  printf(">>> ratio pion %f\n",ratio_pion);
  AddGenerator(gen_pion,"Pion", ratio_pion);
  fTotalRate+=ratio_pion;

// Generating Kaon Physics
  AliGenParam * gen_kaon = new AliGenParam(1, AliGenMUONlib::kKaon, "Vogt", "Kaon");  
  gen_kaon->SetPtRange(0,100.);  
  gen_kaon->SetYRange(-8.,8);
  gen_kaon->SetPhiRange(0.,360.);
  gen_kaon->SetForceDecay(kKaToMu);
  gen_kaon->SetTrackingFlag(1);
  Float_t ratio_kaon; // Ratio with respect to the reaction cross-section for the muonic channel in the kinematics limit of the MUONCocktail
  Float_t sigma_kaon     = 1.8e-4;   // OJO 
  Float_t br_kaon        = 0.6351 ;  // Branching Ratio for Kaon 
  gen_kaon->Init();  // Generating pT and Y parametrsation for the 4pi
  ratio_kaon = sigma_kaon * br_kaon * fNumberOfParticipants/ sigma_reaction * gen_kaon->GetRelativeArea(ptMin,ptMax,yMin,yMax,phiMin,phiMax);
  gen_kaon->SetPtRange(ptMin, ptMax);  
  gen_kaon->SetYRange(yMin, yMax);
  gen_kaon->SetPhiRange(phiMin, phiMax);
  gen_kaon->Init(); // Generating pT and Y parametrsation in the desired kinematic range
  printf(">>> ratio kaon %f\n",ratio_kaon);
  AddGenerator(gen_kaon,"Kaon", ratio_kaon);
  fTotalRate+=ratio_kaon;
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
    Bool_t PrimordialTrigger = kFALSE;

    while(!PrimordialTrigger) {

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
      if (numberOfMuons >= fMuonMultiplicity ) PrimordialTrigger = kTRUE;
    }
    //printf(">>> Trigger Accepted!!! \n");
    fNSucceded++;
    //   Float_t Ratio = (Float_t) fNSucceded/fNGenerated;
    //    printf("Generated %d, Succeded %d and Ratio %f\n",fNGenerated, fNSucceded,Ratio);
}






