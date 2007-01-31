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
//     jpsi, upsilon, non-correlated open and beauty.
// The free parameeters are :
//      pp reaction cross-section
//      production cross-sections in pp collisions and 
//      branching ratios in the muon channel and shadowing
// Hard probes are supposed to scale with Ncoll and hadronic muon production with (0.8Ncoll+0.2*Npart)
// There is a primordial trigger which requires :
//      a minimum muon multiplicity above a pT cut in a theta acceptance cone
//

#include <TObjArray.h>
#include <TParticle.h>
#include <TF1.h>

#include "AliGenCocktailEntry.h"
#include "AliGenMUONCocktail.h"
#include "AliGenMUONlib.h"
#include "AliFastGlauber.h"
#include "AliGenParam.h"
#include "AliMC.h"
#include "AliRun.h"
#include "AliStack.h"
#include "AliLog.h"

ClassImp(AliGenMUONCocktail)
  
  
//________________________________________________________________________
AliGenMUONCocktail::AliGenMUONCocktail()
    :AliGenCocktail(),
     fFastGlauber (0x0),
     fTotalRate(0),  
     fMuonMultiplicity(1),
     fMuonPtCut(1.),
     fMuonThetaMinCut(171.), 
     fMuonThetaMaxCut(178.),
     fNSucceded(0), 
     fNGenerated(0), 
     fLowImpactParameter(0.),
     fHighImpactParameter(5.),
     fAverageImpactParameter(0.),
     fNumberOfCollisions(0.), 
     fNumberOfParticipants(0.),
     fHadronicMuons(kTRUE),
     fInvMassCut (kFALSE),
     fInvMassMinCut (0.),
     fInvMassMaxCut (100.)
{
// Constructor
}
//_________________________________________________________________________
AliGenMUONCocktail::~AliGenMUONCocktail()
{
// Destructor
  if (fFastGlauber) delete fFastGlauber;
}

//_________________________________________________________________________
void AliGenMUONCocktail::Init()
{
  // NN cross section
  Double_t sigmaReaction =   0.072;   // MinBias NN cross section for PbPb LHC energies  http://arxiv.org/pdf/nucl-ex/0302016

  // Initialising Fast Glauber object
  fFastGlauber = new AliFastGlauber();
  fFastGlauber->SetPbPbLHC();
  fFastGlauber->SetNNCrossSection(sigmaReaction*1000.); //Expected NN cross-section in mb at LHC with diffractive part http://arxiv.org/pdf/nucl-ex/0302016 )
  fFastGlauber->Init(1);
 
  // Calculating average number of collisions
  Int_t ib=0;
  Int_t ibmax=10000;
  Double_t b       = 0.;
  fAverageImpactParameter=0.;
  fNumberOfCollisions   = 0.;
  fNumberOfParticipants = 0.;
  for(ib=0; ib<ibmax; ib++) {
    b = fFastGlauber->GetRandomImpactParameter(fLowImpactParameter,fHighImpactParameter);
    fAverageImpactParameter+=b;
    fNumberOfCollisions    += fFastGlauber->GetNumberOfCollisions( b )/(1.-TMath::Exp(-fFastGlauber->GetNumberOfCollisions(b)));
    fNumberOfParticipants  += fFastGlauber->GetNumberOfParticipants( b );
  } 
  fAverageImpactParameter/= ((Double_t) ibmax);
  fNumberOfCollisions    /= ((Double_t) ibmax);
  fNumberOfParticipants  /= ((Double_t) ibmax);;
  AliInfo(Form("<b>=%4.2f, <Ncoll>=%5.1f and and <Npart>=%5.1f",b, fNumberOfCollisions, fNumberOfParticipants));

  // Estimating shadowing on charm a beaty production
  // -----------------------------------------------------
  // Extrapolation of the cross sections from $p-p$ to \mbox{Pb--Pb}
  // interactions
  // is done by means of the Glauber model. For the impact parameter dependence
  // of the shadowing factor we use a simple formula:
  // $C_{sh}(b) = C_{sh}(0) + (1 - C_{sh}(0))(b/16~fm)4$,
  // motivated by the theoretical predictions (see e.g.
  // V. Emelyanov et al., Phys. Rev. C61, 044904 (2000)) and HIJING
  // simulations showing an almost flat behaviour
  // up to 10~$fm$ and a rapid increase to 1 for larger impact parameters.
  // C_{sh}(0)  = 0.60 for Psi and 0.76 for Upsilon (Smba communication). 
  // for open charm and beauty is 0.65 and 0.84
  // ----------------------------------------------------- 
  Double_t charmshadowing       = 0.65 + (1.0-0.65)*TMath::Power(fAverageImpactParameter/16.,4);
  Double_t beautyshadowing      = 0.84 + (1.0-0.84)*TMath::Power(fAverageImpactParameter/16.,4);
  Double_t charmoniumshadowing  = 0.60 + (1.0-0.60)*TMath::Power(fAverageImpactParameter/16.,4);
  Double_t beautoniumshadowing  = 0.76 + (1.0-0.76)*TMath::Power(fAverageImpactParameter/16.,4);
  if (fAverageImpactParameter>16.) {
    charmoniumshadowing = 1.0;
    beautoniumshadowing = 1.0;
    charmshadowing = 1.0;
    beautyshadowing= 1.0;
  }
  AliInfo(Form("Shadowing for charmonium and beautonium production are %4.2f and %4.2f respectively",charmoniumshadowing,beautoniumshadowing));
  AliInfo(Form("Shadowing for charm and beauty production are %4.2f and %4.2f respectively",charmshadowing,beautyshadowing));

  // Defining MUON physics cocktail
  // Kinematical limits for particle generation
  Double_t ptMin  = fPtMin;
  Double_t ptMax  = fPtMax;
  Double_t yMin   = fYMin;;
  Double_t yMax   = fYMax;;
  Double_t phiMin = fPhiMin*180./TMath::Pi();
  Double_t phiMax = fPhiMax*180./TMath::Pi();
  AliInfo(Form("Ranges pT:%4.1f : %4.1f GeV/c, y:%4.2f : %4.2f, Phi:%5.1f : %5.1f degres",ptMin,ptMax,yMin,yMax,phiMin,phiMax));

  // Generating J/Psi Physics 
  // Using CFD scaled distribution (see http://clrwww.in2p3.fr/DIMUON2004/talks/sgrigoryan.pdf )
  AliGenParam * genjpsi = new AliGenParam(1, AliGenMUONlib::kJpsi, "CDF scaled", "Jpsi");
  genjpsi->SetPtRange(0,100.); // 4pi generation
  genjpsi->SetYRange(-8.,8);
  genjpsi->SetPhiRange(0.,360.);
  genjpsi->SetForceDecay(kDiMuon);
  genjpsi->SetTrackingFlag(1);
  // Calculation of the particle multiplicity per event in the muonic channel
  Double_t ratiojpsi; // Ratio with respect to the reaction cross-section for the muonic channel in the kinematics limit of the MUONCocktail
  Double_t sigmajpsi  = 31.0e-6 * charmoniumshadowing;   //   section "6.7 Quarkonia Production" table 6.5 for pp  times shadowing
  Double_t brjpsi     = 0.0588;              // Branching Ratio for JPsi PDG PRC15 (200)
  genjpsi->Init();  // Generating pT and Y parametrsation for the 4pi
  ratiojpsi = sigmajpsi * brjpsi * fNumberOfCollisions / sigmaReaction * genjpsi->GetRelativeArea(ptMin,ptMax,yMin,yMax,phiMin,phiMax);
  AliInfo(Form("Jpsi production cross-section in pp with shadowing %5.3g barns",sigmajpsi));
  AliInfo(Form("Jpsi production probability per collisions in acceptance via the muonic channel %5.3g",ratiojpsi));
  // Generation in the kinematical limits of AliGenMUONCocktail
  genjpsi->SetPtRange(ptMin, ptMax);  
  genjpsi->SetYRange(yMin, yMax);
  genjpsi->SetPhiRange(phiMin, phiMax);
  genjpsi->Init(); // Generating pT and Y parametrsation in the desired kinematic range
  // Adding Generator 
  AddGenerator(genjpsi, "Jpsi", ratiojpsi);
  fTotalRate+=ratiojpsi;

 // Generating Psi prime Physics
 // Using CFD scaled distribution (see http://clrwww.in2p3.fr/DIMUON2004/talks/sgrigoryan.pdf )
  AliGenParam * genpsiP = new AliGenParam(1, AliGenMUONlib::kPsiP, "CDF scaled", "PsiP");
  genpsiP->SetPtRange(0,100.);// 4pi generation
  genpsiP->SetYRange(-8.,8);
  genpsiP->SetPhiRange(0.,360.);
  genpsiP->SetForceDecay(kDiMuon);
  genpsiP->SetTrackingFlag(1);
  // Calculation of the paritcle multiplicity per event in the muonic channel
  Double_t ratiopsiP; // Ratio with respect to the reaction cross-section for the muonic channel in the kinematics limit of the MUONCocktail
  Double_t sigmapsiP     = 4.68e-6 * charmoniumshadowing;   //   section "6.7 Quarkonia Production" table 6.5 for pp  times shadowing
  Double_t brpsiP        = 0.0103;           // Branching Ratio for PsiP
  genpsiP->Init();  // Generating pT and Y parametrsation for the 4pi
  ratiopsiP = sigmapsiP * brpsiP * fNumberOfCollisions / sigmaReaction * genpsiP->GetRelativeArea(ptMin,ptMax,yMin,yMax,phiMin,phiMax);
  AliInfo(Form("Psi prime production cross-section in pp with shadowing %5.3g barns",sigmapsiP));
  AliInfo(Form("Psi prime production probability per collisions in acceptance via the muonic channel %5.3g",ratiopsiP));
  // Generation in the kinematical limits of AliGenMUONCocktail
  genpsiP->SetPtRange(ptMin, ptMax);  
  genpsiP->SetYRange(yMin, yMax);
  genpsiP->SetPhiRange(phiMin, phiMax);
  genpsiP->Init(); // Generating pT and Y parametrsation in the desired kinematic range
  // Adding Generator 
  AddGenerator(genpsiP, "PsiP", ratiopsiP);
  fTotalRate+=ratiopsiP;

  // Generating Upsilon Physics
  // Using CFD scaled distribution (see http://clrwww.in2p3.fr/DIMUON2004/talks/sgrigoryan.pdf )
  AliGenParam * genupsilon = new AliGenParam(1, AliGenMUONlib::kUpsilon, "CDF scaled", "Upsilon");  
  genupsilon->SetPtRange(0,100.);  
  genupsilon->SetYRange(-8.,8);
  genupsilon->SetPhiRange(0.,360.);
  genupsilon->SetForceDecay(kDiMuon);
  genupsilon->SetTrackingFlag(1);
  Double_t ratioupsilon; // Ratio with respect to the reaction cross-section for the muonic channel in the kinematics limit of the MUONCocktail
  Double_t sigmaupsilon     = 0.501e-6 * beautoniumshadowing;   //  section "6.7 Quarkonia Production" table 6.5 for pp  times shadowing
  Double_t brupsilon        = 0.0248;  // Branching Ratio for Upsilon
  genupsilon->Init();  // Generating pT and Y parametrsation for the 4pi
  ratioupsilon = sigmaupsilon * brupsilon * fNumberOfCollisions / sigmaReaction * genupsilon->GetRelativeArea(ptMin,ptMax,yMin,yMax,phiMin,phiMax);
  AliInfo(Form("Upsilon 1S production cross-section in pp with shadowing %5.3g barns",sigmaupsilon));
  AliInfo(Form("Upsilon 1S production probability per collisions in acceptance via the muonic channel %5.3g",ratioupsilon));
  genupsilon->SetPtRange(ptMin, ptMax);  
  genupsilon->SetYRange(yMin, yMax);
  genupsilon->SetPhiRange(phiMin, phiMax);
  genupsilon->Init(); // Generating pT and Y parametrsation in the desired kinematic range
  AddGenerator(genupsilon,"Upsilon", ratioupsilon);
  fTotalRate+=ratioupsilon;

  // Generating UpsilonP Physics
  // Using CFD scaled distribution (see http://clrwww.in2p3.fr/DIMUON2004/talks/sgrigoryan.pdf )
  AliGenParam * genupsilonP = new AliGenParam(1, AliGenMUONlib::kUpsilonP, "CDF Scaled", "UpsilonP");  
  genupsilonP->SetPtRange(0,100.);  
  genupsilonP->SetYRange(-8.,8);
  genupsilonP->SetPhiRange(0.,360.);
  genupsilonP->SetForceDecay(kDiMuon);
  genupsilonP->SetTrackingFlag(1);
  Double_t ratioupsilonP; // Ratio with respect to the reaction cross-section for the muonic channel in the kinematics limit of the MUONCocktail
  Double_t sigmaupsilonP     = 0.246e-6 * beautoniumshadowing;   //  section "6.7 Quarkonia Production" table 6.5 for pp  times shadowing
  Double_t brupsilonP        = 0.0131;  // Branching Ratio for UpsilonP
  genupsilonP->Init();  // Generating pT and Y parametrsation for the 4pi
  ratioupsilonP = sigmaupsilonP * brupsilonP * fNumberOfCollisions / sigmaReaction * genupsilonP->GetRelativeArea(ptMin,ptMax,yMin,yMax,phiMin,phiMax);
  AliInfo(Form("Upsilon 2S production cross-section in pp with shadowing %5.3g barns",sigmaupsilonP));
  AliInfo(Form("Upsilon 2S production probability per collisions in acceptance via the muonic channel %5.3g",ratioupsilonP));
  genupsilonP->SetPtRange(ptMin, ptMax);  
  genupsilonP->SetYRange(yMin, yMax);
  genupsilonP->SetPhiRange(phiMin, phiMax);
  genupsilonP->Init(); // Generating pT and Y parametrsation in the desired kinematic range
  AddGenerator(genupsilonP,"UpsilonP", ratioupsilonP);
  fTotalRate+=ratioupsilonP;

  // Generating UpsilonPP Physics
  // Using CFD scaled distribution (see http://clrwww.in2p3.fr/DIMUON2004/talks/sgrigoryan.pdf )
  AliGenParam * genupsilonPP = new AliGenParam(1, AliGenMUONlib::kUpsilonPP, "CDF Scaled", "UpsilonPP");  
  genupsilonPP->SetPtRange(0,100.);  
  genupsilonPP->SetYRange(-8.,8);
  genupsilonPP->SetPhiRange(0.,360.);
  genupsilonPP->SetForceDecay(kDiMuon);
  genupsilonPP->SetTrackingFlag(1);
  Double_t ratioupsilonPP; // Ratio with respect to the reaction cross-section for the muonic channel in the kinematics limit of the MUONCocktail
  Double_t sigmaupsilonPP     = 0.100e-6 * beautoniumshadowing;  //  section "6.7 Quarkonia Production" table 6.5 for pp  times shadowing
  Double_t brupsilonPP        = 0.0181;  // Branching Ratio for UpsilonPP
  genupsilonPP->Init();  // Generating pT and Y parametrsation for the 4pi
  ratioupsilonPP = sigmaupsilonPP * brupsilonPP * fNumberOfCollisions / sigmaReaction * genupsilonPP->GetRelativeArea(ptMin,ptMax,yMin,yMax,phiMin,phiMax);
  AliInfo(Form("Upsilon 3S production cross-section in pp with shadowing %5.3g barns",sigmaupsilonPP));
  AliInfo(Form("Upsilon 3S production probability per collisions in acceptance via the muonic channel %5.3g",ratioupsilonPP));
  genupsilonPP->SetPtRange(ptMin, ptMax);  
  genupsilonPP->SetYRange(yMin, yMax);
  genupsilonPP->SetPhiRange(phiMin, phiMax);
  genupsilonPP->Init(); // Generating pT and Y parametrsation in the desired kinematic range
  AddGenerator(genupsilonPP,"UpsilonPP", ratioupsilonPP);
  fTotalRate+=ratioupsilonPP;

// Generating non-correlated Charm Physics 
  AliGenParam * gencharm = new AliGenParam(1, AliGenMUONlib::kCharm, "pp", "Charm");  
  gencharm->SetPtRange(0,100.);  
  gencharm->SetYRange(-8.,8);
  gencharm->SetPhiRange(0.,360.);
  gencharm->SetForceDecay(kSemiMuonic);
  gencharm->SetTrackingFlag(1);
  Double_t ratiocharm; // Ratio with respect to the reaction cross-section for the muonic channel in the kinematics limit of the MUONCocktail
  Double_t sigmacharm     = 2. * 6.64e-3 * charmshadowing ;   // 
  Double_t brcharm        = 0.12;  // Branching Ratio for Charm
  gencharm->Init();  // Generating pT and Y parametrsation for the 4pi
  ratiocharm = sigmacharm * brcharm * fNumberOfCollisions / sigmaReaction * gencharm->GetRelativeArea(ptMin,ptMax,yMin,yMax,phiMin,phiMax);
  AliInfo(Form("Charm production cross-section in pp with shadowing %5.3g barns",sigmacharm));
  AliInfo(Form("Charm production probability per collisions in acceptance via the semi-muonic channel %5.3g",ratiocharm));
  gencharm->SetPtRange(ptMin, ptMax);  
  gencharm->SetYRange(yMin, yMax);
  gencharm->SetPhiRange(phiMin, phiMax);
  gencharm->Init(); // Generating pT and Y parametrsation in the desired kinematic range
  AddGenerator(gencharm,"Charm", ratiocharm);
  fTotalRate+=ratiocharm;

// Generating non-correlated Beauty Physics 
  AliGenParam * genbeauty = new AliGenParam(1, AliGenMUONlib::kBeauty, "pp", "Beauty");  
  genbeauty->SetPtRange(0,100.);  
  genbeauty->SetYRange(-8.,8);
  genbeauty->SetPhiRange(0.,360.);
  genbeauty->SetForceDecay(kSemiMuonic);
  genbeauty->SetTrackingFlag(1);
  Double_t ratiobeauty; // Ratio with respect to the reaction cross-section for the muonic channel in the kinematics limit of the MUONCocktail
  Double_t sigmabeauty     = 2. * 0.210e-3 * beautyshadowing;   // 
  Double_t brbeauty        = 0.15;  // Branching Ratio for Beauty
  genbeauty->Init();  // Generating pT and Y parametrsation for the 4pi
  ratiobeauty = sigmabeauty * brbeauty * fNumberOfCollisions / sigmaReaction * 
    genbeauty->GetRelativeArea(ptMin,ptMax,yMin,yMax,phiMin,phiMax);
  AliInfo(Form("Beauty production cross-section in pp with shadowing %5.3g barns",sigmabeauty));
  AliInfo(Form("Beauty production probability per collisions in acceptance via the semi-muonic channel %5.3g",ratiobeauty));
  genbeauty->SetPtRange(ptMin, ptMax);  
  genbeauty->SetYRange(yMin, yMax);
  genbeauty->SetPhiRange(phiMin, phiMax);
  genbeauty->Init(); // Generating pT and Y parametrisation in the desired kinematic range
  AddGenerator(genbeauty,"Beauty", ratiobeauty);
  fTotalRate+=ratiobeauty;

  // Only if hadronic muons are included in the cocktail
  if(fHadronicMuons) { 
    // Generating Pion Physics
    // The scaling with Npart and Ncoll has been obtained to reproduced tha values presented by Valeri lors de presentatation 
    // a Clermont Ferrand http://pcrochet.home.cern.ch/pcrochet/files/valerie.pdf
    //  b range(fm)   Ncoll Npart N_mu pT>0.4 GeV/c
    //    0 -  3      1982   381    3.62
    //    3 -  6      1388   297    3.07
    //    6 -  9       674   177    1.76
    //    9 - 12       188    71    0.655
    //   12 - 16        15    10    0.086
    // We found the hadronic muons scales quite well with the number of participants  
    AliGenParam * genpion = new AliGenParam(1, AliGenMUONlib::kPion, "default", "Pion");  
    genpion->SetPtRange(0,100.);  
    genpion->SetYRange(-8.,8);
    genpion->SetPhiRange(0.,360.);
    genpion->SetForceDecay(kPiToMu);
    genpion->SetTrackingFlag(1);
    Double_t ratiopion; // Ratio with respect to the reaction cross-section for the muonic channel in the kinematics limit of the MUONCocktail
    Double_t sigmapion     = 1.80e-2; // Just for reproducing Valeries's data
    Double_t brpion        = 0.9999;  // Branching Ratio for Pion 
    genpion->Init();  // Generating pT and Y parametrsation for the 4pi
    ratiopion = sigmapion * brpion * (0.93*fNumberOfParticipants+0.07*fNumberOfCollisions) / sigmaReaction * genpion->GetRelativeArea(ptMin,ptMax,yMin,yMax,phiMin,phiMax);
    AliInfo(Form("Pseudo-Pion production cross-section in pp with shadowing %5.3g barns",sigmapion));
    AliInfo(Form("Pion production probability per collisions in acceptance via the muonic channel %5.3g",ratiopion));
    genpion->SetPtRange(ptMin, ptMax);  
    genpion->SetYRange(yMin, yMax);
    genpion->SetPhiRange(phiMin, phiMax);
    genpion->Init(); // Generating pT and Y parametrsation in the desired kinematic range
    AddGenerator(genpion,"Pion", ratiopion);
    fTotalRate+=ratiopion;
    
    // Generating Kaon Physics
    AliGenParam * genkaon = new AliGenParam(1, AliGenMUONlib::kKaon, "default", "Kaon");  
    genkaon->SetPtRange(0,100.);  
    genkaon->SetYRange(-8.,8);
    genkaon->SetPhiRange(0.,360.);
    genkaon->SetForceDecay(kKaToMu);
    genkaon->SetTrackingFlag(1);
    Double_t ratiokaon; // Ratio with respect to the reaction cross-section for the muonic channel in the kinematics limit of the MUONCocktail
    Double_t sigmakaon     = 2.40e-4;   // Valerie presentation Clermont-16-jan-2004 and Alice-int-2002-06 
    Double_t brkaon        = 0.6351 ;  // Branching Ratio for Kaon 
    genkaon->Init();  // Generating pT and Y parametrsation for the 4pi
    ratiokaon = sigmakaon * brkaon * (0.93*fNumberOfParticipants+0.07*fNumberOfCollisions)/ sigmaReaction * genkaon->GetRelativeArea(ptMin,ptMax,yMin,yMax,phiMin,phiMax);
    AliInfo(Form("Pseudo-kaon production cross-section in pp with shadowing %5.3g barns",sigmakaon));
    AliInfo(Form("Kaon production probability per collisions in acceptance via the muonic channel %5.3g",ratiokaon));
    genkaon->SetPtRange(ptMin, ptMax);  
    genkaon->SetYRange(yMin, yMax);
    genkaon->SetPhiRange(phiMin, phiMax);
    genkaon->Init(); // Generating pT and Y parametrsation in the desired kinematic range
    AddGenerator(genkaon,"Kaon", ratiokaon);
    fTotalRate+=ratiokaon;
  }
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

//  Generate the vertex position used by all generators    
    if(fVertexSmear == kPerEvent) Vertex();

    // Loop on primordialTrigger
    Bool_t primordialTrigger = kFALSE;
    while(!primordialTrigger) {
		//Reseting stack
		AliRunLoader * runloader = gAlice->GetRunLoader();
		if (runloader)
			if (runloader->Stack())
				runloader->Stack()->Reset();
		// Loop over generators and generate events
		Int_t igen=0;
		Int_t npart =0;
		while((entry = (AliGenCocktailEntry*)next())) {
			gen = entry->Generator();
			gen->SetVertex(fVertex.At(0), fVertex.At(1), fVertex.At(2));
			if ( (npart = gRandom->Poisson(entry->Rate())) >0 ) {
				igen++;	
				if (igen == 1) entry->SetFirst(0);
				else  entry->SetFirst((partArray->GetEntriesFast())+1);
				// if ( (fHadronicMuons == kFALSE) && ( (gen->GetName() == "Pions") || (gen->GetName() == "Kaons") ) )
				//  { AliInfo(Form("This generator %s is finally not generated. This is option for hadronic muons.",gen->GetName() ) ); }
				// else {
				gen->SetNumberParticles(npart);
				gen->Generate();
				entry->SetLast(partArray->GetEntriesFast());
				preventry = entry;
	    		  // }
			}
		}  
		next.Reset();
		// Testing primordial trigger : Muon  pair in the MUON spectrometer acceptance and pTCut
		Int_t iPart;
		fNGenerated++;
		Int_t numberOfMuons=0;
		//      printf(">>>fNGenerated is %d\n",fNGenerated);
		
		TObjArray GoodMuons; // Used in the Invariant Mass selection cut
		
		for(iPart=0; iPart<partArray->GetEntriesFast(); iPart++){      
			
			if ( (TMath::Abs(gAlice->GetMCApp()->Particle(iPart)->GetPdgCode())==13)  &&
	     		(gAlice->GetMCApp()->Particle(iPart)->Theta()*180./TMath::Pi()>fMuonThetaMinCut) &&
	     		(gAlice->GetMCApp()->Particle(iPart)->Theta()*180./TMath::Pi()<fMuonThetaMaxCut) &&
	     		(gAlice->GetMCApp()->Particle(iPart)->Pt()>fMuonPtCut)                             ) { 
	  				gAlice->GetMCApp()->Particle(iPart)->SetProductionVertex(fVertex.At(0), fVertex.At(1), fVertex.At(2), 0.);   
	  				GoodMuons.AddLast(gAlice->GetMCApp()->Particle(iPart));
					numberOfMuons++;
			}			
		}
		
		// Test the invariant mass of each pair (if cut on Invariant mass is required)
		Bool_t InvMassRangeOK = kTRUE;
		if(fInvMassCut && (numberOfMuons>=2) ){
  			TLorentzVector fV1, fV2, fVtot;
			InvMassRangeOK = kFALSE;
			for(iPart=0; iPart<GoodMuons.GetEntriesFast(); iPart++){      
    			TParticle * mu1 = ((TParticle *)GoodMuons.At(iPart));

				for(int iPart2=iPart+1; iPart2<GoodMuons.GetEntriesFast(); iPart2++){      
    					TParticle * mu2 = ((TParticle *)GoodMuons.At(iPart2));
						
						fV1.SetPxPyPzE(mu1->Px() ,mu1->Py() ,mu1->Pz() ,mu1->Energy() );
						fV2.SetPxPyPzE(mu2->Px() ,mu2->Py() ,mu2->Pz() ,mu2->Energy() );
						fVtot = fV1 + fV2;
						
						if(fVtot.M()>fInvMassMinCut && fVtot.M()<fInvMassMaxCut) {
							InvMassRangeOK = kTRUE;
							break;
						}
				}
				if(InvMassRangeOK) break; // Invariant Mass Cut pass as soon as one pair satisfy the criterion
			}	
		}
		
		
		if ((numberOfMuons >= fMuonMultiplicity) &&  InvMassRangeOK ) primordialTrigger = kTRUE;
	}
    fNSucceded++;

    AliDebug(5,Form("Generated Events are %d and Succeeded Events are %d",fNGenerated,fNSucceded));
}
