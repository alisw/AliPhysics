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

/// \ingroup macros
/// \file fastMUONGen.C
/// \brief An example how to use the AliGenMUONCocktailpp generator
/// without GEANT simulation of the detector
///
/// \author H. Woehri,  A. De Falco, INFN Cagliari, April 2007
///
/// The macro switches on all decay modes for the resonances,
/// while for the minimum bias Pythia event we additionally
/// switch on the muonic decays of pions and Kaons
/// Its outcome can be further processed by the macro
/// fastMUONSim.C

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliGenerator.h"
#include "AliRunLoader.h"
#include "AliRun.h"
#include "AliHeader.h"
#include "AliStack.h"
#include "AliGenCocktail.h"
#include "AliDecayerPythia.h"
#include "AliPDG.h"
#include "AliGenMUONCocktailpp.h"

#include "TParticle.h"
#endif

AliGenerator*  CreateGeneratorMC(Int_t mult);
Int_t SetupOutputDirectory();

/// The third argument \em mult allows to select at a "pre-trigger" level
/// events with "mult" muons in the MUON detector's phase space
/// window, which is defined in the method "CreateGeneratorMC(Int_t mult)".
/// Note that in this routine also a cut on the muon's origin
/// should be placed in order not to trigger on muons from pi/K decays
/// that were decayed by Pythia way inside the absorber or in the muon
/// spectrometer itself
void fastMUONGen(Int_t nev = 1, char* filename = "galice.root", Int_t mult = 2)
{
  Int_t runNumber = SetupOutputDirectory();
  gAlice->SetRunNumber(runNumber);
  printf("\n\n\n\nsetting run number to %d\n\n\n\n", runNumber);

  //  Update data base
  AliPDG::AddParticlesToPdgDataBase();
  //
  AliRunLoader* rl = AliRunLoader::Open(filename,"FASTRUN","recreate");
  rl->SetCompressionLevel(2);
  rl->SetNumberOfEventsPerFile(10000);
  rl->LoadKinematics("RECREATE");
  rl->MakeTree("E");
  gAlice->SetRunLoader(rl);
  //
  rl->MakeStack();
  AliStack* stack      = rl->Stack();
  AliHeader* header = rl->GetHeader();
  //
  //  Create and Initialize Generator
  AliGenerator *gener = CreateGeneratorMC(mult);
  gener->SetStack(stack);
  gener->Init();
  //
  for (int iev = 0; iev < nev; iev++) {

    if(iev%1000 == 0) printf("Event %d\n", iev);
    //  Initialize event
    header->Reset(0,iev);
    rl->SetEventNumber(iev);
    stack->Reset();
    rl->MakeTree("K");
    
    gener->Generate();
/*
    // Printing
    Int_t npart = stack->GetNprimary();
    for (Int_t part=0; part<npart; part++) {
      TParticle *MPart = stack->Particle(part);
      Int_t mpart  = MPart->GetPdgCode();
      printf("Particle %5d %5d  %5d\n", part, mpart, MPart->GetFirstMother());
    }
*/
      
    //  Finish event
    header->SetNprimary(stack->GetNprimary());
    header->SetNtrack(stack->GetNtrack());  
    //      I/O
    // 
    stack->FinishEvent();
    header->SetStack(stack);
    rl->TreeE()->Fill();
    rl->WriteKinematics("OVERWRITE");
  } // event loop
  
  gener->FinishRun();
  rl->WriteHeader("OVERWRITE");
  gener->Write();
  rl->Write();  
}
//===============================================
AliGenerator*  CreateGeneratorMC(Int_t mult)
{
    AliGenMUONCocktailpp* gener = new AliGenMUONCocktailpp();
    gener->SetPtRange(0.,100.);
    gener->SetYRange(-4.,-2.4);
    gener->SetPhiRange(0., 360.);
    gener->SetMuonMultiplicity(mult);
    gener->SetMuonPtCut(0.5);
    gener->SetMuonThetaRange(171.,178.);      
    gener->SetOrigin(0.,0.,0.); 
    gener->SetSigma(0.,0.,5.);
    gener->SetVertexSmear(kPerEvent);
    //for resonances: "kAll" all decay modes are openend
    gener->SetDecayModeResonance(kAll);  
    //for MC Pythia events: all decays, including muonic decays from pion and Kaons are opened
    gener->SetDecayModePythia(kAllMuonic); 
    gener->SetMuonOriginCut(-130.); //prevent the trigger muon(s) from decaying after 130 cm
    AliDecayerPythia* decayer = new AliDecayerPythia();
    gener->SetDecayer(decayer);
    return gener;
}
//===============================================
Int_t SetupOutputDirectory(){

  //Setting up the name for the output directory:
  static Int_t sseed = 0; // Set 0 to use the current time
  gRandom->SetSeed(sseed);
  UInt_t theSeed = gRandom->GetSeed();
  Int_t labIndex = 5; //1 subatech, 2 clermont, 3 Torino, 4 Orsay
                      // 5 Cagliari, etc...
  Int_t runNumber = (theSeed%100000000 + 100000000*labIndex);

  Char_t name[100];
  sprintf(name, "touch run_%d", runNumber);
  gSystem->Exec(name);

  return runNumber;
}
