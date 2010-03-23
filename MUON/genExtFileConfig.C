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
/// \file genExtFileConfig.C
/// \brief Configuration macro for event generator from external file
/// for MUON spectrometer Monte Carlo simulation

// Functions
Float_t EtaToTheta(Float_t arg);
AliGenerator* GeneratorFactory();

void genConfig()
{
  cout << "Running genExtFileConfig.C ... " << endl;

  //=======================================================================
  // ******* GEANT STEERING parameters FOR ALICE SIMULATION *******
  if ( gMC ) {
    gMC->SetProcess("DCAY",1);
    gMC->SetProcess("PAIR",1);
    gMC->SetProcess("COMP",1);
    gMC->SetProcess("PHOT",1);
    gMC->SetProcess("PFIS",0);
    gMC->SetProcess("DRAY",0);
    gMC->SetProcess("ANNI",1);
    gMC->SetProcess("BREM",1);
    gMC->SetProcess("MUNU",1);
    gMC->SetProcess("CKOV",1);
    gMC->SetProcess("HADR",1);
    gMC->SetProcess("LOSS",2);
    gMC->SetProcess("MULS",1);
    gMC->SetProcess("RAYL",1);
    Float_t cut = 1.e-3;        // 1MeV cut by default
    Float_t tofmax = 1.e10;
    gMC->SetCut("CUTGAM", cut);
    gMC->SetCut("CUTELE", cut);
    gMC->SetCut("CUTNEU", cut);
    gMC->SetCut("CUTHAD", cut);
    gMC->SetCut("CUTMUO", cut);
    gMC->SetCut("BCUTE",  cut); 
    gMC->SetCut("BCUTM",  cut); 
    gMC->SetCut("DCUTE",  cut); 
    gMC->SetCut("DCUTM",  cut); 
    gMC->SetCut("PPCUTM", cut);
    gMC->SetCut("TOFMAX", tofmax); 
  }

  //=======================================================================
  // External decayer
  //=======================================================================

  // Set External decayer
  TVirtualMCDecayer *decayer = new AliDecayerPythia();
  decayer->SetForceDecay(kAll);
  decayer->Init();
  if ( gMC ) gMC->SetExternalDecayer(decayer);
  
  //=======================================================================
  // Event generator
  //=======================================================================

  // External generator configuration
  AliGenerator* gener = GeneratorFactory();
  gener->SetOrigin(0, 0, 0);    // vertex position
  //gener->SetSigma(0, 0, 5.3);   // Sigma in (X,Y,Z) (cm) on IP position
  //gener->SetCutVertexZ(1.);     // Truncate at 1 sigma
  //gener->SetVertexSmear(kPerEvent); 
  gener->SetTrackingFlag(1);
  gener->Init();
    
  cout << "Running genExtFileConfig.C finished ... " << endl;
}

Float_t EtaToTheta(Float_t arg){
  return (180./TMath::Pi())*2.*atan(exp(-arg));
}

AliGenerator* GeneratorFactory() {

  AliGenExtFile *gener = new AliGenExtFile(-1);
  AliGenReaderTreeK * reader = new AliGenReaderTreeK();

  reader->SetFileName("galice.root");
  reader->AddDir("$ALICE_ROOT/MUON/gen");
  gener->SetReader(reader);
     
  return gener; 
}

