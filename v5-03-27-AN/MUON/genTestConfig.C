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
/// \file genTestConfig.C
/// \brief Configuration macro for event generator 
/// for MUON spectrometer Monte Carlo simulation

Float_t EtaToTheta(Float_t arg);

AliGenerator* genConfig(char option[6]="param")
{
  cout << "Running genTestConfig.C ... " << endl;

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

  // Set External decayer
  TVirtualMCDecayer *decayer = new AliDecayerPythia();
  decayer->SetForceDecay(kAll);
  decayer->Init();
  if ( gMC ) gMC->SetExternalDecayer(decayer);
  
  //=======================================================================
  // Examples of generators. Only option param is sistematically tested
  if (!strcmp(option,"box")) {
    AliGenBox * gener = new AliGenBox(1);
    gener->SetMomentumRange(20.,20.1);
    gener->SetPhiRange(0., 360.);         
    gener->SetThetaRange(171.000,178.001);
    gener->SetPart(kMuonMinus);           // Muons
    gener->SetOrigin(0.,0., 0.);  //vertex position
    gener->SetSigma(0.0, 0.0, 0.0);         //Sigma in (X,Y,Z) (cm) on IP position
  }
  if (!strcmp(option,"gun")) {
    AliGenFixed *gener = new AliGenFixed(1);
    gener->SetMomentum(10);
    gener->SetPhiRange(0.);
    gener->SetThetaRange(0.);
    gener->SetOrigin(30,30,-1200);//vertex position
    gener->SetPart(kMuonMinus);          //GEANT particle type  13 is muons
  }
  if (!strcmp(option,"scan")) {
    AliGenScan *gener = new AliGenScan(-1);
    gener->SetMomentumRange(10,10);
    gener->SetPhiRange(0, 0);
    gener->SetThetaRange(-180, -180);
    //vertex position
    //gener->SetSigma(1,1,0);           //Sigma in (X,Y,Z) (cm) on IP position
    gener->SetPart(kRootino); 
    gener->SetRange(100, -300., 300., 100, -300., 300., 1, 2000, 2000);
  }  
  if (!strcmp(option,"param")) {
    AliGenParam *gener = new AliGenParam(1, AliGenMUONlib::kUpsilon);
    gener->SetMomentumRange(0,999);
    gener->SetPtRange(0,100.);
    gener->SetPhiRange(0., 360.);
    gener->SetCutOnChild(1);
    gener->SetChildPhiRange(0.,360.);
    gener->SetChildThetaRange(171.0,178.0);
    gener->SetOrigin(0,0,0);          //vertex position    gener->SetSigma(0,0,0);           //Sigma in (X,Y,Z) (cm) on IP position
    gener->SetForceDecay(kDiMuon);
    gener->SetTrackingFlag(1);
    gener->Init();
  }
  if (!strcmp(option,"paramJpsi")) {
    AliGenParam *gener = new AliGenParam(1, AliGenMUONlib::kJpsi);
    gener->SetMomentumRange(0,999);
    gener->SetPtRange(0,100.);
    gener->SetPhiRange(0., 360.);
    gener->SetCutOnChild(1);
    gener->SetChildPhiRange(0.,360.);
    gener->SetChildThetaRange(171.0,178.0);
    gener->SetOrigin(0,0,0);
    gener->SetForceDecay(kDiMuon);
    gener->SetTrackingFlag(1);
    gener->Init();
  }
  if (!strcmp(option,"hijing")) { //Hijing generator from ConfigPPR in macros
    AliGenHijing *gener = new AliGenHijing(-1);
    // centre of mass energy 
    gener->SetEnergyCMS(5500.);
    // reference frame
    gener->SetReferenceFrame("CMS");
    // projectile
    gener->SetProjectile("A", 208, 82);
    gener->SetTarget    ("A", 208, 82);
    // tell hijing to keep the full parent child chain
    gener->KeepFullEvent();
    // enable jet quenching
    gener->SetJetQuenching(1);
    // enable shadowing
    gener->SetShadowing(1);
    // neutral pion and heavy particle decays switched off
    gener->SetDecaysOff(1);
    // Don't track spectators
    gener->SetSpectators(0);
    // kinematic selection
    gener->SetSelectAll(0);
    // impact parameter range
    gener->SetImpactParameterRange(0., 5.); // 0. - 5. fm corresponds to ~10% most central
    gener->Init();
  }
  if (!strcmp(option,"muoncocktail")) { // Muon cocktail for PbPb
    AliGenMUONCocktail * gener = new AliGenMUONCocktail();
    gener->SetPtRange(1.,100.);       // Transverse momentum range  
    gener->SetPhiRange(0.,360.);    // Azimuthal angle range 
    gener->SetYRange(-4.0,-2.5);
    gener->SetMuonPtCut(0.5);
    gener->SetMuonThetaCut(171.,178.);
    gener->SetMuonMultiplicity(2);
    gener->SetImpactParameterRange(0.,5.); // 10% most centra PbPb collisions
    gener->SetVertexSmear(kPerTrack);  
    gener->SetOrigin(0,0,0);        // Vertex position
    gener->SetSigma(0,0,0.0);       // Sigma in (X,Y,Z) (cm) on IP position
    gener->Init();
  }  
  return gener;
  
  cout << "Running genGunConfig.C finished ... " << endl;
}

Float_t EtaToTheta(Float_t arg){
  return (180./TMath::Pi())*2.*atan(exp(-arg));
}
