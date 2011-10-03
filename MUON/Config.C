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
/// \file Config.C
/// \brief Configuration macro for MUON spectormeter simulation
///
/// Remember to define the directory and event generator option:
///
/// gAlice->SetConfigFunction("Config('$HOME','box');");
///
/// Can be compiled by :
/// .L Config.C+
/// (depends on the correct rootlogon.C though...)

#if !defined(__CINT__) || defined(__MAKECINT__)

#include "AliABSOv3.h"
#include "AliBODY.h"
#include "AliConfig.h"
#include "AliDIPOv3.h"
#include "AliDecayerPythia.h"
#include "AliFRAMEv2.h"
#include "AliGenBox.h"
#include "AliGenFixed.h"
#include "AliGenHijing.h"
#include "AliGenMUONCocktail.h"
#include "AliGenMUONlib.h"
#include "AliGenParam.h"
#include "AliGenScan.h"
#include "AliHALLv3.h"
#include "AliITSv11Hybrid.h"
#include "AliMAG.h"
#include "AliMagF.h"
#include "AliMUONv1.h"
#include "AliPIPEv3.h"
#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliSHILv3.h"
#include <Riostream.h>
#include <TDatime.h>
#include <TGeant3TGeo.h>
#include <TGeoGlobalMagField.h>
#include <TPDGCode.h>
#include <TRandom.h>
#include <TSystem.h>
#include <TVirtualMC.h>

//#include "AliZDCv3.h"
//#include "AliFMDv1.h"
//#include "AliPMDv1.h"
//#include "AliVZEROv7.h"

#endif

void Config(const char* directory="", 
            const char* option="param", 
            const char* digitstore="AliMUONDigitStoreV2S",
            bool forEmbedding=kFALSE,
            bool forRealistic=kFALSE)
{
  //=====================================================================
  // Config file for MUON test
  //=====================================================================
  //  Libraries required by geant321
  
  gSystem->Load("liblhapdf.so");      // Parton density functions
  gSystem->Load("libpythia6.so");     // Pythia
  gSystem->Load("libgeant321.so");
  gSystem->Load("libEG");
  gSystem->Load("libEGPythia6");
  gSystem->Load("libAliPythia6.so");  // ALICE specific implementations

  new TGeant3TGeo("C++ Interface to Geant3");

  //=======================================================================
  //  Create the output file    
  Text_t filename[100];
  sprintf(filename,"%sgalice.root",directory);
  cout << ">>> Output file is " << filename << endl;   
  cout << ">>> Config.C: Creating Run Loader ..."<<endl;
  AliRunLoader* rl=0x0;
  rl = AliRunLoader::Open(
			  filename, AliConfig::GetDefaultEventFolderName(), "recreate");
  if (rl == 0x0) 
    { gAlice->Fatal("Config.C","Can not instatiate the Run Loader");
      return; }
  rl->SetCompressionLevel(2);
  rl->SetNumberOfEventsPerFile(10000);
  gAlice->SetRunLoader(rl);
  //=======================================================================
  // For having more debuging messages
  //AliLog::SetModuleDebugLevel("MUON", 1);
  //=======================================================================
  // Set External decayer
  TVirtualMCDecayer *decayer = new AliDecayerPythia();
  decayer->SetForceDecay(kAll);
  decayer->Init();
  gMC->SetExternalDecayer(decayer);
  //=======================================================================
  // ******* GEANT STEERING parameters FOR ALICE SIMULATION *******
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
    gener->SetOrigin(0,0,0);          //vertex position    gener->SetSigma(0,0,0); //Sigma in (X,Y,Z) (cm) on IP position
    if (forRealistic)
    {
      gener->SetVertexSmear(kPerEvent);
    }
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
  //============================================================= 
  // Field (L3 0.5 T)
  TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps", -1., -1, AliMagF::k5kG));
  //============================================================= 
  //=================== Alice BODY parameters =============================
  new AliBODY("BODY","Alice envelop");
  //=================== ABSO parameters ============================
  new AliABSOv3("ABSO", "Muon Absorber");
  //=================== DIPO parameters ============================
  new AliDIPOv3("DIPO", "Dipole version 2");
  //================== HALL parameters ============================
  new AliHALLv3("HALL", "Alice Hall");
  //=================== PIPE parameters ============================
  new AliPIPEv3("PIPE", "Beam Pipe");
  //=================== SHIL parameters ============================
  new AliSHILv3("SHIL", "Shielding Version 2");

  //=================== MUON Subsystem ===========================
  AliMUON *MUON = new AliMUONv1("MUON", "default");

  // The 3 switches below are to be used for the trigger code
  // their default value is set in AliMUON.h
  // activate trigger cluster-size (0=default, 1=cluster-size according to AliMUONResponseTriggerV1
  //  MUON->SetTriggerResponseV1(0);
  // activate 4/4 trigger coincidence (0=default (coinc 3/4), 1=coinc 4/4)
  //  MUON->SetTriggerCoinc44(0);
  // activate trigger chamber efficiency by cells (0=default, 1=trigger efficiency according to AliMUONTriggerEfficiencyCells
  //  MUON->SetTriggerEffCells(0);

  // Activate the following line when running with Fluka
  // MUON->SetIsMaxStep(kFALSE);

  // Use SetDigitStoreClassName() to change the digitStore implementation used by (s)digitizer
  MUON->SetDigitStoreClassName(digitstore);

  if ( forEmbedding ) 
  {
    // Noise-only digits in tracker/trigger (0=no noise, 1=default (noise in tracker), 2=noise in tracker and trigger):
    cout << "****** DISABLING NOISE GENERATION AS WE DO EMBEDDING ******" << endl;
    MUON->SetDigitizerWithNoise(0);
		MUON->SetConvertTrigger(true);
    
    new AliITSv11Hybrid("ITS","ITS v11Hybrid");    
  }
  else if (forRealistic)
  {
    new AliITSv11Hybrid("ITS","ITS v11Hybrid");        
  }
  
  // Use non-high performance raw data decoder 
  //MUON->SetFastTrackerDecoder(kFALSE);  
  //MUON->SetFastTriggerDecoder(kFALSE);  
  
  //
  // If SetAlign, the detection elements transformations
  // are taken from the input file and not from the code
  // MUON->SetAlign("transform.dat");

  // To generate and read scaler trigger events in rawdata
  // MUON->SetTriggerScalerEvent();
  
  // To switch off the tail effect
  // MUON->SetTailEffect(kFALSE);

  // If you want to play with builders, first reset the geometry builder,
  // and then add yours.
  //  MUON->ResetGeometryBuilder();
  //  MUON->AddGeometryBuilder(new AliMUONSt1GeometryBuilderV2(MUON));
  //  MUON->AddGeometryBuilder(new AliMUONSt2GeometryBuilderV2(MUON));
  //  MUON->AddGeometryBuilder(new AliMUONSlatGeometryBuilder(MUON));
  //  MUON->AddGeometryBuilder(new AliMUONTriggerGeometryBuilder(MUON));
}
