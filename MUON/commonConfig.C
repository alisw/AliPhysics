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
/// \file commonConfig.C
/// \brief Configuration macro 
/// for MUON spectrometer Monte Carlo simulation

enum PprTrigConf_t {
    kDefaultPPTrig, kDefaultPbPbTrig
};

const char * pprTrigConfName[] = {
    "p-p","Pb-Pb"
};

// Options 
static AliMagF::BMap_t smag = AliMagF::k5kG;
static PprTrigConf_t strig = kDefaultPPTrig; // default PP trigger configuration
static TString comment;

// Functions
void  LoadPythia();

void commonConfig(Bool_t setRootGeometry = kFALSE,
                  const char* digitstore="AliMUONDigitStoreV2S")
{
  cout << "Running commonConfig.C ... " << endl;

  //=======================================================================
  // Load Pythia libraries
  //=======================================================================

  LoadPythia();

  //=======================================================================
  // ALICE steering object (AliRunLoader)
  //=======================================================================

  // Set Root geometry file
  if ( setRootGeometry ) {
    AliSimulation::Instance()->SetGeometryFile("$ALICE_ROOT/MUON/geometry/geometry.root");
  }

  AliRunLoader* rl 
    = AliRunLoader::Open("galice.root",
			  AliConfig::GetDefaultEventFolderName(),
			  "recreate");
  if ( ! rl ) {
    gAlice->Fatal("Config.C","Can not instatiate the Run Loader");
    return;
  }
  rl->SetCompressionLevel(2);
  rl->SetNumberOfEventsPerFile(100);
  gAlice->SetRunLoader(rl);
  
  //======================================================================
  // Trigger configuration
  //=======================================================================

  gAlice->SetTriggerDescriptor(pprTrigConfName[strig]);
  cout << "Trigger configuration is set to  " << pprTrigConfName[strig] << endl;

  // ============================= 
  // Magnetic field
  // ============================= 

  // Field (L3 0.5 T)
  // Field (L3 0.5 T)
  TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps", -1., -1, AliMagF::k5kG));
 
  //============================================================= 
  //=================== Alice BODY parameters =============================
  AliBODY *BODY = new AliBODY("BODY","Alice envelop");
  //=================== ABSO parameters ============================
  AliABSO *ABSO = new AliABSOv3("ABSO", "Muon Absorber");
  //=================== DIPO parameters ============================
  AliDIPO *DIPO = new AliDIPOv3("DIPO", "Dipole version 2");
  //================== HALL parameters ============================
  AliHALL *HALL = new AliHALLv3("HALL", "Alice Hall");
  //=================== PIPE parameters ============================
  AliPIPE *PIPE = new AliPIPEv3("PIPE", "Beam Pipe");
  //=================== SHIL parameters ============================
  AliSHIL *SHIL = new AliSHILv3("SHIL", "Shielding Version 2");

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

  // Use SetDigitStoreClassName() to change the digitStore implementation used by (s)digitizer
  MUON->SetDigitStoreClassName(digitstore);
  
  cout << "MUON DigitStore is " << MUON->DigitStoreClassName().Data() << endl;
  
  // Noise-only digits in tracker/trigger (0=no noise, 1=default (noise in tracker), 2=noise in tracker and trigger):
  //MUON->SetDigitizerWithNoise(kFALSE);

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

  cout << "Running commonConfig.C finished ... " << endl;
}

void LoadPythia()
{
  // Load Pythia related libraries
  gSystem->Load("liblhapdf.so");      // Parton density functions
  gSystem->Load("libEGPythia6.so");   // TGenerator interface
  gSystem->Load("libpythia6.so");     // Pythia
  gSystem->Load("libAliPythia6.so");  // ALICE specific implementations
}
