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

/// \ingroup macros
/// \file runShuttleOffline
/// \brief Macro to run the MCH Shuttle offline
///
/// This macro is intended to be used in cases where the online Shuttle
/// fails (for whatever reason, but hopefully in very rare cases...).
///
/// - First step is to retrieve the input files for the Shuttle
/// - Then must run them through this macro
/// - Then must check them (using mchview program for most of them but BPEVO
/// which is tested using the MUONBusPatchEvolution.C macro)
/// - Finally upload the resulting root files to the OCDB
///
/// Given it's extremely rare we have to use this macro, it's a bit manual,
/// e.g. you have to specify the input files and timestamps by hand, etc...
///
/// \author Laurent Aphecetche, SUBATECH Nantes

#include "AliCDBManager.h"
#include "AliCDBRunRange.h"
#include "AliMUON2DMap.h"
#include "AliMUONCDB.h"
#include "AliMUONTrackerIO.h"
#include "AliMUONTrackerPreprocessor.h"
#include "AliMpCDB.h"
#include "AliMpDDLStore.h"
#include "AliShuttleInterface.h"
#include "AliTestShuttle.h"
#include "TMap.h"
#include "TString.h"
#include <iostream>

//______________________________________________________________________________
TString fileName(const char* dir, int runNumber, const char* da, int i, const char* type)
{
  const char* format = "%s/run%09d_MCH_%s-%d_%s";
  return Form(format,dir,runNumber,da,i,type);
}

//______________________________________________________________________________
void runMCHShuttleOffline(int runNumber, const char* runType,
        const char* sourceDirectory,
        const char* dcsmapfile,
        int startDCStime,
        int endDCStime) 
{
  AliTestShuttle shuttle(runNumber, startDCStime,endDCStime);

  const char* inputCDB = "local://$ALIROOT_OCDB_ROOT/OCDB";

  AliTestShuttle::SetMainCDB(inputCDB);
  AliTestShuttle::SetMainRefStorage("local://$ALICE_ROOT/OCDB/SHUTTLE/TestShuttle/TestReference");

  printf("Test Shuttle temp dir: %s\n", AliShuttleInterface::GetShuttleTempDir());
  printf("Test Shuttle log dir: %s\n", AliShuttleInterface::GetShuttleLogDir());
  printf("Test OCDB storage Uri: %s\n", AliShuttleInterface::GetMainCDB().Data());
  printf("Test Reference storage Uri: %s\n", AliShuttleInterface::GetMainRefStorage().Data());

  // The shuttle can process files that originate from DCS, DAQ and HLT.
  // To test it, we provide some local files and locations where these would be found when
  // the online machinery would be there.
  // In real life this functions would be produces by the sub-detectors
  // calibration programs in DCS, DAQ or HLT. These files can then be retrieved using the Shuttle.
  //
  // Files are added with the function AliTestShuttle::AddInputFile. The syntax is:
  // AddInputFile(<system>, <detector>, <id>, <source>, <local-file>)
  // In this example we add 4 files originating from different LDCs but with the same id (PEDESTALS)

  if (TString(runType)=="PEDESTAL") {

      TString da = "ldc-MUON_TRK";
      for ( int i = 0; i < 7; ++i ) {
          shuttle.AddInputFile(AliTestShuttle::kDAQ,"MCH","PEDESTALS",Form("LDC%d",i),fileName(sourceDirectory,runNumber,da.Data(),i,"PEDESTALS"));
          shuttle.AddInputFile(AliTestShuttle::kDAQ,"MCH","CONFIG",Form("LDC%d",i),fileName(sourceDirectory,runNumber,da.Data(),i,"CONFIG"));

      }
  } else if (TString(runType)=="PHYSICS") {
      TString da = "mon-DA-MCH";
      shuttle.AddInputFile(AliTestShuttle::kDAQ,"MCH","OCCUPANCY","mon",fileName(sourceDirectory,runNumber,da.Data(),0,"OCCUPANCY"));
      shuttle.AddInputFile(AliTestShuttle::kDAQ,"MCH","BPEVO","mon",fileName(sourceDirectory,runNumber,da.Data(),0,"BPEVO"));
      TFile* f = TFile::Open(dcsmapfile);
      TMap* dcsAliasMap = static_cast<TMap*>(f->Get("DCSmap")->Clone());
      delete f;
      shuttle.SetDCSInput(dcsAliasMap);
  }
  else {
      std::cout << "Invalid runtype " << runType << std::endl;
  }
  shuttle.SetInputRunType(runType);

  shuttle.AddInputRunParameter("totalEvents","20");

  // Create the preprocessor that should be tested, it registers itself automatically to the shuttle
  new AliMUONTrackerPreprocessor(&shuttle);

  shuttle.Print();

  // Test the preprocessor
  shuttle.Process();

}

//______________________________________________________________________________
void runMCHShuttleOffline() {

  runMCHShuttleOffline(272781,"PEDESTAL","$HOME/analysis/2017/LHC17h/MCH_272781/","",0,1);

  runMCHShuttleOffline(272763,"PHYSICS","$HOME/analysis/2017/LHC17h/MCH_272763/","$HOME/analysis/2017/LHC17h/MCH_272763/testDCSMap.root_MCH_1498419754_1498425527_run272763.root",1498419754,1498425527);

  runMCHShuttleOffline(272764,"PHYSICS","$HOME/analysis/2017/LHC17h/MCH_272764/","$HOME/analysis/2017/LHC17h/MCH_272764/testDCSMap.root_MCH_1498425548_1498428314_run272764.root",1498425548,1498428314);

}

