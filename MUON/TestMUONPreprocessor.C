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

#if !defined(__CINT__) || defined(__MAKECINT__)

// This macro runs the test preprocessor
// It uses AliTestShuttle to simulate a full Shuttle process

// The input data is created in the functions
//   CreateDCSAliasMap() creates input that would in the same way come from DCS
//   ReadDCSAliasMap() reads from a file
//   CreateInputFilesMap() creates a list of local files, that can be accessed by the shuttle
//
// According to SHUTTLE/TestShuttle/TestPreprocessor.C
// By Laurent Aphecetche, SUBATECH Nantes

#include "AliCDBManager.h"
#include "AliShuttleInterface.h"
#include "AliCDBId.h"
#include "AliTestShuttle.h"
#include "TMap.h"
#include "Riostream.h"
#include "TSystem.h"
#include "AliMpExMap.h"
#include "TMap.h"
#include "TString.h"
#include "TObjArray.h"
#include "AliMpHelper.h"
#include "AliDCSValue.h"
#include "TObjString.h"
#include "TRandom.h"
#include "AliMUONPreprocessor.h"
#include "AliCDBEntry.h"
#endif

void TestMUONPreprocessor()
{
  // load library
  gSystem->Load("../SHUTTLE/TestShuttle/libTestShuttle.so");
  gSystem->Load("libMUONshuttle.so");
  
  // initialize location of CDB
  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/SHUTTLE/TestShuttle/TestCDB");

  // create AliTestShuttle instance
  // The parameters are run, startTime, endTime
  AliTestShuttle* shuttle = new AliTestShuttle(800, 0, 1);

  // TODO(1)
  //
  // The shuttle can read DCS data, if the preprocessor should be tested to process DCS data,
  // some fake data has to be created.
  //
  // The "fake" input data can be taken using either (a) or (b):
  // (a) data from a file: Use ReadDCSAliasMap()
  //     the format of the file is explained in ReadDCSAliasMap()
  //     To use it uncomment the following line:
  //
// TMap* dcsAliasMap = ReadDCSAliasMap();
  //
  // (b) generated in this macro: Use CreateDCSAliasMap() and its documentation
  //     To use it uncomment the following line:
  //
// TMap* dcsAliasMap = CreateDCSAliasMap();

  // now give the alias map to the shuttle
  // shuttle->SetDCSInput(dcsAliasMap);

  // TODO(2)
  //
  // The shuttle can also process files that originate from DCS, DAQ and HLT.
  // To test it, we provide some local files and locations where these would be found when
  // the online machinery would be there.
  // In real life this functions would be produces by the sub-detectors
  // calibration programs in DCS, DAQ or HLT. These files can then be retrieved using the Shuttle.
  //
  // Files are added with the function AliTestShuttle::AddInputFile. The syntax is:
  // AddInputFile(<system>, <detector>, <id>, <source>, <local-file>)
  // In this example we add a file originating from the GDC with the id PEDESTALS
  // Three files originating from different LDCs but with the same id are also added
//  shuttle->AddInputFile(AliTestShuttle::kDAQ, "DET", "PEDESTALS", "GDC", "file1.root");
//  shuttle->AddInputFile(AliTestShuttle::kDAQ, "DET", "DRIFTVELOCITY", "LDC0", "file2a.root");
//  shuttle->AddInputFile(AliTestShuttle::kDAQ, "DET", "DRIFTVELOCITY", "LDC1", "file2b.root");
//  shuttle->AddInputFile(AliTestShuttle::kDAQ, "DET", "DRIFTVELOCITY", "LDC2", "file2b.root");

  shuttle->AddInputFile(AliTestShuttle::kDAQ,"MCH","PEDESTALS","LDC0","$ALICE_ROOT/MUON/data/LDC0.ped");
  shuttle->AddInputFile(AliTestShuttle::kDAQ,"MCH","PEDESTALS","LDC1","$ALICE_ROOT/MUON/data/LDC1.ped");
  shuttle->AddInputFile(AliTestShuttle::kDAQ,"MCH","PEDESTALS","LDC2","$ALICE_ROOT/MUON/data/LDC2.ped");
  shuttle->AddInputFile(AliTestShuttle::kDAQ,"MCH","PEDESTALS","LDC3","$ALICE_ROOT/MUON/data/LDC3.ped");
  
  shuttle->AddInputFile(AliTestShuttle::kDAQ,"MCH","GMS","GMS","$ALICE_ROOT/MUON/data/GMS.root");

  // TODO(3)
  // Create the preprocessor that should be tested, it registers itself automatically to the shuttle
//  AliPreprocessor* pp = new AliTestPreprocessor("DET", shuttle);
  new AliMUONPreprocessor("MCH", shuttle);

  shuttle->Print();
  
  // Test the preprocessor
  shuttle->Process();

  // TODO(4)
  // In the preprocessor AliShuttleInterface::Store should be called to put the final
  // data to the CDB. To check if all went fine have a look at the files produced in
  // $ALICE_ROOT/SHUTTLE/TestShuttle/TestCDB/<detector>/SHUTTLE/Data
  //
  // Check the file which should have been created
//  AliCDBEntry* entry = AliCDBManager::Instance()->Get("DET/SHUTTLE/Data", 7);
//  if (!entry)
//  {
//    printf("The file is not there. Something went wrong.\n");
//    return;
//  }
//
//  AliTestDataDCS* output = dynamic_cast<AliTestDataDCS*> (entry->GetObject());
//  // If everything went fine, draw the result
//  if (output)
//    output->Draw();
}
