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

// $Id$
//
/// \ingroup macros
/// \file MUONGenerateGeometryData.C
/// \brief Macro for generating the geometry data files:
///  transform.dat, svmap.dat.
///
/// To be run from aliroot:
///
/// .x MUONGenerateGeometryData.C
///
/// The generated files do not replace the existing ones
/// but have different names (with extension ".out").
///
/// \author: I. Hrivnacova, IPN Orsay

#if !defined(__CINT__) || defined(__MAKECINT__)

#include "AliMUON.h"
#include "AliMUONGeometryBuilder.h"
#include "AliMUONGeometryTransformer.h"

#include "AliRun.h"
#include "AliCDBManager.h"
#include "AliMC.h"

#include <Riostream.h>
#include <TROOT.h>
#include <TInterpreter.h>

#endif

void MUONGenerateGeometryData(Bool_t transforms = true, 
                              Bool_t svmaps = true,
                              Bool_t writeEnvelopes = true)
{
/// \param transforms      option to generete transform.dat
/// \param svmaps          option to generete svmap.dat
/// \param writeEnvelope   option to include virtual envelopes
///                        in the volume paths

  // Default CDB and run number
  AliCDBManager* man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE_ROOT");
  man->SetRun(0);

  // Initialize
  TString configFileName = "$ALICE_ROOT/MUON/Config.C";
  gROOT->LoadMacro(configFileName.Data());
  gInterpreter->ProcessLine(gAlice->GetConfigFunction());
  gAlice->GetMCApp()->Init();
  cout << "Init done " << endl;

  // Get MUON detector
  AliMUON* muon = (AliMUON*)gAlice->GetModule("MUON");
  if (!muon) {
    cerr << "MUON detector not defined." << endl;
    return;
  }  

  // Get geometry builder
  AliMUONGeometryBuilder* builder = muon ->GetGeometryBuilder();

  if (transforms) {
    cout << "Generating transformation file ..." << endl;
    builder->GetTransformer()->WriteTransformations("transform.dat.out");
  }  

  if (svmaps) {
    cout << "Generating svmaps file ..." << endl;
    builder->WriteSVMaps("svmap.dat.out", true, writeEnvelopes);
  }  
}  
