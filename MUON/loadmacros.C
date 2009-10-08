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
/// \file loadmacros.C
/// \brief Macro which loads and compiles the MUON macros
///
/// \author I. Hrivnacova, IPN Orsay

#if !defined(__CINT__) || defined(__MAKECINT__)

#include <TROOT.h>
#include <TSystem.h>
#include <TString.h>

#endif

void init() 
{
/// Set include path and load libraries which are not 
/// linked with aliroot

  // Redefine include paths as some macros need
  // to see more than what is define in rootlogon.C
  //
  TString includePath = "-I${ALICE_ROOT}/include ";
  includePath        += "-I${ALICE_ROOT}/RAW ";
  includePath        += "-I${ALICE_ROOT}/FASTSIM ";
  includePath        += "-I${ALICE_ROOT}/EVGEN ";
  includePath        += "-I${ALICE_ROOT}/SHUTTLE/TestShuttle ";
  includePath        += "-I${ALICE_ROOT}/ITS ";
  includePath        += "-I${ALICE_ROOT}/MUON ";
  includePath        += "-I${ALICE_ROOT}/MUON/mapping";
  gSystem->SetIncludePath(includePath.Data());

  // Load libraries not linked with aliroot
  //
  gSystem->Load("$ALICE_ROOT/SHUTTLE/TestShuttle/libTestShuttle.so");
  gSystem->Load("libMUONshuttle.so");
  gSystem->Load("libMUONevaluation.so");
  gSystem->Load("liblhapdf.so");     
  gSystem->Load("libpythia6.so");    
  gSystem->Load("libgeant321.so");
  gSystem->Load("libEG");
  gSystem->Load("libEGPythia6");
  gSystem->Load("libAliPythia6.so");
}  

void loadmacro(const TString& macroName)
{
/// Load the macro with given name

  TString path = "$ALICE_ROOT/MUON/";
  path += macroName;
  path += ".C++";
  gROOT->LoadMacro(path.Data());
}  

void loadmacros () 
{
  init();

  loadmacro("DecodeRecoCocktail");         // Hermine, Alessandro     
  loadmacro("fastMUONGen");                // Hermine, Alessandro
  loadmacro("fastMUONSim");                // Hermine, Alessandro
  loadmacro("loadFromOCDB");               // Philippe P.
  loadmacro("MakeMUONFullMisAlignment");   // Javier, Ivana
  loadmacro("MakeMUONResMisAlignment");    // Javier, Ivana
  loadmacro("MakeMUONZeroMisAlignment");   // Javier, Ivana
  loadmacro("MergeMuonLight");             // Hermine, Alessandro
  loadmacro("MUONAlignment");              // Javier   
  loadmacro("MUONCheck");                  // Frederic          
  loadmacro("MUONCheckDI");                // Artur
  loadmacro("MUONCheckMisAligner");        // Javier
  loadmacro("MUONClusterInfo");            // Philippe P.
  loadmacro("MUONefficiency");             // Christophe
  loadmacro("MUONGenerateBusPatch");       // Christian
  loadmacro("MUONGenerateGeometryData");   // Ivana
  loadmacro("MUONGenerateTestGMS");        // Ivana
  loadmacro("MUONGeometryViewingHelper");  // Ivana
  loadmacro("MUONmassPlot_ESD");           // Christian    
  loadmacro("MUONOfflineShift");           // Laurent
  loadmacro("MUONplotefficiency");         // Christian
  loadmacro("MUONRawStreamTracker");       // Christian   
  loadmacro("MUONRawStreamTrigger");       // Christian 
  loadmacro("MUONReCalcGlobalTrigger");    // Bogdan       
  loadmacro("MUONRecoCheck");              // Hermine, Alessandro
  loadmacro("MUONRefit");                  // Philippe P.
  // loadmacro("MUONResoEffChamber");      // Nicolas
  loadmacro("MUONStatusMap");              // Laurent
  loadmacro("MUONSurveyUtil");             // Javier
  loadmacro("MUONSurveyCh8L");             // Javier
  loadmacro("MUONTimeRawStreamTracker");   // Artur
  loadmacro("MUONTimeRawStreamTrigger");   // Artur 
  loadmacro("MUONTrigger");                // Bogdan
  loadmacro("MUONTriggerChamberEfficiency"); // Diego
  loadmacro("MUONTriggerEfficiency");      // Bogdan
  loadmacro("MUONTriggerEfficiencyPt");    // Bogdan
  loadmacro("ReadRecoCocktail");           // Hermine, Alessandro   
  loadmacro("runReconstruction");          // Laurent
  loadmacro("runSimulation");              // Laurent
  loadmacro("TestMUONPreprocessor");       // Laurent
  loadmacro("TestRecPoints");              // Diego
  loadmacro("TestRecPoints");              // Diego
}
