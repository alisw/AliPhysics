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

// Macro which loads and compiles the MUON macros:
//
// runSimulation.C               - ok, comp,  x;      Laurent
// runReconstruction.C           - ok, comp,  x;      Laurent
// fastMuonGen.C                 - ok, comp,  x;      Hermine, Alessandro
// fastMuonSim.C                 - ok, comp,  x;      Hermine, Alessandro
// DecodeRecoCocktail.C          - ok, comp,  README; Hermine, Alessandro
// ReadRecoCocktail.C            - ok, comp,  README; Hermine, Alessandro
// MergeMuonLight.C              - x,  comp,  README; Hermine, Alessandro
// MakeMUONFullMisAlignment.C    - ok, comp,  README; Javier, Ivana
// MakeMUONResMisAlignment.C     - ok, comp,  README; Javier, Ivana
// MakeMUONZeroMisAlignment.C    - ok, comp,  README; Javier, Ivana
// MUONAlignment.C               - ok, comp,  README; Javier
// MUONCheck.C                   - ok, comp,  x,      Frederic, in test
// MUONCheckDI.C                 - x,  comp,  x       Artur     
// MUONCheckMisAligner.C         - ok, comp,  x,      Javier
// MUONdisplay.C                 - ok, comp,  deprecated
// MUONefficiency.C              - ok, comp,  README; Christophe, in test
// MUONGenerateBusPatch.C        - ok, comp,  x,      Christian
// MUONGenerateGeometryData.C    - ok, comp,  README; Ivana
// MUONGenerateTestGMS.C         - ok, comp,  READMEshuttle;  Ivana
// MUONmassPlot_ESD.C            - ok, comp,  README, Christian
// MUONplotefficiency.C          - ok, comp,  README; Christophe
// MUONRawStreamTracker.C        - x,  comp,  README; Christian
// MUONRawStreamTrigger.C        - x,  comp,  README; Christian
// MUONRecoCheck.C               - ok, comp,  README; Hermine, Alessandro
// MUONResoEffChamber.C          - ok, comp,  x,      Nicolas
// MUONStatusMap.C               - ok, comp,  x,      Laurent
// MUONTrigger.C                 - ok, comp,  README, Philippe C.
// MUONTriggerEfficiency.C       - ok, comp,  x,      Philippe C., in test
// MUONTriggerEfficiencyPt.C     - x,  comp,  README, Philippe C.
// MUONTriggerChamberEfficiency.C- x,  comp,  README, Diego
// TestMUONPreprocessor.C        - ok, comp,  READMEshuttle; Laurent
// 
// 1st item:
// ok/fails/x - if the macro runs; x means that it was not tried either for lack
//              of documentation, or expertise
//              
// 2nd item:
// comp/x     - if the macro can be compiled
//              
// 3rd item:
// README*/x  - if it is documented in README, x means no doxumentation outside the 
//              macro itself
//              
// 4th item:  - author(s), responsible for macro maintenance 
//
// eventually 
// 5th item:  - if the macro is run within test scripts                         
//
// I. Hrivnacova

#if !defined(__CINT__) || defined(__MAKECINT__)

#include <TROOT.h>
#include <TSystem.h>
#include <TString.h>

#endif

void loadmacros () 
{
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

  // Load macros
  //
  gROOT->LoadMacro("$ALICE_ROOT/MUON/runSimulation.C++");      
  gROOT->LoadMacro("$ALICE_ROOT/MUON/runReconstruction.C++");      
  gROOT->LoadMacro("$ALICE_ROOT/MUON/fastMUONGen.C++");      
  gROOT->LoadMacro("$ALICE_ROOT/MUON/fastMUONSim.C++");      
  // gROOT->LoadMacro("$ALICE_ROOT/MUON/DecodeRecoCocktail.C++");      
  gROOT->LoadMacro("$ALICE_ROOT/MUON/ReadRecoCocktail.C++");     
  gROOT->LoadMacro("$ALICE_ROOT/MUON/MergeMuonLight.C++");       
  gROOT->LoadMacro("$ALICE_ROOT/MUON/MakeMUONFullMisAlignment.C++");
  gROOT->LoadMacro("$ALICE_ROOT/MUON/MakeMUONResMisAlignment.C++"); 
  gROOT->LoadMacro("$ALICE_ROOT/MUON/MakeMUONZeroMisAlignment.C++");
  gROOT->LoadMacro("$ALICE_ROOT/MUON/MUONAlignment.C++");           
  gROOT->LoadMacro("$ALICE_ROOT/MUON/MUONCheck.C++");               
  gROOT->LoadMacro("$ALICE_ROOT/MUON/MUONCheckDI.C++");          
  gROOT->LoadMacro("$ALICE_ROOT/MUON/MUONCheckMisAligner.C++");     
  gROOT->LoadMacro("$ALICE_ROOT/MUON/MUONefficiency.C++");          
  gROOT->LoadMacro("$ALICE_ROOT/MUON/MUONGenerateBusPatch.C++");    
  gROOT->LoadMacro("$ALICE_ROOT/MUON/MUONGenerateGeometryData.C++");
  gROOT->LoadMacro("$ALICE_ROOT/MUON/MUONGenerateTestGMS.C++");     
  gROOT->LoadMacro("$ALICE_ROOT/MUON/MUONmassPlot_ESD.C++");        
  gROOT->LoadMacro("$ALICE_ROOT/MUON/MUONplotefficiency.C++");      
  gROOT->LoadMacro("$ALICE_ROOT/MUON/MUONRawStreamTracker.C++");    
  gROOT->LoadMacro("$ALICE_ROOT/MUON/MUONRawStreamTrigger.C++");
  gROOT->LoadMacro("$ALICE_ROOT/MUON/MUONRecoCheck.C++");           
  // gROOT->LoadMacro("$ALICE_ROOT/MUON/MUONResoEffChamber.C++"); 
  gROOT->LoadMacro("$ALICE_ROOT/MUON/MUONStatusMap.C++");        
  gROOT->LoadMacro("$ALICE_ROOT/MUON/MUONTrigger.C++");        
  gROOT->LoadMacro("$ALICE_ROOT/MUON/MUONTriggerEfficiency.C++");   
  gROOT->LoadMacro("$ALICE_ROOT/MUON/MUONTriggerEfficiencyPt.C++"); 
  gROOT->LoadMacro("$ALICE_ROOT/MUON/MUONTriggerChamberEfficiency.C++"); 
  gROOT->LoadMacro("$ALICE_ROOT/MUON/TestMUONPreprocessor.C++");    
}
