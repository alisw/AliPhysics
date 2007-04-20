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
// DecodeRecoCocktail.C          - ok, comp,  README; Hermine, Alessandro
// ReadRecoCocktail.C            - ok, x,     README; Hermine, Alessandro
// MergeMuonLight.C              - x,  x,     README; Hermine, Alessandro
// MakeMUONFullMisAlignment.C    - ok, comp,  README; Javier, Ivana
// MakeMUONResMisAlignment.C     - ok, comp,  README; Javier, Ivana
// MakeMUONZeroMisAlignment.C    - ok, comp,  README; Javier, Ivana
// MUONAlignment.C               - ok, comp,  README; Javier
// MUONCheck.C                   - ok, comp,  x,      Frederic, in test
// MUONCheckDI.C                 - x,  !comp, x       Artur     
// MUONCheckMisAligner.C         - ok, comp,  x,      Javier
// MUONdisplay.C                 - ok, comp,  deprecated
// MUONefficiency.C              - ok, comp,  README; Christophe, in test
// MUONGenerateBusPatch.C        - ok, comp,  x,      Christian
// MUONGenerateGeometryData.C    - ok, comp,  README; Ivana
// MUONGenerateTestGMS.C         - ok, comp,  x,      Ivana
// MUONmassPlot_ESD.C            - ok, comp,  x,      Christian
// MUONplotefficiency.C          - ok, x,     README; Christophe
// MUONRawStreamTracker.C        - x,  comp,  README; Christian
// MUONRawStreamTrigger.C        - x,  comp,  README; Christian
// MUONRecoCheck.C               - ok, comp,  README; Hermine, Alessandro
// MUONResoEffChamber.C          - ok, comp,  x,      Nicolas
// MUONStatusMap.C               - x,  !comp, x,      Laurent
// MUONTracker.C                 - no, !comp, README; Philippe P.
// MUONTrigger.C                 - ok, comp,  README, Philippe C.
// MUONTriggerEfficiency.C       - ok, comp,  x,      Philippe C., in test
// MUONTriggerEfficiencyPt.C     - x,  comp,  README, Philippe C.
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
  includePath        += "-I${ALICE_ROOT}/SHUTTLE/TestShuttle ";
  includePath        += "-I${ALICE_ROOT}/ITS ";
  includePath        += "-I${ALICE_ROOT}/MUON ";
  includePath        += "-I${ALICE_ROOT}/MUON/mapping";
  gSystem->SetIncludePath(includePath.Data());

  // Load libraries not linked with aliroot
  //
  gSystem->Load("../SHUTTLE/TestShuttle/libTestShuttle.so");
  gSystem->Load("libMUONshuttle.so");
  gSystem->Load("libMUONevaluation.so");

  // Load macros
  //
  gROOT->LoadMacro("DecodeRecoCocktail.C++");       
  // gROOT->LoadMacro("ReadRecoCocktail.C++");         
  // gROOT->LoadMacro("MergeMuonLight.C++");           
  gROOT->LoadMacro("MakeMUONFullMisAlignment.C++"); 
  gROOT->LoadMacro("MakeMUONResMisAlignment.C++");  
  gROOT->LoadMacro("MakeMUONZeroMisAlignment.C++"); 
  gROOT->LoadMacro("MUONAlignment.C++");            
  gROOT->LoadMacro("MUONCheck.C++");                
  // gROOT->LoadMacro("MUONCheckDI.C++");              
  gROOT->LoadMacro("MUONCheckMisAligner.C++");      
  gROOT->LoadMacro("MUONdisplay.C++");              
  gROOT->LoadMacro("MUONefficiency.C++");           
  gROOT->LoadMacro("MUONGenerateBusPatch.C++");     
  gROOT->LoadMacro("MUONGenerateGeometryData.C++"); 
  gROOT->LoadMacro("MUONGenerateTestGMS.C++");      
  gROOT->LoadMacro("MUONmassPlot_ESD.C++");         
  // gROOT->LoadMacro("MUONplotefficiency.C++");       
  gROOT->LoadMacro("MUONRawStreamTracker.C++");     
  gROOT->LoadMacro("MUONRawStreamTrigger.C++");     
  gROOT->LoadMacro("MUONRecoCheck.C++");            
  gROOT->LoadMacro("MUONResoEffChamber.C++");       
  // gROOT->LoadMacro("MUONStatusMap.C++");            
  // gROOT->LoadMacro("MUONTracker.C++");              
  gROOT->LoadMacro("MUONTrigger.C++");              
  gROOT->LoadMacro("MUONTriggerEfficiency.C++");    
  gROOT->LoadMacro("MUONTriggerEfficiencyPt.C++");  
  gROOT->LoadMacro("TestMUONPreprocessor.C++");     
}
