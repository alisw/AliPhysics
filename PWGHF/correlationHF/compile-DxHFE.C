//-*- Mode: C++ -*-
// $Id$

//* This file is property of and copyright by the ALICE Project        * 
//* ALICE Experiment at CERN, All rights reserved.                     *
//* See cxx source for full Copyright notice                           *

/// @file   compile-DxHFE.C
/// @author Matthias Richter
/// @date   2012-03-19
/// @brief  Compile classes for the D0 - HFE correlation studies
///
/// Usage: aliroot -l $ALICE_ROOT/PWGHF/correlationHF/compile-DxHFE.C

#if defined(__CINT__) && !defined(__MAKECINT__)
{
  //----------- Loading the required libraries ---------//

  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->AddIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/include");
  gROOT->LoadMacro("$ALICE_ROOT/PWGHF/correlationHF/AliDxHFEParticleSelection.cxx+");
  gROOT->LoadMacro("$ALICE_ROOT/PWGHF/correlationHF/AliDxHFEParticleSelectionD0.cxx+");
  gROOT->LoadMacro("$ALICE_ROOT/PWGHF/correlationHF/AliDxHFEParticleSelectionEl.cxx+");
  gROOT->LoadMacro("$ALICE_ROOT/PWGHF/correlationHF/AliAnalysisTaskDxHFEParticleSelection.cxx+");
  gROOT->LoadMacro("$ALICE_ROOT/PWGHF/correlationHF/AliDxHFECorrelation.cxx+");
  gROOT->LoadMacro("$ALICE_ROOT/PWGHF/correlationHF/AliAnalysisTaskDxHFECorrelation.cxx+");
}
#elif
{
  cerr << "this macro can not be compiled, remove option '+'" << endl;
}
#endif
