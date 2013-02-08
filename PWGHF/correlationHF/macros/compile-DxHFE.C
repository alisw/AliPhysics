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
  gInterpreter->ExecuteMacro("$ALICE_ROOT/PWGHF/vertexingHF/macros/LoadLibraries.C");

  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libPWGHFhfe.so");
  gSystem->Load("libCORRFW");
  gSystem->AddIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_ROOT/PWGHF/vertexingHF -I$ALICE_ROOT/PWGHF/base -I$ALICE_ROOT/PWGHF/hfe ");

  TString dir("$ALICE_ROOT/PWGHF/correlationHF/");
  gROOT->LoadMacro(dir+"AliHFAssociatedTrackCuts.cxx+");
  gROOT->LoadMacro(dir+"AliReducedParticle.cxx+");
  gROOT->LoadMacro(dir+"AliHFCorrelator.cxx+");
  gROOT->LoadMacro(dir+"AliDxHFEParticleSelection.cxx+");
  gROOT->LoadMacro(dir+"AliDxHFEParticleSelectionD0.cxx+");
  gROOT->LoadMacro(dir+"AliDxHFEParticleSelectionEl.cxx+");
  gROOT->LoadMacro(dir+"AliDxHFEToolsMC.cxx+");
  gROOT->LoadMacro(dir+"AliDxHFEParticleSelectionMCD0.cxx+");
  gROOT->LoadMacro(dir+"AliDxHFEParticleSelectionMCEl.cxx+");
  gROOT->LoadMacro(dir+"AliDxHFECorrelation.cxx+");
  gROOT->LoadMacro(dir+"AliDxHFECorrelationMC.cxx+");
  gROOT->LoadMacro(dir+"AliAnalysisTaskDxHFEParticleSelection.cxx+");
  gROOT->LoadMacro(dir+"AliAnalysisTaskDxHFECorrelation.cxx+");

}
#elif
{
  cerr << "this macro can not be compiled, remove option '+'" << endl;
}
#endif
