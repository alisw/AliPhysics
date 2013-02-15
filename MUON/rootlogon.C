/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/// \ingroup macros
/// \file rootlogon.C
/// \brief Macro which is run when starting Root in MUON
///
/// It loads the MUON libraries needed for simulation and reconstruction
/// and sets the include path. 
///
/// \author Laurent Aphecetche

{
  cout << "Loading MUON libraries ..." << endl;
  gROOT->LoadMacro("${ALICE_ROOT}/MUON/loadlibs.C");
  gInterpreter->ProcessLine("loadlibs()");
    
  cout << "Setting include path ..." << endl;
  TString includePath = "-I${ALICE_ROOT}/STEER ";
  includePath        += "-I${ALICE_ROOT}/STEER/STEER ";
  includePath        += "-I${ALICE_ROOT}/STEER/STEERBase ";
  includePath        += "-I${ALICE_ROOT}/STEER/CDB ";
  includePath        += "-I${ALICE_ROOT}/STEER/ESD ";
  includePath        += "-I${ALICE_ROOT}/RAW ";
  includePath        += "-I${ALICE_ROOT}/FASTSIM ";
  includePath        += "-I${ALICE_ROOT}/EVGEN ";
  includePath        += "-I${ALICE_ROOT}/SHUTTLE/TestShuttle ";
  includePath        += "-I${ALICE_ROOT}/ITS ";
  includePath        += "-I${ALICE_ROOT}/MUON ";
  includePath        += "-I${ALICE_ROOT}/MUON/mapping ";
  includePath        += "-I${ALICE_ROOT}/RAW ";
  includePath        += "-I${ALICE_ROOT}/include ";
  
  if ( gSystem->Getenv("ALICE_INSTALL") ) includePath += "-I${ALICE_INSTALL}/include ";
  if ( gSystem->Getenv("ALICE_BUILD") ) includePath += "-I${ALICE_BUILD}/include ";

  gSystem->SetIncludePath(includePath.Data());
}
