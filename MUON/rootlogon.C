/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/// By Laurent Aphecetche

{
  cout << "Loading MUON libraries ..." << endl;
  gROOT->LoadMacro("${ALICE_ROOT}/MUON/loadlibs.C");
  gInterpreter->ProcessLine("loadlibs()");
    
  cout << "Setting include path ..." << endl;
  gSystem->SetIncludePath("-I${ALICE_ROOT}/include -I${ALICE_ROOT}/RAW -I${ALICE_ROOT}/MUON -I${ALICE_ROOT}/MUON/mapping");
}
