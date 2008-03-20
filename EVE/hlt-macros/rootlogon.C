// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

{
  cout << "Setting include path ..." << endl;
  TString includePath = "-I${ALICE_ROOT}/include ";
  includePath        += "-I${ALICE_ROOT}/EVE ";
  includePath        += "-I${ALICE_ROOT}/EVE/EveBase ";
  includePath        += "-I${ALICE_ROOT}/EVE/EveDet ";
  includePath        += "-I${ALICE_ROOT}/EVE/EveHLT ";
  includePath        += "-I${ALICE_ROOT}/HLT/BASE ";
  includePath        += "-I${ALICE_ROOT}/HLT/TPCLib ";
  includePath        += "-I${ALICE_ROOT}/HLT/BASE/HOMER ";
  includePath        += "-I${ALICE_ROOT}/TPC ";
  includePath        += "-I${ALICE_ROOT}/RAW ";
  includePath        += "-I${ALICE_ROOT}/MUON ";
  includePath        += "-I${ALICE_ROOT}/MUON/mapping ";
  gSystem->SetIncludePath(includePath.Data());
}
