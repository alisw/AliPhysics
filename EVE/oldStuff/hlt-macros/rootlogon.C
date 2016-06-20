// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

{
  cout<<"Loading HLT libraries..."<<endl;
  gSystem->Load("libHLTbase");
  gSystem->Load("libAliHLTUtil");
  gSystem->Load("libHLTinterface");
  gSystem->Load("libAliHLTMUON");
  gSystem->Load("libAliHLTTPC");
  gSystem->Load("libAliHLTTRD");
  gSystem->Load("libTRDbase");
  gSystem->Load("libTRDrec");
 
  cout << "Setting include path ..." << endl;
  TString includePath = "-I${ALICE_ROOT}/include ";
  includePath        += "-I${ALICE_ROOT}/EVE ";
  includePath        += "-I${ALICE_ROOT}/EVE/macros ";
  includePath        += "-I${ALICE_ROOT}/EVE/EveBase ";
  includePath        += "-I${ALICE_ROOT}/EVE/EveDet ";
  includePath        += "-I${ALICE_ROOT}/EVE/EveHLT ";
  includePath        += "-I${ALICE_ROOT}/HLT/BASE ";
  includePath        += "-I${ALICE_ROOT}/HLT/TPCLib ";
  includePath        += "-I${ALICE_ROOT}/HLT/PHOS ";
  includePath        += "-I${ALICE_ROOT}/HLT/BASE/HOMER ";
  includePath        += "-I${ALICE_ROOT}/HLT/BASE/util ";
  includePath        += "-I${ALICE_ROOT}/HLT/TRD ";
  includePath        += "-I${ALICE_ROOT}/ITS ";
  includePath        += "-I${ALICE_ROOT}/PHOS ";
  includePath        += "-I${ALICE_ROOT}/TPC ";
  includePath        += "-I${ALICE_ROOT}/TRD ";
  includePath        += "-I${ALICE_ROOT}/RAW ";
  includePath        += "-I${ALICE_ROOT}/MUON ";
  includePath        += "-I${ALICE_ROOT}/MUON/mapping ";
  includePath        += "-I${ALICE_ROOT}/HLT/MUON ";
  includePath        += "-I${ALICE_ROOT}/HLT/TPCLib/tracking-ca ";
  includePath        += "-I${ALICE_ROOT}/PWG0";
  
  gSystem->SetIncludePath(includePath.Data());
}
