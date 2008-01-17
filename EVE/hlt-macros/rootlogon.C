
{
    
  cout << "Setting include path ..." << endl;
  TString includePath = "-I${ALICE_ROOT}/include ";
  includePath        += "-I${ALICE_ROOT}/EVE ";
  includePath        += "-I${ALICE_ROOT}/EVE/Alieve ";
  includePath        += "-I${ALICE_ROOT}/EVE/TEveUtil ";
  includePath        += "-I${ALICE_ROOT}/HLT/BASE ";
  includePath        += "-I${ALICE_ROOT}/HLT/TPCLib ";
  includePath        += "-I${ALICE_ROOT}/HLT/BASE/HOMER ";
  includePath        += "-I${ALICE_ROOT}/TPC ";
  includePath        += "-I${ALICE_ROOT}/RAW ";
  gSystem->SetIncludePath(includePath.Data());
}
