////////////////////////////////////////////////////////////////////////
//
// name: TestMergeC.C
// date: 30.04.2002
// last update: 30.04.2002
// author: Jiri Chudoba
// version: 1.0
// description: 
//          test macro for AliMergeSteer class
//
////////////////////////////////////////////////////////////////////////

TestMergeC(Int_t nEvents = 2, Bool_t flagSim = kFALSE, 
	   Bool_t flagSDigits = kFALSE,
	   Bool_t flagMerge = kFALSE, 
	   Bool_t flagRecoMerged = kFALSE,
	   Bool_t flagCmpMerged = kFALSE) {

  gInterpreter->ProcessLine(".L $(ALICE_ROOT)/STEER/AliMergeSteer.C++");
//  gSystem->Load("$(ALICE_ROOT)/STEER/AliMergeSteer_C.so");

  AliMergeSteer* ms = new AliMergeSteer();

  ms->SetNEvents(nEvents);

  ms->SetFlagSim(flagSim);

  ms->SetFlagSDigits(flagSDigits);

  ms->SetFlagMerge(flagMerge);

  ms->SetFlagRecoMerged(flagRecoMerged);

  ms->SetFlagCmpMerged(flagCmpMerged);

//  ms->SetFlagCmpSignalOnly(kTRUE);

  ms->SetFileNameBgrHits("rfio:bgr.hits.root");
  ms->SetFileNameSDigits("tpc.sdigits.roi.root");
  ms->SetFileNameDigitsMerged("digits.root");

  ms->SetDetectorFlag("TPC",2);

  ms->Exec("");

  delete ms;
}
