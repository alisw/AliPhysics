// Macro for muon reconstruction 
// 
// This macro is a simplified way to call the macro MUONrecoNtuple.C which produces
// an Ntuple with reconstructed tracks in output file "MUONtrackReco.root".
//// 
// Arguments:
//   FirstEvent (default 0)
//   LastEvent (default 0)
//   RecGeantHits (1 to reconstruct GEANT hits) (default 0)
//   FileName (for signal) (default "galice.root")
//   BkgGeantFileName (for background),
//   needed only if RecGeantHits = 1 and background to be added

void MUONrecoNtupleRun (Int_t FirstEvent = 0, Int_t LastEvent = 0, Int_t RecGeantHits = 0, Text_t *FileName = "galice.root", Text_t *BkgGeantFileName = "")
{

  TStopwatch timer;
  timer.Start();

  gSystem->SetIncludePath("-I$ALICE_ROOT/MUON -I$ALICE_ROOT/STEER -I$ROOTSYS/include") ;
  gROOT->ProcessLine(".x loadlibs.C");
  // force macros compilation (++), load the resulting shared library and execute the macro MUONrecoNtuple.C
  gROOT->ProcessLine(".x MUONrecoNtuple.C++(FirstEvent,LastEvent,RecGeantHits,FileName,BkgGeantFileName) ");

  timer.Stop();
  cout << "**** MUONrecoNtupleRun.C , timer" <<endl;
  timer.Print();

}





