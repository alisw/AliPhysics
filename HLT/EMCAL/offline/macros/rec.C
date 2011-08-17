#include "../AliRoot-TRUNK-20110708-build/include/AliReconstruction.h"
void rec() {
  AliReconstruction reco;

  //reco.SetWriteESDfriend();
  //reco.SetWriteAlignmentData();

  reco.SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  reco.SetSpecificStorage("GRP/GRP/Data",
			  Form("local://%s",gSystem->pwd()));

  reco.SetRunQA(":") ;
  
  reco.SetQARefDefaultStorage("local://$ALICE_ROOT/QAref") ;
  
  //for (Int_t det = 0 ; det < AliQA::kNDET ; det++) {
//    reco.SetQACycles(det, 999) ;
    //reco.SetQAWriteExpert(det) ; 
  //}

  TStopwatch timer;
  timer.Start();
  reco.Run();
  timer.Stop();
  timer.Print();
}
