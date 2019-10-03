void CleanOADB(const TString& dest="tmp.root",
	       Bool_t         verb=false,
	       Bool_t         all=false)
{
  const char* fwd = "${ALICE_PHYSICS}/PWGLF/FORWARD/analysis2";
  if (!gROOT->GetClass("AliOADBForward")) {
    Info("", "Loading libraries");
    gROOT->Macro(Form("%s/scripts/LoadLibs.C", fwd));
  }

  AliForwardCorrectionManager& fcm = 
    AliForwardCorrectionManager::Instance();
  fcm.EnableCorrections(0x7F);
  // if (!fcm.Init(0,1,900,+5,false,false,0x7F,false)) { 
  // Warning("CleanOADB", "Failed to init");
  // return;
  // }
  fcm.Print("R");


  Info("", "Cleaning up");
  fcm.CleanUp(dest, verb, all);
  Info("", "done");
}
