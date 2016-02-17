// 
//  A macro to test the Global + Cooked Matrix ITSU tracker 
//

#if !defined(__CINT__) || defined(__MAKECINT__)
   #include <TSystem.h>
   #include <TStopwatch.h>
   #include <TPluginManager.h>

   #include "AliLog.h"
   #include "AliReconstruction.h"
   #include "AliITSURecoParam.h"
#endif

extern TSystem *gSystem;
extern TPluginManager *gPluginMgr;

void recGloCooked() {
  AliLog::SetClassDebugLevel("AliITSUReconstructor",1);

  gSystem->Load("libITSUpgradeBase");
  gSystem->Load("libITSUpgradeSim");
  gSystem->Load("libITSUpgradeRec");

  // Set ITS upgrade reconstructor
  gPluginMgr->AddHandler("AliReconstructor", "ITS",
		 "AliITSUReconstructor","ITS", "AliITSUReconstructor()");
  
  AliReconstruction rec;
  rec.SetRunReconstruction(""); // run cluster finder
  rec.SetRunTracking("ITS TPC"); // Turn on with ITS when tracker is implemented
  
  rec.SetRunVertexFinder(kTRUE); // to be implemented - CreateVertexer
  rec.SetRunMultFinder(kFALSE);   // to be implemented - CreateMultFinder
  rec.SetRunPlaneEff(kFALSE);     // to be implemented - CreateTrackleter

  //  rec.SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  rec.SetSpecificStorage("GRP/GRP/Data",
			 Form("local://%s",gSystem->pwd()));
  rec.SetSpecificStorage("ITS/Align/Data",
			 Form("local://%s",gSystem->pwd()));
  rec.SetSpecificStorage("ITS/Calib/RecoParam",
			 Form("local://%s",gSystem->pwd()));
  

  rec.SetRunQA(":");
  rec.SetRunGlobalQA(0);
  //AliLog::Flush();

  AliITSURecoParam *par=AliITSURecoParam::GetHighFluxParam();
  par->SetTracker(1);     // 1 is the Cooked Matrix tracker  
  par->SetSAonly(kFALSE); // kFALSE is the TPC+ITS mode
  rec.SetRecoParam("ITS",par);

  TStopwatch timer;
  timer.Start();
  rec.Run();
  timer.Stop();
  timer.Print();
}

