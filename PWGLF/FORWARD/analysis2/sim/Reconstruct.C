Int_t getIntEnv(const char* name)
{
  TString env = gSystem->Getenv(name);
  if (env.IsNull()) return 0;
  return env.Atoi();
}
void SetSpecStore(AliCDBManager& s, 
		  const char* key, 
		  const char* sub)
{
  s.SetSpecificStorage(key, Form("alien://Folder=/alice/simulation/%s",sub));
}


void Reconstruct(UInt_t run) 
{
  // -----------------------------------------------------------------
  // 
  // Get GRP parameters.  Defines global "grp" as a pointer to GRPData
  //
  gROOT->Macro(Form("GRP.C(%d)", run));
  
  // -----------------------------------------------------------------
  // 
  // Basic setup 
  //
  AliReconstruction reco;
  reco.SetRunReconstruction("ITS TPC TRD TOF PHOS HMPID "
			    "EMCAL MUON FMD ZDC PMD T0 VZERO");


  // -----------------------------------------------------------------
  //
  // switch off cleanESD, write ESDfriends and Alignment data
  // 
  reco.SetCleanESD(kFALSE);
  reco.SetWriteESDfriend();
  reco.SetFractionFriends(.1);
  reco.SetWriteAlignmentData();

  // -----------------------------------------------------------------
  //
  // ITS Efficiency and tracking errors
  // 
  reco.SetRunPlaneEff(kTRUE);
  reco.SetUseTrackingErrorsForAlignment("ITS");

  
  // -----------------------------------------------------------------
  //
  // Raw OCDB
  //
  AliCDBManager* man = AliCDBManager::Instance();
  man->SetDefaultStorageFromRun(grp->run);

  // --- Get GRP to deduce collision system --------------------------
  Bool_t         isAA  = grp->IsAA();
  Bool_t         is10h = grp->period.EqualTo("LHC10h");

  // --- ITS (2 objects) ---------------------------------------------
  SetSpecStore(*man,"ITS/Align/Data",		"2008/v4-15-Release/Residual");
  SetSpecStore(*man,"ITS/Calib/SPDSparseDead",	"2008/v4-15-Release/Residual");


  // --- MUON objects (1 object) -------------------------------------
  SetSpecStore(*man,"MUON/Align/Data",		"2008/v4-15-Release/Residual");

  // --- TPC (7 objects) ---------------------------------------------
  SetSpecStore(*man,"TPC/Align/Data",		"2008/v4-15-Release/Residual");
  SetSpecStore(*man,"TPC/Calib/ClusterParam",	"2008/v4-15-Release/Residual");
  SetSpecStore(*man,"TPC/Calib/RecoParam",	"2008/v4-15-Release/Residual");
  SetSpecStore(*man,"TPC/Calib/TimeGain",	"2008/v4-15-Release/Residual");
  SetSpecStore(*man,"TPC/Calib/AltroConfig",	"2008/v4-15-Release/Residual");
  SetSpecStore(*man,"TPC/Calib/TimeDrift",	"2008/v4-15-Release/Residual");
  SetSpecStore(*man,"TPC/Calib/Correction",	"2008/v4-15-Release/Residual");


  // --- ZDC ---------------------------------------------------------
  // ZDC for 2010 the following is needed 
  // (https://savannah.cern.ch/task/?func=detailitem&item_id=33180#comment46)
  if (is10h) {
    reco.SetRecoParam("ZDC",AliZDCRecoParamPbPb::GetHighFluxParam(2760));
    SetSpecStore(*man,"ZDC/Align/Data",	"2008/v4-15-Release/Ideal/"); 
  }

  // --- GRP from local OCDB -----------------------------------------
  // man->SetSpecificStorage("GRP/GRP/Data",
  //                         Form("local://%s",gSystem->pwd()));
  

  // --- Override some settings in the ITS reco ----------------------
  if (is10h) {
    printf("Overriding ITS/Calib/RecoParam for run %d to do "
	   "reco even in absence of trigger\n", grp->run);
    man->SetRun(grp->run);
    AliCDBEntry*     entry  = man->Get("ITS/Calib/RecoParam");
    TObjArray*       array  = static_cast<TObjArray*>(entry->GetObject());
    AliITSRecoParam* par    = static_cast<AliITSRecoParam*>(array->RemoveAt(1));
    par->SetSkipSubdetsNotInTriggerCluster(kFALSE);
    reco.SetRecoParam("ITS",par);
  }

  // -----------------------------------------------------------------
  // 
  // Do not run QA in PbPb 
  if (isAA) reco.SetRunQA(":");

  // -------------------------------------------------------
  // 
  // Now run 
  // 
  TStopwatch timer;
  timer.Start();
  reco.Run();
  timer.Stop();
  timer.Print();
}
// 
// EOF
// 
