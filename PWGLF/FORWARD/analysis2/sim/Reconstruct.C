/**
 * @file   Reconstruct.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Wed Oct 15 13:27:45 2014
 * 
 * @brief  Steering script for the reconstruction 
 */
/** 
 * Run the reconstruction 
 * 
 * @param run Run number 
 */
void Reconstruct(UInt_t run) 
{
  // -----------------------------------------------------------------
  // 
  // - Get GRP parameters.  Defines global "grp" as a pointer to GRPData
  // - Load base class definitions in BaseConfig.C
  // - Get which detectors are turned on in "detCfg". 
  // - Create the OCDB configuration object "ocdbCfg"
  // 
  gROOT->Macro(Form("GRP.C(%d)", run));
  gROOT->Macro("BaseConfig.C");
  gROOT->Macro("DetConfig.C"); 
  gROOT->Macro("OCDBConfig.C"); 

  // --- Get GRP to deduce collision system --------------------------
  Bool_t         isAA  = grp->IsAA();
  Bool_t         is10h = grp->period.EqualTo("LHC10h");
 
  // -----------------------------------------------------------------
  // 
  // Basic setup 
  //
  AliReconstruction reco;
  TString enable;  
  detCfg->GetRecoString(enable);
  // if (is10h) enable.ReplaceAll("MUON", "");
  reco.SetRunReconstruction(enable);

  // -----------------------------------------------------------------
  //
  // switch off cleanESD, write ESDfriends and Alignment data, clean
  // up rec-points (except for ITS)
  // 
  reco.SetCleanESD(kFALSE);
  reco.SetWriteESDfriend();
  reco.SetFractionFriends(.1);
  reco.SetWriteAlignmentData();
  TString clean(enable); clean.ReplaceAll("ITS", "");
  reco.SetDeleteRecPoints(clean);

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
  ocdbCfg->Init(false);

  // -----------------------------------------------------------------
  // 
  // Specific reconstruction parameters 
  // 
  // --- ZDC ---------------------------------------------------------
  // ZDC for 2010 the following is needed 
  // (https://savannah.cern.ch/task/?func=detailitem&item_id=33180#comment46)
  if (is10h)
    reco.SetRecoParam("ZDC",AliZDCRecoParamPbPb::GetHighFluxParam(2760));

  // --- Override some settings in the ITS reco ----------------------
  // 
  // This was needed for the ITS stand-alone tracks to be produced in
  // LHC14b8a and LHC14b8b.  However, using the LowMult reco-param
  // object for all events means using AliITSVertexer3D irrespective
  // of multiplicity.  AliITSVertexer3D can use more than 4GB of
  // memory for central Hijing events, which means that the jobs are
  // likely to hit the 8GB hard limit.
  if (false && is10h) {
    printf("Overriding ITS/Calib/RecoParam for run %d to do "
	   "reco even in absence of trigger\n", grp->run);
    man->SetRun(grp->run);
    AliCDBEntry*     entry  = man->Get("ITS/Calib/RecoParam");
    TObjArray*       array  = static_cast<TObjArray*>(entry->GetObject());
    AliITSRecoParam* par    = static_cast<AliITSRecoParam*>(array->RemoveAt(1));
    par->SetSkipSubdetsNotInTriggerCluster(kFALSE);
    reco.SetRecoParam("ITS",par);
  }
  TString host(gSystem->HostName());
  Bool_t  cholm = host.BeginsWith("hehi");
  if (cholm) {
    Warning("Reconstruct", "Disabling check of TPC reco-params");
    reco.RemCheckRecoCDBvsSimuCDB("TPC/Calib/RecoParam");
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
