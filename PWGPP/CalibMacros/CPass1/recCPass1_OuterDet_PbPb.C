/*
   rec.C to be used for pass0
   - reconstruction of raw data
   - QA information switched off
   - store all friends
   - default OCDB storage set to "raw://"

   Example:
   aliroot -b -q 'recCPass1.C("raw.root",100)'
*/

AliITSRecoParam* GetSpecialITSRecoParam()
{
  AliITSRecoParam *itsRecoParam = AliITSRecoParam::GetHighFluxParam();
  itsRecoParam->SetTrackleterPhiWindowL2(0.07);
  itsRecoParam->SetTrackleterZetaWindowL2(0.4);
  itsRecoParam->SetTrackleterPhiWindowL1(0.10);
  itsRecoParam->SetTrackleterZetaWindowL1(0.6);
  itsRecoParam->SetTrackleterPhiWindow(0.06);
  itsRecoParam->SetTrackleterThetaWindow(0.025);
  itsRecoParam->SetTrackleterScaleDThetaBySin2T(kTRUE);
  itsRecoParam->SetTrackleterRemoveClustersFromOverlaps(kTRUE);
  
  itsRecoParam->SetVertexerZ();  
  itsRecoParam->ReconstructOnlySPD();

  return itsRecoParam;
}


void recCPass1_OuterDet(const char *filename="raw.root",Int_t nevents=-1, const char *ocdb="raw://")
{
  // Load some system libs for Grid and monitoring
  // Set the CDB storage location
  AliCDBManager * man = AliCDBManager::Instance();
  man->SetDefaultStorage(ocdb);
  // Reconstruction settings
  AliReconstruction rec;
  // Upload CDB entries from the snapshot (local root file) if snapshot exist
  if (gSystem->AccessPathName("OCDB.root", kFileExists)==0) {        
    //rec.SetFromCDBSnapshot("OCDB.root");
    rec.SetCDBSnapshotMode("OCDB.root");
  }
  if (gSystem->AccessPathName("localOCDBaccessConfig.C", kFileExists)==0) {        
    gROOT->LoadMacro("localOCDBaccessConfig.C");
    localOCDBaccessConfig();
  }

  // only SPD-trackletting will be done
  rec.SetRecoParam("ITS",GetSpecialITSRecoParam());

  // All friends
  rec.SetFractionFriends(1.0);

 // AliReconstruction settings - hardwired MB trigger for calibration

  TString newfilename = filename;
  newfilename += "?Trigger=kCalibOuter";
  rec.SetInput(newfilename.Data());

  // Set protection against too many events in a chunk (should not happen)
  if (nevents>0) rec.SetEventRange(0,nevents);

  // Remove recpoints after each event
  rec.SetDeleteRecPoints("ITS MUON EMCAL PHOS VZERO T0");

  // Switch off the V0 finder - saves time!
  rec.SetRunV0Finder(kFALSE); 

  //
  // QA options - all QA is off
  //
  rec.SetRunQA(":");
  rec.SetRunGlobalQA(kFALSE);

  // AliReconstruction settings
  rec.SetWriteESDfriend(kFALSE);
  rec.SetWriteAlignmentData(kFALSE);
  rec.SetUseTrackingErrorsForAlignment("ITS");
  rec.SetRunReconstruction("ITS MUON EMCAL PHOS VZERO T0");
  rec.SetFillESD("ITS MUON EMCAL PHOS VZERO T0");
  rec.SetCleanESD(kFALSE);

  // Specific reco params for ZDC (why isn't this automatic?)
//  rec.SetRecoParam("ZDC",AliZDCRecoParamPbPb::GetHighFluxParam(2760));

  //Ignore SetStopOnError
  rec.SetStopOnError(kFALSE);

  AliLog::Flush();
  rec.Run();
}

