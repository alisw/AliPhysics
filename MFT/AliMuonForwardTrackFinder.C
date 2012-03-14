//================================================================================================================================

void AliMuonForwardTrackFinder(Int_t run=0,
			       Bool_t bransonCorrection = kTRUE,
			       Double_t zVertexError=0.010,
			       Int_t matching=0,
			       const Char_t *readDir= ".",
			       const Char_t *outDir = ".",
			       Int_t nEventsToAnalyze = -1) {
  
  //  TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps", -1., -1, AliMagF::k5kG));

  // AliLog::SetClassDebugLevel("AliMuonForwardTrackFinder", 1);

  AliCDBManager* man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  man->SetSpecificStorage("GRP/GRP/Data", Form("local://%s",gSystem->pwd()));

  AliMuonForwardTrackFinder *finder = new AliMuonForwardTrackFinder();
  finder->SetDraw(kFALSE);
  finder->Init(run, readDir, outDir, nEventsToAnalyze);

  finder -> SetSigmaSpectrometerCut(4.0);
  finder -> SetSigmaClusterCut(4.0);
  finder -> SetChi2GlobalCut(2.0);
  finder -> SetRAbsorberCut(0.0);
  //  finder -> SetRAbsorberCut(26.4);
  finder -> SetLowPtCut(0.0);
  //  finder -> SetLowPtCut(0.5);
  finder -> SetVertexError(0.015, 0.015, zVertexError);
  finder -> SetMatchingMode(matching);                // 0 -> real matching   1 -> ideal matching
  //  finder -> SetMinResearchRadiusAtPlane(4, 0.0);
  //  finder -> SetMinResearchRadiusAtPlane(3, 0.0);
  finder -> SetMinResearchRadiusAtPlane(4, 0.50);
  finder -> SetMinResearchRadiusAtPlane(3, 0.05);
  finder -> SetBransonCorrection(bransonCorrection);

  while (finder->LoadNextTrack()) continue;

  if (finder->GetNRealTracksAnalyzed()) finder->Terminate();

}

//================================================================================================================================

