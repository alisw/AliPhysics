//================================================================================================================================

void AliMuonForwardTrackFinder(Int_t run=0,
			       const Char_t *readDir= ".",
			       const Char_t *outDir = ".",
			       Int_t nEventsToAnalyze = -1) {
  
  //  TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps", -1., -1, AliMagF::k5kG));

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
  finder -> SetExtrapOriginTransvError(0.05);   // To be imposed if the gen. vertex is fixed in (0,0,0)
  finder -> SetGaussianBlurZVert(5.0);          // To be imposed if the gen. vertex is fixed in (0,0,0) 
  finder -> SetMatchingMode(0);                 // 0 -> real matching   1 -> ideal matching
  finder -> SetMinResearchRadiusAtLastPlane(0.5);

  while (finder->LoadNextTrack()) continue;

  if (finder->GetNRealTracksAnalyzed()) finder->Terminate();

}

//================================================================================================================================

