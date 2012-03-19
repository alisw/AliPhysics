enum {kNoOption, kOpenFlavor, kResonanceOnly};

//================================================================================================================================

void AliMuonForwardTrackAnalysis(const Char_t *readDir= ".",
				 Int_t option = kNoOption,
				 Int_t nMassBin = 100, 
				 Double_t massMin = 0.,
				 Double_t massMax = 10.,
				 const Char_t *outDir = ".",
				 Bool_t singleMuonAnalysis = kTRUE,
				 Bool_t muonPairAnalysis = kTRUE,
				 Int_t firstEvent = -1,
				 Int_t lastEvent = -1, 
				 Int_t myRandom = 0,
				 Int_t maxNWrongClusters = 999,
				 Double_t ptMinSingleMuons = 0.0) {
  
  gROOT -> LoadMacro("./AliMuonForwardTrackAnalysis.cxx+");
  //  AliLog::SetClassDebugLevel("AliMuonForwardTrackPair", 1);
  //  AliLog::SetClassDebugLevel("AliMuonForwardTrackAnalysis", 1);

  //  TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps", -1., -1., AliMagF::k5kG, AliMagF::kBeamTypeAA, 2750.));
  //  TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps", -1., -1., AliMagF::k5kG, AliMagF::kBeamTypepp, 7000.));
  TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps", -1., -1., AliMagF::k5kG));
  
  AliMuonForwardTrackAnalysis *myAnalysis = new AliMuonForwardTrackAnalysis();
  myAnalysis->ReadEvents(firstEvent, lastEvent);
  myAnalysis->SetInputDir(readDir);
  myAnalysis->SetOutputDir(outDir);
  myAnalysis->SetMassRange(nMassBin, massMin, massMax);
  myAnalysis->SetPtDimuRange(10, 0., 5.);
  myAnalysis->SetSingleMuonAnalysis(singleMuonAnalysis);
  myAnalysis->SetMuonPairAnalysis(muonPairAnalysis);
  myAnalysis->SetOption(option);
  myAnalysis->SetMatchTrigger(kTRUE);
  myAnalysis->SetMaxNWrongClustersMC(maxNWrongClusters);
  myAnalysis->SetPtMinSingleMuons(ptMinSingleMuons);

  myAnalysis->UseCutOnOffsetChi2(kFALSE);        // cut on the single muons

  myAnalysis->UseBransonForCut(kFALSE);
  myAnalysis->UseBransonForKinematics(kFALSE);

  myAnalysis->Init("MuonGlobalTracks.root");

  while (myAnalysis->LoadNextEvent()) continue;

  myAnalysis->Terminate(Form("outFiles/outFile.%d.%d.%d.root", myAnalysis->GetFirstEvent(), myAnalysis->GetLastEvent(), myRandom));

}

//================================================================================================================================

