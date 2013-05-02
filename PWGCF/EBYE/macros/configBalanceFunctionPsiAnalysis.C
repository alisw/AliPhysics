//__________________________________________________//
AliBalancePsi *GetBalanceFunctionObject(const char* analysisLevel = "AOD", 
					const char* centralityName = 0x0,
					Double_t centrMin = 0.,
					Double_t centrMax = 100.,
					Bool_t bShuffle = kFALSE,
					Bool_t bResonancesCut = kFALSE,
					Bool_t bHBTCut = kFALSE,
					Bool_t bConversionCut = kFALSE,
					Bool_t bMomentumDifferenceCut = kFALSE,
					TString fArgEventClass = "EventPlane",
					Double_t deltaEtaMax = 2.0,
					Bool_t bVertexBinning = kFALSE) {
  //Function to setup the AliBalance object and return it
  AliBalancePsi *gBalance = new AliBalancePsi();
  gBalance->SetAnalysisLevel(analysisLevel);
  gBalance->SetShuffle(bShuffle);
  if(bResonancesCut) gBalance->UseResonancesCut();
  if(bHBTCut) gBalance->UseHBTCut();
  if(bConversionCut) gBalance->UseConversionCut();
  if(bMomentumDifferenceCut) gBalance->UseMomentumDifferenceCut();
  if(centralityName) gBalance->SetCentralityIdentifier(centralityName);
  if(bVertexBinning) gBalance->SetVertexZBinning();
  gBalance->SetCentralityInterval(centrMin,centrMax);
  gBalance->SetEventClass(fArgEventClass);
  gBalance->SetDeltaEtaMax(deltaEtaMax);

  //Set all analyses separately
  //Rapidity
  //gBalance->SetInterval(AliBalance::kRapidity,-0.8,0.8,32,-1.6,1.6,15.);  
  //Eta
  //gBalance->SetInterval(AliBalance::kEta,-0.8,0.8,32,-1.6,1.6,15);
  //Qlong
  //gBalance->SetInterval(AliBalance::kQlong,-1,1,200,0.0,4.0,15);
  //Qout
  //gBalance->SetInterval(AliBalance::kQout,-1,1,200,0.0,4.0,15);
  //Qside
  //gBalance->SetInterval(AliBalance::kQside,-1,1,200,0.0,4.0,15);
  //Qinv
  //gBalance->SetInterval(AliBalance::kQinv,-1,1,200,0.0,4.0,15);
  //Phi
  //gBalance->SetInterval(AliBalance::kPhi,0.,360.,90,-180.,180.0,15);

  //Init the histograms
  gBalance->InitHistograms();
  
  return gBalance;
}

//__________________________________________________//
AliESDtrackCuts *GetTrackCutsObject(Double_t ptMin, Double_t ptMax, Double_t etaMin, Double_t etaMax, Double_t maxTPCchi2, Double_t maxDCAz, Double_t maxDCAxy, Int_t minNClustersTPC) {

  // only used for ESDs
  // Function to setup the AliESDtrackCuts object and return it
  AliESDtrackCuts *cuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
  cuts->SetRequireTPCStandAlone(kTRUE); // TPC only cuts!  

  // extra TPC cuts (Syst studies)
  if(minNClustersTPC != -1)  cuts->SetMinNClustersTPC(minNClustersTPC);
  else cuts->SetMinNClustersTPC(70); // standard for filter bit 128
  
  if(maxTPCchi2 != -1) cuts->SetMaxChi2PerClusterTPC(maxTPCchi2);

  // extra DCA cuts (Syst studies)  
  if(maxDCAz!=-1 && maxDCAxy != -1){
    cuts->SetMaxDCAToVertexZ(maxDCAz);
    cuts->SetMaxDCAToVertexXY(maxDCAxy);
  }

  cuts->SetPtRange(ptMin,ptMax);
  cuts->SetEtaRange(etaMin,etaMax);
  cuts->DefineHistograms(1);
  //cuts->SaveHistograms("trackCuts");

  return cuts;
}

