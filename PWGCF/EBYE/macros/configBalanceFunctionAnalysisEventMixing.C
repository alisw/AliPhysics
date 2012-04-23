//__________________________________________________//
AliBalanceEventMixing *GetBalanceFunctionObject(const char* analysisLevel = "ESD", 
				     const char* centralityName = 0x0,
				     Double_t centrMin = 0.,
				     Double_t centrMax = 100.,
				     Bool_t bShuffle = kFALSE) {
  //Function to setup the AliBalanceEventMixing object and return it
  AliBalanceEventMixing *gBalance = new AliBalanceEventMixing();
  gBalance->SetAnalysisLevel(analysisLevel);
  gBalance->SetShuffle(bShuffle);
  if(centralityName) gBalance->SetCentralityIdentifier(centralityName);
  gBalance->SetCentralityInterval(centrMin,centrMax);

  //Set all analyses separately
  //Rapidity
  gBalance->SetInterval(AliBalanceEventMixing::kRapidity,-0.8,0.8,64,0.0,1.6);  
  //Eta
  gBalance->SetInterval(AliBalanceEventMixing::kEta,-0.8,0.8,64,0.0,1.6);
  //Qlong
  gBalance->SetInterval(AliBalanceEventMixing::kQlong,-1,1,100,0.0,2.0);
  //Qout
  gBalance->SetInterval(AliBalanceEventMixing::kQout,-1,1,100,0.0,2.0);
  //Qside
  gBalance->SetInterval(AliBalanceEventMixing::kQside,-1,1,100,0.0,2.0);
  //Qinv
  gBalance->SetInterval(AliBalanceEventMixing::kQinv,-1,1,100,0.0,2.0);
  //Phi
  gBalance->SetInterval(AliBalanceEventMixing::kPhi,0.,360.,90,0.,180.0);

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

