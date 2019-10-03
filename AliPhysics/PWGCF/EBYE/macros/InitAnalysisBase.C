//__________________________________________________//
AliPidBFBase *GetBaseBF(const char* analysisLevel = "AOD", 
				     const char* centralityName = 0x0,
				     Double_t centrMin = 0.,
				     Double_t centrMax = 100.,
                                     Bool_t bResonancesCut = kFALSE,
				     Bool_t bHBTcut = kFALSE,
                                     Double_t HBTCutValue = 0.02,
                                     Bool_t bConversionCut = kFALSE,
				     Double_t invMassForConversionCut = 0.04,
                                     Bool_t bMomentumDifferenceCut = kFALSE,
                                     Double_t fQCutMin = 0.0,
                                     Bool_t bVertexBinning = kFALSE,
                                     TString fArgEventClass = "EventPlane",
                                     Bool_t bMomentumOrdering = kTRUE,
                                     Int_t SingleVaribale=3,
                                     Int_t PairVariable=6) {
  //Function to setup the AliPidBFBase object and return it
  AliPidBFBase *gBalance = new AliPidBFBase();
  gBalance->SetAnalysisLevel(analysisLevel);
  gBalance->SetTrackVariableSingle(SingleVaribale);
  gBalance->SetTrackVariablePair(PairVariable);
  if(bHBTcut) gBalance->UseHBTCut(HBTCutValue);
  gBalance->UseMomentumOrdering(bMomentumOrdering);
  if(bResonancesCut) gBalance->UseResonancesCut();
  if(bConversionCut) gBalance->UseConversionCut(invMassForConversionCut);
  if(bMomentumDifferenceCut) gBalance->UseMomentumDifferenceCut(fQCutMin);
  if(bVertexBinning) gBalance->SetVertexZBinning();
  gBalance->SetEventClass(fArgEventClass);
  if(centralityName) gBalance->SetCentralityIdentifier(centralityName);
  gBalance->SetCentralityInterval(centrMin,centrMax);

  //Init the histograms
  gBalance->InitHistograms();
  
  return gBalance;
}

