void MakeGRPRecoParam(AliRecoParam::EventSpecie_t default=AliRecoParam::kLowMult) {
//========================================================================
//
// Steering macro for GRP reconstruction parameters
//
//
//========================================================================


  const char* macroname = "MakeGRPRecoParam.C";

  // Activate CDB storage and load geometry from CDB
  AliCDBManager* cdb = AliCDBManager::Instance();
  if(!cdb->IsDefaultStorageSet()) cdb->SetDefaultStorage("local://OCDB");
  
  TObjArray *recoParamArray = new TObjArray();

  {
    AliGRPRecoParam * param = AliGRPRecoParam::GetCosmicTestParam();
    param->SetEventSpecie(AliRecoParam::kCosmic);
    param->SetVertexerTracksConstraintITS(kFALSE);
    param->SetVertexerTracksConstraintTPC(kFALSE);
    recoParamArray->AddLast(param);
  }
  {
    // new settings for pass 2reco of Dec09 pp data
    AliGRPRecoParam * param = AliGRPRecoParam::GetLowFluxParam();
    param->SetEventSpecie(AliRecoParam::kLowMult);
    param->SetVertexerTracksConstraintITS(kTRUE);
    Double_t cutsITS[12]={0.1,
                          0.1,
			  0.5,
			  4, // minimum 4 clusters (default was 5)
                          1,
                          3.,
                          100.,
                          1000.,
                          3.,
                          30.,
                          1,
                          4};
    param->SetVertexerTracksCutsITS(12,cutsITS);
    param->SetVertexerTracksConstraintTPC(kTRUE);
    recoParamArray->AddLast(param);
  }
  {
    AliGRPRecoParam * param = AliGRPRecoParam::GetHighFluxParam();
    param->SetEventSpecie(AliRecoParam::kHighMult);
    param->SetVertexerTracksConstraintITS(kFALSE);
    param->SetVertexerTracksConstraintTPC(kFALSE);
    recoParamArray->AddLast(param);
  }

  // Set the default
  Bool_t defaultIsSet = kFALSE;
  for(Int_t i =0; i < recoParamArray->GetEntriesFast(); i++) {
    AliDetectorRecoParam *par = (AliDetectorRecoParam *)recoParamArray->UncheckedAt(i);
    if (!par) continue;
    if (default & par->GetEventSpecie()) {
      par->SetAsDefault();
      defaultIsSet = kTRUE;
    }
  }

  if (!defaultIsSet) {
    Error(macroname,"The default reconstruction parameters are not set! Exiting...");
    return;
  }

  // save in CDB storage
  AliCDBMetaData *md= new AliCDBMetaData();
  md->SetResponsible("Cvetan Cheshkov");
  md->SetComment("GRP reconstruction parameters");
  md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  md->SetBeamPeriod(0);
  AliCDBId id("GRP/Calib/RecoParam",0,AliCDBRunRange::Infinity());
  cdb->GetDefaultStorage()->Put(recoParamArray,id, md);

  return;
}
