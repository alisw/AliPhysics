void MakeGRPRecoParam(Bool_t allowCleanup, AliRecoParam::EventSpecie_t default=AliRecoParam::kLowMult) {
//========================================================================
//
// Steering macro for GRP reconstruction parameters
//
//
//========================================================================


  const char* macroname = "MakeGRPRecoParam.C";

  printf("V0s and TPC-only track cleanup : %s\n",allowCleanup ? "ON":"OFF");
  
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
    Double_t cutsITS[24]={0.1,
                          0.1,
			  0.5,
			  //4, // minimum 4 clusters (default was 5)
			  3, // minimum 3 clusters (was 4)
                          1,
                          3.,
                          100.,
                          1000.,
                          3.,
                          30.,
                          6, // 6: MultiVertexer (was 1)
                          4,
			  // multivertexer settings
			  7., 
			  1e3.,
			  5.0,
			  0.05,
			  10e-4,
			  2.,
			  10.,
			  1.,
			  25.,
			  0.,
			  999999.,
			  3.
    };
    param->SetVertexerTracksCutsITS(24,cutsITS);
    param->SetVertexerTracksConstraintTPC(kTRUE);

    if (allowCleanup) {
      // V0 will be validated if it passes at least 1 mass hypothesis
      param->AddV0HypSel( AliV0HypSel("gamma",0.5486e-3, 0.5486e-3, 1.099e-3, 0.001, 20, 0.6, 0.,0.0));
      param->AddV0HypSel( AliV0HypSel("K0",139.570e-3, 139.570e-3, 497.7e-3, 0.003,20,0.07, 1.,0.5));
      param->AddV0HypSel( AliV0HypSel("Lambda",938.272e-3, 139.570e-3, 1115.683e-3, 0.001, 20, 0.07, 1.,0.5));
      param->AddV0HypSel( AliV0HypSel("antiLambda",139.570e-3, 938.272e-3, 1115.683e-3, 0.001, 20, 0.07, 1.,0.5));
      param->AddV0HypSel( AliV0HypSel("HyperTriton",2.8092, 139.570e-3, 2.992, 0.0025, 14, 0.07, 1.,0.5));
      param->AddV0HypSel( AliV0HypSel("antiHyperTriton",139.570e-3, 2.8092, 2.992, 0.0025, 14, 0.07, 1.,0.5));
      
      param->SetFlagsNotToClean(AliESDtrack::kITSin | AliESDtrack::kTRDrefit |
				AliESDtrack::kTOFout | AliESDtrack::kHMPIDout);
      param->SetVertexerV0EtaMax(1.);
      param->SetCleanOfflineV0Prongs(kTRUE); // fill 0s in redundant prongs
      param->SetCleanDCAZCut(30.);
    }
    
    recoParamArray->AddLast(param);
  }
  {
    AliGRPRecoParam * param = AliGRPRecoParam::GetHighFluxParam();
    param->SetEventSpecie(AliRecoParam::kHighMult);
    param->SetVertexerTracksConstraintITS(kTRUE);
    Double_t cutsITS[24]={0.1,
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
                          1,
			  // multivertexer settings
			  7., 
			  1e3.,
			  5.0,
			  0.05,
			  10e-4,
			  2.,
			  10.,
			  1.,
			  25.,
			  0.,
			  999999.,
			  3.
    }; // faster finder algo for Iteration 0
    param->SetVertexerTracksCutsITS(24,cutsITS);
    param->SetVertexerTracksConstraintTPC(kTRUE);

    if (allowCleanup) {
      // V0 will be validated if it passes at least 1 mass hypothesis
      param->AddV0HypSel( AliV0HypSel("gamma",0.5486e-3, 0.5486e-3, 1.099e-3, 0.001, 20, 0.6, 0.,0.0));
      param->AddV0HypSel( AliV0HypSel("K0",139.570e-3, 139.570e-3, 497.7e-3, 0.003,20,0.07, 1.,0.5));
      param->AddV0HypSel( AliV0HypSel("Lambda",938.272e-3, 139.570e-3, 1115.683e-3, 0.001, 20, 0.07, 1.,0.5));
      param->AddV0HypSel( AliV0HypSel("antiLambda",139.570e-3, 938.272e-3, 1115.683e-3, 0.001, 20, 0.07, 1.,0.5));
      param->AddV0HypSel( AliV0HypSel("HyperTriton",2.8092, 139.570e-3, 2.992, 0.0025, 14, 0.07, 1.,0.5));
      param->AddV0HypSel( AliV0HypSel("antiHyperTriton",139.570e-3, 2.8092, 2.992, 0.0025, 14, 0.07, 1.,0.5));
      param->SetVertexerV0EtaMax(1.);
      param->SetCleanOfflineV0Prongs(kTRUE); // fill 0s in redundant prongs
      
      param->SetFlagsNotToClean(AliESDtrack::kITSin | AliESDtrack::kTRDrefit |
				AliESDtrack::kTOFout | AliESDtrack::kHMPIDout);
      param->SetCleanDCAZCut(30.);
    }
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
