void MakeITSRecoParam(AliRecoParam::EventSpecie_t default=AliRecoParam::kLowMult, const char* cdbURI="local://") {
//========================================================================
//
// Steering macro for ITSU reconstruction parameters
//
// Contact: ruben.shahoyan@cern.ch
//
//========================================================================
  const char* macroname = "MakeITSRecoParam.C";
  //
  enum {kBit0=0x1<<0,kBit1=0x1<<1,kBit2=0x1<<2,kBit3=0x1<<3,kBit4=0x1<<4,kBit5=0x1<<5,kBit6=0x1<<6,kBit7=0x7<<2,kBit8=0x1<<8};
  //
  //
  gSystem->Load("libITSUpgradeBase.so");
  gSystem->Load("libITSUpgradeSim.so");
  gSystem->Load("libITSUpgradeRec.so");
  //
  // Activate CDB storage and load geometry from CDB
  AliCDBManager* cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage(cdbURI);
  AliITSUTrackCond* trcond = 0;
  int nLr = 7;
  
  TObjArray *recoParamArray = new TObjArray();
  //
  {
    AliITSURecoParam * itsRecoParam = AliITSURecoParam::GetCosmicTestParam();
    //
    itsRecoParam->SetNLayers(nLr);
    //
    //******************************************************************
    itsRecoParam->SetEventSpecie(AliRecoParam::kCosmic);
    itsRecoParam->SetTitle("Cosmic");
    recoParamArray->AddLast(itsRecoParam);
  }
  //
  int nBranch[7] = {10,10,10,10,10,10,10}; // max branching for the seed on layer
  int nCands[7]  = {10,500,500,500,300,200,100}; // max candidates for the TPC seed
  float tr2clChi2[7] = {50,50,50,50,50,80,100}; // cut on cluster to track chi2 
  float missPen[7] = {2.,2.,2.,2.,2.,2.,2.};    // missing cluster penalty
  //
  {
    AliITSURecoParam * itsRecoParam = AliITSURecoParam::GetLowFluxParam();
    //
    itsRecoParam->SetNLayers(nLr);
    //
    //******************************************************************
    itsRecoParam->SetEventSpecie(AliRecoParam::kLowMult);
    itsRecoParam->SetTitle("LowMult");
    recoParamArray->AddLast(itsRecoParam);
    //******************************************************************
    for (int i=0;i<nLr;i++) itsRecoParam->SetAllowDiagonalClusterization(i,kTRUE);
    //  
    // Add tracking conditions >>>
    trCond = new AliITSUTrackCond();
    trCond->SetNLayers(nLr); 
    //
    for (int i=0;i<nLr;i++) {
      trCond->SetMaxBranches(i,nBranch[i]);    // each seed propagated to given layer can produce max nBranch branches
      trCond->SetMaxCandidates(i,nCands[i]);   // each tpc track may have at most nCands prolongations
      trCond->SetMaxTr2ClChi2(i,tr2clChi2[i]); // cut on cluster to track chi2
      trCond->SetMissPenalty(i,missPen[i]);    // missing cluster penalty
      //
    }
    //
    trCond->AddNewCondition(5); // min hits
    trCond->AddGroupPattern( kBit0|kBit1 );
    trCond->AddGroupPattern( kBit3|kBit4 );
    trCond->AddGroupPattern( kBit5|kBit6 );
    //
    trCond->AddNewCondition(5); // min hits
    trCond->AddGroupPattern( kBit0|kBit1 );
    trCond->AddGroupPattern( kBit3|kBit4 );
    trCond->AddGroupPattern( kBit5|kBit6 );
    //
    trCond->AddNewCondition(5);
    trCond->AddGroupPattern( kBit0|kBit2 );
    trCond->AddGroupPattern( kBit3|kBit4 );
    trCond->AddGroupPattern( kBit5|kBit6 );
    //
    trCond->AddNewCondition(5);
    trCond->AddGroupPattern( kBit1|kBit2 );
    trCond->AddGroupPattern( kBit3|kBit4 );
    trCond->AddGroupPattern( kBit5|kBit6 );
    //
    trCond->Init();
    //
    itsRecoParam->AddTrackingCondition(trCond);
    // Add tracking conditions <<<
  }
  {
    AliITSURecoParam * itsRecoParam = AliITSURecoParam::GetHighFluxParam();
    //
    itsRecoParam->SetNLayers(nLr);
    //
    //******************************************************************
    itsRecoParam->SetEventSpecie(AliRecoParam::kHighMult);
    itsRecoParam->SetTitle("HighMult");
    recoParamArray->AddLast(itsRecoParam);
    //******************************************************************
     for (int i=0;i<nLr;i++) itsRecoParam->SetAllowDiagonalClusterization(i,kTRUE);
    //  
    // Add tracking conditions >>>
    trCond = new AliITSUTrackCond();
    trCond->SetNLayers(nLr); 
    //
    for (int i=0;i<nLr;i++) {
      trCond->SetMaxBranches(i,nBranch[i]);    // each seed propagated to given layer can produce max nBranch branches
      trCond->SetMaxCandidates(i,nCands[i]);   // each tpc track may have at most nCands prolongations
      trCond->SetMaxTr2ClChi2(i,tr2clChi2[i]);   // cut on cluster to track chi2
      trCond->SetMissPenalty(i,missPen[i]);    // missing cluster penalty
      //
    }
    //
    trCond->AddNewCondition(5); // min hits
    trCond->AddGroupPattern( kBit0|kBit1 );
    trCond->AddGroupPattern( kBit3|kBit4 );
    trCond->AddGroupPattern( kBit5|kBit6 );
    //
    trCond->AddNewCondition(5);
    trCond->AddGroupPattern( kBit0|kBit2 );
    trCond->AddGroupPattern( kBit3|kBit4 );
    trCond->AddGroupPattern( kBit5|kBit6 );
    //
    trCond->AddNewCondition(5);
    trCond->AddGroupPattern( kBit1|kBit2 );
    trCond->AddGroupPattern( kBit3|kBit4 );
    trCond->AddGroupPattern( kBit5|kBit6 );
    //
    trCond->Init();
    //
    itsRecoParam->AddTrackingCondition(trCond);
    // Add tracking conditions <<<
    //
  }
  //
  // Set the default
  Bool_t defaultIsSet = kFALSE;
  for(Int_t i =0; i < recoParamArray->GetEntriesFast(); i++) {
    AliDetectorRecoParam *param = (AliDetectorRecoParam *)recoParamArray->UncheckedAt(i);
    if (!param) continue;
    if (default & param->GetEventSpecie()) {
      param->SetAsDefault();
      defaultIsSet = kTRUE;
    }
  }

  if (!defaultIsSet) {
    Error(macroname,"The default reconstruction parameters are not set! Exiting...");
    return;
  }

  // save in CDB storage
  AliCDBMetaData *md= new AliCDBMetaData();
  md->SetResponsible("Andrea Dainese");
  md->SetComment("Reconstruction parameters ITS.");
  md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  md->SetBeamPeriod(0);
  AliCDBId id("ITS/Calib/RecoParam",0,AliCDBRunRange::Infinity());
  cdb->GetDefaultStorage()->Put(recoParamArray,id, md);
  //
  return;
}

