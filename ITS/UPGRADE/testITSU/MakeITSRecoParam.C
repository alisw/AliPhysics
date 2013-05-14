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
  const Bool_t kAllowDiagCl = kFALSE;
  //
  // tuned for 20x20 pixels with x/x0= 0.3,0.3,0.3,0.5,0.5,0.5,0.5%
  int nBranch[7] = {5,10,15,4,6,6,10}; // max branching for the seed on layer
  int nCands[7]  = {10,20,45,20,45,15,10}; // max candidates for the TPC seed
  float tr2clChi2[7] = {20,25,30,40,45,45,70}; // cut on cluster to track chi2 
  float gloChi2[7]   = {10, 15,20,40,70,70,70}; // cut on seed global norm chi2
  float missPen[7] = {2.,2.,2.,2.,2.,2.,2.};    // missing cluster penalty
  float maxChi2Match = 10.;
  float maxChi2SA    = 10.;  
  //
  /*
    // tuned for 20x33 pixels
  int nBranch[7] = {5,10,15,4,6,6,10}; // max branching for the seed on layer
  int nCands[7]  = {10,20,25,20,30,15,20}; // max candidates for the TPC seed
  float tr2clChi2[7] = {20,20,25,25,25,30,40}; // cut on cluster to track chi2 
  float gloChi2[7]   = {9, 10,15,20,30,30,30}; // cut on seed global norm chi2
  float missPen[7] = {2.,2.,2.,2.,2.,2.,2.};    // missing cluster penalty
  float maxChi2Match = 10.;
  float maxChi2SA    = 10.;      
  */

  /*
  // this is for tuning only
  int nBranch[7] = {30,20,20,20,20,20,20}; // max branching for the seed on layer
  int nCands[7]  = {20,100,100,100,100,55,20}; // max candidates for the TPC seed
  float tr2clChi2[7] = {40,50,60,80,80,80,80}; // cut on cluster to track chi2 
  float gloChi2[7]   = {20, 30,40,60,80,80,80}; // cut on seed global norm chi2
  float missPen[7] = {2.,2.,2.,2.,2.,2.,2.};    // missing cluster penalty
  float maxChi2Match = 20.;
  float maxChi2SA    = 20.;  
  */
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
    for (int i=0;i<nLr;i++) itsRecoParam->SetAllowDiagonalClusterization(i,kAllowDiagCl);
    //  
    // Add tracking conditions >>>
    trCond = new AliITSUTrackCond();
    trCond->SetNLayers(nLr); 
    trCond->SetMaxITSTPCMatchChi2(maxChi2Match);
    trCond->SetMaxITSSAChi2(maxChi2SA);
    //
    // to exclude some layer use trCon->ExcludeLayer(lrID);
    //
    for (int i=0;i<nLr;i++) {
      trCond->SetMaxBranches(i,nBranch[i]);    // each seed propagated to given layer can produce max nBranch branches
      trCond->SetMaxCandidates(i,nCands[i]);   // each tpc track may have at most nCands prolongations
      trCond->SetMaxTr2ClChi2(i,tr2clChi2[i]); // cut on cluster to track chi2
      trCond->SetMaxChi2GloNrm(i,gloChi2[i]);  // cut on cluster to track global chi2
      trCond->SetMissPenalty(i,missPen[i]);    // missing cluster penalty
      //
    }
    //
    trCond->AddNewCondition(5); // min hits
    trCond->AddGroupPattern( kBit0|kBit1|kBit2, 2); // at least 2 hits in 3 inner layers
    trCond->AddGroupPattern( kBit3|kBit4      , 1); // at least 1 hit in 2 middle layers
    trCond->AddGroupPattern( kBit5|kBit6      , 1); // at least 1 hit in 2 outer layers
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
    for (int i=0;i<nLr;i++) itsRecoParam->SetAllowDiagonalClusterization(i,kAllowDiagCl);
    //  
    // Add tracking conditions >>>
    trCond = new AliITSUTrackCond();
    trCond->SetNLayers(nLr); 
    trCond->SetMaxITSTPCMatchChi2(maxChi2Match);
    trCond->SetMaxITSSAChi2(maxChi2SA);
    //
    for (int i=0;i<nLr;i++) {
      trCond->SetMaxBranches(i,nBranch[i]);    // each seed propagated to given layer can produce max nBranch branches
      trCond->SetMaxCandidates(i,nCands[i]);   // each tpc track may have at most nCands prolongations
      trCond->SetMaxTr2ClChi2(i,tr2clChi2[i]);   // cut on cluster to track chi2
      trCond->SetMaxChi2GloNrm(i,gloChi2[i]);  // cut on cluster to track global chi2
      trCond->SetMissPenalty(i,missPen[i]);    // missing cluster penalty
      //
    }
    //
    trCond->AddNewCondition(5); // min hits
    trCond->AddGroupPattern( kBit0|kBit1|kBit2, 2); // at least 2 hits in 3 inner layers
    trCond->AddGroupPattern( kBit3|kBit4      , 1); // at least 1 hit in 2 middle layers
    trCond->AddGroupPattern( kBit5|kBit6      , 1); // at least 1 hit in 2 outer layers
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

