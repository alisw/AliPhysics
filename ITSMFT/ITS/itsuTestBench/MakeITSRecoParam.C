void MakeITSRecoParam(AliRecoParam::EventSpecie_t default=AliRecoParam::kLowMult, const char* cdbURI="local://",
		      int defSATracker=2,Bool_t saOnly=kFALSE) {
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
  gSystem->Load("libITSUpgradeBase");
  gSystem->Load("libITSUpgradeSim");
  gSystem->Load("libITSUpgradeRec");
  //
  // Activate CDB storage and load geometry from CDB
  AliCDBManager* cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage(cdbURI);
  AliITSUTrackCond* trcond = 0;
  int nLr = 7;
  
  TObjArray *recoParamArray = new TObjArray();
  //
  {
    AliITSURecoParam * itsRecoParam = AliITSURecoParam::GetCosmicTestParam(kFALSE);
    //
    itsRecoParam->SetNLayers(nLr);
    //
    //******************************************************************
    //    itsRecoParam->SetEventSpecie(AliRecoParam::kCosmic);
    //    itsRecoParam->SetTitle("Cosmic");
    itsRecoParam->SetTracker(defSATracker);
    itsRecoParam->SetSAonly(saOnly);
    itsRecoParam->SetMaxROCycle(126); // AliITSUSimulation::kMaxROCycleAccept
    recoParamArray->AddLast(itsRecoParam);
  }
  //
  const Bool_t kAllowDiagCl = kFALSE;
  const Bool_t kUseLUT[3] = {kTRUE,kTRUE,kFALSE}; // use TGeo mat.queries only for RefitInward
  //
  // long tracks
  /*
  // tuned for 20x20 pixels with x/x0= 0.3,0.3,0.3,0.5,0.5,0.5,0.5% CDR setup
  //
  int   c0nBranch[7]   = {5,10,15,4,6,6,10}; // max branching for the seed on layer
  int   c0nCands[7]    = {10,20,45,20,45,15,10}; // max candidates for the TPC seed
  float c0tr2clChi2[7] = {20,25,30,40,45,45,70}; // cut on cluster to track chi2 
  float c0gloChi2[7]   = {6, 10,20,40,70,70,70}; // cut on seed global norm chi2
  float c0missPen[7]   = {2.,2.,2.,2.,2.,2.,2.};    // missing cluster penalty
  float c0maxChi2SA[14] = {0.,0.,0.,0.,2.,3.,8., 10.,10.,10.,10.,10.,10.,10.};   // chi2SA vs Nclus
  float c0maxChi2Match = 10.;
  */
  //
  //  /*
  // tuned for 20x20 pixels with x/x0= 0.3,0.3,0.3,0.5,0.5,0.5,0.5% TDR5 setup
  int   c0nBranch[7] = {3,9,15,4,5,7,10}; // max branching for the seed on layer
  int   c0nCands[7]  = {10,15,45,20,60,20,10}; // max candidates for the TPC seed
  float c0tr2clChi2[7] = {20,25,30,40,45,45,70}; // cut on cluster to track chi2 
  float c0gloChi2[7]   = {6,10,20,30,60,60,70}; // cut on seed global norm chi2
  float c0missPen[7] = {2.,2.,2.,2.,2.,2.,2.};    // missing cluster penalty
  float c0maxChi2SA[14] = {0.,0.,0.,0.,2.5,5.,10., 20.,20.,20.,20.,20.,20.,20.};   // chi2SA vs Nclus
  float c0maxChi2Match = 10.;
  //  */


  // short tracks from decays
  /*
  // tuned for 20x20 pixels with x/x0= 0.3,0.3,0.3,0.5,0.5,0.5,0.5% CDR setup
  int   c1nBranch[7]   = {0,0,0,4,6,6,10}; // max branching for the seed on layer
  int   c1nCands[7]    = {0,0,0,20,45,15,10}; // max candidates for the TPC seed
  float c1tr2clChi2[7] = {0,0,0,20,20,20,30}; // cut on cluster to track chi2 
  float c1gloChi2[7]   = {0,0,0,16,40,30,30}; // cut on seed global norm chi2
  float c1missPen[7]   = {0.,0.,0.,2.,2.,2.,2.};    // missing cluster penalty
  float c1maxChi2SA[14] = {0.,0.,0.,7.,8.,8.,8., 10.,10.,10.,10.,10.,10.,10.};  // chi2SA vs Nclus
  float c1maxChi2Match = 8.;
  */

  // short tracks from decays
  int   c1nBranch[7]   = {0,0,0,4,6,6,10}; // max branching for the seed on layer
  int   c1nCands[7]    = {0,0,0,5,5,5,8}; // max candidates for the TPC seed
  float c1tr2clChi2[7] = {0,0,0,20,20,20,30}; // cut on cluster to track chi2 
  float c1gloChi2[7]   = {0,0,0,16,40,35,30}; // cut on seed global norm chi2
  float c1missPen[7]   = {0.,0.,0.,2.,2.,2.,2.};    // missing cluster penalty
  float c1maxChi2SA[14] = {0.,0.,0.,5.,13.,13.,18., 10.,10.,10.,10.,10.,10.,10.};  // chi2SA vs Nclus
  float c1maxChi2Match = 10.;

  //
  /*
    // tuned for 20x33 pixels
  int   c0nBranch[7] = {5,10,15,4,6,6,10}; // max branching for the seed on layer
  int   c0nCands[7]  = {10,20,25,20,30,15,20}; // max candidates for the TPC seed
  float c0tr2clChi2[7] = {20,20,25,25,25,30,40}; // cut on cluster to track chi2 
  float c0gloChi2[7]   = {9, 10,15,20,30,30,30}; // cut on seed global norm chi2
  float c0missPen[7] = {2.,2.,2.,2.,2.,2.,2.};    // missing cluster penalty
  float c0maxChi2SA[14] = {0.,0.,0.,0.,12.,13.,18., 20.,20.,20.,20.,20.,20.,20.};   // chi2SA vs Nclus
  float c0maxChi2Match = 10.;
  */
  // very short tracks from decays
  int   c2nBranch[7]   = {0,0,0,0,0,6,10}; // max branching for the seed on layer
  int   c2nCands[7]    = {0,0,0,0,0,5,8}; // max candidates for the TPC seed
  float c2tr2clChi2[7] = {0,0,0,0,0,15,20}; // cut on cluster to track chi2
  float c2gloChi2[7]   = {0,0,0,0,0,15,20}; // cut on seed global norm chi2
  float c2missPen[7]   = {0.,0.,0.,0.,0.,2.,2.};    // missing cluster penalty
  float c2maxChi2SA[14] = {0.,5.,5.,5.,13.,13.,18., 10.,10.,10.,10.,10.,10.,10.};  // chi2SA vs Nclus, meaningless for 2 point tracks 
  float c2maxChi2Match = 6.;


  //
  {
    AliITSURecoParam * itsRecoParam = AliITSURecoParam::GetLowFluxParam(kFALSE);
    //
    itsRecoParam->SetNLayers(nLr);
    //
    //******************************************************************
    //    itsRecoParam->SetEventSpecie(AliRecoParam::kLowMult);
    //    itsRecoParam->SetTitle("LowMult");
    itsRecoParam->SetTracker(defSATracker);
    itsRecoParam->SetSAonly(saOnly);
    itsRecoParam->SetMaxROCycle(126); // AliITSUSimulation::kMaxROCycleAccept
    recoParamArray->AddLast(itsRecoParam);
    //******************************************************************
    for (int i=0;i<nLr;i++) itsRecoParam->SetAllowDiagonalClusterization(i,kAllowDiagCl);
    for (int i=AliITSURecoParam::kNTrackingPhases;i--;) itsRecoParam->SetUseMatLUT(i,kUseLUT[i]);
    //  
    // Add tracking conditions >>>
    trCond = new AliITSUTrackCond();
    trCond->SetNLayers(nLr); 
    trCond->SetMaxITSTPCMatchChi2(c0maxChi2Match);
    //
    // to exclude some layer use trCon->ExcludeLayer(lrID);
    //
    for (int i=0;i<nLr;i++) {
      trCond->SetMaxBranches(i,c0nBranch[i]);    // each seed propagated to given layer can produce max nBranch branches
      trCond->SetMaxCandidates(i,c0nCands[i]);   // each tpc track may have at most nCands prolongations
      trCond->SetMaxTr2ClChi2(i,c0tr2clChi2[i]); // cut on cluster to track chi2
      trCond->SetMaxChi2GloNrm(i,c0gloChi2[i]);  // cut on cluster to track global chi2
      trCond->SetMissPenalty(i,c0missPen[i]);    // missing cluster penalty
    }
    for (int i=1;i<=2*nLr;i++) trCond->SetMaxITSSAChi2(i,c0maxChi2SA[i-1]);
    //
    trCond->AddNewCondition(5); // min hits
    trCond->AddGroupPattern( kBit0|kBit1|kBit2, 2); // at least 2 hits in 3 inner layers
    trCond->AddGroupPattern( kBit3|kBit4      , 1); // at least 1 hit in 2 middle layers
    trCond->AddGroupPattern( kBit5|kBit6      , 1); // at least 1 hit in 2 outer layers
    //
    trCond->Init();
    //
    itsRecoParam->AddTrackingCondition(trCond);
    //-----------------------------------------------------------
    // short tracks
    trCond = new AliITSUTrackCond();
    trCond->SetNLayers(nLr); 
    //
    trCond->ExcludeLayer(0);
    trCond->ExcludeLayer(1);
    trCond->ExcludeLayer(2);
    //
    trCond->SetMaxITSTPCMatchChi2(c1maxChi2Match);
    //
    // to exclude some layer use trCon->ExcludeLayer(lrID);
    //
    for (int i=0;i<nLr;i++) {
      trCond->SetMaxBranches(i,c1nBranch[i]);    // each seed propagated to given layer can produce max nBranch branches
      trCond->SetMaxCandidates(i,c1nCands[i]);   // each tpc track may have at most nCands prolongations
      trCond->SetMaxTr2ClChi2(i,c1tr2clChi2[i]); // cut on cluster to track chi2
      trCond->SetMaxChi2GloNrm(i,c1gloChi2[i]);  // cut on cluster to track global chi2
      trCond->SetMissPenalty(i,c1missPen[i]);    // missing cluster penalty
    }
    for (int i=1;i<=2*nLr;i++) trCond->SetMaxITSSAChi2(i,c1maxChi2SA[i-1]);
    //
    trCond->AddNewCondition(4); // min hits
    trCond->AddGroupPattern( kBit3|kBit4|kBit5|kBit6, 4); // at least 1 hit in 2 outer layers
    //
    trCond->Init();
    //
    itsRecoParam->AddTrackingCondition(trCond); 
    // Add tracking conditions <<<
    //-----------------------------------------------------------
    // very short tracks
    trCond = new AliITSUTrackCond();
    trCond->SetNLayers(nLr); 
    //
    trCond->ExcludeLayer(0);
    trCond->ExcludeLayer(1);
    trCond->ExcludeLayer(2);
    trCond->ExcludeLayer(3);
    trCond->ExcludeLayer(4);
    //
    trCond->SetMaxITSTPCMatchChi2(c2maxChi2Match);
    //
    // to exclude some layer use trCon->ExcludeLayer(lrID);
    //
    for (int i=0;i<nLr;i++) {
      trCond->SetMaxBranches(i,c2nBranch[i]);    // each seed propagated to given layer can produce max nBranch branches
      trCond->SetMaxCandidates(i,c2nCands[i]);   // each tpc track may have at most nCands prolongations
      trCond->SetMaxTr2ClChi2(i,c2tr2clChi2[i]); // cut on cluster to track chi2
      trCond->SetMaxChi2GloNrm(i,c2gloChi2[i]);  // cut on cluster to track global chi2
      trCond->SetMissPenalty(i,c2missPen[i]);    // missing cluster penalty
    }
    for (int i=1;i<=2*nLr;i++) trCond->SetMaxITSSAChi2(i,c2maxChi2SA[i-1]);
    //
    trCond->AddNewCondition(2); // min hits
    trCond->AddGroupPattern( kBit5|kBit6, 2);
    //
    trCond->Init();
    //
    itsRecoParam->AddTrackingCondition(trCond);
    // Add tracking conditions <<<
  }
  {
    AliITSURecoParam * itsRecoParam = AliITSURecoParam::GetHighFluxParam(kFALSE);
    //
    itsRecoParam->SetNLayers(nLr);
    //
    //******************************************************************
    //    itsRecoParam->SetEventSpecie(AliRecoParam::kHighMult);
    //    itsRecoParam->SetTitle("HighMult");
    itsRecoParam->SetTracker(defSATracker);
    itsRecoParam->SetSAonly(saOnly);
    itsRecoParam->SetMaxROCycle(126); // AliITSUSimulation::kMaxROCycleAccept
    recoParamArray->AddLast(itsRecoParam);
    //******************************************************************
    for (int i=0;i<nLr;i++) itsRecoParam->SetAllowDiagonalClusterization(i,kAllowDiagCl);
    for (int i=AliITSURecoParam::kNTrackingPhases;i--;) itsRecoParam->SetUseMatLUT(i,kUseLUT[i]);
    //  
    // Add tracking conditions >>>
    trCond = new AliITSUTrackCond();
    trCond->SetNLayers(nLr); 
    trCond->SetMaxITSTPCMatchChi2(c0maxChi2Match);
    //
    for (int i=0;i<nLr;i++) {
      trCond->SetMaxBranches(i,c0nBranch[i]);    // each seed propagated to given layer can produce max nBranch branches
      trCond->SetMaxCandidates(i,c0nCands[i]);   // each tpc track may have at most nCands prolongations
      trCond->SetMaxTr2ClChi2(i,c0tr2clChi2[i]);   // cut on cluster to track chi2
      trCond->SetMaxChi2GloNrm(i,c0gloChi2[i]);  // cut on cluster to track global chi2
      trCond->SetMissPenalty(i,c0missPen[i]);    // missing cluster penalty
    }
    for (int i=1;i<=2*nLr;i++) trCond->SetMaxITSSAChi2(i,c0maxChi2SA[i-1]);
    //
    trCond->AddNewCondition(5); // min hits
    trCond->AddGroupPattern( kBit0|kBit1|kBit2, 2); // at least 2 hits in 3 inner layers
    trCond->AddGroupPattern( kBit3|kBit4      , 1); // at least 1 hit in 2 middle layers
    trCond->AddGroupPattern( kBit5|kBit6      , 1); // at least 1 hit in 2 outer layers
    //
    trCond->Init();
    //
    itsRecoParam->AddTrackingCondition(trCond);
    //-----------------------------------------------------------
    // short tracks
    trCond = new AliITSUTrackCond();
    trCond->SetNLayers(nLr); 
    //
    trCond->ExcludeLayer(0);
    trCond->ExcludeLayer(1);
    trCond->ExcludeLayer(2);
    //
    trCond->SetMaxITSTPCMatchChi2(c1maxChi2Match);
    //
    // to exclude some layer use trCon->ExcludeLayer(lrID);
    //
    for (int i=0;i<nLr;i++) {
      trCond->SetMaxBranches(i,c1nBranch[i]);    // each seed propagated to given layer can produce max nBranch branches
      trCond->SetMaxCandidates(i,c1nCands[i]);   // each tpc track may have at most nCands prolongations
      trCond->SetMaxTr2ClChi2(i,c1tr2clChi2[i]); // cut on cluster to track chi2
      trCond->SetMaxChi2GloNrm(i,c1gloChi2[i]);  // cut on cluster to track global chi2
      trCond->SetMissPenalty(i,c1missPen[i]);    // missing cluster penalty
    }
    for (int i=1;i<=2*nLr;i++) trCond->SetMaxITSSAChi2(i,c1maxChi2SA[i-1]);
    //
    trCond->AddNewCondition(4); // min hits
    trCond->AddGroupPattern( kBit3|kBit4|kBit5|kBit6, 4);
    //
    trCond->Init();
    //
    itsRecoParam->AddTrackingCondition(trCond);
    // Add tracking conditions <<<
    //
    //-----------------------------------------------------------
    // very short tracks
    trCond = new AliITSUTrackCond();
    trCond->SetNLayers(nLr); 
    //
    trCond->ExcludeLayer(0);
    trCond->ExcludeLayer(1);
    trCond->ExcludeLayer(2);
    trCond->ExcludeLayer(3);
    trCond->ExcludeLayer(4);
    //
    trCond->SetMaxITSTPCMatchChi2(c2maxChi2Match);
    //
    // to exclude some layer use trCon->ExcludeLayer(lrID);
    //
    for (int i=0;i<nLr;i++) {
      trCond->SetMaxBranches(i,c2nBranch[i]);    // each seed propagated to given layer can produce max nBranch branches
      trCond->SetMaxCandidates(i,c2nCands[i]);   // each tpc track may have at most nCands prolongations
      trCond->SetMaxTr2ClChi2(i,c2tr2clChi2[i]); // cut on cluster to track chi2
      trCond->SetMaxChi2GloNrm(i,c2gloChi2[i]);  // cut on cluster to track global chi2
      trCond->SetMissPenalty(i,c2missPen[i]);    // missing cluster penalty
    }
    for (int i=1;i<=2*nLr;i++) trCond->SetMaxITSSAChi2(i,c2maxChi2SA[i-1]);
    //
    trCond->AddNewCondition(2); // min hits
    trCond->AddGroupPattern( kBit5|kBit6, 2);
    //
    trCond->Init();
    //
    itsRecoParam->AddTrackingCondition(trCond);
    // Add tracking conditions <<<
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

