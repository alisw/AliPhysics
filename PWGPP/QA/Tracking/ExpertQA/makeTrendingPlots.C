void SetGraphProperties( TGraphErrors *gr, const char *title, const char *yAxisTitle, Int_t color, Int_t mStyle, Float_t mSize);
void DrawAndSave( TGraphErrors *gr, const char *name );

makeTrendingPlots( const char *TrendingFile )
{
  TFile *f = TFile::Open(TrendingFile);
  if(!f){ std::cout << "ERROR! No Trendingfile given!" <<std::endl; return; }

  TTree *tree = (TTree*) f->Get("trending");
  if(!tree){ std::cout << "ERROR! No Tree available!" <<std::endl; return; }

  //TGraph* MakeGraphSparse(TTree* tree, const char* expr = "Entry", const char* cut = "1", Int_t mstyle = 25, Int_t mcolor = 1, Float_t msize = -1, Float_t offset = 0.0)

  // K0 Res
  // shift
  TGraphErrors *grK0shiftResPosHigh1pt = TStatToolkit::MakeGraphSparse( tree, "shiftK0sResPosHigh1pt:run:eShiftK0sResPosHigh1pt","shiftK0sResPosHigh1pt>-900" );
  SetGraphProperties(grK0shiftResPosHigh1pt,"K0 shift resolution (positive tracks, 1pt = 1)","K0 shift resolution",1,20,.6);
  TGraphErrors *grK0shiftResNegHigh1pt = TStatToolkit::MakeGraphSparse( tree, "shiftK0sResNegHigh1pt:run:eShiftK0sResNegHigh1pt","shiftK0sResNegHigh1pt>-900" );
  SetGraphProperties(grK0shiftResNegHigh1pt,"K0 shift resolution (negative tracks, 1pt = 1)","K0 shift resolution",1,20,.6);
  TGraphErrors *grK0shiftResPosLow1pt  = TStatToolkit::MakeGraphSparse( tree, "shiftK0sResPosLow1pt:run:eShiftK0sResPosLow1pt"  ,"shiftK0sResPosLow1pt>-900" );
  SetGraphProperties(grK0shiftResPosLow1pt,"K0 shift resolution (positive tracks, 1pt = 0)","K0 shift resolution",1,20,.6);
  TGraphErrors *grK0shiftResNegLow1pt  = TStatToolkit::MakeGraphSparse( tree, "shiftK0sResNegLow1pt:run:eShiftK0sResNegLow1pt"  ,"shiftK0sResNegLow1pt>-900" );
  SetGraphProperties(grK0shiftResNegLow1pt,"K0 shift resolution (negative tracks, 1pt = 0)","K0 shift resolution",1,20,.6);
  // sigma
  TGraphErrors *grK0sigmaResPosHigh1pt = TStatToolkit::MakeGraphSparse( tree, "sigmaK0sResPosHigh1pt:run:eSigmaK0sResPosHigh1pt","sigmaK0sResPosHigh1pt>-900" );
  SetGraphProperties(grK0sigmaResPosHigh1pt,"K0 sigma resolution (positive tracks, 1pt = 1)","K0 sigma resolution",1,20,.6);
  TGraphErrors *grK0sigmaResNegHigh1pt = TStatToolkit::MakeGraphSparse( tree, "sigmaK0sResNegHigh1pt:run:eSigmaK0sResNegHigh1pt","sigmaK0sResNegHigh1pt>-900" );
  SetGraphProperties(grK0sigmaResNegHigh1pt,"K0 sigma resolution (negative tracks, 1pt = 1)","K0 sigma resolution",1,20,.6);
  TGraphErrors *grK0sigmaResPosLow1pt  = TStatToolkit::MakeGraphSparse( tree, "sigmaK0sResPosLow1pt:run:eSigmaK0sResPosLow1pt"  ,"sigmaK0sResPosLow1pt>-900" );
  SetGraphProperties(grK0sigmaResPosLow1pt,"K0 sigma resolution (positive tracks, 1pt = 0)","K0 sigma resolution",1,20,.6);
  TGraphErrors *grK0sigmaResNegLow1pt  = TStatToolkit::MakeGraphSparse( tree, "sigmaK0sResNegLow1pt:run:eSigmaK0sResNegLow1pt"  ,"sigmaK0sResNegLow1pt>-900" );
  SetGraphProperties(grK0sigmaResNegLow1pt,"K0 sigma resolution (negative tracks, 1pt = 0)","K0 sigma resolution",1,20,.6);
  // K0 Pull
  // shift
  TGraphErrors *grK0shiftPullPosHigh1pt = TStatToolkit::MakeGraphSparse( tree, "shiftK0sPullPosHigh1pt:run:eShiftK0sPullPosHigh1pt","shiftK0sPullPosHigh1pt>-900" );
  SetGraphProperties(grK0shiftPullPosHigh1pt,"K0 shift pull (positive tracks, 1pt = 1)","K0 shift pull",1,20,.6);
  TGraphErrors *grK0shiftPullNegHigh1pt = TStatToolkit::MakeGraphSparse( tree, "shiftK0sPullNegHigh1pt:run:eShiftK0sPullNegHigh1pt","shiftK0sPullNegHigh1pt>-900" );
  SetGraphProperties(grK0shiftPullNegHigh1pt,"K0 shift pull (negative tracks, 1pt = 1)","K0 shift pull",1,20,.6);
  TGraphErrors *grK0shiftPullPosLow1pt  = TStatToolkit::MakeGraphSparse( tree, "shiftK0sPullPosLow1pt:run:eShiftK0sPullPosLow1pt"  ,"shiftK0sPullPosLow1pt>-900" );
  SetGraphProperties(grK0shiftPullPosLow1pt,"K0 shift pull (positive tracks, 1pt = 0)","K0 shift pull",1,20,.6);
  TGraphErrors *grK0shiftPullNegLow1pt  = TStatToolkit::MakeGraphSparse( tree, "shiftK0sPullNegLow1pt:run:eShiftK0sPullNegLow1pt"  ,"shiftK0sPullNegLow1pt>-900" );
  SetGraphProperties(grK0shiftPullNegLow1pt,"K0 shift pull (negative tracks, 1pt = 0)","K0 shift pull",1,20,.6);
  // sigma
  TGraphErrors *grK0sigmaPullPosHigh1pt = TStatToolkit::MakeGraphSparse( tree, "sigmaK0sPullPosHigh1pt:run:eSigmaK0sPullPosHigh1pt","sigmaK0sPullPosHigh1pt>-900" );
  SetGraphProperties(grK0sigmaPullPosHigh1pt,"K0 sigma pull (positive tracks, 1pt = 1)","K0 sigma pull",1,20,.6);
  TGraphErrors *grK0sigmaPullNegHigh1pt = TStatToolkit::MakeGraphSparse( tree, "sigmaK0sPullNegHigh1pt:run:eSigmaK0sPullNegHigh1pt","sigmaK0sPullNegHigh1pt>-900" );
  SetGraphProperties(grK0sigmaPullNegHigh1pt,"K0 sigma pull (negative tracks, 1pt = 1)","K0 sigma pull",1,20,.6);
  TGraphErrors *grK0sigmaPullPosLow1pt  = TStatToolkit::MakeGraphSparse( tree, "sigmaK0sPullPosLow1pt:run:eSigmaK0sPullPosLow1pt"  ,"sigmaK0sPullPosLow1pt>-900" );
  SetGraphProperties(grK0sigmaPullPosLow1pt,"K0 sigma pull (positive tracks, 1pt = 0)","K0 sigma pull",1,20,.6);
  TGraphErrors *grK0sigmaPullNegLow1pt  = TStatToolkit::MakeGraphSparse( tree, "sigmaK0sPullNegLow1pt:run:eSigmaK0sPullNegLow1pt"  ,"sigmaK0sPullNegLow1pt>-900" );
  SetGraphProperties(grK0sigmaPullNegLow1pt,"K0 sigma pull (negative tracks, 1pt = 0)","K0 sigma pull",1,20,.6);

  //DCAr Res
  // TPC+ITS combined
  TGraphErrors *grDCArResCombinedLow1pt   = TStatToolkit::MakeGraphSparse( tree, "dcaRresCombinedLow1pt:run:edcaRresCombinedLow1pt","dcaRresCombinedLow1pt>-900" );
  SetGraphProperties(grDCArResCombinedLow1pt,"DCAr res (TPC+ITS combined tracking, 1pt = 0)","DCAr res",1,20,.6);
  TGraphErrors *grDCArResCombinedHigh1pt  = TStatToolkit::MakeGraphSparse( tree, "dcaRresCombinedHigh1pt:run:edcaRresCombinedHigh1pt","dcaRresCombinedHigh1pt>-900" );
  SetGraphProperties(grDCArResCombinedHigh1pt,"DCAr res (TPC+ITS combined tracking, 1pt = 1)","DCAr res",1,20,.6);
  // TPC only
  TGraphErrors *grDCArResTPCAsideLow1pt   = TStatToolkit::MakeGraphSparse( tree, "dcaRresTPCAsideLow1pt:run:edcaRresTPCAsideLow1pt","dcaRresTPCAsideLow1pt>-900" );
  SetGraphProperties(grDCArResTPCAsideLow1pt,"DCAr res (TPC only tracking Aside, 1pt = 0)","DCAr res",1,20,.6);
  TGraphErrors *grDCArResTPCAsideHigh1pt  = TStatToolkit::MakeGraphSparse( tree, "dcaRresTPCAsideHigh1pt:run:edcaRresTPCAsideHigh1pt","dcaRresTPCAsideHigh1pt>-900" );
  SetGraphProperties(grDCArResTPCAsideHigh1pt,"DCAr res (TPC only tracking Aside, 1pt = 1)","DCAr res",1,20,.6);
  TGraphErrors *grDCArResTPCCsideLow1pt   = TStatToolkit::MakeGraphSparse( tree, "dcaRresTPCCsideLow1pt:run:edcaRresTPCCsideLow1pt","dcaRresTPCCsideLow1pt>-900" );
  SetGraphProperties(grDCArResTPCCsideLow1pt,"DCAr res (TPC only tracking Cside, 1pt = 0)","DCAr res",1,20,.6);
  TGraphErrors *grDCArResTPCCsideHigh1pt  = TStatToolkit::MakeGraphSparse( tree, "dcaRresTPCCsideHigh1pt:run:edcaRresTPCCsideHigh1pt","dcaRresTPCCsideHigh1pt>-900" );
  SetGraphProperties(grDCArResTPCCsideHigh1pt,"DCAr res (TPC only tracking Cside, 1pt = 1)","DCAr res",1,20,.6);
  //DCAr Pull
  // TPC+ITS combined
  TGraphErrors *grDCArPullCombinedLow1pt  = TStatToolkit::MakeGraphSparse( tree, "dcaRpullCombinedLow1pt:run:edcaRpullCombinedLow1pt","dcaRpullCombinedLow1pt>-900" );
  SetGraphProperties(grDCArPullCombinedLow1pt,"DCAr pull (TPC+ITS combined tracking, 1pt = 0)","DCAr pull",1,20,.6);
  TGraphErrors *grDCArPullCombinedHigh1pt = TStatToolkit::MakeGraphSparse( tree, "dcaRpullCombinedHigh1pt:run:edcaRpullCombinedHigh1pt","dcaRpullCombinedHigh1pt>-900" );
  SetGraphProperties(grDCArPullCombinedHigh1pt,"DCAr pull (TPC+ITS combined tracking, 1pt = 1)","DCAr pull",1,20,.6);
  // TPC only
  TGraphErrors *grDCArPullTPCAsideLow1pt  = TStatToolkit::MakeGraphSparse( tree, "dcaRpullTPCAsideLow1pt:run:edcaRpullTPCAsideLow1pt","dcaRpullTPCAsideLow1pt>-900" );
  SetGraphProperties(grDCArPullTPCAsideLow1pt,"DCAr pull (TPC only tracking Aside, 1pt = 0)","DCAr pull",1,20,.6);
  TGraphErrors *grDCArPullTPCAsideHigh1pt = TStatToolkit::MakeGraphSparse( tree, "dcaRpullTPCAsideHigh1pt:run:edcaRpullTPCAsideHigh1pt","dcaRpullTPCAsideHigh1pt>-900" );
  SetGraphProperties(grDCArPullTPCAsideHigh1pt,"DCAr pull (TPC only tracking Aside, 1pt = 1)","DCAr pull",1,20,.6);
  TGraphErrors *grDCArPullTPCCsideLow1pt  = TStatToolkit::MakeGraphSparse( tree, "dcaRpullTPCCsideLow1pt:run:edcaRpullTPCCsideLow1pt","dcaRpullTPCCsideLow1pt>-900" );
  SetGraphProperties(grDCArPullTPCCsideLow1pt,"DCAr pull (TPC only tracking Cside, 1pt = 0)","DCAr pull",1,20,.6);
  TGraphErrors *grDCArPullTPCCsideHigh1pt = TStatToolkit::MakeGraphSparse( tree, "dcaRpullTPCCsideHigh1pt:run:edcaRpullTPCCsideHigh1pt","dcaRpullTPCCsideHigh1pt>-900" );
  SetGraphProperties(grDCArPullTPCCsideHigh1pt,"DCAr pull (TPC only tracking Cside, 1pt = 1)","DCAr pull",1,20,.6);

  //qoverpt Shift
  // Sin part
  //Aside
  TGraphErrors *grqptShiftCombinedSinAside  = TStatToolkit::MakeGraphSparse( tree, "qptShiftCombinedSinAside:run:eqptShiftCombinedSinAside","qptShiftCombinedSinAside>-900" );
  SetGraphProperties(grqptShiftCombinedSinAside,"qpt shift sin part combined tracking Aside","sin part",1,20,.6);
  TGraphErrors *grqptShiftTPCconstSinAside  = TStatToolkit::MakeGraphSparse( tree, "qptShiftTPCconstSinAside:run:eqptShiftTPCconstSinAside","qptShiftTPCconstSinAside>-900" );
  SetGraphProperties(grqptShiftTPCconstSinAside,"qpt shift sin part TPCconstrained Aside","sin part",1,20,.6);
  TGraphErrors *grqptShiftTPConlySinAside   = TStatToolkit::MakeGraphSparse( tree, "qptShiftTPConlySinAside:run:eqptShiftTPConlySinAside","qptShiftTPConlySinAside>-900" );
  SetGraphProperties(grqptShiftTPConlySinAside,"qpt shift sin part TPConly Aside","sin part",1,20,.6);
  //Cside
  TGraphErrors *grqptShiftCombinedSinCside  = TStatToolkit::MakeGraphSparse( tree, "qptShiftCombinedSinCside:run:eqptShiftCombinedSinCside","qptShiftCombinedSinCside>-900" );
  SetGraphProperties(grqptShiftCombinedSinCside,"qpt shift sin part combined tracking Cside","sin part",1,20,.6);
  TGraphErrors *grqptShiftTPCconstSinCside  = TStatToolkit::MakeGraphSparse( tree, "qptShiftTPCconstSinCside:run:eqptShiftTPCconstSinCside","qptShiftTPCconstSinCside>-900" );
  SetGraphProperties(grqptShiftTPCconstSinCside,"qpt shift sin part TPCconstrained Cside","sin part",1,20,.6);
  TGraphErrors *grqptShiftTPConlySinCside   = TStatToolkit::MakeGraphSparse( tree, "qptShiftTPConlySinCside:run:eqptShiftTPConlySinCside","qptShiftTPConlySinCside>-900" );
  SetGraphProperties(grqptShiftTPConlySinAside,"qpt shift sin part TPConly Aside","sin part",1,20,.6);
  // Cos part
  //Aside
  TGraphErrors *grqptShiftCombinedCosAside  = TStatToolkit::MakeGraphSparse( tree, "qptShiftCombinedCosAside:run:eqptShiftCombinedCosAside","qptShiftCombinedCosAside>-900" );
  SetGraphProperties(grqptShiftCombinedCosAside,"qpt shift cos part combined tracking Aside","cos part",1,20,.6);
  TGraphErrors *grqptShiftTPCconstCosAside  = TStatToolkit::MakeGraphSparse( tree, "qptShiftTPCconstCosAside:run:eqptShiftTPCconstCosAside","qptShiftTPCconstCosAside>-900" );
  SetGraphProperties(grqptShiftTPCconstCosAside,"qpt shift cos part TPCconstrained Aside","cos part",1,20,.6);
  TGraphErrors *grqptShiftTPConlyCosAside   = TStatToolkit::MakeGraphSparse( tree, "qptShiftTPConlyCosAside:run:eqptShiftTPConlyCosAside","qptShiftTPConlyCosAside>-900" );
  SetGraphProperties(grqptShiftTPConlyCosAside,"qpt shift cos part TPConly Aside","cos part",1,20,.6);
  //Cside
  TGraphErrors *grqptShiftCombinedCosCside  = TStatToolkit::MakeGraphSparse( tree, "qptShiftCombinedCosCside:run:eqptShiftCombinedCosCside","qptShiftCombinedCosCside>-900" );
  SetGraphProperties(grqptShiftCombinedCosCside,"qpt shift cos part combined tracking Cside","cos part",1,20,.6);
  TGraphErrors *grqptShiftTPCconstCosCside  = TStatToolkit::MakeGraphSparse( tree, "qptShiftTPCconstCosCside:run:eqptShiftTPCconstCosCside","qptShiftTPCconstCosCside>-900" );
  SetGraphProperties(grqptShiftTPCconstCosCside,"qpt shift cos part TPCconstrained Cside","cos part",1,20,.6);
  TGraphErrors *grqptShiftTPConlyCosCside   = TStatToolkit::MakeGraphSparse( tree, "qptShiftTPConlyCosCside:run:eqptShiftTPConlyCosCside","qptShiftTPConlyCosCside>-900" );
  SetGraphProperties(grqptShiftTPConlyCosCside,"qpt shift cos part TPConly Cside","cos part",1,20,.6);

  //delta Phi
  // res
  //Aside
  TGraphErrors *grdPhiResTPCAsideLow1pt  = TStatToolkit::MakeGraphSparse( tree, "dPhiResTPCAsideLow1pt:run:edPhiResTPCAsideLow1pt","dPhiResTPCAsideLow1pt>-900" );
  SetGraphProperties(grdPhiResTPCAsideLow1pt,"delta Phi resolution (Aside, 1pt = 0)","phi res",1,20,.6);
  TGraphErrors *grdPhiResTPCAsideHigh1pt = TStatToolkit::MakeGraphSparse( tree, "dPhiResTPCAsideHigh1pt:run:edPhiResTPCAsideHigh1pt","dPhiResTPCAsideHigh1pt>-900" );
  SetGraphProperties(grdPhiResTPCAsideHigh1pt,"delta Phi resolution (Aside, 1pt = 1)","phi res",1,20,.6);
  //Cside
  TGraphErrors *grdPhiResTPCCsideLow1pt  = TStatToolkit::MakeGraphSparse( tree, "dPhiResTPCCsideLow1pt:run:edPhiResTPCCsideLow1pt","dPhiResTPCCsideLow1pt>-900" );
  SetGraphProperties(grdPhiResTPCCsideLow1pt,"delta Phi resolution (Cside, 1pt = 0)","phi res",1,20,.6);
  TGraphErrors *grdPhiResTPCCsideHigh1pt = TStatToolkit::MakeGraphSparse( tree, "dPhiResTPCCsideHigh1pt:run:edPhiResTPCCsideHigh1pt","dPhiResTPCCsideHigh1pt>-900" );
  SetGraphProperties(grdPhiResTPCCsideHigh1pt,"delta Phi resolution (Cside, 1pt = 1)","phi res",1,20,.6);
  // pull
  //Aside
  TGraphErrors *grdPhiPullTPCAsideLow1pt  = TStatToolkit::MakeGraphSparse( tree, "dPhiPullTPCAsideLow1pt:run:edPhiPullTPCAsideLow1pt","dPhiPullTPCAsideLow1pt>-900" );
  SetGraphProperties(grdPhiPullTPCAsideLow1pt,"delta Phi pull (Aside, 1pt = 0)","phi pull",1,20,.6);
  TGraphErrors *grdPhiPullTPCAsideHigh1pt = TStatToolkit::MakeGraphSparse( tree, "dPhiPullTPCAsideHigh1pt:run:edPhiPullTPCAsideHigh1pt","dPhiPullTPCAsideHigh1pt>-900" );
  SetGraphProperties(grdPhiPullTPCAsideHigh1pt,"delta Phi pull (Aside, 1pt = 1)","phi pull",1,20,.6);
  //Cside
  TGraphErrors *grdPhiPullTPCCsideLow1pt  = TStatToolkit::MakeGraphSparse( tree, "dPhiPullTPCCsideLow1pt:run:edPhiPullTPCCsideLow1pt","dPhiPullTPCCsideLow1pt>-900" );
  SetGraphProperties(grdPhiPullTPCCsideLow1pt,"delta Phi pull (Cside, 1pt = 0)","phi pull",1,20,.6);
  TGraphErrors *grdPhiPullTPCCsideHigh1pt = TStatToolkit::MakeGraphSparse( tree, "dPhiPullTPCCsideHigh1pt:run:edPhiPullTPCCsideHigh1pt","dPhiPullTPCCsideHigh1pt>-900" );
  SetGraphProperties(grdPhiPullTPCCsideHigh1pt,"delta Phi pull (Cside, 1pt = 1)","phi pull",1,20,.6);

  // TPCITS matching Efficiency
  TGraphErrors *grEfficiencyLowPt   = TStatToolkit::MakeGraphSparse( tree, "EfficiencyLowPt:run:eEfficiencyLowPt","EfficiencyLowPt>-900" );
  SetGraphProperties(grEfficiencyLowPt,"TPCITS matching Efficiency (p_{T} < 1 GeV)","Efficiency",1,20,.6);
  grEfficiencyLowPt->GetYaxis()->SetRangeUser(0.,1.1);
  TGraphErrors *grEfficiencyHighPt  = TStatToolkit::MakeGraphSparse( tree, "EfficiencyHighPt:run:eEfficiencyHighPt","EfficiencyHighPt>-900" );
  SetGraphProperties(grEfficiencyHighPt,"TPCITS matching Efficiency (p_{T} > 4 GeV)","Efficiency",1,20,.6);
  grEfficiencyHighPt->GetYaxis()->SetRangeUser(0.,1.1);



  if(grK0shiftResPosHigh1pt     ) DrawAndSave( grK0shiftResPosHigh1pt,     "K0shiftResPosHigh1pt"    );
  if(grK0shiftResNegHigh1pt     ) DrawAndSave( grK0shiftResNegHigh1pt,     "K0shiftResNegHigh1pt"    );
  if(grK0shiftResPosLow1pt      ) DrawAndSave( grK0shiftResPosLow1pt,      "K0shiftResPosLow1pt"     );
  if(grK0shiftResNegLow1pt      ) DrawAndSave( grK0shiftResNegLow1pt,      "K0shiftResNegLow1pt"     );
  if(grK0sigmaResPosHigh1pt     ) DrawAndSave( grK0sigmaResPosHigh1pt,     "K0sigmaResPosHigh1pt"    );
  if(grK0sigmaResNegHigh1pt     ) DrawAndSave( grK0sigmaResNegHigh1pt,     "K0sigmaResNegHigh1pt"    );
  if(grK0sigmaResPosLow1pt      ) DrawAndSave( grK0sigmaResPosLow1pt,      "K0sigmaResPosLow1pt"     );
  if(grK0sigmaResNegLow1pt      ) DrawAndSave( grK0sigmaResNegLow1pt,      "K0sigmaResNegLow1pt"     );
  if(grK0shiftPullPosHigh1pt    ) DrawAndSave( grK0shiftPullPosHigh1pt,    "K0shiftPullPosHigh1pt"   );
  if(grK0shiftPullNegHigh1pt    ) DrawAndSave( grK0shiftPullNegHigh1pt,    "K0shiftPullNegHigh1pt"   );
  if(grK0shiftPullPosLow1pt     ) DrawAndSave( grK0shiftPullPosLow1pt,     "K0shiftPullPosLow1pt"    );
  if(grK0shiftPullNegLow1pt     ) DrawAndSave( grK0shiftPullNegLow1pt,     "K0shiftPullNegLow1pt"    );
  if(grK0sigmaPullPosHigh1pt    ) DrawAndSave( grK0sigmaPullPosHigh1pt,    "K0sigmaPullPosHigh1pt"   );
  if(grK0sigmaPullNegHigh1pt    ) DrawAndSave( grK0sigmaPullNegHigh1pt,    "K0sigmaPullNegHigh1pt"   );
  if(grK0sigmaPullPosLow1pt     ) DrawAndSave( grK0sigmaPullPosLow1pt,     "K0sigmaPullPosLow1pt"    );
  if(grK0sigmaPullNegLow1pt     ) DrawAndSave( grK0sigmaPullNegLow1pt,     "K0sigmaPullNegLow1pt"    );
  if(grDCArResCombinedLow1pt    ) DrawAndSave( grDCArResCombinedLow1pt,    "DCArResCombinedLow1pt"   );
  if(grDCArResCombinedHigh1pt   ) DrawAndSave( grDCArResCombinedHigh1pt,   "DCArResCombinedHigh1pt"  );
  if(grDCArResTPCAsideLow1pt    ) DrawAndSave( grDCArResTPCAsideLow1pt,    "DCArResTPCAsideLow1pt"   );
  if(grDCArResTPCAsideHigh1pt   ) DrawAndSave( grDCArResTPCAsideHigh1pt,   "DCArResTPCAsideHigh1pt"  );
  if(grDCArResTPCCsideLow1pt    ) DrawAndSave( grDCArResTPCCsideLow1pt,    "DCArResTPCCsideLow1pt"   );
  if(grDCArResTPCCsideHigh1pt   ) DrawAndSave( grDCArResTPCCsideHigh1pt,   "DCArResTPCCsideHigh1pt"  );
  if(grDCArPullCombinedLow1pt   ) DrawAndSave( grDCArPullCombinedLow1pt,   "DCArPullCombinedLow1pt"  );
  if(grDCArPullCombinedHigh1pt  ) DrawAndSave( grDCArPullCombinedHigh1pt,  "DCArPullCombinedHigh1pt" );
  if(grDCArPullTPCAsideLow1pt   ) DrawAndSave( grDCArPullTPCAsideLow1pt,   "DCArPullTPCAsideLow1pt"  );
  if(grDCArPullTPCAsideHigh1pt  ) DrawAndSave( grDCArPullTPCAsideHigh1pt,  "DCArPullTPCAsideHigh1pt" );
  if(grDCArPullTPCCsideLow1pt   ) DrawAndSave( grDCArPullTPCCsideLow1pt,   "DCArPullTPCCsideLow1pt"  );
  if(grDCArPullTPCCsideHigh1pt  ) DrawAndSave( grDCArPullTPCCsideHigh1pt,  "DCArPullTPCCsideHigh1pt" );
  if(grqptShiftCombinedSinAside ) DrawAndSave( grqptShiftCombinedSinAside, "qptShiftCombinedSinAside");
  if(grqptShiftTPCconstSinAside ) DrawAndSave( grqptShiftTPCconstSinAside, "qptShiftTPCconstSinAside");
  if(grqptShiftTPConlySinAside  ) DrawAndSave( grqptShiftTPConlySinAside,  "qptShiftTPConlySinAside" );
  if(grqptShiftCombinedSinCside ) DrawAndSave( grqptShiftCombinedSinCside, "qptShiftCombinedSinCside");
  if(grqptShiftTPCconstSinCside ) DrawAndSave( grqptShiftTPCconstSinCside, "qptShiftTPCconstSinCside");
  if(grqptShiftTPConlySinCside  ) DrawAndSave( grqptShiftTPConlySinCside,  "qptShiftTPConlySinCside" );
  if(grqptShiftCombinedCosAside ) DrawAndSave( grqptShiftCombinedCosAside, "qptShiftCombinedCosAside");
  if(grqptShiftTPCconstCosAside ) DrawAndSave( grqptShiftTPCconstCosAside, "qptShiftTPCconstCosAside");
  if(grqptShiftTPConlyCosAside  ) DrawAndSave( grqptShiftTPConlyCosAside,  "qptShiftTPConlyCosAside" );
  if(grqptShiftCombinedCosCside ) DrawAndSave( grqptShiftCombinedCosCside, "qptShiftCombinedCosCside");
  if(grqptShiftTPCconstCosCside ) DrawAndSave( grqptShiftTPCconstCosCside, "qptShiftTPCconstCosCside");
  if(grqptShiftTPConlyCosCside  ) DrawAndSave( grqptShiftTPConlyCosCside,  "qptShiftTPConlyCosCside" );
  if(grdPhiResTPCAsideLow1pt    ) DrawAndSave( grdPhiResTPCAsideLow1pt,    "dPhiResTPCAsideLow1pt"   );
  if(grdPhiResTPCAsideHigh1pt   ) DrawAndSave( grdPhiResTPCAsideHigh1pt,   "dPhiResTPCAsideHigh1pt"  );
  if(grdPhiResTPCCsideLow1pt    ) DrawAndSave( grdPhiResTPCCsideLow1pt,    "dPhiResTPCCsideLow1pt"   );
  if(grdPhiResTPCCsideHigh1pt   ) DrawAndSave( grdPhiResTPCCsideHigh1pt,   "dPhiResTPCCsideHigh1pt"  );
  if(grdPhiPullTPCAsideLow1pt   ) DrawAndSave( grdPhiPullTPCAsideLow1pt,   "dPhiPullTPCAsideLow1pt"  );
  if(grdPhiPullTPCAsideHigh1pt  ) DrawAndSave( grdPhiPullTPCAsideHigh1pt,  "dPhiPullTPCAsideHigh1pt" );
  if(grdPhiPullTPCCsideLow1pt   ) DrawAndSave( grdPhiPullTPCCsideLow1pt,   "dPhiPullTPCCsideLow1pt"  );
  if(grdPhiPullTPCCsideHigh1pt  ) DrawAndSave( grdPhiPullTPCCsideHigh1pt,  "dPhiPullTPCCsideHigh1pt" );
  if(grEfficiencyLowPt          ) DrawAndSave( grEfficiencyLowPt,          "EfficiencyLowPt"         );
  if(grEfficiencyHighPt         ) DrawAndSave( grEfficiencyHighPt,         "EfficiencyHighPt"        );



}

void SetGraphProperties( TGraphErrors *gr, const char *title, const char *yAxisTitle, Int_t color, Int_t mStyle, Float_t mSize){
  gr->SetTitle(title);
  gr->GetXaxis()->SetTitle("run");
  gr->GetYaxis()->SetTitle(yAxisTitle);
  gr->SetMarkerStyle(mStyle);
  gr->SetMarkerSize(mSize);
  gr->SetLineColor(color);
  gr->SetMarkerColor(color);
}

void DrawAndSave( TGraphErrors *gr, const char *name ){

  TCanvas *can = new TCanvas("can","testing",1200,800);
  can->cd();

  gr->Draw("ap");

  //gSystem->Exec("if [ ! -d ./TrendingPlots ] ; then mkdir -p TrendingPlots ; fi");

  //can->SaveAs( Form("./TrendingPlots/%s.png", name) );
  can->SaveAs( Form("%s.png", name) );

  delete can;

}
