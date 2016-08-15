class LMEECutLib {
  
public:
  
  static enum enCentSel {
    kPbPb2011Central=0,
    kPbPb2011_00to10,
    kPbPb2011MidCentral,
    kPbPb2011_10to20,
    kPbPb2011SemiCentral,
    kPbPb2011_20to50,
    kPbPb2011_10to50,
    kPbPb2011_00to50,
    kCENTSELMAX
  };
  
  static enum enPairCut {
    kPairCut_OFF=0,
    kPairCut_theta20,
    kPairCut_theta50,
    kPairCut_theta100,
    kPairCut_mee10_theta30,
    kPairCut_mee20_theta20, // maybe for prefilter with SPD+SDD
    kPairCut_mee20_theta50,
    kPairCut_mee30_theta60,
    kPairCut_mee40_theta80,
    kPairCut_mee60_theta100,
    kPairCut_mee200_theta300, // for testing
    kPairCut_phiv157_mee40,
    kPairCut_phiv157_mee60,
    kPairCut_phiv157_mee80,
    kPairCut_phiv157_mee100,
    kPairCut_phiv236_mee40,
    kPairCut_phiv236_mee60,
    kPairCut_phiv236_mee80,
    kPairCut_phiv236_mee100
  };
  
  static enum enKineCut {
    kKineCut_p50inf_eta150=0, // for resolution extraction
    kKineCut_pt50_eta080,
    kKineCut_pt50_eta090,
    kKineCut_pt200_eta080,
    kKineCut_pt200_eta090,
    kKineCut_pt300_eta080,
    kKineCut_pt300_eta090,
    kKineCut_switch_PIDcorr,// for all cutSets with global (ana) track pt cuts below 400 MeV, the PID correction made with pt>200 MeV is used.
    kKineCut_pt400_eta080=kKineCut_switch_PIDcorr,
    kKineCut_pt400_eta090
  };
  
  static enum LMMECutSet {
    kCut01=1,
    kCut02,
    kCut03,
    kCut04,
    kCut05,
    kCut06,
    kCut07,
    kCut08,
    kCut09,
    kCut10,
    kCut11,
    kCut12,
    kCut13,
    kCut14,
    kCut15,
    kCut16,
    kCut17,
    kCut18,
    kCut19,
    kCut20,
    kPbPb2011_pidITSTPCTOFif_trkSPDorSDD_7_V0excl,          // syst 8
    kPbPb2011_pidITSTPCTOFif_trkSPDfirst_7_V0excl,          // syst 7
    kPbPb2011MC_pi0Dal_1,
    kPbPb2011_pidITS2gevTPCTOFif_trkSPD5orSDD4cls_6_tight,  // syst 6 (no train run yet)
    kPbPb2011_pidITS2gevTPCTOFif_trkSPDfirst5cls_6_tight,   // syst 5 (no train run yet)
    kPbPb2011_pidITS2gevTPCTOFif_trkSPDorSDD_5_tight,       // syst 4
    kPbPb2011_pidITS2gevTPCTOFif_trkSPDfirst_5_tight,       // syst 2
    kPbPb2011_pidITSTPCTOFif_trkSPD5orSDD4cls_4,            // syst 3
    kPbPb2011_pidITSTPCTOFif_trkSPDfirst5cls_4,             // syst 1
    kPbPb2011_pidITSTPC_trkSPDfirst_3,            // (cutSet w/o pairing)
    kPbPb2011_pidTPC_trkSPDfirst_3,               // (cutSet w/o pairing)
    kPbPb2011_pidITSTPCTOFif_trkSPDorSDD_2_loose, // (cutSet w/o pairing)
    kPbPb2011_pidITSTPCTOFif_trkSPDfirst_2_loose, // (cutSet w/o pairing)
    kPbPb2011_pidTPCTOF_trkSPDorSDD_2_loose,      // (cutSet w/o pairing)
    kPbPb2011_pidTPCTOF_trkSPDfirst_2_loose,      // (cutSet w/o pairing)
    kPbPb2011_pidITSTPCTOFif_trkSPDorSDD_1,
    kPbPb2011_pidITSTPCTOFif_trkSPDfirst_1,       // cutSet for Technical Preliminaries for QM2014 (no prefilter used!)
    kPbPb2011_pidTPCTOF_trkSPDorSDD_1,
    kPbPb2011_pidTPCTOF_trkSPDfirst_1,
    kCutSetMAX,
    //
    // PID
    kPbPb2011PID_ITSTPCTOFif_1=kCutSetMAX,  // ITS+TPC+TOFf
    kPbPb2011PID_ITSTPCTOFif_2,
    kPbPb2011PID_ITSTPCTOFif_3,
    kPbPb2011PID_ITSTPCTOFif_LOOSE,
    kPbPb2011PID_ITSTPCTOFif_looseTPC,
    kPbPb2011PID_TPCITS_1,                  // ITS+TPC
    kPbPb2011PID_TPCITSif_2,
    kPbPb2011PID_TPCITS_3,
    kPbPb2011PID_TPCTOF_1,                  // TPC+TOF
    kPbPb2011PID_TPCTOF_LOOSE,
    kPbPb2011PID_TPC_1,                     // TPC
    kPbPb2011PID_TPC_2,
    kPbPb2011PID_TPC_pre,
    kPbPb2011PID_ITS_1,                     // ITS
    kPbPb2011PID_V0_1,                      // V0
    kPbPb2011PID_V0_1_pionCont,             // additional pion cut in TPC to check pion contamination in ITS dEdx
    kPbPb2011PID_V0_2_TOFif,
    kPbPb2011PID_V0_2_TOFif_pionCont,
    kPbPb2011PID_OFF,
    kPIDMAX,
    //
    // Quality
    kPbPb2011TRK_SPDfirst_1=kPIDMAX,    // main track selection. with SPD first
    kPbPb2011TRK_SPDfirst_2,            // main track selection. >5 ITS cls
    kPbPb2011TRK_SDDfirstSPDnone_1,     // complimentary tracks. strictly without SPD, to be combined with others!
    kPbPb2011TRK_SDDfirstSPDnone_2,     // complimentary tracks. >4 ITS cls
    kPbPb2011TRK_SPDorSDD_1,
    kPbPb2011TRK_SPDorSDD_2,
    //                                  // SPDfirst_1 (12000 tracks in 66 events [pt>50 MeV, no PID cuts])
    //                                  // SPDorSDD_1 (16600 tracks ...)
    kPbPb2011TRK_FilterBit0,            // TPC only + 3 ITS clu (almost flat in phi, ITS PID avail for 28600 tr.). (30900 tracks ...)
    kPbPb2011TRK_FilterBit1,            // ITS standalone, no TPC PID avail. (9100 tracks ...)
    kPbPb2011TRK_FilterBit2,            // TPC + SPDany (holes in phi). (21700 tracks ...)
    kPbPb2011TRK_FilterBit4,
    kPbPb2011TRK_FilterBit6,
    kPbPb2011TRK_V0_1,
    kPbPb2011TRK_V0_2_loose,
    kPbPb2011TRK_V0_3_Arm,
    kQualityMAX
  };
  
  //char* LMEECutNames[kQualityMAX] = { "PbPb2011TPCandTOF","PbPb2011TPCorTOF"};
  
  LMEECutLib();
  
  AliDielectronEventCuts*     GetEventCuts(Int_t cutSet, Bool_t hasMC=kFALSE);
  AliAnalysisCuts*            GetCentralityCuts(Bool_t hasMC=kFALSE);
  AliDielectronTrackRotator*  GetTrackRotator();
  AliDielectronMixingHandler* GetMixingHandler();
  
  AliAnalysisCuts* GetPairCutsAna();
  AliAnalysisCuts* GetPairCutsPre(Int_t cutSet=-1);
  // main track cut functions:
  AliAnalysisCuts* GetTrackCutsAna();
  AliAnalysisCuts* GetTrackCutsPre();
  AliAnalysisCuts* GetESDTrackCutsAna(Int_t cutSet);
  AliAnalysisCuts* GetKineCutsAna(); // needed public for efficiency task
  
  void      SetIsQATask(Bool_t b=kTRUE)         { fIsQATask=b; }
  void      SetIsRandomRejTask(Bool_t b=kTRUE)  { fIsRandomRejTask=b; }
  void      SetDoRejectionStep(Bool_t b=kTRUE)  { fDoRejectionStep=b; }
  void      SetFillPureMC(Bool_t b=kTRUE)       { fFillPureMC=b; }
  void      SetIsESDTask(Bool_t b=kTRUE)        { if (b) fIsESDTask=kYes; else fIsESDTask=kNo; }

  Bool_t    GetDoRejectionStep() { return fDoRejectionStep; }
  
  void      SetITSSigmaEleCorrection(AliDielectron *die, Int_t corrZdim, Int_t corrYdim); //giving default value fails: /* = AliDielectronVarManager::kEta*/
  void      SetTPCSigmaEleCorrection(AliDielectron *die, Int_t corrZdim, Int_t corrYdim);
  void      SetITSSigmaEleCorrectionMC(TNamed* task, Int_t corrZdim, Int_t corrYdim);
  void      SetTPCSigmaEleCorrectionMC(TNamed* task, Int_t corrZdim, Int_t corrYdim);
  
  void      AddMCSignals(TNamed* task, Int_t cutDefinition);
  void      InitHistograms(AliDielectron *die, Int_t cutDefinition);
  void      InitCF(AliDielectron* die, Int_t cutDefinition);
  
  // kept public to avoid writing setters...
  Int_t     selectedCentrality;
  Int_t     selectedPIDAna;
  Int_t     selectedPIDPre;
  Int_t     selectedQualityAna;
  Int_t     selectedQualityPre;
  Int_t     selectedKineCutsAna;
  Int_t     selectedKineCutsPre;
  Int_t     selectedPairCutsAna;
  Int_t     selectedPairCutsPre;
  
  
private:
  static enum enCutType {
    kInclude = 0,
    kExclude = 1
  };
  
  static enum enHistVars {
    kMee=0, kMee500,
    kPtee, kP2D, kRuns,
    kPhiV, kOpAng, kOpAng2, kOpAng3,
    kEta2D, kEta3D, kPhi2D, kY3D,
    kSigmaEle, kSigmaOther, kTPCdEdx,
    kPairDCAsigXY, kPairDCAabsXY, kPairLinDCAsigXY, kPairLinDCAabsXY
  };
  
  static enum enESDTask {
    kNo    = 0,
    kYes   = 1,
    kUnset = 2
  };
  
  // internal track cut functions (called by GetTrackCuts):
  //AliAnalysisCuts* GetKineCutsAna(); // needed public for efficiency task
  AliAnalysisCuts* GetKineCutsPre(Int_t cutSet=-1);
  AliAnalysisCuts* GetPIDCuts(Int_t cutSet=-1);
  AliAnalysisCuts* GetQualityCuts(Int_t cutSet=-1, Int_t doExclusion=kInclude);
  AliAnalysisCuts* GetMCTrackCuts();
  // helper functions
  TVectorD* BinsToVector(Int_t nbins, Double_t min, Double_t max);
  TVectorD* GetVector(Int_t var);
  
  Bool_t    fIsQATask;
  Bool_t    fIsRandomRejTask;
  Bool_t    fDoRejectionStep;
  Bool_t    fFillPureMC;
  Int_t     fIsESDTask;
  
};


//_______________________________________________________________________________________________
LMEECutLib::LMEECutLib() :
selectedCentrality(-1),
selectedPIDAna(-1),
selectedPIDPre(-1),
selectedQualityAna(-1),
selectedQualityPre(-1),
selectedKineCutsAna(LMEECutLib::kKineCut_pt200_eta080),
selectedKineCutsPre(LMEECutLib::kKineCut_pt50_eta090),
selectedPairCutsAna(LMEECutLib::kPairCut_OFF),
selectedPairCutsPre(LMEECutLib::kPairCut_OFF),
fIsQATask(kFALSE),
fIsRandomRejTask(kFALSE),
fDoRejectionStep(kFALSE),
fFillPureMC(kFALSE),
fIsESDTask(kUnset)
{
  // Constructor
}


//_______________________________________________________________________________________________
void LMEECutLib::SetITSSigmaEleCorrection(AliDielectron *die, Int_t corrZdim, Int_t corrYdim) {
  //
  // eta correction for the centroid and width of electron sigmas in the ITS, can be one/two/three-dimensional
  //
  printf("starting LMEECutLib::SetITSSigmaEleCorrection()\n --> for MC task, use ::SetITSSigmaEleCorrectionMC()!\n");
  printf(" correction Zdim = %s\n", AliDielectronVarManager::GetValueName(corrZdim));
  printf(" correction Ydim = %s\n", AliDielectronVarManager::GetValueName(corrYdim));
  
  TF2 *fCntrdCorr=0x0;
  TF2 *fWdthCorr=0x0;
  Double_t fitMinZdim=   0;
  Double_t fitMaxZdim=4000;
  Double_t fitMinEta=-0.9, fitMaxEta=+0.9; // Ydim, usually eta
  
  if (corrYdim!=AliDielectronVarManager::kEta) { printf(" correction only available for Ydim = eta!\n"); return; }
  
  if (corrZdim==AliDielectronVarManager::kRefMultTPConly)
  {
    const int netabins=18; // typically the number of eta bins incl over & underflow.
    TF1* fcnMean_px[netabins];
    TF1* fcnWidth_px[netabins];
    for (int i=0; i<netabins; i++) {
      fcnMean_px[i]  = new TF1(Form("fcnMean_px_%i" , i), "[0]+[1]*x", fitMinZdim, fitMaxZdim);
      fcnWidth_px[i] = new TF1(Form("fcnWidth_px_%i", i), "[0]+[1]*x", fitMinZdim, fitMaxZdim);
    }
    if (selectedKineCutsAna<kKineCut_switch_PIDcorr) { // use PID correction for pt>200 MeV/c
      fcnMean_px[0]->SetParameters(-2.317592e-01, 2.015334e-05);
      fcnMean_px[1]->SetParameters(-2.317592e-01, 2.015334e-05);
      fcnMean_px[2]->SetParameters(-2.388820e-01, 2.100415e-05);
      fcnMean_px[3]->SetParameters(-2.368509e-01, 2.649807e-05);
      fcnMean_px[4]->SetParameters(-2.335289e-01, 2.771517e-05);
      fcnMean_px[5]->SetParameters(-2.203853e-01, 2.299109e-05);
      fcnMean_px[6]->SetParameters(-2.231970e-01, 2.695141e-05);
      fcnMean_px[7]->SetParameters(-2.019306e-01, 2.349839e-05);
      fcnMean_px[8]->SetParameters(-1.999166e-01, 2.257996e-05);
      fcnMean_px[9]->SetParameters(-1.980436e-01, 1.807920e-05);
      fcnMean_px[10]->SetParameters(-1.823919e-01, 9.162104e-06);
      fcnMean_px[11]->SetParameters(-2.069340e-01, 1.779933e-05);
      fcnMean_px[12]->SetParameters(-1.953123e-01, 5.935998e-06);
      fcnMean_px[13]->SetParameters(-2.010407e-01, 1.546144e-05);
      fcnMean_px[14]->SetParameters(-1.978739e-01, 8.981407e-06);
      fcnMean_px[15]->SetParameters(-2.224901e-01, 1.851342e-05);
      fcnMean_px[16]->SetParameters(-2.160992e-01, 1.541774e-05);
      fcnMean_px[17]->SetParameters(-2.160992e-01, 1.541774e-05);
      //
      fcnWidth_px[0]->SetParameters(1.042039e+00, 7.178529e-05);
      fcnWidth_px[1]->SetParameters(1.042039e+00, 7.178529e-05);
      fcnWidth_px[2]->SetParameters(1.042571e+00, 6.053239e-05);
      fcnWidth_px[3]->SetParameters(1.055642e+00, 5.049183e-05);
      fcnWidth_px[4]->SetParameters(1.062871e+00, 4.467255e-05);
      fcnWidth_px[5]->SetParameters(1.077608e+00, 3.696654e-05);
      fcnWidth_px[6]->SetParameters(1.071887e+00, 4.272863e-05);
      fcnWidth_px[7]->SetParameters(1.079705e+00, 4.932264e-05);
      fcnWidth_px[8]->SetParameters(1.091248e+00, 4.373216e-05);
      fcnWidth_px[9]->SetParameters(1.089759e+00, 4.158001e-05);
      fcnWidth_px[10]->SetParameters(1.094837e+00, 3.824771e-05);
      fcnWidth_px[11]->SetParameters(1.075169e+00, 4.479051e-05);
      fcnWidth_px[12]->SetParameters(1.083042e+00, 3.299472e-05);
      fcnWidth_px[13]->SetParameters(1.067151e+00, 4.178259e-05);
      fcnWidth_px[14]->SetParameters(1.052284e+00, 4.735530e-05);
      fcnWidth_px[15]->SetParameters(1.030472e+00, 6.617467e-05);
      fcnWidth_px[16]->SetParameters(1.043461e+00, 7.115885e-05);
      fcnWidth_px[17]->SetParameters(1.043461e+00, 7.115885e-05);
    }
    else {  // use PID correction for pt>400 MeV/c
      fcnMean_px[0]->SetParameters(-3.529055e-01, 1.327113e-05);
      fcnMean_px[1]->SetParameters(-3.529055e-01, 1.327113e-05);
      fcnMean_px[2]->SetParameters(-3.592750e-01, 2.087731e-05);
      fcnMean_px[3]->SetParameters(-3.969130e-01, 5.324868e-05);
      fcnMean_px[4]->SetParameters(-3.743942e-01, 3.691861e-05);
      fcnMean_px[5]->SetParameters(-3.694300e-01, 3.904927e-05);
      fcnMean_px[6]->SetParameters(-3.607261e-01, 1.917498e-05);
      fcnMean_px[7]->SetParameters(-3.600538e-01, 2.700674e-05);
      fcnMean_px[8]->SetParameters(-3.536578e-01, 2.846538e-05);
      fcnMean_px[9]->SetParameters(-3.552292e-01, 2.991948e-05);
      fcnMean_px[10]->SetParameters(-3.186401e-01, 6.020452e-06);
      fcnMean_px[11]->SetParameters(-3.146004e-01, 1.215062e-06);
      fcnMean_px[12]->SetParameters(-3.221104e-01, 5.948753e-06);
      fcnMean_px[13]->SetParameters(-3.221268e-01, 1.084347e-05);
      fcnMean_px[14]->SetParameters(-2.931151e-01, -1.277723e-05);
      fcnMean_px[15]->SetParameters(-3.357687e-01, 1.147554e-05);
      fcnMean_px[16]->SetParameters(-3.835154e-01, 4.926267e-05);
      fcnMean_px[17]->SetParameters(-3.835154e-01, 4.926267e-05);
      //
      fcnWidth_px[0]->SetParameters(1.065834e+00, 4.059001e-05);
      fcnWidth_px[1]->SetParameters(1.065834e+00, 4.059001e-05);
      fcnWidth_px[2]->SetParameters(1.047882e+00, 4.991033e-05);
      fcnWidth_px[3]->SetParameters(1.045476e+00, 5.880931e-05);
      fcnWidth_px[4]->SetParameters(1.062518e+00, 4.287982e-05);
      fcnWidth_px[5]->SetParameters(1.062254e+00, 4.977838e-05);
      fcnWidth_px[6]->SetParameters(1.059952e+00, 4.242825e-05);
      fcnWidth_px[7]->SetParameters(1.059177e+00, 5.430946e-05);
      fcnWidth_px[8]->SetParameters(1.080617e+00, 4.751661e-05);
      fcnWidth_px[9]->SetParameters(1.067084e+00, 5.367871e-05);
      fcnWidth_px[10]->SetParameters(1.097499e+00, 3.357789e-05);
      fcnWidth_px[11]->SetParameters(1.084897e+00, 3.741411e-05);
      fcnWidth_px[12]->SetParameters(1.080744e+00, 3.506669e-05);
      fcnWidth_px[13]->SetParameters(1.079388e+00, 3.235666e-05);
      fcnWidth_px[14]->SetParameters(1.081288e+00, 2.310908e-05);
      fcnWidth_px[15]->SetParameters(1.043360e+00, 5.170997e-05);
      fcnWidth_px[16]->SetParameters(1.010796e+00, 8.893934e-05);
      fcnWidth_px[17]->SetParameters(1.010796e+00, 8.893934e-05);
    }
    Bool_t isFcnZdim=kTRUE; // true if 1D functions go along X direction of 2D histogram.
    Int_t  extrapolate=1;   // fill underflow and overflow bins.
    TH2F *hMean_2Dfit  = new TH2F("hMean_2Dfit" ,"hMean_2Dfit" ,19,400,2300,16,-0.8,0.8);
    TH2F *hWidth_2Dfit = new TH2F("hWidth_2Dfit","hWidth_2Dfit",19,400,2300,16,-0.8,0.8);
    for (int ix=1-extrapolate; ix<=hMean_2Dfit->GetNbinsX()+extrapolate; ix++) {
      for (int iy=1-extrapolate; iy<=hMean_2Dfit->GetNbinsY()+extrapolate; iy++) {
        if (isFcnZdim) {
          //cout << " hMean_2Dfit->GetBinCenter(ix) = " << hMean_2Dfit->GetXaxis()->GetBinCenter(ix) << endl;
          hMean_2Dfit->SetBinContent(ix, iy, fcnMean_px[iy]->Eval(hMean_2Dfit->GetXaxis()->GetBinCenter(ix)));
          hWidth_2Dfit->SetBinContent(ix, iy, fcnWidth_px[iy]->Eval(hWidth_2Dfit->GetXaxis()->GetBinCenter(ix)));
        } else {
          //cout << " hMean_2Dfit->GetBinCenter(iy) = " << hMean_2Dfit->GetYaxis()->GetBinCenter(iy) << endl;
          hMean_2Dfit->SetBinContent(ix, iy, fcnMean_px[ix]->Eval(hMean_2Dfit->GetYaxis()->GetBinCenter(iy)));
          hWidth_2Dfit->SetBinContent(ix, iy, fcnWidth_px[ix]->Eval(hWidth_2Dfit->GetYaxis()->GetBinCenter(iy)));
        }
      }
    }
  } // end: if (kRefMultTPConly)
  else {
    printf(" no correction available for Zdim = %s!\n", AliDielectronVarManager::GetValueName(corrZdim));
    printf(" no correction applied!\n");
    return;
  }
  
//  TCanvas* c1 = new TCanvas("c1","c1",800,750);
//  c1->SetRightMargin(0.22);
//  hMean_2Dfit->DrawCopy("colz");
//  c1->Print(Form("cITScorr_mean_kinecut%i.pdf",selectedKineCutsAna));
//  hWidth_2Dfit->DrawCopy("colz");
//  c1->Print(Form("cITScorr_width_kinecut%i.pdf",selectedKineCutsAna));
  
  die->SetCentroidCorrFunctionITS(hMean_2Dfit, corrZdim, corrYdim);
  die->SetWidthCorrFunctionITS(hWidth_2Dfit, corrZdim, corrYdim);
  printf(" ITS PID eta correction loaded!\n");
}

//_______________________________________________________________________________________________
void LMEECutLib::SetTPCSigmaEleCorrection(AliDielectron *die, Int_t corrZdim, Int_t corrYdim) {
  //
  // eta correction for the centroid and width of electron sigmas in the TPC, can be one/two/three-dimensional
  //
  printf("starting LMEECutLib::SetTPCSigmaEleCorrection()\n --> for MC task, use ::SetTPCSigmaEleCorrectionMC()!\n");
  printf(" correction Zdim = %s\n", AliDielectronVarManager::GetValueName(corrZdim));
  printf(" correction Ydim = %s\n", AliDielectronVarManager::GetValueName(corrYdim));
  //Bool_t hasMC=die->GetHasMC();
  //Bool_t hasTuneOnData=kFALSE; //((AliAnalysisTaskPIDResponse*)AliAnalysisManager::GetAnalysisManager()->GetTasks()->At(0))->GetTuneOnData();
  //printf("tune on data switched: %d \n",hasTuneOnData);
  // printf("name task at 0: %s \n",AliAnalysisManager::GetAnalysisManager()->GetTasks()->At(0)->GetName());
  // printf("input event %p \n", AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  // printf("pid response %p \n",((AliInputEventHandler*)AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler())->GetPIDResponse());
  // printf("pid response task %p \n",AliAnalysisManager::GetAnalysisManager()->GetTasks()->At(0));
  // AliAnalysisManager::GetAnalysisManager()->GetTasks()->At(0)->Dump();;
  
  // AliAnalysisManager* man = AliAnalysisManager::GetAnalysisManager();
  // AliInputEventHandler* inputHandler = dynamic_cast<AliInputEventHandler*>(man->GetInputEventHandler());
  // AliPIDResponse* pidResponse = inputHandler->GetPIDResponse();
  // if(pidResponse) hasTuneOnData = pidResponse->IsTunedOnData();
  // printf("man %p inp %p pid %p ====> %d \n",man,inputHandler,pidResponse,hasTuneOnData);
  
  //TF2 *fCntrdCorr=0x0;
  //TF2 *fWdthCorr=0x0;
  Double_t fitMinZdim=   0;
  Double_t fitMaxZdim=4000;
  //Double_t fitMinEta=-0.9, fitMaxEta=+0.9; // Ydim, usually eta
  
  if (corrYdim!=AliDielectronVarManager::kEta) { printf(" correction only available for Ydim = eta!\n"); return; }
  
  if (corrZdim==AliDielectronVarManager::kRefMultTPConly)
  {
    const int netabins=18; // typically the number of eta bins incl over & underflow.
    TF1* fcnMean_px[netabins];
    TF1* fcnWidth_px[netabins];
    for (int i=0; i<netabins; i++) {
      fcnMean_px[i]  = new TF1(Form("fcnMean_px_%i" , i), "[0]+[1]*x", fitMinZdim, fitMaxZdim);
      fcnWidth_px[i] = new TF1(Form("fcnWidth_px_%i", i), "[0]+[1]*x", fitMinZdim, fitMaxZdim);
    }
    if (selectedKineCutsAna<kKineCut_switch_PIDcorr) { // use PID correction for pt>200 MeV/c
      fcnMean_px[0]->SetParameters(-1.156716e-02, -1.301817e-04);
      fcnMean_px[1]->SetParameters(-1.156716e-02, -1.301817e-04);
      fcnMean_px[2]->SetParameters(-1.312152e-01, -1.587708e-04);
      fcnMean_px[3]->SetParameters(-1.907976e-01, -1.956057e-04);
      fcnMean_px[4]->SetParameters(-1.797112e-01, -2.172959e-04);
      fcnMean_px[5]->SetParameters(-1.161102e-01, -2.236802e-04);
      fcnMean_px[6]->SetParameters(-1.083095e-02, -2.269384e-04);
      fcnMean_px[7]->SetParameters(8.482330e-02, -2.245246e-04);
      fcnMean_px[8]->SetParameters(1.047848e-01, -2.000582e-04);
      fcnMean_px[9]->SetParameters(1.647767e-01, -1.951507e-04);
      fcnMean_px[10]->SetParameters(1.140301e-01, -1.948394e-04);
      fcnMean_px[11]->SetParameters(6.306847e-03, -1.913186e-04);
      fcnMean_px[12]->SetParameters(-1.084654e-01, -1.771152e-04);
      fcnMean_px[13]->SetParameters(-1.757320e-01, -1.619230e-04);
      fcnMean_px[14]->SetParameters(-1.746252e-01, -1.296587e-04);
      fcnMean_px[15]->SetParameters(-1.006578e-01, -8.307905e-05);
      fcnMean_px[16]->SetParameters(1.223413e-02, -5.322799e-05);
      fcnMean_px[17]->SetParameters(1.223413e-02, -5.322799e-05);
      //
      fcnWidth_px[0]->SetParameters(1.086622e+00, 6.606505e-05);
      fcnWidth_px[1]->SetParameters(1.086622e+00, 6.606505e-05);
      fcnWidth_px[2]->SetParameters(1.082054e+00, 5.711400e-05);
      fcnWidth_px[3]->SetParameters(1.077950e+00, 5.009987e-05);
      fcnWidth_px[4]->SetParameters(1.088628e+00, 4.443494e-05);
      fcnWidth_px[5]->SetParameters(1.111643e+00, 4.061893e-05);
      fcnWidth_px[6]->SetParameters(1.146064e+00, 3.688119e-05);
      fcnWidth_px[7]->SetParameters(1.169822e+00, 3.917164e-05);
      fcnWidth_px[8]->SetParameters(1.200732e+00, 3.658010e-05);
      fcnWidth_px[9]->SetParameters(1.213068e+00, 3.089858e-05);
      fcnWidth_px[10]->SetParameters(1.185535e+00, 3.263566e-05);
      fcnWidth_px[11]->SetParameters(1.155687e+00, 2.876342e-05);
      fcnWidth_px[12]->SetParameters(1.118644e+00, 3.407548e-05);
      fcnWidth_px[13]->SetParameters(1.089084e+00, 3.803547e-05);
      fcnWidth_px[14]->SetParameters(1.076255e+00, 4.455558e-05);
      fcnWidth_px[15]->SetParameters(1.074615e+00, 5.149723e-05);
      fcnWidth_px[16]->SetParameters(1.080458e+00, 5.354965e-05);
      fcnWidth_px[17]->SetParameters(1.080458e+00, 5.354965e-05);
    }
    else {  // use PID correction for pt>400 MeV/c
      fcnMean_px[0]->SetParameters(1.069692e-01, -1.853931e-04);
      fcnMean_px[1]->SetParameters(1.069692e-01, -1.853931e-04);
      fcnMean_px[2]->SetParameters(-4.447889e-02, -2.098573e-04);
      fcnMean_px[3]->SetParameters(-8.858894e-02, -2.446114e-04);
      fcnMean_px[4]->SetParameters(-1.794249e-02, -2.697981e-04);
      fcnMean_px[5]->SetParameters(1.061417e-01, -2.766702e-04);
      fcnMean_px[6]->SetParameters(2.621812e-01, -2.810459e-04);
      fcnMean_px[7]->SetParameters(4.075723e-01, -2.616764e-04);
      fcnMean_px[8]->SetParameters(4.933002e-01, -2.473390e-04);
      fcnMean_px[9]->SetParameters(5.797911e-01, -2.468768e-04);
      fcnMean_px[10]->SetParameters(4.790437e-01, -2.384890e-04);
      fcnMean_px[11]->SetParameters(3.266869e-01, -2.333392e-04);
      fcnMean_px[12]->SetParameters(1.541458e-01, -2.282300e-04);
      fcnMean_px[13]->SetParameters(1.676668e-02, -2.212672e-04);
      fcnMean_px[14]->SetParameters(-2.689889e-02, -2.052884e-04);
      fcnMean_px[15]->SetParameters(1.024606e-02, -1.558609e-04);
      fcnMean_px[16]->SetParameters(1.380813e-01, -1.459855e-04);
      fcnMean_px[17]->SetParameters(1.380813e-01, -1.459855e-04);
      //
      fcnWidth_px[0]->SetParameters(1.100951e+00, 5.283585e-05);
      fcnWidth_px[1]->SetParameters(1.100951e+00, 5.283585e-05);
      fcnWidth_px[2]->SetParameters(1.069972e+00, 5.493703e-05);
      fcnWidth_px[3]->SetParameters(1.084033e+00, 4.183144e-05);
      fcnWidth_px[4]->SetParameters(1.097820e+00, 4.153080e-05);
      fcnWidth_px[5]->SetParameters(1.115891e+00, 4.219773e-05);
      fcnWidth_px[6]->SetParameters(1.138289e+00, 4.310334e-05);
      fcnWidth_px[7]->SetParameters(1.162934e+00, 3.682753e-05);
      fcnWidth_px[8]->SetParameters(1.198148e+00, 2.239657e-05);
      fcnWidth_px[9]->SetParameters(1.174714e+00, 3.120157e-05);
      fcnWidth_px[10]->SetParameters(1.153562e+00, 3.546540e-05);
      fcnWidth_px[11]->SetParameters(1.136738e+00, 3.161306e-05);
      fcnWidth_px[12]->SetParameters(1.112814e+00, 3.477281e-05);
      fcnWidth_px[13]->SetParameters(1.089942e+00, 3.794432e-05);
      fcnWidth_px[14]->SetParameters(1.075360e+00, 4.227114e-05);
      fcnWidth_px[15]->SetParameters(1.073305e+00, 4.787609e-05);
      fcnWidth_px[16]->SetParameters(1.108154e+00, 4.165704e-05);
      fcnWidth_px[17]->SetParameters(1.108154e+00, 4.165704e-05);      
    }
    Bool_t isFcnZdim=kTRUE; // true if 1D functions go along X direction of 2D histogram.
    Int_t  extrapolate=1;   // fill underflow and overflow bins.
    TH2F *hMean_2Dfit  = new TH2F("hMean_2Dfit" ,"hMean_2Dfit" ,19,400,2300,16,-0.8,0.8);
    TH2F *hWidth_2Dfit = new TH2F("hWidth_2Dfit","hWidth_2Dfit",19,400,2300,16,-0.8,0.8);
    for (int ix=1-extrapolate; ix<=hMean_2Dfit->GetNbinsX()+extrapolate; ix++) {
      for (int iy=1-extrapolate; iy<=hMean_2Dfit->GetNbinsY()+extrapolate; iy++) {
        if (isFcnZdim) {
          //cout << " hMean_2Dfit->GetBinCenter(ix) = " << hMean_2Dfit->GetXaxis()->GetBinCenter(ix) << endl;
          hMean_2Dfit->SetBinContent(ix, iy, fcnMean_px[iy]->Eval(hMean_2Dfit->GetXaxis()->GetBinCenter(ix)));
          hWidth_2Dfit->SetBinContent(ix, iy, fcnWidth_px[iy]->Eval(hWidth_2Dfit->GetXaxis()->GetBinCenter(ix)));
        } else {
          //cout << " hMean_2Dfit->GetBinCenter(iy) = " << hMean_2Dfit->GetYaxis()->GetBinCenter(iy) << endl;
          hMean_2Dfit->SetBinContent(ix, iy, fcnMean_px[ix]->Eval(hMean_2Dfit->GetYaxis()->GetBinCenter(iy)));
          hWidth_2Dfit->SetBinContent(ix, iy, fcnWidth_px[ix]->Eval(hWidth_2Dfit->GetYaxis()->GetBinCenter(iy)));
        }
      }
    }
  } // end: if (kRefMultTPConly)
  else {
    printf(" no correction available for Zdim = %s!\n", AliDielectronVarManager::GetValueName(corrZdim));
    printf(" no correction applied!\n");
    return;
  }
  
//  TCanvas* c1 = new TCanvas("c1","c1",800,750);
//  c1->SetRightMargin(0.22);
//  hMean_2Dfit->DrawCopy("colz");
//  c1->Print(Form("cTPCcorr_mean_kinecut%i.pdf",selectedKineCutsAna));
//  hWidth_2Dfit->DrawCopy("colz");
//  c1->Print(Form("cTPCcorr_width_kinecut%i.pdf",selectedKineCutsAna));
  
  die->SetCentroidCorrFunction(hMean_2Dfit, corrZdim, corrYdim);
  die->SetWidthCorrFunction(hWidth_2Dfit, corrZdim, corrYdim);
  printf(" TPC PID eta correction loaded!\n");
}


//_______________________________________________________________________________________________
void LMEECutLib::SetITSSigmaEleCorrectionMC(TNamed* task, Int_t corrZdim, Int_t corrYdim) {
  //
  // MC post-correction for the centroid and width of electron sigmas in the ITS, can be one/two/three-dimensional
  //
  printf("starting LMEECutLib::SetITSSigmaEleCorrectionMC()\n");
  printf(" correction Zdim = %s\n", AliDielectronVarManager::GetValueName(corrZdim));
  printf(" correction Ydim = %s\n", AliDielectronVarManager::GetValueName(corrYdim));
  Bool_t use1Dfunctions;  // true if 1D functions are used instead of a fully defined 2D map.
  Bool_t isFcnZdim;       // true if 1D functions go along X direction of 2D histogram.
  Int_t  extrapolate=0;   // =1 to fill underflow and overflow bins.
  TH2F *hMean_2Dfit=0x0;
  TH2F *hWidth_2Dfit=0x0;
  TF1* fcnMean_px[20];
  TF1* fcnWidth_px[20];
  Double_t fitMinEta =-0.8, fitMaxEta =+0.8; // Ydim, usually eta
  Double_t fitMinZdim, fitMaxZdim;
  
  if (corrYdim!=AliDielectronVarManager::kEta) { printf(" correction only available for Ydim = eta!\n"); return; }
  
  // in MC (LHC14a1b) the P-dependence is much stronger than the Nacc-dependence, so we correct for P...
  if (corrZdim==AliDielectronVarManager::kP)
  {
    use1Dfunctions=kTRUE;  isFcnZdim=kTRUE;  extrapolate=1;
    fitMinEta=-0.8, fitMaxEta=+0.8;
    fitMinZdim=0  , fitMaxZdim=5.; // GeV/c
    const int netabins=8;
    const int nZbins=25;
    hMean_2Dfit  = new TH2F("hMean_2Dfit" ,"hMean_2Dfit" ,nZbins,fitMinZdim,fitMaxZdim,netabins,fitMinEta,fitMaxEta);
    hWidth_2Dfit = new TH2F("hWidth_2Dfit","hWidth_2Dfit",nZbins,fitMinZdim,fitMaxZdim,netabins,fitMinEta,fitMaxEta);
    
    for (int i=1-extrapolate; i<=netabins+extrapolate; i++) {
      fcnMean_px[i]  = new TF1(Form("fcnMean_px_%i" , i), "[0]+[1]*x", fitMinZdim, fitMaxZdim);
      fcnWidth_px[i] = new TF1(Form("fcnWidth_px_%i", i), "[0]+[1]*x", fitMinZdim, fitMaxZdim);
    }
    // corrections are done on the 10-50% centrality MC (LHC14a1b).
    // entries [0] and [max] have to stay empty if extrapolate=0.
    fcnMean_px[0]->SetParameters(2.392202e-02, -1.973464e-01);
    fcnMean_px[1]->SetParameters(2.392202e-02, -1.973464e-01);
    fcnMean_px[2]->SetParameters(-7.582297e-02, -1.854471e-01);
    fcnMean_px[3]->SetParameters(-1.338714e-01, -2.262040e-01);
    fcnMean_px[4]->SetParameters(-1.568579e-01, -2.689306e-01);
    fcnMean_px[5]->SetParameters(-1.644852e-01, -2.093024e-01);
    fcnMean_px[6]->SetParameters(-1.260372e-01, -2.248881e-01);
    fcnMean_px[7]->SetParameters(-7.303062e-02, -2.098463e-01);
    fcnMean_px[8]->SetParameters(3.796873e-03, -1.959614e-01);
    fcnMean_px[9]->SetParameters(3.796873e-03, -1.959614e-01);
    //
    fcnWidth_px[0]->SetParameters(1.062272e+00, -1.527472e-02);
    fcnWidth_px[1]->SetParameters(1.062272e+00, -1.527472e-02);
    fcnWidth_px[2]->SetParameters(1.083509e+00, -1.745190e-02);
    fcnWidth_px[3]->SetParameters(1.099541e+00, -2.316445e-02);
    fcnWidth_px[4]->SetParameters(1.123306e+00, -4.434398e-02);
    fcnWidth_px[5]->SetParameters(1.111333e+00, -2.527514e-02);
    fcnWidth_px[6]->SetParameters(1.101949e+00, -2.525758e-02);
    fcnWidth_px[7]->SetParameters(1.089704e+00, -2.994100e-02);
    fcnWidth_px[8]->SetParameters(1.068866e+00, -2.793530e-02);
    fcnWidth_px[9]->SetParameters(1.068866e+00, -2.793530e-02);
  } // end: if (kP)
  else {
    printf(" no correction available for Zdim = %s!\n", AliDielectronVarManager::GetValueName(corrZdim));
    printf(" no correction applied!\n");
    return;
  }
  
  if (use1Dfunctions) {
    for (int ix=1-extrapolate; ix<=hMean_2Dfit->GetNbinsX()+extrapolate; ix++) {
      for (int iy=1-extrapolate; iy<=hMean_2Dfit->GetNbinsY()+extrapolate; iy++) {
        if (isFcnZdim) {
          //cout << " hMean_2Dfit->GetBinCenter(ix) = " << hMean_2Dfit->GetXaxis()->GetBinCenter(ix) << endl;
          hMean_2Dfit->SetBinContent(ix, iy, fcnMean_px[iy]->Eval(hMean_2Dfit->GetXaxis()->GetBinCenter(ix)));
          hWidth_2Dfit->SetBinContent(ix, iy, fcnWidth_px[iy]->Eval(hWidth_2Dfit->GetXaxis()->GetBinCenter(ix)));
        } else {
          //cout << " hMean_2Dfit->GetBinCenter(iy) = " << hMean_2Dfit->GetYaxis()->GetBinCenter(iy) << endl;
          hMean_2Dfit->SetBinContent(ix, iy, fcnMean_px[ix]->Eval(hMean_2Dfit->GetYaxis()->GetBinCenter(iy)));
          hWidth_2Dfit->SetBinContent(ix, iy, fcnWidth_px[ix]->Eval(hWidth_2Dfit->GetYaxis()->GetBinCenter(iy)));
        }
      }
    }
  }
  
  if (task->IsA()==AliAnalysisTaskElectronEfficiency::Class()) {
    (static_cast<AliAnalysisTaskElectronEfficiency*>task)->SetCentroidCorrFunctionITS(hMean_2Dfit, corrZdim, corrYdim);
    (static_cast<AliAnalysisTaskElectronEfficiency*>task)->SetWidthCorrFunctionITS(hWidth_2Dfit, corrZdim, corrYdim);
    printf(" MC ITS PID eta correction loaded (AliAnalysisTaskElectronEfficiency)!\n");
  }
  else if (task->IsA()==AliDielectron::Class()) {
    (static_cast<AliDielectron*>task)->SetCentroidCorrFunctionITS(hMean_2Dfit, corrZdim, corrYdim);
    (static_cast<AliDielectron*>task)->SetWidthCorrFunctionITS(hWidth_2Dfit, corrZdim, corrYdim);
    printf(" MC ITS PID eta correction loaded (AliDielectron)!\n");
  }    
}

//_______________________________________________________________________________________________
void LMEECutLib::SetTPCSigmaEleCorrectionMC(TNamed* task, Int_t corrZdim, Int_t corrYdim) {
  //
  // MC post-correction for the centroid and width of electron sigmas in the TPC, can be one/two/three-dimensional
  //
  printf("TPC PID post-correction not needed for LHC14a1[a,b] (logical when using TuneOnData).\n"); return;
  //  printf("starting LMEECutLib::SetTPCSigmaEleCorrectionMC()\n");
  //  printf(" correction Zdim = %s\n", AliDielectronVarManager::GetValueName(corrZdim));
  //  printf(" correction Ydim = %s\n", AliDielectronVarManager::GetValueName(corrYdim));
  //
  //  printf(" MC TPC PID eta correction loaded!\n");
}



// Note: event cuts are identical for all analysis 'cutDefinition's that run together!
// the selection is hardcoded in the AddTask, currently to 'kPbPb2011_pidITSTPCTOFif_trkSPDfirst_1'
//_______________________________________________________________________________________________
AliDielectronEventCuts* LMEECutLib::GetEventCuts(Int_t cutSet, Bool_t hasMC) {
  AliDielectronEventCuts* eventCuts = 0x0;
  switch (cutSet) {
    case kPbPb2011_pidITSTPCTOFif_trkSPDfirst_1:
      //Basic Event Cuts for pp and Pb-Pb, additional cuts may be in the AddTask
      eventCuts=new AliDielectronEventCuts("eventCuts","Vertex Track && |vtxZ|<10 && ncontrib>0");
      eventCuts->SetRequireVertex();
      eventCuts->SetMinVtxContributors(1);
      eventCuts->SetVertexZ(-10.,10.);
      if (hasMC) {
        //eventCuts->SetVertexType(.....);
        //eventCuts->SetCentralityRange( 0.,90.);
      } else {
        eventCuts->SetVertexType(AliDielectronEventCuts::kVtxSPD); // AOD
        //eventCuts->SetVertexType(AliDielectronEventCuts::kVtxTPC); // kVtxAny // AOD
      }
      break;
    default: cout << "No Event Cut defined" << endl;
  }
  return eventCuts;
}


//Selection of relatively 'flat' centralities is a bit difficult...
//_______________________________________________________________________________________________
AliAnalysisCuts* LMEECutLib::GetCentralityCuts(Bool_t hasMC) {
  AliAnalysisCuts* centCuts = 0x0;
  
  // be careful with MC! the triggers should (must?) not be used.
  
  // Online Trigger MB+Semi
  AliDielectronVarCuts *trgMBorSEMICNT = new AliDielectronVarCuts("trgMBorSEMICNT","trgMBorSEMICNT");
  trgMBorSEMICNT->SetCutType(AliDielectronVarCuts::kAny);
  trgMBorSEMICNT->AddBitCut(AliDielectronVarManager::kTriggerInclONL,  1, kFALSE); // MB
  trgMBorSEMICNT->AddBitCut(AliDielectronVarManager::kTriggerInclONL,  7, kFALSE); // SEMICNT
  //trgMBorSEMICNT->AddBitCut(AliDielectronVarManager::kTriggerInclONL,  4, kFALSE); //CNT
  
  switch (selectedCentrality) {
    case kPbPb2011Central:
    case kPbPb2011_00to10:
      AliDielectronVarCuts* centRange = new AliDielectronVarCuts("centCuts","CentralityPbPb2011_00to10");
      centRange->AddCut(AliDielectronVarManager::kCentrality,0.,10.);
      centCuts = centRange;
      break;
    case kPbPb2011MidCentral:
    case kPbPb2011_10to20:
      AliDielectronVarCuts* centRange = new AliDielectronVarCuts("centCuts","CentralityPbPb2011_10to20");
      centRange->AddCut(AliDielectronVarManager::kCentrality,10.,20.);
      AliDielectronCutGroup* centCutsCG =new AliDielectronCutGroup("centCutsCG","centCutsCG",AliDielectronCutGroup::kCompAND);
      centCutsCG->AddCut(centRange);
      if (!hasMC) centCutsCG->AddCut(trgMBorSEMICNT);
      centCuts = centCutsCG;
      break;
    case kPbPb2011SemiCentral:
    case kPbPb2011_20to50:
      AliDielectronVarCuts* centRange = new AliDielectronVarCuts("centCuts","CentralityPbPb2011_20to50");
      centRange->AddCut(AliDielectronVarManager::kCentrality,20.,50.);
      centCuts = centRange;
      break;
    case kPbPb2011_00to50:
      AliDielectronVarCuts* centRange = new AliDielectronVarCuts("centCuts","CentralityPbPb2011_00to50");
      centRange->AddCut(AliDielectronVarManager::kCentrality,0.,50.);
      ////centRange->AddBitCut(AliDielectronVarManager::kTriggerExclOFF,  4, kTRUE); //this is wrong...
      ////centRange->AddBitCut(AliDielectronVarManager::kTriggerInclONL,  4, kTRUE); //also wrong...
      AliDielectronCutGroup* centCutsCG =new AliDielectronCutGroup("centCutsCG","centCutsCG",AliDielectronCutGroup::kCompAND);
      centCutsCG->AddCut(centRange);
      if (!hasMC) centCutsCG->AddCut(trgMBorSEMICNT);
      centCuts = centCutsCG;
      break;
    case kPbPb2011_10to50:
      AliDielectronVarCuts* centRange = new AliDielectronVarCuts("centCuts","CentralityPbPb2011_10to50");
      centRange->AddCut(AliDielectronVarManager::kCentrality,10.,50.);
      AliDielectronCutGroup* centCutsCG =new AliDielectronCutGroup("centCutsCG","centCutsCG",AliDielectronCutGroup::kCompAND);
      centCutsCG->AddCut(centRange);
      if (!hasMC) centCutsCG->AddCut(trgMBorSEMICNT);
      centCuts = centCutsCG;
      break;
    default: cout << "No Centrality selected" << endl;
  }
  return centCuts;
}


//Basic track rotator settings from J/Psi, more investigation needed
//_______________________________________________________________________________________________
AliDielectronTrackRotator* LMEECutLib::GetTrackRotator() {
  AliDielectronTrackRotator* trackRotator = 0x0;
  Int_t cutSet=-1;
  switch (cutSet) {
    default: cout << "No Rotator defined" << endl;
      //default:
      //  trackRotator = new AliDielectronTrackRotator();
      //  trackRotator->SetIterations(20);
      //  trackRotator->SetConeAnglePhi(TMath::Pi()/180*165);
      //  trackRotator->SetStartAnglePhi(TMath::Pi());
      //  break;
  }
  return trackRotator;
}


//_______________________________________________________________________________________________
AliDielectronMixingHandler* LMEECutLib::GetMixingHandler() {
  AliDielectronMixingHandler* mixingHandler = 0x0;
  Int_t cutSet=1;
  switch (cutSet) {
    case 1:
      mixingHandler = new AliDielectronMixingHandler();
      mixingHandler->AddVariable(AliDielectronVarManager::kZvPrim,"-10,-5,0,5,10");
      mixingHandler->AddVariable(AliDielectronVarManager::kCentrality,"0,5,10,20,30,50,80");
      // now using TPC event plane, uncorrected. (also, the old phi range was wrong, now same effective binning.)
      mixingHandler->AddVariable(AliDielectronVarManager::kTPCrpH2uc, 6, TMath::Pi()/-2., TMath::Pi()/2.);
      mixingHandler->SetDepth(15);
      mixingHandler->SetMixType(AliDielectronMixingHandler::kAll);
      break;
      //[...]
    default: cout << "No Mixer defined" << endl;
  }
  return mixingHandler;
}



//Pair Cuts for Analysis step - take care of logic - inverted compared to other PairCuts!!
// cuts = SELECTION!!!
//_______________________________________________________________________________________________
AliAnalysisCuts* LMEECutLib::GetPairCutsAna() {
  cout << " >>>>>>>>>>>>>>>>>>>>>> GetPairCutsAna() >>>>>>>>>>>>>>>>>>>>>> " << endl;
  AliDielectronVarCuts* pairVarCuts = (AliDielectronVarCuts*) GetPairCutsPre(selectedPairCutsAna);
  if (!pairVarCuts) return 0x0;
  pairVarCuts->InvertCuts();
  return pairVarCuts;
}

//Pair Cuts for PREFILTER step
// cuts = REJECTION!!!
//_______________________________________________________________________________________________
AliAnalysisCuts* LMEECutLib::GetPairCutsPre(Int_t cutSet) {  
  cout << " >>>>>>>>>>>>>>>>>>>>>> GetPairCutsPre() >>>>>>>>>>>>>>>>>>>>>> " << endl;
  // for the default function call, pick the cutSet according to 'selectedPairCutsPre', which is set via the Config file.
  if (cutSet<0) { cutSet = selectedPairCutsPre; }
  
  if (cutSet==kPairCut_OFF) {
    cout << "Pair Cuts disabled " << endl;
    return 0x0;
  }
  AliDielectronVarCuts* pairVarCuts = new AliDielectronVarCuts("pairVarCuts","pairVarCuts");
  switch (cutSet) {
      // Just opening angle
    case kPairCut_theta20:
      pairVarCuts->AddCut(AliDielectronVarManager::kOpeningAngle, 0.0, 0.02);
      break;
    case kPairCut_theta50:
      pairVarCuts->AddCut(AliDielectronVarManager::kOpeningAngle, 0.0, 0.05);
      break;
    case kPairCut_theta100:
      pairVarCuts->AddCut(AliDielectronVarManager::kOpeningAngle, 0.0, 0.10);
      break;
      
      // Mee and opening angle
    case kPairCut_mee10_theta30:
      pairVarCuts->AddCut(AliDielectronVarManager::kM, 0.0, 0.01);
      pairVarCuts->AddCut(AliDielectronVarManager::kOpeningAngle, 0.0, 0.03);
      break;
    case kPairCut_mee20_theta20:
      pairVarCuts->AddCut(AliDielectronVarManager::kM, 0.0, 0.02);
      pairVarCuts->AddCut(AliDielectronVarManager::kOpeningAngle, 0.0, 0.02);
      break;
    case kPairCut_mee20_theta50:
      pairVarCuts->AddCut(AliDielectronVarManager::kM, 0.0, 0.02);
      pairVarCuts->AddCut(AliDielectronVarManager::kOpeningAngle, 0.0, 0.05);
      break;
    case kPairCut_mee30_theta60:
      pairVarCuts->AddCut(AliDielectronVarManager::kM, 0.0, 0.03);
      pairVarCuts->AddCut(AliDielectronVarManager::kOpeningAngle, 0.0, 0.06);
      break;
    case kPairCut_mee40_theta80:
      pairVarCuts->AddCut(AliDielectronVarManager::kM, 0.0, 0.04);
      pairVarCuts->AddCut(AliDielectronVarManager::kOpeningAngle, 0.0, 0.08);
      break;
    case kPairCut_mee60_theta100:
      pairVarCuts->AddCut(AliDielectronVarManager::kM, 0.0, 0.06);
      pairVarCuts->AddCut(AliDielectronVarManager::kOpeningAngle, 0.0, 0.10);
      break;
    case kPairCut_mee200_theta300: // just for testing
      pairVarCuts->AddCut(AliDielectronVarManager::kM, 0.0, 0.20);
      pairVarCuts->AddCut(AliDielectronVarManager::kOpeningAngle, 0.0, 0.30);
      break;
      
      // Mee and phiv
    case kPairCut_phiv157_mee40:
      pairVarCuts->AddCut(AliDielectronVarManager::kM, 0.0, 0.04);
      pairVarCuts->AddCut(AliDielectronVarManager::kPhivPair, 0.5*TMath::Pi(), 3.2); 
      break;
    case kPairCut_phiv157_mee60:
      pairVarCuts->AddCut(AliDielectronVarManager::kM, 0.0, 0.06);
      pairVarCuts->AddCut(AliDielectronVarManager::kPhivPair, 0.5*TMath::Pi(), 3.2); 
      break;
    case kPairCut_phiv157_mee80:
      pairVarCuts->AddCut(AliDielectronVarManager::kM, 0.0, 0.08);
      pairVarCuts->AddCut(AliDielectronVarManager::kPhivPair, 0.5*TMath::Pi(), 3.2); 
      break;
    case kPairCut_phiv157_mee100:
      pairVarCuts->AddCut(AliDielectronVarManager::kM, 0.0, 0.10);
      pairVarCuts->AddCut(AliDielectronVarManager::kPhivPair, 0.5*TMath::Pi(), 3.2); 
      break;
    case kPairCut_phiv236_mee40:
      pairVarCuts->AddCut(AliDielectronVarManager::kM, 0.0, 0.04);
      pairVarCuts->AddCut(AliDielectronVarManager::kPhivPair, 0.75*TMath::Pi(), 3.2); 
      break;
    case kPairCut_phiv236_mee60:
      pairVarCuts->AddCut(AliDielectronVarManager::kM, 0.0, 0.06);
      pairVarCuts->AddCut(AliDielectronVarManager::kPhivPair, 0.75*TMath::Pi(), 3.2); 
      break;
    case kPairCut_phiv236_mee80:
      pairVarCuts->AddCut(AliDielectronVarManager::kM, 0.0, 0.08);
      pairVarCuts->AddCut(AliDielectronVarManager::kPhivPair, 0.75*TMath::Pi(), 3.2); 
      break;
    case kPairCut_phiv236_mee100:
      pairVarCuts->AddCut(AliDielectronVarManager::kM, 0.0, 0.10);
      pairVarCuts->AddCut(AliDielectronVarManager::kPhivPair, 0.75*TMath::Pi(), 3.2); 
      break;
      
    default: cout << "No (Prefilter) Pair Cuts defined " << endl;
  } 
  return pairVarCuts;
}



//_______________________________________________________________________________________________
AliAnalysisCuts* LMEECutLib::GetTrackCutsAna() {
  cout << " >>>>>>>>>>>>>>>>>>>>>> GetTrackCutsAna() >>>>>>>>>>>>>>>>>>>>>> " << endl;
  AliDielectronCutGroup* cgTrackCutsAna = new AliDielectronCutGroup("cgTrackCutsAna","cgTrackCutsAna",AliDielectronCutGroup::kCompAND);
  cgTrackCutsAna->AddCut( GetPIDCuts(selectedPIDAna) );
  cgTrackCutsAna->AddCut( GetQualityCuts(selectedQualityAna, kInclude) );
  cgTrackCutsAna->AddCut( GetKineCutsAna() );
  
  return cgTrackCutsAna;
}

//_______________________________________________________________________________________________
AliAnalysisCuts* LMEECutLib::GetTrackCutsPre() {
  cout << " >>>>>>>>>>>>>>>>>>>>>> GetTrackCutsPre() >>>>>>>>>>>>>>>>>>>>>> " << endl;
  AliDielectronCutGroup* cgTrackCutsPre = new AliDielectronCutGroup("cgTrackCutsPre","cgTrackCutsPre",AliDielectronCutGroup::kCompAND);
  cgTrackCutsPre->AddCut( GetPIDCuts(selectedPIDPre) );
  cgTrackCutsPre->AddCut( GetQualityCuts(selectedQualityPre, kInclude) );
  cgTrackCutsPre->AddCut( GetKineCutsPre() );
  
  // in case the prefilter cuts do not include all needed global tracks, we create an "OR" cutgroup:
  AliDielectronCutGroup* cgInitialTrackFilter = new AliDielectronCutGroup("cgInitialTrackFilter","cgInitialTrackFilter",AliDielectronCutGroup::kCompOR);
  cgInitialTrackFilter->AddCut( cgTrackCutsPre );
  cgInitialTrackFilter->AddCut( GetTrackCutsAna() );
  
  return cgInitialTrackFilter;
}



//_______________________________________________________________________________________________
AliAnalysisCuts* LMEECutLib::GetKineCutsAna() {
  cout << " >>>>>>>>>>>>>>>>>>>>>> GetKineCutsAna() >>>>>>>>>>>>>>>>>>>>>> " << endl;
  AliDielectronVarCuts* kineCuts = (AliDielectronVarCuts*) GetKineCutsPre(selectedKineCutsAna);
  if (!kineCuts) return 0x0;
  return kineCuts;
}

//_______________________________________________________________________________________________
AliAnalysisCuts* LMEECutLib::GetKineCutsPre(Int_t cutSet) {  
  cout << " >>>>>>>>>>>>>>>>>>>>>> GetKineCutsPre() >>>>>>>>>>>>>>>>>>>>>> " << endl;
  // for the default function call, pick the cutSet according to 'selectedKineCutsPre', which is set via the Config file.
  if (cutSet<0) { cutSet = selectedKineCutsPre; }
  
  AliDielectronVarCuts* kineCuts = new AliDielectronVarCuts("kineCuts","kineCuts");
  switch (cutSet) {
    case kKineCut_p50inf_eta150:
      kineCuts->AddCut(AliDielectronVarManager::kP, .05, 100.); // momentum!
      kineCuts->AddCut(AliDielectronVarManager::kEta, -1.50, 1.50);
      break;
    case kKineCut_pt50_eta080:
      kineCuts->AddCut(AliDielectronVarManager::kPt, .05, 3.5);
      kineCuts->AddCut(AliDielectronVarManager::kEta, -0.80, 0.80);
      break;
    case kKineCut_pt50_eta090:
      kineCuts->AddCut(AliDielectronVarManager::kPt, .05, 3.5);
      kineCuts->AddCut(AliDielectronVarManager::kEta, -0.90, 0.90);
      break;
    case kKineCut_pt200_eta080:
      kineCuts->AddCut(AliDielectronVarManager::kPt, .2, 3.5);
      kineCuts->AddCut(AliDielectronVarManager::kEta, -0.80, 0.80);
      break;
    case kKineCut_pt200_eta090:
      kineCuts->AddCut(AliDielectronVarManager::kPt, .2, 3.5);
      kineCuts->AddCut(AliDielectronVarManager::kEta, -0.90, 0.90);
      break;
    case kKineCut_pt300_eta080:
      kineCuts->AddCut(AliDielectronVarManager::kPt, .3, 3.5);
      kineCuts->AddCut(AliDielectronVarManager::kEta, -0.80, 0.80);
      break;
    case kKineCut_pt300_eta090:
      kineCuts->AddCut(AliDielectronVarManager::kPt, .3, 3.5);
      kineCuts->AddCut(AliDielectronVarManager::kEta, -0.90, 0.90);
      break;
    case kKineCut_pt400_eta080:
      kineCuts->AddCut(AliDielectronVarManager::kPt, .4, 3.5);
      kineCuts->AddCut(AliDielectronVarManager::kEta, -0.80, 0.80);
      break;
    case kKineCut_pt400_eta090:
      kineCuts->AddCut(AliDielectronVarManager::kPt, .4, 3.5);
      kineCuts->AddCut(AliDielectronVarManager::kEta, -0.90, 0.90);
      break;
      
    default: cout << "No (Prefilter) Kine Cuts defined " << endl;
  } 
  return kineCuts;
}


//_______________________________________________________________________________________________
AliAnalysisCuts* LMEECutLib::GetPIDCuts(Int_t cutSet) {
  cout << " >>>>>>>>>>>>>>>>>>>>>> GetPIDCuts() >>>>>>>>>>>>>>>>>>>>>> " << endl;
  //-----------------------------------------------
  // PID cuts depend on TPC_inner_p, if not specified
  // PID cut ranges correspond to global momentum P
  // check it again!!!
  //-----------------------------------------------
  
  if (cutSet<0) {
    cout << "Invalid cutSet selected! " << endl;
    return 0x0;
  }
  
  AliDielectronPID *pidCuts = new AliDielectronPID("pidCuts","pidCuts");
  switch (cutSet) {
      
      // made identical to my default PID
    case kCut01: // (was pion<4)
    case kCut05: // (was ITS<0.5, pion<4)
    case kCut18: // (exactly kept)
      // loose TPC
      pidCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron,  -4. ,1. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    );
      pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -1.5,3. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    );
      pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,3. ,0.0, 100., kTRUE ,AliDielectronPID::kRequire    );
      pidCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,  -3. ,3. ,0.0, 100., kFALSE,AliDielectronPID::kIfAvailable);
      return pidCuts;
      break;

    case kCut02: // (was TPC<4)
      // tight TPC, tight TOF
      pidCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron,  -3. ,1. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    );
      pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -0.5,3. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    );
      pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,3. ,0.0, 100., kTRUE ,AliDielectronPID::kRequire    );
      pidCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,  -2. ,2. ,0.0, 100., kFALSE,AliDielectronPID::kIfAvailable);
      return pidCuts;
      break;

    case kCut03: // (exactly kept)
      // tight TOF
      pidCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron,  -4. ,1. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    );
      pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -1. ,3. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    );
      pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,4. ,0.0, 100., kTRUE ,AliDielectronPID::kRequire    );
      pidCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,  -2. ,2. ,0.0, 100., kFALSE,AliDielectronPID::kIfAvailable);
      return pidCuts;
      break;

    case kCut04: // (was ITS -3.5..0.5)
      // loose TPC, loose TOF
      pidCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron,  -3. ,1. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    );
      pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -1.5,3. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    );
      pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,3. ,0.0, 100., kTRUE ,AliDielectronPID::kRequire    );
      pidCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,  -4. ,4. ,0.0, 100., kFALSE,AliDielectronPID::kIfAvailable);
      return pidCuts;
      break;

    case kCut06: // (was ITS -3..1, TPC<4, pion<3.5)
    case kCut19: // (exactly kept)
      // loose ITS, loose TPC
      pidCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron,  -4. ,1.5,0.0, 100., kFALSE,AliDielectronPID::kRequire    );
      pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -1.5,3. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    );
      pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,4. ,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ); // makes sense here
      pidCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,  -3. ,3. ,0.0, 100., kFALSE,AliDielectronPID::kIfAvailable);
      return pidCuts;
      break;

    case kCut07:
      // loose ITS, tight TOF
      pidCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron,  -4. ,1.5,0.0, 100., kFALSE,AliDielectronPID::kRequire    );
      pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -1. ,3. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    );
      pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,3. ,0.0, 100., kTRUE ,AliDielectronPID::kRequire    );
      pidCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,  -2. ,2. ,0.0, 100., kFALSE,AliDielectronPID::kIfAvailable);
      return pidCuts;
      break;

    case kCut08: // (was ITS<1, TPC<4, pion<4)
      // loose ITS
      pidCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron,  -3. ,1.5,0.0, 100., kFALSE,AliDielectronPID::kRequire    );
      pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -1. ,3. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    );
      pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,3. ,0.0, 100., kTRUE ,AliDielectronPID::kRequire    );
      pidCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,  -3. ,3. ,0.0, 100., kFALSE,AliDielectronPID::kIfAvailable);
      return pidCuts;
      break;

    case kCut09: // (was ITS<0, TPC<4, pion<4)
      // tight ITS, loose TPC
      pidCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron,  -3. ,0.5,0.0, 100., kFALSE,AliDielectronPID::kRequire    );
      pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -1.5,3. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    );
      pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,3. ,0.0, 100., kTRUE ,AliDielectronPID::kRequire    );
      pidCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,  -3. ,3. ,0.0, 100., kFALSE,AliDielectronPID::kIfAvailable);
      return pidCuts;
      break;

    case kCut10: // (was ITS<0)
    case kCut20: // (exactly kept)
      // tight ITS, tight TPC
      pidCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron,  -3. ,0.5,0.0, 100., kFALSE,AliDielectronPID::kRequire    );
      pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -0.5,3. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    );
      pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,4. ,0.0, 100., kTRUE ,AliDielectronPID::kRequire    );
      pidCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,  -3. ,3. ,0.0, 100., kFALSE,AliDielectronPID::kIfAvailable);
      return pidCuts;
      break;

    case kCut11: // (was ITS<0.5. set pions to <4, they are cut out anyhow)
      // very tight TPC
      pidCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron,  -3. ,1. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    );
      pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -0.5,2. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    );
      pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,4. ,0.0, 100., kTRUE ,AliDielectronPID::kRequire    );
      pidCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,  -3. ,3. ,0.0, 100., kFALSE,AliDielectronPID::kIfAvailable);
      return pidCuts;
      break;

    case kCut12: // (set pions to <4, they are cut out anyhow)
      // tight ITS, very tight TPC, tight TOF
      pidCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron,  -3. ,0.5,0.0, 100., kFALSE,AliDielectronPID::kRequire    );
      pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -0.5,2. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    );
      pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,4. ,0.0, 100., kTRUE ,AliDielectronPID::kRequire    );
      pidCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,  -2. ,2. ,0.0, 100., kFALSE,AliDielectronPID::kIfAvailable);
      return pidCuts;
      break;

    case kCut13: // (was pion<4)
    case kCut15: // (was pion<4)
      // tight ITS
      pidCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron,  -4. ,0.5,0.0, 100., kFALSE,AliDielectronPID::kRequire    );
      pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -1. ,3. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    );
      pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,3. ,0.0, 100., kTRUE ,AliDielectronPID::kRequire    );
      pidCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,  -3. ,3. ,0.0, 100., kFALSE,AliDielectronPID::kIfAvailable);
      return pidCuts;
      break;

    case kCut14:
      // tight ITS, tight TOF
      pidCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron,  -3. ,0.5,0.0, 100., kFALSE,AliDielectronPID::kRequire    );
      pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -1. ,3. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    );
      pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,4. ,0.0, 100., kTRUE ,AliDielectronPID::kRequire    );
      pidCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,  -2. ,2. ,0.0, 100., kFALSE,AliDielectronPID::kIfAvailable);
      return pidCuts;
      break;

    case kCut16: // (exactly kept)
    case kCut17: // (exactly kept)
      // average PID
      pidCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron,  -4. ,1. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    );
      pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -1. ,3. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    );
      pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,4. ,0.0, 100., kTRUE ,AliDielectronPID::kRequire    );
      pidCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,  -3. ,3. ,0.0, 100., kFALSE,AliDielectronPID::kIfAvailable);
      return pidCuts;
      break;




      ///
      /// ITS + TPC + TOFif
      ///
      // ===== Default PID ITS+TPC+TOFif =====
    case kPbPb2011_pidITSTPCTOFif_trkSPDorSDD_7_V0excl:
    case kPbPb2011_pidITSTPCTOFif_trkSPDfirst_7_V0excl:
    case kPbPb2011MC_pi0Dal_1:
    case kPbPb2011_pidITSTPCTOFif_trkSPD5orSDD4cls_4:
    case kPbPb2011_pidITSTPCTOFif_trkSPDfirst5cls_4:
    case kPbPb2011_pidITSTPCTOFif_trkSPDorSDD_1:
    case kPbPb2011_pidITSTPCTOFif_trkSPDfirst_1:
    case kPbPb2011PID_ITSTPCTOFif_1:
      //TPC: electron inclusion asymmetric
      //     pion     exclusion 3sigma
      //ITS: electron inclusion asymmetric OVER FULL MOMENTUM RANGE
      //TOF: electron inclusion 3sigma - BUT ONLY IF AVAILABLE
      pidCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -4. , 1. , 0. ,100., kFALSE);
      pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -1.5, 3. , 0. ,100., kFALSE);
      pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion,     -3. , 3. , 0. ,100., kTRUE);
      pidCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3. , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);
      return pidCuts;
      break;
      
      // ===== Tighter PID ITS+TPC+TOFif =====
    case kPbPb2011_pidITS2gevTPCTOFif_trkSPD5orSDD4cls_6_tight:
    case kPbPb2011_pidITS2gevTPCTOFif_trkSPDfirst5cls_6_tight:
    case kPbPb2011_pidITS2gevTPCTOFif_trkSPDorSDD_5_tight:
    case kPbPb2011_pidITS2gevTPCTOFif_trkSPDfirst_5_tight:
    case kPbPb2011PID_ITSTPCTOFif_2:
      //ITS: electron inclusion asymmetric in region where p,K cross electrons in TPC
      pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -1.5, 2.5, 0. ,100., kFALSE);
      pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion,     -3. , 3. , 0. ,100., kTRUE);
      pidCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -4. , 0.5, 0. ,  2., kFALSE);
      pidCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -2. , 2. , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);
      return pidCuts;
      break;
      
      // ===== Relaxed PID ITS+TPC+TOFif =====
    case kPbPb2011PID_ITSTPCTOFif_3:
      //ITS: electron inclusion asymmetric in region where p,K cross electrons in TPC
      pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2. , 3., 0. ,100., kFALSE);
      pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion,     -3. , 3., 0. ,100., kTRUE);
      pidCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -4. , 2., 0. ,  2., kFALSE);
      pidCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3., 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);
      return pidCuts;
      break;
      
      // ===== LOOSE PID ITS+TPC+TOFif =====
    case kPbPb2011_pidITSTPCTOFif_trkSPDorSDD_2_loose: // loose "ITSTPCTOFif" PID - for 2D contamination study
    case kPbPb2011_pidITSTPCTOFif_trkSPDfirst_2_loose:
    case kPbPb2011PID_ITSTPCTOFif_LOOSE:
      pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-12. ,20. , 0. ,100., kFALSE);
      pidCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron,-10. ,20. , 0. ,100., kFALSE);
      pidCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3. , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);
      return pidCuts;
      break;
    case kPbPb2011PID_ITSTPCTOFif_looseTPC: // only loose TPC - for Anisas contamination & purity tool
      pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-12. ,20. , 0. ,100., kFALSE);
      pidCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.5, 1.5, 0. ,100., kFALSE);
      pidCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3. , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);
      return pidCuts;
      break;
      
      
      ///
      /// ITS+TPC
      ///
    case kPbPb2011_pidITSTPC_trkSPDfirst_3:
    case kPbPb2011PID_TPCITS_1:
      pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -1.5, 3. , 0. ,100., kFALSE);
      pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion,     -3. , 3. , 0. ,100., kTRUE);
      pidCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -4. , 1. , 0. ,100., kFALSE);
      return pidCuts;
      break;
    case kPbPb2011PID_TPCITSif_2:
      pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2. , 3. , 0. ,100., kFALSE);
      pidCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3. , 2. , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);
      return pidCuts;
      break;
    case kPbPb2011PID_TPCITS_3: // for full efficiency. used in MC for resolution extraction.
      pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -3. , 3. , 0. ,100., kFALSE);
      pidCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -4. , 3. , 0. ,100., kFALSE);
      return pidCuts;
      break;
      
      
      ///
      /// TPC + TOF
      ///
      // ===== Default PID TPC+TOF =====
    case kPbPb2011_pidTPCTOF_trkSPDorSDD_1:
    case kPbPb2011_pidTPCTOF_trkSPDfirst_1:
    case kPbPb2011PID_TPCTOF_1:
      //TPC: electron inclusion asymmetric
      //     pion     exclusion 3sigma
      //TOF: electron inclusion 3sigma in region where p,K cross electrons in TPC
      //     (may be useful to use TOF with pt instead of p, because pt determines its acceptance.)
      pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -1.5, 3., 0. ,100., kFALSE);
      pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion,     -3. , 3., 0. ,100., kTRUE);
      pidCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3., 0. ,1.7 , kFALSE);
      return pidCuts;
      break;
      // ===== LOOSE PID TPC+TOF =====
    case kPbPb2011_pidTPCTOF_trkSPDorSDD_2_loose: // loose "ITSTPCTOFif" PID - for 1D contamination study in TPC
    case kPbPb2011_pidTPCTOF_trkSPDfirst_2_loose:
    case kPbPb2011PID_TPCTOF_LOOSE:
      pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-12. ,20. , 0. ,100., kFALSE);
      pidCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3. , 0. ,1.7 , kFALSE);
      return pidCuts;
      break;
      
      
      ///
      /// TPC (+TOFif)
      ///
    case kPbPb2011_pidTPC_trkSPDfirst_3:
    case kPbPb2011PID_TPC_1:
      pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -1.5, 3., 0. ,100., kFALSE);
      pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion,     -3. , 3., 0. ,100., kTRUE);
      return pidCuts;
      break;
    case kPbPb2011PID_TPC_pre:
    case kPbPb2011PID_TPC_2:
      pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -3. , 3., 0. ,100., kFALSE);
      pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion,     -3. , 3., 0. ,100., kTRUE);
      pidCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3., 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);
      return pidCuts;
      break;
      
      ///
      /// ITS
      ///
    case kPbPb2011PID_ITS_1:
      pidCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3. , 3.);
      return pidCuts;
      break;
      
      ///
      /// Other
      ///
      ///
      /// PID for V0 tasks
    case kPbPb2011PID_V0_1_pionCont:
      pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion,     -3. , 3. , 0. ,100., kFALSE);
      // no break here to add the cuts below...
    case kPbPb2011PID_V0_1:
      pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-12. ,20. , 0. ,100., kFALSE);
      pidCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -1.5, 1.5, 0. ,100., kFALSE);
      return pidCuts;
      break;
    case kPbPb2011PID_V0_2_TOFif_pionCont:
      pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion,     -3. , 3. , 0. ,100., kFALSE);
      // no break here to add the cuts below...
    case kPbPb2011PID_V0_2_TOFif:
      pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-12. ,20. , 0. ,100., kFALSE);
      pidCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3. , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);
      // helps to see what gets rejected...
      return pidCuts;
      break;
      
      // NO PID
    case kPbPb2011PID_OFF:
      return pidCuts;
      break;
      
    default: cout << "No PID Cut defined " << endl;
      return 0x0;
  }
  return 0x0;
  
}


//Make/Tighten track Cuts that are *NOT* already
//done in the AOD production
//**IMPORTANT**: For AODs, select FilterBit
//the method is ignored for ESDs

//_______________________________________________________________________________________________
AliAnalysisCuts* LMEECutLib::GetQualityCuts(Int_t cutSet, Int_t doExclusion) {
  cout << " >>>>>>>>>>>>>>>>>>>>>> GetQualityCuts() >>>>>>>>>>>>>>>>>>>>>> " << endl;
  AliDielectronCutGroup* trackCuts=0x0;
  
  if (cutSet<0) {
    cout << "Invalid cutSet selected! " << endl;
    return 0x0;
  }
  
  
  // define trackCuts
  AliDielectronTrackCuts *trackBit4SPDfirst = new AliDielectronTrackCuts("trackBit4SPDfirst","trackBit4SPDfirst");
  trackBit4SPDfirst->SetAODFilterBit(1<<4); //GetStandardITSTPCTrackCuts2011(kFALSE); loose DCA, SPD any
  trackBit4SPDfirst->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
  //
  AliDielectronTrackCuts *trackBit4SPDboth = new AliDielectronTrackCuts("trackBit4SPDboth","trackBit4SPDboth");
  trackBit4SPDboth->SetAODFilterBit(1<<4); //GetStandardITSTPCTrackCuts2011(kFALSE); loose DCA, SPD any
  trackBit4SPDboth->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kBoth);
  //
  AliDielectronTrackCuts *trackBit6SDDfirst = new AliDielectronTrackCuts("trackBit6SDDfirst","trackBit6SDDfirst");
  trackBit6SDDfirst->SetAODFilterBit(1<<6); // GetStandardITSTPCTrackCuts2011(kTRUE), SPD none, SDD first
  
  // define varCuts and add cuts which are identical for all cutSets
  // main track sample:
  AliDielectronVarCuts* varCutsMain = new AliDielectronVarCuts("varCutsMain","varCutsMain");
  varCutsMain->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
  varCutsMain->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
  // complimentary track sample:
  AliDielectronVarCuts* varCutsCompliment = new AliDielectronVarCuts("varCutsCompliment","varCutsCompliment");
  varCutsCompliment->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
  varCutsCompliment->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
  
  // define cutGroups
  AliDielectronCutGroup* cutGroupMain       = new AliDielectronCutGroup("cutGroupMain","cutGroupMain",AliDielectronCutGroup::kCompAND);
  AliDielectronCutGroup* cutGroupCompliment = new AliDielectronCutGroup("cutGroupCompliment","cutGroupCompliment",AliDielectronCutGroup::kCompAND);
  // the two track selections can be combined with "OR" condition using this cutgroup:
  AliDielectronCutGroup* cutGroupCombineOR  = new AliDielectronCutGroup("cutGroupCombineOR","cutGroupCombineOR",AliDielectronCutGroup::kCompOR);

  
  switch (cutSet) {
      
      // made identical to my default track quality cuts
    case kCut01:
    case kCut10:
      varCutsMain->AddCut(AliDielectronVarManager::kNclsITS,          4.0, 100.0);
      varCutsMain->AddCut(AliDielectronVarManager::kNFclsTPCr,      100.0, 160.0);
      varCutsMain->AddCut(AliDielectronVarManager::kNFclsTPCfCross,   0.8,   1.1); // Theo has 0.5, not possible in AOD
      varCutsMain->AddCut(AliDielectronVarManager::kTPCchi2Cl,        0.0,   4.0);
      //trackCuts->SetMaxChi2PerClusterITS(4.5);
      //trackCuts->SetMinNClustersTPC(80);
      cutGroupMain->AddCut(varCutsMain);
      cutGroupMain->AddCut(trackBit4SPDfirst);
      if (fIsESDTask>0) cutGroupMain->AddCut(GetESDTrackCutsAna(kPbPb2011TRK_FilterBit4));
      trackCuts = cutGroupMain;
      break;
      
    case kCut02: // tight ITS, tight TPC
      varCutsMain->AddCut(AliDielectronVarManager::kNclsITS,          5.0, 100.0);
      varCutsMain->AddCut(AliDielectronVarManager::kNFclsTPCr,      130.0, 160.0);
      varCutsMain->AddCut(AliDielectronVarManager::kNFclsTPCfCross,   0.9,   1.1);
      varCutsMain->AddCut(AliDielectronVarManager::kTPCchi2Cl,        0.0,   3.0);
      //trackCuts->SetMaxChi2PerClusterITS(3.5);
      //trackCuts->SetMinNClustersTPC(100);
      cutGroupMain->AddCut(varCutsMain);
      cutGroupMain->AddCut(trackBit4SPDfirst);
      if (fIsESDTask>0) cutGroupMain->AddCut(GetESDTrackCutsAna(kPbPb2011TRK_FilterBit4));
      trackCuts = cutGroupMain;
      break;
      
    case kCut03:
      varCutsMain->AddCut(AliDielectronVarManager::kNclsITS,          4.0, 100.0);
      varCutsMain->AddCut(AliDielectronVarManager::kNFclsTPCr,       80.0, 160.0);
      varCutsMain->AddCut(AliDielectronVarManager::kNFclsTPCfCross,   0.8,   1.1); // Theo has 0.5, not possible in AOD
      varCutsMain->AddCut(AliDielectronVarManager::kTPCchi2Cl,        0.0,   3.0);
      //trackCuts->SetMaxChi2PerClusterITS(3.5);
      //trackCuts->SetMinNClustersTPC(100);
      cutGroupMain->AddCut(varCutsMain);
      cutGroupMain->AddCut(trackBit4SPDboth);
      if (fIsESDTask>0) cutGroupMain->AddCut(GetESDTrackCutsAna(kPbPb2011TRK_FilterBit4));
      trackCuts = cutGroupMain;
      break;
      
    case kCut04:
      varCutsMain->AddCut(AliDielectronVarManager::kNclsITS,          4.0, 100.0);
      varCutsMain->AddCut(AliDielectronVarManager::kNFclsTPCr,      130.0, 160.0);
      varCutsMain->AddCut(AliDielectronVarManager::kNFclsTPCfCross,   0.9,   1.1); // Theo has 0.7
      varCutsMain->AddCut(AliDielectronVarManager::kTPCchi2Cl,        0.0,   4.0);
      //trackCuts->SetMaxChi2PerClusterITS(4.5);
      //trackCuts->SetMinNClustersTPC(120);
      cutGroupMain->AddCut(varCutsMain);
      cutGroupMain->AddCut(trackBit4SPDfirst);
      if (fIsESDTask>0) cutGroupMain->AddCut(GetESDTrackCutsAna(kPbPb2011TRK_FilterBit4));
      trackCuts = cutGroupMain;
      break;
      
    case kCut05: // Theo has 6 ITS clu. I use SPDboth instead...
      varCutsMain->AddCut(AliDielectronVarManager::kNclsITS,          5.0, 100.0);
      varCutsMain->AddCut(AliDielectronVarManager::kNFclsTPCr,       80.0, 160.0);
      varCutsMain->AddCut(AliDielectronVarManager::kNFclsTPCfCross,   0.9,   1.1);
      varCutsMain->AddCut(AliDielectronVarManager::kTPCchi2Cl,        0.0,   3.0);
      //trackCuts->SetMaxChi2PerClusterITS(2.5);
      //trackCuts->SetMinNClustersTPC(80);
      cutGroupMain->AddCut(varCutsMain);
      cutGroupMain->AddCut(trackBit4SPDboth);
      if (fIsESDTask>0) cutGroupMain->AddCut(GetESDTrackCutsAna(kPbPb2011TRK_FilterBit4));
      trackCuts = cutGroupMain;
      break;
      
    case kCut06:
    case kCut13: // Theo has 6 ITS clu
      varCutsMain->AddCut(AliDielectronVarManager::kNclsITS,          5.0, 100.0);
      varCutsMain->AddCut(AliDielectronVarManager::kNFclsTPCr,      100.0, 160.0);
      varCutsMain->AddCut(AliDielectronVarManager::kNFclsTPCfCross,   0.8,   1.1); // Theo has 0.5, not possible in AOD
      varCutsMain->AddCut(AliDielectronVarManager::kTPCchi2Cl,        0.0,   4.0);
      //trackCuts->SetMaxChi2PerClusterITS(4.5);
      //trackCuts->SetMinNClustersTPC(100);
      cutGroupMain->AddCut(varCutsMain);
      cutGroupMain->AddCut(trackBit4SPDfirst);
      if (fIsESDTask>0) cutGroupMain->AddCut(GetESDTrackCutsAna(kPbPb2011TRK_FilterBit4));
      trackCuts = cutGroupMain;
      break;
      
    case kCut07:
      varCutsMain->AddCut(AliDielectronVarManager::kNclsITS,          5.0, 100.0);
      varCutsMain->AddCut(AliDielectronVarManager::kNFclsTPCr,      100.0, 160.0);
      varCutsMain->AddCut(AliDielectronVarManager::kNFclsTPCfCross,   0.8,   1.1); // Theo has 0.5, not possible in AOD
      varCutsMain->AddCut(AliDielectronVarManager::kTPCchi2Cl,        0.0,   3.0);
      //trackCuts->SetMaxChi2PerClusterITS(4.5);
      //trackCuts->SetMinNClustersTPC(100);
      cutGroupMain->AddCut(varCutsMain);
      cutGroupMain->AddCut(trackBit4SPDfirst);
      if (fIsESDTask>0) cutGroupMain->AddCut(GetESDTrackCutsAna(kPbPb2011TRK_FilterBit4));
      trackCuts = cutGroupMain;
      break;
      
    case kCut08:
      varCutsMain->AddCut(AliDielectronVarManager::kNclsITS,          5.0, 100.0);
      varCutsMain->AddCut(AliDielectronVarManager::kNFclsTPCr,      100.0, 160.0);
      varCutsMain->AddCut(AliDielectronVarManager::kNFclsTPCfCross,   0.9,   1.1);
      varCutsMain->AddCut(AliDielectronVarManager::kTPCchi2Cl,        0.0,   3.0);
      //trackCuts->SetMaxChi2PerClusterITS(3.5);
      //trackCuts->SetMinNClustersTPC(120);
      cutGroupMain->AddCut(varCutsMain);
      cutGroupMain->AddCut(trackBit4SPDboth);
      if (fIsESDTask>0) cutGroupMain->AddCut(GetESDTrackCutsAna(kPbPb2011TRK_FilterBit4));
      trackCuts = cutGroupMain;
      break;
      
    case kCut09:
      varCutsMain->AddCut(AliDielectronVarManager::kNclsITS,          4.0, 100.0);
      varCutsMain->AddCut(AliDielectronVarManager::kNFclsTPCr,       80.0, 160.0);
      varCutsMain->AddCut(AliDielectronVarManager::kNFclsTPCfCross,   0.9,   1.1);
      varCutsMain->AddCut(AliDielectronVarManager::kTPCchi2Cl,        0.0,   3.0);
      //trackCuts->SetMaxChi2PerClusterITS(4.5);
      //trackCuts->SetMinNClustersTPC(120);
      cutGroupMain->AddCut(varCutsMain);
      cutGroupMain->AddCut(trackBit4SPDfirst);
      if (fIsESDTask>0) cutGroupMain->AddCut(GetESDTrackCutsAna(kPbPb2011TRK_FilterBit4));
      trackCuts = cutGroupMain;
      break;
      
    case kCut11:
    case kCut12:
      varCutsMain->AddCut(AliDielectronVarManager::kNclsITS,          5.0, 100.0);
      varCutsMain->AddCut(AliDielectronVarManager::kNFclsTPCr,      120.0, 160.0);
      varCutsMain->AddCut(AliDielectronVarManager::kNFclsTPCfCross,   0.9,   1.1); // Theo has 0.7
      varCutsMain->AddCut(AliDielectronVarManager::kTPCchi2Cl,        0.0,   3.0);
      //trackCuts->SetMinNClustersTPC(100);
      cutGroupMain->AddCut(varCutsMain);
      cutGroupMain->AddCut(trackBit4SPDfirst);
      if (fIsESDTask>0) cutGroupMain->AddCut(GetESDTrackCutsAna(kPbPb2011TRK_FilterBit4));
      trackCuts = cutGroupMain;
      break;

      //
      // cutSets 14-20 use combined track samples:
      //

    case kCut14: // tight ITS, tight TOF
    case kCut16: // average PID (like 17)
      varCutsMain->AddCut(AliDielectronVarManager::kNclsITS,          4.0, 100.0);
      varCutsMain->AddCut(AliDielectronVarManager::kNFclsTPCr,      100.0, 160.0);
      varCutsMain->AddCut(AliDielectronVarManager::kNFclsTPCfCross,   0.8,   1.1);
      varCutsMain->AddCut(AliDielectronVarManager::kTPCchi2Cl,        0.0,   4.0);
      cutGroupMain->AddCut(varCutsMain);
      cutGroupMain->AddCut(trackBit4SPDfirst);
      if (fIsESDTask>0) cutGroupMain->AddCut(GetESDTrackCutsAna(kPbPb2011TRK_FilterBit4));
      
      varCutsCompliment->AddCut(AliDielectronVarManager::kNclsITS,          4.0, 100.0);
      varCutsCompliment->AddCut(AliDielectronVarManager::kNFclsTPCr,      100.0, 160.0);
      varCutsCompliment->AddCut(AliDielectronVarManager::kNFclsTPCfCross,   0.8,   1.1);
      varCutsCompliment->AddCut(AliDielectronVarManager::kTPCchi2Cl,        0.0,   4.0);
      cutGroupCompliment->AddCut(varCutsCompliment);
      cutGroupCompliment->AddCut(trackBit6SDDfirst);
      if (fIsESDTask>0) cutGroupCompliment->AddCut(GetESDTrackCutsAna(kPbPb2011TRK_FilterBit6));
      
      cutGroupCombineOR->AddCut(cutGroupMain);
      cutGroupCombineOR->AddCut(cutGroupCompliment);
      trackCuts = cutGroupCombineOR;
      break;
      
    case kCut15: // tight ITS
    case kCut18: // loose TPC (default PID)
      varCutsMain->AddCut(AliDielectronVarManager::kNclsITS,          4.0, 100.0);
      varCutsMain->AddCut(AliDielectronVarManager::kNFclsTPCr,      100.0, 160.0);
      varCutsMain->AddCut(AliDielectronVarManager::kNFclsTPCfCross,   0.8,   1.1);
      varCutsMain->AddCut(AliDielectronVarManager::kTPCchi2Cl,        0.0,   4.0);
      cutGroupMain->AddCut(varCutsMain);
      cutGroupMain->AddCut(trackBit4SPDfirst);
      if (fIsESDTask>0) cutGroupMain->AddCut(GetESDTrackCutsAna(kPbPb2011TRK_FilterBit4));
      
      varCutsCompliment->AddCut(AliDielectronVarManager::kNclsITS,          3.0, 100.0);
      varCutsCompliment->AddCut(AliDielectronVarManager::kNFclsTPCr,      100.0, 160.0);
      varCutsCompliment->AddCut(AliDielectronVarManager::kNFclsTPCfCross,   0.8,   1.1);
      varCutsCompliment->AddCut(AliDielectronVarManager::kTPCchi2Cl,        0.0,   4.0);
      cutGroupCompliment->AddCut(varCutsCompliment);
      cutGroupCompliment->AddCut(trackBit6SDDfirst);
      if (fIsESDTask>0) cutGroupCompliment->AddCut(GetESDTrackCutsAna(kPbPb2011TRK_FilterBit6));
      
      cutGroupCombineOR->AddCut(cutGroupMain);
      cutGroupCombineOR->AddCut(cutGroupCompliment);
      trackCuts = cutGroupCombineOR;
      break;
      
    case kCut17: // average PID (like 16)
    case kCut19: // loose ITS, loose TPC
    case kCut20: // tight ITS, tight TPC
      varCutsMain->AddCut(AliDielectronVarManager::kNclsITS,          5.0, 100.0);
      varCutsMain->AddCut(AliDielectronVarManager::kNFclsTPCr,      100.0, 160.0);
      varCutsMain->AddCut(AliDielectronVarManager::kNFclsTPCfCross,   0.8,   1.1);
      varCutsMain->AddCut(AliDielectronVarManager::kTPCchi2Cl,        0.0,   4.0);
      cutGroupMain->AddCut(varCutsMain);
      cutGroupMain->AddCut(trackBit4SPDfirst);
      if (fIsESDTask>0) cutGroupMain->AddCut(GetESDTrackCutsAna(kPbPb2011TRK_FilterBit4));
      
      varCutsCompliment->AddCut(AliDielectronVarManager::kNclsITS,          4.0, 100.0);
      varCutsCompliment->AddCut(AliDielectronVarManager::kNFclsTPCr,      120.0, 160.0);
      varCutsCompliment->AddCut(AliDielectronVarManager::kNFclsTPCfCross,   0.9,   1.1);
      varCutsCompliment->AddCut(AliDielectronVarManager::kTPCchi2Cl,        0.0,   3.0);
      cutGroupCompliment->AddCut(varCutsCompliment);
      cutGroupCompliment->AddCut(trackBit6SDDfirst);
      if (fIsESDTask>0) cutGroupCompliment->AddCut(GetESDTrackCutsAna(kPbPb2011TRK_FilterBit6));
      
      cutGroupCombineOR->AddCut(cutGroupMain);
      cutGroupCombineOR->AddCut(cutGroupCompliment);
      trackCuts = cutGroupCombineOR;
      break;




      // all cutGroups objects used below are not properly created with "AliDielectronCutGroup* name = ", but it works...
    case kPbPb2011_pidITSTPCTOFif_trkSPDorSDD_7_V0excl:
      // combine 
      cgTrackCutsAna = new AliDielectronCutGroup("cgTrackCutsAna","cgTrackCutsAna",AliDielectronCutGroup::kCompAND);
      cgTrackCutsAna->AddCut(GetQualityCuts(kPbPb2011_pidITSTPCTOFif_trkSPDorSDD_1));
      cgTrackCutsAna->AddCut(GetQualityCuts(kPbPb2011TRK_V0_2_loose, kExclude));
      trackCuts = cgTrackCutsAna;
      break;
      
    case kPbPb2011_pidITSTPCTOFif_trkSPDfirst_7_V0excl:
      // combine 
      cgTrackCutsAna = new AliDielectronCutGroup("cgTrackCutsAna","cgTrackCutsAna",AliDielectronCutGroup::kCompAND);
      cgTrackCutsAna->AddCut(GetQualityCuts(kPbPb2011_pidITSTPCTOFif_trkSPDfirst_1));
      cgTrackCutsAna->AddCut(GetQualityCuts(kPbPb2011TRK_V0_2_loose, kExclude));
      //// instead, only FOR DOING A COMPARISON: ////
      ////      cgTrackCutsAna->AddCut(GetQualityCuts(kPbPb2011_pidITSTPCTOFif_trkSPDorSDD_1));
      trackCuts = cgTrackCutsAna;
      break;
      
    case kPbPb2011TRK_V0_1:
      // primarily used for selection of V0s for TPC eta correction. with additional track cuts.
      // inspired by "AliDielectronV0Cuts.cxx"
      AliDielectronV0Cuts *gammaV0Cuts = new AliDielectronV0Cuts("gammaV0Cuts","gammaV0Cuts");
      // which V0 finder you want to use
      gammaV0Cuts->SetV0finder(AliDielectronV0Cuts::kOnTheFly);  // kAll(default), kOffline or kOnTheFly // kOnTheFly better for small radii (quote Julian)
      // add some pdg codes (they are used then by the KF package and important for gamma conversions)
      gammaV0Cuts->SetPdgCodes(22,11,11); // mother, daughter1 and 2
      // add default PID cuts (defined in AliDielectronPID)
      // requirement can be set to at least one(kAny) of the tracks or to both(kBoth)
      //gammaV0Cuts->SetDefaultPID(16, AliDielectronV0Cuts::kAny);
      // add the pair cuts for V0 candidates
      // variations from Julian in ( ... )
      gammaV0Cuts->AddCut(AliDielectronVarManager::kCosPointingAngle, TMath::Cos(0.02),   1.0,  kFALSE); // ( 0.02 -- 0.05 )
      gammaV0Cuts->AddCut(AliDielectronVarManager::kChi2NDF,                       0.0,  10.0,  kFALSE);
      gammaV0Cuts->AddCut(AliDielectronVarManager::kLegDist,                       0.0,   0.25, kFALSE);
      gammaV0Cuts->AddCut(AliDielectronVarManager::kR,                             3.0,  90.0,  kFALSE);
      gammaV0Cuts->AddCut(AliDielectronVarManager::kPsiPair,                       0.0,   0.05, kFALSE); // ( 0.05 -- 0.2 )
      gammaV0Cuts->AddCut(AliDielectronVarManager::kM,                             0.0,   0.05, kFALSE); // ( 0.05 -- 0.1 )
      gammaV0Cuts->AddCut(AliDielectronVarManager::kArmPt,                         0.0,   0.05, kFALSE);
      //gammaV0Cuts->AddCut(AliDielectronVarManager::kArmAlpha,                     -0.35,  0.35, kFALSE); // ( |0.35|   ------- 0.3 ) // should increase purity...
      // selection or rejection of V0 tracks
      if (doExclusion==kExclude)  gammaV0Cuts->SetExcludeTracks(kTRUE);
      else                        gammaV0Cuts->SetExcludeTracks(kFALSE);
      // add the V0cuts directly to the track filter or to some cut group of it
      
      AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
      trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,     100.0, 160.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.8, 1.1);
      // some more default cuts automatically set by AliDielectronV0Cuts::InitEvent() !!!
      
      cgTrackCutsV0select = new AliDielectronCutGroup("cgTrackCutsV0select","cgTrackCutsV0select",AliDielectronCutGroup::kCompAND);
      cgTrackCutsV0select->AddCut(gammaV0Cuts);
      cgTrackCutsV0select->AddCut(trackCutsAOD);
      trackCuts = cgTrackCutsV0select;
      break;
      
    case kPbPb2011TRK_V0_2_loose:
      // primarily used for exclusion of V0s from track samples.
      AliDielectronV0Cuts *gammaV0Cuts = new AliDielectronV0Cuts("gammaV0Cuts","gammaV0Cuts");
      gammaV0Cuts->SetV0finder(AliDielectronV0Cuts::kOnTheFly);  // kAll(default), kOffline or kOnTheFly
      gammaV0Cuts->SetPdgCodes(22,11,11); // mother, daughter1 and 2
      gammaV0Cuts->AddCut(AliDielectronVarManager::kCosPointingAngle, TMath::Cos(0.05),   1.0,  kFALSE);
      gammaV0Cuts->AddCut(AliDielectronVarManager::kChi2NDF,                       0.0,  10.0,  kFALSE);
      gammaV0Cuts->AddCut(AliDielectronVarManager::kLegDist,                       0.0,   0.25, kFALSE);
      gammaV0Cuts->AddCut(AliDielectronVarManager::kR,                             3.0,  90.0,  kFALSE);
      gammaV0Cuts->AddCut(AliDielectronVarManager::kPsiPair,                       0.0,   0.20, kFALSE);
      gammaV0Cuts->AddCut(AliDielectronVarManager::kM,                             0.0,   0.10, kFALSE);
      gammaV0Cuts->AddCut(AliDielectronVarManager::kArmPt,                         0.0,   0.05, kFALSE);
      //      gammaV0Cuts->AddCut(AliDielectronVarManager::kArmAlpha,                     -0.35,  0.35, kFALSE);
      if (doExclusion==kExclude)  gammaV0Cuts->SetExcludeTracks(kTRUE);
      else                        gammaV0Cuts->SetExcludeTracks(kFALSE);
      // only default track cuts from AliDielectronV0Cuts::InitEvent() used here...
      cgTrackCutsV0excl = new AliDielectronCutGroup("cgTrackCutsV0excl","cgTrackCutsV0excl",AliDielectronCutGroup::kCompAND);
      cgTrackCutsV0excl->AddCut(gammaV0Cuts);
      trackCuts = cgTrackCutsV0excl;
      break;
      
    case kPbPb2011TRK_V0_3_Arm:
      // primarily meant for inclusion, for quite pure sample...
      AliDielectronV0Cuts *gammaV0Cuts = new AliDielectronV0Cuts("gammaV0Cuts","gammaV0Cuts");
      gammaV0Cuts->SetV0finder(AliDielectronV0Cuts::kOnTheFly);  // kAll(default), kOffline or kOnTheFly
      gammaV0Cuts->SetPdgCodes(22,11,11); // mother, daughter1 and 2
      gammaV0Cuts->AddCut(AliDielectronVarManager::kCosPointingAngle, TMath::Cos(0.02),   1.0,  kFALSE);
      gammaV0Cuts->AddCut(AliDielectronVarManager::kChi2NDF,                       0.0,  10.0,  kFALSE);
      gammaV0Cuts->AddCut(AliDielectronVarManager::kLegDist,                       0.0,   0.25, kFALSE);
      gammaV0Cuts->AddCut(AliDielectronVarManager::kR,                             3.0,  90.0,  kFALSE);
      gammaV0Cuts->AddCut(AliDielectronVarManager::kPsiPair,                       0.0,   0.05, kFALSE);
      gammaV0Cuts->AddCut(AliDielectronVarManager::kM,                             0.0,   0.05, kFALSE);
      gammaV0Cuts->AddCut(AliDielectronVarManager::kArmPt,                         0.0,   0.05, kFALSE);
      gammaV0Cuts->AddCut(AliDielectronVarManager::kArmAlpha,                     -0.35,  0.35, kFALSE); // should increase purity...
      if (doExclusion==kExclude)  gammaV0Cuts->SetExcludeTracks(kTRUE);
      else                        gammaV0Cuts->SetExcludeTracks(kFALSE);
      AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
      trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,     100.0, 160.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.8, 1.1);
      cgTrackCutsV0select = new AliDielectronCutGroup("cgTrackCutsV0select","cgTrackCutsV0select",AliDielectronCutGroup::kCompAND);
      cgTrackCutsV0select->AddCut(gammaV0Cuts);
      cgTrackCutsV0select->AddCut(trackCutsAOD);
      trackCuts = cgTrackCutsV0select;
      break;
      
      //----------
      // MC
      //----------
    case kPbPb2011MC_pi0Dal_1:
      AliDielectronVarCuts* trackCutsMC =new AliDielectronVarCuts("trackCutsMC","trackCutsMC");
      trackCutsMC->SetCutOnMCtruth(kTRUE);
      trackCutsMC->AddCut(AliDielectronVarManager::kPdgCode      , -11.01,  11.01);
      //trackCutsMC->AddCut(AliDielectronVarManager::kPdgCode      , -10.01,  10.01, kTRUE); //excludeRange
      trackCutsMC->AddCut(AliDielectronVarManager::kPdgCodeMother, 111);
      
      cgTrackCutsMC = new AliDielectronCutGroup("cgTrackCutsMC","cgTrackCutsMC",AliDielectronCutGroup::kCompAND);
      cgTrackCutsMC->AddCut(GetQualityCuts(kPbPb2011TRK_SPDfirst_1));
      cgTrackCutsMC->AddCut(trackCutsMC);
      trackCuts = cgTrackCutsMC;
      break;
      
      //----------
      // these MAIN settings have to combine different track selections:
      //----------
    case kPbPb2011_pidITS2gevTPCTOFif_trkSPDorSDD_5_tight:
    case kPbPb2011_pidITSTPCTOFif_trkSPDorSDD_2_loose:
    case kPbPb2011_pidTPCTOF_trkSPDorSDD_2_loose:
    case kPbPb2011_pidITSTPCTOFif_trkSPDorSDD_1:
    case kPbPb2011_pidTPCTOF_trkSPDorSDD_1:
      //----------
    case kPbPb2011TRK_SPDorSDD_1:
      // combine typical and new trackcuts with "kCompOR" condition:
      cgTrackCutsAna = new AliDielectronCutGroup("cgTrackCutsAnaSPDorSDD","cgTrackCutsAnaSPDorSDD",AliDielectronCutGroup::kCompOR);
      cgTrackCutsAna->AddCut(GetQualityCuts(kPbPb2011TRK_SPDfirst_1));         // typical trackcuts with requirement of SPD
      cgTrackCutsAna->AddCut(GetQualityCuts(kPbPb2011TRK_SDDfirstSPDnone_1)); // new additional trackcuts with SDD instead of SPD
      trackCuts = cgTrackCutsAna;
      break;
      
      //----------
      // these MAIN settings just load the main track selection directly below:
      //----------
    case kPbPb2011_pidITS2gevTPCTOFif_trkSPDfirst_5_tight:
    case kPbPb2011_pidITSTPC_trkSPDfirst_3:
    case kPbPb2011_pidTPC_trkSPDfirst_3:
    case kPbPb2011_pidITSTPCTOFif_trkSPDfirst_2_loose:
    case kPbPb2011_pidTPCTOF_trkSPDfirst_2_loose:
    case kPbPb2011_pidITSTPCTOFif_trkSPDfirst_1:
    case kPbPb2011_pidTPCTOF_trkSPDfirst_1:
      //----------
    case kPbPb2011TRK_SPDfirst_1: // main track selection
      AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
      trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,     4.0, 100.0); // means at least 2 with PID
      trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,     100.0, 160.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.8, 1.1); // lower limit 0.8 in most filterbits! // 1.1 since 26.02.2014
      AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
      trackCutsDiel->SetAODFilterBit(1<<4); // (=16) filterbit 4! //GetStandardITSTPCTrackCuts2011(kFALSE); loose DCA, 2D cut
      trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
      
      cgTrackCutsAnaSPDfirst = new AliDielectronCutGroup("cgTrackCutsAnaSPDfirst","cgTrackCutsAnaSPDfirst",AliDielectronCutGroup::kCompAND);
      cgTrackCutsAnaSPDfirst->AddCut(trackCutsDiel);
      cgTrackCutsAnaSPDfirst->AddCut(trackCutsAOD);
      if (fIsESDTask>0) cgTrackCutsAnaSPDfirst->AddCut(GetESDTrackCutsAna(kPbPb2011TRK_SPDfirst_1));
      trackCuts = cgTrackCutsAnaSPDfirst;
      break;
      
      
    case kPbPb2011TRK_SDDfirstSPDnone_1: // complimentary tracks, strictly without SPD, to be combined with others!
      AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
      trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,     3.0, 100.0); // means at least 3 with PID
      trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,     100.0, 160.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.8, 1.1); // lower limit 0.8 in most filterbits! // 1.1 since 26.02.2014
      AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
      trackCutsDiel->SetAODFilterBit(1<<6); // GetStandardITSTPCTrackCuts2011(kTRUE), SPD none, SDD first
      
      cgTrackCutsAnaSDDnoSPD = new AliDielectronCutGroup("cgTrackCutsAnaSDDnoSPD","cgTrackCutsAnaSDDnoSPD",AliDielectronCutGroup::kCompAND);
      cgTrackCutsAnaSDDnoSPD->AddCut(trackCutsDiel);
      cgTrackCutsAnaSDDnoSPD->AddCut(trackCutsAOD);
      if (fIsESDTask>0) cgTrackCutsAnaSDDnoSPD->AddCut(GetESDTrackCutsAna(kPbPb2011TRK_SDDfirstSPDnone_1));
      trackCuts = cgTrackCutsAnaSDDnoSPD;
      break;
      
      //----------
      // MAIN settings - combined trackset - variation 1
      //----------
    case kPbPb2011_pidITS2gevTPCTOFif_trkSPD5orSDD4cls_6_tight:
    case kPbPb2011_pidITSTPCTOFif_trkSPD5orSDD4cls_4:
      //----------
    case kPbPb2011TRK_SPDorSDD_2:
      // combine typical and new trackcuts with "kCompOR" condition:
      cgTrackCutsAna = new AliDielectronCutGroup("cgTrackCutsAnaSPDorSDD","cgTrackCutsAnaSPDorSDD",AliDielectronCutGroup::kCompOR);
      cgTrackCutsAna->AddCut(GetQualityCuts(kPbPb2011TRK_SPDfirst_2));         // typical trackcuts with requirement of SPD
      cgTrackCutsAna->AddCut(GetQualityCuts(kPbPb2011TRK_SDDfirstSPDnone_2)); // new additional trackcuts with SDD instead of SPD
      trackCuts = cgTrackCutsAna;
      break;
      
      //----------
      // MAIN settings - single trackset - variation 1
      //----------
    case kPbPb2011_pidITS2gevTPCTOFif_trkSPDfirst5cls_6_tight:
    case kPbPb2011_pidITSTPCTOFif_trkSPDfirst5cls_4:
      //----------
    case kPbPb2011TRK_SPDfirst_2: // main track selection, 5+ ITS clusters
      AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
      trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,     5.0, 100.0); // means at least 3 with PID
      trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,     100.0, 160.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.8, 1.1);
      AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
      trackCutsDiel->SetAODFilterBit(1<<4); // (=16) filterbit 4! //GetStandardITSTPCTrackCuts2011(kFALSE); loose DCA, 2D cut
      trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
      
      cgTrackCutsAnaSPDfirst = new AliDielectronCutGroup("cgTrackCutsAnaSPDfirst","cgTrackCutsAnaSPDfirst",AliDielectronCutGroup::kCompAND);
      cgTrackCutsAnaSPDfirst->AddCut(trackCutsDiel);
      cgTrackCutsAnaSPDfirst->AddCut(trackCutsAOD);
      if (fIsESDTask>0) cgTrackCutsAnaSPDfirst->AddCut(GetESDTrackCutsAna(kPbPb2011TRK_SPDfirst_2));
      trackCuts = cgTrackCutsAnaSPDfirst;
      break;
      
    case kPbPb2011TRK_SDDfirstSPDnone_2: // complimentary tracks, 4+ ITS clusters, strictly without SPD, to be combined with others!
      AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
      trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,     4.0, 100.0); // means at least 4 with PID
      trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,     100.0, 160.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.8, 1.1);
      AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
      trackCutsDiel->SetAODFilterBit(1<<6); // GetStandardITSTPCTrackCuts2011(kTRUE), SPD none, SDD first
      
      cgTrackCutsAnaSDDnoSPD = new AliDielectronCutGroup("cgTrackCutsAnaSDDnoSPD","cgTrackCutsAnaSDDnoSPD",AliDielectronCutGroup::kCompAND);
      cgTrackCutsAnaSDDnoSPD->AddCut(trackCutsDiel);
      cgTrackCutsAnaSDDnoSPD->AddCut(trackCutsAOD);
      if (fIsESDTask>0) cgTrackCutsAnaSDDnoSPD->AddCut(GetESDTrackCutsAna(kPbPb2011TRK_SDDfirstSPDnone_2));
      trackCuts = cgTrackCutsAnaSDDnoSPD;
      break;
      
      //----------
      // for Prefilter
      //----------
      
    case kPbPb2011TRK_FilterBit0: // TPC only + 3 ITS clusters
      AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
      trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,     3.0, 100.0);
      AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
      trackCutsDiel->SetAODFilterBit(1<<0);
      
      cgTrackCutsPre = new AliDielectronCutGroup("cgTrackCutsPre","cgTrackCutsPre",AliDielectronCutGroup::kCompAND);
      cgTrackCutsPre->AddCut(trackCutsDiel);
      cgTrackCutsPre->AddCut(trackCutsAOD);
      if (fIsESDTask>0) cgTrackCutsPre->AddCut(GetESDTrackCutsAna(kPbPb2011TRK_FilterBit0));
      trackCuts = cgTrackCutsPre;
      break;
    case kPbPb2011TRK_FilterBit1: // ITS standalone
      AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
      trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,     3.0, 100.0);
      AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
      trackCutsDiel->SetAODFilterBit(1<<1);
      
      cgTrackCutsPre = new AliDielectronCutGroup("cgTrackCutsPre","cgTrackCutsPre",AliDielectronCutGroup::kCompAND);
      cgTrackCutsPre->AddCut(trackCutsDiel);
      cgTrackCutsPre->AddCut(trackCutsAOD);
      if (fIsESDTask>0) cgTrackCutsPre->AddCut(GetESDTrackCutsAna(kPbPb2011TRK_FilterBit1));
      trackCuts = cgTrackCutsPre;
      break;
    case kPbPb2011TRK_FilterBit2: // similar to Theos prefilter cuts
      AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
      trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,     3.0, 100.0);
      AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
      trackCutsDiel->SetAODFilterBit(1<<2);
      trackCutsDiel->SetRequireITSRefit(kTRUE); // in addition to FilterBit
      //// FilterBit contains (from bit 0):
      //esdTrackCuts->SetMinNClustersTPC(50);           // different from Theo!
      //esdTrackCuts->SetMaxChi2PerClusterTPC(4);       // different from Theo!
      //// Theo:
      ////global
      //fesdTrackCuts->SetPtRange( 0.08 , 100. );       // different
      //fesdTrackCuts->SetEtaRange( -1.1 , 1.1 );       // different
      //fesdTrackCuts->SetAcceptKinkDaughters(kFALSE);  // ok
      //fesdTrackCuts->SetRequireSigmaToVertex(kFALSE); // not in AliDielectronTrackCuts (for AODs)
      //fesdTrackCuts->SetDCAToVertex2D(kFALSE);        // not in AliDielectronTrackCuts (for AODs)
      //fesdTrackCuts->SetMaxDCAToVertexZ(3.);          // ok
      //fesdTrackCuts->SetMaxDCAToVertexXY(1.);         // ok
      ////ITS cuts
      //fesdTrackCuts->SetRequireITSRefit(kTRUE);       // ok
      //fesdTrackCuts->SetMinNClustersITS(3);           // ok
      //fesdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny); // ok
      
      cgTrackCutsPre = new AliDielectronCutGroup("cgTrackCutsPre","cgTrackCutsPre",AliDielectronCutGroup::kCompAND);
      cgTrackCutsPre->AddCut(trackCutsDiel);
      cgTrackCutsPre->AddCut(trackCutsAOD);
      if (fIsESDTask>0) cgTrackCutsPre->AddCut(GetESDTrackCutsAna(kPbPb2011TRK_FilterBit2));
      trackCuts = cgTrackCutsPre;
      break;
      //----------
      
      //[...]
    default: cout << "No Analysis Track Cut defined " << endl;
  }
  return trackCuts;
} 


//*******************************************************************************
//*******************************************************************************
//** ESD TRACK CUTS TUNED FOR AGREEMENT BETWEEN AODS AND ESDS  ******************
//** NOT NECESSARILY 100% OPTIMIZED FOR DIEL-ANALYSIS          ******************
//*******************************************************************************
//*******************************************************************************

//WHEN RUNNING ON ESDs: LOAD Default Cuts for AODs
//_______________________________________________________________________________________________
AliAnalysisCuts* LMEECutLib::GetESDTrackCutsAna(Int_t cutSet) {
  //cout << " >>>>>>>>>>>>>>>>>>>>>> >>>>>>>>>>>>>>>>>>>>>> >>>>>>>>>>>>>>>>>>>>>> " << endl;
  cout << " >>>>>>>>>>>>>>>>>>>>>>  GetESDTrackCutsAna()  >>>>>>>>>>>>>>>>>>>>>> " << endl;
  //cout << " >>>>>>>>>>>>>>>>>>>>>> Setting ESD Track Cuts >>>>>>>>>>>>>>>>>>>>>> " << endl;
  //cout << " >>>>>>>>>>>>>>>>>>>>>> ( do we run on ESD?! ) >>>>>>>>>>>>>>>>>>>>>> " << endl;
  //cout << " >>>>>>>>>>>>>>>>>>>>>> >>>>>>>>>>>>>>>>>>>>>> >>>>>>>>>>>>>>>>>>>>>> " << endl;
  AliESDtrackCuts* esdTrackCuts = 0x0;
  
  if (fIsESDTask==kUnset) {
    cout << "WARNING: LMEECutLib::SetIsESDTask() was not set to true or false! Adding ESD track cuts now..." << endl;
  }
  
  
  //
  // from http://svn.cern.ch/guest/AliRoot/tags/v5-02-Rev-31/ANALYSIS/macros/AddTaskESDFilter.C
  //
  switch (cutSet) {
      
    case kPbPb2011TRK_FilterBit0: // like filter bit 0
      // Cuts on primary tracks
      esdTrackCuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
      break;
    case kPbPb2011TRK_FilterBit1: // like filter bit 1 // not sure if they work! may reject a track with itself!
      // ITS stand-alone tracks
      esdTrackCuts = new AliESDtrackCuts("ITS stand-alone Track Cuts", "ESD Track Cuts");
      esdTrackCuts->SetRequireITSStandAlone(kTRUE);
      break;
    case kPbPb2011TRK_FilterBit2: // is composed of bit 0 and 2 ("itsStrong->SetFilterMask(1);// AND with Standard track cuts")
      // Cuts on primary tracks
      esdTrackCuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
      // Pixel OR necessary for the electrons
      //esdTrackCuts = new AliESDtrackCuts("ITSorSPD", "pixel requirement for ITS");
      esdTrackCuts->SetTitle("standard TPC only + pixel requirement for ITS");
      esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
      break;
      
    case kPbPb2011TRK_FilterBit4:
    case kPbPb2011TRK_SPDfirst_2:
    case kPbPb2011TRK_SPDfirst_1: // like filter bit 4 (Int: 16), AOD095&115
      // standard cuts with very loose DCA
      esdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE); 
      esdTrackCuts->SetMaxDCAToVertexXY(2.4);
      esdTrackCuts->SetMaxDCAToVertexZ(3.2);
      esdTrackCuts->SetDCAToVertex2D(kTRUE);
      
      //The cuts below should be the onyl ones that are missing
      //explicitely in the TrackCutsAna method
      //To be sure, StandardITSTPCTrackCuts is loaded however
      /* 
       esdTrackCuts = new AliESDtrackCuts();
       esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
       //Not done so far via dielectron cuts:
       */
      /*
       esdTrackCuts->SetDCAToVertex2D(kFALSE);
       esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
       esdTrackCuts->SetMaxChi2PerClusterITS(36);
       */
      break;
      
    case kPbPb2011TRK_FilterBit6:
    case kPbPb2011TRK_SDDfirstSPDnone_2:
    case kPbPb2011TRK_SDDfirstSPDnone_1: // like filter bit 6
      // standard cuts with tight DCA but with requiring the first SDD cluster instead of an SPD cluster
      // tracks selected by this cut are exclusive to those selected by the previous cut
      esdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kTRUE); 
      esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kNone);
      esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSDD, AliESDtrackCuts::kFirst);
      break;
      
    default:
      cout << "No ESD Track Cut defined " << endl;
      cout << "Applying cuts of filter bit 4 " << endl;
      return GetESDTrackCutsAna(kPbPb2011TRK_SPDfirst_1);
  }
  return esdTrackCuts;
}



//_______________________________________________________________________________________________
AliAnalysisCuts* LMEECutLib::GetMCTrackCuts() {
  //cout << " >>>>>>>>>>>>>>>>>>>>>> >>>>>>>>>>>>>>>>>>>>>> >>>>>>>>>>>>>>>>>>>>>> " << endl;
  cout << " >>>>>>>>>>>>>>>>>>>>>>  GetMCTrackCuts()  >>>>>>>>>>>>>>>>>>>>>> " << endl;
  //cout << " >>>>>>>>>>>>>>>>>>>>>> >>>>>>>>>>>>>>>>>>>>>> >>>>>>>>>>>>>>>>>>>>>> " << endl;
  AliAnalysisCuts* trackCuts=0x0;
  switch (selectedQualityAna) {
    default:
      AliDielectronVarCuts* trackCutsMC =new AliDielectronVarCuts("trackCutsMC","trackCutsMC");
      trackCutsMC->SetCutOnMCtruth(kTRUE);
      trackCutsMC->AddCut(AliDielectronVarManager::kPdgCode      , -11.01,  11.01);
      //trackCutsMC->AddCut(AliDielectronVarManager::kPdgCode      , -10.01,  10.01, kTRUE); //excludeRange
      //trackCutsMC->AddCut(AliDielectronVarManager::kPdgCodeMother, 111);
      trackCuts = trackCutsMC;
      break;
  }
  return trackCuts;
}



//_______________________________________________________________________________________________
void LMEECutLib::AddMCSignals(TNamed* task, Int_t cutDefinition) {
  //Do we have an MC handler?
  //if (!die->GetHasMC()) return;
  cout << " >>>>>>>>>>>>>>>>>>>>>>  AddMCSignals()  >>>>>>>>>>>>>>>>>>>>>> " << endl;
  
  AliDielectronSignalMC* eleFinalState = new AliDielectronSignalMC("eleFinalState","eleFinalState");
  eleFinalState->SetFillPureMCStep(fFillPureMC);
  eleFinalState->SetLegPDGs(11,1);  //dummy second leg (never MCtrue)
  eleFinalState->SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleFinalState->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState); // for leg electrons: kFinalState = kPrimary
  //mother
  eleFinalState->SetMotherSources(AliDielectronSignalMC::kPrimary, AliDielectronSignalMC::kPrimary); // has no effect for this leg setting (kFinalState). Still includes injected J/psi.
  
  AliDielectronSignalMC* eleDirect = new AliDielectronSignalMC("eleDirect","eleDirect");
  eleDirect->SetFillPureMCStep(fFillPureMC);
  eleDirect->SetLegPDGs(11,1);  //dummy second leg (never MCtrue)
  eleDirect->SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleDirect->SetLegSources(AliDielectronSignalMC::kDirect, AliDielectronSignalMC::kDirect); // for leg electrons: kDirect is empty
  
  AliDielectronSignalMC* eleNoCocktail = new AliDielectronSignalMC("eleNoCocktail","eleNoCocktail");
  eleNoCocktail->SetFillPureMCStep(fFillPureMC);
  eleNoCocktail->SetLegPDGs(11,1);  //dummy second leg (never MCtrue)
  eleNoCocktail->SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleNoCocktail->SetLegSources(AliDielectronSignalMC::kNoCocktail, AliDielectronSignalMC::kNoCocktail); // for leg electrons: kNoCocktail = kFinalState + kSecondary
  
  AliDielectronSignalMC* eleSecondary = new AliDielectronSignalMC("eleSecondary","eleSecondary");
  eleSecondary->SetFillPureMCStep(fFillPureMC);
  eleSecondary->SetLegPDGs(11,1);  //dummy second leg (never MCtrue)
  eleSecondary->SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleSecondary->SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);

  AliDielectronSignalMC* eleSecondaryWeak = new AliDielectronSignalMC("eleSecondaryWeak","eleSecondaryWeak");
  eleSecondaryWeak->SetFillPureMCStep(fFillPureMC);
  eleSecondaryWeak->SetLegPDGs(11,1);  //dummy second leg (never MCtrue)
  eleSecondaryWeak->SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleSecondaryWeak->SetLegSources(AliDielectronSignalMC::kSecondaryFromWeakDecay, AliDielectronSignalMC::kSecondaryFromWeakDecay); // for leg electrons: kSecondaryFromWeakDecay is empty
  
  
  AliDielectronSignalMC* eleFromBGEvent = new AliDielectronSignalMC("eleFromBGEvent","eleFromBGEvent");
  eleFromBGEvent->SetFillPureMCStep(fFillPureMC);
  eleFromBGEvent->SetLegPDGs(11,1);  //dummy second leg (never MCtrue)
  eleFromBGEvent->SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleFromBGEvent->SetLegSources(AliDielectronSignalMC::kFromBGEvent, AliDielectronSignalMC::kFromBGEvent);

  AliDielectronSignalMC* eleFinalStateFromBGEvent = new AliDielectronSignalMC("eleFinalStateFromBGEvent","eleFinalStateFromBGEvent");
  eleFinalStateFromBGEvent->SetFillPureMCStep(fFillPureMC);
  eleFinalStateFromBGEvent->SetLegPDGs(11,1);  //dummy second leg (never MCtrue)
  eleFinalStateFromBGEvent->SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleFinalStateFromBGEvent->SetLegSources(AliDielectronSignalMC::kFinalStateFromBGEvent, AliDielectronSignalMC::kFinalStateFromBGEvent);
  
  AliDielectronSignalMC* eleFinalStateFromBGEvent_mDir = new AliDielectronSignalMC("eleFinalStateFromBGEvent_mDir","eleFinalStateFromBGEvent_mDir");
  eleFinalStateFromBGEvent_mDir->SetFillPureMCStep(fFillPureMC);
  eleFinalStateFromBGEvent_mDir->SetLegPDGs(11,1);  //dummy second leg (never MCtrue)
  eleFinalStateFromBGEvent_mDir->SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleFinalStateFromBGEvent_mDir->SetLegSources(AliDielectronSignalMC::kFinalStateFromBGEvent, AliDielectronSignalMC::kFinalStateFromBGEvent);
  //mother
  eleFinalStateFromBGEvent_mDir->SetMotherSources(AliDielectronSignalMC::kDirect, AliDielectronSignalMC::kDirect);
  
  
  // In Pythia the mother label of primary hadrons gives the (positive) label of the mother quark,
  //  so mother source kNoCocktail should correspond to non-injected signals.
  // In Hijing however the mother label of primary hadrons is negative (they are not traced back to the quarks),
  //  so mother source kNoCocktail gives hadrons which have a mother, i.e. electrons whith have a grandmother.
  AliDielectronSignalMC* eleWithGrandmother = new AliDielectronSignalMC("eleWithGrandmother","eleWithGrandmother");
  eleWithGrandmother->SetFillPureMCStep(fFillPureMC);
  eleWithGrandmother->SetLegPDGs(11,1);  //dummy second leg (never MCtrue)
  eleWithGrandmother->SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleWithGrandmother->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState); // for leg electrons: kFinalState = kPrimary, also in case of 'SetMotherSources(AliDielectronSignalMC::kNoCocktail, AliDielectronSignalMC::kNoCocktail)'.
  //mother
  eleWithGrandmother->SetMotherSources(AliDielectronSignalMC::kNoCocktail, AliDielectronSignalMC::kNoCocktail); // excludes JPsi, includes conversions.
  eleWithGrandmother->SetMotherPDGs(22,22,kTRUE,kTRUE); // exclude conversion electrons
  
  
  // stuff:
  //mother
  //  ele->SetCheckBothChargesMothers(kTRUE,kTRUE);
  //  ele->SetMotherPDGs(902,902,kTRUE,kTRUE); // exclude open charm,beauty hadrons
  //  ele->SetCheckBothChargesGrandMothers(kTRUE,kTRUE);
  //  ele->SetGrandMotherPDGs(902,902,kTRUE,kTRUE); // exclude open charm,beauty hadrons
  //  ele->SetGrandMotherPDGs(500,500,kTRUE,kTRUE); // exclude non-prompt jpsi eletrons
  
  
  AliDielectronSignalMC* gammaConversion = new AliDielectronSignalMC("gammaConversion","gammaConversion");
  gammaConversion->SetFillPureMCStep(fFillPureMC);
  gammaConversion->SetLegPDGs(11,-11);
  gammaConversion->SetCheckBothChargesLegs(kTRUE,kTRUE);
  gammaConversion->SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
  //mother
  gammaConversion->SetMotherPDGs(22,22);
  gammaConversion->SetMothersRelation(AliDielectronSignalMC::kSame);
  
//  // add direct di lepton resonances
//  AliDielectronSignalMC* directP[7];
//  TParticlePDG *ap;
//  Int_t pdg[] = {111, 113, 221, 223, 331, 333, 443};
//  for(Int_t i=0; i<7; i++) {
//    ap = TDatabasePDG::Instance()->GetParticle(pdg[i]);
//    directP[i] = new AliDielectronSignalMC(Form("direct%s",ap->GetName()),Form("direct%s",ap->GetName()));
//    directP[i]->SetLegPDGs(11,-11);
//    directP[i]->SetMotherPDGs(pdg[i],pdg[i]);
//    directP[i]->SetMothersRelation(AliDielectronSignalMC::kSame);
//    directP[i]->SetFillPureMCStep(fFillPureMC);
//    directP[i]->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
//    directP[i]->SetMotherSources(AliDielectronSignalMC::kDirect, AliDielectronSignalMC::kDirect);
//    // directP[i]->SetCheckBothChargesLegs(kTRUE,kTRUE);
//    // directP[i]->SetCheckBothChargesMothers(kTRUE,kTRUE);
//  }
  
  AliDielectronSignalMC* pair_sameMother = new AliDielectronSignalMC("MCpair_sameMother","MCpairSameMother");
  pair_sameMother->SetFillPureMCStep(fFillPureMC);
  pair_sameMother->SetLegPDGs(11,-11);
  pair_sameMother->SetCheckBothChargesLegs(kTRUE,kTRUE);
  pair_sameMother->SetLegSources(AliDielectronSignalMC::kFinalStateFromBGEvent, AliDielectronSignalMC::kFinalStateFromBGEvent);
  //mother
  pair_sameMother->SetMothersRelation(AliDielectronSignalMC::kSame);
  pair_sameMother->SetMotherPDGs(22,22,kTRUE,kTRUE); // exclude conversion electrons
  
  AliDielectronSignalMC* pair_diffMother = new AliDielectronSignalMC("MCpair_diffMother","MCpairDiffMother");
  pair_diffMother->SetFillPureMCStep(fFillPureMC);
  pair_diffMother->SetLegPDGs(11,-11);
  pair_diffMother->SetCheckBothChargesLegs(kTRUE,kTRUE);
  pair_diffMother->SetLegSources(AliDielectronSignalMC::kFinalStateFromBGEvent, AliDielectronSignalMC::kFinalStateFromBGEvent);
  //mother
  pair_diffMother->SetMothersRelation(AliDielectronSignalMC::kDifferent);
  pair_diffMother->SetMotherPDGs(22,22,kTRUE,kTRUE); // exclude conversion electrons
  
  
  if (task->IsA()==AliAnalysisTaskElectronEfficiency::Class()) {
    // selection. only the first signal will be used.
//    (static_cast<AliAnalysisTaskElectronEfficiency*>task)->AddSignalMC(eleFinalStateFromBGEvent);
    (static_cast<AliAnalysisTaskElectronEfficiency*>task)->AddSignalMC(eleFinalStateFromBGEvent_mDir);
  }
  else if (task->IsA()==AliDielectron::Class()) {
    if (! (static_cast<AliDielectron*>task)->GetHasMC()) return;
    switch (cutDefinition) {
      case -1:
        // selection. multiple signals are possible.
        (static_cast<AliDielectron*>task)->AddSignalMC(eleFromBGEvent);
        (static_cast<AliDielectron*>task)->AddSignalMC(eleFinalStateFromBGEvent);
        (static_cast<AliDielectron*>task)->AddSignalMC(eleFinalStateFromBGEvent_mDir);
        (static_cast<AliDielectron*>task)->AddSignalMC(eleWithGrandmother);
        //die->AddSignalMC(elePrimary);
        //die->AddSignalMC(eleFinalState);    // for leg electrons: kFinalState = kPrimary
        //die->AddSignalMC(eleDirect);        // for leg electrons: kDirect is empty
        //die->AddSignalMC(eleNoCocktail);    // for leg electrons: kNoCocktail = kFinalState + kSecondary
        //die->AddSignalMC(eleSecondary);     // for leg electrons: kSecondary = number of 'gammaConversion'
        //die->AddSignalMC(eleSecondaryWeak); // for leg electrons: kSecondaryFromWeakDecay is empty
        //die->AddSignalMC(gammaConversion); //
        break;
        
//      case 1:
//        for(Int_t i=0; i<7; i++) (static_cast<AliDielectron*>task)->AddSignalMC(directP[i]);
//        break;
        
      case 2:
        (static_cast<AliDielectron*>task)->AddSignalMC(pair_sameMother);
        (static_cast<AliDielectron*>task)->AddSignalMC(pair_diffMother);
        break;
        
      default:
        break;
    }
    for (Int_t i=0; i<(static_cast<AliDielectron*>task)->GetMCSignals()->GetEntriesFast(); ++i) {
      cout << "added MCsignal: " << (static_cast<AliDielectron*>task)->GetMCSignals()->At(i)->GetName() << endl;
    }
  }
}


//_______________________________________________________________________________________________
void LMEECutLib::InitHistograms(AliDielectron *die, Int_t cutDefinition) {
  //
  // Initialise the histograms
  //
  
  //Setup histogram Manager
  AliDielectronHistos *histos=new AliDielectronHistos(die->GetName(),die->GetTitle());
  
  //Initialise histogram classes
  histos->SetReservedWords("Track;Pair;Track_Legs;Pre;RejPair;RejTrack;Random");
  
  //Event class
  histos->AddClass("Event"); // all classes will be stored in 'THashList fHistoList'
  
  //Track classes
  //to fill also track info from 2nd event loop until 3
  // in AliDielectron.cxx: fgkTrackClassNames[4] = {"ev1+","ev1-","ev2+","ev2-"};
  if (die->DoEventProcess()) {
    for (Int_t i=0; i<2; ++i){
      histos->AddClass(Form("Track_%s",AliDielectron::TrackClassName(i)));
    }
  }
  
  //Pair classes
  // to fill also mixed event histograms loop until 10
  // fgkPairClassNames[11] = {
  //  "ev1+_ev1+",  "ev1+_ev1-",  "ev1-_ev1-",  // 0-2 (same event)
  //  "ev1+_ev2+",  "ev1-_ev2+",  "ev2+_ev2+",  // 3-4 (+5)
  //  "ev1+_ev2-",  "ev1-_ev2-",                // 6-7
  //  "ev2+_ev2-",  "ev2-_ev2-",  "ev1+_ev1-_TR"
  // };
  
  if (!fIsQATask && !fIsRandomRejTask) // -> analysis with pairing
  {
    for (Int_t i=0; i<3; ++i){
      histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(i)));
      // Legs of final Pairs. Both charges together. No duplicate entries.
      if (!fIsQATask) {
        histos->AddClass(Form("Track_Legs_%s",AliDielectron::PairClassName(i))); // not TrackClassName, see 'AliDielectron::FillHistograms(...)'
      }
    }
    
    //add MC signal histograms to pair class
    if(die->GetMCSignals()) {
      for (Int_t i=0; i<die->GetMCSignals()->GetEntriesFast(); ++i) {
        //histos->AddClass(Form("Track_Legs_%s",die->GetMCSignals()->At(i)->GetName()));
        histos->AddClass(Form("Pair_%s",die->GetMCSignals()->At(i)->GetName()));
        if (fFillPureMC) histos->AddClass(Form("Pair_%s_MCtruth",die->GetMCSignals()->At(i)->GetName()));
      }
    }
    
    //Mixed event and track rotation
    if (die->GetMixingHandler()) {
      histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(3)));
      histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(4)));
      histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(6)));
      histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(7)));
    }
    if (die->GetTrackRotator()) {
      histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(10)));
    }
    
    if (fDoRejectionStep) 
    {
      //PreFilter Classes
      //to fill also track info from 2nd event loop until 2
//      for (Int_t i=0; i<2; ++i){
//        histos->AddClass(Form("Pre_%s",AliDielectron::TrackClassName(i)));
//        // class 'Pre_%s': "Fill Histogram information for tracks after prefilter"(AliDielectron::FillHistogramsTracks(...)), so it is identical to Track_%s if no paircut is used.
//      }
      
      //Create Classes for Rejected Tracks/Pairs:
      for (Int_t i=0; i<3; ++i){
        histos->AddClass(Form("RejPair_%s",AliDielectron::PairClassName(i)));
        // Legs of rejected Pairs. Both charges together. One track can and will make multiple entries.
        histos->AddClass(Form("RejTrack_%s",AliDielectron::PairClassName(i))); // not TrackClassName, see 'AliDielectron::FillHistogramsPair(...)'
      }
      
      /*
       //track rotation
       histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(AliDielectron::kEv1PMRot)));
       histos->AddClass(Form("Track_Legs_%s",AliDielectron::PairClassName(AliDielectron::kEv1PMRot)));
       */
    }// end: (fDoRejectionStep)
    
  }// end: (!fIsQATask)
  
  
  if (fIsRandomRejTask)
  {
    //
    // _____ histograms for AliAnalysisTaskRandomRejection _____
    //
    histos->AddClass("Rand_Pair");
    histos->AddClass("Rand_RejPair");
    const char* cRandomPairClassNames[4] = { "Testpart", "DataEle", "RejTestpart", "RejDataEle" };
    for (Int_t i=0; i<4; ++i){
      histos->AddClass(Form("Random_%s",cRandomPairClassNames[i]));
    }
    histos->UserHistogram("Random","Pt","",200,0,10.,AliDielectronVarManager::kPt);
    histos->UserHistogram("Random","Eta","",200,-2,2,AliDielectronVarManager::kEta);
    histos->UserHistogram("Random","Phi","",120,0.,TMath::TwoPi(),AliDielectronVarManager::kPhi);
    histos->UserHistogram("Random","Px","",200,0,10.,AliDielectronVarManager::kPx);
    histos->UserHistogram("Random","Py","",200,0,10.,AliDielectronVarManager::kPy);
    histos->UserHistogram("Random","Pz","",200,0,10.,AliDielectronVarManager::kPz);
    histos->UserHistogram("Random","Pt_Eta_Phi","",GetVector(kP2D),GetVector(kEta2D),GetVector(kPhi2D), AliDielectronVarManager::kPt,AliDielectronVarManager::kEta,AliDielectronVarManager::kPhi);
    histos->UserHistogram("Random","dXY_dZ","",100,-1.,1.,150,-3.,3.,AliDielectronVarManager::kImpactParXY,AliDielectronVarManager::kImpactParZ);
    histos->UserHistogram("Random","TPC_dEdx_Pt","",
                          100,0.,2.,100,0.,100., AliDielectronVarManager::kPt,AliDielectronVarManager::kTPCsignal);
  }
  
	//add histograms to event class
	histos->UserHistogram("Event","nEvents","",1,0.,1.,AliDielectronVarManager::kNevents);
  histos->UserHistogram("Event","RunNumber","",AliDielectronHelper::MakeArbitraryBinning("167900, 167987, 167988, 168310, 168311, 168322, 168325, 168341, 168342, 168361, 168362, 168458, 168460, 168464, 168467, 168511, 168512, 168514, 168777, 168826, 168988, 168992, 169035, 169040, 169044, 169045, 169091, 169094, 169099, 169138, 169144, 169145, 169148, 169156, 169160, 169167, 169238, 169411, 169415, 169417, 169418, 169419, 169420, 169475, 169498, 169504, 169506, 169512, 169515, 169550, 169553, 169554, 169555, 169557, 169586, 169587, 169588, 169590, 169591, 169835, 169837, 169838, 169846, 169855, 169858, 169859, 170027, 170040, 170081, 170083, 170084, 170085, 170088, 170089, 170091, 170155, 170159, 170163, 170193, 170203, 170204, 170207, 170228, 170230, 170268, 170269, 170270, 170306, 170308, 170309, 170311, 170312, 170313, 170315, 170387, 170388, 170572, 170593, 170600"),AliDielectronVarManager::kRunNumber);
	histos->UserHistogram("Event","Centrality","","-1,0,10,20,30,40,50,60,70,80,90,100,101;#events",AliDielectronVarManager::kCentrality);
  histos->UserHistogram("Event","centrality","",100,0,100,AliDielectronVarManager::kCentrality);
  histos->UserHistogram("Event","nESDTracks","",8000,0,80000,AliDielectronVarManager::kNTrk);
  //histos->UserHistogram("Event","Nacc","",8000,0,8000,AliDielectronVarManager::kNacc);
  histos->UserHistogram("Event","RefMultTPConly","",8000,0,8000,AliDielectronVarManager::kRefMultTPConly);
  histos->UserHistogram("Event","epTPC","",240,-TMath::Pi(),TMath::Pi(),AliDielectronVarManager::kTPCrpH2uc);
  histos->UserHistogram("Event","epV0AC","",240,-TMath::Pi(),TMath::Pi(),AliDielectronVarManager::kv0ACrpH2);
  histos->UserHistogram("Event","epV0AC_epTPC","",120,-TMath::PiOver2(),TMath::PiOver2(),120,-TMath::PiOver2(),TMath::PiOver2(),AliDielectronVarManager::kTPCrpH2uc,AliDielectronVarManager::kv0ACrpH2);
  
  
  //add histograms to Track classes
  // axis labels are set to the values in 'AliDielectronVarManager.cxx' if the histogram title is empty or starts with ';'! [see AliDielectronHistos::AdaptNameTitle(...)]
  histos->UserHistogram("Track","Pt",";Pt [GeV];#tracks",200,0,10.,AliDielectronVarManager::kPt);
  histos->UserHistogram("Track","Px",";Px [GeV];#tracks",200,0,10.,AliDielectronVarManager::kPx);
  histos->UserHistogram("Track","Py",";Py [GeV];#tracks",200,0,10.,AliDielectronVarManager::kPy);
  histos->UserHistogram("Track","Pz",";Pz [GeV];#tracks",200,0,10.,AliDielectronVarManager::kPz);
  histos->UserHistogram("Track","P_PIn",";p (GeV/c);p_{in} (GeV/c)",
                        GetVector(kP2D), GetVector(kP2D), AliDielectronVarManager::kP,AliDielectronVarManager::kPIn);
  
  // ITS
  histos->UserHistogram("Track","ITS_dEdx_P",";p (GeV/c);ITS signal (arb units)",
                        GetVector(kP2D), BinsToVector(700,0.,700.), AliDielectronVarManager::kP,AliDielectronVarManager::kITSsignal);
  histos->UserHistogram("Track","ITSnSigmaEle_P",";p (GeV/c);n#sigma_{ele}^{ITS}",
                        GetVector(kP2D), GetVector(kSigmaEle), AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaEle);
  histos->UserHistogram("Track","ITSnSigmaEle_Eta",";Eta;n#sigma_{ele}^{ITS}",
                        GetVector(kEta2D), GetVector(kSigmaEle), AliDielectronVarManager::kEta,AliDielectronVarManager::kITSnSigmaEle);
  // 3D
  if (fIsQATask) {
    histos->UserHistogram("Track","ITSnSigmaEle_Eta_P",";Eta;n#sigma_{ele}^{ITS};p (GeV/c)",
                          GetVector(kEta3D), GetVector(kSigmaEle), GetVector(kP2D),
                          AliDielectronVarManager::kEta,AliDielectronVarManager::kITSnSigmaEle,AliDielectronVarManager::kP);
    histos->UserHistogram("Track","ITSnSigmaEle_Eta_RefMultTPConly",";Eta;n#sigma_{ele}^{ITS};N_{TPC ref}",
                          GetVector(kEta3D), GetVector(kSigmaEle), BinsToVector(100,0.,5000.),
                          AliDielectronVarManager::kEta,AliDielectronVarManager::kITSnSigmaEle,AliDielectronVarManager::kRefMultTPConly);
    //histos->UserHistogram("Track","ITSnSigmaPio_P",";p (GeV/c);n#sigma_{pion}^{ITS}",
    //                      GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaPio);
    //histos->UserHistogram("Track","ITSnSigmaKao_P",";p (GeV/c);n#sigma_{kaon}^{ITS}",
    //                      GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaKao);
    //histos->UserHistogram("Track","ITSnSigmaPro_P",";p (GeV/c);n#sigma_{proton}^{ITS}",
    //                      GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaPro);
  }
  // TPC
  histos->UserHistogram("Track","TPC_dEdx_P",";p_{in} (GeV/c);TPC signal (arb units)",
                        GetVector(kP2D), GetVector(kTPCdEdx), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCsignal);
  histos->UserHistogram("Track","TPCnSigmaEle_P",";p_{in} (GeV/c);n#sigma_{ele}^{TPC}",
                        GetVector(kP2D), GetVector(kSigmaEle), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle);
  histos->UserHistogram("Track","TPCnSigmaEle_Eta",";Eta;n#sigma_{ele}^{TPC}",
                        GetVector(kEta2D), GetVector(kSigmaEle), AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEle);
  if (die->DoEventProcess()) {
    // some histograms with SigmaEleRaw to benchmark the dEdx-eta-correction:
    histos->UserHistogram("Track","TPCnSigmaEleRaw_Eta",";Eta;n#sigma_{ele,Raw}^{TPC}",
                          GetVector(kEta2D), GetVector(kSigmaEle), AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEleRaw);
    histos->UserHistogram("Track","TPCnSigmaEle_RefMultTPConly",";N_{TPC ref}; n#sigma_{ele}^{TPC}",
                          BinsToVector(100,0.,5000.), GetVector(kSigmaEle), AliDielectronVarManager::kRefMultTPConly,AliDielectronVarManager::kTPCnSigmaEle);
    histos->UserHistogram("Track","TPCnSigmaEleRaw_RefMultTPConly",";N_{TPC ref}; n#sigma_{ele,Raw}^{TPC}",
                          BinsToVector(100,0.,5000.), GetVector(kSigmaEle), AliDielectronVarManager::kRefMultTPConly,AliDielectronVarManager::kTPCnSigmaEleRaw);
//    histos->UserHistogram("Track","TPCnSigmaEle_Nacc",";N_{acc}; n#sigma_{ele}^{TPC}",
//                          BinsToVector(100,0.,5000.), GetVector(kSigmaEle), AliDielectronVarManager::kNacc,AliDielectronVarManager::kTPCnSigmaEle);
    histos->UserHistogram("Track","TPCnSigmaEle_RunNumber",";run;n#sigma_{ele}^{TPC}",
                          GetVector(kRuns), GetVector(kSigmaEle), AliDielectronVarManager::kRunNumber,AliDielectronVarManager::kTPCnSigmaEle);
  }
  
  if (fIsQATask) {
    histos->UserHistogram("Track","TPC_dEdx_Eta",";Eta;TPC signal (arb units)",
                          GetVector(kEta2D), GetVector(kTPCdEdx), AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCsignal);
    histos->UserHistogram("Track","TPC_dEdx_Eta_P",";Eta;TPC signal (arb units);p_{in} (GeV/c)",
                          GetVector(kEta3D), GetVector(kTPCdEdx), GetVector(kP2D),
                          AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCsignal,AliDielectronVarManager::kPIn);
    histos->UserHistogram("Track","TPCnSigmaEle_P_dEdx",";p_{in} (GeV/c);n#sigma_{ele}^{TPC};TPC signal (arb units)",
                          GetVector(kP2D), GetVector(kSigmaEle), GetVector(kTPCdEdx),
                          AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kTPCsignal);
    
    // 3D
    histos->UserHistogram("Track","TPCnSigmaEle_Eta_P",";Eta;n#sigma_{ele}^{TPC};p_{in} (GeV/c)",
                          GetVector(kEta3D), GetVector(kSigmaEle), GetVector(kP2D),
                          AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kPIn);
    histos->UserHistogram("Track","TPCnSigmaEle_Eta_RefMultTPConly",";Eta;n#sigma_{ele}^{TPC};N_{TPC ref}",
                          GetVector(kEta3D), GetVector(kSigmaEle), BinsToVector(100,0.,5000.),
                          AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kRefMultTPConly);
    histos->UserHistogram("Track","TPCnSigmaEle_Eta_RunNumber",";Eta;n#sigma_{ele}^{TPC};run",
                          GetVector(kEta3D), GetVector(kSigmaEle), GetVector(kRuns),
                          AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kRunNumber);
    histos->UserHistogram("Track","TPCnSigmaEle_RefMultTPConly_RunNumber",";N_{TPC ref};n#sigma_{ele}^{TPC};run",
                          BinsToVector(100,0.,5000.), GetVector(kSigmaEle), GetVector(kRuns),
                          AliDielectronVarManager::kRefMultTPConly,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kRunNumber);
    histos->UserHistogram("Track","TPCnSigmaEle_P_RunNumber",";p_{in} (GeV/c);n#sigma_{ele}^{TPC};run",
                          GetVector(kP2D), GetVector(kSigmaEle), GetVector(kRuns), 
                          AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kRunNumber);
    histos->UserHistogram("Track","TPCnSigmaEle_Eta_Nacc",";Eta;n#sigma_{ele}^{TPC};N_{acc}",
                          GetVector(kEta3D), GetVector(kSigmaEle), BinsToVector(100,0.,5000.),
                          AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kNacc);
    histos->UserHistogram("Track","RefMultTPConly_Nacc",";N_{TPC ref};N_{acc}",
                          BinsToVector(100,0.,5000.), BinsToVector(100,0.,5000.), AliDielectronVarManager::kRefMultTPConly,AliDielectronVarManager::kNacc);
    histos->UserHistogram("Track","TPCnSigmaEleRaw_Eta_RefMultTPConly",";Eta;n#sigma_{ele,Raw}^{TPC};N_{TPC ref}",
                          GetVector(kEta3D), GetVector(kSigmaEle), BinsToVector(100,0.,5000.),
                          AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEleRaw,AliDielectronVarManager::kRefMultTPConly);
    
    histos->UserHistogram("Track","TPCnSigmaPio_P",";p_{in} (GeV/c);n#sigma_{pion}^{TPC}",
                          GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaPio);
    histos->UserHistogram("Track","TPCnSigmaKao_P",";p_{in} (GeV/c);n#sigma_{kaon}^{TPC}",
                          GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaKao);
    histos->UserHistogram("Track","TPCnSigmaPro_P",";p_{in} (GeV/c);n#sigma_{proton}^{TPC}",
                          GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaPro);
    histos->UserHistogram("Track","TPCnSigmaKao_Eta",";Eta;n#sigma_{kaon}^{TPC}",
                          GetVector(kEta2D), GetVector(kSigmaOther), AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaKao);
  }
  
  // TRD
  // TRD variables need lot of computing time. since a performance update by Julian, they will not be computed if not needed! (status 7.3.14)
  //  histos->UserHistogram("Track","TRDpidPobEle_P","TRD PID probability Electrons;p_{in} (GeV/c);TRD prob Electrons",
  //                        GetVector(kP2D), 100,0.,1.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTRDprobEle);
  //  histos->UserHistogram("Track","TRDpidPobPio_P","TRD PID probability Pions;p_{in} (GeV/c);TRD prob Pions",
  //                        GetVector(kP2D), 100,0.,1.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTRDprobPio);
  
  // TOF
  histos->UserHistogram("Track","TOFbeta_P",";p_{in} (GeV/c);TOF beta",
                        GetVector(kP2D), BinsToVector(120,0.,1.2), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFbeta);
  histos->UserHistogram("Track","TOFnSigmaEle_P",";p_{in} (GeV/c);TOF number of sigmas Electrons",
                        GetVector(kP2D), GetVector(kSigmaEle), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaEle);
  if (fIsQATask) {
    histos->UserHistogram("Track","TOFnSigmaPio_P",";p_{in} (GeV/c);TOF number of sigmas Pions",
                          GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaPio);
    histos->UserHistogram("Track","TOFnSigmaKao_P",";p_{in} (GeV/c);TOF number of sigmas Kaons",
                          GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaKao);
    histos->UserHistogram("Track","TOFnSigmaPro_P",";p_{in} (GeV/c);TOF number of sigmas Protons",
                          GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaPro);
  }    
  // 2D-PID
  if (fIsQATask) {
    histos->UserHistogram("Track","PIn_TPCnSigmaEle_ITSnSigmaEle",";p_{in} (GeV/c);n#sigma_{ele}^{TPC};n#sigma_{ele}^{ITS}",
                          GetVector(kP2D), GetVector(kSigmaEle), GetVector(kSigmaEle),
                          AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kITSnSigmaEle);
    histos->UserHistogram("Track","PIn_TPCnSigmaEle_TOFnSigmaEle",";p_{in} (GeV/c);n#sigma_{ele}^{TPC};TOF number of sigmas Electrons",
                          GetVector(kP2D), GetVector(kSigmaEle), BinsToVector(50,-5.,5.),
                          AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kTOFnSigmaEle);
    // for Anisas contamination & purity tool
    histos->UserHistogram("Track","PIn_TPCnSigmaEle_TPCnSigmaPio",";p_{in} (GeV/c);n#sigma_{ele}^{TPC};n#sigma_{pion}^{TPC}",
                          GetVector(kP2D), GetVector(kSigmaEle), GetVector(kSigmaOther),
                          AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kTPCnSigmaPio);
  }
  
  // Eta, Phi, Pt
  histos->UserHistogram("Track","Eta","",200,-2,2,AliDielectronVarManager::kEta);
  histos->UserHistogram("Track","Phi","",120,0.,TMath::TwoPi(),AliDielectronVarManager::kPhi);
  histos->UserHistogram("Track","Eta_Phi","",GetVector(kEta2D),GetVector(kPhi2D), AliDielectronVarManager::kEta,AliDielectronVarManager::kPhi);
  histos->UserHistogram("Track","Eta_Pt","",GetVector(kP2D),GetVector(kEta2D), AliDielectronVarManager::kPt,AliDielectronVarManager::kEta);
  histos->UserHistogram("Track","Phi_Pt","",GetVector(kP2D),GetVector(kPhi2D), AliDielectronVarManager::kPt,AliDielectronVarManager::kPhi);
  
  // DCA
  histos->UserHistogram("Track","dXY","",200,-2.,2.,AliDielectronVarManager::kImpactParXY);
  histos->UserHistogram("Track","dZ" ,"",400,-4.,4.,AliDielectronVarManager::kImpactParZ);
  histos->UserHistogram("Track","dXY_dZ","",100,-1.,1.,150,-3.,3.,AliDielectronVarManager::kImpactParXY,AliDielectronVarManager::kImpactParZ);
  
  // Quality
  histos->UserHistogram("Track","TPCcrossedRowsOverFindable",";TPC crossed rows over findable clusters;#tracks",120,0.,1.2,AliDielectronVarManager::kNFclsTPCfCross);
  histos->UserHistogram("Track","TPCcrossedRows",";TPC crossed rows;#tracks",160,-0.5,159.5,AliDielectronVarManager::kNFclsTPCr);
  histos->UserHistogram("Track","TPCnCls",";TPC number clusters;#tracks",160,-0.5,159.5,AliDielectronVarManager::kNclsTPC);
  histos->UserHistogram("Track","ITSnCls",";ITS number clusters;#tracks",160,-0.5,159.5,AliDielectronVarManager::kNclsITS);
  histos->UserHistogram("Track","TPCchi2",";TPC chi2/Cl;#tracks",100,0.,10.,AliDielectronVarManager::kTPCchi2Cl);
  histos->UserHistogram("Track","ITSchi2",";ITS chi2/Cl;#tracks",100,0.,10.,AliDielectronVarManager::kITSchi2Cl);
  histos->UserHistogram("Track","NclsSFracTPC",";TPC fraction of shared clusters;#tracks",200,0,10.,AliDielectronVarManager::kNclsSFracTPC);
  histos->UserHistogram("Track","TPCclsDiff",";TPC cluster difference;#tracks",200,0,20.,AliDielectronVarManager::kTPCclsDiff);
  histos->UserHistogram("Track","TPCsignalN",";TPC number PID clusters;#tracks",160,-0.5,159.5,AliDielectronVarManager::kTPCsignalN);
  
  if (fIsQATask) {
    histos->UserHistogram("Track","TPCcrossedRows_TPCnCls",";TPC number clusters;TPC crossed rows",
                          160,-0.5,159.5, 160,-0.5,159.5, AliDielectronVarManager::kNclsTPC,AliDielectronVarManager::kNFclsTPCr);
    histos->UserHistogram("Track","TPCcrossedRows_Pt",";Pt [GeV];TPC crossed rows",
                          GetVector(kP2D), BinsToVector(160,-0.5,159.5), AliDielectronVarManager::kPt,AliDielectronVarManager::kNFclsTPCr);
    histos->UserHistogram("Track","TPCcrossedRowsOverFindable_Pt",";Pt [GeV];TPC crossed rows over findable",
                          GetVector(kP2D), BinsToVector(120,0.,1.2), AliDielectronVarManager::kPt,AliDielectronVarManager::kNFclsTPCfCross);
    histos->UserHistogram("Track","TPCcrossedRowsOverFindable_Eta",";Eta;TPC crossed rows over findable",
                          100,-1,1, 120,0.,1.2, AliDielectronVarManager::kEta,AliDielectronVarManager::kNFclsTPCfCross);
    histos->UserHistogram("Track","TPCcrossedRowsOverFindable_Phi",";Phi;TPC crossed rows over findable",
                          120,0.,TMath::TwoPi(), 120,0.,1.2, AliDielectronVarManager::kPhi,AliDielectronVarManager::kNFclsTPCfCross);
    //V0
    // these are pair variables, but only internal in the V0-finder...
    //histos->UserHistogram("Track","ArmenterosAlpha_ArmenterosPt",";kArmPt;kArmAlpha",
    //                      100,0.,0.5, 120,-0.6,0.6, AliDielectronVarManager::kArmPt,AliDielectronVarManager::kArmAlpha);
    //histos->UserHistogram("Track","TPCnSigmaEle_ArmenterosPt",";kArmPt;TPCnSigmaEle",
    //                      100,0.,0.5, 160,-12.,20., AliDielectronVarManager::kArmPt,AliDielectronVarManager::kTPCnSigmaEle);
  }
  
  
  //add histograms to Pair classes
  if (!fIsQATask) 
  {
    histos->UserHistogram("Pair","InvMass","",500,0.,5.,AliDielectronVarManager::kM);
    histos->UserHistogram("Pair","PairPt","",160,0.,8., AliDielectronVarManager::kPt);
    histos->UserHistogram("Pair","Rapidity","",200,-2.,2.,AliDielectronVarManager::kY);
    histos->UserHistogram("Pair","OpeningAngle","",315,0.,3.15,AliDielectronVarManager::kOpeningAngle);
    
    // 2D pair histograms
    histos->UserHistogram("Pair","InvMass_PairPt",";Inv. Mass [GeV];Pair Pt [GeV];#pairs",
                          GetVector(kMee), GetVector(kPtee),
                          AliDielectronVarManager::kM, AliDielectronVarManager::kPt);
    histos->UserHistogram("Pair","Eta_Phi_Pair",";Eta;Phi;#pairs",
                          100,-1.,1., 120,0.,TMath::TwoPi(),
                          AliDielectronVarManager::kEta, AliDielectronVarManager::kPhi);
    // 3D pair histograms
    histos->UserHistogram("Pair","InvMass_PairPt_PhivPair",";Inv. Mass [GeV];Pair Pt [GeV];PhiV",
                          GetVector(kMee), GetVector(kPtee), GetVector(kPhiV), 
                          AliDielectronVarManager::kM, AliDielectronVarManager::kPt, AliDielectronVarManager::kPhivPair);
    histos->UserHistogram("Pair","InvMass_PairPt_OpeningAngle",";Inv. Mass [GeV];Pair Pt [GeV];Opening Angle",
                          GetVector(kMee), GetVector(kPtee), GetVector(kOpAng), 
                          AliDielectronVarManager::kM, AliDielectronVarManager::kPt, AliDielectronVarManager::kOpeningAngle);
//    histos->UserHistogram("Pair","InvMass_PhivPair_OpeningAngle",";Inv. Mass [GeV];PhiV;Opening Angle",
//                          GetVector(kMee500), GetVector(kPhiV), GetVector(kOpAng2), 
//                          AliDielectronVarManager::kM, AliDielectronVarManager::kPhivPair, AliDielectronVarManager::kOpeningAngle);
//    histos->UserHistogram("Pair","InvMass_PairPt_Rapidity",";Inv. Mass [GeV];Pair Pt [GeV];Rapidity",
//                          GetVector(kMee), GetVector(kPtee), GetVector(kY3D), 
//                          AliDielectronVarManager::kM, AliDielectronVarManager::kPt, AliDielectronVarManager::kY);
    
    // more pair histograms
    if (die->DoEventProcess()) {
      //opening angle and PhiV
      histos->UserHistogram("Pair","InvMass_OpeningAngle",";Inv. Mass [GeV];Opening Angle;#pairs",
                            GetVector(kMee), GetVector(kOpAng), 
                            AliDielectronVarManager::kM, AliDielectronVarManager::kOpeningAngle);
      histos->UserHistogram("Pair","InvMass_OpeningAngle_fine",";Inv. Mass [GeV];Opening Angle;#pairs",
                            BinsToVector(40,0.,0.2), GetVector(kOpAng3), 
                            AliDielectronVarManager::kM, AliDielectronVarManager::kOpeningAngle);
      histos->UserHistogram("Pair","InvMass_PhivPair",";Inv. Mass [GeV];PhiV;#pairs",
                            GetVector(kMee), GetVector(kPhiV), 
                            AliDielectronVarManager::kM, AliDielectronVarManager::kPhivPair);
      histos->UserHistogram("Pair","PairPt_OpeningAngle",";Pair Pt [GeV];Opening Angle;#pairs",
                            GetVector(kPtee), GetVector(kOpAng), 
                            AliDielectronVarManager::kPt, AliDielectronVarManager::kOpeningAngle);
      histos->UserHistogram("Pair","PairPt_PhivPair",";Pair Pt [GeV];PhiV;#pairs",
                            GetVector(kPtee), GetVector(kPhiV), 
                            AliDielectronVarManager::kPt, AliDielectronVarManager::kPhivPair);
      histos->UserHistogram("Pair","OpeningAngle_PhivPair",";Opening Angle;PhiV;#pairs",
                            GetVector(kOpAng), GetVector(kPhiV), 
                            AliDielectronVarManager::kOpeningAngle, AliDielectronVarManager::kPhivPair);
      
      // pair DCA
      histos->UserHistogram("Pair","InvMass_PairDCAsigXY","",
                            GetVector(kMee), GetVector(kPairDCAsigXY),
                            AliDielectronVarManager::kM, AliDielectronVarManager::kPairDCAsigXY);
      histos->UserHistogram("Pair","InvMass_PairDCAabsXY","",
                            GetVector(kMee), GetVector(kPairDCAabsXY),
                            AliDielectronVarManager::kM, AliDielectronVarManager::kPairDCAabsXY);
      histos->UserHistogram("Pair","InvMass_PairLinDCAsigXY","",
                            GetVector(kMee), GetVector(kPairLinDCAsigXY),
                            AliDielectronVarManager::kM, AliDielectronVarManager::kPairLinDCAsigXY);
      histos->UserHistogram("Pair","InvMass_PairLinDCAabsXY","",
                            GetVector(kMee), GetVector(kPairLinDCAabsXY),
                            AliDielectronVarManager::kM, AliDielectronVarManager::kPairLinDCAabsXY);
//      // 3D
//      histos->UserHistogram("Pair","InvMass_PairDCAsigXY_PairLinDCAsigXY","",
//                            GetVector(kMee), GetVector(kPairDCAsigXY), GetVector(kPairLinDCAsigXY), 
//                            AliDielectronVarManager::kM, AliDielectronVarManager::kPairDCAsigXY, AliDielectronVarManager::kPairLinDCAsigXY);
    }
    // 3D
    histos->UserHistogram("Pair","InvMass_PairDCAsigXY_PairPt","",
                          GetVector(kMee), GetVector(kPairDCAsigXY), GetVector(kPtee), 
                          AliDielectronVarManager::kM, AliDielectronVarManager::kPairDCAsigXY, AliDielectronVarManager::kPt);
    histos->UserHistogram("Pair","InvMass_PairLinDCAsigXY_PairPt","",
                          GetVector(kMee), GetVector(kPairLinDCAsigXY), GetVector(kPtee), 
                          AliDielectronVarManager::kM, AliDielectronVarManager::kPairLinDCAsigXY, AliDielectronVarManager::kPt);
    
    // pairs vs centrality
    histos->UserHistogram("Pair","InvMass_Centrality",";Inv. Mass [GeV];Centrality;#pairs",
                          GetVector(kMee), BinsToVector(102,-1,101), 
                          AliDielectronVarManager::kM, AliDielectronVarManager::kCentrality);
    histos->UserHistogram("Pair","PairPt_Centrality",";Pair Pt [GeV];Centrality;#pairs",
                          GetVector(kPtee), BinsToVector(102,-1,101), 
                          AliDielectronVarManager::kPt, AliDielectronVarManager::kCentrality);
    
    // the histograms "Pre" may produce a memory leak:
    //    W-TROOT::Append: Replacing existing TH1: Pt (Potential memory leak).
    //    <message repeated 4 times>
    
    //    //add histograms to Track classes
    //    histos->UserHistogram("Pre","Pt",";Pt [GeV];#tracks",200,0,10.,AliDielectronVarManager::kPt);
    //    histos->UserHistogram("Pre","Px",";Px [GeV];#tracks",200,0,10.,AliDielectronVarManager::kPx);
    //    histos->UserHistogram("Pre","Py",";Py [GeV];#tracks",200,0,10.,AliDielectronVarManager::kPy);
    //    histos->UserHistogram("Pre","Pz",";Pz [GeV];#tracks",200,0,10.,AliDielectronVarManager::kPz);
    //    
    //    histos->UserHistogram("Pre","ITS_dEdx_P",";p (GeV/c);ITS signal (arb units)",
    //                          GetVector(kP2D), BinsToVector(700,0.,700.), AliDielectronVarManager::kP,AliDielectronVarManager::kITSsignal);
    //    histos->UserHistogram("Pre","TPC_dEdx_P",";p_{in} (GeV/c);TPC signal (arb units)",
    //                          GetVector(kP2D), BinsToVector(120,0.,120.), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCsignal);
    //    
    //    histos->UserHistogram("Pre","ITSnSigmaEle_P",";p (GeV/c);n#sigma_{ele}^{ITS}",
    //                          GetVector(kP2D), GetVector(kSigmaEle), AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaEle);
    //    histos->UserHistogram("Pre","ITSnSigmaPio_P",";p (GeV/c);n#sigma_{pion}^{ITS}",
    //                          GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaPio);
    //    histos->UserHistogram("Pre","ITSnSigmaKao_P",";p (GeV/c);n#sigma_{kaon}^{ITS}",
    //                          GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaKao);
    //    //histos->UserHistogram("Pre","ITSnSigmaPro_P",";p (GeV/c);n#sigma_{proton}^{ITS}",
    //    //                      GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaPro);
    //    
    //    histos->UserHistogram("Pre","TPCnSigmaEle_P",";p_{in} (GeV/c);n#sigma_{ele}^{TPC}",
    //                          GetVector(kP2D), GetVector(kSigmaEle), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle);
    //    histos->UserHistogram("Pre","TPCnSigmaPio_P",";p_{in} (GeV/c);n#sigma_{pion}^{TPC}",
    //                          GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaPio);
    //    histos->UserHistogram("Pre","TPCnSigmaKao_P",";p_{in} (GeV/c);n#sigma_{kaon}^{TPC}",
    //                          GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaKao);
    //    histos->UserHistogram("Pre","TPCnSigmaPro_P",";p_{in} (GeV/c);n#sigma_{proton}^{TPC}",
    //                          GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaPro);
    
  }// end: (!fIsQATask)
  
  die->SetHistogramManager(histos);
}


//_______________________________________________________________________________________________
TVectorD* LMEECutLib::GetVector(Int_t var) {
  switch (var) {
    case kPhiV:   return AliDielectronHelper::MakeLinBinning(100, 0., TMath::Pi());
    case kOpAng:  return AliDielectronHelper::MakeLinBinning(100, 0., TMath::Pi());
    case kOpAng2: return AliDielectronHelper::MakeLinBinning( 50, 0., TMath::Pi()/2.);
    case kOpAng3: return AliDielectronHelper::MakeLinBinning( 50, 0., 0.5);
    case kEta2D:  return AliDielectronHelper::MakeLinBinning(100,-1,1);
    case kEta3D:  return AliDielectronHelper::MakeLinBinning( 50,-1,1);
    case kPhi2D:  return AliDielectronHelper::MakeLinBinning(120, 0., TMath::TwoPi());
    case kY3D:    return AliDielectronHelper::MakeLinBinning( 20,-1.,1.);
      
    case kSigmaEle:
      if (fIsQATask) return AliDielectronHelper::MakeLinBinning(160,-12.,20.);
      else           return AliDielectronHelper::MakeLinBinning( 50, -5., 5.);
    case kSigmaOther:
      if (fIsQATask) return AliDielectronHelper::MakeLinBinning(100,-20.,20.);
      else           return AliDielectronHelper::MakeLinBinning( 50,-10.,10.);
    case kTPCdEdx:
      if (fIsQATask) return AliDielectronHelper::MakeLinBinning(120,  0.,120.);
      else           return AliDielectronHelper::MakeLinBinning( 50, 50.,100.);
      
    case kMee:    return AliDielectronHelper::MakeArbitraryBinning("0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 
                                                                   0.10, 0.14, 0.18, 0.22, 0.26, 0.30, 0.34, 0.38, 0.42, 0.46, 
                                                                   0.50, 0.54, 0.58, 0.62, 0.66, 0.70, 0.74, 0.78, 0.82, 0.86, 
                                                                   0.90, 0.94, 0.98, 1.02, 1.06, 
                                                                   1.10, 1.30, 1.50, 1.70, 1.90, 2.10, 2.30, 2.50, 2.70, 2.90, 
                                                                   3.10, 3.30, 3.50, 4.00, 4.50, 5.00 
                                                                   ");
    case kMee500: return AliDielectronHelper::MakeArbitraryBinning("0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 
                                                                   0.10, 0.14, 0.18, 0.22, 0.26, 0.30, 0.34, 0.38, 0.42, 0.46, 
                                                                   0.50 
                                                                   ");
    case kPtee:   return AliDielectronHelper::MakeArbitraryBinning("0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 
                                                                   0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 
                                                                   1.00, 1.20, 1.40, 1.60, 1.80, 2.00, 2.20, 2.40, 2.60, 2.80, 
                                                                   3.00, 3.20, 3.40, 3.60, 3.80, 4.00, 4.20, 4.40, 4.60, 4.80, 
                                                                   5.00, 6.00, 7.00, 8.00 
                                                                   ");
    case kP2D:    return AliDielectronHelper::MakeArbitraryBinning("0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 
                                                                   0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 
                                                                   1.00, 1.05, 1.10, 1.15, 1.20, 1.25, 1.30, 1.35, 1.40, 1.45, 
                                                                   1.50, 1.55, 1.60, 1.65, 1.70, 1.75, 1.80, 1.85, 1.90, 1.95, 
                                                                   2.00, 2.20, 2.40, 2.60, 2.80, 3.00, 3.20, 3.40, 3.60, 3.80, 
                                                                   4.00, 4.50, 5.00, 5.50, 6.00, 6.50, 7.00, 7.50, 8.00 
                                                                   ");
      //2.00, 2.05, 2.10, 2.15, 2.20, 2.25, 2.30, 2.35, 2.40, 2.45, 
      //2.50, 2.55, 2.60, 2.65, 2.70, 2.75, 2.80, 2.85, 2.90, 2.95, 
      //3.00, 3.05, 3.10, 3.15, 3.20, 3.25, 3.30, 3.35, 3.40, 3.45, 
      
    case kPairDCAsigXY:
    case kPairLinDCAsigXY: return AliDielectronHelper::MakeLinBinning(50, 0., 10.); // in sigma
    case kPairDCAabsXY:
    case kPairLinDCAabsXY: return AliDielectronHelper::MakeLogBinning(50, 0.0001, 1.); // in cm
      
    case kRuns:   return AliDielectronHelper::MakeArbitraryBinning("167900, 167987, 167988, 168310, 168311, 168322, 168325, 168341, 168342, 168361, 168362, 168458, 168460, 168464, 168467, 168511, 168512, 168514, 168777, 168826, 168988, 168992, 169035, 169040, 169044, 169045, 169091, 169094, 169099, 169138, 169144, 169145, 169148, 169156, 169160, 169167, 169238, 169411, 169415, 169417, 169418, 169419, 169420, 169475, 169498, 169504, 169506, 169512, 169515, 169550, 169553, 169554, 169555, 169557, 169586, 169587, 169588, 169590, 169591, 169835, 169837, 169838, 169846, 169855, 169858, 169859, 170027, 170040, 170081, 170083, 170084, 170085, 170088, 170089, 170091, 170155, 170159, 170163, 170193, 170203, 170204, 170207, 170228, 170230, 170268, 170269, 170270, 170306, 170308, 170309, 170311, 170312, 170313, 170315, 170387, 170388, 170572, 170593, 170600");
      
    default: cout << "ERROR: in 'GetVector(...var)' variable for axis range not defined!" << endl;
      break;
  } 
}


//_______________________________________________________________________________________________
TVectorD* LMEECutLib::BinsToVector(Int_t nbins, Double_t min, Double_t max) {
  return AliDielectronHelper::MakeLinBinning(nbins,min,max);
  //  TVectorD *vec = new TVectorD(nbins+1);
  //  Double_t binwdth = (max-min)/nbins;
  //  for (int i = 0; i < nbins+1; i++) (*vec)[i] = min + i*binwdth;
  //  return vec;
}


//_______________________________________________________________________________________________
void LMEECutLib::InitCF(AliDielectron* die, Int_t cutDefinition) {
  //
  // Setupd the CF Manager if needed
  //
  AliDielectronCF *cf=new AliDielectronCF(die->GetName(),die->GetTitle());
  
  //pair variables
  cf->AddVariable(AliDielectronVarManager::kP,100,0.,5.);
  cf->AddVariable(AliDielectronVarManager::kM,200,-0.01,3.99); //20Mev Steps
  cf->AddVariable(AliDielectronVarManager::kPairType,10,0,10);
  
  cf->AddVariable(AliDielectronVarManager::kCentrality,"0.,5.,10.,20.,30.,50.,80.,100.");
  
  //leg variables
  cf->AddVariable(AliDielectronVarManager::kP,160,0.,8.,kTRUE);
  cf->AddVariable(AliDielectronVarManager::kITSsignal,350,0.,700.,kTRUE);
  cf->AddVariable(AliDielectronVarManager::kTPCsignal,60,0.,120.,kTRUE);
  cf->AddVariable(AliDielectronVarManager::kHaveSameMother,21,-10,10,kTRUE);
  
  //only in this case write MC truth info
  if (die->GetHasMC()) {
    cf->SetStepForMCtruth();
    cf->SetStepsForMCtruthOnly();
    cf->AddVariable(AliDielectronVarManager::kPdgCode,10000,-5000.5,4999.5,kTRUE);
    cf->AddVariable(AliDielectronVarManager::kPdgCodeMother,10000,-5000.5,4999.5,kTRUE);
  }
  
  cf->SetStepsForSignal();
  die->SetCFManagerPair(cf);
}
