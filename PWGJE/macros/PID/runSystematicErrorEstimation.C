enum modes { kPt = 0, kZ = 1, kXi = 2};

//_____________________________________________________________________________________________________________________________
Bool_t CentralityHasDecimalsPlaces(Double_t cent)
{
  // Check, whether the centrality has decimals places or is an integer
  const Double_t temp1 = ((Int_t)cent)*1e6;
  const Double_t temp2 = cent*1e6;
  
  return TMath::Abs(temp1 - temp2) > 1;
}


//_____________________________________________________________________________________________________________________________
void setMeanThreshold(Double_t& setMeanLowerThreshold, Double_t& setMeanUpperThreshold, Int_t mode, Double_t lowerJetPt, Double_t upperJetPt)
{
  const Double_t setMeanLowerThresholdPt = 4;
  const Double_t setMeanUpperThresholdPt = 999.;
  
  // Take effective jetPt for z and xi
  Double_t effectiveJetPt = 0.5 * (lowerJetPt + upperJetPt);
  
  if (mode == kPt) {
    setMeanLowerThreshold = setMeanLowerThresholdPt;
    setMeanUpperThreshold = setMeanUpperThresholdPt;
  }
  else if (mode == kZ) {
    setMeanLowerThreshold = setMeanLowerThresholdPt / effectiveJetPt;
    setMeanUpperThreshold = setMeanUpperThresholdPt / effectiveJetPt;
  }
  else if (mode == kXi) {
    // Thresholds are swapped!
    setMeanLowerThreshold = TMath::Log(effectiveJetPt / setMeanUpperThresholdPt);
    setMeanUpperThreshold = TMath::Log(effectiveJetPt / setMeanLowerThresholdPt);
    if (setMeanUpperThresholdPt > 900)
      setMeanLowerThreshold = 0.;
  }
  else {
    printf("ERROR: Unknown mode!\n");
    exit(-1);
  }
  
}

void runSystematicErrorEstimation(TString chargeString /* "", "_negCharge", "_posCharge"*/)
{
  // NOTE:
  // If mean instead of default is to be used, just set setMean to kTRUE in the following.
  // ALSO: In addUpSystematicErrors the reference must be adapted to the corresponding mean (not the default).
  // Currently, ONLY the graph and the plot with the single outputs for the sys errors are changed.
  // TODO for average instead of default for spline systetematics at high pT > 4 GeV/c:
  // - Seems to make no difference more or less (at least for the old "symmetric" systetematics)
  // - Questionable to set the threshold as function of momentum! Could mean that K and p are move away from the best estimated
  //   because only the different models only affect pions in the rel. rise but for quite some pT not the other species
  // - Also a problem:
  // ATTENTION: Normally, weighted average should be done either for the yields or should be WEIGHTED with the total yield
  // (must be done anyway for the to-pion-ratios!). BUT: The whole approach relies on the assumption that the overall statistics
  // is the same for all compared data points. This must be checked anyway and is true within percent typically!
  //
  // => Looks like much to do and many problems with most likely little gain. Better don't use the average instead of the default!
  gROOT->LoadMacro("SystematicErrorEstimation.C+");
  
  const Double_t nSigma = 0;
  
  const Bool_t ignoreSigmaErrors = kTRUE;
  
  const Int_t numCentralities = 9;
  const Int_t centralities[numCentralities + 1] = { 0, 5, 10, 20, 30, 40, 50, 60, 80, 100 };
  
  const Int_t numCentralitiesRefMultOld = 9;
  const Int_t centralitiesRefMultOldLowEdge[numCentralitiesRefMultOld] = { 0,  7, 13, 20, 29, 40, 50, 60,   0 };
  const Int_t centralitiesRefMultOldUpEdge[numCentralitiesRefMultOld]  = { 7, 13, 20, 29, 40, 50, 60, 95, 125 };
  
  const Int_t numCentralitiesRefMultNew = 18;
  const Int_t centralitiesRefMultNewLowEdge[numCentralitiesRefMultNew] =
    { 1, 4,  7, 10, 15, 20, 25, 30, 40, 50, 60,  70,  100, 15,   25,   40,   60,    0 };
  const Int_t centralitiesRefMultNewUpEdge[numCentralitiesRefMultNew]  =
    { 4, 7, 10, 15, 20, 25, 30, 40, 50, 60, 70, 100, 9999, 25, 9999, 9999, 9999, 9999 };
  
  
  const Int_t numCentralitiesV0Mpp = 15;
  const Double_t centralitiesV0MppLowEdge[numCentralitiesV0Mpp] = {    0, 0.01, 0.1, 1,  5, 10, 15, 20, 30, 40, 50,  70,   0, 0,   0 };
  const Double_t centralitiesV0MppUpEdge[numCentralitiesV0Mpp]  = { 0.01,  0.1,   1, 5, 10, 15, 20, 30, 40, 50, 70, 100, 0.1, 1, 100 };
  
  TString centralityString = "";
  
  const Int_t numJetPtBins = 4;
  const Int_t jetPt[numJetPtBins + 3] = { 5, 10, 15, 20, 30, 10, 40 };
  
  TString jetPtString = "";
  
  Double_t setMeanLowerThreshold = 999.;
  Double_t setMeanUpperThreshold = 999.;
  
  
  /* OLD: All in one go -> Deprecated
  const TString path = "/hera/alice/bhess/analysis/13b.pass3";
  Bool_t useCentralities = kTRUE;
  Bool_t useJetPt = kFALSE;
  
  for (Int_t iCent = 0; iCent < numCentralities; iCent++) {
    centralityString = useCentralities ? Form("_centrality // 10-80 bin
      iJetPt++;
    }%d_%d", centralities[iCent], centralities[iCent + 1]) : "";
    
    const Int_t numFiles = 13;
      
    TString outFileTitle = Form("Combined%s", centralityString.Data());
      
    TString fileNames[numFiles];
    TString histTitles[numFiles];
    
    // Default
    fileNames[0] = 
      Form("%s/bhess_PID_TPCdefaultPriors_ITS_PureGauss_results_LLFit__Pt_2_reg1_regFac1.00_muonToElTunedOnMCStandardTrackCutsPPB_idSpectra%s.root"
            , path.Data(), centralityString.Data());
    
    histTitles[0] = "Default settings";
    
    // Centrality estimator
    fileNames[1] = 
      Form("%s/bhess_PID_TPCdefaultPriors_ITS_PureGauss_SystematicsV0AcentEstimator_results_LLFit__Pt_2_reg1_regFac1.00_muonToElTunedOnMCStandardTrackCutsPPB_idSpectra%s.root"
            , path.Data(), centralityString.Data());
    
    histTitles[1] = "V0A";
    
    // Multiplicity correction
    fileNames[2] = 
      Form("%s/bhess_PID_TPCdefaultPriors_ITS_PureGauss_SystematicsMultDown_results_LLFit__Pt_2_reg1_regFac1.00_muonToElTunedOnMCStandardTrackCutsPPB_idSpectra%s.root"
            , path.Data(), centralityString.Data());
    fileNames[3] = 
      Form("%s/bhess_PID_TPCdefaultPriors_ITS_PureGauss_SystematicsMultUp_results_LLFit__Pt_2_reg1_regFac1.00_muonToElTunedOnMCStandardTrackCutsPPB_idSpectra%s.root"
            , path.Data(), centralityString.Data());
    
    histTitles[2] = "Mult. correction - 0.2%";
    histTitles[3] = "Mult. correction + 0.2%";
    
    // Pre-PID
    fileNames[4] = 
      Form("%s/bhess_PID_TPCdefaultPriors_ITS_PureGauss_results_LLFit__Pt_2_reg1_regFac1.00_muonToElTunedOnMCStandardTrackCutsPPB%s.root"
            , path.Data(), centralityString.Data());
    
    histTitles[4] = "No Pre-PID";
    
    // Shape of det. resp.
    fileNames[5] = 
      Form("%s/bhess_PID_TPCdefaultPriors_ITS_results_LLFit__Pt_2_reg1_regFac1.00_muonToElTunedOnMCStandardTrackCutsPPB_idSpectra%s.root"
            , path.Data(), centralityString.Data());
      
    histTitles[5] = "Asymmetric shape ";
    
    // Sigma parametrisation
    fileNames[6] = 
      Form("%s/bhess_PID_TPCdefaultPriors_ITS_PureGauss_SystematicsSigmaDown_results_LLFit__Pt_2_reg1_regFac1.00_muonToElTunedOnMCStandardTrackCutsPPB_idSpectra%s.root"
            , path.Data(), centralityString.Data());
    fileNames[7] = 
      Form("%s/bhess_PID_TPCdefaultPriors_ITS_PureGauss_SystematicsSigmaUp_results_LLFit__Pt_2_reg1_regFac1.00_muonToElTunedOnMCStandardTrackCutsPPB_idSpectra%s.root"
            , path.Data(), centralityString.Data());
    
    histTitles[6] = "#sigma correction - 5%";
    histTitles[7] = "#sigma correction + 5%";

    // Eta correction
    fileNames[8] = 
      Form("%s/bhess_PID_TPCdefaultPriors_ITS_PureGauss_SystematicsEtaDown_results_LLFit__Pt_2_reg1_regFac1.00_muonToElTunedOnMCStandardTrackCutsPPB_idSpectra%s.root"
            , path.Data(), centralityString.Data());
    fileNames[9] = 
      Form("%s/bhess_PID_TPCdefaultPriors_ITS_PureGauss_SystematicsEtaUp_results_LLFit__Pt_2_reg1_regFac1.00_muonToElTunedOnMCStandardTrackCutsPPB_idSpectra%s.root"
            , path.Data(), centralityString.Data());
    
    histTitles[8] = "#eta correction - 1%, p > 0.45 GeV/c; - %4, p < 0.45 GeV/c";
    histTitles[9] = "#eta correction + 1%, p > 0.45 GeV/c; + %4, p < 0.45 GeV/c";

    // Splines
    fileNames[10] =
      Form("%s/bhess_PID_TPCdefaultPriors_ITS_PureGauss_SystematicsSplinesDown_results_LLFit__Pt_2_reg1_regFac1.00_muonToElTunedOnMCStandardTrackCutsPPB_idSpectra%s.root"
            , path.Data(), centralityString.Data());
    fileNames[11] =
      Form("%s/bhess_PID_TPCdefaultPriors_ITS_PureGauss_SystematicsSplinesUp_results_LLFit__Pt_2_reg1_regFac1.00_muonToElTunedOnMCStandardTrackCutsPPB_idSpectra%s.root"
            , path.Data(), centralityString.Data());
    
    histTitles[10] = "Splines - 0.2%";
    histTitles[11] = "Splines + 0.2%";

    // Muon treatment
    fileNames[12] =
      Form("%s/bhess_PID_TPCdefaultPriors_ITS_PureGauss_results_LLFit__Pt_2_reg1_regFac1.00_muonsEqualElectrons_idSpectra%s.root"
            , path.Data(), centralityString.Data());
    //fileNames[13] =
      //Form("%s/bhess_PID_TPCdefaultPriors_ITS_PureGauss_results_LLFit__Pt_2_reg1_regFac1.00_noMuons_idSpectra%s.root"
      //     , path.Data(), centralityString.Data());
    
    histTitles[12] = "Muon_Electron_Ratio_Is_Unity";
    //histTitles[13] = "NoMuons";

    SystematicErrorEstimation(path, outFileTitle, fileNames, histTitles, numFiles, nSigma, ignoreSigmaErrors, setMean, setMeanLowerThreshold,
                              setMeanUpperThreshold); 
    
    if (!useCentralities)
      break;
  }*/
  
  
  /* pPb
  const TString path = "/hera/alice/bhess/analysis/13b.pass3";
  Bool_t useCentralities = kTRUE;
  Bool_t useJetPt = kFALSE;
  
  // In case of taking the average of the systematics: Only do this within certain range, HERE: pT only
  setMeanThreshold(setMeanLowerThreshold, setMeanUpperThreshold, kPt, -1, -1);
  
  for (Int_t iCent = 0; iCent < numCentralities; iCent++) {
    centralityString = useCentralities ? Form("_centrality%d_%d", centralities[iCent], centralities[iCent + 1]) : "";
    
    {
      Bool_t setMean = kFALSE;
      const Int_t numFiles = 2;
      
      TString outFileTitle = Form("Centrality%s", centralityString.Data());
      
      TString fileNames[numFiles];
      fileNames[0] = 
        Form("%s/bhess_PID_TPCdefaultPriors_ITS_PureGauss_results_LLFit__Pt_2_reg1_regFac1.00_muonToElTunedOnMCStandardTrackCutsPPB_idSpectra%s.root"
             , path.Data(), centralityString.Data());
      fileNames[1] = 
        Form("%s/bhess_PID_TPCdefaultPriors_ITS_PureGauss_SystematicsV0AcentEstimator_results_LLFit__Pt_2_reg1_regFac1.00_muonToElTunedOnMCStandardTrackCutsPPB_idSpectra%s.root"
             , path.Data(), centralityString.Data());
      
      TString histTitles[numFiles];
      histTitles[0] = "V0M (default)";
      histTitles[1] = "V0A";
      
      SystematicErrorEstimation(path, outFileTitle, fileNames, histTitles, numFiles, nSigma, ignoreSigmaErrors, setMean, setMeanLowerThreshold,
                                setMeanUpperThreshold); 
    }
    {
      Bool_t setMean = kFALSE;
      const Int_t numFiles = 3;
      
      TString outFileTitle = Form("Multiplicity%s", centralityString.Data());
      
      TString fileNames[numFiles];
      fileNames[0] = 
        Form("%s/bhess_PID_TPCdefaultPriors_ITS_PureGauss_results_LLFit__Pt_2_reg1_regFac1.00_muonToElTunedOnMCStandardTrackCutsPPB_idSpectra%s.root"
             , path.Data(), centralityString.Data());
      fileNames[1] = 
        Form("%s/bhess_PID_TPCdefaultPriors_ITS_PureGauss_SystematicsMultDown_results_LLFit__Pt_2_reg1_regFac1.00_muonToElTunedOnMCStandardTrackCutsPPB_idSpectra%s.root"
             , path.Data(), centralityString.Data());
      fileNames[2] = 
        Form("%s/bhess_PID_TPCdefaultPriors_ITS_PureGauss_SystematicsMultUp_results_LLFit__Pt_2_reg1_regFac1.00_muonToElTunedOnMCStandardTrackCutsPPB_idSpectra%s.root"
             , path.Data(), centralityString.Data());
      
      TString histTitles[numFiles];
      histTitles[0] = "Standard mult. correction";
      histTitles[1] = "Mult. correction - 0.2%";
      histTitles[2] = "Mult. correction + 0.2%";
      
      SystematicErrorEstimation(path, outFileTitle, fileNames, histTitles, numFiles, nSigma, ignoreSigmaErrors, setMean, setMeanLowerThreshold,
                                setMeanUpperThreshold);  
    }
    {
      Bool_t setMean = kFALSE;
      const Int_t numFiles = 2;
      
      TString outFileTitle = Form("PrePID%s", centralityString.Data());
      
      TString fileNames[numFiles];
      fileNames[0] = 
        Form("%s/bhess_PID_TPCdefaultPriors_ITS_PureGauss_results_LLFit__Pt_2_reg1_regFac1.00_muonToElTunedOnMCStandardTrackCutsPPB_idSpectra%s.root"
             , path.Data(), centralityString.Data());
      fileNames[1] = 
        Form("%s/bhess_PID_TPCdefaultPriors_ITS_PureGauss_results_LLFit__Pt_2_reg1_regFac1.00_muonToElTunedOnMCStandardTrackCutsPPB%s.root"
             , path.Data(), centralityString.Data());
      
      TString histTitles[numFiles];
      histTitles[0] = "PID combined (default)";
      histTitles[1] = "No PID";
      
      SystematicErrorEstimation(path, outFileTitle, fileNames, histTitles, numFiles, nSigma, ignoreSigmaErrors, setMean, setMeanLowerThreshold,
                                setMeanUpperThreshold);  
    }
    {
      Bool_t setMean = kFALSE;
      const Int_t numFiles = 2;
      
      TString outFileTitle = Form("Shape%s", centralityString.Data());
      
      TString fileNames[numFiles];
      fileNames[0] = 
        Form("%s/bhess_PID_TPCdefaultPriors_ITS_PureGauss_results_LLFit__Pt_2_reg1_regFac1.00_muonToElTunedOnMCStandardTrackCutsPPB_idSpectra%s.root"
             , path.Data(), centralityString.Data());
      fileNames[1] = 
        Form("%s/bhess_PID_TPCdefaultPriors_ITS_results_LLFit__Pt_2_reg1_regFac1.00_muonToElTunedOnMCStandardTrackCutsPPB_idSpectra%s.root"
             , path.Data(), centralityString.Data());
        
      TString histTitles[numFiles];
      histTitles[0] = "Pure gauss shape (default)";
      histTitles[1] = "Asymmetric shape ";
      
      SystematicErrorEstimation(path, outFileTitle, fileNames, histTitles, numFiles, nSigma, ignoreSigmaErrors, setMean, setMeanLowerThreshold,
                                setMeanUpperThreshold);  
    }
    {
      Bool_t setMean = kFALSE;
      const Int_t numFiles = 3;
      
      TString outFileTitle = Form("Sigma%s", centralityString.Data());
      
      TString fileNames[numFiles];
      fileNames[0] = 
        Form("%s/bhess_PID_TPCdefaultPriors_ITS_PureGauss_results_LLFit__Pt_2_reg1_regFac1.00_muonToElTunedOnMCStandardTrackCutsPPB_idSpectra%s.root"
             , path.Data(), centralityString.Data());
      fileNames[1] = 
        Form("%s/bhess_PID_TPCdefaultPriors_ITS_PureGauss_SystematicsSigmaDown_results_LLFit__Pt_2_reg1_regFac1.00_muonToElTunedOnMCStandardTrackCutsPPB_idSpectra%s.root"
             , path.Data(), centralityString.Data());
      fileNames[2] = 
        Form("%s/bhess_PID_TPCdefaultPriors_ITS_PureGauss_SystematicsSigmaUp_results_LLFit__Pt_2_reg1_regFac1.00_muonToElTunedOnMCStandardTrackCutsPPB_idSpectra%s.root"
             , path.Data(), centralityString.Data());
      
      TString histTitles[numFiles];
      histTitles[0] = "Standard #sigma map";
      histTitles[1] = "#sigma correction - 5%";
      histTitles[2] = "#sigma correction + 5%";
      
      SystematicErrorEstimation(path, outFileTitle, fileNames, histTitles, numFiles, nSigma, ignoreSigmaErrors, setMean, setMeanLowerThreshold,
                                setMeanUpperThreshold);  
    }
    {
      Bool_t setMean = kFALSE;
      const Int_t numFiles = 3;
      
      TString outFileTitle = Form("Eta%s", centralityString.Data());
      
      TString fileNames[numFiles];
      fileNames[0] = 
        Form("%s/bhess_PID_TPCdefaultPriors_ITS_PureGauss_results_LLFit__Pt_2_reg1_regFac1.00_muonToElTunedOnMCStandardTrackCutsPPB_idSpectra%s.root"
             , path.Data(), centralityString.Data());
      fileNames[1] = 
        Form("%s/bhess_PID_TPCdefaultPriors_ITS_PureGauss_SystematicsEtaDown_results_LLFit__Pt_2_reg1_regFac1.00_muonToElTunedOnMCStandardTrackCutsPPB_idSpectra%s.root"
             , path.Data(), centralityString.Data());
      fileNames[2] = 
        Form("%s/bhess_PID_TPCdefaultPriors_ITS_PureGauss_SystematicsEtaUp_results_LLFit__Pt_2_reg1_regFac1.00_muonToElTunedOnMCStandardTrackCutsPPB_idSpectra%s.root"
             , path.Data(), centralityString.Data());
      
      TString histTitles[numFiles];
      histTitles[0] = "Standard #eta map";
      histTitles[1] = "#eta correction - 1%, p > 0.45 GeV/c; - %4, p < 0.45 GeV/c";
      histTitles[2] = "#eta correction + 1%, p > 0.45 GeV/c; + %4, p < 0.45 GeV/c";

      SystematicErrorEstimation(path, outFileTitle, fileNames, histTitles, numFiles, nSigma, ignoreSigmaErrors, setMean, setMeanLowerThreshold,
                                setMeanUpperThreshold);  
    }
    {
      Bool_t setMean = kFALSE;
      const Int_t numFiles = 3;
        
      TString outFileTitle = Form("Splines%s", centralityString.Data());
      TString fileNames[numFiles];
      fileNames[0] = 
        Form("%s/bhess_PID_TPCdefaultPriors_ITS_PureGauss_results_LLFit__Pt_2_reg1_regFac1.00_muonToElTunedOnMCStandardTrackCutsPPB_idSpectra%s.root"
             , path.Data(), centralityString.Data());
      fileNames[1] =
        Form("%s/bhess_PID_TPCdefaultPriors_ITS_PureGauss_SystematicsSplinesDown_results_LLFit__Pt_2_reg1_regFac1.00_muonToElTunedOnMCStandardTrackCutsPPB_idSpectra%s.root"
             , path.Data(), centralityString.Data());
      fileNames[2] =
        Form("%s/bhess_PID_TPCdefaultPriors_ITS_PureGauss_SystematicsSplinesUp_results_LLFit__Pt_2_reg1_regFac1.00_muonToElTunedOnMCStandardTrackCutsPPB_idSpectra%s.root"
             , path.Data(), centralityString.Data());
      
      TString histTitles[numFiles];
      histTitles[0] = "Standard splines";
      histTitles[1] = "Splines - 0.2%";
      histTitles[2] = "Splines + 0.2%";
      
      SystematicErrorEstimation(path, outFileTitle, fileNames, histTitles, numFiles, nSigma, ignoreSigmaErrors, setMean, setMeanLowerThreshold,
                                setMeanUpperThreshold);
    }
    {
      Bool_t setMean = kFALSE;
      const Int_t numFiles = 2;//3;
      
      TString outFileTitle = Form("Muons%s", centralityString.Data());
      
      TString fileNames[numFiles];
      fileNames[0] =
        Form("%s/bhess_PID_TPCdefaultPriors_ITS_PureGauss_results_LLFit__Pt_2_reg1_regFac1.00_muonToElTunedOnMCStandardTrackCutsPPB_idSpectra%s.root"
             , path.Data(), centralityString.Data());
      fileNames[1] =
        Form("%s/bhess_PID_TPCdefaultPriors_ITS_PureGauss_results_LLFit__Pt_2_reg1_regFac1.00_muonsEqualElectrons_idSpectra%s.root"
             , path.Data(), centralityString.Data());
      //fileNames[2] =
        //Form("%s/bhess_PID_TPCdefaultPriors_ITS_PureGauss_results_LLFit__Pt_2_reg1_regFac1.00_noMuons_idSpectra%s.root"
        //     , path.Data(), centralityString.Data());
      
      TString histTitles[numFiles];
      histTitles[0] = "Muon_Electron_Ratio_Tuned_To_MC";
      histTitles[1] = "Muon_Electron_Ratio_Is_Unity";
      //histTitles[2] = "NoMuons";

      SystematicErrorEstimation(path, outFileTitle, fileNames, histTitles, numFiles, nSigma, ignoreSigmaErrors, setMean, setMeanLowerThreshold,
                                setMeanUpperThreshold);  
    }
    
    if (!useCentralities)
      break;
  }*/
  
  
  // pp Inclusive, no muons
  /*
  //const TString path = "/hera/alice/bhess/analysis/10d_e.pass2_merged/nclCut/differentTrackFilterBits/528";
  //const TString path = "/hera/alice/bhess/analysis/10d_e.pass2_merged/nclCut";
  //const TString path = "/hera/alice/bhess/analysis/10d.pass2/nclCut/";
  //const TString path = "/hera/alice/bhess/analysis/10e.pass2/nclCut/";
  //const TString path = "/hera/alice/bhess/analysis/10f6a/nclCut/noMCidForGen";
  //const TString path = "/hera/alice/bhess/analysis/10f6a/nclCut/noMCidForGen/filterBit528";
  const TString path = "/hera/alice/bhess/analysis/10f6a/nclCut/noMCidForGen/filterBit528/closureCheck";
  Bool_t useCentralities = kFALSE;
  Bool_t useJetPt = kFALSE;
  
  // Inclusive means: Only pT
  setMeanThreshold(setMeanLowerThreshold, setMeanUpperThreshold, kPt, -1, -1);
  

  for (Int_t iJetPt = 0; iJetPt <= numJetPtBins; iJetPt++) {
    if (iJetPt == numJetPtBins) {
      // 10-40 bin
      iJetPt++;
    }
    jetPtString = useJetPt ? Form("_jetPt%d.0_%d.0", jetPt[iJetPt], jetPt[iJetPt + 1]) : "";
    jetPtString.Append(chargeString);
    
    {
      Bool_t setMean = kFALSE;
      const Int_t numFiles = 2;
      
      TString outFileTitle = Form("PrePID%s", jetPtString.Data());
      
      TString fileNames[numFiles];
      fileNames[0] = 
        Form("%s/bhess_PID_Jets_Inclusive_results_LLFit__Pt_2_reg1_regFac1.00_noMuons_idSpectra%s.root"
             , path.Data(), jetPtString.Data());
      fileNames[1] = 
        Form("%s/bhess_PID_Jets_Inclusive_results_LLFit__Pt_2_reg1_regFac1.00_noMuons%s.root"
             , path.Data(), jetPtString.Data());
      
      TString histTitles[numFiles];
      histTitles[0] = "PID combined (default)";
      histTitles[1] = "No PID";
      
      SystematicErrorEstimation(path, outFileTitle, fileNames, histTitles, numFiles, nSigma, ignoreSigmaErrors, setMean, setMeanLowerThreshold,
                                setMeanUpperThreshold);  
    }
    {
      Bool_t setMean = kFALSE;
      const Int_t numFiles = 2;
      
      TString outFileTitle = Form("Shape%s", jetPtString.Data());
      
      TString fileNames[numFiles];
      fileNames[0] = 
        Form("%s/bhess_PID_Jets_Inclusive_results_LLFit__Pt_2_reg1_regFac1.00_noMuons_idSpectra%s.root"
             , path.Data(), jetPtString.Data());
      fileNames[1] = 
        Form("%s/bhess_PID_Jets_Inclusive_PureGauss_results_LLFit__Pt_2_reg1_regFac1.00_noMuons_idSpectra%s.root"
             , path.Data(), jetPtString.Data());
      
        
      TString histTitles[numFiles];
      histTitles[0] = "Asymmetric shape (default)";
      histTitles[1] = "Pure gauss shape ";
      
      SystematicErrorEstimation(path, outFileTitle, fileNames, histTitles, numFiles, nSigma, ignoreSigmaErrors, setMean, setMeanLowerThreshold,
                                setMeanUpperThreshold);  
    }
    {
      Bool_t setMean = kFALSE;
      const Int_t numFiles = 3;
      
      TString outFileTitle = Form("Sigma%s", jetPtString.Data());
      
      TString fileNames[numFiles];
      fileNames[0] = 
        Form("%s/bhess_PID_Jets_Inclusive_results_LLFit__Pt_2_reg1_regFac1.00_noMuons_idSpectra%s.root"
             , path.Data(), jetPtString.Data());
      fileNames[1] = 
        Form("%s/bhess_PID_Jets_Inclusive_SystematicsSigmaDown_results_LLFit__Pt_2_reg1_regFac1.00_noMuons_idSpectra%s.root"
             , path.Data(), jetPtString.Data());
      fileNames[2] = 
        Form("%s/bhess_PID_Jets_Inclusive_SystematicsSigmaUp_results_LLFit__Pt_2_reg1_regFac1.00_noMuons_idSpectra%s.root"
             , path.Data(), jetPtString.Data());
      
      TString histTitles[numFiles];
      histTitles[0] = "Standard #sigma map";
      histTitles[1] = "#sigma correction - 5%";
      histTitles[2] = "#sigma correction + 5%";
      
      SystematicErrorEstimation(path, outFileTitle, fileNames, histTitles, numFiles, nSigma, ignoreSigmaErrors, setMean, setMeanLowerThreshold,
                                setMeanUpperThreshold);  
    }
    {
      Bool_t setMean = kFALSE;
      const Int_t numFiles = 3;
      
      TString outFileTitle = Form("Eta%s", jetPtString.Data());
      
      TString fileNames[numFiles];
      fileNames[0] = 
        Form("%s/bhess_PID_Jets_Inclusive_results_LLFit__Pt_2_reg1_regFac1.00_noMuons_idSpectra%s.root"
             , path.Data(), jetPtString.Data());
      fileNames[1] = 
        Form("%s/bhess_PID_Jets_Inclusive_SystematicsEtaDown_results_LLFit__Pt_2_reg1_regFac1.00_noMuons_idSpectra%s.root"
             , path.Data(), jetPtString.Data());
      fileNames[2] = 
        Form("%s/bhess_PID_Jets_Inclusive_SystematicsEtaUp_results_LLFit__Pt_2_reg1_regFac1.00_noMuons_idSpectra%s.root"
             , path.Data(), jetPtString.Data());
      
      TString histTitles[numFiles];
      histTitles[0] = "Standard #eta map";
      histTitles[1] = "#eta correction - 1%, p > 0.45 GeV/c; - %4, p < 0.45 GeV/c";
      histTitles[2] = "#eta correction + 1%, p > 0.45 GeV/c; + %4, p < 0.45 GeV/c";

      SystematicErrorEstimation(path, outFileTitle, fileNames, histTitles, numFiles, nSigma, ignoreSigmaErrors, setMean, setMeanLowerThreshold,
                                setMeanUpperThreshold);  
    }
    {
      Bool_t setMean = kFALSE;
      const Int_t numFiles = 3;
        
      TString outFileTitle = Form("Splines%s", jetPtString.Data());
      TString fileNames[numFiles];
      fileNames[0] = 
        Form("%s/bhess_PID_Jets_Inclusive_results_LLFit__Pt_2_reg1_regFac1.00_noMuons_idSpectra%s.root"
             , path.Data(), jetPtString.Data());
      fileNames[1] =
        Form("%s/bhess_PID_Jets_Inclusive_SystematicsSplinesDown_results_LLFit__Pt_2_reg1_regFac1.00_noMuons_idSpectra%s.root"
             , path.Data(), jetPtString.Data());
      fileNames[2] =
        Form("%s/bhess_PID_Jets_Inclusive_SystematicsSplinesUp_results_LLFit__Pt_2_reg1_regFac1.00_noMuons_idSpectra%s.root"
             , path.Data(), jetPtString.Data());
      
      TString histTitles[numFiles];
      histTitles[0] = "Standard splines";
      histTitles[1] = "Splines - 0.2%";
      histTitles[2] = "Splines + 0.2%";
      
      SystematicErrorEstimation(path, outFileTitle, fileNames, histTitles, numFiles, nSigma, ignoreSigmaErrors, setMean, setMeanLowerThreshold,
                                setMeanUpperThreshold);
    }
    {
      Bool_t setMean = kFALSE;
      const Int_t numFiles = 3;
        
      TString outFileTitle = Form("TOFmode%s", jetPtString.Data());
      TString fileNames[numFiles];
      fileNames[0] = 
        Form("%s/bhess_PID_Jets_Inclusive_results_LLFit__Pt_2_reg1_regFac1.00_noMuons_idSpectra%s_TOFpatched.root"
            , path.Data(), jetPtString.Data());
      fileNames[1] =
        Form("%s/bhess_PID_Jets_Inclusive_SystematicsTOF0_results_LLFit__Pt_2_reg1_regFac1.00_noMuons_idSpectra%s_TOFpatched.root"
            , path.Data(), jetPtString.Data());
      fileNames[2] =
        Form("%s/bhess_PID_Jets_Inclusive_SystematicsTOF2_results_LLFit__Pt_2_reg1_regFac1.00_noMuons_idSpectra%s_TOFpatched.root"
            , path.Data(), jetPtString.Data());
      
      TString histTitles[numFiles];
      histTitles[0] = "Standard TOF mode 1";
      histTitles[1] = "TOF mode 0";
      histTitles[2] = "TOF mode 2";
      
      SystematicErrorEstimation(path, outFileTitle, fileNames, histTitles, numFiles, nSigma, ignoreSigmaErrors, setMean, setMeanLowerThreshold,
                                setMeanUpperThreshold);
    }
    
//     {
//       Bool_t setMean = kFALSE;
//       const Int_t numFiles = 2;//3;
//       
//       TString outFileTitle = Form("Muons%s", jetPtString.Data());
//       
//       TString fileNames[numFiles];
//       fileNames[0] =
//         Form("%s/bhess_PID_Jets_Inclusive_results_LLFit__Pt_2_reg1_regFac1.00_muonToElTunedOnMCHybridTrackCuts_idSpectra%s.root"
//              , path.Data(), jetPtString.Data());
//       fileNames[1] =
//         Form("%s/bhess_PID_Jets_Inclusive_results_LLFit__Pt_2_reg1_regFac1.00_muonsEqualElectrons_idSpectra%s.root"
//              , path.Data(), jetPtString.Data());
//       //fileNames[2] =
//         //Form("%s/bhess_PID_Jets_Inclusive_results_LLFit__Pt_2_reg1_regFac1.00_noMuons_idSpectra%s.root"
//         //     , path.Data(), jetPtString.Data());
//       
//       TString histTitles[numFiles];
//       histTitles[0] = "Muon_Electron_Ratio_Tuned_To_MC";
//       histTitles[1] = "Muon_Electron_Ratio_Is_Unity";
//       //histTitles[2] = "NoMuons";
// 
//       SystematicErrorEstimation(path, outFileTitle, fileNames, histTitles, numFiles, nSigma, ignoreSigmaErrors, setMean, setMeanLowerThreshold,
//                                setMeanUpperThreshold);  
//     }
    
    if (!useJetPt)
      break;
  }
  */
  
  // pp Jets from direct fit
  /*
  //const TString path = "/hera/alice/bhess/analysis/10d_e.pass2_merged/nclCut/differentTrackFilterBits/528";
  //const TString path = "/hera/alice/bhess/analysis/10d_e.pass2_merged/nclCut";
  //const TString path = "/hera/alice/bhess/analysis/10f6a/nclCut/noMCidForGen";
  const TString path = "/hera/alice/bhess/analysis/10f6a/nclCut/noMCidForGen/filterBit528";
  //const TString path = "/hera/alice/bhess/analysis/10f6a/nclCut/noMCidForGen/filterBit528/closureCheck";
  Bool_t useCentralities = kFALSE;
  Bool_t useJetPt = kTRUE;
  
  TString modeString[3] = { "Pt", "Z", "Xi" };
  for (Int_t mode = 0; mode <= 2; mode++) {
    TString muonString = "muonToElTunedOnMCHybridTrackCutsJets";
    // if (mode >= 1)
      muonString = "noMuons";
    
    for (Int_t iJetPt = 0; iJetPt <= numJetPtBins; iJetPt++) {
      if (iJetPt == numJetPtBins) {
        // 10-40 bin
        iJetPt++;
      }
      jetPtString = useJetPt ? Form("_jetPt%d.0_%d.0", jetPt[iJetPt], jetPt[iJetPt + 1]) : "";
      jetPtString.Append(chargeString);
      
      setMeanThreshold(setMeanLowerThreshold, setMeanUpperThreshold, mode, jetPt[iJetPt], jetPt[iJetPt + 1]);
      
      {
        Bool_t setMean = kFALSE;
        const Int_t numFiles = 2;
        
        TString outFileTitle = Form("PrePID%s", jetPtString.Data());
        
        TString fileNames[numFiles];
        fileNames[0] = 
          Form("%s/bhess_PID_Jets_results_LLFit__%s_2_reg1_regFac1.00_%s_idSpectra%s.root"
              , path.Data(), modeString[mode].Data(), muonString.Data(), jetPtString.Data());
        fileNames[1] = 
          Form("%s/bhess_PID_Jets_results_LLFit__%s_2_reg1_regFac1.00_%s%s.root"
              , path.Data(), modeString[mode].Data(), muonString.Data(), jetPtString.Data());
        
        TString histTitles[numFiles];
        histTitles[0] = "PID combined (default)";
        histTitles[1] = "No PID";
        
        outFileTitle = Form("%s_%s", modeString[mode].Data(), outFileTitle.Data());
        SystematicErrorEstimation(path, outFileTitle, fileNames, histTitles, numFiles, nSigma, ignoreSigmaErrors, setMean, 
                                  setMeanLowerThreshold, setMeanUpperThreshold);  
      }
      {
        Bool_t setMean = kFALSE;
        const Int_t numFiles = 2;
        
        TString outFileTitle = Form("Shape%s", jetPtString.Data());
        
        TString fileNames[numFiles];
        fileNames[0] = 
          Form("%s/bhess_PID_Jets_results_LLFit__%s_2_reg1_regFac1.00_%s_idSpectra%s.root"
              , path.Data(), modeString[mode].Data(), muonString.Data(), jetPtString.Data());
        fileNames[1] = 
          Form("%s/bhess_PID_Jets_PureGauss_results_LLFit__%s_2_reg1_regFac1.00_%s_idSpectra%s.root"
              , path.Data(), modeString[mode].Data(), muonString.Data(), jetPtString.Data());
        
          
        TString histTitles[numFiles];
        histTitles[0] = "Asymmetric shape (default)";
        histTitles[1] = "Pure gauss shape ";
        
        outFileTitle = Form("%s_%s", modeString[mode].Data(), outFileTitle.Data());
        SystematicErrorEstimation(path, outFileTitle, fileNames, histTitles, numFiles, nSigma, ignoreSigmaErrors, setMean, 
                                  setMeanLowerThreshold, setMeanUpperThreshold);  
      }
      {
        Bool_t setMean = kFALSE;
        const Int_t numFiles = 3;
        
        TString outFileTitle = Form("Sigma%s", jetPtString.Data());
        
        TString fileNames[numFiles];
        fileNames[0] = 
          Form("%s/bhess_PID_Jets_results_LLFit__%s_2_reg1_regFac1.00_%s_idSpectra%s.root"
              , path.Data(), modeString[mode].Data(), muonString.Data(), jetPtString.Data());
        fileNames[1] = 
          Form("%s/bhess_PID_Jets_SystematicsSigmaDown_results_LLFit__%s_2_reg1_regFac1.00_%s_idSpectra%s.root"
              , path.Data(), modeString[mode].Data(), muonString.Data(), jetPtString.Data());
        fileNames[2] = 
          Form("%s/bhess_PID_Jets_SystematicsSigmaUp_results_LLFit__%s_2_reg1_regFac1.00_%s_idSpectra%s.root"
              , path.Data(), modeString[mode].Data(), muonString.Data(), jetPtString.Data());
        
        TString histTitles[numFiles];
        histTitles[0] = "Standard #sigma map";
        histTitles[1] = "#sigma correction - 5%";
        histTitles[2] = "#sigma correction + 5%";
        
        outFileTitle = Form("%s_%s", modeString[mode].Data(), outFileTitle.Data());
        SystematicErrorEstimation(path, outFileTitle, fileNames, histTitles, numFiles, nSigma, ignoreSigmaErrors, setMean, 
                                  setMeanLowerThreshold, setMeanUpperThreshold);  
      }
      {
        Bool_t setMean = kFALSE;
        const Int_t numFiles = 3;
        
        TString outFileTitle = Form("Eta%s", jetPtString.Data());
        
        TString fileNames[numFiles];
        fileNames[0] = 
          Form("%s/bhess_PID_Jets_results_LLFit__%s_2_reg1_regFac1.00_%s_idSpectra%s.root"
              , path.Data(), modeString[mode].Data(), muonString.Data(), jetPtString.Data());
        fileNames[1] = 
          Form("%s/bhess_PID_Jets_SystematicsEtaDown_results_LLFit__%s_2_reg1_regFac1.00_%s_idSpectra%s.root"
              , path.Data(), modeString[mode].Data(), muonString.Data(), jetPtString.Data());
        fileNames[2] = 
          Form("%s/bhess_PID_Jets_SystematicsEtaUp_results_LLFit__%s_2_reg1_regFac1.00_%s_idSpectra%s.root"
              , path.Data(), modeString[mode].Data(), muonString.Data(), jetPtString.Data());
        
        TString histTitles[numFiles];
        histTitles[0] = "Standard #eta map";
        histTitles[1] = "#eta correction - 1%, p > 0.45 GeV/c; - %4, p < 0.45 GeV/c";
        histTitles[2] = "#eta correction + 1%, p > 0.45 GeV/c; + %4, p < 0.45 GeV/c";

        outFileTitle = Form("%s_%s", modeString[mode].Data(), outFileTitle.Data());
        SystematicErrorEstimation(path, outFileTitle, fileNames, histTitles, numFiles, nSigma, ignoreSigmaErrors, setMean, 
                                  setMeanLowerThreshold, setMeanUpperThreshold);  
      }
      {
        Bool_t setMean = kFALSE;
        const Int_t numFiles = 3;
          
        TString outFileTitle = Form("Splines%s", jetPtString.Data());
        TString fileNames[numFiles];
        fileNames[0] = 
          Form("%s/bhess_PID_Jets_results_LLFit__%s_2_reg1_regFac1.00_%s_idSpectra%s.root"
              , path.Data(), modeString[mode].Data(), muonString.Data(), jetPtString.Data());
        fileNames[1] =
          Form("%s/bhess_PID_Jets_SystematicsSplinesDown_results_LLFit__%s_2_reg1_regFac1.00_%s_idSpectra%s.root"
              , path.Data(), modeString[mode].Data(), muonString.Data(), jetPtString.Data());
        fileNames[2] =
          Form("%s/bhess_PID_Jets_SystematicsSplinesUp_results_LLFit__%s_2_reg1_regFac1.00_%s_idSpectra%s.root"
              , path.Data(), modeString[mode].Data(), muonString.Data(), jetPtString.Data());
        
        TString histTitles[numFiles];
        histTitles[0] = "Standard splines";
        histTitles[1] = "Splines - 0.2%";
        histTitles[2] = "Splines + 0.2%";
        
        outFileTitle = Form("%s_%s", modeString[mode].Data(), outFileTitle.Data());
        SystematicErrorEstimation(path, outFileTitle, fileNames, histTitles, numFiles, nSigma, ignoreSigmaErrors, setMean, 
                                  setMeanLowerThreshold, setMeanUpperThreshold);
      }
      {
        Bool_t setMean = kFALSE;
        const Int_t numFiles = 3;
          
        TString outFileTitle = Form("TOFmode%s", jetPtString.Data());
        TString fileNames[numFiles];
        fileNames[0] = 
          Form("%s/bhess_PID_Jets_results_LLFit__%s_2_reg1_regFac1.00_%s_idSpectra%s_TOFpatched.root"
              , path.Data(), modeString[mode].Data(), muonString.Data(), jetPtString.Data());
        fileNames[1] =
          Form("%s/bhess_PID_Jets_SystematicsTOF0_results_LLFit__%s_2_reg1_regFac1.00_%s_idSpectra%s_TOFpatched.root"
              , path.Data(), modeString[mode].Data(), muonString.Data(), jetPtString.Data());
        fileNames[2] =
          Form("%s/bhess_PID_Jets_SystematicsTOF2_results_LLFit__%s_2_reg1_regFac1.00_%s_idSpectra%s_TOFpatched.root"
              , path.Data(), modeString[mode].Data(), muonString.Data(), jetPtString.Data());
        
        TString histTitles[numFiles];
        histTitles[0] = "Standard TOF mode 1";
        histTitles[1] = "TOF mode 0";
        histTitles[2] = "TOF mode 2";
        
        outFileTitle = Form("%s_%s", modeString[mode].Data(), outFileTitle.Data());
        SystematicErrorEstimation(path, outFileTitle, fileNames, histTitles, numFiles, nSigma, ignoreSigmaErrors, setMean, 
                                  setMeanLowerThreshold, setMeanUpperThreshold);
      }
      
//      // TODO: NO muons yet for direct fit vs. z,xi
//       if (!useJetPt || mode == 0) {
//         Bool_t setMean = kFALSE;
//         const Int_t numFiles = 2;//3;
//         
//         TString outFileTitle = Form("Muons%s", jetPtString.Data());
//         
//         TString fileNames[numFiles];
//         fileNames[0] =
//           Form("%s/bhess_PID_Jets_results_LLFit__%s_2_reg1_regFac1.00_muonToElTunedOnMCHybridTrackCutsJets_idSpectra%s.root"
//               , path.Data(), modeString[mode].Data(), jetPtString.Data());
//         fileNames[1] =
//           Form("%s/bhess_PID_Jets_results_LLFit__%s_2_reg1_regFac1.00_muonsEqualElectrons_idSpectra%s.root"
//               , path.Data(), modeString[mode].Data(), jetPtString.Data());
//         //fileNames[2] =
//           //Form("%s/bhess_PID_Jets_results_LLFit__%s_2_reg1_regFac1.00_noMuons_idSpectra%s.root"
//           //     , path.Data(), modeString[mode].Data(), jetPtString.Data());
//         
//         TString histTitles[numFiles];
//         histTitles[0] = "Muon_Electron_Ratio_Tuned_To_MC";
//         histTitles[1] = "Muon_Electron_Ratio_Is_Unity";
//         //histTitles[2] = "NoMuons";
// 
//         outFileTitle = Form("%s_%s", modeString[mode].Data(), outFileTitle.Data());
//         SystematicErrorEstimation(path, outFileTitle, fileNames, histTitles, numFiles, nSigma, ignoreSigmaErrors, setMean, 
//                                  setMeanLowerThreshold, setMeanUpperThreshold);  
//       }
      
      if (!useJetPt)
        break;
    }
  }
  */
  
  
  // RefMult - pp Inclusive, no muons
  ///*
  //TString path = "/hera/alice/bhess/analysis/10d_e.pass2_merged/nclCut/multiplicityDependencePP/trackletsMult/pileUpRejectionSPD";
  TString path = "/hera/alice/bhess/analysis/10d_e.pass2_merged/nclCut/multiplicityDependencePP/V0Mmult/pileUpRejectionSPD/10d_e";
  
  //path = path.Append("/results_divided_by_MB");
  //path = path.Append("/MB_with_highPtBinningAsForJets");
  
  const Bool_t multRef = kFALSE;
  const Bool_t newRefBinning = kTRUE;
  
  Int_t numCentralitiesRefMult = newRefBinning ? numCentralitiesRefMultNew : numCentralitiesRefMultOld;
  Int_t* centralitiesRefMultLowEdge = 0x0;
  Int_t* centralitiesRefMultUpEdge = 0x0;
  
  if (multRef) {
    centralitiesRefMultLowEdge = new Int_t[numCentralitiesRefMult];
    centralitiesRefMultUpEdge = new Int_t[numCentralitiesRefMult];
    
    for (Int_t i = 0; i < numCentralitiesRefMult; i++) {
      centralitiesRefMultLowEdge[i] = newRefBinning ? centralitiesRefMultNewLowEdge[i] : centralitiesRefMultOldLowEdge[i];
      centralitiesRefMultUpEdge[i]  = newRefBinning ? centralitiesRefMultNewUpEdge[i]  : centralitiesRefMultOldUpEdge[i];
    }
  }
  
  // Inclusive means: Only pT
  setMeanThreshold(setMeanLowerThreshold, setMeanUpperThreshold, kPt, -1, -1);
  
  const Int_t numMults = multRef ? numCentralitiesRefMult : numCentralitiesV0Mpp;
  for (Int_t iMult = -1; iMult < numMults; iMult++) {
    if (multRef) {
      if (iMult == -1)
        centralityString = "";
      else
        centralityString = Form("_centrality%d_%d", centralitiesRefMultLowEdge[iMult], centralitiesRefMultUpEdge[iMult]);
    }
    else {
      centralityString = "";
      if (iMult >= 0) {
        Bool_t centralityHasDecimalsPlaces = CentralityHasDecimalsPlaces(centralitiesV0MppLowEdge[iMult]) ||
                                             CentralityHasDecimalsPlaces(centralitiesV0MppUpEdge[iMult]);
        
        centralityString = Form(centralityHasDecimalsPlaces ? "_centrality%.0fem2_%.0fem2" : "_centrality%.0f_%.0f", 
                                centralityHasDecimalsPlaces ? 100.*centralitiesV0MppLowEdge[iMult] : centralitiesV0MppLowEdge[iMult],
                                centralityHasDecimalsPlaces ? 100.*centralitiesV0MppUpEdge[iMult] : centralitiesV0MppUpEdge[iMult]);
      }
    }
    centralityString.Append(chargeString);
    
    std::cout << "iMult " << iMult << ": \"" << centralityString.Data() << "\"..." << std::endl;
    {
      Bool_t setMean = kFALSE;
      const Int_t numFiles = 2;
      
      TString outFileTitle = Form("PrePID%s", centralityString.Data());
      
      TString fileNames[numFiles];
      fileNames[0] = 
        Form("%s/bhess_PID_TPCdefaultPriors_ITS_results_LLFit__Pt_2_reg1_regFac1.00_noMuons_idSpectra%s.root"
             , path.Data(), centralityString.Data());
      fileNames[1] = 
        Form("%s/bhess_PID_TPCdefaultPriors_ITS_results_LLFit__Pt_2_reg1_regFac1.00_noMuons%s.root"
             , path.Data(), centralityString.Data());
      
      TString histTitles[numFiles];
      histTitles[0] = "PID combined (default)";
      histTitles[1] = "No PID";
      
      SystematicErrorEstimation(path, outFileTitle, fileNames, histTitles, numFiles, nSigma, ignoreSigmaErrors, setMean, setMeanLowerThreshold,
                                setMeanUpperThreshold);  
    }
    {
      Bool_t setMean = kFALSE;
      const Int_t numFiles = 2;
      
      TString outFileTitle = Form("Shape%s", centralityString.Data());
      
      TString fileNames[numFiles];
      fileNames[0] = 
        Form("%s/bhess_PID_TPCdefaultPriors_ITS_results_LLFit__Pt_2_reg1_regFac1.00_noMuons_idSpectra%s.root"
             , path.Data(), centralityString.Data());
      fileNames[1] = 
        Form("%s/bhess_PID_TPCdefaultPriors_ITS_PureGauss_results_LLFit__Pt_2_reg1_regFac1.00_noMuons_idSpectra%s.root"
             , path.Data(), centralityString.Data());
      
        
      TString histTitles[numFiles];
      histTitles[0] = "Asymmetric shape (default)";
      histTitles[1] = "Pure gauss shape ";
      
      SystematicErrorEstimation(path, outFileTitle, fileNames, histTitles, numFiles, nSigma, ignoreSigmaErrors, setMean, setMeanLowerThreshold,
                                setMeanUpperThreshold);  
    }
    {
      Bool_t setMean = kFALSE;
      const Int_t numFiles = 3;
      
      TString outFileTitle = Form("Sigma%s", centralityString.Data());
      
      TString fileNames[numFiles];
      fileNames[0] = 
        Form("%s/bhess_PID_TPCdefaultPriors_ITS_results_LLFit__Pt_2_reg1_regFac1.00_noMuons_idSpectra%s.root"
             , path.Data(), centralityString.Data());
      fileNames[1] = 
        Form("%s/bhess_PID_TPCdefaultPriors_ITS_SystematicsSigmaDown_results_LLFit__Pt_2_reg1_regFac1.00_noMuons_idSpectra%s.root"
             , path.Data(), centralityString.Data());
      fileNames[2] = 
        Form("%s/bhess_PID_TPCdefaultPriors_ITS_SystematicsSigmaUp_results_LLFit__Pt_2_reg1_regFac1.00_noMuons_idSpectra%s.root"
             , path.Data(), centralityString.Data());
      
      TString histTitles[numFiles];
      histTitles[0] = "Standard #sigma map";
      histTitles[1] = "#sigma correction - 5%";
      histTitles[2] = "#sigma correction + 5%";
      
      SystematicErrorEstimation(path, outFileTitle, fileNames, histTitles, numFiles, nSigma, ignoreSigmaErrors, setMean, setMeanLowerThreshold,
                                setMeanUpperThreshold);  
    }
    {
      Bool_t setMean = kFALSE;
      const Int_t numFiles = 3;
      
      TString outFileTitle = Form("Eta%s", centralityString.Data());
      
      TString fileNames[numFiles];
      fileNames[0] = 
        Form("%s/bhess_PID_TPCdefaultPriors_ITS_results_LLFit__Pt_2_reg1_regFac1.00_noMuons_idSpectra%s.root"
             , path.Data(), centralityString.Data());
      fileNames[1] = 
        Form("%s/bhess_PID_TPCdefaultPriors_ITS_SystematicsEtaDown_results_LLFit__Pt_2_reg1_regFac1.00_noMuons_idSpectra%s.root"
             , path.Data(), centralityString.Data());
      fileNames[2] = 
        Form("%s/bhess_PID_TPCdefaultPriors_ITS_SystematicsEtaUp_results_LLFit__Pt_2_reg1_regFac1.00_noMuons_idSpectra%s.root"
             , path.Data(), centralityString.Data());
      
      TString histTitles[numFiles];
      histTitles[0] = "Standard #eta map";
      histTitles[1] = "#eta correction - 1%, p > 0.45 GeV/c; - %4, p < 0.45 GeV/c";
      histTitles[2] = "#eta correction + 1%, p > 0.45 GeV/c; + %4, p < 0.45 GeV/c";

      SystematicErrorEstimation(path, outFileTitle, fileNames, histTitles, numFiles, nSigma, ignoreSigmaErrors, setMean, setMeanLowerThreshold,
                                setMeanUpperThreshold);  
    }
    {
      Bool_t setMean = kFALSE;
      const Int_t numFiles = 3;
        
      TString outFileTitle = Form("Splines%s", centralityString.Data());
      TString fileNames[numFiles];
      fileNames[0] = 
        Form("%s/bhess_PID_TPCdefaultPriors_ITS_results_LLFit__Pt_2_reg1_regFac1.00_noMuons_idSpectra%s.root"
             , path.Data(), centralityString.Data());
      fileNames[1] =
        Form("%s/bhess_PID_TPCdefaultPriors_ITS_SystematicsSplinesDown_results_LLFit__Pt_2_reg1_regFac1.00_noMuons_idSpectra%s.root"
             , path.Data(), centralityString.Data());
      fileNames[2] =
        Form("%s/bhess_PID_TPCdefaultPriors_ITS_SystematicsSplinesUp_results_LLFit__Pt_2_reg1_regFac1.00_noMuons_idSpectra%s.root"
             , path.Data(), centralityString.Data());
      
      TString histTitles[numFiles];
      histTitles[0] = "Standard splines";
      histTitles[1] = "Splines - 0.2%";
      histTitles[2] = "Splines + 0.2%";
      
      SystematicErrorEstimation(path, outFileTitle, fileNames, histTitles, numFiles, nSigma, ignoreSigmaErrors, setMean, setMeanLowerThreshold,
                                setMeanUpperThreshold);
    }
//     {
//       Bool_t setMean = kFALSE;
//       const Int_t numFiles = 3;
//         
//       TString outFileTitle = Form("TOFmode%s", centralityString.Data());
//       TString fileNames[numFiles];
//       fileNames[0] = 
//         Form("%s/bhess_PID_TPCdefaultPriors_ITS_results_LLFit__Pt_2_reg1_regFac1.00_noMuons_idSpectra%s_TOFpatched.root"
//             , path.Data(), centralityString.Data());
//       fileNames[1] =
//         Form("%s/bhess_PID_TPCdefaultPriors_ITS_SystematicsTOF0_results_LLFit__Pt_2_reg1_regFac1.00_noMuons_idSpectra%s_TOFpatched.root"
//             , path.Data(), centralityString.Data());
//       fileNames[2] =
//         Form("%s/bhess_PID_TPCdefaultPriors_ITS_SystematicsTOF2_results_LLFit__Pt_2_reg1_regFac1.00_noMuons_idSpectra%s_TOFpatched.root"
//             , path.Data(), centralityString.Data());
//       
//       TString histTitles[numFiles];
//       histTitles[0] = "Standard TOF mode 1";
//       histTitles[1] = "TOF mode 0";
//       histTitles[2] = "TOF mode 2";
//       
//       SystematicErrorEstimation(path, outFileTitle, fileNames, histTitles, numFiles, nSigma, ignoreSigmaErrors, setMean, setMeanLowerThreshold,
//                                 setMeanUpperThreshold);
//     }
  }
  //*/
  
  
  // pp Inclusive with muons (old)
  /*
  const TString path = "/hera/alice/bhess/analysis/10d_e.pass2_merged/nclCut";
  //const TString path = "/hera/alice/bhess/analysis/10d.pass2/nclCut/";
  //const TString path = "/hera/alice/bhess/analysis/10f6a/nclCut/noMCidForGen";
  Bool_t useCentralities = kFALSE;
  Bool_t useJetPt = kFALSE;
  
  // Inclusive means: Only pT
  setMeanThreshold(setMeanLowerThreshold, setMeanUpperThreshold, kPt, -1, -1);

  for (Int_t iJetPt = 0; iJetPt <= numJetPtBins; iJetPt++) {
    if (iJetPt == numJetPtBins) {
      // 10-40 bin
      iJetPt++;
    }
    jetPtString = useJetPt ? Form("_jetPt%d.0_%d.0", jetPt[iJetPt], jetPt[iJetPt + 1]) : "";
    
    {
      Bool_t setMean = kFALSE;
      const Int_t numFiles = 2;
      
      TString outFileTitle = Form("PrePID%s", jetPtString.Data());
      
      TString fileNames[numFiles];
      fileNames[0] = 
        Form("%s/bhess_PID_Jets_Inclusive_results_LLFit__Pt_2_reg1_regFac1.00_muonToElTunedOnMCHybridTrackCuts_idSpectra%s.root"
             , path.Data(), jetPtString.Data());
      fileNames[1] = 
        Form("%s/bhess_PID_Jets_Inclusive_results_LLFit__Pt_2_reg1_regFac1.00_muonToElTunedOnMCHybridTrackCuts%s.root"
             , path.Data(), jetPtString.Data());
      
      TString histTitles[numFiles];
      histTitles[0] = "PID combined (default)";
      histTitles[1] = "No PID";
      
      SystematicErrorEstimation(path, outFileTitle, fileNames, histTitles, numFiles, nSigma, ignoreSigmaErrors, setMean, 
                                setMeanLowerThreshold, setMeanUpperThreshold);  
    }
    {
      Bool_t setMean = kFALSE;
      const Int_t numFiles = 2;
      
      TString outFileTitle = Form("Shape%s", jetPtString.Data());
      
      TString fileNames[numFiles];
      fileNames[0] = 
        Form("%s/bhess_PID_Jets_Inclusive_results_LLFit__Pt_2_reg1_regFac1.00_muonToElTunedOnMCHybridTrackCuts_idSpectra%s.root"
             , path.Data(), jetPtString.Data());
      fileNames[1] = 
        Form("%s/bhess_PID_Jets_Inclusive_PureGauss_results_LLFit__Pt_2_reg1_regFac1.00_muonToElTunedOnMCHybridTrackCuts_idSpectra%s.root"
             , path.Data(), jetPtString.Data());
      
        
      TString histTitles[numFiles];
      histTitles[0] = "Asymmetric shape (default)";
      histTitles[1] = "Pure gauss shape ";
      
      SystematicErrorEstimation(path, outFileTitle, fileNames, histTitles, numFiles, nSigma, ignoreSigmaErrors, setMean, 
                                setMeanLowerThreshold, setMeanUpperThreshold);  
    }
    {
      Bool_t setMean = kFALSE;
      const Int_t numFiles = 3;
      
      TString outFileTitle = Form("Sigma%s", jetPtString.Data());
      
      TString fileNames[numFiles];
      fileNames[0] = 
        Form("%s/bhess_PID_Jets_Inclusive_results_LLFit__Pt_2_reg1_regFac1.00_muonToElTunedOnMCHybridTrackCuts_idSpectra%s.root"
             , path.Data(), jetPtString.Data());
      fileNames[1] = 
        Form("%s/bhess_PID_Jets_Inclusive_SystematicsSigmaDown_results_LLFit__Pt_2_reg1_regFac1.00_muonToElTunedOnMCHybridTrackCuts_idSpectra%s.root"
             , path.Data(), jetPtString.Data());
      fileNames[2] = 
        Form("%s/bhess_PID_Jets_Inclusive_SystematicsSigmaUp_results_LLFit__Pt_2_reg1_regFac1.00_muonToElTunedOnMCHybridTrackCuts_idSpectra%s.root"
             , path.Data(), jetPtString.Data());
      
      TString histTitles[numFiles];
      histTitles[0] = "Standard #sigma map";
      histTitles[1] = "#sigma correction - 5%";
      histTitles[2] = "#sigma correction + 5%";
      
      SystematicErrorEstimation(path, outFileTitle, fileNames, histTitles, numFiles, nSigma, ignoreSigmaErrors, setMean, 
                                setMeanLowerThreshold, setMeanUpperThreshold);  
    }
    {
      Bool_t setMean = kFALSE;
      const Int_t numFiles = 3;
      
      TString outFileTitle = Form("Eta%s", jetPtString.Data());
      
      TString fileNames[numFiles];
      fileNames[0] = 
        Form("%s/bhess_PID_Jets_Inclusive_results_LLFit__Pt_2_reg1_regFac1.00_muonToElTunedOnMCHybridTrackCuts_idSpectra%s.root"
             , path.Data(), jetPtString.Data());
      fileNames[1] = 
        Form("%s/bhess_PID_Jets_Inclusive_SystematicsEtaDown_results_LLFit__Pt_2_reg1_regFac1.00_muonToElTunedOnMCHybridTrackCuts_idSpectra%s.root"
             , path.Data(), jetPtString.Data());
      fileNames[2] = 
        Form("%s/bhess_PID_Jets_Inclusive_SystematicsEtaUp_results_LLFit__Pt_2_reg1_regFac1.00_muonToElTunedOnMCHybridTrackCuts_idSpectra%s.root"
             , path.Data(), jetPtString.Data());
      
      TString histTitles[numFiles];
      histTitles[0] = "Standard #eta map";
      histTitles[1] = "#eta correction - 1%, p > 0.45 GeV/c; - %4, p < 0.45 GeV/c";
      histTitles[2] = "#eta correction + 1%, p > 0.45 GeV/c; + %4, p < 0.45 GeV/c";

      SystematicErrorEstimation(path, outFileTitle, fileNames, histTitles, numFiles, nSigma, ignoreSigmaErrors, setMean, 
                                setMeanLowerThreshold, setMeanUpperThreshold);  
    }
    {
      Bool_t setMean = kFALSE;
      const Int_t numFiles = 3;
        
      TString outFileTitle = Form("Splines%s", jetPtString.Data());
      TString fileNames[numFiles];
      fileNames[0] = 
        Form("%s/bhess_PID_Jets_Inclusive_results_LLFit__Pt_2_reg1_regFac1.00_muonToElTunedOnMCHybridTrackCuts_idSpectra%s.root"
             , path.Data(), jetPtString.Data());
      fileNames[1] =
        Form("%s/bhess_PID_Jets_Inclusive_SystematicsSplinesDown_results_LLFit__Pt_2_reg1_regFac1.00_muonToElTunedOnMCHybridTrackCuts_idSpectra%s.root"
             , path.Data(), jetPtString.Data());
      fileNames[2] =
        Form("%s/bhess_PID_Jets_Inclusive_SystematicsSplinesUp_results_LLFit__Pt_2_reg1_regFac1.00_muonToElTunedOnMCHybridTrackCuts_idSpectra%s.root"
             , path.Data(), jetPtString.Data());
      
      TString histTitles[numFiles];
      histTitles[0] = "Standard splines";
      histTitles[1] = "Splines - 0.2%";
      histTitles[2] = "Splines + 0.2%";
      
      SystematicErrorEstimation(path, outFileTitle, fileNames, histTitles, numFiles, nSigma, ignoreSigmaErrors, setMean, 
                                setMeanLowerThreshold, setMeanUpperThreshold);
    }
    {
      Bool_t setMean = kFALSE;
      const Int_t numFiles = 3;
        
      TString outFileTitle = Form("TOFmode%s", jetPtString.Data());
      TString fileNames[numFiles];
      fileNames[0] = 
        Form("%s/bhess_PID_Jets_Inclusive_results_LLFit__Pt_2_reg1_regFac1.00_muonToElTunedOnMCHybridTrackCuts_idSpectra%s_TOFpatched.root"
            , path.Data(), jetPtString.Data());
      fileNames[1] =
        Form("%s/bhess_PID_Jets_Inclusive_SystematicsTOF0_results_LLFit__Pt_2_reg1_regFac1.00_muonToElTunedOnMCHybridTrackCuts_idSpectra%s_TOFpatched.root"
            , path.Data(), jetPtString.Data());
      fileNames[2] =
        Form("%s/bhess_PID_Jets_Inclusive_SystematicsTOF2_results_LLFit__Pt_2_reg1_regFac1.00_muonToElTunedOnMCHybridTrackCuts_idSpectra%s_TOFpatched.root"
            , path.Data(), jetPtString.Data());
      
      TString histTitles[numFiles];
      histTitles[0] = "Standard TOF mode 1";
      histTitles[1] = "TOF mode 0";
      histTitles[2] = "TOF mode 2";
      
      SystematicErrorEstimation(path, outFileTitle, fileNames, histTitles, numFiles, nSigma, ignoreSigmaErrors, setMean, 
                                setMeanLowerThreshold, setMeanUpperThreshold);
    }
    {
      Bool_t setMean = kFALSE;
      const Int_t numFiles = 2;//3;
      
      TString outFileTitle = Form("Muons%s", jetPtString.Data());
      
      TString fileNames[numFiles];
      fileNames[0] =
        Form("%s/bhess_PID_Jets_Inclusive_results_LLFit__Pt_2_reg1_regFac1.00_muonToElTunedOnMCHybridTrackCuts_idSpectra%s.root"
             , path.Data(), jetPtString.Data());
      fileNames[1] =
        Form("%s/bhess_PID_Jets_Inclusive_results_LLFit__Pt_2_reg1_regFac1.00_muonsEqualElectrons_idSpectra%s.root"
             , path.Data(), jetPtString.Data());
      //fileNames[2] =
        //Form("%s/bhess_PID_Jets_Inclusive_results_LLFit__Pt_2_reg1_regFac1.00_noMuons_idSpectra%s.root"
        //     , path.Data(), jetPtString.Data());
      
      TString histTitles[numFiles];
      histTitles[0] = "Muon_Electron_Ratio_Tuned_To_MC";
      histTitles[1] = "Muon_Electron_Ratio_Is_Unity";
      //histTitles[2] = "NoMuons";

      SystematicErrorEstimation(path, outFileTitle, fileNames, histTitles, numFiles, nSigma, ignoreSigmaErrors, setMean, 
                                setMeanLowerThreshold, setMeanUpperThreshold);  
    }
    
    if (!useJetPt)
      break;
  }
  */
  
  
  // pp Jets from extract FFs with muons (old)
  /*
  const TString path = "/hera/alice/bhess/analysis/10d_e.pass2_merged/nclCut";
  Bool_t useCentralities = kFALSE;
  Bool_t useJetPt = kTRUE;

  for (Int_t iJetPt = 0; iJetPt <= numJetPtBins; iJetPt++) {
    if (iJetPt == numJetPtBins) {
      // 10-40 bin
      iJetPt++;
    }
    jetPtString = useJetPt ? Form("_jetPt%d.0_%d.0__centrality_all_jetPt%d.0_%d.0", jetPt[iJetPt], jetPt[iJetPt + 1],
                                  jetPt[iJetPt], jetPt[iJetPt + 1])
                           : "";
    
    // Pt only so far
    setMeanThreshold(setMeanLowerThreshold, setMeanUpperThreshold, kPt, jetPt[iJetPt], jetPt[iJetPt + 1]);
    {
      Bool_t setMean = kFALSE;
      const Int_t numFiles = 2;
      
      TString outFileTitle = Form("PrePID%s", jetPtString.Data());
      
      TString fileNames[numFiles];
      fileNames[0] = 
        Form("%s/output_extractedFFs_bhess_PID_Jets_results_LLFit__Pt_2_reg1_regFac1.00_muonToElTunedOnMCHybridTrackCutsJets_idSpectra%s.root"
             , path.Data(), jetPtString.Data());
      fileNames[1] = 
        Form("%s/output_extractedFFs_bhess_PID_Jets_results_LLFit__Pt_2_reg1_regFac1.00_muonToElTunedOnMCHybridTrackCutsJets%s.root"
             , path.Data(), jetPtString.Data());
      
      TString histTitles[numFiles];
      histTitles[0] = "PID combined (default)";
      histTitles[1] = "No PID";
      
      for (Int_t mode = 0; mode <= 2; mode++)
        SystematicErrorEstimation(path, outFileTitle, fileNames, histTitles, numFiles, nSigma, ignoreSigmaErrors, setMean, 
                                  setMeanLowerThreshold, setMeanUpperThreshold, mode);  
    }
    {
      Bool_t setMean = kFALSE;
      const Int_t numFiles = 2;
      
      TString outFileTitle = Form("Shape%s", jetPtString.Data());
      
      TString fileNames[numFiles];
      fileNames[0] = 
        Form("%s/output_extractedFFs_bhess_PID_Jets_results_LLFit__Pt_2_reg1_regFac1.00_muonToElTunedOnMCHybridTrackCutsJets_idSpectra%s.root"
             , path.Data(), jetPtString.Data());
      fileNames[1] = 
        Form("%s/output_extractedFFs_bhess_PID_Jets_PureGauss_results_LLFit__Pt_2_reg1_regFac1.00_muonToElTunedOnMCHybridTrackCutsJets_idSpectra%s.root"
             , path.Data(), jetPtString.Data());
      
        
      TString histTitles[numFiles];
      histTitles[0] = "Asymmetric shape (default)";
      histTitles[1] = "Pure gauss shape ";
      
      for (Int_t mode = 0; mode <= 2; mode++)
        SystematicErrorEstimation(path, outFileTitle, fileNames, histTitles, numFiles, nSigma, ignoreSigmaErrors, setMean, 
                                  setMeanLowerThreshold, setMeanUpperThreshold, mode);  
    }
    {
      Bool_t setMean = kFALSE;
      const Int_t numFiles = 3;
      
      TString outFileTitle = Form("Sigma%s", jetPtString.Data());
      
      TString fileNames[numFiles];
      fileNames[0] = 
        Form("%s/output_extractedFFs_bhess_PID_Jets_results_LLFit__Pt_2_reg1_regFac1.00_muonToElTunedOnMCHybridTrackCutsJets_idSpectra%s.root"
             , path.Data(), jetPtString.Data());
      fileNames[1] = 
        Form("%s/output_extractedFFs_bhess_PID_Jets_SystematicsSigmaDown_results_LLFit__Pt_2_reg1_regFac1.00_muonToElTunedOnMCHybridTrackCutsJets_idSpectra%s.root"
             , path.Data(), jetPtString.Data());
      fileNames[2] = 
        Form("%s/output_extractedFFs_bhess_PID_Jets_SystematicsSigmaUp_results_LLFit__Pt_2_reg1_regFac1.00_muonToElTunedOnMCHybridTrackCutsJets_idSpectra%s.root"
             , path.Data(), jetPtString.Data());
      
      TString histTitles[numFiles];
      histTitles[0] = "Standard #sigma map";
      histTitles[1] = "#sigma correction - 5%";
      histTitles[2] = "#sigma correction + 5%";
      
      for (Int_t mode = 0; mode <= 2; mode++)
        SystematicErrorEstimation(path, outFileTitle, fileNames, histTitles, numFiles, nSigma, ignoreSigmaErrors, setMean, 
                                  setMeanLowerThreshold, setMeanUpperThreshold, mode);  
    }
    {
      Bool_t setMean = kFALSE;
      const Int_t numFiles = 3;
      
      TString outFileTitle = Form("Eta%s", jetPtString.Data());
      
      TString fileNames[numFiles];
      fileNames[0] = 
        Form("%s/output_extractedFFs_bhess_PID_Jets_results_LLFit__Pt_2_reg1_regFac1.00_muonToElTunedOnMCHybridTrackCutsJets_idSpectra%s.root"
             , path.Data(), jetPtString.Data());
      fileNames[1] = 
        Form("%s/output_extractedFFs_bhess_PID_Jets_SystematicsEtaDown_results_LLFit__Pt_2_reg1_regFac1.00_muonToElTunedOnMCHybridTrackCutsJets_idSpectra%s.root"
             , path.Data(), jetPtString.Data());
      fileNames[2] = 
        Form("%s/output_extractedFFs_bhess_PID_Jets_SystematicsEtaUp_results_LLFit__Pt_2_reg1_regFac1.00_muonToElTunedOnMCHybridTrackCutsJets_idSpectra%s.root"
             , path.Data(), jetPtString.Data());
      
      TString histTitles[numFiles];
      histTitles[0] = "Standard #eta map";
      histTitles[1] = "#eta correction - 1%, p > 0.45 GeV/c; - %4, p < 0.45 GeV/c";
      histTitles[2] = "#eta correction + 1%, p > 0.45 GeV/c; + %4, p < 0.45 GeV/c";

      for (Int_t mode = 0; mode <= 2; mode++)
        SystematicErrorEstimation(path, outFileTitle, fileNames, histTitles, numFiles, nSigma, ignoreSigmaErrors, setMean, 
                                  setMeanLowerThreshold, setMeanUpperThreshold, mode);  
    }
    {
      Bool_t setMean = kFALSE;
      const Int_t numFiles = 3;
        
      TString outFileTitle = Form("Splines%s", jetPtString.Data());
      TString fileNames[numFiles];
      fileNames[0] = 
        Form("%s/output_extractedFFs_bhess_PID_Jets_results_LLFit__Pt_2_reg1_regFac1.00_muonToElTunedOnMCHybridTrackCutsJets_idSpectra%s.root"
             , path.Data(), jetPtString.Data());
      fileNames[1] =
        Form("%s/output_extractedFFs_bhess_PID_Jets_SystematicsSplinesDown_results_LLFit__Pt_2_reg1_regFac1.00_muonToElTunedOnMCHybridTrackCutsJets_idSpectra%s.root"
             , path.Data(), jetPtString.Data());
      fileNames[2] =
        Form("%s/output_extractedFFs_bhess_PID_Jets_SystematicsSplinesUp_results_LLFit__Pt_2_reg1_regFac1.00_muonToElTunedOnMCHybridTrackCutsJets_idSpectra%s.root"
             , path.Data(), jetPtString.Data());
      
      TString histTitles[numFiles];
      histTitles[0] = "Standard splines";
      histTitles[1] = "Splines - 0.2%";
      histTitles[2] = "Splines + 0.2%";
      
      for (Int_t mode = 0; mode <= 2; mode++)
        SystematicErrorEstimation(path, outFileTitle, fileNames, histTitles, numFiles, nSigma, ignoreSigmaErrors, setMean, 
                                  setMeanLowerThreshold, setMeanUpperThreshold, mode);
    }
    {
      Bool_t setMean = kFALSE;
      const Int_t numFiles = 2;//3;
      
      TString outFileTitle = Form("Muons%s", jetPtString.Data());
      
      TString fileNames[numFiles];
      fileNames[0] =
        Form("%s/output_extractedFFs_bhess_PID_Jets_results_LLFit__Pt_2_reg1_regFac1.00_muonToElTunedOnMCHybridTrackCutsJets_idSpectra%s.root"
             , path.Data(), jetPtString.Data());
      fileNames[1] =
        Form("%s/output_extractedFFs_bhess_PID_Jets_results_LLFit__Pt_2_reg1_regFac1.00_muonsEqualElectrons_idSpectra%s.root"
             , path.Data(), jetPtString.Data());
      //fileNames[2] =
        //Form("%s/output_extractedFFs_bhess_PID_Jets_results_LLFit__Pt_2_reg1_regFac1.00_noMuons_idSpectra%s.root"
        //     , path.Data(), jetPtString.Data());
      
      TString histTitles[numFiles];
      histTitles[0] = "Muon_Electron_Ratio_Tuned_To_MC";
      histTitles[1] = "Muon_Electron_Ratio_Is_Unity";
      //histTitles[2] = "NoMuons";

      for (Int_t mode = 0; mode <= 2; mode++)
        SystematicErrorEstimation(path, outFileTitle, fileNames, histTitles, numFiles, nSigma, ignoreSigmaErrors, setMean, 
                                  setMeanLowerThreshold, setMeanUpperThreshold, mode);  
    }
    
    if (!useJetPt)
      break;
  }
  */

  // pp Jets from direct fit with muons (old)
  /*
  //TODO: NO Muons yet => no error on muons yet
  const TString path = "/hera/alice/bhess/analysis/10d_e.pass2_merged/nclCut";
  //const TString path = "/hera/alice/bhess/analysis/10f6a/nclCut/noMCidForGen";
  Bool_t useCentralities = kFALSE;
  Bool_t useJetPt = kTRUE;

  TString modeString[3] = { "Pt", "Z", "Xi" };
  for (Int_t mode = 0; mode <= 2; mode++) {
    TString muonString = "muonToElTunedOnMCHybridTrackCutsJets";
    // TODO: NO muons yet for direct fit
    if (mode >= 1)
      muonString = "noMuons";
    
    for (Int_t iJetPt = 0; iJetPt <= numJetPtBins; iJetPt++) {
      if (iJetPt == numJetPtBins) {
        // 10-40 bin
        iJetPt++;
      }
      jetPtString = useJetPt ? Form("_jetPt%d.0_%d.0", jetPt[iJetPt], jetPt[iJetPt + 1]) : "";
      
      setMeanThreshold(setMeanLowerThreshold, setMeanUpperThreshold, mode, jetPt[iJetPt], jetPt[iJetPt + 1]);
      
      {
        Bool_t setMean = kFALSE;
        const Int_t numFiles = 2;
        
        TString outFileTitle = Form("PrePID%s", jetPtString.Data());
        
        TString fileNames[numFiles];
        fileNames[0] = 
          Form("%s/bhess_PID_Jets_results_LLFit__%s_2_reg1_regFac1.00_%s_idSpectra%s.root"
              , path.Data(), modeString[mode].Data(), muonString.Data(), jetPtString.Data());
        fileNames[1] = 
          Form("%s/bhess_PID_Jets_results_LLFit__%s_2_reg1_regFac1.00_%s%s.root"
              , path.Data(), modeString[mode].Data(), muonString.Data(), jetPtString.Data());
        
        TString histTitles[numFiles];
        histTitles[0] = "PID combined (default)";
        histTitles[1] = "No PID";
        
        outFileTitle = Form("%s_%s", modeString[mode].Data(), outFileTitle.Data());
        SystematicErrorEstimation(path, outFileTitle, fileNames, histTitles, numFiles, nSigma, ignoreSigmaErrors, setMean, 
                                  setMeanLowerThreshold, setMeanUpperThreshold);  
      }
      {
        Bool_t setMean = kFALSE;
        const Int_t numFiles = 2;
        
        TString outFileTitle = Form("Shape%s", jetPtString.Data());
        
        TString fileNames[numFiles];
        fileNames[0] = 
          Form("%s/bhess_PID_Jets_results_LLFit__%s_2_reg1_regFac1.00_%s_idSpectra%s.root"
              , path.Data(), modeString[mode].Data(), muonString.Data(), jetPtString.Data());
        fileNames[1] = 
          Form("%s/bhess_PID_Jets_PureGauss_results_LLFit__%s_2_reg1_regFac1.00_%s_idSpectra%s.root"
              , path.Data(), modeString[mode].Data(), muonString.Data(), jetPtString.Data());
        
          
        TString histTitles[numFiles];
        histTitles[0] = "Asymmetric shape (default)";
        histTitles[1] = "Pure gauss shape ";
        
        outFileTitle = Form("%s_%s", modeString[mode].Data(), outFileTitle.Data());
        SystematicErrorEstimation(path, outFileTitle, fileNames, histTitles, numFiles, nSigma, ignoreSigmaErrors, setMean, 
                                  setMeanLowerThreshold, setMeanUpperThreshold);  
      }
      {
        Bool_t setMean = kFALSE;
        const Int_t numFiles = 3;
        
        TString outFileTitle = Form("Sigma%s", jetPtString.Data());
        
        TString fileNames[numFiles];
        fileNames[0] = 
          Form("%s/bhess_PID_Jets_results_LLFit__%s_2_reg1_regFac1.00_%s_idSpectra%s.root"
              , path.Data(), modeString[mode].Data(), muonString.Data(), jetPtString.Data());
        fileNames[1] = 
          Form("%s/bhess_PID_Jets_SystematicsSigmaDown_results_LLFit__%s_2_reg1_regFac1.00_%s_idSpectra%s.root"
              , path.Data(), modeString[mode].Data(), muonString.Data(), jetPtString.Data());
        fileNames[2] = 
          Form("%s/bhess_PID_Jets_SystematicsSigmaUp_results_LLFit__%s_2_reg1_regFac1.00_%s_idSpectra%s.root"
              , path.Data(), modeString[mode].Data(), muonString.Data(), jetPtString.Data());
        
        TString histTitles[numFiles];
        histTitles[0] = "Standard #sigma map";
        histTitles[1] = "#sigma correction - 5%";
        histTitles[2] = "#sigma correction + 5%";
        
        outFileTitle = Form("%s_%s", modeString[mode].Data(), outFileTitle.Data());
        SystematicErrorEstimation(path, outFileTitle, fileNames, histTitles, numFiles, nSigma, ignoreSigmaErrors, setMean, 
                                  setMeanLowerThreshold, setMeanUpperThreshold);  
      }
      {
        Bool_t setMean = kFALSE;
        const Int_t numFiles = 3;
        
        TString outFileTitle = Form("Eta%s", jetPtString.Data());
        
        TString fileNames[numFiles];
        fileNames[0] = 
          Form("%s/bhess_PID_Jets_results_LLFit__%s_2_reg1_regFac1.00_%s_idSpectra%s.root"
              , path.Data(), modeString[mode].Data(), muonString.Data(), jetPtString.Data());
        fileNames[1] = 
          Form("%s/bhess_PID_Jets_SystematicsEtaDown_results_LLFit__%s_2_reg1_regFac1.00_%s_idSpectra%s.root"
              , path.Data(), modeString[mode].Data(), muonString.Data(), jetPtString.Data());
        fileNames[2] = 
          Form("%s/bhess_PID_Jets_SystematicsEtaUp_results_LLFit__%s_2_reg1_regFac1.00_%s_idSpectra%s.root"
              , path.Data(), modeString[mode].Data(), muonString.Data(), jetPtString.Data());
        
        TString histTitles[numFiles];
        histTitles[0] = "Standard #eta map";
        histTitles[1] = "#eta correction - 1%, p > 0.45 GeV/c; - %4, p < 0.45 GeV/c";
        histTitles[2] = "#eta correction + 1%, p > 0.45 GeV/c; + %4, p < 0.45 GeV/c";

        outFileTitle = Form("%s_%s", modeString[mode].Data(), outFileTitle.Data());
        SystematicErrorEstimation(path, outFileTitle, fileNames, histTitles, numFiles, nSigma, ignoreSigmaErrors, setMean, 
                                  setMeanLowerThreshold, setMeanUpperThreshold);  
      }
      {
        Bool_t setMean = kFALSE;
        const Int_t numFiles = 3;
          
        TString outFileTitle = Form("Splines%s", jetPtString.Data());
        TString fileNames[numFiles];
        fileNames[0] = 
          Form("%s/bhess_PID_Jets_results_LLFit__%s_2_reg1_regFac1.00_%s_idSpectra%s.root"
              , path.Data(), modeString[mode].Data(), muonString.Data(), jetPtString.Data());
        fileNames[1] =
          Form("%s/bhess_PID_Jets_SystematicsSplinesDown_results_LLFit__%s_2_reg1_regFac1.00_%s_idSpectra%s.root"
              , path.Data(), modeString[mode].Data(), muonString.Data(), jetPtString.Data());
        fileNames[2] =
          Form("%s/bhess_PID_Jets_SystematicsSplinesUp_results_LLFit__%s_2_reg1_regFac1.00_%s_idSpectra%s.root"
              , path.Data(), modeString[mode].Data(), muonString.Data(), jetPtString.Data());
        
        TString histTitles[numFiles];
        histTitles[0] = "Standard splines";
        histTitles[1] = "Splines - 0.2%";
        histTitles[2] = "Splines + 0.2%";
        
        outFileTitle = Form("%s_%s", modeString[mode].Data(), outFileTitle.Data());
        SystematicErrorEstimation(path, outFileTitle, fileNames, histTitles, numFiles, nSigma, ignoreSigmaErrors, setMean, 
                                  setMeanLowerThreshold, setMeanUpperThreshold);
      }
      {
        Bool_t setMean = kFALSE;
        const Int_t numFiles = 3;
          
        TString outFileTitle = Form("TOFmode%s", jetPtString.Data());
        TString fileNames[numFiles];
        fileNames[0] = 
          Form("%s/bhess_PID_Jets_results_LLFit__%s_2_reg1_regFac1.00_%s_idSpectra%s_TOFpatched.root"
              , path.Data(), modeString[mode].Data(), muonString.Data(), jetPtString.Data());
        fileNames[1] =
          Form("%s/bhess_PID_Jets_SystematicsTOF0_results_LLFit__%s_2_reg1_regFac1.00_%s_idSpectra%s_TOFpatched.root"
              , path.Data(), modeString[mode].Data(), muonString.Data(), jetPtString.Data());
        fileNames[2] =
          Form("%s/bhess_PID_Jets_SystematicsTOF2_results_LLFit__%s_2_reg1_regFac1.00_%s_idSpectra%s_TOFpatched.root"
              , path.Data(), modeString[mode].Data(), muonString.Data(), jetPtString.Data());
        
        TString histTitles[numFiles];
        histTitles[0] = "Standard TOF mode 1";
        histTitles[1] = "TOF mode 0";
        histTitles[2] = "TOF mode 2";
        
        outFileTitle = Form("%s_%s", modeString[mode].Data(), outFileTitle.Data());
        SystematicErrorEstimation(path, outFileTitle, fileNames, histTitles, numFiles, nSigma, ignoreSigmaErrors, setMean, 
                                  setMeanLowerThreshold, setMeanUpperThreshold);
      }
      
      // TODO: NO muons yet for direct fit vs. z,xi
      if (!useJetPt || mode == 0) {
        Bool_t setMean = kFALSE;
        const Int_t numFiles = 2;//3;
        
        TString outFileTitle = Form("Muons%s", jetPtString.Data());
        
        TString fileNames[numFiles];
        fileNames[0] =
          Form("%s/bhess_PID_Jets_results_LLFit__%s_2_reg1_regFac1.00_muonToElTunedOnMCHybridTrackCutsJets_idSpectra%s.root"
              , path.Data(), modeString[mode].Data(), jetPtString.Data());
        fileNames[1] =
          Form("%s/bhess_PID_Jets_results_LLFit__%s_2_reg1_regFac1.00_muonsEqualElectrons_idSpectra%s.root"
              , path.Data(), modeString[mode].Data(), jetPtString.Data());
        //fileNames[2] =
          //Form("%s/bhess_PID_Jets_results_LLFit__%s_2_reg1_regFac1.00_noMuons_idSpectra%s.root"
          //     , path.Data(), modeString[mode].Data(), jetPtString.Data());
        
        TString histTitles[numFiles];
        histTitles[0] = "Muon_Electron_Ratio_Tuned_To_MC";
        histTitles[1] = "Muon_Electron_Ratio_Is_Unity";
        //histTitles[2] = "NoMuons";

        outFileTitle = Form("%s_%s", modeString[mode].Data(), outFileTitle.Data());
        SystematicErrorEstimation(path, outFileTitle, fileNames, histTitles, numFiles, nSigma, ignoreSigmaErrors, setMean, 
                                  setMeanLowerThreshold, setMeanUpperThreshold);  
      }
      
      if (!useJetPt)
        break;
    }
  }
  */
}