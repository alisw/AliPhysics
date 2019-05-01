AliHistogramManager* gHistManData = 0x0;
AliHistogramManager* gHistManMC = 0x0;
AliHistogramManager* gHistManMCPerRun = 0x0;
AliResonanceFits gFits;

TH1* gOldPointer = 0x0;

const Double_t gkV0ANDcrossSection2015 = 51.2;              // mb  (https://twiki.cern.ch/twiki/bin/viewauth/ALICE/EventNormalization#proton_proton_at_sqrt_s_5_TeV)
const Double_t gkV0ANDcrossSection2015Err = 2.3;            // in percents

const Double_t gkV0ANDcrossSection2017 = 50.77;              // mb (the same value as 2015)
const Double_t gkV0ANDcrossSection2017Err = 2.1;            // in percents (conservative as the vdM scan for 2017 is not ready)

const Double_t gkDeltaY = 1.8;                              //  measurement done in -0.9 -> +0.9
const Double_t gkBranchingRatio = 0.05971;                  //  branching ratio of J/psi -> ee from PDG (2018 version)
const Double_t gkBranchingRatioErr = 0.00032;               //  uncertainty on branching ratio from PDG (2018 version)

const Char_t* gkDirTreesLHC15n = "/home/iarsene/data/2015/LHC15n/dstTrees_1192_20180627/";
const Char_t* gkDirTreesLHC17pq = "/home/iarsene/data/2017/LHC17pq/dstTrees_20180608_LEGO1162_1163_1164_1165/";
const Char_t* gkDirTreesLHC17pq_FAST = "/home/iarsene/data/2017/LHC17pq/dstTrees_20180608_LEGO1162_1164_FAST/";
const Char_t* gkDirTreesLHC17pq_CENT_woSDD = "/home/iarsene/data/2017/LHC17pq/dstTrees_20180608_LEGO1163_1165_CENT_woSDD/";

Char_t* gCurrentDataFilename = "";

const Int_t gkNRuns = 72;
TString gRunListStr = "244340;244343;244351;244355;244359;244364;244377;244411;244416;244418;244421;244453;244456;244480;244481;244482;244483;244484;244531;244540;244542;244617;244618;244619;244626;244627;244628;282008;282016;282021;282025;282030;282031;282050;282051;282078;282098;282099;282118;282119;282120;282122;282123;282125;282126;282127;282146;282147;282189;282206;282224;282227;282229;282230;282247;282302;282303;282304;282305;282306;282307;282309;282312;282313;282314;282340;282341;282342;282343;282365;282366;282367";
const Int_t gRunList[gkNRuns] =  { 
   // LHC15n
   244340,244343,244351,244355,244359,244364,244377,244411,244416,244418,
   244421,244453,244456,244480,244481,244482,244483,244484,244531,244540,
   244542,244617,244618,244619,244626,244627,244628,
   // LHC17p
   282008,282016,282021,282025,282030,282031,282050,282051,282078,282098,
   282099,282118,282119,282120,282122,282123,282125,282126,282127,282146,
   282147,282189,282206,282224,282227,282229,282230,282247,282302,282303,
   282304,282305,282306,282307,282309,282312,282313,282314,282340,282341,
   282342,282343,
   // LHC17q
   282365,282366,282367
};

const Int_t gkNSigCountingVariations = 15;
Double_t gmSigLims[gkNSigCountingVariations][2] = {
   {2.92, 3.16}, {2.88, 3.16}, {2.84, 3.16}, {2.80, 3.16}, {2.76, 3.16}, 
   {2.92, 3.12}, {2.88, 3.12}, {2.84, 3.12}, {2.80, 3.12}, {2.76, 3.12}, 
   {2.92, 3.08}, {2.88, 3.08}, {2.84, 3.08}, {2.80, 3.08}, {2.76, 3.08}
   //{2.92, 3.04}, {2.88, 3.04}, {2.84, 3.04}, {2.80, 3.04}, {2.76, 3.04}
};

const Int_t gkNMassFitVariations = 6;
Double_t gmFitLims[gkNMassFitVariations][2] = {
   {1.20, 4.72}, {1.20, 4.52}, {1.20, 5.0}, 
   {1.24, 4.72}, {1.28, 4.72}, {1.32, 4.72}
};

/*
const Int_t gkNMassFitVariations = 1;
Double_t gmFitLims[gkNMassFitVariations][2] = {
   {1.20, 4.20}
};
*/
/*
const Int_t gkNMassExclVariations = 7;
Double_t gmExclLims[gkNMassExclVariations][2] = {
   {2.5, 3.72}, {2.58, 3.72}, {2.66, 3.72}, {2.74, 3.72}, 
   {2.5, 3.52}, {2.5, 3.32}, {2.5, 3.16}
};
*/

const Int_t gkNMassExclVariations = 1;
Double_t gmExclLims[gkNMassExclVariations][2] = {
   {2.5, 3.72}
};



//________________________________________________________________________________________
Int_t GetPeriodFromRun(Int_t run) {
   // 
   // return an integer containing the period code based on the run number
   // 
   // 0- LHC15n
   // 1- LHC17p
   // 2- LHC17q
   
   if (run >= 244340 && run <= 244628) 
      return 0;
   if (run >= 282008 && run <= 282343) 
      return 1;
   if (run >= 282365 && run <= 282367) 
      return 2;      
}

//________________________________________________________________________________________
void ConfigureSignalExtraction(Double_t ptMin = -1.0, Double_t ptMax = -1.0) {
   // 
   //  configure the AliResonanceFits object
   //
   gFits.AddVariable(AliReducedVarManager::kVtxZ, 2);
//   gFits.AddVariable(AliReducedVarManager::kSPDntracklets, 3);
   gFits.AddVariable(AliReducedVarManager::kPt, 1);
   gFits.AddVariable(AliReducedVarManager::kMass, 0);
   gFits.SetVarRange(AliReducedVarManager::kMass, 1.5, 5.0);
   if (ptMin >= 0.0 || ptMax >= 0.0)
      gFits.SetVarRange(AliReducedVarManager::kPt, ptMin+1.0e-6, ptMax-1.0e-6);
   
   gFits.SetBkgMethod(AliResonanceFits::kBkgMixedEvent);    // kBkgMixedEvent,  kBkgLikeSign, kBkgFitFunction, kBkgMixedEventAndResidualFit
   gFits.SetMEMatchingMethod(AliResonanceFits::kMatchSEOS); // kMatchSEOS, kMatchSELS
   gFits.SetLSmethod(AliResonanceFits::kLSArithmeticMean);  // kLSArithmeticMean, kLSGeometricMean
   gFits.SetUseRfactorCorrection(kTRUE);
   gFits.SetMassFitRange(1.0,4.72);
   gFits.SetMassExclusionRange(2.5,3.72);
   gFits.SetScaleSummedBkg(kTRUE);
   //gFits.SetUseSignificantZero(kFALSE);
   gFits.SetScalingOption(AliResonanceFits::kScaleEntries); // kScaleEntries, kScaleFit, kScaleWeightedAverage
   gFits.SetWeightedAveragePower(2.0);                      // 2 - is default, but it can be modified in some specific situations
   // kMinuitMethodChi2, kMinuitMethodLikelihood
   gFits.SetMinuitFitOption(AliResonanceFits::kMinuitMethodLikelihood);
   //TF1* fitFunc = new TF1("fitFunc", "pol2/pol3(3)", 1.0, 4.6);
   //TF1* fitFunc = new TF1("fitFunc", "pol2/pol3(3)", 1.0, 4.6);
   //TF1* fitFunc = new TF1("fitFunc", "[0]*pol1(1)/x", 0., 5.0);
   //fitFunc->SetParLimits(1, -1.0e+10, 0.); 
   //fitFunc->SetParLimits(2, 6., 1.0e+10);
   //fitFunc->SetParLimits(3, -1.0e+10, 0.);
   //fitFunc->SetParLimits(4, 0.0, 1.0e+10);

   //fitFunc->SetParameter(1, 10.); fitFunc->SetParLimits(1, 5., 100.);
   TF1* fitFunc = new TF1("fitFunc", "pol1", 0., 5.0);
   //fitFunc->SetParameter(0, 10.); fitFunc->SetParameter(1,-0.1);
   gFits.SetBkgFitFunction(fitFunc);
}


//____________________________________________________________________________________________________
void ConfigureSignalExtractionME(Double_t ptMin = -1.0, Double_t ptMax = -1.0, 
                     Double_t minFit = 1.5, Double_t maxFit = 4.72, 
                     Double_t minExcl = 2.5, Double_t maxExcl = 3.72) {
   // 
   //  configure the AliResonanceFits object
   //
   gFits.AddVariable(AliReducedVarManager::kVtxZ, 2);
   //gFits.AddVariable(AliReducedVarManager::kSPDntracklets, 3);
   gFits.AddVariable(AliReducedVarManager::kPt, 1);
   gFits.AddVariable(AliReducedVarManager::kMass, 0);
   gFits.SetVarRange(AliReducedVarManager::kMass, 1.0, 5.0);
   if (ptMin >= 0.0 || ptMax >= 0.0)
      gFits.SetVarRange(AliReducedVarManager::kPt, ptMin+1.0e-6, ptMax-1.0e-6);

   gFits.SetBkgMethod(AliResonanceFits::kBkgMixedEvent);
   gFits.SetMEMatchingMethod(AliResonanceFits::kMatchSEOS); // kMatchSEOS, kMatchSELS
   gFits.SetLSmethod(AliResonanceFits::kLSArithmeticMean);  // kLSArithmeticMean, kLSGeometricMean
   gFits.SetUseRfactorCorrection(kTRUE);
   gFits.SetMassFitRange(minFit,maxFit);
   gFits.SetMassExclusionRange(minExcl,maxExcl);
   gFits.SetScaleSummedBkg(kTRUE);
   //gFits.SetUseSignificantZero(kFALSE);
   gFits.SetScalingOption(AliResonanceFits::kScaleEntries); // kScaleEntries, kScaleFit, kScaleWeightedAverage
   gFits.SetWeightedAveragePower(2.0);                      // 2 - is default, but it can be modified in some specific situations
   // kMinuitMethodChi2, kMinuitMethodLikelihood
   gFits.SetMinuitFitOption(AliResonanceFits::kMinuitMethodChi2);
   
}

//____________________________________________________________________________________________________
void ConfigureSignalExtractionLS(Double_t ptMin = -1.0, Double_t ptMax = -1.0, 
                                 Double_t minFit = 1.5, Double_t maxFit = 4.72, 
                                 Double_t minExcl = 2.5, Double_t maxExcl = 3.72) {
      // 
      //  configure the AliResonanceFits object
      //
      gFits.AddVariable(AliReducedVarManager::kVtxZ, 2);
      //gFits.AddVariable(AliReducedVarManager::kSPDntracklets, 3);
      gFits.AddVariable(AliReducedVarManager::kPt, 1);
      gFits.AddVariable(AliReducedVarManager::kMass, 0);
      gFits.SetVarRange(AliReducedVarManager::kMass, 1.0, 5.0);
      if (ptMin >= 0.0 || ptMax >= 0.0)
         gFits.SetVarRange(AliReducedVarManager::kPt, ptMin+1.0e-6, ptMax-1.0e-6);

      gFits.SetBkgMethod(AliResonanceFits::kBkgLikeSign);
      gFits.SetLSmethod(AliResonanceFits::kLSArithmeticMean);  // kLSArithmeticMean, kLSGeometricMean
      gFits.SetUseRfactorCorrection(kTRUE);
      gFits.SetMassFitRange(minFit,maxFit);
      gFits.SetScaleSummedBkg(kTRUE);
      gFits.SetMassExclusionRange(minExcl,maxExcl);
      //gFits.SetUseSignificantZero(kFALSE);
      gFits.SetScalingOption(AliResonanceFits::kScaleEntries); // kScaleEntries, kScaleFit, kScaleWeightedAverage
      gFits.SetWeightedAveragePower(2.0);                      // 2 - is default, but it can be modified in some specific situations
      // kMinuitMethodChi2, kMinuitMethodLikelihood
      gFits.SetMinuitFitOption(AliResonanceFits::kMinuitMethodChi2);

}

//____________________________________________________________________________________________________
void ConfigureSignalExtractionHybrid(Double_t ptMin = -1.0, Double_t ptMax = -1.0, 
                                    Double_t minFit = 1.5, Double_t maxFit = 4.72) {
      // 
      //  configure the AliResonanceFits object
      //
      gFits.AddVariable(AliReducedVarManager::kVtxZ, 2);
      //gFits.AddVariable(AliReducedVarManager::kSPDntracklets, 3);
      gFits.AddVariable(AliReducedVarManager::kPt, 1);
      gFits.AddVariable(AliReducedVarManager::kMass, 0);
      gFits.SetVarRange(AliReducedVarManager::kMass, 1.0, 5.0);
      if (ptMin >= 0.0 && ptMax >= 0.0)
         gFits.SetVarRange(AliReducedVarManager::kPt, ptMin+1.0e-6, ptMax-1.0e-6);

      gFits.SetBkgMethod(AliResonanceFits::kBkgMixedEventAndResidualFit);
      gFits.SetMEMatchingMethod(AliResonanceFits::kMatchSELS); // kMatchSEOS, kMatchSELS
      gFits.SetLSmethod(AliResonanceFits::kLSArithmeticMean);  // kLSArithmeticMean, kLSGeometricMean
      gFits.SetUseRfactorCorrection(kTRUE);
      if (TMath::Abs(ptMin-0.0)<1.0e-5 && ptMax<2.0 && minFit<1.8)
         gFits.SetMassFitRange(1.8,maxFit);
      else
         gFits.SetMassFitRange(minFit,maxFit);
      gFits.SetMassExclusionRange(2.5,3.2);                 //   used for an initial guess of the bkg function parameters
      gFits.SetUseSignificantZero(kTRUE);
      gFits.SetScalingOption(AliResonanceFits::kScaleEntries); // kScaleEntries, kScaleFit, kScaleWeightedAverage
      gFits.SetWeightedAveragePower(2.0);                      // 2 - is default, but it can be modified in some specific situations
      // kMinuitMethodChi2, kMinuitMethodLikelihood
      gFits.SetMinuitFitOption(AliResonanceFits::kMinuitMethodChi2);
      //TF1* fitFunc = new TF1("fitFunc", "[0]*exp([1]*x)", 0., 5.0);
      //TF1* fitFunc = new TF1("fitFunc", "expo", 0., 5.0);
      //fitFunc->SetParameters(0.1, 0.1);
      TF1* fitFunc = new TF1("fitFunc", "pol2", 0., 5.0);
      gFits.SetBkgFitFunction(fitFunc);
}

//____________________________________________________________________________________________________
void ConfigureSignalExtractionFit(Double_t ptMin = -1.0, Double_t ptMax = -1.0, 
                                 Double_t minFit = 1.5, Double_t maxFit = 4.72) {
      // 
      //  configure the AliResonanceFits object
      //
      gFits.AddVariable(AliReducedVarManager::kVtxZ, 2);
      //gFits.AddVariable(AliReducedVarManager::kSPDntracklets, 3);
      gFits.AddVariable(AliReducedVarManager::kPt, 1);
      gFits.AddVariable(AliReducedVarManager::kMass, 0);
      gFits.SetVarRange(AliReducedVarManager::kMass, 1.0, 5.0);
      if (ptMin >= 0.0 || ptMax >= 0.0)
         gFits.SetVarRange(AliReducedVarManager::kPt, ptMin+1.0e-6, ptMax-1.0e-6);

      gFits.SetBkgMethod(AliResonanceFits::kBkgFitFunction);
      gFits.SetMassFitRange(minFit,maxFit);
      gFits.SetUseSignificantZero(kTRUE);
      gFits.SetScalingOption(AliResonanceFits::kScaleEntries);
      gFits.SetMinuitFitOption(AliResonanceFits::kMinuitMethodChi2);
      TF1* fitFunc = new TF1("fitFunc", "pol2/pol3(3)", 1.5, 5.0);
      gFits.SetBkgFitFunction(fitFunc);
}


//________________________________________________________________________________________
void RunSignalExtraction(const Char_t* filename, const Char_t* settingName, const Char_t* filenameMC="") {
   //
   // run signal extraction
   // 
   /*gHistManData = new AliHistogramManager();
   gHistManData->InitFile(filename, "jpsi2eeHistos");
      
   if(filenameMC[0]!='\0') {
      gHistManMC = new AliHistogramManager();
      gHistManMC->InitFile(filenameMC, "jpsi2eeHistos");
   } */
   
   AliHistogramManager* mcMan1 = new AliHistogramManager();
   mcMan1->InitFile("/home/iarsene/work/ALICE/pp2017Analysis/analysisOutputs/20180725_TrackingSystematics_MC/LHC15n/AnalysisHistograms_jpsi2ee_pp2017_MC.root");
   THnF* sig15n = (THnF*)mcMan1->GetHistogram(Form("PairSEPM_%s_TrueElectron",settingName),"PairInvMass");
   sig15n->GetAxis(1)->SetRangeUser(6.0, 7.0);
   TH1D* proj15n = sig15n->Projection(0);
   proj15n->SetLineColor(2);
   proj15n->SetLineWidth(2);
   
   AliHistogramManager* mcMan2 = new AliHistogramManager();
   mcMan2->InitFile("/home/iarsene/work/ALICE/pp2017Analysis/analysisOutputs/20180725_TrackingSystematics_MC/LHC17pq/AnalysisHistograms_jpsi2ee_pp2017_MC.root");
   THnF* sig17pq = (THnF*)mcMan2->GetHistogram(Form("PairSEPM_%s_TrueElectron",settingName),"PairInvMass");
   sig17pq->GetAxis(1)->SetRangeUser(6.0, 7.0);
   TH1D* proj17pq = sig17pq->Projection(0);
   proj17pq->SetLineColor(4);
   proj17pq->SetLineWidth(2);
   
   proj17pq->Draw();
   proj15n->Draw("same");
   /*
   ConfigureSignalExtractionHybrid();
   SetHistograms(settingName, "TrueElectron");
   gFits.Process();
   gFits.ComputeOutputValues(2.92, 3.16);   
   Draw(&gFits, 2.92, 3.16);
   */
}



//_________________________________________________________________________________________
void BuildSpectra(const Char_t* outfilename="SpectraPP5TeV.root") {
   // 
   // get the data from Results.root files and build spectra
   // 
   //const Char_t* trackingSystDir = "/home/iarsene/work/ALICE/pp2017Analysis/analysisOutputs/20180725_TrackingSystematics/LHC17pq/sigExtrHybrid_final";
   const Char_t* trackingSystDir = "/home/iarsene/work/ALICE/pp2017Analysis/analysisOutputs/20180725_TrackingSystematics/LHC17pq/sigExtrHybrid_finalCS";
   const Char_t* pidSystDir = "/home/iarsene/work/ALICE/pp2017Analysis/analysisOutputs/20180725_PIDSystematicsSPDany/LHC17pq/sigExtrHybrid_final";
   const Char_t* pidSingleEleSystDir = "/home/iarsene/work/ALICE/pp2017Analysis/analysisOutputs/20180922_PIDsystematicsSingleEle/LHC17pq";
   
   
   /*
   const Int_t kNPtBinsSpectra = 10;
   Double_t ptBinLims[kNPtBinsSpectra][2] = {
      {0.0, 1.0}, {1.0, 2.0}, {2.0, 3.0}, {3.0, 4.0}, {4.0, 5.0}, {5.0, 7.0}, {7.0, 10.0}, {1.0, 4.0}, {5.0, 10.0}, {-1.0, -1.0}
   };*/
   const Int_t kNPtBinsSpectra = 14;
   Double_t ptBinLims[kNPtBinsSpectra][2] = {
      {0.0, 0.15}, {0.0, 1.0}, {0.15, 1.3}, {1.0, 2.0}, {1.3, 3.0}, 
      {2.0, 3.0},  {3.0, 4.0}, {3.0, 5.0},  {4.0, 5.0}, {5.0, 7.0}, 
      {7.0, 10.0}, {1.0, 4.0}, {5.0, 10.0}, {-1.0, -1.0}
   };
   const Int_t kLowPtWideBin = 11;                           // 1 - 4 GeV/c,  bin used for the signal shape and tracking cut variation uncertainty in low pt bins (pt<5 GeV)
   const Int_t kHighPtWideBin = 12;                          // 5 - 10 GeV/c,  bin used for the signal shape and tracking cut variation uncertainty in high pt bins (pt>5 GeV/c)
   //Bool_t useDataPoint[kNPtBinsSpectra] = {kTRUE, kTRUE, kTRUE, kTRUE, kTRUE, kTRUE, kTRUE, kFALSE, kFALSE, kTRUE};
   // Used data-points for the pp paper spectrum
   /*Bool_t useDataPoint[kNPtBinsSpectra] = {
      kFALSE, kTRUE, kFALSE, kTRUE, kFALSE, 
      kTRUE, kTRUE, kFALSE, kTRUE, kTRUE, 
      kTRUE, kFALSE, kFALSE, kTRUE
   };*/
   // Used data-points for the Pb-Pb paper RAA
   Bool_t useDataPoint[kNPtBinsSpectra] = {
      kTRUE, kFALSE, kTRUE, kFALSE, kTRUE, 
      kFALSE, kFALSE, kTRUE, kFALSE, kFALSE, 
      kFALSE, kFALSE, kTRUE, kTRUE
   };
   
   Double_t csTracking[kNPtBinsSpectra][5] = {{0.0}};
   Double_t cspol2Tracking[kNPtBinsSpectra][5] = {{0.0}};
   Double_t csPID[kNPtBinsSpectra][5] = {{0.0}};
   Double_t cspol2PID[kNPtBinsSpectra][5] = {{0.0}};
   Double_t csMCkine[kNPtBinsSpectra][2] = {{0.0}};

   for(Int_t ipt = 0;ipt<kNPtBinsSpectra;++ipt) {
      cout << "Syst pt " << ptBinLims[ipt][0] << " -> " << ptBinLims[ipt][1] << endl;
      TString dirName = Form("pt%.2f_%.2f", ptBinLims[ipt][0], ptBinLims[ipt][1]);
      if(ptBinLims[ipt][0] < 0.0 && ptBinLims[ipt][1] < 0.0) dirName = "ptIntegrated";
      if(ptBinLims[ipt][1]>0.0 && ptBinLims[ipt][1] < 1.301) 
         dirName = Form("pt%.2f_%.2f_pol2", ptBinLims[ipt][0], ptBinLims[ipt][1]);
      
      cout << "dirName :: " << dirName.Data() << endl;
      
      // get cross-section and syst from tracking variations
      if(ptBinLims[ipt][1]>0.0 && ptBinLims[ipt][1] < 1.301) 
         GetCrossSection(Form("%s/%s/Results.root", trackingSystDir, dirName.Data()), Form("%s/%s/ResultsRefiltered.root", trackingSystDir, dirName.Data()), cspol2Tracking[ipt]);
      else
         GetCrossSection(Form("%s/%s/Results.root", trackingSystDir, dirName.Data()), Form("%s/%s/ResultsRefiltered.root", trackingSystDir, dirName.Data()), csTracking[ipt]);
         
      // get cross-section and syst from pid variations
      if(ptBinLims[ipt][1]>0.0 && ptBinLims[ipt][1] < 1.301) 
         GetCrossSection(Form("%s/%s/Results.root", pidSystDir, dirName.Data()), Form("%s/%s/ResultsRefiltered.root", pidSystDir, dirName.Data()), cspol2PID[ipt]);
      else
         GetCrossSection(Form("%s/%s/Results.root", pidSystDir, dirName.Data()), Form("%s/%s/ResultsRefiltered.root", pidSystDir, dirName.Data()), csPID[ipt]);
   }
   
   // get systematics from the ITS-TPC single track matching efficiency
   TFile* fileWithoutITSTPCmatch = TFile::Open(Form("%s/JpsiKineAccWithoutMatchingSyst.root", trackingSystDir));
   TH1F* effDiffWO = (TH1F*)fileWithoutITSTPCmatch->Get("hPtEff");
   TProfile* effIntegratedWO = (TProfile*)fileWithoutITSTPCmatch->Get("hEffIntegrated");
   TFile* fileWithITSTPCmatch = TFile::Open(Form("%s/JpsiKineAccWithMatchingSyst.root", trackingSystDir));
   TH1F* effDiffW = (TH1F*)fileWithITSTPCmatch->Get("hPtEff");
   TProfile* effIntegratedW = (TProfile*)fileWithITSTPCmatch->Get("hEffIntegrated");
   
   // get the cross-section shift due to the non-prompt fraction tuning to data
   TFile* nonpromptCorrFile = TFile::Open(Form("%s/NonPromptFractionShift.root", trackingSystDir));
   TH1D* nonpromptShift = (TH1D*)nonpromptCorrFile->Get("correction");
   
   // get efficiencies from single electron PID systematics
   TString standardPidSystName = "Standard";
   TH1F* standardPidSystHist[2] = {0x0};
   TH1F* standardPidSystHistIntegrated[2] = {0x0};
   TFile* pidFile = TFile::Open(Form("%s/data_Standard.root", pidSingleEleSystDir));
   standardPidSystHist[0] = (TH1F*)pidFile->Get("hPtEff")->Clone("data_Standard"); 
   standardPidSystHist[0]->SetDirectory(0x0);
   standardPidSystHistIntegrated[0] = (TH1F*)pidFile->Get("hEffIntegrated")->Clone("data_Standard_ptIntegrated"); 
   standardPidSystHistIntegrated[0]->SetDirectory(0x0);
   pidFile->Close();
   pidFile = TFile::Open(Form("%s/mc_Standard.root", pidSingleEleSystDir));
   standardPidSystHist[1] = (TH1F*)pidFile->Get("hPtEff")->Clone("mc_Standard"); 
   standardPidSystHist[1]->SetDirectory(0x0);
   standardPidSystHistIntegrated[1] = (TH1F*)pidFile->Get("hEffIntegrated")->Clone("mc_Standard_ptIntegrated"); 
   standardPidSystHistIntegrated[1]->SetDirectory(0x0);
   pidFile->Close();
   
   TString pionPidSystName[2] = {"electron0_prot0_pion1", "electron0_prot0_pion2"};
   TH1F* pionPidSystHist[2][2] = {{0x0}};
   TH1F* pionPidSystHistIntegrated[2][2] = {{0x0}};
   for(Int_t i = 0; i<2; i++) {
      TFile* pidFile = TFile::Open(Form("%s/data_%s.root", pidSingleEleSystDir, pionPidSystName[i].Data()));
      pionPidSystHist[i][0] = (TH1F*)pidFile->Get("hPtEff")->Clone(Form("data_%s", pionPidSystName[i].Data())); 
      pionPidSystHist[i][0]->SetDirectory(0x0);
      pionPidSystHistIntegrated[i][0] = (TH1F*)pidFile->Get("hEffIntegrated")->Clone(Form("data_%s_ptIntegrated", pionPidSystName[i].Data())); 
      pionPidSystHistIntegrated[i][0]->SetDirectory(0x0);
      pidFile->Close();
      pidFile = TFile::Open(Form("%s/mc_%s.root", pidSingleEleSystDir, pionPidSystName[i].Data()));
      pionPidSystHist[i][1] = (TH1F*)pidFile->Get("hPtEff")->Clone(Form("mc_%s", pionPidSystName[i].Data())); 
      pionPidSystHist[i][1]->SetDirectory(0x0);
      pionPidSystHistIntegrated[i][1] = (TH1F*)pidFile->Get("hEffIntegrated")->Clone(Form("mc_%s_ptIntegrated", pionPidSystName[i].Data())); 
      pionPidSystHistIntegrated[i][1]->SetDirectory(0x0);
      pidFile->Close();
   }
   
   TString protPidSystName[2] = {"electron0_prot1_pion0", "electron0_prot2_pion0"};
   TH1F* protPidSystHist[2][2] = {{0x0}};
   TH1F* protPidSystHistIntegrated[2][2] = {{0x0}};
   for(Int_t i = 0; i<2; i++) {
      TFile* pidFile = TFile::Open(Form("%s/data_%s.root", pidSingleEleSystDir, protPidSystName[i].Data()));
      protPidSystHist[i][0] = (TH1F*)pidFile->Get("hPtEff")->Clone(Form("data_%s", protPidSystName[i].Data())); 
      protPidSystHist[i][0]->SetDirectory(0x0);
      protPidSystHistIntegrated[i][0] = (TH1F*)pidFile->Get("hEffIntegrated")->Clone(Form("data_%s", protPidSystName[i].Data())); 
      protPidSystHistIntegrated[i][0]->SetDirectory(0x0);
      pidFile->Close();
      pidFile = TFile::Open(Form("%s/mc_%s.root", pidSingleEleSystDir, protPidSystName[i].Data()));
      protPidSystHist[i][1] = (TH1F*)pidFile->Get("hPtEff")->Clone(Form("mc_%s", protPidSystName[i].Data())); 
      protPidSystHist[i][1]->SetDirectory(0x0);
      protPidSystHistIntegrated[i][1] = (TH1F*)pidFile->Get("hEffIntegrated")->Clone(Form("mc_%s", protPidSystName[i].Data())); 
      protPidSystHistIntegrated[i][1]->SetDirectory(0x0);
      pidFile->Close();
   }
   
   TString elePidSystName[4] = {"electron1_prot0_pion0", "electron2_prot0_pion0", "electron3_prot0_pion0", "electron4_prot0_pion0"};
   TH1F* elePidSystHist[4][2] = {{0x0}};
   TH1F* elePidSystHistIntegrated[4][2] = {{0x0}};
   for(Int_t i = 0; i<4; i++) {
      TFile* pidFile = TFile::Open(Form("%s/data_%s.root", pidSingleEleSystDir, elePidSystName[i].Data()));
      elePidSystHist[i][0] = (TH1F*)pidFile->Get("hPtEff")->Clone(Form("data_%s", elePidSystName[i].Data())); 
      elePidSystHist[i][0]->SetDirectory(0x0);
      elePidSystHistIntegrated[i][0] = (TH1F*)pidFile->Get("hEffIntegrated")->Clone(Form("data_%s", elePidSystName[i].Data())); 
      elePidSystHistIntegrated[i][0]->SetDirectory(0x0);
      pidFile->Close();
      pidFile = TFile::Open(Form("%s/mc_%s.root", pidSingleEleSystDir, elePidSystName[i].Data()));
      elePidSystHist[i][1] = (TH1F*)pidFile->Get("hPtEff")->Clone(Form("mc_%s", elePidSystName[i].Data())); 
      elePidSystHist[i][1]->SetDirectory(0x0);
      elePidSystHistIntegrated[i][1] = (TH1F*)pidFile->Get("hEffIntegrated")->Clone(Form("mc_%s", elePidSystName[i].Data())); 
      elePidSystHistIntegrated[i][1]->SetDirectory(0x0);
      pidFile->Close();
   }
      
   enum Values {
      kCS = 0, kCSstat, kCSsyst, kCSfit, kCSsig
   };
   TGraphErrors* spectraStat = new TGraphErrors();                        spectraStat->SetName("StatUncert");
   TGraphErrors* spectraSystTracking = new TGraphErrors();                spectraSystTracking->SetName("SystUncertTracking");
   TGraphErrors* spectraSystTrackingCutVariations = new TGraphErrors();   spectraSystTrackingCutVariations->SetName("SystUncertTrackingCutVariations");
   TGraphErrors* spectraSystTrackingITSTPCmatch = new TGraphErrors();     spectraSystTrackingITSTPCmatch->SetName("SystUncertTrackingITSTPCmatch");
   TGraphErrors* spectraSystPID = new TGraphErrors();                     spectraSystPID->SetName("SystUncertPID");
   TGraphErrors* spectraSystPIDstandard = new TGraphErrors();             spectraSystPIDstandard->SetName("SystUncertPIDstandard");
   TGraphErrors* spectraSystPIDele = new TGraphErrors();                  spectraSystPIDele->SetName("SystUncertPIDele");
   TGraphErrors* spectraSystPIDprot = new TGraphErrors();         spectraSystPIDprot->SetName("SystUncertPIDprot");
   TGraphErrors* spectraSystPIDpion = new TGraphErrors();         spectraSystPIDpion->SetName("SystUncertPIDpion");
   TGraphErrors* spectraSystPIDsingleEle = new TGraphErrors();    spectraSystPIDsingleEle->SetName("SystUncertPIDsingleEle");
   TGraphErrors* spectraSystFit = new TGraphErrors();             spectraSystFit->SetName("SystUncertFit");
   TGraphErrors* spectraSystSig = new TGraphErrors();             spectraSystSig->SetName("SystUncertSignalShape");
   TGraphErrors* spectraSystMCkine = new TGraphErrors();          spectraSystMCkine->SetName("SystUncertMCKine");
   TGraphErrors* spectraSystUncorrelated = new TGraphErrors();    spectraSystUncorrelated->SetName("SystUncertUncorrelated");
   TGraphErrors* spectraSystCorrelated = new TGraphErrors();      spectraSystCorrelated->SetName("SystUncertCorrelated");
   TGraphErrors* spectraSystTotal = new TGraphErrors();           spectraSystTotal->SetName("SystUncertTotal");
   
   TGraphErrors* spectraStatIntegrated = new TGraphErrors();
   spectraStatIntegrated->SetName("StatUncertIntegrated");
   TGraphErrors* spectraSystTrackingIntegrated = new TGraphErrors();        
   spectraSystTrackingIntegrated->SetName("SystUncertTrackingIntegrated");
   TGraphErrors* spectraSystTrackingCutVariationsIntegrated = new TGraphErrors();        
   spectraSystTrackingCutVariationsIntegrated->SetName("SystUncertTrackingCutVariationsIntegrated");
   TGraphErrors* spectraSystTrackingITSTPCmatchIntegrated = new TGraphErrors();        
   spectraSystTrackingITSTPCmatchIntegrated->SetName("SystUncertTrackingITSTPCmatchIntegrated");
   TGraphErrors* spectraSystPIDIntegrated = new TGraphErrors();             
   spectraSystPIDIntegrated->SetName("SystUncertPIDIntegrated");
   TGraphErrors* spectraSystPIDstandardIntegrated = new TGraphErrors();          
   spectraSystPIDstandardIntegrated->SetName("SystUncertPIDstandardIntegrated");
   TGraphErrors* spectraSystPIDeleIntegrated = new TGraphErrors();          
   spectraSystPIDeleIntegrated->SetName("SystUncertPIDeleIntegrated");
   TGraphErrors* spectraSystPIDprotIntegrated = new TGraphErrors();         
   spectraSystPIDprotIntegrated->SetName("SystUncertPIDprotIntegrated");
   TGraphErrors* spectraSystPIDpionIntegrated = new TGraphErrors();         
   spectraSystPIDpionIntegrated->SetName("SystUncertPIDpionIntegrated");
   TGraphErrors* spectraSystPIDsingleEleIntegrated = new TGraphErrors();    
   spectraSystPIDsingleEleIntegrated->SetName("SystUncertPIDsingleEleIntegrated");
   TGraphErrors* spectraSystFitIntegrated = new TGraphErrors();             
   spectraSystFitIntegrated->SetName("SystUncertFitIntegrated");
   TGraphErrors* spectraSystSigIntegrated = new TGraphErrors();             
   spectraSystSigIntegrated->SetName("SystUncertSignalShapeIntegrated");
   TGraphErrors* spectraSystMCkineIntegrated = new TGraphErrors();             
   spectraSystMCkineIntegrated->SetName("SystUncertMCkineIntegrated");
   TGraphErrors* spectraSystUncorrelatedIntegrated = new TGraphErrors();    
   spectraSystUncorrelatedIntegrated->SetName("SystUncertUncorrelatedIntegrated");
   TGraphErrors* spectraSystCorrelatedIntegrated = new TGraphErrors();      
   spectraSystCorrelatedIntegrated->SetName("SystUncertCorrelatedIntegrated");
   TGraphErrors* spectraSystTotalIntegrated = new TGraphErrors();           
   spectraSystTotalIntegrated->SetName("SystUncertTotalIntegrated");
   
   // setup the data points
   Int_t currentDataPoint = 0;
   for(Int_t ipt = 0;ipt<kNPtBinsSpectra;++ipt) {
      if(!useDataPoint[ipt]) continue;                      // skip this data point 
      Double_t pt = 0.5*(ptBinLims[ipt][0]+ptBinLims[ipt][1]);
      Double_t ptErr = 0.5*(ptBinLims[ipt][1]-ptBinLims[ipt][0]);
      if(ptBinLims[ipt][0] < 0.0 && ptBinLims[ipt][1] < 0.0) {
         pt = 0.0;
         ptErr = 1.0;
      }
   
      cout << "***************************************** pt :: " << pt << endl;
         
      // MC kinematics ------------------------------------------------------------------------------
      //  Obtained using the function ComputeInputKineSyst(); each pt interval has its own output
      TString filenameKine = Form("%s/MCinputKineStats_pt%.2f_%.2f.root", trackingSystDir, ptBinLims[ipt][0], ptBinLims[ipt][1]);
      if(ptBinLims[ipt][0] < 0.0 && ptBinLims[ipt][1] < 0.0) 
         filenameKine = Form("%s/MCinputKineStats_pt0.00_10.00.root", trackingSystDir);
      
      TFile* fileKine = TFile::Open(filenameKine.Data());
      TH1F* hEffSpectra = (TH1F*)fileKine->Get("hEffSpectra");
      Double_t invEffSpectra = 1.0/hEffSpectra->GetBinContent(1);
      TH1F* hEffStandard = (TH1F*)fileKine->Get("hEffStandard");
      Double_t invEffStandard = 1.0/hEffStandard->GetBinContent(1);
      // compute the fractional cross-section shift
      csMCkine[ipt][0] = 2.0*(invEffSpectra-invEffStandard)/(invEffSpectra+invEffStandard); 
      
      // get the fractional shift needed to tune the non-prompt fraction to the one observed in data
      Double_t nonPromptCSshift = nonpromptShift->GetBinContent(nonpromptShift->GetXaxis()->FindBin(pt));
      if(ptBinLims[ipt][0] < 0.0 && ptBinLims[ipt][1] < 0.0) nonPromptCSshift = 0.0898;
      
      // compute the fractional syst uncertainty
      TH1D* hInvEffVariations = (TH1D*)fileKine->Get("hInvEff");
      csMCkine[ipt][1] = hInvEffVariations->GetRMS() / hInvEffVariations->GetMean();
      
      // switch off MC kine systematics
      cout << "mc kine shift / syst :: " << 100.0*csMCkine[ipt][0] << " / " << 100.0*csMCkine[ipt][1] << " %" << endl;
      cout << "non-prompt shift :: " << nonPromptCSshift << " %" << endl;
      csMCkine[ipt][0] += 0.01*nonPromptCSshift;
      
      //csMCkine[ipt][0] = 0.0; csMCkine[ipt][1] = 0.0;
         
      // Statistical uncertainty --------------------------------------------------------------------
      if(ptBinLims[ipt][0] < 0.0 && ptBinLims[ipt][1] < 0.0) {
         spectraStatIntegrated->SetPoint(0, pt, csTracking[ipt][kCS]*(1.0+csMCkine[ipt][0]));
         spectraStatIntegrated->SetPointError(0, ptErr, csTracking[ipt][kCSstat]);
         spectraSystMCkineIntegrated->SetPoint(0, pt, csTracking[ipt][kCS]*(1.0+csMCkine[ipt][0]));
         spectraSystMCkineIntegrated->SetPointError(0, ptErr, csTracking[ipt][kCS]*(1.0+csMCkine[ipt][0])*csMCkine[ipt][1]);
      }
      else {
         if(ptBinLims[ipt][1] < 1.301) {
            spectraStat->SetPoint(currentDataPoint, pt, cspol2Tracking[ipt][kCS]*(1.0+csMCkine[ipt][0]));
            spectraStat->SetPointError(currentDataPoint, ptErr, cspol2Tracking[ipt][kCSstat]);
            spectraSystMCkine->SetPoint(currentDataPoint, pt, cspol2Tracking[ipt][kCS]*(1.0+csMCkine[ipt][0]));
            spectraSystMCkine->SetPointError(currentDataPoint, ptErr, cspol2Tracking[ipt][kCS]*(1.0+csMCkine[ipt][0])*csMCkine[ipt][1]);
         }
         else {
            spectraStat->SetPoint(currentDataPoint, pt, csTracking[ipt][kCS]*(1.0+csMCkine[ipt][0]));
            spectraStat->SetPointError(currentDataPoint, ptErr, csTracking[ipt][kCSstat]);
            spectraSystMCkine->SetPoint(currentDataPoint, pt, csTracking[ipt][kCS]*(1.0+csMCkine[ipt][0]));
            spectraSystMCkine->SetPointError(currentDataPoint, ptErr, csTracking[ipt][kCS]*(1.0+csMCkine[ipt][0])*csMCkine[ipt][1]);
         }
      }
      
      // Tracking syst uncertainty ------------------------------------------------------------------
      //   Obtained as the quadrature of the uncertainty on the ITS-TPC matching syst and the uncertainty resulting from 
      //         additional cut variations
      if(ptBinLims[ipt][0] < 0.0 && ptBinLims[ipt][1] < 0.0) {
         spectraSystTrackingCutVariationsIntegrated->SetPoint(0, pt, csTracking[ipt][kCS]*(1.0+csMCkine[ipt][0]));
         spectraSystTrackingCutVariationsIntegrated->SetPointError(0, ptErr, csTracking[ipt][kCSsyst]);
         spectraSystTrackingITSTPCmatchIntegrated->SetPoint(0, pt, csTracking[ipt][kCS]*(1.0+csMCkine[ipt][0]));
         Double_t systITSTPCmatch = TMath::Abs( 2.0*(1.0/effIntegratedWO->GetBinContent(1)-1.0/effIntegratedW->GetBinContent(1)) / (1.0/effIntegratedWO->GetBinContent(1)+1.0/effIntegratedW->GetBinContent(1)) );
         spectraSystTrackingITSTPCmatchIntegrated->SetPointError(0, ptErr, systITSTPCmatch * csTracking[ipt][kCS]*(1.0+csMCkine[ipt][0]));
         spectraSystTrackingIntegrated->SetPoint(0, pt, csTracking[ipt][kCS]*(1.0+csMCkine[ipt][0]));
         spectraSystTrackingIntegrated->SetPointError(0, ptErr, TMath::Sqrt(csTracking[ipt][kCSsyst]*csTracking[ipt][kCSsyst] + systITSTPCmatch*systITSTPCmatch*csTracking[ipt][kCS]   *csTracking[ipt][kCS]*(1.0+csMCkine[ipt][0])*(1.0+csMCkine[ipt][0])));
      }
      // Use the tracking variations from the 2 large pt bins
      else {
         Double_t systITSTPCmatch = TMath::Abs( 2.0*(1.0/effDiffWO->GetBinContent(effDiffWO->FindBin(pt))-1.0/effDiffW->GetBinContent(effDiffW->FindBin(pt))) / (1.0/effDiffW->GetBinContent(effDiffW->FindBin(pt))+1.0/effDiffWO->GetBinContent(effDiffWO->FindBin(pt))));
         
         if(ptBinLims[ipt][1] < 1.301) {
            spectraSystTracking->SetPoint(currentDataPoint, pt, cspol2Tracking[ipt][kCS]*(1.0+csMCkine[ipt][0]));
            spectraSystTrackingCutVariations->SetPoint(currentDataPoint, pt, cspol2Tracking[ipt][kCS]*(1.0+csMCkine[ipt][0]));
            spectraSystTrackingITSTPCmatch->SetPoint(currentDataPoint, pt, cspol2Tracking[ipt][kCS]*(1.0+csMCkine[ipt][0]));
         }
         else {
            spectraSystTracking->SetPoint(currentDataPoint, pt, csTracking[ipt][kCS]*(1.0+csMCkine[ipt][0]));
            spectraSystTrackingCutVariations->SetPoint(currentDataPoint, pt, csTracking[ipt][kCS]*(1.0+csMCkine[ipt][0]));
            spectraSystTrackingITSTPCmatch->SetPoint(currentDataPoint, pt, csTracking[ipt][kCS]*(1.0+csMCkine[ipt][0]));
         }
         
         // if the pt bin is lower than 5 GeV/c, use the syst uncertainty determined for 1-4 GeV/c
         if(ptBinLims[ipt][1] < 5.001) {
            spectraSystTrackingCutVariations->SetPointError(currentDataPoint, ptErr, (csTracking[kLowPtWideBin][kCSsyst]/csTracking[kLowPtWideBin][kCS]) * csTracking[ipt][kCS]*(1.0+csMCkine[ipt][0]));
            if(ptBinLims[ipt][1] < 1.301)
               spectraSystTrackingCutVariations->SetPointError(currentDataPoint, ptErr, (csTracking[kLowPtWideBin][kCSsyst]/csTracking[kLowPtWideBin][kCS]) * cspol2Tracking[ipt][kCS]*(1.0+csMCkine[ipt][0]));
         }
         // otherwise use the one for 5-10 GeV/c
         else
            spectraSystTrackingCutVariations->SetPointError(currentDataPoint, ptErr, (csTracking[kHighPtWideBin][kCSsyst]/csTracking[kHighPtWideBin][kCS]) * csTracking[ipt][kCS]*(1.0+csMCkine[ipt][0]) );
            
         if(ptBinLims[ipt][1] < 1.301)
            spectraSystTrackingITSTPCmatch->SetPointError(currentDataPoint, ptErr, systITSTPCmatch*cspol2Tracking[ipt][kCS]*(1.0+csMCkine[ipt][0]));
         else
            spectraSystTrackingITSTPCmatch->SetPointError(currentDataPoint, ptErr, systITSTPCmatch*csTracking[ipt][kCS]*(1.0+csMCkine[ipt][0]));
         
         spectraSystTracking->SetPointError(currentDataPoint, ptErr, TMath::Sqrt(spectraSystTrackingITSTPCmatch->GetErrorY(currentDataPoint)*spectraSystTrackingITSTPCmatch->GetErrorY(currentDataPoint) +
                                                                                 spectraSystTrackingCutVariations->GetErrorY(currentDataPoint)*spectraSystTrackingCutVariations->GetErrorY(currentDataPoint)));
      }
      
      // Fit syst uncertainty -----------------------------------------------------------------------
      //  Systematics obtained after variations of the fit ranges
      if(ptBinLims[ipt][0] < 0.0 && ptBinLims[ipt][1] < 0.0) {
         spectraSystFitIntegrated->SetPoint(0, pt, csTracking[ipt][kCS]*(1.0+csMCkine[ipt][0]));
         spectraSystFitIntegrated->SetPointError(0, ptErr, csTracking[ipt][kCSfit]);
      }
      else {
         if(ptBinLims[ipt][1] < 1.301) {
            spectraSystFit->SetPoint(currentDataPoint, pt, cspol2Tracking[ipt][kCS]*(1.0+csMCkine[ipt][0]));
            spectraSystFit->SetPointError(currentDataPoint, ptErr, cspol2Tracking[ipt][kCSfit]);
         }
         else {
            spectraSystFit->SetPoint(currentDataPoint, pt, csTracking[ipt][kCS]*(1.0+csMCkine[ipt][0]));
            spectraSystFit->SetPointError(currentDataPoint, ptErr, csTracking[ipt][kCSfit]);
         }
      }
      
      // Signal shape syst uncertainty ---------------------------------------------------------------
      // Assume that the signal shape uncertainty is bin-to-bin correlated,  and to avoid including statistical fluctuations
      //   we use the syst determined from variations in the pt range 1-4 for the pt bins lower than 5 GeV/c.
      //  For higher pt bins we use the uncertainty determined for 5-10 GeV/c
      if(ptBinLims[ipt][0] < 0.0 && ptBinLims[ipt][1] < 0.0) {
         spectraSystSigIntegrated->SetPoint(0, pt, csTracking[ipt][kCS]*(1.0+csMCkine[ipt][0]));
         spectraSystSigIntegrated->SetPointError(0, ptErr, csTracking[ipt][kCSsig]);
      }
      else {
         if(ptBinLims[ipt][1] < 1.301)
            spectraSystSig->SetPoint(currentDataPoint, pt, cspol2Tracking[ipt][kCS]*(1.0+csMCkine[ipt][0]));
         else
            spectraSystSig->SetPoint(currentDataPoint, pt, csTracking[ipt][kCS]*(1.0+csMCkine[ipt][0]));
         // if the pt bin is lower than 5 GeV/c,  use the syst uncertainty determined for 1-4 GeV/c
         if(ptBinLims[ipt][1] < 5.01) {
            if(ptBinLims[ipt][1] < 1.301) 
               spectraSystSig->SetPointError(currentDataPoint, ptErr, (csTracking[kLowPtWideBin][kCSsig]/csTracking[kLowPtWideBin][kCS]) * cspol2Tracking[ipt][kCS] *(1.0+csMCkine[ipt][0]));
            else
               spectraSystSig->SetPointError(currentDataPoint, ptErr, (csTracking[kLowPtWideBin][kCSsig]/csTracking[kLowPtWideBin][kCS]) * csTracking[ipt][kCS]*(1.0+csMCkine[ipt][0]));
         }
         // otherwise use the one for 5-10 GeV/c
         else
            spectraSystSig->SetPointError(currentDataPoint, ptErr, (csTracking[kHighPtWideBin][kCSsig]/csTracking[kHighPtWideBin][kCS]) * csTracking[ipt][kCS] *(1.0+csMCkine[ipt][0]));
      }
      
      // PID syst uncertainty ------------------------------------------------------------------------
      if(ptBinLims[ipt][0] < 0.0 && ptBinLims[ipt][1] < 0.0) {
         spectraSystPIDIntegrated->SetPoint(0, pt, csTracking[ipt][kCS]*(1.0+csMCkine[ipt][0]));
         spectraSystPIDIntegrated->SetPointError(0, ptErr, csPID[ipt][kCSsyst]/csPID[ipt][kCS] * csTracking[ipt][kCS] *(1.0+csMCkine[ipt][0]));
         cout << "PID cut var syst :: " << 100.0*csPID[ipt][kCSsyst]/csPID[ipt][kCS] << " %" << endl;
      }
      else {
         if(ptBinLims[ipt][1] < 1.301) {
            spectraSystPID->SetPoint(currentDataPoint, pt, cspol2Tracking[ipt][kCS] *(1.0+csMCkine[ipt][0]));
            spectraSystPID->SetPointError(currentDataPoint, ptErr, cspol2PID[ipt][kCSsyst]/cspol2PID[ipt][kCS] * cspol2Tracking[ipt][kCS] *(1.0+csMCkine[ipt][0]));
            cout << "PID cut var syst :: " << 100.0*cspol2PID[ipt][kCSsyst]/cspol2PID[ipt][kCS] << " %" << endl;
         }
         else {
            spectraSystPID->SetPoint(currentDataPoint, pt, csTracking[ipt][kCS]*(1.0+csMCkine[ipt][0]));
            spectraSystPID->SetPointError(currentDataPoint, ptErr, csPID[ipt][kCSsyst]/csPID[ipt][kCS]*csTracking[ipt][kCS]*(1.0+csMCkine[ipt][0]));   
            cout << "PID cut var syst :: " << 100.0*csPID[ipt][kCSsyst]/csPID[ipt][kCS] << " %" << endl;
         }
      }
      
      Double_t standardIntegrated = 0.0, standard = 0.0;
      if(ptBinLims[ipt][0] < 0.0 && ptBinLims[ipt][1] < 0.0) {
         Double_t standardDataIntegrated = (1.0/standardPidSystHistIntegrated[0]->GetBinContent(1));
         Double_t standardMcIntegrated = (1.0/standardPidSystHistIntegrated[1]->GetBinContent(1));
         standardIntegrated = TMath::Abs((standardDataIntegrated-standardMcIntegrated)/(standardDataIntegrated+standardMcIntegrated));
         spectraSystPIDstandardIntegrated->SetPoint(0, pt, csTracking[ipt][kCS] *(1.0+csMCkine[ipt][0]));
         spectraSystPIDstandardIntegrated->SetPointError(0, ptErr, TMath::Abs(standardIntegrated)*csTracking[ipt][kCS] *(1.0+csMCkine[ipt][0]));
         cout << "PID standard :: " << 100.0*TMath::Abs(standardIntegrated) << " %" << endl;
      }
      else {
         Double_t standardData = (1.0/standardPidSystHist[0]->GetBinContent(standardPidSystHist[0]->FindBin(pt)));
         Double_t standardMC = (1.0/standardPidSystHist[1]->GetBinContent(standardPidSystHist[1]->FindBin(pt)));
         standard = TMath::Abs((standardData-standardMC)/(standardData+standardMC));
         if(ptBinLims[ipt][1] < 1.301) {
            spectraSystPIDstandard->SetPoint(currentDataPoint, pt, cspol2Tracking[ipt][kCS] *(1.0+csMCkine[ipt][0]));
            spectraSystPIDstandard->SetPointError(currentDataPoint, ptErr, TMath::Abs(standard)*cspol2Tracking[ipt][kCS]*(1.0+csMCkine[ipt][0]));
         }
         else {
            spectraSystPIDstandard->SetPoint(currentDataPoint, pt, csTracking[ipt][kCS] *(1.0+csMCkine[ipt][0]));
            spectraSystPIDstandard->SetPointError(currentDataPoint, ptErr, TMath::Abs(standard)*csTracking[ipt][kCS] *(1.0+csMCkine[ipt][0]));
         }
         cout << "PID standard :: " << 100.0*TMath::Abs(standard) << " %" << endl;
      }
      
      Double_t pionIntegrated = 0.0, pion = 0.0;
      if(ptBinLims[ipt][0] < 0.0 && ptBinLims[ipt][1] < 0.0) {
         Double_t pion1dataIntegrated = (1.0/pionPidSystHistIntegrated[0][0]->GetBinContent(1));
         Double_t pion1mcIntegrated = (1.0/pionPidSystHistIntegrated[0][1]->GetBinContent(1));
         Double_t pion1Integrated = TMath::Abs((pion1dataIntegrated-pion1mcIntegrated)/(pion1dataIntegrated+pion1mcIntegrated));
         Double_t pion2dataIntegrated = (1.0/pionPidSystHistIntegrated[1][0]->GetBinContent(1));
         Double_t pion2mcIntegrated = (1.0/pionPidSystHistIntegrated[1][1]->GetBinContent(1));
         Double_t pion2Integrated = TMath::Abs((pion2dataIntegrated-pion2mcIntegrated)/(pion2dataIntegrated+pion2mcIntegrated));
         pionIntegrated = 0.5*(pion1Integrated+pion2Integrated);
         spectraSystPIDpionIntegrated->SetPoint(0, pt, csTracking[ipt][kCS] *(1.0+csMCkine[ipt][0]));
         spectraSystPIDpionIntegrated->SetPointError(0, ptErr, TMath::Abs(pionIntegrated)*csTracking[ipt][kCS] *(1.0+csMCkine[ipt][0]));
         cout << "pion :: " << pionIntegrated*100.0 << " %" << endl;
      }
      else {
         Double_t pion1data = (1.0/pionPidSystHist[0][0]->GetBinContent(pionPidSystHist[0][0]->FindBin(pt)));
         Double_t pion1mc = (1.0/pionPidSystHist[0][1]->GetBinContent(pionPidSystHist[0][1]->FindBin(pt)));
         Double_t pion1 = TMath::Abs((pion1data-pion1mc)/(pion1data+pion1mc));
         Double_t pion2data = (1.0/pionPidSystHist[1][0]->GetBinContent(pionPidSystHist[1][0]->FindBin(pt)));
         Double_t pion2mc = (1.0/pionPidSystHist[1][1]->GetBinContent(pionPidSystHist[1][1]->FindBin(pt)));
         Double_t pion2 = TMath::Abs((pion2data-pion2mc)/(pion2data+pion2mc));
         pion = 0.5*(pion1+pion2);
         if(ptBinLims[ipt][1] < 1.301) {
            spectraSystPIDpion->SetPoint(currentDataPoint, pt, cspol2Tracking[ipt][kCS] *(1.0+csMCkine[ipt][0]));
            spectraSystPIDpion->SetPointError(currentDataPoint, ptErr, TMath::Abs(pion)*cspol2Tracking[ipt][kCS]*(1.0+csMCkine[ipt][0]));
         }
         else {
            spectraSystPIDpion->SetPoint(currentDataPoint, pt, csTracking[ipt][kCS] *(1.0+csMCkine[ipt][0]));
            spectraSystPIDpion->SetPointError(currentDataPoint, ptErr, TMath::Abs(pion)*csTracking[ipt][kCS] *(1.0+csMCkine[ipt][0]));
         }
         cout << "pion :: " << pion*100.0 << " %" << endl;
      }

      Double_t protIntegrated = 0.0, prot = 0.0;
      if(ptBinLims[ipt][0] < 0.0 && ptBinLims[ipt][1] < 0.0) {
         Double_t prot1dataIntegrated = (1.0/protPidSystHistIntegrated[0][0]->GetBinContent(1));
         Double_t prot1mcIntegrated = (1.0/protPidSystHistIntegrated[0][1]->GetBinContent(1));
         Double_t prot1Integrated = TMath::Abs((prot1dataIntegrated-prot1mcIntegrated)/(prot1dataIntegrated+prot1mcIntegrated));
         Double_t prot2dataIntegrated = (1.0/protPidSystHistIntegrated[1][0]->GetBinContent(1));
         Double_t prot2mcIntegrated = (1.0/protPidSystHistIntegrated[1][1]->GetBinContent(1));
         Double_t prot2Integrated = TMath::Abs((prot2dataIntegrated-prot2mcIntegrated)/(prot2dataIntegrated+prot2mcIntegrated));
         protIntegrated = 0.5*(prot1Integrated+prot2Integrated);
         spectraSystPIDprotIntegrated->SetPoint(0, pt, csTracking[ipt][kCS] *(1.0+csMCkine[ipt][0]));
         spectraSystPIDprotIntegrated->SetPointError(0, ptErr, TMath::Abs(protIntegrated)*csTracking[ipt][kCS] *(1.0+csMCkine[ipt][0]));
         cout << "proton :: " << protIntegrated*100.0 << " %" << endl;
      }
      else {
         Double_t prot1data = (1.0/protPidSystHist[0][0]->GetBinContent(protPidSystHist[0][0]->FindBin(pt)));
         Double_t prot1mc = (1.0/protPidSystHist[0][1]->GetBinContent(protPidSystHist[0][1]->FindBin(pt)));
         Double_t prot1 = TMath::Abs((prot1data-prot1mc)/(prot1data+prot1mc));
         Double_t prot2data = (1.0/protPidSystHist[1][0]->GetBinContent(protPidSystHist[1][0]->FindBin(pt)));
         Double_t prot2mc = (1.0/protPidSystHist[1][1]->GetBinContent(protPidSystHist[1][1]->FindBin(pt)));
         Double_t prot2 = TMath::Abs((prot2data-prot2mc)/(prot2data+prot2mc));
         prot = 0.5*(prot1+prot2);
         if(ptBinLims[ipt][1] < 1.301) {
            spectraSystPIDprot->SetPoint(currentDataPoint, pt, cspol2Tracking[ipt][kCS]*(1.0+csMCkine[ipt][0]));
            spectraSystPIDprot->SetPointError(currentDataPoint, ptErr, TMath::Abs(prot)*cspol2Tracking[ipt][kCS]*(1.0+csMCkine[ipt][0]));
         }
         else {
            spectraSystPIDprot->SetPoint(currentDataPoint, pt, csTracking[ipt][kCS] *(1.0+csMCkine[ipt][0]));
            spectraSystPIDprot->SetPointError(currentDataPoint, ptErr, TMath::Abs(prot)*csTracking[ipt][kCS]*(1.0+csMCkine[ipt][0]));
         }
         cout << "proton :: " << prot*100.0 << " %" << endl;
      }

      Double_t eleIntegrated = 0.0, ele = 0.0;
      if(ptBinLims[ipt][0] < 0.0 && ptBinLims[ipt][1] < 0.0) {
         Double_t ele1dataIntegrated = (1.0/elePidSystHistIntegrated[0][0]->GetBinContent(1));
         Double_t ele1mcIntegrated = (1.0/elePidSystHistIntegrated[0][1]->GetBinContent(1));
         Double_t ele1Integrated = TMath::Abs((ele1dataIntegrated-ele1mcIntegrated)/(ele1dataIntegrated+ele1mcIntegrated));
         Double_t ele2dataIntegrated = (1.0/elePidSystHistIntegrated[1][0]->GetBinContent(1));
         Double_t ele2mcIntegrated = (1.0/elePidSystHistIntegrated[1][1]->GetBinContent(1));
         Double_t ele2Integrated = TMath::Abs((ele2dataIntegrated-ele2mcIntegrated)/(ele2dataIntegrated+ele2mcIntegrated));
         Double_t ele3dataIntegrated = (1.0/elePidSystHistIntegrated[2][0]->GetBinContent(1));
         Double_t ele3mcIntegrated = (1.0/elePidSystHistIntegrated[2][1]->GetBinContent(1));
         Double_t ele3Integrated = TMath::Abs((ele3dataIntegrated-ele3mcIntegrated)/(ele3dataIntegrated+ele3mcIntegrated));
         Double_t ele4dataIntegrated = (1.0/elePidSystHistIntegrated[3][0]->GetBinContent(1)); 
         Double_t ele4mcIntegrated = (1.0/elePidSystHistIntegrated[3][1]->GetBinContent(1));
         Double_t ele4Integrated = TMath::Abs((ele4dataIntegrated-ele4mcIntegrated)/(ele4dataIntegrated+ele4mcIntegrated));
         eleIntegrated = 0.25*(ele1Integrated+ele2Integrated+ele3Integrated+ele4Integrated);
         spectraSystPIDeleIntegrated->SetPoint(0, pt, csTracking[ipt][kCS] *(1.0+csMCkine[ipt][0]));
         spectraSystPIDeleIntegrated->SetPointError(0, ptErr, TMath::Abs(eleIntegrated)*csTracking[ipt][kCS] *(1.0+csMCkine[ipt][0]));
         cout << "electron :: " << eleIntegrated*100.0 << " %" << endl;
      }
      else {
         Double_t ele1data = (1.0/elePidSystHist[0][0]->GetBinContent(elePidSystHist[0][0]->FindBin(pt)));
         Double_t ele1mc = (1.0/elePidSystHist[0][1]->GetBinContent(elePidSystHist[0][1]->FindBin(pt)));
         Double_t ele1 = TMath::Abs((ele1data-ele1mc)/(ele1data+ele1mc));
         Double_t ele2data = (1.0/elePidSystHist[1][0]->GetBinContent(elePidSystHist[1][0]->FindBin(pt)));
         Double_t ele2mc = (1.0/elePidSystHist[1][1]->GetBinContent(elePidSystHist[1][1]->FindBin(pt)));
         Double_t ele2 = TMath::Abs((ele2data-ele2mc)/(ele2data+ele2mc));
         Double_t ele3data = (1.0/elePidSystHist[2][0]->GetBinContent(elePidSystHist[2][0]->FindBin(pt)));
         Double_t ele3mc = (1.0/elePidSystHist[2][1]->GetBinContent(elePidSystHist[2][1]->FindBin(pt)));
         Double_t ele3 = TMath::Abs((ele3data-ele3mc)/(ele3data+ele3mc));
         Double_t ele4data = (1.0/elePidSystHist[3][0]->GetBinContent(elePidSystHist[3][0]->FindBin(pt))); 
         Double_t ele4mc = (1.0/elePidSystHist[3][1]->GetBinContent(elePidSystHist[3][1]->FindBin(pt)));
         Double_t ele4 = TMath::Abs((ele4data-ele4mc)/(ele4data+ele4mc));
         ele = 0.25*(ele1+ele2+ele3+ele4);
         if(ptBinLims[ipt][1] < 1.301) {
            spectraSystPIDele->SetPoint(currentDataPoint, pt, cspol2Tracking[ipt][kCS] *(1.0+csMCkine[ipt][0]));
            spectraSystPIDele->SetPointError(currentDataPoint, ptErr, TMath::Abs(ele)*cspol2Tracking[ipt][kCS] *(1.0+csMCkine[ipt][0]));
         }
         else {
            spectraSystPIDele->SetPoint(currentDataPoint, pt, csTracking[ipt][kCS] *(1.0+csMCkine[ipt][0]));
            spectraSystPIDele->SetPointError(currentDataPoint, ptErr, TMath::Abs(ele)*csTracking[ipt][kCS]*(1.0+csMCkine[ipt][0]));
         }
         cout << "electron :: " << ele*100.0 << " %" << endl;
      }

      if(ptBinLims[ipt][0] < 0.0 && ptBinLims[ipt][1] < 0.0) {
         spectraSystPIDsingleEleIntegrated->SetPoint(0, pt, csTracking[ipt][kCS] *(1.0+csMCkine[ipt][0]));
         // quadrature of pion/proton/electron cuts
         spectraSystPIDsingleEleIntegrated->SetPointError(0, ptErr, TMath::Sqrt(standardIntegrated*standardIntegrated+pionIntegrated*pionIntegrated+protIntegrated*protIntegrated+eleIntegrated*eleIntegrated)*csTracking[ipt][kCS] *(1.0+csMCkine[ipt][0]));
      }
      else {
         if(ptBinLims[ipt][1] < 1.301) {
            spectraSystPIDsingleEle->SetPoint(currentDataPoint, pt, cspol2Tracking[ipt][kCS] *(1.0+csMCkine[ipt][0]));
            spectraSystPIDsingleEle->SetPointError(currentDataPoint, ptErr, TMath::Sqrt(standard*standard+pion*pion+prot*prot+ele*ele)*cspol2Tracking[ipt][kCS] *(1.0+csMCkine[ipt][0]));
         }
         else {
            spectraSystPIDsingleEle->SetPoint(currentDataPoint, pt, csTracking[ipt][kCS] *(1.0+csMCkine[ipt][0]));
            spectraSystPIDsingleEle->SetPointError(currentDataPoint, ptErr, TMath::Sqrt(standard*standard+pion*pion+prot*prot+ele*ele)*csTracking[ipt][kCS] *(1.0+csMCkine[ipt][0]));
         }
      }

      Double_t cs, stat, trackingSyst, fitSyst, sigSyst, pidSyst, systUncorrelated, systCorrelated, systTotal;
      Double_t trackingSystCutVariations, trackingSystITSTPCmatch;
      Double_t mcKineShift, mcKineSyst;
      if(ptBinLims[ipt][0] < 0.0 && ptBinLims[ipt][1] < 0.0) {
         // setup the correlated and uncorrelated syst error: tracking + fit
         cs = spectraStatIntegrated->GetY()[0];
         stat = spectraStatIntegrated->GetErrorY(0);
         trackingSyst = spectraSystTrackingIntegrated->GetErrorY(0);
         trackingSystCutVariations = spectraSystTrackingCutVariationsIntegrated->GetErrorY(0);
         trackingSystITSTPCmatch = spectraSystTrackingITSTPCmatchIntegrated->GetErrorY(0);
         fitSyst = spectraSystFitIntegrated->GetErrorY(0);
         sigSyst = spectraSystSigIntegrated->GetErrorY(0);
         //Double_t pidSyst = spectraSystPID->GetErrorY(currentDataPoint);
         pidSyst = spectraSystPIDsingleEleIntegrated->GetErrorY(0);
         mcKineSyst = spectraSystMCkineIntegrated->GetErrorY(0);
         mcKineShift = csMCkine[ipt][0];
         systUncorrelated = TMath::Sqrt(fitSyst*fitSyst);
         systCorrelated = TMath::Sqrt(trackingSyst*trackingSyst+sigSyst*sigSyst+pidSyst*pidSyst+mcKineSyst*mcKineSyst);
         systTotal = TMath::Sqrt(systCorrelated*systCorrelated + systUncorrelated*systUncorrelated);
         
         spectraSystUncorrelatedIntegrated->SetPoint(0, pt, spectraStatIntegrated->GetY()[0]);
         spectraSystUncorrelatedIntegrated->SetPointError(0, ptErr, systUncorrelated);
         spectraSystCorrelatedIntegrated->SetPoint(0, pt, spectraStatIntegrated->GetY()[0]);
         spectraSystCorrelatedIntegrated->SetPointError(0, ptErr, systCorrelated);
         spectraSystTotalIntegrated->SetPoint(0, pt, spectraStatIntegrated->GetY()[0]);
         spectraSystTotalIntegrated->SetPointError(0, ptErr, systTotal);
      }
      else {
         // setup the correlated and uncorrelated syst error: tracking + fit
         cs = spectraStat->GetY()[currentDataPoint];
         stat = spectraStat->GetErrorY(currentDataPoint);
         trackingSyst = spectraSystTracking->GetErrorY(currentDataPoint);
         trackingSystCutVariations = spectraSystTrackingCutVariations->GetErrorY(currentDataPoint);
         trackingSystITSTPCmatch = spectraSystTrackingITSTPCmatch->GetErrorY(currentDataPoint);
         fitSyst = spectraSystFit->GetErrorY(currentDataPoint);
         sigSyst = spectraSystSig->GetErrorY(currentDataPoint);
         //Double_t pidSyst = spectraSystPID->GetErrorY(currentDataPoint);
         pidSyst = spectraSystPIDsingleEle->GetErrorY(currentDataPoint);
         mcKineSyst = spectraSystMCkine->GetErrorY(currentDataPoint);
         mcKineShift = csMCkine[ipt][0];
         systUncorrelated = TMath::Sqrt(fitSyst*fitSyst);
         systCorrelated = TMath::Sqrt(trackingSyst*trackingSyst+sigSyst*sigSyst+pidSyst*pidSyst);
         systTotal = TMath::Sqrt(systCorrelated*systCorrelated + systUncorrelated*systUncorrelated);
         
         spectraSystUncorrelated->SetPoint(currentDataPoint, pt, spectraStat->GetY()[currentDataPoint]);
         spectraSystUncorrelated->SetPointError(currentDataPoint, ptErr, systUncorrelated);
         spectraSystCorrelated->SetPoint(currentDataPoint, pt, spectraStat->GetY()[currentDataPoint]);
         spectraSystCorrelated->SetPointError(currentDataPoint, ptErr, systCorrelated);
         spectraSystTotal->SetPoint(currentDataPoint, pt, spectraStat->GetY()[currentDataPoint]);
         spectraSystTotal->SetPointError(currentDataPoint, ptErr, systTotal);
      }
      
      cout << "cs / stat / syst.total / syst.uncorr / syst.corr :: "
           << cs << " / " << 100.0*stat/cs << " / " << 100.0*systTotal/cs << " / " << 100.0*systUncorrelated/cs << " / " << 100.0*systCorrelated/cs << endl;
      cout << "tracking / trackingCV / trackingMatch / pid / sig / fit :: " << 100.0*trackingSyst/cs << " / " << 100.0*trackingSystCutVariations/cs << " / " 
           << 100.0*trackingSystITSTPCmatch/cs << " / " << 100.0*pidSyst/cs << " / " << 100.0*sigSyst/cs << " / " << 100.0*fitSyst/cs <<  endl;
      cout << "mcKineShift / mcKineSyst :: " << mcKineShift*100.0 << " / " << 100.0*mcKineSyst/cs << " / " << 100.0*csMCkine[ipt][1] << endl;
            
      if(ptBinLims[ipt][0] > -0.01 && ptBinLims[ipt][1] > -0.01) 
         currentDataPoint++;
   }
   
   TFile* output = new TFile(outfilename, "RECREATE");
   spectraStat->Write();
   spectraSystTracking->Write();
   spectraSystTrackingCutVariations->Write();
   spectraSystTrackingITSTPCmatch->Write();
   spectraSystFit->Write();
   spectraSystSig->Write();
   spectraSystPID->Write();
   spectraSystPIDpion->Write();
   spectraSystPIDprot->Write();
   spectraSystPIDele->Write();
   spectraSystPIDstandard->Write();
   spectraSystPIDsingleEle->Write();
   spectraSystMCkine->Write();
   spectraSystUncorrelated->Write();
   spectraSystCorrelated->Write();
   spectraSystTotal->Write();
   spectraStatIntegrated->Write();
   spectraSystTrackingIntegrated->Write();
   spectraSystTrackingCutVariationsIntegrated->Write();
   spectraSystTrackingITSTPCmatchIntegrated->Write();
   spectraSystFitIntegrated->Write();
   spectraSystSigIntegrated->Write();
   spectraSystPIDIntegrated->Write();
   spectraSystPIDpionIntegrated->Write();
   spectraSystPIDprotIntegrated->Write();
   spectraSystPIDeleIntegrated->Write();
   spectraSystPIDstandardIntegrated->Write();
   spectraSystPIDsingleEleIntegrated->Write();
   spectraSystMCkineIntegrated->Write();
   spectraSystUncorrelatedIntegrated->Write();
   spectraSystCorrelatedIntegrated->Write();
   spectraSystTotalIntegrated->Write();
   output->Close();
}



//________________________________________________________________________________________
void GetCrossSection(const Char_t* filename, const Char_t* filenameOut, Double_t* csReturn = 0x0) {
   
   TFile* resultsFile = TFile::Open(filename);
   TH1D* activeSettings = (TH1D*)resultsFile->Get("ActiveSettings");
   TH1D* activeSettingsRefilter = (TH1D*)activeSettings->Clone("ActiveSettingsRefilter");
   activeSettingsRefilter->SetDirectory(0x0);
   activeSettingsRefilter->Reset();
   resultsFile->Close();
   
   for(Int_t iSetting = 1; iSetting<activeSettingsRefilter->GetXaxis()->GetNbins(); ++iSetting) {
      TString binLabel = activeSettingsRefilter->GetXaxis()->GetBinLabel(iSetting);
      activeSettingsRefilter->SetBinContent(iSetting, 1.0);
      if(binLabel.Contains("StandardCommon")) activeSettingsRefilter->SetBinContent(iSetting, 0.0);
      if(binLabel.Contains("ITS4cls")) activeSettingsRefilter->SetBinContent(iSetting, 0.0);
      if(binLabel.Contains("ITS3cls")) activeSettingsRefilter->SetBinContent(iSetting, 0.0);
      if(binLabel.Contains("SPDboth")) activeSettingsRefilter->SetBinContent(iSetting, 0.0);
   }
   
   TList* outList = RecomputeCrossSection(filename, activeSettingsRefilter);
      
   TH1D* csDistrib = (TH1D*)outList->At(0);
   TH1D* csStatDistrib = (TH1D*)outList->At(1);
   TH1D* fitRMSDistrib = (TH1D*)outList->At(4);
   TH1D* sigRMSDistrib = (TH1D*)outList->At(7);
   
   Double_t cs[5];
   cs[0] = csDistrib->GetMean();
   cs[1] = csStatDistrib->GetMean();
   cs[2] = csDistrib->GetRMS();
   cs[3] = fitRMSDistrib->GetMean();
   cs[4] = sigRMSDistrib->GetMean();
     
   
   cout << "cs = " << cs[0] << " +/- " 
        << cs[1] << " (stat. " << 100.0*cs[1]/cs[0] << "%) +/- " 
        << cs[2] << " (syst. " << 100.0*cs[2]/cs[0] <<"%) +/- " 
        << cs[3] << " (syst.fit " << 100.0*cs[3]/cs[0] << "%)  +/- "
        << cs[4] << " (syst.sig " << 100.0*cs[4]/cs[0] << "%)  nb" <<  endl;
   
   
   if(csReturn) {
      for(Int_t i = 0; i<5; ++i) csReturn[i] = cs[i];
   }
   
   TFile* refilteredResults = new TFile(filenameOut, "RECREATE");
   activeSettingsRefilter->Write();
   outList->Write();
   refilteredResults->Close();
}


//_________________________________________________________________________________________
TList* RecomputeCrossSection(const Char_t* filename, TH1D* activeSettings) {
   // 
   // open the Results.root file and recompute the mean cross-section and systematics according to the activeSettings
   // 
   TFile* resultsFile = TFile::Open(filename);
   TH2D* cs2D = (TH2D*)resultsFile->Get("CrossSection2D");
      
   if(cs2D->GetXaxis()->GetNbins() != activeSettings->GetXaxis()->GetNbins()) {
      cout << "RecomputeCrossSection() ERROR: Active settings histograms for tracking uncert does not have the same number of bins as the histograms in the Result.root file" <<  endl;
      return 0x0;
   }
   
   // make a prerun over settings to find limits for histograms
   Double_t minCS = -1.0; Double_t maxCS = -1.0;            // min and maximum cross-section among the tracking variations
   Double_t minCSstat = -1.0; Double_t maxCSstat = -1.0;    // min and maximum cross-section stat error among the tracking variations
   Double_t minCSFit = -1.0; Double_t maxCSFit = -1.0;      // min and maximum cross-section among all the tracking X fit variations
   Double_t minCSSig = -1.0; Double_t maxCSSig = -1.0;      // min and maximum cross-section among all the tracking X sig variations
   for(Int_t iSetting = 1; iSetting <= cs2D->GetXaxis()->GetNbins(); ++iSetting) {
      if(activeSettings->GetBinContent(iSetting) < 1.0e-5) continue;
      Double_t cs = cs2D->GetBinContent(iSetting, 1);
      if(minCS<0.0) minCS = cs;
      if(maxCS<0.0) maxCS = cs;
      if(cs<minCS) minCS = cs;
      if(cs>maxCS) maxCS = cs;
         
      Double_t csStat = cs2D->GetBinError(iSetting, 1);
      if(minCSstat<0.0) minCSstat = csStat;
      if(maxCSstat<0.0) maxCSstat = csStat;
      if(csStat<minCSstat) minCSstat = csStat;
      if(csStat>maxCSstat) maxCSstat = csStat;
      
      for(Int_t j = 1; j <= cs2D->GetYaxis()->GetNbins(); ++j) {
         cs = cs2D->GetBinContent(iSetting, j);
         if(j == 1 && iSetting == 1) {
            minCSFit = cs;
            maxCSFit = cs;
            minCSSig = cs;
            maxCSSig = cs;
            continue;
         }
         TString binLabel = cs2D->GetYaxis()->GetBinLabel(j);
         if(binLabel.Contains("mSig")) {
            if(cs<minCSSig) minCSSig = cs;
            if(cs>maxCSSig) maxCSSig = cs;
         }
         if(binLabel.Contains("mFit")) {
            if(cs<minCSFit) minCSFit = cs;
            if(cs>maxCSFit) maxCSFit = cs;
         }
      }
   }
   
   TH1D* csDistrib = new TH1D("csDistrib", "Cross section distribution, cut variations", 100, minCS-1*(maxCS-minCS), maxCS+1*(maxCS-minCS));
   csDistrib->GetXaxis()->SetTitle(cs2D->GetZaxis()->GetTitle());
   csDistrib->SetDirectory(0x0);
   TH1D* csStatDistrib = new TH1D("csStatDistrib", "Cross section stat err distribution, cut variations", 100, minCSstat-1*(maxCSstat-minCSstat), maxCSstat+1*(maxCSstat-minCSstat));
   csStatDistrib->GetXaxis()->SetTitle(cs2D->GetZaxis()->GetTitle());
   csStatDistrib->SetDirectory(0x0);
   TH2D* csDistribFit2D = new TH2D("csDistribFit2D", "Cross section distribution, cut x fit variations", cs2D->GetXaxis()->GetNbins(), 0.0, cs2D->GetXaxis()->GetNbins(), 100, minCSFit-1*(maxCSFit-minCSFit), maxCSFit+1*(maxCSFit-minCSFit));
   csDistribFit2D->GetYaxis()->SetTitle(cs2D->GetZaxis()->GetTitle());
   csDistribFit2D->SetDirectory(0x0);
   TH2D* csDistribSig2D = new TH2D("csDistribSig2D", "Cross section distribution, cut x sig variations", cs2D->GetXaxis()->GetNbins(), 0.0, cs2D->GetXaxis()->GetNbins(), 100, minCSSig-1*(maxCSSig-minCSSig), maxCSSig+1*(maxCSSig-minCSSig));
   csDistribSig2D->GetYaxis()->SetTitle(cs2D->GetZaxis()->GetTitle());
   csDistribSig2D->SetDirectory(0x0);
   TH1D* csFitMeanRMS = new TH1D("csFitMeanRMS", "Mean and RMS for fit variations vs cut setting", 
                                         cs2D->GetXaxis()->GetNbins(), 0.0, cs2D->GetXaxis()->GetNbins());
   csFitMeanRMS->SetDirectory(0x0);
   TH1D* csSigMeanRMS = new TH1D("csSigMeanRMS", "Mean and RMS for sig variations vs cut setting", 
      cs2D->GetXaxis()->GetNbins(), 0.0, cs2D->GetXaxis()->GetNbins());
   csSigMeanRMS->SetDirectory(0x0);
                                         
   Double_t minRMS = -1.0; Double_t maxRMS = -1.0;
   Double_t minRMSsig = -1.0; Double_t maxRMSsig = -1.0;
   for(Int_t iSetting = 1; iSetting <= cs2D->GetXaxis()->GetNbins(); ++iSetting) {
      if(activeSettings->GetBinContent(iSetting) < 1.0e-5) continue;
      Double_t cs = cs2D->GetBinContent(iSetting, 1);
      csDistrib->Fill(cs);
      Double_t csStatErr = cs2D->GetBinError(iSetting, 1);
      csStatDistrib->Fill(csStatErr);
      csDistribSig2D->GetXaxis()->SetBinLabel(iSetting, cs2D->GetXaxis()->GetBinLabel(iSetting));
      csSigMeanRMS->GetXaxis()->SetBinLabel(iSetting, cs2D->GetXaxis()->GetBinLabel(iSetting));
      csDistribFit2D->GetXaxis()->SetBinLabel(iSetting, cs2D->GetXaxis()->GetBinLabel(iSetting));
      csFitMeanRMS->GetXaxis()->SetBinLabel(iSetting, cs2D->GetXaxis()->GetBinLabel(iSetting));
      
      TH1D* tempCS = new TH1D("tempCS", "", 100, minCSFit-1*(maxCSFit-minCSFit), maxCSFit+1*(maxCSFit-minCSFit));
      TH1D* tempCSsig = new TH1D("tempCSsig", "", 100, minCSSig-1*(maxCSSig-minCSSig), maxCSSig+1*(maxCSSig-minCSSig));
      for(Int_t j = 1; j <= cs2D->GetYaxis()->GetNbins(); ++j) {
         cs = cs2D->GetBinContent(iSetting, j);
         if(j == 1) {
            csDistribFit2D->Fill(-0.5+iSetting, cs);
            tempCS->Fill(cs);
            csDistribSig2D->Fill(-0.5+iSetting, cs);
            tempCSsig->Fill(cs);
            continue;
         }
         TString binLabel = cs2D->GetYaxis()->GetBinLabel(j);
         if(binLabel.Contains("mSig")) {
            csDistribSig2D->Fill(-0.5+iSetting, cs);
            tempCSsig->Fill(cs);
            continue;
         }
         if(binLabel.Contains("mFit")) {
            csDistribFit2D->Fill(-0.5+iSetting, cs);
            tempCS->Fill(cs);
         }
      }
      if(iSetting == 1) {
         minRMS = tempCS->GetRMS();
         maxRMS = tempCS->GetRMS();
         minRMSsig = tempCSsig->GetRMS();
         maxRMSsig = tempCSsig->GetRMS();
      }
      else {
         if(tempCS->GetRMS()<minRMS) minRMS = tempCS->GetRMS();
         if(tempCS->GetRMS()>maxRMS) maxRMS = tempCS->GetRMS();
         if(tempCSsig->GetRMS()<minRMSsig) minRMSsig = tempCSsig->GetRMS();
         if(tempCSsig->GetRMS()>maxRMSsig) maxRMSsig = tempCSsig->GetRMS();
      }
      
      csFitMeanRMS->SetBinContent(iSetting, tempCS->GetMean());
      csFitMeanRMS->SetBinError(iSetting, tempCS->GetRMS());
      csSigMeanRMS->SetBinContent(iSetting, tempCSsig->GetMean());
      csSigMeanRMS->SetBinError(iSetting, tempCSsig->GetRMS());
      
      delete tempCS;
      delete tempCSsig;
   }
   
   TH1D* csDistribFit = new TH1D("csDistribFit", "Distribution of RMS's from fit variations, from all cut variations", 100, minRMS-1*(maxRMS-minRMS), maxRMS+1*(maxRMS-minRMS));
   csDistribFit->GetXaxis()->SetTitle(cs2D->GetZaxis()->GetTitle());
   csDistribFit->SetDirectory(0x0);
   TH1D* csDistribSig = new TH1D("csDistribSig", "Distribution of RMS's from sig variations, from all cut variations", 100, minRMSsig-1*(maxRMSsig-minRMSsig), maxRMSsig+1*(maxRMSsig-minRMSsig));
   csDistribSig->GetXaxis()->SetTitle(cs2D->GetZaxis()->GetTitle());
   csDistribSig->SetDirectory(0x0);
   for(Int_t iSetting = 1; iSetting <= cs2D->GetXaxis()->GetNbins(); ++iSetting) {
      if(activeSettings->GetBinContent(iSetting) < 1.0e-5) continue;
      csDistribFit->Fill(csFitMeanRMS->GetBinError(iSetting));
      csDistribSig->Fill(csSigMeanRMS->GetBinError(iSetting));
   }
   
   TList* outList = new TList();
   outList->SetOwner(kTRUE);
   outList->Add(csDistrib);
   outList->Add(csStatDistrib);
   outList->Add(csDistribFit2D);
   outList->Add(csFitMeanRMS);
   outList->Add(csDistribFit);
   outList->Add(csDistribSig2D);
   outList->Add(csSigMeanRMS);
   outList->Add(csDistribSig);
      
   resultsFile->Close();
   
   return outList;
}


//_________________________________________________________________________________________
void RunAllSignalExtraction(const Char_t* filename, Int_t sigExtrMethod = 1, 
                            Double_t ptMin = -1.0, Double_t ptMax = -1.0, 
                            Int_t periodCode = 6, 
                            const Char_t* outputDir = "", 
                            const Char_t* filenameMC = "", const Char_t* dirMC = "", 
                            const Char_t* mcSignalElectron = "TrueElectron", 
                            const Char_t* mcSignalJpsi = "mcTruthJpsi", 
                            Bool_t runSigExtrSyst = kFALSE) {
   // 
   // run signal extraction automatically for all settings found in the provided file
   // 
   // Signal extraction method (sigExtrMethod):
   // 1 - hybrid (ME bkg + residual fit)
   // 2 - fit (bkg fit function + signal MC shape)
   // 3 - mixed event
   // 4 - like sign
   // 
   //  Options for periodCode:  (used in the ComputeIntegratedLumi() function)
   // 0 - LHC15n
   // 1 - LHC17p
   // 2 - LHC17q
   // 3 - LHC17p + LHC17q
   // 4 - LHC17p + LHC17q (CENT_woSDD)
   // 5 - LHC17p + LHC17q (FAST)
   // 6 - LHC15n + LHC17p + LHC17q
   
   if (sigExtrMethod<1 || sigExtrMethod>4) {
      cout << "Bad or not implemented signal extraction method! The sigExtrMethod parameter must be an integer from 1 to 4,  check it out !!!" << endl;
      return;
   }
   
   Bool_t isNewFile = kFALSE;
   if(!gHistManData || std::strcmp(gCurrentDataFilename, filename)) {
      if(gHistManData) {
         gHistManData->CloseFile();
         delete gHistManData;
      }
      gHistManData = new AliHistogramManager();
      gHistManData->InitFile(filename, "jpsi2eeHistos");

      if(filenameMC[0]!='\0') {
         gHistManMC = new AliHistogramManager();
         gHistManMC->InitFile(filenameMC, "jpsi2eeHistos");
      }
      gCurrentDataFilename = filename;
      isNewFile = kTRUE;
   }
   gSystem->Exec(Form("mkdir -p %s", outputDir));
   
   
   // compute luminosity: need to access the event statistics histograms in the DST tree files
   Double_t totalLumi = 0.0; Double_t totalLumiErr = 0.0;
   TList* lumiHists = ComputeIntegratedLumi(gkDirTreesLHC17pq, gkDirTreesLHC15n, periodCode, totalLumi, totalLumiErr);
   TH1D* nINT7eventsHist = (TH1D*)lumiHists->FindObject("nINT7evHist");
   
   // find all the setting names in the data file
   THashList* mainDir = gHistManData->GetMainDirectory();
   Int_t nSettings = 0;
   TString settingNames[64];
   for(Int_t i = 0; i < mainDir->GetEntries(); ++i) {
      TString currentSubdir = mainDir->At(i)->GetName();
      //  check that this is a directory which contains SE-PM pairs: name should start with "PairSEPM"; then get the name of the cut setting
      TObjArray* arr = currentSubdir.Tokenize("_");
      TString firstSubString = arr->At(0)->GetName();
      if(firstSubString.CompareTo("PairSEPM"))              //"PairSEPM" must be the first substring 
         continue;
      settingNames[nSettings++] = currentSubdir.Remove(0, 9); // remove the "PairSEPM_" substring; what is left should be the cut setting name
   }

                              
   // Setup histograms -----------------------------------------------------------------------------------------------
   TList outputHistograms;
   outputHistograms.SetOwner(kTRUE);
   
   TH1D* rawCountsHist = new TH1D("RawCounts", "Raw counts", nSettings, 0.0, nSettings);     
   rawCountsHist->Sumw2();
   rawCountsHist->GetYaxis()->SetTitle("counts"); 
   outputHistograms.Add(rawCountsHist);
                               
   TH1D* sOverBHist = new TH1D("SoverB", "S/B", nSettings, 0.0, nSettings);            
   sOverBHist->GetYaxis()->SetTitle("S/B");
   outputHistograms.Add(sOverBHist);
   
   TH1D* signifHist = new TH1D("Significance", "Significance", nSettings, 0.0, nSettings);
   signifHist->GetYaxis()->SetTitle("Significance");
   outputHistograms.Add(signifHist);
   
   TH1D* chi2Hist = new TH1D("Chi2Matching", "#chi^{2} matching", nSettings, 0.0, nSettings); 
   chi2Hist->GetYaxis()->SetTitle("#chi^{2}/NDF");
   outputHistograms.Add(chi2Hist);   
   
   TH1D* chi2MCHist = new TH1D("Chi2MC", "#chi^{2} MC", nSettings, 0.0, nSettings); 
   chi2MCHist->GetYaxis()->SetTitle("#chi^{2}/NDF");
   outputHistograms.Add(chi2MCHist);    
   
   TH1D* effHist = new TH1D("Efficiency", "Total acc X efficiency", nSettings, 0.0, nSettings); 
   effHist->GetYaxis()->SetTitle("acc. #times eff");
   effHist->Sumw2();
   outputHistograms.Add(effHist);
   
   TH1D* sigFracHist = new TH1D("SigFraction", "Signal fraction in counting window", nSettings, 0.0, nSettings); 
   sigFracHist->GetYaxis()->SetTitle("fraction");
   outputHistograms.Add(sigFracHist);
   
   TH1D* corrYieldHist = new TH1D("CorrYield", "Corrected yield, not normalized", nSettings, 0.0, nSettings); 
   corrYieldHist->Sumw2();
   corrYieldHist->GetYaxis()->SetTitle("yield");
   outputHistograms.Add(corrYieldHist);
   
   TH1D* deltaCorrected = new TH1D("DeltaCorrected", "Y_{variation} - Y_{standard}, not normalized", nSettings, 0.0, nSettings); 
   deltaCorrected->Sumw2();
   deltaCorrected->GetYaxis()->SetTitle("deviation");
   outputHistograms.Add(deltaCorrected);
   
   TH1D* deltaCorrectedFrac = new TH1D("DeltaCorrectedFraction", "(Y_{variation} - Y_{standard})/Y_{standard}", nSettings, 0.0, nSettings); 
   deltaCorrectedFrac->Sumw2();
   deltaCorrectedFrac->GetYaxis()->SetTitle("deviation (%)");
   outputHistograms.Add(deltaCorrectedFrac);
   
   TH1D* crossSectionHist = new TH1D("CrossSection", "Cross-section (d#sigma / dy)", nSettings, 0.0, nSettings); 
   crossSectionHist->Sumw2();
   if (ptMin<0 || ptMax<0) crossSectionHist->GetYaxis()->SetTitle("d#sigma / dy (nb)");
   else crossSectionHist->GetYaxis()->SetTitle("d^{2}#sigma / dy / dp_{T} (nb)");
   outputHistograms.Add(crossSectionHist);

   TH1D* activeSettings = new TH1D("ActiveSettings", "Settings used for results", nSettings, 0.0, nSettings);
   activeSettings->GetYaxis()->SetTitle("active / inactive (1/0)");
   outputHistograms.Add(activeSettings);
   
   // define axis labels for the signal extraction variations -----------------------------------------------
   TString axLabels[gkNSigCountingVariations+gkNMassFitVariations*gkNMassExclVariations-1];
   for (Int_t i = 0;i<gkNSigCountingVariations;++i) 
      axLabels[i] = Form("mSig [%.2f ; %.2f]", gmSigLims[i][0], gmSigLims[i][1]);
   for (Int_t i = 0;i<gkNMassFitVariations;++i) {
      for (Int_t j = 0;j<gkNMassExclVariations;++j) {
         if (i == 0 && j == 0) 
            continue;
         axLabels[gkNSigCountingVariations+i*gkNMassExclVariations+j-1] = Form("mFit [%.2f ; %.2f], mExcl [%.2f ; %.2f]", 
                                                                        gmFitLims[i][0], gmFitLims[i][1], gmExclLims[j][0], gmExclLims[j][1]);
      }
   }
      
   Int_t nTotalSigExtrVariations = gkNSigCountingVariations + gkNMassFitVariations*gkNMassExclVariations-1; 
   TH2D* rawCountsHist2D = 0x0; 
   TH2D* sOverBHist2D = 0x0; 
   TH2D* signifHist2D = 0x0; 
   TH2D* chi2Hist2D = 0x0;
   TH2D* chi2MCHist2D = 0x0; 
   TH2D* effHist2D = 0x0; 
   TH2D* sigFracHist2D = 0x0; 
   TH2D* corrYieldHist2D = 0x0;
   TH2D* deltaCorrected2D = 0x0; 
   TH2D* deltaCorrectedFrac2D = 0x0; 
   TH2D* crossSectionHist2D = 0x0;
   if(runSigExtrSyst) {
      rawCountsHist2D = new TH2D("RawCounts2D", "Raw counts", nSettings, 0.0, nSettings, nTotalSigExtrVariations, -0.5, -0.5+nTotalSigExtrVariations);
      rawCountsHist2D->Sumw2();
      rawCountsHist2D->GetZaxis()->SetTitle("counts");
      outputHistograms.Add(rawCountsHist2D);
      
      sOverBHist2D = new TH2D("SoverB2D", "S/B", nSettings, 0.0, nSettings, nTotalSigExtrVariations, -0.5, -0.5+nTotalSigExtrVariations); 
      sOverBHist2D->GetZaxis()->SetTitle("S/B");
      outputHistograms.Add(sOverBHist2D);
      
      signifHist2D = new TH2D("Significance2D", "Significance", nSettings, 0.0, nSettings, nTotalSigExtrVariations, -0.5, -0.5+nTotalSigExtrVariations);
      signifHist2D->GetZaxis()->SetTitle("Significance");
      outputHistograms.Add(signifHist2D);
      
      chi2Hist2D = new TH2D("Chi2Matching2D", "#chi^{2} matching", nSettings, 0.0, nSettings, nTotalSigExtrVariations, -0.5, -0.5+nTotalSigExtrVariations); 
      chi2Hist2D->GetZaxis()->SetTitle("#chi^{2}/NDF");
      outputHistograms.Add(chi2Hist2D);
      
      chi2MCHist2D = new TH2D("Chi2MC2D", "#chi^{2} MC", nSettings, 0.0, nSettings, nTotalSigExtrVariations, -0.5, -0.5+nTotalSigExtrVariations); 
      chi2MCHist2D->GetZaxis()->SetTitle("#chi^{2}/NDF");
      outputHistograms.Add(chi2MCHist2D); 
      
      effHist2D = new TH2D("Efficiency2D", "Total acc X efficiency", nSettings, 0.0, nSettings, nTotalSigExtrVariations, -0.5, -0.5+nTotalSigExtrVariations); 
      effHist2D->Sumw2();
      effHist2D->GetZaxis()->SetTitle("acc #times eff");
      outputHistograms.Add(effHist2D);
      
      sigFracHist2D = new TH2D("SigFraction2D", "Signal fraction in counting window", nSettings, 0.0, nSettings, nTotalSigExtrVariations, -0.5, -0.5+nTotalSigExtrVariations); 
      sigFracHist2D->GetZaxis()->SetTitle("fraction");
      outputHistograms.Add(sigFracHist2D);
      
      corrYieldHist2D = new TH2D("CorrYield2D", "Corrected yield, not normalized", nSettings, 0.0, nSettings, nTotalSigExtrVariations, -0.5, -0.5+nTotalSigExtrVariations); 
      corrYieldHist2D->Sumw2();
      corrYieldHist2D->GetZaxis()->SetTitle("corrected yield");
      outputHistograms.Add(corrYieldHist2D);
      
      deltaCorrected2D = new TH2D("DeltaCorrected2D", "Y_{variation} - Y_{standard}, not normalized", nSettings, 0.0, nSettings, nTotalSigExtrVariations, -0.5, -0.5+nTotalSigExtrVariations); 
      deltaCorrected2D->Sumw2();
      deltaCorrected2D->GetZaxis()->SetTitle("deviation");
      outputHistograms.Add(deltaCorrected2D);
      
      deltaCorrectedFrac2D = new TH2D("DeltaCorrectedFraction2D", "(Y_{variation} - Y_{standard})/Y_{standard}", nSettings, 0.0, nSettings, nTotalSigExtrVariations, -0.5, -0.5+nTotalSigExtrVariations); 
      deltaCorrectedFrac2D->Sumw2();
      deltaCorrectedFrac2D->GetZaxis()->SetTitle("deviation (%)");
      outputHistograms.Add(deltaCorrectedFrac2D);
      
      crossSectionHist2D = new TH2D("CrossSection2D", "Cross-section (d#sigma / dy)", nSettings, 0.0, nSettings, nTotalSigExtrVariations, -0.5, -0.5+nTotalSigExtrVariations); 
      crossSectionHist2D->Sumw2();
      if (ptMin<0 || ptMax<0) 
         crossSectionHist2D->GetZaxis()->SetTitle("d#sigma / dy (nb)");
      else 
         crossSectionHist2D->GetZaxis()->SetTitle("d^{2}#sigma / dy / dp_{T} (nb)");
      outputHistograms.Add(crossSectionHist2D);
      
      for(Int_t i = 0;i<gkNSigCountingVariations+gkNMassFitVariations*gkNMassExclVariations;++i) {
         rawCountsHist2D->GetYaxis()->SetBinLabel(i+1, axLabels[i].Data());
         sOverBHist2D->GetYaxis()->SetBinLabel(i+1, axLabels[i].Data());
         signifHist2D->GetYaxis()->SetBinLabel(i+1, axLabels[i].Data());
         chi2Hist2D->GetYaxis()->SetBinLabel(i+1, axLabels[i].Data());
         chi2MCHist2D->GetYaxis()->SetBinLabel(i+1, axLabels[i].Data());
         effHist2D->GetYaxis()->SetBinLabel(i+1, axLabels[i].Data());
         sigFracHist2D->GetYaxis()->SetBinLabel(i+1, axLabels[i].Data());
         corrYieldHist2D->GetYaxis()->SetBinLabel(i+1, axLabels[i].Data());
         deltaCorrected2D->GetYaxis()->SetBinLabel(i+1, axLabels[i].Data());
         deltaCorrectedFrac2D->GetYaxis()->SetBinLabel(i+1, axLabels[i].Data());
         crossSectionHist2D->GetYaxis()->SetBinLabel(i+1, axLabels[i].Data());
      }
   }                                                        // end if (runSigExtrSyst)

   // set the axis labels for setting names 
   for (Int_t i = 0; i<nSettings; ++i) {
      rawCountsHist->GetXaxis()->SetBinLabel(i+1,  settingNames[i].Data());
      sOverBHist->GetXaxis()->SetBinLabel(i+1,  settingNames[i].Data());
      signifHist->GetXaxis()->SetBinLabel(i+1,  settingNames[i].Data());
      chi2Hist->GetXaxis()->SetBinLabel(i+1,  settingNames[i].Data());
      chi2MCHist->GetXaxis()->SetBinLabel(i+1,  settingNames[i].Data());
      effHist->GetXaxis()->SetBinLabel(i+1, settingNames[i].Data());
      sigFracHist->GetXaxis()->SetBinLabel(i+1, settingNames[i].Data());
      corrYieldHist->GetXaxis()->SetBinLabel(i+1, settingNames[i].Data());
      deltaCorrected->GetXaxis()->SetBinLabel(i+1, settingNames[i].Data());
      deltaCorrectedFrac->GetXaxis()->SetBinLabel(i+1, settingNames[i].Data());
      crossSectionHist->GetXaxis()->SetBinLabel(i+1, settingNames[i].Data());
      activeSettings->GetXaxis()->SetBinLabel(i+1, settingNames[i].Data());
      activeSettings->SetBinContent(i+1, 1.0);              // by default all settings are used
      if (runSigExtrSyst) {
         rawCountsHist2D->GetXaxis()->SetBinLabel(i+1,  settingNames[i].Data());
         sOverBHist2D->GetXaxis()->SetBinLabel(i+1,  settingNames[i].Data());
         signifHist2D->GetXaxis()->SetBinLabel(i+1,  settingNames[i].Data());
         chi2Hist2D->GetXaxis()->SetBinLabel(i+1,  settingNames[i].Data());
         chi2MCHist2D->GetXaxis()->SetBinLabel(i+1,  settingNames[i].Data());
         effHist2D->GetXaxis()->SetBinLabel(i+1, settingNames[i].Data());
         sigFracHist2D->GetXaxis()->SetBinLabel(i+1, settingNames[i].Data());
         corrYieldHist2D->GetXaxis()->SetBinLabel(i+1, settingNames[i].Data());
         deltaCorrected2D->GetXaxis()->SetBinLabel(i+1, settingNames[i].Data());
         deltaCorrectedFrac2D->GetXaxis()->SetBinLabel(i+1, settingNames[i].Data());
         crossSectionHist2D->GetXaxis()->SetBinLabel(i+1, settingNames[i].Data());
      }
   }
      
   // compute run-by-run weighted efficiency
   Double_t eff[64][gkNSigCountingVariations] = {{0.0}}; 
   Double_t effErr[64][gkNSigCountingVariations] = {{0.0}};
   Double_t sigFrac[64][gkNSigCountingVariations] = {{0.0}};
   TList* qaHistsForEfficiency = ComputeEfficiency(settingNames, nSettings, dirMC, nINT7eventsHist, mcSignalElectron, mcSignalJpsi, eff, effErr, sigFrac, runSigExtrSyst, ptMin, ptMax);   
   for (Int_t i = 0;i<nSettings;i++)
      cout << "eff " << settingNames[i].Data() << " " << eff[i][0] << " +/- " << effErr[i][0] << endl;
   
   TList qaHistsForSigExtraction;
   qaHistsForSigExtraction.SetOwner(kTRUE);
   
   for(Int_t i = 0; i < nSettings; ++i) {
      //if (i == 1) break;
     
      cout << "***********************************************************************************" << endl;
      cout << "Running signal extraction for setting " << settingNames[i].Data() << endl;
      if (settingNames[i].Contains("ITS4cls")) continue;
      if (settingNames[i].Contains("SPDboth")) continue;
      if (settingNames[i].Contains("ITS3cls")) continue;
      
      if (isNewFile)  SetHistograms(settingNames[i].Data(), mcSignalElectron);
         
      if(sigExtrMethod == 1) ConfigureSignalExtractionHybrid(ptMin, ptMax, gmFitLims[0][0], gmFitLims[0][1]);
      if(sigExtrMethod == 2) ConfigureSignalExtractionFit(ptMin, ptMax, gmFitLims[0][0], gmFitLims[0][1]);
      if(sigExtrMethod == 3) ConfigureSignalExtractionME(ptMin, ptMax, gmFitLims[0][0], gmFitLims[0][1], gmExclLims[0][0], gmExclLims[0][1]);
      if(sigExtrMethod == 4) ConfigureSignalExtractionLS(ptMin, ptMax, gmFitLims[0][0], gmFitLims[0][1], gmExclLims[0][0], gmExclLims[0][1]);
         
      gFits.Process();
      // add QA histograms to the output list
      AddSigExtrQAtoOutputList(&qaHistsForSigExtraction, settingNames[i].Data(), gmFitLims[0][0], gmFitLims[0][1], gmExclLims[0][0], gmExclLims[0][1]);
      
      for(Int_t isig = 0; isig<(runSigExtrSyst ? gkNSigCountingVariations : 1); ++isig) {
         gFits.ComputeOutputValues(gmSigLims[isig][0], gmSigLims[isig][1]);
         
         if(isig == 0) {
            Draw(&gFits, gmSigLims[isig][0], gmSigLims[isig][1], kTRUE, outputDir, settingNames[i].Data(), kFALSE);
            rawCountsHist->SetBinContent(i+1, gFits.GetFitValues()[AliResonanceFits::kSig]);
            rawCountsHist->SetBinError(i+1, gFits.GetFitValues()[AliResonanceFits::kSigErr]);
            sOverBHist->SetBinContent(i+1,  gFits.GetFitValues()[AliResonanceFits::kSoverB]);
            sOverBHist->SetBinError(i+1,  gFits.GetFitValues()[AliResonanceFits::kSoverBerr]);
            signifHist->SetBinContent(i+1,  gFits.GetFitValues()[AliResonanceFits::kSignif]);
            chi2Hist->SetBinContent(i+1,  gFits.GetFitValues()[AliResonanceFits::kChisqSideBands]);
            chi2MCHist->SetBinContent(i+1,  gFits.GetFitValues()[AliResonanceFits::kChisqMCTotal]);
            effHist->SetBinContent(i+1, eff[i][isig]);
            effHist->SetBinError(i+1, effErr[i][isig]);
            sigFracHist->SetBinContent(i+1, sigFrac[i][isig]);
         }
         if(runSigExtrSyst) {
            rawCountsHist2D->SetBinContent(i+1, isig+1, gFits.GetFitValues()[AliResonanceFits::kSig]);
            rawCountsHist2D->SetBinError(i+1, isig+1, gFits.GetFitValues()[AliResonanceFits::kSigErr]);
            sOverBHist2D->SetBinContent(i+1, isig+1, gFits.GetFitValues()[AliResonanceFits::kSoverB]);
            sOverBHist2D->SetBinError(i+1, isig+1, gFits.GetFitValues()[AliResonanceFits::kSoverBerr]);
            signifHist2D->SetBinContent(i+1, isig+1, gFits.GetFitValues()[AliResonanceFits::kSignif]);
            chi2Hist2D->SetBinContent(i+1, isig+1, gFits.GetFitValues()[AliResonanceFits::kChisqSideBands]);
            chi2MCHist2D->SetBinContent(i+1, isig+1, gFits.GetFitValues()[AliResonanceFits::kChisqMCTotal]);
            effHist2D->SetBinContent(i+1, isig+1, eff[i][isig]);
            effHist2D->SetBinError(i+1, isig+1, effErr[i][isig]);
            sigFracHist2D->SetBinContent(i+1, isig+1, sigFrac[i][isig]);
         }
      }                                                     // end loop over signal counting windows
      if(runSigExtrSyst) {
         for(Int_t ifit = 0;ifit<gkNMassFitVariations;++ifit) {
            for(Int_t iexcl = 0;iexcl<gkNMassExclVariations;++iexcl) {
               if (ifit == 0 && iexcl == 0) continue;
               if (sigExtrMethod == 1) ConfigureSignalExtractionHybrid(ptMin, ptMax, gmFitLims[ifit][0], gmFitLims[ifit][1]);
               if (sigExtrMethod == 2) ConfigureSignalExtractionFit(ptMin, ptMax, gmFitLims[ifit][0], gmFitLims[ifit][1]);
               if (sigExtrMethod == 3) ConfigureSignalExtractionME(ptMin, ptMax, gmFitLims[ifit][0], gmFitLims[ifit][1], gmExclLims[iexcl][0], gmExclLims[iexcl][1]);
               if (sigExtrMethod == 4) ConfigureSignalExtractionLS(ptMin, ptMax, gmFitLims[ifit][0], gmFitLims[ifit][1], gmExclLims[iexcl][0], gmExclLims[iexcl][1]);
               
               gFits.Process(); 
               gFits.ComputeOutputValues(gmSigLims[0][0], gmSigLims[0][1]);
               // add QA histograms to the output list
               //AddSigExtrQAtoOutputList(&qaHistsForSigExtraction, settingNames[i].Data(), gmFitLims[ifit][0], gmFitLims[ifit][1], gmExclLims[iexcl][0], gmExclLims[iexcl][1]);
               
               Int_t ybin = gkNSigCountingVariations + ifit*gkNMassExclVariations+iexcl;
               rawCountsHist2D->SetBinContent(i+1, ybin, gFits.GetFitValues()[AliResonanceFits::kSig]);
               rawCountsHist2D->SetBinError(i+1, ybin, gFits.GetFitValues()[AliResonanceFits::kSigErr]);
               sOverBHist2D->SetBinContent(i+1, ybin, gFits.GetFitValues()[AliResonanceFits::kSoverB]);
               sOverBHist2D->SetBinError(i+1, ybin, gFits.GetFitValues()[AliResonanceFits::kSoverBerr]);
               signifHist2D->SetBinContent(i+1, ybin, gFits.GetFitValues()[AliResonanceFits::kSignif]);
               chi2Hist2D->SetBinContent(i+1, ybin, gFits.GetFitValues()[AliResonanceFits::kChisqSideBands]);
               chi2MCHist2D->SetBinContent(i+1, ybin, gFits.GetFitValues()[AliResonanceFits::kChisqMCTotal]);
               effHist2D->SetBinContent(i+1, ybin, eff[i][0]);
               effHist2D->SetBinError(i+1, ybin, effErr[i][0]);
               sigFracHist2D->SetBinContent(i+1, ybin, sigFrac[i][0]);
            }                                                  // end loop over exclusion range variations
         }                                                     // end loop over fit range variations
      }
   }                                                        // end loop over settings
   cout << "Finished settings " << endl;
      
   //for(Int_t i = 4; i <= nSettings; ++i) activeSettings->SetBinContent(i, 0.0);
   
   corrYieldHist->Divide(rawCountsHist, effHist);
   corrYieldHist->GetXaxis()->SetRange(1, nSettings);
   Double_t min = corrYieldHist->GetMinimum(); 
   Double_t max = corrYieldHist->GetMaximum();
   corrYieldHist->GetXaxis()->SetRange(1, 64);
   TH1D* corrYieldDistrib = new TH1D("CorrYieldDistrib", "Distribution of corrected yields", 50, min-2*(max-min), max+2*(max-min));
   corrYieldDistrib->GetXaxis()->SetTitle("corrected yield");
   corrYieldDistrib->GetYaxis()->SetTitle("counts");
   TH1D* corrYieldDistrib_weighted = new TH1D("CorrYieldDistrib_weighted", "Distribution of corrected yields, weighted", 50, min-2*(max-min), max+2*(max-min));
   corrYieldDistrib_weighted->GetXaxis()->SetTitle("corrected yield");
   corrYieldDistrib_weighted->GetYaxis()->SetTitle("counts");
   for(Int_t i = 1;i <= nSettings;++i) {
      if(activeSettings->GetBinContent(i) < 0.01) continue;
      corrYieldDistrib->Fill(corrYieldHist->GetBinContent(i));
   }
   outputHistograms.Add(corrYieldDistrib);
   outputHistograms.Add(corrYieldDistrib_weighted);
   
   
   crossSectionHist->Divide(rawCountsHist, effHist);
   if (ptMin < 0.0 || ptMax < 0.0) 
      crossSectionHist->Scale(1.0/totalLumi/gkBranchingRatio/gkDeltaY);
   else
      crossSectionHist->Scale(1.0/totalLumi/gkBranchingRatio/gkDeltaY/(ptMax-ptMin));
      
   crossSectionHist->GetXaxis()->SetRange(1, nSettings);
   min = crossSectionHist->GetMinimum(); 
   max = crossSectionHist->GetMaximum();
   crossSectionHist->GetXaxis()->SetRange(1, 64);
   TH1D* crossSectionDistrib = new TH1D("CrossSectionDistrib", "Distribution of cross-section results", 50, min-2*(max-min), max+2*(max-min));
   TH1D* crossSectionDistrib_weighted = new TH1D("CrossSectionDistrib_weighted", "Distribution of cross-section results,  weighted", 50, min-2*(max-min), max+2*(max-min));
   if (ptMin<0 || ptMax<0) {
      crossSectionDistrib->GetXaxis()->SetTitle("d#sigma / dy (nb)");
      crossSectionDistrib_weighted->GetXaxis()->SetTitle("d#sigma / dy (nb)");
   }
   else {
      crossSectionDistrib->GetXaxis()->SetTitle("d^{2}#sigma / dy / dp_{T} (nb)");
      crossSectionDistrib_weighted->GetXaxis()->SetTitle("d^{2}#sigma / dy / dp_{T} (nb)");
   }
   crossSectionDistrib->GetYaxis()->SetTitle("counts");
   crossSectionDistrib_weighted->GetYaxis()->SetTitle("counts");
   for (Int_t i = 1;i <= nSettings;++i) { 
      if(activeSettings->GetBinContent(i) < 0.01) continue;
      crossSectionDistrib->Fill(crossSectionHist->GetBinContent(i));
   }
   outputHistograms.Add(crossSectionDistrib);
   outputHistograms.Add(crossSectionDistrib_weighted);
   
   Double_t minStatErr = 0.0; Double_t maxStatErr = 0.0;
   for(Int_t i = 1; i <= nSettings; i++) {
      if (i == 1) {
         minStatErr = crossSectionHist->GetBinError(i);
         maxStatErr = crossSectionHist->GetBinError(i);
      }
      if(activeSettings->GetBinContent(i) < 0.01) continue;
      if (minStatErr>crossSectionHist->GetBinError(i)) minStatErr = crossSectionHist->GetBinError(i);
      if (maxStatErr<crossSectionHist->GetBinError(i)) maxStatErr = crossSectionHist->GetBinError(i);
   }
   TH1D* crossSectionStatErrDistrib = new TH1D("CrossSectionStatErrDistrib", "Distribution of cross-section stat err results", 50, minStatErr-2*(maxStatErr-minStatErr), maxStatErr+2*(maxStatErr-minStatErr));
   TH1D* crossSectionStatErrDistrib_weighted = new TH1D("CrossSectionStatErrDistrib_weighted", "Distribution of cross-section stat err results,  weighted", 50, minStatErr-2*(maxStatErr-minStatErr), maxStatErr+2*(maxStatErr-minStatErr));
   if (ptMin<0 || ptMax<0) {
      crossSectionStatErrDistrib->GetXaxis()->SetTitle("d#sigma / dy (nb)");
      crossSectionStatErrDistrib_weighted->GetXaxis()->SetTitle("d#sigma / dy (nb)");
   }
   else {
      crossSectionStatErrDistrib->GetXaxis()->SetTitle("d^{2}#sigma / dy / dp_{T} (nb)");
      crossSectionStatErrDistrib_weighted->GetXaxis()->SetTitle("d^{2}#sigma / dy / dp_{T} (nb)");
   }
   crossSectionStatErrDistrib->GetYaxis()->SetTitle("counts");
   crossSectionStatErrDistrib_weighted->GetYaxis()->SetTitle("counts");
   
   for (Int_t i = 1;i <= nSettings;++i)  {
      if(activeSettings->GetBinContent(i) < 0.01) continue;
      crossSectionStatErrDistrib->Fill(crossSectionHist->GetBinError(i));
   }
   outputHistograms.Add(crossSectionStatErrDistrib);
   outputHistograms.Add(crossSectionStatErrDistrib_weighted);
   
   
   if (runSigExtrSyst) {
      corrYieldHist2D->Divide(rawCountsHist2D, effHist2D);
      crossSectionHist2D->Divide(rawCountsHist2D, effHist2D);
      crossSectionHist2D->Scale(1.0/totalLumi/gkBranchingRatio/gkDeltaY);
      if (ptMin >= 0.0 && ptMax >= 0.0) 
         crossSectionHist2D->Scale(1.0/(ptMax-ptMin));
   }
   
   
   for(Int_t i = 1; i <= nSettings; i++) {
      if(!(corrYieldHist->GetBinContent(i)>0.0)) continue;
      deltaCorrected->SetBinContent(i, corrYieldHist->GetBinContent(i)-corrYieldHist->GetBinContent(1));
      deltaCorrected->SetBinError(i, TMath::Sqrt(TMath::Abs(corrYieldHist->GetBinError(i)*corrYieldHist->GetBinError(i) - corrYieldHist->GetBinError(1)*corrYieldHist->GetBinError(1))));
      if (corrYieldHist->GetBinContent(i)>0.0)
         deltaCorrectedFrac->SetBinContent(i, 100.0*deltaCorrected->GetBinContent(i)/corrYieldHist->GetBinContent(i));
      if (corrYieldHist->GetBinContent(i)>0.0)
         deltaCorrectedFrac->SetBinError(i, 100.0*deltaCorrected->GetBinError(i)/corrYieldHist->GetBinContent(i));
      
      cout << "setting " << i << endl;
      Double_t weight = 1.0; 
      if(i != 1) weight = (TMath::Abs(deltaCorrectedFrac->GetBinError(i)) == 0.0 ? 0.0 : TMath::Abs(deltaCorrectedFrac->GetBinContent(i) / deltaCorrectedFrac->GetBinError(i)));
      cout << "setting " << i << " after" << endl;
      corrYieldDistrib_weighted->Fill(corrYieldHist->GetBinContent(i), weight);
      crossSectionDistrib_weighted->Fill(crossSectionHist->GetBinContent(i), weight);
      crossSectionStatErrDistrib_weighted->Fill(crossSectionHist->GetBinError(i), weight);
      
      if(runSigExtrSyst) {
         for(Int_t j = 1; j<gkNSigCountingVariations+gkNMassFitVariations*gkNMassExclVariations; ++j) {
            deltaCorrected2D->SetBinContent(i,j, corrYieldHist2D->GetBinContent(i, j)-corrYieldHist->GetBinContent(i, 1));
            deltaCorrected2D->SetBinError(i,j, TMath::Sqrt(TMath::Abs(corrYieldHist2D->GetBinError(i,j)*corrYieldHist2D->GetBinError(i,j)-
                                                                      corrYieldHist2D->GetBinError(i,1)*corrYieldHist2D->GetBinError(i,1))));
            if (corrYieldHist2D->GetBinContent(i, 1)>0.0) {
               deltaCorrectedFrac2D->SetBinContent(i,j, 100.0*deltaCorrected2D->GetBinContent(i, j)/corrYieldHist2D->GetBinContent(i, 1));
               deltaCorrectedFrac2D->SetBinError(i,j, 100.0*deltaCorrected2D->GetBinError(i, j)/corrYieldHist->GetBinContent(i, 1));
            }
         }
      }
   }
   
   
   // do some postprocessing -> compute central values,  RMS's etc
   TH2D* corrYieldMsigDistrib2D = 0x0; TH2D* corrYieldMfitDistrib2D = 0x0;
   TH1D* corrYieldMsigDistrib2D_meanRMS = 0x0; TH1D* corrYieldMfitDistrib2D_meanRMS = 0x0;
   TH1D* corrYieldMsigDistrib2D_RMSdistrib = 0x0; TH1D* corrYieldMfitDistrib2D_RMSdistrib = 0x0;
   TH2D* crossSectionMsigDistrib2D = 0x0; TH2D* crossSectionMfitDistrib2D = 0x0;
   TH1D* crossSectionMsigDistrib2D_meanRMS = 0x0; TH1D* crossSectionMfitDistrib2D_meanRMS = 0x0;
   TH1D* crossSectionMsigDistrib2D_RMSdistrib = 0x0; TH1D* crossSectionMfitDistrib2D_RMSdistrib = 0x0;
   
   TH2D* corrYieldMsigDistrib2D_weighted = 0x0; TH2D* corrYieldMfitDistrib2D_weighted = 0x0;
   TH1D* corrYieldMsigDistrib2D_meanRMS_weighted = 0x0; TH1D* corrYieldMfitDistrib2D_meanRMS_weighted = 0x0;
   TH1D* corrYieldMsigDistrib2D_RMSdistrib_weighted = 0x0; TH1D* corrYieldMfitDistrib2D_RMSdistrib_weighted = 0x0;
   TH2D* crossSectionMsigDistrib2D_weighted = 0x0; TH2D* crossSectionMfitDistrib2D_weighted = 0x0;
   TH1D* crossSectionMsigDistrib2D_meanRMS_weighted = 0x0; TH1D* crossSectionMfitDistrib2D_meanRMS_weighted = 0x0;
   TH1D* crossSectionMsigDistrib2D_RMSdistrib_weighted = 0x0; TH1D* crossSectionMfitDistrib2D_RMSdistrib_weighted = 0x0;
   
   if (runSigExtrSyst) {
      corrYieldHist2D->GetXaxis()->SetRange(1, nSettings);
      corrYieldHist2D->GetYaxis()->SetRange(1, gkNSigCountingVariations);
      min = corrYieldHist2D->GetMinimum(); max = corrYieldHist2D->GetMaximum();
      corrYieldMsigDistrib2D = new TH2D("CorrYieldDistrib2D_msig", "Corrected yields distributions for sig counting variations", nSettings, 0., nSettings, 50, min-2*(max-min), max+2*(max-min));
      corrYieldMsigDistrib2D->GetYaxis()->SetTitle("corrected yield");
      corrYieldMsigDistrib2D_meanRMS = new TH1D("CorrYieldMsigDistrib2D_meanRMS", "Corrected yields mean and RMS of sig counting window variations per cut setting", nSettings, 0., nSettings);
      corrYieldMsigDistrib2D_meanRMS->GetYaxis()->SetTitle("yield");
      corrYieldMsigDistrib2D_weighted = new TH2D("CorrYieldDistrib2D_msig_weighted", "Corrected yields distributions for sig counting variations,  weighted", nSettings, 0., nSettings, 50, min-2*(max-min), max+2*(max-min));
      corrYieldMsigDistrib2D_weighted->GetYaxis()->SetTitle("corrected yield");
      corrYieldMsigDistrib2D_meanRMS_weighted = new TH1D("CorrYieldMsigDistrib2D_meanRMS_weighted", "Corrected yields mean and RMS of sig counting window variations per cut setting, weighted", nSettings, 0., nSettings);
      corrYieldMsigDistrib2D_meanRMS_weighted->GetYaxis()->SetTitle("yield");
      outputHistograms.Add(corrYieldMsigDistrib2D);
      outputHistograms.Add(corrYieldMsigDistrib2D_meanRMS);
      outputHistograms.Add(corrYieldMsigDistrib2D_weighted);
      outputHistograms.Add(corrYieldMsigDistrib2D_meanRMS_weighted);
      TH1D* temp1 = new TH1D("temp1", "", 50, min-2*(max-min), max+2*(max-min));
      TH1D* temp1_weighted = new TH1D("temp1_weighted", "", 50, min-2*(max-min), max+2*(max-min));
      
      corrYieldHist2D->GetYaxis()->SetRange(gkNSigCountingVariations+1, nTotalSigExtrVariations);
      min = corrYieldHist2D->GetMinimum(); max = corrYieldHist2D->GetMaximum();
      corrYieldMfitDistrib2D = new TH2D("CorrYieldDistrib2D_mfit", "Corrected yields distributions for fit ranges variations", nSettings, 0., nSettings, 50, min-2*(max-min), max+2*(max-min));
      corrYieldMfitDistrib2D->GetYaxis()->SetTitle("corrected yield");
      corrYieldMfitDistrib2D_meanRMS = new TH1D("CorrYieldMfitDistrib2D_meanRMS", "Corrected yields mean and RMS of fit variations per cut setting", nSettings, 0., nSettings);
      corrYieldMfitDistrib2D_meanRMS->GetYaxis()->SetTitle("yield");
      corrYieldMfitDistrib2D_weighted = new TH2D("CorrYieldDistrib2D_mfit_weighted", "Corrected yields distributions for fit ranges variations, weighted", nSettings, 0., nSettings, 50, min-2*(max-min), max+2*(max-min));
      corrYieldMfitDistrib2D_weighted->GetYaxis()->SetTitle("corrected yield");
      corrYieldMfitDistrib2D_meanRMS_weighted = new TH1D("CorrYieldMfitDistrib2D_meanRMS_weighted", "Corrected yields mean and RMS of fit variations per cut setting, weighted", nSettings, 0., nSettings);
      corrYieldMfitDistrib2D_meanRMS_weighted->GetYaxis()->SetTitle("yield");
      outputHistograms.Add(corrYieldMfitDistrib2D);
      outputHistograms.Add(corrYieldMfitDistrib2D_meanRMS);
      outputHistograms.Add(corrYieldMfitDistrib2D_weighted);
      outputHistograms.Add(corrYieldMfitDistrib2D_meanRMS_weighted);
      TH1D* temp2 = new TH1D("temp2", "", 50, min-2*(max-min), max+2*(max-min));
      TH1D* temp2_weighted = new TH1D("temp2_weighted", "", 50, min-2*(max-min), max+2*(max-min));
      corrYieldHist2D->GetXaxis()->SetRange(1, nSettings);
      corrYieldHist2D->GetYaxis()->SetRange(1, nTotalSigExtrVariations);
      
      crossSectionHist2D->GetXaxis()->SetRange(1, nSettings);
      crossSectionHist2D->GetYaxis()->SetRange(1, gkNSigCountingVariations);
      min = crossSectionHist2D->GetMinimum(); max = crossSectionHist2D->GetMaximum();
      crossSectionMsigDistrib2D = new TH2D("CrossSectionDistrib2D_msig", "Cross section distributions for sig counting variations", nSettings, 0., nSettings, 50, min-2*(max-min), max+2*(max-min));
      crossSectionMsigDistrib2D_weighted = new TH2D("CrossSectionDistrib2D_msig_weighted", "Cross section distributions for sig counting variations, weighted", nSettings, 0., nSettings, 50, min-2*(max-min), max+2*(max-min));
      if (ptMin<0 || ptMax<0) {
         crossSectionMsigDistrib2D->GetYaxis()->SetTitle("d#sigma / dy (nb)");
         crossSectionMsigDistrib2D_weighted->GetYaxis()->SetTitle("d#sigma / dy (nb)");
      }
      else {
         crossSectionMsigDistrib2D->GetXaxis()->SetTitle("d^{2}#sigma / dy / dp_{T} (nb)");
         crossSectionMsigDistrib2D_weighted->GetXaxis()->SetTitle("d^{2}#sigma / dy / dp_{T} (nb)");
      }
      crossSectionMsigDistrib2D_meanRMS = new TH1D("CrossSectionMsigDistrib2D_meanRMS", "Cross-section mean and RMS of sig counting window variations per cut setting", nSettings, 0., nSettings);
      crossSectionMsigDistrib2D_meanRMS_weighted = new TH1D("CrossSectionMsigDistrib2D_meanRMS_weighted", "Cross-section mean and RMS of sig counting window variations per cut setting, weighted", nSettings, 0., nSettings);
      if (ptMin<0 || ptMax<0) {
         crossSectionMsigDistrib2D_meanRMS->GetYaxis()->SetTitle("d#sigma / dy (nb)");
         crossSectionMsigDistrib2D_meanRMS_weighted->GetYaxis()->SetTitle("d#sigma / dy (nb)");
      }
      else {
         crossSectionMsigDistrib2D_meanRMS->GetYaxis()->SetTitle("d^{2}#sigma / dy / dp_{T} (nb)");
         crossSectionMsigDistrib2D_meanRMS_weighted->GetYaxis()->SetTitle("d^{2}#sigma / dy / dp_{T} (nb)");
      }
      outputHistograms.Add(crossSectionMsigDistrib2D);
      outputHistograms.Add(crossSectionMsigDistrib2D_meanRMS);
      outputHistograms.Add(crossSectionMsigDistrib2D_weighted);
      outputHistograms.Add(crossSectionMsigDistrib2D_meanRMS_weighted);
      TH1D* temp3 = new TH1D("temp3", "", 50, min-2*(max-min), max+2*(max-min));
      TH1D* temp3_weighted = new TH1D("temp3_weighted", "", 50, min-2*(max-min), max+2*(max-min));
      
      crossSectionHist2D->GetYaxis()->SetRange(gkNSigCountingVariations+1, nTotalSigExtrVariations);
      min = crossSectionHist2D->GetMinimum(); max = crossSectionHist2D->GetMaximum();
      crossSectionMfitDistrib2D = new TH2D("CrossSectionDistrib2D_mfit", "Cross section distributions for fit ranges variations", nSettings, 0., nSettings, 50, min-2*(max-min), max+2*(max-min));
      crossSectionMfitDistrib2D_weighted = new TH2D("CrossSectionDistrib2D_mfit_weighted", "Cross section distributions for fit ranges variations, weighted", nSettings, 0., nSettings, 50, min-2*(max-min), max+2*(max-min));
      if (ptMin<0 || ptMax<0) {
         crossSectionMfitDistrib2D->GetYaxis()->SetTitle("d#sigma / dy (nb)");
         crossSectionMfitDistrib2D_weighted->GetYaxis()->SetTitle("d#sigma / dy (nb)");
      }
      else {
         crossSectionMfitDistrib2D->GetXaxis()->SetTitle("d^{2}#sigma / dy / dp_{T} (nb)");
         crossSectionMfitDistrib2D_weighted->GetXaxis()->SetTitle("d^{2}#sigma / dy / dp_{T} (nb)");
      }
      crossSectionMfitDistrib2D_meanRMS = new TH1D("CrossSectionMfitDistrib2D_meanRMS", "Cross-section mean and RMS of fit range variations per cut setting", nSettings, 0., nSettings);
      crossSectionMfitDistrib2D_meanRMS_weighted = new TH1D("CrossSectionMfitDistrib2D_meanRMS_weighted", "Cross-section mean and RMS of fit range variations per cut setting, weighted", nSettings, 0., nSettings);
      if (ptMin<0 || ptMax<0) {
         crossSectionMfitDistrib2D_meanRMS->GetYaxis()->SetTitle("d#sigma / dy (nb)");
         crossSectionMfitDistrib2D_meanRMS_weighted->GetYaxis()->SetTitle("d#sigma / dy (nb)");
      }
      else {
         crossSectionMfitDistrib2D_meanRMS->GetYaxis()->SetTitle("d^{2}#sigma / dy / dp_{T} (nb)");
         crossSectionMfitDistrib2D_meanRMS_weighted->GetYaxis()->SetTitle("d^{2}#sigma / dy / dp_{T} (nb)");
      }
      outputHistograms.Add(crossSectionMfitDistrib2D);
      outputHistograms.Add(crossSectionMfitDistrib2D_meanRMS);
      outputHistograms.Add(crossSectionMfitDistrib2D_weighted);
      outputHistograms.Add(crossSectionMfitDistrib2D_meanRMS_weighted);
      crossSectionHist2D->GetXaxis()->SetRange(1, nSettings);
      crossSectionHist2D->GetYaxis()->SetRange(1, nTotalSigExtrVariations);
      TH1D* temp4 = new TH1D("temp4", "", 50, min-2*(max-min), max+2*(max-min));
      TH1D* temp4_weighted = new TH1D("temp4_weighted", "", 50, min-2*(max-min), max+2*(max-min));
      
      Double_t minRMS1 = 0.0; Double_t maxRMS1 = 0.0;
      Double_t minRMS2 = 0.0; Double_t maxRMS2 = 0.0;
      Double_t minRMS3 = 0.0; Double_t maxRMS3 = 0.0;
      Double_t minRMS4 = 0.0; Double_t maxRMS4 = 0.0;
      for(Int_t i = 1; i <= nSettings; i++) {
         corrYieldMsigDistrib2D->GetXaxis()->SetBinLabel(i, settingNames[i-1].Data());
         corrYieldMsigDistrib2D_meanRMS->GetXaxis()->SetBinLabel(i, settingNames[i-1].Data());
         corrYieldMfitDistrib2D->GetXaxis()->SetBinLabel(i, settingNames[i-1].Data());
         corrYieldMfitDistrib2D_meanRMS->GetXaxis()->SetBinLabel(i, settingNames[i-1].Data());
         crossSectionMsigDistrib2D->GetXaxis()->SetBinLabel(i, settingNames[i-1].Data());
         crossSectionMsigDistrib2D_meanRMS->GetXaxis()->SetBinLabel(i, settingNames[i-1].Data());
         crossSectionMfitDistrib2D->GetXaxis()->SetBinLabel(i, settingNames[i-1].Data());
         crossSectionMfitDistrib2D_meanRMS->GetXaxis()->SetBinLabel(i, settingNames[i-1].Data());
         
         temp1->Reset(); temp2->Reset(); temp3->Reset(); temp4->Reset();
         temp1_weighted->Reset(); temp2_weighted->Reset(); temp3_weighted->Reset(); temp4_weighted->Reset();
         for (Int_t j = 1;j <= gkNSigCountingVariations;++j) {
            cout << "setting " << i << " sig " << j << endl;
            Double_t weight = 1.0; 
            if(j != 1) weight = (TMath::Abs(deltaCorrectedFrac2D->GetBinError(i, j)) == 0.0 ? 0.0 : TMath::Abs(deltaCorrectedFrac2D->GetBinContent(i, j) / deltaCorrectedFrac2D->GetBinError(i, j)));
            cout << "setting " << i << " sig " << j << "  after" << endl;
            corrYieldMsigDistrib2D->Fill(i-1, corrYieldHist2D->GetBinContent(i, j));
            corrYieldMsigDistrib2D_weighted->Fill(i-1, corrYieldHist2D->GetBinContent(i, j), weight);
            temp1->Fill(corrYieldHist2D->GetBinContent(i, j));
            temp1_weighted->Fill(corrYieldHist2D->GetBinContent(i, j), weight);
            crossSectionMsigDistrib2D->Fill(i-1, crossSectionHist2D->GetBinContent(i, j));
            crossSectionMsigDistrib2D_weighted->Fill(i-1, crossSectionHist2D->GetBinContent(i, j), weight);
            temp3->Fill(crossSectionHist2D->GetBinContent(i, j));
            temp3_weighted->Fill(crossSectionHist2D->GetBinContent(i, j), weight);
            if (j == 1) {
               corrYieldMfitDistrib2D->Fill(i-1, corrYieldHist2D->GetBinContent(i, j));
               temp2->Fill(corrYieldHist2D->GetBinContent(i, j));
               crossSectionMfitDistrib2D->Fill(i-1, crossSectionHist2D->GetBinContent(i, j));
               temp4->Fill(crossSectionHist2D->GetBinContent(i, j));
            }
         }
         for (Int_t j = gkNSigCountingVariations+1;j <= nTotalSigExtrVariations;++j) {
            cout << "setting " << i << " sig " << j << endl;
            Double_t weight = 1.0; 
            if(j != 1) weight = (TMath::Abs(deltaCorrectedFrac2D->GetBinError(i,j)) == 0.0 ? 0.0 : TMath::Abs(deltaCorrectedFrac2D->GetBinContent(i, j) / deltaCorrectedFrac2D->GetBinError(i, j)));
            cout << "setting " << i << " sig " << j << "  after" << endl;
            corrYieldMfitDistrib2D->Fill(i-1, corrYieldHist2D->GetBinContent(i, j));
            corrYieldMfitDistrib2D_weighted->Fill(i-1, corrYieldHist2D->GetBinContent(i, j), weight);
            temp2->Fill(corrYieldHist2D->GetBinContent(i, j));
            temp2_weighted->Fill(corrYieldHist2D->GetBinContent(i, j), weight);
            crossSectionMfitDistrib2D->Fill(i-1, crossSectionHist2D->GetBinContent(i, j));
            crossSectionMfitDistrib2D_weighted->Fill(i-1, crossSectionHist2D->GetBinContent(i, j), weight);
            temp4->Fill(crossSectionHist2D->GetBinContent(i, j));
            temp4_weighted->Fill(crossSectionHist2D->GetBinContent(i, j), weight);
         }
         corrYieldMsigDistrib2D_meanRMS->SetBinContent(i, temp1->GetMean());
         corrYieldMsigDistrib2D_meanRMS->SetBinError(i, temp1->GetRMS());
         corrYieldMfitDistrib2D_meanRMS->SetBinContent(i, temp2->GetMean());
         corrYieldMfitDistrib2D_meanRMS->SetBinError(i, temp2->GetRMS());
         crossSectionMsigDistrib2D_meanRMS->SetBinContent(i, temp3->GetMean());
         crossSectionMsigDistrib2D_meanRMS->SetBinError(i, temp3->GetRMS());
         crossSectionMfitDistrib2D_meanRMS->SetBinContent(i, temp4->GetMean());
         crossSectionMfitDistrib2D_meanRMS->SetBinError(i, temp4->GetRMS());
         corrYieldMsigDistrib2D_meanRMS_weighted->SetBinContent(i, temp1_weighted->GetMean());
         corrYieldMsigDistrib2D_meanRMS_weighted->SetBinError(i, temp1_weighted->GetRMS());
         corrYieldMfitDistrib2D_meanRMS_weighted->SetBinContent(i, temp2_weighted->GetMean());
         corrYieldMfitDistrib2D_meanRMS_weighted->SetBinError(i, temp2_weighted->GetRMS());
         crossSectionMsigDistrib2D_meanRMS_weighted->SetBinContent(i, temp3_weighted->GetMean());
         crossSectionMsigDistrib2D_meanRMS_weighted->SetBinError(i, temp3_weighted->GetRMS());
         crossSectionMfitDistrib2D_meanRMS_weighted->SetBinContent(i, temp4_weighted->GetMean());
         crossSectionMfitDistrib2D_meanRMS_weighted->SetBinError(i, temp4_weighted->GetRMS());
         if (i == 1) {
            minRMS1 = temp1->GetRMS(); maxRMS1 = temp1->GetRMS();
            minRMS2 = temp2->GetRMS(); maxRMS2 = temp2->GetRMS();
            minRMS3 = temp3->GetRMS(); maxRMS3 = temp3->GetRMS();
            minRMS4 = temp4->GetRMS(); maxRMS4 = temp4->GetRMS();
         }
         else{
            if (minRMS1>temp1->GetRMS()) minRMS1 = temp1->GetRMS();
            if (maxRMS1<temp1->GetRMS()) maxRMS1 = temp1->GetRMS();
            if (minRMS2>temp2->GetRMS()) minRMS2 = temp2->GetRMS();
            if (maxRMS2<temp2->GetRMS()) maxRMS2 = temp2->GetRMS();
            if (minRMS3>temp3->GetRMS()) minRMS3 = temp3->GetRMS();
            if (maxRMS3<temp3->GetRMS()) maxRMS3 = temp3->GetRMS();
            if (minRMS4>temp4->GetRMS()) minRMS4 = temp4->GetRMS();
            if (maxRMS4<temp4->GetRMS()) maxRMS4 = temp4->GetRMS();
         }
      }                                                     // end loop over settings
      corrYieldMsigDistrib2D_RMSdistrib = new TH1D("CorrYieldMsigDistrib2D_RMSdistrib", "Distribution of RMS's from signal counting window variations for all cut settings,  corrected yields", 
                  50, minRMS1-2*(maxRMS1-minRMS1), maxRMS1+2*(maxRMS1-minRMS1));
      corrYieldMsigDistrib2D_RMSdistrib_weighted = new TH1D("CorrYieldMsigDistrib2D_RMSdistrib_weighted", "Distribution of RMS's from signal counting window variations for all cut settings,  corrected yields, weighted", 
         50, minRMS1-2*(maxRMS1-minRMS1), maxRMS1+2*(maxRMS1-minRMS1));
      corrYieldMfitDistrib2D_RMSdistrib = new TH1D("CorrYieldMfitDistrib2D_RMSdistrib", "Distribution of RMS's from fit range variations for all cut settings,  corrected yields", 
                  50, minRMS2-2*(maxRMS2-minRMS2), maxRMS2+2*(maxRMS2-minRMS2));
      corrYieldMfitDistrib2D_RMSdistrib_weighted = new TH1D("CorrYieldMfitDistrib2D_RMSdistrib_weighted", "Distribution of RMS's from fit range variations for all cut settings,  corrected yields, weighted", 
         50, minRMS2-2*(maxRMS2-minRMS2), maxRMS2+2*(maxRMS2-minRMS2));
      crossSectionMsigDistrib2D_RMSdistrib = new TH1D("CrossSectionMsigDistrib2D_RMSdistrib", "Distribution of RMS's from signal counting window variations for all cut settings, cross sections", 
         50, minRMS3-2*(maxRMS3-minRMS3), maxRMS3+2*(maxRMS3-minRMS3));
      crossSectionMsigDistrib2D_RMSdistrib_weighted = new TH1D("CrossSectionMsigDistrib2D_RMSdistrib_weighted", "Distribution of RMS's from signal counting window variations for all cut settings, cross sections, weighted", 
         50, minRMS3-2*(maxRMS3-minRMS3), maxRMS3+2*(maxRMS3-minRMS3));
      crossSectionMfitDistrib2D_RMSdistrib = new TH1D("CrossSectionMfitDistrib2D_RMSdistrib", "Distribution of RMS's from fit range variations for all cut settings, cross sections", 
         50, minRMS4-2*(maxRMS4-minRMS4), maxRMS4+2*(maxRMS4-minRMS4));
      crossSectionMfitDistrib2D_RMSdistrib_weighted = new TH1D("CrossSectionMfitDistrib2D_RMSdistrib_weighted", "Distribution of RMS's from fit range variations for all cut settings, cross sections, weighted", 
         50, minRMS4-2*(maxRMS4-minRMS4), maxRMS4+2*(maxRMS4-minRMS4));
      outputHistograms.Add(corrYieldMsigDistrib2D_RMSdistrib);
      outputHistograms.Add(corrYieldMfitDistrib2D_RMSdistrib);
      outputHistograms.Add(crossSectionMsigDistrib2D_RMSdistrib);
      outputHistograms.Add(crossSectionMfitDistrib2D_RMSdistrib);
      outputHistograms.Add(corrYieldMsigDistrib2D_RMSdistrib_weighted);
      outputHistograms.Add(corrYieldMfitDistrib2D_RMSdistrib_weighted);
      outputHistograms.Add(crossSectionMsigDistrib2D_RMSdistrib_weighted);
      outputHistograms.Add(crossSectionMfitDistrib2D_RMSdistrib_weighted);
      for(Int_t i = 1; i <= nSettings; i++) {
         if(activeSettings->GetBinContent(i) < 0.01) continue;
         corrYieldMsigDistrib2D_RMSdistrib->Fill(corrYieldMsigDistrib2D_meanRMS->GetBinError(i));
         corrYieldMfitDistrib2D_RMSdistrib->Fill(corrYieldMfitDistrib2D_meanRMS->GetBinError(i));
         crossSectionMsigDistrib2D_RMSdistrib->Fill(crossSectionMsigDistrib2D_meanRMS->GetBinError(i));
         crossSectionMfitDistrib2D_RMSdistrib->Fill(crossSectionMfitDistrib2D_meanRMS->GetBinError(i));
         corrYieldMsigDistrib2D_RMSdistrib_weighted->Fill(corrYieldMsigDistrib2D_meanRMS_weighted->GetBinError(i));
         corrYieldMfitDistrib2D_RMSdistrib_weighted->Fill(corrYieldMfitDistrib2D_meanRMS_weighted->GetBinError(i));
         crossSectionMsigDistrib2D_RMSdistrib_weighted->Fill(crossSectionMsigDistrib2D_meanRMS_weighted->GetBinError(i));
         crossSectionMfitDistrib2D_RMSdistrib_weighted->Fill(crossSectionMfitDistrib2D_meanRMS_weighted->GetBinError(i));
      }
   }                                                        //  end if (runSigExtrSyst)
   
   // Write results and QA histograms
   TFile* save = new TFile(Form("%s/Results.root", outputDir), "RECREATE");
   outputHistograms.Write();
   save->mkdir("SignalExtractionQA");
   save->mkdir("Efficiency_vs_Run");
   save->mkdir("SignalFraction_vs_Run");
   save->mkdir("Efficiency_vs_Pt");
   save->mkdir("EfficiencyQA_PerRun");
                            
   for (Int_t iSetting=0; iSetting<nSettings; ++iSetting) {
      save->cd(); save->cd("SignalExtractionQA");
      TString strName = Form("%s_mfit%.2f_%.2f_mexcl%.2f_%.2f", settingNames[iSetting].Data(), gmFitLims[0][0], gmFitLims[0][1], gmExclLims[0][0], gmExclLims[0][1]);
      TObject* obj = qaHistsForSigExtraction.FindObject(Form("SplusB_%s", strName.Data()));
      if(obj) obj->Write();
      obj = qaHistsForSigExtraction.FindObject(Form("Bkg_%s", strName.Data()));
      if(obj) obj->Write();
      obj = qaHistsForSigExtraction.FindObject(Form("Signal_%s", strName.Data()));
      if(obj) obj->Write();
      obj = qaHistsForSigExtraction.FindObject(Form("SoverB_%s", strName.Data()));
      if(obj) obj->Write();
      obj = qaHistsForSigExtraction.FindObject(Form("SplusResidualBkg_%s", strName.Data()));
      if(obj) obj->Write();
      obj = qaHistsForSigExtraction.FindObject(Form("BkgCombinatorial_%s", strName.Data()));
      if(obj) obj->Write();
      obj = qaHistsForSigExtraction.FindObject(Form("BkgResidual_%s", strName.Data()));
      if(obj) obj->Write();
      obj = qaHistsForSigExtraction.FindObject(Form("SignalMC_%s", strName.Data()));
      if(obj) obj->Write();
      obj = qaHistsForSigExtraction.FindObject(Form("BkgFunc_%s", strName.Data()));
      if(obj) obj->Write();
      
      save->cd(); save->cd("Efficiency_vs_Run");
      obj = qaHistsForEfficiency->FindObject(Form("effVsRun_%s_%s", settingNames[iSetting].Data(), mcSignalElectron));
      if(obj) obj->Write();
      
      save->cd(); save->cd("SignalFraction_vs_Run");
      obj = qaHistsForEfficiency->FindObject(Form("sigFracVsRun_%s_%s", settingNames[iSetting].Data(), mcSignalElectron));
      if(obj) obj->Write();
      
      save->cd(); save->cd("Efficiency_vs_Pt");
      obj = qaHistsForEfficiency->FindObject(Form("effVsPtAverage_%s_%s", settingNames[iSetting].Data(), mcSignalElectron));
      if(obj) obj->Write();
      
      save->cd(); save->cd("EfficiencyQA_PerRun");
      save->mkdir(Form("EfficiencyQA_PerRun/%s_%s", settingNames[iSetting].Data(), mcSignalElectron));
      save->cd(Form("EfficiencyQA_PerRun/%s_%s", settingNames[iSetting].Data(), mcSignalElectron));
      for(Int_t iRun = 0; iRun<gkNRuns; ++iRun) {
         TList* listForRun = (TList*)qaHistsForEfficiency->FindObject(Form("Run_%d", gRunList[iRun]));
         if(!listForRun) continue;
         obj = listForRun->FindObject(Form("ProjMass_run%d_%s", gRunList[iRun], settingNames[iSetting].Data()));
         if (obj) obj->Write();
         obj = listForRun->FindObject(Form("ProjPt_run%d_%s_msig0", gRunList[iRun], settingNames[iSetting].Data()));
         if (obj) obj->Write();
         obj = listForRun->FindObject(Form("effVsPt_run%d_%s_msig0", gRunList[iRun], settingNames[iSetting].Data()));
         if (obj) obj->Write();
         obj = listForRun->FindObject(Form("PtMCdenominator_run%d_%s", gRunList[iRun], mcSignalElectron));
         if (obj) obj->Write();
      }
   }

   save->mkdir("Luminosity");
   save->cd("Luminosity");
   lumiHists->Write();
   cout <<  "bbb" << endl;
   save->Close();
} 

//____________________________________________________________________________________________________
void AddSigExtrQAtoOutputList(TList* outList, const Char_t* settingName, Double_t mFitLow, Double_t mFitHigh, Double_t mExclLow, Double_t mExclHigh) {
   // 
   // get signal extraction histograms from the current state of the gFits and add them to the output list 
   // 
   TH1* tempHist = gFits.GetSplusB();
   TH1* sPlusBhist = 0x0;
   if(tempHist) sPlusBhist = (TH1*)tempHist->Clone(Form("SplusB_%s_mfit%.2f_%.2f_mexcl%.2f_%.2f", settingName, mFitLow, mFitHigh, mExclLow, mExclHigh));
   if(sPlusBhist) outList->Add(sPlusBhist);

   TH1* bkgHist = 0x0;
   tempHist = gFits.GetBkg();
   if(tempHist) bkgHist = (TH1*)tempHist->Clone(Form("Bkg_%s_mfit%.2f_%.2f_mexcl%.2f_%.2f", settingName, mFitLow, mFitHigh, mExclLow, mExclHigh));
   if(bkgHist) outList->Add(bkgHist);

   TH1* signalHist = 0x0;
   tempHist = gFits.GetSignal();
   if(tempHist) signalHist = (TH1*)tempHist->Clone(Form("Signal_%s_mfit%.2f_%.2f_mexcl%.2f_%.2f", settingName, mFitLow, mFitHigh, mExclLow, mExclHigh));
   if(signalHist) outList->Add(signalHist);

   TH1* sOverBhist = 0x0;
   tempHist = gFits.GetSoverB();
   if(tempHist) sOverBhist = (TH1*)tempHist->Clone(Form("SoverB_%s_mfit%.2f_%.2f_mexcl%.2f_%.2f", settingName, mFitLow, mFitHigh, mExclLow, mExclHigh));
   if(sOverBhist) outList->Add(sOverBhist);     

   TH1* sPlusResidualBkgHist = 0x0;
   tempHist = gFits.GetSplusResidualBkg();
   if(tempHist) sPlusResidualBkgHist = (TH1*)tempHist->Clone(Form("SplusResidualBkg_%s_mfit%.2f_%.2f_mexcl%.2f_%.2f", settingName, mFitLow, mFitHigh, mExclLow, mExclHigh));
   if(sPlusResidualBkgHist) outList->Add(sPlusResidualBkgHist);

   TH1* bkgCombinatorialHist = 0x0; 
   tempHist = gFits.GetBkgCombinatorial();
   if(tempHist) bkgCombinatorialHist = (TH1*)tempHist->Clone(Form("BkgCombinatorial_%s_mfit%.2f_%.2f_mexcl%.2f_%.2f", settingName, mFitLow, mFitHigh, mExclLow, mExclHigh));
   if (bkgCombinatorialHist) outList->Add(bkgCombinatorialHist);

   TH1* residualBkgHist = 0x0;
   tempHist = gFits.GetResidualBkg();
   if(tempHist) residualBkgHist = (TH1*)tempHist->Clone(Form("BkgResidual_%s_mfit%.2f_%.2f_mexcl%.2f_%.2f", settingName, mFitLow, mFitHigh, mExclLow, mExclHigh));
   if(residualBkgHist) outList->Add(residualBkgHist);

   TH1* signalMChist = 0x0;
   tempHist = gFits.GetSignalMC();
   if(tempHist) signalMChist = (TH1*)tempHist->Clone(Form("SignalMC_%s_mfit%.2f_%.2f_mexcl%.2f_%.2f", settingName, mFitLow, mFitHigh, mExclLow, mExclHigh));
   if(signalMChist) outList->Add(signalMChist);
      
   TF1* bkgFunc = 0x0;
   TF1* tempFunc = gFits.GetBkgFitFunction();
   if (tempFunc) bkgFunc = (TF1*)tempFunc->Clone(Form("BkgFunc_%s_mfit%.2f_%.2f_mexcl%.2f_%.2f", settingName, mFitLow, mFitHigh, mExclLow, mExclHigh));
   if (bkgFunc) outList->Add(bkgFunc);
}

//____________________________________________________________________________________________________
TList* ComputeIntegratedLumi(const Char_t* dirLHC17pq = "", const Char_t* dirLHC15n = "", Int_t lumiOption = 6, Double_t& lumi, Double_t& lumiErr) {
   // 
   // compute integrated luminosity
   // 
   // options for lumiOption:
   // 0 - LHC15n
   // 1 - LHC17p
   // 2 - LHC17q
   // 3 - LHC17p + LHC17q
   // 4 - LHC17p + LHC17q (CENT_woSDD)
   // 5 - LHC17p + LHC17q (FAST)
   // 6 - LHC15n + LHC17p + LHC17q

   TH3F* zHist = (TH3F*)gHistManData->GetHistogram("EventTag_BeforeCuts", "EventTags_VtxZ_Run");
   
   TList* outputList = new TList();
   outputList->SetOwner(kTRUE);
   TH1D* nINT7evHist = new TH1D("nINT7evHist", "", zHist->GetXaxis()->GetNbins(), zHist->GetXaxis()->GetXbins()->GetArray());
   TH1D* zFracHist = new TH1D("zFracHist", "Fraction of events in -10 to +10 cm per run", zHist->GetXaxis()->GetNbins(), zHist->GetXaxis()->GetXbins()->GetArray());
   TH1D* chi2ZfitHist = new TH1D("chi2ZfitHist", "#chi^{2}/NDF of the Z vertex fit", zHist->GetXaxis()->GetNbins(), zHist->GetXaxis()->GetXbins()->GetArray());
   TH1D* meanfitHist = new TH1D("meanZfitHist", "Mean of the Z vertex distrib", zHist->GetXaxis()->GetNbins(), zHist->GetXaxis()->GetXbins()->GetArray());
   TH1D* sigmafitHist = new TH1D("sigmaZfitHist", "Width of the Z vertex distrib", zHist->GetXaxis()->GetNbins(), zHist->GetXaxis()->GetXbins()->GetArray());
   TH1D* lumiHist = new TH1D("lumiHist", "", zHist->GetXaxis()->GetNbins(), zHist->GetXaxis()->GetXbins()->GetArray());
   outputList->Add(nINT7evHist);
   outputList->Add(zFracHist);
   outputList->Add(lumiHist);
   outputList->Add(chi2ZfitHist);
   outputList->Add(meanfitHist);
   outputList->Add(sigmafitHist);
   
   Double_t nTotalINT7events = 0.0;
   Double_t nTotalINT7events2015 = 0.0;
   Double_t nTotalINT7events2017 = 0.0;
   lumi = 0.0;
   Double_t lumi2015 = 0.0;
   Double_t lumi2017 = 0.0;
   Double_t lumiErr = 0.0;
   for(Int_t irun = 1; irun <= zHist->GetXaxis()->GetNbins(); ++irun) {
      TString runlabel = zHist->GetXaxis()->GetBinLabel(irun);
      nINT7evHist->GetXaxis()->SetBinLabel(irun, runlabel.Data());
      zFracHist->GetXaxis()->SetBinLabel(irun, runlabel.Data());
      lumiHist->GetXaxis()->SetBinLabel(irun, runlabel.Data());
      chi2ZfitHist->GetXaxis()->SetBinLabel(irun, runlabel.Data());
      meanfitHist->GetXaxis()->SetBinLabel(irun, runlabel.Data());
      sigmafitHist->GetXaxis()->SetBinLabel(irun, runlabel.Data());
      
      Int_t run = runlabel.Atoi();
      Bool_t is2015 = (run < 280000  ? kTRUE : kFALSE);
      Bool_t isLHC17p = (run >= 282008 && run <= 282343  ? kTRUE : kFALSE);
      Bool_t isLHC17q = (run >= 282365 && run <= 282367  ? kTRUE : kFALSE);
      if(is2015 && !(lumiOption == 0 || lumiOption == 6)) continue;
      if(!is2015 && lumiOption == 0) continue;
      if(isLHC17p && (lumiOption == 0 || lumiOption == 2)) continue;
      if(isLHC17q && (lumiOption == 0 || lumiOption == 1)) continue;
      
      Double_t nev = 0.0;
      if (is2015) nev = GetINT7EventsForRun(dirLHC15n, run);
      else        nev = GetINT7EventsForRun(dirLHC17pq, run);
            
      Double_t frac = 0.0; Double_t chi2Zfit = 0.0;
      Double_t mean = 0.0; Double_t meanErr = 0.0;
      Double_t sigma = 0.0; Double_t sigmaErr = 0.0;
      Double_t fracErr = 0.0;
      frac = GetEventFractionInVertexRange(zHist, outputList, run, fracErr, chi2Zfit, mean, meanErr, sigma, sigmaErr);
      zFracHist->SetBinContent(irun, frac);
      zFracHist->SetBinError(irun, fracErr);
      chi2ZfitHist->SetBinContent(irun, chi2Zfit);
      meanfitHist->SetBinContent(irun, mean);
      meanfitHist->SetBinError(irun, meanErr);
      sigmafitHist->SetBinContent(irun, sigma);
      sigmafitHist->SetBinError(irun, sigmaErr);
      nTotalINT7events += nev*frac;
      nINT7evHist->SetBinContent(irun, nev*frac);
      if (is2015) nTotalINT7events2015 += nev*frac;
      else nTotalINT7events2017 += nev*frac;
      
      Double_t tempLumi = nev*frac / (is2015 ? gkV0ANDcrossSection2015 : gkV0ANDcrossSection2017);
      lumi += tempLumi;
      if (is2015) lumi2015 += tempLumi;
      else lumi2017 += tempLumi;
      lumiErr += TMath::Power(tempLumi*fracErr/frac, 2.0);
      lumiHist->SetBinContent(irun, tempLumi);
      lumiHist->SetBinError(irun, tempLumi*fracErr/frac);
   }
   lumiErr = TMath::Sqrt(lumiErr) / 1.0e+6;
   nINT7evHist->SetTitle(Form("Number of INT7 events for each run, total = %.4e events", nTotalINT7events));
   lumi /= 1.0e+6;                                          //  transform to nb^-1
   // Compute the integrated lumi uncertainty
   // NOTE: Assumptions: in 2017 data we currently use the same value of sigma_V0AND
   //                    Uncertaionty on the 2017 sigma_V0AND is 5% and considered fully correlated to the one from 2015
   Double_t uncertRatio = gkV0ANDcrossSection2017Err / gkV0ANDcrossSection2015Err;
   Double_t lumiErrBeam = (nTotalINT7events2015 + uncertRatio*nTotalINT7events2017) * 0.01 * gkV0ANDcrossSection2015Err / gkV0ANDcrossSection2015 / 1.0e+6;
   lumiHist->SetTitle(Form("Integrated luminosity for each run, total = %.3e #pm %.3e (stat) #pm %.3e (beam) nb^{-1}", lumi, lumiErr,  lumiErrBeam));
   return outputList;
}

//____________________________________________________________________________________________________
Int_t GetINT7EventsForRun(const Char_t* dir, Int_t run) {
   // 
   // Get the number of INT7 triggers for a given run number
   // 
   TFile* dstFile = TFile::Open(Form("%s/000%d/chunk1/dstTree.root", dir, run));
   TH2I* evStats = (TH2I*)dstFile->Get("EventStatistics");
   Bool_t is2015 = (run < 280000  ? kTRUE : kFALSE);
   Int_t triggerBin = evStats->GetYaxis()->FindBin("INT7");
   Int_t evCounterBin = (is2015 ? 6 : 2);                   // 2017 - "PS and trigger selected" (old event stats histogram)
   // 2015 - "TS and pileup checked (PC)" (new event stats histogram)
   Int_t nINT7events = evStats->GetBinContent(evCounterBin, triggerBin);   
   return nINT7events;
}

//________________________________________________________________________________________
Double_t GetEventFractionInVertexRange(TH3F* hist, TList* histList, Int_t run, Double_t& fracErr, Double_t& chi2, Double_t& mean, Double_t& meanErr, Double_t& sigma, Double_t& sigmaErr) {
   // 
   // find the vertex distribution for unbiased events in a given run and compute the fraction of them in -10<z<10 cm
   // 
   Int_t runBin = hist->GetXaxis()->FindBin(Form("%d", run));
   Int_t unbiasedEvBin = hist->GetYaxis()->FindBin("unbiased event");
   TH1D* vtxDistrib = hist->ProjectionZ(Form("vtxZdistrib%d", run), runBin, runBin, unbiasedEvBin, unbiasedEvBin);
   Bool_t is2015 = (run < 280000  ? kTRUE : kFALSE);
   if (is2015) {
      Int_t zeroBin = vtxDistrib->GetXaxis()->FindBin(0.0001);
      vtxDistrib->SetBinContent(zeroBin, 0.5*(vtxDistrib->GetBinContent(zeroBin-1)+vtxDistrib->GetBinContent(zeroBin+1)));
      vtxDistrib->SetBinError(zeroBin, 0.5*(vtxDistrib->GetBinError(zeroBin-1)+vtxDistrib->GetBinError(zeroBin+1)));
   }
   vtxDistrib->Fit("gaus", "ME", "Q", -15.0, 15.0);
   histList->Add(vtxDistrib);
   TF1* fit = (TF1*)vtxDistrib->GetListOfFunctions()->At(0);
   Double_t frac = fit->Integral(-10., +10.) / fit->Integral(-200., +200.);
   Double_t selected = fit->Integral(-10., +10.);
   Double_t selectedErr = fit->IntegralError(-10., +10.);
   Double_t total = fit->Integral(-200.0, 200.0);
   Double_t totalErr = fit->IntegralError(-200.0, 200.0);
   // compute "binomial" error for the fraction; assume the correlation factor is selected / total
   fracErr = TMath::Sqrt(TMath::Abs(total*total*selectedErr*selectedErr+selected*selected*totalErr*totalErr-2.0*frac*selected*total*selectedErr*totalErr))/total/total;
   chi2 = fit->GetChisquare() / fit->GetNDF();
   mean = fit->GetParameter(1);
   meanErr = fit->GetParError(1);
   sigma = fit->GetParameter(2);
   sigmaErr = fit->GetParError(2);
   return frac;
}



//________________________________________________________________________________________
void SetHistograms(const Char_t* setting, const Char_t* mcSignal) {
   //
   // signal extraction
   //
   THnF* seos = (THnF*)gHistManData->GetHistogram(Form("PairSEPM_%s",setting),"PairInvMass");
   THnF* meos = (THnF*)gHistManData->GetHistogram(Form("PairMEPM_%s",setting),"PairInvMass");
   THnF* sels1 = (THnF*)gHistManData->GetHistogram(Form("PairSEPP_%s",setting),"PairInvMass");
   THnF* sels2 = (THnF*)gHistManData->GetHistogram(Form("PairSEMM_%s",setting),"PairInvMass");
   THnF* mels1 = (THnF*)gHistManData->GetHistogram(Form("PairMEPP_%s",setting),"PairInvMass");
   THnF* mels2 = (THnF*)gHistManData->GetHistogram(Form("PairMEMM_%s",setting),"PairInvMass");
   
   THnF* seosMC = 0x0;
   if(gHistManMC)
      seosMC = (THnF*)gHistManMC->GetHistogram(Form("PairSEPM_%s_%s",setting, mcSignal),"PairInvMass");
   
   gFits.SetSEOSHistogram(seos);
   gFits.SetMEOSHistogram(meos);
   gFits.SetSELSHistograms(sels1, sels2);
   gFits.SetMELSHistograms(mels1, mels2);
   if(seosMC) gFits.SetSEOSMCHistogram(seosMC);
}


//________________________________________________________________________________________
void ComputeEfficiency(const Char_t* setting, const Char_t* mcSignalElectron, const Char_t* mcSignalJpsi, 
                       Double_t& eff, Double_t& effErr, Double_t& sigFrac, 
                       Double_t minMass, Double_t maxMass, 
                       Double_t ptMin=-1.0, Double_t ptMax=-1.0) {
//
// compute efficiency
//
   //cout << "ComputeEff 0" << endl;
   THnF* nominatorMC = (THnF*)gHistManMC->GetHistogram(Form("PairSEPM_%s_%s",setting, mcSignalElectron),"PairInvMass");
   //cout << "ComputeEff 0.1" << endl;
   TH1F* denominatorMC = (TH1F*)gHistManMC->GetHistogram(Form("%s_PureMCTruth_BeforeSelection", mcSignalJpsi),"PtMC");
   //cout << "ComputeEff 0.2" << endl;
   nominatorMC->Sumw2();
   //cout << "ComputeEff 0.3" << endl;
   //denominatorMC->Sumw2();

   //cout << "ComputeEff 1" << endl;

   if (ptMin >= 0.0 && ptMax >= 0.0)
      nominatorMC->GetAxis(1)->SetRangeUser(ptMin+1.0e-6, ptMax-1.0e-6);
   else 
      nominatorMC->GetAxis(1)->SetRangeUser(0.0+1.0e-6, 15.0-1.0e-6);
      //denominatorMC->GetXaxis()->SetRangeUser(ptMin+1.0e-6, ptMax-1.0e-6);

   Double_t err=0.0;
   sigFrac = 0.0;
   
   TH1D* nomProjMass = nominatorMC->Projection(0); nomProjMass->SetName("nomProjMass");
   TH1D* nomProjPt = nominatorMC->Projection(1); nomProjPt->SetName("nomProjPt");

   Double_t nomYieldErr = 0.0;
   //Double_t nomYield = nomProjPt->IntegralAndError(1, nomProjPt->GetXaxis()->GetNbins(), nomYieldErr);
   Double_t nomYield = nomProjMass->IntegralAndError(nomProjMass->GetXaxis()->FindBin(minMass+0.001), 
      nomProjMass->GetXaxis()->FindBin(maxMass-0.001), nomYieldErr);
   Double_t totalRecYieldErr = 0.0;
   Double_t totalRecYield = nomProjMass->IntegralAndError(1, nomProjMass->GetXaxis()->GetNbins(), totalRecYieldErr);
   sigFrac = nomYield / totalRecYield;

   Double_t denomYieldErr = 0.0;
   //Double_t denomYield = denomProjPt->IntegralAndError(1, nomProjPt->GetXaxis()->GetNbins(), denomYieldErr);
   Double_t denomYield = 0.0;
   if (ptMin >= 0.0 && ptMax >= 0.0)
      denomYield = denominatorMC->IntegralAndError(denominatorMC->GetXaxis()->FindBin(ptMin+1.0e-6), denominatorMC->GetXaxis()->FindBin(ptMax-1.0e-6), denomYieldErr);
   else denomYield = denominatorMC->IntegralAndError(1, denominatorMC->GetXaxis()->GetNbins(), denomYieldErr);
   
   // compute approximate binomial error (efficiency ~ 10%)
   if(denomYield>0.) {
      effErr = nomYield * (denomYield - nomYield) / denomYield / denomYield / denomYield;
      effErr = TMath::Sqrt(effErr);
      eff = nomYield / denomYield;
   }
   else {
      effErr = 0.0;
      eff = 0.0;
   }

   delete nomProjMass; 
   delete nomProjPt; 
}

// ________________________________________________________________________________________
TList* ComputeEfficiency(TString* settings, Int_t nSettings, const Char_t* dirMC, TH1D* hEvents, 
   const Char_t* mcSignalElectron, const Char_t* mcSignalJpsi, 
   Double_t eff[][gkNSigCountingVariations], Double_t effErr[][gkNSigCountingVariations], Double_t sigFrac[][gkNSigCountingVariations], 
   Bool_t runSigExtrSyst, Double_t ptMin=-1.0, Double_t ptMax=-1.0) {
   // 
   //  compute efficiency using the run-by-run weighting
   //  to save on time and CPU, when the MC file for a given run is opened, we compute the efficiency for all 
   //  signal counting window variations and all settings
   // 
   if (!hEvents) return 0x0;
   TList* masterList = new TList();
   masterList->SetOwner(kTRUE);
   masterList->SetName(Form("%s", mcSignalElectron));
   
   // store eff vs run just for the standard signal counting window interval 
   TH1D* effVsRun[64];
   TH1D* sigFracVsRun[64];
   for(Int_t i = 0; i<nSettings; ++i) {
      effVsRun[i] = (TH1D*)hEvents->Clone(Form("effVsRun_%s_%s", settings[i].Data(), mcSignalElectron));
      effVsRun[i]->SetDirectory(0x0);
      masterList->Add(effVsRun[i]);
      sigFracVsRun[i] = (TH1D*)hEvents->Clone(Form("sigFracVsRun_%s_%s", settings[i].Data(), mcSignalElectron));
      sigFracVsRun[i]->SetDirectory(0x0);
      masterList->Add(sigFracVsRun[i]);
   }
      
   Double_t norm = 0.0;
   for(Int_t iSetting = 0; iSetting<64; iSetting++) {
      for(Int_t i=0; i<gkNSigCountingVariations; ++i) {
         eff[iSetting][i] = 0.0;
         effErr[iSetting][i] = 0.0;
         sigFrac[iSetting][i] = 0.0;
      }
   }
   
   TH1D* effVsPtAverage[64] = {0x0};
   
   for(Int_t iRun = 1; iRun <= hEvents->GetXaxis()->GetNbins(); ++iRun) {
      Double_t nevents = hEvents->GetBinContent(iRun);
      if (nevents<1.0) continue;
      TString runStr = hEvents->GetXaxis()->GetBinLabel(iRun);
      Int_t run = runStr.Atoi();
      
      Double_t effRun[64][gkNSigCountingVariations], effRunErr[64][gkNSigCountingVariations], sigFracRun[64][gkNSigCountingVariations]; 
      TList* outListForRun = ComputeEfficienciesForRun(settings, nSettings, dirMC, run, mcSignalElectron, mcSignalJpsi, effRun, effRunErr, sigFracRun, runSigExtrSyst, ptMin, ptMax);
      outListForRun->SetName(Form("Run_%d", run));
      masterList->Add(outListForRun);
      
      for(Int_t iSetting = 0; iSetting<nSettings; ++iSetting) {
         // sum the eff vs pt histograms over runs
         if(!effVsPtAverage[iSetting]) { 
            TH1D* tempHist = (TH1D*)outListForRun->FindObject(Form("effVsPt_run%d_%s_msig0", run, settings[iSetting].Data()));
            effVsPtAverage[iSetting] = (TH1D*)tempHist->Clone(Form("effVsPtAverage_%s_%s", settings[iSetting].Data(), mcSignalElectron));
            effVsPtAverage[iSetting]->Scale(nevents);
         }
         else 
            effVsPtAverage[iSetting]->Add((TH1D*)outListForRun->FindObject(Form("effVsPt_run%d_%s_msig0", run, settings[iSetting].Data())), nevents);

         // set the contents of the eff vs run histograms, just for the standard sig extraction setting
         effVsRun[iSetting]->SetBinContent(iRun, effRun[iSetting][0]);
         effVsRun[iSetting]->SetBinError(iRun, effRunErr[iSetting][0]);
         sigFracVsRun[iSetting]->SetBinContent(iRun, sigFracRun[iSetting][0]);
         
         // sum the eff values over runs for all sig extraction settings
         for(Int_t i=0; i<gkNSigCountingVariations; ++i) {
            eff[iSetting][i] += nevents*effRun[iSetting][i];
            sigFrac[iSetting][i] += nevents*sigFracRun[iSetting][i];
            effErr[iSetting][i] += nevents*nevents*effRunErr[iSetting][i]*effRunErr[iSetting][i];
         }
      }
               
      norm += nevents;
   }
   
   // normalize everything
   for(Int_t iSetting = 0; iSetting<nSettings; iSetting++) {
      for (Int_t i=0; i<gkNSigCountingVariations; ++i) {
         eff[iSetting][i] /= norm;
         effErr[iSetting][i] = TMath::Sqrt(effErr[iSetting][i]);
         effErr[iSetting][i] /= norm;
         sigFrac[iSetting][i] /= norm;
      }
      cout << "eff setting " << settings[iSetting].Data() << " " << eff[iSetting][0] << " +/- " << effErr[iSetting][0] << endl;
      effVsPtAverage[iSetting]->Scale(1./norm);
      masterList->Add(effVsPtAverage[iSetting]);
      
      effVsRun[iSetting]->SetTitle(Form("Efficiency vs run, mass range: %.2f -> %.2f GeV/c^{2}, pt range: %.2f -> %.2f GeV/c, integrated eff = %.4e #pm %.4e", gmSigLims[0][0], gmSigLims[0][1], ptMin, ptMax, eff[iSetting][0], effErr[iSetting][0]));
      sigFracVsRun[iSetting]->SetTitle(Form("Signal fraction vs run, mass range: %.2f -> %.2f GeV/c^{2}, pt range: %.2f -> %.2f GeV/c, integrated sig fraction = %.4e", gmSigLims[0][0], gmSigLims[0][1], ptMin, ptMax, sigFrac[iSetting][0]));
   }
   
   cout << "compute efficiency" << endl;
   return masterList;
}

//________________________________________________________________________________________
TList* ComputeEfficienciesForRun(TString* settings, Int_t nSettings, const Char_t* dirMC, Int_t run, 
                              const Char_t* mcSignalElectron, const Char_t* mcSignalJpsi, 
                              Double_t eff[][gkNSigCountingVariations], Double_t effErr[][gkNSigCountingVariations], Double_t sigFrac[][gkNSigCountingVariations], 
                              Bool_t runSigExtrSyst, Double_t ptMin=-1.0, Double_t ptMax=-1.0) {
   //
   // compute efficiencies for a given run number
   //
   
   cout << "ComputeEfficiencyForRun :: Run " << run << " *** ptMin / ptMax " << ptMin << "/" << ptMax << endl;
   gHistManMCPerRun = new AliHistogramManager();
   gHistManMCPerRun->InitFile(Form("%s/%s/runOutput/%09d/AnalysisHistograms_jpsi2ee_pp2017_MC.root", dirMC, (GetPeriodFromRun(run) == 0 ? "LHC15n" : "LHC17pq"), run), "jpsi2eeHistos");
   gHistManMCPerRun->SetName(Form("histMan_run%d", run));
   
   TList* returnList = new TList();
   returnList->SetOwner(kTRUE);
   
   TH1F* denominatorMC_fromFile = (TH1F*)gHistManMCPerRun->GetHistogram(Form("%s_PureMCTruth_BeforeSelection", mcSignalJpsi),"PtMC");  
   
   // temporary solution : needed because the denominator has a different binning wrt nominator
   TH1D* denominatorMC = new TH1D("denominatorMC", "", 400, 0.0, 20.0);
   for (Int_t i = 1;i <= 400; ++i) {
      denominatorMC->SetBinContent(i, denominatorMC_fromFile->GetBinContent(i));
      denominatorMC->SetBinError(i, denominatorMC_fromFile->GetBinError(i));
   }
   denominatorMC->SetDirectory(0x0);
   
   Double_t denomYieldErr = 0.0;
         
   TH1F* denominatorMCclone = (TH1F*)denominatorMC_fromFile->Clone(Form("PtMCdenominator_run%d_%s", run, mcSignalElectron));
   denominatorMCclone->SetDirectory(0x0);
   returnList->Add(denominatorMCclone);
   
   // function to be used for the pt reweighting -> fitted to the MC anchored in LHC17pq (natural spectrum)
   TF1* fpl= new TF1("fpl","[0]*x/TMath::Power((1+TMath::Power(x/[1],[3])),[2])",0.,30.);
   fpl->SetParameters(1.0, 3.67987, 3.19625, 2.0);       //  these are the parameters of the fit to the 2017 MC pt distribution
      
   for(Int_t iSetting = 0; iSetting<nSettings; ++iSetting) {
      THnF* nominatorMC = (THnF*)gHistManMCPerRun->GetHistogram(Form("PairSEPM_%s_%s",settings[iSetting].Data(), mcSignalElectron),"PairInvMass");
      nominatorMC->Sumw2();
   
      // set a fiducial pt cut of 0 - 15.0 GeV/c if no pt range was specified  
      Int_t originalNbins = nominatorMC->GetAxis(1);
      if(ptMin < 0.0 && ptMax < 0.0)
         nominatorMC->GetAxis(1)->SetRangeUser(0.0+1.0e-6, 15.0-1.0e-6);
      else
         nominatorMC->GetAxis(1)->SetRangeUser(ptMin+1.0e-6, ptMax-1.0e-6);
   
      TH1D* nomProjMass = nominatorMC->Projection(0); 
      // go back to the original bin range
      nominatorMC->GetAxis(1)->UnZoom(); 
      
      nomProjMass->SetName(Form("ProjMass_run%d_%s", run, settings[iSetting].Data()));
      nomProjMass->SetDirectory(0x0);
      Double_t allRecYieldErr = 0.0;
      Double_t allRecYield = nomProjMass->IntegralAndError(1, nomProjMass->GetXaxis()->GetNbins(), allRecYieldErr);
         
      for(Int_t isig = 0; isig<(runSigExtrSyst ? gkNSigCountingVariations : 1); ++isig) {
         // Compute the signal fraction in the mass counting window
         Double_t countedRecYieldErr = 0.0;
         Double_t countedRecYield = nomProjMass->IntegralAndError(nomProjMass->GetXaxis()->FindBin(gmSigLims[isig][0]+0.001), 
                                                                  nomProjMass->GetXaxis()->FindBin(gmSigLims[isig][1]-0.001), countedRecYieldErr);
         sigFrac[iSetting][isig] = countedRecYield / allRecYield;
         // Compute the efficiency
         // +++++++++++ no pt reweighting
         /*Double_t denomYieldErr = 0.0;
         Double_t denomYield = 0.0;
         if(ptMin >= 0.0 && ptMax >= 0.0) 
            denomYield = denominatorMC->IntegralAndError(denominatorMC->GetXaxis()->FindBin(ptMin+1.0e-6), denominatorMC->GetXaxis()->FindBin(ptMax-1.0e-6), denomYieldErr);
         else 
            denomYield = denominatorMC->IntegralAndError(1, denominatorMC->GetXaxis()->FindBin(24.999), denomYieldErr);
         // compute approximate binomial error (efficiency ~ 10%)
         if(denomYield>0.) {
            effErr = countedRecYield * (denomYield - countedRecYield) / denomYield / denomYield / denomYield;
            effErr = TMath::Sqrt(effErr);
            eff = countedRecYield / denomYield;
         }
         else {
            effErr = 0.0;
            eff = 0.0;
         } */
      
         // ++++++++++ pt reweighting
         nominatorMC->GetAxis(0)->SetRangeUser(gmSigLims[isig][0]+1.0e-6, gmSigLims[isig][1]-1.0e-6);
         TH1D* nomProjPt = nominatorMC->Projection(1); 
         nomProjPt->SetName(Form("ProjPt_run%d_%s_msig%d", run, settings[iSetting].Data(), isig));
         nomProjPt->SetDirectory(0x0);
         TH1D* effVsPt = (TH1D*)nomProjPt->Clone(Form("effVsPt_run%d_%s_msig%d", run, settings[iSetting].Data(), isig)); 
         effVsPt->SetDirectory(0x0);
         effVsPt->Reset();
         effVsPt->Divide(nomProjPt, denominatorMC, 1.0, 1.0, "B");
         eff[iSetting][isig] = WeightedAverage(effVsPt, fpl, (ptMin>0.0 ? ptMin : 0.0), (ptMax>0.0 ? ptMax : 15.0), effErr[iSetting][isig]);
         if (isig == 0) {
            returnList->Add(nomProjMass);
            returnList->Add(nomProjPt);
            returnList->Add(effVsPt);
         }
      }                                                     //  end loop over signal counting windows
   }                                                        // end loop over settings
   
   gHistManMCPerRun->CloseFile();
   delete denominatorMC;
   delete gHistManMCPerRun;
   return returnList;
}

//____________________________________________________________________________
Double_t WeightedAverage(TH1D* eff, TF1* weight, Float_t ptMin, Float_t ptMax, Double_t& err) {
   //
   // Takes as arguments the efficiency as a funciton of pt and the input pt spectrum.
   // Calculates the average efficiency in the pt range ptMin - ptMax
   //
   Double_t average = 0.0;
   Double_t norm = 0.0;
   err = 0.0;
   for(Int_t i=1; i<=eff->GetXaxis()->GetNbins(); ++i) {
      Double_t pt=eff->GetXaxis()->GetBinCenter(i);
      if(ptMin>0.0 && pt<ptMin) continue;
      if(ptMax>0.0 && pt>ptMax) continue;
      average += eff->GetBinContent(i) * weight->Eval(pt) * eff->GetXaxis()->GetBinWidth(i);
      err += TMath::Power(eff->GetBinError(i) * weight->Eval(pt) * eff->GetXaxis()->GetBinWidth(i), 2.0);
      norm += weight->Eval(pt) * eff->GetXaxis()->GetBinWidth(i);
   }
   average /= norm;
   err = TMath::Sqrt(err); 
   err /= norm;
   return average;
}

//____________________________________________________________________________
Double_t WeightedAverage(TH2D* eff, TF1* weightPt, Double_t alpha, Float_t ptMin, Float_t ptMax, Bool_t nonPromptFixedPolarization = kTRUE) {
   //
   // Takes as arguments a 2-D efficiency as a function of pt and cos(theta*), the pt spectrum and the cos(theta*) distribution
   // Calculates the average efficiency between ptMin and ptMax
   //
   TF1* fbData = new TF1("fbData", "pol3", 0.0, 10.0);
   fbData->SetParameters(0.0949766, 0.001, 0.001, 0.0001);
   
   // When nonPromptFixedPolarization is kTRUE, then for the non-prompt component use a fixed polarization of alpha=-0.46
   Double_t alphaBABAR = -0.46;
   
   Double_t average = 0.0;
   Double_t norm = 0.0;
   for(Int_t ipt=1; ipt<=eff->GetYaxis()->GetNbins(); ++ipt) {
      Double_t pt=eff->GetYaxis()->GetBinCenter(ipt);
      if(pt<ptMin) continue;
      if(pt>ptMax) continue;
      for(Int_t ict=1; ict<=eff->GetXaxis()->GetNbins(); ++ict) {
         Double_t cTheta=eff->GetXaxis()->GetBinCenter(ict);
         if(TMath::Abs(cTheta)>1.0) continue;
         Double_t wPolariz = 1.0+alpha*cTheta*cTheta;
         if(nonPromptFixedPolarization) 
            wPolariz = 1.0+(alpha*(1.0-fbData->Eval(pt)) + alphaBABAR*fbData->Eval(pt))*cTheta*cTheta;
         average += eff->GetBinContent(ict,ipt) * weightPt->Eval(pt) * wPolariz;
         norm += weightPt->Eval(pt) * wPolariz;
       }
    }
   average /= norm;
   return average;
}


//________________________________________________________________________________________
void BeutifyTH1(TH1* h, const Char_t* title, Double_t lineWidth, Int_t lineColor,
                Int_t markerStyle=1, Int_t markerColor=1) {
   //
   // set drawing options for a TH1
   //
   if(!h) return;
   h->SetTitle(title);
   h->SetLineWidth(lineWidth); h->SetLineColor(lineColor);
   h->SetMarkerStyle(markerStyle);
   h->SetMarkerColor(markerColor);
}
                
//________________________________________________________________________________________
void BeutifyTAxis(TAxis* ax, 
                  const Char_t* title, Bool_t centerTitle, Double_t titleSize, Double_t titleOffset, Int_t titleFont,
                  Double_t labelSize, Int_t labelFont, Int_t nDivisions) {
//
// set drawing options for TAxis
//
   if(!ax) return;
   ax->SetTitle(title);
   if(centerTitle) ax->CenterTitle();
   ax->SetTitleSize(titleSize);
   ax->SetTitleOffset(titleOffset);
   ax->SetTitleFont(titleFont);
   ax->SetLabelSize(labelSize);
   ax->SetLabelFont(labelFont);
   ax->SetNdivisions(nDivisions);
}
                                  
                                  
//________________________________________________________________________________________
void Draw(AliResonanceFits* fits, Double_t massSignalMin, Double_t massSignalMax, 
          Bool_t savePicture=kFALSE, const Char_t* outputDir=".", const Char_t* pictureName="standard", Bool_t preliminaryPlot=kFALSE) {
//
// Make signal extraction plot using the AliResonanceFits output
//
   TH1* histSplusB = fits->GetSplusB(); if(!histSplusB) return;
   TH1* histSig = fits->GetSignal();  if(!histSig) return;
   TH1* histBkg = fits->GetBkg(); // if(!histBkg) return;
   TH1* histSplusResidualBkg = fits->GetSplusResidualBkg();
   TH1* histSigMC = fits->GetSignalMC();
   TF1* bkgFit = ((fits->GetBkgMethod()==AliResonanceFits::kBkgMixedEventAndResidualFit || 
   fits->GetBkgMethod()==AliResonanceFits::kBkgFitFunction) ? fits->GetBkgFitFunction() : 0x0);
   TF1* globalFit = ((fits->GetBkgMethod()==AliResonanceFits::kBkgMixedEventAndResidualFit || 
   fits->GetBkgMethod()==AliResonanceFits::kBkgFitFunction) ? fits->GetGlobalFitFunction() : 0x0);
                            
   Bool_t histogramsAreUnchanged = kFALSE;
   if(histSplusB==gOldPointer) histogramsAreUnchanged = kTRUE;
   else gOldPointer = histSplusB;
                                     
   TLatex* lat=new TLatex();
   lat->SetNDC();
   lat->SetTextSize(0.06);
   lat->SetTextColor(1);
                                     
   TCanvas* oldCanvas = (TCanvas*)gROOT->FindObject("AliResonanceFitsCanvas");
   if(oldCanvas) delete oldCanvas; 
   TCanvas* c1=new TCanvas("AliResonanceFitsCanvas", "AliResonanceFits canvas", 980, 1200);
   c1->SetTopMargin(0.01);
   c1->SetRightMargin(0.005);
   c1->SetBottomMargin(0.01);
   c1->SetLeftMargin(0.005);
   c1->Divide(1,2,0.0,0.0);
                                     
   // Draw the top pad ================================================
   TVirtualPad* pad = c1->cd(1);
   pad->SetTopMargin(0.01);
   pad->SetRightMargin(0.02);
   pad->SetBottomMargin(0.0);
   pad->SetLeftMargin(0.15);
   pad->SetTicks(1,1);
                                     
   histSplusB->SetStats(kFALSE);
   if(!histogramsAreUnchanged)
   histSplusB->GetYaxis()->SetRangeUser(0.01, 1.3*histSplusB->GetMaximum());
   //histSplusB->GetXaxis()->SetRangeUser(fits->GetMassFitRange()[0]-0.1, fits->GetMassFitRange()[1]+0.1);
   BeutifyTH1(histSplusB, "", 3.0, 1);
   BeutifyTAxis(histSplusB->GetXaxis(), "", kFALSE, 0.04, 1., 42, 0.04, 42, 510);
   BeutifyTAxis(histSplusB->GetYaxis(), Form("entries per %.0f MeV/#it{c}^{2}", 1000.*histSplusB->GetXaxis()->GetBinWidth(1)), kTRUE, 0.08, 0.84, 42, 0.065, 42, 507);
   histSplusB->DrawClone("EXY");
   
   if(histBkg) {
      histBkg->SetStats(kFALSE);
      BeutifyTH1(histBkg, "", 2.0, 2);
      //histBkg->GetXaxis()->SetRangeUser(fits->GetMassFitRange()[0]-0.1, fits->GetMassFitRange()[1]+0.1);
      histBkg->DrawClone("sameE");
   }
   
   if(bkgFit && globalFit && fits->GetBkgMethod()==AliResonanceFits::kBkgFitFunction) {
      bkgFit->SetLineColor(4); bkgFit->SetLineWidth(1);
      globalFit->SetLineColor(6); globalFit->SetLineWidth(1);
      globalFit->Draw("same");
      bkgFit->Draw("same");
   }
   
   /*
   cout << "Draw 3.2" << endl;
   TBox* boxL = new TBox(fits->GetMassFitRange()[0], histSplusB->GetMaximum()*0.1, fits->GetMassExclusionRange()[0], histSplusB->GetMaximum()*0.15);
   boxL->SetFillColorAlpha(4, 0.35);
   boxL->SetFillStyle(3003);
   boxL->Draw();
   
   cout << "Draw 3.3" << endl;
   TBox* boxR = new TBox(fits->GetMassExclusionRange()[1], histSplusB->GetMaximum()*0.1, fits->GetMassFitRange()[1], histSplusB->GetMaximum()*0.15);
   boxR->SetFillColorAlpha(4, 0.35);
   boxR->SetFillStyle(3003);
   boxR->Draw();
   cout << "Draw 4" << endl;
   */
   
   TLegend* legend1 = new TLegend(0.18, 0.55, 0.42, 0.73);
   legend1->SetTextSize(0.06);
   legend1->SetTextFont(42);
   legend1->SetFillColor(0);
   legend1->SetBorderSize(0);
   legend1->AddEntry(histSplusB, "Same event", "l");
   if(fits->GetBkgMethod()==AliResonanceFits::kBkgFitFunction) {
      legend1->AddEntry(globalFit, "Global fit", "l");
      legend1->AddEntry(bkgFit, "Bkg fit", "l");
   }
   else if(fits->GetBkgMethod()==AliResonanceFits::kBkgMixedEventAndResidualFit) {
      legend1->AddEntry(histBkg, "ME like sign", "l");
   }
   else 
      legend1->AddEntry(histBkg, (fits->GetBkgMethod()==AliResonanceFits::kBkgMixedEvent ? "Mixed event" : "Like sign"), "l");
   legend1->Draw();
                                     
   TLine line;
   line.SetLineColor(1);
   line.SetLineStyle(1);
   line.SetLineWidth(1.0);
   line.DrawLine(massSignalMin, 0.0, massSignalMin, histSplusB->GetMaximum()*0.75);
   line.DrawLine(massSignalMax, 0.0, massSignalMax, histSplusB->GetMaximum()*0.75);
   lat->SetTextFont(42);
   lat->SetTextColor(1);
   lat->SetTextSize(0.060);
   
   if(preliminaryPlot) {
      lat->DrawLatex(0.2, 0.90, "ALICE Preliminary");
      lat->DrawLatex(0.2, 0.82, "0#minus90% Xe#minusXe #sqrt{#it{s}_{NN}} = 5.44 TeV");
   }
                                     
   if(!preliminaryPlot && (fits->GetBkgMethod()==AliResonanceFits::kBkgMixedEvent || fits->GetBkgMethod()==AliResonanceFits::kBkgLikeSign))
      lat->DrawLatex(0.56, 0.82, Form("#chi^{2}/NDF (%.1f, %.1f) = %.2f", fits->GetMassFitRange()[0], fits->GetMassFitRange()[1], fits->GetFitValues()[AliResonanceFits::kChisqSideBands]));
   //if(fEventVsCent) lat->DrawLatex((noYlabels ? 0.52 : 0.56), 0.72, Form("# events = %.2e", fEventVsCent->Integral(fEventVsCent->GetXaxis()->FindBin(fCentralityRange[0]+0.001), fEventVsCent->GetXaxis()->FindBin(fCentralityRange[1]-0.001))));
   /*lat->DrawLatex(0.63, 0.66, "Fit ranges:");
    *  lat->DrawLatex(0.63, 0.60, Form("%.2f<m_{ee}<%.2f GeV/c^{2}", fFitRange[0], fFitRange[1]));
    *  lat->DrawLatex(0.63, 0.54, Form("%.2f<p_{T}<%.2f GeV/c", fPtFitRange[0], fPtFitRange[1]));
    *  lat->DrawLatex(0.63, 0.45, "Signal:");
    *  lat->DrawLatex(0.63, 0.39, Form("%.2f<m_{ee}<%.2f GeV/c^{2}", fSignalRange[0], fSignalRange[1]));
    *  lat->DrawLatex(0.63, 0.33, Form("%.2f<p_{T}<%.2f GeV/c", fPtRange[0], fPtRange[1])); */
                                  
   // Draw the bottom pad ================================================
   pad = c1->cd(2);
   pad->SetTopMargin(0.0);
   pad->SetRightMargin(0.02);
   pad->SetBottomMargin(0.17);
   pad->SetLeftMargin(0.15);
   pad->SetTicks(1,1);
   
   if(fits->GetBkgMethod()==AliResonanceFits::kBkgMixedEvent ||
      fits->GetBkgMethod()==AliResonanceFits::kBkgLikeSign ||
      fits->GetBkgMethod()==AliResonanceFits::kBkgFitFunction || 
      fits->GetBkgMethod()==AliResonanceFits::kBkgMixedEventAndResidualFit) {
      histSig->SetStats(kFALSE);
      BeutifyTH1(histSig, "", 2.0, 2, 20, 2);
      if(!histogramsAreUnchanged)
         histSig->GetYaxis()->SetRangeUser(histSig->GetMinimum()*1.4, histSig->GetMaximum()*1.8);
         //histSig->GetXaxis()->SetRangeUser(fits->GetMassFitRange()[0]-0.1, fits->GetMassFitRange()[1]+0.1);
         BeutifyTAxis(histSig->GetYaxis(), "", kFALSE, 0.08, 1.2, 42, 0.065, 42, 508);
         BeutifyTAxis(histSig->GetXaxis(), "#it{m}_{e^{#plus}e^{#minus}} (GeV/#it{c}^{2})", kTRUE, 0.08, 0.84, 42, 0.065, 42, 508);
         histSig->DrawClone("XY");
         if(histSigMC) {
         //histSigMC->GetXaxis()->SetRangeUser(fits->GetMassFitRange()[0]-0.1, fits->GetMassFitRange()[1]+0.1);
            BeutifyTH1(histSigMC, "", 2.0, 1);
            TH1* histSigMCclone = (TH1*)histSigMC->Clone("histSigMCclone");
            if(fits->GetBkgMethod()==AliResonanceFits::kBkgFitFunction) 
               histSigMCclone->Scale(globalFit->GetParameter(0));
            histSigMCclone->DrawClone("sameHISTC"); //smooth MC shape drawing
         }
      }
      if(fits->GetBkgMethod()==AliResonanceFits::kBkgMixedEventAndResidualFit) {
         histSplusResidualBkg->SetStats(kFALSE);
         BeutifyTH1(histSplusResidualBkg, "", 2.0, 2, 20, 2);
         if(!histogramsAreUnchanged)
            histSplusResidualBkg->GetYaxis()->SetRangeUser(histSplusResidualBkg->GetMinimum()*1.4, histSplusResidualBkg->GetMaximum()*1.8);
         BeutifyTAxis(histSplusResidualBkg->GetYaxis(), "", kFALSE, 0.08, 1.2, 42, 0.065, 42, 508);
         BeutifyTAxis(histSplusResidualBkg->GetXaxis(), "#it{m}_{e^{#plus}e^{#minus}} (GeV/#it{c}^{2})", kTRUE, 0.08, 0.84, 42, 0.065, 42, 508);
         histSplusResidualBkg->DrawClone("XY");
         bkgFit->SetLineColor(4); bkgFit->SetLineWidth(2);
         globalFit->SetLineColor(6); globalFit->SetLineWidth(2);
                                           
         globalFit->Draw("same");
         bkgFit->Draw("same");
      }
      line.DrawLine(fits->GetMassFitRange()[0],0.0,fits->GetMassFitRange()[1],0.0);
                                        
      Double_t* fitValues = fits->GetFitValues();
                                        
      line.DrawLine(massSignalMin, histSig->GetMinimum()*0.9, massSignalMin, histSig->GetMaximum()*0.8);
      line.DrawLine(massSignalMax, histSig->GetMinimum()*0.9, massSignalMax, histSig->GetMaximum()*0.8);
      lat->DrawLatex(0.18, 0.92, Form("Signal:    %.0f #pm %.0f", fitValues[AliResonanceFits::kSig], fitValues[AliResonanceFits::kSigErr]));
      lat->DrawLatex(0.18, 0.84, Form("#it{S}/#it{B}:        %.3f #pm %.3f", fitValues[AliResonanceFits::kSoverB], fitValues[AliResonanceFits::kSoverBerr]));
      if(fits->GetBkgMethod()==AliResonanceFits::kBkgMixedEvent || fits->GetBkgMethod()==AliResonanceFits::kBkgMixedEventAndResidualFit) 
         lat->DrawLatex(0.18, 0.76, Form("#it{S}/#sqrt{#it{S}+#it{B}}: %.1f", fitValues[AliResonanceFits::kSignif]));
      if(fits->GetBkgMethod()==AliResonanceFits::kBkgLikeSign) 
         lat->DrawLatex(0.18, 0.76, Form("#it{S}/#sqrt{#it{S}+2#it{B}}: %.1f", fitValues[AliResonanceFits::kSignif]));
      if(fits->GetBkgMethod()==AliResonanceFits::kBkgFitFunction) 
         lat->DrawLatex(0.18, 0.76, Form("#it{S}/#delta_{#it{S}}: %.1f", fitValues[AliResonanceFits::kSignif]));
      if(histSigMC && !preliminaryPlot) {
         lat->DrawLatex(0.18, 0.68, Form("#chi^{2}/NDF (%.1f, %.1f) = %.2f", fits->GetMassFitRange()[0], fits->GetMassFitRange()[1], fitValues[AliResonanceFits::kChisqMCTotal]));
         lat->DrawLatex(0.18, 0.60, Form("Probability = %.2f", fitValues[AliResonanceFits::kFitProbability]));
         //lat->DrawLatex(0.56, 0.56, Form("MC yield fraction  = %.2f", fitValues[AliResonanceFits::kMCYieldFraction]));
      }
      TLegend* legend2 = new TLegend(0.70, 0.83, 0.96, 0.97);
      legend2->SetFillColor(0);
      legend2->SetBorderSize(0);
      legend2->SetTextFont(42);
      legend2->SetTextSize(0.06);
      legend2->SetMargin(0.15);
      if(fits->GetBkgMethod()==AliResonanceFits::kBkgMixedEventAndResidualFit)
         legend2->AddEntry(histSplusResidualBkg, "Data", "lp");
      else
         legend2->AddEntry(histSig, "Data", "lp");
      if(fits->GetBkgMethod()==AliResonanceFits::kBkgMixedEventAndResidualFit) {
         legend2->AddEntry(globalFit, "global function", "l");
         legend2->AddEntry(bkgFit, "bkg function", "l");
      }
      else
         if(histSigMC) legend2->AddEntry(histSigMC, "Monte Carlo", "l");
      legend2->Draw();  
                                        
      if(savePicture) {
         TString pictureDetailedName = pictureName;
         if(fits->GetBkgMethod()==AliResonanceFits::kBkgMixedEvent) pictureDetailedName += "_BkgME";
         if(fits->GetBkgMethod()==AliResonanceFits::kBkgLikeSign) pictureDetailedName += "_BkgLS";
         if(fits->GetBkgMethod()==AliResonanceFits::kBkgMixedEventAndResidualFit) pictureDetailedName += "_BkgMEresidual";
         if(fits->GetBkgMethod()==AliResonanceFits::kBkgFitFunction) pictureDetailedName += "_BkgFit";
         if(fits->GetScalingOption() == AliResonanceFits::kScaleEntries) pictureDetailedName += "_ScaleEnt";
         if(fits->GetScalingOption() == AliResonanceFits::kScaleWeightedAverage) pictureDetailedName += "_ScaleWav";
         if(fits->GetScalingOption() == AliResonanceFits::kScaleFit) pictureDetailedName += "_ScaleFit";
         if(fits->GetBkgMethod()==AliResonanceFits::kBkgMixedEvent) {
         if(fits->GetMEMatchingMethod()==AliResonanceFits::kMatchSEOS) pictureDetailedName += "_MatchSEOS";
         if(fits->GetMEMatchingMethod()==AliResonanceFits::kMatchSELS) pictureDetailedName += "_MatchSELS";
      }
      if(fits->GetScalingOption() == AliResonanceFits::kScaleFit) {
         if(fits->GetMinuitFitOption() == AliResonanceFits::kMinuitMethodChi2) pictureDetailedName += "_FitChi2";
         if(fits->GetMinuitFitOption() == AliResonanceFits::kMinuitMethodLikelihood) pictureDetailedName += "_FitLikelihood";
      }
                                           
      c1->SaveAs(Form("%s/%s.png", outputDir, pictureDetailedName.Data()));
      //c1->SaveAs(Form("%s/%s.eps", outputDir, pictureDetailedName.Data()));
   }
}


//____________________________________________________________________________________________________
//TList* CheckPIDsystematicsForSetting(const Char_t* setting, const Char_t* refSetting="pid1") {
TList* CheckPIDsystematicsForSetting(TH1F* data_p_ref, TH1F* data_p, TH1F* mc_p_ref, TH1F* mc_p, const Char_t* setting) {
   //
   //
   //
   const Int_t kNrebinnedBins = 48;
   Double_t rebinEdges[kNrebinnedBins+1];
   Int_t bin = 0;
   for(Double_t p = 0.0; p<2.99; p += 0.1, bin++) rebinEdges[bin] = p;
   for(Double_t p = 3.0; p<4.99; p += 0.2, bin++) rebinEdges[bin] = p;
   for(Double_t p = 5.0; p<7.99; p += 0.5, bin++) rebinEdges[bin] = p;
   for(Double_t p = 8.0; p<10.01; p += 1.0, bin++) rebinEdges[bin] = p;
   cout << "bin  " << bin << endl;
   for(Int_t i=0; i<=kNrebinnedBins; ++i) cout << rebinEdges[i] << " ";
   cout << endl;
   cout << "1" << endl;
   
   //cout << "setting / refSetting  :: " << setting << " / " << refSetting << endl; 
   
   //TH1F* data_p_ref = (TH1F*)gHistManData->GetHistogram(Form("Track_%s", refSetting),"P");
   cout << "1.1" << endl;
   cout << data_p_ref << endl;
   //data_p_ref->SetName(Form("data_p_ref_%s_%.6f", refSetting, gRandom->Rndm()));
   TH1F* data_p_ref_rebinned = dynamic_cast<TH1F*>(data_p_ref->Rebin(kNrebinnedBins, Form("%s_rebinned_%.6f", data_p_ref->GetName(), gRandom->Rndm()), rebinEdges));
   
   cout << "2" << endl;
   //TH1F* mc_p_ref = (TH1F*)gHistManMC->GetHistogram(Form("Track_%s_TrueElectron", refSetting),"P");
   //mc_p_ref->SetName(Form("mc_p_ref_%s_%.6f", refSetting, gRandom->Rndm()));
   //TH1F* mc_eta_ref = (TH1F*)gHistManMC->GetHistogram(Form("Track_%s_TrueElectron", refSetting),"Eta");
   //mc_eta_ref->SetName(Form("mc_eta_ref_%s_%.6f", refSetting, gRandom->Rndm()));
   //TH2F* mc_etap_ref = (TH2F*)gHistManMC->GetHistogram(Form("Track_%s_TrueElectron", refSetting),"Eta_P");
   //mc_etap_ref->SetName(Form("mc_etap_ref_%s_%.6f", refSetting, gRandom->Rndm()));
   TH1F* mc_p_ref_rebinned = dynamic_cast<TH1F*>(mc_p_ref->Rebin(kNrebinnedBins, Form("%s_rebinned_%.6f", mc_p_ref->GetName(), gRandom->Rndm()), rebinEdges));
   //mc_p_ref_rebinned->SetName(Form("mc_p_ref_rebinned%.6f",gRandom->Rndm()));
   //mc_eta_ref->Rebin(2);

   cout << "3" << endl;
   //TH1F* data_p = (TH1F*)gHistManData->GetHistogram(Form("Track_%s", setting),"P");
   //data_p->SetName(Form("data_p_%s", setting));
   //TH1F* data_eta = (TH1F*)gHistManData->GetHistogram(Form("Track_%s", setting),"Eta");
   //data_eta->SetName(Form("data_eta_%s", setting));
   //TH2F* data_etap = (TH2F*)gHistManData->GetHistogram(Form("Track_%s", setting),"Eta_P");
   //data_etap->SetName(Form("data_etap_%s", setting));
   TH1F* data_p_rebinned = dynamic_cast<TH1F*>(data_p->Rebin(kNrebinnedBins, Form("%s_rebinned_%.6f", data_p->GetName(), gRandom->Rndm()), rebinEdges));
   //data_p_rebinned->SetName(Form("data_p_rebinned%.6f",gRandom->Rndm()));
   //data_eta->Rebin(2);

   cout << "4" << endl;
   //TH1F* mc_p = (TH1F*)gHistManMC->GetHistogram(Form("Track_%s_TrueElectron", setting),"P");
   //mc_p->SetName(Form("mc_p_%s", setting));
   //TH1F* mc_eta = (TH1F*)gHistManMC->GetHistogram(Form("Track_%s_TrueElectron", setting),"Eta");
   //mc_eta->SetName(Form("mc_eta_%s", setting));
   //TH2F* mc_etap = (TH2F*)gHistManMC->GetHistogram(Form("Track_%s_TrueElectron", setting),"Eta_P");
   //mc_etap->SetName(Form("mc_etap_%s", setting));
   TH1F* mc_p_rebinned = dynamic_cast<TH1F*>(mc_p->Rebin(kNrebinnedBins, Form("%s_rebinned_%.6f", mc_p->GetName(), gRandom->Rndm()), rebinEdges));
   //mc_p_rebinned->SetName(Form("mc_p_rebinned%.6f",gRandom->Rndm()));
   //mc_eta->Rebin(2);

   cout << "5" << endl;
   TH1F* data_ratio_p = (TH1F*)data_p_rebinned->Clone(Form("data_ratio_p_%s", setting));
   BeutifyTH1(data_ratio_p, "", 2, 4);
   TH1F* mc_ratio_p = (TH1F*)mc_p_rebinned->Clone(Form("mc_ratio_p_%s", setting));
   BeutifyTH1(mc_ratio_p, "", 2, 2);
   data_ratio_p->Divide(data_p_rebinned, data_p_ref_rebinned, 1.,1.,"B");
   mc_ratio_p->Divide(mc_p_rebinned, mc_p_ref_rebinned, 1.,1.,"B");
   
   cout << "6" << endl;
   TH1F* dataMC_ratio_p = (TH1F*)data_ratio_p->Clone(Form("dataMC_ratio_p_%s", setting));
   BeutifyTH1(dataMC_ratio_p, "", 2, 6);
   dataMC_ratio_p->Divide(mc_ratio_p);
   
   TList* returnItems = new TList();
   returnItems->Add(data_ratio_p); returnItems->Add(mc_ratio_p);
   returnItems->Add(dataMC_ratio_p);
   
   return returnItems;
}


//____________________________________________________________________________________________________
void CheckPIDsystematics(const Char_t* dataFilename, const Char_t* mcFilename,
      const Char_t* outputFilename="pidSyst_SingleEle.root") {
   //
   // 
   //
   gHistManData = new AliHistogramManager();
   gHistManData->InitFile(dataFilename, "jpsi2eeHistos");

   gHistManMC = new AliHistogramManager();
   gHistManMC->InitFile(mcFilename, "jpsi2eeHistos");

   const Int_t kNSettings = 9;
   TString settings[kNSettings][2] = {
      "Standard", "electronFullBand",
      "electron0_prot0_pion1", "electronFullBand",
      "electron0_prot0_pion2", "electronFullBand",
      "electron0_prot1_pion0", "electronFullBand",
      "electron0_prot2_pion0", "electronFullBand",
      "electron1_prot0_pion0", "electronFullBand",
      "electron2_prot0_pion0", "electronFullBand", 
      "electron3_prot0_pion0", "electronFullBand", 
      "electron4_prot0_pion0", "electronFullBand"
   };
   
   TList* histLists[kNSettings] = {0x0};
   
   TFile* save=new TFile(outputFilename, "RECREATE");
   TH1F* ptMC = (TH1F*)gHistManMC->GetHistogram("mcTruthJpsi_PureMCTruth_BeforeSelection", "PtMC");
   ptMC->SetName("JpsiPtSpectrum");
   ptMC->Write();
   
   for (Int_t i = 0; i<kNSettings; i++) {
      cout <<  "Check PID systematics for setting " <<  settings[i][0].Data() <<  " with reference " << settings[i][1].Data() <<  endl;
      TH1F* data_p_ref = (TH1F*)gHistManData->GetHistogram(Form("Track_%s", settings[i][1].Data()),"P");
      TH1F* data_p = (TH1F*)gHistManData->GetHistogram(Form("Track_%s", settings[i][0].Data()),"P");
      TH1F* mc_p_ref = (TH1F*)gHistManMC->GetHistogram(Form("Track_%s_TrueElectron", settings[i][1].Data()),"P");
      TH1F* mc_p = (TH1F*)gHistManMC->GetHistogram(Form("Track_%s_TrueElectron", settings[i][0].Data()),"P");
      //histLists[i] = CheckPIDsystematicsForSetting(settings[i][0].Data(), settings[i][1].Data());
      histLists[i] = CheckPIDsystematicsForSetting(data_p_ref, data_p, mc_p_ref, mc_p, settings[i][0].Data());
   }
   
   save->cd();
   for (Int_t i = 0; i<kNSettings; i++) histLists[i]->Write();
   
   save->Close();
}

//____________________________________________________________________________________________________
void MakeSignalExtractionFigure() {
   // 
   // Make the signal extraction performance figure for the paper
   // 
   TCanvas* cv = new TCanvas("cv", "cv", 2200, 1100);
   cv->SetTopMargin(0.0);
   cv->SetRightMargin(0.0);
   cv->SetBottomMargin(0.0);
   cv->SetLeftMargin(0.0);
   cv->Divide(4,2,0.0,0.0);
   
   const Int_t kNBins = 8;
   const Char_t* filenames[kNBins] = {
      "/home/iarsene/work/ALICE/pp2017Analysis/analysisOutputs/20180725_TrackingSystematics/LHC17pq/sigExtrHybrid_final/pt0.00_1.00_pol2/Results.root",
      "/home/iarsene/work/ALICE/pp2017Analysis/analysisOutputs/20180725_TrackingSystematics/LHC17pq/sigExtrHybrid_final/pt1.00_2.00/Results.root",
      "/home/iarsene/work/ALICE/pp2017Analysis/analysisOutputs/20180725_TrackingSystematics/LHC17pq/sigExtrHybrid_final/pt2.00_3.00/Results.root",
      "/home/iarsene/work/ALICE/pp2017Analysis/analysisOutputs/20180725_TrackingSystematics/LHC17pq/sigExtrHybrid_final/pt3.00_4.00/Results.root",
      "/home/iarsene/work/ALICE/pp2017Analysis/analysisOutputs/20180725_TrackingSystematics/LHC17pq/sigExtrHybrid_final/pt4.00_5.00/Results.root",
      "/home/iarsene/work/ALICE/pp2017Analysis/analysisOutputs/20180725_TrackingSystematics/LHC17pq/sigExtrHybrid_final/pt5.00_7.00/Results.root",
      "/home/iarsene/work/ALICE/pp2017Analysis/analysisOutputs/20180725_TrackingSystematics/LHC17pq/sigExtrHybrid_final/pt7.00_10.00/Results.root",
      "/home/iarsene/work/ALICE/pp2017Analysis/analysisOutputs/20180725_TrackingSystematics/LHC17pq/sigExtrHybrid_final/ptIntegrated/Results.root"
   };
   
   for (Int_t iPad = 1; iPad <= kNBins; ++iPad) {
      TVirtualPad* pad = cv->cd(iPad);
      pad->SetTopMargin(0.02);
      pad->SetRightMargin(0.02);
      pad->SetBottomMargin(0.15);
      pad->SetLeftMargin(0.15);
      PlotNiceSigExtraction(filenames[iPad-1], pad, iPad);
   }
}

//____________________________________________________________________________________________________
void PlotNiceSigExtraction(const Char_t* filename, TVirtualPad* pad, Int_t padNumber) {
   // 
   // make a nice sig extraction plot on a provided pad
   //
   pad->cd();
   
   Double_t yRange = 1.0;
   if (padNumber == 1) yRange = 160.;
   if (padNumber == 2) yRange = 180.;
   if (padNumber == 3) yRange = 140.;
   if (padNumber == 4) yRange = 100.;
   if (padNumber == 5) yRange = 70.;
   if (padNumber == 6) yRange = 60.;
   if (padNumber == 7) yRange = 25.;
   if (padNumber == 8) yRange = 700.;
 
   TH2D* range = new TH2D(Form("range%f", gRandom->Rndm()), "", 100, 2.001, 3.99, 100, 0.0, yRange);
   BeutifyTAxis(range->GetXaxis(), "#it{m}_{e^{+}e^{#minus}} (GeV/#it{c}^{2})", kTRUE, 0.05, 1.2, 42, 0.04, 42, 507);
   BeutifyTAxis(range->GetYaxis(), "counts per 40 MeV/#it{c}^{2}", kTRUE, 0.05, 1.2, 42, 0.04, 42, 507);
   range->SetStats(kFALSE);
   range->Draw();
   
   TFile* file = TFile::Open(filename);
   
   TH1D* chi2 = (TH1D*)file->Get("Chi2MC");
   TH1D* counts = (TH1D*)file->Get("RawCounts");
   TH1D* soverb = (TH1D*)file->Get("SoverB");
   
   TString nameSuffix = "Standard_mfit1.20_4.72_mexcl2.50_3.72";
   TH1D* splusb = (TH1D*)file->Get(Form("SignalExtractionQA/SplusB_%s", nameSuffix.Data()));
   BeutifyTH1(splusb, "", 2.0, 1);
   splusb->SetStats(kFALSE);
   splusb->Draw("same");
   
   THStack* hs = new THStack("hs", "");
      
   TH1D* bkg = (TH1D*)file->Get(Form("SignalExtractionQA/Bkg_%s", nameSuffix.Data()));
   BeutifyTH1(bkg, "", 2.0, 7);
   bkg->SetStats(kFALSE);
   bkg->SetFillColor(7);
   bkg->SetFillStyle(3002);
   //bkg->Draw("sameHISTC");
   //hs->Add(bkg);
   
   TF1* bkgRes = (TF1*)file->Get(Form("SignalExtractionQA/BkgFunc_%s", nameSuffix.Data()));
   TH1D* bkgTotal = (TH1D*)bkg->Clone(Form("%s_bkgTotal", bkg->GetName()));
   bkgTotal->Add(bkgRes);
   BeutifyTH1(bkgTotal, "", 2.0, 4);
   bkgTotal->SetFillColor(4);
   bkgTotal->SetFillStyle(3002);
   
   
   TH1D* mcLine = (TH1D*)file->Get(Form("SignalExtractionQA/SignalMC_%s", nameSuffix.Data()));
   mcLine->Scale(counts->GetBinContent(1)/mcLine->Integral(mcLine->FindBin(2.9201), mcLine->FindBin(3.159)));
   TH1D* totalFit = (TH1D*)bkgTotal->Clone(Form("%s_totalFit", bkgTotal->GetName()));
   totalFit->Add(mcLine);
   BeutifyTH1(totalFit, "", 2.0, 2);
   totalFit->SetStats(kFALSE);
   totalFit->SetFillColor(2);
   totalFit->SetFillStyle(3002);
 
   totalFit->Draw("sameHISTC");
   bkgTotal->Draw("sameHISTC");
   bkg->Draw("sameHISTC");
   
      
   TLatex latex;
   latex.SetNDC();
   latex.SetTextFont(42);
   latex.SetTextSize(0.04);
   Double_t ptRanges[8][2] = {
      {0.0, 1.0}, {1.0, 2.0}, {2.0, 3.0}, {3.0, 4.0}, {4.0, 5.0}, {5.0, 7.0}, {7.0, 10.0}, {0.0, 15.0}
    };
   if(padNumber == 1) {
      latex.DrawLatex(0.20, 0.93, "ALICE pp #sqrt{#it{s}}=5.02 TeV, #it{L}_{int} = 19.3 nb^{-1} #pm 5%");
      latex.DrawLatex(0.20, 0.88, "J/#psi #rightarrow e^{+}e^{#minus}, |#it{y}|<0.9");
   }
   if (padNumber>1) {
      latex.DrawLatex(0.61, 0.94, Form("%.1f < #it{p}_{T} < %.1f GeV/#it{c}", ptRanges[padNumber-1][0], ptRanges[padNumber-1][1]));
      latex.DrawLatex(0.61, 0.89, Form("#it{#chi}^{2}/NDF = %.2f", chi2->GetBinContent(1)));
      latex.DrawLatex(0.61, 0.84, Form("#it{N}_{J/#psi} = %.0f #pm %.0f", counts->GetBinContent(1), counts->GetBinError(1)));
   }
   if (padNumber == 1) {
      latex.DrawLatex(0.20, 0.74, Form("%.1f < #it{p}_{T} < %.1f GeV/#it{c}", ptRanges[padNumber-1][0], ptRanges[padNumber-1][1]));
      latex.DrawLatex(0.20, 0.69, Form("#it{#chi}^{2}/NDF = %.2f", chi2->GetBinContent(1)));
      latex.DrawLatex(0.20, 0.64, Form("#it{N}_{J/#psi} = %.0f #pm %.0f", counts->GetBinContent(1), counts->GetBinError(1)));
   }
   //latex.DrawLatex(0.20, 0.70, Form("#it{S}/#it{B} = %.1f #pm %.1f", soverb->GetBinContent(1), soverb->GetBinError(1)));
   
   TLegend* legend = 0x0;
   if (padNumber == 7) {
      legend = new TLegend(0.17, 0.75, 0.5, 0.97);
      legend->SetFillColor(0);
      legend->SetBorderSize(0.);
      legend->SetTextFont(42);
      legend->SetTextSize(0.04);
      legend->AddEntry(totalFit, "Signal", "f");
      legend->AddEntry(bkgTotal, "Correlated bkg.", "f");
      legend->AddEntry(bkg, "Combinatorial bkg.", "f");
      legend->Draw();
   }
   
 //hs->Draw("samePLC");
   //file->Close();
}


//____________________________________________________________________________________________________
void ComputeNonPromptFractionSystematics(Double_t ptMin, Double_t ptMax) {
   // 
   // 
   // 
   
   // parameterized fB vs pt from MC
   TF1* fbMC = new TF1("fbMC", "pol5", 0.0, 10.0);
   fbMC->SetParameters(9.49765e-2, 1.83611e-2, -7.87568e-3, 4.92532e-3, -6.40202e-4, 2.58012e-5);
   
   // parameterized fB vs pt tuned to existing measurements
   TF1* fbData = new TF1("fbData", "pol3", 0.0, 10.0);
   fbData->SetParameters(0.0949766, 0.001, 0.001, 0.0001);
   
   TFile* resultsFileInclusive = TFile::Open("/home/iarsene/work/ALICE/pp2017Analysis/analysisOutputs/20180725_TrackingSystematics/LHC17pq/sigExtrHybrid_final/ptIntegrated/Results.root");
   TFile* resultsFilePrompt = TFile::Open("/home/iarsene/work/ALICE/pp2017Analysis/analysisOutputs/20180912_StandardCutPolarizChecks/sigExtrHybridPrompt/Results.root");
   TFile* resultsFileNonPrompt = TFile::Open("/home/iarsene/work/ALICE/pp2017Analysis/analysisOutputs/20180912_StandardCutPolarizChecks/sigExtrHybridNonPrompt/Results.root");
   TH1D* effVsPtInclusive = resultsFileInclusive->Get("Efficiency_vs_Pt/effVsPtAverage_Standard_TrueElectron");
   TH1D* effVsPtPrompt = resultsFilePrompt->Get("Efficiency_vs_Pt/effVsPtAverage_Standard_TrueElectronPrompt");
   TH1D* effVsPtNonPrompt = resultsFileNonPrompt->Get("Efficiency_vs_Pt/effVsPtAverage_Standard_TrueElectronNonPrompt");
   
   effVsPtInclusive->SetLineColor(1); effVsPtInclusive->SetLineWidth(2);
   effVsPtPrompt->SetLineColor(4); effVsPtPrompt->SetLineWidth(2);
   effVsPtNonPrompt->SetLineColor(2); effVsPtNonPrompt->SetLineWidth(2);
   
   // Get the corrected spectrum as measured at mid-y in pp at 5 TeV (step 1 spectrum)
   TFile* spectraFile = TFile::Open("/home/iarsene/work/ALICE/pp2017Analysis/SpectraPP5TeV.root");
   TGraphErrors* spectra = (TGraphErrors*)spectraFile->Get("StatUncert");

   // fit the measured spectra with a power law function
   TF1* fplSpectra = new TF1("fplSpectra","[0]*x/TMath::Power((1+TMath::Power(x/[1],2.0)),[2])",0.,10.);
   fplSpectra->SetParameters(1.0, 3.67987, 3.19625);
   fplSpectra->SetLineColor(4);
   spectra->Fit(fplSpectra, "", "QMEI", 0.0, 10.0);
   
   return;
   
   Double_t err = 0.0;
   Double_t effInclusive = WeightedAverage(effVsPtInclusive, fplSpectra, ptMin, ptMax, err);
   Double_t effPrompt = WeightedAverage(effVsPtPrompt, fplSpectra, ptMin, ptMax, err);
   Double_t effNonPrompt = WeightedAverage(effVsPtNonPrompt, fplSpectra, ptMin, ptMax, err);
   cout << "eff inclusive / prompt / nonprompt :: " << effInclusive << " / " << effPrompt << " / " << effNonPrompt << endl;
   
   Double_t ptLims[7] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 10.0};
   TH1D* correctionHist = new TH1D("correction", "", 6, ptLims);
   for(Int_t i = 0; i<6; ++i) {
      correctionHist->SetBinContent(i+1, 100.0*WeightedAverageNonPrompt(effVsPtInclusive, effVsPtPrompt, effVsPtNonPrompt, fplSpectra, fbData, 
                     correctionHist->GetXaxis()->GetBinLowEdge(i+1), correctionHist->GetXaxis()->GetBinUpEdge(i+1), err));
   }
   
   Double_t crossSectionShift = WeightedAverageNonPrompt(effVsPtInclusive, effVsPtPrompt, effVsPtNonPrompt, fplSpectra, fbData, ptMin, ptMax, err);
   cout << "shift :: " << 100*crossSectionShift << " %" << endl;
   correctionHist->Draw();
   
   //spectra->Draw("A*"); 
   //return;
   TFile* saveFile = new TFile("NonPromptFractionShift.root", "RECREATE");
   correctionHist->Write();
   saveFile->Close();
   //fbMC->Draw();
   //fbData->Draw("same");
}

//____________________________________________________________________________
Double_t WeightedAverageNonPrompt(TH1D* eff, TH1D* promptEff, TH1D* nonpromptEff, TF1* weight, TF1* fbData, Float_t ptMin, Float_t ptMax, Double_t& err) {
   //
   // Takes as arguments the efficiency using inclusive jpsi,  prompt and non-prompt as a function of pt, the fb distribution as a function of pt
   // from data and MC and the input pt spectrum.
   // Calculates the shift to the cross-section needed for the fB fraction observed in data
   //
   Double_t stdAverage = 0.0;
   Double_t norm = 0.0;
   Double_t corrAverage = 0.0;
   err = 0.0;
   for(Int_t i=1; i<=eff->GetXaxis()->GetNbins(); ++i) {
      Double_t pt=eff->GetXaxis()->GetBinCenter(i);
      if(ptMin>0.0 && pt<ptMin) continue;
      if(ptMax>0.0 && pt>ptMax) continue;
      Double_t stdEff = eff->GetBinContent(i);
      Double_t fb = fbData->Eval(pt);
      Double_t correctEff = fb*nonpromptEff->GetBinContent(i) + (1.0-fb)*promptEff->GetBinContent(i); 
      stdAverage += stdEff * weight->Eval(pt) * eff->GetXaxis()->GetBinWidth(i);
      err += TMath::Power(eff->GetBinError(i) * weight->Eval(pt) * eff->GetXaxis()->GetBinWidth(i), 2.0);
      norm += weight->Eval(pt) * eff->GetXaxis()->GetBinWidth(i);
      corrAverage += correctEff * weight->Eval(pt) * eff->GetXaxis()->GetBinWidth(i);
      //corrHist->SetBinContent(i, 100.0*2.0*(1.0/stdEff-1.0/correctEff)/(1.0/stdEff+1.0/correctEff)); 
   }
   stdAverage /= norm;
   corrAverage /= norm;
   err = TMath::Sqrt(err); 
   err /= norm;
   return 2.0*(1.0/corrAverage-1.0/stdAverage)/(1.0/stdAverage+1.0/corrAverage);
}


//____________________________________________________________________________________________________
void ComputePolarizationSystematics(Double_t ptMin, Double_t ptMax) {
   // 
   // compute syst uncertainty due to unknown jpsi polarization
   //
   
   //  Construct the efficiency as a function of cos theta * and pt
   AliHistogramManager* man = new AliHistogramManager();
   man->InitFile("/home/iarsene/work/ALICE/pp2017Analysis/analysisOutputs/20180912_StandardCutPolarizChecks_MC/LHC17pq/AnalysisHistograms_jpsi2ee_pp2017_MC.root", "jpsi2eeHistos");
   
   TH2F* denomCS = (TH2F*)man->GetHistogram("mcTruthJpsi_PureMCTruth_BeforeSelection", "CosThetaStarCS_ptMC");
   TH2F* denomHE = (TH2F*)man->GetHistogram("mcTruthJpsi_PureMCTruth_BeforeSelection", "CosThetaStarHE_ptMC");
   TH3F* nomCS = (TH3F*)man->GetHistogram("PairSEPM_Standard_TrueElectron", "CosThetaStarCS_pt_mass");
   TH3F* nomHE = (TH3F*)man->GetHistogram("PairSEPM_Standard_TrueElectron", "CosThetaStarHE_pt_mass");
   
   nomCS->GetZaxis()->SetRangeUser(2.921, 3.159);
   TH2D* projNomCS = (TH2D*)nomCS->Project3D("yx");
   projNomCS->SetName("projNomCS");
   nomHE->GetZaxis()->SetRangeUser(2.921, 3.159);
   TH2D* projNomHE = (TH2D*)nomHE->Project3D("yx");
   projNomHE->SetName("projNomHE");
   
   TH2D* effCS = (TH2D*)projNomCS->Clone("effCS");
   effCS->Divide(denomCS);
   TH2D* effHE = (TH2D*)projNomHE->Clone("effHE");
   effHE->Divide(denomHE);
   //denomHE->Draw("colz");
   //return;
   //effHE->Draw("colz");
   //return;
   
   // Get the corrected spectrum as measured at mid-y in pp at 5 TeV (step 1 spectrum)
   TFile* spectraFile = TFile::Open("/home/iarsene/work/ALICE/pp2017Analysis/SpectraPP5TeV.root");
   TGraphErrors* spectra = (TGraphErrors*)spectraFile->Get("StatUncert");

   //effVsPt->Draw();
   // fit the measured spectra with a similar function
   TF1* fplSpectra = new TF1("fplSpectra","[0]*x/TMath::Power((1+TMath::Power(x/[1],2.0)),[2])",0.,10.);
   fplSpectra->SetParameters(1.0, 3.67987, 3.19625);
   fplSpectra->SetLineColor(4);
   spectra->Fit(fplSpectra, "", "QMEI", 0.0, 10.0);
   
   Double_t meanPar1 = fplSpectra->GetParameter(1);
   Double_t par1Err = fplSpectra->GetParError(1);
   Double_t meanPar2 = fplSpectra->GetParameter(2);
   Double_t par2Err = fplSpectra->GetParError(2);
   cout << "par1  :: " << meanPar1 << " +/- " << par1Err << endl;
   cout << "par2  :: " << meanPar2 << " +/- " << par2Err << endl;
   cout << "chi2 / ndf :: " << fplSpectra->GetChisquare() << " / " << fplSpectra->GetNDF() << " = " << fplSpectra->GetChisquare() / fplSpectra->GetNDF() << endl;
   
   
   // construct the polarization profile
   TF1* polariz = new TF1("polariz", "1+[0]*x*x", -1.0, 1.0);
   TGraph* grPolarizCS = new TGraph();
   TGraph* grPolarizHE = new TGraph();
   Int_t currentPoint = 0;
   Double_t effLong = 0.0;
   Double_t effTrans = 0.0;
   Double_t effUnpolariz = 0.0;
   for(Double_t alpha = -1.0; alpha <= 1.01; alpha += 0.05) {
      Double_t eff = WeightedAverage(effCS, fplSpectra, alpha, ptMin, ptMax, kTRUE);
      grPolarizCS->SetPoint(currentPoint, alpha, eff);
      eff = WeightedAverage(effHE, fplSpectra, alpha, ptMin, ptMax, kTRUE);
      if(TMath::Abs(alpha+1.0) <0.001) effLong = 1./eff;
      if(TMath::Abs(alpha-1.0) <0.001) effTrans = 1./eff;
      if(TMath::Abs(alpha-0.0) <0.001) effUnpolariz = 1./eff;
      grPolarizHE->SetPoint(currentPoint, alpha, eff-0.00109);
      //grPolarizHE->SetPoint(currentPoint, alpha, eff);
      ++currentPoint;
   }
   cout << "effUnpolariz, effLong, effTrans :: " << effUnpolariz << " / " << effLong << " / " << effTrans << endl;
   cout << "longitudinal :: " << 100.0*2.0*(effLong-effUnpolariz) / (effLong+effUnpolariz) << " %" << endl;
   cout << "transversal :: " << 100.0*2.0*(effTrans-effUnpolariz) / (effTrans+effUnpolariz) << " %" << endl;
   
   
   grPolarizCS->SetName("Collins-Soper");
   grPolarizHE->SetName("Helicity");
   grPolarizCS->SetLineColor(2); grPolarizCS->SetLineWidth(2);
   grPolarizHE->SetLineColor(4); grPolarizHE->SetLineWidth(2);
   grPolarizCS->Draw("AL");
   grPolarizHE->Draw("sameL");
   return;
   
   
   TFile* resultsFile = TFile::Open("/home/iarsene/work/ALICE/pp2017Analysis/analysisOutputs/20180725_TrackingSystematics/LHC17pq/sigExtrHybrid_final/ptIntegrated/Results.root");
   TH1D* effVsPt = resultsFile->Get("Efficiency_vs_Pt/effVsPtAverage_Standard_TrueElectron");
   
   effVsPt->Draw();
   
   TGraph* grEffPtCSlong = new TGraph(); grEffPtCSlong->SetName("grEffPtCSlong");
   TGraph* grEffPtCStran = new TGraph(); grEffPtCStran->SetName("grEffPtCStran");
   TGraph* grEffPtHElong = new TGraph(); grEffPtHElong->SetName("grEffPtHElong");
   TGraph* grEffPtHEtran = new TGraph(); grEffPtHEtran->SetName("grEffPtHEtran");
   
   Int_t currentPoint = 0;
   for(Int_t i=1; i <= effCS->GetYaxis()->GetNbins(); ++i) {
      Double_t pt = effCS->GetYaxis()->GetBinCenter(i);
      Double_t ptLow = effCS->GetYaxis()->GetBinLowEdge(i);
      Double_t ptHigh = effCS->GetYaxis()->GetBinUpEdge(i);
      if(pt < ptMin) continue;
      if(pt > ptMax) continue;
      //polariz->SetParameter(0, -1.0);
      Double_t eff = WeightedAverage(effCS, fplSpectra, -1.0, ptLow, ptHigh);
      grEffPtCSlong->SetPoint(currentPoint, pt, eff);
      //polariz->SetParameter(0, +1.0);
      eff = WeightedAverage(effCS, fplSpectra, +1.0, ptLow, ptHigh);
      grEffPtCStran->SetPoint(currentPoint, pt, eff);
      //polariz->SetParameter(0, -1.0);
      eff = WeightedAverage(effHE, fplSpectra, -1.0, ptLow, ptHigh);
      grEffPtHElong->SetPoint(currentPoint, pt, eff);
      //polariz->SetParameter(0, +1.0);
      eff = WeightedAverage(effHE, fplSpectra, +1.0, ptLow, ptHigh);
      grEffPtHEtran->SetPoint(currentPoint, pt, eff);
      currentPoint++;
   }
   
   grEffPtCSlong->SetLineStyle(2); grEffPtCStran->SetLineStyle(2);
   grEffPtHElong->SetLineStyle(2); grEffPtHEtran->SetLineStyle(2);
   grEffPtCSlong->SetLineColor(2); grEffPtCStran->SetLineColor(2);
   grEffPtHElong->SetLineColor(4); grEffPtHEtran->SetLineColor(4);
   grEffPtCSlong->SetLineWidth(2); grEffPtCStran->SetLineWidth(2);
   grEffPtHElong->SetLineWidth(2); grEffPtHEtran->SetLineWidth(2);
   grEffPtCSlong->Draw("sameL");
   grEffPtCStran->Draw("sameL");
   grEffPtHElong->Draw("sameL");
   grEffPtHEtran->Draw("sameL");
}


//____________________________________________________________________________________________________
void ComputeInputKineSyst(Double_t ptMin, Double_t ptMax) {
   //
   // computes the systematic uncertainty on the choice of the input MC kinematics
   //
   
   // Get the efficiency vs pt
   TFile* resultsFile = TFile::Open("/home/iarsene/work/ALICE/pp2017Analysis/analysisOutputs/20180725_TrackingSystematics/LHC17pq/sigExtrHybrid_final/ptIntegrated/Results.root");
   //TFile* resultsFile = TFile::Open("/home/iarsene/work/ALICE/pp2017Analysis/analysisOutputs/20180912_StandardCutPolarizChecks/sigExtrHybridPrompt/Results.root");
   //TFile* resultsFile = TFile::Open("/home/iarsene/work/ALICE/pp2017Analysis/analysisOutputs/20180912_StandardCutPolarizChecks/sigExtrHybridNonPrompt/Results.root");
   TH1D* effVsPt = resultsFile->Get("Efficiency_vs_Pt/effVsPtAverage_Standard_TrueElectron");
   //TH1D* effVsPt = resultsFile->Get("Efficiency_vs_Pt/effVsPtAverage_Standard_TrueElectronPrompt");
   //TH1D* effVsPt = resultsFile->Get("Efficiency_vs_Pt/effVsPtAverage_Standard_TrueElectronNonPrompt");
   
   // Function which parametrizes the input MC jpsi distribution used in the LHC17pq anchored simulations
   TF1* fpl= new TF1("fpl","[0]*x/TMath::Power((1+TMath::Power(x/[1],[3])),[2])",0.,30.);
   fpl->SetParameters(1800.0, 3.67987, 3.19625, 2.0);          //  these are the parameters of the fit to the 2017 MC pt distribution
   
   // Get the corrected spectrum as measured at mid-y in pp at 5 TeV (step 1 spectrum)
   TFile* spectraFile = TFile::Open("/home/iarsene/work/ALICE/pp2017Analysis/SpectraPP5TeV.root");
   TGraphErrors* spectra = (TGraphErrors*)spectraFile->Get("StatUncert");
   
   //effVsPt->Draw();
   // fit the measured spectra with a similar function
   TF1* fplSpectra = new TF1("fplSpectra","[0]*x/TMath::Power((1+TMath::Power(x/[1],2.0)),[2])",0.,10.);
   fplSpectra->SetParameters(1.0, 3.67987, 3.19625);
   fplSpectra->SetLineColor(4);
   spectra->Fit(fplSpectra, "", "QMEI", 0.0, 10.0);
   
//return;   
   /*gMinuit->SetErrorDef(1.0);                              // 1 sigmas (2^2)
   TGraph* g12 = (TGraph*)(gMinuit->Contour(200,1,2));
   TGraph* g12Clone = (TGraph*)g12->Clone("g12Clone");
   //g12->Draw("al");
     */
   
   Double_t err = 0.0;
   Double_t effStandard = WeightedAverage(effVsPt, fpl, ptMin, ptMax, err);
   TH1F* hEffStandard = new TH1F("hEffStandard", "", 1, 0., 1.);
   hEffStandard->SetBinContent(1, effStandard);
   cout << "effStandard = " << effStandard << endl;
   Double_t effSpectra = WeightedAverage(effVsPt, fplSpectra, ptMin, ptMax, err);
   TH1F* hEffSpectra = new TH1F("hEffSpectra", "", 1, 0., 1.);
   hEffSpectra->SetBinContent(1, effSpectra);
   cout << "effSpectra = " << effSpectra << endl;
   
   //TCanvas* c1=new TCanvas();   

   Double_t meanPar1 = fplSpectra->GetParameter(1);
   Double_t par1Err = fplSpectra->GetParError(1);
   Double_t meanPar2 = fplSpectra->GetParameter(2);
   Double_t par2Err = fplSpectra->GetParError(2);
   cout << "par1  :: " << meanPar1 << " +/- " << par1Err << endl;
   cout << "par2  :: " << meanPar2 << " +/- " << par2Err << endl;
   cout << "chi2 / ndf :: " << fplSpectra->GetChisquare() << " / " << fplSpectra->GetNDF() << " = " << fplSpectra->GetChisquare() / fplSpectra->GetNDF() << endl;
   
   TF1* parCorrelation = new TF1("parCorrelation", "pol2", 0.0, 10.0);
   parCorrelation->SetParameters(1.06, 0.2, 0.092);
   //parCorrelation->Draw("same");
   //return;
   
   TH2D* fitSurf = new TH2D("fitSurf", "", 100, 0.0, 7.5, 100, 0.0, 7.5);
   TH2D* fitSurface = GetFitSurface(fitSurf, 0.0, TMath::Sqrt(2.3), 200);
   //fitSurface->Draw("colz");
   fitSurface->GetXaxis()->SetTitle("pt0");
   fitSurface->GetYaxis()->SetTitle("n");
   
   TH1D* hEff = new TH1D("hEff", "efficiency after MC input kine variations", 1000, 0.02, 0.16);
   TH1D* hInvEff = new TH1D("hInvEff", "1.0/efficiency after MC input kine variations", 1000, 5.0, 20.0);
   TH2D* paramsSampled = new TH2D("paramsSampled", "", 100, 0.0, 7.5, 100, 0.0, 7.5);
   
   fplSpectra->SetLineWidth(0.5);
   fplSpectra->SetLineColor(2);
   fplSpectra->SetLineStyle(2);
   TCanvas* cv = new TCanvas();
   fplSpectra->Draw();
   for(Int_t i=0; i<1000;++i) {
      Double_t pt0 = 0.0;
      Double_t n = 0.0;
      fitSurface->GetRandom2(pt0, n);
      //pt0 = gRandom->Gaus(meanPar1, par1Err);
      //n = gRandom->Gaus(parCorrelation->Eval(pt0), 0.25);
      
      fplSpectra->SetParameter(1, pt0);
      fplSpectra->SetParameter(2, n);
      paramsSampled->Fill(pt0, n);
      
      Double_t err = 0.0;
      Double_t eff = WeightedAverage(effVsPt, fplSpectra, ptMin, ptMax, err);
      hEff->Fill(eff);
      hInvEff->Fill(1.0/eff);
      fplSpectra->DrawClone("same");
   }
   spectra->Draw("sameE");
   fplSpectra->SetParameter(1, meanPar1);
   fplSpectra->SetParameter(2, meanPar2);
   fplSpectra->SetLineWidth(3);
   fplSpectra->SetLineColor(1); fplSpectra->SetLineStyle(1);
   fplSpectra->Draw("same");
   fpl->SetLineColor(4);
   fpl->Draw("same");
   
      
   TCanvas* c2=new TCanvas();
   hEff->Draw();
   TCanvas* c3=new TCanvas();
   hInvEff->Draw();
   TCanvas* c4=new TCanvas();
   paramsSampled->Draw("colz");
   TCanvas* c5=new TCanvas();
   fitSurface->Draw("colz");
   spectra->Fit(fplSpectra, "", "QMEI", 0.0, 10.0);
   gMinuit->SetErrorDef(1.0);                               // 1 sigmas (2^2)
   TGraph* g12 = (TGraph*)(gMinuit->Contour(200,1,2));
   g12->SetLineWidth(2);
   g12->Draw("sameL");
   TLine line;
   line.SetLineWidth(2); line.SetLineStyle(1);
   line.DrawLine(meanPar1, fitSurface->GetYaxis()->GetXmin(), meanPar1, fitSurface->GetYaxis()->GetXmax());
   line.DrawLine(fitSurface->GetXaxis()->GetXmin(), meanPar2, fitSurface->GetXaxis()->GetXmax(), meanPar2);
   line.SetLineStyle(2.0);
   line.DrawLine(meanPar1-par1Err, meanPar2-par2Err, meanPar1-par1Err, meanPar2+par2Err);
   line.DrawLine(meanPar1+par1Err, meanPar2-par2Err, meanPar1+par1Err, meanPar2+par2Err);
   line.DrawLine(meanPar1-par1Err, meanPar2-par2Err, meanPar1+par1Err, meanPar2-par2Err);
   line.DrawLine(meanPar1-par1Err, meanPar2+par2Err, meanPar1+par1Err, meanPar2+par2Err);
   
   TFile* save = new TFile(Form("MCinputKineStats_pt%.2f_%.2f.root", ptMin, ptMax), "RECREATE");
   fplSpectra->Write();
   hEff->Write();
   hInvEff->Write();
   paramsSampled->Write();
   fitSurface->Write();
   g12->Write();
   hEffSpectra->Write();
   hEffStandard->Write();
   save->Close();
}

// ______________________________________________________________________________________________________________________
TH2D* GetFitSurface(TH2D* h, Double_t sigmaMin, Double_t sigmaMax, Int_t n) {
   
   // build a surface with fit param probabilities
   
   TF1* gausFunc = new TF1("gausFunc", "gaus", -10.0, 10.0);
   gausFunc->SetParameters(1.0, 0.0, 1.0);
   
   for(Int_t i=0; i<n; ++i) {
      //cout << "contour for " << sigmaMax - i*(sigmaMax-sigmaMin)/Double_t(n) << endl;
      gMinuit->SetErrorDef(TMath::Power(sigmaMax - i*(sigmaMax-sigmaMin)/Double_t(n), 2.0));
      Double_t weight = gausFunc->Eval(sigmaMax - i*(sigmaMax-sigmaMin)/Double_t(n));
      TGraph* g12 = (TGraph*)gMinuit->Contour(150,1,2);
      for(Int_t j = 0;j<g12->GetN();++j) {
         h->SetBinContent(h->FindBin(g12->GetX()[j], g12->GetY()[j]), weight);
      }
      delete g12;
   }
   
   TH2D* hClone = (TH2D*)h->Clone("clone");
   
   TF1* fitFunc = new TF1("fitFunc", "gaus", 2.0, 6.0);
   for(Int_t x=1; x<h->GetXaxis()->GetNbins(); ++x) {
      TH1D* proj = h->ProjectionY(Form("projY_%d", x), x, x);
      Int_t nFilledBins = GetFilledBins(proj);
      //cout << "projX " << x << "  entries " << nFilledBins <<  endl;
      if(nFilledBins<7) continue;
      proj->Fit(fitFunc, "Q", "QME", 2.0, 6.0);
      //cout << "amp / mean / width  " << fitFunc->GetParameter(0) << " / " << fitFunc->GetParameter(1) << " / " << fitFunc->GetParameter(2) << endl; 
      
      for(Int_t y=1; y<h->GetYaxis()->GetNbins();++y) {
         if(h->GetBinContent(x,y)>0.0) continue;
         hClone->SetBinContent(x,y, fitFunc->Eval(h->GetYaxis()->GetBinCenter(y)));
      }
      delete proj;
   }
   for(Int_t y=1; y<h->GetYaxis()->GetNbins(); ++y) {
      TH1D* proj = h->ProjectionX(Form("projX_%d", y), y, y);
      Int_t nFilledBins = GetFilledBins(proj);
      //cout << "projY " << y << "  entries " << nFilledBins <<  endl;
      if(nFilledBins<7) continue;
      proj->Fit(fitFunc, "Q", "QME", 2.0, 6.0);
      //cout << "amp / mean / width  " << fitFunc->GetParameter(0) << " / " << fitFunc->GetParameter(1) << " / " << fitFunc->GetParameter(2) << endl; 

      for(Int_t x=1; x<h->GetXaxis()->GetNbins();++x) {
         if(hClone->GetBinContent(x,y)>0.0) continue;
         hClone->SetBinContent(x,y, fitFunc->Eval(h->GetXaxis()->GetBinCenter(x)));
      }
      delete proj;
   }
   /*
   for(Int_t x=1; x<h->GetXaxis()->GetNbins(); ++x) {
      TH1D* proj = hClone->ProjectionY(Form("projY_%d", x), x, x);
      Int_t nFilledBins = GetFilledBins(proj);
      cout << "projX " << x << "  entries " << nFilledBins <<  endl;
      if(nFilledBins<10) continue;
      proj->Fit(fitFunc, "Q", "QME", 2.0, 6.0);
      cout << "amp / mean / width  " << fitFunc->GetParameter(0) << " / " << fitFunc->GetParameter(1) << " / " << fitFunc->GetParameter(2) << endl; 

      for(Int_t y=1; y<h->GetYaxis()->GetNbins();++y) {
         if(hClone->GetBinContent(x,y)>0.0) continue;
         hClone->SetBinContent(x,y, fitFunc->Eval(h->GetYaxis()->GetBinCenter(y)));
       }
    } */
   
   return hClone;
}

// _______________________________________________________________________________________________
Int_t GetFilledBins(TH1D* h) {
   
   Int_t nFilledBins = 0;
   for (Int_t i = 1; i<h->GetXaxis()->GetNbins();++i )
      if (h->GetBinContent(i)>0.0) nFilledBins++;
         
   return nFilledBins;
}

//________________________________________________________________________________________
void CheckMtScalingWithDmeson() {
   
   TF1* fplSpectra = new TF1("fplSpectra","[0]*x/TMath::Power((1+TMath::Power(x/[1],2.0)),[2])",0.,10.);
   fplSpectra->SetParameters(1871., 3.244, 2.6941);
   fplSpectra->SetLineColor(4);
   
  /* Double_t mt = TMath::Sqrt(pt*pt+gkMassJpsi*gkMassJpsi);
   mt -= gkMassJpsi;     // mt - m0
   mt += mass;
   pt = TMath::Sqrt(mt*mt-mass*mass);
*/
   //y = TMath::Sqrt((TMath::Sqrt(x*x+3.1*3.1)-3.1+1.87)*(TMath::Sqrt(x*x+3.1*3.1)-3.1+1.87)-1.87*1.87);   
   
   TF1* mtScaled = new TF1("mtScaled","[0]*TMath::Sqrt((TMath::Sqrt(x*x+3.1*3.1)-3.1+1.87)*(TMath::Sqrt(x*x+3.1*3.1)-3.1+1.87)-1.87*1.87)/TMath::Power((1+TMath::Power(TMath::Sqrt((TMath::Sqrt(x*x+3.1*3.1)-3.1+1.87)*(TMath::Sqrt(x*x+3.1*3.1)-3.1+1.87)-1.87*1.87)/[1],2.0)),[2])",0.,10.);
   mtScaled->SetParameters(1871., 3.244, 2.6941);
   mtScaled->SetLineColor(2);
   
   TF1* ratio = new TF1("ratio","fplSpectra/mtScaled",0.0,10.0);
   ratio->SetParameters(1871.*0.050, 3.244, 2.6941, 1871., 3.244, 2.6941);
   ratio->SetLineColor(6);
   
   ratio->Draw();
   
   TCanvas* c2=new TCanvas();
   
   fplSpectra->Draw();
   mtScaled->Draw("same");
}
