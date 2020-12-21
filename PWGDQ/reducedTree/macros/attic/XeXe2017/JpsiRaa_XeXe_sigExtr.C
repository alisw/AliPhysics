// TAA values for Xe-Xe:       TAA, rms, sys error
//  centrality ranges, respectively, 0-20%, 20-90%, 0-90%
Double_t gTAA_XeXe[3][3] = {              
   {9.896, 2.9, 0.88},
   {1.352, 1.5, 0.2},
   {3.252, 4.0, 0.35}
};

// Ncoll values for Xe-Xe:       ncoll, rms, sys error
//  centrality ranges, respectively, 0-20%, 20-90%, 0-90%
Double_t gNcoll_XeXe[3][3] = {              
   {676.9, 200., 60.},
   {92.51, 110., 14.},
   {222.4, 280., 24.}
};

// Npart values for Xe-Xe:       npart, rms, sys error
//  centrality ranges, respectively, 0-20%, 20-90%, 0-90%
Double_t gNpart_XeXe[3][3] = {              
   {193., 34., 2.5},
   {46.42, 39., 2.5},
   {79.01, 72., 2.5}
};

// J/psi cross-section at 5.02 TeV obtained via interpolation at mid-y
// B.R. x d sigma/ dy 
//Double_t gPPcrossSection5TeVinterpol = 367.0e-9;        // nb

// J/psi cross-section at 5.44 TeV obtained via interpolation at mid-y
// B.R. x d sigma/ dy 
Double_t gPPcrossSection5TeVinterpol = 388.6e-9;        // nb, excluding ln+powerlaw

TH1* gOldPointer = 0x0;

AliHistogramManager* gHistManData = 0x0;
AliHistogramManager* gHistManMC = 0x0;
AliResonanceFits gFits;

//________________________________________________________________________________________
void RunSystematics(const Char_t* filename, const Char_t* filenameMC, Int_t centBin=2, const Char_t* outputDir=".", Bool_t savePictures=kFALSE) {
   //
   // Run systematics
   //
   gHistManData = new AliHistogramManager();
   gHistManData->InitFile(filename, "jpsi2eeHistos");
   
   gHistManMC = new AliHistogramManager();
   gHistManMC->InitFile(filenameMC, "jpsi2eeHistos");
   
   TH1F* eventsVsCent = (TH1F*)gHistManData->GetHistogram("Event_AfterCuts","CentVZERO");
   
   //gFits.AddVariable(AliReducedVarManager::kVtxZ, 3);
   gFits.AddVariable(AliReducedVarManager::kCentVZERO, 2);
   gFits.AddVariable(AliReducedVarManager::kPt, 1);
   gFits.AddVariable(AliReducedVarManager::kMass, 0);
   gFits.SetVarRange(AliReducedVarManager::kMass, 1.5, 4.2);
   
   gFits.SetBkgMethod(AliResonanceFits::kBkgMixedEvent);          // kBkgMixedEvent,  kBkgLikeSign, kBkgFitFunction, kBkgMixedEventAndResidualFit
   gFits.SetMEMatchingMethod(AliResonanceFits::kMatchSEOS);     // kMatchSEOS, kMatchSELS
   gFits.SetLSmethod(AliResonanceFits::kLSArithmeticMean);       // kLSArithmeticMean, kLSGeometricMean
   gFits.SetUseRfactorCorrection(kTRUE);
   gFits.SetMassFitRange(1.5,4.2);
   gFits.SetMassExclusionRange(2.5,4.2);
   gFits.SetScaleSummedBkg(kTRUE);
   //gFits.SetUseSignificantZero(kFALSE);
   gFits.SetScalingOption(AliResonanceFits::kScaleEntries);   // kScaleEntries, kScaleFit, kScaleWeightedAverage
   gFits.SetWeightedAveragePower(2.0);                                 // 2 - is default, but it can be modified in some specific situations
   gFits.SetMinuitFitOption(AliResonanceFits::kMinuitMethodLikelihood);   // kMinuitMethodChi2, kMinuitMethodLikelihood
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
   
   /*
   const Int_t kNsettings = 16;
   const Char_t* settingNames[kNsettings] = {
      "prot36_pion34_spdAny", "prot36_pion34_spdFirst", "prot36_pion34_itsCls3", "prot36_pion34_itsCls4",
      "prot36_pion34_itsChi2_19",  "prot36_pion34_itsChi2_17", "prot36_pion34_itsChi2_13", 
      //"prot36_pion34_itsChi2_11", "prot36_pion34_itsChi2_9", "prot36_pion34_itsChi2_7",
      "prot36_pion34_tpcChi2_60", "prot36_pion34_tpcChi2_50", "prot36_pion34_tpcChi2_35", 
      //"prot36_pion34_tpcChi2_30", "prot36_pion34_tpcChi2_25",
      "prot36_pion34_tpcSegm4", "prot36_pion34_tpcSegm5", "prot36_pion34_tpcSegm7",
      "prot36_pion34_tpcCls90", "prot36_pion34_tpcCls100", "prot36_pion34_tpcCls110"//, "prot36_pion34_tpcCls115", "prot36_pion34_tpcCls120"
      
      
      "prot40_pion36_spdAny", "prot40_pion36_spdFirst", "prot40_pion36_itsCls3", "prot40_pion36_itsCls4",
      "prot40_pion36_itsChi2_19",  "prot40_pion36_itsChi2_17", "prot40_pion36_itsChi2_13", 
      //"prot40_pion36_itsChi2_11", "prot40_pion36_itsChi2_9", "prot40_pion36_itsChi2_7",
      "prot40_pion36_tpcChi2_60", "prot40_pion36_tpcChi2_50", "prot40_pion36_tpcChi2_35", 
      //"prot40_pion36_tpcChi2_30", "prot40_pion36_tpcChi2_25",
      "prot40_pion36_tpcSegm4", "prot40_pion36_tpcSegm5", "prot40_pion36_tpcSegm7",
      "prot40_pion36_tpcCls90", "prot40_pion36_tpcCls100", "prot40_pion36_tpcCls110"//, "prot40_pion36_tpcCls115", "prot40_pion36_tpcCls120"
   };
   const Char_t* settingNamesShort[kNsettings] = {
      "SPD any", "SPD first", "ITS ncls #geq 3", "ITS ncls #geq 4",
      "ITS #chi^2 #GT 19",  "ITS #chi^2 #GT 17", "ITS #chi^2 #GT 13", 
      //"ITS #chi^2 #GT 11", "ITS #chi^2 #GT 9", "ITS #chi^2 #GT 7",
      "TPC #chi^2 #GT 6", "TPC #chi^2 #GT 5", "TPC #chi^2 #GT 3.5", 
      //"TPC #chi^2 #GT 3", "TPC #chi^2 #GT 2.5",
      "TPC nsegm #geq 4", "TPC nsegm #geq 5", "TPC nsegm #geq 7",
      "TPC ncls #geq 90", "TPC ncls #geq 100", "TPC ncls #geq 110"//, "TPC ncls #geq 115", "TPC ncls #geq 120"
   };*/
   
   
   
   const Int_t kNsettings = 48;
   //const Int_t kNsettings = 42;
   const Char_t* settingNames[kNsettings] = {
      "electron30_prot36_pion34", "electron30_prot36_pion36", "electron30_prot36_pion38", "electron30_prot36_pion40",
      "electron30_prot38_pion34", "electron30_prot38_pion36", "electron30_prot38_pion38", "electron30_prot38_pion40",
      "electron30_prot40_pion34", "electron30_prot40_pion36", "electron30_prot40_pion38", "electron30_prot40_pion40",
      "electron30_prot42_pion34", "electron30_prot42_pion36", "electron30_prot42_pion38", "electron30_prot42_pion40",
      
      "electron28_prot36_pion34", "electron28_prot36_pion36", "electron28_prot36_pion38", "electron28_prot36_pion40",
      "electron28_prot38_pion34", "electron28_prot38_pion36", "electron28_prot38_pion38", "electron28_prot38_pion40",
      "electron28_prot40_pion34", "electron28_prot40_pion36", "electron28_prot40_pion38", "electron28_prot40_pion40",
      "electron28_prot42_pion34", "electron28_prot42_pion36", "electron28_prot42_pion38", "electron28_prot42_pion40",
      
      "electron26_prot36_pion34", "electron26_prot36_pion36", "electron26_prot36_pion38", "electron26_prot36_pion40",
      "electron26_prot38_pion34", "electron26_prot38_pion36", "electron26_prot38_pion38", "electron26_prot38_pion40",
      "electron26_prot40_pion34", "electron26_prot40_pion36", "electron26_prot40_pion38", "electron26_prot40_pion40",
      "electron26_prot42_pion34", "electron26_prot42_pion36", "electron26_prot42_pion38", "electron26_prot42_pion40",
   };
   const Char_t* settingNamesShort[kNsettings] = {
      "|n#sigma(e)|<3, n#sigma(p)>3.6, n#sigma(#pi)>3.4", "|n#sigma(e)|<3, n#sigma(p)>3.6, n#sigma(#pi)>3.6",
      "|n#sigma(e)|<3, n#sigma(p)>3.6, n#sigma(#pi)>3.8", "|n#sigma(e)|<3, n#sigma(p)>3.6, n#sigma(#pi)>4.0",
      "|n#sigma(e)|<3, n#sigma(p)>3.8, n#sigma(#pi)>3.4", "|n#sigma(e)|<3, n#sigma(p)>3.8, n#sigma(#pi)>3.6",
      "|n#sigma(e)|<3, n#sigma(p)>3.8, n#sigma(#pi)>3.8", "|n#sigma(e)|<3, n#sigma(p)>3.8, n#sigma(#pi)>4.0",
      "|n#sigma(e)|<3, n#sigma(p)>4.0, n#sigma(#pi)>3.4", "|n#sigma(e)|<3, n#sigma(p)>4.0, n#sigma(#pi)>3.6",
      "|n#sigma(e)|<3, n#sigma(p)>4.0, n#sigma(#pi)>3.8", "|n#sigma(e)|<3, n#sigma(p)>4.0, n#sigma(#pi)>4.0",
      "|n#sigma(e)|<3, n#sigma(p)>4.2, n#sigma(#pi)>3.4", "|n#sigma(e)|<3, n#sigma(p)>4.2, n#sigma(#pi)>3.6",
      "|n#sigma(e)|<3, n#sigma(p)>4.2, n#sigma(#pi)>3.8", "|n#sigma(e)|<3, n#sigma(p)>4.2, n#sigma(#pi)>4.0",
      
      "|n#sigma(e)|<2.8, n#sigma(p)>3.6, n#sigma(#pi)>3.4", "|n#sigma(e)|<2.8, n#sigma(p)>3.6, n#sigma(#pi)>3.6",
      "|n#sigma(e)|<2.8, n#sigma(p)>3.6, n#sigma(#pi)>3.8", "|n#sigma(e)|<2.8, n#sigma(p)>3.6, n#sigma(#pi)>4.0",
      "|n#sigma(e)|<2.8, n#sigma(p)>3.8, n#sigma(#pi)>3.4", "|n#sigma(e)|<2.8, n#sigma(p)>3.8, n#sigma(#pi)>3.6",
      "|n#sigma(e)|<2.8, n#sigma(p)>3.8, n#sigma(#pi)>3.8", "|n#sigma(e)|<2.8, n#sigma(p)>3.8, n#sigma(#pi)>4.0",
      "|n#sigma(e)|<2.8, n#sigma(p)>4.0, n#sigma(#pi)>3.4", "|n#sigma(e)|<2.8, n#sigma(p)>4.0, n#sigma(#pi)>3.6",
      "|n#sigma(e)|<2.8, n#sigma(p)>4.0, n#sigma(#pi)>3.8", "|n#sigma(e)|<2.8, n#sigma(p)>4.0, n#sigma(#pi)>4.0",
      "|n#sigma(e)|<2.8, n#sigma(p)>4.2, n#sigma(#pi)>3.4", "|n#sigma(e)|<2.8, n#sigma(p)>4.2, n#sigma(#pi)>3.6",
      "|n#sigma(e)|<2.8, n#sigma(p)>4.2, n#sigma(#pi)>3.8", "|n#sigma(e)|<2.8, n#sigma(p)>4.2, n#sigma(#pi)>4.0",
      
      "|n#sigma(e)|<2.6, n#sigma(p)>3.6, n#sigma(#pi)>3.4", "|n#sigma(e)|<2.6, n#sigma(p)>3.6, n#sigma(#pi)>3.6",
      "|n#sigma(e)|<2.6, n#sigma(p)>3.6, n#sigma(#pi)>3.8", "|n#sigma(e)|<2.6, n#sigma(p)>3.6, n#sigma(#pi)>4.0",
      "|n#sigma(e)|<2.6, n#sigma(p)>3.8, n#sigma(#pi)>3.4", "|n#sigma(e)|<2.6, n#sigma(p)>3.8, n#sigma(#pi)>3.6",
      "|n#sigma(e)|<2.6, n#sigma(p)>3.8, n#sigma(#pi)>3.8", "|n#sigma(e)|<2.6, n#sigma(p)>3.8, n#sigma(#pi)>4.0",
      "|n#sigma(e)|<2.6, n#sigma(p)>4.0, n#sigma(#pi)>3.4", "|n#sigma(e)|<2.6, n#sigma(p)>4.0, n#sigma(#pi)>3.6",
      "|n#sigma(e)|<2.6, n#sigma(p)>4.0, n#sigma(#pi)>3.8", "|n#sigma(e)|<2.6, n#sigma(p)>4.0, n#sigma(#pi)>4.0",
      "|n#sigma(e)|<2.6, n#sigma(p)>4.2, n#sigma(#pi)>3.4", "|n#sigma(e)|<2.6, n#sigma(p)>4.2, n#sigma(#pi)>3.6",
      "|n#sigma(e)|<2.6, n#sigma(p)>4.2, n#sigma(#pi)>3.8", "|n#sigma(e)|<2.6, n#sigma(p)>4.2, n#sigma(#pi)>4.0"
   };
   
   const Int_t kNSigWindowSettings = 15;
   Double_t sigWindows[kNSigWindowSettings][2] = {
      {2.92, 3.16}, {2.88, 3.16}, {2.84, 3.16}, {2.80, 3.16}, {2.76, 3.16},
      {2.92, 3.20}, {2.88, 3.20}, {2.84, 3.20}, {2.80, 3.20}, {2.76, 3.20},
      {2.92, 3.12}, {2.88, 3.12}, {2.84, 3.12}, {2.80, 3.12}, {2.76, 3.12}
      //{2.92, 3.08}, {2.88, 3.08}, {2.84, 3.08}, {2.80, 3.08}, {2.76, 3.08} 
   };
   const Int_t kNFitRangeSettings = 7;
   Double_t fitRanges[kNFitRangeSettings][2] = {
      {1.5, 4.2}, {1.3, 4.4}, {1.1, 4.6}, {1.7, 4.0}, {1.9, 3.8},
      {1.1, 3.8}, {1.9, 4.6}
     // {1.0, 4.5}, {1.0, 4.3}, {1.0, 4.1}, {1.0, 3.9},
     // {1.2, 4.5}, {1.2, 4.3}, {1.2, 4.1}, {1.2, 3.9}
   };
   const Int_t kNExclusionRangeSettings = 5;
   Double_t exclRanges[kNExclusionRangeSettings][2] = {
      {2.74, 3.2}, {2.62, 3.2}, {2.50, 3.2}, {2.38, 3.2}, {2.26, 3.2}  
      //{2.5, 3.2}
   };
   const Int_t kNTotalSigFitSettings = kNSigWindowSettings+kNFitRangeSettings*kNExclusionRangeSettings-1;  
   // the "-1" is the double counting of the "standard" case
   
   const Char_t* centStrings[3] = {"0_20", "20_90", "0_90"};
   const Double_t centLims[3][2] = {
      {0.,20.}, {20.,90.}, {0.,90.}
   };
   
   TH1D* signalHist = new TH1D(Form("signalHist_%s", centStrings[centBin]), Form("Signal counts, centrality %.0f - %.0f #%", centLims[centBin][0], centLims[centBin][1]), kNsettings, -0.5, -0.5+kNsettings);
   signalHist->GetYaxis()->SetTitle("counts");
   signalHist->GetXaxis()->LabelsOption("v");
   TH2D* signalHist2D = new TH2D(Form("signalHist2D_%s", centStrings[centBin]), Form("Signal counts, centrality %.0f - %.0f #%", centLims[centBin][0], centLims[centBin][1]), kNsettings, -0.5, -0.5+kNsettings, kNTotalSigFitSettings, -0.5, -0.5+kNTotalSigFitSettings);
   signalHist2D->GetZaxis()->SetTitle("counts");
   signalHist2D->GetXaxis()->LabelsOption("v");
   
   TH1D* effHist = new TH1D(Form("effHist_%s", centStrings[centBin]), Form("Efficiency, centrality %.0f - %.0f #%", centLims[centBin][0], centLims[centBin][1]), kNsettings, -0.5, -0.5+kNsettings);
   effHist->GetYaxis()->SetTitle("acc x eff");
   effHist->GetXaxis()->LabelsOption("v");
   TH2D* effHist2D = new TH2D(Form("effHist2D_%s", centStrings[centBin]), Form("Efficiency, centrality %.0f - %.0f #%", centLims[centBin][0], centLims[centBin][1]), kNsettings, -0.5, -0.5+kNsettings, kNTotalSigFitSettings, -0.5, -0.5+kNTotalSigFitSettings);
   effHist2D->GetZaxis()->SetTitle("acc x eff");
   effHist2D->GetXaxis()->LabelsOption("v");
   
   TH1D* sigFracHist = new TH1D(Form("signalFracHist_%s", centStrings[centBin]), Form("Signal fraction in counting window, centrality %.0f - %.0f #%", centLims[centBin][0], centLims[centBin][1]), kNsettings, -0.5, -0.5+kNsettings);
   sigFracHist->GetXaxis()->LabelsOption("v");
   TH2D* sigFracHist2D = new TH2D(Form("signalFracHist2D_%s", centStrings[centBin]), Form("Signal fraction in counting window, centrality %.0f - %.0f #%", centLims[centBin][0], centLims[centBin][1]), kNsettings, -0.5, -0.5+kNsettings, kNTotalSigFitSettings, -0.5, -0.5+kNTotalSigFitSettings);
   sigFracHist2D->GetXaxis()->LabelsOption("v");
   
   TH1D* taaHist = new TH1D(Form("Taa_%s", centStrings[centBin]), "TAA values", kNsettings, -0.5, -0.5+kNsettings);
   taaHist->GetXaxis()->LabelsOption("v");
   TH1D* signifHist = new TH1D(Form("signifHist_%s", centStrings[centBin]), Form("Significance, centrality %.0f - %.0f #%", centLims[centBin][0], centLims[centBin][1]), kNsettings, -0.5, -0.5+kNsettings);
   signifHist->GetYaxis()->SetTitle("Significance");
   signifHist->GetXaxis()->LabelsOption("v");
   
   TH2D* signifHist2D = new TH2D(Form("signifHist2D_%s", centStrings[centBin]), Form("Significance, centrality %.0f - %.0f #%", centLims[centBin][0], centLims[centBin][1]), kNsettings, -0.5, -0.5+kNsettings, kNTotalSigFitSettings, -0.5, -0.5+kNTotalSigFitSettings);
   signifHist2D->GetZaxis()->SetTitle("Significance");
   signifHist2D->GetXaxis()->LabelsOption("v");
   
   TH1D* sOverBHist = new TH1D(Form("sOverBHist_%s", centStrings[centBin]), Form("S/B, centrality %.0f - %.0f #%", centLims[centBin][0], centLims[centBin][1]), kNsettings, -0.5, -0.5+kNsettings);
   sOverBHist->GetYaxis()->SetTitle("S/B");
   sOverBHist->GetXaxis()->LabelsOption("v");
   TH2D* sOverBHist2D = new TH2D(Form("sOverBHist2D_%s", centStrings[centBin]), Form("S/B, centrality %.0f - %.0f #%", centLims[centBin][0], centLims[centBin][1]), kNsettings, -0.5, -0.5+kNsettings, kNTotalSigFitSettings, -0.5, -0.5+kNTotalSigFitSettings);
   sOverBHist2D->GetZaxis()->SetTitle("S/B");
   sOverBHist2D->GetXaxis()->LabelsOption("v");
   
   TH1D* chi2Hist = new TH1D(Form("chi2Hist_%s", centStrings[centBin]), Form("#chi^{2}, centrality %.0f - %.0f #%", centLims[centBin][0], centLims[centBin][1]), kNsettings, -0.5, -0.5+kNsettings);
   chi2Hist->GetYaxis()->SetTitle("#chi^{2}");
   chi2Hist->GetXaxis()->LabelsOption("v");
   TH2D* chi2Hist2D = new TH2D(Form("chi2Hist2D_%s", centStrings[centBin]), Form("#chi^{2}, centrality %.0f - %.0f #%", centLims[centBin][0], centLims[centBin][1]), kNsettings, -0.5, -0.5+kNsettings, kNTotalSigFitSettings, -0.5, -0.5+kNTotalSigFitSettings);
   chi2Hist2D->GetZaxis()->SetTitle("#chi^{2}");
   chi2Hist2D->GetXaxis()->LabelsOption("v");
   
   TH1D* chi2MCHist = new TH1D(Form("chi2MCHist_%s", centStrings[centBin]), Form("#chi^{2} MC, centrality %.0f - %.0f #%", centLims[centBin][0], centLims[centBin][1]), kNsettings, -0.5, -0.5+kNsettings);
   chi2MCHist->GetYaxis()->SetTitle("#chi^{2} MC");
   chi2MCHist->GetXaxis()->LabelsOption("v");
   TH2D* chi2MCHist2D = new TH2D(Form("chi2MCHist2D_%s", centStrings[centBin]), Form("#chi^{2} MC, centrality %.0f - %.0f #%", centLims[centBin][0], centLims[centBin][1]), kNsettings, -0.5, -0.5+kNsettings, kNTotalSigFitSettings, -0.5, -0.5+kNTotalSigFitSettings);
   chi2MCHist2D->GetZaxis()->SetTitle("#chi^{2} MC");
   chi2MCHist2D->GetXaxis()->LabelsOption("v");
   
   for(Int_t i=1; i<=kNsettings; ++i) {
      signalHist->GetXaxis()->SetBinLabel(i, settingNamesShort[i-1]);      
      effHist->GetXaxis()->SetBinLabel(i, settingNamesShort[i-1]);
      sigFracHist->GetXaxis()->SetBinLabel(i, settingNamesShort[i-1]);
      signifHist->GetXaxis()->SetBinLabel(i, settingNamesShort[i-1]);      
      sOverBHist->GetXaxis()->SetBinLabel(i, settingNamesShort[i-1]);      
      chi2Hist->GetXaxis()->SetBinLabel(i, settingNamesShort[i-1]);      
      chi2MCHist->GetXaxis()->SetBinLabel(i, settingNamesShort[i-1]);      
      
      signalHist2D->GetXaxis()->SetBinLabel(i, settingNamesShort[i-1]); 
      effHist2D->GetXaxis()->SetBinLabel(i, settingNamesShort[i-1]); 
      sigFracHist2D->GetXaxis()->SetBinLabel(i, settingNamesShort[i-1]);
      signifHist2D->GetXaxis()->SetBinLabel(i, settingNamesShort[i-1]);      
      sOverBHist2D->GetXaxis()->SetBinLabel(i, settingNamesShort[i-1]);      
      chi2Hist2D->GetXaxis()->SetBinLabel(i, settingNamesShort[i-1]);      
      chi2MCHist2D->GetXaxis()->SetBinLabel(i, settingNamesShort[i-1]);      
   }
   
   for(Int_t i=0; i<kNSigWindowSettings; ++i) {
      signalHist2D->GetYaxis()->SetBinLabel(i+1, Form("msig [%.2f - %.2f]", sigWindows[i][0], sigWindows[i][1])); 
      effHist2D->GetYaxis()->SetBinLabel(i+1, Form("msig [%.2f - %.2f]", sigWindows[i][0], sigWindows[i][1])); 
      sigFracHist2D->GetYaxis()->SetBinLabel(i+1, Form("msig [%.2f - %.2f]", sigWindows[i][0], sigWindows[i][1])); 
      signifHist2D->GetYaxis()->SetBinLabel(i+1, Form("msig [%.2f - %.2f]", sigWindows[i][0], sigWindows[i][1])); 
      sOverBHist2D->GetYaxis()->SetBinLabel(i+1, Form("msig [%.2f - %.2f]", sigWindows[i][0], sigWindows[i][1])); 
      chi2Hist2D->GetYaxis()->SetBinLabel(i+1, Form("msig [%.2f - %.2f]", sigWindows[i][0], sigWindows[i][1])); 
      chi2MCHist2D->GetYaxis()->SetBinLabel(i+1, Form("msig [%.2f - %.2f]", sigWindows[i][0], sigWindows[i][1])); 
   }
   Int_t idxSigExtr=1;
   for(Int_t i=0; i<kNFitRangeSettings; ++i) {
      for(Int_t j=0; j<kNExclusionRangeSettings; ++j) {
         if(i==0 && j==0) continue;     // skip the standard case which is already included
         signalHist2D->GetYaxis()->SetBinLabel(kNSigWindowSettings+idxSigExtr, 
                                               Form("mfit [%.2f - %.2f], mexcl [%.2f - %.2f]", fitRanges[i][0], fitRanges[i][1], exclRanges[j][0], exclRanges[j][1])); 
         effHist2D->GetYaxis()->SetBinLabel(kNSigWindowSettings+idxSigExtr, 
                                               Form("mfit [%.2f - %.2f], mexcl [%.2f - %.2f]", fitRanges[i][0], fitRanges[i][1], exclRanges[j][0], exclRanges[j][1])); 
         sigFracHist2D->GetYaxis()->SetBinLabel(kNSigWindowSettings+idxSigExtr, 
                                            Form("mfit [%.2f - %.2f], mexcl [%.2f - %.2f]", fitRanges[i][0], fitRanges[i][1], exclRanges[j][0], exclRanges[j][1])); 
         signifHist2D->GetYaxis()->SetBinLabel(kNSigWindowSettings+idxSigExtr, 
                                               Form("mfit [%.2f - %.2f], mexcl [%.2f - %.2f]", fitRanges[i][0], fitRanges[i][1], exclRanges[j][0], exclRanges[j][1])); 
         sOverBHist2D->GetYaxis()->SetBinLabel(kNSigWindowSettings+idxSigExtr, 
                                               Form("mfit [%.2f - %.2f], mexcl [%.2f - %.2f]", fitRanges[i][0], fitRanges[i][1], exclRanges[j][0], exclRanges[j][1])); 
         chi2Hist2D->GetYaxis()->SetBinLabel(kNSigWindowSettings+idxSigExtr, 
                                               Form("mfit [%.2f - %.2f], mexcl [%.2f - %.2f]", fitRanges[i][0], fitRanges[i][1], exclRanges[j][0], exclRanges[j][1])); 
         chi2MCHist2D->GetYaxis()->SetBinLabel(kNSigWindowSettings+idxSigExtr, 
                                               Form("mfit [%.2f - %.2f], mexcl [%.2f - %.2f]", fitRanges[i][0], fitRanges[i][1], exclRanges[j][0], exclRanges[j][1])); 
         idxSigExtr++;
      }
   }
   
   for(Int_t iSetting=0; iSetting<kNsettings; ++iSetting) {
   //for(Int_t iSetting=0; iSetting<1; ++iSetting) {
      cout << "Setting " << settingNames[iSetting] << endl;
      
      // standard setting
      //gFits.SetBkgFitFunction(fitFunc);
      ExtractSignal(settingNames[iSetting], 2.92, 3.16, centLims[centBin][0], centLims[centBin][1], kTRUE);
      //gFits.PrintFitValues();
      Double_t* fitValues = gFits.GetFitValues();
      signalHist->SetBinContent(iSetting+1, fitValues[AliResonanceFits::kSig]);
      signalHist->SetBinError(iSetting+1, fitValues[AliResonanceFits::kSigErr]);
      signalHist2D->SetBinContent(iSetting+1, 1, fitValues[AliResonanceFits::kSig]);
      signalHist2D->SetBinError(iSetting+1, 1, fitValues[AliResonanceFits::kSigErr]);
      signifHist->SetBinContent(iSetting+1, fitValues[AliResonanceFits::kSignif]);
      signifHist2D->SetBinContent(iSetting+1,1, fitValues[AliResonanceFits::kSignif]);
      sOverBHist->SetBinContent(iSetting+1, fitValues[AliResonanceFits::kSoverB]);
      sOverBHist2D->SetBinContent(iSetting+1,1, fitValues[AliResonanceFits::kSoverB]);
      chi2Hist->SetBinContent(iSetting+1, fitValues[AliResonanceFits::kChisqSideBands]);
      chi2Hist2D->SetBinContent(iSetting+1,1, fitValues[AliResonanceFits::kChisqSideBands]);
      chi2MCHist->SetBinContent(iSetting+1, fitValues[AliResonanceFits::kChisqMCTotal]);
      chi2MCHist2D->SetBinContent(iSetting+1,1, fitValues[AliResonanceFits::kChisqMCTotal]);
      
      Double_t eff=0.; Double_t effErr=0.; Double_t sigFrac=0.;
      ComputeEfficiency(settingNames[iSetting], eff, effErr, sigFrac, 2.92, 3.16, centLims[centBin][0], centLims[centBin][1]);
      effHist->SetBinContent(iSetting+1, eff);
      effHist->SetBinError(iSetting+1, effErr);
      effHist2D->SetBinContent(iSetting+1,1, eff);
      effHist2D->SetBinError(iSetting+1,1, effErr);
      
      sigFracHist->SetBinContent(iSetting+1, sigFrac);
      sigFracHist2D->SetBinContent(iSetting+1,1, sigFrac);
      
      taaHist->SetBinContent(iSetting+1, 1000.*gTAA_XeXe[centBin][0]);         // factor 1000 to express this in 1/b
      taaHist->SetBinError(iSetting+1, 1000.0*gTAA_XeXe[centBin][2]);
      
      if(savePictures) {
         TString pictureBaseName = settingNames[iSetting];
         pictureBaseName += Form("_%s", centStrings[centBin]);
         pictureBaseName += "_Sig2.92_3.16";
         pictureBaseName += Form("_MFit%.2f_%.2f", gFits.GetMassFitRange()[0], gFits.GetMassFitRange()[1]);
         pictureBaseName += Form("_MExcl%.2f_%.2f", gFits.GetMassExclusionRange()[0], gFits.GetMassExclusionRange()[1]);
         Draw(&gFits, 2.92, 3.16, kTRUE, outputDir, pictureBaseName.Data());
      }
      
      continue;
      // mass counting window variations
      for(Int_t imass=1; imass<kNSigWindowSettings; ++imass) {
         gFits.ComputeOutputValues(sigWindows[imass][0], sigWindows[imass][1]);
         fitValues = gFits.GetFitValues();
         signalHist2D->SetBinContent(iSetting+1, 1+imass, fitValues[AliResonanceFits::kSig]);
         signalHist2D->SetBinError(iSetting+1, 1+imass, fitValues[AliResonanceFits::kSigErr]);
         signifHist2D->SetBinContent(iSetting+1,1+imass, fitValues[AliResonanceFits::kSignif]);
         sOverBHist2D->SetBinContent(iSetting+1,1+imass, fitValues[AliResonanceFits::kSoverB]);
         chi2Hist2D->SetBinContent(iSetting+1,1+imass, fitValues[AliResonanceFits::kChisqSideBands]);
         chi2MCHist2D->SetBinContent(iSetting+1,1+imass, fitValues[AliResonanceFits::kChisqMCTotal]);
         
         ComputeEfficiency(settingNames[iSetting], eff, effErr, sigFrac, sigWindows[imass][0], sigWindows[imass][1], centLims[centBin][0], centLims[centBin][1]);
         effHist2D->SetBinContent(iSetting+1,1+imass, eff);
         effHist2D->SetBinError(iSetting+1,1+imass, effErr);
         sigFracHist2D->SetBinContent(iSetting+1,1+imass, sigFrac);
         
         /*if(savePictures) {
            TString pictureBaseName = settingNames[iSetting];
            pictureBaseName += Form("_%s", centStrings[centBin]);
            pictureBaseName += Form("_Sig%.2f_%.2f", sigWindows[imass][0], sigWindows[imass][1]);
            pictureBaseName += Form("_MFit%.2f_%.2f", gFits.GetMassFitRange()[0], gFits.GetMassFitRange()[1]);
            pictureBaseName += Form("_MExcl%.2f_%.2f", gFits.GetMassExclusionRange()[0], gFits.GetMassExclusionRange()[1]);
            Draw(&gFits, sigWindows[imass][0], sigWindows[imass][1], kTRUE, outputDir, pictureBaseName.Data());
         }*/          
      }  // end loop over mass windows
      
      // recompute efficiency for the standard counting window
      ComputeEfficiency(settingNames[iSetting], eff, effErr, sigFrac, 2.92, 3.16, centLims[centBin][0], centLims[centBin][1]);
      
      // systematic variations on the mass fit and exclusion windows
      idxSigExtr=1;
      for(Int_t ifit=0; ifit<kNFitRangeSettings; ++ifit) {
         gFits.SetMassFitRange(fitRanges[ifit][0], fitRanges[ifit][1]);
         gFits.SetBkgFitFunction(fitFunc);
         for(Int_t iexcl=0; iexcl<kNExclusionRangeSettings; ++iexcl) {
            if(ifit==0 && iexcl==0) continue;
            gFits.SetMassExclusionRange(exclRanges[iexcl][0], exclRanges[iexcl][1]);
            
            ExtractSignal(settingNames[iSetting], 2.92, 3.16, centLims[centBin][0], centLims[centBin][1], kTRUE);
            //gFits.PrintFitValues();
            fitValues = gFits.GetFitValues();
            signalHist2D->SetBinContent(iSetting+1, kNSigWindowSettings+idxSigExtr, fitValues[AliResonanceFits::kSig]);
            signalHist2D->SetBinError(iSetting+1, kNSigWindowSettings+idxSigExtr, fitValues[AliResonanceFits::kSigErr]);
            signifHist2D->SetBinContent(iSetting+1,kNSigWindowSettings+idxSigExtr, fitValues[AliResonanceFits::kSignif]);
            sOverBHist2D->SetBinContent(iSetting+1,kNSigWindowSettings+idxSigExtr, fitValues[AliResonanceFits::kSoverB]);
            chi2Hist2D->SetBinContent(iSetting+1,kNSigWindowSettings+idxSigExtr, fitValues[AliResonanceFits::kChisqSideBands]);
            chi2MCHist2D->SetBinContent(iSetting+1,kNSigWindowSettings+idxSigExtr, fitValues[AliResonanceFits::kChisqMCTotal]);
            
            effHist2D->SetBinContent(iSetting+1, kNSigWindowSettings+idxSigExtr, eff);
            effHist2D->SetBinError(iSetting+1, kNSigWindowSettings+idxSigExtr, effErr);
            sigFracHist2D->SetBinContent(iSetting+1, kNSigWindowSettings+idxSigExtr, sigFrac);
            
            /*if(savePictures) {
               TString pictureBaseName = settingNames[iSetting];
               pictureBaseName += Form("_%s", centStrings[centBin]);
               pictureBaseName += "_Sig2.92_3.16";
               pictureBaseName += Form("_MFit%.2f_%.2f", gFits.GetMassFitRange()[0], gFits.GetMassFitRange()[1]);
               pictureBaseName += Form("_MExcl%.2f_%.2f", gFits.GetMassExclusionRange()[0], gFits.GetMassExclusionRange()[1]);
               Draw(&gFits, 2.92, 3.16, kTRUE, outputDir, pictureBaseName.Data());
            }*/
            
            idxSigExtr++;
         }
      }
      // reset defaults
      gFits.SetMassFitRange(1.5, 4.2);
      gFits.SetMassExclusionRange(2.5, 3.2);
   }   // end loop over settings
   signalHist->Sumw2();
   signalHist2D->Sumw2();
   effHist->Sumw2();
   effHist2D->Sumw2();
   taaHist->Sumw2();
   
   // Compute corrected J/psi yields in Xe-Xe
   TH1D* corrYieldHist = (TH1D*) signalHist->Clone(Form("corrYieldHist_%s", centStrings[centBin]));
   corrYieldHist->SetTitle(Form("Corrected J/#psi yield, centrality %.0f - %.0f #%", centLims[centBin][0], centLims[centBin][1]));
   corrYieldHist->GetYaxis()->SetTitle("B.R. x dN/dy");
   corrYieldHist->GetXaxis()->LabelsOption("v");
   corrYieldHist->Divide(effHist);
   corrYieldHist->Scale(1./1.8);               // normalize to the delta y interval y<|0.9|
   
   TH2D* corrYieldHist2D = (TH2D*) signalHist2D->Clone(Form("corrYieldHist2D_%s", centStrings[centBin]));
   corrYieldHist2D->SetTitle(Form("Corrected J/#psi yield, centrality %.0f - %.0f #%", centLims[centBin][0], centLims[centBin][1]));
   corrYieldHist2D->GetYaxis()->SetTitle("B.R. x dN/dy");
   corrYieldHist2D->GetXaxis()->LabelsOption("v");
   corrYieldHist2D->Divide(effHist2D);
   corrYieldHist2D->Scale(1./1.8);               // normalize to the delta y interval y<|0.9|
   TH2D* corrYieldDistrib2D_sigShape = new TH2D(Form("corrYieldDistrib2D_sigShape_%s", centStrings[centBin]), "", kNsettings, -0.5, -0.5+kNsettings, 50, 0.0, 0.003);
   corrYieldDistrib2D_sigShape->GetYaxis()->SetTitle("B.R. x dN/dy");
   corrYieldDistrib2D_sigShape->GetXaxis()->LabelsOption("v");
   TH2D* corrYieldDistrib2D_bkgMatching = new TH2D(Form("corrYieldDistrib2D_bkgMatching_%s", centStrings[centBin]), "", kNsettings, -0.5, -0.5+kNsettings, 50, 0.0, 0.003);
   corrYieldDistrib2D_bkgMatching->GetYaxis()->SetTitle("B.R. x dN/dy");
   corrYieldDistrib2D_bkgMatching->GetXaxis()->LabelsOption("v");
   TH1D* corrYieldSigExtractionSystematics_sigShape = new TH1D(Form("corrYieldSigExtractionSystematics_sigShape_%s", centStrings[centBin]), "", kNsettings, -0.5, -0.5+kNsettings);
   corrYieldSigExtractionSystematics_sigShape->GetXaxis()->LabelsOption("v");
   TH1D* corrYieldSigExtractionSystematics_bkgMatching = new TH1D(Form("corrYieldSigExtractionSystematics_bkgMatching_%s", centStrings[centBin]), "", kNsettings, -0.5, -0.5+kNsettings);
   corrYieldSigExtractionSystematics_bkgMatching->GetXaxis()->LabelsOption("v");
   for(Int_t ib=1;ib<=kNsettings;++ib) {
      corrYieldDistrib2D_sigShape->GetXaxis()->SetBinLabel(ib, settingNames[ib-1]);
      corrYieldDistrib2D_bkgMatching->GetXaxis()->SetBinLabel(ib, settingNames[ib-1]);
      corrYieldSigExtractionSystematics_sigShape->GetXaxis()->SetBinLabel(ib, settingNames[ib-1]);
      corrYieldSigExtractionSystematics_bkgMatching->GetXaxis()->SetBinLabel(ib, settingNames[ib-1]);
   }
   
   // normalize to the number of events
   Double_t nevents = eventsVsCent->Integral(eventsVsCent->GetXaxis()->FindBin(centLims[centBin][0]+1.0e-6), eventsVsCent->GetXaxis()->FindBin(centLims[centBin][1]-1.0e-6));
   corrYieldHist->Scale(1.0/nevents);
   corrYieldHist2D->Scale(1.0/nevents);
   
   TH1D* corrYieldDistrib = new TH1D("corrYieldDistrib", "", 50, 0.0, 0.003);
   corrYieldDistrib->GetXaxis()->SetTitle("B.R. x dN/dy");
   for(Int_t i=1; i<=corrYieldHist->GetXaxis()->GetNbins(); ++i)
      corrYieldDistrib->Fill(corrYieldHist->GetBinContent(i));
   
   TH1D* absDeviations = (TH1D*) corrYieldHist->Clone(Form("absDeviations_%s", centStrings[centBin]));
   absDeviations->SetTitle(Form("Corrected yield deviations from standard setting #Delta Y= (Y-Y_{0}), centrality %.0f - %.0f #%", centLims[centBin][0], centLims[centBin][1]));
   absDeviations->GetYaxis()->SetTitle("#Delta Y = Y_{0}-Y_{i}");
   absDeviations->GetXaxis()->LabelsOption("v");
   
   TH2D* absDeviations2D = (TH2D*) corrYieldHist2D->Clone(Form("absDeviations2D_%s", centStrings[centBin]));
   absDeviations2D->SetTitle(Form("Corrected yield deviations from standard setting #Delta Y= (Y-Y_{0}), centrality %.0f - %.0f #%", centLims[centBin][0], centLims[centBin][1]));
   absDeviations2D->GetZaxis()->SetTitle("#Delta Y = Y_{0}-Y_{i}");
   absDeviations2D->GetXaxis()->LabelsOption("v");
   
   TH1D* relDeviations = (TH1D*) corrYieldHist->Clone(Form("relDeviations_%s", centStrings[centBin]));
   relDeviations->SetTitle(Form("Relative yield deviations from standard setting, centrality %.0f - %.0f #%", centLims[centBin][0], centLims[centBin][1]));
   relDeviations->GetYaxis()->SetTitle("#Delta Y / Y (percents)");
   relDeviations->GetXaxis()->LabelsOption("v");
   
   TH2D* relDeviations2D = (TH2D*) corrYieldHist2D->Clone(Form("relDeviations2D_%s", centStrings[centBin]));
   relDeviations2D->SetTitle(Form("Relative yield deviations from standard setting, centrality %.0f - %.0f #%", centLims[centBin][0], centLims[centBin][1]));
   relDeviations2D->GetZaxis()->SetTitle("percents");
   relDeviations2D->GetXaxis()->LabelsOption("v");
   
   TH1D* tempDndyHist1 = new TH1D("tempDndyHist1","", 50, 0.0, 0.003);
   TH1D* tempDndyHist2 = new TH1D("tempDndyHist2","", 50, 0.0, 0.003);
   for(Int_t iSetting=0; iSetting<kNsettings; ++iSetting) {
      cout << "iSetting " << iSetting << endl;
      
      absDeviations->SetBinContent(iSetting+1, corrYieldHist->GetBinContent(iSetting+1) - corrYieldHist->GetBinContent(1));
      absDeviations->SetBinError(iSetting+1, TMath::Sqrt(TMath::Abs(corrYieldHist->GetBinError(iSetting+1)*corrYieldHist->GetBinError(iSetting+1) - 
                                                                                                            corrYieldHist->GetBinError(1)*corrYieldHist->GetBinError(1))));
      
      if(absDeviations->GetBinError(iSetting+1)>1.0e-6) {
         relDeviations->SetBinContent(iSetting+1, 100.0*absDeviations->GetBinContent(iSetting+1) / corrYieldHist->GetBinContent(1));
         relDeviations->SetBinError(iSetting+1, 100.0*absDeviations->GetBinError(iSetting+1) / corrYieldHist->GetBinContent(1));
      }
      else
         relDeviations->SetBinContent(iSetting+1, 0.0);
      
      for(Int_t iSigSetting=0; iSigSetting<kNTotalSigFitSettings; ++iSigSetting) {
         absDeviations2D->SetBinContent(iSetting+1, iSigSetting+1, corrYieldHist2D->GetBinContent(iSetting+1,iSigSetting+1)-
                                                                                                         corrYieldHist2D->GetBinContent(iSetting+1,1));
         absDeviations2D->SetBinError(iSetting+1, iSigSetting+1, 
                                      TMath::Sqrt(TMath::Abs(corrYieldHist2D->GetBinError(iSetting+1,iSigSetting+1) * 
                                                                            corrYieldHist2D->GetBinError(iSetting+1,iSigSetting+1) - 
                                                                             corrYieldHist2D->GetBinError(iSetting+1,1)*corrYieldHist2D->GetBinError(iSetting+ 1,1))));
         
         if(absDeviations2D->GetBinError(iSetting+1, iSigSetting+1)>1.0e-6) {
            relDeviations2D->SetBinContent(iSetting+1, iSigSetting+1, 
                                           100.0*absDeviations2D->GetBinContent(iSetting+1, iSigSetting+1) / corrYieldHist->GetBinContent(iSetting+1,1));
            relDeviations2D->SetBinError(iSetting+1, iSigSetting+1, 
                                         100.0*absDeviations2D->GetBinError(iSetting+1, iSigSetting+1) / corrYieldHist->GetBinContent(iSetting+1,1));
         }
         else
            relDeviations2D->SetBinContent(iSetting+1, iSigSetting+1, 0.0);
         
         if(iSigSetting<kNSigWindowSettings) { 
            corrYieldDistrib2D_sigShape->Fill(iSetting, corrYieldHist2D->GetBinContent(iSetting+1,iSigSetting+1));
            tempDndyHist1->Fill(corrYieldHist2D->GetBinContent(iSetting+1,iSigSetting+1));
         }
         else {
            corrYieldDistrib2D_bkgMatching->Fill(iSetting, corrYieldHist2D->GetBinContent(iSetting+1,iSigSetting+1));
            tempDndyHist2->Fill(corrYieldHist2D->GetBinContent(iSetting+1,iSigSetting+1));
         }
      }      // end loop over signal extraction cases
      
      corrYieldSigExtractionSystematics_sigShape->SetBinContent(iSetting+1, tempDndyHist1->GetMean());
      corrYieldSigExtractionSystematics_sigShape->SetBinError(iSetting+1, tempDndyHist1->GetRMS());
      corrYieldSigExtractionSystematics_bkgMatching->SetBinContent(iSetting+1, tempDndyHist2->GetMean());
      corrYieldSigExtractionSystematics_bkgMatching->SetBinError(iSetting+1, tempDndyHist2->GetRMS());
      tempDndyHist1->Reset();
      tempDndyHist2->Reset();
   }
   
   cout << "aaaaaaaaa 1" << endl;
   
   // compute the Raa
   TH1D* raaHist = (TH1D*) corrYieldHist->Clone(Form("raaHist_%s", centStrings[centBin]));
   raaHist->SetTitle(Form("J/#psi R_{AA} in Xe-Xe, centrality %.0f - %.0f #%", centLims[centBin][0], centLims[centBin][1]));
   raaHist->GetYaxis()->SetTitle("J/#psi R_{AA}");
   raaHist->GetXaxis()->LabelsOption("v");
   raaHist->Scale(1.0/1000.0/gTAA_XeXe[centBin][0]);
   raaHist->Scale(1.0/gPPcrossSection5TeVinterpol);
   
   TH2D* raaHist2D = (TH2D*) corrYieldHist2D->Clone(Form("raaHist2D_%s", centStrings[centBin]));
   raaHist2D->SetTitle(Form("J/#psi R_{AA} in Xe-Xe, centrality %.0f - %.0f #%", centLims[centBin][0], centLims[centBin][1]));
   raaHist2D->GetZaxis()->SetTitle("J/#psi R_{AA}");
   raaHist2D->GetXaxis()->LabelsOption("v");
   raaHist2D->Scale(1.0/1000.0/gTAA_XeXe[centBin][0]);
   raaHist2D->Scale(1.0/gPPcrossSection5TeVinterpol);
   
   TH2D* raaDistrib2D_sigShape = new TH2D(Form("raaDistrib2D_sigShape_%s", centStrings[centBin]), "", kNsettings, -0.5, -0.5+kNsettings, 50, 0.0, 2.0);
   raaDistrib2D_sigShape->GetYaxis()->SetTitle("J/#psi R_{AA}");
   raaDistrib2D_sigShape->GetXaxis()->LabelsOption("v");
   TH2D* raaDistrib2D_bkgMatching = new TH2D(Form("raaDistrib2D_bkgMatching_%s", centStrings[centBin]), "", kNsettings, -0.5, -0.5+kNsettings, 50, 0.0, 2.0);
   raaDistrib2D_bkgMatching->GetYaxis()->SetTitle("dN/dy");
   raaDistrib2D_bkgMatching->GetXaxis()->LabelsOption("v");
   TH1D* raaSigExtractionSyst_sigShape = new TH1D(Form("raaSigExtractionSyst_sigShape_%s", centStrings[centBin]), "", kNsettings, -0.5, -0.5+kNsettings);
   raaSigExtractionSyst_sigShape->GetXaxis()->LabelsOption("v");
   TH1D* raaSigExtractionSyst_bkgMatching = new TH1D(Form("raaSigExtractionSyst_bkgMatching_%s", centStrings[centBin]), "", kNsettings, -0.5, -0.5+kNsettings);
   raaSigExtractionSyst_bkgMatching->GetXaxis()->LabelsOption("v");
   for(Int_t ib=1;ib<=kNsettings;++ib) {
      raaDistrib2D_sigShape->GetXaxis()->SetBinLabel(ib, settingNames[ib-1]);
      raaDistrib2D_bkgMatching->GetXaxis()->SetBinLabel(ib, settingNames[ib-1]);
      raaSigExtractionSyst_sigShape->GetXaxis()->SetBinLabel(ib, settingNames[ib-1]);
      raaSigExtractionSyst_bkgMatching->GetXaxis()->SetBinLabel(ib, settingNames[ib-1]);
   }
   
   cout << "aaaaaaaaa 2" << endl;
   
   TH1D* raaDistrib = new TH1D(Form("raaDistrib_%s", centStrings[centBin]), "", 50, 0.0, 2.0);
   TH1D* raaErrDistrib = new TH1D(Form("raaErrDistrib_%s", centStrings[centBin]), "", 25, 0., 0.5);
   TH1D* raaDistribFromMeans_sigShape = new TH1D(Form("raaDistribFromMeans_sigShape_%s", centStrings[centBin]), "", 50, 0.0, 2.0);
   TH1D* raaDistribFromMeans_bkgMatching = new TH1D(Form("raaDistribFromMeans_bkgMatching_%s", centStrings[centBin]), "", 50, 0.0, 2.0);
   TH1D* raaDistribSystErr_sigShape = new TH1D(Form("raaDistribSystErr_sigShape_%s", centStrings[centBin]), "", 20, 0., 20.);
   raaDistribSystErr_sigShape->GetXaxis()->SetTitle("Relative syst uncertainty on J/#psi MC signal shape (percents)");
   TH1D* raaDistribSystErr_bkgMatching = new TH1D(Form("raaDistribSystErr_bkgMatching_%s", centStrings[centBin]), "", 20, 0., 20.);
   raaDistribSystErr_bkgMatching->GetXaxis()->SetTitle("Relative syst uncertainty on bkg matching (percents)");
   
   TH1D* tempHist1 = new TH1D("tempHist1", "", 50., 0.0, 2.0);
   TH1D* tempHist2 = new TH1D("tempHist2", "", 50., 0.0, 2.0);
   for(Int_t iSetting=0; iSetting<kNsettings; ++iSetting) {
      raaDistrib->Fill(raaHist->GetBinContent(iSetting+1));
      raaErrDistrib->Fill(raaHist->GetBinError(iSetting+1));
      
      continue;
      for(Int_t iSigSetting=0;iSigSetting<kNTotalSigFitSettings; ++iSigSetting) {
         if(iSigSetting<kNSigWindowSettings) {
            raaDistrib2D_sigShape->Fill(iSetting, raaHist2D->GetBinContent(iSetting+1, iSigSetting+1));
            tempHist1->Fill(raaHist2D->GetBinContent(iSetting+1, iSigSetting+1));
         }
         else {
            raaDistrib2D_bkgMatching->Fill(iSetting, raaHist2D->GetBinContent(iSetting+1, iSigSetting+1));
            tempHist2->Fill(raaHist2D->GetBinContent(iSetting+1, iSigSetting+1));
         }
      }
      raaSigExtractionSyst_sigShape->SetBinContent(iSetting+1, tempHist1->GetMean());
      raaSigExtractionSyst_sigShape->SetBinError(iSetting+1, tempHist1->GetRMS());
      raaSigExtractionSyst_bkgMatching->SetBinContent(iSetting+1, tempHist2->GetMean());
      raaSigExtractionSyst_bkgMatching->SetBinError(iSetting+1, tempHist2->GetRMS());
      raaDistribFromMeans_sigShape->Fill(tempHist1->GetMean());
      raaDistribFromMeans_bkgMatching->Fill(tempHist2->GetMean());
      raaDistribSystErr_sigShape->Fill(100.0*tempHist1->GetRMS()/tempHist1->GetMean());
      raaDistribSystErr_bkgMatching->Fill(100.0*tempHist2->GetRMS()/tempHist2->GetMean());
      tempHist1->Reset();
      tempHist2->Reset();
   }
   
   cout << "aaaaaaaaa 3" << endl;
   
   TFile* save=new TFile(Form("%s/jpsiRAASyst_%s_selectedFinal.root", outputDir, centStrings[centBin]), "RECREATE");
   signalHist->Write();
   signalHist2D->Write();
   signifHist->Write();
   signifHist2D->Write();
   sOverBHist->Write();
   sOverBHist2D->Write();
   effHist->Write();
   effHist2D->Write();
   sigFracHist->Write();
   sigFracHist2D->Write();
   taaHist->Write();
   corrYieldHist->Write();
   corrYieldDistrib->Write();
   corrYieldHist2D->Write();
   corrYieldDistrib2D_sigShape->Write();
   corrYieldDistrib2D_bkgMatching->Write();
   corrYieldSigExtractionSystematics_sigShape->Write();
   corrYieldSigExtractionSystematics_bkgMatching->Write();
   absDeviations->Write();
   absDeviations2D->Write();
   relDeviations->Write();
   relDeviations2D->Write();
   raaHist->Write();
   raaHist2D->Write();
   raaDistrib->Write();
   raaErrDistrib->Write();
   raaDistribFromMeans_sigShape->Write();
   raaDistribFromMeans_bkgMatching->Write();
   raaDistrib2D_sigShape->Write();
   raaDistrib2D_bkgMatching->Write();
   raaSigExtractionSyst_sigShape->Write();
   raaSigExtractionSyst_bkgMatching->Write();
   raaDistribSystErr_sigShape->Write();
   raaDistribSystErr_bkgMatching->Write();
   chi2Hist->Write();
   chi2Hist2D->Write();
   chi2MCHist->Write();
   chi2MCHist2D->Write();
   save->Close();
}


//________________________________________________________________________________________
void ExtractSignal(const Char_t* setting, Double_t minMass, Double_t maxMass, Double_t minCent, Double_t maxCent, Bool_t savePicture=kFALSE) {
   //
   // test signal extraction
   //
   THnF* seos = (THnF*)gHistManData->GetHistogram(Form("PairSEPM_%s",setting),"PairInvMass");
   THnF* meos = (THnF*)gHistManData->GetHistogram(Form("PairMEPM_%s",setting),"PairInvMass");
   THnF* sels1 = (THnF*)gHistManData->GetHistogram(Form("PairSEPP_%s",setting),"PairInvMass");
   THnF* sels2 = (THnF*)gHistManData->GetHistogram(Form("PairSEMM_%s",setting),"PairInvMass");
   THnF* mels1 = (THnF*)gHistManData->GetHistogram(Form("PairMEPP_%s",setting),"PairInvMass");
   THnF* mels2 = (THnF*)gHistManData->GetHistogram(Form("PairMEMM_%s",setting),"PairInvMass");
   
   THnF* seosMC = 0x0;
   if(savePicture)
      seosMC = (THnF*)gHistManMC->GetHistogram(Form("PairSEPM_%s_TrueElectron",setting),"PairInvMass");
   
   gFits.SetSEOSHistogram(seos);
   gFits.SetMEOSHistogram(meos);
   gFits.SetSELSHistograms(sels1, sels2);
   gFits.SetMELSHistograms(mels1, mels2);
   if(seosMC) gFits.SetSEOSMCHistogram(seosMC);
   
   cout << "minCent / maxCent " << minCent << "   " << maxCent << endl; 
   gFits.SetVarRange(AliReducedVarManager::kCentVZERO, minCent, maxCent);
   
   gFits.Process();
   gFits.ComputeOutputValues(minMass, maxMass);
}

//________________________________________________________________________________________
void ComputeEfficiency(const Char_t* setting, Double_t& eff, Double_t& effErr, Double_t& sigFrac, 
                                      Double_t minMass, Double_t maxMass, 
                                      Double_t centMin=0., Double_t centMax=100.,
                                      Double_t ptMin=0., Double_t ptMax=100.) {
   //
   // compute efficiency
   //
   cout << "ComputeEff 0" << endl;
   THnF* nominatorMC = (THnF*)gHistManMC->GetHistogram(Form("PairSEPM_%s_TrueElectron",setting),"PairInvMass");
   cout << "ComputeEff 0.1" << endl;
   TH3F* denominatorMC = (TH3F*)gHistManMC->GetHistogram("mcTruthJpsi_PureMCTruth_BeforeSelection","MassMC_Pt_CentVZERO");
   cout << "ComputeEff 0.2" << endl;
   nominatorMC->Sumw2();
   cout << "ComputeEff 0.3" << endl;
   //denominatorMC->Sumw2();
   
   cout << "ComputeEff 1" << endl;
   
   nominatorMC->GetAxis(1)->SetRangeUser(ptMin+1.0e-6, ptMax-1.0e-6);
   denominatorMC->GetYaxis()->SetRangeUser(ptMin+1.0e-6, ptMax-1.0e-6);
   nominatorMC->GetAxis(2)->SetRangeUser(centMin+1.0e-6, centMax-1.0e-6);
   denominatorMC->GetZaxis()->SetRangeUser(centMin+1.0e-6, centMax-1.0e-6);
   
   Double_t err=0.0;
   sigFrac = 0.0;
   /*TH3D* nominatorMCProj3D = nominatorMC->Projection(0,1,2);
   sigFrac = nominatorMCProj3D->IntegralAndError(nominatorMCProj3D->GetXaxis()->FindBin(minMass+1.0e-6), 
                                                 nominatorMCProj3D->GetXaxis()->FindBin(maxMass-1.0e-6), 
                                                 nominatorMCProj3D->GetYaxis()->FindBin(ptMin+1.0e-6), 
                                                 nominatorMCProj3D->GetYaxis()->FindBin(ptMax-1.0e-6), 
                                                 nominatorMCProj3D->GetZaxis()->FindBin(centMin+1.0e-6), 
                                                 nominatorMCProj3D->GetZaxis()->FindBin(centMax-1.0e-6),err);
   
   sigFrac /= nominatorMCProj3D->IntegralAndError(1, nominatorMCProj3D->GetXaxis()->GetNbins(), 
                                            nominatorMCProj3D->GetYaxis()->FindBin(ptMin+1.0e-6), nominatorMCProj3D->GetYaxis()->FindBin(ptMax-1.0e-6), 
                                            nominatorMCProj3D->GetZaxis()->FindBin(centMin+1.0e-6), nominatorMCProj3D->GetZaxis()->FindBin(centMax-1.0e-6),err);
   */
   //nominatorMCProj3D->GetXaxis()->SetRangeUser(minMass+1.0e-6, maxMass-1.0e-6);
   
   cout << "ComputeEff 2" << endl;
   
   TH1D* nomProjMass = nominatorMC->Projection(0); nomProjMass->SetName("nomProjMass");
   TH1D* denomProjMass = denominatorMC->ProjectionX("denomProjMass");
   TH1D* nomProjPt = nominatorMC->Projection(1); nomProjPt->SetName("nomProjPt");
   TH1D* denomProjPt = denominatorMC->ProjectionY("denomProjPt");
   TH1D* nomProjCent = nominatorMC->Projection(2); nomProjCent->SetName("nomProjCent");
   TH1D* denomProjCent = denominatorMC->ProjectionZ("denomProjCent");
   
   Double_t nomYieldErr = 0.0;
   //Double_t nomYield = nomProjPt->IntegralAndError(1, nomProjPt->GetXaxis()->GetNbins(), nomYieldErr);
   Double_t nomYield = nomProjMass->IntegralAndError(nomProjMass->GetXaxis()->FindBin(minMass+0.001), 
                                                     nomProjMass->GetXaxis()->FindBin(maxMass-0.001), nomYieldErr);
   Double_t totalRecYieldErr = 0.0;
   Double_t totalRecYield = nomProjMass->IntegralAndError(1, nomProjMass->GetXaxis()->GetNbins(), totalRecYieldErr);
   sigFrac = nomYield / totalRecYield;
   
   Double_t denomYieldErr = 0.0;
   //Double_t denomYield = denomProjPt->IntegralAndError(1, nomProjPt->GetXaxis()->GetNbins(), denomYieldErr);
   Double_t denomYield = denomProjMass->IntegralAndError(1, denomProjMass->GetXaxis()->GetNbins(), denomYieldErr);
   
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
   
   cout << "ComputeEff 3" << endl;
   
   delete nomProjMass; delete denomProjMass;
   delete nomProjPt; delete denomProjPt;
   delete nomProjCent; delete denomProjCent;
}

//________________________________________________________________________________________
void TestResonanceFits(const Char_t* filename, const Char_t* setting, const Char_t* filenameMC="", Bool_t isPreliminaryPlot=kFALSE) {
   //
   // test signal extraction
   //
   gHistManData=new AliHistogramManager();
   gHistManData->InitFile(filename,"jpsi2eeHistos");
   THnF* seos = (THnF*)gHistManData->GetHistogram(Form("PairSEPM_%s",setting),"PairInvMass");
   THnF* meos = (THnF*)gHistManData->GetHistogram(Form("PairMEPM_%s",setting),"PairInvMass");
   THnF* sels1 = (THnF*)gHistManData->GetHistogram(Form("PairSEPP_%s",setting),"PairInvMass");
   THnF* sels2 = (THnF*)gHistManData->GetHistogram(Form("PairSEMM_%s",setting),"PairInvMass");
   THnF* mels1 = (THnF*)gHistManData->GetHistogram(Form("PairMEPP_%s",setting),"PairInvMass");
   THnF* mels2 = (THnF*)gHistManData->GetHistogram(Form("PairMEMM_%s",setting),"PairInvMass");
   
   THnF* seosMC = 0x0;
   TH1F* sigShape = 0x0;
   if(filenameMC[0]!='\0') {
      gHistManMC = new AliHistogramManager();
      gHistManMC->InitFile(filenameMC,"jpsi2eeHistos");
      seosMC = (THnF*)gHistManMC->GetHistogram(Form("PairSEPM_%s_TrueElectron",setting),"PairInvMass");
   }
   else {
      TFile* histFile=TFile::Open("/home/iarsene/work/ALICE/AliResonanceFits_bug/jpsiSignalShape.root");
      sigShape = (TH1F*)histFile->Get("Mass");
   }
   
   gFits.SetSEOSHistogram(seos);
   gFits.SetMEOSHistogram(meos);
   gFits.SetSELSHistograms(sels1, sels2);
   gFits.SetMELSHistograms(mels1, mels2);
   if(seosMC) gFits.SetSEOSMCHistogram(seosMC);
   else gFits.SetSignalMCshape(sigShape);
   
   //gFits.AddVariable(AliReducedVarManager::kCentVZERO, 2);
   //gFits.AddVariable(AliReducedVarManager::kVtxZ, 2);
   //gFits.SetVarRange(AliReducedVarManager::kCentVZERO, 0., 20.);
   gFits.AddVariable(AliReducedVarManager::kPt, 1);
   gFits.AddVariable(AliReducedVarManager::kMass, 0);
   gFits.SetVarRange(AliReducedVarManager::kMass, 1.0, 5.0);
   //gFits.SetVarRange(AliReducedVarManager::kPt, 10.0, 30.0);
   
   gFits.SetBkgMethod(AliResonanceFits::kBkgMixedEvent);          // kBkgMixedEvent,  kBkgLikeSign, kBkgMixedEventAndResidualFit, kBkgFitFunction
   gFits.SetMEMatchingMethod(AliResonanceFits::kMatchSEOS);     // kMatchSEOS, kMatchSELS
   gFits.SetLSmethod(AliResonanceFits::kLSArithmeticMean);    // kLSArithmeticMean, kLSGeometricMean
   gFits.SetUseRfactorCorrection(kTRUE);
   gFits.SetMassFitRange(1.5,4.2);
   gFits.SetMassExclusionRange(2.5, 3.2);
   gFits.SetScaleSummedBkg(kTRUE);
   //gFits.SetUseSignificantZero(kFALSE);
   gFits.SetScalingOption(AliResonanceFits::kScaleEntries);   // kScaleEntries, kScaleFit, kScaleWeightedAverage
   //gFits.SetWeightedAveragePower(2.0);                    // 2 - is default, but it can be modified in some specific situations
   //gFits.SetMinuitFitOption(AliResonanceFits::kMinuitMethodLikelihood);   // kMinuitMethodChi2, kMinuitMethodLikelihood
   //TF1* fitFunc = new TF1("fitFunc", "pol2", 0., 5.0);
   //gFits.SetResidualFitFunction(fitFunc);
   //gFits.SetDebugMode(kTRUE);
   //TF1* fitFunc = new TF1("fitFunc", "([0]*x^2+[1]*x+[2])/([3]*x^3+[4]*x^2+[5]*x+[6])", 1.0, 4.6);
   //TF1* fitFunc = new TF1("fitFunc", "[0]*(x-[1])*(x-[2])/pol3(3)", 0., 5.0);
   //TF1* fitFunc = new TF1("fitFunc", "pol2/pol3(3)", 0., 5.0);
   //fitFunc->SetParameters(1.63610e+04, 7.63652e+04, -1.09040e+04, -3.73453e+02, 1.32915e+03, -7.27885e+02, 1.21654e+02);
   //fitFunc->SetParameter(0, 0.); fitFunc->SetParLimits(0, -100., 0.);
   //fitFunc->SetParameter(1, 10.); fitFunc->SetParLimits(1, 5., 100.);
   TF1* fitFunc = new TF1("fitFunc", "pol1", 0., 5.0);
   //fitFunc->SetParameter(0, 10.); fitFunc->SetParameter(1,-0.1);
   gFits.SetBkgFitFunction(fitFunc);
   
   gFits.Process();
   gFits.ComputeOutputValues(2.92, 3.16);
   gFits.Print();
   gFits.PrintFitValues();
   
   Draw(&gFits, 2.92, 3.16, kTRUE, "testDir/", setting, isPreliminaryPlot);
   /*
   Double_t eff=0.0; Double_t sigFrac = 0.0;
   TH1D* effPt = ComputeEfficiency(0, filenameMC, setting, eff, sigFrac, 2.92, 3.16);
   cout << "eff = " << eff << endl;
   cout << "sigFrac = " << sigFrac << endl;
   
   TCanvas* c2=new TCanvas();
   effPt->Draw();
   */
}


//________________________________________________________________________________________
TH1D* ComputeEfficiency(Int_t option, const Char_t* filenameMC, const Char_t* setting, Double_t& eff, Double_t& sigFrac, 
                        Double_t minMass, Double_t maxMass, 
                        Double_t centMin=0., Double_t centMax=100.,
                        Double_t ptMin=0., Double_t ptMax=100.) {
   //
   // compute efficiency
   //
   AliHistogramManager* histMan=new AliHistogramManager();
   histMan->InitFile(filenameMC,"jpsi2eeHistos");
   THnF* nominatorMC = (THnF*)histMan->GetHistogram(Form("PairSEPM_%s_TrueElectron",setting),"PairInvMass");
   TH3F* denominatorMC = (TH3F*)histMan->GetHistogram("mcTruthJpsi_PureMCTruth_BeforeSelection","MassMC_Pt_CentVZERO");
   nominatorMC->Sumw2();
   //denominatorMC->Sumw2();
   
   /*nominatorMC->GetAxis(1)->SetRangeUser(ptMin+1.0e-6, ptMax-1.0e-6);
   denominatorMC->GetYaxis()->SetRangeUser(ptMin+1.0e-6, ptMax-1.0e-6);
   nominatorMC->GetAxis(2)->SetRangeUser(centMin+1.0e-6, centMax-1.0e-6);
   denominatorMC->GetZaxis()->SetRangeUser(centMin+1.0e-6, centMax-1.0e-6);
   */
   
   Double_t err=0.0;
   /*sigFrac = nominatorMC->IntegralAndError(nominatorMC->GetXaxis()->FindBin(minMass+1.0e-6), nominatorMC->GetXaxis()->FindBin(maxMass-1.0e-6), nominatorMC->GetYaxis()->FindBin(ptMin+1.0e-6), nominatorMC->GetYaxis()->FindBin(ptMax-1.0e-6), 
   nominatorMC->GetZaxis()->FindBin(centMin+1.0e-6), nominatorMC->GetZaxis()->FindBin(centMax-1.0e-6),err);
   
   sigFrac /= nominatorMC->IntegralAndError(1, nominatorMC->GetXaxis()->GetNbins(), 
                                            nominatorMC->GetYaxis()->FindBin(ptMin+1.0e-6), nominatorMC->GetYaxis()->FindBin(ptMax-1.0e-6), 
                                            nominatorMC->GetZaxis()->FindBin(centMin+1.0e-6), nominatorMC->GetZaxis()->FindBin(centMax-1.0e-6),err);
   */
   
   nominatorMC->GetAxis(0)->SetRangeUser(minMass+1.0e-6, maxMass-1.0e-6);
   
   TH1D* nomProjMass = nominatorMC->Projection(0); nomProjMass->SetName("nomProjMass");
   TH1D* denomProjMass = denominatorMC->ProjectionX("denomProjMass");
   TH1D* nomProjPt = nominatorMC->Projection(1); nomProjPt->SetName("nomProjPt");
   TH1D* denomProjPt = denominatorMC->ProjectionY("denomProjPt");
   TH1D* nomProjCent = nominatorMC->Projection(2); nomProjCent->SetName("nomProjCent");
   TH1D* denomProjCent = denominatorMC->ProjectionZ("denomProjCent");
   
   eff = nomProjPt->IntegralAndError(1, nomProjPt->GetXaxis()->GetNbins(), err);
   eff /= denomProjPt->IntegralAndError(1, nomProjPt->GetXaxis()->GetNbins(), err);
   
   nomProjPt->Divide(denomProjPt);
   nomProjCent->Divide(denomProjCent);
   if(option==0) return nomProjPt;
   if(option==1) return nomProjCent;
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
   
   TBox* boxL = new TBox(fits->GetMassFitRange()[0], histSplusB->GetMaximum()*0.1, fits->GetMassExclusionRange()[0], histSplusB->GetMaximum()*0.15);
   boxL->SetFillColorAlpha(4, 0.35);
   //boxL->SetFillStyle(3003);
   //boxL->Draw();
   
   TBox* boxR = new TBox(fits->GetMassExclusionRange()[1], histSplusB->GetMaximum()*0.1, fits->GetMassFitRange()[1], histSplusB->GetMaximum()*0.15);
   boxR->SetFillColorAlpha(4, 0.35);
   //boxR->SetFillStyle(3003);
   //boxR->Draw();
   
   TLegend* legend1 = new TLegend(0.18, 0.05, 0.42, 0.23);
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
   lat->DrawLatex(0.63, 0.60, Form("%.2f<m_{ee}<%.2f GeV/c^{2}", fFitRange[0], fFitRange[1]));
   lat->DrawLatex(0.63, 0.54, Form("%.2f<p_{T}<%.2f GeV/c", fPtFitRange[0], fPtFitRange[1]));
   lat->DrawLatex(0.63, 0.45, "Signal:");
   lat->DrawLatex(0.63, 0.39, Form("%.2f<m_{ee}<%.2f GeV/c^{2}", fSignalRange[0], fSignalRange[1]));
   lat->DrawLatex(0.63, 0.33, Form("%.2f<p_{T}<%.2f GeV/c", fPtRange[0], fPtRange[1])); */

   // Draw the bottom pad ================================================
   pad = c1->cd(2);
   pad->SetTopMargin(0.0);
   pad->SetRightMargin(0.02);
   pad->SetBottomMargin(0.17);
   pad->SetLeftMargin(0.15);
   pad->SetTicks(1,1);
   
   if(fits->GetBkgMethod()==AliResonanceFits::kBkgMixedEvent ||
      fits->GetBkgMethod()==AliResonanceFits::kBkgLikeSign ||
      fits->GetBkgMethod()==AliResonanceFits::kBkgFitFunction) {
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
      
      cout << "global fit " << globalFit << "  par0 " << globalFit->GetParameter(0) << endl;
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
TList* CheckPIDsystematicsForSetting(const Char_t* setting, const Char_t* refSetting="pid1") {
   //
   //
   //
   const Int_t kNrebinnedBins = 16;
   Double_t rebinEdges[kNrebinnedBins+1] = {
      0.0, 0.5, 1.0, 1.1, 1.2, 
      1.3, 1.4, 1.5, 1.6, 1.7, 
      1.8, 2.0, 2.5, 3.0, 4.5, 
      6.0, 10.0};
   TH1F* data_p_ref = (TH1F*)gHistManData->GetHistogram(Form("Track_%s", refSetting),"P");
   data_p_ref->SetName(Form("data_p_ref_%s_%.6f", refSetting, gRandom->Rndm()));
   TH1F* data_eta_ref = (TH1F*)gHistManData->GetHistogram(Form("Track_%s", refSetting),"Eta");
   data_eta_ref->SetName(Form("data_eta_ref_%s_%.6f", refSetting, gRandom->Rndm()));
   TH2F* data_etap_ref = (TH2F*)gHistManData->GetHistogram(Form("Track_%s", refSetting),"Eta_P");
   data_etap_ref->SetName(Form("data_etap_ref_%s_%.6f", refSetting, gRandom->Rndm()));
   TH1F* data_p_ref_rebinned = dynamic_cast<TH1F*>(data_p_ref->Rebin(kNrebinnedBins, Form("data_p_ref_%s_rebinned", refSetting), rebinEdges));
   data_eta_ref->Rebin(2);
   TH1F* mc_p_ref = (TH1F*)gHistManMC->GetHistogram(Form("Track_%s_TrueElectron", refSetting),"P");
   mc_p_ref->SetName(Form("mc_p_ref_%s_%.6f", refSetting, gRandom->Rndm()));
   TH1F* mc_eta_ref = (TH1F*)gHistManMC->GetHistogram(Form("Track_%s_TrueElectron", refSetting),"Eta");
   mc_eta_ref->SetName(Form("mc_eta_ref_%s_%.6f", refSetting, gRandom->Rndm()));
   TH2F* mc_etap_ref = (TH2F*)gHistManMC->GetHistogram(Form("Track_%s_TrueElectron", refSetting),"Eta_P");
   mc_etap_ref->SetName(Form("mc_etap_ref_%s_%.6f", refSetting, gRandom->Rndm()));
   TH1F* mc_p_ref_rebinned = dynamic_cast<TH1F*>(mc_p_ref->Rebin(kNrebinnedBins, Form("mc_p_ref_%s_rebinned", refSetting), rebinEdges));
   mc_eta_ref->Rebin(2);
   
   TH1F* data_p = (TH1F*)gHistManData->GetHistogram(Form("Track_%s", setting),"P");
   data_p->SetName(Form("data_p_%s", setting));
   TH1F* data_eta = (TH1F*)gHistManData->GetHistogram(Form("Track_%s", setting),"Eta");
   data_eta->SetName(Form("data_eta_%s", setting));
   TH2F* data_etap = (TH2F*)gHistManData->GetHistogram(Form("Track_%s", setting),"Eta_P");
   data_etap->SetName(Form("data_etap_%s", setting));
   TH1F* data_p_rebinned = dynamic_cast<TH1F*>(data_p->Rebin(kNrebinnedBins, Form("data_p_%s_rebinned", refSetting), rebinEdges));
   data_eta->Rebin(2);
   
   TH1F* mc_p = (TH1F*)gHistManMC->GetHistogram(Form("Track_%s_TrueElectron", setting),"P");
   mc_p->SetName(Form("mc_p_%s", setting));
   TH1F* mc_eta = (TH1F*)gHistManMC->GetHistogram(Form("Track_%s_TrueElectron", setting),"Eta");
   mc_eta->SetName(Form("mc_eta_%s", setting));
   TH2F* mc_etap = (TH2F*)gHistManMC->GetHistogram(Form("Track_%s_TrueElectron", setting),"Eta_P");
   mc_etap->SetName(Form("mc_etap_%s", setting));
   TH1F* mc_p_rebinned = dynamic_cast<TH1F*>(mc_p->Rebin(kNrebinnedBins, Form("mc_p_%s_rebinned", refSetting), rebinEdges));
   mc_eta->Rebin(2);
   
   TH1F* data_ratio_p = (TH1F*) data_p_rebinned->Clone(Form("data_ratio_p_%s", setting));
   BeutifyTH1(data_ratio_p, "", 2, 4);
   TH1F* mc_ratio_p = (TH1F*) mc_p_rebinned->Clone(Form("mc_ratio_p_%s", setting));
   BeutifyTH1(mc_ratio_p, "", 2, 2);
   data_ratio_p->Divide(data_p_rebinned, data_p_ref_rebinned, 1.,1.,"B");
   mc_ratio_p->Divide(mc_p_rebinned, mc_p_ref_rebinned, 1.,1.,"B");
   TH2F* data_ratio_etap = (TH2F*) data_etap->Clone(Form("data_ratio_etap_%s", setting));
   TH2F* mc_ratio_etap = (TH2F*) mc_etap->Clone(Form("mc_ratio_etap_%s", setting));
   data_ratio_etap->Divide(data_etap, data_etap_ref, 1.,1.,"B");
   mc_ratio_etap->Divide(mc_etap, mc_etap_ref, 1.,1.,"B");
   
   TH1F* dataMC_ratio_p = (TH1F*)data_ratio_p->Clone(Form("dataMC_ratio_p_%s", setting));
   BeutifyTH1(dataMC_ratio_p, "", 2, 6);
   dataMC_ratio_p->Divide(mc_ratio_p);
   TH2F* dataMC_ratio_etap = (TH2F*)data_ratio_etap->Clone(Form("dataMC_ratio_etap_%s", setting));
   dataMC_ratio_etap->Divide(mc_ratio_etap);
   
   TList* returnItems = new TList();
   returnItems->Add(data_ratio_p); returnItems->Add(mc_ratio_p);
   returnItems->Add(data_ratio_etap); returnItems->Add(mc_ratio_etap);
   returnItems->Add(dataMC_ratio_p);
   returnItems->Add(dataMC_ratio_etap);
   
   return returnItems;
}


//____________________________________________________________________________________________________
void CheckPIDsystematics(const Char_t* dataFilename, const Char_t* mcFilename,
      const Char_t* settingName, const Char_t* referenceSettingName,
      const Char_t* outputFilename="pidSyst_SingleEle.root") {
   //
   // 
   //
   gHistManData = new AliHistogramManager();
   gHistManData->InitFile(dataFilename, "jpsi2eeHistos");
   
   gHistManMC = new AliHistogramManager();
   gHistManMC->InitFile(mcFilename, "jpsi2eeHistos");
   
   //TList* electron28Ratios = CheckPIDsystematicsForSetting("electron28","electron30");
   //TList* electron26Ratios = CheckPIDsystematicsForSetting("electron26","electron30");
   //TList* proton42Ratios = CheckPIDsystematicsForSetting("proton42","proton36");
   //TList* proton40Ratios = CheckPIDsystematicsForSetting("proton40","proton36");
   //TList* proton38Ratios = CheckPIDsystematicsForSetting("proton38","proton36");
   //TList* pion40Ratios = CheckPIDsystematicsForSetting("pion40","pion34");
   //TList* pion38Ratios = CheckPIDsystematicsForSetting("pion38","pion34");
   //TList* pion36Ratios = CheckPIDsystematicsForSetting("pion36","pion34");
   //TList* electron28Ratios_uncalibrated = CheckPIDsystematicsForSetting("electron28_uncalibrated","electron30_uncalibrated");
   //TList* electron26Ratios_uncalibrated = CheckPIDsystematicsForSetting("electron26_uncalibrated","electron30_uncalibrated");
   //TList* proton42Ratios_uncalibrated = CheckPIDsystematicsForSetting("proton42_uncalibrated","proton36_uncalibrated");
   //TList* proton40Ratios_uncalibrated = CheckPIDsystematicsForSetting("proton40_uncalibrated","proton36_uncalibrated");
   //TList* proton38Ratios_uncalibrated = CheckPIDsystematicsForSetting("proton38_uncalibrated","proton36_uncalibrated");
   //TList* pion40Ratios_uncalibrated = CheckPIDsystematicsForSetting("pion40_uncalibrated","pion34_uncalibrated");
   //TList* pion38Ratios_uncalibrated = CheckPIDsystematicsForSetting("pion38_uncalibrated","pion34_uncalibrated");
   //TList* pion36Ratios_uncalibrated = CheckPIDsystematicsForSetting("pion36_uncalibrated","pion34_uncalibrated");
   TList* histList = CheckPIDsystematicsForSetting(settingName, referenceSettingName);
   
   
   TH1F* ptMC = (TH1F*)gHistManMC->GetHistogram("mcTruthJpsi_PureMCTruth_BeforeSelection", "PtMC");
   ptMC->SetName("JpsiPtSpectrum");
   
   TFile* save=new TFile(outputFilename, "RECREATE");
   ptMC->Write();
   //electron28Ratios->Write();
   //electron26Ratios->Write();
   //proton42Ratios->Write();
   //proton40Ratios->Write();
   //proton38Ratios->Write();
   //pion40Ratios->Write();
   //pion38Ratios->Write();
   //pion36Ratios->Write();
   //electron28Ratios_uncalibrated->Write();
   //electron26Ratios_uncalibrated->Write();
   //proton42Ratios_uncalibrated->Write();
   //proton40Ratios_uncalibrated->Write();
   //proton38Ratios_uncalibrated->Write();
   //pion40Ratios_uncalibrated->Write();
   //pion38Ratios_uncalibrated->Write();
   //pion36Ratios_uncalibrated->Write();
   histList->Write();
   
   save->Close();
}

//____________________________________________________________________________________________________
void ComputeInputKineSyst(const Char_t* mcFilename, const Char_t* setting) {
   //
   // computes the systematic uncertainty on the choice of the input MC kinematics
   //
   gHistManMC = new AliHistogramManager();
   gHistManMC->InitFile(mcFilename, "jpsi2eeHistos");
   
   TH1F* jpsiPtSpectrumInputMC_coarse = (TH1F*)gHistManMC->GetHistogram("mcTruthJpsi_PureMCTruth_BeforeSelection", "PtMC_coarse");
   TH1F* jpsiPtSpectrumInputMC = (TH1F*)gHistManMC->GetHistogram("mcTruthJpsi_PureMCTruth_BeforeSelection", "PtMC");
   jpsiPtSpectrumInputMC->Rebin(10);
   //TH3F* jpsiInput3D = (TH3F*)gHistManMC->GetHistogram("mcTruthJpsi_PureMCTruth_BeforeSelection", "MassMC_Pt_CentVZERO");
   
   
   //THnF* selectedJpsi = (THnF*)gHistManMC->GetHistogram(Form("PairSEPM_%s_TrueElectron", setting), "PairInvMass");
   //selectedJpsi->GetAxis(0)->SetRangeUser(2.9201, 3.159);
   //TH1D* recYield = (TH1D*)selectedJpsi->Projection(0);
   TH1F* recYield_coarse = (TH1F*)gHistManMC->GetHistogram(Form("PairSEPM_%s_TrueElectron", setting), "Pt_coarse");
   TH1F* recYield = (TH1F*)gHistManMC->GetHistogram(Form("PairSEPM_%s_TrueElectron", setting), "Pt");
   recYield->Rebin(10);
   
   TH1F* effPt = (TH1F*)recYield->Clone("effPt");
   effPt->Divide(recYield, jpsiPtSpectrumInputMC, 1., 1., "B");
   
   TH1F* effPt_coarse = (TH1F*)recYield_coarse->Clone("effPt_coarse");
   effPt_coarse->Divide(recYield_coarse, jpsiPtSpectrumInputMC_coarse, 1., 1., "B");
   
   //effPt->Draw();
   //effPt_coarse->Draw("same");
   
   // get the spectra file from the Pb-Pb analysis
   TFile* fileSpectraTona = TFile::Open("analysisOutputs/StuffForAN/dNdpT_Cent0_20_40_90_24062017.root");
   TH1D* spectraPbPb0_20 = (TH1D*)fileSpectraTona->Get("totalYieldPt020");
   TH1D* spectraPbPb20_40 = (TH1D*)fileSpectraTona->Get("totalYieldPt2040");
   TH1D* spectraPbPb40_90 = (TH1D*)fileSpectraTona->Get("totalYieldPt4090");
   
   TFile* fileMUON = TFile::Open("analysisOutputs/StuffForAN/SpectraPbPbFwdY.root");
   TH1D* spectraJpsi0 = (TH1D*)fileMUON->Get("hPsiAccCor_Cent0");
   TF1* spectraFit = (TF1*)spectraJpsi0->GetListOfFunctions()->At(0);
   TH1D* spectraJpsi1 = (TH1D*)fileMUON->Get("hPsiAccCor_Cent1");
   TH1D* spectraJpsi2 = (TH1D*)fileMUON->Get("hPsiAccCor_Cent2");
   TH1D* spectraJpsi3 = (TH1D*)fileMUON->Get("hPsiAccCor_Cent3");
   TH1D* spectraJpsi4 = (TH1D*)fileMUON->Get("hPsiAccCor_Cent4");
   TH1D* spectraJpsi5 = (TH1D*)fileMUON->Get("hPsiAccCor_Cent5");
   TH1D* spectraJpsi6 = (TH1D*)fileMUON->Get("hPsiAccCor_Cent6");
   
   jpsiPtSpectrumInputMC->Fit(spectraFit, "ME", "Q", 0., 10.);
   
   Double_t efff = WeightedAverage(effPt, spectraFit, 0.0, 10.0);
   cout << "efff = " << efff << endl;
   return;
   //jpsiPtSpectrumInputMC->Draw();
   //spectraFit->Draw("same");
   //cout << spectraFit->GetExpFormula().Data() << endl;
   //return;
   
   TCanvas* c1=new TCanvas();   
   
   spectraPbPb20_40->Fit(spectraFit, "ME","Q",0.,10.);
   spectraFit->SetLineStyle(3);
   spectraFit->SetLineColor(4);
   spectraPbPb0_20->Draw();
   
   /*TCanvas* cc = new TCanvas();
   gMinuit->SetErrorDef(1); // 1 sigmas (2^2)
   TGraph* g12 = static_cast<TGraph*>(gMinuit->Contour(200,1,2));
   g12->Draw("al");
   g12->GetXaxis()->SetTitle("pt0");
   g12->GetYaxis()->SetTitle("n");
   cc->Update();*/
   //return;
   
   
   Double_t meanPar1 = spectraFit->GetParameter(1);
   Double_t par1Err = spectraFit->GetParError(1);
   Double_t meanPar2 = spectraFit->GetParameter(2);
   Double_t par2Err = spectraFit->GetParError(2);
   Double_t meanPar3 = spectraFit->GetParameter(3);
   Double_t par3Err = spectraFit->GetParError(3);
   
   cout << "spectra fit function :: " << spectraFit->GetExpFormula() << endl;
   cout << "par1  :: " << meanPar1 << " +/- " << par1Err << endl;
   cout << "par2  :: " << meanPar2 << " +/- " << par2Err << endl;
   cout << "par3  :: " << meanPar3 << " +/- " << par3Err << endl;
   cout << "chi2 / ndf :: " << spectraFit->GetChisquare() << " / " << spectraFit->GetNDF() << " = " << spectraFit->GetChisquare() / spectraFit->GetNDF() << endl;
   
   c1->cd();
   TH1D* hEff = new TH1D("hEff", "efficiency after MC input kine variations", 1000, 0.03, 0.15);
   TH1D* hInvEff = new TH1D("hInvEff", "1.0/efficiency after MC input kine variations", 1000, 5.0, 15.);
   spectraJpsi0->SetLineColor(2);
   spectraJpsi0->SetMarkerColor(2);
   //spectraJpsi0->Draw();
   for(Int_t i=0; i<1000;++i) {
      Double_t par1 = gRandom->Gaus(meanPar1, par1Err);
      //Double_t par2 = gRandom->Gaus(-0.18 + par1*1.1, 0.3);              // 0-20%
      Double_t par2 = gRandom->Gaus(0.54 + par1*0.77, 0.2);             // 20-40%
      //Double_t par2 = gRandom->Gaus(0.65 + par1*0.636, 0.2);             // 40-90%
      spectraFit->SetParameter(1, par1);
      spectraFit->SetParameter(2, par2);
      //spectraFit->SetParameter(1, gRandom->Gaus(meanPar1, par1Err));
      //spectraFit->SetParameter(2, gRandom->Gaus(meanPar2, par2Err));
      //spectraFit->SetParameter(3, gRandom->Gaus(meanPar3, par3Err));
      Double_t eff = WeightedAverage(effPt, spectraFit, 0.0, 10.0);
      hEff->Fill(eff);
      hInvEff->Fill(1.0/eff);
      spectraFit->DrawClone("same");
   }
   spectraPbPb20_40->Draw("same");
      
   TCanvas* c2=new TCanvas();
   hEff->Draw();
   TCanvas* c3=new TCanvas();
   hInvEff->Draw();
   TCanvas* c4=new TCanvas();
   //effTona->Draw();
   effPt->Draw();
}

//____________________________________________________________________________
Double_t WeightedAverage(TH1F* eff, TF1* weight, Float_t ptMin, Float_t ptMax) {
   //
   // Takes as arguments the efficiency as a funciton of pt and the input pt spectrum.
   // Calculates the average efficiency in the pt range ptMin - ptMax
   //
   Double_t average = 0.0;
   Double_t norm = 0.0;
   for(Int_t i=1; i<=eff->GetXaxis()->GetNbins(); ++i) {
      Double_t pt=eff->GetXaxis()->GetBinCenter(i);
      if(pt<ptMin) continue;
      if(pt>ptMax) continue;
      average += eff->GetBinContent(i) * weight->Eval(pt) * eff->GetXaxis()->GetBinWidth(i);
      norm += weight->Eval(pt) * eff->GetXaxis()->GetBinWidth(i);
   }
   average /= norm;
   return average;
}

//____________________________________________________________________________
Double_t WeightedAverage(TH2F* eff, TF1* weightPt, TF1* weightPolariz, Float_t ptMin, Float_t ptMax) {
   //
   // Takes as arguments a 2-D efficiency as a function of pt and cos(theta*), the pt spectrum and the cos(theta*) distribution
   // Calculates the average efficiency between ptMin and ptMax
   //
   Double_t average = 0.0;
   Double_t norm = 0.0;
   for(Int_t ipt=1; ipt<=eff->GetYaxis()->GetNbins(); ++ipt) {
      Double_t pt=eff->GetYaxis()->GetBinCenter(ipt);
      if(pt<ptMin) continue;
      if(pt>ptMax) continue;
      for(Int_t ict=1; ict<=eff->GetXaxis()->GetNbins(); ++ict) {
         Double_t cTheta=eff->GetXaxis()->GetBinCenter(ict);
         average += eff->GetBinContent(ict,ipt) * weightPt->Eval(pt) * weightPolariz->Eval(cTheta);
         norm += weightPt->Eval(pt) * weightPolariz->Eval(cTheta);
      }
   }
   average /= norm;
   return average;
}


//________________________________________________________________________________________
void ConsistencyCheck() {
   
   Double_t preliminary2040 = 0.7389;
   Double_t preliminary2040stat = 0.065564;
   Double_t preliminary2040syst = 0.05715299;
   
   Double_t current2040 = 0.827;
   Double_t current2040stat = 0.072;
   Double_t current2040syst = 0.065;
   
   Double_t xexe090 = 1.315;
   Double_t xexe090stat = 0.37;
   Double_t xexe090syst = 0.122;
   
   Double_t xexefwd090 = 0.54;
   Double_t xexefwd090stat = 0.11;
   Double_t xexefwd090syst = 0.09;
   
   Double_t refsystRel = 0.17;
   
   // difference wrt preliminary result: 1.36 sigma
   // difference wrt updated result: 1.14 sigma
   // difference wrt XeXe at fwd-y:  1.64 sigma
   // difference wrt XeXe at fwd-y (shifted to Pb-Pb value):  1.4 sigma
   
   TH2F* range = new TH2F("range", "", 10, 0.0, 3.0, 10, 0., 1.);
   range->Draw();
   
   TF1* xexePDFnoRef = new TF1("xexePDFnoRef", "gaus", 0., 3.);
   xexePDFnoRef->SetParameters(1., xexe090, TMath::Sqrt(xexe090stat*xexe090stat + xexe090syst*xexe090syst));
   
   TF1* pbpbPrelPDFnoRef = new TF1("pbpbPrelPDFnoRef", "gaus", 0., 3.);
   pbpbPrelPDFnoRef->SetParameters(1., preliminary2040, TMath::Sqrt(preliminary2040stat*preliminary2040stat + preliminary2040syst*preliminary2040syst));
   
   TH1D* diff = new TH1D("diff", "diff", 600, -2., 4.);
   for(Int_t i=0; i<1000000; ++i) {
      diff->Fill(xexePDFnoRef->GetRandom() - pbpbPrelPDFnoRef->GetRandom());
   }
   diff->Draw();
   
   
   /*xexePDFnoRef->Draw("same");
   pbpbPrelPDFnoRef->Draw("same");*/
}