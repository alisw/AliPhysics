void ExtractJpsiSignal() {
   //
   // Extract jpsi signal
   //
   const Char_t* filename = "/home/iarsene/work/ALICE/treeAnalysis/newdst/development/test/dstAnalysisHistograms.root";
   AliHistogramManager* histMan=new AliHistogramManager();
   histMan->InitFile(filename,"jpsi2eeHistos");
   
   THnF* seos = (THnF*)histMan->GetHistogram("PairSEPM_Pt10","PairInvMass");
   THnF* meos = (THnF*)histMan->GetHistogram("PairMEPM_Pt10","PairInvMass");
   THnF* sepp = (THnF*)histMan->GetHistogram("PairSEPP_Pt10","PairInvMass");
   THnF* semm = (THnF*)histMan->GetHistogram("PairSEMM_Pt10","PairInvMass");
   THnF* mepp = (THnF*)histMan->GetHistogram("PairMEPP_Pt10","PairInvMass");
   THnF* memm = (THnF*)histMan->GetHistogram("PairMEMM_Pt10","PairInvMass");
   
   //TH2F* h2 = (TH2F*)(histMan->GetHistogram("Track_Pt10", "TPCsignal_Pin"))->Clone("cloneHist");
   //h2->DrawClone();
   seos->GetAxis(2)->SetRangeUser(0.1,9.9);
   //seos->GetAxis(0)->SetRangeUser(2.9,3.19);
   seos->Projection(0)->Draw();
   return;
   
   
  AliResonanceFits* jpsiFit = new AliResonanceFits();
  jpsiFit->SetHistograms(seos, meos, sepp, semm, mepp, memm); 
  
  Int_t variables[AliResonanceFits::kNVariables] = {0, 1, -1, 2, 3, 4};
  jpsiFit->SetVars(variables);
  
  jpsiFit->SetFitRange(2.0, 4.0);
  jpsiFit->SetExclusionRange(2.5,3.5);
  jpsiFit->SetMassRange(0.0,5.0);
  jpsiFit->SetSignalRange(2.92,3.16);
  jpsiFit->SetCentralityRange(0.0,10.0);
  
  jpsiFit->SetBkgMethod(1);        // 1-mixed event; 2-like sign
  jpsiFit->SetUse2DMatching(kFALSE);
  jpsiFit->SetMatchingOption(1);       // 1-weighted average, 2-fit, 3-entries
  jpsiFit->SetWeightedAveragePower(2.);
  jpsiFit->SetMinuitFitOption(1.);    // 1.-chi2 minimization; 0.5-log likelihood
  jpsiFit->SetFixBkgScale(kFALSE);    // fix the bkg scale in Minuit (sometimes needed just to compute chi2)
  jpsiFit->SetMEMatchingBkg(2);     // 1-match bkg to the SEOS outside signal region
                                                        // 2-match to the SELS R-factor corrected bkg  (default)
  jpsiFit->SetLSmethod(2);      // 1-arithmetic mean; 2-geometric mean
  jpsiFit->SetPlotingOption(0);   // 0-mass projection, 1-pt projection, 2-(mass,pt projection)
  
  //jpsiFit->DrawSignalExtraction(kFALSE, "", "", kFALSE, 0x0, kFALSE, kFALSE);
}