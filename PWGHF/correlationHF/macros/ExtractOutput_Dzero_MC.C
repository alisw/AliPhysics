//
// On AliRoot shell, call the following before loading the macro:
//
// gSystem->SetIncludePath("-I. -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include -I$ALICE_ROOT/ANALYSIS/macros -I$ROOTSYS/include");
// gROOT->LoadMacro("AliDhCorrelationExtraction.cxx++g"); // set the right path! - if not already in AliPhysics!!
//
// Modify:
// - The arguments of the macro
// - The names and paths in SetInputNames method
// - The SetSBRanges values in the ExtractXPt methods (if you use austoSB=kFALSE, otherwise they are dummy)
//
// NOTE FOR D*: If setting the sideband range externally, the (unique) sideband to be used is LSB (not RSB)
// NOTE FOR D+: If you need to integrate mass-pT bins in a single correlation bin BEFORE extracting the outputs, set plotter->IntegratePtBins(kFALSE) in the ExtractXPt methods 
// 
// setting plotter->SetDebugLevel(1) will print a series of additional plots for debugging; plotter->SetDebugLevel(2) will also be more verbose
//


void ExtractOutput_MC(
   Int_t specie=AliDhCorrelationExtraction::kD0toKpi, //the D-meson decay channel (check the enumerator for the options)
   Bool_t singleBinSB=kTRUE, //kTRUE=a single mass bin is used in the THnSparse for storing the sideband corlelations (used by D0 to save space)
   Int_t npools=9, //number of pools for the event-mixing
   Bool_t poolByPool=kTRUE, //kTRUE=pool-by-pool ME correction; kFALSE=merged-pools ME correction (set the options that you used in the online analysis)
   Double_t deltaEtaMin=-1., Double_t deltaEtaMax=1., //deltaEta ranges for correlation distributions  
   Int_t recomode=AliDhCorrelationExtraction::kReco) //kReco = reco, kKine = kine
{

  //Create and set the correlation plotter class
  AliDhCorrelationExtraction *plotter = new AliDhCorrelationExtraction();
  Bool_t flagSpecie = plotter->SetDmesonSpecie(specie);
  plotter->SetNpools(npools);
  plotter->SetCorrectPoolsSeparately(poolByPool);
  plotter->SetDeltaEtaRange(deltaEtaMin,deltaEtaMax);
  plotter->SetRecoMode(recomode);
  if(!flagSpecie) return;

  plotter->SetDebugLevel(2); //0 = get main results; 1 = get full list of plots; 2 = get debug printouts

  SetInputNames(plotter,recomode);  // check the names in the method!!

  Bool_t read = plotter->ReadInputs();
  if(!read) {
    printf("Error in reading the input file! Exiting...\n");
    return;
  }

//  gSystem->Exec(Form("rm -r Output_Root"));
//  gSystem->Exec(Form("rm -r Output_png")); 
  gSystem->Exec(Form("mkdir Output_Root"));
  gSystem->Exec(Form("mkdir Output_png"));

  ExtractLowPt_MC(plotter);
  ExtractMidPt_MC(plotter);
  ExtractHighPt_MC(plotter); 
  ExtractVeryHighPt_MC(plotter);

  return;
}

//________________________________________
void SetInputNames(AliDhCorrelationExtraction *plotter, Int_t recomode) {

  TString suffMode;
  if(recomode==AliDhCorrelationExtraction::kKine) suffMode = "Kine";
  if(recomode==AliDhCorrelationExtraction::kReco) suffMode = "Reco";

  //D0 filenames and paths
  plotter->SetInputFilenameMass("./AnalysisResults_MCclosure_AllIngr.root");
  plotter->SetInputFilenameSE("./AnalysisResults_MCclosure_AllIngr.root");
  plotter->SetInputFilenameME("./AnalysisResults_MCclosure_AllIngr.root");
  plotter->SetDirNameMass(Form("D0hCorrelOutput_pPb_0_100_%sSE",suffMode.Data()));
  plotter->SetListNameMass(Form("coutputmassD0MassOutput_pPb_0_100_%sSE_0100",suffMode.Data()));
  plotter->SetDirNameSE(Form("D0hCorrelOutput_pPb_0_100_%sSE",suffMode.Data()));
  plotter->SetListNameSE(Form("correlationsOutput_pPb_0_100_%sSE_0100",suffMode.Data()));
  plotter->SetDirNameME(Form("D0hCorrelOutput_pPb_0_100_%sME",suffMode.Data()));
  plotter->SetListNameME(Form("correlationsOutput_pPb_0_100_%sME_0100",suffMode.Data()));

  plotter->SetMassHistoName("histSgn_");
  if(recomode==AliDhCorrelationExtraction::kReco) plotter->SetMassHistoName("histSgn_WeigD0Eff_");
  plotter->SetSECorrelHistoName("hPhi_Charg");
  plotter->SetMECorrelHistoNameSuffix("_EvMix");
  
  plotter->AddOriginType("",AliDhCorrelationExtraction::kOrigAll);
  plotter->AddOriginType("_HF_From_c",AliDhCorrelationExtraction::kOrigC);
  plotter->AddOriginType("_HF_From_b",AliDhCorrelationExtraction::kOrigB);
  plotter->AddOriginType("_NonHF",AliDhCorrelationExtraction::kOrigLF);
  plotter->AddOriginType("_From_c",AliDhCorrelationExtraction::kOrigC);
  plotter->AddOriginType("_From_b",AliDhCorrelationExtraction::kOrigB);

  return;
}

//________________________________________
void ExtractLowPt_MC(AliDhCorrelationExtraction *plotter) {

  plotter->SetNpTbins(2);
  plotter->SetFirstpTbin(4);
  plotter->SetLastpTbin(5);
  
  //extract correlation distributions
  printf("*** Extracting correlations in 3<pT(D)<5 GeV/c, 0.3<pT(assoc)<99 GeV/c ***\n");
  Bool_t corrExtrThr1 = plotter->ExtractCorrelations_MC(0.3,99.);
  printf("*** Extracting correlations in 3<pT(D)<5 GeV/c, 0.3<pT(assoc)<1 GeV/c ***\n");
  Bool_t corrExtrThr2 = plotter->ExtractCorrelations_MC(0.3,1.);
  printf("*** Extracting correlations in 3<pT(D)<5 GeV/c, 1<pT(assoc)<99 GeV/c ***\n");
  Bool_t corrExtrThr3 = plotter->ExtractCorrelations_MC(1.,99.);
  printf("*** Extracting correlations in 3<pT(D)<5 GeV/c, 1<pT(assoc)<2 GeV/c ***\n");
  Bool_t corrExtrThr4 = plotter->ExtractCorrelations_MC(1.,2.);
  printf("*** Extracting correlations in 3<pT(D)<5 GeV/c, 2<pT(assoc)<3 GeV/c ***\n");
  Bool_t corrExtrThr5 = plotter->ExtractCorrelations_MC(2.,3);
  printf("*** Extracting correlations in 3<pT(D)<5 GeV/c, 3<pT(assoc)<99 GeV/c ***\n");
  Bool_t corrExtrThr6 = plotter->ExtractCorrelations_MC(3.,99.);
/*  printf("*** Extracting correlations in 3<pT(D)<5 GeV/c, 2<pT(assoc)<99 GeV/c ***\n");
  Bool_t corrExtrThr4 = plotter->ExtractCorrelations_MC(2.,99.);
  printf("*** Extracting correlations in 3<pT(D)<5 GeV/c, 3<pT(assoc)<99 GeV/c ***\n");
  Bool_t corrExtrThr5 = plotter->ExtractCorrelations_MC(3.,99.);    
  printf("*** Extracting correlations in 3<pT(D)<5 GeV/c, 4<pT(assoc)<99 GeV/c ***\n");
  Bool_t corrExtrThr6 = plotter->ExtractCorrelations_MC(4.,99.);
  printf("*** Extracting correlations in 3<pT(D)<5 GeV/c, 1<pT(assoc)<2 GeV/c ***\n");
  Bool_t corrExtrThr7 = plotter->ExtractCorrelations_MC(1.,2.);
  printf("*** Extracting correlations in 3<pT(D)<5 GeV/c, 1<pT(assoc)<3 GeV/c ***\n");
  Bool_t corrExtrThr8 = plotter->ExtractCorrelations_MC(1.,3);
  printf("*** Extracting correlations in 3<pT(D)<5 GeV/c, 2<pT(assoc)<3 GeV/c ***\n");
  Bool_t corrExtrThr9 = plotter->ExtractCorrelations_MC(2.,3);*/
  if(!corrExtrThr1 || !corrExtrThr2 || !corrExtrThr3 /*|| !corrExtrThr4 || !corrExtrThr5 || !corrExtrThr6 || !corrExtrThr7 || !corrExtrThr8 || !corrExtrThr9*/) {
    printf("Error in the extraction of the correlation distributions! Exiting...\n");
    return;
  }

  plotter->ClearObjects(); //important! Call it after each wide-pT range

}


//________________________________________
void ExtractMidPt_MC(AliDhCorrelationExtraction *plotter) {

  plotter->SetNpTbins(3);
  plotter->SetFirstpTbin(6);
  plotter->SetLastpTbin(8);

  //extract correlation distributions
  printf("*** Extracting correlations in 5<pT(D)<8 GeV/c, 0.3<pT(assoc)<99 GeV/c ***\n");
  Bool_t corrExtrThr1 = plotter->ExtractCorrelations_MC(0.3,99.);
  printf("*** Extracting correlations in 5<pT(D)<8 GeV/c, 0.3<pT(assoc)<1 GeV/c ***\n");
  Bool_t corrExtrThr2 = plotter->ExtractCorrelations_MC(0.3,1.);
  printf("*** Extracting correlations in 5<pT(D)<8 GeV/c, 1<pT(assoc)<99 GeV/c ***\n");
  Bool_t corrExtrThr3 = plotter->ExtractCorrelations_MC(1.,99.);
  printf("*** Extracting correlations in 5<pT(D)<8 GeV/c, 1<pT(assoc)<2 GeV/c ***\n");
  Bool_t corrExtrThr4 = plotter->ExtractCorrelations_MC(1.,2.);
  printf("*** Extracting correlations in 5<pT(D)<8 GeV/c, 2<pT(assoc)<3 GeV/c ***\n");
  Bool_t corrExtrThr5 = plotter->ExtractCorrelations_MC(2.,3);
  printf("*** Extracting correlations in 5<pT(D)<8 GeV/c, 3<pT(assoc)<99 GeV/c ***\n");
  Bool_t corrExtrThr6 = plotter->ExtractCorrelations_MC(3.,99.);
/*  printf("*** Extracting correlations in 5<pT(D)<8 GeV/c, 2<pT(assoc)<99 GeV/c ***\n");
  Bool_t corrExtrThr4 = plotter->ExtractCorrelations_MC(2.,99.);
  printf("*** Extracting correlations in 5<pT(D)<8 GeV/c, 3<pT(assoc)<99 GeV/c ***\n");
  Bool_t corrExtrThr5 = plotter->ExtractCorrelations_MC(3.,99.);    
  printf("*** Extracting correlations in 5<pT(D)<8 GeV/c, 4<pT(assoc)<99 GeV/c ***\n");
  Bool_t corrExtrThr6 = plotter->ExtractCorrelations_MC(4.,99.);
  printf("*** Extracting correlations in 5<pT(D)<8 GeV/c, 1<pT(assoc)<2 GeV/c ***\n");
  Bool_t corrExtrThr7 = plotter->ExtractCorrelations_MC(1.,2.);
  printf("*** Extracting correlations in 5<pT(D)<8 GeV/c, 1<pT(assoc)<3 GeV/c ***\n");
  Bool_t corrExtrThr8 = plotter->ExtractCorrelations_MC(1.,3.);
  printf("*** Extracting correlations in 5<pT(D)<8 GeV/c, 2<pT(assoc)<3 GeV/c ***\n");
  Bool_t corrExtrThr9 = plotter->ExtractCorrelations_MC(2.,3.); */
  if(!corrExtrThr1 || !corrExtrThr2 || !corrExtrThr3 /*|| !corrExtrThr4 || !corrExtrThr5 || !corrExtrThr6 || !corrExtrThr7 || !corrExtrThr8 || !corrExtrThr9*/) {
    printf("Error in the extraction of the correlation distributions! Exiting...\n");
    return;
  }

  plotter->ClearObjects(); //important! Call it after each wide-pT range

}


//________________________________________
void ExtractHighPt_MC(AliDhCorrelationExtraction *plotter) {
	
  plotter->SetNpTbins(2);
  plotter->SetFirstpTbin(9);
  plotter->SetLastpTbin(10);

  plotter->IntegratePtBins(kFALSE); //For D+, set it kTRUE; High pT bin (and only it) has to be evaluated merging mass plots and THnSparses for the 5 pT(D+) bins!

  //extract correlation distributions
  printf("*** Extracting correlations in 8<pT(D)<16 GeV/c, 0.3<pT(assoc)<99 GeV/c ***\n");
  Bool_t corrExtrThr1 = plotter->ExtractCorrelations_MC(0.3,99.);
  printf("*** Extracting correlations in 8<pT(D)<16 GeV/c, 0.3<pT(assoc)<1 GeV/c ***\n");
  Bool_t corrExtrThr2 = plotter->ExtractCorrelations_MC(0.3,1.);
  printf("*** Extracting correlations in 8<pT(D)<16 GeV/c, 1<pT(assoc)<99 GeV/c ***\n");
  Bool_t corrExtrThr3 = plotter->ExtractCorrelations_MC(1.,99.);
  printf("*** Extracting correlations in 8<pT(D)<16 GeV/c, 1<pT(assoc)<2 GeV/c ***\n");
  Bool_t corrExtrThr4 = plotter->ExtractCorrelations_MC(1.,2.);
  printf("*** Extracting correlations in 8<pT(D)<16 GeV/c, 2<pT(assoc)<3 GeV/c ***\n");
  Bool_t corrExtrThr5 = plotter->ExtractCorrelations_MC(2.,3);
  printf("*** Extracting correlations in 8<pT(D)<16 GeV/c, 3<pT(assoc)<99 GeV/c ***\n");
  Bool_t corrExtrThr6 = plotter->ExtractCorrelations_MC(3.,99.);
/*  printf("*** Extracting correlations in 8<pT(D)<16 GeV/c, 2<pT(assoc)<99 GeV/c ***\n");
  Bool_t corrExtrThr4 = plotter->ExtractCorrelations_MC(2.,99.);
  printf("*** Extracting correlations in 8<pT(D)<16 GeV/c, 3<pT(assoc)<99 GeV/c ***\n");
  Bool_t corrExtrThr5 = plotter->ExtractCorrelations_MC(3.,99.);    
  printf("*** Extracting correlations in 8<pT(D)<16 GeV/c, 4<pT(assoc)<99 GeV/c ***\n");
  Bool_t corrExtrThr6 = plotter->ExtractCorrelations_MC(4.,99.);
  printf("*** Extracting correlations in 8<pT(D)<16 GeV/c, 1<pT(assoc)<2 GeV/c ***\n");
  Bool_t corrExtrThr7 = plotter->ExtractCorrelations_MC(1.,2.);
  printf("*** Extracting correlations in 8<pT(D)<16 GeV/c, 1<pT(assoc)<3 GeV/c ***\n");
  Bool_t corrExtrThr8 = plotter->ExtractCorrelations_MC(1.,3.);
  printf("*** Extracting correlations in 8<pT(D)<16 GeV/c, 2<pT(assoc)<3 GeV/c ***\n");
  Bool_t corrExtrThr9 = plotter->ExtractCorrelations_MC(2.,3.); */
  if(!corrExtrThr1 || !corrExtrThr2 || !corrExtrThr3 /*|| !corrExtrThr4 || !corrExtrThr5 || !corrExtrThr6 || !corrExtrThr7 || !corrExtrThr8 || !corrExtrThr9*/) {
    printf("Error in the extraction of the correlation distributions! Exiting...\n");
    return;
  }

  plotter->ClearObjects(); //important! Call it after each wide-pT range

}


//________________________________________
void ExtractVeryHighPt_MC(AliDhCorrelationExtraction *plotter) {
	
  plotter->SetNpTbins(2);
  plotter->SetFirstpTbin(11);
  plotter->SetLastpTbin(12);

  plotter->IntegratePtBins(kFALSE); //For D+, set it kTRUE; High pT bin (and only it) has to be evaluated merging mass plots and THnSparses for the 5 pT(D+) bins!

  //extract correlation distributions
  printf("*** Extracting correlations in 16<pT(D)<24 GeV/c, 0.3<pT(assoc)<99 GeV/c ***\n");
  Bool_t corrExtrThr1 = plotter->ExtractCorrelations_MC(0.3,99.);
  printf("*** Extracting correlations in 16<pT(D)<24 GeV/c, 0.3<pT(assoc)<1 GeV/c ***\n");
  Bool_t corrExtrThr2 = plotter->ExtractCorrelations_MC(0.3,1.);
  printf("*** Extracting correlations in 16<pT(D)<24 GeV/c, 1<pT(assoc)<99 GeV/c ***\n");
  Bool_t corrExtrThr3 = plotter->ExtractCorrelations_MC(1.,99.);
  printf("*** Extracting correlations in 16<pT(D)<24 GeV/c, 1<pT(assoc)<2 GeV/c ***\n");
  Bool_t corrExtrThr4 = plotter->ExtractCorrelations_MC(1.,2.);
  printf("*** Extracting correlations in 16<pT(D)<24 GeV/c, 2<pT(assoc)<3 GeV/c ***\n");
  Bool_t corrExtrThr5 = plotter->ExtractCorrelations_MC(2.,3);
  printf("*** Extracting correlations in 16<pT(D)<24 GeV/c, 3<pT(assoc)<99 GeV/c ***\n");
  Bool_t corrExtrThr6 = plotter->ExtractCorrelations_MC(3.,99.);
/*  printf("*** Extracting correlations in 8<pT(D)<16 GeV/c, 2<pT(assoc)<99 GeV/c ***\n");
  Bool_t corrExtrThr4 = plotter->ExtractCorrelations_MC(2.,99.);
  printf("*** Extracting correlations in 8<pT(D)<16 GeV/c, 3<pT(assoc)<99 GeV/c ***\n");
  Bool_t corrExtrThr5 = plotter->ExtractCorrelations_MC(3.,99.);    
  printf("*** Extracting correlations in 8<pT(D)<16 GeV/c, 4<pT(assoc)<99 GeV/c ***\n");
  Bool_t corrExtrThr6 = plotter->ExtractCorrelations_MC(4.,99.);
  printf("*** Extracting correlations in 8<pT(D)<16 GeV/c, 1<pT(assoc)<2 GeV/c ***\n");
  Bool_t corrExtrThr7 = plotter->ExtractCorrelations_MC(1.,2.);
  printf("*** Extracting correlations in 8<pT(D)<16 GeV/c, 1<pT(assoc)<3 GeV/c ***\n");
  Bool_t corrExtrThr8 = plotter->ExtractCorrelations_MC(1.,3.);
  printf("*** Extracting correlations in 8<pT(D)<16 GeV/c, 2<pT(assoc)<3 GeV/c ***\n");
  Bool_t corrExtrThr9 = plotter->ExtractCorrelations_MC(2.,3.); */
  if(!corrExtrThr1 || !corrExtrThr2 || !corrExtrThr3 /*|| !corrExtrThr4 || !corrExtrThr5 || !corrExtrThr6 || !corrExtrThr7 || !corrExtrThr8 || !corrExtrThr9*/) {
    printf("Error in the extraction of the correlation distributions! Exiting...\n");
    return;
  }

  plotter->ClearObjects(); //important! Call it after each wide-pT range

}
//________________________________________
void DrawMCClosure(Int_t specie=AliDhCorrelationExtraction::kD0toKpi, Int_t nOrig=6) {
	
  //Create and set the correlation plotter class
  AliDhCorrelationExtraction *plotter = new AliDhCorrelationExtraction();
  Bool_t flagSpecie = plotter->SetDmesonSpecie(specie);
  if(!flagSpecie) return;

  plotter->AddOriginType("",AliDhCorrelationExtraction::kOrigAll);
  plotter->AddOriginType("_HF_From_c",AliDhCorrelationExtraction::kOrigC);
  plotter->AddOriginType("_HF_From_b",AliDhCorrelationExtraction::kOrigB);
  plotter->AddOriginType("_NonHF",AliDhCorrelationExtraction::kOrigLF);
  plotter->AddOriginType("_From_c",AliDhCorrelationExtraction::kOrigC);
  plotter->AddOriginType("_From_b",AliDhCorrelationExtraction::kOrigB);
  
  //Build plots
  plotter->DrawMCClosure(nOrig,4,5,0.3,99.);
  plotter->DrawMCClosure(nOrig,4,5,0.3,1.);
  plotter->DrawMCClosure(nOrig,4,5,1.,99.);
  plotter->DrawMCClosure(nOrig,4,5,1.,2.);
  plotter->DrawMCClosure(nOrig,4,5,2.,3.);
  plotter->DrawMCClosure(nOrig,4,5,3.,99.);
  plotter->DrawMCClosure(nOrig,6,8,0.3,99.);
  plotter->DrawMCClosure(nOrig,6,8,0.3,1.);
  plotter->DrawMCClosure(nOrig,6,8,1.,99.);
  plotter->DrawMCClosure(nOrig,6,8,1.,2.);
  plotter->DrawMCClosure(nOrig,6,8,2.,3.);
  plotter->DrawMCClosure(nOrig,6,8,3.,99.);
  plotter->DrawMCClosure(nOrig,9,10,0.3,99.);
  plotter->DrawMCClosure(nOrig,9,10,0.3,1.);
  plotter->DrawMCClosure(nOrig,9,10,1.,99.);
  plotter->DrawMCClosure(nOrig,9,10,1.,2.);
  plotter->DrawMCClosure(nOrig,9,10,2.,3.);
  plotter->DrawMCClosure(nOrig,9,10,3.,99.);
  plotter->DrawMCClosure(nOrig,11,12,0.3,99.);
  plotter->DrawMCClosure(nOrig,11,12,0.3,1.);
  plotter->DrawMCClosure(nOrig,11,12,1.,99.);
  plotter->DrawMCClosure(nOrig,11,12,1.,2.);
  plotter->DrawMCClosure(nOrig,11,12,2.,3.);
  plotter->DrawMCClosure(nOrig,11,12,3.,99.);
  
  return;
  
}
