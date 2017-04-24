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
// - The SetSignRanges values in the ExtractXPt methods [to be created from scratch] (do only if you use autoSign=kFALSE)
//
// NOTE FOR D*: If setting the sideband range externally, the (unique) sideband to be used is LSB (not RSB)
// NOTE FOR D+: If you need to integrate mass-pT bins in a single correlation bin BEFORE extracting the outputs, set plotter->IntegratePtBins(kFALSE) in the ExtractXPt methods 
// 
// setting plotter->SetDebugLevel(1) will print a series of additional plots for debugging; plotter->SetDebugLevel(2) will also be more verbose
//

//Only the most important plots

  Int_t fptBinD0[4][2]={
  {4,5},   // ExtractLowPt(plotter);
  {6,8},   // ExtractMidPt(plotter);
  {9,10},  // ExtractHighPt(plotter);
  {4,10}   // ExtractFullPt(plotter);
  };


//Full chain
/*
Int_t fptBinD0[6][2]={
  {3,5},   // ExtractExtraLowPt(plotter);
  {4,5},   // ExtractLowPt(plotter);
  {6,8},   // ExtractMidPt(plotter);
  {9,10},  // ExtractHighPt(plotter);
  {3,10},  // ExtractFullPtExt(plotter);
  {4,10}   // ExtractFullPt(plotter);
};
*/

const int fnCutSetsD0 = sizeof(fptBinD0)/8;

//Only the most important plots
  Double_t fptEl[3][2]={
  {1.0,4.0},
  {1.0,10.0},
  {1.0,5.0}
  };

/*
//Full chain
Double_t fptEl[23][2]={
  {1.5,2.0},
  {1.5,3.0},
  {1.5,4.0},
  {1.5,5.0},
  {1.5,6.0},
  {1.5,10.0},

  {0.3,1.0},
  {0.3,2.0},
  {0.3,3.0},
  {0.3,4.0},
  {0.3,5.0},
  {0.3,10.0},

  {0.5,1.0},
  {0.5,2.0},
  {0.5,3.0},
  {0.5,4.0},
  {0.5,5.0},
  {0.5,10.0},

  {1.0,2.0},
  {1.0,3.0},
  {1.0,4.0},
  {1.0,5.0},
  {1.0,10.0}
};
*/
const int fnCutSetsEl = sizeof(fptEl)/16;
void ExtractOutput(
   Bool_t treeSE=kFALSE, Bool_t treeME=kFALSE, //read TH2F from offline correlation framewrks (produced from the TTree) for SE, ME analysis instead of THnSparse
   Int_t specie=AliDhCorrelationExtraction::kDxHFE, //kD0toKpi, //the D-meson decay channel (check the enumerator for the options)
   Int_t SandB=AliDhCorrelationExtraction::kBfromBinCount, //how to extract S and B (check the enumerator for the options) - kBfromBinCount is the paper approach
   Int_t SBscale=AliDhCorrelationExtraction::kBinCountScaling, //how to renormalize the sidebands (check the enumerator for the options) - kBfromBinCount is the paper approach
   Int_t rebin=1, //rebin the invariant mass plots - USE WITH CARE! (different bin width w.r.t. THnSparse)
   Double_t leftRng=1.7, Double_t rightRng=2.1, //invariant mass fit range -> use 1.7-2.1 for D0 and D+, 0.14-0.16 for D* (but 1.695-2.1 for D0 in pp for results before 1/5/2015)
   Int_t funcBkg=AliHFMassFitter::kExpo, //background function used for the mass fit -> use kExpo for D0 and D+, kPowEx for D*
   Double_t nsigmaFitter=3, //number of sigma in which to extract S, B, S/B, signficance... (only for the spectra visualization, no influence on the correlations)
   Bool_t autoSign=kTRUE, //kTRUE = define Sign from the fit results, in the range nsigmaS (below); kFALSE = use ranges provided via SetSignRanges
   Double_t nsigmaS=2, //number of sigma in which to define the signal region (for correlations), Valid only if autoSign flag = kTRUE, otherwise dummy
   Bool_t autoSB=kFALSE, //kTRUE = define SB from the fit results, in the range insigma-outsigma (below); kFALSE = use ranges provided via SetSBRanges
   Double_t insigma=0, Double_t outsigma=0, //sideband ranges (in units of sigma). Valid only if autoSB=kTRUES, otherwise dummy
   Bool_t singleBinSB=kFALSE, //kTRUE=a single mass bin is used in the THnSparse for storing the sideband corlelations (used by D0 to save space)
   Int_t npools=6, //number of pools for the event-mixing
   Bool_t poolByPool=kTRUE, //kTRUE=pool-by-pool ME correction; kFALSE=merged-pools ME correction (set the options that you used in the online analysis)
   Double_t deltaEtaMin=-1., Double_t deltaEtaMax=1., //deltaEta ranges for correlation distributions
   Bool_t useMC=kFALSE, Int_t elSource=AliDhCorrelationExtraction::kAll, Int_t D0Source=AliDhCorrelationExtraction::kAll) //Adjust parameters for runs over MC

{
  gSystem->Exec("mkdir -p Output_png"); 
  gSystem->Exec("mkdir -p Output_Root"); 
  //Create and set the correlation plotter class
  AliDhCorrelationExtraction *plotter = new AliDhCorrelationExtraction();
  Bool_t flagSpecie = plotter->SetDmesonSpecie(specie);
  plotter->SetSandBextraction(SandB);
  plotter->SetSBscaling(SBscale);
  plotter->SetRebinMassPlots(rebin);
  plotter->SetFitRanges(leftRng,rightRng); //
  plotter->SetBkgFitFunction(funcBkg); //
  plotter->SetNumberOfSigmasFitter(nsigmaFitter);
  if(autoSign) plotter->SetAutoSignRange(autoSign);
  plotter->SetSignalSigmas(nsigmaS);
  plotter->SetAutoSBRange(autoSB,insigma,outsigma);
  plotter->SetSBSingleBin(singleBinSB);
  plotter->SetNpools(npools);
  plotter->SetCorrectPoolsSeparately(poolByPool);
  plotter->SetDeltaEtaRange(deltaEtaMin,deltaEtaMax);
  plotter->ReadTTreeOutputFiles(treeSE,treeME);
  plotter->SetUseMC(useMC);
  plotter->SetUseElSource(elSource);
  plotter->SetUseD0Source(D0Source);
  if(!flagSpecie) return;
  plotter->SetDebugLevel(1); //0 = get main results; 1 = get full list of plots; 2 = get debug printouts
  SetInputNames(plotter, useMC);  // check the names in the method!!
  
  Bool_t read = plotter->ReadInputs();
  if(!read) {
    printf("Error in reading the input file! Exiting...\n");
    return;
  }
  //  ExtractExtraLowPt(plotter);
  ExtractLowPt(plotter);
  ExtractMidPt(plotter);
  ExtractHighPt(plotter);
  ExtractFullPt(plotter);
  //  ExtractFullPtExt(plotter);
  reflectHistoDxHFE();
  return;
}

//________________________________________
void SetInputNames(AliDhCorrelationExtraction *plotter, Bool_t useMC){

  //D2H mass or from our task?  plotter->SetInputFilenameMass("./AnalysisResults.root");
  plotter->SetInputFilenameMass("./DxHFECorrelation.root");
  plotter->SetInputFilenameSE("./DxHFECorrelation.root");
  plotter->SetInputFilenameME("./DxHFECorrelation.root");
  plotter->SetDirNameMass("DxHFE");
  plotter->SetListNameMass("DxHFElist4ITSkF");
  plotter->SetDirNameSE("DxHFE");
  plotter->SetListNameSE("DxHFElist4ITSkF");
  plotter->SetDirNameME("DxHFEME");
  plotter->SetListNameME("DxHFElist4ITSkFME");

  plotter->SetMassHistoName("D0 info");
  if(!useMC){plotter->SetSECorrelHistoName("AliDxHFECorrelation info");}
  else {plotter->SetSECorrelHistoName("AliDxHFECorrelationMC info");}
  plotter->SetMECorrelHistoNameSuffix("");
  return;
}

//________________________________________
void ExtractExtraLowPt(AliDhCorrelationExtraction *plotter){

  plotter->SetNpTbins(3);
  plotter->SetFirstpTbin(3);
  plotter->SetLastpTbin(5);

  //Limits of sidebands passed from outside. Warning: use limits of the mass bin edges to avoid further adjustments in the normalization!
  Double_t LSBLowLim[3] = {1.7448,1.7448,1.7368}; //IMPORTANT!! to be filled accordingly to those set in the task
  Double_t LSBUppLim[3] = {1.8088,1.8088,1.8008};
  Double_t RSBLowLim[3] = {1.9288,1.9288,1.9288};
  Double_t RSBUppLim[3] = {1.9848,1.9848,1.9928};

  plotter->SetSBRanges(LSBLowLim,LSBUppLim,RSBLowLim,RSBUppLim);
  plotter->IntegratePtBins(kFALSE);

  //fit invariant mass distribution and extract S and B for normalizations
  Bool_t massfit = plotter->FitInvariantMass();
  if(!massfit) {
    printf("Error in the fitting of the mass plots! Exiting...\n");
    return;
  }

  plotter->PrintRanges();
  plotter->PrintSandBForNormal();

  //extract correlation distributions
  for(int i=0; i<fnCutSetsEl; i++){
    printf(Form("*** Extracting correlations in 2<pT(D)<5 GeV/c, %g<pT(assoc)<%g GeV/c ***\n", fptEl[i][0], fptEl[i][1]));
    Bool_t corrExtrThr = plotter->ExtractCorrelations(fptEl[i][0],fptEl[i][1]);
    if(!corrExtrThr ) {
      printf("Error in the extraction of the correlation distributions! Exiting...\n");
      return;
    }  
  }
  plotter->ClearObjects(); //important! Call it after each wide-pT range

}


//________________________________________
void ExtractLowPt(AliDhCorrelationExtraction *plotter){

  plotter->SetNpTbins(2);
  plotter->SetFirstpTbin(4);
  plotter->SetLastpTbin(5);

  //Limits of sidebands passed from outside. Warning: use limits of the mass bin edges to avoid further adjustments in the normalization!
  Double_t LSBLowLim[2] = {1.7448,1.7368}; //IMPORTANT!! to be filled accordingly to those set in the task
  Double_t LSBUppLim[2] = {1.8088,1.8008};
  Double_t RSBLowLim[2] = {1.9288,1.9288};
  Double_t RSBUppLim[2] = {1.9848,1.9928};

  plotter->SetSBRanges(LSBLowLim,LSBUppLim,RSBLowLim,RSBUppLim);
  plotter->IntegratePtBins(kFALSE);

  //fit invariant mass distribution and extract S and B for normalizations
  Bool_t massfit = plotter->FitInvariantMass();
  if(!massfit) {
    printf("Error in the fitting of the mass plots! Exiting...\n");
    return;
  }

  plotter->PrintRanges();
  plotter->PrintSandBForNormal();

  //extract correlation distributions
  for(int i=0; i<fnCutSetsEl; i++){
    printf(Form("*** Extracting correlations in 2<pT(D)<5 GeV/c, %g<pT(assoc)<%g GeV/c ***\n", fptEl[i][0], fptEl[i][1]));
    Bool_t corrExtrThr = plotter->ExtractCorrelations(fptEl[i][0],fptEl[i][1]);
    if(!corrExtrThr ) {
      printf("Error in the extraction of the correlation distributions! Exiting...\n");
      return;
    }  
  }
  plotter->ClearObjects(); //important! Call it after each wide-pT range

}


//________________________________________
void ExtractMidPt(AliDhCorrelationExtraction *plotter){

  plotter->SetNpTbins(3);
  plotter->SetFirstpTbin(6);
  plotter->SetLastpTbin(8);

  //Limits of sidebands passed from outside. Warning: use limits of the mass bin edges to avoid further adjustments in the normalization!
  Double_t LSBLowLim[3] = {1.7048,1.7128,1.7128}; //IMPORTANT!! to be filled accordingly to those set in the task
  Double_t LSBUppLim[3] = {1.7848,1.7928,1.7928};
  Double_t RSBLowLim[3] = {1.9448,1.9448,1.9526};
  Double_t RSBUppLim[3] = {2.0248,2.0248,2.0248};
  plotter->SetSBRanges(LSBLowLim,LSBUppLim,RSBLowLim,RSBUppLim);
  plotter->IntegratePtBins(kFALSE);

  //fit invariant mass distribution and extract S and B for normalizations
  Bool_t massfit = plotter->FitInvariantMass();
  if(!massfit) {
    printf("Error in the fitting of the mass plots! Exiting...\n");
    return;
  }

  plotter->PrintRanges();
  plotter->PrintSandBForNormal();

  //extract correlation distributions
  for(int i=0; i<fnCutSetsEl; i++){
    printf(Form("*** Extracting correlations in 2<pT(D)<5 GeV/c, %g<pT(assoc)<%g GeV/c ***\n", fptEl[i][0], fptEl[i][1]));
    Bool_t corrExtrThr = plotter->ExtractCorrelations(fptEl[i][0],fptEl[i][1]);
    if(!corrExtrThr ) {
      printf("Error in the extraction of the correlation distributions! Exiting...\n");
      return;
    }  
  }
  plotter->ClearObjects(); //important! Call it after each wide-pT range

}


//________________________________________
void ExtractHighPt(AliDhCorrelationExtraction *plotter){

  plotter->SetNpTbins(2);
  plotter->SetFirstpTbin(9);
  plotter->SetLastpTbin(10);

  //Limits of sidebands passed from outside. Warning: use limits of the mass bin edges to avoid bias in the normalization!
  Double_t LSBLowLim[2] = {1.7048,1.7128}; //IMPORTANT!! to be filled accordingly to those set in the task
  Double_t LSBUppLim[2] = {1.7528,1.7688};
  Double_t RSBLowLim[2] = {1.9768,1.9768};
  Double_t RSBUppLim[2] = {2.0888,2.0808};

  plotter->SetSBRanges(LSBLowLim,LSBUppLim,RSBLowLim,RSBUppLim);
  plotter->IntegratePtBins(kFALSE); //For D+, set it kTRUE; High pT bin (and only it) has to be evaluated merging mass plots and THnSparses for the 5 pT(D+) bins!

  //fit invariant mass distribution and extract S and B for normalizations
  Bool_t massfit = plotter->FitInvariantMass();
  if(!massfit) {
    printf("Error in the fitting of the mass plots! Exiting...\n");
    return;
  }

  plotter->PrintRanges();
  plotter->PrintSandBForNormal();

  //extract correlation distributions
  for(int i=0; i<fnCutSetsEl; i++){
    printf(Form("*** Extracting correlations in 2<pT(D)<5 GeV/c, %g<pT(assoc)<%g GeV/c ***\n", fptEl[i][0], fptEl[i][1]));
    Bool_t corrExtrThr = plotter->ExtractCorrelations(fptEl[i][0],fptEl[i][1]);
    if(!corrExtrThr ) {
      printf("Error in the extraction of the correlation distributions! Exiting...\n");
      return;
    }  
  }
  plotter->ClearObjects(); //important! Call it after each wide-pT range

}

//________________________________________
void ExtractFullPt(AliDhCorrelationExtraction *plotter){

  plotter->SetNpTbins(7);
  plotter->SetFirstpTbin(4);
  plotter->SetLastpTbin(10);
  //Limits of sidebands passed from outside. Warning: use limits of the mass bin edges to avoid bias in the normalization!
  Double_t LSBLowLim[7] = {1.7448,1.7368,1.7048,1.7128,1.7128,1.7048,1.7128}; //IMPORTANT!! to be filled accordingly to those set in the task
  Double_t LSBUppLim[7] = {1.8088,1.8008,1.7848,1.7928,1.7928,1.7528,1.7688};
  Double_t RSBLowLim[7] = {1.9288,1.9288,1.9448,1.9448,1.9526,1.9768,1.9768};
  Double_t RSBUppLim[7] = {1.9848,1.9928,2.0248,2.0248,2.0248,2.0888,2.0808};

  plotter->SetSBRanges(LSBLowLim,LSBUppLim,RSBLowLim,RSBUppLim);
  plotter->IntegratePtBins(kFALSE); //For D+, set it kTRUE; High pT bin (and only it) has to be evaluated merging mass plots and THnSparses for the 5 pT(D+) bins!

  //fit invariant mass distribution and extract S and B for normalizations
  Bool_t massfit = plotter->FitInvariantMass();
  if(!massfit) {
    printf("Error in the fitting of the mass plots! Exiting...\n");
    return;
  }

  plotter->PrintRanges();
  plotter->PrintSandBForNormal();
  //extract correlation distributions
  for(int i=0; i<fnCutSetsEl; i++){
    printf(Form("*** Extracting correlations in 2<pT(D)<5 GeV/c, %g<pT(assoc)<%g GeV/c ***\n", fptEl[i][0], fptEl[i][1]));
    Bool_t corrExtrThr = plotter->ExtractCorrelations(fptEl[i][0],fptEl[i][1]);
    if(!corrExtrThr ) {
      printf("Error in the extraction of the correlation distributions! Exiting...\n");
      return;
    }  
  }
  plotter->ClearObjects(); //important! Call it after each wide-pT range

}

void ExtractFullPtExt(AliDhCorrelationExtraction *plotter){

  plotter->SetNpTbins(8);
  plotter->SetFirstpTbin(3);
  plotter->SetLastpTbin(10);
  //Limits of sidebands passed from outside. Warning: use limits of the mass bin edges to avoid bias in the normalization!
  Double_t LSBLowLim[8] = {1.7448,1.7448,1.7368,1.7048,1.7128,1.7128,1.7048,1.7128}; //IMPORTANT!! to be filled accordingly to those set in the task
  Double_t LSBUppLim[8] = {1.8088,1.8088,1.8008,1.7848,1.7928,1.7928,1.7528,1.7688};
  Double_t RSBLowLim[8] = {1.9288,1.9288,1.9288,1.9448,1.9448,1.9526,1.9768,1.9768};
  Double_t RSBUppLim[8] = {1.9848,1.9848,1.9928,2.0248,2.0248,2.0248,2.0888,2.0808};

  plotter->SetSBRanges(LSBLowLim,LSBUppLim,RSBLowLim,RSBUppLim);
  plotter->IntegratePtBins(kFALSE); //For D+, set it kTRUE; High pT bin (and only it) has to be evaluated merging mass plots and THnSparses for the 5 pT(D+) bins!

  //fit invariant mass distribution and extract S and B for normalizations
  Bool_t massfit = plotter->FitInvariantMass();
  if(!massfit) {
    printf("Error in the fitting of the mass plots! Exiting...\n");
    return;
  }

  plotter->PrintRanges();
  plotter->PrintSandBForNormal();
  //extract correlation distributions
  for(int i=0; i<fnCutSetsEl; i++){
    printf(Form("*** Extracting correlations in 2<pT(D)<5 GeV/c, %g<pT(assoc)<%g GeV/c ***\n", fptEl[i][0], fptEl[i][1]));
    Bool_t corrExtrThr = plotter->ExtractCorrelations(fptEl[i][0],fptEl[i][1]);
    if(!corrExtrThr ) {
      printf("Error in the extraction of the correlation distributions! Exiting...\n");
      return;
    }  
  }
  plotter->ClearObjects(); //important! Call it after each wide-pT range

}


void reflectHistoDxHFE(){
 gInterpreter->ExecuteMacro("$ALICE_PHYSICS/../src/PWGHF/correlationHF/macros/LoadLibraries.C");
 gSystem->Exec("mkdir -p Output_Root/Reflect/");
 gSystem->Exec("mkdir -p Output_png/Reflect/");
 for (int i=0;i<fnCutSetsD0;i++){
   for (int j=0;j<fnCutSetsEl;j++){
     createReflectHisto(fptBinD0[i][0],fptBinD0[i][1],fptEl[j][0],fptEl[j][1]);
   }
 }
}
void createReflectHisto(Int_t ptBinMin, Int_t ptBinMax, Double_t thrMin, Double_t thrMax){
  TCanvas *c=new TCanvas("cDraw","cDraw",1400,800);
  c->Divide(2,1);
  f=TFile::Open(Form("Output_Root/AzimCorrDistr_Dzero_Canvas_PtIntBins%dto%d_PoolInt_thr%.1fto%.1f.root",ptBinMin,ptBinMax,thrMin,thrMax));
  TCanvas *ch = (TCanvas*)f->Get(Form("cFinal_%.1fto%.1f",thrMin,thrMax));
  TH1D * h1 = (TH1D*)ch->GetPrimitive("h1D_SubtrNorm");
  //  h1->Sumw2();
  c->cd(1);
  //  h1->Rebin(2);
  h1->Sumw2();
  h1->SetMaximum(h1->GetBinContent(h1->GetMaximumBin())*1.3);
  h1->Draw();
  TH1D *hRef;
  hRef=AliHFCorrelationUtils::ReflectHisto(h1,0.5);
  c->cd(2);
  hRef->Draw();
  c->SaveAs(Form("Output_Root/Reflect/AzimCorrDistr_Dzero_Canvas_PtIntBins%dto%d_PoolInt_thr%.1fto%.1f.root",ptBinMin,ptBinMax,thrMin,thrMax));
  c->SaveAs(Form("Output_png/Reflect/AzimCorrDistr_Dzero_Canvas_PtIntBins%dto%d_PoolInt_thr%.1fto%.1f.png",ptBinMin,ptBinMax,thrMin,thrMax));

}
