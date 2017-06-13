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
// Last updated: Fabio Colamaria on 18/04/2017

void ExtractOutput_Dplus(
                         Bool_t treeSE=kFALSE, Bool_t treeME=kFALSE, //read TH2F from offline correlation framewrks (produced from the TTree) for SE, ME analysis instead of THnSparse
                         Int_t specie=AliDhCorrelationExtraction::kDplusKpipi, //the D-meson decay channel (check the enumerator for the options)
                         Int_t SandB=AliDhCorrelationExtraction::kBfromBinCount, //how to extract S and B (check the enumerator for the options) - kBfromBinCount is the paper approach
                         Int_t SBscale=AliDhCorrelationExtraction::kBinCountScaling, //how to renormalize the sidebands (check the enumerator for the options) - kBfromBinCount is the paper approach
                         Int_t rebin=1, //rebin the invariant mass plots - USE WITH CARE! (different bin width w.r.t. THnSparse)
                         Double_t leftRng=1.7, Double_t rightRng=2.05, //invariant mass fit range -> use 1.7-2.1 for D0 and D+, 0.14-0.16 for D* (but 1.695-2.1 for D0 in pp for results before 1/5/2015)
                         Int_t funcBkg=AliHFMassFitter::kExpo, //background function used for the mass fit -> use kExpo for D0 and D+, kPowEx for D*
                         Double_t nsigmaFitter=2, //number of sigma in which to extract S, B, S/B, signficance... (only for the spectra visualization, no influence on the correlations)
                         Bool_t autoSign=kTRUE, //kTRUE = define Sign from the fit results, in the range nsigmaS (below); kFALSE = use ranges provided via SetSignRanges
                         Double_t nsigmaS=2, //number of sigma in which to define the signal region (for correlations), Valid only if autoSign flag = kTRUE, otherwise dummy
                         Bool_t autoSB=kTRUE, //kTRUE = define SB from the fit results, in the range insigma-outsigma (below); kFALSE = use ranges provided via SetSBRanges
                         Double_t insigma=4, Double_t outsigma=8, //sideband ranges (in units of sigma). Valid only if autoSB=kTRUES, otherwise dummy
                         Bool_t singleBinSB=kFALSE, //kTRUE=a single mass bin is used in the THnSparse for storing the sideband corlelations (used by D0 to save space)
                         Int_t npools=9, //number of pools for the event-mixing
                         Bool_t poolByPool=kTRUE, //kTRUE=pool-by-pool ME correction; kFALSE=merged-pools ME correction (set the options that you used in the online analysis)
                         Double_t deltaEtaMin=-1., Double_t deltaEtaMax=1.,
			 Bool_t use2Dmassplots=kFALSE, Double_t mincent=0., Double_t maxcent=100.) //deltaEta ranges for correlation distributions
{
    Int_t num=0;
    //Create and set the correlation plotter class
    AliDhCorrelationExtraction *plotter = new AliDhCorrelationExtraction();
    Bool_t flagSpecie = plotter->SetDmesonSpecie(specie);
    plotter->SetSandBextraction(SandB);
    plotter->SetSBscaling(SBscale);
    plotter->SetRebinMassPlots(rebin);
    plotter->SetFitRanges(leftRng,rightRng); //use 1.7-2.1 for D0 and D+, 0.14-0.16 for D*
    plotter->SetBkgFitFunction(funcBkg); //use kExpo for D0 and D+, kPowEx for D*
    plotter->SetNumberOfSigmasFitter(nsigmaFitter);
    if(!autoSign) plotter->SetAutoSignRange(autoSign);
    plotter->SetSignalSigmas(nsigmaS);
    plotter->SetAutoSBRange(autoSB,insigma,outsigma); //kTRUE = evaluate SB range automatically (give inner and outer sigma as 2° and 3° args); kFALSE = use ranges provided via SetSBRanges
    plotter->SetSBSingleBin(singleBinSB);
    plotter->SetNpools(npools);
    plotter->SetCorrectPoolsSeparately(poolByPool); //kTRUE = pool.by-pool extraction and correction; kFALSE = merged ME pools
    plotter->SetDeltaEtaRange(deltaEtaMin,deltaEtaMax);

    if(!flagSpecie) return;
    
    plotter->SetDebugLevel(0); //0 = get main results; 1 = get full list of plots; 2 = get debug printouts
    /*
    //-------------setting for online and offline combinations-------------------------------
    std::cout<<""<<std::endl;
    std::cout<< "\033[1;31m Enter 1 for SE and ME Both Online \033[0m\n"<<std::endl;
    std::cout<< "\033[1;31m Enter 2 for SE and ME Both Offline \033[0m\n"<<std::endl;
    std::cout<< "\033[1;31m Enter 3 for SE Online and ME Offline \033[0m\n"<<std::endl;
    std::cout<< "\033[1;31m Enter 4 for SE Offline and ME Online \033[0m\n"<<std::endl;
    std::cin>>num;
    if (num==1)     {treeSE=kFALSE; treeME=kFALSE;}
    else if (num==2){treeSE=kTRUE; treeME=kTRUE;}
    else if (num==3){treeSE=kFALSE; treeME=kTRUE;}
    else if (num==4){treeSE=kTRUE; treeME=kFALSE;}
    else
        std::cout<<" The Choice entered is wrong !!!"<<std::endl;
    */
    plotter->ReadTTreeOutputFiles(treeSE,treeME);
    plotter->SetSubtractSoftPiInMEdistr(kFALSE);
    plotter->SetUseMassVsCentPlots(use2Dmassplots);
    if(use2Dmassplots) plotter->SetCentralitySelection(mincent,maxcent);
    SetInputNames(plotter, treeSE, treeME);  // check the names in the method!!
    
    Bool_t read = plotter->ReadInputs();
    if(!read) {
        printf("Error in reading the input file! Exiting...\n");
        return;
    }
    
    ExtractLowPt(plotter);
    ExtractMidPt(plotter);
    ExtractHighPt(plotter);
    // ExtractNewPt(plotter);
    
    return;
}

//________________________________________
void SetInputNames(AliDhCorrelationExtraction *plotter, Bool_t treeSE, Bool_t treeME){
    
    
    //Dplus paths - paper plots (put merged pools as argument!)
    gSystem->Exec("rm -rf Output_Root/ Output_png/");
    gSystem->Exec("mkdir Output_Root/ Output_png/");
    plotter->SetInputFilenameMass("./AnalysisResults.root");
    //-------------Online SE and ME setting------------
    if (!treeSE && !treeME){
        printf("\033[1;44m Extracting Online SE and ME Correlations!!!!!...\033[0m\n");
        plotter->SetInputFilenameSE("./AnalysisResults.root");
        plotter->SetInputFilenameME("./AnalysisResults.root");
    }
    //-------------Offline SE and ME setting---------------
    else if(treeSE && treeME){
        printf("\033[1;44m Extracting Offline SE and ME Correlations!!!!!...\033[0m\n");
        plotter->SetInputFilenameSE("./OfflineCorrelationsSE.root");
        plotter->SetInputFilenameME("./OfflineCorrelationsME.root");
    }
    //-------------Online SE and Offline ME setting---------------
    else if(!treeSE && treeME){
        printf("\033[1;44m Extracting Online SE and Offline ME Correlations!!!!!...\033[0m\n");
        plotter->SetInputFilenameSE("./AnalysisResults.root");
        plotter->SetInputFilenameME("./OfflineCorrelationsME.root");
    }
    //-------------Offline SE and Onlline ME setting---------------
    else if(treeSE && !treeME){
        printf("\033[1;44m Extracting Offline SE and Online ME Correlations!!!!!...\033[0m\n");
        plotter->SetInputFilenameSE("./OfflineCorrelationsSE.root");
        plotter->SetInputFilenameME("./AnalysisResults.root");
    }
    
    plotter->SetDirNameMass("pPbHadCorr_SE_woTrkEff_woDkEff_SECut2.0FileCutsPlbyPl");
    plotter->SetDirNameSE("pPbHadCorr_SE_woTrkEff_woDkEff_SECut2.0FileCutsPlbyPl");
    plotter->SetListNameMass("coutHistos_pPbHadCorr_SE_woTrkEff_woDkEff_SECut2.0FileCutsPlbyPl"); //old prefix: _PPPass4SE
    plotter->SetListNameSE("coutHistos_pPbHadCorr_SE_woTrkEff_woDkEff_SECut2.0FileCutsPlbyPl");
    plotter->SetDirNameME("pPbHadCorr_ME_woTrkEff_woDkEff_MECut2.0FileCutsPlbyPl");
    plotter->SetListNameME("coutHistos_pPbHadCorr_ME_woTrkEff_woDkEff_MECut2.0FileCutsPlbyPl");
    plotter->SetMassHistoName("hnSparseMa_Hdron_Data_Bin");
    plotter->SetSECorrelHistoName("hnSparseCa_Hdron_Data_Bin");
    plotter->SetMECorrelHistoNameSuffix("_evMix");
    return;
}

//________________________________________
void ExtractLowPt(AliDhCorrelationExtraction *plotter){
    
    plotter->SetNpTbins(2);
    plotter->SetFirstpTbin(3);
    plotter->SetLastpTbin(4);
    
    //Limits of sidebands passed from outside. Warning: use limits of the mass bin edges to avoid further adjustments in the normalization!
    Double_t LSBLowLim[2] = {1.7983,1.7748}; //IMPORTANT!! to be filled accordingly to those set in the task
    Double_t LSBUppLim[2] = {1.8375,1.8297};
    Double_t RSBLowLim[2] = {1.9003,1.9081};
    Double_t RSBUppLim[2] = {1.9395,1.9630};
    // plotter->SetAutoSBRange(kTRUE,4,8);
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
    printf("*** Extracting correlations in 3<pT(D)<5 GeV/c, 0.3<pT(assoc)<99 GeV/c ***\n");
    Bool_t corrExtrThr1 = plotter->ExtractCorrelations(0.3,99.);
    printf("*** Extracting correlations in 3<pT(D)<5 GeV/c, 0.3<pT(assoc)<1 GeV/c ***\n");
    Bool_t corrExtrThr2 = plotter->ExtractCorrelations(0.3,1.);
    printf("*** Extracting correlations in 3<pT(D)<5 GeV/c, 1.0<pT(assoc)<99 GeV/c ***\n");
    Bool_t corrExtrThr3 = plotter->ExtractCorrelations(1.0,99.);
    printf("*** Extracting correlations in 3<pT(D)<5 GeV/c, 2.0<pT(assoc)<99 GeV/c ***\n");
    Bool_t corrExtrThr4 = plotter->ExtractCorrelations(2.0,99.);
    if(!corrExtrThr1 || !corrExtrThr2 || !corrExtrThr3 || !corrExtrThr4) {
        printf("Error in the extraction of the correlation distributions! Exiting...\n");
        return;
    }
    
    plotter->ClearObjects(); //important! Call it after each wide-pT range
    
}


//________________________________________
void ExtractMidPt(AliDhCorrelationExtraction *plotter){
    
    plotter->SetNpTbins(3);
    plotter->SetFirstpTbin(5);
    plotter->SetLastpTbin(7);
    
    //Limits of sidebands passed from outside. Warning: use limits of the mass bin edges to avoid further adjustments in the normalization!
    Double_t LSBLowLim[3] = {1.7905,1.7591,1.7670}; //IMPORTANT!! to be filled accordingly to those set in the task
    Double_t LSBUppLim[3] = {1.8375,1.8219,1.8219};
    Double_t RSBLowLim[3] = {1.9003,1.9160,1.9160};
    Double_t RSBUppLim[3] = {1.9474,1.9787,1.9709};
    //plotter->SetAutoSBRange(kTRUE,4,8);
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
    printf("*** Extracting correlations in 5<pT(D)<8 GeV/c, 0.3<pT(assoc)<99 GeV/c ***\n");
    Bool_t corrExtrThr1 = plotter->ExtractCorrelations(0.3,99.);
    printf("*** Extracting correlations in 5<pT(D)<8 GeV/c, 0.3<pT(assoc)<1 GeV/c ***\n");
    Bool_t corrExtrThr2 = plotter->ExtractCorrelations(0.3,1.);
    printf("*** Extracting correlations in 5<pT(D)<8 GeV/c, 1.0<pT(assoc)<99 GeV/c ***\n");
    Bool_t corrExtrThr3 = plotter->ExtractCorrelations(1.0,99.);
    printf("*** Extracting correlations in 5<pT(D)<8 GeV/c, 2.0<pT(assoc)<99 GeV/c ***\n");
    Bool_t corrExtrThr4 = plotter->ExtractCorrelations(2.,99.);
    if(!corrExtrThr1 || !corrExtrThr2 || !corrExtrThr3 || !corrExtrThr4) {
        printf("Error in the extraction of the correlation distributions! Exiting...\n");
        return;
    }
    
    plotter->ClearObjects(); //important! Call it after each wide-pT range
    
}


//________________________________________
void ExtractHighPt(AliDhCorrelationExtraction *plotter){ // This is for PPb for PP we have first bin=8, LastBin=12 and NPtBins=5 (Shyam)
    
    plotter->SetNpTbins(6);
    plotter->SetFirstpTbin(8);
    plotter->SetLastpTbin(13);
    
    //Limits of sidebands passed from outside. Warning: use limits of the mass bin edges to avoid bias in the normalization!
    Double_t LSBLowLim[6] = {1.7670,1.7670,1.7670,1.7670,1.7670,1.7670}; //IMPORTANT!! to be filled accordingly to those set in the task
    Double_t LSBUppLim[6] = {1.8219,1.8219,1.8219,1.8219,1.8219,1.8219};
    Double_t RSBLowLim[6] = {1.9160,1.9160,1.9160,1.9160,1.9160,1.9160};
    Double_t RSBUppLim[6] = {1.9709,1.9709,1.9709,1.9709,1.9709,1.9709};
    // plotter->SetAutoSBRange(kTRUE,4,8);
    plotter->SetSBRanges(LSBLowLim,LSBUppLim,RSBLowLim,RSBUppLim);
    plotter->IntegratePtBins(kFALSE); //For D+! High pT bin (and only it) has to be evaluated merging mass plots and THnSparses for the 5 pT(D+) bins!
    // Use kTRUE in PP where Stats is less and kFALSE in pPb and PPb where stats is large (Shyam)
    //fit invariant mass distribution and extract S and B for normalizations
    Bool_t massfit = plotter->FitInvariantMass();
    if(!massfit) {
        printf("Error in the fitting of the mass plots! Exiting...\n");
        return;
    }
    
    plotter->PrintRanges();
    plotter->PrintSandBForNormal();
    
    //extract correlation distributions
    printf("*** Extracting correlations in 8<pT(D)<16 GeV/c, 0.3<pT(assoc)<99 GeV/c ***\n");
    Bool_t corrExtrThr1 = plotter->ExtractCorrelations(0.3,99.);
    printf("*** Extracting correlations in 8<pT(D)<16 GeV/c, 0.3<pT(assoc)<1 GeV/c ***\n");
    Bool_t corrExtrThr2 = plotter->ExtractCorrelations(0.3,1.);
    printf("*** Extracting correlations in 8<pT(D)<16 GeV/c, 1.0<pT(assoc)<99 GeV/c ***\n");
    Bool_t corrExtrThr3 = plotter->ExtractCorrelations(1.0,99.);
    printf("*** Extracting correlations in 8<pT(D)<16 GeV/c, 2.0<pT(assoc)<99 GeV/c ***\n");
    Bool_t corrExtrThr4 = plotter->ExtractCorrelations(2.,99.);
    if(!corrExtrThr1 || !corrExtrThr2 || !corrExtrThr3 || !corrExtrThr4) {
        printf("Error in the extraction of the correlation distributions! Exiting...\n");
        return;
    }
    
    plotter->ClearObjects(); //important! Call it after each wide-pT range
    
}
void ExtractNewPt(AliDhCorrelationExtraction *plotter){
    
    plotter->SetNpTbins(1);
    plotter->SetFirstpTbin(14);
    plotter->SetLastpTbin(14);
    
    //Limits of sidebands passed from outside. Warning: use limits of the mass bin edges to avoid further adjustments in the normalization!
    Double_t LSBLowLim[1] = {1.7983}; //IMPORTANT!! to be filled accordingly to those set in the task
    Double_t LSBUppLim[1] = {1.8375};
    Double_t RSBLowLim[1] = {1.9003};
    Double_t RSBUppLim[1] = {1.9395};
    // plotter->SetAutoSBRange(kTRUE,4,8);
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
    printf("*** Extracting correlations in 16<pT(D)<24 GeV/c, 0.3<pT(assoc)<99 GeV/c ***\n");
    Bool_t corrExtrThr1 = plotter->ExtractCorrelations(0.3,99.);
    printf("*** Extracting correlations in 16<pT(D)<24 GeV/c, 0.3<pT(assoc)<1 GeV/c ***\n");
    Bool_t corrExtrThr2 = plotter->ExtractCorrelations(0.3,1.);
    printf("*** Extracting correlations in 16<pT(D)<24 GeV/c, 1.0<pT(assoc)<99 GeV/c ***\n");
    Bool_t corrExtrThr3 = plotter->ExtractCorrelations(1.0,99.);
    printf("*** Extracting correlations in 16<pT(D)<24 GeV/c, 2.0<pT(assoc)<99 GeV/c ***\n");
    Bool_t corrExtrThr4 = plotter->ExtractCorrelations(2.0,99.);
    if(!corrExtrThr1 || !corrExtrThr2 || !corrExtrThr3 || !corrExtrThr4) {
        printf("Error in the extraction of the correlation distributions! Exiting...\n");
        return;
    }
    
    plotter->ClearObjects(); //important! Call it after each wide-pT range
    
}
