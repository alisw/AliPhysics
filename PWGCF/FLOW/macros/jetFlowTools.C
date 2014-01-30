void jetFlowTools() {
    // load and compile the libraries
    // make sure that you have ROOUNFOLD available on your machine,
    // (see http://hepunx.rl.ac.uk/~adye/software/unfold/RooUnfold.html ).
    // and make sure that the Load() function knows where to find
    // the libraries
    Load();

    // read detector response from output of matching task
    // AliAnalysisTaskJetMatching
    // the detector response can also be set manually by
    // calling AliJetFlowTools::SetDetectorResponse(TH2D*)
    TString drInputName = "response.root";
    printf("- Reading file %s ... \n", drInputName.Data());
    TFile drInput(drInputName.Data());          // detector response input matrix
    if(drInput.IsZombie()) {
        printf(" > read error ! < \n");
        return;
    }
    TH2D* detres = (TH2D*)drInput.Get("detector_response");
    if(!detres) {
        printf(" > failed to extract detector respose < \n");
        return;
    } else printf(" > Found detector response < \n");

    // get a TList from the AliAnalysisRhoVnModulation task
    // this will be used as input for the unfolding (jet spectra, dpt distribution)
    // input can also be set manually, by calling
    // AliJetFlowTools::SetRawInput() see the header of AliJetFlowTools.h for a full
    // list of necessary input histograms
    TFile f("AnalysisResults.root");
    if(f.IsZombie()) {
        printf(" > read error ! < \n");
        return;
    }
    TList* l = (TList*)f.Get("RhoVnMod_R03_kCombined_Jet_AKTChargedR030_PicoTracks_pT0150_Rho_TPC_PWGJE");
    if(!l) {
        printf(" > failed to find output list ! \n");
        return;
    }
    // create an instance of the Tools class
    // one instance will do all the unfolding
    AliJetFlowTools* tools = new AliJetFlowTools();
    // set some common variables 
    tools->SetCentralityBin(2); // bin only makes sense when output is taken from AliAnalysisRhoVnModulation
    tools->SetDetectorResponse(detres);

    // set the true (unfolded) bins
    Double_t binsTrue[] = {5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 150};
    tools->SetBinsTrue(new TArrayD(sizeof(binsTrue)/sizeof(binsTrue[0]), binsTrue));
    // set the same binning scheme to be used when chi2 is taken as a prior
    // if not set, binsTrue is used
    tools->SetBinsTruePrior(new TArrayD(sizeof(binsTrue)/sizeof(binsTrue[0]), binsTrue));


    // set the measured (folded) bins
    Double_t binsRec[36];
    for(Int_t i(0); i < 36; i++) binsRec[i] = 2*i+26;
    tools->SetBinsRec(new TArrayD(sizeof(binsRec)/sizeof(binsRec[0]), binsRec));
    // set the same binning scheme to be used when chi2 is taken as a prior
    tools->SetBinsRecPrior(new TArrayD(sizeof(binsRec)/sizeof(binsRec[0]), binsRec));

    // connect input (when using output from AliAnalysisRhoVnModulation)
    tools->SetInputList(l);
    // unfold using different parameters

    // configuration. for all avaialble options, see AliJetFlowTools.h
    tools->SetSmoothenSpectrum(kTRUE, 50, 100, 70);
    tools->SetNormalizeSpectra(10000);
    tools->SetUseDetectorResponse(kTRUE);
    tools->SetSaveFull(kTRUE);
    // set no unfolding
    tools->SetUnfoldingAlgorithm(AliJetFlowTools::kNone);
    tools->CreateOutputList(TString("do_nothing"));
    tools->Make();
    tools->SetTestMode(kTRUE);
 
    // do some chi2 unfolding
    tools->SetUnfoldingAlgorithm(AliJetFlowTools::kChi2);
    Double_t b = 0.05;
    tools->CreateOutputList(TString(Form("beta%.2f", b)));
    tools->SetBeta(b);
    tools->Make();
    b = 0.1;
    tools->CreateOutputList(TString(Form("beta%.2f", b)));
    tools->SetBeta(b);
    tools->Make();

    // do some SVD unfolding
    tools->SetUnfoldingAlgorithm(AliJetFlowTools::kSVD);
    tools->SetPrior(AliJetFlowTools::kPriorChi2);

    // svd unfolding prefers diffefrent binning
    Double_t binsRec2[] = {25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 95};
    tools->SetBinsRec(new TArrayD(sizeof(binsRec2)/sizeof(binsRec2[0]), binsRec2));
    for(Int_t j(3); j < 7; j++) {
        tools->CreateOutputList(TString(Form("SVD_kreg_%i", j)));
        tools->SetSVDReg(j);
        tools->Make();
    }
   
    // bayesian unfolding, different number of iterations
   for(Int_t k(1); k < 5; k++) {
      tools->SetUnfoldingAlgorithm(AliJetFlowTools::kBayesian);
      tools->SetBayesianIter(k);
      tools->CreateOutputList(TString(Form("Bayes iter %i", k)));
      tools->Make();
   }
   // finish the unfolding
   // will write the output to file
   tools->Finish();

   // do some post processing (compares unfolding results from different methods, etc)
   tools->PostProcess(TString("SVD kReg 4"));

}

//_____________________________________________________________________________
void Load() {
    gSystem->Load("libTree.so");
    gSystem->Load("libGeom.so");
    gSystem->Load("libVMC.so");
    gSystem->Load("libPhysics");

    gSystem->Load("libSTEERBase.so");
    gSystem->Load("libESD.so");
    gSystem->Load("libAOD.so");
    gSystem->Load("libANALYSIS.so");
    gSystem->Load("libANALYSISalice.so");

    gSystem->Load("libEMCALUtils.so");
    gSystem->Load("libPHOSUtils.so");
    gSystem->Load("libCGAL.so");
    gSystem->Load("libfastjet.so");
    gSystem->Load("libsiscone.so");
    gSystem->Load("libSISConePlugin.so");

    gSystem->Load("libCORRFW.so");
    gSystem->Load("libPWGTools.so");
    gSystem->Load("libJETAN.so");
    gSystem->Load("libFASTJETAN.so");
    gSystem->Load("libPWGJE.so");

    // include paths, necessary for compilation
    gSystem->AddIncludePath("-Wno-deprecated");
    gSystem->AddIncludePath("-I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_ROOT/EMCAL");
    gSystem->AddIncludePath("-I$ALICE_ROOT/PWGDQ/dielectron -I$ALICE_ROOT/PWGHF/hfe");
    gSystem->AddIncludePath("-I$ALICE_ROOT/JETAN -I$ALICE_ROOT/JETAN/fastjet");

    // attempt to load RooUnfold libraries, 
    gSystem->Load("/home/redmer/Documents/CERN/alice/BUILDS/ROOUNFOLD/RooUnfold-1.1.1/libRooUnfold.so");
    gSystem->AddIncludePath("-I/home/redmer/Documents/CERN/alice/BUILDS/ROOUNFOLD/RooUnfold-1.1.1/src/");
    // compile unfolding class
    
    gROOT->LoadMacro("$ALICE_ROOT/PWG/FLOW/Tasks/AliJetFlowTools.cxx+");
}
//_____________________________________________________________________________
