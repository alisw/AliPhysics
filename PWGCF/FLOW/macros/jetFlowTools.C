void jetFlowTools() {
    // load and compile the libraries
    Load();

   // read detector response from output of matching taks
    // AliAnalysisTaskJetMatching
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
    TFile f("AnalysisResults.root");
    if(f.IsZombie()) {
        printf(" > read error ! < \n");
        return;
    }
    TList* l = (TList*)f.Get("RhoVnMod_R04_kCombined_Jet_AKTChargedR040_PicoTracks_pT0150_Rho_TPC_PWGJE");
    if(!l) {
        printf(" > failed to find output list ! \n");
        return;
    }
    const Double_t ptBins[] = {20, 25, 30, 35,  40,  45,  50, 55, 60, 70, 80, 90, 100};
    BinsTrue = new TArrayD(sizeof(ptBins)/sizeof(ptBins[0]), ptBins);
    Double_t binsY[81];
    for(Int_t i(0); i < 81; i++) binsY[i] = (double)(30+i);
    BinsRec = new TArrayD(sizeof(binsY)/sizeof(binsY[0]), binsY);

    // create an instance of the Tools class
    AliJetFlowTools* tools = new AliJetFlowTools();
    // set some common variables
    tools->SetCentralityBin(2);
    tools->SetDetectorResponse(detres);
    tools->SetBinsTrue(BinsTrue);
    tools->SetBinsRec(BinsRec);
 
    // connect input
    tools->SetInputList(l);


    // unfold using different parameters
    tools->SetUnfoldingAlgorithm(AliJetFlowTools::kSVD);
    tools->SetBeta(0.01);

    tools->CreateOutputList(TString("R04_kCombined_SVD_d3"));
    tools->SetSVDDraw(3);
    tools->Make();

    tools->CreateOutputList(TString("R04_kCombined_SVD_d4"));
    tools->SetSVDDraw(4);
    tools->Make();
    
    tools->CreateOutputList(TString("R04_kCombined_SVD_d5"));
    tools->SetSVDDraw(5);
    tools->Make();

    // unfold using chi2 method
    tools->SetUnfoldingAlgorithm(AliJetFlowTools::kChi2);
    tools->CreateOutputList(TString("beta_01"));
    tools->SetBeta(0.01);
    tools->Make();

    tools->CreateOutputList(TString("beta_005"));
    tools->SetBeta(0.05);
    tools->Make();


    // finish the unfolding
    tools->Finish();
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
    
    gROOT->LoadMacro("$ALICE_ROOT/PWG/FLOW/Tasks/AliJetFlowTools.cxx++g");
}
//_____________________________________________________________________________
