void jetV3tools() {
    // load and compile the libraries
    Load();


    const Int_t SVDbestValueIn = 5;
    const Int_t SVDbestValueOut = 4;
    const Double_t bestBetaIn = 1.25;
    const Double_t bestBetaOut = 1.25;

    Bool_t runUnfolding = 0;
    Bool_t doSystematics = (!runUnfolding);

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

    // get a TList from the AliAnalysisJetV3 task
    TFile f("AnalysisResults.root");
    if(f.IsZombie()) {
        printf(" > read error ! < \n");
        return;
    }
    TList* l = (TList*)f.Get("RhoVnMod_R02_kCombined_Jet_AKTChargedR020_PicoTracks_pT0150_pt_scheme_TpcRho_ExLJ_TPC_PWGJE");
    if(!l) {
        printf(" > failed to find output list ! \n");
        return;
    }
    // create an instance of the Tools class
    AliJetFlowTools* tools = new AliJetFlowTools();

    // set some common variables
    tools->SetCentralityBin(1);
    tools->SetDetectorResponse(detres);
    // set the true (unfolded) bins
    Double_t binsTrue[] = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 170};
    tools->SetBinsTrue(new TArrayD(sizeof(binsTrue)/sizeof(binsTrue[0]), binsTrue));
    // set the measured (folded) bins
    Double_t binsRec[] = {20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90};
    tools->SetBinsRec(new TArrayD(sizeof(binsRec)/sizeof(binsRec[0]), binsRec));
    // connect input
    tools->SetInputList(l);

    if(runUnfolding) {
        tools->SetUnfoldingAlgorithm(AliJetFlowTools::kNone);
        tools->CreateOutputList(TString("DONOTHING"));
        tools->Make();
        // optional: smoothen the spectrum
        tools->SetSaveFull(kTRUE);
        tools->SetSmoothenPrior(kFALSE, 50, 250, 70, kFALSE);
        // optional: normalize the spectrum
        tools->SetUseDetectorResponse(kTRUE);
        // optional: save a lot of raw output
        tools->SetExLJDpt(kTRUE);
        // do some /chi2 unfolding
        tools->SetUnfoldingAlgorithm(AliJetFlowTools::kChi2);
        // first step: fishnet, see what good unfolding regularizations are
        Double_t beta[] = {.001, .01 .1, .25, 1.25};
        Int_t kReg[] = {9, 8, 7, 6, 5, 4, 3, 5};
        // for out 
        Double_t betaOut[] = {.001, .01, .1, .25, 1.25};
        Int_t kRegOut[] = {9, 8, 7, 6, 5, 4, 3, 4};
        for(Int_t b(0); b < sizeof(beta)/sizeof(beta[0]); b++) {
            tools->CreateOutputList(TString(Form("#beta = %.4f", beta[b])));
            tools->SetBeta(beta[b], betaOut[b]);
            tools->Make();
        }
        tools->SetPrior(AliJetFlowTools::kPriorChi2);
        for(Int_t j(0); j < sizeof(kReg)/sizeof(kReg[0]); j++) {
            // do some SVD unfolding
            tools->SetUnfoldingAlgorithm(AliJetFlowTools::kSVD);
            tools->CreateOutputList(TString(Form("SVD kReg %i", kReg[j])));
            tools->SetSVDReg(kReg[j]); 
            tools->Make();
        }
        // after fishnet: 
        // here we change the pt binning, using optimal svd and beta values
        tools->SetBeta(bestBetaIn, bestBetaOut);
        tools->SetSVDReg(SVDbestValueIn, SVDbestValueOut);
        // ###### change the true (unfolded) binning ########
        Double_t binsTrue2[] = {5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 170};
        tools->SetBinsTrue(new TArrayD(sizeof(binsTrue2)/sizeof(binsTrue2[0]), binsTrue2));
        tools->SetUnfoldingAlgorithm(AliJetFlowTools::kSVD);
        tools->CreateOutputList(TString("true bin removed"));
        tools->Make();
        // revert the true bin settings to their default ones
        tools->SetBinsTrue(new TArrayD(sizeof(binsTrue)/sizeof(binsTrue[0]), binsTrue));
        // ####### change the measured binning
        // remove a bin at low pt
        Double_t binsRech[] = {25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90};
        tools->SetBinsRec(new TArrayD(sizeof(binsRech)/sizeof(binsRech[0]), binsRech));
        tools->SetUnfoldingAlgorithm(AliJetFlowTools::kSVD);
        tools->CreateOutputList(TString("measured bin removed"));
        tools->Make();
        // add a bin at low pt
        Double_t binsRecl[] = {15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90};
        tools->SetBinsRec(new TArrayD(sizeof(binsRecl)/sizeof(binsRecl[0]), binsRecl));
        tools->SetUnfoldingAlgorithm(AliJetFlowTools::kSVD);
        tools->CreateOutputList(TString("measured bin added"));
        tools->Make();
        
        // revert to the original binning
        tools->SetBinsTrue(new TArrayD(sizeof(binsTrue)/sizeof(binsTrue[0]), binsTrue));
        tools->SetBinsRec(new TArrayD(sizeof(binsRec)/sizeof(binsRec[0]), binsRec));
 
        //  unfold using different tracking efficiency
        TString drInputName95 = "/home/rbertens/Documents/CERN/jet-flow/RESPONSE/R02/95_pct_efficiency/5gev_leading_track_bias/response.root";
        TFile drInput95(drInputName95.Data());
        if(drInput95.IsZombie()) return;
        TH2D* detres95 = (TH2D*)drInput95.Get("detector_response");
        tools->SetDetectorResponse(detres95);
        tools->CreateOutputList(TString("diced response"));
        tools->Make();

        // switch back to the original detector response
        tools->SetDetectorResponse(detres);

 
        // now do the svd unfolding with a pythia spectrum as a prior
        tools->SetUnfoldingAlgorithm(AliJetFlowTools::kSVD);
        gROOT->LoadMacro("/home/rbertens/Documents/CERN/jet-flow/RESPONSE/pythia.C");
        tools->SetPrior(AliJetFlowTools::kPriorPythia, pythia());
        tools->CreateOutputList(TString("pythia prior"));
        tools->Make();
        
        
        // ########### systematics are done ################
        // write the output to file
        tools->Finish();
        // do some postprocessing
        tools->PostProcess(TString("SVD kReg 4"), 4, 20, 100);
    } // end of run unfolding

    if(doSystematics) {
        /**
         * evaluate systematics              */
        // first element of array should point to the nominal value
        const Int_t reg[] = {9, 4, 8, 10, 16};
        TArrayI* regArray = new TArrayI(sizeof(reg)/sizeof(reg[0]), reg);
        const Int_t rec[] = {9, 12};
        TArrayI* recArray = new TArrayI(sizeof(rec)/sizeof(rec[0]), rec);
        const Int_t tru[] = {9, 13, 14};
        TArrayI* truArray = new TArrayI(sizeof(tru)/sizeof(tru[0]), tru);


        // place holder pointers. these will be assigned by using pointer references in the relevant functions
        //
        // for the nominal points
        TH1D*                   nominalRatio    (0x0);
        TGraphErrors*           nominalV3       (0x0);

        // for the shape uncertainty 
        TGraphAsymmErrors*      shapeRatio      (0x0);
        TGraphAsymmErrors*      shapeV3         (0x0);

        // for the correlated uncertainty
        TGraphAsymmErrors*      corrRatio       (0x0);
        TGraphAsymmErrors*      corrV3          (0x0);
        // get the actual values

        tools->GetNominalValues(
                nominalRatio,
                nominalV3,
                regArray,       // doesn't matter which array is passed, as long as first element points to nominal value
                regArray);

        tools->GetShapeUncertainty(
                shapeRatio,
                shapeV3,
                regArray,        // systematics from regularization      
                regArray,       
                recArray,            // from true spectrum variation
                recArray,    
                truArray,            // from rec spectrum variation
                truArray,
                4, 
                20, 
                100);

        const Int_t cor[] = {9, 15};
        TArrayI* corArray = new TArrayI(sizeof(cor)/sizeof(cor[0]), cor);

        tools->GetCorrelatedUncertainty(
                corrRatio,
                corrV3,
                corArray,        // correlated systematics
                corArray,       
                "diced respose", // name of systematic source
                4, 
                20, 
                100);

        
        
        
        
        using namespace AliJetFlowTools;
        Double_t rangeLow(20.);
        Double_t rangeHigh(90.);

        TFile FinalResults = TFile("FinalResults.root", "RECREATE");
        // combine the final results and write them to a file
        TCanvas* full = new TCanvas("full", "full");
        full->Divide(2);
        full->cd(1);
        AliJetFlowTools::Style(gPad, "PEARSON");
        // shape uncertianty, full boxes
        Style(shapeRatio, kYellow, AliJetFlowTools::kRatio);
        shapeRatio->SetTitle("shape uncertainty");
        shapeRatio->GetXaxis()->SetRangeUser(rangeLow, rangeHigh);
        shapeRatio->GetYaxis()->SetRangeUser(.7, 2.2);
        shapeRatio->DrawClone("a2");

        // correlated uncertainty, open boxes
        Style(corrRatio, kGray, AliJetFlowTools::kRatio);
        corrRatio->SetTitle("correlated uncertainty");
        corrRatio->SetFillStyle(0);
        corrRatio->SetFillColor(kWhite);
        corrRatio->DrawClone("5");
        
        // ratio itself
        Style(nominalRatio, kBlack, AliJetFlowTools::kRatio);
        nominalRatio->DrawCopy("same E1");
        nominalRatio->SetTitle("in / out of plane jet yield");

        AddTPaveText("0-10 \% cc Pb-Pb #sqrt{s_{NN}} = 2.76 TeV");
        AliJetFlowTools::AddLegend(gPad, kTRUE); 
        full->cd(2);
        AliJetFlowTools::Style(gPad, "PEARSON");
        
        // shape uncertainto on v3
        Style(shapeV3, kYellow, AliJetFlowTools::kV3);
        shapeV3->SetTitle("shape uncertainty");
        shapeV3->GetXaxis()->SetRangeUser(rangeLow, rangeHigh);
        shapeV3->GetYaxis()->SetRangeUser(-.5, 1.);
        shapeV3->DrawClone("a2");

        // correlated uncertainty
        Style(corrV3, kGray, AliJetFlowTools::kV3);
        corrV3->SetFillColor(kWhite);
        corrV3->SetLineStyle(0);
        corrV3->SetFillStyle(0);
        corrV3->SetTitle("correlated uncertainty");
        corrV3->DrawClone("5");

        // v3 itself
        Style(nominalV3, kBlack, AliJetFlowTools::kV3);
        nominalV3->SetTitle("jet #it{v}_{2}");
        nominalV3->SetFillColor(kWhite);
        nominalV3->DrawClone("same E1");

        AddTPaveText("0-10 \% cc Pb-Pb #sqrt{s_{NN}} = 2.76 TeV");
        AliJetFlowTools::AddLegend(gPad, kTRUE);
        gStyle->SetTitleStyle(0);
        gStyle->SetStatStyle(0);

        full->Write();
        FinalResults.Close();
    }    
}

//_____________________________________________________________________________
void Load() {
    gSystem->Load("libTree");
    gSystem->Load("libGeom");
    gSystem->Load("libVMC");
    gSystem->Load("libPhysics");

    gSystem->Load("libSTEERBase");
    gSystem->Load("libESD");
    gSystem->Load("libAOD");
    gSystem->Load("libANALYSIS");
    gSystem->Load("libANALYSISalice");
  
    gSystem->Load("libEMCALUtils.so");
    gSystem->Load("libPHOSUtils.so");

    //fastjet 3.0.3 
    gSystem->Load("libCGAL.so");
    gSystem->Load("libfastjet.so");
    gSystem->Load("libsiscone.so");
    gSystem->Load("libsiscone_spherical.so");
    gSystem->Load("libfastjetplugins.so");
    gSystem->Load("libfastjetcontribfragile.so");

    gSystem->Load("libJETAN.so");
    gSystem->Load("libFASTJETAN.so");

    // include paths, necessary for compilation
    gSystem->AddIncludePath("-Wno-deprecated");
    gSystem->AddIncludePath("-I$ALICE_PHYSICS/include");
    gSystem->AddIncludePath("-I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_ROOT/EMCAL");
    gSystem->AddIncludePath("-I$ALICE_ROOT/PWGDQ/dielectron -I$ALICE_ROOT/PWGHF/hfe");
    gSystem->AddIncludePath("-I$ALICE_ROOT/JETAN -I$ALICE_ROOT/JETAN/fastjet");

    // attempt to load RooUnfold libraries, 
    gSystem->Load("/home/rbertens/Documents/CERN/alice/BUILDS/ROOUNFOLD/RooUnfold-1.1.1/libRooUnfold");
    gSystem->AddIncludePath("-I/home/rbertens/Documents/CERN/alice/BUILDS/ROOUNFOLD/RooUnfold-1.1.1/src/");
    // compile unfolding class (only if there are local changes or the .o is not found)
    gROOT->LoadMacro("$PWGJE/EMCALJetTasks/UserTasks/AliJetFlowTools.cxx+");
}
//_____________________________________________________________________________
