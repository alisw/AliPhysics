void readTTree()
{

    // example macro to read data from a ttree and perform a flow analysis using the flow pacakge
    // author: Redmer Alexander Bertens (rbertens@cern.ch)
    // note: this macro can run in ROOT only provided libPWGflowBase is available

    // compile the relevant classes
    // include paths, necessary for compilation
    gSystem->AddIncludePath("-Wno-deprecated");
    gSystem->AddIncludePath("-I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");
 
    // load libraries
    gSystem->Load("libPWGflowBase");

    // comile the encapsulated classes
    gROOT->LoadMacro("../objects/AliFlowTTreeEvent.cxx+");
    gROOT->LoadMacro("../objects/AliFlowTTreeTrack.cxx+");
    gROOT->LoadMacro("../objects/AliFlowEventSimpleFromTTree.cxx+");
    
    TChain* myChain = new TChain("tree");
    myChain->Add("/home/rbertens/Documents/CERN/ALICE_DATA/filtered/000167988.root");
    myChain->Add("/home/rbertens/Documents/CERN/ALICE_DATA/filtered/000168066.root");

    // create pointers for the branches
    AliFlowTTreeEvent* event = 0x0;
    myChain->SetBranchAddress("event", &event);
    TClonesArray* tracks = 0x0;
    myChain->SetBranchAddress("track", &tracks);
    // and an example track
    AliFlowTTreeTrack* firstTrack = 0x0;

    // connection to the flow package
    AliFlowAnalysisWithQCumulants* qc = new AliFlowAnalysisWithQCumulants();
    qc->Init();
    AliFlowTrackSimpleCuts* cutsPOI = new AliFlowTrackSimpleCuts();
    cutsPOI->SetPtMin(0.2);
    cutsPOI->SetPtMin(2.);
    AliFlowTrackSimpleCuts* cutsRP = new AliFlowTrackSimpleCuts();

    // set how many events you want to analyze
    Int_t maxEvents = 10000;
    // event loop
    printf(" > %i events in chain, processing %i of them < \n", myChain->GetEntries(), maxEvents);
    for(Int_t i = 0; i < myChain->GetEntries(); i++) {
        cout << " > Parsing event " << i << "\r"; cout.flush();
        myChain->GetEntry(i);
        // pass info to flow package
        AliFlowEventSimple* flowevent = new AliFlowEventSimpleFromTTree(event, tracks, cutsPOI, cutsRP);
        qc->Make(flowevent);
        delete flowevent;
        maxEvents--;
        if(maxEvents < 1) break;
    }

    // wrap up analysis
    qc->Finish();

    // open a new file which will hold the final results of all methods:
    TString outputFileName = "qc_results.root";
    TFile *outputFile = new TFile(outputFileName.Data(),"RECREATE");
    const Int_t nMethods = 1;
    TString method[nMethods] = {"QC"};
    TDirectoryFile *dirFileFinal[nMethods] = {NULL};
    TString fileName[nMethods];
    for(Int_t i=0; i<nMethods; i++)
    {
        // form a file name for each method:
        fileName[i]+="output";
        fileName[i]+=method[i].Data();
        fileName[i]+="analysis";
        dirFileFinal[i] = new TDirectoryFile(fileName[i].Data(),fileName[i].Data());
    }

    // store the final results
    qc->WriteHistograms(dirFileFinal[0]);
    outputFile->Close();
}

