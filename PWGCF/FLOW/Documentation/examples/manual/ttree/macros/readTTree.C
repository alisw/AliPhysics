void readTTree()
{

    // example macro to read data from a ttree and perform a flow analysis using the flow pacakge
    // author: Redmer Alexander Bertens (rbertens@cern.ch)
    // note: this macro can run in ROOT only provided libPWGflowBase.so is available

    // compile the relevant classes
    // include paths, necessary for compilation
    gSystem->AddIncludePath("-Wno-deprecated");
    gSystem->AddIncludePath("-I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_ROOT/EMCAL");
    gSystem->AddIncludePath("-I$ALICE_ROOT/PWGDQ/dielectron -I$ALICE_ROOT/PWGHF/hfe");
    gSystem->AddIncludePath("-I$ALICE_ROOT/JETAN -I$ALICE_ROOT/JETAN/fastjet");
 
    // load libraries
    gSystem->Load("libCore.so");        
    gSystem->Load("libGeom.so");
    gSystem->Load("libVMC.so");
    gSystem->Load("libPhysics.so");
    gSystem->Load("libTree.so");
    gSystem->Load("libPWGflowBase.so");

    // comile the encapsulated classes
    gROOT->LoadMacro("../objects/AliFlowTTreeEvent.cxx+");
    gROOT->LoadMacro("../objects/AliFlowTTreeTrack.cxx+");
    gROOT->LoadMacro("../objects/AliFlowEventSimpleFromTTree.cxx+");
    
    // read the ttree to get some info
    TFile f("myFilteredTree.root");
    TTree* AliFlowTTreeTree = (TTree*)f.Get("tree");

    // create pointers for the branches
    AliFlowTTreeEvent* event = 0x0;
    AliFlowTTreeTree->SetBranchAddress("event", &event);
    TClonesArray* tracks = 0x0;
    AliFlowTTreeTree->SetBranchAddress("track", &tracks);
    // and an example track
    AliFlowTTreeTrack* firstTrack = 0x0;

    // connection to the flow package
    AliFlowAnalysisWithQCumulants* qc = new AliFlowAnalysisWithQCumulants();
    qc->Init();
    AliFlowTrackSimpleCuts* cutsPOI = new AliFlowTrackSimpleCuts();
    cutsPOI->SetPtMin(0.2);
    cutsPOI->SetPtMin(2.);
    AliFlowTrackSimpleCuts* cutsRP = new AliFlowTrackSimpleCuts();

    // event loop
    for(Int_t i = 0, maxEvents = 1000; i < AliFlowTTreeTree->GetEntries(); i++) {
        cout << " > Parsing event " << i << "\r"; cout.flush();
        AliFlowTTreeTree->GetEntry(i);
        // pass info to flow package
        AliFlowEventSimple* flowevent = new AliFlowEventSimpleFromTTree(event, tracks, cutsPOI, cutsRP);
        qc->Make(flowevent);
        delete flowevent;
        maxEvents--;
        if(maxEvents < 1) break;
    }

    // close file 
    f.Close();

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

