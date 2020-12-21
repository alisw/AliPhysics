SteerAnalysisTaskdNdEtaPP13(const Char_t *inputfilename = NULL, Int_t type = 0, const Char_t *outputfilename = "dNdEtapp13", Int_t maxFiles = kMaxInt, Int_t maxEv = kMaxInt)
{
    
    /* include path for ACLic */
//    gSystem->AddIncludePath("-I$ALICE_ROOT/../build/include");
    gSystem->AddIncludePath("-I$ALICE_ROOT/include");
    gSystem->AddIncludePath("-I$ALICE_PHYSICS/include");
    //    gSystem->AddIncludePath("-I$ALICE_ROOT/include");
    //  gSystem->AddIncludePath("-I$ALICE_ROOT/TOF");
    /* load libraries */
    //    gSystem->Load("libSTEERBase");
    //    gSystem->Load("libESD");
    //    gSystem->Load("libAOD");
    gSystem->Load("libANALYSIS");
    gSystem->Load("libANALYSISalice");
    /* build analysis task class */
    gROOT->LoadMacro("AliITSMultRecBg.cxx+g");
    gROOT->LoadMacro("AliAnalysisTaskdNdEtapp13.cxx+");
    
    if (!inputfilename)
        return;
    
    /* setup input chain */
    TString str = inputfilename;
    const Char_t *filename;
    TChain *chain = new TChain("esdTree");
    if (str.EndsWith(".xml")) {
        TGrid::Connect("alien://");
        Info("", "reading data list from collection:");
        TGridCollection *coll = gGrid->OpenCollection(inputfilename, maxFiles);
        coll->Reset();
        while (coll->Next()) {
            filename = coll->GetTURL();
            Info("", Form("%s", filename));
            chain->Add(filename);
        }
    }
    else if (str.EndsWith(".txt")) {
        Info("", "reading data list from text file:");
        ifstream is(inputfilename);
        Char_t buf[4096];
        while(!is.eof()) {
            is.getline(buf, 4096);
            if (is.eof()) break;
            chain->Add(buf);
            Info("", Form("%s", buf));
        }
        is.close();
    }
    else {
        Info("", "single file:");
        filename = inputfilename;
        Info("", Form("%s", filename));
        chain->Add(filename);
    }
    Info("", Form("chain is ready: %d events", chain->GetEntries()));
    
    /* create analysis manager */
    AliAnalysisManager *mgr = new AliAnalysisManager("dNdEtapp13");
    
    /* define input event handler */
    AliESDInputHandler *esdh = new AliESDInputHandler();
    esdh->SetReadFriends(kFALSE);
    mgr->SetInputEventHandler(esdh);
    
    /* define MC truth event handler */
    if (type == 1) {
        AliMCEventHandler *mch = new AliMCEventHandler();
        mgr->SetMCtruthEventHandler(mch);
    }
    
    /* add Physics Selection */
    gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
    AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection(type != 0);
    
    /* add our task */
    TString outFileName = outputfilename;
    if (type == 0)
      outFileName += ".data.root";
    if (type == 1)
      outFileName += ".mc.root";
    if (type == 2)
      outFileName += ".mcdata.root";
    gROOT->LoadMacro("AddAnalysisTaskdNdEtapp13.C");

    const Char_t *centVar[3] = {"MB", "V0M", "Ref"};
    for (Int_t i = 0; i < 3; i++) {
    
      AliAnalysisTaskdNdEtapp13 *task = AddAnalysisTaskdNdEtapp13(outFileName.Data(), Form("clist_%s", centVar[i]), -5., 5., -10., 10., centVar[i], 1.5, -1., type == 1);
      //      AliAnalysisTaskdNdEtapp13 *task = AddAnalysisTaskdNdEtapp13(outFileName.Data(), Form("clist_%s_L", centVar[i]), -5., 5., 0., 10., centVar[i], 1.5, -1., type == 1);
      //      AliAnalysisTaskdNdEtapp13 *task = AddAnalysisTaskdNdEtapp13(outFileName.Data(), Form("clist_%s_R", centVar[i]), -5., 5., -10., 0., centVar[i], 1.5, -1., type == 1);
      
      //      AliAnalysisTaskdNdEtapp13 *task = AddAnalysisTaskdNdEtapp13(outFileName.Data(), Form("clist_%s_W", centVar[i]), -3., 3., -20., 20., centVar[i], 1.5, -1., type == 2);
      //	    AliAnalysisTaskdNdEtapp13 *task = AddAnalysisTaskdNdEtapp13(outFileName.Data(), Form("clist_%s_LL", centVar[i]), -2.0, -0.5, 10., 15., centVar[i], 1.5, -1., type == 2);
      //	    AliAnalysisTaskdNdEtapp13 *task = AddAnalysisTaskdNdEtapp13(outFileName.Data(), Form("clist_%s_RR", centVar[i]), 0.5, 2.0, -15., -10., centVar[i], 1.5, -1., type == 2);
    
    }
    
    /* start analysis */
    mgr->SetDebugLevel(0);
    if (!mgr->InitAnalysis()) return;
    mgr->PrintStatus();
    mgr->StartAnalysis("local", chain, maxEv);
    
    /* create dummy file to tell we are done */
    gSystem->Exec("touch done");
    
}
