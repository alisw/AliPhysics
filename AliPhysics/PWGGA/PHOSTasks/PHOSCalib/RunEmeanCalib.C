void RunEmeanCalib(
                   TString dataset = "Find;"
                   "BasePath=/alice/data/2015/LHC15h/000233219/muon_calo_pass1/%/;"
                   "FileName=root_archive.zip;"
                   "Anchor=AliESDs.root;"
                   "Tree=/esdTree;"
                   "Mode=remote;",  // <-- much faster dataset creation
                   Bool_t usePhysicsSelection = kTRUE,
                   Int_t numEvents = 999999999,
                   Int_t firstEvent = 0
                   ) {
    
    // Not needed on the VAF
    //gEnv->SetValue("XSec.GSI.DelegProxy","2");
    
    TString extraLibs = "ANALYSIS:ANALYSISalice:PWGGAPHOSTasks"; // extraLibs = "ANALYSIS:OADB:ANALYSISalice:CORRFW:OADB:PWGmuon";
    
    TList *list = new TList();
    list->Add(new TNamed("ALIROOT_EXTRA_LIBS", extraLibs.Data()));
    list->Add(new TNamed("ALIROOT_ENABLE_ALIEN", "1"));  // important: creates token on every PROOF worker
    
    // Not needed on the VAF
    //TProof::Mgr("alice-caf.cern.ch")->SetROOTVersion("VO_ALICE@ROOT::v5-34-08");
    
    // Note the difference between CAF and VAF
    //TProof::Open("alice-caf.cern.ch");
    TProof::Open("pod://");
    
    // Check the dataset before running the analysis!
    gProof->ShowDataSet( dataset.Data() );
    //return;  // <-- uncomment this to test search before running the analysis!
    
    // Not needed on the VAF
    //gProof->EnablePackage("VO_ALICE@AliRoot::v5-04-81-AN", list);
    
    // A single AliRoot package for *all* AliRoot versions: new on VAF
    TString aliceVafPar = "/afs/cern.ch/alice/offline/vaf/AliceVaf.par";
    gProof->UploadPackage(aliceVafPar.Data());
    gProof->EnablePackage(aliceVafPar.Data(), list);  // this "list" is the same as always
    
    AliAnalysisManager *mgr  = new AliAnalysisManager("Analysis Train");
    
    AliESDInputHandler *esdH = new AliESDInputHandler();
    mgr->SetInputEventHandler(esdH);
    
    gProof->Load("AliAnalysisTaskEmeanCalib.cxx+");  // DON'T use double '+' when running multiple times: it uselessly recompiles everything!
    //return;
    
    gROOT->LoadMacro("AddTaskEmeanCalib.C");
    AliAnalysisTaskEmeanCalib *simplePtTask = AddTaskEmeanCalib(usePhysicsSelection);
    
    AliCDBManager::Instance()->SetDefaultStorage("local://OCDB");
    AliPHOSCalibData* c = new AliPHOSCalibData(0);
    simplePtTask->SetCalibrations(c);
    
    //gROOT->LoadMacro("$ALICE_PHYSICS/PWGGA/PHOSTasks/PHOS_pPb/AddTaskPHOSPi0pPb.C");
    //AliAnalysisTaskPi0Flow* simplePtTask = AddTaskPHOSPi0pPb("PHOSPi0pp","",0,"V0M",1,-1,100);
    
    //simplePtTask->SetEnablePHOSModule(2, kTRUE);
    //simplePtTask->SetEnablePHOSModule(4, kTRUE);
    //simplePtTask->EnableTOFCut(kFALSE);
    
    if (usePhysicsSelection) {
        simplePtTask->SelectCollisionCandidates(AliVEvent::kAny);
    }
    
    if (!mgr->InitAnalysis()) return;
    mgr->StartAnalysis("proof", dataset, numEvents, firstEvent);
    
}
