// Macro for analisys task of preselected-central-diffractive events
//------------------------------------------------------------------
// When there is no time to wait for legotrain results,
// please use this macro to run jobs on grid for preselected events
//------------------------------------------------------------------
// Author: Taesoo Kim, Beomkyu Kim
// email:  kimb@cern.ch
//

// good-run list of LHC10b period : 33, 7 TeV
const int LHC10bRuns[] = {117222,117220,117116,117112,117099,117063,117060,117059,117053,117052,117050,117048,116645,116643,116574,116571,116562,116403,116402,115521,115401,115399,115393,115345,115335,115193,115186,114931,114930,114924,114918,114798,114786};
// good-run list of LHC10c period : 31, 7 TeV
const int LHC10cRuns[] = {120829,120825,120824,120823,120822,120821,120820,120758,120750,120741,120671,120617,120616,120505,120503,120244,120079,120076,120073,120072,120069,120067,119849,119846,119845,119844,119842,119841,119163,119161,119159};
// good-run list of LHC10d period : 55, 7 TeV
const int LHC10dRuns[]= {126432,126425,126424,126422,126409,126408,126407,126406,126405,126404,126403,126359,126352,126351,126285,126284,126283,126168,126167,126160,126158,126097,126090,126088,126082,126081,126078,126073,126008,126007,126004,125855,125851,125850,125849,125848,125847,125844,125843,125842,125633,125632,125630,125628,125296,125295,125134,125100,125097,125085,125083,125023,124751,122375,122374};
// good-run list of LHC10e period : 113, 7 TeV
const int LHC10eRuns[] ={130850,130848,130847,130844,130842,130840,130834,130804,130803,130802,130799,130798,130795,130793,130704,130696,130628,130623,130621,130620,130609,130608,130601,130524,130520,130519,130517,130481,130480,130479,130375,130356,130354,130343,130342,130178,130172,130157,130151,130149,129966,129962,129961,129960,129959,129744,129742,129738,129736,129735,129734,129729,129726,129725,129723,129659,129653,129652,129651,129650,129647,129641,129639,129587,129586,129540,129536,129528,129527,129525,129524,129523,129521,129520,129514,129513,129512,129042,128913,128855,128853,128850,128843,128836,128835,128834,128833,128824,128823,128820,128819,128778,128777,128678,128677,128621,128615,128611,128609,128605,128596,128594,128592,128590,128582,128506,128504,128503,128498,128495,128494,128486,128366};
// good-run list of LHC15f period : 64, 13 TeV
const int LHC15fRuns[]={225011,225016,225026,225031,225035,225037,225041,225043,225050,225051,225052,225105,225106,225305,225307,225310,225313,225314,225315,225322,225576,225578,225579,225580,225582,225586,225587,225589,225709,225710,225716,225717,225763,225768,226062,226085,226170,226175,226176,226177,226183,226210,226225,226444,226445,226466,226468,226472,226495,226500,226532,226543,226551,226554,226569,226573,226591,226593,226596,226600,226602,226603,226605,226606};





class AliAnalysisGrid;
void runCD(
      const char *taskname = "CD"
    , const char *option = "LHC10b" // LHC10b or c or d or e or LHC15f
    , const char *gridmode = "full" // or "terminate" to merge
    )
{
    
    // add aliroot indlude path
    gROOT->ProcessLine(Form(".include %s/include",gSystem->ExpandPathName("$ROOTSYS")));
    gROOT->ProcessLine(Form(".include %s/include",gSystem->ExpandPathName("$ALICE_ROOT")));
    gROOT->ProcessLine(Form(".include %s/include",gSystem->ExpandPathName("$ALICE_PHYSICS")));

    // analysis manager
    AliAnalysisManager* mgr = new AliAnalysisManager(Form("%s%s",taskname,option));

    
    // create the alien handler and attach it to the manager
    AliAnalysisAlien *plugin = new AliAnalysisAlien();
    plugin->SetRunMode(gridmode);
    plugin->SetAPIVersion("V1.1x");
    plugin->SetAliPhysicsVersion("vAN-20160126-1");
    plugin->SetDropToShell(0);
    plugin->SetRunPrefix("000");
    plugin->SetNrunsPerMaster(1);
    plugin->SetOutputToRunNo();
    
    TString foption = option;
    if (foption.Contains("LHC10b")){
        for (int i=0; i<33; i++) plugin->AddRunNumber(LHC10bRuns[i]);
        plugin->SetGridDataDir("/alice/cern.ch/user/k/kimb/PreselectionLHC10b/out");
    }
    if (foption.Contains("LHC10c")){
        for (int i=0; i<31; i++) plugin->AddRunNumber(LHC10cRuns[i]);
        plugin->SetGridDataDir("/alice/cern.ch/user/k/kimb/PreselectionLHC10c/out");
    }
    if (foption.Contains("LHC10d")){
        for (int i=0; i<55; i++) plugin->AddRunNumber(LHC10dRuns[i]);
        plugin->SetGridDataDir("/alice/cern.ch/user/k/kimb/PreselectionLHC10d/out");
    }
    if (foption.Contains("LHC10e")){
        for (int i=0; i<113; i++) plugin->AddRunNumber(LHC10eRuns[i]);
        plugin->SetGridDataDir("/alice/cern.ch/user/k/kimb/PreselectionLHC10e/out");
    }
    if (foption.Contains("LHC15f")){
        for (int i=0; i<64; i++) plugin->AddRunNumber(LHC15fRuns[i]);
        plugin->SetGridDataDir("/alice/cern.ch/user/k/kimb/PreselectionLHC15f/out");
    }
    plugin->SetDataPattern("*ESDs.root"); 
    
    plugin->SetGridWorkingDir(Form("%s%s",taskname,option));
    plugin->SetGridOutputDir("out");
    plugin->AddIncludePath("-I$ALICE_ROOT/include  -I$ALICE_ROOT/lib -I$ALICE_PHYSICS/include -I$ALICE_PHYSICS/lib -I$ALICE_PHYSICS/OADB/macros" );
    plugin->SetDefaultOutputs(kFALSE);
    plugin->SetOutputFiles("AnalysisResults.root tree.root");
    plugin->SetSplitMaxInputFileNumber(100);
    plugin->SetMasterResubmitThreshold(90);

    // Optionally set time to live (default 30000 sec)
    plugin->SetTTL(10000);
    // Optionally set input format (default xml-single)
    plugin->SetInputFormat("xml-single");
    // Optionally modify the name of the generated JDL (default analysis.jdl)
    plugin->SetJDLName(Form("%s%s.jdl",taskname,option));
    // Optionally modify the executable name (default analysis.sh)
    plugin->SetExecutable(Form("%s%s.sh",taskname,option));
    // Optionally modify job price (default 1)
    plugin->SetPrice(1);
    // Optionally modify split mode (default 'se')
    plugin->SetSplitMode("se");

    
    mgr->SetGridHandler(plugin);
    AliInputEventHandler* esdH = new AliESDInputHandler();
    esdH->SetNeedField(1);
    mgr->SetInputEventHandler(esdH);
    

    // Add PID task manager if we use PID----------------------------------------
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
    AliAnalysisTaskPIDResponse *pidResponseTask = AddTaskPIDResponse(0);
    if(!pidResponseTask) { Printf("no pidResponseTask"); return; }
    //---------------------------------------------------------------------------
   
    AliAnalysisTaskSE *task;
 
    if (foption.Contains("LHC15f")){
        AliInputEventHandler* hdl = (AliInputEventHandler*)AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();
        if (hdl) hdl->SetNeedField(kTRUE); 
        task = new AliAnalysisTaskCDPWA("ESD");
    } else {
        gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
        AliPhysicsSelectionTask *physSelTask = AddTaskPhysicsSelection(0);
        if(!physSelTask) { Printf("no physSelTask"); return; }

        task = new AliAnalysisTaskCDTree("ESD");
        task->SelectCollisionCandidates(AliVEvent::kMB);
    }

    
    // Create containers for input/output
    AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
    AliAnalysisDataContainer *coutput = mgr->CreateContainer("output", TTree::Class(), AliAnalysisManager::kOutputContainer,"tree.root");
    AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("output2", TList::Class(), AliAnalysisManager::kOutputContainer,"AnalysisResults.root");
    
    mgr->AddTask(task);
    mgr->ConnectInput(task, 0, cinput);
    mgr->ConnectOutput(task, 1, coutput);
    mgr->ConnectOutput(task, 2, coutput2);


    // enable debug printouts
    mgr->SetDebugLevel(2);
    if (!mgr->InitAnalysis()) return;
    mgr->PrintStatus();

    // start analysis
    Printf("Starting Analysis....");
    mgr->StartAnalysis("grid",1234567890,0);
}
