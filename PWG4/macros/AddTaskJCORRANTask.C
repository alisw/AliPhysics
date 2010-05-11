
//---------------------------------------------------
// Macro to do analysis with AliJCORRANTask 
// Can be executed with Root and AliRoot
//
// ALICE Jyvaskyla group 
//
//-------------------------------------------------


const TString kInputData = "ESD";
const TString kJCORRANInputFormat = "ESD"; // ESD, AODout, AODin


AliJCORRANTask* AddTaskJCORRAN()
{
    //--------------------------------------
    // Make the analysis manager
    //-------------------------------------
    AliAnalysisManager *mgr  = AliAnalysisManager::GetAnalysisManager();
    if(!mgr){
      ::Error("AddTaskJCORRAN", "No analysis manager to connect to.");
      return NULL;
    }

    if (!mgr->GetInputEventHandler()) {
       ::Error("AddTaskJets", "This task requires an input event handler");
       return NULL;
    }


    //-------------------------------------------------------------------------
    //Define task, put here any other task that you want to use.
    //-------------------------------------------------------------------------
    Int_t     downscaling     = 20;  //downscaling of normal events
    Double_t  lowerCutOnLPmom =  2;  // 3 GeV
    Double_t  lowerCutOnLeadingCaloClusterE = 1.; //GeV

    //           T R A C K     S E L E C T I O N

    // Apply loose track cuts   //FK// 
    AliESDtrackCuts* esdTrackCutsLoose = new AliESDtrackCuts("AliESDtrackCuts", "Loose");
    esdTrackCutsLoose->SetMinNClustersTPC(60);
    esdTrackCutsLoose->SetMaxChi2PerClusterTPC(4.0);
    esdTrackCutsLoose->SetRequireTPCRefit(kTRUE);
    //esdTrackCutsLoose->SetMaxNsigmaToVertex(3);
    //esdTrackCutsLoose->SetRequireSigmaToVertex(kTRUE);
    esdTrackCutsLoose->SetAcceptKinkDaughters(kFALSE);
    esdTrackCutsLoose->SetMaxDCAToVertexXY(3.5);
    esdTrackCutsLoose->SetMaxDCAToVertexZ(3.5);
    // hard
    //AliESDtrackCuts* esdTrackCutsHard = new AliESDtrackCuts("AliESDtrackCuts", "Hard");
    //esdTrackCutsHard->SetMinNClustersTPC(100);
    //esdTrackCutsHard->SetMaxChi2PerClusterTPC(2.0);
    //esdTrackCutsHard->SetMaxCovDiagonalElements(2,2,0.5,0.5,2);
    //esdTrackCutsHard->SetRequireTPCRefit(kTRUE);
    //esdTrackCutsHard->SetMaxNsigmaToVertex(2);
    //esdTrackCutsHard->SetRequireSigmaToVertex(kTRUE);
    //esdTrackCutsHard->SetAcceptKinkDaughters(kFALSE);

    //---------------------------------------------------------------------------
   
    //       E V E N T     S E L E C T I O N

    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
    AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection();
    AliPhysicsSelection* physSele = physSelTask->GetPhysicsSelection();

    //tag as MB different trigger than CINT1B
    //physSele->AddCollisionTriggerClass("+CSMBB-ABCE-NOPF-ALL"); // Put here trigger to be selected, default is CINT1B 

    AliJCORRANTask *jctask = new AliJCORRANTask("PWG4JCORRANTask",kJCORRANInputFormat, esdTrackCutsLoose, downscaling, lowerCutOnLPmom, lowerCutOnLeadingCaloClusterE); //FK//
    jctask->SetDebugLevel(1);
    mgr->AddTask(jctask);

    jctask->SelectCollisionCandidates();  //Apply offline trigger selection by AliPhysicsSelectionTask

    
    // Create containers for input:
    AliAnalysisDataContainer *cinput0  = mgr->GetCommonInputContainer(); 

    // JCORRAN output containers:
    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("jcorrantree", TTree::Class(),
							      AliAnalysisManager::kOutputContainer, "jcorran.root");
							      
    AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("jcorranQAlist", TList::Class(),
							      AliAnalysisManager::kOutputContainer, "jcorranQA.root");


							      
    // JCORRAN task
    mgr->ConnectInput  (jctask,     0, cinput0  );
    mgr->ConnectOutput (jctask,     1, coutput1 );
    mgr->ConnectOutput (jctask,     2, coutput2 );

    
    
 return jctask;
}


