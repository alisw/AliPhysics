
//---------------------------------------------------
// Macro to do analysis with AliJCORRANTask 
// Can be executed with Root and AliRoot
//
// ALICE Jyvaskyla group 
//
// last change 20th Jun 2010 FK
//-------------------------------------------------


const TString kInputData = "ESD";
const TString kJCORRANInputFormat = "ESD"; // ESD, AOD
const Bool_t  kMC = kFALSE; //With real data kMC = kFALSE, MC data kMC =kTRUE

AliJCORRANTask* AddTaskJCORRAN(const char* aodName="jcorran.root", const char* addPhysSelection="$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C")
{
    //--------------------------------------
    // Make the analysis manager
    //-------------------------------------
    AliAnalysisManager *mgr  = AliAnalysisManager::GetAnalysisManager();
    if(!mgr){
      ::Error("AddTaskJCORRAN", "No analysis manager to connect to.");
      return NULL;
    }

    if(!mgr->GetInputEventHandler()){
       ::Error("AddTaskJets", "This task requires an input event handler");
       return NULL;
    }
 
    //AliAODHandler *aodH = dynamic_cast<AliAODHandler*>(AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler());
    //aodH->SetCreateNonStandardAOD();

    //-------------------------------------------------------------------------
    ///           T R A C K     S E L E C T I O N

    AliESDtrackCuts* esdTrackCuts = new AliESDtrackCuts("AliESDtrackCuts", "Loose");
    esdTrackCuts->SetMinNClustersTPC(70);
    esdTrackCuts->SetMaxChi2PerClusterTPC(4.0);
    esdTrackCuts->SetRequireTPCRefit(kTRUE);
    esdTrackCuts->SetRequireITSRefit(kTRUE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMaxDCAToVertexXY(3.5);
    esdTrackCuts->SetMaxDCAToVertexZ(3.5);
    esdTrackCuts->SetMaxDCAToVertexXYPtDep("0.0350+0.0420/pt^0.9");



    //-------------------------------------------------------------------------
    //       E V E N T     S E L E C T I O N

    Int_t    downscaling     = 1;//20;  //downscaling of normal events
    Double_t lowerCutOnLPmom =  0;//3;  // 3 GeV
    Double_t lowerCutOnLeadingCaloClusterE = 1.; //GeV
    Double_t lowerCutOnCaloClusterE = 0.2; //select only clusters above this energy

    if(addPhysSelection){ 
      gROOT->LoadMacro(addPhysSelection);
      AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection();
      AliPhysicsSelection* physSele = physSelTask->GetPhysicsSelection();
    }
    //---------------------------------------------------------------------------
    //            J C O R R A N    T A S K 

    AliJCORRANTask *jctask = new AliJCORRANTask("PWG4JCORRANTask",kJCORRANInputFormat); 
    jctask->SetESDtrackCuts(esdTrackCuts);
    jctask->SetDownScalingOfMB(downscaling);
    jctask->SetLeadingPaticleMomCut(lowerCutOnLPmom);
    jctask->SetLowerCutOnCaloClusterE(lowerCutOnCaloClusterE);
    jctask->SetLowerCutOnLeadingCaloClusterE(lowerCutOnLeadingCaloClusterE);
    jctask->SetRealOrMC(kMC); 
    jctask->SetOutputAODName(aodName); 
    jctask->SetDebugLevel(1);

    if(!kMC && addPhysSelection){
      //Apply offline trigger selection by AliPhysicsSelectionTask
      jctask->SelectCollisionCandidates(AliVEvent::kMB | AliVEvent::kHighMult);
    }
    
    mgr->AddTask(jctask);

    // Create containers for input:
    AliAnalysisDataContainer *cinput0 = mgr->GetCommonInputContainer(); 

    //connect input to JCORRAN task
    mgr->ConnectInput(jctask, 0, cinput0);

    
    return jctask;
}


