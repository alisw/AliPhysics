
//---------------------------------------------------
// Macro to do analysis with AliJCORRANTask 
// Can be executed with Root and AliRoot
//
// ALICE Jyvaskyla group 
//
// last change 20th Jun 2010 FK
//-------------------------------------------------


const TString kInputData = "ESD";
const TString kJCORRANInputFormat = "ESD"; // ESD, AODout, AODin
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
    //           T R A C K     S E L E C T I O N

    AliESDtrackCuts* esdTrackCuts = new AliESDtrackCuts("AliESDtrackCuts2009","Standard");
  
    /* // Apply loose track cuts   
    AliESDtrackCuts* esdTrackCutsLoose = new AliESDtrackCuts("AliESDtrackCuts", "Loose");
    esdTrackCutsLoose->SetMinNClustersTPC(60);
    esdTrackCutsLoose->SetMaxChi2PerClusterTPC(4.0);
    esdTrackCutsLoose->SetRequireTPCRefit(kTRUE);
    //esdTrackCutsLoose->SetMaxNsigmaToVertex(3);
    //esdTrackCutsLoose->SetRequireSigmaToVertex(kTRUE);
    esdTrackCutsLoose->SetAcceptKinkDaughters(kFALSE);
    esdTrackCutsLoose->SetMaxDCAToVertexXY(3.5);
    esdTrackCutsLoose->SetMaxDCAToVertexZ(3.5);
    */

    //---------------------------------------------------------------------------
    //       E V E N T     S E L E C T I O N

    Int_t     downscaling     = 20;  //downscaling of normal events
    Double_t  lowerCutOnLPmom =  3;  //select all events with a particle above momentum  
    Double_t  lowerCutOnLeadingCaloClusterE = 1.; //GeV   select all events with a calo cluster above the energy   
    Double_t  lowerCutOnCaloClusterE = 0.2; // GeV  store only calo clusters above this energy   

    if(addPhysSelection){ 
      gROOT->LoadMacro(addPhysSelection);
      AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection();
      AliPhysicsSelection* physSele = physSelTask->GetPhysicsSelection();
    
    //tag as MB different trigger than CINT1B
    //physSele->AddCollisionTriggerClass("+CSMBB-ABCE-NOPF-ALL"); // Put here trigger to be selected, default is CINT1B 
    }
    //---------------------------------------------------------------------------
    //            J C O R R A N    T A S K 

    AliJCORRANTask *jctask = new AliJCORRANTask("PWG4JCORRANTask",kJCORRANInputFormat); 
    jctask->SetDebugLevel(1);
    /*jctask->SetAliESDtrackCuts(esdTrackCutsLoose); */
    jctask->SetAliESDtrackCuts(esdTrackCuts->GetStandardITSTPCTrackCuts2009());
    jctask->SetDownscaling(downscaling);
    jctask->SetLowerCutOnLPMom(lowerCutOnLPmom);
    jctask->SetLowerCutOnLeadingCaloClusterE(lowerCutOnLeadingCaloClusterE);
    jctask->SetLowerCutOnCaloClusterE(lowerCutOnCaloClusterE);
    jctask->SetRealOrMC(kMC); //flags whether the input are ESDs from real  exp or MonteCarlo 
    jctask->SetOutputAODName(aodName); 
    mgr->AddTask(jctask);

    jctask->SelectCollisionCandidates();  //Apply offline trigger selection by AliPhysicsSelectionTask

    
    // Create containers for input:
    AliAnalysisDataContainer *cinput0 = mgr->GetCommonInputContainer(); 

    //connect input to JCORRAN task
    mgr->ConnectInput(jctask, 0, cinput0);

    
    return jctask;
}


