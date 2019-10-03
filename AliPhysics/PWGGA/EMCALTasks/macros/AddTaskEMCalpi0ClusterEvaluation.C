//Task to run over AOD EMCal Clusters and tender
//Astrid Morreale 2013
//Esd
//___________________________________________________
void load_libraries( void )
{

    // Root libraries
    gSystem->Load("libTree");
    gSystem->Load("libGeom");
    gSystem->Load("libVMC");
    gSystem->Load("libPhysics");
    gSystem->Load("libMinuit");

    // Analysis framework libraries
    gSystem->Load("libSTEERBase");
    gSystem->Load("libESD");
    gSystem->Load("libAOD");
    gSystem->Load("libANALYSIS");
    gSystem->Load("libOADB");
    gSystem->Load("libANALYSISalice");

    // AliRoot libraries
    gSystem->Load("libGui");
    gSystem->Load("libXMLParser");
    gSystem->Load("libCDB");
    gSystem->Load("libProof");
    gSystem->Load("libRAWDatabase");
    gSystem->Load("libRAWDatarec");
    gSystem->Load("libSTEERBase");
    gSystem->Load("libSTEER");
    gSystem->Load("libTRDbase");
    gSystem->Load("libTOFbase");
    gSystem->Load("libTOFrec");
    gSystem->Load("libVZERObase");
    gSystem->Load("libVZEROrec");
    gSystem->Load("libMinuit");
    gSystem->Load("libEMCALUtils");
    gSystem->Load("libEMCALraw");
    gSystem->Load("libEMCALbase");
    gSystem->Load("libEMCALrec");
    gSystem->Load("libANALYSIS");
    gSystem->Load("libANALYSISalice");
    gSystem->Load("libPHOSUtils");
    gSystem->Load("libTender");
    gSystem->Load("libTenderSupplies");
    gSystem->Load("libPWGflowBase");
    gSystem->Load("libPWGflowTasks");

}

//______________________________________________________
void AddTaskEMCalpi0ClusterEvaluation( UInt_t triggerMaskPbPb = AliVEvent::kEMCEGA )
{

    // load libraries
    load_libraries();

    // Use AliRoot includes to compile our task
    gROOT->ProcessLine(".include $ALICE_PHYSICS");
    gROOT->ProcessLine(".include $ALICE_PHYSICS/include");
    gROOT->ProcessLine(".include $ALICE_PHYSICS/PWG/FLOW/Base");




    // analysis manager
    AliAnalysisManager *analysisManager =  AliAnalysisManager::GetAnalysisManager();

    if (!analysisManager)
    {
    ::Error("AddTaskEMCalpi0ClusterEvaluation", "No analysis manager to connect to.");
    return NULL;
    }

  if (!analysisManager->GetInputEventHandler()) {
    ::Error("AddTaskEMCALpi0ClusterEValuation", "This task requires an input event handler");
    return NULL;
   }

  TString type = analysisManager->GetInputEventHandler()->GetDataType();

  AliAnalysisDataContainer *cinput1  = analysisManager->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 = analysisManager->CreateContainer("MassHistos", TList::Class(), AliAnalysisManager::kOutputContainer, "AllMBLHC11h.root");

   // cluster evaluation
    gROOT->LoadMacro("AliEMCalpi0ClusterEvaluationTask.cxx");

    // create task
    AliEMCalpi0ClusterEvaluationTask*clusterEvaluation = new AliEMCalpi0ClusterEvaluationTask( "clusterEvaluation" );
    clusterEvaluation->SelectCollisionCandidates(triggerMaskPbPb);

    // add task to manager
    analysisManager->AddTask(clusterEvaluation);

    gSystem->AddIncludePath("-I$ALICE_ROOT/ANALYSIS ");


    // I/O
    analysisManager->ConnectInput( clusterEvaluation, 0, cinput1 );
    analysisManager->ConnectOutput( clusterEvaluation, 1, coutput1 );


   return clusterEvaluation;


}
