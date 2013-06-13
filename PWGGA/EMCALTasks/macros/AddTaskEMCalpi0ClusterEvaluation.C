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
    gSystem->Load("libGui.so");
    gSystem->Load("libXMLParser.so");
    gSystem->Load("libCDB.so");
    gSystem->Load("libProof.so");
    gSystem->Load("libRAWDatabase.so");
    gSystem->Load("libRAWDatarec.so");
    gSystem->Load("libSTEERBase.so");
    gSystem->Load("libSTEER.so");
    gSystem->Load("libTRDbase.so");
    gSystem->Load("libTOFbase.so");
    gSystem->Load("libTOFrec.so");
    gSystem->Load("libVZERObase.so");
    gSystem->Load("libVZEROrec.so");
    gSystem->Load("libMinuit.so");
    gSystem->Load("libEMCALUtils.so");
    gSystem->Load("libEMCALraw.so");
    gSystem->Load("libEMCALbase.so");
    gSystem->Load("libEMCALrec.so");
    gSystem->Load("libANALYSIS.so");
    gSystem->Load("libANALYSISalice.so");
    gSystem->Load("libPHOSUtils.so");
    gSystem->Load("libTENDER.so");
    gSystem->Load("libTENDERSupplies.so");
    gSystem->Load("libPWGflowBase.so");
    gSystem->Load("libPWGflowTasks.so");

}

//______________________________________________________
void AddTaskEMCalpi0ClusterEvaluation( UInt_t triggerMaskPbPb = AliVEvent::kEMCEGA )
{

    // load libraries
    load_libraries();

    // Use AliRoot includes to compile our task
    gROOT->ProcessLine(".include $ALICE_ROOT");
    gROOT->ProcessLine(".include $ALICE_ROOT/include");
    gROOT->ProcessLine(".include $ALICE_ROOT/PWG/FLOW/Base");




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
