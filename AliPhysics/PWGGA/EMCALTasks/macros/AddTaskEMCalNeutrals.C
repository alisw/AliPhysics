AliEMCalpi0Task* AddTaskEMCalNeutrals( void )
{

    /// analysis manager

    // Get the pointer to the existing analysis manager via the static access method.
    AliAnalysisManager *analysisManager = AliAnalysisManager::GetAnalysisManager();
    if( !analysisManager )
    {
        ::Error("AddTaskEMCalNeutrals", "No analysis manager to connect to.");
        return NULL;
    }

    // get input event handler type and check
    TString type = analysisManager->GetInputEventHandler()->GetDataType();
    if (!type.Contains( "AOD" ) )
    {
        Error("AddTaskEMCalNeutrals", "AOD input handler needed!");
        return NULL;
    }

    // create task
    AliEMCalpi0Task* clusterEvaluation = new AliEMCalpi0Task( "clusterEvaluation" );
    clusterEvaluation->SelectCollisionCandidates(AliVEvent::kSemiCentral| AliVEvent::kCentral|AliVEvent::kMB );

    // add task to manager
    analysisManager->AddTask(clusterEvaluation);

    //
    analysisManager->ConnectInput( clusterEvaluation, 0, analysisManager->GetCommonInputContainer() );
    analysisManager->ConnectOutput( clusterEvaluation, 0, analysisManager->GetCommonOutputContainer() );

    return clusterEvaluation;

}
