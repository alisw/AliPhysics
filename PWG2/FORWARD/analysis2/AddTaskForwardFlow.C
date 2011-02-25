void AddTaskForwardFlow(TString type = "", Int_t etabins = 20)
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddFMDFlowTask", "No analysis manager to connect to.");
    return NULL;
  }   

  AliAODInputHandler* aodInput = dynamic_cast<AliAODInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
   
  Bool_t aod = kFALSE;
  if (aodInput) aod = kTRUE;
  if (!aod) {
    Error("AddTaskForwardFlow", "No analysis manager to connect to.");
    return NULL;
  }

  // --- Check which harmonics to calculate --- //

  Bool_t v1 = kTRUE;
  Bool_t v2 = kTRUE;
  Bool_t v3 = kTRUE;
  Bool_t v4 = kTRUE;

  if (type.Length() > 0) {
    if (!type.Contains("1")) v1 = kFALSE;
    if (!type.Contains("2")) v2 = kFALSE;
    if (!type.Contains("3")) v3 = kFALSE;
    if (!type.Contains("4")) v4 = kFALSE;
  }

  // --- Create output containers and find input from fmd task --- //

  TString outputFile = AliAnalysisManager::GetCommonFileName();
  outputFile += ":FlowResults";

  AliAnalysisDataContainer* qcout = mgr->CreateContainer("QCumulants", TList::Class(), AliAnalysisManager::kOutputContainer, outputFile);

  // --- For the selected flow tasks the input and output is set --- //
  
  AliForwardFlowTaskQC* qc = new AliForwardFlowTaskQC("QCumulants");

  qc->SetDoHarmonics(v1, v2, v3, v4);
  qc->SetUseNEtaBins(etabins);
  
  mgr->ConnectInput(qc ,0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(qc, 0, mgr->GetCommonOutputContainer());
  mgr->ConnectOutput(qc, 1, qcout);

  return;
}
