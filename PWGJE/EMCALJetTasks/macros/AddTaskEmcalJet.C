// $Id$ 

AliEmcalJetTask* AddTaskEmcalJet(
  const char *nTracks        = "Tracks",
  const char *nClusters      = "CaloClusters",
  const Int_t algo           = 1,
  const Double_t radius      = 0.4,
  const Int_t type           = 0,
  const Double_t minTrPt     = 0.15,
  const Double_t minClPt     = 0.15
)
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskAliEmcalJet", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskAliEmcalJet", "This task requires an input event handler");
    return NULL;
  }
  
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  Double_t minJetPt = 0;

  char *algoString;
  if (algo == 0) {
    algoString = "KT";
    minJetPt = 0.1;
  } else if (algo == 1) {
    algoString = "AKT";
    minJetPt = 1.0;
  }

  char *typeString;
  if (type == 0)
    typeString = "Full";
  else if (type == 1)
    typeString = "Charged";
  else if (type == 2)
    typeString = "Neutral";

  char radiusString[200];
  sprintf(radiusString,"R0%2.0f",radius*100.0);

  char pTString[200];
  if (minTrPt==0)
    sprintf(pTString,"pT0000");
  else if (minTrPt<1.0)
    sprintf(pTString,"pT0%3.0f",minTrPt*1000.0);
  else if (minTrPt>=1.0)
    sprintf(pTString,"pT%4.0f",minTrPt*1000.0);

  char ETString[200];
  if (minClPt==0)
    sprintf(ETString,"ET0000");
  else if (minClPt<1.0)
    sprintf(ETString,"ET0%3.0f",minClPt*1000.0);
  else if (minClPt>=1.0)
    sprintf(ETString,"ET%4.0f",minClPt*1000.0);  

  TString name;
  if (type == 0)
    name = TString(Form("Jet_%s%s%s_%s_%s_%s_%s",
                        algoString,typeString,radiusString,nTracks,pTString,nClusters,ETString));
  else if (type == 1)
    name = TString(Form("Jet_%s%s%s_%s_%s",
                        algoString,typeString,radiusString,nTracks,pTString));
  else if (type == 2)
    name = TString(Form("Jet_%s%s%s_%s_%s",
                        algoString,typeString,radiusString,nClusters,ETString));
 
  AliEmcalJetTask* mgrTask = mgr->GetTask(name.Data());
  if (mgrTask)
    return mgrTask;  

  AliEmcalJetTask* jetTask = new AliEmcalJetTask(name);
  jetTask->SetTracksName(nTracks);
  jetTask->SetClusName(nClusters);
  jetTask->SetJetsName(name);
  jetTask->SetAlgo(algo);
  jetTask->SetMinJetTrackPt(minTrPt);
  jetTask->SetMinJetClusPt(minClPt);
  jetTask->SetMinJetPt(minJetPt);
  jetTask->SetRadius(radius);
  jetTask->SetType(type);

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------

  mgr->AddTask(jetTask);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer()  ;
  mgr->ConnectInput  (jetTask, 0, cinput);

  return jetTask;
}
