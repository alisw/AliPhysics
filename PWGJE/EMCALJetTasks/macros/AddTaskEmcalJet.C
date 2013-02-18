// $Id$ 

AliEmcalJetTask* AddTaskEmcalJet(
  const UInt_t type          = AliEmcalJetTask::kAKT | AliEmcalJetTask::kFullJet | AliEmcalJetTask::kR040Jet,
  const char *nTracks        = "Tracks",
  const char *nClusters      = "CaloClusters",
  const Double_t minTrPt     = 0.15,
  const Double_t minClPt     = 0.30,
  const Double_t ghostArea   = 0.01,
  const Double_t radius      = 0.4,
  const char *tag            = "Jet"
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
  if ((type & AliEmcalJetTask::kKT) != 0) {
    algoString = "KT";
    minJetPt = 0.1;
  } else if ((type & AliEmcalJetTask::kAKT) != 0) {
    algoString = "AKT";
    minJetPt = 1.0;
  }

  char *typeString;
  if ((type & AliEmcalJetTask::kFullJet) != 0)
    typeString = "Full";
  else if ((type & AliEmcalJetTask::kChargedJet) != 0)
    typeString = "Charged";
  else if ((type & AliEmcalJetTask::kNeutralJet) != 0)
    typeString = "Neutral";

  char radiusString[200];
  if ((type & AliEmcalJetTask::kR020Jet) != 0)
    sprintf(radiusString,"R020");
  else if ((type & AliEmcalJetTask::kR030Jet) != 0)
    sprintf(radiusString,"R030");
  else if ((type & AliEmcalJetTask::kR040Jet) != 0)
    sprintf(radiusString,"R040");
  else
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
  if (*nTracks && *nClusters)
    name = TString(Form("%s_%s%s%s_%s_%s_%s_%s",
                        tag,algoString,typeString,radiusString,nTracks,pTString,nClusters,ETString));
  else if (!*nClusters)
    name = TString(Form("%s_%s%s%s_%s_%s",
                        tag,algoString,typeString,radiusString,nTracks,pTString));
  else if (!*nTracks)
    name = TString(Form("%s_%s%s%s_%s_%s",
                        tag,algoString,typeString,radiusString,nClusters,ETString));
 
  AliEmcalJetTask* mgrTask = mgr->GetTask(name.Data());
  if (mgrTask)
    return mgrTask;  

  AliEmcalJetTask* jetTask = new AliEmcalJetTask(name);
  jetTask->SetTracksName(nTracks);
  jetTask->SetClusName(nClusters);
  jetTask->SetJetsName(name);
  jetTask->SetJetType(type);
  jetTask->SetMinJetTrackPt(minTrPt);
  jetTask->SetMinJetClusPt(minClPt);
  jetTask->SetMinJetPt(minJetPt);
  if ((type & (AliEmcalJetTask::kRX1Jet|AliEmcalJetTask::kRX2Jet|AliEmcalJetTask::kRX3Jet)) != 0)
    jetTask->SetRadius(radius);
  jetTask->SetGhostArea(ghostArea);

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------

  mgr->AddTask(jetTask);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer()  ;
  mgr->ConnectInput  (jetTask, 0, cinput);

  return jetTask;
}


AliEmcalJetTask* AddTaskEmcalJet(
  const char *nTracks        = "Tracks",
  const char *nClusters      = "CaloClusters",
  const Int_t algo           = 1,
  const Double_t radius      = 0.4,
  const Int_t type           = 0,
  const Double_t minTrPt     = 0.15,
  const Double_t minClPt     = 0.30,
  const Double_t ghostArea   = 0.01  ,
  const char *tag            = "Jet"
)
{  
  
  UInt_t jetType = 0;

  if (algo == 0) 
    jetType |= AliEmcalJetTask::kKT; 
  else 
    jetType |= AliEmcalJetTask::kAKT;

  if (type==0)
    jetType |= AliEmcalJetTask::kFullJet; 
  else if (type==1) 
    jetType |= AliEmcalJetTask::kChargedJet; 
  else if (type==2) 
    jetType |= AliEmcalJetTask::kNeutralJet;

  if (radius==0.2) 
    jetType |= AliEmcalJetTask::kR020Jet; 
  else if (radius==0.3) 
    jetType |= AliEmcalJetTask::kR030Jet;
  else if (radius==0.4) 
    jetType |= AliEmcalJetTask::kR040Jet;
  else
    jetType |= AliEmcalJetTask::kRX1Jet;

  return AddTaskEmcalJet(jetType, nTracks, nClusters, minTrPt, minClPt, ghostArea, radius, tag);
}
