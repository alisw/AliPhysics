AliAnalysisTaskSE* AddTaskZDCQA()
{
  // Creates a QA task to check ZDC data

  AliCDBEntry* entry = AliCDBManager::Instance()->Get("GRP/GRP/Data");
  if(!entry) {
    ::Error("AddTaskZDCQA", "Cannot get AliCDBEntry for GRP/GRP/Data");
    return NULL;
  }
  const AliGRPObject* grpData = dynamic_cast<AliGRPObject*>(entry->GetObject());
  if(!grpData) {
    ::Error("AddTaskZDCQA", "Cannot get AliGRPObject");
    return NULL;
  }
  TString beamType = grpData->GetBeamType();
  if(beamType==AliGRPObject::GetInvalidString()) {
    ::Error("AddTaskZDCQA", "Missing beam type!!!!");
    return NULL;
  }

  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskZDCQA", "No analysis manager to connect to.");
    return NULL;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskZDCQA", "This task requires an input event handler");
    return NULL;
  }
   TString inputDataType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"

   // Configure analysis
   //===========================================================================
   AliAnalysisTaskSE *task;
   if(((beamType.CompareTo("pp"))==0) || ((beamType.CompareTo("p-p"))==0)
    ||((beamType.CompareTo("PP"))==0) || ((beamType.CompareTo("P-P"))==0)){
     gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/ZDC/AddTaskZDCQApp.C");
     task = reinterpret_cast<AliAnalysisTaskSE*>(gInterpreter->ProcessLine("AddTaskZDCQApp();"));
   }
   else if(((beamType.CompareTo("p-A"))==0) || ((beamType.CompareTo("A-p"))==0)
     ||((beamType.CompareTo("P-A"))==0) || ((beamType.CompareTo("A-P"))==0)){
     gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/ZDC/AddTaskZDCQApA.C");
     task = reinterpret_cast<AliAnalysisTaskSE*>(gInterpreter->ProcessLine("AddTaskZDCQApA();"));
   }
   else if((beamType.CompareTo("A-A")) == 0 || (beamType.CompareTo("AA")) == 0){
     gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/ZDC/AddTaskZDCQAPbPb.C");
     task = reinterpret_cast<AliAnalysisTaskSE*>(gInterpreter->ProcessLine("AddTaskZDCQAPbPb();"));
   }

   return task;
}
