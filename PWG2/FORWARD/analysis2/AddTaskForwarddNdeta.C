/**
 * @file   AddTaskForwarddNdeta.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Fri Jan 28 10:22:26 2011
 * 
 * @brief Script to add a multiplicity task for the central
 *        @f$\eta@f$ region
 * 
 * 
 */
AliAnalysisTask*
AddTaskForwarddNdeta(const char* trig="INEL", Double_t vzMin=-10, Double_t vzMax=10, Float_t centlow = 0, Float_t centhigh = 100, Bool_t cutEdges = false)
{
  // analysis manager
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  
  // Make our object.  2nd argumenent is absolute max Eta 
  // 3rd argument is absolute max Vz
  AliForwarddNdetaTask* task = new AliForwarddNdetaTask("Forward");
  task->SetVertexRange(vzMin, vzMax);
  task->SetTriggerMask(trig);
  task->SetCutEdges(cutEdges);
  task->SetCentLow(centlow);
  task->SetCentHigh(centhigh);
   
  //task->SetTriggerEff(0.93);
  /*if(trig == "INEL") task->SetTriggerEff(0.95);
  if(trig == "NSD")  task->SetTriggerEff(1.04);
  TFile f("/home/canute/ALICE/FMDanalysis/productionData/normalizationHists900GeV.root", "READ");
  //TFile f("/home/canute/ALICE/FMDanalysis/BackgroundCorrection/normalizationHists.root", "READ");
  TH2D* hnorm = 0 ;
  if(trig == "INEL") hnorm = (TH2D*)f.Get("hInelNormalization");
  if(trig == "NSD") hnorm = (TH2D*)f.Get("hNSDNormalization");
  task->SetShapeCorrection(hnorm);*/
  mgr->AddTask(task);

  // create containers for input/output
  AliAnalysisDataContainer *sums = 
    mgr->CreateContainer("ForwardSums", TList::Class(), 
			 AliAnalysisManager::kOutputContainer, 
			 AliAnalysisManager::GetCommonFileName());
  AliAnalysisDataContainer *output = 
    mgr->CreateContainer("ForwardResults", TList::Class(), 
			 AliAnalysisManager::kParamContainer, 
			 AliAnalysisManager::GetCommonFileName());
  
  // connect input/output
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, sums);
  mgr->ConnectOutput(task, 2, output);

  return task;
}

  
//________________________________________________________________________
//
// EOF
// 
