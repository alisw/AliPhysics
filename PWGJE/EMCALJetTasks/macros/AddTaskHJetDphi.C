AliAnalysisTaskHJetDphi *AddTaskHJetDphi(const char *name = "HJetDphi", 
					 const TString period = "lhc11h",
					 const UInt_t pSel = AliVEvent::kAnyINT | AliVEvent::kCentral | AliVEvent::kSemiCentral,
					 const Bool_t isEmbed = kFALSE,
					 const Bool_t ismc = kFALSE, 
					 const Bool_t anaTruth = kFALSE,
					 const char *nonStdFile = "",
					 const char *mcParticleArrayName = "",
					 const char *trkArrayName = "",
					 const char *trkHybridName = "",
					 const char *jetArrayName = "",
					 const char *jetDLArrayName = "",
					 const char *jetPLArrayName = "",
					 const char *rhoName = "",
					 const Double_t radius = 0.4,
					 const UInt_t filterMask = 768,
					 const Bool_t requireITSrefit = kTRUE,
					 const Double_t minTT = 20,
					 const Double_t maxTT = 50,
					 const Double_t ptTrk = 0.15,
					 const Double_t ptJet = 10,
					 const Bool_t runSingleInclHJet = kTRUE,
					 const Bool_t runBkgFlow = kTRUE,
					 const Bool_t runTrkQA = kTRUE,
					 const Bool_t runJetQA = kTRUE,
					 const Bool_t runPLHJet = kFALSE,
					 const Bool_t runDLHJet = kFALSE,
					 const Bool_t runLeadTrkQA = kFALSE,
					 const Int_t aodTrkBit0 = 256,
					 const Int_t aodTrkBit1 = 512)
{
   Printf("Adding h+jet dphi task\n");

   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if(!mgr){
      ::Error("AddTaskHJetDphi", "No analysis manager to connect to.");
      return NULL;
   }
   if(!mgr->GetInputEventHandler()){
      ::Error("AddTaskHJetDphi", "This task requires an input event handler.");
      return NULL;
   }

  if(period.Contains("lhc11a",TString::kIgnoreCase))
    pSel = AliVEvent::kMB | AliVEvent::kEMC1;
  if(period.Contains("lhc12a15e",TString::kIgnoreCase))
    pSel = 0;

  char *tname = Form("%s_Charged_ANTIKT%02d_Filter%d_Cut%05d",name,(Int_t)(radius*10),filterMask,(Int_t)(ptTrk*1e3));
  AliAnalysisTaskHJetDphi *hJetTask = new AliAnalysisTaskHJetDphi(tname);
  hJetTask->SetIsEmbedding(isEmbed);
  hJetTask->SetAnaType(1);
  hJetTask->SetRunPeriod(period.Data());
  if(period.Contains("lhc10h",TString::kIgnoreCase) || period.Contains("lhc11h",TString::kIgnoreCase))
    hJetTask->SetCollisionSystem("PbPb");
  else
    hJetTask->SetCollisionSystem("pp");
  hJetTask->SetIsMC(ismc);
  hJetTask->SetAnalyzeMCTruth(anaTruth);
  hJetTask->SetOfflineTrgMask(pSel);
  hJetTask->SetMaxVtxZ(10);
  hJetTask->SetFilterMask(filterMask);
  hJetTask->SetRequireITSRefit(requireITSrefit);
  hJetTask->SetNonStdFile(nonStdFile);
  hJetTask->SetMcParticleArrName(mcParticleArrayName);
  hJetTask->SetEmbTrkArrName(trkHybridName);
  hJetTask->SetTrackArrName(trkArrayName);
  hJetTask->SetJetArrName(jetArrayName);
  hJetTask->SetPLJetArrName(jetPLArrayName);
  hJetTask->SetDLJetArrName(jetDLArrayName);
  hJetTask->SetRhoName(rhoName);
  hJetTask->SetRadius(radius);
  hJetTask->SetTrkPtRange(ptTrk,1e4);
  hJetTask->SetTrkPhiRange(0,2*TMath::Pi());
  hJetTask->SetTrkEtaRange(-0.9,0.9);
  hJetTask->SetTTRange(minTT,maxTT);
  hJetTask->SetJetPtMin(ptJet);
  hJetTask->SetRunSingleInclHJet(runSingleInclHJet);
  hJetTask->SetRunTrkQA(runTrkQA);
  hJetTask->SetRunJetQA(runJetQA);
  hJetTask->SetRunPLHJet(runPLHJet);
  hJetTask->SetRunDLHJet(runDLHJet);
  hJetTask->SetRunLeadTrkQA(runLeadTrkQA);
  hJetTask->SetRunBkgFlow(runBkgFlow);
  hJetTask->SetAODfilterBits(aodTrkBit0,aodTrkBit1);

  mgr->AddTask(hJetTask);
  TString foutputrmaAliAnalysisTaskHJetDphi = Form("rma_AliHadJetQATask.root");
  AliAnalysisDataContainer *coutput = mgr->CreateContainer(Form("%s_TT%1.0f%1.0f",tname,minTT,maxTT),
							   TList::Class(),
							   AliAnalysisManager::kOutputContainer,Form("%s:%s_TT%1.0f%1.0f",AliAnalysisManager::GetCommonFileName(),tname,minTT,maxTT));
  mgr->ConnectInput (hJetTask, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(hJetTask, 0, mgr->GetCommonOutputContainer());
  mgr->ConnectOutput(hJetTask, 1, coutput);
  
  return hJetTask;
}
