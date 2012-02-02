/////////////////////////////////////////////////////////////////////////////////////////////
//
// AddTaskPhiFlow macro
// Author: Redmer A. Bertens, Utrecht University, 2012
/////////////////////////////////////////////////////////////////////////////////////////////
AliAnalysisTaskPhiFlow* AddTaskPhiFlow( Float_t centrMin = 30.,
                                        Float_t centrMax = 40.,
                                        Float_t vertexZ = 10.,
                                        Float_t bayesThresholdP = 0.5,
                                        Float_t minEtaA = -0.8,
                                        Float_t maxEtaA = 0.,
                                        Float_t minEtaB = 0.,
                                        Float_t maxEtaB = 0.8,
                                        Float_t deltaDip = 0.0,
                                        Float_t deltaDipMaxPt = 0.0,
                                        Float_t eventPlanePtCut = 10,
                                        AliAnalysisTaskPhiFlow::PIDtype pidtype = AliAnalysisTaskPhiFlow::kCombined)
{
   TString centralityName("");
   centralityName += Form("%.0f", centrMin);
   centralityName += "-";
   centralityName += Form("%.0f", centrMax);
   centralityName += "_";
   centralityName += Form("vZ%.1f", vertexZ);
   centralityName += "_";
   centralityName += Form("bayP%.1f", bayesThresholdP);
   centralityName += "_";
   centralityName += "_EtaA";
   centralityName += Form("%.1f", minEtaA);
   centralityName += "-";
   centralityName += Form("%.1f", maxEtaA);
   centralityName += "_EtaB";
   centralityName += Form("%.1f", minEtaB);
   centralityName += "-";
   centralityName += Form("%.1f", maxEtaB);
   centralityName += "-";
   centralityName += Form("dDip%.f", deltaDip);
   centralityName += "-";
   centralityName += Form("dDipPt%.f", deltaDipMaxPt);
   centralityName += "-";
   centralityName += Form("epPtCut%.f", eventPlanePtCut);

   TString fileName = AliAnalysisManager::GetCommonFileName();
   fileName += ":PhiV2FlowEvents";

   // get the manager
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr)
   {
      ::Error("AddTaskPhiV2", "No analysis manager to connect to.");
      return NULL;
   }
   if (!mgr->GetInputEventHandler())
   {
      ::Error("AddTaskPhiV2", "This task requires an input event handler");
      return NULL;
   }

   // create the main task
   AliAnalysisTaskPhiFlow *task = new AliAnalysisTaskPhiFlow("TaskPhiFlow");
   task->SelectCollisionCandidates(AliVEvent::kMB);

   // set RP cuts for EP method
   AliFlowTrackCuts* cutsRP = new AliFlowTrackCuts("GlobalRP");
   cutsRP->GetStandardTPCStandaloneTrackCuts();

   // set POI cuts for EP method
   AliFlowTrackCuts* cutsPOI = new AliFlowTrackCuts("GlobalPOI");
   cutsPOI->GetStandardGlobalTrackCuts2010();

   // set and configure kaon track cuts (for the bayesian PID method)
   // select this method by setting PIDtype to kCombined
   AliFlowTrackCuts* cutsKaon = new AliFlowTrackCuts("kaon_cuts");
   AliFlowTrackCuts::PIDsource pid = AliFlowTrackCuts::kTPCbayesian;
   cutsKaon->SetPID(AliPID::kKaon, pid, bayesThresholdP);
   cutsKaon->SetAllowTOFmismatchFlag(kTRUE);
   cutsKaon->SetRequireStrictTOFTPCagreement(kFALSE);

   // set some more cuts ...
   task->EventPlanePtCut(kTRUE); // set to true to enforce pt cut on event plane tracks
   task->SetEventPlanePtCut(eventPlanePtCut); // set the upper value for Pt
   task->SetMaxDeltaDipAngleAndPt(deltaDip, deltaDipMaxPt); // set an upper value for d-dip angle and pt of tracks on which this should be applied
   task->SetCentralityParameters(centrMin, centrMax, "V0M");
   task->SetVertexZ(vertexZ);
   task->SetRPCuts(cutsRP);
   task->SetPOICuts(cutsPOI);

   // set the kaon cuts, and specify the PID procedure which will be used
   task->SetKaonCuts(cutsKaon);
   task->SetIdentificationType(pidtype);

   //set RP cuts for flow package analysis
   cutsFlowRP = new AliFlowTrackCuts("rp_cuts_111");
   AliFlowTrackCuts::trackParameterType rptype = AliFlowTrackCuts::kGlobal;
   cutsFlowRP->SetParamType(rptype);
   cutsFlowRP->SetPtRange(0.2, 5.0);
   cutsFlowRP->SetEtaRange(-0.8, 0.8);
   cutsFlowRP->SetMinNClustersTPC(70);
   cutsFlowRP->SetMinChi2PerClusterTPC(0.1);
   cutsFlowRP->SetMaxChi2PerClusterTPC(4.0);
   cutsFlowRP->SetRequireTPCRefit(kTRUE);
   cutsFlowRP->SetMaxDCAToVertexXY(0.3);
   cutsFlowRP->SetMaxDCAToVertexZ(0.3);
   cutsFlowRP->SetAcceptKinkDaughters(kFALSE);
   cutsFlowRP->SetMinimalTPCdedx(10.);

   Double_t flowBands[2][30];
   for (Int_t i = 0; i < 30; i++)
   {
      flowBands[0][i] = 0.99 + i * 0.0034;
      flowBands[1][i] = 0.99 + (i + 1) * 0.0034;
   }
   task->SetMassRanges(flowBands);
   task->SetEtaRanges(minEtaA, maxEtaA, minEtaB, maxEtaB);
   task->SetFlowRPCuts(cutsFlowRP);
   mgr->AddTask(task);

   AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
   AliAnalysisDataContainer *coutput = mgr->CreateContainer(Form("PhiV2_%s", centralityName.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, fileName);

   AliAnalysisDataContainer *coutputCandidates[30];

   for (int r = 0; r != 30; ++r)
   {
      coutputCandidates[r] = mgr->CreateContainer(
                                Form("Flow_MassBand_%d", r),
                                AliFlowEventSimple::Class(), AliAnalysisManager::kExchangeContainer);
      mgr->ConnectOutput(task, 2 + r, coutputCandidates[r]);
   }

   // Scalar Product
   AliAnalysisTaskScalarProduct *taskSP[30];
   AliAnalysisDataContainer *coutputSP[30];
   for (int r = 0; r != 30; ++r)
   {
      taskSP[r] = new AliAnalysisTaskScalarProduct(Form("SP_MassBand_%d", r), kFALSE);
      taskSP[r]->SetRelDiffMsub(1.0);
      taskSP[r]->SetApplyCorrectionForNUA(kTRUE);
      mgr->AddTask(taskSP[r]);
      coutputSP[r] = mgr->CreateContainer(Form("cobjSP_MassBand_%d", r), TList::Class(), AliAnalysisManager::kOutputContainer, fileName);
      mgr->ConnectInput(taskSP[r], 0, coutputCandidates[r]);
      mgr->ConnectOutput(taskSP[r], 1, coutputSP[r]);
   }
   // Q-Cumulants
   AliAnalysisTaskQCumulants *taskQC[30];
   AliAnalysisDataContainer *coutputQC[30];
   for (int r = 0; r != 30; ++r)
   {
      taskQC[r] = new AliAnalysisTaskQCumulants(Form("QC_MassBand_%d", r), kFALSE);
      taskQC[r]->SetCalculateCumulantsVsM(kFALSE);
      taskQC[r]->SetnBinsMult(10000);
      taskQC[r]->SetMinMult(0.);
      taskQC[r]->SetMaxMult(10000.);
      taskQC[r]->SetApplyCorrectionForNUA(kTRUE);
      taskQC[r]->SetFillMultipleControlHistograms(kFALSE);
      mgr->AddTask(taskQC[r]);
      coutputQC[r] = mgr->CreateContainer(Form("cobjQC_MassBand_%d", r), TList::Class(), AliAnalysisManager::kOutputContainer, fileName);
      mgr->ConnectInput(taskQC[r], 0, coutputCandidates[r]);
      mgr->ConnectOutput(taskQC[r], 1, coutputQC[r]);
   }

   mgr->ConnectInput(task, 0, cinput);
   mgr->ConnectOutput(task, 1, coutput);

   // return the task
   return task;
};

