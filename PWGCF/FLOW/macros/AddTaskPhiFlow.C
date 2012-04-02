/////////////////////////////////////////////////////////////////////////////////////////////
//
// AddTaskPhiFlow macro
// Author: Redmer A. Bertens, Utrecht University, 2012
/////////////////////////////////////////////////////////////////////////////////////////////
AliAnalysisTaskPhiFlow* AddTaskPhiFlow(Bool_t SP = kTRUE,
                                       Bool_t QC = kTRUE,
                                       Bool_t EP = kTRUE,
                                       Int_t harmonic = 2,
                                       Float_t centrMin = 30.,
                                       Float_t centrMax = 40.,
                                       Float_t bayesThresholdP = 0.5,
                                       TString suffixName = "",
                                       Bool_t bCentralTrigger = kFALSE,
                                       Float_t RPEtaMinA = -0.8,
                                       Float_t RPEtaMaxA = 0.0,
                                       Float_t RPEtaMinB = 0.0,
                                       Float_t RPEtaMaxB = 0.8,
                                       Float_t POIEtaMin = -0.8,
                                       Float_t POIEtaMax = 0.8,
                                       Float_t POIPtMin = 0.15,
                                       Float_t POIPtMax = 10.,
                                       Float_t deltaDip = 0.04,
                                       Float_t deltaDipMaxPt = 1.2,
                                       AliAnalysisTaskPhiFlow::PIDtype pidtype = AliAnalysisTaskPhiFlow::kCombined,
                                       Bool_t strictKaonCuts = kFALSE,
                                       Bool_t TPCStandAloneTracks = kFALSE,
                                       Float_t vertexZ = 10.)
{
   // set up main output container's name
   TString centralityName("");
   centralityName += Form("%.0f", centrMin);
   centralityName += "-";
   centralityName += Form("%.0f", centrMax);
   centralityName += "_";
   centralityName += Form("vZ%.f", vertexZ);
   centralityName += "_";
   centralityName += Form("bayP%.1f", bayesThresholdP);
   centralityName += "_";
   centralityName += "_POIEta";
   centralityName += Form("%.1f", POIEtaMin);
   centralityName += "-";
   centralityName += Form("%.1f", POIEtaMax);
   centralityName += "_RPEta";
   centralityName += Form("%.1f", RPEtaMinA);
   centralityName += "-";
   centralityName += Form("%.1f", RPEtaMaxA);
   centralityName += "-";
   centralityName += Form("%.1f", RPEtaMinB);
   centralityName += "-";
   centralityName += Form("%.1f", RPEtaMaxB);
   centralityName += "-";
   centralityName += Form("dDip%.2f", deltaDip);
   centralityName += "-";
   centralityName += Form("dDipPt%.2f", deltaDipMaxPt);
   if (strictKaonCuts)
   {
      centralityName += "-";
      centralityName += "StrictKaonCuts";
   }
   if (bCentralTrigger)
   {
       centralityName += "-";
       centralityName += "kMBkCkSC";
   }

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

    // set triggers
   if(bCentralTrigger) task->SelectCollisionCandidates(AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral);
   else                task->SelectCollisionCandidates(AliVEvent::kMB);

   // set RP cuts for EP method
   AliFlowTrackCuts* cutsRP = new AliFlowTrackCuts("GlobalRP");
   AliFlowTrackCuts* cutsRP = cutsRP->GetStandardTPCStandaloneTrackCuts();

   // set POI cuts for kaon selection
   AliFlowTrackCuts* cutsPOI = new AliFlowTrackCuts("GlobalPOI");
   AliFlowTrackCuts* cutsPOI = cutsPOI->GetStandardGlobalTrackCuts2010();

   // set some more cuts ...
   task->EventPlanePtCut(kFALSE); // set to true to enforce pt cut on event plane tracks
   task->SetMaxDeltaDipAngleAndPt(deltaDip, deltaDipMaxPt); // set an upper value for d-dip angle and pt of tracks on which this should be applied
   task->SetCentralityParameters(centrMin, centrMax, "V0M");
   task->SetVertexZ(vertexZ);
   task->SetRPCuts(cutsRP);
   task->SetPOICuts(cutsPOI);

   // set the kaon cuts, and specify the PID procedure which will be used
   task->SetIdentificationType(pidtype);
   task->SetBayesianProbability(bayesThresholdP);
   if (strictKaonCuts) task->SetRequireStrictKaonCuts();
   if (TPCStandAloneTracks) task->SetRequireTPCStandAloneKaons();
   task->SetCandidateEtaAndPt(POIEtaMin, POIEtaMax, POIPtMin, POIPtMax);

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
   task->SetEtaRanges(RPEtaMinA, RPEtaMaxA, RPEtaMinB, RPEtaMaxB);
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

   // set, if necessary, a suffixname as a unique identifier for each wagon
   if(suffixName == "") {
       suffixName += Form("%.0f", centrMin);
       suffixName += Form("%.0f", centrMax);
   }

   if (SP) // set up scalar product (SP) tasks
   {
      AliAnalysisTaskScalarProduct *taskSP[30];
      AliAnalysisDataContainer *coutputSP[30];
      for (int r = 0; r != 30; ++r)
      {
         taskSP[r] = new AliAnalysisTaskScalarProduct(Form("SP_MassBand_%d_%s", r, suffixName.Data()), kFALSE);
         taskSP[r]->SetHarmonic(harmonic);
         taskSP[r]->SetRelDiffMsub(1.0);
         taskSP[r]->SetApplyCorrectionForNUA(kTRUE);
         mgr->AddTask(taskSP[r]);
         coutputSP[r] = mgr->CreateContainer(Form("cobjSP_MassBand_%d_%s", r, suffixName.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, fileName);
         mgr->ConnectInput(taskSP[r], 0, coutputCandidates[r]);
         mgr->ConnectOutput(taskSP[r], 1, coutputSP[r]);
      }
   }

   if (EP) // set up event plane (ep) tasks
   {
      AliAnalysisTaskScalarProduct *taskEP[30];
      AliAnalysisDataContainer *coutputEP[30];
      for (int r = 0; r != 30; ++r)
      {
         taskEP[r] = new AliAnalysisTaskScalarProduct(Form("EP_MassBand_%d_%s", r, suffixName.Data()), kFALSE);
         taskEP[r]->SetBehaveAsEP();
         taskEP[r]->SetHarmonic(harmonic);
         taskEP[r]->SetRelDiffMsub(1.0);
         taskEP[r]->SetApplyCorrectionForNUA(kTRUE);
         mgr->AddTask(taskEP[r]);
         coutputEP[r] = mgr->CreateContainer(Form("cobjEP_MassBand_%d_%s", r, suffixName.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, fileName);
         mgr->ConnectInput(taskEP[r], 0, coutputCandidates[r]);
         mgr->ConnectOutput(taskEP[r], 1, coutputEP[r]);
      }
   }

   if (QC) // set up q-cumulants task
   {
      AliAnalysisTaskQCumulants *taskQC[30];
      AliAnalysisDataContainer *coutputQC[30];
      for (int r = 0; r != 30; ++r)
      {
         taskQC[r] = new AliAnalysisTaskQCumulants(Form("QC_MassBand_%d_%s", r, suffixName.Data()), kFALSE);
         taskQC[r]->SetHarmonic(harmonic);
         taskQC[r]->SetCalculateCumulantsVsM(kFALSE);
         taskQC[r]->SetnBinsMult(10000);
         taskQC[r]->SetMinMult(0.);
         taskQC[r]->SetMaxMult(10000.);
         taskQC[r]->SetApplyCorrectionForNUA(kTRUE);
         taskQC[r]->SetFillMultipleControlHistograms(kFALSE);
         mgr->AddTask(taskQC[r]);
         coutputQC[r] = mgr->CreateContainer(Form("cobjQC_MassBand_%d_%s", r, suffixName.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, fileName);
         mgr->ConnectInput(taskQC[r], 0, coutputCandidates[r]);
         mgr->ConnectOutput(taskQC[r], 1, coutputQC[r]);
      }
   }
   mgr->ConnectInput(task, 0, cinput);
   mgr->ConnectOutput(task, 1, coutput);

   // return the task
   return task;
};

