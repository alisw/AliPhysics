/////////////////////////////////////////////////////////////////////////////////////////////
//
// AddTaskPhiFlow macro
// Author: Redmer A. Bertens, Utrecht University, 2012
//         Many thanks to Carlos Perez Lara
//         Commented where necessary
/////////////////////////////////////////////////////////////////////////////////////////////

class AliAnalysisDataContainer;
class AliFlowTrackCuts;
class AliFlowTrackSimpleCuts;
class AliFlowEventCuts;
class AliFlowEventSimpleCuts;
class AliAnalysisDataContainer;

AliAnalysisTaskPhiFlow* AddTaskPhiFlow(Bool_t SP = kTRUE,
                                       Bool_t SPSUB = kTRUE,
                                       Bool_t QC = kTRUE,
                                       Bool_t EP = kTRUE,
                                       Float_t centrMin = 20.,
                                       Float_t centrMax = 30.,
                                       Double_t ITSsigma = 0.,
                                       Double_t ITSrange = 0.,
                                       Double_t TPCcontrol = 1.,
                                       Double_t TPCsigma = 3.,
                                       Double_t TPCrange = 0.,
                                       Double_t ITScontrol = -1.,
                                       Double_t Bpurity = 0.5,
                                       TString suffixName = "UniqueID",
                                       Bool_t bCentralTrigger = kFALSE,
                                       Float_t EtaGap = 0.,
                                       Float_t POIEtaMin = -0.8,
                                       Float_t POIEtaMax = 0.8,
                                       Float_t POIPtMin = 0.15,
                                       Float_t POIPtMax = 10.,
                                       Float_t deltaDip = 0.,
                                       Float_t deltaDipMaxPt = 0.,
                                       TString DCA = "pt",
                                       Int_t harm = 2,
                                       Bool_t TPCStandAloneTracks = kFALSE,
                                       Bool_t useGlobalRPCuts = kTRUE,
                                       Float_t vertexZ = 10.,
                                       Bool_t shrinkSP = kTRUE,
                                       Bool_t debug = kFALSE)
{
   Double_t PIDconfig[] = {ITSsigma, ITSrange, TPCcontrol, TPCsigma, TPCrange, ITScontrol, Bpurity};
   // main function, create and add tasks
   if(debug) cout << " === Adding Task PhiFlow === " << endl;
   // set up main output container's name
   TString fileName = AliAnalysisManager::GetCommonFileName();
   fileName += ":PhiReconstruction";
   suffixName += Form("%.0f", centrMin);
   suffixName += Form("%.0f", centrMax);
   fileName+=suffixName;
   if(debug) cout << "    --> Reconstruction data container: " << fileName << endl;
   // check validity of arguments
   if((!SPSUB)&&(EtaGap > 0.)) {
       if(debug) cout << " --> Fatal error: Eta gap is introduced but method SPSUB is not enabled !!! <-- " << endl;
       return 0x0;
   }
   else if ((QC||SP||EP)&&(EtaGap > 0.)) {
       if(debug) cout << " --> Fatal error: one or more of the flow analyses is not compatible with the Eta Gap !!! <--" << endl;
       return 0x0;
   }
   // get the manager and input event handler
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      if(debug) cout << " Fatal error: no analysis manager found! " << endl;
      return 0x0;
   }
   if (!mgr->GetInputEventHandler()) {
      if(debug) cout << " Fatal error: no imput event handler found!" << endl;
      return 0x0;
   }
   // create the main task3
   AliAnalysisTaskPhiFlow *task = new AliAnalysisTaskPhiFlow("TaskPhiFlow");
   if(debug) cout << " === AliAnalysisTaskPhiFlow === " << task << endl;
   if(!task) {
       if(debug) cout << " --> Unexpected error occurred: NO TASK WAS CREATED! (could be a library problem!) " << endl;
       return 0x0;
   }
   // set triggers
   if (bCentralTrigger) task->SelectCollisionCandidates(AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral);
   else                 task->SelectCollisionCandidates(AliVEvent::kMB);
   if(debug) cout << "    --> Set trigger selection to ";
   if(debug&&bCentralTrigger) cout << " kMB, kCentral, kSemiCentral " << endl;
   if(debug&&(!bCentralTrigger)) cout << " kMB " << endl;
   //set RP cuts for flow package analysis
   cutsRP = new AliFlowTrackCuts("RFPcuts");
   if(!cutsRP) {
       if(debug) cout << " Fatal error: no RP cuts found, could be a library problem! " << endl;
       return 0x0;
   }
   if(useGlobalRPCuts) {
       AliFlowTrackCuts::trackParameterType rptype = AliFlowTrackCuts::kGlobal;
       cutsRP->SetParamType(rptype);
       cutsRP->SetPtRange(0.2, 5.0);
       cutsRP->SetEtaRange(-0.8, 0.8);
       cutsRP->SetMinNClustersTPC(70);
       cutsRP->SetMinChi2PerClusterTPC(0.1);
       cutsRP->SetMaxChi2PerClusterTPC(4.0);
       cutsRP->SetRequireTPCRefit(kTRUE);
       cutsRP->SetMaxDCAToVertexXY(0.3);
       cutsRP->SetMaxDCAToVertexZ(0.3);
       cutsRP->SetAcceptKinkDaughters(kFALSE);
       cutsRP->SetMinimalTPCdedx(10.);
       if(debug) cout << "    --> kGlobal RP's " << cutsRP << endl;
   }
   if(!useGlobalRPCuts) {
       AliFlowTrackCuts::trackParameterType rptype = AliFlowTrackCuts::kTPCstandalone;
       cutsRP->SetParamType(rptype);
       cutsRP->SetPtRange(0.2,5.);
       cutsRP->SetEtaRange(-0.8,0.8);
       cutsRP->SetMinNClustersTPC(70);
       cutsRP->SetMinChi2PerClusterTPC(0.2);
       cutsRP->SetMaxChi2PerClusterTPC(4.0);
       cutsRP->SetMaxDCAToVertexXY(3.0);
       cutsRP->SetMaxDCAToVertexZ(3.0);
       cutsRP->SetDCAToVertex2D(kTRUE);
       cutsRP->SetAcceptKinkDaughters(kFALSE);
       cutsRP->SetMinimalTPCdedx(10.);
       if(debug) cout << "    --> kTPCStandAlone RP's " << cutsRP << endl;
   }
   task->SetRPCuts(cutsRP);
   // set POI cuts for kaon selection
   AliFlowTrackCuts* cutsPOI = new AliFlowTrackCuts("GlobalPOI");
   if(!cutsPOI) {
       if(debug) cout << " Fatal error: no POI cuts (could be a library problem)!" << endl;
       return 0x0;
   }
   AliFlowTrackCuts* cutsPOI = cutsPOI->GetStandardGlobalTrackCuts2010();
   cutsPOI->SetPtRange(0.2, 5); //?
   cutsPOI->SetMaxDCAToVertexXY(0.3);
   cutsPOI->SetMaxDCAToVertexZ(0.3);
   if(debug) cout << "    --> cutsPOI " << cutsPOI << endl;
   task->SetPOICuts(cutsPOI);
   //set POI cuts for aods XY Z - 3 distinct cases. comment out the unwanted case
   Double_t POIDCA[] = {0., 0., 0., 0., 0.};
   // 1 --- do nothing
   if(DCA == "none" ) {
       if (debug) cout << " --> No DCA cut on POI's <-- " << endl;
       for (Int_t i = 0; i < 5; i++) POIDCA[i] = 0.;
   }
   // 2 --- use fixed values for xy z
   if(DCA == "fix" ) {
       if (debug) cout << " --> Fixed DCA cut on POI's <-- " << endl;
       POIDCA[0] = -1.; POIDCA[1] = 0.3; POIDCA[2] = 0.3; POIDCA[3] = 0.; POIDCA[4] = 0.;
   }
   // 3 --- use pt dependent cut
   if(DCA == "pt" ) {
       if (debug) cout << " --> Pt dependent DCA cut on POI's <-- " << endl;
       POIDCA[0] = 1.; POIDCA[1] = 0.0105; POIDCA[2] = 0.0350; POIDCA[3] = 1.1; POIDCA[4] = 2.;
   }
   task->SetPOIDCAXYZ(POIDCA);
   // POI filter cuts, will filter invm mass bands and subevents
   AliFlowTrackSimpleCuts* POIfilterQC[30];
   AliFlowTrackSimpleCuts* POIfilterSP[30][2];
   Double_t flowBands[2][30];
   for (Int_t mb = 0; mb < 30; mb++) {
      flowBands[0][mb] = 0.99 + mb * 0.0034;
      flowBands[1][mb] = 0.99 + (mb + 1) * 0.0034;
      POIfilterSP[mb][0] = new AliFlowTrackSimpleCuts(Form("FilterPOISP_MB%d_ETANEG", mb));
      if(!POIfilterSP[mb][0]) {
          if(debug) cout << " FAILED to initialize POI filters, could be a library problem!" << endl;
          return 0x0;
      }
      POIfilterSP[mb][0]->SetEtaMin(-0.8);
      POIfilterSP[mb][0]->SetEtaMax(0.0);
      POIfilterSP[mb][0]->SetMassMin(flowBands[0][mb]);
      POIfilterSP[mb][0]->SetMassMax(flowBands[1][mb]);
      POIfilterSP[mb][1] = new AliFlowTrackSimpleCuts(Form("FilterPOISP_MB%d_ETAPOS", mb));
      POIfilterSP[mb][1]->SetEtaMin(0.0);
      POIfilterSP[mb][1]->SetEtaMax(+0.8);
      POIfilterSP[mb][1]->SetMassMin(flowBands[0][mb]);
      POIfilterSP[mb][1]->SetMassMax(flowBands[1][mb]);
      POIfilterQC[mb] = new AliFlowTrackSimpleCuts(Form("FilterPOIQC_MB%d", mb));
      POIfilterQC[mb]->SetEtaMin(-0.8);
      POIfilterQC[mb]->SetEtaMax(+0.8);
      POIfilterQC[mb]->SetMassMin(flowBands[0][mb]);
      POIfilterQC[mb]->SetMassMax(flowBands[1][mb]);
   }
   if(debug) cout << "    --> Created poi filters " << endl;
   task->SetCommonConstants(30, flowBands[0][0], flowBands[1][29]);
   if(debug) cout << "    --> Set common constants " << endl;
   // set pair and event cuts
   if((deltaDip>0.005)&&(deltaDipMaxPt>0.005)) task->SetMaxDeltaDipAngleAndPt(deltaDip, deltaDipMaxPt);
   else cout << " --> Disabled Delta-Dip exclusion. <-- " << endl;
   task->SetCandidateEtaAndPt(POIEtaMin, POIEtaMax, POIPtMin, POIPtMax);
   task->SetCentralityParameters(centrMin, centrMax, "TRK");
   task->SetVertexZ(vertexZ);
   if(debug) cout << "    --> Set pair cuts and event cuts" << endl;
   // set the kaon cuts, and specify the PID procedure which will be used
   task->SetPIDConfiguration(PIDconfig);
   if(debug) {
       cout << " ************ Configured PID routine ************ " << endl;
       cout << "    0 < " << PIDconfig[1] << " p_t, ITS || TPC with s < " << PIDconfig[0] << endl;
       if(PIDconfig[2] < 0.) cout << "  --> TPC control disabled " << endl;
       if(PIDconfig[2] > 0.) cout << "  --> TPC control enabled " << endl;
       cout << "    " << PIDconfig[1] << " < " << PIDconfig[4] << " p_t, TPC || ITS with s < " << PIDconfig[3] << endl;
       if(PIDconfig[5] < 0.) cout << "  --> ITS control disabled " << endl;
       if(PIDconfig[5] > 0.) cout << "  --> ITS control enabled " << endl;
       cout << "    " << PIDconfig[4] << " < 7 p_t, TPC / TOF Bayesian with p < " << PIDconfig[6] << endl;
       cout << " ************************************************ " << endl;
   }
   if (TPCStandAloneTracks)
   { // switch to tpc standalone tracks for POI selection
      if(debug) cout << "    --> Switching to TPC standalone analsis " << endl;
      task->SetRequireStrictKaonCuts();
      task->SetRequireTPCStandAloneKaons();
   }
   // create and configure IO containers
   AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
   AliAnalysisDataContainer *coutHist = mgr->CreateContainer(Form("PhiV2_%s", OutputName(centrMin, centrMax,PIDconfig, suffixName, bCentralTrigger, EtaGap, POIEtaMin, POIEtaMax, POIPtMin, POIPtMax, deltaDip, deltaDipMaxPt, DCA, harm, TPCStandAloneTracks, vertexZ, debug, useGlobalRPCuts).Data()), TList::Class(), AliAnalysisManager::kOutputContainer, fileName);
   if(debug) cout << "    --> Created IO containers " << cinput << ", " << coutHist << endl;
   mgr->AddTask(task);
   if(debug) cout << " === ADDING MAIN TASK == " << endl;
   mgr->ConnectInput(task, 0, cinput);
   mgr->ConnectOutput(task, 1, coutHist);
   if(debug) cout << "    --> Connected IO containers " << endl;
   if (SP || EP || QC || SPSUB) // if flow analysis should be done after reconstruction
   {
      if(debug) cout << " === RECEIVED REQUEST FOR FLOW ANALYSIS === " << endl;
      AliAnalysisDataContainer *flowEvent = mgr->CreateContainer(Form("FC%s", suffixName.Data()), AliFlowEventSimple::Class(), AliAnalysisManager::kExchangeContainer);
      mgr->ConnectOutput(task, 2, flowEvent);
      if(debug) cout << "    --> Created IO containers " << flowEvent << endl;
      if(debug) cout << "    --> suffixName " << suffixName << endl;
      for (int mb = 0; mb != 30; ++mb) {
         if (QC) {  // add qc tasks
            AddQCmethod(Form("QCTPCMB_%d_%s", mb, suffixName.Data()), harm, flowEvent, POIfilterQC[mb], debug, 0x0, suffixName.Data());
            if(debug) cout << "    --> Hanging QC task ... " << mb << " succes! "<< endl;
         }
         if (SPSUB) {  // add sp subevent tasks
            AddSPmethod(Form("SPTPCMB_%d_%s", mb, suffixName.Data()), -0.8, -.5*EtaGap, .5*EtaGap, +0.8, "Qa", harm, flowEvent, POIfilterSP[mb][1], NULL, false, shrinkSP, debug, true, suffixName.Data());
            if(debug) cout << "    --> Hanging SP Qa task " << mb << " succes!" << endl;
            AddSPmethod(Form("SPTPCMB_%d_%s", mb, suffixName.Data()), -0.8, -.5*EtaGap, .5*EtaGap, +0.8, "Qb", harm, flowEvent, POIfilterSP[mb][0], NULL, false, shrinkSP, debug, true, suffixName.Data());
            if(debug) cout << "    --> Hanging SP Qb task ..." << mb << " succes!"<< endl;
         }
         if (SP) { // add sp tasks
            AddSPmethod(Form("SPTPCMB_%d_%s", mb, suffixName.Data()), -0.8, -0.0, +0.0, +0.8, "QaQb", harm, flowEvent, POIfilterQC[mb], NULL, false, shrinkSP, debug, 0x0, suffixName.Data());
            if(debug) cout << "    --> Hanging SP task ... " << mb << " succes!" << endl;
         }
         if (EP) { // add ep tasks
            AddSPmethod(Form("EPTPCMB_%d_%s", mb, suffixName.Data()), -0.8, -0.0, +0.0, +0.8, "QaQb", harm, flowEvent, POIfilterQC[mb], NULL, true, shrinkSP, debug, 0x0, suffixName.Data());
            if(debug) cout << "    --> Hanging EP task ... " << mb << " succes!" << endl;
         }
      }
   }
   // print summary to screen
   cout << endl << endl << "       ==== AddTaskPhiFlow launched  ===== " << endl;
   cout << " ************ Configured PID routine ************ " << endl;
   cout << "      0 < " << PIDconfig[1] << " p_t, ITS || TPC with s < " << PIDconfig[0] << endl;
   if(PIDconfig[2] < 0.) cout << "    --> TPC control disabled " << endl;
   if(PIDconfig[2] > 0.) cout << "    --> TPC control enabled " << endl;
   cout << "    " << PIDconfig[1] << " < " << PIDconfig[4] << " p_t, TPC || ITS with s < " << PIDconfig[3] << endl;
   if(PIDconfig[5] < 0.) cout << "    --> ITS control disabled " << endl;
   if(PIDconfig[5] > 0.) cout << "    --> ITS control enabled " << endl;
   cout << "      " << PIDconfig[4] << " < 7 p_t, TPC / TOF Bayesian with p < " << PIDconfig[6] << endl;
   cout << " ************ Configured DCA setup ************** " << endl;
   cout << "  DCA type: " << DCA;
   if (DCA == "") cout << "default";
   cout << endl << " ************* Task statistics ****************** " << endl;
   cout << "   -> Launched PHI reconstruction " << endl;
   if(SP) cout << "   --> Launched 30 QaQb SP filters and corresponding 30 SP tasks " << endl;
   if(EP) cout << "   --> Launched 30 QaQb QC filters and corresponding 30 EP tasks " << endl;
   if(SPSUB) cout << "   --> Launched 30+30 Qa&&Qb SP filters and corresponding 30+30 SP tasks " << endl;
   if(QC) cout << "   --> Launched 30 QaQb QC filters and corresponding 30 QC tasks " << endl;
   cout << " ************************************************ " << endl;
   TString condit = "";
   (task->SetQA(kTRUE)) ? condit+= " --> Enabled QA plots <-- " : condit+= " --> Disabled QA plots <-- ";
   (task->SetIsMC(kFALSE)) ? condit+= " --> MC mode <-- " : condit+= " --> DATA mode <-- ";
   cout << condit << endl;
   cout << "           --> Now go for a coffee! <-- " << endl;
   cout << " ************************************************ " << endl;
   return task;
}
//_____________________________________________________________________________
void AddSPmethod(char *name, double minEtaA, double maxEtaA, double minEtaB, double maxEtaB, char *Qvector, int harmonic, AliAnalysisDataContainer *flowEvent, AliFlowTrackSimpleCuts *cutsPOI = NULL, AliFlowTrackSimpleCuts *cutsRFP = NULL, bool bEP, bool shrink = false, bool debug, bool etagap = false, char *suffixName)
{
   // add sp task and invm filter tasks
   if(debug) (bEP) ? cout << " ****** Reveived request for EP task ****** " << endl : cout << " ******* Switching to SP task ******* " << endl;
   TString fileName = AliAnalysisManager::GetCommonFileName();
   (bEP) ? fileName+=":EP" : fileName+=":SP";
   fileName+=suffixName;
   if(etagap) {
       fileName+="_SUBEVENTS";
       if(debug) cout << "    --> Setting up subevent analysis <-- " << endl;
   }
   if(debug) cout << "    --> fileName " << fileName << endl;
   TString myNameSP;
   (bEP) ? myNameSP = Form("%sEPv%d%s", name, harmonic, Qvector): myNameSP = Form("%sSPv%d%s", name, harmonic, Qvector);
   if(debug) cout << " Task and filter name: " << myNameSP << endl;
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   AliAnalysisDataContainer *flowEvent2 = mgr->CreateContainer(Form("Filter_%s", myNameSP.Data()), AliFlowEventSimple::Class(), AliAnalysisManager::kExchangeContainer);
   AliAnalysisTaskFilterFE *tskFilter = new AliAnalysisTaskFilterFE(Form("TaskFilter_%s", myNameSP.Data()), cutsRFP, cutsPOI);
   tskFilter->SetSubeventEtaRange(minEtaA, maxEtaA, minEtaB, maxEtaB);
   mgr->AddTask(tskFilter);
   mgr->ConnectInput(tskFilter, 0, flowEvent);
   mgr->ConnectOutput(tskFilter, 1, flowEvent2);
   AliAnalysisDataContainer *outSP = mgr->CreateContainer(myNameSP.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, fileName);
   AliAnalysisTaskScalarProduct *tskSP = new AliAnalysisTaskScalarProduct(Form("TaskScalarProduct_%s", myNameSP.Data()), kFALSE);
   tskSP->SetApplyCorrectionForNUA(kTRUE);
   tskSP->SetHarmonic(harmonic);
   tskSP->SetTotalQvector(Qvector);
   if (bEP) tskSP->SetBehaveAsEP();
   if (shrink) tskSP->SetBookOnlyBasicCCH(kTRUE);
   mgr->AddTask(tskSP);
   mgr->ConnectInput(tskSP, 0, flowEvent2);
   mgr->ConnectOutput(tskSP, 1, outSP);
}
//_____________________________________________________________________________
void AddQCmethod(char *name, int harmonic, AliAnalysisDataContainer *flowEvent, AliFlowTrackSimpleCuts *cutsPOI = NULL, Bool_t debug, AliFlowTrackSimpleCuts *cutsRFP = NULL, char *suffixName)
{
   // add qc task and invm filter tasks
   if(debug) cout << " ****** Received request for QC v" << harmonic << " task " << name << ", POI " << cutsPOI << ", IO ****** " << flowEvent << endl;
   TString fileName = AliAnalysisManager::GetCommonFileName();
   fileName+=":QC";
   fileName+=suffixName;
   if(debug) cout << "    --> Common filename: " << fileName << endl;
   TString myName = Form("%s", name);
   if(debug) cout << "    --> myName: " << myName << endl;
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   AliAnalysisDataContainer *flowEvent2 = mgr->CreateContainer(Form("Filter_%s", myName.Data()), AliFlowEventSimple::Class(),AliAnalysisManager::kExchangeContainer);
   AliAnalysisTaskFilterFE *tskFilter = new AliAnalysisTaskFilterFE(Form("TaskFilter_%s", myName.Data()), cutsRFP, cutsPOI);
   mgr->AddTask(tskFilter);
   mgr->ConnectInput(tskFilter, 0, flowEvent);
   mgr->ConnectOutput(tskFilter, 1, flowEvent2);
   AliAnalysisDataContainer *outQC = mgr->CreateContainer(myName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, fileName);
   AliAnalysisTaskQCumulants *tskQC = new AliAnalysisTaskQCumulants(Form("TaskQCumulants_%s", myName.Data()), kFALSE);
   tskQC->SetApplyCorrectionForNUA(kTRUE);
   tskQC->SetHarmonic(harmonic);
   tskQC->SetBookOnlyBasicCCH(kTRUE);
   mgr->AddTask(tskQC);
   mgr->ConnectInput(tskQC, 0, flowEvent2);
   mgr->ConnectOutput(tskQC, 1, outQC);
}
//_____________________________________________________________________________
TString OutputName( Float_t centrMin,
                    Float_t centrMax,
                    Double_t PIDconfig[7],
                    TString suffixName,
                    Bool_t bCentralTrigger,
                    Float_t EtaGap,
                    Float_t POIEtaMin,
                    Float_t POIEtaMax,
                    Float_t POIPtMin,
                    Float_t POIPtMax,
                    Float_t deltaDip,
                    Float_t deltaDipMaxPt,
                    TString DCA,
                    Int_t harm,
                    Bool_t TPCStandAloneTracks,
                    Float_t vertexZ,
                    Bool_t debug,
                    Bool_t useGlobalRPCuts)
{
   // generate output name
   TString centralityName = "";
   centralityName += suffixName;
   centralityName += "_DCA";
   centralityName += DCA;
   centralityName += Form("_vZ%.f", vertexZ);
   centralityName += "_";
   for(Int_t i = 0; i < 7; i++) centralityName += Form("%.1f_", PIDconfig[i]);
   centralityName += "_POIEta";
   centralityName += Form("%.1f", POIEtaMin);
   centralityName += Form("%.1f", POIEtaMax);
   centralityName += "_gap";
   centralityName += Form("%.1f", -0.5*EtaGap);
   centralityName += Form("%.1f", 0.5*EtaGap);
   centralityName += "_";
   centralityName += Form("dDip%.2f", deltaDip);
   centralityName += "-";
   centralityName += Form("dDipPt%.2f", deltaDipMaxPt);
   if (TPCStandAloneTracks) {
      centralityName += "-";
      centralityName += "TPCStandAloneTracks";
   }
   if (bCentralTrigger) {
      centralityName += "-";
      centralityName += "kMBkCkSC";
   }
   if (!useGlobalRPCuts) {
       centralityName += "-";
       centralityName += "TPCRP";
   }
   if(debug) cout << "    --> centralityName " << centralityName << endl;
   return centralityName;
}
//_____________________________________________________________________________

