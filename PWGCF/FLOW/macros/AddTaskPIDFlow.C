class AliAnalysisDataContainer;
class AliFlowTrackCuts;
class AliFlowTrackSimpleCuts;
class AliFlowEventCuts;
class AliFlowEventSimpleCuts;
class AliAnalysisDataContainer;


void AddTaskPIDFlow(Int_t triggerSelectionString=AliVEvent::kMB,
                                   Float_t etamin=-0.8,
                                   Float_t etamax=0.8,
                                   Float_t EtaGap=0.2,
                                   TString fileNameBase="output",
                                   TString uniqueStr="",
                                   Int_t AODfilterBitRP = 272,
                                   Int_t AODfilterBitPOI = 272,
                                   Int_t charge=0,
                                   Bool_t doQA=kTRUE,
                                   Bool_t isPID = kTRUE,
                                   Bool_t is2011 = kFALSE,
                                   AliPID::EParticleType particleType=AliPID::kPion,
                                   AliFlowTrackCuts::PIDsource sourcePID=AliFlowTrackCuts::kTOFbayesian) {
    
Bool_t debug = kTRUE;
 
// Define a range of the detector to exclude
Bool_t ExcludeRegion = kFALSE;
Double_t excludeEtaMin = -0.;
Double_t excludeEtaMax = 0.;
Double_t excludePhiMin = 0.;
Double_t excludePhiMax = 0.;
   
// RUN SETTINGS
// Flow analysis method can be:(set to kTRUE or kFALSE)
Bool_t SPSUB    = kTRUE;
Bool_t SP       = kTRUE;  // scalar product method (similar to eventplane method)
Bool_t QC       = kTRUE;  // cumulants using Q vectors
    
int centrMin[9] = {0,5,10,20,30,40,50,60,70};
int centrMax[9] = {5,10,20,30,40,50,60,70,80};
const int ncentr=6;
    
    
//---------Data selection----------
//kMC, kGlobal, kESD_TPConly, kESD_SPDtracklet
AliFlowTrackCuts::trackParameterType rptype = AliFlowTrackCuts::kTPCstandalone;
AliFlowTrackCuts::trackParameterType poitype = AliFlowTrackCuts::kTPCstandalone;
    
//---------Parameter mixing--------
//kPure - no mixing, kTrackWithMCkine, kTrackWithMCPID, kTrackWithMCpt
AliFlowTrackCuts::trackParameterMix rpmix = AliFlowTrackCuts::kPure;
AliFlowTrackCuts::trackParameterMix poimix = AliFlowTrackCuts::kPure;
    
const char* rptypestr = AliFlowTrackCuts::GetParamTypeName(rptype);
const char* poitypestr = AliFlowTrackCuts::GetParamTypeName(poitype);

    
//===========================================================================
// EVENTS CUTS:
AliFlowEventCuts* cutsEvent[ncentr];
AliFlowTrackCuts* cutsRP[ncentr];
AliFlowTrackCuts* cutsPOI[ncentr];
TString outputSlotName[ncentr][4];
for(int icentr=0;icentr<ncentr;icentr++){
    cutsEvent[icentr] = DefinecutsEvent();
    //cutsEvent[icentr]->SetUsedDataset(is2011);
    cutsEvent[icentr]->SetCentralityPercentileRange(centrMin[icentr],centrMax[icentr]);
    cutsEvent[icentr]->SetCentralityPercentileMethod(AliFlowEventCuts::kV0);
    //  cutsEvent->SetRefMultMethod(AliFlowEventCuts::kVZERO);
    //cutsEvent->SetCentralityPercentileMethod(AliFlowEventCuts::kSPD1tracklets);
    //cutsEvent->SetNContributorsRange(2);
    cutsEvent[icentr]->SetPrimaryVertexZrange(-10.,10.);
    cutsEvent[icentr]->SetQA(doQA);
    cutsEvent[icentr]->SetCutTPCmultiplicityOutliers();
    
    
    // RP TRACK CUTS:
    //  AliFlowTrackCuts* cutsRP2 = AliFlowTrackCuts::GetStandardVZEROOnlyTrackCuts();
    cutsRP[icentr] = DefineRPcuts();
    cutsRP[icentr]->SetParamType(rptype);
    cutsRP[icentr]->SetParamMix(rpmix);
    cutsRP[icentr]->SetPtRange(0.2,5.);
    cutsRP[icentr]->SetEtaRange(etamin,etamax);
    cutsRP[icentr]->SetMinNClustersTPC(70);
    //  cutsRP->SetMinChi2PerClusterTPC(0.1);//
    // cutsRP->SetMaxChi2PerClusterTPC(4.0);//
    cutsRP[icentr]->SetMaxDCAToVertexXY(3.0);
    cutsRP[icentr]->SetMaxDCAToVertexZ(3.0);
    cutsRP[icentr]->SetAcceptKinkDaughters(kFALSE);
    cutsRP[icentr]->SetMinimalTPCdedx(10);
    cutsRP[icentr]->SetAODfilterBit(AODfilterBitRP);
    cutsRP[icentr]->SetQA(doQA);
    
    // POI TRACK CUTS:
    cutsPOI[icentr] = DefinePOIcuts();
    cutsPOI[icentr]->GetBayesianResponse()->ForceOldDedx(); // for 2010 data to use old TPC PID Response instead of the official one
    cutsPOI[icentr]->SetParamType(poitype);
    cutsPOI[icentr]->SetParamMix(poimix);
    cutsPOI[icentr]->SetPtRange(0.2,5.);//
    cutsPOI[icentr]->SetEtaRange(etamin,etamax);
    //cutsPOI->SetRequireCharge(kTRUE);
    //cutsPOI->SetPID(PdgRP);
    cutsPOI[icentr]->SetMinNClustersTPC(70);
    // cutsPOI->SetMinChi2PerClusterTPC(0.1); //
    // cutsPOI->SetMaxChi2PerClusterTPC(4.0); //
    //  cutsPOI->SetRequireITSRefit(kTRUE);
    //  cutsPOI->SetRequireTPCRefit(kTRUE);
    //  cutsPOI->SetMinNClustersITS(2);
    //cutsPOI->SetMaxChi2PerClusterITS(1.e+09);
    cutsPOI[icentr]->SetMaxDCAToVertexXY(3.0);
    cutsPOI[icentr]->SetMaxDCAToVertexZ(3.0);
    //cutsPOI->SetDCAToVertex2D(kTRUE);
    //cutsPOI->SetMaxNsigmaToVertex(1.e+10);
    //cutsPOI->SetRequireSigmaToVertex(kFALSE);
    cutsPOI[icentr]->SetAcceptKinkDaughters(kFALSE);
    if(isPID) cutsPOI[icentr]->SetPID(particleType, sourcePID);//particleType, sourcePID
    if (charge!=0) cutsPOI[icentr]->SetCharge(charge);
    //cutsPOI->SetAllowTOFmismatch(kFALSE);
    cutsPOI[icentr]->SetRequireStrictTOFTPCagreement(kTRUE);
    //iexample: francesco's tunig TPC Bethe Bloch for data:
     cutsPOI[icentr]->SetMinimalTPCdedx(10);
    cutsPOI[icentr]->SetAODfilterBit(AODfilterBitPOI);
    // cutsPOI->SetAODfilterBit(768);
    cutsPOI[icentr]->SetQA(doQA);
    cutsPOI[icentr]->SetPriors((centrMin[icentr]+centrMax[icentr])*0.5); // set priors and PID as a function of the centrality

    //=====================================================================
 
    TString suffixName = "highharmflow";
    suffixName += Form("%i", centrMin[icentr]);
    suffixName += Form("%i", centrMax[icentr]);
    
    for(int harmonic=2;harmonic<6;harmonic++){  //for v2,v3,v4 and v5
        outputSlotName[icentr][harmonic-2] = "";
        outputSlotName[icentr][harmonic-2]+=uniqueStr;
        outputSlotName[icentr][harmonic-2]+=Form("v%i",harmonic);
        outputSlotName[icentr][harmonic-2]+=cutsRP[icentr]->GetName();
        outputSlotName[icentr][harmonic-2]+="_";
        outputSlotName[icentr][harmonic-2]+=cutsPOI[icentr]->GetName();
        outputSlotName[icentr][harmonic-2]+=Form("_%i-",centrMin[icentr]);
        outputSlotName[icentr][harmonic-2]+=Form("%i_",centrMax[icentr]);
        
        
        if(isPID){
            outputSlotName[icentr][harmonic-2]+=AliFlowTrackCuts::PIDsourceName(sourcePID);//sourcePID
            outputSlotName[icentr][harmonic-2]+="_";
            outputSlotName[icentr][harmonic-2]+=AliPID::ParticleName(particleType);//particleType
        }
        else{
            outputSlotName[icentr][harmonic-2]+="AllCharged";
        }
        if (charge<0) outputSlotName[icentr][harmonic-2]+="-";
        if (charge>0) outputSlotName[icentr][harmonic-2]+="+";
            }
}
    
TString fileName(fileNameBase);
fileName.Append(".root");

 
    
//====================================FLOWPACKAGE TASKS=========================//
AliAnalysisDataContainer *cinput1[ncentr];
AliAnalysisDataContainer *coutputFE[ncentr];
AliAnalysisDataContainer* coutputFEQA[ncentr];
AliAnalysisTaskFlowEvent *taskFE[ncentr];
    
AliFlowTrackCuts *QC_POI[ncentr], *SP_POI[ncentr][2], *SP_POI_gap[ncentr][2]; // POIs
    
TString outputQA[ncentr];

for (int icentr=0; icentr<ncentr; icentr++) {
        
    // Get the pointer to the existing analysis manager via the static access method.
        //==============================================================================
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        Error("AddTaskFlowEvent", "No analysis manager to connect to.");
        return NULL;
    }
        
    // Check the analysis type using the event handlers connected to the analysis
    // manager. The availability of MC handler can also be checked here.
    //==============================================================================
    if (!mgr->GetInputEventHandler()) {
        ::Error("AddTaskFlowEvent", "This task requires an input event handler");
        return NULL;
    }

    taskFE[icentr] = new AliAnalysisTaskFlowEvent(Form("TaskFlowEvent_%s",suffixName.Data()),"",doQA);
    taskFE[icentr]->SelectCollisionCandidates(triggerSelectionString);
    taskFE[icentr]->SetSubeventEtaRange(etamin, -0.5*EtaGap, 0.5*EtaGap, etamax);
    mgr->AddTask(taskFE[icentr]);

    // Pass cuts for RPs and POIs to the task:
    taskFE[icentr]->SetCutsEvent(cutsEvent[icentr]);
    taskFE[icentr]->SetCutsRP(cutsRP[icentr]);
    taskFE[icentr]->SetCutsPOI(cutsPOI[icentr]);
    if (cutsRP[icentr]->GetParamType()==AliFlowTrackCuts::kVZERO)
    {
        //TODO: since this is set in a static object all analyses in an analysis train
        //will be affected.
        taskFE[icentr]->SetHistWeightvsPhiMin(0.);
        taskFE[icentr]->SetHistWeightvsPhiMax(200.);
    }
    cinput1[icentr] = mgr->GetCommonInputContainer();
    
    coutputFE[icentr] = mgr->CreateContainer(Form("FlowEventSimple_%s_%d",suffixName.Data(),icentr),AliFlowEventSimple::Class(),AliAnalysisManager::kExchangeContainer);
    mgr->ConnectInput(taskFE[icentr],0,cinput1[icentr]);
    mgr->ConnectOutput(taskFE[icentr],1,coutputFE[icentr]);

    if (taskFE[icentr]->GetQAOn()) {
        outputQA[icentr] = fileName;
        outputQA[icentr] += ":QA";
        coutputFEQA[icentr] = mgr->CreateContainer(Form("QA_%s_%d",suffixName.Data(),icentr), TList::Class(),AliAnalysisManager::kOutputContainer,outputQA[icentr]);
        mgr->ConnectOutput(taskFE[icentr],2,coutputFEQA[icentr]);
    }


    //POIs for SP and QC method
    //===========================================================================
 
    //half window for POIs
    
    QC_POI[icentr] = new AliFlowTrackCuts(Form("Filter_%d",icentr));
    QC_POI[icentr]->SetEtaRange( -0.8,+0.8 );
    QC_POI[icentr]->GetBayesianResponse()->ForceOldDedx(); // for 2010 data to use old TPC PID Response instead of the official one
    QC_POI[icentr]->SetParamType(poitype);
//    QC_POI[icentr]->SetParamMix(poimix);
    QC_POI[icentr]->SetPtRange(0.2,5.);//
    QC_POI[icentr]->SetMinNClustersTPC(70);
    // SP_POI->SetMinChi2PerClusterTPC(0.1); //
    // SP_POI->SetMaxChi2PerClusterTPC(4.0); //
    //  SP_POI->SetRequireITSRefit(kTRUE);
    //  SP_POI->SetRequireTPCRefit(kTRUE);
    //  SP_POI->SetMinNClustersITS(2);
    //  SP_POI->SetMaxChi2PerClusterITS(1.e+09);
    QC_POI[icentr]->SetMaxDCAToVertexXY(3.0);
    QC_POI[icentr]->SetMaxDCAToVertexZ(3.0);
    //SP_POI->SetDCAToVertex2D(kTRUE);
    //SP_POI->SetMaxNsigmaToVertex(1.e+10);
    //SP_POI->SetRequireSigmaToVertex(kFALSE);
    QC_POI[icentr]->SetAcceptKinkDaughters(kFALSE);
    if(isPID) QC_POI[icentr]->SetPID(particleType, sourcePID);//particleType, sourcePID
    if (charge!=0) QC_POI[icentr]->SetCharge(charge);
    //SP_POI->SetAllowTOFmismatch(kFALSE);
    QC_POI[icentr]->SetRequireStrictTOFTPCagreement(kTRUE);
    QC_POI[icentr]->SetMinimalTPCdedx(10);
    QC_POI[icentr]->SetAODfilterBit(AODfilterBitPOI);
    QC_POI[icentr]->SetQA(doQA);
    QC_POI[icentr]->SetPriors((centrMin[icentr]+centrMax[icentr])*0.5);

    if(SP){
        SP_POI[icentr][0] = new AliFlowTrackCuts(Form("Filterhf0_%d",icentr));
        SP_POI[icentr][1] = new AliFlowTrackCuts(Form("Filterhf1_%d",icentr));
        for(int hw=0;hw<2;hw++){
            SP_POI[icentr][hw]->GetBayesianResponse()->ForceOldDedx(); // for 2010 data to use old TPC PID Response instead of the official one
            SP_POI[icentr][hw]->SetParamType(poitype);
            //      SP_POI[icentr][hw]->SetParamMix(poimix);
            SP_POI[icentr][hw]->SetPtRange(0.2,5.);//
            SP_POI[icentr][hw]->SetMinNClustersTPC(70);
            // SP_POI->SetMinChi2PerClusterTPC(0.1); //
            // SP_POI->SetMaxChi2PerClusterTPC(4.0); //
            //  SP_POI->SetRequireITSRefit(kTRUE);
            //  SP_POI->SetRequireTPCRefit(kTRUE);
            //  SP_POI->SetMinNClustersITS(2);
            //  SP_POI->SetMaxChi2PerClusterITS(1.e+09);
            SP_POI[icentr][hw]->SetMaxDCAToVertexXY(3.0);
            SP_POI[icentr][hw]->SetMaxDCAToVertexZ(3.0);
            //SP_POI->SetDCAToVertex2D(kTRUE);
            //SP_POI->SetMaxNsigmaToVertex(1.e+10);
            //SP_POI->SetRequireSigmaToVertex(kFALSE);
            SP_POI[icentr][hw]->SetAcceptKinkDaughters(kFALSE);
            if(isPID) SP_POI[icentr][hw]->SetPID(particleType, sourcePID);//particleType, sourcePID
            if (charge!=0) SP_POI[icentr][hw]->SetCharge(charge);
            //SP_POI->SetAllowTOFmismatch(kFALSE);
            SP_POI[icentr][hw]->SetRequireStrictTOFTPCagreement(kTRUE);
            SP_POI[icentr][hw]->SetMinimalTPCdedx(10);
            SP_POI[icentr][hw]->SetAODfilterBit(AODfilterBitPOI);
            SP_POI[icentr][hw]->SetQA(doQA);
            SP_POI[icentr][hw]->SetPriors((centrMin[icentr]+centrMax[icentr])*0.5);
        }
    }
    if(SPSUB){
        SP_POI_gap[icentr][0] = new AliFlowTrackCuts(Form("Filterhf0_%d",icentr));
        SP_POI_gap[icentr][1] = new AliFlowTrackCuts(Form("Filterhf1_%d",icentr));
        for(int hw=0;hw<2;hw++){
            SP_POI_gap[icentr][hw]->GetBayesianResponse()->ForceOldDedx(); // for 2010 data to use old TPC PID Response instead of the official one
            SP_POI_gap[icentr][hw]->SetParamType(poitype);
            //      SP_POI[icentr][hw]->SetParamMix(poimix);
            SP_POI_gap[icentr][hw]->SetPtRange(0.2,5.);//
            SP_POI_gap[icentr][hw]->SetMinNClustersTPC(70);
            // SP_POI_gap->SetMinChi2PerClusterTPC(0.1); //
            // SP_POI_gap->SetMaxChi2PerClusterTPC(4.0); //
            //  SP_POI_gap->SetRequireITSRefit(kTRUE);
            //  SP_POI_gap->SetRequireTPCRefit(kTRUE);
            //  SP_POI_gap->SetMinNClustersITS(2);
            //  SP_POI_gap->SetMaxChi2PerClusterITS(1.e+09);
            SP_POI_gap[icentr][hw]->SetMaxDCAToVertexXY(3.0);
            SP_POI_gap[icentr][hw]->SetMaxDCAToVertexZ(3.0);
            //SP_POI_gap->SetDCAToVertex2D(kTRUE);
            //SP_POI_gap->SetMaxNsigmaToVertex(1.e+10);
            //SP_POI_gap->SetRequireSigmaToVertex(kFALSE);
            SP_POI_gap[icentr][hw]->SetAcceptKinkDaughters(kFALSE);
            if(isPID) SP_POI_gap[icentr][hw]->SetPID(particleType, sourcePID);//particleType, sourcePID
            if (charge!=0) SP_POI_gap[icentr][hw]->SetCharge(charge);
            //SP_POI_gap->SetAllowTOFmismatch(kFALSE);
            SP_POI_gap[icentr][hw]->SetRequireStrictTOFTPCagreement(kTRUE);
            SP_POI_gap[icentr][hw]->SetMinimalTPCdedx(10);
            SP_POI_gap[icentr][hw]->SetAODfilterBit(AODfilterBitPOI);
            SP_POI_gap[icentr][hw]->SetQA(doQA);
            SP_POI_gap[icentr][hw]->SetPriors((centrMin[icentr]+centrMax[icentr])*0.5);
        }
    }

    

    for(int i=0;i<4;i++){
        int harm = i+2;
        if(QC) {
            AddQCmethod( "QC", fileName, rptypestr, outputSlotName[icentr][i], triggerSelectionString, coutputFE[icentr], NULL, QC_POI[icentr], harm, icentr); // QC TPC
        }
        if(SP) {

            AddSPmethod( "SP", fileName, rptypestr, outputSlotName[icentr][i], triggerSelectionString, coutputFE[icentr], NULL, SP_POI[icentr][0], "Qa", harm, icentr, etamin, etamax,0); // SP Qa
            AddSPmethod( "SP", fileName, rptypestr, outputSlotName[icentr][i], triggerSelectionString, coutputFE[icentr], NULL, SP_POI[icentr][1], "Qb", harm, icentr, etamin, etamax,0); // SP Qb
        }
        if(SPSUB && EtaGap!=0){
            AddSPmethod( "SP", fileName, rptypestr, outputSlotName[icentr][i], triggerSelectionString, coutputFE[icentr], NULL, SP_POI_gap[icentr][0], "Qa", harm, icentr, etamin, etamax, EtaGap); // SP Qa
            AddSPmethod( "SP", fileName, rptypestr, outputSlotName[icentr][i], triggerSelectionString, coutputFE[icentr], NULL, SP_POI_gap[icentr][1], "Qb", harm, icentr, etamin, etamax, EtaGap); // SP Qb
        }
    }


}
}
//AddSPmthod
//AddQCmethod
//===========================================================================
//QC method
void AddQCmethod(char *name="",TString fileName="", TString rptypestr="", TString outputSlotName="", Int_t triggerSelectionString=0,AliAnalysisDataContainer* coutputFE,AliFlowTrackCuts* cutsRP, AliFlowTrackCuts* cutsPOI,Int_t harm=2, Int_t icentr=0) {
    
    AliAnalysisDataContainer *coutputQC;
    AliAnalysisTaskQCumulants *taskQC;
    
    int centrMin[9] = {0,5,10,20,30,40,50,60,70};
    int centrMax[9] = {5,10,20,30,40,50,60,70,80};
    TString folder = "FlowPID_QC_";
    TString myFolder = Form("%sv%d_%d%d",folder.Data(),harm,centrMin[icentr],centrMax[icentr]);
    TString myNameQC = Form("%sv%d%s",name,harm,outputSlotName.Data());
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    AliAnalysisDataContainer *flowEvent = mgr->CreateContainer( Form("TPConly_%s",myFolder.Data()),
                                                               AliFlowEventSimple::Class(),
                                                               AliAnalysisManager::kExchangeContainer );
    
    AliAnalysisDataContainer *flowEvent2 = mgr->CreateContainer( Form("Filter_%s", myNameQC.Data()),
                                                                AliFlowEventSimple::Class(),
                                                                AliAnalysisManager::kExchangeContainer );
    AliAnalysisTaskFilterFE *tskFilter = new AliAnalysisTaskFilterFE( Form("TaskFilter_%s",myNameQC.Data()),cutsRP, cutsPOI);
    mgr->AddTask(tskFilter);
    mgr->ConnectInput( tskFilter,0,coutputFE);
    mgr->ConnectOutput(tskFilter,1,flowEvent2);
    
    taskQC = new AliAnalysisTaskQCumulants(Form("TaskQCumulants_%s",outputSlotName.Data()));
    taskQC->SelectCollisionCandidates(triggerSelectionString);
    taskQC->SetCalculateCumulantsVsM(kFALSE);
    taskQC->SetnBinsMult(10000);
    taskQC->SetMinMult(0.);
    taskQC->SetMaxMult(10000.);
    taskQC->SetHarmonic(harm);
    taskQC->SetApplyCorrectionForNUA(kFALSE);
    taskQC->SetFillMultipleControlHistograms(kFALSE);
    
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

    TString outputQC = fileName;
    outputQC += ":outputQCanalysis";
    outputQC+= rptypestr;
    
    coutputQC = mgr->CreateContainer(Form("QC_%s",outputSlotName.Data()),
                                        TList::Class(),AliAnalysisManager::kOutputContainer,outputQC);
    mgr->AddTask(taskQC);
    mgr->ConnectInput(taskQC,0,flowEvent2);
    mgr->ConnectOutput(taskQC,1,coutputQC);
}
//SPmethod

void AddSPmethod(char *name="", TString fileName="", TString rptypestr="", TString outputSlotName="", Int_t triggerSelectionString=0,AliAnalysisDataContainer* coutputFE, AliFlowTrackCuts* cutsRP,AliFlowTrackCuts* cutsPOI, char* Qvector, Int_t harm=2, Int_t icentr=0, Double_t etamin=-0.8, Double_t etamax=0.8, Double_t gap=0.0) {
 
    AliAnalysisDataContainer *coutputSP;
    AliAnalysisTaskScalarProduct *taskSP;

    int centrMin[9] = {0,5,10,20,30,40,50,60,70};
    int centrMax[9] = {5,10,20,30,40,50,60,70,80};
    
    if(Qvector=="Qa"){
        cutsPOI->SetEtaRange( +0.5*gap, +0.8 );
    }
    if(Qvector=="Qb"){
        cutsPOI->SetEtaRange( -0.8,-0.5*gap );
    }
    
    TString folder = "FlowPID_SP_";
    TString myFolder = Form("%sv%i_%s_%i-%i_%.f",folder.Data(),harm,Qvector,centrMin[icentr],centrMax[icentr],gap*10);
    TString myNameSP = Form("%s_v%d_%s_%s_%.f",name,harm,Qvector,outputSlotName.Data(),gap*10);
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    
    AliAnalysisDataContainer *flowEvent = mgr->CreateContainer( Form("TPConly_%s",myFolder.Data()),
                                                             AliFlowEventSimple::Class(),
                                                             AliAnalysisManager::kExchangeContainer );
    
    AliAnalysisDataContainer *flowEvent2 = mgr->CreateContainer( Form("Filter_%s", myNameSP.Data()),
                                                                AliFlowEventSimple::Class(),
                                                                AliAnalysisManager::kExchangeContainer );
    
    AliAnalysisTaskFilterFE *tskFilter = new AliAnalysisTaskFilterFE( Form("TaskFilter_%s",myNameSP.Data()),
                                                                     cutsRP, cutsPOI);
    tskFilter->SetSubeventEtaRange(etamin, -.5*gap, +.5*gap, etamax);
    mgr->AddTask(tskFilter);
    mgr->ConnectInput( tskFilter,0,coutputFE);
    mgr->ConnectOutput(tskFilter,1,flowEvent2);

    
    taskSP = new AliAnalysisTaskScalarProduct(Form("TaskScalarProduct_%s",outputSlotName.Data()),kFALSE);
    taskSP->SetHarmonic(harm);
    taskSP->SelectCollisionCandidates(triggerSelectionString);
    taskSP->SetRelDiffMsub(1.0);
    taskSP->SetTotalQvector(Qvector);
    taskSP->SetApplyCorrectionForNUA(kTRUE);
    
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();


    TString outputSP = fileName;
    outputSP += ":outputSPanalysis";
    outputSP+= rptypestr;
    coutputSP = mgr->CreateContainer(Form("SP_%s_%s_%.f",outputSlotName.Data(),Qvector,gap*10),
                                        TList::Class(),AliAnalysisManager::kOutputContainer,outputSP);
    mgr->AddTask(taskSP);
    mgr->ConnectInput(taskSP,0,flowEvent2);
    mgr->ConnectOutput(taskSP,1,coutputSP);
}
AliFlowEventCuts* DefinecutsEvent(){
    AliFlowEventCuts* cutsEvent = new AliFlowEventCuts("eventcuts");
    return cutsEvent;
}
AliFlowTrackCuts* DefineRPcuts(){
    AliFlowTrackCuts* cutsRP = new AliFlowTrackCuts("TPConlyRP");
    return cutsRP;
}
AliFlowTrackCuts* DefinePOIcuts(){
    AliFlowTrackCuts* cutsPOI = new AliFlowTrackCuts("TPConlyPOI");
    return cutsPOI;
}


