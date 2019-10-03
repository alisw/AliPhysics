//AddTask macro for higher harmonic flow analysis with QC method only.

class AliAnalysisDataContainer;
class AliFlowTrackCuts;
class AliFlowEventCuts;


void AddTaskPIDFlowQC(Int_t triggerSelectionString=AliVEvent::kMB,
                                   Int_t uptoWhichHarmonics = 2, // 2 --> v2 only, 3 --> v2 and v3, and so on
                                   Float_t etamin=-0.8,
                                   Float_t etamax=0.8,
                                   TString fileNameBase="AnalysisResults",
                                   TString uniqueStr="Pion",
                                   Int_t AODfilterBit = 272,
                                   Int_t charge=0,
                                   Int_t MinTPCdedx = 10,
                                   Int_t ncentralityminlim = 0,//0 start from 0-1cc
                                   Int_t ncentralitymaxlim = 6,//6 run over 0-50%
                                   Bool_t doQA=kTRUE,
                                   Bool_t isPID = kTRUE,
                                   Bool_t is2011 = kFALSE,
                                   Bool_t isAOD = kTRUE,
                                   Bool_t UsePIDParContours = kFALSE,
				   Bool_t UseOldDedx = kTRUE,
                                   AliPID::EParticleType particleType=AliPID::kPion,
                                   AliFlowTrackCuts::PIDsource sourcePID=AliFlowTrackCuts::kTOFbayesian) {

// Define a range of the detector to exclude
Bool_t ExcludeRegion = kFALSE;
Double_t excludeEtaMin = -0.;
Double_t excludeEtaMax = 0.;
Double_t excludePhiMin = 0.;
Double_t excludePhiMax = 0.;
    
int centrMin[8] = {0,0,10,20,30,40,60,60};
int centrMax[8] = {1,2,20,30,40,50,70,80};
const int ncentrminlim = ncentralityminlim;
const int ncentrmaxlim = ncentralitymaxlim;

    
    
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
const int ncentr = ncentrmaxlim - ncentrminlim;
const int nharmonics = uptoWhichHarmonics-1;
AliFlowEventCuts* cutsEvent[ncentr];
AliFlowTrackCuts* cutsRP[ncentr];
AliFlowTrackCuts* cutsPOI[ncentr];
TString outputSlotName[ncentr][nharmonics];
TString suffixName[ncentr];

    
for(int icentr=0;icentr<ncentr;icentr++){
    cutsEvent[icentr] = new AliFlowEventCuts(Form("eventcuts_%d",icentr));
    cutsEvent[icentr]->SetLHC11h(is2011);
    cutsEvent[icentr]->SetCentralityPercentileRange(centrMin[icentr+ncentrminlim],centrMax[icentr+ncentrminlim]);
    cutsEvent[icentr]->SetCentralityPercentileMethod(AliFlowEventCuts::kV0);
    cutsEvent[icentr]->SetPrimaryVertexZrange(-10.,10.);
    cutsEvent[icentr]->SetQA(doQA);
    cutsEvent[icentr]->SetCutTPCmultiplicityOutliers();
    
    
    // RP TRACK CUTS:
    cutsRP[icentr] = new AliFlowTrackCuts(Form("TPConlyRP_%d",icentr));
   // cutsRP[icentr]->SetParamType(rptype);
   // cutsRP[icentr]->SetParamMix(rpmix);
    cutsRP[icentr]->SetPtRange(0.2,5.);
    cutsRP[icentr]->SetEtaRange(etamin,etamax);
    cutsRP[icentr]->SetMinNClustersTPC(70);
    //  cutsRP->SetMinChi2PerClusterTPC(0.1);//
    // cutsRP->SetMaxChi2PerClusterTPC(4.0);//
    cutsRP[icentr]->SetMaxDCAToVertexXY(2.4);
    cutsRP[icentr]->SetMaxDCAToVertexZ(3.0);
    cutsRP[icentr]->SetAcceptKinkDaughters(kFALSE);
    cutsRP[icentr]->SetMinimalTPCdedx(MinTPCdedx);
    cutsRP[icentr]->SetAODfilterBit(AODfilterBit);
    cutsRP[icentr]->SetQA(doQA);

    
    //POIs for QC method
    //===========================================================================
    AliFlowTrackCuts  *QC_POI[ncentr];
    //half window for POIs
    //=======================QC POI Cuts
    QC_POI[icentr] = DefinePOIcuts();

    if(UseOldDedx) QC_POI[icentr]->GetBayesianResponse()->ForceOldDedx(); // for 2010 data to use old TPC PID Response instead of the official one
   // QC_POI[icentr]->SetParamType(poitype);
   // QC_POI[icentr]->SetParamMix(poimix);
    QC_POI[icentr]->SetPtRange(0.2,5.);//
    QC_POI[icentr]->SetMinNClustersTPC(70);
    QC_POI[icentr]->SetEtaRange( etamin,etamax );
  
    // QC_POI[icentr]->SetMinChi2PerClusterTPC(0.1); //
    // QC_POI[icentr]->SetMaxChi2PerClusterTPC(4.0); //
    //  QC_POI[icentr]->SetRequireITSRefit(kTRUE);
    //  QC_POI[icentr]->SetRequireTPCRefit(kTRUE);
    //  QC_POI[icentr]->SetMinNClustersITS(2);
    //  QC_POI[icentr]->SetMaxChi2PerClusterITS(1.e+09);
    QC_POI[icentr]->SetMaxDCAToVertexXY(2.4);
    QC_POI[icentr]->SetMaxDCAToVertexZ(3.0);
    //QC_POI[icentr]->SetDCAToVertex2D(kTRUE);
    //QC_POI[icentr]->SetMaxNsigmaToVertex(1.e+10);
    //QC_POI[icentr]->SetRequireSigmaToVertex(kFALSE);
    QC_POI[icentr]->SetAcceptKinkDaughters(kFALSE);
    if(isPID){
        QC_POI[icentr]->SetPID(particleType, sourcePID);//particleType, sourcePID
        QC_POI[icentr]->SetTPCTOFNsigmaPIDCutContours(UsePIDParContours,centrMin[icentr+ncentrminlim],centrMax[icentr+ncentrminlim]);
    }
    if (charge!=0) QC_POI[icentr]->SetCharge(charge);
    //QC_POI[icentr]->SetAllowTOFmismatch(kFALSE);
    QC_POI[icentr]->SetRequireStrictTOFTPCagreement(kTRUE);
    QC_POI[icentr]->SetMinimalTPCdedx(MinTPCdedx);
    if(isAOD) QC_POI[icentr]->SetAODfilterBit(AODfilterBit);
    QC_POI[icentr]->SetAODfilterBit(AODfilterBit);
    QC_POI[icentr]->SetQA(doQA);
        QC_POI[icentr]->SetPriors((centrMin[icentr+ncentrminlim]+centrMax[icentr+ncentrminlim])*0.5);



    //=====================================================================
 
    suffixName[icentr] = "flow";
    suffixName[icentr] += Form("%i_", centrMin[icentr+ncentrminlim]);
    suffixName[icentr] += Form("%i_", centrMax[icentr+ncentrminlim]);

    if(isPID){
        suffixName[icentr]+=AliFlowTrackCuts::PIDsourceName(sourcePID);
        //suffixName[icentr]+="_";
        //suffixName[icentr]+=AliPID::ParticleName(particleType);//particleType
    }
    else{
        suffixName[icentr]+="AllCharged";
    }
    if (charge<0) suffixName[icentr]+="-";
    if (charge>0) suffixName[icentr]+="+";


    for(int harmonic=2;harmonic<uptoWhichHarmonics+1;harmonic++){  //for v2,v3,v4 and v5
        outputSlotName[icentr][harmonic-2] = "";
        outputSlotName[icentr][harmonic-2]+=uniqueStr;
        outputSlotName[icentr][harmonic-2]+=Form("_v%i_",harmonic);
        outputSlotName[icentr][harmonic-2]+=cutsRP[icentr]->GetName();
        outputSlotName[icentr][harmonic-2]+="_";
        outputSlotName[icentr][harmonic-2]+=QC_POI[icentr]->GetName();
        outputSlotName[icentr][harmonic-2]+=Form("_%i-",centrMin[icentr+ncentrminlim]);
        outputSlotName[icentr][harmonic-2]+=Form("%i_",centrMax[icentr+ncentrminlim]);
        
        
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


AliAnalysisDataContainer *coutputQC[ncentr][nharmonics];
AliAnalysisTaskQCumulants *taskQC[ncentr][nharmonics];

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

    taskFE[icentr] = new AliAnalysisTaskFlowEvent(Form("TaskFlowEvent_%s",suffixName[icentr].Data()),"",doQA);
    taskFE[icentr]->SelectCollisionCandidates(triggerSelectionString);
    mgr->AddTask(taskFE[icentr]);

    // Pass cuts for RPs and POIs to the task:
    taskFE[icentr]->SetCutsEvent(cutsEvent[icentr]);
    taskFE[icentr]->SetCutsRP(cutsRP[icentr]);
    taskFE[icentr]->SetCutsPOI(QC_POI[icentr]);
    if (cutsRP[icentr]->GetParamType()==AliFlowTrackCuts::kVZERO)
    {
        //TODO: since this is set in a static object all analyses in an analysis train
        //will be affected.
        taskFE[icentr]->SetHistWeightvsPhiMin(0.);
        taskFE[icentr]->SetHistWeightvsPhiMax(200.);
    }
    
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    
    cinput1[icentr] = mgr->GetCommonInputContainer();
    
    coutputFE[icentr] = mgr->CreateContainer(Form("FlowEvent_%s",suffixName[icentr].Data()),AliFlowEventSimple::Class(),AliAnalysisManager::kExchangeContainer);

    mgr->ConnectInput(taskFE[icentr],0,cinput1[icentr]);
    mgr->ConnectOutput(taskFE[icentr],1,coutputFE[icentr]);
    //==========================================================

    for(int harm=2;harm<uptoWhichHarmonics+1;harm++){
        
        taskQC[icentr][harm-2] = new AliAnalysisTaskQCumulants(Form("TaskQCumulants_%s",outputSlotName[icentr][harm-2].Data()),kFALSE);
        taskQC[icentr][harm-2]->SelectCollisionCandidates(triggerSelectionString);
      //  taskQC[icentr][harm-2]->SetUsePhiWeights(WEIGHTS[0]);
      //  taskQC[icentr][harm-2]->SetUsePtWeights(WEIGHTS[1]);
      //  taskQC[icentr][harm-2]->SetUseEtaWeights(WEIGHTS[2]);
        taskQC[icentr][harm-2]->SetCalculateCumulantsVsM(kFALSE);
        taskQC[icentr][harm-2]->SetnBinsMult(10000);
        taskQC[icentr][harm-2]->SetMinMult(0.);
        taskQC[icentr][harm-2]->SetMaxMult(10000.);
        taskQC[icentr][harm-2]->SetHarmonic(harm);
        taskQC[icentr][harm-2]->SetApplyCorrectionForNUA(kFALSE);
        taskQC[icentr][harm-2]->SetFillMultipleControlHistograms(kFALSE);
        mgr->AddTask(taskQC[icentr][harm-2]);

        TString outputQC = fileName;
        outputQC += ":outputQCanalysis";
        outputQC+= rptypestr;
        
        coutputQC[icentr][harm-2] = mgr->CreateContainer(Form("QC_%s",outputSlotName[icentr][harm-2].Data()),
                                            TList::Class(),AliAnalysisManager::kOutputContainer,outputQC);
        mgr->ConnectInput(taskQC[icentr][harm-2],0,coutputFE[icentr]);
        mgr->ConnectOutput(taskQC[icentr][harm-2],1,coutputQC[icentr][harm-2]);


    }


    if (taskFE[icentr]->GetQAOn()) {
        outputQA[icentr] = fileName;
        outputQA[icentr] += ":QA";
        coutputFEQA[icentr] = mgr->CreateContainer(Form("QA_%s",suffixName[icentr].Data()), TList::Class(),AliAnalysisManager::kOutputContainer,outputQA[icentr]);
        mgr->ConnectOutput(taskFE[icentr],2,coutputFEQA[icentr]);

    }


}

}

//===========================================================================

AliFlowEventCuts* DefinecutsEvent(Int_t icentr){
    AliFlowEventCuts* cutsEvent = new AliFlowEventCuts(Form("eventcuts_%d",icentr));
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





