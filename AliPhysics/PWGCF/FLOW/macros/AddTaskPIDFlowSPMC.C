
class AliAnalysisDataContainer;
class AliFlowTrackCuts;
class AliFlowTrackSimpleCuts;
class AliFlowEventCuts;
class AliFlowEventSimpleCuts;


void AddTaskPIDFlowSPMC(Int_t triggerSelectionString=AliVEvent::kMB,
                        Int_t uptoWhichHarmonics = 2, // 2 --> v2 only, 3 --> v2 and v3, and so on
                        Float_t etamin=-0.8,
                        Float_t etamax=0.8,
                        Float_t EtaGap=0.2,
                        TString fileNameBase="AnalysisResults",
                        TString Qvector ="Qa",
                        Int_t charge=0,
                        Int_t ncentralitymin = 0,
                        Int_t ncentralitymax = 50,
                        Bool_t isPID = kTRUE,
                        Bool_t useAfterBurner=kFALSE,
                        AliPID::EParticleType particleType=AliPID::kPion) {
    
    
    TF1 *gV2Param = 0x0;
    Bool_t doQA=kTRUE;
    
    // AFTERBURNER
    Double_t v1=0.0;
    Double_t v2=0.0;
    Double_t v3=0.0;
    Double_t v4=0.0;
    Int_t numberOfTrackClones=0; //non-flow
    
    //Define the range for eta subevents (for SP method) with TPC
    Double_t minA = etamin;//
    Double_t maxA = -0.5*EtaGap;//
    Double_t minB = +0.5*EtaGap;//
    Double_t maxB = etamax;//
    
    
    int centrMin[6] = {0,1,10,20,30,40};
    int centrMax[6] = {1,2,20,30,40,50};
    
    for(int i=0;i<6;i++){
        if(ncentralitymin == centrMin[i]) const int ncentrminlim = i;
        if(ncentralitymax == centrMax[i]) const int ncentrmaxlim = i;
    }
    
    //---------Data selection---------- ESD only!!!
    //kMC, kGlobal, kESD_TPConly, kESD_SPDtracklet
    AliFlowTrackCuts::trackParameterType rptype = AliFlowTrackCuts::kMC;
    AliFlowTrackCuts::trackParameterType poitype = AliFlowTrackCuts::kMC;
    
    //---------Parameter mixing--------
    //kPure - no mixing, kTrackWithMCkine, kTrackWithMCPID, kTrackWithMCpt
    AliFlowTrackCuts::trackParameterMix rpmix = AliFlowTrackCuts::kPure;
    AliFlowTrackCuts::trackParameterMix poimix = AliFlowTrackCuts::kPure;
    
    const char* rptypestr = AliFlowTrackCuts::GetParamTypeName(rptype);  //ESD
    const char* poitypestr = AliFlowTrackCuts::GetParamTypeName(poitype); //ESD
    
    
    //===========================================================================
    // EVENTS CUTS:
    const int ncentr = ncentrmaxlim - ncentrminlim + 1;
    
    const int nharmonics = uptoWhichHarmonics-1;
    AliFlowEventCuts* cutsEvent[ncentr];
    AliFlowTrackCuts* cutsRP[ncentr];
    AliFlowTrackCuts* cutsPOI[ncentr];
    TString outputSlotName[ncentr][nharmonics];
    TString suffixName[ncentr];
    
    for(int icentr=0;icentr<ncentr;icentr++){
        cutsEvent[icentr] = new AliFlowEventCuts(Form("eventcuts_%d",icentr));
        
        cutsEvent[icentr]->SetQA(doQA);
        
        cutsRP[icentr] = new AliFlowTrackCuts(Form("RP_%d",icentr));
        cutsRP[icentr]->SetParamType(rptype);
        cutsRP[icentr]->SetParamMix(rpmix);
        cutsRP[icentr]->SetPtRange(0.2,5.);
        cutsRP[icentr]->SetEtaRange(etamin,etamax);
        cutsRP[icentr]->SetQA(doQA);
        //POIs for SP and QC method
        //===========================================================================
        AliFlowTrackCuts  *SP_POI[ncentr];
        //half window for POIs
        //=======================SP POI Cuts
        SP_POI[icentr] = new AliFlowTrackCuts("POI");

        SP_POI[icentr]->SetParamType(poitype);
        SP_POI[icentr]->SetParamMix(poimix);
        SP_POI[icentr]->SetPtRange(0.2,10.);
        if(isPID){
            if(particleType==AliPID::kPion) SP_POI[icentr]->SetMCPID(211);
            if(particleType==AliPID::kKaon) SP_POI[icentr]->SetMCPID(321);
            if(particleType==AliPID::kProton) SP_POI[icentr]->SetMCPID(2212);
        }
        
        if(charge != 0) SP_POI[icentr]->SetCharge(charge);
        
        if(Qvector=="Qa"){
            SP_POI[icentr]->SetEtaRange( +0.5*EtaGap, etamax );
            printf(" > NOTE: Using half TPC (Qb) as POI selection u < \n");
        }
        if(Qvector=="Qb"){
            SP_POI[icentr]->SetEtaRange( etamin,-0.5*EtaGap );
            printf(" > NOTE: Using half TPC (Qa) as POI selection u < \n");
        }
        SP_POI[icentr]->SetQA(doQA);
        //=====================================================================
        
        if(Qvector=="Qa") suffixName[icentr] = "Qa";
        if(Qvector=="Qb") suffixName[icentr] = "Qb";
        suffixName[icentr] += "_flow_";
        suffixName[icentr] += Form("%i_", centrMin[icentr+ncentrminlim]);
        suffixName[icentr] += Form("%i_", centrMax[icentr+ncentrminlim]);
        suffixName[icentr] += Form("%.f_", EtaGap*10);
        
        if(isPID){
            suffixName[icentr]+="_";
            suffixName[icentr]+=AliPID::ParticleName(particleType);
        }
        else{
            suffixName[icentr]+="AllCharged";
        }
        if (charge<0) suffixName[icentr]+="-";
        if (charge>0) suffixName[icentr]+="+";
        
        
        for(int harmonic=2;harmonic<uptoWhichHarmonics+1;harmonic++){  //for v2,v3,v4 and v5
            outputSlotName[icentr][harmonic-2] = "";
            if(isPID){
                outputSlotName[icentr][harmonic-2]+="_";
                outputSlotName[icentr][harmonic-2]+=AliPID::ParticleName(particleType);//particleType
            }else{
                outputSlotName[icentr][harmonic-2]+="AllCharged";
            }
            outputSlotName[icentr][harmonic-2]+=Form("_v%i_",harmonic);
            outputSlotName[icentr][harmonic-2]+=Form("%i-",centrMin[icentr+ncentrminlim]);
            outputSlotName[icentr][harmonic-2]+=Form("%i_",centrMax[icentr+ncentrminlim]);

            if (charge<0) outputSlotName[icentr][harmonic-2]+="-";
            if (charge>0) outputSlotName[icentr][harmonic-2]+="+";
        }
    }//loop over centrality classes
    
    
    TString fileName(fileNameBase);
    fileName.Append(".root");
    
    
    
    //====================================FLOWPACKAGE TASKS=========================//
    AliAnalysisDataContainer *cinput1[ncentr];
    AliAnalysisDataContainer *coutputFE[ncentr];
    AliAnalysisDataContainer* coutputFEQA[ncentr];
    AliAnalysisTaskFlowEvent *taskFE[ncentr];
    
    AliAnalysisDataContainer *flowEvent[ncentr][nharmonics];
    AliAnalysisTaskFilterFE *tskFilter[ncentr][nharmonics];
    
    AliAnalysisDataContainer *coutputSP[ncentr][nharmonics];
    AliAnalysisTaskScalarProduct *taskSP[ncentr][nharmonics];
    
    TString outputQA[ncentr];
    TString myNameSP[ncentr][nharmonics];
    TString slot[ncentr][nharmonics];
    
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
        if(useAfterBurner)
        {
            taskFE[icentr] = new AliAnalysisTaskFlowEvent(Form("TaskFlowEvent_%s",suffixName[icentr].Data()),"",doQA,1);
            //taskFE->SetFlow(v1,v2,v3,v4);
            taskFE[icentr]->SetPtDifferentialV2(gV2Param);
            taskFE[icentr]->SetNonFlowNumberOfTrackClones(numberOfTrackClones);
            taskFE[icentr]->SetAfterburnerOn();
        }else {taskFE[icentr] = new AliAnalysisTaskFlowEvent(Form("TaskFlowEvent_%s",suffixName[icentr].Data()),"",doQA); }
        
        taskFE[icentr]->SetSubeventEtaRange(minA, maxA, minB, maxB);
        mgr->AddTask(taskFE[icentr]);
        
        // Pass cuts for RPs and POIs to the task:
        taskFE[icentr]->SetCutsEvent(cutsEvent[icentr]);
        taskFE[icentr]->SetCutsRP(cutsRP[icentr]);
        taskFE[icentr]->SetCutsPOI(SP_POI[icentr]);
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
        TString Species = "";
        if(isPID) Species += AliPID::ParticleName(particleType);
        else      Species += "Allcharged";
        
        
        for(int harm=2;harm<uptoWhichHarmonics+1;harm++){
            myNameSP[icentr][harm-2] = "SP_";
            myNameSP[icentr][harm-2] += Qvector;
            myNameSP[icentr][harm-2] += Form("_v%i_%s_%.f",harm,outputSlotName[icentr][harm-2].Data(),EtaGap*10);
            
            flowEvent[icentr][harm-2] = mgr->CreateContainer( Form("Filter_%s", myNameSP[icentr][harm-2].Data()),
                                                             AliFlowEventSimple::Class(),
                                                             AliAnalysisManager::kExchangeContainer );
            
            tskFilter[icentr][harm-2] = new AliAnalysisTaskFilterFE( Form("TaskFilter_%s",myNameSP[icentr][harm-2].Data()),cutsRP[icentr], NULL);//SP_POI[icentr]
            tskFilter[icentr][harm-2]->SetSubeventEtaRange(minA, maxA, minB, maxB);
            mgr->AddTask(tskFilter[icentr][harm-2]);
            mgr->ConnectInput( tskFilter[icentr][harm-2],0,coutputFE[icentr]);
            mgr->ConnectOutput(tskFilter[icentr][harm-2],1,flowEvent[icentr][harm-2]);
            
            
            taskSP[icentr][harm-2] = new AliAnalysisTaskScalarProduct(Form("TaskScalarProduct_%s",outputSlotName[icentr][harm-2].Data()),kFALSE);
            taskSP[icentr][harm-2]->SetHarmonic(harm);
            taskSP[icentr][harm-2]->SetRelDiffMsub(1.0);
            taskSP[icentr][harm-2]->SetTotalQvector(Qvector);
            taskSP[icentr][harm-2]->SetApplyCorrectionForNUA(kTRUE);
            
            TString outputSP = fileName;
            outputSP += ":outputSPanalysis";
            outputSP+= rptypestr;
            slot[icentr][harm-2] = "SP_";
            slot[icentr][harm-2] += outputSlotName[icentr][harm-2];
            slot[icentr][harm-2] += "_";
            slot[icentr][harm-2] += Qvector;
            coutputSP[icentr][harm-2] = mgr->CreateContainer(Form("%s_%.f",slot[icentr][harm-2].Data(),EtaGap*10),
                                                             TList::Class(),AliAnalysisManager::kOutputContainer,outputSP);
            mgr->AddTask(taskSP[icentr][harm-2]);
            mgr->ConnectInput(taskSP[icentr][harm-2],0,flowEvent[icentr][harm-2]);
            mgr->ConnectInput(taskSP[icentr][harm-2],0,coutputFE[icentr]);
            mgr->ConnectOutput(taskSP[icentr][harm-2],1,coutputSP[icentr][harm-2]);
        }
        
        
        if (taskFE[icentr]->GetQAOn()) {
            outputQA[icentr] = fileName;
            outputQA[icentr] += ":QA";
            coutputFEQA[icentr] = mgr->CreateContainer(Form("QA_%s",suffixName[icentr].Data()), TList::Class(),AliAnalysisManager::kOutputContainer,outputQA[icentr]);
            mgr->ConnectOutput(taskFE[icentr],2,coutputFEQA[icentr]);
        }
        
        
    }
}
