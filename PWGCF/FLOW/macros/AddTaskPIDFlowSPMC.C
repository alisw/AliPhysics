
class AliAnalysisDataContainer;
class AliFlowTrackCuts;
class AliFlowTrackSimpleCuts;
class AliFlowEventCuts;
class AliFlowEventSimpleCuts;


void AddTaskPIDFlowSPMC(Int_t triggerSelectionString=AliVEvent::kMB,
                        Int_t uptoWhichHarmonics = 2, // 2 --> v2 only, 3 --> v2 and v3, and so on
                        Float_t PurityLevel=0.8,
                        Float_t etamin=-0.8,
                        Float_t etamax=0.8,
                        Float_t EtaGap=0.2,
                        TString fileNameBase="AnalysisResults",
                        TString uniqueStr="Pion_02",
                        TString Qvector ="Qa",
                        Int_t AODfilterBit = 1,
                        Int_t charge=0,
                        Int_t MinTPCdedx = 10,
                        Int_t ncentralitymin = 0,
                        Int_t ncentralitymax = 50,
                        Int_t maxITSCls = 7,
                        Int_t maxChi2ITSCls = 37,
                        Bool_t isPID = kTRUE,
                        Bool_t isVZERO = kFALSE, // use vzero sp method
                        Bool_t isData = kTRUE,
                        Bool_t is2011 = kTRUE,
                        Bool_t isAOD = kTRUE,
                        Bool_t UsePurityPIDmethod = kFALSE,
                        Bool_t useAfterBurner=kTRUE,
                        Bool_t isCheck=kFALSE, //Purity PID pt differential vn extent
                        AliPID::EParticleType particleType=AliPID::kPion,
                        AliFlowTrackCuts::PIDsource sourcePID=AliFlowTrackCuts::kTOFbayesian) {
    
    
    TF1 *gV2Param = 0x0;
    Bool_t doQA=kTRUE;
    
    // AFTERBURNER
    Double_t v1=0.0;
    Double_t v2=0.0;
    Double_t v3=0.0;
    Double_t v4=0.0;
    Int_t numberOfTrackClones=0; //non-flow
    
    // Define a range of the detector to exclude
    Bool_t ExcludeRegion = kFALSE;
    Double_t excludeEtaMin = -0.;
    Double_t excludeEtaMax = 0.;
    Double_t excludePhiMin = 0.;
    Double_t excludePhiMax = 0.;
    
    //Define the range for eta subevents (for SP method) with TPC
    Double_t minA = etamin;//
    Double_t maxA = -0.5*EtaGap;//
    Double_t minB = +0.5*EtaGap;//
    Double_t maxB = etamax;//
    
    
    int centrMin[16] = {0,1,2,3,4,5,6,7,8,9,10,20,30,40,60,70};
    int centrMax[16] = {1,2,3,4,5,6,7,8,9,10,20,30,40,50,70,80};
    
    for(int i=0;i<16;i++){
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
        if(isData){
            cutsEvent[icentr]->SetLHC11h(is2011);
            cutsEvent[icentr]->SetCentralityPercentileRange(centrMin[icentr+ncentrminlim],centrMax[icentr+ncentrminlim]);
            cutsEvent[icentr]->SetCentralityPercentileMethod(AliFlowEventCuts::kV0);
            cutsEvent[icentr]->SetPrimaryVertexZrange(-10.,10.);
            cutsEvent[icentr]->SetCutTPCmultiplicityOutliers();
        }
        cutsEvent[icentr]->SetQA(doQA);
        
        if(isData){
            // RP TRACK CUTS:
            if(!isVZERO){
                cutsRP[icentr] = new AliFlowTrackCuts(Form("RP_%d",icentr));
                cutsRP[icentr]->SetPtRange(0.2,5.);
                cutsRP[icentr]->SetEtaRange(etamin,etamax);
                //cutsRP[icentr]->SetMinNClustersTPC(70);
                //cutsRP[icentr]->SetMinChi2PerClusterTPC(0.1);
                //cutsRP[icentr]->SetMaxChi2PerClusterTPC(4.0);
                //cutsRP[icentr]->SetMaxDCAToVertexXY(2.4);
                //cutsRP[icentr]->SetMaxDCAToVertexZ(3.0);
                //cutsRP[icentr]->SetAcceptKinkDaughters(kFALSE);
                cutsRP[icentr]->SetMinimalTPCdedx(MinTPCdedx);
                cutsRP[icentr]->SetAODfilterBit(AODfilterBit);
            }
            
            if(isVZERO) { // use vzero sub analysis
                if(!is2011) cutsRP[icentr] = AliFlowTrackCuts::GetStandardVZEROOnlyTrackCuts2010(); // select vzero tracks
                if(is2011)  cutsRP[icentr] = AliFlowTrackCuts::GetStandardVZEROOnlyTrackCuts2011(); // select vzero tracks
                
                if(!cutsRP[icentr]) {
                    cout << " Fatal error: no RP cuts found! " << endl;
                    return 0x0;
                }
            }//vzero is not a tracking device it is just a scintillator. so pt range or DCAtoVertex are not set here.
        }
        
        if(!isData){
            cutsRP[icentr] = new AliFlowTrackCuts(Form("RP_%d",icentr));
            cutsRP[icentr]->SetParamType(rptype);
            cutsRP[icentr]->SetParamMix(rpmix);
            cutsRP[icentr]->SetPtRange(0.2,5.);
            cutsRP[icentr]->SetEtaRange(etamin,etamax);
        }
        cutsRP[icentr]->SetQA(doQA);
        //POIs for SP and QC method
        //===========================================================================
        AliFlowTrackCuts  *SP_POI[ncentr];
        //half window for POIs
        //=======================SP POI Cuts
        SP_POI[icentr] = DefinePOIcuts();
        if(isData){
            if(!is2011) SP_POI[icentr]->GetBayesianResponse()->ForceOldDedx(); // for 2010 data to use old TPC PID Response instead of the official one
            if(!isAOD){
                SP_POI[icentr]->SetMaxSharedITSCluster(maxITSCls);
                SP_POI[icentr]->SetMaxChi2perITSCluster(maxChi2ITSCls);
                SP_POI[icentr]->SetMinNClustersITS(2);
                SP_POI[icentr]->SetRequireITSRefit(kTRUE);
                SP_POI[icentr]->SetRequireTPCRefit(kTRUE);
                SP_POI[icentr]->SetMaxDCAToVertexXY(0.3);
                SP_POI[icentr]->SetMaxDCAToVertexZ(0.3);
                SP_POI[icentr]->SetAcceptKinkDaughters(kFALSE);
                SP_POI[icentr]->SetMinimalTPCdedx(10.);
            }
            //SP_POI[icentr]->SetParamMix(poimix);
            SP_POI[icentr]->SetPtRange(0.2,6.);//
            SP_POI[icentr]->SetMinNClustersTPC(70);
            SP_POI[icentr]->SetMinChi2PerClusterTPC(0.1);
            SP_POI[icentr]->SetMaxChi2PerClusterTPC(4.0);
            
            //SP_POI->SetRequireITSRefit(kTRUE);
            //SP_POI->SetRequireTPCRefit(kTRUE);
            //SP_POI->SetMinNClustersITS(2);
            //SP_POI->SetMaxChi2PerClusterITS(1.e+09);
            //SP_POI[icentr]->SetMaxDCAToVertexXY(2.4);
            //SP_POI[icentr]->SetMaxDCAToVertexZ(3.0);
            //SP_POI->SetDCAToVertex2D(kTRUE);
            //SP_POI->SetMaxNsigmaToVertex(1.e+10);
            //SP_POI->SetRequireSigmaToVertex(kFALSE);
            //SP_POI[icentr]->SetAcceptKinkDaughters(kFALSE);
            //SP_POI->SetAllowTOFmismatch(kFALSE);
            SP_POI[icentr]->SetRequireStrictTOFTPCagreement(kTRUE);
            SP_POI[icentr]->SetMinimalTPCdedx(MinTPCdedx);
            if(isAOD) SP_POI[icentr]->SetAODfilterBit(AODfilterBit);
            SP_POI[icentr]->SetPriors((centrMin[icentr+ncentrminlim]+centrMax[icentr+ncentrminlim])*0.5);
            
            if(isPID){
                SP_POI[icentr]->SetPID(particleType, sourcePID);//particleType, sourcePID
                
                if(UsePurityPIDmethod){
                    SP_POI[icentr]->SetCentralityPercentile(centrMin[icentr+ncentrminlim],centrMax[icentr+ncentrminlim]);
                    SP_POI[icentr]->SetTPCTOFNsigmaPIDPurityFunctions(PurityLevel);
                }
            }
            
        }
        if(!isData){
            SP_POI[icentr]->SetParamType(poitype);
            SP_POI[icentr]->SetParamMix(poimix);
            SP_POI[icentr]->SetPtRange(0.2,10.);
            if(isPID){
                if(particleType==AliPID::kPion) SP_POI[icentr]->SetMCPID(211);
                if(particleType==AliPID::kKaon) SP_POI[icentr]->SetMCPID(321);
                if(particleType==AliPID::kProton) SP_POI[icentr]->SetMCPID(2212);
            }
        }
        
        if(charge != 0) SP_POI[icentr]->SetCharge(charge);
        
        if(isCheck){
            SP_POI[icentr]->SetEtaRange( etamin, etamax );
            printf(" > NOTE: Using full TPC as POI selection u < \n");
        }else{
            if(!isVZERO && Qvector=="Qa"){
                SP_POI[icentr]->SetEtaRange( +0.5*EtaGap, etamax );
                printf(" > NOTE: Using half TPC (Qb) as POI selection u < \n");
            }
            if(!isVZERO && Qvector=="Qb"){
                SP_POI[icentr]->SetEtaRange( etamin,-0.5*EtaGap );
                printf(" > NOTE: Using half TPC (Qa) as POI selection u < \n");
            }
            if(isVZERO){
                SP_POI[icentr]->SetEtaRange( etamin,etamax );
                printf(" > NOTE: Using full TPC as POI selection u < \n");
            }
        }
        SP_POI[icentr]->SetQA(doQA);
        //=====================================================================
        
        if(!isVZERO && Qvector=="Qa") suffixName[icentr] = "Qa";
        if(!isVZERO && Qvector=="Qb") suffixName[icentr] = "Qb";
        if(isVZERO) suffixName[icentr] = "vzero";
        suffixName[icentr] += "_flow_";
        suffixName[icentr] += Form("%i_", centrMin[icentr+ncentrminlim]);
        suffixName[icentr] += Form("%i_", centrMax[icentr+ncentrminlim]);
        suffixName[icentr] += Form("%.f_", EtaGap*10);
        
        if(isPID){
            suffixName[icentr]+=AliFlowTrackCuts::PIDsourceName(sourcePID);
            suffixName[icentr]+="_";
            suffixName[icentr]+=uniqueStr;
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
            //outputSlotName[icentr][harmonic-2]+=cutsRP[icentr]->GetName();
            //outputSlotName[icentr][harmonic-2]+="_";
            //outputSlotName[icentr][harmonic-2]+=SP_POI[icentr]->GetName();
            outputSlotName[icentr][harmonic-2]+=Form("%i-",centrMin[icentr+ncentrminlim]);
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
        }
        else {taskFE[icentr] = new AliAnalysisTaskFlowEvent(Form("TaskFlowEvent_%s",suffixName[icentr].Data()),"",doQA); }
        if (ExcludeRegion) {
            taskFE[icentr]->DefineDeadZone(excludeEtaMin, excludeEtaMax, excludePhiMin, excludePhiMax);
        }
        
        if(isData) taskFE[icentr]->SelectCollisionCandidates(triggerSelectionString);
        
        //  if(taskFE[icentr]->SetVZEROSubEvents(EP3sub)) cout << " --> Setting up VZERO subevents method ... " << endl;
        if(!isVZERO) taskFE[icentr]->SetSubeventEtaRange(minA, maxA, minB, maxB);
        if(isVZERO)  taskFE[icentr]->SetSubeventEtaRange(-5,-1.5,+1.5,5);
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
            if(!isVZERO){
                tskFilter[icentr][harm-2]->SetSubeventEtaRange(minA, maxA, minB, maxB);
            }
            if(isVZERO) tskFilter[icentr][harm-2]->SetSubeventEtaRange(-5,-1.5,+1.5,5);
            mgr->AddTask(tskFilter[icentr][harm-2]);
            mgr->ConnectInput( tskFilter[icentr][harm-2],0,coutputFE[icentr]);
            mgr->ConnectOutput(tskFilter[icentr][harm-2],1,flowEvent[icentr][harm-2]);
            
            
            taskSP[icentr][harm-2] = new AliAnalysisTaskScalarProduct(Form("TaskScalarProduct_%s",outputSlotName[icentr][harm-2].Data()),kFALSE);
            taskSP[icentr][harm-2]->SetHarmonic(harm);
            if(isData)taskSP[icentr][harm-2]->SelectCollisionCandidates(triggerSelectionString);
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

//===========================================================================

AliFlowEventCuts* DefinecutsEvent(Int_t icentr){
    AliFlowEventCuts* cutsEvent = new AliFlowEventCuts(Form("eventcuts_%d",icentr));
    return cutsEvent;
}
AliFlowTrackCuts* DefineRPcuts(){
    AliFlowTrackCuts* cutsRP = new AliFlowTrackCuts"RP");
    return cutsRP;
}
AliFlowTrackCuts* DefinePOIcuts(){
    AliFlowTrackCuts* cutsPOI = new AliFlowTrackCuts("POI");
    return cutsPOI;
}




