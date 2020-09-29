AliAnalysisTaskEmcalHfeTagging* AddTaskEmcalHfeTagging(const char * njetsBase,
                                                       const char * njetsUS,
                                                       const char * njetsTrue,
                                                       const char * njetsPartLevel,
                                                       const Double_t R,
                                                       const char * nrhoBase,
                                                       const char * ntracks,
                                                       const char * ntracksUS,
                                                       const char *ntracksPartLevel,
                                                       const char * nclusters,
                                                       const char * ntracksTrue,
                                                       const char *type,
                                                       const char *CentEst,
                                                       Int_t       pSel,
                                                       TString     trigClass      = "",
                                                       TString     kEmcalTriggers = "",
                                                       TString     tag            = "",
                                                       AliAnalysisTaskEmcalHfeTagging::JetShapeType jetShapeType = AliAnalysisTaskEmcalHfeTagging::kData,
                                                       AliAnalysisTaskEmcalHfeTagging::JetShapeSub jetShapeSub = AliAnalysisTaskEmcalHfeTagging::kNoSub,
                                                       AliAnalysisTaskEmcalHfeTagging::JetSelectionType jetSelection = AliAnalysisTaskEmcalHfeTagging::kInclusive,
                                                       Float_t minpTHTrigger =0.,
                                                       Float_t maxpTHTrigger =0.,
                                                       Int_t derivSubtrOrder = 0,
                                                       Int_t MCweight = 0,
                                                       Double_t AssPtCut = 0.1,
                                                       Int_t ITSncut = 3,
                                                       Int_t AssTPCnCut = 60,
                                                       Int_t TPCnCut = 100,
                                                       Double_t SigmaTOFcut = 3.,
                                                       Double_t SigmaTPCcutLowPt = -1.,
                                                       Double_t SigmaTPCcutHighPt = -1.5,
                                                       Double_t SigmTPCcutExcElec = 3.5,
                                                       Double_t DcaXYcut = 1.,
                                                       Double_t DcaZcut = 2.,
                                                       Double_t IMcut = 0.1,
                                                       Double_t EtaCut = 0.7,
                                                       Double_t MinEoPcut = 0.9,
                                                       Double_t MaxEoPcut = 1.3,
                                                       Double_t M20cut = 0.35,
                                                       Double_t MinPtTPC = 0.5,
                                                       Double_t MaxPtTPC = 4.,
                                                       Double_t MinPtEMCal = 4.,
                                                       Double_t MaxPtEMCal = 25.
                                                       
                                                       ) {
    
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr)
    {
        Error("AddTaskEmcalHfeTagging","No analysis manager found.");
        return 0;
    }
    Bool_t ismc=kFALSE;
    ismc = (mgr->GetMCtruthEventHandler())?kTRUE:kFALSE;
    
    // Check the analysis type using the event handlers connected to the analysis manager.
    //==============================================================================
    if (!mgr->GetInputEventHandler())
    {
        ::Error("AddTaskEmcalHfeTagging", "This task requires an input event handler");
        return NULL;
    }
    
    TString wagonName1;
    TString wagonName2;
    
    wagonName1 = Form("JetHfeTaggings_%s_TC%s%sList_All",njetsBase,trigClass.Data(),tag.Data());
    wagonName2 = Form("JetHfeTaggings_%s_TC%s%sTree_All",njetsBase,trigClass.Data(),tag.Data());
    
    //Configure jet tagger task
    AliAnalysisTaskEmcalHfeTagging *task = new AliAnalysisTaskEmcalHfeTagging(wagonName1.Data());
    
    
    //task->SetNCentBins(4);
    task->SetJetShapeType(jetShapeType);
    task->SetJetShapeSub(jetShapeSub);
    task->SetJetSelection(jetSelection);
    task->SetDerivativeSubtractionOrder(derivSubtrOrder);
    task->SetMCweight(MCweight);
    task->SetAssPtCut(AssPtCut);
    task->SetITSncut(ITSncut);
    task->SetAssTPCnCut(AssTPCnCut);
    task->SetTPCnCut(TPCnCut);
    task->SetSigmaTOFcut(SigmaTOFcut);
    task->SetSigmaTPCcutLowPt(SigmaTPCcutLowPt);
    task->SetSigmaTPCcutHighPt(SigmaTPCcutHighPt);
    task->SetSigmTPCcutExcElec(SigmTPCcutExcElec);
    task->SetDcaXYcut(DcaXYcut);
    task->SetDcaZcut(DcaZcut);
    task->SetIMcut(IMcut);
    task->SetEtaCut(EtaCut);
    task->SetMinEoPcut(MinEoPcut);
    task->SetMaxEoPcut(MaxEoPcut);
    task->SetM20cut(M20cut);
    task->SetMinPtTPC(MinPtTPC);
    task->SetMaxPtTPC(MaxPtTPC);
    task->SetMinPtEMCal(MinPtEMCal);
    task->SetMaxPtEMCal(MaxPtEMCal);
    
    
    
    if (jetSelection == AliAnalysisTaskEmcalHfeTagging::kRecoil) task->SetPtTriggerSelections(minpTHTrigger, maxpTHTrigger);
    
    TString thename(njetsBase);
    //if(thename.Contains("Sub")) task->SetIsConstSub(kTRUE);
    //task->SetVzRange(-10.,10.);
    
    AliParticleContainer *trackCont;// = task->AddTrackContainer(ntracks);
    
    if ((jetShapeSub==AliAnalysisTaskEmcalHfeTagging::kConstSub) && ((jetShapeType==AliAnalysisTaskEmcalHfeTagging::kData) || (jetShapeType==AliAnalysisTaskEmcalHfeTagging::kDetEmbPartPythia) || (jetShapeType==AliAnalysisTaskEmcalHfeTagging::kPythiaDef))){
        trackCont = task->AddParticleContainer(ntracks);}
    else trackCont = task->AddTrackContainer(ntracks);
    
    
    //Printf("tracks() = %s, trackCont =%p", ntracks, trackCont);
    AliParticleContainer *trackContUS  = task->AddTrackContainer(ntracksUS);
    //Printf("tracksUS() = %s", ntracksUS);
    AliParticleContainer *trackContTrue = task->AddMCParticleContainer(ntracksTrue);
    //Printf("ntracksTrue() = %s, trackContTrue=%p ", ntracksTrue, trackContTrue);
    
    AliParticleContainer *trackContPartLevel=0;
    
    if ((jetShapeSub==AliAnalysisTaskEmcalHfeTagging::kConstSub) && ((jetShapeType==AliAnalysisTaskEmcalHfeTagging::kMCTrue) || (jetShapeType==AliAnalysisTaskEmcalHfeTagging::kPythiaDef))){
        trackContPartLevel = task->AddParticleContainer(ntracksPartLevel);
    }
    else trackContPartLevel = task->AddMCParticleContainer(ntracksPartLevel);
    
    //Printf("ntracksPartLevel() = %s, trackContPartLevel=%p ", ntracksPartLevel, trackContPartLevel);
    
    
    AliClusterContainer *clusterCont = task->AddClusterContainer(nclusters);
    
    AliJetContainer *jetContBase=0x0;
    AliJetContainer *jetContUS=0x0;
    AliJetContainer *jetContTrue=0x0;
    AliJetContainer *jetContPart=0x0;
    TString strType(type);
    
    if ((jetShapeType==AliAnalysisTaskEmcalHfeTagging::kMCTrue || (jetShapeType==AliAnalysisTaskEmcalHfeTagging::kGenOnTheFly))) {
        jetContBase = task->AddJetContainer(njetsBase,strType,R);
        if(jetContBase) {
            jetContBase->SetRhoName(nrhoBase);
            jetContBase->ConnectParticleContainer(trackContPartLevel);
            jetContBase->ConnectClusterContainer(clusterCont);
            jetContBase->SetPercAreaCut(0.6);
        }
    }
    
    if (jetShapeType==AliAnalysisTaskEmcalHfeTagging::kData){
        jetContBase = task->AddJetContainer(njetsBase,strType,R);
        if(jetContBase) {
            jetContBase->SetRhoName(nrhoBase);
            jetContBase->ConnectParticleContainer(trackCont);
            jetContBase->ConnectClusterContainer(clusterCont);
            jetContBase->SetPercAreaCut(0.6);
            if(jetShapeSub==AliAnalysisTaskEmcalHfeTagging::kConstSub) jetContBase->SetAreaEmcCut(-2);
        }
    }
    
    
    if (jetShapeType==AliAnalysisTaskEmcalHfeTagging::kDetEmbPartPythia){
        jetContBase = task->AddJetContainer(njetsBase,strType,R);
        if(jetContBase) {
            jetContBase->SetRhoName(nrhoBase);
            jetContBase->ConnectParticleContainer(trackCont);
            jetContBase->ConnectClusterContainer(clusterCont);
            jetContBase->SetPercAreaCut(0.6);
            
            if(jetShapeSub==AliAnalysisTaskEmcalHfeTagging::kConstSub) jetContBase->SetAreaEmcCut(-2);
        }
        
        jetContTrue = task->AddJetContainer(njetsTrue,strType,R);
        if(jetContTrue) {
            jetContTrue->SetRhoName(nrhoBase);
            jetContTrue->ConnectParticleContainer(trackContTrue);
            jetContTrue->SetPercAreaCut(0.6);
            
        }
        
        if(jetShapeSub==AliAnalysisTaskEmcalHfeTagging::kConstSub){
            jetContUS=task->AddJetContainer(njetsUS,strType,R);
            if(jetContUS) {
                jetContUS->SetRhoName(nrhoBase);
                jetContUS->ConnectParticleContainer(trackContUS);
                jetContUS->SetPercAreaCut(0.6);
                
            }
        }
        
        jetContPart = task->AddJetContainer(njetsPartLevel,strType,R);
        if(jetContPart) {
            jetContPart->SetRhoName(nrhoBase);
            jetContPart->ConnectParticleContainer(trackContPartLevel);
            jetContPart->SetPercAreaCut(0.6);
            
        }
    }
    
    if (jetShapeType==AliAnalysisTaskEmcalHfeTagging::kPythiaDef){
        
        jetContBase = task->AddJetContainer(njetsBase,strType,R);
        if(jetContBase) {
            jetContBase->ConnectParticleContainer(trackCont);
            jetContBase->ConnectClusterContainer(clusterCont);
            jetContBase->SetPercAreaCut(0.6);
        }
        
        jetContTrue = task->AddJetContainer(njetsTrue,strType,R);
        if(jetContTrue) {
            jetContTrue->SetRhoName(nrhoBase);
            jetContTrue->ConnectParticleContainer(trackContTrue);
            jetContTrue->SetPercAreaCut(0.6);
            
        }
        
        if(jetShapeSub==AliAnalysisTaskEmcalHfeTagging::kConstSub){
            jetContUS=task->AddJetContainer(njetsUS,strType,R);
            if(jetContUS) {
                jetContUS->SetRhoName(nrhoBase);
                jetContUS->ConnectParticleContainer(trackContUS);
                jetContUS->SetPercAreaCut(0.6);
                
            }
        }
        
        jetContPart = task->AddJetContainer(njetsPartLevel,strType,R);
        if(jetContPart) {
            jetContPart->SetRhoName(nrhoBase);
            jetContPart->ConnectParticleContainer(trackContPartLevel);
            jetContPart->SetPercAreaCut(0.6);
        }
        
    }
    
    task->SetCaloTriggerPatchInfoName(kEmcalTriggers.Data());
    task->SetCentralityEstimator(CentEst);
    task->SelectCollisionCandidates(pSel);
    task->SetUseAliAnaUtils(kFALSE);
    
    mgr->AddTask(task);
    
    //Connnect input
    mgr->ConnectInput (task, 0, mgr->GetCommonInputContainer() );
    
    //Connect output
    TString contName1(wagonName1);
    TString contName2(wagonName2);
    
    if (jetShapeType == AliAnalysisTaskEmcalHfeTagging::kMCTrue){
        contName1 += "_MCTrue";
        contName2 += "_MCTrue";
    }
    
    if (jetShapeType == AliAnalysisTaskEmcalHfeTagging::kData){
        contName1 += "_Data";
        contName2 += "_Data";
    }
    
    if (jetShapeType == AliAnalysisTaskEmcalHfeTagging::kPythiaDef){
        contName1 +="_PythiaDef";
        contName2 +="_PythiaDef";
    }
    
    
    if (jetShapeSub == AliAnalysisTaskEmcalHfeTagging::kNoSub){
        contName1 += "_NoSub";
        contName2 += "_NoSub";
    }
    
    
    if (jetShapeSub == AliAnalysisTaskEmcalHfeTagging::kConstSub){
        contName1 += "_ConstSub";
        contName2 += "_ConstSub";
    }
    
    if (jetShapeSub == AliAnalysisTaskEmcalHfeTagging::kDerivSub){
        contName1 += "_DerivSub";
        contName2 += "_DerivSub";
    }
    
    if (jetSelection == AliAnalysisTaskEmcalHfeTagging::kInclusive){
        contName1 += "_Incl";
        contName2 += "_Incl";
    }
    
    
    if (jetSelection == AliAnalysisTaskEmcalHfeTagging::kRecoil) {
        TString recoilTriggerString = Form("_Recoil_%.0f_%0.f", minpTHTrigger, maxpTHTrigger);
        contName1 += recoilTriggerString;
        contName2 += recoilTriggerString;
    }
    
    
    TString outputfile = Form("%s",AliAnalysisManager::GetCommonFileName());
    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contName1.Data(), TList::Class(),AliAnalysisManager::kOutputContainer,outputfile);
    mgr->ConnectOutput(task,1,coutput1);
    
    AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(contName2.Data(), TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile);
    mgr->ConnectOutput(task,2,coutput2);
    
    return task;
    
}

