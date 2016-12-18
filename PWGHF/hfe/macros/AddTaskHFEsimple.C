Bool_t ReadContaminationFunctionsBeauty(TString filename, TF1 **functions, double sigma){
    TFile *in = TFile::Open(Form("ALICE_PHYSICS/PWGHF/hfe/macros/configs/pPb/%s", filename.Data()));
    gROOT->cd();
    int isig = static_cast<int>(sigma * 100.);
    if (isig == -44) isig = -42;
    if (isig == 6) isig = 9;
    printf("Getting hadron background for the sigma cut: %d\n", isig);
    bool status = kTRUE;
    for(int icent = 0; icent < 12; icent++){
        functions[icent] = dynamic_cast<TF1 *>(in->Get(Form("hback_%d_%d", isig, icent)));
        if(functions[icent]) printf("Config for centrality class %d found\n", icent);
        else{
            printf("Config for the centrality class %d not found\n", icent);
            status = kFALSE;
        }
    }
    delete in;
    return status;
}

Bool_t ReadContaminationFunctions(TString filename, TF1 **functions, double sigma){
    TFile *in = TFile::Open(Form("ALICE_PHYSICS/PWGHF/hfe/macros/configs/pPb/%s", filename.Data()));
    gROOT->cd();
    int isig = static_cast<int>(sigma * 100.);
    printf("Getting hadron background for the sigma cut: %d\n", isig);
    bool status = kTRUE;
    for(int icent = 0; icent < 12; icent++){
        functions[icent] = dynamic_cast<TF1 *>(in->Get(Form("hback_%d_%d", isig, icent)));
        if(functions[icent]) printf("Config for centrality class %d found\n", icent);
        else{
            printf("Config for the centrality class %d not found\n", icent);
            status = kFALSE;
        }
    }
    delete in;
    return status;
}

AliAnalysisTask *AddTaskHFEsimple(Bool_t MCthere, Bool_t isAOD, TString str){

    TFile *f = TFile::Open(gSystem->ExpandPathName(str.Data()));
    if(!f || f->IsZombie()){
        printf("Could not read file %s\n",str.Data()); 
        return NULL ;
    }
    if(f->TestBit(TFile::kRecovered)){
        printf("File \"%s\" is corrupt!\n",str.Data());
    }
    gROOT->cd();
    TKey *k;
    TIter next(f->GetListOfKeys());
    while ((k = dynamic_cast<TKey *>(next()))){
        TString s(k->GetClassName());
        if(s.EqualTo("AliHFEparamBag")){
            AliHFEparamBag *b = dynamic_cast<AliHFEparamBag *>(k->ReadObj());
            //set MC and AOD settings (from train)
            b->useMC=MCthere;
            b->isAOD=isAOD;
            RegisterTask(b);
            //maybe printf to state that task was submitted?
            printf("*******************\nTask %s submitted!\n*******************\n",b->appendix.Data());
        } else{
            //report if something else is found
        }

    }
    if(!k){
        printf("No valid cut objects found!\n");
        f->Close(); delete f;
        return NULL;
    } 
    f->Close(); delete f;

    return NULL;
}

AliAnalysisTask *RegisterTask(AliHFEparamBag *abag){

    AliHFEparamBag *bag = new AliHFEparamBag(*abag);


    printf("Add macro appendix %s\n", bag->appendix.Data());

    if(bag->useMC&&!gROOT->GetListOfGlobalFunctions()->FindObject("ConfigWeightFactors")) 
        gROOT->LoadMacro("$TRAIN_ROOT/util/hfe/configs/ConfigWeightFactors.C");

    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

    Bool_t kAnalyseTaggedTracks = kTRUE;
    Bool_t kApplyPreselection = kFALSE;

    //***************************************//
    //        Setting up the HFE cuts        //
    //***************************************//

    AliHFEcuts *hfecuts = new AliHFEcuts(bag->appendix.Data(),"HFE cuts for pPb");
    //hfecuts->SetQAOn();
    hfecuts->CreateStandardCuts();
    hfecuts->SetMinNClustersTPC(bag->TPCcl);
    hfecuts->SetMinNClustersTPCPID(bag->TPCclPID);
    hfecuts->SetMinNClustersITS(bag->ITScl);
    hfecuts->SetMinRatioTPCclusters(0.6);
    hfecuts->SetTPCmodes(AliHFEextraCuts::kFound, AliHFEextraCuts::kFoundOverFindable);
    hfecuts->SetCutITSpixel(bag->itshitpixel);
    hfecuts->SetCheckITSLayerStatus(kFALSE);
    //if (bag->spdcheck>0){
    //hfecuts->SetCheckITSLayerStatus(kTRUE);
    //printf("\n\nWe check the live status of the pixels\n\n");
    //}
    hfecuts->SetEtaRange(bag->etami,bag->etama);
    if(bag->phimi >= 0. && bag->phima >= 0) hfecuts->SetPhiRange(bag->phimi,bag->phima);
    hfecuts->SetRejectKinkDaughters();
    hfecuts->SetAcceptKinkMothers();
    if(bag->isAOD) hfecuts->SetAODFilterBit(4);

    //if((iPixelAny==AliHFEextraCuts::kAny) || (iPixelAny==AliHFEextraCuts::kSecond))     

    hfecuts->SetMaxImpactParam(bag->DCAxy,bag->DCAz);
    hfecuts->SetUseMixedVertex(kTRUE);
    hfecuts->SetVertexRange(10.);
    // New pPb cuts (February 2013)
    hfecuts->SetUseCorrelationVertex();
    hfecuts->SetSPDVtxResolutionCut();
    hfecuts->SetpApileupCut();



    if(bag->isBeauty) hfecuts->SetProductionVertex(0,100,0,100);

    // TOF settings:
    Int_t usetof=0;
    Bool_t kTOFmis=kFALSE;
    if (bag->TOFs>0.){
        usetof = 1;
        printf("CONFIGURATION FILE: TOF is used \n");
        hfecuts->SetTOFPIDStep(kTRUE);
        printf("CONFIGURATION FILE: TOF PID step is requested !!!! \n");
        if (bag->TOFmis>0){
            kTOFmis = kTRUE;
            printf("CONFIGURATION FILE: TOF mismatch rejection is set ON \n");
        }
    }

    //***************************************//
    //        Setting up the task            //
    //***************************************//

    AliAnalysisTaskHFEtemplate *task = new AliAnalysisTaskHFEtemplate(Form("HFEtask%s",bag->appendix.Data()));
    printf("task %p\n", task);
    task->SetpPbAnalysis();
    if(!bag->isAOD) task->SetRemoveFirstEventInChunk();
    task->SetRemovePileUp(kFALSE);
    task->SetHFECuts(hfecuts);
    task->GetPIDQAManager()->SetHighResolutionHistos();
    task->SetRejectKinkMother(kFALSE);

    task->SetParams(bag);

    // Determine the centrality estimator
    task->SetCentralityEstimator("V0A");
    if (bag->icent == 2) task->SetCentralityEstimator("V0M");
    else if (bag->icent == 3) task->SetCentralityEstimator("CL1");
    else if (bag->icent == 4) task->SetCentralityEstimator("ZNA");

    //***************************************//
    //        Prepare preselection           //
    // This mimics the ESD->AOD filter in    //
    // case of the ESD analysis and selects  //
    // only tracks which will be selected in //
    // the AOD analysis with the given filter//
    // bit. Not to be applied for AODS.      //
    // For pPb the cuts used are (bit 4)     //
    // esdTrackCutsHG0 from file $ALICE_ROOT///
    // ANALYSIS/macros/AddTaskESDFilter.C    //
    //***************************************//

    if(kApplyPreselection){    
        AliESDtrackCuts* esdfilter = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE);
        esdfilter->SetMaxDCAToVertexXY(2.4);
        esdfilter->SetMaxDCAToVertexZ(3.2);
        esdfilter->SetDCAToVertex2D(kTRUE);

        task->SetHFECutsPreselect(esdfilter);
        printf("Put a preselection cut\n");
        task->SetFillNoCuts(kTRUE);
    }

    //***************************************//
    //          Variable manager             //
    //***************************************//
    // Define Variables
    Double_t ptbinning[36] = {0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.75, 2., 2.25, 2.5, 2.75, 3., 3.5, 4., 4.5, 5., 5.5, 6., 7., 8., 10., 12., 14., 16., 18., 20.};
    Double_t etabinning[17] = {-0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};

    Int_t sizept=(sizeof(ptbinning)/sizeof(double))-1;
    Int_t sizeeta=(sizeof(etabinning)/sizeof(double))-1;

    AliHFEvarManager *vm = task->GetVarManager();
    vm->AddVariable("pt", sizept, ptbinning);
    vm->AddVariable("eta", sizeeta, -0.8,0.8);
    //vm->AddVariable("phi",18, -0, 2*TMath::Pi());
    vm->AddVariable("phi",3, -0, 2*TMath::Pi());
    vm->AddVariable("charge");
    vm->AddVariable("source");
    vm->AddVariable("centrality");

    // For the moment, remove the part dedicated to the background subtraction.
    // It will be implemented in a different way, reading it from a root file.

    //***************************************//
    //          Configure the PID            //
    //***************************************//

    // Define PID
    AliHFEpid *pid = task->GetPID();
    if(bag->useMC) pid->SetHasMCData(kTRUE);

    if (usetof){
        pid->AddDetector("TOF", 0);
        pid->AddDetector("TPC", 1);
    } else {
        pid->AddDetector("TPC", 0);
    }

    // Configure TPC PID
    // do the identical thing in data and MC
    Double_t paramsTPCdEdxcutlow[12] ={0.0, 0.0, 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
    if(bag->tpcdEdxcutlow) memcpy(paramsTPCdEdxcutlow,bag->tpcdEdxcutlow,sizeof(paramsTPCdEdxcutlow));

    Double_t paramsTPCdEdxcuthigh[12] ={3.0, 3.0, 3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0};
    if(bag->tpcdEdxcuthigh) memcpy(paramsTPCdEdxcuthigh,bag->tpcdEdxcuthigh,sizeof(paramsTPCdEdxcuthigh));

    char *cutmodel;

    if(bag->useMC){ // constant (default) cut for MC
        cutmodel="pol0(0)";
        Double_t params[1];
        params[0]=paramsTPCdEdxcutlow[0];
        pid->ConfigureTPCdefaultCut(cutmodel, params,bag->tpcdEdxcuthigh[0]);
    } else { // correct for mean shift in data
        cutmodel="min(pol1(0),pol0(2))";
        Double_t params[3];
        //params[0]=-0.12; params[1]=0.14; params[2]=0.09;
        params[0]=-0.21 + paramsTPCdEdxcutlow[0];
        params[1]=0.14;
        params[2]=paramsTPCdEdxcutlow[0];
        pid->ConfigureTPCdefaultCut(cutmodel, params,bag->tpcdEdxcuthigh[0]);
    }

    // Configure TOF PID
    if (usetof){
        pid->ConfigureTOF(bag->TOFs);
        AliHFEpidTOF *tofpid = pid->GetDetPID(AliHFEpid::kTOFpid);
        if (kTOFmis){
            tofpid->SetRejectTOFmismatch();
        }
    }

    // Load hadron background
    if(!bag->useMC){
        Bool_t status = kTRUE;
        TF1 *hBackground[12];
        if(bag->isAOD==1) {
            if (usetof)  status = ReadContaminationFunctions("hadroncontamination_AOD139_TOFPID_pPb_eta06.root", hBackground, bag->tpcdEdxcutlow[0]);
            else { 
                if (bag->spdcheck == 0){
                    status = ReadContaminationFunctions("hadroncontamination_AOD139_noTOFPID_pPb_eta06.root", hBackground, bag->tpcdEdxcutlow[0]);
                } else if (bag->spdcheck == 1){ 
                    status = ReadContaminationFunctions("hadroncontamination_noTOFPID_pPb_AOD_eta06_TPCcut0_envelope_minsys.root", hBackground, bag->tpcdEdxcutlow[0]);	    
                } else if (bag->spdcheck == 2){ 
                    status = ReadContaminationFunctions("hadroncontamination_noTOFPID_pPb_AOD_eta06_TPCcut0_envelope_maxsys.root", hBackground, bag->tpcdEdxcutlow[0]);	    
                } else if (bag->spdcheck == 3){ 
                    status = ReadContaminationFunctions("hadroncontamination_AOD139_noTOFPID_pPb_eta06_polyFit.root", hBackground, bag->tpcdEdxcutlow[0]);	    
                }
            }
        }
        else if (bag->isBeauty==1) {
            if (bag->spdcheck == 0){
                printf("hadron cont standard beauty");

                if(bag->TOFs==0){
                    printf("hadron cont no TOF\n");
                    status = ReadContaminationFunctions("hadroncontamination_noTOF_pPb_Beauty_ESD_eta06.root", hBackground, bag->tpcdEdxcutlow[0]);
                }
                else{
                    printf("hadron cont with TOF\n");
                    status = ReadContaminationFunctionsBeauty("hadroncontamination_ESD_Beauty_TOFPID_pPb_eta06_iter2.root", hBackground, bag->tpcdEdxcutlow[0]);
                }
            } else if (bag->spdcheck == 1){
                printf("hadron cont min beauty");
                if(bag->TOFs==0){
                    printf("hadron cont no TOF\n");
                    status = ReadContaminationFunctions("hadroncontamination_noTOF_pPb_Beauty_ESD_eta06_envelopemin.root", hBackground, bag->tpcdEdxcutlow[0]);
                }
                status = ReadContaminationFunctionsBeauty("hadroncontamination_ESD_Beauty_TOFPID_pPb_eta06_iter2_envelope_minsys.root", hBackground, bag->tpcdEdxcutlow[0]);
            } else if (bag->spdcheck == 2){
                printf("hadron cont max beauty");
                if(bag->TOFs==0){
                    printf("hadron cont no TOF\n");
                    status = ReadContaminationFunctions("hadroncontamination_noTOF_pPb_Beauty_ESD_eta06_envelopemax.root", hBackground, bag->tpcdEdxcutlow[0]);
                }
                status = ReadContaminationFunctionsBeauty("hadroncontamination_ESD_Beauty_TOFPID_pPb_eta06_iter2_envelope_maxsys.root", hBackground, bag->tpcdEdxcutlow[0]);
            } else if (bag->spdcheck == 3){
                printf("hadron cont pol3 beauty");
                status = ReadContaminationFunctionsBeauty("hadroncontamination_ESD_Beauty_TOFPID_pPb_eta06_iter2_pol3.root", hBackground, bag->tpcdEdxcutlow[0]);
            }
        }
        else  status = ReadContaminationFunctions("hadroncontamination_TOFTPC_pPb_eta06_newsplines_try3.root", hBackground, bag->tpcdEdxcutlow[0]);
        for(Int_t a=0;a<12;a++) {
            //printf("back %f \n",hBackground[a]);
            if(status) task->SetBackGroundFactorsFunction(hBackground[a],a);
            else printf("not all background functions found\n");
        }
    }

    //***************************************//
    //       Configure NPE plugin            //
    //***************************************//

    AliHFENonPhotonicElectron *backe = new AliHFENonPhotonicElectron(Form("HFEBackGroundSubtractionPID2%s",bag->appendix.Data()),"Background subtraction");  //appendix
    //Setting the Cuts for the Associated electron-pool
    AliHFEcuts *hfeBackgroundCuts = new AliHFEcuts(Form("HFEBackSub%s",bag->appendix.Data()),"Background sub Cuts");
    //  hfeBackgroundCuts->SetEtaRange(assETA);
    hfeBackgroundCuts->SetEtaRange(bag->assETAm,bag->assETAp);
    hfeBackgroundCuts->SetPtRange(bag->assMinPt,20.);

    hfeBackgroundCuts->SetMaxChi2perClusterTPC(4);
    hfeBackgroundCuts->SetMinNClustersITS(bag->assITS);
    hfeBackgroundCuts->SetMinNClustersTPC(bag->assTPCcl);
    hfeBackgroundCuts->SetMinNClustersTPCPID(bag->assTPCPIDcl);
    hfeBackgroundCuts->SetMaxImpactParam(bag->assDCAr,bag->assDCAz);
    if(bag->isAOD) hfeBackgroundCuts->SetAODFilterBit(4);
    hfeBackgroundCuts->SetQAOn();			        // QA

    AliHFEpid *pidbackground = backe->GetPIDBackground();
    if(bag->useMC) pidbackground->SetHasMCData(kTRUE);

    if (bag->etadalwei>0){
        printf("\n\n\n\n\n WEIGHTS FOR ETA DALITZ !!!!!!!!!!!!!!!!!! \n\n\n\n");
        backe->SetEtaDalitzWeightFactor(bag->etadalwei);  
    }   

    if (bag->assTOFs>0.){
        pidbackground->AddDetector("TOF", 0);
        pidbackground->AddDetector("TPC", 1);
    } else {
        pidbackground->AddDetector("TPC", 0);
    }

    Double_t paramsTPCdEdxcutlowAssoc[12] ={-3.0,-3.0,-3.0,-3.0,-3.0,-3.0,-3.0,-3.0,-3.0,-3.0,-3.0,-3.0};
    if(bag->assTPCSminus) memcpy(paramsTPCdEdxcutlowAssoc,bag->assTPCSminus,sizeof(paramsTPCdEdxcutlowAssoc));

    Double_t paramsTPCdEdxcuthighAssoc[12] ={3.0, 3.0, 3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0};
    if(bag->assTPCSplus) memcpy(paramsTPCdEdxcuthighAssoc,bag->assTPCSplus,sizeof(paramsTPCdEdxcuthighAssoc));

    char *cutmodelAssoc;
    cutmodelAssoc="pol0";
    for(Int_t a=0;a<11;a++){
        // Not necessary anymore, since the pPb case is handled similarly to the pp case
        //   cout << a << " " << paramsTPCdEdxcut[a] << endl;
        Double_t tpcparamlow[1]={paramsTPCdEdxcutlowAssoc[a]};
        Float_t tpcparamhigh=paramsTPCdEdxcuthighAssoc[a];
        pidbackground->ConfigureTPCcentralityCut(a,cutmodelAssoc,tpcparamlow,tpcparamhigh);
    }
    pidbackground->ConfigureTPCdefaultCut(cutmodelAssoc,paramsTPCdEdxcutlowAssoc,paramsTPCdEdxcuthighAssoc[0]); // After introducing the pPb flag, pPb is merged with pp and this line defines the cut
    //backe->GetPIDBackgroundQAManager()->SetHighResolutionHistos();

    if (bag->assTOFs>0.){
        pidbackground->ConfigureTOF(bag->TOFs);
    }

    backe->SetHFEBackgroundCuts(hfeBackgroundCuts);

    // Selection of associated tracks for the pool
    if(bag->useCat1Tracks) backe->SelectCategory1Tracks(kTRUE);

    // apply opening angle cut to reduce file size
    backe->SetMaxInvMass(0.3);
    backe->SetPtBinning(sizept, ptbinning);
    backe->SetEtaBinning(sizeeta, etabinning);
    //backe->SetAnaPairGen(kTRUE,2);
    //backe->SetDisplayMCStack();
    // MC weight
    if(bag->useMC) {
        //printf("test put weight %d\n",weightlevelback);
        if((bag->weightlevelback >=0) && (bag->weightlevelback < 3)) backe->SetWithWeights(bag->weightlevelback);
    }
    task->SetHFEBackgroundSubtraction(backe);

    //task->SetWeightHist(); 
    //if(bag->useMC) task->SetDebugStreaming();
    //task->SetCalcContamBeauty(kTRUE);

    // QA
    printf("task %p\n", task);
    //task->SetQAOn(AliAnalysisTaskHFEtemplate::kPIDqa);
    task->SetQAOn(AliAnalysisTaskHFEtemplate::kMCqa);
    task->SwitchOnPlugin(AliAnalysisTaskHFEtemplate::kDEstep);
    if(bag->isBeauty){
        if(nonPhotonicElectronBeauty)task->SwitchOnPlugin(AliAnalysisTaskHFEtemplate::kNonPhotonicElectronBeauty);
    }
    else task->SwitchOnPlugin(AliAnalysisTaskHFEtemplate::kNonPhotonicElectron);

    printf("*************************************\n");
    printf("Configuring standard Task:\n");
    task->PrintStatus();
    pid->PrintStatus();
    printf("*************************************\n");

    if(bag->isAOD)
        task->SetAODAnalysis();
    else
        task->SetESDAnalysis();

    if (bag->useMC)	task->SetHasMCData(kTRUE);
    else		task->SetHasMCData(kFALSE);

    task->SelectCollisionCandidates(AliVEvent::kINT7);

    if(bag->useMC&&(bag->isBeauty || (bag->weightlevelback>=0))) {
        ConfigWeightFactors(task,bag->nonHFEsys,bag->WhichWei,"nonHFEcorrect_pPb.root");
        task->SetNonHFEsystematics(bag->nonHFEsys);
    }

    //create data containers
    AliAnalysisDataContainer *cOutputResults =
        mgr->CreateContainer(Form("HFE_Results_%s",bag->appendix.Data()),
                TList::Class(),
                AliAnalysisManager::kOutputContainer,
                Form("HFE%s.root",bag->appendix.Data()));

    AliAnalysisDataContainer *cOutputQA =
        mgr->CreateContainer(Form("HFE_QA_%s",bag->appendix.Data()),
                TList::Class(),
                AliAnalysisManager::kOutputContainer,
                Form("HFE%s.root",bag->appendix.Data()));

    AliAnalysisDataContainer *cOutputSettings =
        mgr->CreateContainer(Form("HFE_SETTINGS_%s",bag->appendix.Data()),
                TObject::Class(),
                AliAnalysisManager::kParamContainer,
                Form("HFE%s.root",bag->appendix.Data()));


    mgr->ConnectInput(task,  0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task, 1, cOutputResults );
    mgr->ConnectOutput(task, 2, cOutputQA);
    mgr->ConnectOutput(task, 3, cOutputSettings);

    mgr->AddTask(task);

    return NULL;



}
