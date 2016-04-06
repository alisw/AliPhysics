AliAnalysisTask *AddTaskHFEnpepPbMulti(Bool_t MCthere, Bool_t isAOD, Bool_t kNPERef = kTRUE, Bool_t kNPERefTPConly = kFALSE){


    AliHFEparamBag *bag_npepPbMulti = new AliHFEparamBag("mybag");
    // Default settings (TOF-TPC pPb)
    bag_npepPbMulti->useMC=MCthere;
    bag_npepPbMulti->isAOD=isAOD;
    bag_npepPbMulti->isBeauty=0;
    bag_npepPbMulti->TPCcl=100;
    bag_npepPbMulti->TPCclPID=80;
    bag_npepPbMulti->ITScl=4;
    bag_npepPbMulti->DCAxy=1.;
    bag_npepPbMulti->DCAz=2.;
    bag_npepPbMulti->TOFs=3.;
    bag_npepPbMulti->etami=-0.8;
    bag_npepPbMulti->etama=0.8;
    bag_npepPbMulti->phimi=-1;
    bag_npepPbMulti->phima=-1;
    bag_npepPbMulti->itshitpixel=AliHFEextraCuts::kBoth;
    bag_npepPbMulti->spdcheck=0;
    bag_npepPbMulti->assETAm=-0.8;
    bag_npepPbMulti->assETAp=0.8;
    bag_npepPbMulti->assMinPt=0.1;
    bag_npepPbMulti->assITS=2;
    bag_npepPbMulti->assTPCcl=60;
    bag_npepPbMulti->assTPCclPID=60;
    bag_npepPbMulti->assDCAr=1.0;
    bag_npepPbMulti->assDCAz=2.0;
    bag_npepPbMulti->assITSpid=3.0;
    bag_npepPbMulti->assTOFs=0.0;
    bag_npepPbMulti->nonPhotonicElectronBeauty=kFALSE;
    bag_npepPbMulti->nonHFEsys=kFALSE;
    bag_npepPbMulti->ipCharge=kFALSE;
    bag_npepPbMulti->ipOpp=kFALSE;
    bag_npepPbMulti->mcQADebugTree=kFALSE;

    bag_npepPbMulti->RefMulti=19.44;
    bag_npepPbMulti->MultiSystem=AliAnalysisTaskHFEMulti::kNtrk10;

    //bag_npepPbMulti->RefMulti=82.7;
    //bag_npepPbMulti->MultiSystem=AliAnalysisTaskHFEMulti::kVZERO;

    // TPC PID Cuts Inclusive leg:
    // General, if mean=0 and sigma=1: 
    // Above 3 sigma we neglect 0.13%
    // cut in sigma (effective efficiency from cut to 3 sigma)
    // -1 (84%), -0.75 (77.2%), -0.5 (69%), -0.25 (59.7%), -0.129 (55%)
    //  0 (49.9%), 0.122 (45%), 0.25 (40%), 0.5 (30.7%)

    // ESDs: mean 0.06, sigma 1    --> -0.94, -0.69, -0.44,  -0.19, -0.009 ,0.06, 0.182, 0.31, 0.56
    // AODs: mean 0.09, sigma 1    --> -0.91, -0.66, -0.41,  -0.16, -0.039 ,0.09, 0.212, 0.34, 0.59
    // AODs: mean 0.09, sigma 1.03 --> -0.94, -0.68, -0.425, -0.17, -0.043 ,0.09, 0.216, 0.35, 0.60

    // On ESD:
    // mean is actually 0.06 (abs(eta)<0.6)
    Double_t dEdxhm[12] = {3.06,3.06,3.06,3.06,3.06,3.06,3.06,3.06,3.06,3.06,3.06,3.06};  // Above 3 sigma we neglect 0.13%
    Double_t tpcl0[12]  = {-0.94,-0.94,-0.94,-0.94,-0.94,-0.94,-0.94,-0.94,-0.94,-0.94,-0.94,-0.94};  // 84%
    Double_t tpcl1[12]  = {-0.44,-0.44,-0.44,-0.44,-0.44,-0.44,-0.44,-0.44,-0.44,-0.44,-0.44,-0.44};  // 69%
    Double_t tpcl2[12]  = {0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06};  // 50%
    Double_t tpcl3[12]  = {-0.69,-0.69,-0.69,-0.69,-0.69,-0.69,-0.69,-0.69,-0.69,-0.69,-0.69,-0.69};  // 77.2%
    Double_t tpcl4[12]  = {-0.19,-0.19,-0.19,-0.19,-0.19,-0.19,-0.19,-0.19,-0.19,-0.19,-0.19,-0.19};  // 59.7%
    Double_t tpcl5[12]  = {0.186,0.186,0.186,0.186,0.186,0.186,0.186,0.186,0.186,0.186,0.186,0.186};  // 45%
    Double_t tpcl6[12]  = {0.31,0.31,0.31,0.31,0.31,0.31,0.31,0.31,0.31,0.31,0.31,0.31};  // 40%
    Double_t tpcl7[12]  = {0.56,0.56,0.56,0.56,0.56,0.56,0.56,0.56,0.56,0.56,0.56,0.56};  // 30.7%

    // AOD 139 with TPC multi cor; w/o TOF PID
    Double_t dEdxhmAOD[12] = {3.18,3.18,3.18,3.18,3.18,3.18,3.18,3.18,3.18,3.18,3.18,3.18};  // upper cut, 3 sigma
    Double_t dEdxhm1AOD[12] = {2.15,2.15,2.15,2.15,2.15,2.15,2.15,2.15,2.15,2.15,2.15,2.15};  // upper cut, 2 sigma
    Double_t dEdxhm2AOD[12] = {2.67,2.67,2.67,2.67,2.67,2.67,2.67,2.67,2.67,2.67,2.67,2.67};  // upper cut, 2.5 sigma
    Double_t dEdxhm3AOD[12] = {4.21,4.21,4.21,4.21,4.21,4.21,4.21,4.21,4.21,4.21,4.21,4.21};  // upper cut, 4 sigma
    Double_t tpcl0AOD[12]  = {-0.94,-0.94,-0.94,-0.94,-0.94,-0.94,-0.94,-0.94,-0.94,-0.94,-0.94,-0.94};  // 84% <-- hadron contamination
    Double_t tpcl1AOD[12]  = {-0.68,-0.68,-0.68,-0.68,-0.68,-0.68,-0.68,-0.68,-0.68,-0.68,-0.68,-0.68};  // 77.2%
    Double_t tpcl2AOD[12]  = {-0.425,-0.425,-0.425,-0.425,-0.425,-0.425,-0.425,-0.425,-0.425,-0.425,-0.425,-0.425};  // 69% <-- had cont
    Double_t tpcl3AOD[12]  = {-0.17,-0.17,-0.17,-0.17,-0.17,-0.17,-0.17,-0.17,-0.17,-0.17,-0.17,-0.17};  // 59.7%           <-- had cont
    Double_t tpcl4AOD[12]  = {-0.043,-0.043,-0.043,-0.043,-0.043,-0.043,-0.043,-0.043,-0.043,-0.043,-0.043,-0.043};  // 55%
    Double_t tpcl5AOD[12]  = {0.09,0.09,0.09,0.09,0.09,0.09,0.09,0.09,0.09,0.09,0.09,0.09};  // 49.9%                       <-- had cont
    Double_t tpcl6AOD[12]  = {0.216,0.216,0.216,0.216,0.216,0.216,0.216,0.216,0.216,0.216,0.216,0.216};  // 45%
    Double_t tpcl7AOD[12]  = {0.35,0.35,0.35,0.35,0.35,0.35,0.35,0.35,0.35,0.35,0.35,0.35};  // 40%
    Double_t tpcl8AOD[12]  = {0.60,0.60,0.60,0.60,0.60,0.60,0.60,0.60,0.60,0.60,0.60,0.60};  // 30.7%


    // Default setting for the associated electron for the NonPhotonic Analysis
    const double	kassTPCSminus	= -3.0;
    const double	kassTPCSplus	=  3.0;

    enum{
        kHFEV0A = 1,
        kHFEV0M = 2,
        kHFECL1 = 3,
        kHFEZNA = 4
    };

    const Bool_t isBeauty = kFALSE;

    //get the current analysis manager
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        Error("AddTask_hfe_HFE", "No analysis manager found.");
        return 0;
    }

    AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();

    //@@ 0 @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    Double_t dEdxaclm[12], dEdxachm[12],dEdxaclm1[12], dEdxachm1[12],dEdxaclm2[12], dEdxachm2[12];
    for(int icent = 0; icent < 12; icent++){
        dEdxaclm[icent] = kassTPCSminus;
        dEdxachm[icent] = kassTPCSplus;
        dEdxaclm1[icent] = -2;
        dEdxachm1[icent] = 2;
        dEdxaclm2[icent] = -4;
        dEdxachm2[icent] = 4;
    }

    Int_t kWei = -1;
    if (MCthere) kWei = 0;
    enum {
        // updated after changes by Jan on 17.04, commit 21648
        kb2weiData = 3,  // DPMJET weights 
        kd3weiData = 4,  // HIJING weights 
        kd3weib2   = 5,  // HIJING weights for DPMJET
    };

    memcpy(bag_npepPbMulti->tpcdEdxcuthigh, dEdxhmAOD, sizeof(Double_t)*12);
    memcpy(bag_npepPbMulti->assTPCSminus,   dEdxaclm,  sizeof(Double_t)*12);
    memcpy(bag_npepPbMulti->assTPCSplus,    dEdxachm, sizeof(Double_t)*12);
    bag_npepPbMulti->icent=kHFEV0A;
    bag_npepPbMulti->useCat1Tracks=kTRUE;
    bag_npepPbMulti->useCat2Tracks=kFALSE;
    bag_npepPbMulti->weightlevelback=kWei;
    bag_npepPbMulti->WhichWei=kd3weiData;
    bag_npepPbMulti->etadalwei=0;


    if(kNPERef){
        // **************************************************************
        // 
        // Reference task
        //
        // **************************************************************


        if(isAOD==1){
    //bag_npepPbMulti->RefMulti=82.7;
    //bag_npepPbMulti->MultiSystem=AliAnalysisTaskHFEMulti::kVZERO;
            //TPC1
            AliHFEparamBag *bag1 = new AliHFEparamBag(*bag_npepPbMulti);
            bag1->appendix = "TPCTOF1";
            memcpy(bag1->tpcdEdxcutlow,  tpcl5AOD,  sizeof(Double_t)*12);
            RegisterTaskNPEpPbMulti(bag1);

            AliHFEparamBag *bagX = new AliHFEparamBag(*bag_npepPbMulti);
            bagX->RefMulti=82.7;
            bagX->MultiSystem=AliAnalysisTaskHFEMulti::kVZERO;
            bagX->appendix = "TPCTOF_NV0A";
            RegisterTaskNPEpPbMulti(bagX);
//DCA
//ITS
//SPD
//TPC cluster
/*
            //TPC2
            AliHFEparamBag *bag2 = new AliHFEparamBag(*bag_npepPbMulti);
            bag2->TPCcl=90;
            bag2->appendix = new TString("TPCTOF_TPCcl1");
            RegisterTaskNPEpPbMulti(bag2);

            //TPC3
            AliHFEparamBag *bag3 = new AliHFEparamBag(*bag_npepPbMulti);
            bag3->TPCcl=95;
            bag3->appendix = new TString("TPCTOF_TPCcl2");
            RegisterTaskNPEpPbMulti(bag3);

            //TPC4
            AliHFEparamBag *bag4 = new AliHFEparamBag(*bag_npepPbMulti);
            bag4->TPCcl=105;
            bag4->appendix = new TString("TPCTOF_TPCcl3");
            RegisterTaskNPEpPbMulti(bag4);

            //TPC5
            AliHFEparamBag *bag5 = new AliHFEparamBag(*bag_npepPbMulti);
            bag5->TPCcl=110;
            bag5->appendix = new TString("TPCTOF_TPCcl4");
            RegisterTaskNPEpPbMulti(bag5);

            //TPC6
            AliHFEparamBag *bag6 = new AliHFEparamBag(*bag_npepPbMulti);
            bag6->TPCcl=115;
            bag6->appendix = new TString("TPCTOF_TPCcl5");
            RegisterTaskNPEpPbMulti(bag6);

            //TPC7
            AliHFEparamBag *bag7 = new AliHFEparamBag(*bag_npepPbMulti);
            bag7->TPCcl=120;
            bag7->appendix = new TString("TPCTOF_TPCcl6");
            RegisterTaskNPEpPbMulti(bag7);

//TPC PID cluster
            //TPC2
            AliHFEparamBag *bag8 = new AliHFEparamBag(*bag_npepPbMulti);
            bag8->TPCclPID=60;
            bag8->appendix = new TString("TPCTOF_TPCclPID1");
            RegisterTaskNPEpPbMulti(bag8);

            //TPC3
            AliHFEparamBag *bag9 = new AliHFEparamBag(*bag_npepPbMulti);
            bag9->TPCclPID=70;
            bag9->appendix = new TString("TPCTOF_TPCclPID2");
            RegisterTaskNPEpPbMulti(bag9);

            //TPC4
            AliHFEparamBag *bag10 = new AliHFEparamBag(*bag_npepPbMulti);
            bag10->TPCclPID=90;
            bag10->appendix = new TString("TPCTOF_TPCclPID3");
            RegisterTaskNPEpPbMulti(bag10);

            //TPC5
            AliHFEparamBag *bag11 = new AliHFEparamBag(*bag_npepPbMulti);
            bag11->TPCclPID=100;
            bag11->appendix = new TString("TPCTOF_TPCclPID4");
            RegisterTaskNPEpPbMulti(bag11);
//TOF PID
//TPC PID

 */


        }
        else {
            // Reference
            memcpy(bag_npepPbMulti->tpcdEdxcutlow,  tpcl1,  sizeof(Double_t)*12);
            RegisterTaskNPEpPbMulti(bag_npepPbMulti);
        }
    }

    if(kNPERefTPConly){
        // **************************************************************
        // 
        // Reference task for TPC-only on the inclusive leg
        // From July 21, without ITS SA
        // **************************************************************
        if(isAOD==1){
            // Reference: 50% PID cut. For the moment, old contamination file

            AliHFEparamBag *bagY = new AliHFEparamBag(*bag_npepPbMulti);
            bagY->TOFs=0;
            bagY->appendix = "TPConly1";
            memcpy(bagY->tpcdEdxcutlow,  tpcl5AOD,  sizeof(Double_t)*12);
            RegisterTaskNPEpPbMulti(bagY);

            AliHFEparamBag *bagY2 = new AliHFEparamBag(*bag_npepPbMulti);
            bagY2->TOFs=0;
            bagY2->appendix = "TPConly_NV0A";
            bagY2->RefMulti=82.7;
            bagY2->MultiSystem=AliAnalysisTaskHFEMulti::kVZERO;
            RegisterTaskNPEpPbMulti(bagY2);
/*
 *
            //TPC2
            AliHFEparamBag *bag9 = new AliHFEparamBag(*bag_npepPbMulti);
            bag9->TPCcl=120;
            bag9->appendix = new TString("TPConly2");
            RegisterTaskNPEpPbMulti(bag9);

            //TPC3
            AliHFEparamBag *bag10 = new AliHFEparamBag(*bag_npepPbMulti);
            bag10->TPCcl=80;
            bag10->appendix = new TString("TPConly3");
            RegisterTaskNPEpPbMulti(bag10);

            //TPC4
            AliHFEparamBag *bag11 = new AliHFEparamBag(*bag_npepPbMulti);
            bag11->appendix = new TString("TPConly4");
            memcpy(bag11->tpcdEdxcutlow,  tpcl2AOD,  sizeof(Double_t)*12);
            RegisterTaskNPEpPbMulti(bag11);

            //TPC5
            AliHFEparamBag *bag12 = new AliHFEparamBag(*bag_npepPbMulti);
            bag12->appendix = new TString("TPConly5");
            memcpy(bag12->tpcdEdxcutlow,  tpcl7AOD,  sizeof(Double_t)*12);
            RegisterTaskNPEpPbMulti(bag12);

            //TPC6
            AliHFEparamBag *bag13 = new AliHFEparamBag(*bag_npepPbMulti);
            bag13->TPCclPID=100;
            bag13->appendix = new TString("TPConly6");
            RegisterTaskNPEpPbMulti(bag13);

            //TPC7
            AliHFEparamBag *bag14 = new AliHFEparamBag(*bag_npepPbMulti);
            bag14->TPCclPID=60;
            bag14->appendix = new TString("TPConly7");
            RegisterTaskNPEpPbMulti(bag14);
 */
        }
    }
    return NULL;
}

AliAnalysisTask RegisterTaskNPEpPbMulti(AliHFEparamBag *abag){

    AliHFEparamBag *bag = new AliHFEparamBag(*abag);


    printf("Add macro appendix %s\n", bag->appendix.Data());

 if(!gROOT->GetListOfGlobalFunctions()->FindObject("ConfigWeightFactors")) 
     gROOT->LoadMacro("$ALICE_PHYSICS/PWGHF/hfe/macros/configs/ConfigWeightFactors.C");
 if(!gROOT->GetListOfGlobalFunctions()->FindObject("ConfigHFEnpepPbTOFworkaround")) 
     gROOT->LoadMacro("$ALICE_PHYSICS/PWGHF/hfe/macros/configs/pPb/ConfigHFEnpepPbMulti.C");

    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

    AliAnalysisTaskHFEMulti *task = ConfigHFEnpepPbMulti(bag);

    if(bag->isAOD)
        task->SetAODAnalysis();
    else
        task->SetESDAnalysis();

    if (bag->useMC)	task->SetHasMCData(kTRUE);
    else		task->SetHasMCData(kFALSE);

    task->SelectCollisionCandidates(AliVEvent::kINT7);

    if(bag->useMC&&(bag->isBeauty || (bag->weightlevelback>=0))) {
        ConfigWeightFactors(task,bag->nonHFEsys,bag->WhichWei);
        task->SetNonHFEsystematics(bag->nonHFEsys);
    }

    //Set multiplicity estimator
    task->SetMultiEstimator(bag->MultiSystem);

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
