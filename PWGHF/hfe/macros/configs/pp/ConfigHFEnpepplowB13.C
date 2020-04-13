AliAnalysisTaskHFE* ConfigHFEnpepplowB13(Bool_t useMC, Bool_t isAOD, TString appendix,
                                   UChar_t TPCcl=70, UChar_t TPCclPID = 80,
                                   UChar_t ITScl=3, Double_t DCAxy=1000., Double_t DCAz=1000.,
                                   Double_t* tpcdEdxcutlow=NULL, Double_t* tpcdEdxcuthigh=NULL,
                                   Double_t TOFs=3., Int_t TOFmis=0,
                                   Int_t itshitpixel = 0, Int_t ikink = 0,
                                   Double_t etami=-0.8, Double_t etama=0.8,
                                   Double_t assETAm=-0.8, Double_t assETAp=0.8,
                                   Double_t assMinPt=0.2, Int_t assITS=2,
                                   Int_t assTPCcl=100, Int_t assTPCPIDcl=80,
                                   Double_t assDCAr=1.0, Double_t assDCAz=2.0,
                                   Double_t *assTPCSminus=NULL, Double_t *assTPCSplus=NULL,
                                   Double_t assITSpid=-3., Double_t assTOFs=3.,
                                   Bool_t useCat1Tracks = kTRUE, Bool_t useCat2Tracks = kTRUE,
                                   Int_t weightlevelback = -1, Int_t HadronContFunc=0 , Int_t Chi2perTPCcluster=4)
{
    
    //***************************************//
    //        Setting up the HFE cuts        //
    //***************************************//
    
    AliHFEcuts *hfecuts = new AliHFEcuts(appendix,"HFE cuts for pp");
    //hfecuts->SetQAOn();
    hfecuts->CreateStandardCuts();
    if(isAOD) hfecuts->SetAODFilterBit(4);
    
    hfecuts->SetMinNClustersTPC(TPCcl);
    hfecuts->SetMinNClustersTPCPID(TPCclPID);
    hfecuts->SetMinNClustersITS(ITScl);
    hfecuts->SetMinRatioTPCclusters(0.6);
    hfecuts->SetTPCmodes(AliHFEextraCuts::kFound, AliHFEextraCuts::kFoundOverFindable);
    hfecuts->SetCutITSpixel(itshitpixel);
    hfecuts->SetCheckITSLayerStatus(kFALSE);
    hfecuts->SetEtaRange(etami,etama);
    hfecuts->SetRejectKinkDaughters();
    if (ikink == 0){
        hfecuts->SetRejectKinkMothers();
    } else {
        hfecuts->SetAcceptKinkMothers();
    }
    hfecuts->SetMaxImpactParam(DCAxy,DCAz);
    
    hfecuts->SetUseMixedVertex(kTRUE);
    hfecuts->SetUseSPDVertex(kFALSE);
    hfecuts->SetUseTrackVertex(kFALSE);
    
    hfecuts->SetVertexRange(10.);
    if (hfecuts->GetUseSPDVertex() || hfecuts->GetUseMixedVertex() ) {
        //if SPD vertices are used acivate resolution cut
        printf("########## CONFIGURATION FILE: Switch is working ######### \n");
        hfecuts->SetSPDVtxResolutionCut();
    }
    
    
    // TOF settings:
    Int_t usetof=0;
    Bool_t kTOFmis=kFALSE;
    if (TOFs>0.){
        usetof = 1;
        printf("CONFIGURATION FILE: TOF is used \n");
        hfecuts->SetTOFPIDStep(kTRUE);
        printf("CONFIGURATION FILE: TOF PID step is requested !!!! \n");
        if (TOFmis>0){
            kTOFmis = kTRUE;
            printf("CONFIGURATION FILE: TOF mismatch rejection is set ON \n");
        }
    }
    
    //***************************************//
    //        Setting up the task            //
    //***************************************//
    
    AliAnalysisTaskHFE *task = new AliAnalysisTaskHFE(Form("HFEtask%s",appendix.Data()));
    printf("task %p\n", task);
    task->SetppAnalysis();
    //if(!isAOD) task->SetRemoveFirstEventInChunk();
    task->SetRemovePileUp(kTRUE);
    task->SetHFECuts(hfecuts);
    task->GetPIDQAManager()->SetHighResolutionHistos();
    
    //***************************************//
    //          Variable manager             //
    //***************************************//
    // Define Variables
    //standard binning
    Double_t ptbinning[36] = {0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.75, 2., 2.25, 2.5, 2.75, 3., 3.5, 4., 4.5, 5., 5.5, 6., 7., 8., 10., 12., 14., 16., 18., 20.};
    Double_t etabinning[17] = {-0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};
    
    Double_t phibinning[5] = {0, 0.5*TMath::Pi(), TMath::Pi(), 1.5*TMath::Pi(),2*TMath::Pi()};
    
    Int_t sizept=(sizeof(ptbinning)/sizeof(double))-1;
    Int_t sizeeta=(sizeof(etabinning)/sizeof(double))-1;
    Int_t sizephi=(sizeof(phibinning)/sizeof(double))-1;
    
    AliHFEvarManager *vm = task->GetVarManager();
    vm->AddVariable("pt", sizept, ptbinning);
    vm->AddVariable("eta", sizeeta, -0.8,0.8);
    vm->AddVariable("phi", sizephi, phibinning);
    vm->AddVariable("charge");
    vm->AddVariable("source");
    
    // The part dedicated to the background subtraction
    // should be implemented in a different way, reading it from a root file.
    
    if(!useMC){
      
      	// For eta < 0.8
       TF1 *hBackground;   // 0 = Landau + 2 Gaus, 1 =  Error + 2 Gaus
        // First hadron contamination fit for 13TeV  by Andrea
        // relative to the case of a TPC+TOF PID cut at -1 sigma // flipped on Dec 10 2019 by Sudhir
        if(HadronContFunc == 2){
        hBackground = new TF1("hadronicBackgroundFunction", "[0]*TMath::Landau(x,[1],[2]) + [3]*TMath::Gaus(x, [4], [5]) + [6]*TMath::Gaus(x, [7], [8])", 0. ,60.);
    	hBackground->SetParameter(0, 3.98051e+00);
   	hBackground->SetParameter(1, 1.07724e+01);
  	hBackground->SetParameter(2, 2.55463e+00);
 	hBackground->SetParameter(3, 4.16569e-03);
	hBackground->SetParameter(4, 8.60415e-01);
	hBackground->SetParameter(5, 4.66792e-02);
 	hBackground->SetParameter(6, 6.21460e-04);
        hBackground->SetParameter(7, 4.46974e-01);
        hBackground->SetParameter(8, 2.54008e-02);

        }
        if(HadronContFunc == 3){
        hBackground = new TF1("hadronicBackgroundFunction", "[0]+[1]*TMath::Erf([2]*x+[3]) + [4]*TMath::Gaus(x, [5], [6]) + [7]*TMath::Gaus(x, [8], [9])", 0. ,60.);
        hBackground->SetParameter(0, 2.20105e-01);
        hBackground->SetParameter(1, 2.20105e-01);
        hBackground->SetParameter(2, 6.48606e-01);
        hBackground->SetParameter(3,-3.81881e+00);
        hBackground->SetParameter(4, 4.16489e-03);
        hBackground->SetParameter(5, 8.60414e-01);
        hBackground->SetParameter(6, 4.66569e-02);
        hBackground->SetParameter(7, 6.18557e-04);
        hBackground->SetParameter(8, 4.47756e-01);
        hBackground->SetParameter(9, 2.54008e-02);
        }
        // For eta < 0.5
        if(HadronContFunc == 0){
        hBackground = new TF1("hadronicBackgroundFunction", "[0]*TMath::Landau(x,[1],[2]) + [3]*TMath::Gaus(x, [4], [5]) + [6]*TMath::Gaus(x, [7], [8]) + [9]*TMath::Gaus(x, [10], [11]) + [12]*TMath::Gaus(x, [13], [14])", 0. ,60.);
    	hBackground->SetParameter(0, 5.19953e+02);
   	hBackground->SetParameter(1, 1.67295e+01);
  	hBackground->SetParameter(2, 3.93376e+00);
 	hBackground->SetParameter(3, 9.61806e-05);
        hBackground->SetParameter(4, 0.25);
        hBackground->SetParameter(5, 3.14008e-02);
        hBackground->SetParameter(6, 3.2e-03);
        hBackground->SetParameter(7, 5.00000e-01);
        hBackground->SetParameter(8, 5.52508e-02);
        hBackground->SetParameter(9, 1.2e-02);
        hBackground->SetParameter(10, 0.85);
        hBackground->SetParameter(11, 7.84008e-02);
        hBackground->SetParameter(12, 2e-02);
        hBackground->SetParameter(13, 1.05);
        hBackground->SetParameter(14, 2.14008e-02);

        }
        if(HadronContFunc == 1){
        hBackground = new TF1("hadronicBackgroundFunction", "[0]+[1]*TMath::Erf([2]*x+[3]) + [4]*TMath::Gaus(x, [5], [6]) + [7]*TMath::Gaus(x, [8], [9]) + [10]*TMath::Gaus(x, [11], [12]) + [13]*TMath::Gaus(x, [14], [15])", 0. ,60.);
        hBackground->SetParameter(0, 2.20155e-01);
        hBackground->SetParameter(1, 2.20165e-01);
        hBackground->SetParameter(2, 6.80038e-01);
        hBackground->SetParameter(3,-3.62570e+00);
        hBackground->SetParameter(4, 1.00000e-04);
        hBackground->SetParameter(5, 0.25);
        hBackground->SetParameter(6, 3.14008e-02);
        hBackground->SetParameter(7, 3.2e-03);
        hBackground->SetParameter(8, 5.00000e-01);
        hBackground->SetParameter(9, 5.52508e-02);
        hBackground->SetParameter(10, 1.2e-02);
        hBackground->SetParameter(11, 0.85);
        hBackground->SetParameter(12, 7.64008e-02);
        hBackground->SetParameter(13, 2e-02);
        hBackground->SetParameter(14, 1.05);
        hBackground->SetParameter(15, 2.14008e-02);
        }
        
       /* if(HadronContFunc == 2){
            hBackground= new TF1("hadronicBackgroundFunction", "[0]*TMath::Landau(x,[1],[2])", 0. ,60.);
            hBackground->SetParameter(0, 7.61936e+03);
            hBackground->SetParameter(1, 2.46441e+01);
            hBackground->SetParameter(2, 5.82430e+00);
            //error function
        }
        if(HadronContFunc == 3){
            hBackground = new TF1("hadronicBackgroundFunction", "[0]+[1]*TMath::Erf([2]*x+[3])", 0. ,60.);
            hBackground->SetParameter(0, 9.22700e-01);
            hBackground->SetParameter(1, 9.22697e-01);
            hBackground->SetParameter(2, 5.96595e-01);
            hBackground->SetParameter(3,-3.85842e+00);
        }*/
        
        task->SetBackGroundFactorsFunction(hBackground);
        
    }
    
    //***************************************//
    //          Configure the PID            //
    //***************************************//
    
    // Define PID
    AliHFEpid *pid = task->GetPID();
    if(useMC) pid->SetHasMCData(kTRUE);
    
    if (usetof){
        pid->AddDetector("TOF", 0);
        pid->AddDetector("TPC", 1);
    } else {
        pid->AddDetector("TPC", 0);
    }
    
    // Configure TPC PID
    // do the identical thing in data and MC
    Double_t paramsTPCdEdxcutlow[12] ={0.0, 0.0, 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
    if(tpcdEdxcutlow) memcpy(paramsTPCdEdxcutlow,tpcdEdxcutlow,sizeof(paramsTPCdEdxcutlow));
    
    Double_t paramsTPCdEdxcuthigh[12] ={3.0, 3.0, 3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0};
    if(tpcdEdxcuthigh) memcpy(paramsTPCdEdxcuthigh,tpcdEdxcuthigh,sizeof(paramsTPCdEdxcuthigh));
    
    char *cutmodel;
    
    //  if(useMC){ // constant (default) cut for MC
    cutmodel="pol0(0)";
    Double_t params[1];
    params[0]=paramsTPCdEdxcutlow[0];
    pid->ConfigureTPCdefaultCut(cutmodel, params,tpcdEdxcuthigh[0]);
    /*
     } else { // correct for mean shift in data
     cutmodel="min(pol1(0),pol0(2))";
     Double_t params[3];
     //params[0]=-0.12; params[1]=0.14; params[2]=0.09;
     params[0]=-0.21 + paramsTPCdEdxcutlow[0];
     params[1]=0.14;
     params[2]=paramsTPCdEdxcutlow[0];
     pid->ConfigureTPCdefaultCut(cutmodel, params,tpcdEdxcuthigh[0]);
     }
     */
    
    // Configure TOF PID
    if (usetof){
        pid->ConfigureTOF(TOFs);
        AliHFEpidTOF *tofpid = pid->GetDetPID(AliHFEpid::kTOFpid);
        if (kTOFmis){
            tofpid->SetRejectTOFmismatch();
        }
    }
    
    // To make different upper TOF cut to see contamination effect
    // The below two lines should be removed after this check
    //AliHFEpidTOF *tofpid = pid->GetDetPID(AliHFEpid::kTOFpid);
    //if(TOFs<3.) tofpid->SetTOFnSigmaBand(-3,TOFs); //only to check the assymmetric tof cut
    
    
    //***************************************//
    //       Configure NPE plugin            //
    //***************************************//
    
    AliHFENonPhotonicElectron *backe = new AliHFENonPhotonicElectron(Form("HFEBackGroundSubtractionPID2%s",appendix.Data()),"Background subtraction");  //appendix
    //Setting the Cuts for the Associated electron-pool
    AliHFEcuts *hfeBackgroundCuts = new AliHFEcuts(Form("HFEBackSub%s",appendix.Data()),"Background sub Cuts");
    hfeBackgroundCuts->SetEtaRange(assETAm,assETAp);
    hfeBackgroundCuts->SetPtRange(assMinPt,20.);
    
    hfeBackgroundCuts->SetMaxChi2perClusterTPC(Chi2perTPCcluster);
    
    hfeBackgroundCuts->SetMinNClustersITS(assITS);
    hfeBackgroundCuts->SetMinNClustersTPC(assTPCcl);
    hfeBackgroundCuts->SetMinNClustersTPCPID(assTPCPIDcl);
    hfeBackgroundCuts->SetMaxImpactParam(assDCAr,assDCAz);
    if(isAOD) hfeBackgroundCuts->SetAODFilterBit(0); // Standard TPC only tracks
    hfeBackgroundCuts->SetQAOn();			        // QA
    
    AliHFEpid *pidbackground = backe->GetPIDBackground();
    if(useMC) pidbackground->SetHasMCData(kTRUE);
    
    if (assTOFs>0.){
        pidbackground->AddDetector("TOF", 0);
        pidbackground->AddDetector("TPC", 1);
    } else {
        pidbackground->AddDetector("TPC", 0);
    }
    
    Double_t paramsTPCdEdxcutlowAssoc[12] ={-3.0,-3.0,-3.0,-3.0,-3.0,-3.0,-3.0,-3.0,-3.0,-3.0,-3.0,-3.0};
    if(assTPCSminus) memcpy(paramsTPCdEdxcutlowAssoc,assTPCSminus,sizeof(paramsTPCdEdxcutlowAssoc));
    
    Double_t paramsTPCdEdxcuthighAssoc[12] ={3.0, 3.0, 3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0};
    if(assTPCSplus) memcpy(paramsTPCdEdxcuthighAssoc,assTPCSplus,sizeof(paramsTPCdEdxcuthighAssoc));
    
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
    
    if (assTOFs>0.){
        pidbackground->ConfigureTOF(TOFs);
    }
    
    backe->SetHFEBackgroundCuts(hfeBackgroundCuts);
    
    // Selection of associated tracks for the pool
    if(useCat1Tracks) backe->SelectCategory1Tracks(kTRUE);
    if(useCat2Tracks){
        backe->SelectCategory2Tracks(kTRUE);
        backe->SetITSMeanShift(-0.5);
        backe->SetITSnSigmaHigh(assITSpid);
        Double_t assITSminus = -1.0 * assITSpid;
        backe->SetITSnSigmaLow(assITSminus);
        //backe->SetminPt(assMinPt);
    }
    
    // apply opening angle cut to reduce file size
    backe->SetMaxInvMass(0.6);
    
    Double_t massbinning[32] = {0.0,0.005,0.01,0.015,0.02,0.025,0.03,0.035,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.16,0.18,0.2,0.24,0.28,0.32,0.36,0.4,0.44,0.48,0.52,0.56,0.6};
    Int_t sizemass=(sizeof(massbinning)/sizeof(double))-1;
    
    backe->SetPtBinning(sizept, ptbinning);
    backe->SetEtaBinning(sizeeta, etabinning);
    backe->SetInvMassBinning(sizemass, massbinning);
    backe->SetPhiBinning(sizephi,phibinning);
    
    if(useMC) {
        if((weightlevelback >=0) && (weightlevelback <3)) backe->SetWithWeights(weightlevelback);
    }
    
    task->SetHFEBackgroundSubtraction(backe);
    task->SetWeightHist();
    
    // Also May 30: try to add tagged tracks (photon conversions)
    // not used since long ...
    AliHFEcuts *v0trackCuts = new AliHFEcuts("V0trackCuts", "Track Cuts for tagged track Analysis");
    v0trackCuts->CreateStandardCuts();
    v0trackCuts->SetMinNClustersTPC(TPCcl);
    v0trackCuts->SetMinNClustersTPCPID(TPCclPID);
    v0trackCuts->SetMinRatioTPCclusters(0.6);
    v0trackCuts->SetTPCmodes(AliHFEextraCuts::kFound, AliHFEextraCuts::kFoundOverFindable);
    v0trackCuts->SetMinNClustersITS(1);
    v0trackCuts->SetCutITSpixel(AliHFEextraCuts::kFirst);
    v0trackCuts->SetCheckITSLayerStatus(kFALSE);
    v0trackCuts->UnsetVertexRequirement();
    //hfecuts->SetSigmaToVertex(10);
    //v0trackCuts->SetTOFPIDStep(kTRUE);
    v0trackCuts->SetQAOn();
    
    task->SwitchOnPlugin(AliAnalysisTaskHFE::kTaggedTrackAnalysis);
    task->SetTaggedTrackCuts(v0trackCuts);
    task->SetCleanTaggedTrack(kTRUE);
    // end tagged tracks
    
    // QA
    printf("task %p\n", task);
    task->SetQAOn(AliAnalysisTaskHFE::kPIDqa);
    task->SetQAOn(AliAnalysisTaskHFE::kMCqa);
    task->SwitchOnPlugin(AliAnalysisTaskHFE::kDEstep);
    task->SwitchOnPlugin(AliAnalysisTaskHFE::kNonPhotonicElectron);
    
    printf("*************************************\n");
    printf("Configuring standard Task:\n");
    task->PrintStatus();
    pid->PrintStatus();
    printf("*************************************\n");
    return task;
}
