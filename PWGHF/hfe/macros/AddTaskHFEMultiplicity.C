AliAnalysisTask* AddTaskHFEMultiplicity(TString suffixName = "",
					Bool_t readMC	=	kFALSE,
                    Bool_t SwitchPi0EtaWeightCalc = kTRUE,
					Bool_t PhysSelINT7 =  kTRUE,
					Bool_t useTender   =  kTRUE,
					Bool_t ClsTypeEMC  =  kTRUE,
					Bool_t ClsTypeDCAL =  kTRUE,
					TString estimatorFilename="",
                    Bool_t zvtxcut =  kTRUE,
                    Bool_t zvtxQA =  kTRUE,
					Double_t refMult=61.2,
					Int_t NcontV =2,
					Int_t MaxTPCclus = 100.,
					Int_t TPCNclus =80.,
  					Int_t ITSNclus =3.,
  					Double_t DCAxyCut =2.4,
  					Double_t DCAzCut =3.2,
  					Double_t Etarange = 0.7,
  					Double_t TrackPtMin = 1.,
  					Double_t EopEMin = 0.8,
 					Double_t EopEMax = 1.2,
  					Double_t TPCNSigMin = -1.,
  					Double_t TPCNSigMax = 3.,
  					Double_t M20Min =0.02,
  					Double_t M20Max =0.35,
  					Int_t AssoTPCCluster = 80.,
  					Int_t AssoITSCluster =3.,
  					Double_t AssoEPt =0.1,
  					Double_t AssoEEtarange =0.9,
  					Double_t AssoENsigma =3.,
  					Bool_t AssoITSRefit =kTRUE,
  					Double_t InvmassCut = 0.14					)

{ // Get the pointer to the existing analysis manager via the static access method.
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    return 0x0;
  }
  
  if (!mgr->GetInputEventHandler()) {
    return 0x0;
  }
  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"

  Bool_t  fReadMC = kTRUE;
  AliMCEventHandler *mcH = dynamic_cast<AliMCEventHandler*>(mgr->GetMCtruthEventHandler());
  if(!mcH){
    fReadMC=kFALSE;
  }


  TString fileName = AliAnalysisManager::GetCommonFileName();
  fileName += ":MyTask";
  
  if(PhysSelINT7){
    AliAnalysisTaskHFEMultiplicity* HFEtaskINT7 = new AliAnalysisTaskHFEMultiplicity("");
    mgr->AddTask(HFEtaskINT7); //HFEtask is my task
    HFEtaskINT7->SetReadMC(readMC);
    HFEtaskINT7->SelectCollisionCandidates(AliVEvent::kINT7);
    HFEtaskINT7->SetTenderSwitch(useTender);
    HFEtaskINT7->SwitchPi0EtaWeightCalc(SwitchPi0EtaWeightCalc);
    HFEtaskINT7->SetClusterTypeEMC(ClsTypeEMC);
    HFEtaskINT7->SetClusterTypeDCAL(ClsTypeDCAL);
    HFEtaskINT7->Setzvtxcut(zvtxcut);
    HFEtaskINT7->SetzvtxQA(zvtxQA);
    HFEtaskINT7->SetNcontVCut(NcontV);
    HFEtaskINT7->SetEtaRange(Etarange);
    HFEtaskINT7->SetMaxTPCCluster(MaxTPCclus);
    HFEtaskINT7->SetNTPCCluster(TPCNclus);
    HFEtaskINT7->SetNITSCluster(ITSNclus);
    HFEtaskINT7->SetTrackpTMin(TrackPtMin);
    HFEtaskINT7->SetDCACut(DCAxyCut,DCAzCut);
    HFEtaskINT7->SetTPCnsigma(TPCNSigMin,TPCNSigMax);
    HFEtaskINT7->SetEopE(EopEMin,EopEMax);
    HFEtaskINT7->SetShowerShapeEM20(M20Min,M20Max);
    HFEtaskINT7->SetInvMassCut(InvmassCut);
    HFEtaskINT7->SetAssoTPCclus(AssoTPCCluster);
    HFEtaskINT7->SetAssoITSclus(AssoITSCluster);
    HFEtaskINT7->SetAssoITSrefit(AssoITSRefit);
    HFEtaskINT7->SetAssopTMin(AssoEPt);
    HFEtaskINT7->SetAssoEtarange(AssoEEtarange);
    HFEtaskINT7->SetAssoTPCnsig(AssoENsigma);
 
    TFile* fileEstimator=TFile::Open(estimatorFilename.Data());
    if(!fileEstimator)  {
        Printf("File with multiplicity estimator not found\n");
      return NULL;
    }
    HFEtaskINT7->SetReferenceMultiplicity(refMult);
    const Char_t* profilebasename="SPDmult";
    
    const Char_t* periodNames[2] = {"LHC16s", "LHC16r"};
    TProfile* multEstimatorAvg[2];
    for(Int_t ip=0; ip<2; ip++) {
      cout<< " Trying to get "<<Form("%s_%s",profilebasename,periodNames[ip])<<endl;
      multEstimatorAvg[ip] = (TProfile*)(fileEstimator->Get(Form("%s_%s",profilebasename,periodNames[ip]))->Clone(Form("%s_%s_clone",profilebasename,periodNames[ip])));
      if (!multEstimatorAvg[ip]) {
	Printf("Multiplicity estimator for %s not found! Please check your estimator file",periodNames[ip]);
	return NULL;
      }
    }

    HFEtaskINT7->SetMultiProfileLHC16s(multEstimatorAvg[0]);
    HFEtaskINT7->SetMultiProfileLHC16r(multEstimatorAvg[1]);

    // Create containers for input/output
    TString finDirname         = "_INT7";
    TString outBasicname       = "EID";
    TString profname       = "coutputProf";
      
    finDirname 	      +=   suffixName.Data();
    outBasicname      +=   finDirname.Data();
    profname          +=   finDirname.Data();
      
      
      
    mgr->ConnectInput(HFEtaskINT7,0,mgr->GetCommonInputContainer());
    mgr->ConnectOutput(HFEtaskINT7,1,mgr->CreateContainer(outBasicname, TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
     mgr->ConnectOutput(HFEtaskINT7,2,mgr->CreateContainer(profname, TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));

   

 }
    
   if(!PhysSelINT7){
     

// EMCal EGA GA1
      
     AliAnalysisTaskHFEMultiplicity* HFEtaskGA1 = new AliAnalysisTaskHFEMultiplicity("");
     mgr->AddTask(HFEtaskGA1);
     HFEtaskGA1->SetReadMC(readMC);
     HFEtaskGA1->SelectCollisionCandidates(AliVEvent::kEMCEGA);
     HFEtaskGA1->SetEMCalTriggerEG1(kTRUE);
     HFEtaskGA1->SetEMCalTriggerDG1(kTRUE);
     HFEtaskGA1->SetTenderSwitch(useTender);
     HFEtaskGA1->SwitchPi0EtaWeightCalc(SwitchPi0EtaWeightCalc);
     HFEtaskGA1->SetClusterTypeEMC(ClsTypeEMC);
     HFEtaskGA1->SetClusterTypeDCAL(ClsTypeDCAL);
     HFEtaskGA1->Setzvtxcut(zvtxcut);
     HFEtaskGA1->SetzvtxQA(zvtxQA);
     HFEtaskGA1->SetNcontVCut(NcontV);
     HFEtaskGA1->SetEtaRange(Etarange);
     HFEtaskGA1->SetMaxTPCCluster(MaxTPCclus);
     HFEtaskGA1->SetNTPCCluster(TPCNclus);
     HFEtaskGA1->SetNITSCluster(ITSNclus);
     HFEtaskGA1->SetTrackpTMin(TrackPtMin);
     HFEtaskGA1->SetDCACut(DCAxyCut,DCAzCut);
     HFEtaskGA1->SetTPCnsigma(TPCNSigMin,TPCNSigMax);
     HFEtaskGA1->SetEopE(EopEMin,EopEMax);
     HFEtaskGA1->SetShowerShapeEM20(M20Min,M20Max);
     HFEtaskGA1->SetInvMassCut(InvmassCut);
     HFEtaskGA1->SetAssoTPCclus(AssoTPCCluster);
     HFEtaskGA1->SetAssoITSclus(AssoITSCluster);
     HFEtaskGA1->SetAssoITSrefit(AssoITSRefit);
     HFEtaskGA1->SetAssopTMin(AssoEPt);
     HFEtaskGA1->SetAssoEtarange(AssoEEtarange);
     HFEtaskGA1->SetAssoTPCnsig(AssoENsigma);

    TFile* fileEstimatorGA1=TFile::Open(estimatorFilename.Data());
    if(!fileEstimatorGA1)  {
      Printf("File with multiplicity estimator not found\n");
      return NULL;
    }
    HFEtaskGA1->SetReferenceMultiplicity(refMult);
    const Char_t* profilebasenameEG1="SPDmultEG1";
       const Char_t* periodNamesGA1[2] = {"LHC16s", "LHC16r"};
    TProfile* multEstimatorAvgGA1[2];
    for(Int_t ip=0; ip<2; ip++) {
      cout<< " Trying to get "<<Form("%s_%s",profilebasenameEG1,periodNamesGA1[ip])<<endl;
      multEstimatorAvgGA1[ip] = (TProfile*)(fileEstimatorGA1->Get(Form("%s_%s",profilebasenameEG1,periodNamesGA1[ip]))->Clone(Form("%s_%s_clone",profilebasenameEG1,periodNamesGA1[ip])));
      if (!multEstimatorAvgGA1[ip]) {
	Printf("Multiplicity estimator for %s not found! Please check your estimator file",periodNamesGA1[ip]);
	return NULL;
      }
    }

    HFEtaskGA1->SetMultiProfileLHC16s(multEstimatorAvgGA1[0]);
    HFEtaskGA1->SetMultiProfileLHC16r(multEstimatorAvgGA1[1]);
       
       
     TString finDirnameGA1         = "_TrigGA1";
     TString outBasicnameGA1       = "EID";
     TString profnameGA1       = "coutputProf";
       
     finDirnameGA1 	      +=   suffixName.Data();
     outBasicnameGA1      +=   finDirnameGA1.Data();
     profnameGA1          +=   finDirnameGA1.Data();
       
       
     mgr->ConnectInput(HFEtaskGA1,0,mgr->GetCommonInputContainer());
     mgr->ConnectOutput(HFEtaskGA1,1,mgr->CreateContainer(outBasicnameGA1, TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
     mgr->ConnectOutput(HFEtaskGA1,2,mgr->CreateContainer(profnameGA1, TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
     

// EMCal EGA GA2
       
     AliAnalysisTaskHFEMultiplicity* HFEtaskGA2 = new AliAnalysisTaskHFEMultiplicity("");
     mgr->AddTask(HFEtaskGA2);
     HFEtaskGA2->SetReadMC(readMC);
     HFEtaskGA2->SelectCollisionCandidates(AliVEvent::kEMCEGA);
     HFEtaskGA2->SetEMCalTriggerEG2(kTRUE);
     HFEtaskGA2->SetEMCalTriggerDG2(kTRUE);
     HFEtaskGA2->SetTenderSwitch(useTender);
     HFEtaskGA2->SwitchPi0EtaWeightCalc(SwitchPi0EtaWeightCalc);
     HFEtaskGA2->SetClusterTypeEMC(ClsTypeEMC);
     HFEtaskGA2->SetClusterTypeDCAL(ClsTypeDCAL);
     HFEtaskGA2->Setzvtxcut(zvtxcut);
     HFEtaskGA2->SetzvtxQA(zvtxQA);
     HFEtaskGA2->SetNcontVCut(NcontV);
     HFEtaskGA2->SetEtaRange(Etarange);
     HFEtaskGA2->SetNTPCCluster(TPCNclus);
     HFEtaskGA2->SetNITSCluster(ITSNclus);
     HFEtaskGA2->SetTrackpTMin(TrackPtMin);
     HFEtaskGA2->SetDCACut(DCAxyCut,DCAzCut);
     HFEtaskGA2->SetTPCnsigma(TPCNSigMin,TPCNSigMax);
     HFEtaskGA2->SetEopE(EopEMin,EopEMax);
     HFEtaskGA2->SetShowerShapeEM20(M20Min,M20Max);
     HFEtaskGA2->SetInvMassCut(InvmassCut);
     HFEtaskGA2->SetAssoTPCclus(AssoTPCCluster);
     HFEtaskGA2->SetAssoITSclus(AssoITSCluster);
     HFEtaskGA2->SetAssoITSrefit(AssoITSRefit);
     HFEtaskGA2->SetAssopTMin(AssoEPt);
     HFEtaskGA2->SetAssoEtarange(AssoEEtarange);
     HFEtaskGA2->SetAssoTPCnsig(AssoENsigma);

    TFile* fileEstimatorGA2=TFile::Open(estimatorFilename.Data());
    if(!fileEstimatorGA2)  {
      Printf("File with multiplicity estimator not found\n");
      return NULL;
    }
    HFEtaskGA2->SetReferenceMultiplicity(refMult);
    const Char_t* profilebasenameEG2="SPDmultEG2";
       const Char_t* periodNamesGA2[2] = {"LHC16s", "LHC16r"};
    TProfile* multEstimatorAvgGA2[2];
    for(Int_t ip=0; ip<2; ip++) {
      cout<< " Trying to get "<<Form("%s_%s",profilebasenameEG2,periodNamesGA2[ip])<<endl;
      multEstimatorAvgGA2[ip] = (TProfile*)(fileEstimatorGA2->Get(Form("%s_%s",profilebasenameEG2,periodNamesGA2[ip]))->Clone(Form("%s_%s_clone",profilebasenameEG2,periodNamesGA2[ip])));
      if (!multEstimatorAvgGA2[ip]) {
	Printf("Multiplicity estimator for %s not found! Please check your estimator file",periodNamesGA2[ip]);
	return NULL;
      }
    }

    HFEtaskGA2->SetMultiProfileLHC16s(multEstimatorAvgGA2[0]);
    HFEtaskGA2->SetMultiProfileLHC16r(multEstimatorAvgGA2[1]);
       
       
     TString finDirnameGA2         = "_TrigGA2";
     TString outBasicnameGA2       = "EID";
     TString profnameGA2       = "coutputProf";
       
     finDirnameGA2 	      +=   suffixName.Data();
     outBasicnameGA2      +=   finDirnameGA2.Data();
     profnameGA2          +=   finDirnameGA2.Data();
       
       
     mgr->ConnectInput(HFEtaskGA2,0,mgr->GetCommonInputContainer());
     mgr->ConnectOutput(HFEtaskGA2,1,mgr->CreateContainer(outBasicnameGA2, TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
        mgr->ConnectOutput(HFEtaskGA2,2,mgr->CreateContainer(profnameGA2, TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
     
       
    
     
}
    
  
    
  return NULL;
}
