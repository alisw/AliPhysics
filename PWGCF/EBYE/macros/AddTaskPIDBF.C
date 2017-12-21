
//Task for  BF with PID --------------------------------------
//[Noor Alam : sk.noor.alam@cern.ch/noor1989phyalam@gmail.com]
//[Varibale Energy Cyclotron Centre, Kolkata, India]
//------------------------------------------------------------ 

Bool_t kUseNSigmaPID = kTRUE;

//_________________________________________________________//
AliAnalysisTaskPIDBF *AddTaskPIDBF(Double_t centrMin=0.,
	                           Double_t centrMax=100.,
			           Bool_t gRunMixing=kTRUE,
				   Bool_t gRunMixingWithEventPlane=kFALSE,
				   TString centralityEstimator="V0M",
				   Double_t vertexZ=10.,
				   Double_t DCAxy=-1,
				   Double_t DCAz=-1,
				   Double_t ptMin=0.3,
				   Double_t ptMax=1.5,
                                   Double_t kTPCPtMax=0.0,
			           Double_t etaMin=-0.8,
				   Double_t etaMax=0.8,
				   Double_t maxTPCchi2 = -1, 
				   Int_t minNClustersTPC = -1,
				   Bool_t kUsePID = kTRUE,
                                   Int_t ParticleType_ = 0,
                                   Int_t DetectorType_ =0,
                                   Bool_t bResonancesCut = kTRUE,
				   Bool_t bHBTcut = kTRUE,
				   Double_t HBTCutValue = 0.02,
				   Bool_t bConversionCut = kTRUE,
			           Double_t invMassForConversionCut = 0.04,
				   Bool_t bMomentumDifferenceCut = kTRUE,
				   Double_t fQCutMin = 0.0,
				   Int_t AODfilterBit = 128,
				   TString fileNameBase="AnalysisResults",
				   TString dirNameExtra="",
				   TString fArgEventClass="Centrality",
				   TString analysisTypeUser="AOD",
				   Bool_t bVertexBinning=kTRUE,
				   Double_t sigmaElectronRejection=3,
                                   Double_t nSigmaMax = 3.0,
				   Bool_t electronExclusiveRejection=kFALSE,
				   TString correctionFileName = "",
			           Int_t nCentralityArrayBinsForCorrection = -1,
	                           Double_t *gCentralityArrayForCorrections = 0x0,
                                   Bool_t bMomentumOrdering = kTRUE,
                                   Bool_t bCentralTrigger=kFALSE) {
  // Creates a balance function analysis task and adds it to the analysis manager.
  // Get the pointer to the existing analysis manager via the static access method.
  TString outputFileName(fileNameBase);
  outputFileName.Append(".root");

  //===========================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskBF", "No analysis manager to connect to.");
    return NULL;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  //===========================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskBF", "This task requires an input event handler");
    return NULL;
  }
  TString analysisType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  if(dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler())) analysisType = "MC";
  
  // to set the analysis type manually
  if(analysisTypeUser != ""){
    analysisType = analysisTypeUser;
    ::Info("AddTaskBF",Form("Analysis Type manually set to %s",analysisType.Data()));
  }

  gROOT->LoadMacro("$ALICE_PHYSICS/PWGCF/EBYE/macros/InitAnalysisBase.C");
  AliPidBFBase *bf  = 0;  // Balance Function object
  AliPidBFBase *bfm = 0;  // mixing Balance function object

  //maximum Delta eta range
  Double_t deltaEtaMax=TMath::Abs(etaMax-etaMin);

  if (analysisType=="AOD"){
   bf  = GetBaseBF("AOD",centralityEstimator,centrMin,centrMax,bResonancesCut,bHBTcut,HBTCutValue,bConversionCut,invMassForConversionCut,bMomentumDifferenceCut,fQCutMin,bVertexBinning,fArgEventClass,bMomentumOrdering,3,6);
  if(gRunMixing)  bfm= GetBaseBF("AOD",centralityEstimator,centrMin,centrMax,bResonancesCut,bHBTcut,HBTCutValue,bConversionCut,invMassForConversionCut,bMomentumDifferenceCut,fQCutMin,bVertexBinning,fArgEventClass,bMomentumOrdering,3,6);
  }

  else if (analysisType=="MCAOD"){
 bf  = GetBaseBF("MCAOD",centralityEstimator,centrMin,centrMax,bResonancesCut,bHBTcut,HBTCutValue,bConversionCut,invMassForConversionCut,bMomentumDifferenceCut,fQCutMin,bVertexBinning,fArgEventClass,bMomentumOrdering,3,6);
  if(gRunMixing)  bfm= GetBaseBF("MCAOD",centralityEstimator,centrMin,centrMax,bResonancesCut,bHBTcut,HBTCutValue,bConversionCut,invMassForConversionCut,bMomentumDifferenceCut,fQCutMin,bVertexBinning,fArgEventClass,bMomentumOrdering,3,6);
}
  else if (analysisType=="MCAODrec"){
 bf  = GetBaseBF("MCAODrec",centralityEstimator,centrMin,centrMax,bResonancesCut,bHBTcut,HBTCutValue,bConversionCut,invMassForConversionCut,bMomentumDifferenceCut,fQCutMin,bVertexBinning,fArgEventClass,bMomentumOrdering,3,6);
  if(gRunMixing)  bfm= GetBaseBF("MCAODrec",centralityEstimator,centrMin,centrMax,bResonancesCut,bHBTcut,HBTCutValue,bConversionCut,invMassForConversionCut,bMomentumDifferenceCut,fQCutMin,bVertexBinning,fArgEventClass,bMomentumOrdering,3,6);
  
}
 else{
    ::Error("AddTaskBF", "analysis type NOT known.");
    return NULL;
  }
  
  // Task to add manager 
 
  AliAnalysisTaskPIDBF *taskBF = new AliAnalysisTaskPIDBF(Form("TaskBFPsi_%.0f-%.0f_Bit%d_%s%s",centrMin,centrMax,AODfilterBit,centralityEstimator.Data(),dirNameExtra.Data()));
  
  //Event characteristics scheme
  taskBF->SetEventClass(fArgEventClass);
  taskBF->UseOfflineTrigger(bCentralTrigger);
  taskBF->UseMultSelection(bCentralTrigger);
/*
  taskBF->SetRapidityUse(kFALSE);
  taskBF->SetMisMatchTOFProb(.01,kTRUE);*/

// Set Check Pile up : 
//   taskBF->CheckPileUp();

 
   if(fArgEventClass == "Multiplicity") {
    taskBF->SetMultiplicityRange(centrMin,centrMax);
    taskBF->SetMultiplicityEstimator(centralityEstimator);
    //cout<<"Multiplicity estimator "<<centralityEstimator.Data()<<endl;
  }
  else if(fArgEventClass == "Centrality") {
    if(analysisType == "MC")
      taskBF->SetImpactParameterRange(centrMin,centrMax);
    else {
      taskBF->SetCentralityPercentileRange(centrMin,centrMax);
      // centrality estimator (default = V0M)
      taskBF->SetCentralityEstimator(centralityEstimator);
     // cout<<"Centrality estimator "<<centralityEstimator.Data()<<endl;
    }
  }
 

  //  correction factor from MC data : 

   if(correctionFileName != "")
    taskBF->SetInputCorrection(Form("$ALICE_PHYSICS/PWGCF/EBYE/BalanceFunctions/Corrections/%s",correctionFileName.Data()),nCentralityArrayBinsForCorrection,gCentralityArrayForCorrections);

  taskBF->SetAnalysisObject(bf);
  if(gRunMixing){
    taskBF->SetMixingObject(bfm);
    taskBF->SetMixingTracks(50000);
    if(gRunMixingWithEventPlane){
      taskBF->SetMixingWithEventPlane(gRunMixingWithEventPlane);
    }
  }

 
  if(analysisType == "AOD") {
    taskBF->SetAODtrackCutBit(AODfilterBit);
    taskBF->SetKinematicsCutsAOD(ptMin,ptMax,etaMin,etaMax);

    // set extra DCA cuts (-1 no extra cut)
    taskBF->SetExtraDCACutsAOD(DCAxy,DCAz);

    // set extra TPC chi2 / nr of clusters cut
    taskBF->SetExtraTPCCutsAOD(maxTPCchi2, minNClustersTPC);

    // electron rejection (so far only for AOD), <0 --> no rejection
    if(sigmaElectronRejection > 0){
      if(electronExclusiveRejection) taskBF->SetElectronOnlyRejection(sigmaElectronRejection); // no other particle in nsigma 
      else                           taskBF->SetElectronRejection(sigmaElectronRejection); // check only if electrons in nsigma
    }

    //++++++++++++++++//
    if(kUsePID) {
       if(kUseNSigmaPID){
 	taskBF->SetUseNSigmaPID(nSigmaMax);
        taskBF->SetParticleType(ParticleType_);   // Added  Noor Alam
        taskBF->SetDetectorPID(DetectorType_);   // Added  Noor Alam
        taskBF->SetTPCPtMinMax(ptMin,kTPCPtMax); // TPC Pt min and max (N.Alam )
        taskBF->SetTOFPtMinMax(kTPCPtMax,ptMax); // TOF Pt min and max (N.Alam)
      //   N.A taskBF->SetParticleOfInterest(AliAnalysisTaskPIDBF::kKaon);
     // taskBF->SetDetectorUsedForPID(AliAnalysisTaskPIDBF::kTPCTOF); //TOFpid,TPCpid
    }
}
    //++++++++++++++++//

  }
 

 
  else if(analysisType == "MCAOD") {
    // pt and eta cut (pt_min, pt_max, eta_min, eta_max)
    taskBF->SetAODtrackCutBit(AODfilterBit);
    taskBF->SetKinematicsCutsAOD(ptMin,ptMax,etaMin,etaMax);
    taskBF->ExcludeElectronsInMC();
   if(kUsePID) {
    taskBF->SetUseNSigmaPID(nSigmaMax);
    taskBF->SetParticleType(ParticleType_);
    }   
  }

  else if(analysisType == "MCAODrec") {     //++++++++++++++++
    taskBF->SetAODtrackCutBit(AODfilterBit);
    taskBF->SetKinematicsCutsAOD(ptMin,ptMax,etaMin,etaMax); 
    taskBF->ExcludeElectronsInMC();

    // set extra DCA cuts (-1 no extra cut)
    taskBF->SetExtraDCACutsAOD(DCAxy,DCAz);

    // set extra TPC chi2 / nr of clusters cut
    taskBF->SetExtraTPCCutsAOD(maxTPCchi2, minNClustersTPC);

    // electron rejection (so far only for AOD), <0 --> no rejection
    if(sigmaElectronRejection > 0){
      if(electronExclusiveRejection) taskBF->SetElectronOnlyRejection(sigmaElectronRejection); // no other particle in nsigma 
      else                           taskBF->SetElectronRejection(sigmaElectronRejection); // check only if electrons in nsigma
    }
 
    if(kUsePID) {
       if(kUseNSigmaPID){
        taskBF->SetUseNSigmaPID(nSigmaMax);
        taskBF->SetParticleType(ParticleType_);   // Added  Noor Alam
        taskBF->SetDetectorPID(DetectorType_);   // Added  Noor Alam
        taskBF->SetTPCPtMinMax(ptMin,kTPCPtMax); // TPC Pt min and max (N.Alam )
        taskBF->SetTOFPtMinMax(kTPCPtMax,ptMax); // TOF Pt min and max (N.Alam)
    }
}

  }//++++++++++++++++

     if(bCentralTrigger) taskBF->SelectCollisionCandidates(AliVEvent::kINT7);

     else taskBF->SelectCollisionCandidates(AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral);



  // centrality estimator (default = V0M)
  taskBF->SetCentralityEstimator(centralityEstimator);
  
  // vertex cut (x,y,z)
  taskBF->SetVertexDiamond(3.,3.,vertexZ);

  mgr->AddTask(taskBF);
  
  //==============================================================================
  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  outputFileName += ":PWGCFEbyE.outputBalanceFunctionPsiAnalysis";
  AliAnalysisDataContainer *coutQA = mgr->CreateContainer(Form("listQAPsi_%.0f-%.0f_Bit%d_%s%s",centrMin,centrMax,AODfilterBit,centralityEstimator.Data(),dirNameExtra.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName.Data());
  AliAnalysisDataContainer *coutBF = mgr->CreateContainer(Form("listBFPsi_%.0f-%.0f_Bit%d_%s%s",centrMin,centrMax,AODfilterBit,centralityEstimator.Data(),dirNameExtra.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName.Data());
  if(gRunMixing) AliAnalysisDataContainer *coutBFM = mgr->CreateContainer(Form("listBFPsiMixed_%.0f-%.0f_Bit%d_%s%s",centrMin,centrMax,AODfilterBit,centralityEstimator.Data(),dirNameExtra.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName.Data());
  if(kUsePID || sigmaElectronRejection > 0) AliAnalysisDataContainer *coutQAPID = mgr->CreateContainer(Form("listQAPIDPsi_%.0f-%.0f_Bit%d_%s%s",centrMin,centrMax,AODfilterBit,centralityEstimator.Data(),dirNameExtra.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName.Data());

   AliAnalysisDataContainer *coutQA_AliEventCuts = mgr->CreateContainer(Form("listQA_AliEventCuts_%.0f-%.0f_Bit%d_%s%s",centrMin,centrMax,AODfilterBit,centralityEstimator.Data(),dirNameExtra.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName.Data());


  mgr->ConnectInput(taskBF, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskBF, 1, coutQA);
  mgr->ConnectOutput(taskBF, 2, coutBF);
  if(gRunMixing) mgr->ConnectOutput(taskBF, 3, coutBFM);
  if(kUsePID||sigmaElectronRejection > 0) mgr->ConnectOutput(taskBF, 4, coutQAPID);
  mgr->ConnectOutput(taskBF, 5, coutQA_AliEventCuts);
  //if((kUsePID && analysisType == "AOD")||sigmaElectronRejection > 0) mgr->ConnectOutput(taskBF, 5, coutQAPID);
  //if((kUsePID && analysisType == "ESD")||sigmaElectronRejection > 0) mgr->ConnectOutput(taskBF, 5, coutQAPID);

  return taskBF;
}
