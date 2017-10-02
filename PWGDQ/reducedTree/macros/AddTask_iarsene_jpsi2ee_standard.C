//void Setup(AliReducedAnalysisJpsi2ee* processor, TString prod="LHC10h");
//AliHistogramManager* SetupHistogramManager(TString prod="LHC10h");
//void DefineHistograms(AliHistogramManager* man, TString prod="LHC10h");

//__________________________________________________________________________________________
AliAnalysisTask* AddTask_iarsene_jpsi2ee_standard(Bool_t isAliRoot=kTRUE, Int_t runMode=1, TString prod="LHC10h"){    
   //
   // isAliRoot=kTRUE for ESD/AOD analysis in AliROOT, kFALSE for root analysis on reduced trees
   // runMode=1 (AliAnalysisTaskReducedEventProcessor::kUseOnTheFlyReducedEvents)
   //               =2 (AliAnalysisTaskReducedEventProcessor::kUseEventsFromTree)
   //
   //get the current analysis manager

  printf("INFO on AddTask_iarsene_jpsi2ee(): (isAliRoot, runMode) :: (%d,%d) \n", isAliRoot, runMode);

  AliReducedAnalysisJpsi2ee* jpsi2eeAnalysis = new AliReducedAnalysisJpsi2ee("Jpsi2eeAnalysis","Jpsi->ee analysis");
  jpsi2eeAnalysis->Init();
  //jpsi2eeAnalysis->SetLoopOverTracks(kFALSE);
  //jpsi2eeAnalysis->SetRunEventMixing(kFALSE);
  //jpsi2eeAnalysis->SetRunPairing(kFALSE);
  //jpsi2eeAnalysis->SetRunOverMC(kTRUE);
  //jpsi2eeAnalysis->SetRunLikeSignPairing(kFALSE);
  Setup(jpsi2eeAnalysis, prod);
  // initialize an AliAnalysisTask which will wrapp the AliReducedAnalysisJpsi2ee such that it can be run in an aliroot analysis train (e.g. LEGO, local analysis etc)
  AliAnalysisTaskReducedEventProcessor* task = new AliAnalysisTaskReducedEventProcessor("ReducedEventAnalysisManager", runMode);
  task->AddTask(jpsi2eeAnalysis);
  
  if(isAliRoot){
     AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
     if (!mgr) {
        Error("AddTask_iarsene_dst", "No analysis manager found.");
        return 0;
     }
     
     AliAnalysisDataContainer* cReducedEvent = NULL;
     if(runMode==AliAnalysisTaskReducedEventProcessor::kUseOnTheFlyReducedEvents) {
       printf("INFO on AddTask_iarsene_jpsi2ee(): use on the fly events ");
       cReducedEvent = (AliAnalysisDataContainer*)mgr->GetContainers()->FindObject("ReducedEventDQ");
       if(!cReducedEvent) {
         printf("ERROR: In AddTask_iarsene_jpsi2ee(), couldn't find exchange container with ReducedEvent! ");
         printf("             You need to use AddTask_iarsene_dst() in the on-the-fly reduced events mode!");
         return 0x0;
       }
     }
            
     mgr->AddTask(task);
      
     if(runMode==AliAnalysisTaskReducedEventProcessor::kUseEventsFromTree) 
        mgr->ConnectInput(task,  0, mgr->GetCommonInputContainer());
      
     if(runMode==AliAnalysisTaskReducedEventProcessor::kUseOnTheFlyReducedEvents) 
       mgr->ConnectInput(task, 0, cReducedEvent);
  
     AliAnalysisDataContainer *cOutputHist = mgr->CreateContainer("jpsi2eeHistos", THashList::Class(),
                                                                  AliAnalysisManager::kOutputContainer, "dstAnalysisHistograms.root");
     mgr->ConnectOutput(task, 1, cOutputHist );
  }
  else {
    // nothing at the moment   
  }
  
  return task;
}


//_________________________________________________________________
void Setup(AliReducedAnalysisJpsi2ee* processor, TString prod /*="LHC10h"*/) {
  //
  // Configure the analysis task
  // Setup histograms, handlers, cuts, etc.
  //
  if(prod.Contains("LHC15o_highIR") && !processor->GetRunOverMC()) {
    TFile* corrFile = TFile::Open("/home/iarsene/work/ALICE/treeAnalysis/newdst/development/test/CorrMaps_CentralityNew_LHC15o_pass1.root");
    AliReducedVarManager::SetTPCelectronCorrectionMaps((TH2F*)corrFile->Get("CorrMapCentroid"), (TH2F*)corrFile->Get("CorrMapWidth"), AliReducedVarManager::kEta, AliReducedVarManager::kCentVZERO);
    corrFile->Close();
    AliReducedVarManager::SetVZEROCalibrationPath("/home/iarsene/data/2015/LHC15o/DQ_PbPb_LEGO287/analysisOutputs/2017_01_09_VZEROepRecentering/");
    AliReducedVarManager::SetCalibrateVZEROqVector(kTRUE);
    AliReducedVarManager::SetRecenterVZEROqVector(kTRUE);
  }
  if((prod.Contains("LHC15o_pidfix") || prod.Contains("LHC15o_lowIR")) && !processor->GetRunOverMC()) {
     TFile* corrFile = TFile::Open("/home/iarsene/work/ALICE/treeAnalysis/newdst/development/test/CorrMaps_CentralityNew_LHC15o_pass1_pidfix.root");
     AliReducedVarManager::SetTPCelectronCorrectionMaps((TH2F*)corrFile->Get("CorrMapCentroid"), (TH2F*)corrFile->Get("CorrMapWidth"), AliReducedVarManager::kEta, AliReducedVarManager::kCentVZERO);
     corrFile->Close();
     if(prod.Contains("LHC15o_pidfix")) {
       AliReducedVarManager::SetVZEROCalibrationPath("/home/iarsene/data/2015/LHC15o/DQ_PbPb_LEGO288_pidfix/analysisOutputs/2017_01_09_VZEROepRecentering/");
       AliReducedVarManager::SetCalibrateVZEROqVector(kTRUE);
       AliReducedVarManager::SetRecenterVZEROqVector(kTRUE);
     }
  }
  
  if(prod.Contains("LHC15o")) {
    TFile* grpFile = TFile::Open("/home/iarsene/work/ALICE/treeAnalysis/newdst/development/test/grpData_LHC15o_allRuns.root");
    AliReducedVarManager::SetLHCDataInfo((TH1F*)grpFile->Get("lumiTotal"), (TH1F*)grpFile->Get("intTotal0"), 
                                       (TH1F*)grpFile->Get("intTotal1"), (TH1I*)grpFile->Get("fillNumber"));
    AliReducedVarManager::SetGRPDataInfo((TH1I*)grpFile->Get("dipolePolarity"), (TH1I*)grpFile->Get("l3Polarity"), 
                                       (TH1I*)grpFile->Get("timeStart"), (TH1I*)grpFile->Get("timeStop"));
    grpFile->Close();
  }
  
  
  // Set event cuts
  AliReducedEventCut* evCut1 = new AliReducedEventCut("Centrality","Centrality selection");
  //if(prod.Contains("LHC15o")) evCut1->AddCut(AliReducedVarManager::kCentVZERO, 30., 90.);
  evCut1->AddCut(AliReducedVarManager::kVZEROTotalMult, 0.0, 37500.0);
  evCut1->AddCut(AliReducedVarManager::kVtxZ, -10.0, 10.0);
  evCut1->AddCut(AliReducedVarManager::kIsPhysicsSelection, 0.1, 2.);   // request physics selection
  //evCut1->AddCut(AliReducedVarManager::kDeltaVtxZ, -0.2, 0.2);   // select based on the difference between the vtxZ and the tPC vtxZ
  //evCut1->AddCut(AliReducedVarManager::kIRIntClosestIntMap+1, 0.99, 5000., kTRUE);   // exclude out-of-bunch pileup
  //evCut1->EnableVertexDistanceCut();
  //evCut1->AddCut(AliReducedVarManager::kTZEROpileup, -0.1, 0.1);
  TF1* cutCorrTPCoutVZEROmult = new TF1("cutCorrTPCoutVZEROmult", "[0]+[1]*x+[2]*x*x", 0., 1.e+5);
  cutCorrTPCoutVZEROmult->SetParameters(-2200., 2.5, 1.2e-5);
  if(prod.Contains("LHC15o") && !processor->GetRunOverMC()) 
      evCut1->AddCut(AliReducedVarManager::kVZEROTotalMult, cutCorrTPCoutVZEROmult, 99999., kFALSE, AliReducedVarManager::kNTracksPerTrackingStatus+AliReducedVarManager::kTPCout, 0., 99998.);
  
  TF1* cutCorrVZEROmultVtxContr = new TF1("cutCorrVZEROmultVtxContr", "[0]+[1]*x", 0.0, 40000.);
  cutCorrVZEROmultVtxContr->SetParameters(-150., 0.09);
  if(prod.Contains("LHC15o") && !processor->GetRunOverMC()) 
     evCut1->AddCut(AliReducedVarManager::kNVtxContributors, cutCorrVZEROmultVtxContr, 99999., kFALSE, AliReducedVarManager::kVZEROTotalMult, 0., 99998.);
  
  
  processor->AddEventCut(evCut1);
  
  // Set track cuts
  AliReducedTrackCut* standardCut = new AliReducedTrackCut("standard","");
  standardCut->AddCut(AliReducedVarManager::kPt, 1.0,50.0);
  standardCut->AddCut(AliReducedVarManager::kEta, -0.9,0.9);
  standardCut->AddCut(AliReducedVarManager::kDcaXY, -1.0,1.0);
  standardCut->AddCut(AliReducedVarManager::kDcaZ, -3.0,3.0);
  standardCut->AddCut(AliReducedVarManager::kTPCncls, 70.,160.0);
  if(prod.Contains("LHC15o_highIR")) {
     AddCentralityDependentProtRej(standardCut, kFALSE, processor->GetRunOverMC(), -0.50);
     AddCentralityDependentPionRej(standardCut, kFALSE, processor->GetRunOverMC(), -1.50);
  }
  if(prod.Contains("LHC15o_pidfix") || prod.Contains("LHC15o_lowIR")) {
     AddCentralityDependentProtRej(standardCut, kTRUE, processor->GetRunOverMC(), -0.50);
     AddCentralityDependentPionRej(standardCut, kTRUE, processor->GetRunOverMC(), -1.50);
  }
  if(!processor->GetRunOverMC())
    standardCut->AddCut(AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron, -3.0, 3.0);
  else
    standardCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron, -3.0, 3.0); 
  //standardCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron, -3.0, 3.0);
  //standardCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kProton, 3.0, 30000.0);
  //standardCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kPion, 3.0, 30000.0);
  standardCut->SetRejectKinks();
  standardCut->SetRequestITSrefit();
  standardCut->SetRequestTPCrefit();
  //standardCut->SetRequestTOFout();
  standardCut->SetRequestSPDany();
  TF1* chi2Cut=new TF1("chi2Cut","[0]+[1]*x",0.,15000.);
  chi2Cut->SetParameters(1.9, 1.1e-4);
  //standardCut->AddCut(AliReducedVarManager::kTPCchi2, 0., chi2Cut, kFALSE, AliReducedVarManager::kNTracksPerTrackingStatus+AliReducedVarManager::kTPCout, 0., 99998.);
  standardCut->AddCut(AliReducedVarManager::kTPCchi2, 0., 4.0);
  standardCut->AddCut(AliReducedVarManager::kITSchi2, 0., 10.);
  standardCut->AddCut(AliReducedVarManager::kTPCNclusBitsFired, 6., 9.);
  standardCut->AddCut(AliReducedVarManager::kTPCnclsSharedRatio, 0.3, 2., kTRUE);
  standardCut->AddCut(AliReducedVarManager::kTPCcrossedRowsOverFindableClusters, 0.8, 2.);
  //standardCut->AddCut(AliReducedVarManager::kChi2TPCConstrainedVsGlobal, 36., 1.0e+10, kTRUE);
  processor->AddTrackCut(standardCut); 
  
  
  AliReducedTrackCut* standardCutNoKink = new AliReducedTrackCut("standardNoKinkCut","");
  standardCutNoKink->AddCut(AliReducedVarManager::kPt, 1.0,50.0);
  standardCutNoKink->AddCut(AliReducedVarManager::kEta, -0.9,0.9);
  standardCutNoKink->AddCut(AliReducedVarManager::kDcaXY, -1.0,1.0);
  standardCutNoKink->AddCut(AliReducedVarManager::kDcaZ, -3.0,3.0);
  standardCutNoKink->AddCut(AliReducedVarManager::kTPCncls, 70.,160.0);
  if(prod.Contains("LHC15o_highIR")) {
     AddCentralityDependentProtRej(standardCutNoKink, kFALSE, processor->GetRunOverMC(), -0.50);
     AddCentralityDependentPionRej(standardCutNoKink, kFALSE, processor->GetRunOverMC(), -1.50);
  }
  if(prod.Contains("LHC15o_pidfix") || prod.Contains("LHC15o_lowIR")) {
     AddCentralityDependentProtRej(standardCutNoKink, kTRUE, processor->GetRunOverMC(), -0.50);
     AddCentralityDependentPionRej(standardCutNoKink, kTRUE, processor->GetRunOverMC(), -1.50);
  }
  if(!processor->GetRunOverMC())
     standardCutNoKink->AddCut(AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron, -3.0, 3.0);
  else
     standardCutNoKink->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron, -3.0, 3.0); 
  //standardCutNoKink->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron, -3.0, 3.0);
  //standardCutNoKink->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kProton, 3.0, 30000.0);
  //standardCutNoKink->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kPion, 3.0, 30000.0);
  //standardCutNoKink->SetRejectKinks();
  standardCutNoKink->SetRequestITSrefit();
  standardCutNoKink->SetRequestTPCrefit();
  //standardCutNoKink->SetRequestTOFout();
  standardCutNoKink->SetRequestSPDany();
  TF1* chi2Cut=new TF1("chi2Cut","[0]+[1]*x",0.,15000.);
  chi2Cut->SetParameters(1.9, 1.1e-4);
  //standardCutNoKink->AddCut(AliReducedVarManager::kTPCchi2, 0., chi2Cut, kFALSE, AliReducedVarManager::kNTracksPerTrackingStatus+AliReducedVarManager::kTPCout, 0., 99998.);
  standardCutNoKink->AddCut(AliReducedVarManager::kTPCchi2, 0., 4.0);
  standardCutNoKink->AddCut(AliReducedVarManager::kITSchi2, 0., 10.);
  standardCutNoKink->AddCut(AliReducedVarManager::kTPCNclusBitsFired, 6., 9.);
  standardCutNoKink->AddCut(AliReducedVarManager::kTPCnclsSharedRatio, 0.3, 2., kTRUE);
  standardCutNoKink->AddCut(AliReducedVarManager::kTPCcrossedRowsOverFindableClusters, 0.8, 2.);
  //standardCutNoKink->AddCut(AliReducedVarManager::kChi2TPCConstrainedVsGlobal, 36., 1.0e+10, kTRUE);
  processor->AddTrackCut(standardCutNoKink); 
    
  // set track prefilter cuts
  AliReducedTrackCut* prefTrackCut1 = new AliReducedTrackCut("prefCutPt09","prefilter P selection");
  prefTrackCut1->AddCut(AliReducedVarManager::kP, 10.7,100.0);
  prefTrackCut1->SetRequestTPCrefit();
  processor->AddPrefilterTrackCut(prefTrackCut1);  
  
  // set pair prefilter cuts
  AliReducedVarCut* prefPairCut = new AliReducedVarCut("prefCutM50MeV","prefilter pair cuts");
  prefPairCut->AddCut(AliReducedVarManager::kMass, 0.0, 0.05, kTRUE);
  processor->AddPrefilterPairCut(prefPairCut);
  
  // Set pair cuts
  AliReducedTrackCut* pairCut1 = new AliReducedTrackCut("Ptpair","Pt pair selection");
  pairCut1->AddCut(AliReducedVarManager::kPt, 0.0,100.0);
  pairCut1->AddCut(AliReducedVarManager::kRap, -0.9,0.9);
  processor->AddPairCut(pairCut1);
  
  //SetupHistogramManager(processor->GetHistogramManager(), prod);
  SetupHistogramManager(processor, prod);
  SetupMixingHandler(processor);
}

//______________________________________________________________________________________
void AddCentralityDependentProtRej(AliReducedTrackCut* cut, Bool_t isPidFix=kFALSE, Bool_t isMC=kFALSE, Double_t nSigOffset=0.0) {
   //
   // setup the centrality dependent proton rejection
   //
   Double_t centMin[8] = {0., 5., 10., 15., 20., 30., 40., 60.};
   Double_t centMax[8] = {5., 10., 15., 20., 30., 40., 60., 90.};
   TF1* fprot[8];
   for(Int_t i=0;i<8;++i) {
      fprot[i] = new TF1(Form("fprot%d",i), "pol3", 0.,50.);
      if(isPidFix) 
         SetupProtonRejFunctionPidfix(fprot[i],i, nSigOffset);
      else       
         SetupProtonRejFunction(fprot[i],i, nSigOffset);
      if(isMC)
         cut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron, fprot[i], 1.0e+30, kFALSE, AliReducedVarManager::kPin, 0., 3.);
      else
         cut->AddCut(AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron, fprot[i], 1.0e+30, kFALSE, AliReducedVarManager::kPin, 0., 3.);
   }
}

//______________________________________________________________________________________
void AddCentralityDependentPionRej(AliReducedTrackCut* cut, Bool_t isPidFix=kFALSE, Bool_t isMC=kFALSE, Double_t nSigOffset=0.0) {
   //
   // setup the centrality dependent proton rejection
   //
   Double_t centMin[8] = {0., 5., 10., 15., 20., 30., 40., 60.};
   Double_t centMax[8] = {5., 10., 15., 20., 30., 40., 60., 90.};
   TF1* fpion[8];
   for(Int_t i=0;i<8;++i) {
      fpion[i] = new TF1(Form("fpion%d",i), "pol3", 0.,50.);
      if(isPidFix)
         SetupPionRejFunctionPidfix(fpion[i],i, nSigOffset);
      else
         SetupPionRejFunction(fpion[i],i, nSigOffset);
      if(isMC)
         cut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron, fpion[i], 1.0e+30, kFALSE, AliReducedVarManager::kPin, 1.5, 20.);
      else
         cut->AddCut(AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron, fpion[i], 1.0e+30, kFALSE, AliReducedVarManager::kPin, 1.5, 20.);
   }
}


//______________________________________________________________________________________
void SetupProtonRejFunction(TF1* func, Int_t centBin, Double_t nSigOffset=0.0) {
   //
   // setup the proton rejection function
   // Here use the numbers which Dennis provided for 3.5 sigma proton rejection
   //
   switch(centBin) {
      case 0:                     // cent 0-5%
         func->SetParameters(49.4868+nSigOffset,-85.1973,49.0375,-9.88178);
         break;
      case 1:                     // cent 5-10%
         func->SetParameters(44.6979+nSigOffset,-74.2031,40.815,-7.89076);
         break;
      case 2:                     // cent 10-15%
         func->SetParameters(48.6892+nSigOffset,-81.7643,45.4575,-8.83387);
         break;
      case 3:                     // cent 15-20%
         func->SetParameters(53.3788+nSigOffset,-92.5709,53.6544,-10.8983);
         break;
      case 4:                     // cent 20-30%
         func->SetParameters(50.1662+nSigOffset,-84.1989,46.6445,-9.03298);
         break;
      case 5:                     // cent 30-40%
         func->SetParameters(61.5201+nSigOffset,-109.668,65.4643,-13.6647);
         break;
      case 6:                     // cent 40-60%
         func->SetParameters(66.5285+nSigOffset,-120.471,73.0039,-15.4409);
         break;
      case 7:                     // cent 60-90%
         func->SetParameters(60.668+nSigOffset,-108.065,64.4816,-13.6468);
         break;
   };
}

//______________________________________________________________________________________
void SetupProtonRejFunctionPidfix(TF1* func, Int_t centBin, Double_t nSigOffset=0.0) {
   //
   // setup the proton rejection function
   // Here use the numbers which Dennis provided for 3.5 sigma proton rejection
   //
   switch(centBin) {
      case 0:                     // cent 0-5%
         func->SetParameters(46.7265+nSigOffset,-81.1999,47.4381,-9.80234);
         break;
      case 1:                     // cent 5-10%
         func->SetParameters(48.5281+nSigOffset,-84.0613,48.8339,-10.0214);
         break;
      case 2:                     // cent 10-15%
         func->SetParameters(46.9739+nSigOffset,-79.9369,45.2782,-9.04988);
         break;
      case 3:                     // cent 15-20%
         func->SetParameters(53.7956+nSigOffset,-93.7956,54.5638,-11.1288);
         break;
      case 4:                     // cent 20-30%
         func->SetParameters(57.0128+nSigOffset,-100.898,59.774,-12.4282);
         break;
      case 5:                     // cent 30-40%
         func->SetParameters(56.1116+nSigOffset,-97.7581,56.5004,-11.4338);
         break;
      case 6:                     // cent 40-60%
         func->SetParameters(62.3849+nSigOffset,-112.061,67.2837,-14.2049);
         break;
      case 7:                     // cent 60-90%
         func->SetParameters(80.1836+nSigOffset,-154.34,100.532,-22.9687);
         break;
   };
}


//______________________________________________________________________________________
void SetupPionRejFunction(TF1* func, Int_t centBin, Double_t nSigOffset=0.0) {
   //
   // setup the pion rejection function
   // Here use the numbers which Dennis provided for 3.5 sigma pion rejection
   //
   switch(centBin) {
      case 0:                     // cent 0-5%
         func->SetParameters(-4.39039+nSigOffset,1.41342,-0.149678,0.00622121);
         break;
      case 1:                     // cent 5-10%
         func->SetParameters(-4.51303+nSigOffset,1.43655,-0.151552,0.00631887);
         break;
      case 2:                     // cent 10-15%
         func->SetParameters(-4.59046+nSigOffset,1.4011,-0.139758,0.00556735);
         break;
      case 3:                     // cent 15-20%
         func->SetParameters(-4.59451+nSigOffset,1.35289,-0.127082,0.00473946);
         break;
      case 4:                     // cent 20-30%
         func->SetParameters(-4.88737+nSigOffset,1.44344,-0.140034,0.00539735);
         break;
      case 5:                     // cent 30-40%
         func->SetParameters(-5.39416+nSigOffset,1.63295,-0.170718,0.00705499);
         break;
      case 6:                     // cent 40-60%
         func->SetParameters(-5.99412+nSigOffset,1.78627,-0.189964,0.00792811);
         break;
      case 7:                     // cent 60-90%
         func->SetParameters(-6.48072+nSigOffset,1.88094,-0.200417,0.00835834);
         break;
   };
}


//______________________________________________________________________________________
void SetupPionRejFunctionPidfix(TF1* func, Int_t centBin, Double_t nSigOffset=0.0) {
   //
   // setup the pion rejection function
   // Here use the numbers which Dennis provided for 3.5 sigma pion rejection
   //
   switch(centBin) {
      case 0:                     // cent 0-5%
         func->SetParameters(-4.39414+nSigOffset,1.30979,-0.123276,0.00446862);
         break;
      case 1:                     // cent 5-10%
         func->SetParameters(-4.82105+nSigOffset,1.53861,-0.167065,0.00708157);
         break;
      case 2:                     // cent 10-15%
         func->SetParameters(-4.98141+nSigOffset,1.5569,-0.167242,0.00710832);
         break;
      case 3:                     // cent 15-20%
         func->SetParameters(-5.09913+nSigOffset,1.56274,-0.165177,0.00698126);
         break;
      case 4:                     // cent 20-30%
         func->SetParameters(-5.30148+nSigOffset,1.61432,-0.17017,0.00705404);
         break;
      case 5:                     // cent 30-40%
         func->SetParameters(-5.68263+nSigOffset,1.68404,-0.175564,0.00720695);
         break;
      case 6:                     // cent 40-60%
         func->SetParameters(-6.20825+nSigOffset,1.79933,-0.189546,0.00780421);
         break;
      case 7:                     // cent 60-90%
         func->SetParameters(-6.20573+nSigOffset,1.54408,-0.132792,0.00458051);
         break;
   };
}


//_________________________________________________________________
void SetupMixingHandler(AliReducedAnalysisJpsi2ee* task) {
   //
   // setup the mixing handler
   //
   AliMixingHandler* handler = task->GetMixingHandler();
   handler->SetPoolDepth(50);
   handler->SetMixingThreshold(1.0);
   handler->SetDownscaleEvents(1);
   handler->SetDownscaleTracks(1);
   //handler->SetEventVariables(AliReducedVarManager::kCentVZERO,AliReducedVarManager::kTimeRelativeSOR,AliReducedVarManager::kVZERORP+6+1);
   //handler->SetEventVariables(AliReducedVarManager::kCentVZERO,AliReducedVarManager::kVtxZ,AliReducedVarManager::kNTracksTPCoutFromPileup);
   //handler->SetEventVariables(AliReducedVarManager::kCentVZERO, AliReducedVarManager::kVtxZ, AliReducedVarManager::kVZERORP+6+1);
   //handler->SetEventVariables(AliReducedVarManager::kCentVZERO,AliReducedVarManager::kVtxZ,AliReducedVarManager::kNTracksPerTrackingStatus+AliReducedVarManager::kTPCout);
   handler->SetEventVariables(AliReducedVarManager::kCentVZERO, AliReducedVarManager::kNTracksPerTrackingStatus+AliReducedVarManager::kTPCout, AliReducedVarManager::kVZERORP+6+1);
   
   Float_t centLims[27] = {
     0.0, 1.0, 2.0, 3.0, 4.0, 
     5.0, 6.0, 7.0, 8.0, 9.0, 
     10., 12.0, 14.0, 16.0, 18.0, 
     20., 22.5, 25.0, 27.5, 30.0, 
     35., 40., 50., 60., 70., 
     80., 90.
   };
   Float_t ntpcOutLims[26] = {
      0.0, 1000., 2000., 3000., 4000., 
      5000., 6000., 7000., 8000., 8500., 
      9000., 9500., 10000., 10500., 11000., 
      11500., 12000., 12250., 12500., 12750., 
      13000., 13250., 13500., 13750., 14000., 
      14250.
   };
   
   Float_t ntpcOutLimsMC[2];
   ntpcOutLimsMC[0] = 0.0; ntpcOutLimsMC[1] = 2500.0; 
   //Float_t zLims[14] = {-10.,-9.,-8.,-7.,-5.,-3.,-1.,1.,3.,5.,7.,8.,9.,10.};
   Float_t zLims[2] = {-10.,10.};

   Float_t timeLims[2];
   for(Int_t i=0;i<2;++i) timeLims[i] = 0.0+i*450.;
   
   Float_t epLims[11];
   //epLims[0] = -0.5*TMath::Pi();
   //epLims[1] = +0.5*TMath::Pi();
   for(Int_t i=0;i<=10;++i) epLims[i] = -0.5*TMath::Pi()+i*TMath::Pi()/10.;
   
   handler->SetCentralityLimits(17,centLims);
   handler->SetEventPlaneLimits(11,epLims);
   handler->SetEventVertexLimits(16,ntpcOutLims);
   /*if(task->GetRunOverMC())
      handler->SetEventVertexLimits(2,ntpcOutLimsMC);
   else
     handler->SetEventVertexLimits(2,zLims);*/
}


//_________________________________________________________________
void SetupHistogramManager(AliReducedAnalysisJpsi2ee* task, TString prod /*="LHC10h"*/) {
  //
  // setup the histograms manager
  //
  AliReducedVarManager::SetDefaultVarNames();
  
  DefineHistograms(task, prod);
  
  AliReducedVarManager::SetUseVars(task->GetHistogramManager()->GetUsedVars());
}


//_________________________________________________________________
void DefineHistograms(AliReducedAnalysisJpsi2ee* task, TString prod /*="LHC10h"*/) {
  //
  // define histograms
  // NOTE: The DefineHistograms need to be called after the track cuts are defined because the name of the cuts
  //           are used in the histogram lists
   // NOTE: The name of the pair histogram lists for event mixing need to contain "PairMEPP", "PairMEPM", "PairMEMM", in their name, see below.
   //  TODO: make needed changes such that this becomes less prone to mistakes
   
  AliHistogramManager* man = task->GetHistogramManager(); 
   
  TString histClasses = "";
  histClasses += "Event_BeforeCuts;";
  histClasses += "Event_AfterCuts;";
  if(!task->GetRunOverMC()) {
    histClasses += "EventTag_BeforeCuts;";   
    histClasses += "EventTag_AfterCuts;";   
    histClasses += "EventTriggers_BeforeCuts;";
    histClasses += "EventTriggers_AfterCuts;";   
  }
  if(task->GetLoopOverTracks())
    histClasses += "Track_BeforeCuts;";
  if(!task->GetRunOverMC()) {
    histClasses += "TrackStatusFlags_BeforeCuts;";
  }
  if(task->GetLoopOverTracks()) {
    histClasses += "TrackITSclusterMap_BeforeCuts;";
    histClasses += "TrackTPCclusterMap_BeforeCuts;";
  }
  if(task->GetRunOverMC()) {
    histClasses += "MCTruth_BeforeSelection;";
    histClasses += "MCTruth_AfterSelection;";
  }
  if(task->GetLoopOverTracks()) {
    for(Int_t i=0; i<task->GetNTrackCuts(); ++i) {
      TString cutName = task->GetTrackCutName(i);
      histClasses += Form("Track_%s;", cutName.Data());
      if(!task->GetRunOverMC()) {
        histClasses += Form("TrackStatusFlags_%s;", cutName.Data());
      }
      histClasses += Form("TrackITSclusterMap_%s;", cutName.Data());
      histClasses += Form("TrackTPCclusterMap_%s;", cutName.Data());
      if(task->GetRunOverMC()) {
        histClasses += Form("Track_%s_MCTruth;", cutName.Data());
        histClasses += Form("TrackStatusFlags_%s_MCTruth;", cutName.Data());
        histClasses += Form("TrackITSclusterMap_%s_MCTruth;", cutName.Data());
        histClasses += Form("TrackTPCclusterMap_%s_MCTruth;", cutName.Data());
      }
      if(!task->GetRunLikeSignPairing())
         histClasses += Form("PairSEPM_%s;", cutName.Data());
      else
         histClasses += Form("PairSEPP_%s;PairSEPM_%s;PairSEMM_%s;", cutName.Data(), cutName.Data(), cutName.Data());
      if(!task->GetRunOverMC())
        histClasses += Form("PairMEPP_%s;PairMEPM_%s;PairMEMM_%s;", cutName.Data(), cutName.Data(), cutName.Data());
      if(task->GetRunOverMC())
        histClasses += Form("PairSEPM_%s_MCTruth;", cutName.Data());
    }
  }
  
  Int_t runNBins = 0;
  Double_t runHistRange[2] = {0.0,0.0};
  
  // Pb-Pb from 2010 run range is default: LHC10h
  runNBins = 2500;
  runHistRange[0] = 137100.;
  runHistRange[1] = 139600.;
  
  // Pb-Pb of 2011
  if(prod.Contains("LHC11h")) {
    runNBins = 2700;
    runHistRange[0] = 167900.;
    runHistRange[1] = 170600.;
  }
  
  // Pb-Pb of 2015
  if(prod.Contains("LHC15o")) {
     runNBins = 2100;
     runHistRange[0] = 244900.;
     runHistRange[1] = 247000.;
  }
  
  // p-Pb of 2013
  if(prod.Contains("LHC13b") || prod.Contains("LHC13c")) {
    runNBins = 400;
    runHistRange[0] = 195300.;
    runHistRange[1] = 195700.;
  }
  
  // pp at 13 TeV
  if(prod.Contains("LHC16l")) {
     runNBins = 1140;
     runHistRange[0] = 258880.;
     runHistRange[1] = 260020.;
  }
  
  // p-Pb at 8.16 TeV
  if(prod.Contains("LHC16r")) {
     runNBins = 1000;
     runHistRange[0] = 265400.;
     runHistRange[1] = 266400.;
  }
  
  const Int_t kNBCBins = 3600;
  Double_t bcHistRange[2] = {-0.5,3599.5};
  
  TString classesStr(histClasses);
  TObjArray* arr=classesStr.Tokenize(";");
  
  cout << "Histogram classes included in the Histogram Manager" << endl;
  for(Int_t iclass=0; iclass<arr->GetEntries(); ++iclass) {
    TString classStr = arr->At(iclass)->GetName();
    
    if(classStr.Contains("MCTruth_")) {
       man->AddHistClass(classStr.Data());
       cout << classStr.Data() << endl;
       const Int_t kNPtBins = 18;
       Double_t ptBins[kNPtBins] = {
          0.0, 0.02, 0.04, 0.06, 0.08, 
          0.1, 0.12, 0.14, 0.16, 0.18,
          0.3, 0.5, 1.0, 3.0, 5.0,
          7.0, 10.0, 20.0
       };
       const Int_t kNMassBins = 126;
       Double_t massBins[kNMassBins];
       for(Int_t i=0; i<kNMassBins;++i) massBins[i] = 0.0+i*0.04;
       const Int_t kNCentBins = 7;
       Double_t centBins[kNCentBins] = {30.,40.,50.,60.,70.,80.,90.};
       
       man->AddHistogram(classStr.Data(), "MassMC_Pt_CentVZERO", "", kFALSE, kNMassBins-1, massBins, AliReducedVarManager::kMassMC, kNPtBins-1, ptBins, AliReducedVarManager::kPtMC, kNCentBins-1, centBins, AliReducedVarManager::kCentVZERO);
       
       man->AddHistogram(classStr.Data(), "MassMC", "MC mass", kFALSE, 200, 0., 5.0, AliReducedVarManager::kMassMC);
       man->AddHistogram(classStr.Data(), "RapidityMC", "MC rapidity", kFALSE, 48, -1.2, 1.2, AliReducedVarManager::kRapMC);
       man->AddHistogram(classStr.Data(), "PtMC", "p_{T} MC", kFALSE, 1000, 0., 10.0, AliReducedVarManager::kPtMC);
       man->AddHistogram(classStr.Data(), "PtMC_coarse", "p_{T} MC", kFALSE, 20, 0., 20.0, AliReducedVarManager::kPtMC);
       man->AddHistogram(classStr.Data(), "PhiMC", "MC #varphi", kFALSE, 100, 0., 6.3, AliReducedVarManager::kPhiMC);
       man->AddHistogram(classStr.Data(), "EtaMC", "MC #eta", kFALSE, 100, -1.5, 1.5, AliReducedVarManager::kEtaMC);
       man->AddHistogram(classStr.Data(), "PtMC_RapMC", "", kFALSE, 100, -1.2, 1.2, AliReducedVarManager::kRapMC, 100, 0., 15., AliReducedVarManager::kPtMC);
       man->AddHistogram(classStr.Data(), "CosThetaStarCS", "cos(#theta^{*})_{CS}", kFALSE, 22, -1.1, 1.1, AliReducedVarManager::kPairThetaCS);
       man->AddHistogram(classStr.Data(), "CosThetaStarCS_ptMC", "cos(#theta^{*})_{CS} vs MC pt", kFALSE, 22, -1.1, 1.1, AliReducedVarManager::kPairThetaCS, 50, 0., 1., AliReducedVarManager::kPtMC);
       man->AddHistogram(classStr.Data(), "CosThetaStarHE", "cos(#theta^{*})_{HE}", kFALSE, 22, -1.1, 1.1, AliReducedVarManager::kPairThetaHE);
       man->AddHistogram(classStr.Data(), "CosThetaStarHE_ptMC", "cos(#theta^{*})_{HE} vs MC pt", kFALSE, 22, -1.1, 1.1, AliReducedVarManager::kPairThetaHE, 100, 0., 1., AliReducedVarManager::kPtMC);
       man->AddHistogram(classStr.Data(), "CosThetaStarHE_ptMC_coarse", "cos(#theta^{*})_{HE} vs MC pt", kFALSE, 10, -1.0, 1.0, AliReducedVarManager::kPairThetaHE, 20, 0., 20., AliReducedVarManager::kPtMC);
       man->AddHistogram(classStr.Data(), "PhiStarCS", "#varphi^{*}_{CS}", kFALSE, 22, -3.3, 3.3, AliReducedVarManager::kPairPhiCS);
       man->AddHistogram(classStr.Data(), "PhiStarHE", "#varphi^{*}_{HE}", kFALSE, 22, -3.3, 3.3, AliReducedVarManager::kPairPhiHE);
       
       continue;
    }
    
    // Event wise histograms
    if(classStr.Contains("EventTag_")) {
       man->AddHistClass(classStr.Data());
       cout << classStr.Data() << endl;
       TString tagNames = "";
       tagNames += "AliAnalysisUtils 2013 selection;";
       tagNames += "AliAnalysisUtils MV pileup;";
       tagNames += "AliAnalysisUtils MV pileup, no BC check;";
       tagNames += "AliAnalysisUtils MV pileup, min wght dist 10;";
       tagNames += "AliAnalysisUtils MV pileup, min wght dist 5;";
       tagNames += "IsPileupFromSPD(3,0.6,3.,2.,5.);";
       tagNames += "IsPileupFromSPD(4,0.6,3.,2.,5.);";
       tagNames += "IsPileupFromSPD(5,0.6,3.,2.,5.);";
       tagNames += "IsPileupFromSPD(6,0.6,3.,2.,5.);";
       tagNames += "IsPileupFromSPD(3,0.8,3.,2.,5.);";
       tagNames += "IsPileupFromSPD(4,0.8,3.,2.,5.);";
       tagNames += "IsPileupFromSPD(5,0.8,3.,2.,5.);";
       tagNames += "IsPileupFromSPD(6,0.8,3.,2.,5.);";
       tagNames += "vtx distance selected;";
       man->AddHistogram(classStr.Data(), "EventTags", "Event tags", kFALSE,
                         20, -0.5, 19.5, AliReducedVarManager::kEventTag, 0, 0.0, 0.0, AliReducedVarManager::kNothing, 0, 0.0, 0.0, AliReducedVarManager::kNothing, tagNames.Data());
       man->AddHistogram(classStr.Data(), "EventTags_CentVZERO", "Event tags vs VZERO centrality", kFALSE,
                         20, -0.5, 19.5, AliReducedVarManager::kEventTag, 9, 0.0, 90.0, AliReducedVarManager::kCentVZERO, 0, 0.0, 0.0, AliReducedVarManager::kNothing, tagNames.Data());
       continue;
    }
    
    if(classStr.Contains("EventTriggers_")) {
       man->AddHistClass(classStr.Data());
       cout << classStr.Data() << endl;
       TString triggerNames = "";
       for(Int_t i=0; i<64; ++i) {triggerNames += AliReducedVarManager::fgkOfflineTriggerNames[i]; triggerNames+=";";}
       
       man->AddHistogram(classStr.Data(), "Triggers2", "", kFALSE,
                         64, -0.5, 63.5, AliReducedVarManager::kOnlineTriggerFired2, 0, 0.0, 0.0, AliReducedVarManager::kNothing, 0, 0.0, 0.0, AliReducedVarManager::kNothing, triggerNames.Data());
       man->AddHistogram(classStr.Data(), "CentVZERO_Triggers2", "", kFALSE,
                         64, -0.5, 63.5, AliReducedVarManager::kOnlineTriggerFired2, 9, 0.0, 90.0, AliReducedVarManager::kCentVZERO, 0, 0.0, 0.0, AliReducedVarManager::kNothing, triggerNames.Data());
       continue;
    }
    
    if(classStr.Contains("Event_")) {
      man->AddHistClass(classStr.Data());
      cout << classStr.Data() << endl;
      man->AddHistogram(classStr.Data(),"RunNo","Run numbers",kFALSE, runNBins, runHistRange[0], runHistRange[1], AliReducedVarManager::kRunNo);
      man->AddHistogram(classStr.Data(),"TimeFromSOR","Events vs time",kFALSE, 450, 0., 450., AliReducedVarManager::kTimeRelativeSOR);
      man->AddHistogram(classStr.Data(),"VtxX","Vtx X",kFALSE,300,-0.4,0.4,AliReducedVarManager::kVtxX);
      man->AddHistogram(classStr.Data(),"VtxX_TimeFromSOR_prof","<Vtx X> vs time from SOR",kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 300,-0.4,0.4,AliReducedVarManager::kVtxX);
      man->AddHistogram(classStr.Data(),"VtxY","Vtx Y",kFALSE,300,-0.4,0.4,AliReducedVarManager::kVtxY);
      man->AddHistogram(classStr.Data(),"VtxY_TimeFromSOR_prof","<Vtx Y> vs time from SOR",kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 300,-0.4,0.4,AliReducedVarManager::kVtxY);
      man->AddHistogram(classStr.Data(),"VtxZ","Vtx Z",kFALSE,300,-15.,15.,AliReducedVarManager::kVtxZ);
      man->AddHistogram(classStr.Data(),"VtxZ_TimeFromSOR_prof","<Vtx Z> vs time from SOR",kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 300,-15.0,15.0,AliReducedVarManager::kVtxZ);
      man->AddHistogram(classStr.Data(),"CentVZERO","Centrality(VZERO)",kFALSE, 100, 0.0, 100.0, AliReducedVarManager::kCentVZERO);
      man->AddHistogram(classStr.Data(),"CentVZERO_VtxZ","Centrality(VZERO) vs vtxZ",kFALSE, 50, 0.0, 100.0, AliReducedVarManager::kCentVZERO,
         50, -10., 10., AliReducedVarManager::kVtxZ);
      man->AddHistogram(classStr.Data(),"CentVZERO_TimeFromSOR","Centrality(VZERO) vs time from SOR",kFALSE,
                        90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 50, 0.0, 100.0, AliReducedVarManager::kCentVZERO);
      man->AddHistogram(classStr.Data(),"CentSPD","Centrality(SPD)",kFALSE, 100, 0.0, 100.0, AliReducedVarManager::kCentSPD);
      man->AddHistogram(classStr.Data(),"CentSPD_VtxZ","Centrality(SPD) vs vtxZ",kFALSE, 50, 0.0, 100.0, AliReducedVarManager::kCentSPD,
         50, -10., 10., AliReducedVarManager::kVtxZ);
      man->AddHistogram(classStr.Data(),"CentTPC","Centrality(TPC)",kFALSE, 100, 0.0, 100.0, AliReducedVarManager::kCentTPC);
      man->AddHistogram(classStr.Data(),"CentTPC_VtxZ","Centrality(TPC) vs vtxZ",kFALSE, 100, 0.0, 100.0, AliReducedVarManager::kCentTPC,
         50, -10., 10., AliReducedVarManager::kVtxZ);
      man->AddHistogram(classStr.Data(),"CentQuality","Centrality quality",kFALSE, 500, 0., 500., AliReducedVarManager::kCentQuality);
      man->AddHistogram(classStr.Data(),"IRIntClosestInt1Map","Closest out of bunch interactions Int1",kFALSE,7000,-3500.,3500.,AliReducedVarManager::kIRIntClosestIntMap);
      man->AddHistogram(classStr.Data(),"IRIntClosestInt2Map","Closest out of bunch interactions Int2",kFALSE,7000,-3500.,3500.,AliReducedVarManager::kIRIntClosestIntMap+1);
      man->AddHistogram(classStr.Data(),"NTracksTotal","Number of total tracks per event",kFALSE,500,0.,30000.,AliReducedVarManager::kNtracksTotal);
      man->AddHistogram(classStr.Data(),"NTracksTotal_BeamIntensity0_prof","Number of total tracks per event",kTRUE,100, 3.0e+12, 9.0e+12, AliReducedVarManager::kBeamIntensity0, 500,0.,20000.,AliReducedVarManager::kNtracksTotal);
      man->AddHistogram(classStr.Data(),"NTracksSelected","Number of selected tracks per event",kFALSE,300,0.,300.,AliReducedVarManager::kNtracksSelected);
      man->AddHistogram(classStr.Data(),"NTracksSelected_TimeFromSOR","Averaged number of selected tracks per event vs time from SOR",kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 100,0.,100.,AliReducedVarManager::kNtracksSelected);
      man->AddHistogram(classStr.Data(),"NTracksSelected_CentVZERO_TimeFromSOR","Averaged number of selected tracks per event per centrality vs time from SOR",kTRUE, 20, 0.0, 100.0, AliReducedVarManager::kCentVZERO, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 100,0.,100.,AliReducedVarManager::kNtracksSelected);
      man->AddHistogram(classStr.Data(),"DeltaVtxZ","Z_{global}-Z_{TPC}",kFALSE,300,-1.,1.,AliReducedVarManager::kDeltaVtxZ);
      man->AddHistogram(classStr.Data(),"DeltaVtxZ_TimeFromSOR_prof","Z_{global}-Z_{TPC} vs time from SOR",kTRUE,90, 0.0, 450.,  AliReducedVarManager::kTimeRelativeSOR, 300,-1.,1.,AliReducedVarManager::kDeltaVtxZ);
      man->AddHistogram(classStr.Data(),"NVtxTPCContributors","",kFALSE,500,0.,10000.,AliReducedVarManager::kNVtxTPCContributors);
      man->AddHistogram(classStr.Data(),"NVtxTPCContributors_BeamIntensity0","",kTRUE, 100, 3.0e+12, 9.0e+12, AliReducedVarManager::kBeamIntensity0, 500,0.,10000.,AliReducedVarManager::kNVtxTPCContributors);
      man->AddHistogram(classStr.Data(),"SPDntracklets", "SPD #tracklets in |#eta|<1.0", kFALSE, 200, 0., 5000., AliReducedVarManager::kSPDntracklets);
      man->AddHistogram(classStr.Data(),"SPDntracklets_TimeFromSOR_prof", "SPD <#tracklets> in |#eta|<1.0 vs time from SOR", kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 200, 0., 5000., AliReducedVarManager::kSPDntracklets);
      for(Int_t il=0; il<2; ++il) {
        man->AddHistogram(classStr.Data(), Form("SPDfiredChips_layer%d",il+1), Form("SPD fired chips in layer %d",il+1), 
			  kFALSE, 200, 0., 1000., AliReducedVarManager::kSPDFiredChips+il);
        man->AddHistogram(classStr.Data(), Form("SPDfiredChips_layer%d_TimeFromSOR_prof",il+1), Form("SPD <#fired chips> in layer %d vs time from SOR",il+1), 
                          kTRUE, 90, 0.0, 450., AliReducedVarManager::kTimeRelativeSOR, 200, 0., 1000., AliReducedVarManager::kSPDFiredChips+il);
      }
      for(Int_t il=0; il<6; ++il) {
        man->AddHistogram(classStr.Data(), Form("ITSclusters_layer%d",il+1), Form("ITS clusters in layer %d",il+1), 
			  kFALSE, 300, 0., 15000., AliReducedVarManager::kITSnClusters+il);
        man->AddHistogram(classStr.Data(), Form("ITSclusters_layer%d_TimeFromSOR_prof",il+1), Form("ITS <clusters> in layer %d vs time from SOR",il+1), 
                          kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 300, 0., 15000., AliReducedVarManager::kITSnClusters+il);
      }
      man->AddHistogram(classStr.Data(), "SPDnSingleClusters", "SPD single clusters", 
			kFALSE, 200, 0., 10000., AliReducedVarManager::kSPDnSingleClusters);	
      man->AddHistogram(classStr.Data(), "SPDnSingleClusters_TimeFromSOR_prof", "SPD single <#clusters> vs time from SOR", 
                        kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 200, 0., 10000., AliReducedVarManager::kSPDnSingleClusters);	
      man->AddHistogram(classStr.Data(), "SPDnSingleClusters_BeamInt0_prof", "SPD single <#clusters> vs beam intensity", 
                        kTRUE, 100, 3.0e+12, 9.0e+12., AliReducedVarManager::kBeamIntensity0, 200, 0., 10000., AliReducedVarManager::kSPDnSingleClusters);	
      man->AddHistogram(classStr.Data(),"VZEROmult", "", kFALSE, 500, 0.0, 50000., AliReducedVarManager::kVZEROTotalMult);
      man->AddHistogram(classStr.Data(),"VZEROmult_VtxContributors", "", kFALSE, 
		   200, 0.0, 40000., AliReducedVarManager::kVZEROTotalMult, 200, 0.0, 10000., AliReducedVarManager::kNVtxContributors);
      man->AddHistogram(classStr.Data(),"VZEROmult_SPDntracklets", "", kFALSE, 
                        200, 0.0, 5000., AliReducedVarManager::kSPDntracklets, 200, 0.0, 40000., AliReducedVarManager::kVZEROTotalMult);
      man->AddHistogram(classStr.Data(),"IsPhysicsSelection","Physics selection flag",kFALSE,
                        2,-0.5,1.5,AliReducedVarManager::kIsPhysicsSelection, 0,0.0,0.0,AliReducedVarManager::kNothing, 0,0.0,0.0,AliReducedVarManager::kNothing, "off;on");
      man->AddHistogram(classStr.Data(),"IsPhysicsSelection_TimeFromSOR_prof","Fraction of physics selection selected",kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 2,-0.5,1.5,AliReducedVarManager::kIsPhysicsSelection);
      man->AddHistogram(classStr.Data(),"NPosTracksAnalyzed","#positrons per event",kFALSE,100,0.,100.,AliReducedVarManager::kNtracksPosAnalyzed);
      man->AddHistogram(classStr.Data(),"NNegTracksAnalyzed","#electrons per event",kFALSE,100,0.,100.,AliReducedVarManager::kNtracksNegAnalyzed);
      man->AddHistogram(classStr.Data(),"NTotalTracksAnalyzed","#leg candidates per event",kFALSE,100,0.,100.,AliReducedVarManager::kNtracksAnalyzed);
      man->AddHistogram(classStr.Data(),"NTotalPairsAnalyzed","#dielectron candidates per event",kFALSE,100,0.,100.,AliReducedVarManager::kNpairsSelected);
      man->AddHistogram(classStr.Data(),"NPosTracksAnalyzed_TimeFromSOR_prof","<#positrons> per event vs time from SOR",kTRUE,90, 0.0, 450., AliReducedVarManager::kTimeRelativeSOR, 100,0.,100.,AliReducedVarManager::kNtracksPosAnalyzed);
      man->AddHistogram(classStr.Data(),"NNegTracksAnalyzed_TimeFromSOR_prof","<#electrons> per event vs time from SOR",kTRUE,90, 0.0, 450., AliReducedVarManager::kTimeRelativeSOR,100,0.,100.,AliReducedVarManager::kNtracksNegAnalyzed);
      man->AddHistogram(classStr.Data(),"NTotalTracksAnalyzed_TimeFromSOR_prof","<#leg candidates> per event vs time from SOR",kTRUE,90, 0.0, 450., AliReducedVarManager::kTimeRelativeSOR,100,0.,100.,AliReducedVarManager::kNtracksAnalyzed);
      man->AddHistogram(classStr.Data(),"NTotalTracksAnalyzed_BeamInt0_prof","<#leg candidates> per event vs beam intensity",kTRUE,100, 3.0e+12, 9.0e+12, AliReducedVarManager::kBeamIntensity0,100,0.,100.,AliReducedVarManager::kNtracksAnalyzed);
      man->AddHistogram(classStr.Data(),"NTotalPairsAnalyzed_TimeFromSOR_prof","<#dielectron> candidates per event vs time from SOR",kTRUE,90, 0.0, 450., AliReducedVarManager::kTimeRelativeSOR,100,0.,100.,AliReducedVarManager::kNpairsSelected);
      man->AddHistogram(classStr.Data(),"NTotalPairsAnalyzed_BeamInt0_prof","<#dielectron> candidates vs beam intensity",kTRUE,100, 3.0e+12, 9.0e+12, AliReducedVarManager::kBeamIntensity0,100,0.,100.,AliReducedVarManager::kNpairsSelected);
      
      /*for(Int_t i=0;i<36;++i)
         man->AddHistogram(classStr.Data(),Form("NTotalTracksAnalyzedInPhiBins_phiSec%d_%s_TimeFromSOR_prof", i%18, (i<18 ? "negEta" : "posEta")), Form("<#leg candidates> per event in #varphi sector %d and %s vs time from SOR", i%18, (i<18 ? "negative #eta" : "positive #eta")),kTRUE,90, 0.0, 450., AliReducedVarManager::kTimeRelativeSOR,100,0.,100.,AliReducedVarManager::kNtracksAnalyzedInPhiBins+i);
      */
      
      man->AddHistogram(classStr.Data(),"TZEROpileup", "TZERO pileup", kFALSE, 2, -0.5, 1.5, AliReducedVarManager::kTZEROpileup);
      man->AddHistogram(classStr.Data(),"TZEROpileup_TimeFromSOR_prof", "TZERO pileup vs time from SOR", kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 2, -0.5, 1.5, AliReducedVarManager::kTZEROpileup);
      man->AddHistogram(classStr.Data(),"TZEROsatellite", "TZERO satellite", kFALSE, 2, -0.5, 1.5, AliReducedVarManager::kTZEROsatellite);
      man->AddHistogram(classStr.Data(),"TZEROsatellite_TimeFromSOR_prof", "TZERO satellite vs time from SOR", kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 2, -0.5, 1.5, AliReducedVarManager::kTZEROsatellite);
      
      /*
      man->AddHistogram(classStr.Data(),"NTracksTPCoutVsITSout","",kFALSE,100,0.,15.,AliReducedVarManager::kNTracksTPCoutVsITSout);
      man->AddHistogram(classStr.Data(),"NTracksTPCoutVsITSout_TimeFromSOR_prof","",kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 100,0.,15.,AliReducedVarManager::kNTracksTPCoutVsITSout);
      man->AddHistogram(classStr.Data(),"NTracksTRDoutVsITSout","",kFALSE,100,0.,5.,AliReducedVarManager::kNTracksTRDoutVsITSout);
      man->AddHistogram(classStr.Data(),"NTracksTRDoutVsITSout_TimeFromSOR_prof","",kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 100,0.,15.,AliReducedVarManager::kNTracksTRDoutVsITSout);
      man->AddHistogram(classStr.Data(),"NTracksTOFoutVsITSout","",kFALSE,100,0.,1.,AliReducedVarManager::kNTracksTOFoutVsITSout);
      man->AddHistogram(classStr.Data(),"NTracksTOFoutVsITSout_TimeFromSOR_prof","",kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 100,0.,15.,AliReducedVarManager::kNTracksTOFoutVsITSout);
      man->AddHistogram(classStr.Data(),"NTracksTRDoutVsTPCout","",kFALSE,100,0.,1.,AliReducedVarManager::kNTracksTRDoutVsTPCout);
      man->AddHistogram(classStr.Data(),"NTracksTRDoutVsTPCout_TimeFromSOR_prof","",kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 100,0.,15.,AliReducedVarManager::kNTracksTRDoutVsTPCout);
      man->AddHistogram(classStr.Data(),"NTracksTOFoutVsTPCout","",kFALSE,100,0.,0.6,AliReducedVarManager::kNTracksTOFoutVsTPCout);
      man->AddHistogram(classStr.Data(),"NTracksTOFoutVsTPCout_TimeFromSOR_prof","",kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 100,0.,15.,AliReducedVarManager::kNTracksTOFoutVsTPCout);
      man->AddHistogram(classStr.Data(),"NTracksTOFoutVsTRDout","",kFALSE,100,0.,2.,AliReducedVarManager::kNTracksTOFoutVsTRDout);
      man->AddHistogram(classStr.Data(),"NTracksTOFoutVsTRDout_TimeFromSOR_prof","",kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 100,0.,15.,AliReducedVarManager::kNTracksTOFoutVsTRDout);
      man->AddHistogram(classStr.Data(),"NTracksITSoutVsSPDtracklets","",kFALSE,100,0.,50.,AliReducedVarManager::kNTracksITSoutVsSPDtracklets);
      man->AddHistogram(classStr.Data(),"NTracksITSoutVsSPDtracklets_TimeFromSOR_prof","",kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 100,0.,15.,AliReducedVarManager::kNTracksITSoutVsSPDtracklets);
      man->AddHistogram(classStr.Data(),"NTracksTPCoutVsSPDtracklets","",kFALSE,100,0.,50.,AliReducedVarManager::kNTracksTPCoutVsSPDtracklets);
      man->AddHistogram(classStr.Data(),"NTracksTPCoutVsSPDtracklets_TimeFromSOR_prof","",kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 100,0.,15.,AliReducedVarManager::kNTracksTPCoutVsSPDtracklets);
      man->AddHistogram(classStr.Data(),"NTracksTRDoutVsSPDtracklets","",kFALSE,100,0.,30.,AliReducedVarManager::kNTracksTRDoutVsSPDtracklets);
      man->AddHistogram(classStr.Data(),"NTracksTRDoutVsSPDtracklets_TimeFromSOR_prof","",kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 100,0.,15.,AliReducedVarManager::kNTracksTRDoutVsSPDtracklets);
      man->AddHistogram(classStr.Data(),"NTracksTOFoutVsSPDtracklets","",kFALSE,100,0.,10.,AliReducedVarManager::kNTracksTOFoutVsSPDtracklets);
      man->AddHistogram(classStr.Data(),"NTracksTOFoutVsSPDtracklets_TimeFromSOR_prof","",kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 100,0.,15.,AliReducedVarManager::kNTracksTOFoutVsSPDtracklets);*/
      man->AddHistogram(classStr.Data(),"NTracksITSin","",kFALSE,300,0.,15000.,AliReducedVarManager::kNTracksPerTrackingStatus+AliReducedVarManager::kITSin);
      man->AddHistogram(classStr.Data(),"NTracksITSout","",kFALSE,300,0.,15000.,AliReducedVarManager::kNTracksPerTrackingStatus+AliReducedVarManager::kITSout);
      man->AddHistogram(classStr.Data(),"NTracksTPCin","",kFALSE,300,0.,30000.,AliReducedVarManager::kNTracksPerTrackingStatus+AliReducedVarManager::kTPCin);
      man->AddHistogram(classStr.Data(),"NTracksTPCout","",kFALSE,1800,0.,30000.,AliReducedVarManager::kNTracksPerTrackingStatus+AliReducedVarManager::kTPCout);
      man->AddHistogram(classStr.Data(),"NTracksTRDout","",kFALSE,300,0.,15000.,AliReducedVarManager::kNTracksPerTrackingStatus+AliReducedVarManager::kTRDout);
      man->AddHistogram(classStr.Data(),"NTracksTRDin","",kFALSE,300,0.,15000.,AliReducedVarManager::kNTracksPerTrackingStatus+AliReducedVarManager::kTRDin);
      man->AddHistogram(classStr.Data(),"NTracksTPCout_VZEROmult","",kFALSE,900,0.,30000.,AliReducedVarManager::kNTracksPerTrackingStatus+AliReducedVarManager::kTPCout, 300, 0., 40000., AliReducedVarManager::kVZEROTotalMult);
      man->AddHistogram(classStr.Data(),"NTracksTPCrefit_VZEROmult","",kFALSE,900,0.,30000.,AliReducedVarManager::kNTracksPerTrackingStatus+AliReducedVarManager::kTPCrefit, 300, 0., 40000., AliReducedVarManager::kVZEROTotalMult);
      man->AddHistogram(classStr.Data(),"NTracksTPCoutFromPileup_VZEROmult","",kFALSE,1200,-10000.,30000.,AliReducedVarManager::kNTracksTPCoutFromPileup, 400, 0., 40000., AliReducedVarManager::kVZEROTotalMult);
      man->AddHistogram(classStr.Data(),"NTracksTPCoutFromPileup_VZEROmult_nAnalyzedPairs_prof","",kTRUE,40,-10000.,30000.,AliReducedVarManager::kNTracksTPCoutFromPileup, 40, 0., 40000., AliReducedVarManager::kVZEROTotalMult, 10, 0., 100., AliReducedVarManager::kNpairsSelected);
      man->AddHistogram(classStr.Data(),"NTracksTPCoutFromPileup_TimeFromSOR_nAnalyzedPairs_prof","",kTRUE, 45, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 40,-10000.,30000.,AliReducedVarManager::kNTracksTPCoutFromPileup, 10, 0., 100., AliReducedVarManager::kNpairsSelected);
      man->AddHistogram(classStr.Data(),"NTracksTPCout_VZEROmult_nAnalyzedPairs_prof","",kTRUE,60,0.,30000.,AliReducedVarManager::kNTracksPerTrackingStatus+AliReducedVarManager::kTPCout, 40, 0., 40000., AliReducedVarManager::kVZEROTotalMult, 10, 0., 100., AliReducedVarManager::kNpairsSelected);
      man->AddHistogram(classStr.Data(),"NTracksTPCoutFromPileup","",kFALSE,2000,-10000.,30000.,AliReducedVarManager::kNTracksTPCoutFromPileup);
      man->AddHistogram(classStr.Data(),"NTracksTRDout_VZEROmult","",kFALSE,300,0.,15000.,AliReducedVarManager::kNTracksPerTrackingStatus+AliReducedVarManager::kTRDout, 300, 0., 40000., AliReducedVarManager::kVZEROTotalMult);
      man->AddHistogram(classStr.Data(),"NTracksTRDrefit_VZEROmult","",kFALSE,300,0.,15000.,AliReducedVarManager::kNTracksPerTrackingStatus+AliReducedVarManager::kTRDrefit, 300, 0., 40000., AliReducedVarManager::kVZEROTotalMult);
      man->AddHistogram(classStr.Data(),"NTracksITSout_VZEROmult","",kFALSE,300,0.,30000.,AliReducedVarManager::kNTracksPerTrackingStatus+AliReducedVarManager::kITSout, 300, 0., 40000., AliReducedVarManager::kVZEROTotalMult);
      man->AddHistogram(classStr.Data(),"NTracksITSrefit_VZEROmult","",kFALSE,300,0.,30000.,AliReducedVarManager::kNTracksPerTrackingStatus+AliReducedVarManager::kITSrefit, 300, 0., 40000., AliReducedVarManager::kVZEROTotalMult);
      man->AddHistogram(classStr.Data(),"NTracksTOFout_VZEROmult","",kFALSE,300,0.,3000.,AliReducedVarManager::kNTracksPerTrackingStatus+AliReducedVarManager::kTOFout, 300, 0., 40000., AliReducedVarManager::kVZEROTotalMult);
      man->AddHistogram(classStr.Data(),"NTracksTOFrefit_VZEROmult","",kFALSE,300,0.,3000.,AliReducedVarManager::kNTracksPerTrackingStatus+AliReducedVarManager::kTOFrefit, 300, 0., 40000., AliReducedVarManager::kVZEROTotalMult);
      man->AddHistogram(classStr.Data(),"NTracksTotal_VZEROmult","",kFALSE,300,0.,30000.,AliReducedVarManager::kNtracksTotal, 300, 0., 40000., AliReducedVarManager::kVZEROTotalMult);
      man->AddHistogram(classStr.Data(),"NTracksSelected_VZEROmult","",kFALSE,300,0.,300.,AliReducedVarManager::kNtracksSelected, 300, 0., 40000., AliReducedVarManager::kVZEROTotalMult);
      man->AddHistogram(classStr.Data(),"NTracksAnalyzed_VZEROmult","",kFALSE,100,0.,100.,AliReducedVarManager::kNtracksAnalyzed, 300, 0., 40000., AliReducedVarManager::kVZEROTotalMult);
      man->AddHistogram(classStr.Data(),"NTracksTPCoutVsVZEROmult","",kFALSE,1000,0.,10.,AliReducedVarManager::kNTracksTPCoutVsVZEROTotalMult);
      man->AddHistogram(classStr.Data(),"NTracksTPCoutVsVZEROmult_CentVZERO","",kFALSE, 20, 0., 100., AliReducedVarManager::kCentVZERO, 1000,0.,10.,AliReducedVarManager::kNTracksTPCoutVsVZEROTotalMult);
      man->AddHistogram(classStr.Data(),"NTracksAnalyzed_RunNo_prof","",kTRUE, runNBins, runHistRange[0], runHistRange[1], AliReducedVarManager::kRunNo, 10, 0., 10., AliReducedVarManager::kNtracksAnalyzed);
      man->AddHistogram(classStr.Data(),"NPairsAnalyzed_RunNo_prof","",kTRUE, runNBins, runHistRange[0], runHistRange[1], AliReducedVarManager::kRunNo, 10, 0., 10., AliReducedVarManager::kNpairsSelected);
      man->AddHistogram(classStr.Data(),"EventAverageTPCchi2","Average TPC chi2 per track",kFALSE,100,0.,5.,AliReducedVarManager::kEvAverageTPCchi2);
      man->AddHistogram(classStr.Data(),"EventAverageTPCchi2_run_prof","Average TPC chi2 per track vs run",kTRUE,runNBins, runHistRange[0], runHistRange[1], AliReducedVarManager::kRunNo, 100,0.,5.,AliReducedVarManager::kEvAverageTPCchi2);
      man->AddHistogram(classStr.Data(),"EventAverageTPCchi2_NTPCout","Average TPC chi2 per track vs n TPCout",kFALSE,100,0.,5.,AliReducedVarManager::kEvAverageTPCchi2, 50, 0.0, 15000., AliReducedVarManager::kNTracksPerTrackingStatus+AliReducedVarManager::kTPCout);
      man->AddHistogram(classStr.Data(),"EventAverageTPCchi2_NTPCout_prof","Average TPC chi2 per track vs n TPCout",kTRUE,50, 0.0, 15000., AliReducedVarManager::kNTracksPerTrackingStatus+AliReducedVarManager::kTPCout, 100,0.,5.,AliReducedVarManager::kEvAverageTPCchi2);
      
      
      Bool_t isCalibrated = kTRUE;
      /*for(Int_t i=0; i<64; ++i) {
         man->AddHistogram(classStr.Data(), Form("VZEROmult_ch%d",i), Form("Multiplicity in VZERO channel %d",i), kFALSE, 
                      200, 0.0, (isCalibrated ? 10.0 : 1000.0), AliReducedVarManager::kVZEROChannelMult+i);
         man->AddHistogram(classStr.Data(), Form("VZEROmult_ch%d_zoomLowMult",i), Form("Multiplicity in VZERO channel %d",i), kFALSE, 
                           100, 0.0, (isCalibrated ? 0.1 : 10.0), AliReducedVarManager::kVZEROChannelMult+i);
         man->AddHistogram(classStr.Data(),Form("VZEROmult_ch%d_Run_prof",i),Form("Multiplicity in VZERO channel %d vs run",i), kTRUE,
                      runNBins, runHistRange[0], runHistRange[1], AliReducedVarManager::kRunNo, 500, 0.0, 1000.0, AliReducedVarManager::kVZEROChannelMult+i);
         man->AddHistogram(classStr.Data(),Form("VZEROmult_ch%d_VtxZ_prof",i),Form("Multiplicity in VZERO channel %d vs event vertex Z",i), kTRUE,
                      20, -10., 10., AliReducedVarManager::kVtxZ, 500, 0.0, 1000.0, AliReducedVarManager::kVZEROChannelMult+i);
         man->AddHistogram(classStr.Data(),Form("VZEROmult_ch%d_VtxCent_prof",i),Form("Multiplicity in VZERO channel %d vs (vtxZ,centrality SPD);",i),kTRUE,
                      20, -10.,10., AliReducedVarManager::kVtxZ, 18, 0.0, 90.0, AliReducedVarManager::kCentSPD, 100, 0., 1000., AliReducedVarManager::kVZEROChannelMult+i);
         man->AddHistogram(classStr.Data(),Form("VZEROmultOccupancy_ch%d_VtxCent",i),Form("Multiplicity in VZERO channel %d vs (vtxZ,centrality SPD);",i),kFALSE,
                      20, -10., 10., AliReducedVarManager::kVtxZ, 18, 0., 90., AliReducedVarManager::kCentSPD, 100, 0., 1000., AliReducedVarManager::kVZEROChannelMult+i);
      }  // end loop over channels
      */
      man->AddHistogram(classStr.Data(), "VZEROA_NEmptyChannels_VtxCent_prof", "No. VZERO-A empty channels per event vs. centrality SPD and vertex Z", kTRUE,
                   24, -12.0, +12.0, AliReducedVarManager::kVtxZ, 20, 0.0, 100., AliReducedVarManager::kCentSPD, 32, 0.0, 32.0, AliReducedVarManager::kVZEROAemptyChannels);
      man->AddHistogram(classStr.Data(), "VZEROC_NEmptyChannels_VtxCent_prof", "No. VZERO-C empty channels per event vs. centrality SPD and vertex Z", kTRUE,
                   24, -12.0, +12.0, AliReducedVarManager::kVtxZ, 20, 0.0, 100., AliReducedVarManager::kCentSPD, 32, 0.0, 32.0, AliReducedVarManager::kVZEROCemptyChannels);
      
      const Char_t sname[2] = {'A','C'};
      for(Int_t ih=1; ih<2; ++ih) {
         man->AddHistogram(classStr.Data(), Form("QvecX_sideA_vs_sideC_h%d", ih+1), Form("Q_{x}, side A vs side C, harmonic %d",ih+1), kFALSE, 
                      100, (isCalibrated ? -20.0 : -1500.0), (isCalibrated ? 20.0 : 1500.0), AliReducedVarManager::kVZEROQvecX+0*6+ih,
                      100, (isCalibrated ? -20.0 : -1500.0), (isCalibrated ? 20.0 : 1500.0), AliReducedVarManager::kVZEROQvecX+1*6+ih);
         man->AddHistogram(classStr.Data(), Form("QvecY_sideA_vs_sideC_h%d", ih+1), Form("Q_{y}, side A vs side C, harmonic %d",ih+1), kFALSE, 
                           100, (isCalibrated ? -20.0 : -1500.0), (isCalibrated ? 20.0 : 1500.0), AliReducedVarManager::kVZEROQvecY+0*6+ih,
                           100, (isCalibrated ? -20.0 : -1500.0), (isCalibrated ? 20.0 : 1500.0), AliReducedVarManager::kVZEROQvecY+1*6+ih);
         man->AddHistogram(classStr.Data(), Form("RP_sideA_vs_sideC_h%d",ih+1), Form("VZERO RP-Aside vs RP-Cside, harmonic %d",ih+1), kFALSE, 
                           100, -4.0/Double_t(ih+1), 4.0/Double_t(ih+1), AliReducedVarManager::kVZERORP+0*6+ih,
                           100, -4.0/Double_t(ih+1), 4.0/Double_t(ih+1), AliReducedVarManager::kVZERORP+1*6+ih);
         man->AddHistogram(classStr.Data(), Form("VZERORPres_h%d_Cent",ih+1), Form("VZERO RP resolution vs centrality, harmonic %d", ih+1),
                           kTRUE, 100, 0.0, 100.0, AliReducedVarManager::kCentSPD, 100, -1., +1., AliReducedVarManager::kVZERORPres+ih);
         man->AddHistogram(classStr.Data(), Form("VZERORPres_h%d_VtxZ",ih+1), Form("VZERO RP resolution vs vertex Z, harmonic %d", ih+1),
                           kTRUE, 48, -12.0, 12.0, AliReducedVarManager::kVtxZ, 100, -1., +1., AliReducedVarManager::kVZERORPres+ih);
         man->AddHistogram(classStr.Data(), Form("VZERORPres_h%d_CentVtxZ",ih+1), Form("VZERO RP resolution vs (centrality,vtx. Z), harmonic %d", ih+1),
                           kTRUE, 100, 0.0, 100.0, AliReducedVarManager::kCentSPD, 48, -12.0, 12.0, AliReducedVarManager::kVtxZ, 100, -1., +1., AliReducedVarManager::kVZERORPres+ih);
         /*
         man->AddHistogram(classStr.Data(), Form("VZERO_XaXc_h%d_centSPD", ih+1), 
                      Form("VZERO Q_{x}^{A}#timesQ_{x}^{C}, harmonic %d, vs centSPD", ih+1), kTRUE,
                           20, 0.0, 100.0, AliReducedVarManager::kCentSPD, 500, (isCalibrated ? -20.0 : -2000.0), (isCalibrated ? 20.0 : 2000.0), AliReducedVarManager::kVZEROXaXc+ih);
         man->AddHistogram(classStr.Data(), Form("VZERO_XaYa_h%d_centSPD", ih+1), 
                      Form("VZERO Q_{x}^{A}#timesQ_{y}^{A}, harmonic %d, vs centSPD", ih+1), kTRUE,
                           20, 0.0, 100.0, AliReducedVarManager::kCentSPD, 500, (isCalibrated ? -20.0 : -2000.0), (isCalibrated ? 20.0 : 2000.0), AliReducedVarManager::kVZEROXaYa+ih);
         man->AddHistogram(classStr.Data(), Form("VZERO_XaYc_h%d_centSPD", ih+1), 
                      Form("VZERO Q_{x}^{A}#timesQ_{y}^{C}, harmonic %d, vs centSPD", ih+1), kTRUE,
                           20, 0.0, 100.0, AliReducedVarManager::kCentSPD, 500, (isCalibrated ? -20.0 : -2000.0), (isCalibrated ? 20.0 : 2000.0), AliReducedVarManager::kVZEROXaYc+ih);
         man->AddHistogram(classStr.Data(), Form("VZERO_YaXc_h%d_centSPD", ih+1), 
                      Form("VZERO Q_{y}^{A}#timesQ_{x}^{C}, harmonic %d, vs centSPD", ih+1), kTRUE,
                           20, 0.0, 100.0, AliReducedVarManager::kCentSPD, 500, (isCalibrated ? -20.0 : -2000.0), (isCalibrated ? 20.0 : 2000.0), AliReducedVarManager::kVZEROYaXc+ih);
         man->AddHistogram(classStr.Data(), Form("VZERO_YaYc_h%d_centSPD", ih+1), 
                      Form("VZERO Q_{y}^{A}#timesQ_{y}^{C}, harmonic %d, vs centSPD", ih+1), kTRUE,
                           20, 0.0, 100.0, AliReducedVarManager::kCentSPD, 500, (isCalibrated ? -20.0 : -2000.0), (isCalibrated ? 20.0 : 2000.0), AliReducedVarManager::kVZEROYaYc+ih);
         man->AddHistogram(classStr.Data(), Form("VZERO_XcYc_h%d_centSPD", ih+1), 
                      Form("VZERO Q_{x}^{C}#timesQ_{y}^{C}, harmonic %d, vs centSPD", ih+1), kTRUE,
                           20, 0.0, 100.0, AliReducedVarManager::kCentSPD, 500, (isCalibrated ? -20.0 : -2000.0), (isCalibrated ? 20.0 : 2000.0), AliReducedVarManager::kVZEROXcYc+ih);
         
         man->AddHistogram(classStr.Data(), Form("VZERO_XaXc_h%d_VtxZ", ih+1), 
                      Form("VZERO Q_{x}^{A}#timesQ_{x}^{C}, harmonic %d, vs VtxZ", ih+1), kTRUE,
                           24, -12.0, 12.0, AliReducedVarManager::kVtxZ, 500, (isCalibrated ? -20.0 : -2000.0), (isCalibrated ? 20.0 : 2000.0), AliReducedVarManager::kVZEROXaXc+ih);
         man->AddHistogram(classStr.Data(), Form("VZERO_XaYa_h%d_VtxZ", ih+1), 
                      Form("VZERO Q_{x}^{A}#timesQ_{y}^{A}, harmonic %d, vs VtxZ", ih+1), kTRUE,
                           24, -12.0, 12.0, AliReducedVarManager::kVtxZ, 500, (isCalibrated ? -20.0 : -2000.0), (isCalibrated ? 20.0 : 2000.0), AliReducedVarManager::kVZEROXaYa+ih);
         man->AddHistogram(classStr.Data(), Form("VZERO_XaYc_h%d_VtxZ", ih+1), 
                      Form("VZERO Q_{x}^{A}#timesQ_{y}^{C}, harmonic %d, vs VtxZ", ih+1), kTRUE,
                           24, -12.0, 12.0, AliReducedVarManager::kVtxZ, 500, (isCalibrated ? -20.0 : -2000.0), (isCalibrated ? 20.0 : 2000.0), AliReducedVarManager::kVZEROXaYc+ih);
         man->AddHistogram(classStr.Data(), Form("VZERO_YaXc_h%d_VtxZ", ih+1), 
                      Form("VZERO Q_{y}^{A}#timesQ_{x}^{C}, harmonic %d, vs VtxZ", ih+1), kTRUE,
                           24, -12.0, 12.0, AliReducedVarManager::kVtxZ, 500, (isCalibrated ? -20.0 : -2000.0), (isCalibrated ? 20.0 : 2000.0), AliReducedVarManager::kVZEROYaXc+ih);
         man->AddHistogram(classStr.Data(), Form("VZERO_YaYc_h%d_VtxZ", ih+1), 
                      Form("VZERO Q_{y}^{A}#timesQ_{y}^{C}, harmonic %d, vs VtxZ", ih+1), kTRUE,
                           24, -12.0, 12.0, AliReducedVarManager::kVtxZ, 500, (isCalibrated ? -20.0 : -2000.0), (isCalibrated ? 20.0 : 2000.0), AliReducedVarManager::kVZEROYaYc+ih);
         man->AddHistogram(classStr.Data(), Form("VZERO_XcYc_h%d_VtxZ", ih+1), 
                      Form("VZERO Q_{x}^{C}#timesQ_{y}^{C}, harmonic %d, vs VtxZ", ih+1), kTRUE,
                           24, -12.0, 12.0, AliReducedVarManager::kVtxZ, 500, (isCalibrated ? -20.0 : -2000.0), (isCalibrated ? 20.0 : 2000.0), AliReducedVarManager::kVZEROXcYc+ih);
         */
         /*for(Int_t iS=0; iS<2; ++iS) {
            man->AddHistogram(classStr.Data(), Form("QvecX_side%c_h%d_CentSPD",sname[iS],ih+1), 
                         Form("Q_{x}, side %c, harmonic %d vs CentSPD",sname[iS],ih+1), kFALSE, 
                              20, 0.0, 100.0, AliReducedVarManager::kCentSPD, 100, (isCalibrated ? -10.0 : -1500.0), (isCalibrated ? 10.0 : 1500.0), AliReducedVarManager::kVZEROQvecX+iS*6+ih);
            man->AddHistogram(classStr.Data(), Form("QvecY_side%c_h%d_CentSPD",sname[iS],ih+1), 
                         Form("Q_{y}, side %c, harmonic %d vs CentSPD",sname[iS],ih+1), kFALSE, 
                              20, 0.0, 100.0, AliReducedVarManager::kCentSPD, 100, (isCalibrated ? -10.0 : -1500.0), (isCalibrated ? 10.0 : 1500.0), AliReducedVarManager::kVZEROQvecY+iS*6+ih);
            man->AddHistogram(classStr.Data(), Form("QvecX_side%c_h%d_Run_prof", sname[iS], ih+1), 
                         Form("<Q_{x}>, VZERO side %c, harmonic %d, vs run", sname[iS], ih+1), kTRUE,
                              runNBins, runHistRange[0], runHistRange[1], AliReducedVarManager::kRunNo, 100, -100., 100., AliReducedVarManager::kVZEROQvecX+iS*6+ih);
            man->AddHistogram(classStr.Data(), Form("QvecY_side%c_h%d_Run_prof", sname[iS], ih+1), 
                         Form("<Q_{y}>, VZERO side %c, harmonic %d, vs run", sname[iS], ih+1), kTRUE,
                              runNBins, runHistRange[0], runHistRange[1], AliReducedVarManager::kRunNo, 100, -100., 100., AliReducedVarManager::kVZEROQvecY+iS*6+ih);
            man->AddHistogram(classStr.Data(), Form("QvecX_side%c_h%d_CentSPDVtxZ_prof",sname[iS],ih+1), 
                         Form("Q_{x}, side %c, harmonic %d, vs centSPD and vtxZ",sname[iS],ih+1), kTRUE, 
                              18, 0.,90., AliReducedVarManager::kCentSPD, 20, -10.,10., AliReducedVarManager::kVtxZ, 10, 0.,100., AliReducedVarManager::kVZEROQvecX+iS*6+ih);
            man->AddHistogram(classStr.Data(), Form("QvecY_side%c_h%d_CentSPDVtxZ_prof",sname[iS],ih+1), 
                         Form("Q_{y}, side %c, harmonic %d, vs centSPD and vtxZ",sname[iS],ih+1), kTRUE, 
                              18, 0.,90., AliReducedVarManager::kCentSPD, 20, -10.,10., AliReducedVarManager::kVtxZ, 10, 0.,100., AliReducedVarManager::kVZEROQvecY+iS*6+ih);
            man->AddHistogram(classStr.Data(), Form("RP_side%c_h%d_CentSPD",sname[iS],ih+1), 
                         Form("VZERO reaction plane, side %c, harmonic %d, vs centrality SPD",sname[iS],ih+1), kFALSE, 
                              100, -4.0/Double_t(ih+1), 4.0/Double_t(ih+1), AliReducedVarManager::kVZERORP+iS*6+ih, 20, 0.0, 100.0, AliReducedVarManager::kCentSPD);
            man->AddHistogram(classStr.Data(), Form("RP_side%c_h%d_VtxZ",sname[iS],ih+1), 
                         Form("VZERO reaction plane, side %c, harmonic %d, vs vtxZ",sname[iS],ih+1), kFALSE, 
                              100, -4.0/Double_t(ih+1), 4.0/Double_t(ih+1), AliReducedVarManager::kVZERORP+iS*6+ih, 24, -12.0, +12.0, AliReducedVarManager::kVtxZ);
            
         }   // end loop over VZERO sides
         */
      }   // end loop over harmonics
      
      continue;
    }  // end if className contains "Event"    
    
    // Track histograms
    if(classStr.Contains("TrackITSclusterMap_")) {
       man->AddHistClass(classStr.Data());
       cout << classStr.Data() << endl;
       man->AddHistogram(classStr.Data(), "ITSlayerHit", "Hits in the ITS layers", kFALSE,
                         6, 0.5, 6.5, AliReducedVarManager::kITSlayerHit);
       man->AddHistogram(classStr.Data(), "ITSlayerHit_Phi", "Hits in the ITS layers vs #varphi", kFALSE,
                         180, 0.0, 6.29, AliReducedVarManager::kPhi, 6, 0.5, 6.5, AliReducedVarManager::kITSlayerHit);
       man->AddHistogram(classStr.Data(), "ITSlayerHit_Eta", "Hits in the ITS layers vs #eta", kFALSE,
                         100, -1.0, 1.0, AliReducedVarManager::kEta, 6, 0.5, 6.5, AliReducedVarManager::kITSlayerHit);
       //man->AddHistogram(classStr.Data(), "ITSlayerHit_DCAxy", "Hits in the ITS layers vs DCA_{xy}", kFALSE,
       //                  1000, -0.5, 0.5, AliReducedVarManager::kDcaXY, 6, 0.5, 6.5, AliReducedVarManager::kITSlayerHit);
       //man->AddHistogram(classStr.Data(), "ITSlayerHit_DCAz", "Hits in the ITS layers vs DCA_{z}", kFALSE,
       //                  1800, -1.0, 1.0, AliReducedVarManager::kDcaZ, 6, 0.5, 6.5, AliReducedVarManager::kITSlayerHit);
       continue;
    }  // end of ITSclusterMap histogram definitions
    
    if(classStr.Contains("TrackTPCclusterMap_")) {
       man->AddHistClass(classStr.Data());
       cout << classStr.Data() << endl;
       man->AddHistogram(classStr.Data(), "TPCclusterMap", "TPC cluster map", kFALSE,
                         8, -0.5, 7.5, AliReducedVarManager::kTPCclusBitFired);
       man->AddHistogram(classStr.Data(), "TPCclusterMap_Phi", "TPC cluster map vs #varphi", kFALSE,
                         180, 0.0, 6.29, AliReducedVarManager::kPhi, 8, -0.5, 7.5, AliReducedVarManager::kTPCclusBitFired);
       man->AddHistogram(classStr.Data(), "TPCclusterMap_Eta", "TPC cluster map vs #eta", kFALSE,
                         100, -1.0, 1.0, AliReducedVarManager::kEta, 8, -0.5, 7.5, AliReducedVarManager::kTPCclusBitFired);
       man->AddHistogram(classStr.Data(), "TPCclusterMap_Pt", "TPC cluster map vs p_{T}", kFALSE,
                         100, 0.0, 10.0, AliReducedVarManager::kPt, 8, -0.5, 7.5, AliReducedVarManager::kTPCclusBitFired);
       continue;
    }  // end of TPCclusterMap histogram definitions
    
    TString trkStatusNames = "";
    for(Int_t iflag=0; iflag<AliReducedVarManager::kNTrackingStatus; ++iflag) {
       trkStatusNames += AliReducedVarManager::fgkTrackingStatusNames[iflag];
       trkStatusNames += ";";
    }
    if(classStr.Contains("TrackStatusFlags_")) {
       man->AddHistClass(classStr.Data());
       cout << classStr.Data() << endl;
       man->AddHistogram(classStr.Data(), "TrackingFlags", "Tracking flags;;", kFALSE,
                         AliReducedVarManager::kNTrackingStatus, -0.5, AliReducedVarManager::kNTrackingStatus-0.5, AliReducedVarManager::kTrackingFlag, 
                         0, 0.0, 0.0, AliReducedVarManager::kNothing, 0, 0.0, 0.0, AliReducedVarManager::kNothing, trkStatusNames.Data());
       /*man->AddHistogram(classStr.Data(),"TPCnsigElectronCorrected_Pin_TrackingFlags","Corrected TPC N_{#sigma} electron vs. inner param P vs tracking flags;;",kFALSE,
                         50,0.0,5.0,AliReducedVarManager::kPin,50,-5.0,5.0,AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron, AliReducedVarManager::kNTrackingStatus, -0.5, AliReducedVarManager::kNTrackingStatus-0.5, AliReducedVarManager::kTrackingFlag, "", "", trkStatusNames.Data());
       man->AddHistogram(classStr.Data(),"TPCncls_Eta_NTracksTPCoutFromPileup_TrackingFlags_prof","",kTRUE, 20,-1.0,1.0,AliReducedVarManager::kEta,
                         43,-1500., 20000., AliReducedVarManager::kNTracksTPCoutFromPileup, AliReducedVarManager::kNTrackingStatus, -0.5, AliReducedVarManager::kNTrackingStatus-0.5, AliReducedVarManager::kTrackingFlag, "", "", trkStatusNames.Data(),  AliReducedVarManager::kTPCncls);*/
       continue;
    }
    
    if(classStr.Contains("Track_")) {
      man->AddHistClass(classStr.Data());
      cout << classStr.Data() << endl;
      man->AddHistogram(classStr.Data(), "P_Pin", "P vs Pin", kFALSE, 200, 0.0, 10.0, AliReducedVarManager::kP, 200, 0.0, 10.0, AliReducedVarManager::kPin);
      man->AddHistogram(classStr.Data(), "Pt", "p_{T} distribution", kFALSE, 1000, 0.0, 50.0, AliReducedVarManager::kPt);
      man->AddHistogram(classStr.Data(), "Pt_TimeFromSOR", "<p_{T}> vs time from SOR", kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 1000, 0.0, 50.0, AliReducedVarManager::kPt);
      man->AddHistogram(classStr.Data(), "Eta", "", kFALSE, 1000, -1.5, 1.5, AliReducedVarManager::kEta);
      man->AddHistogram(classStr.Data(), "Eta_Phi", "", kFALSE, 100, -1.0, 1.0, AliReducedVarManager::kEta, 100, 0., 6.3, AliReducedVarManager::kPhi);
      man->AddHistogram(classStr.Data(), "Eta_TimeFromSOR", "<#eta> vs time from SOR", kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 1000, -1.5, 1.5, AliReducedVarManager::kEta);
      man->AddHistogram(classStr.Data(), "Phi", "", kFALSE, 1000, 0.0, 6.3, AliReducedVarManager::kPhi);
      man->AddHistogram(classStr.Data(), "Phi_TimeFromSOR", "<#varphi> vs time from SOR", kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 1000, 0.0, 6.3, AliReducedVarManager::kPhi);
      man->AddHistogram(classStr.Data(), "DCAxy_Pt", "DCAxy", kFALSE, 100, -2.0, 2.0, AliReducedVarManager::kDcaXY, 50, 0.0, 5.0, AliReducedVarManager::kPt);
      man->AddHistogram(classStr.Data(), "DCAxy_TPCchi2", "DCAxy vs TPC chi2", kFALSE, 100, -2.0, 2.0, AliReducedVarManager::kDcaXY, 50, 0.0, 5.0, AliReducedVarManager::kTPCchi2);
      man->AddHistogram(classStr.Data(), "DCAz_TPCchi2", "DCAz vs TPC chi2", kFALSE, 100, -5.0, 5.0, AliReducedVarManager::kDcaZ, 50, 0.0, 5.0, AliReducedVarManager::kTPCchi2);
      man->AddHistogram(classStr.Data(), "DCAxy_ITSchi2", "DCAxy vs ITS chi2", kFALSE, 100, -2.0, 2.0, AliReducedVarManager::kDcaXY, 100, 0.0, 25.0, AliReducedVarManager::kITSchi2);
      man->AddHistogram(classStr.Data(), "DCAz_ITSchi2", "DCAz vs ITS chi2", kFALSE, 100, -5.0, 5.0, AliReducedVarManager::kDcaZ, 100, 0.0, 25.0, AliReducedVarManager::kITSchi2);
      man->AddHistogram(classStr.Data(), "DCAxy_goldenChi2", "DCAxy vs golden chi2", kFALSE, 100, -2.0, 2.0, AliReducedVarManager::kDcaXY, 100, 0.0, 100.0, AliReducedVarManager::kChi2TPCConstrainedVsGlobal);
      man->AddHistogram(classStr.Data(), "DCAz_goldenChi2", "DCAz vs golden chi2", kFALSE, 100, -5.0, 5.0, AliReducedVarManager::kDcaZ, 100, 0.0, 100.0, AliReducedVarManager::kChi2TPCConstrainedVsGlobal);
      man->AddHistogram(classStr.Data(), "DCAxy", "DCAxy", kFALSE, 200, -5.0, 5.0, AliReducedVarManager::kDcaXY);
      man->AddHistogram(classStr.Data(), "DCAxy_TimeFromSOR", "<DCAxy> vs time from SOR", kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 200, -5.0, 5.0, AliReducedVarManager::kDcaXY);
      man->AddHistogram(classStr.Data(), "DCAz", "DCAz", kFALSE, 200, -5.0, 5.0, AliReducedVarManager::kDcaZ);
      man->AddHistogram(classStr.Data(), "DCAz_Pt", "DCAz", kFALSE, 100, -3.0, 3.0, AliReducedVarManager::kDcaZ, 50, 0.0, 5.0, AliReducedVarManager::kPt);
      man->AddHistogram(classStr.Data(), "DCAz_TimeFromSOR", "<DCAz> vs time from SOR", kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 200, -5.0, 5.0, AliReducedVarManager::kDcaZ);
      man->AddHistogram(classStr.Data(),"DCAxy_Eta_MultVZERO_NTracksTPCoutFromPileup_prof","",kTRUE, 10, 0., 40000., AliReducedVarManager::kVZEROTotalMult, 22,-2000., 20000., AliReducedVarManager::kNTracksTPCoutFromPileup, 20,-1.0,1.0, AliReducedVarManager::kEta, "","","", AliReducedVarManager::kDcaXY);
      man->AddHistogram(classStr.Data(),"DCAz_Eta_MultVZERO_NTracksTPCoutFromPileup_prof","",kTRUE, 10, 0., 40000., AliReducedVarManager::kVZEROTotalMult, 22,-2000., 20000., AliReducedVarManager::kNTracksTPCoutFromPileup, 20,-1.0,1.0, AliReducedVarManager::kEta, "","","", AliReducedVarManager::kDcaZ);
      
        man->AddHistogram(classStr.Data(),"ITSncls", "ITS nclusters", kFALSE, 7,-0.5,6.5,AliReducedVarManager::kITSncls);
        man->AddHistogram(classStr.Data(),"ITSnclsShared", "ITS nclusters shared", kFALSE, 7,-0.5,6.5,AliReducedVarManager::kITSnclsShared);
        man->AddHistogram(classStr.Data(),"ITSncls_TimeFromSOR", "<ITS nclusters> vs time from SOR", kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 7,-0.5,6.5,AliReducedVarManager::kITSncls);
        man->AddHistogram(classStr.Data(),"ITSnclsShared_TimeFromSOR", "<ITS nclusters shared> vs time from SOR", kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 7,-0.5,6.5,AliReducedVarManager::kITSnclsShared);
        man->AddHistogram(classStr.Data(),"Eta_Phi_ITSncls_prof","ITS <nclusters> vs (#eta,#phi)",kTRUE,
                     36, -0.9, 0.9, AliReducedVarManager::kEta, 90, 0.0, 2.0*TMath::Pi(), AliReducedVarManager::kPhi, 7, -0.5, 6.5, AliReducedVarManager::kITSncls);
        man->AddHistogram(classStr.Data(),"Eta_Phi_ITSnclsShared_prof","ITS <nclusters-shared> vs (#eta,#phi)",kTRUE,
                          36, -0.9, 0.9, AliReducedVarManager::kEta, 90, 0.0, 2.0*TMath::Pi(), AliReducedVarManager::kPhi, 7, -0.5, 6.5, AliReducedVarManager::kITSnclsShared);
	man->AddHistogram(classStr.Data(),"ITSchi2", "ITS #chi^{2}", kFALSE, 200,0.0,50.0, AliReducedVarManager::kITSchi2);
        man->AddHistogram(classStr.Data(),"ITSchi2_RunNo_prof","",kTRUE, runNBins, runHistRange[0], runHistRange[1], AliReducedVarManager::kRunNo, 10, 0., 10., AliReducedVarManager::kITSchi2);
        man->AddHistogram(classStr.Data(),"ITSchi2_ITSncls", "ITS #chi^{2} vs ITS clusters", kFALSE, 100,0.0,100.0, AliReducedVarManager::kITSchi2,7,-0.5,6.5, AliReducedVarManager::kITSncls);
        man->AddHistogram(classStr.Data(),"ITSchi2_TimeFromSOR", "<ITS #chi^{2}> vs time from SOR", kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 200,0.0,20.0, AliReducedVarManager::kITSchi2);
	man->AddHistogram(classStr.Data(),"Eta_Phi_ITSchi2_prof","ITS <#chi^{2}> vs (#eta,#phi)",kTRUE,
                     36, -0.9, 0.9, AliReducedVarManager::kEta, 90, 0.0, 2.0*TMath::Pi(), AliReducedVarManager::kPhi, 100, 0.0, 1000, AliReducedVarManager::kITSchi2);
        man->AddHistogram(classStr.Data(),"ITSncls_Eta_MultVZERO_NTracksTPCoutFromPileup_prof","",kTRUE, 10, 0., 40000., AliReducedVarManager::kVZEROTotalMult, 22,-2000., 20000., AliReducedVarManager::kNTracksTPCoutFromPileup, 20,-1.0,1.0, AliReducedVarManager::kEta, "","","", AliReducedVarManager::kITSncls);
        man->AddHistogram(classStr.Data(),"ITSchi2_Eta_MultVZERO_NTracksTPCoutFromPileup_prof","",kTRUE, 10, 0., 40000., AliReducedVarManager::kVZEROTotalMult, 22,-2000., 20000., AliReducedVarManager::kNTracksTPCoutFromPileup, 20,-1.0,1.0, AliReducedVarManager::kEta, "","","", AliReducedVarManager::kITSchi2);
        man->AddHistogram(classStr.Data(),"ITSnSharedCls_Eta_MultVZERO_NTracksTPCoutFromPileup_prof","",kTRUE, 10, 0., 40000., AliReducedVarManager::kVZEROTotalMult, 22,-2000., 20000., AliReducedVarManager::kNTracksTPCoutFromPileup, 20,-1.0,1.0, AliReducedVarManager::kEta, "","","", AliReducedVarManager::kITSnclsShared);
        
        man->AddHistogram(classStr.Data(),"Chi2TPCConstrainedVsGlobal","",kFALSE,200,0.,200.,AliReducedVarManager::kChi2TPCConstrainedVsGlobal);
        man->AddHistogram(classStr.Data(),"Chi2TPCConstrainedVsGlobal_TimeFromSOR_prof","<Chi2 TPC constrained vs global> vs time from SOR",kTRUE,90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 100,0.,100.,AliReducedVarManager::kChi2TPCConstrainedVsGlobal);
        man->AddHistogram(classStr.Data(),"Chi2TPCConstrainedVsGlobal_Eta_MultVZERO_NTracksTPCoutFromPileup_prof","",kTRUE, 10, 0., 40000., AliReducedVarManager::kVZEROTotalMult, 22,-2000., 20000., AliReducedVarManager::kNTracksTPCoutFromPileup, 20,-1.0,1.0, AliReducedVarManager::kEta, "","","", AliReducedVarManager::kChi2TPCConstrainedVsGlobal);
        man->AddHistogram(classStr.Data(),"Chi2TPCConstrainedVsGlobal_Eta_MultVZERO_NTracksTPCout_prof","",kTRUE, 10, 0., 40000., AliReducedVarManager::kVZEROTotalMult, 20,0., 30000., AliReducedVarManager::kNTracksPerTrackingStatus+AliReducedVarManager::kTPCout, 20,-1.0,1.0, AliReducedVarManager::kEta, "","","", AliReducedVarManager::kChi2TPCConstrainedVsGlobal);
        
        man->AddHistogram(classStr.Data(),"TPCncls","",kFALSE,160,-0.5,159.5,AliReducedVarManager::kTPCncls);
        man->AddHistogram(classStr.Data(),"TPCncls_RunNo_prof","",kTRUE, runNBins, runHistRange[0], runHistRange[1], AliReducedVarManager::kRunNo, 10, 0., 10., AliReducedVarManager::kTPCncls);
        man->AddHistogram(classStr.Data(),"TPCncls_TimeFromSOR","<TPC #cls> vs time from SOR",kTRUE,90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 160,-0.5,159.5,AliReducedVarManager::kTPCncls);
	man->AddHistogram(classStr.Data(),"TPCcrossedRows","", kFALSE, 160,-0.5,159.5,AliReducedVarManager::kTPCcrossedRows);
        man->AddHistogram(classStr.Data(),"TPCfindableClusters","", kFALSE, 160,-0.5,159.5,AliReducedVarManager::kTPCnclsF);
        //man->AddHistogram(classStr.Data(),"TPCfindableClusters_TPCcrossedRows","", kFALSE, 160,-0.5,159.5,AliReducedVarManager::kTPCnclsF, 160, -0.5, 159.5, AliReducedVarManager::kTPCcrossedRows);
        man->AddHistogram(classStr.Data(),"TPCcrossedRowsOverFindableClusters","", kFALSE, 200,0.0,1.5,AliReducedVarManager::kTPCcrossedRowsOverFindableClusters);
        man->AddHistogram(classStr.Data(),"TPCcrossedRowsOverFindableClusters_TimeFromSOR","",kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 100,0.,1.,AliReducedVarManager::kTPCcrossedRowsOverFindableClusters);
        //man->AddHistogram(classStr.Data(),"TPCcrossedRowsOverFindableClusters_TPCnclsSharedRatio","", kFALSE, 200,0.0,1.5,AliReducedVarManager::kTPCcrossedRowsOverFindableClusters, 200, 0.0, 1.0, AliReducedVarManager::kTPCnclsSharedRatio);
	man->AddHistogram(classStr.Data(),"TPCnclsShared","",kFALSE, 160,-0.5,159.5,AliReducedVarManager::kTPCnclsShared);
        man->AddHistogram(classStr.Data(),"TPCnclsSharedRatio","",kFALSE, 200,0.,1.,AliReducedVarManager::kTPCnclsSharedRatio);
        //man->AddHistogram(classStr.Data(),"TPCnclsShared_TPCncls","",kFALSE, 160,-0.5,159.5,AliReducedVarManager::kTPCnclsShared, 160, -0.5, 159.5, AliReducedVarManager::kTPCncls);
        man->AddHistogram(classStr.Data(),"TPCnclsShared_TimeFromSOR","",kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 160,-0.5,159.5,AliReducedVarManager::kTPCnclsShared);
        man->AddHistogram(classStr.Data(),"TPCnclsSharedRatio_TimeFromSOR","",kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 100,0.,1.,AliReducedVarManager::kTPCnclsSharedRatio);
        man->AddHistogram(classStr.Data(),"TPCnclsShared_Eta_MultVZERO_NTracksTPCoutFromPileup_prof","",kTRUE, 10, 0., 40000., AliReducedVarManager::kVZEROTotalMult, 22,-2000., 20000., AliReducedVarManager::kNTracksTPCoutFromPileup, 20,-1.0,1.0, AliReducedVarManager::kEta, "","","", AliReducedVarManager::kTPCnclsShared);
        man->AddHistogram(classStr.Data(),"TPCnclsShared_Eta_MultVZERO_NTracksTPCout_prof","",kTRUE, 10, 0., 40000., AliReducedVarManager::kVZEROTotalMult, 20,0., 30000., AliReducedVarManager::kNTracksPerTrackingStatus+AliReducedVarManager::kTPCout, 20,-1.0,1.0, AliReducedVarManager::kEta, "","","", AliReducedVarManager::kTPCnclsShared);
        man->AddHistogram(classStr.Data(),"TPCncls_Eta_NTracksTPCoutFromPileup_prof","",kTRUE, 20,-1.0,1.0,AliReducedVarManager::kEta,
          43,-1500., 20000., AliReducedVarManager::kNTracksTPCoutFromPileup, 10, 0., 160., AliReducedVarManager::kTPCncls);
        man->AddHistogram(classStr.Data(),"TPCncls_Eta_MultVZERO_NTracksTPCoutFromPileup_prof","",kTRUE, 10, 0., 40000., AliReducedVarManager::kVZEROTotalMult, 22,-2000., 20000., AliReducedVarManager::kNTracksTPCoutFromPileup, 20,-1.0,1.0, AliReducedVarManager::kEta, "","","", AliReducedVarManager::kTPCncls);
        man->AddHistogram(classStr.Data(),"TPCncls_Eta_MultVZERO_NTracksTPCout_prof","",kTRUE, 10, 0., 40000., AliReducedVarManager::kVZEROTotalMult, 20,0., 30000., AliReducedVarManager::kNTracksPerTrackingStatus+AliReducedVarManager::kTPCout, 20,-1.0,1.0, AliReducedVarManager::kEta, "","","", AliReducedVarManager::kTPCncls);
        
        
        man->AddHistogram(classStr.Data(),"TPCchi2","",kFALSE, 200,0.0,10.0,AliReducedVarManager::kTPCchi2);
        man->AddHistogram(classStr.Data(),"TPCchi2_nTPCout","",kFALSE, 100, 0.0, 15000., AliReducedVarManager::kNTracksPerTrackingStatus+AliReducedVarManager::kTPCout, 200,0.0,10.0,AliReducedVarManager::kTPCchi2);
        man->AddHistogram(classStr.Data(),"TPCchi2_Pt","",kFALSE, 200,0.0,10.0,AliReducedVarManager::kTPCchi2, 100, 0., 10.0, AliReducedVarManager::kPt);
        man->AddHistogram(classStr.Data(),"TPCchi2_RunNo_prof","",kTRUE, runNBins, runHistRange[0], runHistRange[1], AliReducedVarManager::kRunNo, 10, 0., 10., AliReducedVarManager::kTPCchi2);
        man->AddHistogram(classStr.Data(),"TPCchi2_TPCncls", "TPC #chi^{2} vs TPC clusters", kFALSE, 100,0.0,10.0, AliReducedVarManager::kTPCchi2,160,0.,160., AliReducedVarManager::kTPCncls);
        man->AddHistogram(classStr.Data(),"TPCsegments","",kFALSE, 9,0.0,9.0,AliReducedVarManager::kTPCNclusBitsFired);
        man->AddHistogram(classStr.Data(),"TPCsegments_TimeFromSOR","",kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 8,0.0,8.0,AliReducedVarManager::kTPCNclusBitsFired);
        man->AddHistogram(classStr.Data(),"TPCsegments_Eta_MultVZERO_NTracksTPCoutFromPileup_prof","",kTRUE, 10, 0., 40000., AliReducedVarManager::kVZEROTotalMult, 22,-2000., 20000., AliReducedVarManager::kNTracksTPCoutFromPileup, 20,-1.0,1.0, AliReducedVarManager::kEta, "","","", AliReducedVarManager::kTPCNclusBitsFired);
        man->AddHistogram(classStr.Data(),"TPCchi2_TimeFromSOR","TPC <#chi^{2}> vs time from SOR",kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 200,0.0,10.0,AliReducedVarManager::kTPCchi2);
        man->AddHistogram(classStr.Data(),"TPCchi2_Eta_MultVZERO_NTracksTPCoutFromPileup_prof","",kTRUE, 10, 0., 40000., AliReducedVarManager::kVZEROTotalMult, 22,-2000., 20000., AliReducedVarManager::kNTracksTPCoutFromPileup, 20,-1.0,1.0, AliReducedVarManager::kEta, "","","", AliReducedVarManager::kTPCchi2);
        man->AddHistogram(classStr.Data(),"TPCchi2_Eta_MultVZERO_NTracksTPCout_prof","",kTRUE, 10, 0., 40000., AliReducedVarManager::kVZEROTotalMult, 30,0., 30000., AliReducedVarManager::kNTracksPerTrackingStatus+AliReducedVarManager::kTPCout, 20,-1.0,1.0, AliReducedVarManager::kEta, "","","", AliReducedVarManager::kTPCchi2);
        man->AddHistogram(classStr.Data(),"Eta_Phi_TPCncls_prof","TPC <nclusters> vs (#eta,#phi)",kTRUE,
                     36, -0.9, 0.9, AliReducedVarManager::kEta, 90, 0.0, 2.0*TMath::Pi(), AliReducedVarManager::kPhi, 160, -0.5, 159.5, AliReducedVarManager::kTPCncls);    
	man->AddHistogram(classStr.Data(),"Eta_Phi_TPCcrossedRows_prof","TPC <n crossed rows> vs (#eta,#phi)",kTRUE,
                     36, -0.9, 0.9, AliReducedVarManager::kEta, 90, 0.0, 2.0*TMath::Pi(), AliReducedVarManager::kPhi, 160, -0.5, 159.5, AliReducedVarManager::kTPCcrossedRows);
	man->AddHistogram(classStr.Data(),"Eta_Phi_TPCnclsF_prof","TPC <nclusters findable> vs (#eta,#phi)",kTRUE,
                     36, -0.9, 0.9, AliReducedVarManager::kEta, 90, 0.0, 2.0*TMath::Pi(), AliReducedVarManager::kPhi, 160, -0.5, 159.5, AliReducedVarManager::kTPCnclsF);
	man->AddHistogram(classStr.Data(),"Eta_Phi_TPCchi2_prof","TPC <#chi^{2}> vs (#eta,#phi)",kTRUE,
                     36, -0.9, 0.9, AliReducedVarManager::kEta, 90, 0.0, 2.0*TMath::Pi(), AliReducedVarManager::kPhi, 160, -0.5, 159.5, AliReducedVarManager::kTPCchi2);
        man->AddHistogram(classStr.Data(),"Eta_Phi_TPCchi2","TPC #chi^{2} vs (#eta,#phi)",kFALSE,
                          36, -0.9, 0.9, AliReducedVarManager::kEta, 90, 0.0, 2.0*TMath::Pi(), AliReducedVarManager::kPhi, 100, 0., 10.0, AliReducedVarManager::kTPCchi2);
        man->AddHistogram(classStr.Data(),"TPCsignal_Pin","TPC dE/dx vs. inner param P",kFALSE,
                     300,0.0,30.0,AliReducedVarManager::kPin,200,-0.5,199.5,AliReducedVarManager::kTPCsignal);
	man->AddHistogram(classStr.Data(),"TPCsignalN","",kFALSE,160,-0.5,159.5,AliReducedVarManager::kTPCsignalN);
        man->AddHistogram(classStr.Data(),"TPCsignalN_TimeFromSOR","TPC <#cls pid> vs time from SOR",kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 160,-0.5,159.5,AliReducedVarManager::kTPCsignalN);
	man->AddHistogram(classStr.Data(),"Eta_Phi_TPCsignalN_prof","TPC <nclusters pid> vs (#eta,#phi)",kTRUE,
                     36, -0.9, 0.9, AliReducedVarManager::kEta, 90, 0.0, 2.0*TMath::Pi(), AliReducedVarManager::kPhi, 160, -0.5, 159.5, AliReducedVarManager::kTPCsignalN);    
        man->AddHistogram(classStr.Data(),"TPCnsigElectron_TimeFromSOR","TPC N_{#sigma} electron vs. time from SOR",kFALSE,
                          90, 0., 450.,AliReducedVarManager::kTimeRelativeSOR,100,-5.0,5.0,AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron);
        man->AddHistogram(classStr.Data(),"TPCnsigElectron_TimeFromSOR_prof","<TPC N_{#sigma} electron> vs. time from SOR",kTRUE,
                          90, 0., 450.,AliReducedVarManager::kTimeRelativeSOR,100,-5.0,5.0,AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron);
        man->AddHistogram(classStr.Data(),"TPCnsigElectron_Pin","TPC N_{#sigma} electron vs. inner param P",kFALSE,
                     100,0.0,10.0,AliReducedVarManager::kPin,100,-5.0,5.0,AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron);
        /*man->AddHistogram(classStr.Data(),"TPCnsigElectron_Pin_TimeFromSOR","TPC N_{#sigma} electron vs. inner param P vs time from SOR",kFALSE,
                          100,0.0,10.0,AliReducedVarManager::kPin,100,-5.0,5.0,AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron, 15, 0., 450., AliReducedVarManager::kTimeRelativeSOR);*/
        man->AddHistogram(classStr.Data(),"TPCnsigElectron_Eta","TPC N_{#sigma} electron vs. #eta",kFALSE,
                          36,-0.9,0.9,AliReducedVarManager::kEta,100,-5.0,5.0,AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron);
        /*man->AddHistogram(classStr.Data(),"TPCnsigElectron_Eta_TimeFromSOR","TPC N_{#sigma} electron vs. #eta vs time from SOR",kFALSE,
                          36,-0.9,0.9,AliReducedVarManager::kEta,100,-5.0,5.0,AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron, 15, 0., 450., AliReducedVarManager::kTimeRelativeSOR);*/
        man->AddHistogram(classStr.Data(),"TPCnSigEle_Eta_MultVZERO_NTracksTPCoutFromPileup_prof","",kTRUE, 10, 0., 40000., AliReducedVarManager::kVZEROTotalMult, 22,-2000., 20000., AliReducedVarManager::kNTracksTPCoutFromPileup, 20,-1.0,1.0, AliReducedVarManager::kEta, "","","", AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron);
        man->AddHistogram(classStr.Data(),"TPCnSigEle_Eta_MultVZERO_NTracksTPCout_prof","",kTRUE, 10, 0., 40000., AliReducedVarManager::kVZEROTotalMult, 20,0., 30000., AliReducedVarManager::kNTracksPerTrackingStatus+AliReducedVarManager::kTPCout, 20,-1.0,1.0, AliReducedVarManager::kEta, "","","", AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron);
        man->AddHistogram(classStr.Data(),"TPCnsigElectronCorrected_Pin","Corrected TPC N_{#sigma} electron vs. inner param P",kFALSE,
                          100,0.0,10.0,AliReducedVarManager::kPin,100,-5.0,5.0,AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron);
        
        /*
        man->AddHistogram(classStr.Data(),"TPCnsigElectronCorrected_Pin_DCAxy","DCAxy vs corrected TPC N_{#sigma} electron vs. inner param P",kFALSE,
                          50,0.0,5.0,AliReducedVarManager::kPin,50,-5.0,5.0,AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron, 80, -2., 2., AliReducedVarManager::kDcaXY);
        man->AddHistogram(classStr.Data(),"TPCnsigElectronCorrected_Pin_DCAz","DCAz vs corrected TPC N_{#sigma} electron vs. inner param P",kFALSE,
                          50,0.0,5.0,AliReducedVarManager::kPin,50,-5.0,5.0,AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron, 120, -3., 3., AliReducedVarManager::kDcaZ);
        man->AddHistogram(classStr.Data(),"TPCnsigElectronCorrected_Pin_ITSchi2","ITS chi2 vs corrected TPC N_{#sigma} electron vs. inner param P",kFALSE,
                          50,0.0,5.0,AliReducedVarManager::kPin,50,-5.0,5.0,AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron, 50, 0., 50., AliReducedVarManager::kITSchi2);
        man->AddHistogram(classStr.Data(),"TPCnsigElectronCorrected_Pin_ITSchi2_prof","ITS <chi2> vs corrected TPC N_{#sigma} electron vs. inner param P",kTRUE,
                          50,0.0,5.0,AliReducedVarManager::kPin,50,-5.0,5.0,AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron, 50, 0., 50., AliReducedVarManager::kITSchi2);
        man->AddHistogram(classStr.Data(),"TPCnsigElectronCorrected_Pin_TPCchi2","TPC chi2 vs corrected TPC N_{#sigma} electron vs. inner param P",kFALSE,
                          50,0.0,5.0,AliReducedVarManager::kPin,50,-5.0,5.0,AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron, 50, 0., 5., AliReducedVarManager::kTPCchi2);
        man->AddHistogram(classStr.Data(),"TPCnsigElectronCorrected_Pin_TPCchi2_prof","TPC <chi2> vs corrected TPC N_{#sigma} electron vs. inner param P",kTRUE,
                          50,0.0,5.0,AliReducedVarManager::kPin,50,-5.0,5.0,AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron, 50, 0., 5., AliReducedVarManager::kTPCchi2);
        man->AddHistogram(classStr.Data(),"TPCnsigElectronCorrected_Pin_Chi2TPCConstrainedVsGlobal","TPC constrained vs global chi2 vs corrected TPC N_{#sigma} electron vs. inner param P",kFALSE,
                          50,0.0,5.0,AliReducedVarManager::kPin,50,-5.0,5.0,AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron, 50, 0., 200., AliReducedVarManager::kChi2TPCConstrainedVsGlobal);
        man->AddHistogram(classStr.Data(),"TPCnsigElectronCorrected_Pin_Chi2TPCConstrainedVsGlobal_prof","TPC constrained vs global <chi2> vs corrected TPC N_{#sigma} electron vs. inner param P",kTRUE,
                          50,0.0,5.0,AliReducedVarManager::kPin,50,-5.0,5.0,AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron, 50, 0., 200., AliReducedVarManager::kChi2TPCConstrainedVsGlobal);
        man->AddHistogram(classStr.Data(),"TPCnsigElectronCorrected_Pin_NclsITS_prof","<#ITScls> vs corrected TPC N_{#sigma} electron vs. inner param P",kTRUE,
                          50,0.0,5.0,AliReducedVarManager::kPin,50,-5.0,5.0,AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron, 7, 0., 7., AliReducedVarManager::kITSncls);
        man->AddHistogram(classStr.Data(),"TPCnsigElectronCorrected_Pin_NclsTPC_prof","<#TPCcls> vs corrected TPC N_{#sigma} electron vs. inner param P",kTRUE,
                          50,0.0,5.0,AliReducedVarManager::kPin,50,-5.0,5.0,AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron, 80, 0., 160., AliReducedVarManager::kTPCncls);
        */
        
        man->AddHistogram(classStr.Data(),"TPCnsigElectronCorrected_Eta","TPC N_{#sigma} electron vs. #eta",kFALSE,
                          36,-0.9,0.9,AliReducedVarManager::kEta,100,-5.0,5.0,AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron);
        man->AddHistogram(classStr.Data(),"TPCnsigElectronCorrected_TimeFromSOR_prof","<corrected TPC N_{#sigma} electron> vs. time from SOR",kTRUE,
                          90, 0., 450.,AliReducedVarManager::kTimeRelativeSOR,80,-4.0,4.0,AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron);
        
        /*man->AddHistogram(classStr.Data(),"TPCnsigElectronCorrected_Pin_MassUsedForTracking_prof","<mass used for tracking> vs corrected TPC N_{#sigma} electron vs. inner param P",kTRUE,
                          100,0.0,10.0,AliReducedVarManager::kPin,100,-5.0,5.0,AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron, 40, 0., 4., AliReducedVarManager::kMassUsedForTracking);
        man->AddHistogram(classStr.Data(),"TPCnsigElectronCorrected_Eta_MassUsedForTracking_prof","<mass used for tracking> vs corrected TPC N_{#sigma} electron vs. eta",kTRUE,
                          36,-0.9,0.9,AliReducedVarManager::kEta,100,-5.0,5.0,AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron, 40, 0., 4., AliReducedVarManager::kMassUsedForTracking);*/
	man->AddHistogram(classStr.Data(),"TPCnsigElectron_Run","TPC N_{#sigma} electron vs. run",kTRUE,
                     runNBins, runHistRange[0], runHistRange[1], AliReducedVarManager::kRunNo,100,-5.0,5.0,AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron);
      
        man->AddHistogram(classStr.Data(),"TOFbeta_P","TOF #beta vs P",kFALSE,
                     200,0.0,20.0,AliReducedVarManager::kP, 220,0.0,1.1,AliReducedVarManager::kTOFbeta);
        
        man->AddHistogram(classStr.Data(),"TPCchi2_ITSchi2","",kFALSE, 100,0.0,10.0,AliReducedVarManager::kTPCchi2, 100,0.0,40.0,AliReducedVarManager::kITSchi2);
        man->AddHistogram(classStr.Data(),"ITSchi2_Chi2TPCconstrainedVsGlobal","",kFALSE, 500,0.0,1000.0,AliReducedVarManager::kChi2TPCConstrainedVsGlobal, 100,0.0,40.0,AliReducedVarManager::kITSchi2);
        man->AddHistogram(classStr.Data(),"TPCchi2_Chi2TPCconstrainedVsGlobal","",kFALSE, 500,0.0,1000.0,AliReducedVarManager::kChi2TPCConstrainedVsGlobal, 100,0.0,10.0,AliReducedVarManager::kTPCchi2);
        
        if(classStr.Contains("MCTruth")) {
          man->AddHistogram(classStr.Data(), "PtMC", "p_{T} MC", kFALSE, 150, 0., 15.0, AliReducedVarManager::kPtMC);
          man->AddHistogram(classStr.Data(), "PtRec_PtMC", "p_{T} MC vs p_{T} reconstructed", kFALSE, 150, 0., 15.0, AliReducedVarManager::kPtMC, 150, 0., 15.0, AliReducedVarManager::kPt);
          man->AddHistogram(classStr.Data(), "PhiMC", "#varphi MC", kFALSE, 180, 0., 6.3, AliReducedVarManager::kPhiMC);
          man->AddHistogram(classStr.Data(), "PhiRec_PhiMC", "#varphi MC vs #varphi reconstructed", kFALSE, 180, 0., 6.3, AliReducedVarManager::kPhiMC, 180, 0., 6.3, AliReducedVarManager::kPhi);
          man->AddHistogram(classStr.Data(), "EtaMC", "#eta MC", kFALSE, 100, -1.0, 1.0, AliReducedVarManager::kEtaMC);
          man->AddHistogram(classStr.Data(), "EtaRec_EtaMC", "#eta MC vs #eta reconstructed", kFALSE, 100, -1.0, 1.0, AliReducedVarManager::kEtaMC, 100, -1.0, 1.0, AliReducedVarManager::kEta);          
          man->AddHistogram(classStr.Data(), "PDGcode0", "PDG code of the track", kFALSE, 12000, -6000., 6000.0, AliReducedVarManager::kPdgMC);
          man->AddHistogram(classStr.Data(), "PDGcode1", "PDG code of the track's mother", kFALSE, 12000, -6000., 6000.0, AliReducedVarManager::kPdgMC+1);
          man->AddHistogram(classStr.Data(), "PDGcode2", "PDG code of the track's grand-mother", kFALSE, 12000, -6000., 6000.0, AliReducedVarManager::kPdgMC+2);
          man->AddHistogram(classStr.Data(), "PDGcode3", "PDG code of the track's grand-grand mother", kFALSE, 12000, -6000., 6000.0, AliReducedVarManager::kPdgMC+3);
        }
      continue;
    }  // end if "TrackQA"
        
    const Int_t kNMassBins = 101;
    Double_t massBins[kNMassBins];
    for(Int_t i=0; i<kNMassBins; ++i) massBins[i] = 1.0 + i*0.04; 
    
    //const Int_t kNPtBins = 55;
    const Int_t kNPtBins = 7;
    //const Int_t kNPtBins = 2;
    //Double_t ptBins[kNPtBins];
    /*Double_t ptBins[kNPtBins] = {
      0.0, 0.02, 0.04, 0.06, 0.08, 
      0.1, 0.12, 0.14, 0.16, 0.18,
      0.2, 0.25, 0.30, 0.35, 0.40,
      0.45, 0.5, 1.0, 3.0, 5.0,
      7.0, 10.0, 20.0
    };*/
    /*Double_t ptBins[kNPtBins] = {
       0.0, 0.02, 0.04, 0.06, 0.08, 
       0.1, 0.12, 0.14, 0.16, 0.18,
       0.3, 0.5, 1.0, 3.0, 5.0,
       7.0, 10.0, 20.0
    };*/
    Double_t ptBins[kNPtBins] = {
       0.0, 0.15, 1.3, 3.0, 5.0, 10.0,50. 
    };
    //for(Int_t i=0; i<=25; ++i) ptBins[i] = 0.0+ i*0.02;
    //ptBins[26] = 1.5; ptBins[27] = 4.5; ptBins[28] = 10.0; ptBins[29] = 20.0;
    
   const Int_t kNCentBins = 27;
    
   Double_t centBins2[kNCentBins] = {
      0.0, 1.0, 2.0, 3.0, 4.0, 
      5.0, 6.0, 7.0, 8.0, 9.0, 
      10., 12.0, 14.0, 16.0, 18.0, 
      20., 22.5, 25.0, 27.5, 30.0, 
      35., 40., 50., 60., 70., 
      80., 90.
   };
   
   const Int_t kNtpcOutBins = 26;
   const Int_t kNtpcOutBinsMC = 2;
   Double_t ntpcOutLims[kNtpcOutBins] = {
      0.0, 1000., 2000., 3000., 4000., 
      5000., 6000., 7000., 8000., 8500., 
      9000., 9500., 10000., 10500., 11000., 
      11500., 12000., 12250., 12500., 12750., 
      13000., 13250., 13500., 13750., 14000., 
      14250.
   };
   //for(Int_t i=0; i<kNtpcOutBins;++i) ntpcOutLims[i] = 0.0 + i*200.;
   Double_t ntpcOutLimsMC[kNtpcOutBinsMC] = {0., 15000.};
   
    const Int_t kNVtxBins = 2;
    //Double_t vtxBins[kNVtxBins] = {-10.,-9.,-8.,-7.,-5.,-3.,-1.,1.,3.,5.,7.,8.,9.,10.};
    Double_t vtxBins[kNVtxBins] = {-10.,10.};
    const Int_t kNTimeBins = 2;
    Double_t timeBins[kNTimeBins];
    for(Int_t i=0; i<kNTimeBins; ++i) timeBins[i] = 0.0 + i*450.0;
    
    const Int_t kNEPbins = 11;
    //Double_t epBins[kNEPbins] = {-0.5*TMath::Pi(),0.5*TMath::Pi()};
    Double_t epBins[kNEPbins];
    for(Int_t i=0;i<=10;++i) epBins[i] = -0.5*TMath::Pi()+i*TMath::Pi()/10.;
    
    //const Int_t kNTPCoutLims = 2;
    //Double_t tpcOutLims[kNTPCoutLims] = {-50000., 50000.};
    /*   -2000.,-700.,-400., -200., -100.,
       -50., 0.,   100.,   200.,   350.,   
       500., 750., 1000., 1500., 2000., 
       2500., 3000., 3500., 4000., 5000., 
       6000., 8000., 10000., 12000., 15000., 
       20000., 25000.};*/
    
    //Int_t vars[5] = {AliReducedVarManager::kMass, AliReducedVarManager::kPt, 
    //   AliReducedVarManager::kCentVZERO, AliReducedVarManager::kVtxZ, AliReducedVarManager::kVZERORP+6+1};
    //Int_t vars[5] = {AliReducedVarManager::kMass, AliReducedVarManager::kPt, 
    //   AliReducedVarManager::kCentVZERO, AliReducedVarManager::kVtxZ, AliReducedVarManager::kNTracksPerTrackingStatus+AliReducedVarManager::kTPCout};
    //Int_t vars[5] = {AliReducedVarManager::kMass, AliReducedVarManager::kPt, 
     //     AliReducedVarManager::kCentVZERO, AliReducedVarManager::kVtxZ, AliReducedVarManager::kVZERORP+6+1};
     Int_t vars[5] = {AliReducedVarManager::kMass, AliReducedVarManager::kPt, AliReducedVarManager::kCentVZERO,
        AliReducedVarManager::kNTracksPerTrackingStatus+AliReducedVarManager::kTPCout, AliReducedVarManager::kVZERORP+6+1};
       
    //TArrayD pairHistBinLimits[5] = {TArrayD(kNMassBins,massBins), TArrayD(kNPtBins,ptBins), TArrayD(kNCentBins,centBins), TArrayD(kNVtxBins,vtxBins), TArrayD(kNEPbins,epBins)};
    TArrayD pairHistBinLimits[5];
    pairHistBinLimits[0] = TArrayD(kNMassBins,massBins);
    pairHistBinLimits[1] = TArrayD(kNPtBins,ptBins);
    pairHistBinLimits[2] = TArrayD(kNCentBins,centBins2);
    //pairHistBinLimits[3] = TArrayD(kNVtxBins,vtxBins);
    //pairHistBinLimits[3] = TArrayD(kNTimeBins,timeBins);
    //pairHistBinLimits[4] = TArrayD(kNEPbins,epBins);
    pairHistBinLimits[3] = TArrayD(kNtpcOutBins,ntpcOutLims);
    /*if(task->GetRunOverMC())
       pairHistBinLimits[4] = TArrayD(kNtpcOutBinsMC,ntpcOutLimsMC);
    else
      pairHistBinLimits[4] = TArrayD(kNtpcOutBins,ntpcOutLims);*/
    pairHistBinLimits[4] = TArrayD(kNEPbins,epBins);
    
    // Histograms for pairs
    if(classStr.Contains("Pair")) {
      man->AddHistClass(classStr.Data());
      cout << classStr.Data() << endl;
      man->AddHistogram(classStr.Data(), "PairInvMass", "Differential pair inv. mass", 5, vars, pairHistBinLimits);
      man->AddHistogram(classStr.Data(), "MeanPt_massCent", "<p_{T}> vs (mass,cent)", kTRUE, 100, 1.0, 5.0, AliReducedVarManager::kMass, 9, 0.0, 90., AliReducedVarManager::kCentVZERO, 10, 0., 100., AliReducedVarManager::kPt);
      man->AddHistogram(classStr.Data(), "MeanPt2_massCent", "<p^{2}_{T}> vs (mass,cent)", kTRUE, 100, 1.0, 5.0, AliReducedVarManager::kMass, 9, 0.0, 90., AliReducedVarManager::kCentVZERO, 10, 0., 100., AliReducedVarManager::kPtSquared);
      if(classStr.Contains("PairSE") || classStr.Contains("PairPrefilterSE")) {
        man->AddHistogram(classStr.Data(), "PairType", "Pair type", kFALSE, 4, -0.5, 3.5, AliReducedVarManager::kPairType);
	man->AddHistogram(classStr.Data(), "Mass", "Invariant mass", kFALSE, kNMassBins, 0.0, 5.0, AliReducedVarManager::kMass);
        man->AddHistogram(classStr.Data(), "Mass_TimeFromSOR", "Invariant mass vs time from SOR", kFALSE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, kNMassBins, 0.0, 5.0, AliReducedVarManager::kMass);
        man->AddHistogram(classStr.Data(), "Mass_TimeFromSOR_prof", "<Invariant mass> vs time from SOR", kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, kNMassBins, 0.0, 5.0, AliReducedVarManager::kMass);
        man->AddHistogram(classStr.Data(), "Mass_Cent_TimeFromSOR_prof", "<Invariant mass> vs centrality and time from SOR", kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 18, 0.0, 90., AliReducedVarManager::kCentVZERO, kNMassBins, 0.0, 5.0, AliReducedVarManager::kMass);
        man->AddHistogram(classStr.Data(), "Pt", "", kFALSE, 1000, 0.0, 10.0, AliReducedVarManager::kPt);
        man->AddHistogram(classStr.Data(), "Pt_coarse", "", kFALSE, 20, 0.0, 20.0, AliReducedVarManager::kPt);
	man->AddHistogram(classStr.Data(), "Rapidity", "Rapidity", kFALSE, 240, -1.2, 1.2, AliReducedVarManager::kRap);
	man->AddHistogram(classStr.Data(), "Phi", "Azimuthal distribution", kFALSE, 315, 0.0, 6.3, AliReducedVarManager::kPhi);
        man->AddHistogram(classStr.Data(), "Leg1TPCchi2_Leg2TPCchi2", "", kFALSE, 100, 0.0, 6.0, AliReducedVarManager::kPairLegTPCchi2, 100, 0.0, 6.0, AliReducedVarManager::kPairLegTPCchi2+1);
        man->AddHistogram(classStr.Data(), "Leg1ITSchi2_Leg2ITSchi2", "", kFALSE, 100, 0.0, 6.0, AliReducedVarManager::kPairLegITSchi2, 100, 0.0, 6.0, AliReducedVarManager::kPairLegITSchi2+1);
        man->AddHistogram(classStr.Data(), "CosThetaStarCS", "cos(#theta^{*})_{CS}", kFALSE, 22, -1.1, 1.1, AliReducedVarManager::kPairThetaCS);
        man->AddHistogram(classStr.Data(), "CosThetaStarCS_ptMC", "cos(#theta^{*})_{CS} vs MC pt", kFALSE, 22, -1.1, 1.1, AliReducedVarManager::kPairThetaCS, 50, 0., 1., AliReducedVarManager::kPtMC);
        man->AddHistogram(classStr.Data(), "CosThetaStarHE", "cos(#theta^{*})_{HE}", kFALSE, 22, -1.1, 1.1, AliReducedVarManager::kPairThetaHE);
        man->AddHistogram(classStr.Data(), "CosThetaStarHE_pt", "cos(#theta^{*})_{HE} vs pt", kFALSE, 22, -1.1, 1.1, AliReducedVarManager::kPairThetaHE, 100, 0., 1., AliReducedVarManager::kPt);
        man->AddHistogram(classStr.Data(), "CosThetaStarHE_pt_coarse", "cos(#theta^{*})_{HE} vs pt", kFALSE, 10, -1.0, 1.0, AliReducedVarManager::kPairThetaHE, 20, 0., 20., AliReducedVarManager::kPt);
        man->AddHistogram(classStr.Data(), "PhiStarCS", "#varphi^{*}_{CS}", kFALSE, 22, -3.3, 3.3, AliReducedVarManager::kPairPhiCS);
        man->AddHistogram(classStr.Data(), "PhiStarHE", "#varphi^{*}_{HE}", kFALSE, 22, -3.3, 3.3, AliReducedVarManager::kPairPhiHE);         
        Double_t v2PtLims[9] = {0.,0.15,0.3,0.5,1.0,3.0,5.0,7.0,10.0};
        Double_t v2MassLims[8] = {2.4,2.6,2.8,2.92,3.04,3.16,3.50,4.0};
        Double_t v2CentLims[4] = {30.,50.,70.,90.};
        //man->AddHistogram(classStr.Data(), "v2VZEROA_massPtCent", "v_{2}^{VZERO-A} (mass,p_{T})", kTRUE, 7, v2MassLims, AliReducedVarManager::kMass, 8, v2PtLims, AliReducedVarManager::kPt, 3, v2CentLims, AliReducedVarManager::kCentVZERO, "", "", "", AliReducedVarManager::kVZEROFlowVn+0*6+1);
        //man->AddHistogram(classStr.Data(), "v2VZEROC_massPtCent", "v_{2}^{VZERO-C} (mass,p_{T})", kTRUE, 7, v2MassLims, AliReducedVarManager::kMass, 8, v2PtLims, AliReducedVarManager::kPt, 3, v2CentLims, AliReducedVarManager::kCentVZERO, "", "", "", AliReducedVarManager::kVZEROFlowVn+1*6+1);
      }   // end if "QA"
      
      if(classStr.Contains("MCTruth")) {
         man->AddHistogram(classStr.Data(), "PtMC", "p_{T} MC", kFALSE, 1000, 0., 10.0, AliReducedVarManager::kPtMC);
         man->AddHistogram(classStr.Data(), "PtMC_coarse", "p_{T} MC", kFALSE, 20, 0., 20.0, AliReducedVarManager::kPtMC);
         man->AddHistogram(classStr.Data(), "MassMC_Mass", "Invariant mass, MC vs reconstructed", kFALSE, 150, 2.0, 3.5, AliReducedVarManager::kMass,
            150, 2.0, 3.5, AliReducedVarManager::kMassMC);
         man->AddHistogram(classStr.Data(), "PtMC_Pt", "pair pT, MC vs reconstructed", kFALSE, 150, 0.0, 15., AliReducedVarManager::kPt,
                           150, 0.0, 15.0, AliReducedVarManager::kPtMC);
         man->AddHistogram(classStr.Data(), "PtMC_MassMC", "", kFALSE, 100, 0.0, 1., AliReducedVarManager::kPtMC,
                           30, 2.7, 3.3, AliReducedVarManager::kMassMC);
         man->AddHistogram(classStr.Data(), "Pt_Mass", "", kFALSE, 100, 0.0, 1., AliReducedVarManager::kPt,
                           30, 2.7, 3.3, AliReducedVarManager::kMass);
         man->AddHistogram(classStr.Data(), "PtMC_Mass", "", kFALSE, 100, 0.0, 1., AliReducedVarManager::kPtMC,
                           30, 2.7, 3.3, AliReducedVarManager::kMass);
         man->AddHistogram(classStr.Data(), "Pt_MassMC", "", kFALSE, 100, 0.0, 1., AliReducedVarManager::kPt,
                           30, 2.7, 3.3, AliReducedVarManager::kMassMC);
         man->AddHistogram(classStr.Data(), "PtMC_Pt_lowPtZoom", "pair pT, MC vs reconstructed", kFALSE, 100, 0.0, 0.5, AliReducedVarManager::kPt,
                           100, 0.0, 0.5, AliReducedVarManager::kPtMC);
         man->AddHistogram(classStr.Data(), "PtMC_Pt_Mass", "pair pT, MC vs reconstructed, vs mass", kFALSE, 100, 0.0, 0.5, AliReducedVarManager::kPt,
                           100, 0.0, 0.5, AliReducedVarManager::kPtMC, 15, 2.72, 3.32, AliReducedVarManager::kMass);
         man->AddHistogram(classStr.Data(), "PtMC_Pt_Mass_coarse", "pair pT, MC vs reconstructed, vs mass", kFALSE, 25, 0.0, 0.5, AliReducedVarManager::kPt,
                           25, 0.0, 0.5, AliReducedVarManager::kPtMC, 15, 2.72, 3.32, AliReducedVarManager::kMass);
         
         man->AddHistogram(classStr.Data(), "Mass_Pt_CentVZERO", "", kFALSE, kNMassBins-1, massBins, AliReducedVarManager::kMass, kNPtBins-1, ptBins, AliReducedVarManager::kPt, kNCentBins-1, centBins, AliReducedVarManager::kCentVZERO);
      }
      
      continue;
    }   // end if for Pair classes of histograms
  }  // end loop over histogram classes
}
