//
//
void SetDefaults(AliAnalysisTaskEbyeIterPID *defaultTask);
TTree *GetLookUpTable(Bool_t runOnGrid, Int_t index);
//
//
AliAnalysisTaskEbyeIterPID* Config_marsland_EbyeIterPID(Bool_t getFromAlien, Int_t settingType, Int_t lhcPeriod, Int_t lookUpTableIndex, TString combinedName) {
  //
  // Configuration file for the AliAnalysisTaskEbyeIterPID.cxx class
  //
  AliAnalysisTaskEbyeIterPID *task = new AliAnalysisTaskEbyeIterPID(combinedName);
  SetDefaults(task);
  if (lhcPeriod==1) task->SelectCollisionCandidates(AliVEvent::kMB);   // select minimum bias events for LHC10h
  if (lhcPeriod==2) task->SelectCollisionCandidates(AliVEvent::kINT7); // select minimum bias events for LHC15o
  //
  // Get the lookup table
  TTree *lookUpTree=NULL;
  if (lookUpTableIndex>0 && settingType>199) lookUpTree = GetLookUpTable(getFromAlien,lookUpTableIndex);

  std::cout << " Info::marsland: ===== In the Config --> Running with lhcPeriod =  ";
  std::cout << lhcPeriod << " --- lookUpTableIndex = " << lookUpTableIndex << " --- settingType = " << settingType << std::endl;
  // Other Specific settings

  switch (settingType) {
    //
    // ====================================================================================
    // ============================= Real Data Settings ===================================
    // ====================================================================================
    //
    case 0:{
      std::cout << settingType << " Info::marsland: (Default event & track cuts) + allCuts + ArmPodTree filled " << std::endl;
      task->SetDefaultTrackCuts(kTRUE);
      task->SetDefaultEventCuts(kTRUE);
      task->SetFillAllCutVariables(kTRUE);
      task->SetFillArmPodTree(kTRUE);
      task->SetRunOnGrid(kTRUE);
      task->fEventCuts.fUseVariablesCorrelationCuts = true;
    }
    break;
    case 1:{
      std::cout << settingType << " Info::marsland: (Open event & track cuts) + allCuts + ArmPodTree filled  " << std::endl;
      task->SetDefaultTrackCuts(kFALSE);
      task->SetDefaultEventCuts(kFALSE);
      task->SetFillAllCutVariables(kTRUE);
      task->SetFillArmPodTree(kTRUE);
      task->SetRunOnGrid(kTRUE);
      task->fEventCuts.fUseVariablesCorrelationCuts = true;
    }
    break;
    case 2:{
      std::cout << settingType << " Info::marsland: (REFERENCE settings) + centBinning 5 " << std::endl;
      task->SetUseCouts(kFALSE);
      task->SetIsMCtrue(kFALSE);
      task->SetFillAllCutVariables(kTRUE);
      task->SetFillArmPodTree(kTRUE);
      task->fEventCuts.fUseVariablesCorrelationCuts = true;
      const Int_t tmpCentbins = 17;
      Float_t tmpfxCentBins[tmpCentbins] = {0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80};
      task->SetCentralityBinning(tmpCentbins,tmpfxCentBins);
    }
    break;
    case 3:{
      std::cout << settingType << " Info::marsland: (REFERENCE settings) + centBinning 10 " << std::endl;
      task->SetUseCouts(kFALSE);
      task->SetFillAllCutVariables(kTRUE);
      task->SetFillArmPodTree(kTRUE);
      task->fEventCuts.fUseVariablesCorrelationCuts = true;
      const Int_t tmpCentbins = 9;
      Float_t tmpfxCentBins[tmpCentbins] = {0,10,20,30,40,50,60,70,80};
      task->SetCentralityBinning(tmpCentbins,tmpfxCentBins);
    }
    break;
    case 4:{
      std::cout << settingType << " Info::marsland: Fill hists for all syst setting " << std::endl;
      // Real data settings
      task->SetNEtabins(16);
      task->SetEtaLowerEdge(-0.8.);
      task->SetEtaUpperEdge( 0.8.);
      task->SetNMomBins(600);
      task->SetMomLowerEdge(0.);
      task->SetMomUpperEdge(12.);
      task->SetDeDxBinWidth(1);
      task->SetDeDxLowerEdge(20.);
      task->SetDeDxUpperEdge(1020.);
      task->SetFillAllCutVariables(kTRUE);
      //
      task->SetDefaultTrackCuts(kFALSE);
      task->SetFillArmPodTree(kFALSE);
      task->SetFillOnlyHists(kTRUE);
      task->SetRunOnGrid(kTRUE);
      task->fEventCuts.fUseVariablesCorrelationCuts = true;
    }
    break;
    case 5:{
      std::cout << settingType << " Info::marsland: Marians event tree " << std::endl;
      task->SetDefaultTrackCuts(kFALSE);
      task->SetRunOnGrid(kTRUE);
      task->SetFillEventInfo(kTRUE);
      task->SetUseCouts(kTRUE);
    }
    break;
    //
    // ====================================================================================
    // =========================== MC Closure on Lego train  ==============================
    // ====================================================================================
    //
    case 50:{
      std::cout << settingType << " Info::marsland: MC full on lego train --> eff matrix is not filled " << std::endl;
      task->SetEffMatrix(kTRUE);  task->SetIsMCtrue(kTRUE);  task->SetFillAllCutVariables(kTRUE);  // conditions to enter FillMCFull_NetParticles()
      //
      task->SetFillTreeMC(kTRUE);
      task->SetFillGenDistributions(kFALSE);
      task->SetTrackOriginType(0);   // 0:prim, 1: prim+weak, 2:prim+material, 3:prim+material+weak, 4:full scan
      task->SetUsePtCut(2); // 0: tpc momcut, 1: vertex momcut, 2: pT cut
      // acceptance
      const Int_t tmpEtaBinsMC = 8;
      const Int_t tmpMomBinsMC = 2;
      Float_t tmpetaDownArr[tmpEtaBinsMC] = {-0.1,-0.2,-0.3,-0.4,-0.5,-0.6,-0.7,-0.8};
      Float_t tmpetaUpArr[tmpEtaBinsMC]   = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};
      Float_t tmppDownArr[tmpMomBinsMC] = { 0.4, 0.6};
      Float_t tmppUpArr[tmpMomBinsMC]   = { 1.0, 1.5};
      task->SetMCEtaScanArray(tmpEtaBinsMC, tmpetaDownArr, tmpetaUpArr);
      task->SetMCMomScanArray(tmpMomBinsMC, tmppDownArr,   tmppUpArr);
      // resonances to exclude
      const Int_t tmpNresonances = 1;
      TString tmpResArr[tmpNresonances] = {"xxx"};
      task->SetMCResonanceArray(tmpNresonances,tmpResArr);
    }
    break;
    case 51:{
      std::cout << settingType << " Info::marsland: MC full on lego train --> eff matrix is not filled " << std::endl;
      task->SetEffMatrix(kTRUE);  task->SetIsMCtrue(kTRUE);  task->SetFillAllCutVariables(kTRUE);  // conditions to enter FillMCFull_NetParticles()
      //
      task->SetFillTreeMC(kTRUE);
      task->SetFillGenDistributions(kFALSE);
      task->SetTrackOriginType(0);   // 0:prim, 1: prim+weak, 2:prim+material, 3:prim+material+weak, 4:full scan
      task->SetUsePtCut(2); // 0: tpc momcut, 1: vertex momcut, 2: pT cut
      // acceptance
      const Int_t tmpEtaBinsMC = 8;
      const Int_t tmpMomBinsMC = 1;
      Float_t tmpetaDownArr[tmpEtaBinsMC] = {-0.1,-0.2,-0.3,-0.4,-0.5,-0.6,-0.7,-0.8};
      Float_t tmpetaUpArr[tmpEtaBinsMC]   = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};
      Float_t tmppDownArr[tmpMomBinsMC] = { 0.6};
      Float_t tmppUpArr[tmpMomBinsMC]   = { 1.5};
      task->SetMCEtaScanArray(tmpEtaBinsMC, tmpetaDownArr, tmpetaUpArr);
      task->SetMCMomScanArray(tmpMomBinsMC, tmppDownArr,   tmppUpArr);
      // resonances to exclude
      const Int_t tmpNresonances = 1;
      TString tmpResArr[tmpNresonances] = {"xxx"};
      task->SetMCResonanceArray(tmpNresonances,tmpResArr);
    }
    break;
    case 52:{
      std::cout << settingType << " Info::marsland: MC full on lego train --> eff matrix is not filled " << std::endl;
      task->SetEffMatrix(kTRUE);  task->SetIsMCtrue(kTRUE);  task->SetFillAllCutVariables(kTRUE);  // conditions to enter FillMCFull_NetParticles()
      //
      task->SetFillTreeMC(kTRUE);
      task->SetFillGenDistributions(kFALSE);
      task->SetTrackOriginType(0);   // 0:prim, 1: prim+weak, 2:prim+material, 3:prim+material+weak, 4:full scan
      task->SetUsePtCut(2); // 0: tpc momcut, 1: vertex momcut, 2: pT cut
      //
      // acceptance
      const Int_t tmpEtaBinsMC = 8;
      const Int_t tmpMomBinsMC = 1;
      Float_t tmpetaDownArr[tmpEtaBinsMC] = {-0.1,-0.2,-0.3,-0.4,-0.5,-0.6,-0.7,-0.8};
      Float_t tmpetaUpArr[tmpEtaBinsMC]   = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};
      Float_t tmppDownArr[tmpMomBinsMC] = { 0.6};
      Float_t tmppUpArr[tmpMomBinsMC]   = { 2.0};
      task->SetMCEtaScanArray(tmpEtaBinsMC, tmpetaDownArr, tmpetaUpArr);
      task->SetMCMomScanArray(tmpMomBinsMC, tmppDownArr,   tmppUpArr);
      // resonances to exclude
      const Int_t tmpNresonances = 1;
      TString tmpResArr[tmpNresonances] = {"xxx"};
      task->SetMCResonanceArray(tmpNresonances,tmpResArr);
    }
    break;
    case 53:{
      std::cout << settingType << " Info::marsland: MC full on lego train --> eff matrix is not filled " << std::endl;
      task->SetEffMatrix(kTRUE);  task->SetIsMCtrue(kTRUE);  task->SetFillAllCutVariables(kTRUE);  // conditions to enter FillMCFull_NetParticles()
      //
      task->SetFillTreeMC(kTRUE);
      task->SetFillGenDistributions(kFALSE);
      task->SetTrackOriginType(0);   // 0:prim, 1: prim+weak, 2:prim+material, 3:prim+material+weak, 4:full scan
      task->SetUsePtCut(2); // 0: tpc momcut, 1: vertex momcut, 2: pT cut
      // acceptance
      const Int_t tmpEtaBinsMC = 8;
      const Int_t tmpMomBinsMC = 1;
      Float_t tmpetaDownArr[tmpEtaBinsMC] = {-0.1,-0.2,-0.3,-0.4,-0.5,-0.6,-0.7,-0.8};
      Float_t tmpetaUpArr[tmpEtaBinsMC]   = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};
      Float_t tmppDownArr[tmpMomBinsMC] = { 0.8};
      Float_t tmppUpArr[tmpMomBinsMC]   = { 2.0};
      task->SetMCEtaScanArray(tmpEtaBinsMC, tmpetaDownArr, tmpetaUpArr);
      task->SetMCMomScanArray(tmpMomBinsMC, tmppDownArr,   tmppUpArr);
      // resonances to exclude
      const Int_t tmpNresonances = 1;
      TString tmpResArr[tmpNresonances] = {"xxx"};
      task->SetMCResonanceArray(tmpNresonances,tmpResArr);
    }
    break;
    //
    // ====================================================================================
    // ======================== MC to run at GSI or runGrid.C  ============================
    // ====================================================================================
    //
    case 60:{
      std::cout << settingType << " Info::marsland: fTreeMC + mcFull + EffMatrix + + dist + full cout  --> At GSI or with RunGrid.C " << std::endl;
      task->SetEffMatrix(kTRUE);  task->SetIsMCtrue(kTRUE);  task->SetFillAllCutVariables(kTRUE);  // conditions to enter FillMCFull_NetParticles()
      //
      task->SetFillGenDistributions(kTRUE);
      task->SetUseCouts(kTRUE);
      task->SetFillTreeMC(kTRUE);
      task->SetTrackOriginType(4);   // 0:prim, 1: prim+weak, 2:prim+material, 3:prim+material+weak, 4:full scan
      task->SetUsePtCut(1); // 0: tpc momcut, 1: vertex momcut, 2: pT cut
      // acceptance
      const Int_t tmpEtaBinsMC = 8;
      const Int_t tmpMomBinsMC = 4;
      Float_t tmpetaDownArr[tmpEtaBinsMC] = {-0.1,-0.2,-0.3,-0.4,-0.5,-0.6,-0.7,-0.8};
      Float_t tmpetaUpArr[tmpEtaBinsMC]   = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};
      Float_t tmppDownArr[tmpMomBinsMC] = { 0.4, 0.6, 0.6, 0.8};
      Float_t tmppUpArr[tmpMomBinsMC]   = { 1.0, 1.5, 2.0, 2.0};
      task->SetMCEtaScanArray(tmpEtaBinsMC, tmpetaDownArr, tmpetaUpArr);
      task->SetMCMomScanArray(tmpMomBinsMC, tmppDownArr,   tmppUpArr);
      // resonances to exclude
      const Int_t tmpNresonances = 1;
      TString tmpResArr[tmpNresonances] = {"xxx"};
      task->SetMCResonanceArray(tmpNresonances,tmpResArr);
    }
    break;
    case 61:{
      std::cout << settingType << " Info::marsland: mcFull + EffMatrix " << std::endl;
      task->SetEffMatrix(kTRUE);  task->SetIsMCtrue(kTRUE);  task->SetFillAllCutVariables(kTRUE);  // conditions to enter FillMCFull_NetParticles()
      //
      task->SetUseCouts(kFALSE);
      task->SetTrackOriginType(4);   // 0:prim, 1: prim+weak, 2:prim+material, 3:prim+material+weak, 4:full scan
      task->SetUsePtCut(1); // 0: tpc momcut, 1: vertex momcut, 2: pT cut
      // acceptance
      const Int_t tmpEtaBinsMC = 1;
      const Int_t tmpMomBinsMC = 1;
      Float_t tmpetaDownArr[tmpEtaBinsMC] = {-0.8};
      Float_t tmpetaUpArr[tmpEtaBinsMC]   = { 0.8};
      Float_t tmppDownArr[tmpMomBinsMC] = { 0.6};
      Float_t tmppUpArr[tmpMomBinsMC]   = { 1.5};
      task->SetMCEtaScanArray(tmpEtaBinsMC, tmpetaDownArr, tmpetaUpArr);
      task->SetMCMomScanArray(tmpMomBinsMC, tmppDownArr,   tmppUpArr);
      // resonances to exclude
      const Int_t tmpNresonances = 1;
      TString tmpResArr[tmpNresonances] = {"xxx"};
      task->SetMCResonanceArray(tmpNresonances,tmpResArr);
    }
    break;
    case 62:{
      std::cout << settingType << " Info::marsland: MC full on lego train --> full acceptance " << std::endl;
      task->SetEffMatrix(kTRUE);  task->SetIsMCtrue(kTRUE);  task->SetFillAllCutVariables(kTRUE);  // conditions to enter FillMCFull_NetParticles()
      //
      // task->SetRunOnGrid(kTRUE); // do not fill eff matrix
      task->SetUseCouts(kTRUE);
      task->SetRunFastSimulation(kTRUE);
      task->SetFillTreeMC(kFALSE);
      task->SetFillGenDistributions(kFALSE);
      task->SetDefaultTrackCuts(kTRUE);
      task->SetFillNudynFastGen(kTRUE);
      task->SetTrackOriginType(0);   // 0:prim, 1: prim+weak, 2:prim+material, 3:prim+material+weak, 4:full scan
      task->SetRapidityType(1);      // 0:pseudorapidity, 1: rapidity
      task->SetUsePtCut(1); // 0: tpc momcut, 1: vertex momcut, 2: pT cut
      //
      // acceptance
      const Int_t tmpEtaBinsMC = 24;
      const Int_t tmpMomBinsMC = 3;
      Float_t tmpetaDownArr[tmpEtaBinsMC] = {-0.5,-1.,-1.5,-2.,-2.5,-3.,-3.5,-4.,-4.5,-5.,-5.5,-6.,-6.5,-7.,-7.5,-8,-8.5,-9,-9.5,-10.,-10.5,-11.,-11.5,-12.};
      Float_t tmpetaUpArr[tmpEtaBinsMC]   = { 0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5., 5.5, 6., 6.5, 7., 7.5, 8, 8.5, 9, 9.5, 10., 10.5, 11., 11.5, 12.};
      Float_t tmppDownArr[tmpMomBinsMC] = { 0.2, 0.6, 0.    };
      Float_t tmppUpArr[tmpMomBinsMC]   = { 1.5, 1.5, 10000.};
      task->SetMCEtaScanArray(tmpEtaBinsMC, tmpetaDownArr, tmpetaUpArr);
      task->SetMCMomScanArray(tmpMomBinsMC, tmppDownArr,   tmppUpArr);
      // resonances to exclude
      const Int_t tmpNresonances = 1;
      TString tmpResArr[tmpNresonances] = {"xxx"};
      task->SetMCResonanceArray(tmpNresonances,tmpResArr);
      //
      // baryons to be included for netbaryon analysis
      const Int_t tmpNbaryons = 7;
      Int_t tmpBaryonArr[tmpNbaryons] = {2212,2112,2224,2214,2114,1114,3122};  // {p,n,delta++,delta+,delta0,delta-,Lambda,}
      task->SetMCBaryonArray(tmpNbaryons,tmpBaryonArr);
      //
      // // baryons to be included for netbaryon analysis
      // const Int_t tmpNbaryons = 2;
      // Int_t tmpBaryonArr[tmpNbaryons] = {2212,2112};  // {p,n}
      // task->SetMCBaryonArray(tmpNbaryons,tmpBaryonArr);
    }
    break;
    case 63:{
      std::cout << settingType << " Info::marsland: fTreeMC and mcFull for the net-particle study " << std::endl;
      task->SetEffMatrix(kTRUE);  task->SetIsMCtrue(kTRUE);  task->SetFillAllCutVariables(kTRUE);  // conditions to enter FillMCFull_NetParticles()
      //
      // task->SetRunOnGrid(kTRUE); // do not fill eff matrix
      task->SetUseCouts(kTRUE);
      task->SetFillTreeMC(kTRUE);         // fills fTreeMC
      task->SetFillArmPodTree(kTRUE);         // fills fTreeMC
      task->SetDefaultEventCuts(kFALSE);
      task->SetDefaultTrackCuts(kTRUE);
      task->SetTrackOriginType(0);   // 0:prim, 1: prim+weak, 2:prim+material, 3:prim+material+weak, 4:full scan
      task->SetRapidityType(0);      // 0:pseudorapidity, 1: rapidity
      task->SetUsePtCut(1);          // 0: tpc momcut, 1: vertex momcut, 2: pT cut
      //
      // acceptance
      const Int_t tmpEtaBinsMC = 8;
      const Int_t tmpMomBinsMC = 4;
      Float_t tmpetaDownArr[tmpEtaBinsMC] = {-0.1,-0.2,-0.3,-0.4,-0.5,-0.6,-0.7,-0.8};
      Float_t tmpetaUpArr[tmpEtaBinsMC]   = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};
      Float_t tmppDownArr[tmpMomBinsMC] = { 0.4, 0.6, 0.6, 0.8};
      Float_t tmppUpArr[tmpMomBinsMC]   = { 1.0, 1.5, 2.0, 2.0};
      task->SetMCEtaScanArray(tmpEtaBinsMC, tmpetaDownArr, tmpetaUpArr);
      task->SetMCMomScanArray(tmpMomBinsMC, tmppDownArr,   tmppUpArr);
      // resonances to exclude
      const Int_t tmpNresonances = 1;
      TString tmpResArr[tmpNresonances] = {"xxx"};
      task->SetMCResonanceArray(tmpNresonances,tmpResArr);
      //
      // baryons to be included for netbaryon analysis
      const Int_t tmpNbaryons = 7;
      Int_t tmpBaryonArr[tmpNbaryons] = {2212,2112,2224,2214,2114,1114,3122};  // {p,n,delta++,delta+,delta0,delta-,Lambda,}
      task->SetMCBaryonArray(tmpNbaryons,tmpBaryonArr);
    }
    break;
    //
    // ====================================================================================
    // ============================== MC to run Local  ====================================
    // ====================================================================================
    //
    case 70:{
      std::cout << settingType << " Info::marsland: too see distributions " << std::endl;
      task->SetEffMatrix(kTRUE);  task->SetIsMCtrue(kTRUE);  task->SetFillAllCutVariables(kTRUE);  // conditions to enter FillMCFull_NetParticles()
      //
      task->SetFillGenDistributions(kTRUE);
      task->SetRunOnGrid(kTRUE); // do not fill eff matrix
      task->SetTrackOriginType(0);   // 0:prim, 1: prim+weak, 2:prim+material, 3:prim+material+weak, 4:full scan
      task->SetUsePtCut(1); // 0: tpc momcut, 1: vertex momcut, 2: pT cut
      //
      task->SetMomLowerEdge(0.);
      task->SetMomUpperEdge(10000.);
      task->SetEtaLowerEdge(-10000.);
      task->SetEtaUpperEdge( 10000.);
      // acceptance
      const Int_t tmpEtaBinsMC = 1;
      const Int_t tmpMomBinsMC = 1;
      Float_t tmpetaDownArr[tmpEtaBinsMC] = {-0.8};
      Float_t tmpetaUpArr[tmpEtaBinsMC]   = { 0.8};
      Float_t tmppDownArr[tmpMomBinsMC] = { 0.6};
      Float_t tmppUpArr[tmpMomBinsMC]   = { 1.5};
      task->SetMCEtaScanArray(tmpEtaBinsMC, tmpetaDownArr, tmpetaUpArr);
      task->SetMCMomScanArray(tmpMomBinsMC, tmppDownArr,   tmppUpArr);
      // resonances to exclude
      const Int_t tmpNresonances = 1;
      TString tmpResArr[tmpNresonances] = {"xxx"};
      task->SetMCResonanceArray(tmpNresonances,tmpResArr);
    }
    break;
    //
    // ====================================================================================
    // ============================== MC Eff Matrix =======================================
    // ====================================================================================
    //
    case 80:{
      std::cout << settingType << " Info::marsland: fill only eff matrix " << std::endl;
      task->SetEffMatrix(kTRUE);  task->SetIsMCtrue(kTRUE);  task->SetFillAllCutVariables(kTRUE);  // conditions to enter FillMCFull_NetParticles()
      //
      task->SetTrackOriginType(0);   // 0:prim, 1: prim+weak, 2:prim+material, 3:prim+material+weak, 4:full scan
      task->SetUsePtCut(1); // 0: tpc momcut, 1: vertex momcut, 2: pT cut
      //
      // acceptance
      const Int_t tmpEtaBinsMC = 1;
      const Int_t tmpMomBinsMC = 1;
      Float_t tmpetaDownArr[tmpEtaBinsMC] = {-0.8};
      Float_t tmpetaUpArr[tmpEtaBinsMC]   = { 0.8};
      Float_t tmppDownArr[tmpMomBinsMC] = { 0.6};
      Float_t tmppUpArr[tmpMomBinsMC]   = { 1.5};
      task->SetMCEtaScanArray(tmpEtaBinsMC, tmpetaDownArr, tmpetaUpArr);
      task->SetMCMomScanArray(tmpMomBinsMC, tmppDownArr,   tmppUpArr);
      // resonances to exclude
      const Int_t tmpNresonances = 1;
      TString tmpResArr[tmpNresonances] = {"xxx"};
      task->SetMCResonanceArray(tmpNresonances,tmpResArr);
    }
    break;
    //
    // ====================================================================================
    // =================================== MC FastGen  ====================================
    // ====================================================================================
    //
    case 100:{
      std::cout << settingType << " Info::marsland: FastSimul: StandardTPCITScuts 8EtaBin_150pBins_9centBins (REFERENCE settings) " << std::endl;
      task->SetIsMCtrue(kTRUE);
      task->SetNEtabins(8);
      task->SetRunFastSimulation(kTRUE);
    }
    break;
    case 101:{
      std::cout << settingType << " Info::marsland: FastSimul: StandardTPCITScuts 8EtaBin_150pBins_9centBins ETA DEPENDENCE " << std::endl;
      task->SetIsMCtrue(kTRUE);
      task->SetNEtabins(8);
      task->SetRunFastSimulation(kTRUE);
      const Int_t tmpCentbins  = 10;
      const Int_t tmpEtaBinsMC = 8;
      const Int_t tmpMomBinsMC = 2;
      Float_t tmpfxCentBins[tmpCentbins] = {0,5,10,20,30,40,50,60,70,80};
      Float_t tmpetaDownArr[tmpEtaBinsMC] = {-0.1,-0.2,-0.3,-0.4,-0.5,-0.6,-0.7,-0.8};
      Float_t tmpetaUpArr[tmpEtaBinsMC]   = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};
      Float_t tmppDownArr[tmpMomBinsMC] = { 0.2, 0.6};
      Float_t tmppUpArr[tmpMomBinsMC]   = { 1.5, 1.5};
      task->SetMCEtaScanArray(tmpEtaBinsMC, tmpetaDownArr, tmpetaUpArr);
      task->SetMCMomScanArray(tmpMomBinsMC, tmppDownArr,   tmppUpArr);
      task->SetCentralityBinning(tmpCentbins,tmpfxCentBins);
    }
    break;
    case 102:{
      std::cout << settingType << " Info::marsland: FastSimul: StandardTPCITScuts 8EtaBin_150pBins_9centBins Momentum DEPENDENCE " << std::endl;
      task->SetIsMCtrue(kTRUE);
      task->SetNEtabins(8);
      task->SetRunFastSimulation(kTRUE);
      const Int_t tmpCentbins  = 10;
      const Int_t tmpEtaBinsMC = 1;
      const Int_t tmpMomBinsMC = 9;
      Float_t tmpfxCentBins[tmpCentbins] = {0,5,10,20,30,40,50,60,70,80};
      Float_t tmpetaDownArr[tmpEtaBinsMC] = {-0.8};
      Float_t tmpetaUpArr[tmpEtaBinsMC]   = { 0.8};
      Float_t tmppDownArr[tmpMomBinsMC] = { 0.2, 0.2, 0.2, 0.3, 0.3, 0.3, 0.5, 0.5, 0.5};
      Float_t tmppUpArr[tmpMomBinsMC]   = { 1.5, 1.8, 2.0, 1.5, 1.8, 2.0, 1.5, 1.8, 2.0};
      task->SetMCEtaScanArray(tmpEtaBinsMC, tmpetaDownArr, tmpetaUpArr);
      task->SetMCMomScanArray(tmpMomBinsMC, tmppDownArr,   tmppUpArr);
      task->SetCentralityBinning(tmpCentbins,tmpfxCentBins);
    }
    break;
    case 103:{
      std::cout << settingType << " Info::marsland: FullSimul: StandardTPCITScuts 8EtaBin_150pBins_9centBins EffMatrix " << std::endl;
      task->SetEffMatrix(kTRUE);
      task->SetIsMCtrue(kTRUE);
      task->SetNEtabins(8);
      task->SetRunFastSimulation(kFALSE);
      task->SetFillArmPodTree(kFALSE);
      task->SetDeDxCheck(kFALSE);
      const Int_t tmpCentbins  = 10;
      const Int_t tmpEtaBinsMC = 10;
      const Int_t tmpMomBinsMC = 2;
      const Int_t tmpNresonances = 2;
      Float_t tmpfxCentBins[tmpCentbins] = {0,5,10,20,30,40,50,60,70,80};
      Float_t tmpetaDownArr[tmpEtaBinsMC] = {-0.1,-0.2,-0.3,-0.4,-0.5,-0.6,-0.7,-0.8,-0.9,-1.};
      Float_t tmpetaUpArr[tmpEtaBinsMC]   = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.};
      Float_t tmppDownArr[tmpMomBinsMC] = { 0.2, 0.6};
      Float_t tmppUpArr[tmpMomBinsMC]   = { 1.5, 1.5};
      //       TString tmpResArr[tmpNresonances] = {"rho","phi","eta","omega","Delta","Lambda"};
      //       TString tmpResArr[tmpNresonances] = {"rho","phi","Delta","omega"};
      TString tmpResArr[tmpNresonances] = {"rho","phi"};
      task->SetMCResonanceArray(tmpNresonances,tmpResArr);
      task->SetMCEtaScanArray(tmpEtaBinsMC, tmpetaDownArr, tmpetaUpArr);
      task->SetMCMomScanArray(tmpMomBinsMC, tmppDownArr,   tmppUpArr);
      task->SetCentralityBinning(tmpCentbins,tmpfxCentBins);
    }
    break;
    case 104:{
      std::cout << settingType << " Info::marsland: FullSimul: StandardTPCITScuts 8EtaBin_150pBins_9centBins EffMatrix " << std::endl;
      task->SetEffMatrix(kTRUE);
      task->SetIsMCtrue(kTRUE);
      task->SetNEtabins(8);
      task->SetRunFastSimulation(kFALSE);
      task->SetFillArmPodTree(kFALSE);
      task->SetDeDxCheck(kFALSE);
      const Int_t tmpCentbins  = 10;
      const Int_t tmpEtaBinsMC = 1;
      const Int_t tmpMomBinsMC = 1;
      Float_t tmpfxCentBins[tmpCentbins] = {0,5,10,20,30,40,50,60,70,80};
      Float_t tmpetaDownArr[tmpEtaBinsMC] = {-0.8};
      Float_t tmpetaUpArr[tmpEtaBinsMC]   = { 0.8};
      Float_t tmppDownArr[tmpMomBinsMC] = { 0.2};
      Float_t tmppUpArr[tmpMomBinsMC]   = { 3.2};
      task->SetMCEtaScanArray(tmpEtaBinsMC, tmpetaDownArr, tmpetaUpArr);
      task->SetMCMomScanArray(tmpMomBinsMC, tmppDownArr,   tmppUpArr);
      task->SetCentralityBinning(tmpCentbins,tmpfxCentBins);
    }
    break;
    case 105:{
      std::cout << settingType << " Info::marsland: FastSimul: StandardTPCITScuts 8EtaBin_150pBins_9centBins (REFERENCE settings): Fill only DnchDeta " << std::endl;
      task->SetIsMCtrue(kTRUE);
      task->SetNEtabins(8);
      task->SetRunFastSimulation(kTRUE);
      task->SetFillDnchDeta(kTRUE);
    }
    break;
    case 106:{
      std::cout << settingType << " Info::marsland: THnSparse: StandardTPCITScuts 8EtaBin_150pBins_9centBins (REFERENCE settings) MC CLOSURE with TOF cut " << std::endl;
      task->SetIncludeTOF(kTRUE);
      task->SetIsMCtrue(kTRUE);
      task->SetEffMatrix(kTRUE);
      task->SetNEtabins(8);
      const Int_t tmpCentbins  = 10;
      const Int_t tmpEtaBinsMC = 2;
      const Int_t tmpMomBinsMC = 4;
      Float_t tmpfxCentBins[tmpCentbins] = {0,5,10,20,30,40,50,60,70,80};
      Float_t tmpetaDownArr[tmpEtaBinsMC] = {-0.5,-0.8};
      Float_t tmpetaUpArr[tmpEtaBinsMC]   = { 0.5, 0.8};
      Float_t tmppDownArr[tmpMomBinsMC] = { 0.2, 0.2, 0.6, 0.6};
      Float_t tmppUpArr[tmpMomBinsMC]   = { 1.5, 1.8, 1.5, 1.8};
      task->SetMCEtaScanArray(tmpEtaBinsMC, tmpetaDownArr, tmpetaUpArr);
      task->SetMCMomScanArray(tmpMomBinsMC, tmppDownArr,   tmppUpArr);
      task->SetCentralityBinning(tmpCentbins,tmpfxCentBins);
    }
    break;
    case 107:{
      std::cout << settingType << " Info::marsland: THnSparse: StandardTPCITScuts 8EtaBin_150pBins_9centBins (REFERENCE settings) MC CLOSURE with NO ITS cut " << std::endl;
      task->SetIncludeTOF(kFALSE);
      task->SetIsMCtrue(kTRUE);
      task->SetEffMatrix(kTRUE);
      task->SetNEtabins(8);
      const Int_t tmpCentbins  = 10;
      const Int_t tmpEtaBinsMC = 2;
      const Int_t tmpMomBinsMC = 4;
      Float_t tmpfxCentBins[tmpCentbins] = {0,5,10,20,30,40,50,60,70,80};
      Float_t tmpetaDownArr[tmpEtaBinsMC] = {-0.5,-0.8};
      Float_t tmpetaUpArr[tmpEtaBinsMC]   = { 0.5, 0.8};
      Float_t tmppDownArr[tmpMomBinsMC] = { 0.2, 0.2, 0.6, 0.6};
      Float_t tmppUpArr[tmpMomBinsMC]   = { 1.5, 1.8, 1.5, 1.8};
      task->SetMCEtaScanArray(tmpEtaBinsMC, tmpetaDownArr, tmpetaUpArr);
      task->SetMCMomScanArray(tmpMomBinsMC, tmppDownArr,   tmppUpArr);
      task->SetCentralityBinning(tmpCentbins,tmpfxCentBins);
    }
    break;
    case 108:{
      std::cout << settingType << " Info::marsland: FULL MC at GSI --> eta and mom scan + removal of resonances + EffMatrix [0.2,0.6]<p<[1.5,3.2] GeV/c  " << std::endl;
      // FULL MC settings
      task->SetIsMCtrue(kTRUE);
      task->SetFillAllCutVariables(kTRUE);
      task->SetEffMatrix(kTRUE);
      task->SetFillArmPodTree(kTRUE);
      //
      task->SetUseCouts(kTRUE);
      task->SetNEtabins(10);
      task->SetEtaLowerEdge(-1.);
      task->SetEtaUpperEdge(1.);
      // eta bin scan
      const Int_t tmpEtaBinsMC = 2;
      Float_t tmpetaDownArr[tmpEtaBinsMC] = {-0.8,-1.0};
      Float_t tmpetaUpArr[tmpEtaBinsMC]   = { 0.8, 1.0};
      task->SetMCEtaScanArray(tmpEtaBinsMC, tmpetaDownArr, tmpetaUpArr);
      // mom bin scan
      const Int_t tmpMomBinsMC = 2;
      Float_t tmppDownArr[tmpMomBinsMC] = { 0.6, 0.2};
      Float_t tmppUpArr[tmpMomBinsMC]   = { 1.5, 1.5};
      task->SetMCMomScanArray(tmpMomBinsMC, tmppDownArr,   tmppUpArr);
      // cent bins
      const Int_t tmpCentbins  = 10;
      Float_t tmpfxCentBins[tmpCentbins] = {0,5,10,20,30,40,50,60,70,80};
      task->SetCentralityBinning(tmpCentbins,tmpfxCentBins);
      // resonances to exclude
      const Int_t tmpNresonances = 1;
      TString tmpResArr[tmpNresonances] = {"xxx"};
      task->SetMCResonanceArray(tmpNresonances,tmpResArr);
      //       const Int_t tmpNresonances = 5;
      //       TString tmpResArr[tmpNresonances] = {"rho","phi","omega","eta","Delta"};
      //       task->SetMCResonanceArray(tmpNresonances,tmpResArr);
    }
    break;
    case 109:{
      std::cout << settingType << " Info::marsland: Higher Moments MC closure " << std::endl;
      task->SetRunOnGrid(getFromAlien);
      std::cout << settingType << " Info::marsland: runOnGrid = " << getFromAlien << std::endl;
      task->SetUseCouts(kTRUE);
      task->SetIsMCtrue(kTRUE);
      task->SetFillHigherMomentsMCclosure(kTRUE);
      // eta cent and mom scan
      const Int_t tmpEtaBinsMC = 8;
      const Int_t tmpMomBinsMC = 3;
      const Int_t tmpCentbins  = 10;
      Float_t tmpfxCentBins[tmpCentbins]  = {0,5,10,20,30,40,50,60,70,80};
      Float_t tmpetaDownArr[tmpEtaBinsMC] = {-0.1,-0.2,-0.3,-0.4,-0.5,-0.6,-0.7,-0.8};
      Float_t tmpetaUpArr[tmpEtaBinsMC]   = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};
      Float_t tmppDownArr[tmpMomBinsMC]   = { 0.4, 0.2, 0.6};
      Float_t tmppUpArr[tmpMomBinsMC]     = { 1.0, 1.5, 1.5};
      task->SetMCEtaScanArray(tmpEtaBinsMC, tmpetaDownArr, tmpetaUpArr);
      task->SetMCMomScanArray(tmpMomBinsMC, tmppDownArr,   tmppUpArr);
      task->SetCentralityBinning(tmpCentbins,tmpfxCentBins);
      // Read First Moments for higher moments
      task->SetLookUpTableEfficiencyCorrection(lookUpTree, 0, tmppDownArr, tmppUpArr, tmpfxCentBins, tmpetaUpArr, tmpMomBinsMC, tmpCentbins, tmpEtaBinsMC);
      task->SetLookUpTableEfficiencyCorrection(lookUpTree, 1, tmppDownArr, tmppUpArr, tmpfxCentBins, tmpetaUpArr, tmpMomBinsMC, tmpCentbins, tmpEtaBinsMC);
      task->SetLookUpTableEfficiencyCorrection(lookUpTree, 2, tmppDownArr, tmppUpArr, tmpfxCentBins, tmpetaUpArr, tmpMomBinsMC, tmpCentbins, tmpEtaBinsMC);
      // void SetLookUpTableEfficiencyCorrection(lookUpTree,partType,pDownArr[],pUpArr[],centArr[],etaArr[],tmpMomBinsMC,tmpCentbins,tmpEtaBinsMC)
    }
    break;
    //
    // ====================================================================================
    // ================================ FastGen  Settings =================================
    // ====================================================================================
    //
    case 200:{
      std::cout << settingType << " Info::marsland: Fast Gen on LEGO --> Calculate higher moments using look-up table [0.2,0.6]<p<1.5 GeV/c " << std::endl;
      task->SetRunFastHighMomentCal(kTRUE);
      task->SetRunOnGrid(kTRUE);
      task->SetUseCouts(kFALSE);
      task->SetIsMCtrue(kTRUE);
      // eta bin scan
      const Int_t tmpEtaBinsMC = 20;
      Float_t tmpetaDownArr[tmpEtaBinsMC] = {-0.1,-0.2,-0.3,-0.4,-0.5,-0.6,-0.7,-0.8,-0.9,-1.,-1.1,-1.2,-1.3,-1.4,-1.5,-1.6,-1.7,-1.8,-1.9,-2};
      Float_t tmpetaUpArr[tmpEtaBinsMC]   = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2};
      task->SetMCEtaScanArray(tmpEtaBinsMC, tmpetaDownArr, tmpetaUpArr);
      // mom bin scan
      const Int_t tmpMomBinsMC = 2;
      Float_t tmppDownArr[tmpMomBinsMC] = { 0.2, 0.6};
      Float_t tmppUpArr[tmpMomBinsMC]   = { 1.5, 1.5};
      task->SetMCMomScanArray(tmpMomBinsMC, tmppDownArr,   tmppUpArr);
      // cent bins
      const Int_t tmpCentbins  = 10;
      Float_t tmpfxCentBins[tmpCentbins] = {0,5,10,20,30,40,50,60,70,80};
      task->SetCentralityBinning(tmpCentbins,tmpfxCentBins);
      // resonances to exclude
      const Int_t tmpNresonances = 1;
      TString tmpResArr[tmpNresonances] = {"xxx"};
      task->SetMCResonanceArray(tmpNresonances,tmpResArr);
      // Read First Moments
      task->SetLookUpTableFirstMoments(lookUpTree,0, tmppDownArr,tmpfxCentBins,tmpetaUpArr,tmpMomBinsMC,tmpCentbins,tmpEtaBinsMC);
      task->SetLookUpTableFirstMoments(lookUpTree,1, tmppDownArr,tmpfxCentBins,tmpetaUpArr,tmpMomBinsMC,tmpCentbins,tmpEtaBinsMC);
      task->SetLookUpTableFirstMoments(lookUpTree,2, tmppDownArr,tmpfxCentBins,tmpetaUpArr,tmpMomBinsMC,tmpCentbins,tmpEtaBinsMC);
      task->SetLookUpTableFirstMoments(lookUpTree,9, tmppDownArr,tmpfxCentBins,tmpetaUpArr,tmpMomBinsMC,tmpCentbins,tmpEtaBinsMC);
      task->SetLookUpTableFirstMoments(lookUpTree,11,tmppDownArr,tmpfxCentBins,tmpetaUpArr,tmpMomBinsMC,tmpCentbins,tmpEtaBinsMC);
    }
    break;
    case 201:{
      std::cout << settingType << " Info::marsland: Fast Gen at GSI --> eta and mom scan + removal of resonances [0.2,0.6]<p<1.5 GeV/c " << std::endl;
      task->SetRunFastSimulation(kTRUE);
      task->SetRunOnGrid(kTRUE);
      task->SetPercentageOfEvents(0);
      task->SetUseCouts(kFALSE);
      task->SetIsMCtrue(kTRUE);
      // eta bin scan
      const Int_t tmpEtaBinsMC = 20;
      Float_t tmpetaDownArr[tmpEtaBinsMC] = {-0.1,-0.2,-0.3,-0.4,-0.5,-0.6,-0.7,-0.8,-0.9,-1.,-1.1,-1.2,-1.3,-1.4,-1.5,-1.6,-1.7,-1.8,-1.9,-2};
      Float_t tmpetaUpArr[tmpEtaBinsMC]   = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2};
      task->SetMCEtaScanArray(tmpEtaBinsMC, tmpetaDownArr, tmpetaUpArr);
      // mom bin scan
      const Int_t tmpMomBinsMC = 2;
      Float_t tmppDownArr[tmpMomBinsMC] = { 0.2,0.6};
      Float_t tmppUpArr[tmpMomBinsMC]   = { 1.5,1.5};
      task->SetMCMomScanArray(tmpMomBinsMC, tmppDownArr,   tmppUpArr);
      // cent bins
      const Int_t tmpCentbins  = 10;
      Float_t tmpfxCentBins[tmpCentbins] = {0,5,10,20,30,40,50,60,70,80};
      task->SetCentralityBinning(tmpCentbins,tmpfxCentBins);
      // resonances to exclude
      const Int_t tmpNresonances = 1;
      TString tmpResArr[tmpNresonances] = {"xxx"};
      task->SetMCResonanceArray(tmpNresonances,tmpResArr);
      //       const Int_t tmpNresonances = 5;
      //       TString tmpResArr[tmpNresonances] = {"rho","phi","omega","eta","Delta"};
      //       task->SetMCResonanceArray(tmpNresonances,tmpResArr);
    }
    break;
    case 202:{
      std::cout << settingType << " Info::marsland: Fast Gen at GSI --> eta and mom scan + removal of resonances in full acceptance [0.2,0.6]<p<1.5 GeV/c  " << std::endl;
      task->SetRunFastSimulation(kTRUE);
      task->SetRunOnGrid(kTRUE);
      task->SetPercentageOfEvents(0);
      task->SetIsMCtrue(kTRUE);
      // eta bin scan
      const Int_t tmpEtaBinsMC = 20;
      Float_t tmpetaDownArr[tmpEtaBinsMC] = {-0.5,-1.,-1.5,-2.,-2.5,-3.,-3.5,-4.,-4.5,-5.,-5.5,-6.,-6.5,-7.,-7.5,-8,-8.5,-9,-9.5,-10.};
      Float_t tmpetaUpArr[tmpEtaBinsMC]   = { 0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5., 5.5, 6., 6.5, 7., 7.5, 8, 8.5, 9, 9.5, 10.};
      task->SetMCEtaScanArray(tmpEtaBinsMC, tmpetaDownArr, tmpetaUpArr);
      // mom bin scan
      const Int_t tmpMomBinsMC = 2;
      Float_t tmppDownArr[tmpMomBinsMC] = { 0.2,0.6};
      Float_t tmppUpArr[tmpMomBinsMC]   = { 1.5,1.5};
      task->SetMCMomScanArray(tmpMomBinsMC, tmppDownArr,   tmppUpArr);
      // cent bins
      const Int_t tmpCentbins  = 10;
      Float_t tmpfxCentBins[tmpCentbins] = {0,5,10,20,30,40,50,60,70,80};
      task->SetCentralityBinning(tmpCentbins,tmpfxCentBins);
      // resonances to exclude
      const Int_t tmpNresonances = 1;
      TString tmpResArr[tmpNresonances] = {"xxx"};
      task->SetMCResonanceArray(tmpNresonances,tmpResArr);
    }
    break;
    case 203:{
      std::cout << settingType << " Info::marsland: Fast Gen on LEGO --> eta and mom scan + removal of resonances 0.2<p<1.5 GeV/c  " << std::endl;
      task->SetRunFastSimulation(kTRUE);
      task->SetRunOnGrid(kTRUE);
      task->SetIsMCtrue(kTRUE);
      // eta bin scan
      const Int_t tmpEtaBinsMC = 20;
      Float_t tmpetaDownArr[tmpEtaBinsMC] = {-0.1,-0.2,-0.3,-0.4,-0.5,-0.6,-0.7,-0.8,-0.9,-1.,-1.1,-1.2,-1.3,-1.4,-1.5,-1.6,-1.7,-1.8,-1.9,-2};
      Float_t tmpetaUpArr[tmpEtaBinsMC]   = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2};
      task->SetMCEtaScanArray(tmpEtaBinsMC, tmpetaDownArr, tmpetaUpArr);
      // mom bin scan
      const Int_t tmpMomBinsMC = 1;
      Float_t tmppDownArr[tmpMomBinsMC] = { 0.2};
      Float_t tmppUpArr[tmpMomBinsMC]   = { 1.5};
      task->SetMCMomScanArray(tmpMomBinsMC, tmppDownArr,   tmppUpArr);
      // cent bins
      const Int_t tmpCentbins  = 10;
      Float_t tmpfxCentBins[tmpCentbins] = {0,5,10,20,30,40,50,60,70,80};
      task->SetCentralityBinning(tmpCentbins,tmpfxCentBins);
      // resonances to exclude
      const Int_t tmpNresonances = 1;
      TString tmpResArr[tmpNresonances] = {"xxx"};
      task->SetMCResonanceArray(tmpNresonances,tmpResArr);
    }
    break;
    case 204:{
      std::cout << settingType << " Info::marsland:  Fast Gen on LEGO --> eta and mom scan + removal of resonances 0.6<p<1.5 GeV/c  " << std::endl;
      task->SetRunFastSimulation(kTRUE);
      task->SetRunOnGrid(kTRUE);
      task->SetIsMCtrue(kTRUE);
      // eta bin scan
      const Int_t tmpEtaBinsMC = 20;
      Float_t tmpetaDownArr[tmpEtaBinsMC] = {-0.1,-0.2,-0.3,-0.4,-0.5,-0.6,-0.7,-0.8,-0.9,-1.,-1.1,-1.2,-1.3,-1.4,-1.5,-1.6,-1.7,-1.8,-1.9,-2};
      Float_t tmpetaUpArr[tmpEtaBinsMC]   = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2};
      task->SetMCEtaScanArray(tmpEtaBinsMC, tmpetaDownArr, tmpetaUpArr);
      // mom bin scan
      const Int_t tmpMomBinsMC = 1;
      Float_t tmppDownArr[tmpMomBinsMC] = { 0.6};
      Float_t tmppUpArr[tmpMomBinsMC]   = { 1.5};
      task->SetMCMomScanArray(tmpMomBinsMC, tmppDownArr,   tmppUpArr);
      // cent bins
      const Int_t tmpCentbins  = 10;
      Float_t tmpfxCentBins[tmpCentbins] = {0,5,10,20,30,40,50,60,70,80};
      task->SetCentralityBinning(tmpCentbins,tmpfxCentBins);
      // resonances to exclude
      const Int_t tmpNresonances = 1;
      TString tmpResArr[tmpNresonances] = {"xxx"};
      task->SetMCResonanceArray(tmpNresonances,tmpResArr);
    }
    break;
    case 205:{
      std::cout << settingType << " Info::marsland: Fast Gen at GSI --> eta and mom scan + removal of resonances in full acceptance [0.2,0.6]<p<1.5 GeV/c  " << std::endl;
      task->SetRunFastSimulation(kTRUE);
      task->SetRunOnGrid(kTRUE);
      task->SetPercentageOfEvents(0);
      task->SetIsMCtrue(kTRUE);
      // eta bin scan
      const Int_t tmpEtaBinsMC = 20;
      Float_t tmpetaDownArr[tmpEtaBinsMC] = {-0.5,-1.,-1.5,-2.,-2.5,-3.,-3.5,-4.,-4.5,-5.,-5.5,-6.,-6.5,-7.,-7.5,-8,-8.5,-9,-9.5,-10.};
      Float_t tmpetaUpArr[tmpEtaBinsMC]   = { 0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5., 5.5, 6., 6.5, 7., 7.5, 8, 8.5, 9, 9.5, 10.};
      task->SetMCEtaScanArray(tmpEtaBinsMC, tmpetaDownArr, tmpetaUpArr);
      // mom bin scan
      const Int_t tmpMomBinsMC = 2;
      Float_t tmppDownArr[tmpMomBinsMC] = { 0.2,0.2};
      Float_t tmppUpArr[tmpMomBinsMC]   = { 3. ,10.};
      task->SetMCMomScanArray(tmpMomBinsMC, tmppDownArr,   tmppUpArr);
      // cent bins
      const Int_t tmpCentbins  = 10;
      Float_t tmpfxCentBins[tmpCentbins] = {0,5,10,20,30,40,50,60,70,80};
      task->SetCentralityBinning(tmpCentbins,tmpfxCentBins);
      // resonances to exclude
      const Int_t tmpNresonances = 1;
      TString tmpResArr[tmpNresonances] = {"xxx"};
      task->SetMCResonanceArray(tmpNresonances,tmpResArr);
    }
    break;
    case 206:{
      std::cout << settingType << " Info::marsland: Fast Gen at GSI --> eta and mom scan + removal of resonances in full acceptance 0<p<1000 GeV/c  " << std::endl;
      task->SetRunFastSimulation(kTRUE);
      task->SetRunOnGrid(kTRUE);
      task->SetPercentageOfEvents(0);
      task->SetIsMCtrue(kTRUE);
      task->SetUseCouts(kTRUE);
      // eta bin scan
      const Int_t tmpEtaBinsMC = 20;
      Float_t tmpetaDownArr[tmpEtaBinsMC] = {-0.5,-1.,-1.5,-2.,-2.5,-3.,-3.5,-4.,-4.5,-5.,-5.5,-6.,-6.5,-7.,-7.5,-8,-8.5,-9,-9.5,-10.};
      Float_t tmpetaUpArr[tmpEtaBinsMC]   = { 0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5., 5.5, 6., 6.5, 7., 7.5, 8, 8.5, 9, 9.5, 10.};
      task->SetMCEtaScanArray(tmpEtaBinsMC, tmpetaDownArr, tmpetaUpArr);
      // mom bin scan
      const Int_t tmpMomBinsMC = 2;
      Float_t tmppDownArr[tmpMomBinsMC] = { 0.};
      Float_t tmppUpArr[tmpMomBinsMC]   = { 10000.};
      task->SetMCMomScanArray(tmpMomBinsMC, tmppDownArr,   tmppUpArr);
      // cent bins
      const Int_t tmpCentbins  = 10;
      Float_t tmpfxCentBins[tmpCentbins] = {0,5,10,20,30,40,50,60,70,80};
      task->SetCentralityBinning(tmpCentbins,tmpfxCentBins);
      // resonances to exclude
      const Int_t tmpNresonances = 1;
      TString tmpResArr[tmpNresonances] = {"xxx"};
      task->SetMCResonanceArray(tmpNresonances,tmpResArr);
      // baryons to be included for netbaryon analysis
      //       const Int_t tmpNbaryons = 7;
      //       Int_t tmpBaryonArr[tmpNbaryons] = {2212,2112,2224,2214,2114,1114,3122};
      //       // {p,n,delta++,delta+,delta0,delta-,Lambda,}
      //       task->SetMCBaryonArray(tmpNbaryons,tmpBaryonArr);
      const Int_t tmpNbaryons = 2;
      Int_t tmpBaryonArr[tmpNbaryons] = {2212,2112};
      // {p,n,delta++,delta+,delta0,delta-,Lambda,}
      task->SetMCBaryonArray(tmpNbaryons,tmpBaryonArr);

    }
    break;
    case 207:{
      std::cout << settingType << " Info::marsland: Fast Gen on LEGO --> Calculate higher moments using look-up table [0.2,0.6]<p<1.5 GeV/c " << std::endl;
      task->SetRunFastHighMomentCal(kTRUE);
      task->SetRunOnGrid(kTRUE);
      task->SetUseCouts(kFALSE);
      task->SetIsMCtrue(kTRUE);
      // eta bin scan
      const Int_t tmpEtaBinsMC = 20;
      Float_t tmpetaDownArr[tmpEtaBinsMC] = {-0.5,-1.,-1.5,-2.,-2.5,-3.,-3.5,-4.,-4.5,-5.,-5.5,-6.,-6.5,-7.,-7.5,-8,-8.5,-9,-9.5,-10.};
      Float_t tmpetaUpArr[tmpEtaBinsMC]   = { 0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5., 5.5, 6., 6.5, 7., 7.5, 8, 8.5, 9, 9.5, 10.};
      task->SetMCEtaScanArray(tmpEtaBinsMC, tmpetaDownArr, tmpetaUpArr);
      // mom bin scan
      const Int_t tmpMomBinsMC = 2;
      Float_t tmppDownArr[tmpMomBinsMC] = { 0.};
      Float_t tmppUpArr[tmpMomBinsMC]   = { 1000.};
      task->SetMCMomScanArray(tmpMomBinsMC, tmppDownArr,   tmppUpArr);
      // cent bins
      const Int_t tmpCentbins  = 10;
      Float_t tmpfxCentBins[tmpCentbins] = {0,5,10,20,30,40,50,60,70,80};
      task->SetCentralityBinning(tmpCentbins,tmpfxCentBins);
      // resonances to exclude
      const Int_t tmpNresonances = 1;
      TString tmpResArr[tmpNresonances] = {"xxx"};
      task->SetMCResonanceArray(tmpNresonances,tmpResArr);
      // Read First Moments
      task->SetLookUpTableFirstMoments(lookUpTree,0, tmppDownArr,tmpfxCentBins,tmpetaUpArr,tmpMomBinsMC,tmpCentbins,tmpEtaBinsMC);
      task->SetLookUpTableFirstMoments(lookUpTree,1, tmppDownArr,tmpfxCentBins,tmpetaUpArr,tmpMomBinsMC,tmpCentbins,tmpEtaBinsMC);
      task->SetLookUpTableFirstMoments(lookUpTree,2, tmppDownArr,tmpfxCentBins,tmpetaUpArr,tmpMomBinsMC,tmpCentbins,tmpEtaBinsMC);
      task->SetLookUpTableFirstMoments(lookUpTree,9, tmppDownArr,tmpfxCentBins,tmpetaUpArr,tmpMomBinsMC,tmpCentbins,tmpEtaBinsMC);
      task->SetLookUpTableFirstMoments(lookUpTree,11,tmppDownArr,tmpfxCentBins,tmpetaUpArr,tmpMomBinsMC,tmpCentbins,tmpEtaBinsMC);
    }
    break;
    case 208:{
      std::cout << settingType << " Info::marsland: Higher Moments MC closure " << std::endl;
      task->SetRunOnGrid(getFromAlien);
      std::cout << settingType << " Info::marsland: runOnGrid = " << getFromAlien << std::endl;
      task->SetUseCouts(kTRUE);
      task->SetIsMCtrue(kTRUE);
      task->SetFillHigherMomentsMCclosure(kTRUE);
      // eta cent and mom scan
      const Int_t tmpEtaBinsMC = 8;
      const Int_t tmpMomBinsMC = 4;
      const Int_t tmpCentbins  = 10;
      Float_t tmpfxCentBins[tmpCentbins]  = {0,5,10,20,30,40,50,60,70,80};
      Float_t tmpetaDownArr[tmpEtaBinsMC] = {-0.1,-0.2,-0.3,-0.4,-0.5,-0.6,-0.7,-0.8};
      Float_t tmpetaUpArr[tmpEtaBinsMC]   = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};
      Float_t tmppDownArr[tmpMomBinsMC]   = { 0.4, 0.4, 0.6, 0.6};
      Float_t tmppUpArr[tmpMomBinsMC]     = { 0.8, 1.0, 1.5, 2.0};
      task->SetMCEtaScanArray(tmpEtaBinsMC, tmpetaDownArr, tmpetaUpArr);
      task->SetMCMomScanArray(tmpMomBinsMC, tmppDownArr,   tmppUpArr);
      task->SetCentralityBinning(tmpCentbins,tmpfxCentBins);
      // Read First Moments for higher moments
      task->SetLookUpTableEfficiencyCorrection(lookUpTree, 0, tmppDownArr, tmppUpArr, tmpfxCentBins, tmpetaUpArr, tmpMomBinsMC, tmpCentbins, tmpEtaBinsMC);
      task->SetLookUpTableEfficiencyCorrection(lookUpTree, 1, tmppDownArr, tmppUpArr, tmpfxCentBins, tmpetaUpArr, tmpMomBinsMC, tmpCentbins, tmpEtaBinsMC);
      task->SetLookUpTableEfficiencyCorrection(lookUpTree, 2, tmppDownArr, tmppUpArr, tmpfxCentBins, tmpetaUpArr, tmpMomBinsMC, tmpCentbins, tmpEtaBinsMC);
    }

  }

  // Finally initialize the task after physics selection
  task->Initialize();
  return task;

}
// ____________________________________________________________________________________________
void SetDefaults(AliAnalysisTaskEbyeIterPID *defaultTask)
{

  // Setters for the eta momentum dEdx and centrality bins
  std::cout << " Info::marsland: ------------------------------------------------------------------------------------- " << std::endl;
  std::cout << " Info::marsland: ------------------------------------------------------------------------------------- " << std::endl;
  std::cout << " Info::marsland: ------------------- Set default settings for the task object ------------------------ " << std::endl;
  std::cout << " Info::marsland: ------------------------------------------------------------------------------------- " << std::endl;
  std::cout << " Info::marsland: ------------------------------------------------------------------------------------- " << std::endl;

  defaultTask->SetSampleDeDxUpperEdge(400.);
  defaultTask->SetDeDxBinWidth(2.5);
  defaultTask->SetDeDxLowerEdge(20.);
  defaultTask->SetDeDxUpperEdge(2020.);
  defaultTask->SetNEtabins(20);
  defaultTask->SetEtaLowerEdge(-1.);
  defaultTask->SetEtaUpperEdge( 1.);
  defaultTask->SetNMomBins(150);
  defaultTask->SetMomLowerEdge(0.1);
  defaultTask->SetMomUpperEdge(3.1);
  defaultTask->SetNGenprotonBins(100);

  // DEFAULT SETTINGS
  const Int_t tmpCentbins  = 10;
  const Int_t tmpEtaBinsMC = 3;
  const Int_t tmpMomBinsMC = 4;
  Float_t tmpfxCentBins[tmpCentbins] = {0,5,10,20,30,40,50,60,70,80};
  Float_t tmpetaDownArr[tmpEtaBinsMC] = {-0.5,-0.8,-1.};
  Float_t tmpetaUpArr[tmpEtaBinsMC]   = { 0.5, 0.8, 1.};
  Float_t tmppDownArr[tmpMomBinsMC] = { 0.2, 0.6, 0.2, 0.6};
  Float_t tmppUpArr[tmpMomBinsMC]   = { 1.5, 1.5, 1.8, 1.8};

  // centrality binning and Eta Momentum Scans for MC
  defaultTask->SetCentralityBinning(tmpCentbins,tmpfxCentBins);
  defaultTask->SetMCEtaScanArray(tmpEtaBinsMC, tmpetaDownArr, tmpetaUpArr);
  defaultTask->SetMCMomScanArray(tmpMomBinsMC, tmppDownArr,   tmppUpArr);

  // Boolians which are by default === ON ===
  defaultTask->SetRunOnGrid(kFALSE);
  defaultTask->SetIsMCtrue(kFALSE);
  defaultTask->SetIncludeITScuts(kFALSE);
  defaultTask->SetFillArmPodTree(kFALSE);
  defaultTask->SetUsePtCut(1);
  defaultTask->SetTrackOriginType(0);   // 0:prim, 1: prim+weak, 2:prim+material, 3:prim+material+weak, 4:full scan
  defaultTask->SetRapidityType(0);      // 0:pseudorapidity, 1: rapidity

  // Extra Boolians which are by default === OFF ===
  defaultTask->SetDeDxCheck(kFALSE);
  defaultTask->SetEffMatrix(kFALSE);
  defaultTask->SetFillAllCutVariables(kFALSE);
  defaultTask->SetFillOnlyHists(kFALSE);
  defaultTask->SetRunFastSimulation(kFALSE);
  defaultTask->SetFillHigherMomentsMCclosure(kFALSE);
  defaultTask->SetFillDnchDeta(kFALSE);
  defaultTask->SetIncludeTOF(kFALSE);
  defaultTask->SetUseThnSparse(kFALSE);
  defaultTask->SetUseCouts(kFALSE);
  defaultTask->SetWeakAndMaterial(kFALSE);
  defaultTask->SetFillEventInfo(kFALSE);
  defaultTask->SetFillTreeMC(kFALSE);
  defaultTask->SetFillAllCutVariables(kFALSE);
  defaultTask->SetFillNudynFastGen(kFALSE);
  defaultTask->SetDefaultTrackCuts(kTRUE);
  defaultTask->SetDefaultEventCuts(kFALSE);

  // Setters for the systematic uncertainty checks
  defaultTask->SetSystCentEstimator(0);
  defaultTask->SetSystDCAxy(0);
  defaultTask->SetSystNCrossedRows(0);
  defaultTask->SetSystTPCChi2(0);
  defaultTask->SetSystVz(0);

}
// ____________________________________________________________________________________________
TTree *GetLookUpTable(Bool_t runOnGrid, Int_t index)
{

  std::cout << " Info::marsland: Copy LookUp table from alien " << std::endl;
  // Define the lookup table file name
  TString fileName  ="";
  if (index==1) fileName="MomentsTree_AccCan_HIJING.root";
  if (index==2) fileName="MomentsTree_AccCan_LHC13f3a.root";
  if (index==3) fileName="MomentsTree_AccCan_LHC13f3b.root";
  if (index==4) fileName="MomentsTree_AccCan_LHC13f3c.root";
  if (index==5) fileName="MomentsTree_AccCan_HIJING_FullAcc.root";
  if (index==6) fileName="MomentsTree_AccCan_LHC13f3a_FullAcc.root";
  if (index==7) fileName="MomentsTree_AccCan_LHC13f3b_FullAcc.root";
  if (index==8) fileName="MomentsTree_AccCan_LHC13f3c_FullAcc.root";
  if (index==9) fileName="LookupTable_HighOrderCumulants_10h_HIJING.root";
  if (index==10) fileName="LookupTable_res_0_cent_HighOrderCumulants_10h_HIJING.root";
  if (index==11) fileName="LookupTable_res_0_centimp_HighOrderCumulants_10h_HIJING.root";

  //
  TTree *tree=NULL;
  TFile *fInputLookUp=NULL;
  TString lookUpDir ="";
  TString lookUpPath="";
  //
  // connect to alien for the lookup table
  std::cout << " Info::marsland: Connecting to GRID for the lookup table " << std::endl;

  if (runOnGrid){
    TGrid * alien = TGrid::Connect("alien://",0,0,"t");
    lookUpDir = "alien:///alice/cern.ch/user/m/marsland/PWGCF/EBYE/IdentityMethodEbyeFluctuations/macros";
    gSystem->Exec(Form("alien_cp %s/%s .",lookUpDir.Data(),fileName.Data()));
    lookUpPath = Form("%s/%s",gSystem->pwd(),fileName.Data());
  } else {
    lookUpDir = "/lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/MC/LookUpTables";
    lookUpPath = Form("%s/%s",lookUpDir.Data(),fileName.Data());
  }

  //
  // retrieve the ttree
  std::cout << " Info::marsland: LookUp table used is = " << lookUpPath << std::endl;
  fInputLookUp = new TFile(lookUpPath);
  tree = (TTree*)fInputLookUp->Get("mcMoments");
  //
  // return the lookup tree for further processing
  if(tree) return tree;
  else { std::cout << " Error::marsland: There is no lookUp table" << std::endl; return 0;}

}
