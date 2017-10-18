//
//
void SetDefaults(AliAnalysisTaskEbyeIterPID *defaultTask);
//
//
//
AliAnalysisTaskEbyeIterPID* Config_marsland_EbyeIterPID(Int_t settingType, Int_t lhcPeriod) {
  
  // 
  // Configuration file for the AliAnalysisTaskEbyeIterPID.cxx class
  //
  /*
   *    Real data --> settingType = 
   *    0.)  THnSparse is used: StandardTPCITScuts              16EtaBin_150pBins_9centBins (REFERENCE settings)
   *    1.)  THnSparse is used: StandardTPCITScuts + TIGHTCUTS  16EtaBin_150pBins_9centBins (active length cut)
   *    2.)  THnSparse is used: StandardTPCITScuts + DCAxySmall 16EtaBin_150pBins_9centBins  
   *    3.)  THnSparse is used: StandardTPCITScuts + DCAxyLarge 16EtaBin_150pBins_9centBins  
   *    4.)  THnSparse is used: StandardTPCITScuts + cRows60    16EtaBin_150pBins_9centBins  
   *    5.)  THnSparse is used: StandardTPCITScuts + cRows100   16EtaBin_150pBins_9centBins  
   *    6.)  THnSparse is used: StandardTPCITScuts + centEstCL1 16EtaBin_150pBins_9centBins  
   *    7.)  THnSparse is used: StandardTPCITScuts + centEstTRK 16EtaBin_150pBins_9centBins  
   *    8.)  THnSparse is used: StandardTPCITScuts + Vz8        16EtaBin_150pBins_9centBins  
   *    9.)  THnSparse is used: StandardTPCITScuts + Vz12       16EtaBin_150pBins_9centBins  
   *    10.) THnSparse is used: StandardTPCITScuts + Chi2Small  16EtaBin_150pBins_9centBins  
   *    11.) THnSparse is used: StandardTPCITScuts + Chi2Large  16EtaBin_150pBins_9centBins  
   *    12.) THnSparse is used: (REFERENCE settings) + allCuts are filled  
   *    13.) THnSparse is used: (REFERENCE settings) + Bayesian Probabilities are filled  
   *    14.) THnSparse is used: (REFERENCE settings) + dEdxTree is filled
   *    15.) THnSparse is used: (REFERENCE settings) + 18centBins  
   *    16.) THnSparse is used: (REFERENCE settings) + centBin 10  
   *    17.) THnSparse is used: (REFERENCE settings) + centBin 5  
   *    18.) THnSparse is used: ITS is OFF 
   *    19.) THnSparse is used: (REFERENCE settings) + THnSparse is used: number of eta bins = 32 
   *    20.) THnSparse is used: StandardTPCITScuts + TightCuts TPC dEdx preliminary plot 
   *    21.) THnSparse is used: (REFERENCE settings) + dEdx + allCuts + ArmPodTree filled + eta range extended 
   * 
   *    MC data --> settingType = 
   *    100.) THnSparse is used: StandardTPCITScuts 8EtaBin_150pBins_9centBins (REFERENCE settings) MC CLOSURE
   *    101.) FastSimul is used: StandardTPCITScuts 8EtaBin_150pBins_9centBins (REFERENCE settings)
   *    102.) FastSimul is used: StandardTPCITScuts 8EtaBin_150pBins_9centBins ETA DEPENDENCE
   *    103.) FastSimul is used: StandardTPCITScuts 8EtaBin_150pBins_9centBins Momentum DEPENDENCE
   *    104.) FullSinul is used: StandardTPCITScuts 16EtaBin_150pBins_9centBins EffMatrix
   *    105.) FullSinul is used: Tight Cuts         16EtaBin_150pBins_9centBins EffMatrix
   *    106.) FastSimul is used: StandardTPCITScuts 8EtaBin_150pBins_9centBins (REFERENCE settings): Fill only DnchDeta
   *    107.) THnSparse is used: StandardTPCITScuts 8EtaBin_150pBins_9centBins (REFERENCE settings) MC CLOSURE with TOF cut
   *    108.) THnSparse is used: StandardTPCITScuts 8EtaBin_150pBins_9centBins (REFERENCE settings) MC CLOSURE with NO ITS cut 
   *    109.) THnSparse is used: (REFERENCE settings) + dEdx + allCuts + ArmPodTree filled + eta range extended
   */
  AliAnalysisTaskEbyeIterPID *task = new AliAnalysisTaskEbyeIterPID("marslandEbyeIterPID");
  SetDefaults(task);
  if (lhcPeriod==1) task->SelectCollisionCandidates(AliVEvent::kMB);   // select minimum bias events for LHC10h
  if (lhcPeriod==2) task->SelectCollisionCandidates(AliVEvent::kINT7); // select minimum bias events for LHC15o
  std::cout << " ===== In the Config --> Running with lhcPeriod =  " << lhcPeriod << " ===== " << std::endl;
  // Other Specific settings

  switch (settingType) {
    //     
    // ====================================================================================
    // ============================= Real Data Settings ===================================
    // ====================================================================================
    // 
    // THnSparse is used: StandardTPCITScuts 16EtaBin_150pBins_9centBins (REFERENCE settings)
    case 0:{ 
      std::cout << settingType << " THnSparse is used: StandardTPCITScuts 16EtaBin_150pBins_9centBins (REFERENCE settings) " << std::endl;
      task->SetFillDeDxTree(kFALSE); 
    }
    break;
    case 1:{ 
      std::cout << settingType << " THnSparse is used: (REFERENCE settings)  + TIGHTCUTS  (active length cut) " << std::endl;
      task->SetFillDeDxTree(kFALSE);  
      task->SetTightCuts(kTRUE); 
    }
    break;
    case 2:{ 
      std::cout << settingType << " THnSparse is used: (REFERENCE settings)  + DCAxySmall " << std::endl;
      task->SetFillDeDxTree(kFALSE);  
      task->SetSystDCAxy(-1); 
    }
    break;
    case 3:{ 
      std::cout << settingType << " THnSparse is used: (REFERENCE settings)  + DCAxyLarge " << std::endl;
      task->SetFillDeDxTree(kFALSE);  
      task->SetSystDCAxy(1); 
    }
    break;
    case 4:{ 
      std::cout << settingType << " THnSparse is used: (REFERENCE settings)  + cRows60 " << std::endl;
      task->SetFillDeDxTree(kFALSE);  
      task->SetSystNCrossedRows(-1); 
    }
    break;
    case 5:{ 
      std::cout << settingType << " THnSparse is used: (REFERENCE settings)  + cRows100 " << std::endl;
      task->SetFillDeDxTree(kFALSE);  
      task->SetSystNCrossedRows(1); 
    }
    break;
    case 6:{ 
      std::cout << settingType << " THnSparse is used: (REFERENCE settings)  + centEstCL1 " << std::endl;
      task->SetFillDeDxTree(kFALSE);  
      task->SetSystCentEstimator(-1); 
    }
    break;
    case 7:{ 
      std::cout << settingType << " THnSparse is used: (REFERENCE settings)  + centEstTRK " << std::endl;
      task->SetFillDeDxTree(kFALSE);  
      task->SetSystCentEstimator(1); 
    }
    break;
    case 8:{ 
      std::cout << settingType << " THnSparse is used: (REFERENCE settings)  + Vz8 " << std::endl;
      task->SetFillDeDxTree(kFALSE);  
      task->SetSystVz(-1); 
    }
    break;
    case 9:{ 
      std::cout << settingType << " THnSparse is used: (REFERENCE settings)  + Vz12 " << std::endl;
      task->SetFillDeDxTree(kFALSE);  
      task->SetSystVz(1); 
    }
    break;
    case 10:{ 
      std::cout << settingType << " THnSparse is used: (REFERENCE settings)  + Chi2Small " << std::endl;
      task->SetFillDeDxTree(kFALSE); 
      task->SetSystTPCChi2(-1); 
    }
    break;
    case 11:{ 
      std::cout << settingType << " THnSparse is used: (REFERENCE settings)  + Chi2Large " << std::endl;
      task->SetFillDeDxTree(kFALSE); 
      task->SetSystTPCChi2(1); 
    }
    break;
    case 12:{ 
      std::cout << settingType << " THnSparse is used: (REFERENCE settings) + allCuts are filled " << std::endl;
      task->SetFillDeDxTree(kFALSE); 
      task->SetFillAllCutVariables(kTRUE); 
      task->SetFillArmPodTree(kTRUE); 
    }
    break;
    case 13:{ 
      std::cout << settingType << " THnSparse is used: (REFERENCE settings) + Bayesian Probabilities are filled " << std::endl;
      task->SetFillDeDxTree(kFALSE); 
      task->SetFillBayesianProb(kTRUE); 
    }
    break;
    case 14:{ 
      std::cout << settingType << " (REFERENCE settings) + dEdx Tree is filled " << std::endl;
      task->SetFillDeDxTree(kTRUE); 
    }
    break;
    case 15:{ 
      std::cout << settingType << " THnSparse is used: (REFERENCE settings) + 18centBins " << std::endl;
      task->SetFillDeDxTree(kFALSE); 
      const Int_t tmpCentbins = 19; 
      Float_t tmpfxCentBins[tmpCentbins] = {0,2.5,5,7.5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80}; 
      task->SetCentralityBinning(tmpCentbins,tmpfxCentBins); 
    }
    break;
    case 16:{ 
      std::cout << settingType << " THnSparse is used: (REFERENCE settings) + centBin 10 " << std::endl;
      task->SetFillDeDxTree(kFALSE); 
      const Int_t tmpCentbins = 9;  
      Float_t tmpfxCentBins[tmpCentbins] = {0,10,20,30,40,50,60,70,80}; 
      task->SetCentralityBinning(tmpCentbins,tmpfxCentBins); 
    }
    break;
    case 17:{ 
      std::cout << settingType << " THnSparse is used: (REFERENCE settings) + centBin 5 " << std::endl;
      task->SetFillDeDxTree(kFALSE); 
      const Int_t tmpCentbins = 17; 
      Float_t tmpfxCentBins[tmpCentbins] = {0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80}; 
      task->SetCentralityBinning(tmpCentbins,tmpfxCentBins); 
    }
    break;
    case 18:{ 
      std::cout << settingType << " THnSparse is used: ITS is OFF " << std::endl;
      task->SetFillDeDxTree(kFALSE);
      task->SetSystDCAxy(1); 
      task->SetSystVz(1); 
    }
    break;
    case 19:{ 
      std::cout << settingType << " THnSparse is used: (REFERENCE settings) + THnSparse is used: number of eta bins = 32 " << std::endl;
      task->SetFillDeDxTree(kFALSE);
      task->SetNEtabins(32); 
      const Int_t tmpCentbins  = 10;   
      Float_t tmpfxCentBins[tmpCentbins] = {0,5,10,20,30,40,50,60,70,80};
      task->SetCentralityBinning(tmpCentbins,tmpfxCentBins); 
    }
    break;
    case 20:{ 
      std::cout << settingType << " THnSparse is used: StandardTPCITScuts + TightCuts TPC dEdx preliminary plot  " << std::endl;
      task->SetDeDxBinWidth(2.5);
      task->SetDeDxLowerEdge(20.);
      task->SetDeDxUpperEdge(1020.);
      task->SetNEtabins(8);
      task->SetEtaLowerEdge(-0.8);
      task->SetEtaUpperEdge(0.8);
      task->SetNMomBins(200);
      task->SetMomLowerEdge(0.);
      task->SetMomUpperEdge(200.);
      task->SetFillDeDxTree(kTRUE); 
      task->SetTightCuts(kTRUE);
    }
    break;
    case 21:{ 
      std::cout << settingType << " THnSparse is used: (REFERENCE settings) + allCuts + ArmPodTree filled + eta range extended  " << std::endl;
      task->SetFillAllCutVariables(kTRUE); 
      task->SetFillTIdenTrees(kTRUE);
      task->SetFillArmPodTree(kTRUE); 
      task->SetUseCouts(kTRUE);  
      task->SetNEtabins(10);
      task->SetEtaLowerEdge(-1.);
      task->SetEtaUpperEdge(1.);
    }
    break;
    //     
    // ====================================================================================
    // =================================== MC Settings ====================================
    // ====================================================================================
    // 
    case 100:{ 
      std::cout << settingType << " THnSparse is used: StandardTPCITScuts 8EtaBin_150pBins_9centBins (REFERENCE settings) MC CLOSURE " << std::endl;
      task->SetFillDeDxTree(kFALSE);   
      task->SetIsMCtrue(kTRUE); 
      task->SetIncludeTOF(kFALSE); 
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
    case 101:{ 
      std::cout << settingType << " FastSimul is used: StandardTPCITScuts 8EtaBin_150pBins_9centBins (REFERENCE settings) " << std::endl;
      task->SetIsMCtrue(kTRUE); 
      task->SetNEtabins(8);
      task->SetRunFastSimulation(kTRUE);
      task->SetCentralityBinning(tmpCentbins,tmpfxCentBins); 
    }
    break;
    case 102:{ 
      std::cout << settingType << " FastSimul is used: StandardTPCITScuts 8EtaBin_150pBins_9centBins ETA DEPENDENCE " << std::endl;
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
    case 103:{ 
      std::cout << settingType << " FastSimul is used: StandardTPCITScuts 8EtaBin_150pBins_9centBins Momentum DEPENDENCE " << std::endl;
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
    case 104:{ 
      std::cout << settingType << " FullSimul is used: StandardTPCITScuts 8EtaBin_150pBins_9centBins EffMatrix " << std::endl;
      task->SetEffMatrix(kTRUE);
      task->SetIsMCtrue(kTRUE); 
      task->SetTightCuts(kFALSE);
      task->SetNEtabins(8);
      task->SetRunFastSimulation(kFALSE);
      task->SetFillArmPodTree(kFALSE);
      task->SetFillDeDxTree(kFALSE);
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
    case 105:{ 
      std::cout << settingType << " FullSimul is used: StandardTPCITScuts 8EtaBin_150pBins_9centBins EffMatrix " << std::endl;
      task->SetEffMatrix(kTRUE);
      task->SetIsMCtrue(kTRUE); 
      task->SetTightCuts(kFALSE);
      task->SetNEtabins(8);
      task->SetRunFastSimulation(kFALSE);
      task->SetFillArmPodTree(kFALSE);
      task->SetFillDeDxTree(kFALSE);
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
    case 106:{ 
      std::cout << settingType << " FullSimul is used: Tight Cuts 8EtaBin_150pBins_9centBins EffMatrix " << std::endl;
      task->SetEffMatrix(kTRUE);
      task->SetIsMCtrue(kTRUE); 
      task->SetNEtabins(8);
      task->SetRunFastSimulation(kFALSE);
      task->SetTightCuts(kTRUE);
      task->SetFillArmPodTree(kFALSE);
      task->SetFillDeDxTree(kFALSE);
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
    case 107:{ 
      std::cout << settingType << " FastSimul is used: StandardTPCITScuts 8EtaBin_150pBins_9centBins (REFERENCE settings): Fill only DnchDeta " << std::endl;
      task->SetIsMCtrue(kTRUE); 
      task->SetNEtabins(8);
      task->SetRunFastSimulation(kTRUE);
      task->SetFillDnchDeta(kTRUE);
    }
    break;
    case 108:{ 
      std::cout << settingType << " THnSparse is used: StandardTPCITScuts 8EtaBin_150pBins_9centBins (REFERENCE settings) MC CLOSURE with TOF cut " << std::endl;
      task->SetFillDeDxTree(kFALSE);  
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
    case 109:{ 
      std::cout << settingType << " THnSparse is used: StandardTPCITScuts 8EtaBin_150pBins_9centBins (REFERENCE settings) MC CLOSURE with NO ITS cut " << std::endl;
      task->SetFillDeDxTree(kFALSE);  
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
    case 110:{ 
      std::cout << settingType << " FULL MC --> eta and mom scan + removal of resonances + EffMatrix  " << std::endl;
      task->SetIsMCtrue(kTRUE); 
      task->SetIncludeTOF(kFALSE); 
      task->SetEffMatrix(kTRUE);
      task->SetNEtabins(10);
      task->SetEtaLowerEdge(-1.);
      task->SetEtaUpperEdge(1.);
      // eta bin scan
      const Int_t tmpEtaBinsMC = 2;
      Float_t tmpetaDownArr[tmpEtaBinsMC] = {-0.8,-1.0};
      Float_t tmpetaUpArr[tmpEtaBinsMC]   = { 0.8, 1.0};
      task->SetMCEtaScanArray(tmpEtaBinsMC, tmpetaDownArr, tmpetaUpArr);
      // mom bin scan
      const Int_t tmpMomBinsMC = 1;
      Float_t tmppDownArr[tmpMomBinsMC] = {0.2};
      Float_t tmppUpArr[tmpMomBinsMC]   = {2.2};
      task->SetMCMomScanArray(tmpMomBinsMC, tmppDownArr,   tmppUpArr); 
      // cent bins
      const Int_t tmpCentbins  = 10;   
      Float_t tmpfxCentBins[tmpCentbins] = {0,5,10,20,30,40,50,60,70,80};
      task->SetCentralityBinning(tmpCentbins,tmpfxCentBins);
      // resonances to exclude
      const Int_t tmpNresonances = 3;
      TString tmpResArr[tmpNresonances] = {"rho","phi","Delta"};
      task->SetMCResonanceArray(tmpNresonances,tmpResArr);
    }
    break;
    case 111:{ 
      std::cout << settingType << " Fast Gen --> eta and mom scan + removal of resonances  " << std::endl;
      task->SetRunFastSimulation(kTRUE);
      task->SetIsMCtrue(kTRUE); 
      task->SetNEtabins(20);
      task->SetEtaLowerEdge(-2.);
      task->SetEtaUpperEdge(2.);
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
      //       const Int_t tmpNresonances = 1;
      //       TString tmpResArr[tmpNresonances] = {"xxx"};
      //       task->SetMCResonanceArray(tmpNresonances,tmpResArr);
      const Int_t tmpNresonances = 5;
      TString tmpResArr[tmpNresonances] = {"rho","phi","omega","eta","Delta"};
      task->SetMCResonanceArray(tmpNresonances,tmpResArr); 
    }
    break;
    case 112:{ 
      std::cout << settingType << " Fast Gen --> eta and mom scan + removal of resonances  " << std::endl;
      task->SetRunFastSimulation(kTRUE);
      task->SetIsMCtrue(kTRUE);
      task->SetNEtabins(20);
      task->SetEtaLowerEdge(-2.);
      task->SetEtaUpperEdge(2.);
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
      //       const Int_t tmpNresonances = 1;
      //       TString tmpResArr[tmpNresonances] = {"xxx"};
      //       task->SetMCResonanceArray(tmpNresonances,tmpResArr);
      const Int_t tmpNresonances = 5;
      TString tmpResArr[tmpNresonances] = {"rho","phi","omega","eta","Delta"};
      task->SetMCResonanceArray(tmpNresonances,tmpResArr);
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
  std::cout << " ------------------------------------------------------------------------------------- " << std::endl;
  std::cout << " ------------------------------------------------------------------------------------- " << std::endl;
  std::cout << " ------------------- Set default settings for the task object ------------------------ " << std::endl;
  std::cout << " ------------------------------------------------------------------------------------- " << std::endl;
  std::cout << " ------------------------------------------------------------------------------------- " << std::endl;

  defaultTask->SetSampleDeDxUpperEdge(140.);
  defaultTask->SetDeDxBinWidth(2.5);
  defaultTask->SetDeDxLowerEdge(20.);
  defaultTask->SetDeDxUpperEdge(1020.);
  defaultTask->SetNEtabins(16);
  defaultTask->SetEtaLowerEdge(-0.8);
  defaultTask->SetEtaUpperEdge(0.8);
  defaultTask->SetNMomBins(150);
  defaultTask->SetMomLowerEdge(0.2);
  defaultTask->SetMomUpperEdge(3.2);
  
  // DEFAULT SETTINGS
  const Int_t tmpCentbins  = 10;   
  const Int_t tmpEtaBinsMC = 2;
  const Int_t tmpMomBinsMC = 4;
  Float_t tmpfxCentBins[tmpCentbins] = {0,5,10,20,30,40,50,60,70,80};
  Float_t tmpetaDownArr[tmpEtaBinsMC] = {-0.5,-0.8};
  Float_t tmpetaUpArr[tmpEtaBinsMC]   = { 0.5, 0.8};
  Float_t tmppDownArr[tmpMomBinsMC] = { 0.2, 0.2, 0.3, 0.3};
  Float_t tmppUpArr[tmpMomBinsMC]   = { 1.5, 1.8, 1.5, 1.8};
  
  // centrality binning and Eta Momentum Scans for MC
  defaultTask->SetCentralityBinning(tmpCentbins,tmpfxCentBins);
  defaultTask->SetMCEtaScanArray(tmpEtaBinsMC, tmpetaDownArr, tmpetaUpArr);
  defaultTask->SetMCMomScanArray(tmpMomBinsMC, tmppDownArr,   tmppUpArr);  
  
  // Boolians which are by default === ON ===
  defaultTask->SetIsMCtrue(kFALSE);
  defaultTask->SetIncludeITScuts(kFALSE);
  defaultTask->SetFillArmPodTree(kFALSE);
  defaultTask->SetFillDeDxTree(kFALSE);
  
  // Extra Boolians which are by default === OFF ===
  defaultTask->SetTightCuts(kFALSE);
  defaultTask->SetDeDxCheck(kFALSE);
  defaultTask->SetEffMatrix(kFALSE);
  defaultTask->SetCleanSamplesOnly(kFALSE);
  defaultTask->SetFillBayesianProb(kFALSE);
  defaultTask->SetFillAllCutVariables(kFALSE);
  defaultTask->SetRunFastSimulation(kFALSE);
  defaultTask->SetFillDnchDeta(kFALSE);
  defaultTask->SetIncludeTOF(kFALSE);   
  defaultTask->SetUseThnSparse(kFALSE);  
  defaultTask->SetUseCouts(kFALSE);  
  defaultTask->SetWeakAndMaterial(kFALSE);
  defaultTask->SetFillTIdenTrees(kFALSE);
  
  // Setters for the systematic uncertainty checks
  defaultTask->SetSystCentEstimator(0);
  defaultTask->SetSystDCAxy(0);
  defaultTask->SetSystNCrossedRows(0);
  defaultTask->SetSystTPCChi2(0);
  defaultTask->SetSystVz(0);
  
}
