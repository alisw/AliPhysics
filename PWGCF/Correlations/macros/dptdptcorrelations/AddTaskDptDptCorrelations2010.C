//
// Macro designed for use with the AliAnalysisTaskDptDptCorrelations task.
//
// Author: Prabhat Pujahari & Claude Pruneau, Wayne State
// 
//           system:  0: PbPb                 1: pp
//      singlesOnly:  0: full correlations    1: singles only
//       useWeights:  0: no                   1: yes
// centralityMethod:  3: track count  4: V0 centrality
//        chargeSet:  0: ++    1: +-    2: -+    3: --
/////////////////////////////////////////////////////////////////////////////////
AliAnalysisTaskDptDptCorrelations *AddTaskDptDptCorrelations2010(int    system                 = 0,
                                                                 int    singlesOnly            = 0,
                                                                 int    useWeights             = 1,
                                                                 int    centralityMethod       = 4,
                                                                 int    chargeSet              = 1)

{
  // Set Default Configuration of this analysis
  // ==========================================
  int    debugLevel             = 0;
  int    rejectPileup           = 1;
  int    rejectPairConversion   = 1;
  int    sameFilter             = 1;

  
  int    nCentrality;
  double minCentrality[10];
  double maxCentrality[10];

  if (system==0) // PbPb
    {
    if (centralityMethod == 4)
      {
      nCentrality = 10;
      minCentrality[0] = 0.0; maxCentrality[0] = 5.0;
      minCentrality[1] = 5.0; maxCentrality[1] = 10.;
      minCentrality[2] = 10.; maxCentrality[2] = 20.;
      minCentrality[3] = 20.; maxCentrality[3] = 30.;
      minCentrality[4] = 30.; maxCentrality[4] = 40.;
      minCentrality[5] = 40.; maxCentrality[5] = 50.;
      minCentrality[6] = 50.; maxCentrality[6] = 60.;
      minCentrality[7] = 60.; maxCentrality[7] = 70.;
      minCentrality[8] = 70.; maxCentrality[8] = 80.;
      minCentrality[9] = 80.; maxCentrality[9] = 90.;
      
      }
    else
      {
      cout << "-F- AddTaskDptDptCorrelations() system:" << system << ". centralityMethod:" << centralityMethod << " Option NOT AVAILABLE. ABORT."
      return 0;
      }
    }
  else if (system==1) // pp
    {
    if (centralityMethod == 3)
      {
      nCentrality = 4;
      minCentrality[0] = 2;   maxCentrality[0] = 100.0;
      minCentrality[1] = 2;   maxCentrality[1] = 20.;
      minCentrality[2] = 20.; maxCentrality[2] = 50.;
      minCentrality[3] = 50.; maxCentrality[3] = 100.;
      }
    else
      {
      cout << "-F- AddTaskDptDptCorrelations() system:" << system << ". centralityMethod:" << centralityMethod << " Option NOT AVAILABLE. ABORT."
      return 0;
      }
    }
  else
    {
    cout << "-F- AddTaskDptDptCorrelations() system:" << system << ". Option NOT CURRENTLY AVAILABLE. ABORT."
    return 0;
    }

  double zMin                   = -10.;
  double zMax                   =  10.;
  double ptMin                  =  0.2;
  double ptMax                  =  2.0;
  double etaMin                 = -1.0;
  double etaMax                 =  1.0;
  double dcaZMin                = -3.0;
  double dcaZMax                =  3.0;
  double dcaXYMin               = -2.4;
  double dcaXYMax               =  2.4;
  double dedxMin                =  0.0;
  double dedxMax                =  20000.0;
  int    nClusterMin            =   70;
  int    trackFilterBit         =  128;
  
  int    requestedCharge1       =  1; //default
  int    requestedCharge2       = -1; //default
  
  
  // Get the pointer to the existing analysis manager via the static access method.
  // ==============================================================================
  AliAnalysisManager *analysisManager = AliAnalysisManager::GetAnalysisManager();
  
  if (!analysisManager) 
    {
    ::Error("AddTaskDptDptCorrelations", "No analysis manager to connect to.");
    return NULL;
    }  
  
  TString part1Name;
  TString part2Name;
  TString eventName;
  TString prefixName        = "Corr_";
  TString pileupRejecSuffix = "_PileupRejec";
  TString pairRejecSuffix   = "_PairRejec";
  TString calibSuffix       = "_calib";
  TString singlesOnlySuffix = "_SO";
  TString suffix;
  
  TString inputPath         = ".";
  TString outputPath        = ".";
  TString baseName;
  TString listName;
  TString taskName;
  TString inputHistogramFileName;
  TString outputHistogramFileName;
  
  // Create the task and add subtask.
  // ===========================================================================
  int iTask = 0; // task counter
  AliAnalysisDataContainer *taskInputContainer;
  AliAnalysisDataContainer *taskOutputContainer;
  AliAnalysisTaskDptDptCorrelations* task;
  
  for (int iCentrality=0; iCentrality < nCentrality; ++iCentrality)
    {
      switch (chargeSet)
        {
          case 0: part1Name = "P_"; part2Name = "P_"; requestedCharge1 =  1; requestedCharge2 =  1; sameFilter = 1; break;
          case 1: part1Name = "P_"; part2Name = "M_"; requestedCharge1 =  1; requestedCharge2 = -1; sameFilter = 0;   break;
          case 2: part1Name = "M_"; part2Name = "P_"; requestedCharge1 = -1; requestedCharge2 =  1; sameFilter = 0;   break;
          case 3: part1Name = "M_"; part2Name = "M_"; requestedCharge1 = -1; requestedCharge2 = -1; sameFilter = 1;   break;
        }
      //part1Name += int(1000*etaMin);
      part1Name += "eta";
      part1Name += int(1000*etaMax);
      part1Name += "_";
      part1Name += int(1000*ptMin);
      part1Name += "pt";
      part1Name += int(1000*ptMax);
      part1Name += "_";
      //part2Name += int(1000*etaMin);
      part2Name += "eta";
      part2Name += int(1000*etaMax);
      part2Name += "_";
      part2Name += int(1000*ptMin);
      part2Name += "pt";
      part2Name += int(1000*ptMax);
      part2Name += "_";
      eventName =  "";
      eventName += int(10.*minCentrality[iCentrality] );
      eventName += "Vo";
      eventName += int(10.*maxCentrality[iCentrality] );
      //eventName += "_";
      //eventName += int(10*zMin ); 
      //eventName += "Z";
      //eventName += int(10*zMax ); 
      //if (rejectPileup)         eventName += pileupRejecSuffix;
      //if (rejectPairConversion) eventName += pairRejecSuffix;
      baseName     =   prefixName;
      baseName     +=  part1Name;
      baseName     +=  part2Name;
      baseName     +=  eventName;
      listName     =   baseName;
      taskName     =   baseName;

      //inputHistogramFileName =  "/Users/prabhat/Archive/calib/PbPb276Calibration_PP_2010.root";
      inputHistogramFileName =  "alien:///alice/cern.ch/user/p/prabhat/CalibFiles/PbPb276Calibration_2010.root";

      outputHistogramFileName = baseName;
      if (singlesOnly) outputHistogramFileName += singlesOnlySuffix;
      outputHistogramFileName += ".root";
      
    cout << "============================================================" << endl;
    cout << "                   iTask: " << iTask << endl;
      cout << "               Task Name: " << taskName << endl;
      cout << "               List Name: " << listName << endl;
      cout << "  inputHistogramFileName: " << inputHistogramFileName  << endl;
      cout << " outputHistogramFileName: " << outputHistogramFileName << endl;
      cout << "                  system: " << system << endl;
      cout << "             singlesOnly: " << singlesOnly << endl;
      cout << "           using weights: " << useWeights << endl;
      cout << "        centralityMethod: " << centralityMethod << endl;
      cout << "               chargeSet: " << chargeSet    << endl;
      cout << "              debugLevel: " << debugLevel   << endl;
      cout << "            rejectPileup: " << rejectPileup << endl;
      cout << "    rejectPairConversion: " << rejectPairConversion  << endl;
      cout << "              sameFilter: " << sameFilter  << endl;
      cout << "                    zMin: " << zMin        << endl;
      cout << "                    zMax: " << zMax        << endl;
      cout << "                   ptMin: " << ptMin       << endl;
      cout << "                   ptMax: " << ptMax       << endl;
      cout << "                  etaMin: " << etaMin   << endl;
      cout << "                  etaMax: " << etaMax   << endl;
      cout << "                 dcaZMin: " << dcaZMin  << endl;
      cout << "                 dcaZMax: " << dcaZMax  << endl;
      cout << "                dcaXYMin: " << dcaXYMin << endl;
      cout << "                dcaXYMax: " << dcaXYMax << endl;
      cout << "                 dedxMin: " << dedxMin  << endl;
      cout << "                 dedxMax: " << dedxMax  << endl;
      cout << "             nClusterMin: " << nClusterMin      << endl;
      cout << "          trackFilterBit: " << trackFilterBit   << endl;
      cout << "        requestedCharge1: " << requestedCharge1 << endl;
      cout << "        requestedCharge2: " << requestedCharge2 << endl;
    cout << "============================================================" << endl;
    
    TFile  * inputFile  = 0;
      TList  * histoList  = 0;
      TH3F   * weight_1   = 0;
      TH3F   * weight_2   = 0;
      if (useWeights)
        {
        TGrid::Connect("alien:");
        inputFile = TFile::Open(inputHistogramFileName,"OLD");
        if (!inputFile)
          {
          cout << "Requested file:" << inputHistogramFileName << " was not opened. ABORT." << endl;
          return;
          }
        TString nameHistoBase = "correction_";
        TString nameHisto;
        nameHistoBase += eventName;
        if (requestedCharge1 == 1)
          {
          nameHisto = nameHistoBase + "_p";
          cout << "Input Histogram named: " << nameHisto << endl;
          weight_1 = (TH3F *) inputFile->Get(nameHisto);
          }
        else
          {
          nameHisto = nameHistoBase + "_m";
          cout << "Input Histogram named: " << nameHisto << endl;
          weight_1 = (TH3F *) inputFile->Get(nameHisto);
          }
        if (!weight_1) 
          {
          cout << "Requested histogram 'correction_p/m' was not found. ABORT." << endl;
          return 0;
          }
        
        if (!sameFilter)
          {
          weight_2 = 0;
          if (requestedCharge2 == 1)
            {
            nameHisto = nameHistoBase + "_p";
            cout << "Input Histogram named: " << nameHisto << endl;
            weight_2 = (TH3F *) inputFile->Get(nameHisto);
            }
          else
            {
            nameHisto = nameHistoBase + "_m";
            cout << "Input Histogram named: " << nameHisto << endl;
            weight_2 = (TH3F *) inputFile->Get(nameHisto);
            }
          if (!weight_2) 
            {
            cout << "Requested histogram 'correction_p/m' was not found. ABORT." << endl;
            return 0;
            }
          }  
        }
      task = new  AliAnalysisTaskDptDptCorrelations(taskName);
      //configure my task
      task->SetDebugLevel(          debugLevel      ); 
      task->SetSameFilter(          sameFilter      );
      task->SetSinglesOnly(         singlesOnly     ); 
      task->SetUseWeights(          useWeights      ); 
      task->SetRejectPileup(        rejectPileup    ); 
      task->SetRejectPairConversion(rejectPairConversion); 
      task->SetVertexZMin(          zMin            ); 
      task->SetVertexZMax(          zMax            ); 
      task->SetVertexXYMin(         -1.            ); 
      task->SetVertexXYMax(          1.            ); 
      task->SetCentralityMethod(    centralityMethod);
      task->SetCentrality(          minCentrality[iCentrality], maxCentrality[iCentrality]);
      task->SetPtMin1(              ptMin           ); 
      task->SetPtMax1(              ptMax           ); 
      task->SetEtaMin1(             etaMin          ); 
      task->SetEtaMax1(             etaMax          ); 
      task->SetPtMin2(              ptMin           ); 
      task->SetPtMax2(              ptMax           ); 
      task->SetEtaMin2(             etaMin          ); 
      task->SetEtaMax2(             etaMax          ); 
      task->SetDcaZMin(             dcaZMin         ); 
      task->SetDcaZMax(             dcaZMax         ); 
      task->SetDcaXYMin(            dcaXYMin        ); 
      task->SetDcaXYMax(            dcaXYMax        ); 
      task->SetDedxMin(             dedxMin         ); 
      task->SetDedxMax(             dedxMax         ); 
      task->SetNClusterMin(         nClusterMin     ); 
      task->SetTrackFilterBit(      trackFilterBit  );
      task->SetRequestedCharge_1(   requestedCharge1); 
      task->SetRequestedCharge_2(   requestedCharge2); 
      task->SetWeigth_1(            weight_1        );
      task->SetWeigth_2(            weight_2        );
      
      
      cout << "Creating task output container" << endl;
      taskOutputContainer = analysisManager->CreateContainer(listName, 
                                                             TList::Class(),    
                                                             AliAnalysisManager::kOutputContainer, 
                                                             Form("%s:Histos", AliAnalysisManager::GetCommonFileName()));
      cout << "Add task to analysis manager and connect it to input and output containers" << endl;
      analysisManager->AddTask(task);
      analysisManager->ConnectInput( task,  0, analysisManager->GetCommonInputContainer());
      analysisManager->ConnectOutput(task,  0, taskOutputContainer );
      cout << "Task added ...." << endl;
      
      iTask++;
    
    }
  
  
  
  return task;
}
