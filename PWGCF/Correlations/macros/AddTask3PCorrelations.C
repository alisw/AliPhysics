//
// Macro designed for use with the AliAnalysisTaskDptDptCorrelations task.
//
// Author: Claude Pruneau, Wayne State
// 
//           system:  0: PbPb                 1: pp
//      singlesOnly:  0: full correlations    1: singles only
//       useWeights:  0: no                   1: yes
// centralityMethod:  3: track count  4: V0 centrality
//////////////////////////////////////////////////////////////////////////////
AliAnalysisTask3PCorrelations * AddTask3PCorrelations(int    etaOption              = 0,
                                                      int    system                 = 0,
                                                      int    singlesOnly            = 0,
                                                      int    useWeights             = 0,
                                                      int    centralityMethod       = 4)

{
  // Set Default Configuration of this analysis
  // ==========================================
  int    debugLevel             = 0;
  int    rejectPileup           = 1;
  int    rejectPairConversion   = 1;
  int    nCentrality;
  double minCentrality[10];
  double maxCentrality[10];

  if (system==0) // PbPb
    {
    if (centralityMethod == 4)
      {
      nCentrality = 10;
      minCentrality[0] = 0.5; maxCentrality[0] = 5.0;
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
      cout << "-F- AddTask3PCorrelations() system:" << system << ". centralityMethod:" << centralityMethod << " Option NOT AVAILABLE. ABORT."
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
      cout << "-F- AddTask3PCorrelations() system:" << system << ". centralityMethod:" << centralityMethod << " Option NOT AVAILABLE. ABORT."
      return 0;
      }
    }
  else
    {
    cout << "-F- AddTask3PCorrelations() system:" << system << ". Option NOT CURRENTLY AVAILABLE. ABORT."
    return 0;
    }

  double zMin                   = -10.;
  double zMax                   =  10.;
  double dcaZMin                = -3.0;
  double dcaZMax                =  3.0;
  double dcaXYMin               = -3.0;
  double dcaXYMax               =  3.0;
  double dedxMin                =  0.0;
  double dedxMax                =  20000.0;
  int    nClusterMin            =   70;
  int    trackFilterBit         =  128;
  double ptMin1                 =   2.0;
  double ptMax1                 =  10.0;
  double etaMin1                = -0.25;
  double etaMax1                =  0.25;
  double ptMin2                 =  1.0;
  double ptMax2                 =  2.0;
  double etaMin2                = -1.0;
  double etaMax2                = -0.7;
  double ptMin3                 =  1.0;
  double ptMax3                 =  2.0;
  double etaMin3                =  0.7;
  double etaMax3                =  1.0;
  
  if (etaOption == 0)
    {
    etaMin2 = -1.0;
    etaMax2 = -0.7;
    etaMin3 =  0.7;
    etaMax3 =  1.0;
    }
  else
    {
    etaMin2 = etaMin1;
    etaMax2 = etaMax1;
    etaMin3 = etaMin1;
    etaMax3 = etaMax1;
    
    }
  
  // Get the pointer to the existing analysis manager via the static access method.
  // ==============================================================================
  AliAnalysisManager *analysisManager = AliAnalysisManager::GetAnalysisManager();
  
  if (!analysisManager) 
    {
    ::Error("AddTask3PCorrelations", "No analysis manager to connect to.");
    return NULL;
    }  
  
  TString part1Name;
  TString part2Name;
  TString part3Name;
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
  AliAnalysisTask3PCorrelations * task;
  
  for (int iCentrality=0; iCentrality < nCentrality; ++iCentrality)
    {
    part1Name = "";
    part1Name += int(1000*ptMin1);
      part1Name += "pt";
      part1Name += int(1000*ptMax1);
      part1Name += "_";

    part2Name = "";
      part2Name += int(1000*ptMin2);
      part2Name += "pt";
      part2Name += int(1000*ptMax2);
      part2Name += "_";

    part3Name = "";
    part3Name += int(1000*ptMin3);
    part3Name += "pt";
    part3Name += int(1000*ptMax3);
    part3Name += "_";

    eventName =  "";
      eventName += int(10.*minCentrality[iCentrality] );
      eventName += "Vo";
      eventName += int(10.*maxCentrality[iCentrality] );
      //if (rejectPileup)         eventName += pileupRejecSuffix;
      //if (rejectPairConversion) eventName += pairRejecSuffix;
      baseName     =   prefixName;
      baseName     +=  part1Name;
    baseName     +=  part2Name;
    baseName     +=  part3Name;
      baseName     +=  eventName;
      listName     =   baseName;
      taskName     =   baseName;
    inputHistogramFileName =  "alien:///alice/cern.ch/user/c/cpruneau/PbPb273Calibration.root"; //TFile::Open();
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
      cout << "              debugLevel: " << debugLevel   << endl;
      cout << "            rejectPileup: " << rejectPileup << endl;
      cout << "    rejectPairConversion: " << rejectPairConversion  << endl;
      cout << "                    zMin: " << zMin        << endl;
      cout << "                    zMax: " << zMax        << endl;
    
    cout << "                   ptMin1: " << ptMin1       << endl;
    cout << "                   ptMax1: " << ptMax1       << endl;
    cout << "                  etaMin1: " << etaMin1   << endl;
    cout << "                  etaMax1: " << etaMax1   << endl;
      
    cout << "                   ptMin2: " << ptMin2       << endl;
    cout << "                   ptMax2: " << ptMax2       << endl;
    cout << "                  etaMin2: " << etaMin2   << endl;
    cout << "                  etaMax2: " << etaMax2   << endl;

    cout << "                   ptMin3: " << ptMin3       << endl;
    cout << "                   ptMax3: " << ptMax3       << endl;
    cout << "                  etaMin3: " << etaMin3   << endl;
    cout << "                  etaMax3: " << etaMax3   << endl;

    
    cout << "                 dcaZMin: " << dcaZMin  << endl;
      cout << "                 dcaZMax: " << dcaZMax  << endl;
      cout << "                dcaXYMin: " << dcaXYMin << endl;
      cout << "                dcaXYMax: " << dcaXYMax << endl;
      cout << "                 dedxMin: " << dedxMin  << endl;
      cout << "                 dedxMax: " << dedxMax  << endl;
      cout << "             nClusterMin: " << nClusterMin      << endl;
      cout << "          trackFilterBit: " << trackFilterBit   << endl;
    cout << "============================================================" << endl;
    
    TFile  * inputFile  = 0;
      TList  * histoList  = 0;
    TH3F   * weight_1P   = 0;
    TH3F   * weight_1M   = 0;
    TH3F   * weight_2P   = 0;
    TH3F   * weight_2M   = 0;
    TH3F   * weight_3P   = 0;
    TH3F   * weight_3M   = 0;
    TString nameHistoBase = "correction_";
    TString nameHisto;
      if (useWeights)
        {
        TGrid::Connect("alien:");
        inputFile = TFile::Open(inputHistogramFileName,"OLD");
        if (!inputFile)
          {
          cout << "Requested file:" << inputHistogramFileName << " was not opened. ABORT." << endl;
          return;
          }
        nameHistoBase += eventName;
        nameHisto = nameHistoBase + "_1p"; cout << "Input Histogram named: " << nameHisto << endl;
        weight_1P = (TH3F *) inputFile->Get(nameHisto);
        nameHisto = nameHistoBase + "_1m"; cout << "Input Histogram named: " << nameHisto << endl;
        weight_1M = (TH3F *) inputFile->Get(nameHisto);
        nameHisto = nameHistoBase + "_2p"; cout << "Input Histogram named: " << nameHisto << endl;
        weight_2P = (TH3F *) inputFile->Get(nameHisto);
        nameHisto = nameHistoBase + "_2m"; cout << "Input Histogram named: " << nameHisto << endl;
        weight_2M = (TH3F *) inputFile->Get(nameHisto);
        nameHisto = nameHistoBase + "_3p"; cout << "Input Histogram named: " << nameHisto << endl;
        weight_3P = (TH3F *) inputFile->Get(nameHisto);
        nameHisto = nameHistoBase + "_3m"; cout << "Input Histogram named: " << nameHisto << endl;
        weight_3M = (TH3F *) inputFile->Get(nameHisto);

        }
      task = new  AliAnalysisTask3PCorrelations(taskName);
      //configure my task
      task->SetDebugLevel(          debugLevel      ); 
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
      task->SetPtMin1(              ptMin1           ); 
      task->SetPtMax1(              ptMax1           ); 
      task->SetEtaMin1(             etaMin1          ); 
      task->SetEtaMax1(             etaMax1          ); 
    task->SetPtMin2(              ptMin2           ); 
    task->SetPtMax2(              ptMax2           ); 
    task->SetEtaMin2(             etaMin2          ); 
    task->SetEtaMax2(             etaMax2          ); 
    task->SetPtMin3(              ptMin3           ); 
    task->SetPtMax3(              ptMax3           ); 
    task->SetEtaMin3(             etaMin3          ); 
    task->SetEtaMax3(             etaMax3          ); 
      task->SetDcaZMin(             dcaZMin         ); 
      task->SetDcaZMax(             dcaZMax         ); 
      task->SetDcaXYMin(            dcaXYMin        ); 
      task->SetDcaXYMax(            dcaXYMax        ); 
      task->SetDedxMin(             dedxMin         ); 
      task->SetDedxMax(             dedxMax         ); 
      task->SetNClusterMin(         nClusterMin     ); 
      task->SetTrackFilterBit(      trackFilterBit  );
    task->SetWeigth_1P(           weight_1P       );
    task->SetWeigth_1M(           weight_1M       );
    task->SetWeigth_2P(           weight_2P       );
    task->SetWeigth_2M(           weight_2M       );
    task->SetWeigth_3P(           weight_3P       );
    task->SetWeigth_3M(           weight_3M       );

      
      
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
