//
// Macro designed for use with the AliAnalysisTaskDptDptCorrelations task.
//
// Author: Claude Pruneau, Wayne State
/////////////////////////////////////////////////////////////////////////////////
AliAnalysisTaskDptDptCorrelations *AddTaskDptDptCorrelations(int    singlesOnly            = 0,
                                                             int    useWeights             = 1,
                                                             int    centralityMethod       = 4,
                                                             int    chargeSet              = 1)

{
  // Set Default Configuration of this analysis
  // ==========================================
  int    debugLevel             = 0;
  //int    singlesOnly            = 1;
  //int    useWeights             = 0;
  int    rejectPileup           = 1;
  int    rejectPairConversion   = 1;
  int    sameFilter             = 1;
  
  //int    centralityMethod       =  4;
  int    nCentrality            =  10;
  double minCentrality[] = { 0.5,   5., 10., 20., 30., 40., 50., 60., 70., 80. };
  double maxCentrality[] = { 5.0,  10., 20., 30., 40., 50., 60., 70., 80., 90. }; 
  
  double zMin                   = -10.;
  double zMax                   =  10.;
  double ptMin                  =  0.2;
  double ptMax                  =  2.0;
  double etaMin                 = -1.0;
  double etaMax                 =  1.0;
  double dcaZMin                = -3.0;
  double dcaZMax                =  3.0;
  double dcaXYMin               = -3.0;
  double dcaXYMax               =  3.0;
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
      //inputHistogramFileName = inputPath;
      //inputHistogramFileName += "/";
      //inputHistogramFileName =  baseName;
      //inputHistogramFileName += calibSuffix;
      //inputHistogramFileName += ".root";
      //outputHistogramFileName = outputPath;
      //outputHistogramFileName += "/";
      //inputHistogramFileName =  "/home/pruneau//wrk/Alice/PbPb/2730GeV/DptDpt/Calib/PbPb273Calibration.root";
    inputHistogramFileName =  "alien:///alice/cern.ch/user/c/cpruneau/PbPb273Calibration.root"; //TFile::Open();
    outputHistogramFileName = baseName;
      if (singlesOnly) outputHistogramFileName += singlesOnlySuffix;
      outputHistogramFileName += ".root";
      
      cout << "                   iTask: " << iTask << endl;
      cout << "               Task Name: " << taskName << endl;
      cout << "               List Name: " << listName << endl;
      cout << "  inputHistogramFileName: " << inputHistogramFileName  << endl;
      cout << " outputHistogramFileName: " << outputHistogramFileName << endl;
      cout << "           using weights: " << useWeights              << endl;
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
