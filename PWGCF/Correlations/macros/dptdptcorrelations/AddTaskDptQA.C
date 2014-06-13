//           system:  0: PbPb                 1: pp
//      singlesOnly:  0: full correlations    1: singles only
//       useWeights:  0: no                   1: yes
// centralityMethod:  3: track count  4: V0 centrality
/////////////////////////////////////////////////////////////////////////////////
AliAnalysisTaskDptDptQA *AddTaskDptQA
(int    system                 = 0,
 int    singlesOnly            = 0,
 int    useWeights             = 0,
 int    centralityMethod       = 4,
 int    centralitySelected     = 4,
 double etaMin                 = -0.8,
 double etaMax                 =  0.8,
 int    trackFilterBit         =  128,
 char *inputHistogramFileName  = "alien:///alice/cern.ch/user/p/prabhat/CalibFiles/PbPbCalib_dca1.root"
 )
  
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
	  
	  return 0;
	}
    }
  
  double zMin                   = -10.;
  double zMax                   =  10.;
  double ptMin                  =  0.2;
  double ptMax                  =  2.0;
  
  double dcaZMin                = -3.0;
  double dcaZMax                =  3.0;
  double dcaXYMin               = -3.0;
  double dcaXYMax               =  3.0;
  double dedxMin                =  0.0;
  double dedxMax                =  20000.0;
  int    nClusterMin            =   70;
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
  TString outputHistogramFileName;
  
  // Create the task and add subtask.
  // ===========================================================================
  int iTask = 0; // task counter
  AliAnalysisDataContainer *taskInputContainer;
  AliAnalysisDataContainer *taskOutputContainer;
  AliAnalysisTaskDptDptQA* task;

  TFile  * inputFile  = 0;
  TList  * histoList  = 0;
  TH3F   * weight_1   = 0;
  TH3F   * weight_2   = 0;
  
  int iCentrality = centralitySelected;
  outputHistogramFileName = baseName;
  if (singlesOnly) outputHistogramFileName += singlesOnlySuffix;
  outputHistogramFileName += ".root";

  if (useWeights)
  {
    TGrid::Connect("alien:");
    inputFile = TFile::Open(inputHistogramFileName,"OLD");
    if (!inputFile)
      {
      cout << "Requested file:" << inputHistogramFileName << " was not opened. ABORT." << endl;
      return;
      }
    }
  
  //============================
  // (+,+)
  //============================
  requestedCharge1 = 1;
  requestedCharge2 = 1;
  sameFilter = 1;
  part1Name = "P_";
  part1Name += trackFilterBit;
  part1Name += "_";
  part1Name += int(1000*etaMax);
  part1Name += "_";
  part1Name += int(1000*ptMin);
  part1Name += "pt";
  part1Name += int(1000*ptMax);
  part1Name += "_";
  part2Name = "P_"; 
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
  baseName  =   prefixName;
  baseName  +=  part1Name;
  baseName  +=  part2Name;
  baseName  +=  eventName;
  listName  =   baseName;
  taskName  =   baseName;
  if (useWeights)
    {
    TString nameHistoBase = "correction_";
    TString nameHisto;
    nameHistoBase += eventName;
    nameHisto = nameHistoBase + "_p";
    cout << "Input Histogram named: " << nameHisto << endl;
    weight_1 = (TH3F *) inputFile->Get(nameHisto);
    weight_2 = weight_1;
    }
  else
    {
    weight_1 = 0;
    weight_2 = 0;
    }
  task = new  AliAnalysisTaskDptDptQA(taskName);
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
  cout << "(+,+) Task added ...." << endl;
  
  //============================
  // (+,-)
  //============================
  requestedCharge1 =  1;
  requestedCharge2 = -1;
  sameFilter = 0;
  part1Name = "P_"; 
  part1Name += trackFilterBit;
  part1Name += "_";    
  part1Name += int(1000*etaMax);
  part1Name += "_";
  part1Name += int(1000*ptMin);
  part1Name += "pt";
  part1Name += int(1000*ptMax);
  part1Name += "_";
  part2Name = "M_"; 
  part2Name += trackFilterBit;
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
  baseName  =   prefixName;
  baseName  +=  part1Name;
  baseName  +=  part2Name;
  baseName  +=  eventName;
  listName  =   baseName;
  taskName  =   baseName;
  if (useWeights)
    {
    TString nameHistoBase = "correction_";
    TString nameHisto;
    nameHistoBase += eventName;
    nameHisto = nameHistoBase + "_p";
    cout << "Input Histogram named: " << nameHisto << endl;
    weight_1 = (TH3F *) inputFile->Get(nameHisto);
    nameHisto = nameHistoBase + "_m";
    cout << "Input Histogram named: " << nameHisto << endl;
    weight_2 = (TH3F *) inputFile->Get(nameHisto);
    }
  else
    {
    weight_1 = 0;
    weight_2 = 0;
    }
  task = new  AliAnalysisTaskDptDptQA(taskName);
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

  //============================
  // (-,-)
  //============================
  requestedCharge1 = -1;
  requestedCharge2 = -1;
  sameFilter = 1;
  part1Name = "M_"; 
  part1Name += trackFilterBit;
  part1Name += "_";
  part1Name += int(1000*etaMax);
  part1Name += "_";
  part1Name += int(1000*ptMin);
  part1Name += "pt";
  part1Name += int(1000*ptMax);
  part1Name += "_";
  part2Name = "M_"; 
  part2Name += trackFilterBit;
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
  baseName  =   prefixName;
  baseName  +=  part1Name;
  baseName  +=  part2Name;
  baseName  +=  eventName;
  listName  =   baseName;
  taskName  =   baseName;
  if (useWeights)
    {
    TString nameHistoBase = "correction_";
    TString nameHisto;
    nameHistoBase += eventName;
    nameHisto = nameHistoBase + "_m";
    cout << "Input Histogram named: " << nameHisto << endl;
    weight_1 = (TH3F *) inputFile->Get(nameHisto);
    weight_2 = weight_1;
    }
  else
    {
    weight_1 = 0;
    weight_2 = 0;
    }
  task = new  AliAnalysisTaskDptDptQA(taskName);
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
  analysisManager->AddTask(task);
  analysisManager->ConnectInput( task,  0, analysisManager->GetCommonInputContainer());
  analysisManager->ConnectOutput(task,  0, taskOutputContainer );
  return task;
}
