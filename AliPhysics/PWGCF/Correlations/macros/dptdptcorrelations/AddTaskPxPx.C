// Macro designed for use with the AliAnalysisTaskDptDptCorrelations task.
// Author: Prabhat Pujahari & Claude Pruneau, Wayne State
//           system:  0: PbPb                 1: pPb
//      singlesOnly:  0: full correlations    1: singles only
//       useWeights:  0: no                   1: yes
// centralityMethod:  3: track count  4: V0 centrality  7: V0A centrality for pPb
//        chargeSet:  0: ++    1: +-    2: -+    3: --
/////////////////////////////////////////////////////////////////////////////////
AliAnalysisTaskpxpx *AddTaskPxPx
(int    system                  = 0,
 int    singlesOnly             = 0,
 int    useWeights              = 1, 
 int    centralityMethod        = 4,
 int    chargeSet               = 1,
 double zMin                   = -10.,
 double zMax                   =  10.,
 int    trackFilterBit          = 128,
 int    nClusterMin             = 80, 
 double eta1Min                 = -0.8,
 double eta1Max                 = 0.8,
 double eta2Min                 = -0.8,
 double eta2Max                 = 0.8,
 double dcaZMin                 = -3.2,
 double dcaZMax                 =  3.2,
 double dcaXYMin                = -2.4,
 double dcaXYMax                =  2.4,
 int nCentrality                =  1,
 Bool_t trigger                 = kFALSE,
 const char* taskname           = "dcaz2",
 char *inputHistogramFileName   = "alien:///alice/cern.ch/user/p/prabhat/CalibFiles/PbPbCalib_dca1.root")
  
{
  // Set Default Configuration of this analysis
  // ==========================================
  int    debugLevel             = 0;
  int    rejectPileup           = 1;
  int    rejectPairConversion   = 1;
  int    sameFilter             = 1;
  
  //int    nCentrality;
  double minCentrality[10];
  double maxCentrality[10];

  if (system==0) // PbPb
    {
    if (centralityMethod == 4)
      {
	minCentrality[0] = 0.0;  maxCentrality[0]  = 5.0;
        minCentrality[1] = 30.;  maxCentrality[1]  = 40.;
        minCentrality[2] = 60.;  maxCentrality[2]  = 70.;
        
      }
    else
      {

      return 0;
      }
    }
  else if (system==1) // pPb
    {
    if (centralityMethod == 7)
      {
	minCentrality[0] = 0;   maxCentrality[0] = 20.0;
	minCentrality[1] = 20.; maxCentrality[1] = 40.;
	minCentrality[2] = 40.; maxCentrality[2] = 60.;
	minCentrality[3] = 60.; maxCentrality[3] = 80.;
      }
    else
      {
      return 0;
      }
    }
  else
    {
    return 0;
    }

  double ptMin                  =  0.2;
  double ptMax                  =  2.0;
  double dedxMin                =  0.0;
  double dedxMax                =  20000.0;
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
  //TString inputHistogramFileName;
  TString outputHistogramFileName;
  
  // Create the task and add subtask.
  // ===========================================================================
  int iTask = 0; // task counter
  AliAnalysisDataContainer *taskInputContainer;
  AliAnalysisDataContainer *taskOutputContainer;
  AliAnalysisTaskpxpx* task;
  
  for (int iCentrality=0; iCentrality < nCentrality; ++iCentrality)
    {
      switch (chargeSet)
        {
          case 0: part1Name = "P_"; part2Name = "P_"; requestedCharge1 =  1; requestedCharge2 =  1; sameFilter = 1; break;
          case 1: part1Name = "P_"; part2Name = "M_"; requestedCharge1 =  1; requestedCharge2 = -1; sameFilter = 0;   break;
          case 2: part1Name = "M_"; part2Name = "P_"; requestedCharge1 = -1; requestedCharge2 =  1; sameFilter = 0;   break;
          case 3: part1Name = "M_"; part2Name = "M_"; requestedCharge1 = -1; requestedCharge2 = -1; sameFilter = 1;   break;
        }

      part1Name += "eta";
      part1Name += int(1000*eta1Max);
      part1Name += "_";
      part1Name += int(1000*ptMin);
      part1Name += "pt";
      part1Name += int(1000*ptMax);
      part1Name += "_";
      part1Name += int(1000*dcaZMin);
      part1Name += "DCA";
      part1Name += int(1000*dcaZMax);
      part1Name += "_";

      part2Name += "eta";
      part2Name += int(1000*eta2Max);
      part2Name += "_";
      part2Name += int(1000*ptMin);
      part2Name += "pt";
      part2Name += int(1000*ptMax);
      part2Name += "_";
      part2Name += int(1000*dcaZMin);
      part2Name += "DCA";
      part2Name += int(1000*dcaZMax);
      part2Name += "_";

      eventName =  "";
      eventName += int(10.*minCentrality[iCentrality] );
      eventName += "Vo";
      eventName += int(10.*maxCentrality[iCentrality] );

      baseName     =   prefixName;
      baseName     +=  part1Name;
      baseName     +=  part2Name;
      baseName     +=  eventName;
      listName     =   baseName;
      taskName     =   baseName;


      outputHistogramFileName = baseName;
      if (singlesOnly) outputHistogramFileName += singlesOnlySuffix;
      outputHistogramFileName += ".root";
      
     
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
	    //cout << "Requested file:" << inputHistogramFileName << " was not opened. ABORT." << endl;
          return;
          }
        TString nameHistoBase = "correction_";
        TString nameHisto;
        nameHistoBase += eventName;
        if (requestedCharge1 == 1)
          {
          nameHisto = nameHistoBase + "_p";
          //cout << "Input Histogram named: " << nameHisto << endl;
          weight_1 = (TH3F *) inputFile->Get(nameHisto);
          }
        else
          {
          nameHisto = nameHistoBase + "_m";
          //cout << "Input Histogram named: " << nameHisto << endl;
          weight_1 = (TH3F *) inputFile->Get(nameHisto);
          }
        if (!weight_1) 
          {
	    //cout << "Requested histogram 'correction_p/m' was not found. ABORT." << endl;
          return 0;
          }
        
        if (!sameFilter)
          {
          weight_2 = 0;
          if (requestedCharge2 == 1)
            {
            nameHisto = nameHistoBase + "_p";
            //cout << "Input Histogram named: " << nameHisto << endl;
            weight_2 = (TH3F *) inputFile->Get(nameHisto);
            }
          else
            {
            nameHisto = nameHistoBase + "_m";
            //cout << "Input Histogram named: " << nameHisto << endl;
            weight_2 = (TH3F *) inputFile->Get(nameHisto);
            }
          if (!weight_2) 
            {
	      //cout << "Requested histogram 'correction_p/m' was not found. ABORT." << endl;
            return 0;
            }
          }  
        }
      task = new  AliAnalysisTaskpxpx(taskName);
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
      task->SetEtaMin1(             eta1Min          ); 
      task->SetEtaMax1(             eta1Max          ); 
      task->SetPtMin2(              ptMin           ); 
      task->SetPtMax2(              ptMax           ); 
      task->SetEtaMin2(             eta2Min          ); 
      task->SetEtaMax2(             eta2Max          ); 
      task->SetDcaZMin(             dcaZMin         ); 
      task->SetDcaZMax(             dcaZMax         ); 
      task->SetDcaXYMin(            dcaXYMin        ); 
      task->SetDcaXYMax(            dcaXYMax        ); //checking by prp
      task->SetDedxMin(             dedxMin         ); 
      task->SetDedxMax(             dedxMax         ); 
      task->SetNClusterMin(         nClusterMin     ); 
      task->SetTrackFilterBit(      trackFilterBit  );
      task->SetRequestedCharge_1(   requestedCharge1); 
      task->SetRequestedCharge_2(   requestedCharge2); 
      task->SetWeigth_1(            weight_1        );
      task->SetWeigth_2(            weight_2        );
      

      if(trigger) task->SelectCollisionCandidates(AliVEvent::kINT7);
      else task->SelectCollisionCandidates(AliVEvent::kMB);

      cout << "Creating task output container" << endl;

      taskOutputContainer = analysisManager->CreateContainer(listName,
                                                             TList::Class(),
                                                             AliAnalysisManager::kOutputContainer,
                                                             Form("%s:%s", AliAnalysisManager::GetCommonFileName(),taskname));
      cout << "Add task to analysis manager and connect it to input and output containers" << endl;

      analysisManager->AddTask(task);
      analysisManager->ConnectInput( task,  0, analysisManager->GetCommonInputContainer());
      analysisManager->ConnectOutput(task,  0, taskOutputContainer );
      //cout << "Task added ...." << endl;
      
      iTask++;
    
    }
  
  
  
  return task;
}
