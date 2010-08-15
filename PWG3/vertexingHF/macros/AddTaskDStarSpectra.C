//DEFINITION OF A FEW CONSTANTS
const Int_t    numberOfSigmasPID  = 3;
// ANALYSIS TYPE 
const Bool_t anaType   = 1;//0 HD; 1 UU;
const Bool_t usePID = kTRUE;
//----------------------------------------------------

AliAnalysisTaskSEDStarSpectra *AddTaskDStarSpectra(Bool_t theMCon=kFALSE)
{
  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskDStarSpectra", "No analysis manager to connect to.");
    return NULL;
  }  
  
  TString cutobjname="Dstar";
  
  // D0 daughters pre-selections
  //
  AliRDHFCutsDStartoKpipi* RDHFDStartoKpipi=new AliRDHFCutsDStartoKpipi(); 
  RDHFDStartoKpipi->SetName("DstartoKpipi");
  RDHFDStartoKpipi->SetTitle("cuts for D* analysis");

  AliESDtrackCuts* esdTrackCuts=new AliESDtrackCuts();
  esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
  //default
  esdTrackCuts->SetRequireTPCRefit(kTRUE);
  esdTrackCuts->SetRequireITSRefit(kTRUE);
  esdTrackCuts->SetMinNClustersITS(4); // default is 5
  esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
  					 AliESDtrackCuts::kAny); //test d0 asimmetry
                                                                 // default is kBoth, otherwise kAny
  esdTrackCuts->SetMinDCAToVertexXY(0.);
  esdTrackCuts->SetPtRange(0.1,1.e10);
  
  //
  // soft pion pre-selections 
  // 
  AliESDtrackCuts* esdSoftPicuts=new AliESDtrackCuts();
  esdSoftPicuts->SetRequireSigmaToVertex(kFALSE);
  //default
  esdSoftPicuts->SetRequireTPCRefit(kFALSE);
  esdSoftPicuts->SetRequireITSRefit(kFALSE);
  esdSoftPicuts->SetMinNClustersITS(4); // default is 4
  esdSoftPicuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
					  AliESDtrackCuts::kAny); //test d0 asimmetry
  esdSoftPicuts->SetPtRange(0.0,1.e10);

  // set pre selections
  RDHFDStartoKpipi->AddTrackCuts(esdTrackCuts);
  RDHFDStartoKpipi->AddTrackCutsSoftPi(esdSoftPicuts);
  
  const Int_t nvars=14;
  
  const Int_t nptbins=7;
  Float_t* ptbins;
  ptbins=new Float_t[nptbins+1];
  ptbins[0]=0.7;
  ptbins[1]=1.;
  ptbins[2]=2.;
  ptbins[3]=3.;
  ptbins[4]=5.;
  ptbins[5]=8.; 
  ptbins[6]=12.;
  ptbins[7]=18.;
  
  RDHFDStartoKpipi->SetPtBins(nptbins+1,ptbins);
  

  Float_t** rdcutsvalmine;
  rdcutsvalmine=new Float_t*[nvars];
  for(Int_t iv=0;iv<nvars;iv++){
    rdcutsvalmine[iv]=new Float_t[nptbins];
  }

  if(anaType==1){ // UU cuts

    rdcutsvalmine[0][0]=0.7;
    rdcutsvalmine[1][0]=0.022;
    rdcutsvalmine[2][0]=0.7;
    rdcutsvalmine[3][0]=0.21;
    rdcutsvalmine[4][0]=0.21;
    rdcutsvalmine[5][0]=0.05;
    rdcutsvalmine[6][0]=0.05;
    rdcutsvalmine[7][0]=-0.00005;
    rdcutsvalmine[8][0]=0.85;
    rdcutsvalmine[9][0]=0.3.;
    rdcutsvalmine[10][0]=0.1;
    rdcutsvalmine[11][0]=0.05;
    rdcutsvalmine[12][0]=100.;
    rdcutsvalmine[13][0]=0.5;

    rdcutsvalmine[0][1]=0.7;
    rdcutsvalmine[1][1]=0.03;
    rdcutsvalmine[2][1]=0.8;
    rdcutsvalmine[3][1]=0.45;
    rdcutsvalmine[4][1]=0.45;
    rdcutsvalmine[5][1]=0.09;
    rdcutsvalmine[6][1]=0.09;
    rdcutsvalmine[7][1]=-0.00029;
    rdcutsvalmine[8][1]=0.8;
    rdcutsvalmine[9][1]=0.3;
    rdcutsvalmine[10][1]=0.1;
    rdcutsvalmine[11][1]=0.05;
    rdcutsvalmine[12][1]=100.;
    rdcutsvalmine[13][1]=1.0;
    
    rdcutsvalmine[0][2]=0.7;
    rdcutsvalmine[1][2]=0.02;
    rdcutsvalmine[2][2]=0.8;
    rdcutsvalmine[3][2]=0.7;
    rdcutsvalmine[4][2]=0.7;
    rdcutsvalmine[5][2]=0.08;
    rdcutsvalmine[6][2]=0.08;
    rdcutsvalmine[7][2]=-0.00018;
    rdcutsvalmine[8][2]=0.90;
    rdcutsvalmine[9][2]=0.3;
    rdcutsvalmine[10][2]=0.1;
    rdcutsvalmine[11][2]=0.05;
    rdcutsvalmine[12][2]=100.;
    rdcutsvalmine[13][2]=0.5;
    
    rdcutsvalmine[0][3]=0.7;
    rdcutsvalmine[1][3]=0.05;
    rdcutsvalmine[2][3]=0.8;
    rdcutsvalmine[3][3]=1.;
    rdcutsvalmine[4][3]=1.;
    rdcutsvalmine[5][3]=0.042;
    rdcutsvalmine[6][3]=0.056;
    rdcutsvalmine[7][3]=-0.000065;
    rdcutsvalmine[8][3]=0.9;
    rdcutsvalmine[9][3]=0.3;
    rdcutsvalmine[10][3]=0.1;
    rdcutsvalmine[11][3]=0.05;
    rdcutsvalmine[12][3]=100.;
    rdcutsvalmine[13][3]=0.5;
    
    rdcutsvalmine[0][4]=0.7;
    rdcutsvalmine[1][4]=0.08;
    rdcutsvalmine[2][4]=0.9;
    rdcutsvalmine[3][4]=1.2;
    rdcutsvalmine[4][4]=1.2;
    rdcutsvalmine[5][4]=0.07;
    rdcutsvalmine[6][4]=0.07;
    rdcutsvalmine[7][4]=0.0001;
    rdcutsvalmine[8][4]=0.9;
    rdcutsvalmine[9][4]=0.3;
    rdcutsvalmine[10][4]=0.1;
    rdcutsvalmine[11][4]=0.05;
    rdcutsvalmine[12][4]=100.;
    rdcutsvalmine[13][4]=0.5;
 

    rdcutsvalmine[0][5]=0.7;
    rdcutsvalmine[1][5]=0.1;
    rdcutsvalmine[2][5]=1.0;
    rdcutsvalmine[3][5]=1.;
    rdcutsvalmine[4][5]=1.;
    rdcutsvalmine[5][5]=0.08;
    rdcutsvalmine[6][5]=0.08;
    rdcutsvalmine[7][5]=0.0004;
    rdcutsvalmine[8][5]=0.9;
    rdcutsvalmine[9][5]=0.3;
    rdcutsvalmine[10][5]=0.1;
    rdcutsvalmine[11][5]=0.05;
    rdcutsvalmine[12][5]=100000.;
    rdcutsvalmine[13][5]=0.5;

    rdcutsvalmine[0][6]=0.7;
    rdcutsvalmine[1][6]=0.1;
    rdcutsvalmine[2][6]=1.0;
    rdcutsvalmine[3][6]=1.;
    rdcutsvalmine[4][6]=1.;
    rdcutsvalmine[5][6]=0.1;
    rdcutsvalmine[6][6]=0.1;
    rdcutsvalmine[7][6]=0.0005;
    rdcutsvalmine[8][6]=0.9;
    rdcutsvalmine[9][6]=0.3;
    rdcutsvalmine[10][6]=0.1;
    rdcutsvalmine[11][6]=0.05;
    rdcutsvalmine[12][6]=100.;
    rdcutsvalmine[13][6]=0.5;

    rdcutsvalmine[0][7]=0.7;
    rdcutsvalmine[1][7]=0.1;
    rdcutsvalmine[2][7]=1.0;
    rdcutsvalmine[3][7]=1.;
    rdcutsvalmine[4][7]=1.;
    rdcutsvalmine[5][7]=0.1;
    rdcutsvalmine[6][7]=0.1;
    rdcutsvalmine[7][7]=0.001;
    rdcutsvalmine[8][7]=0.9;
    rdcutsvalmine[9][7]=0.3;
    rdcutsvalmine[10][7]=0.1;
    rdcutsvalmine[11][7]=0.05;
    rdcutsvalmine[12][7]=100.;
    rdcutsvalmine[13][7]=0.5;

   
  }else{ // HD cuts - To be setted
    
    rdcutsvalmine[0][0]=0.450;
    rdcutsvalmine[1][0]=0.02;
    rdcutsvalmine[2][0]=0.7;
    rdcutsvalmine[3][0]=0.8;
    rdcutsvalmine[4][0]=0.8;
    rdcutsvalmine[5][0]=0.1;
    rdcutsvalmine[6][0]=0.1;
    rdcutsvalmine[7][0]=0.000002;
    rdcutsvalmine[8][0]=0.9;
    rdcutsvalmine[9][0]=200.;
    rdcutsvalmine[10][0]=200.;
    rdcutsvalmine[11][0]=0.021;
    rdcutsvalmine[12][0]=100.;
    rdcutsvalmine[13][0]=0.5;

    rdcutsvalmine[0][1]=0.45;
    rdcutsvalmine[1][1]=0.03;
    rdcutsvalmine[2][1]=0.7;
    rdcutsvalmine[3][1]=0.8;
    rdcutsvalmine[4][1]=0.8;
    rdcutsvalmine[5][1]=0.1;
    rdcutsvalmine[6][1]=0.1;
    rdcutsvalmine[7][1]=-0.00002;
    rdcutsvalmine[8][1]=0.9;
    rdcutsvalmine[9][1]=200.;
    rdcutsvalmine[10][1]=200.;
    rdcutsvalmine[11][1]=0.021;
    rdcutsvalmine[12][1]=100.;
    rdcutsvalmine[13][1]=0.5;
    
    rdcutsvalmine[0][2]=0.450;
    rdcutsvalmine[1][2]=0.03;
    rdcutsvalmine[2][2]=0.7;
    rdcutsvalmine[3][2]=0.8;
    rdcutsvalmine[4][2]=0.8;
    rdcutsvalmine[5][2]=0.1;
    rdcutsvalmine[6][2]=0.1;
    rdcutsvalmine[7][2]=-0.00002;
    rdcutsvalmine[8][2]=0.9;
    rdcutsvalmine[9][2]=200.;
    rdcutsvalmine[10][2]=200.;
    rdcutsvalmine[11][2]=0.021;
    rdcutsvalmine[12][2]=100.;
    rdcutsvalmine[13][2]=0.5;
    
    rdcutsvalmine[0][3]=0.45;
    rdcutsvalmine[1][3]=0.03;
    rdcutsvalmine[2][3]=0.7;
    rdcutsvalmine[3][3]=0.9;
    rdcutsvalmine[4][3]=0.9;
    rdcutsvalmine[5][3]=0.1;
    rdcutsvalmine[6][3]=0.1;
    rdcutsvalmine[7][3]=0.000002;
    rdcutsvalmine[8][3]=0.8;
    rdcutsvalmine[9][3]=200.;
    rdcutsvalmine[10][3]=200.;
    rdcutsvalmine[11][3]=0.021;
    rdcutsvalmine[12][3]=100.;
    rdcutsvalmine[13][3]=0.5;
    
    rdcutsvalmine[0][4]=0.45;
    rdcutsvalmine[1][4]=0.03;
    rdcutsvalmine[2][4]=0.7;
    rdcutsvalmine[3][4]=1.;
    rdcutsvalmine[4][4]=1.;
    rdcutsvalmine[5][4]=0.1;
    rdcutsvalmine[6][4]=0.1;
    rdcutsvalmine[7][4]=0.000002;
    rdcutsvalmine[8][4]=0.8;
    rdcutsvalmine[9][4]=200.;
    rdcutsvalmine[10][4]=200.;
    rdcutsvalmine[11][4]=0.021;
    rdcutsvalmine[12][4]=100.;
    rdcutsvalmine[13][4]=0.5;
    
    rdcutsvalmine[0][5]=0.45;
    rdcutsvalmine[1][5]=0.03;
    rdcutsvalmine[2][5]=0.7;
    rdcutsvalmine[3][5]=1.;
    rdcutsvalmine[4][5]=1.;
    rdcutsvalmine[5][5]=0.1;
    rdcutsvalmine[6][5]=0.1;
    rdcutsvalmine[7][5]=0.000002;
    rdcutsvalmine[8][5]=0.8;
    rdcutsvalmine[9][5]=200.;
    rdcutsvalmine[10][5]=200.;
    rdcutsvalmine[11][5]=0.021;
    rdcutsvalmine[12][5]=100.;
    rdcutsvalmine[13][5]=0.5;
    
    rdcutsvalmine[0][6]=0.45;
    rdcutsvalmine[1][6]=0.03;
    rdcutsvalmine[2][6]=0.7;
    rdcutsvalmine[3][6]=1.;
    rdcutsvalmine[4][6]=1.;
    rdcutsvalmine[5][6]=0.1;
    rdcutsvalmine[6][6]=0.1;
    rdcutsvalmine[7][6]=0.000002;
    rdcutsvalmine[8][6]=0.8;
    rdcutsvalmine[9][6]=200.;
    rdcutsvalmine[10][6]=200.;
    rdcutsvalmine[11][6]=0.021;
    rdcutsvalmine[12][6]=100.;
    rdcutsvalmine[13][6]=0.5;
    
    rdcutsvalmine[0][7]=0.45;
    rdcutsvalmine[1][7]=0.03;
    rdcutsvalmine[2][7]=0.7;
    rdcutsvalmine[3][7]=1.;
    rdcutsvalmine[4][7]=1.;
    rdcutsvalmine[5][7]=0.1;
    rdcutsvalmine[6][7]=0.1;
    rdcutsvalmine[7][7]=0.000002;
    rdcutsvalmine[8][7]=0.8;
    rdcutsvalmine[9][7]=200.;
    rdcutsvalmine[10][7]=200.;
    rdcutsvalmine[11][7]=0.021;
    rdcutsvalmine[12][7]=100.;
    rdcutsvalmine[13][7]=0.5;

  }
  
  RDHFDStartoKpipi->SetCuts(nvars,nptbins,rdcutsvalmine);
  RDHFDStartoKpipi->PrintAll();
  
 
  //CREATE THE TASK
  printf("CREATE TASK\n");
  // create the task
  AliAnalysisTaskSEDStarSpectra *task = new AliAnalysisTaskSEDStarSpectra("AliAnalysisTaskSEDStarSpectra",RDHFDStartoKpipi);
  task->SetAnalysisType(anaType);
  task->SetNSigmasPID(numberOfSigmasPID);
  task->SetMC(theMCon);
  task->SetPID(usePID);
  task->SetDebugLevel(0);

  mgr->AddTask(task);

  // Create and connect containers for input/output
  
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  if (anaType == 0) outputfile += ":PWG3_D2H_DStarSpectraHD";
  if (anaType == 1) outputfile += ":PWG3_D2H_DStarSpectraUU";
  
  // ------ input data ------

  //AliAnalysisDataContainer *cinput0  = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *cinput0  =  mgr->CreateContainer("indstar",TChain::Class(), 
							     AliAnalysisManager::kInputContainer);
  // ----- output data -----
  
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("chist1",TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
  AliAnalysisDataContainer *coutputDStar1 = mgr->CreateContainer("DStarSpectrum",TList::Class(),AliAnalysisManager::kOutputContainer, outputfile.Data());
  AliAnalysisDataContainer *coutputDStar2 = mgr->CreateContainer("DStarAll",TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
  AliAnalysisDataContainer *coutputDStar3 = mgr->CreateContainer("DStarPID3",TList::Class(),AliAnalysisManager::kOutputContainer, outputfile.Data());
  AliAnalysisDataContainer *coutputDStar4 = mgr->CreateContainer("DStarPID2",TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
  AliAnalysisDataContainer *coutputDStar5 = mgr->CreateContainer("DStarPID1",TList::Class(),AliAnalysisManager::kOutputContainer, outputfile.Data());
  AliAnalysisDataContainer *coutputDStar6 = mgr->CreateContainer("cuts",AliRDHFCutsDStartoKpipi::Class(),AliAnalysisManager::kOutputContainer, outputfile.Data()); //cuts

  mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());

  mgr->ConnectOutput(task,1,coutput1);
  mgr->ConnectOutput(task,2,coutputDStar1);
  mgr->ConnectOutput(task,3,coutputDStar2);
  mgr->ConnectOutput(task,4,coutputDStar3);
  mgr->ConnectOutput(task,5,coutputDStar4);
  mgr->ConnectOutput(task,6,coutputDStar5);
  mgr->ConnectOutput(task,7,coutputDStar6);
  
  return task;
}

