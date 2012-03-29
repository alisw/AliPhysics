//#include "exception.h"
//For running on PbPb data 0-50% most central
AliAnalysisTask *AddTask_DiHadron(Int_t IncludeLowPtBins=0){

  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_sma_PriVtx", "No analysis manager found.");
    return 0;
  }

 
  //========= SetInitial Parameters =====
  //Int_t IncludeLowPtBins=0;//Set to 1 to include low pt triggers

  //Track Quality Cuts
  Int_t MinimumClustersTPC=70;
  Float_t MinClusterRatio=0.51;//Must have at least this ratio not shared
  Float_t MaxTPCchi2=4;
  Int_t MinimumClustersITS=0;
  Float_t EtaCut=0.9;//Tracks in +/- Eta are used
  Float_t TriggerEtaCut=0.5;//Trigger particle restriction so flat area in acceptanc can be created
  Float_t NearPhiCut=1.5;//Cut used to seperate near and away side for delta eta plots
  Float_t XECut=NearPhiCut;//For XE distribution near and away seperation
  Float_t MaxDCA=3;//Total DCA Cut
  Float_t MaxDCAXY=2.4;
  Float_t MaxDCAZ=3.2;
  Int_t DCAMethod=1;//0 MaxDCA used, 1 MaxDCAXY and MaxDCAZ used 2 pT dependent DCA cut
  Int_t TPCRefit=1;
  Int_t ITSRefit=0;//1 for all particles, 2 for particles above 5 GeV/c
  Int_t SPDCut=0;//check for a point in 1 of first 2 layers of the its
 Float_t MinimumPt=0.25;//Minimum Pt considered by the code
  Float_t MaximumPt=50;
  Float_t ZVertexCut=10;//in cm

  //Options
  Int_t RunOnAOD=1;
  Int_t EfficiencyCorrection=1;//do efficiency corrections in this code
  Int_t MakeMCHistos=0;//if 0 MC histograms are not made (will be empty if 1 and ran on real data)
  Int_t DEBUG=0;//for debugging
 
  
 
  //Binning
  Int_t nBinPhi=60;//Number of bins for #Delta#phi histograms
  Int_t nBinEta=54;//Number of bins for #Delta#eta histograms
  Int_t nBinsPhiEtaPhi=20;//Number of bins for #Delta#phi-#Delta#eta in #Delta#phi
  Int_t nBinsPhiEtaEta=18;//Number of bins for #Delta#phi-#Delta#eta in #Delta#phi
  Int_t nBinsPhiPhi=30;//Number of bins for #Delta#phi-#Delta#phi
  Int_t nBinsEtaEta=27;//Number of bins for #Delta#eta-#Delta#eta  
  const Float_t fPi=3.1415926535898;
  Float_t PhiPlotMin=-fPi/3;//Min bin edge in #Delta#phi
  Float_t PhiPlotMax=2*fPi+PhiPlotMin;//Max bin edge
  
  //Size of some arrays change array contents below
  const Int_t NTriggerPtBins2=11;//max20
  const Int_t NTriggerPtBins1=8;//max=20
  Int_t NTriggerPtBins=NTriggerPtBins1;
  if(IncludeLowPtBins)NTriggerPtBins=NTriggerPtBins2;
  const Int_t NEventsToMix=10;//max=100
  const Int_t NCentralityBins=7;//max=10  //6
  const Int_t PercentageCentralityBins=1;//0 or 1
  const Int_t NAssociatedPtBins=25;//max=50
  const Int_t N3ParticleAssociatedPtBins=4;//max=50
  const Int_t NZVertexBinsForMixing=1;//max=20
  const Int_t NXEBins=1;//max=20
  const Int_t NumberOfTriggerIDs=1;
  Float_t EffFitPtCut=3;

  TF1 *EfficiencyFitLow=new TF1("EfficiencyFitLow","[0]/[1]*exp(-0.5*pow(x/[1],2))+[2]+[3]*x+[4]*x**2+[5]*x**3",MinimumPt,EffFitPtCut);
  TF1 *EfficiencyFitHigh=new TF1("EfficiencyFitHigh","[0]+[1]*(x-3)",EffFitPtCut,MaximumPt);
  const Int_t NParamFitLow=6;
  const Int_t NParamFitHigh=2;
  
  //Not high enough occupancy to worry about the centrality in pp
  //For overlapping centrality bins efficiencies from first bin are used
  //For 1% bins from pass2 
Float_t FitLowParam[NCentralityBins*NParamFitLow]={
  -0.0227546, 0.21379, 0.916013, 0.0586031, -0.0429047, 0.0064962,
  -0.0224257, 0.215105, 0.916696, 0.0593096, -0.0429726, 0.00650312,
  -0.0220935, 0.208194, 0.915207, 0.0613371, -0.0416018, 0.00576639,
  -0.0229694, 0.208645, 0.919995, 0.0527738, -0.0344433, 0.00430404,
  -0.0235676, 0.206099, 0.922875, 0.0476261, -0.0298983, 0.0032692,
  -0.0238139, 0.206366, 0.925183, 0.0430703, -0.0264516, 0.00250388,
  -0.024815, 0.204409, 0.928448, 0.03423, -0.019894, 0.00121746
};
 Float_t FitHighParam[NCentralityBins*NParamFitHigh]={
   0.881078, 0.00158217,
   0.883456, 0.00144851,
   0.880495, 0.00369636,
   0.884535, 0.00297087,
   0.884937, 0.00166201,
   0.883934, 0.00290625,
   0.884964, -0.000784607
 };
 
  Float_t V2FitPtCut=2.75;
  Float_t V3FitPtCut=2.75;
  Float_t V4FitPtCut=2.75;
  TF1 *V2FitLow=new TF1("V2FitLow","[0]*(x+[1]*x**2+[2]*x**3+[3]*x**4)",MinimumPt,V2FitPtCut);
  TF1 *V2FitHigh=new TF1("V2FitHigh","[0]*exp(-0.5*((x-[1])/[2])**2)+[3]",V2FitPtCut,MaximumPt);
  TF1 *V3FitLow=new TF1("V3FitLow","[0]*(x+[1]*x**2+[2]*x**3+[3]*x**4)",MinimumPt,V2FitPtCut);
  TF1 *V3FitHigh=new TF1("V3FitHigh","[0]*exp(-0.5*((x-[1])/[2])**2)+[3]",V2FitPtCut,MaximumPt);
 TF1 *V4FitLow=new TF1("V4FitLow","[0]*(x+[1]*x**2+[2]*x**3+[3]*x**4)",MinimumPt,V2FitPtCut);
  TF1 *V4FitHigh=new TF1("V4FitHigh","[0]*exp(-0.5*((x-[1])/[2])**2)+[3]",V2FitPtCut,MaximumPt);

  const Int_t NParamV2FitLow=4;
  const Int_t NParamV2FitHigh=4;
  const Int_t NParamV3FitLow=4;
  const Int_t NParamV3FitHigh=4;
  const Int_t NParamV4FitLow=4;
  const Int_t NParamV4FitHigh=4;

  //From FlowScale2 macro
Float_t FitLowParamV2[NCentralityBins*NParamV2FitLow]={
0.046931, -0.504895, 0.186345, -0.030504,
0.060953, -0.503882, 0.197361, -0.030849,
0.088939, -0.396433, 0.127467, -0.016969,
0.122602, -0.347615, 0.098626, -0.012576,
0.155403, -0.296297, 0.065681, -0.007492,
0.177546, -0.271958, 0.049388, -0.005644,
0.187634, -0.223967, 0.011281, 0.001462};
Float_t FitHighParamV2[NCentralityBins*NParamV2FitHigh]={
0.000000, 4.478095, 0.500000, 0.049007,
0.012659, 3.603740, 0.656049, 0.071827,
0.043836, 3.623140, 1.414293, 0.089569,
0.071277, 3.490653, 1.454361, 0.112755,
0.082163, 3.324474, 1.429585, 0.145809,
0.085348, 3.194132, 1.562318, 0.163847,
0.061056, 2.994950, 1.442946, 0.194451};
Float_t FitLowParamV3[NCentralityBins*NParamV3FitLow]={
0.019078, 0.793599, -0.436079, 0.075946,
0.025818, 0.088555, 0.070194, -0.026702,
0.028562, 0.297268, -0.112696, 0.013257,
0.034723, 0.137282, -0.008578, -0.009524,
0.041567, -0.002427, 0.093642, -0.032998,
0.045574, 0.059023, 0.029342, -0.019022,
0.052432, -0.067332, 0.087218, -0.028239};
Float_t FitHighParamV3[NCentralityBins*NParamV3FitHigh]={
0.079111, 3.802577, 3.182298, 0.000000,
0.105446, 4.216029, 2.377002, 0.000000,
0.106553, 4.500000, 2.825647, 0.009652,
0.099684, 3.875928, 1.968600, 0.022106,
0.103187, 3.597519, 1.980556, 0.022403,
0.123473, 3.773586, 2.099909, 0.011130,
0.140945, 3.618698, 1.852780, 0.000000};
Float_t FitLowParamV4[NCentralityBins*NParamV4FitLow]={
0.005289, 1.000000, 0.398581, -0.095251,
0.018611, -1.000000, 0.837275, -0.163058,
0.022107, -1.000000, 0.847701, -0.168469,
0.025434, -1.000000, 0.854238, -0.176690,
0.030768, -1.000000, 0.811683, -0.164335,
0.037000, -1.000000, 0.781250, -0.158122,
0.043484, -1.000000, 0.752546, -0.153147};
  Float_t FitHighParamV4[NCentralityBins*NParamV4FitHigh]={
0.061365, 3.669890, 1.724126, 0.014456,
0.086516, 4.194769, 1.723464, 0.000000,
0.030263, 4.008413, 0.563914, 0.067831,
0.073456, 4.255665, 1.551659, 0.025494,
0.058771, 4.102898, 1.191926, 0.048796,
0.054682, 4.221963, 1.119594, 0.064515,
0.094564, 4.500000, 1.568750, 0.038027};



  if(IncludeLowPtBins){ Float_t TriggerPtBins[(NTriggerPtBins2+1)]={0.75,1,2,2.5,3,4,6,8,10,15,20,25};}
  else{  Float_t TriggerPtBins[(NTriggerPtBins1+1)]={2.5,3,4,6,8,10,15,20,25};}

  Float_t AssociatedPtBins[(NAssociatedPtBins+1)]={0.25,0.5,0.75,1,1.5,2,2.5,3,3.5,4,4.5,5,6,7,8,9,10,12,15,20,25,30,40,50,70,100};
 
  Float_t AssociatedPtBins31[N3ParticleAssociatedPtBins]={0.5,0.75,1,2};
  Float_t AssociatedPtBins32[N3ParticleAssociatedPtBins]={0.75,1,2,3};
 
  Int_t CentralityBins1[NCentralityBins]={0,0,5,10,20,30,40};
  Int_t CentralityBins2[NCentralityBins]={2,5,10,20,30,40,50};

  Float_t XEBins[(NXEBins+1)]={0,0.01};
  //char *TriggerIDArray="CINT1B";//seperate multiple with ,
  //char *TriggerIDArray="CSH1-B";//seperate multiple with ,
  char *TriggerIDArray="C";//PbPb test, using SelectCollisionsCanidates instead

  ////////////////////////
  //Add the task
  ////////////////////////
  AliAnalysisTaskDiHadron *task = new AliAnalysisTaskDiHadron("julery_DiHadron");
  task->SetCuts(MinimumClustersTPC,MinClusterRatio,MaxTPCchi2,MinimumClustersITS, EtaCut,TriggerEtaCut,NearPhiCut,XECut,MaxDCA,MaxDCAXY,MaxDCAZ, DCAMethod, TPCRefit,ITSRefit,SPDCut,MinimumPt,MaximumPt,ZVertexCut,NumberOfTriggerIDs,TriggerIDArray);
  task->SetOptions(RunOnAOD,EfficiencyCorrection,DEBUG,MakeMCHistos);
  task->SetBins(nBinPhi,nBinEta,nBinsPhiEtaPhi,nBinsPhiEtaEta,nBinsPhiPhi,nBinsEtaEta,PhiPlotMin,PhiPlotMax,NTriggerPtBins,NEventsToMix,NCentralityBins,PercentageCentralityBins,NAssociatedPtBins,N3ParticleAssociatedPtBins,NZVertexBinsForMixing,NXEBins,TriggerPtBins,AssociatedPtBins,AssociatedPtBins31,AssociatedPtBins32,CentralityBins1,CentralityBins2,XEBins);
  task->SetEfficiencies(EffFitPtCut,EfficiencyFitLow,EfficiencyFitHigh,NParamFitLow,NParamFitHigh,FitLowParam,FitHighParam);
  task->SetFlow(V2FitPtCut,V3FitPtCut,V4FitPtCut,V2FitLow,V2FitHigh,V3FitLow,V3FitHigh,V4FitLow,V4FitHigh,NParamV2FitLow,NParamV2FitHigh,NParamV3FitLow,NParamV3FitHigh,NParamV4FitLow,NParamV4FitHigh,FitLowParamV2,FitHighParamV2,FitLowParamV3,FitHighParamV3,FitLowParamV4,FitHighParamV4);

// physics selection
Int_t isMC=0;//1 for MC 0 for DATA
//gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
//AliPhysicsSelectionTask *PhysicsTask=AddTaskPhysicsSelection(isMC, 0); //isMC is true when processing monte carlo, the second 0 disables the cluster vs tracklets
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
  //AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection(isMC,0);
  task->SelectCollisionCandidates(AliVEvent::kMB);
  //AliCentralitySelectionTask *centSelTask = AliCentralitySelectionTask("CentralitySelection");
 
  mgr->AddTask(task);
 

  //================================================
  //              data containers
  //================================================
  //            find input container
  //below the trunk version
  //   AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer();
  //this is the old way!!!
  AliAnalysisDataContainer *cinput  = (AliAnalysisDataContainer*)mgr->GetContainers()->FindObject("cAUTO_INPUT");
  
  //            define output containers, please use 'username'_'somename'
  AliAnalysisDataContainer *coutput1 = 
    mgr->CreateContainer("julery_DiHadron", TList::Class(),
			     AliAnalysisManager::kOutputContainer,"julery_DiHadron.root");

  //           connect containers
  mgr->ConnectInput  (task,  0, cinput );
  mgr->ConnectOutput (task,  0, coutput1);

  return task;
  
}
