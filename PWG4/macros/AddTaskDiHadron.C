//#include "exception.h"
//For running at CERN
AliAnalysisTask *AddTaskDiHadron(){
  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_DiHadron", "No analysis manager found.");
    return 0;
  }

  
  //=========  Set initial parameters=====

  //Track Quality Cuts
  Int_t MinimumClustersTPC=70;
  Float_t MinClusterRatio=0.1;//1/2 would remove split tracks if not for sharing of clusters
  Float_t MaxTPCchi2=4;
  Int_t MinimumClustersITS=0;
  Float_t EtaCut=0.8;//Tracks in +/- Eta are used
  Float_t TriggerEtaCut=0.5;//Smaller trigger for flat acceptance on near-side
  Float_t NearPhiCut=1.5;//Cut used to seperate near and away side for delta eta plots
  Float_t XECut=NearPhiCut;//For XE distribution near and away seperation
  Float_t MaxDCA=3;//Total DCA Cut
  Float_t MaxDCAXY=2.4;
  Float_t MaxDCAZ=3.2;
  Int_t DCAMethod=2;//0 MaxDCA used, 1 MaxDCAXY and MaxDCAZ used 2 pT dependent DCA cut
  Int_t TPCRefit=1;
  Int_t ITSRefit=1;//1 for all particles, 2 for particles above 5 GeV/c
  Int_t SPDCut=1;//check for a point in 1 of first 2 layers of the its
 Float_t MinimumPt=0.25;//Minimum Pt considered by the code
  Float_t MaximumPt=50;
  Float_t ZVertexCut=10;//in cm

  Int_t EfficiencyCorrection=1;//do efficiency corrections in this code
  Int_t MakeMCHistos=1;//if 0 MC histograms are not made (will be empty if 1 and ran on real data)
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
  const Int_t NTriggerPtBins=12;//max=20
  const Int_t NEventsToMix=100;//max=100
  const Int_t NCentralityBins=4;//max=10
  const Int_t NAssociatedPtBins=25;//max=50
  const Int_t N3ParticleAssociatedPtBins=10;//max=50
  const Int_t NZVertexBinsForMixing=7;//max=20
  const Int_t NXEBins=11;//max=20
  const Int_t NumberOfTriggerIDs=1;
  Float_t EffFitPtCut=3;

  TF1 *EfficiencyFitLow=new TF1("EfficiencyFitLow","[0]/[1]*exp(-0.5*pow(x/[1],2))+[2]+[3]*x",MinimumPt,EffFitPtCut);
   TF1 *EfficiencyFitHigh=new TF1("EfficiencyFitHigh","[0]",EffFitPtCut,MaximumPt);
   const Int_t NParamFitLow=4;
   const Int_t NParamFitHigh=1;

   //Not high enough occupancy to worry about the centrality in pp
   //For overlapping centrality bins efficiencies from first bin are used
   //7Pythia_LHC10b5
   Float_t FitLowParam[NCentralityBins*NParamFitLow]={
     -0.030749, 0.254311, 0.858824, -0.0323708,
     -0.0304332, 0.252195, 0.851405, -0.03164,
     -0.0295618, 0.248594, 0.869159, -0.0354148,
     -0.0300529, 0.236693, 0.875875, -0.0370379};
   
   Float_t FitHighParam[NCentralityBins*NParamFitHigh]={
     0.75813,
     0.750521,
     0.751902,
     0.68658};

   /*
   //LHC10c6_900Pythia
   Float_t FitLowParam[NCentralityBins*NParamFitLow]={
     -0.027393, 0.236723, 0.814427, -0.023897,
     -0.0271116, 0.232689, 0.809889, -0.0246341,
     -0.0284658, 0.245163, 0.856043, -0.0318309,
     -0.117114, 0.0355117, 0.828829, -0.0207492}
   
   Float_t FightHighParam[NCentralityBins*NParamFitHigh]={
     0.729888,
     0.719402,
     0.708409,
     0.829268}
   */


  Float_t TriggerPtBins[(NTriggerPtBins+1)]={2,2.5,3,4,6,8,10,15,20,30,40,50};
  Float_t AssociatedPtBins[(NAssociatedPtBins+1)]={0.25,0.5,0.75,1,1.5,2,2.5,3,3.5,4,4.5,5,6,7,8,9,10,12,15,20,25,30,40,50,70,100};
  Float_t AssociatedPtBins31[N3ParticleAssociatedPtBins]={0.5,1.0,1.5,2.0,3,4,1};
  Float_t AssociatedPtBins32[N3ParticleAssociatedPtBins]={1.0,1.5,2.0,3.0,4,5,2};
  Int_t CentralityBins1[NCentralityBins]={0,     0, 20, 40};
  Int_t CentralityBins2[NCentralityBins]={500,20,40,500};
  Float_t XEBins[(NXEBins+1)]={0,0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.6,0.7,0.8,1};
  char *TriggerIDArray[NumberOfTriggerIDs]={"CINT1B"};

  //==================================
  //Add the task
  //===================================
  AliAnalysisTaskDiHadron *task = new AliAnalysisTaskDiHadron("DiHadron");
  task->SetCuts(MinimumClustersTPC,MinClusterRatio,MaxTPCchi2,MinimumClustersITS, EtaCut,TriggerEtaCut,NearPhiCut,XECut,MaxDCA,MaxDCAXY,MaxDCAZ, DCAMethod, TPCRefit,ITSRefit,SPDCut,MinimumPt,MaximumPt,ZVertexCut,NumberOfTriggerIDs,TriggerIDArray);
  task->SetOptions(EfficiencyCorrection,DEBUG,MakeMCHistos);
  task->SetBins(nBinPhi,nBinEta,nBinsPhiEtaPhi,nBinsPhiEtaEta,nBinsPhiPhi,nBinsEtaEta,PhiPlotMin,PhiPlotMax,NTriggerPtBins,NEventsToMix,NCentralityBins,NAssociatedPtBins,N3ParticleAssociatedPtBins,NZVertexBinsForMixing,NXEBins,TriggerPtBins,AssociatedPtBins,AssociatedPtBins31,AssociatedPtBins32,CentralityBins1,CentralityBins2,XEBins);
  task->SetEfficiencies(EffFitPtCut,EfficiencyFitLow,EfficiencyFitHigh,NParamFitLow,NParamFitHigh,FitLowParam,FitHighParam);
 
  mgr->AddTask(task);


  //================================================
  //              data containers
  //================================================
 AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  mgr->ConnectInput(task, 0, cinput); 

  
  //            define output containers, please use 'username'_'somename'
  AliAnalysisDataContainer *coutput1 = 
    mgr->CreateContainer("DiHadron", TList::Class(),
			     AliAnalysisManager::kOutputContainer,"DiHadron.root");

  //           connect containers
  mgr->ConnectInput  (task,  0, cinput );
  mgr->ConnectOutput (task,  0, coutput1);

  return task;
}
