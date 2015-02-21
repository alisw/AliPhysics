void InitHistograms(AliDielectron *die, Int_t cutDefinition);
void InitCF(AliDielectron* die, Int_t cutDefinition);

void SetupEventCuts(AliDielectron *die, Int_t cutDefinition);
void SetupTrackCuts(AliDielectron *die, Int_t cutDefinition);
void SetupPairCuts(AliDielectron *die,  Int_t cutDefinition);
//QAtask

void SetupV0cuts(AliDielectron *die, Int_t cutDefinition);
void SetupV0Add(AliDielectron *die, Int_t cutDefinition);
TVectorD *GetRunNumbers();
/*
namespace ConfDef {
	
	void InitHistograms(AliDielectron *die, Int_t cutDefinition);
	void InitCF(AliDielectron* die, Int_t cutDefinition);
	void SetupV0Cuts(AliDielectron *die, Int_t cutDefinition);
	void SetupV0Add(AliDielectron *die, Int_t cutDefinition);
	void SetupTrackCuts(AliDielectron *die, Int_t cutDefinition);
	void SetupPairCuts(AliDielectron *die, Int_t cutDefinition);
	void SetEtaCorrection();	
	enum {kDefault,kHF,kLmee,kDefault_activevolume,kDefault_conversions,kDefault_conversions_wPID};
}
*/

enum ConfDef {kDefault,kHF,kLmee,kDefault_activevolume,kDefault_conversions,kDefault_conversions_wPID};
void SetupMCsignals(AliDielectron *die);




TVectorD *GetRunNumbers() {
// returns a vector with the runnumber used in the period                                                                                                                                                                                                                                                            
  Double_t first=0;
  Double_t last =1;

 
  switch(iPeriod) {
  case k10b: first=114737; last=117223; break;
  case k10c: first=117777; last=121417; break;
  case k10d: first=121692; last=126437; break;
  case k10e: first=127102; last=130850; break;
  case k10f: first=130931; last=135031; break;
  case k10h: first=136831; last=139517; break;
  case k11a: first=141052; last=146974; break;
  case k11d: first=155838; last=159649; break;
  case k11h: first=165772; last=170718; break;
  case k12h: first=188720; last=192738; break;
  }
  // printf("iPeriod: %d \t %.0f-%.0f \n",iPeriod,first,last);                                                                                                                                                                                                                                                        
  return (AliDielectronHelper::MakeLinBinning(last-first, first, last));
}




//TString run_numbers;
//if( list.Contains("LHC10b")){

TString run_numbers = "114786,114798,114918,114920,114924,114930,114931,115186,115193,115310,115318,115322,115328,115335,115345,115393,115399,115401,115414,115521,116079,116081,116102,116288,116402,116403,116562,116571,116574,116643,116645,117048,117050,117052,117053,117059,117060,117063,117092,117099,117109,117112,117116,117220,117222,118506,118507,118512,118518,118556,118558,118560,118561,119159,119161,119163,119841,119842,119844,119845,119846,119849,119853,119856,119859,119863,120067,120069,120072,120073,120076,120079,120244,120503,120504,120505,120616,120617,120671,120741,120750,120758,120820,120821,120822,120823,120824,120825,120829,121039,121040,122370,122374,122375,124751,125023,125085,125097,125100,125101,125134,125296,125628,125630,125632,125633,125842,125843,125844,125847,125848,125849,125850,125851,125855,125941,125942,125943,125944,125945,125946,125947,125948,125949,125950,125951,125952,125953,125954,125955,125956,125957,125958,125959,125960,125961,125962,125963,125964,125965,125966,125969,125970,125976,125978,125981,125982,125983,125984,125985,125986,125997,125998,125999,126004,126007,126008,126073,126078,126081,126082,126088,126090,126097,126158,126160,126168,126283,126284,126285,126351,126352,126359,126403,126404,126405,126406,126407,126408,126409,126422,126424,126425,126432";

/*
TVectorD *GetRunNumbers() {
  Double_t runLHC10b[]={114786,114920,126422,126424,126425,126432};
  }
  
  Double_t runLHC10h[] = {139510, 139507, 139505, 139503, 139465, 139438, 139437, 139360, 139329, 139328, 139314, 139310, 139309, 139173, 139107, 139105, 139038, 139037, 139036, 139029, 139028, 138872, 138871, 138870, 138837, 138732, 138730, 138666, 138662, 138653, 138652, 138638, 138624, 138621, 138583, 138582, 138579, 138578, 138534, 138469, 138442, 138439, 138438, 138396, 138364, 138275, 138225, 138201, 138197, 138192, 138190, 137848, 137844, 137752, 137751, 137724, 137722, 137718, 137704, 137693, 137692, 137691, 137686, 137685, 137639, 137638, 137608, 137595, 137549, 137546, 137544, 137541, 137539, 137531, 137530, 137443, 137441,137440, 137439, 137434, 137432, 137431, 137430, 137366, 137243, 137236, 137235, 137232, 137231, 137230, 137162, 137161, 137135};
  Double_t runLHC11h[] = {167915, 167920, 167985, 167987, 167988, 168069, 168076, 168105, 168107, 168108, 168115, 168310, 168311, 168322, 168325, 168341, 168342, 168361, 168362, 168458, 168460, 168464, 168467, 168511, 168512, 168514, 168777, 168826, 168988, 168992, 169035, 169040, 169044, 169045, 169091, 169094, 169099, 169138, 169144, 169145, 169148, 169156, 169160, 169167, 169238, 169411, 169415, 169417, 169418, 169419, 169420, 169475, 169498, 169504, 169506, 169512, 169515, 169550, 169553, 169554, 169555, 169557, 169586, 169587, 169588, 169590, 169591, 169835, 169837, 169838, 169846, 169855, 169858, 169859, 169923, 169965, 170027,170040, 170081, 170083, 170084, 170085, 170088, 170089, 170091, 170155, 170159, 170163, 170193, 170203, 170204, 170207, 170228, 170230, 170268, 170269, 170270, 170306, 170308, 170309, 170311, 170312, 170313, 170315, 170387, 170388, 170572, 170593,170593+1};

  Int_t size = (int) (sizeof(runLHC10b)/sizeof(Double_t));
  TVectorD *vec = new TVectorD(size,runLHC10b);
  //vec->Print("");                                                                                                                                                                                                                                                                                                   
  return vec;
}
*/


  //TString run_numbers = "114786,114798,114918,114920,114924,114930,114931,115186,115193,115310,115318,115322,115328,115335,115345,115393,115399,115401,115414,115521,116079,116081,116102,116288,116402,116403,116562,116571,116574,116643,116645,117048,117050,117052,117053,117059,117060,117063,117092,117099,117109,117112,117116,117220,117222,118506,118507,118512,118518,118556,118558,118560,118561,119159,119161,119163,119841,119842,119844,119845,119846,119849,119853,119856,119859,119863,120067,120069,120072,120073,120076,120079,120244,120503,120504,120505,120616,120617,120671,120741,120750,120758,120820,120821,120822,120823,120824,120825,120829,121039,121040,122370,122374,122375,124751,125023,125085,125097,125100,125101,125134,125296,125628,125630,125632,125633,125842,125843,125844,125847,125848,125849,125850,125851,125855,125941,125942,125943,125944,125945,125946,125947,125948,125949,125950,125951,125952,125953,125954,125955,125956,125957,125958,125959,125960,125961,125962,125963,125964,125965,125966,125969,125970,125976,125978,125981,125982,125983,125984,125985,125986,125997,125998,125999,126004,126007,126008,126073,126078,126081,126082,126088,126090,126097,126158,126160,126168,126283,126284,126285,126351,126352,126359,126403,126404,126405,126406,126407,126408,126409,126422,126424,126425,126432";
//}


  //TString run_numbers = "114785,114786,114798,114918,114920,114924,114930,114931,115056,115186,115193,115310,115312,115318,115322,115327,115328,115335,115345,115369,115393,115399,115401,115406,115414,115521,116079,116081,116102,116288,116358,116360,116372,116402,116403,116429,116431,116432,116562,116571,116574,116643,116645,116787,117048,117050,117052,117053,117059,117060,117063,117092,117099,117109,117112,117116,117220,117222,118506,118507,118512,118518,118556,118558,118560,118561,119159,119161,119163,119841,119842,119844,119845,119846,119849,119853,119856,119859,119862,120067,120069,120072,120073,120076,120079,120244,120503,120504,120505,120616,120617,120671,120741,120750,120758,120820,120821,120822,120823,120824,120825,120829,121039,121040,122374,122375,124751,125023,125083,125085,125097,125100,125101,125133,125134,125139,125140,125156,125186,125295,125296,125628,125630,125632,125633,125842,125843,125844,125847,125848,125849,125850,125851,125855,126004,126007,126008,126073,126078,126081,126082,126088,126090,126097,126158,126160,126167,126168,126283,126284,126285,126350126351,126352,126359,126403,126404,126405,126406,126407,126408,126409,126422,126424,126425,




//TVectorD *GetRunNumbers();
//*TVectorD *GetRunNumbers() {
// if (period="LHC10d"){
      // return AliDielectronHelper::MakeArbitraryBinning ("122374,122375,124751,125023,125085,125097,125100,125101,125134,125296,125628,125630,125632,125633,125842,125843,125844,125847,125848,125849,125850,125851,125855,125941,125942,125943,125944,125945,125946,125947,125948,125949,125950,125951,125952,125953,125954,125955,125956,125957,125958,125959,125960,125961,125962,125963,125964,125965,125966,125969,125970,125976,125978,125981,125982,125983,125984,125985,125986,125997,125998,125999,126004,126007,126008,126073,126078,126081,126082,126088,126090,126097,126158,126160,126168,126283,126284,126285,126351,126352,126359,126403,126404,126405,126406,126407,126408,126409,126422,126424,126425,126432");
    //   }
//}
  // increasing order of run array is important

//TString names=("default");
TString names=("JPsi;kHFe;Lmee;Jpsi_activevolume;Jpsi_conversions;Jpsi_conversions_wPID");
TObjArray *arrNames=names.Tokenize(";");

const Int_t nDie=arrNames->GetEntries();
Bool_t isAOD=kFALSE;
Bool_t hasMC=kFALSE;
Int_t iPeriod=-1;
enum { k10b=0, k10c, k10d, k10e, k10f, k10h, k11a, k11d, k11h, k12h, k13b, k13c, k13d, k13e, k13f };

//______________________________________________________________________________________
AliAnalysisTask* AddTask_jpsi_Default(TString prod="", Bool_t isMC=kFALSE)
{
  //get the current analysis manager

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_jpsi_JPsi", "No analysis manager found.");
    return 0;
  }

  //Do we have an MC handler?
  hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);

  //AOD input?
  isAOD=mgr->GetInputEventHandler()->IsA()==AliAODInputHandler::Class();
  if(isAOD) hasMC=isMC;

  //Get the current train configuration
  //  TString trainConfig=gSystem->Getenv("CONFIG_FILE");
  //TString list=gSystem->Getenv("LIST");
  //if( list.IsNull()) list=prod;

  //printf("list************* %s",list.Data());
  // selected period
  if(      !prod.CompareTo("LHC10b") ) iPeriod = k10b;
  else if( !prod.CompareTo("LHC10c") ) iPeriod = k10c;
  else if( !prod.CompareTo("LHC10d") ) iPeriod = k10d;
  else if( !prod.CompareTo("LHC10e") ) iPeriod = k10e;
  else if( !prod.CompareTo("LHC10f") ) iPeriod = k10f;
  else if( !prod.CompareTo("LHC10h") ) iPeriod = k10h;
  else if( !prod.CompareTo("LHC11a") ) iPeriod = k11a;
  else if( !prod.CompareTo("LHC11d") ) iPeriod = k11d;
  else if( !prod.CompareTo("LHC11h") ) iPeriod = k11h;
  else if( !prod.CompareTo("LHC12h") ) iPeriod = k12h;
  else if( !prod.CompareTo("LHC13b") ) iPeriod = k13b;
  else if( !prod.CompareTo("LHC13c") ) iPeriod = k13c;
  else if( !prod.CompareTo("LHC12d") ) iPeriod = k13d;
  else if( !prod.CompareTo("LHC12e") ) iPeriod = k13e;
  else if( !prod.CompareTo("LHC12f") ) iPeriod = k13f;

  //  printf("prod***********%s",prod.Data());
  ///
  /*  if(!prod.CompareTo("LHC10b")){
    //TString run_numbers= "117220"                                                                                                                                       
    run_numbers="114786,114798,114918,114824,114930,114931,115186,115193,115310,115318,115322,115328,115335,115345,115393,115399,115401,115414,115521,116079,116081,116102,16288,116402,116403,116562,116571,116574,116643,116645,117048,117050,117052,117053,117059,117060,117063,117092,117099,117109,117112,117116,117220,117222";
  }
  else if (!prod.CompareTo("LHC10d")){
    run_numbers = "122370,122374,122375,124751,125023,125085,125097,125100,125101,125134,125296,125628,125630,125632,125633,125842,125843,125844,125847,125848,125849,125850,125851,125855,125941,125942,125943,125944,125945,125946,125947,125948,125949,125950,125951,125952,125953,125954,125955,125956,125957,125958,125959,125960,125961,125962,125963,125964,125965,125966,125969,125970,125976,125978,125981,125982,125983,125984,125985,125986,125997,125998,125999,126004,126007,126008,126073,126078,126081,126082,126088,126090,126097,126158,126160,126168,126283,126284,126285,126351,126352,126359,126403,126404,126405,126406,126407,126408,126409,126422,126424,126425,126432";
  }
  */
  //
  //////                                                                                                                                                                                                                                                                

  //Do 
 

  // // aod monte carlo
  // if( list.Contains("LHC11a10") ||
  //     list.Contains("LHC11b10") ||
  //     list.Contains("LHC12a17") ||
  //     list.Contains("fix")
  //     ) hasMC=kTRUE;

  //create task and add it to the manager
  AliAnalysisTaskMultiDielectron *task=new AliAnalysisTaskMultiDielectron("JpsiDefault");
  //  task->SetBeamEnergy(1380.); // not neeeded since we are not looking at helicity and Collins-Soper coordinates
  if (!hasMC) task->UsePhysicsSelection();

  // add special triggers
  switch(iPeriod) {
  case k11d: task->SetTriggerMask(AliVEvent::kEMCEJE+AliVEvent::kEMC7+AliVEvent::kEMCEGA);     break;
  case k11h: task->SetTriggerMask(AliVEvent::kMB+AliVEvent::kCentral+AliVEvent::kSemiCentral); break;
  case k12h: task->SetTriggerMask(AliVEvent::kAnyINT); break;                                      
  case k13b: task->SetTriggerMask(AliVEvent::kINT7); break;
  case k13c: task->SetTriggerMask(AliVEvent::kINT7); break;
  case k13d: task->SetTriggerMask(AliVEvent::kAnyINT); break;
  case k13e: task->SetTriggerMask(AliVEvent::kAnyINT); break;
  case k13f: task->SetTriggerMask(AliVEvent::kAnyINT); break;
  }
  mgr->AddTask(task);

  //add dielectron analysis with different cuts to the task
  for (Int_t i=0; i<nDie; ++i){ //nDie defined in config file
    AliDielectron *jpsi=ConfigDefault(i);
    if (!jpsi) continue;
    jpsi->SetHasMC(hasMC);
    task->AddDielectron(jpsi);
  }

  //   task->SetTriggerOnV0AND();
  //   if ( trainConfig=="pp" ) task->SetRejectPileup();
  
  //create output container
  AliAnalysisDataContainer *coutput1 =
    mgr->CreateContainer("jpsi_Default_tree",
                         TTree::Class(),
                         AliAnalysisManager::kExchangeContainer,
                         "jpsi_Default_default");
  
  AliAnalysisDataContainer *cOutputHist1 =
    mgr->CreateContainer("jpsi_Default_QA",
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         "jpsi_Default.root");

  AliAnalysisDataContainer *cOutputHist2 =
    mgr->CreateContainer("jpsi_Default_CF",
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         "jpsi_Default.root");
  
  AliAnalysisDataContainer *cOutputHist3 =
    mgr->CreateContainer("jpsi_Default_EventStat",
                         TH1D::Class(),
                         AliAnalysisManager::kOutputContainer,
                         "jpsi_Default.root");
  
  mgr->ConnectInput(task,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 0, coutput1 );
  mgr->ConnectOutput(task, 1, cOutputHist1);
  mgr->ConnectOutput(task, 2, cOutputHist2);
  mgr->ConnectOutput(task, 3, cOutputHist3);
  
  return task;
}


//______________________________________________________________________________________
//______________________________________________________________________________________
//______________________________________________________________________________________
//
// Here the configuration part starts
//
AliDielectron* ConfigDefault(Int_t cutDefinition)
{
  //
  // Setup the instance of AliDielectron
  //
  
  // create the actual framework object
  TString name=Form("%02d",cutDefinition);
  if (cutDefinition<arrNames->GetEntriesFast()){
    name=arrNames->At(cutDefinition)->GetName();
  }
  AliDielectron *die =
  new AliDielectron(Form("%s",name.Data()),
                    Form("Track cuts: %s",name.Data()));
  
  /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv CUTS vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
  SetupEventCuts(die,cutDefinition);
  SetupTrackCuts(die,cutDefinition);
  if (cutDefinition ==ConfDef::kDefault){
    SetupPairCuts(die,cutDefinition);
  }

  if (cutDefinition ==ConfDef::kDefault || cutDefinition ==ConfDef::kLmee || cutDefinition ==ConfDef::kDefault_activevolume){
		SetupV0Cuts(die,cutDefinition);
	}
		
//V0s to have a pure e+e- sample to check the TPC nsigma	
  if (cutDefinition ==ConfDef::kDefault_conversions || cutDefinition ==ConfDef::kDefault_conversions_wPID){
    SetupV0add(die,cutDefinition);
  }

  
  if (cutDefinition !=ConfDef::kDefault){
  die->SetNoPairing();
  }	
  
  //SetupV0Cuts(die,cutDefinition); 


  /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv MISC vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
  // Monte Carlo Signals
  if (hasMC) SetupMCsignals(die);
  // prefilter settings
  // die->SetPreFilterUnlikeOnly();//  die->SetNoPairing();//  die->SetPreFilterAllSigns();
  // cut QA
  // die->SetCutQA();

  /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv OUTPUT vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
  InitHistograms(die,cutDefinition);
  //No CF container for now
   if (cutDefinition ==0){
   InitCF(die,cutDefinition);
   }
  //   AliDielectronMixingHandler *mix=new AliDielectronMixingHandler;
  //   mix->AddVariable(AliDielectronVarManager::kZvPrim,"-10,-8,-5,0,5,8,10");
  //   mix->AddVariable(AliDielectronVarManager::kPhi ,"0,3.14,6.3");
  //   mix->SetDepth(10);
  //  die->SetMixingHandler(mix);
  //
  
  return die;
}

//______________________________________________________________________________________
void SetupEventCuts(AliDielectron *die, Int_t cutDefinition)
{
  //
  // Setup the event cuts
  //

  AliDielectronEventCuts *eventCuts=new AliDielectronEventCuts("eventCuts","eventCuts");
  if(isAOD) eventCuts->SetVertexType(AliDielectronEventCuts::kVtxAny);
  eventCuts->SetRequireVertex();
  eventCuts->SetMinVtxContributors(1);
  eventCuts->SetVertexZ(-10.,+10.);
  eventCuts->Print();
  die->GetEventFilter().AddCuts(eventCuts);

}

//______________________________________________________________________________________
void SetupTrackCuts(AliDielectron *die, Int_t cutDefinition)
{
  //
  // Setup the track cuts
  //

  AliDielectronCutGroup* cuts = new AliDielectronCutGroup("cuts","cuts",AliDielectronCutGroup::kCompAND);
  die->GetTrackFilter().AddCuts(cuts);
  
  //default quality cuts
  AliDielectronTrackCuts *refit=new AliDielectronTrackCuts("refit","refit");
  
  if (cutDefinition==ConfDef::kDefault){
    refit->SetRequireITSRefit(kTRUE);
    refit->SetRequireTPCRefit(kTRUE);
    refit->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
  }else if(cutDefinition==ConfDef::kHF){
    refit->SetRequireITSRefit(kTRUE);
    refit->SetRequireTPCRefit(kTRUE);
    refit->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
  }else if(cutDefinition==ConfDef::kLmee){
    refit->SetRequireITSRefit(kTRUE);
    refit->SetRequireTPCRefit(kTRUE);
    refit->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
  }else if (cutDefinition ==ConfDef::kDefault_activevolume){
    refit->SetRequireITSRefit(kTRUE);
    refit->SetRequireTPCRefit(kTRUE);
    refit->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
  }else if (cutDefinition ==ConfDef::kDefault_conversions){
    refit->SetRequireITSRefit(kTRUE);
    refit->SetRequireTPCRefit(kTRUE);
    //    refit->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
  }else if(cutDefinition ==ConfDef::kDefault_conversions_wPID){
    refit->SetRequireITSRefit(kTRUE);
    refit->SetRequireTPCRefit(kTRUE);
    //  refit->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
    
  }
  cuts->AddCut(refit);
  
	
	//pt and kink mother
  AliDielectronVarCuts *pt = new AliDielectronVarCuts("ptCut","pt cut");
  
  if (cutDefinition==kDefault){
    pt->AddCut(AliDielectronVarManager::kPt,1.,1.e30);
    pt->AddCut(AliDielectronVarManager::kEta,-0.9,0.9);
    pt->AddCut(AliDielectronVarManager::kKinkIndex0,0.);	
    pt->AddCut(AliDielectronVarManager::kTPCnSigmaEle,-3.,3.);
    pt->AddCut(AliDielectronVarManager::kTPCnSigmaPio,3.5,1000.);
    pt->AddCut(AliDielectronVarManager::kTPCnSigmaPro,3.,1000.);
    pt->AddCut(AliDielectronVarManager::kNclsTPC,70.,160.);  
    pt->AddCut(AliDielectronVarManager::kTPCchi2Cl,0.,4.);
    //impact parameter
    pt->AddCut(AliDielectronVarManager::kImpactParXY,-1.,1.);
    pt->AddCut(AliDielectronVarManager::kImpactParZ,-3.,3.);
  }
  else if (cutDefinition==kHF) {
    pt->AddCut(AliDielectronVarManager::kPt,1.,1.e30);
    pt->AddCut(AliDielectronVarManager::kEta,-0.8,0.8);
    //to be checked
    pt->AddCut(AliDielectronVarManager::kKinkIndex0,0.);
    pt->AddCut(AliDielectronVarManager::kTPCnSigmaEle,-0.5.,3.);
    if (AliDielectronVarManager::kTOFPIDBit >0.8){
      pt->AddCut(AliDielectronVarManager::kTOFnSigmaEle,-3.,3.);
    }
    //TPC PID ClusteR
    pt->AddCut(AliDielectronVarManager::kTPCsignalN,80.,160.);
    //rjim NTPCcluster cut
    pt->AddCut(AliDielectronVarManager::kNclsTPC,120.,160.);  
    pt->AddCut(AliDielectronVarManager::kNFclsTPCfCross,0.6,1.1);    
    pt->AddCut(AliDielectronVarManager::kTPCchi2Cl,0.,4.);
    //    pt->AddCut(AliDielectronVarManager::kNClsITS,4.,200); 
    pt->AddCut(AliDielectronVarManager::kImpactParXY,-1.,1.);
    pt->AddCut(AliDielectronVarManager::kImpactParZ,-2.,2.);  
  }
  else if (cutDefinition==kLmee) {
    pt->AddCut(AliDielectronVarManager::kPt,0.2,1.e30);
    pt->AddCut(AliDielectronVarManager::kEta,-0.8,0.8);
    pt->AddCut(AliDielectronVarManager::kKinkIndex0,0.);
    pt->AddCut(AliDielectronVarManager::kTPCnSigmaEle,-3.,3.);
    pt->AddCut(AliDielectronVarManager::kTPCnSigmaPio,3.,1000.);
    pt->AddCut(AliDielectronVarManager::kTOFnSigmaEle,-3.,3.); 
    //rjim Nitsclusters cut
    //  pt->AddCut(AliDielectronVarManager::kNclsITS,4.,160.);
    //add ncrossed rows tpc instead cluster instead Ncluster
    pt->AddCut(AliDielectronVarManager::kNclsTPC,70.,160.);
    pt->AddCut(AliDielectronVarManager::kNFclsTPCfCross,0.8,1.0);    
    pt->AddCut(AliDielectronVarManager::kTPCchi2Cl,0.,4.); 
    pt->AddCut(AliDielectronVarManager::kImpactParXY,-1.,1.);
    pt->AddCut(AliDielectronVarManager::kImpactParZ,-3.,3.);
  }else if(cutDefinition==kDefault_activevolume){
    pt->AddCut(AliDielectronVarManager::kPt,1.,1.e30);
    pt->AddCut(AliDielectronVarManager::kEta,-0.9,0.9);
    pt->AddCut(AliDielectronVarManager::kKinkIndex0,0.);
    pt->AddCut(AliDielectronVarManager::kTPCnSigmaEle,-3.,3.);
    pt->AddCut(AliDielectronVarManager::kTPCnSigmaPio,3.5,1000.);
    pt->AddCut(AliDielectronVarManager::kTPCnSigmaPro,3.,1000.);
    //    pt->AddCut(AliDielectronVarManager::kTPCactvol,120.,200.);
    // NTPCclusters
    pt->AddCut(AliDielectronVarManager::kNclsTPC,70.,160.);  
    pt->AddCut(AliDielectronVarManager::kTPCchi2Cl,0.,4.);
    //impact parameter
    pt->AddCut(AliDielectronVarManager::kImpactParXY,-1.,1.);
    pt->AddCut(AliDielectronVarManager::kImpactParZ,-3.,3.);
    //  pt->AddCut(AliDielectronVarManager::kTOFPIDBit,0.8,2.0);
  }else if(cutDefinition==kDefault_conversions){
    pt->AddCut(AliDielectronVarManager::kPt,1.,1.e30);
    pt->AddCut(AliDielectronVarManager::kEta,-0.9,0.9);
    pt->AddCut(AliDielectronVarManager::kKinkIndex0,0.);
    //    pt->AddCut(AliDielectronVarManager::kTPCnSigmaEle,-3.,3.);
    pt->AddCut(AliDielectronVarManager::kTPCnSigmaPio,3.5,1000.);
    pt->AddCut(AliDielectronVarManager::kTPCnSigmaPro,3.,1000.);
    pt->AddCut(AliDielectronVarManager::kNclsTPC,70.,160.);  
    pt->AddCut(AliDielectronVarManager::kTPCchi2Cl,0.,4.);
    pt->AddCut(AliDielectronVarManager::kImpactParXY,-1.,1.);
    pt->AddCut(AliDielectronVarManager::kImpactParZ,-3.,3.);
  }else if(cutDefinition==kDefault_conversions_wPID){
    pt->AddCut(AliDielectronVarManager::kPt,1.,1.e30);
    pt->AddCut(AliDielectronVarManager::kEta,-0.9,0.9);
    pt->AddCut(AliDielectronVarManager::kKinkIndex0,0.);
    pt->AddCut(AliDielectronVarManager::kTPCnSigmaEle,-3.,3.);
    pt->AddCut(AliDielectronVarManager::kTPCnSigmaPio,3.5,1000.);
    pt->AddCut(AliDielectronVarManager::kTPCnSigmaPro,3.,1000.);
    pt->AddCut(AliDielectronVarManager::kNclsTPC,70.,160.);  
    pt->AddCut(AliDielectronVarManager::kTPCchi2Cl,0.,4.);
    pt->AddCut(AliDielectronVarManager::kImpactParXY,-1.,1.);
    pt->AddCut(AliDielectronVarManager::kImpactParZ,-3.,3.);
  }
  cuts->AddCut(pt);
  
	
	
	/*
	
  // vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv TRACK CUTS vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv /
  AliDielectronVarCuts *varAccCuts   = new AliDielectronVarCuts("acc","acc");
  varAccCuts->AddCut(AliDielectronVarManager::kPt,           0.8, 1e30);
  varAccCuts->AddCut(AliDielectronVarManager::kEta,         -0.9,  0.9);
  die->GetTrackFilter().AddCuts(varAccCuts);
  varAccCuts->Print();


  AliDielectronVarCuts   *varRecCuts = new AliDielectronVarCuts("VarRecCuts","VarRecCuts");
  varRecCuts->AddCut(AliDielectronVarManager::kNclsTPC,      70.,   160.);
  varRecCuts->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
  varRecCuts->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
  varRecCuts->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
  //  varRecCuts->AddCut(AliDielectronVarManager::kITSchi2Cl,     0.  ,   36.     ); // not defined in AOD
  varRecCuts->AddCut(AliDielectronVarManager::kKinkIndex0,.000001,1e30,kTRUE);

  AliDielectronTrackCuts *trkRecCuts = new AliDielectronTrackCuts("TrkRecCuts","TrkRecCuts");
  trkRecCuts->SetITSclusterCut(AliDielectronTrackCuts::kOneOf, 3); // SPD any
  trkRecCuts->SetRequireITSRefit(kTRUE);
  trkRecCuts->SetRequireTPCRefit(kTRUE);

  AliDielectronCutGroup  *grpRecCuts = new AliDielectronCutGroup("rec","rec",AliDielectronCutGroup::kCompAND);
  grpRecCuts->AddCut(trkRecCuts);
  grpRecCuts->AddCut(varRecCuts);
  die->GetTrackFilter().AddCuts(grpRecCuts);
  grpRecCuts->Print();


  /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv PID CUTS vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
  AliDielectronVarCuts *pidVarCuts = new AliDielectronVarCuts("varPIDCuts","varPIDCuts");
  pidVarCuts->AddCut(AliDielectronVarManager::kTPCnSigmaEle,  -2. ,    3.     ); //-3.0
  pidVarCuts->AddCut(AliDielectronVarManager::kTPCnSigmaPio,   3.5, 1000.     );
  pidVarCuts->AddCut(AliDielectronVarManager::kTPCnSigmaPro,   4. , 1000.     ); //3.0
  //  AliDielectronPID *pidCuts        = new AliDielectronPID("PIDCuts","PIDCuts"); //not used
  AliDielectronCutGroup *grpPIDCuts = new AliDielectronCutGroup("PID","PID",AliDielectronCutGroup::kCompAND);
  // grpPIDCuts->AddCut(pidVarCuts);
  //grpPIDCuts->AddCut(pidCuts);
  //  die->GetTrackFilter().AddCuts(grpPIDCuts);
  // grpPIDCuts->Print();

  //


  //exclude conversion electrons selected by the tender
  //   AliDielectronTrackCuts *noconv=new AliDielectronTrackCuts("noConv","noConv");
  //   noconv->SetV0DaughterCut(AliPID::kElectron,kTRUE);
  //   cuts->AddCut(noconv);

  //	cuts->Print();
}
void SetupV0Cuts(AliDielectron *die, Int_t cutDefinition)
{
	//
	// Setup the V0 cuts
	//
	AliDielectronV0Cuts *gammaV0Cuts = new AliDielectronV0Cuts("IsGamma","IsGamma");
	gammaV0Cuts->SetPdgCodes(22,11,11);
	gammaV0Cuts->SetDefaultPID(13); // TPC+-10, TOF+-3 if available
	gammaV0Cuts->AddCut(AliDielectronVarManager::kCosPointingAngle,TMath::Cos(0.02),1.0, kFALSE);
	gammaV0Cuts->AddCut(AliDielectronVarManager::kChi2NDF, 0.0,  10.0, kFALSE);//to be checked, if properly filled
	gammaV0Cuts->AddCut(AliDielectronVarManager::kLegDist, 0.0,   0.25, kFALSE);
	gammaV0Cuts->AddCut(AliDielectronVarManager::kR, 3.0,  90.0, kFALSE);
	gammaV0Cuts->AddCut(AliDielectronVarManager::kPsiPair, 0.0,   0.05, kFALSE);
	gammaV0Cuts->AddCut(AliDielectronVarManager::kM, 0.0,   0.05, kFALSE);
	// gammaV0Cuts->AddCut(AliDielectronVarManager::kOpeningAngle, 0.0,   0.1, kFALSE);
	gammaV0Cuts->AddCut(AliDielectronVarManager::kArmPt, 0.0,   0.05, kFALSE);
	// gammaV0Cuts->AddCut(AliDielectronVarManager::kArmAlpha, -0.35,0.35, kFALSE); // not sure if it works as expected
	gammaV0Cuts->SetExcludeTracks(kTRUE);//ktrue excludes tracks v0s, 
	//kfalse 
	// gammaV0Cuts->Print();
	//  const Double_t |cutAlphaG| < 0.35; &&  const Double_t cutQTG < 0.05;
	//  const Double_t |cutAlphaG2|[2] = {0.6, 0.8}; && const Double_t cutQTG2 < 0.04;
	die->GetTrackFilter().AddCuts(gammaV0Cuts);
	gammaV0Cuts->Print();	// 
}

void SetupV0add(AliDielectron *die, Int_t cutDefinition)
{
	//
	// Setup the V0 cuts
	//
	AliDielectronV0Cuts *gammaV0Add = new AliDielectronV0Cuts("IsGamma2","IsGamma2");
	gammaV0Add->SetPdgCodes(22,11,11);
	gammaV0Add->SetDefaultPID(13); // TPC+-10, TOF+-3 if available
	gammaV0Add->AddCut(AliDielectronVarManager::kCosPointingAngle,TMath::Cos(0.02),1.0, kFALSE);
	gammaV0Add->AddCut(AliDielectronVarManager::kChi2NDF, 0.0,  10.0, kFALSE);//to be checked, if properly filled
	gammaV0Add->AddCut(AliDielectronVarManager::kLegDist, 0.0,   0.25, kFALSE);
	gammaV0Add->AddCut(AliDielectronVarManager::kR, 3.0,  90.0, kFALSE);
	gammaV0Add->AddCut(AliDielectronVarManager::kPsiPair, 0.0,   0.05, kFALSE);
	gammaV0Add->AddCut(AliDielectronVarManager::kM, 0.0,   0.05, kFALSE);
	// gammaV0Cuts->AddCut(AliDielectronVarManager::kOpeningAngle, 0.0,   0.1, kFALSE);
	gammaV0Add->AddCut(AliDielectronVarManager::kArmPt, 0.0,   0.05, kFALSE);
	// gammaV0Cuts->AddCut(AliDielectronVarManager::kArmAlpha, -0.35,0.35, kFALSE); // not sure if it works as expected
	
	gammaV0Add->SetExcludeTracks(kFALSE);//ktrue excludes tracks v0s, 
	
	//kfalse 
	// gammaV0Cuts->Print();
	//  const Double_t |cutAlphaG| < 0.35; &&  const Double_t cutQTG < 0.05;
	//  const Double_t |cutAlphaG2|[2] = {0.6, 0.8}; && const Double_t cutQTG2 < 0.04;
	die->GetTrackFilter().AddCuts(gammaV0Add);
	gammaV0Add->Print();	// 
}
//______________________________________________________________________________________
void SetupPairCuts(AliDielectron *die, Int_t cutDefinition)
{
  //
  // Setup the pair cuts
  //

  // add conversion rejection
  AliDielectronVarCuts *gammaCut=new AliDielectronVarCuts("gammaCut","gammaCut");
  gammaCut->AddCut(AliDielectronVarManager::kM,0.,.05);
  die->GetPairPreFilter().AddCuts(gammaCut);

}

//______________________________________________________________________________________
void InitHistograms(AliDielectron *die, Int_t cutDefinition)
{
  //
  // Initialise the histograms
  //

  //Setup histogram Manager
  AliDielectronHistos *histos=
  new AliDielectronHistos(die->GetName(),
                          die->GetTitle());

  //Initialise histogram classes
  histos->SetReservedWords("Track;Pair");

  //Track classes
  //to fill also track info from 2nd event loop until 2
  for (Int_t i=0; i<2; ++i){
    histos->AddClass(Form("Track_%s",AliDielectron::TrackClassName(i)));
  }
  
  if(cutDefinition==0){
  //Pair classes
  // to fill also mixed event histograms loop until 10
  for (Int_t i=0; i<3; ++i){
    histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(i)));
  }
  }
  //add MC signal histograms to track and pair class
  if(die->GetMCSignals()) {
    for(Int_t isig=0; isig<die->GetMCSignals()->GetEntriesFast(); isig++) {
      TString sigMCname = die->GetMCSignals()->At(isig)->GetName(); 

      // mc truth
      histos->AddClass(Form("Pair_%s_MCtruth",       sigMCname.Data()));
      histos->AddClass(Form("Track_Legs_%s_MCtruth", sigMCname.Data())); 
      // mc reconstructed
      histos->AddClass(Form("Pair_%s",               sigMCname.Data()));
      histos->AddClass(Form("Track_Legs_%s",         sigMCname.Data())); 
    }
  }
  if(cutDefinition==0){

  //add histograms to event class
    histos->AddClass("Event");
    histos->UserHistogram("Event","","",300,-15.,15.,AliDielectronVarManager::kZvPrim);
    histos->UserHistogram("Event","","",
                          100,-2.,2.,100,-2.,2.,AliDielectronVarManager::kXvPrim,AliDielectronVarManager::kYvPrim);
    // histos->UserHistogram("Event","","",AliDielectronHelper::MakeArbitraryBinning(run_numbers),AliDielectronVarManager::kRunNumber);
    histos->UserProfile("Event","","",AliDielectronVarManager::kPairs, GetRunNumbers(), AliDielectronVarManager::kRunNumber);
    //    histos->UserHistogram("Events","","",100,0.,10.,AliDielectronVarManager::kPairs);
  }

  //add histograms to Track classes
  //  histos->UserProfile("Track","","",200,0,20.,AliDielectronVarManager::kPIn;AliDielectronHelper::MakeArbitraryBinning(run_numbers),AliDielectronVarManager::kRunNumber);
  //  histos->UserProfile("Track","","", AliDielectronVarManager::kTPCnSigmaEle,AliDielectronHelper::MakeArbitraryBinning(run_numbers),AliDielectronVarManager::kRunNumber);
  // histos->UserHistogram("Track","","",AliDielectronHelper::MakeLinBinning(4745, 121692, 126437),AliDielectronHelper::MakeLinBinning(100,-5.,+5.),AliDielectronVarManager::kRunNumber,AliDielectronVarManager::kTPCnSigmaEle);  
  histos->UserHistogram("Track","","",GetRunNumbers(),AliDielectronHelper::MakeLinBinning(100,-5.,+5.),AliDielectronVarManager::kRunNumber,AliDielectronVarManager::kTPCnSigmaEle);
  //  histos->UserProfile("Track","","", AliDielectronVarManager::kTOFnSigmaEle,AliDielectronHelper::MakeArbitraryBinning(run_numbers), AliDielectronVarManager::kRunNumber);
  histos->UserHistogram("Track","","",GetRunNumbers(),AliDielectronHelper::MakeLinBinning(100,-5.,+5.),AliDielectronVarManager::kRunNumber,AliDielectronVarManager::kTOFnSigmaEle);
   //histos->UserHistogram("Track","","", run_numbers,AliDielectronHelper::MakeLinBinning(100,-5.,+5.),
   //                   AliDielectronVarManager::kRunNumber,AliDielectronVarManager::kTOFnSigmaEle);

  //  histos->UserProfile("Track","","", AliDielectronVarManager::kITSnSigmaEle,AliDielectronHelper::MakeArbitraryBinning(run_numbers), AliDielectronVarManager::kRunNumber);
   histos->UserHistogram("Track","","",GetRunNumbers(),AliDielectronHelper::MakeLinBinning(100,-5.,+5.),AliDielectronVarManager::kRunNumber,AliDielectronVarManager::kITSnSigmaEle);
  
   //histos->UserProfile("Track","","", AliDielectronVarManager::kNclsTRD,AliDielectronHelper::MakeArbitraryBinning(run_numbers), AliDielectronVarManager::kRunNumber);
   histos->UserHistogram("Track","","",GetRunNumbers(),AliDielectronHelper::MakeLinBinning(300,0.,+300.),AliDielectronVarManager::kRunNumber,AliDielectronVarManager::kNclsTRD);

   //histos->UserProfile("Track","","",AliDielectronVarManager::kEta,AliDielectronHelper::MakeArbitraryBinning(run_numbers),AliDielectronVarManager::kRunNumber);
   histos->UserHistogram("Track","","",GetRunNumbers(),AliDielectronHelper::MakeLinBinning(20,-1.0,1.0),AliDielectronVarManager::kRunNumber,AliDielectronVarManager::kEta);

   //histos->UserProfile("Track","","", AliDielectronVarManager::kPhi,AliDielectronHelper::MakeArbitraryBinning(run_numbers),AliDielectronVarManager::kRunNumber);
   histos->UserHistogram("Track","","",GetRunNumbers(),AliDielectronHelper::MakeLinBinning(63,0.,6.32),AliDielectronVarManager::kRunNumber,AliDielectronVarManager::kPhi);

   histos->UserProfile("Track","","", AliDielectronVarManager::kNclsTPC, 200.,0.,20., AliDielectronVarManager::kPt);
   histos->UserHistogram("Track","","", 200.,0.,200.,AliDielectronVarManager::kNclsTPC, 200.,0.,20., AliDielectronVarManager::kPt);
   histos->UserHistogram("Track","","",AliDielectronHelper::MakeLinBinning(300,0.5,300.5),AliDielectronVarManager::kNclsTRD);
   histos->UserHistogram("Track","","",AliDielectronHelper::MakeLinBinning(150,0.,150.),AliDielectronVarManager::kNaccTrckltsEsd10Corr);
   histos->UserHistogram("Track","","",AliDielectronHelper::MakeLinBinning(5,-0.5,2.0),AliDielectronVarManager::kTOFPIDBit);
   //  histos->UserHistogram("Track","","",GetRunNumbers(),
   //	        400,0.2,20.,200,0.,200.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCsignal,kTRUE);
   histos->UserHistogram("Track","","",
			 100,0.2,20.,100,-10.,10.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kITSnSigmaEle);
  histos->UserHistogram("Track","","",
                        100,0.2,20.,100,-10.,10.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,kTRUE);
  histos->UserHistogram("Track","","",
                        100,0.2,20.,100,-10.,10.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaPio,kTRUE);
  histos->UserHistogram("Track","","",
                        100,0.2,20.,100,-10.,10.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaPro,kTRUE);
  histos->UserHistogram("Track","","",
			200,0.2,20.,100,-15.,15.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaKao,kTRUE);    	
  histos->UserHistogram("Track","","",
			200,0.2,20.,100,0.,1.2,AliDielectronVarManager::kP,AliDielectronVarManager::kTOFbeta,kTRUE);
  histos->UserHistogram("Track","","",
			200,0.2,20.,200,-20.,20.,AliDielectronVarManager::kP,AliDielectronVarManager::kTOFnSigmaEle,kTRUE);    
  histos->UserHistogram("Track","","",
			144,0.0,6.285,100,0.0,200,AliDielectronVarManager::kPhi,AliDielectronVarManager::kTPCsignal);
  histos->UserHistogram("Track","","",
			40,-1.0,1.0,100,0.0,200,AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCsignal);
  histos->UserHistogram("Track","","",
			100,-5,5.,100,-5.,5.,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kITSnSigmaEle);
  histos->UserHistogram("Track","","",
			100,-5,5.,100,-5.,5.,AliDielectronVarManager::kTOFnSigmaEle,AliDielectronVarManager::kTPCnSigmaEle);
  //histos->UserHistogram("Track","","",
  //                    100,-2,2,144,0,6.285,AliDielectronVarManager::kEta,AliDielectronVarManager::kPhi);
  histos->UserHistogram("Track","","",
                        160,-0.5,159.5,AliDielectronVarManager::kNclsTPC);
  histos->UserHistogram("Track","","",
                        100,0.,10.,AliDielectronVarManager::kTPCchi2Cl);
  histos->UserHistogram("Track","","",
                        150,-15,15,160,-0.5,159.5,AliDielectronVarManager::kYsignedIn,AliDielectronVarManager::kTPCsignalN);
   //rjim findable cluster vs pt
  histos->UserHistogram("Track","","",
			200.,0.0.,20.0, 161,-0.5,161.5,AliDielectronVarManager::kPt,AliDielectronVarManager::kNFclsTPCr);
  //rjim frac find vs pt
  histos->UserHistogram("Track","","",
			200.0,0.0,20.0,160,0,1.1,AliDielectronVarManager::kPt,AliDielectronVarManager::kNFclsTPCrFrac);
  
  // histos rjimenez 2nd part
  //check tofbit
  
  histos->UserHistogram("Track","","",
			3.,-0.5,2.5, 200,0.,20.,AliDielectronVarManager::kTOFPIDBit,AliDielectronVarManager::kPt);
  
  histos->UserHistogram("Track","","",
			160,-0.5,159.5,AliDielectronVarManager::kNFclsTPCr);
  histos->UserHistogram("Track","","",
			160,0,1.1,AliDielectronVarManager::kNFclsTPCrFrac);
  
  //rjimenez nsigma vseta
  histos->UserHistogram("Track","","",
			20,-1.0,1.0,100,-10.,10.,AliDielectronVarManager::kEta,AliDielectronVarManager::kITSnSigmaEle);
  
  histos->UserHistogram("Track","","",
			20,-1.0,1.0,100,-10.,10.,AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEle);
  
  histos->UserHistogram("Track","","",
			20,-1.0,1.0,100,-10.,10.,AliDielectronVarManager::kEta,AliDielectronVarManager::kTOFnSigmaEle); 
  //rjimenez nsigma vs phi
  histos->UserHistogram("Track","","",
			20,0.0,7.0,100,-10.,10.,AliDielectronVarManager::kPhi,AliDielectronVarManager::kITSnSigmaEle);
  
  histos->UserHistogram("Track","","",
			20,0.0,7.0,100,-10.,10.,AliDielectronVarManager::kPhi,AliDielectronVarManager::kTPCnSigmaEle);
  
  histos->UserHistogram("Track","","",
			20,0.0,7.0,100,-10.,10.,AliDielectronVarManager::kPhi,AliDielectronVarManager::kTOFnSigmaEle);
  
  //for now SPD tracklets but multiplicity should be implemented
  histos->UserHistogram("Track","","",
			100,0.0,100.,100,-10.,10.,AliDielectronVarManager::kNaccTrckltsEsd10Corr,AliDielectronVarManager::kTPCnSigmaEle,kTRUE);
  histos->UserHistogram("Track","","",
			100,0.0,100.,200,-20.,20.,AliDielectronVarManager::kNaccTrckltsEsd10Corr,AliDielectronVarManager::kTOFnSigmaEle,kTRUE);      
  //
	
  histos->UserHistogram("Track","","",
			200,-20.,20.,100,-10.,10.,AliDielectronVarManager::kYsignedIn,AliDielectronVarManager::kTPCnSigmaEle);
  
  histos->UserHistogram("Track","","",
			200,-20.,20.,200,0.2,20.,AliDielectronVarManager::kYsignedIn,AliDielectronVarManager::kPIn);

  //inner and outer read out TPC clusters
  histos->UserHistogram("Track","","",
			70,0.0,7.0,160.,0.0,160.0,AliDielectronVarManager::kPhi,AliDielectronVarManager::kTPCclsIRO);
  histos->UserHistogram("Track","","",
			70,0.0,7.0,160.,0.0,160.0,AliDielectronVarManager::kPhi,AliDielectronVarManager::kTPCclsORO);
  histos->UserHistogram("Track","","",
			10000,0.0,1.0,AliDielectronVarManager::kM);
  histos->UserHistogram("Track","","",
			200,0.0,20.0,10000,0.0,1.0,AliDielectronVarManager::kPt,AliDielectronVarManager::kM);
  histos->UserHistogram("Track","","",
			11,-0.5,10.5.,AliDielectronVarManager::kTRDntracklets);
  
  histos->UserHistogram("Track","","",
			200.,0.0.,20.0,11,-0.5,10.5,AliDielectronVarManager::kPt,AliDielectronVarManager::kTRDntracklets);
  
  histos->UserHistogram("Track","","",
			150.,-15.0.,15.0,20,0.,1.,AliDielectronVarManager::kYsignedIn,AliDielectronVarManager::kTRDprobEle);
  
  histos->UserHistogram("Track","","",
			150.,-15.0,15.0,20,0.,1.,AliDielectronVarManager::kYsignedIn,AliDielectronVarManager::kTRDprobPio);
  
  histos->UserHistogram("Track","","",
			200.,0.2,20.0,20,0.,1.,AliDielectronVarManager::kP,AliDielectronVarManager::kTRDprob2DEle);
  
  histos->UserHistogram("Track","","",
			200.,0.2,20.0,20,0.,1.,AliDielectronVarManager::kP,AliDielectronVarManager::kTRDprob2DPio);
  
  histos->UserHistogram("Track","","",
			8.,-0.5,7.5,AliDielectronVarManager::kTRDpidQuality);
  
  histos->UserHistogram("Track","","",
			100.,0.0,1.0,AliDielectronVarManager::kTRDpidEffLeg);
  
  histos->UserHistogram("Track","","",
			20,-1.0,1.0,70,-3.5,3.5,AliDielectronVarManager::kEta,AliDielectronVarManager::kTRDphi);
  
  histos->UserHistogram("Track","","",500,-1.,1.,AliDielectronVarManager::kImpactParXY);
  histos->UserHistogram("Track","","",600,-3.,3.,AliDielectronVarManager::kImpactParZ);

  //add histograms to Pair classes

  if(cutDefinition ==0){
  histos->UserHistogram("Pair","","",
                        301,-.01,6.01,AliDielectronVarManager::kM);
  histos->UserHistogram("Pair","","",
                        100,-1.,1.,AliDielectronVarManager::kY);
  histos->UserHistogram("Pair","","",
                        100,0.,3.15,AliDielectronVarManager::kOpeningAngle);
  histos->UserHistogram("Pair","","",
                        301,-.01,6.01,100,0.,3.15,AliDielectronVarManager::kM,AliDielectronVarManager::kOpeningAngle);
  histos->UserHistogram("Pair","","",
                        301,-.01,6.01,200,0.2,20.,AliDielectronVarManager::kM,AliDielectronVarManager::kPt);

  }
  die->SetHistogramManager(histos);

}

//______________________________________________________________________________________
void InitCF(AliDielectron* die, Int_t cutDefinition)
{
  //
  // Setupd the CF Manager if needed
  //

  AliDielectronCF *cf=new AliDielectronCF(die->GetName(),die->GetTitle());

  //pair variables
  cf->AddVariable(AliDielectronVarManager::kPt,"0.0, 1.0, 1.3, 2.0, 3.0, 5., 7.0, 10.0, 100.0");

  cf->AddVariable(AliDielectronVarManager::kZvPrim,"-10,-5,0,5,10");
  cf->AddVariable(AliDielectronVarManager::kY,"-1,-0.9,-0.8,-0.3,0,0.3,0.9,1.0");
  cf->AddVariable(AliDielectronVarManager::kM,125,0.,125*.04); //40Mev Steps
  cf->AddVariable(AliDielectronVarManager::kPairType,11,0,11);

  //leg variables
  cf->AddVariable(AliDielectronVarManager::kPt,"0.0, 0.8, 1.0, 1.1, 1.2, 1.3, 100.0",kTRUE);
  cf->AddVariable(AliDielectronVarManager::kEta,"-1.,-0.9,-0.8,0,0.8,0.9,1.0",kTRUE);

  //cf->AddVariable(AliDielectronVarManager::kITSLayerFirstCls,6,0,6,kTRUE);

  if (hasMC && 0){ //ATTENTION SWITCHED OFF
    cf->AddVariable(AliDielectronVarManager::kPdgCode,10000,-5000.5,4999.5,kTRUE);
    cf->AddVariable(AliDielectronVarManager::kPdgCodeMother,10000,-5000.5,4999.5,kTRUE);
    cf->AddVariable(AliDielectronVarManager::kPdgCodeGrandMother,10000,-5000.5,4999.5,kTRUE);
  }

  if(hasMC) {
    //only steps for efficiencies
    cf->SetStepsForMCtruthOnly();
    //only in this case write MC truth info
    if (cutDefinition==0){
      cf->SetStepForMCtruth();
    }
  }
  // cf->SetStepsForSignal();

  die->SetCFManagerPair(cf);
}

//______________________________________________________________________________________
void SetupMCsignals(AliDielectron *die){
  AliDielectronSignalMC* promptJpsi = new AliDielectronSignalMC("promptJpsi","Prompt J/psi");   // prompt J/psi (not from beauty decays)
  promptJpsi->SetLegPDGs(11,-11);
  promptJpsi->SetMotherPDGs(443,443);
  promptJpsi->SetGrandMotherPDGs(503,503,kTRUE,kTRUE);   // not from beauty hadrons
  promptJpsi->SetMothersRelation(AliDielectronSignalMC::kSame);
  promptJpsi->SetFillPureMCStep(kTRUE);
  promptJpsi->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  promptJpsi->SetCheckBothChargesLegs(kTRUE,kTRUE);
  promptJpsi->SetCheckBothChargesMothers(kTRUE,kTRUE);
  promptJpsi->SetCheckBothChargesGrandMothers(kTRUE,kTRUE);
  die->AddSignalMC(promptJpsi);

  AliDielectronSignalMC* promptJpsiNonRad = new AliDielectronSignalMC("promptJpsiNonRad","Prompt J/psi non-Radiative");   // prompt J/psi (not from beauty decays)
  promptJpsiNonRad->SetLegPDGs(11,-11);
  promptJpsiNonRad->SetMotherPDGs(443,443);
  promptJpsiNonRad->SetGrandMotherPDGs(503,503,kTRUE,kTRUE);   // not from beauty hadrons
  promptJpsiNonRad->SetMothersRelation(AliDielectronSignalMC::kSame);
  promptJpsiNonRad->SetFillPureMCStep(kTRUE);
  promptJpsiNonRad->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  promptJpsiNonRad->SetCheckBothChargesLegs(kTRUE,kTRUE);
  promptJpsiNonRad->SetCheckBothChargesMothers(kTRUE,kTRUE);
  promptJpsiNonRad->SetCheckBothChargesGrandMothers(kTRUE,kTRUE);
  promptJpsiNonRad->SetJpsiRadiative(AliDielectronSignalMC::kIsNotRadiative);
  die->AddSignalMC(promptJpsiNonRad);

  AliDielectronSignalMC* promptJpsiRad = new AliDielectronSignalMC("promptJpsiRad","Prompt J/psi Radiative");   // prompt J/psi (not from beauty decays)
  promptJpsiRad->SetLegPDGs(11,-11);
  promptJpsiRad->SetMotherPDGs(443,443);
  promptJpsiRad->SetGrandMotherPDGs(503,503,kTRUE,kTRUE);   // not from beauty hadrons
  promptJpsiRad->SetMothersRelation(AliDielectronSignalMC::kSame);
  promptJpsiRad->SetFillPureMCStep(kTRUE);
  promptJpsiRad->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  promptJpsiRad->SetCheckBothChargesLegs(kTRUE,kTRUE);
  promptJpsiRad->SetCheckBothChargesMothers(kTRUE,kTRUE);
  promptJpsiRad->SetCheckBothChargesGrandMothers(kTRUE,kTRUE);
  promptJpsiRad->SetJpsiRadiative(AliDielectronSignalMC::kIsRadiative);
  die->AddSignalMC(promptJpsiRad);

}


