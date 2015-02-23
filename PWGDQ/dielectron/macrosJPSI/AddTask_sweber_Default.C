
void InitHistograms(AliDielectron *die, Int_t cutDefinition);
void SetupEventCuts(AliDielectron *die, Int_t cutDefinition);
void SetupTrackCuts(AliDielectron *die, Int_t cutDefinition);
void SetupPairCuts(AliDielectron *die,  Int_t cutDefinition);
void SetupV0cuts(AliDielectron *die, Int_t cutDefinition);

enum  {kDefault,kFirst, kStrict,kLoose};


#if 0
void InitCF(AliDielectron* die, Int_t cutDefinition);
void SetupMCsignals(AliDielectron *die);
void SetupV0Add(AliDielectron *die, Int_t cutDefinition);

#endif


#if 0
TVectorD *GetRunNumbers() {
    if (period="LHC10d"){
    return AliDielectronHelper::MakeArbitraryBinning ("122374,122375,124751,125023,125085,125097,125100,125101,125134,125296,125628,125630,125632,125633,125842,125843,125844,125847,125848,125849,125850,125851,125855,125941,125942,125943,125944,125945,125946,125947,125948,125949,125950,125951,125952,125953,125954,125955,125956,125957,125958,125959,125960,125961,125962,125963,125964,125965,125966,125969,125970,125976,125978,125981,125982,125983,125984,125985,125986,125997,125998,125999,126004,126007,126008,126073,126078,126081,126082,126088,126090,126097,126158,126160,126168,126283,126284,126285,126351,126352,126359,126403,126404,126405,126406,126407,126408,126409,126422,126424,126425,126432");
     }
}
#endif


  // increasing order of run array is important

//TString names=("default");
//TString names=("default;SPDfirst;Strict;TOF");
TString names=("default;SPDfirst;Strict;Loose");
TObjArray *arrNames=names.Tokenize(";");

const Int_t nDie=arrNames->GetEntries();
Bool_t isAOD=kFALSE;
Bool_t hasMC=kFALSE;
Int_t iPeriod=-1;
enum { k10b=0, k10c, k10d, k10e, k10f, k10h, k11a, k11d, k11h, k12h, k13b, k13c, k13d, k13e, k13f };

//______________________________________________________________________________________
AliAnalysisTask* AddTask_sweber_Default(TString prod="", Bool_t isMC=kFALSE)
{
  //get the current analysis manager
  
  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_sweber_Default", "No analysis manager found.");
    return 0;
  }

  //Do we have an MC handler?
  hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);

  //AOD input?
  isAOD=mgr->GetInputEventHandler()->IsA()==AliAODInputHandler::Class();
  if(isAOD) hasMC=isMC;

  //Get the current train configuration
  TString trainConfig=gSystem->Getenv("CONFIG_FILE");
  TString list=gSystem->Getenv("LIST");
  if( list.IsNull()) list=prod;

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

  
	AliDielectronEventCuts *eventCuts=new AliDielectronEventCuts("eventCuts","Vertex Track && |vtxZ|<10 && ncontrib>0");

  
	eventCuts->SetRequireVertex();
	eventCuts->SetMinVtxContributors(1);
	eventCuts->SetVertexZ(-10.,10.);
#if 0
  
//	eventCuts->SetRunRejection( AliDielectronHelper::MakeArbitraryBinning("114737,114740,114743,114744,114745,114746,114747,114750,114751,114753,114757,114778,114783,114785,114916,114917,114919,115165,115173,115237,115315,115325,115338,115413,115514,115516,115880,115881,115882,115887,115888,115889,115890,115892,115987,116111,116112,116118,116123,116130,116134,116197,116198,116203,116204,116287,116401,116559,116561,116572,116609,116610,116611,116640,116642,116644,116681,116682,116684,116948,117033,117034,117035,117039,117041,117042,117045,117046,117049,117051,117054,117061,117064,117065,117077,117078,117082,117086,117098,117100,117110,117113,117117,117118,117119,117120,117121,117221,117223,118359,118503,118504,118557,118742,118903,119022,119033,119034,119035,119037,119038,119039,119040,119041,119042,119043,119045,119047,119048,119055,119057,119061,119067,119068,119077,119079,119084,119085,119086,119156,119158,119160,119162,119164,119837,119838,119840,119843,119854,119865,119902,119903,119904,119907,119909,119912,119913,119915,119917,119923,119924,119926,119932,119934,119935,119941,119946,119948,119952,119954,119959,119961,119963,119965,119967,119969,119970,119971,119973,119978,119979,119994,119998,119999,120000,120001,120002,120003,120004,120064,120065,120066,120068,120071,120075,120078,120241,120242,120243,120275,120281,120476,120479,120480,120481,120482,120483,120485,120486,120487,120488,120489,120490,120492,120493,120494,120495,120496,120498,120500,120611,120612,120613,120614,120615,120618,120619,120622,120624,120625,120670,120672,120673,120675,120677,120678,120679,120681,120682,120683,120685,120687,120688,120689,120691,120692,120693,120742,120743,120744,120746,120818,120819,120828,120967,121694,121720,122195,122372,122373,124183,124186,124187,124191,124355,124358,124359,124360,124362,124364,124367,124371,124374,124378,124380,124381,124383,124385,124388,124400,124600,124602,124603,124604,124605,124606,124607,124608,124693,124702,124712,124713,124714,124738,124739,124743,124744,124745,124746,124750,124850,124886,125022,125076,125077,125078,125131,125292,125294,125634,125841,125941,125942,125943,125944,125945,125946,125947,125948,125949,125950,125951,125952,125953,125954,125955,125956,125957,125958,125959,125960,125961,125962,125963,125964,125965,125966,125969,125970,125976,125978,125981,125982,125983,125984,125985,125986,125997,125998,125999,126086,126087,126162,126169,126170,126171,126172,126173,126174,126175,126177,126411,126416,126419,126420,126421,126437,127102,127105,127112,127115,127123,127129,127133,127138,127712,127714,127715,127718,127719,127723,127724,127729,127730,127813,127814,127815,127817,127819,127820,127822,127930,127931,127932,127933,127935,127936,127937,127940,127941,127942,128050,128053,128175,128180,128182,128183,128185,128186,128189,128190,128191,128192,128256,128257,128260,128262,128263,128357,128358,128359,128361,128363,128364,128367,128369,128371,128372,128373,128451,128453,128454,128455,128456,128457,128458,128459,128463,128468,128469,128470,128471,128472,128473,128474,128475,128483,128507,128581,128589,128610,128776,128813,128814,128817,128818,128849,128910,128911,128912,129041,129505,129506,129508,129509,129510,129511,129518,129522,129526,129529,129541,129589,129597,129598,129648,129649,129654,129655,129665,129667,129731,129745,129747,129748,129750,129760,129763,130148,130156,130170,130179,130348,130353,130365,130369,130529,130612,130619,130627,130640,130694,130831,130833,133004,133005,133282,133328,133415,133417,133418,133419,133671,133672,133673,133674,133675,133680,133761,133924,133979,133981,133985,134094,134199,134200,134201,134202,134203,134204,134298,134299,134301,134302,134303,134304,134305,134306,134489,134497,134657,134660,134666,134667,134670,134671,134672,134674,134679,134681,134685,134690,134691,134692,134778,134779,134780,134781,134840,134841,134844,134905,134907,134908,134909,134910,134911,134912,134913,134914,134916,134919,134920,134921,134922,134924,134925,134926,134928,134929,134930,134931,135025,135029,135031")  );
#endif
  //eventCuts->Print();
	task->SetEventFilter(eventCuts);
  
  
  
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
    mgr->CreateContainer("sweber_Default_tree",
                         TTree::Class(),
                         AliAnalysisManager::kExchangeContainer,
                         "sweber_Default_default");
  
  AliAnalysisDataContainer *cOutputHist1 =
    mgr->CreateContainer("sweber_Default_QA",
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         "JPSI.root");

  AliAnalysisDataContainer *cOutputHist2 =
    mgr->CreateContainer("sweber_Default_CF",
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         "JPSI.root");
	
  
  AliAnalysisDataContainer *cOutputHist3 =
    mgr->CreateContainer("sweber_Default_EventStat",
                         TH1D::Class(),
                         AliAnalysisManager::kOutputContainer,
                         "JPSI.root");
  
	
	
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
  
  SetupV0Cuts(die,cutDefinition);

  SetupPairCuts(die,cutDefinition);

  /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv MISC vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
  // Monte Carlo Signals
  //if (hasMC) SetupMCsignals(die);
  // prefilter settings
  // die->SetPreFilterUnlikeOnly();//  die->SetNoPairing();//  die->SetPreFilterAllSigns();
  // cut QA
  // die->SetCutQA();

  /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv OUTPUT vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
  InitHistograms(die,cutDefinition);
  //InitCF(die,cutDefinition);


     AliDielectronMixingHandler *mix=new AliDielectronMixingHandler;
     mix->AddVariable(AliDielectronVarManager::kZvPrim,"-10,-8,-5,0,5,8,10");
     mix->AddVariable(AliDielectronVarManager::kPhi ,"0,3.14,6.3");
     mix->SetDepth(10);
    die->SetMixingHandler(mix);
  //
    AliDielectronTrackRotator *rot=new AliDielectronTrackRotator;
   rot->SetConeAnglePhi(TMath::Pi()/180.*135.);
  rot->SetIterations(20);
  die->SetTrackRotator(rot);
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
	
	if (cutDefinition==kDefault|| cutDefinition==kLoose){
		refit->SetRequireITSRefit(kTRUE);
		refit->SetRequireTPCRefit(kTRUE);
		refit->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
	}else if( cutDefinition==kFirst  || cutDefinition==kStrict ){
		refit->SetRequireITSRefit(kTRUE);
		refit->SetRequireTPCRefit(kTRUE);
		refit->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
	}
	cuts->AddCut(refit);
	
	
	
	AliDielectronVarCuts *pt = new AliDielectronVarCuts("ptCut","pt cut");
	
	
	
	if ( cutDefinition==kDefault || cutDefinition==kFirst ){	
		pt->AddCut(AliDielectronVarManager::kPt,1.,1.e30);
		pt->AddCut(AliDielectronVarManager::kEta,-0.9,0.9);
		pt->AddCut(AliDielectronVarManager::kKinkIndex0,0.);	
		//rjim NTPCclusters cut
		pt->AddCut(AliDielectronVarManager::kNclsTPC,70.,160.);  
		pt->AddCut(AliDielectronVarManager::kTPCchi2Cl,0.,4.);
		//impact parameter
		pt->AddCut(AliDielectronVarManager::kImpactParXY,-1.,1.);
		pt->AddCut(AliDielectronVarManager::kImpactParZ,-3.,3.);
		
		
		pt->AddCut(AliDielectronVarManager::kTPCnSigmaEle,-3.,3.);
		pt->AddCut(AliDielectronVarManager::kTPCnSigmaPio,3.5,1000.);
		pt->AddCut(AliDielectronVarManager::kTPCnSigmaPro,3.,1000.);	
	}
	
	else if( cutDefinition==kLoose ){
		
		
		pt->AddCut(AliDielectronVarManager::kPt,.5 ,1.e30);
		pt->AddCut(AliDielectronVarManager::kEta,-0.9,0.9);
		pt->AddCut(AliDielectronVarManager::kKinkIndex0,0.);	
		//rjim NTPCclusters cut
		pt->AddCut(AliDielectronVarManager::kNclsTPC,70.,160.);  
		pt->AddCut(AliDielectronVarManager::kTPCchi2Cl,0.,4.);
		//impact parameter
		pt->AddCut(AliDielectronVarManager::kImpactParXY,-1.,1.);
		pt->AddCut(AliDielectronVarManager::kImpactParZ,-3.,3.);
		
		pt->AddCut(AliDielectronVarManager::kTPCnSigmaEle,-3.,3.);
		pt->AddCut(AliDielectronVarManager::kTPCnSigmaPio,4.,1000.);

    	AliDielectronCutGroup *pTOFproton = new AliDielectronCutGroup("pTOFproton", "pTOFproton", AliDielectronCutGroup::kCompOR);//NOTE taken out for testing
  //1. electron candidates inside proton band and with electron TOFpid, rejected if outside 3sigma TOF
    	AliDielectronVarCuts *pidTOF = new AliDielectronVarCuts("PIDTOFthere ","TOF nSigma |e|<3 + TPC nSigma |p|<3.5");
  //TOF bit already internally required
    	pidTOF->AddCut(AliDielectronVarManager::kTOFnSigmaEle,-3.,3.);
		pidTOF->AddCut(AliDielectronVarManager::kTPCnSigmaPro,-3.5,3.5);
    	pTOFproton->AddCut(pidTOF);

    //2. electron candidates outside proton band accepted if no matching inside with TOF
    	AliDielectronVarCuts *pidproton = new AliDielectronVarCuts("PIDTPCproton ","TPC nSigma |P|>3.5");
    	pidproton->AddCut(AliDielectronVarManager::kTPCnSigmaPro,-3.5,3.5,kTRUE);
    	pTOFproton->AddCut(pidproton);
		cuts->AddCut(pTOFproton);
		
		
		
		
		
		// same for kaons
		
    	AliDielectronCutGroup *pTOFkaon = new AliDielectronCutGroup("pTOFkaon", "pTOFkaon", AliDielectronCutGroup::kCompOR);
    	AliDielectronVarCuts *pidTOF2 = new AliDielectronVarCuts("PIDTOFthere2 ","TOF nSigma |e|<3 + TPC nSigma |k|<3.5");
    	pidTOF2->AddCut(AliDielectronVarManager::kTOFnSigmaEle,-3.,3.);
		pidTOF2->AddCut(AliDielectronVarManager::kTPCnSigmaKao,-3.5,3.5);
    	pTOFkaon->AddCut(pidTOF2);
    	AliDielectronVarCuts *pidkaon = new AliDielectronVarCuts("pidkaon ","TPC nSigma |k|>3.5");
    	pidkaon->AddCut(AliDielectronVarManager::kTPCnSigmaKao,-3.5,3.5,kTRUE);
    	pTOFkaon->AddCut(pidkaon);
	
		cuts->AddCut(pTOFkaon);
		
		/*
		// same for pions
		
    	AliDielectronCutGroup *pTOFpion = new AliDielectronCutGroup("pTOFpion", "pTOFpion", AliDielectronCutGroup::kCompOR);
    	AliDielectronVarCuts *pidTOF3 = new AliDielectronVarCuts("PIDTOFthere3 ","TOF nSigma |e|<3 + TPC nSigma |pi|<3.5");
    	pidTOF3->AddCut(AliDielectronVarManager::kTOFnSigmaEle,-3.,3.);
		pidTOF3->AddCut(AliDielectronVarManager::kTPCnSigmaPio,-3.5,3.5);
    	pTOFpion->AddCut(pidTOF3);
    	AliDielectronVarCuts *pidpion = new AliDielectronVarCuts("pidpion ","TPC nSigma |pi|>3.5");
    	pidpion->AddCut(AliDielectronVarManager::kTPCnSigmaPio,-3.5,3.5,kTRUE);
    	pTOFpion->AddCut(pidpion);
	
		cuts->AddCut(pTOFpion);

		*/
		
	}
	else if( cutDefinition==kStrict ){
		
		
		pt->AddCut(AliDielectronVarManager::kPt,.8 ,1.e30);
		pt->AddCut(AliDielectronVarManager::kEta,-0.9,0.9);
		pt->AddCut(AliDielectronVarManager::kKinkIndex0,0.);	
		//rjim NTPCclusters cut
		pt->AddCut(AliDielectronVarManager::kNclsTPC,70.,160.);  
		pt->AddCut(AliDielectronVarManager::kTPCchi2Cl,0.,4.);
		//impact parameter
		pt->AddCut(AliDielectronVarManager::kImpactParXY,-1.,1.);
		pt->AddCut(AliDielectronVarManager::kImpactParZ,-3.,3.);
		
		pt->AddCut(AliDielectronVarManager::kTPCnSigmaEle,-3.,3.);
		pt->AddCut(AliDielectronVarManager::kTPCnSigmaPio,4.,1000.);
		pt->AddCut(AliDielectronVarManager::kTPCnSigmaPro,3.5,1000.);
		
		
	}
	
	cuts->AddCut(pt);



  //exclude conversion electrons selected by the tender
  //   AliDielectronTrackCuts *noconv=new AliDielectronTrackCuts("noConv","noConv");
  //   noconv->SetV0DaughterCut(AliPID::kElectron,kTRUE);
  //   cuts->AddCut(noconv);

	cuts->Print();
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
	//gammaV0Cuts->Print();	// 
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

  //Pair classes
  // to fill also mixed event histograms loop until 10
  Int_t pairClasses [5] = {0,1,2,4, 10};
  
  
  for (Int_t i=0; i<5; ++i){
    histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(pairClasses[i])));
  }

  //add histograms to event class
    histos->AddClass("Event");
  histos->UserHistogram("Event","nEvents","Number of processed events after cuts;Number events",
                        4,0.,4.,AliDielectronVarManager::kNevents);
    histos->UserHistogram("Event","VtxZ","Vertex Z;Z[cm]",300,-15.,15.,AliDielectronVarManager::kZvPrim);
    histos->UserHistogram("Event","VtxX_VtxY","Vertex Y vs. X;X (cm);Y (cm)",
                          100,-2.,2.,100,-2.,2.,AliDielectronVarManager::kXvPrim,AliDielectronVarManager::kYvPrim);
    // }
    //rjimenez tracks/event
    histos->UserHistogram("Event","Ntrk","tracks/event; Ntracks; #Entries",
			  700,0.0,700.,AliDielectronVarManager::kNTrk);
  
  
  
  histos->UserHistogram("Track","Pt","Pt;Pt [GeV];#tracks",200,0,20.,AliDielectronVarManager::kPt);
  histos->UserHistogram("Track","dEdx_P","dEdx;P [GeV];TPC signal (arb units);#tracks",
                        400,0.2,20.,200,0.,200.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCsignal,kTRUE);
  
  histos->UserHistogram("Track","TPCnSigmaEle_P","TPC number of sigmas Electrons;P [GeV];TPC number of sigmas Electrons",
                        100,0.2,20.,100,-10.,10.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,kTRUE);
  histos->UserHistogram("Track","TPCnSigmaPio_P","TPC number of sigmas Pions;P [GeV];TPC number of sigmas Pions",
                        100,0.2,20.,100,-15.,15.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaPio,kTRUE);
  histos->UserHistogram("Track","TPCnSigmaPro_P","TPC number of sigmas Protons;P [GeV];TPC number of sigmas Protons",
                        100,0.2,20.,100,-15.,15.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaPro,kTRUE);
  histos->UserHistogram("Track","TPCnSigmaKao_P","TPC number of sigmas Kaons;P [GeV];TPC number of sigmas Kaons",
                        100,0.2,20.,100,-15.,15.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaKao,kTRUE);
  
  histos->UserHistogram("Track","TOFnSigmaEle_P","TOF number of sigmas Electrons;P [GeV];TOF number of sigmas Electrons",
                        100,0.2,20.,100,-10.,10.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaEle,kTRUE);
  histos->UserHistogram("Track","TOFnSigmaPio_P","TOF number of sigmas Pions;P [GeV];TOF number of sigmas Pions",
                        100,0.2,20.,100,-15.,15.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaPio,kTRUE);
  histos->UserHistogram("Track","TOFnSigmaPro_P","TOF number of sigmas Protons;P [GeV];TOF number of sigmas Protons",
                        100,0.2,20.,100,-15.,15.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaPro,kTRUE);
  histos->UserHistogram("Track","TOFnSigmaKao_P","TOF number of sigmas Kaons;P [GeV];TOF number of sigmas Kaons",
                        100,0.2,20.,100,-15.,15.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaKao,kTRUE);

  
	histos->UserHistogram("Track","Eta","Eta;Eta;#tracks",200,-2,2.,AliDielectronVarManager::kEta);
  histos->UserHistogram("Track","TPCnCls","Number of Clusters TPC;TPC number clusteres;#tracks",162,-1,161,AliDielectronVarManager::kNclsTPC);
  										
  histos->UserHistogram("Track","dXY","dXY;dXY [cm];#tracks",200,-2.,2.,AliDielectronVarManager::kImpactParXY);
  histos->UserHistogram("Track","dZ","dZ;dZ [cm];#tracks",200,-10.,10.,AliDielectronVarManager::kImpactParZ);

  
  //add histograms to Pair classes
  histos->UserHistogram("Pair","InvMass","Inv.Mass;Inv. Mass [GeV];#pairs",
                        125,1.,125.*.04,AliDielectronVarManager::kM);
  histos->UserHistogram("Pair","Rapidity","Rapidity;Rapidity;#pairs",
                        100,-1.,1.,AliDielectronVarManager::kY);
  histos->UserHistogram("Pair","OpeningAngle","Opening angle;angle",
                        100,1.,3.15,AliDielectronVarManager::kOpeningAngle);
  //
  
   histos->UserHistogram("Pair","Eta","eta;Eta;#pairs", 200,-2.,2.,AliDielectronVarManager::kEta);
 // histos->UserHistogram("Pair","OpeningAngle","Opening angle;angle",
 //                       100,0.,3.15,AliDielectronVarManager::kOpeningAngle);
  histos->UserHistogram("Pair","InvMass_OpeningAngle","Opening angle:Inv.Mass;Inv. Mass [GeV];angle",
                        125,1.,125.*.04,100,0.,3.15,AliDielectronVarManager::kM,AliDielectronVarManager::kOpeningAngle);

  
  

  die->SetHistogramManager(histos);

}


#if 0
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
  //  printf("iPeriod: %d \t %.0f-%.0f \n",iPeriod,first,last);
  return (AliDielectronHelper::MakeLinBinning(last-first, first, last));
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



#endif
