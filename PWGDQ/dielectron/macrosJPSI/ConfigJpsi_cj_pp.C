void SetupEventCutsDieleFilter(AliDielectron *diele, Int_t cutDefinition, Bool_t isAOD, Int_t MultSel);
void SetupTrackCutsDieleData(AliDielectron *diele, Int_t cutDefinition, Bool_t isAOD, Bool_t isMC);
void SetupPairCutsDieleData(AliDielectron *diele, Int_t cutDefinition, Bool_t isAOD, Bool_t isMC);
void InitHistogramsDieleData(AliDielectron *diele, Int_t cutDefinition, Bool_t isAOD);
void InitCFDieleData(AliDielectron *diele, Int_t cutDefinition, Bool_t isAOD, Bool_t isMC);
void AddMCSignals(AliDielectron *diele);


AliESDtrackCuts *SetupESDtrackCutsDieleData(Int_t cutDefinition);

TString namesDieleData=("TPC;EMCal;EMCal2");//EMCal2 with loose cuts only for CFContainers!

TObjArray *arrNamesDieleData=namesDieleData.Tokenize("; ");

const Int_t nDie=arrNamesDieleData->GetEntries();

AliDielectron* ConfigJpsi_cj_pp(Int_t cutDefinition, Bool_t isAOD=kFALSE, Int_t trigger_index=0, Bool_t isMC, Int_t MultSel=0)
{
  //
  // Setup the instance of AliDielectron
  //
  
  // create the actual framework object
  TString name=Form("%02d",cutDefinition);
  if (cutDefinition<arrNamesDieleData->GetEntriesFast()){
    name=arrNamesDieleData->At(cutDefinition)->GetName();
  }
  AliDielectron *diele = new AliDielectron(Form("%s",name.Data()),
                             Form("Track cuts: %s",name.Data()));
	
	
	//MC signals
	//from Marcel December 13th
	if(isMC) {
		AddMCSignals(diele);
		printf(" Add %d MC signals \n",diele->GetMCSignals()->GetEntriesFast());
	}
  // cut setup
  SetupEventCutsDieleFilter(diele, cutDefinition, isAOD, MultSel);
  SetupTrackCutsDieleData(diele, cutDefinition, isAOD, isMC);
  SetupPairCutsDieleData(diele, cutDefinition, isAOD, trigger_index, isMC);

  
  //

  // histogram setup
  // only if an AliDielectronHistos object is attached to the
  // dielelectron framework histograms will be filled
  //
  InitHistogramsDieleData(diele, cutDefinition, isAOD);
		// InitCFDieleData(diele, cutDefinition, isAOD, isMC);
	

//   AliDielectronTrackRotator *rot=new AliDielectronTrackRotator;
//   rot->SetIterations(20);
//   rot->SetConeAnglePhi(TMath::Pi()/2);
//   rot->SetStartAnglePhi(TMath::Pi()/2);
//   diele->SetTrackRotator(rot);

	
// mixing
/*
  AliDielectronMixingHandler *mix=new AliDielectronMixingHandler;
  mix->AddVariable(AliDielectronVarManager::kZvPrim, 100,-10.,10.);
  mix->SetMixType(AliDielectronMixingHandler::kAll);
  mix->SetDepth(20);
  diele->SetMixingHandler(mix);
 */
	
  
  return diele;
}

//______________________________________________________________________________________
void SetupEventCutsDieleFilter(AliDielectron *diele, Int_t cutDefinition, Bool_t isAOD, Int_t MultSel)
{

//   // Setup the event cuts
	/*
	AliDielectronEventCuts *eventCuts=new AliDielectronEventCuts("eventCuts","Vertex Track && |vtxZ|<10 && ncontrib>0");
	if(isAOD) eventCuts->SetVertexType(AliDielectronEventCuts::kVtxAny);
	eventCuts->SetRequireVertex();
	eventCuts->SetMinVtxContributors(1);
	eventCuts->SetVertexZ(-10.,10.);
		//eventCuts->SetCentralityRange(0.0,80.0);
		//task->SetEventFilter(eventCuts);
	 */
	

	AliDielectronVarCuts *Mult = new AliDielectronVarCuts("Mult","kNaccTrcklts10Corr cut");
	if(MultSel==1){
		Mult->AddCut(AliDielectronVarManager::AliDielectronVarManager::kNaccTrcklts10Corr,0.,29.99); 
	}
	if(MultSel==2){
		Mult->AddCut(AliDielectronVarManager::AliDielectronVarManager::kNaccTrcklts10Corr,30.,59.99); 
	}
	if(MultSel==3){
		Mult->AddCut(AliDielectronVarManager::AliDielectronVarManager::kNaccTrcklts10Corr,60.,99.99); 
	}
	if(MultSel==4){
		Mult->AddCut(AliDielectronVarManager::AliDielectronVarManager::kNaccTrcklts10Corr,100.,1000.); 
	}
	
	diele->GetEventFilter().AddCuts(Mult);
  
	
}





//______________________________________________________________________________________
void SetupTrackCutsDieleData(AliDielectron *diele, Int_t cutDefinition, Bool_t isAOD, Bool_t isMC)
{
  //
  // Setup the track cuts
  //
  
  //ESD quality cuts DielectronTrackCuts
  if (!isAOD) {
    diele->GetTrackFilter().AddCuts(SetupESDtrackCutsDieleData(cutDefinition));
  } 
  else {
    AliDielectronTrackCuts *trackCuts=new AliDielectronTrackCuts("trackCuts","trackCuts");
    trackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
    trackCuts->SetRequireTPCRefit(kTRUE);
    trackCuts->SetRequireITSRefit(kTRUE);
    diele->GetTrackFilter().AddCuts(trackCuts);
  }

  //Pt cut ----------------------------------------------------------
  AliDielectronVarCuts *pt = new AliDielectronVarCuts("ptCut","pt cut");
  pt->AddCut(AliDielectronVarManager::kPt,1.0,1e30);
  pt->AddCut(AliDielectronVarManager::kKinkIndex0,0.);
  
  //AOD additions since there are no AliESDtrackCuts -----------------
  //
  if (isAOD){
	pt->AddCut(AliDielectronVarManager::kNclsTPC,85.,160.);
    pt->AddCut(AliDielectronVarManager::kEta,-0.9.,0.9);
		  //pt->AddCut(AliDielectronVarManager::kPhi,0.,4.);
    pt->AddCut(AliDielectronVarManager::kImpactParXY,-1.,1.);
    pt->AddCut(AliDielectronVarManager::kImpactParZ,-3.,3.);	  
    pt->AddCut(AliDielectronVarManager::kITSLayerFirstCls,0.,4.);
		  //new
	pt->AddCut(AliDielectronVarManager::kTPCchi2Cl,0.,4.);
	  
	if(cutDefinition==2){
		pt->AddCut(AliDielectronVarManager::kNclsTPC,70.,160.);
		pt->AddCut(AliDielectronVarManager::kEta,-0.9.,0.9);
			//pt->AddCut(AliDielectronVarManager::kPhi,0.,4.);
		pt->AddCut(AliDielectronVarManager::kImpactParXY,-2.,2.);
		pt->AddCut(AliDielectronVarManager::kImpactParZ,-4.,4.);	  
		pt->AddCut(AliDielectronVarManager::kITSLayerFirstCls,0.,6.);
			//new
		pt->AddCut(AliDielectronVarManager::kTPCchi2Cl,0.,10.);
	  
	  
	}
	  	  
	pt->AddCut(AliDielectronVarManager::kTPCnSigmaEle,-2.25,3.0);		  
	if(cutDefinition==2)pt->AddCut(AliDielectronVarManager::kTPCnSigmaEle,-3.0,4.0);
	  	
		 	  
  }
	
  diele->GetTrackFilter().AddCuts(pt);

}

//______________________________________________________________________________________
void SetupPairCutsDieleData(AliDielectron *diele, Int_t cutDefinition, Bool_t isAOD, Int_t trigger_index, Bool_t isMC)
{
  // Setup the pair cuts
  //Invariant mass and rapidity selection
  Double_t gCut  = 0.050;             // default
  

	//AliDielectronVarCuts *gammaCut=new AliDielectronVarCuts("gammaCut","gammaCut");
	//if(cutDefinition==1)	gammaCut->AddCut(AliDielectronVarManager::kM,0.,gCut);
	//if(cutDefinition>0)	diele->GetPairPreFilter().AddCuts(gammaCut);
	//if(cutDefinition>0)	diele->SetPreFilterUnlikeOnly();

  //Invariant mass and rapidity selection
  AliDielectronVarCuts *pairCut=new AliDielectronVarCuts("0<M<5+|Y|<.9","0<M<5 + |Y|<.9");
  pairCut->AddCut(AliDielectronVarManager::kM,0.,15.);
  pairCut->AddCut(AliDielectronVarManager::kY,-0.9,0.9);
  pairCut->AddCut(AliDielectronVarManager::kPt,1,50.);
	
  if(cutDefinition==1){
		  //EMC7 trigger
	    if(trigger_index == 1)pairCut->AddCut(AliDielectronVarManager::kPt,3,50.);
	  
		if(trigger_index == 3 )pairCut->AddCut(AliDielectronVarManager::kPt,7,50.);
		if(trigger_index == 4 || trigger_index == 40)pairCut->AddCut(AliDielectronVarManager::kPt,5,50.);
	    if(trigger_index == 6 || trigger_index == 60 )pairCut->AddCut(AliDielectronVarManager::kPt,11,50.);
	  
		
  }
  
  diele->GetPairFilter().AddCuts(pairCut);
		
  if(cutDefinition==1){
		AliDielectronVarCuts *mycut = new AliDielectronVarCuts("CutEMCAL","cut for EMCal");
	    mycut->AddCut(AliDielectronVarManager::kEMCALEoverP,0.8,1.3);
		  //trigger_index == 2 for both thresholds together, without thresholds separation (just cross check)
	    if(trigger_index == 3 || trigger_index == 2) mycut->AddCut(AliDielectronVarManager::kEMCALE,7,50.);
	    else if(trigger_index == 4 || trigger_index == 40) mycut->AddCut(AliDielectronVarManager::kEMCALE,5,50.);	
		  //for 16k period, threshold at 10 GeV:
		  //Dcal have both thresholds at 10 GeV!!
		  //trigger 20 do not separate thresholds for dcal
	    else if(trigger_index == 6 || trigger_index == 60  ) mycut->AddCut(AliDielectronVarManager::kEMCALE,11,50.);
		  //EMC7:
	    else if(trigger_index == 1) mycut->AddCut(AliDielectronVarManager::kEMCALE,3,50.);
        else mycut->AddCut(AliDielectronVarManager::kEMCALE,1,50.);
		 
		//dcal cut if using EG triggers
		//if using EMCal trigger, it excludes if the tracks matches DCAL:
		//(we don't want energy cut on DCAL if the trigger is on EMCal)
	    if(trigger_index == 2 ||trigger_index == 3 || trigger_index == 4 || trigger_index == 6) mycut->AddCut(AliDielectronVarManager::kPhi,4.377,5.7071, kTRUE);//kTRUE means to exclude!
	  
	  
		//emcal cut if using DG triggers
	    if(trigger_index == 40 || trigger_index == 60) mycut->AddCut(AliDielectronVarManager::kPhi,1.396,3.2637, kTRUE);

		AliDielectronPairLegCuts *varpair=new AliDielectronPairLegCuts();	
		varpair->GetLeg1Filter().AddCuts(mycut);
		varpair->GetLeg2Filter().AddCuts(mycut);
		varpair->SetCutType(AliDielectronPairLegCuts::kAnyLeg);	  
		diele->GetPairFilter().AddCuts(varpair);
  }
	
	if(cutDefinition==2){
		AliDielectronVarCuts *mycut = new AliDielectronVarCuts("CutEMCAL","cut for EMCal");
	    mycut->AddCut(AliDielectronVarManager::kEMCALEoverP,0.7,1.3);
        mycut->AddCut(AliDielectronVarManager::kEMCALE,4,50.);
		
			//dcal cut if using EG triggers
			//if using EMCal trigger, it excludes if the tracks matches DCAL:
			//(we don't want energy cut on DCAL if the trigger is on EMCal)
	    if(trigger_index == 2 ||trigger_index == 3 || trigger_index == 4 || trigger_index == 6) mycut->AddCut(AliDielectronVarManager::kPhi,4.377,5.7071, kTRUE);//kTRUE means to exclude!
		
		
			//emcal cut if using DG triggers
	    if(trigger_index == 40 || trigger_index == 60) mycut->AddCut(AliDielectronVarManager::kPhi,1.396,3.2637, kTRUE);
		
		AliDielectronPairLegCuts *varpair=new AliDielectronPairLegCuts();	
		varpair->GetLeg1Filter().AddCuts(mycut);
		varpair->GetLeg2Filter().AddCuts(mycut);
		varpair->SetCutType(AliDielectronPairLegCuts::kAnyLeg);	  
		diele->GetPairFilter().AddCuts(varpair);
	}

	
	

}

//______________________________________________________________________________________
AliESDtrackCuts *SetupESDtrackCutsDieleData(Int_t cutDefinition)
{
  //
  // Setup default AliESDtrackCuts
  //
  AliESDtrackCuts *esdTrackCuts = new AliESDtrackCuts;

  // basic track quality cuts  (basicQ)
  esdTrackCuts->SetMaxDCAToVertexZ(3.0);
  esdTrackCuts->SetMaxDCAToVertexXY(1.0);

  esdTrackCuts->SetEtaRange( -0.9 , 0.9 );

  esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
  esdTrackCuts->SetRequireITSRefit(kTRUE);
  esdTrackCuts->SetRequireTPCRefit(kTRUE);

  esdTrackCuts->SetPtRange(1,1e30);

  esdTrackCuts->SetMinNClustersTPC(85);
  esdTrackCuts->SetMaxChi2PerClusterTPC(4);

  // default SPD any
  esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);

  return esdTrackCuts;
}


//______________________________________________________________________________________
void InitHistogramsDieleData(AliDielectron *diele, Int_t cutDefinition, Bool_t isAOD)
{
  //
  // Initialise the histograms
  //
  
  //Setup histogram Manager
  AliDielectronHistos *histos=new AliDielectronHistos(diele->GetName(),diele->GetTitle());
  
  //Initialise histogram classes
  histos->SetReservedWords("Track;Pair");
  
  //Track classes
  //EMCal tracks are same as TPC tracks... only changes the pair cuts
	//to fill also track info from 2nd event loop until 2
	if(cutDefinition!=2){
			for (Int_t i=0; i<2; ++i){
					histos->AddClass(Form("Track_%s",AliDielectron::TrackClassName(i)));
			}
		
	}
		
  
  //Pair classes
  // to fill also mixed event histograms loop until 10

	if(cutDefinition!=2){	
	
		for (Int_t i=0; i<3; ++i){
				histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(i)));

		}
  
			//legs from pair
		for (Int_t i=0; i<3; ++i){
				histos->AddClass(Form("Track_Legs_%s",AliDielectron::PairClassName(i)));
		}
		
	}
	
	
  //track rotation
  //histos->AddClass(Form("Pair_%s",PairClassName(AliDielectron::kEv1PMRot)));
  //histos->AddClass(Form("Track_Legs_%s",PairClassName(AliDielectron::kEv1PMRot)));
  

    //add histograms to event class
if(cutDefinition==0){
    histos->AddClass("Event");
    histos->UserHistogram("Event","VtxZ","Vertex Z;Z[cm]",500,-40.,40.,AliDielectronVarManager::kZvPrim);
    histos->UserHistogram("Event","VtxYxVtxZ","Vertexyz;Z[cm];Y[cm]",500,-40.,40.,400,-0.5,0.5,AliDielectronVarManager::kZvPrim,AliDielectronVarManager::kYvPrim);
    histos->UserHistogram("Event","VtxXxVtxZ","Vertexxz;Z[cm];X[cm]",500,-40.,40.,400,-0.5,0.5,AliDielectronVarManager::kZvPrim,AliDielectronVarManager::kXvPrim);
    histos->UserHistogram("Event","VtxYxVtxX","Vertexxz;Z[cm];X[cm]",400,-0.5,0.5,400,-0.5,0.5,AliDielectronVarManager::kXvPrim,AliDielectronVarManager::kYvPrim);

		//kNaccTrcklts10Corr
		//histos->UserHistogram("Event","SPDTracklets","SPDTracklets;SPDTracklets;Entries",300,0.,300.,AliDielectronVarManager::kCentralitySPDTracklets);	

	
	histos->UserHistogram("Event","kNaccTrcklts10Corr","kNaccTrcklts10Corr;kNaccTrcklts10Corr;Entries",300,0.,300.,AliDielectronVarManager::kNaccTrcklts10Corr);	
}
 
  //add histograms to Track classes
	
 

	
  histos->UserHistogram("Track","Pt","Pt;Pt [GeV];#tracks",100,0,20.,AliDielectronVarManager::kPt,kTRUE);
	/*
  histos->UserHistogram("Track","TPCnCls","Number of Clusters TPC;TPC number clusteres;#tracks",160,-0.5,159.5,AliDielectronVarManager::kNclsTPC,kTRUE);
  histos->UserHistogram("Track","TPCchi2Cl","Chi-2/Clusters TPC;Chi2/ncls number clusteres;#tracks",100,0,10,AliDielectronVarManager::kTPCchi2Cl,kTRUE);
  histos->UserHistogram("Track","TPCnFCls","Number of findable Clusters TPC;TPC number clusteres;#tracks",160,-0.5,159.5,AliDielectronVarManager::kNFclsTPC,kTRUE);
  histos->UserHistogram("Track","TPCnFClsfCross","fraction crossed rows/findable;TPC number clusteres;#tracks",200,0.,2.,AliDielectronVarManager::kNFclsTPCfCross,kTRUE);
  histos->UserHistogram("Track","TPCnFClsr","Number of findable Clusters(crossed rows) TPC;TPC number clusteres;#tracks",160,-0.5,159.5,AliDielectronVarManager::kNFclsTPCr,kTRUE);
  histos->UserHistogram("Track","TPCnFClsrFrac","Number of found/findable Clusters TPC;TPC number clusteres;#tracks",200,0.,2.,AliDielectronVarManager::kNFclsTPCrFrac,kTRUE);
  histos->UserHistogram("Track","TPCnFClsTPCfCross","Fraction of rows/findable Clusters TPC;TPC number clusteres;#tracks",200,0.,2.,AliDielectronVarManager::kNFclsTPCfCross,kTRUE);  
  histos->UserHistogram("Track","TPCsignalN","Number of points for TPC Signal;TPC Npoints dEdx;#tracks",160,-0.5,159.5,AliDielectronVarManager::kTPCsignalN,kTRUE);    
 */
	
  //histos->UserHistogram("Track","SPDTracklets","SPDTracklets;SPDTracklets;Entries",300,0.,300.,AliDielectronVarManager::kCentralitySPDTracklets);
  // histos->UserHistogram("Track","kNaccTrcklts10Corr","kNaccTrcklts10Corr;kNaccTrcklts10Corr;Entries",300,0.,300.,AliDielectronVarManager::kNaccTrcklts10Corr);
	

  histos->UserHistogram("Track","dXY","dXY;dXY [cm];#tracks",200,-1.,1.,AliDielectronVarManager::kImpactParXY,kTRUE);
  
  histos->UserHistogram("Track","dZ","dZ;dZ [cm];#tracks",200,-3.,3.,AliDielectronVarManager::kImpactParZ,kTRUE);
  
  histos->UserHistogram("Track","Eta_Phi","Eta Phi Map; Eta; Phi;#tracks",
                        100,-1,1,144,0,6.285,AliDielectronVarManager::kEta,AliDielectronVarManager::kPhi,kTRUE);

  histos->UserHistogram("Track","dEdx_P","dEdx;P [GeV];TPC signal (arb units);#tracks",
                        200,0.2,20.,800,20.,200.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCsignal,kTRUE);

  histos->UserHistogram("Track","dEdx_Pt","dEdx;Pt [GeV];TPC signal (arb units);#tracks",
                        200,0.2,20.,800,20.,200.,AliDielectronVarManager::kPt,AliDielectronVarManager::kTPCsignal,kTRUE);


  histos->UserHistogram("Track","TPCnSigmaEle_P","TPCnSigmaEle;P [GeV];TPCnSigmaEle;#tracks",
                        200,0.2,20.,800,-12.,12.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,kTRUE);
  
  histos->UserHistogram("Track","TPCnSigmaEle_Pt","TPCnSigmaEle;Pt [GeV];TPCnSigmaEle;#tracks",
                        200,0.2,20.,800,-12.,12.,AliDielectronVarManager::kPt,AliDielectronVarManager::kTPCnSigmaEle,kTRUE);

  histos->UserHistogram("Track","TPCnSigmaEle_Phi","TPCnSigmaEle;#phi [rad];TPCnSigmaEle;#tracks",
                        200,0.,2*TMath::Pi(),800,-12.,12.,AliDielectronVarManager::kPhi,AliDielectronVarManager::kTPCnSigmaEle,kTRUE);
  
  histos->UserHistogram("Track","TPCnSigmaEle_Eta","TPCnSigmaEle;#eta;TPCnSigmaEle;#tracks",
                        200,-1.,1.,800,-12.,12.,AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEle,kTRUE);
    

  histos->UserHistogram("Track","dEdx_Phi","dEdx;#phi [rad];TPC signal (arb units);#tracks",
                        200,0.,2*TMath::Pi(),800,20.,200.,AliDielectronVarManager::kPhi,AliDielectronVarManager::kTPCsignal,kTRUE);

  histos->UserHistogram("Track","dEdx_Eta","dEdx;#eta;TPC signal (arb units);#tracks",
                        200,-1.,1.,800,20.,200.,AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCsignal,kTRUE);


  histos->UserHistogram("Track","dEdx_nSigmaEMCal","dEdx;NsigmaEmcal;TPC signal (arb units);NSigmaEMCAL",
                        200,-5.,5.,800,20.,200.,AliDielectronVarManager::kEMCALnSigmaEle,AliDielectronVarManager::kTPCsignal,kTRUE);

  histos->UserHistogram("Track","dEdx_TPCnSigmaEle","dEdx;TPC signal (arbunits);TPC number of sigmas Electrons;TPC signal (a.u.);#tracks",
                        100,-10.,10.,800,20.,200.,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kTPCsignal,kTRUE);

  histos->UserHistogram("Track","dEdx_EoverP","dEdx;EoverP;TPC signal (arbunits);E/P",
						100,0.,5.,800,20.,200.,AliDielectronVarManager::kEMCALEoverP,AliDielectronVarManager::kTPCsignal,kTRUE);
  
  histos->UserHistogram("Track","nSigmaEMCal_EoverP","NsigmaEmcal;EoverP;NSigmaEMCAL;E/P"
						,100,0.,5.,200,-5.,5.,AliDielectronVarManager::kEMCALEoverP,AliDielectronVarManager::kEMCALnSigmaEle,kTRUE);
 
  histos->UserHistogram("Track","EMCal_E","EmcalE;Cluster Energy [GeV];#Clusters",
						200,0.,40.,AliDielectronVarManager::kEMCALE,kTRUE);

  histos->UserHistogram("Track","ITS_FirstCls","ITS First Layer;ITS First Layer;#Entries",
						6,0.,6.,AliDielectronVarManager::kITSLayerFirstCls,kTRUE);
 
	
  //add histograms to Pair classes
  //new Cris
	
  //Ecluster versus Phi to separate EMCal and DCal
  histos->UserHistogram("Track","EMCal_E_Phi","Energy vs. Phi; EMCal_E;Phi",
						200,0.,40.,20,0.,2*TMath::Pi(),AliDielectronVarManager::kEMCALE,AliDielectronVarManager::kPhi,kTRUE);

	
  histos->UserHistogram("Track","TPCnSigmaEle_EoverP","TPCnSigmaEle;EoverP;TPC signal (arbunits);E/P",
						100,0.,5.,100,-10.,10.,AliDielectronVarManager::kEMCALEoverP,AliDielectronVarManager::kTPCnSigmaEle,kTRUE);
  histos->UserHistogram("Track","EoverP_pt","Pt;EoverP;p_{T} (GeV/c);E/p",
						3000,0.,30.,200,0.,2.,AliDielectronVarManager::kPt,AliDielectronVarManager::kEMCALEoverP,kTRUE);

  histos->UserHistogram("Pair","OpeningAngle2D","Opening angle vs p_{T} ;p_{T} (GeV/c); angle",
						300,0.,30., 2000,-10.,10.,AliDielectronVarManager::kPt, AliDielectronVarManager::kOpeningAngle);
	
  histos->UserHistogram("Pair","InvMass","Inv.Mass;Inv. Mass [GeV];#pairs",
                        375,0.0,15.0,AliDielectronVarManager::kM);
	
  histos->UserHistogram("Pair","InvMass2D","Inv.Mass;Pt [GeV]; Inv. Mass [GeV]",
                        20,0.,20.,375,0,15.0,AliDielectronVarManager::kPt,AliDielectronVarManager::kM);
	
  
  histos->UserHistogram("Pair","Rapidity","Rapidity;Rapidity;#pairs",
                        50,-1.,1.,AliDielectronVarManager::kY);
  histos->UserHistogram("Pair","OpeningAngle","Opening angle;angle",
                        50,0.,3.15,AliDielectronVarManager::kOpeningAngle);
	
  histos->UserHistogram("Pair","PseudoProperTime","Pseudoproper decay length; pseudoproper-decay-length[#mum];Entries/40#mum",
                          150,-0.3.,0.3,AliDielectronVarManager::kPseudoProperTime);
	
	
 //histos->UserHistogram("Pair","SPDTracklets","SPDTracklets;SPDTracklets;Entries",300,0.,300.,AliDielectronVarManager::kCentralitySPDTracklets);
 //histos->UserHistogram("Pair","kNaccTrcklts10Corr","kNaccTrcklts10Corr;kNaccTrcklts10Corr;Entries",300,0.,300.,AliDielectronVarManager::kNaccTrcklts10Corr);

	
	
	
  
/*
 histos->UserHistogram("Pair","Chi2/NDF","#Chi^{2}/NDF;#Chi^{2}/NDF",
                        100, 0., 20., AliDielectronVarManager::kChi2NDF);
	
*/
	
  
  diele->SetHistogramManager(histos);
}


void InitCFDieleData(AliDielectron *diele, Int_t cutDefinition, Bool_t isAOD, Bool_t isMC)
{
  //
  // Setupd the CF Manager if needed
  //
  if(cutDefinition==1||cutDefinition==2){
	
	
  AliDielectronCF *cf=new AliDielectronCF(diele->GetName(),diele->GetTitle());
  //for centrality selection
  // cf->AddVariable(AliDielectronVarManager::kCentralitySPDTracklets, "0,20,40,60,100, 200, 300");
  //cf->AddVariable(AliDielectronVarManager::kNaccTrcklts10Corr, "0,30,60,100,1000");
  //pair variables
	
  cf->AddVariable(AliDielectronVarManager::kPt,"5.0,7.0,9.0,11.0,15.0,20.0");
  cf->AddVariable(AliDielectronVarManager::kM,375,0.,375*.04); //40Mev Steps
  
  cf->AddVariable(AliDielectronVarManager::kPairType,3,0,3);
	
  cf->AddVariable(AliDielectronVarManager::kOpeningAngle,16,0,3.2);
	
  cf->AddVariable(AliDielectronVarManager::kEta,"-1.0,-0.9,-0.8,0,0.8,0.9,1.0");
		// cf->AddVariable(AliDielectronVarManager::kY,40,-1.,1.);
  cf->AddVariable(AliDielectronVarManager::kPhi,"0,1.0,2.0,2.5,3.0,3.5,4.0,4.5,5.0,6.0,6.5");
	
		//cf->AddVariable(AliDielectronVarManager::kPseudoProperTime,300,-0.3,0.3);
		//cf->AddVariable(AliDielectronVarManager::kPseudoProperTimeErr,200,0.,0.1);
		//cf->AddVariable(AliDielectronVarManager::kPseudoProperTimeResolution,400,-0.1,0.1);
		//cf->AddVariable(AliDielectronVarManager::kPseudoProperTimePull,400,-0.1,0.1); 
		//cf->AddVariable(AliDielectronVarManager::kChi2NDF,100, 0., 20.);
	
   
  //leg variables
  cf->AddVariable(AliDielectronVarManager::kPt,"1.0,3.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,15.0,20.0",kTRUE);
  
  cf->AddVariable(AliDielectronVarManager::kNclsTPC,"70, 80, 85, 90, 100,160",kTRUE);
  cf->AddVariable(AliDielectronVarManager::kTPCchi2Cl,5, 0., 10.,kTRUE);
		  //cf->AddVariable(AliDielectronVarManager::kTPCsignalN,160,-0.5,159.5,kTRUE);
	
  cf->AddVariable(AliDielectronVarManager::kEta,"-1.0,-0.9,-0.8,0,0.8,0.9,1.0",kTRUE);
  cf->AddVariable(AliDielectronVarManager::kPhi,"0,1.0,2.0,2.5,3.0,3.5,4.0,4.5,5.0,6.0,6.5",kTRUE);
   
	
  cf->AddVariable(AliDielectronVarManager::kEMCALE,"4.5, 5.0, 5.5, 6.5, 7.0, 7.5,  10.0, 10.5, 11.0, 11.5,  15.0, 20.0, 30.0",kTRUE); 
		//cf->AddVariable(AliDielectronVarManager::kEMCALnSigmaEle,"-4.5,-4.,-3.75,-3.5,-3.25,-3.0,-2.75,-2.5,-2.0,-1.0,1.0,2.0,2.5,3.0,3.25,3.5,3.75,4.0,5.0",kTRUE);
		// cf->AddVariable(AliDielectronVarManager::kEMCALNCells,25,0,25,kTRUE);
  cf->AddVariable(AliDielectronVarManager::kEMCALEoverP,"0,0.4,0.6,0.7,0.8,0.9,1.0,1.1, 1.2,1.3,1.4,2",kTRUE);
		// cf->AddVariable(AliDielectronVarManager::kTPCsignal,"40.,50.,55.,60.,65.,68.,70.,72.,75.,80.,90.,100.,110.,200.",kTRUE);
	
  cf->AddVariable(AliDielectronVarManager::kTPCnSigmaEle," -3.0,-2.5, -2.0, -1.0, 0, 1.0, 2.0, 2.5, 3.0, 4.0",kTRUE);
	
  if(!isMC){
	 cf->AddVariable(AliDielectronVarManager::kTPCnSigmaPio,"-10, -3.0, -2.0,-1.5, -1.0, 0, 1.0,1.5, 2.0, 3.0, 4.0, 10",kTRUE);
	 cf->AddVariable(AliDielectronVarManager::kTPCnSigmaPro,"-10, -3.0, -2.0,-1.5, -1.0, 0, 1.0,1.5, 2.0, 3.0, 4.0, 10",kTRUE);
	 cf->AddVariable(AliDielectronVarManager::kTPCnSigmaKao,"-10, -3.0, -2.0,-1.5, -1.0, 0, 1.0,1.5, 2.0, 3.0, 4.0, 10",kTRUE);
  }
	
//cf->AddVariable(AliDielectronVarManager::kTOFnSigmaEle,",-3.5,-3.0,-2.75,-2.5,-2.25,-2.0,-1.75,-1.5,-1.0,-0.5,0.0,1.0,2.0,2.25,2.5,2.75,3.0,3.25,3.5,4.0",kTRUE);
//cf->AddVariable(AliDielectronVarManager::kTOFnSigmaPio,"1.,2.,2.5,3.0,3.5,4.0,4.5,100",kTRUE);
	
  cf->AddVariable(AliDielectronVarManager::kITSLayerFirstCls,6,0.,6.,kTRUE);
		  //cf->AddVariable(AliDielectronVarManager::kZvPrim,20,-20.,20.);
  cf->AddVariable(AliDielectronVarManager::kImpactParXY,"-2.0, -1.0,-0.5,0.5, 1.0, 2.0",kTRUE);
  cf->AddVariable(AliDielectronVarManager::kImpactParZ,"-4.0,-3.0, -2.0, 2.0,3.0, 4.0",kTRUE);
    
	
	// if (!isAOD){
	// Bool_t hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);
	/*
    if (isMC){
			// printf("Is MC in the containers!!!\n");
      cf->AddVariable(AliDielectronVarManager::kPdgCode,10000,-5000.5,4999.5,kTRUE);
      cf->AddVariable(AliDielectronVarManager::kPdgCodeMother,10000,-5000.5,4999.5,kTRUE);
      cf->AddVariable(AliDielectronVarManager::kPdgCodeGrandMother,10000,-5000.5,4999.5,kTRUE);
		// }
	}
	 */
	
	 
	
    //only in this case write MC truth info
	if((cutDefinition==1) && isMC){
		cf->SetStepForMCtruth();
	}

   diele->SetCFManagerPair(cf);
		
  }
  
}
void AddMCSignals(AliDielectron *diele){
	
	AliDielectronSignalMC* inclusiveJpsi = new AliDielectronSignalMC("inclusiveJpsi","Inclusive J/psi");
	inclusiveJpsi->SetLegPDGs(11,-11);
	inclusiveJpsi->SetMotherPDGs(443,443);
	inclusiveJpsi->SetMothersRelation(AliDielectronSignalMC::kSame);
	inclusiveJpsi->SetFillPureMCStep(kTRUE);
	inclusiveJpsi->SetCheckBothChargesLegs(kTRUE,kTRUE);
	inclusiveJpsi->SetCheckBothChargesMothers(kTRUE,kTRUE);
	diele->AddSignalMC(inclusiveJpsi);
	
	
	
	/*
		// all prompt J/psi
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
	diele->AddSignalMC(promptJpsi);
*/
}

