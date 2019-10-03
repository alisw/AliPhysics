void InitHistograms(AliDielectron *die, Int_t cutDefinition);
void InitCF(AliDielectron* die, Int_t cutDefinition);

void SetupTrackCuts(AliDielectron *die, Int_t cutDefinition);
void SetupPairCuts(AliDielectron *die, Int_t cutDefinition);
void SetupV0Cuts(AliDielectron *die, Int_t cutDefinition);

void SetEtaCorrection();



// define different cuts


UInt_t maskSPDany = 1<<0;
UInt_t maskTrackPt23 = (1<<1) ;
UInt_t maskTrackPt35 = (1<<2);
UInt_t maskTrackPt = maskTrackPt23  | maskTrackPt35;
UInt_t maskTriggeredQU = (1<<3) ;
UInt_t maskTriggeredSE = (1<<4);
UInt_t maskTriggered = maskTriggeredQU | maskTriggeredSE;


Bool_t hasMC = (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);



AliDielectron* ConfigJpsi(AliAnalysisTaskMultiDielectron *task)
{
  //
  // Setup the instance of AliDielectron
  //
  vector<Int_t> cutsToUse;
  cutsToUse.push_back(0);
  cutsToUse.push_back(maskSPDany);
   cutsToUse.push_back(maskTrackPt23);
   cutsToUse.push_back(maskTrackPt23 | maskSPDany);
  cutsToUse.push_back(maskTriggeredQU);
  cutsToUse.push_back(maskTriggeredQU | maskSPDany);
  cutsToUse.push_back(maskTriggeredQU | maskTrackPt23);
  cutsToUse.push_back(maskTriggeredQU | maskTrackPt23 | maskSPDany);
  TString names = "default+SPDfirst;default+SPDany;Pt>2.3+SPDfirst;Pt>2.3+SPDany;QUtrig+SPDfirst;QUtrig+SPDany;QUtrig+Pt>2.3+SPDfirst;QUtrig+Pt>2.3+SPDany";
  TObjArray *arrNames=names.Tokenize(";");
  
/**
*
* event cuts    
*
**/
    AliDielectronEventCuts *eventCuts=new AliDielectronEventCuts("eventCuts","Vertex Track && |vtxZ|<10 && ncontrib>0");
    eventCuts->SetRequireVertex();
    eventCuts->SetMinVtxContributors(1);
    eventCuts->SetVertexZ(-10.,10.);
//   eventCuts->SetCentralityRange(0.,80.);
    task->SetEventFilter(eventCuts);
  
  
  for(Int_t cutNumber=0; cutNumber < cutsToUse.size() ; cutNumber++){
	Int_t cutDefinition = cutsToUse[cutNumber];
	// create the actual framework object
	TString name=Form("%02d",cutDefinition);
	if (cutNumber < arrNames->GetEntriesFast()){
		name=arrNames->At(cutNumber)->GetName();
	}
	AliDielectron *die =
		new AliDielectron(Form("%s",name.Data()),
						Form("Track cuts: %s",name.Data()));

	if (hasMC){
			SetupMCsignals(die);
			InitCF(die,cutDefinition);
		}
		
	else{
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
	}
			
			// cut setup
	SetupTrackCuts(die,cutDefinition);  
	//
	
	
		SetupPairCuts(die,cutDefinition);

		
		
		
	// ConfDef::SetupV0Cuts(die, cutDefinition);
		
	//
	InitHistograms(die,cutDefinition);
	//


	
	//
	// ConfDef::SetEtaCorrection();
	task->AddDielectron(die);
  }
  
}

	
/**
*	
*	Setup the track cuts according to given cut definition
*	
**/
void SetupTrackCuts(AliDielectron *die, Int_t cutDefinition)
{
	
	
  AliDielectronCutGroup* cuts = new AliDielectronCutGroup("cuts","cuts",AliDielectronCutGroup::kCompAND);
  die->GetTrackFilter().AddCuts(cuts);

  //default quality cuts
  AliDielectronTrackCuts *cut1=new AliDielectronTrackCuts("cut1","cut1");
  cut1->SetRequireITSRefit(kTRUE);
  cut1->SetRequireTPCRefit(kTRUE);
  if( cutDefinition & maskSPDany ){
		cut1->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
	}
	else{
		cut1->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
	}
	cuts->AddCut(cut1);
	
	

	AliDielectronVarCuts *trackQuality = new AliDielectronVarCuts("trackQuality1","trackQuality1");
	trackQuality->AddCut(AliDielectronVarManager::kImpactParZ, -3. ,  3. );
	trackQuality->AddCut(AliDielectronVarManager::kImpactParXY,-1. ,   1. );
	trackQuality->AddCut(AliDielectronVarManager::kEta,        - .9,   .9);
	trackQuality->AddCut(AliDielectronVarManager::kKinkIndex0,   0.       );
	trackQuality->AddCut(AliDielectronVarManager::kNclsTPC,     70., 160. );
	trackQuality->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.,   4. );
	trackQuality->AddCut(AliDielectronVarManager::kPt,         1. ,  1e30);
	cuts->AddCut(trackQuality);
	
	
// PID cuts 
	
	
	AliDielectronPID *pid=new AliDielectronPID("PID","PID");
	pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,		-3.  ,3.0);
	//pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,		3.5, 1000);
	//pid->AddCut(AliDielectronPID::kTPC,AliPID::kProton, 		3. ,1000);

	cuts->AddCut(pid);

  //exclude conversion electrons selected by the tender
	//AliDielectronTrackCuts *noconv=new AliDielectronTrackCuts("noConv","noConv");
//	noconv->SetV0DaughterCut(AliPID::kElectron,kTRUE);
	//cuts->AddCut(noconv);
}


/**
*
* Setup the pair cuts according to the given cut definition
*
**/


void SetupPairCuts(AliDielectron *die, Int_t cutDefinition)
{ 
	
	
	
  // add conversion rejection
  AliDielectronVarCuts *gammaCut=new AliDielectronVarCuts("gammaCut","gammaCut");
  gammaCut->AddCut(AliDielectronVarManager::kM,0.,.1);
  die->GetPairPreFilter().AddCuts(gammaCut);
  die->SetPreFilterUnlikeOnly();
	
	
		
		
//		AliDielectronVarCuts *invMassCut=new AliDielectronVarCuts("InvMass","1.6<M");
//		invMassCut->AddCut(AliDielectronVarManager::kM,1.6,5.0);
//		die->GetPairFilter().AddCuts(invMassCut);


//		AliDielectronVarCuts *yCut=new AliDielectronVarCuts("y0.9Cut","y0.9Cut");
//		yCut->AddCut(AliDielectronVarManager::kY,-0.9,0.9);
//		die->GetPairFilter().AddCuts(yCut); 

		AliDielectronVarCuts *y09Cut=new AliDielectronVarCuts("y0.9Cut","y0.9Cut");
		y09Cut->AddCut(AliDielectronVarManager::kY,-.9, .9	);
		die->GetPairFilter().AddCuts(y09Cut); 

// require that one track p_T is above the trigger threshold
	if( cutDefinition & maskTrackPt  || cutDefinition & maskTriggered ){
		Double_t offlineThreshold = (cutDefinition & maskTrackPt23 ) ? 2.3 : 3.5;
		Double_t onlineThreshold = (cutDefinition & maskTriggeredQU ) ? 2. : 3.;
		Double_t pidThreshold = (cutDefinition & maskTriggeredQU ) ? 163.9 : 143.9;
		
		
		AliDielectronVarCuts* reqtrigger = new AliDielectronVarCuts("tracktriggered","tracktriggered");
		if(cutDefinition & maskTriggered){
			reqtrigger->AddCut(AliDielectronVarManager::kTRDonlinePt, -1 * onlineThreshold, onlineThreshold, kTRUE);
			reqtrigger->AddCut(AliDielectronVarManager::kTRDonlinePID, pidThreshold, 999.);
			reqtrigger->AddCut(AliDielectronVarManager::kTRDonlineNTracklets, 4.9,999.);
			reqtrigger->AddCut(AliDielectronVarManager::kTRDonlineFirstLayer, 0.1,999.,kTRUE);
		}
		if(cutDefinition & maskTrackPt){
			reqtrigger->AddCut(AliDielectronVarManager::kPt,offlineThreshold,999.);
		}
		AliDielectronPairLegCuts *triggerCut=new AliDielectronPairLegCuts("RequireTriggeredTrack","RequireTriggeredTrack");
		triggerCut->GetLeg1Filter().AddCuts(reqtrigger);
    	triggerCut->GetLeg2Filter().AddCuts(reqtrigger);
    	triggerCut->SetCutType(AliDielectronPairLegCuts::kAnyLeg);
		die->GetPairFilter().AddCuts(triggerCut);
	}



}

//______________________________________________________________________________________
void SetupV0Cuts(AliDielectron *die, Int_t cutDefinition)
{
  //
  // Setup the V0 cuts
  //

  
  AliDielectronV0Cuts *gammaV0Cuts = new AliDielectronV0Cuts("IsGamma","IsGamma");
  gammaV0Cuts->SetPdgCodes(22,11,11);
  gammaV0Cuts->SetDefaultPID(13); // TPC+-10, TOF+-3 if available
  gammaV0Cuts->AddCut(AliDielectronVarManager::kCosPointingAngle, TMath::Cos(0.02),   1.0, kFALSE);
  gammaV0Cuts->AddCut(AliDielectronVarManager::kChi2NDF,                       0.0,  10.0, kFALSE);//to be checked, if properly filled
  gammaV0Cuts->AddCut(AliDielectronVarManager::kLegDist,                       0.0,   0.25, kFALSE);
  gammaV0Cuts->AddCut(AliDielectronVarManager::kR,                             3.0,  90.0, kFALSE);
  gammaV0Cuts->AddCut(AliDielectronVarManager::kPsiPair,                       0.0,   0.05, kFALSE);
  gammaV0Cuts->AddCut(AliDielectronVarManager::kM,                             0.0,   0.05, kFALSE);
  //  gammaV0Cuts->AddCut(AliDielectronVarManager::kOpeningAngle,              0.0,   0.1, kFALSE);
  gammaV0Cuts->AddCut(AliDielectronVarManager::kArmPt,                         0.0,   0.05, kFALSE);
  //  gammaV0Cuts->AddCut(AliDielectronVarManager::kArmAlpha,                     -0.35,  0.35, kFALSE); // not sure if it works as expected
  gammaV0Cuts->SetExcludeTracks(kTRUE);
 // gammaV0Cuts->Print();
  
  //  const Double_t |cutAlphaG| < 0.35; &&  const Double_t cutQTG < 0.05;
 //  const Double_t |cutAlphaG2|[2] = {0.6, 0.8}; &&  const Double_t cutQTG2 < 0.04;

  die->GetTrackFilter().AddCuts(gammaV0Cuts);
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
  
  //Event class
  histos->AddClass("Event");

  //Track classes
  for (Int_t i=0; i<2; ++i){
    histos->AddClass(Form("Track_%s",AliDielectron::TrackClassName(i)));
  }


  //Pair classes
  // to fill also mixed event histograms loop until 10
 Int_t pairClasses [11] = {0,1,2,3,4,5,6,7,8,9, 10};
  
  
  for (Int_t i=0; i<10; ++i){
    histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(pairClasses[i])));
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
  
  
  //add histograms to event class

    
    //add histograms to event class
    histos->AddClass("Event");
    histos->UserHistogram("Event","","",2000,0.,2000.,AliDielectronVarManager::kNTrk);
    histos->UserHistogram("Event","","",100,0.,10.,AliDielectronVarManager::kPairs);
    histos->UserHistogram("Event","","",4,0.,4.,AliDielectronVarManager::kNevents);    
    histos->UserHistogram("Event","","",300,-15.,15.,AliDielectronVarManager::kZvPrim);
    histos->UserHistogram("Event","","",
                          100,-2.,2.,100,-2.,2.,AliDielectronVarManager::kXvPrim,AliDielectronVarManager::kYvPrim);
    // histos->UserHistogram("Event","","",AliDielectronHelper::MakeArbitraryBinning(run_numbers),AliDielectronVarManager::kRunNumber);
   

  //nsigma
  histos->UserHistogram("Track","","",
            200,0.2,20.,AliDielectronVarManager::kPt); 
  histos->UserHistogram("Track","","",80,-1.0,1.0,AliDielectronVarManager::kEta); 
  histos->UserHistogram("Track","","",400,0.0,4.0,AliDielectronVarManager::kPhi);  
  
  histos->UserHistogram("Track","","",400,0.0,4.0,80,-1.0,1.0,AliDielectronVarManager::kPhi,AliDielectronVarManager::kEta);  
  histos->UserHistogram("Track","","",
            200,0.2,20.,100,-10.,10.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kITSnSigmaEle);
  histos->UserHistogram("Track","","",
                        200,0.2,20.,100,-10.,10.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,kTRUE);
  histos->UserHistogram("Track","","",
                        200,0.2,20.,100,-10.,10.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaPio,kTRUE);
  histos->UserHistogram("Track","","",
                        200,0.2,20.,100,-10.,10.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaPro,kTRUE);
  histos->UserHistogram("Track","","",
            200,0.2,20.,100,-15.,15.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaKao,kTRUE);       
  histos->UserHistogram("Track","","",
                        100,0.2,20.,200,0.,200.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCsignal); 
  histos->UserHistogram("Track","","",
                        100,0.2,20.,200,0.,200.,AliDielectronVarManager::kP,AliDielectronVarManager::kTPCsignal); 
  histos->UserHistogram("Track","","",
                        100,0.2,20.,200,0.,200.,AliDielectronVarManager::kPt,AliDielectronVarManager::kTPCsignal); 
  histos->UserHistogram("Track","","",
            144,0.0,6.285,100,0.0,200,AliDielectronVarManager::kPhi,AliDielectronVarManager::kTPCsignal);
  histos->UserHistogram("Track","","",
            40,-1.0,1.0,100,0.0,200,AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCsignal);
  histos->UserHistogram("Track","","",
                        160,-0.5,159.5,AliDielectronVarManager::kNclsTPC);
  histos->UserHistogram("Track","","",
                        100,0.,10.,AliDielectronVarManager::kTPCchi2Cl);
  histos->UserHistogram("Track","","",
            160,-0.5,159.5,AliDielectronVarManager::kNFclsTPCr);
  histos->UserHistogram("Track","","",
            160,0,1.1,AliDielectronVarManager::kNFclsTPCrFrac);
 

// TRD online variable
  histos->UserHistogram("Track","","",256,0.,256.,AliDielectronVarManager::kTRDonlinePID );
  histos->UserHistogram("Track","","",200,-20.,20.,AliDielectronVarManager::kTRDonlinePt );
  histos->UserHistogram("Track","","",10,-1,9., AliDielectronVarManager::kTRDonlineNTracklets );
  histos->UserHistogram("Track","","",10,-1,9., AliDielectronVarManager::kTRDonlineFirstLayer );

 
  histos->UserHistogram("Track","","",
            20,-1.0,1.0,100,-10.,10.,AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEle);
  
  
  histos->UserHistogram("Track","","",
            20,0.0,7.0,100,-10.,10.,AliDielectronVarManager::kPhi,AliDielectronVarManager::kTPCnSigmaEle);
  
  
  histos->UserHistogram("Track","","",500,-1.,1.,AliDielectronVarManager::kImpactParXY);
histos->UserHistogram("Track","","",600,-3.,3.,AliDielectronVarManager::kImpactParZ);

  //add histograms to Pair classes

  histos->UserHistogram("Pair","","", 350,1.5,5.,AliDielectronVarManager::kM);
  histos->UserHistogram("Pair","","", 200,0.1,20.,AliDielectronVarManager::kPt);
  histos->UserHistogram("Pair","","", 100,-1.,1.,AliDielectronVarManager::kY);
  histos->UserHistogram("Pair","","", 100,0.,3.15,AliDielectronVarManager::kOpeningAngle);
  histos->UserHistogram("Pair","","", 350,1.5,5.,100,0.,3.15,AliDielectronVarManager::kM,AliDielectronVarManager::kOpeningAngle);
  histos->UserHistogram("Pair","","", 359,1.5,5.,200,0.2,20.,AliDielectronVarManager::kM,AliDielectronVarManager::kPt);


  die->SetHistogramManager(histos);
  
}

void InitCF(AliDielectron* die, Int_t cutDefinition)
{
  //
  // Setupd the CF Manager if needed
  //
  
  AliDielectronCF *cf=new AliDielectronCF(die->GetName(),die->GetTitle());
  
  //pair variables
  cf->AddVariable(AliDielectronVarManager::kPt,"0., .5,  1., 1.5, 2., 2.5,  3., 5., 7., 10., 100.");
  
  cf->AddVariable(AliDielectronVarManager::kZvPrim,"-10,-5,0,5,10");
  cf->AddVariable(AliDielectronVarManager::kY,"-1,-0.9,-0.8,-0.3,0.,0.3,0.8,0.9,1.0");
  cf->AddVariable(AliDielectronVarManager::kM,125,0.,125*.04); //40Mev Steps
  cf->AddVariable(AliDielectronVarManager::kPairType,11,0,11);

  //leg variables
  cf->AddVariable(AliDielectronVarManager::kPt,"0.0, 1.0, 2.3, 3.5, 100.0",kTRUE);
  cf->AddVariable(AliDielectronVarManager::kEta,"-1.,-0.9,-0.8,0,0.8,0.9,1.0",kTRUE);
  
  //cf->AddVariable(AliDielectronVarManager::kITSLayerFirstCls,6,0,6,kTRUE);

  
  if (hasMC){
    cf->AddVariable(AliDielectronVarManager::kPdgCode,10000,-5000.5,4999.5,kTRUE);
    cf->AddVariable(AliDielectronVarManager::kPdgCodeMother,10000,-5000.5,4999.5,kTRUE);
    cf->AddVariable(AliDielectronVarManager::kPdgCodeGrandMother,10000,-5000.5,4999.5,kTRUE);

    //only steps for efficiencies
    cf->SetStepsForMCtruthOnly();
  }
  
    cf->SetStepForMCtruth();
  
  cf->SetStepsForSignal();
  
  die->SetCFManagerPair(cf);
}

void SetupMCsignals(AliDielectron *die){

	
  AliDielectronSignalMC* inclusiveJpsi = new AliDielectronSignalMC("inclusiveJpsi","Inclusive J/psi");
  inclusiveJpsi->SetLegPDGs(11,-11);
  inclusiveJpsi->SetMotherPDGs(443,443);
  inclusiveJpsi->SetMothersRelation(AliDielectronSignalMC::kSame);
  inclusiveJpsi->SetFillPureMCStep(kTRUE);
  inclusiveJpsi->SetCheckBothChargesLegs(kTRUE,kTRUE);
  inclusiveJpsi->SetCheckBothChargesMothers(kTRUE,kTRUE);
  die->AddSignalMC(inclusiveJpsi);
	
	
	
  AliDielectronSignalMC* promptJpsi = new AliDielectronSignalMC("promptJpsi","Prompt J/psi");   // prompt J/psi (not from beauty decays)
  promptJpsi->SetLegPDGs(11,-11);
  promptJpsi->SetMotherPDGs(443,443);
  promptJpsi->SetGrandMotherPDGs(503,503,kTRUE,kTRUE);   // not from beauty hadrons
  promptJpsi->SetMothersRelation(AliDielectronSignalMC::kSame);
  promptJpsi->SetFillPureMCStep(kTRUE);
  promptJpsi->SetCheckBothChargesLegs(kTRUE,kTRUE);
  promptJpsi->SetCheckBothChargesMothers(kTRUE,kTRUE);
  promptJpsi->SetCheckBothChargesGrandMothers(kTRUE,kTRUE);
  die->AddSignalMC(promptJpsi);
	
	
  AliDielectronSignalMC* nonmpromptJpsi = new AliDielectronSignalMC("nonmpromptJpsi","Nonmprompt J/psi");   // nonmprompt J/psi (from beauty decays)
  nonmpromptJpsi->SetLegPDGs(11,-11);
  nonmpromptJpsi->SetMotherPDGs(443,443);
  nonmpromptJpsi->SetGrandMotherPDGs(503,503);   // from beauty hadrons
  nonmpromptJpsi->SetMothersRelation(AliDielectronSignalMC::kSame);
  nonmpromptJpsi->SetFillPureMCStep(kTRUE);
  nonmpromptJpsi->SetCheckBothChargesLegs(kTRUE,kTRUE);
  nonmpromptJpsi->SetCheckBothChargesMothers(kTRUE,kTRUE);
  nonmpromptJpsi->SetCheckBothChargesGrandMothers(kTRUE,kTRUE);
  die->AddSignalMC(nonmpromptJpsi);
  
}

void SetEtaCorrection()
{
  if (AliDielectronPID::GetEtaCorrFunction()) return;
  
  TString list=gSystem->Getenv("LIST");
  
  TFile f("$TRAIN_ROOT/jpsi_JPSI/EtaCorrMaps.root");
  if (!f.IsOpen()) return;
  TList *keys=f.GetListOfKeys();
  
  for (Int_t i=0; i<keys->GetEntries(); ++i){
    TString kName=keys->At(i)->GetName();
    TPRegexp reg(kName);
    if (reg.MatchB(list)){
      printf("Using Eta Correction Function: %s\n",kName.Data());
      AliDielectronPID::SetEtaCorrFunction((TF1*)f.Get(kName.Data()));
    }
  }
}
