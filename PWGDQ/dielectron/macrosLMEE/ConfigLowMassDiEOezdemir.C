
   void InitHistograms(AliDielectron *die, Int_t cutDefinition);
   void InitCF(AliDielectron* die, Int_t cutDefinition);

   void InitMCCuts(Bool_t mcStatus, Int_t nPDGC, TArrayI* pdgC);

   void SetupTrackCuts(AliDielectron *die, Int_t cutDefinition);

   void SetupPairCuts(AliDielectron *die, Int_t cutDefinition);

   void EnableMC();
   AliESDtrackCuts *SetupESDtrackCuts(Int_t cutDefinition);
   AliESDtrackCuts *SetupESDtrackCutsITSSA(Int_t cutDefinition);

   TString names=("Base050");
   TObjArray *arrNames=names.Tokenize(";");
   const Int_t nDie=arrNames->GetEntries();
   Bool_t MCenabled=kFALSE;

   Bool_t kMix = 0; //Event mixing (M.Koehler)

   AliDielectron* ConfigLowMassDiEOezdemir(Int_t cutDefinition,Bool_t hasMC=kFALSE)
   {
	 //
	 // Setup the instance of AliDielectron
	 //
	 
	 MCenabled=hasMC;
	 
	 // create the actual framework object
	
	TString name=Form("%02d",cutDefinition);
	 if ((cutDefinition)<arrNames->GetEntriesFast()){
	   name=arrNames->At((cutDefinition))->GetName();
	 }

	 //thisCut only relevant for MC:
	 AliDielectron *die =
	   new AliDielectron(Form
	   ("%s",name.Data()),
		   Form("Track cuts: %s",name.Data()));
	//Event mixing (M.Koehler)
	if(kMix){
 	AliDielectronMixingHandler *mix = new AliDielectronMixingHandler;
 	mix->SetMixType(AliDielectronMixingHandler::kAll);
 	mix->AddVariable(AliDielectronVarManager::kZvPrim,"-10., -7.5, -5., -2.5 , 0., 2.5, 5., 7.5 , 10.");
	mix->AddVariable(AliDielectronVarManager::kNacc,"0,10000");
 	mix->SetDepth(10);
 	die->SetMixingHandler(mix);
 	}
	 
	 // cut setup
	 SetupTrackCuts(die,cutDefinition);
	 SetupPairCuts(die,cutDefinition);
	 
	 // histogram setup
	 // only if an AliDielectronHistos object is attached to the
	 // dielectron framework histograms will be filled
	 //
	 InitHistograms(die,cutDefinition);

	 // the last definition uses no cuts and only the QA histograms should be filled!
	 InitCF(die,cutDefinition);
/*
	 AliDielectronTrackRotator *rot=new AliDielectronTrackRotator;
	 rot->SetIterations(1);
	 rot->SetConeAnglePhi(TMath::Pi());
	 die->SetTrackRotator(rot);
*/
	 return die;
   }

   //______________________________________________________________________________________
   void SetupTrackCuts(AliDielectron *die, Int_t cutDefinition)
   {
	 //
	 // Setup the track cuts
	 //
	 //


	 //
	 // ============ DALITZ REJECTION CUTS ====================
	 //
	 AliDielectronCutGroup* cgITSTPC = new AliDielectronCutGroup("cgITSTPC","cgITSTPC",AliDielectronCutGroup::kCompAND);
	 AliDielectronCutGroup* cgITSSA = new AliDielectronCutGroup("cgITSSA","cgITSSA",AliDielectronCutGroup::kCompAND);

	 cgITSTPC->AddCut(SetupESDtrackCuts(cutDefinition));
	 
	 AliDielectronPID *pidITSTPC = new AliDielectronPID("TPCpre","TPCpre");
	 pidITSTPC->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-3.,3.);
	 cgITSTPC->AddCut(pidITSTPC);

	 //---

	 cgITSSA->AddCut(SetupESDtrackCutsITSSA(cutDefinition));

	 AliDielectronPID *pidITSSA = new  AliDielectronPID("pidITSSA","pidITSSA");
	 //pidITSSA->AddCut(AliDielectronPID::kITS,AliPID::kPion,-999999.,-2.5,0.0,0.15,kFALSE ); //18.01.11: electron ID@lowpT
	 //pidITSSA->AddCut(AliDielectronPID::kITS,AliPID::kElectron,-3.,3.);
	 cgITSSA->AddCut(pidITSSA);

	 //
	 // ================Setup PID-Cuts (for final acceptance)===============
	 //
	 AliDielectronVarCuts *pt = new AliDielectronVarCuts("P>.2","P>.2");
	 pt->AddCut(AliDielectronVarManager::kP,.2,1e30);

	 // + [2] added for Dec2010 Cut! 
	 TF1 *lowerCut = new TF1("lowerCut", "[0] * TMath::Exp([1]*x) + [2]", 0, 20);
     /* until Nov2010
	 lowerCut->SetParameter(0, -2.7);
	 lowerCut->SetParameter(1, -0.4357);
	 */
	 /* 18.01.2011 ALiHFEpid.cxx */
	 lowerCut->SetParameter(0,-3.7);
	 lowerCut->SetParameter(1,-0.8);
	 lowerCut->SetParameter(2,-0.35);
/*
	 if (MCenabled) { //overwrite parameters
	   lowerCut->SetParameter(0,-2.5);
	   lowerCut->SetParameter(2,-2.2);
	 }
*/
	 //___________________________________________

	 AliDielectronPID *pidTPChardTOF = new AliDielectronPID("TPC-TOF-HFE","TPC-TOF-HFE");
	 pidTPChardTOF->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,lowerCut,3.,0.0,100.,kFALSE);
	 pidTPChardTOF->AddCut(AliDielectronPID::kTOF ,AliPID::kElectron , -3. , 3. , 0.0 , 100., kFALSE );

	 //___________________________________________
	 AliDielectronPID *pidTT = new AliDielectronPID("TPC-TOF-HFE","TPC-TOF-HFE");
	 pidTT->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,lowerCut,3.,0.0,100.,kFALSE);
	 pidTT->AddCut( AliDielectronPID::kTOF ,AliPID::kElectron , -3. , 3. , 0.4 , 100., kFALSE );
	 pidTT->AddCut(AliDielectronPID::kTPC,AliPID::kPion,-3,3.,0.,0.4,kTRUE);
	 pidTT->AddCut(AliDielectronPID::kTPC,AliPID::kProton,-3,3.,0.,0.4,kTRUE); 
	 pidTT->AddCut(AliDielectronPID::kTPC,AliPID::kKaon,-3,3.,0.,0.4,kTRUE); 

	//___________________________________________
	 AliDielectronPID *pidTTnew = new AliDielectronPID("TPC-TOF-HFE","TPC-TOF-HFE");
	 pidTTnew->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-1.5,3.,0.2,100.,kFALSE);
	 pidTTnew->AddCut( AliDielectronPID::kTOF ,AliPID::kElectron , -3. , 3. , 0.4 , 5., kFALSE );
	 pidTTnew->AddCut(AliDielectronPID::kTPC,AliPID::kPion,-100.,4.,0.2,100.,kTRUE);
	 //pidTTnew->AddCut(AliDielectronPID::kTPC,AliPID::kProton,-3.,3.,0.2,0.4,kTRUE); 
	 //pidTTnew->AddCut(AliDielectronPID::kTPC,AliPID::kKaon,-3.,3.,0.2,0.4,kTRUE); 

	 //___________________________________________
	 AliDielectronPID *pidTTonlyKP = new AliDielectronPID("TPC-TOF-HFE","TPC-TOF-HFE");
	 pidTTonlyKP->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,lowerCut,3.,0.0,100.,kFALSE);
	 pidTTonlyKP->AddCut( AliDielectronPID::kTOF ,AliPID::kElectron , -3. , 3. , 0.4 , 100., kFALSE );
	 pidTTonlyKP->AddCut(AliDielectronPID::kTPC,AliPID::kProton,-3,3.,0.,0.4,kTRUE); 
	 pidTTonlyKP->AddCut(AliDielectronPID::kTPC,AliPID::kKaon,-3,3.,0.,0.4,kTRUE); 

	 //___________________________________________

	 AliDielectronPID *pidTTupp = new AliDielectronPID("TPC-TOF-HFE","TPC-TOF-HFE");
	 pidTTupp->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-1.,3.,0.0,100.,kFALSE);
	 pidTTupp->AddCut( AliDielectronPID::kTOF ,AliPID::kElectron , -1. , 3. , 0.4 , 100., kFALSE );
	 pidTTupp->AddCut(AliDielectronPID::kTPC,AliPID::kPion,-1,3.,0.,0.4,kTRUE);
	 pidTTupp->AddCut(AliDielectronPID::kTPC,AliPID::kProton,-1,3.,0.,0.4,kTRUE); 
	 pidTTupp->AddCut(AliDielectronPID::kTPC,AliPID::kKaon,-1,3.,0.,0.4,kTRUE); 
	 //___________________________________________

	 AliDielectronPID *pidITS = new AliDielectronPID("ITS PID","ITS PID");
	 //pidITS->AddCut(AliDielectronPID::kITS,AliPID::kElectron,-1.,1.,0.0,100.,kFALSE );
	 pidITS->AddCut(AliDielectronPID::kITS,AliPID::kPion,-999999.,-2.5,0.0,0.15,kFALSE ); //18.01.11: electron ID@lowpT
	 pidITS->AddCut(AliDielectronPID::kITS,AliPID::kPion,-2.,2.,0.15,0.3,kTRUE ); //18.01.11: electron ID@lowpT
	 //move to config section, to give additional point in CF
	 pidITS->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,lowerCut,3.,0.0,100.,kFALSE);
	 pidITS->AddCut(AliDielectronPID::kTPC,AliPID::kPion,-3,3.,0.0,100.,kTRUE);
	 pidITS->AddCut(AliDielectronPID::kTOF ,AliPID::kElectron , -3. , 3. , 0.3 , 100., kFALSE );

	 AliDielectronVarCuts *ptITS = new AliDielectronVarCuts("Pt>.05","Pt>.05");
	 ptITS->AddCut(AliDielectronVarManager::kPt,.05,1e30);
	 //___________________________________________

	 AliDielectronPID *pidJpsi = new AliDielectronPID("TPC Jpsi-like","TPC Jpsi-like");
	 pidJpsi->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,lowerCut,3.,kFALSE); 
	 pidJpsi->AddCut(AliDielectronPID::kTPC,AliPID::kProton,-3,3.,kTRUE); //Difference: HFE-like inclusion
	 pidJpsi->AddCut(AliDielectronPID::kTPC,AliPID::kPion,-3,3.,kTRUE); //Difference: HFE-like inclusion
	 pidJpsi->AddCut(AliDielectronPID::kTPC,AliPID::kKaon,-3,3.,kTRUE); //Difference: HFE-like inclusion

         //___________________________________________
	 AliDielectronPID *pidTPConlyE = new AliDielectronPID("TPC PID","TPC PID");
	 pidTPConlyE->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-1.5,3.,0.,100.,kFALSE);

	//___________________________________________
	 AliDielectronPID *pidTPC = new AliDielectronPID("TPC PID","TPC PID");
	 pidTPC->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-1.5,3.,0.,100.,kFALSE);
	 pidTPC->AddCut(AliDielectronPID::kTPC,AliPID::kProton,-1.5,1.5,0.,100.,kTRUE); 
	 pidTPC->AddCut(AliDielectronPID::kTPC,AliPID::kKaon,-1.5,1.5,0.,100.,kTRUE);
	 pidTPC->AddCut(AliDielectronPID::kTPC,AliPID::kPion,-1.5,1.5,0.,100.,kTRUE); 

	 // 
	 // =============CREATE GROUPS FOR ==============
	 //

	 //Soft cut conditions for initial single-track analysis
	 AliDielectronCutGroup* cgInitialTrackFilter = new AliDielectronCutGroup("cgInitialTrackFilter","cgInitialTrackFilter",AliDielectronCutGroup::kCompOR);
	 cgInitialTrackFilter->AddCut(cgITSTPC);
	 cgInitialTrackFilter->AddCut(cgITSSA);


	 //Setup possible CutGroups for Step2
	 //HFE-Like:
	 //
	 //ITS-TPC-TOF
	 AliDielectronCutGroup* cgSecondTrackFilterITS = new AliDielectronCutGroup("cgSecondTrackFilterITS","cgSecondTrackFilterITS",AliDielectronCutGroup::kCompAND);
	 cgSecondTrackFilterITS->AddCut(cgITSTPC);
	 cgSecondTrackFilterITS->AddCut(pidITS);
	 cgSecondTrackFilterITS->AddCut(ptITS);


	 //Config 6: TPC-TOF classic rej. only kaon & proton
	 AliDielectronCutGroup* cgSecondTrackFilterPIDTPConlyKP = new AliDielectronCutGroup("cgPIDTPConlyKP","cgPIDTPConlyKP",AliDielectronCutGroup::kCompAND);
	 cgSecondTrackFilterPIDTPConlyKP->AddCut(pt);
	 cgSecondTrackFilterPIDTPConlyKP->AddCut(SetupESDtrackCuts(cutDefinition));
	 cgSecondTrackFilterPIDTPConlyKP->AddCut(pidTTonlyKP);


	 //Config 0: TPC-TOF classic 
	 AliDielectronCutGroup* cgSecondTrackFilterPIDTPCclassic = new AliDielectronCutGroup("cgPIDTPCclassic","cgPIDTPCclassic",AliDielectronCutGroup::kCompAND);
	 cgSecondTrackFilterPIDTPCclassic->AddCut(pt);
	 cgSecondTrackFilterPIDTPCclassic->AddCut(SetupESDtrackCuts(cutDefinition));
	 cgSecondTrackFilterPIDTPCclassic->AddCut(pidTT);

	 //Config 1: TPC-TOF Classic with -1.5/3 sigma (soft TOF, TPC exclusion)
	 AliDielectronCutGroup* cgSecondTrackFilterPIDTPC = new AliDielectronCutGroup("cgPIDTPC","cgPIDTPC",AliDielectronCutGroup::kCompAND);
	 //cgSecondTrackFilterPIDTPC->AddCut(pt);
	 cgSecondTrackFilterPIDTPC->AddCut(SetupESDtrackCuts(cutDefinition));
	 cgSecondTrackFilterPIDTPC->AddCut(pidTTnew);

	//Config 2: TPC only electrons
	 AliDielectronCutGroup* cgSecondTrackFilterPIDTPConlyE = new AliDielectronCutGroup("cgPIDTPConlyE","cgPIDTPConlyE",AliDielectronCutGroup::kCompAND);
	 cgSecondTrackFilterPIDTPConlyE->AddCut(SetupESDtrackCuts(cutDefinition));
	 cgSecondTrackFilterPIDTPConlyE->AddCut(pidTPConlyE);

	//TPC (kaon, proton, pion exclusion)
	AliDielectronCutGroup* cgSecondTrackFilterPIDTPCSA = new AliDielectronCutGroup("cgPIDTPCSA","cgPIDTPCSA",AliDielectronCutGroup::kCompAND);
	 cgSecondTrackFilterPIDTPCSA->AddCut(SetupESDtrackCuts(cutDefinition));
	 cgSecondTrackFilterPIDTPCSA->AddCut(pidTPC);

	 //Not in use
	 AliDielectronCutGroup* cgSecondTrackFilterJpsi = new AliDielectronCutGroup("cgSecondTrackFilterJpsi","cgSecondTrackFilterJpsi",AliDielectronCutGroup::kCompAND);
	 cgSecondTrackFilterJpsi->AddCut(pt);
	 cgSecondTrackFilterJpsi->AddCut(pidJpsi);
	 cgSecondTrackFilterJpsi->AddCut(cgITSTPC);

	 AliDielectronCutGroup* cgSecondTrackFilterCombinedTPC=  new AliDielectronCutGroup("cgSecondTrackFilterCombinedTPC","cgSecondTrackFilterCombinedTPC",AliDielectronCutGroup::kCompOR);
	 cgSecondTrackFilterCombinedTPC->AddCut(cgSecondTrackFilterPIDTPC);
	 cgSecondTrackFilterCombinedTPC->AddCut(cgSecondTrackFilterJpsi);



	 //Correct BG scheme for Combinatorial BG 
	 //die->SetPreFilterUnlikeOnly();
	 die->SetPreFilterAllSigns();
	 die->SetUseKF(kFALSE);

	 //New Cut-scheme
	 //
	 //
	/*
	 if ( cutDefinition==0 ) {
	   //Step1: TrackFilter, nur grob
	   die->GetTrackFilter().AddCuts(cgInitialTrackFilter);
	   //Step2: "Final" PID
	   die->GetPairPreFilterLegs().AddCuts(cgSecondTrackFilterPIDTPC);
	 
	 }
	 
	 if ( cutDefinition==1 ) {
	   //Step1: TrackFilter, nur grob
	   die->GetTrackFilter().AddCuts(cgInitialTrackFilter);
	   //Step2: "Final" PID
	   die->GetPairPreFilterLegs().AddCuts(cgSecondTrackFilterPIDTPConlyKP);
	 
	 }
	
	if ( cutDefinition==2 ) {
	   //Step1: TrackFilter, nur grob
	   die->GetTrackFilter().AddCuts(cgInitialTrackFilter);
	   //Step2: "Final" PID
	   die->GetPairPreFilterLegs().AddCuts(cgSecondTrackFilterPIDTPConlyE);
	 
	 }
	 
	 if ( cutDefinition==3 ) {
	   //Step1: TrackFilter, nur grob
	   die->GetTrackFilter().AddCuts(cgInitialTrackFilter);
	   //Step2: "Final" PID
	   die->GetPairPreFilterLegs().AddCuts(cgSecondTrackFilterPIDTPCSA);
	 
	 }
        

	if ( cutDefinition==4 ) {
	   //Step1: TrackFilter, nur grob
	   die->GetTrackFilter().AddCuts(cgInitialTrackFilter);
	   //Step2: "Final" PID
	   die->GetPairPreFilterLegs().AddCuts(cgSecondTrackFilterPIDTPC);
	 
	 }
         */
	
	if ( cutDefinition ==0 ) {
	   //Step1: TrackFilter, nur grob
	   die->GetTrackFilter().AddCuts(cgInitialTrackFilter);
	   //Step2: "Final" PID
	   die->GetPairPreFilterLegs().AddCuts(cgSecondTrackFilterPIDTPC);	
	}
	/*
	 if ( ((cutDefinition>=0) && (cutDefinition <=2) ) || (cutDefinition ==4)|| (cutDefinition ==7)|| (cutDefinition ==8)|| (cutDefinition ==9)) {
	   AliDielectronVarCuts *vtxZ = new AliDielectronVarCuts("vtxZ","Vertex z cut");
	   vtxZ->AddCut(AliDielectronVarManager::kZvPrim,-15.,15.);
	   die->GetEventFilter().AddCuts(vtxZ);
	 }
	*/

   }
   //______________________________________________________________________________________

   void SetupPairCuts(AliDielectron *die, Int_t cutDefinition)
   {
	 //
	 // Setup the pair cuts
	 //
	 
	 /*
	 //Invariant mass selection
	 if (cutDefinition==1) {
		 AliDielectronVarCuts *dalitzpairReject=new AliDielectronVarCuts("InvMass","InvMass > 150 MeV");
		 dalitzpairReject->AddCut(AliDielectronVarManager::kM,0.15,100.,kTRUE);
		 die->GetPairPreFilter().AddCuts(dalitzpairReject);
	   }
	*/
	/*
	if (cutDefinition==0){
		AliDielectronVarCuts *openingAnglecut = new AliDielectronVarCuts("OpeningAngle","OpeningAngle > 0.035 rad");
		openingAnglecut->AddCut(AliDielectronVarManager::kOpeningAngle, 0.,0.035);
		die->GetPairPreFilter().AddCuts(openingAnglecut);
	}

	if (cutDefinition==4) {
		 AliDielectronVarCuts *dalitzpairReject=new AliDielectronVarCuts("InvMass","InvMass > 150 MeV");
		 dalitzpairReject->AddCut(AliDielectronVarManager::kM,0.15,100.,kTRUE);
		 die->GetPairPreFilter().AddCuts(dalitzpairReject);
	   }
	*/
	
	if(cutDefinition == 0){

	AliDielectronTrackCuts *noconv=new AliDielectronTrackCuts("noConv","conversion tagging");
	noconv->SetV0DaughterCut(AliPID::kElectron,kTRUE);
 	die->GetTrackFilter().AddCuts(noconv);
	            
	//pair rapidity
        AliDielectronVarCuts *RapidityCut=new AliDielectronVarCuts("RapidityCut","RapidityCut");
	RapidityCut->AddCut(AliDielectronVarManager::kY, -0.8 , 0.8);
        die->GetPairFilter().AddCuts(RapidityCut);

        AliDielectronVarCuts *PhiV = new AliDielectronVarCuts("PhiV","PhiV");//mass and Phiv together
        PhiV->AddCut(AliDielectronVarManager::kM, 0.0 , 0.05);
        PhiV->AddCut(AliDielectronVarManager::kPhivPair, 2.5 , 3.2 );
	die->GetPairPreFilter().AddCuts(PhiV);
	/*
        AliDielectronVarCuts *PhiV2 = new AliDielectronVarCuts("PhiV2","PhiV2");//mass and Phiv together
        PhiV2->AddCut(AliDielectronVarManager::kM, 0.05 , 1000.);

        AliDielectronCutGroup* pairCutsCG2 =new AliDielectronCutGroup("pairCutsCG2","pairCutsCG2",AliDielectronCutGroup::kCompOR);
        pairCutsCG2->AddCut(PhiV);
        pairCutsCG2->AddCut(PhiV2);
        die->GetPairFilter().AddCuts(pairCutsCG2);
	*/
		}


   }

   //______________________________________________________________________________________
   AliESDtrackCuts *SetupESDtrackCuts(Int_t cutDefinition)
   {
	 //
	 // Setup default AliESDtrackCuts
	 //
	 AliESDtrackCuts *esdTrackCuts = new AliESDtrackCuts;
	  //esdTrackCuts->SetDCAToVertex2D(kTRUE);
	 
	 esdTrackCuts->SetMaxDCAToVertexZ(3.0);
	 esdTrackCuts->SetMaxDCAToVertexXY(1.0); 
	 //esdTrackCuts->SetMaxNsigmaToVertex(3.0); 

	 esdTrackCuts->SetPtRange(  0.2 , 100.0);
	 esdTrackCuts->SetEtaRange( -0.8 , 0.8 );
	 esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
	 esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
	 esdTrackCuts->SetDCAToVertex2D(kFALSE);

	 esdTrackCuts->SetRequireITSRefit(kTRUE);
	 //esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
	 esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
	 esdTrackCuts->SetMinNClustersITS(3);
	/*
	 if(cutDefinition == 1) esdTrackCuts->SetMinNClustersITS(4); 
	 else esdTrackCuts->SetMinNClustersITS(3);
	*/
	 esdTrackCuts->SetRequireTPCRefit(kTRUE);
	 //esdTrackCuts->SetMinNClustersTPC(60);
	 esdTrackCuts->SetMinNClustersTPC(80);
	/*
	 if(cutDefinition == 2) esdTrackCuts->SetMinNClustersTPC(100);
 	 else esdTrackCuts->SetMinNClustersTPC(80);
	*/
	 //esdTrackCuts->SetMinNClustersTPC(110);
	 //esdTrackCuts->SetMinNClustersTPC(70);
	 esdTrackCuts->SetMinNCrossedRowsTPC(100);
	/*
	 if(cutDefinition == 3) esdTrackCuts->SetMinNCrossedRowsTPC(120);
         else esdTrackCuts->SetMinNCrossedRowsTPC(100);
	*/
	 esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.5);
	 esdTrackCuts->SetMaxChi2PerClusterTPC(4);

	 return esdTrackCuts;
   }

   //______________________________________________________________________________________

   AliESDtrackCuts *SetupESDtrackCutsITSSA(Int_t cutDefinition)
   {
	 //
	 // Setup default AliESDtrackCuts
	 //
	 AliESDtrackCuts *esdTrackCuts = new AliESDtrackCuts;
	  esdTrackCuts->SetDCAToVertex2D(kTRUE);

	 esdTrackCuts->SetMaxDCAToVertexZ(3.0);
	 esdTrackCuts->SetMaxDCAToVertexXY(1.0); 

	 esdTrackCuts->SetEtaRange( -0.84 , 0.84 );
	 
	 esdTrackCuts->SetPtRange(  0.0, 200.);

	 esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
	 esdTrackCuts->SetRequireITSRefit(kTRUE);
	 esdTrackCuts->SetRequireITSStandAlone(kTRUE);

	 esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
	 esdTrackCuts->SetMinNClustersITS(3); //PhotonGroup-Dalitz: 2?!
	 

	 return esdTrackCuts;
   }

   //______________________________________________________________________________________

void InitHistograms(AliDielectron *die, Int_t cutDefinition)
   {
/*
	 //
	 // Initialise the histograms
	   //
	 
	 //Setup histogram Manager
	 AliDielectronHistos *histos=
	   new AliDielectronHistos(die->GetName(),
							   die->GetTitle());
	 
	 //Initialise histogram classes
	 histos->SetReservedWords("Track;Pair;Pre;RejTrack;RejPair");
	 
	 //Track classes
	 //to fill also track info from 2nd event loop until 2
	 for (Int_t i=0; i<2; ++i){
	   histos->AddClass(Form("Track_%s",AliDielectron::TrackClassName(i)));
	 }
	 
	 //Pair classes
	 // to fill also mixed event histograms loop until 10
	 for (Int_t i=0; i<3; ++i){
	   histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(i)));
	 }
		 
	 //PreFilter Classes
	 //to fill also track info from 2nd event loop until 2
	 for (Int_t i=0; i<2; ++i){
	   histos->AddClass(Form("Pre_%s",AliDielectron::TrackClassName(i)));
	 }


	 //Create Classes for Rejected Tracks/Pairs:
	 for (Int_t i=0; i<3; ++i){
	   histos->AddClass(Form("RejTrack_%s",AliDielectron::PairClassName(i)));
	 }
	 for (Int_t i=0; i<3; ++i){
	   histos->AddClass(Form("RejPair_%s",AliDielectron::PairClassName(i)));
	 }
	
	 
	 //Event class
	 histos->AddClass("Event");
	 //if (cutDefinition==nDie-1) histos->AddClass("Event");
	*/ 
         /*	
	 //track rotation
	 histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(AliDielectron::kEv1PMRot)));
	 histos->AddClass(Form("Track_Legs_%s",AliDielectron::PairClassName(AliDielectron::kEv1PMRot)));
	 */
/*
	 //Event mixing (M.Koehler)
	 if(kMix){
 	  histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(3))); //ME ++
 	  histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(4)));//ME -+
 	  histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(6)));//ME +-
 	  histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(7))); // ME --
 	  }	
	
	 if (cutDefinition==nDie-1){
	   //add histograms to event class
	 }

	 histos->UserHistogram("Event","nEvents","Number of processed events after cuts;Number events",
							 1,0.,1.,AliDielectronVarManager::kNevents);
	   histos->UserHistogram("Event","Nacc","Number of accepted tracks;Number tracks after cuts",
							 1,0.,1.,AliDielectronVarManager::kNacc);
	   histos->UserHistogram("Event","NTrk","Number of accepted ;Number tracks ",
							 1,0.,1.,AliDielectronVarManager::kNTrk);
  
  //add histograms to Track classes
  histos->UserHistogram("Track","Pt","Pt;Pt [GeV];#tracks",200,0,20.,AliDielectronVarManager::kPt);
  

  histos->UserHistogram("Track","ITS_dEdx_P","ITS_dEdx;P [GeV];ITS signal (arb units);#tracks",
                        400,0.0,20.,1000,0.,1000.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kITSsignal,kTRUE);

  histos->UserHistogram("Track","dEdx_P","dEdx;P [GeV];TPC signal (arb units);#tracks",
                        400,0.0,20.,200,0.,200.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCsignal,kTRUE);


  histos->UserHistogram("Track","TPCnSigmaEle_P","TPC number of sigmas Electrons;P [GeV];TPC number of sigmas Electrons;#tracks",
                        400,0.0,20.,100,-5.,5.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,kTRUE);
  histos->UserHistogram("Track","TPCnSigmaKao_P","TPC number of sigmas Kaons;P [GeV];TPC number of sigmas Kaons;#tracks",
                        400,0.0,20.,100,-5.,5.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaKao,kTRUE);
  histos->UserHistogram("Track","TPCnSigmaPio_P","TPC number of sigmas Pions;P [GeV];TPC number of sigmas Pions;#tracks",
                        400,0.0,20.,100,-5.,5.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaPio,kTRUE);
  
  histos->UserHistogram("Track","TRDpidPobEle_P","TRD PID probability Electrons;P [GeV];TRD prob Electrons;#tracks",
                        400,0.0,20.,100,0.,1.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTRDprobEle,kTRUE);
  histos->UserHistogram("Track","TRDpidPobPio_P","TRD PID probability Pions;P [GeV];TRD prob Pions;#tracks",
                        400,0.0,20.,100,0.,1.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTRDprobPio,kTRUE);
  
  histos->UserHistogram("Track","TOFnSigmaKao_P","TOF number of sigmas Kaons;P [GeV];TOF number of sigmas Kaons;#tracks",
                        400,0.0,20.,100,-5.,5.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaKao,kTRUE);
  histos->UserHistogram("Track","TOFnSigmaPro_P","TOF number of sigmas Protons;P [GeV];TOF number of sigmas Protons;#tracks",
                        400,0.0,20.,100,-5.,5.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaPro,kTRUE);

  histos->UserHistogram("Track","Eta_Phi","Eta Phi Map; Eta; Phi;#tracks",
                        200,-2,2,200,0,3.15,AliDielectronVarManager::kEta,AliDielectronVarManager::kPhi,kTRUE);
  
  histos->UserHistogram("Track","dXY","dXY;dXY [cm];#tracks",200,-2.,2.,AliDielectronVarManager::kImpactParXY,kTRUE);
  
  histos->UserHistogram("Track","TPCnCls","Number of Clusters TPC;TPC number clusteres;#tracks",159,0.,159.,AliDielectronVarManager::kNclsTPC,kTRUE);

  histos->UserHistogram("Track","TPCnCls_kNFclsTPCr","nTPC vs nTPCr;nTPC vs nTPCr;#tracks",159,0.,159.,159,0.,159.,AliDielectronVarManager::kNclsTPC,AliDielectronVarManager::kNFclsTPCr);

  histos->UserHistogram("Track","kNFclsTPCr_pT","nTPCr vs pt;nTPCr vs pt;#tracks",159,0.,159.,200,0.,20.,AliDielectronVarManager::kNFclsTPCr,AliDielectronVarManager::kPt);
  
  //add histograms to Pair classes
  histos->UserHistogram("Pair","InvMass","Inv.Mass;Inv. Mass [GeV];#pairs",
                         500,0.0,5.00,AliDielectronVarManager::kM);
  histos->UserHistogram("Pair","Rapidity","Rapidity;Rapidity;#pairs",
                        100,-2.,2.,AliDielectronVarManager::kY);
  histos->UserHistogram("Pair","OpeningAngle","Opening angle;angle",
                        100,0.,3.15,AliDielectronVarManager::kOpeningAngle);
  //2D Histo Plot
  histos->UserHistogram("Pair","InvMassPairPt","Inv.Mass vs PairPt;Inv. Mass [GeV], pT [GeV];#pairs",
                         500,0.0,5.0,500,0.,50.,AliDielectronVarManager::kM,AliDielectronVarManager::kPt);

  histos->UserHistogram("Pair","InvMassOpeningAngle","Opening Angle vs Inv.Mass;Inv. Mass [GeV];#pairs",
                         500,0.0,5.0,200,0.,6.3,AliDielectronVarManager::kM,AliDielectronVarManager::kOpeningAngle);
 
  
  die->SetHistogramManager(histos);
*/

 //
  // Initialise the histograms
  //
  
  //Setup histogram classes
  AliDielectronHistos *histos=
    new AliDielectronHistos(die->GetName(),
                            die->GetTitle());
  


  //Initialise histogram classes
  histos->SetReservedWords("Track;Pair");
  
  //Event class
  histos->AddClass("Event");
  

  //Track classes
  //to fill also track info from 2nd event loop until 2
  for (Int_t i=0; i<2; ++i){
    histos->AddClass(Form("Track_%s",AliDielectron::TrackClassName(i)));
  }
  
  //Pair classes
  // to fill also mixed event histograms loop until 10

  for (Int_t i=0; i<3; ++i){
    histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(i)));

  }

  if(kMix){
  histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(3))); //ME ++
  histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(4)));//ME -+
  histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(6)));//ME +-
  histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(7))); // ME --
  }

  //if(kRot)histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(10)));//Rot


  //add histograms to event class
  histos->UserHistogram("Event","nEvents","Number of processed events after cuts;Number events",1,0,1,AliDielectronVarManager::kNevents);
  histos->UserHistogram("Event","ZVertex","ZVertex;ZVertex/cm",480,-12.,12.,AliDielectronVarManager::kZvPrim);
  histos->UserHistogram("Event","nESDTracks","ESD tracks;ESD tracks;Number events",100,0,200,AliDielectronVarManager::kTracks);
  histos->UserHistogram("Event","Nacc","Number of accepted tracks;Number events",100,0,200,AliDielectronVarManager::kNacc);
  histos->UserHistogram("Event","NVtxContrib","Number of Vertex Contributor;NVtx;Number events",100,0,100,AliDielectronVarManager::kNVtxContrib);


  //add histograms to Track classes
  histos->UserHistogram("Track","P","P;P [GeV];#tracks",500,0.,10.,AliDielectronVarManager::kPIn);
  histos->UserHistogram("Track","Pt","Pt;Pt [GeV];#tracks",500,0.,10.,AliDielectronVarManager::kPt);
  histos->UserHistogram("Track","Eta","Eta; Eta ;#tracks",100,-1.,1.,AliDielectronVarManager::kEta);
  histos->UserHistogram("Track","Phi","Phi; Phi ;#tracks",640,0.,6.4,AliDielectronVarManager::kPhi);


  histos->UserHistogram("Track","Eta_phi","Eta vs Phi;Eta;Phi",100,-1.,1.,320,0.,6.4,AliDielectronVarManager::kEta,AliDielectronVarManager::kPhi);
  histos->UserHistogram("Track","Eta_pt","Eta vs Pt;Eta;Pt",100,-1.,1.,500,0.,10.,AliDielectronVarManager::kEta,AliDielectronVarManager::kPt);
  histos->UserHistogram("Track","ImpParXY_Pt","ImpParXY_Pt; Pt; ÍmpParXY",500,0.,10.,500,-5.,5.,AliDielectronVarManager::kPt,AliDielectronVarManager::kImpactParXY);
  histos->UserHistogram("Track","ImpParZ_Pt","ImpParZ_Pt; Pt; ÍmpParZ",500,0.,10.,500,-5.,5.,AliDielectronVarManager::kPt,AliDielectronVarManager::kImpactParZ);
 

  //track checks (ITS)
  histos->UserHistogram("Track","ITSchi2Cl_Mom","ITS Chi2 vs Momentum;Mom;ITS chi2",500,0.,5.,50,0.,5.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kITSchi2Cl);
  histos->UserHistogram("Track","NclsITS_Mom",";Mom;kNclsITS",500,0.,5.,7,0,7,AliDielectronVarManager::kPIn,AliDielectronVarManager::kNclsITS);

  //track checks (TPC)
  histos->UserHistogram("Track","TPCsignalNfrac_Mom",";fraction TPCSignalN/TPCncls vs Momentum;Mom;TPCSignalN/TPCncls",500,0.,5.,60,0.,1.2,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCsignalNfrac);
  histos->UserHistogram("Track","TPCchi2Cl_Mom","TPC Chi2 vs Momentum;Mom;TPC Chi2",500,0.,10.,100,0,5,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCchi2Cl);
  histos->UserHistogram("Track","TPCclsDiff_Mom","kTPCclsDiff vs Momentum;Mom;kTPCclsDiff",500,0.,10.,100,-10,10,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCclsDiff);
  histos->UserHistogram("Track","FindableTPCcls_Mom","kNFclsTPC vs Momentum;Mom;kNFclsTPC",500,0.,10.,200,0.,200.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kNFclsTPC);
  histos->UserHistogram("Track","TPCcls_Mom","kNclsTPC vs Momentum;Mom;kNclsTPC",500,0.,10.,200,0.,200.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kNclsTPC);
  histos->UserHistogram("Track","kNclsSFracTPC_Mom","kNclsSFracTPC vs Momentum;Mom;kTPCclsSFrac",500,0.,10.,1000,0.,1.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kNclsSFracTPC);
  histos->UserHistogram("Track","kNFclsTPCrFrac_Mom","kNFclsTPCrFrac vs Momentum;Mom;kNFclsTPCrFrac",500,0.,10.,60,0.,1.2.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kNFclsTPCrFrac);

  //track checks (TOF)
  histos->UserHistogram("Track","TOFbeta_Mom","kTOFbeta vs Momentum;Mom;TOFbeta",500,0.,5.,120,0.,1.2,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFbeta);
  histos->UserHistogram("Track","TOFPIDBit_Mom","kTOFPIDBit vs Momentum;Mom;TOFPIDbit",500,0.,5.,2,0.,1.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFPIDBit);

  //track checks (PID)
  histos->UserHistogram("Track","ITSnSigma_MomPio","ITS number of sigmas Pion vs Momentum;Mom;ITSsigmaPion",500,0.,5.,1000,-20,20,AliDielectronVarManager::kPIn,AliDielectronVarManager::kITSnSigmaPio);
  histos->UserHistogram("Track","TPCnSigma_MomPio","TPC number of sigmas Pions vs Momentum;Mom;TPCsigmaPion",500,0.,5.,1000,-20,20,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaPio);

  histos->UserHistogram("Track","ITSnSigma_MomEle","ITS number of sigmas Electrons vs Momentum;Mom;ITSsigmaEle",500,0.,5.,800,-40,40,AliDielectronVarManager::kPIn,AliDielectronVarManager::kITSnSigmaEle);
  histos->UserHistogram("Track","TPCnSigma_MomEle","TPC number of sigmas Electrons vs Momentum;Mom;TPCsigmaEle",1000,0.,10.,800,-40,40,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle);
  histos->UserHistogram("Track","TOFnSigma_MomEle","TOF number of sigmas Electrons vs Momentum;Mom;TOFsigmaEle",500,0.,5.,800,-40,40,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaEle);

  histos->UserHistogram("Track","ITSdEdx_P","dEdx;P [GeV];ITS signal (arb units) vs Momentum;Mom;ITSsignal", 500,0,5,800,0,200,AliDielectronVarManager::kPIn,AliDielectronVarManager::kITSsignal);
  histos->UserHistogram("Track","TPCdEdx_P","dEdx;P [GeV];TPC signal (arb units) vs Momentum;Mom;TPCsignal", 500,0,10.,800,0,200,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCsignal);

  //
  //add histograms to Pair classes
  //

  histos->UserHistogram("Pair","InvMass","Inv.Mass;Inv. Mass [GeV];#pairs",
                        500,0,4,AliDielectronVarManager::kM);

  histos->UserHistogram("Pair","InvMass_low","Inv.Mass;Inv. Mass [GeV];#pairs",
                        500,0,0.5,AliDielectronVarManager::kM);

  histos->UserHistogram("Pair","InvMass10","Inv.Mass;Inv. Mass [GeV];#pairs",
                        500,0.,5.,AliDielectronVarManager::kM);

  histos->UserHistogram("Pair","kDeltaEta","kDeltaEta;kDeltaEta;#pairs",
                        160,0.,1.6,AliDielectronVarManager::kDeltaEta);

  histos->UserHistogram("Pair","kDeltaEta_low","kDeltaEta;kDeltaEta;#pairs",
                        500,0.,0.5,AliDielectronVarManager::kDeltaEta);

  histos->UserHistogram("Pair","kDeltaPhi","kDeltaPhi;kDeltaPhi;#pairs",
                        320,0.,6.4,AliDielectronVarManager::kDeltaPhi);

  histos->UserHistogram("Pair",
                        "kDeltaEta_kDeltaPhi","kDeltaEta_kDeltaPhi;kDeltaEta;kDeltaPhi",
                        160, 0. , 1.6, 320 , 0., 6.4 ,
                         AliDielectronVarManager::kDeltaEta , AliDielectronVarManager::kDeltaPhi );

  histos->UserHistogram("Pair","PhiV",";PhiV;#pairs",
                        320,0.,6.4,AliDielectronVarManager::kPhivPair);


   histos->UserHistogram("Pair","PhiV_Pt",";Pt;PhiV",
			 100,0.,10.,320,0.,6.4,AliDielectronVarManager::kPt,AliDielectronVarManager::kPhivPair);

   histos->UserHistogram("Pair","InvMass_bin2","Inv.Mass;Inv. Mass [GeV];#pairs",
	"0. , 0.025, 0.05 , 0.075 ,0.1 ,0.15 , 0.2 , 0.25 , 0.3 , 
 	0.4 ,  0.5 , 0.6, 0.65 , 0.688 , 0.725, 0.75, 0.775, 0.8 , 0.85 ,
	 0.95,  0.975 , 1.0 , 1.025 , 1.05, 1.125 , 1.25 , 1.5 , 1.75 , 2.0 , 
	2.25, 2.5 , 2.75 , 2.85, 2.95,3.05, 3.1 , 3.15 , 
	3.3 , 3.5, 3.75 , 4.0",AliDielectronVarManager::kM);


  histos->UserHistogram("Pair",
                        "InvMass_Pt","InvMass_Pt;InvMass;Pt",
                        500, 0. , 4., 100 , 0., 5. ,
                         AliDielectronVarManager::kM , AliDielectronVarManager::kPt );

  histos->UserHistogram("Pair",
                        "InvMass_PhivPair","InvMass_PhivPair;InvMass;PhivPair",
                        500, 0. , 4., 320 , 0., 3.2 ,
                         AliDielectronVarManager::kM , AliDielectronVarManager::kPhivPair );

  histos->UserHistogram("Pair",
                        "OpAngle","Opening angle;Opening angle;#pairs",
                        320, 0. , 3.2, 
                         AliDielectronVarManager::kOpeningAngle);

  histos->UserHistogram("Pair",
                        "OpAngle_InvMass","OpAngle_InvMass;Opening angle;Invariant Mass",
                        320, 0. , 3.2, 500 , 0. , 4. ,
                         AliDielectronVarManager::kOpeningAngle,AliDielectronVarManager::kM);


  histos->UserHistogram("Pair",
                        "Phi","Phi;counts;Phi",
                        320, 0. , 6.4, AliDielectronVarManager::kPhi);

  histos->UserHistogram("Pair",
                        "Y","Y;counts;Y",
                        120, -1.2 , 1.2, AliDielectronVarManager::kY);

  die->SetHistogramManager(histos);

}


void InitCF(AliDielectron* die, Int_t cutDefinition)
{
  //
  // Setupd the CF Manager if needed
  //
  AliDielectronCF *cf=new AliDielectronCF(die->GetName(),die->GetTitle());
 
  //pair variables
  cf->AddVariable(AliDielectronVarManager::kPt,100,0,10);
  //cf->AddVariable(AliDielectronVarManager::kP,200,0,20);
  //cf->AddVariable(AliDielectronVarManager::kPhi,64, -3.2, 3.2);
  cf->AddVariable(AliDielectronVarManager::kY,40,-2,2);
  cf->AddVariable(AliDielectronVarManager::kM,500,0.,4.); 
  //cf->AddVariable(AliDielectronVarManager::kPairType,10,0,10);
  //cf->AddVariable(AliDielectronVarManager::kOpeningAngle,315,0,3.15);
  //cf->AddVariable(AliDielectronVarManager::kDeltaEta,200,-2,2);
  //cf->AddVariable(AliDielectronVarManager::kDeltaPhi,100,0,3.15);
  //cf->AddVariable(AliDielectronVarManager::kHaveSameMother,21,-10,10);
  cf->AddVariable(AliDielectronVarManager::kNumberOfDaughters,5,0,5);
  //leg variables
  cf->AddVariable(AliDielectronVarManager::kPt,100,0.,10.,kTRUE);
  //cf->AddVariable(AliDielectronVarManager::kP,200,0.,20.,kTRUE);
  //cf->AddVariable(AliDielectronVarManager::kY,40,-2.,2.,kTRUE);
  cf->AddVariable(AliDielectronVarManager::kEta,20,-1.,1.,kTRUE);
  cf->AddVariable(AliDielectronVarManager::kPhi,100,0.,3.15,kTRUE);
//  cf->AddVariable(AliDielectronVarManager::kITSsignal,1000,0.0.,1000.,kTRUE);
//  cf->AddVariable(AliDielectronVarManager::kHaveSameMother,21,-10,10,kTRUE);
//  cf->AddVariable(AliDielectronVarManager::kNumberOfDaughters,11,0,10,kTRUE);
 
  //only in this case write MC truth info
 
  if (MCenabled) {
	cf->SetStepForMCtruth();
	//cf->SetStepsForMCtruthOnly();
	//cf->SetStepForNoCutsMCmotherPid();
	cout << "MC ENABLED ------------------------------------------------------" << endl;
    cf->AddVariable(AliDielectronVarManager::kPdgCode,10000,-5000.5,4999.5,kTRUE);
    cf->AddVariable(AliDielectronVarManager::kPdgCodeMother,10000,-5000.5,4999.5,kTRUE);
  }

  //cf->SetStepsForEachCut();
  //cf->SetStepForPreFilter();
  cf->SetStepForAfterAllCuts();
  //cf->SetStepsForBackground();
  cf->SetStepsForSignal();
  
  die->SetCFManagerPair(cf);

/*

  AliDielectronSignalMC* lowMassDiele=new
	AliDielectronSignalMC("lowMassDiele","low mass dielectron pairs");
  lowMassDiele->SetLegPDGs(11,-11);
  lowMassDiele->SetCheckBothChargesLegs(kTRUE,kTRUE);
  lowMassDiele->SetLegSources(AliDielectronSignalMC::kPrimary,
	  AliDielectronSignalMC::kPrimary);
  lowMassDiele->SetFillPureMCStep(kTRUE);
  die->AddSignalMC(lowMassDiele);

  AliDielectronSignalMC* secondary=new
	AliDielectronSignalMC("secondary","secondary electrons pairs");
  secondary->SetLegPDGs(11,-11);
  secondary->SetCheckBothChargesLegs(kTRUE,kTRUE);
  secondary->SetLegSources(AliDielectronSignalMC::kSecondary,
	  AliDielectronSignalMC::kSecondary);
  die->AddSignalMC(secondary);

  AliDielectronSignalMC* eleFromConversions=new
	AliDielectronSignalMC("eleFromConversions","conversion electrons");
  eleFromConversions->SetLegPDGs(11,-11);
  eleFromConversions->SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleFromConversions->SetLegSources(AliDielectronSignalMC::kSecondary,
	  AliDielectronSignalMC::kSecondary);
  eleFromConversions->SetMotherPDGs(22,22);    // 22- photon
  die->AddSignalMC(eleFromConversions);

  AliDielectronSignalMC* misIdPions=new
	AliDielectronSignalMC("misIdPions","mis id. pion pairs");
  misIdPions->SetLegPDGs(211,-211);
  misIdPions->SetCheckBothChargesLegs(kTRUE,kTRUE);
  die->AddSignalMC(misIdPions);

  AliDielectronSignalMC* dalitzDecays=new
	AliDielectronSignalMC("dalitzDecays","dalitz Pairs");
  dalitzDecays->SetLegPDGs(11,-11);
  dalitzDecays->SetCheckBothChargesLegs(kTRUE,kTRUE);
  dalitzDecays->SetLegSources(AliDielectronSignalMC::kSecondary,
	        AliDielectronSignalMC::kSecondary);
  dalitzDecays->SetMotherPDGs(111,111);
  dalitzDecays->SetFillPureMCStep(kTRUE);
  die->AddSignalMC(dalitzDecays);
*/

  AliDielectronSignalMC* PhiDecays=new
	AliDielectronSignalMC("PhiDecays","Phi Pairs");
  PhiDecays->SetLegPDGs(11,-11);
  PhiDecays->SetCheckBothChargesLegs(kTRUE,kTRUE);
  PhiDecays->SetLegSources(AliDielectronSignalMC::kFinalState,
	        AliDielectronSignalMC::kFinalState);
  PhiDecays->SetMotherPDGs(333,333);
  PhiDecays->SetMothersRelation(AliDielectronSignalMC::kSame); 
  PhiDecays->SetFillPureMCStep(kTRUE);
  die->AddSignalMC(PhiDecays);

  AliDielectronSignalMC* OmegaDecays=new
	AliDielectronSignalMC("OmegaDecays","Omega Pairs");
  OmegaDecays->SetLegPDGs(11,-11);
  OmegaDecays->SetCheckBothChargesLegs(kTRUE,kTRUE);
  OmegaDecays->SetLegSources(AliDielectronSignalMC::kFinalState,
	        AliDielectronSignalMC::kFinalState);
  OmegaDecays->SetMotherPDGs(223,223);
  OmegaDecays->SetMothersRelation(AliDielectronSignalMC::kSame); 
  OmegaDecays->SetFillPureMCStep(kTRUE);
  die->AddSignalMC(OmegaDecays);
/*
AliDielectronSignalMC* DieleConti= new AliDielectronSignalMC("DieleConti","low mass dielectron pairs");
DieleConti->SetLegPDGs(11,-11);
DieleConti->SetMotherPDGs(0,0,22,22);
DieleConti->SetCheckBothChargesLegs(kTRUE,kTRUE);
DieleConti->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
DieleConti->SetMothersRelation(AliDielectronSignalMC::kSame);
DieleConti->SetFillPureMCStep(kTRUE);
die->AddSignalMC(DieleConti);
*/

}

//--------------------------------------
   void EnableMC() {
	MCenabled=kTRUE;
   }

