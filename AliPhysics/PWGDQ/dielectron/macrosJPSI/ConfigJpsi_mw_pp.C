void InitHistograms(AliDielectron *die, Int_t cutDefinition);
void InitCF(AliDielectron* die, Int_t cutDefinition);

void SetupTrackCuts(AliDielectron *die, Int_t cutDefinition);
void SetupPairCuts(AliDielectron *die, Int_t cutDefinition);

AliESDtrackCuts *SetupESDtrackCuts(Int_t cutDefinition);

TString names=("trackQ+highpT+11dE/p(oneleg)+55<dEdx<120");
TObjArray *arrNames=names.Tokenize(";");

const Int_t nDie=arrNames->GetEntries();

AliDielectron* ConfigJpsi_mw_pp(Int_t cutDefinition)
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
  
  // cut setup
  SetupTrackCuts(die,cutDefinition);
  SetupPairCuts(die,cutDefinition);
  
  //
  // histogram setup
  // only if an AliDielectronHistos object is attached to the
  //  dielectron framework histograms will be filled
  //
  InitHistograms(die,cutDefinition);

  // init the debug tree. Use this only if you know what you are doing!!!
  InitDebugTree(die,cutDefinition);
  
  return die;
}

//______________________________________________________________________________________
void SetupTrackCuts(AliDielectron *die, Int_t cutDefinition)
{
  //
  // Setup the track cuts
  //

  AliDielectronCutGroup* cuts = new AliDielectronCutGroup("cuts","cuts",AliDielectronCutGroup::kCompAND);
  die->GetTrackFilter().AddCuts(cuts);
  
  //Pt + nsigma cut
  AliDielectronVarCuts *pt = new AliDielectronVarCuts("Pt>1.5+55<dEdx<120","Pt>1.5");
	
	  pt->AddCut(AliDielectronVarManager::kPt,1.5,1e30);
   //pt->AddCut(AliDielectronVarManager::kTPCsignal,70.,100.);//was before 60!!
	 pt->AddCut(AliDielectronVarManager::kNclsTPC,50.,160.);
		//here also cut for 0.9 eta range of legs
	 pt->AddCut(AliDielectronVarManager::kEta,-0.9,0.9);
	 
	 /*
  pt->AddCut(AliDielectronVarManager::kPt,.6,1e30);
  //  pt->AddCut(AliDielectronVarManager::kTPCsignal,60.,100.);
  pt->AddCut(AliDielectronVarManager::kTPCnSigmaEle,-3.,5.);
  */
	pt->AddCut(AliDielectronVarManager::kTPCsignal,55.,120.);//to be checked, usually rather 70, but after check in reduced trees and look at QA plots from GSI and consultation with Jens from 55 to 120 chosen, lower (seems to be also logical due to relative pull down of electron compared to pion line (MIP usually set to 50 by default)  by lower TPC gain )
  cuts->AddCut(pt);
  
  //ESD quality cuts
  cuts->AddCut(SetupESDtrackCuts(cutDefinition));
  
  //remove conversions tagged by the V0 tender supply
  AliDielectronTrackCuts *noconv=new AliDielectronTrackCuts("noConv","noConv");
  noconv->SetV0DaughterCut(AliPID::kElectron,kTRUE);
  cuts->AddCut(noconv);
	
  
}

//______________________________________________________________________________________
void SetupPairCuts(AliDielectron *die, Int_t cutDefinition)
{
  //
  // Setup the pair cuts
  //
  
   //
  // Setup the pair cuts
  //
  //Invarian mass selection
  AliDielectronVarCuts *invMassCut=new AliDielectronVarCuts("InvMass","1.6<M<5.0, pT>3.");
  // Minv > 1.8
  invMassCut->AddCut(AliDielectronVarManager::kM,1.6,5.0);
//invMassCut->AddCut(AliDielectronVarManager::kPairType,1.);
  // ptJpsi > 1GeV
  invMassCut->AddCut(AliDielectronVarManager::kPt,3.,1e30);//pT in the end also more stringent 5 or something like this
  die->GetPairFilter().AddCuts(invMassCut);

  AliDielectronVarCuts *gammaCut=new AliDielectronVarCuts("gammaCut","gammaCut");
  gammaCut->AddCut(AliDielectronVarManager::kM,0.,.05);// was first at 0.05
  die->GetPairPreFilter().AddCuts(gammaCut);

	//EMCAL cut only on one leg!
  AliDielectronPID *EMCAL = new AliDielectronPID("EMCAL","EMCAL nSigma e>-3. +e<4.");
	EMCAL->AddCut(AliDielectronPID::kEMCAL, AliPID::kElectron,-3.,4.);

	AliDielectronVarCuts *ptoneleg = new AliDielectronVarCuts("Pt>2.","Pt>2.");//finally perhaps 3.5, just more open for playing
  // pt > 0.7GeV
  ptoneleg->AddCut(AliDielectronVarManager::kPt,2.0,1e30);

	AliDielectronPairLegCuts *EMCALAny = new AliDielectronPairLegCuts("EMCALpid any","EMCALpid any");
	EMCALAny->GetLeg1Filter().AddCuts(EMCAL);
	EMCALAny->GetLeg1Filter().AddCuts(ptoneleg);
	EMCALAny->GetLeg2Filter().AddCuts(EMCAL);
	EMCALAny->GetLeg2Filter().AddCuts(ptoneleg);
	EMCALAny->SetCutType(AliDielectronPairLegCuts::kAnyLeg);
	die->GetPairFilter().AddCuts(EMCALAny);



  die->SetPreFilterUnlikeOnly();
}

//______________________________________________________________________________________
AliESDtrackCuts *SetupESDtrackCuts(Int_t cutDefinition)
{
  //
  // Setup default AliESDtrackCuts
  //
  AliESDtrackCuts *esdTrackCuts = new AliESDtrackCuts;
  
  esdTrackCuts->SetMaxDCAToVertexZ(3.0);
  esdTrackCuts->SetMaxDCAToVertexXY(1.0); 
  esdTrackCuts->SetRequireTPCRefit(kTRUE);
  //esdTrackCuts->SetRequireITSRefit(kTRUE);
  esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
//   esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
  
  esdTrackCuts->SetMinNClustersTPC(50);
  esdTrackCuts->SetMaxChi2PerClusterTPC(4);

  return esdTrackCuts;
}


//______________________________________________________________________________________
void InitHistograms(AliDielectron *die, Int_t cutDefinition)
{
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
  
  //add histograms to event class
  histos->UserHistogram("Event","VtxZ","Vertex Z;Z[cm]",300,-15.,15.,AliDielectronVarManager::kZvPrim);
  
  //add histograms to Track classes
  histos->UserHistogram("Track","Pt","Pt;Pt [GeV];#tracks",200,0,20.,AliDielectronVarManager::kPt);
  histos->UserHistogram("Track","dEdx_P","dEdx;P [GeV];TPC signal (arb units);#tracks",
                        200,0.2,20.,100,0.,200.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCsignal,kTRUE);
  histos->UserHistogram("Track","TPCnSigmaEle_P","TPC number of sigmas Electrons;P [GeV];TPC number of sigmas Electrons;#tracks",
                        200,0.2,20.,100,-10.,10.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,kTRUE);
  histos->UserHistogram("Track","TPCnSigmaKao_P","TPC number of sigmas Kaons;P [GeV];TPC number of sigmas Kaons;#tracks",
                        200,0.2,20.,100,-10.,10.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaKao,kTRUE);
  histos->UserHistogram("Track","TPCnSigmaPio_P","TPC number of sigmas Pions;P [GeV];TPC number of sigmas Pions;#tracks",
                        200,0.2,20.,100,-10.,10.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaPio,kTRUE);
  
  
  histos->UserHistogram("Track","TOFnSigmaEle_P","TOF number of sigmas Electrons;P [GeV];TOF number of sigmas Electrons;#tracks",
                        200,0.2,20.,100,-10.,10.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaEle,kTRUE);
  histos->UserHistogram("Track","TOFnSigmaKao_P","TOF number of sigmas Kaons;P [GeV];TOF number of sigmas Kaons;#tracks",
                        200,0.2,20.,100,-10.,10.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaKao,kTRUE);
  histos->UserHistogram("Track","TOFnSigmaPro_P","TOF number of sigmas Protons;P [GeV];TOF number of sigmas Protons;#tracks",
                        200,0.2,20.,100,-10.,10.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaPro,kTRUE);
												
histos->UserHistogram("Track","EmCalnSigmaEle_P","EMCAL number of sigmas Electrons;P [GeV];EMCAL number of sigmas Electrons;#tracks",
                        200,0.2,20.,100,-10.,10.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kEMCALnSigmaEle,kTRUE);
  
  histos->UserHistogram("Track","dXY","dXY;dXY [cm];#tracks",200,-2.,2.,AliDielectronVarManager::kImpactParXY);
  
  histos->UserHistogram("Track","TPCnCls","Number of Clusters TPC;TPC number clusteres;#tracks",159,0.,159.,AliDielectronVarManager::kNclsTPC);
  histos->UserHistogram("Track","TPCnClsIter1","Number of Clusters TPC;TPC number clusteres;#tracks",160,-0.5,159.5,AliDielectronVarManager::kNclsTPCiter1);
  histos->UserHistogram("Track","TPCsignalN","Number of Clusters TPC;TPC number clusteres;#tracks",160,-0.5,159.5,AliDielectronVarManager::kTPCsignalN);
  
  //add histograms to Pair classes
  histos->UserHistogram("Pair","InvMass","Inv.Mass;Inv. Mass [GeV];#pairs",
                        201,-.01,4.01,AliDielectronVarManager::kM);
  histos->UserHistogram("Pair","Rapidity","Rapidity;Rapidity;#pairs",
                        100,-2.,2.,AliDielectronVarManager::kY);
  histos->UserHistogram("Pair","OpeningAngle","Opening angle;angle",
                        100,0.,3.15,AliDielectronVarManager::kOpeningAngle);
  histos->UserHistogram("Pair","Chi2/NDF","#Chi^{2}/NDF;#Chi^{2}/NDF",
                        100, 0., 20., AliDielectronVarManager::kChi2NDF);
  
  die->SetHistogramManager(histos);
}

//______________________________________________________________________________________
void InitDebugTree(AliDielectron *die, Int_t cutDefinition)
{
  //
  // Create a debug tree manager
  // it should be only defined for ONE AliDielectron manager!!!
  
  //follwing lines taken out for grid compatibility and replaced by lines after this!
  /*
  AliDielectronDebugTree *tree=new AliDielectronDebugTree(Form("%sDebug",die->GetName()), "DielectronDebugTree");
  tree->SetOutputFileName(Form("jpsi_debug%02d.root",cutDefinition));
  TString addoutput=gSystem->Getenv("ADD_OUTPUT_FILES");
  if (addoutput.Length()) addoutput+=",";
  addoutput+=Form("jpsi_debug%02d.root",cutDefinition);
  gSystem->Setenv("ADD_OUTPUT_FILES",addoutput.Data());
  */

  AliDielectronDebugTree *tree=new AliDielectronDebugTree(Form("%sDebug",die->GetName()), "DielectronDebugTree");
  tree->SetOutputFileName("Jpsi_debugTree.root"); 
  //  tree->AddPairVariable(AliDielectronVarManager::kPx);
  //  tree->AddPairVariable(AliDielectronVarManager::kPy);
  //  tree->AddPairVariable(AliDielectronVarManager::kPz);
  tree->AddPairVariable(AliDielectronVarManager::kZvPrim);
  tree->AddPairVariable(AliDielectronVarManager::kPt);
  tree->AddPairVariable(AliDielectronVarManager::kY);
  tree->AddPairVariable(AliDielectronVarManager::kPhi);
  tree->AddPairVariable(AliDielectronVarManager::kM);
  tree->AddPairVariable(AliDielectronVarManager::kPairType);
  //  tree->AddPairVariable(AliDielectronVarManager::kLegDist);     // needs fix in KF package
  //  tree->AddPairVariable(AliDielectronVarManager::kLegDistXY);
  tree->AddPairVariable(AliDielectronVarManager::kOpeningAngle);
  tree->AddPairVariable(AliDielectronVarManager::kChi2NDF);
  tree->AddPairVariable(AliDielectronVarManager::kThetaCS);
  tree->AddPairVariable(AliDielectronVarManager::kThetaHE);
  
  //   leg variables
  //  tree->AddLegVariable(AliDielectronVarManager::kPx);
  //  tree->AddLegVariable(AliDielectronVarManager::kPy);
  //  tree->AddLegVariable(AliDielectronVarManager::kPz);
  tree->AddLegVariable(AliDielectronVarManager::kPt);
  tree->AddLegVariable(AliDielectronVarManager::kP);
  //  tree->AddLegVariable(AliDielectronVarManager::kE);
  //  tree->AddLegVariable(AliDielectronVarManager::kM);
  tree->AddLegVariable(AliDielectronVarManager::kEta);
  tree->AddLegVariable(AliDielectronVarManager::kPhi);
  tree->AddLegVariable(AliDielectronVarManager::kImpactParXY);
  tree->AddLegVariable(AliDielectronVarManager::kImpactParZ);
  
  tree->AddLegVariable(AliDielectronVarManager::kPIn);
  tree->AddLegVariable(AliDielectronVarManager::kTPCsignal);
  tree->AddLegVariable(AliDielectronVarManager::kTPCsignalN);
  tree->AddLegVariable(AliDielectronVarManager::kNclsTPC);
//   tree->AddLegVariable(AliDielectronVarManager::kNFclsTPCr);
//   tree->AddLegVariable(AliDielectronVarManager::kNFclsTPCrFrac);
  tree->AddLegVariable(AliDielectronVarManager::kTPCchi2Cl);
  tree->AddLegVariable(AliDielectronVarManager::kITSchi2Cl);
  tree->AddLegVariable(AliDielectronVarManager::kTrackStatus);
  
  tree->AddLegVariable(AliDielectronVarManager::kTPCnSigmaEle);
  tree->AddLegVariable(AliDielectronVarManager::kTPCnSigmaPio);
//   tree->AddLegVariable(AliDielectronVarManager::kTPCnSigmaMuo);
  tree->AddLegVariable(AliDielectronVarManager::kTPCnSigmaKao);
  tree->AddLegVariable(AliDielectronVarManager::kTPCnSigmaPro);

  tree->AddLegVariable(AliDielectronVarManager::kEMCALnSigmaEle);
  
  tree->AddLegVariable(AliDielectronVarManager::kKinkIndex0);
  
//   tree->AddLegVariable(AliDielectronVarManager::kITSclusterMap);
  tree->AddLegVariable(AliDielectronVarManager::kITSLayerFirstCls);
  
  tree->AddLegVariable(AliDielectronVarManager::kTOFnSigmaEle);
  //  tree->AddLegVariable(AliDielectronVarManager::kTOFnSigmaPio);
  //  tree->AddLegVariable(AliDielectronVarManager::kTOFnSigmaMuo);
  //  tree->AddLegVariable(AliDielectronVarManager::kTOFnSigmaKao);
  //  tree->AddLegVariable(AliDielectronVarManager::kTOFnSigmaPro);

  //  tree->AddLegVariable(AliDielectronVarManager::kTRDpidQuality);
  //  tree->AddLegVariable(AliDielectronVarManager::kTRDprobEle);
  //  tree->AddLegVariable(AliDielectronVarManager::kTRDprobPio);
  
  Bool_t hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);
  if (hasMC){
    tree->AddLegVariable(AliDielectronVarManager::kPdgCode);
    tree->AddLegVariable(AliDielectronVarManager::kPdgCodeMother);
    tree->AddLegVariable(AliDielectronVarManager::kPdgCodeGrandMother);
  }
  
  die->SetDebugTree(tree);
}


