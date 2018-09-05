void InitHistograms(AliDielectron *die, Int_t cutDefinition);
void SetupCuts(AliDielectron *die, Int_t cutDefinition, TString userPathWeightFile, Double_t userTMVACutValue);
void SetupAODtrackCutsTMVAPIDFirst(AliDielectron *die, Bool_t bMCPID, Bool_t bTMVAPID, TString userPathWeightFile, Double_t userTMVACutValue);
void SetTPCCorr(AliDielectron *die);
const AliDielectronEventCuts *GetEventCuts();
void SetupMCsignals(AliDielectron* die);

Bool_t isRandomRejTask=kFALSE;//needed for InitHistograms() //dont change!!!
Bool_t kRot = kFALSE;
Bool_t kMix = kTRUE;
Bool_t kNoPairing   = kTRUE;
Bool_t randomizeDau = kTRUE;

// available cut defintions
const Int_t nMax = 3; 
const Int_t nPF  = 999; // use prefiltering for cuts > nPF

AliDielectron* Config_miweber_LMEE_TMVA(Double_t userTMVACutValue = 0.0,
					Int_t cutDefinition=1,
					TString userPathWeightFile="",
					Bool_t bESDANA = kFALSE,
					Bool_t bCutQA = kFALSE,
					Bool_t isRandomRej=kFALSE,
					Bool_t useTPCCorr=kFALSE)
{
  //
  // Setup the instance of AliDielectron
  // (inspired from simpler task of Theo Broeker w/o the usage of a CutLib file)
  //

  isRandomRejTask=isRandomRej;

  // First check if cut definition is defined
  if(!(cutDefinition >= 0 && cutDefinition < nMax))){
    Printf("=====================================");
    Printf("Cut set = %d is not available.",cutDefinition);
    Printf("Use [%d,%d[ for cuts",0,nMax);
    Printf("=====================================");
    return NULL;
  }
  else{
    Printf("=====================================");
    Printf("Configuring cut set %d",cutDefinition);
    Printf("=====================================");
  }
    

  // create the actual framework object
  TString name=Form("TMVACuts%d",cutDefinition);
  
  AliDielectron *die =
    new AliDielectron(Form("%s",name.Data()),
                      Form("Track cuts: %s",name.Data()));

  if(bCutQA){
    die->SetCutQA(bCutQA);
  }
  
  if(kRot){
    AliDielectronTrackRotator *rot = new AliDielectronTrackRotator;
    rot->SetConeAnglePhi(TMath::Pi());
    rot->SetIterations(10);
    die->SetTrackRotator(rot);
  }//kRot
  
  
  if(kMix && !(die->GetHasMC()) ){ // need second since there is a problem when mixing MC events (TRef?)
    AliDielectronMixingHandler *mix = new AliDielectronMixingHandler;

    mix->AddVariable(AliDielectronVarManager::kZvPrim,"-10,-5,0,5,10");
    mix->AddVariable(AliDielectronVarManager::kCentrality,"0,5,10,20,30,50,80");
    mix->SetDepth(15);
    mix->SetMixType(AliDielectronMixingHandler::kAll);
    
    // using TPC event plane, uncorrected. (also, the old phi range was wrong, now same effective binning.)
    // mix->AddVariable(AliDielectronVarManager::kQnTPCrpH2, 6, TMath::Pi()/-2., TMath::Pi()/2.);
        
    die->SetMixingHandler(mix);
  }//kMix


  // set track cuts
  SetupCuts(die,cutDefinition,bESDANA,userPathWeightFile,userTMVACutValue);
  
  if(useTPCCorr)  SetTPCCorr(die);
 
 //
 // histogram setup
 // only if an AliDielectronHistos object is attached to the
 // dielectron framework histograms will be filled
 //

 //
 // MC signals setup
 // SetupMCsignals(die);
 //

 InitHistograms(die,cutDefinition);
 //  InitCF(die,cutDefinition);
 
 die->SetNoPairing(kNoPairing);
 
 return die;

}

//______________________________________________________________________________________
void SetupCuts(AliDielectron *die, Int_t cutDefinition, Bool_t bESDANA = kFALSE, TString userPathWeightFile = "", Double_t userTMVACutValue)
{
  // Setup the track cuts

  //pairing with TLorentzVector
  die->SetUseKF(kFALSE);

  if( cutDefinition == 0 ){
  // first set of cuts to test TMVA usage for PID (only MC PID)
    if(!bESDANA){
      SetupAODtrackCutsTMVAPIDFirst(die,kTRUE,kFALSE,userPathWeightFile,userTMVACutValue);
    }
  }
  else if( cutDefinition == 1 ){
  // first set of cuts to test TMVA usage for PID (MC + TMVA PID)
    if(!bESDANA){      
      SetupAODtrackCutsTMVAPIDFirst(die,kTRUE,kTRUE,userPathWeightFile,userTMVACutValue);
    }
  }
  else if( cutDefinition == 2 ){
  // first set of cuts to test TMVA usage for PID (only TMVA PID)
    if(!bESDANA){      
      SetupAODtrackCutsTMVAPIDFirst(die,kFALSE,kTRUE,userPathWeightFile,userTMVACutValue);
    }
  }
}


//______________________________________________________________________________________
void SetupAODtrackCutsTMVAPIDFirst(AliDielectron *die, Bool_t bMCPID, Bool_t bTMVAPID, TString userPathWeightFile, Double_t userTMVACutValue)
{
  //
  // Setup the track cuts
  // - these are similar cuts used for nano AOD creation LHC15o (ConfigLMEE_nano_PbPb2015.C), but tighter
  // - Add TMVA usage for PID (first try)
  // 

  Bool_t hasMC=die->GetHasMC();

  // Quality cuts
  AliDielectronCutGroup* cuts = new AliDielectronCutGroup("cuts","cuts",AliDielectronCutGroup::kCompAND);

  /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv FILTER CUTS vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
  // AOD track filter (needs to be first cut to speed up)
  AliDielectronTrackCuts *trkFilter = new AliDielectronTrackCuts("TrkFilter","TrkFilter");
  trkFilter->SetAODFilterBit(AliDielectronTrackCuts::kTPCqualSPDany);

  /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv TRACK CUTS vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
  AliDielectronVarCuts *varCuts    = new AliDielectronVarCuts("VarCuts","VarCuts");
  AliDielectronVarCuts *varCuts2   = new AliDielectronVarCuts("VarCuts2","VarCuts2");
  AliDielectronTrackCuts *trkCuts  = new AliDielectronTrackCuts("TrkCuts","TrkCuts");
  
  // specific cuts
  trkCuts->SetITSclusterCut(AliDielectronTrackCuts::kOneOf, 3); // SPD any
  trkCuts->SetRequireITSRefit(kTRUE);
  trkCuts->SetRequireTPCRefit(kTRUE); // not useful when using prefilter

  // standard cuts
  varCuts->AddCut(AliDielectronVarManager::kNclsTPC,      80.0, 160.0);
  varCuts->AddCut(AliDielectronVarManager::kNclsITS,      3.0, 100.0);
  varCuts->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   15.0);
  varCuts->AddCut(AliDielectronVarManager::kNclsSITS,     0.0,   3.1); // means 0 and 1 shared Cluster
  varCuts->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   8.0);
  varCuts->AddCut(AliDielectronVarManager::kNFclsTPCr,    80.0, 160.0);
  varCuts->AddCut(AliDielectronVarManager::kPt,           0.2, 8.);
  varCuts->AddCut(AliDielectronVarManager::kEta,         -0.8,   0.8);
  varCuts->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
  varCuts->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
  varCuts->AddCut(AliDielectronVarManager::kKinkIndex0,   0.);

  // additional 
  varCuts2->AddCut(AliDielectronVarManager::kNFclsTPCfCross, 0.95, 1.05);

  // /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv PID CUTS (TMVA) vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */

  // Copy xml file to local directory first
  TString weightFile="alien:///alice/cern.ch/user/m/miweber/TMVA_weights/TMVA_e_0_05_05_05_electron0555first.weights.xml";
  if(userPathWeightFile!="")
    weightFile = userPathWeightFile;
  Printf("Use TMVA weight input file from Alien: %s",weightFile.Data());
  
  AliDielectronTMVACuts *pidCuts = new AliDielectronTMVACuts("PIDCutsTMVA","PIDCutsTMVA");
  pidCuts->AddTMVAInput("pt", AliDielectronVarManager::kPt);
  pidCuts->AddTMVAInput("nsigITSe", AliDielectronVarManager::kITSnSigmaEle);
  pidCuts->AddTMVAInput("nsigTPCe", AliDielectronVarManager::kTPCnSigmaEle);
  pidCuts->AddTMVAInput("nsigTOFe", AliDielectronVarManager::kTOFnSigmaEle);
  pidCuts->AddTMVAInput("nsigITSpi", AliDielectronVarManager::kITSnSigmaPio);
  pidCuts->AddTMVAInput("nsigTPCpi", AliDielectronVarManager::kTPCnSigmaPio);
  pidCuts->AddTMVAInput("nsigTOFpi", AliDielectronVarManager::kTOFnSigmaPio);
  pidCuts->AddTMVASpectator("pdg", AliDielectronVarManager::kPdgCode);
  pidCuts->SetTMVAWeights("BDT method", weightFile.Data());

  Printf("Use TMVA cut value = %f",userTMVACutValue);
  pidCuts->SetTMVACutValue(userTMVACutValue);

  /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv MC PID CUTS vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
  AliDielectronVarCuts *pidCutsMC = new AliDielectronVarCuts("PIDCutsMC","PIDCutsMC");
  pidCutsMC->SetCutType(AliDielectronVarCuts::kAny);
  //pidCutsMC->SetCutOnMCtruth(kTRUE);//un-comment to use only positive labels
  pidCutsMC->AddCut(AliDielectronVarManager::kPdgCode, -11.1, -10.9);
  pidCutsMC->AddCut(AliDielectronVarManager::kPdgCode, +10.9, +11.1);

  // activate the cut sets (order might be CPU timewise important)
  if(bMCPID){
    cuts->AddCut(pidCutsMC); 
    die->GetTrackFilter().AddCuts(pidCutsMC);
  }
  cuts->AddCut(trkFilter);
  die->GetTrackFilter().AddCuts(trkFilter);
  cuts->AddCut(varCuts);
  die->GetTrackFilter().AddCuts(varCuts);
  cuts->AddCut(trkCuts);
  die->GetTrackFilter().AddCuts(trkCuts);
  cuts->AddCut(varCuts2);
  die->GetTrackFilter().AddCuts(varCuts2);
  if(bTMVAPID){
    cuts->AddCut(pidCuts);
    die->GetTrackFilter().AddCuts(pidCuts);
  }
  cuts->Print();

}

//______________________________________________________________________________________
void SetTPCCorr(AliDielectron *die){
  ::Info("Config_miweber_LMEE_TMVA","starting LMEECutLib::SetEtaCorrectionTPC()\n");
  TString path="alien:///alice/cern.ch/user/s/selehner/recal/recalib_data_tpc_nsigmaele.root";
  gSystem->Exec(TString::Format("alien_cp %s .",path.Data()));
  ::Info("Config_miweber_LMEE_TMVA","Copy TPC correction from Alien: %s",path.Data());
  _file = TFile::Open("recalib_data_tpc_nsigmaele.root");
  
  TH3D* mean = dynamic_cast<TH3D*>(_file->Get("sum_mean_correction"));
  TH3D* width= dynamic_cast<TH3D*>(_file->Get("sum_width_correction"));

  if(mean)   ::Info("Config_miweber_LMEE_TMVA","Mean Correction Histo loaded, entries: %f",mean->GetEntries());
  else {
    ::Info("Config_miweber_LMEE_TMVA","Mean Correction Histo not loaded! entries: %f",mean->GetEntries());
    return 0;
  }

    die->SetCentroidCorrFunction(mean, AliDielectronVarManager::kP, AliDielectronVarManager::kEta, AliDielectronVarManager::kRefMultTPConly);
    die->SetWidthCorrFunction(width, AliDielectronVarManager::kP, AliDielectronVarManager::kEta, AliDielectronVarManager::kRefMultTPConly);
       
}

//______________________________________________________________________________________
void InitHistograms(AliDielectron *die, Int_t cutDefinition)
{
  
  //Setup histogram classes
  AliDielectronHistos *histos=
    new AliDielectronHistos(die->GetName(),
                           die->GetTitle());
  
  //Initialise histogram classes
  //histos->SetReservedWords("Track;Pair");
  histos->SetReservedWords("Track;Pair;Track_Legs;Pre;RejPair;RejTrack;Random");
  //histos->SetReservedWords("Track");  

  //Event class
  histos->AddClass("Event");

  if(!isRandomRejTask){
    //Track classes
    //to fill also track info from 2nd event loop until 2
    for (Int_t i=0; i<2; ++i){
      histos->AddClass(Form("Track_%s",AliDielectron::TrackClassName(i)));
    }
    //Pair classes
    for (Int_t i=0; i<3; ++i){
      histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(i)));
      // Legs of final Pairs. Both charges together. No duplicate entries.
      if(cutDefinition > nPF) 
	histos->AddClass(Form("Track_Legs_%s",AliDielectron::PairClassName(i))); // not TrackClassName, see 'AliDielectron::FillHistograms(...)'
    }
    //ME and track rot
    if (die->GetMixingHandler()) {
      histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(3)));
      histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(4)));
      histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(6)));
      histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(7)));
    }
    if (die->GetTrackRotator()) {
      histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(10)));
    }
    if(cutDefinition > nPF){  
      //PreFilter Classes
      //to fill also track info from 2nd event loop until 2
      for (Int_t i=0; i<2; ++i){
        histos->AddClass(Form("Pre_%s",AliDielectron::TrackClassName(i)));
      }
      
      //Create Classes for Rejected Tracks/Pairs:
      for (Int_t i=0; i<3; ++i){
        histos->AddClass(Form("RejPair_%s",AliDielectron::PairClassName(i)));
        // Legs of rejected Pairs. Both charges together. One track can and will make multiple entries.
        histos->AddClass(Form("RejTrack_%s",AliDielectron::PairClassName(i))); // not TrackClassName, see 'AliDielectron::FillHistogramsPair(...)'
      }
    }
  }
  if(isRandomRejTask){
    //
    // _____ histograms for AliAnalysisTaskMultiDielectronPR _____
    //
    //    histos->AddClass("Rand_Pair");
    //    histos->AddClass("Rand_RejPair");
    const char* cRandomPairClassNames[2] = { "Testpart", "RejTestpart" };
    for (Int_t i=0; i<2; ++i){
      histos->AddClass(Form("Random_%s",cRandomPairClassNames[i]));
    }
    histos->UserHistogram("Random","Pt","",200,0,10.,AliDielectronVarManager::kPt);
    histos->UserHistogram("Random","Eta","",200,-2,2,AliDielectronVarManager::kEta);
    histos->UserHistogram("Random","Phi","",120,0.,TMath::TwoPi(),AliDielectronVarManager::kPhi);
    histos->UserHistogram("Random","Px","",200,0,10.,AliDielectronVarManager::kPx);
    histos->UserHistogram("Random","Py","",200,0,10.,AliDielectronVarManager::kPy);
    histos->UserHistogram("Random","Pz","",200,0,10.,AliDielectronVarManager::kPz);
    histos->UserHistogram("Random","Pt_Eta_Phi","",
                          500,0.,10.,16,-0.8,0.8,30,0.,2*TMath::Pi(),
                          AliDielectronVarManager::kPt,AliDielectronVarManager::kEta,AliDielectronVarManager::kPhi);
  }

  //add MC signal histograms to pair class
  if(die->GetMCSignals()) {
    for (Int_t i=0; i<die->GetMCSignals()->GetEntriesFast(); ++i) {
      histos->AddClass(Form("Pair_%s",die->GetMCSignals()->At(i)->GetName()));
      // histos->AddClass(Form("Track_%s",die->GetMCSignals()->At(i)->GetName()));
      // histos->AddClass(Form("Track_Legs_%s",die->GetMCSignals()->At(i)->GetName()));
      // histos->AddClass(Form("Pair_%s_MCtruth",die->GetMCSignals()->At(i)->GetName()));
    }
  }
  
  //add histograms to event class
  histos->UserHistogram("Event","nEvents","Number of processed events after cuts;Number events",1,0,1,AliDielectronVarManager::kNevents);
  histos->UserHistogram("Event","ZVertex","ZVertex;ZVertex/cm",120,-12.,12.,AliDielectronVarManager::kZvPrim);
  histos->UserHistogram("Event","Centrality","Centrality;Centrality/%",202,-1.,100.,AliDielectronVarManager::kCentralityNew);
  histos->UserHistogram("Event","nEvTPC_eventplaneents",";;ev plane;",AliDielectronHelper::MakeLinBinning(180,  TMath::Pi()/-2.,TMath::Pi()/2.),AliDielectronVarManager::kQnTPCrpH2);


  //add histograms to track class
  histos->UserHistogram("Track","Pt","Pt;Pt [GeV];#tracks",500,0.,10.,AliDielectronVarManager::kPt);
  histos->UserHistogram("Track","P","P;P [GeV];#tracks",500,0.,10.,AliDielectronVarManager::kP);
  histos->UserHistogram("Track","PIn","PIn;PIn [GeV];#tracks",500,0.,10.,AliDielectronVarManager::kPIn);
  histos->UserHistogram("Track","Eta_phi","Eta vs Phi;Eta;Phi",90,-0.9,0.9,160,0.,6.4,AliDielectronVarManager::kEta,AliDielectronVarManager::kPhi);
  histos->UserHistogram("Track","Pt_phi","Pt vs Phi;Pt;Phi [GeV];#tracks",500,0.,5.,320,0.,6.4,AliDielectronVarManager::kPt,AliDielectronVarManager::kPhi);
  histos->UserHistogram("Track","ImpParXY","ImpParXY; ImpParXY ;#tracks",100,-5.,5.,AliDielectronVarManager::kImpactParXY);
  histos->UserHistogram("Track","ImpParZ","ImpParZ; ImpParZ ;#tracks",100,-5.,5.,AliDielectronVarManager::kImpactParZ);
  
  histos->UserHistogram("Track","NClusterTPC","NClusterTPC; NClusterTPC ;#tracks",200,-0.5,199.5,AliDielectronVarManager::kNclsTPC);
  histos->UserHistogram("Track","NClusterTPCShared","NClusterTPCShared; NClusterTPCShared ;#tracks",200,-0.5,199.5,AliDielectronVarManager::kNclsSTPC);
  histos->UserHistogram("Track","CrossedRows","CrossedRows; CrossedRows ;#tracks",200,-0.5,199.5,AliDielectronVarManager::kNFclsTPCr);
  histos->UserHistogram("Track","CrossedRowsOverFindable","CrRowsOverFindable; CrRows/FindableCls ;#tracks",120,0.,1.2,AliDielectronVarManager::kNFclsTPCfCross);
  histos->UserHistogram("Track","TPCchi2perCls","TPCchi2perCls; TPCchi2perCls ;#tracks",100,0.,10.,AliDielectronVarManager::kTPCchi2Cl);
  
  histos->UserHistogram("Track","NClusterITS","NClusterITS; NClusterITS ;#tracks",8,-0.5,7.5,AliDielectronVarManager::kNclsITS);
  histos->UserHistogram("Track","ITSchi2perCls","ITSchi2perCls; ITSchi2perCls ;#tracks",100,0.,10.,AliDielectronVarManager::kITSchi2Cl);

  histos->UserHistogram("Track","TPCdEdx_P","dEdx;P [GeV];TPC signal (arb units) vs Momentum;Mom;TPCsignal",     200,0.,10.,150,  0.,150. ,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCsignal);  
  histos->UserHistogram("Track","TPCnSigma_MomEle","TPC number of sigmas Electrons vs Momentum;Mom;TPCsigmaEle", 200,0.,10.,300,-30., 30. ,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle);
  histos->UserHistogram("Track","ITSdEdx_P","dEdx;P [GeV];ITS signal (arb units) vs Momentum;Mom;ITSsignal",     200,0.,10.,150,  0.,150. ,AliDielectronVarManager::kPIn,AliDielectronVarManager::kITSsignal);
  histos->UserHistogram("Track","ITSnSigma_MomEle","ITS number of sigmas Electrons vs Momentum;Mom;ITSsigmaEle"                           ,     200,0.,10.,300,-30., 30. ,AliDielectronVarManager::kPIn,AliDielectronVarManager::kITSnSigmaEle);
  histos->UserHistogram("Track","TOFbeta_Mom","kTOFbeta vs Momentum;Mom;TOFbeta"                           ,     200,0.,10.,120,  0.,  1.2,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFbeta);
  histos->UserHistogram("Track","TOFnSigma_MomEle","TOF number of sigmas Electrons vs Momentum;Mom;TOFsigmaEle"                           ,     200,0.,10.,300,-30., 30. ,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaEle);

  //add histograms to pair classes
  histos->UserHistogram("Pair",
                        "InvMass_pPt","Inv.Mass:PairPt;Inv. Mass (GeV/c^{2});Pair Pt (GeV/c)",
                        500,0.,5.,250,0.,5.,
                        AliDielectronVarManager::kM, AliDielectronVarManager::kPt);

  histos->UserHistogram("Pair","InvMass_PairPt_PhivPair","InvMass:PairPt:PhivPair;Inv. Mass [GeV];Pair Pt [GeV];PhiV",
                        600,0.,6., 600,0.,6., 20,0.,TMath::Pi(),
                        AliDielectronVarManager::kM, AliDielectronVarManager::kPt, AliDielectronVarManager::kPhivPair);  
  						
  histos->UserHistogram("Pair",
                        "Eta_phi_pair","Eta vs Phi (pair);Eta;Phi",
                        50,-1.,1.,80,0.,6.4,
                        AliDielectronVarManager::kEta,AliDielectronVarManager::kPhi);
			
  histos->UserHistogram("Pair",
                        "InvMass_PhivPair","InvMass_PhivPair;InvMass;PhivPair",
                         50, 0. , 0.5, 160 , 0., 3.2,
                         AliDielectronVarManager::kM , AliDielectronVarManager::kPhivPair );

  histos->UserHistogram("Pair",
		            	"InvMass_OpAngle","InvMass_OpAngle;Invariant Mass;Opening angle",
		            	100, 0., 0.5 ,160, 0., 3.2,
		            	AliDielectronVarManager::kM,AliDielectronVarManager::kOpeningAngle);

  histos->UserHistogram("Pair",
                        "Y","Y;counts;Y",
                        60, -1.2 , 1.2, 
                        AliDielectronVarManager::kY);

  if(cutDefinition > nPF){ 
    histos->UserHistogram("Pre","Pt",";Pt [GeV];#tracks",200,0,10.,AliDielectronVarManager::kPt);

    histos->UserHistogram("RejPair","Pt",";Pt [GeV];#tracks",200,0,10.,AliDielectronVarManager::kPt);

    histos->UserHistogram("RejPair",
                          "OpAngle_InvMass","InvMass_openingAngle;Invariant Mass;opening angle",
                          100, 0., 0.2, 100, 0. ,0.2,
                          AliDielectronVarManager::kM, AliDielectronVarManager::kOpeningAngle);
                        
    histos->UserHistogram("RejTrack","Pt",";Pt [GeV];#tracks",200,0,10.,AliDielectronVarManager::kPt); 
  
    histos->UserHistogram("Track_Legs","Pt",";Pt [GeV];#tracks",200,0,10.,AliDielectronVarManager::kPt); 
  }
  die->SetHistogramManager(histos);

}

const AliDielectronEventCuts *GetEventCuts(){

  AliDielectronEventCuts *eventCuts=new AliDielectronEventCuts("eventCuts","Vertex SPD && |vtxZ|<10 && ncontrib>0");
  eventCuts->SetRequireVertex();
  eventCuts->SetVertexType(AliDielectronEventCuts::kVtxAny); // AOD
  eventCuts->SetVertexZ(-10.,10.);
  eventCuts->SetMinVtxContributors(1); 

  //no centrality cuts for the moment
  //Bool_t isRun2 = kTRUE;
  //eventCuts->SetCentralityRange(0,80,isRun2);

  return eventCuts;
}


void SetupMCsignals(AliDielectron* die){

  Printf("Setting up MC signals...");

  // ##################### "real" pairs from signals (pi0,eta,eta',rho, omega, phi) ##############################
   AliDielectronSignalMC* pi0Sig = new AliDielectronSignalMC("pi0", "pi0Signal"); ///pi0 dalitz pairs 
  pi0Sig->SetLegPDGs(11,-11);
  pi0Sig->SetMotherPDGs(111,111);
  pi0Sig->SetMothersRelation(AliDielectronSignalMC::kSame);
  pi0Sig->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  pi0Sig->SetMotherSources(AliDielectronSignalMC::kPrimary, AliDielectronSignalMC::kPrimary);
  pi0Sig->SetCheckBothChargesLegs(kTRUE,kTRUE);
  pi0Sig->SetCheckBothChargesMothers(kTRUE,kTRUE);
  pi0Sig->SetFillPureMCStep(kTRUE);
  die->AddSignalMC(pi0Sig);

  AliDielectronSignalMC* pi0All = new AliDielectronSignalMC("pi0", "pi0All"); ///pi0 dalitz pairs (also from secondary)
  pi0All->SetLegPDGs(11,-11);
  pi0All->SetMotherPDGs(111,111);
  pi0All->SetMothersRelation(AliDielectronSignalMC::kSame);
  pi0All->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  pi0All->SetCheckBothChargesLegs(kTRUE,kTRUE);
  pi0All->SetCheckBothChargesMothers(kTRUE,kTRUE);
  pi0All->SetFillPureMCStep(kTRUE);
  die->AddSignalMC(pi0All);


  AliDielectronSignalMC* etaSig = new AliDielectronSignalMC("Eta", "etaSignal"); ///eta dalitz pairs 
  etaSig->SetLegPDGs(11,-11);
  etaSig->SetMotherPDGs(221,221);
  etaSig->SetMothersRelation(AliDielectronSignalMC::kSame);
  etaSig->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  etaSig->SetMotherSources(AliDielectronSignalMC::kPrimary, AliDielectronSignalMC::kPrimary);
  etaSig->SetCheckBothChargesLegs(kTRUE,kTRUE);
  etaSig->SetCheckBothChargesMothers(kTRUE,kTRUE);
  etaSig->SetFillPureMCStep(kTRUE);
  die->AddSignalMC(etaSig);


  AliDielectronSignalMC* etaprimeSig = new AliDielectronSignalMC("Etaprime", "etaprimeSignal"); ///etaprime pairs 
  etaprimeSig->SetLegPDGs(11,-11);
  etaprimeSig->SetMotherPDGs(331,331);
  etaprimeSig->SetMothersRelation(AliDielectronSignalMC::kSame);
  etaprimeSig->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  etaprimeSig->SetMotherSources(AliDielectronSignalMC::kPrimary, AliDielectronSignalMC::kPrimary);
  etaprimeSig->SetCheckBothChargesLegs(kTRUE,kTRUE);
  etaprimeSig->SetCheckBothChargesMothers(kTRUE,kTRUE);
  etaprimeSig->SetFillPureMCStep(kTRUE);
  die->AddSignalMC(etaprimeSig);


  AliDielectronSignalMC* rhoSig = new AliDielectronSignalMC("Rho", "rhoSignal"); ///rho pairs 
  rhoSig->SetLegPDGs(11,-11);
  rhoSig->SetMotherPDGs(113,113);
  rhoSig->SetMothersRelation(AliDielectronSignalMC::kSame);
  rhoSig->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  rhoSig->SetMotherSources(AliDielectronSignalMC::kPrimary, AliDielectronSignalMC::kPrimary);
  rhoSig->SetCheckBothChargesLegs(kTRUE,kTRUE);
  rhoSig->SetCheckBothChargesMothers(kTRUE,kTRUE);
  rhoSig->SetFillPureMCStep(kTRUE);
  die->AddSignalMC(rhoSig);

  AliDielectronSignalMC* omegaSig = new AliDielectronSignalMC("Omega", "omegaSignal"); ///omega pairs 
  omegaSig->SetLegPDGs(11,-11);
  omegaSig->SetMotherPDGs(223,223);
  omegaSig->SetMothersRelation(AliDielectronSignalMC::kSame);
  omegaSig->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  omegaSig->SetMotherSources(AliDielectronSignalMC::kPrimary, AliDielectronSignalMC::kPrimary);
  omegaSig->SetCheckBothChargesLegs(kTRUE,kTRUE);
  omegaSig->SetCheckBothChargesMothers(kTRUE,kTRUE);
  omegaSig->SetFillPureMCStep(kTRUE);
  die->AddSignalMC(omegaSig);
  
  AliDielectronSignalMC* phiSig = new AliDielectronSignalMC("Phi", "phiSignal"); ///phi pairs 
  phiSig->SetLegPDGs(11,-11);
  phiSig->SetMotherPDGs(333,333);
  phiSig->SetMothersRelation(AliDielectronSignalMC::kSame);
  phiSig->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  phiSig->SetMotherSources(AliDielectronSignalMC::kPrimary, AliDielectronSignalMC::kPrimary);
  phiSig->SetCheckBothChargesLegs(kTRUE,kTRUE);
  phiSig->SetCheckBothChargesMothers(kTRUE,kTRUE);
  phiSig->SetFillPureMCStep(kTRUE);
  die->AddSignalMC(phiSig);
  
  // ##################### "real" pairs from photon conversions in the detector material ##############################
  AliDielectronSignalMC* signalFromResonance_ULS_gammaConv = new AliDielectronSignalMC("signalFromResonance_ULS_gammaConv", "signalFromResonance_ULS_gammaConv");
  signalFromResonance_ULS_gammaConv->SetLegPDGs(11,-11);
  signalFromResonance_ULS_gammaConv->SetMotherPDGs(22,22);
  signalFromResonance_ULS_gammaConv->SetMothersRelation(AliDielectronSignalMC::kSame);
  signalFromResonance_ULS_gammaConv->SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary); // kSecondary means decays in the detector
  signalFromResonance_ULS_gammaConv->SetCheckBothChargesLegs(kTRUE,kTRUE);
  die->AddSignalMC(signalFromResonance_ULS_gammaConv);

  // ##################### combinatorial pairs ##############################
  AliDielectronSignalMC* diEleCombinatiorial = new AliDielectronSignalMC("diEleCombinatiorial", "diEleCombinatiorial");
  diEleCombinatiorial->SetLegPDGs(11,-11);
  diEleCombinatiorial->SetMothersRelation(AliDielectronSignalMC::kDifferent);
  diEleCombinatiorial->SetCheckBothChargesLegs(kTRUE,kTRUE);
  die->AddSignalMC(diEleCombinatiorial);

  AliDielectronSignalMC* diEleCombinatiorialConversion = new AliDielectronSignalMC("diEleCombinatiorialConversion", "diEleCombinatiorialConversion");
  diEleCombinatiorialConversion->SetLegPDGs(11,-11);
  diEleCombinatiorialConversion->SetMotherPDGs(22,0);// 1 leg from photons + 1 leg from everything
  diEleCombinatiorialConversion->SetMothersRelation(AliDielectronSignalMC::kDifferent);
  diEleCombinatiorialConversion->SetCheckBothChargesLegs(kTRUE,kTRUE);
  die->AddSignalMC(diEleCombinatiorialConversion);

    // ##################### HF pairs ##############################
  AliDielectronSignalMC* diEleHF = new AliDielectronSignalMC("diEleHF", "diEleHF");
  diEleHF->SetLegPDGs(11,-11);
  diEleHF->SetMotherPDGs(401,401);
  diEleHF->SetMothersRelation(AliDielectronSignalMC::kDifferent);
  diEleHF->SetCheckBothChargesLegs(kTRUE,kTRUE);
  die->AddSignalMC(diEleHF);

 
  // #################### D-Mesons
  AliDielectronSignalMC* diEleOpenCharmCharged = new AliDielectronSignalMC("DmesonsCharged","di-electrons from open charm D+- mesons no B grandmother");  // dielectrons originating from open charm hadrons
  diEleOpenCharmCharged->SetLegPDGs(11,-11);
  diEleOpenCharmCharged->SetMotherPDGs(401,401);
  diEleOpenCharmCharged->SetLegSources(AliDielectronSignalMC::kDontCare, AliDielectronSignalMC::kDontCare);
  diEleOpenCharmCharged->SetMothersRelation(AliDielectronSignalMC::kDifferent);
  // diEleOpenCharmCharged->SetFillPureMCStep(kTRUE);
  diEleOpenCharmCharged->SetCheckStackForPDG(kTRUE);
  diEleOpenCharmCharged->SetPDGforStack(503);
  diEleOpenCharmCharged->SetCheckBothChargesLegs(kTRUE,kTRUE);
  diEleOpenCharmCharged->SetCheckBothChargesMothers(kTRUE,kTRUE);
  die->AddSignalMC(diEleOpenCharmCharged);

  AliDielectronSignalMC* diEleOpenCharmNeutral = new AliDielectronSignalMC("DmesonsNeutral","di-electrons from open charm D0 mesons no B grandmother");  // dielectrons originating from open charm hadrons
  diEleOpenCharmNeutral->SetLegPDGs(11,-11);
  diEleOpenCharmNeutral->SetMotherPDGs(405,405);
  diEleOpenCharmNeutral->SetLegSources(AliDielectronSignalMC::kDontCare, AliDielectronSignalMC::kDontCare);
  diEleOpenCharmNeutral->SetMothersRelation(AliDielectronSignalMC::kDifferent);
  diEleOpenCharmNeutral->SetCheckStackForPDG(kTRUE);
  diEleOpenCharmNeutral->SetPDGforStack(503);
  diEleOpenCharmNeutral->SetCheckBothChargesLegs(kTRUE,kTRUE);
  diEleOpenCharmNeutral->SetCheckBothChargesMothers(kTRUE,kTRUE);
  die->AddSignalMC(diEleOpenCharmNeutral);

  //B meson (3)
  AliDielectronSignalMC* diEleOneOpenB = new AliDielectronSignalMC("B2ee","di-electrons from one B meson");  // dielectrons originating from open charm hadrons
  diEleOneOpenB->SetLegPDGs(11,-11);
  diEleOneOpenB->SetMotherPDGs(401,501);
  diEleOneOpenB->SetLegSources(AliDielectronSignalMC::kDontCare, AliDielectronSignalMC::kDontCare);
  diEleOneOpenB->SetGrandMotherPDGs(501,0);
  diEleOneOpenB->SetCheckMotherGrandmotherRelation(kTRUE,kTRUE);
  diEleOneOpenB->SetCheckBothChargesLegs(kTRUE,kTRUE);
  diEleOneOpenB->SetCheckBothChargesMothers(kTRUE,kTRUE);
  diEleOneOpenB->SetCheckBothChargesGrandMothers(kTRUE,kTRUE);
  die->AddSignalMC(diEleOneOpenB);

  //B meson (1)(1)
  AliDielectronSignalMC* diEleOpenB = new AliDielectronSignalMC("BMesons","di-electrons from B mesons");  // dielectrons originating from open charm hadrons
  diEleOpenB->SetLegPDGs(11,-11);
  diEleOpenB->SetMotherPDGs(501,501);
  diEleOpenB->SetLegSources(AliDielectronSignalMC::kDontCare, AliDielectronSignalMC::kDontCare);
  diEleOpenB->SetMothersRelation(AliDielectronSignalMC::kDifferent);
  diEleOpenB->SetCheckBothChargesLegs(kTRUE,kTRUE);
  diEleOpenB->SetCheckBothChargesMothers(kTRUE,kTRUE);
  die->AddSignalMC(diEleOpenB);

  //B meson (2)(2)
  AliDielectronSignalMC* diEleOpenBtoD = new AliDielectronSignalMC("B2D2ee","di-electrons from B->D-> e");  // dielectrons originating from open charm hadrons
  diEleOpenBtoD->SetLegPDGs(11,-11);
  diEleOpenBtoD->SetMotherPDGs(401,401);
  diEleOpenBtoD->SetLegSources(AliDielectronSignalMC::kDontCare, AliDielectronSignalMC::kDontCare);
  diEleOpenBtoD->SetGrandMotherPDGs(501,501);
  diEleOpenBtoD->SetGrandMothersRelation(AliDielectronSignalMC::kDifferent);
  diEleOpenBtoD->SetCheckBothChargesLegs(kTRUE,kTRUE);
  diEleOpenBtoD->SetCheckBothChargesMothers(kTRUE,kTRUE);
  diEleOpenBtoD->SetCheckBothChargesGrandMothers(kTRUE,kTRUE);
  die->AddSignalMC(diEleOpenBtoD);

  //B meson (1)(2)
  AliDielectronSignalMC* diEleOpenBandBtoD = new AliDielectronSignalMC("B2eAndB2D2e","di-electrons from B->e and B->D->e");  // dielectrons originating from open charm hadrons
  diEleOpenBandBtoD->SetLegPDGs        (11,11);
  diEleOpenBandBtoD->SetMotherPDGs     (401,501);
  diEleOpenBandBtoD->SetGrandMotherPDGs(501,0);
  diEleOpenBandBtoD->SetLegSources(AliDielectronSignalMC::kDontCare, AliDielectronSignalMC::kDontCare);
  diEleOpenBandBtoD->SetCheckBothChargesLegs(kTRUE,kTRUE);
  diEleOpenBandBtoD->SetCheckBothChargesMothers(kTRUE,kTRUE);
  diEleOpenBandBtoD->SetCheckBothChargesGrandMothers(kTRUE,kTRUE);
  //do i need this?
  diEleOpenBandBtoD->SetCheckMotherGrandmotherRelation(kTRUE,kFALSE);
  die->AddSignalMC(diEleOpenBandBtoD);
}

