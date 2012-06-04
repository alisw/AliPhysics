void InitHistograms(AliDielectron *die, Int_t cutDefinition);
void InitCF(AliDielectron* die, Int_t cutDefinition);

void SetupTrackCuts(AliDielectron *die, Int_t cutDefinition);
void SetupElectronCuts(AliDielectron *die, Int_t cutDefinition);

TVectorD *GetRunNumbers();
TVectorD *BinsToVector(Int_t nbins, Double_t min, Double_t max);

TString names=("CutStats;noPID;TPC;TOF;TRD;TOFTRD;TRDeff;e_mcPid;jpsi_mcPid;conv_mcPid");
enum {kCutStats=0, knoPID, kTPC, kTOF, kTRD, kTOFTRD, kTRDeff, kEleMC, kEleJPsiMC, kEleConvMC};
TObjArray *arrNames=names.Tokenize(";");

const Int_t nDie=arrNames->GetEntries();

Bool_t hasMC=kFALSE;

AliDielectron* ConfigJpsiQA_jb_PbPb(Int_t cutDefinition, Bool_t isMC=kFALSE)
{
  //
  // Setup the instance of AliDielectron
  //
  
  // MC event handler?
  hasMC=isMC;
  //hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);    

  //ESD handler?
  Bool_t isESD=(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()->IsA()==AliESDInputHandler::Class());
  
  // switch off some configurations
  switch(cutDefinition) {
    case kCutStats:   
      //case knoPID:
    case kTPC:
    case kTOF:
    case kTOFTRD:
    case kTRDeff:
      if(hasMC) return 0x0;
      break;
    case kTRD: return 0x0; break;
    case kEleMC: return 0x0; break;
    case kEleJPsiMC:
    case kEleConvMC:
      if(!hasMC) return 0x0;
      break;
  }
  
  // create the actual framework object
  TString name=Form("%02d",cutDefinition);
  if (cutDefinition<arrNames->GetEntriesFast()){
    name=arrNames->At(cutDefinition)->GetName();
  }
  AliDielectron *die = new AliDielectron(Form("%s",name.Data()),
                                         Form("Track cuts: %s",name.Data()));

  //only track QA, no Pairing
  if(cutDefinition<kTRDeff) die->SetNoPairing();
  
  // cut setup
  if(cutDefinition!=kTRDeff) SetupTrackCuts(die,cutDefinition);
  if(cutDefinition==kTRDeff) SetupElectronCuts(die,cutDefinition);
  
  
  
  // histogram setup
  InitHistograms(die,cutDefinition);

  return die;
}
//______________________________________________________________________________________
void SetupElectronCuts(AliDielectron *die, Int_t cutDefinition)
{
  // Quality cuts
  AliDielectronCutGroup* cuts = new AliDielectronCutGroup("cuts","cuts",AliDielectronCutGroup::kCompAND);
  die->GetTrackFilter().AddCuts(cuts);
  
  // track cuts ESD and AOD
  AliDielectronVarCuts *varCuts = new AliDielectronVarCuts("VarCuts","VarCuts");
  varCuts->AddCut(AliDielectronVarManager::kEta,         -0.9,   0.9);
  varCuts->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
  varCuts->AddCut(AliDielectronVarManager::kNclsTPC,     60.0, 160.0);
  varCuts->AddCut(AliDielectronVarManager::kKinkIndex0,   0.0);
  varCuts->AddCut(AliDielectronVarManager::kTRDpidQuality,1.0,   6.0);
  varCuts->AddCut(AliDielectronVarManager::kPt,           0.8,  1e30);
  cuts->AddCut(varCuts);
  
  AliDielectronTrackCuts *trkCuts = new AliDielectronTrackCuts("TrkCuts","TrkCuts");
  trkCuts->SetRequireTPCRefit(kTRUE);
  cuts->AddCut(trkCuts);
  
  // PID
  AliDielectronPID *pid = new AliDielectronPID("PID","PID cut");
  pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-3.0, 3.0, 0.0, 0.0,  kFALSE, AliDielectronPID::kRequire); 
  cuts->AddCut(pid);
  
  // Pair inclusion
  AliDielectronVarCuts *gammaCuts = new AliDielectronVarCuts("GammaCuts","GammaCuts");
  gammaCuts->AddCut(AliDielectronVarManager::kOpeningAngle, 0.0,   0.035); // 0.1
//  gammaCuts->AddCut(AliDielectronVarManager::kLegDist,      0.0,   0.25);
  gammaCuts->AddCut(AliDielectronVarManager::kR,            3.0,   90.0);
//  gammaCuts->AddCut(AliDielectronVarManager::kPsiPair,      0.0,   0.05);
//  gammaCuts->AddCut(AliDielectronVarManager::kChi2NDF,      0.0,   10.0);
  gammaCuts->AddCut(AliDielectronVarManager::kM,            0.0,   0.05);
  die->GetPairFilter().AddCuts(gammaCuts);

  // Pair exclusion
  //  AliDielectronVarCuts *lambdaCuts = new AliDielectronVarCuts("LambdaCuts","LambdaCuts");
  //  lambdaCuts->AddCut(AliDielectronVarManager::kOpeningAngle, 0.0,   0.035, kTRUE); // 0.035
  //  lambdaCuts->AddCut(AliDielectronVarManager::kR,            3.0,   40.0, kTRUE);
  //  lambdaCuts->AddCut(AliDielectronVarManager::kChi2NDF,      0.0,   0.5,  kTRUE);
  //  lambdaCuts->AddCut(AliDielectronVarManager::kM,            1.01,  1.03, kTRUE);
  //  lambdaCuts->AddCut(AliDielectronVarManager::kY,            -0.9,  0.9,  kTRUE);
  //  die->GetPairFilter().AddCuts(lambdaCuts);
  //
  //  AliDielectronVarCuts *k0sCuts = new AliDielectronVarCuts("K0sCuts","K0sCuts");
  //  k0sCuts->AddCut(AliDielectronVarManager::kR,            3.0,   40.0,  kTRUE);
  //  k0sCuts->AddCut(AliDielectronVarManager::kChi2NDF,      0.0,   0.5,   kTRUE);
  //  k0sCuts->AddCut(AliDielectronVarManager::kM,            4.90,  0.504, kTRUE);
  //  k0sCuts->AddCut(AliDielectronVarManager::kY,            -0.9,  0.9,   kTRUE);
  //  die->GetPairFilter().AddCuts(k0sCuts);
  
  
}
//______________________________________________________________________________________
void SetupTrackCuts(AliDielectron *die, Int_t cutDefinition)
{
  //
  // Setup the track cuts
  //
  //ESD handler?
  Bool_t isESD=(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()->IsA()==AliESDInputHandler::Class());

  // Quality cuts
  AliDielectronCutGroup* cuts = new AliDielectronCutGroup("cuts","cuts",AliDielectronCutGroup::kCompAND);
  die->GetTrackFilter().AddCuts(cuts);
  
  //Pt cut, should make execution a bit faster
  AliDielectronVarCuts *pt = new AliDielectronVarCuts("Pt>.8","Pt>.8");
  if(cutDefinition<kEleMC) 
    pt->AddCut(AliDielectronVarManager::kPt,1.1,1e30);
  else   
    pt->AddCut(AliDielectronVarManager::kPt,0.8,1e30);
  if(cutDefinition!=kCutStats && cutDefinition!=knoPID) cuts->AddCut(pt);
  
  // track cuts ESD and AOD
  AliDielectronVarCuts *varCuts = new AliDielectronVarCuts("VarCuts","VarCuts");
  varCuts->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
  varCuts->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
  varCuts->AddCut(AliDielectronVarManager::kEta,         -0.9,   0.9);
  varCuts->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
  varCuts->AddCut(AliDielectronVarManager::kNclsTPC,     70.0, 160.0);
  varCuts->AddCut(AliDielectronVarManager::kKinkIndex0,   0.0);
  if(cutDefinition!=kCutStats) cuts->AddCut(varCuts);
  
  AliDielectronTrackCuts *trkCuts = new AliDielectronTrackCuts("TrkCuts","TrkCuts");
  trkCuts->SetClusterRequirementITS(AliDielectronTrackCuts::kSPD,AliDielectronTrackCuts::kAny);  
  trkCuts->SetRequireITSRefit(kTRUE);
  trkCuts->SetRequireTPCRefit(kTRUE);
  if(cutDefinition!=kCutStats) cuts->AddCut(trkCuts);

  /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv PID CUTS vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
  AliDielectronPID *pid = new AliDielectronPID("PID","PID");
  
  ////////////////////////////////// DATA
  if(!hasMC) {
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,-100.,3.5,0.,0.,kTRUE);
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kProton,-100.,3.5,0.,0.,kTRUE);
  }
  
  ////////////////////////////////// MC
  if(hasMC) {
    if (cutDefinition>=kEleMC) {
      AliDielectronVarCuts *pdgCuts=new AliDielectronVarCuts("pdgCuts","pdgCuts");
      pdgCuts->SetCutType(AliDielectronVarCuts::kAny);
      pdgCuts->AddCut(AliDielectronVarManager::kPdgCode,11.);
      pdgCuts->AddCut(AliDielectronVarManager::kPdgCode,-11.);
      cuts->AddCut(pdgCuts);
      
      AliDielectronVarCuts *pdgMotherCutsT=new AliDielectronVarCuts("pdgMotherCutsT","pdgMotherCutsT");
      AliDielectronVarCuts *pdgMotherCutsP=new AliDielectronVarCuts("pdgMotherCutsP","pdgMotherCutsP");
      if (cutDefinition==kEleJPsiMC){
        //        pdgMotherCutsT->AddCut(AliDielectronVarManager::kPdgCodeMother,443.);
        //        cuts->AddCut(pdgMotherCutsT);        
        pdgMotherCutsP->AddCut(AliDielectronVarManager::kPdgCode,443.);
        pdgMotherCutsP->AddCut(AliDielectronVarManager::kHaveSameMother,1.);
        die->GetPairFilter().AddCuts(pdgMotherCutsP);
      }
      if (cutDefinition==kEleConvMC){
        pdgMotherCutsT->AddCut(AliDielectronVarManager::kPdgCodeMother,22.);
        cuts->AddCut(pdgMotherCutsT);        
        
        pdgMotherCutsP->AddCut(AliDielectronVarManager::kHaveSameMother,1.);
        die->GetPairFilter().AddCuts(pdgMotherCutsP);
      }
      
    }
  }
	////////////////////////////////// DATA + MC
  // pid cuts TPC + TOF & TRD
  pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-3.,3.);
  if(cutDefinition==kTOF || cutDefinition==kTOFTRD) 
    pid->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,-3,3.,0.,0.,kFALSE,AliDielectronPID::kIfAvailable);
  if(cutDefinition==kTRD || cutDefinition==kTOFTRD) 
    pid->AddCut(AliDielectronPID::kTRDeleEff,AliPID::kElectron,.8,1.,3.5.,6.,kFALSE,
                AliDielectronPID::kIfAvailable,AliDielectronVarManager::kTRDpidQuality);

    if(cutDefinition!=knoPID    && 
       cutDefinition!=kCutStats && 
       cutDefinition!=kTRDeff   && 
       !hasMC ) cuts->AddCut(pid);
  /* ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ PID CUTS ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/

  // exclude conversion electrons selected by the tender
  //AliDielectronTrackCuts *noconv=new AliDielectronTrackCuts("noConv","noConv");
  //noconv->SetV0DaughterCut(AliPID::kElectron,kTRUE);
  //cuts->AddCut(noconv);
  
}

//______________________________________________________________________________________
void InitHistograms(AliDielectron *die, Int_t cutDefinition)
{
  //
  // Initialise the histograms
  //
  
  //Setup histogram Manager
  AliDielectronHistos *histos=new AliDielectronHistos(die->GetName(), die->GetTitle());
  
  //Initialise histogram classes
  histos->SetReservedWords("Track;Pair");

  // booleans for histo selection
  Bool_t bHistEvts = kFALSE, bHistPair = kFALSE, bHistCuts = kFALSE, bHistPID = kFALSE, bHistEff=kFALSE, bHistRunQA=kFALSE;
  switch (cutDefinition) {
    case kCutStats:  bHistCuts=kTRUE; break;
    case knoPID:     bHistEvts=kTRUE; bHistPID=kTRUE; bHistRunQA=kTRUE; break;
    case kTPC:       
    case kTOF:       
    case kTRD:       
    case kTOFTRD: 
      bHistPID=kTRUE; break;
    case kTRDeff:    bHistEff=kTRUE; break;
    case kEleMC:
    case kEleJPsiMC: 
    case kEleConvMC: 
      bHistPair=kTRUE; break;
  }
  
  
  //Track classes
  //to fill also track info from 2nd event loop until 2
  for (Int_t i=0; i<2; ++i){
    histos->AddClass(Form("Track_%s",AliDielectron::TrackClassName(i)));
  }
  
  //add histograms to event class
  if (bHistEvts) {
    histos->AddClass("Event");
    histos->UserHistogram("Event","RunNumber",";run;#events",GetRunNumbers(),AliDielectronVarManager::kRunNumber);
    histos->UserHistogram("Event","VtxZ",";z_{vtx} (cm)",300,-15.,15.,AliDielectronVarManager::kZvPrim);
    if(hasMC) {
      histos->AddClass("MCEvent");
      histos->UserHistogram("MCEvent","NumberOfJPsis",";N_{J/#psi};#events",20,0.,20.,AliDielectronVarManager::kNumberOfJPsis);
    }
  }
  
  if (bHistRunQA && !hasMC) {
    // Flow QA
    histos->UserHistogram("Event","TPCrpH2uc_RunNumber",";run;#Psi^{TPC} (rad.)",
                          GetRunNumbers(), BinsToVector(100,-2.,2.) ,
                          AliDielectronVarManager::kRunNumber, AliDielectronVarManager::kTPCrpH2uc);
    histos->UserHistogram("Event","vOArpH2_RunNumber",";run;#Psi_{2}^{V0A} (rad.)",
                          GetRunNumbers(), BinsToVector(100,-2.,2.) ,
                          AliDielectronVarManager::kRunNumber, AliDielectronVarManager::kv0ArpH2);
    histos->UserHistogram("Event","vOCrpH2_RunNumber",";run;#Psi_{2}^{V0C} (rad.)",
                          GetRunNumbers(), BinsToVector(100,-2.,2.) ,
                          AliDielectronVarManager::kRunNumber, AliDielectronVarManager::kv0CrpH2);
    
    histos->UserHistogram("Event","TPCrpH2uc_Cent_RunNumber",";centrality (%);#Psi^{TPC} (rad.);run",
                          BinsToVector(10,0.,100.), BinsToVector(100,-2.,2.), GetRunNumbers(),
                          AliDielectronVarManager::kCentrality,
                          AliDielectronVarManager::kTPCrpH2uc,
                          AliDielectronVarManager::kRunNumber );
    histos->UserHistogram("Event","vOArpH2_Cent_RunNumber",";centrality (%);#Psi_{2}^{V0A} (rad.);run",
                          BinsToVector(10,0.,100.), BinsToVector(100,-2.,2.), GetRunNumbers(),
                          AliDielectronVarManager::kCentrality,
                          AliDielectronVarManager::kv0ArpH2,
                          AliDielectronVarManager::kRunNumber );
    histos->UserHistogram("Event","vOCrpH2_Cent_RunNumber",";centrality (%);#Psi_{2}^{V0C} (rad.);run",
                          BinsToVector(10,0.,100.), BinsToVector(100,-2.,2.), GetRunNumbers(),
                          AliDielectronVarManager::kCentrality,
                          AliDielectronVarManager::kv0CrpH2,
                          AliDielectronVarManager::kRunNumber );                      

    // PID QA
    histos->UserHistogram("Track","TPCnSigmaPio_Nacc_RunNumber",";N_{acc};n#sigma_{pio}^{TPC};run",
                          BinsToVector(60,0.,3000.), BinsToVector(40,-5.,5.), GetRunNumbers(),
                          AliDielectronVarManager::kNacc,
                          AliDielectronVarManager::kTPCnSigmaPio,
                          AliDielectronVarManager::kRunNumber);
    histos->UserHistogram("Track","TPCnSigmaPio_Nacc",";N_{acc};n#sigma_{pio}^{TPC}",
                          BinsToVector(60,0.,3000.), BinsToVector(40,-5.,5.),
                          AliDielectronVarManager::kNacc,
                          AliDielectronVarManager::kTPCnSigmaPio);
    
    histos->UserProfile("Track","TPCnSigmaPio-Nacc-RunNumber",";N_{acc};run;n#sigma_{pio}^{TPC}",
                        AliDielectronVarManager::kTPCnSigmaPio,
                        BinsToVector(60,0.,3000.),      GetRunNumbers(),
                        AliDielectronVarManager::kNacc, AliDielectronVarManager::kRunNumber, "s;-5;5" );
    histos->UserProfile("Track","TPCnSigmaPio-RunNumber",";run;n#sigma_{pio}^{TPC}",
                        AliDielectronVarManager::kTPCnSigmaPio,
                        GetRunNumbers(), AliDielectronVarManager::kRunNumber, "s;-5;5");
    histos->UserProfile("Track","TPCnSigmaPio-Nacc",";N_{acc};n#sigma_{pio}^{TPC}",
                        AliDielectronVarManager::kTPCnSigmaPio,
                        BinsToVector(60,0.,3000.), AliDielectronVarManager::kNacc , "s;-5;5");
    
  }
  
  if (bHistPID) {	
    
    if(cutDefinition==kTPC) {     // centrality dependence
      histos->UserHistogram("Track","dEdx_P_Cent",";p (GeV/c);TPC signal (arb units);centrality (%)",
                            400,0.2,20., 200,0.,200., 10,0.,100.,
                            AliDielectronVarManager::kPIn,
                            AliDielectronVarManager::kTPCsignal,
                            AliDielectronVarManager::kCentrality, 
                            kTRUE);
      
      histos->UserHistogram("Track","TPCnSigmaEle_P_Cent",";p (GeV/c);n#sigma_{ele}^{TPC};centrality (%)",
                            100,0.2,20., 100,-10.,10., 10,0.,100.,
                            AliDielectronVarManager::kPIn,
                            AliDielectronVarManager::kTPCnSigmaEle,
                            AliDielectronVarManager::kCentrality, 
                            kTRUE);
      
      histos->UserHistogram("Track","TPCnSigmaPio_P_Cent",";p (GeV/c);n#sigma_{pio}^{TPC};centrality (%)",
                            100,0.2,20., 100,-10.,10., 10,0.,100.,
                            AliDielectronVarManager::kPIn,
                            AliDielectronVarManager::kTPCnSigmaPio,
                            AliDielectronVarManager::kCentrality, 
                            kTRUE);
      
      histos->UserHistogram("Track","TPCnSigmaPro_P_Cent",";p (GeV/c);n#sigma_{pro}^{TPC};centrality (%)",
                            100,0.2,20., 100,-10.,10., 10,0.,100.,
                            AliDielectronVarManager::kPIn,
                            AliDielectronVarManager::kTPCnSigmaPro,
                            AliDielectronVarManager::kCentrality, 
                            kTRUE);
      
      histos->UserHistogram("Track","TOFbeta_P_Cent",";p (GeV/c);#beta;centrality (%);#tracks",
                            250,0.0,5., 300,0.,1.2, 10,0.,100.,
                            AliDielectronVarManager::kPIn,
                            AliDielectronVarManager::kTOFbeta,
                            AliDielectronVarManager::kCentrality);
      
      histos->UserHistogram("Track","TOFnSigmaEle_P_Cent","dEdxTOF;p (GeV/c);n#sigma_{ele}^{TOF};centrality (%);#tracks",
                            250,0.0,5., 100,-10.,10., 10,0.,100.,
                            AliDielectronVarManager::kPIn,
                            AliDielectronVarManager::kTOFnSigmaEle,
                            AliDielectronVarManager::kCentrality);
    }
    
    // main pid spectra
    histos->UserHistogram("Track","dEdx_P",";p (GeV/c);TPC signal (arb units);#tracks",
                          400,0.2,20.,200,0.,200.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCsignal,kTRUE);
    histos->UserHistogram("Track","TPCnSigmaEle_P",";p (GeV/c);n#sigma_{ele}^{TPC}",
                          100,0.2,20.,100,-10.,10.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,kTRUE);
    histos->UserHistogram("Track","TPCnSigmaPio_P",";p (GeV/c);n#sigma_{pio}^{TPC}",
                        100,0.2,20.,100,-10.,10.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaPio,kTRUE);
    histos->UserHistogram("Track","TPCnSigmaPro_P",";p (GeV/c);n#sigma_{pro}^{TPC}",
                          100,0.2,20.,100,-10.,10.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaPro,kTRUE);
    
    histos->UserHistogram("Track","TOFbeta_P",";p (GeV/c);#beta;#tracks",
                          250,0.0,5.0,300,0.,1.2,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFbeta);
    histos->UserHistogram("Track","TOFnSigmaEle_P",";p (GeV/c);n#sigma_{ele}^{TOF}",
                          100,0.2,20.,50,-10.,10.,AliDielectronVarManager::kP,AliDielectronVarManager::kTOFnSigmaEle,kTRUE);
  }
  
  if (bHistPair) {
    //Pair classes
    for (Int_t i=1; i<=1; ++i){ // only +- pairs
      histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(i)));
    }
    //add histograms to Pair classes
    histos->UserHistogram("Pair","InvMass",";m_{ee} (GeV/c^{2});#pairs",
                          100,.0,100*0.005, AliDielectronVarManager::kM); // 5MeV bins, 0.5GeV/c2
    histos->UserHistogram("Pair","OpeningAngle",";angle (rad.)",
                          100,0.,3.15, AliDielectronVarManager::kOpeningAngle);
    histos->UserHistogram("Pair","Chi2NDF",";#chi^{2}/NDF",
                          100,0.,20, AliDielectronVarManager::kChi2NDF);
    histos->UserHistogram("Pair","PsiPair","PsiPair;#psi",
                          100,0.,3.15, AliDielectronVarManager::kPsiPair);
    histos->UserHistogram("Pair","R","Radius;r (cm)",
                          200,0.,100., AliDielectronVarManager::kR);
    histos->UserHistogram("Pair","LegDist",";dca (cm)",
                          100,0.,1., AliDielectronVarManager::kLegDist);
    histos->UserHistogram("Pair","LegDistXY",";dca_{xy} (cm)",
                          100,0.,1., AliDielectronVarManager::kLegDistXY);
    histos->UserHistogram("Pair","PdgCode",";mother PDG code;#tracks",
                          10000,-5000.5,4999.5,AliDielectronVarManager::kPdgCode);
    
    // ITS tracks
    histos->UserHistogram("Track","NclsITS",";N_{cls}^{ITS};#tracks", 
                          "0,1,2,3,4,5,6,7",  AliDielectronVarManager::kNclsITS);
    histos->UserHistogram("Track","ITSLayerFirstCls",";ITS layer first cls;#tracks", 
                          "-1,0,1,2,3,4,5,6", AliDielectronVarManager::kITSLayerFirstCls);
  }

  // TRD efficiency
  if(bHistEff) {
    for (Int_t i=1; i<=1; ++i){ // only +- pairs
      //      histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(i)));
      //legs from pair
      histos->AddClass(Form("Track_Legs_%s",AliDielectron::PairClassName(i)));
    }
        
    // purity
    histos->UserHistogram("Track","TPCnSigmaEle",";n#sigma_{ele}^{TPC};#tracks",
                          100,-10.,10.,AliDielectronVarManager::kTPCnSigmaEle);
    // TRD
    histos->UserHistogram("Track","TRDpidQuality",";N_{trkl}^{TRD};#tracks", 
                          "0,1,2,3,4,5,6,7",    AliDielectronVarManager::kTRDpidQuality);
    histos->UserHistogram("Track","TRDprobEle",";P_{ele}^{TRD};#tracks", 
                          100,0.,1.,            AliDielectronVarManager::kTRDprobEle);
    histos->UserHistogram("Track","TRDprobPio",";P_{pio}^{TRD};#tracks", 
                          100,0.,1.,            AliDielectronVarManager::kTRDprobPio);

    histos->UserHistogram("Track","TRDprobEle_TRDpidQuality",";N_{trkl}^{TRD};P_{ele}^{TRD}", 
                          7,0.,7., 20,0.,1.,
                          AliDielectronVarManager::kTRDpidQuality, AliDielectronVarManager::kTRDprobEle);
    histos->UserHistogram("Track","TRDprobPio_TRDpidQuality",";N_{trkl}^{TRD};P_{pio}^{TRD}", 
                          7,0.,7., 20,0.,1.,
                          AliDielectronVarManager::kTRDpidQuality, AliDielectronVarManager::kTRDprobPio);
    
    histos->UserProfile("Track","TRDprobEle-TRDpidQuality",";N_{trkl}^{TRD};P_{ele}^{TRD}", 
                        AliDielectronVarManager::kTRDprobEle,
                        "0,1,2,3,4,5,6,7",    AliDielectronVarManager::kTRDpidQuality);
    histos->UserProfile("Track","TRDprobPio-TRDpidQuality",";N_{trkl}^{TRD};P_{pio}^{TRD}", 
                        AliDielectronVarManager::kTRDprobPio,
                        "0,1,2,3,4,5,6,7",    AliDielectronVarManager::kTRDpidQuality);
    
  }
  
  if(0) {    
    histos->UserHistogram("Track","Pt",";p_{T} (GeV/c);#tracks",200,0,20.,AliDielectronVarManager::kPt);
    histos->UserHistogram("Track","Eta",";#eta;#tracks",100,-2.,2.,AliDielectronVarManager::kEta);
    histos->UserHistogram("Track","Phi",";#phi;#tracks",360,0.,6.285,AliDielectronVarManager::kPhi);  
    if(hasMC)
      histos->UserHistogram("Track","PdgCodeMother",";mother PDG code;#tracks",10000,-5000.5,4999.5,AliDielectronVarManager::kPdgCodeMother);
    //    histos->UserHistogram("Track","PdgCode",";tracks PDG code;#tracks",10000,-5000.5,4999.5,AliDielectronVarManager::kPdgCode);
  }
  
  //add histograms to get cut statistics
  if (bHistCuts) {	
    histos->UserHistogram("Track","Eta",";Eta;#tracks",
                          "-5.,-0.9,-0.8,0.8,0.9,5.",          AliDielectronVarManager::kEta);
    histos->UserHistogram("Track","Pt",";p_{T} (GeV/c);#tracks", 
                          "0.0,0.8, 1.0, 1.2, 1.5, 100.0",    AliDielectronVarManager::kPt);
    histos->UserHistogram("Track","ImpactParXY",";dXY (cm);#tracks",  
                          500,-1.,1.,               AliDielectronVarManager::kImpactParXY);
    histos->UserHistogram("Track","ImpactParZ",";dZ (cm);#tracks",	  
                          600,-3.,3.,               AliDielectronVarManager::kImpactParZ);
    // ITS
    histos->UserHistogram("Track","NclsITS",";N_{cls}^{ITS};#tracks", 
                          "0,1,2,3,4,5,6,7",          AliDielectronVarManager::kNclsITS);
    histos->UserHistogram("Track","ITSLayerFirstCls",";ITS layer first cls;#tracks", 
                          "-1,0,1,2,3,4,5,6",          AliDielectronVarManager::kITSLayerFirstCls);
    // TPC
    histos->UserHistogram("Track","NclsTPC",";N_{cls}^{TPC};#tracks", 
                          "70, 90, 100, 120, 160",   AliDielectronVarManager::kNclsTPC);
    histos->UserHistogram("Track","TPCchi2PerCluster",";#chi^{2}/N_{cls}^{TPC};#tracks", 
                          10,0,10,AliDielectronVarManager::kTPCchi2Cl);
    histos->UserHistogram("Track","TPC_nSigma_Electrons",";n#sigma_{ele}^{TPC};#tracks", 
                          "-100,-4,-3,-2,-1.5,-1,1,1.5,2,3,4,100",AliDielectronVarManager::kTPCnSigmaEle);
    histos->UserHistogram("Track","TPC_nSigma_Pions",";n#sigma_{pio}^{TPC};#tracks", 
                          "-100,3.5,4.0,4.5,5.0,5.5,100",AliDielectronVarManager::kTPCnSigmaPio);
    histos->UserHistogram("Track","TPC_nSigma_Protons",";n#sigma_{pro}^{TPC};#tracks", 
                          "-100,3.5,4.0,4.5,5.0,5.5,100",AliDielectronVarManager::kTPCnSigmaPro);
    // TRD
    histos->UserHistogram("Track","TRDpidQuality",";N_{trkl}^{TRD};#tracks", 
                          "0,1,2,3,4,5,6,7",            AliDielectronVarManager::kTRDpidQuality);
    histos->UserHistogram("Track","TRDprobEle",";P_{ele}^{TRD};#tracks", 
                          100,0.,1.,            AliDielectronVarManager::kTRDprobEle);
    //    histos->UserHistogram("Track","TRD_PIDbit",";TRD pid bit;#tracks",
    //                          "-.5,.5,1.5",              AliDielectronVarManager::kTRDPIDBit);

    // TOF
    histos->UserHistogram("Track","TOF_PIDbit",";TOF pid bit;#tracks",
                          "-.5,.5,1.5",              AliDielectronVarManager::kTOFPIDBit);
    histos->UserHistogram("Track","TOF_nSigma_Electrons",";n#sigma_{ele}^{TOF};#tracks", 
                          "-100,-3,-2,2,3,100",AliDielectronVarManager::kTOFnSigmaEle);
  }
  
  
  // track histos
  //  histos->UserHistogram("Track","Pt",";p_{T} (GeV/c);#tracks",200,0,20.,AliDielectronVarManager::kPt);
  //  histos->UserHistogram("Track","Eta",";#eta;#tracks",100,-2.,2.,AliDielectronVarManager::kEta);
  //  histos->UserHistogram("Track","Phi",";#phi;#tracks",360,0.,6.285,AliDielectronVarManager::kPhi);  
  //  histos->UserHistogram("Track","Eta_Phi",";#eta; #phi;#tracks",
  //                        100,-2,2,144,0,6.285,AliDielectronVarManager::kEta,AliDielectronVarManager::kPhi);
  //  
  //  histos->UserHistogram("Track","Ncl",";Number clusters TPC;Number clusters TPC",
  //                        160,-0.5,159.5,AliDielectronVarManager::kNclsTPC);
  //  histos->UserHistogram("Track","NclFrFrac",";found/findable clusters (robust);#tracks",
  //                        160,0,1.1,AliDielectronVarManager::kNFclsTPCrFrac);
  //  histos->UserHistogram("Track","TPCsignalN","Number of Clusters TPC PID;#clusteres TPC PID;#tracks",
  //                        160,-0.5,159.5,AliDielectronVarManager::kTPCsignalN);
  //
  //  histos->UserHistogram("Track","Ncl_NclPid","Number clusters TPC vs. number of clusters PID;#clusters TPC; #clusters TPC PID",
  //                        160,-0.5,159.5,160,-0.5,159.5,AliDielectronVarManager::kNclsTPC,AliDielectronVarManager::kTPCsignalN);
  //  histos->UserHistogram("Track","TPCchi2Cl","TPC #chi^{2}/cluster;TPC #chi^{2}/cluster;#tracks",
  //                        100,0.,10.,AliDielectronVarManager::kTPCchi2Cl);
  
  die->SetHistogramManager(histos);
}

TVectorD *GetRunNumbers() {
  
  Double_t runLHC10h[] = { // all good runs based on RCT 29.Mai
    139510, 139507, 139505, 139503, 139465, 139438, 139437, 139360, 139329, 139328, 139314, 139310, 139309, 139173, 139107, 139105, 139038, 139037, 139036, 139029, 139028, 138872, 138871, 138870, 138837, 138732, 138730, 138666, 138662, 138653, 138652, 138638, 138624, 138621, 138583, 138582, 138579, 138578, 138534, 138469, 138442, 138439, 138438, 138396, 138364, 138275, 138225, 138201, 138197, 138192, 138190, 137848, 137844, 137752, 137751, 137724, 137722, 137718, 137704, 137693, 137692, 137691, 137686, 137685, 137639, 137638, 137608, 137595, 137549, 137546, 137544, 137541, 137539, 137531, 137530, 137443, 137441, 137440, 137439, 137434, 137432, 137431, 137430, 137366, 137243, 137236, 137235, 137232, 137231, 137230, 137162, 137161, 137135
  };
  
  Double_t runLHC11h[] = { // all good runs based on RCT 29.Mai
    170593, 170572, 170388, 170387, 170315, 170313, 170312, 170311, 170309, 170308, 170306, 170270, 170269, 170268, 170230, 170228, 170207, 170204, 170203, 170193, 170163, 170159, 170155, 170091, 170089, 170088, 170085, 170084, 170083, 170081, 170040, 170027, 169965, 169923, 169859, 169858, 169855, 169846, 169838, 169837, 169835, 169591, 169588, 169586, 169557, 169554, 169550, 169512, 169504, 169498, 169475, 169419, 169418, 169417, 169415, 169411, 169238, 169167, 169160, 169156, 169148, 169145, 169144, 169138, 169099, 169094, 169091, 169044, 169035, 168992, 168988, 168826, 168777, 168514, 168512, 168511, 168467, 168464, 168460, 168458, 168362, 168361, 168342, 168341, 168325, 168322, 168311, 168310, 168115, 168108, 168107, 168105, 168076, 168069, 167988, 167987, 167985, 167920, 167915
  };
  
  // selection via environement variable (works only for gsi trains)
  TString list=gSystem->Getenv("LIST");
  
  if(list.Contains("10h") || list.Contains("11a10b")) {
    Int_t size = (int) (sizeof(runLHC10h)/sizeof(Double_t));
    TVectorD *vec = new TVectorD(size+1);
    
    (*vec)[size] = runLHC10h[0] + 1;
    for (int i = 0; i < size; i++) {
      (*vec)[i] = runLHC10h[size-1-i];
    }
//    vec->Print("");
    return vec;
  }
  
  if( list.IsNull() || list.Contains("11h") || list.Contains("12a17") ) {
    
    Int_t size = (int) (sizeof(runLHC11h)/sizeof(Double_t));
    TVectorD *vec = new TVectorD(size+1);
    
    (*vec)[size] = runLHC11h[0] + 1;
    for (int i = 0; i < size; i++) {
      (*vec)[i] = runLHC11h[size-1-i];
    }
//    vec->Print("");
    return vec;
  }
  
}

TVectorD *BinsToVector(Int_t nbins, Double_t min, Double_t max) {
  return AliDielectronHelper::MakeLinBinning(nbins,min,max);
  //  TVectorD *vec = new TVectorD(nbins+1);
  //
  //  Double_t binwdth = (max-min)/nbins;
  //  for (int i = 0; i < nbins+1; i++) (*vec)[i] = min + i*binwdth;
  //  
  //  return vec;
}


