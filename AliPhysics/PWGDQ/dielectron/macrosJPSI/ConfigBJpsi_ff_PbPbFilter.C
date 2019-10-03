void InitHistograms(AliDielectron *die, Int_t cutDefinition);

void SetupTrackCuts(Bool_t isESD, AliDielectron *die, Int_t cutDefinition);
void SetupPairCuts(AliDielectron *die, Int_t cutDefinition);

void AddMCSignals(AliDielectron *die);
void SetEtaCorrection();

TString names=("TOFTRDany");
enum { kTOFTRD};

TObjArray *arrNames=names.Tokenize(";");
const Int_t nDie=arrNames->GetEntries();

Bool_t hasMC=kFALSE;

AliDielectron* ConfigBJpsi_ff_PbPbFilter(Int_t cutDefinition, Bool_t isMC=kFALSE)
{
  //
  // Setup the instance of AliDielectron
  //
  

  // MC event handler?
  hasMC=isMC;
    //(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);    

  //ESD handler?
  Bool_t isESD=(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()->IsA()==AliESDInputHandler::Class());
  
  
  // create the actual framework object
  TString name=Form("%02d",cutDefinition);
  if (cutDefinition<arrNames->GetEntriesFast()){
    name=arrNames->At(cutDefinition)->GetName();
  }
  AliDielectron *die = new AliDielectron(Form("%s",name.Data()),
                                         Form("Track cuts: %s",name.Data()));
  
  // Monte Carlo Signals and TRD efficiency tables
  if(hasMC) {
    AddMCSignals(die);
    
    // trd tables
    TString pidTab="$TRAIN_ROOT/util/dielectron/dielectron/TRDpidEff_eleProb07_TRDntr4_6.root";
    TString trainRoot=gSystem->Getenv("TRAIN_ROOT");
    if (trainRoot.IsNull()) pidTab="$ALICE_PHYSICS/PWGDQ/dielectron/files/TRDpidEff_eleProb07_TRDntr4_6.root";

    if (gSystem->AccessPathName(gSystem->ExpandPathName(pidTab.Data())))
      Error("ConfigPbPb","PID table not found: %s",pidTab.Data());
    else 
      die->SetTRDcorrectionFilename(pidTab.Data());
  }
  
  // cut setup
  SetupTrackCuts(isESD,die,cutDefinition);
  SetupPairCuts(die,cutDefinition);
  
  // histogram setup
  if(cutDefinition == kTOFTRD) 
    InitHistograms(die,cutDefinition);
  
  // setup eta correction
  SetEtaCorrection();
  
  return die;
}

//______________________________________________________________________________________
void SetupTrackCuts(Bool_t isESD, AliDielectron *die, Int_t cutDefinition)
{
  //
  // Setup the track cuts
  //
  
  // Quality cuts
  AliDielectronCutGroup* cuts = new AliDielectronCutGroup("cuts","cuts",AliDielectronCutGroup::kCompAND);
  die->GetTrackFilter().AddCuts(cuts);
 
  //Pt cut, should make execution a bit faster
  AliDielectronVarCuts *pt = new AliDielectronVarCuts("Pt>.85","Pt>.85");
  pt->AddCut(AliDielectronVarManager::kPt,0.85,1e30);
  cuts->AddCut(pt);
  
  // track cuts ESD and AOD
  AliDielectronVarCuts *varCuts = new AliDielectronVarCuts("VarCuts","VarCuts");
  varCuts->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
  varCuts->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
  varCuts->AddCut(AliDielectronVarManager::kEta,         -0.9,   0.9);
  varCuts->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
  varCuts->AddCut(AliDielectronVarManager::kNclsTPC,     70.0, 160.0);
  varCuts->AddCut(AliDielectronVarManager::kKinkIndex0,   0.0);

   AliDielectronTrackCuts *trkCuts = new AliDielectronTrackCuts("TrkCuts","TrkCuts");
   varCuts->AddCut(AliDielectronVarManager::kITSLayerFirstCls,-0.01,1.5); //ITS(0-1) = SPDany
 
  cuts->AddCut(varCuts);
 
  trkCuts->SetRequireITSRefit(kTRUE);
  trkCuts->SetRequireTPCRefit(kTRUE);
  cuts->AddCut(trkCuts);
  
  //Do we have an MC handler?
  //  Bool_t hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);
  
  /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv PID CUTS vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
  AliDielectronPID *pid = new AliDielectronPID("PID","PID");
  
  ////////////////////////////////// DATA
  if(!hasMC) {
    
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,-3.0,3.0,0.,0.,kTRUE);
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kProton,-3.0,3.0,0.,0.,kTRUE);
  }
  
  ////////////////////////////////// MC
  if(hasMC) {
    
    // electron
    Double_t nSigmaPi = 3.5; Double_t nSigmaP = 3.5;
    Double_t resolution=0.0549;
    Double_t BBpro[5] = {0};
    Double_t BBpio[5] = {0};
    
    for(Int_t icent=0; icent<8; icent++) {
      
      switch (icent) {
        case 0:  // 0-10%
          BBpro[0] = 0.031555;  BBpro[1] = 26.0595; BBpro[2] = 3.02422e-11;  BBpro[3] = 2.05594; BBpro[4] = 5.99848;
          BBpio[0] = 0.0252122; BBpio[1] = 38.8991; BBpio[2] = 4.0901e-11;   BBpio[3] = 5.27988; BBpio[4] = 4.3108;
          break;
        case 1:  // 10-20%
          BBpro[0] = 0.0315171; BBpro[1] = 25.8656; BBpro[2] = 3.03896e-11;  BBpro[3] = 2.05802; BBpro[4] = 5.99999;
          BBpio[0] = 0.0252127; BBpio[1] = 33.8617; BBpio[2] = 3.56866e-11;  BBpio[3] = 5.24831; BBpio[4] = 4.31093;
          break;
        case 2:  // 20-30%
          BBpro[0] = 0.0315171; BBpro[1] = 25.8656; BBpro[2] = 3.03896e-11;  BBpro[3] = 2.05802; BBpro[4] = 5.99999;
          BBpio[0] = 0.0263205; BBpio[1] = 37.9307; BBpio[2] = 4.29724e-11;  BBpio[3] = 5.74458; BBpio[4] = 4.32459;
          break;
        case 3:  // 30-40%
          BBpro[0] = 0.0315171; BBpro[1] = 25.8656; BBpro[2] = 3.03896e-11;  BBpro[3] = 2.05802; BBpro[4] = 5.99999;
          BBpio[0] = 0.026294;  BBpio[1] = 39.0346; BBpio[2] = 4.12261e-11;  BBpio[3] = 5.28808; BBpio[4] = 4.31301;
          break;
        case 4:  // 40-50%
          BBpro[0] = 0.0315171; BBpro[1] = 25.8656; BBpro[2] = 3.03896e-11;  BBpro[3] = 2.05802; BBpro[4] = 5.99999;
          BBpio[0] = 0.0263134; BBpio[1] = 38.2084; BBpio[2] = 3.75159e-11;  BBpio[3] = 5.78125; BBpio[4] = 4.31363;
          break;
        case 5:  // 50-60%
          BBpro[0] = 0.0315171; BBpro[1] = 25.8656; BBpro[2] = 3.03896e-11;  BBpro[3] = 2.05802; BBpro[4] = 5.99999;
          BBpio[0] = 0.0263134; BBpio[1] = 38.2084; BBpio[2] = 3.75159e-11;  BBpio[3] = 5.78125; BBpio[4] = 4.31363;
          break;
        case 6:  // 60-70%
          BBpro[0] = 0.031555;  BBpro[1] = 26.0595; BBpro[2] = 3.02422e-11;  BBpro[3] = 2.05594; BBpro[4] = 5.99848;
          BBpio[0] = 0.026302;  BBpio[1] = 38.6888; BBpio[2] = 3.56792e-11;  BBpio[3] = 5.2465;  BBpio[4] = 4.31094;
          break;
        case 7:  // 70-80%
          BBpro[0] = 0.0315171; BBpro[1] = 25.8656; BBpro[2] = 3.03896e-11;  BBpro[3] = 2.05802; BBpro[4] = 5.99999;
          BBpio[0] = 0.0263134; BBpio[1] = 38.2084; BBpio[2] = 3.75159e-11;  BBpio[3] = 5.78125; BBpio[4] = 4.31363;
          break;
        case 8:  // 80-90%
          BBpro[0] = 0.0313438; BBpro[1] = 25.8666; BBpro[2] = 4.5457e-11;   BBpro[3] = 2.07912; BBpro[4] = 5.99986;
          BBpio[0] = 0.0252127; BBpio[1] = 33.8617; BBpio[2] = 3.56866e-11;  BBpio[3] = 5.24831; BBpio[4] = 4.31093;
          break;
        case 9:  // 90-100%
          BBpro[0] = 0.0319126; BBpro[1] = 36.8784; BBpro[2] = 3.4274e-11;   BBpro[3] = 3.2431;  BBpro[4] = 5.93388;
          BBpio[0] = 0.027079;  BBpio[1] = 67.5936; BBpio[2] = 9.72548e-11;  BBpio[3] = 9.61382; BBpio[4] = 5.99372;
          break;
      }
      
      
      TF1 *ffPro=new TF1(Form("fBethe%d_c%d",AliPID::kProton,icent), Form("(%f*%f+(AliExternalTrackParam::BetheBlochAleph(x/%f,[0],[1],[2],[3],[4])-AliExternalTrackParam::BetheBlochAleph(x/%f,[0],[1],[2],[3],[4])))/%f", nSigmaP,resolution, AliPID::ParticleMass(AliPID::kProton), AliPID::ParticleMass(AliPID::kElectron), resolution), 0.05,200.);
      
      TF1 *ffPio=new TF1(Form("fBethe%d_c%d",AliPID::kPion,icent), Form("(%f*%f+(AliExternalTrackParam::BetheBlochAleph(x/%f,[0],[1],[2],[3],[4])-AliExternalTrackParam::BetheBlochAleph(x/%f,[0],[1],[2],[3],[4])))/%f", nSigmaPi,resolution, AliPID::ParticleMass(AliPID::kPion), AliPID::ParticleMass(AliPID::kElectron), resolution), 0.05,200.);
      
      TString list=gSystem->Getenv("LIST");
      
      //LHC11a10b
      if (list.Contains("LHC11a10b") || list.IsNull()) {
        printf("LHC11a10b parameters\n");
        ffPro->SetParameters(BBpro[0],BBpro[1],BBpro[2],BBpro[3],BBpro[4]);
        ffPio->SetParameters(BBpio[0],BBpio[1],BBpio[2],BBpio[3],BBpio[4]);
        
        // proton cut
        pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,ffPro,10,((double)icent)*10.,((double)icent+1)*10,
                    kFALSE,AliDielectronPID::kRequire,AliDielectronVarManager::kCentrality);
        // pion cut
        pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,ffPio,10,((double)icent)*10.,((double)icent+1)*10,
                    kFALSE,AliDielectronPID::kRequire,AliDielectronVarManager::kCentrality);
      }
    }
    
    // shifts for the nSigma electrons
    TGraph* nSigmaCorrection = new TGraph();
    // LHC11a10b
    if (list.Contains("LHC11a10b") || list.IsNull()) {
      nSigmaCorrection->SetPoint(0, 137161., -0.50-(0.28));
      nSigmaCorrection->SetPoint(1, 139510., -0.50-(0.28));
      pid->SetCorrGraph(nSigmaCorrection);
    }
    
  } //hasMC
  
  ////////////////////////////////// DATA + MC
  // pid cuts TPC + TOF 
  pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-4.,4.);
    pid->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,-3,3.,0.,0.,kFALSE,AliDielectronPID::kIfAvailable);
  
  cuts->AddCut(pid);  
  /* ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ PID CUTS ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
  
  
  // exclude conversion electrons selected by the tender
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
  
  // conversion rejection
  Double_t gCut = 0.05;             // default
  
  AliDielectronVarCuts *gammaCut=new AliDielectronVarCuts("gammaCut","gammaCut");
  gammaCut->AddCut(AliDielectronVarManager::kM,0.,gCut);
  die->GetPairPreFilter().AddCuts(gammaCut);
  
  
}

//______________________________________________________________________________________
void InitHistograms(AliDielectron *die, Int_t cutDefinition)
{
  //
  // Initialise the histograms
  //
  //  Bool_t hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);
  
  //Setup histogram Manager
  AliDielectronHistos *histos=new AliDielectronHistos(die->GetName(),die->GetTitle());
  
  
  //add histograms to event class
  histos->AddClass("Event");
  histos->AddClass("Event_noCuts");
  histos->UserHistogram("Event","VtxZ","Vertex Z;z (cm)",
                        300,-15.,15.,AliDielectronVarManager::kZvPrim);
 histos->UserHistogram("Event_noCuts","VtxZ","Vertex Z;z (cm)",
                        300,-15.,15.,AliDielectronVarManager::kZvPrim);
 
  histos->UserHistogram("Event","Centrality","Centrality;centrality (%)",
                        20,0.,100.,AliDielectronVarManager::kCentrality);

  histos->UserHistogram("Event","Multiplicity","Multiplicity V0;Multiplicity V0",
                        500,0.,25000., AliDielectronVarManager::kMultV0);
  histos->UserHistogram("Event","Multiplicity_nTracks","Multiplicity V0 vs #tracks; #tracks; Multiplicity V0;",
                        500,0.,25000.,500,0.,25000.,AliDielectronVarManager::kNTrk, AliDielectronVarManager::kMultV0);

  histos->UserHistogram("Event","Cent_Mult","Centrality vs. Multiplicity;centrality (%);Multiplicity V0",
                        10,0.,100., 500,0.,500.,AliDielectronVarManager::kCentrality,AliDielectronVarManager::kMultV0);
  histos->UserProfile("Event","Cent_Nacc",
                      "accepted tracks;centrality (%)",
                      AliDielectronVarManager::kNacc,
                      "0.,5.,10.,20.,40.,50.,60.,80.,100.",
                      AliDielectronVarManager::kCentrality);
  histos->UserProfile("Event","Cent_NVtxContrib",
                      "number of vertex contributors;centrality (%)",
                      AliDielectronVarManager::kNVtxContrib,
                      "0.,5.,10.,20.,40.,50.,60.,80.,100.",
                      AliDielectronVarManager::kCentrality);
  
  
 
  
  
  //Initialise histogram classes
  histos->SetReservedWords("Track;Pair");
  
  
  //Pair classes
  // to fill also mixed event histograms loop until 10
  for (Int_t i=0; i<3; ++i){
    histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(i)));
  }
  
  if(cutDefinition <= kTOFTRD) {
    
    //legs from pair
    for (Int_t i=0; i<3; ++i){
      histos->AddClass(Form("Track_Legs_%s",AliDielectron::PairClassName(i)));
    }
    
    //Track classes
    //to fill also track info from 2nd event loop until 2
    for (Int_t i=0; i<2; ++i){
      histos->AddClass(Form("Track_%s",AliDielectron::TrackClassName(i)));
    }
    
    //add histograms to Track classes
    //histos->UserHistogram("Track","TOFbit","TOFbit;bit;#tracks",19,-9.5,9.5,AliDielectronVarManager::kTOFPIDBit);
    histos->UserHistogram("Track","Pt","Pt;Pt [GeV];#tracks",
                          400,0,20.,
                          AliDielectronVarManager::kPt);
   
    histos->UserHistogram("Track","TPCnCls","Number of Clusters TPC;TPC number clusteres;#tracks",
                          160,-0.5,159.5,
                          AliDielectronVarManager::kNclsTPC);
    histos->UserHistogram("Track","TPCsignalN","Number of Clusters TPC;TPC number clusteres;#tracks",
                          160,-0.5,159.5,
                          AliDielectronVarManager::kTPCsignalN);
    
    histos->UserHistogram("Track","dXY","dXY;dXY [cm];#tracks",
                          500,-1.,1.,
                          AliDielectronVarManager::kImpactParXY);
    histos->UserHistogram("Track","dZ","dZ;dZ [cm];#tracks",
                          600,-3.,3.,
                          AliDielectronVarManager::kImpactParZ);

    histos->UserHistogram("Track","Centrality_Eta_Nsigma","Cent_Eta_nSigma; Centrality; #eta; TPCnSigma_Electron",20,0.,100.,200,-1.,1.,200,-10.,10.,AliDielectronVarManager::kCentrality,AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEle);

    histos->UserHistogram("Track","Eta_Phi","Eta Phi Map; Eta; Phi;#tracks",
                          200,-1,1,200,0,6.285,
                          AliDielectronVarManager::kEta,AliDielectronVarManager::kPhi);
    
    histos->UserHistogram("Track","dEdx_P","dEdx;P [GeV];TPC signal (arb units);#tracks",
                          400,0.2,20.,200,0.,200.,
                          AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCsignal,kTRUE);

    histos->UserHistogram("Track","TPCnSigmaEle_P","TPC number of sigmas Electrons;P [GeV];TPC number of sigmas Electrons;#tracks",400,0.2,20.,200,-10.,10., AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,kTRUE);

    histos->UserHistogram("Track","nSigmaTOF_P","TOF sigmas Electrons;P [GeV];TOF number of sigmas Electrons;#tracks", 400,0.2,20.,200,-10.,10., AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaEle,kTRUE);
    
    histos->UserHistogram("Track","Ncl",";Number clusters TPC;Number clusters TPC",
                          160,-0.5,159.5,
                          AliDielectronVarManager::kNclsTPC);
    histos->UserHistogram("Track","NclFr",";Number of findable clusters (robust);Number findable clusters TPC",
                          160,-0.5,159.5,
                          AliDielectronVarManager::kNFclsTPCr);
    histos->UserHistogram("Track","Ncl_NclFr","Number of (findable) clusters TPC;found clusters;findable clusters",
                          160,-0.5,159.5,160,-0.5,159.5,
                          AliDielectronVarManager::kNclsTPC,AliDielectronVarManager::kNFclsTPCr);
    histos->UserHistogram("Track","NtrklTRD",";Number tracklets TRD for pid;Number tracklets TRD",
                          8,-0.5,7.5,
                          AliDielectronVarManager::kTRDpidQuality);
    
    //add histograms to Pair classes
    histos->UserHistogram("Pair","InvMass","Inv.Mass;Inv. Mass [GeV];#pairs",
                          300,.0,300*0.04,
                          AliDielectronVarManager::kM); // 40MeV bins, 12GeV/c2
    histos->UserHistogram("Pair","Rapidity","Rapidity;Rapidity;#pairs",
                          100,-1.,1.,
                          AliDielectronVarManager::kY);
    histos->UserHistogram("Pair","OpeningAngle","Opening angle;angle",
                          100,0.,3.15,
                          AliDielectronVarManager::kOpeningAngle);
    histos->UserHistogram("Pair","Chi2NDF","#chi^{2}/NDF;#chi^{2}/NDF",
                          100,0.,20,
                          AliDielectronVarManager::kChi2NDF);
  
   histos->UserHistogram("Pair","PseudoProperTime","Pseudoproper decay length; pseudoproper-decay-length[#mum];Entries/40#mum",
                          150,-0.3.,0.3,AliDielectronVarManager::kPseudoProperTime);


   }
  
  die->SetHistogramManager(histos);
}


void AddMCSignals(AliDielectron *die){
  //Do we have an MC handler?
  //Bool_t hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);
  if (!hasMC) return;
  
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
  promptJpsi->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  promptJpsi->SetCheckBothChargesLegs(kTRUE,kTRUE);
  promptJpsi->SetCheckBothChargesMothers(kTRUE,kTRUE);
  promptJpsi->SetCheckBothChargesGrandMothers(kTRUE,kTRUE);
  die->AddSignalMC(promptJpsi);
  
  AliDielectronSignalMC* beautyJpsi = new AliDielectronSignalMC("beautyJpsi","Beauty J/psi");
  beautyJpsi->SetLegPDGs(11,-11);
  beautyJpsi->SetMotherPDGs(443,443);
  beautyJpsi->SetMothersRelation(AliDielectronSignalMC::kSame);
  beautyJpsi->SetGrandMotherPDGs(500,500);
  beautyJpsi->SetFillPureMCStep(kTRUE);
  beautyJpsi->SetCheckBothChargesLegs(kTRUE,kTRUE);
  beautyJpsi->SetCheckBothChargesMothers(kTRUE,kTRUE);
  beautyJpsi->SetCheckBothChargesGrandMothers(kTRUE,kTRUE);
  die->AddSignalMC(beautyJpsi);
  
  AliDielectronSignalMC* directJpsi = new AliDielectronSignalMC("directJpsi","Direct J/psi");   // embedded J/psi
  directJpsi->SetLegPDGs(11,-11);
  directJpsi->SetMotherPDGs(443,443);
  directJpsi->SetMothersRelation(AliDielectronSignalMC::kSame);
  directJpsi->SetFillPureMCStep(kTRUE);
  directJpsi->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  directJpsi->SetMotherSources(AliDielectronSignalMC::kDirect, AliDielectronSignalMC::kDirect);
  directJpsi->SetCheckBothChargesLegs(kTRUE,kTRUE);
  directJpsi->SetCheckBothChargesMothers(kTRUE,kTRUE);
  die->AddSignalMC(directJpsi);
}

void SetEtaCorrection()
{

  if (AliDielectronPID::GetEtaCorrFunction()) return;
  
  //TString list=gSystem->Getenv("LIST");
  TString list="LHC11h.pass2";

  //TString etaMap="$TRAIN_ROOT/jpsi_JPSI/EtaCorrMaps.root";
  //TString trainRoot=gSystem->Getenv("TRAIN_ROOT");
  //if (trainRoot.IsNull()) 
  TString etaMap="$ALICE_PHYSICS/PWGDQ/dielectron/files/EtaCorrMaps.root";
  if (gSystem->AccessPathName(gSystem->ExpandPathName(etaMap.Data()))){
    Error("ConfigPbPb","Eta map not found: %s",etaMap.Data());
    return;
  }

  TFile f(etaMap.Data());
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

