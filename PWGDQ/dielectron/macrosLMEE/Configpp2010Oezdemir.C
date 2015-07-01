void InitHistograms(AliDielectron *die, Int_t cutDefinition);
void InitCF(AliDielectron* die, Int_t cutDefinition);
void AddMCSignals(AliDielectron* die);
void SetupCuts(AliDielectron *die, Int_t cutDefinition);
void InitHF(AliDielectron* die, Int_t cutDefinition);
AliDielectronEventCuts *GetEventCuts();
TVectorD *BinsToVector(Int_t nbins, Double_t min, Double_t max);
TVectorD *GetVector(Int_t var);
enum {kMee=0, kMee500, kPtee, kP2D, kPhiV, kOpAng, kOpAng2, kEta2D, kEta3D, kPhi2D, kY3D, kSigmaEle, kSigmaOther, kTPCdEdx};

AliESDtrackCuts *SetupESDtrackCutsAna(Int_t cutDefinition);
AliDielectronPID *SetPIDcutsAna(Int_t cutDefinition);
AliESDtrackCuts *SetupESDtrackCutsPre(Int_t cutDefinition);
AliDielectronPID *SetPIDcutsPre(Int_t cutDefinition);

TString names= ("PairPre60MeV50mrad;PairPre100MeV50mrad");

Bool_t kMix = kFALSE;
Bool_t isQAtask=kFALSE;//needed for InitHistograms() and GetVector()
Bool_t doRejectionStep=kFALSE;//needed for SetupCuts() and InitHistograms()
// _____ for Random Rejection task
Bool_t isRandomRejTask=kFALSE;//needed for InitHistograms()

TObjArray *arrNames=names.Tokenize(";");
const Int_t nDie=arrNames->GetEntriesFast();

AliDielectron* Configpp2010(Int_t cutDefinition, Bool_t isRandomRej=kFALSE)
{
  //
  // Setup the instance of AliDielectron
  //

  isRandomRejTask=isRandomRej;
  kMix=kTRUE;
  doRejectionStep=kTRUE;
  isQAtask=kFALSE;

  // create the actual framework object
  TString name=Form("%02d",cutDefinition);
  if (cutDefinition<arrNames->GetEntriesFast()){
    name=arrNames->At(cutDefinition)->GetName();
  }

  Bool_t hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);
  AliDielectron *die = new AliDielectron(Form("%s",name.Data()),Form("Track cuts: %s",name.Data()));
  if(hasMC) AddMCSignals(die);
  
  // set track-PID and pair cuts
  SetupCuts(die,cutDefinition);

  if(kMix){
        AliDielectronMixingHandler *mix = new AliDielectronMixingHandler;
	mix->SetMixType(AliDielectronMixingHandler::kAll);
	mix->AddVariable(AliDielectronVarManager::kZvPrim,"-10., -7.5, -5., -2.5 , 0., 2.5, 5., 7.5 , 10.");
        mix->AddVariable(AliDielectronVarManager::kNacc,"0,10000");
	mix->SetDepth(17);
	die->SetMixingHandler(mix); 
        }//kMix

  //
  // histogram setup
  // only if an AliDielectronHistos object is attached to the
  // dielectron framework histograms will be filled
  //

  InitHistograms(die,cutDefinition);
  //InitCF(die,cutDefinition);
  //InitHF(die,cutDefinition);
return die;

}

//______________________________________________________________________________________
void SetupCuts(AliDielectron *die, Int_t cutDefinition)
{
  //
  // Setup the track cuts
  //
        //options
        //die->SetPreFilterAllSigns();
        die->SetPreFilterUnlikeOnly();    
       
        //pairing with TLorentzVector
        die->SetUseKF(kFALSE);
              
 
        if (doRejectionStep){
        //trackPre cuts
        die->GetTrackFilter().AddCuts(SetupESDtrackCutsPre(cutDefinition));
        //pidPre cuts
//        die->GetTrackFilter().AddCuts(SetPIDcutsPre(cutDefinition));
        //pairPre cuts
        AliDielectronVarCuts *PairPre = new AliDielectronVarCuts("PairPre","PairPre");//mass and Phiv together
          if(cutDefinition ==0){
             PairPre->AddCut(AliDielectronVarManager::kM, 0.0 , 0.06);
             PairPre->AddCut(AliDielectronVarManager::kOpeningAngle, 0., 0.05 );
             die->GetPairPreFilter().AddCuts(PairPre);
          }
          if(cutDefinition ==1){
             PairPre->AddCut(AliDielectronVarManager::kM, 0.0 , 0.1);
             PairPre->AddCut(AliDielectronVarManager::kOpeningAngle, 0., 0.05 );
             die->GetPairPreFilter().AddCuts(PairPre);
          }      
        //trackAna cuts
        die->GetPairPreFilterLegs().AddCuts(SetupESDtrackCutsAna(cutDefinition));
        //pidAna cuts
        die->GetPairPreFilterLegs().AddCuts(SetPIDcutsAna(cutDefinition));
        //pair cuts
        AliDielectronTrackCuts *noconv=new AliDielectronTrackCuts("noConv","conversion tagging");
        noconv->SetV0DaughterCut(AliPID::kElectron,kTRUE);//real cut
        //      noconv->SetV0DaughterCut(AliPID::kElectron,kFALSE);//V0 analysis
        die->GetPairPreFilterLegs().AddCuts(noconv);
       
        AliDielectronVarCuts *PhiV = new AliDielectronVarCuts("PhiV","PhiV");//mass and Phiv together
        PhiV->AddCut(AliDielectronVarManager::kM, 0.0 , 0.1);
        PhiV->AddCut(AliDielectronVarManager::kPhivPair, 0. , 2. );
        
        AliDielectronVarCuts *PhiV2 = new AliDielectronVarCuts("PhiV2","PhiV2");//mass and Phiv together
        PhiV2->AddCut(AliDielectronVarManager::kM, 0.1 , 1000.);

        AliDielectronCutGroup* pairCutsCG2 =new AliDielectronCutGroup("pairCutsCG2","pairCutsCG2",AliDielectronCutGroup::kCompOR);
        pairCutsCG2->AddCut(PhiV);
        pairCutsCG2->AddCut(PhiV2);
        die->GetPairFilter().AddCuts(pairCutsCG2);
        }
        else{
        
        //trackAna cuts
        die->GetTrackFilter().AddCuts(SetupESDtrackCutsAna(cutDefinition));
        //pidAna cuts
        die->GetTrackFilter().AddCuts(SetPIDcutsAna(cutDefinition));
        //pair cuts
        AliDielectronTrackCuts *noconv=new AliDielectronTrackCuts("noConv","conversion tagging");
        noconv->SetV0DaughterCut(AliPID::kElectron,kTRUE);//real cut
        //      noconv->SetV0DaughterCut(AliPID::kElectron,kFALSE);//V0 analysis
        die->GetTrackFilter.AddCuts(noconv);

        AliDielectronVarCuts *PhiV = new AliDielectronVarCuts("PhiV","PhiV");//mass and Phiv together
        PhiV->AddCut(AliDielectronVarManager::kM, 0.0 , 0.1);
        PhiV->AddCut(AliDielectronVarManager::kPhivPair, 0. , 2.);

        AliDielectronVarCuts *PhiV2 = new AliDielectronVarCuts("PhiV2","PhiV2");//mass and Phiv together
        PhiV2->AddCut(AliDielectronVarManager::kM, 0.1 , 1000.);

        AliDielectronCutGroup* pairCutsCG2 =new AliDielectronCutGroup("pairCutsCG2","pairCutsCG2",AliDielectronCutGroup::kCompOR);
        pairCutsCG2->AddCut(PhiV);
        pairCutsCG2->AddCut(PhiV2);
        die->GetPairFilter().AddCuts(pairCutsCG2);        
        }
}
//______________________________________________________________________________________
//-----------------------------------pid------------------------------------------------

AliDielectronPID *SetPIDcutsPre(Int_t cutDefinition){

  AliDielectronPID *pidPre = new AliDielectronPID();
  //stronger electron inclusion

    //pidPre->AddCut(AliDielectronPID::kITS,AliPID::kElectron,-5.,1.5,0.0,100.,kFALSE,AliDielectronPID::kRequire,AliDielectronVarManager::kPt);
//    pidPre->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,-4.,4.,0., 100., kFALSE, AliDielectronPID::kIfAvailable,AliDielectronVarManager::kPt);
    //pidPre->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-5.,5.,0.,100.,kFALSE,AliDielectronPID::kRequire,AliDielectronVarManager::kPt);
    //pidPre->AddCut(AliDielectronPID::kTPC,AliPID::kPion,-100.,4.,0.2,10.,kTRUE,AliDielectronPID::kRequire,AliDielectronVarManager::kPt);
   
  
    //ITS+TPC+TOFifAvailable

  return pidPre;

}

AliDielectronPID *SetPIDcutsAna(Int_t cutDefinition){

  AliDielectronPID *pidAna = new AliDielectronPID();
  //stronger electron inclusion
    pidAna->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,4. ,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pidAna->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,  -3. ,3. ,0.0, 100., kFALSE,AliDielectronPID::kIfAvailable,AliDielectronVarManager::kPt);
    pidAna->AddCut(AliDielectronPID::kITS,AliPID::kElectron,  -5. ,1. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pidAna->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -1.5,3. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
   
  //ITS+TPC+TOFifAvailable

  return pidAna;

}


//______________________________________________________________________________________
AliESDtrackCuts *SetupESDtrackCutsPre(Int_t cutDefinition)
{

  AliESDtrackCuts *fesdTrackCutsPre = new AliESDtrackCuts;


  //global
  fesdTrackCutsPre->SetPtRange( 0.08 , 100. );
  fesdTrackCutsPre->SetEtaRange( -1.1 , 1.1 );
  fesdTrackCutsPre->SetAcceptKinkDaughters(kFALSE);
  fesdTrackCutsPre->SetRequireSigmaToVertex(kFALSE);
  fesdTrackCutsPre->SetDCAToVertex2D(kFALSE);
  fesdTrackCutsPre->SetMaxDCAToVertexZ(3.);
  fesdTrackCutsPre->SetMaxDCAToVertexXY(1.);

  //ITS
  fesdTrackCutsPre->SetRequireITSRefit(kTRUE);
  fesdTrackCutsPre->SetMinNClustersITS(3);
  fesdTrackCutsPre->SetMaxChi2PerClusterITS(5);
 
  return fesdTrackCutsPre;


}


AliESDtrackCuts *SetupESDtrackCutsAna(Int_t cutDefinition)
{

  AliESDtrackCuts *fesdTrackCutsAna = new AliESDtrackCuts;


  //global
  fesdTrackCutsAna->SetPtRange( 0.2 , 100. );
  fesdTrackCutsAna->SetEtaRange( -0.8 , 0.8 );
  fesdTrackCutsAna->SetAcceptKinkDaughters(kFALSE);
  fesdTrackCutsAna->SetRequireSigmaToVertex(kFALSE);
  fesdTrackCutsAna->SetDCAToVertex2D(kFALSE);
  fesdTrackCutsAna->SetMaxDCAToVertexZ(3.);
  fesdTrackCutsAna->SetMaxDCAToVertexXY(1.);

  //ITS
  fesdTrackCutsAna->SetRequireITSRefit(kTRUE);
  fesdTrackCutsAna->SetMinNClustersITS(4);
  fesdTrackCutsAna->SetMaxChi2PerClusterITS(5);
  fesdTrackCutsAna->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);

  //TPC
  fesdTrackCutsAna->SetRequireTPCRefit(kTRUE);
  fesdTrackCutsAna->SetMinNClustersTPC(80);
  fesdTrackCutsAna->SetMinNCrossedRowsTPC(100);
  fesdTrackCutsAna->SetMinRatioCrossedRowsOverFindableClustersTPC(0.5);
  fesdTrackCutsAna->SetMaxChi2PerClusterTPC(4);   
 
  return fesdTrackCutsAna;


}

//______________________________________________________________________________________
void InitHistograms(AliDielectron *die, Int_t cutDefinition)
{
  //
  // Initialise the histograms
  //
  
  //Setup histogram Manager
  AliDielectronHistos *histos=new AliDielectronHistos(die->GetName(),die->GetTitle());
  
  //Initialise histogram classes
  histos->SetReservedWords("Track;Pair;Track_Legs;Pre;RejPair;RejTrack;Random");
  
  //Event class
  histos->AddClass("Event"); // all classes will be stored in 'THashList fHistoList'
  
  //Track classes
  //to fill also track info from 2nd event loop until 3
  // in AliDielectron.cxx: fgkTrackClassNames[4] = {"ev1+","ev1-","ev2+","ev2-"};
  for (Int_t i=0; i<2; ++i){
    histos->AddClass(Form("Track_%s",AliDielectron::TrackClassName(i)));
  }
  
  //Pair classes
  // to fill also mixed event histograms loop until 10
  // fgkPairClassNames[11] = {
  //  "ev1+_ev1+",  "ev1+_ev1-",  "ev1-_ev1-",  // 0-2 (same event)
  //  "ev1+_ev2+",  "ev1-_ev2+",  "ev2+_ev2+",  // 3-4 (+5)
  //  "ev1+_ev2-",  "ev1-_ev2-",                // 6-7
  //  "ev2+_ev2-",  "ev2-_ev2-",  "ev1+_ev1-_TR"
  // };
  
  if (!isQAtask && !isRandomRejTask) // -> analysis with pairing
  {
    for (Int_t i=0; i<3; ++i){
      histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(i)));
      // Legs of final Pairs. Both charges together. No duplicate entries.
      if (!isQAtask) {
        histos->AddClass(Form("Track_Legs_%s",AliDielectron::PairClassName(i))); // not TrackClassName, see 'AliDielectron::FillHistograms(...)'
      }
    }
    
    //Mixed event and track rotation
    if (die->GetMixingHandler()) {
      histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(3)));
      histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(4)));
      histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(6)));
      histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(7)));
    }
    if (die->GetTrackRotator()) {
      histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(10)));
    }
    
    if (doRejectionStep) 
    {
      //PreFilter Classes
      //to fill also track info from 2nd event loop until 2
      for (Int_t i=0; i<2; ++i){
        histos->AddClass(Form("Pre_%s",AliDielectron::TrackClassName(i)));
        // class 'Pre_%s': "Fill Histogram information for tracks after prefilter"(AliDielectron::FillHistogramsTracks(...)), so it is identical to Track_%s if no paircut is used.
      }
      
      //Create Classes for Rejected Tracks/Pairs:
      for (Int_t i=0; i<3; ++i){
        histos->AddClass(Form("RejPair_%s",AliDielectron::PairClassName(i)));
        // Legs of rejected Pairs. Both charges together. One track can and will make multiple entries.
        histos->AddClass(Form("RejTrack_%s",AliDielectron::PairClassName(i))); // not TrackClassName, see 'AliDielectron::FillHistogramsPair(...)'
      }
      
      /*
       //track rotation
       histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(AliDielectron::kEv1PMRot)));
       histos->AddClass(Form("Track_Legs_%s",AliDielectron::PairClassName(AliDielectron::kEv1PMRot)));
       */
    }// end: (doRejectionStep)
    
  }// end: (!isQAtask)
  
   
  if (isRandomRejTask)
  {
    //
    // _____ histograms for AliAnalysisTaskMultiDielectronPR _____
    //
    histos->AddClass("Rand_Pair");
    histos->AddClass("Rand_RejPair");
    const char* cRandomPairClassNames[4] = { "Testpart", "DataEle", "RejTestpart", "RejDataEle" };
    for (Int_t i=0; i<4; ++i){
      histos->AddClass(Form("Random_%s",cRandomPairClassNames[i]));
    }
    histos->UserHistogram("Random","Pt","",200,0,10.,AliDielectronVarManager::kPt);
    histos->UserHistogram("Random","Eta","",200,-2,2,AliDielectronVarManager::kEta);
    histos->UserHistogram("Random","Phi","",120,0.,TMath::TwoPi(),AliDielectronVarManager::kPhi);
    histos->UserHistogram("Random","Px","",200,0,10.,AliDielectronVarManager::kPx);
    histos->UserHistogram("Random","Py","",200,0,10.,AliDielectronVarManager::kPy);
    histos->UserHistogram("Random","Pz","",200,0,10.,AliDielectronVarManager::kPz);
    histos->UserHistogram("Random","Pt_Eta_Phi","",GetVector(kP2D),GetVector(kEta2D),GetVector(kPhi2D), AliDielectronVarManager::kPt,AliDielectronVarManager::kEta,AliDielectronVarManager::kPhi);
    histos->UserHistogram("Random","dXY_dZ","",100,-1.,1.,150,-3.,3.,AliDielectronVarManager::kImpactParXY,AliDielectronVarManager::kImpactParZ);
    histos->UserHistogram("Random","TPC_dEdx_Pt","",
                          100,0.,2.,100,0.,100., AliDielectronVarManager::kPt,AliDielectronVarManager::kTPCsignal);
  }
  
	//add histograms to event class
	histos->UserHistogram("Event","nEvents","",1,0.,1.,AliDielectronVarManager::kNevents);
  histos->UserHistogram("Event","nESDTracks","",8000,0,80000,AliDielectronVarManager::kNTrk);
  histos->UserHistogram("Event","Nacc","",8000,0,8000,AliDielectronVarManager::kNacc);
  histos->UserHistogram("Event","RefMultTPConly","",8000,0,8000,AliDielectronVarManager::kRefMultTPConly);
  histos->UserHistogram("Event","epTPC","",240,-TMath::Pi(),TMath::Pi(),AliDielectronVarManager::kTPCrpH2uc);
  histos->UserHistogram("Event","epV0AC","",240,-TMath::Pi(),TMath::Pi(),AliDielectronVarManager::kv0ACrpH2);
  histos->UserHistogram("Event","epV0AC_epTPC","",120,-TMath::PiOver2(),TMath::PiOver2(),120,-TMath::PiOver2(),TMath::PiOver2(),AliDielectronVarManager::kTPCrpH2uc,AliDielectronVarManager::kv0ACrpH2);
  
  
  //add histograms to Track classes
  // axis labels are set to the values in 'AliDielectronVarManager.cxx' if the histogram title is empty or starts with ';'! [see AliDielectronHistos::AdaptNameTitle(...)]
  histos->UserHistogram("Track","Pt",";Pt [GeV];#tracks",200,0,10.,AliDielectronVarManager::kPt);
  histos->UserHistogram("Track","Px",";Px [GeV];#tracks",200,0,10.,AliDielectronVarManager::kPx);
  histos->UserHistogram("Track","Py",";Py [GeV];#tracks",200,0,10.,AliDielectronVarManager::kPy);
  histos->UserHistogram("Track","Pz",";Pz [GeV];#tracks",200,0,10.,AliDielectronVarManager::kPz);
  histos->UserHistogram("Track","P_PIn",";p (GeV/c);p_{in} (GeV/c)",
                        GetVector(kP2D), GetVector(kP2D), AliDielectronVarManager::kP,AliDielectronVarManager::kPIn);
  
  // ITS
  histos->UserHistogram("Track","ITS_dEdx_P",";p (GeV/c);ITS signal (arb units)",
                        GetVector(kP2D), BinsToVector(700,0.,700.), AliDielectronVarManager::kP,AliDielectronVarManager::kITSsignal);
  histos->UserHistogram("Track","ITSnSigmaEle_P",";p (GeV/c);n#sigma_{ele}^{ITS}",
                        GetVector(kP2D), GetVector(kSigmaEle), AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaEle);
  histos->UserHistogram("Track","ITSnSigmaPio_P",";p (GeV/c);n#sigma_{pion}^{ITS}",
                        GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaPio);
  histos->UserHistogram("Track","ITSnSigmaKao_P",";p (GeV/c);n#sigma_{kaon}^{ITS}",
                        GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaKao);
  //histos->UserHistogram("Track","ITSnSigmaPro_P",";p (GeV/c);n#sigma_{proton}^{ITS}",
  //                      GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaPro);
  // TPC
  histos->UserHistogram("Track","TPC_dEdx_P",";p_{in} (GeV/c);TPC signal (arb units)",
                        GetVector(kP2D), GetVector(kTPCdEdx), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCsignal);
  histos->UserHistogram("Track","TPCnSigmaEle_P",";p_{in} (GeV/c);n#sigma_{ele}^{TPC}",
                        GetVector(kP2D), GetVector(kSigmaEle), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle);
  histos->UserHistogram("Track","TPCnSigmaEle_Eta",";Eta;n#sigma_{ele}^{TPC}",
                        GetVector(kEta2D), GetVector(kSigmaEle), AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEle);
  // some histograms with SigmaEleRaw to benchmark the dEdx-eta-correction:
  histos->UserHistogram("Track","TPCnSigmaEleRaw_Eta",";Eta;n#sigma_{ele,Raw}^{TPC}",
                        GetVector(kEta2D), GetVector(kSigmaEle), AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEleRaw);
  histos->UserHistogram("Track","TPCnSigmaEle_Nacc",";N_{acc}; n#sigma_{ele}^{TPC}",
                        BinsToVector(100,0.,5000.), GetVector(kSigmaEle), AliDielectronVarManager::kNacc,AliDielectronVarManager::kTPCnSigmaEle);
  histos->UserHistogram("Track","TPCnSigmaEle_RefMultTPConly",";N_{TPC ref}; n#sigma_{ele}^{TPC}",
                        BinsToVector(100,0.,5000.), GetVector(kSigmaEle), AliDielectronVarManager::kRefMultTPConly,AliDielectronVarManager::kTPCnSigmaEle);
  histos->UserHistogram("Track","TPCnSigmaEleRaw_RefMultTPConly",";N_{TPC ref}; n#sigma_{ele,Raw}^{TPC}",
                        BinsToVector(100,0.,5000.), GetVector(kSigmaEle), AliDielectronVarManager::kRefMultTPConly,AliDielectronVarManager::kTPCnSigmaEleRaw);
  // 3D
  histos->UserHistogram("Track","TPCnSigmaEle_Eta_P",";Eta;n#sigma_{ele}^{TPC};p_{in} (GeV/c)",
                        GetVector(kEta3D), GetVector(kSigmaEle), GetVector(kP2D),
                        AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kPIn);
  histos->UserHistogram("Track","TPCnSigmaEle_Eta_RefMultTPConly",";Eta;n#sigma_{ele}^{TPC};N_{TPC ref}",
                        GetVector(kEta3D), GetVector(kSigmaEle), BinsToVector(100,0.,5000.),
                        AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kRefMultTPConly);
  
  if (isQAtask) {
    histos->UserHistogram("Track","TPC_dEdx_Eta",";Eta;TPC signal (arb units)",
                          GetVector(kEta2D), GetVector(kTPCdEdx), AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCsignal);
    histos->UserHistogram("Track","TPC_dEdx_Eta_P",";Eta;TPC signal (arb units);p_{in} (GeV/c)",
                          GetVector(kEta3D), GetVector(kTPCdEdx), GetVector(kP2D),
                          AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCsignal,AliDielectronVarManager::kPIn);
    histos->UserHistogram("Track","TPCnSigmaEle_P_dEdx",";p_{in} (GeV/c);n#sigma_{ele}^{TPC};TPC signal (arb units)",
                          GetVector(kP2D), GetVector(kSigmaEle), GetVector(kTPCdEdx),
                          AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kTPCsignal);
    
    //histos->UserHistogram("Track","TPCnSigmaEle_Eta_P"," ... // see above
    //histos->UserHistogram("Track","TPCnSigmaEle_Eta_RefMultTPConly", ... // see above
    histos->UserHistogram("Track","TPCnSigmaEle_Eta_Nacc",";Eta;n#sigma_{ele}^{TPC};N_{acc}",
                          GetVector(kEta3D), GetVector(kSigmaEle), BinsToVector(100,0.,5000.),
                          AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kNacc);
    histos->UserHistogram("Track","RefMultTPConly_Nacc",";N_{TPC ref};N_{acc}",
                          BinsToVector(100,0.,5000.), BinsToVector(100,0.,5000.), AliDielectronVarManager::kRefMultTPConly,AliDielectronVarManager::kNacc);
    histos->UserHistogram("Track","TPCnSigmaEleRaw_Eta_RefMultTPConly",";Eta;n#sigma_{ele,Raw}^{TPC};N_{TPC ref}",
                          GetVector(kEta3D), GetVector(kSigmaEle), BinsToVector(100,0.,5000.),
                          AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEleRaw,AliDielectronVarManager::kRefMultTPConly);
    
    histos->UserHistogram("Track","TPCnSigmaPio_P",";p_{in} (GeV/c);n#sigma_{pion}^{TPC}",
                          GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaPio);
    histos->UserHistogram("Track","TPCnSigmaKao_P",";p_{in} (GeV/c);n#sigma_{kaon}^{TPC}",
                          GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaKao);
    histos->UserHistogram("Track","TPCnSigmaPro_P",";p_{in} (GeV/c);n#sigma_{proton}^{TPC}",
                          GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaPro);
    histos->UserHistogram("Track","TPCnSigmaKao_Eta",";Eta;n#sigma_{kaon}^{TPC}",
                          GetVector(kEta2D), GetVector(kSigmaOther), AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaKao);
  }
  
  
  // TOF
  histos->UserHistogram("Track","TOFbeta_P",";p_{in} (GeV/c);TOF beta",
                        GetVector(kP2D), BinsToVector(120,0.,1.2) ,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFbeta);
  histos->UserHistogram("Track","TOFnSigmaEle_P",";p_{in} (GeV/c);TOF number of sigmas Electrons",
                        GetVector(kP2D), GetVector(kSigmaEle), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaEle);
  if (isQAtask) {
    histos->UserHistogram("Track","TOFnSigmaPio_P",";p_{in} (GeV/c);TOF number of sigmas Pions",
                          GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaPio);
    histos->UserHistogram("Track","TOFnSigmaKao_P",";p_{in} (GeV/c);TOF number of sigmas Kaons",
                          GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaKao);
    histos->UserHistogram("Track","TOFnSigmaPro_P",";p_{in} (GeV/c);TOF number of sigmas Protons",
                          GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaPro);
  }    
  // 2D-PID
  if (isQAtask) {
    histos->UserHistogram("Track","PIn_TPCnSigmaEle_ITSnSigmaEle",";p_{in} (GeV/c);n#sigma_{ele}^{TPC};n#sigma_{ele}^{ITS}",
                          50,0.,2.5, 160,-12.,20., 150,-10.,20.,
                          AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kITSnSigmaEle);
    histos->UserHistogram("Track","PIn_TPCnSigmaEle_TOFnSigmaEle",";p_{in} (GeV/c);n#sigma_{ele}^{TPC};TOF number of sigmas Electrons",
                          50,0.,2.5, 160,-12.,20., 50,-5.,5.,
                          AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kTOFnSigmaEle);
  }
  
  // Eta, Phi, Pt
  histos->UserHistogram("Track","Eta","",200,-2,2,AliDielectronVarManager::kEta);
  histos->UserHistogram("Track","Phi","",120,0.,TMath::TwoPi(),AliDielectronVarManager::kPhi);
  histos->UserHistogram("Track","Eta_Phi","",GetVector(kEta2D),GetVector(kPhi2D), AliDielectronVarManager::kEta,AliDielectronVarManager::kPhi);
  histos->UserHistogram("Track","Eta_Pt","",GetVector(kP2D),GetVector(kEta2D), AliDielectronVarManager::kPt,AliDielectronVarManager::kEta);
  histos->UserHistogram("Track","Phi_Pt","",GetVector(kP2D),GetVector(kPhi2D), AliDielectronVarManager::kPt,AliDielectronVarManager::kPhi);
  
  // DCA
  histos->UserHistogram("Track","dXY","",200,-2.,2.,AliDielectronVarManager::kImpactParXY);
  histos->UserHistogram("Track","dZ" ,"",400,-4.,4.,AliDielectronVarManager::kImpactParZ);
  histos->UserHistogram("Track","dXY_dZ","",100,-1.,1.,150,-3.,3.,AliDielectronVarManager::kImpactParXY,AliDielectronVarManager::kImpactParZ);
  
  // Quality
  histos->UserHistogram("Track","TPCcrossedRowsOverFindable",";TPC crossed rows over findable clusters;#tracks",120,0.,1.2,AliDielectronVarManager::kNFclsTPCfCross);
  histos->UserHistogram("Track","TPCcrossedRows",";TPC crossed rows;#tracks",160,-0.5,159.5,AliDielectronVarManager::kNFclsTPCr);
  histos->UserHistogram("Track","TPCnCls",";TPC number clusters;#tracks",160,-0.5,159.5,AliDielectronVarManager::kNclsTPC);
  histos->UserHistogram("Track","ITSnCls",";ITS number clusters;#tracks",160,-0.5,159.5,AliDielectronVarManager::kNclsITS);
  histos->UserHistogram("Track","TPCchi2",";TPC chi2/Cl;#tracks",100,0.,10.,AliDielectronVarManager::kTPCchi2Cl);
  histos->UserHistogram("Track","ITSchi2",";ITS chi2/Cl;#tracks",100,0.,10.,AliDielectronVarManager::kITSchi2Cl);
  histos->UserHistogram("Track","NclsSFracTPC",";TPC fraction of shared clusters;#tracks",200,0,10.,AliDielectronVarManager::kNclsSFracTPC);
  histos->UserHistogram("Track","TPCclsDiff",";TPC cluster difference;#tracks",200,0,20.,AliDielectronVarManager::kTPCclsDiff);
  histos->UserHistogram("Track","TPCsignalN",";TPC number PID clusters;#tracks",160,-0.5,159.5,AliDielectronVarManager::kTPCsignalN);
  
  if (isQAtask) {
    histos->UserHistogram("Track","TPCcrossedRows_TPCnCls",";TPC number clusters;TPC crossed rows",
                          160,-0.5,159.5, 160,-0.5,159.5, AliDielectronVarManager::kNclsTPC,AliDielectronVarManager::kNFclsTPCr);
    histos->UserHistogram("Track","TPCcrossedRows_Pt",";Pt [GeV];TPC crossed rows",
                          GetVector(kP2D), BinsToVector(160,-0.5,159.5), AliDielectronVarManager::kPt,AliDielectronVarManager::kNFclsTPCr);
    histos->UserHistogram("Track","TPCcrossedRowsOverFindable_Pt",";Pt [GeV];TPC crossed rows over findable",
                          GetVector(kP2D), BinsToVector(120,0.,1.2), AliDielectronVarManager::kPt,AliDielectronVarManager::kNFclsTPCfCross);
    histos->UserHistogram("Track","TPCcrossedRowsOverFindable_Eta",";Eta;TPC crossed rows over findable",
                          100,-1,1, 120,0.,1.2, AliDielectronVarManager::kEta,AliDielectronVarManager::kNFclsTPCfCross);
    histos->UserHistogram("Track","TPCcrossedRowsOverFindable_Phi",";Phi;TPC crossed rows over findable",
                          120,0.,TMath::TwoPi(), 120,0.,1.2, AliDielectronVarManager::kPhi,AliDielectronVarManager::kNFclsTPCfCross);
  }
  
  if (!isQAtask) 
  {
    //add histograms to Pair classes
    histos->UserHistogram("Pair","InvMass","",500,0.,5.,AliDielectronVarManager::kM);
    histos->UserHistogram("Pair","PairPt","",160,0.,8., AliDielectronVarManager::kPt);
    histos->UserHistogram("Pair","Rapidity","",200,-2.,2.,AliDielectronVarManager::kY);
    histos->UserHistogram("Pair","OpeningAngle","",240,0.,TMath::Pi(),AliDielectronVarManager::kOpeningAngle);
    
    //2D and 3D histograms
    histos->UserHistogram("Pair","InvMass_PairPt",";Inv. Mass [GeV];Pair Pt [GeV];#pairs",
                          GetVector(kMee), GetVector(kPtee),
                          AliDielectronVarManager::kM, AliDielectronVarManager::kPt);
    histos->UserHistogram("Pair","Eta_Phi_Pair",";Eta;Phi;#pairs",
                          100,-1.,1., 120,0.,TMath::TwoPi(),
                          AliDielectronVarManager::kEta, AliDielectronVarManager::kPhi);
    histos->UserHistogram("Pair","InvMass_PairPt_PhivPair",";Inv. Mass [GeV];Pair Pt [GeV];PhiV",
                          GetVector(kMee), GetVector(kPtee), GetVector(kPhiV), 
                          AliDielectronVarManager::kM, AliDielectronVarManager::kPt, AliDielectronVarManager::kPhivPair);
    histos->UserHistogram("Pair","InvMass_PairPt_OpeningAngle",";Inv. Mass [GeV];Pair Pt [GeV];Opening Angle",
                          GetVector(kMee), GetVector(kPtee), GetVector(kOpAng), 
                          AliDielectronVarManager::kM, AliDielectronVarManager::kPt, AliDielectronVarManager::kOpeningAngle);
    histos->UserHistogram("Pair","InvMass_PhivPair_OpeningAngle",";Inv. Mass [GeV];PhiV;Opening Angle",
                          GetVector(kMee500), GetVector(kPhiV), GetVector(kOpAng2), 
                          AliDielectronVarManager::kM, AliDielectronVarManager::kPhivPair, AliDielectronVarManager::kOpeningAngle);
    histos->UserHistogram("Pair","InvMass_PairPt_Rapidity",";Inv. Mass [GeV];Pair Pt [GeV];Rapidity",
                          GetVector(kMee), GetVector(kPtee), GetVector(kY3D), 
                          AliDielectronVarManager::kM, AliDielectronVarManager::kPt, AliDielectronVarManager::kY);
    
    //opening angle and PhiV
    histos->UserHistogram("Pair","InvMass_OpeningAngle",";Inv. Mass [GeV];Opening Angle;#pairs",
                          GetVector(kMee), GetVector(kOpAng), 
                          AliDielectronVarManager::kM, AliDielectronVarManager::kOpeningAngle);
    histos->UserHistogram("Pair","InvMass_PhivPair",";Inv. Mass [GeV];PhiV;#pairs",
                          GetVector(kMee), GetVector(kPhiV), 
                          AliDielectronVarManager::kM, AliDielectronVarManager::kPhivPair);
    histos->UserHistogram("Pair","PairPt_OpeningAngle",";Pair Pt [GeV];Opening Angle;#pairs",
                          GetVector(kPtee), GetVector(kOpAng), 
                          AliDielectronVarManager::kPt, AliDielectronVarManager::kOpeningAngle);
    histos->UserHistogram("Pair","PairPt_PhivPair",";Pair Pt [GeV];PhiV;#pairs",
                          GetVector(kPtee), GetVector(kPhiV), 
                          AliDielectronVarManager::kPt, AliDielectronVarManager::kPhivPair);
    histos->UserHistogram("Pair","OpeningAngle_PhivPair",";Opening Angle;PhiV;#pairs",
                          GetVector(kOpAng), GetVector(kPhiV), 
                          AliDielectronVarManager::kOpeningAngle, AliDielectronVarManager::kPhivPair);
    
  }// end: (!isQAtask)
  
  die->SetHistogramManager(histos);
}


//______________________________________________________________________________________
TVectorD *GetVector(Int_t var) 
{
  switch (var) 
  {
    case kPhiV:   return AliDielectronHelper::MakeLinBinning(100, 0., TMath::Pi());
    case kOpAng:  return AliDielectronHelper::MakeLinBinning(100, 0., TMath::Pi());
    case kOpAng2: return AliDielectronHelper::MakeLinBinning( 50, 0., TMath::Pi()/2.);
    case kEta2D:  return AliDielectronHelper::MakeLinBinning(100,-1,1);
    case kEta3D:  return AliDielectronHelper::MakeLinBinning( 50,-1,1);
    case kPhi2D:  return AliDielectronHelper::MakeLinBinning(120, 0., TMath::TwoPi());
    case kY3D:    return AliDielectronHelper::MakeLinBinning( 50,-2.5,2.5);
      
    case kSigmaEle:
      if (isQAtask) return AliDielectronHelper::MakeLinBinning(100,-10.,10.);
      else          return AliDielectronHelper::MakeLinBinning( 50, -5., 5.);
    case kSigmaOther:
      if (isQAtask) return AliDielectronHelper::MakeLinBinning(100,-20.,20.);
      else          return AliDielectronHelper::MakeLinBinning( 50,-10.,10.);
    case kTPCdEdx:
      if (isQAtask) return AliDielectronHelper::MakeLinBinning(120,  0.,120.);
      else          return AliDielectronHelper::MakeLinBinning( 50, 50.,100.);
      
    case kMee:    return AliDielectronHelper::MakeArbitraryBinning("0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 
                                                                   0.10, 0.14, 0.18, 0.22, 0.26, 0.30, 0.34, 0.38, 0.42, 0.46, 
                                                                   0.50, 0.54, 0.58, 0.62, 0.66, 0.70, 0.74, 0.78, 0.82, 0.86, 
                                                                   0.90, 0.94, 0.98, 1.02, 1.06, 
                                                                   1.10, 1.30, 1.50, 1.70, 1.90, 2.10, 2.30, 2.50, 2.70, 2.90, 
                                                                   3.10, 3.30, 3.50, 4.00, 4.50, 5.00 
                                                                   ");
    case kMee500: return AliDielectronHelper::MakeArbitraryBinning("0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 
                                                                   0.10, 0.14, 0.18, 0.22, 0.26, 0.30, 0.34, 0.38, 0.42, 0.46, 
                                                                   0.50 
                                                                   ");
    case kPtee:   return AliDielectronHelper::MakeArbitraryBinning("0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 
                                                                   0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 
                                                                   1.00, 1.20, 1.40, 1.60, 1.80, 2.00, 2.20, 2.40, 2.60, 2.80, 
                                                                   3.00, 3.20, 3.40, 3.60, 3.80, 4.00, 4.20, 4.40, 4.60, 4.80, 
                                                                   5.00, 6.00, 7.00, 8.00 
                                                                   ");
    case kP2D:    return AliDielectronHelper::MakeArbitraryBinning("0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 
                                                                   0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 
                                                                   1.00, 1.05, 1.10, 1.15, 1.20, 1.25, 1.30, 1.35, 1.40, 1.45, 
                                                                   1.50, 1.55, 1.60, 1.65, 1.70, 1.75, 1.80, 1.85, 1.90, 1.95, 
                                                                   2.00, 2.20, 2.40, 2.60, 2.80, 3.00, 3.20, 3.40, 3.60, 3.80, 
                                                                   4.00, 4.50, 5.00, 5.50, 6.00, 6.50, 7.00, 7.50, 8.00 
                                                                   ");
      
      
    default: cout << "ERROR: in 'GetVector(...var)' variable for axis range not defined!" << endl;
      break;
  } 
}


//______________________________________________________________________________________
TVectorD *BinsToVector(Int_t nbins, Double_t min, Double_t max) {
  return AliDielectronHelper::MakeLinBinning(nbins,min,max);
  //  TVectorD *vec = new TVectorD(nbins+1);
  //
  //  Double_t binwdth = (max-min)/nbins;
  //  for (int i = 0; i < nbins+1; i++) (*vec)[i] = min + i*binwdth;
  //  
  //  return vec;
}



void InitCF(AliDielectron* die, Int_t cutDefinition)
{
  //
  // Setupd the CF Manager if needed
  //
  AliDielectronCF *cf=new AliDielectronCF(die->GetName(),die->GetTitle());
 
  //pair variables
  cf->AddVariable(AliDielectronVarManager::kPt,100,0.,5.);
  //cf->AddVariable(AliDielectronVarManager::kP,200,0,20);
  //cf->AddVariable(AliDielectronVarManager::kPhi,64, -3.2, 3.2);
  cf->AddVariable(AliDielectronVarManager::kY,20,-1.,1.);
  cf->AddVariable(AliDielectronVarManager::kM,500,0.,4.); 
  //cf->AddVariable(AliDielectronVarManager::kPairType,10,0,10);
  //cf->AddVariable(AliDielectronVarManager::kOpeningAngle,315,0,3.15);
  //cf->AddVariable(AliDielectronVarManager::kDeltaEta,200,-2,2);
  //cf->AddVariable(AliDielectronVarManager::kDeltaPhi,100,0,3.15);
  //cf->AddVariable(AliDielectronVarManager::kHaveSameMother,21,-10,10);
  //cf->AddVariable(AliDielectronVarManager::kNumberOfDaughters,5,0,5);
  //leg variables
  cf->AddVariable(AliDielectronVarManager::kPt,100,0.,5.,kTRUE);
  //cf->AddVariable(AliDielectronVarManager::kP,200,0.,20.,kTRUE);
  cf->AddVariable(AliDielectronVarManager::kY,40,-2.,2.,kTRUE);
  //cf->AddVariable(AliDielectronVarManager::kEta,20,-1.,1.,kTRUE);
  cf->AddVariable(AliDielectronVarManager::kPhi,100,0.,3.15,kTRUE);
//  cf->AddVariable(AliDielectronVarManager::kITSsignal,1000,0.0.,1000.,kTRUE);
//  cf->AddVariable(AliDielectronVarManager::kHaveSameMother,21,-10,10,kTRUE);
//  cf->AddVariable(AliDielectronVarManager::kNumberOfDaughters,11,0,10,kTRUE);
  
  cf->SetStepForMCtruth();
  //cf->SetStepsForEachCut();
  //cf->SetStepForPreFilter();
  cf->SetStepForAfterAllCuts();
  //cf->SetStepsForBackground();
  cf->SetStepsForSignal();
  
  die->SetCFManagerPair(cf);

}

void AddMCSignals(AliDielectron* die){
  
  /*
  AliDielectronSignalMC* lowMassDiele = new AliDielectronSignalMC("lowMassDiele","low mass dielectron pairs");
  lowMassDiele->SetLegPDGs(11,-11);
  lowMassDiele->SetCheckBothChargesLegs(kTRUE,kTRUE);
  lowMassDiele->SetLegSources(AliDielectronSignalMC::kPrimary,AliDielectronSignalMC::kPrimary);
  lowMassDiele->SetFillPureMCStep(kTRUE);
  die->AddSignalMC(lowMassDiele);
  */
  
  //AliDielectronSignalMC* secondary = new AliDielectronSignalMC("secondary","secondary electrons pairs");
  //secondary->SetLegPDGs(11,-11);
  //secondary->SetCheckBothChargesLegs(kTRUE,kTRUE);
  //secondary->SetLegSources(AliDielectronSignalMC::kSecondary,AliDielectronSignalMC::kSecondary);
  //die->AddSignalMC(secondary);
  
  //AliDielectronSignalMC* eleFromConversions = new AliDielectronSignalMC("eleFromConversions","conversion electrons");
  //eleFromConversions->SetLegPDGs(11,-11);
  //eleFromConversions->SetCheckBothChargesLegs(kTRUE,kTRUE);
  //eleFromConversions->SetLegSources(AliDielectronSignalMC::kSecondary,AliDielectronSignalMC::kSecondary);
  //eleFromConversions->SetMotherPDGs(22,22);    // 22- photon
  //die->AddSignalMC(eleFromConversions);
  
  //AliDielectronSignalMC* misIdPions = new AliDielectronSignalMC("misIdPions","mis id. pion pairs");
  //misIdPions->SetLegPDGs(211,-211);
  //misIdPions->SetCheckBothChargesLegs(kTRUE,kTRUE);
  //misIdPions->SetMotherPDGs(0,0);
  //misIdPions->SetLegSources(AliDielectronSignalMC::kDontCare,AliDielectronSignalMC::kDontCare);
  //misIdPions->SetMothersRelation(AliDielectronSignalMC::kUndefined);
  //die->AddSignalMC(misIdPions);
  
  //AliDielectronSignalMC* misIdPionEle = new AliDielectronSignalMC("misIdPionEle","mis id. pion ele");
  //misIdPionEle->SetLegPDGs(11,-211);
  //misIdPionEle->SetCheckBothChargesLegs(kTRUE,kTRUE);
  //misIdPionEle->SetMotherPDGs(0,0);
  //misIdPionEle->SetLegSources(AliDielectronSignalMC::kDontCare,AliDielectronSignalMC::kDontCare);
  //misIdPionEle->SetMothersRelation(AliDielectronSignalMC::kUndefined);
  //die->AddSignalMC(misIdPionEle);
 
  //AliDielectronSignalMC* misIdElePion = new AliDielectronSignalMC("misIdElePion","mis id. ele pion");
  //misIdElePion->SetLegPDGs(211,-11);
  //misIdElePion->SetCheckBothChargesLegs(kTRUE,kTRUE);
  //misIdElePion->SetMotherPDGs(0,0);
  //misIdElePion->SetLegSources(AliDielectronSignalMC::kDontCare,AliDielectronSignalMC::kDontCare);
  //misIdElePion->SetMothersRelation(AliDielectronSignalMC::kUndefined);
  //die->AddSignalMC(misIdElePion);
  
  //AliDielectronSignalMC* misIdKaons = new AliDielectronSignalMC("misIdKaons","mis id. kaon pairs");
  //misIdKaons->SetLegPDGs(321,-321);
  //misIdKaons->SetCheckBothChargesLegs(kTRUE,kTRUE);
  //misIdKaons->SetMotherPDGs(0,0);
  //misIdKaons->SetLegSources(AliDielectronSignalMC::kDontCare,AliDielectronSignalMC::kDontCare);
  //misIdKaons->SetMothersRelation(AliDielectronSignalMC::kUndefined);
  //die->AddSignalMC(misIdKaons);
  
  //AliDielectronSignalMC* misIdKaonEle = new AliDielectronSignalMC("misIdKaonEle","mis id. kaon ele");
  //misIdKaonEle->SetLegPDGs(11,-321);
  //misIdKaonEle->SetCheckBothChargesLegs(kTRUE,kTRUE);
  //misIdKaonEle->SetMotherPDGs(0,0);
  //misIdKaonEle->SetLegSources(AliDielectronSignalMC::kDontCare,AliDielectronSignalMC::kDontCare);
  //misIdKaonEle->SetMothersRelation(AliDielectronSignalMC::kUndefined);
  //die->AddSignalMC(misIdKaonEle);

  //AliDielectronSignalMC* misIdEleKaon = new AliDielectronSignalMC("misIdEleKaon","mis id. ele kaon");
  //misIdEleKaon->SetLegPDGs(321,-11);
  //misIdEleKaon->SetCheckBothChargesLegs(kTRUE,kTRUE);
  //misIdEleKaon->SetMotherPDGs(0,0);
  //misIdEleKaon->SetLegSources(AliDielectronSignalMC::kDontCare,AliDielectronSignalMC::kDontCare);
  //misIdEleKaon->SetMothersRelation(AliDielectronSignalMC::kUndefined);
  //die->AddSignalMC(misIdEleKaon);

  //AliDielectronSignalMC* misIdProtons = new AliDielectronSignalMC("misIdProtons","mis id. proton pairs");
  //misIdProtons->SetLegPDGs(321,-321);
  //misIdProtons->SetCheckBothChargesLegs(kTRUE,kTRUE);
  //misIdProtons->SetMotherPDGs(0,0);
  //misIdProtons->SetLegSources(AliDielectronSignalMC::kDontCare,AliDielectronSignalMC::kDontCare);
  //misIdProtons->SetMothersRelation(AliDielectronSignalMC::kUndefined);
  //die->AddSignalMC(misIdProtons);
   
  //AliDielectronSignalMC* misIdProtonEle = new AliDielectronSignalMC("misIdProtonEle","mis id. proton ele");
  //misIdProtonEle->SetLegPDGs(11,-2212);
  //misIdProtonEle->SetCheckBothChargesLegs(kTRUE,kTRUE);
  //misIdProtonEle->SetMotherPDGs(0,0);
  //misIdProtonEle->SetLegSources(AliDielectronSignalMC::kDontCare,AliDielectronSignalMC::kDontCare);
  //misIdProtonEle->SetMothersRelation(AliDielectronSignalMC::kUndefined);
  //die->AddSignalMC(misIdProtonEle);

  //AliDielectronSignalMC* misIdEleProton = new  AliDielectronSignalMC("misIdEleProton","mis id. ele proton");
  //misIdEleProton->SetLegPDGs(2212,-11);
  //misIdEleProton->SetCheckBothChargesLegs(kTRUE,kTRUE);
  //misIdEleProton->SetMotherPDGs(0,0);
  //misIdEleProton->SetLegSources(AliDielectronSignalMC::kDontCare,AliDielectronSignalMC::kDontCare);
  //misIdEleProton->SetMothersRelation(AliDielectronSignalMC::kUndefined);
  //die->AddSignalMC(misIdEleProton);   
  
  //AliDielectronSignalMC* dalitzDecays = new AliDielectronSignalMC("dalitzDecays","dalitz Pairs");
  //dalitzDecays->SetLegPDGs(11,-11);
  //dalitzDecays->SetCheckBothChargesLegs(kTRUE,kTRUE);
  //dalitzDecays->SetLegSources(AliDielectronSignalMC::kSecondary,AliDielectronSignalMC::kSecondary);
  //dalitzDecays->SetMotherPDGs(111,111);
  //dalitzDecays->SetFillPureMCStep(kTRUE);
  //die->AddSignalMC(dalitzDecays);

  
    AliDielectronSignalMC* PhiDecays= new AliDielectronSignalMC("PhiDecays","Phi Pairs");
    PhiDecays->SetLegPDGs(11,-11);
    PhiDecays->SetCheckBothChargesLegs(kTRUE,kTRUE);
    PhiDecays->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
    PhiDecays->SetMotherPDGs(333,333);
    PhiDecays->SetMothersRelation(AliDielectronSignalMC::kSame); 
    PhiDecays->SetFillPureMCStep(kTRUE);
    die->AddSignalMC(PhiDecays);
  
  
    AliDielectronSignalMC* OmegaDecays= new AliDielectronSignalMC("OmegaDecays","Omega Pairs");
    OmegaDecays->SetLegPDGs(11,-11);
    OmegaDecays->SetCheckBothChargesLegs(kTRUE,kTRUE);
    OmegaDecays->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
    OmegaDecays->SetMotherPDGs(223,223);
    OmegaDecays->SetMothersRelation(AliDielectronSignalMC::kSame);
    OmegaDecays->SetDalitz(AliDielectronSignalMC::kIsNotDalitz); 
    OmegaDecays->SetFillPureMCStep(kTRUE);
    die->AddSignalMC(OmegaDecays);
  
    /*
    AliDielectronSignalMC* RhoDecays= new AliDielectronSignalMC("RhoDecays","Rho Pairs");
    RhoDecays->SetLegPDGs(11,-11);
    RhoDecays->SetCheckBothChargesLegs(kTRUE,kTRUE);
    RhoDecays->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
    RhoDecays->SetMotherPDGs(113,113);
    RhoDecays->SetMothersRelation(AliDielectronSignalMC::kSame); 
    RhoDecays->SetFillPureMCStep(kTRUE);
    die->AddSignalMC(RhoDecays);
    */

  //AliDielectronSignalMC* charm = new AliDielectronSignalMC("charm","charm electrons pairs");
  //charm->SetLegPDGs(11,-11);
  //charm->SetMotherPDGs(402,402);
  //charm->SetCheckBothChargesLegs(kTRUE,kTRUE);
  //charm->SetCheckBothChargesMothers(kTRUE,kTRUE);
  //charm->SetLegSources(AliDielectronSignalMC::kFinalState,AliDielectronSignalMC::kFinalState);
  //charm->SetMotherPDGs(402,402);  
  //die->AddSignalMC(charm);

  //AliDielectronSignalMC* DieleConti= new AliDielectronSignalMC("DieleConti","low mass ee pairs");
  //DieleConti->SetLegPDGs(11,-11);
  //DieleConti->SetMotherPDGs(0,0,22,22);
  //DieleConti->SetCheckBothChargesLegs(kTRUE,kTRUE);
  //DieleConti->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  //DieleConti->SetMothersRelation(AliDielectronSignalMC::kSame);
  //DieleConti->SetFillPureMCStep(kTRUE);
  //die->AddSignalMC(DieleConti);
  

}


void InitHF(AliDielectron* die, Int_t cutDefinition)
{
  //
  // Setup the HF arrays
  //

  AliDielectronHF *hf=new AliDielectronHF(die->GetName(),die->GetTitle());
  hf->SetPairTypes(AliDielectronHF::kSeMeAll);
  //hf->SetStepForMCGenerated(); 
  hf->AddCutVariable(AliDielectronVarManager::kM, 100,0.0,4.0);
  hf->UserHistogram("Pair", AliDielectronHelper::MakeLinBinning(100,0.0,4.0), AliDielectronVarManager::kM);
  die->SetHistogramArray(hf);
}

AliDielectronEventCuts* GetEventCuts(){

  AliDielectronEventCuts *eventCuts=new AliDielectronEventCuts("eventCuts","Vertex Track && |vtxZ|<10 && ncontrib>0");
  eventCuts->SetRequireVertex();
  eventCuts->SetVertexZ(-10.,10.);
  eventCuts->SetMinVtxContributors(1); 
  
  return eventCuts;
}
