void InitHistograms(AliDielectron *die, Int_t cutDefinition=0);
void SetupCuts(AliDielectron *die, Int_t cutDefinition, TString userPathWeightFile, Double_t userTMVACutValue);
void SetTPCCorr(AliDielectron *die);
const AliDielectronEventCuts *GetEventCuts();
void SetupMCsignals(AliDielectron* die);
TVectorD* getBinVec(int val);
AliDielectronCutGroup* SetupTrackCutsAndSettings(Int_t selTr, Int_t selPID, Int_t MVACut=0., Bool_t useAODFilterCuts, TString TMVAweight);
        
Bool_t isRandomRejTask=kFALSE;//needed for InitHistograms() //dont change!!!
Bool_t kRot = kFALSE;
Bool_t kMix = kTRUE;
Bool_t randomizeDau = kTRUE;
Bool_t bUsePileUpCutsTPCClusters=kTRUE;

enum {kMee = 0, kPtee, kPhi, kCent};

// available cut defintions
const Int_t nMax = 3; 
const Int_t nPF  = 999; // use prefiltering for cuts > nPF

void Config_slehner_diele_TMVA(AliAnalysisTaskMultiDielectron *task,Bool_t usePIDCorr=kFALSE,Bool_t hasMC=kFALSE, Bool_t useAODFilterCuts,  TString TMVAweight)
{
  // create the actual framework object
  Int_t trackCut=0;
  Int_t PIDCut=0;
  Int_t MVACut=0;
  Int_t pairCut=-1;
  
  for(int glcut = 0; glcut <=30; ++glcut){
//  for(int glcut = 0; glcut <=0; glcut++){
    ////////DEFINE THE CUTS AS FUNCTION OF GLCUT//////
    if(glcut>0 && glcut<21) continue;    
    if(glcut==0){
      trackCut=-1;
      PIDCut=0;
    }
    else{
      PIDCut=glcut-10;
      trackCut=glcut;
    }
    //event plane cuts for pairs
//    if(glcut==0) trackCut=30;    //in event plane
//    if(glcut==0) PIDCut=20;
//    if(glcut==0) pairCut=1;
//    if(glcut==1) trackCut=30;    //out event plane
//    if(glcut==1) PIDCut=20;
//    if(glcut==1) pairCut=0;    
//    if(glcut==2) trackCut=22;     //in of event plane
//    if(glcut==2) PIDCut=12;
//    if(glcut==2) pairCut=1;    
//    if(glcut==3) trackCut=22;     //out event plane
//    if(glcut==3) PIDCut=12;
//    if(glcut==3) pairCut=0;
    
    for(MVACut = 0; MVACut<10;MVACut++){
//      if(MVACut!=0 && MVACut!=3) continue;
      TString name=TString::Format("DieleTr%d_PID%d_MVA%d",trackCut,PIDCut, MVACut);
//      TString name=TString::Format("DieleTr%d_PID%d_Pair%d_MVA%d",trackCut,PIDCut, pairCut, MVACut);
      cout<<"Adding Diele Task: "<<name.Data()<<endl;    
      AliDielectron * diel_low = new AliDielectron(Form("%s",name.Data()), Form("Name: %s",name.Data()));
      if(!diel_low){
        Printf("=======================================");
        Printf("No AliDielectron object loaded -> EXIT ");
        Printf("=======================================");
        return NULL; 
      }  
      if(hasMC) SetupMCsignals(diel_low);
      if(kMix && !hasMC ){ // need second since there is a problem when mixing MC events (TRef?)
        AliDielectronMixingHandler *mix = new AliDielectronMixingHandler;

        mix->AddVariable(AliDielectronVarManager::kZvPrim,"-10,-5,0,5,10");
        mix->AddVariable(AliDielectronVarManager::kCentrality,"0,5,10,20,30,50,80");
        mix->SetDepth(15);
        mix->SetMixType(AliDielectronMixingHandler::kAll);
        diel_low->SetMixingHandler(mix);
      }

      if(usePIDCorr){
       SetITSCorr(diel_low,hasMC);
       SetTPCCorr(diel_low,hasMC);
//       SetTOFCorr(diel_low,hasMC);
      }

      diel_low->SetUseKF(kFALSE);   //keep this one, otherwise masses are slightly wrong and R factors very wrong!
      InitHistograms(diel_low,0,hasMC);

      if(bUsePileUpCutsTPCClusters){
        Double_t pileUpCutsTPCClustersMin = 2.0e+06;
        Double_t pileUpCutsTPCClustersMax = 3.1e+06;        
        
        printf("Adding TPC-Cluster based Pile-up rejection! \n");
        TF1* fFitMin = new TF1("fFit","pol6",0,90);
        fFitMin->SetParameters(pileUpCutsTPCClustersMin,-95678.946999,2152.010478,-50.119000,0.780528,-0.006150,0.000019);
        TF1* fFitMax = new TF1("fFit","pol6",0,90);
        fFitMax->SetParameters(pileUpCutsTPCClustersMax,-95678.946999,2152.010478,-50.119000,0.780528,-0.006150,0.000019);
               
        AliDielectronEventCuts *pileUpCuts = new AliDielectronEventCuts("pileUpCuts","pileUpCuts");
        pileUpCuts->SetMinCorrCutFunction(fFitMin, AliDielectronVarManager::kCentralityNew, AliDielectronVarManager::kNTPCclsEvent);
        pileUpCuts->SetMaxCorrCutFunction(fFitMax, AliDielectronVarManager::kCentralityNew, AliDielectronVarManager::kNTPCclsEvent);

        pileUpCuts->Print();
        diel_low->GetEventFilter().AddCuts(pileUpCuts);
    }      
      
      std::cout << "CutTr: "<<trackCut<<" CutPID: "<<PIDCut<<" MVAcut: "<<-1+MVACut*0.2<<" being added"<< std::endl;
      diel_low->GetTrackFilter().AddCuts(SetupTrackCutsAndSettings(trackCut, PIDCut, MVACut, useAODFilterCuts,TMVAweight));   
      if(pairCut!=-1) diel_low->GetPairFilter().AddCuts(SetupPairCuts(pairCut));   
      
      task->AddDielectron(diel_low);
      printf("successfully added AliDielectron: %s\n",diel_low->GetName());           
      }

  }
 return;

}

AliDielectronCutGroup* SetupTrackCutsAndSettings(Int_t selTr, Int_t selPID, Int_t MVACut, Bool_t useAODFilterCuts,TString TMVAweight)
{
  std::cout<<"SetupTrackCutsAndSettings: "<<selTr<<","<<selPID<<std::endl;
//  AliAnalysisFilter *anaFilter = new AliAnalysisFilter("anaFilter","anaFilter"); // named constructor seems mandatory!

  LMEECutLib* LMcutlib = new LMEECutLib();

  std::cout << "Get CutTr: "<<selTr<<" CutPID: "<<selPID<<std::endl;
  return (LMcutlib->GetTrackCuts(selTr, selPID, MVACut ,useAODFilterCuts, TMVAweight));
}

AliDielectronCutGroup* SetupPairCuts(Int_t sel)
{
  std::cout<<"SetupPairCuts for Pair Cut Nr: "<<sel<<std::endl;
//  AliAnalysisFilter *anaFilter = new AliAnalysisFilter("anaFilter","anaFilter"); // named constructor seems mandatory!

  LMEECutLib* LMcutlib = new LMEECutLib();

  return (LMcutlib->GetPairCuts(sel));
}


//______________________________________________________________________________________
void SetTPCCorr(AliDielectron *die, Bool_t hasMC){
  ::Info("Config_slehner_LMEE_TMVA","starting LMEECutLib::SetEtaCorrection for TPC\n");
  TString path="alien:///alice/cern.ch/user/s/selehner/recal/";
  if(hasMC) TString fName= "recalib_mc_tpc_nsigmaele.root";
  else      TString fName= "recalib_data_tpc_nsigmaele.root";
  TFile* _file = TFile::Open(fName.Data());
//  _file = TFile::Open(fName.Data());
  if(!_file){
    gSystem->Exec(TString::Format("alien_cp %s .",(path+fName).Data()));
    ::Info("Config_slehner_LMEE_TMVA","Get TPC correction from Alien: %s",(path+fName).Data());
    _file = TFile::Open(fName.Data());
    if(!_file)  ::Error("Config_slehner_LMEE_TMVA","Cannot get TPC correction from Alien: %s",(path+fName).Data());
  }
  
  TH3D* mean = dynamic_cast<TH3D*>(_file->Get("sum_mean_correction"));
  TH3D* width= dynamic_cast<TH3D*>(_file->Get("sum_width_correction"));

  if(mean)   ::Info("Config_slehner_LMEE_TMVA","Mean Correction Histo loaded, entries: %f",mean->GetEntries());
  else {
    ::Info("Config_slehner_LMEE_TMVA","Mean Correction Histo not loaded! entries: %f",mean->GetEntries());
    return 0;
  }
    die->SetCentroidCorrFunction(mean, AliDielectronVarManager::kP, AliDielectronVarManager::kEta, AliDielectronVarManager::kRefMultTPConly);
    die->SetWidthCorrFunction(width, AliDielectronVarManager::kP, AliDielectronVarManager::kEta, AliDielectronVarManager::kRefMultTPConly);     
}

//______________________________________________________________________________________
void SetITSCorr(AliDielectron *die, Bool_t hasMC){
  ::Info("Config_slehner_LMEE_TMVA","starting SetITSCorr\n");
  TString path="alien:///alice/cern.ch/user/s/selehner/recal/";
  if(hasMC) TString fName= "recalib_mc_its_nsigmaele.root";
  else      TString fName= "recalib_data_its_nsigmaele.root";
  gSystem->Exec(TString::Format("alien_cp %s .",(path+fName).Data()));
  ::Info("Config_slehner_LMEE_TMVA","Get ITS correction from Alien: %s",(path+fName).Data());
  _file = TFile::Open(fName.Data());
  
  TH3D* mean = dynamic_cast<TH3D*>(_file->Get("sum_mean_correction"));
  TH3D* width= dynamic_cast<TH3D*>(_file->Get("sum_width_correction"));

  if(mean)   ::Info("Config_slehner_LMEE_TMVA","Mean Correction Histo loaded, entries: %f",mean->GetEntries());
  else {
    ::Info("Config_slehner_LMEE_TMVA","Mean Correction Histo not loaded! entries: %f",mean->GetEntries());
    return 0;
  }
    die->SetCentroidCorrFunctionITS(mean, AliDielectronVarManager::kP, AliDielectronVarManager::kEta, AliDielectronVarManager::kRefMultTPConly);
    die->SetWidthCorrFunctionITS(width, AliDielectronVarManager::kP, AliDielectronVarManager::kEta, AliDielectronVarManager::kRefMultTPConly);     
}

//______________________________________________________________________________________
void SetTOFCorr(AliDielectron *die, Bool_t hasMC){
  ::Info("Config_slehner_LMEE_TMVA","starting SetTOFCorr for TOF\n");
  TString path="alien:///alice/cern.ch/user/s/selehner/recal/";
  if(hasMC) TString fName= "recalib_mc_tof_nsigmaele.root";
  else      TString fName= "recalib_data_tof_nsigmaele.root";
  gSystem->Exec(TString::Format("alien_cp %s .",(path+fName).Data()));
  ::Info("Config_slehner_LMEE_TMVA","Get TOF correction from Alien: %s",(path+fName).Data());
  _file = TFile::Open(fName.Data());
  
  TH3D* mean = dynamic_cast<TH3D*>(_file->Get("sum_mean_correction"));
  TH3D* width= dynamic_cast<TH3D*>(_file->Get("sum_width_correction"));

  if(mean)   ::Info("Config_slehner_LMEE_TMVA","Mean Correction Histo loaded, entries: %f",mean->GetEntries());
  else {
    ::Info("Config_slehner_LMEE_TMVA","Mean Correction Histo not loaded! entries: %f",mean->GetEntries());
    return 0;
  }
    die->SetCentroidCorrFunctionTOF(mean, AliDielectronVarManager::kP, AliDielectronVarManager::kEta, AliDielectronVarManager::kRefMultTPConly);
    die->SetWidthCorrFunctionTOF(width, AliDielectronVarManager::kP, AliDielectronVarManager::kEta, AliDielectronVarManager::kRefMultTPConly);     
}

//______________________________________________________________________________________
void InitHistograms(AliDielectron *die, Int_t cutDefinition, Bool_t isMC)
{
  
  //Setup histogram classes
  AliDielectronHistos *histos= new AliDielectronHistos(die->GetName(), die->GetTitle());
  histos->SetReservedWords("Track;Pair;Track_Legs;Pre;RejPair;RejTrack;Random");

  //Event class
  histos->AddClass("Event");

  TH1::AddDirectory(kFALSE);
  
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
//  if(isRandomRejTask){
//    //
//    // _____ histograms for AliAnalysisTaskMultiDielectronPR _____
//    //
//    //    histos->AddClass("Rand_Pair");
//    //    histos->AddClass("Rand_RejPair");
//    const char* cRandomPairClassNames[2] = { "Testpart", "RejTestpart" };
//    for (Int_t i=0; i<2; ++i){
//      histos->AddClass(Form("Random_%s",cRandomPairClassNames[i]));
//    }
//    histos->UserHistogram("Random","Pt","",200,0,10.,AliDielectronVarManager::kPt);
//    histos->UserHistogram("Random","Eta","",200,-2,2,AliDielectronVarManager::kEta);
//    histos->UserHistogram("Random","Phi","",120,0.,TMath::TwoPi(),AliDielectronVarManager::kPhi);
//    histos->UserHistogram("Random","Px","",200,0,10.,AliDielectronVarManager::kPx);
//    histos->UserHistogram("Random","Py","",200,0,10.,AliDielectronVarManager::kPy);
//    histos->UserHistogram("Random","Pz","",200,0,10.,AliDielectronVarManager::kPz);
//    histos->UserHistogram("Random","Pt_Eta_Phi","",
//                          500,0.,10.,16,-0.8,0.8,30,0.,2*TMath::Pi(),
//                          AliDielectronVarManager::kPt,AliDielectronVarManager::kEta,AliDielectronVarManager::kPhi);
//  }

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
//  histos->UserHistogram("Event","ZVertex","ZVertex;ZVertex/cm",120,-12.,12.,AliDielectronVarManager::kZvPrim);
  histos->UserHistogram("Event","Centrality","Centrality;Centrality/%",202,-1.,100.,AliDielectronVarManager::kCentralityNew);
//  histos->UserHistogram("Event","nEvTPC",";;ev plane;",AliDielectronHelper::MakeLinBinning(180,  -TMath::Pi()/2.,TMath::Pi()/2.),AliDielectronVarManager::kQnTPCrpH2);
//  histos->UserHistogram("Event","TPCrpH2:Centrality", "", 180, -TMath::Pi()/2.,TMath::Pi()/2.,202,-1.,100.,AliDielectronVarManager::kQnTPCrpH2,AliDielectronVarManager::kCentralityNew);


  //add histograms to track class

//  histos->UserHistogram("Track","nITS","nITS;nITS;nTracks",10,0,10.,AliDielectronVarManager::kNclsITS); 
//  histos->UserHistogram("Track", "ITS1Shared", "ITS1Shared", 10,0,10., AliDielectronVarManager::kClsS1ITS);
//  histos->UserHistogram("Track", "ITS2Shared", "ITS2Shared", 10,0,10., AliDielectronVarManager::kClsS2ITS);
//  histos->UserHistogram("Track", "ITS3Shared", "ITS3Shared", 10,0,10., AliDielectronVarManager::kClsS3ITS);
//  histos->UserHistogram("Track", "ITS4Shared", "ITS4Shared", 10,0,10., AliDielectronVarManager::kClsS4ITS);
//  histos->UserHistogram("Track", "ITS5Shared", "ITS5Shared", 10,0,10., AliDielectronVarManager::kClsS5ITS);
//  histos->UserHistogram("Track", "ITS6Shared", "ITS6Shared", 10,0,10., AliDielectronVarManager::kClsS6ITS);
//  histos->UserHistogram("Track","nITSshared_frac", "nITSshared_frac",100,0,1.1,AliDielectronVarManager::kNclsSFracITS);
//  histos->UserHistogram("Track","NCrossedRowsTPC","",200,0,200,  AliDielectronVarManager::kNclsCrTPC);
//  histos->UserHistogram("Track","NClustersTPC", "",200,0,200, AliDielectronVarManager::kNclsTPC);
//  histos->UserHistogram("Track","NTPCSignal","" ,200,0,200,  AliDielectronVarManager::kTPCsignalN);
//  histos->UserHistogram("Track","log(abs(DCAxy))","", 100,-20,2,  AliDielectronVarManager::kLogDCAXY);
//  histos->UserHistogram("Track","log(abs(DCAz))" ,"", 100,-20,2,  AliDielectronVarManager::kLogDCAZ);  
//  histos->UserHistogram("Track","chi2GlobalPerNDF","", 100,0,10,  AliDielectronVarManager::kChi2GlobalNDF);
//  histos->UserHistogram("Track","chi2ITS","" , 100,0,100,  AliDielectronVarManager::kITSchi2);
//  histos->UserHistogram("Track","eta","" , 100,-0.8,0.8, AliDielectronVarManager::kEta);
//  histos->UserHistogram("Track","phi","" , 100,0,7, AliDielectronVarManager::kPhi);
  histos->UserHistogram("Track","pt", "", 100,0,10,AliDielectronVarManager::kPt);  
  if(isMC)  histos->UserHistogram("Track","PdgCodeMother", "", 1,0,0,AliDielectronVarManager::kPdgCodeMother);  
  histos->UserHistogram("Track","nSigmaITSEl:pt", "", 100,0,10,100,-5,5,AliDielectronVarManager::kPt,AliDielectronVarManager::kITSnSigmaEle);  
  histos->UserHistogram("Track","nSigmaTPCEl:pt", "", 100,0,10,100,-5,5,AliDielectronVarManager::kPt,AliDielectronVarManager::kTPCnSigmaEle);  
  histos->UserHistogram("Track","nSigmaTOFEl:pt", "", 100,0,10,100,-5,5,AliDielectronVarManager::kPt,AliDielectronVarManager::kTOFnSigmaEle);  
  histos->UserHistogram("Track","nSigmaITSEl:eta", "", 100,-1,1,100,-5,5,AliDielectronVarManager::kEta,AliDielectronVarManager::kITSnSigmaEle);  
  histos->UserHistogram("Track","nSigmaTPCEl:eta", "", 100,-1,1,100,-5,5,AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEle);  
  histos->UserHistogram("Track","nSigmaTOFEl:eta", "", 100,-1,1,100,-5,5,AliDielectronVarManager::kEta,AliDielectronVarManager::kTOFnSigmaEle);  
//  histos->UserHistogram("Track","DeltaPhiTrackTPCrpH2", "", 100,-TMath::Pi(),TMath::Pi(),AliDielectronVarManager::kQnDeltaPhiTrackTPCrpH2);  
//  histos->UserHistogram("Track","Phi", "", 100,0,2*TMath::Pi(),AliDielectronVarManager::kPhi);  
//  histos->UserHistogram("Track","TPCrpH2", "", 100,-TMath::Pi()/2.,TMath::Pi()/2.,AliDielectronVarManager::kQnTPCrpH2);  
 
//lmee mass spectrum
//  TVectorD* mbins=  AliDielectronHelper::MakeArbitraryBinning(" 0.00, 0.02 ,0.04 ,0.08 ,0.14 ,0.22 ,0.38 ,0.54 ,1.1 ,1.7 ,2.5 ,2.9 ,3.0 ,3.1 ,3.3 ,3.5 ,4.0 ,5.0"); //Carsten's binning
//  TVectorD* ptbins= AliDielectronHelper::MakeArbitraryBinning("0.0,0.4,0.6,1,2.5,8");
//  histos->UserHistogram("Pair","InvMass_pPt","Inv.Mass:PairPt;Inv. Mass (GeV/c^{2});Pair Pt (GeV/c); ",
//                        mbins, ptbins,
//                        AliDielectronVarManager::kM, AliDielectronVarManager::kPt);
  
////low ptee
  TVectorD* mbins=  AliDielectronHelper::MakeArbitraryBinning(" 0.0,0.1,0.4,0.5 ,0.6 , 0.7 ,1.1, 1.5, 2.0 ,2.7 , 2.9, 3.1 , 5.0"); // for low ptee
  TVectorD* ptbins= AliDielectronHelper::MakeArbitraryBinning("0.0, 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 1, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0");
  TVectorD* centbins= AliDielectronHelper::MakeLinBinning(20,0,100);
  TVectorD* phibins= AliDielectronHelper::MakeLinBinning(20,0,TMath::Pi());
  
  histos->UserHistogram("Pair","InvMass_pPt_cent",";Inv. Mass (GeV/c^{2});Pair Pt (GeV/c); Centrality (V0M)",
                        getBinVec(kMee), getBinVec(kPtee), getBinVec(kCent),
                        AliDielectronVarManager::kM, AliDielectronVarManager::kPt, AliDielectronVarManager::kCentrality);

  histos->UserHistogram("Pair","InvMass_pPt_DeltaPhiSumDiff",";Inv. Mass (GeV/c^{2});Pair Pt (GeV/c);DeltaPhiSumDiff (rad)",
                        getBinVec(kMee), getBinVec(kPtee), getBinVec(kPhi),
                        AliDielectronVarManager::kM, AliDielectronVarManager::kPt, AliDielectronVarManager::kDeltaPhiSumDiff);

  histos->UserHistogram("Pair","InvMass_pPt_DeltaPhiSumPos",";Inv. Mass (GeV/c^{2});Pair Pt (GeV/c);DeltaPhiSumPos (rad)",
                        getBinVec(kMee), getBinVec(kPtee), getBinVec(kPhi),
                        AliDielectronVarManager::kM, AliDielectronVarManager::kPt, AliDielectronVarManager::kDeltaPhiSumPos);

  histos->UserHistogram("Pair","InvMass_pPt_DeltaPhiSumNeg",";Inv. Mass (GeV/c^{2});Pair Pt (GeV/c);DeltaPhiSumNeg (rad)",
                        getBinVec(kMee), getBinVec(kPtee), getBinVec(kPhi),
                        AliDielectronVarManager::kM, AliDielectronVarManager::kPt, AliDielectronVarManager::kDeltaPhiSumNeg);

  histos->UserHistogram("Pair","DeltaPhiSumDiff_DeltaPhiSumPos",";DeltaPhiSumDiff (rad);DeltaPhiSumPos (rad)",
                        getBinVec(kPhi), getBinVec(kPhi),
                        AliDielectronVarManager::kDeltaPhiSumDiff, AliDielectronVarManager::kDeltaPhiSumPos);

  histos->UserHistogram("Pair","DeltaPhiSumDiff_DeltaPhiSumNeg",";DeltaPhiSumDiff (rad);DeltaPhiSumNeg (rad)",
                        getBinVec(kPhi), getBinVec(kPhi),
                        AliDielectronVarManager::kDeltaPhiSumDiff, AliDielectronVarManager::kDeltaPhiSumNeg);

  histos->UserHistogram("Pair","DeltaPhiSumPos_DeltaPhiSumNeg",";DeltaPhiSumPos (rad);DeltaPhiSumNeg (rad)",
                        getBinVec(kPhi), getBinVec(kPhi),
                        AliDielectronVarManager::kDeltaPhiSumPos, AliDielectronVarManager::kDeltaPhiSumNeg);


  
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

//  // ##################### "real" pairs from signals (pi0,eta,eta',rho, omega, phi) ##############################
//   AliDielectronSignalMC* pi0Sig = new AliDielectronSignalMC("pi0", "pi0Signal"); ///pi0 dalitz pairs 
//  pi0Sig->SetLegPDGs(11,-11);
//  pi0Sig->SetMotherPDGs(111,111);
//  pi0Sig->SetMothersRelation(AliDielectronSignalMC::kSame);
//  pi0Sig->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
//  pi0Sig->SetMotherSources(AliDielectronSignalMC::kPrimary, AliDielectronSignalMC::kPrimary);
//  pi0Sig->SetCheckBothChargesLegs(kTRUE,kTRUE);
//  pi0Sig->SetCheckBothChargesMothers(kTRUE,kTRUE);
//  pi0Sig->SetFillPureMCStep(kTRUE);
//  die->AddSignalMC(pi0Sig);
//
//  AliDielectronSignalMC* pi0All = new AliDielectronSignalMC("pi0", "pi0All"); ///pi0 dalitz pairs (also from secondary)
//  pi0All->SetLegPDGs(11,-11);
//  pi0All->SetMotherPDGs(111,111);
//  pi0All->SetMothersRelation(AliDielectronSignalMC::kSame);
//  pi0All->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
//  pi0All->SetCheckBothChargesLegs(kTRUE,kTRUE);
//  pi0All->SetCheckBothChargesMothers(kTRUE,kTRUE);
//  pi0All->SetFillPureMCStep(kTRUE);
//  die->AddSignalMC(pi0All);
//
//
//  AliDielectronSignalMC* etaSig = new AliDielectronSignalMC("Eta", "etaSignal"); ///eta dalitz pairs 
//  etaSig->SetLegPDGs(11,-11);
//  etaSig->SetMotherPDGs(221,221);
//  etaSig->SetMothersRelation(AliDielectronSignalMC::kSame);
//  etaSig->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
//  etaSig->SetMotherSources(AliDielectronSignalMC::kPrimary, AliDielectronSignalMC::kPrimary);
//  etaSig->SetCheckBothChargesLegs(kTRUE,kTRUE);
//  etaSig->SetCheckBothChargesMothers(kTRUE,kTRUE);
//  etaSig->SetFillPureMCStep(kTRUE);
//  die->AddSignalMC(etaSig);
//
//
//  AliDielectronSignalMC* etaprimeSig = new AliDielectronSignalMC("Etaprime", "etaprimeSignal"); ///etaprime pairs 
//  etaprimeSig->SetLegPDGs(11,-11);
//  etaprimeSig->SetMotherPDGs(331,331);
//  etaprimeSig->SetMothersRelation(AliDielectronSignalMC::kSame);
//  etaprimeSig->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
//  etaprimeSig->SetMotherSources(AliDielectronSignalMC::kPrimary, AliDielectronSignalMC::kPrimary);
//  etaprimeSig->SetCheckBothChargesLegs(kTRUE,kTRUE);
//  etaprimeSig->SetCheckBothChargesMothers(kTRUE,kTRUE);
//  etaprimeSig->SetFillPureMCStep(kTRUE);
//  die->AddSignalMC(etaprimeSig);
//
//
//  AliDielectronSignalMC* rhoSig = new AliDielectronSignalMC("Rho", "rhoSignal"); ///rho pairs 
//  rhoSig->SetLegPDGs(11,-11);
//  rhoSig->SetMotherPDGs(113,113);
//  rhoSig->SetMothersRelation(AliDielectronSignalMC::kSame);
//  rhoSig->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
//  rhoSig->SetMotherSources(AliDielectronSignalMC::kPrimary, AliDielectronSignalMC::kPrimary);
//  rhoSig->SetCheckBothChargesLegs(kTRUE,kTRUE);
//  rhoSig->SetCheckBothChargesMothers(kTRUE,kTRUE);
//  rhoSig->SetFillPureMCStep(kTRUE);
//  die->AddSignalMC(rhoSig);
//
//  AliDielectronSignalMC* omegaSig = new AliDielectronSignalMC("Omega", "omegaSignal"); ///omega pairs 
//  omegaSig->SetLegPDGs(11,-11);
//  omegaSig->SetMotherPDGs(223,223);
//  omegaSig->SetMothersRelation(AliDielectronSignalMC::kSame);
//  omegaSig->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
//  omegaSig->SetMotherSources(AliDielectronSignalMC::kPrimary, AliDielectronSignalMC::kPrimary);
//  omegaSig->SetCheckBothChargesLegs(kTRUE,kTRUE);
//  omegaSig->SetCheckBothChargesMothers(kTRUE,kTRUE);
//  omegaSig->SetFillPureMCStep(kTRUE);
//  die->AddSignalMC(omegaSig);
//  
//  AliDielectronSignalMC* phiSig = new AliDielectronSignalMC("Phi", "phiSignal"); ///phi pairs 
//  phiSig->SetLegPDGs(11,-11);
//  phiSig->SetMotherPDGs(333,333);
//  phiSig->SetMothersRelation(AliDielectronSignalMC::kSame);
//  phiSig->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
//  phiSig->SetMotherSources(AliDielectronSignalMC::kPrimary, AliDielectronSignalMC::kPrimary);
//  phiSig->SetCheckBothChargesLegs(kTRUE,kTRUE);
//  phiSig->SetCheckBothChargesMothers(kTRUE,kTRUE);
//  phiSig->SetFillPureMCStep(kTRUE);
//  die->AddSignalMC(phiSig);
//  
//  // ##################### "real" pairs from photon conversions in the detector material ##############################
//  AliDielectronSignalMC* signalFromResonance_ULS_gammaConv = new AliDielectronSignalMC("signalFromResonance_ULS_gammaConv", "signalFromResonance_ULS_gammaConv");
//  signalFromResonance_ULS_gammaConv->SetLegPDGs(11,-11);
//  signalFromResonance_ULS_gammaConv->SetMotherPDGs(22,22);
//  signalFromResonance_ULS_gammaConv->SetMothersRelation(AliDielectronSignalMC::kSame);
//  signalFromResonance_ULS_gammaConv->SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary); // kSecondary means decays in the detector
//  signalFromResonance_ULS_gammaConv->SetCheckBothChargesLegs(kTRUE,kTRUE);
//  die->AddSignalMC(signalFromResonance_ULS_gammaConv);
//
//  // ##################### combinatorial pairs ##############################
//  AliDielectronSignalMC* diEleCombinatiorial = new AliDielectronSignalMC("diEleCombinatiorial", "diEleCombinatiorial");
//  diEleCombinatiorial->SetLegPDGs(11,-11);
//  diEleCombinatiorial->SetMothersRelation(AliDielectronSignalMC::kDifferent);
//  diEleCombinatiorial->SetCheckBothChargesLegs(kTRUE,kTRUE);
//  die->AddSignalMC(diEleCombinatiorial);
//
//  AliDielectronSignalMC* diEleCombinatiorialConversion = new AliDielectronSignalMC("diEleCombinatiorialConversion", "diEleCombinatiorialConversion");
//  diEleCombinatiorialConversion->SetLegPDGs(11,-11);
//  diEleCombinatiorialConversion->SetMotherPDGs(22,0);// 1 leg from photons + 1 leg from everything
//  diEleCombinatiorialConversion->SetMothersRelation(AliDielectronSignalMC::kDifferent);
//  diEleCombinatiorialConversion->SetCheckBothChargesLegs(kTRUE,kTRUE);
//  die->AddSignalMC(diEleCombinatiorialConversion);
//
//    // ##################### HF pairs ##############################
//  AliDielectronSignalMC* diEleHF = new AliDielectronSignalMC("diEleHF", "diEleHF");
//  diEleHF->SetLegPDGs(11,-11);
//  diEleHF->SetMotherPDGs(401,401);
//  diEleHF->SetMothersRelation(AliDielectronSignalMC::kDifferent);
//  diEleHF->SetCheckBothChargesLegs(kTRUE,kTRUE);
//  die->AddSignalMC(diEleHF);
//
// 
//  // #################### D-Mesons
//  AliDielectronSignalMC* diEleOpenCharmCharged = new AliDielectronSignalMC("DmesonsCharged","di-electrons from open charm D+- mesons no B grandmother");  // dielectrons originating from open charm hadrons
//  diEleOpenCharmCharged->SetLegPDGs(11,-11);
//  diEleOpenCharmCharged->SetMotherPDGs(401,401);
//  diEleOpenCharmCharged->SetLegSources(AliDielectronSignalMC::kDontCare, AliDielectronSignalMC::kDontCare);
//  diEleOpenCharmCharged->SetMothersRelation(AliDielectronSignalMC::kDifferent);
//  // diEleOpenCharmCharged->SetFillPureMCStep(kTRUE);
//  diEleOpenCharmCharged->SetCheckStackForPDG(kTRUE);
//  diEleOpenCharmCharged->SetPDGforStack(503);
//  diEleOpenCharmCharged->SetCheckBothChargesLegs(kTRUE,kTRUE);
//  diEleOpenCharmCharged->SetCheckBothChargesMothers(kTRUE,kTRUE);
//  die->AddSignalMC(diEleOpenCharmCharged);
//
//  AliDielectronSignalMC* diEleOpenCharmNeutral = new AliDielectronSignalMC("DmesonsNeutral","di-electrons from open charm D0 mesons no B grandmother");  // dielectrons originating from open charm hadrons
//  diEleOpenCharmNeutral->SetLegPDGs(11,-11);
//  diEleOpenCharmNeutral->SetMotherPDGs(405,405);
//  diEleOpenCharmNeutral->SetLegSources(AliDielectronSignalMC::kDontCare, AliDielectronSignalMC::kDontCare);
//  diEleOpenCharmNeutral->SetMothersRelation(AliDielectronSignalMC::kDifferent);
//  diEleOpenCharmNeutral->SetCheckStackForPDG(kTRUE);
//  diEleOpenCharmNeutral->SetPDGforStack(503);
//  diEleOpenCharmNeutral->SetCheckBothChargesLegs(kTRUE,kTRUE);
//  diEleOpenCharmNeutral->SetCheckBothChargesMothers(kTRUE,kTRUE);
//  die->AddSignalMC(diEleOpenCharmNeutral);
//
//  //B meson (3)
//  AliDielectronSignalMC* diEleOneOpenB = new AliDielectronSignalMC("B2ee","di-electrons from one B meson");  // dielectrons originating from open charm hadrons
//  diEleOneOpenB->SetLegPDGs(11,-11);
//  diEleOneOpenB->SetMotherPDGs(401,501);
//  diEleOneOpenB->SetLegSources(AliDielectronSignalMC::kDontCare, AliDielectronSignalMC::kDontCare);
//  diEleOneOpenB->SetGrandMotherPDGs(501,0);
//  diEleOneOpenB->SetCheckMotherGrandmotherRelation(kTRUE,kTRUE);
//  diEleOneOpenB->SetCheckBothChargesLegs(kTRUE,kTRUE);
//  diEleOneOpenB->SetCheckBothChargesMothers(kTRUE,kTRUE);
//  diEleOneOpenB->SetCheckBothChargesGrandMothers(kTRUE,kTRUE);
//  die->AddSignalMC(diEleOneOpenB);
//
//  //B meson (1)(1)
//  AliDielectronSignalMC* diEleOpenB = new AliDielectronSignalMC("BMesons","di-electrons from B mesons");  // dielectrons originating from open charm hadrons
//  diEleOpenB->SetLegPDGs(11,-11);
//  diEleOpenB->SetMotherPDGs(501,501);
//  diEleOpenB->SetLegSources(AliDielectronSignalMC::kDontCare, AliDielectronSignalMC::kDontCare);
//  diEleOpenB->SetMothersRelation(AliDielectronSignalMC::kDifferent);
//  diEleOpenB->SetCheckBothChargesLegs(kTRUE,kTRUE);
//  diEleOpenB->SetCheckBothChargesMothers(kTRUE,kTRUE);
//  die->AddSignalMC(diEleOpenB);
//
//  //B meson (2)(2)
//  AliDielectronSignalMC* diEleOpenBtoD = new AliDielectronSignalMC("B2D2ee","di-electrons from B->D-> e");  // dielectrons originating from open charm hadrons
//  diEleOpenBtoD->SetLegPDGs(11,-11);
//  diEleOpenBtoD->SetMotherPDGs(401,401);
//  diEleOpenBtoD->SetLegSources(AliDielectronSignalMC::kDontCare, AliDielectronSignalMC::kDontCare);
//  diEleOpenBtoD->SetGrandMotherPDGs(501,501);
//  diEleOpenBtoD->SetGrandMothersRelation(AliDielectronSignalMC::kDifferent);
//  diEleOpenBtoD->SetCheckBothChargesLegs(kTRUE,kTRUE);
//  diEleOpenBtoD->SetCheckBothChargesMothers(kTRUE,kTRUE);
//  diEleOpenBtoD->SetCheckBothChargesGrandMothers(kTRUE,kTRUE);
//  die->AddSignalMC(diEleOpenBtoD);
//
//  //B meson (1)(2)
//  AliDielectronSignalMC* diEleOpenBandBtoD = new AliDielectronSignalMC("B2eAndB2D2e","di-electrons from B->e and B->D->e");  // dielectrons originating from open charm hadrons
//  diEleOpenBandBtoD->SetLegPDGs        (11,11);
//  diEleOpenBandBtoD->SetMotherPDGs     (401,501);
//  diEleOpenBandBtoD->SetGrandMotherPDGs(501,0);
//  diEleOpenBandBtoD->SetLegSources(AliDielectronSignalMC::kDontCare, AliDielectronSignalMC::kDontCare);
//  diEleOpenBandBtoD->SetCheckBothChargesLegs(kTRUE,kTRUE);
//  diEleOpenBandBtoD->SetCheckBothChargesMothers(kTRUE,kTRUE);
//  diEleOpenBandBtoD->SetCheckBothChargesGrandMothers(kTRUE,kTRUE);
//  //do i need this?
//  diEleOpenBandBtoD->SetCheckMotherGrandmotherRelation(kTRUE,kFALSE);
//  die->AddSignalMC(diEleOpenBandBtoD);
  
  
  AliDielectronSignalMC* diElePhoto = new AliDielectronSignalMC("sameMotherEle","dielectrons from photo production");  
  diElePhoto->SetLegPDGs(11,-11);
//  diElePhoto->SetMotherPDGs(-9999,-9999);
  diElePhoto->SetMothersRelation(AliDielectronSignalMC::kSame);  
  diElePhoto->SetLegSources(AliDielectronSignalMC::kDontCare, AliDielectronSignalMC::kDontCare);
//  diElePhoto->SetGrandMotherPDGs(501,501);
//  diElePhoto->SetGrandMothersRelation(AliDielectronSignalMC::kDifferent);
  diElePhoto->SetCheckBothChargesLegs(kTRUE,kTRUE);
//  diElePhoto->SetCheckBothChargesMothers(kTRUE,kTRUE);
//  diElePhoto->SetCheckBothChargesGrandMothers(kTRUE,kTRUE);
  die->AddSignalMC(diElePhoto); 
  
}

TVectorD* getBinVec(int val){
  switch(val){
      case kMee:   return AliDielectronHelper::MakeArbitraryBinning(" 0.0,0.1,0.4,0.5 ,0.6 , 0.7 ,1.1, 1.5, 2.0 ,2.7 , 2.9, 3.1 , 5.0");
      case kPtee:  return AliDielectronHelper::MakeArbitraryBinning("0.0, 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 1, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0");
      case kPhi:   return AliDielectronHelper::MakeLinBinning(20,0,TMath::Pi());
      case kCent:  return AliDielectronHelper::MakeLinBinning(20,0,100);
      default: std::cout << "ERROR: in 'GetVector(...var)' variable for axis range not defined!" << std::endl;
      break;
  }
  return 0x0;      
}