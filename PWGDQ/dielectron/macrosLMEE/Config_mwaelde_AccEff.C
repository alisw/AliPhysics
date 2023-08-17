#include <AliAnalysisManager.h>
#include <AliAODInputHandler.h>

// cuts
// TString names("TrackCut_1;TrackCut_2");
// TString names("TrackCut_1;TrackCut_2;TrackCut_3");
TString names("TrackCut_1;TrackCut_2;TrackCut_3;TrackCut_4;TrackCut_5;TrackCut_6;TrackCut_7;TrackCut_8;TrackCut_9;TrackCut_10");
TObjArray *arrNames=names.Tokenize(";");
const Int_t nDie=arrNames->GetEntriesFast();
Int_t GetN(){return nDie;}

// generator index
const UInt_t genIDD[1] = {0}; // main generator
const Int_t ngenIDD =  0; // off
Int_t GetGenID(){return ngenIDD;}
UInt_t GetGenID(Int_t i) {return genIDD[i];}


// Binning of resolution histograms
const int    NbinsDeltaMom    = 1200; //12
Int_t GetNbinsDeltaMom() {return NbinsDeltaMom;}
const int    NbinsRelMom      = 400; //4
Int_t GetNbinsRelMom(){return NbinsRelMom;}
const int    NbinsDeltaEta    = 200; //2
Int_t GetNbinsDeltaEta(){return NbinsDeltaEta;}
const int    NbinsDeltaTheta  = 200; //2
Int_t GetNbinsDeltaTheta(){return NbinsDeltaTheta;}
const int    NbinsDeltaPhi    = 200; // 2
Int_t GetNbinsDeltaPhi(){return NbinsDeltaPhi;}

AliAnalysisCuts* SetupTrackCuts(Int_t cutDefinition);
AliAnalysisCuts* SetupPIDcuts(Int_t cutDefinition);
AliAnalysisFilter* Config_mwaelde_AccEff(Int_t cutDefinition);
AliAnalysisFilter* Config_mwaelde_AccEff_PID(Int_t cutDefinition);


// #########################################################
// #########################################################

//________________________________________________________________
AliAnalysisFilter* Config_mwaelde_AccEff(Int_t cutDefinition)
{
  std::cout << "SetupTrackCutsAndSettings()" <<std::endl;
  // cutDefinition = cutDefinition+11;

  AliAnalysisFilter *anaFilter = new AliAnalysisFilter(Form("anaFilter_Track_%d",cutDefinition),Form("anaFilter_Track_%d",cutDefinition)); // named constructor seems mandatory!

  anaFilter->AddCuts( SetupTrackCuts(cutDefinition) );
  std::cout << "...Trackcuts added!" <<std::endl;

  return anaFilter;
}


//________________________________________________________________
AliAnalysisCuts* SetupTrackCuts(Int_t cutDefinition)
{

   AliESDtrackCuts *fesdTrackCuts = new AliESDtrackCuts;

  // List of cuts
  fesdTrackCuts->SetPtRange( 0.2 , 10. );
  fesdTrackCuts->SetEtaRange( -0.8 , 0.8 );
  fesdTrackCuts->SetAcceptKinkDaughters(kFALSE);
  fesdTrackCuts->SetRequireSigmaToVertex(kFALSE);
  fesdTrackCuts->SetDCAToVertex2D(kFALSE);

  fesdTrackCuts->SetMaxDCAToVertexZ(3.);
  fesdTrackCuts->SetMaxDCAToVertexXY(1.);

  fesdTrackCuts->SetRequireITSRefit(kTRUE);
  fesdTrackCuts->SetRequireTPCRefit(kTRUE);

  // ITS number of clusters
  if((cutDefinition==1) || (cutDefinition==3) || (cutDefinition==5) || (cutDefinition==6) || (cutDefinition==8) || (cutDefinition==9) || (cutDefinition==10) || (cutDefinition==14) || (cutDefinition==15) || (cutDefinition==23) || (cutDefinition==24) || (cutDefinition==27)) fesdTrackCuts->SetMinNClustersITS(3);
  if((cutDefinition==2) || (cutDefinition==4) || (cutDefinition==7) || (cutDefinition==16) || (cutDefinition==18) || (cutDefinition==19) || (cutDefinition==21) || (cutDefinition==22) || (cutDefinition==25) || (cutDefinition==30)) fesdTrackCuts->SetMinNClustersITS(2);
  if((cutDefinition==11) || (cutDefinition==12) || (cutDefinition==13) || (cutDefinition==17) || (cutDefinition==20) || (cutDefinition==26) || (cutDefinition==28) || (cutDefinition==29)) fesdTrackCuts->SetMinNClustersITS(4);


  // ITS chi2
  if((cutDefinition==2) || (cutDefinition==5) || (cutDefinition==9) || (cutDefinition==11) || (cutDefinition==17) || (cutDefinition==24) || (cutDefinition==25) || (cutDefinition==29)) fesdTrackCuts->SetMaxChi2PerClusterITS(4.5);
  if((cutDefinition==1) || (cutDefinition==7) || (cutDefinition==8) || (cutDefinition==10) || (cutDefinition==12) || (cutDefinition==13) || (cutDefinition==16) || (cutDefinition==21) || (cutDefinition==26)) fesdTrackCuts->SetMaxChi2PerClusterITS(5.5);
  if((cutDefinition==3) || (cutDefinition==4) || (cutDefinition==6) || (cutDefinition==14) || (cutDefinition==15) || (cutDefinition==18) || (cutDefinition==19) || (cutDefinition==20) || (cutDefinition==22) || (cutDefinition==23) || (cutDefinition==27) || (cutDefinition==28) || (cutDefinition==30)) fesdTrackCuts->SetMaxChi2PerClusterITS(6.5);


  // ITS refit
  fesdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);

  // Min TPC number of clusters
  fesdTrackCuts->SetMinNClustersTPC(80);

  // Min TPC number of cross rows
  if((cutDefinition==1) || (cutDefinition==4) || (cutDefinition==11) || (cutDefinition==13) || (cutDefinition==15) || (cutDefinition==23) || (cutDefinition==24) || (cutDefinition==27) || (cutDefinition==29) || (cutDefinition==30)) fesdTrackCuts->SetMinNCrossedRowsTPC(100);
  if((cutDefinition==2) || (cutDefinition==3) || (cutDefinition==5) || (cutDefinition==8) || (cutDefinition==16) || (cutDefinition==17) || (cutDefinition==20) || (cutDefinition==25)) fesdTrackCuts->SetMinNCrossedRowsTPC(120);
  if((cutDefinition==6) || (cutDefinition==7) || (cutDefinition==9) || (cutDefinition==10) || (cutDefinition==12) || (cutDefinition==14) || (cutDefinition==18) || (cutDefinition==19) || (cutDefinition==21) || (cutDefinition==22) || (cutDefinition==26) || (cutDefinition==28)) fesdTrackCuts->SetMinNCrossedRowsTPC(140);

  // cross row/findable
  if((cutDefinition==1) || (cutDefinition==2) || (cutDefinition==3) || (cutDefinition==8) || (cutDefinition==9) || (cutDefinition==12) || (cutDefinition==14) || (cutDefinition==17) || (cutDefinition==18) || (cutDefinition==20) || (cutDefinition==21) || (cutDefinition==25) || (cutDefinition==30)) fesdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.6);
  if((cutDefinition==4) || (cutDefinition==6) || (cutDefinition==7) || (cutDefinition==10) || (cutDefinition==13) || (cutDefinition==15) || (cutDefinition==19) || (cutDefinition==22) || (cutDefinition==24) || (cutDefinition==26) || (cutDefinition==27) || (cutDefinition==29)) fesdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.7);
  if((cutDefinition==5) || (cutDefinition==11) || (cutDefinition==16) || (cutDefinition==23) || (cutDefinition==28)) fesdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.5);

  // max TPC chi2
  if((cutDefinition==1) || (cutDefinition==5) || (cutDefinition==8) || (cutDefinition==10) || (cutDefinition==11) || (cutDefinition==16) || (cutDefinition==19) || (cutDefinition==21) || (cutDefinition==27) || (cutDefinition==28)) fesdTrackCuts->SetMaxChi2PerClusterTPC(4);
  if((cutDefinition==2) || (cutDefinition==3) || (cutDefinition==4) || (cutDefinition==6) || (cutDefinition==9) || (cutDefinition==23)) fesdTrackCuts->SetMaxChi2PerClusterTPC(3);
  if((cutDefinition==7) || (cutDefinition==12) || (cutDefinition==13) || (cutDefinition==14) || (cutDefinition==15) || (cutDefinition==17) || (cutDefinition==18) || (cutDefinition==20) || (cutDefinition==22) || (cutDefinition==24) || (cutDefinition==25) || (cutDefinition==26) || (cutDefinition==29) || (cutDefinition==30)) fesdTrackCuts->SetMaxChi2PerClusterTPC(5);


  // shared TPC cls
  if(cutDefinition<11) fesdTrackCuts->SetMaxFractionSharedTPCClusters(0.4);
  if((cutDefinition>10) && (cutDefinition < 21)) fesdTrackCuts->SetMaxFractionSharedTPCClusters(0.6);
  if(cutDefinition>20) fesdTrackCuts->SetMaxFractionSharedTPCClusters(0.8);


   return fesdTrackCuts;

}
AliAnalysisFilter* Config_mwaelde_AccEff_PID(Int_t cutDefinition)
{
  std::cout << "SetupTrackCutsAndSettings()" <<std::endl;
  // cutDefinition = cutDefinition+11;

  AliAnalysisFilter *anaFilter_PID = new AliAnalysisFilter(Form("anaFilter_PID_%d",cutDefinition),Form("anaFilter_PID_%d",cutDefinition)); // named constructor seems mandatory!

  anaFilter_PID->AddCuts( SetupPIDcuts(cutDefinition) );
  std::cout << "...PIDcuts added!" <<std::endl;

  return anaFilter_PID;
}


//________________________________________________________________
AliAnalysisCuts* SetupPIDcuts(Int_t cutDefinition)
{
  AliAnalysisCuts* fancyCut=0x0;

  AliDielectronPID *mastermind_TPC = new AliDielectronPID("mastermind_TPC","mastermind_TPC");
  AliDielectronPID *mastermind_TOF = new AliDielectronPID("mastermind_TOF","mastermind_TOF");

  //TPC electrons: includes electrons and exclude all possible other contributions using the TPC

  // TPC kaon rejection
  if(cutDefinition==1) mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kKaon,    -3. ,3.,0.0, 100., kTRUE, AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==2) mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kKaon,    -2.5 ,3.5,0.0, 100., kTRUE, AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==3) mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kKaon,    -3. ,3.,0.0, 100., kTRUE, AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==4) mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kKaon,    -2.5 ,3.5,0.0, 100., kTRUE, AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==5) mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kKaon,    -3. ,3.,0.0, 100., kTRUE, AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==6) mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kKaon,    -3. ,3.,0.0, 100., kTRUE, AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==7) mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kKaon,    -3. ,3.,0.0, 100., kTRUE, AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==8) mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kKaon,    -3. ,3.,0.0, 100., kTRUE, AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==9) mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kKaon,    -3.5 ,2.5,0.0, 100., kTRUE, AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==10) mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kKaon,    -3.5 ,2.5,0.0, 100., kTRUE, AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==11) mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kKaon,    -3.5 ,2.5,0.0, 100., kTRUE, AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==12) mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kKaon,    -2.5 ,3.5,0.0, 100., kTRUE, AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==13) mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kKaon,    -3.5 ,2.5,0.0, 100., kTRUE, AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==14) mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kKaon,    -3. ,3.,0.0, 100., kTRUE, AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==15) mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kKaon,    -3. ,3.,0.0, 100., kTRUE, AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==16) mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kKaon,    -3.5 ,2.5,0.0, 100., kTRUE, AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==17) mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kKaon,    -3.5 ,2.5,0.0, 100., kTRUE, AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==18) mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kKaon,    -2.5 ,3.5,0.0, 100., kTRUE, AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==19) mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kKaon,    -3.5 ,2.5,0.0, 100., kTRUE, AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==20) mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kKaon,    -3. ,3.,0.0, 100., kTRUE, AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==21) mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kKaon,    -3.5 ,2.5,0.0, 100., kTRUE, AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==22) mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kKaon,    -2.5 ,3.5,0.0, 100., kTRUE, AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==23) mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kKaon,    -3.5 ,2.5,0.0, 100., kTRUE, AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==24) mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kKaon,    -2.5 ,3.5,0.0, 100., kTRUE, AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==25) mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kKaon,    -3.5 ,2.5,0.0, 100., kTRUE, AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==26) mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kKaon,    -2.5 ,3.5,0.0, 100., kTRUE, AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==27) mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kKaon,    -3. ,3.,0.0, 100., kTRUE, AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==28) mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kKaon,    -2.5 ,3.5,0.0, 100., kTRUE, AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==29) mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kKaon,    -3.5 ,2.5,0.0, 100., kTRUE, AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==30) mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kKaon,    -3.5 ,2.5,0.0, 100., kTRUE, AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);



  // TPC proton rejection
  if(cutDefinition==1) mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kProton,    -3. ,3.,0.0, 100., kTRUE, AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==2) mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kProton,    -2.5 ,3.5,0.0, 100., kTRUE, AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==3) mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kProton,    -2.5 ,3.5,0.0, 100., kTRUE, AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==4) mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kProton,    -3. ,3.,0.0, 100., kTRUE, AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==5) mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kProton,    -2.5 ,3.5,0.0, 100., kTRUE, AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==6) mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kProton,    -3. ,3.,0.0, 100., kTRUE, AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==7) mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kProton,    -2.5 ,3.5,0.0, 100., kTRUE, AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==8) mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kProton,    -3. ,3.,0.0, 100., kTRUE, AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==9) mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kProton,    -2.5 ,3.5,0.0, 100., kTRUE, AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==10) mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kProton,    -2.5 ,3.5,0.0, 100., kTRUE, AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==11) mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kProton,    -3. ,3.,0.0, 100., kTRUE, AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==12) mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kProton,    -3. ,3.,0.0, 100., kTRUE, AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==13) mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kProton,    -2.5 ,3.5,0.0, 100., kTRUE, AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==14) mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kProton,    -3.5 ,2.5,0.0, 100., kTRUE, AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==15) mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kProton,    -2.5 ,3.5,0.0, 100., kTRUE, AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==16) mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kProton,    -3. ,3.,0.0, 100., kTRUE, AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==17) mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kProton,    -2.5 ,3.5,0.0, 100., kTRUE, AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==18) mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kProton,    -3. ,3.,0.0, 100., kTRUE, AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==19) mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kProton,    -3. ,3.,0.0, 100., kTRUE, AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==20) mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kProton,    -3. ,3.,0.0, 100., kTRUE, AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==21) mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kProton,    -3.5 ,2.5,0.0, 100., kTRUE, AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==22) mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kProton,    -2.5 ,3.5,0.0, 100., kTRUE, AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==23) mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kProton,    -3.5 ,2.5,0.0, 100., kTRUE, AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==24) mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kProton,    -3.5 ,2.5,0.0, 100., kTRUE, AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==25) mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kProton,    -3.5 ,2.5,0.0, 100., kTRUE, AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==26) mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kProton,    -3. , 3.,0.0, 100., kTRUE, AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==27) mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kProton,    -3.  ,3.,0.0, 100., kTRUE, AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==28) mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kProton,    -3.5 ,2.5,0.0, 100., kTRUE, AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==29) mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kProton,    -3. ,3.,0.0, 100., kTRUE, AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==30) mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kProton,    -3. ,3.,0.0, 100., kTRUE, AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);

  //TPC electron
  if(cutDefinition==1) {
    Double_t sigmae = 3.;
    mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -sigmae , sigmae, 0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    mastermind_TOF->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -sigmae, sigmae, 0. ,100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
  else if(cutDefinition==2) {
    Double_t sigmae = 3.;
    mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -sigmae , sigmae, 0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    mastermind_TOF->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -sigmae, sigmae, 0. ,100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
  else if(cutDefinition==3) {
    Double_t sigmae = 3.;
    mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -sigmae , sigmae, 0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    mastermind_TOF->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -sigmae, sigmae, 0. ,100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
  else if(cutDefinition==4) {
    Double_t sigmae = 3.;
    mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -sigmae , sigmae, 0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    mastermind_TOF->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -sigmae, sigmae, 0. ,100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
  else if(cutDefinition==5) {
    Double_t sigmae = 3.;
    mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -sigmae , sigmae, 0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    mastermind_TOF->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -sigmae, sigmae, 0. ,100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
  else if(cutDefinition==6) {
    Double_t sigmae = 3.;
    mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -sigmae , sigmae, 0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    mastermind_TOF->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -sigmae, sigmae, 0. ,100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
  else if(cutDefinition==7) {
    Double_t sigmae = 3.;
    mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -sigmae , sigmae, 0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    mastermind_TOF->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -sigmae, sigmae, 0. ,100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
  else if(cutDefinition==8) {
    Double_t sigmae = 3.;
    mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -sigmae , sigmae, 0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    mastermind_TOF->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -sigmae, sigmae, 0. ,100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
  else if(cutDefinition==9) {
    Double_t sigmae = 3.;
    mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -sigmae , sigmae, 0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    mastermind_TOF->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -sigmae, sigmae, 0. ,100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
  else if(cutDefinition==10) {
    Double_t sigmae = 3.;
    mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -sigmae , sigmae, 0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    mastermind_TOF->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -sigmae, sigmae, 0. ,100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
  else if(cutDefinition==11) {
    Double_t sigmae = 2.5;
    mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -sigmae , sigmae, 0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    mastermind_TOF->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -sigmae, sigmae, 0. ,100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
  else if(cutDefinition==12) {
    Double_t sigmae = 2.5;
    mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -sigmae , sigmae, 0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    mastermind_TOF->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -sigmae, sigmae, 0. ,100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
  else if(cutDefinition==13) {
    Double_t sigmae = 2.5;
    mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -sigmae , sigmae, 0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    mastermind_TOF->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -sigmae, sigmae, 0. ,100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
  else if(cutDefinition==14) {
    Double_t sigmae = 2.5;
    mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -sigmae , sigmae, 0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    mastermind_TOF->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -sigmae, sigmae, 0. ,100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
  else if(cutDefinition==15) {
    Double_t sigmae = 2.5;
    mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -sigmae , sigmae, 0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    mastermind_TOF->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -sigmae, sigmae, 0. ,100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
  else if(cutDefinition==16) {
    Double_t sigmae = 2.5;
    mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -sigmae , sigmae, 0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    mastermind_TOF->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -sigmae, sigmae, 0. ,100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
  else if(cutDefinition==17) {
    Double_t sigmae = 2.5;
    mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -sigmae , sigmae, 0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    mastermind_TOF->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -sigmae, sigmae, 0. ,100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
  else if(cutDefinition==18) {
    Double_t sigmae = 2.5;
    mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -sigmae , sigmae, 0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    mastermind_TOF->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -sigmae, sigmae, 0. ,100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
  else if(cutDefinition==19) {
    Double_t sigmae = 2.5;
    mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -sigmae , sigmae, 0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    mastermind_TOF->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -sigmae, sigmae, 0. ,100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
  else if(cutDefinition==20) {
    Double_t sigmae = 2.5;
    mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -sigmae , sigmae, 0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    mastermind_TOF->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -sigmae, sigmae, 0. ,100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
  else if(cutDefinition==21) {
    Double_t sigmae = 3.5;
    mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -sigmae , sigmae, 0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    mastermind_TOF->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -sigmae, sigmae, 0. ,100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
  else if(cutDefinition==22) {
    Double_t sigmae = 3.5;
    mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -sigmae , sigmae, 0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    mastermind_TOF->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -sigmae, sigmae, 0. ,100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
  else if(cutDefinition==23) {
    Double_t sigmae = 3.5;
    mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -sigmae , sigmae, 0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    mastermind_TOF->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -sigmae, sigmae, 0. ,100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
  else if(cutDefinition==24) {
    Double_t sigmae = 3.5;
    mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -sigmae , sigmae, 0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    mastermind_TOF->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -sigmae, sigmae, 0. ,100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
  else if(cutDefinition==25) {
    Double_t sigmae = 3.5;
    mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -sigmae , sigmae, 0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    mastermind_TOF->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -sigmae, sigmae, 0. ,100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
  else if(cutDefinition==26) {
    Double_t sigmae = 3.5;
    mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -sigmae , sigmae, 0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    mastermind_TOF->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -sigmae, sigmae, 0. ,100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
  else if(cutDefinition==27) {
    Double_t sigmae = 3.5;
    mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -sigmae , sigmae, 0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    mastermind_TOF->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -sigmae, sigmae, 0. ,100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
  else if(cutDefinition==28) {
    Double_t sigmae = 3.5;
    mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -sigmae , sigmae, 0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    mastermind_TOF->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -sigmae, sigmae, 0. ,100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
  else if(cutDefinition==29) {
    Double_t sigmae = 3.5;
    mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -sigmae , sigmae, 0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    mastermind_TOF->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -sigmae, sigmae, 0. ,100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
  else if(cutDefinition==30) {
    Double_t sigmae = 3.5;
    mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -sigmae , sigmae, 0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    mastermind_TOF->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -sigmae, sigmae, 0. ,100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }


  // TPC pions rejection
  if(cutDefinition==1) {
    Double_t sigmapion = 3.5;
    mastermind_TOF->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -100, sigmapion, 0.0 ,100., kTRUE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,sigmapion,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
  else if(cutDefinition==2) {
    Double_t sigmapion = 3.5;
    mastermind_TOF->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -100, sigmapion, 0.0 ,100., kTRUE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,sigmapion,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
   else if(cutDefinition==3) {
    Double_t sigmapion = 3.5;
    mastermind_TOF->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -100, sigmapion, 0.0 ,100., kTRUE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,sigmapion,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
   else if(cutDefinition==4) {
    Double_t sigmapion = 4.;
    mastermind_TOF->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -100, sigmapion, 0.0 ,100., kTRUE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,sigmapion,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
   else if(cutDefinition==5) {
    Double_t sigmapion = 4.;
    mastermind_TOF->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -100, sigmapion, 0.0 ,100., kTRUE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,sigmapion,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
   else if(cutDefinition==6) {
    Double_t sigmapion = 4.;
    mastermind_TOF->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -100, sigmapion, 0.0 ,100., kTRUE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,sigmapion,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
   else if(cutDefinition==7) {
    Double_t sigmapion = 3.5;
    mastermind_TOF->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -100, sigmapion, 0.0 ,100., kTRUE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,sigmapion,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
   else if(cutDefinition==8) {
    Double_t sigmapion = 3.5;
    mastermind_TOF->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -100, sigmapion, 0.0 ,100., kTRUE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,sigmapion,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
   else if(cutDefinition==9) {
    Double_t sigmapion = 3.;
    mastermind_TOF->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -100, sigmapion, 0.0 ,100., kTRUE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,sigmapion,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
   else if(cutDefinition==10) {
    Double_t sigmapion = 3.5;
    mastermind_TOF->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -100, sigmapion, 0.0 ,100., kTRUE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,sigmapion,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
   else if(cutDefinition==11) {
    Double_t sigmapion = 3.0;
    mastermind_TOF->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -100, sigmapion, 0.0 ,100., kTRUE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,sigmapion,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
   else if(cutDefinition==12) {
    Double_t sigmapion = 3.5;
    mastermind_TOF->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -100, sigmapion, 0.0 ,100., kTRUE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,sigmapion,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
   else if(cutDefinition==13) {
    Double_t sigmapion = 3.0;
    mastermind_TOF->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -100, sigmapion, 0.0 ,100., kTRUE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,sigmapion,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
   else if(cutDefinition==14) {
    Double_t sigmapion = 3.0;
    mastermind_TOF->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -100, sigmapion, 0.0 ,100., kTRUE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,sigmapion,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
   else if(cutDefinition==15) {
    Double_t sigmapion = 4.0;
    mastermind_TOF->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -100, sigmapion, 0.0 ,100., kTRUE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,sigmapion,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
   else if(cutDefinition==16) {
    Double_t sigmapion = 3.5;
    mastermind_TOF->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -100, sigmapion, 0.0 ,100., kTRUE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,sigmapion,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
   else if(cutDefinition==17) {
    Double_t sigmapion = 4.0;
    mastermind_TOF->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -100, sigmapion, 0.0 ,100., kTRUE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,sigmapion,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
   else if(cutDefinition==18) {
    Double_t sigmapion = 3.0;
    mastermind_TOF->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -100, sigmapion, 0.0 ,100., kTRUE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,sigmapion,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
   else if(cutDefinition==19) {
    Double_t sigmapion = 4.0;
    mastermind_TOF->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -100, sigmapion, 0.0 ,100., kTRUE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,sigmapion,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
   else if(cutDefinition==20) {
    Double_t sigmapion = 3.0;
    mastermind_TOF->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -100, sigmapion, 0.0 ,100., kTRUE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,sigmapion,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
   else if(cutDefinition==21) {
    Double_t sigmapion = 3.5;
    mastermind_TOF->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -100, sigmapion, 0.0 ,100., kTRUE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,sigmapion,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
   else if(cutDefinition==22) {
    Double_t sigmapion = 4.0;
    mastermind_TOF->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -100, sigmapion, 0.0 ,100., kTRUE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,sigmapion,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
   else if(cutDefinition==23) {
    Double_t sigmapion = 4.0;
    mastermind_TOF->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -100, sigmapion, 0.0 ,100., kTRUE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,sigmapion,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
   else if(cutDefinition==24) {
    Double_t sigmapion = 3.0;
    mastermind_TOF->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -100, sigmapion, 0.0 ,100., kTRUE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,sigmapion,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
   else if(cutDefinition==25) {
    Double_t sigmapion = 3.0;
    mastermind_TOF->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -100, sigmapion, 0.0 ,100., kTRUE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,sigmapion,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
   else if(cutDefinition==26) {
    Double_t sigmapion = 4.0;
    mastermind_TOF->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -100, sigmapion, 0.0 ,100., kTRUE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,sigmapion,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
   else if(cutDefinition==27) {
    Double_t sigmapion = 4.0;
    mastermind_TOF->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -100, sigmapion, 0.0 ,100., kTRUE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,sigmapion,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
   else if(cutDefinition==28) {
    Double_t sigmapion = 3.0;
    mastermind_TOF->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -100, sigmapion, 0.0 ,100., kTRUE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,sigmapion,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
   else if(cutDefinition==29) {
    Double_t sigmapion = 3.5;
    mastermind_TOF->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -100, sigmapion, 0.0 ,100., kTRUE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,sigmapion,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
   else if(cutDefinition==30) {
    Double_t sigmapion = 3.0;
    mastermind_TOF->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -100, sigmapion, 0.0 ,100., kTRUE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,sigmapion,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }

  // TOF e sigma
  if(cutDefinition==1) mastermind_TOF->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3. , 0. ,100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==2) mastermind_TOF->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -2. , 2. , 0. ,100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==3) mastermind_TOF->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -2. , 2. , 0. ,100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==4) mastermind_TOF->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -2. , 2. , 0. ,100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==5) mastermind_TOF->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3. , 0. ,100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==6) mastermind_TOF->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -2. , 2. , 0. ,100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==7) mastermind_TOF->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3. , 0. ,100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==8) mastermind_TOF->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -2. , 2. , 0. ,100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==9) mastermind_TOF->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3. , 0. ,100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==10) mastermind_TOF->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3. , 0. ,100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==11) mastermind_TOF->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -2. , 2. , 0. ,100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==12) mastermind_TOF->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3. , 0. ,100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==13) mastermind_TOF->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3. , 0. ,100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==14) mastermind_TOF->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3. , 0. ,100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==15) mastermind_TOF->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -2. , 2. , 0. ,100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==16) mastermind_TOF->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -2. , 2. , 0. ,100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==17) mastermind_TOF->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3. , 0. ,100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==18) mastermind_TOF->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3. , 0. ,100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==19) mastermind_TOF->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -2. , 2. , 0. ,100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==20) mastermind_TOF->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -2. , 2. , 0. ,100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==21) mastermind_TOF->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -2. , 2. , 0. ,100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==22) mastermind_TOF->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -2. , 2. , 0. ,100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==23) mastermind_TOF->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3. , 0. ,100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==24) mastermind_TOF->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3. , 0. ,100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==25) mastermind_TOF->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -2. , 2. , 0. ,100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==26) mastermind_TOF->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -2. , 2. , 0. ,100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==27) mastermind_TOF->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -2. , 2. , 0. ,100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==28) mastermind_TOF->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3. , 0. ,100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==29) mastermind_TOF->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -2. , 2. , 0. ,100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  else if(cutDefinition==30) mastermind_TOF->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -2. , 2. , 0. ,100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);

  // Combine
  AliDielectronCutGroup* mastermind_cg = new AliDielectronCutGroup("mastermind_cg","mastermind_cg",AliDielectronCutGroup::kCompOR);
  mastermind_cg->AddCut(mastermind_TPC);
  mastermind_cg->AddCut(mastermind_TOF);
  fancyCut = mastermind_cg;


  return fancyCut;

}
