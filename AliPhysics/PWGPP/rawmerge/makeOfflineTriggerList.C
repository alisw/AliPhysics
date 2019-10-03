/*
   gSystem->AddIncludePath("-I$AliPhysics_SRC/PWGPP/ -I$AliPhysics_SRC/OADB/");  // ? why not in the alienv,  why not available ?
  .L $AliPhysics_SRC/PWGPP/rawmerge/makeOfflineTriggerList.C+
  makeOfflineTriggerList(gSystem->ExpandPathName("$NOTES/JIRA/PWGPP-227/data/2016/LHC16t/000267161/pass1_CENT_wSDD/filtered.list"))
  //
  // makeOfflineTriggerList for selected filtered trees
  // configuration parameters are taken fro env variables  defined in makeOfflineTriggerListSetupXXX.sh
  // 
   
  aliroot -b -q $AliPhysics_SRC/PWGPP/rawmerge/makeOfflineTriggerList.C+\(\"filtered.list\"\) |tee makeOfflineTriggerList.log

*/


#include "TFile.h"
#include "TCut.h"
#include "TH1.h"
#include "TSystem.h"
#include "TMath.h"
#include "TObjString.h"
#include "TTree.h"
#include "AliTreePlayer.h"
#include "AliXRDPROOFtoolkit.h"
#include "TTreeStream.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include "TPad.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TStatToolkit.h"
#include "AliExternalInfo.h"
#include "TStopwatch.h"
#include "AliDrawStyle.h"
#include "AliAnalysisTaskFilteredTree.h"

using namespace std;
using std::cout;
using std::setw;



void TriggerHighMultiplicity( const char * chinput,  const char * filter="ntracks<1000", Long64_t ntracksPrim=200000, Float_t fractionTracks=0.8, Long64_t nEvents=200000);
void TriggerHighPt( const char * chinput,  const char * filter="(esdTrack.fFlags&0x4401)>0", Long64_t nEvents=200000, Double_t dcaCut=4);
void TriggerHighPtV0s( const char * chinput,  const char * filter, Long64_t nEvents, Double_t zCut, Double_t covarQPtCut, Int_t filterMask);
void TriggerCosmicPairs( const char * chinput,  const char * filter="abs(t0.fP[4])<0.33", Long64_t nEvents=200000);
void TriggerLaser( const char * chinput, Long64_t nEvents=200000);
void triggerCalibHighPt( const char * chinput,  const char * filter, Long64_t nEvents, Double_t dcaCut, Double_t multRatioFraction, Double_t multFraction, Int_t maxTracks);
void triggerCalibV0( const char * chinput,  const char * filter, Long64_t nEvents, Double_t multRatioFraction, Double_t multFraction, Int_t maxTracks, Double_t zCut, Double_t covarQPtCut, Int_t filterMask);
void GetRawSummary();

 

void makeOfflineTriggerList(const char * chinput){
  //
  // 
  // Make high mulitplicity event selection for the calibration
  // TString highMultiplicityFilter="ntracks<1000";  
  if (TString(chinput).Contains("GetRawSummary")) {
    ::Info("makeOfflineTriggerList","makeOfflineTriggerList");
    return GetRawSummary();
  }
  AliDrawStyle::SetDefaults();
  AliDrawStyle::ApplyStyle("figTemplate");

  Long64_t  highMultiplicityNTracksPrim=500000;   // default 500000 tracks for the calibration
  Float_t highMultiplicityFractionTracks=0.8;   //
  Long64_t  highMultiplicityNEvents=10000000;
  // if (gSystem->Getenv("highMultiplicityFilter"))        highMultiplicityFilter=gSystem->Getenv("highMultiplicityFilter");
  if (gSystem->Getenv("highMultiplicityNTracksPrim"))   highMultiplicityNTracksPrim=TString(gSystem->Getenv("highMultiplicityNTracksPrim")).Atof();
  if (gSystem->Getenv("highMultiplicityFractionTracks"))highMultiplicityFractionTracks=TString(gSystem->Getenv("highMultiplicityFractionTracks")).Atof();
  if (gSystem->Getenv("highMultiplicityNEvents"))       highMultiplicityNEvents=TString(gSystem->Getenv("highMultiplicityNEvents")).Atof();
  // TriggerHighMultiplicity(chinput, highMultiplicityFilter.Data(),highMultiplicityNTracksPrim,  highMultiplicityFractionTracks, highMultiplicityNEvents);
  //
  // High pt filters
  TString highPtFilter="(esdTrack.fFlags&0x4401)>0&&esdTrack.Pt()>8";
  Long64_t   highPtNEvents=1000000000;
  Double_t highPtDCAcut=4;  
  if (gSystem->Getenv("highPtFilter"))        highPtFilter=gSystem->Getenv("highPtFilter");
  if (gSystem->Getenv("highPtNEvents"))       highPtNEvents=TString(gSystem->Getenv("highPtNEvents")).Atoi();
  if (gSystem->Getenv("highPtDCAcut"))        highPtDCAcut =TString(gSystem->Getenv("highPtDCAcut")).Atoi();
  TriggerHighPt( chinput, highPtFilter.Data(), highPtNEvents,  highPtDCAcut);
  // HighPt V0s
  //
  TString  highPtV0Filter="((type==1)*max(track0.Pt(),track1.Pt())>2)||(max(track0.Pt(),track1.Pt())>4)"; 
  Int_t    highPtV0NEvents=100000000; 
  Double_t highPtV0ZCut=20.;
  Double_t highPtV0CovarQPtCut=0.03;
  Int_t    highPtV0FilterMask=0x4401;
  if (gSystem->Getenv("highPtV0Filter"))        highPtV0Filter=gSystem->Getenv("highPtV0Filter");
  if (gSystem->Getenv("highPtV0NEvents"))       highPtV0NEvents=TString(gSystem->Getenv("highPtV0NEvents")).Atoi();
  if (gSystem->Getenv("highPtV0ZCut"))          highPtV0ZCut=TString(gSystem->Getenv("highPtV0ZCut")).Atoi();
  if (gSystem->Getenv("highPtV0CovarQPtCut"))   highPtV0CovarQPtCut=TString(gSystem->Getenv("highPtV0CovarQPtCut")).Atof();
  if (gSystem->Getenv("highPtV0CovarFilterMask"))   highPtV0FilterMask=TString(gSystem->Getenv("highPtV0FilterMask")).Atoi();
  TriggerHighPtV0s(chinput, highPtV0Filter.Data(),  highPtV0NEvents, highPtV0ZCut,  highPtV0CovarQPtCut,  highPtV0FilterMask);  
  // CosmicPairs
  //
  TString  cosmicsFilter="abs(t0.fP[4])<0.33"; 
  Int_t    cosmicsNEvents=100000000; 
  if (gSystem->Getenv("cosmicsFilter"))        cosmicsFilter=gSystem->Getenv("cosmicsFilter");
  if (gSystem->Getenv("cosmicsNEvents"))       cosmicsNEvents=TString(gSystem->Getenv("cosmicsNEvents")).Atoi();
  TriggerCosmicPairs(chinput, cosmicsFilter.Data(),  cosmicsNEvents);  
  // Laser
  //
  Int_t    laserNEvents=100000000; 
  if (gSystem->Getenv("laserNEvents"))       laserNEvents=TString(gSystem->Getenv("cosmicsNEvents")).Atoi();
  TriggerLaser(chinput,  laserNEvents);  
  //
  // calibration trigger && highpt - high multiplicity low backgorund, enhanced high pt
  TString  calibhighPtFilter="(esdTrack.fFlags&0x4401)>0&&esdTrack.Pt()>2";
  Long64_t   calibhighPtNEvents=1000000000;
  Double_t calibhighPtDCAcut=4;  
  if (gSystem->Getenv("calibhighPtFilter"))        highPtFilter=gSystem->Getenv("calibhighPtFilter");
  if (gSystem->Getenv("calibhighPtNEvents"))       highPtNEvents=TString(gSystem->Getenv("calibhighPtNEvents")).Atoi();
  if (gSystem->Getenv("calibhighPtDCAcut"))        highPtDCAcut =TString(gSystem->Getenv("calibhighPtDCAcut")).Atoi();
  triggerCalibHighPt(chinput, calibhighPtFilter,calibhighPtNEvents, calibhighPtDCAcut, highMultiplicityFractionTracks, 0.5, highMultiplicityNTracksPrim); 
  // calibration trigger && highpt V0 - high multiplicity low backgorund, enhanced high pt
  TString  calibhighPtV0Filter="((type==1)*max(track0.Pt(),track1.Pt())>1)||(max(track0.Pt(),track1.Pt())>2)"; 
  triggerCalibV0(chinput,calibhighPtV0Filter, calibhighPtNEvents, highMultiplicityFractionTracks, 0.5, highMultiplicityNTracksPrim, highPtV0ZCut,  highPtV0CovarQPtCut,  highPtV0FilterMask);

}


void TriggerHighMultiplicity( const char * chinput,  const char * filter, Long64_t ntracksPrim, Float_t fractionTracks, Long64_t nEvents){
  /*
    chinput=gSystem->ExpandPathName("$NOTES/JIRA/PWGPP-227/data/2016/LHC16t/000267161/pass1_CENT_wSDD/filtered.list");
    filter="ntracks<10000"; 
    Int_t ntracksPrim=50000;
    Float_t fractionTracks=0.8; 
  */
  //Int_t nEvents=100000;
  const char * treeName="eventInfoV0";
  TTree *tree=NULL; 
  if (fractionTracks>1) {
    ::Error("TriggerHighMultiplicity","Invalid fraction of tracks - %f should be in <0,1",fractionTracks);
  }
  if (TString(chinput).Contains(".root")){
    TFile * fInputFile = TFile::Open(chinput);
    tree               = (TTree*)fInputFile->Get(treeName);
  } else {
    tree=(TTree*)AliXRDPROOFtoolkit::MakeChain(chinput, treeName, 0, 1000000000,0,1);
  }
  tree->SetEstimate(tree->GetEntries());
  TTreeSRedirector * pcstream = new TTreeSRedirector("filteredMult.root","recreate");
  //
  //  step1 - get cumulat fractionTracks from mult/ntracks 
  Long64_t events = tree->Draw("mult/ntracks",filter,"",nEvents);
  Float_t fractionTracksCut = TMath::KOrdStat(events, tree->GetV1(), Long64_t(Double_t(events)*fractionTracks));
  tree->GetHistogram()->Write("multNTracks");
  //
  // step 2 select events  
  TCut weightFraction=TString::Format("mult*((mult/ntracks>%f)&&%s)",fractionTracksCut,filter).Data();  
  tree->Draw("mult",weightFraction,"",nEvents);
  tree->GetHistogram()->Write("multWithWeight");
  Int_t nBins = tree->GetHistogram()->GetXaxis()->GetNbins();
  Double_t cutMult=0;
  Double_t tracksAll=tree->GetHistogram()->Integral();
  Double_t *integral=tree->GetHistogram()->GetIntegral();
  for (Int_t iBin=nBins-1; iBin>0;iBin--){
    if ((1-integral[iBin])*tracksAll<ntracksPrim){
      cutMult=tree->GetHistogram()->GetBinLowEdge(iBin);
      continue;
    }
  }
  // Dump input parameters
  ::Info("KeyValue.TriggerHighMultiplicity.Input.Filter","%s",filter);
  ::Info("KeyValue.TriggerHighMultiplicity.Input.ntracksPrim","%llu",ntracksPrim);
  ::Info("KeyValue.TriggerHighMultiplicity.Input.fractionTracks","%f",fractionTracks);
  ::Info("KeyValue.TriggerHighMultiplicity.Input.nEvent","%llu",nEvents);
  //
  // Dump selected Events
  tree->GetUserInfo()->AddLast(new TObjString("highMult"));
  Int_t triggerIndex = tree->GetUserInfo()->GetEntries()-1;  
  tree->SetAlias("cutMultSelect",TString::Format("((mult/ntracks>%f)&&%s)&&mult>%f",fractionTracksCut,filter,cutMult).Data());
  ::Info("TriggerHighMultiplicity","Selected events cut = %s", tree->GetAlias("cutMultSelect"));
  TString outputString="";
  outputString+="fileName.GetString();1;fname;/C:";
  outputString+="run;1;run;/i:";
  outputString+="evtNumberInFile;1.20;eventID;/i:";
  outputString+="time;1.20;timeStamp;/i:";
  outputString+="gid;1.30;gid;/l:";
  outputString+=TString::Format("This->GetUserInfo()->At(%d)->GetName();1;trigger;/C:",triggerIndex);
  outputString+="triggerClass.GetName();20.20;triggerClass;/C:";
  //
  AliTreePlayer::selectWhatWhereOrderBy(tree, outputString.Data(), "cutMultSelect","",0,nEvents, "csv","filteredMult.list");
  // Export counter - histograms in filteredMult.root file and log file
  Int_t nEventsAll = tree->Draw("mult","","",nEvents);
  tree->GetHistogram()->Write("multAll");
  tree->Draw("mult","mult","",nEvents);
  tree->GetHistogram()->Write("multPrimAllWeighted");
  Int_t nTracksPrimAll= tree->GetHistogram()->Integral();
  tree->Draw("ntracks","ntracks","",nEvents);
  tree->GetHistogram()->Write("multTracksAllWeighted");
  Double_t nTracksAll= tree->GetHistogram()->Integral();
  //
  Int_t nHighMultiplicitySelected = tree->Draw("mult", "cutMultSelect","",nEvents);
  tree->GetHistogram()->Write("multSelected");
  tree->Draw("mult", "mult*(cutMultSelect)","",nEvents);
  tree->GetHistogram()->Write("multPrimSelectedWeighted");
  Int_t nTracksPrimSelected= tree->GetHistogram()->Integral();
  tree->Draw("ntracks", "ntracks*(cutMultSelect)","",nEvents);
  tree->GetHistogram()->Write("multSelectedWeighted");
  Int_t nTracksSelected= tree->GetHistogram()->Integral();
  // 
  tree->Draw("run");
  // Dump counters
  ::Info("KeyValue.TriggerHighMultiplicity.Selection","%s",tree->GetAlias("cutMultSelect"));
  ::Info("KeyValue.TriggerHighMultiplicity.NEventsAll","%d",nEventsAll);
  ::Info("KeyValue.TriggerHighMultiplicity.NHighMultiplicitySelected","%d",nHighMultiplicitySelected);
  ::Info("KeyValue.TriggerHighMultiplicity.NTracksAll","%.0f",nTracksAll);
  ::Info("KeyValue.TriggerHighMultiplicity.NTracksSelected","%d",nTracksSelected);
  ::Info("KeyValue.TriggerHighMultiplicity.NTracksPrimAll","%d",nTracksPrimAll);
  ::Info("KeyValue.TriggerHighMultiplicity.NTracksPrimSelected","%d",nTracksPrimSelected);
  delete pcstream;
}


void TriggerHighPt( const char * chinput,  const char * filter, Long64_t nEvents, Double_t dcaCut){
  /*
    chinput=gSystem->ExpandPathName("$NOTES/JIRA/PWGPP-227/data/2016/LHC16t/000267161/pass1_CENT_wSDD/filtered.list");
    filter="(esdTrack.fFlags&0x4401)>0&&esdTrack.Pt()>8"; 
    Int_t ntracksPrim=50000;
    Int_t nEvents=10000; 
    Double_t dcaCut=3;    
  */
  const char * treeName="highPt";
  TTree *tree=NULL; 
  TTree *treeEvent=NULL; 

  if (TString(chinput).Contains(".root")){
    TFile * fInputFile = TFile::Open(chinput);
    tree               = (TTree*)fInputFile->Get(treeName);
    treeEvent          = (TTree*)fInputFile->Get("eventInfoV0");
  } else {
    tree=(TTree*)AliXRDPROOFtoolkit::MakeChain(chinput, treeName, 0, 1000000000,0,1);
    treeEvent=(TTree*)AliXRDPROOFtoolkit::MakeChain(chinput, "eventInfoV0", 0, 1000000000,0,1);
  }
  
  tree->SetEstimate(tree->GetEntries());
  AliAnalysisTaskFilteredTree::SetDefaultAliasesHighPt(tree);
  TTreeSRedirector * pcstream = new TTreeSRedirector("filteredHighPt.root","recreate");
  // Dump selected Events
  tree->GetUserInfo()->AddLast(new TObjString("highPt"));
  Int_t triggerIndex = tree->GetUserInfo()->GetEntries()-1;  
  tree->SetAlias("cutHighPt",TString::Format("abs(esdTrack.fdTPC)<%.3f&&abs(esdTrack.fzTPC)<%.3f&&esdTrack.fzTPC!=0&&sqrt(esdTrack.fC[14])<0.01&&%s",dcaCut,dcaCut, filter).Data());
  //
  TString outputString="";
  outputString+="fileName.GetString();1;fname;/C:";
  outputString+="runNumber;1;run;/i:";
  outputString+="evtNumberInFile;1.20;eventID;/i:";
  outputString+="evtTimeStamp;1.20;timeStamp;/i:";
  outputString+="gid;1.30;gid;/l:";
  outputString+=TString::Format("This->GetUserInfo()->At(%d)->GetName();1;trigger;/C:",triggerIndex);
  outputString+="triggerClass.GetName();20.20;triggerClass;/C:";
  AliTreePlayer::selectWhatWhereOrderBy(tree, outputString.Data(), "cutHighPt","",0,nEvents, "csv","filteredHighPt.list");
  //
  tree->Draw("-esdTrack.Pt():mult:ntracks","cutHighPt&&esdTrack.Pt()<100","goff",nEvents);
  Double_t meanMult=TMath::Mean(tree->GetSelectedRows(),tree->GetV2());
  Double_t meanNTracks=TMath::Mean(tree->GetSelectedRows(),tree->GetV3());
  Double_t mult95=TMath::KOrdStat(tree->GetSelectedRows(),tree->GetV2(),Long64_t(tree->GetSelectedRows()*0.95));
  Double_t ntracks95=TMath::KOrdStat(tree->GetSelectedRows(),tree->GetV3(),Long64_t(tree->GetSelectedRows()*0.95)); 
  Double_t ptQuant95=-TMath::KOrdStat(tree->GetSelectedRows(),tree->GetV1(),Long64_t(tree->GetSelectedRows()*0.05));
  //
  // Make QA plots
  //
  // QA plots
  TString qaHistograms="";
  qaHistograms+=TString::Format("mult:#1>>hisMultAll(%.0f,0,%.0f);",TMath::Min(2*mult95,200.),2*mult95);
  qaHistograms+=TString::Format("mult:#cutHighPt>>hisMultSelected(%.0f,0,%.0f);",TMath::Min(2*mult95,200.),2*mult95);
  qaHistograms+=TString::Format("ntracks:#1>>hisNtracksAll(%.0f,0,%.0f);",TMath::Min(ntracks95,200.),2*ntracks95);
  qaHistograms+=TString::Format("ntracks:#cutHighPt>>hisNtracksSelected(%.0f,0,%.0f);",TMath::Min(ntracks95,200.),2*ntracks95);
  qaHistograms+=TString::Format("mult/ntracks:#1>>hisMultNtracksAll(100,0,1);");
  qaHistograms+=TString::Format("mult/ntracks:#cutHighPt>>hisMultNtracksSelected(100,0,1);");
  qaHistograms+=TString::Format("qPt:#1>>hisQPtAll(200,-5,5);",2.*ptQuant95);
  qaHistograms+=TString::Format("qPt:#cutHighPt>>hisQPtSelected(200,-5,5);",2.*ptQuant95);
  TStopwatch timer; TObjArray *hisArray = AliTreePlayer::MakeHistograms(tree, qaHistograms, "cutHighPt",0,nEvents,nEvents,15); timer.Print();
  //
  TString drawExpression="";
  drawExpression="[1,1,1,1]:";
  drawExpression+="%Ogridx,gridy;hisMultAll(0,10000)(0)(err);hisMultSelected(0,10000)(0)(err):";
  drawExpression+="%Ogridx,gridy;hisMultNtracksAll(0,10000)(0)(err);hisMultNtracksSelected(0,10000)(0)(err):";
  drawExpression+="%Ogridx,gridy;hisNtracksAll(0,10000)(0)(err);hisNtracksSelected(0,10000)(0)(err):";
  drawExpression+="%Ogridx,gridy;hisQPtAll(0,200)(0)(err);hisQPtSelected(0,200)(0)(err):";
  TObjArray *keepArray = new TObjArray(100);
  TPad * pad = AliTreePlayer::DrawHistograms(0,hisArray,drawExpression,keepArray, 15); 
  ((TCanvas*)pad)->SetWindowSize(1800,1000);
  pad->SaveAs("triggerHighPt.png");
  pad->SaveAs("triggerHighPt.C");
  pcstream->GetFile()->cd();
  hisArray->Write("hisArray",TObjArray::kSingleKey);
  keepArray->Write("keepArray",TObjArray::kSingleKey);
  //
  //
 
  TString eventListName = "filteredHighPt.list";
  AliTreePlayer::selectWhatWhereOrderBy(tree, outputString.Data(), "cutHighPt","",0,nEvents, "csv",eventListName);
  gSystem->Exec(Form("{ rm %s && uniq > %s; } < %s ",eventListName.Data(),eventListName.Data(),eventListName.Data()));
  //
  Int_t nHighPtSelected = tree->Draw("esdTrack.fTPCncls","cutHighPt");
  Int_t nEventsAll=treeEvent->GetEntries();
   // Dump input parameters
  
  ::Info("KeyValue.TriggerHighPt.Input.Filter","%s",filter);
  ::Info("KeyValue.TriggerHighPt.Input.dcaCut","%f",dcaCut);
  ::Info("KeyValue.TriggerHighPt.Input.nEvents","%llu",nEvents);
  // Dump counters
  ::Info("KeyValue.TriggerHighPt.Selection","%s",tree->GetAlias("cutHighPt"));
  ::Info("KeyValue.TriggerHighPt.NEventsAll","%d",nEventsAll);
  ::Info("KeyValue.TriggerHighPt.NHighPtSelected","%d",nHighPtSelected);
  ::Info("KeyValue.TriggerHighPt.MeanMult","%f",meanMult);  
  ::Info("KeyValue.TriggerHighPt.Mult95","%f", mult95);  
  ::Info("KeyValue.TriggerHighPt.NTracks95","%f",ntracks95);  
  delete pcstream;
}


void TriggerHighPtV0s( const char * chinput,  const char * filter, Long64_t nEvents, Double_t zCut, Double_t covarQPtCut, Int_t filterMask){
  /*
    chinput=gSystem->ExpandPathName("$NOTES/JIRA/PWGPP-227/data/2016/LHC16t/000267161/pass1_CENT_wSDD/filtered.list");
    filter="((type==1)*max(track0.Pt(),track1.Pt())>2)||(max(track0.Pt(),track1.Pt())>4)"; 
    Int_t nEvents=10000; 
    Double_t zCut=20.;
    Double_t covarQPtCut=0.03;
    filterMask=0x4401;
    TriggerHighPtV0s(chinput, filter, nEvents, zCut,covarQPtCut,filterMask);

  */
  const char * treeName="V0s";
  TTree *tree=NULL; 
  TTree *treeEvent=NULL; 
  if (TString(chinput).Contains(".root")){
    TFile * fInputFile = TFile::Open(chinput);
    tree               = (TTree*)fInputFile->Get(treeName);
    treeEvent          = (TTree*)fInputFile->Get("eventInfoV0");
  } else {
    tree=(TTree*)AliXRDPROOFtoolkit::MakeChain(chinput, treeName, 0, 1000000000,0,1);
    treeEvent=(TTree*)AliXRDPROOFtoolkit::MakeChain(chinput, "eventInfoV0", 0, 1000000000,0,1);
  }
  tree->SetEstimate(tree->GetEntries()); 
  AliAnalysisTaskFilteredTree::SetDefaultAliasesV0(tree);
  TTreeSRedirector * pcstream = new TTreeSRedirector("filteredHighPtV0s.root","recreate");
  //
  tree->SetAlias("filterMask",TString::Format("((track0.fFlags&%x)+(track1.fFlags&%x)>0)",filterMask,filterMask).Data());
  tree->SetAlias("dcaZcut",TString::Format("abs(track0.fzTPC)<%.3f&&abs(track1.fzTPC)<%.3f",zCut,zCut).Data());
  tree->SetAlias("covarQPtCut",TString::Format("sqrt(track0.fC[14])<%.3f&&sqrt(track1.fC[14])<%0.3f",covarQPtCut,covarQPtCut).Data());
  tree->SetAlias("userCut",filter);
  TString cutAll=TString::Format("filterMask&&dcaZcut&&covarQPtCut&&%s",filter);
  // Dump input parameters
  ::Info("KeyValue.TriggerHighPtV0s.Input.filter","%s",filter);
  ::Info("KeyValue.TriggerHighPtV0s.Input.nEvents","%llu",nEvents);
  ::Info("KeyValue.TriggerHighPtV0s.Input.zCut","%f",zCut);
  ::Info("KeyValue.TriggerHighPtV0s.Input.covarQPtCut","%f",covarQPtCut);
  ::Info("KeyValue.TriggerHighPtV0s.Input.filterMask","%x",filterMask);
  //
  // Dump selected Events
  tree->GetUserInfo()->AddLast(new TObjString("highPtV0s"));
  Int_t triggerIndex = tree->GetUserInfo()->GetEntries()-1;  
  //
  TString outputString="";
  outputString+="fileName.GetString();1;fname;/C:";
  outputString+="runNumber;1;run;/i:";
  outputString+="evtNumberInFile;1.20;eventID;/i:";
  outputString+="evtTimeStamp;1.20;timeStamp;/i:";
  outputString+="gid;1.30;gid;/l:";
  outputString+=TString::Format("This->GetUserInfo()->At(%d)->GetName();1;trigger;/C:",triggerIndex);
  outputString+="type;10.20;triggerType;/i:";
  outputString+="triggerClass.GetName();20.20;triggerClass;/C:";
  TString eventListName = "filteredHighPtV0s.list";
  AliTreePlayer::selectWhatWhereOrderBy(tree, outputString.Data(), cutAll.Data(),"",0,nEvents, "csv",eventListName);
  gSystem->Exec(Form("{ rm %s && uniq > %s; } < %s ",eventListName.Data(),eventListName.Data(),eventListName.Data()));
  //
  //
  Int_t nV0Selected = tree->Draw("1",cutAll.Data());
  Int_t nEventsAll=treeEvent->GetEntries();
  // Dump counters
  ::Info("KeyValue.TriggerHighPtV0s.Selection","%s",cutAll.Data());
  ::Info("KeyValue.TriggerHighPtV0s.NEventsAll","%d",nEventsAll);
  ::Info("KeyValue.TriggerHighPtV0s.NV0Selected","%d",nV0Selected);
  delete pcstream;
}


void TriggerCosmicPairs( const char * chinput,  const char * filter, Long64_t nEvents){
  /*
    chinput=gSystem->ExpandPathName("$NOTES/JIRA/PWGPP-227/data/2016/LHC16t/000267161/pass1_CENT_wSDD/filtered.list");
    filter=""; 
    Int_t nEvents=10000; 
    TriggerCosmicPairs(chinput, filter, nEvents);
  */
  const char * treeName="CosmicPairs";
  TTree *tree=NULL; 
  TTree *treeEvent=NULL; 
  if (TString(chinput).Contains(".root")){
    TFile * fInputFile = TFile::Open(chinput);
    tree               = (TTree*)fInputFile->Get(treeName);
    treeEvent          = (TTree*)fInputFile->Get("eventInfoV0");
  } else {
    tree=(TTree*)AliXRDPROOFtoolkit::MakeChain(chinput, treeName, 0, 1000000000,0,1);
    treeEvent=(TTree*)AliXRDPROOFtoolkit::MakeChain(chinput, "eventInfoV0", 0, 1000000000,0,1);
  }
  tree->SetEstimate(tree->GetEntries());
  TTreeSRedirector * pcstream = new TTreeSRedirector("filteredCosmicPairs.root","recreate");
  //
  // Additional cuts for cosmics
  TString cutAll="1";
  if (gSystem->Getenv("isCosmic")!=NULL &&strstr(gSystem->Getenv("isCosmic"),"0")!=0){
      tree->SetAlias("alphaPrime","abs(t0.fAlpha)-pi/2");
      tree->SetAlias("alphaPrimeFit","TMath::Gaus(alphaPrime,0,0.573+0)");
      tree->SetAlias("alphaPrimeDownscale","alphaPrimeFit*rndm<0.05");  
      tree->SetAlias("itsFiducial","min(abs(t0.fP[0]),abs(t1.fP[0]))<16&&min(abs(t0.fP[1]),abs(t1.fP[1]))<30");
      cutAll="itsFiducial || alphaPrimeDownscale";
  }
  else{
      TCut cutDCA="abs(0.5*(t0.fD-t1.fD))>5&&abs(0.5*(t0.fD-t1.fD))<80"; //tracks crossing the inner field cage (80cm)
      TCut cutCross="t0.fOp.fP[1]*t1.fOp.fP[1]<0"; //tracks crossing central electrode
      cutAll= cutDCA && cutCross;
  }
  cutAll=TString::Format("%s&&%s",cutAll.Data(),filter);
  //
  // Dump input parameters
  ::Info("KeyValue.TriggerCosmicPairs.Input.filter","%s",filter);
  ::Info("KeyValue.TriggerCosmicPairs.Input.nEvents","%llu",nEvents);
  //
  // Dump selected Events
  tree->GetUserInfo()->AddLast(new TObjString("CosmicPairs"));
  Int_t triggerIndex = tree->GetUserInfo()->GetEntries()-1;  
  //
  TString outputString="";
  outputString+="fileName.GetString();1;fname;/C:";
  outputString+="runNumber;1;run;/i:";
  outputString+="evtNumberInFile;1.20;eventID;/i:";
  outputString+="evtTimeStamp;1.20;timeStamp;/i:";
  outputString+="gid;1.30;gid;/l:";
  outputString+=TString::Format("This->GetUserInfo()->At(%d)->GetName();1;trigger;/C:",triggerIndex);
  //   outputString+="type;10.20;triggerType;/i:";
  outputString+="triggerClass.GetName();20.20;triggerClass;/C:";
  TString eventListName = "filteredCosmicPairs.list";
  AliTreePlayer::selectWhatWhereOrderBy(tree, outputString.Data(), cutAll.Data(),"",0,nEvents, "csv",eventListName);
  gSystem->Exec(Form("{ rm %s && uniq > %s; } < %s ",eventListName.Data(),eventListName.Data(),eventListName.Data()));
  //
  //
  Int_t nCosmicPairsSelected = tree->Draw("1",cutAll.Data());
  Int_t nEventsAll=treeEvent->GetEntries();
  // Dump counters
  ::Info("KeyValue.TriggerCosmicPairs.Selection","%s",cutAll.Data());
  ::Info("KeyValue.TriggerCosmicPairs.NEventsAll","%d",nEventsAll);
  ::Info("KeyValue.TriggerCosmicPairs.NCosmicPairsSelected","%d",nCosmicPairsSelected);
  delete pcstream;
}


void TriggerLaser( const char * chinput, Long64_t nEvents){
  /*
    chinput=gSystem->ExpandPathName("$NOTES/JIRA/PWGPP-227/data/2016/LHC16t/000267161/pass1_CENT_wSDD/filtered.list");
    Int_t nEvents=10000; 
    TriggerLaser(chinput, filter, nEvents);
  */
  const char * treeName="Laser";
  TTree *tree=NULL; 
  TTree *treeEvent=NULL; 
  if (TString(chinput).Contains(".root")){
    TFile * fInputFile = TFile::Open(chinput);
    tree               = (TTree*)fInputFile->Get(treeName);
    treeEvent          = (TTree*)fInputFile->Get("eventInfoV0");
  } else {
    tree=(TTree*)AliXRDPROOFtoolkit::MakeChain(chinput, treeName, 0, 1000000000,0,1);
    treeEvent=(TTree*)AliXRDPROOFtoolkit::MakeChain(chinput, "eventInfoV0", 0, 1000000000,0,1);
  }
  tree->SetEstimate(tree->GetEntries());
  TTreeSRedirector * pcstream = new TTreeSRedirector("filteredLaser.root","recreate");
  
  TString cutAll = "1";
  
  //
  // Dump input parameters
  ::Info("KeyValue.TriggerLaser.Input.nEvents","%llu",nEvents);
  //
  // Dump selected Events
  tree->GetUserInfo()->AddLast(new TObjString("Laser"));
  Int_t triggerIndex = tree->GetUserInfo()->GetEntries()-1;  
  //
  TString outputString="";
  //   outputString+="fileName.GetString();1;fname;/C:";
  outputString+="runNumber;1;run;/i:";
  outputString+="evtNumberInFile;1.20;eventID;/i:";
  outputString+="evtTimeStamp;1.20;timeStamp;/i:";
  outputString+="gid;1.30;gid;/l:";
  outputString+=TString::Format("This->GetUserInfo()->At(%d)->GetName();1;trigger;/C:",triggerIndex);
  outputString+="triggerClass.GetName();20.20;triggerClass;/C:";
  TString eventListName = "filteredLaser.list";
  AliTreePlayer::selectWhatWhereOrderBy(tree, outputString.Data(), cutAll.Data(),"",0,nEvents, "csv",eventListName);
  gSystem->Exec(Form("{ rm %s && uniq > %s; } < %s ",eventListName.Data(),eventListName.Data(),eventListName.Data()));

  //
  //
  Int_t nLaserSelected = tree->Draw("1",cutAll.Data());
  Int_t nEventsAll=treeEvent->GetEntries();
  // Dump counters
  ::Info("KeyValue.TriggerLaser.Selection","%s",cutAll.Data());
  ::Info("KeyValue.TriggerLaser.NEventsAll","%d",nEventsAll);
  ::Info("KeyValue.TriggerLaser.NLaserSelected","%d",nLaserSelected);
  delete pcstream;
}


void triggerCalibHighPt( const char * chinput,  const char * filter, Long64_t nEvents, Double_t dcaCut, Double_t multRatioFraction, Double_t multFraction, Int_t maxTracks){
  // Calibration events:
  //   * multiplicty cumulant cut
  //   * pilaup ratio cumulant cut
  //   * pt cumulant cut
  //
  /*
    TString chinput="filtered1.list";
    filter="(esdTrack.fFlags&0x4401)>0&&esdTrack.Pt()>2&&mult<1000"; 
    Int_t nEvents=100000; 
    Double_t dcaCut=3;    
    Double_t multRatioFraction=0.8;
    Double_t multFraction=0.5;
    Int_t maxTracks=100000;   
    triggerCalibHighPt("filtered1.list", "(esdTrack.fFlags&0x4401)>0&&esdTrack.Pt()>2&&mult<1000", 10000000, 3,0.8,0.5, 10000);
  */

  TTree *treeHighPt=NULL; 
  TTree *treeEvent=NULL; 
  if (TString(chinput).Contains(".root")){
    TFile * fInputFile = TFile::Open(chinput);
    treeHighPt               = (TTree*)fInputFile->Get("highPt");
    treeEvent          = (TTree*)fInputFile->Get("eventInfoV0");
  } else {
    treeHighPt=(TTree*)AliXRDPROOFtoolkit::MakeChain(chinput, "highPt", 0, 1000000000,0,1);
    treeEvent=(TTree*)AliXRDPROOFtoolkit::MakeChain(chinput, "eventInfoV0", 0, 1000000000,0,1);
  }
  TTreeSRedirector * pcstream = new TTreeSRedirector("filteredCalibHighPt.root","recreate");
  //
  // fill high pt calibration events
  //
  treeHighPt->SetEstimate(treeHighPt->GetEntries());
  AliAnalysisTaskFilteredTree::SetDefaultAliasesHighPt(treeHighPt);
  // Dump selected Events
  Int_t triggerIndex = treeHighPt->GetUserInfo()->GetEntries()-1;  
  treeHighPt->Draw("mult/ntracks:mult",filter,"goff",nEvents);
  Double_t fractionCut=TMath::KOrdStat(treeHighPt->GetSelectedRows(),treeHighPt->GetV1(),Long64_t(treeHighPt->GetSelectedRows()*multRatioFraction));
  Double_t multCut=TMath::KOrdStat(treeHighPt->GetSelectedRows(),treeHighPt->GetV2(),Long64_t(treeHighPt->GetSelectedRows()*multFraction));
  
  treeHighPt->SetAlias("cutHighPt",TString::Format("abs(esdTrack.fdTPC)<%.3f&&abs(esdTrack.fzTPC)<%.3f&&esdTrack.fzTPC!=0&&sqrt(esdTrack.fC[14])<0.01&&%s&&mult/ntracks>%f&&mult>%f",dcaCut,dcaCut, filter,fractionCut, multCut).Data());
  treeHighPt->Draw("-esdTrack.Pt():mult:ntracks","cutHighPt","goff",nEvents);
  Double_t meanMult=TMath::Mean(treeHighPt->GetSelectedRows(),treeHighPt->GetV2());
  Double_t mult95=TMath::KOrdStat(treeHighPt->GetSelectedRows(),treeHighPt->GetV2(),Long64_t(treeHighPt->GetSelectedRows()*0.95));
  Double_t ntracks95=TMath::KOrdStat(treeHighPt->GetSelectedRows(),treeHighPt->GetV3(),Long64_t(treeHighPt->GetSelectedRows()*0.95));
  Double_t ptQuant95=-TMath::KOrdStat(treeHighPt->GetSelectedRows(),treeHighPt->GetV1(),Long64_t(treeHighPt->GetSelectedRows()*0.05));
  Double_t ptQuant=-TMath::KOrdStat(treeHighPt->GetSelectedRows(),treeHighPt->GetV1(),TMath::Min(treeHighPt->GetSelectedRows()-1,Long64_t(maxTracks/meanMult)-1));
  treeHighPt->SetAlias("cutHighPt",TString::Format("abs(esdTrack.fdTPC)<%.3f&&abs(esdTrack.fzTPC)<%.3f&&esdTrack.fzTPC!=0&&sqrt(esdTrack.fC[14])<0.01&&%s&&mult/ntracks>%f&&esdTrack.Pt()>%f&&mult>%f",dcaCut,dcaCut, filter,fractionCut,ptQuant,multCut).Data());
  //
  TString outputString="";
  outputString+="fileName.GetString();1;fname;/C:";
  outputString+="runNumber;1;run;/i:";
  outputString+="evtNumberInFile;1.20;eventID;/i:";
  outputString+="evtTimeStamp;1.20;timeStamp;/i:";
  outputString+="gid;1.30;gid;/l:";
  outputString+="esdTrack.Pt();1.30;trackPt;/l:";
  outputString+="mult;1.30;mult;/l:";
  outputString+=TString::Format("This->GetUserInfo()->At(%d)->GetName();1;trigger;/C:",triggerIndex);
  outputString+="triggerClass.GetName();20.20;triggerClass;/C:";
  TString eventListName = "filteredCalibHighPt.list";
  AliTreePlayer::selectWhatWhereOrderBy(treeHighPt, outputString.Data(), "cutHighPt","",0,nEvents, "csv",eventListName);  
  gSystem->Exec(Form("{ rm %s && uniq > %s; } < %s ",eventListName.Data(),eventListName.Data(),eventListName.Data()));
  //
  // QA plots
  TString qaHistograms="";
  qaHistograms+=TString::Format("mult:#1>>hisMultAll(%.0f,0,%.0f);",TMath::Min(2*mult95,200.),2*mult95);
  qaHistograms+=TString::Format("mult:#cutHighPt>>hisMultSelected(%.0f,0,%.0f);",TMath::Min(2*mult95,200.),2*mult95);
  qaHistograms+=TString::Format("ntracks:#1>>hisNtracksAll(%.0f,0,%.0f);",TMath::Min(ntracks95,200.),2*ntracks95);
  qaHistograms+=TString::Format("ntracks:#cutHighPt>>hisNtracksSelected(%.0f,0,%.0f);",TMath::Min(ntracks95,200.),2*ntracks95);
  qaHistograms+=TString::Format("mult/ntracks:#1>>hisMultNtracksAll(100,0,1);");
  qaHistograms+=TString::Format("mult/ntracks:#cutHighPt>>hisMultNtracksSelected(100,0,1);");
  qaHistograms+=TString::Format("esdTrack.Pt():#1>>hisPtAll(100,0,%.0f);",2.*ptQuant95);
  qaHistograms+=TString::Format("esdTrack.Pt():#cutHighPt>>hisPtSelected(100,0,%.0f);",2.*ptQuant95);
  TStopwatch timer; TObjArray *hisArray = AliTreePlayer::MakeHistograms(treeHighPt, qaHistograms, filter,0,nEvents,nEvents,15); timer.Print();
  //
  TString drawExpression="";
  drawExpression="[1,1,1,1]:";
  drawExpression+="%Ogridx,gridy;hisMultAll(0,200)(0)(err);hisMultSelected(0,200)(0)(err):";
  drawExpression+="%Ogridx,gridy;hisMultNtracksAll(0,10000)(0)(err);hisMultNtracksSelected(0,10000)(0)(err):";
  drawExpression+="%Ogridx,gridy;hisNtracksAll(0,10000)(0)(err);hisNtracksSelected(0,10000)(0)(err):";
  drawExpression+="%Ogridx,gridy;hisPtAll(0,100)(0)(err);hisPtSelected(0,100)(0)(err):";
  TObjArray *keepArray = new TObjArray(100);
  TPad * pad = AliTreePlayer::DrawHistograms(0,hisArray,drawExpression,keepArray, 15); 
  ((TCanvas*)pad)->SetWindowSize(1800,1000);
  pad->SaveAs("triggerCalibHighPt.png");
  pad->SaveAs("triggerCalibHighPt.C");
  pcstream->GetFile()->cd();
  hisArray->Write("hisArray",TObjArray::kSingleKey);
  keepArray->Write("keepArray",TObjArray::kSingleKey);
  //
  //
  Int_t nHighPtSelected = treeHighPt->Draw("esdTrack.fTPCncls","cutHighPt","goff",nEvents);
  Int_t nEventsAll=treeEvent->GetEntries();
  // Dump input parameters
  ::Info("SKeyValue.TriggerCalibHighPt.Input.Filter","%s",filter);
  ::Info("KeyValue.TriggerCalibHighPt.Input.dcaCut","%f",dcaCut);
  ::Info("KeyValue.TriggerCalibHighPt.Input.nEvents","%llu",nEvents);
  //
  ::Info("SKeyValue.TriggerCalibHighPt.Cut.Selection","%s",treeHighPt->GetAlias("cutHighPt"));
  ::Info("KeyValue.TriggerCalibHighPt.Cut.multCut","%f",multCut);
  ::Info("KeyValue.TriggerCalibHighPt.Cut.fractionCut","%f",fractionCut);
  ::Info("KeyValue.TriggerCalibHighPt.Cut.ptQuant","%f",ptQuant);
  ::Info("KeyValue.TriggerCalibHighPt.Cut.meanMult","%f",meanMult);
  // Dump counters
  ::Info("KeyValue.TriggerCalibHighPt.NEventsAll","%d",nEventsAll);
  ::Info("KeyValue.TriggerCalibHighPt.NHighPtSelected","%d",nHighPtSelected);
  delete pcstream;
}


void triggerCalibV0( const char * chinput,  const char * filter, Long64_t nEvents, Double_t multRatioFraction, Double_t multFraction, Int_t maxTracks, Double_t zCut, Double_t covarQPtCut, Int_t filterMask){
  // Calibration events:
  //   * multiplicty cumulant cut
  //   * pilaup ratio cumulant cut
  //   * pt cumulant cut
  //
  /*
    TString chinput="filtered1.list";
    filter="((type==1)*max(track0.Pt(),track1.Pt())>1)||(max(track0.Pt(),track1.Pt())>2)"; 
    Int_t nEvents=100000; 
    Double_t dcaCut=3;    
    Double_t multRatioFraction=0.8;
    Double_t multFraction=0.5;
    Int_t maxTracks=100000;   
    Double_t zCut=20.;
    Double_t covarQPtCut=0.03;
    filterMask=0x4401;

  */

  TTree *treeV0=NULL; 
  TTree *treeEvent=NULL; 
  if (TString(chinput).Contains(".root")){
    TFile * fInputFile = TFile::Open(chinput);
    treeV0               = (TTree*)fInputFile->Get("V0s");
    treeEvent          = (TTree*)fInputFile->Get("eventInfoV0");
  } else {
    treeV0=(TTree*)AliXRDPROOFtoolkit::MakeChain(chinput, "V0s", 0, 1000000000,0,1);
    treeEvent=(TTree*)AliXRDPROOFtoolkit::MakeChain(chinput, "eventInfoV0", 0, 1000000000,0,1);
  }
  treeEvent->BuildIndex("gid");
  treeV0->BuildIndex("gid");
  treeV0->AddFriend(treeEvent,"E.");
  treeV0->SetAlias("filterMask",TString::Format("((track0.fFlags&%x)+(track1.fFlags&%x)>0)",filterMask,filterMask).Data());
  treeV0->SetAlias("dcaZcut",TString::Format("abs(track0.fzTPC)<%.3f&&abs(track1.fzTPC)<%.3f",zCut,zCut).Data());
  treeV0->SetAlias("covarQPtCut",TString::Format("sqrt(track0.fC[14])<%.3f&&sqrt(track1.fC[14])<%0.3f",covarQPtCut,covarQPtCut).Data());
  //
  TTreeSRedirector * pcstream = new TTreeSRedirector("filteredCalibV0.root","recreate");
  //
  // fill high pt calibration events
  //
  treeV0->SetEstimate(treeV0->GetEntries()); 
  AliAnalysisTaskFilteredTree::SetDefaultAliasesV0(treeV0);
  // Dump selected Events
  Int_t triggerIndex = treeV0->GetUserInfo()->GetEntries()-1;  
  treeV0->Draw("mult/ntracks:mult",filter,"goff",nEvents);
  Double_t fractionCut=TMath::KOrdStat(treeV0->GetSelectedRows(),treeV0->GetV1(),Long64_t(treeV0->GetSelectedRows()*multRatioFraction));
  Double_t multCut=TMath::KOrdStat(treeV0->GetSelectedRows(),treeV0->GetV2(),Long64_t(treeV0->GetSelectedRows()*multFraction));
  //
  treeV0->SetAlias("cutHighPt",TString::Format("filterMask&&dcaZcut&&covarQPtCut&&%s&&mult/ntracks>%f&&mult>%f",filter,fractionCut, multCut).Data());
  treeV0->Draw("-v0.Pt():mult:ntracks","cutHighPt","goff",nEvents);
  Double_t meanMult=TMath::Mean(treeV0->GetSelectedRows(),treeV0->GetV2());
  Double_t mult95=TMath::KOrdStat(treeV0->GetSelectedRows(),treeV0->GetV2(),Long64_t(treeV0->GetSelectedRows()*0.95));
  Double_t ntracks95=TMath::KOrdStat(treeV0->GetSelectedRows(),treeV0->GetV3(),Long64_t(treeV0->GetSelectedRows()*0.95));
  Double_t ptQuant95=-TMath::KOrdStat(treeV0->GetSelectedRows(),treeV0->GetV1(),Long64_t(treeV0->GetSelectedRows()*0.05));
  Double_t ptQuant=-TMath::KOrdStat(treeV0->GetSelectedRows(),treeV0->GetV1(),TMath::Min(treeV0->GetSelectedRows()-1,Long64_t(maxTracks/meanMult))-1);
  treeV0->SetAlias("cutHighPt",TString::Format("filterMask&&dcaZcut&&covarQPtCut&&%s&&mult/ntracks>%f&&v0.Pt()>%f&&mult>%f",filter,fractionCut,ptQuant,multCut).Data());
  //
  TString outputString="";
  outputString+="fileName.GetString();1;fname;/C:";
  outputString+="runNumber;1;run;/i:";
  outputString+="evtNumberInFile;1.20;eventID;/i:";
  outputString+="evtTimeStamp;1.20;timeStamp;/i:";
  outputString+="gid;1.30;gid;/l:";
  outputString+="v0.Pt();1.30;trackPt;/l:";
  outputString+="mult;1.30;mult;/l:";
  outputString+=TString::Format("This->GetUserInfo()->At(%d)->GetName();1;trigger;/C:",triggerIndex);
  outputString+="triggerClass.GetName();20.20;triggerClass;/C:";
  TString eventListName = "filteredCalibV0.list";
  AliTreePlayer::selectWhatWhereOrderBy(treeV0, outputString.Data(), "cutHighPt","",0,nEvents, "csv",eventListName);  
  gSystem->Exec(Form("{ rm %s && uniq > %s; } < %s ",eventListName.Data(),eventListName.Data(),eventListName.Data()));
  //
  // QA plots
  TString qaHistograms="";
  qaHistograms+=TString::Format("mult:#1>>hisMultAll(%.0f,0,%.0f);",TMath::Min(2*mult95,200.),2*mult95);
  qaHistograms+=TString::Format("mult:#cutHighPt>>hisMultSelected(%.0f,0,%.0f);",TMath::Min(2*mult95,200.),2*mult95);
  qaHistograms+=TString::Format("ntracks:#1>>hisNtracksAll(%.0f,0,%.0f);",TMath::Min(ntracks95,200.),2*ntracks95);
  qaHistograms+=TString::Format("ntracks:#cutHighPt>>hisNtracksSelected(%.0f,0,%.0f);",TMath::Min(ntracks95,200.),2*ntracks95);
  qaHistograms+=TString::Format("mult/ntracks:#1>>hisMultNtracksAll(100,0,1);");
  qaHistograms+=TString::Format("mult/ntracks:#cutHighPt>>hisMultNtracksSelected(100,0,1);");
  qaHistograms+=TString::Format("v0.Pt():#1>>hisPtAll(100,0,%.0f);",2.*ptQuant95);
  qaHistograms+=TString::Format("v0.Pt():#cutHighPt>>hisPtSelected(100,0,%.0f);",2.*ptQuant95);
  TStopwatch timer; TObjArray *hisArray = AliTreePlayer::MakeHistograms(treeV0, qaHistograms, filter,0,nEvents,nEvents,15); timer.Print();
  //
  TString drawExpression="";
  drawExpression="[1,1,1,1]:";
  drawExpression+="%Ogridx,gridy;hisMultAll(0,200)(0)(err);hisMultSelected(0,200)(0)(err):";
  drawExpression+="%Ogridx,gridy;hisMultNtracksAll(0,10000)(0)(err);hisMultNtracksSelected(0,10000)(0)(err):";
  drawExpression+="%Ogridx,gridy;hisNtracksAll(0,10000)(0)(err);hisNtracksSelected(0,10000)(0)(err):";
  drawExpression+="%Ogridx,gridy;hisPtAll(0,100)(0)(err);hisPtSelected(0,100)(0)(err):";
  TObjArray *keepArray = new TObjArray(100);
  TPad * pad = AliTreePlayer::DrawHistograms(0,hisArray,drawExpression,keepArray, 15); 
  ((TCanvas*)pad)->SetWindowSize(1400,800);
  pad->SaveAs("triggerCalibV0.png");
  pad->SaveAs("triggerCalibV0.C");
  pcstream->GetFile()->cd();
  hisArray->Write("hisArray",TObjArray::kSingleKey);
  keepArray->Write("keepArray",TObjArray::kSingleKey);
  //
  //
  Int_t nHighPtSelected = treeV0->Draw("track0.fTPCncls","cutHighPt","goff",nEvents);
  Int_t nEventsAll=treeEvent->GetEntries();
  // Dump input parameters
  ::Info("KeyValue.TriggerCalibV0.Input.filter","%s",filter);
  ::Info("KeyValue.TriggerCalibV0.Input.nEvents","%llu",nEvents);
  ::Info("KeyValue.TriggerCalibV0.Input.zCut","%f",zCut);
  ::Info("KeyValue.TriggerCalibV0.Input.covarQPtCut","%f",covarQPtCut);
  ::Info("KeyValue.TriggerCalibV0.Input.filterMask","%x",filterMask);
  //
  ::Info("SKeyValue.TriggerCalibV0.Cut.Selection","%s",treeV0->GetAlias("cutHighPt"));
  ::Info("KeyValue.TriggerCalibV0.Cut.multCut","%f",multCut);
  ::Info("KeyValue.TriggerCalibV0.Cut.fractionCut","%f",fractionCut);
  ::Info("KeyValue.TriggerCalibV0.Cut.ptQuant","%f",ptQuant);
  ::Info("KeyValue.TriggerCalibV0.Cut.meanMult","%f",meanMult);
  // Dump counters
  ::Info("KeyValue.TriggerCalibV0.NEventsAll","%d",nEventsAll);
  ::Info("KeyValue.TriggerCalibV0.NV0Selected","%d",nHighPtSelected);
  delete pcstream;
}


void GetRawSummary(){
  //
  // Compare content of the raw data files with the raw gids
  // Input:
  //   gidraw.tree  - raw 
  //   ref.list         - list of files with refernce triggers 
  // Assumption - all files are in local directory
  // 
  // Output: 
  //   log file with KeyValue of events summary
  TObjArray *refArray=(gSystem->GetFromPipe("cat ref.list")).Tokenize("\n");
  TTree  *treeRaw = new TTree;
  treeRaw->ReadFile("gidraw.tree","",'\t');
  treeRaw->BuildIndex("gid");
  ::Info("GetRawSummary.KeyValue.RawAll","%d",Int_t(treeRaw->GetEntries()));
  for  (Int_t iRef=0; iRef<refArray->GetEntries(); iRef++){
    TTree * tree = new TTree;
    tree->ReadFile(refArray->At(iRef)->GetName());
    TString tName = TString(refArray->At(iRef)->GetName()).ReplaceAll(".list","").ReplaceAll("filtered","");
    tree->SetName(tName);
    tree->BuildIndex("gid");
    tree->AddFriend(treeRaw,"RAW");    
    treeRaw->AddFriend(tree,TString::Format("T%d",iRef).Data());
    Int_t nFiltered=tree->GetEntries();
    Int_t nRaws0=treeRaw->Draw("gid%10000",TString::Format("gid==T%d.gid",iRef),"goff");
    Int_t nRaws1=tree->Draw("gid%10000","gid==RAW.gid","goff");
    ::Info(TString::Format("GetRawSummary.KeyValue.Filter%sAll",tName.Data()),"%d", nFiltered);
    ::Info(TString::Format("GetRawSummary.KeyValue.Filter%sRaw0",tName.Data()),"%d",nRaws0);
    ::Info(TString::Format("GetRawSummary.KeyValue.Filter%sRaw1",tName.Data()),"%d",nRaws1);
    ::Info(TString::Format("GetRawSummary.KeyValue.Filter%sRatio",tName.Data()),"%f",nRaws1/Double_t(nFiltered));
  }
}

void SummarizeLogs(const char * logTreeFile="log.tree", const char * pass=0){
  //
  // Input data assumed to be in fromatted csv file with header
  //     year/d:period/C:run/d:name/C:key/C:value/I
  // Input parsed log files  obtained using script - where the logPath is the prefix path to the log directories  
  // ( source $AliPhysics_SRC/PWGPP/rawmerge/makeOfflineTriggerList.sh; ProcessfilterLog $logPath )

  //
  TTree logTree;
  logTree.ReadFile(logTreeFile,"",'\t');  
  if (logTree.GetEntries()<=0){
    ::Error("makeOfflineTriggerList.SummarizeLogs","Invalid input file %s. No entries. Exiting",logTreeFile);
    return;
  }
  char cperiod[1000];
  TBranch * br = logTree.GetBranch("period");
  if (br==NULL){
    ::Error("makeOfflineTriggerList.SummarizeLogs","Invalid input file %s. Missing branch period. Exiting",logTreeFile);
    return;
  }
  br->SetAddress(cperiod);
  br->GetEntry(0);
  AliExternalInfo info; info.fVerbose=0;
  TTree *logbookTree =info.GetTree("Logbook",cperiod,"","Logbook.detector:TPC:detector==\"TPC\"");
  logTree.BuildIndex("run");
  logTree.AddFriend(logbookTree,"Logbook");
  logTree.AddFriend(logbookTree->GetFriend("Logbook.detector_TPC"),"LogbookTPC");
  TTree *evsTree=NULL;
  if (pass!=NULL){
    evsTree =info.GetTree("QA.EVS",cperiod,pass,"");
    evsTree->BuildIndex("run");
    logTree.AddFriend(evsTree,"EVS");
  }
  
  //
  TGraph *grLogbook= TStatToolkit::MakeGraphSparse(&logTree, "Logbook.totalEventsPhysics:run","Logbook.run==run",25,1,1);
  TGraph *grLogbookTPC= TStatToolkit::MakeGraphSparse(&logTree, "LogbookTPC.eventCountPhysics:run","run==LogbookTPC.run",21,2,1);
  TGraph *grAll    = TStatToolkit::MakeGraphSparse(&logTree, "value:run","(strstr(key,\"TriggerHighPt.NEventsAll\")!=0)",21,4,1);
  TGraph *grHighPt = TStatToolkit::MakeGraphSparse(&logTree, "100*value:run","(strstr(key,\"PtSelected\")!=0)",22,6,1);
  TGraph *grV0s    = TStatToolkit::MakeGraphSparse(&logTree, "100*value:run","(strstr(key,\"V0Selected\")!=0)",23,7,1);
  TGraph *grMult   = TStatToolkit::MakeGraphSparse(&logTree, "100*value:run","(strstr(key,\"HighMultiplicitySelected\")!=0)",24,8,1);
  TGraph * graphs[6]={grLogbook,grLogbookTPC,  grAll, grHighPt, grV0s, grMult};
  Double_t minY=grAll->GetY()[0], maxY=grAll->GetY()[0];

  for (Int_t igr=0; igr<6; igr++){
    if ( graphs[igr]==NULL) continue;
    Double_t gmin=TMath::MinElement(graphs[igr]->GetN(),graphs[igr]->GetY());
    Double_t gmax=TMath::MaxElement(graphs[igr]->GetN(),graphs[igr]->GetY());
    if (minY>gmin) minY=gmin;
    if (maxY<gmax) maxY=gmax;
  }
  grLogbook->SetMaximum(maxY+(maxY-minY)*0.5);
  grLogbook->SetMinimum(0);
  //
  TCanvas *canvasStat= new TCanvas("canvasStat","canvasStat",1400,600);
  canvasStat->SetRightMargin(0.05);
  TLegend *legend = new TLegend(0.11,0.70,0.6,0.88,TString::Format("Offline trigger counters period %s  (makeOfflineTriggerList:SummarizeLogs)",cperiod).Data());
  legend->SetNColumns(3);
  legend->SetBorderSize(0);
  grLogbook->GetYaxis()->SetTitle("# Events");
  grLogbook->GetXaxis()->SetTitle("run");
  for (Int_t igr=0; igr<6; igr++){
    if ( graphs[igr]==NULL) continue;
    if (igr==0) graphs[igr]->Draw("alp");
    graphs[igr]->Draw("lp");
  }
  legend->AddEntry(grLogbook,"Logbook #totalEventsPhysics","p");
  legend->AddEntry(grLogbookTPC,"TPC #totalEventsPhysics","p");
  legend->AddEntry(grAll,"ESD. #Events processed ","p");
  legend->AddEntry(grHighPt,"# High pt x 100","p");
  legend->AddEntry(grV0s,"# V0s x 100","p");
  if (grMult) legend->AddEntry(grMult,"#High mult x 100","p");
  //
  legend->Draw();
  canvasStat->SaveAs("makeOfflineTriggerListEventSummary.png");
  canvasStat->SaveAs("makeOfflineTriggerListEventSummary.pdf");
  canvasStat->SaveAs("makeOfflineTriggerListEventSummary.C");
  // make period counter summary
  gStyle->SetOptStat(0);
  logTree.SetMarkerStyle(21);
  logTree.Draw("key","( (strstr(key,\"TriggerHighPt.NEventsAll\")>0)||(strstr(key,\"elected\")>0))*value/1000000.");
  logTree.GetHistogram()->GetXaxis()->SetTitle("Trigger");
  logTree.GetHistogram()->GetYaxis()->SetTitle("#Events/1000000");
  logTree.GetHistogram()->Draw("hist TEXT0");
  canvasStat->SaveAs("numberOfEventsTriggered.png");
  canvasStat->SaveAs("numberOfEventsTriggered.C");
  // make period fraction summary
  logTree.GetHistogram()->Scale(1/logTree.GetHistogram()->GetBinContent(1));
  logTree.GetHistogram()->GetXaxis()->SetRangeUser(1,10000);
  logTree.GetHistogram()->GetYaxis()->SetTitle("fraction");
  logTree.GetHistogram()->Draw("hist TEXT0");
  canvasStat->SaveAs("fractionOfEventsTriggered.png");
  canvasStat->SaveAs("fractionOfEventsTriggered.C");

}


void checkReconstuctedTrigger(){
  /*
    TFile *f = TFile::Open("/lustre/nyx/alice/users/marsland/alice-tpc-notes/JIRA/PWGPP-126/filtering/2016/LHC16t/benchmark/test2/2016/LHC16t/000267165/cpass0/211_rawSelected0/AliESDs.root");
    TTree * esdTree= (TTree*)f->Get("esdTree");
  */
  TTree* esdTree = AliXRDPROOFtoolkit::MakeChain("esd.list","esdTree",0, 2000);
  esdTree->SetAlias("gid","((AliESDHeader.fPeriodNumber<<36)+(AliESDHeader.fOrbitNumber<<12)+AliESDHeader.fBunchCrossNumber)"); 
  esdTree->BuildIndex("gid");
  esdTree->Scan("gid%1000","","",3);
  //
  //
  //
  TTree treeHighPt,treeHighPtV0,treeHighMilt;
  treeHighPt.ReadFile("filteredHighPt.list","",'\t');
  treeHighPt.BuildIndex("gid");
  treeHighPt.AddFriend(esdTree,"E");
  esdTree->Scan("gid%1000","","",3);
  treeHighPt.Draw("Tracks[].Pt()>>hisPt(100, 2,30)","gid==E.gid&&Tracks[].fFlags&&0x44==0x44&&Tracks[].fTPCncls>100","",10000);
  gPad->SaveAs("TracksPt.png");
  gPad->SaveAs("TracksPt.C");
  //
  // High mult
  //
  TTree treeHighMult;
  treeHighMult.ReadFile("filteredMult.list","",'\t');
  treeHighMult.BuildIndex("gid");
  treeHighMult.AddFriend(esdTree,"E");  
  treeHighMult.Draw("E.SPDVertex.fNContributors>>hisContibutors","gid==E.gid","");
  gPad->SaveAs("SPDVertexfNContributors.png");

  
  
}
