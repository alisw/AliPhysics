/*
  .L $ALICE_PHYSICS/../src/PWGPP/rawmerge/makeOfflineTriggerList.C+
  makeOfflineTriggerList(gSystem->ExpandPathName("$NOTES/JIRA/PWGPP-227/data/2016/LHC16t/000267161/pass1_CENT_wSDD/filtered.list"))
  //
  // makeOfflineTriggerList for selected filtered trees
  // configuration parameters are taken fro env variables  defined in makeOfflineTriggerListSetupXXX.sh
  // 
   
  aliroot -b -q $ALICE_PHYSICS/../src/PWGPP/rawmerge/makeOfflineTriggerList.C+\(\"$NOTES/JIRA/PWGPP-227/data/2016/LHC16t/000267161/pass1_CENT_wSDD/filtered.list\"\) >makeOfflineTriggerList.log

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
using namespace std;
using std::cout;
using std::setw;



void TriggerHighMultiplicity( const char * chinput,  const char * filter="ntracks<1000", Long_t ntracksPrim=200000, Float_t fractionTracks=0.8, Long_t nEvents=200000);
void TriggerHighPt( const char * chinput,  const char * filter="(esdTrack.fFlags&0x4401)>0", Long_t nEvents=200000, Double_t dcaCut=4);
void TriggerHighPtV0( const char * chinput,  const char * filter, Long_t nEvents, Double_t zCut, Double_t covarQPtCut, Int_t filterMask);


void makeOfflineTriggerList(const char * chinput){  
  //
  // 
  // Make high mulitplicity event selection for the calibration
  TString highMultiplicityFilter="ntracks<1000";  
  Long_t  highMultiplicityNTracksPrim=500000;   // default 500000 tracks for the calibration
  Float_t highMultiplicityFractionTracks=0.8;   //
  Long_t  highMultiplicityNEvents=10000000;
  if (gSystem->Getenv("highMultiplicityFilter"))        highMultiplicityFilter=gSystem->Getenv("highMultiplicityFilter");
  if (gSystem->Getenv("highMultiplicityNTracksPrim"))   highMultiplicityNTracksPrim=TString(gSystem->Getenv("highMultiplicityNTracksPrim")).Atof();
  if (gSystem->Getenv("highMultiplicityFractionTracks"))highMultiplicityFractionTracks=TString(gSystem->Getenv("highMultiplicityFractionTracks")).Atof();
  if (gSystem->Getenv("highMultiplicityNEvents"))       highMultiplicityNEvents=TString(gSystem->Getenv("highMultiplicityNEvents")).Atof();
  TriggerHighMultiplicity(chinput, highMultiplicityFilter.Data(),highMultiplicityNTracksPrim,  highMultiplicityFractionTracks, highMultiplicityNEvents);
  // High pt filters
  TString highPtFilter="(esdTrack.fFlags&0x4401)>0&&esdTrack.Pt()>8";
  Long_t   highPtNEvents=1000000000;
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
  TriggerHighPtV0(chinput, highPtV0Filter.Data(),  highPtV0NEvents, highPtV0ZCut,  highPtV0CovarQPtCut,  highPtV0FilterMask);  
}

void TriggerHighMultiplicity( const char * chinput,  const char * filter, Long_t ntracksPrim, Float_t fractionTracks, Long_t nEvents){
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
    tree=(TTree*)AliXRDPROOFtoolkit::MakeChain(chinput, treeName, 0, 1000000000,0);
  }
  tree->SetEstimate(tree->GetEntries());
  TTreeSRedirector * pcstream = new TTreeSRedirector("filteredMultQA.root","recreate");
  //
  //  step1 - get cumulat fractionTracks from mult/ntracks 
  Long_t events = tree->Draw("mult/ntracks",filter,"",nEvents);
  Float_t fractionTracksCut = TMath::KOrdStat(events, tree->GetV1(), Long_t(Double_t(events)*fractionTracks));
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
  ::Info("KeyValue.TriggerHighMultiplicity.Input.ntracksPrim","%ld",ntracksPrim);
  ::Info("KeyValue.TriggerHighMultiplicity.Input.fractionTracks","%f",fractionTracks);
  ::Info("KeyValue.TriggerHighMultiplicity.Input.nEvent","%ld",nEvents);
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




void TriggerHighPt( const char * chinput,  const char * filter, Long_t nEvents, Double_t dcaCut){
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
    tree=(TTree*)AliXRDPROOFtoolkit::MakeChain(chinput, treeName, 0, 1000000000,0);
    treeEvent=(TTree*)AliXRDPROOFtoolkit::MakeChain(chinput, "eventInfoV0", 0, 1000000000,0);
  }

  tree->SetEstimate(tree->GetEntries());
  TTreeSRedirector * pcstream = new TTreeSRedirector("filteredHighPtQA.root","recreate");
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
  Int_t nHighPtSelected = tree->Draw("esdTrack.fTPCncls","cutHighPt");
  Int_t nEventsAll=treeEvent->GetEntries();
   // Dump input parameters
  
  ::Info("KeyValue.TriggerHighPt.Input.Filter","%s",filter);
  ::Info("KeyValue.TriggerHighPt.Input.dcaCut","%f",dcaCut);
  ::Info("KeyValue.TriggerHighPt.Input.nEvents","%ld",nEvents);
  // Dump counters
  ::Info("KeyValue.TriggerHighPt.Selection","%s",tree->GetAlias("cutHighPt"));
  ::Info("KeyValue.TriggerHighPt.NEventsAll","%d",nEventsAll);
  ::Info("KeyValue.TriggerHighPt.NHighPtSelected","%d",nHighPtSelected);
  delete pcstream;
}


void TriggerHighPtV0( const char * chinput,  const char * filter, Long_t nEvents, Double_t zCut, Double_t covarQPtCut, Int_t filterMask){
  /*
    chinput=gSystem->ExpandPathName("$NOTES/JIRA/PWGPP-227/data/2016/LHC16t/000267161/pass1_CENT_wSDD/filtered.list");
    filter="((type==1)*max(track0.Pt(),track1.Pt())>2)||(max(track0.Pt(),track1.Pt())>4)"; 
    Int_t nEvents=10000; 
    Double_t zCut=20.;
    Double_t covarQPtCut=0.03;
    filterMask=0x4401;
    TriggerHighPtV0(chinput, filter, nEvents, zCut,covarQPtCut,filterMask);

  */
  const char * treeName="V0s";
  TTree *tree=NULL; 
  TTree *treeEvent=NULL; 
  if (TString(chinput).Contains(".root")){
    TFile * fInputFile = TFile::Open(chinput);
    tree               = (TTree*)fInputFile->Get(treeName);
    treeEvent          = (TTree*)fInputFile->Get("eventInfoV0");
  } else {
    tree=(TTree*)AliXRDPROOFtoolkit::MakeChain(chinput, treeName, 0, 1000000000,0);
    treeEvent=(TTree*)AliXRDPROOFtoolkit::MakeChain(chinput, "eventInfoV0", 0, 1000000000,0);
  }
  tree->SetEstimate(tree->GetEntries());
  TTreeSRedirector * pcstream = new TTreeSRedirector("filteredHighPtV0.root","recreate");
  //
  tree->SetAlias("filterMask",TString::Format("((track0.fFlags&%x)+(track1.fFlags&%x)>0)",filterMask,filterMask).Data());
  tree->SetAlias("dcaZcut",TString::Format("abs(track0.fzTPC)<%.3f&&abs(track1.fzTPC)<%.3f",zCut,zCut).Data());
  tree->SetAlias("covarQPtCut",TString::Format("sqrt(track0.fC[14])<%.3f&&sqrt(track1.fC[14])<%0.3f",covarQPtCut,covarQPtCut).Data());
  tree->SetAlias("userCut",filter);
  TString cutAll=TString::Format("filterMask&&dcaZcut&&covarQPtCut&&%s",filter);
  // Dump input parameters
  ::Info("KeyValue.TriggerHighPtV0s.Input.filter","%s",filter);
  ::Info("KeyValue.TriggerHighPtV0s.Input.nEvents","%ld",nEvents);
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
  AliTreePlayer::selectWhatWhereOrderBy(tree, outputString.Data(), cutAll.Data(),"",0,nEvents, "csv","filteredHighPtV0s.list");
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


void SummarizeLogs(const char * logTreeFile="log.tree"){
  //
  // Input parsed log files 
  // obtained as        egrep KeyValue 0002*/lists/*.log | sed s_":I-KeyValue."_"\t"_ >  log.tree
  //
  TTree logTree;
  logTree.ReadFile(logTreeFile,"name/C:key/C:value/d",'\t');  
  //
  logTree->SetAlias("nEvents","value");
  logTree.SetMarkerStyle(21);
  logTree.SetMarkerColor(1);
  logTree.Draw("value:name","(strstr(key,\"TriggerHighPt.NEventsAll\")!=0)");
  gPad->SetBottomeMargin(0.25);
  logTree.SetMarkerColor(2);
  logTree.Draw("value*100:name","(strstr(key,\"PtSelected\")!=0)","same");
  logTree.SetMarkerColor(4);
  logTree.Draw("value*100:name","(strstr(key,\"V0Selected\")!=0)","same");
  //
  //Int_t sumAll=

}
