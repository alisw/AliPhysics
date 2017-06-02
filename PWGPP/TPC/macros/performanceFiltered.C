/*
  gSystem->AddIncludePath("-I$AliPhysics_SRC/PWGPP/ -I$AliPhysics_SRC/OADB/");  // ? why not in the alienv,  why not available ?

  AliDrawStyle::SetDefaults()
  AliDrawStyle::PrintStyles(0,TPRegexp("."));
  AliDrawStyle::ApplyStyle("figTemplate2");

  .L $AliPhysics_SRC/PWGPP/TPC/macros/performanceFiltered.C+ 
  gStyle->SetOptStat(0);
  //
  //
  //  performanceFiltered(20000000);
  //  InitAnalysis();
  //  AnalyzeHistograms();
  //    AnalyzeHistograms()
  //  MakeResidualDistortionMaps()
  Combined tracking performance 
  //
  aliroot -b -q $HOME/rootlogon.C $AliPhysics_SRC/PWGPP/TPC/macros/performanceFiltered.C+\(2000000\)
//
  aliroot -b -q $HOME/rootlogon.C $AliPhysics_SRC/PWGPP/TPC/macros/performanceFiltered.C+\(2000000,1\) 

*/
#include "TSystem.h"
#include "TChain.h"
#include "TProof.h"
#include "TFile.h"
#include "TCut.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TString.h"
#include "TH2.h"
#include "THn.h"
#include "stdio.h"
#include "TVectorF.h"
#include "Riostream.h"
#include <iostream>
#include "TPRegexp.h"
#include "TTreeStream.h"
#include "AliLumiTools.h"
#include "AliXRDPROOFtoolkit.h"
#include "AliExternalInfo.h"
#include "TStatToolkit.h"
#include "TLatex.h"
#include "TLegend.h"
#include "AliTreePlayer.h"
#include "TStyle.h"
#include "TKey.h"
#include "AliDrawStyle.h"
#include "AliAnalysisTaskFilteredTree.h"
//
TChain * chain=0;
TChain * chainV0=0;
TTreeSRedirector *pcstream = 0;
TObjArray fitSlicesArray(3);
TObject * toStore =0;
TGraph *lumiGraph=0;
TObjArray * hisArray=0;
TObjArray * hisArrayV0=0;
TObjArray * keepArray=0;
// Ranges
Double_t timeStart,timeEnd;
Int_t    timeBins=1;
Int_t    ntracksEnd, multEnd;
Int_t    chainEntries=0;
//   parameters
Int_t run=246272;
TString period="LHC15o";
Double_t deltaT=300;  // 5 minutes binning


void InitAnalysis();
TObjArray * FillPerfomanceHisto(Int_t maxEntries);
void SetMetadata();
void MakeResidualDistortionMaps();
//
void GetNclReport(TObjArray * hisArray,  TObjArray *keepArray );
void GetDCAReport(TObjArray * hisArray,  TObjArray *keepArray );
void GetChi2Report();
void makeP2Report();
void makeP4Report();


void performanceFiltered(Int_t maxEvents, Int_t action=0){
  //
  //   .L $NOTES/JIRA/PWGPP-221/code/performanceFiltered.C+
  //
  if (action==1) {
    MakeResidualDistortionMaps();
    return;
  }
  InitAnalysis();
  if (chain==NULL) {
    ::Error("performanceFiltered","Empty input chain");
    return;
  }
  hisArray = FillPerfomanceHisto(maxEvents);
  keepArray=new TObjArray();
  //GetNclReport(hisArray,keepArray);
  //GetDCAReport(hisArray,keepArray);
  keepArray->Write("keepArray",TObjArray::kSingleKey);
  (*pcstream)<<"perf"<<"\n";
  delete pcstream;
}


void AnalyzeHistograms(){
  //
  //
  //
  TFile *finput = TFile::Open("performanceHisto.root","read");
  hisArray=new TObjArray();
  TList * keys = finput->GetListOfKeys();
  for (Int_t iKey=0; iKey<keys->GetEntries(); iKey++){    
    TObject * o = finput->Get(TString::Format("%s;%d",keys->At(iKey)->GetName(),((TKey*)keys->At(iKey))->GetCycle()).Data());
    hisArray->AddLast(o);
  }

  pcstream = new TTreeSRedirector("performanceSummary.root","recreate");
  keepArray=new TObjArray();
  AliDrawStyle::ApplyStyle("figTemplate2");
  GetNclReport(hisArray,keepArray);
  GetDCAReport(hisArray,keepArray);
  GetChi2Report();
  makeP2Report();
  makeP4Report();
}



void SetMetadata(){
  //
  chain->SetAlias("logTracks5","log(1+ntracks/5.)");
  chain->SetAlias("esdTrackPt","esdTrack.Pt()");
  chain->SetAlias("esdTrackQPt","esdTrack.fP[4]");
  chain->SetAlias("esdTrackfIpPt","esdTrack.fIp.Pt()");
  chain->SetAlias("phiInner","atan2(esdTrack.fIp.Py(),esdTrack.fIp.Px()+0)");
  chain->SetAlias("secInner","9*(atan2(esdTrack.fIp.Py(),esdTrack.fIp.Px()+0)/pi)+18*(esdTrack.fIp.Py()<0)");
  chain->SetAlias("dalphaQ","sign(esdTrack.fP[4])*(esdTrack.fIp.fP[0]/esdTrack.fIp.fX)");
  //
  chain->SetAlias("nclTPC","esdTrack.fTPCncls");
  chain->SetAlias("nclTRD","esdTrack.fTRDncls");
  chain->SetAlias("nclITS","esdTrack.fITSncls");
  chain->SetAlias("mdEdx","40./max(esdTrack.fTPCsignal,40.)");
  chain->SetAlias("smdEdx","sqrt(40./max(esdTrack.fTPCsignal,40.))");
  chain->SetAlias("nclROCA","esdTrack.GetTPCClusterInfo(3,1,0,159)");
  chain->SetAlias("nclROC0","esdTrack.GetTPCClusterInfo(3,1,0,62)");
  chain->SetAlias("nclROC1","esdTrack.GetTPCClusterInfo(3,1,63,126)");
  chain->SetAlias("nclROC2","esdTrack.GetTPCClusterInfo(3,1,127,159)");
  chain->SetAlias("nclFROCA","esdTrack.GetTPCClusterInfo(3,0,0,159)");
  chain->SetAlias("nclFROC0","esdTrack.GetTPCClusterInfo(3,0,0,62)");
  chain->SetAlias("nclFROC1","esdTrack.GetTPCClusterInfo(3,0,63,126)");
  chain->SetAlias("nclFROC2","esdTrack.GetTPCClusterInfo(3,0,127,159)");

  //
  chain->SetAlias("TPCASide","esdTrack.fIp.fP[1]>0");
  chain->SetAlias("TPCCSide","esdTrack.fIp.fP[1]<0");
  chain->SetAlias("TPCCross","esdTrack.fIp.fP[1]*esdTrack.fIp.fP[3]<0");
  //
  // V0 aliases
  //  
  chainV0->SetAlias("track0mP","1/track0.fIp.P()");
  chainV0->SetAlias("track1mP","1/track1.fIp.P()");
  chainV0->SetAlias("track0tgl","abs(track0.fP[3]+0)");
  chainV0->SetAlias("track1tgl","abs(track1.fP[3]+0)");
  chainV0->SetAlias("track0chi2","(track0.fTPCchi2/track0.fTPCncls+0)");
  chainV0->SetAlias("track1chi2","(track1.fTPCchi2/track1.fTPCncls+0)");
  chainV0->SetAlias("track0NclF","track0.GetTPCClusterInfo(3,0)");
  chainV0->SetAlias("track1NclF","track0.GetTPCClusterInfo(3,0)");
  chainV0->SetAlias("dedx0PionNorm","track0.fTPCsignal/dEdx0DPion");
  chainV0->SetAlias("dedx1PionNorm","track1.fTPCsignal/dEdx1DPion");
  chainV0->SetAlias("dedx0ProtonNorm","track0.fTPCsignal/dEdx0DProton");
  chainV0->SetAlias("dedx1ProtonNorm","track1.fTPCsignal/dEdx1DProton");

  chainV0->SetAlias("cleanPion0FromK0","K0Like0>0.05&&K0Like0>(LLike0+ALLike0+ELike0)*3&&abs(K0Delta)<0.006&&V0Like>0.1&&abs(track1.fTPCsignal/dEdx1DPion-50)<8&&v0.PtArmV0()>0.06");
  chainV0->SetAlias("cleanPion1FromK0","K0Like0>0.05&&K0Like0>(LLike0+ALLike0+ELike0)*3&&abs(K0Delta)<0.006&&V0Like>0.1&&abs(track0.fTPCsignal/dEdx0DPion-50)<8&&v0.PtArmV0()>0.06");
  chainV0->SetAlias("cleanPion0FromLambda","ALLike>0.05&&ALLike0>(K0Like0+LLike0+ELike0)*3&&abs(ALDelta)<0.006&&V0Like>0.1&&abs(track1.fTPCsignal/dEdx1DProton-50)<8");
  chainV0->SetAlias("cleanPion1FromLambda","LLike>0.05&&LLike0>(K0Like0+ALLike0+ELike0)*3&&abs(LDelta)<0.006&&V0Like>0.1&&abs(track0.fTPCsignal/dEdx0DProton-50)<8");
  //
  chainV0->SetAlias("cleanProton0FromLambda","LLike>0.05&&LLike0>(K0Like0+ALLike0+ELike0)*3&&abs(LDelta)<0.001&&V0Like>0.1&&abs(track1.fTPCsignal/dEdx1DPion-50)<8");
  chainV0->SetAlias("cleanProton1FromLambda","ALLike>0.05&&ALLike0>(K0Like0+LLike0+ELike0)*3&&abs(ALDelta)<0.001&&V0Like>0.1&&abs(track0.fTPCsignal/dEdx0DPion-50)<8");

}


void InitAnalysis(){
  //
  // get parameters
  // Init analysis
  ::Info("InitAnalysis()","START");
  pcstream = new TTreeSRedirector("performanceHisto.root","recreate");
  pcstream->GetFile()->cd();

  run=TString(gSystem->Getenv("run")).Atoi();
  period=gSystem->Getenv("period");
  if (gSystem->Getenv("deltaT")!=NULL) deltaT=TString(gSystem->Getenv("deltaT")).Atof();
  //
  // get chain
  chain = AliXRDPROOFtoolkit::MakeChainRandom("filtered.list","highPt",0,40000,0,1);
  chainV0 = AliXRDPROOFtoolkit::MakeChainRandom("filtered.list","V0s",0,40000,0,1);
  AliAnalysisTaskFilteredTree::SetDefaultAliasesHighPt(chain);
  AliAnalysisTaskFilteredTree::SetDefaultAliasesV0(chainV0);
  SetMetadata();
  //
  Int_t selected = chain->Draw("ntracks:mult:evtTimeStamp","","goff",100000);
  if (selected<=0){
    ::Error("performanceFiltered.InitAnalysis","Empty or corrupted input list");
    return;
  }
  ntracksEnd=TMath::KOrdStat(selected,chain->GetV1(),Int_t(selected*0.98))*1.02; // max Ntracks
  multEnd=TMath::KOrdStat(selected,chain->GetV2(),Int_t(selected*0.98))*1.02; // max Mult (primary)
  timeStart=TMath::KOrdStat(selected,chain->GetV3(),Int_t(selected*0.005)); 
  timeEnd=TMath::KOrdStat(selected,chain->GetV3(),Int_t(selected*0.995)); 
  timeBins=(timeEnd-timeStart)/deltaT+1;
  (*pcstream)<<"perf"<<
    "run="<<run<<
    "timeStart="<<timeStart<<
    "timeEnd="<<timeEnd;


  AliExternalInfo info;
  TTree* treeLogbook = info.GetTree("Logbook",period.Data(),"");

  if (treeLogbook==NULL) {
    //::Error("performanceFiltered.InitAnalysis","logbook tree not avaiable for period %s",period.Data()); 
    return;
  }
  Int_t entries = treeLogbook->Draw("DAQ_time_start:DAQ_time_end",TString::Format("run==%d",run).Data(),"goff");
  if (entries>0) {timeStart=treeLogbook->GetV1()[0];  timeEnd=treeLogbook->GetV2()[0];}
  timeBins=(timeEnd-timeStart)/deltaT+1;
  AliLumiTools lumiTool;
  lumiGraph = lumiTool.GetLumiFromCTP(run,"local:///cvmfs/alice.cern.ch/calibration/data/2015/OCDB/");
  //
  TVectorF vecX(timeBins), vecLumi(timeBins);
  for (Int_t itime=0; itime<timeBins; itime++){
    Double_t ctime= timeStart+(timeEnd-timeStart)*(itime+0.5)/timeBins;
    vecX[itime]=ctime;
    vecLumi[itime]=lumiGraph->Eval(ctime);    
  }
  TGraph* grLumiBin=new TGraph(timeBins, vecX.GetMatrixArray(), vecLumi.GetMatrixArray());
  (*pcstream)<<"perf"<<
    "lumiGraph.="<<lumiGraph<<
    "lumiBin.="<<grLumiBin;
}


TObjArray * FillPerfomanceHisto(Int_t maxEntries){
  //
  // Fill perfomance histograms
  //      return array of histograms
  //
  /*
    Int_t maxEntries=200000; 
  */
  Int_t nTracks=chain->GetEntries();
  chain->SetEstimate(chain->GetEntries());
  TString timeRange=TString::Format( "%d,%.0f,%.0f",timeBins,timeStart, timeEnd);
  //
  TString defaultCut="esdTrack.GetTPCClusterInfo(3,1)>60&&esdTrack.IsOn(0x1)>0";
  TString defaultCutMatch="esdTrack.GetTPCClusterInfo(3,1)>60";
  const Int_t nqaHistos=23;
  const char * qaHistos[nqaHistos]={"nclITS","nclTPC","nclTRD",		\
				    "normChi2ITS","normChi2TPC","normChi2TRD", \
				    "nclROC0","nclROC1","nclROC2", "nclROCA", \
				    "nclFROC0","nclFROC1","nclFROC2","nclFROCA", \
				    "deltaPC2Norm", "deltaPC3Norm", "deltaPC4Norm", \
				    "pullPC2", "pullPC3", "pullPC4", \
				    "covarPC2Norm", "covarPC3Norm", "covarPC4Norm"};

  const Int_t histosBins[nqaHistos]={8,80,80,	\
				     50,50,50,	\
				     64,64,32,160,	\
				     55,55,55,55,\
				     60,60,50, \
				     50,50,50, \
				     50,50,50 };
  const Double_t histosMin[nqaHistos]={0,0,0,	\
				       0,0,0,	\
				       0,0,0,0,	\
				       0,0,0,0, \
				       -0.015,-0.015,-0.1, \
				       -10,-10,-10, \
				       0.00,0.00,0.0};
  const Double_t histosMax[nqaHistos]={8,160,160,	\
				       10,10,10,	\
				       64,64,32,160,	\
				       1.1,1.1,1.1,1.1, \
				       0.015,0.015,0.1, \
				       10,10,10, \
				       0.01,0.01,0.1};
  

  TString hisString="";
  {
    // Standard kinematic histograms
    hisString+="esdTrackPt:#esdTrack.fTPCncls>60>>hisPtAll(100,0,30);"; 
    hisString+="esdTrackPt:#(esdTrack.fFlags&0x4)>0>>hisPtITS(100,1,30);";    
    hisString+="esdTrackfIpPt:#(esdTrack.fFlags&0x4)>0>>hisPtTPCOnly(100,1,30);";  
    hisString+="esdTrackQPt:tgl:secInner:#esdTrack.fTPCncls>60>>hisQptTglSecAll(40,-2,2,10,-1,1,90,0,18);"; 
  }
  // QA variables histograms
  for (Int_t iPar=0; iPar<nqaHistos; iPar++){
    // 
    hisString+=TString::Format("%s:qPt:tgl:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn01&&nclCut>>qahis%s_v_qPt_tgl(%d,%f,%f,200,-5,5,10,-1,1);", \
			       qaHistos[iPar],qaHistos[iPar],histosBins[iPar],histosMin[iPar],histosMax[iPar]);
    hisString+=TString::Format("%s:qPt:tgl:logTracks5:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn01&&nclCut>>qahis%s_v_qPt_tgl_LogTracks5(%d,%f,%f,50,-5,5,10,-1,1,10,0,10.);", \
			       qaHistos[iPar],qaHistos[iPar],histosBins[iPar],histosMin[iPar],histosMax[iPar]);
    hisString+=TString::Format("%s:qPt:tgl:smdEdx:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn01&&nclCut>>qahis%s_v_qPt_tgl_smdEdx(%d,%f,%f,50,-5,5,10,-1,1,10,0,1);", \
			       qaHistos[iPar],qaHistos[iPar],histosBins[iPar],histosMin[iPar],histosMax[iPar]);
    hisString+=TString::Format("%s:qPt:tgl:smdEdx:logTracks5:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn01&&nclCut>>qahis%s_v_qPt_tgl_smdEdx_LogTracks5(%d,%f,%f,50,-5,5,10,-1,1,10,0,1,10,0,10.);", \
			       qaHistos[iPar],qaHistos[iPar],histosBins[iPar],histosMin[iPar],histosMax[iPar]);
    hisString+=TString::Format("%s:qPt:tgl:dalphaQ:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn01&&nclCut>>qahis%s_v_qPt_tgl_dalphaQ(%d,%f,%f,48,-3,3,10,-1,1,50,-0.18,0.18);", \
			       qaHistos[iPar],qaHistos[iPar],histosBins[iPar],histosMin[iPar],histosMax[iPar]);
    hisString+=TString::Format("%s:qPt:tgl:alphaV:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn01&&nclCut>>qahis%s_v_qPt_tgl_alphaV(%d,%f,%f,48,-3,3,10,-1,1,90,-3.145,3.145);",\
			       qaHistos[iPar],qaHistos[iPar],histosBins[iPar],histosMin[iPar],histosMax[iPar]);
  }
  // Matchin efficiency histograms for primary +-4 sigma tracks
  //
  const Int_t nmatchHistos=4;
  const char * matchHistos[nmatchHistos]={"ITSOn","ITSRefit","TPCRefit","TRDOn"};
   // QA variables histograms 
  TString hisMatch="";
  for (Int_t iPar=0; iPar<nmatchHistos; iPar++){
    // 
    hisMatch+=TString::Format("%s:qPt:tgl:#TPCOn&&TOFOn&&IsPrim4&&IsPrim4TPC>>matchhis%s_v_qPt_tgl(2,-0.5,1.5,200,-5,5,10,-1,1);", \
			       matchHistos[iPar],matchHistos[iPar]);
    hisMatch+=TString::Format("%s:qPt:tgl:dalphaQ:#TPCOn&&TOFOn&&IsPrim4&&IsPrim4TPC>>matchhis%s_v_qPt_tgl_dalphaQ(2,-0.5,1.5,48,-3,3,10,-1,1,50,-0.18,0.18);", \
			       matchHistos[iPar],matchHistos[iPar]);
    hisMatch+=TString::Format("%s:qPt:tgl:alphaV:#TPCOn&&TOFOn&&IsPrim4&&IsPrim4TPC>>matchhis%s_v_qPt_tgl_alphaV(2,-0.5,1.5,48,-3,3,10,-1,1,90,-3.145,3.145);", \
			       matchHistos[iPar],matchHistos[iPar]);
    }
  



  {
    // N clusters per time;
    hisString+=TString::Format("esdTrack.fTPCncls:secInner:evtTimeStamp:#TPCASide>>hisTPCNclSecTimeA(100,60,160,180,0,18,%s);",timeRange.Data());
    hisString+=TString::Format("esdTrack.fTPCncls:secInner:evtTimeStamp:#TPCCSide>>hisTPCNclSecTimeC(100,60,160,180,0,18,%s);",timeRange.Data());
    hisString+=TString::Format("esdTrack.fTPCncls:secInner:evtTimeStamp:#TPCCross>>hisTPCNclSecTimeCross(100,60,160,180,0,18,%s);",timeRange.Data());
    // N assigned clusters fraction
    hisString+=TString::Format("esdTrack.GetTPCClusterInfo(3,0,0,63):secInner:evtTimeStamp:#TPCASide>>hisTPCClFracSecTimeA(20,0.1,1,180,0,18,%s);",timeRange.Data());
    hisString+=TString::Format("esdTrack.GetTPCClusterInfo(3,0,0,63):secInner:evtTimeStamp:#TPCCSide>>hisTPCClFracSecTimeC(20,0.1,1,180,0,18,%s);",timeRange.Data());
    hisString+=TString::Format("esdTrack.GetTPCClusterInfo(3,0,0,63):secInner:evtTimeStamp:#TPCCross>>hisTPCClFracSecTimeCross(20,0.1,1,180,0,18,%s);",timeRange.Data());
    // Chi2 histograms ITS
    hisString+=TString::Format("normChi2ITS:secInner:evtTimeStamp:#TPCASide>>hisChi2ITSSecTimeA(50,0,10,180,0,18,%s);",timeRange.Data());
    hisString+=TString::Format("normChi2ITS:secInner:evtTimeStamp:#TPCCSide>>hisChi2ITSSecTimeC(50,0,10,180,0,18,%s);",timeRange.Data());
    hisString+=TString::Format("normChi2ITS:secInner:evtTimeStamp:#TPCCross>>hisChi2ITSSecTimeCross(50,0,10,180,0,18,%s);",timeRange.Data());
    hisString+=TString::Format("normChi2ITS:tgl:ntracks>>hisChi2ITSP3NTracks(50,0,10,10,-1,1,10,0,%d);",ntracksEnd);
    hisString+=TString::Format("normChi2ITS:tgl:mult>>hisChi2ITSP3Mult(50,0,10,10,-1,1,10,0,%d);",multEnd);
    // Chi2 histograms TPC    
    hisString+=TString::Format("normChi2TPC:secInner:evtTimeStamp:#TPCASide>>hisChi2TPCSecTimeA(50,0,10,180,0,18,%s);",timeRange.Data());
    hisString+=TString::Format("normChi2TPC:secInner:evtTimeStamp:#TPCCSide>>hisChi2TPCSecTimeC(50,0,10,180,0,18,%s);",timeRange.Data());
    hisString+=TString::Format("normChi2TPC:secInner:evtTimeStamp:#TPCCross>>hisChi2TPCSecTimeCross(50,0,10,180,0,18,%s);",timeRange.Data());
    hisString+=TString::Format("normChi2TPC:tgl:ntracks>>hisChi2TPCP3NTracks(50,0,10,10,-1,1,10,0,%d);",ntracksEnd);
    hisString+=TString::Format("normChi2TPC:tgl:mult>>hisChi2TPCP3Mult(50,0,10,10,-1,1,10,0,%d);",multEnd);
    // DCAr histograms normalized to Pt
    hisString+=TString::Format("normDCAR:secInner:evtTimeStamp:#TPCASide&&abs(esdTrack.fD)<0.2>>hisTPCDCARSecTimeA(100,-3,3,180,0,18,%s);",timeRange.Data());
    hisString+=TString::Format("normDCAR:secInner:evtTimeStamp:#TPCCSide&&abs(esdTrack.fD)<0.2>>hisTPCDCARSecTimeC(100,-3,3,180,0,18,%s);",timeRange.Data());
    hisString+=TString::Format("normDCAR:secInner:evtTimeStamp:#TPCCross&&abs(esdTrack.fD)<0.2>>hisTPCDCARSecTimeCross(100,-3,3,180,0,18,%s);",timeRange.Data());
    // DCAz histograms normalized to Pt
    hisString+=TString::Format("normDCAZ:secInner:evtTimeStamp:#TPCASide&&abs(esdTrack.fD)<0.2>>hisTPCDCAZSecTimeA(100,-3,3,180,0,18,%s);",timeRange.Data());
    hisString+=TString::Format("normDCAZ:secInner:evtTimeStamp:#TPCCSide&&abs(esdTrack.fD)<0.2>>hisTPCDCAZSecTimeC(100,-3,3,180,0,18,%s);",timeRange.Data());
    hisString+=TString::Format("normDCAZ:secInner:evtTimeStamp:#TPCCross&&abs(esdTrack.fD)<0.2>>hisTPCDCAZSecTimeCross(100,-3,3,180,0,18,%s);",timeRange.Data());
  }
//
  // Kinematics matching
  Double_t range[5]={3,3,0.01,0.01,0.05};
  Double_t rangeP[5]={8,8,8,8,8};
  Double_t rangeCITS[5]={0.2,0.2,0.01,0.01,0.05};
  Double_t fnull=0;
  for (Int_t iPar=0; iPar<5; iPar++){
    // 
    hisString+=TString::Format("deltaP%d:qPt:tgl:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn01&&nclCut>>hisDeltaP%d_Allv_qPt_tgl(400,%f,%f,200,-5,5,10,-1,1);",iPar,iPar,-range[iPar],range[iPar]);
    hisString+=TString::Format("deltaP%d:qPt:tgl:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn01&&nclCut&&TRDOn>>hisDeltaP%d_TRDv_qPt_tgl(400,%f,%f,200,-5,5,10,-1,1);",iPar,iPar,-range[iPar],range[iPar]);
    hisString+=TString::Format("pullP%d:qPt:tgl:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn01&&nclCut>>hisPullP%d_Allv_qPt_tgl(400,-8,8,200,-5,5,10,-1,1);",iPar,iPar);
    hisString+=TString::Format("pullP%d:qPt:tgl:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn01&&nclCut&&TRDOn>>hisPullP%d_TRDv_qPt_tgl(400,-8,8,200,-5,5,10,-1,1);",iPar,iPar);
    hisString+=TString::Format("covarP%d:qPt:tgl:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn01&&nclCut>>hisCovarP%d_Allv_qPt_tgl(400,%f,%f,200,-5,5,10,-1,1);",iPar,iPar,fnull,range[iPar]);
    hisString+=TString::Format("covarP%d:qPt:tgl:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn01&&nclCut&&TRDOn>>hisCovarP%d_TRDv_qPt_tgl(400,%f,%f,200,-5,5,10,-1,1);",iPar,iPar,fnull,range[iPar]);
    hisString+=TString::Format("covarP%dITS:qPt:tgl:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn01&&nclCut>>hisCovarP%dITS_Allv_qPt_tgl(400,%f,%f,200,-5,5,10,-1,1);",iPar,iPar,fnull,rangeCITS[iPar]);
    hisString+=TString::Format("covarP%dITS:qPt:tgl:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn01&&nclCut&&TRDOn>>hisCovarP%dITS_TRDv_qPt_tgl(400,%f,%f,200,-5,5,10,-1,1);",iPar,iPar,fnull,rangeCITS[iPar]);
    // Residual calibration histogramming
    hisString+=TString::Format("deltaP%d:qPt:tgl:alphaV:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn01&&nclCut>>hisDeltaP%d_Allv_qPt_tgl_alphaV(100,%f,%f,48,-3,3,10,-1,1,90,-3.145,3.145);",iPar,iPar, -range[iPar],range[iPar]);
    hisString+=TString::Format("deltaP%d:qPt:tgl:alphaV:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn01&&nclCut&&TRDOn>>hisDeltaP%d_TRDv_qPt_tgl_alphaV(100,%f,%f,48,-3,3,10,-1,1,90,-3.145,3.145);",iPar,iPar, -range[iPar],range[iPar]);
    hisString+=TString::Format("pullP%d:qPt:tgl:alphaV:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn01&&nclCut>>hisPullP%d_Allv_qPt_tgl_alphaV(100,%f,%f,48,-3,3,10,-1,1,90,-3.145,3.145);",iPar,iPar, -rangeP[iPar],rangeP[iPar]);
    hisString+=TString::Format("pullP%d:qPt:tgl:alphaV:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn01&&nclCut&&TRDOn>>hisPullP%d_TRDv_qPt_tgl_alphaV(100,%f,%f,48,-3,3,10,-1,1,90,-3.145,3.145);",iPar,iPar, -rangeP[iPar],rangeP[iPar]);
    hisString+=TString::Format("covarP%d:qPt:tgl:alphaV:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn01&&nclCut>>hisCovarP%d_Allv_qPt_tgl_alphaV(100,%f,%f,48,-3,3,10,-1,1,90,-3.145,3.145);",iPar,iPar, fnull,range[iPar]);
    hisString+=TString::Format("covarP%d:qPt:tgl:alphaV:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn01&&nclCut&&TRDOn>>hisCovarP%d_TRDv_qPt_tgl_alphaV(100,%f,%f,48,-3,3,10,-1,1,90,-3.145,3.145);",iPar,iPar, fnull,range[iPar]);
    hisString+=TString::Format("covarP%dITS:qPt:tgl:alphaV:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn01&&nclCut>>hisCovarP%dITS_Allv_qPt_tgl_alphaV(100,%f,%f,48,-3,3,10,-1,1,90,-3.145,3.145);",iPar,iPar, fnull,rangeCITS[iPar]);
    hisString+=TString::Format("covarP%dITS:qPt:tgl:alphaV:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn01&&nclCut&&TRDOn>>hisCovarP%dITS_TRDv_qPt_tgl_alphaV(100,%f,%f,48,-3,3,10,-1,1,90,-3.145,3.145);",iPar,iPar, fnull,rangeCITS[iPar]);
    // Edge Effect histogramming
    hisString+=TString::Format("deltaP%d:qPt:tgl:dalphaQ:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn01&&nclCut>>hisDeltaP%d_Allv_qPt_tgl_dalphaQ(100,%f,%f,48,-3,3,10,-1,1,50,-0.18,0.18);",iPar,iPar, -range[iPar],range[iPar]);
    hisString+=TString::Format("deltaP%d:qPt:tgl:dalphaQ:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn01&&nclCut&&TRDOn>>hisDeltaP%d_TRDv_qPt_tgl_dalphaQ(100,%f,%f,48,-3,3,10,-1,1,50,-0.18,0.18);",iPar,iPar, -range[iPar],range[iPar]);
    hisString+=TString::Format("pullP%d:qPt:tgl:dalphaQ:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn01&&nclCut>>hisPullP%d_Allv_qPt_tgl_dalphaQ(100,%f,%f,48,-3,3,10,-1,1,50,-0.18,0.18);",iPar,iPar, -rangeP[iPar],rangeP[iPar]);
    hisString+=TString::Format("pullP%d:qPt:tgl:dalphaQ:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn01&&nclCut&&TRDOn>>hisPullP%d_TRDv_qPt_tgl_dalphaQ(100,%f,%f,48,-3,3,10,-1,1,50,-0.18,0.18);",iPar,iPar, -rangeP[iPar],rangeP[iPar]);
    hisString+=TString::Format("covarP%d:qPt:tgl:dalphaQ:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn01&&nclCut>>hisCovarP%d_Allv_qPt_tgl_dalphaQ(100,%f,%f,48,-3,3,10,-1,1,50,-0.18,0.18);",iPar,iPar, fnull,range[iPar]);
    hisString+=TString::Format("covarP%d:qPt:tgl:dalphaQ:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn01&&nclCut&&TRDOn>>hisCovarP%d_TRDv_qPt_tgl_dalphaQ(100,%f,%f,48,-3,3,10,-1,1,50,-0.18,0.18);",iPar,iPar, fnull,range[iPar]);
    hisString+=TString::Format("covarP%dITS:qPt:tgl:dalphaQ:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn01&&nclCut>>hisCovarP%dITS_Allv_qPt_tgl_dalphaQ(100,%f,%f,48,-3,3,10,-1,1,50,-0.18,0.18);",iPar,iPar, fnull,rangeCITS[iPar]);
    hisString+=TString::Format("covarP%dITS:qPt:tgl:dalphaQ:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn01&&nclCut&&TRDOn>>hisCovarP%dITS_TRDv_qPt_tgl_dalphaQ(100,%f,%f,48,-3,3,10,-1,1,50,-0.18,0.18);",iPar,iPar, fnull,rangeCITS[iPar]);
    //
    // multiplicity 
    //
    hisString+=TString::Format("deltaP%d:qPt:tgl:logTracks5:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn01&&nclCut>>hisDeltaP%d_Allv_qPt_tgl_logTracks5(400,%f,%f,200,-5,5,10,-1,1,10,0,10);",iPar,iPar,-range[iPar],range[iPar]);
    hisString+=TString::Format("deltaP%d:qPt:tgl:logTracks5:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn01&&nclCut&&TRDOn>>hisDeltaP%d_TRDv_qPt_tgl_logTracks5(400,%f,%f,200,-5,5,10,-1,1,10,0,10);",iPar,iPar,-range[iPar],range[iPar]);
    hisString+=TString::Format("pullP%d:qPt:tgl:logTracks5:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn01&&nclCut>>hisPullP%d_Allv_qPt_tgl_logTracks5(400,-8,8,200,-5,5,10,-1,1,10,0,10);",iPar,iPar);
    hisString+=TString::Format("pullP%d:qPt:tgl:logTracks5:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn01&&nclCut&&TRDOn>>hisPullP%d_TRDv_qPt_tgl_logTracks5(400,-8,8,200,-5,5,10,-1,1,10,0,10);",iPar,iPar);

  }

  TString hisV0String="";
  chainV0->SetAlias("K0cut","abs(K0PIDPull)<2&&(abs(LPull)>2.&&abs(ALPull)>2.)");
  {
    // K0 performance
    hisV0String+="K0Delta:mpt:tglV0:#K0cut>>hisK0DMassQPtTgl(100,-0.03,0.03,80,0,2,10,-1,1);";  
    hisV0String+="K0Pull:mpt:tglV0:#K0cut>>hisK0PullQPtTgl(100,-6.0,6.0,80,0,2,10,-1,1);";
    // K0 resolution/maps - in respec to sector edge
    hisV0String+="K0Delta:mpt:tglV0:dalphaV0:#K0cut>>hisK0DMassQPtTglDSec(100,-0.03,0.03,10,0,1,10,-1,1,10,0.0,0.35);";  
    hisV0String+="K0Pull:mpt:tglV0:dalphaV0:#K0cut>>hisK0PullQPtTglDSec(100,-6.0,6.0,10,0,1,10,-1,1,10,0.0,0.35);";
    // K0 resolution/maps - 
    hisV0String+="K0Delta:mpt:tglV0:alphaV0:#K0cut>>hisK0DMassQPtTglAlpha(100,-0.03,0.03,10,0,1,5,-1,1,18,-3.1415,3.1415);";  
    hisV0String+="K0Pull:mpt:tglV0:alphaV0:#K0cut>>hisK0PullQPtTglAlpha(100,-6.0,6.0,10,0,1,5,-1,1,18,-3.1415,3.1415);";
  }  
  {
    // dEdx clean pions
    hisV0String+="dedx0PionNorm:track0mP:track0tgl:#cleanPion0FromK0>>v0hisdEdxPion0FromK0_v_qPt_tgl(100,20,100,60,0,3,10,0,1);";
    hisV0String+="dedx1PionNorm:track1mP:track1tgl:#cleanPion1FromK0>>v0hisdEdxPion1FromK0_v_qPt_tgl(100,20,100,60,0,3,10,0,1);";
    hisV0String+="dedx0PionNorm:track0mP:track0tgl:#cleanPion0FromLambda>>v0hisdEdxPion0FromLambda_v_qPt_tgl(100,20,100,60,0,3,10,0,1);";
    hisV0String+="dedx1PionNorm:track1mP:track1tgl:#cleanPion1FromLambda>>v0hisdEdxPion1FromLambda_v_qPt_tgl(100,20,100,60,0,3,10,0,1);";
    // dEdxclean protons
    hisV0String+="dedx0ProtonNorm:track0mP:track0tgl:#cleanProton0FromLambda>>v0hisdEdxProton0FromLambda_v_qPt_tgl(100,20,100,60,0,3,10,0,1);";
    hisV0String+="dedx1ProtonNorm:track1mP:track1tgl:#cleanProton1FromLambda>>v0hisdEdxProton1FromLambda_v_qPt_tgl(100,20,100,60,0,3,10,0,1);";
    // chi2 clean sample
    hisV0String+="track0chi2:track0mP:track0tgl:#(cleanPion0FromK0||cleanPion0FromLambda)>>v0hisChi2Pion0_v_qPt_tgl(100,0,10,60,0,3,10,0,1);";
    hisV0String+="track1chi2:track1mP:track1tgl:#(cleanPion1FromK0||cleanPion1FromLambda)>>v0hisChi2Pion1_v_qPt_tgl(100,0,10,60,0,3,10,0,1);";
    hisV0String+="track0chi2:track0mP:track0tgl:#cleanProton0FromLambda>>v0hisChi2Proton0_v_qPt_tgl(100,0,10,60,0,3,10,0,1);";
    hisV0String+="track1chi2:track1mP:track1tgl:#cleanProton1FromLambda>>v0hisChi2Proton1_v_qPt_tgl(100,0,10,60,0,3,10,0,1);";
    // nclf
    hisV0String+="track0NclF:track0mP:track0tgl:#(cleanPion0FromK0||cleanPion0FromLambda)>>v0hisNclFPion0_v_qPt_tgl(60,0.5,1.1,60,0,3,10,0,1);";
    hisV0String+="track1NclF:track1mP:track1tgl:#(cleanPion1FromK0||cleanPion1FromLambda)>>v0hisNclFPion1_v_qPt_tgl(60,0.5,1.1,60,0,3,10,0,1);";
    hisV0String+="track0NclF:track0mP:track0tgl:#cleanProton0FromLambda>>v0hisNclFProton0_v_qPt_tgl(60,0.5,1.1,60,0,3,10,0,1);";
    hisV0String+="track1NclF:track1mP:track1tgl:#cleanProton1FromLambda>>v0hisNclFProton1_v_qPt_tgl(60,0.5,1.1,60,0,3,10,0,1);";

  }

  //
  TStopwatch timer;
  //
  timer.Start();
  hisArrayV0 = AliTreePlayer::MakeHistograms(chainV0, hisV0String, "1",0,maxEntries,200000,15);
  timer.Print();
  (*pcstream).GetFile()->cd();
  for (Int_t iKey=0; iKey<hisArrayV0->GetEntries(); iKey++){
    hisArrayV0->At(iKey)->Write( hisArrayV0->At(iKey)->GetName());
  }
  delete  hisArrayV0;
  //
  timer.Start();
  hisArray = AliTreePlayer::MakeHistograms(chain, hisString, defaultCut,0,maxEntries,200000,15);
  timer.Print();
  (*pcstream).GetFile()->cd();
  for (Int_t iKey=0; iKey<hisArray->GetEntries(); iKey++){
    hisArray->At(iKey)->Write( hisArray->At(iKey)->GetName());
  }
  delete hisArray;
  timer.Start();
  hisArray = AliTreePlayer::MakeHistograms(chain, hisMatch, defaultCutMatch,0,maxEntries,200000,15);
  timer.Print();
  (*pcstream).GetFile()->cd();
  //  hisArray->Write("perfArray",  TObjArray::kSingleKey);
  for (Int_t iKey=0; iKey<hisArray->GetEntries(); iKey++){
    hisArray->At(iKey)->Write( hisArray->At(iKey)->GetName());
  }
  delete hisArray;
  //
  //  hisArray->Write("perfArray");
  //  hisArrayV0->Write("perfArrayV0");
  (*pcstream).GetFile()->Flush();
  return 0;
}




void GetNclReport(TObjArray * hisArray,  TObjArray *keepArray ){
  //
  // NCl report
  // 
  TString drawExpressionNcl="";
  drawExpressionNcl="[3,3,3]:";
  drawExpressionNcl+="%Ogridx,gridy;hisTPCNclSecTimeA(0,100,0,180,0,10000)(0)(err):";
  drawExpressionNcl+="%Ogridx,gridy;hisTPCNclSecTimeC(0,100,0,180,0,10000)(0)(err):";
  drawExpressionNcl+="%Ogridx,gridy;hisTPCNclSecTimeC(0,100,0,180,0,10000)(0)(err):";
  drawExpressionNcl+="hisTPCNclSecTimeA(0,100,0,179,0,10000)(0,1)(f-mean p):";
  drawExpressionNcl+="hisTPCNclSecTimeC(0,100,0,179,0,10000)(0,1)(f-mean p):";
  drawExpressionNcl+="hisTPCNclSecTimeCross(0,100,0,179,0,10000)(0,1)(f-mean p):";
  drawExpressionNcl+="%Otimex,gridx,gridy;";
  drawExpressionNcl+="hisTPCNclSecTimeA(0,100,0,179,0,10000)(0,2)(f-mean p);hisTPCNclSecTimeC(0,100,0,179,0,10000)(0,2)(f-mean p);:";
  drawExpressionNcl+="%Otimex,gridx,gridy;";
  drawExpressionNcl+="hisTPCNclSecTimeA(0,100,U0,U18,0,10000)(0,2)(f-mean p);hisTPCNclSecTimeA(0,100,U1,U9,0,10000)(0,2)(f-mean p);hisTPCNclSecTimeA(0,100,U10,U18,0,10000)(0,2)(f-mean p);:";
  drawExpressionNcl+="%Otimex,gridx,gridy;";
  drawExpressionNcl+="hisTPCNclSecTimeC(0,100,U0,U18,0,10000)(0,2)(f-mean p);hisTPCNclSecTimeC(0,100,U1,U3,0,10000)(0,2)(f-mean p);hisTPCNclSecTimeC(0,100,U11,U13,0,10000)(0,2)(f-mean p);:";  
  TPad * padNcl = AliTreePlayer::DrawHistograms(0,hisArray,drawExpressionNcl,keepArray, 1+2+4+8);
  ((TCanvas*)padNcl)->SetWindowSize(1800,1000);
  padNcl->SaveAs("nclReport.pdf");
  padNcl->SaveAs("nclReport.png");
  {(*pcstream)<<"perf"<<
      "grNclASideTime.="<<keepArray->FindObject("hisTPCNclSecTimeA(0,100,U0,U18,0,10000)(0,2)(f-mean p)")<<
      "grNclCSideTime.="<<keepArray->FindObject("hisTPCNclSecTimeC(0,100,U0,U18,0,10000)(0,2)(f-mean p)")<<
      "grNclASideSector.="<<keepArray->FindObject("hisTPCNclSecTimeA(0,100,0,179,0,10000)(0,1)(f-mean p)")<<
      "grNclCSidesector.="<<keepArray->FindObject("hisTPCNclSecTimeA(0,100,0,179,0,10000)(0,1)(f-mean p)");      
  }
  // ClFrac findable report  
  TString drawExpressionClFrac="";
  drawExpressionClFrac="[3,3,3]:";
  drawExpressionClFrac+="%Ogridx,gridy;hisTPCClFracSecTimeA(0,100,0,180,0,10000)(0)(err):";
  drawExpressionClFrac+="%Ogridx,gridy;hisTPCClFracSecTimeC(0,100,0,180,0,10000)(0)(err):";
  drawExpressionClFrac+="%Ogridx,gridy;hisTPCClFracSecTimeC(0,100,0,180,0,10000)(0)(err):";
  drawExpressionClFrac+="hisTPCClFracSecTimeA(0,100,0,179,0,10000)(0,1)(f-mean p):";
  drawExpressionClFrac+="hisTPCClFracSecTimeC(0,100,0,179,0,10000)(0,1)(f-mean p):";
  drawExpressionClFrac+="hisTPCClFracSecTimeCross(0,100,0,179,0,10000)(0,1)(f-mean p):";
  drawExpressionClFrac+="%Otimex,gridx,gridy;";
  drawExpressionClFrac+="hisTPCClFracSecTimeA(0,100,0,179,0,10000)(0,2)(f-mean p);hisTPCClFracSecTimeC(0,100,0,179,0,10000)(0,2)(f-mean p);:";
  drawExpressionClFrac+="%Otimex,gridx,gridy;";
  drawExpressionClFrac+="hisTPCClFracSecTimeA(0,100,U0,U18,0,10000)(0,2)(f-mean p);hisTPCClFracSecTimeA(0,100,U1,U9,0,10000)(0,2)(f-mean p);hisTPCClFracSecTimeA(0,100,U10,U18,0,10000)(0,2)(f-mean p);:";
  drawExpressionClFrac+="%Otimex,gridx,gridy;";
  drawExpressionClFrac+="hisTPCClFracSecTimeC(0,100,U0,U18,0,10000)(0,2)(f-mean p);hisTPCClFracSecTimeC(0,100,U1,U3,0,10000)(0,2)(f-mean p);hisTPCClFracSecTimeC(0,100,U11,U13,0,10000)(0,2)(f-mean p);:";
  
  TPad * padClFrac = AliTreePlayer::DrawHistograms(0,hisArray,drawExpressionClFrac,keepArray, 1+2+4+8);
  ((TCanvas*)padClFrac)->SetWindowSize(1800,1000);
  padClFrac->SaveAs("tpcclFractionReport.pdf");
  padClFrac->SaveAs("tpcclfractionReport.png");
  padClFrac->SaveAs("tpcclfractionReport.C");
  {(*pcstream)<<"perf"<<
      "grClFracASideTime.="<<keepArray->FindObject("hisTPCClFracSecTimeA(0,100,U0,U18,0,10000)(0,2)(f-mean p)")<<
      "grClFracCSideTime.="<<keepArray->FindObject("hisTPCClFracSecTimeC(0,100,U0,U18,0,10000)(0,2)(f-mean p)")<<
      "grClFracASideSector.="<<keepArray->FindObject("hisTPCClFracSecTimeA(0,100,0,179,0,10000)(0,1)(f-mean p)")<<
      "grClFracCSidesector.="<<keepArray->FindObject("hisTPCClFracSecTimeA(0,100,0,179,0,10000)(0,1)(f-mean p)");      
  }
}

void GetDCAReport(TObjArray * hisArray,  TObjArray *keepArray ){
  //
  // DCA run report
  //
  // 1.) TPC - DCAR
  TString drawExpressionTPCDCAR="";
  drawExpressionTPCDCAR="[1,2,2,1,1,1]:";
  drawExpressionTPCDCAR+="%Ogridx,gridy;";
  drawExpressionTPCDCAR+="hisTPCDCARSecTimeA(0,100,0,180,0,10000)(0)(err);hisTPCDCARSecTimeC(0,100,0,180,0,10000)(0)(err):";
  //
  drawExpressionTPCDCAR+="%Ogridx,gridy;hisTPCDCARSecTimeA(0,100,0,180,0,10000)(0,1)(f-rms p);:";
  drawExpressionTPCDCAR+="%Ogridx,gridy;hisTPCDCARSecTimeC(0,100,0,180,0,10000)(0,1)(f-rms p):";
  drawExpressionTPCDCAR+="%Ogridx,gridy;hisTPCDCARSecTimeA(0,100,0,180,0,10000)(0,1)(f-mean p);:";
  drawExpressionTPCDCAR+="%Ogridx,gridy;hisTPCDCARSecTimeC(0,100,0,180,0,10000)(0,1)(f-mean p):";
  drawExpressionTPCDCAR+="%Otimex,gridx,gridy;";
  drawExpressionTPCDCAR+="hisTPCDCARSecTimeA(0,100,U0,U18,0,10000)(0,2)(f-rms p);hisTPCDCARSecTimeA(0,100,U1,U9,0,10000)(0,2)(f-rms p);hisTPCDCARSecTimeA(0,100,U10,U18,0,10000)(0,2)(f-rms p);:";
  drawExpressionTPCDCAR+="%Otimex,gridx,gridy;";
  drawExpressionTPCDCAR+="hisTPCDCARSecTimeC(0,100,U0,U18,0,10000)(0,2)(f-rms p);hisTPCDCARSecTimeC(0,100,U1.5,U2.5,0,10000)(0,2)(f-rms p);hisTPCDCARSecTimeC(0,100,U11.5,U12.5,0,10000)(0,2)(f-rms p);:";
  drawExpressionTPCDCAR+="%Otimex,gridx,gridy;";
  drawExpressionTPCDCAR+="hisTPCDCARSecTimeA(0,100,U0,U18,0,10000)(0,2)(f-mean p);hisTPCDCARSecTimeC(0,100,U0,U18,0,10000)(0,2)(f-mean p);:";
  TPad * padTPCDCAR = AliTreePlayer::DrawHistograms(0,hisArray,drawExpressionTPCDCAR,keepArray,1+2+4+8);
  ((TCanvas*)padTPCDCAR)->SetWindowSize(1600,1000);
  padTPCDCAR->SaveAs("tpcDCARReport.png");
  padTPCDCAR->SaveAs("tpcDCARReport.pdf");
  padTPCDCAR->SaveAs("tpcDCARReport.C");
  {(*pcstream)<<"perf"<<    // export DCAR trending variables
      "grRMSTPCDCARSecA.="<<keepArray->FindObject("hisTPCDCARSecTimeA(0,100,0,180,0,10000)(0)(err)")<<
      "grRMSTPCDCARSecC.="<<keepArray->FindObject("hisTPCDCARSecTimeC(0,100,0,180,0,10000)(0)(err)")<<
      "grRMSTPCDCARTimeA.="<<keepArray->FindObject("hisTPCDCARSecTimeA(0,100,U0,U18,0,10000)(0,2)(f-rms p)")<<
      "grRMSTPCDCARTimeA.="<<keepArray->FindObject("hisTPCDCARSecTimeC(0,100,U0,U18,0,10000)(0,2)(f-rms p)")<<      
      "grMeanTPCDCARTimeA.="<<keepArray->FindObject("hisTPCDCARSecTimeA(0,100,U0,U18,0,10000)(0,2)(f-mean p)")<<
      "grMeanTPCDCARTimeA.="<<keepArray->FindObject("hisTPCDCARSecTimeC(0,100,U0,U18,0,10000)(0,2)(f-mean p)");      
  }
  // 2.) TPC - DCAZ
  TString drawExpressionTPCDCAZ="";
  drawExpressionTPCDCAZ="[1,2,1,1,1]:";
  drawExpressionTPCDCAZ+="%Ogridx,gridy;hisTPCDCAZSecTimeA(30,70,0,180,0,10000)(0)(err);hisTPCDCAZSecTimeC(0,100,0,180,0,10000)(0)(err):";
  drawExpressionTPCDCAZ+="%Ogridx,gridy;hisTPCDCAZSecTimeA(30,70,0,180,0,10000)(0,1)(f-rms p):";
  drawExpressionTPCDCAZ+="%Ogridx,gridy;hisTPCDCAZSecTimeC(30,70,0,180,0,10000)(0,1)(f-rms p):";
  drawExpressionTPCDCAZ+="%Otimex,gridx,gridy;";
  drawExpressionTPCDCAZ+="hisTPCDCAZSecTimeA(0,100,U0,U18,0,10000)(0,2)(f-rms p);hisTPCDCAZSecTimeA(0,100,U1,U9,0,10000)(0,2)(f-rms p);hisTPCDCAZSecTimeA(0,100,U10,U18,0,10000)(0,2)(f-rms p);:";
  drawExpressionTPCDCAZ+="%Otimex,gridx,gridy;";
  drawExpressionTPCDCAZ+="hisTPCDCAZSecTimeC(0,100,U0,U18,0,10000)(0,2)(f-rms p);hisTPCDCAZSecTimeC(0,100,U1.5,U2.5,0,10000)(0,2)(f-rms p);hisTPCDCAZSecTimeC(0,100,U11.5,U12.5,0,10000)(0,2)(f-rms p);:";
  drawExpressionTPCDCAZ+="%Otimex,gridx,gridy;";
  drawExpressionTPCDCAZ+="hisTPCDCAZSecTimeA(0,100,U0,U18,0,10000)(0,2)(f-mean p);hisTPCDCAZSecTimeC(0,100,U1.5,U2.5,0,10000)(0,2)(f-mean p);:";
  TPad * padTPCDCAZ = AliTreePlayer::DrawHistograms(0,hisArray,drawExpressionTPCDCAZ,keepArray,1+2+4+8);
  ((TCanvas*)padTPCDCAZ)->SetWindowSize(1600,1000);
  padTPCDCAR->SaveAs("tpcDCAZReport.png");
  padTPCDCAR->SaveAs("tpcDCAZReport.pdf");
  {(*pcstream)<<"perf"<<     // export DCAR trending variables
      "grRMSTPCDCAZSecA.="<<keepArray->FindObject("hisTPCDCAZSecTimeA(0,100,0,180,0,10000)(0)(err)")<<
      "grRMSTPCDCAZSecC.="<<keepArray->FindObject("hisTPCDCAZSecTimeC(0,100,0,180,0,10000)(0)(err)")<< 
      "grRMSTPCDCAZTimeA.="<<keepArray->FindObject("hisTPCDCAZSecTimeA(0,100,U0,U18,0,10000)(0,2)(f-rms p)")<<
      "grRMSTPCDCAZTimeA.="<<keepArray->FindObject("hisTPCDCAZSecTimeC(0,100,U0,U18,0,10000)(0,2)(f-rms p)")<<      
      "grMeanTPCDCAZTimeA.="<<keepArray->FindObject("hisTPCDCAZSecTimeA(0,100,U0,U18,0,10000)(0,2)(f-mean p)")<<
      "grMeanTPCDCAZTimeA.="<<keepArray->FindObject("hisTPCDCAZSecTimeC(0,100,U1.5,U2.5,0,10000)(0,2)(f-mean p)");     
 }

}


void GetChi2Report(){
  // Run propertes
  TString drawExpressionChi2TPC="";
  drawExpressionChi2TPC="[2,2,1]:";
  drawExpressionChi2TPC+="%Ogridx,gridy;hisChi2TPCSecTimeA(0,40,0,180,0,10)(0,1)(f-mean p);:";
  drawExpressionChi2TPC+="%Ogridx,gridy;hisChi2TPCSecTimeC(0,40,0,180,0,10)(0,1)(f-mean p);:";
  drawExpressionChi2TPC+="%Ogridx,gridy;hisChi2TPCP3NTracks(0,100,1,3,0,10)(0,2)(f-mean-p);hisChi2TPCP3NTracks(0,100,3,5,0,10)(0,2)(f-mean-p):";
  drawExpressionChi2TPC+="%Ogridx,gridy;hisChi2TPCP3NTracks(0,100,6,8,0,10)(0,2)(f-mean-p);hisChi2TPCP3NTracks(0,100,8,10,0,10)(0,2)(f-mean-p):";
  drawExpressionChi2TPC+="%Ogridx,gridy;hisChi2TPCP3NTracks(0,100,0,10,1,3)(0,1)(f-mean-p);hisChi2TPCP3NTracks(0,100,0,10,4,6)(0,1)(f-mean-p);hisChi2TPCP3NTracks(0,100,0,10,7,9)(0,1)(f-mean-p):";
  TPad * padChi2TPC = AliTreePlayer::DrawHistograms(0,hisArray,drawExpressionChi2TPC,keepArray, 1+2+4+8);
  ((TCanvas*)padChi2TPC)->SetWindowSize(1600,1000);
  padChi2TPC->SaveAs("tpcChi2Report.png");
  padChi2TPC->SaveAs("tpcChi2Report.pdf");
  padChi2TPC->SaveAs("tpcChi2Report.C");
  // Run propertes
  TString drawExpressionChi2ITS="";
  drawExpressionChi2ITS="[2,2,1]:";
  drawExpressionChi2ITS+="%Ogridx,gridy;hisChi2ITSSecTimeA(0,40,0,180,0,10)(0,1)(f-mean p);:";
  drawExpressionChi2ITS+="%Ogridx,gridy;hisChi2ITSSecTimeC(0,40,0,180,0,10)(0,1)(f-mean p);:";
  drawExpressionChi2ITS+="%Ogridx,gridy;hisChi2ITSP3NTracks(0,100,1,3,0,10)(0,2)(f-mean-p);hisChi2ITSP3NTracks(0,100,3,5,0,10)(0,2)(f-mean-p):";
  drawExpressionChi2ITS+="%Ogridx,gridy;hisChi2ITSP3NTracks(0,100,6,8,0,10)(0,2)(f-mean-p);hisChi2ITSP3NTracks(0,100,8,10,0,10)(0,2)(f-mean-p):";
  drawExpressionChi2ITS+="%Ogridx,gridy;hisChi2ITSP3NTracks(0,100,0,10,1,3)(0,1)(f-mean-p);hisChi2ITSP3NTracks(0,100,0,10,4,6)(0,1)(f-mean-p);hisChi2ITSP3NTracks(0,100,0,10,7,9)(0,1)(f-mean-p):";
  TPad * padChi2ITS = AliTreePlayer::DrawHistograms(0,hisArray,drawExpressionChi2ITS,keepArray, 1+2+4+8);
  ((TCanvas*)padChi2ITS)->SetWindowSize(1600,1000);
  padChi2ITS->SaveAs("itsChi2Report.png");
  padChi2ITS->SaveAs("itsChi2Report.pdf");
  padChi2ITS->SaveAs("itsChi2Report.C");
  // Time properties
  TString drawExpressionChi2TPCTime="";
  drawExpressionChi2TPCTime="[1,1,1]:";
  drawExpressionChi2TPCTime+="%Otimex,gridx,gridy;";
  drawExpressionChi2TPCTime+="hisChi2TPCSecTimeA(0,40,0,180,0,10000)(0,2)(f-mean p);hisChi2TPCSecTimeC(0,40,0,180,0,10000)(0,2)(f-mean p);:";
  drawExpressionChi2TPCTime+="%Otimex,gridx,gridy;";
  drawExpressionChi2TPCTime+="hisChi2TPCSecTimeA(0,40,U0,U18,0,10000)(0,2)(f-mean p);hisChi2TPCSecTimeA(0,40,U1,U9,0,10000)(0,2)(f-mean p);hisChi2TPCSecTimeA(0,40,U10,U18,0,10000)(0,2)(f-mean p);:";
  drawExpressionChi2TPCTime+="%Otimex,gridx,gridy;";
  drawExpressionChi2TPCTime+="hisChi2TPCSecTimeC(0,40,U0,U18,0,10000)(0,2)(f-mean p);hisChi2TPCSecTimeC(0,40,U1.5,U2.5,0,10000)(0,2)(f-mean p);hisChi2TPCSecTimeC(0,40,U11.5,U12.5,0,10000)(0,2)(f-mean p);:";
  TPad * padChi2TPCTime = AliTreePlayer::DrawHistograms(0,hisArray,drawExpressionChi2TPCTime,keepArray, 1+2+4+8);
  ((TCanvas*)padChi2TPCTime)->SetWindowSize(1600,1000);
  padChi2TPC->SaveAs("tpcChi2TimeReport.png");
  padChi2TPC->SaveAs("tpcChi2TimeReport.pdf");
  padChi2TPC->SaveAs("tpcChi2TimeReport.C");
  //
  Int_t hotSpots[8]={2,4,6,7,9, 20,30,35};
  TString drawExpressionChi2TPCSectorTime="[2,2,2,2]:";
  TString drawExpressionChi2ITSSectorTime="[2,2,2,2]:";
  for (Int_t iHot=0; iHot<8; iHot++){
    if (hotSpots[iHot]<18){
      drawExpressionChi2TPCSectorTime+=TString::Format("%Otimex,gridx,gridy;hisChi2TPCSecTimeA(0,40,U%.1f,U%.1f,0,10000)(0,2)(f-mean p):",\
						       hotSpots[iHot]-0.5,hotSpots[iHot]+0.5);
      drawExpressionChi2ITSSectorTime+=TString::Format("%Otimex,gridx,gridy;hisChi2ITSSecTimeA(0,40,U%.1f,U%.1f,0,10000)(0,2)(f-mean p):",\
						       hotSpots[iHot]-0.5,hotSpots[iHot]+0.5);
    }
    if (hotSpots[iHot]>18){
      drawExpressionChi2TPCSectorTime+=TString::Format("%Otimex,gridx,gridy;hisChi2TPCSecTimeC(0,40,U%.1f,U%.1f,0,10000)(0,2)(f-mean p):",\
						       hotSpots[iHot]-0.5-18,hotSpots[iHot]+0.5-18);
      drawExpressionChi2ITSSectorTime+=TString::Format("%Otimex,gridx,gridy;hisChi2ITSSecTimeC(0,40,U%.1f,U%.1f,0,10000)(0,2)(f-mean p):",\
						       hotSpots[iHot]-0.5-18,hotSpots[iHot]+0.5-18);
    }
  }
  TPad * padChi2TPCSectorTime = AliTreePlayer::DrawHistograms(0,hisArray,drawExpressionChi2TPCSectorTime,keepArray, 1+2+4+8);
  ((TCanvas*)padChi2TPCSectorTime)->SetWindowSize(1600,1000);
  padChi2TPCSectorTime->SaveAs("tpcChi2SectorTimeReport.png");
  padChi2TPCSectorTime->SaveAs("tpcChi2SectorTimeReport.pdf");
  padChi2TPCSectorTime->SaveAs("tpcChi2SectorTimeReport.C");
  //
  
  TPad * padChi2ITSSectorTime = AliTreePlayer::DrawHistograms(0,hisArray,drawExpressionChi2ITSSectorTime,keepArray, 1+2+4+8);
  ((TCanvas*)padChi2ITSSectorTime)->SetWindowSize(1600,1000);
  padChi2ITSSectorTime->SaveAs("itsChi2SectorTimeReport.png");
  padChi2ITSSectorTime->SaveAs("itsChi2SectorTimeReport.pdf");
  padChi2ITSSectorTime->SaveAs("itsChi2SectorTimeReport.C");



}

void makeP2Report(){
  TString drawExpressionP2="";
  drawExpressionP2="[1,1,2]:";
  //drawExpressionP2+="%Ogridx,gridy;hisDeltaP2CQPtTglTRD(0,400,80,120,0,10)(0,1)(f-mean p);:";
  drawExpressionP2+="%Ogridx,gridy;hisDeltaP2CQPtTglTRD(0,400,80,120,0,10)(0,1)(f-grms p);:";
  //drawExpressionP2+="%Ogridx,gridy;hisPullP2CQPtTglTRD(0,400,80,120,0,10)(0,1)(f-mean p);:";
  drawExpressionP2+="%Ogridx,gridy;hisPullP2CQPtTglTRD(0,400,80,120,0,10)(0,1)(f-grms p);:";
  drawExpressionP2+="%Ogridx,gridy;hisDeltaP2CQPtTglTRD(100,300,100,101,0,10)(0)(err);:";
  drawExpressionP2+="%Ogridx,gridy;hisPullP2CQPtTglTRD(10,90,100,101,0,10)(0)(err);:";
  TPad * padP2 = AliTreePlayer::DrawHistograms(0,hisArray,drawExpressionP2,keepArray, 1+2+4+8);
  padP2->SaveAs("deltaP2TPCtoFull.png");
  padP2->SaveAs("deltaP2TPCtoFull.C");

  TString drawExpressionDeltaP2Const="";
  drawExpressionDeltaP2Const="[1,1,2]:";
  //drawExpressionDeltaP2Const+="%Ogridx,gridy;hisDeltaP2ConstCQPtTglTRD(0,400,80,120,0,10)(0,1)(f-mean p);:";
  drawExpressionDeltaP2Const+="%Ogridx,gridy;hisDeltaP2ConstCQPtTglTRD(0,400,80,120,0,10)(0,1)(f-grms p);:";
  //drawExpressionDeltaP2Const+="%Ogridx,gridy;hisPullP2CQPtTglTRD(0,400,80,120,0,10)(0,1)(f-mean p);:";
  drawExpressionDeltaP2Const+="%Ogridx,gridy;hisPullP2ConstCQPtTglTRD(0,400,80,120,0,10)(0,1)(f-grms p);:";
  drawExpressionDeltaP2Const+="%Ogridx,gridy;hisDeltaP2ConstCQPtTglTRD(100,300,99,101,0,10)(0)(err);:";
  drawExpressionDeltaP2Const+="%Ogridx,gridy;hisPullP2ConstCQPtTglTRD(10,90,97,102,0,10)(0)(err);:";
  TPad * padDeltaP2const = AliTreePlayer::DrawHistograms(0,hisArray,drawExpressionDeltaP2Const,keepArray, 1+2+4+8);
  padDeltaP2const->SaveAs("deltaP2TPCConsttoFullConst.png");
  padDeltaP2const->SaveAs("deltaP2TPCConsttoFullConst.C");
}


void makeP4Report(){
  TString drawExpressionP4="";
  drawExpressionP4="[1,1,2]:";
  //drawExpressionP4+="%Ogridx,gridy;hisDeltaP4CQPtTglTRD(0,400,80,120,0,10)(0,1)(f-mean p);:";
  drawExpressionP4+="%Ogridx,gridy;hisDeltaP4CQPtTglTRD(0,400,80,120,0,10)(0,1)(f-grms p);:";
  //drawExpressionP4+="%Ogridx,gridy;hisPullP4CQPtTglTRD(0,400,80,120,0,10)(0,1)(f-mean p);:";
  drawExpressionP4+="%Ogridx,gridy;hisPullP4CQPtTglTRD(0,400,80,120,0,10)(0,1)(f-grms p);:";
  drawExpressionP4+="%Ogridx,gridy;hisDeltaP4CQPtTglTRD(100,300,99,102,0,10)(0)(err);:";
  drawExpressionP4+="%Ogridx,gridy;hisPullP4CQPtTglTRD(10,90,99,102,0,10)(0)(err);:";
  TPad * padP4 = AliTreePlayer::DrawHistograms(0,hisArray,drawExpressionP4,keepArray, 1+2+4+8);
  padP4->SaveAs("deltaP4TPCtoFull.png");
  padP4->SaveAs("deltaP4TPCtoFull.C");

  TString drawExpressionDeltaP4Const="";
  drawExpressionDeltaP4Const="[1,1,2]:";
  //drawExpressionDeltaP4Const+="%Ogridx,gridy;hisDeltaP4ConstCQPtTglTRD(0,400,80,120,0,10)(0,1)(f-mean p);:";
  drawExpressionDeltaP4Const+="%Ogridx,gridy;hisDeltaP4ConstCQPtTglTRD(0,400,80,120,0,10)(0,1)(f-grms p);:";
  //drawExpressionDeltaP4Const+="%Ogridx,gridy;hisPullP4CQPtTglTRD(0,400,80,120,0,10)(0,1)(f-mean p);:";
  drawExpressionDeltaP4Const+="%Ogridx,gridy;hisPullP4ConstCQPtTglTRD(0,400,80,120,0,10)(0,1)(f-grms p);:";
  drawExpressionDeltaP4Const+="%Ogridx,gridy;hisDeltaP4ConstCQPtTglTRD(100,300,99,101,0,10)(0)(err);:";
  drawExpressionDeltaP4Const+="%Ogridx,gridy;hisPullP4ConstCQPtTglTRD(10,90,97,102,0,10)(0)(err);:";
  TPad * padDeltaP4const = AliTreePlayer::DrawHistograms(0,hisArray,drawExpressionDeltaP4Const,keepArray, 1+2+4+8);
  padDeltaP4const->SaveAs("deltaP4TPCConsttoFullConst.png");
  padDeltaP4const->SaveAs("deltaP4TPCConsttoFullConst.C");
}


void MakeResidualDistortionMaps(){
  //
  // MakeResidualDistortionMaps
  //    Input: performanceHisto.root with sets of histograms
  //    Ouput: residualMap.root
  TFile *finput = TFile::Open("performanceHisto.root","read");
  hisArray=new TObjArray();
  TList * keys = finput->GetListOfKeys();
  for (Int_t iKey=0; iKey<keys->GetEntries(); iKey++){    
    TObject * object = finput->Get(TString::Format("%s;%d",keys->At(iKey)->GetName(),((TKey*)keys->At(iKey))->GetCycle()).Data());
    THnBase * his  = dynamic_cast<THnBase*>(object);
    if (his) hisArray->AddLast(his);
  }
  TTreeSRedirector * pcstream = new TTreeSRedirector("residualMap.root","recreate");
  // Residual histogram -> maps creation
  TPRegexp regexpHis("^(his|matchhis|qahis)");    // make residual maps for each delta histogram  
  TPRegexp regexpMatch("^matchhis");  
  TPRegexp regexpK0("hisK0");   
  //
  //
  TMatrixD projectionInfo(5,5);
  projectionInfo(0,0)=0;  projectionInfo(0,1)=0;  projectionInfo(0,2)=0;   
  projectionInfo(1,0)=1;  projectionInfo(1,1)=1;  projectionInfo(1,2)=0; 
  projectionInfo(2,0)=2;  projectionInfo(2,1)=0;  projectionInfo(2,2)=0;  
  projectionInfo(3,0)=3;  projectionInfo(3,1)=1;  projectionInfo(3,2)=0;    
  projectionInfo(4,0)=4;  projectionInfo(4,1)=0;  projectionInfo(4,2)=0;    
  for (Int_t iHis=0; iHis<hisArray->GetEntries(); iHis++){
    Int_t proj[6]={0,1,2,3,4,5};
    THn * hisInput=(THn*)hisArray->At(iHis);    
    if (hisInput->GetNdimensions()<2) continue;
    THnBase *hisProj=0;
    ::Info("MakeResidualDistortionMaps","%s\t%d\t%d\t%d",hisInput->GetName(),hisInput->GetNdimensions(), hisInput->GetNbins(), Int_t(hisInput->GetEntries()));
    if (regexpHis.Match(TString(hisInput->GetName())) && regexpK0.Match(TString(hisInput->GetName()))==0){
      Double_t fraction=(regexpMatch.Match(TString(hisInput->GetName()))>0)?0.0:0.1;
      hisInput->Print(); 
      TStatToolkit::MakeDistortionMapFast(hisInput,pcstream,projectionInfo,0,fraction);
      Int_t nDim=hisInput->GetNdimensions();
      if (nDim<2) continue;
      hisProj=hisInput->ProjectionND(nDim-1,proj);
      hisProj->SetName(TString::Format("%sProj%d",hisInput->GetName(),nDim).Data());
      TStatToolkit::MakeDistortionMapFast(hisProj,pcstream,projectionInfo,0,fraction);
      delete hisProj;
    }
  }
  // Track performance maps
  TPRegexp regexpPerf("_qPt_tgl$");  
  for (Int_t iHis=0; iHis<hisArray->GetEntries(); iHis++){
    THnBase *hisProj=0;
    Int_t proj[5]={0,1};
    if ( (regexpPerf.Match(TString(hisArray->At(iHis)->GetName()))>0) && (regexpK0.Match(TString(hisArray->At(iHis)->GetName()))==0) ){
      hisArray->At(iHis)->Print(); 
      THn * hisInput=(THn*)hisArray->At(iHis);
      Double_t fraction=(regexpMatch.Match(TString(hisInput->GetName()))>0)?0.0:0.1;
      //A side
      hisInput->GetAxis(2)->SetRangeUser(0,1);
      hisProj=hisInput->ProjectionND(2,proj);
      hisProj->SetName(TString::Format("%sASide",hisInput->GetName()).Data());
      TStatToolkit::MakeDistortionMapFast(hisProj,pcstream,projectionInfo,0,fraction);
      delete hisProj;
      //C side
      hisInput->GetAxis(2)->SetRangeUser(-1,0);
      hisProj=hisInput->ProjectionND(2,proj);
      hisProj->SetName(TString::Format("%sCSide",hisInput->GetName()).Data());
      TStatToolkit::MakeDistortionMapFast(hisProj,pcstream,projectionInfo,0,fraction);
      delete hisProj;
    }
  }
  //
  // V0 performance maps
  //
  //  TPRegexp regexpK0("hisK0");  
  projectionInfo(0,0)=0;  projectionInfo(0,1)=0;  projectionInfo(0,2)=0;   // merge pt bins
  projectionInfo(1,0)=1;  projectionInfo(1,1)=0;  projectionInfo(1,2)=0; 
  projectionInfo(2,0)=2;  projectionInfo(2,1)=1;  projectionInfo(2,2)=0;  
  //
  for (Int_t iHis=0; iHis<hisArray->GetEntries(); iHis++){
    if ( (regexpK0.Match(TString(hisArray->At(iHis)->GetName()))>0) ){
      hisArray->At(iHis)->Print(); 
      THn * hisInput=(THn*)hisArray->At(iHis);
      TStatToolkit::MakeDistortionMapFast(hisInput,pcstream,projectionInfo,0,0.1);
      Int_t proj[2]={0,1};
      THnBase * hisProj2D=hisInput->ProjectionND(2,proj);
      TStatToolkit::MakeDistortionMapFast(hisProj2D,pcstream,projectionInfo,0,0.1);
      if (hisInput->GetNdimensions()>3){ // for 4 diminsional histogram skip eta dapendence we use just A side c side
	Int_t etaBins=hisInput->GetAxis(2)->GetNbins();
	Int_t rebinEta[4]={1,1,etaBins/2,1};
	THnBase * hisAC=hisInput->Rebin(rebinEta);
	projectionInfo(2,1)=0;   
	TStatToolkit::MakeDistortionMapFast(hisAC,pcstream,projectionInfo,0,0.1);
	projectionInfo(2,1)=0; 
      }
    }
  }
  delete pcstream;
}




void DrawSummaryAndDump(TObjArray * hisArray){
  //
  //
  gStyle->SetOptStat(0);
  TObjArray *keepArray = new TObjArray(100);
  TString drawExpressionPt="";
  drawExpressionPt="[1,1,1]:";
  drawExpressionPt+="hisPtAll(0,100)(0)(errpl);hisPtITS(0,100)(0)(err);hisPtTPCOnly(0,100)(0)(er./000196528/pass4/AOD/AOD/001/FilterEvents_Trees.rootr):";
  TPad * padPt = AliTreePlayer::DrawHistograms(0,hisArray,drawExpressionPt,keepArray,1+2+4+8);
  // 
  //
  TString drawExpressionChi2TPC="";
  drawExpressionChi2TPC="[1,1,1,1]:";
  drawExpressionChi2TPC+="%Ogridx,gridy;";
  drawExpressionChi2TPC+="hisChi2TPCSecTimeA(0,40,0,180,0,10000)(0,1)(f-mean p);hisChi2TPCSecTimeC(0,40,0,180,0,10000)(0,1)(f-mean p):";
  drawExpressionChi2TPC+="%Otimex,gridx,gridy;";
  drawExpressionChi2TPC+="hisChi2TPCSecTimeA(0,40,0,180,0,10000)(0,2)(f-mean p);hisChi2TPCSecTimeC(0,40,0,180,0,10000)(0,2)(f-mean p);:";
  drawExpressionChi2TPC+="%Otimex,gridx,gridy;";
  drawExpressionChi2TPC+="hisChi2TPCSecTimeA(0,40,U0,U18,0,10000)(0,2)(f-mean p);hisChi2TPCSecTimeA(0,40,U1,U9,0,10000)(0,2)(f-mean p);hisChi2TPCSecTimeA(0,40,U10,U18,0,10000)(0,2)(f-mean p);:";
  drawExpressionChi2TPC+="%Otimex,gridx,gridy;";
  drawExpressionChi2TPC+="hisChi2TPCSecTimeC(0,40,U0,U18,0,10000)(0,2)(f-mean p);hisChi2TPCSecTimeC(0,40,U1.5,U2.5,0,10000)(0,2)(f-mean p);hisChi2TPCSecTimeC(0,40,U11.5,U12.5,0,10000)(0,2)(f-mean p);:";
  TPad * padChi2TPC = AliTreePlayer::DrawHistograms(0,hisArray,drawExpressionChi2TPC,keepArray, 1+2+4+8);
  ((TCanvas*)padChi2TPC)->SetWindowSize(1600,1000);
  //
  TString drawExpressionChi2ITS="";
  drawExpressionChi2ITS="[1,1,1,1]:";
  drawExpressionChi2ITS+="%Ogridx,gridy;";
  drawExpressionChi2ITS+="hisChi2ITSSecTimeA(0,40,0,180,0,10000)(0,1)(f-mean p);hisChi2ITSSecTimeC(0,40,0,180,0,10000)(0,1)(f-mean p):";
  drawExpressionChi2ITS+="%Otimex,gridx,gridy;";
  drawExpressionChi2ITS+="hisChi2ITSSecTimeA(0,40,0,180,0,10000)(0,2)(f-mean p);hisChi2ITSSecTimeC(0,40,0,180,0,10000)(0,2)(f-mean p);:";
  drawExpressionChi2ITS+="%Otimex,gridx,gridy;";
  drawExpressionChi2ITS+="hisChi2ITSSecTimeA(0,40,U0,U18,0,10000)(0,2)(f-mean p);hisChi2ITSSecTimeA(0,40,U1,U9,0,10000)(0,2)(f-mean p);hisChi2ITSSecTimeA(0,40,U10,U18,0,10000)(0,2)(f-mean p);:";
  drawExpressionChi2ITS+="%Otimex,gridx,gridy;";
  drawExpressionChi2ITS+="hisChi2ITSSecTimeC(0,40,U0,U18,0,10000)(0,2)(f-mean p);hisChi2ITSSecTimeC(0,40,U1.5,U2.5,0,10000)(0,2)(f-mean p);hisChi2ITSSecTimeC(0,40,U11.5,U12.5,0,10000)(0,2)(f-mean p);:";
  TPad * padChi2ITS = AliTreePlayer::DrawHistograms(0,hisArray,drawExpressionChi2ITS,keepArray,1+2+4+8);
  padChi2ITS->SetRightMargin(0.05);
  ((TCanvas*)padChi2ITS)->SetWindowSize(1600,1000);
 
}




void DrawSummary(){
  //
  // Draw summary trending information
  //
  TChain *chain = AliXRDPROOFtoolkit::MakeChainRandom("perf.list","perf",0,20,0);
  chain->SetAlias("normResolP4","pull2Time2.fY*profC14time.fY");
  chain->SetAlias("covarResolP4","profC14time.fY");
  TStatToolkit::AddMetadata(chain,"lumiBin.fY.AxisTitle","Inst. lumi (Hz/mb)");
  TStatToolkit::AddMetadata(chain,"lumiBin.fX.AxisTitle","time");
  TStatToolkit::AddMetadata(chain,"covarResolP4.AxisTitle","#sigma_{q/pt} (1/Gev/c)");
  TStatToolkit::AddMetadata(chain,"normResolP4.AxisTitle","#sigma_{q/pt} (1/Gev/c)");
  //
  TStatToolkit::AddMetadata(chain,"covarResolP4.Legend","#sigma_{q/pt} from covariance");
  TStatToolkit::AddMetadata(chain,"normResolP4.Legend","#sigma_{q/pt} x pull  ");
  //
  TCanvas * canvasDraw = new TCanvas("canvasDraw","canvasDraw",1100,600);
  canvasDraw->SetLeftMargin(0.15);
  //
  TGraph * grNorm = TStatToolkit::MakeGraphErrors(chain, "normResolP4:lumiBin.fY:profC14time.fEY","profC14time.fEY>0",25,2,0.5);
  TGraph * grCovar = TStatToolkit::MakeGraphErrors(chain, "covarResolP4:lumiBin.fY:profC14time.fEY","profC14time.fEY>0",21,4,0.5);
  grNorm->SetMinimum(0.0005);
  grNorm->SetMaximum(0.0035);
  grNorm->GetYaxis()->SetTitleOffset(1.5);
  {
    grNorm->Fit("pol1");
    grCovar->Fit("pol1");
    grNorm->Draw("ap");
    grCovar->Draw("p");
    TLegend * legend = new TLegend(0.16,0.6,0.5,0.88,"TPC+ITS q/pt resolution (|q/pt|<0.15)");
    legend->SetBorderSize(0);
    legend->AddEntry(grCovar, TStatToolkit::GetMetadata(chain,"covarResolP4.Legend")->GetTitle(),"p");
    legend->AddEntry(grNorm, TStatToolkit::GetMetadata(chain,"normResolP4.Legend")->GetTitle(),"p");
    legend->Draw();
  }
  canvasDraw->SaveAs("LHC15opass1_qptResol_qpt015_vslumi.png");
  canvasDraw->SaveAs("LHC15opass1_qptResol_qpt015_vslumi.ps");
  //
  TGraph * grNormT = TStatToolkit::MakeGraphErrors(chain, "normResolP4:lumiBin.fX:profC14time.fEY","",25,2,0.5);
  TGraph * grCovarT = TStatToolkit::MakeGraphErrors(chain, "covarResolP4:lumiBin.fX:profC14time.fEY","",21,4,0.5);
  grNormT->SetMinimum(grCovar->GetMinimum());
  grNormT->GetYaxis()->SetTitleOffset(1.5);
  //  grNormT->GetXaxis()->SetTimeDisplay();

  {
    grNormT->Draw("ap");
    grCovarT->Draw("p");
    TLegend * legend = new TLegend(0.16,0.6,0.5,0.88,"TPC+ITS q/pt resolution (|q/pt|<0.15)");
    legend->SetBorderSize(0);
    legend->AddEntry(grCovar, TStatToolkit::GetMetadata(chain,"covarResolP4.Legend")->GetTitle(),"p");
    legend->AddEntry(grNorm, TStatToolkit::GetMetadata(chain,"normResolP4.Legend")->GetTitle(),"p");
    legend->Draw();
  }
  canvasDraw->SaveAs("LHC15opass1_qptResol_qpt015_vstime.png");
  canvasDraw->SaveAs("LHC15opass1_qptResol_qpt015_vstime.ps");
  //
  //
  gPad->SetRightMargin(0.15);
  chain->Draw("profC14P4.fY:pull2P4.fX:interactionRate","","colz");
  chain->GetHistogram()->GetXaxis()->SetTitle("q/pt (1/GeV)");
  chain->GetHistogram()->GetYaxis()->SetTitle("#sigma_{q/pt} (1/GeV)");
  chain->GetHistogram()->GetYaxis()->SetTitleOffset(1.5);
  chain->GetHistogram()->Draw("colz");
  TLatex latex;
  latex.DrawTextNDC(0.2,0.8,"LHC15o pass1");
  canvasDraw->SaveAs("LHC15opass1_qptResol_qpt.png");
  canvasDraw->SaveAs("LHC15opass1_qptResol_qpt.ps");
  //
  gPad->SetRightMargin(0.15);
  chain->Draw("pull2P4.fY*profC14P4.fY:pull2P4.fX:interactionRate","","colz");
  chain->GetHistogram()->GetXaxis()->SetTitle("q/pt (1/GeV)");
  chain->GetHistogram()->GetYaxis()->SetTitle("#sigma_{q/pt} (1/GeV)");
  chain->GetHistogram()->GetYaxis()->SetTitleOffset(1.5);
  chain->GetHistogram()->Draw("colz");
  latex.DrawTextNDC(0.2,0.8,"LHC15o pass1");
  canvasDraw->SaveAs("LHC15opass1_qptResolxPull_qpt.png");
  canvasDraw->SaveAs("LHC15opass1_qptResolxPull_qpt.ps");
}


