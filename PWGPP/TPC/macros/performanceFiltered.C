/*
  gSystem->AddIncludePath("-I$AliPhysics_SRC/PWGPP/ -I$AliPhysics_SRC/OADB/");  // ? why not in the alienv,  why not available ?

  .x $NOTES/aux/NimStyle.C  
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


void SetMetadata();
void InitAnalysis();
TObjArray * FillPerfomanceHisto(Int_t maxEntries);
void SetMetadata();
void MakeResidualDistortionMaps();
//
void GetNclReport(TObjArray * hisArray,  TObjArray *keepArray );
void GetDCAReport(TObjArray * hisArray,  TObjArray *keepArray );


void performanceFiltered(Int_t maxEvents, Int_t action=0){
  //
  //   .L $NOTES/JIRA/PWGPP-221/code/performanceFiltered.C+
  //
  if (action==1) {
    MakeResidualDistortionMaps();
    return;
  }
  InitAnalysis(); 
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

}



void SetMetadata(){
  //
  chain->SetAlias("phiInner","atan2(esdTrack.fIp.Py(),esdTrack.fIp.Px()+0)");
  chain->SetAlias("secInner","9*(atan2(esdTrack.fIp.Py(),esdTrack.fIp.Px()+0)/pi)+18*(esdTrack.fIp.Py()<0)");
  //
  chain->SetAlias("deltaP0","(extInnerParamV.fP[0]-esdTrack.fP[0])");
  chain->SetAlias("pullP0","(extInnerParamV.fP[0]-esdTrack.fP[0])/sqrt(extInnerParamV.fC[0]+esdTrack.fC[0])");
  chain->SetAlias("deltaP0C","(extInnerParamC.fP[0]-esdTrack.fCp.fP[0])");
  chain->SetAlias("pullP0C","(extInnerParamC.fP[0]-esdTrack.fCp.fP[0])/sqrt(extInnerParamC.fC[0]+esdTrack.fCp.fC[0])");
  chain->SetAlias("deltaP2","(extInnerParamV.fP[2]-esdTrack.fP[2])");
  chain->SetAlias("pullP2","(extInnerParamV.fP[2]-esdTrack.fP[2])/sqrt(extInnerParamV.fC[5]+esdTrack.fC[5])");
  chain->SetAlias("deltaP2C","(extInnerParamC.fP[2]-esdTrack.fCp.fP[2])");
  chain->SetAlias("pullP2C","(extInnerParamC.fP[2]-esdTrack.fCp.fP[2])/sqrt(extInnerParamC.fC[5]+esdTrack.fCp.fC[5])");
  chain->SetAlias("deltaP4","(extInnerParam.fP[4]-esdTrack.fP[4])");
  chain->SetAlias("pullP4","(extInnerParam.fP[4]-esdTrack.fP[4])/sqrt(extInnerParam.fC[14]+esdTrack.fC[14])");
  chain->SetAlias("deltaP4C","(extInnerParamC.fP[4]-esdTrack.fCp.fP[4])");
  chain->SetAlias("pullP4C","(extInnerParamC.fP[4]-esdTrack.fCp.fP[4])/sqrt(extInnerParamC.fC[14]+esdTrack.fCp.fC[14])");
  //
  chain->SetAlias("normChi2ITS","sqrt(esdTrack.fITSchi2/esdTrack.fITSncls)");
  chain->SetAlias("normChi2TPC","esdTrack.fTPCchi2/esdTrack.fTPCncls");
  chain->SetAlias("normDCAR","esdTrack.fdTPC/sqrt(1+esdTrack.fP[4]**2)");
  chain->SetAlias("normDCAZ","esdTrack.fzTPC/sqrt(1+esdTrack.fP[4]**2)");
  chain->SetAlias("TPCASide","esdTrack.fIp.fP[1]>0");
  chain->SetAlias("TPCCSide","esdTrack.fIp.fP[1]<0");
  chain->SetAlias("TPCCross","esdTrack.fIp.fP[1]*esdTrack.fIp.fP[3]<0");
  chain->SetAlias("qPt","esdTrack.fP[4]");
  chain->SetAlias("tgl","esdTrack.fP[3]");
  chain->SetAlias("alphaV","esdTrack.fAlpha");
  //
  chain->SetAlias("ITSOn","((esdTrack.fFlags&0x1)>0)");
  chain->SetAlias("TPCOn","((esdTrack.fFlags&0x10)>0)");
  chain->SetAlias("ITSRefit","((esdTrack.fFlags&0x4)>0)");
  chain->SetAlias("TPCRefit","((esdTrack.fFlags&0x40)>0)");
  chain->SetAlias("TOFOn","((esdTrack.fFlags&0x2000)>0)");
  chain->SetAlias("TRDOn","((esdTrack.fFlags&0x400)>0)");
  chain->SetAlias("ITSOn0","esdTrack.fITSncls>4&&esdTrack.HasPointOnITSLayer(0)&&esdTrack.HasPointOnITSLayer(1)");
  chain->SetAlias("nclCut","(esdTrack.GetTPCClusterInfo(3,1)+esdTrack.fTRDncls)>140-5*(abs(esdTrack.fP[4]))");
  chain->SetAlias("IsPrim4","abs(esdTrack.fD/sqrt(esdTrack.fCdd))<4");
  

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
  chain = AliXRDPROOFtoolkit::MakeChainRandom("filtered.list","highPt",0,40000);
  chainV0 = AliXRDPROOFtoolkit::MakeChainRandom("filtered.list","V0s",0,40000);
  AliAnalysisTaskFilteredTree::SetDefaultAliasesHighPt(chain);
  AliAnalysisTaskFilteredTree::SetDefaultAliasesV0(chainV0);
  SetMetadata();
  //
  Int_t selected = chain->Draw("ntracks:mult:evtTimeStamp","","goff",100000);
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
  TString defaultCut="esdTrack.fTPCncls>60&&esdTrack.IsOn(0x1)>0";
  TString hisString="";
  {
    hisString+="esdTrack.Pt():#esdTrack.fTPCncls>60>>hisPtAll(100,0,30);"; 
    hisString+="esdTrack.Pt():#(esdTrack.fFlags&0x4)>0>>hisPtITS(100,1,30);";    
    hisString+="esdTrack.fIp.Pt():#(esdTrack.fFlags&0x4)>0>>hisPtTPCOnly(100,1,30);";  
    // Kinematic histograms
    hisString+="esdTrack.fP[4]:esdTrack.fP[3]:secInner:#esdTrack.fTPCncls>60>>hisQptTglSecAll(40,-2,2,10,-1,1,90,0,18);"; 
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
    hisString+=TString::Format("normChi2ITS:esdTrack.fP[3]:ntracks>>hisChi2ITSP3NTracks(50,0,10,10,-1,1,10,0,%d);",ntracksEnd);
    hisString+=TString::Format("normChi2ITS:esdTrack.fP[3]:mult>>hisChi2ITSP3Mult(50,0,10,10,-1,1,10,0,%d);",multEnd);
    // Chi2 histograms TPC    
    hisString+=TString::Format("normChi2TPC:secInner:evtTimeStamp:#TPCASide>>hisChi2TPCSecTimeA(50,0,10,180,0,18,%s);",timeRange.Data());
    hisString+=TString::Format("normChi2TPC:secInner:evtTimeStamp:#TPCCSide>>hisChi2TPCSecTimeC(50,0,10,180,0,18,%s);",timeRange.Data());
    hisString+=TString::Format("normChi2TPC:secInner:evtTimeStamp:#TPCCross>>hisChi2TPCSecTimeCross(50,0,10,180,0,18,%s);",timeRange.Data());
    hisString+=TString::Format("normChi2TPC:esdTrack.fP[3]:ntracks>>hisChi2TPCP3NTracks(50,0,10,10,-1,1,10,0,%d);",ntracksEnd);
    hisString+=TString::Format("normChi2TPC:esdTrack.fP[3]:mult>>hisChi2TPCP3Mult(50,0,10,10,-1,1,10,0,%d);",multEnd);
    // DCAr histograms normalized to Pt
    hisString+=TString::Format("normDCAR:secInner:evtTimeStamp:#TPCASide&&abs(esdTrack.fD)<0.2>>hisTPCDCARSecTimeA(100,-3,3,180,0,18,%s);",timeRange.Data());
    hisString+=TString::Format("normDCAR:secInner:evtTimeStamp:#TPCCSide&&abs(esdTrack.fD)<0.2>>hisTPCDCARSecTimeC(100,-3,3,180,0,18,%s);",timeRange.Data());
    hisString+=TString::Format("normDCAR:secInner:evtTimeStamp:#TPCCross&&abs(esdTrack.fD)<0.2>>hisTPCDCARSecTimeCross(100,-3,3,180,0,18,%s);",timeRange.Data());
    // DCAz histograms normalized to Pt
    hisString+=TString::Format("normDCAZ:secInner:evtTimeStamp:#TPCASide&&abs(esdTrack.fD)<0.2>>hisTPCDCAZSecTimeA(100,-3,3,180,0,18,%s);",timeRange.Data());
    hisString+=TString::Format("normDCAZ:secInner:evtTimeStamp:#TPCCSide&&abs(esdTrack.fD)<0.2>>hisTPCDCAZSecTimeC(100,-3,3,180,0,18,%s);",timeRange.Data());
    hisString+=TString::Format("normDCAZ:secInner:evtTimeStamp:#TPCCross&&abs(esdTrack.fD)<0.2>>hisTPCDCAZSecTimeCross(100,-3,3,180,0,18,%s);",timeRange.Data());
  }

  {
    // deltaP0 and pullP0 histograms  (TPC+TRD - ITS+TPF+TRD)
    hisString+=TString::Format("deltaP0:qPt:tgl:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn0&&nclCut>>hisDeltaP0CQPtTglAll(400,-3.0,3.0,200,-5,5,10,-1,1);");
    hisString+=TString::Format("deltaP0:qPt:tgl:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn0&&TRDOn&&nclCut>>hisDeltaP0CQPtTglTRD(400,-3.0,3.0,200,-5,5,10,-1,1);");
    hisString+=TString::Format("pullP0:qPt:tgl:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn0&&nclCut>>hisPullP0CQPtTglAll(100,-6,6,200,-5,5,10,-1,1);");
    hisString+=TString::Format("pullP0:qPt:tgl:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn0&&TRDOn&&nclCut>>hisPullP0CQPtTglTRD(100,-6,6,200,-5,5,10,-1,1);");
    hisString+=TString::Format("deltaP0C:qPt:tgl:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn0&&nclCut>>hisDeltaP0ConstCQPtTglAll(400,-3.0,3.0,200,-5,5,10,-1,1);");
    hisString+=TString::Format("deltaP0C:qPt:tgl:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn0&&TRDOn&&nclCut>>hisDeltaP0ConstCQPtTglTRD(400,-3.0,3.0,200,-5,5,10,-1,1);");
    hisString+=TString::Format("pullP0C:qPt:tgl:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn0&&nclCut>>hisPullP0ConstCQPtTglAll(100,-6,6,200,-5,5,10,-1,1);");
    hisString+=TString::Format("pullP0C:qPt:tgl:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn0&&TRDOn&&nclCut>>hisPullP0ConstCQPtTglTRD(100,-6,6,200,-5,5,10,-1,1);");
    // deltaP2 and pullP2 histograms  (TPC+TRD - ITS+TPF+TRD)
    hisString+=TString::Format("deltaP2:qPt:tgl:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn0&&nclCut>>hisDeltaP2CQPtTglAll(400,-0.01,0.01,200,-5,5,10,-1,1);");
    hisString+=TString::Format("deltaP2:qPt:tgl:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn0&&TRDOn&&nclCut>>hisDeltaP2CQPtTglTRD(400,-0.01,0.01,200,-5,5,10,-1,1);");
    hisString+=TString::Format("pullP2:qPt:tgl:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn0&&nclCut>>hisPullP2CQPtTglAll(100,-6,6,200,-5,5,10,-1,1);");
    hisString+=TString::Format("pullP2:qPt:tgl:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn0&&TRDOn&&nclCut>>hisPullP2CQPtTglTRD(100,-6,6,200,-5,5,10,-1,1);");
    hisString+=TString::Format("deltaP2C:qPt:tgl:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn0&&nclCut>>hisDeltaP2ConstCQPtTglAll(400,-0.01,0.01,200,-5,5,10,-1,1);");
    hisString+=TString::Format("deltaP2C:qPt:tgl:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn0&&TRDOn&&nclCut>>hisDeltaP2ConstCQPtTglTRD(400,-0.01,0.01,200,-5,5,10,-1,1);");
    hisString+=TString::Format("pullP2C:qPt:tgl:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn0&&nclCut>>hisPullP2ConstCQPtTglAll(100,-6,6,200,-5,5,10,-1,1);");
    hisString+=TString::Format("pullP2C:qPt:tgl:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn0&&TRDOn&&nclCut>>hisPullP2ConstCQPtTglTRD(100,-6,6,200,-5,5,10,-1,1);");
    // deltaP4 and pullP4 histograms (TPC+TRD - ITS+TPF+TRD)
    hisString+=TString::Format("deltaP4:qPt:tgl:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn0&&nclCut>>hisDeltaP4CQPtTglAll(400,-0.05,0.05,200,-5,5,10,-1,1);");
    hisString+=TString::Format("deltaP4:qPt:tgl:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn0&&TRDOn&&nclCut>>hisDeltaP4CQPtTglTRD(400,-0.05,0.05,200,-5,5,10,-1,1);");
    hisString+=TString::Format("pullP4:qPt:tgl:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn0&&nclCut>>hisPullP4CQPtTglAll(100,-6,6,200,-5,5,10,-1,1);");
    hisString+=TString::Format("pullP4:qPt:tgl:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn0&&TRDOn&&nclCut>>hisPullP4CQPtTglTRD(100,-6,6,200,-5,5,10,-1,1);");
    hisString+=TString::Format("deltaP4C:qPt:tgl:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn0&&nclCut>>hisDeltaP4ConstCQPtTglAll(400,-0.05,0.05,200,-5,5,10,-1,1);");
    hisString+=TString::Format("deltaP4C:qPt:tgl:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn0&&TRDOn&&nclCut>>hisDeltaP4ConstCQPtTglTRD(400,-0.05,0.05,200,-5,5,10,-1,1);");
    hisString+=TString::Format("pullP4C:qPt:tgl:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn0&&nclCut>>hisPullP4ConstCQPtTglAll(100,-6,6,200,-5,5,10,-1,1);");
    hisString+=TString::Format("pullP4C:qPt:tgl:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn0&&TRDOn&&nclCut>>hisPullP4ConstCQPtTglTRD(100,-6,6,200,-5,5,10,-1,1);");
    // deltaP4 and pullP4 histograms (TPC+TRD - ITS+TPF+TRD)
    hisString+=TString::Format("sqrt(esdTrack.fC[14]):qPt:tgl:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn0&&nclCut>>hisCovarP4CQPtTglAll(400,-0.05,0.05,200,-2.5,2.5,10,-1,1);");
    hisString+=TString::Format("sqrt(esdTrack.fC[14]):qPt:tgl:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn0&&TRDOn&&nclCut>>hisCovarP4CQPtTglTRD(400,-0.05,0.05,200,-2.5,2.5,10,-1,1);");
    hisString+=TString::Format("sqrt(esdTrack.fCp.fC[14]):qPt:tgl:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn0&&nclCut>>hisCovarP4ConstCQPtTglAll(400,-0.05,0.05,200,-2.5,2.5,10,-1,1);");
    hisString+=TString::Format("sqrt(esdTrack.fCp.fC[14]):qPt:tgl:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn0&&TRDOn&&nclCut>>hisCovarP4ConstCQPtTglTRD(400,-0.05,0.05,200,-2.5,2.5,10,-1,1);");
  }
  // Residual calibration histogramming
  {
    // deltaP0 and pullP0 histograms  (TPC+TRD - ITS+TPF+TRD)
    hisString+=TString::Format("deltaP0:qPt:tgl:alphaV:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn0&&nclCut>>hisDeltaP0CQPtTglAlphaVAll(100,-3.0,3.0,12,-3,3,10,-1,1,90,-3.145,3.145);");
    hisString+=TString::Format("deltaP0:qPt:tgl:alphaV:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn0&&TRDOn&&nclCut>>hisDeltaP0CQPtTglAlphaVTRD(100,-3.0,3.0,12,-3,3,10,-1,1,90,-3.145,3.145);");
    hisString+=TString::Format("pullP0:qPt:tgl:alphaV:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn0&&nclCut>>hisPullP0CQPtTglAlphaVAll(100,-6,6,12,-3,3,10,-1,1,90,-3.145,3.145);");
    hisString+=TString::Format("pullP0:qPt:tgl:alphaV:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn0&&TRDOn&&nclCut>>hisPullP0CQPtTglAlphaVTRD(100,-6,6,12,-3,3,10,-1,1,90,-3.145,3.145);");
    hisString+=TString::Format("deltaP0C:qPt:tgl:alphaV:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn0&&nclCut>>hisDeltaP0ConstCQPtTglAlphaVAll(100,-3.0,3.0,12,-3,3,10,-1,1,90,-3.145,3.145);");
    hisString+=TString::Format("deltaP0C:qPt:tgl:alphaV:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn0&&TRDOn&&nclCut>>hisDeltaP0ConstCQPtTglAlphaVTRD(100,-3.0,3.0,12,-3,3,10,-1,1,90,-3.145,3.145);");
    hisString+=TString::Format("pullP0C:qPt:tgl:alphaV:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn0&&nclCut>>hisPullP0ConstCQPtTglAlphaVAll(100,-6,6,12,-3,3,10,-1,1,90,-3.145,3.145);");
    hisString+=TString::Format("pullP0C:qPt:tgl:alphaV:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn0&&TRDOn&&nclCut>>hisPullP0ConstCQPtTglAlphaVTRD(100,-6,6,12,-3,3,10,-1,1,90,-3.145,3.145);");
    // deltaP2 and pullP2 histograms  (TPC+TRD - ITS+TPF+TRD)
    hisString+=TString::Format("deltaP2:qPt:tgl:alphaV:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn0&&nclCut>>hisDeltaP2CQPtTglAlphaVAll(100,-0.01,0.01,12,-3,3,10,-1,1,90,-3.145,3.145);");
    hisString+=TString::Format("deltaP2:qPt:tgl:alphaV:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn0&&TRDOn&&nclCut>>hisDeltaP2CQPtTglAlphaVTRD(100,-0.01,0.01,12,-3,3,10,-1,1,90,-3.145,3.145);");
    hisString+=TString::Format("pullP2:qPt:tgl:alphaV:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn0&&nclCut>>hisPullP2CQPtTglAlphaVAll(100,-6,6,12,-3,3,10,-1,1,90,-3.145,3.145);");
    hisString+=TString::Format("pullP2:qPt:tgl:alphaV:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn0&&TRDOn&&nclCut>>hisPullP2CQPtTglAlphaVTRD(100,-6,6,12,-3,3,10,-1,1,90,-3.145,3.145);");
    hisString+=TString::Format("deltaP2C:qPt:tgl:alphaV:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn0&&nclCut>>hisDeltaP2ConstCQPtTglAlphaVAll(100,-0.01,0.01,12,-3,3,10,-1,1,90,-3.145,3.145);");
    hisString+=TString::Format("deltaP2C:qPt:tgl:alphaV:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn0&&TRDOn&&nclCut>>hisDeltaP2ConstCQPtTglAlphaVTRD(100,-0.01,0.01,12,-3,3,10,-1,1,90,-3.145,3.145);");
    hisString+=TString::Format("pullP2C:qPt:tgl:alphaV:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn0&&nclCut>>hisPullP2ConstCQPtTglAlphaVAll(100,-6,6,12,-3,3,10,-1,1,90,-3.145,3.145);");
    hisString+=TString::Format("pullP2C:qPt:tgl:alphaV:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn0&&TRDOn&&nclCut>>hisPullP2ConstCQPtTglAlphaVTRD(100,-6,6,12,-3,3,10,-1,1,90,-3.145,3.145);");
    // deltaP4 and pullP4 histograms (TPC+TRD - ITS+TPF+TRD)
    hisString+=TString::Format("deltaP4:qPt:tgl:alphaV:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn0&&nclCut>>hisDeltaP4CQPtTglAlphaVAll(100,-0.05,0.05,12,-3,3,10,-1,1,90,-3.145,3.145);");
    hisString+=TString::Format("deltaP4:qPt:tgl:alphaV:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn0&&TRDOn&&nclCut>>hisDeltaP4CQPtTglAlphaVTRD(100,-0.05,0.05,12,-3,3,10,-1,1,90,-3.145,3.145);");
    hisString+=TString::Format("pullP4:qPt:tgl:alphaV:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn0&&nclCut>>hisPullP4CQPtTglAlphaVAll(100,-6,6,12,-3,3,10,-1,1,90,-3.145,3.145);");
    hisString+=TString::Format("pullP4:qPt:tgl:alphaV:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn0&&TRDOn&&nclCut>>hisPullP4CQPtTglAlphaVTRD(100,-6,6,12,-3,3,10,-1,1,90,-3.145,3.145);");
    hisString+=TString::Format("deltaP4C:qPt:tgl:alphaV:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn0&&nclCut>>hisDeltaP4ConstCQPtTglAlphaVAll(100,-0.05,0.05,12,-3,3,10,-1,1,90,-3.145,3.145);");
    hisString+=TString::Format("deltaP4C:qPt:tgl:alphaV:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn0&&TRDOn&&nclCut>>hisDeltaP4ConstCQPtTglAlphaVTRD(100,-0.05,0.05,12,-3,3,10,-1,1,90,-3.145,3.145);");
    hisString+=TString::Format("pullP4C:qPt:tgl:alphaV:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn0&&nclCut>>hisPullP4ConstCQPtTglAlphaVAll(100,-6,6,12,-3,3,10,-1,1,90,-3.145,3.145);");
    hisString+=TString::Format("pullP4C:qPt:tgl:alphaV:#IsPrim4&&TPCOn&&ITSRefit&&ITSOn0&&TRDOn&&nclCut>>hisPullP4ConstCQPtTglAlphaVTRD(100,-6,6,12,-3,3,10,-1,1,90,-3.145,3.145);");
  }
  TString hisV0String="";
  chainV0->SetAlias("K0cut","abs(K0PIDPull)<2&&(abs(LPull)>2.&&abs(ALPull)>2.)");
  {
    // K0 performance
    hisV0String+="K0Delta:1/v0.Pt():tglV0:#K0cut>>hisK0DMassQPtTgl(100,-0.03,0.03,80,0,2,5,10,-1,1);";  
    hisV0String+="K0Pull:1/v0.Pt():tglV0:#K0cut>>hisK0PullQPtTgl(100,-6.0,6.0,80,0,2,5,10,-1,1);";
    // K0 resolution/maps - in respec to sector edge
    hisV0String+="K0Delta:1/v0.Pt():tglV0:dalphaV0:#K0cut>>hisK0DMassQPtTglDSec(100,-0.03,0.03,10,0,1,10,-1,1,10,0.0,0.35);";  
    hisV0String+="K0Pull:1/v0.Pt():tglV0:dalphaV0:#K0cut>>hisK0PullQPtTglDSec(100,-6.0,6.0,10,0,1,10,-1,1,10,0.0,0.35);";
    // K0 resolution/maps - 
    hisV0String+="K0Delta:1/v0.Pt():tglV0:alphaV0:#K0cut>>hisK0DMassQPtTglAlpha(100,-0.03,0.03,10,0,1,5,-1,1,18,-3.1415,3.1415);";  
    hisV0String+="K0Pull:1/v0.Pt():tglV0:alphaV0:#K0cut>>hisK0PullQPtTglAlpha(100,-6.0,6.0,10,0,1,5,-1,1,18,-3.1415,3.1415);";
  }  
  //
  TStopwatch timer;

  timer.Start();
  hisArrayV0 = AliTreePlayer::MakeHistograms(chainV0, hisV0String, "",0,maxEntries,10000000,15);
  timer.Print();

  timer.Start();
  hisArray = AliTreePlayer::MakeHistograms(chain, hisString, defaultCut,0,maxEntries,10000000,15);
  timer.Print();


  (*pcstream).GetFile()->cd();
  //  hisArray->Write("perfArray",  TObjArray::kSingleKey);
  hisArray->Write("perfArray");
  hisArrayV0->Write("perfArrayV0");
  (*pcstream).GetFile()->Flush();
  return hisArray;
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
  //
  
  TPad * padChi2ITSSectorTime = AliTreePlayer::DrawHistograms(0,hisArray,drawExpressionChi2ITSSectorTime,keepArray, 1+2+4+8);
  ((TCanvas*)padChi2ITSSectorTime)->SetWindowSize(1600,1000);
  padChi2ITSSectorTime->SaveAs("itsChi2SectorTimeReport.png");
  padChi2ITSSectorTime->SaveAs("itsChi2SectorTimeReport.pdf");



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
    TObject * o = finput->Get(TString::Format("%s;%d",keys->At(iKey)->GetName(),((TKey*)keys->At(iKey))->GetCycle()).Data());
    hisArray->AddLast(o);
  }
  TTreeSRedirector * pcstream = new TTreeSRedirector("residualMap.root","recreate");
  // Residual histogram -> maps creation
  TPRegexp regexpDelta("(eltaP|ullP|covar).*phaV");    // make residual maps for each delta histogram  
  TMatrixD projectionInfo(4,5);
  projectionInfo(0,0)=0;  projectionInfo(0,1)=0;  projectionInfo(0,2)=0;   
  projectionInfo(1,0)=1;  projectionInfo(1,1)=1;  projectionInfo(1,2)=0; 
  projectionInfo(2,0)=2;  projectionInfo(2,1)=0;  projectionInfo(2,2)=0;  
  projectionInfo(3,0)=3;  projectionInfo(3,1)=1;  projectionInfo(3,2)=0;    
  for (Int_t iHis=0; iHis<hisArray->GetEntries(); iHis++){
    if (regexpDelta.Match(TString(hisArray->At(iHis)->GetName()))){
      hisArray->At(iHis)->Print(); 
      THn * hisInput=(THn*)hisArray->At(iHis);
      TStatToolkit::MakeDistortionMapFast(hisInput,pcstream,projectionInfo,1,0.1);
    }
  }
  // Track performance maps
  TPRegexp regexpPerf("(eltaP|ullP|hisCovar).*");  
  for (Int_t iHis=0; iHis<hisArray->GetEntries(); iHis++){
    if ( (regexpPerf.Match(TString(hisArray->At(iHis)->GetName()))>0) && (regexpDelta.Match(TString(hisArray->At(iHis)->GetName()))==0) ){
      hisArray->At(iHis)->Print(); 
      THn * hisInput=(THn*)hisArray->At(iHis);
      TStatToolkit::MakeDistortionMapFast(hisInput,pcstream,projectionInfo,0,0.1);
    }
  }
  //
  // V0 performance maps
  //
  TPRegexp regexpK0("hisK0");  
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


