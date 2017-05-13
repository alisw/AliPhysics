/*
  .x $NOTES/aux/NimStyle.C(1) 
  NimStyleBig()
  .L $AliPhysics_SRC/PWGPP/TPC/macros/comparePerformanceMaps.C+
   InitTrees();  
  
 

 */
#include "TFile.h"
#include "TTree.h"
#include "TPRegexp.h"
#include "TSystem.h"
#include "TStatToolkit.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TStyle.h"

TTree *  InitMapTree(TPRegexp regExp,  TPRegexp notReg, TString axisAlias,  TString axisTitle);
TObjArray * samples=new TObjArray(1000);


TPRegexp dummy("^!");
TPRegexp regTreeK0Qpt("hisK0.*proj_0_1");
TPRegexp regTreeNotK0Qpt("hisK0.*(Alpha|DSec)");
TPRegexp regTreeK0QptDSec("hisK0.*QPtTglDSec_1_1_5_1Dist");
TPRegexp regTreeK0Alpha("hisK0.*AlphaDist");
//TPRegexp regTreeDelta("(his|qahis).*(Delta|Pull|Covar|ncl|Chi2).*_tglDist");  // regular expression for standard trees
TPRegexp regTreeDelta("^(his|qahis).*(Delta|Pull|Covar|ncl|Chi2|covar|delta|pull).*_tglDist$");  // regular expression for standard trees
TPRegexp regTreeNotDeltaInt("his.*(lpha|DSec)");
//
TPRegexp regCovar("hisCovar");
TPRegexp regTreeDeltaAlpha("(his|qahis)(Delta|Pull|Covar|ncl|Chi2|covar|delta|pull).*alphaVDist$");
TPRegexp regTreeDeltaDAlphaQ("(his|qahis)(Delta|Pull|Covarr|ncl|Chi2|covar|delta|pull).*_dalphaQDist$");

TTree * treeDelta0,  *treeK0proj_0_1,  *treeK0QptDSec, * treeCovar, * treeDeltaAlpha, *treeDeltaDAlphaQ, *treeK0Alpha;
TCanvas *canvasDraw=0;


void makeCanvas(){
  canvasDraw = new TCanvas("xxx","xxx",1100,500);
  canvasDraw->SetRightMargin(0.15);
  canvasDraw->SetBottomMargin(0.15);
  canvasDraw->SetLeftMargin(0.1);
  gStyle->SetTitleOffset(0.8,"Z");
  gStyle->SetTitleOffset(1.0,"X");
  gStyle->SetTitleOffset(0.6,"Y");
}

void InitTrees(){
  treeDelta0= InitMapTree(regTreeDelta,regTreeNotDeltaInt,"qPt:tgl","q/p_{t}(1/GeV):unit");
  treeK0proj_0_1= InitMapTree(regTreeK0Qpt,regTreeNotK0Qpt , "mpt:tgl","1/p_{t} (1/GeV): tan(#lambda)");
  treeK0QptDSec= InitMapTree(regTreeK0QptDSec,dummy , "mpt:side:dsec","1/p_{t} (1/GeV): side:dsec");
  treeCovar= InitMapTree(regCovar,dummy , "qPt:tgl","q/p_{t} (1/GeV): tan(#lambda)");
  treeDeltaAlpha= InitMapTree(regTreeDeltaAlpha,dummy , "qPt:tgl:alpha","q/p_{t} (1/GeV):tan(#lambda):alpha");
  treeDeltaDAlphaQ= InitMapTree(regTreeDeltaDAlphaQ,dummy , "qPt:tgl:alpha","q/p_{t} (1/GeV):tan(#lambda):dalphaQ");
  treeK0Alpha= InitMapTree(regTreeK0Alpha,dummy , "qPt:tgl:alpha","q/p_{t} (1/GeV):tan(#lambda):alpha");

}


TTree *  InitMapTree(TPRegexp regExp, TPRegexp notReg,  TString axisAlias,  TString axisTitle){
  //
  TObjArray * residualMapList=  gSystem->GetFromPipe("cat residualMap.list").Tokenize("\n");
  Int_t nFiles=residualMapList->GetEntries()/2;
  TTree * treeBase =0; 
  TObjArray *arrayAxisAlias=axisAlias.Tokenize(":");
  TObjArray *arrayAxisTitle=axisTitle.Tokenize(":");
  for (Int_t iFile=0; iFile<nFiles; iFile++){
    TString name0=residualMapList->At(iFile*2)->GetName();
    samples->AddAt(new TObjString(name0),iFile);
    TFile * finput = TFile::Open(residualMapList->At(iFile*2+1)->GetName());
    if (finput==NULL){
      ::Error("MakeResidualDistortionReport","Invalid file name %s",residualMapList->At(iFile*2+1)->GetName()); 
      continue;
    }
    TList * keys = finput->GetListOfKeys();
    for (Int_t iKey=0; iKey<keys->GetEntries(); iKey++){   
      if (regExp.Match(keys->At(iKey)->GetName())==0) continue;
      if (notReg.Match(keys->At(iKey)->GetName())!=0) continue;
      TTree * tree = (TTree*)finput->Get(keys->At(iKey)->GetName()); // better to use dynamic cast
      if (treeBase==NULL) {
	TFile * finput2 = TFile::Open(residualMapList->At(iFile*2+1)->GetName());
	treeBase= (TTree*)finput2->Get(keys->At(iKey)->GetName());
      }
      treeBase->AddFriend(tree,TString::Format("%s.%s",name0.Data(),keys->At(iKey)->GetName()).Data());
      Int_t entriesF=tree->GetEntries();
      Int_t entriesB=treeBase->GetEntries();
      if (entriesB==entriesF){
	::Info("InitMapTree", "%s\t%s.%s:\t%d\t%d",treeBase->GetName(),name0.Data(),keys->At(iKey)->GetName(),entriesB, entriesF);
      }else{
	::Error("InitMapTree","%s\t%s.%s:\t%d\t%d",treeBase->GetName(),name0.Data(),keys->At(iKey)->GetName(),entriesB, entriesF);
      }
    }
  }
  if (treeBase){
    for (Int_t iAxis=0; iAxis<arrayAxisAlias->GetEntries(); iAxis++){
      treeBase->SetAlias(arrayAxisAlias->At(iAxis)->GetName(), TString::Format("%s%dCenter",arrayAxisAlias->At(iAxis)->GetName(),iAxis+1).Data());
      TStatToolkit::AddMetadata(treeBase,TString::Format("%s.AxisTitle",arrayAxisAlias->At(iAxis)->GetName()).Data(),  arrayAxisTitle->At(iAxis)->GetName());
    }
    treeBase->SetMarkerStyle(25);
    treeBase->SetMarkerSize(0.5);
  }
  return treeBase;
}



void MakeK0DefaultReport(){
  //
  // 0. Init Variables
  TCut selection="entries>100&&mpt<1";
  gStyle->SetOptTitle(0);
  TMultiGraph* mGraph =0;
  TCanvas *canvas =  new TCanvas("canvas","canvas",1200, 800);
  TLegend *legend=0;
  TString drawString="";
  canvas->SetGrid(1,1);
  // 1. Set metadata
  Int_t nSamples=samples->GetEntries();
  for (Int_t i=0; i<nSamples; i++){
    treeK0proj_0_1->SetAlias(Form("%s.dMassK0.rms",samples->At(i)->GetName()),Form("1000*%s.hisK0DMassQPtTgl_proj_0_1Dist.rmsG+0",samples->At(i)->GetName()));
    treeK0proj_0_1->SetAlias(Form("%s.dMassK0.err",samples->At(i)->GetName()),Form("1000*%s.hisK0DMassQPtTgl_proj_0_1Dist.rmsG/sqrt(%s.hisK0DMassQPtTgl_proj_0_1Dist.entries+0)",samples->At(i)->GetName(),samples->At(i)->GetName()));
    treeK0proj_0_1->SetAlias(Form("%s.dMassK0.mean",samples->At(i)->GetName()),Form("1000*%s.hisK0DMassQPtTgl_proj_0_1Dist.meanG+0",samples->At(i)->GetName()));
    treeK0proj_0_1->SetAlias(Form("%s.dMassK0.meanErr",samples->At(i)->GetName()),Form("1000*%s.hisK0DMassQPtTgl_proj_0_1Dist.rmsG/sqrt(%s.hisK0DMassQPtTgl_proj_0_1Dist.entries+0)",samples->At(i)->GetName(),samples->At(i)->GetName()));
    treeK0proj_0_1->SetAlias(Form("%s.dMassK0.pull",samples->At(i)->GetName()),Form("%s.hisK0PullQPtTgl_proj_0_1Dist.rmsG+0",samples->At(i)->GetName()));
    treeK0proj_0_1->SetAlias(Form("%s.dMassK0.pullErr",samples->At(i)->GetName()),Form("%s.hisK0PullQPtTgl_proj_0_1Dist.rmsG/sqrt(%s.hisK0PullQPtTgl_proj_0_1Dist.entries+0)",samples->At(i)->GetName(),samples->At(i)->GetName()));
  }
  for (Int_t i=0; i<nSamples; i++){
    TStatToolkit::AddMetadata(treeK0proj_0_1,TString::Format("%s.dMassK0.rms.AxisTitle",samples->At(i)->GetName()).Data(),"#sigma_{M} (MeV/c)"); 
    TStatToolkit::AddMetadata(treeK0proj_0_1,TString::Format("%s.dMassK0.rms.Legend",samples->At(i)->GetName()).Data(),samples->At(i)->GetName()); 
    TStatToolkit::AddMetadata(treeK0proj_0_1,TString::Format("%s.dMassK0.pull.AxisTitle",samples->At(i)->GetName()).Data(),"pull  = RMS(#frac{#Delta_{M}}{#sigma^{exp}_{M}})"); 
    TStatToolkit::AddMetadata(treeK0proj_0_1,TString::Format("%s.dMassK0.pull.Legend",samples->At(i)->GetName()).Data(),samples->At(i)->GetName()); 
    TStatToolkit::AddMetadata(treeK0proj_0_1,TString::Format("%s.dMassK0.mean.AxisTitle",samples->At(i)->GetName()).Data(),"#delta_{M} (MeV/c)"); 
    TStatToolkit::AddMetadata(treeK0proj_0_1,TString::Format("%s.dMassK0.mean.Legend",samples->At(i)->GetName()).Data(),samples->At(i)->GetName()); 
  }
  // Fig 1.1) Inv mss resol
  {
    gPad->Clear(); 
    legend = new TLegend(0.2,0.15,0.8,0.3,"K^{0} Inv.mass resolution (gauss fit)");
    legend->SetBorderSize(0);legend->SetNColumns(2);
    drawString="";
    for (Int_t i=0; i<nSamples; i++) {drawString+=Form("%s.dMassK0.rms",samples->At(i)->GetName()); drawString+=(i==nSamples-1)?":":";";}
    drawString+="mpt:";
    for (Int_t i=0; i<nSamples; i++) {drawString+=Form("%s.dMassK0.err",samples->At(i)->GetName()); drawString+=(i==nSamples-1)?"":";";}
    mGraph = TStatToolkit::MakeMultGraph(treeK0proj_0_1,"K0",drawString.Data(),selection,"25;21;20;24","1;2;4;6",kFALSE,1,10,legend);
    mGraph->Draw();
    legend->Draw();
    canvas->SaveAs("k0RMSMassMPt.png");
    canvas->SaveAs("k0RMSMassMPt.pdf");
  }
  // Inv mass delta
  {
    gPad->Clear(); 
    legend = new TLegend(0.2,0.15,0.6,0.3,"K^{0} #Delta Inv.mass");
    legend->SetFillStyle(4050);
    legend->SetBorderSize(0);legend->SetNColumns(2);
    drawString="";
    for (Int_t i=0; i<nSamples; i++) {drawString+=Form("%s.dMassK0.mean",samples->At(i)->GetName()); drawString+=(i==nSamples-1)?":":";";}
    drawString+="mpt:";
    for (Int_t i=0; i<nSamples; i++) {drawString+=Form("%s.dMassK0.meanErr",samples->At(i)->GetName()); drawString+=(i==nSamples-1)?"":";";}
    mGraph = TStatToolkit::MakeMultGraph(treeK0proj_0_1,"K0",drawString.Data(),selection,"25;21;20;24","1;2;4;6",kFALSE,1,8,legend);
    mGraph->Draw();
    legend->Draw();
    canvas->SaveAs("k0DeltaMassMPt.png");
    canvas->SaveAs("k0DeltaMassMPt.pdf");
  }
  //
  // Inv mass pull - rms 
  {
    gPad->Clear(); 
    legend = new TLegend(0.2,0.15,0.6,0.3,"K^{0} #Delta Inv.mass");
    legend->SetFillStyle(4050);
    legend->SetBorderSize(0);legend->SetNColumns(2);
    drawString="";
    for (Int_t i=0; i<nSamples; i++) {drawString+=Form("%s.dMassK0.pull",samples->At(i)->GetName()); drawString+=(i==nSamples-1)?":":";";}
    drawString+="mpt:";
    for (Int_t i=0; i<nSamples; i++) {drawString+=Form("%s.dMassK0.pullErr",samples->At(i)->GetName()); drawString+=(i==nSamples-1)?"":";";}
    mGraph = TStatToolkit::MakeMultGraph(treeK0proj_0_1,"K0",drawString.Data(),selection,"25;21;20;24","1;2;4;6",kFALSE,1,8,legend);
    mGraph->Draw();
    legend->Draw();
    canvas->SaveAs("k0PullMassMPt.png");
    canvas->SaveAs("k0PullMassMPt.pdf");
  }
  //
}

void MakeK0DSectorReport(){
  //
  // 0. Init Variables
  TCut selection="entries>100&&abs(mpt-0.11)<0.05&&side>0";
  gStyle->SetOptTitle(0);
  TMultiGraph* mGraph =0;
  TCanvas *canvas =  new TCanvas("canvas","canvas",1200, 800);
  TLegend *legend=0;
  TString drawString="";
  canvas->SetGrid(1,1);
  // 1. Set metadata
  Int_t nSamples=samples->GetEntries();
  for (Int_t i=0; i<nSamples; i++){
    treeK0QptDSec->SetAlias(Form("%s.dMassK0.rms",samples->At(i)->GetName()),Form("1000*%s.hisK0DMassQPtTglDSec_1_1_5_1Dist.rmsG+0",samples->At(i)->GetName()));
    treeK0QptDSec->SetAlias(Form("%s.dMassK0.err",samples->At(i)->GetName()),Form("1000*%s.hisK0DMassQPtTglDSec_1_1_5_1Dist.rmsG/sqrt(%s.hisK0DMassQPtTglDSec_1_1_5_1Dist.entries+0)",samples->At(i)->GetName(),samples->At(i)->GetName()));
    treeK0QptDSec->SetAlias(Form("%s.dMassK0.mean",samples->At(i)->GetName()),Form("1000*%s.hisK0DMassQPtTglDSec_1_1_5_1Dist.meanG+0",samples->At(i)->GetName()));
    treeK0QptDSec->SetAlias(Form("%s.dMassK0.meanErr",samples->At(i)->GetName()),Form("1000*%s.hisK0DMassQPtTglDSec_1_1_5_1Dist.rmsG/sqrt(%s.hisK0DMassQPtTglDSec_1_1_5_1Dist.entries+0)",samples->At(i)->GetName(),samples->At(i)->GetName()));
    treeK0QptDSec->SetAlias(Form("%s.dMassK0.pull",samples->At(i)->GetName()),Form("%s.hisK0DMassQPtTglDSec_1_1_5_1Dist.rmsG+0",samples->At(i)->GetName()));
    treeK0QptDSec->SetAlias(Form("%s.dMassK0.pullErr",samples->At(i)->GetName()),Form("%s.hisK0DMassQPtTglDSec_1_1_5_1Dist.rmsG/sqrt(%s.hisK0DMassQPtTglDSec_1_1_5_1Dist.entries+0)",samples->At(i)->GetName(),samples->At(i)->GetName()));
  }
  for (Int_t i=0; i<nSamples; i++){
    TStatToolkit::AddMetadata(treeK0QptDSec,TString::Format("%s.dMassK0.rms.AxisTitle",samples->At(i)->GetName()).Data(),"#sigma_{M} (MeV/c)"); 
    TStatToolkit::AddMetadata(treeK0QptDSec,TString::Format("%s.dMassK0.rms.Legend",samples->At(i)->GetName()).Data(),samples->At(i)->GetName()); 
    TStatToolkit::AddMetadata(treeK0QptDSec,TString::Format("%s.dMassK0.pull.AxisTitle",samples->At(i)->GetName()).Data(),"pull  = RMS(#frac{#Delta_{M}}{#sigma^{exp}_{M}})"); 
    TStatToolkit::AddMetadata(treeK0QptDSec,TString::Format("%s.dMassK0.pull.Legend",samples->At(i)->GetName()).Data(),samples->At(i)->GetName()); 
    TStatToolkit::AddMetadata(treeK0QptDSec,TString::Format("%s.dMassK0.mean.AxisTitle",samples->At(i)->GetName()).Data(),"#delta_{M} (MeV/c)"); 
    TStatToolkit::AddMetadata(treeK0QptDSec,TString::Format("%s.dMassK0.mean.Legend",samples->At(i)->GetName()).Data(),samples->At(i)->GetName()); 
  }
  // Fig 1.1) Inv mass resol
  {
    gPad->Clear(); 
    legend = new TLegend(0.2,0.15,0.8,0.3,"K^{0} Inv.mass resolution (p_{t} 10GeV/c)");
    legend->SetBorderSize(0);legend->SetNColumns(2);
    drawString="";
    for (Int_t i=0; i<nSamples; i++) {drawString+=Form("%s.dMassK0.rms",samples->At(i)->GetName()); drawString+=(i==nSamples-1)?":":";";}
    drawString+="dsec:";
    for (Int_t i=0; i<nSamples; i++) {drawString+=Form("%s.dMassK0.err",samples->At(i)->GetName()); drawString+=(i==nSamples-1)?"":";";}
    mGraph = TStatToolkit::MakeMultGraph( treeK0QptDSec,"K0",drawString.Data(),selection,"25;21;20;24","1;2;4;6",kFALSE,1,4,legend);
    mGraph->Draw();
    legend->Draw();
    canvas->SaveAs("k0RMSMassMPtDSec.png");
    canvas->SaveAs("k0RMSMassMPtDSec.pdf");
  }
  // Fig 1.2) Inv mass bias
  {
    gPad->Clear(); 
    legend = new TLegend(0.2,0.15,0.8,0.3,"K^{0} Inv.mass bias (p_{t} 10GeV/c)");
    legend->SetBorderSize(0);legend->SetNColumns(2);
    drawString="";
    for (Int_t i=0; i<nSamples; i++) {drawString+=Form("%s.dMassK0.mean",samples->At(i)->GetName()); drawString+=(i==nSamples-1)?":":";";}
    drawString+="dsec:";
    for (Int_t i=0; i<nSamples; i++) {drawString+=Form("%s.dMassK0.meanErr",samples->At(i)->GetName()); drawString+=(i==nSamples-1)?"":";";}
    mGraph = TStatToolkit::MakeMultGraph( treeK0QptDSec,"K0",drawString.Data(),selection,"25;21;20;24","1;2;4;6",kFALSE,1,2.5,legend);
    mGraph->Draw();
    legend->Draw();
    canvas->SaveAs("k0BiasMassMPtDSec.png");
    canvas->SaveAs("k0BiasMassMPtDSec.pdf");
  }
}


void PreliminaryPlot(){
  //
  //
  //
  treeDelta0->SetMarkerStyle(21); treeDelta0->SetMarkerColor(2); 
  TLatex  latex;
  //
  //
  treeDelta0->Draw("WithoutTRD.hisDeltaP0_TRDv_qPt_tglDist.rms/WithTRD2.hisDeltaP0_TRDv_qPt_tglDist.rms:qPtCenter:abs(tglCenter)","entries>50&&abs(qPtCenter)<1","colz");
  treeDelta0->GetHistogram()->GetXaxis()->SetTitle("q/pt (1/GeV)");
  treeDelta0->GetHistogram()->GetYaxis()->SetTitle("#sigma_{r#phi}_{WithoutTRD}/#sigma_{r#phi}_{TRD2}");
  latex.DrawLatexNDC(0.2,0.15,"#sigma(r#phi_{(TPC+(TRD))}/#sigma(r#phi_{(ITS+TPC+(TRD))})");
  gPad->SaveAs("deltaP0RatioWithoutTo2.png");

  treeDelta0->Draw("WithoutTRD.hisDeltaP2_TRDv_qPt_tglDist.rms/WithTRD2.hisDeltaP2_TRDv_qPt_tglDist.rms:qPtCenter:abs(tglCenter)","entries>50&&abs(qPtCenter)<1","colz");
  treeDelta0->GetHistogram()->GetXaxis()->SetTitle("q/pt (1/GeV)");
  treeDelta0->GetHistogram()->GetYaxis()->SetTitle("#sigma_{#phi}_{WithoutTRD}/#sigma_{#phi}_{TRD2}");
  latex.DrawLatexNDC(0.2,0.15,"#sigma(#phi_{(TPC+(TRD))}/#sigma(#phi_{(ITS+TPC+(TRD))})");
  gPad->SaveAs("deltaP2RatioWithoutTo2.png");
  //
  treeDelta0->Draw("WithoutTRD.hisDeltaP4_TRDv_qPt_tglDist.rms/WithTRD2.hisDeltaP4_TRDv_qPt_tglDist.rms:qPtCenter:abs(tglCenter)","entries>50&&abs(qPtCenter)<1","colz");
  treeDelta0->GetHistogram()->GetXaxis()->SetTitle("q/pt (1/GeV)");
  treeDelta0->GetHistogram()->GetYaxis()->SetTitle("#sigma_{#phi}_{WithoutTRD}/#sigma_{#phi}_{TRD2}");
  latex.DrawLatexNDC(0.2,0.15,"#sigma(qpt_{(TPC+(TRD))}/#sigma(qpt_{(ITS+TPC+(TRD))})");
  gPad->SaveAs("deltaP4RatioWithoutTo2.png");

  //
  //
  treeDelta0->Draw("WithoutTRD.hisDeltaP0_TRDv_qPt_tglDist.rms/WithTRD2.hisDeltaP0_TRDv_qPt_tglDist.rms:qPtCenter:abs(tglCenter)","entries>50&&abs(qPtCenter)<1","colz");
  treeDelta0->GetHistogram()->GetXaxis()->SetTitle("q/pt (1/GeV)");
  treeDelta0->GetHistogram()->GetYaxis()->SetTitle("#sigma_{r#phi}_{WithoutTRD}/#sigma_{r#phi}_{TRD2}");
  latex.DrawLatexNDC(0.2,0.15,"#sigma(r#phi_{(TPC+(TRD))}/#sigma(r#phi_{(ITS+TPC+(TRD))})");
  gPad->SaveAs("deltaP0RatioWithoutTo2.png");



  treeDelta0->Draw("WithoutTRD.hisDeltaP0CQPtTglTRDDist.rmsG/WithTRD2.hisDeltaP0CQPtTglTRDDist.rmsG:qPtCenter:abs(tglCenter)","entries>200&&abs(qPtCenter)<1&&abs(rmsG/rms-1.2)<0.3","colz");
  treeDelta0->GetHistogram()->GetXaxis()->SetTitle("sin(#phi)");
  treeDelta0->GetHistogram()->GetYaxis()->SetTitle("#sigma_{r#phi}_{WithoutTRD}/#sigma_{r#phi}_{TRD4}");
  latex.DrawLatexNDC(0.2,0.15,"#sigma(r#phi_{(TPC+(TRD))}-r#phi_{(ITS+TPC+(TRD))})");
  gPad->SaveAs("deltaP0RatioWithoutTo4.png");
  //
  treeDelta0->Draw("WithoutTRD.hisPullP0CQPtTglTRDDist.rmsG/WithTRD2.hisPullP0CQPtTglTRDDist.rmsG:qPtCenter:abs(tglCenter)","entries>200&&abs(qPtCenter)<1&&abs(rmsG/rms-1.2)<0.3","colz");
  treeDelta0->GetHistogram()->GetXaxis()->SetTitle("q/pt (GeV)");
  treeDelta0->GetHistogram()->GetYaxis()->SetTitle("pull_{q/pt}_{WithoutTRD}/pull_{q/pt}_{TRD4}");
  latex.DrawLatexNDC(0.2,0.15,"pull(r#phi_{(TPC+(TRD))}-r#phi_{(ITS+TPC+(TRD))})");
  gPad->SaveAs("pullP0RatioWithoutTo4.png");

  //
  treeDelta0->Draw("WithoutTRD.hisDeltaP4CQPtTglTRDDist.rmsG/WithTRD2.hisDeltaP4CQPtTglTRDDist.rmsG:qPtCenter:abs(tglCenter)","entries>200&&abs(qPtCenter)<1&&abs(rmsG/rms-1.2)<0.3","colz");
  treeDelta0->GetHistogram()->GetXaxis()->SetTitle("q/pt (GeV)");
  treeDelta0->GetHistogram()->GetYaxis()->SetTitle("#sigma_{q/pt}_{WithoutTRD}/#sigma_{q/pt}_{TRD4}");
  latex.DrawLatexNDC(0.2,0.15,"#sigma(q/pt_{(TPC+(TRD))}-q/pt_{(ITS+TPC+(TRD))})");
  gPad->SaveAs("deltaP4RatioWithoutTo4.png");
  //
  treeDelta0->Draw("WithoutTRD.hisPullP4CQPtTglTRDDist.rmsG/WithTRD2.hisPullP4CQPtTglTRDDist.rmsG:qPtCenter:abs(tglCenter)","entries>200&&abs(qPtCenter)<1&&abs(rmsG/rms-1.2)<0.3","colz");
  treeDelta0->GetHistogram()->GetXaxis()->SetTitle("q/pt (GeV)");
  treeDelta0->GetHistogram()->GetYaxis()->SetTitle("pull_{q/pt}_{WithoutTRD}/pull_{q/pt}_{TRD4}");
  latex.DrawLatexNDC(0.2,0.15,"pull(q/pt_{(TPC+(TRD))}-q/pt_{(ITS+TPC+(TRD))})");
  gPad->SaveAs("pullP4RatioWithoutTo4.png");
  //
  //
  treeDelta0->Draw("WithTRD1.hisDeltaP2_TRDv_qPt_tglDist.rmsG/WithTRD2.hisDeltaP2_TRDv_qPt_tglDist.rmsG:qPtCenter:abs(tglCenter)","entries>200&&abs(qPtCenter)<1&&abs(rmsG/rms-1.2)<0.3","colz");
  treeDelta0->GetHistogram()->GetXaxis()->SetTitle("sin(#phi)");
  treeDelta0->GetHistogram()->GetYaxis()->SetTitle("#sigma_{#phi}_{WithTRD1}/#sigma_{#phi}_{TRD4}");
  latex.DrawLatexNDC(0.2,0.15,"#sigma(#phi_{(TPC+(TRD))}-#phi_{(ITS+TPC+(TRD))})");
  gPad->SaveAs("deltaP2RatioWithTRD1To4.png");
  //
  treeDelta0->Draw("WithTRD1.hisPullP2CQPtTglTRDDist.rmsG/WithTRD2.hisPullP2CQPtTglTRDDist.rmsG:qPtCenter:abs(tglCenter)","entries>200&&abs(qPtCenter)<1&&abs(rmsG/rms-1.2)<0.3","colz");
  treeDelta0->GetHistogram()->GetXaxis()->SetTitle("q/pt (GeV)");
  treeDelta0->GetHistogram()->GetYaxis()->SetTitle("pull_{q/pt}_{TRD1}/pull_{q/pt}_{TRD4}");
  latex.DrawLatexNDC(0.2,0.15,"pull(#phi_{(TPC+(TRD))}-#phi_{(ITS+TPC+(TRD))})");
  gPad->SaveAs("pullP2RatioWithTRD1To4.png");
  //
}


void DrawResolcomparison(){
  TLatex  latex;
  //
  TCanvas *canvas =  new TCanvas("canvas","canvas",1200, 800);
  treeCovar->SetMarkerStyle(21); treeCovar->SetMarkerSize(0.5);
  //
  treeCovar->Draw("WithTRD2.hisCovarP4ConstCQPtTglTRDDist.mean:qPtCenter:tgl","abs(qPtCenter)<0.3&tgl>0","colz");
  treeCovar->GetHistogram()->GetXaxis()->SetTitle("q/pt (GeV)");
  treeCovar->GetHistogram()->GetYaxis()->SetTitle("exp. #sigma_{q/pt}_{TRD2}");
  latex.DrawLatexNDC(0.2,0.15,"ITS+TPC+TRD setting4");
  gPad->SaveAs("covarP4WithTRD2.png");
  //
  treeCovar->Draw("WithoutTRD.hisCovarP4ConstCQPtTglTRDDist.mean/WithTRD2.hisCovarP4ConstCQPtTglTRDDist.mean:qPtCenter:tgl","entries>10&&abs(qPtCenter)<2&tgl>0","colz");
  treeCovar->GetHistogram()->GetXaxis()->SetTitle("q/pt (GeV)");
  treeCovar->GetHistogram()->GetYaxis()->SetTitle("exp. #sigma_{q/pt}_{WithoutTRD} #sigma_{q/pt}_{TRD2}");
  latex.DrawLatexNDC(0.2,0.15,"ITS+TPC+TRD setting4");
  gPad->SaveAs("covarRatioP4WithoutToTRD2.png");
  

}  


void errorComp(){

  TObjArray *hisArray= new TObjArray(10);   
  treeDelta0->SetMarkerStyle(22); treeDelta0->SetMarkerColor(2);
  treeDelta0->Draw("WithTRD2.hisDeltaP4CQPtTglTRDDist.rmsG/WithTRD1.hisDeltaP4CQPtTglTRDDist.rmsG-1:qpt>>hisRationRMS4(20,-0.5,0.5)","entries>50","prof");
  treeDelta0->GetHistogram()->GetXaxis()->SetTitle("q/p_{t} (GeV/c)");
  treeDelta0->GetHistogram()->SetMaximum(0.1);
  treeDelta0->GetHistogram()->SetMinimum(-0.2);

  hisArray->AddLast(treeDelta0->GetHistogram());
  treeDelta0->SetMarkerStyle(21); treeDelta0->SetMarkerColor(1);
  treeDelta0->Draw("WithTRD2.hisDeltaP0CQPtTglTRDDist.rmsG/WithTRD1.hisDeltaP0CQPtTglTRDDist.rmsG-1:qpt>>hisRationRMS0(20,-0.5,0.5)","entries>50","profsame");
  hisArray->AddLast(treeDelta0->GetHistogram());
  treeDelta0->SetMarkerStyle(25); treeDelta0->SetMarkerColor(4);
  treeDelta0->Draw("WithTRD2.hisDeltaP2_TRDv_qPt_tglDist.rmsG/WithTRD1.hisDeltaP2_TRDv_qPt_tglDist.rmsG-1:qpt>>hisRationRMS2(20,-0.5,0.5)","entries>50","profsame");
  hisArray->AddLast(treeDelta0->GetHistogram());


}



void xxx(){
}
