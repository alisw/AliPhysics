/*
  .L $AliPhysics_SRC/PWGPP/TPC/macros/TPCQAWebpage/tpcCalibValidation.C+
  AliDrawStyle::SetDefaults();
  AliDrawStyle::ApplyStyle("figTemplate");
  // InitTPCMCValidation("LHC15o","pass1",3,1);
  // InitTPCMCValidation("LHC17*","cpass1_pass1",3,1);
  RegisterDefaultCalibFitters();
  MakeGainFitsMinuit();

 */
#include <TError.h>
#include "TCanvas.h"
#include "TLatex.h"
#include "TTree.h"
#include "TObjArray.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TTreeFormula.h"
#include "AliExternalInfo.h"
#include "AliTreeTrending.h"
#include "TStatToolkit.h"
#include "TPRegexp.h"
#include "TGaxis.h"
#include "TStyle.h"
#include "AliDrawStyle.h"
#include "AliPainter.h"
#include "AliTreePlayer.h"
#include "AliTMinuitToolkit.h"
#include "TPaletteAxis.h"

AliExternalInfo *externalInfo = 0;
AliTreeTrending *trendingDraw = 0;
TTree *treeCalib=0;
TString period, pass;
AliTMinuitToolkit * fitter1D,* fitter2D, *fitter3D, *fitter4D;
TFormula *formula1D, *formula2D, *formula3D, *formula4D;
TF1 likeGausCachy("likeGausCachy", AliTMinuitToolkit::GaussCachyLogLike,-10,10,2);
TF1 likeAbs("likeAbs", "abs(x)",-10,10);

void RegisterDefaultCalibFitters();
Bool_t InitTPCMCValidation(TString period, TString pass, Int_t verbose,  Int_t doCheck);
void AddMetadata();

TVectorF rocGainIROC(36), rocGainOROCMedium(36), rocGainOROCLong(36);


/// Initialization of the calibration validation
/// \param period
/// \param pass
/// \param verbose
/// \param doCheck
/// \return
Bool_t InitTPCMCValidation(TString pPeriod, TString pPass, Int_t verbose,  Int_t doCheck) {
  period = pPeriod;
  pass = pPass;
  externalInfo = new AliExternalInfo(".", "", verbose);
  trendingDraw = new AliTreeTrending("mcAnchor", "mcAnchor");
  trendingDraw->SetDefaultStyle();
  gStyle->SetOptTitle(0);
  if (period.Contains("*")) {
    treeCalib = externalInfo->GetChain("QA.rawTPC", period, pass, TString("QA.TPC;QA.EVS;Logbook;Logbook.detector:TPC:detector==\"TPC\""));
  } else {
    treeCalib = externalInfo->GetTree("QA.rawTPC", period, pass, TString("QA.TPC;QA.EVS;Logbook;Logbook.detector:TPC:detector==\"TPC\""));
  }
  treeCalib->SetMarkerStyle(21);
  treeCalib->SetMarkerSize(0.4);
  RegisterDefaultCalibFitters();
  AddMetadata();
}

///

/// Cache selected production trees
/// List of production from MonaLISA
/// \param select      - selection mask
/// \param reject      - rejection mask
/// \param sourceList  - list of detectors to cache
/*!
   Example usage:
        MakeProductionCache(TPRegexp("LHC17.*"),TPRegexp("cpass0"),"QA.TPC;QA.EVS;QA.TRD;QA.rawTPC;QA.ITS;Logbook;QA.TOF;Logbook.detector");
*/
void CacheProduction(TPRegexp select, TPRegexp reject, TString sourceList){
  AliExternalInfo info;
  TTree* treeProd = info.GetTreeProdCycle();
  Int_t entries=treeProd->GetEntries();
  TObjArray * detectorArray=sourceList.Tokenize(";");
  for (Int_t i=0; i<entries; i++){
    treeProd->GetEntry(i);
    char * productionTag= (char*)treeProd->GetLeaf("Tag")->GetValuePointer();
    if (select.Match(productionTag)==0) continue;
    if (reject.Match(productionTag)==1) continue;
    printf("Caching\t%s\n",productionTag);
    TString production(productionTag);
    Int_t pos=production.First('_');
    if (pos<0) continue;
    if (pos>production.Length()-4) continue;
    printf("Caching\t%s\n",productionTag);
    TString period( production(0,pos));
    TString pass(production(pos+1, production.Length()-pos-1));
    printf("Caching\t%s\t%s\t%s\n",productionTag,period.Data(),pass.Data());
    for (Int_t iDet=0;iDet<detectorArray->GetEntries(); iDet++) {
      info.Cache(detectorArray->At(iDet)->GetName(), period.Data(), pass.Data());
    }
  }
}


void AddMetadata(){
  //
  treeCalib->SetAlias("isSelected","(QA.TPC.run==run&&QA.EVS.run==run)");
  treeCalib->SetAlias("Logbook.Bz","0.5*L3_magnetCurrent/30000.");
  treeCalib->SetAlias("rate","QA.EVS.interactionRate/1000000.");
  treeCalib->SetAlias("co2", "gasValues.fElements[2]/100");
  treeCalib->SetAlias("n2", "gasValues.fElements[3]/100");

  treeCalib->SetAlias("ptrel", "1+ptrel0");                     /// TODO  check definition in the code
  treeCalib->SetAlias("CGainMIP","gainMIP*(1+(run>278800)*0.3626)");
  TStatToolkit::AddMetadata(treeCalib,"rate.AxisTitle","rate (MHz)");
  TStatToolkit::AddMetadata(treeCalib,"rate.Title","IR");
  TStatToolkit::AddMetadata(treeCalib,"Logbook.Bz.AxisTitle","Bz (T)");
  TStatToolkit::AddMetadata(treeCalib,"Logbook.Bz.Title","B_{z}");
  //
  Int_t selected = treeCalib->Draw("co2:CGainMIP:ptrel:vdriftITS:n2",     "co2>0&&run==QA.TPC.run","goffpara");
  treeCalib->SetAlias("co2_RobustMean",TString::Format("(%f+0)",TMath::Median(selected,treeCalib->GetV1())).Data());
  treeCalib->SetAlias("CGainMIP_RobustMean",TString::Format("(%f+0)",TMath::Median(selected,treeCalib->GetV2())).Data());
  treeCalib->SetAlias("ptrel_RobustMean",TString::Format("(%f+0)",TMath::Median(selected,treeCalib->GetV3())).Data());
  treeCalib->SetAlias("vdriftITS_RobustMean",TString::Format("(%f+0)",TMath::Median(selected,treeCalib->GetV4())).Data());
  treeCalib->SetAlias("n2_RobustMean",TString::Format("(%f+0)",TMath::Median(selected,treeCalib->GetVal(4))).Data());
  treeCalib->SetAlias("dco2","(co2-co2_RobustMean)");
  treeCalib->SetAlias("dn2","(n2-n2_RobustMean)");
  treeCalib->SetAlias("co2_Norm","(co2/co2_RobustMean)");
  treeCalib->SetAlias("n2_Norm","(n2/n2_RobustMean)");
  treeCalib->SetAlias("CGainMIP_Norm","(CGainMIP/CGainMIP_RobustMean)");
  treeCalib->SetAlias("ptrel_Norm","(ptrel/ptrel_RobustMean)");
  TStatToolkit::AddMetadata(treeCalib,"CGainMIP.AxisTitle","Gain (a.u.)");
  TStatToolkit::AddMetadata(treeCalib,"CGainMIP.Legend","Gain");
  TStatToolkit::AddMetadata(treeCalib,"ptrel0.Title","#Delta_{P/T}");
  TStatToolkit::AddMetadata(treeCalib,"dco2.Title","(f_{CO2}-0.12)");
  TStatToolkit::AddMetadata(treeCalib,"dco2.AxisTitle","#delta_f_{CO2})");
  TStatToolkit::AddMetadata(treeCalib,"co2.AxisTitle","f_{CO2}");
  //
  TStatToolkit::AddMetadata(treeCalib,"co2_Norm.AxisTitle","#frac{f_{CO2}}{<f_{CO2}>}");
  TStatToolkit::AddMetadata(treeCalib,"co2_Norm.Title","#frac{f_{CO2}}{<f_{CO2}>}");
  TStatToolkit::AddMetadata(treeCalib,"co2_Norm.Legend","#frac{f_{CO2}}{<f_{CO2}>}");
  TStatToolkit::AddMetadata(treeCalib,"CGainMIP_Norm.AxisTitle","#frac{G}{<G>}");
  TStatToolkit::AddMetadata(treeCalib,"n2_Norm.AxisTitle","#frac{f_{N2}}{<f_{N2}>}");
  TStatToolkit::AddMetadata(treeCalib,"n2_Norm.Title","#frac{f_{N2}}{<f_{N2}>}");
  TStatToolkit::AddMetadata(treeCalib,"n2_Norm.Legend","#frac{f_{N2}}{<f_{N2}>}");
  //
  TStatToolkit::AddMetadata(treeCalib,"CGainMIP_Norm.Title","#frac{G}{<G>}");
  TStatToolkit::AddMetadata(treeCalib,"CGainMIP_Norm.Legend","#frac{G}{<G>}");
  TStatToolkit::AddMetadata(treeCalib,"ptrel_Norm.AxisTitle","#frac{P/T}{<P/T>}");
  TStatToolkit::AddMetadata(treeCalib,"ptrel_Norm.Legend","#frac{P/T}{<P/T>}");
  TStatToolkit::AddMetadata(treeCalib,"ptrel_Norm.Title","#frac{{P/T}{<P/T>}");
  //
  //
  treeCalib->Draw("rocGainIROC.fElements:Iteration$>>rocGainIROC(36,0,36)","isSelected&&rocGainIROC.fElements>0","profgoff");
  for (Int_t i=0; i<36; i++) rocGainIROC[i]=treeCalib->GetHistogram()->GetBinContent(i+1);
  treeCalib->Draw("rocGainOROCMedium.fElements:Iteration$>>rocGainOROCMedium(36,0,36)","isSelected&&rocGainIROC.fElements>0","profgoff");
  for (Int_t i=0; i<36; i++) rocGainOROCMedium[i]=treeCalib->GetHistogram()->GetBinContent(i+1);
  treeCalib->Draw("rocGainOROCLong.fElements:Iteration$>>rocGainOROCMedium(36,0,36)","isSelected&&rocGainIROC.fElements>0","profgoff");
  for (Int_t i=0; i<36; i++) rocGainOROCLong[i]=treeCalib->GetHistogram()->GetBinContent(i+1);

}

/// register default fitters
void RegisterDefaultCalibFitters(){
  AliTMinuitToolkit::RegisterDefaultFitters();
  fitter1D  = new AliTMinuitToolkit("AliTMinuitToolkitTest1D.root");
  fitter2D  = new AliTMinuitToolkit("AliTMinuitToolkitTest2D.root");
  fitter3D  = new AliTMinuitToolkit("AliTMinuitToolkitTest3D.root");
  fitter4D  = new AliTMinuitToolkit("AliTMinuitToolkitTest4D.root");
  formula1D = new TFormula("formula1D","[0]*x[0]");
  formula2D = new TFormula("formula2D","[0]*x[0]+[1]*x[1]");
  formula3D = new TFormula("formula3D","[0]*x[0]+[1]*x[1]+[2]*x[2]");
  formula4D = new TFormula("formula4D","[0]*x[0]+[1]*x[1]+[2]*x[2]+[3]*x[3]");
  fitter1D->SetVerbose(0x1); fitter2D->SetVerbose(0x1); fitter3D->SetVerbose(0x1);
  fitter1D->SetFitFunction((TF1*)formula1D,kTRUE);fitter2D->SetFitFunction((TF1*)formula2D,kTRUE);
  fitter3D->SetFitFunction((TF1*)formula3D,kTRUE);fitter4D->SetFitFunction((TF1*)formula4D,kTRUE);
  // set likelihood function
  likeGausCachy.SetParameters(0.8,1);
  fitter1D->SetLogLikelihoodFunction(&likeGausCachy); fitter2D->SetLogLikelihoodFunction(&likeGausCachy);
  fitter3D->SetLogLikelihoodFunction(&likeGausCachy); fitter4D->SetLogLikelihoodFunction(&likeGausCachy);
  fitter1D->SetInitialParam(new TMatrixD(1,4)); fitter2D->SetInitialParam(new TMatrixD(2,4));
  fitter3D->SetInitialParam(new TMatrixD(3,4)); fitter4D->SetInitialParam(new TMatrixD(4,4));
}


void MakeGainFitsMinuit(){
  treeCalib->SetAlias("gainFitCut","QA.EVS.run==run&&Logbook.run==run&&QA.TPC.run==run&&abs(dco2)<0.05&&abs(dn2)<0.05&&abs(attachMIP)<0.1");
  // make fits and set fit description metadata
  Int_t entries = fitter4D->FillFitter(treeCalib,"CGainMIP:0.03","1:dco2:rate:ptrel0", "gainFitCut", 0,10000000);
  TMatrixD& initParam=(*fitter4D->GetInitialParam());
  initParam(2,0)=0.1; initParam(2,1)=1;
  initParam(3,0)=0.1; initParam(3,1)=1;
  fitter4D->Fit();
   treeCalib->SetAlias("CGainMIP_Fit4",fitter4D->GetFitFunctionAsAlias().Data());
  TStatToolkit::AddMetadata(treeCalib,"CGainMIP_Fit4.AxisTitle","Fit_{LL} (a.u.)");
  TStatToolkit::AddMetadata(treeCalib,"CGainMIP_Fit4.Legend",TString::Format("Fit_{LL} %s",fitter4D->GetFitFunctionAsAlias("latex",treeCalib).Data()).Data());
  fitter4D->Bootstrap(30,"bootstrapChi2Norm");
  treeCalib->SetAlias("CGainMIP_Bootstrap4",fitter4D->GetFitFunctionAsAlias().Data());
  treeCalib->SetAlias("CGainMIP_FitRatio","CGainMIP/CGainMIP_Bootstrap4");
  TStatToolkit::AddMetadata(treeCalib,"CGainMIP_Bootstrap4.AxisTitle","Fit_{bs30} (a.u.)");
  TStatToolkit::AddMetadata(treeCalib,"CGainMIP_Bootstrap4.Legend",TString::Format("Fit_{bs30} %s",fitter4D->GetFitFunctionAsAlias("latex",treeCalib).Data()).Data());
   TStatToolkit::AddMetadata(treeCalib,"CGainMIP_FitRatio.AxisTitle","#frac{G}{G_{fit}}");
  TStatToolkit::AddMetadata(treeCalib,"CGainMIP_FitRatio.Legend","#frac{G}{G_{fit}}");

  // Draw results
  TCanvas * canvasGainFit = new TCanvas("canvasGaiFit","canvasGaiFit", 1600,1000);
  AliPainter::DivideTPad(canvasGainFit,"[1,1,1,2]",0);
  canvasGainFit->cd(1);
  TLegend* legend = new TLegend(0.11,0.65,0.5,0.89, "TPC gain");
  legend->SetMargin(0.03); legend->SetBorderSize(0);
  legend->SetEntrySeparation(0.2);
  TMultiGraph *graph = TStatToolkit::MakeMultGraph(treeCalib, "",  "CGainMIP;CGainMIP_Fit4;CGainMIP_Bootstrap4:time:0.003;0.;0","isSelected", "25;21;21","1;2;4", 0,0.6,6,legend);
  AliPainter::SetMultiGraphTimeAxis(graph,"X");
  TStatToolkit::DrawMultiGraph(graph,"ap");
  legend->Draw();
  //
  canvasGainFit->cd(2);
  legend = new TLegend(0.11,0.60,0.3,0.89, "TPC normalized gain");
  legend->SetMargin(0.2); legend->SetBorderSize(0); legend->SetNColumns(4); legend->SetEntrySeparation(0.5);
  graph = TStatToolkit::MakeMultGraph(treeCalib, "",  "CGainMIP_Norm;CGainMIP_FitRatio;co2_Norm;ptrel_Norm:time:0.001;0.001;0.001;0.001","QA.EVS.run==run", "25;25;21;21","1;2;3;4", 0,0.6,3,legend);
  AliPainter::SetMultiGraphTimeAxis(graph,"X");
  TStatToolkit::DrawMultiGraph(graph,"ap");
  legend->Draw();
  //
  canvasGainFit->cd(3);
  legend = new TLegend(0.11,0.60,0.3,0.89, "TPC normalized gain");
  legend->SetMargin(0.1); legend->SetBorderSize(0); //legend->SetNColumns(3);
  graph = TStatToolkit::MakeMultGraph(treeCalib, "",  "ptrel0;attachMIP:time:0.001;0.001","QA.EVS.run==run", "25;21;21","1;2;4", 0,0.75,10,legend);
  AliPainter::SetMultiGraphTimeAxis(graph,"X");
  TStatToolkit::DrawMultiGraph(graph,"ap");
  legend->Draw();
  //
  canvasGainFit->cd(4);
  treeCalib->Draw("CGainMIP:CGainMIP_Bootstrap4:co2","gainFitCut","colz");
  canvasGainFit->cd(4)->Update();  // needed in order to create palette
  TStatToolkit::AdaptHistoMetadata(treeCalib,0,"colz");
  canvasGainFit->cd(5);
  treeCalib->Draw("rate:co2:CGainMIP_Bootstrap4","gainFitCut","colz");
  canvasGainFit->cd(5)->Update(); // needed in order to create palette
  TStatToolkit::AdaptHistoMetadata(treeCalib,0,"colz");
  canvasGainFit->SaveAs(period+"_"+pass+"GainFit.png");
}

///
///
void makeGaindEdxReport(){
  //
  treeCalib->SetAlias("runCut","QA.TPC.run==run&&QA.EVS.run==run");
  TCanvas *canvasTime=new TCanvas("canvasTime","canvasTime",1400,500);
  TCanvas *canvasDep=new TCanvas("canvasDep","canvasDep",800,500);
  TLegend * legend=0;AddMetadata();
  TMultiGraph *graph=0;
  // 1.) Gain/CO2/N2 report
  // 1.a) time evolution
  canvasTime->cd(0);
  legend = new TLegend(0.11,0.60,0.3,0.89, "TPC normalized gain");legend->SetNColumns(3); legend->SetEntrySeparation(0.5); legend->SetBorderSize(0);
  graph = TStatToolkit::MakeMultGraph(treeCalib, "",  "CGainMIP_Norm;co2_Norm;n2_Norm:time","QA.EVS.run==run&&CGainMIP>0", "25;21;21","1;2;4", 0,0.6,6,legend);
  AliPainter::SetMultiGraphTimeAxis(graph,"X");
  TStatToolkit::DrawMultiGraph(graph,"ap");legend->Draw();
  gPad->SaveAs("trendingGain_CO2_N2.png");
  // 1.b) dependence (for different N2 concentration -different slope (2017)
  canvasDep->cd(0);
  treeCalib->Draw("CGainMIP_Norm:co2_Norm:n2_Norm","QA.TPC.run==run&&CGainMIP_Norm>0.5&&abs(n2_Norm-1)<0.03","colz");
  TStatToolkit::AdaptHistoMetadata(treeCalib,0,"colz");
  gPad->SaveAs("gain_VsCO2_N2.png");
  // 2.) Vdrift/CO2/N2
  // 2.a) time evolution
  canvasTime->cd(0);
  legend = new TLegend(0.11,0.60,0.3,0.89, "TPC vdrift correction");legend->SetNColumns(3); legend->SetEntrySeparation(0.5); legend->SetBorderSize(0);
  graph = TStatToolkit::MakeMultGraph(treeCalib, "",  "1+vdriftITS;co2_Norm;n2_Norm:time","QA.EVS.run==run&&CGainMIP>0.5", "25;21;21","1;2;4", 0,0.6,6,legend);
  AliPainter::SetMultiGraphTimeAxis(graph,"X");
  TStatToolkit::DrawMultiGraph(graph,"ap"); legend->Draw();
  gPad->SaveAs("trendingDriftCorr_CO2_N2.png");
  // 2.b) dependence (for different N2 concentration -different slope (2017)
  canvasDep->cd(0);
  treeCalib->Draw("vdriftITS:co2_Norm:n2_Norm","QA.TPC.run==run&&CGainMIP_Norm>0.5&&abs(n2_Norm-1)<0.03","colz");
  TStatToolkit::AdaptHistoMetadata(treeCalib,0,"colz");
  gPad->SaveAs("vdriftCorr_VsCO2_N2.png");
  //
  canvasTime->cd(0);
  legend = new TLegend(0.11,0.60,0.3,0.89, "TPC PID");legend->SetNColumns(1); legend->SetEntrySeparation(0.5); legend->SetBorderSize(0);
  graph = TStatToolkit::MakeMultGraph(treeCalib, "",  "(meanMIPele/meanMIP)/1.6;(resolutionMIP/resolutionMIPele)/1.2;co2_Norm;CGainMIP_Norm:time","runCut&&CGainMIP>0&&interactionRate<300000", "25;25;21;21","1;2;4;6", 0,0.6,3,legend);
  AliPainter::SetMultiGraphTimeAxis(graph,"X");
  TStatToolkit::DrawMultiGraph(graph,"ap");legend->Draw();
  gPad->SaveAs("trendingPIDSeparation.png");
  //
  canvasDep->cd(0);
  treeCalib->Draw("resolutionMIP/resolutionMIPele:CGainMIP_Norm:interactionRate","runCut&&meanMIP>40&&CGainMIP_Norm>0.5&&interactionRate<300000&&resolutionMIP/resolutionMIPele>1","colz");
  TStatToolkit::AdaptHistoMetadata(treeCalib,0,"colz");
  gPad->SaveAs("resolMIPToEl_vsGain_Rate.png");

  canvasDep->cd(0);
  treeCalib->Draw("meanMIPele/meanMIP:CGainMIP_Norm:interactionRate","runCut&&meanMIP>40&&CGainMIP_Norm>0.5&&interactionRate<300000&&resolutionMIP/resolutionMIPele>1","colz");
  TStatToolkit::AdaptHistoMetadata(treeCalib,0,"colz");
  gPad->SaveAs("meanEltoMIP_vsGain_Rate.png");

  canvasDep->cd(0);
  treeCalib->Draw("resolutionMIP/resolutionMIPele:CGainMIP_Norm:1+attachMIP","runCut&&meanMIP>40&&CGainMIP_Norm>0.5&&interactionRate<300000&&resolutionMIP/resolutionMIPele>1&&abs(attachMIP)<0.02","colz");
  TStatToolkit::AdaptHistoMetadata(treeCalib,0,"colz");
  gPad->SaveAs("resolMIPToEl_vsGain_Attach.png");


}

/// Logbook data volume studies
void FitDataVolume(TTree * treeLogbook) {
  // treeLogbook = externalInfo->GetChain("Logbook", "LHC17*", "cpass1_pass1", TString("QA.TPC;QA.EVS;Logbook;Logbook.detector:TPC:detector==\"TPC\";Logbook"));
  // AliTMinuitToolkit::RegisterDefaultFitters();
  treeLogbook->SetMarkerSize(0.5);treeLogbook->SetMarkerStyle(25);
  gStyle->SetStatX(0.88); gStyle->SetStatY(0.88);
  Int_t entries= treeLogbook->Draw("(bytesInjectedPhysics/eventCountPhysics):QA.EVS.interactionRate:run","Logbook.detector_TPC.run==run&&QA.EVS.run==run&&QA.EVS.interactionRate>10","colz");
  TGraphErrors *gr = new TGraphErrors(entries, treeLogbook->GetV2(), treeLogbook->GetV1());
  AliTMinuitToolkit::Fit(gr,"pol1R","bootstrap50", "","funOption(2,2,2)");  // fit robust linear - bootstrap20
  TStatToolkit::AdaptHistoMetadata(treeLogbook,0,"colz");
  gr->Draw("p");
  //
  treeLogbook->SetAlias("sizeFitRatio","(bytesInjectedPhysics/eventCountPhysics)/(400000+QA.EVS.interactionRate*40.69)"); // TODO get from Fit
  treeLogbook->Draw("sizeFitRatio:time","Logbook.detector_TPC.run==run&&QA.EVS.run==run&&QA.EVS.interactionRate>10","",10000);
  treeLogbook->GetHistogram()->GetXaxis()->SetTimeFormat("%d/%m");
  TMultiGraph * graph = TStatToolkit::MakeMultGraph(treeLogbook, "",  "sizeFitRatio:time","Logbook.detector_TPC.run==run&&QA.EVS.run==run&&QA.EVS.interactionRate>10", "25;25;21;21","1;2;4;6", 0,0.6,3,0);
  AliPainter::SetMultiGraphTimeAxis(graph,"X");
  TStatToolkit::DrawMultiGraph(graph,"ap");

}
