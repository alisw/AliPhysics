/*
  .L $AliPhysics_SRC/PWGPP/TPC/macros/TPCQAWebpage/tpcCalibValidation.C+
  AliDrawStyle::SetDefaults();
  AliDrawStyle::ApplyStyle("figTemplate");
  InitTPCMCValidation("LHC15o","pass1",3,1);

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
  period=pPeriod;
  pass=pPass;
  externalInfo = new AliExternalInfo(".", "", verbose);
  trendingDraw = new AliTreeTrending("mcAnchor","mcAnchor");
  trendingDraw->SetDefaultStyle();
  gStyle->SetOptTitle(0);
  treeCalib = externalInfo->GetTree("QA.rawTPC",period,pass,"QA.TPC;QA.EVS;Logbook");
  treeCalib->SetMarkerStyle(21);
  treeCalib->SetMarkerSize(0.4);
  RegisterDefaultCalibFitters();
  AddMetadata();
}

void AddMetadata(){
  //
  treeCalib->SetAlias("isSelected","QA.TPC.run==run&&QA.EVS.run==run");
  treeCalib->SetAlias("Logbook.Bz","0.5*L3_magnetCurrent/30000.");
  treeCalib->SetAlias("rate","QA.EVS.interactionRate/1000000.");
  treeCalib->SetAlias("co2", "gasValues.fElements[2]/100");
  treeCalib->SetAlias("ptrel", "1+ptrel0");   /// TODO  check definition in the code
  TStatToolkit::AddMetadata(treeCalib,"rate.AxisTitle","rate (MHz)");
  TStatToolkit::AddMetadata(treeCalib,"rate.Title","IR");
  TStatToolkit::AddMetadata(treeCalib,"Logbook.Bz.AxisTitle","Bz (T)");
  TStatToolkit::AddMetadata(treeCalib,"Logbook.Bz.Title","B_{z}");
  //
  TStatToolkit::SetStatusAlias(treeCalib, "co2",     "isSelected&&co2>0", Form("varname_RobustMean:(MeanEF+0):%f",0.98));
  TStatToolkit::SetStatusAlias(treeCalib, "gainMIP", "isSelected&&co2>0", Form("varname_RobustMean:(MeanEF+0):%f",0.98));
  TStatToolkit::SetStatusAlias(treeCalib, "ptrel", "isSelected&&co2>0", Form("varname_RobustMean:(MeanEF+0):%f",0.98));
  treeCalib->SetAlias("dco2","(co2-co2_RobustMean)");
  treeCalib->SetAlias("co2_Norm","(co2/co2_RobustMean)");
  treeCalib->SetAlias("gainMIP_Norm","(gainMIP/gainMIP_RobustMean)");
  treeCalib->SetAlias("ptrel_Norm","(ptrel/ptrel_RobustMean)");
  TStatToolkit::AddMetadata(treeCalib,"gainMIP.AxisTitle","Gain (a.u.)");
  TStatToolkit::AddMetadata(treeCalib,"gainMIP.Legend","Gain");
  TStatToolkit::AddMetadata(treeCalib,"ptrel0.Title","#Delta_{P/T}");
  TStatToolkit::AddMetadata(treeCalib,"dco2.Title","(f_{CO2}-0.12)");
  TStatToolkit::AddMetadata(treeCalib,"dco2.AxisTitle","#delta_f_{CO2})");
  TStatToolkit::AddMetadata(treeCalib,"co2.AxisTitle","f_{CO2}");
  //
  TStatToolkit::AddMetadata(treeCalib,"co2_Norm.AxisTitle","#frac{f_{CO2}}{<f_{CO2}>}");
  TStatToolkit::AddMetadata(treeCalib,"co2_Norm.Title","#frac{f_{CO2}}{<f_{CO2}>}");
  TStatToolkit::AddMetadata(treeCalib,"co2_Norm.Legend","#frac{f_{CO2}}{<f_{CO2}>}");
  TStatToolkit::AddMetadata(treeCalib,"gainMIP_Norm.AxisTitle","#frac{G}{<G>}");
  TStatToolkit::AddMetadata(treeCalib,"gainMIP_Norm.Title","#frac{G}{<G>}");
  TStatToolkit::AddMetadata(treeCalib,"gainMIP_Norm.Legend","#frac{G}{<G>}");
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
  treeCalib->SetAlias("gainFitCut","QA.EVS.run==run&&Logbook.run==run&&abs(dco2)<0.05");
  // make fits and set fit description metadata
  Int_t entries = fitter4D->FillFitter(treeCalib,"gainMIP:0.03","1:dco2:rate:ptrel0", "gainFitCut", 0,10000000);
  TMatrixD& initParam=(*fitter4D->GetInitialParam());
  initParam(0,0)=2; initParam(0,1)=1;
  initParam(1,0)=0.1; initParam(1,1)=1;
  initParam(2,0)=0.1; initParam(2,1)=1;
  initParam(3,0)=0.1; initParam(3,1)=1;
  fitter4D->Fit();
   treeCalib->SetAlias("gainMIP_Fit4",fitter4D->GetFitFunctionAsAlias().Data());
  TStatToolkit::AddMetadata(treeCalib,"gainMIP_Fit4.AxisTitle","Fit_{LL} (a.u.)");
  TStatToolkit::AddMetadata(treeCalib,"gainMIP_Fit4.Legend",TString::Format("Fit_{LL} %s",fitter4D->GetFitFunctionAsAlias("latex",treeCalib).Data()).Data());
  fitter4D->Bootstrap(30,"bootstrapChi2Norm");
  treeCalib->SetAlias("gainMIP_Bootstrap4",fitter4D->GetFitFunctionAsAlias().Data());
  treeCalib->SetAlias("gainMIP_FitRatio","gainMIP/gainMIP_Bootstrap4");
  TStatToolkit::AddMetadata(treeCalib,"gainMIP_Bootstrap4.AxisTitle","Fit_{bs30} (a.u.)");
  TStatToolkit::AddMetadata(treeCalib,"gainMIP_Bootstrap4.Legend",TString::Format("Fit_{bs30} %s",fitter4D->GetFitFunctionAsAlias("latex",treeCalib).Data()).Data());
   TStatToolkit::AddMetadata(treeCalib,"gainMIP_FitRatio.AxisTitle","#frac{G}{G_{fit}}");
  TStatToolkit::AddMetadata(treeCalib,"gainMIP_FitRatio.Legend","#frac{G}{G_{fit}}");

  // Draw results
  TCanvas * canvasGainFit = new TCanvas("canvasGaiFit","canvasGaiFit", 1600,1000);
  AliDrawStyle::DivideTPad(canvasGainFit,"1,1,1,2");
  canvasGainFit->cd(1);
  TLegend* legend = new TLegend(0.11,0.65,0.5,0.89, "TPC gain");
  legend->SetMargin(0.03); legend->SetBorderSize(0);
  legend->SetEntrySeparation(0.2);
  TMultiGraph *graph = TStatToolkit::MakeMultGraph(treeCalib, "",  "gainMIP;gainMIP_Fit4;gainMIP_Bootstrap4:time:0.003;0.;0","isSelected", "25;21;21","1;2;4", 0,0.6,3,legend);
  AliDrawStyle::SetMultiGraphTimeAxis(graph,"X");
  TStatToolkit::DrawMultiGraph(graph,"ap");
  legend->Draw();
  //
  canvasGainFit->cd(2);
  legend = new TLegend(0.11,0.60,0.3,0.89, "TPC normalized gain");
  legend->SetMargin(0.2); legend->SetBorderSize(0); legend->SetNColumns(4); legend->SetEntrySeparation(0.5);
  graph = TStatToolkit::MakeMultGraph(treeCalib, "",  "gainMIP_Norm;gainMIP_FitRatio;co2_Norm;ptrel_Norm:time:0.001;0.001;0.001;0.001","QA.EVS.run==run", "25;25;21;21","1;2;3;4", 0,0.6,3,legend);
  AliDrawStyle::SetMultiGraphTimeAxis(graph,"X");
  TStatToolkit::DrawMultiGraph(graph,"ap");
  legend->Draw();
  //
  canvasGainFit->cd(3);
  legend = new TLegend(0.11,0.60,0.3,0.89, "TPC normalized gain");
  legend->SetMargin(0.1); legend->SetBorderSize(0); //legend->SetNColumns(3);
  graph = TStatToolkit::MakeMultGraph(treeCalib, "",  "ptrel0;attachMIP:time:0.001;0.001","QA.EVS.run==run", "25;21;21","1;2;4", 0,0.75,10,legend);
  AliDrawStyle::SetMultiGraphTimeAxis(graph,"X");
  TStatToolkit::DrawMultiGraph(graph,"ap");
  legend->Draw();
  //
  canvasGainFit->cd(4);
  treeCalib->Draw("gainMIP:gainMIP_Bootstrap4:co2","gainFitCut","colz");
  canvasGainFit->cd(4)->Update();  // needed in order to create palette
  TStatToolkit::AdaptHistoMetadata(treeCalib,0,"colz");
  canvasGainFit->cd(5);
  treeCalib->Draw("rate:co2:gainMIP_Bootstrap4","gainFitCut","colz");
  canvasGainFit->cd(5)->Update(); // needed in order to create palette
  TStatToolkit::AdaptHistoMetadata(treeCalib,0,"colz");
}


