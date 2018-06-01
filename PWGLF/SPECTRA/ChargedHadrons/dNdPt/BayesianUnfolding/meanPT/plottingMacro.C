#include "TStyle.h"
#include "TFile.h"
#include "TList.h"
#include "TH1D.h"
#include "TH2D.h"
#include "THn.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TGaxis.h"
#include "TGraph2D.h"
#include "TObjArray.h"
#include "TGFrame.h"
#include "TGFileBrowser.h"
#include "TBrowser.h"



#include <string>
#include <iostream>

#include "Plotting.h"

using namespace std;
TList* makeIndividualPlots(string fileName, TList* histos);
TList* makeComparisons(TList* histos);

string outputFilename("OutputCanvae.root");

Bool_t includeSimulation = kFALSE;
Bool_t includeSystematics = kFALSE;
Bool_t includeClosureTests = kTRUE;


void plottingMacro(){
  TGaxis::SetMaxDigits(3);
  //==================================================================================================
  // Input:
  //==================================================================================================

  // Get Input histograms and store them in List (with respective sublists)
  TList* histos = new TList();
  histos->SetOwner();

  if(includeSystematics){
    histos->Add(getHistosFromFile("pp_5TeV_08_Sys"));
    histos->Add(getHistosFromFile("pp_7TeV_08_Sys"));
    histos->Add(getHistosFromFile("pp_8TeV_08_Sys"));
    histos->Add(getHistosFromFile("pp_13TeV_08_Sys"));
  }

  histos->Add(getHistosFromFile("pp_5TeV_08"));
  histos->Add(getHistosFromFile("pp_5TeV_08_pass3"));
  histos->Add(getHistosFromFile("pp_5TeV_08_pass3_Pythia6"));
  histos->Add(getHistosFromFile("pp_5TeV_08_pass3_secondMC"));
  histos->Add(getHistosFromFile("pp_7TeV_08"));
  histos->Add(getHistosFromFile("pp_8TeV_08"));
  histos->Add(getHistosFromFile("pp_13TeV_08"));
  histos->Add(getHistosFromFile("pp_13TeV_08_EPOS"));
  histos->Add(getHistosFromFile("pp_13TeV_08_secondMC"));
  histos->Add(getHistosFromFile("pPb_5TeV_08"));
  histos->Add(getHistosFromFile("PbPb_5TeV_08_rebin"));
  histos->Add(getHistosFromFile("XeXe_5TeV_08_rebin"));
  histos->Add(getHistosFromFile("pp_5TeV_08_iterationCheck"));

  // obtain plots from publication
  histos->Add(getHistosFromMacro());


//==================================================================================================
// Plotting:
//==================================================================================================

  TObjArray* plots = new TObjArray();
  plots->SetOwner();

  // individual plots
  plots->Add(makeIndividualPlots("pp_5TeV_08", histos));

//  plots->Add(makeIndividualPlots("pp_5TeV_08", histos));
//  plots->Add(makeIndividualPlots("pp_5TeV_08_pass3", histos));
//  plots->Add(makeIndividualPlots("pp_5TeV_08_pass3_Pythia6", histos));
//  plots->Add(makeIndividualPlots("pp_5TeV_08_pass3_secondMC", histos));
//  plots->Add(makeIndividualPlots("pp_13TeV_08", histos));
//  plots->Add(makeIndividualPlots("pp_13TeV_08_EPOS", histos));
//  plots->Add(makeIndividualPlots("pp_13TeV_08_secondMC", histos));
//  plots->Add(makeIndividualPlots("pp_8TeV_08", histos));


//  ((TCanvas*)((TList*)plots->At(0))->At(0))->ls();//Print("basdf.pdf");//tempCanvas->Print("mario.pdf");

/*
  plots->Add(makeIndividualPlots("pp_7TeV_08", histos));
  plots->Add(makeIndividualPlots("pp_8TeV_08", histos));

  plots->Add(makeIndividualPlots("pPb_5TeV_08", histos));
  plots->Add(makeIndividualPlots("PbPb_5TeV_08_rebin", histos));
  plots->Add(makeIndividualPlots("XeXe_5TeV_08_rebin", histos));
*/
  // combined plots
  plots->Add(makeComparisons(histos));


//==================================================================================================
// Output:
//==================================================================================================

  TFile* outputFile =  new TFile(outputFilename.c_str(),"RECREATE");
  outputFile->cd();
  histos->Write("histos", TObject::kSingleKey);

  // Write Plots in file
  for(Int_t i = 0; i < plots->GetEntries(); i++){

    TList* currList = (TList*)plots->At(i);
    string currListName = currList->GetName();

    outputFile->mkdir(currListName.c_str());
    outputFile->cd(currListName.c_str());
    currList->Write();
  }

  TBrowser* outputInspector = new TBrowser();
}











//==================================================================================================
//==================================================================================================
// Functions
//==================================================================================================
//==================================================================================================

TList* makeIndividualPlots(string fileName, TList* histos){

  TList* list = new TList();
  list->SetOwner();
  list->SetName((fileName + "_plots").c_str());

  string alice = "#bf{ALICE work in progress}";
  string mc = "#bf{ALICE simulation}";
  string charged = "Charged Particles";

  string colsys = "pp";
  string ptReach = "0.15 GeV/#it{c} < #it{p}_{T} < 10 GeV/#it{c}";
  string erg = "#sqrt{#it{s}} = 5.02 TeV";
  string eta = "|#it{#eta}| < 0.8";
  string nEvents;
  string nEventsMC;
  string nIter = "";//""after 10 iterations";

  Int_t viewingRangeFull = 100;

  Int_t viewingRangeMoment1Measured = 50;
  Int_t viewingRangeMoment2Measured = 40;
  Int_t viewingRangeMoment3Measured = 30;

  Int_t viewingRangeMoment1 = 60;//60
  Int_t viewingRangeMoment2 = 60; //60
  Int_t viewingRangeMoment3 = 60;//50

  Double_t yRangeMoments[] = {0.41, 5.3};
  Double_t yRangeMomentsMeasured[] = {0.41, 5.3};

  TString qualifier = fileName.c_str();

  if(qualifier.Contains("5TeV")) erg = "#sqrt{#it{s}} = 5.02 TeV";
  if(qualifier.Contains("7TeV")) erg = "#sqrt{#it{s}} = 7 TeV";
  if(qualifier.Contains("8TeV")) erg = "#sqrt{#it{s}} = 8 TeV";
  if(qualifier.Contains("13TeV")) erg = "#sqrt{#it{s}} = 13 TeV";
  if(qualifier.Contains("XeXe_5TeV")) erg = "#sqrt{#it{s}} = 5.44 TeV";

  if(qualifier.Contains("pp")) {colsys = "pp";}
  if(qualifier.Contains("pPb")) {
    colsys = "p-Pb"; erg.insert(12,"_{NN}");
    viewingRangeMoment1 = 200; viewingRangeMoment2 = 200; viewingRangeMoment3 = 200;
    viewingRangeMoment1Measured = 200; viewingRangeMoment2Measured = 200; viewingRangeMoment3Measured = 200;
  }
  if(qualifier.Contains("PbPb")) {
    colsys = "Pb-Pb"; erg.insert(12,"_{NN}");
    viewingRangeMoment1 = 4500; viewingRangeMoment2 = 4500; viewingRangeMoment3 = 4500;
    viewingRangeMoment1Measured = 4500; viewingRangeMoment2Measured = 4500; viewingRangeMoment3Measured = 4500;
  }
  if(qualifier.Contains("XeXe")) {
    colsys = "Xe-Xe"; erg = "#sqrt{#it{s_{NN}}} = 5.44 TeV";
    viewingRangeMoment1 = 4500; viewingRangeMoment2 = 4500; viewingRangeMoment3 = 4500;
    viewingRangeMoment1Measured = 4500; viewingRangeMoment2Measured = 4500; viewingRangeMoment3Measured = 4500;
  }

  if(qualifier.Contains("08")) eta = "|#it{#eta}| < 0.8";
  if(qualifier.Contains("03")) eta = "|#it{#eta}| < 0.3";


  nEvents = getIntegral("multDistMeasured", fileName, histos);
  nEventsMC = getIntegral("multDistMeasuredMC", fileName, histos);


  string tempText;
  TObjArray* tempArray = new TObjArray(1);
  tempArray->SetOwner(kFALSE);
  TObjArray* tempArrayRatios = new TObjArray(1);
  tempArrayRatios->SetOwner(kFALSE);

  string tempTitles[15];
  TCanvas* tempCanvas = NULL;
  string tempCanvasName;

  Short_t defaultColors[14]={4,2,kGreen+3,kMagenta+2,8,9,11,12,13,14,15,16,17,18};
  Short_t defaultMarker[14]={20,21,34,33,24,25,26,27,28,29,30,2,3,5};
  Short_t tempColors[14]={kGray,kGreen+2,28,7,8,9,11,kRed,kBlue,13,kBlack,16,17,12};
  Short_t tempMarker[14]={20,21,34,33,24,25,26,27,28,29,30,2,3,5};

  Short_t colors[14]={kGray,kGreen+2,28,7,8,9,11,kRed,kBlue,13,kBlack,16,17,12};
  Short_t markers[14]={20,21,22,23,24,25,26,27,28,29,30,2,3,5};


  //-------------------- Unfolding Ratio  ------------------------------------------------------
  tempCanvasName = "ratioUnfoldedNotUnfolded";
  TH2D* tempStandard = (TH2D*)getClone("multPtMeasuredMC", fileName, histos);
  TH2D* tempUnfolded = (TH2D*)getClone("multPtMeasuredPtResUnfoldedMC", fileName, histos);
  TH2D* tempGenerated = (TH2D*)getClone("multPtGeneratedMC", fileName, histos);

  TH1D* withUnf = tempUnfolded->ProjectionY();
  TH1D* woUnf = tempStandard->ProjectionY();
  TH1D* gen = tempGenerated->ProjectionY();

  TH1D* ratioWith = withUnf->Clone();
  ratioWith->Divide(gen);
  TH1D* ratioWo = woUnf->Clone();
  ratioWo->Divide(gen);

  TH1D* ratioWwo = withUnf->Clone();
  ratioWwo->Divide(woUnf);

  tempArray->Add(withUnf);
  tempArray->Add(woUnf);
//  tempArray->Add(gen);

  tempText = mc + "<|>" + colsys + ", " + erg + ", " +  eta + "<|>" + ptReach;
//  tempArray->Add(makeText(0.4, 0.3, tempText, 3));

  tempArrayRatios->Add(ratioWwo);
//  tempArrayRatios->Add(ratioWith);
//  tempArrayRatios->Add(ratioWo);
  tempArrayRatios->Add(new TF1("line", "1", 0.15, 10));

  tempTitles[0] = "unfolded";
  tempTitles[1] = "not unfolded";
//  tempTitles[2] = "truth";
  tempArray->Add(makeLegend(0.15, 0.75, tempArray, tempTitles));


  tempCanvas = makeCanvas(appendFileName(tempCanvasName, fileName), tempArray, tempArrayRatios, "logX logY square");

  list->Add((tempCanvas->Clone()));
  tempCanvas->Close();
  tempArray->Clear();

  //-------------------- ptRes  ------------------------------------------------------
  tempCanvasName = "ptResolution";
  tempArray->Add(getClone("ptResolution", fileName, histos));

  tempCanvas = makeCanvas(appendFileName(tempCanvasName, fileName), tempArray, 0, "logZ square");

  list->Add((tempCanvas->Clone()));
  tempCanvas->Close();
  tempArray->Clear();

  tempCanvasName = "ptResponseMatrix";
  tempArray->Add(getClone("ptResponseMatrix", fileName, histos));

  tempCanvas = makeCanvas(appendFileName(tempCanvasName, fileName), tempArray, 0, "logZ square");

  list->Add((tempCanvas->Clone()));
  tempCanvas->Close();
  tempArray->Clear();

/*
  //-------------------- Efficiency  ------------------------------------------------------
  tempCanvasName = "efficiencyPt";
  tempArray->Add(getClone("efficiencyPt", fileName, histos));
  setTitle("Y", "#it{#epsilon}^{MC}_{prim}", (TH1D*)tempArray->At(0));

  tempText = mc + "<|>" + colsys + ", " + erg + ", " +  eta + "<|>" + ptReach;
  tempArray->Add(makeText(0.4, 0.3, tempText, 3));

  tempCanvas = makeCanvas(appendFileName(tempCanvasName, fileName), tempArray, 0, "logX square");

  list->Add((tempCanvas->Clone()));
  tempCanvas->Close();
  tempArray->Clear();

  //-------------------- Secondary contamination  ------------------------------------------------------
  tempCanvasName = "secContPt";
  tempArray->Add(getClone("secContPt", fileName, histos));
  setTitle("Y", "#it{#xi}^{MC}_{sec}", (TH1D*)tempArray->At(0));

  tempText = mc + "<|>" + colsys + ", " + erg + ", " +  eta + "<|>" + ptReach;
  tempArray->Add(makeText(0.4, 0.3, tempText, 3));

  tempCanvas = makeCanvas(appendFileName(tempCanvasName, fileName), tempArray, 0, "logX square");

  list->Add((tempCanvas->Clone()));
  tempCanvas->Close();
  tempArray->Clear();
*/
  //-------------------- Response Matrix ------------------------------------------------------
  tempCanvasName = "responseMatrix";
  tempArray->Add(getClone("responseMatrix", fileName, histos));
  setTitle("Z", "#it{P}(#it{N}_{acc}|#it{N}_{ch})", (TH1D*)tempArray->At(0));
  setRange("Z", 1e-3, 1, (TH1D*)tempArray->At(0));

  tempText = mc + "<|>" + colsys + ", " + erg + ", " +  eta + "<|>" + ptReach;
  tempArray->Add(makeText(0.4, 0.3, tempText, 3));

  tempCanvas = makeCanvas(appendFileName(tempCanvasName, fileName), tempArray, 0, "logz square");

  list->Add((tempCanvas->Clone()));
  tempCanvas->Close();
  tempArray->Clear();

  //-------------------- Response Matrix Tracks -----------------------------------------------
  tempCanvasName = "responseMatrixTracks";
  tempArray->Add(getClone("responseMatrixTracks", fileName, histos));
  setTitle("Z", "#it{P}_{part}(#it{N}_{acc}|#it{N}_{ch})", (TH1D*)tempArray->At(0));
  setRange("Z", 1e-3, 1, (TH1D*)tempArray->At(0));

  tempText = mc + "<|>" + colsys + ", " + erg + ", " +  eta + "<|>" + ptReach;
  tempArray->Add(makeText(0.4, 0.3, tempText, 3));

  tempCanvas = makeCanvas(appendFileName(tempCanvasName, fileName), tempArray, 0, "logz square");

  list->Add((tempCanvas->Clone()));
  tempCanvas->Close();
  tempArray->Clear();
  //-------------------- RespMat tracks projY ------------------------------------------------------
  tempCanvasName = "initialDist";
  tempArray->Add(((TH2D*)getClone("responseMatrixTracksOrig", fileName, histos))->ProjectionY());

  setRange("X", 0, viewingRangeMoment1, (TH1D*)tempArray->At(0));

  tempText = alice + "<|>" + colsys + ", " + erg + ", " +  eta + "<|>" + ptReach;
  tempArray->Add(makeText(0.5, 0.9, tempText, 3));

  tempCanvas = makeCanvas(appendFileName(tempCanvasName, fileName), tempArray, 0, "thick square");

  list->Add((tempCanvas->Clone()));
  tempCanvas->Close();
  tempArray->Clear();

  //-------------------- Unfolding Matrix ------------------------------------------------------
  tempCanvasName = "unfoldingMatrix";
  tempArray->Add(getClone("unfoldingMatrix", fileName, histos));
  setTitle("Z", "#it{P}(#it{N}_{ch}|#it{N}_{acc})", (TH1D*)tempArray->At(0));
  setRange("Z", 1e-3, 1, (TH1D*)tempArray->At(0));

  tempText = alice + "<|>" + colsys + ", " + erg + ", " +  eta + "<|>" + ptReach;
  tempArray->Add(makeText(0.4, 0.3, tempText, 3));

  tempCanvas = makeCanvas(appendFileName(tempCanvasName, fileName), tempArray, 0, "logz square");

  list->Add((tempCanvas->Clone()));
  tempCanvas->Close();
  tempArray->Clear();

  //-------------------- Unfolding Matrix ------------------------------------------------------
  tempCanvasName = "unfoldingMatrixMC";
  tempArray->Add(getClone("unfoldingMatrixMC", fileName, histos));
  setTitle("Z", "#it{P}(#it{N}_{ch}|#it{N}_{acc})", (TH1D*)tempArray->At(0));
  setRange("Z", 1e-3, 1, (TH1D*)tempArray->At(0));

  tempText = mc + "<|>" + colsys + ", " + erg + ", " +  eta + "<|>" + ptReach;
  tempArray->Add(makeText(0.4, 0.3, tempText, 3));

  tempCanvas = makeCanvas(appendFileName(tempCanvasName, fileName), tempArray, 0, "logz square");

  list->Add((tempCanvas->Clone()));
  tempCanvas->Close();
  tempArray->Clear();

  //----------------- Multiplicity Distributions ---------------------------------------------
  tempCanvasName = "multDistributions";
  tempArray->Add(getClone("multDistMeasured", fileName, histos));
  tempArray->Add(getClone("multDistUnfolded", fileName, histos));
  tempArray->Add(getClone("multDistGeneratedMC", fileName, histos));

  setTitle("X", "Multiplicity", (TH1D*)tempArray->At(0));
  setTitle("Y", "# Events", (TH1D*)tempArray->At(0));

  tempTitles[0] = "#it{n}_{evt}(#it{N}_{acc})";
  tempTitles[1] = "#hat{#it{n}}_{evt}(#it{N}_{ch}) " + nIter;
  tempTitles[2] = "#it{n}_{evt}(#it{N}_{ch})_{MC}";

  tempLeg = makeLegend(0.44, 0.92, tempArray, tempTitles);
  tempArray->Add(tempLeg);

  tempText = alice + "<|>" + colsys + ", " + erg + ", " +  eta + "<|>" + ptReach;
  tempArray->Add(makeText(0.16, 0.32, tempText));

  tempCanvas = makeCanvas(appendFileName(tempCanvasName, fileName), tempArray, 0, "logY square thick");

  list->Add((tempCanvas->Clone()));
  tempCanvas->Close();
  tempArray->Clear();

  //-------------------- Pt Spectrum Measured ------------------------------------------------------
  tempCanvasName = "multPtMeasured";
  tempArray->Add(getClone("multPtMeasured", fileName, histos));
  setTitle("Z", "1/#it{N}_{evt} 1/(2#pi #it{p}_{T}) (d^{3}#it{N})/(d#it{p}_{T}d#it{#eta}d#it{N}_{acc}) (GeV/#it{c})^{-2}", (TH1D*)tempArray->At(0));

  setRange("X", 0, 90, (TH1D*)tempArray->At(0));
  setRange("Y", 0.15, 10, (TH1D*)tempArray->At(0));
  setRange("Z", 1e-10, 5, (TH1D*)tempArray->At(0));

  tempText = alice + "<|>" + colsys + ", " + erg + ", " +  eta + "<|>" + ptReach;
  tempArray->Add(makeText(0.15, 0.5, tempText));

  tempCanvas = makeCanvas(appendFileName(tempCanvasName, fileName), tempArray, 0, "logz square 2D");

  list->Add((tempCanvas->Clone()));
  tempCanvas->Close();
  tempArray->Clear();

  //-------------------- Pt Spectrum Measured MC --------------------------------------------------
  tempCanvasName = "multPtMeasuredMC";
  tempArray->Add(getClone("multPtMeasuredMC", fileName, histos));
  setTitle("Z", "1/#it{N}_{evt} 1/(2#pi #it{p}_{T}) (d^{3}#it{N})/(d#it{p}_{T}d#it{#eta}d#it{N}_{acc}) (GeV/#it{c})^{-2}", (TH1D*)tempArray->At(0));
  setTitle("Y", "#it{p}_{T} (GeV/#it{c})", (TH1D*)tempArray->At(0));
  setRange("X", 0, 90, (TH1D*)tempArray->At(0));
  setRange("Y", 0.15, 10, (TH1D*)tempArray->At(0));
  setRange("Z", 1e-10, 5, (TH1D*)tempArray->At(0));

  tempText = mc + "<|>" + colsys + ", " + erg + ", " +  eta + "<|>" + ptReach;
  tempArray->Add(makeText(0.14, 0.28, tempText));

  tempCanvas = makeCanvas(appendFileName(tempCanvasName, fileName), tempArray, 0, "logz square 2D");

  list->Add((tempCanvas->Clone()));
  tempCanvas->Close();
  tempArray->Clear();

  //-------------------- Pt Spectrum Unfolded ------------------------------------------------------
  tempCanvasName = "multPtUnfolded";
  tempArray->Add(getClone("multPtUnfolded", fileName, histos));
  setTitle("Z", "1/#it{N}_{evt} 1/(2#pi #it{p}_{T}) (d^{3}#it{N})/(d#it{p}_{T}d#it{#eta}d#it{N}_{ch}) (GeV/#it{c})^{-2}", (TH1D*)tempArray->At(0));
  setRange("X", 0, 90, (TH1D*)tempArray->At(0));
  setRange("Y", 0.15, 10, (TH1D*)tempArray->At(0));
  setRange("Z", 1e-10, 5, (TH1D*)tempArray->At(0));

  tempText = alice + "<|>" + colsys + ", " + erg + ", " +  eta + "<|>" + ptReach;
  tempArray->Add(makeText(0.15, 0.5, tempText));

  tempCanvas = makeCanvas(appendFileName(tempCanvasName, fileName), tempArray, 0, "logz square 2D");

  list->Add((tempCanvas->Clone()));
  tempCanvas->Close();
  tempArray->Clear();

  //-------------------- Pt Spectrum Unfolded MC -----------------------------------------------
  tempCanvasName = "multPtUnfoldedMC";
  tempArray->Add(getClone("multPtUnfoldedMC", fileName, histos));
  setTitle("Z", "1/#it{N}_{evt} 1/(2#pi #it{p}_{T}) (d^{3}#it{N})/(d#it{p}_{T}d#it{#eta}d#it{N}_{ch}) (GeV/#it{c})^{-2}", (TH1D*)tempArray->At(0));
  setTitle("Y", "#it{p}_{T} (GeV/#it{c})", (TH1D*)tempArray->At(0));
  setRange("X", 0, 90, (TH1D*)tempArray->At(0));
  setRange("Y", 0.15, 10, (TH1D*)tempArray->At(0));
  setRange("Z", 1e-10, 5, (TH1D*)tempArray->At(0));

  tempText = mc + "<|>" + colsys + ", " + erg + ", " +  eta + "<|>" + ptReach;
  tempArray->Add(makeText(0.14, 0.28, tempText));

  tempCanvas = makeCanvas(appendFileName(tempCanvasName, fileName), tempArray, 0, "logz square 2D");

  list->Add((tempCanvas->Clone()));
  tempCanvas->Close();
  tempArray->Clear();

  //-------------------- Moments Measured ----------------------------------------------------
  tempCanvasName = "momentsMeasured";
  tempArray->Add(getClone("momentMeasured1", fileName, histos));
  tempArray->Add(getClone("momentMeasured2", fileName, histos));
  tempArray->Add(getClone("momentMeasured3", fileName, histos));

  setTitle("Y", "#LT#it{p}^{n}_{T}#GT (GeV/#it{c})^{n}", (TH1D*)tempArray->At(0));
  setRange("X", 0, viewingRangeFull, (TH1D*)tempArray->At(0));
  setRange("Y", yRangeMomentsMeasured[0], yRangeMomentsMeasured[1], (TH1D*)tempArray->At(0));

  cutHist((TH1D*)tempArray->At(0), viewingRangeMoment1Measured);
  cutHist((TH1D*)tempArray->At(1), viewingRangeMoment2Measured);
  cutHist((TH1D*)tempArray->At(2), viewingRangeMoment3Measured);

  tempTitles[0] = "n = 1";
  tempTitles[1] = "n = 2";
  tempTitles[2] = "n = 3";

  tempLeg = makeLegend(0.15, 0.7, tempArray, tempTitles);
  tempArray->Add(tempLeg);

  tempText = alice + "<|>" + colsys + ", " + erg + ", " +  eta + "<|>" + ptReach + "<|><|>Measured Moments";
  tempArray->Add(makeText(0.15, 0.9, tempText, 5));

  tempCanvas = makeCanvas(appendFileName(tempCanvasName, fileName), tempArray, 0, "thick square");

  list->Add((tempCanvas->Clone()));
  tempCanvas->Close();
  tempArray->Clear();

  //-------------------- Moments Reweighted ----------------------------------------------------
  tempCanvasName = "momentsReweighted";
  tempArray->Add(getClone("momentReweighted1", fileName, histos));
  tempArray->Add(getClone("momentReweighted2", fileName, histos));
  tempArray->Add(getClone("momentReweighted3", fileName, histos));
  setRange("X", 0, viewingRangeFull, (TH1D*)tempArray->At(0));
  setRange("Y", yRangeMoments[0], yRangeMoments[1], (TH1D*)tempArray->At(0));
  setTitle("Y", "#LT#it{p}^{n}_{T}#GT (GeV/#it{c})^{n}", (TH1D*)tempArray->At(0));

  cutHist((TH1D*)tempArray->At(0), viewingRangeMoment1);
  cutHist((TH1D*)tempArray->At(1), viewingRangeMoment2);
  cutHist((TH1D*)tempArray->At(2), viewingRangeMoment3);

  tempTitles[0] = "n = 1";
  tempTitles[1] = "n = 2";
  tempTitles[2] = "n = 3";

  tempLeg = makeLegend(0.15, 0.7, tempArray, tempTitles);
  tempArray->Add(tempLeg);

  tempText = alice + "<|>" + colsys + ", " + erg + ", " +  eta + "<|>" + ptReach + "<|><|>Re-weighted Moments";
  tempArray->Add(makeText(0.15, 0.9, tempText, 5));

  tempCanvas = makeCanvas(appendFileName(tempCanvasName, fileName), tempArray, 0, "thick square");

  list->Add((tempCanvas->Clone()));
  tempCanvas->Close();
  tempArray->Clear();

  //-------------------- Moments Unfolded ------------------------------------------------------
  tempCanvasName = "momentsUnfolded";
  tempArray->Add(getClone("momentUnfolded1", fileName, histos));
  tempArray->Add(getClone("momentUnfolded2", fileName, histos));
  tempArray->Add(getClone("momentUnfolded3", fileName, histos));
  setRange("X", 0, viewingRangeFull, (TH1D*)tempArray->At(0));
  setRange("Y", yRangeMoments[0], yRangeMoments[1], (TH1D*)tempArray->At(0));
  setTitle("Y", "#LT#it{p}^{n}_{T}#GT (GeV/#it{c})^{n}", (TH1D*)tempArray->At(0));

  cutHist((TH1D*)tempArray->At(0), viewingRangeMoment1);
  cutHist((TH1D*)tempArray->At(1), viewingRangeMoment2);
  cutHist((TH1D*)tempArray->At(2), viewingRangeMoment3);

  tempTitles[0] = "n = 1";
  tempTitles[1] = "n = 2";
  tempTitles[2] = "n = 3";

  tempLeg = makeLegend(0.15, 0.7, tempArray, tempTitles);
  tempArray->Add(tempLeg);

  tempText = alice + "<|>" + colsys + ", " + erg + ", " +  eta + "<|>" + ptReach + "<|><|>Unfolded Moments";
  tempArray->Add(makeText(0.15, 0.9, tempText, 5));

  tempCanvas = makeCanvas(appendFileName(tempCanvasName, fileName), tempArray, 0, "thick square");

  list->Add((tempCanvas->Clone()));
  tempCanvas->Close();
  tempArray->Clear();

  //-------------------- Mean + RMS Measured ------------------------------------------------------
  tempCanvasName = "meanRMSMeasured";
  tempArray->Add(getClone("momentMeasured1", fileName, histos));
  tempArray->Add(getClone("momentMeasuredRMS", fileName, histos));

  tempTitles[0] = "mean";
  tempTitles[1] = "rms";


  cutHist((TH1D*)tempArray->At(0), 50);
  cutHist((TH1D*)tempArray->At(1), 40);
  setRange("Y", 0.52, 1.65, (TH1D*)tempArray->At(0));
  setRange("X", 0, viewingRangeFull, (TH1D*)tempArray->At(0));
  setTitle("Y", "#LT#it{p}_{T}#GT, #sqrt{#LT#it{p}^{2}_{T}#GT} (GeV/#it{c})", (TH1D*)tempArray->At(0));

  tempArray->Add(makeLegend(0.15, 0.7, tempArray, tempTitles));

  tempText = alice + "<|>" + colsys + ", " + erg + ", " +  eta + "<|>" + ptReach + "<|><|>Measured";
  tempArray->Add(makeText(0.15, 0.9, tempText, 5));

  tempCanvas = makeCanvas(appendFileName(tempCanvasName, fileName), tempArray, 0, "thick square");

  list->Add((tempCanvas->Clone()));
  tempCanvas->Close();
  tempArray->Clear();
  //-------------------- Mean + RMS Reweighted ------------------------------------------------------
  tempCanvasName = "meanRMSReweighted";
  tempArray->Add(getClone("momentReweighted1", fileName, histos));
  tempArray->Add(getClone("momentReweightedRMS", fileName, histos));

  tempTitles[0] = "mean";
  tempTitles[1] = "rms";

  cutHist((TH1D*)tempArray->At(0), viewingRangeMoment1);
  cutHist((TH1D*)tempArray->At(1), viewingRangeMoment1);
  setRange("Y", 0.52, 1.65, (TH1D*)tempArray->At(0));
  setRange("X", 0, viewingRangeMoment1, (TH1D*)tempArray->At(0));
  setTitle("Y", "#LT#it{p}_{T}#GT, #sqrt{#LT#it{p}^{2}_{T}#GT} (GeV/#it{c})", (TH1D*)tempArray->At(0));

  tempArray->Add(makeLegend(0.15, 0.7, tempArray, tempTitles));

  tempText = alice + "<|>" + colsys + ", " + erg + ", " +  eta + "<|>" + ptReach + "<|><|>Re-weighted";
  tempArray->Add(makeText(0.15, 0.9, tempText, 5));

  tempCanvas = makeCanvas(appendFileName(tempCanvasName, fileName), tempArray, 0, "thick square");

  list->Add((tempCanvas->Clone()));
  tempCanvas->Close();
  tempArray->Clear();
  //-------------------- Mean + RMS Unfolded ------------------------------------------------------
  tempCanvasName = "meanRMSUnfolded";
  tempArray->Add(getClone("momentUnfolded1", fileName, histos));
  tempArray->Add(getClone("momentUnfoldedRMS", fileName, histos));

  tempTitles[0] = "mean";
  tempTitles[1] = "rms";


  cutHist((TH1D*)tempArray->At(0), viewingRangeMoment1);
  cutHist((TH1D*)tempArray->At(1), viewingRangeMoment1);
  setRange("Y", 0.48, 1.32, (TH1D*)tempArray->At(0));
  setRange("X", 0, viewingRangeMoment1, (TH1D*)tempArray->At(0));
  setTitle("Y", "#LT#it{p}_{T}#GT, #sqrt{#LT#it{p}^{2}_{T}#GT} (GeV/#it{c})", (TH1D*)tempArray->At(0));

  tempArray->Add(makeLegend(0.15, 0.7, tempArray, tempTitles));

  tempText = alice + "<|>" + colsys + ", " + erg + ", " +  eta + "<|>" + ptReach + "<|><|>Unfolded";
  tempArray->Add(makeText(0.15, 0.9, tempText, 5));

  tempCanvas = makeCanvas(appendFileName(tempCanvasName, fileName), tempArray, 0, "thick square");

  list->Add((tempCanvas->Clone()));
  tempCanvas->Close();
  tempArray->Clear();

  //-------------------- Mean + RMS Measured ------------------------------------------------------
  tempCanvasName = "meanRMSMeasuredMC";
  tempArray->Add(getClone("momentMeasuredMC1", fileName, histos));
  tempArray->Add(getClone("momentMeasuredMCRMS", fileName, histos));

  tempTitles[0] = "mean";
  tempTitles[1] = "rms";


  cutHist((TH1D*)tempArray->At(0), 50);
  cutHist((TH1D*)tempArray->At(1), 40);
  setRange("Y", 0.52, 1.65, (TH1D*)tempArray->At(0));
  setRange("X", 0, viewingRangeFull, (TH1D*)tempArray->At(0));
  setTitle("Y", "#LT#it{p}_{T}#GT, #sqrt{#LT#it{p}^{2}_{T}#GT} (GeV/#it{c})", (TH1D*)tempArray->At(0));

  tempArray->Add(makeLegend(0.15, 0.7, tempArray, tempTitles));

  tempText = alice + "<|>" + colsys + ", " + erg + ", " +  eta + "<|>" + ptReach + "<|><|>Measured";
  tempArray->Add(makeText(0.15, 0.9, tempText, 5));

  tempCanvas = makeCanvas(appendFileName(tempCanvasName, fileName), tempArray, 0, "thick square");

  list->Add((tempCanvas->Clone()));
  tempCanvas->Close();
  tempArray->Clear();
  //-------------------- Mean + RMS Reweighted ------------------------------------------------------
  tempCanvasName = "meanRMSReweightedMC";
  tempArray->Add(getClone("momentReweightedMC1", fileName, histos));
  tempArray->Add(getClone("momentReweightedMCRMS", fileName, histos));

  tempTitles[0] = "mean";
  tempTitles[1] = "rms";

  cutHist((TH1D*)tempArray->At(0), viewingRangeMoment1);
  cutHist((TH1D*)tempArray->At(1), viewingRangeMoment1);
  setRange("Y", 0.52, 1.65, (TH1D*)tempArray->At(0));
  setRange("X", 0, viewingRangeMoment1, (TH1D*)tempArray->At(0));
  setTitle("Y", "#LT#it{p}_{T}#GT, #sqrt{#LT#it{p}^{2}_{T}#GT} (GeV/#it{c})", (TH1D*)tempArray->At(0));

  tempArray->Add(makeLegend(0.15, 0.7, tempArray, tempTitles));

  tempText = alice + "<|>" + colsys + ", " + erg + ", " +  eta + "<|>" + ptReach + "<|><|>Re-weighted";
  tempArray->Add(makeText(0.15, 0.9, tempText, 5));

  tempCanvas = makeCanvas(appendFileName(tempCanvasName, fileName), tempArray, 0, "thick square");

  list->Add((tempCanvas->Clone()));
  tempCanvas->Close();
  tempArray->Clear();

  //-------------------- Mean + RMS Unfolded ------------------------------------------------------
  tempCanvasName = "meanRMSUnfoldedMC";
  tempArray->Add(getClone("momentUnfoldedMC1", fileName, histos));
  tempArray->Add(getClone("momentUnfoldedMCRMS", fileName, histos));

  tempTitles[0] = "mean";
  tempTitles[1] = "rms";


  cutHist((TH1D*)tempArray->At(0), viewingRangeMoment1);
  cutHist((TH1D*)tempArray->At(1), viewingRangeMoment1);
  setRange("Y", 0.48, 1.32, (TH1D*)tempArray->At(0));
  setRange("X", 0, viewingRangeMoment1, (TH1D*)tempArray->At(0));
  setTitle("Y", "#LT#it{p}_{T}#GT, #sqrt{#LT#it{p}^{2}_{T}#GT} (GeV/#it{c})", (TH1D*)tempArray->At(0));

  tempArray->Add(makeLegend(0.15, 0.7, tempArray, tempTitles));

  tempText = alice + "<|>" + colsys + ", " + erg + ", " +  eta + "<|>" + ptReach + "<|><|>Unfolded";
  tempArray->Add(makeText(0.15, 0.9, tempText, 5));

  tempCanvas = makeCanvas(appendFileName(tempCanvasName, fileName), tempArray, 0, "thick square");

  list->Add((tempCanvas->Clone()));
  tempCanvas->Close();
  tempArray->Clear();

  //-------------------- Mean + RMS Generated ------------------------------------------------------
  tempCanvasName = "meanRMSGeneratedMC";
  tempArray->Add(getClone("momentGeneratedMC1", fileName, histos));
  tempArray->Add(getClone("momentGeneratedMCRMS", fileName, histos));

  tempTitles[0] = "mean";
  tempTitles[1] = "rms";


  cutHist((TH1D*)tempArray->At(0), viewingRangeMoment1);
  cutHist((TH1D*)tempArray->At(1), viewingRangeMoment1);
  setRange("Y", 0.48, 1.32, (TH1D*)tempArray->At(0));
  setRange("X", 0, viewingRangeMoment1, (TH1D*)tempArray->At(0));
  setTitle("Y", "#LT#it{p}_{T}#GT, #sqrt{#LT#it{p}^{2}_{T}#GT} (GeV/#it{c})", (TH1D*)tempArray->At(0));

  tempArray->Add(makeLegend(0.15, 0.7, tempArray, tempTitles));

  tempText = alice + "<|>" + colsys + ", " + erg + ", " +  eta + "<|>" + ptReach + "<|><|>Generated";
  tempArray->Add(makeText(0.15, 0.9, tempText, 5));

  tempCanvas = makeCanvas(appendFileName(tempCanvasName, fileName), tempArray, 0, "thick square");

  list->Add((tempCanvas->Clone()));
  tempCanvas->Close();
  tempArray->Clear();

  //-------------------- Model Comparison Moment 1 ------------------------------------------------------
  if(includeSimulation){
    Short_t colorsModelComp[14]={kBlue+1,kBlack,kRed,kGreen+4,8,9,11,kRed,kBlue,13,kBlack,16,17,12};
    Short_t markersModelComp[14]={20,24,4,28,24,25,26,27,28,29,30,2,3,5};

    tempCanvasName = "modelCompMoment1";
    tempArray->Add(getClone("momentUnfolded1", fileName, histos));
    tempArray->Add(getClone("momentUnfolded1SystErr", fileName, histos));
    tempArray->Add(getClone("momentGeneratedMC1", fileName, histos));
    tempArray->Add(getClone("momentSimulatedMC1", fileName, histos));
    tempArray->Add(getClone("momentSimulatedMPI1MC1", fileName, histos));
    tempArray->Add(getClone("momentSimulatedMPI2MC1", fileName, histos));
    tempArray->Add(getClone("momentSimulatedMPI3MC1", fileName, histos));

    setRange("X", 0, 75, (TH1D*)tempArray->At(0));
    setRange("Y", 0.4, 0.95, (TH1D*)tempArray->At(0));
    cutHist((TH1D*)tempArray->At(0), 75);

    tempTitles[0] = "data";
    tempTitles[1] = "dummy";
    tempTitles[2] = "PYTHIA8";
    tempTitles[3] = "PYTHIA8 no CR";
    tempTitles[4] = "PYTHIA8 nMPI = 1";
    tempTitles[5] = "PYTHIA8 nMPI = 2";
    tempTitles[6] = "PYTHIA8 nMPI = 3";
    tempArray->Add(makeLegend(0.15, 0.75, tempArray, tempTitles));

    tempText = alice + "<|>" + colsys + ", " + erg + ", " +  eta + "<|>" + ptReach;
    tempArray->Add(makeText(0.15, 0.9, tempText,4));

    tempCanvas = makeCanvas(appendFileName(tempCanvasName, fileName), tempArray, 0, "thick square", colorsModelComp, markersModelComp);

    list->Add((tempCanvas->Clone()));
    tempCanvas->Close();
    tempArray->Clear();

    //-------------------- Model Comparison Moment 2 ------------------------------------------------------
    tempCanvasName = "modelCompMoment2";
    tempArray->Add(getClone("momentUnfolded2", fileName, histos));
    tempArray->Add(getClone("momentUnfolded2SystErr", fileName, histos));
    tempArray->Add(getClone("momentGeneratedMC2", fileName, histos));
    tempArray->Add(getClone("momentSimulatedMC2", fileName, histos));

    setRange("X", 0, 75, (TH1D*)tempArray->At(0));
    cutHist((TH1D*)tempArray->At(0), 75);
    setRange("Y", 0.25, 1.425, (TH1D*)tempArray->At(0));

    tempTitles[0] = "data";
    tempTitles[1] = "dummy";
    tempTitles[2] = "PYTHIA8";
    tempTitles[3] = "PYTHIA8 no CR";
    tempArray->Add(makeLegend(0.15, 0.75, tempArray, tempTitles));

    tempText = alice + "<|>" + colsys + ", " + erg + ", " +  eta + "<|>" + ptReach;
    tempArray->Add(makeText(0.15, 0.9, tempText,4));


    tempCanvas = makeCanvas(appendFileName(tempCanvasName, fileName), tempArray, 0, "thick square", colorsModelComp, markersModelComp);

    list->Add((tempCanvas->Clone()));
    tempCanvas->Close();
    tempArray->Clear();

    //-------------------- Model Comparison Moment 3 ------------------------------------------------------
    tempCanvasName = "modelCompMoment3";
    tempArray->Add(getClone("momentUnfolded3", fileName, histos));
    tempArray->Add(getClone("momentUnfolded3SystErr", fileName, histos));
    tempArray->Add(getClone("momentGeneratedMC3", fileName, histos));
    tempArray->Add(getClone("momentSimulatedMC3", fileName, histos));

    setRange("X", 0, 75, (TH1D*)tempArray->At(0));
    cutHist((TH1D*)tempArray->At(0), 75);

    tempTitles[0] = "data";
    tempTitles[1] = "dummy";
    tempTitles[2] = "PYTHIA8";
    tempTitles[3] = "PYTHIA8 no CR";
    tempArray->Add(makeLegend(0.15, 0.75, tempArray, tempTitles));

    tempText = alice + "<|>" + colsys + ", " + erg + ", " +  eta + "<|>" + ptReach;
    tempArray->Add(makeText(0.15, 0.9, tempText,4));

    tempCanvas = makeCanvas(appendFileName(tempCanvasName, fileName), tempArray, 0, "thick square", colorsModelComp, markersModelComp);

    list->Add((tempCanvas->Clone()));
    tempCanvas->Close();
    tempArray->Clear();
  }

  //-------------------- Systematics ------------------------------------------------------
  if(includeSystematics){
    Short_t colorsSyst[14]={kPink+8,kGreen+2,28,kOrange+2,8,kCyan-6,kMagenta+3,kRed,kBlue,13,kBlack,16,17,12};
    Short_t markersSyst[14]={20,21,22,23,24,25,26,27,28,29,30,2,3,5};

    tempCanvasName = "trackSystematicsMean";

    tempArray->Add(getClone("momentUnfolded1_DCAtoVertexZ", (fileName + "_Sys").c_str(), histos));
    tempArray->Add(getClone("momentUnfolded1_DCAtoVertexXYPtDep", (fileName + "_Sys").c_str(), histos));
    tempArray->Add(getClone("momentUnfolded1_RatioCrossedRowsOverFindableClustersTPC", (fileName + "_Sys").c_str(), histos));
    tempArray->Add(getClone("momentUnfolded1_FractionSharedClustersTPC", (fileName + "_Sys").c_str(), histos));
    tempArray->Add(getClone("momentUnfolded1_Maxchi2perTPCclu", (fileName + "_Sys").c_str(), histos));
    tempArray->Add(getClone("momentUnfolded1_Maxchi2perITSclu", (fileName + "_Sys").c_str(), histos));
    tempArray->Add(getClone("momentUnfolded1_ClusterReqITS", (fileName + "_Sys").c_str(), histos));
    tempArray->Add(getClone("momentUnfolded1_DeadzoneWidth", (fileName + "_Sys").c_str(), histos));
    tempArray->Add(getClone("momentUnfolded1_Ncrnclgeomlength", (fileName + "_Sys").c_str(), histos));
    tempArray->Add(getClone("momentUnfolded1_MaxChi2TPCConstrained", (fileName + "_Sys").c_str(), histos));
    tempArray->Add(getClone("momentUnfolded1_total", (fileName + "_Sys").c_str(), histos));

    tempTitles[0] = "  1  ";
    tempTitles[1] = "  2  ";
    tempTitles[2] = "  3  ";
    tempTitles[3] = "  4  ";
    tempTitles[4] = "  5  ";
    tempTitles[5] = "  6  ";
    tempTitles[6] = "  7  ";
    tempTitles[7] = "  8  ";
    tempTitles[8] = "  9  ";
    tempTitles[9] = "  10 ";
    tempTitles[10] = "total";


    setRange("X", 0, viewingRangeMoment1, (TH1D*)tempArray->At(0));
    setRange("Y", 0.0, 0.0215, (TH1D*)tempArray->At(0));
    setTitle("Y", "(variation_{ID} - nominal)/nominal", (TH1D*)tempArray->At(0));

    tempArray->Add(makeLegend(0.27, 0.8, tempArray, tempTitles,0.6, kTRUE));

    tempText = "Systematic uncertainty of #LT#it{p}_{T}#GT<|>for track selection criterion with ID:<|>";
    tempArray->Add(makeText(0.27, 0.9, tempText, 3));

    tempText = alice + "<|>" + colsys + ", " + erg + ", " +  eta + "<|>" + ptReach;
    tempArray->Add(makeText(0.27, 0.6, tempText, 5));


    tempCanvas = makeCanvas(appendFileName(tempCanvasName, fileName), tempArray, 0, "2xthick square", colorsSyst, markersSyst);

    list->Add((tempCanvas->Clone()));
    tempCanvas->Close();
    tempArray->Clear();
  }

  //-------------------- Closure Tests ------------------------------------------------------
  if(includeClosureTests){

    TList* closureList = new TList();
    closureList->SetOwner();
    closureList->SetName((fileName + "_closure").c_str());

    //----------------- Closure Test MultDist --------------------------------------
    tempCanvasName = "closureTestMultDist";
    tempArray->Add(getClone("multDistMeasuredClosureTest", fileName, histos));
    tempArray->Add(getClone("multDistUnfoldedClosureTest", fileName, histos));
    tempArray->Add(getClone("multDistGeneratedClosureTest", fileName, histos));

    setTitle("X", "Multiplicity", (TH1D*)tempArray->At(0));
    setTitle("Y", "# Events", (TH1D*)tempArray->At(0));
    setRange("Y", 0.2, 3e7, (TH1D*)tempArray->At(0));

    tempArrayRatios->Add(getRatio("multDistUnfoldedClosureTest", "multDistGeneratedClosureTest", "unf./gen.", fileName, histos));
    tempArrayRatios->Add(new TF1("line", "1", -0.5, 99.5));
    setRange("Y", 0.8, 1.2, (TH1D*)tempArrayRatios->At(0));
    setTitle("X", "Multiplicity", (TH1D*)tempArrayRatios->At(0));

    tempTitles[0] = "#it{n}_{evt}(#it{N}_{acc})_{MC} - reconstructed";
    tempTitles[1] = "#hat{#it{n}}_{evt}(#it{N}_{ch})_{MC}   - unfolded " + nIter;
    tempTitles[2] = "#it{n}_{evt}(#it{N}_{ch})_{MC}   - generated";

    tempArray->Add(makeLegend(0.33, 0.94, tempArray, tempTitles));

    tempText = mc + "<|>" + colsys + ", " + erg + ", " +  eta + "<|>" + ptReach;
    tempArray->Add(makeText(0.15, 0.23, tempText,4));

    tempCanvas = makeCanvas(appendFileName(tempCanvasName, fileName), tempArray, tempArrayRatios, "thick square logY");

    closureList->Add((tempCanvas->Clone()));
    tempCanvas->Close();
    tempArray->Clear();
    tempArrayRatios->Clear();
/*
    //----------------- Closure Test MultDist Flat --------------------------------------
    tempCanvasName = "closureTestMultDistFlat";
    tempArray->Add(getClone("multDistMeasuredClosureTest", fileName, histos));
    tempArray->Add(getClone("multDistUnfoldedClosureTestFlat", fileName, histos));
    tempArray->Add(getClone("multDistGeneratedClosureTest", fileName, histos));
    tempArray->Add(getClone("multDistInitialClosureTestFlat", fileName, histos));

    setTitle("X", "Multiplicity", (TH1D*)tempArray->At(0));
    setTitle("Y", "# Events", (TH1D*)tempArray->At(0));
    setRange("Y", 0.2, 3e7, (TH1D*)tempArray->At(0));

    tempArrayRatios->Add(getRatio("multDistUnfoldedClosureTestFlat", "multDistGeneratedClosureTest", "unf./gen.", fileName, histos));
    tempArrayRatios->Add(new TF1("line", "1", -0.5, 99.5));
    setRange("Y", 0.8, 1.2, (TH1D*)tempArrayRatios->At(0));
    setTitle("X", "Multiplicity", (TH1D*)tempArrayRatios->At(0));

    tempTitles[0] = "#it{n}_{evt}(#it{N}_{acc})_{MC} - reconstructed";
    tempTitles[1] = "#hat{#it{n}}_{evt}(#it{N}_{ch})_{MC}   - unfolded " + nIter;
    tempTitles[2] = "#it{n}_{evt}(#it{N}_{ch})_{MC}   - generated";
    tempTitles[3] = "#it{n}_{evt}(#it{N}_{ch})_{MC}   - initial";

    tempArray->Add(makeLegend(0.33, 0.94, tempArray, tempTitles));

    tempText = mc + "<|>" + colsys + ", " + erg + ", " +  eta + "<|>" + ptReach;
    tempArray->Add(makeText(0.15, 0.23, tempText,4));

    tempCanvas = makeCanvas(appendFileName(tempCanvasName, fileName), tempArray, tempArrayRatios, "thick square logY");

    closureList->Add((tempCanvas->Clone()));
    tempCanvas->Close();
    tempArray->Clear();
    tempArrayRatios->Clear();
*/
    //----------------- Closure Test pt MultDistTracks --------------------------------------
    tempCanvasName = "closureTestPtDistTracks";
    tempArray->Add(getClone("ptDistMeasuredMC", fileName, histos));
    tempArray->Add(getClone("ptDistGeneratedMC", fileName, histos));
    setTitle("Y", "# Tracks", (TH1D*)tempArray->At(0));

    tempArrayRatios->Add(getRatio("ptDistMeasuredMC", "ptDistGeneratedMC", "meas./gen.", fileName, histos));
    setTitle("X", "Multiplicity", (TH1D*)tempArrayRatios->At(0));
    tempArrayRatios->Add(new TF1("line", "1", -0.5, 99.5));

    tempTitles[0] = "meas";
    tempTitles[1] = "gen";

    tempArray->Add(makeLegend(0.33, 0.94, tempArray, tempTitles));

    tempText = mc + "<|>" + colsys + ", " + erg + ", " +  eta + "<|>" + ptReach;
    tempArray->Add(makeText(0.15, 0.2, tempText,4));

    tempCanvas = makeCanvas(appendFileName(tempCanvasName, fileName), tempArray, tempArrayRatios, "thick square logY");

    closureList->Add((tempCanvas->Clone()));
    tempCanvas->Close();
    tempArray->Clear();
    tempArrayRatios->Clear();

    //----------------- Closure Test MultDistTracks --------------------------------------
    tempCanvasName = "closureTestMultDistTracks";
    tempArray->Add(getClone("multDistTracksMeasuredClosureTest", fileName, histos));
    tempArray->Add(getClone("multDistTracksUnfoldedClosureTest", fileName, histos));
    tempArray->Add(getClone("multDistTracksGeneratedClosureTest", fileName, histos));
    tempArray->Add(getClone("multDistTracksInitialClosureTest", fileName, histos));

    setTitle("X", "Multiplicity", (TH1D*)tempArray->At(0));
    setTitle("Y", "# Tracks", (TH1D*)tempArray->At(0));

    tempArrayRatios->Add(getRatio("multDistTracksUnfoldedClosureTest", "multDistTracksGeneratedClosureTest", "unf./gen.", fileName, histos));
    setTitle("X", "Multiplicity", (TH1D*)tempArrayRatios->At(0));
    tempArrayRatios->Add(new TF1("line", "1", -0.5, 99.5));

    tempTitles[0] = "#it{n}_{evt}(#it{N}_{acc})_{MC} - reconstructed";
    tempTitles[1] = "#hat{#it{n}}_{evt}(#it{N}_{ch})_{MC}   - unfolded " + nIter;
    tempTitles[2] = "#it{n}_{evt}(#it{N}_{ch})_{MC}   - generated";
    tempTitles[3] = "#it{n}_{evt}(#it{N}_{ch})_{MC}   - initial";
    tempArray->Add(makeLegend(0.33, 0.94, tempArray, tempTitles));

    tempText = mc + "<|>" + colsys + ", " + erg + ", " +  eta + "<|>" + ptReach;
    tempArray->Add(makeText(0.15, 0.2, tempText,4));
    tempCanvas = makeCanvas(appendFileName(tempCanvasName, fileName), tempArray, tempArrayRatios, "thick square logY");

    closureList->Add((tempCanvas->Clone()));
    tempCanvas->Close();
    tempArray->Clear();
    tempArrayRatios->Clear();


/*
    //----------------- Closure Test MultPt --------------------------------------
    tempCanvasName = "closureTestMultPt";

    TH2D* efficiency2D = (TH2D*)getClone("multPtUnfoldedMC", fileName, histos);
    TH2D* temp2DHist = (TH2D*)getClone("multPtGeneratedMC", fileName, histos);
    efficiency2D->Divide(temp2DHist);
    tempArray->Add(efficiency2D);
    delete temp2DHist;

    setRange("X", -0.5, 72, (TH1D*)tempArray->At(0));
    setRange("Y", 0.15, 10, (TH1D*)tempArray->At(0));

    tempCanvas = makeCanvas(appendFileName(tempCanvasName, fileName), tempArray, 0, "thick square logY");

    closureList->Add((tempCanvas->Clone()));
    tempCanvas->Close();
    tempArray->Clear();
    tempArrayRatios->Clear();


    //----------------- Closure Test MultPt --------------------------------------
    tempCanvasName = "closureTestMultPt_ProjMult";

    tempArray->Add(efficiency2D->ProjectionX());
    tempCanvas = makeCanvas(appendFileName(tempCanvasName, fileName), tempArray, 0, "thick square");

    closureList->Add((tempCanvas->Clone()));
    tempCanvas->Close();
    tempArray->Clear();

    //----------------- Closure Test MultPt --------------------------------------
    tempCanvasName = "closureTestMultPt_ProjPt";

    tempArray->Add(efficiency2D->ProjectionY());
    tempCanvas = makeCanvas(appendFileName(tempCanvasName, fileName), tempArray, 0, "thick square logX");

    closureList->Add((tempCanvas->Clone()));
    tempCanvas->Close();
    tempArray->Clear();


    //----------------- Closure Test MultPt staterr --------------------------------------
    tempCanvasName = "closureTestMultPtErrors";

    TH2D* relativeError = getRelativeErrorRatio(efficiency2D);
    tempArray->Add(relativeError);

    setRange("X", -0.5, 72, (TH1D*)tempArray->At(0));
    setRange("Y", 0.15, 10, (TH1D*)tempArray->At(0));

    tempCanvas = makeCanvas(appendFileName(tempCanvasName, fileName), tempArray, 0, "thick square logY");

    closureList->Add((tempCanvas->Clone()));
    tempCanvas->Close();
    tempArray->Clear();
    tempArrayRatios->Clear();
    delete efficiency2D;
    delete relativeError;
*/
    //----------------- Closure Test Moments ----------------------------------------
    Short_t colorsMoments[14]={kBlack,kRed,kBlue+1,kOrange+2,8,kCyan-6,kMagenta+3,kRed,kBlue,13,kBlack,16,17,12};
    Short_t markersMoments[14]={20,28,4,23,24,25,26,27,28,29,30,2,3,5};

    Short_t colorsMomentsRatios[14]={kRed,kBlack,kBlue,kOrange+2,8,kCyan-6,kMagenta+3,kRed,kBlue,13,kBlack,16,17,12};
    Short_t markersMomentsRatios[14]={28,21,4,23,24,25,26,27,28,29,30,2,3,5};


    //----------------- Closure Test Moment 1 ----------------------------------------
    tempCanvasName = "closureTestMoment1";
    tempArray->Add(getClone("momentGeneratedMC1", fileName, histos));
    tempArray->Add(getClone("momentUnfoldedMC1", fileName, histos));
    tempArray->Add(getClone("momentReweightedMC1", fileName, histos));
    setRange("X", 0, viewingRangeMoment1, (TH1D*)tempArray->At(0));


    cutHist((TH1D*)tempArray->At(0), viewingRangeMoment1);
    cutHist((TH1D*)tempArray->At(1), viewingRangeMoment1);
    cutHist((TH1D*)tempArray->At(2), viewingRangeMoment1);

    tempArrayRatios->Add(getRatio("momentUnfoldedMC1", "momentGeneratedMC1", "ratio", fileName, histos));
    tempArrayRatios->Add(new TF1("line", "1", -0.5, viewingRangeFull));
    tempArrayRatios->Add(getRatio("momentReweightedMC1", "momentGeneratedMC1", "ratio", fileName, histos));
    setRange("X", 0, viewingRangeMoment1, (TH1D*)tempArrayRatios->At(0));

    setRange("Y", 0.46, 0.89, (TH1D*)tempArray->At(0));
    setRange("Y", 0.94, 1.055, (TH1D*)tempArrayRatios->At(0));

    tempTitles[0] = "truth";
    tempTitles[1] = "unfolding method";
    tempTitles[2] = "re-weighting method";

    tempArray->Add(makeLegend(0.55, 0.27, tempArray, tempTitles));

    tempText = "";
    tempText = tempText + mc + "<|>" + colsys + ", " + erg + ", " +  eta + "<|>" + ptReach;
    tempArray->Add(makeText(0.15, 0.93, tempText,4));

    tempCanvas = makeCanvas(appendFileName(tempCanvasName, fileName), tempArray, tempArrayRatios, "thick square", colorsMoments, markersMoments, colorsMomentsRatios, markersMomentsRatios);

    closureList->Add((tempCanvas->Clone()));
    tempCanvas->Close();
    tempArray->Clear();
    tempArrayRatios->Clear();

    //----------------- Closure Test Moment 2 ----------------------------------------
    tempCanvasName = "closureTestMoment2";
    tempArray->Add(getClone("momentGeneratedMC2", fileName, histos));
    tempArray->Add(getClone("momentUnfoldedMC2", fileName, histos));
    tempArray->Add(getClone("momentReweightedMC2", fileName, histos));
    setRange("X", 0, viewingRangeMoment1, (TH1D*)tempArray->At(0));


    cutHist((TH1D*)tempArray->At(0), viewingRangeMoment1);
    cutHist((TH1D*)tempArray->At(1), viewingRangeMoment1);
    cutHist((TH1D*)tempArray->At(2), viewingRangeMoment1);

    tempArrayRatios->Add(getRatio("momentUnfoldedMC2", "momentGeneratedMC2", "ratio", fileName, histos));
    tempArrayRatios->Add(new TF1("line", "1", -0.5, viewingRangeFull));
    tempArrayRatios->Add(getRatio("momentReweightedMC2", "momentGeneratedMC2", "ratio", fileName, histos));
    setRange("X", 0, viewingRangeMoment1, (TH1D*)tempArrayRatios->At(0));

    setRange("Y", 0.25, 1.425, (TH1D*)tempArray->At(0));
    setRange("Y", 0.75, 1.25, (TH1D*)tempArrayRatios->At(0));

    tempTitles[0] = "truth";
    tempTitles[1] = "unfolding method";
    tempTitles[2] = "re-weighting method";

    tempArray->Add(makeLegend(0.55, 0.27, tempArray, tempTitles));

    tempText = "";
    tempText = tempText + mc + "<|>" + colsys + ", " + erg + ", " +  eta + "<|>" + ptReach;
    tempArray->Add(makeText(0.15, 0.93, tempText,4));

    tempCanvas = makeCanvas(appendFileName(tempCanvasName, fileName), tempArray, tempArrayRatios, "thick square", colorsMoments, markersMoments, colorsMomentsRatios, markersMomentsRatios);

    closureList->Add((tempCanvas->Clone()));
    tempCanvas->Close();
    tempArray->Clear();
    tempArrayRatios->Clear();

    //----------------- Closure Test Moment 3  ----------------------------------------
   tempCanvasName = "closureTestMoment3";
   tempArray->Add(getClone("momentGeneratedMC3", fileName, histos));
   tempArray->Add(getClone("momentUnfoldedMC3", fileName, histos));
   tempArray->Add(getClone("momentReweightedMC3", fileName, histos));
   setRange("X", 0, viewingRangeMoment1, (TH1D*)tempArray->At(0));
   setRange("Y", 0.2, 4.4, (TH1D*)tempArray->At(0));


   cutHist((TH1D*)tempArray->At(0), viewingRangeMoment1);
   cutHist((TH1D*)tempArray->At(1), viewingRangeMoment1);
   cutHist((TH1D*)tempArray->At(2), viewingRangeMoment1);

   tempArrayRatios->Add(getRatio("momentUnfoldedMC3", "momentGeneratedMC3", "ratio", fileName, histos));
   tempArrayRatios->Add(new TF1("line", "1", -0.5, viewingRangeFull));
   tempArrayRatios->Add(getRatio("momentReweightedMC3", "momentGeneratedMC3", "ratio", fileName, histos));
   setRange("X", 0, viewingRangeMoment1, (TH1D*)tempArrayRatios->At(0));
    setRange("Y", 0.75, 1.25, (TH1D*)tempArrayRatios->At(0));

    tempTitles[0] = "truth";
    tempTitles[1] = "unfolding method";
    tempTitles[2] = "re-weighting method";

   tempArray->Add(makeLegend(0.55, 0.27, tempArray, tempTitles));

   tempText = "";
   tempText = tempText + mc + "<|>" + colsys + ", " + erg + ", " +  eta + "<|>" + ptReach;
   tempArray->Add(makeText(0.15, 0.93, tempText,4));


   tempCanvas = makeCanvas(appendFileName(tempCanvasName, fileName), tempArray, tempArrayRatios, "thick square", colorsMoments, markersMoments, colorsMomentsRatios, markersMomentsRatios);

   closureList->Add((tempCanvas->Clone()));
   tempCanvas->Close();
   tempArray->Clear();
   tempArrayRatios->Clear();

/*
    //----------------- Closure Test Pt bins ----------------------------------------
    TH2D* multPtMeasuredMC = ((TH2D*) getClone("multPtMeasuredMC", fileName, histos));
    TH2D* multPtUnfoldedMC = ((TH2D*) getClone("multPtUnfoldedMC", fileName, histos));
    TH2D* multPtGeneratedMC = ((TH2D*) getClone("multPtGeneratedMC", fileName, histos));
    TH2D* responseMatrixTracksOrig = ((TH2D*) getClone("responseMatrixTracksOrig", fileName, histos));
    TH1D* startDist = responseMatrixTracksOrig->ProjectionY("d");

    for(Int_t ptBin = multPtMeasuredMC->GetYaxis()->FindBin(0.151); ptBin <= multPtMeasuredMC->GetYaxis()->FindBin(9.9); ptBin++)
    {
      multPtMeasuredMC->GetYaxis()->SetRange(ptBin,ptBin);
      multPtUnfoldedMC->GetYaxis()->SetRange(ptBin,ptBin);
      multPtGeneratedMC->GetYaxis()->SetRange(ptBin,ptBin);


      TH1D* tempProjectionMeasured = multPtMeasuredMC->ProjectionX("dummyMeasured");
      TH1D* tempProjectionUnfolded = multPtUnfoldedMC->ProjectionX("dummyUnfolded");
      TH1D* tempProjectionGenerated = multPtGeneratedMC->ProjectionX("dummyGenerated");

      TH1D* tempRatio = tempProjectionUnfolded->Clone("dummyRatio");
      tempRatio->Divide(tempProjectionGenerated);

      ostringstream tempID;
      tempID <<  "closureTestPtBin_";
      tempID << ptBin << "_" << multPtMeasuredMC->GetYaxis()->GetBinCenter(ptBin) << "GeV" ;

      tempCanvasName = tempID.str();
      tempArray->Add(tempProjectionMeasured->Clone());
      tempArray->Add(tempProjectionUnfolded->Clone());
      tempArray->Add(tempProjectionGenerated->Clone());
      setRange("X", 0, viewingRangeMoment1, (TH1D*)tempArray->At(0));

      tempTitles[0] = "measured";
      tempTitles[1] = "unfolded";
      tempTitles[2] = "generated";

      tempArray->Add(makeLegend(0.5, 0.4, tempArray, tempTitles));

      tempArrayRatios->Add(tempRatio->Clone());
      tempArrayRatios->Add(new TF1("line", "1", 0, 100));
      setRange("X", 0, viewingRangeMoment1, (TH1D*)tempArrayRatios->At(0));
      setRange("Y", 0.8, 1.2, (TH1D*)tempArrayRatios->At(0));

      ostringstream ptReachTemp;
      ptReachTemp << multPtMeasuredMC->GetYaxis()->GetBinLowEdge(ptBin) << " GeV/#it{c} - " << multPtMeasuredMC->GetYaxis()->GetBinUpEdge(ptBin) << " GeV/#it{c}";

      tempText = "";
      tempText = tempText + mc + "<|>" + colsys + ", " + erg + ", " +  eta + "<|>" + ptReachTemp.str();;
      tempArray->Add(makeText(0.5, 0.93, tempText,4));



      tempCanvas = makeCanvas(appendFileName(tempCanvasName, fileName), tempArray, tempArrayRatios, "logz square");

      closureList->Add((tempCanvas->Clone()));
      tempCanvas->Close();
      tempArray->Clear();
      tempArrayRatios->Clear();

      delete tempProjectionMeasured;
      delete tempProjectionUnfolded;
      delete tempProjectionGenerated;
      delete tempRatio;
    }
*/

    list->Add(closureList);

  }



/*
//======================= DPG playground =========================================================
  Short_t colorsTemp[14]={kBlue+2,kBlue+2,kRed+1,kRed+1,8,kCyan-6,kMagenta+3,kRed,kBlue,13,kBlack,16,17,12};
  Short_t markersTemp[14]={20,24,20,24,24,25,26,27,28,29,30,2,3,5};

  THnD* multPtMultGen = ((THnD*) getClone("multPtMultGen", fileName, histos));

  Int_t nchBin = 16;

  Int_t ptBin = multPtMultGen->GetAxis(0)->FindBin(1.05);

  multPtMultGen->GetAxis(0)->SetRange(ptBin, ptBin);
  TH1D* measuredMultDistAllNch = multPtMultGen->Projection(1);
  measuredMultDistAllNch->SetName("a");
  multPtMultGen->GetAxis(2)->SetRange(nchBin, nchBin);
  TH1D* measuredMultDistParticularNch = multPtMultGen->Projection(1);
  measuredMultDistParticularNch->SetName("b");

  tempCanvasName = "temp";
  tempArray->Add(measuredMultDistAllNch->Clone());
  tempArray->Add(measuredMultDistParticularNch->Clone());
  setTitle("Y", "1/#it{N}_{evt} 1/(2#pi #it{p}_{T}) (d^{3}#it{N})/(d#it{p}_{T}d#it{#eta}d#it{N}_{acc}) (GeV/#it{c})^{-2}", (TH1D*)tempArray->At(0));
  setRange("X", 0, 50, (TH1D*)tempArray->At(0));

  delete measuredMultDistAllNch;
  delete measuredMultDistParticularNch;

  multPtMultGen->GetAxis(2)->SetRange();
  ptBin = multPtMultGen->GetAxis(0)->FindBin(2.1);
  multPtMultGen->GetAxis(0)->SetRange(ptBin, ptBin);
  measuredMultDistAllNch = multPtMultGen->Projection(1);
  measuredMultDistAllNch->SetName("a");
  multPtMultGen->GetAxis(2)->SetRange(nchBin, nchBin);
  measuredMultDistParticularNch = multPtMultGen->Projection(1);
  measuredMultDistParticularNch->SetName("b");


//  tempArray->Add(measuredMultDistAllNch->Clone());
//  tempArray->Add(measuredMultDistParticularNch->Clone());

  tempTitles[0] = "1.0 GeV/#it{c} < #it{p}_{T} < 1.1 GeV/#it{c}";
  tempTitles[1] = "2.0 GeV/#it{c} < #it{p}_{T} < 2.2 GeV/#it{c}";

//  tempArray->Add(makeLegend(0.5, 0.9, tempArray, tempTitles));

  tempText = "1.0 GeV/#it{c} < #it{p}_{T} < 1.1 GeV/#it{c}";
  tempArray->Add(makeText(0.5, 0.9, tempText,4));


  tempCanvas = makeCanvas(appendFileName(tempCanvasName, fileName), tempArray, 0, "thick square logY dpg", colorsTemp, 0);

  list->Add((tempCanvas->Clone()));
  tempCanvas->Close();
  tempArray->Clear();
*/
  cout << "-> created " << list->GetEntries() << " plots in " << list->GetName() << endl;
  return list;
}


TList* makeComparisons(TList* histos){

  string alice = "#bf{ALICE work in progress}";
  string mc = "#bf{ALICE simulation}";

  TList* list = new TList();
  list->SetOwner();
  list->SetName("combinedPlots");

  TObjArray* tempArray = new TObjArray(1);
  tempArray->SetOwner(kFALSE);
  TObjArray* tempArrayRatios = new TObjArray(1);
  tempArrayRatios->SetOwner(kFALSE);

  string tempText;
  string tempTitles[15];
  TCanvas* tempCanvas = NULL;
  string tempCanvasName;

  //-------------------- System Size Comparison MC ------------------------------------------
  Short_t colorsSystSize[14]={0};
  Short_t markersSystSize[14]={0};
  Short_t tempColors[7] = {kWhite, kBlue+1, kMagenta+1, kRed+1, kRed, kRed, kRed};
  Short_t tempMarkers[7] = {20, 20, 24, 21, 25, 34, 28};

  colorsSystSize[0] = tempColors[0];
  colorsSystSize[1] = tempColors[1];
  colorsSystSize[2] = tempColors[1];
  colorsSystSize[3] = tempColors[2];
  colorsSystSize[4] = tempColors[2];
  colorsSystSize[5] = tempColors[3];
  colorsSystSize[6] = tempColors[3];

  markersSystSize[0] = tempMarkers[0];
  markersSystSize[1] = tempMarkers[1];
  markersSystSize[2] = tempMarkers[2];
  markersSystSize[3] = tempMarkers[3];
  markersSystSize[4] = tempMarkers[4];
  markersSystSize[5] = tempMarkers[5];
  markersSystSize[6] = tempMarkers[6];

  TH1D* dummyHist  = (TH1D*)getClone("momentUnfoldedMC1", "PbPb_5TeV_08_rebin", histos);
  dummyHist->SetName("dummyHist");
  dummyHist->Reset();

  tempCanvasName = "systemSizeMoment1MC";
  tempArray->Add(dummyHist->Clone());
  tempArray->Add(getClone("momentUnfoldedMC1", "pp_5TeV_08", histos));
  tempArray->Add(getClone("momentGeneratedMC1", "pp_5TeV_08", histos));
  tempArray->Add(getClone("momentUnfoldedMC1", "pPb_5TeV_08", histos));
  tempArray->Add(getClone("momentGeneratedMC1", "pPb_5TeV_08", histos));
  tempArray->Add(getClone("momentUnfoldedMC1", "PbPb_5TeV_08_rebin", histos));
  tempArray->Add(getClone("momentGeneratedMC1", "PbPb_5TeV_08_rebin", histos));

  setRange("X", 0, 250, (TH1D*)tempArray->At(0));
  setRange("Y", 0.45, 0.92, (TH1D*)tempArray->At(0));

  for(Int_t i = 2; i <= 8; i++){
    ((TH1D*)tempArray->At(5))->SetBinContent(i, 0);
    ((TH1D*)tempArray->At(5))->SetBinError(i, 0);
    ((TH1D*)tempArray->At(6))->SetBinContent(i, 0);
    ((TH1D*)tempArray->At(6))->SetBinError(i, 0);
  }

  cutHist((TH1D*)tempArray->At(1), 60);
  cutHist((TH1D*)tempArray->At(2), 60);
  cutHist((TH1D*)tempArray->At(3), 100);
  cutHist((TH1D*)tempArray->At(4), 100);
  cutHist((TH1D*)tempArray->At(5), 250);
  cutHist((TH1D*)tempArray->At(6), 250);

  tempTitles[1] = "pp (PYTHIA8 Monash13)";
  tempTitles[2] = "dummy";
  tempTitles[3] = "p-Pb (DPMJET)";
  tempTitles[4] = "dummy";
  tempTitles[5] = "Pb-Pb 0-80% (HIJING)";
  tempTitles[6] = "dummy";
//  tempTitles[1] = "Xe-Xe (#sqrt{#it{s}_{NN}} = 5.44 TeV)";
  tempArray->Add(makeLegend(0.43, 0.63, tempArray, tempTitles));
  tempText = "";
  tempText = tempText + "closed: unfolded<|>open:   truth";
  tempArray->Add(makeText(0.45, 0.72, tempText,2));

  tempText = "";
  tempText = tempText + mc + "<|>" +  "#sqrt{#it{s}} = 5.02 TeV , |#it{#eta}| < 0.8" + "<|>" + "0.15 GeV/#it{c} < #it{p}_{T} < 10 GeV/#it{c}";
  tempArray->Add(makeText(0.45, 0.9, tempText,3));

  tempCanvas = makeCanvas(tempCanvasName.c_str(), tempArray, 0, "thick square", colorsSystSize, markersSystSize);

  list->Add((tempCanvas->Clone()));
  tempCanvas->Close();
  tempArray->Clear();

  //-------------------- System Size Comparison Moment 2 MC ------------------------------------------
  delete dummyHist;
  TH1D* dummyHist  = (TH1D*)getClone("momentUnfoldedMC2", "PbPb_5TeV_08_rebin", histos);
  dummyHist->SetName("dummyHist");
  dummyHist->Reset();

  tempCanvasName = "systemSizeMoment2MC";
  tempArray->Add(dummyHist->Clone());
  tempArray->Add(getClone("momentUnfoldedMC2", "pp_5TeV_08", histos));
  tempArray->Add(getClone("momentGeneratedMC2", "pp_5TeV_08", histos));
  tempArray->Add(getClone("momentUnfoldedMC2", "pPb_5TeV_08", histos));
  tempArray->Add(getClone("momentGeneratedMC2", "pPb_5TeV_08", histos));
  tempArray->Add(getClone("momentUnfoldedMC2", "PbPb_5TeV_08_rebin", histos));
  tempArray->Add(getClone("momentGeneratedMC2", "PbPb_5TeV_08_rebin", histos));

  setRange("X", 0, 250, (TH1D*)tempArray->At(0));
  setRange("Y", 0.3, 1.4, (TH1D*)tempArray->At(0));

  for(Int_t i = 2; i <= 8; i++){
    ((TH1D*)tempArray->At(5))->SetBinContent(i, 0);
    ((TH1D*)tempArray->At(5))->SetBinError(i, 0);
    ((TH1D*)tempArray->At(6))->SetBinContent(i, 0);
    ((TH1D*)tempArray->At(6))->SetBinError(i, 0);
  }

  cutHist((TH1D*)tempArray->At(1), 60);
  cutHist((TH1D*)tempArray->At(2), 60);
  cutHist((TH1D*)tempArray->At(3), 100);
  cutHist((TH1D*)tempArray->At(4), 100);
  cutHist((TH1D*)tempArray->At(5), 250);
  cutHist((TH1D*)tempArray->At(6), 250);

  tempTitles[1] = "pp";
  tempTitles[2] = "dummy";
  tempTitles[3] = "p-Pb";
  tempTitles[4] = "dummy";
  tempTitles[5] = "Pb-Pb";
  tempTitles[6] = "dummy";
//  tempTitles[1] = "Xe-Xe (#sqrt{#it{s}_{NN}} = 5.44 TeV)";
  tempArray->Add(makeLegend(0.6, 0.63, tempArray, tempTitles));
  tempText = "";
  tempText = tempText + "closed: unfolded<|>open:   truth";
  tempArray->Add(makeText(0.6, 0.7, tempText,2));

  tempCanvas = makeCanvas(tempCanvasName.c_str(), tempArray, 0, "thick square", colorsSystSize, markersSystSize);

  list->Add((tempCanvas->Clone()));
  tempCanvas->Close();
  tempArray->Clear();

  //-------------------- System Size Comparison MC ------------------------------------------
  delete dummyHist;
  TH1D* dummyHist  = (TH1D*)getClone("momentUnfoldedMC3", "PbPb_5TeV_08_rebin", histos);
  dummyHist->SetName("dummyHist");
  dummyHist->Reset();

  tempCanvasName = "systemSizeMoment3MC";
  tempArray->Add(dummyHist->Clone());
  tempArray->Add(getClone("momentUnfoldedMC3", "pp_5TeV_08", histos));
  tempArray->Add(getClone("momentUnfoldedMC3", "pPb_5TeV_08", histos));
  tempArray->Add(getClone("momentUnfoldedMC3", "PbPb_5TeV_08_rebin", histos));
//  tempArray->Add(getClone("momentUnfoldedMC3", "XeXe_5TeV_08_rebin", histos));

  setRange("X", 0, 250, (TH1D*)tempArray->At(0));
  setRange("Y", 0.4, 3.5, (TH1D*)tempArray->At(0));

  for(Int_t i = 2; i <= 8; i++){
    ((TH1D*)tempArray->At(3))->SetBinContent(i, 0);
    ((TH1D*)tempArray->At(3))->SetBinError(i, 0);
  }

  cutHist((TH1D*)tempArray->At(1), 60);
  cutHist((TH1D*)tempArray->At(2), 100);
  cutHist((TH1D*)tempArray->At(3), 250);

  tempTitles[1] = "pp";
  tempTitles[2] = "p-Pb";
  tempTitles[3] = "Pb-Pb";
//  tempTitles[1] = "Xe-Xe (#sqrt{#it{s}_{NN}} = 5.44 TeV)";
  tempArray->Add(makeLegend(0.7, 0.35, tempArray, tempTitles));

  tempText = "";
  tempText = tempText + mc + "<|>" +  "#sqrt{#it{s}} = 5.02 TeV , |#it{#eta}| < 0.8" + "<|>" + "0.15 GeV/#it{c} < #it{p}_{T} < 10 GeV/#it{c}" + "<|>"  + "PbPb: 0\% - 80\% cent";
  tempArray->Add(makeText(0.5, 0.9, tempText,4));

  tempCanvas = makeCanvas(tempCanvasName.c_str(), tempArray, 0, "thick square", colorsSystSize, markersSystSize);

  list->Add((tempCanvas->Clone()));
  tempCanvas->Close();
  tempArray->Clear();

  //-------------------- System Size Comparison data ------------------------------------------

  colorsSystSize[0] = tempColors[0];
  colorsSystSize[1] = tempColors[1];
  colorsSystSize[2] = tempColors[2];
  colorsSystSize[3] = tempColors[3];
  colorsSystSize[4] = tempColors[4];
  colorsSystSize[5] = tempColors[3];
  colorsSystSize[6] = tempColors[3];

  markersSystSize[0] = tempMarkers[0];
  markersSystSize[1] = tempMarkers[1];
  markersSystSize[2] = tempMarkers[3];
  markersSystSize[3] = tempMarkers[5];
  markersSystSize[4] = tempMarkers[4];
  markersSystSize[5] = tempMarkers[5];
  markersSystSize[6] = tempMarkers[6];

  delete dummyHist;
  TH1D* dummyHist  = (TH1D*)getClone("momentUnfolded1", "PbPb_5TeV_08_rebin", histos);
  dummyHist->SetName("dummyHist");
  dummyHist->Reset();

  tempCanvasName = "systemSizeMoment1";
  tempArray->Add(dummyHist->Clone());
  tempArray->Add(getClone("momentUnfolded1", "pp_5TeV_08", histos));
  tempArray->Add(getClone("momentUnfolded1", "pPb_5TeV_08", histos));
  tempArray->Add(getClone("momentUnfolded1", "PbPb_5TeV_08_rebin", histos));

  setRange("X", 0, 250, (TH1D*)tempArray->At(0));
  setRange("Y", 0.45, 0.92, (TH1D*)tempArray->At(0));

//  ((TH1D*)tempArray->At(0))->SetBinContent(2, 0);
//  ((TH1D*)tempArray->At(0))->SetBinError(2, 0);

  cutHist((TH1D*)tempArray->At(1), 60);
  cutHist((TH1D*)tempArray->At(2), 100);
  cutHist((TH1D*)tempArray->At(3), 250);

  tempTitles[0] = "dummy";
  tempTitles[1] = "pp";
  tempTitles[2] = "p-Pb";
  tempTitles[3] = "Pb-Pb";
//  tempTitles[1] = "Xe-Xe (#sqrt{#it{s}_{NN}} = 5.44 TeV)";
  tempArray->Add(makeLegend(0.7, 0.35, tempArray, tempTitles));

  tempText = "";
  tempText = tempText + alice + "<|>" +  "#sqrt{#it{s}_{NN}} = 5.02 TeV , |#it{#eta}| < 0.8" + "<|>" + "0.15 GeV/#it{c} < #it{p}_{T} < 10 GeV/#it{c}" + "<|>"  + "Pb-Pb: 0\% - 80\% cent";
  tempArray->Add(makeText(0.5, 0.9, tempText,4));

  tempCanvas = makeCanvas(tempCanvasName.c_str(), tempArray, 0, "thick square", colorsSystSize, markersSystSize);

  list->Add((tempCanvas->Clone()));
  tempCanvas->Close();
  tempArray->Clear();

  //-------------------- System Size Comparison data ------------------------------------------
  delete dummyHist;
  TH1D* dummyHist  = (TH1D*)getClone("momentUnfolded2", "PbPb_5TeV_08_rebin", histos);
  dummyHist->SetName("dummyHist");
  dummyHist->Reset();

  tempCanvasName = "systemSizeMoment2";
  tempArray->Add(dummyHist->Clone());
  tempArray->Add(getClone("momentUnfolded2", "pp_5TeV_08", histos));
  tempArray->Add(getClone("momentUnfolded2", "pPb_5TeV_08", histos));
  tempArray->Add(getClone("momentUnfolded2", "PbPb_5TeV_08_rebin", histos));
//  tempArray->Add(getClone("momentUnfolded2", "XeXe_5TeV_08_rebin", histos));

  setRange("X", 0, 250, (TH1D*)tempArray->At(0));
  setRange("Y", 0.3, 1.4, (TH1D*)tempArray->At(0));

//  ((TH1D*)tempArray->At(0))->SetBinContent(2, 0);
//  ((TH1D*)tempArray->At(0))->SetBinError(2, 0);
//  ((TH1D*)tempArray->At(1))->SetBinContent(2, 0);
//  ((TH1D*)tempArray->At(1))->SetBinError(2, 0);

  cutHist((TH1D*)tempArray->At(1), 60);
  cutHist((TH1D*)tempArray->At(2), 100);
  cutHist((TH1D*)tempArray->At(3), 250);

  tempTitles[0] = "dummy";
  tempTitles[1] = "pp";
  tempTitles[2] = "p-Pb";
  tempTitles[3] = "Pb-Pb";
//  tempTitles[1] = "Xe-Xe (#sqrt{#it{s}_{NN}} = 5.44 TeV)";
  tempArray->Add(makeLegend(0.7, 0.35, tempArray, tempTitles));

  tempText = "";
  tempText = tempText + alice + "<|>" +  "#sqrt{#it{s}_{NN}} = 5.02 TeV , |#it{#eta}| < 0.8" + "<|>" + "0.15 GeV/#it{c} < #it{p}_{T} < 10 GeV/#it{c}" + "<|>"  + "PbPb: 0\% - 80\% cent";
  tempArray->Add(makeText(0.5, 0.9, tempText,4));

  tempCanvas = makeCanvas(tempCanvasName.c_str(), tempArray, 0, "thick square", colorsSystSize, markersSystSize);

  list->Add((tempCanvas->Clone()));
  tempCanvas->Close();
  tempArray->Clear();

  //-------------------- System Size Comparison data ------------------------------------------
  delete dummyHist;
  TH1D* dummyHist  = (TH1D*)getClone("momentUnfolded3", "PbPb_5TeV_08_rebin", histos);
  dummyHist->SetName("dummyHist");
  dummyHist->Reset();

  tempCanvasName = "systemSizeMoment3";
  tempArray->Add(dummyHist->Clone());
  tempArray->Add(getClone("momentUnfolded3", "pp_5TeV_08", histos));
  tempArray->Add(getClone("momentUnfolded3", "pPb_5TeV_08", histos));
  tempArray->Add(getClone("momentUnfolded3", "PbPb_5TeV_08_rebin", histos));
//  tempArray->Add(getClone("momentUnfolded3", "XeXe_5TeV_08_rebin", histos));

  setRange("X", 0, 250, (TH1D*)tempArray->At(0));
  setRange("Y", 0.4, 3.5, (TH1D*)tempArray->At(0));

//  ((TH1D*)tempArray->At(0))->SetBinContent(2, 0);
//  ((TH1D*)tempArray->At(0))->SetBinError(2, 0);
//  ((TH1D*)tempArray->At(1))->SetBinContent(2, 0);
//  ((TH1D*)tempArray->At(1))->SetBinError(2, 0);

  cutHist((TH1D*)tempArray->At(1), 60);
  cutHist((TH1D*)tempArray->At(2), 100);
  cutHist((TH1D*)tempArray->At(3), 250);

  tempTitles[0] = "dummy";
  tempTitles[1] = "pp";
  tempTitles[2] = "p-Pb";
  tempTitles[3] = "Pb-Pb";
//  tempTitles[1] = "Xe-Xe (#sqrt{#it{s}_{NN}} = 5.44 TeV)";
  tempArray->Add(makeLegend(0.7, 0.35, tempArray, tempTitles));

  tempText = "";
  tempText = tempText + alice + "<|>" +  "#sqrt{#it{s}_{NN}} = 5.02 TeV , |#it{#eta}| < 0.8" + "<|>" + "0.15 GeV/#it{c} < #it{p}_{T} < 10 GeV/#it{c}" + "<|>"  + "PbPb: 0\% - 80\% cent";
  tempArray->Add(makeText(0.5, 0.9, tempText,4));

  tempCanvas = makeCanvas(tempCanvasName.c_str(), tempArray, 0, "thick square", colorsSystSize, markersSystSize);

  list->Add((tempCanvas->Clone()));
  tempCanvas->Close();
  tempArray->Clear();



  //-------------------- System Size Comparison data ------------------------------------------
  delete dummyHist;
  TH1D* dummyHist  = (TH1D*)getClone("momentUnfolded1", "PbPb_5TeV_08_rebin", histos);
  dummyHist->SetName("dummyHist");
  dummyHist->Reset();

  tempCanvasName = "systemSizeMoment1_withPub";
  tempArray->Add(dummyHist->Clone());
  tempArray->Add(getClone("momentUnfolded1", "pp_5TeV_08", histos));
  tempArray->Add(getClone("momentReweighted1Sys_pp_7TeV_03", "pub", histos));
  tempArray->Add(getClone("momentUnfolded1", "pPb_5TeV_08", histos));
  tempArray->Add(getClone("momentReweighted1Sys_pPb_5TeV_03", "pub", histos));
  tempArray->Add(getClone("momentUnfolded1", "PbPb_5TeV_08_rebin", histos));
  tempArray->Add(getClone("momentReweighted1Sys_PbPb_2TeV_03", "pub", histos));

  setRange("X", 0, 100, (TH1D*)tempArray->At(0));
  setRange("Y", 0.45, 0.92, (TH1D*)tempArray->At(0));

//  ((TH1D*)tempArray->At(0))->SetBinContent(2, 0);
//  ((TH1D*)tempArray->At(0))->SetBinError(2, 0);

  cutHist((TH1D*)tempArray->At(1), 60);
  cutHist((TH1D*)tempArray->At(3), 100);
  cutHist((TH1D*)tempArray->At(4), 100);

  tempTitles[0] = "dummy";
  tempTitles[1] = "pp #sqrt{#it{s}} = 5.02 TeV , |#it{#eta}| < 0.8";
  tempTitles[2] = "pp #sqrt{#it{s}} = 7 TeV , |#it{#eta}| < 0.3";
  tempTitles[3] = "p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV , |#it{#eta}| < 0.8";
  tempTitles[4] = "p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV , |#it{#eta}| < 0.3";
  tempTitles[5] = "Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV , |#it{#eta}| < 0.8";
  tempTitles[6] = "Pb-Pb #sqrt{#it{s}_{NN}} = 2.76 TeV , |#it{#eta}| < 0.3";
//  tempTitles[1] = "Xe-Xe (#sqrt{#it{s}_{NN}} = 5.44 TeV)";
  tempArray->Add(makeLegend(0.5, 0.35, tempArray, tempTitles));

  tempText = "";
  tempText = tempText + alice + "<|>" + "0.15 GeV/#it{c} < #it{p}_{T} < 10 GeV/#it{c}";
  tempArray->Add(makeText(0.5, 0.9, tempText,2));

  tempCanvas = makeCanvas(tempCanvasName.c_str(), tempArray, 0, "thick square", colorsSystSize, markersSystSize);

  list->Add((tempCanvas->Clone()));
  tempCanvas->Close();
  tempArray->Clear();


  //-------------------- Energy Comparison MC ------------------------------------------
  Short_t colorsErg[14]={0};
  Short_t markersErg[14]={0};
  Short_t tempColorsErg[4] = {kRed, kBlue, kGreen+2, kMagenta+1};
  Short_t tempMarkersErg[4] = {20, 21, 34, 20};

  tempCanvasName = "energyDependenceMoment1MC";
  tempArray->Add(getClone("momentUnfoldedMC1", "pp_5TeV_08", histos));
//  tempArray->Add(getClone("momentUnfoldedMC1", "pp_7TeV_08", histos));
  tempArray->Add(getClone("momentUnfoldedMC1", "pp_8TeV_08", histos));
  tempArray->Add(getClone("momentUnfoldedMC1", "pp_13TeV_08", histos));


  colorsErg[0] = tempColorsErg[0];
  colorsErg[1] = tempColorsErg[2];
  colorsErg[2] = tempColorsErg[3];
  colorsErg[3] = tempColorsErg[3];

  markersErg[0] = tempMarkersErg[0];
  markersErg[1] = tempMarkersErg[2];
  markersErg[2] = tempMarkersErg[3];
  markersErg[3] = tempMarkersErg[3];

  setRange("X", 0, 60, (TH1D*)tempArray->At(0));
  setRange("Y", 0.44, 0.95, (TH1D*)tempArray->At(0));

  cutHist((TH1D*)tempArray->At(0), 60);
  cutHist((TH1D*)tempArray->At(1), 60);

  tempTitles[0] = "#sqrt{#it{s}} = 5.02 TeV";
  tempTitles[1] = "#sqrt{#it{s}} = 8 TeV";
  tempTitles[2] = "#sqrt{#it{s}} = 13 TeV";

//  tempTitles[0] = "#sqrt{#it{s}} = 5.02 TeV PYTHIA8 Monash13";
//  tempTitles[1] = "#sqrt{#it{s}} = 7 TeV PYTHIA6 Perugia11";
//  tempTitles[2] = "#sqrt{#it{s}} = 8 TeV PYTHIA8 Monash13";
//  tempTitles[3] = "#sqrt{#it{s}} = 13 TeV PYTHIA8 Monash13";
  tempArray->Add(makeLegend(0.15, 0.75, tempArray, tempTitles));

  tempText = "";
  tempText = tempText + mc + "<|>" +  "pp , |#it{#eta}| < 0.8" + "<|>" + "0.15 GeV/#it{c} < #it{p}_{T} < 10 GeV/#it{c}";
  tempArray->Add(makeText(0.15, 0.9, tempText,3));

  tempCanvas = makeCanvas(tempCanvasName.c_str(), tempArray, 0, "thick square", colorsErg, markersErg);

  list->Add((tempCanvas->Clone()));
  tempCanvas->Close();
  tempArray->Clear();

  //-------------------- Energy Comparison Data ------------------------------------------

  colorsErg[0] = tempColorsErg[0];
  colorsErg[1] = tempColorsErg[0];
  colorsErg[2] = tempColorsErg[1];
  colorsErg[3] = tempColorsErg[1];
  colorsErg[4] = tempColorsErg[2];
  colorsErg[5] = tempColorsErg[2];
  colorsErg[6] = tempColorsErg[3];
  colorsErg[7] = tempColorsErg[3];

  markersErg[0] = tempMarkersErg[0];
  markersErg[1] = tempMarkersErg[0];
  markersErg[2] = tempMarkersErg[1];
  markersErg[3] = tempMarkersErg[1];
  markersErg[4] = tempMarkersErg[2];
  markersErg[5] = tempMarkersErg[2];
  markersErg[6] = tempMarkersErg[3];
  markersErg[7] = tempMarkersErg[3];

  tempCanvasName = "energyDependenceMoment1";
  tempArray->Add(getClone("momentUnfolded1", "pp_5TeV_08", histos));
  tempArray->Add(getClone("momentUnfolded1SystErr", "pp_5TeV_08", histos));
  tempArray->Add(getClone("momentUnfolded1", "pp_7TeV_08", histos));
  tempArray->Add(getClone("momentUnfolded1SystErr", "pp_7TeV_08", histos));
  tempArray->Add(getClone("momentUnfolded1", "pp_8TeV_08", histos));
  tempArray->Add(getClone("momentUnfolded1SystErr", "pp_8TeV_08", histos));
  tempArray->Add(getClone("momentUnfolded1", "pp_13TeV_08", histos));
  tempArray->Add(getClone("momentUnfolded1SystErr", "pp_13TeV_08", histos));

  setRange("X", 0, 75, (TH1D*)tempArray->At(0));
  setRange("Y", 0.44, 0.95, (TH1D*)tempArray->At(0));
  cutHist((TH1D*)tempArray->At(0), 75);
  cutHist((TH1D*)tempArray->At(1), 75);

  tempTitles[0] = "#sqrt{#it{s}} = 5.02 TeV";
  tempTitles[1] = "dummy";
  tempTitles[2] = "#sqrt{#it{s}} = 7 TeV";
  tempTitles[3] = "dummy";
  tempTitles[4] = "#sqrt{#it{s}} = 8 TeV";
  tempTitles[6] = "dummy";
  tempTitles[6] = "#sqrt{#it{s}} = 13 TeV";
  tempTitles[7] = "dummy";
  tempArray->Add(makeLegend(0.2, 0.5, tempArray, tempTitles));

  tempText = "";
  tempText = tempText + alice + "<|>" +  "pp , |#it{#eta}| < 0.8" + "<|>" + "0.15 GeV/#it{c} < #it{p}_{T} < 10 GeV/#it{c}";
  tempArray->Add(makeText(0.15, 0.9, tempText,3));

  tempCanvas = makeCanvas(tempCanvasName.c_str(), tempArray, 0, "thick square", colorsErg, markersErg);

  list->Add((tempCanvas->Clone()));
  tempCanvas->Close();
  tempArray->Clear();


  tempCanvasName = "energyDependenceMoment2";
  tempArray->Add(getClone("momentUnfolded2", "pp_5TeV_08", histos));
  tempArray->Add(getClone("momentUnfolded2SystErr", "pp_5TeV_08", histos));
  tempArray->Add(getClone("momentUnfolded2", "pp_7TeV_08", histos));
  tempArray->Add(getClone("momentUnfolded2SystErr", "pp_7TeV_08", histos));
  tempArray->Add(getClone("momentUnfolded2", "pp_8TeV_08", histos));
  tempArray->Add(getClone("momentUnfolded2SystErr", "pp_8TeV_08", histos));
  tempArray->Add(getClone("momentUnfolded2", "pp_13TeV_08", histos));
  tempArray->Add(getClone("momentUnfolded2SystErr", "pp_13TeV_08", histos));

  setRange("X", 0, 60, (TH1D*)tempArray->At(0));
  setRange("Y", 0.3, 1.45, (TH1D*)tempArray->At(0));
  cutHist((TH1D*)tempArray->At(0), 60);
  cutHist((TH1D*)tempArray->At(1), 60);


  tempTitles[0] = "#sqrt{#it{s}} = 5.02 TeV";
  tempTitles[1] = "dummy";
  tempTitles[2] = "#sqrt{#it{s}} = 7 TeV";
  tempTitles[3] = "dummy";
  tempTitles[4] = "#sqrt{#it{s}} = 8 TeV";
  tempTitles[6] = "dummy";
  tempTitles[6] = "#sqrt{#it{s}} = 13 TeV";
  tempTitles[7] = "dummy";
  tempArray->Add(makeLegend(0.15, 0.75, tempArray, tempTitles));

  tempText = "";
  tempText = tempText + alice + "<|>" +  "pp , |#it{#eta}| < 0.8" + "<|>" + "0.15 GeV/#it{c} < #it{p}_{T} < 10 GeV/#it{c}";
  tempArray->Add(makeText(0.15, 0.9, tempText,3));

  tempCanvas = makeCanvas(tempCanvasName.c_str(), tempArray, 0, "thick square", colorsErg, markersErg);

  list->Add((tempCanvas->Clone()));
  tempCanvas->Close();
  tempArray->Clear();


  tempCanvasName = "energyDependenceMoment3";
  tempArray->Add(getClone("momentUnfolded3", "pp_5TeV_08", histos));
  tempArray->Add(getClone("momentUnfolded3SystErr", "pp_5TeV_08", histos));
  tempArray->Add(getClone("momentUnfolded3", "pp_7TeV_08", histos));
  tempArray->Add(getClone("momentUnfolded3SystErr", "pp_7TeV_08", histos));
  tempArray->Add(getClone("momentUnfolded3", "pp_8TeV_08", histos));
  tempArray->Add(getClone("momentUnfolded3SystErr", "pp_8TeV_08", histos));
  tempArray->Add(getClone("momentUnfolded3", "pp_13TeV_08", histos));
  tempArray->Add(getClone("momentUnfolded3SystErr", "pp_13TeV_08", histos));

  setRange("X", 0, 60, (TH1D*)tempArray->At(0));
  setRange("Y", 0.2, 4.2, (TH1D*)tempArray->At(0));
  cutHist((TH1D*)tempArray->At(0), 60);
  cutHist((TH1D*)tempArray->At(1), 60);


  tempTitles[0] = "#sqrt{#it{s}} = 5.02 TeV";
  tempTitles[1] = "dummy";
  tempTitles[2] = "#sqrt{#it{s}} = 7 TeV";
  tempTitles[3] = "dummy";
  tempTitles[4] = "#sqrt{#it{s}} = 8 TeV";
  tempTitles[6] = "dummy";
  tempTitles[6] = "#sqrt{#it{s}} = 13 TeV";
  tempTitles[7] = "dummy";
  tempArray->Add(makeLegend(0.15, 0.75, tempArray, tempTitles));

  tempText = "";
  tempText = tempText + alice + "<|>" +  "pp , |#it{#eta}| < 0.8" + "<|>" + "0.15 GeV/#it{c} < #it{p}_{T} < 10 GeV/#it{c}";
  tempArray->Add(makeText(0.15, 0.9, tempText,3));

  tempCanvas = makeCanvas(tempCanvasName.c_str(), tempArray, 0, "thick square", colorsErg, markersErg);

  list->Add((tempCanvas->Clone()));
  tempCanvas->Close();
  tempArray->Clear();


// ratio to 8 tev
tempCanvasName = "energyRatios";
tempArray->Add(getClone("momentUnfolded1", "pp_5TeV_08", histos));
tempArray->Add(getClone("momentUnfolded1", "pp_8TeV_08", histos));
tempArray->Add(getClone("momentUnfolded1", "pp_13TeV_08", histos));

setRange("X", 0, 60, (TH1D*)tempArray->At(0));
setRange("Y", 0.44, 0.95, (TH1D*)tempArray->At(0));
cutHist((TH1D*)tempArray->At(0), 60);
cutHist((TH1D*)tempArray->At(1), 60);

TH1D* base5 = getClone("momentUnfolded1", "pp_5TeV_08", histos);
TH1D* ratio8 = getClone("momentUnfolded1", "pp_8TeV_08", histos);
ratio8->Divide(base5);
TH1D* ratio13 = getClone("momentUnfolded1", "pp_13TeV_08", histos);
ratio13->Divide(base5);

tempArrayRatios->Add(ratio8);
tempArrayRatios->Add(ratio13);
//tempArrayRatios->Add(new TF1("line", "1", -0.5, 99.5));
setRange("X", 0, 60, (TH1D*)tempArrayRatios->At(0));


tempTitles[0] = "#sqrt{#it{s}} = 5.02 TeV";
tempTitles[1] = "#sqrt{#it{s}} = 8 TeV";
tempTitles[2] = "#sqrt{#it{s}} = 13 TeV";
tempArray->Add(makeLegend(0.15, 0.75, tempArray, tempTitles));

tempText = "";
tempText = tempText + alice + "<|>" +  "pp , |#it{#eta}| < 0.8" + "<|>" + "0.15 GeV/#it{c} < #it{p}_{T} < 10 GeV/#it{c}";
tempArray->Add(makeText(0.15, 0.9, tempText,3));

tempCanvas = makeCanvas(tempCanvasName.c_str(), tempArray, tempArrayRatios, "thick square", colorsErg, markersErg);

list->Add((tempCanvas->Clone()));
tempCanvas->Close();
tempArray->Clear();
tempArrayRatios->Clear();

  //_________________________________________________________________________________

  Short_t colorsSyst[14]={kPink+8,kGreen+2,28,kOrange+2,8,kCyan-6,kMagenta+3,kRed,kBlue,13,kBlack,16,17,12};
  Short_t markersSyst[14]={20,21,22,23,24,25,26,27,28,29,30,2,3,5};

  tempCanvasName = "iterationDependence";

  tempArray->Add(getClone("momentUnfolded1_nIter_1", "pp_5TeV_08_iterationCheck", histos));
  tempArray->Add(getClone("momentUnfolded1_nIter_2", "pp_5TeV_08_iterationCheck", histos));
  tempArray->Add(getClone("momentUnfolded1_nIter_3", "pp_5TeV_08_iterationCheck", histos));
  tempArray->Add(getClone("momentUnfolded1_nIter_4", "pp_5TeV_08_iterationCheck", histos));
  tempArray->Add(getClone("momentUnfolded1_nIter_5", "pp_5TeV_08_iterationCheck", histos));
  tempArray->Add(getClone("momentUnfolded1_nIter_6", "pp_5TeV_08_iterationCheck", histos));
  tempArray->Add(getClone("momentUnfolded1_nIter_7", "pp_5TeV_08_iterationCheck", histos));
  tempArray->Add(getClone("momentUnfolded1_nIter_8", "pp_5TeV_08_iterationCheck", histos));
  tempArray->Add(getClone("momentUnfolded1_nIter_9", "pp_5TeV_08_iterationCheck", histos));
  tempArray->Add(getClone("momentUnfolded1_nIter_10", "pp_5TeV_08_iterationCheck", histos));

  tempTitles[0] = "  1  ";
  tempTitles[1] = "  2  ";
  tempTitles[2] = "  3  ";
  tempTitles[3] = "  4  ";
  tempTitles[4] = "  5  ";
  tempTitles[5] = "  6  ";
  tempTitles[6] = "  7  ";
  tempTitles[7] = "  8  ";
  tempTitles[8] = "  9  ";
  tempTitles[9] = "  10 ";

//  setRange("X", 0, 60, (TH1D*)tempArray->At(0));
  setRange("Y", 0.48, 0.83, (TH1D*)tempArray->At(0));

  tempArray->Add(makeLegend(0.5, 0.5, tempArray, tempTitles,0.6));

  tempText = "";
  tempText = tempText + alice + "<|>" +  "pp , #sqrt{#it{s}} = 5.02 TeV , |#it{#eta}| < 0.8" + "<|>" + "0.15 GeV/#it{c} < #it{p}_{T} < 10 GeV/#it{c}";
  tempArray->Add(makeText(0.5, 0.8, tempText, 3));

  tempText = "n_{iter} = ";
  tempArray->Add(makeText(0.5, 0.6, tempText, 1));

  tempCanvas = makeCanvas(tempCanvasName.c_str(), tempArray, 0, "2xthick square", colorsSyst, markersSyst);

  list->Add((tempCanvas->Clone()));
  tempCanvas->Close();
  tempArray->Clear();
  //_________________________________________________________________________________


  TFile* publFile =  new TFile("Histograms/paperMeanPt.root","READ");
  TGraphErrors* publMoment1FromSpectraDnchDeta = (TGraphErrors*)publFile->FindObjectAny("publMoment1FromSpectraDnchDeta");

  //-------------------- Mean Unfolded ------------------------------------------------------
  tempCanvasName = "moment1dnchdeta";
  tempArray->Add(compressAxis((TH1D*)getClone("momentUnfolded1", "PbPb_5TeV_08_rebin", histos), 1/1.8));
  tempArray->Add(publMoment1FromSpectraDnchDeta->Clone());

  setRange("X", 0, 2000, (TH1D*)tempArray->At(0));
  setRange("Y", 0.48, 0.75, (TH1D*)tempArray->At(0));

  tempTitles[0] = "Mario";
  tempTitles[1] = "Julius";
  tempArray->Add(makeLegend(0.5, 0.5, tempArray, tempTitles,0.6));

  tempCanvas = makeCanvas(tempCanvasName.c_str(), tempArray, 0, "thick square");

  list->Add((tempCanvas->Clone()));
  tempCanvas->Close();
  tempArray->Clear();


//------------------ QM playground -----------------------------------------------
//----------------------secondMC -------------------------------------------------
tempCanvasName = "ratioToDataMoment1";
tempArray->Add(getClone("momentGeneratedMC1", "pp_5TeV_08_pass3", histos));
tempArray->Add(getClone("momentGeneratedMC1", "pp_5TeV_08_pass3_Pythia6", histos));
tempArray->Add(getClone("momentUnfolded1", "pp_5TeV_08_pass3", histos));
setRange("X", 0, 75, (TH1D*)tempArray->At(0));

tempTitles[0] = "Monash13 truth";
tempTitles[1] = "Perugia11 truth";
tempTitles[2] = "Data unfolded";
tempArray->Add(makeLegend(0.5, 0.35, tempArray, tempTitles));

tempText = "";
tempText = tempText + alice + "<|>" +  "pp , #sqrt{#it{s}} = 5.02 TeV , |#it{#eta}| < 0.8" + "<|>" + "0.15 GeV/#it{c} < #it{p}_{T} < 10 GeV/#it{c}";
tempArray->Add(makeText(0.5, 0.8, tempText, 3));

tempArrayRatios->Add(getRatio("momentGeneratedMC1", "momentUnfolded1", "mc/data", "pp_5TeV_08_pass3", histos, "pp_5TeV_08_pass3"));
tempArrayRatios->Add(getRatio("momentGeneratedMC1", "momentUnfolded1", "mc/data", "pp_5TeV_08_pass3_Pythia6", histos, "pp_5TeV_08_pass3"));
tempArrayRatios->Add(new TF1("line", "1", -0.5, 99.5));

setRange("X", 0, 60, (TH1D*)tempArrayRatios->At(0));


tempCanvas = makeCanvas(tempCanvasName.c_str(), tempArray, tempArrayRatios, "thick square");

list->Add((tempCanvas->Clone()));
tempCanvas->Close();
tempArray->Clear();
tempArrayRatios->Clear();

tempCanvasName = "ratioToDataMoment2";
tempArray->Add(getClone("momentGeneratedMC2", "pp_5TeV_08_pass3", histos));
tempArray->Add(getClone("momentGeneratedMC2", "pp_5TeV_08_pass3_Pythia6", histos));
tempArray->Add(getClone("momentUnfolded2", "pp_5TeV_08_pass3", histos));
setRange("X", 0, 75, (TH1D*)tempArray->At(0));

tempTitles[0] = "Monash13 truth";
tempTitles[1] = "Perugia11 truth";
tempTitles[2] = "Data unfolded";
tempArray->Add(makeLegend(0.5, 0.35, tempArray, tempTitles));

tempText = "";
tempText = tempText + alice + "<|>" +  "pp , #sqrt{#it{s}} = 5.02 TeV , |#it{#eta}| < 0.8" + "<|>" + "0.15 GeV/#it{c} < #it{p}_{T} < 10 GeV/#it{c}";
tempArray->Add(makeText(0.5, 0.8, tempText, 3));

tempArrayRatios->Add(getRatio("momentGeneratedMC2", "momentUnfolded2", "mc/data", "pp_5TeV_08_pass3", histos, "pp_5TeV_08_pass3"));
tempArrayRatios->Add(getRatio("momentGeneratedMC2", "momentUnfolded2", "mc/data", "pp_5TeV_08_pass3_Pythia6", histos, "pp_5TeV_08_pass3"));
tempArrayRatios->Add(new TF1("line", "1", -0.5, 99.5));

setRange("X", 0, 60, (TH1D*)tempArrayRatios->At(0));


tempCanvas = makeCanvas(tempCanvasName.c_str(), tempArray, tempArrayRatios, "thick square");

list->Add((tempCanvas->Clone()));
tempCanvas->Close();
tempArray->Clear();
tempArrayRatios->Clear();

tempCanvasName = "ratioToDataMoment3";
tempArray->Add(getClone("momentGeneratedMC3", "pp_5TeV_08_pass3", histos));
tempArray->Add(getClone("momentGeneratedMC3", "pp_5TeV_08_pass3_Pythia6", histos));
tempArray->Add(getClone("momentUnfolded3", "pp_5TeV_08_pass3", histos));
setRange("X", 0, 75, (TH1D*)tempArray->At(0));

tempTitles[0] = "Monash13 truth";
tempTitles[1] = "Perugia11 truth";
tempTitles[2] = "Data unfolded";
tempArray->Add(makeLegend(0.5, 0.35, tempArray, tempTitles));

tempText = "";
tempText = tempText + alice + "<|>" +  "pp , #sqrt{#it{s}} = 5.02 TeV , |#it{#eta}| < 0.8" + "<|>" + "0.15 GeV/#it{c} < #it{p}_{T} < 10 GeV/#it{c}";
tempArray->Add(makeText(0.5, 0.8, tempText, 3));

tempArrayRatios->Add(getRatio("momentGeneratedMC3", "momentUnfolded3", "mc/data", "pp_5TeV_08_pass3", histos, "pp_5TeV_08_pass3"));
tempArrayRatios->Add(getRatio("momentGeneratedMC3", "momentUnfolded3", "mc/data", "pp_5TeV_08_pass3_Pythia6", histos, "pp_5TeV_08_pass3"));
tempArrayRatios->Add(new TF1("line", "1", -0.5, 99.5));

setRange("X", 0, 60, (TH1D*)tempArrayRatios->At(0));


tempCanvas = makeCanvas(tempCanvasName.c_str(), tempArray, tempArrayRatios, "thick square");

list->Add((tempCanvas->Clone()));
tempCanvas->Close();
tempArray->Clear();
tempArrayRatios->Clear();




tempCanvasName = "multDistMeasuredMC_MonashPerugia";
tempArray->Add(getClone("multDistMeasuredMC", "pp_5TeV_08_pass3_Pythia6", histos));
tempArray->Add(getClone("multDistMeasuredMC", "pp_5TeV_08_pass3", histos));
setRange("X", 0, 60, (TH1D*)tempArray->At(0));

tempTitles[0] = "Perugia11";
tempTitles[1] = "Monash13";
tempArray->Add(makeLegend(0.5, 0.35, tempArray, tempTitles));

tempArrayRatios->Add(getRatio("multDistMeasuredMC", "multDistMeasuredMC", "Peru/Mona", "pp_5TeV_08_pass3_Pythia6", histos, "pp_5TeV_08_pass3"));
tempArrayRatios->Add(new TF1("line", "1", -0.5, 99.5));

setRange("X", 0, 60, (TH1D*)tempArrayRatios->At(0));


tempCanvas = makeCanvas(tempCanvasName.c_str(), tempArray, tempArrayRatios, "thick square logY");

list->Add((tempCanvas->Clone()));
tempCanvas->Close();
tempArray->Clear();
tempArrayRatios->Clear();


tempCanvasName = "multDistGeneratedMC_MonashPerugia";
tempArray->Add(getClone("multDistGeneratedMC", "pp_5TeV_08_pass3_Pythia6", histos));
tempArray->Add(getClone("multDistGeneratedMC", "pp_5TeV_08_pass3", histos));
setRange("X", 0, 80, (TH1D*)tempArray->At(0));

tempTitles[0] = "Perugia11";
tempTitles[1] = "Monash13";
tempArray->Add(makeLegend(0.5, 0.35, tempArray, tempTitles));

tempArrayRatios->Add(getRatio("multDistGeneratedMC", "multDistGeneratedMC", "Peru/Mona", "pp_5TeV_08_pass3_Pythia6", histos, "pp_5TeV_08_pass3"));
tempArrayRatios->Add(new TF1("line", "1", -0.5, 99.5));
setRange("X", 0, 80, (TH1D*)tempArrayRatios->At(0));


tempCanvas = makeCanvas(tempCanvasName.c_str(), tempArray, tempArrayRatios, "thick square logY");

list->Add((tempCanvas->Clone()));
tempCanvas->Close();
tempArray->Clear();
tempArrayRatios->Clear();


tempCanvasName = "ptDistMeasuredMC_MonashPerugia";
tempArray->Add(getClone("ptDistMeasuredMC", "pp_5TeV_08_pass3_Pythia6", histos));
tempArray->Add(getClone("ptDistMeasuredMC", "pp_5TeV_08_pass3", histos));

tempTitles[0] = "Perugia11";
tempTitles[1] = "Monash13";
tempArray->Add(makeLegend(0.5, 0.35, tempArray, tempTitles));

tempArrayRatios->Add(getRatio("ptDistMeasuredMC", "ptDistMeasuredMC", "Peru/Mona", "pp_5TeV_08_pass3_Pythia6", histos, "pp_5TeV_08_pass3"));
tempArrayRatios->Add(new TF1("line", "1", -0.5, 99.5));


tempCanvas = makeCanvas(tempCanvasName.c_str(), tempArray, tempArrayRatios, "thick square logY logX");

list->Add((tempCanvas->Clone()));
tempCanvas->Close();
tempArray->Clear();
tempArrayRatios->Clear();


tempCanvasName = "ptDistGeneratedMC_MonashPerugia";
tempArray->Add(getClone("ptDistGeneratedMC", "pp_5TeV_08_pass3_Pythia6", histos));
tempArray->Add(getClone("ptDistGeneratedMC", "pp_5TeV_08_pass3", histos));

tempTitles[0] = "Perugia11";
tempTitles[1] = "Monash13";
tempArray->Add(makeLegend(0.5, 0.35, tempArray, tempTitles));

tempArrayRatios->Add(getRatio("ptDistGeneratedMC", "ptDistGeneratedMC", "Peru/Mona", "pp_5TeV_08_pass3_Pythia6", histos, "pp_5TeV_08_pass3"));
tempArrayRatios->Add(new TF1("line", "1", -0.5, 99.5));

tempCanvas = makeCanvas(tempCanvasName.c_str(), tempArray, tempArrayRatios, "thick square logY logX");

list->Add((tempCanvas->Clone()));
tempCanvas->Close();
tempArray->Clear();
tempArrayRatios->Clear();


tempCanvasName = "moment1Ratio_MonashPerugia";
tempArray->Add(getClone("momentGeneratedMC1", "pp_5TeV_08_pass3_Pythia6", histos));
tempArray->Add(getClone("momentGeneratedMC1", "pp_5TeV_08_pass3", histos));
setRange("X", 0, 60, (TH1D*)tempArray->At(0));

tempTitles[0] = "Perugia11";
tempTitles[1] = "Monash13";
tempArray->Add(makeLegend(0.5, 0.35, tempArray, tempTitles));

tempArrayRatios->Add(getRatio("momentGeneratedMC1", "momentGeneratedMC1", "Peru/Mona", "pp_5TeV_08_pass3_Pythia6", histos, "pp_5TeV_08_pass3"));
tempArrayRatios->Add(new TF1("line", "1", -0.5, 99.5));
setRange("X", 0, 60, (TH1D*)tempArrayRatios->At(0));

tempCanvas = makeCanvas(tempCanvasName.c_str(), tempArray, tempArrayRatios, "thick square");

list->Add((tempCanvas->Clone()));
tempCanvas->Close();
tempArray->Clear();
tempArrayRatios->Clear();

tempCanvasName = "moment2Ratio_MonashPerugia";
tempArray->Add(getClone("momentGeneratedMC2", "pp_5TeV_08_pass3_Pythia6", histos));
tempArray->Add(getClone("momentGeneratedMC2", "pp_5TeV_08_pass3", histos));
setRange("X", 0, 60, (TH1D*)tempArray->At(0));

tempTitles[0] = "Perugia11";
tempTitles[1] = "Monash13";
tempArray->Add(makeLegend(0.5, 0.35, tempArray, tempTitles));

tempArrayRatios->Add(getRatio("momentGeneratedMC2", "momentGeneratedMC2", "Peru/Mona", "pp_5TeV_08_pass3_Pythia6", histos, "pp_5TeV_08_pass3"));
tempArrayRatios->Add(new TF1("line", "1", -0.5, 99.5));
setRange("X", 0, 60, (TH1D*)tempArrayRatios->At(0));

tempCanvas = makeCanvas(tempCanvasName.c_str(), tempArray, tempArrayRatios, "thick square");

list->Add((tempCanvas->Clone()));
tempCanvas->Close();
tempArray->Clear();
tempArrayRatios->Clear();


// _EPOS

tempCanvasName = "multDistMeasuredMC_MonashEPOS";
tempArray->Add(getClone("multDistMeasuredMC", "pp_13TeV_08_EPOS", histos));
tempArray->Add(getClone("multDistMeasuredMC", "pp_13TeV_08", histos));
setRange("X", 0, 60, (TH1D*)tempArray->At(0));

tempTitles[0] = "EPOS";
tempTitles[1] = "Monash13";
tempArray->Add(makeLegend(0.5, 0.35, tempArray, tempTitles));

tempArrayRatios->Add(getRatio("multDistMeasuredMC", "multDistMeasuredMC", "EPOS/Mona", "pp_13TeV_08_EPOS", histos, "pp_13TeV_08"));
tempArrayRatios->Add(new TF1("line", "1", -0.5, 99.5));

setRange("X", 0, 60, (TH1D*)tempArrayRatios->At(0));


tempCanvas = makeCanvas(tempCanvasName.c_str(), tempArray, tempArrayRatios, "thick square logY");

list->Add((tempCanvas->Clone()));
tempCanvas->Close();
tempArray->Clear();
tempArrayRatios->Clear();


tempCanvasName = "multDistGeneratedMC_MonashEPOS";
tempArray->Add(getClone("multDistGeneratedMC", "pp_13TeV_08_EPOS", histos));
tempArray->Add(getClone("multDistGeneratedMC", "pp_13TeV_08", histos));
setRange("X", 0, 80, (TH1D*)tempArray->At(0));

tempTitles[0] = "EPOS";
tempTitles[1] = "Monash13";
tempArray->Add(makeLegend(0.5, 0.35, tempArray, tempTitles));

tempArrayRatios->Add(getRatio("multDistGeneratedMC", "multDistGeneratedMC", "EPOS/Mona", "pp_13TeV_08_EPOS", histos, "pp_13TeV_08"));
tempArrayRatios->Add(new TF1("line", "1", -0.5, 99.5));
setRange("X", 0, 80, (TH1D*)tempArrayRatios->At(0));


tempCanvas = makeCanvas(tempCanvasName.c_str(), tempArray, tempArrayRatios, "thick square logY");

list->Add((tempCanvas->Clone()));
tempCanvas->Close();
tempArray->Clear();
tempArrayRatios->Clear();


tempCanvasName = "ptDistMeasuredMC_MonashEPOS";
tempArray->Add(getClone("ptDistMeasuredMC", "pp_13TeV_08_EPOS", histos));
tempArray->Add(getClone("ptDistMeasuredMC", "pp_13TeV_08", histos));

tempTitles[0] = "EPOS";
tempTitles[1] = "Monash13";
tempArray->Add(makeLegend(0.5, 0.35, tempArray, tempTitles));

tempArrayRatios->Add(getRatio("ptDistMeasuredMC", "ptDistMeasuredMC", "EPOS/Mona", "pp_13TeV_08_EPOS", histos, "pp_13TeV_08"));
tempArrayRatios->Add(new TF1("line", "1", -0.5, 99.5));


tempCanvas = makeCanvas(tempCanvasName.c_str(), tempArray, tempArrayRatios, "thick square logY logX");

list->Add((tempCanvas->Clone()));
tempCanvas->Close();
tempArray->Clear();
tempArrayRatios->Clear();


tempCanvasName = "ptDistGeneratedMC_MonashEPOS";
tempArray->Add(getClone("ptDistGeneratedMC", "pp_13TeV_08_EPOS", histos));
tempArray->Add(getClone("ptDistGeneratedMC", "pp_13TeV_08", histos));

tempTitles[0] = "EPOS";
tempTitles[1] = "Monash13";
tempArray->Add(makeLegend(0.5, 0.35, tempArray, tempTitles));

tempArrayRatios->Add(getRatio("ptDistGeneratedMC", "ptDistGeneratedMC", "EPOS/Mona", "pp_13TeV_08_EPOS", histos, "pp_13TeV_08"));
tempArrayRatios->Add(new TF1("line", "1", -0.5, 99.5));

tempCanvas = makeCanvas(tempCanvasName.c_str(), tempArray, tempArrayRatios, "thick square logY logX");

list->Add((tempCanvas->Clone()));
tempCanvas->Close();
tempArray->Clear();
tempArrayRatios->Clear();


tempCanvasName = "moment1Ratio_MonashEPOS";
tempArray->Add(getClone("momentGeneratedMC1", "pp_13TeV_08_EPOS", histos));
tempArray->Add(getClone("momentGeneratedMC1", "pp_13TeV_08", histos));
setRange("X", 0, 80, (TH1D*)tempArray->At(0));

tempTitles[0] = "EPOS";
tempTitles[1] = "Monash13";
tempArray->Add(makeLegend(0.5, 0.35, tempArray, tempTitles));

tempArrayRatios->Add(getRatio("momentGeneratedMC1", "momentGeneratedMC1", "EPOS/Mona", "pp_13TeV_08_EPOS", histos, "pp_13TeV_08"));
tempArrayRatios->Add(new TF1("line", "1", -0.5, 99.5));
setRange("X", 0, 80, (TH1D*)tempArrayRatios->At(0));

tempCanvas = makeCanvas(tempCanvasName.c_str(), tempArray, tempArrayRatios, "thick square");

list->Add((tempCanvas->Clone()));
tempCanvas->Close();
tempArray->Clear();
tempArrayRatios->Clear();

tempCanvasName = "moment2Ratio_MonashEPOS";
tempArray->Add(getClone("momentGeneratedMC2", "pp_13TeV_08_EPOS", histos));
tempArray->Add(getClone("momentGeneratedMC2", "pp_13TeV_08", histos));
setRange("X", 0, 80, (TH1D*)tempArray->At(0));

tempTitles[0] = "EPOS";
tempTitles[1] = "Monash13";
tempArray->Add(makeLegend(0.5, 0.35, tempArray, tempTitles));

tempArrayRatios->Add(getRatio("momentGeneratedMC2", "momentGeneratedMC2", "EPOS/Mona", "pp_13TeV_08_EPOS", histos, "pp_13TeV_08"));
tempArrayRatios->Add(new TF1("line", "1", -0.5, 99.5));
setRange("X", 0, 80, (TH1D*)tempArrayRatios->At(0));

tempCanvas = makeCanvas(tempCanvasName.c_str(), tempArray, tempArrayRatios, "thick square");

list->Add((tempCanvas->Clone()));
tempCanvas->Close();
tempArray->Clear();
tempArrayRatios->Clear();




// end qm playground
  cout << "-> created " << list->GetEntries() << " plots in " << list->GetName() << endl;
  return list;
}
