/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/
//
// Class for spectrum correction
// Subtraction of hadronic background, Unfolding of the data and
// Renormalization done here
// The following containers have to be set:
//  - Correction framework container for real data
//  - Correction framework container for MC (Efficiency Map)
//  - Correction framework container for background coming from data
//  - Correction framework container for background coming from MC
//
//  Author: 
//            Raphaelle Bailhache <R.Bailhache@gsi.de>
//            Markus Fasel <M.Fasel@gsi.de>
//

#include <TArrayD.h>
#include <TH1.h>
#include <TList.h>
#include <TObjArray.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TMath.h>
#include <TAxis.h>
#include <TGraphErrors.h>
#include <TFile.h>
#include <TPad.h>
#include <TH2D.h>


#include "AliCFContainer.h"
#include "AliCFDataGrid.h"
#include "AliCFEffGrid.h"
#include "AliCFGridSparse.h"
#include "AliCFUnfolding.h"
#include "AliLog.h"

#include "AliHFEspectrum.h"

ClassImp(AliHFEspectrum)

//____________________________________________________________
AliHFEspectrum::AliHFEspectrum(const char *name):
  TNamed(name, ""),
  fCFContainers(NULL),
  fTemporaryObjects(NULL),
  fCorrelation(NULL),
  fBackground(NULL),
  fBackgroundSource(kMCbackground),
  fDumpToFile(kFALSE),
  fNEvents(0),
  fStepMC(-1),
  fStepTrue(-1),
  fStepData(-1),
  fStepGuessedUnfolding(-1),
  fNumberOfIterations(5)
{
  //
  // Default constructor
  //
}

//____________________________________________________________
AliHFEspectrum::~AliHFEspectrum(){
  //
  // Destructor
  //
  if(fCFContainers) delete fCFContainers;
  if(fTemporaryObjects){
    fTemporaryObjects->Clear();
    delete fTemporaryObjects;
  }
}

//____________________________________________________________
void AliHFEspectrum::Correct(AliCFContainer *datacontainer,AliCFContainer *mccontainer,THnSparseF *mccorrelation, AliCFContainer *contaminationcontainer){
  //
  // Correct the spectrum for efficiency and unfolding
  // with both method and compare
  //
  
  gStyle->SetPalette(1);
  gStyle->SetOptStat(1111);
  gStyle->SetPadBorderMode(0);
  gStyle->SetCanvasColor(10);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadRightMargin(0.13);

  SetMCEffStep(8);
  SetMCTruthStep(0);
  SetNumberOfIteration(5);
  SetStepGuessedUnfolding(16);
  SetNumberOfIteration(5);
  SetStepToCorrect(20);
  
  //////////////////
  // Take only pt
  /////////////////
  
  AliCFDataGrid *dataspectrumbeforesubstraction = (AliCFDataGrid *) ((AliCFDataGrid *)GetSpectrum(datacontainer,20))->Clone();
  dataspectrumbeforesubstraction->SetName("dataspectrumbeforesubstraction");

  //
  AliCFDataGrid *dataspectrumaftersubstraction = 0x0;
  SetCorrelation(mccorrelation);
  SetContainer(datacontainer,AliHFEspectrum::kDataContainer);
  SetContainer(mccontainer,AliHFEspectrum::kMCContainer);
  if(contaminationcontainer) {
    SetContainer(contaminationcontainer,AliHFEspectrum::kBackgroundData);
    SetBackgroundSource(AliHFEspectrum::kDataBackground);
    dataspectrumaftersubstraction = SubtractBackground(1,kTRUE);
  }
  
  // Unfold
  
  TList *listunfolded = Unfold(1,dataspectrumaftersubstraction);
  if(!listunfolded){
    printf("Unfolded failed\n");
    return;
  }
  THnSparse *correctedspectrum = (THnSparse *) listunfolded->At(0);
  THnSparse *residualspectrum = (THnSparse *) listunfolded->At(1);
  if(!correctedspectrum){
    AliError("No corrected spectrum\n");
    return;
  }
  if(!residualspectrum){
    AliError("No residul spectrum\n");
    return;
  }   
  
  // Simply correct
  SetMCEffStep(20);  
  AliCFDataGrid *alltogetherCorrection = CorrectForEfficiency(1,dataspectrumaftersubstraction);
  
  //////////
  // Plot
  //////////
  AliCFDataGrid *dataspectrum = 0x0;
  TH1D *measuredTH1Daftersubstraction = 0x0;
  TH1D *measuredTH1Dbeforesubstraction = 0x0;
  TH1D *measuredTH1background = 0x0;
  AliCFDataGrid *contaminationspectrum = 0x0;
  if(dataspectrumaftersubstraction) dataspectrum = dataspectrumaftersubstraction;
  else dataspectrum = dataspectrumbeforesubstraction;

  if(dataspectrumaftersubstraction) {
    
    TCanvas * cbackgroundsubtraction = new TCanvas("backgroundsubtraction","backgroundsubtraction",1000,700);
    if(fBackground) cbackgroundsubtraction->Divide(3,1);
    cbackgroundsubtraction->cd(1);
    measuredTH1Daftersubstraction = dataspectrumaftersubstraction->Project(0);
    measuredTH1Dbeforesubstraction = dataspectrumbeforesubstraction->Project(0);
    CorrectFromTheWidth(measuredTH1Daftersubstraction);
    CorrectFromTheWidth(measuredTH1Dbeforesubstraction);
    measuredTH1Daftersubstraction->SetStats(0);
    measuredTH1Daftersubstraction->SetTitle("");
    measuredTH1Daftersubstraction->GetYaxis()->SetTitle("dN/dp_{T} [(GeV/c)^{-1}]");
    measuredTH1Daftersubstraction->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    measuredTH1Daftersubstraction->SetMarkerStyle(25);
    measuredTH1Daftersubstraction->SetMarkerColor(kBlack);
    measuredTH1Daftersubstraction->SetLineColor(kBlack);
    measuredTH1Dbeforesubstraction->SetStats(0);
    measuredTH1Dbeforesubstraction->SetTitle("");
    measuredTH1Dbeforesubstraction->GetYaxis()->SetTitle("dN/dp_{T} [(GeV/c)^{-1}]");
    measuredTH1Dbeforesubstraction->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    measuredTH1Dbeforesubstraction->SetMarkerStyle(24);
    measuredTH1Dbeforesubstraction->SetMarkerColor(kBlue);
    measuredTH1Dbeforesubstraction->SetLineColor(kBlue);
    measuredTH1Daftersubstraction->Draw();
    measuredTH1Dbeforesubstraction->Draw("same");
    TLegend *legsubstraction = new TLegend(0.4,0.6,0.89,0.89);
    legsubstraction->AddEntry(measuredTH1Dbeforesubstraction,"before substraction","p");
    legsubstraction->AddEntry(measuredTH1Daftersubstraction,"after substraction","p");
    legsubstraction->Draw("same");
    cbackgroundsubtraction->cd(2);
    TH1D* ratiomeasuredcontamination = (TH1D*)measuredTH1Dbeforesubstraction->Clone();
    ratiomeasuredcontamination->SetName("ratiomeasuredcontamination");
    ratiomeasuredcontamination->SetTitle("");
    ratiomeasuredcontamination->GetYaxis()->SetTitle("(with contamination - without contamination) / with contamination");
    ratiomeasuredcontamination->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    ratiomeasuredcontamination->Add(measuredTH1Daftersubstraction,-1.0);
    ratiomeasuredcontamination->Divide(measuredTH1Dbeforesubstraction);
    ratiomeasuredcontamination->SetStats(0);
    ratiomeasuredcontamination->SetMarkerStyle(26);
    ratiomeasuredcontamination->SetMarkerColor(kBlack);
    ratiomeasuredcontamination->SetLineColor(kBlack);
    ratiomeasuredcontamination->Draw();
    if(fBackground) {
      cbackgroundsubtraction->cd(3);
      measuredTH1background = fBackground->Project(0);
      CorrectFromTheWidth(measuredTH1background);
      measuredTH1background->SetStats(0);
      measuredTH1background->SetTitle("");
      measuredTH1background->GetYaxis()->SetTitle("dN/dp_{T} [(GeV/c)^{-1}]");
      measuredTH1background->GetXaxis()->SetTitle("p_{T} [GeV/c]");
      measuredTH1background->SetMarkerStyle(26);
      measuredTH1background->SetMarkerColor(kBlack);
      measuredTH1background->SetLineColor(kBlack);
      measuredTH1background->Draw();
    }
    
    contaminationspectrum = (AliCFDataGrid *) ((AliCFDataGrid *)GetSpectrum(contaminationcontainer,1))->Clone();
    contaminationspectrum->SetName("contaminationspectrum");
    TCanvas * ccontaminationspectrum = new TCanvas("contaminationspectrum","contaminationspectrum",1000,700);
    ccontaminationspectrum->Divide(3,1);
    ccontaminationspectrum->cd(1);
    TH2D * contaminationspectrum2dpteta = contaminationspectrum->Project(1,0);
    TH2D * contaminationspectrum2dptphi = contaminationspectrum->Project(2,0);
    TH2D * contaminationspectrum2detaphi = contaminationspectrum->Project(1,2);
    contaminationspectrum2dpteta->SetStats(0);
    contaminationspectrum2dpteta->SetTitle("");
    contaminationspectrum2dpteta->GetXaxis()->SetTitle("#eta");
    contaminationspectrum2dpteta->GetYaxis()->SetTitle("p_{T} [GeV/c]");
    contaminationspectrum2dptphi->SetStats(0);
    contaminationspectrum2dptphi->SetTitle("");
    contaminationspectrum2dptphi->GetXaxis()->SetTitle("#phi [rad]");
    contaminationspectrum2dptphi->GetYaxis()->SetTitle("p_{T} [GeV/c]");
    contaminationspectrum2detaphi->SetStats(0);
    contaminationspectrum2detaphi->SetTitle("");
    contaminationspectrum2detaphi->GetXaxis()->SetTitle("#eta");
    contaminationspectrum2detaphi->GetYaxis()->SetTitle("#phi [rad]");
    contaminationspectrum2dptphi->Draw("colz");
    ccontaminationspectrum->cd(2);
    contaminationspectrum2dpteta->Draw("colz");
    ccontaminationspectrum->cd(3);
    contaminationspectrum2detaphi->Draw("colz");

  }
  
  
  TCanvas * cresidual = new TCanvas("residual","residual",1000,700);
  cresidual->Divide(2,1);
  cresidual->cd(1);
  gPad->SetLogy();
  TGraphErrors* residualspectrumD = Normalize(residualspectrum);
  if(!residualspectrumD) {
    AliError("Number of Events not set for the normalization");
    return;
  }
  residualspectrumD->SetTitle("");
  residualspectrumD->GetYaxis()->SetTitleOffset(1.5);
  residualspectrumD->GetYaxis()->SetRangeUser(0.000000001,1.0);
  residualspectrumD->SetMarkerStyle(26);
  residualspectrumD->SetMarkerColor(kBlue);
  residualspectrumD->SetLineColor(kBlue);
  residualspectrumD->Draw("AP");
  TGraphErrors* measuredspectrumD = Normalize(dataspectrum);
  measuredspectrumD->SetTitle("");  
  measuredspectrumD->GetYaxis()->SetTitleOffset(1.5);
  measuredspectrumD->GetYaxis()->SetRangeUser(0.000000001,1.0);
  measuredspectrumD->SetMarkerStyle(25);
  measuredspectrumD->SetMarkerColor(kBlack);
  measuredspectrumD->SetLineColor(kBlack);
  measuredspectrumD->Draw("P");
  TLegend *legres = new TLegend(0.4,0.6,0.89,0.89);
  legres->AddEntry(residualspectrumD,"Residual","p");
  legres->AddEntry(measuredspectrumD,"Measured","p");
  legres->Draw("same");
  cresidual->cd(2);
  TH1D *residualTH1D = residualspectrum->Projection(0);
  TH1D *measuredTH1D = dataspectrum->Project(0);
  TH1D* ratioresidual = (TH1D*)residualTH1D->Clone();
  ratioresidual->SetName("ratioresidual");
  ratioresidual->SetTitle("");
  ratioresidual->GetYaxis()->SetTitle("Folded/Measured");
  ratioresidual->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  ratioresidual->Divide(residualTH1D,measuredTH1D,1,1);
  ratioresidual->SetStats(0);
  ratioresidual->Draw();

  //

  
  TCanvas * ccorrected = new TCanvas("corrected","corrected",1000,700);
  ccorrected->Divide(2,1);
  ccorrected->cd(1);
  gPad->SetLogy();
  TGraphErrors* correctedspectrumD = Normalize(correctedspectrum);
  correctedspectrumD->SetTitle("");
  correctedspectrumD->GetYaxis()->SetTitleOffset(1.5);
  correctedspectrumD->GetYaxis()->SetRangeUser(0.000000001,1.0);
  correctedspectrumD->SetMarkerStyle(26);
  correctedspectrumD->SetMarkerColor(kBlue);
  correctedspectrumD->SetLineColor(kBlue);
  correctedspectrumD->Draw("AP");
  TGraphErrors* alltogetherspectrumD = Normalize(alltogetherCorrection);
  alltogetherspectrumD->SetTitle("");
  alltogetherspectrumD->GetYaxis()->SetTitleOffset(1.5);
  alltogetherspectrumD->GetYaxis()->SetRangeUser(0.000000001,1.0);
  alltogetherspectrumD->SetMarkerStyle(25);
  alltogetherspectrumD->SetMarkerColor(kBlack);
  alltogetherspectrumD->SetLineColor(kBlack);
  alltogetherspectrumD->Draw("P");
  TLegend *legcorrected = new TLegend(0.4,0.6,0.89,0.89);
  legcorrected->AddEntry(correctedspectrumD,"Corrected","p");
  legcorrected->AddEntry(alltogetherspectrumD,"Alltogether","p");
  legcorrected->Draw("same");
  ccorrected->cd(2);
  TH1D *correctedTH1D = correctedspectrum->Projection(0);
  TH1D *alltogetherTH1D = alltogetherCorrection->Project(0);
  TH1D* ratiocorrected = (TH1D*)correctedTH1D->Clone();
  ratiocorrected->SetName("ratiocorrected");
  ratiocorrected->SetTitle("");
  ratiocorrected->GetYaxis()->SetTitle("Unfolded/DirectCorrected");
  ratiocorrected->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  ratiocorrected->Divide(correctedTH1D,alltogetherTH1D,1,1);
  ratiocorrected->SetStats(0);
  ratiocorrected->Draw();
  
  // Efficiency correction

  AliCFEffGrid  *efficiencymcPID = (AliCFEffGrid*)  GetEfficiency(mccontainer,8,7);
  AliCFEffGrid  *efficiencymctrackinggeo = (AliCFEffGrid*)  GetEfficiency(mccontainer,7,0);
  AliCFEffGrid  *efficiencymcall = (AliCFEffGrid*)  GetEfficiency(mccontainer,8,0);
  
  AliCFEffGrid  *efficiencyesdall = (AliCFEffGrid*)  GetEfficiency(mccontainer,20,0);
  
  TCanvas * cefficiency = new TCanvas("efficiency","efficiency",1000,700);
  cefficiency->cd(1);
  TH1D* efficiencymcPIDD = efficiencymcPID->Project(0);
  efficiencymcPIDD->SetTitle("");
  efficiencymcPIDD->SetStats(0);
  efficiencymcPIDD->SetMarkerStyle(25);
  efficiencymcPIDD->Draw();
  TH1D* efficiencymctrackinggeoD = efficiencymctrackinggeo->Project(0);
  efficiencymctrackinggeoD->SetTitle("");
  efficiencymctrackinggeoD->SetStats(0);
  efficiencymctrackinggeoD->SetMarkerStyle(26);
  efficiencymctrackinggeoD->Draw("same");
  TH1D* efficiencymcallD = efficiencymcall->Project(0);
  efficiencymcallD->SetTitle("");
  efficiencymcallD->SetStats(0);
  efficiencymcallD->SetMarkerStyle(27);
  efficiencymcallD->Draw("same");
  TH1D* efficiencyesdallD = efficiencyesdall->Project(0);
  efficiencyesdallD->SetTitle("");
  efficiencyesdallD->SetStats(0);
  efficiencyesdallD->SetMarkerStyle(24);
  efficiencyesdallD->Draw("same");
  TLegend *legeff = new TLegend(0.4,0.6,0.89,0.89);
  legeff->AddEntry(efficiencymcPIDD,"PID efficiency","p");
  legeff->AddEntry(efficiencymctrackinggeoD,"Tracking geometry efficiency","p");
  legeff->AddEntry(efficiencymcallD,"Overall efficiency","p");
  legeff->AddEntry(efficiencyesdallD,"Overall efficiency ESD","p");
  legeff->Draw("same");

  // Dump to file if needed

  if(fDumpToFile) {
    
    TFile *out = new TFile("finalSpectrum.root","recreate");
    out->cd();
    //
    residualspectrumD->SetName("UnfoldingResidualSpectrum");
    residualspectrumD->Write();
    measuredspectrumD->SetName("MeasuredSpectrum");
    measuredspectrumD->Write();
    ratioresidual->SetName("RatioResidualSpectrum");
    ratioresidual->Write();
    //
    correctedspectrumD->SetName("UnfoldingCorrectedSpectrum");
    correctedspectrumD->Write();
    alltogetherspectrumD->SetName("AlltogetherSpectrum");
    alltogetherspectrumD->Write();
    ratiocorrected->SetName("RatioUnfoldingAlltogetherSpectrum");
    ratiocorrected->Write();
    //
    correctedspectrum->SetName("UnfoldingCorrectedNotNormalizedSpectrum");
    correctedspectrum->Write();
    alltogetherCorrection->SetName("AlltogetherCorrectedNotNormalizedSpectrum");
    alltogetherCorrection->Write();
    //
    if(measuredTH1Daftersubstraction && measuredTH1Dbeforesubstraction) {
      measuredTH1Daftersubstraction->SetName("Rawptspectrummeasuredaftersubtraction");
      measuredTH1Daftersubstraction->Write();
      measuredTH1Dbeforesubstraction->SetName("Rawptspectrummeasuredbeforesubtraction");
      measuredTH1Dbeforesubstraction->Write();
    }
    if(measuredTH1background) {
      measuredTH1background->SetName("Rawptspectrummeasuredbackground");
      measuredTH1background->Write();
    }
    //
    out->Close(); delete out;

    
  }
}
//____________________________________________________________
AliCFDataGrid* AliHFEspectrum::SubtractBackground(Int_t dimensions, Bool_t setBackground){
  //
  // Apply background subtraction
  //

  AliCFContainer *dataContainer = GetContainer(kDataContainer);
  if(!dataContainer){
    AliError("Data Container not available");
    return NULL;
  }
 
  // Get the data grid in the requested format
  Bool_t initContainer = kFALSE;
  Int_t dims[3];
  switch(dimensions){
    case 1:   dims[0] = 0;
              initContainer = kTRUE;
              break;
    case 3:   for(Int_t i = 0; i < 3; i++) dims[i] = i;
              initContainer = kTRUE;
              break;
    default:
              AliError("Container with this number of dimensions not foreseen (yet)");
  };
  if(!initContainer){
    AliError("Creation of the sliced containers failed. Background estimation not possible");
    return NULL;
  }
  AliCFContainer *spectrumSliced = GetSlicedContainer(dataContainer, dimensions, dims);

  // Make Background Estimate
  AliCFDataGrid *backgroundGrid = NULL;
  if(fBackgroundSource == kDataBackground)
    backgroundGrid = TakeBackgroundFromData(dimensions);
  else
    backgroundGrid = MakeBackgroundEstimateFromMC(dimensions);
  if(!backgroundGrid){
    AliError("Background Estimate not created");
    ClearObject(spectrumSliced);
    return NULL;
  }


  AliCFDataGrid *spectrumSubtracted = new AliCFDataGrid("spectrumSubtracted", "Data Grid for spectrum after Background subtraction", *spectrumSliced,fStepData);
  //spectrumSubtracted->ApplyBGCorrection(*backgroundGrid);
  spectrumSubtracted->Add(backgroundGrid,-1.0);
  if(setBackground){
    if(fBackground) delete fBackground;
    fBackground = backgroundGrid;
  } else delete backgroundGrid;

  return spectrumSubtracted;
}

//____________________________________________________________
TList *AliHFEspectrum::Unfold(Int_t dimensions,AliCFDataGrid* const bgsubpectrum){
  
  //
  // Unfold and eventually correct for efficiency the bgsubspectrum
  //

  AliCFContainer *mcContainer = GetContainer(kMCContainer);
  if(!mcContainer){
    AliError("MC Container not available");
    return NULL;
  }

  if(!fCorrelation){
    AliError("No Correlation map available");
    return NULL;
  }

  // Get the mc grid in the requested format
  Bool_t initContainer = kFALSE;
  Int_t dims[3];
  switch(dimensions){
    case 1:   dims[0] = 0;
              initContainer = kTRUE;
              break;
    case 3:   for(Int_t i = 0; i < 3; i++) dims[i] = i;
              initContainer = kTRUE;
              break;
    default:
              AliError("Container with this number of dimensions not foreseen (yet)");
  };
  if(!initContainer){
    AliError("Creation of the sliced containers failed. Background estimation not possible");
    return NULL;
  }
  AliCFContainer *mcContainerD = GetSlicedContainer(mcContainer, dimensions, dims);
  THnSparse *correlationD = GetSlicedCorrelation(dimensions, dims);

  // Data in the right format
  AliCFDataGrid *dataGrid = 0x0;  
  if(bgsubpectrum) {
    if(bgsubpectrum->GetNVar()!= dimensions) {
      AliError("Not the expected number of dimensions for the AliCFDataGrid\n");
      return NULL;
    }
    dataGrid = bgsubpectrum;
    
  }
  else {

    AliCFContainer *dataContainer = GetContainer(kDataContainer);
    if(!dataContainer){
      AliError("Data Container not available");
      return NULL;
    }

    AliCFContainer *dataContainerD = GetSlicedContainer(dataContainer, dimensions, dims);
    dataGrid = new AliCFDataGrid("dataGrid","dataGrid",*dataContainerD, fStepData);
    
  } 
  
  
  // Guessed
  AliCFDataGrid* guessedGrid = new AliCFDataGrid("guessed","",*mcContainerD, fStepGuessedUnfolding);
  THnSparse* guessedTHnSparse = ((AliCFGridSparse*)guessedGrid->GetData())->GetGrid();

  // Efficiency
  AliCFEffGrid* efficiencyD = new AliCFEffGrid("efficiency","",*mcContainerD);
  efficiencyD->CalculateEfficiency(fStepMC,fStepTrue);

  // Unfold 
  
  AliCFUnfolding unfolding("unfolding","",dimensions,correlationD,efficiencyD->GetGrid(),dataGrid->GetGrid(),guessedTHnSparse);
  unfolding.SetMaxNumberOfIterations(fNumberOfIterations);
  unfolding.UseSmoothing();
  unfolding.Unfold();

  // Results
  THnSparse* result = unfolding.GetUnfolded();
  THnSparse* residual = unfolding.GetEstMeasured();
  TList *listofresults = new TList;
  listofresults->SetOwner();
  listofresults->AddAt((THnSparse*)result->Clone(),0);
  listofresults->AddAt((THnSparse*)residual->Clone(),1);

  return listofresults;

}

//____________________________________________________________
AliCFDataGrid *AliHFEspectrum::CorrectForEfficiency(Int_t dimensions,AliCFDataGrid* const bgsubpectrum){
  
  //
  // Apply unfolding and efficiency correction together to bgsubspectrum
  //

  AliCFContainer *mcContainer = GetContainer(kMCContainer);
  if(!mcContainer){
    AliError("MC Container not available");
    return NULL;
  }

  // Get the data grid in the requested format
  Bool_t initContainer = kFALSE;
  Int_t dims[3];
  switch(dimensions){
    case 1:   dims[0] = 0;
              initContainer = kTRUE;
              break;
    case 3:   for(Int_t i = 0; i < 3; i++) dims[i] = i;
              initContainer = kTRUE;
              break;
    default:
              AliError("Container with this number of dimensions not foreseen (yet)");
  };
  if(!initContainer){
    AliError("Creation of the sliced containers failed. Background estimation not possible");
    return NULL;
  }
  AliCFContainer *mcContainerD = GetSlicedContainer(mcContainer, dimensions, dims);

  // Efficiency
  AliCFEffGrid* efficiencyD = new AliCFEffGrid("efficiency","",*mcContainerD);
  efficiencyD->CalculateEfficiency(fStepMC,fStepTrue);

  // Data in the right format
  AliCFDataGrid *dataGrid = 0x0;  
  if(bgsubpectrum) {
    if(bgsubpectrum->GetNVar()!= dimensions) {
      AliError("Not the expected number of dimensions for the AliCFDataGrid\n");
      return NULL;
    }
    dataGrid = bgsubpectrum;
    
  }
  else {

    AliCFContainer *dataContainer = GetContainer(kDataContainer);
    if(!dataContainer){
      AliError("Data Container not available");
      return NULL;
    }

    AliCFContainer *dataContainerD = GetSlicedContainer(dataContainer, dimensions, dims);
    dataGrid = new AliCFDataGrid("dataGrid","dataGrid",*dataContainerD, fStepData);
    
  } 

  // Correct
  AliCFDataGrid *result = (AliCFDataGrid *) dataGrid->Clone();
  result->ApplyEffCorrection(*efficiencyD);
  

  return result;

}
//__________________________________________________________________________________
TGraphErrors *AliHFEspectrum::Normalize(THnSparse * const spectrum) const {
  //
  // Normalize the spectrum to 1/(2*Pi*p_{T})*dN/dp_{T} (GeV/c)^{-2}
  // Give the final pt spectrum to be compared
  //
 
  if(fNEvents > 0) {
    
    TH1D* projection = spectrum->Projection(0);
    CorrectFromTheWidth(projection);
    TGraphErrors *graphError = NormalizeTH1(projection);
    return graphError;
  
  }
    
  return 0x0;
  

}
//__________________________________________________________________________________
TGraphErrors *AliHFEspectrum::Normalize(AliCFDataGrid * const spectrum) const {
  //
  // Normalize the spectrum to 1/(2*Pi*p_{T})*dN/dp_{T} (GeV/c)^{-2}
  // Give the final pt spectrum to be compared
  //
  if(fNEvents > 0) {

    TH1D* projection = spectrum->Project(0);
    CorrectFromTheWidth(projection);
    TGraphErrors *graphError = NormalizeTH1(projection);
    return graphError;
    
  }

  return 0x0;
  

}
//__________________________________________________________________________________
TGraphErrors *AliHFEspectrum::NormalizeTH1(TH1 *input) const {
  //
  // Normalize the spectrum to 1/(2*Pi*p_{T})*dN/dp_{T} (GeV/c)^{-2}
  // Give the final pt spectrum to be compared
  //
  if(fNEvents > 0) {

    TGraphErrors *spectrumNormalized = new TGraphErrors(input->GetNbinsX());
    Double_t p = 0, dp = 0; Int_t point = 1;
    Double_t n = 0, dN = 0;
    Double_t nCorr = 0, dNcorr = 0;
    Double_t errdN = 0, errdp = 0;
    for(Int_t ibin = input->GetXaxis()->GetFirst(); ibin <= input->GetXaxis()->GetLast(); ibin++){
      point = ibin - input->GetXaxis()->GetFirst();
      p = input->GetXaxis()->GetBinCenter(ibin);
      dp = input->GetXaxis()->GetBinWidth(ibin)/2.;
      n = input->GetBinContent(ibin);
      dN = input->GetBinError(ibin);

      // New point
      nCorr = 0.5 * 1/1.6 * 1./(Double_t)(fNEvents) * 1/(2. * TMath::Pi() * p) * n;
      errdN = 1./(2. * TMath::Pi() * p);
      errdp = 1./(2. * TMath::Pi() * p*p) * n;
      dNcorr = 0.5 * 1/1.6 * 1./(Double_t)(fNEvents) * TMath::Sqrt(errdN * errdN * dN *dN + errdp *errdp * dp *dp);
      
      spectrumNormalized->SetPoint(point, p, nCorr);
      spectrumNormalized->SetPointError(point, dp, dNcorr);
    }
    spectrumNormalized->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    spectrumNormalized->GetYaxis()->SetTitle("#frac{1}{2 #pi p_{T}} #frac{dN}{dp_{T}} / [GeV/c]^{-2}");
    spectrumNormalized->SetMarkerStyle(22);
    spectrumNormalized->SetMarkerColor(kBlue);
    spectrumNormalized->SetLineColor(kBlue);
    
    return spectrumNormalized;
    
  }

  return 0x0;
  

}
//____________________________________________________________
void AliHFEspectrum::SetContainer(AliCFContainer *cont, AliHFEspectrum::CFContainer_t type){
  //
  // Set the container for a given step to the 
  //
  if(!fCFContainers) fCFContainers = new TList;
  fCFContainers->AddAt(cont, type);
}

//____________________________________________________________
AliCFContainer *AliHFEspectrum::GetContainer(AliHFEspectrum::CFContainer_t type){
  //
  // Get Correction Framework Container for given type
  //
  if(!fCFContainers) return NULL;
  return dynamic_cast<AliCFContainer *>(fCFContainers->At(type));
}

//____________________________________________________________
AliCFDataGrid *AliHFEspectrum::TakeBackgroundFromData(Int_t nDim) {
  // 
  // Take Background Estimate from Data
  //
 
  AliCFContainer *backgroundContainer = GetContainer(kBackgroundData);
  if(!backgroundContainer){
    AliError("MC background container not found");
    return NULL;
  }
  AliInfo(Form("Making background Estimate from Data in %d Dimensions", nDim));

  Bool_t initContainer = kFALSE;
  Int_t dimensions[3];
  switch(nDim){
    case 1:   dimensions[0] = 0;
              initContainer = kTRUE;
              break;
    case 3:   for(Int_t i = 0; i < 3; i++) dimensions[i] = i;
              initContainer = kTRUE;
              break;
    default:
              AliError("Container with this number of dimensions not foreseen (yet)");
  };
  if(!initContainer){
    AliError("Creation of the sliced containers failed. Background estimation not possible");
    return NULL;
  }
  AliCFContainer *slicedBackground = GetSlicedContainer(backgroundContainer, nDim, dimensions);
  AliCFDataGrid *contaminationspectrum = new AliCFDataGrid("ContaminationGrid","ContaminationGrid",*slicedBackground,1);
  
  return contaminationspectrum;
}

//____________________________________________________________
AliCFDataGrid *AliHFEspectrum::MakeBackgroundEstimateFromMC(Int_t nDim){
  //
  // Make Background Estimate using MC
  // Calculate ratio of hadronic background from MC and 
  // apply this on data
  //
  AliCFContainer *backgroundContainer = GetContainer(kBackgroundMC);
  AliCFContainer *dataContainer = GetContainer(kDataContainer);
  if(!backgroundContainer){
    AliError("MC background container not found");
    return NULL;
  }
  if(!dataContainer){
    AliError("Data container not found");
    return NULL;
  }
  AliInfo(Form("Making background Estimate from MC in %d Dimensions", nDim));

  Bool_t initContainer = kFALSE;
  Int_t dimensions[3];
  switch(nDim){
    case 1:   dimensions[0] = 1;
              initContainer = kTRUE;
              break;
    case 3:   for(Int_t i = 0; i < 3; i++) dimensions[i] = i;
              initContainer = kTRUE;
              break;
    default:
              AliError("Container with this number of dimensions not foreseen (yet)");
  };
  if(!initContainer){
    AliError("Creation of the sliced containers failed. Background estimation not possible");
    return NULL;
  }
  AliCFContainer *slicedBackground = GetSlicedContainer(backgroundContainer, nDim, dimensions);
  AliCFContainer *slicedData = GetSlicedContainer(dataContainer, nDim, dimensions);
  
  // Create Efficiency Grid and data grid
  AliCFEffGrid backgroundRatio("backgroundRatio", "BackgroundRatio", *slicedBackground);
  backgroundRatio.CalculateEfficiency(1, 0);
  AliCFDataGrid *backgroundEstimate = new AliCFDataGrid("backgroundEstimate", "Grid for Background Estimate", *slicedData, fStepData);
  backgroundEstimate->ApplyEffCorrection(backgroundRatio);

  return backgroundEstimate;
}

//____________________________________________________________
AliCFContainer *AliHFEspectrum::GetSlicedContainer(AliCFContainer *container, Int_t nDim, Int_t *dimensions) {
  //
  // Slice Pt bin
  //
  
  Double_t *varMin = new Double_t[container->GetNVar()],
           *varMax = new Double_t[container->GetNVar()];

  Double_t *binLimits;
  for(Int_t ivar = 0; ivar < container->GetNVar(); ivar++){
    
    binLimits = new Double_t[container->GetNBins(ivar)+1];
    container->GetBinLimits(ivar,binLimits);
    varMin[ivar] = binLimits[0];
    varMax[ivar] = binLimits[container->GetNBins(ivar)];
    delete[] binLimits;
  }
  
  AliCFContainer *k = container->MakeSlice(nDim, dimensions, varMin, varMax);
  AddTemporaryObject(k);
  delete[] varMin; delete[] varMax;

  return k;

}

//_________________________________________________________________________
THnSparse *AliHFEspectrum::GetSlicedCorrelation(Int_t nDim, Int_t *dimensions) const {
  //
  // Slice Pt correlation
  //

  Int_t ndimensions = fCorrelation->GetNdimensions();
  printf("Number of dimension %d correlation map\n",ndimensions);
  if(ndimensions < (2*nDim)) {
    AliError("Problem in the dimensions");
    return NULL;
  }
  Int_t ndimensionsContainer = (Int_t) ndimensions/2;
  printf("Number of dimension %d container\n",ndimensionsContainer);

  Int_t *dim = new Int_t[nDim*2];
  for(Int_t iter=0; iter < nDim; iter++){
    dim[iter] = dimensions[iter];
    dim[iter+nDim] = ndimensionsContainer + dimensions[iter];
    printf("For iter %d: %d and iter+nDim %d: %d\n",iter,dimensions[iter],iter+nDim,ndimensionsContainer + dimensions[iter]);
  }
    
  THnSparse *k = fCorrelation->Projection(nDim*2,dim);

  delete[] dim; 
  return k;
  
}
//___________________________________________________________________________
void AliHFEspectrum::CorrectFromTheWidth(TH1D *h1) const {
  //
  // Correct from the width of the bins --> dN/dp_{T} (GeV/c)^{-1}
  //

  TAxis *axis = h1->GetXaxis();
  Int_t nbinX = h1->GetNbinsX();

  for(Int_t i = 1; i <= nbinX; i++) {

    Double_t width = axis->GetBinWidth(i);
    Double_t content = h1->GetBinContent(i);
    Double_t error = h1->GetBinError(i); 
    h1->SetBinContent(i,content/width); 
    h1->SetBinError(i,error/width);
  }

}
//____________________________________________________________
void AliHFEspectrum::AddTemporaryObject(TObject *o){
  // 
  // Emulate garbage collection on class level
  // 
  if(!fTemporaryObjects) {
    fTemporaryObjects= new TList;
    fTemporaryObjects->SetOwner();
  }
  fTemporaryObjects->Add(o);
}

//____________________________________________________________
void AliHFEspectrum::ClearObject(TObject *o){
  //
  // Do a safe deletion
  //
  if(fTemporaryObjects){
    if(fTemporaryObjects->FindObject(o)) fTemporaryObjects->Remove(o);
    delete o;
  } else{
    // just delete
    delete o;
  }
}
//_________________________________________________________________________
TObject* AliHFEspectrum::GetSpectrum(AliCFContainer * const c, Int_t step) {
  AliCFDataGrid* data = new AliCFDataGrid("data","",*c, step);
  return data;
}
//_________________________________________________________________________
TObject* AliHFEspectrum::GetEfficiency(AliCFContainer * const c, Int_t step, Int_t step0) {
  // 
  // Create efficiency grid and calculate efficiency
  // of step to step0
  //
  TString name("eff");
  name += step;
  name+= step0;
  AliCFEffGrid* eff = new AliCFEffGrid((const char*)name,"",*c);
  eff->CalculateEfficiency(step,step0);
  return eff;
}
