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

/* $Id$ */

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

#include "AliPID.h"
#include "AliCFContainer.h"
#include "AliCFDataGrid.h"
#include "AliCFEffGrid.h"
#include "AliCFGridSparse.h"
#include "AliCFUnfolding.h"
#include "AliLog.h"

#include "AliHFEspectrum.h"
#include "AliHFEcuts.h"
#include "AliHFEcontainer.h"

ClassImp(AliHFEspectrum)

//____________________________________________________________
AliHFEspectrum::AliHFEspectrum(const char *name):
  TNamed(name, ""),
  fCFContainers(NULL),
  fTemporaryObjects(NULL),
  fCorrelation(NULL),
  fBackground(NULL),
  fInclusiveSpectrum(kTRUE),
  fDumpToFile(kFALSE),
  fNbDimensions(1),
  fNEvents(0),
  fStepMC(-1),
  fStepTrue(-1),
  fStepData(-1),
  fStepBeforeCutsV0(-1),
  fStepAfterCutsV0(-1),
  fStepGuessedUnfolding(-1),
  fNumberOfIterations(5),
  fDebugLevel(0)
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
Bool_t AliHFEspectrum::Init(const AliHFEcontainer *datahfecontainer, const AliHFEcontainer *mchfecontainer, const AliHFEcontainer *v0hfecontainer){
  //
  // Init what we need for the correction:
  //
  // Raw spectrum, hadron contamination
  // MC efficiency maps, correlation matrix
  // V0 efficiency if wanted
  //
  // This for a given dimension.
  // If no inclusive spectrum, then take only efficiency map for beauty electron
  // and the appropriate correlation matrix
  //

  // Get the requested format
  Int_t dims[3];
  switch(fNbDimensions){
  case 1:   dims[0] = 0;
    break;
  case 2:   for(Int_t i = 0; i < 2; i++) dims[i] = i;
    break;
  case 3:   for(Int_t i = 0; i < 3; i++) dims[i] = i;
    break;
  default:
    AliError("Container with this number of dimensions not foreseen (yet)");
    return kFALSE;
  };
  
  // Data container: raw spectrum + hadron contamination  
  AliCFContainer *datacontainer = datahfecontainer->GetCFContainer("recTrackContReco");
  AliCFContainer *contaminationcontainer = datahfecontainer->GetCFContainer("hadronicBackground");
  if((!datacontainer) || (!contaminationcontainer)) return kFALSE;

  AliCFContainer *datacontainerD = GetSlicedContainer(datacontainer, fNbDimensions, dims);
  AliCFContainer *contaminationcontainerD = GetSlicedContainer(contaminationcontainer, fNbDimensions, dims);
  if((!datacontainerD) || (!contaminationcontainerD)) return kFALSE;
  SetContainer(datacontainerD,AliHFEspectrum::kDataContainer);
  SetContainer(contaminationcontainerD,AliHFEspectrum::kBackgroundData);

  // MC container: ESD/MC efficiency maps + MC/MC efficiency maps 
  AliCFContainer *mccontaineresd = 0x0;
  AliCFContainer *mccontainermc = 0x0;
  if(fInclusiveSpectrum) {
    mccontaineresd = mchfecontainer->MakeMergedCFContainer("sumesd","sumesd","MCTrackCont:recTrackContReco");
    mccontainermc = mchfecontainer->MakeMergedCFContainer("summc","summc","MCTrackCont:recTrackContMC");
  }
  else {
    mccontaineresd = mchfecontainer->MakeMergedCFContainer("sumesd","sumesd","MCTrackCont:recTrackContReco:recTrackContDEReco");
    mccontainermc = mchfecontainer->MakeMergedCFContainer("summc","summc","MCTrackCont:recTrackContMC:recTrackContDEMC");
  }
  if((!mccontaineresd) || (!mccontainermc)) return kFALSE;  

  Int_t source = -1;
  if(!fInclusiveSpectrum) source = 1;
  AliCFContainer *mccontaineresdD = GetSlicedContainer(mccontaineresd, fNbDimensions, dims, source);
  AliCFContainer *mccontainermcD = GetSlicedContainer(mccontainermc, fNbDimensions, dims, source);
  if((!mccontaineresdD) || (!mccontainermcD)) return kFALSE;
  SetContainer(mccontainermcD,AliHFEspectrum::kMCContainerMC);
  SetContainer(mccontaineresdD,AliHFEspectrum::kMCContainerESD);

  // MC container: correlation matrix
  THnSparseF *mccorrelation = 0x0;
  if(fInclusiveSpectrum) {
    if((fStepMC==(AliHFEcuts::kNcutStepsMCTrack + AliHFEcuts::kStepHFEcutsTRD + 2))) mccorrelation = mchfecontainer->GetCorrelationMatrix("correlationstepafterPID");
    if((fStepMC==(AliHFEcuts::kNcutStepsMCTrack + AliHFEcuts::kStepHFEcutsTRD + 1))) mccorrelation = mchfecontainer->GetCorrelationMatrix("correlationstepafterTOF");
    if((fStepMC==(AliHFEcuts::kNcutStepsMCTrack + AliHFEcuts::kStepHFEcutsTRD))) mccorrelation = mchfecontainer->GetCorrelationMatrix("correlationstepbeforePID");
  }
  else mccorrelation = mchfecontainer->GetCorrelationMatrix("correlationstepafterDE");
  if(!mccorrelation) return kFALSE;
  THnSparseF *mccorrelationD = GetSlicedCorrelation(mccorrelation, fNbDimensions, dims);
  if(!mccorrelationD) {
    printf("No correlation\n");
    return kFALSE;
  }
  SetCorrelation(mccorrelationD);

  // V0 container Electron, pt eta phi
  if(v0hfecontainer) {
    AliCFContainer *containerV0 = v0hfecontainer->GetCFContainer("taggedTrackContainerReco");
    if(!containerV0) return kFALSE;
    AliCFContainer *containerV0Electron = GetSlicedContainer(containerV0, fNbDimensions, dims, AliPID::kElectron+1);
    if(!containerV0Electron) return kFALSE;
    SetContainer(containerV0Electron,AliHFEspectrum::kDataContainerV0);
  } 

  if(fDebugLevel>0){
   
    AliCFDataGrid *contaminationspectrum = (AliCFDataGrid *) ((AliCFDataGrid *)GetSpectrum(contaminationcontainer,1))->Clone();
    contaminationspectrum->SetName("contaminationspectrum");
    TCanvas * ccontaminationspectrum = new TCanvas("contaminationspectrum","contaminationspectrum",1000,700);
    ccontaminationspectrum->Divide(3,1);
    ccontaminationspectrum->cd(1);
    TH2D * contaminationspectrum2dpteta = (TH2D *) contaminationspectrum->Project(1,0);
    TH2D * contaminationspectrum2dptphi = (TH2D *) contaminationspectrum->Project(2,0);
    TH2D * contaminationspectrum2detaphi = (TH2D *) contaminationspectrum->Project(1,2);
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


  return kTRUE;
 
  
}

//____________________________________________________________
Bool_t AliHFEspectrum::Correct(Bool_t subtractcontamination){
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

  ///////////////////////////
  // Check initialization
  ///////////////////////////

  if((!GetContainer(kDataContainer)) || (!GetContainer(kMCContainerMC)) || (!GetContainer(kMCContainerESD))){
    AliInfo("You have to init before");
    return kFALSE;
  }
  
  if((fStepTrue == 0) && (fStepMC == 0) && (fStepData == 0)) {
    AliInfo("You have to set the steps before: SetMCTruthStep, SetMCEffStep, SetStepToCorrect");
    return kFALSE;
  }
 
  SetNumberOfIteration(50);
  SetStepGuessedUnfolding(AliHFEcuts::kStepRecKineITSTPC + AliHFEcuts::kNcutStepsMCTrack);
    
  AliCFDataGrid *dataGridAfterFirstSteps = 0x0;
  //////////////////////////////////
  // Subtract hadron background
  /////////////////////////////////
  AliCFDataGrid *dataspectrumaftersubstraction = 0x0;
  if(subtractcontamination) {
    dataspectrumaftersubstraction = SubtractBackground(kTRUE);
    dataGridAfterFirstSteps = dataspectrumaftersubstraction;
  }

  ////////////////////////////////////////////////
  // Correct for TPC efficiency from V0
  ///////////////////////////////////////////////
  AliCFDataGrid *dataspectrumafterV0efficiencycorrection = 0x0;
  AliCFContainer *dataContainerV0 = GetContainer(kDataContainerV0);
  if(dataContainerV0){
    dataspectrumafterV0efficiencycorrection = CorrectV0Efficiency(dataspectrumaftersubstraction);
    dataGridAfterFirstSteps = dataspectrumafterV0efficiencycorrection;  
  }
  
  ///////////////
  // Unfold
  //////////////
  TList *listunfolded = Unfold(dataGridAfterFirstSteps);
  if(!listunfolded){
    printf("Unfolded failed\n");
    return kFALSE;
  }
  THnSparse *correctedspectrum = (THnSparse *) listunfolded->At(0);
  THnSparse *residualspectrum = (THnSparse *) listunfolded->At(1);
  if(!correctedspectrum){
    AliError("No corrected spectrum\n");
    return kFALSE;
  }
  if(!residualspectrum){
    AliError("No residul spectrum\n");
    return kFALSE;
  }   
  
  /////////////////////
  // Simply correct
  ////////////////////
  AliCFDataGrid *alltogetherCorrection = CorrectForEfficiency(dataGridAfterFirstSteps);
  

  //////////
  // Plot
  //////////
  if(fDebugLevel > 0.0) {
  
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
    TH1D *alltogetherTH1D = (TH1D *) alltogetherCorrection->Project(0);
    TH1D* ratiocorrected = (TH1D*)correctedTH1D->Clone();
    ratiocorrected->SetName("ratiocorrected");
    ratiocorrected->SetTitle("");
    ratiocorrected->GetYaxis()->SetTitle("Unfolded/DirectCorrected");
    ratiocorrected->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    ratiocorrected->Divide(correctedTH1D,alltogetherTH1D,1,1);
    ratiocorrected->SetStats(0);
    ratiocorrected->Draw();


    // Dump to file if needed

    if(fDumpToFile) {
      
      TFile *out = new TFile("finalSpectrum.root","recreate");
      out->cd();
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
      out->Close(); delete out;
    }
  

  }


  

  return kTRUE;
}
//____________________________________________________________
AliCFDataGrid* AliHFEspectrum::SubtractBackground(Bool_t setBackground){
  //
  // Apply background subtraction
  //
  
  // Raw spectrum
  AliCFContainer *dataContainer = GetContainer(kDataContainer);
  if(!dataContainer){
    AliError("Data Container not available");
    return NULL;
  }
  AliCFDataGrid *spectrumSubtracted = new AliCFDataGrid("spectrumSubtracted", "Data Grid for spectrum after Background subtraction", *dataContainer,fStepData);

  AliCFDataGrid *dataspectrumbeforesubstraction = (AliCFDataGrid *) ((AliCFDataGrid *)GetSpectrum(GetContainer(kDataContainer),fStepData))->Clone();
  dataspectrumbeforesubstraction->SetName("dataspectrumbeforesubstraction"); 
 
  // Background Estimate
  AliCFContainer *backgroundContainer = GetContainer(kBackgroundData);
  if(!backgroundContainer){
    AliError("MC background container not found");
    return NULL;
  }
 
  Int_t stepbackground = 1;
  if(!fInclusiveSpectrum) stepbackground = 2;
  AliCFDataGrid *backgroundGrid = new AliCFDataGrid("ContaminationGrid","ContaminationGrid",*backgroundContainer,stepbackground);

  // Subtract
  spectrumSubtracted->Add(backgroundGrid,-1.0);
  if(setBackground){
    if(fBackground) delete fBackground;
    fBackground = backgroundGrid;
  } else delete backgroundGrid;


  if(fDebugLevel > 0) {
    
    TCanvas * cbackgroundsubtraction = new TCanvas("backgroundsubtraction","backgroundsubtraction",1000,700);
    cbackgroundsubtraction->Divide(3,1);
    cbackgroundsubtraction->cd(1);
    gPad->SetLogy();
    TH1D *measuredTH1Daftersubstraction = (TH1D *) spectrumSubtracted->Project(0);
    TH1D *measuredTH1Dbeforesubstraction = (TH1D *) dataspectrumbeforesubstraction->Project(0);
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
    legsubstraction->AddEntry(measuredTH1Dbeforesubstraction,"With hadron contamination","p");
    legsubstraction->AddEntry(measuredTH1Daftersubstraction,"Without hadron contamination ","p");
    legsubstraction->Draw("same");
    cbackgroundsubtraction->cd(2);
    gPad->SetLogy();
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
    cbackgroundsubtraction->cd(3);
    TH1D *measuredTH1background = (TH1D *) backgroundGrid->Project(0);
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
  

  return spectrumSubtracted;
}
//____________________________________________________________
AliCFDataGrid *AliHFEspectrum::CorrectV0Efficiency(AliCFDataGrid* const bgsubpectrum){
  
  //
  // Apply TPC pid efficiency correction from V0
  //

  AliCFContainer *v0Container = GetContainer(kDataContainerV0);
  if(!v0Container){
    AliError("V0 Container not available");
    return NULL;
  }

  // Efficiency
  AliCFEffGrid* efficiencyD = new AliCFEffGrid("efficiency","",*v0Container);
  efficiencyD->CalculateEfficiency(fStepAfterCutsV0,fStepBeforeCutsV0);

  // Data in the right format
  AliCFDataGrid *dataGrid = 0x0;  
  if(bgsubpectrum) {
    dataGrid = bgsubpectrum;
  }
  else {
    
    AliCFContainer *dataContainer = GetContainer(kDataContainer);
    if(!dataContainer){
      AliError("Data Container not available");
      return NULL;
    }

    dataGrid = new AliCFDataGrid("dataGrid","dataGrid",*dataContainer, fStepData);
  } 

  // Correct
  AliCFDataGrid *result = (AliCFDataGrid *) dataGrid->Clone();
  result->ApplyEffCorrection(*efficiencyD);

  if(fDebugLevel > 0) {
    
    TCanvas * cV0Efficiency = new TCanvas("V0Efficiency","V0Efficiency",1000,700);
    cV0Efficiency->Divide(2,1);
    cV0Efficiency->cd(1);
    TH1D *afterE = (TH1D *) result->Project(0);
    TH1D *beforeE = (TH1D *) dataGrid->Project(0);
    afterE->SetStats(0);
    afterE->SetTitle("");
    afterE->GetYaxis()->SetTitle("dN/dp_{T} [(GeV/c)^{-1}]");
    afterE->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    afterE->SetMarkerStyle(25);
    afterE->SetMarkerColor(kBlack);
    afterE->SetLineColor(kBlack);
    beforeE->SetStats(0);
    beforeE->SetTitle("");
    beforeE->GetYaxis()->SetTitle("dN/dp_{T} [(GeV/c)^{-1}]");
    beforeE->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    beforeE->SetMarkerStyle(24);
    beforeE->SetMarkerColor(kBlue);
    beforeE->SetLineColor(kBlue);
    afterE->Draw();
    beforeE->Draw("same");
    TLegend *legV0efficiency = new TLegend(0.4,0.6,0.89,0.89);
    legV0efficiency->AddEntry(beforeE,"Before Efficiency correction","p");
    legV0efficiency->AddEntry(afterE,"After Efficiency correction","p");
    legV0efficiency->Draw("same");
    cV0Efficiency->cd(2);
    TH1D* efficiencyDproj = (TH1D *) efficiencyD->Project(0);
    efficiencyDproj->SetTitle("");
    efficiencyDproj->SetStats(0);
    efficiencyDproj->SetMarkerStyle(25);
    efficiencyDproj->Draw();


  }

  
  return result;

}
//____________________________________________________________
TList *AliHFEspectrum::Unfold(AliCFDataGrid* const bgsubpectrum){
  
  //
  // Unfold and eventually correct for efficiency the bgsubspectrum
  //

  AliCFContainer *mcContainer = GetContainer(kMCContainerMC);
  if(!mcContainer){
    AliError("MC Container not available");
    return NULL;
  }

  if(!fCorrelation){
    AliError("No Correlation map available");
    return NULL;
  }

  // Data 
  AliCFDataGrid *dataGrid = 0x0;  
  if(bgsubpectrum) {
    dataGrid = bgsubpectrum;
  }
  else {

    AliCFContainer *dataContainer = GetContainer(kDataContainer);
    if(!dataContainer){
      AliError("Data Container not available");
      return NULL;
    }

    dataGrid = new AliCFDataGrid("dataGrid","dataGrid",*dataContainer, fStepData);
  } 
  
  // Guessed
  AliCFDataGrid* guessedGrid = new AliCFDataGrid("guessed","",*mcContainer, fStepGuessedUnfolding);
  THnSparse* guessedTHnSparse = ((AliCFGridSparse*)guessedGrid->GetData())->GetGrid();

  // Efficiency
  AliCFEffGrid* efficiencyD = new AliCFEffGrid("efficiency","",*mcContainer);
  efficiencyD->CalculateEfficiency(fStepMC,fStepTrue);

  // Unfold 
  
  AliCFUnfolding unfolding("unfolding","",fNbDimensions,fCorrelation,efficiencyD->GetGrid(),dataGrid->GetGrid(),guessedTHnSparse);
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

  if(fDebugLevel > 0) {
    
    TCanvas * cresidual = new TCanvas("residual","residual",1000,700);
    cresidual->Divide(2,1);
    cresidual->cd(1);
    gPad->SetLogy();
    TGraphErrors* residualspectrumD = Normalize(residual);
    if(!residualspectrumD) {
      AliError("Number of Events not set for the normalization");
      return kFALSE;
    }
    residualspectrumD->SetTitle("");
    residualspectrumD->GetYaxis()->SetTitleOffset(1.5);
    residualspectrumD->GetYaxis()->SetRangeUser(0.000000001,1.0);
    residualspectrumD->SetMarkerStyle(26);
    residualspectrumD->SetMarkerColor(kBlue);
    residualspectrumD->SetLineColor(kBlue);
    residualspectrumD->Draw("AP");
    AliCFDataGrid *dataGridBis = (AliCFDataGrid *) (dataGrid->Clone());
    dataGridBis->SetName("dataGridBis"); 
    TGraphErrors* measuredspectrumD = Normalize(dataGridBis);
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
    TH1D *residualTH1D = residual->Projection(0);
    TH1D *measuredTH1D = (TH1D *) dataGridBis->Project(0);
    TH1D* ratioresidual = (TH1D*)residualTH1D->Clone();
    ratioresidual->SetName("ratioresidual");
    ratioresidual->SetTitle("");
    ratioresidual->GetYaxis()->SetTitle("Folded/Measured");
    ratioresidual->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    ratioresidual->Divide(residualTH1D,measuredTH1D,1,1);
    ratioresidual->SetStats(0);
    ratioresidual->Draw();
    
  }

  return listofresults;

}

//____________________________________________________________
AliCFDataGrid *AliHFEspectrum::CorrectForEfficiency(AliCFDataGrid* const bgsubpectrum){
  
  //
  // Apply unfolding and efficiency correction together to bgsubspectrum
  //

  AliCFContainer *mcContainer = GetContainer(kMCContainerESD);
  if(!mcContainer){
    AliError("MC Container not available");
    return NULL;
  }

  // Efficiency
  AliCFEffGrid* efficiencyD = new AliCFEffGrid("efficiency","",*mcContainer);
  efficiencyD->CalculateEfficiency(fStepMC,fStepTrue);

  // Data in the right format
  AliCFDataGrid *dataGrid = 0x0;  
  if(bgsubpectrum) {
    dataGrid = bgsubpectrum;
  }
  else {
    
    AliCFContainer *dataContainer = GetContainer(kDataContainer);
    if(!dataContainer){
      AliError("Data Container not available");
      return NULL;
    }

    dataGrid = new AliCFDataGrid("dataGrid","dataGrid",*dataContainer, fStepData);
  } 

  // Correct
  AliCFDataGrid *result = (AliCFDataGrid *) dataGrid->Clone();
  result->ApplyEffCorrection(*efficiencyD);

  if(fDebugLevel > 0) {
    
    printf("Step MC: %d\n",fStepMC);
    printf("Step tracking: %d\n",AliHFEcuts::kStepHFEcutsTRD + AliHFEcuts::kNcutStepsMCTrack);
    printf("Step MC true: %d\n",fStepTrue);
    AliCFEffGrid  *efficiencymcPID = (AliCFEffGrid*)  GetEfficiency(GetContainer(kMCContainerMC),fStepMC,AliHFEcuts::kStepHFEcutsTRD + AliHFEcuts::kNcutStepsMCTrack);
    AliCFEffGrid  *efficiencymctrackinggeo = (AliCFEffGrid*)  GetEfficiency(GetContainer(kMCContainerMC),AliHFEcuts::kStepHFEcutsTRD + AliHFEcuts::kNcutStepsMCTrack,fStepTrue);
    AliCFEffGrid  *efficiencymcall = (AliCFEffGrid*)  GetEfficiency(GetContainer(kMCContainerMC),fStepMC,fStepTrue);
    
    AliCFEffGrid  *efficiencyesdall = (AliCFEffGrid*)  GetEfficiency(GetContainer(kMCContainerESD),fStepMC,fStepTrue);
    
    TCanvas * cefficiency = new TCanvas("efficiency","efficiency",1000,700);
    cefficiency->cd(1);
    TH1D* efficiencymcPIDD = (TH1D *) efficiencymcPID->Project(0);
    efficiencymcPIDD->SetTitle("");
    efficiencymcPIDD->SetStats(0);
    efficiencymcPIDD->SetMarkerStyle(25);
    efficiencymcPIDD->Draw();
    TH1D* efficiencymctrackinggeoD = (TH1D *) efficiencymctrackinggeo->Project(0);
    efficiencymctrackinggeoD->SetTitle("");
    efficiencymctrackinggeoD->SetStats(0);
    efficiencymctrackinggeoD->SetMarkerStyle(26);
    efficiencymctrackinggeoD->Draw("same");
    TH1D* efficiencymcallD = (TH1D *) efficiencymcall->Project(0);
    efficiencymcallD->SetTitle("");
    efficiencymcallD->SetStats(0);
    efficiencymcallD->SetMarkerStyle(27);
    efficiencymcallD->Draw("same");
    TH1D* efficiencyesdallD = (TH1D *) efficiencyesdall->Project(0);
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
  }
  
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

    TH1D* projection = (TH1D *) spectrum->Project(0);
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
AliCFContainer *AliHFEspectrum::GetSlicedContainer(AliCFContainer *container, Int_t nDim, Int_t *dimensions,Int_t source) {
  //
  // Slice bin for a given source of electron
  //
  
  Double_t *varMin = new Double_t[container->GetNVar()],
           *varMax = new Double_t[container->GetNVar()];

  Double_t *binLimits;
  for(Int_t ivar = 0; ivar < container->GetNVar(); ivar++){
    
    binLimits = new Double_t[container->GetNBins(ivar)+1];
    container->GetBinLimits(ivar,binLimits);
    if((ivar == 4) && ((source>= 0) && (source<container->GetNBins(ivar)))) {
      varMin[ivar] = binLimits[source];
      varMax[ivar] = binLimits[source];
    }
    else {
      varMin[ivar] = binLimits[0];
      varMax[ivar] = binLimits[container->GetNBins(ivar)];
    }

    delete[] binLimits;
  }
  
  AliCFContainer *k = container->MakeSlice(nDim, dimensions, varMin, varMax);
  AddTemporaryObject(k);
  delete[] varMin; delete[] varMax;

  return k;

}

//_________________________________________________________________________
THnSparseF *AliHFEspectrum::GetSlicedCorrelation(THnSparseF *correlationmatrix, Int_t nDim, Int_t *dimensions) const {
  //
  // Slice correlation
  //

  Int_t ndimensions = correlationmatrix->GetNdimensions();
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
    
  THnSparseF *k = (THnSparseF *) correlationmatrix->Projection(nDim*2,dim);

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
TObject* AliHFEspectrum::GetEfficiency(AliCFContainer * const c, Int_t step, Int_t step0){
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
