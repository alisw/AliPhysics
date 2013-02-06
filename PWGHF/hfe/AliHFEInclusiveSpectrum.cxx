
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
#include <TF1.h>

#include "AliPID.h"
#include "AliCFContainer.h"
#include "AliCFDataGrid.h"
#include "AliCFEffGrid.h"
#include "AliCFGridSparse.h"
#include "AliCFUnfolding.h"
#include "AliLog.h"

#include "AliHFEInclusiveSpectrum.h"
#include "AliHFEInclusiveSpectrumQA.h"
#include "AliHFEcuts.h"
#include "AliHFEcontainer.h"
#include "AliHFEtools.h"

ClassImp(AliHFEInclusiveSpectrum)

//____________________________________________________________
AliHFEInclusiveSpectrum::AliHFEInclusiveSpectrum(const char *name):
  AliHFECorrectSpectrumBase(name),
  fQA(NULL)
{
  //
  // Default constructor
  //

  fQA = new AliHFEInclusiveSpectrumQA("AliHFEInclusiveSpectrumQA");
 
}
//____________________________________________________________
AliHFEInclusiveSpectrum::AliHFEInclusiveSpectrum(const AliHFEInclusiveSpectrum &ref):
  AliHFECorrectSpectrumBase(ref),
  fQA(ref.fQA)
{
  //
  // Copy constructor
  //
  ref.Copy(*this);

}
//____________________________________________________________
AliHFEInclusiveSpectrum &AliHFEInclusiveSpectrum::operator=(const AliHFEInclusiveSpectrum &ref){
  //
  // Assignment operator
  //
  if(this == &ref) 
    ref.Copy(*this);
  return *this;
}
//____________________________________________________________
void AliHFEInclusiveSpectrum::Copy(TObject &o) const {
  // 
  // Copy into object o
  //
  AliHFEInclusiveSpectrum &target = dynamic_cast<AliHFEInclusiveSpectrum &>(o);
  target.fQA = fQA;


}
//____________________________________________________________
AliHFEInclusiveSpectrum::~AliHFEInclusiveSpectrum(){
  //
  // Destructor
  //
  if(fQA) delete fQA;
  
}
//____________________________________________________________
Bool_t AliHFEInclusiveSpectrum::Init(const AliHFEcontainer *datahfecontainer, const AliHFEcontainer *mchfecontainer, const AliHFEcontainer */*bghfecontainer*/, const AliHFEcontainer *v0hfecontainer,AliCFContainer *photoniccontainerD){
  //
  // Init what we need for the correction:
  //
  // Raw spectrum, hadron contamination
  // MC efficiency maps, correlation matrix
  // V0 efficiency if wanted
  //
  // This for a given dimension.
  //
  //

  Int_t kNdim = 3;

  Int_t dims[kNdim];
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

  //
  // Data container: raw spectrum + hadron contamination
  //
  AliCFContainer *datacontainer = datahfecontainer->GetCFContainer("recTrackContReco");
  AliCFContainer *contaminationcontainer = datahfecontainer->GetCFContainer("hadronicBackground");
  if((!datacontainer) || (!contaminationcontainer)) return kFALSE;
  AliCFContainer *datacontainerD = GetSlicedContainer(datacontainer, fNbDimensions, dims, -1, fChargeChoosen,fTestCentralityLow,fTestCentralityHigh);
  AliCFContainer *contaminationcontainerD = GetSlicedContainer(contaminationcontainer, fNbDimensions, dims, -1, fChargeChoosen,fTestCentralityLow,fTestCentralityHigh);
  if((!datacontainerD) || (!contaminationcontainerD)) return kFALSE;
  SetContainer(datacontainerD,AliHFECorrectSpectrumBase::kDataContainer);
  SetContainer(contaminationcontainerD,AliHFECorrectSpectrumBase::kBackgroundData);

  // Photonic Background
  SetContainer(photoniccontainerD,AliHFECorrectSpectrumBase::kPhotonicBackground);

  // QA 
  Int_t dimqa = datacontainer->GetNVar();
  Int_t dimsqa[dimqa];
  for(Int_t i = 0; i < dimqa; i++) dimsqa[i] = i;
  AliCFContainer *datacontainerDQA = GetSlicedContainer(datacontainer, dimqa, dimsqa, -1, fChargeChoosen,fTestCentralityLow,fTestCentralityHigh);
  fQA->AddResultAt(datacontainerDQA,AliHFEInclusiveSpectrumQA::kDataProjection);

  //
  // MC container: ESD/MC efficiency maps + MC/MC efficiency maps 
  //
  AliCFContainer *mccontaineresd = 0x0;
  AliCFContainer *mccontainermc = 0x0;
  mccontaineresd = mchfecontainer->MakeMergedCFContainer("sumesd","sumesd","MCTrackCont:recTrackContReco");
  mccontainermc = mchfecontainer->MakeMergedCFContainer("summc","summc","MCTrackCont:recTrackContMC");
  if((!mccontaineresd) || (!mccontainermc)) return kFALSE;  
  Int_t source = -1;
  AliCFContainer *mccontaineresdD = GetSlicedContainer(mccontaineresd, fNbDimensions, dims, source, fChargeChoosen,fTestCentralityLow,fTestCentralityHigh);
  AliCFContainer *mccontainermcD = GetSlicedContainer(mccontainermc, fNbDimensions, dims, source, fChargeChoosen,fTestCentralityLow,fTestCentralityHigh);
  if((!mccontaineresdD) || (!mccontainermcD)) return kFALSE;
  SetContainer(mccontainermcD,AliHFECorrectSpectrumBase::kMCContainerMC);
  SetContainer(mccontaineresdD,AliHFECorrectSpectrumBase::kMCContainerESD);

  //
  // Correlation matrix
  //  
  THnSparseF *mccorrelation = mchfecontainer->GetCorrelationMatrix("correlationstepbeforePID");
  if(!mccorrelation) mccorrelation = mchfecontainer->GetCorrelationMatrix("correlationstepafterPID");
  if(!mccorrelation) return kFALSE;
  THnSparseF *mccorrelationD = GetSlicedCorrelation(mccorrelation, fNbDimensions, dims,fChargeChoosen,fTestCentralityLow,fTestCentralityHigh);
  if(!mccorrelationD) {
    printf("No correlation\n");
    return kFALSE;
  }
  SetCorrelation(mccorrelationD);
  // QA
  fQA->AddResultAt(mccorrelation,AliHFEInclusiveSpectrumQA::kCMProjection);
 
  //
  // V0 container Electron, pt eta phi
  //
  if(v0hfecontainer) {
    AliCFContainer *containerV0 = v0hfecontainer->GetCFContainer("taggedTrackContainerReco");
    if(!containerV0) return kFALSE;
    AliCFContainer *containerV0Electron = GetSlicedContainer(containerV0, fNbDimensions, dims, AliPID::kElectron,fChargeChoosen,fTestCentralityLow,fTestCentralityHigh);
    if(!containerV0Electron) return kFALSE;
    SetContainer(containerV0Electron,AliHFECorrectSpectrumBase::kDataContainerV0);
  }

  //
  fQA->DrawProjections();

  
  return kTRUE;
}
//____________________________________________________________
Bool_t AliHFEInclusiveSpectrum::Correct(Bool_t subtractcontamination, Bool_t subtractphotonic){
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

  printf("Steps are: stepdata %d, stepMC %d, steptrue %d, stepV0after %d, stepV0before %d\n",fStepData,fStepMC,fStepTrue,fStepAfterCutsV0,fStepBeforeCutsV0);

  ///////////////////////////
  // Check initialization
  ///////////////////////////

  if((!GetContainer(kDataContainer)) || (!GetContainer(kMCContainerMC)) || (!GetContainer(kMCContainerESD))){
    AliInfo("You have to init before");
    return kFALSE;
  }

  if((fStepTrue < 0) && (fStepMC < 0) && (fStepData < 0)) {
    AliInfo("You have to set the steps before: SetMCTruthStep, SetMCEffStep, SetStepToCorrect");
    return kFALSE;
  }

  SetStepGuessedUnfolding(AliHFEcuts::kStepRecKineITSTPC + AliHFEcuts::kNcutStepsMCTrack);

  AliCFDataGrid *dataGridAfterFirstSteps = 0x0;

  //////////////////////////////////
  // Subtract hadron background
  /////////////////////////////////
  AliCFDataGrid *dataspectrumaftersubstraction = 0x0;
  if(subtractcontamination) {
    dataspectrumaftersubstraction = SubtractBackground();
    dataGridAfterFirstSteps = dataspectrumaftersubstraction;
  }

  //////////////////////////////////
  // Subtract Photonic background
  /////////////////////////////////
  AliCFDataGrid *dataspectrumafterphotonicsubstraction = 0x0;
  if(subtractphotonic) {
    dataspectrumafterphotonicsubstraction = SubtractPhotonicBackground();
    dataGridAfterFirstSteps = dataspectrumafterphotonicsubstraction;
  }

  ////////////////////////////////////////////////
  // Correct for TPC efficiency from V0 if any
  ///////////////////////////////////////////////
  AliCFDataGrid *dataspectrumafterV0efficiencycorrection = 0x0;
  AliCFContainer *dataContainerV0 = GetContainer(kDataContainerV0);
  if(dataContainerV0){
    dataspectrumafterV0efficiencycorrection = CorrectV0Efficiency(dataspectrumaftersubstraction);
    dataGridAfterFirstSteps = dataspectrumafterV0efficiencycorrection;
  }

  //////////////////////////////////////////////////////////////////////////////
  // Correct for efficiency parametrized (if TPC efficiency is parametrized)
  /////////////////////////////////////////////////////////////////////////////
  AliCFDataGrid *dataspectrumafterefficiencyparametrizedcorrection = 0x0;
  if(fEfficiencyFunction){
    dataspectrumafterefficiencyparametrizedcorrection = CorrectParametrizedEfficiency(dataGridAfterFirstSteps);
    dataGridAfterFirstSteps = dataspectrumafterefficiencyparametrizedcorrection;  
  }
    
  ///////////////
  // Unfold
  //////////////
  THnSparse *correctedspectrum = Unfold(dataGridAfterFirstSteps);
  if(!correctedspectrum){
    AliError("No corrected spectrum\n");
    return kFALSE;
  }
  
  /////////////////////
  // Simply correct
  ////////////////////
  AliCFDataGrid *alltogetherCorrection = CorrectForEfficiency(dataGridAfterFirstSteps);

  // QA final results
  TGraphErrors* correctedspectrumD = Normalize(correctedspectrum);
  TGraphErrors* alltogetherspectrumD = Normalize(alltogetherCorrection);
  fQA->AddResultAt(correctedspectrumD,AliHFEInclusiveSpectrumQA::kFinalResultUnfolded);
  fQA->AddResultAt(alltogetherspectrumD,AliHFEInclusiveSpectrumQA::kFinalResultDirectEfficiency);
  fQA->DrawResult();

  return kTRUE;
}

//____________________________________________________________
AliCFDataGrid* AliHFEInclusiveSpectrum::SubtractBackground(){
  //
  // Apply background subtraction
  //

  // Raw spectrum
  AliCFContainer *dataContainer = GetContainer(kDataContainer);
  if(!dataContainer){
    AliError("Data Container not available");
    return NULL;
  }
  printf("Step data: %d\n",fStepData);
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
  AliCFDataGrid *backgroundGrid = new AliCFDataGrid("ContaminationGrid","ContaminationGrid",*backgroundContainer,stepbackground);

  // Subtract 
  spectrumSubtracted->Add(backgroundGrid,-1.0);

  // QA
  TH1D *subtractedspectrum = (TH1D *) spectrumSubtracted->Project(0);
  CorrectFromTheWidth(subtractedspectrum);
  TH1D *rawspectrum = (TH1D *) dataspectrumbeforesubstraction->Project(0);
  CorrectFromTheWidth(rawspectrum);
  fQA->AddResultAt(subtractedspectrum,AliHFEInclusiveSpectrumQA::kAfterSC);
  fQA->AddResultAt(rawspectrum,AliHFEInclusiveSpectrumQA::kBeforeSC);
  fQA->DrawSubtractContamination();

  return spectrumSubtracted;
}

//____________________________________________________________
AliCFDataGrid* AliHFEInclusiveSpectrum::SubtractPhotonicBackground(){
  //
  // Apply Photonic background subtraction
  //

  printf("Photonic Background Subtraction \n");

  // Raw spectrum
  AliCFContainer *dataContainer = GetContainer(kDataContainer);
  if(!dataContainer){
    AliError("Data Container not available");
    return NULL;
  }
  printf("Step data: %d\n",fStepData);
  AliCFDataGrid *spectrumPhotonicSubtracted = new AliCFDataGrid("spectrumPhotonicSubtracted", "Data Grid for spectrum after Photonic Background subtraction", *dataContainer,fStepData);

  AliCFDataGrid *dataSpectrumBeforePhotonicSubstraction = (AliCFDataGrid *) ((AliCFDataGrid *)GetSpectrum(GetContainer(kDataContainer),fStepData))->Clone();
  dataSpectrumBeforePhotonicSubstraction->SetName("dataSpectrumBeforePhotonicSubstraction");


  // Background Estimate
  AliCFContainer *photonicContainer = GetContainer(kPhotonicBackground);
  if(!photonicContainer){
    AliError("Photonic background container not found");
    return NULL;
  }

  Int_t stepbackground = 0;
  AliCFDataGrid *photonicGrid = new AliCFDataGrid("ContaminationGrid","ContaminationGrid",*photonicContainer,stepbackground);

  // Subtract
  spectrumPhotonicSubtracted->Add(photonicGrid,-1.0);

  // QA
  TH1D *photonicsubtractedspectrum = (TH1D *) spectrumPhotonicSubtracted->Project(0);
  CorrectFromTheWidth(photonicsubtractedspectrum);
  TH1D *newrawspectrum = (TH1D *) dataSpectrumBeforePhotonicSubstraction->Project(0);
  CorrectFromTheWidth(newrawspectrum);
  fQA->AddResultAt(photonicsubtractedspectrum,AliHFEInclusiveSpectrumQA::kAfterSPB);
  fQA->AddResultAt(newrawspectrum,AliHFEInclusiveSpectrumQA::kBeforeSPB);
  fQA->DrawSubtractPhotonicBackground();

  return spectrumPhotonicSubtracted;
}


//____________________________________________________________
AliCFDataGrid *AliHFEInclusiveSpectrum::CorrectParametrizedEfficiency(AliCFDataGrid* const bgsubpectrum){
  
  //
  // Apply TPC pid efficiency correction from parametrisation
  //

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
  AliCFDataGrid *result = (AliCFDataGrid *) dataGrid->Clone();
  result->SetName("ParametrizedEfficiencyBefore");
  THnSparse *h = result->GetGrid();
  Int_t nbdimensions = h->GetNdimensions();
  //printf("CorrectParametrizedEfficiency::We have dimensions %d\n",nbdimensions);
  AliCFContainer *dataContainer = GetContainer(kDataContainer);
  if(!dataContainer){
    AliError("Data Container not available");
    return NULL;
  }
  AliCFContainer *dataContainerbis = (AliCFContainer *) dataContainer->Clone();
  dataContainerbis->Add(dataContainerbis,-1.0);


  Int_t* coord = new Int_t[nbdimensions];
  memset(coord, 0, sizeof(Int_t) * nbdimensions);
  Double_t* points = new Double_t[nbdimensions];

  ULong64_t nEntries = h->GetNbins();
  for (ULong64_t i = 0; i < nEntries; ++i) {

    Double_t value = h->GetBinContent(i, coord);
    //Double_t valuecontainer = dataContainerbis->GetBinContent(coord,fStepData);
    //printf("Value %f, and valuecontainer %f\n",value,valuecontainer);

    // Get the bin co-ordinates given an coord
    for (Int_t j = 0; j < nbdimensions; ++j)
      points[j] = h->GetAxis(j)->GetBinCenter(coord[j]);

    if (!fEfficiencyFunction->IsInside(points))
         continue;
    TF1::RejectPoint(kFALSE);

    // Evaulate function at points
    Double_t valueEfficiency = fEfficiencyFunction->EvalPar(points, NULL);
    //printf("Value efficiency is %f\n",valueEfficiency);

    if(valueEfficiency > 0.0) {
      h->SetBinContent(coord,value/valueEfficiency);
      dataContainerbis->SetBinContent(coord,fStepData,value/valueEfficiency);
    }
    Double_t error = h->GetBinError(i);
    h->SetBinError(coord,error/valueEfficiency);
    dataContainerbis->SetBinError(coord,fStepData,error/valueEfficiency);

   
  } 

  delete[] coord;
  delete[] points;

  AliCFDataGrid *resultt = new AliCFDataGrid("spectrumEfficiencyParametrized", "Data Grid for spectrum after Efficiency parametrized", *dataContainerbis,fStepData);

 // QA
  TH1D *afterE = (TH1D *) resultt->Project(0);
  CorrectFromTheWidth(afterE);
  TH1D *beforeE = (TH1D *) dataGrid->Project(0);
  CorrectFromTheWidth(beforeE);
  fQA->AddResultAt(afterE,AliHFEInclusiveSpectrumQA::kAfterPE);
  fQA->AddResultAt(beforeE,AliHFEInclusiveSpectrumQA::kBeforePE);
  fQA->AddResultAt(fEfficiencyFunction,AliHFEInclusiveSpectrumQA::kPEfficiency);
  fQA->DrawCorrectWithEfficiency(AliHFEInclusiveSpectrumQA::kParametrized);
  
  return resultt;

}
//____________________________________________________________
AliCFDataGrid *AliHFEInclusiveSpectrum::CorrectV0Efficiency(AliCFDataGrid* const bgsubpectrum){
  
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

  // QA
  TH1D *afterE = (TH1D *) result->Project(0);
  CorrectFromTheWidth(afterE);
  TH1D *beforeE = (TH1D *) dataGrid->Project(0);
  CorrectFromTheWidth(beforeE);
  TH1D* efficiencyDproj = (TH1D *) efficiencyD->Project(0);
  fQA->AddResultAt(afterE,AliHFEInclusiveSpectrumQA::kAfterV0);
  fQA->AddResultAt(beforeE,AliHFEInclusiveSpectrumQA::kBeforeV0);
  fQA->AddResultAt(efficiencyDproj,AliHFEInclusiveSpectrumQA::kV0Efficiency);
  fQA->DrawCorrectWithEfficiency(AliHFEInclusiveSpectrumQA::kV0);
 

  return result;

}
//____________________________________________________________
THnSparse *AliHFEInclusiveSpectrum::Unfold(AliCFDataGrid* const bgsubpectrum){
  
  //
  // Return the unfolded spectrum
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

  AliCFUnfolding unfolding("unfolding","",fNbDimensions,fCorrelation,efficiencyD->GetGrid(),dataGrid->GetGrid(),guessedTHnSparse,1.e-06,0,fNumberOfIterations);
  if(fSetSmoothing) unfolding.UseSmoothing();
  unfolding.Unfold();

  // Results
  THnSparse* result = unfolding.GetUnfolded();
  THnSparse* residual = unfolding.GetEstMeasured();

  // QA
  TH1D *residualh = (TH1D *) residual->Projection(0);
  TH1D *beforeE = (TH1D *) dataGrid->Project(0);
  TH1D* efficiencyDproj = (TH1D *) efficiencyD->Project(0);
  TH1D *afterE = (TH1D *) result->Projection(0);
  CorrectFromTheWidth(residualh);
  CorrectFromTheWidth(beforeE);
  CorrectFromTheWidth(afterE);
  fQA->AddResultAt(residualh,AliHFEInclusiveSpectrumQA::kResidualU);
  fQA->AddResultAt(afterE,AliHFEInclusiveSpectrumQA::kAfterU);
  fQA->AddResultAt(beforeE,AliHFEInclusiveSpectrumQA::kBeforeU);
  fQA->AddResultAt(efficiencyDproj,AliHFEInclusiveSpectrumQA::kUEfficiency);
  fQA->DrawUnfolding();

  return (THnSparse *) result->Clone();

}
//____________________________________________________________
AliCFDataGrid *AliHFEInclusiveSpectrum::CorrectForEfficiency(AliCFDataGrid* const bgsubpectrum){

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

  // QA
  TH1D *afterE = (TH1D *) result->Project(0);
  CorrectFromTheWidth(afterE);
  TH1D *beforeE = (TH1D *) dataGrid->Project(0);
  CorrectFromTheWidth(beforeE);
  TH1D* efficiencyDproj = (TH1D *) efficiencyD->Project(0);
  fQA->AddResultAt(afterE,AliHFEInclusiveSpectrumQA::kAfterMCE);
  fQA->AddResultAt(beforeE,AliHFEInclusiveSpectrumQA::kBeforeMCE);
  fQA->AddResultAt(efficiencyDproj,AliHFEInclusiveSpectrumQA::kMCEfficiency);
  fQA->DrawCorrectWithEfficiency(AliHFEInclusiveSpectrumQA::kMC);

  return result;

}
//____________________________________________________________
void AliHFEInclusiveSpectrum::WriteResults(const char *filename)
{
  //
  // Write results
  //

  AliCFContainer *dataContainer = GetContainer(kDataContainer);
  AliCFContainer *mcContainer = GetContainer(kMCContainerMC);
  TObject *unfolded = 0x0;
  TObject *correctedspectrum = 0x0;
  if(fQA) {
    unfolded = fQA->GetResult(AliHFEInclusiveSpectrumQA::kFinalResultUnfolded);
    correctedspectrum = fQA->GetResult(AliHFEInclusiveSpectrumQA::kFinalResultDirectEfficiency);
  }

  TFile *file = TFile::Open(filename,"recreate");
  if(dataContainer) dataContainer->Write("data");
  if(mcContainer) mcContainer->Write("mcefficiency");
  if(fCorrelation) fCorrelation->Write("correlationmatrix");
  if(unfolded) unfolded->Write("unfoldedspectrum");
  if(correctedspectrum) correctedspectrum->Write("correctedspectrum");
  if(fQA) fQA->Write("QAResults");
  file->Close();

}

