
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

#include "AliHFEspectrum.h"
#include "AliHFEcuts.h"
#include "AliHFEcontainer.h"
#include "AliHFEtools.h"

ClassImp(AliHFEspectrum)

//____________________________________________________________
AliHFEspectrum::AliHFEspectrum(const char *name):
  TNamed(name, ""),
  fCFContainers(new TObjArray(kDataContainerV0+1)),
  fTemporaryObjects(NULL),
  fCorrelation(NULL),
  fBackground(NULL),
  fEfficiencyFunction(NULL),
  fWeightCharm(NULL),
  fInclusiveSpectrum(kTRUE),
  fDumpToFile(kFALSE),
  fEtaSelected(kFALSE),
  fUnSetCorrelatedErrors(kTRUE),
  fSetSmoothing(kFALSE),
  fIPanaHadronBgSubtract(kFALSE),
  fIPanaCharmBgSubtract(kFALSE),
  fIPanaConversionBgSubtract(kFALSE),
  fIPanaNonHFEBgSubtract(kFALSE),
  fIPParameterizedEff(kFALSE),
  fNonHFEsyst(0),
  fBeauty2ndMethod(kFALSE),
  fIPEffCombinedSamples(kTRUE),
  fNbDimensions(1),
  fNMCEvents(0),
  fStepMC(-1),
  fStepTrue(-1),
  fStepData(-1),
  fStepBeforeCutsV0(-1),
  fStepAfterCutsV0(-1),
  fStepGuessedUnfolding(-1),
  fNumberOfIterations(5),
  fNRandomIter(0),
  fChargeChoosen(kAllCharge),
  fNCentralityBinAtTheEnd(0),
  fTestCentralityLow(-1),
  fTestCentralityHigh(-1),
  fFillMoreCorrelationMatrix(kFALSE),
  fHadronEffbyIPcut(NULL),
  fConversionEffbgc(NULL),
  fNonHFEEffbgc(NULL),      
  fBSpectrum2ndMethod(NULL),
  fkBeauty2ndMethodfilename(""),
  fBeamType(0),
  fEtaSyst(kTRUE),
  fDebugLevel(0),
  fWriteToFile(kFALSE),
  fUnfoldBG(kFALSE)
{
  //
  // Default constructor
  //

  for(Int_t k = 0; k < 20; k++){
      fNEvents[k] = 0;
      fNMCbgEvents[k] = 0;
      fLowBoundaryCentralityBinAtTheEnd[k] = 0;
      fHighBoundaryCentralityBinAtTheEnd[k] = 0;
      if(k<kCentrality)
      {
          fEfficiencyTOFPIDD[k] = 0;
          fEfficiencyesdTOFPIDD[k] = 0;
          fEfficiencyIPCharmD[k] = 0;     
          fEfficiencyIPBeautyD[k] = 0;    
          fEfficiencyIPBeautyesdD[k] = 0;
          fEfficiencyIPConversionD[k] = 0;
          fEfficiencyIPNonhfeD[k] = 0;   

	  fConversionEff[k] = 0;
	  fNonHFEEff[k] = 0;
	  fCharmEff[k] = 0;
	  fBeautyEff[k] = 0;
	  fEfficiencyCharmSigD[k] = 0;
          fEfficiencyBeautySigD[k] = 0;
          fEfficiencyBeautySigesdD[k] = 0;
      }
  }
  memset(fEtaRange, 0, sizeof(Double_t) * 2);
  memset(fEtaRangeNorm, 0, sizeof(Double_t) * 2);
  memset(fConvSourceContainer, 0, sizeof(AliCFContainer*) * kElecBgSources * kBgLevels);
  memset(fNonHFESourceContainer, 0, sizeof(AliCFContainer*) * kElecBgSources * kBgLevels);
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
Bool_t AliHFEspectrum::Init(const AliHFEcontainer *datahfecontainer, const AliHFEcontainer *mchfecontainer, const AliHFEcontainer *v0hfecontainer, const AliHFEcontainer *bghfecontainer){
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

  
  if(fBeauty2ndMethod) CallInputFileForBeauty2ndMethod();

  Int_t kNdim = 3;
  Int_t kNcentr =1;
  //Int_t ptpr =0;
  if(fBeamType==0) kNdim=3;
  if(fBeamType==1)
  {
      kNdim=4;
      kNcentr=11;
      //ptpr=1;
  }

  Int_t dims[kNdim];
  // Get the requested format
  if(fBeamType==0){
    // pp analysis
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
  }

  if(fBeamType==1){
    // PbPb analysis; centrality as first dimension
    Int_t nbDimensions = fNbDimensions;
    fNbDimensions = fNbDimensions + 1;
    switch(nbDimensions){
      case 1: dims[0] = 5;
              dims[1] = 0;
              break;
      case 2: dims[0] = 5;
              for(Int_t i = 0; i < 2; i++) dims[(i+1)] = i;
              break;
      case 3: dims[0] = 5;
              for(Int_t i = 0; i < 3; i++) dims[(i+1)] = i;
              break;
      default:
              AliError("Container with this number of dimensions not foreseen (yet)");
              return kFALSE;
    };
  }

  // Data container: raw spectrum + hadron contamination  
  AliCFContainer *datacontainer = 0x0; 
  if(fInclusiveSpectrum) {
    datacontainer = datahfecontainer->GetCFContainer("recTrackContReco");
  }
  else{

      datacontainer = datahfecontainer->MakeMergedCFContainer("sumreco","sumreco","recTrackContReco:recTrackContDEReco");
  }
  AliCFContainer *contaminationcontainer = datahfecontainer->GetCFContainer("hadronicBackground");
  if((!datacontainer) || (!contaminationcontainer)) return kFALSE;

  AliCFContainer *datacontainerD = GetSlicedContainer(datacontainer, fNbDimensions, dims, -1, fChargeChoosen,fTestCentralityLow,fTestCentralityHigh);
  AliCFContainer *contaminationcontainerD = GetSlicedContainer(contaminationcontainer, fNbDimensions, dims, -1, fChargeChoosen,fTestCentralityLow,fTestCentralityHigh);
  if((!datacontainerD) || (!contaminationcontainerD)) return kFALSE;

  SetContainer(datacontainerD,AliHFEspectrum::kDataContainer);
  SetContainer(contaminationcontainerD,AliHFEspectrum::kBackgroundData);

  // MC container: ESD/MC efficiency maps + MC/MC efficiency maps 
  AliCFContainer *mccontaineresd = 0x0;
  AliCFContainer *mccontaineresdbg = 0x0;
  AliCFContainer *mccontainermc = 0x0;
  AliCFContainer *mccontainermcbg = 0x0;
  AliCFContainer *nonHFEweightedContainer = 0x0;
  AliCFContainer *convweightedContainer = 0x0;
  AliCFContainer *nonHFEtempContainer = 0x0;//temporary container to be sliced for the fnonHFESourceContainers
  AliCFContainer *convtempContainer = 0x0;//temporary container to be sliced for the fConvSourceContainers
   
  if(fInclusiveSpectrum) {
    mccontaineresd = mchfecontainer->MakeMergedCFContainer("sumesd","sumesd","MCTrackCont:recTrackContReco");
    mccontainermc = mchfecontainer->MakeMergedCFContainer("summc","summc","MCTrackCont:recTrackContMC");
  }
  else {
    mccontaineresd = mchfecontainer->MakeMergedCFContainer("sumesd","sumesd","MCTrackCont:recTrackContReco:recTrackContDEReco");
    mccontaineresdbg = bghfecontainer->MakeMergedCFContainer("sumesd","sumesd","MCTrackCont:recTrackContReco:recTrackContDEReco");
    mccontainermc = mchfecontainer->MakeMergedCFContainer("summc","summc","MCTrackCont:recTrackContMC:recTrackContDEMC");
    mccontainermcbg = bghfecontainer->MakeMergedCFContainer("summcbg","summcbg","MCTrackCont:recTrackContMC:recTrackContDEMC");

    if(fNonHFEsyst){   
      const Char_t *sourceName[kElecBgSources]={"Pion","Eta","Omega","Phi","EtaPrime","Rho"};
      const Char_t *levelName[kBgLevels]={"Best","Lower","Upper"};
      for(Int_t iSource = 0; iSource < kElecBgSources; iSource++){
	for(Int_t iLevel = 0; iLevel < kBgLevels; iLevel++){
          nonHFEtempContainer =  bghfecontainer->GetCFContainer(Form("mesonElecs%s%s",sourceName[iSource],levelName[iLevel]));
	  convtempContainer =  bghfecontainer->GetCFContainer(Form("conversionElecs%s%s",sourceName[iSource],levelName[iLevel]));
	  for(Int_t icentr=0;icentr<kNcentr;icentr++)
	  {
	      if(fBeamType==0)
	      {
		  fConvSourceContainer[iSource][iLevel][icentr] = GetSlicedContainer(convtempContainer, fNbDimensions, dims, -1, fChargeChoosen);
		  fNonHFESourceContainer[iSource][iLevel][icentr] = GetSlicedContainer(nonHFEtempContainer, fNbDimensions, dims, -1, fChargeChoosen);
	      }
	      if(fBeamType==1)
	      {
              
		fConvSourceContainer[iSource][iLevel][icentr] = GetSlicedContainer(convtempContainer, fNbDimensions, dims, -1, fChargeChoosen,icentr,icentr);
		fNonHFESourceContainer[iSource][iLevel][icentr] = GetSlicedContainer(nonHFEtempContainer, fNbDimensions, dims, -1, fChargeChoosen,icentr,icentr);
	      }
//	      if((!fConvSourceContainer[iSource][iLevel][icentr])||(!fNonHFESourceContainer[iSource][iLevel])) return kFALSE;
	  }
          if(fBeamType == 1)break;
	}
      }
    }
    // else{      
      nonHFEweightedContainer = bghfecontainer->GetCFContainer("mesonElecs");
      convweightedContainer = bghfecontainer->GetCFContainer("conversionElecs");
      if((!convweightedContainer)||(!nonHFEweightedContainer)) return kFALSE;  
      //}
  }
  if((!mccontaineresd) || (!mccontainermc)) return kFALSE;  
  
  Int_t source = -1;
  if(!fInclusiveSpectrum) source = 1; //beauty
  AliCFContainer *mccontaineresdD = GetSlicedContainer(mccontaineresd, fNbDimensions, dims, source, fChargeChoosen,fTestCentralityLow,fTestCentralityHigh);
  AliCFContainer *mccontainermcD = GetSlicedContainer(mccontainermc, fNbDimensions, dims, source, fChargeChoosen,fTestCentralityLow,fTestCentralityHigh);
  if((!mccontaineresdD) || (!mccontainermcD)) return kFALSE;
  SetContainer(mccontainermcD,AliHFEspectrum::kMCContainerMC);
  SetContainer(mccontaineresdD,AliHFEspectrum::kMCContainerESD);

  // set charm, nonHFE container to estimate BG
  if(!fInclusiveSpectrum){
   source = 0; //charm
   mccontainermcD = GetSlicedContainer(mccontainermcbg, fNbDimensions, dims, source, fChargeChoosen);
   SetContainer(mccontainermcD,AliHFEspectrum::kMCContainerCharmMC);

   //if(!fNonHFEsyst){
     AliCFContainer *nonHFEweightedContainerD = GetSlicedContainer(nonHFEweightedContainer, fNbDimensions, dims, -1, fChargeChoosen);
     SetContainer(nonHFEweightedContainerD,AliHFEspectrum::kMCWeightedContainerNonHFEESD);
     AliCFContainer *convweightedContainerD = GetSlicedContainer(convweightedContainer, fNbDimensions, dims, -1, fChargeChoosen);
     SetContainer(convweightedContainerD,AliHFEspectrum::kMCWeightedContainerConversionESD);
     //}

     SetParameterizedEff(mccontainermc, mccontainermcbg, mccontaineresd, mccontaineresdbg, dims);

  }
  // MC container: correlation matrix
  THnSparseF *mccorrelation = 0x0;
  if(fInclusiveSpectrum) {
    if(fStepMC==(AliHFEcuts::kNcutStepsMCTrack + AliHFEcuts::kStepHFEcutsTRD + 2)) mccorrelation = mchfecontainer->GetCorrelationMatrix("correlationstepbeforePID");
    else if(fStepMC==(AliHFEcuts::kNcutStepsMCTrack + AliHFEcuts::kStepHFEcutsTRD + 1)) mccorrelation = mchfecontainer->GetCorrelationMatrix("correlationstepbeforePID");
    else if(fStepMC==(AliHFEcuts::kNcutStepsMCTrack + AliHFEcuts::kStepHFEcutsTRD)) mccorrelation = mchfecontainer->GetCorrelationMatrix("correlationstepbeforePID");
    else if(fStepMC==(AliHFEcuts::kNcutStepsMCTrack + AliHFEcuts::kStepHFEcutsTRD - 1)) mccorrelation = mchfecontainer->GetCorrelationMatrix("correlationstepbeforePID");
    else mccorrelation = mchfecontainer->GetCorrelationMatrix("correlationstepafterPID");
    
    if(!mccorrelation) mccorrelation = mchfecontainer->GetCorrelationMatrix("correlationstepafterPID");
  }
  else mccorrelation = mchfecontainer->GetCorrelationMatrix("correlationstepafterPID"); // we confirmed that we get same result by using it instead of correlationstepafterDE
  //else mccorrelation = mchfecontainer->GetCorrelationMatrix("correlationstepafterDE");
  if(!mccorrelation) return kFALSE;
  Int_t testCentralityLow = fTestCentralityLow;
  Int_t testCentralityHigh = fTestCentralityHigh;
  if(fFillMoreCorrelationMatrix) {
    testCentralityLow = fTestCentralityLow-1;
    testCentralityHigh = fTestCentralityHigh+1;
  }
  THnSparseF *mccorrelationD = GetSlicedCorrelation(mccorrelation, fNbDimensions, dims,testCentralityLow,testCentralityHigh);
  if(!mccorrelationD) {
    printf("No correlation\n");
    return kFALSE;
  }
  SetCorrelation(mccorrelationD);

  // V0 container Electron, pt eta phi
  if(v0hfecontainer) {
    AliCFContainer *containerV0 = v0hfecontainer->GetCFContainer("taggedTrackContainerReco");
    if(!containerV0) return kFALSE;
    AliCFContainer *containerV0Electron = GetSlicedContainer(containerV0, fNbDimensions, dims, AliPID::kElectron,fChargeChoosen,fTestCentralityLow,fTestCentralityHigh);
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

    TCanvas * ccorrelation = new TCanvas("ccorrelation","ccorrelation",1000,700);
    ccorrelation->cd(1);
    if(mccorrelationD) {
      if(fBeamType==0){
        TH2D * ptcorrelation = (TH2D *) mccorrelationD->Projection(fNbDimensions,0);
        ptcorrelation->Draw("colz");
      }
      if(fBeamType==1){
        TH2D * ptcorrelation = (TH2D *) mccorrelationD->Projection(fNbDimensions+1,1);
        ptcorrelation->Draw("colz");
      }
    }
    if(fWriteToFile){ 
      ccontaminationspectrum->SaveAs("contaminationspectrum.eps");
      ccorrelation->SaveAs("correlationMatrix.eps");
    }
  }

  TFile *file = TFile::Open("tests.root","recreate");
  datacontainerD->Write("data");
  mccontainermcD->Write("mcefficiency");
  mccorrelationD->Write("correlationmatrix");
  file->Close();

  return kTRUE;
}


//____________________________________________________________
void AliHFEspectrum::CallInputFileForBeauty2ndMethod(){
  //
  // get spectrum for beauty 2nd method
  //
  //
    TFile *inputfile2ndmethod=TFile::Open(fkBeauty2ndMethodfilename);
    fBSpectrum2ndMethod = new TH1D(*static_cast<TH1D*>(inputfile2ndmethod->Get("BSpectrum")));
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
 
  SetNumberOfIteration(10);
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

  printf("cloning spectrum\n");
  AliCFDataGrid *rawsave(NULL);
  if(dataspectrumaftersubstraction){
     rawsave = (AliCFDataGrid *)dataspectrumaftersubstraction->Clone("rawdata");
  } else {
    AliCFContainer *dataContainer = GetContainer(kDataContainer);
    if(!dataContainer){
      AliError("Data Container not available");
    }
    rawsave = new AliCFDataGrid("rawsave", "raw spectrum after subtraction",*dataContainer, fStepData);
  }
  printf("cloned: %p\n", rawsave);

  ////////////////////////////////////////////////
  // Correct for TPC efficiency from V0
  ///////////////////////////////////////////////
  AliCFDataGrid *dataspectrumafterV0efficiencycorrection = 0x0;
  AliCFContainer *dataContainerV0 = GetContainer(kDataContainerV0);
  if(dataContainerV0){printf("Got the V0 container\n");
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

    Int_t ptpr = 0;
    if(fBeamType==0) ptpr=0;
    if(fBeamType==1) ptpr=1;
  
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
    TH1D *correctedTH1D = correctedspectrum->Projection(ptpr);
    TH1D *alltogetherTH1D = (TH1D *) alltogetherCorrection->Project(ptpr);
    TH1D* ratiocorrected = (TH1D*)correctedTH1D->Clone();
    ratiocorrected->SetName("ratiocorrected");
    ratiocorrected->SetTitle("");
    ratiocorrected->GetYaxis()->SetTitle("Unfolded/DirectCorrected");
    ratiocorrected->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    ratiocorrected->Divide(correctedTH1D,alltogetherTH1D,1,1);
    ratiocorrected->SetStats(0);
    ratiocorrected->Draw();
    if(fWriteToFile)ccorrected->SaveAs("CorrectedPbPb.eps");

    //TH1D unfoldingspectrac[fNCentralityBinAtTheEnd];
    //TGraphErrors unfoldingspectracn[fNCentralityBinAtTheEnd];
    //TH1D correctedspectrac[fNCentralityBinAtTheEnd];
    //TGraphErrors correctedspectracn[fNCentralityBinAtTheEnd];

    TH1D *unfoldingspectrac = new TH1D[fNCentralityBinAtTheEnd];
    TGraphErrors *unfoldingspectracn = new TGraphErrors[fNCentralityBinAtTheEnd];
    TH1D *correctedspectrac = new TH1D[fNCentralityBinAtTheEnd];
    TGraphErrors *correctedspectracn = new TGraphErrors[fNCentralityBinAtTheEnd];

    if(fBeamType==1) {

      TCanvas * ccorrectedallspectra = new TCanvas("correctedallspectra","correctedallspectra",1000,700);
      ccorrectedallspectra->Divide(2,1);
      TLegend *legtotal = new TLegend(0.4,0.6,0.89,0.89);
      TLegend *legtotalg = new TLegend(0.4,0.6,0.89,0.89);
      
      THnSparseF* sparsesu = (THnSparseF *) correctedspectrum;
      TAxis *cenaxisa = sparsesu->GetAxis(0);
      THnSparseF* sparsed = (THnSparseF *) alltogetherCorrection->GetGrid();
      TAxis *cenaxisb = sparsed->GetAxis(0);
      Int_t nbbin = cenaxisb->GetNbins();
      Int_t stylee[20] = {20,21,22,23,24,25,26,27,28,30,4,5,7,29,29,29,29,29,29,29};
      Int_t colorr[20] = {2,3,4,5,6,7,8,9,46,38,29,30,31,32,33,34,35,37,38,20};
      for(Int_t binc = 0; binc < fNCentralityBinAtTheEnd; binc++){
        TString titlee("corrected_centrality_bin_");
        titlee += "[";
        titlee += fLowBoundaryCentralityBinAtTheEnd[binc];
        titlee += "_";
        titlee += fHighBoundaryCentralityBinAtTheEnd[binc];
        titlee += "[";
        TString titleec("corrected_check_projection_bin_");
        titleec += "[";
        titleec += fLowBoundaryCentralityBinAtTheEnd[binc];
        titleec += "_";
        titleec += fHighBoundaryCentralityBinAtTheEnd[binc];
        titleec += "[";
        TString titleenameunotnormalized("Unfolded_Notnormalized_centrality_bin_");
        titleenameunotnormalized += "[";
        titleenameunotnormalized += fLowBoundaryCentralityBinAtTheEnd[binc];
        titleenameunotnormalized += "_";
        titleenameunotnormalized += fHighBoundaryCentralityBinAtTheEnd[binc];
        titleenameunotnormalized += "[";
        TString titleenameunormalized("Unfolded_normalized_centrality_bin_");
        titleenameunormalized += "[";
        titleenameunormalized += fLowBoundaryCentralityBinAtTheEnd[binc];
        titleenameunormalized += "_";
        titleenameunormalized += fHighBoundaryCentralityBinAtTheEnd[binc];	
        titleenameunormalized += "[";
        TString titleenamednotnormalized("Dirrectcorrected_Notnormalized_centrality_bin_");
        titleenamednotnormalized += "[";
        titleenamednotnormalized += fLowBoundaryCentralityBinAtTheEnd[binc];
        titleenamednotnormalized += "_";
        titleenamednotnormalized += fHighBoundaryCentralityBinAtTheEnd[binc];
        titleenamednotnormalized += "[";
        TString titleenamednormalized("Dirrectedcorrected_normalized_centrality_bin_");
        titleenamednormalized += "[";
        titleenamednormalized += fLowBoundaryCentralityBinAtTheEnd[binc];
        titleenamednormalized += "_";
        titleenamednormalized += fHighBoundaryCentralityBinAtTheEnd[binc];	
        titleenamednormalized += "[";
        Int_t nbEvents = 0;
        for(Int_t k = fLowBoundaryCentralityBinAtTheEnd[binc]; k < fHighBoundaryCentralityBinAtTheEnd[binc]; k++) {
          printf("Number of events %d in the bin %d added!!!\n",fNEvents[k],k);
          nbEvents += fNEvents[k];
        }
        Double_t lowedgega = cenaxisa->GetBinLowEdge(fLowBoundaryCentralityBinAtTheEnd[binc]+1);
        Double_t upedgega = cenaxisa->GetBinUpEdge(fHighBoundaryCentralityBinAtTheEnd[binc]);
        printf("Bin Low edge %f, up edge %f for a\n",lowedgega,upedgega);
        Double_t lowedgegb = cenaxisb->GetBinLowEdge(fLowBoundaryCentralityBinAtTheEnd[binc]+1);
        Double_t upedgegb = cenaxisb->GetBinUpEdge(fHighBoundaryCentralityBinAtTheEnd[binc]);
        printf("Bin Low edge %f, up edge %f for b\n",lowedgegb,upedgegb);
        cenaxisa->SetRange(fLowBoundaryCentralityBinAtTheEnd[binc]+1,fHighBoundaryCentralityBinAtTheEnd[binc]);
        cenaxisb->SetRange(fLowBoundaryCentralityBinAtTheEnd[binc]+1,fHighBoundaryCentralityBinAtTheEnd[binc]);
        TCanvas * ccorrectedcheck = new TCanvas((const char*) titleec,(const char*) titleec,1000,700);
        ccorrectedcheck->cd(1);
        TH1D *aftersuc = (TH1D *) sparsesu->Projection(0);
        TH1D *aftersdc = (TH1D *) sparsed->Projection(0);
        aftersuc->Draw();
        aftersdc->Draw("same");
        TCanvas * ccorrectede = new TCanvas((const char*) titlee,(const char*) titlee,1000,700);
        ccorrectede->Divide(2,1);
        ccorrectede->cd(1);
        gPad->SetLogy();
        TH1D *aftersu = (TH1D *) sparsesu->Projection(1);
        CorrectFromTheWidth(aftersu);
        aftersu->SetName((const char*)titleenameunotnormalized);
        unfoldingspectrac[binc] = *aftersu;
        ccorrectede->cd(1);
        TGraphErrors* aftersun = NormalizeTH1N(aftersu,nbEvents);
        aftersun->SetTitle("");
        aftersun->GetYaxis()->SetTitleOffset(1.5);
        aftersun->GetYaxis()->SetRangeUser(0.000000001,1.0);
        aftersun->SetMarkerStyle(26);
        aftersun->SetMarkerColor(kBlue);
        aftersun->SetLineColor(kBlue);
        aftersun->Draw("AP");
        aftersun->SetName((const char*)titleenameunormalized);
        unfoldingspectracn[binc] = *aftersun;
        ccorrectede->cd(1);
        TH1D *aftersd = (TH1D *) sparsed->Projection(1);
        CorrectFromTheWidth(aftersd);
        aftersd->SetName((const char*)titleenamednotnormalized);
        correctedspectrac[binc] = *aftersd;
        ccorrectede->cd(1);
        TGraphErrors* aftersdn = NormalizeTH1N(aftersd,nbEvents);
        aftersdn->SetTitle("");
        aftersdn->GetYaxis()->SetTitleOffset(1.5);
        aftersdn->GetYaxis()->SetRangeUser(0.000000001,1.0);
        aftersdn->SetMarkerStyle(25);
        aftersdn->SetMarkerColor(kBlack);
        aftersdn->SetLineColor(kBlack);
        aftersdn->Draw("P");
        aftersdn->SetName((const char*)titleenamednormalized);
        correctedspectracn[binc] = *aftersdn;
        TLegend *legcorrectedud = new TLegend(0.4,0.6,0.89,0.89);
        legcorrectedud->AddEntry(aftersun,"Corrected","p");
        legcorrectedud->AddEntry(aftersdn,"Alltogether","p");
        legcorrectedud->Draw("same");
        ccorrectedallspectra->cd(1);
        gPad->SetLogy();
        TH1D *aftersunn = (TH1D *) aftersun->Clone();
        aftersunn->SetMarkerStyle(stylee[binc]);
        aftersunn->SetMarkerColor(colorr[binc]);
        if(binc==0) aftersunn->Draw("AP");
        else aftersunn->Draw("P");
        legtotal->AddEntry(aftersunn,(const char*) titlee,"p");
        ccorrectedallspectra->cd(2);
        gPad->SetLogy();
        TH1D *aftersdnn = (TH1D *) aftersdn->Clone();
        aftersdnn->SetMarkerStyle(stylee[binc]);
        aftersdnn->SetMarkerColor(colorr[binc]);
        if(binc==0) aftersdnn->Draw("AP");
        else aftersdnn->Draw("P");
        legtotalg->AddEntry(aftersdnn,(const char*) titlee,"p");
        ccorrectede->cd(2);
        TH1D* ratiocorrectedbinc = (TH1D*)aftersu->Clone();
        TString titleee("ratiocorrected_bin_");
        titleee += binc;
        ratiocorrectedbinc->SetName((const char*) titleee);
        ratiocorrectedbinc->SetTitle("");
        ratiocorrectedbinc->GetYaxis()->SetTitle("Unfolded/DirectCorrected");
        ratiocorrectedbinc->GetXaxis()->SetTitle("p_{T} [GeV/c]");
        ratiocorrectedbinc->Divide(aftersu,aftersd,1,1);
        ratiocorrectedbinc->SetStats(0);
        ratiocorrectedbinc->Draw();
      }

      ccorrectedallspectra->cd(1);
      legtotal->Draw("same");
      ccorrectedallspectra->cd(2);
      legtotalg->Draw("same");
      
      cenaxisa->SetRange(0,nbbin);
      cenaxisb->SetRange(0,nbbin);
      if(fWriteToFile) ccorrectedallspectra->SaveAs("CorrectedPbPb.eps");
    }

    // Dump to file if needed
    if(fDumpToFile) {
      TFile *out = new TFile("finalSpectrum.root","recreate");
      correctedspectrumD->SetName("UnfoldingCorrectedSpectrum");
      correctedspectrumD->Write();
      alltogetherspectrumD->SetName("AlltogetherSpectrum");
      alltogetherspectrumD->Write();
      ratiocorrected->SetName("RatioUnfoldingAlltogetherSpectrum");
      ratiocorrected->Write();
      correctedspectrum->SetName("UnfoldingCorrectedNotNormalizedSpectrum");
      correctedspectrum->Write();
      alltogetherCorrection->SetName("AlltogetherCorrectedNotNormalizedSpectrum");
      alltogetherCorrection->Write();
      rawsave->Write();
      for(Int_t binc = 0; binc < fNCentralityBinAtTheEnd; binc++){
	      unfoldingspectrac[binc].Write();
	      unfoldingspectracn[binc].Write();
	      correctedspectrac[binc].Write();
	      correctedspectracn[binc].Write();
      }
      out->Close(); delete out;
    }

    if (unfoldingspectrac) delete[] unfoldingspectrac;
    if (unfoldingspectracn)  delete[] unfoldingspectracn;
    if (correctedspectrac) delete[] correctedspectrac;
    if (correctedspectracn) delete[] correctedspectracn;

  }

  return kTRUE;
}

//____________________________________________________________
Bool_t AliHFEspectrum::CorrectBeauty(Bool_t subtractcontamination){
  //
  // Correct the spectrum for efficiency and unfolding for beauty analysis
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
 
  SetNumberOfIteration(10);
  SetStepGuessedUnfolding(AliHFEcuts::kStepRecKineITSTPC + AliHFEcuts::kNcutStepsMCTrack);
    
  AliCFDataGrid *dataGridAfterFirstSteps = 0x0;
  //////////////////////////////////
  // Subtract hadron background
  /////////////////////////////////
  AliCFDataGrid *dataspectrumaftersubstraction = 0x0;
  AliCFDataGrid *unnormalizedRawSpectrum = 0x0;
  TGraphErrors *gNormalizedRawSpectrum = 0x0;
  if(subtractcontamination) {
      if(!fBeauty2ndMethod) dataspectrumaftersubstraction = SubtractBackground(kTRUE);
      else dataspectrumaftersubstraction = GetRawBspectra2ndMethod();
      unnormalizedRawSpectrum = (AliCFDataGrid*)dataspectrumaftersubstraction->Clone();
      dataGridAfterFirstSteps = dataspectrumaftersubstraction;
      gNormalizedRawSpectrum = Normalize(unnormalizedRawSpectrum);
  }

  printf("after normalize getting IP \n");

  /////////////////////////////////////////////////////////////////////////////////////////
  // Correct for IP efficiency for beauty electrons after subtracting all the backgrounds
  /////////////////////////////////////////////////////////////////////////////////////////

  AliCFDataGrid *dataspectrumafterefficiencyparametrizedcorrection = 0x0;
  AliCFDataGrid *dataspectrumafterV0efficiencycorrection = 0x0;
  AliCFContainer *dataContainerV0 = GetContainer(kDataContainerV0);

  if(fEfficiencyFunction){
    dataspectrumafterefficiencyparametrizedcorrection = CorrectParametrizedEfficiency(dataGridAfterFirstSteps);
    dataGridAfterFirstSteps = dataspectrumafterefficiencyparametrizedcorrection;
  }
  else if(dataContainerV0){
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
    AliError("No residual spectrum\n");
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

    Int_t ptpr = 0;
    if(fBeamType==0) ptpr=0;
    if(fBeamType==1) ptpr=1;
  
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
    TH1D *correctedTH1D = correctedspectrum->Projection(ptpr);
    TH1D *alltogetherTH1D = (TH1D *) alltogetherCorrection->Project(ptpr);
    TH1D* ratiocorrected = (TH1D*)correctedTH1D->Clone();
    ratiocorrected->SetName("ratiocorrected");
    ratiocorrected->SetTitle("");
    ratiocorrected->GetYaxis()->SetTitle("Unfolded/DirectCorrected");
    ratiocorrected->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    ratiocorrected->Divide(correctedTH1D,alltogetherTH1D,1,1);
    ratiocorrected->SetStats(0);
    ratiocorrected->Draw();
    if(fWriteToFile) ccorrected->SaveAs("CorrectedBeauty.eps");

    if(fBeamType == 0){
	if(fNonHFEsyst){
	    CalculateNonHFEsyst(0);
	}
    }

    // Dump to file if needed

    if(fDumpToFile) {
        // to do centrality dependent

      TFile *out;
      out = new TFile("finalSpectrum.root","recreate");
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
      if(unnormalizedRawSpectrum) {
      unnormalizedRawSpectrum->SetName("beautyAfterIP");
      unnormalizedRawSpectrum->Write();
      }

      if(gNormalizedRawSpectrum){
        gNormalizedRawSpectrum->SetName("normalizedBeautyAfterIP");
        gNormalizedRawSpectrum->Write();
      }

      if(fBeamType==0) {
          Int_t countpp=0;
	  fEfficiencyCharmSigD[countpp]->SetTitle(Form("IPEfficiencyForCharmSigCent%i",countpp));
	  fEfficiencyCharmSigD[countpp]->SetName(Form("IPEfficiencyForCharmSigCent%i",countpp));
	  fEfficiencyCharmSigD[countpp]->Write();
	  fEfficiencyBeautySigD[countpp]->SetTitle(Form("IPEfficiencyForBeautySigCent%i",countpp));
	  fEfficiencyBeautySigD[countpp]->SetName(Form("IPEfficiencyForBeautySigCent%i",countpp));
	  fEfficiencyBeautySigD[countpp]->Write();
	  fCharmEff[countpp]->SetTitle(Form("IPEfficiencyForCharmCent%i",countpp));
	  fCharmEff[countpp]->SetName(Form("IPEfficiencyForCharmCent%i",countpp));
	  fCharmEff[countpp]->Write();
	  fBeautyEff[countpp]->SetTitle(Form("IPEfficiencyForBeautyCent%i",countpp));
	  fBeautyEff[countpp]->SetName(Form("IPEfficiencyForBeautyCent%i",countpp));
	  fBeautyEff[countpp]->Write();
	  fConversionEff[countpp]->SetTitle(Form("IPEfficiencyForConversionCent%i",countpp));
	  fConversionEff[countpp]->SetName(Form("IPEfficiencyForConversionCent%i",countpp));
	  fConversionEff[countpp]->Write();
	  fNonHFEEff[countpp]->SetTitle(Form("IPEfficiencyForNonHFECent%i",countpp));
	  fNonHFEEff[countpp]->SetName(Form("IPEfficiencyForNonHFECent%i",countpp));
	  fNonHFEEff[countpp]->Write();
      }

      if(fBeamType==1) {

	  TGraphErrors* correctedspectrumDc[kCentrality];
          TGraphErrors* alltogetherspectrumDc[kCentrality];
	  for(Int_t i=0;i<kCentrality-2;i++)
	  {
	      correctedspectrum->GetAxis(0)->SetRange(i+1,i+1);
	      correctedspectrumDc[i] = Normalize(correctedspectrum,i);
              if(correctedspectrumDc[i]){
                correctedspectrumDc[i]->SetTitle(Form("UnfoldingCorrectedSpectrum_%i",i));
                correctedspectrumDc[i]->SetName(Form("UnfoldingCorrectedSpectrum_%i",i));
                correctedspectrumDc[i]->Write();
              }
	      alltogetherCorrection->GetAxis(0)->SetRange(i+1,i+1);
	      alltogetherspectrumDc[i] = Normalize(alltogetherCorrection,i);
              if(alltogetherspectrumDc[i]){
                alltogetherspectrumDc[i]->SetTitle(Form("AlltogetherSpectrum_%i",i));
                alltogetherspectrumDc[i]->SetName(Form("AlltogetherSpectrum_%i",i));
                alltogetherspectrumDc[i]->Write();
              }
	    
	      TH1D *centrcrosscheck = correctedspectrum->Projection(0);
	      centrcrosscheck->SetTitle(Form("centrality_%i",i));
	      centrcrosscheck->SetName(Form("centrality_%i",i));
	      centrcrosscheck->Write();

	      TH1D *correctedTH1Dc = correctedspectrum->Projection(ptpr);
	      TH1D *alltogetherTH1Dc = (TH1D *) alltogetherCorrection->Project(ptpr);

	      TH1D *centrcrosscheck2 =  (TH1D *) alltogetherCorrection->Project(0);
	      centrcrosscheck2->SetTitle(Form("centrality2_%i",i));
	      centrcrosscheck2->SetName(Form("centrality2_%i",i));
	      centrcrosscheck2->Write();

	      TH1D* ratiocorrectedc = (TH1D*)correctedTH1D->Clone();
	      ratiocorrectedc->Divide(correctedTH1Dc,alltogetherTH1Dc,1,1);
	      ratiocorrectedc->SetTitle(Form("RatioUnfoldingAlltogetherSpectrum_%i",i));
	      ratiocorrectedc->SetName(Form("RatioUnfoldingAlltogetherSpectrum_%i",i));
              ratiocorrectedc->Write();

              fEfficiencyCharmSigD[i]->SetTitle(Form("IPEfficiencyForCharmSigCent%i",i));
	      fEfficiencyCharmSigD[i]->SetName(Form("IPEfficiencyForCharmSigCent%i",i));
              fEfficiencyCharmSigD[i]->Write();
	      fEfficiencyBeautySigD[i]->SetTitle(Form("IPEfficiencyForBeautySigCent%i",i));
	      fEfficiencyBeautySigD[i]->SetName(Form("IPEfficiencyForBeautySigCent%i",i));
              fEfficiencyBeautySigD[i]->Write();
	      fCharmEff[i]->SetTitle(Form("IPEfficiencyForCharmCent%i",i));
	      fCharmEff[i]->SetName(Form("IPEfficiencyForCharmCent%i",i));
              fCharmEff[i]->Write();
	      fBeautyEff[i]->SetTitle(Form("IPEfficiencyForBeautyCent%i",i));
	      fBeautyEff[i]->SetName(Form("IPEfficiencyForBeautyCent%i",i));
              fBeautyEff[i]->Write();
	      fConversionEff[i]->SetTitle(Form("IPEfficiencyForConversionCent%i",i));
	      fConversionEff[i]->SetName(Form("IPEfficiencyForConversionCent%i",i));
              fConversionEff[i]->Write();
	      fNonHFEEff[i]->SetTitle(Form("IPEfficiencyForNonHFECent%i",i));
	      fNonHFEEff[i]->SetName(Form("IPEfficiencyForNonHFECent%i",i));
              fNonHFEEff[i]->Write();
	  }

      }

      out->Close(); 
      delete out;
    }
  }

  return kTRUE;
}

//____________________________________________________________
AliCFDataGrid* AliHFEspectrum::SubtractBackground(Bool_t setBackground){
  //
  // Apply background subtraction
  //

    Int_t ptpr = 0;
    Int_t nbins=1;
    if(fBeamType==0)
    {
	ptpr=0;
        nbins=1;
    }
    if(fBeamType==1)
    {
	ptpr=1;
	nbins=2;
    }

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
 
  Int_t stepbackground = 1; // 2 for !fInclusiveSpectrum analysis(old method)
  AliCFDataGrid *backgroundGrid = new AliCFDataGrid("ContaminationGrid","ContaminationGrid",*backgroundContainer,stepbackground);

  if(!fInclusiveSpectrum){
    //Background subtraction for IP analysis

    TH1D *incElecCent[kCentrality-1];
    TH1D *charmCent[kCentrality-1];
    TH1D *convCent[kCentrality-1];
    TH1D *nonHFECent[kCentrality-1]; 
    TH1D *subtractedCent[kCentrality-1];

    for(Int_t initia=0; initia<(kCentrality-1); initia++) {
      incElecCent[initia]=0x0;
      charmCent[initia]=0x0;
      convCent[initia]=0x0;
      nonHFECent[initia]=0x0; 
      subtractedCent[initia]=0x0;
    }

    TH1D *measuredTH1Draw = (TH1D *) dataspectrumbeforesubstraction->Project(ptpr);
    CorrectFromTheWidth(measuredTH1Draw);
    if(fBeamType==1){
      THnSparseF* sparseIncElec = (THnSparseF *) dataspectrumbeforesubstraction->GetGrid();
      for(Int_t icent =  1; icent < kCentrality-1; icent++){
        sparseIncElec->GetAxis(0)->SetRange(icent,icent);
        incElecCent[icent-1] = (TH1D *) sparseIncElec->Projection(ptpr);
        CorrectFromTheWidth(incElecCent[icent-1]);
      }
    }
    TCanvas *rawspectra = new TCanvas("rawspectra","rawspectra",500,400);
    rawspectra->cd();
    rawspectra->SetLogy();
    gStyle->SetOptStat(0);
    TLegend *lRaw = new TLegend(0.55,0.55,0.85,0.85);
    measuredTH1Draw->SetMarkerStyle(20);
    measuredTH1Draw->Draw();
    measuredTH1Draw->GetXaxis()->SetRangeUser(0.0,7.9);
    lRaw->AddEntry(measuredTH1Draw,"measured raw spectrum");
    TH1D* htemp;
    Int_t* bins=new Int_t[2];
    if(fIPanaHadronBgSubtract){
      // Hadron background
	printf("Hadron background for IP analysis subtracted!\n");
	if(fBeamType==0)
	{

	    htemp  = (TH1D *) fHadronEffbyIPcut->Projection(0);
	    bins[0]=htemp->GetNbinsX();
	}
	if(fBeamType==1)
	{
	    htemp  = (TH1D *) fHadronEffbyIPcut->Projection(0);
	    bins[0]=htemp->GetNbinsX();
	    htemp  = (TH1D *) fHadronEffbyIPcut->Projection(1);
	    bins[1]=htemp->GetNbinsX();
	}
      AliCFDataGrid *hbgContainer = new AliCFDataGrid("hbgContainer","hadron bg after IP cut",nbins,bins);
      hbgContainer->SetGrid(fHadronEffbyIPcut);
      backgroundGrid->Multiply(hbgContainer,1);
      // draw raw hadron bg spectra
      TH1D *hadronbg= (TH1D *) backgroundGrid->Project(ptpr);
      CorrectFromTheWidth(hadronbg);
      hadronbg->SetMarkerColor(7);
      hadronbg->SetMarkerStyle(20);
      rawspectra->cd();
      hadronbg->Draw("samep");
      lRaw->AddEntry(hadronbg,"hadrons");
      // subtract hadron contamination
      spectrumSubtracted->Add(backgroundGrid,-1.0);
    }
    if(fIPanaCharmBgSubtract){
      // Charm background
      printf("Charm background for IP analysis subtracted!\n");
      AliCFDataGrid *charmbgContainer = (AliCFDataGrid *) GetCharmBackground();
      // draw charm bg spectra
      TH1D *charmbg= (TH1D *) charmbgContainer->Project(ptpr);
      CorrectFromTheWidth(charmbg);
      charmbg->SetMarkerColor(3);
      charmbg->SetMarkerStyle(20);
      rawspectra->cd();
      charmbg->Draw("samep");
      lRaw->AddEntry(charmbg,"charm elecs");
      // subtract charm background
      spectrumSubtracted->Add(charmbgContainer,-1.0);
      if(fBeamType==1){
        THnSparseF* sparseCharmElec = (THnSparseF *) charmbgContainer->GetGrid();
        for(Int_t icent =  1; icent < kCentrality-1; icent++){ 
          sparseCharmElec->GetAxis(0)->SetRange(icent,icent);
          charmCent[icent-1] = (TH1D *) sparseCharmElec->Projection(ptpr);
          CorrectFromTheWidth(charmCent[icent-1]);
        }
      }
    }
    if(fIPanaConversionBgSubtract){
      // Conversion background
      AliCFDataGrid *conversionbgContainer = (AliCFDataGrid *) GetConversionBackground(); 
      // draw conversion bg spectra
      TH1D *conversionbg= (TH1D *) conversionbgContainer->Project(ptpr);
      CorrectFromTheWidth(conversionbg);
      conversionbg->SetMarkerColor(4);
      conversionbg->SetMarkerStyle(20);
      rawspectra->cd();
      conversionbg->Draw("samep");
      lRaw->AddEntry(conversionbg,"conversion elecs");
      // subtract conversion background
      spectrumSubtracted->Add(conversionbgContainer,-1.0);
      if(fBeamType==1){
        THnSparseF* sparseconvElec = (THnSparseF *) conversionbgContainer->GetGrid();
        for(Int_t icent =  1; icent < kCentrality-1; icent++){
          sparseconvElec->GetAxis(0)->SetRange(icent,icent);
          convCent[icent-1] = (TH1D *) sparseconvElec->Projection(ptpr);
          CorrectFromTheWidth(convCent[icent-1]);
        }
      }
    }
    if(fIPanaNonHFEBgSubtract){
      // NonHFE background
      AliCFDataGrid *nonHFEbgContainer = (AliCFDataGrid *) GetNonHFEBackground();
      // draw Dalitz/dielectron bg spectra
      TH1D *nonhfebg= (TH1D *) nonHFEbgContainer->Project(ptpr);
      CorrectFromTheWidth(nonhfebg);
      nonhfebg->SetMarkerColor(6);
      nonhfebg->SetMarkerStyle(20);
      rawspectra->cd();
      nonhfebg->Draw("samep");
      lRaw->AddEntry(nonhfebg,"non-HF elecs");
      // subtract Dalitz/dielectron background
      spectrumSubtracted->Add(nonHFEbgContainer,-1.0);
      if(fBeamType==1){
        THnSparseF* sparseNonHFEElec = (THnSparseF *) nonHFEbgContainer->GetGrid();
        for(Int_t icent =  1; icent < kCentrality-1; icent++){
          sparseNonHFEElec->GetAxis(0)->SetRange(icent,icent);
          nonHFECent[icent-1] = (TH1D *) sparseNonHFEElec->Projection(ptpr);
          CorrectFromTheWidth(nonHFECent[icent-1]);
        }
      }
    }

    TH1D *rawbgsubtracted = (TH1D *) spectrumSubtracted->Project(ptpr);
    CorrectFromTheWidth(rawbgsubtracted);
    rawbgsubtracted->SetMarkerStyle(24);
    rawspectra->cd();
    lRaw->AddEntry(rawbgsubtracted,"subtracted raw spectrum");
    rawbgsubtracted->Draw("samep");
    lRaw->Draw("SAME");
    gPad->SetGrid();
    //rawspectra->SaveAs("rawspectra.eps");
    
    if(fBeamType==1){
      THnSparseF* sparseSubtracted = (THnSparseF *) spectrumSubtracted->GetGrid();
      for(Int_t icent =  1; icent < kCentrality-1; icent++){
        sparseSubtracted->GetAxis(0)->SetRange(icent,icent);
        subtractedCent[icent-1] = (TH1D *) sparseSubtracted->Projection(ptpr);
        CorrectFromTheWidth(subtractedCent[icent-1]);
      }
    
      TLegend *lCentRaw = new TLegend(0.55,0.55,0.85,0.85);
      TCanvas *centRaw = new TCanvas("centRaw","centRaw",1000,800);
      centRaw->Divide(3,3);
      for(Int_t icent = 1; icent < kCentrality-1; icent++){
        centRaw->cd(icent);
        gPad->SetLogx();
        gPad->SetLogy();
	if(!incElecCent[icent-1] || !charmCent[icent-1] || !convCent[icent-1] || !nonHFECent[icent-1] || !subtractedCent[icent-1]) continue;
        incElecCent[icent-1]->GetXaxis()->SetRangeUser(0.4,8.);
        incElecCent[icent-1]->Draw("p");
        incElecCent[icent-1]->SetMarkerColor(1);
        incElecCent[icent-1]->SetMarkerStyle(20);
        charmCent[icent-1]->Draw("samep");
        charmCent[icent-1]->SetMarkerColor(3);
        charmCent[icent-1]->SetMarkerStyle(20);
        convCent[icent-1]->Draw("samep");
        convCent[icent-1]->SetMarkerColor(4);
        convCent[icent-1]->SetMarkerStyle(20);
        nonHFECent[icent-1]->Draw("samep");
        nonHFECent[icent-1]->SetMarkerColor(6);
        nonHFECent[icent-1]->SetMarkerStyle(20);
        subtractedCent[icent-1]->Draw("samep");
        subtractedCent[icent-1]->SetMarkerStyle(24);
        if(icent == 1){
          lCentRaw->AddEntry(incElecCent[0],"inclusive electron spectrum");
          lCentRaw->AddEntry(charmCent[0],"charm elecs");
          lCentRaw->AddEntry(convCent[0],"conversion elecs");
          lCentRaw->AddEntry(nonHFECent[0],"non-HF elecs");
          lCentRaw->AddEntry(subtractedCent[0],"subtracted electron spectrum");
          lCentRaw->Draw("SAME");
        }
      }
    }

    delete[] bins; 

  }
  else{
    // Subtract 
    spectrumSubtracted->Add(backgroundGrid,-1.0);
  }

  if(setBackground){
    if(fBackground) delete fBackground;
    fBackground = backgroundGrid;
  } else delete backgroundGrid;


  if(fDebugLevel > 0) {

    Int_t ptprd=0;
    if(fBeamType==0) ptprd=0;
    if(fBeamType==1) ptprd=1;
    
    TCanvas * cbackgroundsubtraction = new TCanvas("backgroundsubtraction","backgroundsubtraction",1000,700);
    cbackgroundsubtraction->Divide(3,1);
    cbackgroundsubtraction->cd(1);
    gPad->SetLogy();
    TH1D *measuredTH1Daftersubstraction = (TH1D *) spectrumSubtracted->Project(ptprd);
    TH1D *measuredTH1Dbeforesubstraction = (TH1D *) dataspectrumbeforesubstraction->Project(ptprd);
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
    ratiomeasuredcontamination->Sumw2();
    ratiomeasuredcontamination->Add(measuredTH1Daftersubstraction,-1.0);
    ratiomeasuredcontamination->Divide(measuredTH1Dbeforesubstraction);
    ratiomeasuredcontamination->SetStats(0);
    ratiomeasuredcontamination->SetMarkerStyle(26);
    ratiomeasuredcontamination->SetMarkerColor(kBlack);
    ratiomeasuredcontamination->SetLineColor(kBlack);
    for(Int_t k=0; k < ratiomeasuredcontamination->GetNbinsX(); k++){
      ratiomeasuredcontamination->SetBinError(k+1,0.0);
    }
    ratiomeasuredcontamination->Draw("P");
    cbackgroundsubtraction->cd(3);
    TH1D *measuredTH1background = (TH1D *) backgroundGrid->Project(ptprd);
    CorrectFromTheWidth(measuredTH1background);
    measuredTH1background->SetStats(0);
    measuredTH1background->SetTitle("");
    measuredTH1background->GetYaxis()->SetTitle("dN/dp_{T} [(GeV/c)^{-1}]");
    measuredTH1background->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    measuredTH1background->SetMarkerStyle(26);
    measuredTH1background->SetMarkerColor(kBlack);
    measuredTH1background->SetLineColor(kBlack);
    measuredTH1background->Draw();
    if(fWriteToFile) cbackgroundsubtraction->SaveAs("BackgroundSubtracted.eps");

    if(fBeamType==1) {

      TCanvas * cbackgrounde = new TCanvas("BackgroundSubtraction_allspectra","BackgroundSubtraction_allspectra",1000,700);
      cbackgrounde->Divide(2,1);
      TLegend *legtotal = new TLegend(0.4,0.6,0.89,0.89);
      TLegend *legtotalg = new TLegend(0.4,0.6,0.89,0.89);
     
      THnSparseF* sparsesubtracted = (THnSparseF *) spectrumSubtracted->GetGrid();
      TAxis *cenaxisa = sparsesubtracted->GetAxis(0);
      THnSparseF* sparsebefore = (THnSparseF *) dataspectrumbeforesubstraction->GetGrid();
      TAxis *cenaxisb = sparsebefore->GetAxis(0);
      Int_t nbbin = cenaxisb->GetNbins();
      Int_t stylee[20] = {20,21,22,23,24,25,26,27,28,30,4,5,7,29,29,29,29,29,29,29};
      Int_t colorr[20] = {2,3,4,5,6,7,8,9,46,38,29,30,31,32,33,34,35,37,38,20};
      for(Int_t binc = 0; binc < nbbin; binc++){
        TString titlee("BackgroundSubtraction_centrality_bin_");
        titlee += binc;
        TCanvas * cbackground = new TCanvas((const char*) titlee,(const char*) titlee,1000,700);
        cbackground->Divide(2,1);
        cbackground->cd(1);
        gPad->SetLogy();
        cenaxisa->SetRange(binc+1,binc+1);
        cenaxisb->SetRange(binc+1,binc+1);
        TH1D *aftersubstraction = (TH1D *) sparsesubtracted->Projection(1);
        TH1D *beforesubstraction = (TH1D *) sparsebefore->Projection(1);
        CorrectFromTheWidth(aftersubstraction);
        CorrectFromTheWidth(beforesubstraction);
        aftersubstraction->SetStats(0);
        aftersubstraction->SetTitle((const char*)titlee);
        aftersubstraction->GetYaxis()->SetTitle("dN/dp_{T} [(GeV/c)^{-1}]");
        aftersubstraction->GetXaxis()->SetTitle("p_{T} [GeV/c]");
        aftersubstraction->SetMarkerStyle(25);
        aftersubstraction->SetMarkerColor(kBlack);
        aftersubstraction->SetLineColor(kBlack);
        beforesubstraction->SetStats(0);
        beforesubstraction->SetTitle((const char*)titlee);
        beforesubstraction->GetYaxis()->SetTitle("dN/dp_{T} [(GeV/c)^{-1}]");
        beforesubstraction->GetXaxis()->SetTitle("p_{T} [GeV/c]");
        beforesubstraction->SetMarkerStyle(24);
        beforesubstraction->SetMarkerColor(kBlue);
        beforesubstraction->SetLineColor(kBlue);
        aftersubstraction->Draw();
        beforesubstraction->Draw("same");
        TLegend *lega = new TLegend(0.4,0.6,0.89,0.89);
        lega->AddEntry(beforesubstraction,"With hadron contamination","p");
        lega->AddEntry(aftersubstraction,"Without hadron contamination ","p");
        lega->Draw("same");
        cbackgrounde->cd(1);
        gPad->SetLogy();
        TH1D *aftersubtractionn = (TH1D *) aftersubstraction->Clone();
        aftersubtractionn->SetMarkerStyle(stylee[binc]);
        aftersubtractionn->SetMarkerColor(colorr[binc]);
        if(binc==0) aftersubtractionn->Draw();
        else aftersubtractionn->Draw("same");
        legtotal->AddEntry(aftersubtractionn,(const char*) titlee,"p");
        cbackgrounde->cd(2);
        gPad->SetLogy();
        TH1D *aftersubtractionng = (TH1D *) aftersubstraction->Clone();
        aftersubtractionng->SetMarkerStyle(stylee[binc]);
        aftersubtractionng->SetMarkerColor(colorr[binc]);
        if(fNEvents[binc] > 0.0) aftersubtractionng->Scale(1/(Double_t)fNEvents[binc]);
        if(binc==0) aftersubtractionng->Draw();
        else aftersubtractionng->Draw("same");
        legtotalg->AddEntry(aftersubtractionng,(const char*) titlee,"p");
        cbackground->cd(2);
        TH1D* ratiocontamination = (TH1D*)beforesubstraction->Clone();
        ratiocontamination->SetName("ratiocontamination");
        ratiocontamination->SetTitle((const char*)titlee);
        ratiocontamination->GetYaxis()->SetTitle("(with contamination - without contamination) / with contamination");
        ratiocontamination->GetXaxis()->SetTitle("p_{T} [GeV/c]");
        ratiocontamination->Add(aftersubstraction,-1.0);
        ratiocontamination->Divide(beforesubstraction);
        Int_t totalbin = ratiocontamination->GetXaxis()->GetNbins();
        for(Int_t nbinpt = 0; nbinpt < totalbin; nbinpt++) {
          ratiocontamination->SetBinError(nbinpt+1,0.0);
        }
        ratiocontamination->SetStats(0);
        ratiocontamination->SetMarkerStyle(26);
        ratiocontamination->SetMarkerColor(kBlack);
        ratiocontamination->SetLineColor(kBlack);
        ratiocontamination->Draw("P");
      }
     
      cbackgrounde->cd(1);
      legtotal->Draw("same");
      cbackgrounde->cd(2);
      legtotalg->Draw("same");

      cenaxisa->SetRange(0,nbbin);
      cenaxisb->SetRange(0,nbbin);
      if(fWriteToFile) cbackgrounde->SaveAs("BackgroundSubtractedPbPb.eps");
    }
  }
  
  return spectrumSubtracted;
}

//____________________________________________________________
AliCFDataGrid* AliHFEspectrum::GetCharmBackground(){
  //
  // calculate charm background
  //
    Int_t ptpr = 0;
    Int_t nDim = 1;
    if(fBeamType==0)
    {
	ptpr=0;
    }
    if(fBeamType==1)
    {
	ptpr=1;
	nDim=2;
    }

  Double_t evtnorm=0;
  if(fNMCbgEvents[0]) evtnorm= double(fNEvents[0])/double(fNMCbgEvents[0]);

  AliCFContainer *mcContainer = GetContainer(kMCContainerCharmMC);
  if(!mcContainer){
    AliError("MC Container not available");
    return NULL;
  }

  if(!fCorrelation){
    AliError("No Correlation map available");
    return NULL;
  }

  AliCFDataGrid *charmBackgroundGrid= 0x0;
  charmBackgroundGrid = new AliCFDataGrid("charmBackgroundGrid","charmBackgroundGrid",*mcContainer, fStepMC-1); // use MC eff. up to right before PID

  TH1D *charmbgaftertofpid = (TH1D *) charmBackgroundGrid->Project(0);
  Int_t* bins=new Int_t[2];
  bins[0]=charmbgaftertofpid->GetNbinsX();

  if(fBeamType==1)
  {
      bins[0]=12;
      charmbgaftertofpid = (TH1D *) charmBackgroundGrid->Project(1);
      bins[1]=charmbgaftertofpid->GetNbinsX();

  }

  AliCFDataGrid *ipWeightedCharmContainer = new AliCFDataGrid("ipWeightedCharmContainer","ipWeightedCharmContainer",nDim,bins);
  ipWeightedCharmContainer->SetGrid(GetPIDxIPEff(0)); // get charm efficiency
  TH1D* parametrizedcharmpidipeff = (TH1D*)ipWeightedCharmContainer->Project(ptpr);

  charmBackgroundGrid->Multiply(ipWeightedCharmContainer,1.);

  Int_t *nBinpp=new Int_t[1];
  Int_t* binspp=new Int_t[1];
  binspp[0]=charmbgaftertofpid->GetNbinsX();  // number of pt bins

  Int_t *nBinPbPb=new Int_t[2];
  Int_t* binsPbPb=new Int_t[2];
  binsPbPb[1]=charmbgaftertofpid->GetNbinsX();  // number of pt bins
  binsPbPb[0]=12;

  Int_t looppt=binspp[0];
  if(fBeamType==1) looppt=binsPbPb[1];

  for(Long_t iBin=1; iBin<= looppt;iBin++){
      if(fBeamType==0)
      {
          nBinpp[0]=iBin;
          charmBackgroundGrid->SetElementError(nBinpp, charmBackgroundGrid->GetElementError(nBinpp)*evtnorm);
          charmBackgroundGrid->SetElement(nBinpp,charmBackgroundGrid->GetElement(nBinpp)*evtnorm);
      }
      if(fBeamType==1)
      {
          // loop over centrality
          for(Long_t iiBin=1; iiBin<= binsPbPb[0];iiBin++){
              nBinPbPb[0]=iiBin;
              nBinPbPb[1]=iBin;
              Double_t evtnormPbPb=0;
              if(fNMCbgEvents[iiBin-1]) evtnormPbPb= double(fNEvents[iiBin-1])/double(fNMCbgEvents[iiBin-1]);
              charmBackgroundGrid->SetElementError(nBinPbPb, charmBackgroundGrid->GetElementError(nBinPbPb)*evtnormPbPb);
              charmBackgroundGrid->SetElement(nBinPbPb,charmBackgroundGrid->GetElement(nBinPbPb)*evtnormPbPb);
          }
      }
  }

  TH1D* charmbgafteripcut = (TH1D*)charmBackgroundGrid->Project(ptpr);

  AliCFDataGrid *weightedCharmContainer = new AliCFDataGrid("weightedCharmContainer","weightedCharmContainer",nDim,bins);
  weightedCharmContainer->SetGrid(GetCharmWeights()); // get charm weighting factors
  TH1D* charmweightingfc = (TH1D*)weightedCharmContainer->Project(ptpr);
  charmBackgroundGrid->Multiply(weightedCharmContainer,1.);
  TH1D* charmbgafterweight = (TH1D*)charmBackgroundGrid->Project(ptpr);

  // Efficiency (set efficiency to 1 for only folding) 
  AliCFEffGrid* efficiencyD = new AliCFEffGrid("efficiency","",*mcContainer);
  efficiencyD->CalculateEfficiency(0,0);

  // Folding
  if(fBeamType==0)nDim = 1;
  if(fBeamType==1)nDim = 2;
  AliCFUnfolding folding("unfolding","",nDim,fCorrelation,efficiencyD->GetGrid(),charmBackgroundGrid->GetGrid(),charmBackgroundGrid->GetGrid());
  folding.SetMaxNumberOfIterations(1);
  folding.Unfold();

  // Results
  THnSparse* result1= folding.GetEstMeasured(); // folded spectra
  THnSparse* result=(THnSparse*)result1->Clone();
  charmBackgroundGrid->SetGrid(result);
  TH1D* charmbgafterfolding = (TH1D*)charmBackgroundGrid->Project(ptpr);

  //Charm background evaluation plots

  TCanvas *cCharmBgEval = new TCanvas("cCharmBgEval","cCharmBgEval",1000,600);
  cCharmBgEval->Divide(3,1);

  cCharmBgEval->cd(1);

  if(fBeamType==0)charmbgaftertofpid->Scale(evtnorm);
  if(fBeamType==1)
  {
      Double_t evtnormPbPb=0;
      for(Int_t kCentr=0;kCentr<bins[0];kCentr++)
      {
	  if(fNMCbgEvents[kCentr]) evtnormPbPb= evtnormPbPb+double(fNEvents[kCentr])/double(fNMCbgEvents[kCentr]);
      }
      charmbgaftertofpid->Scale(evtnormPbPb);
  }

  CorrectFromTheWidth(charmbgaftertofpid);
  charmbgaftertofpid->SetMarkerStyle(25);
  charmbgaftertofpid->Draw("p");
  charmbgaftertofpid->GetYaxis()->SetTitle("yield normalized by # of data events");
  charmbgaftertofpid->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  gPad->SetLogy();

  CorrectFromTheWidth(charmbgafteripcut);
  charmbgafteripcut->SetMarkerStyle(24);
  charmbgafteripcut->Draw("samep");

  CorrectFromTheWidth(charmbgafterweight);
  charmbgafterweight->SetMarkerStyle(24);
  charmbgafterweight->SetMarkerColor(4);
  charmbgafterweight->Draw("samep");

  CorrectFromTheWidth(charmbgafterfolding);
  charmbgafterfolding->SetMarkerStyle(24);
  charmbgafterfolding->SetMarkerColor(2);
  charmbgafterfolding->Draw("samep");

  cCharmBgEval->cd(2);
  parametrizedcharmpidipeff->SetMarkerStyle(24);
  parametrizedcharmpidipeff->Draw("p");
  parametrizedcharmpidipeff->GetXaxis()->SetTitle("p_{T} (GeV/c)");

  cCharmBgEval->cd(3);
  charmweightingfc->SetMarkerStyle(24);
  charmweightingfc->Draw("p");
  charmweightingfc->GetYaxis()->SetTitle("weighting factor for charm electron");
  charmweightingfc->GetXaxis()->SetTitle("p_{T} (GeV/c)");

  cCharmBgEval->cd(1);
  TLegend *legcharmbg = new TLegend(0.3,0.7,0.89,0.89);
  legcharmbg->AddEntry(charmbgaftertofpid,"After TOF PID","p");
  legcharmbg->AddEntry(charmbgafteripcut,"After IP cut","p");
  legcharmbg->AddEntry(charmbgafterweight,"After Weighting","p");
  legcharmbg->AddEntry(charmbgafterfolding,"After Folding","p");
  legcharmbg->Draw("same");

  cCharmBgEval->cd(2);
  TLegend *legcharmbg2 = new TLegend(0.3,0.7,0.89,0.89);
  legcharmbg2->AddEntry(parametrizedcharmpidipeff,"PID + IP cut eff.","p");
  legcharmbg2->Draw("same");

  CorrectStatErr(charmBackgroundGrid);
  if(fWriteToFile) cCharmBgEval->SaveAs("CharmBackground.eps");

  delete[] bins;
  delete[] nBinpp;
  delete[] binspp;
  delete[] nBinPbPb;
  delete[] binsPbPb;

  if(fUnfoldBG) UnfoldBG(charmBackgroundGrid);

  return charmBackgroundGrid;
}

//____________________________________________________________
AliCFDataGrid* AliHFEspectrum::GetConversionBackground(){
  //
  // calculate conversion background
  //
  
  Double_t evtnorm[1] = {0.0};
  if(fNMCbgEvents[0]) evtnorm[0]= double(fNEvents[0])/double(fNMCbgEvents[0]);
  printf("check event!!! %lf \n",evtnorm[0]);
  
  AliCFContainer *backgroundContainer = 0x0;
  
  if(fNonHFEsyst){
    backgroundContainer = (AliCFContainer*)fConvSourceContainer[0][0][0]->Clone();
    for(Int_t iSource = 1; iSource < kElecBgSources; iSource++){
      backgroundContainer->Add(fConvSourceContainer[iSource][0][0]); // make centrality dependent
    }  
  }
  else{    
    // Background Estimate
    backgroundContainer = GetContainer(kMCWeightedContainerConversionESD);    
  } 
  if(!backgroundContainer){
    AliError("MC background container not found");
    return NULL;
  }
  
  Int_t stepbackground = 3;

  AliCFDataGrid *backgroundGrid = new AliCFDataGrid("ConversionBgGrid","ConversionBgGrid",*backgroundContainer,stepbackground);
  Int_t *nBinpp=new Int_t[1];
  Int_t* binspp=new Int_t[1];
  binspp[0]=fConversionEff[0]->GetNbinsX();  // number of pt bins

  Int_t *nBinPbPb=new Int_t[2];
  Int_t* binsPbPb=new Int_t[2];
  binsPbPb[1]=fConversionEff[0]->GetNbinsX();  // number of pt bins
  binsPbPb[0]=12;

  Int_t looppt=binspp[0];
  if(fBeamType==1) looppt=binsPbPb[1];

  for(Long_t iBin=1; iBin<= looppt;iBin++){
      if(fBeamType==0)
      {
	  nBinpp[0]=iBin;
	  backgroundGrid->SetElementError(nBinpp, backgroundGrid->GetElementError(nBinpp)*evtnorm[0]);
	  backgroundGrid->SetElement(nBinpp,backgroundGrid->GetElement(nBinpp)*evtnorm[0]);
      }
      if(fBeamType==1)
      {
          // loop over centrality
	  for(Long_t iiBin=1; iiBin<= binsPbPb[0];iiBin++){
	      nBinPbPb[0]=iiBin;
	      nBinPbPb[1]=iBin;
              Double_t evtnormPbPb=0;
              if(fNMCbgEvents[iiBin-1]) evtnormPbPb= double(fNEvents[iiBin-1])/double(fNMCbgEvents[iiBin-1]);
	      backgroundGrid->SetElementError(nBinPbPb, backgroundGrid->GetElementError(nBinPbPb)*evtnormPbPb);
	      backgroundGrid->SetElement(nBinPbPb,backgroundGrid->GetElement(nBinPbPb)*evtnormPbPb);
	  }
      }
  }
  //end of workaround for statistical errors

  AliCFDataGrid *weightedConversionContainer;
  if(fBeamType==0) weightedConversionContainer = new AliCFDataGrid("weightedConversionContainer","weightedConversionContainer",1,binspp);
  else weightedConversionContainer = new AliCFDataGrid("weightedConversionContainer","weightedConversionContainer",2,binsPbPb);
  weightedConversionContainer->SetGrid(GetPIDxIPEff(2));
  backgroundGrid->Multiply(weightedConversionContainer,1.0);

  delete[] nBinpp;
  delete[] binspp;
  delete[] nBinPbPb;
  delete[] binsPbPb; 

  return backgroundGrid;
}


//____________________________________________________________
AliCFDataGrid* AliHFEspectrum::GetNonHFEBackground(){
  //
  // calculate non-HFE background
  //

  Double_t evtnorm[1] = {0.0};
  if(fNMCbgEvents[0]) evtnorm[0]= double(fNEvents[0])/double(fNMCbgEvents[0]);
  printf("check event!!! %lf \n",evtnorm[0]);
  
  AliCFContainer *backgroundContainer = 0x0;
  if(fNonHFEsyst){
    backgroundContainer = (AliCFContainer*)fNonHFESourceContainer[0][0][0]->Clone();
    for(Int_t iSource = 1; iSource < kElecBgSources; iSource++){
      if(iSource == 1)
        backgroundContainer->Add(fNonHFESourceContainer[iSource][0][0],1.41);//correction for the eta Dalitz decay branching ratio in PYTHIA
      else
        backgroundContainer->Add(fNonHFESourceContainer[iSource][0][0]);
    } 
  }
  else{    
    // Background Estimate 
    backgroundContainer = GetContainer(kMCWeightedContainerNonHFEESD);
  }
  if(!backgroundContainer){
    AliError("MC background container not found");
    return NULL;
  }
  
  
  Int_t stepbackground = 3;

  AliCFDataGrid *backgroundGrid = new AliCFDataGrid("NonHFEBgGrid","NonHFEBgGrid",*backgroundContainer,stepbackground);
  Int_t *nBinpp=new Int_t[1];
  Int_t* binspp=new Int_t[1];
  binspp[0]=fConversionEff[0]->GetNbinsX();  // number of pt bins

  Int_t *nBinPbPb=new Int_t[2];
  Int_t* binsPbPb=new Int_t[2];
  binsPbPb[1]=fConversionEff[0]->GetNbinsX();  // number of pt bins
  binsPbPb[0]=12;

  Int_t looppt=binspp[0];
  if(fBeamType==1) looppt=binsPbPb[1];


  for(Long_t iBin=1; iBin<= looppt;iBin++){
      if(fBeamType==0)
      {
	  nBinpp[0]=iBin;
	  backgroundGrid->SetElementError(nBinpp, backgroundGrid->GetElementError(nBinpp)*evtnorm[0]);
	  backgroundGrid->SetElement(nBinpp,backgroundGrid->GetElement(nBinpp)*evtnorm[0]);
      }
      if(fBeamType==1)
      {
	  for(Long_t iiBin=1; iiBin<=binsPbPb[0];iiBin++){
	      nBinPbPb[0]=iiBin;
	      nBinPbPb[1]=iBin;
	      Double_t evtnormPbPb=0;
              if(fNMCbgEvents[iiBin-1]) evtnormPbPb= double(fNEvents[iiBin-1])/double(fNMCbgEvents[iiBin-1]);
	      backgroundGrid->SetElementError(nBinPbPb, backgroundGrid->GetElementError(nBinPbPb)*evtnormPbPb);
	      backgroundGrid->SetElement(nBinPbPb,backgroundGrid->GetElement(nBinPbPb)*evtnormPbPb);
	  }
      }
  }
  //end of workaround for statistical errors
  AliCFDataGrid *weightedNonHFEContainer;
  if(fBeamType==0) weightedNonHFEContainer = new AliCFDataGrid("weightedNonHFEContainer","weightedNonHFEContainer",1,binspp);
  else weightedNonHFEContainer = new AliCFDataGrid("weightedNonHFEContainer","weightedNonHFEContainer",2,binsPbPb);
  weightedNonHFEContainer->SetGrid(GetPIDxIPEff(3));
  backgroundGrid->Multiply(weightedNonHFEContainer,1.0);

  delete[] nBinpp;
  delete[] binspp;
  delete[] nBinPbPb;
  delete[] binsPbPb; 

  return backgroundGrid;
}

//____________________________________________________________
AliCFDataGrid *AliHFEspectrum::CorrectParametrizedEfficiency(AliCFDataGrid* const bgsubpectrum){
  
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

  if(fDebugLevel > 0) {

    TCanvas * cEfficiencyParametrized = new TCanvas("EfficiencyParametrized","EfficiencyParametrized",1000,700);
    cEfficiencyParametrized->Divide(2,1);
    cEfficiencyParametrized->cd(1);
    TH1D *afterE = (TH1D *) resultt->Project(0);
    TH1D *beforeE = (TH1D *) dataGrid->Project(0);
    CorrectFromTheWidth(afterE);
    CorrectFromTheWidth(beforeE);
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
    gPad->SetLogy();
    afterE->Draw();
    beforeE->Draw("same");
    TLegend *legefficiencyparametrized = new TLegend(0.4,0.6,0.89,0.89);
    legefficiencyparametrized->AddEntry(beforeE,"Before Efficiency correction","p");
    legefficiencyparametrized->AddEntry(afterE,"After Efficiency correction","p");
    legefficiencyparametrized->Draw("same");
    cEfficiencyParametrized->cd(2);
    fEfficiencyFunction->Draw();
    //cEfficiencyParametrized->cd(3);
    //TH1D *ratioefficiency = (TH1D *) beforeE->Clone();
    //ratioefficiency->Divide(afterE);
    //ratioefficiency->Draw();

    if(fWriteToFile) cEfficiencyParametrized->SaveAs("efficiency.eps");

  }

  return resultt;

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

    Int_t ptpr=0;
    if(fBeamType==0) ptpr=0;
    if(fBeamType==1) ptpr=1;
    
    TCanvas * cV0Efficiency = new TCanvas("V0Efficiency","V0Efficiency",1000,700);
    cV0Efficiency->Divide(2,1);
    cV0Efficiency->cd(1);
    TH1D *afterE = (TH1D *) result->Project(ptpr);
    TH1D *beforeE = (TH1D *) dataGrid->Project(ptpr);
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
    TH1D* efficiencyDproj = (TH1D *) efficiencyD->Project(ptpr);
    efficiencyDproj->SetTitle("");
    efficiencyDproj->SetStats(0);
    efficiencyDproj->SetMarkerStyle(25);
    efficiencyDproj->Draw();

    if(fBeamType==1) {

      TCanvas * cV0Efficiencye = new TCanvas("V0Efficiency_allspectra","V0Efficiency_allspectra",1000,700);
      cV0Efficiencye->Divide(2,1);
      TLegend *legtotal = new TLegend(0.4,0.6,0.89,0.89);
      TLegend *legtotalg = new TLegend(0.4,0.6,0.89,0.89);
     
      THnSparseF* sparseafter = (THnSparseF *) result->GetGrid();
      TAxis *cenaxisa = sparseafter->GetAxis(0);
      THnSparseF* sparsebefore = (THnSparseF *) dataGrid->GetGrid();
      TAxis *cenaxisb = sparsebefore->GetAxis(0);
      THnSparseF* efficiencya = (THnSparseF *) efficiencyD->GetGrid();
      TAxis *cenaxisc = efficiencya->GetAxis(0);
      Int_t nbbin = cenaxisb->GetNbins();
      Int_t stylee[20] = {20,21,22,23,24,25,26,27,28,30,4,5,7,29,29,29,29,29,29,29};
      Int_t colorr[20] = {2,3,4,5,6,7,8,9,46,38,29,30,31,32,33,34,35,37,38,20};
      for(Int_t binc = 0; binc < nbbin; binc++){
        TString titlee("V0Efficiency_centrality_bin_");
        titlee += binc;
        TCanvas * ccV0Efficiency = new TCanvas((const char*) titlee,(const char*) titlee,1000,700);
        ccV0Efficiency->Divide(2,1);
        ccV0Efficiency->cd(1);
        gPad->SetLogy();
        cenaxisa->SetRange(binc+1,binc+1);
        cenaxisb->SetRange(binc+1,binc+1);
        cenaxisc->SetRange(binc+1,binc+1);
        TH1D *aftere = (TH1D *) sparseafter->Projection(1);
        TH1D *beforee = (TH1D *) sparsebefore->Projection(1);
        CorrectFromTheWidth(aftere);
        CorrectFromTheWidth(beforee);
        aftere->SetStats(0);
        aftere->SetTitle((const char*)titlee);
        aftere->GetYaxis()->SetTitle("dN/dp_{T} [(GeV/c)^{-1}]");
        aftere->GetXaxis()->SetTitle("p_{T} [GeV/c]");
        aftere->SetMarkerStyle(25);
        aftere->SetMarkerColor(kBlack);
        aftere->SetLineColor(kBlack);
        beforee->SetStats(0);
        beforee->SetTitle((const char*)titlee);
        beforee->GetYaxis()->SetTitle("dN/dp_{T} [(GeV/c)^{-1}]");
        beforee->GetXaxis()->SetTitle("p_{T} [GeV/c]");
        beforee->SetMarkerStyle(24);
        beforee->SetMarkerColor(kBlue);
        beforee->SetLineColor(kBlue);
        aftere->Draw();
        beforee->Draw("same");
        TLegend *lega = new TLegend(0.4,0.6,0.89,0.89);
        lega->AddEntry(beforee,"Before correction","p");
        lega->AddEntry(aftere,"After correction","p");
        lega->Draw("same");
        cV0Efficiencye->cd(1);
        gPad->SetLogy();
        TH1D *afteree = (TH1D *) aftere->Clone();
        afteree->SetMarkerStyle(stylee[binc]);
        afteree->SetMarkerColor(colorr[binc]);
        if(binc==0) afteree->Draw();
        else afteree->Draw("same");
        legtotal->AddEntry(afteree,(const char*) titlee,"p");
        cV0Efficiencye->cd(2);
        gPad->SetLogy();
        TH1D *aftereeu = (TH1D *) aftere->Clone();
        aftereeu->SetMarkerStyle(stylee[binc]);
        aftereeu->SetMarkerColor(colorr[binc]);
        if(fNEvents[binc] > 0.0) aftereeu->Scale(1/(Double_t)fNEvents[binc]);
        if(binc==0) aftereeu->Draw();
        else aftereeu->Draw("same");
        legtotalg->AddEntry(aftereeu,(const char*) titlee,"p");
        ccV0Efficiency->cd(2);
        TH1D* efficiencyDDproj = (TH1D *) efficiencya->Projection(1);
        efficiencyDDproj->SetTitle("");
        efficiencyDDproj->SetStats(0);
        efficiencyDDproj->SetMarkerStyle(25);
        efficiencyDDproj->Draw();
      }
     
      cV0Efficiencye->cd(1);
      legtotal->Draw("same");
      cV0Efficiencye->cd(2);
      legtotalg->Draw("same");

      cenaxisa->SetRange(0,nbbin);
      cenaxisb->SetRange(0,nbbin);
      cenaxisc->SetRange(0,nbbin);

      if(fWriteToFile) cV0Efficiencye->SaveAs("V0efficiency.eps");
    }

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

  if(!fBeauty2ndMethod)
  {
      // Consider parameterized IP cut efficiency
      if(!fInclusiveSpectrum){
	  Int_t* bins=new Int_t[1];
	  bins[0]=35;

	  AliCFEffGrid *beffContainer = new AliCFEffGrid("beffContainer","beffContainer",1,bins);
	  beffContainer->SetGrid(GetBeautyIPEff(kTRUE));
	  efficiencyD->Multiply(beffContainer,1);
      }
  }
  

  // Unfold 
  
  AliCFUnfolding unfolding("unfolding","",fNbDimensions,fCorrelation,efficiencyD->GetGrid(),dataGrid->GetGrid(),guessedTHnSparse);
  if(fUnSetCorrelatedErrors) unfolding.UnsetCorrelatedErrors();
  unfolding.SetMaxNumberOfIterations(fNumberOfIterations);
  if(fNRandomIter > 0) unfolding.SetNRandomIterations(fNRandomIter);
  if(fSetSmoothing) unfolding.UseSmoothing();
  unfolding.Unfold();

  // Results
  THnSparse* result = unfolding.GetUnfolded();
  THnSparse* residual = unfolding.GetEstMeasured();
  TList *listofresults = new TList;
  listofresults->SetOwner();
  listofresults->AddAt((THnSparse*)result->Clone(),0);
  listofresults->AddAt((THnSparse*)residual->Clone(),1);

  if(fDebugLevel > 0) {

    Int_t ptpr=0;
    if(fBeamType==0) ptpr=0;
    if(fBeamType==1) ptpr=1;
    
    TCanvas * cresidual = new TCanvas("residual","residual",1000,700);
    cresidual->Divide(2,1);
    cresidual->cd(1);
    gPad->SetLogy();
    TGraphErrors* residualspectrumD = Normalize(residual);
    if(!residualspectrumD) {
      AliError("Number of Events not set for the normalization");
      return NULL;
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
    TH1D *residualTH1D = residual->Projection(ptpr);
    TH1D *measuredTH1D = (TH1D *) dataGridBis->Project(ptpr);
    TH1D* ratioresidual = (TH1D*)residualTH1D->Clone();
    ratioresidual->SetName("ratioresidual");
    ratioresidual->SetTitle("");
    ratioresidual->GetYaxis()->SetTitle("Folded/Measured");
    ratioresidual->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    ratioresidual->Divide(residualTH1D,measuredTH1D,1,1);
    ratioresidual->SetStats(0);
    ratioresidual->Draw();

    if(fWriteToFile) cresidual->SaveAs("Unfolding.eps");
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


  if(!fBeauty2ndMethod)
  {
      // Consider parameterized IP cut efficiency
      if(!fInclusiveSpectrum){
	  Int_t* bins=new Int_t[1];
	  bins[0]=35;

	  AliCFEffGrid *beffContainer = new AliCFEffGrid("beffContainer","beffContainer",1,bins);
	  beffContainer->SetGrid(GetBeautyIPEff(kFALSE));
	  efficiencyD->Multiply(beffContainer,1);
      }
  }

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

    Int_t ptpr=0;
    if(fBeamType==0) ptpr=0;
    if(fBeamType==1) ptpr=1;
    
    printf("Step MC: %d\n",fStepMC);
    printf("Step tracking: %d\n",AliHFEcuts::kStepHFEcutsTRD + AliHFEcuts::kNcutStepsMCTrack);
    printf("Step MC true: %d\n",fStepTrue);
    AliCFEffGrid  *efficiencymcPID = (AliCFEffGrid*)  GetEfficiency(GetContainer(kMCContainerMC),fStepMC,AliHFEcuts::kStepHFEcutsTRD + AliHFEcuts::kNcutStepsMCTrack);
    AliCFEffGrid  *efficiencymctrackinggeo = (AliCFEffGrid*)  GetEfficiency(GetContainer(kMCContainerMC),AliHFEcuts::kStepHFEcutsTRD + AliHFEcuts::kNcutStepsMCTrack,fStepTrue);
    AliCFEffGrid  *efficiencymcall = (AliCFEffGrid*)  GetEfficiency(GetContainer(kMCContainerMC),fStepMC,fStepTrue);
    
    AliCFEffGrid  *efficiencyesdall = (AliCFEffGrid*)  GetEfficiency(GetContainer(kMCContainerESD),fStepMC,fStepTrue);
    
    TCanvas * cefficiency = new TCanvas("efficiency","efficiency",1000,700);
    cefficiency->cd(1);
    TH1D* efficiencymcPIDD = (TH1D *) efficiencymcPID->Project(ptpr);
    efficiencymcPIDD->SetTitle("");
    efficiencymcPIDD->SetStats(0);
    efficiencymcPIDD->SetMarkerStyle(25);
    efficiencymcPIDD->Draw();
    TH1D* efficiencymctrackinggeoD = (TH1D *) efficiencymctrackinggeo->Project(ptpr);
    efficiencymctrackinggeoD->SetTitle("");
    efficiencymctrackinggeoD->SetStats(0);
    efficiencymctrackinggeoD->SetMarkerStyle(26);
    efficiencymctrackinggeoD->Draw("same");
    TH1D* efficiencymcallD = (TH1D *) efficiencymcall->Project(ptpr);
    efficiencymcallD->SetTitle("");
    efficiencymcallD->SetStats(0);
    efficiencymcallD->SetMarkerStyle(27);
    efficiencymcallD->Draw("same");
    TH1D* efficiencyesdallD = (TH1D *) efficiencyesdall->Project(ptpr);
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

    if(fBeamType==1) {

      THnSparseF* sparseafter = (THnSparseF *) result->GetGrid();
      TAxis *cenaxisa = sparseafter->GetAxis(0);
      THnSparseF* sparsebefore = (THnSparseF *) dataGrid->GetGrid();
      TAxis *cenaxisb = sparsebefore->GetAxis(0);
      THnSparseF* efficiencya = (THnSparseF *) efficiencyD->GetGrid();
      TAxis *cenaxisc = efficiencya->GetAxis(0);
      //Int_t nbbin = cenaxisb->GetNbins();
      //Int_t stylee[20] = {20,21,22,23,24,25,26,27,28,30,4,5,7,29,29,29,29,29,29,29};
      //Int_t colorr[20] = {2,3,4,5,6,7,8,9,46,38,29,30,31,32,33,34,35,37,38,20};
      for(Int_t binc = 0; binc < fNCentralityBinAtTheEnd; binc++){
        TString titlee("Efficiency_centrality_bin_");
        titlee += fLowBoundaryCentralityBinAtTheEnd[binc];
        titlee += "_";
        titlee += fHighBoundaryCentralityBinAtTheEnd[binc];
        TCanvas * cefficiencye = new TCanvas((const char*) titlee,(const char*) titlee,1000,700);
        cefficiencye->Divide(2,1);
        cefficiencye->cd(1);
        gPad->SetLogy();
        cenaxisa->SetRange(fLowBoundaryCentralityBinAtTheEnd[binc]+1,fHighBoundaryCentralityBinAtTheEnd[binc]);
        cenaxisb->SetRange(fLowBoundaryCentralityBinAtTheEnd[binc]+1,fHighBoundaryCentralityBinAtTheEnd[binc]);
        TH1D *afterefficiency = (TH1D *) sparseafter->Projection(ptpr);
        TH1D *beforeefficiency = (TH1D *) sparsebefore->Projection(ptpr);
        CorrectFromTheWidth(afterefficiency);
        CorrectFromTheWidth(beforeefficiency);
        afterefficiency->SetStats(0);
        afterefficiency->SetTitle((const char*)titlee);
        afterefficiency->GetYaxis()->SetTitle("dN/dp_{T} [(GeV/c)^{-1}]");
        afterefficiency->GetXaxis()->SetTitle("p_{T} [GeV/c]");
        afterefficiency->SetMarkerStyle(25);
        afterefficiency->SetMarkerColor(kBlack);
        afterefficiency->SetLineColor(kBlack);
        beforeefficiency->SetStats(0);
        beforeefficiency->SetTitle((const char*)titlee);
        beforeefficiency->GetYaxis()->SetTitle("dN/dp_{T} [(GeV/c)^{-1}]");
        beforeefficiency->GetXaxis()->SetTitle("p_{T} [GeV/c]");
        beforeefficiency->SetMarkerStyle(24);
        beforeefficiency->SetMarkerColor(kBlue);
        beforeefficiency->SetLineColor(kBlue);
        afterefficiency->Draw();
        beforeefficiency->Draw("same");
        TLegend *lega = new TLegend(0.4,0.6,0.89,0.89);
        lega->AddEntry(beforeefficiency,"Before efficiency correction","p");
        lega->AddEntry(afterefficiency,"After efficiency correction","p");
        lega->Draw("same");
        cefficiencye->cd(2);
        cenaxisc->SetRange(fLowBoundaryCentralityBinAtTheEnd[binc]+1,fLowBoundaryCentralityBinAtTheEnd[binc]+1);
        TH1D* efficiencyDDproj = (TH1D *) efficiencya->Projection(ptpr);
        efficiencyDDproj->SetTitle("");
        efficiencyDDproj->SetStats(0);
        efficiencyDDproj->SetMarkerStyle(25);
        efficiencyDDproj->SetMarkerColor(2);
        efficiencyDDproj->Draw();
        cenaxisc->SetRange(fHighBoundaryCentralityBinAtTheEnd[binc],fHighBoundaryCentralityBinAtTheEnd[binc]);
        TH1D* efficiencyDDproja = (TH1D *) efficiencya->Projection(ptpr);
        efficiencyDDproja->SetTitle("");
        efficiencyDDproja->SetStats(0);
        efficiencyDDproja->SetMarkerStyle(26);
        efficiencyDDproja->SetMarkerColor(4);
        efficiencyDDproja->Draw("same");
      }
    }

    if(fWriteToFile) cefficiency->SaveAs("efficiencyCorrected.eps");
  }
  
  return result;

}

//__________________________________________________________________________________
TGraphErrors *AliHFEspectrum::Normalize(THnSparse * const spectrum,Int_t i) const {
  //
  // Normalize the spectrum to 1/(2*Pi*p_{T})*dN/dp_{T} (GeV/c)^{-2}
  // Give the final pt spectrum to be compared
  //

  if(fNEvents[i] > 0) {

    Int_t ptpr = 0;
    if(fBeamType==0) ptpr=0;
    if(fBeamType==1) ptpr=1;
    
    TH1D* projection = spectrum->Projection(ptpr);
    CorrectFromTheWidth(projection);
    TGraphErrors *graphError = NormalizeTH1(projection,i);
    return graphError;
  
  }
    
  return 0x0;
  

}
//__________________________________________________________________________________
TGraphErrors *AliHFEspectrum::Normalize(AliCFDataGrid * const spectrum,Int_t i) const {
  //
  // Normalize the spectrum to 1/(2*Pi*p_{T})*dN/dp_{T} (GeV/c)^{-2}
  // Give the final pt spectrum to be compared
  //
  if(fNEvents[i] > 0) {

    Int_t ptpr=0;
    if(fBeamType==0) ptpr=0;
    if(fBeamType==1) ptpr=1;
    
    TH1D* projection = (TH1D *) spectrum->Project(ptpr);
    CorrectFromTheWidth(projection);
    TGraphErrors *graphError = NormalizeTH1(projection,i);

    return graphError;
    
  }

  return 0x0;
  

}
//__________________________________________________________________________________
TGraphErrors *AliHFEspectrum::NormalizeTH1(TH1 *input,Int_t i) const {
  //
  // Normalize the spectrum to 1/(2*Pi*p_{T})*dN/dp_{T} (GeV/c)^{-2}
  // Give the final pt spectrum to be compared
  //
  Double_t chargecoefficient = 0.5;
  if(fChargeChoosen != 0) chargecoefficient = 1.0;

  Double_t etarange = fEtaSelected ? fEtaRangeNorm[1] - fEtaRangeNorm[0] : 1.6;
  printf("Normalizing Eta Range %f\n", etarange);
  if(fNEvents[i] > 0) {

    TGraphErrors *spectrumNormalized = new TGraphErrors(input->GetNbinsX());
    Double_t p = 0, dp = 0; Int_t point = 1;
    Double_t n = 0, dN = 0;
    Double_t nCorr = 0, dNcorr = 0;
    //Double_t errdN = 0, errdp = 0;
    Double_t errdN = 0;
    for(Int_t ibin = input->GetXaxis()->GetFirst(); ibin <= input->GetXaxis()->GetLast(); ibin++){
      point = ibin - input->GetXaxis()->GetFirst();
      p = input->GetXaxis()->GetBinCenter(ibin);
      dp = input->GetXaxis()->GetBinWidth(ibin)/2.;
      n = input->GetBinContent(ibin);
      dN = input->GetBinError(ibin);
      // New point
      nCorr = chargecoefficient * 1./etarange * 1./(Double_t)(fNEvents[i]) * 1./(2. * TMath::Pi() * p) * n;
      errdN = 1./(2. * TMath::Pi() * p);
      //errdp = 1./(2. * TMath::Pi() * p*p) * n;
      //dNcorr = chargecoefficient * 1./etarange * 1./(Double_t)(fNEvents[i]) * TMath::Sqrt(errdN * errdN * dN *dN + errdp *errdp * dp *dp);
      dNcorr = chargecoefficient * 1./etarange * 1./(Double_t)(fNEvents[i]) * TMath::Sqrt(errdN * errdN * dN *dN);
      
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
//__________________________________________________________________________________
TGraphErrors *AliHFEspectrum::NormalizeTH1N(TH1 *input,Int_t normalization) const {
  //
  // Normalize the spectrum to 1/(2*Pi*p_{T})*dN/dp_{T} (GeV/c)^{-2}
  // Give the final pt spectrum to be compared
  //
  Double_t chargecoefficient = 0.5;
  if(fChargeChoosen != kAllCharge) chargecoefficient = 1.0;

  Double_t etarange = fEtaSelected ? fEtaRangeNorm[1] - fEtaRangeNorm[0] : 1.6;
  printf("Normalizing Eta Range %f\n", etarange);
  if(normalization > 0) {

    TGraphErrors *spectrumNormalized = new TGraphErrors(input->GetNbinsX());
    Double_t p = 0, dp = 0; Int_t point = 1;
    Double_t n = 0, dN = 0;
    Double_t nCorr = 0, dNcorr = 0;
    //Double_t errdN = 0, errdp = 0;
    Double_t errdN = 0;
    for(Int_t ibin = input->GetXaxis()->GetFirst(); ibin <= input->GetXaxis()->GetLast(); ibin++){
      point = ibin - input->GetXaxis()->GetFirst();
      p = input->GetXaxis()->GetBinCenter(ibin);
      //dp = input->GetXaxis()->GetBinWidth(ibin)/2.;
      n = input->GetBinContent(ibin);
      dN = input->GetBinError(ibin);
      // New point
      nCorr = chargecoefficient * 1./etarange * 1./(Double_t)(normalization) * 1./(2. * TMath::Pi() * p) * n;
      errdN = 1./(2. * TMath::Pi() * p);
      //errdp = 1./(2. * TMath::Pi() * p*p) * n;
      //dNcorr = chargecoefficient * 1./etarange * 1./(Double_t)(normalization) * TMath::Sqrt(errdN * errdN * dN *dN + errdp *errdp * dp *dp);
      dNcorr = chargecoefficient * 1./etarange * 1./(Double_t)(normalization) * TMath::Sqrt(errdN * errdN * dN *dN);
      
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
  if(!fCFContainers) fCFContainers = new TObjArray(kDataContainerV0+1);
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
AliCFContainer *AliHFEspectrum::GetSlicedContainer(AliCFContainer *container, Int_t nDim, Int_t *dimensions,Int_t source,Chargetype_t charge, Int_t centralitylow, Int_t centralityhigh) {
  //
  // Slice bin for a given source of electron
  // nDim is the number of dimension the corrections are done
  // dimensions are the definition of the dimensions
  // source is if we want to keep only one MC source (-1 means we don't cut on the MC source)
  // positivenegative if we want to keep positive (1) or negative (0) or both (-1)
  // centrality (-1 means we do not cut on centrality)
  //
  
  Double_t *varMin = new Double_t[container->GetNVar()],
           *varMax = new Double_t[container->GetNVar()];

  Double_t *binLimits;
  for(Int_t ivar = 0; ivar < container->GetNVar(); ivar++){
    
    binLimits = new Double_t[container->GetNBins(ivar)+1];
    container->GetBinLimits(ivar,binLimits);
    varMin[ivar] = binLimits[0];
    varMax[ivar] = binLimits[container->GetNBins(ivar)];
    // source
    if(ivar == 4){
      if((source>= 0) && (source<container->GetNBins(ivar))) {
              varMin[ivar] = container->GetAxis(4,0)->GetBinLowEdge(container->GetAxis(4,0)->FindBin(binLimits[source]));
              varMax[ivar] = container->GetAxis(4,0)->GetBinUpEdge(container->GetAxis(4,0)->FindBin(binLimits[source]));
      }     
    }
    // charge
    if(ivar == 3) {
      if(charge != kAllCharge){
        varMin[ivar] = container->GetAxis(3,0)->GetBinLowEdge(container->GetAxis(3,0)->FindBin(charge));
        varMax[ivar] = container->GetAxis(3,0)->GetBinUpEdge(container->GetAxis(3,0)->FindBin(charge));
      }
    }
    // eta
    if(ivar == 1){
      for(Int_t ic = 1; ic <= container->GetAxis(1,0)->GetLast(); ic++) 
        AliDebug(1, Form("eta bin %d, min %f, max %f\n", ic, container->GetAxis(1,0)->GetBinLowEdge(ic), container->GetAxis(1,0)->GetBinUpEdge(ic))); 
      if(fEtaSelected){
        varMin[ivar] = fEtaRange[0];
        varMax[ivar] = fEtaRange[1];
      }
    }
    if(fEtaSelected){
      fEtaRangeNorm[0] = container->GetAxis(1,0)->GetBinLowEdge(container->GetAxis(1,0)->FindBin(fEtaRange[0]));
      fEtaRangeNorm[1] = container->GetAxis(1,0)->GetBinUpEdge(container->GetAxis(1,0)->FindBin(fEtaRange[1]));
      AliInfo(Form("Normalization done in eta range [%f,%f]\n", fEtaRangeNorm[0], fEtaRangeNorm[0]));
    }
    if(ivar == 5){
	if((centralitylow>= 0) && (centralitylow<container->GetNBins(ivar)) && (centralityhigh>= 0) && (centralityhigh<container->GetNBins(ivar))) {
	    varMin[ivar] = binLimits[centralitylow];
	    varMax[ivar] = binLimits[centralityhigh];

	    TAxis *axistest = container->GetAxis(5,0);
	    printf("GetSlicedContainer: Number of bin in centrality direction %d\n",axistest->GetNbins());
	    printf("Project from %f to %f\n",binLimits[centralitylow],binLimits[centralityhigh]);
	    Double_t lowcentrality = axistest->GetBinLowEdge(axistest->FindBin(binLimits[centralitylow]));
	    Double_t highcentrality = axistest->GetBinUpEdge(axistest->FindBin(binLimits[centralityhigh]));
	    printf("GetSlicedContainer: Low centrality %f and high centrality %f\n",lowcentrality,highcentrality);
	
	}

    }


    delete[] binLimits;
  }
  
  AliCFContainer *k = container->MakeSlice(nDim, dimensions, varMin, varMax);
  AddTemporaryObject(k);
  delete[] varMin; delete[] varMax;

  return k;

}

//_________________________________________________________________________
THnSparseF *AliHFEspectrum::GetSlicedCorrelation(THnSparseF *correlationmatrix, Int_t nDim, Int_t *dimensions, Int_t centralitylow, Int_t centralityhigh) const {
  //
  // Slice correlation
  //

  Int_t ndimensions = correlationmatrix->GetNdimensions();
  //printf("Number of dimension %d correlation map\n",ndimensions);
  if(ndimensions < (2*nDim)) {
    AliError("Problem in the dimensions");
    return NULL;
  }
  
  // Cut in centrality is centrality > -1
  if((centralitylow >=0) && (centralityhigh >=0)) {

    TAxis *axiscentrality0 = correlationmatrix->GetAxis(5);
    TAxis *axiscentrality1 = correlationmatrix->GetAxis(5+((Int_t)(ndimensions/2.)));

    Int_t bins0 = axiscentrality0->GetNbins();
    Int_t bins1 = axiscentrality1->GetNbins();
    
    printf("Number of centrality bins: %d and %d\n",bins0,bins1);
    if(bins0 != bins1) {
      AliError("Problem in the dimensions");
      return NULL;
    }
    
    if((centralitylow>= 0) && (centralitylow<bins0) && (centralityhigh>= 0) && (centralityhigh<bins0)) {
      axiscentrality0->SetRangeUser(centralitylow,centralityhigh);
      axiscentrality1->SetRangeUser(centralitylow,centralityhigh);

      Double_t lowcentrality0 = axiscentrality0->GetBinLowEdge(axiscentrality0->FindBin(centralitylow));
      Double_t highcentrality0 = axiscentrality0->GetBinUpEdge(axiscentrality0->FindBin(centralityhigh));
      Double_t lowcentrality1 = axiscentrality1->GetBinLowEdge(axiscentrality1->FindBin(centralitylow));
      Double_t highcentrality1 = axiscentrality1->GetBinUpEdge(axiscentrality1->FindBin(centralityhigh));
      printf("GetCorrelation0: Low centrality %f and high centrality %f\n",lowcentrality0,highcentrality0);
      printf("GetCorrelation1: Low centrality %f and high centrality %f\n",lowcentrality1,highcentrality1);

      TCanvas * ctestcorrelation = new TCanvas("testcorrelationprojection","testcorrelationprojection",1000,700);
      ctestcorrelation->cd(1);
      TH2D* projection = (TH2D *) correlationmatrix->Projection(5,5+((Int_t)(ndimensions/2.)));
      projection->Draw("colz");

    }
    
  }


  Int_t ndimensionsContainer = (Int_t) ndimensions/2;
  //printf("Number of dimension %d container\n",ndimensionsContainer);

  Int_t *dim = new Int_t[nDim*2];
  for(Int_t iter=0; iter < nDim; iter++){
    dim[iter] = dimensions[iter];
    dim[iter+nDim] = ndimensionsContainer + dimensions[iter];
    //printf("For iter %d: %d and iter+nDim %d: %d\n",iter,dimensions[iter],iter+nDim,ndimensionsContainer + dimensions[iter]);
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

//___________________________________________________________________________
void AliHFEspectrum::CorrectStatErr(AliCFDataGrid *backgroundGrid) const { 
  //
  // Correct statistical error
  //

  TH1D *h1 = (TH1D*)backgroundGrid->Project(0);
  Int_t nbinX = h1->GetNbinsX();
  Int_t bins[1];
  for(Long_t i = 1; i <= nbinX; i++) {
    bins[0] = i;
    Float_t content = h1->GetBinContent(i);
    if(content>0){
      Float_t error = TMath::Sqrt(content);
      backgroundGrid->SetElementError(bins, error);
    }
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
TObject* AliHFEspectrum::GetSpectrum(const AliCFContainer * const c, Int_t step) {
  AliCFDataGrid* data = new AliCFDataGrid("data","",*c, step);
  return data;
}
//_________________________________________________________________________
TObject* AliHFEspectrum::GetEfficiency(const AliCFContainer * const c, Int_t step, Int_t step0){
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

//____________________________________________________________________________
THnSparse* AliHFEspectrum::GetCharmWeights(){
  
  //
  // Measured D->e based weighting factors
  //

  const Int_t nDimpp=1;
  Int_t nBinpp[nDimpp] = {35};
  Double_t ptbinning1[36] = {0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.75, 2., 2.25, 2.5, 2.75, 3., 3.5, 4., 4.5, 5., 5.5, 6., 7., 8., 10., 12., 14., 16., 18., 20.};
  const Int_t nDimPbPb=2;
  Int_t nBinPbPb[nDimPbPb] = {11,35};
  Double_t kCentralityRange[12] = {0.,1.,2., 3., 4., 5., 6., 7.,8.,9., 10., 11.};
  Int_t loopcentr=1;
  Int_t looppt=nBinpp[0];
  if(fBeamType==0)
  {
      fWeightCharm = new THnSparseF("weightHisto", "weighting factor; pt[GeV/c]", nDimpp, nBinpp);
      fWeightCharm->SetBinEdges(0, ptbinning1);
  }
  if(fBeamType==1)
  {
      fWeightCharm = new THnSparseF("weightHisto", "weighting factor; centrality bin; pt[GeV/c]", nDimPbPb, nBinPbPb);
      fWeightCharm->SetBinEdges(1, ptbinning1);
      fWeightCharm->SetBinEdges(0, kCentralityRange);
      loopcentr=nBinPbPb[0];
  }

  // Weighting factor for pp
  Double_t weight[35]={0.859260, 0.872552, 0.847475, 0.823631, 0.839386, 0.874024, 0.916755, 0.942801, 0.965856, 0.933905, 0.933414, 0.931936, 0.847826, 0.810902, 0.796608, 0.727002, 0.659227, 0.583610, 0.549956, 0.512633, 0.472254, 0.412364, 0.353191, 0.319145, 0.305412, 0.290334, 0.269863, 0.254646, 0.230245, 0.200859, 0.275953, 0.276271, 0.227332, 0.197004, 0.474385};
  
  // Weighting factor for PbPb (0-20%)
  //Double_t weight[35]={0.641897,  0.640472,  0.615228,  0.650469,  0.737762,  0.847867,  1.009317,  1.158594,  1.307482,  1.476973,  1.551131,  1.677131,  1.785478,  1.888933,  2.017957,  2.074757,  1.926700,  1.869495,  1.546558,  1.222873,  1.160313,  0.903375,  0.799642,  0.706244,  0.705449,  0.599947,  0.719570,  0.499422,  0.703978,  0.477452,  0.325057,  0.093391,  0.096675,  0.000000,  0.000000};

  // Weighting factor for PbPb (40-80%)
  //Double_t weight[35]={0.181953,  0.173248,  0.166799,  0.182558,  0.206581,  0.236955,  0.279390,  0.329129,  0.365260,  0.423059,  0.452057,  0.482726,  0.462627,  0.537770,  0.584663,  0.579452,  0.587194,  0.499498,  0.443299,  0.398596,  0.376695,  0.322331,  0.260890,  0.374834,  0.249114,  0.310330,  0.287326,  0.243174,  0.758945,  0.138867,  0.170576,  0.107797,  0.011390,  0.000000,  0.000000};

  //points
  Double_t pt[1];
  Double_t contents[2];

  for(int icentr=0; icentr<loopcentr; icentr++)
  {
      for(int i=0; i<looppt; i++){
	  pt[0]=(ptbinning1[i]+ptbinning1[i+1])/2.;
	  if(fBeamType==1)
	  {
	      contents[0]=icentr;
	      contents[1]=pt[0];
	  }
	  if(fBeamType==0)
	  {
	      contents[0]=pt[0];
	  }

	  fWeightCharm->Fill(contents,weight[i]);
      }
  }

  Int_t nDimSparse = fWeightCharm->GetNdimensions();
  Int_t* binsvar = new Int_t[nDimSparse]; // number of bins for each variable
  Long_t nBins = 1; // used to calculate the total number of bins in the THnSparse

  for (Int_t iVar=0; iVar<nDimSparse; iVar++) {
      binsvar[iVar] = fWeightCharm->GetAxis(iVar)->GetNbins();
      nBins *= binsvar[iVar];
  }

  Int_t *binfill = new Int_t[nDimSparse]; // bin to fill the THnSparse (holding the bin coordinates)
  // loop that sets 0 error in each bin
  for (Long_t iBin=0; iBin<nBins; iBin++) {
    Long_t bintmp = iBin ;
    for (Int_t iVar=0; iVar<nDimSparse; iVar++) {
      binfill[iVar] = 1 + bintmp % binsvar[iVar] ;
      bintmp /= binsvar[iVar] ;
    }
    fWeightCharm->SetBinError(binfill,0.); // put 0 everywhere
  }

  delete[] binsvar;
  delete[] binfill;

  return fWeightCharm;
}

//____________________________________________________________________________
void AliHFEspectrum::SetParameterizedEff(AliCFContainer *container, AliCFContainer *containermb, AliCFContainer *containeresd, AliCFContainer *containeresdmb, Int_t *dimensions){

   // TOF PID efficiencies
   Int_t ptpr=0;
   if(fBeamType==0) ptpr=0;
   if(fBeamType==1) ptpr=1;

   Int_t loopcentr=1;
   const Int_t nCentralitybinning=11; //number of centrality bins
   if(fBeamType==1)
   {
     loopcentr=nCentralitybinning;
   }

   TF1 *fittofpid = new TF1("fittofpid","[0]*([1]+[2]*log(x)+[3]*log(x)*log(x))*tanh([4]*x-[5])",0.9,8.);
   TF1 *fipfit = new TF1("fipfit","[0]*([1]+[2]*log(x)+[3]*log(x)*log(x))*tanh([4]*x-[5])",0.9,8.);
   TF1 *fipfitnonhfe = new TF1("fipfitnonhfe","[0]*([1]+[2]*log(x)+[3]*log(x)*log(x))*tanh([4]*x-[5])",0.3,10.0);

   TCanvas * cefficiencyParamtof = new TCanvas("efficiencyParamtof","efficiencyParamtof",600,600);
   cefficiencyParamtof->cd();

   AliCFContainer *mccontainermcD = 0x0;
   AliCFContainer *mccontaineresdD = 0x0;
   TH1D* efficiencysigTOFPIDD[nCentralitybinning];
   TH1D* efficiencyTOFPIDD[nCentralitybinning];
   TH1D* efficiencysigesdTOFPIDD[nCentralitybinning];
   TH1D* efficiencyesdTOFPIDD[nCentralitybinning];
   Int_t source = -1; //get parameterized TOF PID efficiencies

   for(int icentr=0; icentr<loopcentr; icentr++) {
      // signal sample
      if(fBeamType==0) mccontainermcD = GetSlicedContainer(container, fNbDimensions, dimensions, source, fChargeChoosen);
      else mccontainermcD = GetSlicedContainer(container, fNbDimensions, dimensions, source, fChargeChoosen,icentr+1);
      AliCFEffGrid* efficiencymcsigParamTOFPID= new AliCFEffGrid("efficiencymcsigParamTOFPID","",*mccontainermcD);
      efficiencymcsigParamTOFPID->CalculateEfficiency(fStepMC,fStepMC-1); // TOF PID efficiencies

      // mb sample for double check
      if(fBeamType==0)  mccontainermcD = GetSlicedContainer(containermb, fNbDimensions, dimensions, source, fChargeChoosen);
      else mccontainermcD = GetSlicedContainer(containermb, fNbDimensions, dimensions, source, fChargeChoosen,icentr+1);
      AliCFEffGrid* efficiencymcParamTOFPID= new AliCFEffGrid("efficiencymcParamTOFPID","",*mccontainermcD);
      efficiencymcParamTOFPID->CalculateEfficiency(fStepMC,fStepMC-1); // TOF PID efficiencies

      // mb sample with reconstructed variables
      if(fBeamType==0)  mccontainermcD = GetSlicedContainer(containeresdmb, fNbDimensions, dimensions, source, fChargeChoosen);
      else mccontainermcD = GetSlicedContainer(containeresdmb, fNbDimensions, dimensions, source, fChargeChoosen,icentr+1);
      AliCFEffGrid* efficiencyesdParamTOFPID= new AliCFEffGrid("efficiencyesdParamTOFPID","",*mccontainermcD);
      efficiencyesdParamTOFPID->CalculateEfficiency(fStepMC,fStepMC-1); // TOF PID efficiencies

      // mb sample with reconstructed variables
      if(fBeamType==0)  mccontainermcD = GetSlicedContainer(containeresd, fNbDimensions, dimensions, source, fChargeChoosen);
      else mccontainermcD = GetSlicedContainer(containeresd, fNbDimensions, dimensions, source, fChargeChoosen,icentr+1);
      AliCFEffGrid* efficiencysigesdParamTOFPID= new AliCFEffGrid("efficiencysigesdParamTOFPID","",*mccontainermcD);
      efficiencysigesdParamTOFPID->CalculateEfficiency(fStepMC,fStepMC-1); // TOF PID efficiencies

      //fill histo
      efficiencysigTOFPIDD[icentr] = (TH1D *) efficiencymcsigParamTOFPID->Project(ptpr);
      efficiencyTOFPIDD[icentr] = (TH1D *) efficiencymcParamTOFPID->Project(ptpr);
      efficiencysigesdTOFPIDD[icentr] = (TH1D *) efficiencysigesdParamTOFPID->Project(ptpr);
      efficiencyesdTOFPIDD[icentr] = (TH1D *) efficiencyesdParamTOFPID->Project(ptpr);
      efficiencysigTOFPIDD[icentr]->SetName(Form("efficiencysigTOFPIDD%d",icentr));
      efficiencyTOFPIDD[icentr]->SetName(Form("efficiencyTOFPIDD%d",icentr));
      efficiencysigesdTOFPIDD[icentr]->SetName(Form("efficiencysigesdTOFPIDD%d",icentr));
      efficiencyesdTOFPIDD[icentr]->SetName(Form("efficiencyesdTOFPIDD%d",icentr));

      //fit (mc enhenced sample)
      fittofpid->SetParameters(0.5,0.319,0.0157,0.00664,6.77,2.08);
      efficiencysigTOFPIDD[icentr]->Fit(fittofpid,"R");
      efficiencysigTOFPIDD[icentr]->GetYaxis()->SetTitle("Efficiency");
      fEfficiencyTOFPIDD[icentr] = efficiencysigTOFPIDD[icentr]->GetFunction("fittofpid");

      //fit (esd enhenced sample)
      efficiencysigesdTOFPIDD[icentr]->Fit(fittofpid,"R");
      efficiencysigesdTOFPIDD[icentr]->GetYaxis()->SetTitle("Efficiency");
      fEfficiencyesdTOFPIDD[icentr] = efficiencysigesdTOFPIDD[icentr]->GetFunction("fittofpid");

   }

   // draw (for PbPb, only 1st bin)
   //sig mc
   efficiencysigTOFPIDD[0]->SetTitle("");
   efficiencysigTOFPIDD[0]->SetStats(0);
   efficiencysigTOFPIDD[0]->SetMarkerStyle(25);
   efficiencysigTOFPIDD[0]->SetMarkerColor(2);
   efficiencysigTOFPIDD[0]->SetLineColor(2);
   efficiencysigTOFPIDD[0]->Draw();

   //mb mc
   efficiencyTOFPIDD[0]->SetTitle("");
   efficiencyTOFPIDD[0]->SetStats(0);
   efficiencyTOFPIDD[0]->SetMarkerStyle(24);
   efficiencyTOFPIDD[0]->SetMarkerColor(4);
   efficiencyTOFPIDD[0]->SetLineColor(4);
   efficiencyTOFPIDD[0]->Draw("same");

   //sig esd
   efficiencysigesdTOFPIDD[0]->SetTitle("");
   efficiencysigesdTOFPIDD[0]->SetStats(0);
   efficiencysigesdTOFPIDD[0]->SetMarkerStyle(25);
   efficiencysigesdTOFPIDD[0]->SetMarkerColor(3);
   efficiencysigesdTOFPIDD[0]->SetLineColor(3);
   efficiencysigesdTOFPIDD[0]->Draw("same");

   //mb esd
   efficiencyesdTOFPIDD[0]->SetTitle("");
   efficiencyesdTOFPIDD[0]->SetStats(0);
   efficiencyesdTOFPIDD[0]->SetMarkerStyle(25);
   efficiencyesdTOFPIDD[0]->SetMarkerColor(1);
   efficiencyesdTOFPIDD[0]->SetLineColor(1);
   efficiencyesdTOFPIDD[0]->Draw("same");

   //signal mc fit
   if(fEfficiencyTOFPIDD[0]){
     fEfficiencyTOFPIDD[0]->SetLineColor(2);
     fEfficiencyTOFPIDD[0]->Draw("same");
   }
   //mb esd fit
   if(fEfficiencyesdTOFPIDD[0]){
       fEfficiencyesdTOFPIDD[0]->SetLineColor(3);
       fEfficiencyesdTOFPIDD[0]->Draw("same");
     }

   TLegend *legtofeff = new TLegend(0.3,0.15,0.79,0.44);
   legtofeff->AddEntry(efficiencysigTOFPIDD[0],"TOF PID Step Efficiency","");
   legtofeff->AddEntry(efficiencysigTOFPIDD[0],"vs MC p_{t} for enhenced samples","p");
   legtofeff->AddEntry(efficiencyTOFPIDD[0],"vs MC p_{t} for mb samples","p");
   legtofeff->AddEntry(efficiencysigesdTOFPIDD[0],"vs esd p_{t} for enhenced samples","p");
   legtofeff->AddEntry(efficiencyesdTOFPIDD[0],"vs esd p_{t} for mb samples","p");
   legtofeff->Draw("same");


   TCanvas * cefficiencyParamIP = new TCanvas("efficiencyParamIP","efficiencyParamIP",500,500);
   cefficiencyParamIP->cd();
   gStyle->SetOptStat(0);

   // IP cut efficiencies
   for(int icentr=0; icentr<loopcentr; icentr++)  {

     AliCFContainer *charmCombined = 0x0; 
     AliCFContainer *beautyCombined = 0x0;
     AliCFContainer *beautyCombinedesd = 0x0;

     printf("centrality printing %i \n",icentr);

     source = 0; //charm enhenced
     if(fBeamType==0) mccontainermcD = GetSlicedContainer(container, fNbDimensions, dimensions, source, fChargeChoosen);
     else mccontainermcD = GetSlicedContainer(container, fNbDimensions, dimensions, source, fChargeChoosen, icentr+1);
     AliCFEffGrid* efficiencyCharmSig = new AliCFEffGrid("efficiencyCharmSig","",*mccontainermcD);
     efficiencyCharmSig->CalculateEfficiency(AliHFEcuts::kNcutStepsMCTrack + fStepData,AliHFEcuts::kNcutStepsMCTrack + fStepData-1); // ip cut efficiency. 

     charmCombined= (AliCFContainer*)mccontainermcD->Clone("charmCombined");  

     source = 1; //beauty enhenced
     if(fBeamType==0) mccontainermcD = GetSlicedContainer(container, fNbDimensions, dimensions, source, fChargeChoosen);
     else mccontainermcD = GetSlicedContainer(container, fNbDimensions, dimensions, source, fChargeChoosen, icentr+1);
     AliCFEffGrid* efficiencyBeautySig = new AliCFEffGrid("efficiencyBeautySig","",*mccontainermcD);
     efficiencyBeautySig->CalculateEfficiency(AliHFEcuts::kNcutStepsMCTrack + fStepData,AliHFEcuts::kNcutStepsMCTrack + fStepData-1); // ip cut efficiency. 

     beautyCombined = (AliCFContainer*)mccontainermcD->Clone("beautyCombined"); 

     if(fBeamType==0) mccontaineresdD = GetSlicedContainer(containeresd, fNbDimensions, dimensions, source, fChargeChoosen);
     else mccontaineresdD = GetSlicedContainer(containeresd, fNbDimensions, dimensions, source, fChargeChoosen, icentr+1);
     AliCFEffGrid* efficiencyBeautySigesd = new AliCFEffGrid("efficiencyBeautySigesd","",*mccontaineresdD);
     efficiencyBeautySigesd->CalculateEfficiency(AliHFEcuts::kNcutStepsMCTrack + fStepData,AliHFEcuts::kNcutStepsMCTrack + fStepData-1); // ip cut efficiency.

     beautyCombinedesd = (AliCFContainer*)mccontaineresdD->Clone("beautyCombinedesd");

     source = 0; //charm mb
     if(fBeamType==0) mccontainermcD = GetSlicedContainer(containermb, fNbDimensions, dimensions, source, fChargeChoosen);
     else mccontainermcD = GetSlicedContainer(containermb, fNbDimensions, dimensions, source, fChargeChoosen, icentr+1);
     AliCFEffGrid* efficiencyCharm = new AliCFEffGrid("efficiencyCharm","",*mccontainermcD);
     efficiencyCharm->CalculateEfficiency(AliHFEcuts::kNcutStepsMCTrack + fStepData,AliHFEcuts::kNcutStepsMCTrack + fStepData-1); // ip cut efficiency. 

     charmCombined->Add(mccontainermcD); 
     AliCFEffGrid* efficiencyCharmCombined = new AliCFEffGrid("efficiencyCharmCombined","",*charmCombined); 
     efficiencyCharmCombined->CalculateEfficiency(AliHFEcuts::kNcutStepsMCTrack + fStepData,AliHFEcuts::kNcutStepsMCTrack + fStepData-1); 

     source = 1; //beauty mb
     if(fBeamType==0) mccontainermcD = GetSlicedContainer(containermb, fNbDimensions, dimensions, source, fChargeChoosen);
     else mccontainermcD = GetSlicedContainer(containermb, fNbDimensions, dimensions, source, fChargeChoosen, icentr+1);
     AliCFEffGrid* efficiencyBeauty = new AliCFEffGrid("efficiencyBeauty","",*mccontainermcD);
     efficiencyBeauty->CalculateEfficiency(AliHFEcuts::kNcutStepsMCTrack + fStepData,AliHFEcuts::kNcutStepsMCTrack + fStepData-1); // ip cut efficiency. 

     beautyCombined->Add(mccontainermcD);
     AliCFEffGrid* efficiencyBeautyCombined = new AliCFEffGrid("efficiencyBeautyCombined","",*beautyCombined); 
     efficiencyBeautyCombined->CalculateEfficiency(AliHFEcuts::kNcutStepsMCTrack + fStepData,AliHFEcuts::kNcutStepsMCTrack + fStepData-1); 

     if(fBeamType==0) mccontaineresdD = GetSlicedContainer(containeresdmb, fNbDimensions, dimensions, source, fChargeChoosen);
     else mccontaineresdD = GetSlicedContainer(containeresdmb, fNbDimensions, dimensions, source, fChargeChoosen, icentr+1);
     AliCFEffGrid* efficiencyBeautyesd = new AliCFEffGrid("efficiencyBeautyesd","",*mccontaineresdD);
     efficiencyBeautyesd->CalculateEfficiency(AliHFEcuts::kNcutStepsMCTrack + fStepData,AliHFEcuts::kNcutStepsMCTrack + fStepData-1); // ip cut efficiency.

     beautyCombinedesd->Add(mccontaineresdD);
     AliCFEffGrid* efficiencyBeautyCombinedesd = new AliCFEffGrid("efficiencyBeautyCombinedesd","",*beautyCombinedesd);
     efficiencyBeautyCombinedesd->CalculateEfficiency(AliHFEcuts::kNcutStepsMCTrack + fStepData,AliHFEcuts::kNcutStepsMCTrack + fStepData-1);

     source = 2; //conversion mb
     if(fBeamType==0) mccontainermcD = GetSlicedContainer(containermb, fNbDimensions, dimensions, source, fChargeChoosen);
     else mccontainermcD = GetSlicedContainer(containermb, fNbDimensions, dimensions, source, fChargeChoosen, icentr+1);
     AliCFEffGrid* efficiencyConv = new AliCFEffGrid("efficiencyConv","",*mccontainermcD);
     efficiencyConv->CalculateEfficiency(AliHFEcuts::kNcutStepsMCTrack + fStepData,AliHFEcuts::kNcutStepsMCTrack + fStepData-1); // ip cut efficiency. 

     source = 3; //non HFE except for the conversion mb
     if(fBeamType==0) mccontainermcD = GetSlicedContainer(containermb, fNbDimensions, dimensions, source, fChargeChoosen);
     else mccontainermcD = GetSlicedContainer(containermb, fNbDimensions, dimensions, source, fChargeChoosen, icentr+1);
     AliCFEffGrid* efficiencyNonhfe= new AliCFEffGrid("efficiencyNonhfe","",*mccontainermcD);
     efficiencyNonhfe->CalculateEfficiency(AliHFEcuts::kNcutStepsMCTrack + fStepData,AliHFEcuts::kNcutStepsMCTrack + fStepData-1); // ip cut efficiency.

     if(fIPEffCombinedSamples){
       fEfficiencyCharmSigD[icentr] = (TH1D*)efficiencyCharmCombined->Project(ptpr); //signal enhenced + mb 
       fEfficiencyBeautySigD[icentr] = (TH1D*)efficiencyBeautyCombined->Project(ptpr); //signal enhenced + mb
       fEfficiencyBeautySigesdD[icentr] = (TH1D*)efficiencyBeautyCombinedesd->Project(ptpr); //signal enhenced + mb
     }
     else{
       fEfficiencyCharmSigD[icentr] = (TH1D*)efficiencyCharmSig->Project(ptpr); //signal enhenced only
       fEfficiencyBeautySigD[icentr] = (TH1D*)efficiencyBeautySig->Project(ptpr); //signal enhenced only
       fEfficiencyBeautySigesdD[icentr] = (TH1D*)efficiencyBeautySigesd->Project(ptpr); //signal enhenced only
     }
     fCharmEff[icentr] = (TH1D*)efficiencyCharm->Project(ptpr); //mb only
     fBeautyEff[icentr] = (TH1D*)efficiencyBeauty->Project(ptpr); //mb only
     fConversionEff[icentr] = (TH1D*)efficiencyConv->Project(ptpr); //mb only
     fNonHFEEff[icentr] = (TH1D*)efficiencyNonhfe->Project(ptpr); //mb only

   }

   if(fBeamType==0){
     AliCFEffGrid  *nonHFEEffGrid = (AliCFEffGrid*)  GetEfficiency(GetContainer(kMCWeightedContainerNonHFEESD),1,0);
     fNonHFEEffbgc = (TH1D *) nonHFEEffGrid->Project(0);

     AliCFEffGrid  *conversionEffGrid = (AliCFEffGrid*)  GetEfficiency(GetContainer(kMCWeightedContainerConversionESD),1,0);
     fConversionEffbgc = (TH1D *) conversionEffGrid->Project(0);
   }

   for(int icentr=0; icentr<loopcentr; icentr++)  {
     fipfit->SetParameters(0.5,0.319,0.0157,0.00664,6.77,2.08);
     fipfit->SetLineColor(2);
     fEfficiencyBeautySigD[icentr]->Fit(fipfit,"R");
     fEfficiencyBeautySigD[icentr]->GetYaxis()->SetTitle("Efficiency");
     if(fBeauty2ndMethod)fEfficiencyIPBeautyD[icentr] = fEfficiencyBeautySigD[0]->GetFunction("fipfit"); //why do we need this line?
     else fEfficiencyIPBeautyD[icentr] = fEfficiencyBeautySigD[icentr]->GetFunction("fipfit");

     fipfit->SetParameters(0.5,0.319,0.0157,0.00664,6.77,2.08);
     fipfit->SetLineColor(6);
     fEfficiencyBeautySigesdD[icentr]->Fit(fipfit,"R");
     fEfficiencyBeautySigesdD[icentr]->GetYaxis()->SetTitle("Efficiency");
     if(fBeauty2ndMethod)fEfficiencyIPBeautyesdD[icentr] = fEfficiencyBeautySigesdD[0]->GetFunction("fipfit"); //why do we need this line?
     else fEfficiencyIPBeautyesdD[icentr] = fEfficiencyBeautySigesdD[icentr]->GetFunction("fipfit");

     fipfit->SetParameters(0.5,0.319,0.0157,0.00664,6.77,2.08);
     fipfit->SetLineColor(1);
     fEfficiencyCharmSigD[icentr]->Fit(fipfit,"R");
     fEfficiencyCharmSigD[icentr]->GetYaxis()->SetTitle("Efficiency");
     fEfficiencyIPCharmD[icentr] = fEfficiencyCharmSigD[icentr]->GetFunction("fipfit");
     
     if(fIPParameterizedEff){
       fipfitnonhfe->SetParameters(0.5,0.319,0.0157,0.00664,6.77,2.08);
       fipfitnonhfe->SetLineColor(3);
       fConversionEff[icentr]->Fit(fipfitnonhfe,"R");
       fConversionEff[icentr]->GetYaxis()->SetTitle("Efficiency");
       fEfficiencyIPConversionD[icentr] = fConversionEff[icentr]->GetFunction("fipfitnonhfe");

       fipfitnonhfe->SetParameters(0.5,0.319,0.0157,0.00664,6.77,2.08);
       fipfitnonhfe->SetLineColor(4);
       fNonHFEEff[icentr]->Fit(fipfitnonhfe,"R");
       fNonHFEEff[icentr]->GetYaxis()->SetTitle("Efficiency");
       fEfficiencyIPNonhfeD[icentr] = fNonHFEEff[icentr]->GetFunction("fipfitnonhfe");
     }
   }

   // draw (for PbPb, only 1st bin)
   fEfficiencyCharmSigD[0]->SetMarkerStyle(21);
   fEfficiencyCharmSigD[0]->SetMarkerColor(1);
   fEfficiencyCharmSigD[0]->SetLineColor(1);
   fEfficiencyBeautySigD[0]->SetMarkerStyle(21);
   fEfficiencyBeautySigD[0]->SetMarkerColor(2);
   fEfficiencyBeautySigD[0]->SetLineColor(2);
   fEfficiencyBeautySigesdD[0]->SetStats(0);
   fEfficiencyBeautySigesdD[0]->SetMarkerStyle(21);
   fEfficiencyBeautySigesdD[0]->SetMarkerColor(6);
   fEfficiencyBeautySigesdD[0]->SetLineColor(6);
   fCharmEff[0]->SetMarkerStyle(24);
   fCharmEff[0]->SetMarkerColor(1);
   fCharmEff[0]->SetLineColor(1);
   fBeautyEff[0]->SetMarkerStyle(24);
   fBeautyEff[0]->SetMarkerColor(2);
   fBeautyEff[0]->SetLineColor(2);
   fConversionEff[0]->SetMarkerStyle(24);
   fConversionEff[0]->SetMarkerColor(3);
   fConversionEff[0]->SetLineColor(3);
   fNonHFEEff[0]->SetMarkerStyle(24);
   fNonHFEEff[0]->SetMarkerColor(4);
   fNonHFEEff[0]->SetLineColor(4);

   fEfficiencyCharmSigD[0]->Draw();
   fEfficiencyCharmSigD[0]->GetXaxis()->SetRangeUser(0.0,7.9);
   fEfficiencyCharmSigD[0]->GetYaxis()->SetRangeUser(0.0,0.5);

   fEfficiencyBeautySigD[0]->Draw("same");
   fEfficiencyBeautySigesdD[0]->Draw("same");
   //fCharmEff[0]->Draw("same");
   //fBeautyEff[0]->Draw("same");

   if(fBeamType==0){
     fConversionEffbgc->SetMarkerStyle(25);
     fConversionEffbgc->SetMarkerColor(3);
     fConversionEffbgc->SetLineColor(3);
     fNonHFEEffbgc->SetMarkerStyle(25);
     fNonHFEEffbgc->SetMarkerColor(4);
     fNonHFEEffbgc->SetLineColor(4);
     fConversionEffbgc->Draw("same");
     fNonHFEEffbgc->Draw("same");
   }
   else{
     fConversionEff[0]->Draw("same");
     fNonHFEEff[0]->Draw("same");
   }
   if(fEfficiencyIPBeautyD[0])
      fEfficiencyIPBeautyD[0]->Draw("same");
   if(fEfficiencyIPBeautyesdD[0])
     fEfficiencyIPBeautyesdD[0]->Draw("same");
   if( fEfficiencyIPCharmD[0])
     fEfficiencyIPCharmD[0]->Draw("same");
   if(fIPParameterizedEff){
     fEfficiencyIPConversionD[0]->Draw("same");
     fEfficiencyIPNonhfeD[0]->Draw("same");
   }
   TLegend *legipeff = new TLegend(0.58,0.2,0.88,0.39);
   legipeff->AddEntry(fEfficiencyBeautySigD[0],"IP Step Efficiency","");
   legipeff->AddEntry(fEfficiencyBeautySigD[0],"beauty e","p");
   legipeff->AddEntry(fEfficiencyBeautySigesdD[0],"beauty e(esd pt)","p");
   legipeff->AddEntry(fEfficiencyCharmSigD[0],"charm e","p");
   legipeff->AddEntry(fConversionEffbgc,"conversion e(esd pt)","p");
   legipeff->AddEntry(fNonHFEEffbgc,"Dalitz e(esd pt)","p");
   //legipeff->AddEntry(fConversionEff[0],"conversion e","p");
   //legipeff->AddEntry(fNonHFEEff[0],"Dalitz e","p");
   legipeff->Draw("same");
   gPad->SetGrid();
   //cefficiencyParamIP->SaveAs("efficiencyParamIP.eps");
}

//____________________________________________________________________________
THnSparse* AliHFEspectrum::GetBeautyIPEff(Bool_t isMCpt){
  //
  // Return beauty electron IP cut efficiency
  //

  const Int_t nPtbinning1 = 35;//number of pt bins, according to new binning
  const Int_t nCentralitybinning=11;//number of centrality bins
  Double_t kPtRange[nPtbinning1+1] = { 0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.75, 2., 2.25, 2.5, 2.75, 3., 3.5, 4., 4.5, 5., 5.5, 6., 7., 8., 10., 12., 14., 16., 18., 20.};//pt bin limits
  Double_t kCentralityRange[nCentralitybinning+1] = {0.,1.,2., 3., 4., 5., 6., 7.,8.,9., 10., 11.};
  //Int_t ptpr = 0;
  Int_t nDim=1;  //dimensions of the efficiency weighting grid
  if(fBeamType==1)
  {
    //ptpr=1;
    nDim=2; //dimensions of the efficiency weighting grid
  }
  Int_t nBin[1] = {nPtbinning1};
  Int_t nBinPbPb[2] = {nCentralitybinning,nPtbinning1};


  THnSparseF *ipcut;
  if(fBeamType==0) ipcut = new THnSparseF("beff", "b IP efficiency; p_{t}(GeV/c)", nDim, nBin);
  else ipcut = new THnSparseF("beff", "b IP efficiency; centrality bin; p_{t}(GeV/c)", nDim, nBinPbPb);
 
  for(Int_t idim = 0; idim < nDim; idim++)
  {
    if(nDim==1) ipcut->SetBinEdges(idim, kPtRange);
    if(nDim==2)
      {
        ipcut->SetBinEdges(0, kCentralityRange);
        ipcut->SetBinEdges(1, kPtRange);
      }
  }
  Double_t pt[1];
  Double_t weight[nCentralitybinning];
  Double_t weightErr[nCentralitybinning];
  Double_t contents[2];

  for(Int_t a=0;a<11;a++)
  {
      weight[a] = 1.0;
      weightErr[a] = 1.0;
  }


  Int_t looppt=nBin[0];
  Int_t loopcentr=1;
  Int_t ibin[2];
  if(fBeamType==1)
  {
      loopcentr=nBinPbPb[0];
  }


  for(int icentr=0; icentr<loopcentr; icentr++)
  {
      for(int i=0; i<looppt; i++)
      {
	  pt[0]=(kPtRange[i]+kPtRange[i+1])/2.;
          if(isMCpt){
            if(fEfficiencyIPBeautyD[icentr]){
              weight[icentr]=fEfficiencyIPBeautyD[icentr]->Eval(pt[0]);
              weightErr[icentr] = 0;
            }
            else{
              printf("Fit failed on beauty IP cut efficiency for centrality %d. Contents in histo used!\n",icentr);
              weight[icentr] = fEfficiencyBeautySigD[icentr]->GetBinContent(i+1); 
              weightErr[icentr] = fEfficiencyBeautySigD[icentr]->GetBinError(i+1);
            }
          }
          else{
            if(fEfficiencyIPBeautyesdD[icentr]){
              weight[icentr]=fEfficiencyIPBeautyesdD[icentr]->Eval(pt[0]);
              weightErr[icentr] = 0;
            }
            else{
              printf("Fit failed on beauty IP cut efficiency for centrality %d. Contents in histo used!\n",icentr);
              weight[icentr] = fEfficiencyBeautySigesdD[icentr]->GetBinContent(i+1);
              weightErr[icentr] = fEfficiencyBeautySigD[icentr]->GetBinError(i+1);
            }
          }

          if(fBeamType==1){
              contents[0]=icentr;
              contents[1]=pt[0];
              ibin[0]=icentr;
              ibin[1]=i+1;
          }
          if(fBeamType==0){
              contents[0]=pt[0];
              ibin[0]=i+1;
          }
          ipcut->Fill(contents,weight[icentr]);
          ipcut->SetBinError(ibin,weightErr[icentr]);
      }
  } 

  Int_t nDimSparse = ipcut->GetNdimensions();
  Int_t* binsvar = new Int_t[nDimSparse]; // number of bins for each variable
  Long_t nBins = 1; // used to calculate the total number of bins in the THnSparse

  for (Int_t iVar=0; iVar<nDimSparse; iVar++) {
      binsvar[iVar] = ipcut->GetAxis(iVar)->GetNbins();
      nBins *= binsvar[iVar];
  }

  Int_t *binfill = new Int_t[nDimSparse]; // bin to fill the THnSparse (holding the bin coordinates)
  // loop that sets 0 error in each bin
  for (Long_t iBin=0; iBin<nBins; iBin++) {
    Long_t bintmp = iBin ;
    for (Int_t iVar=0; iVar<nDimSparse; iVar++) {
      binfill[iVar] = 1 + bintmp % binsvar[iVar] ;
      bintmp /= binsvar[iVar] ;
    }
    //ipcut->SetBinError(binfill,0.); // put 0 everywhere
  }

  delete[] binsvar;
  delete[] binfill;

  return ipcut;
}

//____________________________________________________________________________
THnSparse* AliHFEspectrum::GetPIDxIPEff(Int_t source){
  //
  // Return PID x IP cut efficiency
  //
    const Int_t nPtbinning1 = 35;//number of pt bins, according to new binning
    const Int_t nCentralitybinning=11;//number of centrality bins
    Double_t kPtRange[nPtbinning1+1] = { 0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.75, 2., 2.25, 2.5, 2.75, 3., 3.5, 4., 4.5, 5., 5.5, 6., 7., 8., 10., 12., 14., 16., 18., 20.};//pt bin limits
    Double_t kCentralityRange[nCentralitybinning+1] = {0.,1.,2., 3., 4., 5., 6., 7.,8.,9., 10., 11.};
    //Int_t ptpr = 0;
    Int_t nDim=1;  //dimensions of the efficiency weighting grid
    if(fBeamType==1)
    {
      //ptpr=1;
	nDim=2; //dimensions of the efficiency weighting grid
    }
    Int_t nBin[1] = {nPtbinning1};
    Int_t nBinPbPb[2] = {nCentralitybinning,nPtbinning1};

    THnSparseF *pideff;
    if(fBeamType==0) pideff = new THnSparseF("pideff", "PID efficiency; p_{t}(GeV/c)", nDim, nBin);
    else pideff = new THnSparseF("pideff", "PID efficiency; centrality bin; p_{t}(GeV/c)", nDim, nBinPbPb);
    for(Int_t idim = 0; idim < nDim; idim++)
    {

	if(nDim==1) pideff->SetBinEdges(idim, kPtRange);
	if(nDim==2)
	{
	    pideff->SetBinEdges(0, kCentralityRange);
	    pideff->SetBinEdges(1, kPtRange);
	}
    }

  Double_t pt[1];
  Double_t weight[nCentralitybinning];
  Double_t weightErr[nCentralitybinning];
  Double_t contents[2];

  for(Int_t a=0;a<nCentralitybinning;a++)
  {
      weight[a] = 1.0;
      weightErr[a] = 1.0;
  }

  Int_t looppt=nBin[0];
  Int_t loopcentr=1;
  Int_t ibin[2];
  if(fBeamType==1)
  {
      loopcentr=nBinPbPb[0];
  }

  for(int icentr=0; icentr<loopcentr; icentr++)
  {
      Double_t trdtpcPidEfficiency = fEfficiencyFunction->Eval(0); // assume we have constant TRD+TPC PID efficiency
      for(int i=0; i<looppt; i++)
      {
	  pt[0]=(kPtRange[i]+kPtRange[i+1])/2.;

          Double_t tofpideff = 0.;
          Double_t tofpideffesd = 0.;
          if(fEfficiencyTOFPIDD[icentr])
            tofpideff = fEfficiencyTOFPIDD[icentr]->Eval(pt[0]); 
          else{
            printf("TOF PID fit failed on conversion for centrality %d. The result is wrong!\n",icentr);
          }  
          if(fEfficiencyesdTOFPIDD[icentr])
            tofpideffesd = fEfficiencyesdTOFPIDD[icentr]->Eval(pt[0]);
          else{
            printf("TOF PID fit failed on conversion for centrality %d. The result is wrong!\n",icentr);
          }

          //tof pid eff x tpc pid eff x ip cut eff
          if(fIPParameterizedEff){
            if(source==0) {
              if(fEfficiencyIPCharmD[icentr]){
                weight[icentr] = tofpideff*trdtpcPidEfficiency*fEfficiencyIPCharmD[icentr]->Eval(pt[0]);
                weightErr[icentr] = 0; 
              }
              else{
                printf("Fit failed on charm IP cut efficiency for centrality %d\n",icentr);
                weight[icentr] = tofpideff*trdtpcPidEfficiency*fEfficiencyCharmSigD[icentr]->GetBinContent(i+1);
                weightErr[icentr] = tofpideff*trdtpcPidEfficiency*fEfficiencyCharmSigD[icentr]->GetBinError(i+1); 
              }
            } 
	    else if(source==2) {
              if(fEfficiencyIPConversionD[icentr]){
                weight[icentr] = tofpideffesd*trdtpcPidEfficiency*fEfficiencyIPConversionD[icentr]->Eval(pt[0]); 
                weightErr[icentr] = 0; 
              }
              else{
                printf("Fit failed on conversion IP cut efficiency for centrality %d\n",icentr);
                weight[icentr] = tofpideffesd*trdtpcPidEfficiency*fConversionEff[icentr]->GetBinContent(i+1);
                weightErr[icentr] = tofpideffesd*trdtpcPidEfficiency*fConversionEff[icentr]->GetBinError(i+1);
              }
            }
	    else if(source==3) {
              if(fEfficiencyIPNonhfeD[icentr]){
                weight[icentr] = tofpideffesd*trdtpcPidEfficiency*fEfficiencyIPNonhfeD[icentr]->Eval(pt[0]); 
                weightErr[icentr] = 0; 
              }
              else{
                printf("Fit failed on dalitz IP cut efficiency for centrality %d\n",icentr);
                weight[icentr] = tofpideffesd*trdtpcPidEfficiency*fNonHFEEff[icentr]->GetBinContent(i+1);
                weightErr[icentr] = tofpideffesd*trdtpcPidEfficiency*fNonHFEEff[icentr]->GetBinError(i+1);
              }  
            }
          }
          else{
            if(source==0){ 
              if(fEfficiencyIPCharmD[icentr]){
                weight[icentr] = tofpideff*trdtpcPidEfficiency*fEfficiencyIPCharmD[icentr]->Eval(pt[0]);
                weightErr[icentr] = 0;
              }
              else{
                printf("Fit failed on charm IP cut efficiency for centrality %d\n",icentr);
                weight[icentr] = tofpideff*trdtpcPidEfficiency*fEfficiencyCharmSigD[icentr]->GetBinContent(i+1);
                weightErr[icentr] = tofpideff*trdtpcPidEfficiency*fEfficiencyCharmSigD[icentr]->GetBinError(i+1);
              }
            }
	    else if(source==2){
              if(fBeamType==0){
                weight[icentr] = tofpideffesd*trdtpcPidEfficiency*fConversionEffbgc->GetBinContent(i+1); // conversion
                weightErr[icentr] = tofpideffesd*trdtpcPidEfficiency*fConversionEffbgc->GetBinError(i+1);
              }
              else{
                weight[icentr] = tofpideffesd*trdtpcPidEfficiency*fConversionEff[icentr]->GetBinContent(i+1); // conversion
                weightErr[icentr] = tofpideffesd*trdtpcPidEfficiency*fConversionEff[icentr]->GetBinError(i+1);
              }
            }
	    else if(source==3){
              if(fBeamType==0){
                weight[icentr] = tofpideffesd*trdtpcPidEfficiency*fNonHFEEffbgc->GetBinContent(i+1); // conversion
                weightErr[icentr] = tofpideffesd*trdtpcPidEfficiency*fNonHFEEffbgc->GetBinError(i+1);
              }
              else{ 
                weight[icentr] = tofpideffesd*trdtpcPidEfficiency*fNonHFEEff[icentr]->GetBinContent(i+1); // Dalitz
                weightErr[icentr] = tofpideffesd*trdtpcPidEfficiency*fNonHFEEff[icentr]->GetBinError(i+1);
              }
            }
          }

	  if(fBeamType==1){
	      contents[0]=icentr;
	      contents[1]=pt[0];
              ibin[0]=icentr;
              ibin[1]=i+1;
	  }
	  if(fBeamType==0){
	      contents[0]=pt[0];
              ibin[0]=i+1;
	  }

	  pideff->Fill(contents,weight[icentr]);
          pideff->SetBinError(ibin,weightErr[icentr]);
      }
  }

  Int_t nDimSparse = pideff->GetNdimensions();
  Int_t* binsvar = new Int_t[nDimSparse]; // number of bins for each variable
  Long_t nBins = 1; // used to calculate the total number of bins in the THnSparse

  for (Int_t iVar=0; iVar<nDimSparse; iVar++) {
      binsvar[iVar] = pideff->GetAxis(iVar)->GetNbins();
      nBins *= binsvar[iVar];
  }

  Int_t *binfill = new Int_t[nDimSparse]; // bin to fill the THnSparse (holding the bin coordinates)
  // loop that sets 0 error in each bin
  for (Long_t iBin=0; iBin<nBins; iBin++) {
    Long_t bintmp = iBin ;
    for (Int_t iVar=0; iVar<nDimSparse; iVar++) {
      binfill[iVar] = 1 + bintmp % binsvar[iVar] ;
      bintmp /= binsvar[iVar] ;
    }
  }

  delete[] binsvar;
  delete[] binfill;


  return pideff;
}

//__________________________________________________________________________
AliCFDataGrid *AliHFEspectrum::GetRawBspectra2ndMethod(){
 //
 // retrieve AliCFDataGrid for raw beauty spectra obtained from fit method
    //
    Int_t ptpr = 0;
    Int_t nDim = 1;
    if(fBeamType==0)
    {
	ptpr=0;
    }
    if(fBeamType==1)
    {
	ptpr=1;
	nDim=2;
    }

    const Int_t nPtbinning1 = 18;//number of pt bins, according to new binning
    const Int_t nCentralitybinning=11;//number of centrality bins
    Double_t kPtRange[19] = {0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 12, 16, 20};
   
    Double_t kCentralityRange[nCentralitybinning+1] = {0.,1.,2., 3., 4., 5., 6., 7.,8.,9., 10., 11.};
    Int_t nBin[1] = {nPtbinning1};
    Int_t nBinPbPb[2] = {nCentralitybinning,nPtbinning1};

    AliCFDataGrid *rawBeautyContainer;
    if(fBeamType==0)  rawBeautyContainer = new AliCFDataGrid("rawBeautyContainer","rawBeautyContainer",nDim,nBin);
    else rawBeautyContainer = new AliCFDataGrid("rawBeautyContainer","rawBeautyContainer",nDim,nBinPbPb);
    //  printf("number of bins= %d\n",bins[0]);


    
    
    THnSparseF *brawspectra;
    if(fBeamType==0) brawspectra= new THnSparseF("brawspectra", "beauty yields ; p_{t}(GeV/c)", nDim, nBin);
    else brawspectra= new THnSparseF("brawspectra", "beauty yields ; p_{t}(GeV/c)", nDim, nBinPbPb);
    if(fBeamType==0) brawspectra->SetBinEdges(0, kPtRange);
    if(fBeamType==1)
      {
	//      brawspectra->SetBinEdges(0, centralityBins);
	brawspectra->SetBinEdges(0, kCentralityRange);
	brawspectra->SetBinEdges(1, kPtRange);
      }
    
    Double_t pt[1];
    Double_t yields= 0.;
    Double_t valuesb[2];
    
    //Int_t looppt=nBin[0];
    Int_t loopcentr=1;
    if(fBeamType==1)
      {
	loopcentr=nBinPbPb[0];
      }
    
    for(int icentr=0; icentr<loopcentr; icentr++)
      {
	
	for(int i=0; i<fBSpectrum2ndMethod->GetNbinsX(); i++){
	  pt[0]=(kPtRange[i]+kPtRange[i+1])/2.;
	  
	  yields = fBSpectrum2ndMethod->GetBinContent(i+1);
	  
	  if(fBeamType==1)
	    {
	      valuesb[0]=icentr;
	      valuesb[1]=pt[0];
	    }
	  if(fBeamType==0) valuesb[0]=pt[0];
	  brawspectra->Fill(valuesb,yields);
	}
      }
    
    
    
    Int_t nDimSparse = brawspectra->GetNdimensions();
    Int_t* binsvar = new Int_t[nDimSparse]; // number of bins for each variable
    Long_t nBins = 1; // used to calculate the total number of bins in the THnSparse
    
    for (Int_t iVar=0; iVar<nDimSparse; iVar++) {
      binsvar[iVar] = brawspectra->GetAxis(iVar)->GetNbins();
      nBins *= binsvar[iVar];
    }
    
    Int_t *binfill = new Int_t[nDimSparse]; // bin to fill the THnSparse (holding the bin coordinates)
    // loop that sets 0 error in each bin
    for (Long_t iBin=0; iBin<nBins; iBin++) {
      Long_t bintmp = iBin ;
      for (Int_t iVar=0; iVar<nDimSparse; iVar++) {
	binfill[iVar] = 1 + bintmp % binsvar[iVar] ;
	bintmp /= binsvar[iVar] ;
      }
      brawspectra->SetBinError(binfill,0.); // put 0 everywhere
    }
    
    
    rawBeautyContainer->SetGrid(brawspectra); // get charm efficiency
    TH1D* hRawBeautySpectra = (TH1D*)rawBeautyContainer->Project(ptpr);
    
    new TCanvas;
    fBSpectrum2ndMethod->SetMarkerStyle(24);
    fBSpectrum2ndMethod->Draw("p");
    hRawBeautySpectra->SetMarkerStyle(25);
    hRawBeautySpectra->Draw("samep");

    delete[] binfill;
    delete[] binsvar; 

    return rawBeautyContainer;
}

//__________________________________________________________________________
void AliHFEspectrum::CalculateNonHFEsyst(Int_t centrality){
  //
  // Calculate non HFE sys
  //
  //

  if(!fNonHFEsyst)
    return;

  Double_t evtnorm[1] = {0.0};
  if(fNMCbgEvents[0]>0) evtnorm[0]= double(fNEvents[0])/double(fNMCbgEvents[0]);
  
  AliCFDataGrid *convSourceGrid[kElecBgSources][kBgLevels];
  AliCFDataGrid *nonHFESourceGrid[kElecBgSources][kBgLevels];

  AliCFDataGrid *bgLevelGrid[2][kBgLevels];//for pi0 and eta based errors
  AliCFDataGrid *bgNonHFEGrid[kBgLevels];
  AliCFDataGrid *bgConvGrid[kBgLevels];

  Int_t stepbackground = 3;
  Int_t* bins=new Int_t[1];
  const Char_t *bgBase[2] = {"pi0","eta"};
 
  bins[0]=fConversionEff[centrality]->GetNbinsX();
   
  AliCFDataGrid *weightedConversionContainer = new AliCFDataGrid("weightedConversionContainer","weightedConversionContainer",1,bins);
  AliCFDataGrid *weightedNonHFEContainer = new AliCFDataGrid("weightedNonHFEContainer","weightedNonHFEContainer",1,bins);

  for(Int_t iLevel = 0; iLevel < kBgLevels; iLevel++){   
    for(Int_t iSource = 0; iSource < kElecBgSources; iSource++){
      convSourceGrid[iSource][iLevel] = new AliCFDataGrid(Form("convGrid_%d_%d_%d",iSource,iLevel,centrality),Form("convGrid_%d_%d_%d",iSource,iLevel,centrality),*fConvSourceContainer[iSource][iLevel][centrality],stepbackground);
      weightedConversionContainer->SetGrid(GetPIDxIPEff(2));
      convSourceGrid[iSource][iLevel]->Multiply(weightedConversionContainer,1.0);
      
      nonHFESourceGrid[iSource][iLevel] = new AliCFDataGrid(Form("nonHFEGrid_%d_%d_%d",iSource,iLevel,centrality),Form("nonHFEGrid_%d_%d_%d",iSource,iLevel,centrality),*fNonHFESourceContainer[iSource][iLevel][centrality],stepbackground);
      weightedNonHFEContainer->SetGrid(GetPIDxIPEff(3));
      nonHFESourceGrid[iSource][iLevel]->Multiply(weightedNonHFEContainer,1.0);
    }
    
    bgConvGrid[iLevel] = (AliCFDataGrid*)convSourceGrid[0][iLevel]->Clone();
    for(Int_t iSource = 2; iSource < kElecBgSources; iSource++){
      bgConvGrid[iLevel]->Add(convSourceGrid[iSource][iLevel]);
    }
    if(!fEtaSyst)
      bgConvGrid[iLevel]->Add(convSourceGrid[1][iLevel]);
    
    bgNonHFEGrid[iLevel] = (AliCFDataGrid*)nonHFESourceGrid[0][iLevel]->Clone(); 
    for(Int_t iSource = 2; iSource < kElecBgSources; iSource++){//add other sources to pi0, to get overall background from all meson decays, exception: eta (independent error calculation)
      bgNonHFEGrid[iLevel]->Add(nonHFESourceGrid[iSource][iLevel]);
    }
    if(!fEtaSyst)
      bgNonHFEGrid[iLevel]->Add(nonHFESourceGrid[1][iLevel]);
    
    bgLevelGrid[0][iLevel] = (AliCFDataGrid*)bgConvGrid[iLevel]->Clone();
    bgLevelGrid[0][iLevel]->Add(bgNonHFEGrid[iLevel]);
    if(fEtaSyst){
      bgLevelGrid[1][iLevel] = (AliCFDataGrid*)nonHFESourceGrid[1][iLevel]->Clone();//background for eta source
      bgLevelGrid[1][iLevel]->Add(convSourceGrid[1][iLevel]);
    }
  }
 
  
  //Now subtract the mean from upper, and lower from mean container to get the error based on the pion yield uncertainty (-> this error sums linearly, since its contribution to all meson yields is correlated; exception: eta errors in pp 7 TeV sum with others the gaussian way, as they are independent from pi0) 
  AliCFDataGrid *bgErrorGrid[2][2];//for pions/eta error base, for lower/upper
  TH1D* hBaseErrors[2][2];//pi0/eta and lower/upper
  for(Int_t iErr = 0; iErr < 2; iErr++){//errors for pi0 and eta base
    bgErrorGrid[iErr][0] = (AliCFDataGrid*)bgLevelGrid[iErr][1]->Clone();
    bgErrorGrid[iErr][0]->Add(bgLevelGrid[iErr][0],-1.);
    bgErrorGrid[iErr][1] = (AliCFDataGrid*)bgLevelGrid[iErr][2]->Clone();    
    bgErrorGrid[iErr][1]->Add(bgLevelGrid[iErr][0],-1.);

  //plot absolute differences between limit yields (upper/lower limit, based on pi0 and eta errors) and best estimate
 
    hBaseErrors[iErr][0] = (TH1D*)bgErrorGrid[iErr][0]->Project(0);
    hBaseErrors[iErr][0]->Scale(-1.);
    hBaseErrors[iErr][0]->SetTitle(Form("Absolute %s-based systematic errors from non-HF meson decays and conversions",bgBase[iErr]));
    hBaseErrors[iErr][1] = (TH1D*)bgErrorGrid[iErr][1]->Project(0);
    if(!fEtaSyst)break;
  }
  
 
  //Calculate the scaling errors for electrons from all mesons except for pions (and in pp 7 TeV case eta): square sum of (0.3 * best yield estimate), where 0.3 is the error generally assumed for m_t scaling
  TH1D *hSpeciesErrors[kElecBgSources-1];
  for(Int_t iSource = 1; iSource < kElecBgSources; iSource++){
    if(fEtaSyst && (iSource == 1))continue;
    hSpeciesErrors[iSource-1] = (TH1D*)convSourceGrid[iSource][0]->Project(0);
    TH1D *hNonHFEtemp = (TH1D*)nonHFESourceGrid[iSource][0]->Project(0);
    hSpeciesErrors[iSource-1]->Add(hNonHFEtemp);
    hSpeciesErrors[iSource-1]->Scale(0.3);   
  }
  
  //Int_t firstBgSource = 0;//if eta systematics are not from scaling
  //if(fEtaSyst){firstBgSource = 1;}//source 0 histograms are not filled if eta errors are independently determined!
  TH1D *hOverallSystErrLow = (TH1D*)hSpeciesErrors[1]->Clone();
  TH1D *hOverallSystErrUp = (TH1D*)hSpeciesErrors[1]->Clone();
  TH1D *hScalingErrors = (TH1D*)hSpeciesErrors[1]->Clone();

  TH1D *hOverallBinScaledErrsUp = (TH1D*)hOverallSystErrUp->Clone();
  TH1D *hOverallBinScaledErrsLow = (TH1D*)hOverallSystErrLow->Clone();

  for(Int_t iBin = 1; iBin <= kBgPtBins; iBin++){
    Double_t pi0basedErrLow,pi0basedErrUp,etaErrLow,etaErrUp;    
    pi0basedErrLow = hBaseErrors[0][0]->GetBinContent(iBin); 
    pi0basedErrUp = hBaseErrors[0][1]->GetBinContent(iBin);
    if(fEtaSyst){
      etaErrLow = hBaseErrors[1][0]->GetBinContent(iBin); 
      etaErrUp = hBaseErrors[1][1]->GetBinContent(iBin);
    }
    else{ etaErrLow = etaErrUp = 0.;}

    Double_t sqrsumErrs= 0;
    for(Int_t iSource = 1; iSource < kElecBgSources; iSource++){
      if(fEtaSyst && (iSource == 1))continue;
      Double_t scalingErr=hSpeciesErrors[iSource-1]->GetBinContent(iBin);
      sqrsumErrs+=(scalingErr*scalingErr);
    }
    for(Int_t iErr = 0; iErr < 2; iErr++){
      for(Int_t iLevel = 0; iLevel < 2; iLevel++){
        hBaseErrors[iErr][iLevel]->SetBinContent(iBin,hBaseErrors[iErr][iLevel]->GetBinContent(iBin)/hBaseErrors[iErr][iLevel]->GetBinWidth(iBin));
      }
      if(!fEtaSyst)break;
    }
    hOverallSystErrUp->SetBinContent(iBin, TMath::Sqrt((pi0basedErrUp*pi0basedErrUp)+(etaErrUp*etaErrUp)+sqrsumErrs));
    hOverallSystErrLow->SetBinContent(iBin, TMath::Sqrt((pi0basedErrLow*pi0basedErrLow)+(etaErrLow*etaErrLow)+sqrsumErrs));
    hScalingErrors->SetBinContent(iBin, TMath::Sqrt(sqrsumErrs)/hScalingErrors->GetBinWidth(iBin));

    hOverallBinScaledErrsUp->SetBinContent(iBin,hOverallSystErrUp->GetBinContent(iBin)/hOverallBinScaledErrsUp->GetBinWidth(iBin));
    hOverallBinScaledErrsLow->SetBinContent(iBin,hOverallSystErrLow->GetBinContent(iBin)/hOverallBinScaledErrsLow->GetBinWidth(iBin));		  
  }
   
  
  TCanvas *cPiErrors = new TCanvas("cPiErrors","cPiErrors",1000,600);
  cPiErrors->cd();
  cPiErrors->SetLogx();
  cPiErrors->SetLogy();
  hBaseErrors[0][0]->Draw();
  //hBaseErrors[0][1]->SetMarkerColor(kBlack);
  //hBaseErrors[0][1]->SetLineColor(kBlack);
  //hBaseErrors[0][1]->Draw("SAME");
  if(fEtaSyst){
    hBaseErrors[1][0]->Draw("SAME");
    hBaseErrors[1][0]->SetMarkerColor(kBlack);
    hBaseErrors[1][0]->SetLineColor(kBlack);
  //hBaseErrors[1][1]->SetMarkerColor(13);
  //hBaseErrors[1][1]->SetLineColor(13);
  //hBaseErrors[1][1]->Draw("SAME");
  }
  //hOverallBinScaledErrsUp->SetMarkerColor(kBlue);
  //hOverallBinScaledErrsUp->SetLineColor(kBlue);
  //hOverallBinScaledErrsUp->Draw("SAME");
  hOverallBinScaledErrsLow->SetMarkerColor(kGreen);
  hOverallBinScaledErrsLow->SetLineColor(kGreen);
  hOverallBinScaledErrsLow->Draw("SAME");
  hScalingErrors->SetLineColor(kBlue);
  hScalingErrors->Draw("SAME");

  TLegend *lPiErr = new TLegend(0.6,0.6, 0.95,0.95);
  lPiErr->AddEntry(hBaseErrors[0][0],"Lower error from pion error");
  //lPiErr->AddEntry(hBaseErrors[0][1],"Upper error from pion error");
  if(fEtaSyst){
  lPiErr->AddEntry(hBaseErrors[1][0],"Lower error from eta error");
  //lPiErr->AddEntry(hBaseErrors[1][1],"Upper error from eta error");
  }
  lPiErr->AddEntry(hScalingErrors, "scaling error");
  lPiErr->AddEntry(hOverallBinScaledErrsLow, "overall lower systematics");
  //lPiErr->AddEntry(hOverallBinScaledErrsUp, "overall upper systematics");
  lPiErr->Draw("SAME");

  //Normalize errors
  TH1D *hUpSystScaled = (TH1D*)hOverallSystErrUp->Clone();
  TH1D *hLowSystScaled = (TH1D*)hOverallSystErrLow->Clone();
  hUpSystScaled->Scale(evtnorm[0]);//scale by N(data)/N(MC), to make data sets comparable to saved subtracted spectrum (calculations in separate macro!)
  hLowSystScaled->Scale(evtnorm[0]);
  TH1D *hNormAllSystErrUp = (TH1D*)hUpSystScaled->Clone();
  TH1D *hNormAllSystErrLow = (TH1D*)hLowSystScaled->Clone();
  //histograms to be normalized to TGraphErrors
  CorrectFromTheWidth(hNormAllSystErrUp);
  CorrectFromTheWidth(hNormAllSystErrLow);

  TCanvas *cNormOvErrs = new TCanvas("cNormOvErrs","cNormOvErrs");
  cNormOvErrs->cd();
  cNormOvErrs->SetLogx();
  cNormOvErrs->SetLogy();

  TGraphErrors* gOverallSystErrUp = NormalizeTH1(hNormAllSystErrUp);
  TGraphErrors* gOverallSystErrLow = NormalizeTH1(hNormAllSystErrLow);
  gOverallSystErrUp->SetTitle("Overall Systematic non-HFE Errors");
  gOverallSystErrUp->SetMarkerColor(kBlack);
  gOverallSystErrUp->SetLineColor(kBlack);
  gOverallSystErrLow->SetMarkerColor(kRed);
  gOverallSystErrLow->SetLineColor(kRed);
  gOverallSystErrUp->Draw("AP");
  gOverallSystErrLow->Draw("PSAME");
  TLegend *lAllSys = new TLegend(0.4,0.6,0.89,0.89);
  lAllSys->AddEntry(gOverallSystErrLow,"lower","p");
  lAllSys->AddEntry(gOverallSystErrUp,"upper","p");
  lAllSys->Draw("same");


  AliCFDataGrid *bgYieldGrid;
  if(fEtaSyst){
    bgLevelGrid[0][0]->Add(bgLevelGrid[1][0]);//Addition of the eta background best estimate to the rest. Needed to be separated for error treatment - now overall background necessary! If no separate eta systematics exist, the corresponding grid has already been added before.
  }
  bgYieldGrid = (AliCFDataGrid*)bgLevelGrid[0][0]->Clone();

  TH1D *hBgYield = (TH1D*)bgYieldGrid->Project(0);
  TH1D* hRelErrUp = (TH1D*)hOverallSystErrUp->Clone();
  hRelErrUp->Divide(hBgYield);
  TH1D* hRelErrLow = (TH1D*)hOverallSystErrLow->Clone();
  hRelErrLow->Divide(hBgYield);

  TCanvas *cRelErrs = new TCanvas("cRelErrs","cRelErrs");
  cRelErrs->cd();
  cRelErrs->SetLogx();
  hRelErrUp->SetTitle("Relative error of non-HFE background yield");
  hRelErrUp->Draw();
  hRelErrLow->SetLineColor(kBlack);
  hRelErrLow->Draw("SAME");

  TLegend *lRel = new TLegend(0.6,0.6,0.95,0.95);
  lRel->AddEntry(hRelErrUp, "upper");
  lRel->AddEntry(hRelErrLow, "lower");
  lRel->Draw("SAME");

  //CorrectFromTheWidth(hBgYield);
  //hBgYield->Scale(evtnorm[0]);
 
 
  //write histograms/TGraphs to file
  TFile *output = new TFile("systHists.root","recreate");

  hBgYield->SetName("hBgYield");  
  hBgYield->Write();
  hRelErrUp->SetName("hRelErrUp");
  hRelErrUp->Write();
  hRelErrLow->SetName("hRelErrLow");
  hRelErrLow->Write();
  hUpSystScaled->SetName("hOverallSystErrUp");
  hUpSystScaled->Write();
  hLowSystScaled->SetName("hOverallSystErrLow");
  hLowSystScaled->Write();
  gOverallSystErrUp->SetName("gOverallSystErrUp");
  gOverallSystErrUp->Write();
  gOverallSystErrLow->SetName("gOverallSystErrLow");
  gOverallSystErrLow->Write(); 

  output->Close(); 
  delete output;
  
}

//____________________________________________________________
void AliHFEspectrum::UnfoldBG(AliCFDataGrid* const bgsubpectrum){

  //
  // Unfold backgrounds to check its sanity
  //

  AliCFContainer *mcContainer = GetContainer(kMCContainerCharmMC);
  //AliCFContainer *mcContainer = GetContainer(kMCContainerMC);
  if(!mcContainer){
    AliError("MC Container not available");
    return;
  }

  if(!fCorrelation){
    AliError("No Correlation map available");
    return;
  }

  // Data 
  AliCFDataGrid *dataGrid = 0x0;
  dataGrid = bgsubpectrum;

  // Guessed
  AliCFDataGrid* guessedGrid = new AliCFDataGrid("guessed","",*mcContainer, fStepGuessedUnfolding);
  THnSparse* guessedTHnSparse = ((AliCFGridSparse*)guessedGrid->GetData())->GetGrid();

  // Efficiency
  AliCFEffGrid* efficiencyD = new AliCFEffGrid("efficiency","",*mcContainer);
  efficiencyD->CalculateEfficiency(fStepMC+2,fStepTrue);

  // Unfold background spectra
  Int_t nDim=1;
  if(fBeamType==0)nDim = 1;
  if(fBeamType==1)nDim = 2;
  AliCFUnfolding unfolding("unfolding","",nDim,fCorrelation,efficiencyD->GetGrid(),dataGrid->GetGrid(),guessedTHnSparse);
  if(fUnSetCorrelatedErrors) unfolding.UnsetCorrelatedErrors();
  unfolding.SetMaxNumberOfIterations(fNumberOfIterations);
  if(fSetSmoothing) unfolding.UseSmoothing();
  unfolding.Unfold();

  // Results
  THnSparse* result = unfolding.GetUnfolded();
  TCanvas *ctest = new TCanvas("yvonnetest","yvonnetest",1000,600);
  if(fBeamType==1)
  {
      ctest->Divide(2,2);
      ctest->cd(1);
      result->GetAxis(0)->SetRange(1,1);
      TH1D* htest1=(TH1D*)result->Projection(0);
      htest1->Draw();
      ctest->cd(2);
      result->GetAxis(0)->SetRange(1,1);
      TH1D* htest2=(TH1D*)result->Projection(1);
      htest2->Draw();
      ctest->cd(3);
      result->GetAxis(0)->SetRange(6,6);
      TH1D* htest3=(TH1D*)result->Projection(0);
      htest3->Draw();
      ctest->cd(4);
      result->GetAxis(0)->SetRange(6,6);
      TH1D* htest4=(TH1D*)result->Projection(1);
      htest4->Draw();

  }





  TGraphErrors* unfoldedbgspectrumD = Normalize(result);
  if(!unfoldedbgspectrumD) {
    AliError("Unfolded background spectrum doesn't exist");
  }
  else{
    TFile *file = TFile::Open("unfoldedbgspectrum.root","recreate");
    if(fBeamType==0)unfoldedbgspectrumD->Write("unfoldedbgspectrum");

    if(fBeamType==1)
    {
        Int_t centr=1;
	result->GetAxis(0)->SetRange(centr,centr);
	unfoldedbgspectrumD = Normalize(result,centr-1);
	unfoldedbgspectrumD->Write("unfoldedbgspectrum_centr0_20");
        centr=6;
	result->GetAxis(0)->SetRange(centr,centr);
	unfoldedbgspectrumD = Normalize(result,centr-1);
        unfoldedbgspectrumD->Write("unfoldedbgspectrum_centr40_80");
    }

    file->Close();
  }
}
// this comment serves no purpose
