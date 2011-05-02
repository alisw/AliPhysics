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
  fCFContainers(NULL),
  fTemporaryObjects(NULL),
  fCorrelation(NULL),
  fBackground(NULL),
  fEfficiencyFunction(NULL),
  fWeightCharm(NULL),
  fInclusiveSpectrum(kTRUE),
  fDumpToFile(kFALSE),
  fUnSetCorrelatedErrors(kTRUE),
  fIPanaHadronBgSubtract(kFALSE),
  fIPanaCharmBgSubtract(kFALSE),
  fIPanaConversionBgSubtract(kFALSE),
  fIPanaNonHFEBgSubtract(kFALSE),
  fNonHFEbgMethod2(kFALSE),
  fNbDimensions(1),
  fNMCEvents(0),
  fNMCbgEvents(0),
  fStepMC(-1),
  fStepTrue(-1),
  fStepData(-1),
  fStepBeforeCutsV0(-1),
  fStepAfterCutsV0(-1),
  fStepGuessedUnfolding(-1),
  fNumberOfIterations(5),
  fChargeChoosen(-1),
  fNCentralityBinAtTheEnd(0),
  fHadronEffbyIPcut(NULL),
  fBeamType(0),
  fDebugLevel(0)
{
  //
  // Default constructor
  //

  for(Int_t k = 0; k < 20; k++){
    fNEvents[k] = 0;
    fLowBoundaryCentralityBinAtTheEnd[k] = 0;
    fHighBoundaryCentralityBinAtTheEnd[k] = 0;
  }
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
    Int_t kNdim;
    if(fBeamType==0) kNdim=3;
    if(fBeamType==1) kNdim=4;

    Int_t dims[kNdim];
    // Get the requested format
    if(fBeamType==0)
    {
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

    if(fBeamType==1)
    {
        // PbPb analysis; centrality as first dimension
      Int_t nbDimensions = fNbDimensions;
      fNbDimensions = fNbDimensions + 1;
	switch(nbDimensions){
	case 1:
	  dims[0] = 5;
	  dims[1] = 0;
	  break;
	case 2:
	  dims[0] = 5;
	  for(Int_t i = 0; i < 2; i++) dims[(i+1)] = i;
	  break;
	case 3:
	  dims[0] = 5;
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

  AliCFContainer *datacontainerD = GetSlicedContainer(datacontainer, fNbDimensions, dims, -1, fChargeChoosen);
  AliCFContainer *contaminationcontainerD = GetSlicedContainer(contaminationcontainer, fNbDimensions, dims, -1, fChargeChoosen);
  if((!datacontainerD) || (!contaminationcontainerD)) return kFALSE;

  SetContainer(datacontainerD,AliHFEspectrum::kDataContainer);
  SetContainer(contaminationcontainerD,AliHFEspectrum::kBackgroundData);

  if(bghfecontainer){
    // set nonHFE backgrounds using temporary stored in hadronicBackground, will be replaced with proper container
    AliCFContainer *bgcontainer = bghfecontainer->GetCFContainer("hadronicBackground");
    if(!bgcontainer) return kFALSE;
    AliCFContainer *bgcontainerD = GetSlicedContainer(bgcontainer, fNbDimensions, dims, -1, fChargeChoosen);
    if(!bgcontainerD) return kFALSE;
    SetContainer(bgcontainerD,AliHFEspectrum::kNonHFEBackgroundData);
  }

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
  AliCFContainer *mccontaineresdD = GetSlicedContainer(mccontaineresd, fNbDimensions, dims, source, fChargeChoosen);
  AliCFContainer *mccontainermcD = GetSlicedContainer(mccontainermc, fNbDimensions, dims, source, fChargeChoosen);
  if((!mccontaineresdD) || (!mccontainermcD)) return kFALSE;
  SetContainer(mccontainermcD,AliHFEspectrum::kMCContainerMC);
  SetContainer(mccontaineresdD,AliHFEspectrum::kMCContainerESD);

  // set charm, nonHFE container to estimate BG
  if(!fInclusiveSpectrum){
   source = 0; //charm
   mccontainermcD = GetSlicedContainer(mccontainermc, fNbDimensions, dims, source, fChargeChoosen);
   SetContainer(mccontainermcD,AliHFEspectrum::kMCContainerCharmMC);
   source = 2; //conversion
   mccontainermcD = GetSlicedContainer(mccontainermc, fNbDimensions, dims, source, fChargeChoosen);
   SetContainer(mccontainermcD,AliHFEspectrum::kMCContainerConversionMC);
   source = 3; //non HFE except for the conversion
   mccontainermcD = GetSlicedContainer(mccontainermc, fNbDimensions, dims, source, fChargeChoosen);
   SetContainer(mccontainermcD,AliHFEspectrum::kMCContainerNonHFEMC);
  }

  // MC container: correlation matrix
  THnSparseF *mccorrelation = 0x0;
  if(fInclusiveSpectrum) {
    if((fStepMC==(AliHFEcuts::kNcutStepsMCTrack + AliHFEcuts::kStepHFEcutsTRD + 2))) mccorrelation = mchfecontainer->GetCorrelationMatrix("correlationstepafterPID");
    else if((fStepMC==(AliHFEcuts::kNcutStepsMCTrack + AliHFEcuts::kStepHFEcutsTRD + 1))) mccorrelation = mchfecontainer->GetCorrelationMatrix("correlationstepafterPID");
    else if((fStepMC==(AliHFEcuts::kNcutStepsMCTrack + AliHFEcuts::kStepHFEcutsTRD))) mccorrelation = mchfecontainer->GetCorrelationMatrix("correlationstepbeforePID");
    else if((fStepMC==(AliHFEcuts::kNcutStepsMCTrack + AliHFEcuts::kStepHFEcutsTRD - 1))) mccorrelation = mchfecontainer->GetCorrelationMatrix("correlationstepbeforePID");
    else mccorrelation = mchfecontainer->GetCorrelationMatrix("correlationstepafterPID");
  }
  else mccorrelation = mchfecontainer->GetCorrelationMatrix("correlationstepafterPID"); // we confirmed that we get same result by using it instead of correlationstepafterDE
  //else mccorrelation = mchfecontainer->GetCorrelationMatrix("correlationstepafterDE");
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

  printf("Steps are: stepdata %d, stepMC %d, steptrue %d, stepV0after %d, stepV0before %d\n",fStepData,fStepMC,fStepTrue,fStepAfterCutsV0,fStepBeforeCutsV0);

  ///////////////////////////
  // Check initialization
  ///////////////////////////

  if((!GetContainer(kDataContainer)) || (!GetContainer(kMCContainerMC)) || (!GetContainer(kMCContainerESD))){
    AliInfo("You have to init before");
    return kFALSE;
  }
  
  if((fStepTrue == -1) && (fStepMC == -1) && (fStepData == -1)) {
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

  ////////////////////////////////////////////////
  // Correct for TPC efficiency from V0
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

   
   

    Int_t ptpr;
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

    TH1D unfoldingspectrac[fNCentralityBinAtTheEnd];
    TGraphErrors unfoldingspectracn[fNCentralityBinAtTheEnd];
    TH1D correctedspectrac[fNCentralityBinAtTheEnd];
    TGraphErrors correctedspectracn[fNCentralityBinAtTheEnd];

    

    if(fBeamType==1) {

      TCanvas * ccorrected_allspectra = new TCanvas("correctedallspectra","correctedallspectra",1000,700);
      ccorrected_allspectra->Divide(2,1);
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
	  //printf("Number of events %d in the bin %d added!!!\n",fNEvents[k],k);
	  nbEvents += fNEvents[k];
	}
	//Double_t lowedgega = cenaxisa->GetBinLowEdge(fLowBoundaryCentralityBinAtTheEnd[binc]+1);
	//Double_t upedgega = cenaxisa->GetBinUpEdge(fHighBoundaryCentralityBinAtTheEnd[binc]);
	//printf("Bin Low edge %f, up edge %f for a\n",lowedgega,upedgega);
	//Double_t lowedgegb = cenaxisb->GetBinLowEdge(fLowBoundaryCentralityBinAtTheEnd[binc]+1);
	//Double_t upedgegb = cenaxisb->GetBinUpEdge(fHighBoundaryCentralityBinAtTheEnd[binc]);
	//printf("Bin Low edge %f, up edge %f for b\n",lowedgegb,upedgegb);
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
	ccorrected_allspectra->cd(1);
	gPad->SetLogy();
	TH1D *aftersunn = (TH1D *) aftersun->Clone();
	aftersunn->SetMarkerStyle(stylee[binc]);
	aftersunn->SetMarkerColor(colorr[binc]);
	if(binc==0) aftersunn->Draw("AP");
	else aftersunn->Draw("P");
	legtotal->AddEntry(aftersunn,(const char*) titlee,"p");
	ccorrected_allspectra->cd(2);
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

      ccorrected_allspectra->cd(1);
      legtotal->Draw("same");
      ccorrected_allspectra->cd(2);
      legtotalg->Draw("same");
      
      cenaxisa->SetRange(0,nbbin);
      cenaxisb->SetRange(0,nbbin);

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
      for(Int_t binc = 0; binc < fNCentralityBinAtTheEnd; binc++){
	unfoldingspectrac[binc].Write();
	unfoldingspectracn[binc].Write();
	correctedspectrac[binc].Write();
	correctedspectracn[binc].Write();
      }
      out->Close(); delete out;
    }

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
  if(subtractcontamination) {
    dataspectrumaftersubstraction = SubtractBackground(kTRUE);
    dataGridAfterFirstSteps = dataspectrumaftersubstraction;
  }

  /////////////////////////////////////////////////////////////////////////////////////////
  // Correct for IP efficiency for beauty electrons after subtracting all the backgrounds
  /////////////////////////////////////////////////////////////////////////////////////////


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
    if(fIPanaHadronBgSubtract){
      // Hadron background
      printf("Hadron background for IP analysis subtracted!\n");
      Int_t* bins=new Int_t[1];
      TH1D* htemp  = (TH1D *) fHadronEffbyIPcut->Projection(0);
      bins[0]=htemp->GetNbinsX();
      AliCFDataGrid *hbgContainer = new AliCFDataGrid("hbgContainer","hadron bg after IP cut",1,bins);
      hbgContainer->SetGrid(fHadronEffbyIPcut);
      backgroundGrid->Multiply(hbgContainer,1);
      spectrumSubtracted->Add(backgroundGrid,-1.0);
    }
    if(fIPanaCharmBgSubtract){
      // Charm background
      printf("Charm background for IP analysis subtracted!\n");
      AliCFDataGrid *charmbgContainer = (AliCFDataGrid *) GetCharmBackground();
      spectrumSubtracted->Add(charmbgContainer,-1.0);
    }
    if(fIPanaConversionBgSubtract){
      // Conversion background
      AliCFDataGrid *conversionbgContainer = (AliCFDataGrid *) GetConversionBackground();
      spectrumSubtracted->Add(conversionbgContainer,-1.0);
      printf("Conversion background subtraction is preliminary!\n");
    }
    if(fIPanaNonHFEBgSubtract){
      // NonHFE background
      AliCFDataGrid *nonHFEbgContainer = (AliCFDataGrid *) GetNonHFEBackground();
      spectrumSubtracted->Add(nonHFEbgContainer,-1.0);
      printf("Non HFE background subtraction is preliminary!\n");
    }
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

    Int_t ptpr;
    if(fBeamType==0) ptpr=0;
    if(fBeamType==1) ptpr=1;
    
    TCanvas * cbackgroundsubtraction = new TCanvas("backgroundsubtraction","backgroundsubtraction",1000,700);
    cbackgroundsubtraction->Divide(3,1);
    cbackgroundsubtraction->cd(1);
    gPad->SetLogy();
    TH1D *measuredTH1Daftersubstraction = (TH1D *) spectrumSubtracted->Project(ptpr);
    TH1D *measuredTH1Dbeforesubstraction = (TH1D *) dataspectrumbeforesubstraction->Project(ptpr);
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
    TH1D *measuredTH1background = (TH1D *) backgroundGrid->Project(ptpr);
    CorrectFromTheWidth(measuredTH1background);
    measuredTH1background->SetStats(0);
    measuredTH1background->SetTitle("");
    measuredTH1background->GetYaxis()->SetTitle("dN/dp_{T} [(GeV/c)^{-1}]");
    measuredTH1background->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    measuredTH1background->SetMarkerStyle(26);
    measuredTH1background->SetMarkerColor(kBlack);
    measuredTH1background->SetLineColor(kBlack);
    measuredTH1background->Draw();

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
      
    }
   
       
  }
  

  return spectrumSubtracted;
}

//____________________________________________________________
AliCFDataGrid* AliHFEspectrum::GetCharmBackground(){
  //
  // calculate charm background
  //

  Double_t evtnorm=0;
  if(fNMCEvents) evtnorm= double(fNEvents[0])/double(fNMCEvents);

  AliCFContainer *mcContainer = GetContainer(kMCContainerCharmMC);
  if(!mcContainer){
    AliError("MC Container not available");
    return NULL;
  }

  if(!fCorrelation){
    AliError("No Correlation map available");
    return NULL;
  }

  Int_t nDim = 1;
  AliCFDataGrid *charmBackgroundGrid= 0x0;
  charmBackgroundGrid = new AliCFDataGrid("charmBackgroundGrid","charmBackgroundGrid",*mcContainer, fStepMC+2);
  TH1D *tmphisto = (TH1D *) charmBackgroundGrid->Project(0);
  Int_t* bins=new Int_t[1];
  bins[0]=tmphisto->GetNbinsX();
  AliCFDataGrid *weightedCharmContainer = new AliCFDataGrid("weightedCharmContainer","weightedCharmContainer",1,bins);
  weightedCharmContainer->SetGrid(GetWeights());
  charmBackgroundGrid->Multiply(weightedCharmContainer,evtnorm);

  // Efficiency (set efficiency to 1 for only folding) 
  AliCFEffGrid* efficiencyD = new AliCFEffGrid("efficiency","",*mcContainer);
  efficiencyD->CalculateEfficiency(0,0);

  // Folding 
  AliCFUnfolding folding("unfolding","",nDim,fCorrelation,efficiencyD->GetGrid(),charmBackgroundGrid->GetGrid(),charmBackgroundGrid->GetGrid());
  folding.SetMaxNumberOfIterations(1);
  folding.Unfold();

  // Results
  THnSparse* result1= folding.GetEstMeasured(); // folded spectra
  THnSparse* result=(THnSparse*)result1->Clone();
  charmBackgroundGrid->SetGrid(result);

  return charmBackgroundGrid;
}

//____________________________________________________________
AliCFDataGrid* AliHFEspectrum::GetConversionBackground(){
  //
  // calculate conversion background
  //

  Double_t evtnorm[1];
  if(fNMCbgEvents) evtnorm[0]= double(fNEvents[0])/double(fNMCbgEvents);
  printf("check event!!! %lf \n",evtnorm[0]);

  // Background Estimate
  AliCFContainer *backgroundContainer = GetContainer(kNonHFEBackgroundData);
  if(!backgroundContainer){
    AliError("MC background container not found");
    return NULL;
  }

  Int_t stepbackground = 2;
  AliCFDataGrid *backgroundGrid = new AliCFDataGrid("ConversionBgGrid","ConversionBgGrid",*backgroundContainer,stepbackground);
  backgroundGrid->Scale(evtnorm);

  new TCanvas;
  TH1D *histodraw111 = (TH1D *) backgroundGrid->Project(0);
  CorrectFromTheWidth(histodraw111);
  histodraw111->SetMarkerStyle(20);
  histodraw111->SetLineColor(4);
  histodraw111->SetMarkerColor(4);
  histodraw111->Draw("p");


  return backgroundGrid;
}


//____________________________________________________________
AliCFDataGrid* AliHFEspectrum::GetNonHFEBackground(){
  //
  // calculate non HFE background
  //

  if(fNonHFEbgMethod2){
    Double_t evtnorm=0;
    if(fNMCEvents) evtnorm= double(fNEvents[0])/double(fNMCEvents);

    AliCFContainer *mcContainer = GetContainer(kMCContainerNonHFEMC);
    if(!mcContainer){
      AliError("MC Container not available");
      return NULL;
    }

    if(!fCorrelation){
      AliError("No Correlation map available");
      return NULL;
    }

    Int_t nDim = 1;
    AliCFDataGrid *nonHFEBackgroundGrid= 0x0;
    nonHFEBackgroundGrid = new AliCFDataGrid("nonHFEBackgroundGrid","nonHFEBackgroundGrid",*mcContainer, fStepMC+2);
    TH1D *tmphisto = (TH1D *) nonHFEBackgroundGrid->Project(0);
    Int_t* bins=new Int_t[1];
    bins[0]=tmphisto->GetNbinsX();
    AliCFDataGrid *weightedNonHFEContainer = new AliCFDataGrid("weightedNonHFEContainer","weightedNonHFEContainer",1,bins);
    weightedNonHFEContainer->SetGrid(GetWeights()); // we have to use different weithing function
    nonHFEBackgroundGrid->Multiply(weightedNonHFEContainer,evtnorm);

    // Efficiency (set efficiency to 1 for only folding) 
    AliCFEffGrid* efficiencyD = new AliCFEffGrid("efficiency","",*mcContainer);
    efficiencyD->CalculateEfficiency(0,0);

    // Folding 
    AliCFUnfolding folding("unfolding","",nDim,fCorrelation,efficiencyD->GetGrid(),nonHFEBackgroundGrid->GetGrid(),nonHFEBackgroundGrid->GetGrid());
    folding.SetMaxNumberOfIterations(1);
    folding.Unfold();

    // Results
    THnSparse* result1= folding.GetEstMeasured(); // folded spectra
    THnSparse* result=(THnSparse*)result1->Clone();
    nonHFEBackgroundGrid->SetGrid(result);

    return nonHFEBackgroundGrid;
  }
  else{
    Double_t evtnorm[1];
    if(fNMCbgEvents) evtnorm[0]= double(fNEvents[0])/double(fNMCbgEvents);

    // Background Estimate
    AliCFContainer *backgroundContainer = GetContainer(kNonHFEBackgroundData);
    if(!backgroundContainer){
      AliError("MC background container not found");
      return NULL;
    }

    Int_t stepbackground = 3;
    AliCFDataGrid *backgroundGrid = new AliCFDataGrid("NonHFEBgGrid","NonHFEBgGrid",*backgroundContainer,stepbackground);
    backgroundGrid->Scale(evtnorm);

    new TCanvas;
    TH1D *histodraw222 = (TH1D *) backgroundGrid->Project(0);
    CorrectFromTheWidth(histodraw222);
    histodraw222->SetMarkerStyle(20);
    histodraw222->SetLineColor(4);
    histodraw222->SetMarkerColor(4);
    histodraw222->Draw("p");

    return backgroundGrid;
  }

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

    Int_t ptpr;
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

  // Unfold 
  
  AliCFUnfolding unfolding("unfolding","",fNbDimensions,fCorrelation,efficiencyD->GetGrid(),dataGrid->GetGrid(),guessedTHnSparse);
  //unfolding.SetUseCorrelatedErrors();
  if(fUnSetCorrelatedErrors) unfolding.UnsetCorrelatedErrors();
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

    Int_t ptpr;
    if(fBeamType==0) ptpr=0;
    if(fBeamType==1) ptpr=1;
    
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

    Int_t ptpr;
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

    Int_t ptpr;
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

    Int_t ptpr;
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
  if(fNEvents[i] > 0) {

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
      nCorr = 0.5 * 1/1.6 * 1./(Double_t)(fNEvents[i]) * 1/(2. * TMath::Pi() * p) * n;
      errdN = 1./(2. * TMath::Pi() * p);
      errdp = 1./(2. * TMath::Pi() * p*p) * n;
      dNcorr = 0.5 * 1/1.6 * 1./(Double_t)(fNEvents[i]) * TMath::Sqrt(errdN * errdN * dN *dN + errdp *errdp * dp *dp);
      
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
  if(normalization > 0) {

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
      nCorr = 0.5 * 1/1.6 * 1./(Double_t)(normalization) * 1/(2. * TMath::Pi() * p) * n;
      errdN = 1./(2. * TMath::Pi() * p);
      errdp = 1./(2. * TMath::Pi() * p*p) * n;
      dNcorr = 0.5 * 1/1.6 * 1./(Double_t)(normalization) * TMath::Sqrt(errdN * errdN * dN *dN + errdp *errdp * dp *dp);
      
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
AliCFContainer *AliHFEspectrum::GetSlicedContainer(AliCFContainer *container, Int_t nDim, Int_t *dimensions,Int_t source,Int_t positivenegative) {
  //
  // Slice bin for a given source of electron
  // nDim is the number of dimension the corrections are done
  // dimensions are the definition of the dimensions
  // source is if we want to keep only one MC source (-1 means we don't cut on the MC source)
  // positivenegative if we want to keep positive (1) or negative (0) or both (-1)
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
	varMin[ivar] = binLimits[source];
	varMax[ivar] = binLimits[source];
      }     
    }
    // charge
    if(ivar == 3) {
      if((positivenegative>= 0) && (positivenegative<container->GetNBins(ivar))) {
	varMin[ivar] = binLimits[positivenegative];
	varMax[ivar] = binLimits[positivenegative];
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
THnSparseF *AliHFEspectrum::GetSlicedCorrelation(THnSparseF *correlationmatrix, Int_t nDim, Int_t *dimensions) const {
  //
  // Slice correlation
  //

  Int_t ndimensions = correlationmatrix->GetNdimensions();
  //printf("Number of dimension %d correlation map\n",ndimensions);
  if(ndimensions < (2*nDim)) {
    AliError("Problem in the dimensions");
    return NULL;
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
THnSparse* AliHFEspectrum::GetWeights(){
  
  //
  // Measured D->e based weighting factors
  //

  const Int_t nDim=1;
  Int_t nBin[nDim] = {44};
  const Double_t kPtbound[2] = {0.1, 20.};

  Double_t* binEdges[nDim];
  binEdges[0] =  AliHFEtools::MakeLogarithmicBinning(nBin[0], kPtbound[0], kPtbound[1]);

  fWeightCharm = new THnSparseF("weightHisto", "weiting factor; pt[GeV/c]", nDim, nBin);
  for(Int_t idim = 0; idim < nDim; idim++)
     fWeightCharm->SetBinEdges(idim, binEdges[idim]);
  /* //fit  
  Double_t weight[44]={ 2.20939,2.1742,2.13569,2.09369,2.04806,1.99869,1.94552,1.88854,1.82782,1.76349,1.69577,1.62497,1.55149,1.4758,1.39848,1.32016,1.24153,1.16328,1.08615,1.01081,0.937916,0.868033,0.801642,0.739117,0.680724,0.626621,0.576859,0.5314,0.490125,0.452854,0.419361,0.389388,0.362661,0.338898,0.317822,0.299165,0.282674,0.268115,0.255271,0.243946,0.233964,0.225166,0.217413,0.210578};*/
  //points
  Double_t weight[44]={ 0.443336,0.437069,0.424299,0.41365,0.410649,0.385012,0.364252,0.360181,0.347073,0.354445,0.369578,0.358253,0.381115,0.412863,0.451792,0.47061,0.544857,0.566261,0.656843,0.669739,0.725247,0.779511,0.750385,0.777815,0.673391,0.694613,0.598251,0.513847,0.45729,0.468098,0.379422,0.358718,0.403793,0.392902,0.349981,0.27205,0.302447,0.279313,0.258986,0.14738,0.414487,0.0552302,0.0184101,0};
  Double_t pt[44]={0.106398,0.120014,0.135371,0.152694,0.172234,0.194274,0.219135,0.247177,0.278807,0.314485,0.354729,0.400122,0.451324,0.509078,0.574223,0.647704,0.730589,0.824079,0.929534,1.04848,1.18265,1.33399,1.5047,1.69725,1.91444,2.15943,2.43576,2.74745,3.09904,3.49561,3.94293,4.44749,5.01662,5.65858,6.38268,7.19945,8.12074,9.15992,10.3321,11.6542,13.1456,14.8278,16.7252,18.8655};
  Double_t dataE[1];
  for(int i=0; i<44; i++){
    dataE[0]=pt[i];
    fWeightCharm->Fill(dataE,weight[i]);
  }
  Int_t* ibins[nDim];
  for(Int_t ivar = 0; ivar < nDim; ivar++)
    ibins[ivar] = new Int_t[nBin[ivar] + 1];
  fWeightCharm->SetBinError(ibins[0],0);

  return fWeightCharm;
}
