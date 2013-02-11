 
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

#include "AliHFEBeautySpectrum.h"
#include "AliHFEBeautySpectrumQA.h"
#include "AliHFEcuts.h"
#include "AliHFEcontainer.h"
#include "AliHFEtools.h"

ClassImp(AliHFEBeautySpectrum)

//____________________________________________________________
AliHFEBeautySpectrum::AliHFEBeautySpectrum(const char *name):
AliHFECorrectSpectrumBase(name),
  fQA(NULL),
  fTemporaryObjects(NULL),
  fBackground(NULL),
  fEfficiencyTOFPIDD(NULL),
  fEfficiencyesdTOFPIDD(NULL),
  fEfficiencyIPCharmD(NULL),
  fEfficiencyIPBeautyD(NULL),
  fEfficiencyIPBeautyesdD(NULL),
  fEfficiencyIPConversionD(NULL),
  fEfficiencyIPNonhfeD(NULL), 
  fNonHFEbg(NULL),
  fWeightCharm(NULL),
  fInclusiveSpectrum(kFALSE),
  fDumpToFile(kFALSE),
  fUnSetCorrelatedErrors(kTRUE),
  fIPanaHadronBgSubtract(kFALSE),
  fIPanaCharmBgSubtract(kFALSE),
  fIPanaNonHFEBgSubtract(kFALSE),
  fIPParameterizedEff(kFALSE),
  fNonHFEsyst(0),
  fBeauty2ndMethod(kFALSE),
  fIPEffCombinedSamples(kTRUE),
  fNMCEvents(0),
  fNMCbgEvents(0),
  fNCentralityBinAtTheEnd(0),
  fFillMoreCorrelationMatrix(kFALSE),
  fHadronEffbyIPcut(NULL),
  fEfficiencyCharmSigD(NULL),
  fEfficiencyBeautySigD(NULL),
  fEfficiencyBeautySigesdD(NULL),
  fConversionEff(NULL),
  fNonHFEEff(NULL),
  fCharmEff(NULL),
  fBeautyEff(NULL),
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

  fQA = new AliHFEBeautySpectrumQA("AliHFEBeautySpectrumQA");

  for(Int_t k = 0; k < 20; k++){
    fLowBoundaryCentralityBinAtTheEnd[k] = 0;
    fHighBoundaryCentralityBinAtTheEnd[k] = 0;
  }
  fNMCbgEvents = 0;  
  fEfficiencyTOFPIDD = 0;
  fEfficiencyesdTOFPIDD = 0;
  fEfficiencyIPCharmD = 0;     
  fEfficiencyIPBeautyD = 0;    
  fEfficiencyIPBeautyesdD = 0;
  fEfficiencyIPConversionD = 0;
  fEfficiencyIPNonhfeD = 0;   
  
  fConversionEff = 0;
  fNonHFEEff = 0;
  fCharmEff = 0;
  fBeautyEff = 0;
  fEfficiencyCharmSigD = 0;
  fEfficiencyBeautySigD = 0;
  fEfficiencyBeautySigesdD = 0;
  
  
  memset(fConvSourceContainer, 0, sizeof(AliCFContainer*) * kElecBgSources * kBgLevels);
  memset(fNonHFESourceContainer, 0, sizeof(AliCFContainer*) * kElecBgSources * kBgLevels);
}
//____________________________________________________________
AliHFEBeautySpectrum::AliHFEBeautySpectrum(const AliHFEBeautySpectrum &ref):
  AliHFECorrectSpectrumBase(ref),
  fQA(ref.fQA),
  fTemporaryObjects(ref.fTemporaryObjects),
  fBackground(ref.fBackground),
  fEfficiencyTOFPIDD(ref.fEfficiencyTOFPIDD),
  fEfficiencyesdTOFPIDD(ref.fEfficiencyesdTOFPIDD),
  fEfficiencyIPCharmD(ref.fEfficiencyIPCharmD),
  fEfficiencyIPBeautyD(ref.fEfficiencyIPBeautyD),
  fEfficiencyIPBeautyesdD(ref.fEfficiencyIPBeautyesdD),
  fEfficiencyIPConversionD(ref.fEfficiencyIPConversionD),
  fEfficiencyIPNonhfeD(ref.fEfficiencyIPNonhfeD), 
  fNonHFEbg(ref.fNonHFEbg),
  fWeightCharm(ref.fWeightCharm),
  fInclusiveSpectrum(ref.fInclusiveSpectrum),
  fDumpToFile(ref.fDumpToFile),
  fUnSetCorrelatedErrors(ref.fUnSetCorrelatedErrors),
  fIPanaHadronBgSubtract(ref.fIPanaHadronBgSubtract),
  fIPanaCharmBgSubtract(ref.fIPanaCharmBgSubtract),
  fIPanaNonHFEBgSubtract(ref.fIPanaNonHFEBgSubtract),
  fIPParameterizedEff(ref.fIPParameterizedEff),
  fNonHFEsyst(ref.fNonHFEsyst),
  fBeauty2ndMethod(ref.fBeauty2ndMethod),
  fIPEffCombinedSamples(ref.fIPEffCombinedSamples),
  fNMCEvents(ref.fNMCEvents),
  fNMCbgEvents(ref.fNMCbgEvents),
  fNCentralityBinAtTheEnd(ref.fNCentralityBinAtTheEnd),
  fFillMoreCorrelationMatrix(ref.fFillMoreCorrelationMatrix),
  fHadronEffbyIPcut(ref.fHadronEffbyIPcut),
  fEfficiencyCharmSigD(ref.fEfficiencyCharmSigD),
  fEfficiencyBeautySigD(ref.fEfficiencyBeautySigD),
  fEfficiencyBeautySigesdD(ref.fEfficiencyBeautySigesdD),
  fConversionEff(ref.fConversionEff),
  fNonHFEEff(ref.fNonHFEEff),
  fCharmEff(ref.fCharmEff),
  fBeautyEff(ref.fBeautyEff),
  fConversionEffbgc(ref.fConversionEffbgc),
  fNonHFEEffbgc(ref.fNonHFEEffbgc),      
  fBSpectrum2ndMethod(ref.fBSpectrum2ndMethod),
  fkBeauty2ndMethodfilename(ref.fkBeauty2ndMethodfilename),
  fBeamType(ref.fBeamType),
  fEtaSyst(ref.fEtaSyst),
  fDebugLevel(ref.fDebugLevel),
  fWriteToFile(ref.fWriteToFile),
  fUnfoldBG(ref.fUnfoldBG)
{
  //
  // Copy constructor
  //

  ref.Copy(*this);  
}
//____________________________________________________________
AliHFEBeautySpectrum &AliHFEBeautySpectrum::operator=(const AliHFEBeautySpectrum &ref){
 //
  // Assignment operator
  //
  if(this == &ref) 
    ref.Copy(*this);
  return *this;
}
//____________________________________________________________
void AliHFEBeautySpectrum::Copy(TObject &o) const {
  // 
  // Copy into object o
  //
  AliHFEBeautySpectrum &target = dynamic_cast<AliHFEBeautySpectrum &>(o);
  target.fQA = fQA;

  target.fQA = fQA;
  target.fTemporaryObjects = fTemporaryObjects;
  target.fBackground = fBackground;
  target.fWeightCharm = fWeightCharm;
  target.fDumpToFile = fDumpToFile;
  target.fUnSetCorrelatedErrors = fUnSetCorrelatedErrors;
  target.fIPanaHadronBgSubtract = fIPanaHadronBgSubtract;
  target.fIPanaCharmBgSubtract = fIPanaCharmBgSubtract;
  target.fIPanaNonHFEBgSubtract = fIPanaNonHFEBgSubtract;
  target.fIPParameterizedEff = fIPParameterizedEff;
  target.fNonHFEsyst = fNonHFEsyst;
  target.fBeauty2ndMethod = fBeauty2ndMethod;
  target.fIPEffCombinedSamples = fIPEffCombinedSamples;
  target.fNMCEvents = fNMCEvents;
  target.fNCentralityBinAtTheEnd = fNCentralityBinAtTheEnd;
  target.fFillMoreCorrelationMatrix = fFillMoreCorrelationMatrix;
  target.fHadronEffbyIPcut = fHadronEffbyIPcut;
  target.fConversionEffbgc = fConversionEffbgc;
  target.fNonHFEEffbgc = fNonHFEEffbgc;      
  target.fBSpectrum2ndMethod = fBSpectrum2ndMethod;
  target.fkBeauty2ndMethodfilename = fkBeauty2ndMethodfilename;
  target.fBeamType = fBeamType;
  target.fEtaSyst = fEtaSyst;
  target.fDebugLevel = fDebugLevel;
  target.fWriteToFile = fWriteToFile;
  target.fUnfoldBG = fUnfoldBG;
}
//____________________________________________________________
AliHFEBeautySpectrum::~AliHFEBeautySpectrum(){
  //
  // Destructor
  //
  if(fTemporaryObjects){
    fTemporaryObjects->Clear();
    delete fTemporaryObjects;
  }
  if(fQA) delete fQA;
}
//____________________________________________________________
Bool_t AliHFEBeautySpectrum::Init(const AliHFEcontainer *datahfecontainer, const AliHFEcontainer *mchfecontainer, const AliHFEcontainer *bghfecontainer, const AliHFEcontainer */*v0hfecontainer*/,AliCFContainer */*photoniccontainerD*/){
  //
  // Init what we need for the correction:
  //
  // Raw spectrum, hadron contamination
  // MC efficiency maps, correlation matrix
  //
  // This for a given dimension.
  // If no inclusive spectrum, then take only efficiency map for beauty electron
  // and the appropriate correlation matrix
  //

  
  if(fBeauty2ndMethod) CallInputFileForBeauty2ndMethod();

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
  AliCFContainer *datacontainer = datahfecontainer->MakeMergedCFContainer("sumreco","sumreco","recTrackContReco:recTrackContDEReco");
  AliCFContainer *contaminationcontainer = datahfecontainer->GetCFContainer("hadronicBackground");
  if((!datacontainer) || (!contaminationcontainer)) return kFALSE;
  AliCFContainer *datacontainerD = GetSlicedContainer(datacontainer, fNbDimensions, dims, -1, fChargeChoosen,fTestCentralityLow,fTestCentralityHigh);
  AliCFContainer *contaminationcontainerD = GetSlicedContainer(contaminationcontainer, fNbDimensions, dims, -1, fChargeChoosen,fTestCentralityLow,fTestCentralityHigh);
  if((!datacontainerD) || (!contaminationcontainerD)) return kFALSE;
  SetContainer(datacontainerD,AliHFECorrectSpectrumBase::kDataContainer);
  SetContainer(contaminationcontainerD,AliHFECorrectSpectrumBase::kBackgroundData);

  // QA : see with MJ if it is necessary for Beauty
  Int_t dimqa = datacontainer->GetNVar();
  Int_t dimsqa[dimqa];
  for(Int_t i = 0; i < dimqa; i++) dimsqa[i] = i;
  AliCFContainer *datacontainerDQA = GetSlicedContainer(datacontainer, dimqa, dimsqa, -1, fChargeChoosen,fTestCentralityLow,fTestCentralityHigh);
  fQA->AddResultAt(datacontainerDQA,AliHFEBeautySpectrumQA::kDataProjection);

  //
  // MC container: ESD/MC efficiency maps + MC/MC efficiency maps
  // 
  AliCFContainer *mccontaineresd = 0x0;
  AliCFContainer *mccontaineresdbg = 0x0;
  AliCFContainer *mccontainermc = 0x0;
  AliCFContainer *mccontainermcbg = 0x0;
  AliCFContainer *nonHFEweightedContainer = 0x0;
  AliCFContainer *convweightedContainer = 0x0;
  AliCFContainer *nonHFEweightedContainerSig = 0x0;//mjmj
  AliCFContainer *convweightedContainerSig = 0x0;//mjmj
  AliCFContainer *nonHFEtempContainer = 0x0;//temporary container to be sliced for the fnonHFESourceContainers
  AliCFContainer *convtempContainer = 0x0;//temporary container to be sliced for the fConvSourceContainers
   
  mccontaineresd = mchfecontainer->MakeMergedCFContainer("sumesd","sumesd","MCTrackCont:recTrackContReco:recTrackContDEReco");
  mccontaineresdbg = bghfecontainer->MakeMergedCFContainer("sumesd","sumesd","MCTrackCont:recTrackContReco:recTrackContDEReco");
  mccontainermc = mchfecontainer->MakeMergedCFContainer("summc","summc","MCTrackCont:recTrackContMC:recTrackContDEMC");
  mccontainermcbg = bghfecontainer->MakeMergedCFContainer("summcbg","summcbg","MCTrackCont:recTrackContMC:recTrackContDEMC");
  
  if(fNonHFEsyst){   
    const Char_t *sourceName[kElecBgSources]={"Pion","Eta","Omega","Phi","EtaPrime","Rho","K0sSec","OtherSecPi0","Else"};
    const Char_t *levelName[kBgLevels]={"Best","Lower","Upper"};
    for(Int_t iSource = 0; iSource < kElecBgSources; iSource++){
      for(Int_t iLevel = 0; iLevel < kBgLevels; iLevel++){
        if(!(nonHFEtempContainer = bghfecontainer->GetCFContainer(Form("mesonElecs%s%s",sourceName[iSource],levelName[iLevel]))))continue;
        convtempContainer = bghfecontainer->GetCFContainer(Form("conversionElecs%s%s",sourceName[iSource],levelName[iLevel]));     
        if(!(fConvSourceContainer[iSource][iLevel] = GetSlicedContainer(convtempContainer, fNbDimensions, dims, -1, fChargeChoosen,fTestCentralityLow,fTestCentralityHigh)))continue;
        if(!(fNonHFESourceContainer[iSource][iLevel] = GetSlicedContainer(nonHFEtempContainer, fNbDimensions, dims, -1, fChargeChoosen,fTestCentralityLow,fTestCentralityHigh)))continue;
        //	      if((!fConvSourceContainer[iSource][iLevel][icentr])||(!fNonHFESourceContainer[iSource][iLevel])) return kFALSE;
        if(fBeamType == 1)break;//no systematics yet for PbPb non-HFE backgrounds
      }
    }
  }
  // else{      
  nonHFEweightedContainer = bghfecontainer->GetCFContainer("mesonElecs");
  convweightedContainer = bghfecontainer->GetCFContainer("conversionElecs");
  if((!convweightedContainer)||(!nonHFEweightedContainer)) return kFALSE; 
  nonHFEweightedContainerSig = mchfecontainer->GetCFContainer("mesonElecs");//mjmj
  convweightedContainerSig = mchfecontainer->GetCFContainer("conversionElecs");//mjmj
  if((!convweightedContainerSig)||(!nonHFEweightedContainerSig)) return kFALSE;
  //}
  
  if((!mccontaineresd) || (!mccontainermc)) return kFALSE;  
  
  Int_t source = 1; //beauty
  AliCFContainer *mccontaineresdD = GetSlicedContainer(mccontaineresd, fNbDimensions, dims, source, fChargeChoosen,fTestCentralityLow,fTestCentralityHigh);
  AliCFContainer *mccontainermcD = GetSlicedContainer(mccontainermc, fNbDimensions, dims, source, fChargeChoosen,fTestCentralityLow,fTestCentralityHigh);
  if((!mccontaineresdD) || (!mccontainermcD)) return kFALSE;
  SetContainer(mccontainermcD,AliHFECorrectSpectrumBase::kMCContainerMC);
  SetContainer(mccontaineresdD,AliHFECorrectSpectrumBase::kMCContainerESD);

  // set charm, nonHFE container to estimate BG
  source = 0; //charm
  mccontainermcD = GetSlicedContainer(mccontainermc, fNbDimensions, dims, source, fChargeChoosen,fTestCentralityLow,fTestCentralityHigh); //mjmjmj
  //mccontainermcD = GetSlicedContainer(mccontainermcbg, fNbDimensions, dims, source, fChargeChoosen,fTestCentralityLow,fTestCentralityHigh);
  SetContainer(mccontainermcD,AliHFECorrectSpectrumBase::kMCContainerCharmMC);

  //if(!fNonHFEsyst){
  AliCFContainer *nonHFEweightedContainerD = GetSlicedContainer(nonHFEweightedContainer, fNbDimensions, dims, -1, fChargeChoosen,fTestCentralityLow,fTestCentralityHigh);
  SetContainer(nonHFEweightedContainerD,AliHFECorrectSpectrumBase::kMCWeightedContainerNonHFEESD);
  AliCFContainer *convweightedContainerD = GetSlicedContainer(convweightedContainer, fNbDimensions, dims, -1, fChargeChoosen,fTestCentralityLow,fTestCentralityHigh);
  SetContainer(convweightedContainerD,AliHFECorrectSpectrumBase::kMCWeightedContainerConversionESD);
  
  //mjmj
  AliCFContainer *nonHFEweightedContainerSigD = GetSlicedContainer(nonHFEweightedContainerSig, fNbDimensions, dims, -1, fChargeChoosen,fTestCentralityLow,fTestCentralityHigh);
  SetContainer(nonHFEweightedContainerSigD,AliHFECorrectSpectrumBase::kMCWeightedContainerNonHFEESDSig);
  if(!GetContainer(AliHFECorrectSpectrumBase::kMCWeightedContainerNonHFEESDSig)){printf("merged non-HFE container not found!\n");return kFALSE;};
  AliCFContainer *convweightedContainerSigD = GetSlicedContainer(convweightedContainerSig, fNbDimensions, dims, -1, fChargeChoosen,fTestCentralityLow,fTestCentralityHigh);
  SetContainer(convweightedContainerSigD,AliHFECorrectSpectrumBase::kMCWeightedContainerConversionESDSig); 
  if(!GetContainer(AliHFECorrectSpectrumBase::kMCWeightedContainerConversionESDSig)){printf(" merged conversion container not found!\n");return kFALSE;}; 
  //}
  
  SetParameterizedEff(mccontainermc, mccontainermcbg, mccontaineresd, mccontaineresdbg, dims);
  
  //
  // MC container: correlation matrix
  //
  THnSparseF *mccorrelation = mchfecontainer->GetCorrelationMatrix("correlationstepafterPID"); // we confirmed that we get same result by using it instead of correlationstepafterDE
  //mccorrelation = mchfecontainer->GetCorrelationMatrix("correlationstepafterDE");
  if(!mccorrelation) return kFALSE;
  Int_t testCentralityLow = fTestCentralityLow;
  Int_t testCentralityHigh = fTestCentralityHigh;
  if(fFillMoreCorrelationMatrix) {
    testCentralityLow = fTestCentralityLow-1;
    testCentralityHigh = fTestCentralityHigh+1;
  }
  THnSparseF *mccorrelationD = GetSlicedCorrelation(mccorrelation, fNbDimensions, dims, fChargeChoosen, testCentralityLow,testCentralityHigh);
  if(!mccorrelationD) {
    printf("No correlation\n");
    return kFALSE;
  }
  SetCorrelation(mccorrelationD);

  // QA
  fQA->AddResultAt(mccorrelation,AliHFEBeautySpectrumQA::kCMProjection); 
  //
  fQA->DrawProjections();


  return kTRUE;
}


//____________________________________________________________
void AliHFEBeautySpectrum::CallInputFileForBeauty2ndMethod(){
  //
  // get spectrum for beauty 2nd method
  //
  //
    TFile *inputfile2ndmethod=TFile::Open(fkBeauty2ndMethodfilename);
    fBSpectrum2ndMethod = new TH1D(*static_cast<TH1D*>(inputfile2ndmethod->Get("BSpectrum")));
}

//____________________________________________________________
Bool_t AliHFEBeautySpectrum::Correct(Bool_t subtractcontamination, Bool_t /*subtractphotonic*/){
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

  printf("Steps are: stepdata %d, stepMC %d, steptrue %d\n",fStepData,fStepMC,fStepTrue);

  ///////////////////////////
  // Check initialization
  ///////////////////////////

  if((!GetContainer(kDataContainer)) || (!GetContainer(kMCContainerMC)) || (!GetContainer(kMCContainerESD))){
    AliInfo("You have to init before");
    return kFALSE;
  }
  
  if((fStepTrue <= 0) && (fStepMC <= 0) && (fStepData <= 0)) {
    AliInfo("You have to set the steps before: SetMCTruthStep, SetMCEffStep, SetStepToCorrect");
    return kFALSE;
  }
 
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

  /////////////////////////////////////////////////////////////////////////////////////////
  // Correct for IP efficiency for beauty electrons after subtracting all the backgrounds
  /////////////////////////////////////////////////////////////////////////////////////////

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
  

  //////////
  // Plot
  //////////

  // QA final results
  TGraphErrors* correctedspectrumD = Normalize(correctedspectrum);
  TGraphErrors* alltogetherspectrumD = Normalize(alltogetherCorrection);
  fQA->AddResultAt(correctedspectrumD,AliHFEBeautySpectrumQA::kFinalResultUnfolded);
  fQA->AddResultAt(alltogetherspectrumD,AliHFEBeautySpectrumQA::kFinalResultDirectEfficiency);
  fQA->AddResultAt(correctedspectrum,AliHFEBeautySpectrumQA::kFinalResultUnfSparse);
  fQA->AddResultAt(alltogetherCorrection,AliHFEBeautySpectrumQA::kFinalResultDirectEffSparse);
  fQA->DrawResult();

 
  if(fBeamType == 0){
    if(fNonHFEsyst){
      CalculateNonHFEsyst();
    }
  }

  // Dump to file if needed
  
  if(fDumpToFile) {
    TFile *out;
    out = new TFile("finalSpectrum.root","update");
    out->cd();
    // to do centrality dependent
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
        
    fEfficiencyCharmSigD->SetTitle("IPEfficiencyForCharmSig");
    fEfficiencyCharmSigD->SetName("IPEfficiencyForCharmSig");
    fEfficiencyCharmSigD->Write();
    fEfficiencyBeautySigD->SetTitle("IPEfficiencyForBeautySig");
    fEfficiencyBeautySigD->SetName("IPEfficiencyForBeautySig");
    fEfficiencyBeautySigD->Write();
    fCharmEff->SetTitle("IPEfficiencyForCharm");
    fCharmEff->SetName("IPEfficiencyForCharm");
    fCharmEff->Write();
    fBeautyEff->SetTitle("IPEfficiencyForBeauty");
    fBeautyEff->SetName("IPEfficiencyForBeauty");
    fBeautyEff->Write();
    fConversionEff->SetTitle("IPEfficiencyForConversion");
    fConversionEff->SetName("IPEfficiencyForConversion");
    fConversionEff->Write();
    fNonHFEEff->SetTitle("IPEfficiencyForNonHFE");
    fNonHFEEff->SetName("IPEfficiencyForNonHFE");
    fNonHFEEff->Write();
    
   
    out->Close(); 
    delete out;
   }
  return kTRUE;
}
//____________________________________________________________
AliCFDataGrid* AliHFEBeautySpectrum::SubtractBackground(Bool_t setBackground){
  //
  // Apply background subtraction for IP analysis
  //
  
  Int_t nbins=1;
  printf("subtracting background\n");
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
  
  const Int_t bgPlots = 8;
  TH1D *incElec = 0x0;
  TH1D *hadron = 0x0;
  TH1D *charm = 0x0;
  TH1D *nonHFE[bgPlots];
  TH1D *subtracted = 0x0;
  
  THnSparseF* sparseIncElec = (THnSparseF *) dataspectrumbeforesubstraction->GetGrid();  
  incElec = (TH1D *) sparseIncElec->Projection(0);
  CorrectFromTheWidth(incElec);

  TH1D* htemp;
  Int_t* bins=new Int_t[2];

  if(fIPanaHadronBgSubtract){
    // Hadron background
    printf("Hadron background for IP analysis subtracted!\n");
    htemp  = (TH1D *) fHadronEffbyIPcut->Projection(0);
    bins[0]=htemp->GetNbinsX();
    
    AliCFDataGrid *hbgContainer = new AliCFDataGrid("hbgContainer","hadron bg after IP cut",nbins,bins);
    hbgContainer->SetGrid(fHadronEffbyIPcut);
    backgroundGrid->Multiply(hbgContainer,1);
    // subtract hadron contamination
    spectrumSubtracted->Add(backgroundGrid,-1.0);
    // draw raw hadron bg spectra     
    THnSparseF* sparseHbg = (THnSparseF *) hbgContainer->GetGrid();
    hadron = (TH1D *) sparseHbg->Projection(0);
    CorrectFromTheWidth(hadron);
  }

  if(fIPanaCharmBgSubtract){
    // Charm background
    printf("Charm background for IP analysis subtracted!\n");
    AliCFDataGrid *charmbgContainer = (AliCFDataGrid *) GetCharmBackground();
    spectrumSubtracted->Add(charmbgContainer,-1.0);
    THnSparseF *sparseCharmElec = (THnSparseF *) charmbgContainer->GetGrid();  
    charm = (TH1D *) sparseCharmElec->Projection(0);
    CorrectFromTheWidth(charm); 
  }

  const Char_t *sourceName[bgPlots]={"Conversion","Dalitz","K0s secondary pion Dalitz","K0s secondary pion conversions","Other secondary pion Dalitz","Other secondary pion Dalitz conversions","Other Dalitz", "Other conversions"};
  const Int_t style[2] = {21,22};
  const Int_t color[3] = {7,12,8};
  if(fIPanaNonHFEBgSubtract){ 
    Int_t iPlot = 0;
    for(Int_t iSource = 0; iSource < kElecBgSources; iSource++){
      if((iSource>0)&&(iSource<6))continue;//standard sources are summed up in case of iSource == 0; not treated separately in plots
      for(Int_t decay = 0; decay < 2; decay++){
        AliCFDataGrid *nonHFEbgContainer = (AliCFDataGrid *) GetNonHFEBackground(decay,iSource);//conversions/Dalitz decays from different sources
        if(iSource == 0)
          spectrumSubtracted->Add(nonHFEbgContainer,-1.0);
        THnSparseF* sparseNonHFEelecs = (THnSparseF *) nonHFEbgContainer->GetGrid();
        nonHFE[iPlot] = (TH1D *) sparseNonHFEelecs->Projection(0);
        CorrectFromTheWidth(nonHFE[iPlot]);
        iPlot++;
      }
      if(!(fNonHFEsyst && (fBeamType == 1)))break;
    } 
  }

  TLegend *lRaw = new TLegend(0.55,0.55,0.85,0.85);
  TCanvas *cRaw = new TCanvas("cRaw","cRaw",1000,800);     
  cRaw->cd();
  gPad->SetLogx();
  gPad->SetLogy();
  incElec->GetXaxis()->SetRangeUser(0.4,8.);
  incElec->Draw("p");
  incElec->SetMarkerColor(1);
  incElec->SetMarkerStyle(20);
  lRaw->AddEntry(incElec,"inclusive electron spectrum");
  
  if(fIPanaCharmBgSubtract){
    charm->Draw("samep");
    charm->SetMarkerColor(3);
    charm->SetMarkerStyle(20);
    charm->GetXaxis()->SetRangeUser(0.4,8.);
    lRaw->AddEntry(charm,"charm elecs");
  }
  
  if(fIPanaNonHFEBgSubtract){
    Int_t iPlot = 0;
    for(Int_t iSource = 0; iSource < kElecBgSources; iSource++){
      if((iSource>0)&&(iSource<6))continue;//standard sources are summed up in case of iSource == 0; not treated separately in plots
      for(Int_t decay = 0; decay < 2; decay++){
        nonHFE[iPlot]->Draw("samep");
        if(iSource == 0){
          nonHFE[iPlot]->SetMarkerStyle(20);
          if(decay == 0)
            nonHFE[iPlot]->SetMarkerColor(4);
          else
            nonHFE[iPlot]->SetMarkerColor(6);
        }
        else{     
          nonHFE[iPlot]->SetMarkerColor(color[iSource-6]);      
          nonHFE[iPlot]->SetMarkerStyle(style[decay]);
        }
        nonHFE[iPlot]->GetXaxis()->SetRangeUser(0.4,8.);
        lRaw->AddEntry(nonHFE[iPlot],sourceName[iPlot]);
        iPlot++;
      }
      if(!(fNonHFEsyst && (fBeamType == 1)))break;//only draw standard sources
    }
  }

  THnSparseF* sparseSubtracted = (THnSparseF *) spectrumSubtracted->GetGrid();
  subtracted = (TH1D *) sparseSubtracted->Projection(0);
  CorrectFromTheWidth(subtracted);
  subtracted->Draw("samep");
  subtracted->SetMarkerStyle(24);
  lRaw->AddEntry(subtracted,"subtracted electron spectrum");
  lRaw->Draw("SAME");
 
  delete[] bins; 

  TH1D *measuredTH1background = (TH1D *) backgroundGrid->Project(0);
  CorrectFromTheWidth(measuredTH1background);

  if(setBackground){
    if(fBackground) delete fBackground;
    fBackground = backgroundGrid;
  } else delete backgroundGrid;
  
  fQA->AddResultAt(subtracted,AliHFEBeautySpectrumQA::kAfterSC);
  fQA->AddResultAt(incElec,AliHFEBeautySpectrumQA::kBeforeSC);
  fQA->AddResultAt(measuredTH1background, AliHFEBeautySpectrumQA::kMeasBG);
  fQA->DrawSubtractContamination();
  
  return spectrumSubtracted;
}
//____________________________________________________________________________________________
AliCFDataGrid* AliHFEBeautySpectrum::GetNonHFEBackground(Int_t decay, Int_t source){
  //
  // calculate non-HF electron background
  // Arguments: decay = 0 for conversions, decay = 1 for meson decays to electrons
  // source: gives mother either of the electron (for meson->electron) or of the gamma (conversions);
  // source numbers as given by array in AliAnalysisTaskHFE: 
  // {"Pion","Eta","Omega","Phi","EtaPrime","Rho","K0sSec","OtherSecPi0","Else"};
  //
  
  Double_t evtnorm[1] = {1.};
  if(fNMCbgEvents) evtnorm[0]= double(fNEvents)/double(fNMCbgEvents);
  
  Int_t stepbackground = 1;//for non-standard background, step after IP cut is taken directly
  if(source == 0) stepbackground = 3;//for standard sources, take step before PID -> correction of PID x IP efficiencies from enhanced samples
  AliCFContainer *backgroundContainer = 0x0;
  if(fNonHFEsyst){
    if(decay == 0){//conversions
      backgroundContainer = (AliCFContainer*)fConvSourceContainer[source][0]->Clone();
      if(source == 0)//standard conversion containers
        for(Int_t iSource = 1; iSource < kElecBgSources-3; iSource++){//add all standard sources
          backgroundContainer->Add(fConvSourceContainer[iSource][0]);     
        }
    }
    else if(decay == 1){//Dalitz decays, same procedure as above
      backgroundContainer = (AliCFContainer*)fNonHFESourceContainer[source][0]->Clone();
      if(source == 0)
        for(Int_t iSource = 1; iSource < kElecBgSources-3; iSource++){
          backgroundContainer->Add(fNonHFESourceContainer[iSource][0]);   
        }
    }   
  }
  else{//careful: these containers include the non-standard electron sources (K0s, etc.). If necessary, use SetRange in species axis to fix it?
    if(source == 0){
      // Background Estimate
      if(decay == 0)
        backgroundContainer = GetContainer(kMCWeightedContainerConversionESD);
      else if(decay == 1)
        backgroundContainer = GetContainer(kMCWeightedContainerNonHFEESD);
    } 
  }
  if(!backgroundContainer){
    AliError("MC background container not found");
    return NULL;
  }
  AliCFDataGrid *backgroundGrid = new AliCFDataGrid("backgroundGrid","backgroundGrid",*backgroundContainer,stepbackground);
 
  Double_t rangelow = 1.;
  Double_t rangeup = 6.;
  if(decay == 1) rangelow = 0.9;


  const Char_t *dmode[2]={"Conversions","Dalitz decays"};
  TF1 *fitHagedorn = new TF1("fitHagedorn", "[0]/TMath::Power(([1]+x/[2]), [3])", rangelow, rangeup);
  fNonHFEbg = 0x0;
  TH1D *h1 = (TH1D*)backgroundGrid->Project(0);
  TAxis *axis = h1->GetXaxis();
  CorrectFromTheWidth(h1);
  if(source == 0){
    fitHagedorn->SetParameter(0, 0.15);
    fitHagedorn->SetParameter(1, 0.09);
    fitHagedorn->SetParameter(2, 8.4);
    fitHagedorn->SetParameter(3, 6.3);
    // TCanvas *ch1conversion = new TCanvas("ch1conversion","ch1conversion",500,400);
    // ch1conversion->cd();
    //fHagedorn->SetLineColor(2);
    h1->Fit(fitHagedorn,"R");
    // h1->Draw();
    if(!(fNonHFEbg = h1->GetFunction("fitHagedorn")))printf("electron background fit failed for %s\n",dmode[decay]);
  }
  
  Int_t *nBinpp=new Int_t[1];
  Int_t *binspp=new Int_t[1];
  binspp[0]=kSignalPtBins;// number of pt bins
  
  Int_t looppt=binspp[0];
  
  for(Long_t iBin=1; iBin<= looppt;iBin++){       
    Double_t iipt= h1->GetBinCenter(iBin);
    Double_t iiwidth= axis->GetBinWidth(iBin);
    nBinpp[0]=iBin;
    Double_t fitFactor = backgroundGrid->GetElement(nBinpp);//if no fit available, just take bin-by-bin information
    if(fNonHFEbg)fitFactor = fNonHFEbg->Eval(iipt)*iiwidth;
    backgroundGrid->SetElementError(nBinpp, backgroundGrid->GetElementError(nBinpp)*evtnorm[0]);
    backgroundGrid->SetElement(nBinpp,fitFactor*evtnorm[0]);
  }
  //end of workaround for statistical errors
  if(source == 0){  
    AliCFDataGrid *weightedBGContainer = 0x0;   
    weightedBGContainer = new AliCFDataGrid("weightedBGContainer","weightedBGContainer",1,binspp);
    if(decay == 0)
      weightedBGContainer->SetGrid(GetPIDxIPEff(2));
    else if(decay == 1)
      weightedBGContainer->SetGrid(GetPIDxIPEff(3));
    if(stepbackground == 3){    
      backgroundGrid->Multiply(weightedBGContainer,1.0);
    }  
  }
  delete[] nBinpp;
  delete[] binspp;  
  return backgroundGrid;
}
//____________________________________________________________
AliCFDataGrid* AliHFEBeautySpectrum::GetCharmBackground(){
  //
  // calculate charm background; when using centrality binning 2: centBin = 0 for 0-20%, centBin = 5 for 40-80%
  //

  Int_t nDim = 1;
 
  Double_t evtnorm=0;
  if(fNMCEvents) evtnorm= double(fNEvents)/double(fNMCEvents); //mjmjmj
  //if(fNMCbgEvents) evtnorm= double(fNEvents)/double(fNMCbgEvents);
  printf("events for data: %d",fNEvents);
  printf("events for MC: %d",fNMCEvents);
  printf("normalization factor for charm: %f",evtnorm);
  
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
   
  AliCFDataGrid *ipWeightedCharmContainer = new AliCFDataGrid("ipWeightedCharmContainer","ipWeightedCharmContainer",nDim,bins);
  ipWeightedCharmContainer->SetGrid(GetPIDxIPEff(0)); // get charm efficiency
  TH1D* parametrizedcharmpidipeff = (TH1D*)ipWeightedCharmContainer->Project(0);
  
  charmBackgroundGrid->Multiply(ipWeightedCharmContainer,1.);
  Int_t *nBinpp=new Int_t[1];
  Int_t* binspp=new Int_t[1];
  binspp[0]=charmbgaftertofpid->GetNbinsX();  // number of pt bins
  
  Int_t looppt=binspp[0];
  
  for(Long_t iBin=1; iBin<= looppt;iBin++){
    nBinpp[0]=iBin;
    charmBackgroundGrid->SetElementError(nBinpp, charmBackgroundGrid->GetElementError(nBinpp)*evtnorm);
    charmBackgroundGrid->SetElement(nBinpp,charmBackgroundGrid->GetElement(nBinpp)*evtnorm);
  }
  
  TH1D* charmbgafteripcut = (TH1D*)charmBackgroundGrid->Project(0);
  
  AliCFDataGrid *weightedCharmContainer = new AliCFDataGrid("weightedCharmContainer","weightedCharmContainer",nDim,bins);
  weightedCharmContainer->SetGrid(GetCharmWeights(fTestCentralityLow)); // get charm weighting factors
  TH1D* charmweightingfc = (TH1D*)weightedCharmContainer->Project(0);
  charmBackgroundGrid->Multiply(weightedCharmContainer,1.);
  TH1D* charmbgafterweight = (TH1D*)charmBackgroundGrid->Project(0);
  
  // Efficiency (set efficiency to 1 for only folding) 
  AliCFEffGrid* efficiencyD = new AliCFEffGrid("efficiency","",*mcContainer);
  efficiencyD->CalculateEfficiency(0,0);
  
  // Folding
  if(fBeamType==0)nDim = 1;
  AliCFUnfolding folding("unfolding","",nDim,fCorrelation,efficiencyD->GetGrid(),charmBackgroundGrid->GetGrid(),charmBackgroundGrid->GetGrid());
  folding.SetMaxNumberOfIterations(1);
  folding.Unfold();
  
  // Results
  THnSparse* result1= folding.GetEstMeasured(); // folded spectra
  THnSparse* result=(THnSparse*)result1->Clone();
  charmBackgroundGrid->SetGrid(result);
  TH1D* charmbgafterfolding = (TH1D*)charmBackgroundGrid->Project(0);

  //Charm background evaluation plots
  
  TCanvas *cCharmBgEval = new TCanvas("cCharmBgEval","cCharmBgEval",1000,600);
  cCharmBgEval->Divide(3,1);
  
  cCharmBgEval->cd(1);
  
  if(fBeamType==0)charmbgaftertofpid->Scale(evtnorm);
    
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
  
  if(fUnfoldBG) UnfoldBG(charmBackgroundGrid);

  return charmBackgroundGrid;
}
//____________________________________________________________
AliCFDataGrid *AliHFEBeautySpectrum::CorrectParametrizedEfficiency(AliCFDataGrid* const bgsubpectrum){  
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
  fQA->AddResultAt(afterE,AliHFEBeautySpectrumQA::kAfterPE);
  fQA->AddResultAt(beforeE,AliHFEBeautySpectrumQA::kBeforePE);
  fQA->AddResultAt(fEfficiencyFunction,AliHFEBeautySpectrumQA::kPEfficiency);
  fQA->DrawCorrectWithEfficiency(AliHFEBeautySpectrumQA::kParametrized);
  
  return resultt;
}
//____________________________________________________________
THnSparse *AliHFEBeautySpectrum::Unfold(AliCFDataGrid* const bgsubpectrum){
  
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
    Int_t* bins=new Int_t[1];
    bins[0]=kSignalPtBins;
    
    AliCFEffGrid *beffContainer = new AliCFEffGrid("beffContainer","beffContainer",1,bins);
    beffContainer->SetGrid(GetBeautyIPEff(kTRUE));
    efficiencyD->Multiply(beffContainer,1);
  }
  

  // Unfold 
  
  AliCFUnfolding unfolding("unfolding","",fNbDimensions,fCorrelation,efficiencyD->GetGrid(),dataGrid->GetGrid(),guessedTHnSparse,1.e-06,0,fNumberOfIterations);//check with MinJung if last arguments are correct (taken from inclusive analysis...
  if(fUnSetCorrelatedErrors) unfolding.UnsetCorrelatedErrors();
  unfolding.SetMaxNumberOfIterations(fNumberOfIterations);
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
  fQA->AddResultAt(residualh,AliHFEBeautySpectrumQA::kResidualU);
  fQA->AddResultAt(afterE,AliHFEBeautySpectrumQA::kAfterU);
  fQA->AddResultAt(beforeE,AliHFEBeautySpectrumQA::kBeforeU);
  fQA->AddResultAt(efficiencyDproj,AliHFEBeautySpectrumQA::kUEfficiency);
  fQA->DrawUnfolding();

  return (THnSparse *) result->Clone();
}

//____________________________________________________________
AliCFDataGrid *AliHFEBeautySpectrum::CorrectForEfficiency(AliCFDataGrid* const bgsubpectrum){
  
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
    Int_t* bins=new Int_t[1];
    bins[0]=kSignalPtBins;
    
    AliCFEffGrid *beffContainer = new AliCFEffGrid("beffContainer","beffContainer",1,bins);
    beffContainer->SetGrid(GetBeautyIPEff(kFALSE));
    efficiencyD->Multiply(beffContainer,1);
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

  // QA
  TH1D *afterE = (TH1D *) result->Project(0);
  CorrectFromTheWidth(afterE);
  TH1D *beforeE = (TH1D *) dataGrid->Project(0);
  CorrectFromTheWidth(beforeE);
  TH1D* efficiencyDproj = (TH1D *) efficiencyD->Project(0);
  fQA->AddResultAt(afterE,AliHFEBeautySpectrumQA::kAfterMCE);
  fQA->AddResultAt(beforeE,AliHFEBeautySpectrumQA::kBeforeMCE);
  fQA->AddResultAt(efficiencyDproj,AliHFEBeautySpectrumQA::kMCEfficiency);
  fQA->DrawCorrectWithEfficiency(AliHFEBeautySpectrumQA::kMC);
  
  return result;

}
//____________________________________________________________
void AliHFEBeautySpectrum::AddTemporaryObject(TObject *o){
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
void AliHFEBeautySpectrum::ClearObject(TObject *o){
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
THnSparse* AliHFEBeautySpectrum::GetCharmWeights(Int_t centBin){
 
  //
  // Measured D->e based weighting factors
  //

  const Int_t nDimpp=1;
  Int_t nBinpp[nDimpp] = {kSignalPtBins};
  Double_t ptbinning1[kSignalPtBins+1] = {0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.75, 2., 2.25, 2.5, 2.75, 3., 3.5, 4., 4.5, 5., 5.5, 6., 7., 8., 10., 12., 14., 16., 18., 20.};
  Int_t looppt=nBinpp[0];
  
  fWeightCharm = new THnSparseF("weightHisto", "weighting factor; pt[GeV/c]", nDimpp, nBinpp);
  fWeightCharm->SetBinEdges(0, ptbinning1);
 
  // //if(fBeamType == 0){// Weighting factor for pp
  //Double_t weightpp[kSignalPtBins]={0.859260, 0.872552, 0.847475, 0.823631, 0.839386, 0.874024, 0.916755, 0.942801, 0.965856, 0.933905, 0.933414, 0.931936, 0.847826, 0.810902, 0.796608, 0.727002, 0.659227, 0.583610, 0.549956, 0.512633, 0.472254, 0.412364, 0.353191, 0.319145, 0.305412, 0.290334, 0.269863, 0.254646, 0.230245, 0.200859, 0.275953, 0.276271, 0.227332, 0.197004, 0.474385};
//TPC+TOF standard cut 4800
  Double_t weightpp[kSignalPtBins]={0.050326, 0.045826, 0.042043, 0.039641, 0.039872, 0.041796, 0.044320, 0.046273, 0.047376, 0.047657, 0.047973, 0.048307, 0.045325, 0.044067, 0.043367, 0.040417, 0.037048, 0.033695, 0.032192, 0.029270, 0.027270, 0.024451, 0.020846, 0.019032, 0.018210, 0.017554, 0.015604, 0.015194, 0.013542, 0.013447, 0.015160, 0.015507, 0.014989, 0.012533, 0.025603}; 
  //}
  //else{
  //if(centBin == 0){
      // Weighting factor for PbPb (0-20%)
  Double_t weightPbPb1[kSignalPtBins]={0.641897,  0.640472,  0.615228,  0.650469,  0.737762,  0.847867,  1.009317,  1.158594,  1.307482,  1.476973,  1.551131,  1.677131,  1.785478,  1.888933,  2.017957,  2.074757,  1.926700,  1.869495,  1.546558,  1.222873,  1.160313,  0.903375,  0.799642,  0.706244,  0.705449,  0.599947,  0.719570,  0.499422,  0.703978,  0.477452,  0.325057,  0.093391,  0.096675,  0.000000,  0.000000};
  //}
  //else{
  // Weighting factor for PbPb (40-80%)
  Double_t weightPbPb2[kSignalPtBins]={0.181953,  0.173248,  0.166799,  0.182558,  0.206581,  0.236955,  0.279390,  0.329129,  0.365260,  0.423059,  0.452057,  0.482726,  0.462627,  0.537770,  0.584663,  0.579452,  0.587194,  0.499498,  0.443299,  0.398596,  0.376695,  0.322331,  0.260890,  0.374834,  0.249114,  0.310330,  0.287326,  0.243174,  0.758945,  0.138867,  0.170576,  0.107797,  0.011390,  0.000000,  0.000000};
  //  }
  //}
  Double_t weight[kSignalPtBins];
  for(Int_t ipt = 0; ipt < kSignalPtBins; ipt++){
    if(fBeamType == 0)weight[ipt] = weightpp[ipt];
    else if(centBin == 0)weight[ipt] = weightPbPb1[ipt];
    else weight[ipt] = weightPbPb2[ipt];
  }  

  //points
  Double_t pt[1];
  Double_t contents[2];
  
  for(int i=0; i<looppt; i++){
    pt[0]=(ptbinning1[i]+ptbinning1[i+1])/2.;
    contents[0]=pt[0];
    fWeightCharm->Fill(contents,weight[i]);
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
void AliHFEBeautySpectrum::SetParameterizedEff(AliCFContainer *container, AliCFContainer *containermb, AliCFContainer *containeresd, AliCFContainer *containeresdmb, Int_t *dimensions){

   // TOF PID efficiencies
  
   TF1 *fittofpid = new TF1("fittofpid","[0]*([1]+[2]*log(x)+[3]*log(x)*log(x))*tanh([4]*x-[5])",0.95,8.);
   TF1 *fipfit = new TF1("fipfit","[0]*([1]+[2]*log(x)+[3]*log(x)*log(x))*tanh([4]*x-[5])",0.95,6.);
   TF1 *fipfitnonhfe = new TF1("fipfitnonhfe","[0]*([1]+[2]*log(x)+[3]*log(x)*log(x))*tanh([4]*x-[5])",0.5,8.0);

   if(fBeamType == 1){//not in use yet - adapt as soon as possible!
     fittofpid = new TF1("fittofpid","[0]*([1]+[2]*log(x)+[3]*log(x)*log(x))*tanh([4]*x-[5])",0.95,8.);
     fipfit = new TF1("fipfit","[0]*([1]+[2]*log(x)+[3]*log(x)*log(x))*tanh([4]*x-[5])",0.95,8.);
     fipfitnonhfe = new TF1("fipfitnonhfe","[0]*([1]+[2]*log(x)+[3]*log(x)*log(x))*tanh([4]*x-[5])",0.5,8.);
   }

   TCanvas * cefficiencyParamtof = new TCanvas("efficiencyParamtof","efficiencyParamtof",600,600);
   cefficiencyParamtof->cd();

   AliCFContainer *mccontainermcD = 0x0;
   AliCFContainer *mccontaineresdD = 0x0;
   TH1D* efficiencysigTOFPIDD;
   TH1D* efficiencyTOFPIDD;
   TH1D* efficiencysigesdTOFPIDD;
   TH1D* efficiencyesdTOFPIDD;
   Int_t source = -1; //get parameterized TOF PID efficiencies

   // signal sample
   mccontainermcD = GetSlicedContainer(container, fNbDimensions, dimensions, source, fChargeChoosen,fTestCentralityLow,fTestCentralityHigh);
   AliCFEffGrid* efficiencymcsigParamTOFPID= new AliCFEffGrid("efficiencymcsigParamTOFPID","",*mccontainermcD);
   efficiencymcsigParamTOFPID->CalculateEfficiency(fStepMC,fStepMC-1); // TOF PID efficiencies
   
   // mb sample for double check
   if(fBeamType==0)  mccontainermcD = GetSlicedContainer(containermb, fNbDimensions, dimensions, source, fChargeChoosen,fTestCentralityLow,fTestCentralityHigh);
   AliCFEffGrid* efficiencymcParamTOFPID= new AliCFEffGrid("efficiencymcParamTOFPID","",*mccontainermcD);
   efficiencymcParamTOFPID->CalculateEfficiency(fStepMC,fStepMC-1); // TOF PID efficiencies
   
   // mb sample with reconstructed variables
   if(fBeamType==0)  mccontainermcD = GetSlicedContainer(containeresdmb, fNbDimensions, dimensions, source, fChargeChoosen,fTestCentralityLow,fTestCentralityHigh);
   AliCFEffGrid* efficiencyesdParamTOFPID= new AliCFEffGrid("efficiencyesdParamTOFPID","",*mccontainermcD);
   efficiencyesdParamTOFPID->CalculateEfficiency(fStepMC,fStepMC-1); // TOF PID efficiencies
   
   // mb sample with reconstructed variables
   if(fBeamType==0)  mccontainermcD = GetSlicedContainer(containeresd, fNbDimensions, dimensions, source, fChargeChoosen,fTestCentralityLow,fTestCentralityHigh);
   AliCFEffGrid* efficiencysigesdParamTOFPID= new AliCFEffGrid("efficiencysigesdParamTOFPID","",*mccontainermcD);
   efficiencysigesdParamTOFPID->CalculateEfficiency(fStepMC,fStepMC-1); // TOF PID efficiencies
   
   //fill histo
   efficiencysigTOFPIDD = (TH1D *) efficiencymcsigParamTOFPID->Project(0);
   efficiencyTOFPIDD = (TH1D *) efficiencymcParamTOFPID->Project(0);
   efficiencysigesdTOFPIDD = (TH1D *) efficiencysigesdParamTOFPID->Project(0);
   efficiencyesdTOFPIDD = (TH1D *) efficiencyesdParamTOFPID->Project(0);
   efficiencysigTOFPIDD->SetName("efficiencysigTOFPIDD");
   efficiencyTOFPIDD->SetName("efficiencyTOFPIDD");
   efficiencysigesdTOFPIDD->SetName("efficiencysigesdTOFPIDD");
   efficiencyesdTOFPIDD->SetName("efficiencyesdTOFPIDD");
   
   //fit (mc enhenced sample)
   fittofpid->SetParameters(0.5,0.319,0.0157,0.00664,6.77,2.08);
   efficiencysigTOFPIDD->Fit(fittofpid,"R");
   efficiencysigTOFPIDD->GetYaxis()->SetTitle("Efficiency");
   fEfficiencyTOFPIDD = efficiencysigTOFPIDD->GetFunction("fittofpid");
   
   //fit (esd enhenced sample)
   efficiencysigesdTOFPIDD->Fit(fittofpid,"R");
   efficiencysigesdTOFPIDD->GetYaxis()->SetTitle("Efficiency");
   fEfficiencyesdTOFPIDD = efficiencysigesdTOFPIDD->GetFunction("fittofpid");
   
   // draw (for PbPb, only 1st bin)
   //sig mc
   efficiencysigTOFPIDD->SetTitle("");
   efficiencysigTOFPIDD->SetStats(0);
   efficiencysigTOFPIDD->SetMarkerStyle(25);
   efficiencysigTOFPIDD->SetMarkerColor(2);
   efficiencysigTOFPIDD->SetLineColor(2);
   efficiencysigTOFPIDD->Draw();

   //mb mc
   efficiencyTOFPIDD->SetTitle("");
   efficiencyTOFPIDD->SetStats(0);
   efficiencyTOFPIDD->SetMarkerStyle(24);
   efficiencyTOFPIDD->SetMarkerColor(4);
   efficiencyTOFPIDD->SetLineColor(4);
   efficiencyTOFPIDD->Draw("same");

   //sig esd
   efficiencysigesdTOFPIDD->SetTitle("");
   efficiencysigesdTOFPIDD->SetStats(0);
   efficiencysigesdTOFPIDD->SetMarkerStyle(25);
   efficiencysigesdTOFPIDD->SetMarkerColor(3);
   efficiencysigesdTOFPIDD->SetLineColor(3);
   efficiencysigesdTOFPIDD->Draw("same");

   //mb esd
   efficiencyesdTOFPIDD->SetTitle("");
   efficiencyesdTOFPIDD->SetStats(0);
   efficiencyesdTOFPIDD->SetMarkerStyle(25);
   efficiencyesdTOFPIDD->SetMarkerColor(1);
   efficiencyesdTOFPIDD->SetLineColor(1);
   efficiencyesdTOFPIDD->Draw("same");

   //signal mc fit
   if(fEfficiencyTOFPIDD){
     fEfficiencyTOFPIDD->SetLineColor(2);
     fEfficiencyTOFPIDD->Draw("same");
   }
   //mb esd fit
   if(fEfficiencyesdTOFPIDD){
       fEfficiencyesdTOFPIDD->SetLineColor(3);
       fEfficiencyesdTOFPIDD->Draw("same");
     }

   TLegend *legtofeff = new TLegend(0.3,0.15,0.79,0.44);
   legtofeff->AddEntry(efficiencysigTOFPIDD,"TOF PID Step Efficiency","");
   legtofeff->AddEntry(efficiencysigTOFPIDD,"vs MC p_{t} for enhenced samples","p");
   legtofeff->AddEntry(efficiencyTOFPIDD,"vs MC p_{t} for mb samples","p");
   legtofeff->AddEntry(efficiencysigesdTOFPIDD,"vs esd p_{t} for enhenced samples","p");
   legtofeff->AddEntry(efficiencyesdTOFPIDD,"vs esd p_{t} for mb samples","p");
   legtofeff->Draw("same");


   TCanvas * cefficiencyParamIP = new TCanvas("efficiencyParamIP","efficiencyParamIP",500,500);
   cefficiencyParamIP->cd();
   gStyle->SetOptStat(0);

   // IP cut efficiencies
   AliCFContainer *charmCombined = 0x0; 
   AliCFContainer *beautyCombined = 0x0;
   AliCFContainer *beautyCombinedesd = 0x0;

   source = 0; //charm enhenced
   mccontainermcD = GetSlicedContainer(container, fNbDimensions, dimensions, source, fChargeChoosen,fTestCentralityLow,fTestCentralityHigh);
   AliCFEffGrid* efficiencyCharmSig = new AliCFEffGrid("efficiencyCharmSig","",*mccontainermcD);
   efficiencyCharmSig->CalculateEfficiency(AliHFEcuts::kNcutStepsMCTrack + fStepData,AliHFEcuts::kNcutStepsMCTrack + fStepData-1); // ip cut efficiency. 

   charmCombined= (AliCFContainer*)mccontainermcD->Clone("charmCombined");  

   source = 1; //beauty enhenced
   mccontainermcD = GetSlicedContainer(container, fNbDimensions, dimensions, source, fChargeChoosen,fTestCentralityLow,fTestCentralityHigh);
   AliCFEffGrid* efficiencyBeautySig = new AliCFEffGrid("efficiencyBeautySig","",*mccontainermcD);
   efficiencyBeautySig->CalculateEfficiency(AliHFEcuts::kNcutStepsMCTrack + fStepData,AliHFEcuts::kNcutStepsMCTrack + fStepData-1); // ip cut efficiency. 

   beautyCombined = (AliCFContainer*)mccontainermcD->Clone("beautyCombined"); 

   mccontaineresdD = GetSlicedContainer(containeresd, fNbDimensions, dimensions, source, fChargeChoosen,fTestCentralityLow,fTestCentralityHigh);
   AliCFEffGrid* efficiencyBeautySigesd = new AliCFEffGrid("efficiencyBeautySigesd","",*mccontaineresdD);
   efficiencyBeautySigesd->CalculateEfficiency(AliHFEcuts::kNcutStepsMCTrack + fStepData,AliHFEcuts::kNcutStepsMCTrack + fStepData-1); // ip cut efficiency.

   beautyCombinedesd = (AliCFContainer*)mccontaineresdD->Clone("beautyCombinedesd");

   source = 0; //charm mb
   mccontainermcD = GetSlicedContainer(containermb, fNbDimensions, dimensions, source, fChargeChoosen,fTestCentralityLow,fTestCentralityHigh);
   AliCFEffGrid* efficiencyCharm = new AliCFEffGrid("efficiencyCharm","",*mccontainermcD);
   efficiencyCharm->CalculateEfficiency(AliHFEcuts::kNcutStepsMCTrack + fStepData,AliHFEcuts::kNcutStepsMCTrack + fStepData-1); // ip cut efficiency. 
   
   charmCombined->Add(mccontainermcD); 
   AliCFEffGrid* efficiencyCharmCombined = new AliCFEffGrid("efficiencyCharmCombined","",*charmCombined); 
   efficiencyCharmCombined->CalculateEfficiency(AliHFEcuts::kNcutStepsMCTrack + fStepData,AliHFEcuts::kNcutStepsMCTrack + fStepData-1); 

   source = 1; //beauty mb
   mccontainermcD = GetSlicedContainer(containermb, fNbDimensions, dimensions, source, fChargeChoosen,fTestCentralityLow,fTestCentralityHigh);
   AliCFEffGrid* efficiencyBeauty = new AliCFEffGrid("efficiencyBeauty","",*mccontainermcD);
   efficiencyBeauty->CalculateEfficiency(AliHFEcuts::kNcutStepsMCTrack + fStepData,AliHFEcuts::kNcutStepsMCTrack + fStepData-1); // ip cut efficiency. 
   
   beautyCombined->Add(mccontainermcD);
   AliCFEffGrid* efficiencyBeautyCombined = new AliCFEffGrid("efficiencyBeautyCombined","",*beautyCombined); 
   efficiencyBeautyCombined->CalculateEfficiency(AliHFEcuts::kNcutStepsMCTrack + fStepData,AliHFEcuts::kNcutStepsMCTrack + fStepData-1); 
   
   mccontaineresdD = GetSlicedContainer(containeresdmb, fNbDimensions, dimensions, source, fChargeChoosen,fTestCentralityLow,fTestCentralityHigh);
   AliCFEffGrid* efficiencyBeautyesd = new AliCFEffGrid("efficiencyBeautyesd","",*mccontaineresdD);
   efficiencyBeautyesd->CalculateEfficiency(AliHFEcuts::kNcutStepsMCTrack + fStepData,AliHFEcuts::kNcutStepsMCTrack + fStepData-1); // ip cut efficiency.

   beautyCombinedesd->Add(mccontaineresdD);
   AliCFEffGrid* efficiencyBeautyCombinedesd = new AliCFEffGrid("efficiencyBeautyCombinedesd","",*beautyCombinedesd);
   efficiencyBeautyCombinedesd->CalculateEfficiency(AliHFEcuts::kNcutStepsMCTrack + fStepData,AliHFEcuts::kNcutStepsMCTrack + fStepData-1);
   
   source = 2; //conversion mb
   mccontainermcD = GetSlicedContainer(containermb, fNbDimensions, dimensions, source, fChargeChoosen,fTestCentralityLow,fTestCentralityHigh);
   AliCFEffGrid* efficiencyConv = new AliCFEffGrid("efficiencyConv","",*mccontainermcD);
   efficiencyConv->CalculateEfficiency(AliHFEcuts::kNcutStepsMCTrack + fStepData,AliHFEcuts::kNcutStepsMCTrack + fStepData-1); // ip cut efficiency. 

   source = 3; //non HFE except for the conversion mb
   mccontainermcD = GetSlicedContainer(containermb, fNbDimensions, dimensions, source, fChargeChoosen,fTestCentralityLow,fTestCentralityHigh);
   AliCFEffGrid* efficiencyNonhfe= new AliCFEffGrid("efficiencyNonhfe","",*mccontainermcD);
   efficiencyNonhfe->CalculateEfficiency(AliHFEcuts::kNcutStepsMCTrack + fStepData,AliHFEcuts::kNcutStepsMCTrack + fStepData-1); // ip cut efficiency.

   if(fIPEffCombinedSamples){printf("combined samples taken for beauty and charm\n");
     fEfficiencyCharmSigD = (TH1D*)efficiencyCharmCombined->Project(0); //signal enhenced + mb 
     fEfficiencyBeautySigD = (TH1D*)efficiencyBeautyCombined->Project(0); //signal enhenced + mb
     fEfficiencyBeautySigesdD = (TH1D*)efficiencyBeautyCombinedesd->Project(0); //signal enhenced + mb
     }
     else{printf("signal enhanced samples taken for beauty and charm\n");
       fEfficiencyCharmSigD = (TH1D*)efficiencyCharmSig->Project(0); //signal enhenced only
       fEfficiencyBeautySigD = (TH1D*)efficiencyBeautySig->Project(0); //signal enhenced only
       fEfficiencyBeautySigesdD = (TH1D*)efficiencyBeautySigesd->Project(0); //signal enhenced only
     }
     fCharmEff = (TH1D*)efficiencyCharm->Project(0); //mb only
     fBeautyEff = (TH1D*)efficiencyBeauty->Project(0); //mb only
     fConversionEff = (TH1D*)efficiencyConv->Project(0); //mb only
     fNonHFEEff = (TH1D*)efficiencyNonhfe->Project(0); //mb only

   if(fBeamType==0){
     //AliCFEffGrid  *nonHFEEffGrid = (AliCFEffGrid*)  GetEfficiency(GetContainer(kMCWeightedContainerNonHFEESD),1,0);
     AliCFContainer *cfcontainer = GetContainer(AliHFECorrectSpectrumBase::kMCWeightedContainerNonHFEESDSig);
     if(!cfcontainer) return;
     AliCFEffGrid  *nonHFEEffGrid = (AliCFEffGrid*)  GetEfficiency(cfcontainer,1,0); //mjmj
     fNonHFEEffbgc = (TH1D *) nonHFEEffGrid->Project(0);
     
     //AliCFEffGrid  *conversionEffGrid = (AliCFEffGrid*)  GetEfficiency(GetContainer(kMCWeightedContainerConversionESD),1,0);
     AliCFContainer *cfcontainerr = GetContainer(AliHFECorrectSpectrumBase::kMCWeightedContainerConversionESDSig);
     if(!cfcontainerr) return;
     AliCFEffGrid  *conversionEffGrid = (AliCFEffGrid*)  GetEfficiency(cfcontainerr,1,0); //mjmj
     fConversionEffbgc = (TH1D *) conversionEffGrid->Project(0);
     
     //MHMH
     if(fBeamType == 0){
       fipfitnonhfe->SetParameters(0.5,0.319,0.0157,0.00664,6.77,2.08);
       fNonHFEEffbgc->Fit(fipfitnonhfe,"R"); 
       fEfficiencyIPNonhfeD = fNonHFEEffbgc->GetFunction("fipfitnonhfe");
       
       fipfitnonhfe = new TF1("fipfitnonhfe","[0]*([1]+[2]*log(x)+[3]*log(x)*log(x))*tanh([4]*x-[5])",0.5,8.);
       fipfitnonhfe->SetParameters(0.5,0.319,0.0157,0.00664,6.77,2.08);      
       fConversionEffbgc->Fit(fipfitnonhfe,"R");
       fEfficiencyIPConversionD = fConversionEffbgc->GetFunction("fipfitnonhfe");
     }
     //MHMH
   }
   
   fipfit->SetParameters(0.5,0.319,0.0157,0.00664,6.77,2.08);
   fipfit->SetLineColor(2);
   fEfficiencyBeautySigD->Fit(fipfit,"R");
   fEfficiencyBeautySigD->GetYaxis()->SetTitle("Efficiency");
   fEfficiencyIPBeautyD = fEfficiencyBeautySigD->GetFunction("fipfit");
   
   fipfit->SetParameters(0.5,0.319,0.0157,0.00664,6.77,2.08);
   fipfit->SetLineColor(6);
   fEfficiencyBeautySigesdD->Fit(fipfit,"R");
   fEfficiencyBeautySigesdD->GetYaxis()->SetTitle("Efficiency");
   fEfficiencyIPBeautyesdD = fEfficiencyBeautySigesdD->GetFunction("fipfit");
   
   fipfit->SetParameters(0.5,0.319,0.0157,0.00664,6.77,2.08);
   fipfit->SetLineColor(1);
   fEfficiencyCharmSigD->Fit(fipfit,"R");
   fEfficiencyCharmSigD->GetYaxis()->SetTitle("Efficiency");
   fEfficiencyIPCharmD = fEfficiencyCharmSigD->GetFunction("fipfit");
   
   //if(0){
   if(fIPParameterizedEff){
     // fipfitnonhfe->SetParameters(0.5,0.319,0.0157,0.00664,6.77,2.08);
     fipfitnonhfe->SetLineColor(3);
     if(fBeamType==1){
       fConversionEff->Fit(fipfitnonhfe,"R");
       fConversionEff->GetYaxis()->SetTitle("Efficiency");
       fEfficiencyIPConversionD = fConversionEff->GetFunction("fipfitnonhfe");
     }
     else{
       fConversionEffbgc->Fit(fipfitnonhfe,"R");
       fConversionEffbgc->GetYaxis()->SetTitle("Efficiency");
       fEfficiencyIPConversionD = fConversionEffbgc->GetFunction("fipfitnonhfe");
     }       
     // fipfitnonhfe->SetParameters(0.5,0.319,0.0157,0.00664,6.77,2.08);
     fipfitnonhfe->SetLineColor(4);
     if(fBeamType==1){
       fNonHFEEff->Fit(fipfitnonhfe,"R");
       fNonHFEEff->GetYaxis()->SetTitle("Efficiency");
       fEfficiencyIPNonhfeD = fNonHFEEff->GetFunction("fipfitnonhfe");
     }
     else{
       fNonHFEEffbgc->Fit(fipfitnonhfe,"R");
       fNonHFEEffbgc->GetYaxis()->SetTitle("Efficiency");
       fEfficiencyIPNonhfeD = fNonHFEEffbgc->GetFunction("fipfitnonhfe");
     }
   }
   
   // draw
   fEfficiencyCharmSigD->SetMarkerStyle(21);
   fEfficiencyCharmSigD->SetMarkerColor(1);
   fEfficiencyCharmSigD->SetLineColor(1);
   fEfficiencyBeautySigD->SetMarkerStyle(21);
   fEfficiencyBeautySigD->SetMarkerColor(2);
   fEfficiencyBeautySigD->SetLineColor(2);
   fEfficiencyBeautySigesdD->SetStats(0);
   fEfficiencyBeautySigesdD->SetMarkerStyle(21);
   fEfficiencyBeautySigesdD->SetMarkerColor(6);
   fEfficiencyBeautySigesdD->SetLineColor(6);
   fCharmEff->SetMarkerStyle(24);
   fCharmEff->SetMarkerColor(1);
   fCharmEff->SetLineColor(1);
   fBeautyEff->SetMarkerStyle(24);
   fBeautyEff->SetMarkerColor(2);
   fBeautyEff->SetLineColor(2);
   fConversionEff->SetMarkerStyle(24);
   fConversionEff->SetMarkerColor(3);
   fConversionEff->SetLineColor(3);
   fNonHFEEff->SetMarkerStyle(24);
   fNonHFEEff->SetMarkerColor(4);
   fNonHFEEff->SetLineColor(4);

   fEfficiencyCharmSigD->Draw();
   fEfficiencyCharmSigD->GetXaxis()->SetRangeUser(0.0,7.9);
   fEfficiencyCharmSigD->GetYaxis()->SetRangeUser(0.0,0.5);

   fEfficiencyBeautySigD->Draw("same");
   fEfficiencyBeautySigesdD->Draw("same");
   if(fBeamType == 1){
     fNonHFEEff->Draw("same");
     fConversionEff->Draw("same");
   }
   //fCharmEff->Draw("same");
   //fBeautyEff->Draw("same");

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
 
   if(fEfficiencyIPBeautyD)
      fEfficiencyIPBeautyD->Draw("same");
   if(fEfficiencyIPBeautyesdD)
     fEfficiencyIPBeautyesdD->Draw("same");
   if( fEfficiencyIPCharmD)
     fEfficiencyIPCharmD->Draw("same");
   if(fIPParameterizedEff){
     if(fEfficiencyIPConversionD)
       fEfficiencyIPConversionD->Draw("same");
     if(fEfficiencyIPNonhfeD)
       fEfficiencyIPNonhfeD->Draw("same");
   }
   TLegend *legipeff = new TLegend(0.58,0.2,0.88,0.39);
   legipeff->AddEntry(fEfficiencyBeautySigD,"IP Step Efficiency","");
   legipeff->AddEntry(fEfficiencyBeautySigD,"beauty e","p");
   legipeff->AddEntry(fEfficiencyBeautySigesdD,"beauty e(esd pt)","p");
   legipeff->AddEntry(fEfficiencyCharmSigD,"charm e","p");
   if(fBeamType == 0){
     legipeff->AddEntry(fConversionEffbgc,"conversion e(esd pt)","p");
     legipeff->AddEntry(fNonHFEEffbgc,"Dalitz e(esd pt)","p");
   }
   else{
     legipeff->AddEntry(fConversionEff,"conversion e","p");
     legipeff->AddEntry(fNonHFEEff,"Dalitz e","p");
   }
   legipeff->Draw("same");
   gPad->SetGrid();
   //cefficiencyParamIP->SaveAs("efficiencyParamIP.eps");
}

//____________________________________________________________________________
THnSparse* AliHFEBeautySpectrum::GetBeautyIPEff(Bool_t isMCpt){
  //
  // Return beauty electron IP cut efficiency
  //

  const Int_t nPtbinning1 = kSignalPtBins;//number of pt bins, according to new binning
  Double_t kPtRange[nPtbinning1+1] = { 0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.75, 2., 2.25, 2.5, 2.75, 3., 3.5, 4., 4.5, 5., 5.5, 6., 7., 8., 10., 12., 14., 16., 18., 20.};//pt bin limits
 
  Int_t nDim=1;  //dimensions of the efficiency weighting grid
  
  Int_t nBin[1] = {nPtbinning1};
 
  THnSparseF *ipcut = new THnSparseF("beff", "b IP efficiency; p_{t}(GeV/c)", nDim, nBin);
  
  ipcut->SetBinEdges(0, kPtRange);
  
  Double_t pt[1];
  Double_t weight;
  Double_t weightErr;
  Double_t contents[2];

  weight = 1.0;
  weightErr = 1.0;
  
  Int_t looppt=nBin[0];
 
  Int_t ibin[2];
  
  for(int i=0; i<looppt; i++)
    {
      pt[0]=(kPtRange[i]+kPtRange[i+1])/2.;
      if(isMCpt){
        if(fEfficiencyIPBeautyD){
          weight=fEfficiencyIPBeautyD->Eval(pt[0]);
          weightErr = 0;
        }
        else{
          printf("Fit failed on beauty IP cut efficiency. Contents in histo used!\n");
          weight = fEfficiencyBeautySigD->GetBinContent(i+1); 
          weightErr = fEfficiencyBeautySigD->GetBinError(i+1);
        }
      }
      else{
        if(fEfficiencyIPBeautyesdD){
          weight=fEfficiencyIPBeautyesdD->Eval(pt[0]);
          weightErr = 0;
        }
        else{
          printf("Fit failed on beauty IP cut efficiency. Contents in histo used!\n");
          weight = fEfficiencyBeautySigesdD->GetBinContent(i+1);
          weightErr = fEfficiencyBeautySigD->GetBinError(i+1);
        }
      }
      
      contents[0]=pt[0];
      ibin[0]=i+1;
      
      ipcut->Fill(contents,weight);
      ipcut->SetBinError(ibin,weightErr);
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
THnSparse* AliHFEBeautySpectrum::GetPIDxIPEff(Int_t source){
  //
  // Return PID x IP cut efficiency
  //
  const Int_t nPtbinning1 = kSignalPtBins;//number of pt bins, according to new binning
  Double_t kPtRange[nPtbinning1+1] = { 0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.75, 2., 2.25, 2.5, 2.75, 3., 3.5, 4., 4.5, 5., 5.5, 6., 7., 8., 10., 12., 14., 16., 18., 20.};//pt bin limits
  Int_t nDim=1;  //dimensions of the efficiency weighting grid
  
  Int_t nBin[1] = {nPtbinning1};
  
  THnSparseF *pideff;
  pideff = new THnSparseF("pideff", "PID efficiency; p_{t}(GeV/c)", nDim, nBin);
  pideff->SetBinEdges(0, kPtRange);
  
  Double_t pt[1];
  Double_t weight;
  Double_t weightErr;
  Double_t contents[2];
  
  weight = 1.0;
  weightErr = 1.0;
  
  Int_t looppt=nBin[0];
  Int_t ibin[2];
      
  Double_t trdtpcPidEfficiency = fEfficiencyFunction->Eval(0); // assume we have constant TRD+TPC PID efficiency
  for(int i=0; i<looppt; i++)
    {
      pt[0]=(kPtRange[i]+kPtRange[i+1])/2.;
      
      Double_t tofpideff = 1.;
      Double_t tofpideffesd = 1.;
      if(fEfficiencyTOFPIDD)
        tofpideff = fEfficiencyTOFPIDD->Eval(pt[0]); 
      else{
        printf("TOF PID fit failed on conversion. The result is wrong!\n");
      }  
      if(fEfficiencyesdTOFPIDD)
        tofpideffesd = fEfficiencyesdTOFPIDD->Eval(pt[0]);
      else{
        printf("TOF esd PID fit failed on conversion. The result is wrong!\n");
      }
      
      //tof pid eff x tpc pid eff x ip cut eff
      if(fIPParameterizedEff){
        if(source==0) {
          if(fEfficiencyIPCharmD){
            weight = tofpideff*trdtpcPidEfficiency*fEfficiencyIPCharmD->Eval(pt[0]);
            weightErr = 0; 
          }
          else{
            printf("Fit failed on charm IP cut efficiency\n");
            weight = tofpideff*trdtpcPidEfficiency*fEfficiencyCharmSigD->GetBinContent(i+1);
            weightErr = tofpideff*trdtpcPidEfficiency*fEfficiencyCharmSigD->GetBinError(i+1); 
          }
        } 
        else if(source==2) {
          if(fEfficiencyIPConversionD){
            weight = tofpideffesd*trdtpcPidEfficiency*fEfficiencyIPConversionD->Eval(pt[0]); 
            weightErr = 0; 
          }
          else{
            printf("Fit failed on conversion IP cut efficiency\n");
            weight = tofpideffesd*trdtpcPidEfficiency*fConversionEff->GetBinContent(i+1);
            weightErr = tofpideffesd*trdtpcPidEfficiency*fConversionEff->GetBinError(i+1);
          }
        }
        else if(source==3) {
          if(fEfficiencyIPNonhfeD){
            weight = tofpideffesd*trdtpcPidEfficiency*fEfficiencyIPNonhfeD->Eval(pt[0]); 
            weightErr = 0; 
          }
          else{
            printf("Fit failed on dalitz IP cut efficiency\n");
            weight = tofpideffesd*trdtpcPidEfficiency*fNonHFEEff->GetBinContent(i+1);
            weightErr = tofpideffesd*trdtpcPidEfficiency*fNonHFEEff->GetBinError(i+1);
          }  
        }
      }
      else{
        if(source==0){ 
          if(fEfficiencyIPCharmD){
            weight = tofpideff*trdtpcPidEfficiency*fEfficiencyIPCharmD->Eval(pt[0]);
            weightErr = 0;
          }
          else{
            printf("Fit failed on charm IP cut efficiency\n");
            weight = tofpideff*trdtpcPidEfficiency*fEfficiencyCharmSigD->GetBinContent(i+1);
            weightErr = tofpideff*trdtpcPidEfficiency*fEfficiencyCharmSigD->GetBinError(i+1);
          }
        }
        else if(source==2){
          if(fBeamType==0){
            weight = tofpideffesd*trdtpcPidEfficiency*fConversionEffbgc->GetBinContent(i+1); // conversion
            weightErr = tofpideffesd*trdtpcPidEfficiency*fConversionEffbgc->GetBinError(i+1);
          }
          else{
            weight = tofpideffesd*trdtpcPidEfficiency*fConversionEff->GetBinContent(i+1); // conversion
            weightErr = tofpideffesd*trdtpcPidEfficiency*fConversionEff->GetBinError(i+1);
          }
        }
        else if(source==3){
          if(fBeamType==0){
            weight = tofpideffesd*trdtpcPidEfficiency*fNonHFEEffbgc->GetBinContent(i+1); // conversion
            weightErr = tofpideffesd*trdtpcPidEfficiency*fNonHFEEffbgc->GetBinError(i+1);
          }
          else{ 
            weight = tofpideffesd*trdtpcPidEfficiency*fNonHFEEff->GetBinContent(i+1); // Dalitz
            weightErr = tofpideffesd*trdtpcPidEfficiency*fNonHFEEff->GetBinError(i+1);
          }
        }
      }
      
      contents[0]=pt[0];
      ibin[0]=i+1;
      
      pideff->Fill(contents,weight);
      pideff->SetBinError(ibin,weightErr);
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
AliCFDataGrid *AliHFEBeautySpectrum::GetRawBspectra2ndMethod(){
 //
 // retrieve AliCFDataGrid for raw beauty spectra obtained from fit method
    //
  Int_t nDim = 1;
  
  const Int_t nPtbinning1 = 18;//number of pt bins, according to new binning
  Double_t kPtRange[19] = {0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 12, 16, 20};
  
  Int_t nBin[1] = {nPtbinning1};
  
  AliCFDataGrid *rawBeautyContainer = new AliCFDataGrid("rawBeautyContainer","rawBeautyContainer",nDim,nBin);
      
  THnSparseF *brawspectra;
  brawspectra= new THnSparseF("brawspectra", "beauty yields ; p_{t}(GeV/c)", nDim, nBin);
  
  brawspectra->SetBinEdges(0, kPtRange);
      
  Double_t pt[1];
  Double_t yields= 0.;
  Double_t valuesb[2];
            
  for(int i=0; i<fBSpectrum2ndMethod->GetNbinsX(); i++){
    pt[0]=(kPtRange[i]+kPtRange[i+1])/2.;
    
    yields = fBSpectrum2ndMethod->GetBinContent(i+1);
        
    valuesb[0]=pt[0];
    brawspectra->Fill(valuesb,yields);
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
  TH1D* hRawBeautySpectra = (TH1D*)rawBeautyContainer->Project(0);
  
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
void AliHFEBeautySpectrum::CalculateNonHFEsyst(){
  //
  // Calculate non HFE sys
  //
  //

  if(!fNonHFEsyst)
    return;

  Double_t evtnorm[1] = {0.0};
  if(fNMCbgEvents>0) evtnorm[0]= double(fNEvents)/double(fNMCbgEvents);
  
  AliCFDataGrid *convSourceGrid[kElecBgSources-3][kBgLevels];
  AliCFDataGrid *nonHFESourceGrid[kElecBgSources-3][kBgLevels];

  AliCFDataGrid *bgLevelGrid[2][kBgLevels];//for pi0 and eta based errors
  AliCFDataGrid *bgNonHFEGrid[kBgLevels];
  AliCFDataGrid *bgConvGrid[kBgLevels];

  Int_t stepbackground = 3;
  Int_t* bins=new Int_t[1];
  const Char_t *bgBase[2] = {"pi0","eta"};
 
  bins[0]=kSignalPtBins;//fConversionEff[centrality]->GetNbinsX();
   
  AliCFDataGrid *weightedConversionContainer = new AliCFDataGrid("weightedConversionContainer","weightedConversionContainer",1,bins);
  AliCFDataGrid *weightedNonHFEContainer = new AliCFDataGrid("weightedNonHFEContainer","weightedNonHFEContainer",1,bins);

  for(Int_t iLevel = 0; iLevel < kBgLevels; iLevel++){   
    for(Int_t iSource = 0; iSource < kElecBgSources-3; iSource++){
      convSourceGrid[iSource][iLevel] = new AliCFDataGrid(Form("convGrid_%d_%d",iSource,iLevel),Form("convGrid_%d_%d",iSource,iLevel),*fConvSourceContainer[iSource][iLevel],stepbackground);
      weightedConversionContainer->SetGrid(GetPIDxIPEff(2));
      convSourceGrid[iSource][iLevel]->Multiply(weightedConversionContainer,1.0);
      
      nonHFESourceGrid[iSource][iLevel] = new AliCFDataGrid(Form("nonHFEGrid_%d_%d",iSource,iLevel),Form("nonHFEGrid_%d_%d",iSource,iLevel),*fNonHFESourceContainer[iSource][iLevel],stepbackground);
      weightedNonHFEContainer->SetGrid(GetPIDxIPEff(3));
      nonHFESourceGrid[iSource][iLevel]->Multiply(weightedNonHFEContainer,1.0);
    }
    
    bgConvGrid[iLevel] = (AliCFDataGrid*)convSourceGrid[0][iLevel]->Clone();
    for(Int_t iSource = 2; iSource < kElecBgSources-3; iSource++){
      bgConvGrid[iLevel]->Add(convSourceGrid[iSource][iLevel]);
    }
    if(!fEtaSyst)
      bgConvGrid[iLevel]->Add(convSourceGrid[1][iLevel]);
    
    bgNonHFEGrid[iLevel] = (AliCFDataGrid*)nonHFESourceGrid[0][iLevel]->Clone(); 
    for(Int_t iSource = 2; iSource < kElecBgSources-3; iSource++){//add other sources to pi0, to get overall background from all meson decays, exception: eta (independent error calculation)
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
  TH1D *hSpeciesErrors[kElecBgSources-4];
  for(Int_t iSource = 1; iSource < kElecBgSources-3; iSource++){
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
    for(Int_t iSource = 1; iSource < kElecBgSources-3; iSource++){
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
void AliHFEBeautySpectrum::UnfoldBG(AliCFDataGrid* const bgsubpectrum){
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
      ctest->Divide(2);
      ctest->cd(1);
      TH1D* htest1=(TH1D*)result->Projection(0);
      htest1->Draw();
      ctest->cd(2);
      TH1D* htest2=(TH1D*)result->Projection(1);
      htest2->Draw();
  }


  TGraphErrors* unfoldedbgspectrumD = Normalize(result);
  if(!unfoldedbgspectrumD) {
    AliError("Unfolded background spectrum doesn't exist");
  }
  else{
    TFile *file = TFile::Open("unfoldedbgspectrum.root","recreate");
    if(fBeamType==0)unfoldedbgspectrumD->Write("unfoldedbgspectrum");

    file->Close();
  }
}
