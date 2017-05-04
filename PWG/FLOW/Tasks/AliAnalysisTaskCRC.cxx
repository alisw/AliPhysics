/*************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
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

/********************************
 * analysis task for CRC        *
 *                              *
 * author: Jacopo Margutti      *
 *         (margutti@nikhef.nl) *
 * ******************************/

class TFile;
class TString;
class TList;
class AliAnalysisTaskSE;

#include "Riostream.h"
#include "AliAODEvent.h"
#include "AliAODHeader.h"
#include "AliCentrality.h"
#include "AliFlowVector.h"
#include "AliFlowEvent.h"
#include "AliFlowEventSimple.h"
#include "AliAnalysisTaskCRC.h"
#include "AliFlowAnalysisCRC.h"
#include "AliLog.h"

class AliFlowVector;
class TVector;

using std::cout;
using std::endl;
ClassImp(AliAnalysisTaskCRC)

//================================================================================================================

AliAnalysisTaskCRC::AliAnalysisTaskCRC(const char *name, Bool_t useParticleWeights):
AliAnalysisTaskSE(name),
fEvent(NULL),
fQC(NULL),
fListHistos(NULL),
fBookOnlyBasicCCH(kTRUE),
fFillMultipleControlHistograms(kFALSE),
fHarmonic(1),
fApplyCorrectionForNUA(kFALSE),
fApplyCorrectionForNUAVsM(kFALSE),
fPropagateErrorAlsoFromNIT(kFALSE),
fCalculateDiffFlow(kTRUE),
fCalculate2DDiffFlow(kFALSE),
fCalculateDiffFlowVsEta(kTRUE),
fStoreDistributions(kFALSE),
fCalculateCumulantsVsM(kFALSE),
fCalculateAllCorrelationsVsM(kFALSE),
fCalculateMixedHarmonics(kFALSE),
fCalculateMixedHarmonicsVsM(kFALSE),
fStoreControlHistograms(kFALSE),
fMinimumBiasReferenceFlow(kTRUE),
fForgetAboutCovariances(kFALSE),
fStoreVarious(kFALSE),
fExactNoRPs(0),
fUse2DHistograms(kFALSE),
fFillProfilesVsMUsingWeights(kTRUE),
fUseQvectorTerms(kFALSE),
fnBinsMult(10000),
fMinMult(0.),
fMaxMult(10000.),
fUseParticleWeights(useParticleWeights),
fUsePhiWeights(kFALSE),
fUsePtWeights(kFALSE),
fUseEtaWeights(kFALSE),
fUseTrackWeights(kFALSE),
fUsePhiEtaWeights(kFALSE),
fUsePhiEtaWeightsChDep(kFALSE),
fUsePhiEtaWeightsVtxDep(kFALSE),
fUsePhiEtaCuts(kFALSE),
fUseZDCESEMulWeights(kFALSE),
fUseZDCESESpecWeights(kFALSE),
fWeightsList(NULL),
fMultiplicityWeight(NULL),
fMultiplicityIs(AliFlowCommonConstants::kRP),
fnBinsForCorrelations(10000),
fUseBootstrap(kFALSE),
fUseBootstrapVsM(kFALSE),
fnSubsamples(10),
fCalculateCRC(kTRUE),
fCalculateCRCPt(kFALSE),
fCalculateCME(kFALSE),
fCalculateCRCInt(kFALSE),
fCalculateCRC2(kFALSE),
fCalculateCRCVZ(kFALSE),
fCalculateCRCZDC(kFALSE),
fCalculateEbEFlow(kFALSE),
fStoreZDCQVecVtxPos(kFALSE),
fCRC2nEtaBins(6),
fCalculateFlowQC(kFALSE),
fCalculateFlowZDC(kFALSE),
fCalculateFlowVZ(kFALSE),
fUseVZERO(kFALSE),
fUseZDC(kFALSE),
fRecenterZDC(kFALSE),
fDivSigma(kTRUE),
fInvertZDC(kFALSE),
fCRCTestSin(kFALSE),
fVtxRbR(kFALSE),
fUseNUAforCRC(kFALSE),
fUseCRCRecenter(kFALSE),
fCRCEtaMin(0.),
fCRCEtaMax(0.),
fnCenBin(10),
fFlowQCCenBin(100),
fFlowQCDeltaEta(0.4),
fCenBinWidth(10.),
fDataSet(""),
fInteractionRate(""),
fSelectCharge(""),
fPOIExtraWeights(""),
fCorrWeight("TPCuVZuZDCu"),
fQVecList(NULL),
fCRCZDCCalibList(NULL),
fCRCVZEROCalibList(NULL),
fCRCZDCResList(NULL),
fZDCESEList(NULL),
fCenWeightsHist(NULL),
fPhiExclZoneHist(NULL),
fQAZDCCuts(kFALSE),
fMinMulZN(1),
fMaxDevZN(5.),
fZDCGainAlpha(0.395)
{
  // constructor
  AliDebug(2,"AliAnalysisTaskCRC::AliAnalysisTaskCRC(const char *name, Bool_t useParticleWeights)");
  
  // Define input and output slots here
  // Input slot #0 works with an AliFlowEventSimple
  DefineInput(0, AliFlowEventSimple::Class());
  
  // Output slot #0 is reserved
  // Output slot #1 writes into a TList container
  DefineOutput(1, TList::Class());
  
  // Event weights:
  fMultiplicityWeight = new TString("combinations");
  
  // b) Initialize default min and max values of correlations:
  //    (Remark: The default values bellow were chosen for v2=5% and M=500)
  fMinValueOfCorrelation[0] = -0.015; // <2>_min
  fMaxValueOfCorrelation[0] = 0.03; // <2>_max
  fMinValueOfCorrelation[1] = -0.6e-3; // <4>_min
  fMaxValueOfCorrelation[1] = 0.07; // <4>_max
  fMinValueOfCorrelation[2] = -0.08e-3; // <6>_min
  fMaxValueOfCorrelation[2] = 0.015; // <6>_max
  fMinValueOfCorrelation[3] = -20.e-6; // <8>_min
  fMaxValueOfCorrelation[3] = 0.003; // <8>_max
  
  // c) Initialize default min and max values of correlation products:
  //    (Remark: The default values bellow were chosen for v2=5% and M=500)
  fMinValueOfCorrelationProduct[0] = -15.e-6; // <2><4>_min
  fMaxValueOfCorrelationProduct[0] = 0.02; // <2><4>_max
  
  // d) Initialize default min and max values of q-vector terms:
  fMinValueOfQvectorTerms[0] = 0.;
  fMaxValueOfQvectorTerms[0] = 30.;
  fMinValueOfQvectorTerms[1] = 0.;
  fMaxValueOfQvectorTerms[1] = 20.;
  fMinValueOfQvectorTerms[2] = 0.;
  fMaxValueOfQvectorTerms[2] = 200.;
  fMinValueOfQvectorTerms[3] = -30.;
  fMaxValueOfQvectorTerms[3] = 80.;
  
  for(Int_t c=0; c<10; c++) {
    fPtWeightsHist[c] = NULL;
    for(Int_t b=0; b<21; b++) {
      for(Int_t k=0; k<2; k++) {
        fEtaWeightsHist[c][b][k] = NULL;
      }
    }
  }
  for(Int_t c=0; c<2; c++) {
    for(Int_t k=0; k<2; k++) {
      fNvsCenCut[c][k] = NULL;
    }
  }
  for(Int_t k=0; k<5; k++) {
    fZDCESEMultWeightsHist[k] = NULL;
    fZDCESESpecWeightsHist[k] = NULL;
  }
  
}

//================================================================================================================

AliAnalysisTaskCRC::AliAnalysisTaskCRC():
AliAnalysisTaskSE(),
fEvent(NULL),
fQC(NULL),
fListHistos(NULL),
fBookOnlyBasicCCH(kFALSE),
fFillMultipleControlHistograms(kFALSE),
fHarmonic(1),
fApplyCorrectionForNUA(kFALSE),
fApplyCorrectionForNUAVsM(kFALSE),
fPropagateErrorAlsoFromNIT(kFALSE),
fCalculateDiffFlow(kFALSE),
fCalculate2DDiffFlow(kFALSE),
fCalculateDiffFlowVsEta(kTRUE),
fStoreDistributions(kFALSE),
fCalculateCumulantsVsM(kFALSE),
fCalculateAllCorrelationsVsM(kFALSE),
fCalculateMixedHarmonics(kFALSE),
fCalculateMixedHarmonicsVsM(kFALSE),
fStoreControlHistograms(kFALSE),
fMinimumBiasReferenceFlow(kFALSE),
fForgetAboutCovariances(kFALSE),
fStoreVarious(kFALSE),
fExactNoRPs(0),
fUse2DHistograms(kFALSE),
fFillProfilesVsMUsingWeights(kTRUE),
fUseQvectorTerms(kFALSE),
fnBinsMult(0),
fMinMult(0.),
fMaxMult(0.),
fUseParticleWeights(kFALSE),
fUsePhiWeights(kFALSE),
fUsePtWeights(kFALSE),
fUseEtaWeights(kFALSE),
fUseTrackWeights(kFALSE),
fUsePhiEtaWeights(kFALSE),
fUsePhiEtaWeightsChDep(kFALSE),
fUsePhiEtaWeightsVtxDep(kFALSE),
fUsePhiEtaCuts(kFALSE),
fUseZDCESEMulWeights(kFALSE),
fUseZDCESESpecWeights(kFALSE),
fWeightsList(NULL),
fMultiplicityWeight(NULL),
fMultiplicityIs(AliFlowCommonConstants::kRP),
fnBinsForCorrelations(0),
fUseBootstrap(kFALSE),
fUseBootstrapVsM(kFALSE),
fnSubsamples(10),
fCalculateCRC(kTRUE),
fCalculateCRCPt(kFALSE),
fCalculateCME(kFALSE),
fCalculateCRCInt(kFALSE),
fCalculateCRC2(kFALSE),
fCalculateCRCVZ(kFALSE),
fCalculateCRCZDC(kFALSE),
fCalculateEbEFlow(kFALSE),
fStoreZDCQVecVtxPos(kFALSE),
fCRC2nEtaBins(6),
fCalculateFlowQC(kFALSE),
fCalculateFlowZDC(kFALSE),
fCalculateFlowVZ(kFALSE),
fUseVZERO(kFALSE),
fUseZDC(kFALSE),
fRecenterZDC(kFALSE),
fDivSigma(kTRUE),
fInvertZDC(kFALSE),
fCRCTestSin(kFALSE),
fVtxRbR(kFALSE),
fUseNUAforCRC(kFALSE),
fUseCRCRecenter(kFALSE),
fCRCEtaMin(0.),
fCRCEtaMax(0.),
fnCenBin(10),
fFlowQCCenBin(100),
fFlowQCDeltaEta(0.4),
fCenBinWidth(10.),
fDataSet(""),
fInteractionRate(""),
fSelectCharge(""),
fPOIExtraWeights(""),
fCorrWeight("TPCuVZuZDCu"),
fQVecList(NULL),
fCRCZDCCalibList(NULL),
fCRCVZEROCalibList(NULL),
fCRCZDCResList(NULL),
fZDCESEList(NULL),
fCenWeightsHist(NULL),
fPhiExclZoneHist(NULL),
fQAZDCCuts(kFALSE),
fMinMulZN(1),
fMaxDevZN(5.),
fZDCGainAlpha(0.395)
{
  // Dummy constructor
  AliDebug(2,"AliAnalysisTaskCRC::AliAnalysisTaskCRC()");
  
  // b) Initialize default min and max values of correlations:
  //    (Remark: The default values bellow were chosen for v2=5% and M=500)
  fMinValueOfCorrelation[0] = -0.015; // <2>_min
  fMaxValueOfCorrelation[0] = 0.03; // <2>_max
  fMinValueOfCorrelation[1] = -0.6e-3; // <4>_min
  fMaxValueOfCorrelation[1] = 0.07; // <4>_max
  fMinValueOfCorrelation[2] = -0.08e-3; // <6>_min
  fMaxValueOfCorrelation[2] = 0.015; // <6>_max
  fMinValueOfCorrelation[3] = -20.e-6; // <8>_min
  fMaxValueOfCorrelation[3] = 0.003; // <8>_max
  
  // c) Initialize default min and max values of correlation products:
  //    (Remark: The default values bellow were chosen for v2=5% and M=500)
  fMinValueOfCorrelationProduct[0] = -15.e-6; // <2><4>_min
  fMaxValueOfCorrelationProduct[0] = 0.02; // <2><4>_max
  
  // d) Initialize default min and max values of q-vector terms:
  fMinValueOfQvectorTerms[0] = 0.;
  fMaxValueOfQvectorTerms[0] = 30.;
  fMinValueOfQvectorTerms[1] = 0.;
  fMaxValueOfQvectorTerms[1] = 20.;
  fMinValueOfQvectorTerms[2] = 0.;
  fMaxValueOfQvectorTerms[2] = 200.;
  fMinValueOfQvectorTerms[3] = -30.;
  fMaxValueOfQvectorTerms[3] = 80.;
  
  for(Int_t c=0; c<10; c++) {
    fPtWeightsHist[c] = NULL;
    for(Int_t b=0; b<21; b++) {
      for(Int_t k=0; k<2; k++) {
        fEtaWeightsHist[c][b][k] = NULL;
      }
    }
  }
  for(Int_t c=0; c<2; c++) {
    for(Int_t k=0; k<2; k++) {
      fNvsCenCut[c][k] = NULL;
    }
  }
  for(Int_t k=0; k<5; k++) {
    fZDCESEMultWeightsHist[k] = NULL;
    fZDCESESpecWeightsHist[k] = NULL;
  }
  
}

//==========================================================================================================

void AliAnalysisTaskCRC::UserCreateOutputObjects()
{
  // Called at every worker node to initialize
  AliDebug(2,"AliAnalysisTaskCRC::UserCreateOutputObjects()");
  
  // Analyser:
  fQC = new AliFlowAnalysisCRC("AliFlowAnalysisCRC",fnCenBin,fCenBinWidth);
  
  // Common:
  fQC->SetBookOnlyBasicCCH(fBookOnlyBasicCCH);
  fQC->SetFillMultipleControlHistograms(fFillMultipleControlHistograms);
  fQC->SetHarmonic(fHarmonic);
  fQC->SetApplyCorrectionForNUA(fApplyCorrectionForNUA);
  fQC->SetApplyCorrectionForNUAVsM(fApplyCorrectionForNUAVsM);
  fQC->SetPropagateErrorAlsoFromNIT(fPropagateErrorAlsoFromNIT);
  fQC->SetCalculateDiffFlow(fCalculateDiffFlow);
  fQC->SetCalculate2DDiffFlow(fCalculate2DDiffFlow);
  fQC->SetCalculateDiffFlowVsEta(fCalculateDiffFlowVsEta);
  fQC->SetStoreDistributions(fStoreDistributions);
  fQC->SetCalculateCumulantsVsM(fCalculateCumulantsVsM);
  fQC->SetCalculateAllCorrelationsVsM(fCalculateAllCorrelationsVsM);
  fQC->SetCalculateMixedHarmonics(fCalculateMixedHarmonics);
  fQC->SetCalculateMixedHarmonicsVsM(fCalculateMixedHarmonicsVsM);
  fQC->SetStoreControlHistograms(fStoreControlHistograms);
  fQC->SetMinimumBiasReferenceFlow(fMinimumBiasReferenceFlow);
  fQC->SetForgetAboutCovariances(fForgetAboutCovariances);
  fQC->SetExactNoRPs(fExactNoRPs);
  if(fDataSet.EqualTo("2010")) fQC->SetDataSet(AliFlowAnalysisCRC::k2010);
  if(fDataSet.EqualTo("2011")) fQC->SetDataSet(AliFlowAnalysisCRC::k2011);
  if(fDataSet.EqualTo("2015")) fQC->SetDataSet(AliFlowAnalysisCRC::k2015);
  if(fDataSet.EqualTo("2015v6")) fQC->SetDataSet(AliFlowAnalysisCRC::k2015v6);
  if(fInteractionRate.EqualTo("high")) fQC->SetInteractionRate(AliFlowAnalysisCRC::kHigh);
  if(fInteractionRate.EqualTo("low"))  fQC->SetInteractionRate(AliFlowAnalysisCRC::kLow);
  if(fInteractionRate.EqualTo("pos"))  fQC->SetInteractionRate(AliFlowAnalysisCRC::kPos);
  if(fInteractionRate.EqualTo("neg"))  fQC->SetInteractionRate(AliFlowAnalysisCRC::kNeg);
  if(fSelectCharge.EqualTo("pos")) fQC->SetSelectCharge(AliFlowAnalysisCRC::kPosCh);
  if(fSelectCharge.EqualTo("neg")) fQC->SetSelectCharge(AliFlowAnalysisCRC::kNegCh);
  fQC->SetCalculateCRC(fCalculateCRC);
  fQC->SetCalculateCRCPt(fCalculateCRCPt);
  fQC->SetCalculateCME(fCalculateCME);
  fQC->SetCalculateCRCInt(fCalculateCRCInt),
  fQC->SetCalculateCRC2(fCalculateCRC2);
  fQC->SetCalculateCRCVZ(fCalculateCRCVZ);
  fQC->SetCalculateCRCZDC(fCalculateCRCZDC);
  fQC->SetCalculateEbEFlow(fCalculateEbEFlow);
  fQC->SetStoreZDCQVecVtxPos(fStoreZDCQVecVtxPos);
  fQC->SetCRC2nEtaBins(fCRC2nEtaBins);
  fQC->SetCalculateFlowQC(fCalculateFlowQC);
  fQC->SetFlowQCCenBin(fFlowQCCenBin);
  fQC->SetFlowQCDeltaEta(fFlowQCDeltaEta);
  fQC->SetCalculateFlowZDC(fCalculateFlowZDC);
  fQC->SetCalculateFlowVZ(fCalculateFlowVZ);
  fQC->SetUseVZERO(fUseVZERO);
  fQC->SetUseZDC(fUseZDC);
  fQC->SetRecenterZDC(fRecenterZDC);
  fQC->SetDivSigma(fDivSigma);
  fQC->SetInvertZDC(fInvertZDC);
  fQC->SetQAZDCCuts(fQAZDCCuts);
  fQC->SetMinMulZN(fMinMulZN);
  fQC->SetMaxDevZN(fMaxDevZN);
  fQC->SetTestSin(fCRCTestSin);
  fQC->SetRecenterZDCVtxRbR(fVtxRbR);
  fQC->SetNUAforCRC(fUseNUAforCRC);
  fQC->SetUseCRCRecenter(fUseCRCRecenter);
  fQC->SetCRCEtaRange(fCRCEtaMin,fCRCEtaMax);
  fQC->SetUsePtWeights(fUsePtWeights);
  fQC->SetUseEtaWeights(fUseEtaWeights);
  if(fCorrWeight.Contains("TPCu")) fQC->SetCorrWeightTPC(AliFlowAnalysisCRC::kUnit);
  else if(fCorrWeight.Contains("TPCm")) fQC->SetCorrWeightTPC(AliFlowAnalysisCRC::kMultiplicity);
  if(fCorrWeight.Contains("VZu"))  fQC->SetCorrWeightVZ(AliFlowAnalysisCRC::kUnit);
  else if(fCorrWeight.Contains("VZm"))  fQC->SetCorrWeightVZ(AliFlowAnalysisCRC::kMultiplicity);
  if(fCorrWeight.Contains("ZDCu")) fQC->SetCorrWeightZDC(AliFlowAnalysisCRC::kUnit);
  else if(fCorrWeight.Contains("ZDCm")) fQC->SetCorrWeightZDC(AliFlowAnalysisCRC::kMultiplicity);
  // Multiparticle correlations vs multiplicity:
  fQC->SetnBinsMult(fnBinsMult);
  fQC->SetMinMult(fMinMult);
  fQC->SetMaxMult(fMaxMult);
  // Particle weights:
  if(fUseParticleWeights) {
    // Pass the flags to class:
    if(fUsePhiWeights){fQC->SetUsePhiWeights(fUsePhiWeights);}
    if(fUsePtWeights){fQC->SetUsePtWeights(fUsePtWeights);}
    if(fUseEtaWeights){fQC->SetUseEtaWeights(fUseEtaWeights);}
    if(fUseTrackWeights){fQC->SetUseTrackWeights(fUseTrackWeights);}
    if(fPOIExtraWeights.EqualTo("EtaPhi")) fQC->SetPOIExtraWeights(AliFlowAnalysisCRC::kEtaPhi);
    if(fPOIExtraWeights.EqualTo("EtaPhiCh")) fQC->SetPOIExtraWeights(AliFlowAnalysisCRC::kEtaPhiCh);
    if(fPOIExtraWeights.EqualTo("EtaPhiVtx")) fQC->SetPOIExtraWeights(AliFlowAnalysisCRC::kEtaPhiVtx);
    if(fPOIExtraWeights.EqualTo("EtaPhiChPt")) fQC->SetPOIExtraWeights(AliFlowAnalysisCRC::kEtaPhiChPt);
    if(fPOIExtraWeights.EqualTo("EtaPhiRbR")) fQC->SetPOIExtraWeights(AliFlowAnalysisCRC::kEtaPhiRbR);
    if(fPOIExtraWeights.EqualTo("EtaPhiChRbR")) fQC->SetPOIExtraWeights(AliFlowAnalysisCRC::kEtaPhiChRbR);
    // Pass the list with weights to class:
    if(fWeightsList) fQC->SetWeightsList(fWeightsList);
  }
  if(fUsePhiEtaCuts) fQC->SetUsePhiEtaCuts(fUsePhiEtaCuts);
  // Event weights:
  if(!fMultiplicityWeight->Contains("combinations")) {
    fQC->SetMultiplicityWeight(fMultiplicityWeight->Data());
  }
  // Q Vector weights:
  if(fUseCRCRecenter) {
    if(fQVecList) fQC->SetCRCQVecWeightsList(fQVecList);
  }
  if (fRecenterZDC) {
    if(fCRCZDCCalibList) fQC->SetCRCZDCCalibList(fCRCZDCCalibList);
    if(fCRCZDCResList) fQC->SetCRCZDCResList(fCRCZDCResList);
  }
  if(fCRCVZEROCalibList) fQC->SetCRCVZEROCalibList(fCRCVZEROCalibList);
  if (fQAZDCCuts) {
    if(fZDCESEList) fQC->SetZDCESEList(fZDCESEList);
  }
  if(fCenWeightsHist) fQC->SetCenWeightsHist(fCenWeightsHist);
  if(fPhiExclZoneHist) fQC->SetPhiExclZoneHist(fPhiExclZoneHist);
  if(fUsePtWeights){
    for(Int_t c=0; c<10; c++) {
      if(fPtWeightsHist[c]) fQC->SetPtWeightsHist(fPtWeightsHist[c],c);
    }
  }
  if(fUseEtaWeights){
    for(Int_t h=0; h<10; h++) {
      for(Int_t b=0; b<21; b++) {
        for(Int_t c=0; c<2; c++) {
          if(fEtaWeightsHist[h][b][c]) fQC->SetEtaWeightsHist(fEtaWeightsHist[h][b][c],h,b,c);
        }
      }
    }
  }
  if(fMinMulZN>1){
    for(Int_t c=0; c<2; c++) {
      for(Int_t k=0; k<2; k++) {
        if(fNvsCenCut[c][k]) fQC->SetNvsCenCut(fNvsCenCut[c][k],c,k);
      }
    }
  }
  if(fUseZDCESEMulWeights) {
    fQC->SetUseZDCESEMulWeights(fUseZDCESEMulWeights);
    for(Int_t k=0; k<5; k++) {
      if(fZDCESEMultWeightsHist[k]) fQC->SetZDCESEMultWeightsHist(fZDCESEMultWeightsHist[k],k);
    }
  }
  if(fUseZDCESESpecWeights) {
    fQC->SetUseZDCESESpecWeights(fUseZDCESESpecWeights);
    for(Int_t k=0; k<5; k++) {
      if(fZDCESESpecWeightsHist[k]) fQC->SetZDCESESpecWeightsHist(fZDCESESpecWeightsHist[k],k);
    }
  }
  fQC->SetZDCGainAlpha(fZDCGainAlpha);
  fQC->SetMultiplicityIs(fMultiplicityIs);
  fQC->SetnBinsForCorrelations(fnBinsForCorrelations);
  fQC->SetUse2DHistograms(fUse2DHistograms);
  fQC->SetFillProfilesVsMUsingWeights(fFillProfilesVsMUsingWeights);
  fQC->SetUseQvectorTerms(fUseQvectorTerms);
  
  // Store phi distribution for one event to illustrate flow:
  fQC->SetStoreVarious(fStoreVarious);
  
  // Initialize default min and max values of correlations:
  for(Int_t ci=0;ci<4;ci++) {
    fQC->SetMinValueOfCorrelation(ci,fMinValueOfCorrelation[ci]);
    fQC->SetMaxValueOfCorrelation(ci,fMaxValueOfCorrelation[ci]);
  }
  
  // Initialize default min and max values of correlation products:
  for(Int_t cpi=0;cpi<1;cpi++) {
    fQC->SetMinValueOfCorrelationProduct(cpi,fMinValueOfCorrelationProduct[cpi]);
    fQC->SetMaxValueOfCorrelationProduct(cpi,fMaxValueOfCorrelationProduct[cpi]);
  }
  
  // Initialize default min and max values of Q-vector terms:
  for(Int_t ci=0;ci<4;ci++) {
    fQC->SetMinValueOfQvectorTerms(ci,fMinValueOfQvectorTerms[ci]);
    fQC->SetMaxValueOfQvectorTerms(ci,fMaxValueOfQvectorTerms[ci]);
  }
  
  // Bootstrap:
  fQC->SetUseBootstrap(fUseBootstrap);
  fQC->SetUseBootstrapVsM(fUseBootstrapVsM);
  fQC->SetnSubsamples(fnSubsamples);
  
  fQC->Init();
  
  if(fQC->GetHistList()) {
    fListHistos = fQC->GetHistList();
  } else {
    Printf("ERROR: Could not retrieve histogram list (QC, Task::UserCreateOutputObjects()) !!!!");
  }
  
  PostData(1,fListHistos);
  
} // end of void AliAnalysisTaskCRC::UserCreateOutputObjects()

//================================================================================================================

void AliAnalysisTaskCRC::UserExec(Option_t *)
{
  // main loop (called for each event)
  fEvent = dynamic_cast<AliFlowEvent*>(GetInputData(0));
  
  // Q-cumulants
  if(fEvent) {
    fQC->SetRunNumber(fEvent->GetRun());
    fQC->Make(fEvent);
  } else {
    cout<<"WARNING: No input data (QC, Task::UserExec()) !!!!"<<endl;
    cout<<endl;
  }
  
  PostData(1,fListHistos);
}

//================================================================================================================

void AliAnalysisTaskCRC::Terminate(Option_t *)
{
  //accessing the merged output list:
  fListHistos = (TList*)GetOutputData(1);
  
  fQC = new AliFlowAnalysisCRC("AliFlowAnalysisCRC",fnCenBin,fCenBinWidth);
  if(fDataSet.EqualTo("2010")) fQC->SetDataSet(AliFlowAnalysisCRC::k2010);
  if(fDataSet.EqualTo("2011")) fQC->SetDataSet(AliFlowAnalysisCRC::k2011);
  if(fDataSet.EqualTo("2015")) fQC->SetDataSet(AliFlowAnalysisCRC::k2015);
  if(fDataSet.EqualTo("2015v6")) fQC->SetDataSet(AliFlowAnalysisCRC::k2015v6);
  if(fInteractionRate.EqualTo("high")) fQC->SetInteractionRate(AliFlowAnalysisCRC::kHigh);
  if(fInteractionRate.EqualTo("low"))  fQC->SetInteractionRate(AliFlowAnalysisCRC::kLow);
  if(fInteractionRate.EqualTo("pos"))  fQC->SetInteractionRate(AliFlowAnalysisCRC::kPos);
  if(fInteractionRate.EqualTo("neg"))  fQC->SetInteractionRate(AliFlowAnalysisCRC::kNeg);
  fQC->SetRunList();
  
  if(fListHistos) {
    fQC->GetOutputHistograms(fListHistos);
    fQC->Finish();
    PostData(1,fListHistos);
  } else {
    cout<<" WARNING: histogram list pointer is empty (QC, Task::Terminate()) !!!!"<<endl;
    cout<<endl;
  }
  
} // end of void AliAnalysisTaskCRC::Terminate(Option_t *)





















