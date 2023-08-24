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

#include "TList.h"
#include "TFile.h"
#include "Riostream.h"
#include "AliAODEvent.h"
#include "AliAODHeader.h"
#include "AliCentrality.h"
#include "AliFlowVector.h"
#include "AliFlowEvent.h"
#include "TProfile2D.h"
#include "AliFlowEventSimple.h"
#include "AliAnalysisTaskQvec.h"
#include "AliFlowAnalysisQvec.h"
#include "AliLog.h"

class AliFlowVector;
class TVector;

using std::cout;
using std::endl;
ClassImp(AliAnalysisTaskQvec)

//================================================================================================================

AliAnalysisTaskQvec::AliAnalysisTaskQvec(const char *name, Bool_t useParticleWeights):
AliAnalysisTaskSE(name),
fEvent(NULL),
fQC(NULL),
fListHistos(NULL),
fExactNoRPs(0),
fUseParticleWeights(useParticleWeights),
fUsePtWeights(kFALSE),
fUsePhiEtaCuts(kFALSE),
fUseZDCESEMulWeights(kFALSE),
fUseZDCESESpecWeights(kFALSE),
fCutMultiplicityOutliers(kFALSE),
fWeightsList(NULL),
fMultiplicityWeight(NULL),
fCalculateCRC(kTRUE),
fCalculateCME(kFALSE),
fUseVZERO(kFALSE),
fUseZDC(kFALSE),
fRemoveSplitMergedTracks(kFALSE),
fRecenterZDC(kFALSE),
fDivSigma(kTRUE),
fInvertZDC(kFALSE),
fCRCEtaMin(0.),
fCRCEtaMax(0.),
fStoreZDCQVecVtxPos(kFALSE),
fCRCTestSin(kFALSE),
fnCenBin(10),
fCenBinWidth(10.),
fDataSet(""),
fInteractionRate(""),
fSelectCharge(""),
fPOIExtraWeights(""),
fCRCZDCCalibList(NULL),
fCRCZDC2DCutList(NULL),
fCRCVZEROCalibList(NULL),
fCRCZDCResList(NULL),
fZDCESEList(NULL),
fZDCCalibListFinalCommonPart(NULL),
fCenWeightsHist(NULL),
fRefMultRbRPro(NULL),
fAvEZDCCRbRPro(NULL),
fAvEZDCARbRPro(NULL),
fQAZDCCuts(kFALSE),
fUseTracklets(kFALSE),
fMinMulZN(1)
{
  // constructor
  AliDebug(2,"AliAnalysisTaskQvec::AliAnalysisTaskQvec(const char *name, Bool_t useParticleWeights)");

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

  // c) Initialize default min and max values of correlation products:
  //    (Remark: The default values bellow were chosen for v2=5% and M=500)

  // d) Initialize default min and max values of q-vector terms:

  for(Int_t c=0; c<10; c++) {
    fPtWeightsHist[c] = NULL;
  }
  for(Int_t k=0; k<5; k++) {
    fZDCESEMultWeightsHist[k] = NULL;
    fZDCESESpecWeightsHist[k] = NULL;
  }

}

//================================================================================================================

AliAnalysisTaskQvec::AliAnalysisTaskQvec():
AliAnalysisTaskSE(),
fEvent(NULL),
fQC(NULL),
fListHistos(NULL),
fExactNoRPs(0),
fUseParticleWeights(kFALSE),
fUsePtWeights(kFALSE),
fUsePhiEtaCuts(kFALSE),
fUseZDCESEMulWeights(kFALSE),
fUseZDCESESpecWeights(kFALSE),
fCutMultiplicityOutliers(kFALSE),
fWeightsList(NULL),
fMultiplicityWeight(NULL),
fCalculateCRC(kTRUE),
fCalculateCME(kFALSE),
fUseVZERO(kFALSE),
fUseZDC(kFALSE),
fRemoveSplitMergedTracks(kFALSE),
fRecenterZDC(kFALSE),
fDivSigma(kTRUE),
fInvertZDC(kFALSE),
fCRCEtaMin(0.),
fCRCEtaMax(0.),
fStoreZDCQVecVtxPos(kFALSE),
fCRCTestSin(kFALSE),
fnCenBin(10),
fCenBinWidth(10.),
fDataSet(""),
fInteractionRate(""),
fSelectCharge(""),
fPOIExtraWeights(""),
fCRCZDCCalibList(NULL),
fCRCZDC2DCutList(NULL),
fCRCVZEROCalibList(NULL),
fCRCZDCResList(NULL),
fZDCESEList(NULL),
fZDCCalibListFinalCommonPart(NULL),
fCenWeightsHist(NULL),
fRefMultRbRPro(NULL),
fAvEZDCCRbRPro(NULL),
fAvEZDCARbRPro(NULL),
fQAZDCCuts(kFALSE),
fUseTracklets(kFALSE),
fMinMulZN(1)
{
  // Dummy constructor
  AliDebug(2,"AliAnalysisTaskQvec::AliAnalysisTaskQvec()");

  // b) Initialize default min and max values of correlations:
  //    (Remark: The default values bellow were chosen for v2=5% and M=500)

  // c) Initialize default min and max values of correlation products:
  //    (Remark: The default values bellow were chosen for v2=5% and M=500)

  // d) Initialize default min and max values of q-vector terms:

  for(Int_t c=0; c<10; c++) {
    fPtWeightsHist[c] = NULL;
  }
  for(Int_t k=0; k<5; k++) {
    fZDCESEMultWeightsHist[k] = NULL;
    fZDCESESpecWeightsHist[k] = NULL;
  }

}

//==========================================================================================================

void AliAnalysisTaskQvec::UserCreateOutputObjects()
{
  // Called at every worker node to initialize
  AliDebug(2,"AliAnalysisTaskQvec::UserCreateOutputObjects()");

  // Analyser:
  fQC = new AliFlowAnalysisQvec("AliFlowAnalysisQvec",fnCenBin,fCenBinWidth);

  // Common:
  fQC->SetExactNoRPs(fExactNoRPs);
  if(fDataSet.EqualTo("2010")) fQC->SetDataSet(AliFlowAnalysisQvec::k2010);
  if(fDataSet.EqualTo("2011")) fQC->SetDataSet(AliFlowAnalysisQvec::k2011);
  if(fDataSet.EqualTo("2015")) fQC->SetDataSet(AliFlowAnalysisQvec::k2015);
  if(fDataSet.EqualTo("2015v6")) fQC->SetDataSet(AliFlowAnalysisQvec::k2015v6);
  if(fDataSet.EqualTo("2015pidfix")) fQC->SetDataSet(AliFlowAnalysisQvec::k2015pidfix);
  if(fInteractionRate.EqualTo("high")) fQC->SetInteractionRate(AliFlowAnalysisQvec::kHigh);
  if(fInteractionRate.EqualTo("low"))  fQC->SetInteractionRate(AliFlowAnalysisQvec::kLow);
  if(fInteractionRate.EqualTo("pos"))  fQC->SetInteractionRate(AliFlowAnalysisQvec::kPos);
  if(fInteractionRate.EqualTo("neg"))  fQC->SetInteractionRate(AliFlowAnalysisQvec::kNeg);
  if(fSelectCharge.EqualTo("pos")) fQC->SetSelectCharge(AliFlowAnalysisQvec::kPosCh);
  if(fSelectCharge.EqualTo("neg")) fQC->SetSelectCharge(AliFlowAnalysisQvec::kNegCh);
  fQC->SetCalculateCRC(fCalculateCRC);
  fQC->SetCalculateCME(fCalculateCME);
  fQC->SetStoreZDCQVecVtxPos(fStoreZDCQVecVtxPos);
  fQC->SetUseVZERO(fUseVZERO);
  fQC->SetUseZDC(fUseZDC);
  fQC->SetRemoveSplitMergedTracks(fRemoveSplitMergedTracks);
  fQC->SetRecenterZDC(fRecenterZDC);
  fQC->SetDivSigma(fDivSigma);
  fQC->SetInvertZDC(fInvertZDC);
  fQC->SetQAZDCCuts(fQAZDCCuts);
  fQC->SetUseTracklets(fUseTracklets);
  fQC->SetMinMulZN(fMinMulZN);
  fQC->SetTestSin(fCRCTestSin);
  fQC->SetUsePtWeights(fUsePtWeights);
  fQC->SetCRCEtaRange(fCRCEtaMin,fCRCEtaMax);
  
  // Multiparticle correlations vs multiplicity:
  
  // Particle weights:
  if(fUseParticleWeights) {
    // Pass the flags to class:
    if(fUsePtWeights){fQC->SetUsePtWeights(fUsePtWeights);}
    if(fPOIExtraWeights.EqualTo("EtaPhi")) fQC->SetPOIExtraWeights(AliFlowAnalysisQvec::kEtaPhi);
    if(fPOIExtraWeights.EqualTo("EtaPhiCh")) fQC->SetPOIExtraWeights(AliFlowAnalysisQvec::kEtaPhiCh);
    if(fPOIExtraWeights.EqualTo("EtaPhiVtx")) fQC->SetPOIExtraWeights(AliFlowAnalysisQvec::kEtaPhiVtx);
    if(fPOIExtraWeights.EqualTo("EtaPhiChPt")) fQC->SetPOIExtraWeights(AliFlowAnalysisQvec::kEtaPhiChPt);
    if(fPOIExtraWeights.EqualTo("EtaPhiRbR")) fQC->SetPOIExtraWeights(AliFlowAnalysisQvec::kEtaPhiRbR);
    if(fPOIExtraWeights.EqualTo("EtaPhiChRbR")) fQC->SetPOIExtraWeights(AliFlowAnalysisQvec::kEtaPhiChRbR);
    if(fPOIExtraWeights.EqualTo("EtaPhiVtxRbR")) fQC->SetPOIExtraWeights(AliFlowAnalysisQvec::kEtaPhiVtxRbR);
    // Pass the list with weights to class:
    if(fWeightsList) fQC->SetWeightsList(fWeightsList);
  }
  if(fUsePhiEtaCuts) fQC->SetUsePhiEtaCuts(fUsePhiEtaCuts);
  // Event weights:
  if(!fMultiplicityWeight->Contains("combinations")) {
    fQC->SetMultiplicityWeight(fMultiplicityWeight->Data());
  }
  // Q Vector weights:

  if (fRecenterZDC) {
    if(fCRCZDCCalibList) fQC->SetCRCZDCCalibList(fCRCZDCCalibList);
    if(fCRCZDC2DCutList) fQC->SetCRCZDC2DCutList(fCRCZDC2DCutList);
    if(fCRCZDCResList) fQC->SetCRCZDCResList(fCRCZDCResList);
    //@Shi set my ZDC calib file
    if(fZDCCalibListFinalCommonPart) {
		fQC->SetZDCCalibListFinalCommonPart(fZDCCalibListFinalCommonPart);
	}
  }
  if(fCRCVZEROCalibList) fQC->SetCRCVZEROCalibList(fCRCVZEROCalibList);
  if (fQAZDCCuts) {
    if(fZDCESEList) fQC->SetZDCESEList(fZDCESEList);
  }
  if(fCenWeightsHist) fQC->SetCenWeightsHist(fCenWeightsHist);
  if(fRefMultRbRPro) fQC->SetRefMultRbRPro(fRefMultRbRPro);
  if(fAvEZDCCRbRPro && fAvEZDCARbRPro) {
    fQC->SetAvEZDCRbRPro(fAvEZDCCRbRPro,fAvEZDCARbRPro);
  }
  cout<<"===> fUsePtWeights === "<<fUsePtWeights<<endl;
  if(fUsePtWeights){
	cout<<"===> setting pT weight Hist"<<endl;
    for(Int_t c=0; c<10; c++) {
      if(fPtWeightsHist[c]) fQC->SetPtWeightsHist(fPtWeightsHist[c],c);
      cout<<"===> fPtWeightsHist["<<c<<"]->GetNbinsX() = "<<fPtWeightsHist[c]->GetNbinsX()<<endl;
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
	  cout<<"===>  fZDCESESpecWeightsHist[k]->GetNbinsX() == "<<fZDCESESpecWeightsHist[k]->GetNbinsX()<<endl;
      if(fZDCESESpecWeightsHist[k]) fQC->SetZDCESESpecWeightsHist(fZDCESESpecWeightsHist[k],k);
    }
  }
  fQC->SetCutMultiplicityOutliers(fCutMultiplicityOutliers);

  // Store phi distribution for one event to illustrate flow:

  // Initialize default min and max values of correlations:

  // Initialize default min and max values of correlation products:

  // Initialize default min and max values of Q-vector terms:

  // Bootstrap:
  fQC->Init();

  if(fQC->GetHistList()) {
    fListHistos = fQC->GetHistList();
  } else {
    Printf("ERROR: Could not retrieve histogram list (QC, Task::UserCreateOutputObjects()) !!!!");
  }
  PostData(1,fListHistos);

} // end of void AliAnalysisTaskQvec::UserCreateOutputObjects()

//================================================================================================================

void AliAnalysisTaskQvec::UserExec(Option_t *)
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

void AliAnalysisTaskQvec::Terminate(Option_t *)
{
  //accessing the merged output list:
  fListHistos = (TList*)GetOutputData(1);

  fQC = new AliFlowAnalysisQvec("AliFlowAnalysisQvec",fnCenBin,fCenBinWidth);
  if(fDataSet.EqualTo("2010")) fQC->SetDataSet(AliFlowAnalysisQvec::k2010);
  if(fDataSet.EqualTo("2011")) fQC->SetDataSet(AliFlowAnalysisQvec::k2011);
  if(fDataSet.EqualTo("2015")) fQC->SetDataSet(AliFlowAnalysisQvec::k2015);
  if(fDataSet.EqualTo("2015v6")) fQC->SetDataSet(AliFlowAnalysisQvec::k2015v6);
  if(fDataSet.EqualTo("2015pidfix")) fQC->SetDataSet(AliFlowAnalysisQvec::k2015pidfix);
  if(fInteractionRate.EqualTo("high")) fQC->SetInteractionRate(AliFlowAnalysisQvec::kHigh);
  if(fInteractionRate.EqualTo("low"))  fQC->SetInteractionRate(AliFlowAnalysisQvec::kLow);
  if(fInteractionRate.EqualTo("pos"))  fQC->SetInteractionRate(AliFlowAnalysisQvec::kPos);
  if(fInteractionRate.EqualTo("neg"))  fQC->SetInteractionRate(AliFlowAnalysisQvec::kNeg);
  fQC->SetRunList();

  if(fListHistos) {
    fQC->GetOutputHistograms(fListHistos);
    fQC->Finish();
    PostData(1,fListHistos);
  } else {
    cout<<" WARNING: histogram list pointer is empty (QC, Task::Terminate()) !!!!"<<endl;
    cout<<endl;
  }

} // end of void AliAnalysisTaskQvec::Terminate(Option_t *)

//================================================================================================================

void AliAnalysisTaskQvec::NotifyRun()
{
  //open file
  TGrid::Connect("alien://");
  TString ZDCRecenterFileName = Form("alien:///alice/cern.ch/user/s/sqiu/15o_ZDCRunByRunCalib/15o_ZDCcalibVar_%d.root",fCurrentRunNumber);
  TFile* ZDCRecenterFileRunByRun = TFile::Open(ZDCRecenterFileName, "READ");
  
  if(ZDCRecenterFileRunByRun) {
    TList* ZDCRecenterListRunByRun = (TList*)(ZDCRecenterFileRunByRun->FindObjectAny("Q Vectors")); // hardcoded TList Q Vectors
    if(ZDCRecenterListRunByRun) {
	  fQC->SetZDCCalibListFinalRunByRun(ZDCRecenterListRunByRun);
	} else {
      std::cout << "ERROR: ZDCRecenterList do not exist!" << std::endl;
      exit(1);
    }
  } else {
	std::cout << "ERROR: if fStepZDCRecenter larger than 0, ZDCRecenterFile should exist!" << std::endl;
    exit(1);
  }

  delete ZDCRecenterFileRunByRun;

}
