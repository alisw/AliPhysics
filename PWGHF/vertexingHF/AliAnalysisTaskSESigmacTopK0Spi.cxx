/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appeuear in the supporting documentation. The authors make no claims   *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id: */

//
//
//               Base class for Lc2V0 Analysis to be used with TMVA
//
//

//------------------------------------------------------------------------------------------
//
//  Author: C. Zampolli, A. Alici
//
//------------------------------------------------------------------------------------------

#include <TSystem.h>
#include <TParticle.h>
#include <TParticlePDG.h>
#include <TH1F.h>
#include <TH1I.h>
#include <TProfile.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TTree.h>
#include "TROOT.h"
#include <TDatabasePDG.h>
#include <TMath.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <THnSparse.h>
#include <TObjArray.h>

#include <AliAnalysisDataSlot.h>
#include <AliAnalysisDataContainer.h>
#include "AliMCEvent.h"
#include "AliAnalysisManager.h"
#include "AliAODMCHeader.h"
#include "AliAODHandler.h"
#include "AliLog.h"
#include "AliAODVertex.h"
#include "AliAODRecoDecay.h"
#include "AliAODRecoDecayHF.h"
#include "AliAODRecoCascadeHF.h"
#include "AliAnalysisVertexingHF.h"
#include "AliESDtrack.h"
#include "AliAODTrack.h"
#include "AliAODv0.h"
#include "AliAODMCParticle.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskSESigmacTopK0Spi.h"
#include "AliNormalizationCounter.h"
#include "AliAODPidHF.h"
#include "AliPIDResponse.h"
#include "AliPIDCombined.h"
#include "AliTOFPIDResponse.h"
#include "AliInputEventHandler.h"
#include "AliAODRecoDecayHF3Prong.h"
#include "AliKFParticle.h"
#include "AliKFVertex.h"
#include "AliExternalTrackParam.h"
#include "AliESDUtils.h"
#include "AliDataFile.h"
#include <TMVA/Tools.h>
#include <TMVA/Reader.h>
#include <TMVA/MethodCuts.h>

#include "IClassifierReader.h"

using std::cout;
using std::endl;

#include <dlfcn.h>

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskSESigmacTopK0Spi);
/// \endcond

//__________________________________________________________________________
AliAnalysisTaskSESigmacTopK0Spi::AliAnalysisTaskSESigmacTopK0Spi():
  AliAnalysisTaskSE(),
  fUseMCInfo(kFALSE),
  fOutput(0),
  fPIDResponse(0),
  fPIDCombined(0),
  fIsK0sAnalysis(kFALSE),
  fCounter(0),
  fCounterC(0),
  fAnalCuts(0),
  fListCuts(0),
  fListWeight(0),
  fListCounters(0),
  fListProfiles(0),
  fUseOnTheFlyV0(kFALSE),
  fIsEventSelected(kFALSE),
  fVariablesTreeSgn(0),
  fVariablesTreeBkg(0),
  fCandidateVariables(),
  fHistoCentrality(0),
  fHistoEvents(0),
  fHistoTracklets_1(0),
  fHistoTracklets_1_cent(0),
  fHistoTracklets_All(0),
  fHistoTracklets_All_cent(0),
  fHistoLc(0),
  fHistoLcOnTheFly(0),
  fFillOnlySgn(kFALSE),
  fHistoLcBeforeCuts(0),
  fHistoFiducialAcceptance(0),
  fHistoCodesSgn(0),
  fHistoCodesBkg(0),
  fVtx1(0),
  fmcLabelLc(-1),
  fKeepingOnlyHIJINGBkg(kFALSE),
  fUtils(0),
  fHistoBackground(0),
  fCurrentEvent(-1),
  fBField(0),
  fKeepingOnlyPYTHIABkg(kFALSE),
  fHistoMCLcK0SpGen(0),
  fHistoMCLcK0SpGenAcc(0),
  fHistoMCLcK0SpGenLimAcc(0),
  fTriggerMask(0),
  fFuncWeightPythia(0),
  fFuncWeightFONLL5overLHC13d3(0),
  fFuncWeightFONLL5overLHC13d3Lc(0),
  fHistoMCNch(0),
  fNTracklets_1(0),
  fNTracklets_All(0),
  fCentrality(0),
  fFillTree(0),
  fUseWeightsLibrary(kFALSE),
  fBDTReader(0),
  fTMVAlibName(""),
  fTMVAlibPtBin(""),
  fNamesTMVAVar(""),
  fBDTHisto(0),
  fBDTHistoVsMassK0S(0),
  fBDTHistoVstImpParBach(0),
  fBDTHistoVstImpParV0(0),
  fBDTHistoVsBachelorPt(0),
  fBDTHistoVsCombinedProtonProb(0),
  fBDTHistoVsCtau(0),
  fBDTHistoVsCosPAK0S(0),
  fBDTHistoVsSignd0(0),
  fBDTHistoVsCosThetaStar(0),
  fBDTHistoVsnSigmaTPCpr(0),
  fBDTHistoVsnSigmaTOFpr(0),
  fBDTHistoVsnSigmaTPCpi(0),
  fBDTHistoVsnSigmaTPCka(0),
  fBDTHistoVsBachelorP(0),
  fBDTHistoVsBachelorTPCP(0),
  fHistoNsigmaTPC(0),
  fHistoNsigmaTOF(0),
  fDebugHistograms(kFALSE),
  fAODProtection(1),
  fUsePIDresponseForNsigma(kFALSE),
  fNVars(14),
  fTimestampCut(0),
  fUseXmlWeightsFile(kTRUE),
  fReader(0),
  fVarsTMVA(0),
  fNVarsSpectators(0),
  fVarsTMVASpectators(0),
  fNamesTMVAVarSpectators(""),
  fXmlWeightsFile(""),
  fBDTHistoTMVA(0),  
  fRefMult(9.26),
  fYearNumber(16),
  fHistoNtrUnCorr(0),
  fHistoNtrCorr(0),
  fHistoVzVsNtrUnCorr(0),
  fHistoVzVsNtrCorr(0),
  fUseMultCorrection(kFALSE),
  fMultiplicityEstimator(kNtrk10),
  fDoVZER0ParamVertexCorr(1),
  fUseMultiplicityCut(kFALSE),
  fMultiplicityCutMin(0.),
  fMultiplicityCutMax(99999.),
  fUseXmlFileFromCVMFS(kFALSE),
  fXmlFileFromCVMFS(""),
  ftrackArraySelSoftPi(0x0),
  fnSelSoftPi(),
  fESDtrackCutsSoftPion(0x0),
  fhistMCSpectrumAccLc(0x0),
  fhistMCSpectrumAccLcFromSc(0x0),
  fhistMCSpectrumAccSc(0x0),
  fhSparseAnalysisSigma(0x0),
  fLcMassWindowForSigmaC(0.20),
  fSigmaCDeltaMassWindow(0.230),
  fNRotations(0),
  fMinAngleForRot(5*TMath::Pi()/6),
  fMaxAngleForRot(7*TMath::Pi()/6),
  fMinPtSigmacCand(0),
  fMaxPtSigmacCand(999),
  fOutputSparse(0),
  isLcAnalysis(0)
{
  /// Default ctor
  //
  for(Int_t i=0; i<14; i++) fMultEstimatorAvg[i]=0;
}
//___________________________________________________________________________
AliAnalysisTaskSESigmacTopK0Spi::AliAnalysisTaskSESigmacTopK0Spi(const Char_t* name,
									     AliRDHFCutsLctoV0* analCuts, Bool_t useOnTheFly) :
  AliAnalysisTaskSE(name),
  fUseMCInfo(kFALSE),
  fOutput(0),
  fPIDResponse(0),
  fPIDCombined(0),
  fIsK0sAnalysis(kFALSE),
  fCounter(0),
  fCounterC(0),
  fAnalCuts(analCuts),
  fListCuts(0),
  fListWeight(0),
  fListCounters(0),
  fListProfiles(0),
  fUseOnTheFlyV0(useOnTheFly),
  fIsEventSelected(kFALSE),
  fVariablesTreeSgn(0),
  fVariablesTreeBkg(0),
  fCandidateVariables(),
  fHistoCentrality(0),
  fHistoEvents(0),
  fHistoTracklets_1(0),
  fHistoTracklets_1_cent(0),
  fHistoTracklets_All(0),
  fHistoTracklets_All_cent(0),
  fHistoLc(0),
  fHistoLcOnTheFly(0),
  fFillOnlySgn(kFALSE),
  fHistoLcBeforeCuts(0),
  fHistoFiducialAcceptance(0),
  fHistoCodesSgn(0),
  fHistoCodesBkg(0),
  fVtx1(0),
  fmcLabelLc(-1),
  fKeepingOnlyHIJINGBkg(kFALSE),
  fUtils(0),
  fHistoBackground(0),
  fCurrentEvent(-1),
  fBField(0),
  fKeepingOnlyPYTHIABkg(kFALSE),
  fHistoMCLcK0SpGen(0),
  fHistoMCLcK0SpGenAcc(0),
  fHistoMCLcK0SpGenLimAcc(0),
  fTriggerMask(0),
  fFuncWeightPythia(0),
  fFuncWeightFONLL5overLHC13d3(0),
  fFuncWeightFONLL5overLHC13d3Lc(0),
  fHistoMCNch(0),
  fNTracklets_1(0),
  fNTracklets_All(0),
  fCentrality(0),  
  fFillTree(0),
  fUseWeightsLibrary(kFALSE),
  fBDTReader(0),
  fTMVAlibName(""),
  fTMVAlibPtBin(""),
  fNamesTMVAVar(""),
  fBDTHisto(0),
  fBDTHistoVsMassK0S(0),
  fBDTHistoVstImpParBach(0),
  fBDTHistoVstImpParV0(0),
  fBDTHistoVsBachelorPt(0),
  fBDTHistoVsCombinedProtonProb(0),
  fBDTHistoVsCtau(0),
  fBDTHistoVsCosPAK0S(0),
  fBDTHistoVsSignd0(0),
  fBDTHistoVsCosThetaStar(0),
  fBDTHistoVsnSigmaTPCpr(0),
  fBDTHistoVsnSigmaTOFpr(0),
  fBDTHistoVsnSigmaTPCpi(0),
  fBDTHistoVsnSigmaTPCka(0),
  fBDTHistoVsBachelorP(0),
  fBDTHistoVsBachelorTPCP(0),
  fHistoNsigmaTPC(0),
  fHistoNsigmaTOF(0),
  fDebugHistograms(kFALSE),
  fAODProtection(1),
  fUsePIDresponseForNsigma(kFALSE),
  fNVars(14),
  fTimestampCut(0),
  fUseXmlWeightsFile(kTRUE),
  fReader(0),
  fVarsTMVA(0),
  fNVarsSpectators(0),
  fVarsTMVASpectators(0),
  fNamesTMVAVarSpectators(""),
  fXmlWeightsFile(""),
  fBDTHistoTMVA(0),  
  fRefMult(9.26),
  fYearNumber(16),
  fHistoNtrUnCorr(0),
  fHistoNtrCorr(0),
  fHistoVzVsNtrUnCorr(0),
  fHistoVzVsNtrCorr(0),
  fUseMultCorrection(kFALSE),
  fMultiplicityEstimator(kNtrk10),
  fDoVZER0ParamVertexCorr(1),
  fUseMultiplicityCut(kFALSE),
  fMultiplicityCutMin(0.),
  fMultiplicityCutMax(99999.),
  fUseXmlFileFromCVMFS(kFALSE),
  fXmlFileFromCVMFS(""),
  ftrackArraySelSoftPi(0x0),
  fnSelSoftPi(),
  fESDtrackCutsSoftPion(0x0),
  fhistMCSpectrumAccLc(0x0),
  fhistMCSpectrumAccLcFromSc(0x0),
  fhistMCSpectrumAccSc(0x0),
  fhSparseAnalysisSigma(0x0),
  fLcMassWindowForSigmaC(0.20),
  fSigmaCDeltaMassWindow(0.230),
  fNRotations(0),
  fMinAngleForRot(5*TMath::Pi()/6),
  fMaxAngleForRot(7*TMath::Pi()/6),
  fMinPtSigmacCand(0),
  fMaxPtSigmacCand(999),
  fOutputSparse(0),
  isLcAnalysis(0)
{
  //
  /// Constructor. Initialization of Inputs and Outputs
  //
  Info("AliAnalysisTaskSESigmacTopK0Spi","Calling Constructor");

  for(Int_t i=0; i<14; i++) fMultEstimatorAvg[i]=0;

  DefineOutput(1, TList::Class());  // Tree signal + Tree Bkg + histoEvents
  DefineOutput(2, TList::Class()); // normalization counter object
  DefineOutput(3, TList::Class());  // Cuts
  DefineOutput(4, TTree::Class());  // Tree signal + Tree Bkg + histoEvents
  DefineOutput(5, TTree::Class());  // Tree signal + Tree Bkg + histoEvents
  DefineOutput(6, TList::Class());  // Sparse
  DefineOutput(7, TList::Class());  // weights
}

//___________________________________________________________________________
AliAnalysisTaskSESigmacTopK0Spi::~AliAnalysisTaskSESigmacTopK0Spi() {
  //
  /// destructor
  //
  Info("~AliAnalysisTaskSESigmacTopK0Spi","Calling Destructor");

  if (fOutput) {
    delete fOutput;
    fOutput = 0;
  }

  if (fPIDResponse) {
    delete  fPIDResponse;
  }

  if (fPIDCombined) {
    delete  fPIDCombined;
  }

  if (fCounter) {
    delete fCounter;
    fCounter = 0;
  }

  if (fCounterC) {
    delete fCounterC;
    fCounterC = 0;
  }

  if (fAnalCuts) {
    delete fAnalCuts;
    fAnalCuts = 0;
  }

  if (fListCuts) {
    delete fListCuts;
    fListCuts = 0;
  }

  if (fListCounters) {
    delete fListCounters;
    fListCounters = 0;
  }

  if (fListWeight) {
    delete fListWeight;
    fListWeight = 0;
  }

  if(fVariablesTreeSgn){
    delete fVariablesTreeSgn;
    fVariablesTreeSgn = 0;
  }

  if(fVariablesTreeBkg){
    delete fVariablesTreeBkg;
    fVariablesTreeBkg = 0;
  }

  if (fUtils) {
    delete fUtils;
    fUtils = 0;
  }
  
  if (fBDTReader) {
    //delete fBDTReader;
    fBDTReader = 0;
  }

  if (fReader) {
    delete fReader;
    fReader = 0;
  }

  if (fVarsTMVA) {
    delete fVarsTMVA;
    fVarsTMVA = 0;
  }

  if (fVarsTMVASpectators) {
    delete fVarsTMVASpectators;
    fVarsTMVASpectators = 0;
  }

  for(Int_t i=0; i<14; i++){
    if (fMultEstimatorAvg[i]) delete fMultEstimatorAvg[i];
  }

  if (fOutputSparse) {
    delete fOutputSparse;
    fOutputSparse = 0;
  }

  if(ftrackArraySelSoftPi) delete ftrackArraySelSoftPi;
  if(fESDtrackCutsSoftPion) delete fESDtrackCutsSoftPion;
  if(fhistMCSpectrumAccLc) delete fhistMCSpectrumAccLc;
  if(fhistMCSpectrumAccLcFromSc) delete fhistMCSpectrumAccLcFromSc;
  if(fhistMCSpectrumAccSc) delete fhistMCSpectrumAccSc;
  if(fhSparseAnalysisSigma) delete fhSparseAnalysisSigma;
  
}
//_________________________________________________
void AliAnalysisTaskSESigmacTopK0Spi::Init() {
  //
  /// Initialization
  //

  fIsEventSelected = kFALSE;

  if (fDebug > 1) AliInfo("Init");

  fListCuts = new TList();
  fListCuts->SetOwner();
  fListCuts->Add(new AliRDHFCutsLctoV0(*fAnalCuts));
  PostData(3, fListCuts);

  // Save the weight functions or histograms
  fListWeight = new TList();
  fListWeight->SetOwner();
  fListWeight->Add(fHistoMCNch);
  PostData(7, fListWeight);
  
  if (fUseMCInfo && (fKeepingOnlyHIJINGBkg || fKeepingOnlyPYTHIABkg)) fUtils = new AliVertexingHFUtils();

  if(!fESDtrackCutsSoftPion){
    fESDtrackCutsSoftPion = new AliESDtrackCuts("AliESDtrackCuts", "default");
    fESDtrackCutsSoftPion->SetMinNClustersITS(3);
    fESDtrackCutsSoftPion->SetMaxDCAToVertexXY(0.065);
    fESDtrackCutsSoftPion->SetPtRange(0.05, 1.e10);
    fESDtrackCutsSoftPion->SetMaxDCAToVertexZ(0.15);
    fESDtrackCutsSoftPion->SetEtaRange(-0.9, +0.9);    
  } 
  
  return;
}

//________________________________________ terminate ___________________________
void AliAnalysisTaskSESigmacTopK0Spi::Terminate(Option_t*)
{
  /// The Terminate() function is the last function to be called during
  /// a query. It always runs on the client, it can be used to present
  /// the results graphically or save the results to file.

  AliInfo("Terminate");
  AliAnalysisTaskSE::Terminate();

  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutput) {
    AliError("fOutput not available");
    return;
  }

  if(fHistoMCLcK0SpGen) {
    AliInfo(Form("At MC level, %f Lc --> K0S + p were found", fHistoMCLcK0SpGen->GetEntries()));
  } else {
    AliInfo("fHistoMCLcK0SpGen not available");
  }
  if(fHistoMCLcK0SpGenAcc) {
    AliInfo(Form("At MC level, %f Lc --> K0S + p were found in the acceptance", fHistoMCLcK0SpGenAcc->GetEntries()));
  } else {
    AliInfo("fHistoMCLcK0SpGenAcc not available");
  }
  if(fVariablesTreeSgn) {
    AliInfo(Form("At Reco level, %lld Lc --> K0S + p were found", fVariablesTreeSgn->GetEntries()));
  } else {
    AliInfo("fVariablesTreeSgn not available");
  }

  fOutputSparse = dynamic_cast<TList*> (GetOutputData(6));
  if (!fOutputSparse) {
    AliError("fOutputSparse not available");
    return;
  }

  return;
}

//___________________________________________________________________________
void AliAnalysisTaskSESigmacTopK0Spi::UserCreateOutputObjects() {

  AliInfo(Form("CreateOutputObjects of task %s\n", GetName()));

  fOutput = new TList();
  fOutput->SetOwner();
  fOutput->SetName("listTrees");

  // Output slot 1: list of 2 trees (Sgn + Bkg) of the candidate variables
  const char* nameoutput = GetOutputSlot(1)->GetContainer()->GetName();
  fVariablesTreeSgn = new TTree(Form("%s_Sgn", nameoutput), "Candidates variables tree, Signal");
  fVariablesTreeBkg = new TTree(Form("%s_Bkg", nameoutput), "Candidates variables tree, Background");

  Int_t nVar; 
  if (fUseMCInfo)  nVar = 54; //"full" tree if MC
  else nVar = 38; //"reduced" tree if data
  
  fCandidateVariables = new Float_t [nVar];
  TString * fCandidateVariableNames = new TString[nVar];
  
  if (fUseMCInfo) { // "full tree" for MC
    fCandidateVariableNames[0] = "massLc2K0Sp";
    fCandidateVariableNames[1] = "massLc2Lambdapi";
    fCandidateVariableNames[2] = "massK0S";
    fCandidateVariableNames[3] = "massLambda";
    fCandidateVariableNames[4] = "massLambdaBar";
    fCandidateVariableNames[5] = "cosPAK0S";
    fCandidateVariableNames[6] = "dcaV0";
    fCandidateVariableNames[7] = "tImpParBach";
    fCandidateVariableNames[8] = "tImpParV0";
    fCandidateVariableNames[9] = "nSigmaTPCpr";
    fCandidateVariableNames[10] = "nSigmaTOFpr";
    fCandidateVariableNames[11] = "bachelorPt";
    fCandidateVariableNames[12] = "V0positivePt";
    fCandidateVariableNames[13] = "V0negativePt";
    fCandidateVariableNames[14] = "dcaV0pos";
    fCandidateVariableNames[15] = "dcaV0neg";
    fCandidateVariableNames[16] = "v0Pt";
    fCandidateVariableNames[17] = "massGamma";
    fCandidateVariableNames[18] = "LcPt";
    fCandidateVariableNames[19] = "combinedProtonProb";
    fCandidateVariableNames[20] = "LcEta";
    fCandidateVariableNames[21] = "V0positiveEta";
    fCandidateVariableNames[22] = "V0negativeEta";
    fCandidateVariableNames[23] = "TPCProtonProb";
    fCandidateVariableNames[24] = "TOFProtonProb";
    fCandidateVariableNames[25] = "bachelorEta";
    fCandidateVariableNames[26] = "LcP";
    fCandidateVariableNames[27] = "bachelorP";
    fCandidateVariableNames[28] = "v0P";
    fCandidateVariableNames[29] = "V0positiveP";
    fCandidateVariableNames[30] = "V0negativeP";
    fCandidateVariableNames[31] = "v0Eta";
    fCandidateVariableNames[32] = "LcPtMC";
    fCandidateVariableNames[33] = "DecayLengthK0S";
    fCandidateVariableNames[34] = "bachCode";
    fCandidateVariableNames[35] = "k0SCode";
    fCandidateVariableNames[36] = "alphaArm";
    fCandidateVariableNames[37] = "ptArm";
    fCandidateVariableNames[38] = "CosThetaStar";
    fCandidateVariableNames[39] = "weightPtFlat";
    //fCandidateVariableNames[40] = "weightFONLL5overLHC13d3";
    //fCandidateVariableNames[41] = "weightFONLL5overLHC13d3Lc";
    //fCandidateVariableNames[40] = "CosPALc";
    fCandidateVariableNames[40] = "SigmacPt";
    fCandidateVariableNames[41] = "CosThetaStarSoftPi";
    fCandidateVariableNames[42] = "weightNch";
    fCandidateVariableNames[43] = "NtrkRaw";
    fCandidateVariableNames[44] = "NtrkCorr";
    fCandidateVariableNames[45] = "signd0";
    fCandidateVariableNames[46] = "centrality";
    fCandidateVariableNames[47] = "NtrkAll";
    fCandidateVariableNames[48] = "origin";
    fCandidateVariableNames[49] = "nSigmaTPCpi";
    fCandidateVariableNames[50] = "nSigmaTPCka";
    //fCandidateVariableNames[51] = "bachTPCmom";
    fCandidateVariableNames[51] = "deltaM";
    fCandidateVariableNames[52] = "ptArmLc";
    fCandidateVariableNames[53] = "alphaArmLc";
  }
  else {   // "light mode"
    fCandidateVariableNames[0] = "massLc2K0Sp";
    fCandidateVariableNames[1] = "alphaArm";
    fCandidateVariableNames[2] = "massK0S";
    fCandidateVariableNames[3] = "massLambda";
    fCandidateVariableNames[4] = "massLambdaBar";
    fCandidateVariableNames[5] = "cosPAK0S";
    fCandidateVariableNames[6] = "dcaV0";
    fCandidateVariableNames[7] = "tImpParBach";
    fCandidateVariableNames[8] = "tImpParV0";
    fCandidateVariableNames[9] = "nSigmaTPCpr";
    fCandidateVariableNames[10] = "nSigmaTOFpr";
    fCandidateVariableNames[11] = "bachelorPt";
    fCandidateVariableNames[12] = "V0positivePt";
    fCandidateVariableNames[13] = "V0negativePt";
    fCandidateVariableNames[14] = "dcaV0pos";
    fCandidateVariableNames[15] = "dcaV0neg";
    fCandidateVariableNames[16] = "v0Pt";
    //fCandidateVariableNames[17] = "bachTPCmom";
    fCandidateVariableNames[17] = "SigmacPt";
    fCandidateVariableNames[18] = "LcPt";
    fCandidateVariableNames[19] = "combinedProtonProb";
    fCandidateVariableNames[20] = "V0positiveEta";
    //fCandidateVariableNames[21] = "bachelorP"; // we replaced the V0negativeEta with the bachelor P as this is more useful (for PID) while the V0 daughters' eta we don't use... And are practically the same (positive and negative)
    fCandidateVariableNames[21] = "ptArmLc";
    fCandidateVariableNames[22] = "bachelorEta";
    fCandidateVariableNames[23] = "v0P";
    fCandidateVariableNames[24] = "DecayLengthK0S";
    fCandidateVariableNames[25] = "nSigmaTPCpi";
    fCandidateVariableNames[26] = "nSigmaTPCka";
    fCandidateVariableNames[27] = "NtrkRaw";
    fCandidateVariableNames[28] = "NtrkCorr";
    fCandidateVariableNames[29] = "CosThetaStar";
    fCandidateVariableNames[30] = "signd0";        
    fCandidateVariableNames[31] = "centrality"; 
    fCandidateVariableNames[32] = "NtrkAll";
    fCandidateVariableNames[33] = "origin";
    fCandidateVariableNames[34] = "ptArm";
    fCandidateVariableNames[35] = "deltaM";
    //fCandidateVariableNames[36] = "CosPALc";
    fCandidateVariableNames[36] = "CosThetaStarSoftPi";
    fCandidateVariableNames[37] = "alphaArmLc";
  }
  
  for(Int_t ivar=0; ivar < nVar; ivar++){
    fVariablesTreeSgn->Branch(fCandidateVariableNames[ivar].Data(), &fCandidateVariables[ivar], Form("%s/f",fCandidateVariableNames[ivar].Data()));
    fVariablesTreeBkg->Branch(fCandidateVariableNames[ivar].Data(), &fCandidateVariables[ivar], Form("%s/f",fCandidateVariableNames[ivar].Data()));
  }
  
  fHistoCentrality = new TH1F("fHistoCentrality", "fHistoCentrality", 100, 0., 100.);

  fHistoEvents = new TH1F("fHistoEvents", "fHistoEvents", 6, -0.5, 5.5);
  TString labelEv[6] = {"RejectedDeltaMismatch", "AcceptedDeltaMismatch", "NotSelected", "TimeStampCut", "Selected", "AcceptedMultCut"};
  for (Int_t ibin = 1; ibin <= fHistoEvents->GetNbinsX(); ibin++){
    fHistoEvents->GetXaxis()->SetBinLabel(ibin, labelEv[ibin-1].Data());
  }

  fHistoTracklets_1 = new TH1F("fHistoTracklets_1", "fHistoTracklets_1", 1000, 0, 5000);
  fHistoTracklets_1_cent = new TH2F("fHistoTracklets_1_cent", "fHistoTracklets_1_cent; centrality; SPD tracklets [-1, 1]", 100, 0., 100., 1000, 0, 5000);
  fHistoTracklets_All = new TH1F("fHistoTracklets_All", "fHistoTracklets_All", 1000, 0, 5000);
  fHistoTracklets_All_cent = new TH2F("fHistoTracklets_All_cent", "fHistoTracklets_All_cent; centrality; SPD tracklets [-999, 999]", 100, 0., 100., 1000, 0, 5000);

  fHistoLc = new TH1F("fHistoLc", "fHistoLc", 2, -0.5, 1.5);

  fHistoLcOnTheFly = new TH1F("fHistoLcOnTheFly", "fHistoLcOnTheFly", 4, -0.5, 3.5);
  TString labelOnTheFly[4] = {"OnTheFly-Bkg", "OfflineBkg", "OnTheFly-Sgn", "OfflineSgn"};
  for (Int_t ibin = 1; ibin <= fHistoLcOnTheFly->GetNbinsX(); ibin++){
    fHistoLcOnTheFly->GetXaxis()->SetBinLabel(ibin, labelOnTheFly[ibin-1].Data());
  }

  fHistoLcBeforeCuts = new TH1F("fHistoLcBeforeCuts", "fHistoLcBeforeCuts", 2, -0.5, 1.5);
  TString labelBeforeCuts[2] = {"Bkg", "Sgn"};
  for (Int_t ibin = 1; ibin <= fHistoLcBeforeCuts->GetNbinsX(); ibin++){
    fHistoLcBeforeCuts->GetXaxis()->SetBinLabel(ibin, labelBeforeCuts[ibin-1].Data());
    fHistoLc->GetXaxis()->SetBinLabel(ibin, labelBeforeCuts[ibin-1].Data());
  }

  fHistoFiducialAcceptance = new TH1F("fHistoFiducialAcceptance", "fHistoFiducialAcceptance", 4, -0.5, 3.5);
  TString labelAcc[4] = {"NotOk-Bkg", "Ok-Bkg", "NotOk-Sgn", "Ok-Sgn"};
  for (Int_t ibin = 1; ibin <= fHistoFiducialAcceptance->GetNbinsX(); ibin++){
    fHistoFiducialAcceptance->GetXaxis()->SetBinLabel(ibin, labelAcc[ibin-1].Data());
  }

  fHistoCodesSgn = new TH2F("fHistoCodesSgn", "fHistoCodes for Signal; bachelor; K0S", 7, -1.5, 5.5, 9, -1.5, 7.5);
  fHistoCodesBkg = new TH2F("fHistoCodesBkg", "fHistoCodes for Background; bachelor; K0S", 7, -1.5, 5.5, 9, -1.5, 7.5);

  TString labelx[7] = { "kBachInvalid", "kBachFake", "kBachNoProton", "kBachPrimary", "kBachNoLambdaMother",
			"kBachDifferentLambdaMother",	"kBachCorrectLambdaMother"};
  TString labely[9] = { "kK0SInvalid", "kK0SFake", "kK0SNoK0S", "kK0SWithoutMother", "kK0SNotFromK0",
			"kK0Primary", "kK0NoLambdaMother", "kK0DifferentLambdaMother", "kK0CorrectLambdaMother"};

  for (Int_t ibin = 1; ibin <= fHistoCodesSgn->GetNbinsX(); ibin++){
    fHistoCodesSgn->GetXaxis()->SetBinLabel(ibin, labelx[ibin-1].Data());
    fHistoCodesBkg->GetXaxis()->SetBinLabel(ibin, labelx[ibin-1].Data());
  }
  for (Int_t ibin = 1; ibin <= fHistoCodesSgn->GetNbinsY(); ibin++){
    fHistoCodesSgn->GetYaxis()->SetBinLabel(ibin, labely[ibin-1].Data());
    fHistoCodesBkg->GetYaxis()->SetBinLabel(ibin, labely[ibin-1].Data());
  }

  fHistoBackground = new TH1F("fHistoBackground", "fHistoBackground", 4, -0.5, 3.5);
  TString labelBkg[4] = {"Injected", "Non-injected", "Non-PYTHIA", "PYTHIA"};
  for (Int_t ibin = 1; ibin <= fHistoBackground->GetNbinsX(); ibin++){
    fHistoBackground->GetXaxis()->SetBinLabel(ibin, labelBkg[ibin-1].Data());
  }

  const Float_t ptbins[15] = {0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 12., 17., 25., 35.};

  fHistoMCLcK0SpGen = new TH1F("fHistoMCLcK0SpGen", "fHistoMCLcK0SpGen", 14, ptbins);
  fHistoMCLcK0SpGenAcc = new TH1F("fHistoMCLcK0SpGenAcc", "fHistoMCLcK0SpGenAcc", 14, ptbins);
  fHistoMCLcK0SpGenLimAcc = new TH1F("fHistoMCLcK0SpGenLimAcc", "fHistoMCLcK0SpGenLimAcc", 14, ptbins);

  fBDTHisto = new TH2D("fBDTHisto", "Lc inv mass vs bdt output; bdt; m_{inv}(pK^{0}_{S})[GeV/#it{c}^{2}]", 10000, -1, 1, 1000, 2.05, 2.55);
  fBDTHistoTMVA = new TH2D("fBDTHistoTMVA", "Lc inv mass vs bdt output; bdt; m_{inv}(pK^{0}_{S})[GeV/#it{c}^{2}]", 10000, -1, 1, 1000, 2.05, 2.55);
  if (fDebugHistograms) {    
    fBDTHistoVsMassK0S = new TH2D("fBDTHistoVsMassK0S", "K0S inv mass vs bdt output; bdt; m_{inv}(#pi^{+}#pi^{#minus})[GeV/#it{c}^{2}]", 1000, -1, 1, 1000, 0.485, 0.51);
    fBDTHistoVstImpParBach = new TH2D("fBDTHistoVstImpParBach", "d0 bachelor vs bdt output; bdt; d_{0, bachelor}[cm]", 1000, -1, 1, 100, -1, 1);
    fBDTHistoVstImpParV0 = new TH2D("fBDTHistoVstImpParV0", "d0 K0S vs bdt output; bdt; d_{0, V0}[cm]", 1000, -1, 1, 100, -1, 1);
    fBDTHistoVsBachelorPt = new TH2D("fBDTHistoVsBachelorPt", "bachelor pT vs bdt output; bdt; p_{T, bachelor}[GeV/#it{c}]", 1000, -1, 1, 100, 0, 20);
    fBDTHistoVsCombinedProtonProb = new TH2D("fBDTHistoVsCombinedProtonProb", "combined proton probability vs bdt output; bdt; Bayesian PID_{bachelor}", 1000, -1, 1, 100, 0, 1);
    fBDTHistoVsCtau = new TH2D("fBDTHistoVsCtau", "K0S ctau vs bdt output; bdt; c#tau_{V0}[cm]",  1000, -1, 1, 1000, 0, 100);
    fBDTHistoVsCosPAK0S = new TH2D("fBDTHistoVsCosPAK0S", "V0 cosine pointing angle vs bdt output; bdt; CosPAK^{0}_{S}", 1000, -1, 1, 100, 0.9, 1);
    fBDTHistoVsCosThetaStar = new TH2D("fBDTHistoVsCosThetaStar", "proton emission angle in pK0s pair rest frame vs bdt output; bdt; Cos#Theta*", 1000, -1, 1, 100, -1, 1);
    fBDTHistoVsSignd0 = new TH2D("fBDTHistoVsSignd0", "signed d0 bachelor vs bdt output; bdt; signd_{0, bachelor}[cm]", 1000, -1, 1, 100, -1, 1);
    fBDTHistoVsnSigmaTPCpr = new TH2D("fBDTHistoVsnSigmaTPCpr", "nSigmaTPCpr vs bdt output; bdt; n_{#sigma}^{TPC}_{pr}", 1000, -1, 1, 1000, -10, 10);
    fBDTHistoVsnSigmaTOFpr = new TH2D("fBDTHistoVsnSigmaTOFpr", "nSigmaTOFpr vs bdt output; bdt; n_{#sigma}^{TOF}_{pr}", 1000, -1, 1, 1000, -10, 10);
    fBDTHistoVsnSigmaTPCpi = new TH2D("fBDTHistoVsnSigmaTPCpi", "nSigmaTPCpi vs bdt output; bdt; n_{#sigma}^{TPC}_{pi}", 1000, -1, 1, 1000, -10, 10);
    fBDTHistoVsnSigmaTPCka = new TH2D("fBDTHistoVsnSigmaTPCka", "nSigmaTPCka vs bdt output; bdt; n_{#sigma}^{TPC}_{ka}", 1000, -1, 1, 1000, -10, 10);
    fBDTHistoVsBachelorP = new TH2D("fBDTHistoVsBachelorP", "bachelor p vs bdt output; bdt; p_{bachelor}[GeV/#it{c}]", 1000, -1, 1, 100, 0, 20);
    fBDTHistoVsBachelorTPCP = new TH2D("fBDTHistoVsBachelorTPCP", "bachelor TPC momentum vs bdt output; bdt; p_{TPC, bachelor}[GeV/#it{c}]", 1000, -1, 1, 100, 0, 20);
    fHistoNsigmaTPC = new TH2D("fHistoNsigmaTPC", "; #it{p} (GeV/#it{c}); n_{#sigma}^{TPC} (proton hypothesis)", 500, 0, 5, 1000, -5, 5);
    fHistoNsigmaTOF = new TH2D("fHistoNsigmaTOF", "; #it{p} (GeV/#it{c}); n_{#sigma}^{TOF} (proton hypothesis)", 500, 0, 5, 1000, -5, 5);
  }
  
  fOutput->Add(fHistoEvents);
  fOutput->Add(fHistoTracklets_1);
  fOutput->Add(fHistoTracklets_1_cent);
  fOutput->Add(fHistoTracklets_All);
  fOutput->Add(fHistoTracklets_All_cent);
  fOutput->Add(fHistoLc);
  fOutput->Add(fHistoLcOnTheFly);
  fOutput->Add(fHistoLcBeforeCuts);
  fOutput->Add(fHistoFiducialAcceptance);
  fOutput->Add(fHistoCodesSgn);
  fOutput->Add(fHistoCodesBkg);
  fOutput->Add(fHistoBackground);
  fOutput->Add(fHistoMCLcK0SpGen);
  fOutput->Add(fHistoMCLcK0SpGenAcc);
  fOutput->Add(fHistoMCLcK0SpGenLimAcc);
  fOutput->Add(fHistoCentrality);
  fOutput->Add(fBDTHisto);
  fOutput->Add(fBDTHistoTMVA);
  if (fDebugHistograms) {    
    fOutput->Add(fBDTHistoVsMassK0S);
    fOutput->Add(fBDTHistoVstImpParBach);
    fOutput->Add(fBDTHistoVstImpParV0);
    fOutput->Add(fBDTHistoVsBachelorPt);
    fOutput->Add(fBDTHistoVsCombinedProtonProb);
    fOutput->Add(fBDTHistoVsCtau);
    fOutput->Add(fBDTHistoVsCosPAK0S);
    fOutput->Add(fBDTHistoVsCosThetaStar);
    fOutput->Add(fBDTHistoVsSignd0);
    fOutput->Add(fBDTHistoVsnSigmaTPCpr);
    fOutput->Add(fBDTHistoVsnSigmaTOFpr);
    fOutput->Add(fBDTHistoVsnSigmaTPCpi);
    fOutput->Add(fBDTHistoVsnSigmaTPCka);
    fOutput->Add(fBDTHistoVsBachelorP);
    fOutput->Add(fBDTHistoVsBachelorTPCP);
    fOutput->Add(fHistoNsigmaTPC);
    fOutput->Add(fHistoNsigmaTOF);
  }
  
  if(fUseMultCorrection){

    Int_t nMultBins      = 200;
    Float_t firstMultBin = -0.5;
    Float_t lastMultBin  = 199.5;
    const char *estimatorName="tracklets";
    if(fMultiplicityEstimator==kVZERO || fMultiplicityEstimator==kVZEROA || fMultiplicityEstimator==kVZEROEq || fMultiplicityEstimator==kVZEROAEq) {
      nMultBins = 400;
      lastMultBin = 799.5;
      estimatorName = "vzero";
    }

    fHistoNtrUnCorr = new TH1F("fHistoNtrUnCorr", Form("Uncorrected %s mulitplicity; %s; Entries",estimatorName,estimatorName),nMultBins, firstMultBin, lastMultBin);
    fHistoNtrCorr   = new TH1F("fHistoNtrCorr", Form("Corrected %s mulitplicity; %s; Entries",estimatorName,estimatorName),nMultBins, firstMultBin, lastMultBin);

    fHistoVzVsNtrUnCorr = new TH2F("fHistoVzVsNtrUnCorr", Form("VtxZ vs Uncorrected %s mulitplicity; VtxZ; %s; Entries",estimatorName,estimatorName), 300,-15.,15., nMultBins, firstMultBin, lastMultBin);
    fHistoVzVsNtrCorr   = new TH2F("fHistoVzVsNtrCorr", Form("VtxZ vs Corrected %s mulitplicity; VtxZ; %s; Entries",estimatorName,estimatorName), 300,-15.,15., nMultBins, firstMultBin, lastMultBin);

    fOutput->Add(fHistoNtrUnCorr);
    fOutput->Add(fHistoNtrCorr);
    fOutput->Add(fHistoVzVsNtrUnCorr);
    fOutput->Add(fHistoVzVsNtrCorr);
  }

  PostData(1, fOutput);
  PostData(4, fVariablesTreeSgn);
  PostData(5, fVariablesTreeBkg);
  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  fPIDResponse = inputHandler->GetPIDResponse();

  // Setting properties of PID
  fPIDCombined = new AliPIDCombined;
  fPIDCombined->SetDefaultTPCPriors();
  fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTPC+AliPIDResponse::kDetTOF);
  //fPIDCombined->SetPriorDistribution((AliPID::EParticleType)ispec,fPriors[ispec]);

  fCounter = new AliNormalizationCounter("NormalizationCounter");
  fCounter->Init();

  fCounterC = new AliNormalizationCounter("NormalizationCounterCorrMult");
  fCounterC->SetStudyMultiplicity(kTRUE,1.);
  fCounterC->Init();

  fListCounters = new TList();
  fListCounters->SetOwner();
  fListCounters->SetName("ListCounters");
  fListCounters->Add(fCounter);
  fListCounters->Add(fCounterC);
  
  PostData(2, fListCounters);
    
  // weight function from ratio of flat value (1/30) to pythia
  // use to normalise to flat distribution (should lead to flat pT distribution)
  fFuncWeightPythia = new TF1("funcWeightPythia","1./(30. *[0]*x/TMath::Power(1.+(TMath::Power((x/[1]),[3])),[2]))", 0.15, 30);
  fFuncWeightPythia->SetParameters(0.36609, 1.94635, 1.40463,2.5);
  
  // weight function from the ratio of the LHC13d3 MC
  // and FONLL calculations for pp data
  //fFuncWeightFONLL5overLHC13d3 = new TF1("funcWeightFONLL5overLHC13d3","([0]*x)/TMath::Power([2],(1+TMath::Power([3],x/[1])))+[4]*TMath::Exp([5]+[6]*x)+[7]*TMath::Exp([8]*x)", 0.15, 30.);
  //fFuncWeightFONLL5overLHC13d3->SetParameters(2.94999e+00, 3.47032e+00, 2.81278e+00, 2.5, 1.93370e-02, 3.86865e+00, -1.54113e-01, 8.86944e-02, 2.56267e-02);

  //fFuncWeightFONLL5overLHC13d3Lc = new TF1("funcWeightFONLL5overLHC13d3Lc","([0]*x)/TMath::Power([2],(1+TMath::Power([3],x/[1])))+[4]*TMath::Exp([5]+[6]*x)+[7]*TMath::Exp([8]*x)", 0.15, 20.);
  //fFuncWeightFONLL5overLHC13d3Lc->SetParameters(5.94428e+01, 1.63585e+01, 9.65555e+00, 6.71944e+00, 8.88338e-02, 2.40477e+00, -4.88649e-02, -6.78599e-01, -2.10951e-01);

  //PostData(6, fOutputKF);
 
  PostData(7, fListWeight);

  if (!fFillTree) {
    Printf("Booking methods");
    // creating the BDT and TMVA reader
    fVarsTMVA = new Float_t[fNVars];
    fVarsTMVASpectators = new Float_t[fNVarsSpectators];
    fReader = new TMVA::Reader( "!Color:!Silent" );
    std::vector<std::string> inputNamesVec;
    TObjArray *tokens = fNamesTMVAVar.Tokenize(",");
    for(Int_t i = 0; i < tokens->GetEntries(); i++){
      TString variable = ((TObjString*)(tokens->At(i)))->String();
      std::string tmpvar = variable.Data();
      inputNamesVec.push_back(tmpvar);
      if (fUseXmlWeightsFile || fUseXmlFileFromCVMFS) fReader->AddVariable(variable.Data(), &fVarsTMVA[i]);
    }      
    delete tokens;
    TObjArray *tokensSpectators = fNamesTMVAVarSpectators.Tokenize(",");
    for(Int_t i = 0; i < tokensSpectators->GetEntries(); i++){
      TString variable = ((TObjString*)(tokensSpectators->At(i)))->String();
      if (fUseXmlWeightsFile || fUseXmlFileFromCVMFS) fReader->AddSpectator(variable.Data(), &fVarsTMVASpectators[i]);
    }
    delete tokensSpectators;
    if (fUseWeightsLibrary) {
      void* lib = dlopen(fTMVAlibName.Data(), RTLD_NOW);
      void* p = dlsym(lib, Form("%s", fTMVAlibPtBin.Data()));
      IClassifierReader* (*maker1)(std::vector<std::string>&) = (IClassifierReader* (*)(std::vector<std::string>&)) p;
      fBDTReader = maker1(inputNamesVec);
    }
    
    if (fUseXmlWeightsFile) fReader->BookMVA("BDT method", fXmlWeightsFile);

    if (fUseXmlFileFromCVMFS){
       TString pathToFileCVMFS = AliDataFile::GetFileName(fXmlFileFromCVMFS.Data());
       if (pathToFileCVMFS.IsNull()){
          AliFatal("Cannot access data files from CVMFS");
       }
       fReader->BookMVA("BDT method", pathToFileCVMFS); 
    }
  }

  ftrackArraySelSoftPi = new TArrayI(10000);

  fhistMCSpectrumAccLc = new TH3F("fhistMCSpectrumAccLc", "fhistMCSpectrumAccLc", 250, 0, 50, 20, -0.5, 19.5, 2, 3.5, 5.5); // 
  fhistMCSpectrumAccSc = new TH3F("fhistMCSpectrumAccSc", "fhistMCSpectrumAccSc", 250, 0, 50, 20, -0.5, 19.5, 2, 3.5, 5.5); // 

  const Int_t nbinsAccLcFromSc = 6;
  Int_t binsAccLcFromSc[nbinsAccLcFromSc] = {250, 20, 2, 20, 250, 40};
  Double_t lowedgesAccLcFromSc[nbinsAccLcFromSc] = {0, -0.5, 3.5, -1, 0, -2};
  Double_t upedgesAccLcFromSc[nbinsAccLcFromSc] = {50, 19.5, 5.5, 1, 50, 2};
  fhistMCSpectrumAccLcFromSc = new THnSparseF("fhistMCSpectrumAccLcFromSc", "fhistMCSpectrumAccLcFromSc; ptLc; codeLc; Qorigin; yLc; ptSc; ySc", nbinsAccLcFromSc, binsAccLcFromSc, lowedgesAccLcFromSc, upedgesAccLcFromSc); // 
  
  Int_t nbinsSparseSigma[9] = {25, 400, 400, 20, 25, 2, 1, 2000, 2};
  Double_t lowEdgesSigma[9] = {0, 0.130, 2.100, -1, 0, 3.5, 0.5, -1, 0};
  Double_t upEdgesSigma[9] = {25, 0.330, 2.500, 1, 25, 5.5, 1.5, 1, 2};
  if(!fFillTree)  fhSparseAnalysisSigma = new THnSparseF("fhSparseAnalysisSigma", "fhSparseAnalysis; pt; deltamass; LcMass; CosThetaStarSoftPion; ptsigmac; checkorigin; isRotated; bdtresp; softPiITSrefit", 9, nbinsSparseSigma, lowEdgesSigma, upEdgesSigma);

  fOutputSparse = new TList();
  fOutputSparse->SetOwner();
  fOutputSparse->SetName("listOutputSparse");

  fOutputSparse->Add(fhistMCSpectrumAccLc);
  fOutputSparse->Add(fhistMCSpectrumAccSc);
  fOutputSparse->Add(fhistMCSpectrumAccLcFromSc);
  if(fhSparseAnalysisSigma)  fOutputSparse->Add(fhSparseAnalysisSigma);
  
  PostData(6, fOutputSparse);
  
  return;
}

//_________________________________________________
void AliAnalysisTaskSESigmacTopK0Spi::UserExec(Option_t *)
{
  /// user exec
  if (!fInputEvent) {
    AliError("NO EVENT FOUND!");
    return;
  }

  fCurrentEvent++;
  AliDebug(2, Form("Processing event = %d", fCurrentEvent));
  AliAODEvent* aodEvent = dynamic_cast<AliAODEvent*>(fInputEvent);

  if(fAODProtection >= 0){
    //   Protection against different number of events in the AOD and deltaAOD
    //   In case of discrepancy the event is rejected.
    Int_t matchingAODdeltaAODlevel = AliRDHFCuts::CheckMatchingAODdeltaAODevents();
    if (matchingAODdeltaAODlevel < 0 || (matchingAODdeltaAODlevel == 0 && fAODProtection == 1)) {
      // AOD/deltaAOD trees have different number of entries || TProcessID do not match while it was required
      fHistoEvents->Fill(0);
      return;
    }
    fHistoEvents->Fill(1);
  }
  
  TClonesArray *arrayLctopKos=0;

  if (!aodEvent && AODEvent() && IsStandardAOD()) {
    // In case there is an AOD handler writing a standard AOD, use the AOD
    // event in memory rather than the input (ESD) event.
    aodEvent = dynamic_cast<AliAODEvent*> (AODEvent());
    // in this case the braches in the deltaAOD (AliAOD.VertexingHF.root)
    // have to taken from the AOD event hold by the AliAODExtension
    AliAODHandler* aodHandler = (AliAODHandler*)
      ((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());

    if (aodHandler->GetExtensions()) {
      AliAODExtension *ext = (AliAODExtension*)aodHandler->GetExtensions()->FindObject("AliAOD.VertexingHF.root");
      AliAODEvent *aodFromExt = ext->GetAOD();
      arrayLctopKos = (TClonesArray*)aodFromExt->GetList()->FindObject("CascadesHF");
    }
  } else {
    arrayLctopKos = (TClonesArray*)aodEvent->GetList()->FindObject("CascadesHF");
  }
  
  if (!fUseMCInfo) {
    fAnalCuts->SetTriggerClass("");
    fAnalCuts->SetTriggerMask(fTriggerMask);
  }
  
  Int_t runnumber = aodEvent->GetRunNumber();
  if (aodEvent->GetTriggerMask() == 0 && (runnumber >= 195344 && runnumber <= 195677)){
    AliDebug(3, "Event rejected because of null trigger mask");
    return;
  }

  fCounter->StoreEvent(aodEvent, fAnalCuts, fUseMCInfo);

  // mc analysis
  TClonesArray *mcArray = 0;
  AliAODMCHeader *mcHeader = 0;

  // multiplicity definition with tracklets
  fNTracklets_1 = static_cast<Int_t>(AliVertexingHFUtils::GetNumberOfTrackletsInEtaRange(aodEvent, -1., 1.));
  fNTracklets_All = static_cast<Int_t>(AliVertexingHFUtils::GetNumberOfTrackletsInEtaRange(aodEvent, -999., 999.));
  if (fUseMCInfo) {
    // MC array need for matching
    mcArray = dynamic_cast<TClonesArray*>(aodEvent->FindListObject(AliAODMCParticle::StdBranchName()));
    if (!mcArray) {
      AliError("Could not find Monte-Carlo in AOD");
      return;
    }
    // load MC header
    mcHeader = (AliAODMCHeader*)aodEvent->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    if (!mcHeader) {
      AliError("AliAnalysisTaskSESigmacTopK0Spi::UserExec: MC header branch not found!\n");
      return;
    }
    
    Double_t zMCVertex = mcHeader->GetVtxZ();
    if (TMath::Abs(zMCVertex) > fAnalCuts->GetMaxVtxZ()){
      AliDebug(3, Form("z coordinate of MC vertex = %f, it was required to be within [-%f, +%f], skipping event", zMCVertex, fAnalCuts->GetMaxVtxZ(), fAnalCuts->GetMaxVtxZ()));
      return;
    }
    
    FillMCHisto(mcArray);
    LoopOverGenParticles(mcArray);

  }
  
  //centrality
  fCentrality = fAnalCuts->GetCentrality(aodEvent, AliRDHFCuts::kCentV0M);
  
  // AOD primary vertex
  fVtx1 = (AliAODVertex*)aodEvent->GetPrimaryVertex();
  if (!fVtx1) {
    AliDebug(2, "No primary vertex found, returning");
    return;
  }
  if (fVtx1->GetNContributors()<1) {
    AliDebug(2, "Number of contributors to vertex < 1, returning");
    return;
  }
  
  Int_t fNTracklets_03=AliVertexingHFUtils::GetNumberOfTrackletsInEtaRange(aodEvent,-0.3,0.3);
  Int_t fNTracklets_05=AliVertexingHFUtils::GetNumberOfTrackletsInEtaRange(aodEvent,-0.5,0.5);
  Int_t fNTracklets_16=AliVertexingHFUtils::GetNumberOfTrackletsInEtaRange(aodEvent,-1.6,1.6);

  Int_t vzeroMult=0, vzeroMultA=0, vzeroMultC=0;
  Int_t vzeroMultEq=0, vzeroMultAEq=0, vzeroMultCEq=0;
  AliAODVZERO *vzeroAOD = (AliAODVZERO*)aodEvent->GetVZEROData();
  if(vzeroAOD) {
    vzeroMultA = static_cast<Int_t>(vzeroAOD->GetMTotV0A());
    vzeroMultC = static_cast<Int_t>(vzeroAOD->GetMTotV0C());
    vzeroMult = vzeroMultA + vzeroMultC;
    vzeroMultAEq = static_cast<Int_t>(AliVertexingHFUtils::GetVZEROAEqualizedMultiplicity(aodEvent));
    vzeroMultCEq = static_cast<Int_t>(AliVertexingHFUtils::GetVZEROCEqualizedMultiplicity(aodEvent));
    vzeroMultEq = vzeroMultAEq + vzeroMultCEq;
  }

  Int_t countMult = fNTracklets_1;
  if(fMultiplicityEstimator==kNtrk03) { countMult = fNTracklets_03; }
  else if(fMultiplicityEstimator==kNtrk05) { countMult = fNTracklets_05; }
  else if(fMultiplicityEstimator==kNtrk10to16) { countMult = fNTracklets_16  - fNTracklets_1; }
  else if(fMultiplicityEstimator==kVZERO) { countMult = vzeroMult; }
  else if(fMultiplicityEstimator==kVZEROA) { countMult = vzeroMultA; }
  else if(fMultiplicityEstimator==kVZEROEq) { countMult = vzeroMultEq; }
  else if(fMultiplicityEstimator==kVZEROAEq) { countMult = vzeroMultAEq; }

  // Double_t countTreta1corr = fNTracklets_1;
  Double_t countCorr=countMult;

  if(fUseMultCorrection){
    // In case of VZERO multiplicity, consider the zvtx correction flag
    //  fDoVZER0ParamVertexCorr: 0= none, 1= usual d2h, 2=AliESDUtils
    Bool_t isDataDrivenZvtxCorr=kTRUE;
    Int_t vzeroMultACorr=vzeroMultA, vzeroMultCCorr=vzeroMultC, vzeroMultCorr=vzeroMult;
    Int_t vzeroMultAEqCorr=vzeroMultAEq, vzeroMultCEqCorr=vzeroMultCEq, vzeroMultEqCorr=vzeroMultEq;

    if( (fMultiplicityEstimator==kVZERO) || (fMultiplicityEstimator==kVZEROA) ||
        (fMultiplicityEstimator==kVZEROEq) || (fMultiplicityEstimator==kVZEROAEq)  )
    {
      if(fDoVZER0ParamVertexCorr==0){
        // do not correct
        isDataDrivenZvtxCorr=kFALSE;
      }
      else if (fDoVZER0ParamVertexCorr==2){
        // use AliESDUtils correction
        Float_t zvtx = fVtx1->GetZ();
        isDataDrivenZvtxCorr=kFALSE;
        vzeroMultACorr = static_cast<Int_t>(AliESDUtils::GetCorrV0A(vzeroMultA,zvtx));
        vzeroMultCCorr = static_cast<Int_t>(AliESDUtils::GetCorrV0C(vzeroMultC,zvtx));
        vzeroMultCorr = vzeroMultACorr + vzeroMultCCorr;
        vzeroMultAEqCorr = static_cast<Int_t>(AliESDUtils::GetCorrV0A(vzeroMultAEq,zvtx));
        vzeroMultCEqCorr =static_cast<Int_t>( AliESDUtils::GetCorrV0C(vzeroMultCEq,zvtx));
        vzeroMultEqCorr = vzeroMultAEqCorr + vzeroMultCEqCorr;
        if(fMultiplicityEstimator==kVZERO) { countCorr = vzeroMultCorr; }
        else if(fMultiplicityEstimator==kVZEROA) { countCorr = vzeroMultACorr; }
        else if(fMultiplicityEstimator==kVZEROEq) { countCorr = vzeroMultEqCorr; }
        else if(fMultiplicityEstimator==kVZEROAEq) { countCorr = vzeroMultAEqCorr; }
      }
    }
    if(isDataDrivenZvtxCorr){
      TProfile* estimatorAvg = GetEstimatorHistogram(aodEvent);
      if(estimatorAvg){
        // countTreta1corr = static_cast<Int_t>(AliVertexingHFUtils::GetCorrectedNtracklets(estimatorAvg,fNTracklets_1,vtx1->GetZ(),fRefMult));
        countCorr       = static_cast<Int_t>(AliVertexingHFUtils::GetCorrectedNtracklets(estimatorAvg,countMult,fVtx1->GetZ(),fRefMult));
      }
    }
  }

  fCounterC->StoreEvent(aodEvent, fAnalCuts, fUseMCInfo,countCorr);
  
  Bool_t isSelectedMultCut = kTRUE;
  if(fUseMultiplicityCut){
    // Multiplicity cut used
    if(countCorr >= fMultiplicityCutMin && countCorr < fMultiplicityCutMax){
      // Within multiplicity range
      fHistoEvents->Fill(5);
    }
    else{
      // Outside multiplicity range
      isSelectedMultCut = kFALSE;
    }      
  }

  fIsEventSelected = fAnalCuts->IsEventSelected(aodEvent);

  if ( !fIsEventSelected || !isSelectedMultCut) {
    fHistoEvents->Fill(2);
    return; // don't take into account not selected events
  }

  // check on the timestamp
  AliVHeader* h = aodEvent->GetHeader();
  UInt_t timestamp = h->GetTimeStamp();
  //Printf("timestamp = %d, cut = %u", timestamp, fTimestampCut);
  if (fTimestampCut != 0) {
    //Printf("timestamp = %d, cut = %u", timestamp, fTimestampCut);
    if (timestamp > fTimestampCut) {
      fHistoEvents->Fill(3);
      return;
    }
  }

  fHistoEvents->Fill(4);

  fHistoTracklets_1->Fill(fNTracklets_1);
  fHistoTracklets_All->Fill(fNTracklets_All);
  fHistoTracklets_1_cent->Fill(fCentrality, fNTracklets_1);
  fHistoTracklets_All_cent->Fill(fCentrality, fNTracklets_All);

  fHistoCentrality->Fill(fCentrality);
    
  if(fUseMultCorrection){
    fHistoNtrUnCorr->Fill(countMult);
    fHistoNtrCorr->Fill(countCorr);
    fHistoVzVsNtrUnCorr->Fill(fVtx1->GetZ(),countMult);
    fHistoVzVsNtrCorr->Fill(fVtx1->GetZ(),countCorr);
  }
    
  //Setting magnetic field for KF vertexing
  fBField = aodEvent->GetMagneticField();
  AliKFParticle::SetField(fBField);


  PrepareTracks(aodEvent); //needed for SigmaC loop
  
  Int_t nSelectedAnal = 0;
  if (fIsK0sAnalysis) {
    MakeAnalysisForLc2prK0S(aodEvent, arrayLctopKos, mcArray,
			    nSelectedAnal, fAnalCuts, mcHeader);
  }
  fCounter->StoreCandidates(aodEvent, nSelectedAnal, kFALSE);
  
  Double_t nchWeight = 1.;
  if (fNTracklets_1 > 0) {
    
    if(!fHistoMCNch) AliInfo("Input histos to evaluate Nch weights missing"); 
    if(fHistoMCNch) nchWeight *= fHistoMCNch->GetBinContent(fHistoMCNch->FindBin(fNTracklets_1));
  }


  PostData(1, fOutput);
  PostData(2, fListCounters);
  PostData(4, fVariablesTreeSgn);
  PostData(5, fVariablesTreeBkg);
  PostData(6, fOutputSparse);
  PostData(7, fListWeight);
}
//-------------------------------------------------------------------------------
void AliAnalysisTaskSESigmacTopK0Spi::FillMCHisto(TClonesArray *mcArray){

  /// method to fill MC histo: how many Lc --> K0S + p are there at MC level
  
  for (Int_t iPart = 0; iPart < mcArray->GetEntriesFast(); iPart++) {
    AliAODMCParticle* mcPart = dynamic_cast<AliAODMCParticle*>(mcArray->At(iPart));
    if (!mcPart){
      AliError("Failed casting particle from MC array!, Skipping particle");
      continue;
    }
    Int_t pdg = mcPart->GetPdgCode();
    if (TMath::Abs(pdg) != 4122){
      AliDebug(2, Form("MC particle %d is not a Lc: its pdg code is %d", iPart, pdg));
      continue;
    }
    AliDebug(2, Form("Step 0 ok: MC particle %d is a Lc: its pdg code is %d", iPart, pdg));
    Int_t labeldaugh0 = mcPart->GetDaughterLabel(0);
    Int_t labeldaugh1 = mcPart->GetDaughterLabel(1);
    if (labeldaugh0 <= 0 || labeldaugh1 <= 0){
      AliDebug(2, Form("The MC particle doesn't have correct daughters, skipping!!"));
      continue;
    }
    else if (labeldaugh1 - labeldaugh0 == 1){
      AliDebug(2, Form("Step 1 ok: The MC particle has correct daughters!!"));
      AliAODMCParticle* daugh0 = dynamic_cast<AliAODMCParticle*>(mcArray->At(labeldaugh0));
      AliAODMCParticle* daugh1 = dynamic_cast<AliAODMCParticle*>(mcArray->At(labeldaugh1));
      if(!daugh0 || !daugh1){
	AliDebug(2,"Particle daughters not properly retrieved!");
	return;
      }
      Int_t pdgCodeDaugh0 = TMath::Abs(daugh0->GetPdgCode());
      Int_t pdgCodeDaugh1 = TMath::Abs(daugh1->GetPdgCode());
      AliAODMCParticle* bachelorMC = daugh0;
      AliAODMCParticle* v0MC = daugh1;
      AliDebug(2, Form("pdgCodeDaugh0 = %d, pdgCodeDaugh1 = %d", pdgCodeDaugh0, pdgCodeDaugh1));
      if ((pdgCodeDaugh0 == 311 && pdgCodeDaugh1 == 2212) || (pdgCodeDaugh0 == 2212 && pdgCodeDaugh1 == 311)){
	// we are in the case of Lc --> K0 + p; now we have to check if the K0 decays in K0S, and if this goes in pi+pi-
	/// first, we set the bachelor and the v0: above we assumed first proton and second V0, but we could have to change it:
	if (pdgCodeDaugh0 == 311 && pdgCodeDaugh1 == 2212) {
	  bachelorMC = daugh1;
	  v0MC = daugh0;
	}
	AliDebug(2, Form("Number of Daughters of v0 = %d", v0MC->GetNDaughters()));
	if (v0MC->GetNDaughters() != 1) {
	  AliDebug(2, "The K0 does not decay in 1 body only! Impossible... Continuing...");
	  continue;
	}
	else { // So far: Lc --> K0 + p, K0 with 1 daughter
	  AliDebug(2, "Step 2 ok: The K0 does decay in 1 body only! ");
	  Int_t labelK0daugh = v0MC->GetDaughterLabel(0);
	  AliAODMCParticle* partK0S = dynamic_cast<AliAODMCParticle*>(mcArray->At(labelK0daugh));
	  if(!partK0S){
	    AliError("Error while casting particle! returning a NULL array");
	    continue;
	  }
	  else { // So far: Lc --> K0 + p, K0 with 1 daughter that we can access
	    if (partK0S->GetNDaughters() != 2 || TMath::Abs(partK0S->GetPdgCode() != 310)){
	      AliDebug(2, "The K0 daughter is not a K0S or does not decay in 2 bodies");
	      continue;
	    }
	    else { // So far: Lc --> K0 + p, K0 --> K0S, K0S in 2 bodies
	      AliDebug(2, "Step 3 ok: The K0 daughter is a K0S and does decay in 2 bodies");
	      Int_t labelK0Sdaugh0 = partK0S->GetDaughterLabel(0);
	      Int_t labelK0Sdaugh1 = partK0S->GetDaughterLabel(1);
	      AliAODMCParticle* daughK0S0 = dynamic_cast<AliAODMCParticle*>(mcArray->At(labelK0Sdaugh0));
	      AliAODMCParticle* daughK0S1 = dynamic_cast<AliAODMCParticle*>(mcArray->At(labelK0Sdaugh1));
	      if (!daughK0S0 || ! daughK0S1){
		AliDebug(2, "Could not access K0S daughters, continuing...");
		continue;
	      }
	      else { // So far: Lc --> K0 + p, K0 --> K0S, K0S in 2 bodies that we can access
		AliDebug(2, "Step 4 ok: Could access K0S daughters, continuing...");
		Int_t pdgK0Sdaugh0 = daughK0S0->GetPdgCode();
		Int_t pdgK0Sdaugh1 = daughK0S1->GetPdgCode();
		if (TMath::Abs(pdgK0Sdaugh0) != 211 || TMath::Abs(pdgK0Sdaugh1) != 211){
		  AliDebug(2, "The K0S does not decay in pi+pi-, continuing");
		  //AliInfo("The K0S does not decay in pi+pi-, continuing");
		}
		else { // Full chain: Lc --> K0 + p, K0 --> K0S, K0S --> pi+pi-
		  if (fAnalCuts->IsInFiducialAcceptance(mcPart->Pt(), mcPart->Y())) {
		    AliDebug(2, Form("----> Filling histo with pt = %f", mcPart->Pt()));
		    if(TMath::Abs(mcPart->Y()) < 0.5) fHistoMCLcK0SpGenLimAcc->Fill(mcPart->Pt());
		    //AliInfo(Form("\nparticle = %d, Filling MC Gen histo\n", iPart));
		    fHistoMCLcK0SpGen->Fill(mcPart->Pt());
		    if(!(TMath::Abs(bachelorMC->Eta()) > 0.9 || bachelorMC->Pt() < 0.1 ||
			 TMath::Abs(daughK0S0->Eta()) > 0.9 || daughK0S0->Pt() < 0.1 ||
			 TMath::Abs(daughK0S1->Eta()) > 0.9 || daughK0S1->Pt() < 0.1)) {
		      fHistoMCLcK0SpGenAcc->Fill(mcPart->Pt());
		    }
		  }
		  else {
		    AliDebug(2, "not in fiducial acceptance! Skipping");
		    continue;
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  } // closing loop over mcArray
  
  return;
  
}

//-------------------------------------------------------------------------------
void AliAnalysisTaskSESigmacTopK0Spi::MakeAnalysisForLc2prK0S(AliAODEvent *aodEvent,
							      TClonesArray *arrayLctopKos,
							      TClonesArray *mcArray,
							      Int_t &nSelectedAnal,
							      AliRDHFCutsLctoV0 *cutsAnal,
							      AliAODMCHeader* aodheader){

  /// Lc prong needed to MatchToMC method
  
  Int_t pdgCand = 4122;
  Int_t pdgDgLctoV0bachelor[2] = {2212, 310};
  Int_t pdgDgV0toDaughters[2] = {211, 211};

  // loop to search for candidates Lc->K0sp
  Int_t nCascades= arrayLctopKos->GetEntriesFast();

  // loop over cascades to search for candidates Lc->p+K0S

  Int_t mcLabel = -1;

  AliAnalysisVertexingHF *vHF = new AliAnalysisVertexingHF();
  for (Int_t iLctopK0s = 0; iLctopK0s < nCascades; iLctopK0s++) {

    // Lc candidates and K0s from Lc
    AliAODRecoCascadeHF* lcK0spr = dynamic_cast<AliAODRecoCascadeHF*>(arrayLctopKos->At(iLctopK0s));
    if (!lcK0spr) {
      AliDebug(2, Form("Cascade %d doens't exist, skipping",iLctopK0s));
      continue;
    }

    if (!(lcK0spr->CheckCascadeFlags())) {
      AliDebug(2, Form("Cascade %d is not flagged as Lc candidate",iLctopK0s));
      continue;
    }

    // use Preselect to filter out the tracks according to the pT
    if (cutsAnal->GetUsePreselect()){
      TObjArray arrTracks(2);
      for(Int_t ipr = 0; ipr < 2; ipr++){
	AliAODTrack *tr;
	if (ipr == 0) tr = vHF->GetProng(aodEvent, lcK0spr, ipr);
	else tr = (AliAODTrack*)(aodEvent->GetV0(lcK0spr->GetProngID(1)));
	arrTracks.AddAt(tr, ipr);
      }
      Int_t preSelectLc = cutsAnal->PreSelect(arrTracks);
      if (preSelectLc == 0) continue;
    }

    if(!vHF->FillRecoCasc(aodEvent, lcK0spr, kFALSE)){ //Fill the data members of the candidate only if they are empty.
      continue;
    }
    //if (!(vHF->RecoSecondaryVertexForCascades(aodEvent, lcK0spr))) continue;
    
    
    if (!lcK0spr->GetSecondaryVtx()) {
      AliInfo("No secondary vertex");
      continue;
    }

    if (lcK0spr->GetNDaughters()!=2) {
      AliDebug(2, Form("Cascade %d has not 2 daughters (nDaughters=%d)", iLctopK0s, lcK0spr->GetNDaughters()));
      continue;
    }

    AliAODv0 * v0part = dynamic_cast<AliAODv0*>(lcK0spr->Getv0());
    AliAODTrack * bachPart = dynamic_cast<AliAODTrack*>(lcK0spr->GetBachelor());
    if (!v0part || !bachPart) {
      AliDebug(2, Form("Cascade %d has no V0 or no bachelor object", iLctopK0s));
      continue;
    }

    if (!v0part->GetSecondaryVtx()) {
      AliDebug(2, Form("No secondary vertex for V0 by cascade %d", iLctopK0s));
      continue;
    }

    if (v0part->GetNDaughters()!=2) {
      AliDebug(2, Form("current V0 has not 2 daughters (onTheFly=%d, nDaughters=%d)", v0part->GetOnFlyStatus(), v0part->GetNDaughters()));
      continue;
    }

    AliAODTrack * v0Pos = dynamic_cast<AliAODTrack*>(lcK0spr->Getv0PositiveTrack());
    AliAODTrack * v0Neg = dynamic_cast<AliAODTrack*>(lcK0spr->Getv0NegativeTrack());
    if (!v0Neg || !v0Pos) {
      AliDebug(2, Form("V0 by cascade %d has no V0positive of V0negative object", iLctopK0s));
      continue;
    }

    if (v0Pos->Charge() == v0Neg->Charge()) {
      AliDebug(2, Form("V0 by cascade %d has daughters with the same sign: IMPOSSIBLE!", iLctopK0s));
      continue;
    }

    //Filling a control histogram with no cuts
    if (fUseMCInfo) {

      Int_t pdgCode=-2;
      
      // find associated MC particle for Lc -> p+K0 and K0S->pi+pi
      fmcLabelLc = lcK0spr->MatchToMC(pdgCand, pdgDgLctoV0bachelor[1], pdgDgLctoV0bachelor, pdgDgV0toDaughters, mcArray, kTRUE);
      //if (fmcLabelLc>=0) {
      if (fmcLabelLc != -1) {
	AliDebug(2, Form("----> cascade number %d (total cascade number = %d) is a Lc!", iLctopK0s, nCascades));

	AliAODMCParticle *partLc = dynamic_cast<AliAODMCParticle*>(mcArray->At(fmcLabelLc));
	if(partLc){
	  pdgCode = partLc->GetPdgCode();
	  if (pdgCode<0) AliDebug(2,Form(" ¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯ MClabel=%d ~~~~~~~~~~ pdgCode=%d", fmcLabelLc, pdgCode));
	  pdgCode = TMath::Abs(pdgCode);
	  fHistoLcBeforeCuts->Fill(1);
	}
      }
      else {
	fHistoLcBeforeCuts->Fill(0);
	pdgCode = -1;
      }

    }
    
    Int_t isLc = 0;
    
    if (fUseMCInfo) {

      Int_t pdgCode = -2;

      // find associated MC particle for Lc -> p+K0 and K0S->pi+pi
      mcLabel = lcK0spr->MatchToMC(pdgCand, pdgDgLctoV0bachelor[1], pdgDgLctoV0bachelor, pdgDgV0toDaughters, mcArray, kTRUE);
      if (mcLabel >= 0) {
	AliDebug(2, Form(" ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~cascade number %d (total cascade number = %d)", iLctopK0s, nCascades));

	AliAODMCParticle *partLc = dynamic_cast<AliAODMCParticle*>(mcArray->At(mcLabel));
	if(partLc){
	  pdgCode = partLc->GetPdgCode();
	  if (pdgCode < 0) AliDebug(2, Form(" ¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯ MClabel=%d ~~~~~~~~~~ pdgCode=%d", mcLabel, pdgCode));
	  pdgCode = TMath::Abs(pdgCode);
	  isLc = 1;
	  fHistoLc->Fill(1);
	}
      }
      else {
	fHistoLc->Fill(0);
	pdgCode = -1;
      }
    }
    AliDebug(2, Form("\n\n\n Analysing candidate %d\n", iLctopK0s));
    AliDebug(2, Form(">>>>>>>>>> Candidate is background, fFillOnlySgn = %d --> SKIPPING", fFillOnlySgn));
    if (!isLc) {
      if (fFillOnlySgn) { // if it is background, and we want only signal, we do not fill the tree
	continue;
      }
      else { // checking if we want to fill the background
	if (fKeepingOnlyHIJINGBkg){
	  // we have decided to fill the background only when the candidate has the daugthers that all come from HIJING underlying event!
	  Bool_t isCandidateInjected = fUtils->HasCascadeCandidateAnyDaughInjected(lcK0spr, aodheader, mcArray);
	  if (!isCandidateInjected){
	    AliDebug(2, "The candidate is from HIJING (i.e. not injected), keeping it to fill background");
	    fHistoBackground->Fill(1);
	  }
	  else {
	    AliDebug(2, "The candidate is NOT from HIJING, we skip it when filling background");
	    fHistoBackground->Fill(0);
	    continue;
	  }
	}
	else if (fKeepingOnlyPYTHIABkg){
	  // we have decided to fill the background only when the candidate has the daugthers that all come from PYTHIA underlying event!
	  AliAODTrack *bachelor = (AliAODTrack*)lcK0spr->GetBachelor();
	  AliAODTrack *v0pos = (AliAODTrack*)lcK0spr->Getv0PositiveTrack();
	  AliAODTrack *v0neg = (AliAODTrack*)lcK0spr->Getv0NegativeTrack();
	  if (!bachelor || !v0pos || !v0neg) {
	    AliDebug(2, "Cannot retrieve one of the tracks while checking origin, continuing");
	    continue;
	  }
	  else {
	    Int_t labelbachelor = TMath::Abs(bachelor->GetLabel());
	    Int_t labelv0pos = TMath::Abs(v0pos->GetLabel());
	    Int_t labelv0neg = TMath::Abs(v0neg->GetLabel());
	    AliAODMCParticle* MCbachelor =  (AliAODMCParticle*)mcArray->At(labelbachelor);
	    AliAODMCParticle* MCv0pos =  (AliAODMCParticle*)mcArray->At(labelv0pos);
	    AliAODMCParticle* MCv0neg =  (AliAODMCParticle*)mcArray->At(labelv0neg);
	    if (!MCbachelor || !MCv0pos || !MCv0neg) {
	      AliDebug(2, "Cannot retrieve MC particle for one of the tracks while checking origin, continuing");
	      continue;
	    }
	    else {
	      Int_t isBachelorFromPythia = fUtils->CheckOrigin(mcArray, MCbachelor, kTRUE);
	      Int_t isv0posFromPythia = fUtils->CheckOrigin(mcArray, MCv0pos, kTRUE);
	      Int_t isv0negFromPythia = fUtils->CheckOrigin(mcArray, MCv0neg, kTRUE);
	      if (isBachelorFromPythia != 0 && isv0posFromPythia != 0 && isv0negFromPythia != 0){
		AliDebug(2, "The candidate is from PYTHIA (i.e. all daughters originate from a quark), keeping it to fill background");
		fHistoBackground->Fill(2);
	      }
	      else {
		AliDebug(2, "The candidate is NOT from PYTHIA, we skip it when filling background");
		fHistoBackground->Fill(3);
		continue;
	      }
	    }
	  }
	}
      }
    }

    //FillLc2pK0Sspectrum(lcK0spr, isLc, nSelectedAnal, cutsAnal, mcArray, iLctopK0s);
    //FillLc2pK0Sspectrum(lcK0spr, isLc, nSelectedAnal, cutsAnal, mcArray, mcLabel);
    FillLc2pK0Sspectrum(lcK0spr, isLc, nSelectedAnal, cutsAnal, mcArray, mcLabel, aodEvent);
  }
  
  delete vHF;

  return;

}
//________________________________________________________________________
void AliAnalysisTaskSESigmacTopK0Spi::FillLc2pK0Sspectrum(AliAODRecoCascadeHF *part,
							  Int_t isLc,
							  Int_t &nSelectedAnal,
							  AliRDHFCutsLctoV0 *cutsAnal,
							  TClonesArray *mcArray, Int_t iLctopK0s, AliAODEvent *aod){
  //
  /// Fill histos for Lc -> K0S+proton
  //
  /*
    if (!part->GetOwnPrimaryVtx()) {
    //Printf("No primary vertex for Lc found!!");
    part->SetOwnPrimaryVtx(fVtx1);
    }
    else {
    //Printf("Yu-huuuu!!! primary vertex for Lc found!!");
    }
  */
  Double_t invmassLc = part->InvMassLctoK0sP();

  AliAODv0 * v0part = part->Getv0();
  Bool_t onFlyV0 = v0part->GetOnFlyStatus(); // on-the-flight V0s
  if (onFlyV0){ // on-the-fly V0
    if (isLc) { // Lc
      fHistoLcOnTheFly->Fill(2.);
    }
    else { // not Lc
      fHistoLcOnTheFly->Fill(0.);
    }
  }
  else { // offline V0
    if (isLc) { // Lc
      fHistoLcOnTheFly->Fill(3.);
    }
    else { // not Lc
      fHistoLcOnTheFly->Fill(1.);
    }
  }

  Double_t dcaV0 = v0part->GetDCA();
  Double_t invmassK0s = v0part->MassK0Short();

  if ( (cutsAnal->IsInFiducialAcceptance(part->Pt(), part->Y(4122))) ) {
    if (isLc) {
      fHistoFiducialAcceptance->Fill(3.);
    }
    else {
      fHistoFiducialAcceptance->Fill(1.);
    }
  }
  else {
    if (isLc) {
      fHistoFiducialAcceptance->Fill(2.);
    }
    else {
      fHistoFiducialAcceptance->Fill(0.);
    }
  }

  Int_t isInV0window = (((cutsAnal->IsSelectedSingleCut(part, AliRDHFCuts::kCandidate, 2, aod)) & (AliRDHFCutsLctoV0::kLcToK0Spr)) == (AliRDHFCutsLctoV0::kLcToK0Spr)); // cut on V0 invMass

  if (isInV0window == 0) {
    AliDebug(2, "No: The candidate has NOT passed the V0 window cuts!");
    if (isLc) AliDebug(2, "SIGNAL candidate rejected: V0 window cuts");
    return;
  }
  else AliDebug(2, "Yes: The candidate has passed the mass cuts!");  

  Bool_t isInCascadeWindow = (((cutsAnal->IsSelectedSingleCut(part, AliRDHFCuts::kCandidate, 0, aod)) & (AliRDHFCutsLctoV0::kLcToK0Spr)) == (AliRDHFCutsLctoV0::kLcToK0Spr)); // cut on Lc->p+K0S invMass

  if (!isInCascadeWindow) {
    AliDebug(2, "No: The candidate has NOT passed the cascade window cuts!");
    if (isLc) AliDebug(2, "SIGNAL candidate rejected: cascade window cuts");
    return;
  }
  else AliDebug(2, "Yes: The candidate has passed the cascade window cuts!");

  Bool_t isCandidateSelectedCuts = (((cutsAnal->IsSelected(part, AliRDHFCuts::kCandidate, aod)) & (AliRDHFCutsLctoV0::kLcToK0Spr)) == (AliRDHFCutsLctoV0::kLcToK0Spr)); // kinematic/topological cuts
  AliDebug(2, Form("recoAnalysisCuts = %d", cutsAnal->IsSelected(part, AliRDHFCuts::kCandidate, aod) & (AliRDHFCutsLctoV0::kLcToK0Spr)));
  if (!isCandidateSelectedCuts){
    AliDebug(2, "No: Analysis cuts kCandidate level NOT passed");
    if (isLc) AliDebug(2, "SIGNAL candidate rejected");
    return;
  }
  else {
    AliDebug(2, "Yes: Analysis cuts kCandidate level passed");
  }

  AliAODTrack *bachelor = (AliAODTrack*)part->GetBachelor();
  if (!bachelor) {
    AliDebug(2, Form("Very weird, the bachelor is not there... returning for this candidate"));
    return;
  }

  //Bool_t isBachelorID = (((cutsAnal->IsSelected(part,AliRDHFCuts::kPID, aod))&(AliRDHFCutsLctoV0::kLcToK0Spr))==(AliRDHFCutsLctoV0::kLcToK0Spr)); // ID x bachelor
  Double_t probTPCTOF[AliPID::kSPECIES] = {-1.};

  UInt_t detUsed = fPIDCombined->ComputeProbabilities(bachelor, fPIDResponse, probTPCTOF);
  AliDebug(2, Form("detUsed (TPCTOF case) = %d", detUsed));
  Double_t probProton = -1.;
  //  Double_t probPion = -1.;
  //  Double_t probKaon = -1.;
  if (detUsed == (UInt_t)fPIDCombined->GetDetectorMask() ) {
    AliDebug(2, Form("We have found the detector mask for TOF + TPC: probProton will be set to %f", probTPCTOF[AliPID::kProton]));
    probProton = probTPCTOF[AliPID::kProton];
    // probPion = probTPCTOF[AliPID::kPion];
    // probKaon = probTPCTOF[AliPID::kKaon];
  }
  else { // if you don't have both TOF and TPC, try only TPC
    fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTPC);
    AliDebug(2, "We did not find the detector mask for TOF + TPC, let's see only TPC");
    detUsed = fPIDCombined->ComputeProbabilities(bachelor, fPIDResponse, probTPCTOF);
    AliDebug(2, Form(" detUsed (TPC case) = %d", detUsed));
    if (detUsed == (UInt_t)fPIDCombined->GetDetectorMask()) {
      probProton = probTPCTOF[AliPID::kProton];
      // probPion = probTPCTOF[AliPID::kPion];
      // probKaon = probTPCTOF[AliPID::kKaon];
      AliDebug(2, Form("TPC only worked: probProton will be set to %f", probTPCTOF[AliPID::kProton]));
    }
    else {
      AliDebug(2, "Only TPC did not work...");
    }
    // resetting mask to ask for both TPC+TOF
    fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTPC+AliPIDResponse::kDetTOF);
  }
  AliDebug(2, Form("probProton = %f", probProton));

  // now we get the TPC and TOF single PID probabilities (only for Proton, or the tree will explode :) )
  Double_t probProtonTPC = -1.;
  Double_t probProtonTOF = -1.;
  Double_t pidTPC[AliPID::kSPECIES]={-1.};
  Double_t pidTOF[AliPID::kSPECIES]={-1.};
  Int_t respTPC = fPIDResponse->ComputePIDProbability(AliPIDResponse::kDetTPC, bachelor, AliPID::kSPECIES, pidTPC);
  Int_t respTOF = fPIDResponse->ComputePIDProbability(AliPIDResponse::kDetTOF, bachelor, AliPID::kSPECIES, pidTOF);
  if (respTPC == AliPIDResponse::kDetPidOk) probProtonTPC = pidTPC[AliPID::kProton];
  if (respTOF == AliPIDResponse::kDetPidOk) probProtonTOF = pidTOF[AliPID::kProton];

  // checking V0 status (on-the-fly vs offline)
  if ( !( !onFlyV0 || (onFlyV0 && fUseOnTheFlyV0) ) ) {
    AliDebug(2, "On-the-fly discarded");
    return;
  }

  if ( (((cutsAnal->IsSelected(part,AliRDHFCuts::kAll, aod)) & (AliRDHFCutsLctoV0::kLcToK0Spr)) == (AliRDHFCutsLctoV0::kLcToK0Spr)) ) {
    nSelectedAnal++;
  }

  if ( !(cutsAnal->IsInFiducialAcceptance(part->Pt(), part->Y(4122))) ) return;

  if ( !( ( (cutsAnal->IsSelected(part, AliRDHFCuts::kTracks, aod)) & (AliRDHFCutsLctoV0::kLcToK0Spr)) == (AliRDHFCutsLctoV0::kLcToK0Spr) ) ) { // esd track cuts
    if (isLc) AliDebug(2, "SIGNAL candidate rejected");
    AliDebug(2, "No: Analysis cuts kTracks level NOT passed");
    return;
  }
  else {
    AliDebug(2, "Yes: Analysis cuts kTracks level passed");
  }

  Int_t pdgCand = 4122;
  Int_t pdgDgLctoV0bachelor[2] = {211, 3122}; // case of wrong decay! Lc --> L + pi
  Int_t pdgDgV0toDaughters[2] = {2212, 211}; // case of wrong decay! Lc --> L + pi
  Int_t isLc2LBarpi = 0, isLc2Lpi = 0;
  Int_t currentLabel = part->GetLabel();
  Int_t mcLabel = 0;
  if (fUseMCInfo) {
    mcLabel = part->MatchToMC(pdgCand, pdgDgLctoV0bachelor[1], pdgDgLctoV0bachelor, pdgDgV0toDaughters, mcArray, kTRUE);
    if (mcLabel >= 0) {
      if (bachelor->Charge() == -1) isLc2LBarpi=1;
      if (bachelor->Charge() == +1) isLc2Lpi=1;
    }
  }

  Int_t pdgDg2prong[2] = {211, 211};
  Int_t labelK0S = 0;
  Int_t isK0S = 0;
  if (fUseMCInfo) {
    labelK0S = v0part->MatchToMC(310, mcArray, 2, pdgDg2prong);
    if (labelK0S>=0) isK0S = 1;
  }

  pdgDg2prong[0] = 211;
  pdgDg2prong[1] = 2212;
  Int_t isLambda = 0;
  Int_t isLambdaBar = 0;
  Int_t lambdaLabel = 0;
  if (fUseMCInfo) {
    lambdaLabel = v0part->MatchToMC(3122, mcArray, 2, pdgDg2prong);
    if (lambdaLabel >= 0) {
      AliAODMCParticle *lambdaTrack = (AliAODMCParticle*)mcArray->At(lambdaLabel);
      if (lambdaTrack->GetPdgCode() == 3122) isLambda = 1;
      else if (lambdaTrack->GetPdgCode() == -3122) isLambdaBar = 1;
    }
  }

  pdgDg2prong[0] = 11;
  pdgDg2prong[1] = 11;
  Int_t isGamma = 0;
  Int_t gammaLabel = 0;
  if (fUseMCInfo) {
    gammaLabel = v0part->MatchToMC(22, mcArray, 2, pdgDg2prong);
    if (gammaLabel>=0) {
      AliAODMCParticle *gammaTrack = (AliAODMCParticle*)mcArray->At(gammaLabel);
      if (gammaTrack->GetPdgCode() == 22) isGamma = 1;
    }
  }

  Int_t pdgTemp = -1;
  if (currentLabel != -1){
    AliAODMCParticle *tempPart = (AliAODMCParticle*)mcArray->At(currentLabel);
    pdgTemp = tempPart->GetPdgCode();
  }
  if (isLc) AliDebug(2, Form("Signal: Candidate is a Lc in K0s+p"));
  else if (isLc2LBarpi) AliDebug(2, Form("Background: Candidate is a Lc in Lbar + pi"));
  else if (isLc2Lpi) AliDebug(2, Form("Background: Candidate is a Lc in L + pi"));
  else AliDebug(2, Form("Pure bkg: Candidate has pdg = %d", pdgTemp));
  if (isK0S) AliDebug(2, Form("V0 is a K0S"));
  else if (isLambda)  AliDebug(2, Form("V0 is a Lambda"));
  else if (isLambdaBar)  AliDebug(2, Form("V0 is a LambdaBar"));
  else if (isGamma)  AliDebug(2, Form("V0 is a Gamma"));
  //else AliDebug(2, Form("V0 is something else!!"));

  Double_t invmassLc2Lpi = part->InvMassLctoLambdaPi();
  Double_t invmassLambda = v0part->MassLambda();
  Double_t invmassLambdaBar = v0part->MassAntiLambda();

  //Double_t nSigmaITSpr = -999.;
  Double_t nSigmaTPCpr = -999.;
  Double_t nSigmaTOFpr = -999.;

  //Double_t nSigmaITSpi = -999.;
  Double_t nSigmaTPCpi = -999.;
  Double_t nSigmaTOFpi = -999.;

  //Double_t nSigmaITSka = -999.;
  Double_t nSigmaTPCka = -999.;
  Double_t nSigmaTOFka = -999.;

  /*
    cutsAnal->GetPidHF()->GetnSigmaITS(bachelor, 4, nSigmaITSpr);
    cutsAnal->GetPidHF()->GetnSigmaTPC(bachelor, 4, nSigmaTPCpr);
    cutsAnal->GetPidHF()->GetnSigmaTOF(bachelor, 4, nSigmaTOFpr);
    cutsAnal->GetPidHF()->GetnSigmaITS(bachelor, 2, nSigmaITSpi);
    cutsAnal->GetPidHF()->GetnSigmaTPC(bachelor, 2, nSigmaTPCpi);
    cutsAnal->GetPidHF()->GetnSigmaTOF(bachelor, 2, nSigmaTOFpi);
    cutsAnal->GetPidHF()->GetnSigmaITS(bachelor, 3, nSigmaITSka);
    cutsAnal->GetPidHF()->GetnSigmaTPC(bachelor, 3, nSigmaTPCka);
    cutsAnal->GetPidHF()->GetnSigmaTOF(bachelor, 3, nSigmaTOFka);
  */

  if (fUsePIDresponseForNsigma) {
    nSigmaTPCpi = fPIDResponse->NumberOfSigmasTPC(bachelor, (AliPID::kPion));
    nSigmaTPCka = fPIDResponse->NumberOfSigmasTPC(bachelor, (AliPID::kKaon));
    nSigmaTPCpr = fPIDResponse->NumberOfSigmasTPC(bachelor, (AliPID::kProton));
    nSigmaTOFpi = fPIDResponse->NumberOfSigmasTOF(bachelor, (AliPID::kPion));
    nSigmaTOFka = fPIDResponse->NumberOfSigmasTOF(bachelor, (AliPID::kKaon));
    nSigmaTOFpr = fPIDResponse->NumberOfSigmasTOF(bachelor, (AliPID::kProton));
  }
  else {
    cutsAnal->GetPidHF()->GetnSigmaTPC(bachelor, (AliPID::kPion), nSigmaTPCpi);
    cutsAnal->GetPidHF()->GetnSigmaTPC(bachelor, (AliPID::kKaon), nSigmaTPCka);
    cutsAnal->GetPidHF()->GetnSigmaTPC(bachelor, (AliPID::kProton), nSigmaTPCpr);
    cutsAnal->GetPidHF()->GetnSigmaTOF(bachelor, (AliPID::kPion), nSigmaTOFpi);
    cutsAnal->GetPidHF()->GetnSigmaTOF(bachelor, (AliPID::kKaon), nSigmaTOFka);
    cutsAnal->GetPidHF()->GetnSigmaTOF(bachelor, (AliPID::kProton), nSigmaTOFpr);    
  }
  
  Double_t ptLcMC = -1;
  Double_t weightPythia = -1, weight5LHC13d3 = -1, weight5LHC13d3Lc = -1; 

  AliAODMCParticle *partLcMC = 0x0;
   
  if (fUseMCInfo) {
    if (iLctopK0s >= 0) {
      partLcMC = (AliAODMCParticle*)mcArray->At(iLctopK0s);
      ptLcMC = partLcMC->Pt();
      //Printf("--------------------- Reco pt = %f, MC particle pt = %f", part->Pt(), ptLcMC);
      weightPythia = fFuncWeightPythia->Eval(ptLcMC);
      //weight5LHC13d3 = fFuncWeightFONLL5overLHC13d3->Eval(ptLcMC);
      //weight5LHC13d3Lc = fFuncWeightFONLL5overLHC13d3Lc->Eval(ptLcMC);
    }
  }

  Double_t weightNch = 1;
  if (fUseMCInfo) {
    //Int_t nChargedMCPhysicalPrimary=AliVertexingHFUtils::GetGeneratedPhysicalPrimariesInEtaRange(mcArray,-1.0,1.0);
    //  if(nChargedMCPhysicalPrimary > 0)

    if(fNTracklets_1 > 0){
      if(!fHistoMCNch) AliDebug(2, "Input histos to evaluate Nch weights missing"); 
      if(fHistoMCNch) weightNch *= fHistoMCNch->GetBinContent(fHistoMCNch->FindBin(fNTracklets_1));
    }
  }

  // Fill candidate variable Tree (track selection, V0 invMass selection)
  // if (!onFlyV0 && isInV0window && isInCascadeWindow && part->CosV0PointingAngle()>0.99 && TMath::Abs(nSigmaTPCpr) <= 3 && 
  //    TMath::Abs(nSigmaTOFpr) <= 3 && v0part->Getd0Prong(0) < 20 && v0part->Getd0Prong(1) < 20) {
  if (!onFlyV0 && isInV0window && isInCascadeWindow){
    
    EBachelor bachCode = kBachInvalid;
    EK0S k0SCode = kK0SInvalid;
    if (fUseMCInfo) {
      bachCode = CheckBachelor(part, bachelor, mcArray);
      k0SCode = CheckK0S(part, v0part, mcArray);
    }
    
    AliAODTrack *v0pos = (AliAODTrack*)part->Getv0PositiveTrack();
    AliAODTrack *v0neg = (AliAODTrack*)part->Getv0NegativeTrack();
    //Int_t kfResult;
    TVector3 mom1(bachelor->Px(), bachelor->Py(), bachelor->Pz());
    TVector3 mom2(v0part->Px(), v0part->Py(), v0part->Pz());
    TVector3 momTot(part->Px(), part->Py(), part->Pz());
    
    Double_t Ql1 = mom1.Dot(momTot)/momTot.Mag();
    Double_t Ql2 = mom2.Dot(momTot)/momTot.Mag();
    
    Double_t alphaArmLc = (Ql1 - Ql2)/(Ql1 + Ql2);
    Double_t alphaArmLcCharge = ( bachelor->Charge() > 0 ? (Ql1 - Ql2)/(Ql1 + Ql2) : (Ql2 - Ql1)/(Ql1 + Ql2) );
    Double_t ptArmLc = mom1.Perp(momTot);
    
    Double_t massK0SPDG = TDatabasePDG::Instance()->GetParticle(310)->Mass();    // mass K0S PDG
    Double_t massPrPDG = TDatabasePDG::Instance()->GetParticle(2212)->Mass();    // mass Proton PDG
    Double_t massLcPDG = TDatabasePDG::Instance()->GetParticle(4122)->Mass();    // mass Lc PDG
    
    // Cosine of proton emission angle (theta*) in the rest frame of the mother particle
    // for prong ip (0 or 1) with mass hypotheses massLcPDG for mother particle (from AliAODRecoDecay)
    Double_t pStar = TMath::Sqrt((massLcPDG*massLcPDG - massPrPDG*massPrPDG - massK0SPDG*massK0SPDG)*(massLcPDG*massLcPDG - massPrPDG*massPrPDG - massK0SPDG*massK0SPDG) - 4.*massPrPDG*massPrPDG*massK0SPDG*massK0SPDG)/(2.*massLcPDG);
    Double_t e = part->E(4122);
    Double_t beta = part->P()/e;
    Double_t gamma = e/massLcPDG;
    //Double_t cts = (Ql1/gamma-beta*TMath::Sqrt(pStar*pStar+massPrPDG*massPrPDG))/pStar;
    
    // Cosine of proton emission angle (theta*) in the rest frame of the mother particle
    // (from AliRDHFCutsLctoV0)
    TLorentzVector vpr, vk0s,vlc;
    vpr.SetXYZM(part->PxProng(0), part->PyProng(0), part->PzProng(0), massPrPDG);
    vk0s.SetXYZM(part->PxProng(1), part->PyProng(1), part->PzProng(1), massK0SPDG);
    vlc = vpr + vk0s;
    TVector3 vboost = vlc.BoostVector();
    vpr.Boost(-vboost);
    Double_t cts = TMath::Cos(vpr.Angle(vlc.Vect()));   
    
    Double_t countTreta1corr = fNTracklets_1; 
    
    AliAODEvent* aodEvent = dynamic_cast<AliAODEvent*>(fInputEvent);
    
    if(fUseMultCorrection){
      TProfile* estimatorAvg = GetEstimatorHistogram(aodEvent);
      if(estimatorAvg) {
        countTreta1corr = static_cast<Int_t>(AliVertexingHFUtils::GetCorrectedNtracklets(estimatorAvg, fNTracklets_1, fVtx1->GetZ(), fRefMult));
      }
    }    
    
    //AliAODVertex *primvert = dynamic_cast<AliAODVertex*>(part->GetPrimaryVtx());
    
    Double_t d0z0bach[2], covd0z0bach[3];
    //bachelor->PropagateToDCA(primvert, fBField, kVeryBig, d0z0bach, covd0z0bach);
    bachelor->PropagateToDCA(fVtx1, fBField, kVeryBig, d0z0bach, covd0z0bach);
    Double_t tx[3];
    bachelor->GetXYZ(tx);
    //tx[0] -= primvert->GetX();
    //tx[1] -= primvert->GetY();
    //tx[2] -= primvert->GetZ();
    tx[0] -= fVtx1->GetX();
    tx[1] -= fVtx1->GetY();
    tx[2] -= fVtx1->GetZ();
    Double_t innerpro = tx[0]*part->Px() + tx[1]*part->Py();
    Double_t signd0 = 1.;
    if(innerpro < 0.) signd0 = -1.;
    
    signd0 = signd0*TMath::Abs(d0z0bach[0]);

    Double_t BDTResponse = -1;
    Double_t tmva = -1;
    
    if(!fFillTree){
      std::vector<Double_t> inputVars(fNVars);
      if (fNVars == 11) {
	inputVars[0] = invmassK0s;
	inputVars[1] = part->Getd0Prong(0);
	inputVars[2] = part->Getd0Prong(1);
	inputVars[3] = (part->DecayLengthV0())*0.497/(v0part->P());
	inputVars[4] = part->CosV0PointingAngle();
	inputVars[5] = cts;
	inputVars[6] = signd0;
	inputVars[7] = nSigmaTOFpr;
	inputVars[8] = nSigmaTPCpr;
	inputVars[9] = nSigmaTPCpi;
	inputVars[10] = nSigmaTPCka;
      }
      else if (fNVars == 8) {
	inputVars[0] = invmassK0s;
	inputVars[1] = part->Getd0Prong(0);
	inputVars[2] = part->Getd0Prong(1);
	inputVars[3] = (part->DecayLengthV0())*0.497/(v0part->P());
	inputVars[4] = part->CosV0PointingAngle();
	inputVars[5] = cts;
	inputVars[6] = signd0;
	inputVars[7] = probProton;
      }
      
      for (Int_t i = 0; i < fNVars; i++) {
	fVarsTMVA[i] = inputVars[i];
      }
      
      if (fUseXmlWeightsFile || fUseXmlFileFromCVMFS) tmva = fReader->EvaluateMVA("BDT method");
      if (fUseWeightsLibrary) BDTResponse = fBDTReader->GetMvaValue(inputVars);
      //Printf("BDTResponse = %f, invmassLc = %f", BDTResponse, invmassLc);
      //Printf("tmva = %f", tmva); 
      fBDTHisto->Fill(BDTResponse, invmassLc); 
      fBDTHistoTMVA->Fill(tmva, invmassLc); 
      if (fDebugHistograms) {
	if (fUseXmlWeightsFile || fUseXmlFileFromCVMFS) BDTResponse = tmva; // we fill the debug histogram with the output from the xml file
	fBDTHistoVsMassK0S->Fill(BDTResponse, invmassK0s);
	fBDTHistoVstImpParBach->Fill(BDTResponse, part->Getd0Prong(0));
	fBDTHistoVstImpParV0->Fill(BDTResponse, part->Getd0Prong(1));
	fBDTHistoVsBachelorPt->Fill(BDTResponse, bachelor->Pt());
	fBDTHistoVsCombinedProtonProb->Fill(BDTResponse, probProton);
	fBDTHistoVsCtau->Fill(BDTResponse, (part->DecayLengthV0())*0.497/(v0part->P()));
	fBDTHistoVsCosPAK0S->Fill(BDTResponse, part->CosV0PointingAngle());
	fBDTHistoVsSignd0->Fill(BDTResponse, signd0);
	fBDTHistoVsCosThetaStar->Fill(BDTResponse, cts);
	fBDTHistoVsnSigmaTPCpr->Fill(BDTResponse, nSigmaTPCpr);
	fBDTHistoVsnSigmaTOFpr->Fill(BDTResponse, nSigmaTOFpr);
	fBDTHistoVsnSigmaTPCpi->Fill(BDTResponse, nSigmaTPCpi);
	fBDTHistoVsnSigmaTPCka->Fill(BDTResponse, nSigmaTPCka);
	fBDTHistoVsBachelorP->Fill(BDTResponse, bachelor->P());
	fBDTHistoVsBachelorTPCP->Fill(BDTResponse, bachelor->GetTPCmomentum());
	fHistoNsigmaTPC->Fill(bachelor->P(), nSigmaTPCpr);
	fHistoNsigmaTOF->Fill(bachelor->P(), nSigmaTOFpr);
      }
    }
    
    
    Int_t checkOrigin = -1;
    Bool_t isFromSigmaC = kFALSE;
    AliAODMCParticle *mcpartMum = 0x0;
    
    Double_t sigmaCpt = -1;
    Double_t cosThetaStarSoftPi = -1.1;
    Int_t labelSoftPi = -1;
    Double_t ptsigmacMC = -1;
    Double_t ptlambdacMC = -1;
    Double_t ysigmacMC = -9;
    Double_t ylambdacMC = -9;
    Double_t pointlcsc[6];
    
    if(partLcMC){
      
      checkOrigin = AliVertexingHFUtils::CheckOrigin(mcArray, partLcMC, kTRUE);
      
      fhistMCSpectrumAccLc->Fill(partLcMC->Pt(), kRecoPID, checkOrigin);     

      Int_t indSc = partLcMC->GetMother();
      if(indSc >= 0){
	mcpartMum = (AliAODMCParticle*)mcArray->At(indSc); // Sigmac candidate 
	Int_t pdgLcMum = TMath::Abs(mcpartMum->GetPdgCode());
	if(pdgLcMum == 4112 || pdgLcMum == 4222) isFromSigmaC = kTRUE;
      }
 
      if(isFromSigmaC){
	pointlcsc[0] = partLcMC->Pt();
	pointlcsc[1] = kRecoLcPID;
	pointlcsc[2] = checkOrigin;
	pointlcsc[3] = partLcMC->Y();
	pointlcsc[4] = mcpartMum->Pt();
	pointlcsc[5] = mcpartMum->Y();
	  
	fhistMCSpectrumAccLcFromSc->Fill(pointlcsc); // !!!! really needed?
	
      }// end is from sigmaC
    }// end if partLcMC

    if(mcpartMum && isFromSigmaC){
    //if(mcpartMum){
      ptsigmacMC = mcpartMum->Pt();
      //if(mcpartMum->GetNDaughters() != 2) return;
      if(mcpartMum->GetNDaughters() == 2) {
	for(Int_t k = mcpartMum->GetDaughterLabel(0); k <= mcpartMum->GetDaughterLabel(1); k++){
	  if(k >= 0){
	    AliAODMCParticle *mcpartScdau = (AliAODMCParticle*)mcArray->At(k);
	    if(TMath::Abs(mcpartScdau->GetPdgCode()) == 211){ // daughter is soft pion
	      labelSoftPi = k;
	    }
	    else if(TMath::Abs(mcpartScdau->GetPdgCode()) == 4122){ // daughter is LambdaC
	      //mcpartLc = mcpartScdau;
	      //partLcMC = mcpartScdau;
	      //ptlambdacMC = partLcMC->Pt();
	      //ylambdacMC = partLcMC->Y();
	      ptlambdacMC = mcpartScdau->Pt();
	      ylambdacMC = mcpartScdau->Y();
	    }
	  }
	}
      }
    }


    if(isLcAnalysis){ //fill the tree with only Lc candidates
      if(fFillTree && isLcAnalysis){
	if (fUseMCInfo) {   //  save full tree if on MC
	  fCandidateVariables[0] = invmassLc;
	  fCandidateVariables[1] = invmassLc2Lpi;
	  fCandidateVariables[2] = invmassK0s;
	  fCandidateVariables[3] = invmassLambda;
	  fCandidateVariables[4] = invmassLambdaBar;
	  fCandidateVariables[5] = part->CosV0PointingAngle();
	  fCandidateVariables[6] = dcaV0;
	  fCandidateVariables[7] = part->Getd0Prong(0);
	  fCandidateVariables[8] = part->Getd0Prong(1);
	  fCandidateVariables[9] = nSigmaTPCpr;
	  fCandidateVariables[10] = nSigmaTOFpr;
	  fCandidateVariables[11] = bachelor->Pt();
	  fCandidateVariables[12] = v0pos->Pt();
	  fCandidateVariables[13] = v0neg->Pt();
	  fCandidateVariables[14] = v0part->Getd0Prong(0);
	  fCandidateVariables[15] = v0part->Getd0Prong(1);
	  fCandidateVariables[16] = v0part->Pt();
	  fCandidateVariables[17] = v0part->InvMass2Prongs(0,1,11,11);
	  fCandidateVariables[18] = part->Pt();
	  fCandidateVariables[19] = probProton;
	  fCandidateVariables[20] = part->Eta();
	  fCandidateVariables[21] = v0pos->Eta();
	  fCandidateVariables[22] = v0neg->Eta();
	  fCandidateVariables[23] = probProtonTPC;
	  fCandidateVariables[24] = probProtonTOF;
	  fCandidateVariables[25] = bachelor->Eta();      
	  fCandidateVariables[26] = part->P();
	  fCandidateVariables[27] = bachelor->P();
	  fCandidateVariables[28] = v0part->P();
	  fCandidateVariables[29] = v0pos->P();
	  fCandidateVariables[30] = v0neg->P();
	  fCandidateVariables[31] = v0part->Eta();
	  fCandidateVariables[32] = ptLcMC;
	  fCandidateVariables[33] = part->DecayLengthV0();
	  fCandidateVariables[34] = bachCode;
	  fCandidateVariables[35] = k0SCode;
	  fCandidateVariables[36] = v0part->AlphaV0();
	  fCandidateVariables[37] = v0part->PtArmV0();	
	  fCandidateVariables[38] = cts;
	  fCandidateVariables[39] = weightPythia;
	  fCandidateVariables[40] = sigmaCpt;
	  fCandidateVariables[41] = cosThetaStarSoftPi;
	  fCandidateVariables[42] = weightNch;
	  fCandidateVariables[43] = fNTracklets_1;      
	  fCandidateVariables[44] = countTreta1corr;
	  fCandidateVariables[45] = signd0;
	  fCandidateVariables[46] = fCentrality;
	  fCandidateVariables[47] = fNTracklets_All;
	  fCandidateVariables[48] = checkOrigin;
	  fCandidateVariables[49] = nSigmaTPCpi;
	  fCandidateVariables[50] = nSigmaTPCka;
	  fCandidateVariables[51] = 0;
	  fCandidateVariables[52] = ptArmLc;
	  fCandidateVariables[53] = alphaArmLc;
	}      
	else { //remove MC-only variables from tree if data
	  fCandidateVariables[0] = invmassLc;
	  fCandidateVariables[1] = v0part->AlphaV0();
	  fCandidateVariables[2] = invmassK0s;
	  fCandidateVariables[3] = invmassLambda;
	  fCandidateVariables[4] = invmassLambdaBar;
	  fCandidateVariables[5] = part->CosV0PointingAngle();
	  fCandidateVariables[6] = dcaV0;
	  fCandidateVariables[7] = part->Getd0Prong(0);
	  fCandidateVariables[8] = part->Getd0Prong(1);
	  fCandidateVariables[9] = nSigmaTPCpr;
	  fCandidateVariables[10] = nSigmaTOFpr;
	  fCandidateVariables[11] = bachelor->Pt();
	  fCandidateVariables[12] = v0pos->Pt();
	  fCandidateVariables[13] = v0neg->Pt();
	  fCandidateVariables[14] = v0part->Getd0Prong(0);
	  fCandidateVariables[15] = v0part->Getd0Prong(1);
	  fCandidateVariables[16] = v0part->Pt();
	  fCandidateVariables[17] = sigmaCpt;	      
	  fCandidateVariables[18] = part->Pt();
	  fCandidateVariables[19] = probProton;
	  fCandidateVariables[20] = v0pos->Eta();
	  fCandidateVariables[21] = ptArmLc;
	  fCandidateVariables[22] = bachelor->Eta();
	  fCandidateVariables[23] = v0part->P();
	  fCandidateVariables[24] = part->DecayLengthV0();
	  fCandidateVariables[25] = nSigmaTPCpi;
	  fCandidateVariables[26] = nSigmaTPCka;
	  fCandidateVariables[27] = fNTracklets_1;
	  fCandidateVariables[28] = countTreta1corr;
	  fCandidateVariables[29] = cts;
	  fCandidateVariables[30] = signd0;       
	  fCandidateVariables[31] = fCentrality;
	  fCandidateVariables[32] = fNTracklets_All;
	  fCandidateVariables[33] = -1;
	  fCandidateVariables[34] = v0part->PtArmV0();
	  fCandidateVariables[35] = 0;
	  fCandidateVariables[36] = cosThetaStarSoftPi;
	  fCandidateVariables[37] = alphaArmLc;
	}
      }
      
      if (fUseMCInfo) {
	if (isLc){
	  AliDebug(2, Form("Reco particle %d --> Filling Sgn", iLctopK0s));
	  if(fFillTree) fVariablesTreeSgn->Fill();
	  fHistoCodesSgn->Fill(bachCode, k0SCode);
	}
	else {
	  if (fFillOnlySgn == kFALSE){
	    AliDebug(2, "Filling Bkg");
	    if(fFillTree) fVariablesTreeBkg->Fill();
	    fHistoCodesBkg->Fill(bachCode, k0SCode);
	  }
	}
      }
      else {
	if(fFillTree) fVariablesTreeSgn->Fill();   
      }
    }
    
    else { //isLcAnalysis  

      //if(TMath::Abs(mass - 2.28646) > fLcMassWindowForSigmaC) return; //Lc mass window selection  
      //AliDebug(2, Form("Good Lc candidate , will loop over %d pions",fnSelSoftPi));
      
      Double_t pointSigma[9]; //Lcpt, deltam, Lcmass, cosThetaStarSoftPi, Sigmapt, origin, isRotated, bdtresp, softPiITSrefit
      
      
      pointSigma[0] = part->Pt();
      //pointSigma[2] = part->CosPointingAngle();
      pointSigma[2] = invmassLc;
      pointSigma[5] = checkOrigin;
      pointSigma[6] = 1;
      pointSigma[7] = -1;
      pointSigma[8] = 0;
      // Bool_t arrayVariableIsFilled = kFALSE;
      // if(pointS){    
      //   for(Int_t k = 0; k < 8; k++){
      // 	pointSigma[k] = pointS[k];
      //   }
      //   arrayVariableIsFilled = kTRUE;
      // }
      
      // Loop over soft pions
      Double_t p2 = part->P2();    
      for(Int_t isoft = 0; isoft < fnSelSoftPi; isoft++){
	Int_t indsof = ftrackArraySelSoftPi->At(isoft);
	AliAODTrack *tracksoft = (AliAODTrack*)aodEvent->GetTrack(indsof);    		
	if(mcpartMum){
	  if(TMath::Abs(tracksoft->GetLabel()) != labelSoftPi) continue;
	}
	
	Bool_t skip = kFALSE;
	for(Int_t k = 0; k < 2; k++){
	  if((Int_t)(part->GetProngID(k)) == tracksoft->GetID()){
	    //Printf("Skipping Lc candidate with itself");
	    skip = kTRUE;
	    break;
	  }
	}
	if(skip) continue;
	
	Double_t psoft[3], psoftOrig[3];
	tracksoft->PxPyPz(psoftOrig);
	psoft[0] = psoftOrig[0];
	psoft[1] = psoftOrig[1];
	psoft[2] = psoftOrig[2];
	Double_t pcand[3];
	part->PxPyPz(pcand);
	Double_t rotStep = 0.;
	
	if(tracksoft->TestFilterBit(AliAODTrack::kITSrefit)) pointSigma[8] = 1;
	
	//pointSigma[7] = 1;
	//if(fNRotations>1) rotStep=(fMaxAngleForRot-fMinAngleForRot)/(fNRotations-1); // -1 is to ensure that the last rotation is done with angle=fMaxAngleForRot     
	//for(Int_t irot = -1; irot < fNRotations; irot++){
	// tracks are rotated to provide further background, if required
	// ASSUMPTIONS: there is no need to repeat single track selection after rotation, because we just rotate in the transverse plane (-> pt and eta does not change; the aspects related to the detector are, like potential intersection of dead modules in the roated direction are not considered)
	
	//if(irot >= 0){
	//Double_t phirot = fMinAngleForRot + rotStep*irot;	
	//psoft[0] = psoftOrig[0]*TMath::Cos(phirot) - psoftOrig[1]*TMath::Sin(phirot);
	//psoft[1] = psoftOrig[0]*TMath::Sin(phirot) + psoftOrig[1]*TMath::Cos(phirot);
	//pointSigma[7] = 0;
	//}
	
	Double_t psigma[3] = {pcand[0]+psoft[0], pcand[1]+psoft[1], pcand[2]+psoft[2]};
	Double_t e1, e2;	      
	//Double_t cosThetaStarSoftPi = -1.1;
	
	if(TMath::Abs(invmassLc - 2.28646) < fLcMassWindowForSigmaC){// here we may be more restrictive and check also resp_only_pid, given that later is  done before filling the histogram
	  
	  e1 = TMath::Sqrt(invmassLc*invmassLc + p2);
	  e2 = TMath::Sqrt(0.019479785 + psoft[0]*psoft[0] + psoft[1]*psoft[1] + psoft[2]*psoft[2]);// 0.019479785 =  0.13957*0.13957
	  TLorentzVector lsum(psoft[0] + part->Px(), psoft[1] + part->Py(), psoft[2] + part->Pz(), e1 + e2);
	  sigmaCpt = lsum.Pt();
	  pointSigma[4] = sigmaCpt;
	  
	  if(sigmaCpt < fMinPtSigmacCand || sigmaCpt > fMaxPtSigmacCand) continue;
	  
	  pointlcsc[0] = ptlambdacMC;
	  pointlcsc[1] = kRecoPID;
	  pointlcsc[2] = checkOrigin;
	  pointlcsc[3] = ylambdacMC;
	  pointlcsc[4] = ptsigmacMC;
	  pointlcsc[5] = ysigmacMC;
	  if(mcpartMum && isFromSigmaC){
	    fhistMCSpectrumAccSc->Fill(ptsigmacMC, kRecoPID, checkOrigin);	      
	    fhistMCSpectrumAccLcFromSc->Fill(pointlcsc);
	  }
	  
	  cosThetaStarSoftPi = CosThetaStar(psigma, psoft, TDatabasePDG::Instance()->GetParticle(4222)->Mass(), TDatabasePDG::Instance()->GetParticle(211)->Mass());
	  pointSigma[3] = cosThetaStarSoftPi;
	  Double_t deltaM = lsum.M() - invmassLc;
	  pointSigma[1] = deltaM;
	  
	  if(deltaM < fSigmaCDeltaMassWindow){ // good candidate
	    
	    if(fFillTree){
	      if (fUseMCInfo) {   //  save full tree if on MC
		fCandidateVariables[0] = invmassLc;
		fCandidateVariables[1] = invmassLc2Lpi;
		fCandidateVariables[2] = invmassK0s;
		fCandidateVariables[3] = invmassLambda;
		fCandidateVariables[4] = invmassLambdaBar;
		fCandidateVariables[5] = part->CosV0PointingAngle();
		fCandidateVariables[6] = dcaV0;
		fCandidateVariables[7] = part->Getd0Prong(0);
		fCandidateVariables[8] = part->Getd0Prong(1);
		fCandidateVariables[9] = nSigmaTPCpr;
		fCandidateVariables[10] = nSigmaTOFpr;
		fCandidateVariables[11] = bachelor->Pt();
		fCandidateVariables[12] = v0pos->Pt();
		fCandidateVariables[13] = v0neg->Pt();
		fCandidateVariables[14] = v0part->Getd0Prong(0);
		fCandidateVariables[15] = v0part->Getd0Prong(1);
		fCandidateVariables[16] = v0part->Pt();
		fCandidateVariables[17] = v0part->InvMass2Prongs(0,1,11,11);
		fCandidateVariables[18] = part->Pt();
		fCandidateVariables[19] = probProton;
		fCandidateVariables[20] = part->Eta();
		fCandidateVariables[21] = v0pos->Eta();
		fCandidateVariables[22] = v0neg->Eta();
		fCandidateVariables[23] = probProtonTPC;
		fCandidateVariables[24] = probProtonTOF;
		fCandidateVariables[25] = bachelor->Eta();      
		fCandidateVariables[26] = part->P();
		fCandidateVariables[27] = bachelor->P();
		fCandidateVariables[28] = v0part->P();
		fCandidateVariables[29] = v0pos->P();
		fCandidateVariables[30] = v0neg->P();
		fCandidateVariables[31] = v0part->Eta();
		fCandidateVariables[32] = ptLcMC;
		fCandidateVariables[33] = part->DecayLengthV0();
		fCandidateVariables[34] = bachCode;
		fCandidateVariables[35] = k0SCode;
		fCandidateVariables[36] = v0part->AlphaV0();
		fCandidateVariables[37] = v0part->PtArmV0();
		
		AliDebug(2, Form("v0pos->GetStatus() & AliESDtrack::kITSrefit= %d, v0neg->GetStatus() & AliESDtrack::kITSrefit = %d, v0pos->GetTPCClusterInfo(2, 1)= %f, v0neg->GetTPCClusterInfo(2, 1) = %f", (Int_t)(v0pos->GetStatus() & AliESDtrack::kITSrefit), (Int_t)(v0pos->GetStatus() & AliESDtrack::kITSrefit), v0pos->GetTPCClusterInfo(2, 1), v0neg->GetTPCClusterInfo(2, 1)));
		fCandidateVariables[38] = cts;
		fCandidateVariables[39] = weightPythia;
		fCandidateVariables[40] = sigmaCpt;
		fCandidateVariables[41] = cosThetaStarSoftPi;
		fCandidateVariables[42] = weightNch;
		fCandidateVariables[43] = fNTracklets_1;      
		fCandidateVariables[44] = countTreta1corr;
		fCandidateVariables[45] = signd0;
		fCandidateVariables[46] = fCentrality;
		fCandidateVariables[47] = fNTracklets_All;
		fCandidateVariables[48] = checkOrigin;
		fCandidateVariables[49] = nSigmaTPCpi;
		fCandidateVariables[50] = nSigmaTPCka;
		//fCandidateVariables[51] = bachelor->GetTPCmomentum();
		fCandidateVariables[51] = deltaM;
		fCandidateVariables[52] = ptArmLc;
		fCandidateVariables[53] = alphaArmLc;
	      }      
	      else { //remove MC-only variables from tree if data
		fCandidateVariables[0] = invmassLc;
		fCandidateVariables[1] = v0part->AlphaV0();
		fCandidateVariables[2] = invmassK0s;
		fCandidateVariables[3] = invmassLambda;
		fCandidateVariables[4] = invmassLambdaBar;
		fCandidateVariables[5] = part->CosV0PointingAngle();
		fCandidateVariables[6] = dcaV0;
		fCandidateVariables[7] = part->Getd0Prong(0);
		fCandidateVariables[8] = part->Getd0Prong(1);
		fCandidateVariables[9] = nSigmaTPCpr;
		fCandidateVariables[10] = nSigmaTOFpr;
		fCandidateVariables[11] = bachelor->Pt();
		fCandidateVariables[12] = v0pos->Pt();
		fCandidateVariables[13] = v0neg->Pt();
		fCandidateVariables[14] = v0part->Getd0Prong(0);
		fCandidateVariables[15] = v0part->Getd0Prong(1);
		fCandidateVariables[16] = v0part->Pt();
		fCandidateVariables[17] = sigmaCpt;	      
		fCandidateVariables[18] = part->Pt();
		fCandidateVariables[19] = probProton;
		fCandidateVariables[20] = v0pos->Eta();
		//fCandidateVariables[21] = bachelor->P();
		fCandidateVariables[21] = ptArmLc;
		fCandidateVariables[22] = bachelor->Eta();
		fCandidateVariables[23] = v0part->P();
		fCandidateVariables[24] = part->DecayLengthV0();
		fCandidateVariables[25] = nSigmaTPCpi;
		fCandidateVariables[26] = nSigmaTPCka;
		fCandidateVariables[27] = fNTracklets_1;
		fCandidateVariables[28] = countTreta1corr;
		fCandidateVariables[29] = cts;
		fCandidateVariables[30] = signd0;       
		fCandidateVariables[31] = fCentrality;
		fCandidateVariables[32] = fNTracklets_All;
		fCandidateVariables[33] = -1;
		fCandidateVariables[34] = v0part->PtArmV0();
		fCandidateVariables[35] = deltaM;
		fCandidateVariables[36] = cosThetaStarSoftPi;
		fCandidateVariables[37] = alphaArmLc;
	      }
	    }
	    
	    // fill multiplicity histograms for events with a candidate   
	    //fHistNtrUnCorrEvWithCand->Fill(fNTracklets_1, weightNch);
	    //fHistNtrCorrEvWithCand->Fill(countTreta1corr, weightNch);
	    if (fUseMCInfo) {
	      if (isLc){
		AliDebug(2, Form("Reco particle %d --> Filling Sgn", iLctopK0s));
		if(fFillTree) fVariablesTreeSgn->Fill();
		fHistoCodesSgn->Fill(bachCode, k0SCode);
	      }
	      else {
		if (fFillOnlySgn == kFALSE){
		  AliDebug(2, "Filling Bkg");
		  if(fFillTree) fVariablesTreeBkg->Fill();
		  fHistoCodesBkg->Fill(bachCode, k0SCode);
		}
	      }
	    }
	    else {
	      if(fFillTree) fVariablesTreeSgn->Fill();   
	    }
	    
	    
	    if(!fFillTree){
	      if (fUseXmlWeightsFile || fUseXmlFileFromCVMFS) pointSigma[7] = tmva;
	      if (fUseWeightsLibrary) pointSigma[7] = BDTResponse;
	      if(fhSparseAnalysisSigma)  {
		if(!mcpartMum) fhSparseAnalysisSigma->Fill(pointSigma);
		else {
		  AliAODTrack *trkd = (AliAODTrack*)part->GetDaughter(0); // Daughter(0) of Cascade is always a proton
		  AliAODMCParticle* pProt = (AliAODMCParticle*)mcArray->At(TMath::Abs(trkd->GetLabel()));
		  if(TMath::Abs(pProt->GetPdgCode()) == 2212){
		    pointSigma[5] = ptsigmacMC;
		    pointSigma[0] = ptlambdacMC;
		    fhSparseAnalysisSigma->Fill(pointSigma);
		    fhistMCSpectrumAccSc->Fill(ptsigmacMC, kRecoPID, checkOrigin);	      
		    pointlcsc[0] = ptlambdacMC;
		    pointlcsc[1] = kRecoPID;
		    pointlcsc[2] = checkOrigin;
		    pointlcsc[3] = ylambdacMC;
		    pointlcsc[4] = ptsigmacMC;
		    pointlcsc[5] = ysigmacMC;
		    fhistMCSpectrumAccLcFromSc->Fill(pointlcsc);
		  }
		}
	      }
	    }
	  } 
	}
      }
    }
    
  }
  
  return;
  
  
}
//________________________________________________________________________
AliAnalysisTaskSESigmacTopK0Spi::EBachelor AliAnalysisTaskSESigmacTopK0Spi::CheckBachelor( AliAODRecoCascadeHF *part,
											   AliAODTrack* bachelor,
											   TClonesArray *mcArray ){
  
  //Printf("In CheckBachelor");

  /// function to check if the bachelor is a p, if it is a p but not from Lc
  /// to be filled for background candidates

  Int_t label = bachelor->GetLabel();
  if (label == -1) {
    return kBachFake;
  }

  AliAODMCParticle *mcpart = dynamic_cast<AliAODMCParticle*>(mcArray->At(TMath::Abs(label)));
  if (!mcpart) return kBachInvalid;
  Int_t pdg = mcpart->PdgCode();
  if (TMath::Abs(pdg) != 2212) {
    AliDebug(2, Form("Bachelor is not a p, but a particle with pdg code =  %d", pdg));
    return kBachNoProton;
  }
  else { // it is a proton
    //Int_t labelLc = part->GetLabel();
    Int_t labelLc = FindLcLabel(part, mcArray);
    //Printf(">>>>>>>>>>>>> label for Lc = %d", labelLc);
    Int_t bachelorMotherLabel = mcpart->GetMother();
    //Printf(">>>>>>>>>>>>> label for bachelorMother = %d", bachelorMotherLabel);
    if (bachelorMotherLabel == -1) {
      return kBachPrimary;
    }
    AliAODMCParticle *bachelorMother = dynamic_cast<AliAODMCParticle*>(mcArray->At(bachelorMotherLabel));
    if (!bachelorMother) return kBachInvalid;
    Int_t pdgMother = bachelorMother->PdgCode();
    if (TMath::Abs(pdgMother) != 4122) {
      AliDebug(2, Form("The proton does not come from a Lc, but from a particle with pdgcode = %d", pdgMother));
      return kBachNoLambdaMother;
    }
    else { // it comes from Lc
      if (labelLc != bachelorMotherLabel){
	//AliInfo(Form("The proton comes from a Lc, but it is not the candidate we are analyzing (label Lc = %d, label p mother = %d", labelLc, bachelorMotherLabel));
	AliDebug(2, Form("The proton comes from a Lc, but it is not the candidate we are analyzing"));
	return kBachDifferentLambdaMother;
      }
      else { // it comes from the correct Lc
	AliDebug(2, Form("The proton comes from a Lc, and it is EXACTLY the candidate we are analyzing"));
	return kBachCorrectLambdaMother;
      }
    }
  }

  return kBachInvalid;

}

//________________________________________________________________________
AliAnalysisTaskSESigmacTopK0Spi::EK0S AliAnalysisTaskSESigmacTopK0Spi::CheckK0S( AliAODRecoCascadeHF *part,
										 AliAODv0* v0part,
										 //AliAODTrack* v0part,
										 TClonesArray *mcArray ){
  
  /// function to check if the K0Spart is a p, if it is a p but not from Lc
  /// to be filled for background candidates

  //Printf(" CheckK0S");

  //Int_t labelMatchToMC = v0part->MatchToMC(310, mcArray);
  //Int_t label = v0part->GetLabel();
  Int_t labelFind = FindV0Label(v0part, mcArray);
  //Printf("\n >>>>>>>>>>>>> label for V0 = %d, from MatchToMC = %d, from Find = %d", label, labelMatchToMC, labelFind);
  if (labelFind == -1) {
    return kK0SFake;
  }

  AliAODMCParticle *mcpart = dynamic_cast<AliAODMCParticle*>(mcArray->At(labelFind));
  if (!mcpart) return kK0SInvalid;
  Int_t pdg = mcpart->PdgCode();
  if (TMath::Abs(pdg) != 310) {
    AliDebug(2, Form("K0Spart is not a K0S, but a particle with pdg code =  %d", pdg));
    //AliInfo(Form("K0Spart is not a K0S, but a particle with pdg code =  %d", pdg));
    return kK0SNoK0S;
  }
  else { // it is a K0S
    //Int_t labelLc = part->GetLabel();
    Int_t labelLc = FindLcLabel(part, mcArray);
    Int_t K0SpartMotherLabel = mcpart->GetMother();
    if (K0SpartMotherLabel == -1) {
      return kK0SWithoutMother;
    }
    AliAODMCParticle *K0SpartMother = dynamic_cast<AliAODMCParticle*>(mcArray->At(K0SpartMotherLabel));  // should be a K0 in case of signal
    if (!K0SpartMother) return kK0SInvalid;
    Int_t pdgMotherK0S = K0SpartMother->PdgCode();
    if (TMath::Abs(pdgMotherK0S) != 311) {
      AliDebug(2, Form("The K0S does not come from a K0, but from a particle with pdgcode = %d", pdgMotherK0S));
      //AliInfo(Form("The K0S does not come from a K0, but from a particle with pdgcode = %d", pdgMotherK0S));
      return kK0SNotFromK0;
    }
    else { // the K0S comes from a K0
      Int_t K0MotherLabel = K0SpartMother->GetMother(); // mother of K0 --> Lc in case of signal
      if (K0MotherLabel == -1) {
	return kK0Primary;
      }
      AliAODMCParticle *K0Mother = dynamic_cast<AliAODMCParticle*>(mcArray->At(K0MotherLabel));
      if (!K0Mother) return kK0SInvalid;
      Int_t pdgK0Mother = K0Mother->PdgCode();
      if (TMath::Abs(pdgK0Mother) != 4122) { // the K0 does not come from a Lc
	AliDebug(2, Form("The K0 does not come from a Lc, but from a particle with pdgcode = %d", pdgK0Mother));
	//AliInfo(Form("The K0 does not come from a Lc, but from a particle with pdgcode = %d", pdgK0Mother));
	return kK0NoLambdaMother;
      }
      else { // the K0 comes from Lc
	if (labelLc != K0MotherLabel){ // The K0 comes from a different Lc
	  AliDebug(2, Form("The K0S comes from a Lc, but it is not the candidate we are analyzing"));
	  //AliInfo(Form("The K0S comes from a Lc, but it is not the candidate we are analyzing"));
	  return kK0DifferentLambdaMother;
	}
	else { // it comes from the correct Lc
	  AliDebug(2, Form("The K0S comes from a Lc, and it is EXACTLY the candidate we are analyzing"));
	  //AliInfo(Form("The K0S comes from a Lc, and it is EXACTLY the candidate we are analyzing"));
	  return kK0CorrectLambdaMother;
	}
      }
    }
  }

  return kK0SInvalid;

}

//----------------------------------------------------------------------------
Int_t AliAnalysisTaskSESigmacTopK0Spi::FindV0Label(AliAODRecoDecay* v0part, TClonesArray *mcArray) const
{

  //Printf(" FindV0Label");

  /// finding the label of teh V0; inspired from AliAODRecoDecay::MatchToMC

  Int_t labMother[2] = {-1, -1};
  AliAODMCParticle *part = 0;
  AliAODMCParticle *mother = 0;
  Int_t dgLabels = -1;

  Int_t ndg = v0part->GetNDaughters();
  if (ndg != 2) {
    AliFatal(Form("IMPOSSIBLE!! There are %d daughters, but they should be 2!", ndg));
  }

  for(Int_t i = 0; i < 2; i++) {
    AliAODTrack *trk = (AliAODTrack*)v0part->GetDaughter(i);
    dgLabels = trk->GetLabel();
    if (dgLabels == -1) {
      //printf("daughter with negative label %d\n", dgLabels);
      return -1;
    }
    part = (AliAODMCParticle*)mcArray->At(TMath::Abs(dgLabels));
    if (!part) {
      //printf("no MC particle\n");
      return -1;
    }
    labMother[i] = part->GetMother();
    if (labMother[i] != -1){
      mother = (AliAODMCParticle*)mcArray->At(labMother[i]);
      if(!mother) {
	//printf("no MC mother particle\n");
	return -1;
      }
    }
    else {
      return -1;
    }
  }

  if (labMother[0] == labMother[1]) return labMother[0];

  else return -1;

}

//_____________________________________________________________________________
Int_t AliAnalysisTaskSESigmacTopK0Spi::FindLcLabel(AliAODRecoCascadeHF* cascade, TClonesArray *mcArray) const
{

  /// finding the label of teh V0; inspired from AliAODRecoDecay::MatchToMC

  //Printf(" FindLcLabel");

  AliAODMCParticle *part = 0;
  AliAODMCParticle *mother = 0;
  AliAODMCParticle *grandmother = 0;
  Int_t dgLabels = -1;

  Int_t ndg = cascade->GetNDaughters();
  if (ndg != 2) {
    AliFatal(Form("IMPOSSIBLE!! There are %d daughters, but they should be 2!", ndg));
  }

  // daughter 0 --> bachelor, daughter 1 --> V0
  AliAODTrack *trk = (AliAODTrack*)cascade->GetDaughter(0); // bachelor
  dgLabels = trk->GetLabel();
  if (dgLabels == -1 ) {
    //printf("daughter with negative label %d\n", dgLabels);
    return -1;
  }
  part = (AliAODMCParticle*)mcArray->At(TMath::Abs(dgLabels));
  if (!part) {
    //printf("no MC particle\n");
    return -1;
  }
  Int_t labMotherBach = part->GetMother();
  if (labMotherBach == -1){
    return -1;
  }
  mother = (AliAODMCParticle*)mcArray->At(labMotherBach);
  if(!mother) {
    //printf("no MC mother particle\n");
    return -1;
  }

  AliAODv0 *v0 = (AliAODv0*)cascade->GetDaughter(1); // V0
  dgLabels = FindV0Label(v0, mcArray);
  if (dgLabels == -1 ) {
    //printf("daughter with negative label (v0 was a fake) %d\n", dgLabels);
    return -1;
  }
  part = (AliAODMCParticle*)mcArray->At(TMath::Abs(dgLabels));
  if (!part) {
    //printf("no MC particle\n");
    return -1;
  }
  Int_t labMotherv0 = part->GetMother();
  if (labMotherv0 == -1){
    return -1;
  }
  mother = (AliAODMCParticle*)mcArray->At(labMotherv0);
  if(!mother) {
    //printf("no MC mother particle\n");
    return -1;
  }
  Int_t labGrandMotherv0 = mother->GetMother();
  if (labGrandMotherv0 == -1){
    return -1;
  }
  grandmother = (AliAODMCParticle*)mcArray->At(labGrandMotherv0);
  if(!grandmother) {
    //printf("no MC mother particle\n");
    return -1;
  }

  //Printf("labMotherBach = %d, labMotherv0 = %d, labGrandMotherv0 = %d", labMotherBach, labMotherv0, labGrandMotherv0);
  if (labGrandMotherv0 == labMotherBach) return labMotherBach;

  else return -1;
  
}

//____________________________________________________________________________
TProfile* AliAnalysisTaskSESigmacTopK0Spi::GetEstimatorHistogram(const AliVEvent* event){
  /// Get Estimator Histogram from period event->GetRunNumber();
  ///
  /// If you select SPD tracklets in |eta|<1 you should use type == 1
    
  Int_t runNo  = event->GetRunNumber();
  Int_t period = -1;   // pp: 0-LHC10b, 1-LHC10c, 2-LHC10d, 3-LHC10e
  
  if(fYearNumber==10){
    if(runNo>114930 && runNo<117223) period = 0;//10b
    if(runNo>119158 && runNo<120830) period = 1;//10c
    if(runNo>122373 && runNo<126438) period = 2;//10d
    if(runNo>127711 && runNo<130851) period = 3;//10e
    if(period<0 || period>3) return 0;
  }else if(fYearNumber==16){
    if(runNo>=252235 && runNo<=252375)period = 0;//16d
    if(runNo>=252603 && runNo<=253591)period = 1;//16e
    if(runNo>=254124 && runNo<=254332)period = 2;//16g
    if(runNo>=254378  && runNo<=255469 )period = 3;//16h_1
    if(runNo>=254418  && runNo<=254422 )period = 4;//16h_2 negative mag
    if(runNo>=256146  && runNo<=256420 )period = 5;//16j
    if(runNo>=256504  && runNo<=258537 )period = 6;//16k
    if(runNo>=258883  && runNo<=260187)period = 7;//16l
    if(runNo>=262395  && runNo<=264035 )period = 8;//16o
    if(runNo>=264076  && runNo<=264347 )period = 9;//16p
  }else if(fYearNumber==17){
    if(runNo>=270822 && runNo<=270830)period = 0;//17e
    if(runNo>=270854 && runNo<=270865)period = 1;//17f
    if(runNo>=271868 && runNo<=273103)period = 2;//17h
    if(runNo>=273591  && runNo<=274442)period = 3;//17i
    if(runNo>=274593  && runNo<=274671)period = 4;//17j 
    if(runNo>=274690  && runNo<=276508)period = 5;//17k
    if(runNo>=276551  && runNo<=278216)period = 6;//17l
    if(runNo>=278914  && runNo<=280140)period = 7;//17m
    if(runNo>=280282   && runNo<=281961)period = 8;//17o
    if(runNo>=282504  && runNo<=282704)period = 9;//17r
  }else if(fYearNumber==18){     
    if(runNo>=285008 && runNo<=285447)period = 0;//18b
    if(runNo>=285978 && runNo<=286350)period = 1;//18d
    if(runNo>=286380 && runNo<=286937)period = 2;//18e
    if(runNo>=287000  && runNo<=287977)period = 3;//18f
    if(runNo>=288619  && runNo<=288750)period = 4;//18g
    if(runNo>=288804  && runNo<=288806)period = 5;//18h
    if(runNo>=288861  && runNo<=288909 )period = 6;//18i
    if(runNo==288943)period = 7;//18j
    if(runNo>=289165   && runNo<=289201)period = 8;//18k
    if(runNo>=289240  && runNo<=289971)period = 9;//18l
    if(runNo>=290222  && runNo<=292839)period = 10;//18m
    if(runNo>=293357   && runNo<=293359)period = 11;//18n
    if(runNo>=293368   && runNo<=293898)period = 12;//18o
    if(runNo>=294009  && runNo<=294925)period = 13;//18p
  }
  return fMultEstimatorAvg[period];

}
//________________________________________________________________________
void AliAnalysisTaskSESigmacTopK0Spi::PrepareTracks(AliAODEvent *aod){
  
  // SELECT TRACKS and flag them: could consider to include common tracks (dca) to reduce cpu time (call once propagatetodca)
  ftrackArraySelSoftPi->Reset();
  fnSelSoftPi=0;
  
  for(Int_t itrack = 0; itrack < aod->GetNumberOfTracks(); itrack++){

    Int_t iSelSoftPionCuts = -1;

    AliAODTrack *track = dynamic_cast<AliAODTrack*>(aod->GetTrack(itrack));
    if(!track) AliFatal("Not a standard AOD");

    //Printf("selecting track");
    AliESDtrack *trackESD = SelectTrack(track, iSelSoftPionCuts, fESDtrackCutsSoftPion); //select soft-pion candidates
    
    if(!trackESD) continue;
    //Printf("good track");    
    
    if(iSelSoftPionCuts>=0){
      ftrackArraySelSoftPi->AddAt(itrack, fnSelSoftPi);
      fnSelSoftPi++;
    }    
    
    delete trackESD;
  }
}
///_______________________________________________________________________________________
AliESDtrack* AliAnalysisTaskSESigmacTopK0Spi::SelectTrack(AliAODTrack *aodtr, Int_t &isSelSoftPion, AliESDtrackCuts *cutsSoftPion){
  
  isSelSoftPion = -1;

  if(aodtr->GetID() < 0) return 0x0;

  if(TMath::Abs(aodtr->Charge()) != 1) return 0x0;

  AliESDtrack *esdTrack=new AliESDtrack(aodtr);
  // set the TPC cluster info
  esdTrack->SetTPCClusterMap(aodtr->GetTPCClusterMap());
  esdTrack->SetTPCSharedMap(aodtr->GetTPCSharedMap());
  esdTrack->SetTPCPointsF(aodtr->GetTPCNclsF());
  // needed to calculate the impact parameters
  Double_t pos[3], cov[6];
  fVtx1->GetXYZ(pos);
  fVtx1->GetCovarianceMatrix(cov);
  const AliESDVertex vESD(pos, cov, 100., 100);

  esdTrack->RelateToVertex(&vESD, 0., 3.);
  if(cutsSoftPion){
    if(cutsSoftPion->IsSelected(esdTrack)){
      isSelSoftPion = 0;
    }
  }

  if(isSelSoftPion < 0){
    delete esdTrack;
    return 0x0;   
  }
  
  return esdTrack;
}
//________________________________________________________________________
Double_t AliAnalysisTaskSESigmacTopK0Spi::CosThetaStar(Double_t mumVector[3],Double_t daughtVector[3],Double_t massMum,Double_t massDaught){
  
  Double_t mumP2=mumVector[0]*mumVector[0]+mumVector[1]*mumVector[1]+mumVector[2]*mumVector[2];
  Double_t mumP=TMath::Sqrt(mumP2);
  Double_t eMum=TMath::Sqrt(mumP2+massMum*massMum);
  Double_t daughtP2=daughtVector[0]*daughtVector[0]+daughtVector[1]*daughtVector[1]+daughtVector[2]*daughtVector[2];
  Double_t eDaugh=TMath::Sqrt(daughtP2+massDaught*massDaught);
  Double_t plLab=(mumVector[0]*daughtVector[0]+mumVector[1]*daughtVector[1]+mumVector[2]*daughtVector[2])/mumP;
  Double_t beta = mumP/eMum;
  Double_t gamma = eMum/massMum;
  Double_t plStar=gamma*(plLab-beta*eDaugh);
  Double_t daughtpT2=daughtP2-plLab*plLab;
  return plStar/TMath::Sqrt(plStar*plStar+daughtpT2);  
  
}
//______________________________________________________
void AliAnalysisTaskSESigmacTopK0Spi::LoopOverGenParticles(TClonesArray *mcArray){
  
  for(Int_t kmc = 0; kmc < mcArray->GetEntries(); kmc++){
    AliAODMCParticle *mcpart = (AliAODMCParticle*)mcArray->At(kmc);
    
    Int_t pdg = mcpart->GetPdgCode();
    Int_t arrayDauLab[3];
    Double_t pointLcSc[6];
    Double_t ptpartSc;
    Double_t ypartSc;
    
    if(TMath::Abs(pdg) == 4122){

      if(AliVertexingHFUtils::CheckLcV0bachelorDecay(mcArray, mcpart, arrayDauLab) >= 1) {
	Int_t checkOrigin = AliVertexingHFUtils::CheckOrigin(mcArray, mcpart, kTRUE);
	if(checkOrigin == 0)continue;
	
	Double_t ptpart = mcpart->Pt();
	Double_t ypart = mcpart->Y();
	
	Bool_t isFromSigmaC = kFALSE;
	Int_t indSc = mcpart->GetMother();
	AliAODMCParticle *mcpartMum = 0x0;
	if(indSc >= 0){
	  mcpartMum = (AliAODMCParticle*)mcArray->At(indSc); 
	  Int_t pdgLcMum = TMath::Abs(mcpartMum->GetPdgCode());
	  if(pdgLcMum == 4112 || pdgLcMum == 4222) isFromSigmaC = kTRUE;
	}
	
	if(isFromSigmaC){
	  ptpartSc = mcpartMum->Pt();
	  ypartSc = mcpartMum->Y();
	  pointLcSc[0] = ptpart;
	  pointLcSc[1] = -1;
	  pointLcSc[2] = checkOrigin;
	  pointLcSc[3] = ypart;
	  pointLcSc[4] = ptpartSc;
	  pointLcSc[5] = ypartSc;
	}
	if(TMath::Abs(ypart) < 0.5){
	  fhistMCSpectrumAccLc->Fill(ptpart, kGenLimAcc, checkOrigin);// Gen Level
	  
	  if(isFromSigmaC){
	    pointLcSc[1] = kGenLimAcc;
	    fhistMCSpectrumAccLcFromSc->Fill(pointLcSc);
	  }
	}
	
	Bool_t isInAcc = kTRUE;
	// check GenAcc level
	if(fAnalCuts){
	  if(!fAnalCuts->IsInFiducialAcceptance( ptpart, ypart)){
	    isInAcc = kFALSE;
	  }
	}
	else {
	  if(TMath::Abs(ypart) > 0.8){
	    isInAcc = kFALSE;
	  }
	}
	if(isInAcc){
	  fhistMCSpectrumAccLc->Fill(ptpart, kGenAccMother, checkOrigin);// Gen Acc Mother
	  if(isFromSigmaC){
	    pointLcSc[1] = kGenAccMother;
	    fhistMCSpectrumAccLcFromSc->Fill(pointLcSc);
	  }
	  
	  for(Int_t k = 0; k < 2; k++){
	    AliAODMCParticle *mcpartdau=(AliAODMCParticle*)mcArray->At(arrayDauLab[k]);
	    if(TMath::Abs(mcpartdau->Eta()) > 0.9){
	      isInAcc = kFALSE;
	    }	    
	  }
	  if(isInAcc){
	    fhistMCSpectrumAccLc->Fill(mcpart->Pt(), kGenAcc, checkOrigin);// Gen Acc
	    if(isFromSigmaC){
	      pointLcSc[1] = kGenAcc;
	      fhistMCSpectrumAccLcFromSc->Fill(pointLcSc);
	    }
	  }
	}
	
	if(isFromSigmaC){
	  // LimAcc level
	  if(TMath::Abs(ypartSc) < 0.5){
	    fhistMCSpectrumAccSc->Fill(ptpartSc, kGenLimAcc, checkOrigin);// Gen Level
	  }
	  // check GenAcc level
	  Bool_t isInAccSc = kTRUE;
	  if(fAnalCuts){
	    if(!fAnalCuts->IsInFiducialAcceptance(ptpartSc, ypartSc)){
	      isInAccSc = kFALSE;
	    }
	  }
	  else {
	    if(TMath::Abs(mcpartMum->Y()) > 0.8){
	      isInAccSc = kFALSE;
	    }
	  }
	  if(isInAccSc){
	    fhistMCSpectrumAccSc->Fill(ptpartSc, kGenAccMother, checkOrigin);// Gen Acc Mother
	    
	    if(isInAcc){// both Sc and Lc in fiducial acceptance + Lc daughter in Acc
	      for(Int_t k = mcpartMum->GetDaughterLabel(0); k < mcpartMum->GetDaughterLabel(1); k++){
		if(k >= 0){
		  AliAODMCParticle *mcpartMumdau = (AliAODMCParticle*)mcArray->At(k);
		  if(TMath::Abs(mcpartMumdau->GetPdgCode() == 211) && TMath::Abs(mcpartMumdau->Eta()) > 0.9){
		    isInAccSc = kFALSE;
		  }	    
		}
	      }
	      if(isInAccSc){
		fhistMCSpectrumAccSc->Fill(ptpartSc, kGenAcc, checkOrigin);// Gen Acc
	      }		
	    }
	  }	    
	}
      }
    }
  }
}

