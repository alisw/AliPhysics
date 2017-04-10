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

/* $Id: AliAnalysisTaskSELc2V0bachelorTMVA.cxx 62231 2013-04-29 09:47:26Z fprino $ */

//
//
//               Base class for Lc2V0 Analysis to be used with TMVA
//
//

//------------------------------------------------------------------------------------------
//
//  Author: C. Zampolli (from AliAnalysisTaskSELc2V0bachelor class by A.De Caro, P. Pagano)
//
//------------------------------------------------------------------------------------------

#include <TSystem.h>
#include <TParticle.h>
#include <TParticlePDG.h>
#include <TH1F.h>
#include <TH1I.h>
#include <TH2F.h>
#include <TTree.h>
#include "TROOT.h"
#include <TDatabasePDG.h>
#include "TVector3.h"

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
#include "AliAnalysisTaskSELc2V0bachelorTMVA.h"
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
#include "AliESDtrack.h"

using std::cout;
using std::endl;

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskSELc2V0bachelorTMVA);
/// \endcond

//__________________________________________________________________________
AliAnalysisTaskSELc2V0bachelorTMVA::AliAnalysisTaskSELc2V0bachelorTMVA():
AliAnalysisTaskSE(),
  fUseMCInfo(kFALSE),
  fOutput(0),
  fCEvents(0),
  fPIDResponse(0),
  fPIDCombined(0),
  fIsK0sAnalysis(kFALSE),
  fCounter(0),
  fAnalCuts(0),
  fListCuts(0),
  fListWeight(0),
  fHistoMCNch(0x0),
  fUseOnTheFlyV0(kFALSE),
  fIsEventSelected(kFALSE),
  fVariablesTreeSgn(0),
  fVariablesTreeBkg(0),
  fCandidateVariables(),
  fIspA(kFALSE),
  fHistoEvents(0),
  fHistoLc(0),
  fHistoLcOnTheFly(0),
  fFillOnlySgn(kFALSE),
  fHistoLcBeforeCuts(0),
  fHistoFiducialAcceptance(0),
  fHistoCodesSgn(0),
  fHistoCodesBkg(0),
  fHistoLcpKpiBeforeCuts(0),
  fVtx1(0),

  fHistoDistanceLcToPrimVtx(0),
  fHistoDistanceV0ToPrimVtx(0),
  fHistoDistanceV0ToLc(0),

  fHistoDistanceLcToPrimVtxSgn(0),
  fHistoDistanceV0ToPrimVtxSgn(0),
  fHistoDistanceV0ToLcSgn(0),

  fHistoVtxLcResidualToPrimVtx(0),
  fHistoVtxV0ResidualToPrimVtx(0),
  fHistoVtxV0ResidualToLc(0),

  fHistoMassV0All(0),
  fHistoDecayLengthV0All(0),
  fHistoLifeTimeV0All(0),

  fHistoMassV0True(0),
  fHistoDecayLengthV0True(0),
  fHistoLifeTimeV0True(0),

  fHistoMassV0TrueFromAOD(0),

  fHistoMassV0TrueK0S(0),
  fHistoDecayLengthV0TrueK0S(0),
  fHistoLifeTimeV0TrueK0S(0),

  fHistoMassV0TrueK0SFromAOD(0),

  fHistoMassLcAll(0),
  fHistoDecayLengthLcAll(0),
  fHistoLifeTimeLcAll(0),

  fHistoMassLcTrue(0),
  fHistoDecayLengthLcTrue(0),
  fHistoLifeTimeLcTrue(0),

  fHistoMassLcTrueFromAOD(0),

  fHistoMassV0fromLcAll(0),
  fHistoDecayLengthV0fromLcAll(0),
  fHistoLifeTimeV0fromLcAll(0),

  fHistoMassV0fromLcTrue(0),
  fHistoDecayLengthV0fromLcTrue(0),
  fHistoLifeTimeV0fromLcTrue(0),

  fHistoMassLcSgn(0),
  fHistoMassLcSgnFromAOD(0),
  fHistoDecayLengthLcSgn(0),
  fHistoLifeTimeLcSgn(0),

  fHistoMassV0fromLcSgn(0),
  fHistoDecayLengthV0fromLcSgn(0),
  fHistoLifeTimeV0fromLcSgn(0),

  fHistoKF(0),
  fHistoKFV0(0),
  fHistoKFLc(0),

  fHistoMassKFV0(0),
  fHistoDecayLengthKFV0(0),
  fHistoLifeTimeKFV0(0),

  fHistoMassKFLc(0),
  fHistoDecayLengthKFLc(0),
  fHistoLifeTimeKFLc(0),

  fHistoArmenterosPodolanskiV0KF(0),
  fHistoArmenterosPodolanskiV0KFSgn(0),
  fHistoArmenterosPodolanskiV0AOD(0),
  fHistoArmenterosPodolanskiV0AODSgn(0),

  fOutputKF(0),
  fmcLabelLc(-1),
  ftopoConstraint(kTRUE),
  fCallKFVertexing(kFALSE),
  fKeepingOnlyHIJINGBkg(kFALSE),
  fUtils(0),
  fHistoBackground(0),
  fCutKFChi2NDF(999999.),
  fCutKFDeviationFromVtx(999999.),
  fCutKFDeviationFromVtxV0(0.),
  fCurrentEvent(-1),
  fBField(0),
  fKeepingOnlyPYTHIABkg(kFALSE),
  fHistoMCLcK0SpGen(0x0),
  fHistoMCLcK0SpGenAcc(0x0),
  fHistoMCLcK0SpGenLimAcc(0x0),
  fTriggerMask(0),
  fFuncWeightPythia(0),
  fFuncWeightFONLL5overLHC13d3(0),
  fFuncWeightFONLL5overLHC13d3Lc(0),
  fNTracklets(0)
{
  //
  /// Default ctor
  //
}
//___________________________________________________________________________
AliAnalysisTaskSELc2V0bachelorTMVA::AliAnalysisTaskSELc2V0bachelorTMVA(const Char_t* name,
								       AliRDHFCutsLctoV0* analCuts, Bool_t useOnTheFly) :
  AliAnalysisTaskSE(name),
  fUseMCInfo(kFALSE),
  fOutput(0),
  fCEvents(0),
  fPIDResponse(0),
  fPIDCombined(0),
  fIsK0sAnalysis(kFALSE),
  fCounter(0),
  fAnalCuts(analCuts),
  fListCuts(0),
  fListWeight(0),
  fHistoMCNch(0x0),
  fUseOnTheFlyV0(useOnTheFly),
  fIsEventSelected(kFALSE),
  fVariablesTreeSgn(0),
  fVariablesTreeBkg(0),
  fCandidateVariables(),
  fIspA(kFALSE),
  fHistoEvents(0),
  fHistoLc(0),
  fHistoLcOnTheFly(0),
  fFillOnlySgn(kFALSE),
  fHistoLcBeforeCuts(0),
  fHistoFiducialAcceptance(0),
  fHistoCodesSgn(0),
  fHistoCodesBkg(0),
  fHistoLcpKpiBeforeCuts(0),
  fVtx1(0),

  fHistoDistanceLcToPrimVtx(0),
  fHistoDistanceV0ToPrimVtx(0),
  fHistoDistanceV0ToLc(0),

  fHistoDistanceLcToPrimVtxSgn(0),
  fHistoDistanceV0ToPrimVtxSgn(0),
  fHistoDistanceV0ToLcSgn(0),

  fHistoVtxLcResidualToPrimVtx(0),
  fHistoVtxV0ResidualToPrimVtx(0),
  fHistoVtxV0ResidualToLc(0),

  fHistoMassV0All(0),
  fHistoDecayLengthV0All(0),
  fHistoLifeTimeV0All(0),

  fHistoMassV0True(0),
  fHistoDecayLengthV0True(0),
  fHistoLifeTimeV0True(0),

  fHistoMassV0TrueFromAOD(0),

  fHistoMassV0TrueK0S(0),
  fHistoDecayLengthV0TrueK0S(0),
  fHistoLifeTimeV0TrueK0S(0),

  fHistoMassV0TrueK0SFromAOD(0),

  fHistoMassLcAll(0),
  fHistoDecayLengthLcAll(0),
  fHistoLifeTimeLcAll(0),

  fHistoMassLcTrue(0),
  fHistoDecayLengthLcTrue(0),
  fHistoLifeTimeLcTrue(0),

  fHistoMassLcTrueFromAOD(0),

  fHistoMassV0fromLcAll(0),
  fHistoDecayLengthV0fromLcAll(0),
  fHistoLifeTimeV0fromLcAll(0),

  fHistoMassV0fromLcTrue(0),
  fHistoDecayLengthV0fromLcTrue(0),
  fHistoLifeTimeV0fromLcTrue(0),

  fHistoMassLcSgn(0),
  fHistoMassLcSgnFromAOD(0),
  fHistoDecayLengthLcSgn(0),
  fHistoLifeTimeLcSgn(0),

  fHistoMassV0fromLcSgn(0),
  fHistoDecayLengthV0fromLcSgn(0),
  fHistoLifeTimeV0fromLcSgn(0),

  fHistoKF(0),
  fHistoKFV0(0),
  fHistoKFLc(0),

  fHistoMassKFV0(0),
  fHistoDecayLengthKFV0(0),
  fHistoLifeTimeKFV0(0),

  fHistoMassKFLc(0),
  fHistoDecayLengthKFLc(0),
  fHistoLifeTimeKFLc(0),

  fHistoArmenterosPodolanskiV0KF(0),
  fHistoArmenterosPodolanskiV0KFSgn(0),
  fHistoArmenterosPodolanskiV0AOD(0),
  fHistoArmenterosPodolanskiV0AODSgn(0),

  fOutputKF(0),
  fmcLabelLc(-1),
  ftopoConstraint(kTRUE),
  fCallKFVertexing(kFALSE),
  fKeepingOnlyHIJINGBkg(kFALSE),
  fUtils(0),
  fHistoBackground(0),
  fCutKFChi2NDF(999999.),
  fCutKFDeviationFromVtx(999999.),
  fCutKFDeviationFromVtxV0(0.),
  fCurrentEvent(-1),
  fBField(0),
  fKeepingOnlyPYTHIABkg(kFALSE),
  fHistoMCLcK0SpGen(0x0),
  fHistoMCLcK0SpGenAcc(0x0),
  fHistoMCLcK0SpGenLimAcc(0x0),
  fTriggerMask(0),

  fFuncWeightPythia(0),
  fFuncWeightFONLL5overLHC13d3(0),
  fFuncWeightFONLL5overLHC13d3Lc(0),
  fNTracklets(0)
{
  //
  /// Constructor. Initialization of Inputs and Outputs
  //
  Info("AliAnalysisTaskSELc2V0bachelorTMVA","Calling Constructor");

  DefineOutput(1, TList::Class());  // Tree signal + Tree Bkg + histoEvents
  DefineOutput(2, AliNormalizationCounter::Class()); // normalization counter object
  DefineOutput(3, TList::Class());  // Cuts
  DefineOutput(4, TTree::Class());  // Tree signal + Tree Bkg + histoEvents
  DefineOutput(5, TTree::Class());  // Tree signal + Tree Bkg + histoEvents
  DefineOutput(6, TList::Class());  // Tree signal + Tree Bkg + histoEvents
  DefineOutput(7, TList::Class());  // weights
}

//___________________________________________________________________________
AliAnalysisTaskSELc2V0bachelorTMVA::~AliAnalysisTaskSELc2V0bachelorTMVA() {
  //
  /// destructor
  //
  Info("~AliAnalysisTaskSELc2V0bachelorTMVA","Calling Destructor");

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

  if (fAnalCuts) {
    delete fAnalCuts;
    fAnalCuts = 0;
  }

  if (fListCuts) {
    delete fListCuts;
    fListCuts = 0;
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

  if (fOutputKF) {
    delete fOutputKF;
    fOutputKF = 0;
  }

  if (fUtils) {
    delete fUtils;
    fUtils = 0;
  }

}
//_________________________________________________
void AliAnalysisTaskSELc2V0bachelorTMVA::Init() {
  //
  /// Initialization
  //

  fIsEventSelected=kFALSE;

  if (fDebug > 1) AliInfo("Init");

  fListCuts = new TList();
  fListCuts->SetOwner();
  fListCuts->Add(new AliRDHFCutsLctoV0(*fAnalCuts));
  PostData(3,fListCuts);

  // Save the weight functions or histograms
  fListWeight = new TList();
  fListWeight->Add(fHistoMCNch);
  PostData(7,fListWeight);

  if (fUseMCInfo && (fKeepingOnlyHIJINGBkg || fKeepingOnlyPYTHIABkg)) fUtils = new AliVertexingHFUtils();

  return;
}

//________________________________________ terminate ___________________________
void AliAnalysisTaskSELc2V0bachelorTMVA::Terminate(Option_t*)
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


  //AliDebug(2, Form("At MC level, %f Lc --> K0S + p were found", fHistoMCLcK0SpGen->GetEntries()));
  //AliDebug(2, Form("At MC level, %f Lc --> K0S + p were found in the acceptance", fHistoMCLcK0SpGenAcc->GetEntries()));
  //AliDebug(2, Form("At Reco level, %lld Lc --> K0S + p were found", fVariablesTreeSgn->GetEntries()));
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

  fOutputKF = dynamic_cast<TList*> (GetOutputData(6));
  if (!fOutputKF) {
    AliError("fOutputKF not available");
    return;
  }

  return;
}

//___________________________________________________________________________
void AliAnalysisTaskSELc2V0bachelorTMVA::UserCreateOutputObjects() {
  /// output
  AliInfo(Form("CreateOutputObjects of task %s\n", GetName()));

  //slot #1
  //OpenFile(1);
  fOutput = new TList();
  fOutput->SetOwner();
  fOutput->SetName("listTrees");

  // Output slot 1: list of 2 trees (Sgn + Bkg) of the candidate variables
  const char* nameoutput = GetOutputSlot(1)->GetContainer()->GetName();
  fVariablesTreeSgn = new TTree(Form("%s_Sgn", nameoutput), "Candidates variables tree, Signal");
  fVariablesTreeBkg = new TTree(Form("%s_Bkg", nameoutput), "Candidates variables tree, Background");
  Int_t nVar = 90;
  fCandidateVariables = new Float_t [nVar];
  TString * fCandidateVariableNames = new TString[nVar];
  fCandidateVariableNames[0]="massLc2K0Sp";
  fCandidateVariableNames[1]="massLc2Lambdapi";
  fCandidateVariableNames[2]="massK0S";
  fCandidateVariableNames[3]="massLambda";
  fCandidateVariableNames[4]="massLambdaBar";
  fCandidateVariableNames[5]="cosPAK0S";
  fCandidateVariableNames[6]="dcaV0";
  fCandidateVariableNames[7]="tImpParBach";
  fCandidateVariableNames[8]="tImpParV0";
  fCandidateVariableNames[9]="nSigmaTPCpr";
  fCandidateVariableNames[10]="nSigmaTPCpi";
  fCandidateVariableNames[11]="nSigmaTPCka";
  fCandidateVariableNames[12]="nSigmaTOFpr";
  fCandidateVariableNames[13]="nSigmaTOFpi";
  fCandidateVariableNames[14]="nSigmaTOFka";
  fCandidateVariableNames[15]="bachelorPt";
  fCandidateVariableNames[16]="V0positivePt";
  fCandidateVariableNames[17]="V0negativePt";
  fCandidateVariableNames[18]="dcaV0pos";
  fCandidateVariableNames[19]="dcaV0neg";
  fCandidateVariableNames[20]="v0Pt";
  fCandidateVariableNames[21]="massGamma";
  fCandidateVariableNames[22]="LcPt";
  fCandidateVariableNames[23]="combinedProtonProb";
  fCandidateVariableNames[24]="LcEta";
  fCandidateVariableNames[25]="V0positiveEta";
  fCandidateVariableNames[26]="V0negativeEta";
  fCandidateVariableNames[27]="TPCProtonProb";
  fCandidateVariableNames[28]="TOFProtonProb";
  fCandidateVariableNames[29]="bachelorEta";
  fCandidateVariableNames[30]="LcP";
  fCandidateVariableNames[31]="bachelorP";
  fCandidateVariableNames[32]="v0P";
  fCandidateVariableNames[33]="V0positiveP";
  fCandidateVariableNames[34]="V0negativeP";
  fCandidateVariableNames[35]="LcY";
  fCandidateVariableNames[36]="v0Y";
  fCandidateVariableNames[37]="bachelorY";
  fCandidateVariableNames[38]="V0positiveY";
  fCandidateVariableNames[39]="V0negativeY";
  fCandidateVariableNames[40]="v0Eta";
  fCandidateVariableNames[41]="DecayLengthLc";
  fCandidateVariableNames[42]="DecayLengthK0S";
  fCandidateVariableNames[43]="CtLc";
  fCandidateVariableNames[44]="CtK0S";
  fCandidateVariableNames[45]="bachCode";
  fCandidateVariableNames[46]="k0SCode";

  fCandidateVariableNames[47]="V0KFmass";
  fCandidateVariableNames[48]="V0KFdecayLength";
  fCandidateVariableNames[49]="V0KFlifeTime";

  fCandidateVariableNames[50]="V0KFmassErr";
  fCandidateVariableNames[51]="V0KFdecayTimeErr";
  fCandidateVariableNames[52]="V0KFlifeTimeErr";

  fCandidateVariableNames[53]="LcKFmass";
  fCandidateVariableNames[54]="LcKFdecayLength";
  fCandidateVariableNames[55]="LcKFlifeTime";

  fCandidateVariableNames[56]="LcKFmassErr";
  fCandidateVariableNames[57]="LcKFdecayTimeErr";
  fCandidateVariableNames[58]="LcKFlifeTimeErr";

  fCandidateVariableNames[59]="LcKFDistToPrimVtx";
  fCandidateVariableNames[60]="V0KFDistToPrimVtx";
  fCandidateVariableNames[61]="V0KFDistToLc";
  fCandidateVariableNames[62]="alphaArmKF";
  fCandidateVariableNames[63]="ptArmKF";
  fCandidateVariableNames[64]="alphaArm";
  fCandidateVariableNames[65]="ptArm";

  fCandidateVariableNames[66]="ITSrefitV0pos";
  fCandidateVariableNames[67]="ITSrefitV0neg";

  fCandidateVariableNames[68]="TPCClV0pos";
  fCandidateVariableNames[69]="TPCClV0neg";

  fCandidateVariableNames[70]="v0Xcoord";
  fCandidateVariableNames[71]="v0Ycoord";
  fCandidateVariableNames[72]="v0Zcoord";
  fCandidateVariableNames[73]="primVtxX";
  fCandidateVariableNames[74]="primVtxY";
  fCandidateVariableNames[75]="primVtxZ";

  fCandidateVariableNames[76]="ITSclBach";
  fCandidateVariableNames[77]="SPDclBach";

  fCandidateVariableNames[78]="ITSclV0pos";
  fCandidateVariableNames[79]="SPDclV0pos";
  fCandidateVariableNames[80]="ITSclV0neg";
  fCandidateVariableNames[81]="SPDclV0neg";

  fCandidateVariableNames[82]="alphaArmLc";
  fCandidateVariableNames[83]="alphaArmLcCharge";
  fCandidateVariableNames[84]="ptArmLc";

  fCandidateVariableNames[85]="CosThetaStar";

  fCandidateVariableNames[86]="weightPtFlat";
  fCandidateVariableNames[87]="weightFONLL5overLHC13d3";
  fCandidateVariableNames[88]="weightFONLL5overLHC13d3Lc";
  fCandidateVariableNames[89]="weightNch";


  for(Int_t ivar=0; ivar<nVar; ivar++){
    fVariablesTreeSgn->Branch(fCandidateVariableNames[ivar].Data(),&fCandidateVariables[ivar],Form("%s/f",fCandidateVariableNames[ivar].Data()));
    fVariablesTreeBkg->Branch(fCandidateVariableNames[ivar].Data(),&fCandidateVariables[ivar],Form("%s/f",fCandidateVariableNames[ivar].Data()));
  }

  fHistoEvents = new TH1F("fHistoEvents", "fHistoEvents", 2, -0.5, 1.5);
  TString labelEv[2] = {"NotSelected", "Selected"};
  for (Int_t ibin = 1; ibin <= fHistoEvents->GetNbinsX(); ibin++){
    fHistoEvents->GetXaxis()->SetBinLabel(ibin, labelEv[ibin-1].Data());
  }

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

  fHistoLcpKpiBeforeCuts = new TH1F("fHistoLcpKpiBeforeCuts", "fHistoLcpKpiBeforeCuts", 2, -0.5, 1.5);
  for (Int_t ibin = 1; ibin <= fHistoLcpKpiBeforeCuts->GetNbinsX(); ibin++){
    fHistoLcpKpiBeforeCuts->GetXaxis()->SetBinLabel(ibin, labelBeforeCuts[ibin-1].Data());
  }

  fHistoBackground = new TH1F("fHistoBackground", "fHistoBackground", 4, -0.5, 3.5);
  TString labelBkg[4] = {"Injected", "Non-injected", "Non-PYTHIA", "PYTHIA"};
  for (Int_t ibin = 1; ibin <= fHistoBackground->GetNbinsX(); ibin++){
    fHistoBackground->GetXaxis()->SetBinLabel(ibin, labelBkg[ibin-1].Data());
  }

  //fOutput->Add(fVariablesTreeSgn);
  //fOutput->Add(fVariablesTreeBkg);

  const Float_t ptbins[15] = {0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 12., 17., 25., 35.};

  fHistoMCLcK0SpGen = new TH1F("fHistoMCLcK0SpGen", "fHistoMCLcK0SpGen", 14, ptbins);
  fHistoMCLcK0SpGenAcc = new TH1F("fHistoMCLcK0SpGenAcc", "fHistoMCLcK0SpGenAcc", 14, ptbins);
  fHistoMCLcK0SpGenLimAcc = new TH1F("fHistoMCLcK0SpGenLimAcc", "fHistoMCLcK0SpGenLimAcc", 14, ptbins);

  fOutput->Add(fHistoEvents);
  fOutput->Add(fHistoLc);
  fOutput->Add(fHistoLcOnTheFly);
  fOutput->Add(fHistoLcBeforeCuts);
  fOutput->Add(fHistoFiducialAcceptance);
  fOutput->Add(fHistoCodesSgn);
  fOutput->Add(fHistoCodesBkg);
  fOutput->Add(fHistoLcpKpiBeforeCuts);
  fOutput->Add(fHistoBackground);
  fOutput->Add(fHistoMCLcK0SpGen);
  fOutput->Add(fHistoMCLcK0SpGenAcc);
  fOutput->Add(fHistoMCLcK0SpGenLimAcc);

  PostData(1, fOutput);
  PostData(4, fVariablesTreeSgn);
  PostData(5, fVariablesTreeBkg);
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  fPIDResponse = inputHandler->GetPIDResponse();

  if (fAnalCuts->GetIsUsePID()){
    /*
      fAnalCuts->GetPidHF()->SetPidResponse(fPIDResponse);
      fAnalCuts->GetPidV0pos()->SetPidResponse(fPIDResponse);
      fAnalCuts->GetPidV0neg()->SetPidResponse(fPIDResponse);
      fAnalCuts->GetPidHF()->SetOldPid(kFALSE);
      fAnalCuts->GetPidV0pos()->SetOldPid(kFALSE);
      fAnalCuts->GetPidV0neg()->SetOldPid(kFALSE);
    */
    fAnalCuts->SetUsePID(kFALSE); // I don't want to use the PID through the cut object, but I will use the PID response directly!!!
  }

  // Setting properties of PID
  fPIDCombined=new AliPIDCombined;
  fPIDCombined->SetDefaultTPCPriors();
  fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTPC+AliPIDResponse::kDetTOF);
  //fPIDCombined->SetPriorDistribution((AliPID::EParticleType)ispec,fPriors[ispec]);

  fCounter = new AliNormalizationCounter("NormalizationCounter");
  fCounter->Init();
  PostData(2, fCounter);

  // Histograms from KF

  if (fCallKFVertexing){
    fHistoDistanceLcToPrimVtx    = new TH1D("fHistoDistanceLcToPrimVtx", "Lc distance to Prim Vertex from KF; distance [cm]", 1000, 0., 1);
    fHistoDistanceV0ToPrimVtx    = new TH1D("fHistoDistanceV0ToPrimVtx", "V0 distance to Prim Vertex from KF; distance [cm]", 1000, 0., 100.);
    fHistoDistanceV0ToLc         = new TH1D("fHistoDistanceV0ToLc", "V0 distance to Lc Vertex from KF; distance [cm]", 1000, 0., 100.);

    fHistoDistanceLcToPrimVtxSgn = new TH1D("fHistoDistanceLcToPrimVtxSgn", "Lc Sgn distance to Prim Vertex from KF; distance [cm]", 1000, 0., 1);
    fHistoDistanceV0ToPrimVtxSgn = new TH1D("fHistoDistanceV0ToPrimVtxSgn", "V0 Sgn distance to Prim Vertex from KF; distance [cm]", 1000, 0., 100.);
    fHistoDistanceV0ToLcSgn      = new TH1D("fHistoDistanceV0ToLcSgn", "V0 Sgn distance to Lc Vertex from KF; distance [cm]", 1000, 0., 100.);

    fHistoVtxLcResidualToPrimVtx = new TH1D("fHistoVtxLcResidualToPrimVtx", "Residual between MC and KF (MC - KF): Lc to Prim Vtx; distance [cm]", 1000, -5., 5.);
    fHistoVtxV0ResidualToPrimVtx = new TH1D("fHistoVtxV0ResidualToPrimVtx", "Residual between MC and KF (MC - KF): V0 to Prim Vtx; distance [cm]", 1000, -5., 5.);
    fHistoVtxV0ResidualToLc      = new TH1D("fHistoVtxV0ResidualToLc", "Residual between MC and KF: V0 to Lc (MC - KF); distance [cm]", 1000, -5., 5.);

    fHistoMassV0All              = new TH1D("fHistoMassV0All", "V0 Mass; mass", 500, 0.4, 0.6);
    fHistoDecayLengthV0All       = new TH1D("fHistoDecayLengthV0All", "V0 Decay Length; decayLength", 500, -10, 10.0);
    fHistoLifeTimeV0All          = new TH1D("fHistoLifeTimeV0All", "V0 Life Time; lifeTime", 500, -10.0, 10.0);

    fHistoMassV0True             = new TH1D("fHistoMassV0True", "True V0 Mass; mass", 500, 0.4, 0.6);
    fHistoDecayLengthV0True      = new TH1D("fHistoDecayLengthV0True", "True V0 Decay Length; decayLength", 500, -10, 10.0);
    fHistoLifeTimeV0True         = new TH1D("fHistoLifeTimeV0True", "True V0 Life Time; lifeTime", 500, -10.0, 10.0);

    fHistoMassV0TrueFromAOD      = new TH1D("fHistoMassV0TrueFormAOD", "True V0 Mass (AOD); mass", 500, 0.4, 0.6);

    fHistoMassV0TrueK0S          = new TH1D("fHistoMassV0TrueK0S", "True V0-K0S Mass; mass", 500, 0.4, 0.6);
    fHistoDecayLengthV0TrueK0S   = new TH1D("fHistoDecayLengthV0TrueK0S", "True V0-K0S Decay Length; decayLength", 500, -10, 10.0);
    fHistoLifeTimeV0TrueK0S      = new TH1D("fHistoLifeTimeV0TrueK0S", "True V0-K0S Life Time; lifeTime", 500, -10.0, 10.0);

    fHistoMassV0TrueK0SFromAOD   = new TH1D("fHistoMassV0TrueK0SFormAOD", "True V0-K0S Mass (AOD); mass", 500, 0.4, 0.6);

    fHistoMassLcAll              = new TH1D("fHistoMassLcAll", "Lc Mass; mass", 500, 2.0, 3.0);
    fHistoDecayLengthLcAll       = new TH1D("fHistoDecayLengthLcAll", "Lc Decay Length; decayLenght", 100000, -0.1, 0.1);
    fHistoLifeTimeLcAll          = new TH1D("fHistoLifeTimeLcAll", "Lc Life Time; lifeTime", 100000, -0.1, 0.1);

    fHistoMassLcTrue             = new TH1D("fHistoMassLcTrue", "True Lc Mass; mass", 500, 2.0, 3.0);
    fHistoDecayLengthLcTrue      = new TH1D("fHistoDecayLengthLcTrue", "True Lc Decay Length; decayLength", 100000, -0.1, 0.1);
    fHistoLifeTimeLcTrue         = new TH1D("fHistoLifeTimeLcTrue", "True Lc Life Time; lifeTime", 100000, -0.1, 0.1);

    fHistoMassLcTrueFromAOD      = new TH1D("fHistoMassLcTrueFromAOD", "True Lc Mass (AOD); mass", 500, 2.0, 3.0);

    fHistoMassV0fromLcAll        = new TH1D("fHistoMassV0fromLcAll", "V0 mass from Lc built in KF; mass", 500, 0.4, 0.6);
    fHistoDecayLengthV0fromLcAll = new TH1D("fHistoDecayLengthV0fromLcAll", "V0 Decay Length from Lc built in KF; decayLength", 500, 0, 10.0);
    fHistoLifeTimeV0fromLcAll    = new TH1D("fHistoLifeTimeV0fromLcAll", "V0 Life Time from Lc built in KF; lifeTime", 500, 0.0, 3.0);

    fHistoMassV0fromLcTrue       = new TH1D("fHistoMassV0fromLcTrue", "V0 mass from true Lc built in KF; mass", 500, 0.4, 0.6);
    fHistoDecayLengthV0fromLcTrue= new TH1D("fHistoDecayLengthV0fromLcTrue", "V0 Decay Length from true Lc built in KF; decayLength", 500, 0, 10.0);
    fHistoLifeTimeV0fromLcTrue   = new TH1D("fHistoLifeTimeV0fromLcTrue", "V0 Life Time from true Lc built in KF; lifeTime", 500, 0.0, 3.0);

    fHistoMassLcSgn              = new TH1D("fHistoMassLcSgn", "True Lc Signal Mass; mass", 500, 2.0, 3.0);
    fHistoMassLcSgnFromAOD       = new TH1D("fHistoMassLcSgnFromAOD", "True Lc Signal Mass (AOD); mass", 500, 2.0, 3.0);
    fHistoDecayLengthLcSgn       = new TH1D("fHistoDecayLengthLcSgn", "True Lc Signal Decay Length; decayLength", 100000, -0.1, 0.1);
    fHistoLifeTimeLcSgn          = new TH1D("fHistoLifeTimeLcSgn", "True Lc Signal Life Time; lifeTime", 100000, -0.1, 0.1);

    fHistoMassV0fromLcSgn        = new TH1D("fHistoMassV0fromLcSgn", "V0 from True Lc Signal Mass; mass", 500, 0.4, 0.6);
    fHistoDecayLengthV0fromLcSgn = new TH1D("fHistoDecayLengthV0fromLcSgn", "V0 True Lc Signal Decay Length; decayLength", 500, 0, 10.0);
    fHistoLifeTimeV0fromLcSgn    = new TH1D("fHistoLifeTimeV0fromLcSgn", "V0 True Lc Signal Life Time; lifeTime", 500, 0.0, 3.0);

    fHistoKF = new TH2D("fHistoKF", "Summary from KF; V0 KF; Lc KF", 16, -0.5, 15.5, 16, -0.5, 15.5);
    fHistoKFV0 = new TH1D("fHistoKFV0", "Summary from KF; V0 KF", 16, -0.5, 15.5);
    fHistoKFLc = new TH1D("fHistoKFLc", "Summary from KF; V0 KF", 16, -0.5, 15.5);
    TString axisLabel[16] = {"AllOk", "M_NotOk", "Sm_NotOk", "Dl_NotOk",
			     "Lt_NotOk", "M_Sm_NotOk", "M_Dl_NotOk", "M_Lt_NotOk",
			     "Dl_Sm_NotOk", "Dl_Lt_NotOk", "Sm_Lt_NotOk", "M_Sm_Dl_NotOk",
			     "M_Sm_Lt_NotOk", "Sm_Dl_Lt_NotOk", "M_Dl_Lt_NotOk", "All_NotOk"};

    for (Int_t ibin = 1; ibin <=16; ibin++){
      fHistoKF->GetXaxis()->SetBinLabel(ibin, axisLabel[ibin-1].Data());
      fHistoKF->GetYaxis()->SetBinLabel(ibin, axisLabel[ibin-1].Data());
      fHistoKFV0->GetXaxis()->SetBinLabel(ibin, axisLabel[ibin-1].Data());
      fHistoKFLc->GetXaxis()->SetBinLabel(ibin, axisLabel[ibin-1].Data());
    }

    fHistoMassKFV0 = new TH2D("fHistoMassKFV0", "mass vs sigmaMass for V0; mass; sigmaMass", 500, 0.4, 0.6, 500, 0., 10);
    fHistoDecayLengthKFV0 = new TH2D("fHistoDecayLengthKFV0", "decayLength vs sigmaDecayLength for V0; decayLength; sigmaDecayLength", 500, -10, 10, 500, 0., 10);
    fHistoLifeTimeKFV0 = new TH2D("fHistoLifeTimeKFV0", "lifeTime vs sigmalifeTime for V0; lifeTime; sigmaLifeTime", 500, -10, 10, 500, 0., 10);

    fHistoMassKFLc = new TH2D("fHistoMassKFLc", "mass vs sigmaMass for Lc; mass; sigmaMass", 500, 0.4, 0.6, 500, 0., 10);
    fHistoDecayLengthKFLc = new TH2D("fHistoDecayLengthKFLc", "decayLength vs sigmaDecayLength for Lc; decayLength; sigmaDecayLength", 500, -10, 10, 500, 0., 10);
    fHistoLifeTimeKFLc = new TH2D("fHistoLifeTimeKFLc", "lifeTime vs sigmalifeTime for Lc; lifeTime; sigmaLifeTime", 500, -10, 10, 500, 0., 10);

    fHistoArmenterosPodolanskiV0KF = new TH2D("fHistoArmenterosPodolanskiV0KF", "V0 ArmenterosPodolanski from KF; #alpha; Qt", 1000, -1, 1, 1000, 0, 1);
    fHistoArmenterosPodolanskiV0KFSgn = new TH2D("fHistoArmenterosPodolanskiV0KFSgn", "V0 (signal) ArmenterosPodolanski from KF; #alpha; Qt", 1000, -1, 1, 1000, 0, 1);
    fHistoArmenterosPodolanskiV0AOD = new TH2D("fHistoArmenterosPodolanskiV0AOD", "V0 ArmenterosPodolanski from AOD; #alpha; Qt", 1000, -1, 1, 1000, 0, 1);
    fHistoArmenterosPodolanskiV0AODSgn = new TH2D("fHistoArmenterosPodolanskiV0AODSgn", "V0 (signal) ArmenterosPodolanski from AOD; #alpha; Qt", 1000, -1, 1, 1000, 0, 1);
  }

  fOutputKF = new TList();
  fOutputKF->SetOwner();
  fOutputKF->SetName("listHistoKF");

  if (fCallKFVertexing){
    fOutputKF->Add(fHistoDistanceLcToPrimVtx);
    fOutputKF->Add(fHistoDistanceV0ToPrimVtx);
    fOutputKF->Add(fHistoDistanceV0ToLc);

    fOutputKF->Add(fHistoDistanceLcToPrimVtxSgn);
    fOutputKF->Add(fHistoDistanceV0ToPrimVtxSgn);
    fOutputKF->Add(fHistoDistanceV0ToLcSgn);

    fOutputKF->Add(fHistoVtxLcResidualToPrimVtx);
    fOutputKF->Add(fHistoVtxV0ResidualToPrimVtx);
    fOutputKF->Add(fHistoVtxV0ResidualToLc);

    fOutputKF->Add(fHistoMassV0All);
    fOutputKF->Add(fHistoDecayLengthV0All);
    fOutputKF->Add(fHistoLifeTimeV0All);

    fOutputKF->Add(fHistoMassV0True);
    fOutputKF->Add(fHistoDecayLengthV0True);
    fOutputKF->Add(fHistoLifeTimeV0True);

    fOutputKF->Add(fHistoMassV0TrueFromAOD);

    fOutputKF->Add(fHistoMassV0TrueK0S);
    fOutputKF->Add(fHistoDecayLengthV0TrueK0S);
    fOutputKF->Add(fHistoLifeTimeV0TrueK0S);

    fOutputKF->Add(fHistoMassV0TrueK0SFromAOD);

    fOutputKF->Add(fHistoMassLcAll);
    fOutputKF->Add(fHistoDecayLengthLcAll);
    fOutputKF->Add(fHistoLifeTimeLcAll);

    fOutputKF->Add(fHistoMassLcTrue);
    fOutputKF->Add(fHistoDecayLengthLcTrue);
    fOutputKF->Add(fHistoLifeTimeLcTrue);

    fOutputKF->Add(fHistoMassLcTrueFromAOD);

    fOutputKF->Add(fHistoMassV0fromLcAll);
    fOutputKF->Add(fHistoDecayLengthV0fromLcAll);
    fOutputKF->Add(fHistoLifeTimeV0fromLcAll);

    fOutputKF->Add(fHistoMassV0fromLcTrue);
    fOutputKF->Add(fHistoDecayLengthV0fromLcTrue);
    fOutputKF->Add(fHistoLifeTimeV0fromLcTrue);

    fOutputKF->Add(fHistoMassLcSgn);
    fOutputKF->Add(fHistoMassLcSgnFromAOD);
    fOutputKF->Add(fHistoDecayLengthLcSgn);
    fOutputKF->Add(fHistoLifeTimeLcSgn);

    fOutputKF->Add(fHistoMassV0fromLcSgn);
    fOutputKF->Add(fHistoDecayLengthV0fromLcSgn);
    fOutputKF->Add(fHistoLifeTimeV0fromLcSgn);

    fOutputKF->Add(fHistoKF);
    fOutputKF->Add(fHistoKFV0);
    fOutputKF->Add(fHistoKFLc);

    fOutputKF->Add(fHistoMassKFV0);
    fOutputKF->Add(fHistoDecayLengthKFV0);
    fOutputKF->Add(fHistoLifeTimeKFV0);

    fOutputKF->Add(fHistoMassKFLc);
    fOutputKF->Add(fHistoDecayLengthKFLc);
    fOutputKF->Add(fHistoLifeTimeKFLc);

    fOutputKF->Add(fHistoArmenterosPodolanskiV0KF);
    fOutputKF->Add(fHistoArmenterosPodolanskiV0KFSgn);
    fOutputKF->Add(fHistoArmenterosPodolanskiV0AOD);
    fOutputKF->Add(fHistoArmenterosPodolanskiV0AODSgn);
  }

  // weight function from ratio of flat value (1/30) to pythia
  // use to normalise to flat distribution (should lead to flat pT distribution)
  fFuncWeightPythia = new TF1("funcWeightPythia","1./(30. *[0]*x/TMath::Power(1.+(TMath::Power((x/[1]),[3])),[2]))",0.15,30);
  fFuncWeightPythia->SetParameters(0.36609,1.94635,1.40463,2.5);

  // weight function from the ratio of the LHC13d3 MC
  // and FONLL calculations for pp data
  fFuncWeightFONLL5overLHC13d3 = new TF1("funcWeightFONLL5overLHC13d3","([0]*x)/TMath::Power([2],(1+TMath::Power([3],x/[1])))+[4]*TMath::Exp([5]+[6]*x)+[7]*TMath::Exp([8]*x)",0.15,30.);
  fFuncWeightFONLL5overLHC13d3->SetParameters(2.94999e+00,3.47032e+00,2.81278e+00,2.5,1.93370e-02,3.86865e+00,-1.54113e-01,8.86944e-02,2.56267e-02);

  fFuncWeightFONLL5overLHC13d3Lc = new TF1("funcWeightFONLL5overLHC13d3Lc","([0]*x)/TMath::Power([2],(1+TMath::Power([3],x/[1])))+[4]*TMath::Exp([5]+[6]*x)+[7]*TMath::Exp([8]*x)",0.15,20.);
  fFuncWeightFONLL5overLHC13d3Lc->SetParameters(5.94428e+01,1.63585e+01,9.65555e+00,6.71944e+00,8.88338e-02,2.40477e+00,-4.88649e-02,-6.78599e-01,-2.10951e-01);

  PostData(6, fOutputKF);

  return;
}

//_________________________________________________
void AliAnalysisTaskSELc2V0bachelorTMVA::UserExec(Option_t *)
{
  /// user exec
  if (!fInputEvent) {
    AliError("NO EVENT FOUND!");
    return;
  }

  fCurrentEvent++;
  AliDebug(2, Form("Processing event = %d", fCurrentEvent));
  AliAODEvent* aodEvent = dynamic_cast<AliAODEvent*>(fInputEvent);
  TClonesArray *arrayLctopKos=0;

  TClonesArray *array3Prong = 0;

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
      arrayLctopKos=(TClonesArray*)aodFromExt->GetList()->FindObject("CascadesHF");

      array3Prong=(TClonesArray*)aodFromExt->GetList()->FindObject("Charm3Prong");
    }
  } else {
    arrayLctopKos=(TClonesArray*)aodEvent->GetList()->FindObject("CascadesHF");

    array3Prong=(TClonesArray*)aodEvent->GetList()->FindObject("Charm3Prong");
  }
  
  if (!fUseMCInfo && fIspA) {
    fAnalCuts->SetTriggerClass("");
    fAnalCuts->SetTriggerMask(fTriggerMask);
  }

  Int_t runnumber = aodEvent->GetRunNumber();
  if (aodEvent->GetTriggerMask() == 0 && (runnumber >= 195344 && runnumber <= 195677)){
    AliDebug(3,"Event rejected because of null trigger mask");
    return;
  }

  fCounter->StoreEvent(aodEvent,fAnalCuts,fUseMCInfo);

  // mc analysis
  TClonesArray *mcArray = 0;
  AliAODMCHeader *mcHeader=0;

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
      AliError("AliAnalysisTaskSELc2V0bachelorTMVA::UserExec: MC header branch not found!\n");
      return;
    }

    Double_t zMCVertex = mcHeader->GetVtxZ();
    if (TMath::Abs(zMCVertex) > fAnalCuts->GetMaxVtxZ()){
      AliDebug(3,Form("z coordinate of MC vertex = %f, it was required to be within [-%f, +%f], skipping event", zMCVertex, fAnalCuts->GetMaxVtxZ(), fAnalCuts->GetMaxVtxZ()));
      AliInfo(Form("z coordinate of MC vertex = %f, it was required to be within [-%f, +%f], skipping event", zMCVertex, fAnalCuts->GetMaxVtxZ(), fAnalCuts->GetMaxVtxZ()));
      return;
    }
    
    // multiplicity definition with tracklets
    fNTracklets = static_cast<Int_t>(AliVertexingHFUtils::GetNumberOfTrackletsInEtaRange(aodEvent,-1.,1.));

    // multiplicity definition with tracklets
    fNTracklets = static_cast<Int_t>(AliVertexingHFUtils::GetNumberOfTrackletsInEtaRange(aodEvent,-1.,1.));

    //Printf("Filling MC histo");
    FillMCHisto(mcArray);
  }

  // AOD primary vertex
  fVtx1 = (AliAODVertex*)aodEvent->GetPrimaryVertex();
  if (!fVtx1) return;
  if (fVtx1->GetNContributors()<1) return;

  fIsEventSelected = fAnalCuts->IsEventSelected(aodEvent);

  if ( !fIsEventSelected ) {
    fHistoEvents->Fill(0);
    return; // don't take into account not selected events
  }
  fHistoEvents->Fill(1);

  // Setting magnetic field for KF vertexing
  fBField = aodEvent->GetMagneticField();
  AliKFParticle::SetField(fBField);

  Int_t nSelectedAnal = 0;
  if (fIsK0sAnalysis) {
    MakeAnalysisForLc2prK0S(arrayLctopKos, mcArray,
			    nSelectedAnal, fAnalCuts,
			    array3Prong, mcHeader);
  }
  fCounter->StoreCandidates(aodEvent,nSelectedAnal,kFALSE);

  PostData(1, fOutput);
  PostData(2, fCounter);
  PostData(4, fVariablesTreeSgn);
  PostData(5, fVariablesTreeBkg);
  PostData(6, fOutputKF);

}
//-------------------------------------------------------------------------------
void AliAnalysisTaskSELc2V0bachelorTMVA::FillMCHisto(TClonesArray *mcArray){

  /// method to fill MC histo: how many Lc --> K0S + p are there at MC level
  for (Int_t iPart=0; iPart<mcArray->GetEntriesFast(); iPart++) {
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
    Int_t labeldaugh0 = mcPart->GetDaughter(0);
    Int_t labeldaugh1 = mcPart->GetDaughter(1);
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
	  Int_t labelK0daugh = v0MC->GetDaughter(0);
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
	      Int_t labelK0Sdaugh0 = partK0S->GetDaughter(0);
	      Int_t labelK0Sdaugh1 = partK0S->GetDaughter(1);
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
void AliAnalysisTaskSELc2V0bachelorTMVA::MakeAnalysisForLc2prK0S(TClonesArray *arrayLctopKos,
								 TClonesArray *mcArray,
								 Int_t &nSelectedAnal,
								 AliRDHFCutsLctoV0 *cutsAnal, TClonesArray *array3Prong,
								 AliAODMCHeader* aodheader){
  /// Lc prong needed to MatchToMC method

  Int_t pdgCand = 4122;
  Int_t pdgDgLctoV0bachelor[2]={2212, 310};
  Int_t pdgDgV0toDaughters[2]={211, 211};

  Int_t pdgDgLctopKpi[3]={2212, 321, 211};

  // loop to search for candidates Lc->p+K+pi
  Int_t n3Prong = array3Prong->GetEntriesFast();
  Int_t nCascades= arrayLctopKos->GetEntriesFast();

  //AliInfo(Form("\n\n\n\n3 prong candidates = %d, ncascades = %d \n\n\n\n\n", n3Prong, nCascades));
  for (Int_t i3Prong = 0; i3Prong < n3Prong; i3Prong++) {
    AliAODRecoDecayHF3Prong *d = (AliAODRecoDecayHF3Prong*)array3Prong->UncheckedAt(i3Prong);
    //Filling a control histogram with no cuts
    if (fUseMCInfo) {

      // find associated MC particle for Lc -> p+K+pi
      Int_t mcLabel = d->MatchToMC(4122, mcArray, 3, pdgDgLctopKpi);
      //Int_t mcLabelTemp = d->MatchToMC(4122, mcArray);
      //Printf("mcLabel = %d, mcLabelTemp = %d", mcLabel, mcLabelTemp);
      if (mcLabel >= 0) {

	AliAODMCParticle *partLcpKpi = dynamic_cast<AliAODMCParticle*>(mcArray->At(mcLabel));
	if(partLcpKpi){
	  Int_t pdgCode = partLcpKpi->GetPdgCode();
	  AliDebug(2,Form(" ¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯ MClabel=%d ~~~~~~~~~~ pdgCode=%d", mcLabel, pdgCode));
	  //AliInfo(Form(" ¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯ MClabel=%d ~~~~~~~~~~ pdgCode=%d", mcLabel, pdgCode));
	  fHistoLcpKpiBeforeCuts->Fill(1);

	}
      }
      else {
	//AliInfo(Form(" ¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯ MClabel=%d ~~~~~~~~~~", mcLabel));
	fHistoLcpKpiBeforeCuts->Fill(0);
      }
    }
  }

  // loop over cascades to search for candidates Lc->p+K0S

  Int_t mcLabel = -1;
  for (Int_t iLctopK0s = 0; iLctopK0s < nCascades; iLctopK0s++) {

    // Lc candidates and K0s from Lc
    AliAODRecoCascadeHF* lcK0spr = dynamic_cast<AliAODRecoCascadeHF*>(arrayLctopKos->At(iLctopK0s));
    if (!lcK0spr) {
      AliDebug(2,Form("Cascade %d doens't exist, skipping",iLctopK0s));
      continue;
    }

    if (!(lcK0spr->CheckCascadeFlags())) {
      AliDebug(2,Form("Cascade %d is not flagged as Lc candidate",iLctopK0s));
      continue;
    }

    //Filling a control histogram with no cuts
    if (fUseMCInfo) {

      Int_t pdgCode=-2;

      // find associated MC particle for Lc -> p+K0 and K0S->pi+pi
      fmcLabelLc = lcK0spr->MatchToMC(pdgCand, pdgDgLctoV0bachelor[1], pdgDgLctoV0bachelor, pdgDgV0toDaughters, mcArray, kTRUE);
      if (fmcLabelLc>=0) {
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
	pdgCode=-1;
      }
    }

    //if (!lcK0spr->GetSecondaryVtx()) {
    //	AliInfo("No secondary vertex");
    //continue;
    //}

    if (lcK0spr->GetNDaughters()!=2) {
      AliDebug(2,Form("Cascade %d has not 2 daughters (nDaughters=%d)",iLctopK0s,lcK0spr->GetNDaughters()));
      continue;
    }

    AliAODv0 * v0part = dynamic_cast<AliAODv0*>(lcK0spr->Getv0());
    AliAODTrack * bachPart = dynamic_cast<AliAODTrack*>(lcK0spr->GetBachelor());
    if (!v0part || !bachPart) {
      AliDebug(2,Form("Cascade %d has no V0 or no bachelor object",iLctopK0s));
      continue;
    }


    if (!v0part->GetSecondaryVtx()) {
      AliDebug(2,Form("No secondary vertex for V0 by cascade %d",iLctopK0s));
      continue;
    }

    if (v0part->GetNDaughters()!=2) {
      AliDebug(2,Form("current V0 has not 2 daughters (onTheFly=%d, nDaughters=%d)",v0part->GetOnFlyStatus(),v0part->GetNDaughters()));
      continue;
    }

    AliAODTrack * v0Pos = dynamic_cast<AliAODTrack*>(lcK0spr->Getv0PositiveTrack());
    AliAODTrack * v0Neg = dynamic_cast<AliAODTrack*>(lcK0spr->Getv0NegativeTrack());
    if (!v0Neg || !v0Pos) {
      AliDebug(2,Form("V0 by cascade %d has no V0positive of V0negative object",iLctopK0s));
      continue;
    }


    if (v0Pos->Charge() == v0Neg->Charge()) {
      AliDebug(2,Form("V0 by cascade %d has daughters with the same sign: IMPOSSIBLE!",iLctopK0s));
      continue;
    }

    Int_t isLc = 0;
    
    if (fUseMCInfo) {

      Int_t pdgCode = -2;

      // find associated MC particle for Lc -> p+K0 and K0S->pi+pi
      mcLabel = lcK0spr->MatchToMC(pdgCand, pdgDgLctoV0bachelor[1], pdgDgLctoV0bachelor, pdgDgV0toDaughters, mcArray, kTRUE);
      if (mcLabel>=0) {
	AliDebug(2,Form(" ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~cascade number %d (total cascade number = %d)", iLctopK0s, nCascades));

	AliAODMCParticle *partLc = dynamic_cast<AliAODMCParticle*>(mcArray->At(mcLabel));
	if(partLc){
	  pdgCode = partLc->GetPdgCode();
	  if (pdgCode<0) AliDebug(2,Form(" ¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯ MClabel=%d ~~~~~~~~~~ pdgCode=%d", mcLabel, pdgCode));
	  pdgCode = TMath::Abs(pdgCode);
	  isLc = 1;
	  fHistoLc->Fill(1);
	}
      }
      else {
	fHistoLc->Fill(0);
	pdgCode=-1;
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
	  // we have decided to fill the background only when the candidate has the daugthers that all come from HIJING underlying event!
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
    FillLc2pK0Sspectrum(lcK0spr, isLc, nSelectedAnal, cutsAnal, mcArray, mcLabel);
  }

  return;

}
//________________________________________________________________________
void AliAnalysisTaskSELc2V0bachelorTMVA::FillLc2pK0Sspectrum(AliAODRecoCascadeHF *part,
							     Int_t isLc,
							     Int_t &nSelectedAnal,
							     AliRDHFCutsLctoV0 *cutsAnal,
							     TClonesArray *mcArray, Int_t iLctopK0s){
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

  if ( (cutsAnal->IsInFiducialAcceptance(part->Pt(),part->Y(4122))) ) {
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

  Int_t isInV0window = (((cutsAnal->IsSelectedSingleCut(part, AliRDHFCuts::kCandidate, 2)) & (AliRDHFCutsLctoV0::kLcToK0Spr)) == (AliRDHFCutsLctoV0::kLcToK0Spr)); // cut on V0 invMass

  if (isInV0window == 0) {
    AliDebug(2, "No: The candidate has NOT passed the V0 window cuts!");
    if (isLc) Printf("SIGNAL candidate rejected: V0 window cuts");
    return;
  }
  else AliDebug(2, "Yes: The candidate has passed the mass cuts!");

  Bool_t isInCascadeWindow = (((cutsAnal->IsSelectedSingleCut(part, AliRDHFCuts::kCandidate, 0)) & (AliRDHFCutsLctoV0::kLcToK0Spr)) == (AliRDHFCutsLctoV0::kLcToK0Spr)); // cut on Lc->p+K0S invMass

  if (!isInCascadeWindow) {
    AliDebug(2, "No: The candidate has NOT passed the cascade window cuts!");
    if (isLc) Printf("SIGNAL candidate rejected: cascade window cuts");
    return;
  }
  else AliDebug(2, "Yes: The candidate has passed the cascade window cuts!");

  Bool_t isCandidateSelectedCuts = (((cutsAnal->IsSelected(part, AliRDHFCuts::kCandidate)) & (AliRDHFCutsLctoV0::kLcToK0Spr)) == (AliRDHFCutsLctoV0::kLcToK0Spr)); // kinematic/topological cuts
  AliDebug(2, Form("recoAnalysisCuts = %d", cutsAnal->IsSelected(part, AliRDHFCuts::kCandidate) & (AliRDHFCutsLctoV0::kLcToK0Spr)));
  if (!isCandidateSelectedCuts){
    AliDebug(2, "No: Analysis cuts kCandidate level NOT passed");
    if (isLc) Printf("SIGNAL candidate rejected");
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

  //Bool_t isBachelorID = (((cutsAnal->IsSelected(part,AliRDHFCuts::kPID))&(AliRDHFCutsLctoV0::kLcToK0Spr))==(AliRDHFCutsLctoV0::kLcToK0Spr)); // ID x bachelor
  Double_t probTPCTOF[AliPID::kSPECIES]={-1.};

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
    AliDebug(2,Form(" detUsed (TPC case) = %d", detUsed));
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

  if ( (((cutsAnal->IsSelected(part,AliRDHFCuts::kAll))&(AliRDHFCutsLctoV0::kLcToK0Spr))==(AliRDHFCutsLctoV0::kLcToK0Spr)) ) {
    nSelectedAnal++;
  }

  if ( !(cutsAnal->IsInFiducialAcceptance(part->Pt(),part->Y(4122))) ) return;

  if ( !( ( (cutsAnal->IsSelected(part, AliRDHFCuts::kTracks)) & (AliRDHFCutsLctoV0::kLcToK0Spr)) == (AliRDHFCutsLctoV0::kLcToK0Spr) ) ) { // esd track cuts
    if (isLc) Printf("SIGNAL candidate rejected");
    AliDebug(2, "No: Analysis cuts kTracks level NOT passed");
    return;
  }
  else {
    AliDebug(2, "Yes: Analysis cuts kTracks level passed");
  }

  Int_t pdgCand = 4122;
  Int_t pdgDgLctoV0bachelor[2]={211, 3122}; // case of wrong decay! Lc --> L + pi
  Int_t pdgDgV0toDaughters[2]={2212, 211}; // case of wrong decay! Lc --> L + pi
  Int_t isLc2LBarpi=0, isLc2Lpi=0;
  Int_t currentLabel = part->GetLabel();
  Int_t mcLabel = 0;
  if (fUseMCInfo) {
    mcLabel = part->MatchToMC(pdgCand, pdgDgLctoV0bachelor[1], pdgDgLctoV0bachelor, pdgDgV0toDaughters, mcArray, kTRUE);
    if (mcLabel>=0) {
      if (bachelor->Charge()==-1) isLc2LBarpi=1;
      if (bachelor->Charge()==+1) isLc2Lpi=1;
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
    if (lambdaLabel>=0) {
      AliAODMCParticle *lambdaTrack = (AliAODMCParticle*)mcArray->At(lambdaLabel);
      if (lambdaTrack->GetPdgCode()==3122) isLambda = 1;
      else if (lambdaTrack->GetPdgCode()==-3122) isLambdaBar = 1;
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
      if (gammaTrack->GetPdgCode()==22) isGamma = 1;
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

  //Double_t nSigmaITSpr=-999.;
  Double_t nSigmaTPCpr=-999.;
  Double_t nSigmaTOFpr=-999.;

  //Double_t nSigmaITSpi=-999.;
  Double_t nSigmaTPCpi=-999.;
  Double_t nSigmaTOFpi=-999.;

  //Double_t nSigmaITSka=-999.;
  Double_t nSigmaTPCka=-999.;
  Double_t nSigmaTOFka=-999.;

  /*
    cutsAnal->GetPidHF()->GetnSigmaITS(bachelor,4,nSigmaITSpr);
    cutsAnal->GetPidHF()->GetnSigmaTPC(bachelor,4,nSigmaTPCpr);
    cutsAnal->GetPidHF()->GetnSigmaTOF(bachelor,4,nSigmaTOFpr);
    cutsAnal->GetPidHF()->GetnSigmaITS(bachelor,2,nSigmaITSpi);
    cutsAnal->GetPidHF()->GetnSigmaTPC(bachelor,2,nSigmaTPCpi);
    cutsAnal->GetPidHF()->GetnSigmaTOF(bachelor,2,nSigmaTOFpi);
    cutsAnal->GetPidHF()->GetnSigmaITS(bachelor,3,nSigmaITSka);
    cutsAnal->GetPidHF()->GetnSigmaTPC(bachelor,3,nSigmaTPCka);
    cutsAnal->GetPidHF()->GetnSigmaTOF(bachelor,3,nSigmaTOFka);
  */

  nSigmaTPCpi = fPIDResponse->NumberOfSigmasTPC(bachelor,(AliPID::kPion));
  nSigmaTPCka = fPIDResponse->NumberOfSigmasTPC(bachelor,(AliPID::kKaon));
  nSigmaTPCpr = fPIDResponse->NumberOfSigmasTPC(bachelor,(AliPID::kProton));
  nSigmaTOFpi = fPIDResponse->NumberOfSigmasTOF(bachelor,(AliPID::kPion));
  nSigmaTOFka = fPIDResponse->NumberOfSigmasTOF(bachelor,(AliPID::kKaon));
  nSigmaTOFpr = fPIDResponse->NumberOfSigmasTOF(bachelor,(AliPID::kProton));
  
  Double_t ptLcMC = -1;
  Double_t weightPythia = -1, weight5LHC13d3 = -1, weight5LHC13d3Lc = -1; 

  if (fUseMCInfo) {
    if (iLctopK0s >= 0) {
      AliAODMCParticle *partLcMC = (AliAODMCParticle*)mcArray->At(iLctopK0s);
      ptLcMC = partLcMC->Pt();
      //Printf("--------------------- Reco pt = %f, MC particle pt = %f", part->Pt(), ptLcMC);
      weightPythia = fFuncWeightPythia->Eval(ptLcMC);
      weight5LHC13d3 = fFuncWeightFONLL5overLHC13d3->Eval(ptLcMC);
      weight5LHC13d3Lc = fFuncWeightFONLL5overLHC13d3Lc->Eval(ptLcMC);
    }
  }

  Double_t weightNch = 1;
  if (fUseMCInfo) {
    //Int_t nChargedMCPhysicalPrimary=AliVertexingHFUtils::GetGeneratedPhysicalPrimariesInEtaRange(mcArray,-1.0,1.0);
    //  if(nChargedMCPhysicalPrimary > 0){
    if(fNTracklets > 0){
      if(!fHistoMCNch) AliInfo("Input histos to evaluate Nch weights missing"); 
      if(fHistoMCNch) weightNch *= fHistoMCNch->GetBinContent(fHistoMCNch->FindBin(fNTracklets));
    }
  }
  

  // Fill candidate variable Tree (track selection, V0 invMass selection)
  if (!onFlyV0 && isInV0window && isInCascadeWindow && part->CosV0PointingAngle()>0.99 && TMath::Abs(nSigmaTPCpr) <= 3 && v0part->Getd0Prong(0) < 20 && v0part->Getd0Prong(1) < 20) {

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
    fCandidateVariables[10] = nSigmaTPCpi;
    fCandidateVariables[11] = nSigmaTPCka;
    fCandidateVariables[12] = nSigmaTOFpr;
    fCandidateVariables[13] = nSigmaTOFpi;
    fCandidateVariables[14] = nSigmaTOFka;
    fCandidateVariables[15] = bachelor->Pt();
    AliAODTrack *v0pos = (AliAODTrack*)part->Getv0PositiveTrack();
    fCandidateVariables[16] = v0pos->Pt();
    AliAODTrack *v0neg = (AliAODTrack*)part->Getv0NegativeTrack();
    fCandidateVariables[17] = v0neg->Pt();
    fCandidateVariables[18] = v0part->Getd0Prong(0);
    fCandidateVariables[19] = v0part->Getd0Prong(1);
    fCandidateVariables[20] = v0part->Pt();
    fCandidateVariables[21] = v0part->InvMass2Prongs(0,1,11,11);
    fCandidateVariables[22] = part->Pt();
    fCandidateVariables[23] = probProton;
    fCandidateVariables[24] = part->Eta();
    fCandidateVariables[25] = v0pos->Eta();
    fCandidateVariables[26] = v0neg->Eta();
    fCandidateVariables[27] = probProtonTPC;
    fCandidateVariables[28] = probProtonTOF;
    fCandidateVariables[29] = bachelor->Eta();

    fCandidateVariables[30] = part->P();
    fCandidateVariables[31] = bachelor->P();
    fCandidateVariables[32] = v0part->P();
    fCandidateVariables[33] = v0pos->P();
    fCandidateVariables[34] = v0neg->P();

    fCandidateVariables[35] = part->Y(4122);
    fCandidateVariables[36] = bachelor->Y(2212);
    fCandidateVariables[37] = v0part->Y(310);
    fCandidateVariables[38] = v0pos->Y(211);
    fCandidateVariables[39] = v0neg->Y(211);

    fCandidateVariables[40] = v0part->Eta();

    fCandidateVariables[41] = part->DecayLength();
    fCandidateVariables[42] = part->DecayLengthV0();
    fCandidateVariables[43] = part->Ct(4122);
    fCandidateVariables[44] = v0part->Ct(310, v0part->GetSecondaryVtx());

    EBachelor bachCode = kBachInvalid;
    EK0S k0SCode = kK0SInvalid;
    if (fUseMCInfo) {
      bachCode = CheckBachelor(part, bachelor, mcArray);
      k0SCode = CheckK0S(part, v0part, mcArray);
    }

    fCandidateVariables[45] = bachCode;
    fCandidateVariables[46] = k0SCode;

    Double_t V0KF[3] = {-999999, -999999, -999999};
    Double_t errV0KF[3] = {-999999, -999999, -999999};
    Double_t LcKF[3] = {-999999, -999999, -999999};
    Double_t errLcKF[3] = {-999999, -999999, -999999};
    Double_t distances[3] = {-999999, -999999, -999999};
    Double_t armPolKF[2] = {-999999, -999999};

    if (fCallKFVertexing){
      Int_t kfResult = CallKFVertexing(part, v0part, bachelor, mcArray, &V0KF[0], &errV0KF[0], &LcKF[0], &errLcKF[0], &distances[0], &armPolKF[0]);
      AliDebug(2, Form("Result from KF = %d", kfResult));
    }

    /*
      for (Int_t i = 0; i< 3; i++){
      Printf("i = %d, V0KF = %f, errV0KF = %f, LcKF = %f, errLcKF = %f", V0KF[i], errV0KF[i], LcKF[i], errLcKF[i]);
      }
    */

    fCandidateVariables[47] = V0KF[0];
    fCandidateVariables[48] = V0KF[1];
    fCandidateVariables[49] = V0KF[2];

    fCandidateVariables[50] = errV0KF[0];
    fCandidateVariables[51] = errV0KF[1];
    fCandidateVariables[52] = errV0KF[2];

    fCandidateVariables[53] = LcKF[0];
    fCandidateVariables[54] = LcKF[1];
    fCandidateVariables[55] = LcKF[2];

    fCandidateVariables[56] = errLcKF[0];
    fCandidateVariables[57] = errLcKF[1];
    fCandidateVariables[58] = errLcKF[2];

    fCandidateVariables[59] = distances[0];
    fCandidateVariables[60] = distances[1];
    fCandidateVariables[61] = distances[2];
    fCandidateVariables[62] = armPolKF[0];
    fCandidateVariables[63] = armPolKF[1];
    fCandidateVariables[64] = v0part->AlphaV0();
    fCandidateVariables[65] = v0part->PtArmV0();

    AliDebug(2, Form("v0pos->GetStatus() & AliESDtrack::kITSrefit= %d, v0neg->GetStatus() & AliESDtrack::kITSrefit = %d, v0pos->GetTPCClusterInfo(2, 1)= %f, v0neg->GetTPCClusterInfo(2, 1) = %f", (Int_t)(v0pos->GetStatus() & AliESDtrack::kITSrefit), (Int_t)(v0pos->GetStatus() & AliESDtrack::kITSrefit), v0pos->GetTPCClusterInfo(2, 1), v0neg->GetTPCClusterInfo(2, 1)));
    fCandidateVariables[66] = v0pos->GetStatus() & AliESDtrack::kITSrefit;
    fCandidateVariables[67] = v0neg->GetStatus() & AliESDtrack::kITSrefit;
    fCandidateVariables[68] = v0pos->GetTPCClusterInfo(2, 1);
    fCandidateVariables[69] = v0neg->GetTPCClusterInfo(2, 1);

    fCandidateVariables[70] = v0part->Xv();
    fCandidateVariables[71] = v0part->Yv();
    fCandidateVariables[72] = v0part->Zv();

    fCandidateVariables[73] = fVtx1->GetX();
    fCandidateVariables[74] = fVtx1->GetY();
    fCandidateVariables[75] = fVtx1->GetZ();

    fCandidateVariables[76] = bachelor->GetITSNcls();
    fCandidateVariables[77] = bachelor->HasPointOnITSLayer(0) + bachelor->HasPointOnITSLayer(1);

    fCandidateVariables[78] = v0pos->GetITSNcls();
    fCandidateVariables[79] = v0pos->HasPointOnITSLayer(0) + v0pos->HasPointOnITSLayer(1);

    fCandidateVariables[80] = v0neg->GetITSNcls();
    fCandidateVariables[81] = v0neg->HasPointOnITSLayer(0) + v0neg->HasPointOnITSLayer(1);

    TVector3 mom1(bachelor->Px(), bachelor->Py(), bachelor->Pz());
    TVector3 mom2(v0part->Px(), v0part->Py(), v0part->Pz());
    TVector3 momTot(part->Px(), part->Py(), part->Pz());

    Double_t Ql1 = mom1.Dot(momTot)/momTot.Mag();
    Double_t Ql2 = mom2.Dot(momTot)/momTot.Mag();

    Double_t alphaArmLc = (Ql1 - Ql2)/(Ql1 + Ql2);
    Double_t alphaArmLcCharge = ( bachelor->Charge() > 0 ? (Ql1 - Ql2)/(Ql1 + Ql2) : (Ql2 - Ql1)/(Ql1 + Ql2) );
    Double_t ptArmLc = mom1.Perp(momTot);

    fCandidateVariables[82] = alphaArmLc;
    fCandidateVariables[83] = alphaArmLcCharge;
    fCandidateVariables[84] = ptArmLc;

    Double_t massK0SPDG = TDatabasePDG::Instance()->GetParticle(310)->Mass();    // mass K0S PDG
    Double_t massPrPDG = TDatabasePDG::Instance()->GetParticle(2212)->Mass();    // mass Proton PDG
    Double_t massLcPDG = TDatabasePDG::Instance()->GetParticle(4122)->Mass();    // mass Lc PDG

    Double_t pStar = TMath::Sqrt((massLcPDG*massLcPDG-massPrPDG*massPrPDG-massK0SPDG*massK0SPDG)*(massLcPDG*massLcPDG-massPrPDG*massPrPDG-massK0SPDG*massK0SPDG)-4.*massPrPDG*massPrPDG*massK0SPDG*massK0SPDG)/(2.*massLcPDG);
    Double_t e = part->E(4122);
    Double_t beta = part->P()/e;
    Double_t gamma = e/massLcPDG;

    Double_t cts = (Ql1/gamma-beta*TMath::Sqrt(pStar*pStar+massPrPDG*massPrPDG))/pStar;

    fCandidateVariables[85] = cts;

    fCandidateVariables[86] = weightPythia;
    fCandidateVariables[87] = weight5LHC13d3;
    fCandidateVariables[88] = weight5LHC13d3Lc;
    fCandidateVariables[89] = weightNch;


    if (fUseMCInfo) {
      if (isLc){
	AliDebug(2, Form("Reco particle %d --> Filling Sgn", iLctopK0s));
	fVariablesTreeSgn->Fill();
	fHistoCodesSgn->Fill(bachCode, k0SCode);
      }
      else {
	if (fFillOnlySgn == kFALSE){
	  AliDebug(2, "Filling Bkg");
	  fVariablesTreeBkg->Fill();
	  fHistoCodesBkg->Fill(bachCode, k0SCode);
	}
      }
    }
    else {
      fVariablesTreeSgn->Fill();
    }
  }

  return;

}

//________________________________________________________________________
Int_t AliAnalysisTaskSELc2V0bachelorTMVA::CallKFVertexing(AliAODRecoCascadeHF *cascade, AliAODv0* v0part, AliAODTrack* bach, TClonesArray *mcArray,
							  Double_t* V0KF, Double_t* errV0KF, Double_t* LcKF, Double_t* errLcKF,
							  Double_t* distances, Double_t* armPolKF) {

  //
  /// method to perform KF vertexing
  /// elements: [0] = mass, [1] = DecayLength, [2] = lifeTime
  //

  Int_t codeKFV0 = -1, codeKFLc = -1;

  AliKFVertex primVtxCopy;
  Int_t nt = 0, ntcheck = 0;
  Double_t pos[3] = {0., 0., 0.};

  fVtx1->GetXYZ(pos);
  Int_t contr = fVtx1->GetNContributors();
  Double_t covmatrix[6] = {0.};
  fVtx1->GetCovarianceMatrix(covmatrix);
  Double_t chi2 = fVtx1->GetChi2();
  AliESDVertex primaryESDVtxCopy(pos, covmatrix, chi2, contr, "Vertex");

  // topological constraint
  primVtxCopy = AliKFVertex(primaryESDVtxCopy);
  nt = primaryESDVtxCopy.GetNContributors();
  ntcheck = nt;

  Int_t pdg[2] = {211, -211};
  Int_t pdgLc[2] = {2212, 310};

  Int_t pdgDgV0toDaughters[2] = {211, 211};

  Int_t mcLabelV0 = v0part->MatchToMC(310, mcArray, 2, pdgDgV0toDaughters);

  // the KF vertex for the V0 has to be built with the prongs of the V0!
  Bool_t isMCokV0 = kTRUE, isBkgV0 = kFALSE;
  AliKFParticle  V0, positiveV0KF, negativeV0KF;
  Int_t labelsv0daugh[2] = {-1, -1};
  Int_t idv0daugh[2] = {-1, -1};
  AliExternalTrackParam* esdv0Daugh1 = 0x0;
  AliExternalTrackParam* esdv0Daugh2 = 0x0;
  for(Int_t ipr= 0; ipr < 2; ipr++){ // 0 is positive, 1 is negative
    AliAODTrack *aodTrack = (AliAODTrack*)v0part->GetDaughter(ipr);
    if(!aodTrack) {
      AliDebug(2, "No V0 daughters available");
      return -1;
    }
    Double_t xyz[3], pxpypz[3], cv[21];
    Short_t sign;
    aodTrack->GetXYZ(xyz);
    aodTrack->PxPyPz(pxpypz);
    aodTrack->GetCovarianceXYZPxPyPz(cv);
    sign = aodTrack->Charge();
    AliExternalTrackParam tmp1( xyz, pxpypz, cv, sign);

    if (ipr == 0) esdv0Daugh1 = new AliExternalTrackParam( xyz, pxpypz, cv, sign);
    else esdv0Daugh2 = new AliExternalTrackParam( xyz, pxpypz, cv, sign);
    labelsv0daugh[ipr] = TMath::Abs(aodTrack->GetLabel());
    idv0daugh[ipr] = aodTrack->GetID();
    if (labelsv0daugh[ipr] == -1) isBkgV0 = kTRUE;

    //Printf("v0 daughter %d has label %d", ipr, labelsv0daugh[ipr]);

    AliKFParticle daughterKF(*aodTrack, pdg[ipr]); // we assume that the PDG is correct
    if (aodTrack->Charge() > 0) { // assigning positive and negative track to KF V0 for Armenteros-Podolanski plot
      positiveV0KF = daughterKF;
    }
    else {
      negativeV0KF = daughterKF;
    }
  }

  Double_t xn=0., xp=0.;//, dca;
  AliDebug(2, Form("bField = %f, esdv0Daugh1 = %p, esdv0Daugh2 = %p", fBField, esdv0Daugh1, esdv0Daugh2));
  //  dca = esdv0Daugh1->GetDCA(esdv0Daugh2, fBField, xn, xp);

  AliExternalTrackParam tr1(*esdv0Daugh1);
  AliExternalTrackParam tr2(*esdv0Daugh2);
  tr1.PropagateTo(xn, fBField);
  tr2.PropagateTo(xp, fBField);

  AliKFParticle daughterKF1(tr1, 211);
  AliKFParticle daughterKF2(tr2, 211);
  V0.AddDaughter(positiveV0KF);
  V0.AddDaughter(negativeV0KF);
  //V0.AddDaughter(daughterKF1);
  //V0.AddDaughter(daughterKF2);

  delete esdv0Daugh1;
  delete esdv0Daugh2;
  esdv0Daugh1=0;
  esdv0Daugh2=0;
  // Checking the quality of the KF V0 vertex
  if( V0.GetNDF() < 1 ) {
    //Printf("Number of degrees of freedom < 1, continuing");
    return -1;
  }
  if( TMath::Sqrt(TMath::Abs(V0.GetChi2()/V0.GetNDF())) > fCutKFChi2NDF ) {
    //Printf("Chi2 per DOF too big, continuing");
    return -1;
  }

  if(ftopoConstraint && nt > 0){
    for(Int_t ipr = 0; ipr < 2; ipr++){ // 0 is positive, 1 is negative
      AliAODTrack *aodTrack = (AliAODTrack*)v0part->GetDaughter(ipr);
      //* subtruct daughters from primary vertex
      if(!aodTrack->GetUsedForPrimVtxFit()) {
	//Printf("Track %d was not used for primary vertex, continuing", i);
	continue;
      }
      AliKFParticle daughterKF(*aodTrack, pdg[ipr]); // we assume that the PDG is correct
      primVtxCopy -= daughterKF;
      ntcheck--;
    }
  }

  // Check V0 Chi^2 deviation from primary vertex // not needed for V0 for Lc decay!!
  /*
  if( V0.GetDeviationFromVertex( primVtxCopy ) < fCutKFDeviationFromVtxV0) {
  //Printf("Deviation from vertex too big, continuing");
    return -1;
  }
  */

  //* Get V0 invariant mass
  Double_t massV0 = 999999, sigmaMassV0 = 999999;
  Int_t retMV0 = V0.GetMass( massV0, sigmaMassV0 );
  if( retMV0 ) {
    if (massV0 < 0) {
      codeKFV0 = 1; // Mass not ok
      if (sigmaMassV0 > 1e19) codeKFV0 = 5;  // Mass and SigmaMass not ok
    }
    else if (sigmaMassV0 > 1e19) codeKFV0 = 2; // SigmaMass not ok
  }
  fHistoMassKFV0->Fill(massV0, sigmaMassV0);

  if (massV0 < 0.4) Printf("\n\n>>>>>>>>>> Found the Funny V0 (mass = %f, sigma = %f, AOD mass = %f): labels of the tracks = %d, %d, id = %d and %d", massV0, sigmaMassV0, v0part->MassK0Short(), labelsv0daugh[0], labelsv0daugh[1], idv0daugh[0], idv0daugh[1]);
  if (massV0 > 0.55) Printf("\n\n>>>>>>>>>> Found the Funny V0 (mass = %f, , sigma = %f, AOD mass = %f): labels of the tracks = %d, %d, id = %d and %d", massV0, sigmaMassV0, v0part->MassK0Short(), labelsv0daugh[0], labelsv0daugh[1], idv0daugh[0], idv0daugh[1]);

  Printf("Vertices: KF:  x = %f, y = %f, z = %f", V0.GetX(), V0.GetY(), V0.GetZ());
  Printf("Vertices: AOD: x = %f, y = %f, z = %f", v0part->Xv(), v0part->Yv(), v0part->Zv());

  //Printf("Got MC vtx for V0");
  if (fUseMCInfo && TMath::Abs(labelsv0daugh[0] - labelsv0daugh[1]) == 1) {
    AliAODMCParticle* tmpdaughv01 = dynamic_cast<AliAODMCParticle*>(mcArray->At(labelsv0daugh[0]));
    AliAODMCParticle* tmpdaughv02 = dynamic_cast<AliAODMCParticle*>(mcArray->At(labelsv0daugh[1]));
    if (!tmpdaughv01 && labelsv0daugh[0] > 0){
      AliDebug(2, "Could not access MC info for first daughter of V0, continuing");
    }
    if (!tmpdaughv02 && labelsv0daugh[1] > 0){
      AliDebug(2, "Could not access MC info for second daughter of V0, continuing");
    }
    if(tmpdaughv01){
      Double_t xPionMC = tmpdaughv01->Xv(); //Production vertex of Pion --> Where K0S decays
      Double_t yPionMC = tmpdaughv01->Yv();
      Double_t zPionMC = tmpdaughv01->Zv();
      //Printf("Got MC vtx for Pion");
      Printf("Vertices: MC:  x = %f, y = %f, z = %f", xPionMC, yPionMC, zPionMC);
    }
  }
  else {
    Printf("Not a true V0");
  }
  //massV0=-1;//return -1;// !!!!

  // now use what we just try with the bachelor, to build the Lc

  // topological constraint
  nt = primVtxCopy.GetNContributors();
  ntcheck = nt;

  Bool_t isMCokLc = kTRUE, isBkgLc = kFALSE;
  AliKFParticle  Lc;
  Int_t labelsLcdaugh[2] = {-1, -1};
  labelsLcdaugh[0] = TMath::Abs(bach->GetLabel());
  labelsLcdaugh[1] = mcLabelV0;

  if (bach->Charge() < 0) pdgLc[0] = -pdgLc[0];
  AliKFParticle daughterKFLc(*bach, pdgLc[0]);
  Lc.AddDaughter(daughterKFLc);
  TParticlePDG* particlePDG = TDatabasePDG::Instance()->GetParticle(310);
  Double_t massPDGK0S = particlePDG->Mass();
  V0.SetMassConstraint(massPDGK0S);
  Lc.AddDaughter(V0);
  if( Lc.GetNDF() < 1 ) {
    AliDebug(2, Form("Lc: Number of degrees of freedom < 1 (%d), continuing", Lc.GetNDF()));
    return -1;
  }
  if( TMath::Sqrt(TMath::Abs(Lc.GetChi2()/Lc.GetNDF())) > fCutKFChi2NDF) {
    AliDebug(2, Form("Lc: Chi2 per DOF too big, continuing (%f)", TMath::Sqrt(TMath::Abs(Lc.GetChi2()/Lc.GetNDF()))));
    return -1;
  }

  if(ftopoConstraint && nt > 0){
    //* subtruct daughters from primary vertex
    if(!bach->GetUsedForPrimVtxFit()) {
      AliDebug(3, "Lc: Bachelor was not used for primary vertex, not subtracting it from primary vertex");
    }
    else{
      primVtxCopy -= daughterKFLc;
      ntcheck--;
    }
    /* the V0 was added above, so it is ok to remove it without checking
       if(!V0->GetUsedForPrimVtxFit()) {
       Printf("Lc: V0 was not used for primary vertex, continuing");
       continue;
       }
    */
    //primVtxCopy -= V0;
    //ntcheck--;
  }

  // Check Lc Chi^2 deviation from primary vertex
  /*
  if( Lc.GetDeviationFromVertex( primVtxCopy ) > fCutKFDeviationFromVtx) {
    AliDebug(2, Form("Lc: Deviation from vertex too big, continuing (%f)", Lc.GetDeviationFromVertex( primVtxCopy )));
    return -1;
  }

  if(ftopoConstraint){
    if(ntcheck>0) {
      // Add Lc to primary vertex to improve the primary vertex resolution
      primVtxCopy += Lc;
      Lc.SetProductionVertex(primVtxCopy);
    }
  }
  */
  //* Check chi^2
  if( TMath::Sqrt( TMath::Abs(Lc.GetChi2()/Lc.GetNDF())) > fCutKFChi2NDF) {
    AliDebug(2, Form("Lc: Final Chi2 per DOF too big, continuing (%f)", TMath::Sqrt( TMath::Abs(Lc.GetChi2()/Lc.GetNDF()))));
    return -1;
  }

  if(ftopoConstraint){
    V0.SetProductionVertex(Lc);
  }

  // After setting the vertex of the V0, getting/filling some info

  //* Get V0 decayLength
  Double_t decayLengthV0 = 999999, sigmaDecayLengthV0 = 999999;
  Int_t retDLV0 = V0.GetDecayLength( decayLengthV0, sigmaDecayLengthV0 );
  if( retDLV0 ) {
    if (sigmaDecayLengthV0 > 1e19) {
      codeKFV0 = 3; // DecayLength not ok
      if (massV0 < 0) {
	codeKFV0 = 6; // DecayLength and Mass not ok
	if (sigmaMassV0 > 1e19) codeKFV0 = 11;  // DecayLength and Mass and SigmaMass not ok
      }
      else if (sigmaMassV0 > 1e19) codeKFV0 = 8;  // DecayLength and SigmaMass not ok
    }
  }
  fHistoDecayLengthKFV0->Fill(decayLengthV0, sigmaDecayLengthV0);

  //* Get V0 life time
  Double_t lifeTimeV0 = 999999, sigmaLifeTimeV0 = 999999;
  Int_t retTLV0 = V0.GetLifeTime( lifeTimeV0, sigmaLifeTimeV0 );
  if( retTLV0 ) {
    if (sigmaLifeTimeV0 > 1e19) {
      codeKFV0 = 4;  // LifeTime not ok
      if (sigmaDecayLengthV0 > 1e19) {
	codeKFV0 = 9;  // LifeTime and DecayLength not ok
	if (massV0 < 0) {
	  codeKFV0 = 14; // LifeTime and DecayLength and Mass not ok
	  if (sigmaMassV0 > 1e19) codeKFV0 = 15; // LifeTime and DecayLength and Mass and SigmaMass not ok
	}
	else if (sigmaMassV0 > 1e19) codeKFV0 = 13; // LifeTime and DecayLength and SigmaMass not ok
      }
      else if (massV0 < 0) { // LifeTime and Mass and SigmaMass not ok
	codeKFV0 = 7; // LifeTime and Mass not ok
	if (sigmaMassV0 > 1e19) codeKFV0 = 12; // LifeTime and Mass and SigmaMass not ok
      }
      else if (sigmaMassV0 > 1e19) codeKFV0 = 10;  // LifeTime and SigmaMass not ok
    }
  }
  fHistoLifeTimeKFV0->Fill(lifeTimeV0, sigmaLifeTimeV0);

  if (codeKFV0 == -1) codeKFV0 = 0;
  fHistoKFV0->Fill(codeKFV0);

  AliDebug(2, Form("V0: mass = %f, decay length = %f, life time = %f", massV0, decayLengthV0, lifeTimeV0 ));

  fHistoMassV0All->Fill(massV0);
  fHistoDecayLengthV0All->Fill(decayLengthV0);
  fHistoLifeTimeV0All->Fill(lifeTimeV0);

  Double_t qtAlphaV0[2] = {0., 0.};
  Double_t vtxV0KF[3] = {V0.GetX(), V0.GetY(), V0.GetZ()};
  positiveV0KF.TransportToPoint(vtxV0KF);
  negativeV0KF.TransportToPoint(vtxV0KF);
  V0.GetArmenterosPodolanski(positiveV0KF, negativeV0KF, qtAlphaV0);
  AliDebug(2, Form("Armenteros-Podolanski variables: alpha = %f, qt = %f", qtAlphaV0[1], qtAlphaV0[0]));
  fHistoArmenterosPodolanskiV0KF->Fill(qtAlphaV0[1], qtAlphaV0[0]);
  fHistoArmenterosPodolanskiV0AOD->Fill(v0part->AlphaV0(), v0part->PtArmV0());
  armPolKF[0] = qtAlphaV0[1];
  armPolKF[1] = qtAlphaV0[0];

  // Checking MC info for V0

  AliAODMCParticle *motherV0 = 0x0;
  AliAODMCParticle *daughv01 = 0x0;
  AliAODMCParticle *daughv02 = 0x0;

  if (fUseMCInfo) {
    daughv01 = dynamic_cast<AliAODMCParticle*>(mcArray->At(labelsv0daugh[0]));
    daughv02 = dynamic_cast<AliAODMCParticle*>(mcArray->At(labelsv0daugh[1]));
    if (!daughv01 && labelsv0daugh[0] > 0){
      AliDebug(2, "Could not access MC info for first daughter of V0, continuing");
      isMCokV0 = kFALSE;
    }
    if (!daughv02 && labelsv0daugh[1] > 0){
      AliDebug(2, "Could not access MC info for second daughter of V0, continuing");
      isMCokV0 = kFALSE;
    }
    if (isMCokV0){
      if( daughv01->GetMother() ==  daughv02->GetMother() && daughv01->GetMother()>=0 ){
	AliDebug(3, Form("The mother has label %d", daughv01->GetMother()));
	motherV0 = dynamic_cast<AliAODMCParticle*>(mcArray->At(daughv01->GetMother()));
	if( motherV0 && TMath::Abs(motherV0->GetPdgCode()) != 21 ){ // These are all V0 that are truly V0, not only K0S, but no gluons
	  if( motherV0->GetNDaughters() == 2 ){
	    fHistoMassV0True->Fill(massV0);
	    fHistoDecayLengthV0True->Fill(decayLengthV0);
	    fHistoLifeTimeV0True->Fill(lifeTimeV0);
	    fHistoMassV0TrueFromAOD->Fill(v0part->MassK0Short());
	    if (TMath::Abs(motherV0->GetPdgCode()) == 310){ // These are true V0 that are also K0S
	      fHistoMassV0TrueK0S->Fill(massV0);
	      fHistoDecayLengthV0TrueK0S->Fill(decayLengthV0);
	      fHistoLifeTimeV0TrueK0S->Fill(lifeTimeV0);
	      fHistoMassV0TrueK0SFromAOD->Fill(v0part->MassK0Short());
	    }
	  }
	  AliDebug(2, Form("PDG V0 = %d, pi = %d, pj = %d, ndaughters = %d, mc mass = %f, reco mass = %f, v0 mass = %f", motherV0->GetPdgCode(), daughv01->GetPdgCode(), daughv02->GetPdgCode(), motherV0->GetNDaughters(), motherV0->GetCalcMass(), massV0, v0part->MassK0Short()));
	}
	else if (!motherV0){
	  AliDebug(3, "could not access MC info for mother, continuing");
	}
	else {
	  AliDebug(3, "MC mother is a gluon, continuing");
	}
      }
      else {
	AliDebug(3, "Background V0!");
	isBkgV0 = kTRUE;
      }
    }
  }

  AliDebug(2, Form("isMCokV0 = %d, isBkgV0 = %d", (Int_t)isMCokV0, (Int_t)isBkgV0));

  // Going back to Lc

  //* Get Lc invariant mass
  Double_t massLc = 999999, sigmaMassLc= 999999;
  Int_t retMLc = Lc.GetMass( massLc, sigmaMassLc );
  if( retMLc ) {
    AliDebug(3, Form("----> Could not get mass (%e), and sigma(%e) for Lc, continuing", massLc, sigmaMassLc));
    if (massLc < 0) {
      codeKFLc = 1; // Mass not ok
      if (sigmaMassLc > 1e19) codeKFLc = 5;  // Mass and SigmaMass not ok
    }
    else if (sigmaMassLc > 1e19) codeKFLc = 2; // SigmaMass not ok
  }
  fHistoMassKFLc->Fill(massLc, sigmaMassLc);

  //* Get Lc Decay length
  Double_t decayLengthLc = 999999, sigmaDecayLengthLc = 999999;
  Int_t retDLLc = Lc.GetDecayLength( decayLengthLc, sigmaDecayLengthLc );
  if( retDLLc ) {
    AliDebug(3, "----> Lc: Could not get decay length, and sigma");
    if (sigmaDecayLengthLc > 1e19) {
      codeKFLc = 3; // DecayLength not ok
      if (massLc < 0) {
	codeKFLc = 6; // DecayLength and Mass not ok
	if (sigmaMassLc > 1e19) codeKFLc = 11;  // DecayLength and Mass and SigmaMass not ok
      }
      else if (sigmaMassLc > 1e19) codeKFLc = 8;  // DecayLength and SigmaMass not ok
    }
  }
  AliDebug(3, Form("retDLLc = %d, with decayLength = %f and error = %e", retDLLc, decayLengthLc, sigmaDecayLengthLc));

  fHistoDecayLengthKFLc->Fill(decayLengthLc, sigmaDecayLengthLc);

  //* Get Lc life time
  Double_t lifeTimeLc = 999999, sigmaLifeTimeLc = 999999;
  Int_t retTLLc = Lc.GetLifeTime( lifeTimeLc, sigmaLifeTimeLc );
  if( retTLLc ) {
    AliDebug(3, "----> Lc: Could not get lifeTime, and sigma");
    if (sigmaLifeTimeLc > 1e19) {
      codeKFLc = 4;  // LifeTime not ok
      if (sigmaDecayLengthLc > 1e19) {
	codeKFLc = 9;  // LifeTime and DecayLength not ok
	if (massLc < 0) {
	  codeKFLc = 14; // LifeTime and DecayLength and Mass not ok
	  if (sigmaMassLc > 1e19) codeKFLc = 15; // LifeTime and DecayLength and Mass and SigmaMass not ok
	}
	else if (sigmaMassLc > 1e19) codeKFLc = 13; // LifeTime and DecayLength and SigmaMass not ok
      }
      else if (massLc < 0) { // LifeTime and Mass and SigmaMass not ok
	codeKFLc = 7; // LifeTime and Mass not ok
	if (sigmaMassLc > 1e19) codeKFLc = 12; // LifeTime and Mass and SigmaMass not ok
      }
      else if (sigmaMassLc > 1e19) codeKFLc = 10;  // LifeTime and SigmaMass not ok
    }
  }

  fHistoLifeTimeKFLc->Fill(lifeTimeLc, sigmaLifeTimeLc);

  AliDebug(2, Form("Lc: mass = %f (error = %e), decay length = %f (error = %e), life time = %f (error = %e) --> codeKFLc = %d", massLc, sigmaMassLc, decayLengthLc, sigmaDecayLengthLc, lifeTimeLc, sigmaLifeTimeLc, codeKFLc));

  if (codeKFLc == -1) codeKFLc = 0;
  fHistoKFLc->Fill(codeKFLc);

  fHistoKF->Fill(codeKFV0, codeKFLc);

  // here we fill the histgrams for all the reconstructed KF vertices for the cascade
  fHistoMassLcAll->Fill(massLc);
  fHistoDecayLengthLcAll->Fill(decayLengthLc);
  fHistoLifeTimeLcAll->Fill(lifeTimeLc);

  fHistoMassV0fromLcAll->Fill(massV0);
  fHistoDecayLengthV0fromLcAll->Fill(decayLengthV0);
  fHistoLifeTimeV0fromLcAll->Fill(lifeTimeV0);

  Double_t xV0 = V0.GetX();
  Double_t yV0 = V0.GetY();
  Double_t zV0 = V0.GetZ();

  Double_t xLc = Lc.GetX();
  Double_t yLc = Lc.GetY();
  Double_t zLc = Lc.GetZ();

  Double_t xPrimVtx = primVtxCopy.GetX();
  Double_t yPrimVtx = primVtxCopy.GetY();
  Double_t zPrimVtx = primVtxCopy.GetZ();

  Double_t distanceLcToPrimVtx = TMath::Sqrt((xPrimVtx - xLc) * (xPrimVtx - xLc) +
					     (yPrimVtx - yLc) * (yPrimVtx - yLc) +
					     (zPrimVtx - zLc) * (zPrimVtx - zLc));

  Double_t distanceV0ToPrimVtx = TMath::Sqrt((xPrimVtx - xV0) * (xPrimVtx - xV0) +
					     (yPrimVtx - yV0) * (yPrimVtx - yV0) +
					     (zPrimVtx - zV0) * (zPrimVtx - zV0));

  Double_t distanceV0ToLc = TMath::Sqrt((xLc - xV0)*(xLc - xV0) +
					(yLc - yV0)*(yLc - yV0) +
					(zLc - zV0)*(zLc - zV0));

  //Printf("distanceLcToPrimVtx = %e, distanceV0ToPrimVtx= %f, distanceV0ToLc = %f", distanceLcToPrimVtx, distanceV0ToPrimVtx, distanceV0ToLc);

  fHistoDistanceLcToPrimVtx->Fill(distanceLcToPrimVtx);
  fHistoDistanceV0ToPrimVtx->Fill(distanceV0ToPrimVtx);
  fHistoDistanceV0ToLc->Fill(distanceV0ToLc);

  distances[0] = distanceLcToPrimVtx;
  distances[1] = distanceV0ToPrimVtx;
  distances[2] = distanceV0ToLc;

  if (fUseMCInfo) {

    AliAODMCParticle *daughv01Lc = 0x0;
    AliAODMCParticle *K0S = 0x0;
    AliAODMCParticle *daughv02Lc = 0x0;

    if (labelsLcdaugh[0] >= 0) {
      //      Printf("Getting Bachelor from label %d", labelsLcdaugh[1]);
      daughv01Lc = dynamic_cast<AliAODMCParticle*>(mcArray->At(labelsLcdaugh[0]));
      if (!daughv01Lc){
	AliDebug(3, "Could not access MC info for first daughter of Lc");
	isMCokLc = kFALSE;
      }
      else {
	AliDebug(2, Form("The bachelor has label = %d", daughv01Lc->GetLabel()));
	if (TMath::Abs(daughv01Lc->GetPdgCode()) != 2212) isBkgLc = kTRUE;
      }
    }
    else { // The bachelor is a fake
      isBkgLc = kTRUE;
    }

    if (labelsLcdaugh[1] >= 0) {
      //Printf("Getting K0S");
      K0S = dynamic_cast<AliAODMCParticle*>(mcArray->At(labelsLcdaugh[1]));
      if (!K0S) {
	AliDebug(3, "Could not access MC info for second daughter of Lc");
	isMCokLc = kFALSE;
      }
      else{
	if (TMath::Abs(K0S->GetPdgCode()) != 310) isBkgLc = kTRUE;
      }
    }
    else{
      AliDebug(2, "The K0S is not true --> it does not have a label, continuing...");
      isBkgLc = kTRUE;
    }

    if (!isBkgLc){ // so far, we only checked that the V0 and the bachelor are not fake, and in particular, we know that the V0 is a K0S since we used the MatchToMC method
      if (isMCokLc) { // We can then access its MC info, and it might then be that also the Lc is a true Lc
	Int_t iK0 = K0S->GetMother();
	if (iK0 < 0) {
	  Printf("The K0S has no mother... IMPOSSIBLE"); // the K0S MUST have a mother!
	}
	else { // The K0S has a mother
	  daughv02Lc = dynamic_cast<AliAODMCParticle*>(mcArray->At(iK0));
	  if (!daughv02Lc){
	    AliDebug(3, "Could not access MC info for second daughter of Lc");
	  }
	  else { // we can access safely the K0S mother in the MC
	    if( daughv01Lc && (daughv01Lc->GetMother() ==  daughv02Lc->GetMother()) && (daughv01Lc->GetMother()>=0) ){  // This is a true cascade! bachelor and V0 come from the same mother
	      //Printf("Lc: The mother has label %d", daughv01Lc->GetMother());
	      AliAODMCParticle *motherLc = dynamic_cast<AliAODMCParticle*>(mcArray->At(daughv01Lc->GetMother()));
	      Int_t pdgMum = 0, pdgBach = 0, pdgV0 = 0;
	      if (motherLc) pdgMum = motherLc->GetPdgCode();
	      if (daughv01Lc) pdgBach = daughv01Lc->GetPdgCode();
	      if (daughv02Lc) pdgV0 = daughv02Lc->GetPdgCode();
	      AliDebug(2, Form("pdgMum = %d, pdgBach = %d, pdgV0 = %d", pdgMum, pdgBach, pdgV0));

	      if( motherLc && TMath::Abs(motherLc->GetPdgCode()) != 21 ){ // These are true cascades, we don't know yet if they are Lc
		fHistoMassLcTrue->Fill(massLc);
		fHistoDecayLengthLcTrue->Fill(decayLengthLc);
		fHistoLifeTimeLcTrue->Fill(lifeTimeLc);
		fHistoMassLcTrueFromAOD->Fill(cascade->InvMassLctoK0sP());

		fHistoMassV0fromLcTrue->Fill(massV0);
		fHistoDecayLengthV0fromLcTrue->Fill(decayLengthV0);
		fHistoLifeTimeV0fromLcTrue->Fill(lifeTimeV0);

		if (TMath::Abs(motherLc->GetPdgCode()) == 4122 && TMath::Abs(motherV0->GetPdgCode()) == 310 && TMath::Abs(daughv01Lc->GetPdgCode()) == 2212){  // This is Lc --> K0S + p (the check on the PDG code of the V0 is useless, since we used MathcToMC with it, but fine...
		  AliDebug(2, Form("IT IS SIGNAL!!! with label = %d", motherLc->GetLabel()));

		  fHistoArmenterosPodolanskiV0KFSgn->Fill(qtAlphaV0[1], qtAlphaV0[0]);
		  fHistoArmenterosPodolanskiV0AODSgn->Fill(v0part->AlphaV0(), v0part->PtArmV0());

		  fHistoDistanceLcToPrimVtxSgn->Fill(distanceLcToPrimVtx);
		  fHistoDistanceV0ToPrimVtxSgn->Fill(distanceV0ToPrimVtx);
		  fHistoDistanceV0ToLcSgn->Fill(distanceV0ToLc);

		  fHistoMassLcSgn->Fill(massLc);
		  fHistoMassLcSgnFromAOD->Fill(cascade->InvMassLctoK0sP());

		  fHistoDecayLengthLcSgn->Fill(decayLengthLc);
		  fHistoLifeTimeLcSgn->Fill(lifeTimeLc);

		  fHistoMassV0fromLcSgn->Fill(massV0);
		  fHistoDecayLengthV0fromLcSgn->Fill(decayLengthV0);
		  fHistoLifeTimeV0fromLcSgn->Fill(lifeTimeV0);
		}
		else {
		  isBkgLc = kTRUE;  // it is not a Lc, since the pdg != 4122
		}

		// if we got here, we can compare with MC information; this part is done also for cases where the cascade is not a Lc!
		// First, we compare the vtx of the Lc
		Double_t xLcMC = motherLc->Xv();
		Double_t yLcMC = motherLc->Yv();
		Double_t zLcMC = motherLc->Zv();
		//Printf("Got MC vtx for Lc");
		Double_t xProtonMC = daughv01Lc->Xv();
		Double_t yProtonMC = daughv01Lc->Yv();
		Double_t zProtonMC = daughv01Lc->Zv();
		//Printf("Got MC vtx for Proton");

		Double_t vtxLcResidualToPrimVtx = TMath::Sqrt((xLcMC - xProtonMC) * (xLcMC - xProtonMC) +
							      (yLcMC - yProtonMC) * (yLcMC - yProtonMC) +
							      (zLcMC - zProtonMC) * (zLcMC - zProtonMC)) - distanceLcToPrimVtx;

		// Then, we compare the vtx of the K0s
		Double_t xV0MC = motherV0->Xv();     //Production vertex of K0S --> Where Lc decays
		Double_t yV0MC = motherV0->Yv();
		Double_t zV0MC = motherV0->Zv();

		//Printf("Got MC vtx for V0");
		Double_t xPionMC = daughv01->Xv(); //Production vertex of Pion --> Where K0S decays
		Double_t yPionMC = daughv01->Yv();
		Double_t zPionMC = daughv01->Zv();
		//Printf("Got MC vtx for Pion");
		Printf("Vertices: MC:  x = %f, y = %f, z = %f", xPionMC, yPionMC, zPionMC);

		Double_t vtxV0ResidualToLc = TMath::Sqrt((xV0MC - xPionMC) * (xV0MC - xPionMC) +
							 (yV0MC - yPionMC) * (yV0MC - yPionMC) +
							 (zV0MC - zPionMC) * (zV0MC - zPionMC)) - distanceV0ToLc;

		// Then, the K0S vertex wrt primary vertex

		Double_t vtxV0ResidualToPrimVtx = TMath::Sqrt((xPionMC - xLcMC) * (xPionMC - xLcMC) +
							      (yPionMC - yLcMC) * (yPionMC - yLcMC) +
							      (zPionMC - zLcMC) * (zPionMC - zLcMC)) - distanceV0ToPrimVtx;

		fHistoVtxLcResidualToPrimVtx->Fill(vtxLcResidualToPrimVtx);
		fHistoVtxV0ResidualToLc->Fill(vtxV0ResidualToLc);
		fHistoVtxV0ResidualToPrimVtx->Fill(vtxV0ResidualToPrimVtx);

	      } //closing: if( motherLc && TMath::Abs(motherLc->GetPdgCode()) != 21 )
	      else if (!motherLc){
		AliDebug(2, "We could not access MC info for Lc mother, so we did nothing");
	      }
	      else {
		AliDebug(2, "MC Lc mother is a gluon, so we did nothing");
	      }
	    } //closing: if( daughv01Lc->GetMother() ==  daughv02Lc->GetMother() && daughv01Lc->GetMother()>=0 )
	    else {
	      isBkgLc = kTRUE; // it cannot be a Lc, since the daughters do not have the same mother
	    }
	  } // closing: else { // we can access safely the K0S mother in the MC
	} // closing: else { // The K0S has a mother
      } // closing isMCLcok
    }  // closing !isBkgLc
  } // closing fUseMCInfo

  //Printf("retMV0 = %d, retMLc = %d", retMV0, retMLc);
  if ( retMV0 == 0 && retMLc == 0){
    V0KF[0] = massV0;
    errV0KF[0] = sigmaMassV0;
    V0KF[1] = decayLengthV0;
    errV0KF[1] = sigmaDecayLengthV0;
    V0KF[2] = lifeTimeV0;
    errV0KF[2] = sigmaLifeTimeV0;
    LcKF[0] = massLc;
    errLcKF[0] = sigmaMassLc;
    LcKF[1] = decayLengthLc;
    errLcKF[1] = sigmaDecayLengthLc;
    LcKF[2] = lifeTimeLc;
    errLcKF[2] = sigmaLifeTimeLc;
  }

  return 1;

}
//________________________________________________________________________
AliAnalysisTaskSELc2V0bachelorTMVA::EBachelor AliAnalysisTaskSELc2V0bachelorTMVA::CheckBachelor( AliAODRecoCascadeHF *part,
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
AliAnalysisTaskSELc2V0bachelorTMVA::EK0S AliAnalysisTaskSELc2V0bachelorTMVA::CheckK0S( AliAODRecoCascadeHF *part,
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
Int_t AliAnalysisTaskSELc2V0bachelorTMVA::FindV0Label(AliAODRecoDecay* v0part, TClonesArray *mcArray) const
{

  //Printf(" FindV0Label");

  /// finding the label of teh V0; inspired from AliAODRecoDecay::MatchToMC

  Int_t labMother[2]={-1, -1};
  AliAODMCParticle *part=0;
  AliAODMCParticle *mother=0;
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

//----------------------------------------------------------------------------
Int_t AliAnalysisTaskSELc2V0bachelorTMVA::FindLcLabel(AliAODRecoCascadeHF* cascade, TClonesArray *mcArray) const
{

  /// finding the label of teh V0; inspired from AliAODRecoDecay::MatchToMC

  //Printf(" FindLcLabel");

  AliAODMCParticle *part=0;
  AliAODMCParticle *mother=0;
  AliAODMCParticle *grandmother=0;
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






