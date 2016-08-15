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

//////////////////////////////////////////////////////////////////////////
// Author: M. Nicassio m.nicassio@gsi.de maria.nicassio@cern.ch 
// First version 
// October 2014
//
// 
// Parts of this code are taken and adapted from the code below 
//
// PWGCF/FEMTOSCOPY/AliFemto and AliFemtoUser (AliFemto framework)
// PWGCF/FEMTOSCOPY/V0LamAnalysis (Jai)
// PWGCF/FEMTOSCOPY/Chaoticity (Dhevan)
// PWGCF/FEMTOSCOPY/K0S/PLamAnalysis (Hans) 
//  
// PWGLF/STRANGENESS 
//
// Possibility to run on
// - first particle: pion (0), proton (1)
// - second particle: pion (0) only when first particle pion (0), 
//   proton (1) only when first particle proton (1), Xi (2), Omega (3), 
// - AODs
// - ESDs (not mantained)
//
// November 2015
// - first attempt of shared daughter cut with selection criterion on inv mass (code optimization January 2016) 
// December 2015/January 2016
// - TH2D to run pipi (new binning for pions)
// - for identical particle case: _randomization_ code optimized 
// - small bug fixes (doSkipOver for protons not initialized, doPickOne for Xi not initialized, dphis calculation with analytical method)
// February 2016
// - dont use global tracks IP cut for same PID femto analysis, problems in CF (enhancement at low k* seen clearly in MC and data)
// - PropagateToDCA for global tracks to get IP changed according to AliAODTrack::PropagateToDCA code
// - new method for track propagation at fixed radius tested --> seems to work better for pi and p pairs than Hans propagation but to be further tested,
//   still no shift applied,
// - histos for momentum resolution correction, prim MC particle flag set for p and Xi and k*gen calculation
// May 2016
// - possibility to set Xi topological cuts (4 sets: default, tight, loose or otherwise reco)
// - added histos for p purity studies (data and MC) 
/////////////////////////////////////////////////////////////////////////

class AliESDVertex;

class AliAODVertex;

#include "TChain.h"

#include "TList.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TMath.h"
#include "TRandom.h"
#include "TDatabasePDG.h"

#include "AliAnalysisManager.h"
#include "AliAnalysisTask.h"
#include "AliInputEventHandler.h"

#include "AliCentrality.h"
#include "AliPIDResponse.h"

#include "AliVTrack.h"

#include "AliESDEvent.h"
#include "AliESDtrackCuts.h"
#include "AliESDcascade.h"
#include "AliESDv0.h"

#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODcascade.h"
#include "AliAODv0.h"
#include "AliAODMCParticle.h"


#include "AliStack.h"
#include "AliMCEvent.h"

#include "AliVVertex.h"
#include "AliAnalysisTaskhCascadeFemto.h"

#define PI 3.1415927
 
using namespace std;


ClassImp(AliAnalysisTaskhCascadeFemto)

AliAnalysisTaskhCascadeFemto::AliAnalysisTaskhCascadeFemto():AliAnalysisTaskSE(),
  fESDevent(NULL),
  fAODevent(NULL),
  fAnalysisType("ESD"),
  fFirstpart(kPion),
  fSecondpart(kXi),
  fXicutSet(kdefault),
  fMassWindowCascades(0.008),
  fnSigmaTPCPIDfirstParticle(3.),
  fnSigmaTPCTOFPIDfirstParticle(3.),
  fnSigmaTPCPIDsecondParticleDau(4.),
  fkCascadeSideBands(kFALSE),
  
  fReadMCTruth(kFALSE),
  fUseContainer(kFALSE),
  fUseStandardCuts(0),
  fkApplyTtc(kFALSE),
  fDphisMin(0.),
  fDetasMin(0.),
  fMomemtumLimitForTOFPID(0.8),
  fkApplyRatioCrRnFindCut(kFALSE),
  fkCutOnCrossedRows(kFALSE),
  fkCutOnTPCIP(kFALSE),  
  fIPCutxy(0.1),
  fIPCutz(0.15),
  fMinPtForPrim(0.7),
  fMaxPtForPrim(4.),
  fMinPtForCasc(0.8),
  fMaxPtForCasc(10.),
  
  // Cascade and V0 topological and kinematical selection
  fIPBac(0.03),
  fDcaXiDaughters(0.3),
  fXiCosineOfPointingAngle(0.999),
  fXiRadiusMin(0.9),
  fXiRadiusMax(100.),
  fCtau(999999.), 
  fIPV0ToPrimVertexXi(0.05),
  fMassWindowL(0.008),
  fV0CosineOfPointingAngle(0.998),
  fV0RadiusXiMin(0.9),
  fV0RadiusXiMax(100.),
  fDcaV0Daughters(1.),
  fIPPosToPrimVertexXi(0.1),
  fIPNegToPrimVertexXi(0.1),
  fkApplyYcutCasc(kFALSE),
  fkPropagateGlobal(kFALSE),
  fkPropagateAtFixedR(kFALSE),
  fkCutOnTtcProp(kFALSE),
  fkApplySharedDaughterCutXi(kTRUE),

  fESDtrackCuts(0),
  fPIDResponse(0),

  fCentrLowLim(0),
  fCentrUpLim(0),

  fNbinsMass(0),
  fLowMassLim(0),
  fUpMassLim(0),

  fNbinsMassGm(0),
  fLowMassLimGm(0),
  fUpMassLimGm(0),

  farrGT(0),fTrackBufferSize(20200), // was 18k
  fEventColl(0x0),
  fEvt(0x0),
  fMaxPMult(3000), // 1000 for protons 
  fMaxXiMult(20),  // was 100
  fPDGXi(1.32171),
  fPDGOmega(1.6724),
  fPDGp(0.938272046),
  fPDGL(1.115683),
  fPDGpi(0.13957),
  fPDGfirst(0.),
  fPDGsecond(0.),
  fzVertexBins(10),
  fnCentBins(20),
  fnEventsToMix(7),


  fHistCentrality(0x0),
  fHistVertexDistribution(0x0),
  fHistMultiplicityOfMixedEvent(0x0),

  fHistnTPCCrossedR(0x0),
  fHistRationTPCCrossedRnFind(0x0),
  fHistSharedFrTPCcl(0x0), 
  fHistprimpTPCdEdx(0x0),

  fHistyptProtons(0x0),
  fHistphietaProtons(0x0),
  fHistIPtoPVxyzTPC(0x0),
  fHistIPtoPVxyzGlobal(0x0),

  fHistpTOFmisvsp(0x0),
  fHistpTOFnsigmavsp(0x0),
  fHistpTOFsignalvsp(0x0),
  fHistpnsigTOFnsigTPC(0x0),
  fHistpsignalTOFsignalTPC(0x0),
  fHistProtonMultvsCent(0x0),

  fHistMCPrimProtons(0x0),
  fHistMCFromWdecayProtons(0x0),
  fHistMCFromMaterialProtons(0x0),
  fHistMCOtherProtons(0x0),
  fHistMCPrimAProtons(0x0),
  fHistMCFromWdecayAProtons(0x0),
  fHistMCFromMaterialAProtons(0x0),
  fHistMCOtherAProtons(0x0),

  fHistpTPCdEdx(0x0),
  fHistnTPCdEdx(0x0),
  fHistbTPCdEdx(0x0),

  fHistPosV0TPCClusters(0x0),
  fHistNegV0TPCClusters(0x0),
  fHistBachTPCClusters(0x0),

  fHistSharedFrTPCclPos(0x0),
  fHistSharedFrTPCclNeg(0x0),
  fHistSharedFrTPCclBach(0x0),

  fHistyptXi(0x0),
  fHistphietaXi(0x0),
  fHistptL(0x0),
  fHistptpL(0x0),
  fHistptnL(0x0),
  fHistptbac(0x0),

  fHistInvMassXiMinus(0x0),
  fHistInvMassL(0x0),
  fHistInvMassXiPlus(0x0),
  fHistInvMassAntiL(0x0),
  fHistInvMassXiInPairs(0x0),
  fHistXiMultvsCent(0x0),
  fHistIPtoPVxyGlobalvspt(0x0),
//  fCFContCascadeCuts(0x0),

  fHistpXiSignalRealKstar(0x0),
  fHistapXiSignalRealKstar(0x0),
  fHistpaXiSignalRealKstar(0x0),
  fHistapaXiSignalRealKstar(0x0),
  fHistpXiSignalMixedKstargenvsrec(0x0),
  fHistapXiSignalMixedKstargenvsrec(0x0),
  fHistpaXiSignalMixedKstargenvsrec(0x0),
  fHistapaXiSignalMixedKstargenvsrec(0x0),

  fHistpXiSignalBkgKstar(0x0),
  fHistapXiSignalBkgKstar(0x0),
  fHistpaXiSignalBkgKstar(0x0),
  fHistapaXiSignalBkgKstar(0x0),
  fHistFractionOfXiWithSharedDaughters(0x0),
  fHistTotMaxFractionOfXiWithSharedDaughters(0x0),

  fHistpXibacDetaSDphiS(0x0),
  fHistpXiposDetaSDphiS(0x0),
  fHistpXinegDetaSDphiS(0x0),
  fHistpaXibacDetaSDphiS(0x0),
  fHistpaXiposDetaSDphiS(0x0),
  fHistpaXinegDetaSDphiS(0x0),
  fHistapXibacDetaSDphiS(0x0),
  fHistapXiposDetaSDphiS(0x0),
  fHistapXinegDetaSDphiS(0x0),
  fHistapaXibacDetaSDphiS(0x0),
  fHistapaXiposDetaSDphiS(0x0),
  fHistapaXinegDetaSDphiS(0x0),

  fHistpXibacDetaSDphiSBkg(0x0),
  fHistpXiposDetaSDphiSBkg(0x0),
  fHistpXinegDetaSDphiSBkg(0x0),
  fHistpaXibacDetaSDphiSBkg(0x0),
  fHistpaXiposDetaSDphiSBkg(0x0),
  fHistpaXinegDetaSDphiSBkg(0x0),
  fHistapXibacDetaSDphiSBkg(0x0),
  fHistapXiposDetaSDphiSBkg(0x0),
  fHistapXinegDetaSDphiSBkg(0x0),
  fHistapaXibacDetaSDphiSBkg(0x0),
  fHistapaXiposDetaSDphiSBkg(0x0),
  fHistapaXinegDetaSDphiSBkg(0x0),

  fHistTrackBufferOverflow(0x0),

  fOutputContainer(NULL)

{
  fPDGp = TDatabasePDG::Instance()->GetParticle(2212)->Mass();
  fPDGpi = TDatabasePDG::Instance()->GetParticle(211)->Mass();

  fPDGXi = TDatabasePDG::Instance()->GetParticle(3312)->Mass();
  fPDGOmega = TDatabasePDG::Instance()->GetParticle(3334)->Mass();
  fPDGL = TDatabasePDG::Instance()->GetParticle(3122)->Mass();
  

  // default Constructor
  // Define input and output slots here
}

//_______________________________________________________________

AliAnalysisTaskhCascadeFemto::AliAnalysisTaskhCascadeFemto(const char *name):AliAnalysisTaskSE(name),
  fESDevent(NULL),
  fAODevent(NULL),
  fAnalysisType("ESD"),
  fFirstpart(kPion),
  fSecondpart(kXi),
  fXicutSet(kdefault),
  fMassWindowCascades(0.008),
  fnSigmaTPCPIDfirstParticle(3.),
  fnSigmaTPCTOFPIDfirstParticle(3.),
  fnSigmaTPCPIDsecondParticleDau(4.),
  fkCascadeSideBands(kFALSE),
 
  fReadMCTruth(kFALSE),
  fUseContainer(kFALSE),
  fUseStandardCuts(0),
  fkApplyTtc(kFALSE),
  fDphisMin(0.),
  fDetasMin(0.),
  fMomemtumLimitForTOFPID(0.8),
  fkApplyRatioCrRnFindCut(kFALSE),
  fkCutOnCrossedRows(kFALSE),
  fkCutOnTPCIP(kFALSE),
  fIPCutxy(0.1),
  fIPCutz(0.15),
  fMinPtForPrim(0.7),
  fMaxPtForPrim(4.),
  fMinPtForCasc(0.8),
  fMaxPtForCasc(10.),
    
  // Cascade and V0 topological and kinematical selection
  fIPBac(0.03),
  fDcaXiDaughters(0.3),
  fXiCosineOfPointingAngle(0.999),
  fXiRadiusMin(0.9),
  fXiRadiusMax(100.),
  fCtau(999999.), 
  fIPV0ToPrimVertexXi(0.05),
  fMassWindowL(0.008),
  fV0CosineOfPointingAngle(0.998),
  fV0RadiusXiMin(0.9),
  fV0RadiusXiMax(100.),
  fDcaV0Daughters(1.),
  fIPPosToPrimVertexXi(0.1),
  fIPNegToPrimVertexXi(0.1),

  fkApplyYcutCasc(kFALSE),
  fkPropagateGlobal(kFALSE),
  fkPropagateAtFixedR(kFALSE),
  fkCutOnTtcProp(kFALSE),
  fkApplySharedDaughterCutXi(kTRUE),

  fESDtrackCuts(0),
  fPIDResponse(0),

  fCentrLowLim(0),
  fCentrUpLim(0),                        
  
  fNbinsMass(0),                          
  fLowMassLim(0),
  fUpMassLim(0),
  
  fNbinsMassGm(0),                                             
  fLowMassLimGm(0),
  fUpMassLimGm(0),

  farrGT(0),fTrackBufferSize(20200), // was 18k
  fEventColl(0x0),
  fEvt(0x0),
  fMaxPMult(3000), // 1000 for protons 
  fMaxXiMult(20),  // was 100
  fPDGXi(1.32171),
  fPDGOmega(1.6724), 
  fPDGp(0.938272046),
  fPDGL(1.115683),
  fPDGpi(0.13957),
  fPDGfirst(0.),
  fPDGsecond(0.),
  fzVertexBins(10), 
  fnCentBins(20),
  fnEventsToMix(7),


  fHistCentrality(0x0),
  fHistVertexDistribution(0x0),
  fHistMultiplicityOfMixedEvent(0x0),

  fHistnTPCCrossedR(0x0),
  fHistRationTPCCrossedRnFind(0x0),
  fHistSharedFrTPCcl(0x0),
  fHistprimpTPCdEdx(0x0),
 
  fHistyptProtons(0x0),
  fHistphietaProtons(0x0),
  fHistIPtoPVxyzTPC(0x0),
  fHistIPtoPVxyzGlobal(0x0),

  fHistpTOFmisvsp(0x0),
  fHistpTOFnsigmavsp(0x0),
  fHistpTOFsignalvsp(0x0),
  fHistpnsigTOFnsigTPC(0x0),
  fHistpsignalTOFsignalTPC(0x0),
  fHistProtonMultvsCent(0x0),

  fHistMCPrimProtons(0x0),
  fHistMCFromWdecayProtons(0x0),
  fHistMCFromMaterialProtons(0x0),
  fHistMCOtherProtons(0x0),
  fHistMCPrimAProtons(0x0),
  fHistMCFromWdecayAProtons(0x0),
  fHistMCFromMaterialAProtons(0x0),
  fHistMCOtherAProtons(0x0),

  fHistpTPCdEdx(0x0),
  fHistnTPCdEdx(0x0),
  fHistbTPCdEdx(0x0),

  fHistPosV0TPCClusters(0x0),
  fHistNegV0TPCClusters(0x0),
  fHistBachTPCClusters(0x0),

  fHistSharedFrTPCclPos(0x0),
  fHistSharedFrTPCclNeg(0x0),
  fHistSharedFrTPCclBach(0x0),

  fHistyptXi(0x0),
  fHistphietaXi(0x0),
  fHistptL(0x0),
  fHistptpL(0x0),
  fHistptnL(0x0),
  fHistptbac(0x0),

  fHistInvMassXiMinus(0x0),
  fHistInvMassL(0x0),
  fHistInvMassXiPlus(0x0),
  fHistInvMassAntiL(0x0),
  fHistInvMassXiInPairs(0x0),
  fHistXiMultvsCent(0x0),
  fHistIPtoPVxyGlobalvspt(0x0),
//  fCFContCascadeCuts(0x0),

  fHistpXiSignalRealKstar(0x0),
  fHistapXiSignalRealKstar(0x0),
  fHistpaXiSignalRealKstar(0x0),
  fHistapaXiSignalRealKstar(0x0),
  fHistpXiSignalMixedKstargenvsrec(0x0),
  fHistapXiSignalMixedKstargenvsrec(0x0),
  fHistpaXiSignalMixedKstargenvsrec(0x0),
  fHistapaXiSignalMixedKstargenvsrec(0x0),

  fHistpXiSignalBkgKstar(0x0),
  fHistapXiSignalBkgKstar(0x0),
  fHistpaXiSignalBkgKstar(0x0),
  fHistapaXiSignalBkgKstar(0x0),
  fHistFractionOfXiWithSharedDaughters(0x0),
  fHistTotMaxFractionOfXiWithSharedDaughters(0x0),

  fHistpXibacDetaSDphiS(0x0),
  fHistpXiposDetaSDphiS(0x0),
  fHistpXinegDetaSDphiS(0x0),
  fHistpaXibacDetaSDphiS(0x0),
  fHistpaXiposDetaSDphiS(0x0),
  fHistpaXinegDetaSDphiS(0x0),
  fHistapXibacDetaSDphiS(0x0),
  fHistapXiposDetaSDphiS(0x0),
  fHistapXinegDetaSDphiS(0x0),
  fHistapaXibacDetaSDphiS(0x0),
  fHistapaXiposDetaSDphiS(0x0),
  fHistapaXinegDetaSDphiS(0x0),

  fHistpXibacDetaSDphiSBkg(0x0),
  fHistpXiposDetaSDphiSBkg(0x0),
  fHistpXinegDetaSDphiSBkg(0x0),
  fHistpaXibacDetaSDphiSBkg(0x0),
  fHistpaXiposDetaSDphiSBkg(0x0),
  fHistpaXinegDetaSDphiSBkg(0x0),
  fHistapXibacDetaSDphiSBkg(0x0),
  fHistapXiposDetaSDphiSBkg(0x0),
  fHistapXinegDetaSDphiSBkg(0x0),
  fHistapaXibacDetaSDphiSBkg(0x0),
  fHistapaXiposDetaSDphiSBkg(0x0),
  fHistapaXinegDetaSDphiSBkg(0x0),

  fHistTrackBufferOverflow(0x0),

  fOutputContainer(NULL)
{
  fPDGp = TDatabasePDG::Instance()->GetParticle(2212)->Mass(); 
  fPDGpi = TDatabasePDG::Instance()->GetParticle(211)->Mass();
  fPDGXi = TDatabasePDG::Instance()->GetParticle(3312)->Mass();
  fPDGOmega = TDatabasePDG::Instance()->GetParticle(3334)->Mass();
  fPDGL = TDatabasePDG::Instance()->GetParticle(3122)->Mass(); 
  // Constructor
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TH1 container
  DefineOutput(1, TList::Class());
//  DefineOutput(2, AliCFContainer::Class());

}

//_______________________________________________________________

AliAnalysisTaskhCascadeFemto::~AliAnalysisTaskhCascadeFemto() {
  //
  // Destructor
  //
  if (fOutputContainer && !AliAnalysisManager::GetAnalysisManager()->IsProofMode())   { delete fOutputContainer;     fOutputContainer = 0x0;    }
//  if (fCFContCascadeCuts && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) { delete fCFContCascadeCuts; fCFContCascadeCuts = 0x0; }
  if (fESDtrackCuts) delete fESDtrackCuts;


  if (farrGT)
    delete[] farrGT;
  farrGT=0;


  for (unsigned short i=0; i<fzVertexBins; i++) {
    for (unsigned short j=0; j<fnCentBins; j++) {
      delete fEventColl[i][j];
    }
    delete[] fEventColl[i];
  }
  delete[] fEventColl;

}

//_______________________________________________________________

void AliAnalysisTaskhCascadeFemto::UserCreateOutputObjects() {

  // PID object
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  fPIDResponse = inputHandler->GetPIDResponse();

  // Event mixing part

  fEventColl = new AliAnalysishCascadeEventCollection **[fzVertexBins]; 

  for (unsigned short i=0; i<fzVertexBins; i++) {
    fEventColl[i] = new AliAnalysishCascadeEventCollection *[fnCentBins];
    for (unsigned short j=0; j<fnCentBins; j++) {
      fEventColl[i][j] = new AliAnalysishCascadeEventCollection(fnEventsToMix+1, fMaxXiMult, fMaxPMult);
    }
  }

  // Set PDG for first and second particle
  if (fFirstpart == kPion) fPDGfirst = fPDGpi; 
  else if (fFirstpart == kProton) fPDGfirst = fPDGp;
  else { cout<<" First particle not known or not possible to deal with in this task for the moment! "<<endl; return;}
  if (fSecondpart == kPion) fPDGsecond = fPDGpi;
  else if (fSecondpart == kProton) fPDGsecond = fPDGp;
  else if (fSecondpart == kXi) fPDGsecond = fPDGXi;
  else if (fSecondpart == kOmega) fPDGsecond = fPDGOmega;
  else { cout<<" Second particle not known or not possible to deal with in this task for the moment! "<<endl; return;}

    // set of topological cuts for cascades (default or for systematic studies)
  fIPBac = 0.03; // >
  fDcaXiDaughters = 0.3; //<
  fXiRadiusMax = 1000.; //<
  fIPPosToPrimVertexXi = 0.1;  // >
  fIPNegToPrimVertexXi = 0.1;  // >
  fV0CosineOfPointingAngle = 0.998; // >
  fV0RadiusXiMax = 1000.; //<
  fCtau = 99999.; // >

  //cout<<"Set of cuts: "<<fXicutSet<<endl;

  if (fXicutSet == kdefault) {

    fXiCosineOfPointingAngle = 0.9992;// >
    if (fSecondpart == kXi) {
      fXiRadiusMin = 1.5; // >
      fCtau = 14.9; // <
    } else if (fSecondpart == kOmega) {
      fXiRadiusMin = 1.;
      fCtau = 7.9;
    }
    fIPV0ToPrimVertexXi = 0.1; // >
    fMassWindowL = 0.005; // <>
    fV0RadiusXiMin = 3.; // >
    fDcaV0Daughters = 0.8; // <
    // window xi mass for omega 0.008
    
  } else if (fXicutSet == kloose) { 

    fXiCosineOfPointingAngle = 0.9991;
    if (fSecondpart == kXi) {
      fXiRadiusMin = 1.3;
      fCtau = 17.5;
    } else if (fSecondpart == kOmega) {
      fXiRadiusMin = 0.9; 
      fCtau = 8.9;
    }
    fIPV0ToPrimVertexXi = 0.09;
    fMassWindowL = 0.006;
    fV0RadiusXiMin = 2.5;
    fDcaV0Daughters = 0.9;

  } else if (fXicutSet == ktight) {

    fXiCosineOfPointingAngle = 0.9993; 
    if (fSecondpart == kXi) {
      fXiRadiusMin = 1.7;
      fCtau = 12.5;
    } else if (fSecondpart == kOmega) {
      fXiRadiusMin = 1.2;  
      fCtau = 6.9;
    }
    fIPV0ToPrimVertexXi = 0.12;
    fMassWindowL = 0.004;
    fV0RadiusXiMin = 4.;
    fDcaV0Daughters = 0.7;

  } else { // kreco

    fXiCosineOfPointingAngle = 0.999;
    fXiRadiusMin = 0.9;
    fIPV0ToPrimVertexXi = 0.05;
    fMassWindowL = 0.008;
    // 7 in V0 vertexer (-chi2) + mass
    fV0RadiusXiMin = 0.9;
    fDcaV0Daughters = 1.;

  } 

  /*cout<<"Cut "<<fDcaV0Daughters<<endl;
  cout<<"Cut "<<fV0RadiusXiMin<<endl;
  cout<<"Cut "<<fMassWindowL<<endl;
  cout<<"Cut "<<fIPV0ToPrimVertexXi<<endl;
  cout<<"Cut "<<fXiRadiusMin<<endl;
  cout<<"Cut "<<fXiCosineOfPointingAngle<<endl;
  cout<<"Cut "<<fIPBac<<endl;
  cout<<"Cut "<<fCtau<<endl;*/


  // Store pointer to global tracks
  farrGT = new Int_t[fTrackBufferSize];
  
  //Define and fill the OutputContainer
  fOutputContainer = new TList();
  fOutputContainer->SetOwner(kTRUE);


  // Create histograms 
  fHistCentrality = new TH1F("fHistCentrality", "Number of events", 20, 0., 100.);
  fHistCentrality ->GetXaxis()->SetTitle("Centrality");
  fHistCentrality ->GetYaxis()->SetTitle("Entries");
  fOutputContainer->Add(fHistCentrality);

  fHistVertexDistribution = new TH1F("fHistVertexDistribution", "Primary vertex distribution", 40, -20., 20.);
  fHistVertexDistribution ->GetXaxis()->SetTitle("z_{v} (cm)");
  fHistVertexDistribution ->GetYaxis()->SetTitle("Entries");
  fOutputContainer->Add(fHistVertexDistribution); 

  fHistMultiplicityOfMixedEvent = new TH1F("fHistMultiplicityOfMixedEvent","",100, 0.,100.);
  fOutputContainer->Add(fHistMultiplicityOfMixedEvent);

  fHistnTPCCrossedR = new TH1F("fHistnTPCCrossedR","",200, 0.,200.);
  fOutputContainer->Add(fHistnTPCCrossedR);
  fHistRationTPCCrossedRnFind = new TH1F("fHistRationTPCCrossedRnFind","",2000, 0.,170.);
  fOutputContainer->Add(fHistRationTPCCrossedRnFind);
  fHistSharedFrTPCcl = new TH1F("fHistSharedFrTPCcl","",400, 0.,2.);
  fOutputContainer->Add(fHistSharedFrTPCcl);

  fHistprimpTPCdEdx = new TH2F("fHistprimpTPCdEdx", "", 400, -6.0, 6.0, 500, 0.0, 2000);
  fOutputContainer->Add(fHistprimpTPCdEdx);

  fHistyptProtons = new TH3F("fHistyptProtons","",100,0.,10.,40,-2.,2.,10,0.,100.);
  fOutputContainer->Add(fHistyptProtons);
  fHistphietaProtons = new TH2F("fHistphietaProtons","",400,0.,2*TMath::Pi(),100, -2, 2.);
  fOutputContainer->Add(fHistphietaProtons);
  fHistIPtoPVxyzTPC = new TH2F("fHistIPtoPVxyzTPC","",200, -5.,5., 200, -5., 5.);
  fOutputContainer->Add(fHistIPtoPVxyzTPC);
  fHistIPtoPVxyzGlobal = new TH2F("fHistIPtoPVxyzGlobal","",200, -5.,5., 200, -5., 5.);
  fOutputContainer->Add(fHistIPtoPVxyzGlobal);

  
  fHistpTOFmisvsp = new TH2F("fHistpTOFmisvsp", "", 200, 0., 10., 101, 0.0, 1.01);
  fOutputContainer->Add(fHistpTOFmisvsp);
  fHistpTOFnsigmavsp = new TH2F("fHistpTOFnsigmavsp", "", 200, 0., 10., 100, -5., 5.);
  fOutputContainer->Add(fHistpTOFnsigmavsp);
  fHistpTOFsignalvsp = new TH2F("fHistpTOFsignalvsp", "", 200, 0., 10, 100,-3000,3000);//1000 , 0., 2.);
  fHistpTOFsignalvsp->GetYaxis()->SetTitle("t_{meas}-t_{0}-t_{expected} (ps)");
  fHistpTOFsignalvsp->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
  fOutputContainer->Add(fHistpTOFsignalvsp);
  fHistpnsigTOFnsigTPC = new TH3F("fHistpnsigTOFnsigTPC","",100, 0., 5., 100, -10., 10., 100, -10., 10);
  fOutputContainer->Add(fHistpnsigTOFnsigTPC);
  fHistpsignalTOFsignalTPC = new TH3F("fHistpsignalTOFsignalTPC","", 100, 0., 5., 100, -3000, 3000, 500, 0.0, 2000);
  fHistpsignalTOFsignalTPC->GetZaxis()->SetTitle("TPC dE/dx (a.u.)");
  fHistpsignalTOFsignalTPC->GetYaxis()->SetTitle("t_{meas}-t_{0}-t_{expected} (ps)");
  fHistpsignalTOFsignalTPC->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
  fOutputContainer->Add(fHistpsignalTOFsignalTPC);

  fHistProtonMultvsCent = new TH2F("fHistProtonMultvsCent","",400,0.,2000.,10,0.,100.);
  fOutputContainer->Add(fHistProtonMultvsCent);

  fHistMCPrimProtons = new TH2F("fHistMCPrimProtons", "", 30, 0., 5.,200, -3.,3.);
  fOutputContainer->Add(fHistMCPrimProtons);
  fHistMCFromWdecayProtons = new TH2F("fHistMCFromWdecayProtons", "", 30, 0., 5.,200, -3.,3.);
  fOutputContainer->Add(fHistMCFromWdecayProtons);
  fHistMCFromMaterialProtons = new TH2F("fHistMCFromMaterialProtons", "", 30, 0., 5.,200, -3.,3.);
  fOutputContainer->Add(fHistMCFromMaterialProtons);
  fHistMCOtherProtons = new TH2F("fHistMCOtherProtons", "", 30, 0., 5.,30, -3.,3.);
  fOutputContainer->Add(fHistMCOtherProtons);
  fHistMCPrimAProtons = new TH2F("fHistMCPrimAProtons", "", 30, 0., 5.,30, -3.,3.);
  fOutputContainer->Add(fHistMCPrimAProtons);
  fHistMCFromWdecayAProtons = new TH2F("fHistMCFromWdecayAProtons", "", 30, 0., 5.,200, -3.,3.);
  fOutputContainer->Add(fHistMCFromWdecayAProtons);
  fHistMCFromMaterialAProtons = new TH2F("fHistMCFromMaterialAProtons", "", 30, 0., 5.,200, -3.,3.);
  fOutputContainer->Add(fHistMCFromMaterialAProtons);
  fHistMCOtherAProtons = new TH2F("fHistMCOtherAProtons", "", 30, 0., 5.,30, -3.,3.);
  fOutputContainer->Add(fHistMCOtherAProtons);

  fHistpTPCdEdx = new TH2F("fHistpTPCdEdx", "", 400, -6.0, 6.0, 500, 0.0, 2000);
  fOutputContainer->Add(fHistpTPCdEdx);
  fHistnTPCdEdx = new TH2F("fHistnTPCdEdx", "", 400, -6.0, 6.0, 500, 0.0, 2000);
  fOutputContainer->Add(fHistnTPCdEdx);
  fHistbTPCdEdx = new TH2F("fHistbTPCdEdx", "", 400, -6.0, 6.0, 500, 0.0, 2000);
  fOutputContainer->Add(fHistbTPCdEdx);

  fHistPosV0TPCClusters = new TH1F("fHistPosV0TPCClusters", "", 100, 0.,200.);
  fOutputContainer->Add(fHistPosV0TPCClusters);
  fHistNegV0TPCClusters = new TH1F("fHistNegV0TPCClusters", "", 100, 0.,200.);
  fOutputContainer->Add(fHistNegV0TPCClusters);
  fHistBachTPCClusters = new TH1F("fHistBachTPCClusters", "", 100, 0.,200.);
  fOutputContainer->Add(fHistBachTPCClusters);

  fHistSharedFrTPCclPos = new TH1F("fHistSharedFrTPCclPos","",400, 0.,2.);
  fOutputContainer->Add(fHistSharedFrTPCclPos);
  fHistSharedFrTPCclNeg = new TH1F("fHistSharedFrTPCclNeg","",400, 0.,2.); 
  fOutputContainer->Add(fHistSharedFrTPCclNeg);
  fHistSharedFrTPCclBach = new TH1F("fHistSharedFrTPCclBach","",400, 0.,2.);
  fOutputContainer->Add(fHistSharedFrTPCclBach);

  fHistyptXi = new TH2F("fHistyptXi","",100,0.,10.,40,-2.,2.);
  fOutputContainer->Add(fHistyptXi);
  fHistphietaXi = new TH2F("fHistphietaXi","",400,0.,2*TMath::Pi(),100, -2, 2.);
  fOutputContainer->Add(fHistphietaXi);

  fHistptL = new TH1F("fHistptL","",100,0.,10.);
  fOutputContainer->Add(fHistptL);
  fHistptpL = new TH1F("fHistptpL","",100,0.,10.);
  fOutputContainer->Add(fHistptpL);
  fHistptnL = new TH1F("fHistptnL","",100,0.,10.);
  fOutputContainer->Add(fHistptnL);
  fHistptbac = new TH1F("fHistptbac","",100,0.,10.);
  fOutputContainer->Add(fHistptbac);

  fHistInvMassXiMinus = new TH2F("fHistInvMassXiMinus","", fNbinsMassGm, fLowMassLimGm, fUpMassLimGm,100,0.,10.);
  fOutputContainer->Add(fHistInvMassXiMinus);

  fHistInvMassL = new TH2F("fHistInvMassL","", fNbinsMass, fLowMassLim, fUpMassLim,100, 0., 100.);
  fOutputContainer->Add(fHistInvMassL);
  fHistInvMassXiPlus = new TH2F("fHistInvMassXiPlus","", fNbinsMassGm, fLowMassLimGm, fUpMassLimGm,100,0.,10.);
  fOutputContainer->Add(fHistInvMassXiPlus);
  fHistInvMassAntiL = new TH2F("fHistInvMassAntiL","", fNbinsMass, fLowMassLim, fUpMassLim,100, 0., 100.);
  fOutputContainer->Add(fHistInvMassAntiL);
  fHistXiMultvsCent = new TH2F("fHistXiMultvsCent","",200,0.,200.,10,0.,100.);
  fOutputContainer->Add(fHistXiMultvsCent);

  fHistInvMassXiInPairs = new TH2F("fHistInvMassXiInPairs","", fNbinsMassGm, fLowMassLimGm, fUpMassLimGm,20,0.,100.);
  fOutputContainer->Add(fHistInvMassXiInPairs);

  fHistIPtoPVxyGlobalvspt = new TH2F("fHistIPtoPVxyGlobalvspt", "", 30, 0., 5.,200, -3.,3.);
  fOutputContainer->Add(fHistIPtoPVxyGlobalvspt);

/*  if(! fCFContCascadeCuts) {

   // Container meant to store all the relevant distributions corresponding to the cut variables.
   // - NB overflow/underflow of variables on which we want to cut later should be 0!!!
   const Int_t lNbSteps = 2 ;
   const Int_t lNbVariables = 22 ;
   //array for the number of bins in each dimension :
   Int_t lNbBinsPerVar[lNbVariables] = {0};
   lNbBinsPerVar[0] = 100;
   lNbBinsPerVar[1] = 126;
   lNbBinsPerVar[2] = 24;
   lNbBinsPerVar[3] = 220;
   lNbBinsPerVar[4] = 30;
   lNbBinsPerVar[5] = 50;
  
   lNbBinsPerVar[6] = 101;
  
   lNbBinsPerVar[7] = 102;
   lNbBinsPerVar[8] = 101;
   lNbBinsPerVar[9] = 26;
   lNbBinsPerVar[10] = 26;

   lNbBinsPerVar[11] = 150; // 75 2-MeV/c2 bins
   lNbBinsPerVar[12] = 120; // 60 2-MeV/c2 bins

   lNbBinsPerVar[13] = 100;

   lNbBinsPerVar[14] = 44; // 0.05 in rapidity units
   lNbBinsPerVar[15] = 100; // 44 0.05 in rapidity units

   lNbBinsPerVar[16] = 20;

   lNbBinsPerVar[17] = 11;
   lNbBinsPerVar[18] = 100;
   lNbBinsPerVar[19] = 112; // Proper time of cascade
   lNbBinsPerVar[20] = 112; // Proper time of V0
   lNbBinsPerVar[21] = 112; // Distance V0-Xi in transverse plane

   fCFContCascadeCuts = new AliCFContainer("fCFContCascadeCuts","Container for Cascade cuts", lNbSteps, lNbVariables, lNbBinsPerVar );


   //setting the bin limits

   //0
   fCFContCascadeCuts->SetBinLimits(0, 0., 2.); // DcaXiDaughters : 0.0 to 2.0
   //1
   Double_t *lBinLim1 = new Double_t[ lNbBinsPerVar[1]+1 ];
   lBinLim1[0] = 0.0;
   lBinLim1[1] = 0.03;
   for(Int_t i=2; i< lNbBinsPerVar[1];i++) lBinLim1[i] = (Double_t)0.03 + (5. - 0.03 )/(lNbBinsPerVar[1]-2) * (Double_t)(i-1) ;
   lBinLim1[ lNbBinsPerVar[1] ] = 100.0;
   fCFContCascadeCuts -> SetBinLimits(1, lBinLim1 );
   delete [] lBinLim1; // DcaBachToPrimVertexXi : 0.0 to 0.5
   //2
   fCFContCascadeCuts->SetBinLimits(2, .9988, 1.); // XiCosineOfPointingAngle : 0.99 to 1.0
   //3
   fCFContCascadeCuts -> SetBinLimits(3, 0., 110.); // XiRadius : 0.0 to 110.0
   //4
   fCFContCascadeCuts->SetBinLimits(4, 1.1, 1.13 ); // InvMassLambdaAsCascDghter
   //5
   fCFContCascadeCuts->SetBinLimits(5, 0., 2.); // DcaV0DaughtersXi : 0.0 to 2.0
   //6
   fCFContCascadeCuts->SetBinLimits(6, .98, 1.0002); // V0CosineOfPointingAngleXi : 0.99 to 1.0
   //7
   Double_t *lBinLim7 = new Double_t[ lNbBinsPerVar[7]+1 ];
   for(Int_t i=0; i< lNbBinsPerVar[7]-1;i++) lBinLim7[i] = (Double_t)0.0 + (100. - 0.0 )/(lNbBinsPerVar[7]-2) * (Double_t)i ;
   lBinLim7[ lNbBinsPerVar[7]-1] = 200.0;
   lBinLim7[ lNbBinsPerVar[7]] = 1000.0;
   fCFContCascadeCuts -> SetBinLimits(7, lBinLim7 );
   delete [] lBinLim7; // V0RadiusXi : 0.0 to 100.0
   //8
   Double_t *lBinLim8 = new Double_t[ lNbBinsPerVar[8]+1 ];
   for(Int_t i=0; i< lNbBinsPerVar[8];i++) lBinLim8[i] = (Double_t)0.0 + (0.4 - 0.0 )/(lNbBinsPerVar[8]-1) * (Double_t)i ;
   lBinLim8[ lNbBinsPerVar[8] ] = 100.0;
   fCFContCascadeCuts -> SetBinLimits(8, lBinLim8 );
   delete [] lBinLim8; // DcaV0ToPrimVertexXi : 0.0 to 0.4
   //9
   Double_t *lBinLim9 = new Double_t[ lNbBinsPerVar[9]+1 ];
   for(Int_t i=0; i< lNbBinsPerVar[9];i++) lBinLim9[i] = (Double_t)0.0 + (0.25 - 0.0 )/(lNbBinsPerVar[9]-1) * (Double_t)i ;
   lBinLim9[ lNbBinsPerVar[9] ] = 100.0;
   fCFContCascadeCuts -> SetBinLimits(9, lBinLim9 );
   delete [] lBinLim9; // DcaPosToPrimVertexXi : 0.0 to 0.25
   //10
   Double_t *lBinLim10 = new Double_t[ lNbBinsPerVar[10]+1 ];
   for(Int_t i=0; i< lNbBinsPerVar[10];i++) lBinLim10[i] = (Double_t)0.0 + (0.25 - 0.0 )/(lNbBinsPerVar[10]-1) * (Double_t)i ;
   lBinLim10[ lNbBinsPerVar[10] ] = 100.0;
   fCFContCascadeCuts -> SetBinLimits(10, lBinLim10 );
   delete [] lBinLim10; // DcaNegToPrimVertexXi : 0.0 to 0.25
   //11
   fCFContCascadeCuts->SetBinLimits(11, 1.25, 1.40); // InvMassXi
   fCFContCascadeCuts->SetBinLimits(12, 1.62, 1.74); // InvMassOmega
   fCFContCascadeCuts->SetBinLimits(13, 0.0, 10.0); // XiTransvMom
   fCFContCascadeCuts->SetBinLimits(14, -1.1, 1.1); // Y(Xi)
   fCFContCascadeCuts->SetBinLimits(15, 0.,10.);//-1.1, 1.1); // Y(Omega)
   fCFContCascadeCuts->SetBinLimits(16, -10.0, 10.0); // BestPrimaryVtxPosZ
   Double_t *lBinLim17 = new Double_t[ lNbBinsPerVar[17]+1 ];
   for(Int_t i=3; i< lNbBinsPerVar[17]+1;i++) lBinLim17[i] = (Double_t)(i-1)*10.;
   lBinLim17[0] = 0.0;
   lBinLim17[1] = 5.0;
   lBinLim17[2] = 10.0;
   fCFContCascadeCuts -> SetBinLimits(17, lBinLim17 ); // Centrality
   delete [] lBinLim17;
   fCFContCascadeCuts->SetBinLimits(18, 0.0, 6000.0); // Track multiplicity not used now
  
   Double_t *lBinLim19 = new Double_t[ lNbBinsPerVar[19]+1 ];
   for(Int_t i=0; i< lNbBinsPerVar[19];i++) lBinLim19[i] = (Double_t)-1. + (110. + 1.0 )/(lNbBinsPerVar[19]-1) * (Double_t)i ;
   lBinLim19[ lNbBinsPerVar[19] ] = 2000.0;
   fCFContCascadeCuts->SetBinLimits(19, lBinLim19); // Proper time cascade
  
   fCFContCascadeCuts->SetBinLimits(20, lBinLim19); // Proper time V0
  
   fCFContCascadeCuts->SetBinLimits(21, lBinLim19); // Distance V0-Xi in tansverse plane
   delete [] lBinLim19;
   // Setting the number of steps : one for each cascade species (Xi-, Xi+ and Omega-, Omega+)
   fCFContCascadeCuts->SetStepTitle(0, "#Xi^{-} candidates");
   fCFContCascadeCuts->SetStepTitle(1, "#bar{#Xi}^{+} candidates");
 //  fCFContCascadeCuts->SetStepTitle(2, "#Omega^{-} candidates");
//   fCFContCascadeCuts->SetStepTitle(3, "#bar{#Omega}^{+} candidates");
  
   // Setting the variable title, per axis
   fCFContCascadeCuts->SetVarTitle(0, "DCA(#Xi daughters) (cm)");
   fCFContCascadeCuts->SetVarTitle(1, "IP(Bachelor to prim. vertex) (cm)");
   fCFContCascadeCuts->SetVarTitle(2, "cos(#Xi pointing angle)");
   fCFContCascadeCuts->SetVarTitle(3, "R_{xy}(#Xi decay) (cm)");
   fCFContCascadeCuts->SetVarTitle(4, "M_{#Lambda}(as casc. daughter) (GeV/#it{c}^{2})");
   fCFContCascadeCuts->SetVarTitle(5, "DCA(#Lambda daughters) (cm)");
  
   fCFContCascadeCuts->SetVarTitle(6, "cos(#Lambda pointing angle)");
   fCFContCascadeCuts->SetVarTitle(7, "R_{xy}(#Lambda decay) (cm)");
   fCFContCascadeCuts->SetVarTitle(8, "IP(#Lambda to prim. vertex) (cm)");
   fCFContCascadeCuts->SetVarTitle(9, "IP(Pos to prim. vertex) (cm)");
   fCFContCascadeCuts->SetVarTitle(10, "IP(Neg to prim. vertex) (cm)");
  
   fCFContCascadeCuts->SetVarTitle(11, "M_{#Lambda#pi} (GeV/#it{c}^{2})");
   fCFContCascadeCuts->SetVarTitle(12, "M_{#LambdaK} (GeV/#it{c}^{2})");
  
   fCFContCascadeCuts->SetVarTitle(13, "p_{T} (GeV/#it{c})");

   fCFContCascadeCuts->SetVarTitle(14, "y_{#Xi}");
   fCFContCascadeCuts->SetVarTitle(15, "#chi^{2} KFparticle fit");//Y(Omega)"); // not implemented 

   fCFContCascadeCuts->SetVarTitle(16, "Z-position of prim. vertex (cm)");

   fCFContCascadeCuts->SetVarTitle(17, "Centrality");

   fCFContCascadeCuts->SetVarTitle(18, "Event track multiplicity");

   fCFContCascadeCuts->SetVarTitle(19, "Proper time cascade");
   fCFContCascadeCuts->SetVarTitle(20, "Proper time #Lambda");
   fCFContCascadeCuts->SetVarTitle(21, "Distance xy #Lambda-cascade ");
  
 }
*/  

  // Histos for CF in k*
  int kstarbins = 0;
  float kstarmax = 5.;
  float kstarmin = 0.;
  if (fSecondpart == kPion) {kstarbins = 500; kstarmax = 2.5; kstarmin = 0.0;}
  else {kstarbins = 500; kstarmax=5; }

  fHistpXiSignalRealKstar = new TH2D("fHistpXiSignalRealKstar","",kstarbins,kstarmin,kstarmax,20,0.,100.);
  fOutputContainer->Add(fHistpXiSignalRealKstar);
  fHistapXiSignalRealKstar = new TH2D("fHistapXiSignalRealKstar","",kstarbins,kstarmin,kstarmax,20,0.,100.);
  fOutputContainer->Add(fHistapXiSignalRealKstar);
  fHistpaXiSignalRealKstar = new TH2D("fHistpaXiSignalRealKstar","",kstarbins,kstarmin,kstarmax,20,0.,100.);
  fOutputContainer->Add(fHistpaXiSignalRealKstar);
  fHistapaXiSignalRealKstar = new TH2D("fHistapaXiSignalRealKstar","",kstarbins,kstarmin,kstarmax,20,0.,100.);
  fOutputContainer->Add(fHistapaXiSignalRealKstar);
  fHistpXiSignalMixedKstargenvsrec = new TH2D("fHistpXiSignalMixedKstargenvsrec","",999,0.,0.9995,999,0.,0.9995);
  fOutputContainer->Add(fHistpXiSignalMixedKstargenvsrec);
  fHistapXiSignalMixedKstargenvsrec = new TH2D("fHistapXiSignalMixedKstargenvsrec","",999,0.,0.9995,999,0.,0.9995);
  fOutputContainer->Add(fHistapXiSignalMixedKstargenvsrec);
  fHistpaXiSignalMixedKstargenvsrec = new TH2D("fHistpaXiSignalMixedKstargenvsrec","",999,0.,0.9995,999,0.,0.9995);
  fOutputContainer->Add(fHistpaXiSignalMixedKstargenvsrec);
  fHistapaXiSignalMixedKstargenvsrec = new TH2D("fHistapaXiSignalMixedKstargenvsrec","",999,0.,0.9995,999,0.,0.9995);
  fOutputContainer->Add(fHistapaXiSignalMixedKstargenvsrec);

  fHistpXiSignalBkgKstar = new TH2D("fHistpXiSignalBkgKstar","",kstarbins,kstarmin,kstarmax,20,0.,100.);
  fOutputContainer->Add(fHistpXiSignalBkgKstar);
  fHistapXiSignalBkgKstar = new TH2D("fHistapXiSignalBkgKstar","",kstarbins,kstarmin,kstarmax,20,0.,100.);
  fOutputContainer->Add(fHistapXiSignalBkgKstar);
  fHistpaXiSignalBkgKstar = new TH2D("fHistpaXiSignalBkgKstar","",kstarbins,kstarmin,kstarmax,20,0.,100.);
  fOutputContainer->Add(fHistpaXiSignalBkgKstar);
  fHistapaXiSignalBkgKstar = new TH2D("fHistapaXiSignalBkgKstar","",kstarbins,kstarmin,kstarmax,20,0.,100.);
  fOutputContainer->Add(fHistapaXiSignalBkgKstar);

  fHistFractionOfXiWithSharedDaughters = new TH2F("fHistFractionOfXiWithSharedDaughters","",200,0.,200.,201,0.,1.005);
  fOutputContainer->Add(fHistFractionOfXiWithSharedDaughters);
  fHistTotMaxFractionOfXiWithSharedDaughters = new TH2F("fHistTotMaxFractionOfXiWithSharedDaughters","",200,0.,200.,201,0.,1.005);
  fOutputContainer->Add(fHistTotMaxFractionOfXiWithSharedDaughters);

  // Histos for CF in 
  int detabins = 0;
  float detamin = -2.;
  float detamax = 2.;
  int dphibins = 1000;
  float dphimin = -1.;
  float dphimax = 1.;
  if (fSecondpart == kPion) {detabins = 200; detamin = -0.5; detamax = 0.5; dphibins = 500; dphimin = -.5; dphimax = .5; }
  else if (fSecondpart == kProton) detabins = 400;
  else detabins = 200;

  fHistpXibacDetaSDphiS = new TH2D("fHistpXibacDetaSDphiS","",dphibins, dphimin, dphimax, detabins, detamin, detamax);    
  fOutputContainer->Add(fHistpXibacDetaSDphiS);
  fHistpXiposDetaSDphiS = new TH2D("fHistpXiposDetaSDphiS","",dphibins, dphimin, dphimax, detabins, detamin, detamax);
  fOutputContainer->Add(fHistpXiposDetaSDphiS);
  fHistpXinegDetaSDphiS = new TH2D("fHistpXinegDetaSDphiS","",dphibins, dphimin, dphimax, detabins, detamin, detamax);
  fOutputContainer->Add(fHistpXinegDetaSDphiS);
  fHistpaXibacDetaSDphiS = new TH2D("fHistpaXibacDetaSDphiS","",dphibins, dphimin, dphimax, detabins, detamin, detamax);
  fOutputContainer->Add(fHistpaXibacDetaSDphiS);
  fHistpaXiposDetaSDphiS = new TH2D("fHistpaXiposDetaSDphiS","",dphibins, dphimin, dphimax, detabins, detamin, detamax);
  fOutputContainer->Add(fHistpaXiposDetaSDphiS);
  fHistpaXinegDetaSDphiS = new TH2D("fHistpaXinegDetaSDphiS","",dphibins, dphimin, dphimax, detabins, detamin, detamax);
  fOutputContainer->Add(fHistpaXinegDetaSDphiS);
  fHistapXibacDetaSDphiS = new TH2D("fHistapXibacDetaSDphiS","",dphibins, dphimin, dphimax, detabins, detamin, detamax);
  fOutputContainer->Add(fHistapXibacDetaSDphiS);
  fHistapXiposDetaSDphiS = new TH2D("fHistapXiposDetaSDphiS","",dphibins, dphimin, dphimax, detabins, detamin, detamax);
  fOutputContainer->Add(fHistapXiposDetaSDphiS);
  fHistapXinegDetaSDphiS = new TH2D("fHistapXinegDetaSDphiS","",dphibins, dphimin, dphimax, detabins, detamin, detamax);
  fOutputContainer->Add(fHistapXinegDetaSDphiS);
  fHistapaXibacDetaSDphiS = new TH2D("fHistapaXibacDetaSDphiS","",dphibins, dphimin, dphimax, detabins, detamin, detamax);
  fOutputContainer->Add(fHistapaXibacDetaSDphiS);
  fHistapaXiposDetaSDphiS = new TH2D("fHistapaXiposDetaSDphiS","",dphibins, dphimin, dphimax, detabins, detamin, detamax);
  fOutputContainer->Add(fHistapaXiposDetaSDphiS);
  fHistapaXinegDetaSDphiS = new TH2D("fHistapaXinegDetaSDphiS","",dphibins, dphimin, dphimax, detabins, detamin, detamax);
  fOutputContainer->Add(fHistapaXinegDetaSDphiS);

  fHistpXibacDetaSDphiSBkg = new TH2D("fHistpXibacDetaSDphiSBkg","",dphibins, dphimin, dphimax, detabins, detamin, detamax);
  fOutputContainer->Add(fHistpXibacDetaSDphiSBkg);
  fHistpXiposDetaSDphiSBkg = new TH2D("fHistpXiposDetaSDphiSBkg","",dphibins, dphimin, dphimax, detabins, detamin, detamax);
  fOutputContainer->Add(fHistpXiposDetaSDphiSBkg);
  fHistpXinegDetaSDphiSBkg = new TH2D("fHistpXinegDetaSDphiSBkg","",dphibins, dphimin, dphimax, detabins, detamin, detamax);
  fOutputContainer->Add(fHistpXinegDetaSDphiSBkg);
  fHistpaXibacDetaSDphiSBkg = new TH2D("fHistpaXibacDetaSDphiSBkg","",dphibins, dphimin, dphimax, detabins, detamin, detamax);
  fOutputContainer->Add(fHistpaXibacDetaSDphiSBkg);
  fHistpaXiposDetaSDphiSBkg = new TH2D("fHistpaXiposDetaSDphiSBkg","",dphibins, dphimin, dphimax, detabins, detamin, detamax);
  fOutputContainer->Add(fHistpaXiposDetaSDphiSBkg);
  fHistpaXinegDetaSDphiSBkg = new TH2D("fHistpaXinegDetaSDphiSBkg","",dphibins, dphimin, dphimax, detabins, detamin, detamax);
  fOutputContainer->Add(fHistpaXinegDetaSDphiSBkg);
  fHistapXibacDetaSDphiSBkg = new TH2D("fHistapXibacDetaSDphiSBkg","",dphibins, dphimin, dphimax, detabins, detamin, detamax);
  fOutputContainer->Add(fHistapXibacDetaSDphiSBkg);
  fHistapXiposDetaSDphiSBkg = new TH2D("fHistapXiposDetaSDphiSBkg","",dphibins, dphimin, dphimax, detabins, detamin, detamax);
  fOutputContainer->Add(fHistapXiposDetaSDphiSBkg);
  fHistapXinegDetaSDphiSBkg = new TH2D("fHistapXinegDetaSDphiSBkg","",dphibins, dphimin, dphimax, detabins, detamin, detamax);
  fOutputContainer->Add(fHistapXinegDetaSDphiSBkg);
  fHistapaXibacDetaSDphiSBkg = new TH2D("fHistapaXibacDetaSDphiSBkg","",dphibins, dphimin, dphimax, detabins, detamin, detamax);
  fOutputContainer->Add(fHistapaXibacDetaSDphiSBkg);
  fHistapaXiposDetaSDphiSBkg = new TH2D("fHistapaXiposDetaSDphiSBkg","",dphibins, dphimin, dphimax, detabins, detamin, detamax);
  fOutputContainer->Add(fHistapaXiposDetaSDphiSBkg);
  fHistapaXinegDetaSDphiSBkg = new TH2D("fHistapaXinegDetaSDphiSBkg","",dphibins, dphimin, dphimax, detabins, detamin, detamax);
  fOutputContainer->Add(fHistapaXinegDetaSDphiSBkg);

  fHistTrackBufferOverflow = new TH1F("fHistTrackBufferOverflow","",2,0,2);
  fOutputContainer->Add(fHistTrackBufferOverflow);

  PostData(1, fOutputContainer );
//  PostData(2, fCFContCascadeCuts);

}

//_______________________________________________________________

void AliAnalysisTaskhCascadeFemto::UserExec(Option_t *) {
  // Main loop
  // Called for each event

  if(!fPIDResponse) {
    AliError("Cannot get pid response");
    return;
  }

  AliVVertex *vertexmain =0x0;
  AliCentrality* centrality = 0x0;
  Double_t lBestPrimaryVtxPos[3] = {-100.0, -100.0, -100.0};


  AliMCEvent  *lMCevent  = 0x0;
  AliStack    *lMCstack  = 0x0;
  TClonesArray *arrayMC = 0x0;

  Int_t ncascades = 0;
  Int_t ntracks = 0;

  if (fAnalysisType == "ESD") {
    fESDevent = dynamic_cast<AliESDEvent*>( InputEvent() );
    if (!fESDevent) {
      AliWarning("ERROR: ESDevent not available \n");
      return;
    }

    if (fReadMCTruth) {
      lMCevent = MCEvent();
      if (!lMCevent) {
        Printf("ERROR: Could not retrieve MC event \n");
        //cout << "Name of the file with pb :" <<  CurrentFileName() << endl;
        return;
      }
      lMCstack = lMCevent->Stack();
      if (!lMCstack) {
        Printf("ERROR: Could not retrieve MC stack \n");
        //cout << "Name of the file with pb :" <<  CurrentFileName() << endl;
        return;
      }
    }

    if((TMath::Abs(lBestPrimaryVtxPos[2])) > 10.) return;
    ncascades = fESDevent->GetNumberOfCascades();
    ntracks = fESDevent->GetNumberOfTracks(); 
    centrality = fESDevent->GetCentrality();

    const AliESDVertex *lPrimaryBestESDVtx = fESDevent->GetPrimaryVertex();
    if (!lPrimaryBestESDVtx){
      AliWarning("No prim. vertex in ESD... return!");
      return;
    }
    vertexmain = (AliVVertex*) lPrimaryBestESDVtx;
    lPrimaryBestESDVtx->GetXYZ( lBestPrimaryVtxPos );

  } else if (fAnalysisType == "AOD") {
    fAODevent = dynamic_cast<AliAODEvent*>( InputEvent() );
    if (!fAODevent) {
      AliWarning("ERROR: AODevent not available \n");
      return;
    }

    ncascades  = fAODevent->GetNumberOfCascades();
    ntracks = fAODevent->GetNumberOfTracks();
    centrality = fAODevent->GetCentrality();

    const AliAODVertex *lPrimaryBestAODVtx = fAODevent->GetPrimaryVertex();
    if (!lPrimaryBestAODVtx){
      AliWarning("No prim. vertex in AOD... return!");
      return;
    }
    vertexmain = (AliVVertex*) lPrimaryBestAODVtx;
    lPrimaryBestAODVtx->GetXYZ( lBestPrimaryVtxPos );

    if (fReadMCTruth) {
//      Printf("Reading MC truth!!! \n");
      arrayMC = (TClonesArray*) fAODevent->GetList()->FindObject(AliAODMCParticle::StdBranchName());

      if (!arrayMC) AliFatal("Error: MC particles branch not found!\n");

    }

  } else {

    Printf("Analysis type (ESD or AOD) not specified \n");
    return;

  }

  Float_t lcentrality = centrality->GetCentralityPercentile("V0M");
  if (lcentrality<fCentrLowLim||lcentrality>=fCentrUpLim)  return;
  if((TMath::Abs(lBestPrimaryVtxPos[2])) > 10.) return;

  fHistCentrality->Fill(lcentrality);
  fHistVertexDistribution->Fill(lBestPrimaryVtxPos[2]); 

  const Float_t bfield = (InputEvent())->GetMagneticField();
  //cout<<"B field"<<bfield<<endl;;
  int fieldsign;
  if (bfield >=0.) fieldsign = 1;
  else fieldsign = -1;

  // Store the event in the buffer to do mixing
  // ... find vertex... 
  int zBin=0;

  double zStep=2*10/double(fzVertexBins), zStart=-10.;

  for (int i=0; i<fzVertexBins; i++) {
    if ((lBestPrimaryVtxPos[2] > zStart+i*zStep) && (lBestPrimaryVtxPos[2] < zStart+(i+1)*zStep)) {
      zBin=i;
      break;
    }
  }

  // ... and centrality
  int centralityBin=0;

  if(lcentrality < 5.) centralityBin=19;  // changed <= with < to be consistent with histogram binning, except last bin 
  else if(lcentrality < 10.) centralityBin=18;
  else if(lcentrality < 15.) centralityBin=17;
  else if(lcentrality < 20.) centralityBin=16;
  else if(lcentrality < 25.) centralityBin=15;
  else if(lcentrality < 30.) centralityBin=14;
  else if(lcentrality < 35.) centralityBin=13;
  else if(lcentrality < 40.) centralityBin=12;
  else if(lcentrality < 45.) centralityBin=11;
  else if(lcentrality < 50.) centralityBin=10;
  else if(lcentrality < 55.) centralityBin=9; 
  else if(lcentrality < 60.) centralityBin=8;
  else if(lcentrality < 65.) centralityBin=7;
  else if(lcentrality < 70.) centralityBin=6;
  else if(lcentrality < 75.) centralityBin=5;
  else if(lcentrality < 80.) centralityBin=4;
  else if(lcentrality < 85.) centralityBin=3;
  else if(lcentrality < 90.) centralityBin=2;   // from here on, wont be filled because the range selected in AddTask is 0-90
  else if(lcentrality < 95.) centralityBin=1;
  else if(lcentrality <= 100.) centralityBin=0;

  fEventColl[zBin][centralityBin]->FifoShift();

  fEvt = fEventColl[zBin][centralityBin]->fEvt;

//  printf("buffer size: %d\n",fTrackBufferSize);
//  printf("ntracks: %d\n",ntracks);

  for (Int_t igt = 0; igt < fTrackBufferSize; igt++) farrGT[igt] = -1;

  AliAODTrack *globaltrack = 0x0;
 
  // Read and store global tracks to retrieve PID information for TPC only tracks
  for (Int_t igt = 0; igt < ntracks; igt++) {
    globaltrack = (AliAODTrack*) fAODevent->GetTrack(igt);

    if (!globaltrack) continue; 
    if (globaltrack->GetID()<0 ) continue;
    //cout<<" Global track!!  "<<igt<<endl;

    if ( globaltrack->GetTPCNcls()<=0   ) { // such tracks do not have the TPC refit either, filter map is 2 --> ITS constrained
      //cout<<" No TPC cl for this global track!!  "<<igt<<endl;
      //bool statusTPC = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC,globaltrack);
      //if (AliPIDResponse::kDetPidOk == statusTPC) cout<<"But has TPC PID!!!!!"<<endl; 
      //else cout<<"And no TPC PID info"<<endl; // always when no TPC clusters are present
      //if (!globaltrack->IsOn(AliAODTrack::kTPCrefit)) cout<<" ... and also no tpc refit  "<<globaltrack->GetFilterMap()<<endl;
      continue;
    }

    //if (!globaltrack->IsOn(AliAODTrack::kTPCrefit)) { cout<<" No tpc refit  "<<globaltrack->GetFilterMap()<<endl; continue;} // there are such tracks with no TPC clusters

//      if ( !(globaltrack->GetFilterMap()) ) { 
//        cout<<" No filter map for this global track!!  "<<globaltrack->GetFilterMap()<<endl;
//        continue;
//    }


    // Check id is not too big for buffer
    if (globaltrack->GetID()>=fTrackBufferSize) {
      printf("Warning: track ID too big for buffer: ID: %d, buffer %d\n",globaltrack->GetID(),fTrackBufferSize);
      fHistTrackBufferOverflow->Fill(1);
      continue;
    }

    //cout<<" The array will contain "<<igt<<" , it contains "<<farrGT[globaltrack->GetID()]<<endl; 

    // Warn if we overwrite a track
    if (farrGT[globaltrack->GetID()]>=0) { // Two tracks same id  --> checked that never happens
      cout<<" Array already filled "<<farrGT[globaltrack->GetID()]<<endl; // Warning but we dont overwrite
    } else { 
      farrGT[globaltrack->GetID()] = igt;           // solution adopted in the femto framework, Hans uses pointers
//    cout<<" Set array, now it contains  "<<farrGT[globaltrack->GetID()]<<endl;
    }
  }
 
  globaltrack = 0x0; 

  // Read, plot and store pions or protons
  
  AliAODTrack *track = 0x0;
  AliVTrack *vtrackg = 0x0;
  AliVTrack *vtrack = 0x0;

  Bool_t isTOFPIDok = kFALSE;
  Float_t nsigmaTOFp = 0.;
  Float_t nsigmaTPCp = 0.;
  Float_t probMis = 0.;
  Double32_t tTOF = 0.;
  Double_t beta = 0.; 
  Double_t length = 0.;

  Float_t nTPCCrossedRows = 0.; 
  Float_t rationCrnFind = 0.;

  UShort_t nTPCclSp = 0;
  UShort_t nTPCclp = 0;

  Float_t sharedFractionTPCcl = 0.;
  AliPIDResponse::EDetPidStatus statusTOF;
  AliPIDResponse::EDetPidStatus statusTPC;
 
  Double_t expectedTimes[AliPID::kSPECIES];

  Double_t dz[2]= {-999.,-999.}; //Double_t covar[3]={-999.,-999.,-999.};
  Double_t dzg[2]= {-999.,-999.}; Double_t covarg[3]={-999.,-999.,-999.};
  AliExternalTrackParam etp1;
  Double_t xyz[3] = {0.,0.,0.};
 

  Float_t rapidity = 0.;
  Short_t charge = -2;
  Float_t trackmomentum = 0.;
  Float_t tracktransversemomentum = 0.;
  Float_t tpcSignal = 0.;


  int pCount = 0;

  double pMomentumTruth[3] = {0.,0.,0.};
  AliReconstructedProton::MCProtonOrigin_t mcProtonOrigin = AliReconstructedProton::kUnassigned;

  Bool_t isP = kFALSE;  // particle
  Bool_t isaP = kFALSE; // anti-particle 

  for (Int_t ip = 0; ip < ntracks; ip++) {

    vtrack = fAODevent->GetTrack(ip);
    if (!vtrack) continue;
    track = dynamic_cast<AliAODTrack*>(vtrack);
    if(!track) AliFatal("Not a standard AOD");

    // Take TPC only tracks constrained to the SPD PV during filtering 
    // Filter Bit 7 (1<<7) (Filter Mask 128) but PID is not stored --> find the corresponding global track 
    if(!track->TestFilterBit(AliAODTrack::kTrkTPCOnlyConstrained)) continue; // Bit 7 corresponds to the following cuts (applied to original ESD track, see Ruben presentation 13.02.2014) : 
                                                                             // ncl 50  
                                                                             // dcaxy 2.4 dcaz 3.2 // this cut however is on the global track!! PWGCF meeting (17.02.2014)
                                                                             // no kinks 
                                                                             // chi2percl 4 
//  if (track->TestFilterBit(1 << 7)) cout<<" Filter bit of those traks is 7 "<<endl;
//  if (track->TestFilterBit(128)) cout<<" Filter bit of those traks is 7  "<<endl;
//  if (track->P()==0.) cout<<" track has momentum zero!!"<<endl; // never

//    nTPCCrossedRows = track->GetTPCClusterInfo(2, 1); // The method below is equivalent
    nTPCCrossedRows = track->GetTPCNCrossedRows();
    nTPCclSp =  track->GetTPCnclsS();
    nTPCclp = track->GetTPCncls();

    //cout<<"no of crossed rows "<<nTPCCrossedRows<<endl;
   
    if (fkCutOnCrossedRows) {
      if(nTPCCrossedRows<70) continue;  // 0 in MC 2010 AOD
    } else {
      if (nTPCclp<70) continue; 
    }

    rationCrnFind = 0.;
    //if (!track->GetTPCNclsF()) {cout<<" findable n cl 0 !!"<<endl; continue;}
    if (track->GetTPCNclsF()) { rationCrnFind = nTPCCrossedRows/track->GetTPCNclsF(); /*cout<<" ratio cr find "<<rationCrnFind<<endl;*/  } // not there in 2010 AOD
    if (fkApplyRatioCrRnFindCut) if(rationCrnFind< .8) continue;  

    sharedFractionTPCcl = (Float_t) nTPCclSp/nTPCclp; // can be > 1.
    //cout<<"Shared fraction ncl "<<sharedFractionTPCcl<<endl;
 
    if (sharedFractionTPCcl>0.) continue; // for protons it is around 4% (after all cuts) // FIXME Hans has a different way of retrieving it

    // Get the corresponding global track to use PID --> stored only for global tracks
    // 
    // Check that the array fGTI isn't too small
    // for the track id
    if (-track->GetID()-1 >= fTrackBufferSize) {
      printf ("Exceeding buffer size!!");
      continue;
    }

    vtrackg = fAODevent->GetTrack(farrGT[-track->GetID()-1]);
    if (!vtrackg) {
      printf ("No global info! iTrack %d, ID %d\n",ip,track->GetID());
      continue;
    }
    if (farrGT[-track->GetID()-1]>=ntracks || farrGT[-track->GetID()-1]<0) { 
      //cout<<"This index is out of range!!"<<farrGT[-track->GetID()-1]<<endl; 
      continue;
    }
    globaltrack = dynamic_cast<AliAODTrack*>(vtrackg); 
    if(!globaltrack) AliFatal("Not a standard AOD");

//    cout<<" Filter map for the global track "<<globaltrack->GetFilterMap()<<endl;
 

    // IP to PV of tracks
    dz[0]= -999.; dz[1] = -999.;
    dzg[0]= -999.; dzg[1] = -999.;

    if (track->TestBit(AliAODTrack::kIsDCA)) {  // always

      dz[0] = track->DCA();    // the TPC one should be applied the other biases the CF --> from Maciejs note, checked YES!  
      dz[1] = track->ZAtDCA(); // for those two lines check AliAODTrack.h // FIXME shifted distributions, known problem, asked Mac and Marian. 
    } //else cout<<"Bit is not DCA!!!!!!!!!!!"<<endl;                                                                      
//    Double_t p[3]; track->GetXYZ(p); // first two elements of p give the same as above, ignore the third, those are original DCA of TPC only tracks (with or w/o constraint?)
//    cout<<"d xy "<<dz[0]<<"while value is for other "<<p[0]<<endl;
//    cout<<"d z "<<dz[1]<<"while value is for other "<<p[1]<<endl;


    //if (dz[0]==0. || dz[1]==0. ) cout<<"dcaxy "<<dz[0]<<" dcaz "<<dz[1]<<endl; // never

    isP = kFALSE;
    isaP = kFALSE;
    
    charge = globaltrack->Charge();
    trackmomentum = track->P();
    tracktransversemomentum = track->Pt();
 
    if (TMath::Abs(track->Eta())> 0.8) continue;
    if (track->Chi2perNDF()>4.) continue; // should be redundant already applied in the TPC only track ESD filtering
    if ((tracktransversemomentum<fMinPtForPrim) || (tracktransversemomentum>fMaxPtForPrim)) continue;  

/*    if (fkCutOnTPCIP) {
      if (TMath::Abs(dz[0])>fIPCutxy ) continue;  // 2.4 proton 1. pion 
      if (TMath::Abs(dz[1])>fIPCutz ) continue;   // 3.2  proton 1. pion
    } else {
      // in AliAODTrack::PropagateToDCA() the lines below are iplemented 
      globaltrack->GetPosition(xyz);
      if (xyz[0]*xyz[0]+xyz[1]*xyz[1] >3.*3.) {
        //cout<<" Allowed only for tracks inside the beam pipe!!"<<endl;
         continue;
      } 
      etp1.CopyFromVTrack(vtrackg);
      if (!etp1.PropagateToDCA(vertexmain,bfield,100.,dzg,covarg)) { // Does not work for TPC constrained
        //cout<<"wrong propagation!!! "<<dzg[0]<<" and  "<<dzg[1]<<endl;
        continue;
      } 
      if (TMath::Abs(dzg[0])>fIPCutxy ) continue;  // 0.1 proton/2 pion
      if (TMath::Abs(dzg[1])>fIPCutz ) continue;   // 0.15 proton/pion
    }
*/

    if (fFirstpart == kPion) { 
      if (charge>0) isaP = kTRUE;
      else if (charge<0) isP = kTRUE;
    } else if (fFirstpart == kProton) {
      if (charge>0) isP = kTRUE;
      else if (charge<0) isaP = kTRUE;
    }

    rapidity = 0.5*TMath::Log( (track->E(fPDGfirst) + track->Pz()) / (track->E(fPDGfirst) - track->Pz() +1.e-13) );

    nsigmaTOFp = 0.; // be careful with those initialization
    nsigmaTPCp = 0.;
    probMis = 0.;

    statusTPC = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC,globaltrack);

    if (AliPIDResponse::kDetPidOk != statusTPC) { 
      //printf ("TPC PID not there"); 
      continue; 
    }

    tpcSignal = globaltrack->GetTPCsignal();

    if (fFirstpart == kPion) nsigmaTPCp = fPIDResponse->NumberOfSigmasTPC(globaltrack, AliPID::kPion); 
    else if (fFirstpart == kProton) nsigmaTPCp = fPIDResponse->NumberOfSigmasTPC(globaltrack, AliPID::kProton); 
   
    if (fFirstpart == kProton) { // For protons use TPC and TOF

      statusTOF = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF,globaltrack); // this checks kTOFout and kTIME https://twiki.cern.ch/twiki/bin/viewauth/ALICE/TOF

      tTOF = 0.0;
      beta = -1000.;
      length = 0.;
      isTOFPIDok = kFALSE;

      if ( (statusTOF == AliPIDResponse::kDetPidOk) ) {  
        isTOFPIDok = kTRUE;
        probMis = fPIDResponse->GetTOFMismatchProbability(globaltrack);
        nsigmaTOFp = fPIDResponse->NumberOfSigmasTOF(globaltrack, AliPID::kProton); 
        tTOF = globaltrack->GetTOFsignal(); 
        tTOF -= fPIDResponse->GetTOFResponse().GetStartTime(globaltrack->P());      
        globaltrack->GetIntegratedTimes(expectedTimes);
        tTOF -= expectedTimes[AliPID::kProton];  

        // length = globaltrack->GetIntegratedLength();  // FIXME length is zero!! from a mail february 2014: this info is not available for AODs 115
                                                         // use AODs 145
        // if (tTOF > 0.) {//&&length>350) {
        // beta = length/(tTOF*2.99792457999999984e-02); 
        // cout<<" rack length  "<<length<<" beta "<<beta<<endl; 
        // gamma = 1/TMath::Sqrt(1 - beta*beta);
        // mass = ptot/TMath::Sqrt(gamma*gamma - 1); // using inner TPC mom. as approx. 
        // }
      } else { 
        probMis = 0.; nsigmaTOFp = 0.; //cout<<"The corresponding global track has no tof pid!"<<endl;
      } 
    }     
    
    if (!fkCutOnTPCIP) {   
      // in AliAODTrack::PropagateToDCA() the lines below are iplemented 
      globaltrack->GetPosition(xyz);
      if (xyz[0]*xyz[0]+xyz[1]*xyz[1] >3.*3.) {
        //cout<<" Allowed only for tracks inside the beam pipe!!"<<endl;
         continue;
      }
      etp1.CopyFromVTrack(vtrackg);
      if (!etp1.PropagateToDCA(vertexmain,bfield,100.,dzg,covarg)) { // Does not work for TPC constrained
        //cout<<"wrong propagation!!! "<<dzg[0]<<" and  "<<dzg[1]<<endl;
        continue;
      }
    }

    // Fill histos for purity studies, after that cut on IP
    // Preselection on DCA z
    if (fkCutOnTPCIP) { 
      if (TMath::Abs(dz[1])>fIPCutz ) continue;   // 3.2  proton 1. pion
    } else {
      if (TMath::Abs(dzg[1])>fIPCutz ) continue;   // 0.15 proton/pion
    }

    // Fill here histos for purity studies only for tracks with IP within cuts
    if (lcentrality>0. && lcentrality<=10.) {
      if (fkCutOnTPCIP) {
        if (TMath::Abs(dz[0])<fIPCutxy) {
          if (trackmomentum < fMomemtumLimitForTOFPID) fHistpsignalTOFsignalTPC->Fill(trackmomentum,0.,tpcSignal);
          else {
            if (isTOFPIDok&&(probMis>0.01)) fHistpsignalTOFsignalTPC->Fill(trackmomentum,tTOF,tpcSignal);
          }
        }
      } else {
        if (TMath::Abs(dzg[0])<fIPCutxy) {
          if (trackmomentum < fMomemtumLimitForTOFPID) fHistpsignalTOFsignalTPC->Fill(trackmomentum,0.,tpcSignal);
          else {
            if (isTOFPIDok&&(probMis>0.01)) fHistpsignalTOFsignalTPC->Fill(trackmomentum,tTOF,tpcSignal);
          }
        }
      }
    }

    // Cut on PID
    if (fFirstpart == kProton) {
      if (trackmomentum< fMomemtumLimitForTOFPID) {
        if (TMath::Abs(nsigmaTPCp)> fnSigmaTPCPIDfirstParticle) continue;
      } else {
        if (isTOFPIDok) {
          if (probMis > 0.01) continue;
          if (TMath::Sqrt(nsigmaTOFp*nsigmaTOFp+nsigmaTPCp*nsigmaTPCp)> fnSigmaTPCTOFPIDfirstParticle) continue;   // this cleans the TOF corrected time plot vs p       
        } else {
          continue; // if no TOF we dont use this track 
          // if (TMath::Abs(nsigmaTPCp)> 3.) continue; 
        }
      }
    } else if (fFirstpart == kPion) { // for pions only TPC
        if (TMath::Abs(nsigmaTPCp)> fnSigmaTPCPIDfirstParticle) continue;
    }

    if (lcentrality>0. && lcentrality<=10.) {
      fHistIPtoPVxyGlobalvspt->Fill(tracktransversemomentum,dzg[0]);
      if (fReadMCTruth) { 
        mcProtonOrigin = ProtonOrigin( TMath::Abs(globaltrack->GetLabel()), arrayMC, pMomentumTruth); //  for purity and momentum resolution correction
        if (isP) {
          if (mcProtonOrigin == AliReconstructedProton::kPrimaryP) fHistMCPrimProtons->Fill(tracktransversemomentum,dzg[0]);   
          else if (mcProtonOrigin == AliReconstructedProton::kSecondaryWeak) fHistMCFromWdecayProtons->Fill(tracktransversemomentum,dzg[0]);
          else if (mcProtonOrigin == AliReconstructedProton::kMaterial) fHistMCFromMaterialProtons->Fill(tracktransversemomentum,dzg[0]);
          else if (mcProtonOrigin == AliReconstructedProton::kUnassigned) fHistMCOtherProtons->Fill(tracktransversemomentum,dzg[0]);
        } else if (isaP) {
          if (mcProtonOrigin == AliReconstructedProton::kPrimaryP) fHistMCPrimAProtons->Fill(tracktransversemomentum,dzg[0]);
          else if (mcProtonOrigin == AliReconstructedProton::kSecondaryWeak) fHistMCFromWdecayAProtons->Fill(tracktransversemomentum,dzg[0]);
          else if (mcProtonOrigin == AliReconstructedProton::kMaterial) fHistMCFromMaterialAProtons->Fill(tracktransversemomentum,dzg[0]);
          else if (mcProtonOrigin == AliReconstructedProton::kUnassigned) fHistMCOtherAProtons->Fill(tracktransversemomentum,dzg[0]);
        }
      }
    }
    if (fkCutOnTPCIP) {
      if (TMath::Abs(dz[0])>fIPCutxy ) continue;  // 2.4 proton 1. pion 
    } else {
      if (TMath::Abs(dzg[0])>fIPCutxy ) continue;  // 0.1 proton/2 pion
    }

    
    // 
    fHistyptProtons->Fill(tracktransversemomentum,rapidity,lcentrality);
    fHistphietaProtons->Fill(track->Phi(),track->Eta());
    fHistIPtoPVxyzTPC->Fill(dz[0],dz[1]); 
    fHistIPtoPVxyzGlobal->Fill(dzg[0],dzg[1]);
   
//       cout<<"TOF signal "<<globaltrack->GetTOFsignal()<<endl;
    if (lcentrality>0. && lcentrality<=10.) {
      if (isTOFPIDok) {
        fHistpTOFmisvsp->Fill(trackmomentum,probMis);
        fHistpTOFnsigmavsp->Fill(trackmomentum,nsigmaTOFp);
        fHistpTOFsignalvsp->Fill(trackmomentum,tTOF);  
        fHistpnsigTOFnsigTPC->Fill(trackmomentum,nsigmaTOFp,nsigmaTPCp);
      }

      fHistprimpTPCdEdx->Fill(globaltrack->GetTPCmomentum()*charge, tpcSignal);
    } 
   
    fHistnTPCCrossedR->Fill(nTPCCrossedRows);
    fHistRationTPCCrossedRnFind->Fill(rationCrnFind);
    fHistSharedFrTPCcl->Fill(sharedFractionTPCcl);

    if (fReadMCTruth)  {
      for (int pIndex = 0; pIndex <3; pIndex++) {
        fEvt->fReconstructedProton[pCount].pMomentumTruth[pIndex]  = pMomentumTruth[pIndex];
      }
    }

    //Save first particle information
    fEvt->fReconstructedProton[pCount].mcProtonOriginType = mcProtonOrigin;
    fEvt->fReconstructedProton[pCount].pCharge = charge;
    //cout<<"Charge of pion "<<fEvt->fReconstructedProton[pCount].pCharge<<endl;
    
    fEvt->fReconstructedProton[pCount].pMomentum[0]  = track->Px();
    fEvt->fReconstructedProton[pCount].pMomentum[1]  = track->Py();
    fEvt->fReconstructedProton[pCount].pMomentum[2]  = track->Pz();

    fEvt->fReconstructedProton[pCount].pPt     = tracktransversemomentum;
    fEvt->fReconstructedProton[pCount].pEta    = track->Eta();
    fEvt->fReconstructedProton[pCount].pPhi    = track->Phi();
    fEvt->fReconstructedProton[pCount].pRap    = rapidity; 
    fEvt->fReconstructedProton[pCount].doSkipOver = kFALSE;

    if (isP) { 
      fEvt->fReconstructedProton[pCount].isP = kTRUE; 
      fEvt->fReconstructedProton[pCount].isaP = kFALSE;
    } else if (isaP) {
      fEvt->fReconstructedProton[pCount].isaP = kTRUE; 
      fEvt->fReconstructedProton[pCount].isP = kFALSE;
    }

    fEvt->fReconstructedProton[pCount].index = TMath::Abs(globaltrack->GetID());

//    cout<<" Pos 125 cm before set "<<fEvt->fReconstructedProton[pCount].pShiftedGlobalPosition[0]<<" "<<fEvt->fReconstructedProton[pCount].pShiftedGlobalPosition[1]<<" "<<fEvt->fReconstructedProton[pCount].pShiftedGlobalPosition[2]<<endl;
    if (fkPropagateGlobal) {
      if (!fkPropagateAtFixedR) SetSftPosR125(vtrackg, bfield, lBestPrimaryVtxPos, fEvt->fReconstructedProton[pCount].pShiftedGlobalPosition); // Hans popagates the TPC tracks, I try both, flag to choose 
      else                      SetPosR125(vtrackg, bfield, fEvt->fReconstructedProton[pCount].pShiftedGlobalPosition );

    } else {
      if (!fkPropagateAtFixedR) SetSftPosR125(vtrack, bfield, lBestPrimaryVtxPos, fEvt->fReconstructedProton[pCount].pShiftedGlobalPosition);
      else                      SetPosR125(vtrack, bfield, fEvt->fReconstructedProton[pCount].pShiftedGlobalPosition );   
    }
//    cout<<" Pos set 125 cm "<<fEvt->fReconstructedProton[pCount].pShiftedGlobalPosition[0]<<" "<<fEvt->fReconstructedProton[pCount].pShiftedGlobalPosition[1]<<" "<<fEvt->fReconstructedProton[pCount].pShiftedGlobalPosition[2]<<endl; 
    fEvt->fReconstructedProton[pCount].pEtaS = EtaS(fEvt->fReconstructedProton[pCount].pShiftedGlobalPosition); 
//    cout<<" Etas "<<fEvt->fReconstructedProton[pCount].pEtaS<<endl;

//    cout<<" This proton in the stack with pt "<<fEvt->fReconstructedProton[pCount].pPt<<endl;

    pCount++;
    if (fMaxPMult <= pCount){
      cerr<<"Proton counts exceeded "<<fMaxPMult<<"!"<<endl;
      break;
    }

  }

  
  fEvt->fNumberCandidateProton = pCount;
  fHistProtonMultvsCent->Fill(pCount,lcentrality);

  if (fSecondpart == kXi || fSecondpart == kOmega) {

  // Read and plot Xi both ESDs and AODs 
    Double_t lInvMassXiMinus = 0.;
    Double_t lInvMassXiPlus = 0.;
    Double_t lInvMassOmegaMinus = 0.;
    Double_t lInvMassOmegaPlus = 0.;

    Double_t lInvMassXi = 0.;
    Double_t lInvMassLambdaAsCascDghter = 0.;
    Double_t lInvMassAntiLambdaAsCascDghter = 0.;
    Double_t lInvMassLambda = 0.;

    Int_t lPosTPCClusters = -1; 
    Int_t lNegTPCClusters = -1;  
    Int_t lBachTPCClusters = -1; 

    Int_t lPosTPCClustersS;
    Int_t lNegTPCClustersS;
    Int_t lBachTPCClustersS;

    Float_t sharedFractionTPCclPos;
    Float_t sharedFractionTPCclNeg;
    Float_t sharedFractionTPCclBach;

    AliVTrack *pTrackXi = 0x0;
    AliVTrack *nTrackXi = 0x0;
    AliVTrack *bachTrackXi = 0x0;

    Bool_t lIsBachelorPionForTPC = kFALSE;
    Bool_t lIsBachelorKaonForTPC = kFALSE;
    Bool_t lIsNegPionForTPC      = kFALSE;
    Bool_t lIsNegProtonForTPC    = kFALSE;
    Bool_t lIsPosPionForTPC      = kFALSE;
    Bool_t lIsPosProtonForTPC    = kFALSE;

    Double_t lXiMomX       = 0., lXiMomY = 0., lXiMomZ = 0.;
    Double_t lXiTransvMom  = 0., lXiMom  = 0. ;
    Double_t lV0PMomX       = 0., lV0PMomY = 0., lV0PMomZ = 0.;
    Double_t lV0NMomX       = 0., lV0NMomY = 0., lV0NMomZ = 0.;
    Double_t lV0TransvMom     = 0., lV0Mom   = 0. ; 

    Double_t lV0MomX = 0;
    Double_t lV0MomY = 0;
    Double_t lV0MomZ = 0;

    Short_t lChargeXi = -2;

    AliESDcascade *xi = 0x0;
    UInt_t lIdxPosXi  = 0;
    UInt_t lIdxNegXi  = 0;
    UInt_t lBachIdx   = 0;
    Double_t lV0quality = 0;

    ULong_t pStatus    = 0;
    ULong_t nStatus    = 0;
    ULong_t bachStatus = 0;

    const AliAODcascade *xiaod = 0x0;
    AliAODTrack* aodpTrackXi = 0x0;
    AliAODTrack* aodnTrackXi = 0x0;
    AliAODTrack* aodbachTrackXi = 0x0;

    Float_t lRapXi;  
    Float_t lEtaXi;
    Float_t lPhiXi;

    Float_t etaPos;
    Float_t etaNeg;
    Float_t etaBach;

    double xiMomentumTruth[3];
    AliReconstructedXi::MCXiOrigin_t mcXiOrigin = AliReconstructedXi::kUnassigned;
    int xiCount = 0;

    Bool_t isXi = kFALSE;
    Bool_t isaXi = kFALSE;
    Bool_t isOmega = kFALSE;
    Bool_t isaOmega = kFALSE;

    int indexB, indexP, indexN;
    int labelB, labelP, labelN;

    // Topological variables
    Double_t lContainerCutVars[22] = {0.};

    Double_t lChi2Xi = 0. ; // not implemented 
    Double_t lDcaXiDaughters = -1. ;  
    Double_t lXiCosineOfPointingAngle = -1. ; 
    Double_t lPosXi[3] = { -1000.0, -1000.0, -1000.0 }; 
    Double_t lXiRadius = -1000. ; 
    Float_t lctau = -1.;
    Float_t lctauV0 = -1.;
    Float_t distV0Xi = -1.;
    Float_t ldistTV0Xi = -1.;

    // V0 part in cascades
    Double_t lDcaV0DaughtersXi = -1.; 
    Double_t lV0RadiusXi = -1.;  
    Double_t lPosV0Xi[3] = {-1.,-1.,-1.}; 

    Double_t lDcaBachToPrimVertexXi = -1., lDcaV0ToPrimVertexXi = -1.; 
    Double_t lDcaPosToPrimVertexXi = -1.; 
    Double_t lDcaNegToPrimVertexXi = -1.; 
    Double_t lV0CosineOfPointingAngle = -1. ; 
    Double_t lV0toXiCosineOfPointingAngle = 0. ; // not used

    for (Int_t iXi = 0; iXi < ncascades; iXi++) {

      lInvMassXiMinus = 0.;
      lInvMassXiPlus = 0.;
      lInvMassOmegaMinus = 0.;
      lInvMassOmegaPlus = 0.; 
      lInvMassXi = 0.;
      lInvMassLambdaAsCascDghter = 0.;
      lInvMassAntiLambdaAsCascDghter = 0.;
      lPosTPCClusters = 0;
      lNegTPCClusters = 0;
      lBachTPCClusters = 0;

      pTrackXi = 0x0;
      nTrackXi = 0x0;
      bachTrackXi = 0x0;

      lIsBachelorPionForTPC = kFALSE;
      lIsBachelorKaonForTPC = kFALSE;
      lIsNegPionForTPC      = kFALSE;
      lIsNegProtonForTPC    = kFALSE;
      lIsPosPionForTPC      = kFALSE;
      lIsPosProtonForTPC    = kFALSE;

      lXiMomX       = 0.; lXiMomY  = 0.; lXiMomZ  = 0.;
      lXiTransvMom  = 0.;
      lV0PMomX      = 0.; lV0PMomY = 0.; lV0PMomZ = 0.;
      lV0NMomX      = 0.; lV0NMomY = 0.; lV0NMomZ = 0.;
      lV0TransvMom  = 0.; lV0Mom   = 0. ;


      lChargeXi = -2;

      isXi = kFALSE;
      isaXi = kFALSE;
      isOmega = kFALSE;
      isaOmega = kFALSE;

     
      if (fAnalysisType == "ESD") { // not mantained 

        xi = fESDevent->GetCascade(iXi);
        if (!xi) continue;

        lIdxPosXi  = (UInt_t) TMath::Abs( xi->GetPindex() );
        lIdxNegXi  = (UInt_t) TMath::Abs( xi->GetNindex() );
        lBachIdx   = (UInt_t) TMath::Abs( xi->GetBindex() );
        if (lBachIdx == lIdxNegXi) { // prevented already in cascadevertexer
          AliWarning("Pb / Idx(Bach. track) = Idx(Neg. track) ... continue!"); continue;
        }
        if (lBachIdx == lIdxPosXi) {
          AliWarning("Pb / Idx(Bach. track) = Idx(Pos. track) ... continue!"); continue;
        }

        pTrackXi             = (AliVTrack*) fESDevent->GetTrack( lIdxPosXi ); // FIXME check if conversion to VTrack works
        nTrackXi             = (AliVTrack*) fESDevent->GetTrack( lIdxNegXi );
        bachTrackXi          = (AliVTrack*) fESDevent->GetTrack( lBachIdx );
        if (!pTrackXi || !nTrackXi || !bachTrackXi ) {
          AliWarning("ERROR: Could not retrieve one of the 3 ESD daughter tracks of the cascade ...");
          continue;
        }

        lV0quality = 0;
        // Charge
        if ( bachTrackXi->Charge() < 0 )  {  
          lV0quality = 0.;
          xi->ChangeMassHypothesis(lV0quality , 3312);   // Xi- 
          lInvMassXiMinus = xi->GetEffMassXi();
          lInvMassLambdaAsCascDghter = xi->GetEffMass(); // This method without parameters returns fEffMass (set in CascadeVertexer)
        } else if ( bachTrackXi->Charge() >  0 ) {
          lV0quality = 0.;
          xi->ChangeMassHypothesis(lV0quality , -3312);  // Xi+, no need to go back to def hyp 
          lInvMassXiPlus = xi->GetEffMassXi();
          lInvMassAntiLambdaAsCascDghter = xi->GetEffMass(); // This method without parameters returns fEffMass (set in CascadeVertexer) 
        }
        lChargeXi = xi->Charge();

        pStatus    = pTrackXi->GetStatus();
        nStatus    = nTrackXi->GetStatus();
        bachStatus = bachTrackXi->GetStatus();

//        if ((pStatus&AliESDtrack::kTPCrefit)    == 0) { AliWarning("Pb / V0 Pos. track has no TPCrefit ... continue!"); continue; } // already in V0 finder
//        if ((nStatus&AliESDtrack::kTPCrefit)    == 0) { AliWarning("Pb / V0 Neg. track has no TPCrefit ... continue!"); continue; } 
        if ((bachStatus&AliESDtrack::kTPCrefit) == 0) { AliWarning("Pb / Bach.   track has no TPCrefit ... continue!"); continue; }

        xi->GetPxPyPz( lXiMomX, lXiMomY, lXiMomZ );

        xi->GetNPxPyPz(lV0NMomX,lV0NMomY,lV0NMomZ);
        xi->GetPPxPyPz(lV0PMomX,lV0PMomY,lV0PMomZ);

        lV0TransvMom = TMath::Sqrt(TMath::Power(lV0NMomX+lV0PMomX,2)+TMath::Power(lV0NMomY+lV0PMomY,2));


      } else if (fAnalysisType == "AOD") {

        xiaod = fAODevent->GetCascade(iXi);
        if (!xiaod) continue;

        pTrackXi    = dynamic_cast<AliVTrack*>( xiaod->GetDaughter(0) );
        nTrackXi    = dynamic_cast<AliVTrack*>( xiaod->GetDaughter(1) );
        bachTrackXi = dynamic_cast<AliVTrack*>( xiaod->GetDecayVertexXi()->GetDaughter(0) );
        if (!pTrackXi || !nTrackXi || !bachTrackXi ) {
          AliWarning("ERROR: Could not retrieve one of the 3 AOD daughter tracks of the cascade ...");
          continue;
        }

//      lIdxPosXi  = (UInt_t) TMath::Abs( aodpTrackXi->GetID() );  
//      lIdxNegXi  = (UInt_t) TMath::Abs( aodnTrackXi->GetID() );
//      lBachIdx   = (UInt_t) TMath::Abs( aodbachTrackXi->GetID() );

//      if (lBachIdx == lIdxNegXi) {  // Redundand, already checked in Cascade vertexer 
//        AliWarning("Pb / Idx(Bach. track) = Idx(Neg. track) ... continue!"); continue;
//      }
//      if (lBachIdx == lIdxPosXi) {
//        AliWarning("Pb / Idx(Bach. track) = Idx(Pos. track) ... continue!"); continue;
//      }

        aodpTrackXi = dynamic_cast<AliAODTrack*>(pTrackXi); // to convert aod in vtracks
        aodnTrackXi = dynamic_cast<AliAODTrack*>(nTrackXi);
        aodbachTrackXi = dynamic_cast<AliAODTrack*>(bachTrackXi);

        //if (aodpTrackXi->Charge()<0) cout<<"Wrong Sign for pos!!!!!!!!!!!!"<<endl;
        //if (aodnTrackXi->Charge()>0) cout<<"Wrong Sign for neg!!!!!!!!!!!!"<<endl;

//      if ((aodpTrackXi->IsOn(AliAODTrack::kITSpureSA)))    {  cout<<" ITSSAtrack POS!! "<<endl; }
//      if ((aodnTrackXi->IsOn(AliAODTrack::kITSpureSA)))    {  cout<<" ITSSAtrack NEG!! "<<endl; }
//      if ((aodbachTrackXi->IsOn(AliAODTrack::kITSpureSA))) {  cout<<" ITSSAtrack BACH!!"<<endl; }


//        if (!(aodpTrackXi->IsOn(AliAODTrack::kTPCrefit)))    { AliWarning(" V0 Pos. track has no TPCrefit ... continue!"); continue; } // already in V0 finder
//        if (!(aodnTrackXi->IsOn(AliAODTrack::kTPCrefit)))    { AliWarning(" V0 Neg. track has no TPCrefit ... continue!"); continue; }
        if (!(aodbachTrackXi->IsOn(AliAODTrack::kTPCrefit))) { /*AliWarning(" Bach.   track has no TPCrefit ... continue!");*/ continue; }



        lChargeXi = xiaod->ChargeXi();
        if ( lChargeXi < 0 )  { 
          lInvMassXiMinus = xiaod->MassXi(); 
          lInvMassLambdaAsCascDghter = xiaod->MassLambda(); 
          lInvMassOmegaMinus = xiaod->MassOmega();
        }
        if ( lChargeXi > 0 )  { 
          lInvMassXiPlus = xiaod->MassXi(); 
          lInvMassAntiLambdaAsCascDghter  = xiaod->MassAntiLambda();
          lInvMassOmegaPlus = xiaod->MassOmega();
        }

        indexB = bachTrackXi->GetID();
        indexP = pTrackXi->GetID();
        indexN = nTrackXi->GetID();
        // cout<<" ID bac  "<<indexB<<" n ID  "<<indexN<<" p ID "<<indexP<<endl;      

        lXiMomX = xiaod->MomXiX();
        lXiMomY = xiaod->MomXiY();
        lXiMomZ = xiaod->MomXiZ();
        lV0MomX = xiaod->MomV0X();
        lV0MomY = xiaod->MomV0Y();
        lV0MomZ = xiaod->MomV0Z();
 
        lDcaXiDaughters = xiaod->DcaXiDaughters();

        lDcaV0DaughtersXi = xiaod->DcaV0Daughters();
        lDcaV0ToPrimVertexXi = xiaod->DcaV0ToPrimVertex();

        lDcaBachToPrimVertexXi = xiaod->DcaBachToPrimVertex();
        lPosV0Xi[0] = xiaod->DecayVertexV0X();
        lPosV0Xi[1] = xiaod->DecayVertexV0Y();
        lPosV0Xi[2] = xiaod->DecayVertexV0Z();
        lV0RadiusXi = TMath::Sqrt( lPosV0Xi[0]*lPosV0Xi[0] + lPosV0Xi[1]*lPosV0Xi[1] );

        lV0CosineOfPointingAngle = xiaod->CosPointingAngle( lBestPrimaryVtxPos );

        lV0toXiCosineOfPointingAngle = xiaod->CosPointingAngle( xiaod->GetDecayVertexXi() );

        lXiCosineOfPointingAngle = xiaod->CosPointingAngleXi( lBestPrimaryVtxPos[0], lBestPrimaryVtxPos[1], lBestPrimaryVtxPos[2] );
 
        lDcaPosToPrimVertexXi = xiaod->DcaPosToPrimVertex();
        lDcaNegToPrimVertexXi = xiaod->DcaNegToPrimVertex();
 
        //cout<<" dca xi to PV "<<xiaod->DcaXiToPrimVertex()<<endl;; // gives back -999 FIXME

        lPosXi[0] = xiaod->DecayVertexXiX();
        lPosXi[1] = xiaod->DecayVertexXiY();
        lPosXi[2] = xiaod->DecayVertexXiZ();

        if (fSecondpart == kXi) lRapXi = xiaod->RapXi();
        else if (fSecondpart == kOmega) lRapXi = xiaod->RapOmega();
      
        // those are for the V0 
        //lEta      = xiaod->Eta();                             
        //lTheta    = xiaod->Theta() *180.0/TMath::Pi();       
        //lPhi      = xiaod->Phi()   *180.0/TMath::Pi();       

        lPosTPCClustersS   = aodpTrackXi->GetTPCnclsS();
        lNegTPCClustersS   = aodnTrackXi->GetTPCnclsS();
        lBachTPCClustersS  = aodbachTrackXi->GetTPCnclsS();

        etaPos = aodpTrackXi->Eta();
        etaNeg = aodnTrackXi->Eta();
        etaBach = aodbachTrackXi->Eta();

      }

      lPosTPCClusters   = pTrackXi->GetTPCNcls();
      lNegTPCClusters   = nTrackXi->GetTPCNcls();
      lBachTPCClusters  = bachTrackXi->GetTPCNcls();

      if (lPosTPCClusters  < 70) { //AliWarning("Pb / V0 Pos. track has less than minn TPC clusters ... continue!");
        continue;
      }
      if (lNegTPCClusters  < 70) { //AliWarning("Pb / V0 Neg. track has less than minn TPC clusters ... continue!");
        continue;
      }
      if (lBachTPCClusters < 70) { //AliWarning("Pb / Bach.   track has less than minn TPC clusters ... continue!");
        continue;
      }

//    if (lXiTransvMom<0.8||lXiTransvMom>8.) continue;
//    if (TMath::Abs(lRapXi)>0.5) continue;

      if (TMath::Abs(etaBach)>0.8) continue;
      if (TMath::Abs(etaPos)>0.8) continue;
      if (TMath::Abs(etaNeg)>0.8) continue;
      if (fkApplyYcutCasc) if (TMath::Abs(lRapXi)>0.5) continue; // basically contained in the previous cut

      sharedFractionTPCclPos  = (Float_t) lPosTPCClustersS/lPosTPCClusters;
      sharedFractionTPCclNeg  = (Float_t) lNegTPCClustersS/lNegTPCClusters;
      sharedFractionTPCclBach = (Float_t) lBachTPCClustersS/lBachTPCClusters;
//      cout<<" shared pos "<<sharedFractionTPCclPos<<" shared neg "<<sharedFractionTPCclNeg<<" shared bach "<<sharedFractionTPCclBach<<endl;
      if (sharedFractionTPCclPos>0||sharedFractionTPCclNeg>0||sharedFractionTPCclBach>0) continue; // around 5 % for all tracks


      lXiRadius = TMath::Sqrt( lPosXi[0]*lPosXi[0] + lPosXi[1]*lPosXi[1] );

      lV0TransvMom = TMath::Sqrt(TMath::Power(lV0MomX,2)+TMath::Power(lV0MomY,2));
      lV0Mom = TMath::Sqrt(TMath::Power(lV0MomX,2)+TMath::Power(lV0MomY,2)+TMath::Power(lV0MomZ,2));
      lV0RadiusXi = TMath::Sqrt( lPosV0Xi[0]*lPosV0Xi[0] + lPosV0Xi[1]*lPosV0Xi[1] );

      lXiTransvMom = TMath::Sqrt( lXiMomX*lXiMomX   + lXiMomY*lXiMomY );
      lXiMom = TMath::Sqrt( lXiMomX*lXiMomX   + lXiMomY*lXiMomY + lXiMomZ*lXiMomZ);

      lEtaXi = 0.5*TMath::Log((lXiMom+lXiMomZ)/(lXiMom-lXiMomZ+1.e-13));
      lPhiXi = TMath::Pi()+TMath::ATan2(-lXiMomY,-lXiMomX); // from AliAODRecoDecay.h

      lctau = -1.;
      lctau = TMath::Sqrt(TMath::Power((lPosXi[0]-lBestPrimaryVtxPos[0]),2)+TMath::Power((lPosXi[1]-lBestPrimaryVtxPos[1]),2)+TMath::Power(( lPosXi[2]-lBestPrimaryVtxPos[2]),2));
      if (fSecondpart == kXi) {
        if (lXiMom!=0) lctau = lctau*fPDGXi/lXiMom;
        else lctau = -1.;
      } else if (fSecondpart == kOmega ) {
        if (lXiMom!=0) lctau = lctau*fPDGOmega/lXiMom;
        else lctau = -1.;
      }

      distV0Xi = TMath::Sqrt(TMath::Power((lPosV0Xi[0]-lPosXi[0]),2)+TMath::Power((lPosV0Xi[1]-lPosXi[1]),2)+TMath::Power((lPosV0Xi[2]-lPosXi[2]),2));
      lctauV0 = -1.;
      if (lV0Mom!=0) lctauV0 = distV0Xi*fPDGL/lV0Mom;
      ldistTV0Xi = TMath::Sqrt(TMath::Power((lPosV0Xi[0]-lPosXi[0]),2)+TMath::Power((lPosV0Xi[1]-lPosXi[1]),2));


      // Bachelor
      if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC( bachTrackXi,AliPID::kKaon)) < fnSigmaTPCPIDsecondParticleDau) lIsBachelorKaonForTPC = kTRUE;
      if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC( bachTrackXi,AliPID::kPion)) < fnSigmaTPCPIDsecondParticleDau) lIsBachelorPionForTPC = kTRUE;
      // Negative V0 daughter
      if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC( nTrackXi,AliPID::kPion   )) < fnSigmaTPCPIDsecondParticleDau) lIsNegPionForTPC   = kTRUE;
      if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC( nTrackXi,AliPID::kProton )) < fnSigmaTPCPIDsecondParticleDau) lIsNegProtonForTPC = kTRUE;
      // Positive V0 daughter
      if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC( pTrackXi,AliPID::kPion   )) < fnSigmaTPCPIDsecondParticleDau) lIsPosPionForTPC   = kTRUE;
      if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC( pTrackXi,AliPID::kProton )) < fnSigmaTPCPIDsecondParticleDau) lIsPosProtonForTPC = kTRUE;


      if (fSecondpart == kXi) {

        if (lChargeXi < 0 && lIsBachelorPionForTPC && lIsNegPionForTPC && lIsPosProtonForTPC) {

          isXi = kTRUE;
          lInvMassLambda = lInvMassLambdaAsCascDghter;

          lInvMassXi = lInvMassXiMinus;

        } else if (lChargeXi > 0 && lIsBachelorPionForTPC && lIsPosPionForTPC && lIsNegProtonForTPC) {

          isaXi = kTRUE;
          lInvMassLambda = lInvMassAntiLambdaAsCascDghter;

          lInvMassXi = lInvMassXiPlus;

        } else  continue;

      } else if (fSecondpart == kOmega) {

        if (lChargeXi < 0 && lIsBachelorKaonForTPC && lIsNegPionForTPC && lIsPosProtonForTPC)  {
          isXi = kTRUE; // this flag is used to set the reconstructed second paticle collection 
          isOmega = kTRUE;
          lInvMassLambda = lInvMassLambdaAsCascDghter;

          lInvMassXi = lInvMassOmegaMinus;

        } else if (lChargeXi > 0 && lIsBachelorKaonForTPC && lIsPosPionForTPC && lIsNegProtonForTPC) {

          isaXi = kTRUE;
          isaOmega = kTRUE;
          lInvMassLambda = lInvMassAntiLambdaAsCascDghter;

          lInvMassXi = lInvMassOmegaPlus;

        } else continue;
      }

      // Select on topological cuts
      if (lDcaBachToPrimVertexXi<fIPBac) continue;  
      if (lDcaXiDaughters>fDcaXiDaughters) continue;
      if (lXiCosineOfPointingAngle<fXiCosineOfPointingAngle) continue;
      if (lXiRadius<fXiRadiusMin || lXiRadius>fXiRadiusMax) continue;
      if (lctau>fCtau) continue;
      if (fSecondpart == kOmega) {
        if (TMath::Abs(xiaod->MassXi()-fPDGXi)<0.008) continue; // reject omega candiates around the xi mass value
      }
      if (!fkCascadeSideBands) if ((lInvMassLambda<fPDGL-fMassWindowL)||(lInvMassLambda>fPDGL+fMassWindowL)) continue; // if I want side bands and tight cuts I should not apply the cut        
      if (lV0CosineOfPointingAngle<fV0CosineOfPointingAngle) continue;//lV0toXiCosineOfPointingAngle;
      if (lV0RadiusXi<fV0RadiusXiMin||lV0RadiusXi>fV0RadiusXiMax) continue;
      if (lDcaV0ToPrimVertexXi<fIPV0ToPrimVertexXi) continue;
      if (lDcaPosToPrimVertexXi<fIPPosToPrimVertexXi) continue;
      if (lDcaNegToPrimVertexXi<fIPNegToPrimVertexXi) continue;
      if (lDcaV0DaughtersXi>fDcaV0Daughters) continue;
      // end of selections
 
      // Fill here inv mass histos for purity estimation before cutting on onvariant mass. 
        // NB Some Xi will be removed in the pair loop but still included here
      if (lcentrality>0. && lcentrality<=10.) {
        if (isXi)           fHistInvMassXiMinus->Fill(lInvMassXiMinus, lXiTransvMom);
        else if (isaXi)     fHistInvMassXiPlus->Fill(lInvMassXiPlus, lXiTransvMom);
      }
      if (!fkCascadeSideBands) { // Peak
        if (TMath::Abs(lInvMassXi-fPDGsecond) > fMassWindowCascades ) continue; // save only xis in the selected mass range
      } else { // off-peak: fMassWindowCascades is in this case larger than the peak width to avoid the tails of the peak 
        if ((TMath::Abs(lInvMassXi-fPDGsecond) < fMassWindowCascades)|| (lInvMassXi >= fPDGsecond+4*fMassWindowCascades) || (lInvMassXi <= fPDGsecond-4*fMassWindowCascades) )
          continue;
        if (TMath::Abs(lInvMassLambda-fPDGL)<0.005) continue;
      }
      if (lXiTransvMom<fMinPtForCasc||lXiTransvMom>fMaxPtForCasc) continue;
      if (isXi) {

        fHistInvMassL->Fill(lInvMassLambdaAsCascDghter, lV0TransvMom);
        fHistpTPCdEdx->Fill(pTrackXi->GetTPCmomentum()*pTrackXi->Charge(), pTrackXi->GetTPCsignal());
        fHistnTPCdEdx->Fill(nTrackXi->GetTPCmomentum()*nTrackXi->Charge(), nTrackXi->GetTPCsignal());
        fHistbTPCdEdx->Fill(bachTrackXi->GetTPCmomentum()*bachTrackXi->Charge(), bachTrackXi->GetTPCsignal());


      } else if (isaXi) {

        fHistInvMassAntiL->Fill(lInvMassAntiLambdaAsCascDghter, lV0TransvMom);
        fHistpTPCdEdx->Fill(pTrackXi->GetTPCmomentum()*pTrackXi->Charge(), pTrackXi->GetTPCsignal());
        fHistnTPCdEdx->Fill(nTrackXi->GetTPCmomentum()*nTrackXi->Charge(), nTrackXi->GetTPCsignal());
        fHistbTPCdEdx->Fill(bachTrackXi->GetTPCmomentum()*bachTrackXi->Charge(), bachTrackXi->GetTPCsignal());

      } 
 

      fHistPosV0TPCClusters->Fill(lPosTPCClusters);  
      fHistNegV0TPCClusters->Fill(lNegTPCClusters);
      fHistBachTPCClusters->Fill(lBachTPCClusters); 

      fHistSharedFrTPCclPos->Fill(sharedFractionTPCclPos);
      fHistSharedFrTPCclNeg->Fill(sharedFractionTPCclNeg);
      fHistSharedFrTPCclBach->Fill(sharedFractionTPCclBach);

      fHistyptXi->Fill(lXiTransvMom,lRapXi); 
      fHistphietaXi->Fill(lPhiXi,lEtaXi);
      fHistptL->Fill(lV0TransvMom);
      fHistptpL->Fill(pTrackXi->Pt());
      fHistptnL->Fill(nTrackXi->Pt());
      fHistptbac->Fill(bachTrackXi->Pt());

      // Fill the container for topological cuts optimization
/*      if (fUseContainer) {

        lContainerCutVars[0] = lDcaXiDaughters;
        lContainerCutVars[1] = lDcaBachToPrimVertexXi;
        lContainerCutVars[2] = lXiCosineOfPointingAngle;
        lContainerCutVars[3] = lXiRadius;
        lContainerCutVars[4] = lInvMassLambda;
        lContainerCutVars[5] = lDcaV0DaughtersXi;
        lContainerCutVars[6] = lV0CosineOfPointingAngle;//lV0toXiCosineOfPointingAngle;
        lContainerCutVars[7] = lV0RadiusXi;
        lContainerCutVars[8] = lDcaV0ToPrimVertexXi;
        lContainerCutVars[9] = lDcaPosToPrimVertexXi;
        lContainerCutVars[10] = lDcaNegToPrimVertexXi;

        lContainerCutVars[13] = lXiTransvMom;
        lContainerCutVars[16] = lBestPrimaryVtxPos[2];
        lContainerCutVars[17] = lcentrality;
        lContainerCutVars[18] = 0.;
        // lContainerCutVars[19] = lBachTPCClusters;
        lContainerCutVars[19] = lctau;
        lContainerCutVars[20] = lctauV0;
        lContainerCutVars[21] = ldistTV0Xi;
        if ( lChargeXi < 0 ) {
          lContainerCutVars[11] = lInvMassXiMinus;
          lContainerCutVars[12] = lInvMassOmegaMinus;//1.63;
          lContainerCutVars[14] = lRapXi;
          lContainerCutVars[15] = lChi2Xi;
          if ( lIsBachelorPionForTPC && lIsPosProtonForTPC && lIsNegPionForTPC ) 
            fCFContCascadeCuts->Fill(lContainerCutVars,0); // for Xi-
      
        } else {
          lContainerCutVars[11] = lInvMassXiPlus;
          lContainerCutVars[12] = lInvMassOmegaPlus;//1.63;
          lContainerCutVars[14] = lRapXi;
          lContainerCutVars[15] = lChi2Xi;
          if ( lIsBachelorPionForTPC && lIsNegProtonForTPC && lIsPosPionForTPC )
            fCFContCascadeCuts->Fill(lContainerCutVars,1); // for Xi+

        }
      }
*/
    
      fEvt->fReconstructedXi[xiCount].xiMomentum[0]  = lXiMomX;
      fEvt->fReconstructedXi[xiCount].xiMomentum[1]  = lXiMomY;
      fEvt->fReconstructedXi[xiCount].xiMomentum[2]  = lXiMomZ;

      if (fReadMCTruth)  {
        labelB = TMath::Abs(bachTrackXi->GetLabel()); labelP = TMath::Abs(pTrackXi->GetLabel()); labelN = TMath::Abs(nTrackXi->GetLabel());
        mcXiOrigin = XiOrigin(labelB, labelP, labelN, arrayMC, xiMomentumTruth); //  for purity and momentum resolution correction   
        for (int xiIndex = 0; xiIndex <3; xiIndex++) {
          fEvt->fReconstructedXi[xiCount].xiMomentumTruth[xiIndex]  = xiMomentumTruth[xiIndex];
        }
      }
      fEvt->fReconstructedXi[xiCount].mcXiOriginType = mcXiOrigin;
      fEvt->fReconstructedXi[xiCount].xiPt     = lXiTransvMom;
      fEvt->fReconstructedXi[xiCount].xiEta    = lEtaXi;
      fEvt->fReconstructedXi[xiCount].xiPhi    = lPhiXi;
      fEvt->fReconstructedXi[xiCount].xiRap    = lRapXi;
      fEvt->fReconstructedXi[xiCount].xiMass   = lInvMassXi;


      //fEvt->fReconstructedXi[xiCount].xiDCA   = ->DcaV0ToPrimVertex();
      fEvt->fReconstructedXi[xiCount].isXi = isXi; 
      fEvt->fReconstructedXi[xiCount].isaXi = isaXi; 
      fEvt->fReconstructedXi[xiCount].indexB = TMath::Abs(indexB);
      fEvt->fReconstructedXi[xiCount].indexP = TMath::Abs(indexP);
      fEvt->fReconstructedXi[xiCount].indexN = TMath::Abs(indexN);
      fEvt->fReconstructedXi[xiCount].doSkipOver = kFALSE;
      fEvt->fReconstructedXi[xiCount].doPickOne  = kFALSE; 

      // Store for each track the position at 1.2 radius
      if (!fkPropagateAtFixedR) {
        SetSftPosR125( bachTrackXi, bfield, lBestPrimaryVtxPos, fEvt->fReconstructedXi[xiCount].daughterBacShiftedGlobalPosition);
        SetSftPosR125( pTrackXi, bfield, lBestPrimaryVtxPos, fEvt->fReconstructedXi[xiCount].daughterPosShiftedGlobalPosition);
        SetSftPosR125( nTrackXi, bfield, lBestPrimaryVtxPos, fEvt->fReconstructedXi[xiCount].daughterNegShiftedGlobalPosition);
      } else {
        SetPosR125( bachTrackXi, bfield, fEvt->fReconstructedXi[xiCount].daughterBacShiftedGlobalPosition);
        SetPosR125( pTrackXi, bfield, fEvt->fReconstructedXi[xiCount].daughterPosShiftedGlobalPosition);
        SetPosR125( nTrackXi, bfield, fEvt->fReconstructedXi[xiCount].daughterNegShiftedGlobalPosition);
      }
      fEvt->fReconstructedXi[xiCount].daughterBacEtaS = EtaS(fEvt->fReconstructedXi[xiCount].daughterBacShiftedGlobalPosition); 
      fEvt->fReconstructedXi[xiCount].daughterPosEtaS = EtaS(fEvt->fReconstructedXi[xiCount].daughterPosShiftedGlobalPosition);
      fEvt->fReconstructedXi[xiCount].daughterNegEtaS = EtaS(fEvt->fReconstructedXi[xiCount].daughterNegShiftedGlobalPosition);

//    cout<<"Array is xi or axi "<<fEvt->fReconstructedXi[xiCount].isXi<<" "<<fEvt->fReconstructedXi[xiCount].isaXi<<" from flag s xi or axi "<<isXi<<"  "<<isaXi<<endl;
//    cout<<"PDGXi mass "<<fPDGXi<<endl; 
//    cout<<"mass "<<lInvMassXi<<endl;
      fEvt->fReconstructedXi[xiCount].isXiMass = kTRUE;  // selection done earlier 
//    else fEvt->fReconstructedXi[xiCount].isXiMass = kFALSE; // modified to reduce the size of the array

//    cout<<"This Xi in the stack with pT "<<fEvt->fReconstructedXi[xiCount].xiPt<<endl; 
//    cout<<"and the pt of the xi was "<<lXiTransvMom<<endl;
      xiCount++;
      if (fMaxXiMult <= xiCount) {
        cerr<<"Xi counts exceeded "<<fMaxXiMult<<"!"<<endl;
        break;
      }
    } 

    fEvt->fNumberCandidateXi = xiCount;
    //cout<<"Xi mult "<<xiCount<<endl; 
    fHistXiMultvsCent->Fill(xiCount,lcentrality);

    if (fkApplySharedDaughterCutXi) {
      // Check if Xi share daughters!!
      Int_t nXiWithSharedDaughtersWithXiMass = 0.;
      Int_t nXiWithXiMass = xiCount; //0.; 
      Int_t daughterIndex = -1;
      Int_t ncheckeddaughters = 0;  
      Int_t checkeddaughters[(fEvt->fNumberCandidateXi)*3];

      for (int i=0; i < fEvt->fNumberCandidateXi*3; i++) checkeddaughters[i] = -99999;

      for (int i=0; i < fEvt->fNumberCandidateXi-1; i++) {
        nXiWithSharedDaughtersWithXiMass = 0;

        if (fEvt->fReconstructedXi[i].isXiMass == kTRUE && !(fEvt->fReconstructedXi[i].doSkipOver) ) {
          //cout<<"Checking Xi: "<<i<<endl;
          //cout<<"Checking bachelor  "<<endl;
          daughterIndex = fEvt->fReconstructedXi[i].indexB;
          nXiWithSharedDaughtersWithXiMass = CheckDaughterTrack ( i, daughterIndex, checkeddaughters);
          if (nXiWithSharedDaughtersWithXiMass!=-1) {
            fHistFractionOfXiWithSharedDaughters->Fill(nXiWithXiMass, nXiWithSharedDaughtersWithXiMass/(Float_t)  nXiWithXiMass);
            //if (nXiWithSharedDaughtersWithXiMass>nXiWithXiMass) cout<<"Problem xitot<xishared anal bac "<<nXiWithXiMass<<" "<<nXiWithSharedDaughtersWithXiMass<<endl;
            //if (nXiWithSharedDaughtersWithXiMass>0) cout<<"Xi sharing this daughter "<<nXiWithSharedDaughtersWithXiMass<<" out of tot Xi with mass "<<nXiWithXiMass<<endl;

            checkeddaughters[ncheckeddaughters] = daughterIndex; 
            ncheckeddaughters++;
          }
        
          //cout<<"Checking pos daughter "<<endl;
          daughterIndex = fEvt->fReconstructedXi[i].indexP;
          nXiWithSharedDaughtersWithXiMass = CheckDaughterTrack ( i, daughterIndex, checkeddaughters);

          if (nXiWithSharedDaughtersWithXiMass!=-1) {
            fHistFractionOfXiWithSharedDaughters->Fill(nXiWithXiMass, nXiWithSharedDaughtersWithXiMass/(Float_t)  nXiWithXiMass);
            //if (nXiWithSharedDaughtersWithXiMass>nXiWithXiMass) cout<<"Problem xitot<xishared anal p "<<nXiWithXiMass<<" "<<nXiWithSharedDaughtersWithXiMass<<endl;
            //if (nXiWithSharedDaughtersWithXiMass>0) cout<<"Xi sharing this daughter "<<nXiWithSharedDaughtersWithXiMass<<" out of tot Xi with mass "<<nXiWithXiMass<<endl;

            checkeddaughters[ncheckeddaughters] = daughterIndex;
            ncheckeddaughters++;
          }

          //cout<<"Checking neg daughter  "<<endl;
          daughterIndex = fEvt->fReconstructedXi[i].indexN;
          nXiWithSharedDaughtersWithXiMass = CheckDaughterTrack ( i, daughterIndex, checkeddaughters);
          if (nXiWithSharedDaughtersWithXiMass!=-1) {
            fHistFractionOfXiWithSharedDaughters->Fill(nXiWithXiMass, nXiWithSharedDaughtersWithXiMass/(Float_t)  nXiWithXiMass);
            //if (nXiWithSharedDaughtersWithXiMass>nXiWithXiMass) cout<<"Problem xitot<xishared anal n "<<nXiWithXiMass<<" "<<nXiWithSharedDaughtersWithXiMass<<endl;
            //if (nXiWithSharedDaughtersWithXiMass>0) cout<<"Xi sharing this daughter "<<nXiWithSharedDaughtersWithXiMass<<" out of tot Xi with mass "<<nXiWithXiMass<<endl;
            checkeddaughters[ncheckeddaughters] = daughterIndex;
            ncheckeddaughters++;
          }
        }   
      }
      // cout<<"number of daughters actually checked "<<ncheckeddaughters<<endl;  
      // cout<<"Fraction of xis with shared daughters "<< nXiWithSharedDaughtersWithXiMass/ (Float_t) nXiWithXiMass<<endl;
      // cout<<"Fraction of antixis with shared daughters "<< naXiWithSharedDaughtersWithXiMass/(Float_t)  naXiWithXiMass<<endl;
      nXiWithSharedDaughtersWithXiMass = 0;
      for (int i=0; i < fEvt->fNumberCandidateXi; i++) { // never in this loop if nXi in inv mass window is 0
        if ((fEvt->fReconstructedXi[i].isXiMass == kTRUE) && (fEvt->fReconstructedXi[i].doSkipOver == kTRUE)) nXiWithSharedDaughtersWithXiMass++ ;
      }
      fHistTotMaxFractionOfXiWithSharedDaughters->Fill(nXiWithXiMass, nXiWithSharedDaughtersWithXiMass/(Float_t)  nXiWithXiMass);
    } 
 
      // Remove protons anx Xi when the proton is one of the daughter of the Xi, 
      // NB done after shared daughter flagging to avoid unflagging if the iterative procedure for shared daughter is implemented
      //cout<<"Proton mult "<<fEvt->fNumberCandidateProton<<endl;
    for (int i=0; i < fEvt->fNumberCandidateXi; i++) {
      for (int j=0; j<fEvt->fNumberCandidateProton; j++) {
        if (fEvt->fReconstructedXi[i].indexB == fEvt->fReconstructedProton[j].index ||
            fEvt->fReconstructedXi[i].indexP == fEvt->fReconstructedProton[j].index ||
            fEvt->fReconstructedXi[i].indexN == fEvt->fReconstructedProton[j].index   ) {
            //cout<<"Proton is also daughter of Xi!! Skipping proton and Xi!"<<endl;
          fEvt->fReconstructedXi[i].doSkipOver = kTRUE;
          fEvt->fReconstructedProton[j].doSkipOver = kTRUE;

        }
      }
    }


    DoPairshCasc(lcentrality); 
  } else DoPairshh(lcentrality, fieldsign);  

// Post output data
  PostData(1, fOutputContainer);
//  PostData(2, fCFContCascadeCuts);

}

//_______________________________________________________________

void AliAnalysisTaskhCascadeFemto::DoPairshCasc (const Float_t lcentrality) {

 
  bool ismassXi;
  bool isXi;
  bool isaXi;
  bool isP;
  bool isaP;
  double dphisbac;
  double detasbac;
  double dphispos;
  double detaspos;
  double dphisneg;
  double detasneg;

//  cout<<"Number of xi candates in this event "<< event->fNumberCandidateXi<<endl;


  int evmultmixed = 0;
  bool multmixedcounted = kFALSE;
  double pairKstar = 0.;
  double pairKstargen = 0.;

  for (int i=0; i < fEvt->fNumberCandidateXi; i++) { //Start looping over reconstructed Xis in this event 

   if (fReadMCTruth) {
     if (fSecondpart == kXi) if ((fEvt->fReconstructedXi[i].mcXiOriginType != AliReconstructedXi::kPrimaryXi) && (fEvt->fReconstructedXi[i].mcXiOriginType!=AliReconstructedXi::kPrimaryAntiXi)) continue;
     else if (fSecondpart == kOmega) if ((fEvt->fReconstructedXi[i].mcXiOriginType != AliReconstructedXi::kPrimaryOmega) && (fEvt->fReconstructedXi[i].mcXiOriginType!=AliReconstructedXi::kPrimaryAntiOmega)) continue;
    }

    ismassXi  = fEvt->fReconstructedXi[i].isXiMass;  
    

    if (!ismassXi) continue;
    if (fEvt->fReconstructedXi[i].doSkipOver) continue;
//    cout<<"Found Xi with correct mass"<<endl;
    isXi = fEvt->fReconstructedXi[i].isXi;
    isaXi = fEvt->fReconstructedXi[i].isaXi;    

    fHistInvMassXiInPairs->Fill(fEvt->fReconstructedXi[i].xiMass,lcentrality);

    for (int eventNumber=0; eventNumber<fnEventsToMix+1; eventNumber++) { // Event buffer loop: eventNumber=0 is the current event, all other eventNumbers are past events 
                                                                          

      // For same event pairs

      if (!multmixedcounted && eventNumber!=0 && ((fEvt+eventNumber)->fNumberCandidateProton)!=0.) evmultmixed++; 

      for (int j=0; j<(fEvt+eventNumber)->fNumberCandidateProton; j++) { // Second particle loop (from past or current event)

        if (fReadMCTruth)
          if ((fEvt+eventNumber)->fReconstructedProton[j].mcProtonOriginType != AliReconstructedProton::kPrimaryP) continue;

        if ((fEvt+eventNumber)->fReconstructedProton[j].doSkipOver) continue; // remove also from mixing, if they were used in Xi we are not sure about the quality  
        isP = (fEvt+eventNumber)->fReconstructedProton[j].isP;    
        isaP = (fEvt+eventNumber)->fReconstructedProton[j].isaP;  

        //Calculate k* for the pair
        pairKstar = CalculateKstar(fEvt->fReconstructedXi[i].xiMomentum, (fEvt+eventNumber)->fReconstructedProton[j].pMomentum, fPDGsecond,fPDGfirst);
        if (fReadMCTruth) pairKstargen = CalculateKstar(fEvt->fReconstructedXi[i].xiMomentumTruth, (fEvt+eventNumber)->fReconstructedProton[j].pMomentumTruth, fPDGsecond,fPDGfirst);          
        // Pair histogramming
        dphisbac = CalculateDphiSatR12m( fEvt->fReconstructedXi[i].daughterBacShiftedGlobalPosition, (fEvt+eventNumber)->fReconstructedProton[j].pShiftedGlobalPosition );
        detasbac = fEvt->fReconstructedXi[i].daughterBacEtaS - (fEvt+eventNumber)->fReconstructedProton[j].pEtaS ;

        dphispos = CalculateDphiSatR12m( fEvt->fReconstructedXi[i].daughterPosShiftedGlobalPosition, (fEvt+eventNumber)->fReconstructedProton[j].pShiftedGlobalPosition );
        detaspos = fEvt->fReconstructedXi[i].daughterPosEtaS - (fEvt+eventNumber)->fReconstructedProton[j].pEtaS ;

        dphisneg = CalculateDphiSatR12m( fEvt->fReconstructedXi[i].daughterNegShiftedGlobalPosition, (fEvt+eventNumber)->fReconstructedProton[j].pShiftedGlobalPosition );
        detasneg = fEvt->fReconstructedXi[i].daughterNegEtaS - (fEvt+eventNumber)->fReconstructedProton[j].pEtaS ;

  //      if (fEvt->fReconstructedXi[i].daughterBacShiftedGlobalPosition[0]==-9999. || fEvt->fReconstructedXi[i].daughterPosShiftedGlobalPosition[0]==-9999. || fEvt->fReconstructedXi[i].daughterNegShiftedGlobalPosition[0]==-9999.) {
//          cout<<" pion coord "<<(fEvt+eventNumber)->fReconstructedProton[j].pShiftedGlobalPosition[0]<<" xibac "<<fEvt->fReconstructedXi[i].daughterBacShiftedGlobalPosition[0]<<" vo dau "<<fEvt->fReconstructedXi[i].daughterPosShiftedGlobalPosition[0]<<" v0 dau "<<fEvt->fReconstructedXi[i].daughterNegShiftedGlobalPosition[0]<<endl; 

//          cout<<" dphibac "<<dphisbac<<" detasbac "<<detasbac<<" dphispos "<<dphispos<<" detaspos "<<detaspos<<" dphisneg "<<dphisneg<<" detasneg "<<detasneg<<endl;

  //     }
        // Apply two-track cuts
        if (fkApplyTtc) {  // dont cut cases in which propagation failed for both tracks and all coordinates are set to 9999, there dphi deta are 0
                           // happens seldom, mainly for primary pions of very low momentum (0.14 GeV/c)
          if (isP&&isXi) {
            if (fFirstpart==kPion && (fSecondpart==kXi || fSecondpart==kOmega)) {
              if (dphisneg!=0. && detasneg!=0. && TMath::Abs(dphisneg)<fDphisMin && TMath::Abs(detasneg)<fDetasMin) continue;
              if (dphisbac!=0. && detasbac!=0. && TMath::Abs(dphisbac)<fDphisMin && TMath::Abs(detasbac)<fDetasMin) continue;
            } else if (fFirstpart==kProton && (fSecondpart==kXi || fSecondpart==kOmega)) {
              if (dphispos!=0. && detaspos!=0. && TMath::Abs(dphispos)<fDphisMin && TMath::Abs(detaspos)<fDetasMin) continue;
            }
          } else if (isaP&&isaXi) {
            if (fFirstpart==kPion && (fSecondpart==kXi || fSecondpart==kOmega)) {
              if (dphispos!=0. && detaspos!=0. && TMath::Abs(dphispos)<fDphisMin && TMath::Abs(detaspos)<fDetasMin) continue;
              if (dphisbac!=0. && detasbac!=0. && TMath::Abs(dphisbac)<fDphisMin && TMath::Abs(detasbac)<fDetasMin) continue;
            } else if (fFirstpart==kProton && (fSecondpart==kXi || fSecondpart==kOmega)) {
              if (dphisneg!=0. && detasneg!=0. && TMath::Abs(dphisneg)<fDphisMin && TMath::Abs(detasneg)<fDetasMin) continue;
            }
          } else if (isaP&&isXi) {
     
            if (fFirstpart==kPion && (fSecondpart==kXi || fSecondpart==kOmega)) {
              if (dphispos!=0. && detaspos!=0. && TMath::Abs(dphispos)<fDphisMin && TMath::Abs(detaspos)<fDetasMin) continue;
            } else if (fFirstpart==kProton && (fSecondpart==kXi || fSecondpart==kOmega)) {
              if (dphisbac!=0. && detasbac!=0. && TMath::Abs(dphisbac)<fDphisMin && TMath::Abs(detasbac)<fDetasMin) continue;
              if (dphisneg!=0. && detasneg!=0. && TMath::Abs(dphisneg)<fDphisMin && TMath::Abs(detasneg)<fDetasMin) continue;
            }

          } else if (isP&&isaXi) {
            if (fFirstpart==kPion && (fSecondpart==kXi || fSecondpart==kOmega)) {
              if (dphisneg!=0. && detasneg!=0. && TMath::Abs(dphisneg)<fDphisMin && TMath::Abs(detasneg)<fDetasMin) continue;
            } else if (fFirstpart==kProton && (fSecondpart==kXi || fSecondpart==kOmega)) {
              if (dphisbac!=0. && detasbac!=0. && TMath::Abs(dphisbac)<fDphisMin && TMath::Abs(detasbac)<fDetasMin) continue;
              if (dphispos!=0. && detaspos!=0. && TMath::Abs(dphispos)<fDphisMin && TMath::Abs(detaspos)<fDetasMin) continue;
            }

          }
           
        }



        if (eventNumber==0) {//Same event pair histogramming
          if (isXi&&isP) {
  //          if (TMath::Abs(dphispos)>0.015&&TMath::Abs(detaspos)>0.1) {
            fHistpXiSignalRealKstar->Fill(pairKstar,lcentrality);//}
//            if (pairKstar<0.1) {
            if (dphisbac!=0.&& detasbac!=0.) fHistpXibacDetaSDphiS->Fill(dphisbac,detasbac); 
            if (dphispos!=0.&& detaspos!=0.) fHistpXiposDetaSDphiS->Fill(dphispos,detaspos);
            if (dphisneg!=0.&& detasneg!=0.) fHistpXinegDetaSDphiS->Fill(dphisneg,detasneg); 
//            }
          } else if (isXi&&isaP) {
            fHistapXiSignalRealKstar->Fill(pairKstar,lcentrality);
//            if (pairKstar<0.1) {
            if (dphisbac!=0.&& detasbac!=0.) fHistapXibacDetaSDphiS->Fill(dphisbac,detasbac);
            if (dphispos!=0.&& detaspos!=0.) fHistapXiposDetaSDphiS->Fill(dphispos,detaspos);
            if (dphisneg!=0.&& detasneg!=0.) fHistapXinegDetaSDphiS->Fill(dphisneg,detasneg);
//            }
          } else if (isaXi&&isP) {
            fHistpaXiSignalRealKstar->Fill(pairKstar,lcentrality);
//            if (pairKstar<0.1) {
            if (dphisbac!=0.&& detasbac!=0.) fHistpaXibacDetaSDphiS->Fill(dphisbac,detasbac); 
            if (dphispos!=0.&& detaspos!=0.) fHistpaXiposDetaSDphiS->Fill(dphispos,detaspos); 
            if (dphisneg!=0.&& detasneg!=0.) fHistpaXinegDetaSDphiS->Fill(dphisneg,detasneg); 
//            }
          } else if (isaXi&&isaP) {
         //   if (TMath::Abs(dphisneg)>0.015&&TMath::Abs(detasneg)>0.1) { 
            fHistapaXiSignalRealKstar->Fill(pairKstar,lcentrality);//}
//            if (pairKstar<0.1) {
            if (dphisbac!=0.&& detasbac!=0.) fHistapaXibacDetaSDphiS->Fill(dphisbac,detasbac); 
            if (dphispos!=0.&& detaspos!=0.) fHistapaXiposDetaSDphiS->Fill(dphispos,detaspos); 
            if (dphisneg!=0.&& detasneg!=0.) fHistapaXinegDetaSDphiS->Fill(dphisneg,detasneg);
//            }
          }


        } else {//Mixed-event pair histogramming
/*         if (dphisbac==0.&& detasbac==0.) cout<<"Mixing and Deta nad phi bac are 0!!!!"<<endl;
          if (dphispos==0.&& detaspos==0.) cout<<"Mixing and Deta nad phi pos are 0!!!!"<<endl;
          if (dphisneg==0.&& detasneg==0.) cout<<"Mixing and Deta nad phi neg are 0!!!!"<<endl;
*/ 

          if (isXi&&isP) {
  //          if (TMath::Abs(dphispos)>0.015&&TMath::Abs(detaspos)>0.1) {
            fHistpXiSignalBkgKstar->Fill(pairKstar,lcentrality);//}
            if (fReadMCTruth) if (lcentrality>0. && lcentrality<=10.) fHistpXiSignalMixedKstargenvsrec->Fill(pairKstargen,pairKstar);            
//            if (pairKstar<0.1) {
            if (dphisbac!=0.&& detasbac!=0.) fHistpXibacDetaSDphiSBkg->Fill(dphisbac,detasbac);
            if (dphispos!=0.&& detaspos!=0.) fHistpXiposDetaSDphiSBkg->Fill(dphispos,detaspos);
            if (dphisneg!=0.&& detasneg!=0.) fHistpXinegDetaSDphiSBkg->Fill(dphisneg,detasneg);
//            }
          } else if (isXi&&isaP) {
            fHistapXiSignalBkgKstar->Fill(pairKstar,lcentrality);
            if (fReadMCTruth) if (lcentrality>0. && lcentrality<=10.) fHistapXiSignalMixedKstargenvsrec->Fill(pairKstargen,pairKstar); 
//            if (pairKstar<0.1) {
            if (dphisbac!=0.&& detasbac!=0.) fHistapXibacDetaSDphiSBkg->Fill(dphisbac,detasbac);
            if (dphispos!=0.&& detaspos!=0.) fHistapXiposDetaSDphiSBkg->Fill(dphispos,detaspos);
            if (dphisneg!=0.&& detasneg!=0.) fHistapXinegDetaSDphiSBkg->Fill(dphisneg,detasneg);
//            }
          } else if (isaXi&&isP) {
            fHistpaXiSignalBkgKstar->Fill(pairKstar,lcentrality);
            if (fReadMCTruth) if (lcentrality>0. && lcentrality<=10.) fHistpaXiSignalMixedKstargenvsrec->Fill(pairKstargen,pairKstar);
//            if (pairKstar<0.1) {
            if (dphisbac!=0.&& detasbac!=0.) fHistpaXibacDetaSDphiSBkg->Fill(dphisbac,detasbac);
            if (dphispos!=0.&& detaspos!=0.) fHistpaXiposDetaSDphiSBkg->Fill(dphispos,detaspos);
            if (dphisneg!=0.&& detasneg!=0.) fHistpaXinegDetaSDphiSBkg->Fill(dphisneg,detasneg);
//            }

          } else if (isaXi&&isaP) {
  //          if (TMath::Abs(dphisneg)>0.015&&TMath::Abs(detasneg)>0.1) {
            fHistapaXiSignalBkgKstar->Fill(pairKstar,lcentrality);//}
            if (fReadMCTruth) if (lcentrality>0. && lcentrality<=10.) fHistapaXiSignalMixedKstargenvsrec->Fill(pairKstargen,pairKstar);
//            if (pairKstar<0.1) {
            if (dphisbac!=0.&& detasbac!=0.) fHistapaXibacDetaSDphiSBkg->Fill(dphisbac,detasbac);
            if (dphispos!=0.&& detaspos!=0.) fHistapaXiposDetaSDphiSBkg->Fill(dphispos,detaspos);
            if (dphisneg!=0.&& detasneg!=0.) fHistapaXinegDetaSDphiSBkg->Fill(dphisneg,detasneg);
//            }
          }

        }//end mixed event pair histogramming

      }


     }//end event loop
     if (evmultmixed!=0) multmixedcounted = kTRUE; // we count only the first time we go through all events and check if there are protons to mix, then the ncounts will always be the same

   }

   if  (multmixedcounted) fHistMultiplicityOfMixedEvent->Fill(evmultmixed);


}

//_______________________________________________________________

void AliAnalysisTaskhCascadeFemto::DoPairshh (const Float_t lcentrality, int fieldsign) {

  bool is1P;
  bool is1aP;
  bool is2P;
  bool is2aP;
  double dphis;
  Double_t dphis2;
  Double_t dphisprop;
  Double_t detasprop;
  double detas;
  double deta;

  int evmultmixed = 0;
  bool multmixedcounted = kFALSE;
  double pairKstar = 0.;

  int im;
  int jn;

  for (int i=0; i<fEvt->fNumberCandidateProton; i++) {

    if (fEvt->fReconstructedProton[i].doSkipOver) continue;

    is1P  = fEvt->fReconstructedProton[i].isP;
    is1aP = fEvt->fReconstructedProton[i].isaP;

    for (int eventNumber=0; eventNumber<fnEventsToMix+1; eventNumber++) { 
      // For same event pairs
      if (!multmixedcounted && eventNumber!=0 && ((fEvt+eventNumber)->fNumberCandidateProton)!=0.) evmultmixed++; 
      for (int j=0; j<(fEvt+eventNumber)->fNumberCandidateProton; j++) {
        //if(eventNumber!=0) cout<<" fEvt+ev number "<<fEvt+eventNumber<<" event number "<<eventNumber<<" fevt "<<fEvt<<endl;

        if ((fEvt+eventNumber)->fReconstructedProton[j].doSkipOver) continue;
        //cout<<" event number "<<eventNumber<<endl;
        if ( (eventNumber == 0) && (j<=i)) continue; 
        is2P = (fEvt+eventNumber)->fReconstructedProton[j].isP;
        is2aP = (fEvt+eventNumber)->fReconstructedProton[j].isaP;

        // Instead of loop variables use im and jn to do the swapping in a clean way  
        im = i;
        jn = j;
        // Pair ramdomization 
        if (eventNumber == 0) {
          if (gRandom->Rndm()>=0.5) { 
            im = j; jn = i;     
          }  
        }
        //Calculate k* for the pair
        pairKstar = CalculateKstar(fEvt->fReconstructedProton[im].pMomentum, (fEvt+eventNumber)->fReconstructedProton[jn].pMomentum, fPDGsecond,fPDGfirst);

        /*Double_t tot1Mom = TMath::Sqrt(
                                      fEvt->fReconstructedProton[im].pMomentum[0]*
                                      fEvt->fReconstructedProton[im].pMomentum[0] +
                                      fEvt->fReconstructedProton[im].pMomentum[1]*
                                      fEvt->fReconstructedProton[im].pMomentum[1] +
                                      fEvt->fReconstructedProton[im].pMomentum[2]*
                                      fEvt->fReconstructedProton[im].pMomentum[2]);
 
        Double_t tot2Mom = TMath::Sqrt(
                                      (fEvt+eventNumber)->fReconstructedProton[jn].pMomentum[0]*
                                      (fEvt+eventNumber)->fReconstructedProton[jn].pMomentum[0] +
                                      (fEvt+eventNumber)->fReconstructedProton[jn].pMomentum[1]*
                                      (fEvt+eventNumber)->fReconstructedProton[jn].pMomentum[1] +
                                      (fEvt+eventNumber)->fReconstructedProton[jn].pMomentum[2]*
                                      (fEvt+eventNumber)->fReconstructedProton[jn].pMomentum[2]);

        Double_t ediff = TMath::Sqrt(tot1Mom*tot1Mom+fPDGsecond*fPDGsecond)-TMath::Sqrt(tot2Mom*tot2Mom+fPDGsecond*fPDGsecond);

        Double_t qinv = TMath::Sqrt(TMath::Abs(
                       (fEvt->fReconstructedProton[im].pMomentum[0]-(fEvt+eventNumber)->fReconstructedProton[jn].pMomentum[0])*
                       (fEvt->fReconstructedProton[im].pMomentum[0]-(fEvt+eventNumber)->fReconstructedProton[jn].pMomentum[0]) +
                       (fEvt->fReconstructedProton[im].pMomentum[1]-(fEvt+eventNumber)->fReconstructedProton[jn].pMomentum[1])*
                       (fEvt->fReconstructedProton[im].pMomentum[1]-(fEvt+eventNumber)->fReconstructedProton[jn].pMomentum[1]) +
                       (fEvt->fReconstructedProton[im].pMomentum[2]-(fEvt+eventNumber)->fReconstructedProton[jn].pMomentum[2])*
                       (fEvt->fReconstructedProton[im].pMomentum[2]-(fEvt+eventNumber)->fReconstructedProton[jn].pMomentum[2])-
                       (ediff*ediff)));
*/        //cout<<"qinv from kstar "<<pairKstar*2<<" qinv direct "<<qinv<<endl;
        //if (pairKstar==0.) cout<<" K* is zero!!!! "<<endl;


        dphis = CalculateDphiSatR12mAnal((fEvt+eventNumber)->fReconstructedProton[jn].pCharge,
                                          fEvt->fReconstructedProton[im].pCharge, fieldsign, 
                                         (fEvt+eventNumber)->fReconstructedProton[jn].pPt, 
                                          fEvt->fReconstructedProton[im].pPt , 
                                         (fEvt+eventNumber)->fReconstructedProton[jn].pPhi, 
                                          fEvt->fReconstructedProton[im].pPhi, &dphis2);

        deta = fEvt->fReconstructedProton[im].pEta-(fEvt+eventNumber)->fReconstructedProton[jn].pEta;

        dphisprop = CalculateDphiSatR12m( fEvt->fReconstructedProton[im].pShiftedGlobalPosition, 
                                         (fEvt+eventNumber)->fReconstructedProton[jn].pShiftedGlobalPosition );
        detasprop = fEvt->fReconstructedProton[im].pEtaS-(fEvt+eventNumber)->fReconstructedProton[jn].pEtaS;


        //if ((fEvt+eventNumber)->fReconstructedProton[jn].pShiftedGlobalPosition[0]==-9999. && 
        //     fEvt->fReconstructedProton[im].pShiftedGlobalPosition[0]==-9999.) 
        //  cout<<" dphisprop  "<<dphisprop<<" detasprop "<<detasprop<<endl;

        // Apply two-track cuts
        if (fkApplyTtc) { 
          if ( (is1P && is2P) || (is1aP && is2aP) ) {
            if (fSecondpart == kProton) {
              if (fkCutOnTtcProp) {
                if ((dphisprop!=0.) && (detasprop!=0.) && (TMath::Abs(dphisprop)<fDphisMin) && (TMath::Abs(detasprop)<fDetasMin)) continue;
              } else {
                if (TMath::Abs(dphis)<fDphisMin && TMath::Abs(deta)<fDetasMin) continue;  
              } 
            } else if (fSecondpart == kPion) {
              if (fkCutOnTtcProp) {
                if ((dphisprop!=0.) && (detasprop!=0.) && TMath::Sqrt(dphisprop*dphisprop+detasprop*detasprop)<fDphisMin && TMath::Abs(detasprop)<fDetasMin ) continue; 
              } else {
                if (TMath::Sqrt(dphis*dphis+deta*deta)<fDphisMin && TMath::Abs(deta)<fDetasMin) continue; // cut paper
              }
            }
          } 
        }

        if (eventNumber==0) {//Same event pair histogramming
           // simplest pair cut // FIXME not needed maybe
          if ((fEvt+eventNumber)->fReconstructedProton[jn].index == fEvt->fReconstructedProton[im].index) { 
            cout<<"In the same event the two particles have the same index!"<<endl; 
            continue;
          }

          if (is1P && is2P) {
            fHistpXiSignalRealKstar->Fill(pairKstar,lcentrality);
            fHistpXibacDetaSDphiS->Fill(dphis,deta);
            fHistpXiposDetaSDphiS->Fill(dphis2,deta);
            if ( (dphisprop!=0.) && (detasprop!=0.))
              fHistpXinegDetaSDphiS->Fill(dphisprop,detasprop);

          } else if (is1aP && is2aP) {

            fHistapaXiSignalRealKstar->Fill(pairKstar,lcentrality);
            fHistapaXibacDetaSDphiS->Fill(dphis,deta);
            fHistapaXiposDetaSDphiS->Fill(dphis2,deta);
            if ((dphisprop!=0.) && (detasprop!=0.))
              fHistapaXinegDetaSDphiS->Fill(dphisprop,detasprop); 
          } else if ((is1P && is2aP) || (is1aP && is2P)) { // only because we are combining same particles

            fHistpaXiSignalRealKstar->Fill(pairKstar,lcentrality);
            fHistpaXibacDetaSDphiS->Fill(dphis,deta);
            fHistpaXiposDetaSDphiS->Fill(dphis2,deta);
            if ((dphisprop!=0.) && (detasprop!=0.))
              fHistpaXinegDetaSDphiS->Fill(dphisprop,detasprop);

          } 


        } else {//Mixed-event pair histogramming

          if (is1P && is2P) {

            fHistpXiSignalBkgKstar->Fill(pairKstar,lcentrality);
            fHistpXibacDetaSDphiSBkg->Fill(dphis,deta);
            fHistpXiposDetaSDphiSBkg->Fill(dphis2,deta);
            if ((dphisprop!=0.) && (detasprop!=0.))
              fHistpXinegDetaSDphiSBkg->Fill(dphisprop,detasprop);            


          } else if (is1aP && is2aP) {

            fHistapaXiSignalBkgKstar->Fill(pairKstar,lcentrality);
            fHistapaXibacDetaSDphiSBkg->Fill(dphis,deta);
            fHistapaXiposDetaSDphiSBkg->Fill(dphis2,deta);
            if ((dphisprop!=0.) && (detasprop!=0.))
              fHistapaXinegDetaSDphiSBkg->Fill(dphisprop,detasprop);

          } else if ((is1P && is2aP) || (is1aP && is2P)) {

            fHistpaXiSignalBkgKstar->Fill(pairKstar,lcentrality);
            fHistpaXibacDetaSDphiSBkg->Fill(dphis,deta);
            fHistpaXiposDetaSDphiSBkg->Fill(dphis2,deta);
            if ((dphisprop!=0.) && (detasprop!=0.))
              fHistpaXinegDetaSDphiSBkg->Fill(dphisprop,detasprop);
          } 

        }
      } // second part

    }//end event loop

    if (evmultmixed!=0) multmixedcounted = kTRUE;
    
  } // first part

  if  (multmixedcounted) fHistMultiplicityOfMixedEvent->Fill(evmultmixed);
  
}

//_______________________________________________________________

AliReconstructedProton::MCProtonOrigin_t AliAnalysisTaskhCascadeFemto::ProtonOrigin(int trackLabel, TClonesArray *arrayMC, double* pMomentumTruth) {

// for proton purity: take DCA from global tracks (because for TPC only prim and sec are too similar)
// plot the dcaxy for primary, sec from weak decays and material in pt bins. from these templates one fits 
// data with the ROOT TFractionFitter. From fits fractions vs pt for each contribution can be derived

   AliAODMCParticle *partMC = (AliAODMCParticle*) arrayMC->At(trackLabel);
   if (!partMC) {
     //Printf("Pointer to = 0x0 ! Skip ...\n");
     return AliReconstructedProton::kUnassigned;
   }
   if (!(partMC->PxPyPz(pMomentumTruth))) return AliReconstructedProton::kUnassigned;
   if ( partMC->IsPhysicalPrimary()) {
     if ( (fFirstpart == kProton && TMath::Abs(partMC->GetPdgCode()) == 2212) ||
          (fFirstpart == kPion && TMath::Abs(partMC->GetPdgCode()) == 211 ) ) return AliReconstructedProton::kPrimaryP;   
   } else if (partMC->IsSecondaryFromMaterial()) {
     if ( (fFirstpart == kProton && TMath::Abs(partMC->GetPdgCode()) == 2212) ||
          (fFirstpart == kPion && TMath::Abs(partMC->GetPdgCode()) == 211 ) ) { return AliReconstructedProton::kMaterial;  } 
   } else if (partMC->IsSecondaryFromWeakDecay()) { 
     if ( (fFirstpart == kProton &&TMath::Abs(partMC->GetPdgCode()) == 2212 ) ||
          (fFirstpart == kPion && TMath::Abs(partMC->GetPdgCode()) == 211 ) ) { return AliReconstructedProton::kSecondaryWeak; }   
   }

   return AliReconstructedProton::kUnassigned;

}

//_______________________________________________________________

AliReconstructedXi::MCXiOrigin_t AliAnalysisTaskhCascadeFemto::XiOrigin(int labelB, int labelP, int labelN, TClonesArray *arrayMC, double *xiMomentumTruth) {

        // From daughter indeces look for the mother...
        AliAODMCParticle* mcBach  = (AliAODMCParticle*) arrayMC->At( labelB );
        AliAODMCParticle* mcPosV0Dghter   = (AliAODMCParticle*) arrayMC->At( labelP );
        AliAODMCParticle* mcNegV0Dghter   = (AliAODMCParticle*) arrayMC->At( labelN );

        // ID of daughter tracks 
        if( !(mcBach->PdgCode() == -211) && !(mcBach->PdgCode() ==  211) && !(mcBach->PdgCode() == -321) && !(mcBach->PdgCode() ==  321) ) {
         //cout<<" Bach not ok !!"<<endl;
         return AliReconstructedXi::kFake;  
        } 
        if( !(mcPosV0Dghter->PdgCode() == 211) && !(mcPosV0Dghter->PdgCode() == 2212 )) {
         //cout<<" Pos not ok !!"<<endl;
         return AliReconstructedXi::kFake;
        }
        if( !(mcNegV0Dghter->PdgCode() ==-211) && !(mcPosV0Dghter->PdgCode() == -2212)) {
          //cout<<" Neg not ok !!"<<endl;
          return AliReconstructedXi::kFake; 
        }
 
        // Check the lambda
        Int_t lblMotherPosV0Dghter = mcPosV0Dghter->GetMother();
        Int_t lblMotherNegV0Dghter = mcNegV0Dghter->GetMother();

        if( lblMotherPosV0Dghter != lblMotherNegV0Dghter) return AliReconstructedXi::kFake; // same mother 
        if( lblMotherPosV0Dghter < 0 ) return AliReconstructedXi::kFake; // this particle is primary, no mother

        //cout<<"V0 same mother "<<lblMotherPosV0Dghter<<"  "<<lblMotherNegV0Dghter<<endl; 
        AliAODMCParticle* mcMotherPosV0Dghter = (AliAODMCParticle*) arrayMC->At( lblMotherPosV0Dghter );

        Int_t lblGdMotherPosV0Dghter = mcMotherPosV0Dghter->GetMother() ;
        
        //if( lblGdMotherPosV0Dghter != lblGdMotherNegV0Dghter ) {  cout<<"this never happens if we are here... "<<endl; return AliReconstructedXi::kFake;}   
        if( lblGdMotherPosV0Dghter < 0 ) { /*cout<<"Label gm pos v0 is negative! "<<endl;*/ return AliReconstructedXi::kFake;} // primary lambda ...

        AliAODMCParticle* mcGdMotherPosV0Dghter = (AliAODMCParticle*) arrayMC->At( lblGdMotherPosV0Dghter );

        Int_t lblMotherBach = (Int_t) TMath::Abs( mcBach->GetMother() );
        if (lblMotherBach<0) { /*cout<<"bach has negative label !! "<<endl;*/ return AliReconstructedXi::kFake;} 

        if( lblMotherBach != lblGdMotherPosV0Dghter ) { /*cout<<"V0 and bach have different mother !! "<<endl;*/ return AliReconstructedXi::kFake;} //same mother for bach and V0 daughters

        AliAODMCParticle* mcMotherBach = (AliAODMCParticle*) arrayMC->At( lblMotherBach ); // redundant because at this point this is the same mother of the V0 
        //cout<<"all checks ok!"<<endl;

        xiMomentumTruth[0] = mcGdMotherPosV0Dghter->Px();
        xiMomentumTruth[1] = mcGdMotherPosV0Dghter->Py();
        xiMomentumTruth[2] = mcGdMotherPosV0Dghter->Pz();

        //if (mcMotherBach                ->GetPdgCode() != mcGdMotherPosV0Dghter       ->GetPdgCode()) cout<<"This never happens "<<endl; 
        if( mcMotherBach->GetPdgCode() == 3312 ) { 
          if (mcMotherBach->IsPhysicalPrimary())  { //cout<<"PRIMARY XI!! "<<endl;    
            return AliReconstructedXi::kPrimaryXi;
          }
          else { //cout<<"SECONDARY XI!! "<<endl; 
            return AliReconstructedXi::kSecondaryXi; 
          }
        } else if( mcMotherBach->GetPdgCode() == -3312 ) {
          if (mcMotherBach->IsPhysicalPrimary()) { //cout<<"PRIMARY ANTI XI!! "<<endl;
            return AliReconstructedXi::kPrimaryAntiXi;
          } else { //cout<<"SECONDARY ANTI XI!! "<<endl; 
            return AliReconstructedXi::kSecondaryAntiXi; 
          }
        } else if( mcMotherBach->GetPdgCode() == 3334 ) {
          if (mcMotherBach->IsPhysicalPrimary())  { //cout<<"PRIMARY OMEGA!! "<<endl;
            return AliReconstructedXi::kPrimaryOmega;
          } else { //cout<<"SECONDARY OMEGA!! "<<endl; 
            return AliReconstructedXi::kSecondaryOmega;
          }
        } else if( mcMotherBach->GetPdgCode() == -3334 ) {
          if (mcMotherBach->IsPhysicalPrimary())  { //cout<<"SECONDARY ANTI OMEGA!! "<<endl;
            return AliReconstructedXi::kPrimaryAntiOmega;
          } else { //cout<<"SECONDARY ANTI OMEGA!! "<<endl;
            return AliReconstructedXi::kSecondaryAntiOmega; 
          }
        }

        return AliReconstructedXi::kFake;

}

//_______________________________________________________________

Int_t AliAnalysisTaskhCascadeFemto::CheckDaughterTrack (Int_t xiIndex, Int_t daughterIndex, Int_t* checkeddaughters ) {

   //cout<<"In check method "<<endl;
   Int_t nXiWithSharedDaughtersWithXiMass = 0;
   Bool_t kShared = kFALSE;
   for (int j=xiIndex+1; j < fEvt->fNumberCandidateXi; j++) { 
     if (fEvt->fReconstructedXi[j].isXiMass == kTRUE && !(fEvt->fReconstructedXi[j].doSkipOver)) {
       // check if this track was already checked, if so skip this track (return)
       for (int icd = 0; icd < (fEvt->fNumberCandidateXi)*3; icd++) { 
         if (daughterIndex == checkeddaughters[icd]) { 
           //cout<<"This daughter was already checked exiting with -1 and label is "<<daughterIndex<<endl; 
           return -1;
         } 
       }  
       //cout<<"in loop checking"<<endl; 
       // compare with all daughters of all the other Xi
       Int_t indexB2 = fEvt->fReconstructedXi[j].indexB;
       Int_t indexP2 = fEvt->fReconstructedXi[j].indexP;
       Int_t indexN2 = fEvt->fReconstructedXi[j].indexN;
       if ( daughterIndex == indexB2) { //cout<<"multiple use of this now is bac ! "<<endl; 
         fEvt->fReconstructedXi[j].doPickOne = kTRUE; 
         fEvt->fReconstructedXi[xiIndex].doPickOne = kTRUE;
         kShared = kTRUE;
       } else if ( daughterIndex == indexP2) { //cout<<"multiple use of this now is pos ! "<<endl; 
         fEvt->fReconstructedXi[j].doPickOne = kTRUE; 
         fEvt->fReconstructedXi[xiIndex].doPickOne = kTRUE;
         kShared = kTRUE;
       } else if ( daughterIndex == indexN2) { //cout<<"multiple use of this now is neg! "<<endl; 
         fEvt->fReconstructedXi[j].doPickOne = kTRUE; 
         fEvt->fReconstructedXi[xiIndex].doPickOne = kTRUE;
         kShared = kTRUE;
       } 
   
     }
   }

   if (kShared) {  // at least one Xi with one common daughter was found
     SelectBestCandidate (); 
     for (int j=0; j < fEvt->fNumberCandidateXi; j++) {
       if (fEvt->fReconstructedXi[j].isXiMass == kTRUE) {
         if (fEvt->fReconstructedXi[j].doPickOne == kTRUE) {
           nXiWithSharedDaughtersWithXiMass++;
           fEvt->fReconstructedXi[j].doPickOne = kFALSE; // reset
         }
       }
     }
   }
   

   return nXiWithSharedDaughtersWithXiMass;

}

//_______________________________________________________________

void AliAnalysisTaskhCascadeFemto::SelectBestCandidate () {

 // select for each track the best candidate and flag the worst cascades (doSkipOver)
 // to be improved 
 //    - checking bad Xi with good Xi to unflag them if they do not share daughters any more with the good ones (Jai AN Notes/node/64)
 //      (note that the same flag is used to tag Xi that have the primary proton as daughter but this is done after this check in DoPairs so should be safe)     
 //    - consider other criteria
 //    - validate criteria using MC

 Float_t deltamassi = 0.;
 Float_t deltamassj = 0.;
 for (int i=0; i < fEvt->fNumberCandidateXi; i++) { 

   if (fEvt->fReconstructedXi[i].doPickOne && !(fEvt->fReconstructedXi[i].doSkipOver)) {
   
     deltamassi = TMath::Abs(fEvt->fReconstructedXi[i].xiMass - fPDGsecond);
     //cout<<" Delta mass is for first "<<deltamassi<<endl;
     for (int j=i+1; j < fEvt->fNumberCandidateXi; j++) {
       if (fEvt->fReconstructedXi[j].doPickOne) {
         deltamassj = TMath::Abs(fEvt->fReconstructedXi[j].xiMass - fPDGsecond);
         //cout<<" Delta mass is for second "<<deltamassj<<endl;
         if (deltamassj > deltamassi) { 
         fEvt->fReconstructedXi[j].doSkipOver = kTRUE;
         //cout<<"Taking first "<<endl;
         } else  { fEvt->fReconstructedXi[i].doSkipOver = kTRUE; /*cout<<"Taking second "<<endl;*/ i = j;break;} // continue the check from j 
       }// check which is best and set flag to skip
     } 
   }
 }
 return;

}

//_______________________________________________________________

double AliAnalysisTaskhCascadeFemto::CalculateKstar(double momentum1[3], double momentum2[3], double mass1, double mass2) { // Jai S

   // Calculate k* for any pair of particles, regardless of whether the

   // particles have the same mass.

   double kstar = 0.;

   double e1 = 0.;

   double e2 = 0.;

   for(int i = 0; i < 3; i++){

     kstar -= pow(momentum1[i]-momentum2[i],2);

     e1 += pow(momentum1[i],2);

     e2 += pow(momentum2[i],2);

   }

   e1 += pow(mass1,2);

   e1 = sqrt(e1);

   e2 += pow(mass2,2);

   e2 = sqrt(e2);

   

   kstar += pow(e1-e2,2);

 

   double totalMomentumSquared = 0;

   for(int i = 0; i < 3; i++){

     totalMomentumSquared -= pow(momentum1[i]+momentum2[i],2);

   }

   totalMomentumSquared += pow(e1+e2,2);

   kstar -= pow((pow(mass1,2)-pow(mass2,2)),2)/totalMomentumSquared;

 

   kstar *= -1.;

   kstar = sqrt(kstar); //At this point, we've actually calculated Qinv

   kstar *= 0.5; // kstar is 0.5*Qinv

   return kstar;
}

//_______________________________________________________________

double AliAnalysisTaskhCascadeFemto::CalculateDphiSatR12mAnal(Short_t chg1, Short_t chg2, Int_t magSign, Double_t ptv1, Double_t ptv2, Double_t phi1, Double_t phi2, Double_t *dps2) { 

  double rad = 1.2;
  double afsi1b = 0.075*chg1*magSign*rad/ptv1; // 0.07510020733 = - 0.3 (= e in Heaviside-Lorentz units) *0.5 (= B in T) /2 (see later for the -), pT in GeV/c
  double afsi2b = 0.075*chg2*magSign*rad/ptv2;

  if (fabs(afsi1b) >=1.) return 9999.; // angle is pi/2 or not defined --> dont cut
  if (fabs(afsi2b) >=1.) return 9999.; // MN modified these two lines returning 9999 and not kTRUE
  double dps = phi2 - phi1 -TMath::ASin(afsi2b) + TMath::ASin(afsi1b); // - sign of e is outside
  *dps2 = phi2 - phi1 +TMath::ASin(afsi2b) - TMath::ASin(afsi1b);

  dps = TVector2::Phi_mpi_pi(dps);
 
  *dps2 = TVector2::Phi_mpi_pi(*dps2); 
 
  // The following is equivalent
/*  double phi1bis =0.;
  double phi2bis =0.;
  phi1bis = phi1-TMath::ASin(afsi1b); // the minus sign is outside ASin // which procedure is now correct?
  if(phi1bis > 2*PI) phi1bis -= 2*PI;  // alice conv is phi in 0,2 pi
  if(phi1bis < 0) phi1bis += 2*PI;
  phi2bis = phi2 - TMath::ASin(afsi2b);
  if(phi2bis > 2*PI) phi2bis -= 2*PI;
  if(phi2bis < 0) phi2bis += 2*PI;
  double deltaphi = phi2bis - phi1bis;
  if(deltaphi > PI) deltaphi -= 2*PI;  // NB in Dhevan code on git the 2 is omitted! 
//    cout<<"dphi dev >pi!! "<<deltaphi<<endl; 
  if(deltaphi < -PI) deltaphi += 2*PI;
//    cout<<"dphi dev <-pi!! "<<deltaphi<<endl;

  cout<<"dphi bis "<<deltaphi<<" dphi used "<<dps<<endl; 
*/
  return dps; //deltaphi;

}

//_______________________________________________________________

double AliAnalysisTaskhCascadeFemto::CalculateDphiSatR12m(Double_t pos1SftR125[3], Double_t pos2SftR125[3]) { // Hans B

  // Returns delta phi star at R = 1.2 m

  const Float_t distSft = TMath::Sqrt(TMath::Power(pos1SftR125[0] - pos2SftR125[0],2)
                          + TMath::Power(pos1SftR125[1] - pos2SftR125[1],2));
  return 2.0 * TMath::ATan(distSft/2./(125.)); 

}

//_______________________________________________________________

void AliAnalysisTaskhCascadeFemto::SetSftPosR125(AliVTrack *track, const Float_t bfield, Double_t priVtx[3], Double_t posSftR125[3] ) {  // Hans B
  // Sets the spatial position of the track at the radius R=1.25m in the shifted coordinate system
  
  
  // Initialize the array to something indicating there was no propagation
  posSftR125[0]=-9999.;
  posSftR125[1]=-9999.;
  posSftR125[2]=-9999.;
   // Make a copy of the track to not change parameters of the track
  AliExternalTrackParam etp;
  etp.CopyFromVTrack(track);
  
  // The global position of the track
  Double_t xyz[3]={-9999.,-9999.,-9999.};  

  // The radius we want to propagate to, squared, for faster code
  const Float_t rSquared = 125.*125.;


  // Propagation is done in local x of the track
  for (Float_t x = 58.;x<247.;x+=1.){
    // Starts at 83 / Sqrt(2) and goes outwards. 85/Sqrt(2) is the smallest local x
    // for global radius 85 cm. x = 245 is the outer radial limit of the TPC when
    // the track is straight, i.e. has inifinite pt and doesn't get bent. 
    // If the track's momentum is smaller than infinite, it will develop a y-component,
    // which adds to the global radius
    // We don't change the propagation steps to not mess up things!

    // Stop if the propagation was not succesful. This can happen for low pt tracks
    // that don't reach outer radii
    if (!etp.PropagateTo(x,bfield)) { 
      //cout<<"propagation failed!! and coords are -9999. Pt of thetrack "<<track->Pt()<<endl; // very low pT
      break;
    }
    etp.GetXYZ(xyz); // GetXYZ returns global coordinates  
    // Calculate the shifted radius we are at, squared. 
    // Compare squared radii for faster code
    Float_t shiftedRadiusSquared = (xyz[0]-priVtx[0])*(xyz[0]-priVtx[0])
                                 + (xyz[1]-priVtx[1])*(xyz[1]-priVtx[1]);

    // Roughly reached the radius we want
    if(shiftedRadiusSquared > rSquared){
      
      // Bigger loop has bad precision, we're nearly one centimeter too far, 
      // go back in small steps.
      while (shiftedRadiusSquared>rSquared) {
        // Propagate a mm inwards
        x-=.1;
        if (!etp.PropagateTo(x,bfield)){
          // Propagation failed but we're already with a
          // cm precision at R=1.25m so we only break the 
          // inner loop
          //cout<<"propagation failed!! and etss is "<<EtaS(posSftR125)<<endl; 

          break;
        }
        // Get the global position
        etp.GetXYZ(xyz);
        // Calculate shifted radius, squared
        shiftedRadiusSquared = (xyz[0]-priVtx[0])*(xyz[0]-priVtx[0])
                           + (xyz[1]-priVtx[1])*(xyz[1]-priVtx[1]);
      }
      // We reached R=1.25m with a precission of a cm to a mm,
      // set the spatial position
      posSftR125[0]=xyz[0]-priVtx[0];
      posSftR125[1]=xyz[1]-priVtx[1];
      posSftR125[2]=xyz[2]-priVtx[2];
      //cout<<" Pos 125 cm in function end "<<posSftR125[0]<<" "<<posSftR125[1]<<" "<<posSftR125[2]<<endl;
      // Done
      return;
    } // End of if roughly reached radius
  } // End of coarse propagation loop
}

//_______________________________________________________________

void AliAnalysisTaskhCascadeFemto::SetPosR125(AliVTrack *track, const Float_t bfield, Double_t posSftR125[3] ) {  
  // Sets the spatial position of the track at the radius R=1.25m in the shifted coordinate system

  // New method to propagate to a radius, seems fails less and is more precise at large R. 
  // Initialize the array to something indicating there was no propagation
  posSftR125[0]=-9999.;
  posSftR125[1]=-9999.;
  posSftR125[2]=-9999.;
   // Make a copy of the track to not change parameters of the track
  AliExternalTrackParam etp;
  etp.CopyFromVTrack(track);

  // The global position of the track
  Double_t xyz[3]={-9999.,-9999.,-9999.};

  // The radius we want to propagate to, squared, for faster code
  const double r = 125.;
  const double rSquared = r*r;
   
  // FIXME how to shift to the primary vertex position? Should it be done iteratively?
  if (!etp.GetXYZatR(r,bfield,xyz,NULL)) return;  
//  Float_t shiftedRadiusSquared = (xyz[0]-priVtx[0])*(xyz[0]-priVtx[0])
//                                 + (xyz[1]-priVtx[1])*(xyz[1]-priVtx[1]);

//  Float_t unshiftedRadiusSquared = (xyz[0]*xyz[0])+ (xyz[1]*xyz[1]); // this is exactly 125 
  //cout<<"Coord at r "<<r<<" x "<<xyz[0]<<" y "<<xyz[1]<<" z "<<xyz[2]<<" shifted r "<<TMath::Sqrt(shiftedRadiusSquared)<<" unshifted "<<TMath::Sqrt(unshiftedRadiusSquared)<<endl;
  posSftR125[0]=xyz[0];//-priVtx[0];
  posSftR125[1]=xyz[1];//-priVtx[1];
  posSftR125[2]=xyz[2];//-priVtx[2];


  return ;
  
}

//_______________________________________________________________

Double_t AliAnalysisTaskhCascadeFemto::ThetaS( Double_t posSftR125[3] ) const { // Hans B
  // Returns the longitudinal angle of the particles propagated
  // position at R=1.25m. See
  // https://edms.cern.ch/file/406391/2/ALICE-INT-2003-038.pdf
  // for the ALICE coordinate system. Theta is zero at positive z,
  // pi/2 at z = 0 aka the xy plane and pi at negative z 

  // R^    ^  
  //  |   /
  //  |?'/
  //  | / ?
  //  |/____>z
  // 
  // Let's compute ?' and ? = pi/2 - ?'
  // where ?' can even be and should 
  // sometimes be negative
  // tan(?') = z/R
  // ?' = arctan(z/R)
  // ? = pi/2 - ?'
  //   = pi/2 - arctan(z/R)
  // Note that in the doc above theta
  // is calculated as arccos(z/sqrt(x^2+y^2+z^2))

  // Array of positions is 85,105,125,..cm,
  // we take the z position at R=1.25m
  // return TMath::Pi()/2. - TMath::ATan(fXshifted[2][2]/125.);
  return TMath::Pi()/2. - TMath::ATan(posSftR125[2]/125.); // ok here R is really there --> transverse plane 
}

//_______________________________________________________________

Double_t AliAnalysisTaskhCascadeFemto::EtaS( Double_t posSftR125[3] ) const {  // Hans B
  // Returns the corresponding eta of a pri. part. 
  // with this particles pos at R=1.25m

  // http://en.wikipedia.org/wiki/Pseudorapidity
  // ? = -ln[ tan(?/2)]
  // printf("z: %+04.0f, thetaS %+03.2f etaS %+1.2f\n"
  // ,fXshifted[2][2],ThetaS(),-TMath::Log( TMath::Tan(ThetaS()/2.) ));
  return -TMath::Log( TMath::Tan(ThetaS(posSftR125 )/2.) );
}

//_______________________________________________________________

void AliAnalysisTaskhCascadeFemto::Terminate(const Option_t *) {
  // Draw result to the screen
  // Called once at the end of the query

  if (!GetOutputData(0)) return;
}

