/************************************************************************* 
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. * 
 *                                                                        * 
 * Author: X. Sanchez Castro                                              * 
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

//git test

#include <TCanvas.h>
#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <THnSparse.h>
#include <TH3F.h>
#include <TPDGCode.h>
#include <TDatabasePDG.h>
#include <TClonesArray.h>
#include <TROOT.h>

#include "AliOADBContainer.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliGenHijingEventHeader.h"

#include "AliAODEvent.h"
#include "AliAODv0.h"
#include "AliAODcascade.h"

#include "AliCFContainer.h"
#include "AliCentrality.h"

#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliAODPid.h"

#include "AliInputEventHandler.h"

#include "AliAnalysisTaskLambdaOverK0sJets.h"

//extern TROOT *gROOT;


ClassImp(AliAnalysisTaskLambdaOverK0sJets)
ClassImp(AliMiniParticle)

// Global variables:
static Int_t    nbins = 100;                 // Number of bins for l, pt, mass for V0
static Int_t    nbinsPhi = 120;              // Number of bins for Phi
static Int_t    nbinsdPhi = 20;              // Number of bins for dPhi
static Int_t    nbinsdEta = 30;              // Number of bins for dEta
static Int_t    nbinPtLP = 200;
static Int_t    nbinsVtx = 20;

static Float_t pMin = 0.0;                   // Lower cut for transverse momentum
static Float_t pMax = 10.;                   // Max cut for transverse momentum for V0
static Float_t ptMaxLP = 200.;               // Max cut for transverse momentum for trigger particle

static Float_t lMin = 0.0;                  // Limits in the histo for fidutial volume
static Float_t lMax = 100.;                 // Limits in the fidutial volume

static Int_t   nMaxEvMix = 250;

//
//  
//

AliAnalysisTaskLambdaOverK0sJets::AliAnalysisTaskLambdaOverK0sJets(const char *name) :
AliAnalysisTaskSE(name),

fAOD(0),  fCollision("PbPb2010"), fIsMC(kFALSE), fUsePID(kFALSE), fCentMin(0.), fCentMax(90.), fDoQA(kFALSE), fDoMixEvt(kFALSE), fTriggerFB(768), fTrigPtMin(5.), fTrigPtMax(10.), fTrigPtMCMin(5.), fTrigPtMCMax(10000.), fTrigEtaMax(0.8), fTriggerNCls(10), fCheckIDTrig(kFALSE), fSeparateInjPart(kTRUE), fEndOfHijingEvent(-1),  fPIDResponse(0),

  fMinPtDaughter(0.160), fMaxEtaDaughter(0.8), fMaxDCADaughter(1.0), fUseEtaCut(kFALSE), fYMax(0.7), fDCAToPrimVtx(0.1), fMinCPA(0.998), fNSigma(3.0), fDaugNClsTPC(70.), fMinCtau(0.), fMaxCtau(3.), fIdTrigger(-1), fIsV0LP(0), fPtV0LP(0.), fIsSndCheck(0),

fTPCRadius(125.), fFracTPCcls(1.0), fDiffTrigDaugFracTPCSharedCls(0.06),

fOutput(0), fOutputQA(0), fOutputME(0), fMEList(0x0), fTriggerParticles(0x0), fTriggerPartMC(0x0), fAssocParticles(0x0), fAssocPartMC(0x0), fEvents(0), fEvtPerCent(0), fCentrality(0),  fCentrality2(0), fCentralityTrig(0), fPrimayVtxGlobalvsSPD(0), fPrimaryVertexX(0), fPrimaryVertexY(0), fPrimaryVertexZ(0), fChargedMultiplicity(0),

fTriggerEventPlane(0),  fTriggerMCPtCent(0), fTriggerMCResPt(0), fTriggerMCResEta(0), fTriggerMCResPhi(0), fTriggerPtCent(0),  fTriggerPtCentCh(0), fNTrigPerEvt(0), fTriggerWiSPDHit(0), fTriggerEtaPhi(0), fTrigFracShTPCcls(0), fTriggerDCA(0), fCheckTriggerFromV0Daug(0), fTriggerComingFromDaug(0), fTriggerIsV0(0), fCheckIDTrigPtK0s(0), fCheckIDTrigPhiK0s(0), fCheckIDTrigEtaK0s(0), fCheckIDTrigNclsK0s(0), fCheckIDTrigPtLambda(0), fCheckIDTrigPhiLambda(0), fCheckIDTrigEtaLambda(0),  fCheckIDTrigNclsLambda(0), fCheckIDTrigPtAntiLambda(0), fCheckIDTrigPhiAntiLambda(0), fCheckIDTrigEtaAntiLambda(0), fCheckIDTrigNclsAntiLambda(0), 

  fInjectedParticles(0),

  fK0sMCPt(0), fK0sMCPtRap(0), fK0sMCPtRap2(0),  fK0sMCPtRapEmbeded(0), fK0sAssocPt(0), fK0sAssocPtArm(0),  fK0sAssocPtRap(0), fK0sAssocPtRapEmbeded(0), fK0sMCResEta(0), fK0sMCResPhi(0), fK0sMCResPt(0), fK0sPosMCResEta(0), fK0sPosMCResPhi(0), fK0sPosMCResPt(0), fK0sNegMCResEta(0), fK0sNegMCResPhi(0), fK0sNegMCResPt(0), fLambdaMCPt(0), fLambdaMCPtRap(0), fLambdaMCPtRap2(0),  fLambdaMCPtRapEmbeded(0),  fLambdaMCFromXi(0), fLambdaAssocPt(0), fLambdaAssocPtRap(0), fLambdaAssocFromXi(0), fLambdaMCResEta(0), fLambdaMCResPhi(0), fLambdaMCResPt(0), fLambdaPosMCResEta(0), fLambdaPosMCResPhi(0), fLambdaPosMCResPt(0), fLambdaNegMCResEta(0), fLambdaNegMCResPhi(0), fLambdaNegMCResPt(0), fAntiLambdaMCPt(0), fAntiLambdaMCPtRap(0), fAntiLambdaMCPtRap2(0), fAntiLambdaMCPtRapEmbeded(0), fAntiLambdaMCFromXi(0), fAntiLambdaAssocPt(0), fAntiLambdaAssocPtRap(0), fAntiLambdaAssocFromXi(0), fAntiLambdaMCResEta(0), fAntiLambdaMCResPhi(0), fAntiLambdaMCResPt(0), fAntiLambdaPosMCResEta(0), fAntiLambdaPosMCResPhi(0), fAntiLambdaPosMCResPt(0), fAntiLambdaNegMCResEta(0), fAntiLambdaNegMCResPhi(0), fAntiLambdaNegMCResPt(0), 

  fHistArmenterosPodolanski(0), fHistArmPodBckg(0),
  
fK0sMass(0), fK0sMassEmbeded(0), fK0sMassPtEta(0), fK0sMassPtPhi(0), fK0sPosDaugFracShTPCcls(0), fK0sNegDaugFracShTPCcls(0), fK0sDaughtersPt(0), fK0sPosDaugFracShTPCclsTrig(0), fK0sNegDaugFracShTPCclsTrig(0),  fK0sDCADaugToPrimVtx(0), fK0sSpatialRes(0), fK0sBckgDecLength(0), fK0sBckgDCADaugToPrimVtx(0), fK0sBckgEtaPhi(0), fK0sBckgPhiRadio(0), fK0sBckgDCANegDaugToPrimVtx(0), fK0sBckgDCAPosDaugToPrimVtx(0), fV0MassCascade(0),
  
fLambdaMass(0), fLambdaMassEmbeded(0), fLambdaMass2(0), fLambdaMass2Embeded(0), fLambdaMassPtEta(0), fLambdaMassPtPhi(0), fLambdaPosDaugFracShTPCcls(0), fLambdaNegDaugFracShTPCcls(0), fLambdaDaughtersPt(0), fLambdaPosDaugFracShTPCclsTrig(0), fLambdaNegDaugFracShTPCclsTrig(0), fLambdaDCADaugToPrimVtx(0), fLambdaSpatialRes(0), fLambdaBckgDecLength(0), fLambdaBckgDCADaugToPrimVtx(0), fLambdaBckgEtaPhi(0), fLambdaBckgPhiRadio(0), fLambdaBckgDCANegDaugToPrimVtx(0), fLambdaBckgDCAPosDaugToPrimVtx(0), 

fAntiLambdaMass(0), fAntiLambdaMassEmbeded(0), fAntiLambdaMass2(0), fAntiLambdaMass2Embeded(0), fAntiLambdaMassPtEta(0), fAntiLambdaMassPtPhi(0), fAntiLambdaPosDaugFracShTPCcls(0), fAntiLambdaNegDaugFracShTPCcls(0), fAntiLambdaDaughtersPt(0),  fAntiLambdaPosDaugFracShTPCclsTrig(0), fAntiLambdaNegDaugFracShTPCclsTrig(0), fAntiLambdaDCADaugToPrimVtx(0), fAntiLambdaSpatialRes(0), fAntiLambdaBckgDecLength(0), fAntiLambdaBckgDCADaugToPrimVtx(0), fAntiLambdaBckgEtaPhi(0), fAntiLambdaBckgPhiRadio(0), fAntiLambdaBckgDCANegDaugToPrimVtx(0), fAntiLambdaBckgDCAPosDaugToPrimVtx(0), 

fK0sPtPosDaug(0), fK0sPtNegDaug(0), fK0sBckgPtPosDaug(0), fK0sBckgPtNegDaug(0), fK0sPhiEtaPosDaug(0), fK0sPhiEtaNegDaug(0), fK0sBckgPhiEtaPosDaug(0), fK0sBckgPhiEtaNegDaug(0), fK0sDCAPosDaug(0), fK0sDCANegDaug(0), fK0sBckgDCAPosDaug(0), fK0sBckgDCANegDaug(0), fK0sDecayPos(0), fK0sBckgDecayPos(0), fK0sDecayVertex(0), fK0sBckgDecayVertex(0), fK0sCPA(0), fK0sBckgCPA(0), fK0sDCAV0Daug(0), fK0sBckgDCAV0Daug(0), fK0sNClustersTPC(0), fK0sBckgNClustersTPC(0), fK0sNClustersITSPos(0), fK0sNClustersITSNeg(0), fK0sBckgNClustersITSPos(0), fK0sBckgNClustersITSNeg(0), fK0sCTau(0), fK0sBckgCTau(0),

  fLambdaPtPosDaug(0), fLambdaPtNegDaug(0), fLambdaBckgPtPosDaug(0), fLambdaBckgPtNegDaug(0), fLambdaPhiEtaPosDaug(0),fLambdaPhiEtaNegDaug(0), fLambdaBckgPhiEtaPosDaug(0), fLambdaBckgPhiEtaNegDaug(0), fLambdaDCAPosDaug(0),fLambdaDCANegDaug(0), fLambdaBckgDCAPosDaug(0), fLambdaBckgDCANegDaug(0), fLambdaDecayPos(0), fLambdaBckgDecayPos(0), fLambdaDecayVertex(0), fLambdaBckgDecayVertex(0), fLambdaCPA(0), fLambdaBckgCPA(0), fLambdaDCAV0Daug(0), fLambdaBckgDCAV0Daug(0), fLambdaNClustersTPC(0), fLambdaBckgNClustersTPC(0), fLambdaNClustersITSPos(0), fLambdaNClustersITSNeg(0), fLambdaBckgNClustersITSPos(0), fLambdaBckgNClustersITSNeg(0), fLambdaCTau(0), fLambdaBckgCTau(0),

fAntiLambdaPtPosDaug(0), fAntiLambdaPtNegDaug(0), fAntiLambdaBckgPtPosDaug(0), fAntiLambdaBckgPtNegDaug(0), fAntiLambdaPhiEtaPosDaug(0),fAntiLambdaPhiEtaNegDaug(0), fAntiLambdaBckgPhiEtaPosDaug(0),fAntiLambdaBckgPhiEtaNegDaug(0), fAntiLambdaDCAPosDaug(0),fAntiLambdaDCANegDaug(0), fAntiLambdaBckgDCAPosDaug(0), fAntiLambdaBckgDCANegDaug(0), fAntiLambdaDecayPos(0), fAntiLambdaBckgDecayPos(0), fAntiLambdaDecayVertex(0), fAntiLambdaBckgDecayVertex(0), fAntiLambdaCPA(0), fAntiLambdaBckgCPA(0), fAntiLambdaDCAV0Daug(0), fAntiLambdaBckgDCAV0Daug(0), fAntiLambdaNClustersTPC(0), fAntiLambdaBckgNClustersTPC(0), fAntiLambdaNClustersITSPos(0), fAntiLambdaNClustersITSNeg(0), fAntiLambdaBckgNClustersITSPos(0),  fAntiLambdaBckgNClustersITSNeg(0), fAntiLambdaCTau(0), fAntiLambdaBckgCTau(0)
  
{
  // Dummy Constructor

  // variables for track splitting:
  // shifted positionf for thw tracks
  for(Int_t i=0; i<3; i++){
    fTrigSftR125[i] = -9999.;
    fDaugSftR125[i] = -9999.;     
  }

  // Particles properties in MC
  for (Int_t i=0; i<kNCent; i++){ 
    
    // K0s
    fK0sMCPtRapVtx[i] = 0;
    fK0sMCPtRapVtxEmbeded[i] = 0;
    fK0sMCPtRapPtDaugPt[i] = 0x0;
    fK0sMCPtRapPtDaugPtEmbeded[i] = 0x0;
    fK0sMCPtPhiEta[i] = 0;
    fK0sAssocPtPhiEta[i] = 0;
    // -- Natural particles
    fK0sAssocPtMassArm[i] = 0x0;
    fK0sAssocMassPtVtx[i] = 0x0;
    fK0sAssocMassPtDCADaug[i] = 0x0;
    fK0sAssocMassPtCPA[i] = 0x0;
    fK0sAssocMassPtDCAPV[i] = 0x0;
    fK0sAssocMassPtDaugNClsTPC[i] = 0x0;
    fK0sAssocMassPtShTPCcls[i] = 0x0;
    fK0sAssocMassPtDaugPt[i] = 0x0;
    fK0sAssocMassPtCtau[i] = 0x0;
    fK0sAssocMassPtFidVolume[i] = 0x0;
    // -- Embeded particles
    fK0sAssocPtMassArmEmbeded[i] = 0x0;
    fK0sAssocMassPtVtxEmbeded[i] = 0x0;
    fK0sAssocMassPtDCADaug[i] = 0x0;
    fK0sAssocMassPtCPAEmbeded[i] = 0x0;
    fK0sAssocMassPtDCAPVEmbeded[i] = 0x0;
    fK0sAssocMassPtDaugNClsTPCEmbeded[i] = 0x0;
    fK0sAssocMassPtShTPCclsEmbeded[i] = 0x0;
    fK0sAssocMassPtDaugPtEmbeded[i] = 0x0;
    fK0sAssocMassPtCtauEmbeded[i] = 0x0;
    fK0sAssocMassPtFidVolumeEmbeded[i] = 0x0;
    // -- Mass vs rapidity vs pt vs centrlaity
    fK0sMassPtRap[i] = 0;
    // -- Splitting checks
    fK0sPosDaugSplCheckCovMat[i] = 0x0;
    fK0sNegDaugSplCheckCovMat[i] = 0x0;
    fK0sPosDaugdPhiSdEtaS[i] = 0x0;   
    fK0sNegDaugdPhiSdEtaS[i] = 0x0;
    fK0sPosMCResdEtaSdPhiS[i] = 0x0;
    fK0sNegMCResdEtaSdPhiS[i] = 0x0;

    // Lambda
    fLambdaMCPtRapVtx[i] = 0;
    fLambdaMCPtRapVtxEmbeded[i] = 0;
    fLambdaMCPtRapPtDaugPt[i] = 0x0;
    fLambdaMCPtRapPtDaugPtEmbeded[i] = 0x0;
    fLambdaMCPtPhiEta[i] = 0;
    fLambdaAssocPtPhiEta[i] = 0;
    // -- Natural particles
    fLambdaAssocMassPtRap[i] = 0x0;
    fLambdaAssocMassPtRap2[i] = 0x0;
    fLambdaAssocMassPtVtx[i] = 0x0;
    fLambdaAssocMassPtDCADaug[i] = 0x0;
    fLambdaAssocMassPtCPA[i] = 0x0;
    fLambdaAssocMassPtDCAPV[i] = 0x0;
    fLambdaAssocMassPtDaugNClsTPC[i] = 0x0;
    fLambdaAssocMassPtShTPCcls[i] = 0x0;
    fLambdaAssocMassPtDaugPt[i] = 0x0;
    fLambdaAssocMassPtCtau[i] = 0x0;
    fLambdaAssocMassPtFidVolume[i] = 0x0;
    // -- Embeded particles
    fLambdaAssocMassPtRapEmbeded[i] = 0x0;
    fLambdaAssocMassPtRapEmbeded2[i] = 0x0;
    fLambdaAssocMassPtVtxEmbeded[i] = 0x0;
    fLambdaAssocMassPtDCADaug[i] = 0x0;
    fLambdaAssocMassPtCPAEmbeded[i] = 0x0;
    fLambdaAssocMassPtDCAPVEmbeded[i] = 0x0;
    fLambdaAssocMassPtDaugNClsTPCEmbeded[i] = 0x0;
    fLambdaAssocMassPtShTPCclsEmbeded[i] = 0x0;
    fLambdaAssocMassPtDaugPtEmbeded[i] = 0x0;
    fLambdaAssocMassPtCtauEmbeded[i] = 0x0;
    fLambdaAssocMassPtFidVolumeEmbeded[i] = 0x0;
    // -- Mass vs rapidity vs pt vs centrlaity
    fLambdaMassPtRap[i] = 0;
    // -- Splitting checks
    fLambdaPosDaugSplCheckCovMat[i] = 0x0;
    fLambdaNegDaugSplCheckCovMat[i] =0x0;
    fLambdaPosDaugdPhiSdEtaS[i] = 0x0;
    fLambdaNegDaugdPhiSdEtaS[i] = 0x0;
    fLambdaPosMCResdEtaSdPhiS[i] = 0x0;
    fLambdaNegMCResdEtaSdPhiS[i] = 0x0;

    // AntiLambda
    fAntiLambdaMCPtRapVtx[i] = 0;
    fAntiLambdaMCPtRapVtxEmbeded[i] = 0;
    fAntiLambdaMCPtRapPtDaugPt[i] = 0x0;
    fAntiLambdaMCPtRapPtDaugPtEmbeded[i] = 0x0;
    fAntiLambdaMCPtPhiEta[i] = 0;
    fAntiLambdaAssocPtPhiEta[i] = 0;
    // -- Natural particles
    fAntiLambdaAssocMassPtRap[i] = 0x0;
    fAntiLambdaAssocMassPtRap2[i] = 0x0;
    fAntiLambdaAssocMassPtVtx[i] = 0x0;
    fAntiLambdaAssocMassPtDCADaug[i] = 0x0;
    fAntiLambdaAssocMassPtCPA[i] = 0x0;
    fAntiLambdaAssocMassPtDCAPV[i] = 0x0;
    fAntiLambdaAssocMassPtDaugNClsTPC[i] = 0x0;
    fAntiLambdaAssocMassPtShTPCcls[i] = 0x0;
    fAntiLambdaAssocMassPtDaugPt[i] = 0x0;
    fAntiLambdaAssocMassPtCtau[i] = 0x0;
    fAntiLambdaAssocMassPtFidVolume[i] = 0x0;
    // -- Embeded particles
    fAntiLambdaAssocMassPtRapEmbeded[i] = 0x0;
    fAntiLambdaAssocMassPtRapEmbeded2[i] = 0x0;
    fAntiLambdaAssocMassPtVtxEmbeded[i] = 0x0;
    fAntiLambdaAssocMassPtDCADaug[i] = 0x0;
    fAntiLambdaAssocMassPtCPAEmbeded[i] = 0x0;
    fAntiLambdaAssocMassPtDCAPVEmbeded[i] = 0x0;
    fAntiLambdaAssocMassPtDaugNClsTPCEmbeded[i] = 0x0;
    fAntiLambdaAssocMassPtShTPCclsEmbeded[i] = 0x0;
    fAntiLambdaAssocMassPtDaugPtEmbeded[i] = 0x0;
    fAntiLambdaAssocMassPtCtauEmbeded[i] = 0x0;
    fAntiLambdaAssocMassPtFidVolumeEmbeded[i] = 0x0;
    // -- Mass vs rapidity vs pt vs centrlaity
    fAntiLambdaMassPtRap[i] = 0;
    // -- Splitting checks
    fAntiLambdaPosDaugSplCheckCovMat[i] = 0x0;
    fAntiLambdaNegDaugSplCheckCovMat[i] = 0x0;
    fAntiLambdaPosDaugdPhiSdEtaS[i] = 0x0;
    fAntiLambdaNegDaugdPhiSdEtaS[i] = 0x0;
    fAntiLambdaPosMCResdEtaSdPhiS[i] = 0x0;
    fAntiLambdaNegMCResdEtaSdPhiS[i] = 0x0;

  }

  // Correlations in MC
  for (Int_t i=0; i<kNCent*kN1; i++){     
    // K0s
    fK0sdPhidEtaMC[i] = 0;
    // Lambda
    fLambdadPhidEtaMC[i] = 0;
    // AntiLambda
    fAntiLambdadPhidEtaMC[i] = 0;
  }

  // Correlations
  for (Int_t i=0; i<(kNCent*kN1*kNVtxZ); i++){     
    // K0s
    fK0sdPhidEtaPtL[i] = 0;
    // Lambda
    fLambdadPhidEtaPtL[i] = 0;
    // AntiLambda
    fAntiLambdadPhidEtaPtL[i] = 0;  
  }

  // Mixed events distributions
  for (Int_t i=0; i<(kN1*kNVtxZ*kNCent); i++){ 
    fK0sdPhidEtaME[i] = 0;
    fLambdadPhidEtaME[i] = 0;
    fAntiLambdadPhidEtaME[i] = 0;
  }
 
  // Constructor. Initialization of pointers
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
  DefineOutput(3, TList::Class());
 

}

//___________________________________________________________________________________________

AliAnalysisTaskLambdaOverK0sJets::~AliAnalysisTaskLambdaOverK0sJets() 
{

  // Destructor
  if(fMEList){
    
    for(Int_t icent=0; icent<kNCent; icent++){
      for(Int_t iz=0; iz<kNVtxZ; iz++){
	fMEList[icent*kNVtxZ+iz]->Delete();  delete fMEList[icent*kNVtxZ+iz];
      }
    }
    delete[] fMEList; fMEList=0x0;
  }
  
  if(fTriggerParticles) {
    delete fTriggerParticles;
    fTriggerParticles=0x0;
  }

  if(fTriggerPartMC) {
    delete fTriggerPartMC;
    fTriggerPartMC=0x0;
  }

  if(fAssocParticles) {
    delete fAssocParticles;
    fAssocParticles=0x0;
  }

  if(fAssocPartMC) {
    delete fAssocPartMC;
    fAssocPartMC=0x0;
  }
 
    
}

//___________________________________________________________________________________________

void AliAnalysisTaskLambdaOverK0sJets::UserCreateOutputObjects()
{ 
  // Creating the histograms that are needed for the output 
  
  fOutput = new TList(); 
  fOutput->SetOwner();

  fOutputQA = new TList(); 
  fOutputQA->SetOwner();

  fOutputME = new TList(); 
  fOutputME->SetOwner();

  fMEList = new TList*[kNCent*kNVtxZ];
  for(Int_t icent=0; icent<kNCent; icent++){
    for(Int_t iz=0; iz<kNVtxZ; iz++){
      fMEList[icent*kNVtxZ+iz] = new TList();
      fMEList[icent*kNVtxZ+iz]->SetOwner(kFALSE);
    }
  }

  char hNameHist[100];

  // ====== General characteristics of the event and tracks ====== //

  // Counter for the number of events in each step:
  fEvents =new TH1F("fEvents","Number of events",14,-0.5,13.5);
  fEvents->GetXaxis()->SetBinLabel(1,"calls to UserExec()");
  fEvents->GetXaxis()->SetBinLabel(2,"AOD available");
  fEvents->GetXaxis()->SetBinLabel(3,"CINT1B");
  fEvents->GetXaxis()->SetBinLabel(4,"V0M Cent");
  fEvents->GetXaxis()->SetBinLabel(5,"Global Vtx > 3 part");
  fEvents->GetXaxis()->SetBinLabel(6,"SPD Vtx > 3 part");
  fEvents->GetXaxis()->SetBinLabel(7,"|ZVtx Global - Zvtx SPD| < 0.5 cm");
  fEvents->GetXaxis()->SetBinLabel(8,"|VtxZ| < 10 cm");
  fEvents->GetXaxis()->SetBinLabel(9,"V0 is LP");
  fEvents->GetXaxis()->SetBinLabel(10," ");
  fEvents->GetXaxis()->SetBinLabel(11," ");
  fEvents->GetXaxis()->SetBinLabel(12,"Triggered");
  fEvents->GetXaxis()->SetBinLabel(13,"NOT Triggered");
  fEvents->GetXaxis()->SetBinLabel(14," ");
  fEvents->GetYaxis()->SetTitle("Counts"); 
  fOutput->Add(fEvents);

  fEvtPerCent = new TH2F("fEvtPerCent","Events per centrality bin;Step;Centrality bin",5,-0.5,4.5,4,-0.5,3.5);
  fOutput->Add(fEvtPerCent);

  // Centrality:
  fCentrality = new TH1F("fCentrality","Centrality;Centrality (%);Events",500,0.,100.);
  fOutput->Add(fCentrality);

  fCentrality2 = new TH1F("fCentrality2","Centrality in events with |VtxZ|<10 cm;Centrality (%);Events",500,0.,100.);
  fOutput->Add(fCentrality2);

  fCentralityTrig = new TH2F("fCentralityTrig","Centrality in events per trigger selection;Centrality (%);Triger Selection",100,0.,100.,3,0.5,3.5);
  fCentralityTrig->GetYaxis()->SetBinLabel(1,"kCentral");
  fCentralityTrig->GetYaxis()->SetBinLabel(2,"kSemiCentral");
  fCentralityTrig->GetYaxis()->SetBinLabel(3,"kMB");
  fOutput->Add(fCentralityTrig);

  // Primary Vertex:
  fPrimayVtxGlobalvsSPD = new TH2F("fPrimayVtxGlobalvsSPD",";Z_{vtx,tr} (cm);Z_{SPD,tr} (cm)",200,-20,20,200,-20,20);
  fOutput->Add(fPrimayVtxGlobalvsSPD);

  fPrimaryVertexX = new TH1F("fPrimaryVertexX", "Primary Vertex Position X;Primary Vertex Position X (cm);Events",100,-0.5,0.5);
  fOutput->Add(fPrimaryVertexX);
  
  fPrimaryVertexY = new TH1F("fPrimaryVertexY", "Primary Vertex Position Y;Primary Vertex Position Y (cm);Events",100,-0.5,0.5);
  fOutput->Add(fPrimaryVertexY);
  
  fPrimaryVertexZ = new TH1F("fPrimaryVertexZ", "Primary Vertex Position Z;Primary Vertex Position Z (cm);Events",200,-20,20);
  fOutput->Add(fPrimaryVertexZ);
  
  fChargedMultiplicity  = new TH2F("fChargedMultiplicity","Multiplicity;Multiplicity;centrality bin",1225,0,2500.,4,0.,4.);
  fOutput->Add(fChargedMultiplicity);

  // ====== Trigger Particle characteristics ====== //
  
  // Difference between Event plane and the Trigger particles:
  fTriggerEventPlane = new TH1F("fTriggerEventPlane", ";#varphi_{EP}-#varphi_{Trig};Events",50,0.,TMath::Pi());
  fOutput->Add(fTriggerEventPlane);

  // MC: Pt Trigger particle vs centrality:
  if(fIsMC){
    fTriggerMCPtCent = new TH2F("fTriggerMCPtCent","Trigger particle MC;#it{p}_{T} (GeV/#it{c});centrality (%)",2*nbinPtLP,pMin,2*ptMaxLP,100,0.,100.);
    fOutput->Add(fTriggerMCPtCent);

    fTriggerMCResPt = new TH3F("fTriggerMCResPt","Trigger particle MC: #it{p}_{T} resolution;(p_{t,MC}-p_{t,Rec})/p_{t,Rec};#it{p}_{T} (GeV/#it{c});centrality",60,-0.3,0.3,2*nbinPtLP,pMin,ptMaxLP,100,0.,100.);
    fOutput->Add(fTriggerMCResPt);

    fTriggerMCResEta = new TH3F("fTriggerMCResEta","Trigger particle MC: #eta resolution; #eta_{MC}-#eta_{Rec};#it{p}_{T} (GeV/#it{c}); centrality",40,-0.1,0.1,2*nbinPtLP,pMin,ptMaxLP,100,0.,100.);
    fOutput->Add(fTriggerMCResEta);

    fTriggerMCResPhi = new TH3F("fTriggerMCResPhi","Trigger particle MC: #varphi resolution; #varphi_{MC}-#varphi_{Rec};#it{p}_{T} (GeV/#it{c}); centrality",40,-0.1,0.1,2*nbinPtLP,pMin,ptMaxLP,100,0.,100.);
    fOutput->Add(fTriggerMCResPhi);
  }

  // Pt Trigger particle vs centrality:
  fTriggerPtCent = new TH3F("fTriggerPtCent","Trigger particle;#it{p}_{T} (GeV/#it{c});centrality (%);Vertex Z (cm)",nbinPtLP,pMin,ptMaxLP,100,0.,100.,nbinsVtx,-10.,10.);
  fOutput->Add(fTriggerPtCent);

  fTriggerPtCentCh = new TH3F("fTriggerPtCentCh","Trigger particle;#it{p}_{T} (GeV/#it{c});centrality (%);Vertex Z (cm)",nbinPtLP,pMin,ptMaxLP,100,0.,100.,nbinsVtx,-10.,10.);
  fOutput->Add(fTriggerPtCentCh);

  fNTrigPerEvt = new TH2F("fNTrigPerEvt","Number of Trigger Particles Per Event;Counts;Centrality",51,-0.5,50.5,100,0.,100);
  fOutput->Add(fNTrigPerEvt);

  fTriggerWiSPDHit = new TH1F("fTriggerWiSPDHit","Number of Trigger Particles wi SPD Hits",3,0.,3.);
  fOutput->Add(fTriggerWiSPDHit);

  // Phi vs pseudorapidity:
  fTriggerEtaPhi = new TH2F("fTriggerEtaPhi","Trigger particle;#varphi (rad);#eta",nbinsPhi,0.,2.*TMath::Pi(),100,-1.,1.);
  fOutput->Add(fTriggerEtaPhi);
  
  // DCA to primary vertex:
  fTriggerDCA = new TH2F("fTriggerDCA","Trigger particle;DCA (cm);",32,0.,3.2,2,0.5,2.5);
  fTriggerDCA->GetYaxis()->SetBinLabel(1,"XY");
  fTriggerDCA->GetYaxis()->SetBinLabel(2,"Z");
  fOutput->Add(fTriggerDCA);

  // Fraction of shared TPC cls
  fTrigFracShTPCcls =
    new TH2F("fTrigFracShTPCcls","Trigger particle: vs #it{p}_{T} vs fraction Shared TPC cls;#it{p}_{T} (GeV/#it{c});fraction Shared TPC cls",
             nbinPtLP,pMin,ptMaxLP,50,0,1.);
  fOutput->Add(fTrigFracShTPCcls);

  // Check if Trigger particle comes from a V0 daughter:
  fCheckTriggerFromV0Daug = 
    new TH1F("fCheckTriggerFromV0Daug","Trigger particle from a V0 daughter;;Counts",4,-0.5,3.5);
  fCheckTriggerFromV0Daug->GetXaxis()->SetTitle("Flag"); 
  fCheckTriggerFromV0Daug->GetXaxis()->SetBinLabel(1,"NOT V0 daug");
  fCheckTriggerFromV0Daug->GetXaxis()->SetBinLabel(2,"V0 daug");
  fCheckTriggerFromV0Daug->GetXaxis()->SetBinLabel(3,"V0 daug & V0 LP");
  fOutput->Add(fCheckTriggerFromV0Daug);
  
  fTriggerComingFromDaug = new TH1F("fTriggerComingFromDaug","Trigger particle from a V0 daughter;#it{p}_{T} (GeV/#it{c});Counts",240, 0, 12);
  fOutput->Add(fTriggerComingFromDaug);

  fTriggerIsV0 = new TH1F("fTriggerIsV0","V0 candidate is a LP;#it{p}_{T} (GeV/#it{c});Counts",nbinPtLP,pMin,ptMaxLP);
  fOutput->Add(fTriggerIsV0);

  // ------------------- > Comparing properties of this trigger with the daughters
  //   K0s
  fCheckIDTrigPtK0s = new TH3F("fCheckIDTrigPtK0s","K^{0}_{S};#deltap/p_{tri};;p_{V0}",120,-0.6,0.6,7,-0.5,6.5,100,1.,11.);
  fCheckIDTrigPtK0s->GetYaxis()->SetBinLabel(1,"Pos Daug X");
  fCheckIDTrigPtK0s->GetYaxis()->SetBinLabel(2,"Pos Daug Y");
  fCheckIDTrigPtK0s->GetYaxis()->SetBinLabel(3,"Pos Daug Z");
  fCheckIDTrigPtK0s->GetYaxis()->SetBinLabel(4,"Neg Daug X");
  fCheckIDTrigPtK0s->GetYaxis()->SetBinLabel(5,"Neg Daug Y");
  fCheckIDTrigPtK0s->GetYaxis()->SetBinLabel(6,"Neg Daug Z");
  fOutput->Add(fCheckIDTrigPtK0s);

  fCheckIDTrigPhiK0s = new TH3F("fCheckIDTrigPhiK0s","K^{0}_{S};#delta#varphi;;p_{V0}",150,-TMath::Pi(),TMath::Pi(),3,-0.5,2.5,100,1.,11.);
  fCheckIDTrigPhiK0s->GetYaxis()->SetBinLabel(1,"Pos Daug");
  fCheckIDTrigPhiK0s->GetYaxis()->SetBinLabel(2,"Neg Daug");
  fOutput->Add(fCheckIDTrigPhiK0s);

  fCheckIDTrigEtaK0s = new TH3F("fCheckIDTrigEtaK0s","K^{0}_{S};#delta#eta;;p_{V0}",200,-2.,2.,3,-0.5,2.5,100,1.,11.);
  fCheckIDTrigEtaK0s->GetYaxis()->SetBinLabel(1,"Pos Daug");
  fCheckIDTrigEtaK0s->GetYaxis()->SetBinLabel(2,"Neg Daug");
  fOutput->Add(fCheckIDTrigEtaK0s);

  fCheckIDTrigNclsK0s = new TH3F("fCheckIDTrigNclsK0s","K^{0}_{S};NCls TPC;;p_{V0}",181,0.5,180.5,3,-0.5,2.5,100,1.,11.);
  fCheckIDTrigNclsK0s->GetYaxis()->SetBinLabel(1,"Pos Daug");
  fCheckIDTrigNclsK0s->GetYaxis()->SetBinLabel(2,"Neg Daug");
  fOutput->Add(fCheckIDTrigNclsK0s);

  //   Lambda
  fCheckIDTrigPtLambda = new TH3F("fCheckIDTrigPtLambda","#Lambda",120,-0.6,0.6,7,-0.5,6.5,100,1.,11.);
  fCheckIDTrigPtLambda->GetYaxis()->SetBinLabel(1,"Pos Daug X");
  fCheckIDTrigPtLambda->GetYaxis()->SetBinLabel(2,"Pos Daug Y");
  fCheckIDTrigPtLambda->GetYaxis()->SetBinLabel(3,"Pos Daug Z");
  fCheckIDTrigPtLambda->GetYaxis()->SetBinLabel(4,"Neg Daug X");
  fCheckIDTrigPtLambda->GetYaxis()->SetBinLabel(5,"Neg Daug Y");
  fCheckIDTrigPtLambda->GetYaxis()->SetBinLabel(6,"Neg Daug Z");
  fOutput->Add(fCheckIDTrigPtLambda);

  fCheckIDTrigPhiLambda  = new TH3F("fCheckIDTrigPhiLambda","#Lambda",150,-TMath::Pi(),TMath::Pi(),3,-0.5,2.5,100,1.,11.);
  fCheckIDTrigPhiLambda->GetYaxis()->SetBinLabel(1,"Pos Daug");
  fCheckIDTrigPhiLambda->GetYaxis()->SetBinLabel(2,"Neg Daug");
  fOutput->Add(fCheckIDTrigPhiLambda);

  fCheckIDTrigEtaLambda  = new TH3F("fCheckIDTrigEtaLambda","#Lambda",200,-2.,2.,3,-0.5,2.5,100,1.,11.);
  fCheckIDTrigEtaLambda->GetYaxis()->SetBinLabel(1,"Pos Daug");
  fCheckIDTrigEtaLambda->GetYaxis()->SetBinLabel(2,"Neg Daug");
  fOutput->Add(fCheckIDTrigEtaLambda);

  fCheckIDTrigNclsLambda = new TH3F("fCheckIDTrigNclsLambda","#Lambda;NCls TPC;;p_{V0}",181,0.5,180.5,3,-0.5,2.5,100,1.,11.);
  fCheckIDTrigNclsLambda->GetYaxis()->SetBinLabel(1,"Pos Daug");
  fCheckIDTrigNclsLambda->GetYaxis()->SetBinLabel(2,"Neg Daug");
  fOutput->Add(fCheckIDTrigNclsLambda);

  //   AntiLambda
  fCheckIDTrigPtAntiLambda = new TH3F("fCheckIDTrigPtAntiLambda","#bar{#Lambda}",120,-0.8,0.8,7,-0.5,6.5,100,1.,11.);
  fCheckIDTrigPtAntiLambda->GetYaxis()->SetBinLabel(1,"Pos Daug X");
  fCheckIDTrigPtAntiLambda->GetYaxis()->SetBinLabel(2,"Pos Daug Y");
  fCheckIDTrigPtAntiLambda->GetYaxis()->SetBinLabel(3,"Pos Daug Z");
  fCheckIDTrigPtAntiLambda->GetYaxis()->SetBinLabel(4,"Neg Daug X");
  fCheckIDTrigPtAntiLambda->GetYaxis()->SetBinLabel(5,"Neg Daug Y");
  fCheckIDTrigPtAntiLambda->GetYaxis()->SetBinLabel(6,"Neg Daug Z");
  fOutput->Add(fCheckIDTrigPtAntiLambda);

  fCheckIDTrigPhiAntiLambda  = new TH3F("fCheckIDTrigPhiAntiLambda","#bar{#Lambda}",150,-TMath::Pi(),TMath::Pi(),3,-0.5,2.5,100,1.,11.);
  fCheckIDTrigPhiAntiLambda->GetYaxis()->SetBinLabel(1,"Pos Daug");
  fCheckIDTrigPhiAntiLambda->GetYaxis()->SetBinLabel(2,"Neg Daug");
  fOutput->Add(fCheckIDTrigPhiAntiLambda);

  fCheckIDTrigEtaAntiLambda  = new TH3F("fCheckIDTrigEtaAntiLambda","#bar{#Lambda}",200,-2.,2.,3,-0.5,2.5,100,1.,11.);
  fCheckIDTrigEtaAntiLambda->GetYaxis()->SetBinLabel(1,"Pos Daug");
  fCheckIDTrigEtaAntiLambda->GetYaxis()->SetBinLabel(2,"Neg Daug");
  fOutput->Add(fCheckIDTrigEtaAntiLambda);

  fCheckIDTrigNclsAntiLambda = new TH3F("fCheckIDTrigNclsAntiLambda","#bar{#Lambda};NCls TPC;;p_{V0}",181,0.5,180.5,3,-0.5,2.5,100,1.,11.);
  fCheckIDTrigNclsAntiLambda->GetYaxis()->SetBinLabel(1,"Pos Daug");
  fCheckIDTrigNclsAntiLambda->GetYaxis()->SetBinLabel(2,"Neg Daug");
  fOutput->Add(fCheckIDTrigNclsAntiLambda);

  // ====== MC-true and  MC-Association information ====== //
  if(fIsMC){

    fInjectedParticles = new TH1F("fInjectedParticles","Injected particles;;Counts",2,0.,2.);
    fInjectedParticles->GetXaxis()->SetBinLabel(1,"Injected");
    fInjectedParticles->GetXaxis()->SetBinLabel(2,"Natural");
    fOutput->Add(fInjectedParticles);
    
    // K0s MC-true:
    fK0sMCPt       = new TH1F("fK0sMCPt", "K^{0}_{S} MC;#it{p}_{T} (GeV/#it{c});Counts",nbins,pMin,pMax);
    fOutput->Add(fK0sMCPt);

    fK0sMCPtRap    = new TH3F("fK0sMCPtRap", "K^{0}_{S} MC;#it{p}_{T} (GeV/#it{c});y;centrality",nbins,pMin,pMax,20,-1.0,1.0,100,0.,100.);
    fOutput->Add(fK0sMCPtRap);

    fK0sMCPtRap2   = new TH3F("fK0sMCPtRap2", "K^{0}_{S} MC;#it{p}_{T} (GeV/#it{c});y;centrality",nbins,pMin,pMax,20,-1.0,1.0,100,0.,100.);
    fOutput->Add(fK0sMCPtRap2);

    for(Int_t jj=0;jj<kNCent;jj++){
      snprintf(hNameHist,100, "fK0sMCPtRapVtx_Cent_%d",jj);
      fK0sMCPtRapVtx[jj] = new TH3F(hNameHist, "K^{0}_{S} MC  |VtxZ|;#it{p}_{T} (GeV/#it{c});y;VtxZ",nbins,pMin,pMax,20,-1.0,1.0,20,-10.,10.);
      fOutput->Add(fK0sMCPtRapVtx[jj]);
    }

    fK0sMCPtRapEmbeded   = new TH3F("fK0sMCPtRapEmbeded", "K^{0}_{S} Embeded MC;#it{p}_{T} (GeV/#it{c});y;centrality",nbins,pMin,pMax,20,-1.,1.,100,0.,100.);
    fOutput->Add(fK0sMCPtRapEmbeded);

    for(Int_t jj=0;jj<kNCent;jj++){
      snprintf(hNameHist,100, "fK0sMCPtRapVtxEmbeded_Cent_%d",jj);
      fK0sMCPtRapVtxEmbeded[jj] = new TH3F(hNameHist, "K^{0}_{S} Embeded MC |VtxZ|;#it{p}_{T} (GeV/#it{c});y;VtxZ",nbins,pMin,pMax,20,-1.0,1.0,20,-10.,10.);
      fOutput->Add(fK0sMCPtRapVtxEmbeded[jj]);
    }
  

    Int_t binsK0sMC[4] = {nbins,20,100,100};   Double_t xminK0sMC[4] = {pMin,-1.0,0.,0.};   Double_t xmaxK0sMC[4] = {pMax,1.0,10.,10.}; // gral efficiency
    for(Int_t jj=0;jj<kNCent;jj++){
       snprintf(hNameHist,100, "fK0sMCPtRapPtDaugPt_Cent_%d",jj);
       fK0sMCPtRapPtDaugPt[jj] = new THnSparseD(hNameHist, "K^{0}_{S} MC;#it{p}_{T} (GeV/#it{c});y;#it{p}_{T,Pos Daug} (GeV/#it{c});#it{p}_{T,Neg Daug} (GeV/#it{c});",
						4,binsK0sMC,xminK0sMC,xmaxK0sMC);
       fOutput->Add(fK0sMCPtRapPtDaugPt[jj]);


       snprintf(hNameHist,100, "fK0sMCPtRapPtDaugPtEmbeded_Cent_%d",jj);
       fK0sMCPtRapPtDaugPtEmbeded[jj] = new THnSparseD(hNameHist, "K^{0}_{S} MC Embeded;#it{p}_{T} (GeV/#it{c});y;#it{p}_{T,Pos Daug} (GeV/#it{c});#it{p}_{T,Neg Daug} (GeV/#it{c});",
						       4,binsK0sMC,xminK0sMC,xmaxK0sMC);
       fOutput->Add(fK0sMCPtRapPtDaugPtEmbeded[jj]);
    }


    for(Int_t jj=0;jj<kNCent;jj++){
      snprintf(hNameHist,100, "fK0sMCPtPhiEta_Cent_%d",jj);
      fK0sMCPtPhiEta[jj]    = new TH3F(hNameHist, "K^{0}_{S} MC;#varphi (rad);#eta;#it{p}_{T} (GeV/#it{c})",nbinsPhi,0.,2.*TMath::Pi(),20,-1.,1.,nbins,pMin,pMax);
      fOutput->Add(fK0sMCPtPhiEta[jj]);
    }
  
    // K0s MC-Association:
    fK0sAssocPt = 
      new TH1F("fK0sAssocPt","K^{0}_{S} Assoc: L_{T} vs #it{p}_{T};#it{p}_{T} (GeV/#it{c});Counts",nbins,pMin,pMax);
    fOutput->Add(fK0sAssocPt);

    fK0sAssocPtArm = 
      new TH3F("fK0sAssocPtArm","K^{0}_{S} Assoc: #it{p}_{T} vs y vs centrality;#it{p}_{T} (GeV/#it{c});y;centrality",nbins,pMin,pMax,20,-1.0,1.0,100,0.,100.);
    fOutput->Add(fK0sAssocPtArm);

    fK0sAssocPtRap    = new TH3F("fK0sAssocPtRap","K^{0}_{S} Assoc;#it{p}_{T} (GeV/#it{c});y;centrality",nbins,pMin,pMax,20,-1.0,1.0,100,0.,100.);
    fOutput->Add(fK0sAssocPtRap);

    fK0sAssocPtRapEmbeded    = new TH3F("fK0sAssocPtRapEmbeded","K^{0}_{S} Assoc  - Embeded MC;#it{p}_{T} (GeV/#it{c});y;centrality",nbins,pMin,pMax,20,-1.0,1.0,100,0.,100.);
    fOutput->Add(fK0sAssocPtRapEmbeded);

    for(Int_t jj=0;jj<kNCent;jj++){
      snprintf(hNameHist,100, "fK0sAssocPtPhiEta_Cent_%d",jj);
      fK0sAssocPtPhiEta[jj]    = new TH3F(hNameHist,"K^{0}_{S} Assoc;#varphi;#eta;#it{p}_{T} (GeV/#it{c})",nbinsPhi,0.,2.*TMath::Pi(),20,-1.0,1.0,nbins,pMin,pMax);
      fOutput->Add(fK0sAssocPtPhiEta[jj]);
    }


    // Histogramas para estudios sistematicos de la eficiencia
    Int_t binsEff1[3] = {nbins,nbins,20};              Double_t xminEff1[3] = {0.398,pMin,-1.0};             Double_t xmaxEff1[3] = {0.598,pMax,1.0};             // gral efficiency
    Int_t binsEff2[4] = {nbins,nbins,20,20};           Double_t xminEff2[4] = {0.398,pMin,-1.0,-10.};         Double_t xmaxEff2[4] = {0.598,pMax,1.0,10.};         // vtx cut
    Int_t binsEff3[4] = {nbins,nbins,20,60};           Double_t xminEff3[4] = {0.398,pMin,-1.0,0.};           Double_t xmaxEff3[4] = {0.598,pMax,1.0,1.2};         // dca between daughters
    Int_t binsEff4[4] = {nbins,nbins,20,50};           Double_t xminEff4[4] = {0.398,pMin,-1.0,0.9975};       Double_t xmaxEff4[4] = {0.598,pMax,1.0,1.};          // cpa
    Int_t binsEff5[5] = {nbins,nbins,20,99,99};        Double_t xminEff5[5] = {0.398,pMin,-1.0,0.,0.};        Double_t xmaxEff5[5] = {0.598,pMax,1.0,3.3,3.3};     // dca to prim. vtx
    Int_t binsEff6[5] = {nbins,nbins,20,170,170};        Double_t xminEff6[5] = {0.398,pMin,-1.0,0.5,0.5};      Double_t xmaxEff6[5] = {0.598,pMax,1.0,170.5,170};   // No. TPC Cls
    Int_t binsEffKsh[5] = {nbins,nbins,20,50,50};        Double_t xminEffKsh[5] = {0.398,pMin,-1.0,0.,0.};      Double_t xmaxEffKsh[5] = {0.598,pMax,1.0,1.0,1.0};   //shared TPC cls
    Int_t binsEffKPtDaug[5] = {nbins,nbins,20,100,100};  Double_t xminEffKPtDaug[5] = {0.398,pMin,-1.0,0.,0.};  Double_t xmaxEffKPtDaug[5] = {0.598,pMax,1.0,10.,10.};     //PtDaug
    Int_t binsEffKCtau[4] = {nbins,nbins,20,100};         Double_t xminEffKCtau[4] = {0.398,pMin,-1.0,0.};       Double_t xmaxEffKCtau[4] = {0.598,pMax,1.0,50.0};          //CTau
    Int_t binsEffKFidVol[4] = {nbins,nbins,20,100};      Double_t xminEffKFidVol[4] = {0.398,pMin,-1.0,0.};     Double_t xmaxEffKFidVol[4] = {0.598,pMax,1.0,100.};        //Fiducial volume

    for(Int_t i=0; i<kNCent; i++){
     
      /// ------- Natural particles
      snprintf(hNameHist,100, "fK0sAssocPtMassArm_Cent_%d",i);
      fK0sAssocPtMassArm[i]    = new THnSparseD(hNameHist,"K^{0}_{S} Assoc;Mass (GeV/c^{2});#it{p}_{T} (GeV/#it{c});rap;",3,binsEff1,xminEff1,xmaxEff1);
      fOutput->Add(fK0sAssocPtMassArm[i]);

      snprintf(hNameHist,100, "fK0sAssocMassPtVtx_Cent_%d",i);
      fK0sAssocMassPtVtx[i]  = new THnSparseD(hNameHist, "K^{0}_{S}; Mass (GeV/c^{2}); #it{p}_{T}; rap; VtxZ;",4,binsEff2,xminEff2,xmaxEff2);
      fOutput->Add(fK0sAssocMassPtVtx[i]);      

      snprintf(hNameHist,100, "fK0sAssocMassPtDCADaug_Cent_%d",i);
      fK0sAssocMassPtDCADaug[i]  = new THnSparseD(hNameHist, "K^{0}_{S}; Mass (GeV/c^{2}); #it{p}_{T}; rap; DCADaug;",4,binsEff3,xminEff3,xmaxEff3);
      fOutput->Add(fK0sAssocMassPtDCADaug[i]); 

      snprintf(hNameHist,100, "fK0sAssocMassPtCPA_Cent_%d",i);
      fK0sAssocMassPtCPA[i]  = new THnSparseD(hNameHist, "K^{0}_{S}; Mass (GeV/c^{2}); #it{p}_{T}; rap; CPA;",4,binsEff4,xminEff4,xmaxEff4);
      fOutput->Add(fK0sAssocMassPtCPA[i]);  
      
      snprintf(hNameHist,100, "fK0sAssocMassPtDCAPV_Cent_%d",i);
      fK0sAssocMassPtDCAPV[i]  = new THnSparseD(hNameHist, "K^{0}_{S}; Mass (GeV/c^{2}); #it{p}_{T}; rap; Pos DCA to Prim. Vtx; Neg DCA to Prim. Vtx;",5,binsEff5,xminEff5,xmaxEff5);
      fOutput->Add(fK0sAssocMassPtDCAPV[i]);  
     
      snprintf(hNameHist,100, "fK0sAssocMassPtDaugNClsTPC_Cent_%d",i);
      fK0sAssocMassPtDaugNClsTPC[i]  = new THnSparseD(hNameHist, "K^{0}_{S}; Mass (GeV/c^{2}); #it{p}_{T}; rap; Pos # TPC Cls; Neg # TPC Cls;",5,binsEff6,xminEff6,xmaxEff6);
      fOutput->Add(fK0sAssocMassPtDaugNClsTPC[i]); 

      snprintf(hNameHist,100, "fK0sAssocMassPtShTPCcls_Cent_%d",i);
      fK0sAssocMassPtShTPCcls[i]  = new THnSparseD(hNameHist, "K^{0}_{S}; Mass (GeV/c^{2}); #it{p}_{T}; rap; Pos fraction shared TPC Cls; Neg fraction shared TPC Cls;",5,binsEffKsh,xminEffKsh,xmaxEffKsh);
      fOutput->Add(fK0sAssocMassPtShTPCcls[i]); 

      snprintf(hNameHist,100, "fK0sAssocMassPtDaugPt_Cent_%d",i);
      fK0sAssocMassPtDaugPt[i]  = new THnSparseD(hNameHist, "K^{0}_{S}; Mass (GeV/c^{2}); #it{p}_{T}; rap; #it{p}_{T,Pos Daug} (GeV/#it{c}); #it{p}_{T,Neg Daug} (GeV/#it{c});",5,binsEffKPtDaug,xminEffKPtDaug,xmaxEffKPtDaug);
      fOutput->Add(fK0sAssocMassPtDaugPt[i]); 

      snprintf(hNameHist,100, "fK0sAssocMassPtCtau_Cent_%d",i);
      fK0sAssocMassPtCtau[i]  = new THnSparseD(hNameHist, "K^{0}_{S}; Mass (GeV/c^{2}); #it{p}_{T}; rap; c#tau;",4,binsEffKCtau,xminEffKCtau,xmaxEffKCtau);
      fOutput->Add(fK0sAssocMassPtCtau[i]); 

      snprintf(hNameHist,100, "fK0sAssocMassPtFidVolume_Cent_%d",i);
      fK0sAssocMassPtFidVolume[i]  = new THnSparseD(hNameHist, "K^{0}_{S}; Mass (GeV/c^{2}); #it{p}_{T}; rap; l_{T} (cm);",4,binsEffKFidVol,xminEffKFidVol,xmaxEffKFidVol);
      fOutput->Add(fK0sAssocMassPtFidVolume[i]); 


      /// ----- Embeded particles 
      snprintf(hNameHist,100, "fK0sAssocPtMassArmEmbeded_Cent_%d",i);
      fK0sAssocPtMassArmEmbeded[i]    = new THnSparseD(hNameHist,"K^{0}_{S} Assoc Embeded;Mass (GeV/c^{2});#it{p}_{T} (GeV/#it{c});rap;",3,binsEff1,xminEff1,xmaxEff1);
      fOutput->Add(fK0sAssocPtMassArmEmbeded[i]);

      snprintf(hNameHist,100, "fK0sAssocMassPtVtxEmbeded_Cent_%d",i);
      fK0sAssocMassPtVtxEmbeded[i]  = new THnSparseD(hNameHist, "K^{0}_{S} Embeded; Mass (GeV/c^{2}); #it{p}_{T}; rap; VtxZ;",4,binsEff2,xminEff2,xmaxEff2);
      fOutput->Add(fK0sAssocMassPtVtxEmbeded[i]);      

      snprintf(hNameHist,100, "fK0sAssocMassPtDCADaugEmbeded_Cent_%d",i);
      fK0sAssocMassPtDCADaugEmbeded[i]  = new THnSparseD(hNameHist, "K^{0}_{S}; Mass (GeV/c^{2}); #it{p}_{T}; rap; DCADaug;",4,binsEff3,xminEff3,xmaxEff3);
      fOutput->Add(fK0sAssocMassPtDCADaugEmbeded[i]); 

      snprintf(hNameHist,100, "fK0sAssocMassPtCPAEmbeded_Cent_%d",i);
      fK0sAssocMassPtCPAEmbeded[i]  = new THnSparseD(hNameHist, "K^{0}_{S}; Mass (GeV/c^{2}); #it{p}_{T}; rap; CPA;",4,binsEff4,xminEff4,xmaxEff4);
      fOutput->Add(fK0sAssocMassPtCPAEmbeded[i]);  

      snprintf(hNameHist,100, "fK0sAssocMassPtDCAPVEmbeded_Cent_%d",i);
      fK0sAssocMassPtDCAPVEmbeded[i]  = new THnSparseD(hNameHist, "K^{0}_{S}; Mass (GeV/c^{2}); #it{p}_{T}; rap; Pos DCA to Prim. Vtx; Neg DCA to Prim. Vtx;",5,binsEff5,xminEff5,xmaxEff5);
      fOutput->Add(fK0sAssocMassPtDCAPVEmbeded[i]);  

      snprintf(hNameHist,100, "fK0sAssocMassPtDaugNClsTPCEmbeded_Cent_%d",i);
      fK0sAssocMassPtDaugNClsTPCEmbeded[i]  = new THnSparseD(hNameHist, "K^{0}_{S}; Mass (GeV/c^{2}); #it{p}_{T}; rap; Pos # TPC Cls; Neg # TPC Cls;",5,binsEff6,xminEff6,xmaxEff6);
      fOutput->Add(fK0sAssocMassPtDaugNClsTPCEmbeded[i]); 

      snprintf(hNameHist,100, "fK0sAssocMassPtShTPCclsEmbeded_Cent_%d",i);
      fK0sAssocMassPtShTPCclsEmbeded[i]  = new THnSparseD(hNameHist, "K^{0}_{S}; Mass (GeV/c^{2}); #it{p}_{T}; rap; Pos fraction shared TPC Cls; Neg fraction shared TPC Cls;",5,binsEffKsh,xminEffKsh,xmaxEffKsh);
      fOutput->Add(fK0sAssocMassPtShTPCclsEmbeded[i]); 

      snprintf(hNameHist,100, "fK0sAssocMassPtDaugPtEmbeded_Cent_%d",i);
      fK0sAssocMassPtDaugPtEmbeded[i]  = new THnSparseD(hNameHist, "K^{0}_{S}; Mass (GeV/c^{2}); #it{p}_{T}; rap; #it{p}_{T,Pos Daug} (GeV/#it{c}); #it{p}_{T,Neg Daug} (GeV/#it{c});",5,binsEffKPtDaug,xminEffKPtDaug,xmaxEffKPtDaug);
      fOutput->Add(fK0sAssocMassPtDaugPtEmbeded[i]); 

      snprintf(hNameHist,100, "fK0sAssocMassPtCtauEmbeded_Cent_%d",i);
      fK0sAssocMassPtCtauEmbeded[i]  = new THnSparseD(hNameHist, "K^{0}_{S}; Mass (GeV/c^{2}); #it{p}_{T}; rap; c#tau;",4,binsEffKCtau,xminEffKCtau,xmaxEffKCtau);
      fOutput->Add(fK0sAssocMassPtCtauEmbeded[i]); 

      snprintf(hNameHist,100, "fK0sAssocMassPtFidVolumeEmbeded_Cent_%d",i);
      fK0sAssocMassPtFidVolumeEmbeded[i]  = new THnSparseD(hNameHist, "K^{0}_{S}; Mass (GeV/c^{2}); #it{p}_{T}; rap; l_{T} (cm);",4,binsEffKFidVol,xminEffKFidVol,xmaxEffKFidVol);
      fOutput->Add(fK0sAssocMassPtFidVolumeEmbeded[i]); 

    }
    
    fK0sMCResEta     = new TH3F("fK0sMCResEta","K^{0}_{S} Assoc: #eta resolution; #eta_{MC}-#eta_{Rec};#it{p}_{T} (GeV/#it{c}); centrality",40,-0.1,0.1,nbins,pMin,pMax,100,0.,100.);
    fOutput->Add(fK0sMCResEta);

    fK0sMCResPhi     = new TH3F("fK0sMCResPhi","K^{0}_{S} Assoc: #varphi resolution; #varphi_{MC}-#varphi_{Rec};#it{p}_{T} (GeV/#it{c}); centrality",40,-0.1,0.1,nbins,pMin,pMax,100,0.,100.);
    fOutput->Add(fK0sMCResPhi);

    fK0sMCResPt     = new TH3F("fK0sMCResPt","K^{0}_{S} Assoc: pt resolution; #it{p}_{T,MC}-#it{p]_{T,Rec};#it{p}_{T} (GeV/#it{c}); centrality",60,-0.3,0.3,nbins,pMin,pMax,100,0.,100.);
    fOutput->Add(fK0sMCResPt);

    fK0sPosMCResEta     = new TH3F("fK0sPosMCResEta","K^{0}_{S} Pos. Daug.: #eta resolution; #eta_{MC}-#eta_{Rec};#it{p}_{T} (GeV/#it{c}); centrality",40,-0.1,0.1,nbins,pMin,pMax,100,0.,100.);
    fOutput->Add(fK0sPosMCResEta);

    fK0sPosMCResPhi     = new TH3F("fK0sPosMCResPhi","K^{0}_{S}  Pos. Daug.: #varphi resolution; #varphi_{MC}-#varphi_{Rec};#it{p}_{T} (GeV/#it{c}); centrality",40,-0.1,0.1,nbins,pMin,pMax,100,0.,100.);
    fOutput->Add(fK0sPosMCResPhi);

    fK0sPosMCResPt     = new TH3F("fK0sPosMCResPt","K^{0}_{S}  Pos. Daug.: pt resolution; #it{p}_{T,MC}-#it{p]_{T,Rec};#it{p}_{T} (GeV/#it{c}); centrality",60,-0.3,0.3,nbins,pMin,pMax,100,0.,100.);
    fOutput->Add(fK0sPosMCResPt);  

    fK0sNegMCResEta     = new TH3F("fK0sNegMCResEta","K^{0}_{S} Neg. Daug.: #eta resolution; #eta_{MC}-#eta_{Rec};#it{p}_{T} (GeV/#it{c}); centrality",40,-0.1,0.1,nbins,pMin,pMax,100,0.,100.);
    fOutput->Add(fK0sNegMCResEta);

    fK0sNegMCResPhi     = new TH3F("fK0sNegMCResPhi","K^{0}_{S}  Neg. Daug.: #varphi resolution; #varphi_{MC}-#varphi_{Rec};#it{p}_{T} (GeV/#it{c}); centrality",40,-0.1,0.1,nbins,pMin,pMax,100,0.,100.);
    fOutput->Add(fK0sNegMCResPhi);

    fK0sNegMCResPt     = new TH3F("fK0sNegMCResPt","K^{0}_{S}  Neg. Daug.: pt resolution; #it{p}_{T,MC}-#it{p]_{T,Rec};#it{p}_{T} (GeV/#it{c}); centrality",60,-0.3,0.3,nbins,pMin,pMax,100,0.,100.);
    fOutput->Add(fK0sNegMCResPt);  

    // Lambda MC-true: 
    fLambdaMCPt = new TH1F("fLambdaMCPt","#Lambda MC;#it{p}_{T} (GeV/#it{c});Counts",nbins,pMin,pMax);
    fOutput->Add(fLambdaMCPt);

    fLambdaMCPtRap = new TH3F("fLambdaMCPtRap","#Lambda MC;#it{p}_{T} (GeV/#it{c});y;centrality",nbins,pMin,pMax,20,-1.0,1.0,100,0.,100.);
    fOutput->Add(fLambdaMCPtRap);

    fLambdaMCPtRap2 = new TH3F("fLambdaMCPtRap2","#Lambda MC;#it{p}_{T} (GeV/#it{c});y;centrality",nbins,pMin,pMax,20,-1.0,1.0,100,0.,100.);
    fOutput->Add(fLambdaMCPtRap2);

    for(Int_t jj=0;jj<kNCent;jj++){
      snprintf(hNameHist,100, "fLambdaMCPtRapVtx_Cent_%d",jj);
      fLambdaMCPtRapVtx[jj] = new TH3F(hNameHist,"#Lambda MC  |VtxZ|<3 cm;#it{p}_{T} (GeV/#it{c});y;zv",nbins,pMin,pMax,20,-1.0,1.0,20,-10.,10.);
      fOutput->Add(fLambdaMCPtRapVtx[jj]);
    }

    fLambdaMCPtRapEmbeded = new TH3F("fLambdaMCPtRapEmbeded","#Lambda Embeded MC;#it{p}_{T} (GeV/#it{c});y;centrality",nbins,pMin,pMax,20,-1.0,1.0,100,0.,100.);
    fOutput->Add(fLambdaMCPtRapEmbeded);
  
    for(Int_t jj=0;jj<kNCent;jj++){
      snprintf(hNameHist,100, "fLambdaMCPtRapVtxEmbeded_Cent_%d",jj);
      fLambdaMCPtRapVtxEmbeded[jj] = new TH3F(hNameHist,"#Lambda Embeded MC |VtxZ|<3 cm;#it{p}_{T} (GeV/#it{c});y;zv",nbins,pMin,pMax,20,-1.0,1.0,20,-10.,10.);
      fOutput->Add(fLambdaMCPtRapVtxEmbeded[jj]);
    }

    Int_t binsLambdaMC[4] = {nbins,20,100,100};   Double_t xminLambdaMC[4] = {pMin,-1.0,0.,0.};   Double_t xmaxLambdaMC[4] = {pMax,1.0,10.,10.}; // gral efficiency
    for(Int_t jj=0;jj<kNCent;jj++){
       snprintf(hNameHist,100, "fLambdaMCPtRapPtDaugPt_Cent_%d",jj);
       fLambdaMCPtRapPtDaugPt[jj] = new THnSparseD(hNameHist, "#Lambda MC;#it{p}_{T} (GeV/#it{c});y;#it{p}_{T,Pos Daug} (GeV/#it{c});#it{p}_{T,Neg Daug} (GeV/#it{c});",
						   4,binsLambdaMC,xminLambdaMC,xmaxLambdaMC);
       fOutput->Add(fLambdaMCPtRapPtDaugPt[jj]);


       snprintf(hNameHist,100, "fLambdaMCPtRapPtDaugPtEmbeded_Cent_%d",jj);
       fLambdaMCPtRapPtDaugPtEmbeded[jj] = new THnSparseD(hNameHist, "#Lambda MC Embeded;#it{p}_{T} (GeV/#it{c});y;#it{p}_{T,Pos Daug} (GeV/#it{c});#it{p}_{T,Neg Daug} (GeV/#it{c});",
							  4,binsLambdaMC,xminLambdaMC,xmaxLambdaMC);
       fOutput->Add(fLambdaMCPtRapPtDaugPtEmbeded[jj]);
    }

    fLambdaMCFromXi  = new TH2F("fLambdaMCFromXi", "#Lambda from Xi MC;#it{p}_{T} (GeV/#it{c});centrality",nbins,pMin,pMax,100,0.,100.);
    fOutput->Add(fLambdaMCFromXi);

    for(Int_t jj=0;jj<kNCent;jj++){
      snprintf(hNameHist,100, "fLambdaMCPtPhiEta_Cent_%d",jj);
      fLambdaMCPtPhiEta[jj] = new TH3F(hNameHist,"#Lambda MC;#varphi (rad);#eta;#it{p}_{T} (GeV/#it{c})",nbinsPhi,0.,2.*TMath::Pi(),20,-1.0,1.0,nbins,pMin,pMax);
      fOutput->Add(fLambdaMCPtPhiEta[jj]);
    }

    // Lambda MC-Association:
    fLambdaAssocPt = 
      new TH1F("fLambdaAssocPt","#Lambda Assoc: L_{T} vs #it{p}_{T};#it{p}_{T} (GeV/#it{c});Counts",nbins,pMin,pMax);
    fOutput->Add(fLambdaAssocPt);

    fLambdaAssocPtRap = new TH3F("fLambdaAssocPtRap", "#Lambda Assoc;#it{p}_{T} (GeV/#it{c});y;centrality",nbins,pMin,pMax,20,-1.0,1.0,100,0.,100.);
    fOutput->Add(fLambdaAssocPtRap);
    
    fLambdaAssocFromXi  = new TH2F("fLambdaAssocFromXi", "#Lambda from Xi Assoc;#it{p}_{T} (GeV/#it{c});centrality",nbins,pMin,pMax,100,0.,100.);
    fOutput->Add(fLambdaAssocFromXi);

    for(Int_t jj=0;jj<kNCent;jj++){
      snprintf(hNameHist,100, "fLambdaAssocPtPhiEta_Cent_%d",jj);
      fLambdaAssocPtPhiEta[jj] = new TH3F(hNameHist, "#Lambda Assoc;#varphi (rad);#eta;#it{p}_{T} (GeV/#it{c})",nbinsPhi,0.,2.*TMath::Pi(),20,-1.0,1.0,nbins,pMin,pMax);
      fOutput->Add(fLambdaAssocPtPhiEta[jj]);
    }
    
    // Histogramas para estudios sistematicos de la eficiencia
    Int_t binsEff7[3] = {nbins,nbins,20};          Double_t xminEff7[3] = {1.065,pMin,-1.0};           Double_t xmaxEff7[3] = {1.165,pMax,1.0};               // gral efficiency
    Int_t binsEff8[4] = {nbins,nbins,20,20};       Double_t xminEff8[4] = {1.065,pMin,-1.0,-10.};      Double_t xmaxEff8[4] = {1.165,pMax,1.0,10.};            // vtx
    Int_t binsEff9[4] = {nbins,nbins,20,60};       Double_t xminEff9[4] = {1.065,pMin,-1.0,0.};        Double_t xmaxEff9[4] = {1.165,pMax,1.0,1.2};            // dca between daughters
    Int_t binsEff10[4] = {nbins,nbins,20,50};      Double_t xminEff10[4] = {1.065,pMin,-1.0,0.9975};   Double_t xmaxEff10[4] = {1.165,pMax,1.0,1.};            // cpa
    Int_t binsEff11[5] = {nbins,nbins,20,99,99};   Double_t xminEff11[5] = {1.065,pMin,-1.0,0.,0.};    Double_t xmaxEff11[5] = {1.165,pMax,1.0,3.3,3.3};       // dca to prim. vtx
    Int_t binsEff12[5] = {nbins,nbins,20,170,170}; Double_t xminEff12[5] = {1.065,pMin,-1.0,0.5,0.5};  Double_t xmaxEff12[5] = {1.165,pMax,1.0,170.5,170.5};   // No. TPC Cls
    Int_t binsEffLsh[5] = {nbins,nbins,20,50,50};  Double_t xminEffLsh[5] = {1.065,pMin,-1.0,0.,0.};   Double_t xmaxEffLsh[5] = {1.165,pMax,1.0,1.0,1.0};      // shared TPC cls
    Int_t binsEffLPtDaug[5] = {nbins,nbins,20,100,100};  Double_t xminEffLPtDaug[5] = {1.065,pMin,-1.0,0.,0.};  Double_t xmaxEffLPtDaug[5] = {1.165,pMax,1.0,10.,10.};     //PtDaug
    Int_t binsEffLCtau[4]   = {nbins,nbins,20,100};       Double_t xminEffLCtau[4]   = {1.065,pMin,-1.0,0.};     Double_t xmaxEffLCtau[4]   = {1.165,pMax,1.0,50.0};         //CTau
    Int_t binsEffLFidVol[4] = {nbins,nbins,20,100};      Double_t xminEffLFidVol[4] = {1.065,pMin,-1.0,0.};     Double_t xmaxEffLFidVol[4] = {1.165,pMax,1.0,100.};        //Fiducial volume


    for(Int_t i=0; i<kNCent; i++){

      // --------- Natural particles
      snprintf(hNameHist,100, "fLambdaAssocMassPtRap_Cent_%d",i);
      fLambdaAssocMassPtRap[i]  = new THnSparseD(hNameHist, "#Lambda; Mass (GeV/c^{2}); #it{p}_{T}; rap;",3,binsEff7,xminEff7,xmaxEff7);
      fOutput->Add(fLambdaAssocMassPtRap[i]);      

      snprintf(hNameHist,100, "fLambdaAssocMassPtRap2_Cent_%d",i);
      fLambdaAssocMassPtRap2[i]  = new THnSparseD(hNameHist, "#Lambda; Mass (GeV/c^{2}); #it{p}_{T}; rap;",3,binsEff7,xminEff7,xmaxEff7);
      fOutput->Add(fLambdaAssocMassPtRap2[i]);     
      
      snprintf(hNameHist,100, "fLambdaAssocMassPtVtx_Cent_%d",i);
      fLambdaAssocMassPtVtx[i]  = new THnSparseD(hNameHist, "#Lambda; Mass (GeV/c^{2}); #it{p}_{T}; rap; VtxZ;",4,binsEff8,xminEff8,xmaxEff8);
      fOutput->Add(fLambdaAssocMassPtVtx[i]);      
     
      snprintf(hNameHist,100, "fLambdaAssocMassPtDCADaug_Cent_%d",i);
      fLambdaAssocMassPtDCADaug[i]  = new THnSparseD(hNameHist, "#Lambda; Mass (GeV/c^{2}); #it{p}_{T}; rap; DCADaug;",4,binsEff9,xminEff9,xmaxEff9);
      fOutput->Add(fLambdaAssocMassPtDCADaug[i]); 
     
      snprintf(hNameHist,100, "fLambdaAssocMassPtCPA_Cent_%d",i);
      fLambdaAssocMassPtCPA[i]  = new THnSparseD(hNameHist, "#Lambda; Mass (GeV/c^{2}); #it{p}_{T}; rap; CPA;",4,binsEff10,xminEff10,xmaxEff10);
      fOutput->Add(fLambdaAssocMassPtCPA[i]);  
    
      snprintf(hNameHist,100, "fLambdaAssocMassPtDCAPV_Cent_%d",i);
      fLambdaAssocMassPtDCAPV[i]  = new THnSparseD(hNameHist, "#Lambda; Mass (GeV/c^{2}); #it{p}_{T}; rap; Pos DCA to Prim. Vtx; Neg DCA to Prim. Vtx;",5,binsEff11,xminEff11,xmaxEff11);
      fOutput->Add(fLambdaAssocMassPtDCAPV[i]);  

      snprintf(hNameHist,100, "fLambdaAssocMassPtDaugNClsTPC_Cent_%d",i);
      fLambdaAssocMassPtDaugNClsTPC[i]  = new THnSparseD(hNameHist, "#Lambda; Mass (GeV/c^{2}); #it{p}_{T}; rap; Pos # TPC Cls; Neg # TPC Cls;",5,binsEff12,xminEff12,xmaxEff12);
      fOutput->Add(fLambdaAssocMassPtDaugNClsTPC[i]); 

      snprintf(hNameHist,100, "fLambdaAssocMassPtShTPCcls_Cent_%d",i);
      fLambdaAssocMassPtShTPCcls[i]  = new THnSparseD(hNameHist, "#Lambda; Mass (GeV/c^{2}); #it{p}_{T}; rap; Pos fraction shared TPC Cls; Neg fraction shared TPC Cls;",5,binsEffLsh,xminEffLsh,xmaxEffLsh);
      fOutput->Add(fLambdaAssocMassPtShTPCcls[i]); 

      snprintf(hNameHist,100, "fLambdaAssocMassPtDaugPt_Cent_%d",i);
      fLambdaAssocMassPtDaugPt[i]  = new THnSparseD(hNameHist, "#Lambda; Mass (GeV/c^{2}); #it{p}_{T}; rap; #it{p}_{T,Pos Daug} (GeV/#it{c}); #it{p}_{T,Neg Daug} (GeV/#it{c});",5,binsEffLPtDaug,xminEffLPtDaug,xmaxEffLPtDaug);
      fOutput->Add(fLambdaAssocMassPtDaugPt[i]); 

      snprintf(hNameHist,100, "fLambdaAssocMassPtCtau_Cent_%d",i);
      fLambdaAssocMassPtCtau[i]  = new THnSparseD(hNameHist, "#Lambda; Mass (GeV/c^{2}); #it{p}_{T}; rap; c#tau;",4,binsEffLCtau,xminEffLCtau,xmaxEffLCtau);
      fOutput->Add(fLambdaAssocMassPtCtau[i]); 

      snprintf(hNameHist,100, "fLambdaAssocMassPtFidVolume_Cent_%d",i);
      fLambdaAssocMassPtFidVolume[i]  = new THnSparseD(hNameHist, "#Lambda; Mass (GeV/c^{2}); #it{p}_{T}; rap; l_{T} (cm);",4,binsEffLFidVol,xminEffLFidVol,xmaxEffLFidVol);
      fOutput->Add(fLambdaAssocMassPtFidVolume[i]);  

      // ------------ Embeded particles
      snprintf(hNameHist,100, "fLambdaAssocMassPtRapEmbeded_Cent_%d",i);
      fLambdaAssocMassPtRapEmbeded[i]  = new THnSparseD(hNameHist, "#Lambda Embeded; Mass (GeV/c^{2}); #it{p}_{T}; rap;",3,binsEff7,xminEff7,xmaxEff7);
      fOutput->Add(fLambdaAssocMassPtRapEmbeded[i]);  

      snprintf(hNameHist,100, "fLambdaAssocMassPtRapEmbeded2_Cent_%d",i);
      fLambdaAssocMassPtRapEmbeded2[i]  = new THnSparseD(hNameHist, "#Lambda Embeded; Mass (GeV/c^{2}); #it{p}_{T}; rap;",3,binsEff7,xminEff7,xmaxEff7);
      fOutput->Add(fLambdaAssocMassPtRapEmbeded2[i]);    

      snprintf(hNameHist,100, "fLambdaAssocMassPtVtxEmbeded_Cent_%d",i);
      fLambdaAssocMassPtVtxEmbeded[i]  = new THnSparseD(hNameHist, "#Lambda Embeded; Mass (GeV/c^{2}); #it{p}_{T}; rap; VtxZ;",4,binsEff8,xminEff8,xmaxEff8);
      fOutput->Add(fLambdaAssocMassPtVtxEmbeded[i]);      

      snprintf(hNameHist,100, "fLambdaAssocMassPtDCADaugEmbeded_Cent_%d",i);
      fLambdaAssocMassPtDCADaugEmbeded[i]  = new THnSparseD(hNameHist, "#Lambda; Mass (GeV/c^{2}); #it{p}_{T}; rap; DCADaug;",4,binsEff9,xminEff9,xmaxEff9);
      fOutput->Add(fLambdaAssocMassPtDCADaugEmbeded[i]); 
 
      snprintf(hNameHist,100, "fLambdaAssocMassPtCPAEmbeded_Cent_%d",i);
      fLambdaAssocMassPtCPAEmbeded[i]  = new THnSparseD(hNameHist, "#Lambda; Mass (GeV/c^{2}); #it{p}_{T}; rap; CPA;",4,binsEff10,xminEff10,xmaxEff10);
      fOutput->Add(fLambdaAssocMassPtCPAEmbeded[i]);  

      snprintf(hNameHist,100, "fLambdaAssocMassPtDCAPVEmbeded_Cent_%d",i);
      fLambdaAssocMassPtDCAPVEmbeded[i]  = new THnSparseD(hNameHist, "#Lambda; Mass (GeV/c^{2}); #it{p}_{T}; rap; Pos DCA to Prim. Vtx; Neg DCA to Prim. Vtx;",5,binsEff11,xminEff11,xmaxEff11);
      fOutput->Add(fLambdaAssocMassPtDCAPVEmbeded[i]);  

      snprintf(hNameHist,100, "fLambdaAssocMassPtDaugNClsTPCEmbeded_Cent_%d",i);
      fLambdaAssocMassPtDaugNClsTPCEmbeded[i]  = new THnSparseD(hNameHist, "#Lambda; Mass (GeV/c^{2}); #it{p}_{T}; rap;  Pos # TPC Cls; Neg # TPC Cls;",5,binsEff12,xminEff12,xmaxEff12);
      fOutput->Add(fLambdaAssocMassPtDaugNClsTPCEmbeded[i]);

      snprintf(hNameHist,100, "fLambdaAssocMassPtShTPCclsEmbeded_Cent_%d",i);
      fLambdaAssocMassPtShTPCclsEmbeded[i]  = new THnSparseD(hNameHist, "#Lambda; Mass (GeV/c^{2}); #it{p}_{T}; rap; Pos fraction shared TPC Cls; Neg fraction shared TPC Cls;",5,binsEffLsh,xminEffLsh,xmaxEffLsh);
      fOutput->Add(fLambdaAssocMassPtShTPCclsEmbeded[i]); 

      snprintf(hNameHist,100, "fLambdaAssocMassPtDaugPtEmbeded_Cent_%d",i);
      fLambdaAssocMassPtDaugPtEmbeded[i]  = new THnSparseD(hNameHist, "#Lambda; Mass (GeV/c^{2}); #it{p}_{T}; rap; #it{p}_{T,Pos Daug} (GeV/#it{c}); #it{p}_{T,Neg Daug} (GeV/#it{c});",5,binsEffLPtDaug,xminEffLPtDaug,xmaxEffLPtDaug);
      fOutput->Add(fLambdaAssocMassPtDaugPtEmbeded[i]); 

      snprintf(hNameHist,100, "fLambdaAssocMassPtCtauEmbeded_Cent_%d",i);
      fLambdaAssocMassPtCtauEmbeded[i]  = new THnSparseD(hNameHist, "#Lambda; Mass (GeV/c^{2}); #it{p}_{T}; rap;  c#tau ;",4,binsEffLCtau,xminEffLCtau,xmaxEffLCtau);
      fOutput->Add(fLambdaAssocMassPtCtauEmbeded[i]); 

      snprintf(hNameHist,100, "fLambdaAssocMassPtFidVolumeEmbeded_Cent_%d",i);
      fLambdaAssocMassPtFidVolumeEmbeded[i]  = new THnSparseD(hNameHist, "#Lambda; Mass (GeV/c^{2}); #it{p}_{T}; rap; l_{T} (cm);",4,binsEffLFidVol,xminEffLFidVol,xmaxEffLFidVol);
      fOutput->Add(fLambdaAssocMassPtFidVolumeEmbeded[i]);  

    } 

    fLambdaMCResEta     = new TH3F("fLambdaMCResEta","#Lambda Assoc: #eta resolution; #eta_{MC}-#eta_{Rec};#it{p}_{T} (GeV/#it{c}); centrality",40,-0.1,0.1,nbins,pMin,pMax,100,0.,100.);
    fOutput->Add(fLambdaMCResEta);

    fLambdaMCResPhi     = new TH3F("fLambdaMCResPhi","#Lambda Assoc: #varphi resolution; #varphi_{MC}-#varphi_{Rec};#it{p}_{T} (GeV/#it{c}); centrality",40,-0.1,0.1,nbins,pMin,pMax,100,0.,100.);
    fOutput->Add(fLambdaMCResPhi);

    fLambdaMCResPt     = new TH3F("fLambdaMCResPt","#Lambda Assoc: pt resolution; #it{p}_{T,MC}-#it{p]_{T,Rec};#it{p}_{T} (GeV/#it{c}); centrality",60,-0.3,0.3,nbins,pMin,pMax,100,0.,100.);
    fOutput->Add(fLambdaMCResPt);

    fLambdaPosMCResEta     = new TH3F("fLambdaPosMCResEta","#Lambda Pos. Daug.: #eta resolution; #eta_{MC}-#eta_{Rec};#it{p}_{T} (GeV/#it{c}); centrality",40,-0.1,0.1,nbins,pMin,pMax,100,0.,100.);
    fOutput->Add(fLambdaPosMCResEta);

    fLambdaPosMCResPhi     = new TH3F("fLambdaPosMCResPhi","#Lambda  Pos. Daug.: #varphi resolution; #varphi_{MC}-#varphi_{Rec};#it{p}_{T} (GeV/#it{c}); centrality",40,-0.1,0.1,nbins,pMin,pMax,100,0.,100.);
    fOutput->Add(fLambdaPosMCResPhi);

    fLambdaPosMCResPt     = new TH3F("fLambdaPosMCResPt","#Lambda  Pos. Daug.: pt resolution; #it{p}_{T,MC}-#it{p]_{T,Rec};#it{p}_{T} (GeV/#it{c}); centrality",60,-0.3,0.3,nbins,pMin,pMax,100,0.,100.);
    fOutput->Add(fLambdaPosMCResPt);  

    fLambdaNegMCResEta     = new TH3F("fLambdaNegMCResEta","#Lambda Neg. Daug.: #eta resolution; #eta_{MC}-#eta_{Rec};#it{p}_{T} (GeV/#it{c}); centrality",40,-0.1,0.1,nbins,pMin,pMax,100,0.,100.);
    fOutput->Add(fLambdaNegMCResEta);

    fLambdaNegMCResPhi     = new TH3F("fLambdaNegMCResPhi","#Lambda  Neg. Daug.: #varphi resolution; #varphi_{MC}-#varphi_{Rec};#it{p}_{T} (GeV/#it{c}); centrality",40,-0.1,0.1,nbins,pMin,pMax,100,0.,100.);
    fOutput->Add(fLambdaNegMCResPhi);

    fLambdaNegMCResPt     = new TH3F("fLambdaNegMCResPt","#Lambda  Neg. Daug.: pt resolution; #it{p}_{T,MC}-#it{p]_{T,Rec};#it{p}_{T} (GeV/#it{c}); centrality",60,-0.3,0.3,nbins,pMin,pMax,100,0.,100.);
    fOutput->Add(fLambdaNegMCResPt);  

    // AntiLambda MC-true: 
    fAntiLambdaMCPt = new TH1F("fAntiLambdaMCPt","#bar{#Lambda} MC;#it{p}_{T} (GeV/#it{c});Counts",nbins,pMin,pMax);
    fOutput->Add(fAntiLambdaMCPt);
  
    fAntiLambdaMCPtRap = new TH3F("fAntiLambdaMCPtRap","#bar{#Lambda} MC;#it{p}_{T} (GeV/#it{c});y;centrality",nbins,pMin,pMax,20,-1.0,1.0,100,0.,100.);
    fOutput->Add(fAntiLambdaMCPtRap);
  
    fAntiLambdaMCPtRap2 = new TH3F("fAntiLambdaMCPtRap2","#bar{#Lambda} MC;#it{p}_{T} (GeV/#it{c});y;centrality",nbins,pMin,pMax,20,-1.0,1.0,100,0.,100.);
    fOutput->Add(fAntiLambdaMCPtRap2);

    for(Int_t jj=0;jj<kNCent;jj++){
      snprintf(hNameHist,100, "fAntiLambdaMCPtRapVtx_Cent_%d",jj);
      fAntiLambdaMCPtRapVtx[jj] = new TH3F(hNameHist,"#bar{#Lambda} MC |VtxZ|;#it{p}_{T} (GeV/#it{c});y;zv",nbins,pMin,pMax,20,-1.0,1.0,20,-10.,10.);
      fOutput->Add(fAntiLambdaMCPtRapVtx[jj]);  
    }

    fAntiLambdaMCPtRapEmbeded = new TH3F("fAntiLambdaMCPtRapEmbeded","#bar{#Lambda} Embeded MC;#it{p}_{T} (GeV/#it{c});y;centrality",nbins,pMin,pMax,20,-1.0,1.0,100,0.,100.);
    fOutput->Add(fAntiLambdaMCPtRapEmbeded);
    
    for(Int_t jj=0;jj<kNCent;jj++){
      snprintf(hNameHist,100, "fAntiLambdaMCPtRapVtxEmbeded_Cent_%d",jj);
      fAntiLambdaMCPtRapVtxEmbeded[jj] = new TH3F(hNameHist,"#bar{#Lambda} Embeded MC |VtxZ|;#it{p}_{T} (GeV/#it{c});y;zv",nbins,pMin,pMax,20,-1.0,1.0,20,-10.,10.);
      fOutput->Add(fAntiLambdaMCPtRapVtxEmbeded[jj]); 
    }

    Int_t binsAntiLambdaMC[4] = {nbins,20,100,100};   Double_t xminAntiLambdaMC[4] = {pMin,-1.0,0.,0.};   Double_t xmaxAntiLambdaMC[4] = {pMax,1.0,10.,10.}; // gral efficiency
    for(Int_t jj=0;jj<kNCent;jj++){
       snprintf(hNameHist,100, "fAntiLambdaMCPtRapPtDaugPt_Cent_%d",jj);
       fAntiLambdaMCPtRapPtDaugPt[jj] = new THnSparseD(hNameHist, "#bar{#Lambda} MC;#it{p}_{T} (GeV/#it{c});y;#it{p}_{T,Pos Daug} (GeV/#it{c});#it{p}_{T,Neg Daug} (GeV/#it{c});",
						       4,binsAntiLambdaMC,xminAntiLambdaMC,xmaxAntiLambdaMC);
       fOutput->Add(fAntiLambdaMCPtRapPtDaugPt[jj]);


       snprintf(hNameHist,100, "fLambdaMCPtRapPtDaugPtEmbeded_Cent_%d",jj);
       fAntiLambdaMCPtRapPtDaugPtEmbeded[jj] = new THnSparseD(hNameHist, "#bar{#Lambda} MC Embeded;#it{p}_{T} (GeV/#it{c});y;#it{p}_{T,Pos Daug} (GeV/#it{c});#it{p}_{T,Neg Daug} (GeV/#it{c});",
							      4,binsAntiLambdaMC,xminAntiLambdaMC,xmaxAntiLambdaMC);
       fOutput->Add(fAntiLambdaMCPtRapPtDaugPtEmbeded[jj]);
    }

    fAntiLambdaMCFromXi  = new TH2F("fAntiLambdaMCFromXi", "#bar{#Lambda} from Xi MC;#it{p}_{T} (GeV/#it{c});centrality",nbins,pMin,pMax,100,0.,100.);
    fOutput->Add(fAntiLambdaMCFromXi);

    for(Int_t jj=0;jj<kNCent;jj++){
      snprintf(hNameHist,100, "fAntiLambdaMCPtPhiEta_Cent_%d",jj);
      fAntiLambdaMCPtPhiEta[jj] = new TH3F(hNameHist,"#bar{#Lambda} MC;#varphi (rad);#eta;#it{p}_{T} (GeV/#it{c})",nbinsPhi,0.,2.*TMath::Pi(),20,-1.0,1.0,nbins,pMin,pMax);
      fOutput->Add(fAntiLambdaMCPtPhiEta[jj]);
    }
  
    // AntiLambda MC-Association:
    fAntiLambdaAssocPt = 
      new TH1F("fAntiLambdaAssocPt","#bar{#Lambda} Assoc: L_{T} vs #it{p}_{T};#it{p}_{T} (GeV/#it{c})",nbins,pMin,pMax);
    fOutput->Add(fAntiLambdaAssocPt);
  
    fAntiLambdaAssocPtRap = new TH3F("fAntiLambdaAssocPtRap", "#bar{#Lambda} Assoc;#it{p}_{T} (GeV/#it{c});y;centrality",nbins,pMin,pMax,20,-1.0,1.0,100,0.,100.);
    fOutput->Add(fAntiLambdaAssocPtRap);
  
    fAntiLambdaAssocFromXi  = new TH2F("fAntiLambdaAssocFromXi", "#bar{#Lambda} from Xi MC;#it{p}_{T} (GeV/#it{c});centrality",nbins,pMin,pMax,100,0.,100.);
    fOutput->Add(fAntiLambdaAssocFromXi);

    for(Int_t jj=0;jj<kNCent;jj++){
      snprintf(hNameHist,100, "fAntiLambdaAssocPtPhiEta_Cent_%d",jj);
      fAntiLambdaAssocPtPhiEta[jj] = new TH3F(hNameHist, "#Lambda Assoc;#varphi (rad);#eta;#it{p}_{T} (GeV/#it{c})",nbinsPhi,0.,2.*TMath::Pi(),20,-1.0,1.0,nbins,pMin,pMax);
      fOutput->Add(fAntiLambdaAssocPtPhiEta[jj]);
    }

    // Histogramas para estudios sistematicos de la eficiencia
    Int_t binsEff13[3] = {nbins,nbins,20};          Double_t xminEff13[3] = {1.065,pMin,-1.0};          Double_t xmaxEff13[3] = {1.165,pMax,1.0};              // gral efficiency
    Int_t binsEff14[4] = {nbins,nbins,20,20};       Double_t xminEff14[4] = {1.065,pMin,-1.0,-10.};     Double_t xmaxEff14[4] = {1.165,pMax,1.0,10.};          // vtx
    Int_t binsEff15[4] = {nbins,nbins,20,60};       Double_t xminEff15[4] = {1.065,pMin,-1.0,0.};       Double_t xmaxEff15[4] = {1.165,pMax,1.0,1.2};          // dca between daug
    Int_t binsEff16[4] = {nbins,nbins,20,50};       Double_t xminEff16[4] = {1.065,pMin,-1.0,0.9975};   Double_t xmaxEff16[4] = {1.165,pMax,1.0,1.};           // cpa
    Int_t binsEff17[5] = {nbins,nbins,20,99,99};    Double_t xminEff17[5] = {1.065,pMin,-1.0,0.,0.};    Double_t xmaxEff17[5] = {1.165,pMax,1.0,3.3,3.3};      // dca to prim. vtx
    Int_t binsEff18[5] = {nbins,nbins,20,170,170};  Double_t xminEff18[5] = {1.065,pMin,-1.0,0.5,0.5};  Double_t xmaxEff18[5] = {1.165,pMax,1.0,170.5,170.5};  // No. TPC Cls
    Int_t binsEffALsh[5] = {nbins,nbins,20,50,50};  Double_t xminEffALsh[5] = {1.065,pMin,-1.0,0.,0.};  Double_t xmaxEffALsh[5] = {1.165,pMax,1.0,1.0,1.0};    // shared TPC cls
    Int_t binsEffALPtDaug[5] = {nbins,nbins,20,100,100};  Double_t xminEffALPtDaug[5] = {1.065,pMin,-1.0,0.,0.};  Double_t xmaxEffALPtDaug[5] = {1.065,pMax,1.0,10.,10.};     //PtDaug
    Int_t binsEffALCtau[4]   = {nbins,nbins,20,100};       Double_t xminEffALCtau[4]   = {1.065,pMin,-1.0,0.};     Double_t xmaxEffALCtau[4]   = {1.065,pMax,1.0,50.0};         //CTau 
    Int_t binsEffALFidVol[4] = {nbins,nbins,20,100};      Double_t xminEffALFidVol[4] = {1.065,pMin,-1.0,0.};     Double_t xmaxEffALFidVol[4] = {1.065,pMax,1.0,100.};        //Fiducial volume

    for(Int_t i=0; i<kNCent; i++){
      // --------- Natural particles
      snprintf(hNameHist,100, "fAntiLambdaAssocMassPtRap_Cent_%d",i);
      fAntiLambdaAssocMassPtRap[i]  = new THnSparseD(hNameHist, "#bar{#Lambda}; Mass (GeV/c^{2}); #it{p}_{T}; rap;",3,binsEff13,xminEff13,xmaxEff13);
      fOutput->Add(fAntiLambdaAssocMassPtRap[i]);      
  
      snprintf(hNameHist,100, "fAntiLambdaAssocMassPtRap2_Cent_%d",i);
      fAntiLambdaAssocMassPtRap2[i]  = new THnSparseD(hNameHist, "#bar{#Lambda}; Mass (GeV/c^{2}); #it{p}_{T}; rap;",3,binsEff13,xminEff13,xmaxEff13);
      fOutput->Add(fAntiLambdaAssocMassPtRap2[i]); 
      
      snprintf(hNameHist,100, "fAntiLambdaAssocMassPtVtx_Cent_%d",i);
      fAntiLambdaAssocMassPtVtx[i]  = new THnSparseD(hNameHist, "#bar{#Lambda}; Mass (GeV/c^{2}); #it{p}_{T}; rap; VtxZ;",4,binsEff14,xminEff14,xmaxEff14);
      fOutput->Add(fAntiLambdaAssocMassPtVtx[i]);      

      snprintf(hNameHist,100, "fAntiLambdaAssocMassPtDCADaug_Cent_%d",i);
      fAntiLambdaAssocMassPtDCADaug[i]  = new THnSparseD(hNameHist, "#bar{#Lambda}; Mass (GeV/c^{2}); #it{p}_{T}; rap; DCADaug;",4,binsEff15,xminEff15,xmaxEff15);
      fOutput->Add(fAntiLambdaAssocMassPtDCADaug[i]); 

      snprintf(hNameHist,100, "fAntiLambdaAssocMassPtCPA_Cent_%d",i);
      fAntiLambdaAssocMassPtCPA[i]  = new THnSparseD(hNameHist, "#bar{#Lambda}; Mass (GeV/c^{2}); #it{p}_{T}; rap; CPA;",4,binsEff16,xminEff16,xmaxEff16);
      fOutput->Add(fAntiLambdaAssocMassPtCPA[i]);  

      snprintf(hNameHist,100, "fAntiLambdaAssocMassPtDCAPV_Cent_%d",i);
      fAntiLambdaAssocMassPtDCAPV[i]  = new THnSparseD(hNameHist, "#bar{#Lambda}; Mass (GeV/c^{2}); #it{p}_{T}; rap; Pos DCA to Prim. Vtx; Neg DCA to Prim. Vtx;",5,binsEff17,xminEff17,xmaxEff17);
      fOutput->Add(fAntiLambdaAssocMassPtDCAPV[i]);  

      snprintf(hNameHist,100, "fAntiLambdaAssocMassPtDaugNClsTPC_Cent_%d",i);
      fAntiLambdaAssocMassPtDaugNClsTPC[i]  = new THnSparseD(hNameHist, "#bar{#Lambda}; Mass (GeV/c^{2}); #it{p}_{T}; rap;  Pos # TPC Cls; Neg # TPC Cls;",5,binsEff18,xminEff18,xmaxEff18);
      fOutput->Add(fAntiLambdaAssocMassPtDaugNClsTPC[i]); 

      snprintf(hNameHist,100, "fAntiLambdaAssocMassPtShTPCcls_Cent_%d",i);
      fAntiLambdaAssocMassPtShTPCcls[i]  = new THnSparseD(hNameHist, "#bar{#Lambda}; Mass (GeV/c^{2}); #it{p}_{T}; rap; Pos fraction shared TPC Cls; Neg fraction shared TPC Cls;",5,binsEffALsh,xminEffALsh,xmaxEffALsh);
      fOutput->Add(fAntiLambdaAssocMassPtShTPCcls[i]); 

      snprintf(hNameHist,100, "fAntiLambdaAssocMassPtDaugPt_Cent_%d",i);
      fAntiLambdaAssocMassPtDaugPt[i]  = new THnSparseD(hNameHist, "#bar{#Lambda}; Mass (GeV/c^{2}); #it{p}_{T}; rap; #it{p}_{T,Pos Daug} (GeV/#it{c}); #it{p}_{T,Neg Daug} (GeV/#it{c});",5,binsEffALPtDaug,xminEffALPtDaug,xmaxEffALPtDaug);
      fOutput->Add(fAntiLambdaAssocMassPtDaugPt[i]); 

      snprintf(hNameHist,100, "fAntiLambdaAssocMassPtCtau_Cent_%d",i);
      fAntiLambdaAssocMassPtCtau[i]  = new THnSparseD(hNameHist, "#bar{#Lambda}; Mass (GeV/c^{2}); #it{p}_{T}; rap; c#tau;",4,binsEffALCtau,xminEffALCtau,xmaxEffALCtau);
      fOutput->Add(fAntiLambdaAssocMassPtCtau[i]); 

      snprintf(hNameHist,100, "fAntiLambdaAssocMassPtFidVolume_Cent_%d",i);
      fAntiLambdaAssocMassPtFidVolume[i]  = new THnSparseD(hNameHist, "#bar{#Lambda}; Mass (GeV/c^{2}); #it{p}_{T}; rap; l_{T} (cm);",4,binsEffALFidVol,xminEffALFidVol,xmaxEffALFidVol);
      fOutput->Add(fAntiLambdaAssocMassPtFidVolume[i]);  

      // ------------ Embeded particles
      snprintf(hNameHist,100, "fAntiLambdaAssocMassPtRapEmbeded_Cent_%d",i);
      fAntiLambdaAssocMassPtRapEmbeded[i]  = new THnSparseD(hNameHist, "#bar{#Lambda} Embeded; Mass (GeV/c^{2}); #it{p}_{T}; rap;",3,binsEff13,xminEff13,xmaxEff13);
      fOutput->Add(fAntiLambdaAssocMassPtRapEmbeded[i]);    

      snprintf(hNameHist,100, "fAntiLambdaAssocMassPtRapEmbeded2_Cent_%d",i);
      fAntiLambdaAssocMassPtRapEmbeded2[i]  = new THnSparseD(hNameHist, "#bar{#Lambda} Embeded; Mass (GeV/c^{2}); #it{p}_{T}; rap;",3,binsEff13,xminEff13,xmaxEff13);
      fOutput->Add(fAntiLambdaAssocMassPtRapEmbeded2[i]);    

      snprintf(hNameHist,100, "fAntiLambdaAssocMassPtVtxEmbeded_Cent_%d",i);
      fAntiLambdaAssocMassPtVtxEmbeded[i]  = new THnSparseD(hNameHist, "#bar{#Lambda} Embeded; Mass (GeV/c^{2}); #it{p}_{T}; rap; VtxZ;",4,binsEff14,xminEff14,xmaxEff14);
      fOutput->Add(fAntiLambdaAssocMassPtVtxEmbeded[i]);      

      snprintf(hNameHist,100, "fAntiLambdaAssocMassPtDCADaugEmbeded_Cent_%d",i);
      fAntiLambdaAssocMassPtDCADaugEmbeded[i]  = new THnSparseD(hNameHist, "#bar{#Lambda}; Mass (GeV/c^{2}); #it{p}_{T}; rap; DCADaug;",4,binsEff15,xminEff15,xmaxEff15);
      fOutput->Add(fAntiLambdaAssocMassPtDCADaugEmbeded[i]); 

      snprintf(hNameHist,100, "fAntiLambdaAssocMassPtCPAEmbeded_Cent_%d",i);
      fAntiLambdaAssocMassPtCPAEmbeded[i]  = new THnSparseD(hNameHist, "#bar{#Lambda}; Mass (GeV/c^{2}); #it{p}_{T}; rap; CPA;",4,binsEff16,xminEff16,xmaxEff16);
      fOutput->Add(fAntiLambdaAssocMassPtCPAEmbeded[i]);  

      snprintf(hNameHist,100, "fAntiLambdaAssocMassPtDCAPVEmbeded_Cent_%d",i);
      fAntiLambdaAssocMassPtDCAPVEmbeded[i]  = new THnSparseD(hNameHist, "#bar{#Lambda}; Mass (GeV/c^{2}); #it{p}_{T}; rap; Pos DCA to Prim. Vtx; Neg DCA to Prim. Vtx;",5,binsEff17,xminEff17,xmaxEff17);
      fOutput->Add(fAntiLambdaAssocMassPtDCAPVEmbeded[i]);  

      snprintf(hNameHist,100, "fAntiLambdaAssocMassPtDaugNClsTPCEmbeded_Cent_%d",i);
      fAntiLambdaAssocMassPtDaugNClsTPCEmbeded[i]  = new THnSparseD(hNameHist, "#bar{#Lambda}; Mass (GeV/c^{2}); #it{p}_{T}; rap;  Pos # TPC Cls; Neg # TPC Cls;",5,binsEff18,xminEff18,xmaxEff18);
      fOutput->Add(fAntiLambdaAssocMassPtDaugNClsTPCEmbeded[i]);
  
      snprintf(hNameHist,100, "fAntiLambdaAssocMassPtShTPCclsEmbeded_Cent_%d",i);
      fAntiLambdaAssocMassPtShTPCclsEmbeded[i]  = new THnSparseD(hNameHist, "#bar{#Lambda}; Mass (GeV/c^{2}); #it{p}_{T}; rap; Pos fraction shared TPC Cls; Neg fraction shared TPC Cls;",5,binsEffALsh,xminEffALsh,xmaxEffALsh);
      fOutput->Add(fAntiLambdaAssocMassPtShTPCclsEmbeded[i]); 

      snprintf(hNameHist,100, "fAntiLambdaAssocMassPtDaugPtEmbeded_Cent_%d",i);
      fAntiLambdaAssocMassPtDaugPtEmbeded[i]  = new THnSparseD(hNameHist, "#bar{#Lambda}; Mass (GeV/c^{2}); #it{p}_{T}; rap; #it{p}_{T,Pos Daug} (GeV/#it{c}); #it{p}_{T,Neg Daug} (GeV/#it{c});",5,binsEffALPtDaug,xminEffALPtDaug,xmaxEffALPtDaug);
      fOutput->Add(fAntiLambdaAssocMassPtDaugPtEmbeded[i]); 

      snprintf(hNameHist,100, "fAntiLambdaAssocMassPtCtauEmbeded_Cent_%d",i);
      fAntiLambdaAssocMassPtCtauEmbeded[i]  = new THnSparseD(hNameHist, "#bar{#Lambda}; Mass (GeV/c^{2}); #it{p}_{T}; rap; c#tau;",4,binsEffALCtau,xminEffALCtau,xmaxEffALCtau);
      fOutput->Add(fAntiLambdaAssocMassPtCtauEmbeded[i]); 

      snprintf(hNameHist,100, "fAntiLambdaAssocMassPtFidVolumeEmbeded_Cent_%d",i);
      fAntiLambdaAssocMassPtFidVolumeEmbeded[i]  = new THnSparseD(hNameHist, "#bar{#Lambda}; Mass (GeV/c^{2}); #it{p}_{T}; rap; l_{T} (cm);",4,binsEffALFidVol,xminEffALFidVol,xmaxEffALFidVol);
      fOutput->Add(fAntiLambdaAssocMassPtFidVolumeEmbeded[i]);  


  } 

    fAntiLambdaMCResEta     = new TH3F("fAntiLambdaMCResEta","#bar{#Lambda} Assoc: #eta resolution; #eta_{MC}-#eta_{Rec};#it{p}_{T} (GeV/#it{c}); centrality",40,-0.1,0.1,nbins,pMin,pMax,100,0.,100.);
    fOutput->Add(fAntiLambdaMCResEta);

    fAntiLambdaMCResPhi     = new TH3F("fAntiLambdaMCResPhi","#bar{#Lambda} Assoc: #varphi resolution; #varphi_{MC}-#varphi_{Rec};#it{p}_{T} (GeV/#it{c}); centrality",40,-0.1,0.1,nbins,pMin,pMax,100,0.,100.);
    fOutput->Add(fAntiLambdaMCResPhi);

    fAntiLambdaMCResPt     = new TH3F("fAntiLambdaMCResPt","#bar{#Lambda} Assoc: pt resolution; #it{p}_{T,MC}-#it{p]_{T,Rec};#it{p}_{T} (GeV/#it{c}); centrality",60,-0.3,0.3,nbins,pMin,pMax,100,0.,100.);
    fOutput->Add(fAntiLambdaMCResPt);

    fAntiLambdaPosMCResEta     = new TH3F("fAntiLambdaPosMCResEta","#bar{#Lambda} Pos. Daug.: #eta resolution; #eta_{MC}-#eta_{Rec};#it{p}_{T} (GeV/#it{c}); centrality",40,-0.1,0.1,nbins,pMin,pMax,100,0.,100.);
    fOutput->Add(fAntiLambdaPosMCResEta);

    fAntiLambdaPosMCResPhi     = new TH3F("fAntiLambdaPosMCResPhi","#bar{#Lambda}  Pos. Daug.: #varphi resolution; #varphi_{MC}-#varphi_{Rec};#it{p}_{T} (GeV/#it{c}); centrality",40,-0.1,0.1,nbins,pMin,pMax,100,0.,100.);
    fOutput->Add(fAntiLambdaPosMCResPhi);

    fAntiLambdaPosMCResPt     = new TH3F("fAntiLambdaPosMCResPt","#bar{#Lambda}  Pos. Daug.: pt resolution; #it{p}_{T,MC}-#it{p]_{T,Rec};#it{p}_{T} (GeV/#it{c}); centrality",60,-0.3,0.3,nbins,pMin,pMax,100,0.,100.);
    fOutput->Add(fAntiLambdaPosMCResPt);  

    fAntiLambdaNegMCResEta     = new TH3F("fAntiLambdaNegMCResEta","#bar{#Lambda} Neg. Daug.: #eta resolution; #eta_{MC}-#eta_{Rec};#it{p}_{T} (GeV/#it{c}); centrality",40,-0.1,0.1,nbins,pMin,pMax,100,0.,100.);
    fOutput->Add(fAntiLambdaNegMCResEta);

    fAntiLambdaNegMCResPhi     = new TH3F("fAntiLambdaNegMCResPhi","#bar{#Lambda}  Neg. Daug.: #varphi resolution; #varphi_{MC}-#varphi_{Rec};#it{p}_{T} (GeV/#it{c}); centrality",40,-0.1,0.1,nbins,pMin,pMax,100,0.,100.);
    fOutput->Add(fAntiLambdaNegMCResPhi);

    fAntiLambdaNegMCResPt     = new TH3F("fAntiLambdaNegMCResPt","#bar{#Lambda}  Neg. Daug.: pt resolution; #it{p}_{T,MC}-#it{p]_{T,Rec};#it{p}_{T} (GeV/#it{c}); centrality",60,-0.3,0.3,nbins,pMin,pMax,100,0.,100.);
    fOutput->Add(fAntiLambdaNegMCResPt);  

  } //End MC

  // ======================================================== //
  // ========== Reconstruction information in AOD =========== //
  fHistArmenterosPodolanski  =
    new TH3F("fHistArmenterosPodolanski","Armenteros-Podolanski phase space;#alpha;p_{t} arm",
             100,-1.0,1.0,50,0,0.5,7,-0.5,6.5);
  fHistArmenterosPodolanski->GetZaxis()->SetBinLabel(1,"K^{0}_{S} Inv. Mass Peak");
  fHistArmenterosPodolanski->GetZaxis()->SetBinLabel(2,"K^{0}_{S} Bckg");
  fHistArmenterosPodolanski->GetZaxis()->SetBinLabel(3,"#Lambda Inv. Mass Peak");
  fHistArmenterosPodolanski->GetZaxis()->SetBinLabel(4,"#Lambda Bckg");
  fHistArmenterosPodolanski->GetZaxis()->SetBinLabel(5,"#bar{#Lambda} Inv. Mass Peak");
  fHistArmenterosPodolanski->GetZaxis()->SetBinLabel(6,"#bar{#Lambda} Bckg");
  fOutput->Add(fHistArmenterosPodolanski);
 
  fHistArmPodBckg =
    new TH3F("fHistArmPodBckg","Armenteros-Podolanski phase space for correlations;#alpha;p_{t} arm",
             100,-1.0,1.0,50,0,0.5,6,-0.5,5.5);
  fHistArmPodBckg->GetZaxis()->SetBinLabel(1,"K^{0}_{S} SigBck: Trig events");
  fHistArmPodBckg->GetZaxis()->SetBinLabel(2,"K^{0}_{S} Bck: Trig events");
  fHistArmPodBckg->GetZaxis()->SetBinLabel(3,"#Lambda SigBck: Trig events");
  fHistArmPodBckg->GetZaxis()->SetBinLabel(4,"#Lambda Bck: Trig events");
  fHistArmPodBckg->GetZaxis()->SetBinLabel(5,"#bar{#Lambda} SigBck: Trig events");
  fHistArmPodBckg->GetZaxis()->SetBinLabel(6,"#bar{#Lambda} Bck: Trig events");
  fOutput->Add(fHistArmPodBckg);
 
  // ****** K0s ******
  fK0sMass =
    new TH3F("fK0sMass", "K^{0}_{s}: mass vs #it{p}_{T};Mass (GeV/#it{c}^2);#it{p}_{T} (GeV/#it{c});centrality",nbins,0.398,0.598,nbins,pMin,pMax,100,0.,100.);

  fOutput->Add(fK0sMass);
 
  fK0sMassEmbeded =
    new TH3F("fK0sMassEmbeded", "K^{0}_{s} Embeded: mass vs #it{p}_{T}",nbins,0.398,0.598,nbins,pMin,pMax,100,0.,100.);
  fK0sMassEmbeded->GetXaxis()->SetTitle("Mass (GeV/c^2)");
  fK0sMassEmbeded->GetYaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fK0sMassEmbeded->GetZaxis()->SetTitle("centrality");
  fOutput->Add(fK0sMassEmbeded);


  fK0sMassPtEta =
    new TH3F("fK0sMassPtEta","K^{0}_{s}: Mass vs #it{p}_{T} vs #eta;Mass (GeV/C^{2});#it{p}_{T} (GeV/#it{c});#eta",
	     nbins,0.398,0.598,nbins,pMin,pMax,20,-1.0,1.0);
  fOutput->Add(fK0sMassPtEta);
 
  for(Int_t i=0; i<kNCent; i++){
    fK0sMassPtRap[i] =
      new TH3F(Form("fK0sMassPtRap_cent_%.0lf_%.0lf",kBinCent[i],kBinCent[i+1]),"K^{0}_{s}: mass vs #it{p}_{T} vs y;Mass (GeV/C^{2});#it{p}_{T} (GeV/#it{c});y",
	       nbins,0.398,0.598,nbins,pMin,pMax,20,-1.0,1.0);
    fOutput->Add(fK0sMassPtRap[i]);
  } 
 
  fK0sMassPtPhi  =
    new TH3F("fK0sMassPtPhi","K^{0}_{s}: mass vs #it{p}_{T} vs #varphi;Mass (GeV/c^2);#it{p}_{T} (GeV/#it{c});#varphi (rad)",
             nbins,0.398,0.598,nbins,pMin,pMax,nbinsPhi,0.,2.*TMath::Pi());
  fOutput->Add(fK0sMassPtPhi);
  
  // ----- Fraction of shared TPC cls
  fK0sPosDaugFracShTPCcls =
    new TH3F("fK0sPosDaugFracShTPCcls","K^{0}_{s}: mass vs #it{p}_{T} vs fraction Shared TPC cls;Mass (GeV/c^2);#it{p}_{T} (GeV/#it{c});fraction Shared TPC cls",
             nbins,0.398,0.598,nbins,pMin,pMax,50,0,1.);
  fOutput->Add(fK0sPosDaugFracShTPCcls);

  fK0sNegDaugFracShTPCcls =
    new TH3F("fK0sNegDaugFracShTPCcls","K^{0}_{s}: mass vs #it{p}_{T} vs fraction Shared TPC cls;Mass (GeV/c^2);#it{p}_{T} (GeV/#it{c});fraction Shared TPC cls",
             nbins,0.398,0.598,nbins,pMin,pMax,50,0,1.);
  fOutput->Add(fK0sNegDaugFracShTPCcls);


  // ================== Correlations =================

  // ------------------------  Splitting:  
  Double_t binsDev[121];
  binsDev[0] = 0;

  for (Int_t k=-7;k<=4;k++)
    for (Int_t j=1;j<=10;j++)
      binsDev[(k+7)*10+j] = j*TMath::Power(10,k);
     
  Int_t binsSplit[9] = {100,nbins,100,2,301,101,101,120,9};   Double_t xminSplit[9] = {pMin,0.398,pMin,-0.5,-0.001,-0.005,-0.005,0,-0.5}; Double_t xmaxSplit[9] = {pMax,0.598,pMax,1.5,0.3,1.005,1.005,10e+4,8.5};    

  Int_t binsSplit2[14] = {100,nbins,100,2,10,20,101,101,100,120,99,99,99,2};   
  Double_t xminSplit2[14] = {pMin,0.398,pMin,-0.5,0.,-0.1,-0.005,-0.005,-1.,0,0.,0.,0.,-0.5}; 
  Double_t xmaxSplit2[14] = {pMax,0.598,pMax,1.5,0.1,0.1,1.005,1.005,1.,10e+4,3.3,3.3,3.3,1.5};    

  Int_t binsSplit3[6] = {100,nbins,100,46,46,2};
  Double_t xminSplit3[6] = {pMin,0.398,pMin,-0.16,-0.16,-0.5};
  Double_t xmaxSplit3[6] = {pMax,0.598,pMax,0.16,0.16,1.5};

  for(Int_t j=0; j<kNCent; j++){

    // positive daughter
    fK0sPosDaugSplCheckCovMat[j]   = new THnSparseD( Form("fK0sPosDaugSplCheckCovMat_%d",j), "K^{0}_{S} Pos. daughter; #it{p}_{T,V0} (GeV/#it{c}); Mass (GeV/c^2); #it{p}_{Daug} (GeV/#it{c}); Same Sign as Trigger Particle;  R#Delta#varphi*_{max}; Trigger: fraction of TPC shared cls; Daughter: fraction of TPC shared cls;  (X-X')^{2}/( #sigma^{2} + #sigma'^{2} ); Variables;",9,binsSplit,xminSplit,xmaxSplit);
    fK0sPosDaugSplCheckCovMat[j]->SetBinEdges(7,binsDev);
    fOutput->Add(fK0sPosDaugSplCheckCovMat[j]);  

    // negative daughter
    fK0sNegDaugSplCheckCovMat[j]   = new THnSparseD( Form("fK0sNegDaugSplCheckCovMat_%d",j), "K^{0}_{S} Neg. daughter; #it{p}_{T,V0} (GeV/#it{c}); Mass (GeV/c^2); #it{p}_{Daug} (GeV/#it{c}); Same Sign as Trigger Particle;  R#Delta#varphi*_{max}; Trigger: fraction of TPC shared cls; Daughter: fraction of TPC shared cls;  (X-X')^{2}/( #sigma^{2} + #sigma'^{2} ); Variables;",9,binsSplit,xminSplit,xmaxSplit);
    fK0sNegDaugSplCheckCovMat[j]->SetBinEdges(7,binsDev);
    fOutput->Add(fK0sNegDaugSplCheckCovMat[j]); 

    // Positive daughter:
    fK0sPosDaugdPhiSdEtaS[j]  = new THnSparseD(Form("fK0sPosDaugdPhiSdEtaS_%d",j), "K^{0}_{S} Pos. daughter; #it{p}_{T,V0} (GeV/#it{c}); Mass (GeV/c^2); #it{p}_{Daug} (GeV/#it{c}); Same Sign as Trigger Particle; #Delta#varphi*; #Delta#eta*; Trigger: fraction of TPC shared cls; Daughter: fraction of TPC shared cls; Correlation fraction of shared cls: Trigger - Daughter; #sum_{x,y,z}(#it{p}_{i,Trig}-#it{p}_{i,Daug})^{2}/( #sigma_{i,Trig}^{2} + #sigma_{i,Daug}^{2} ); DCA to prim. vtx; Trigger: DCA_{XY}; Trigger: DCA_{Z};same MC label;",14,binsSplit2,xminSplit2,xmaxSplit2);
    fK0sPosDaugdPhiSdEtaS[j]->SetBinEdges(9,binsDev);    
    fOutput->Add(fK0sPosDaugdPhiSdEtaS[j]);  
    
    // Negative daughter:
    fK0sNegDaugdPhiSdEtaS[j]  = new THnSparseD(Form("fK0sNegDaugdPhiSdEtaS_%d",j), "K^{0}_{S} Neg. daughter; #it{p}_{T,V0} (GeV/#it{c}); Mass (GeV/c^2); #it{p}_{Daug} (GeV/#it{c}); Same Sign as Trigger Particle; #Delta#varphi*; #Delta#eta*; Trigger: fraction of TPC shared cls; Daughter: fraction of TPC shared cls; Correlation fraction of shared cls: Trigger - Daughter;  #sum_{x,y,z}(#it{p}_{i,Trig}-#it{p}_{i,Daug})^{2}/( #sigma_{i,Trig}^{2} + #sigma_{i,Daug}^{2} ); DCA to prim. vtx;  Trigger: DCA_{XY}; Trigger: DCA_{Z}; same MC label;",14,binsSplit2,xminSplit2,xmaxSplit2);
    fK0sNegDaugdPhiSdEtaS[j]->SetBinEdges(9,binsDev);
    fOutput->Add(fK0sNegDaugdPhiSdEtaS[j]);  

    if(fIsMC){
      // Positive daughter:
      fK0sPosMCResdEtaSdPhiS[j]  = new THnSparseD(Form("fK0sPosMCResdEtaSdPhiS_%d",j), "K^{0}_{S} Pos. daughter; #it{p}_{T,V0} (GeV/#it{c}); Mass (GeV/c^2); #it{p}_{Daug} (GeV/#it{c}); #Delta#varphi*; #Delta#eta*; Same Sign as Trigger Particle;",6,binsSplit3,xminSplit3,xmaxSplit3);
      fOutput->Add(fK0sPosMCResdEtaSdPhiS[j]);  
    
      // Negative daughter:
      fK0sNegMCResdEtaSdPhiS[j]  = new THnSparseD(Form("fK0sNegMCResdEtaSdPhiS_%d",j), "K^{0}_{S} Neg. daughter; #it{p}_{T,V0} (GeV/#it{c}); Mass (GeV/c^2); #it{p}_{Daug} (GeV/#it{c}); #Delta#varphi*; #Delta#eta*;  Same Sign as Trigger Particle;",6,binsSplit3,xminSplit3,xmaxSplit3);
      fOutput->Add(fK0sNegMCResdEtaSdPhiS[j]);  
    }

  }


  // ----- Fraction of shared TPC cls
  fK0sPosDaugFracShTPCclsTrig =
    new TH3F("fK0sPosDaugFracShTPCclsTrig","K^{0}_{s}: mass vs #it{p}_{T} vs fraction Shared TPC cls;Mass (GeV/c^2);#it{p}_{T} (GeV/#it{c});fraction Shared TPC cls",
             nbins,0.398,0.598,nbins,pMin,pMax,50,0,1.);
  fOutput->Add(fK0sPosDaugFracShTPCclsTrig);

  fK0sNegDaugFracShTPCclsTrig =
    new TH3F("fK0sNegDaugFracShTPCclsTrig","K^{0}_{s}: mass vs #it{p}_{T} vs fraction Shared TPC cls;Mass (GeV/c^2);#it{p}_{T} (GeV/#it{c});fraction Shared TPC cls",
             nbins,0.398,0.598,nbins,pMin,pMax,50,0,1.);
  fOutput->Add(fK0sNegDaugFracShTPCclsTrig);


  //    DCA to prim vertex
  fK0sDCADaugToPrimVtx  
    = new TH3F("fK0sDCADaugToPrimVtx","K^{0}_{S} Bckg: dca daughter vs. #it{p}_{T,l};DCA Pos daug (cm);DCA Neg daug (cm);#it{p}_{T,l} (GeV/#it{c})",
	       90,0.,3.3,90,0.,3.3,nbinPtLP,pMin,ptMaxLP);
  fOutput->Add(fK0sDCADaugToPrimVtx);

  //    Spatial Resoltuion between trigger- and asosciated- particles
  fK0sSpatialRes = new TH3F("fK0sSpatialRes","K^{0}_{S}: Spatial resolution;#Delta#varphi (rad);trig-assoc. resolution (cm);dec. length (cm)",
			    20,-0.1,0.1,100,0.,10,2*nbins,lMin,lMax);
  fOutput->Add(fK0sSpatialRes);

  for(Int_t jj=0;jj<kNCent;jj++){
    for(Int_t k=0;k<kN1;k++){

      // Monte-Carlo level:
      if(fIsMC){
	snprintf(hNameHist,100, "fK0sdPhidEtaMC_%.2f_%.2f_Cent_%.0f_%.0f",kPtBinV0[k],kPtBinV0[k+1],kBinCent[jj],kBinCent[jj+1]); 
	fK0sdPhidEtaMC[jj*kN1+k] = new TH3F(hNameHist,"K^{0}_{S} MC: #Delta#varphi vs #Delta#eta vs #it{p}_{T,l}",
					    nbinsdPhi,-TMath::PiOver2(),3*TMath::PiOver2(),
					    nbinsdEta,-1.5,1.5,
					    nbinsVtx,-10.,10.);
	fK0sdPhidEtaMC[jj*kN1+k]->GetXaxis()->SetTitle("#Delta#varphi (rad)"); 
	fK0sdPhidEtaMC[jj*kN1+k]->GetYaxis()->SetTitle("#Delta#eta"); 
	fK0sdPhidEtaMC[jj*kN1+k]->GetZaxis()->SetTitle("Vertex Z (cm)"); 
	fOutput->Add(fK0sdPhidEtaMC[jj*kN1+k]);
      }
  
      // Reconstruction level:
      for(Int_t ll=0;ll<kNVtxZ;ll++){
	snprintf(hNameHist,100, "fK0sdPhidEtaPtL_%.2f_%.2f_Cent_%.0f_%.0f_%d",kPtBinV0[k],kPtBinV0[k+1],kBinCent[jj],kBinCent[jj+1],ll); 
	fK0sdPhidEtaPtL[jj*kN1*kNVtxZ  + k*kNVtxZ + ll] = new TH3F(hNameHist,"K^{0}_{S}: #Delta#varphi vs #Delta#eta vs Inv. Mass",
								   nbinsdPhi,-TMath::PiOver2(),3*TMath::PiOver2(),
								   nbinsdEta,-1.5,1.5,
								   nbins,0.398,0.598);
	fK0sdPhidEtaPtL[jj*kN1*kNVtxZ  + k*kNVtxZ + ll]->GetXaxis()->SetTitle("#Delta#varphi (rad)"); 
	fK0sdPhidEtaPtL[jj*kN1*kNVtxZ  + k*kNVtxZ + ll]->GetYaxis()->SetTitle("#Delta#eta"); 
	fK0sdPhidEtaPtL[jj*kN1*kNVtxZ  + k*kNVtxZ + ll]->GetZaxis()->SetTitle("Inv. Mass"); 
	fOutput->Add(fK0sdPhidEtaPtL[jj*kN1*kNVtxZ  + k*kNVtxZ + ll]);
      }
    }
  }

  // Correlations (side-band):
  fK0sBckgDecLength
    = new TH2F("fK0sBckgDecLength","K^{0}_{S} Bckg: c#tau vs. #it{p}_{T,l}",
	       100,0.,15.,nbinPtLP,pMin,ptMaxLP);
  fK0sBckgDecLength->GetXaxis()->SetTitle("c#tau (cm)"); 
  fK0sBckgDecLength->GetYaxis()->SetTitle("#it{p}_{T,l} (GeV/#it{c})"); 
  fOutput->Add(fK0sBckgDecLength);

  fK0sBckgDCADaugToPrimVtx  
    = new TH3F("fK0sBckgDCADaugToPrimVtx","K^{0}_{S} Bckg: dca daughter vs. #it{p}_{T,l}",
	       90,0.,3.3,90,0.,3.3,nbinPtLP,pMin,ptMaxLP);
  fK0sBckgDCADaugToPrimVtx->GetXaxis()->SetTitle("DCA Pos daug (cm)"); 
  fK0sBckgDCADaugToPrimVtx->GetYaxis()->SetTitle("DCA Neg daug (cm)"); 
  fK0sBckgDCADaugToPrimVtx->GetZaxis()->SetTitle("#it{p}_{T,l} (GeV/#it{c})"); 
  fOutput->Add(fK0sBckgDCADaugToPrimVtx);
  
  fK0sBckgEtaPhi = 
    new TH2F("fK0sBckgEtaPhi","K^{0}_{s} Bckg: #varphi vs #eta",
	     nbinsPhi,0.,2.*TMath::Pi(),100,-1.,1.);
  fK0sBckgEtaPhi->GetXaxis()->SetTitle("#varphi (rad)"); 
  fK0sBckgEtaPhi->GetYaxis()->SetTitle("#eta"); 
  fOutput->Add(fK0sBckgEtaPhi);

  fK0sBckgPhiRadio
    = new TH2F("fK0sBckgPhiRadio","K^{0}_{S} Bckg: #varphi vs l_{T}",
	       nbinsPhi,0.,2.*TMath::Pi(),2*nbins,lMin,lMax);
  fK0sBckgPhiRadio->GetXaxis()->SetTitle("#varphi (rad)"); 
  fK0sBckgPhiRadio->GetYaxis()->SetTitle("l_{T} (cm)"); 
  fOutput->Add(fK0sBckgPhiRadio);
 
  fK0sBckgDCANegDaugToPrimVtx  
    = new TH2F("fK0sBckgDCANegDaugToPrimVtx","K^{0}_{S} Bckg: dca NegDaughter",
	       7,-0.5,6.5,90,0.,3.3);
  fK0sBckgDCANegDaugToPrimVtx->GetXaxis()->SetTitle("MC Production"); 
  fK0sBckgDCANegDaugToPrimVtx->GetXaxis()->SetBinLabel(1,"Rec");
  fK0sBckgDCANegDaugToPrimVtx->GetXaxis()->SetBinLabel(2,"Primary");
  fK0sBckgDCANegDaugToPrimVtx->GetXaxis()->SetBinLabel(3,"V0's");
  fK0sBckgDCANegDaugToPrimVtx->GetXaxis()->SetBinLabel(4,"Cascades");
  fK0sBckgDCANegDaugToPrimVtx->GetXaxis()->SetBinLabel(5,"Gamma conv.");
  fK0sBckgDCANegDaugToPrimVtx->GetXaxis()->SetBinLabel(6,"Unidentified mother");
  fK0sBckgDCANegDaugToPrimVtx->GetXaxis()->SetBinLabel(7,"Other");
  fK0sBckgDCANegDaugToPrimVtx->GetYaxis()->SetTitle("DCA Neg Daug (cm)"); 
  fOutput->Add(fK0sBckgDCANegDaugToPrimVtx);

  fK0sBckgDCAPosDaugToPrimVtx  
    = new TH2F("fK0sBckgDCAPosDaugToPrimVtx","K^{0}_{S} Bckg: dca PosDaughter",
	       7,-0.5,6.5,90,0.,3.3);
  fK0sBckgDCAPosDaugToPrimVtx->GetXaxis()->SetTitle("MC Production"); 
  fK0sBckgDCAPosDaugToPrimVtx->GetXaxis()->SetBinLabel(1,"Rec");
  fK0sBckgDCAPosDaugToPrimVtx->GetXaxis()->SetBinLabel(2,"Primary");
  fK0sBckgDCAPosDaugToPrimVtx->GetXaxis()->SetBinLabel(3,"V0's");
  fK0sBckgDCAPosDaugToPrimVtx->GetXaxis()->SetBinLabel(4,"Cascades");
  fK0sBckgDCAPosDaugToPrimVtx->GetXaxis()->SetBinLabel(5,"Gamma conv.");
  fK0sBckgDCAPosDaugToPrimVtx->GetXaxis()->SetBinLabel(6,"Unidentified mother");
  fK0sBckgDCAPosDaugToPrimVtx->GetXaxis()->SetBinLabel(7,"Other");
  fK0sBckgDCAPosDaugToPrimVtx->GetYaxis()->SetTitle("DCA Pos Daug (cm)"); 
  fOutput->Add(fK0sBckgDCAPosDaugToPrimVtx);
        
  fV0MassCascade
    = new TH2F("fV0MassCascade","Cascade Reconstruction wi V0's candiates;Invariant Mass (GeV/c^{2});Cascade type",650, 1.2, 2.5,12,0.5,12.5);
  fOutput->Add(fV0MassCascade);


  // ****** Lambda ******
  fLambdaMass = 
    new TH3F("fLambdaMass","Mass vs #it{p}_{T} for \\Lambda",nbins,1.065,1.165,nbins,pMin,pMax,100,0.,100.);
  fLambdaMass->GetXaxis()->SetTitle("Mass (GeV/c^2)");
  fLambdaMass->GetYaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})"); 
  fLambdaMass->GetZaxis()->SetTitle("centrality"); 
  fOutput->Add(fLambdaMass);
  
  fLambdaMassEmbeded =
    new TH3F("fLambdaMassEmbeded","Mass vs #it{p}_{T} for \\Lambda Embeded",nbins,1.065,1.165,nbins,pMin,pMax,100,0.,100.);
  fLambdaMassEmbeded->GetXaxis()->SetTitle("Mass (GeV/c^2)");
  fLambdaMassEmbeded->GetYaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fLambdaMassEmbeded->GetZaxis()->SetTitle("centrality");
  fOutput->Add(fLambdaMassEmbeded);

  fLambdaMass2 =
    new TH3F("fLambdaMass2","Mass vs #it{p}_{T} for \\Lambda",nbins,1.065,1.165,nbins,pMin,pMax,100,0.,100.);
  fLambdaMass2->GetXaxis()->SetTitle("Mass (GeV/c^2)");
  fLambdaMass2->GetYaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fLambdaMass2->GetZaxis()->SetTitle("centrality");
  fOutput->Add(fLambdaMass2);

  fLambdaMass2Embeded =
    new TH3F("fLambdaMass2Embeded","Mass vs #it{p}_{T} for \\Lambda Embeded",nbins,1.065,1.165,nbins,pMin,pMax,100,0.,100.);
  fLambdaMass2Embeded->GetXaxis()->SetTitle("Mass (GeV/c^2)");
  fLambdaMass2Embeded->GetYaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fLambdaMass2Embeded->GetZaxis()->SetTitle("centrality");
  fOutput->Add(fLambdaMass2Embeded);

  fLambdaMassPtEta =
    new TH3F("fLambdaMassPtEta","\\Lambda: mass vs #it{p}_{T} vs #eta;Mass (GeV/#it{c}^2);#it{p}_{T} (GeV/#it{c});#eta",
	     nbins,1.065,1.165,nbins,pMin,pMax,20,-1.0,1.0);
  fOutput->Add(fLambdaMassPtEta);

  for(Int_t i=0; i<kNCent; i++){
    fLambdaMassPtRap[i] =
      new TH3F(Form("fLambdaMassPtRap_cent_%.0lf_%.0lf",kBinCent[i],kBinCent[i+1]),"\\Lambda: mass vs #it{p}_{T} vs y;Mass (GeV/#it{c}^2);#it{p}_{T} (GeV/#it{c});y",
	       nbins,1.065,1.165,nbins,pMin,pMax,20,-1.0,1.0);
    fOutput->Add(fLambdaMassPtRap[i]);
  }

  fLambdaMassPtPhi  = 
    new TH3F("fLambdaMassPtPhi","#Lambda: mass vs #it{p}_{T} vs #varphi;Mass (GeV/#it{c}^2);#it{p}_{T} (GeV/#it{c});#varphi (rad)",
	     nbins,1.065,1.165,nbins,pMin,pMax,nbinsPhi,0.,2.*TMath::Pi());
  fOutput->Add(fLambdaMassPtPhi);

  
  // ----- Fraction of shared TPC cls
  fLambdaPosDaugFracShTPCcls =
    new TH3F("fLambdaPosDaugFracShTPCcls","#Lambda: mass vs #it{p}_{T} vs fraction Shared TPC cls;Mass (GeV/c^2);#it{p}_{T} (GeV/#it{c});fraction Shared TPC cls",
         	     nbins,1.065,1.165,nbins,pMin,pMax,50,0,1.);
  fOutput->Add(fLambdaPosDaugFracShTPCcls);

  fLambdaNegDaugFracShTPCcls =
    new TH3F("fLambdaNegDaugFracShTPCcls","#Lambda: mass vs #it{p}_{T} vs fraction Shared TPC cls;Mass (GeV/c^2);#it{p}_{T} (GeV/#it{c});fraction Shared TPC cls",
          	     nbins,1.065,1.165,nbins,pMin,pMax,50,0,1.);
  fOutput->Add(fLambdaNegDaugFracShTPCcls);



  // ================== Correlations =================

  // ----------------Splitting:
  xminSplit[1] = 1.065;   xmaxSplit[1] = 1.165;    
  xminSplit2[1] = 1.065;  xmaxSplit2[1] = 1.165;
  xminSplit3[1] = 1.065;  xmaxSplit3[1] = 1.165;

  for(Int_t j=0; j<kNCent; j++){

    // positive daughter:
    fLambdaPosDaugSplCheckCovMat[j]   = new THnSparseD( Form("fLambdaPosDaugSplCheckCovMat_%d",j), "#Lambda Pos. daughter;   #it{p}_{T,V0} (GeV/#it{c}); Mass (GeV/c^2); #it{p}_{Daug} (GeV/#it{c}); Same Sign as Trigger Particle;  R#Delta#varphi*_{max}; Trigger: fraction of TPC shared cls; Daughter: fraction of TPC shared cls;  (X-X')^{2}/( #sigma^{2} + #sigma'^{2} ); Variables;",9,binsSplit,xminSplit,xmaxSplit);   
    fLambdaPosDaugSplCheckCovMat[j]->SetBinEdges(7,binsDev);
    fOutput->Add(fLambdaPosDaugSplCheckCovMat[j]);  

    // negative daughter:
    fLambdaNegDaugSplCheckCovMat[j]   = new THnSparseD( Form("fLambdaNegDaugSplCheckCovMat_%d",j), "#Lambda Neg. daughter;   #it{p}_{T,V0} (GeV/#it{c}); Mass (GeV/c^2); #it{p}_{Daug} (GeV/#it{c}); Same Sign as Trigger Particle;  R#Delta#varphi*_{max}; Trigger: fraction of TPC shared cls; Daughter: fraction of TPC shared cls;  (X-X')^{2}/( #sigma^{2} + #sigma'^{2} ); Variables;",9,binsSplit,xminSplit,xmaxSplit);   
    fLambdaNegDaugSplCheckCovMat[j]->SetBinEdges(7,binsDev);
    fOutput->Add(fLambdaNegDaugSplCheckCovMat[j]); 

    // Positive daughter:
    fLambdaPosDaugdPhiSdEtaS[j]  = new THnSparseD(Form("fLambdaPosDaugdPhiSdEtaS_%d",j), "#Lambda Pos. daughter; #it{p}_{T,V0} (GeV/#it{c}); Mass (GeV/c^2); #it{p}_{Daug} (GeV/#it{c}); Same Sign as Trigger Particle; #Delta#varphi*; #Delta#eta*; Trigger: fraction of TPC shared cls; Daughter: fraction of TPC shared cls;  Correlation fraction of shared cls: Trigger - Daughter; #sum_{x,y,z}(#it{p}_{i,Trig}-#it{p}_{i,Daug})^{2}/( #sigma_{i,Trig}^{2} + #sigma_{i,Daug}^{2} ); DCA to prim. vtx;  Trigger: DCA_{XY}; Trigger: DCA_{Z}; same MC label;",14,binsSplit2,xminSplit2,xmaxSplit2);
    fLambdaPosDaugdPhiSdEtaS[j]->SetBinEdges(9,binsDev);
    fOutput->Add(fLambdaPosDaugdPhiSdEtaS[j]);  
    
    // Negative daughter:
    fLambdaNegDaugdPhiSdEtaS[j]  = new THnSparseD(Form("fLambdaNegDaugdPhiSdEtaS_%d",j), "#Lambda Neg. daughter; #it{p}_{T,V0} (GeV/#it{c}); Mass (GeV/c^2); #it{p}_{Daug} (GeV/#it{c}); Same Sign as Trigger Particle; #Delta#varphi*; #Delta#eta*; Trigger: fraction of TPC shared cls; Daughter: fraction of TPC shared cls;  Correlation fraction of shared cls: Trigger - Daughter; #sum_{x,y,z}(#it{p}_{i,Trig}-#it{p}_{i,Daug})^{2}/( #sigma_{i,Trig}^{2} + #sigma_{i,Daug}^{2} ); DCA to prim. vtx;  Trigger: DCA_{XY}; Trigger: DCA_{Z}; same MC label;",14,binsSplit2,xminSplit2,xmaxSplit2);
    fLambdaNegDaugdPhiSdEtaS[j]->SetBinEdges(9,binsDev);
    fOutput->Add(fLambdaNegDaugdPhiSdEtaS[j]);  

    if(fIsMC){
      // Positive daughter:
      fLambdaPosMCResdEtaSdPhiS[j]  = new THnSparseD(Form("fLambdaPosMCResdEtaSdPhiS_%d",j), "#Lambda Pos. daughter; #it{p}_{T,V0} (GeV/#it{c}); Mass (GeV/c^2); #it{p}_{Daug} (GeV/#it{c}); #Delta#varphi*; #Delta#eta*; Same Sign as Trigger Particle;",6,binsSplit3,xminSplit3,xmaxSplit3);
      fOutput->Add(fLambdaPosMCResdEtaSdPhiS[j]);  
    
      // Negative daughter:
      fLambdaNegMCResdEtaSdPhiS[j]  = new THnSparseD(Form("fLambdaNegMCResdEtaSdPhiS_%d",j), "#Lambda Neg. daughter; #it{p}_{T,V0} (GeV/#it{c}); Mass (GeV/c^2); #it{p}_{Daug} (GeV/#it{c}); #Delta#varphi*; #Delta#eta*;  Same Sign as Trigger Particle;",6,binsSplit3,xminSplit3,xmaxSplit3);
      fOutput->Add(fLambdaNegMCResdEtaSdPhiS[j]);  
    }

  }


  // ----- Fraction of shared TPC cls
  fLambdaPosDaugFracShTPCclsTrig =
    new TH3F("fLambdaPosDaugFracShTPCclsTrig","#Lambda: mass vs #it{p}_{T} vs fraction Shared TPC cls;Mass (GeV/c^2);#it{p}_{T} (GeV/#it{c});fraction Shared TPC cls",
         	     nbins,1.065,1.165,nbins,pMin,pMax,50,0,1.);
  fOutput->Add(fLambdaPosDaugFracShTPCclsTrig);

  fLambdaNegDaugFracShTPCclsTrig =
    new TH3F("fLambdaNegDaugFracShTPCclsTrig","#Lambda: mass vs #it{p}_{T} vs fraction Shared TPC cls;Mass (GeV/c^2);#it{p}_{T} (GeV/#it{c});fraction Shared TPC cls",
          	     nbins,1.065,1.165,nbins,pMin,pMax,50,0,1.);
 fOutput->Add(fLambdaNegDaugFracShTPCclsTrig);


  //    DCA to prim vertex
  fLambdaDCADaugToPrimVtx  
    = new TH3F("fLambdaDCADaugToPrimVtx","#Lambda Bckg: dca daughter vs. #it{p}_{T,l}",
	       90,0.,3.3,90,0.,3.3,nbinPtLP,pMin,ptMaxLP);
  fLambdaDCADaugToPrimVtx->GetXaxis()->SetTitle("DCA Pos daug (cm)"); 
  fLambdaDCADaugToPrimVtx->GetYaxis()->SetTitle("DCA Neg daug (cm)"); 
  fLambdaDCADaugToPrimVtx->GetZaxis()->SetTitle("#it{p}_{T,l} (GeV/#it{c})"); 
  fOutput->Add(fLambdaDCADaugToPrimVtx);

  //    Spatial Resoltuion between trigger- and asosciated- particles
  fLambdaSpatialRes = new TH3F("fLambdaSpatialRes","#Lambda: Spatial resolution;#Delta#varphi (rad);trig-assoc. resolution (cm);dec. length (cm)",
			       20,-0.1,0.1,100,0.,10,2*nbins,lMin,lMax);
  fOutput->Add(fLambdaSpatialRes);


  for(Int_t jj=0;jj<kNCent;jj++){
    for(Int_t k=0;k<kN1;k++){

      // Monte-Carlo level:
      if(fIsMC){
	snprintf(hNameHist,100, "fLambdadPhidEtaMC_%.2f_%.2f_Cent_%.0f_%.0f",kPtBinV0[k],kPtBinV0[k+1],kBinCent[jj],kBinCent[jj+1]); 
	fLambdadPhidEtaMC[jj*kN1+k] = new TH3F(hNameHist,"#Lambda MC: #Delta#varphi vs #Delta#eta vs #it{p}_{T,l}",
					       nbinsdPhi,-TMath::PiOver2(),3*TMath::PiOver2(),
					       nbinsdEta,-1.5,1.5,
					       nbinsVtx,-10.,10.);
	fLambdadPhidEtaMC[jj*kN1+k]->GetXaxis()->SetTitle("#Delta#varphi (rad)"); 
	fLambdadPhidEtaMC[jj*kN1+k]->GetYaxis()->SetTitle("#Delta#eta"); 
	fLambdadPhidEtaMC[jj*kN1+k]->GetZaxis()->SetTitle("Vertex Z (cm)"); 
	fOutput->Add(fLambdadPhidEtaMC[jj*kN1+k]);
      }

      // Reconstruction level:
      for(Int_t ll=0;ll<kNVtxZ;ll++){
	snprintf(hNameHist,100, "fLambdadPhidEtaPtL_%.2f_%.2f_Cent_%.0f_%.0f_%d",kPtBinV0[k],kPtBinV0[k+1],kBinCent[jj],kBinCent[jj+1],ll); 
	fLambdadPhidEtaPtL[jj*kN1*kNVtxZ  + k*kNVtxZ + ll] = new TH3F(hNameHist,"#Lambda: #Delta#varphi vs #Delta#eta vs #it{p}_{T,l}",
								      nbinsdPhi,-TMath::PiOver2(),3*TMath::PiOver2(),
								      nbinsdEta,-1.5,1.5,
								      nbins,1.065,1.165);
	fLambdadPhidEtaPtL[jj*kN1*kNVtxZ  + k*kNVtxZ + ll]->GetXaxis()->SetTitle("#Delta#varphi (rad)"); 
	fLambdadPhidEtaPtL[jj*kN1*kNVtxZ  + k*kNVtxZ + ll]->GetYaxis()->SetTitle("#Delta#eta"); 
	fLambdadPhidEtaPtL[jj*kN1*kNVtxZ  + k*kNVtxZ + ll]->GetZaxis()->SetTitle("Inv. Mass");
	fOutput->Add(fLambdadPhidEtaPtL[jj*kN1*kNVtxZ  + k*kNVtxZ + ll]);
      }
    }
  }

  // Correlations (side-band):
  fLambdaBckgDecLength
    = new TH2F("fLambdaBckgDecLength","#Lambda Bckg: c#tau vs. #it{p}_{T,l}",
	       100,0.,25.,nbinPtLP,pMin,ptMaxLP);
  fLambdaBckgDecLength->GetXaxis()->SetTitle("c#tau (cm)"); 
  fLambdaBckgDecLength->GetYaxis()->SetTitle("#it{p}_{T,l} (GeV/#it{c})"); 
  fOutput->Add(fLambdaBckgDecLength);
  
  fLambdaBckgDCADaugToPrimVtx  
    = new TH3F("fLambdaBckgDCADaugToPrimVtx","#Lambda Bckg: dca daughter vs. #it{p}_{T,l}",
	       90,0.,3.3,90,0.,3.3,nbinPtLP,pMin,ptMaxLP);
  fLambdaBckgDCADaugToPrimVtx->GetXaxis()->SetTitle("DCA Pos daug (cm)"); 
  fLambdaBckgDCADaugToPrimVtx->GetYaxis()->SetTitle("DCA Neg daug (cm)"); 
  fLambdaBckgDCADaugToPrimVtx->GetZaxis()->SetTitle("#it{p}_{T,l} (GeV/#it{c})"); 
  fOutput->Add(fLambdaBckgDCADaugToPrimVtx);
  
  fLambdaBckgEtaPhi = 
    new TH2F("fLambdaBckgEtaPhi","#Lambda Bckg: #varphi vs #eta",
	     nbinsPhi,0.,2.*TMath::Pi(),100,-1.,1.);
  fLambdaBckgEtaPhi->GetXaxis()->SetTitle("#varphi (rad)"); 
  fLambdaBckgEtaPhi->GetYaxis()->SetTitle("#eta"); 
  fOutput->Add(fLambdaBckgEtaPhi);
    
  fLambdaBckgPhiRadio
    = new TH2F("fLambdaBckgPhiRadio","#Lambda Bckg: #varphi vs l_{T}",
	       nbinsPhi,0.,2.*TMath::Pi(),2*nbins,lMin,lMax);
  fLambdaBckgPhiRadio->GetXaxis()->SetTitle("#varphi (rad)"); 
  fLambdaBckgPhiRadio->GetYaxis()->SetTitle("l_{T} (cm)"); 
  fOutput->Add(fLambdaBckgPhiRadio);


  fLambdaBckgDCANegDaugToPrimVtx  
    = new TH2F("fLambdaBckgDCANegDaugToPrimVtx","#Lambda Bckg: dca NegDaughter",
	       7,-0.5,6.5,90,0.,3.3);
  fLambdaBckgDCANegDaugToPrimVtx->GetXaxis()->SetTitle("MC Production"); 
  fLambdaBckgDCANegDaugToPrimVtx->GetXaxis()->SetBinLabel(1,"Rec");
  fLambdaBckgDCANegDaugToPrimVtx->GetXaxis()->SetBinLabel(2,"Primary");
  fLambdaBckgDCANegDaugToPrimVtx->GetXaxis()->SetBinLabel(3,"V0's");
  fLambdaBckgDCANegDaugToPrimVtx->GetXaxis()->SetBinLabel(4,"Cascades");
  fLambdaBckgDCANegDaugToPrimVtx->GetXaxis()->SetBinLabel(5,"Gamma conv.");
  fLambdaBckgDCANegDaugToPrimVtx->GetXaxis()->SetBinLabel(6,"Unidentified mother");
  fLambdaBckgDCANegDaugToPrimVtx->GetXaxis()->SetBinLabel(7,"Other");
  fLambdaBckgDCANegDaugToPrimVtx->GetYaxis()->SetTitle("DCA Neg Daug (cm)"); 
  fOutput->Add(fLambdaBckgDCANegDaugToPrimVtx);


  fLambdaBckgDCAPosDaugToPrimVtx  
    = new TH2F("fLambdaBckgDCAPosDaugToPrimVtx","#Lambda Bckg: dca PosDaughter",
	       7,-0.5,6.5,90,0.,3.3);
  fLambdaBckgDCAPosDaugToPrimVtx->GetXaxis()->SetTitle("MC Production"); 
  fLambdaBckgDCAPosDaugToPrimVtx->GetXaxis()->SetBinLabel(1,"Rec");
  fLambdaBckgDCAPosDaugToPrimVtx->GetXaxis()->SetBinLabel(2,"Primary");
  fLambdaBckgDCAPosDaugToPrimVtx->GetXaxis()->SetBinLabel(3,"V0's");
  fLambdaBckgDCAPosDaugToPrimVtx->GetXaxis()->SetBinLabel(4,"Cascades");
  fLambdaBckgDCAPosDaugToPrimVtx->GetXaxis()->SetBinLabel(5,"Gamma conv.");
  fLambdaBckgDCAPosDaugToPrimVtx->GetXaxis()->SetBinLabel(6,"Unidentified mother");
  fLambdaBckgDCAPosDaugToPrimVtx->GetXaxis()->SetBinLabel(7,"Other");
  fLambdaBckgDCAPosDaugToPrimVtx->GetYaxis()->SetTitle("DCA Pos Daug (cm)"); 
  fOutput->Add(fLambdaBckgDCAPosDaugToPrimVtx);


  // ****** AntiLambda ******
  fAntiLambdaMass = 
    new TH3F("fAntiLambdaMass","Mass vs #it{p}_{T} for #bar{#Lambda}",nbins,1.065,1.165,nbins,pMin,pMax,100,0.,100.);
  fAntiLambdaMass->GetXaxis()->SetTitle("Mass (GeV/c^2)");
  fAntiLambdaMass->GetYaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})"); 
  fAntiLambdaMass->GetZaxis()->SetTitle("centrality"); 
  fOutput->Add(fAntiLambdaMass);
  
  fAntiLambdaMassEmbeded =
    new TH3F("fAntiLambdaMassEmbeded","Mass vs #it{p}_{T} for #bar{#Lambda} Embeded",nbins,1.065,1.165,nbins,pMin,pMax,100,0.,100.);
  fAntiLambdaMassEmbeded->GetXaxis()->SetTitle("Mass (GeV/c^2)");
  fAntiLambdaMassEmbeded->GetYaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fAntiLambdaMassEmbeded->GetZaxis()->SetTitle("centrality");
  fOutput->Add(fAntiLambdaMassEmbeded);

  fAntiLambdaMass2 =
    new TH3F("fAntiLambdaMass2","Mass vs #it{p}_{T} for #bar{#Lambda}",nbins,1.065,1.165,nbins,pMin,pMax,100,0.,100.);
  fAntiLambdaMass2->GetXaxis()->SetTitle("Mass (GeV/c^2)");
  fAntiLambdaMass2->GetYaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fAntiLambdaMass2->GetZaxis()->SetTitle("centrality");
  fOutput->Add(fAntiLambdaMass2);  

  fAntiLambdaMass2Embeded =
    new TH3F("fAntiLambdaMass2Embeded","Mass vs #it{p}_{T} for #bar{#Lambda} Embeded",nbins,1.065,1.165,nbins,pMin,pMax,100,0.,100.);
  fAntiLambdaMass2Embeded->GetXaxis()->SetTitle("Mass (GeV/c^2)");
  fAntiLambdaMass2Embeded->GetYaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fAntiLambdaMass2Embeded->GetZaxis()->SetTitle("centrality");
  fOutput->Add(fAntiLambdaMass2Embeded);  

  fAntiLambdaMassPtEta =
    new TH3F("fAntiLambdaMassPtEta","#bar{#Lambda}: mass vs #it{p}_{T} vs #eta;Mass (GeV/#it{c}^2);#it{p}_{T} (GeV/#it{c});#eta",nbins,1.065,1.165,nbins,pMin,pMax,20,-1.0,1.0);
  fOutput->Add(fAntiLambdaMassPtEta);

  for(Int_t i=0; i<kNCent; i++){
    fAntiLambdaMassPtRap[i] =
      new TH3F(Form("fAntiLambdaMassPtRap_cent_%.0lf_%.0lf",kBinCent[i],kBinCent[i+1]),"#bar{#Lambda}: mass vs #it{p}_{T} vs y;Mass (GeV/#it{c}^2);#it{p}_{T} (GeV/#it{c});y",nbins,1.065,1.165,nbins,pMin,pMax,20,-1.0,1.0);
    fOutput->Add(fAntiLambdaMassPtRap[i]);
  }

  fAntiLambdaMassPtPhi  = 
    new TH3F("fAntiLambdaMassPtPhi","#bar{#Lambda}: mass vs #it{p}_{T} vs #varphi;Mass (GeV/#it{c}^2);#it{p}_{T} (GeV/#it{c});#varphi (rad)",
	     nbins,1.065,1.165,nbins,pMin,pMax,nbinsPhi,0.,2.*TMath::Pi());
  fOutput->Add(fAntiLambdaMassPtPhi);

  fAntiLambdaPosDaugFracShTPCcls =
    new TH3F("fAntiLambdaPosDaugFracShTPCcls","#bar{#Lambda}: mass vs #it{p}_{T} vs fraction Shared TPC cls;Mass (GeV/c^2);#it{p}_{T} (GeV/#it{c});fraction Shared TPC cls",
         	     nbins,1.065,1.165,nbins,pMin,pMax,50,0,1.);
  fOutput->Add(fAntiLambdaPosDaugFracShTPCcls);

  fAntiLambdaNegDaugFracShTPCcls =
    new TH3F("fAntiLambdaNegDaugFracShTPCcls","#bar{#Lambda}: mass vs #it{p}_{T} vs fraction Shared TPC cls;Mass (GeV/c^2);#it{p}_{T} (GeV/#it{c});fraction Shared TPC cls",
          	     nbins,1.065,1.165,nbins,pMin,pMax,50,0,1.);
  fOutput->Add(fAntiLambdaNegDaugFracShTPCcls);


  // ================== Correlations =================

  // ----------------Splitting:
  for(Int_t j=0; j<kNCent; j++){

    // positive daughter:
    fAntiLambdaPosDaugSplCheckCovMat[j]  = new THnSparseD(Form("fAntiLambdaPosDaugSplCheckCovMat_%d",j), "#bar{#Lambda} Pos. daughter;  #it{p}_{T,V0} (GeV/#it{c}); Mass (GeV/c^2); #it{p}_{Daug} (GeV/#it{c}); Same Sign as Trigger Particle;  R#Delta#varphi*_{max}; Trigger: fraction of TPC shared cls; Daughter: fraction of TPC shared cls;  (X-X')^{2}/( #sigma^{2} + #sigma'^{2} ); Variables;",9,binsSplit,xminSplit,xmaxSplit); 
    fAntiLambdaPosDaugSplCheckCovMat[j]->SetBinEdges(7,binsDev);
    fOutput->Add(fAntiLambdaPosDaugSplCheckCovMat[j]);  

    // negative daughter:
    fAntiLambdaNegDaugSplCheckCovMat[j]  = new THnSparseD(Form("fAntiLambdaNegDaugSplCheckCovMat_%d",j), "#bar{#Lambda} Neg. daughter;  #it{p}_{T,V0} (GeV/#it{c}); Mass (GeV/c^2); #it{p}_{Daug} (GeV/#it{c}); Same Sign as Trigger Particle;  R#Delta#varphi*_{max}; Trigger: fraction of TPC shared cls; Daughter: fraction of TPC shared cls;  (X-X')^{2}/( #sigma^{2} + #sigma'^{2} ); Variables;",9,binsSplit,xminSplit,xmaxSplit);       
    fAntiLambdaNegDaugSplCheckCovMat[j]->SetBinEdges(7,binsDev);
    fOutput->Add(fAntiLambdaNegDaugSplCheckCovMat[j]); 

    // Positive daughter:
    fAntiLambdaPosDaugdPhiSdEtaS[j]  = new THnSparseD(Form("fAntiLambdaPosDaugdPhiSdEtaS_%d",j), "#bar{#Lambda} Pos. daughter; #it{p}_{T,V0} (GeV/#it{c}); Mass (GeV/c^2); #it{p}_{Daug} (GeV/#it{c}); Same Sign as Trigger Particle; #Delta#varphi*; #Delta#eta*; Trigger: fraction of TPC shared cls; Daughter: fraction of TPC shared cls; Correlation fraction of shared cls: Trigger - Daughter;  Trigger: DCA_{XY}; Trigger: DCA_{Z}; same MC label;",14,binsSplit2,xminSplit2,xmaxSplit2);
    fAntiLambdaPosDaugdPhiSdEtaS[j]->SetBinEdges(9,binsDev);  
      fOutput->Add(fAntiLambdaPosDaugdPhiSdEtaS[j]);  
    
    // Negative daughter:
    fAntiLambdaNegDaugdPhiSdEtaS[j]  = new THnSparseD(Form("fAntiLambdaNegDaugdPhiSdEtaS_%d",j), "#bar{#Lambda} Neg. daughter; #it{p}_{T,V0} (GeV/#it{c}); Mass (GeV/c^2); #it{p}_{Daug} (GeV/#it{c}); Same Sign as Trigger Particle; #Delta#varphi*; #Delta#eta*; Trigger: fraction of TPC shared cls; Daughter: fraction of TPC shared cls; Correlation fraction of shared cls: Trigger - Daughter;  Trigger: DCA_{XY}; Trigger: DCA_{Z}; same MC label;",14,binsSplit2,xminSplit2,xmaxSplit2);
    fAntiLambdaNegDaugdPhiSdEtaS[j]->SetBinEdges(9,binsDev);    
      fOutput->Add(fAntiLambdaNegDaugdPhiSdEtaS[j]);  

   if(fIsMC){
      // Positive daughter:
      fAntiLambdaPosMCResdEtaSdPhiS[j]  = new THnSparseD(Form("fAntiLambdaPosMCResdEtaSdPhiS_%d",j), "#bar{#Lambda} Pos. daughter; #it{p}_{T,V0} (GeV/#it{c}); Mass (GeV/c^2); #it{p}_{Daug} (GeV/#it{c}); #Delta#varphi*; #Delta#eta*; Same Sign as Trigger Particle;",6,binsSplit3,xminSplit3,xmaxSplit3);
      fOutput->Add(fAntiLambdaPosMCResdEtaSdPhiS[j]);  
    
      // Negative daughter:
      fAntiLambdaNegMCResdEtaSdPhiS[j]  = new THnSparseD(Form("fAntiLambdaNegMCResdEtaSdPhiS_%d",j), "#bar{#Lambda} Neg. daughter; #it{p}_{T,V0} (GeV/#it{c}); Mass (GeV/c^2); #it{p}_{Daug} (GeV/#it{c}); #Delta#varphi*; #Delta#eta*;  Same Sign as Trigger Particle;",6,binsSplit3,xminSplit3,xmaxSplit3);
      fOutput->Add(fAntiLambdaNegMCResdEtaSdPhiS[j]);  
    }

  }

  
  // ----- Fraction of shared TPC cls
  fAntiLambdaPosDaugFracShTPCclsTrig =
    new TH3F("fAntiLambdaPosDaugFracShTPCclsTrig","#bar{#Lambda}: mass vs #it{p}_{T} vs fraction Shared TPC cls;Mass (GeV/c^2);#it{p}_{T} (GeV/#it{c});fraction Shared TPC cls",
         	     nbins,1.065,1.165,nbins,pMin,pMax,50,0,1.);
  fOutput->Add(fAntiLambdaPosDaugFracShTPCclsTrig);

  fAntiLambdaNegDaugFracShTPCclsTrig =
    new TH3F("fAntiLambdaNegDaugFracShTPCclsTrig","#bar{#Lambda}: mass vs #it{p}_{T} vs fraction Shared TPC cls;Mass (GeV/c^2);#it{p}_{T} (GeV/#it{c});fraction Shared TPC cls",
          	     nbins,1.065,1.165,nbins,pMin,pMax,50,0,1.);
 fOutput->Add(fAntiLambdaNegDaugFracShTPCclsTrig);


  //    DCA to prim vertex
  fAntiLambdaDCADaugToPrimVtx  
    = new TH3F("fAntiLambdaDCADaugToPrimVtx","#bar{#Lambda} Bckg: dca daughter vs. #it{p}_{T,l}",
	       90,0.,3.3,90,0.,3.3,nbinPtLP,pMin,ptMaxLP);
  fAntiLambdaDCADaugToPrimVtx->GetXaxis()->SetTitle("DCA Pos daug (cm)"); 
  fAntiLambdaDCADaugToPrimVtx->GetYaxis()->SetTitle("DCA Neg daug (cm)"); 
  fAntiLambdaDCADaugToPrimVtx->GetZaxis()->SetTitle("#it{p}_{T,l} (GeV/#it{c})"); 
  fOutput->Add(fAntiLambdaDCADaugToPrimVtx);

  //    Spatial Resoltuion between trigger- and asosciated- particles
  fAntiLambdaSpatialRes = new TH3F("fAntiLambdaSpatialRes","#bar{#Lambda}: Spatial resolution;#Delta#varphi (rad);trig-assoc. resolution (cm);dec. length (cm)",
				   20,-0.1,0.1,100,0.,10,2*nbins,lMin,lMax);
  fOutput->Add(fAntiLambdaSpatialRes);

  for(Int_t jj=0;jj<kNCent;jj++){
    for(Int_t k=0;k<kN1;k++){

      // Monte-Carlo level:
      if(fIsMC){
	snprintf(hNameHist,100, "fAntiLambdadPhidEtaMC_%.2f_%.2f_Cent_%.0f_%.0f",kPtBinV0[k],kPtBinV0[k+1],kBinCent[jj],kBinCent[jj+1]); 
	fAntiLambdadPhidEtaMC[jj*kN1+k] = new TH3F(hNameHist,"#bar{#Lambda} MC: #Delta#varphi vs #Delta#eta vs #it{p}_{T,l}",
						   nbinsdPhi,-TMath::PiOver2(),3*TMath::PiOver2(),
						   nbinsdEta,-1.5,1.5,
						   nbinsVtx,-10.,10.);
	fAntiLambdadPhidEtaMC[jj*kN1+k]->GetXaxis()->SetTitle("#Delta#varphi (rad)"); 
	fAntiLambdadPhidEtaMC[jj*kN1+k]->GetYaxis()->SetTitle("#Delta#eta"); 
	fAntiLambdadPhidEtaMC[jj*kN1+k]->GetZaxis()->SetTitle("Vertex Z (cm)"); 
	fOutput->Add(fAntiLambdadPhidEtaMC[jj*kN1+k]);
      }

      // Reconstruction level:
      for(Int_t ll=0;ll<kNVtxZ;ll++){
	snprintf(hNameHist,100, "fAntiLambdadPhidEtaPtL_%.2f_%.2f_Cent_%.0f_%.0f_%d",kPtBinV0[k],kPtBinV0[k+1],kBinCent[jj],kBinCent[jj+1],ll); 
	fAntiLambdadPhidEtaPtL[jj*kN1*kNVtxZ  + k*kNVtxZ + ll] = new TH3F(hNameHist,"#bar{#Lambda}: #Delta#varphi vs #Delta#eta vs #it{p}_{T,l}",
									  nbinsdPhi,-TMath::PiOver2(),3*TMath::PiOver2(),
									  nbinsdEta,-1.5,1.5,
									  nbins,1.065,1.165);
	fAntiLambdadPhidEtaPtL[jj*kN1*kNVtxZ  + k*kNVtxZ + ll]->GetXaxis()->SetTitle("#Delta#varphi (rad)"); 
	fAntiLambdadPhidEtaPtL[jj*kN1*kNVtxZ  + k*kNVtxZ + ll]->GetYaxis()->SetTitle("#Delta#eta"); 
	fAntiLambdadPhidEtaPtL[jj*kN1*kNVtxZ  + k*kNVtxZ + ll]->GetZaxis()->SetTitle("Inv. Mass");
	fOutput->Add(fAntiLambdadPhidEtaPtL[jj*kN1*kNVtxZ  + k*kNVtxZ + ll]);
      }
    }
  }

  // Correlations (side-band):
  fAntiLambdaBckgDecLength
    = new TH2F("fAntiLambdaBckgDecLength","#bar{#Lambda} Bckg: c#tau vs. #it{p}_{T,l}",
	       100,0.,25.,nbinPtLP,pMin,ptMaxLP);
  fAntiLambdaBckgDecLength->GetXaxis()->SetTitle("c#tau (cm)"); 
  fAntiLambdaBckgDecLength->GetYaxis()->SetTitle("#it{p}_{T,l} (GeV/#it{c})"); 
  fOutput->Add(fAntiLambdaBckgDecLength);
  
  fAntiLambdaBckgDCADaugToPrimVtx  
    = new TH3F("fAntiLambdaBckgDCADaugToPrimVtx","#bar{#Lambda} Bckg: dca daughter vs. #it{p}_{T,l}",
	       90,0.,3.3,90,0.,3.3,nbinPtLP,pMin,ptMaxLP);
  fAntiLambdaBckgDCADaugToPrimVtx->GetXaxis()->SetTitle("DCA Pos daug (cm)"); 
  fAntiLambdaBckgDCADaugToPrimVtx->GetYaxis()->SetTitle("DCA Neg daug (cm)"); 
  fAntiLambdaBckgDCADaugToPrimVtx->GetZaxis()->SetTitle("#it{p}_{T,l} (GeV/#it{c})"); 
  fOutput->Add(fAntiLambdaBckgDCADaugToPrimVtx);
  
  fAntiLambdaBckgEtaPhi = 
    new TH2F("fAntiLambdaBckgEtaPhi","#bar{#Lambda} Bckg: #varphi vs #eta;#varphi (rad);l_{T} (cm)",
	     nbinsPhi,0.,2.*TMath::Pi(),100,-1.,1.);
  fOutput->Add(fAntiLambdaBckgEtaPhi);
    
  fAntiLambdaBckgPhiRadio
    = new TH2F("fAntiLambdaBckgPhiRadio","#bar{#Lambda} Bckg: #varphi vs l_{T};#varphi (rad);l_{T} (cm)",
	       nbinsPhi,0.,2.*TMath::Pi(),2*nbins,lMin,lMax);
  fOutput->Add(fAntiLambdaBckgPhiRadio);


  fAntiLambdaBckgDCANegDaugToPrimVtx  
    = new TH2F("fAntiLambdaBckgDCANegDaugToPrimVtx","#bar{#Lambda} Bckg: dca NegDaughter",
	       7,-0.5,6.5,90,0.,3.3);
  fAntiLambdaBckgDCANegDaugToPrimVtx->GetXaxis()->SetTitle("MC Production"); 
  fAntiLambdaBckgDCANegDaugToPrimVtx->GetXaxis()->SetBinLabel(1,"Rec");
  fAntiLambdaBckgDCANegDaugToPrimVtx->GetXaxis()->SetBinLabel(2,"Primary");
  fAntiLambdaBckgDCANegDaugToPrimVtx->GetXaxis()->SetBinLabel(3,"V0's");
  fAntiLambdaBckgDCANegDaugToPrimVtx->GetXaxis()->SetBinLabel(4,"Cascades");
  fAntiLambdaBckgDCANegDaugToPrimVtx->GetXaxis()->SetBinLabel(5,"Gamma conv.");
  fAntiLambdaBckgDCANegDaugToPrimVtx->GetXaxis()->SetBinLabel(6,"Unidentified mother");
  fAntiLambdaBckgDCANegDaugToPrimVtx->GetXaxis()->SetBinLabel(7,"Other");
  fAntiLambdaBckgDCANegDaugToPrimVtx->GetYaxis()->SetTitle("DCA Neg Daug (cm)"); 
  fOutput->Add(fAntiLambdaBckgDCANegDaugToPrimVtx);


  fAntiLambdaBckgDCAPosDaugToPrimVtx  
    = new TH2F("fAntiLambdaBckgDCAPosDaugToPrimVtx","#bar{#Lambda} Bckg: dca PosDaughter",
	       7,-0.5,6.5,90,0.,3.3);
  fAntiLambdaBckgDCAPosDaugToPrimVtx->GetXaxis()->SetTitle("MC Production"); 
  fAntiLambdaBckgDCAPosDaugToPrimVtx->GetXaxis()->SetBinLabel(1,"Rec");
  fAntiLambdaBckgDCAPosDaugToPrimVtx->GetXaxis()->SetBinLabel(2,"Primary");
  fAntiLambdaBckgDCAPosDaugToPrimVtx->GetXaxis()->SetBinLabel(3,"V0's");
  fAntiLambdaBckgDCAPosDaugToPrimVtx->GetXaxis()->SetBinLabel(4,"Cascades");
  fAntiLambdaBckgDCAPosDaugToPrimVtx->GetXaxis()->SetBinLabel(5,"Gamma conv.");
  fAntiLambdaBckgDCAPosDaugToPrimVtx->GetXaxis()->SetBinLabel(6,"Unidentified mother");
  fAntiLambdaBckgDCAPosDaugToPrimVtx->GetXaxis()->SetBinLabel(7,"Other");
  fAntiLambdaBckgDCAPosDaugToPrimVtx->GetYaxis()->SetTitle("DCA Pos Daug (cm)"); 
  fOutput->Add(fAntiLambdaBckgDCAPosDaugToPrimVtx);


  // ============================================================= //

  // K0s in ME:  
  for(Int_t ll=0;ll<kNCent;ll++){
    for(Int_t k=0;k<kN1;k++){
      for(Int_t j=0;j<kNVtxZ;j++){
      
	snprintf(hNameHist,100,"fK0sdPhidEtaME_%.2f_%.2f_%.0f_%.0f_%d",kPtBinV0[k],kPtBinV0[k+1],kBinCent[ll],kBinCent[ll+1],j);                  
	fK0sdPhidEtaME[ll*kN1*kNVtxZ + k*kNVtxZ + j] = new TH3F(hNameHist,"K^{0}_{S}: #Delta#varphi vs #Delta#eta in ME;#Delta#varphi (rad);#Delta#eta",
								nbinsdPhi,-TMath::PiOver2(),3*TMath::PiOver2(),
								nbinsdEta,-1.5,1.5,nbins,0.398,0.598);
	fOutputME->Add(fK0sdPhidEtaME[ll*kN1*kNVtxZ + k*kNVtxZ + j]);
      }
    }
  }

  // Lambda in ME:  
  for(Int_t ll=0;ll<kNCent;ll++){
    for(Int_t k=0;k<kN1;k++){
      for(Int_t j=0;j<kNVtxZ;j++){

	snprintf(hNameHist,100,"fLambdadPhidEtaME_%.2f_%.2f_%.0lf_%.0lf_%d",kPtBinV0[k],kPtBinV0[k+1],kBinCent[ll],kBinCent[ll+1],j);
	fLambdadPhidEtaME[ll*kN1*kNVtxZ + k*kNVtxZ + j] = new TH3F(hNameHist,"#Lambda: #Delta#varphi vs #Delta#eta in ME;#Delta#varphi (rad);#Delta#eta",
								   nbinsdPhi,-TMath::PiOver2(),3*TMath::PiOver2(),
								   nbinsdEta,-1.5,1.5,nbins,1.065,1.165);
	fOutputME->Add(fLambdadPhidEtaME[ll*kN1*kNVtxZ + k*kNVtxZ + j]);
      }
    }
  }

  // AntiLambda in ME:
  for(Int_t ll=0;ll<kNCent;ll++){
    for(Int_t k=0;k<kN1;k++){
      for(Int_t j=0;j<kNVtxZ;j++){

	snprintf(hNameHist,100,"fAntiLambdadPhidEtaME_%.2f_%.2f_%.0lf_%.0lf_%d",kPtBinV0[k],kPtBinV0[k+1],kBinCent[ll],kBinCent[ll+1],j);
	fAntiLambdadPhidEtaME[ll*kN1*kNVtxZ + k*kNVtxZ + j] = new TH3F(hNameHist,"#bar{#Lambda}: #Delta#varphi vs #Delta#eta in ME;#Delta#varphi (rad);#Delta#eta",
								       nbinsdPhi,-TMath::PiOver2(),3*TMath::PiOver2(),
								       nbinsdEta,-1.5,1.5,nbins,1.065,1.165);
	fOutputME->Add(fAntiLambdadPhidEtaME[ll*kN1*kNVtxZ + k*kNVtxZ + j]);
      }
    }
  }

  
  // ============================================================= //

  if(fDoQA){

    // ----------------------------
    // Quality Assurance K0s:

    // Transverse momentum:
    //     --- signal ---
    fK0sPtPosDaug =
      new TH2F("fK0sPtPosDaug","K^{0}_{S}: #it{p}_{T};#it{p}_{T};#it{p}_{T} V0",nbins,pMin,pMax,nbins,pMin,pMax);
    fOutputQA->Add(fK0sPtPosDaug);

    fK0sPtNegDaug =
      new TH2F("fK0sPtNegDaug","K^{0}_{S}: #it{p}_{T};#it{p}_{T};#it{p}_{T} V0",nbins,pMin,pMax,nbins,pMin,pMax);
    fOutputQA->Add(fK0sPtNegDaug);

    //     --- background ---
    fK0sBckgPtPosDaug =
      new TH2F("fK0sBckgPtPosDaug","K^{0}_{S}: #it{p}_{T};#it{p}_{T};#it{p}_{T} V0",nbins,pMin,pMax,nbins,pMin,pMax);
    fOutputQA->Add(fK0sBckgPtPosDaug);

    fK0sBckgPtNegDaug =
      new TH2F("fK0sBckgPtNegDaug","K^{0}_{S}: #it{p}_{T};#it{p}_{T};#it{p}_{T} V0",nbins,pMin,pMax,nbins,pMin,pMax);
    fOutputQA->Add(fK0sBckgPtNegDaug);

    // Phi Eta
    //     --- signal ---
    fK0sPhiEtaPosDaug = 
      new TH3F("fK0sPhiEtaPosDaug","K^{0}_{S}: #varphi vs #eta Pos. Daug.;#varphi;#eta;#it{p}_{T} V0",nbinsPhi,0.,2.*TMath::Pi(),100,-1.,1.,nbins,pMin,pMax);
    fOutputQA->Add(fK0sPhiEtaPosDaug);

    fK0sPhiEtaNegDaug  = 
      new TH3F("fK0sPhiEtaNegDaug","K^{0}_{S}: #varphi vs #eta Neg. Daug.;#varphi;#eta;#it{p}_{T} V0",nbinsPhi,0.,2.*TMath::Pi(),100,-1.,1.,nbins,pMin,pMax);
    fOutputQA->Add(fK0sPhiEtaNegDaug);

    //     --- background ---
    fK0sBckgPhiEtaPosDaug = 
      new TH3F("fK0sBckgPhiEtaPosDaug","K^{0}_{S} Bckg: #varphi vs #eta Pos. Daug.;#varphi;#eta;#it{p}_{T} V0",nbinsPhi,0.,2.*TMath::Pi(),100,-1.,1.,nbins,pMin,pMax);
    fOutputQA->Add(fK0sBckgPhiEtaPosDaug);

    fK0sBckgPhiEtaNegDaug  = 
      new TH3F("fK0sBckgPhiEtaNegDaug","K^{0}_{S} Bckg: #varphi vs #eta Neg. Daug.;#varphi;#eta;#it{p}_{T} V0",nbinsPhi,0.,2.*TMath::Pi(),100,-1.,1.,nbins,pMin,pMax);
    fOutputQA->Add(fK0sBckgPhiEtaNegDaug);

    // Distance of closest approach:
    //     --- signal ---
    fK0sDCAPosDaug = 
      new TH2F("fK0sDCAPosDaug","K^{0}_{S}: dca Pos;dca;#it{p}_{T} V0",66,0.,3.3,nbins,pMin,pMax);
    fOutputQA->Add(fK0sDCAPosDaug);

    fK0sDCANegDaug =  
      new TH2F("fK0sDCANegDaug","K^{0}_{S}: dca Neg;dca;#it{p}_{T} V0",66,0.,3.3,nbins,pMin,pMax);
    fOutputQA->Add(fK0sDCANegDaug);
    
    //     --- background ---
    fK0sBckgDCAPosDaug = 
      new TH2F("fK0sBckgDCAPosDaug","K^{0}_{S} Bckg: dca Pos;dca;#it{p}_{T} V0",66,0.,3.3,nbins,pMin,pMax);
    fOutputQA->Add(fK0sBckgDCAPosDaug);

    fK0sBckgDCANegDaug =  
      new TH2F("fK0sBckgDCANegDaug","K^{0}_{S} Bckg: dca Neg;dca;#it{p}_{T} V0",66,0.,3.3,nbins,pMin,pMax);
    fOutputQA->Add(fK0sBckgDCANegDaug);

    // Decay vertex reconstruction:
    //     --- signal ---
    fK0sDecayPos  =  
      new TH3F("fK0sDecayPos","K^{0}_{S}: Position of Dec. Vtx",200,-100.,100.,200,-100.,100.,nbins,pMin,pMax);
    fK0sDecayPos->GetXaxis()->SetTitle("Pos. X"); 
    fK0sDecayPos->GetYaxis()->SetTitle("Pos. Y"); 
    fK0sDecayPos->GetZaxis()->SetTitle("#it{p}_{T} V0"); 
    fOutputQA->Add(fK0sDecayPos);

    fK0sDecayVertex  =  
      new TH2F("fK0sDecayVertex","K^{0}_{S}: decay length",100,0.,100.,nbins,pMin,pMax);
    fK0sDecayVertex->GetXaxis()->SetTitle("l_{T}"); 
    fK0sDecayVertex->GetYaxis()->SetTitle("#it{p}_{T} V0"); 
    fOutputQA->Add(fK0sDecayVertex);

    //     --- background ---
    fK0sBckgDecayPos  =  
      new TH3F("fK0sBckgDecayPos","K^{0}_{S}: Position of Dec. Vtx",200,-100.,100.,200,-100.,100.,nbins,pMin,pMax);
    fK0sBckgDecayPos->GetXaxis()->SetTitle("Pos. X"); 
    fK0sBckgDecayPos->GetYaxis()->SetTitle("Pos. Y"); 
    fK0sBckgDecayPos->GetZaxis()->SetTitle("#it{p}_{T} V0"); 
    fOutputQA->Add(fK0sBckgDecayPos);

    fK0sBckgDecayVertex  =  
      new TH2F("fK0sBckgDecayVertex","K^{0}_{S} Bckg: decay vertex",100,0.,100.,nbins,pMin,pMax);
    fK0sBckgDecayVertex->GetXaxis()->SetTitle("l_{T}"); 
    fK0sBckgDecayVertex->GetYaxis()->SetTitle("#it{p}_{T} V0"); 
    fOutputQA->Add(fK0sBckgDecayVertex);

    // Cosine of the Pointing Angle:
    //     --- signal ---
    fK0sCPA  =  
      new TH2F("fK0sCPA","K^{0}_{S}: cosine of the pointing angle",100,0.98,1.,nbins,pMin,pMax);
    fK0sCPA->GetXaxis()->SetTitle("cpa"); 
    fK0sCPA->GetYaxis()->SetTitle("#it{p}_{T} V0"); 
    fOutputQA->Add(fK0sCPA);
    //     --- background ---
    fK0sBckgCPA  =  
      new TH2F("fK0sBckgCPA","K^{0}_{S} Bckg: cosine of the pointing angle",100,0.98,1.,nbins,pMin,pMax);
    fK0sBckgCPA->GetXaxis()->SetTitle("cpa"); 
    fK0sBckgCPA->GetYaxis()->SetTitle("#it{p}_{T} V0"); 
    fOutputQA->Add(fK0sBckgCPA);

    // DCA between daughters:
    //     --- signal ---
    fK0sDCAV0Daug  =  
      new TH2F("fK0sDCAV0Daug","K^{0}_{S}: DCA daughters",60,0,1.2,nbins,pMin,pMax);
    fK0sDCAV0Daug->GetXaxis()->SetTitle("dca between daughters"); 
    fK0sDCAV0Daug->GetYaxis()->SetTitle("#it{p}_{T} V0"); 
    fOutputQA->Add(fK0sDCAV0Daug);
    //     --- background ---
    fK0sBckgDCAV0Daug  =  
      new TH2F("fK0sBckgDCAV0Daug","K^{0}_{S} Bckg: DCA daughters",60,0,1.2,nbins,pMin,pMax);
    fK0sBckgDCAV0Daug->GetXaxis()->SetTitle("dca between daughters"); 
    fK0sBckgDCAV0Daug->GetYaxis()->SetTitle("#it{p}_{T} V0"); 
    fOutputQA->Add(fK0sBckgDCAV0Daug);

    // Number of TPC clusters:
    //     --- signal ---
    fK0sNClustersTPC =  // Positive momentum to positive daugther - Negative momentum to negative daugther 
      new TH3F("fK0sNClustersTPC","K^{0}_{S};#varphi;Num. TPC Clusters; #it{p}_{T} (GeV/#it{c})",nbinsPhi,0.,2.*TMath::Pi(),181,0.5,180.5,nbins,-pMax,pMax); 
    fOutputQA->Add(fK0sNClustersTPC);
    //     --- background ---
    fK0sBckgNClustersTPC =  // Positive momentum to positive daugther - Negative momentum to negative daugther 
      new TH3F("fK0sBckgNClustersTPC","K^{0}_{S} Bckg;#varphi;Num. TPC Clusters; #it{p}_{T} (GeV/#it{c})",nbinsPhi,0.,2.*TMath::Pi(),181,0.5,180.5,nbins,-pMax,pMax); 
    fOutputQA->Add(fK0sBckgNClustersTPC);
 
    // Number of ITS clusters:
    //     --- signal ---
    fK0sNClustersITSPos = 
      new TH3F("fK0sNClustersITSPos","K^{0}_{S}: Pos. Daug;#varphi;Num. ITS Clusters;#it{p}_{T} (GeV/#it{c})",nbinsPhi,0.,2.*TMath::Pi(),7,-0.5,6.5,nbins,pMin,pMax); 
    fOutputQA->Add(fK0sNClustersITSPos);

    fK0sNClustersITSNeg = 
      new TH3F("fK0sNClustersITSNeg","K^{0}_{S}: Neg. Daug;#varphi;Num. ITS Clusters;#it{p}_{T} (GeV/#it{c})",nbinsPhi,0.,2.*TMath::Pi(),7,-0.5,6.5,nbins,pMin,pMax); 
    fOutputQA->Add(fK0sNClustersITSNeg);
    //     --- background ---
    fK0sBckgNClustersITSPos = 
      new TH3F("fK0sBckgNClustersITSPos","K^{0}_{S} Bckg: Pos. Daug;#varphi;Num. ITS Clusters;;#it{p}_{T} (GeV/#it{c})",nbinsPhi,0.,2.*TMath::Pi(),7,-0.5,6.5,nbins,pMin,pMax); 
    fOutputQA->Add(fK0sBckgNClustersITSPos);

    fK0sBckgNClustersITSNeg = 
      new TH3F("fK0sBckgNClustersITSNeg","K^{0}_{S} Bckg: Neg. Daug;#varphi;Num. ITS Clusters;;#it{p}_{T} (GeV/#it{c})",nbinsPhi,0.,2.*TMath::Pi(),7,-0.5,6.5,nbins,pMin,pMax); 
    fOutputQA->Add(fK0sBckgNClustersITSNeg);
  
    // ctau:
    //     --- signal ---
    fK0sCTau = 
      new TH2F("fK0sCTau","K^{0}_{S}: #it{c}#tau;c#tau (cm);#it{p}_{T} (GeV/#it{c})",100,0.,50.,nbins,pMin,pMax); 
    fOutputQA->Add(fK0sCTau);
    //     --- background ---
    fK0sBckgCTau = 
      new TH2F("fK0sBckgCTau","K^{0}_{S} Bckg: #it{c}#tau;c#tau (cm);#it{p}_{T} (GeV/#it{c})",100,0.,50,nbins,pMin,pMax); 
    fOutputQA->Add(fK0sBckgCTau);


    // ----------------------------
    // Quality Assurance Lambda:

    // Transverse momentum:
    //     --- signal ---
    fLambdaPtPosDaug =
      new TH2F("fLambdaPtPosDaug","#Lambda: #it{p}_{T};#it{p}_{T};#it{p}_{T} V0",nbins,pMin,pMax,nbins,pMin,pMax);
    fOutputQA->Add(fLambdaPtPosDaug);

    fLambdaPtNegDaug =
      new TH2F("fLambdaPtNegDaug","#Lambda: #it{p}_{T};#it{p}_{T};#it{p}_{T} V0",nbins,pMin,pMax,nbins,pMin,pMax);
    fOutputQA->Add(fLambdaPtNegDaug);

    //     --- background ---
    fLambdaBckgPtPosDaug =
      new TH2F("fLambdaBckgPtPosDaug","#Lambda: #it{p}_{T};#it{p}_{T};#it{p}_{T} V0",nbins,pMin,pMax,nbins,pMin,pMax);
    fOutputQA->Add(fLambdaBckgPtPosDaug);

    fLambdaBckgPtNegDaug =
      new TH2F("fLambdaBckgPtNegDaug","#Lambda: #it{p}_{T};#it{p}_{T};#it{p}_{T} V0",nbins,pMin,pMax,nbins,pMin,pMax);
    fOutputQA->Add(fLambdaBckgPtNegDaug);

    // Phi Eta
    //     --- signal ---
    fLambdaPhiEtaPosDaug = 
      new TH3F("fLambdaPhiEtaPosDaug","#Lambda: #varphi vs #eta Pos. Daug.",nbinsPhi,0.,2.*TMath::Pi(),100,-1.,1.,nbins,pMin,pMax);
    fLambdaPhiEtaPosDaug->GetXaxis()->SetTitle("#varphi"); 
    fLambdaPhiEtaPosDaug->GetYaxis()->SetTitle("#eta"); 
    fLambdaPhiEtaPosDaug->GetZaxis()->SetTitle("#it{p}_{T} V0"); 
    fOutputQA->Add(fLambdaPhiEtaPosDaug);

    fLambdaPhiEtaNegDaug  = 
      new TH3F("fLambdaPhiEtaNegDaug","#Lambda: #varphi vs #eta Neg. Daug.",nbinsPhi,0.,2.*TMath::Pi(),100,-1.,1.,nbins,pMin,pMax);
    fLambdaPhiEtaNegDaug->GetXaxis()->SetTitle("#varphi"); 
    fLambdaPhiEtaNegDaug->GetYaxis()->SetTitle("#eta"); 
    fLambdaPhiEtaNegDaug->GetZaxis()->SetTitle("#it{p}_{T} V0"); 
    fOutputQA->Add(fLambdaPhiEtaNegDaug);

    //     --- background ---
    fLambdaBckgPhiEtaPosDaug = 
      new TH3F("fLambdaBckgPhiEtaPosDaug","#Lambda: #varphi vs #eta Pos. Daug.",nbinsPhi,0.,2.*TMath::Pi(),100,-1.,1.,nbins,pMin,pMax);
    fLambdaBckgPhiEtaPosDaug->GetXaxis()->SetTitle("#varphi"); 
    fLambdaBckgPhiEtaPosDaug->GetYaxis()->SetTitle("#eta"); 
    fLambdaBckgPhiEtaPosDaug->GetZaxis()->SetTitle("#it{p}_{T} V0"); 
    fOutputQA->Add(fLambdaBckgPhiEtaPosDaug);

    fLambdaBckgPhiEtaNegDaug  = 
      new TH3F("fLambdaBckgPhiEtaNegDaug","#Lambda: #varphi vs #eta Neg. Daug.",nbinsPhi,0.,2.*TMath::Pi(),100,-1.,1.,nbins,pMin,pMax);
    fLambdaBckgPhiEtaNegDaug->GetXaxis()->SetTitle("#varphi"); 
    fLambdaBckgPhiEtaNegDaug->GetYaxis()->SetTitle("#eta"); 
    fLambdaBckgPhiEtaNegDaug->GetZaxis()->SetTitle("#it{p}_{T} V0"); 
    fOutputQA->Add(fLambdaBckgPhiEtaNegDaug);

    // Distance of closest approach
    //     --- signal ---
    fLambdaDCAPosDaug = 
      new TH2F("fLambdaDCAPosDaug","#Lambda: dca Pos",66,0.,3.3,nbins,pMin,pMax);
    fLambdaDCAPosDaug->GetXaxis()->SetTitle("dca"); 
    fLambdaDCAPosDaug->GetYaxis()->SetTitle("#it{p}_{T} V0"); 
    fOutputQA->Add(fLambdaDCAPosDaug);

    fLambdaDCANegDaug =  
      new TH2F("fLambdaDCANegDaug","#Lambda: dca Neg",66,0.,3.3,nbins,pMin,pMax);
    fLambdaDCANegDaug->GetXaxis()->SetTitle("dca"); 
    fLambdaDCANegDaug->GetYaxis()->SetTitle("#it{p}_{T} V0"); 
    fOutputQA->Add(fLambdaDCANegDaug);
    
    //     --- background ---
    fLambdaBckgDCAPosDaug = 
      new TH2F("fLambdaBckgDCAPosDaug","#Lambda Bckg: dca Pos",66,0.,3.3,nbins,pMin,pMax);
    fLambdaBckgDCAPosDaug->GetXaxis()->SetTitle("dca"); 
    fLambdaBckgDCAPosDaug->GetYaxis()->SetTitle("#it{p}_{T} V0"); 
    fOutputQA->Add(fLambdaBckgDCAPosDaug);

    fLambdaBckgDCANegDaug =  
      new TH2F("fLambdaBckgDCANegDaug","#Lambda Bckg: dca Neg",66,0.,3.3,nbins,pMin,pMax);
    fLambdaBckgDCANegDaug->GetXaxis()->SetTitle("dca"); 
    fLambdaBckgDCANegDaug->GetYaxis()->SetTitle("#it{p}_{T} V0"); 
    fOutputQA->Add(fLambdaBckgDCANegDaug);


    // Decay vertex reconstruction
    //     --- signal ---
    fLambdaDecayPos  =  
      new TH3F("fLambdaDecayPos","#Lambda: Position of Dec. Vtx",200,-100.,100.,200,-100.,100.,nbins,pMin,pMax);
    fLambdaDecayPos->GetXaxis()->SetTitle("Pos. X"); 
    fLambdaDecayPos->GetYaxis()->SetTitle("Pos. Y"); 
    fLambdaDecayPos->GetZaxis()->SetTitle("#it{p}_{T} V0"); 
    fOutputQA->Add(fLambdaDecayPos);

    fLambdaDecayVertex  =  
      new TH2F("fLambdaDecayVertex","#Lambda: decay length",100,0.,100.,nbins,pMin,pMax);
    fLambdaDecayVertex->GetXaxis()->SetTitle("l_{T}"); 
    fLambdaDecayVertex->GetYaxis()->SetTitle("#it{p}_{T} V0"); 
    fOutputQA->Add(fLambdaDecayVertex);

    //     --- background ---
    fLambdaBckgDecayPos  =  
      new TH3F("fLambdaBckgDecayPos","#Lambda Bckg: Position of Dec. Vtx",200,-100.,100.,200,-100.,100.,nbins,pMin,pMax);
    fLambdaBckgDecayPos->GetXaxis()->SetTitle("Pos. X"); 
    fLambdaBckgDecayPos->GetYaxis()->SetTitle("Pos. Y"); 
    fLambdaBckgDecayPos->GetZaxis()->SetTitle("#it{p}_{T} V0"); 
    fOutputQA->Add(fLambdaBckgDecayPos);

    fLambdaBckgDecayVertex  =  
      new TH2F("fLambdaBckgDecayVertex","#Lambda Bckg: decay length",100,0.,100.,nbins,pMin,pMax);
    fLambdaBckgDecayVertex->GetXaxis()->SetTitle("l_{T}"); 
    fLambdaBckgDecayVertex->GetYaxis()->SetTitle("#it{p}_{T} V0"); 
    fOutputQA->Add(fLambdaBckgDecayVertex);

    // Cosine of the Pointing Angle
    //     --- signal ---
    fLambdaCPA  =  
      new TH2F("fLambdaCPA","#Lambda: cosine of the pointing angle",100,0.98,1.,nbins,pMin,pMax);
    fLambdaCPA->GetXaxis()->SetTitle("cpa"); 
    fLambdaCPA->GetYaxis()->SetTitle("#it{p}_{T} V0"); 
    fOutputQA->Add(fLambdaCPA);
    //     --- background ---
    fLambdaBckgCPA  =  
      new TH2F("fLambdaBckgCPA","#Lambda Bckg: cosine of the pointing angle",100,0.98,1.,nbins,pMin,pMax);
    fLambdaBckgCPA->GetXaxis()->SetTitle("cpa"); 
    fLambdaBckgCPA->GetYaxis()->SetTitle("#it{p}_{T} V0"); 
    fOutputQA->Add(fLambdaBckgCPA);

    // DCA between daughters
    //     --- signal ---
    fLambdaDCAV0Daug  =  
      new TH2F("fLambdaDCAV0Daug","#Lambda: DCA daughters",60,0,1.2,nbins,pMin,pMax);
    fLambdaDCAV0Daug->GetXaxis()->SetTitle("dca between daughters"); 
    fLambdaDCAV0Daug->GetYaxis()->SetTitle("#it{p}_{T} V0"); 
    fOutputQA->Add(fLambdaDCAV0Daug);
    //     --- background ---
    fLambdaBckgDCAV0Daug  =  
      new TH2F("fLambdaBckgDCAV0Daug","#Lambda Bckg: DCA daughters",60,0,1.2,nbins,pMin,pMax);
    fLambdaBckgDCAV0Daug->GetXaxis()->SetTitle("dca between daughters"); 
    fLambdaBckgDCAV0Daug->GetYaxis()->SetTitle("#it{p}_{T} V0"); 
    fOutputQA->Add(fLambdaBckgDCAV0Daug);
  
    // Number of TPC clusters:
    //     --- signal ---
    fLambdaNClustersTPC = 
      new TH3F("fLambdaNClustersTPC","#Lambda;#varphi;Num. TPC Clusters;#it{p}_{T} (GeV/#it{c})",nbinsPhi,0.,2.*TMath::Pi(),181,0.5,180.5,nbins,-pMax,pMax); 
    fOutputQA->Add(fLambdaNClustersTPC);
    //     --- background ---
    fLambdaBckgNClustersTPC = 
      new TH3F("fLambdaBckgNClustersTPC","#Lambda Bckg;#varphi;Num. TPC Clusters;#it{p}_{T} (GeV/#it{c})",nbinsPhi,0.,2.*TMath::Pi(),181,0.5,180.5,nbins,-pMax,pMax); 
    fOutputQA->Add(fLambdaBckgNClustersTPC);
 
    // Number of ITS clusters:
    //     --- signal ---
    fLambdaNClustersITSPos = 
      new TH3F("fLambdaNClustersITSPos","#Lambda: Pos. Daug;#varphi;Num. ITS Clusters;#it{p}_{T} (GeV/#it{c})",nbinsPhi,0.,2.*TMath::Pi(),7,-0.5,6.5,nbins,pMin,pMax); 
    fOutputQA->Add(fLambdaNClustersITSPos);

    fLambdaNClustersITSNeg = 
      new TH3F("fLambdaNClustersITSNeg","#Lambda: Neg. Daug;#varphi;Num. ITS Clusters;#it{p}_{T} (GeV/#it{c})",nbinsPhi,0.,2.*TMath::Pi(),7,-0.5,6.5,nbins,pMin,pMax); 
    fOutputQA->Add(fLambdaNClustersITSNeg);
    //     --- background ---
    fLambdaBckgNClustersITSPos = 
      new TH3F("fLambdaBckgNClustersITSPos","#Lambda Bckg: Pos. Daug;#varphi;Num. ITS Clusters;;#it{p}_{T} (GeV/#it{c})",nbinsPhi,0.,2.*TMath::Pi(),7,-0.5,6.5,nbins,pMin,pMax); 
    fOutputQA->Add(fLambdaBckgNClustersITSPos);

    fLambdaBckgNClustersITSNeg = 
      new TH3F("fLambdaBckgNClustersITSNeg","#Lambda Bckg: Neg. Daug;#varphi;Num. ITS Clusters;;#it{p}_{T} (GeV/#it{c})",nbinsPhi,0.,2.*TMath::Pi(),7,-0.5,6.5,nbins,pMin,pMax); 
    fOutputQA->Add(fLambdaBckgNClustersITSNeg);


    // ctau:
    //     --- signal ---
    fLambdaCTau = 
      new TH2F("fLambdaCTau","#Lambda: #it{c}#tau;c#tau (cm);#it{p}_{T} (GeV/#it{c})",100,0.,50.,nbins,pMin,pMax); 
    fOutputQA->Add(fLambdaCTau);
    //     --- background ---
    fLambdaBckgCTau = 
      new TH2F("fLambdaBckgCTau","#Lambda Bckg: #it{c}#tau;c#tau (cm);#it{p}_{T} (GeV/#it{c})",100,0.,50.,nbins,pMin,pMax); 
    fOutputQA->Add(fLambdaBckgCTau);


    // ----------------------------
    // Quality Assurance AntiLambda:
    // Transverse momentum:
    //     --- signal ---
    fAntiLambdaPtPosDaug =
      new TH2F("fAntiLambdaPtPosDaug","#bar{#Lambda}: #it{p}_{T};#it{p}_{T};#it{p}_{T} V0",nbins,pMin,pMax,nbins,pMin,pMax);
    fOutputQA->Add(fAntiLambdaPtPosDaug);

    fAntiLambdaPtNegDaug =
      new TH2F("fAntiLambdaPtNegDaug","#bar{#Lambda}: #it{p}_{T};#it{p}_{T};#it{p}_{T} V0",nbins,pMin,pMax,nbins,pMin,pMax);
    fOutputQA->Add(fAntiLambdaPtNegDaug);

    //     --- background ---
    fAntiLambdaBckgPtPosDaug =
      new TH2F("fAntiLambdaBckgPtPosDaug","#bar{#Lambda}: #it{p}_{T};#it{p}_{T};#it{p}_{T} V0",nbins,pMin,pMax,nbins,pMin,pMax);
    fOutputQA->Add(fAntiLambdaBckgPtPosDaug);

    fAntiLambdaBckgPtNegDaug =
      new TH2F("fAntiLambdaBckgPtNegDaug","#bar{#Lambda}: #it{p}_{T};#it{p}_{T};#it{p}_{T} V0",nbins,pMin,pMax,nbins,pMin,pMax);
    fOutputQA->Add(fAntiLambdaBckgPtNegDaug);

    // Phi Eta
    //     --- signal ---
    fAntiLambdaPhiEtaPosDaug = 
      new TH3F("fAntiLambdaPhiEtaPosDaug","#bar{#Lambda}: #varphi vs #eta Pos. Daug.",nbinsPhi,0.,2.*TMath::Pi(),100,-1.,1.,nbins,pMin,pMax);
    fAntiLambdaPhiEtaPosDaug->GetXaxis()->SetTitle("#varphi"); 
    fAntiLambdaPhiEtaPosDaug->GetYaxis()->SetTitle("#eta"); 
    fAntiLambdaPhiEtaPosDaug->GetZaxis()->SetTitle("#it{p}_{T} V0"); 
    fOutputQA->Add(fAntiLambdaPhiEtaPosDaug);

    fAntiLambdaPhiEtaNegDaug  = 
      new TH3F("fAntiLambdaPhiEtaNegDaug","#bar{#Lambda}: #varphi vs #eta Neg. Daug.",nbinsPhi,0.,2.*TMath::Pi(),100,-1.,1.,nbins,pMin,pMax);
    fAntiLambdaPhiEtaNegDaug->GetXaxis()->SetTitle("#varphi"); 
    fAntiLambdaPhiEtaNegDaug->GetYaxis()->SetTitle("#eta"); 
    fAntiLambdaPhiEtaNegDaug->GetZaxis()->SetTitle("#it{p}_{T} V0"); 
    fOutputQA->Add(fAntiLambdaPhiEtaNegDaug);

    //     --- background ---
    fAntiLambdaBckgPhiEtaPosDaug = 
      new TH3F("fAntiLambdaBckgPhiEtaPosDaug","#bar{#Lambda}: #varphi vs #eta Pos. Daug.",nbinsPhi,0.,2.*TMath::Pi(),100,-1.,1.,nbins,pMin,pMax);
    fAntiLambdaBckgPhiEtaPosDaug->GetXaxis()->SetTitle("#varphi"); 
    fAntiLambdaBckgPhiEtaPosDaug->GetYaxis()->SetTitle("#eta"); 
    fAntiLambdaBckgPhiEtaPosDaug->GetZaxis()->SetTitle("#it{p}_{T} V0"); 
    fOutputQA->Add(fAntiLambdaBckgPhiEtaPosDaug);

    fAntiLambdaBckgPhiEtaNegDaug  = 
      new TH3F("fAntiLambdaBckgPhiEtaNegDaug","#bar{#Lambda}: #varphi vs #eta Neg. Daug.",nbinsPhi,0.,2.*TMath::Pi(),100,-1.,1.,nbins,pMin,pMax);
    fAntiLambdaBckgPhiEtaNegDaug->GetXaxis()->SetTitle("#varphi"); 
    fAntiLambdaBckgPhiEtaNegDaug->GetYaxis()->SetTitle("#eta"); 
    fAntiLambdaBckgPhiEtaNegDaug->GetZaxis()->SetTitle("#it{p}_{T} V0"); 
    fOutputQA->Add(fAntiLambdaBckgPhiEtaNegDaug);

    // Distance of closest approach
    //     --- signal ---
    fAntiLambdaDCAPosDaug = 
      new TH2F("fAntiLambdaDCAPosDaug","#bar{#Lambda}: dca Pos",66,0.,3.3,nbins,pMin,pMax);
    fAntiLambdaDCAPosDaug->GetXaxis()->SetTitle("dca"); 
    fAntiLambdaDCAPosDaug->GetYaxis()->SetTitle("#it{p}_{T} V0"); 
    fOutputQA->Add(fAntiLambdaDCAPosDaug);

    fAntiLambdaDCANegDaug =  
      new TH2F("fAntiLambdaDCANegDaug","#bar{#Lambda}: dca Neg",66,0.,3.3,nbins,pMin,pMax);
    fAntiLambdaDCANegDaug->GetXaxis()->SetTitle("dca"); 
    fAntiLambdaDCANegDaug->GetYaxis()->SetTitle("#it{p}_{T} V0"); 
    fOutputQA->Add(fAntiLambdaDCANegDaug);
    
    //     --- background ---
    fAntiLambdaBckgDCAPosDaug = 
      new TH2F("fAntiLambdaBckgDCAPosDaug","#bar{#Lambda} Bckg: dca Pos",66,0.,3.3,nbins,pMin,pMax);
    fAntiLambdaBckgDCAPosDaug->GetXaxis()->SetTitle("dca"); 
    fAntiLambdaBckgDCAPosDaug->GetYaxis()->SetTitle("#it{p}_{T} V0"); 
    fOutputQA->Add(fAntiLambdaBckgDCAPosDaug);

    fAntiLambdaBckgDCANegDaug =  
      new TH2F("fAntiLambdaBckgDCANegDaug","#bar{#Lambda} Bckg: dca Neg",66,0.,3.3,nbins,pMin,pMax);
    fAntiLambdaBckgDCANegDaug->GetXaxis()->SetTitle("dca"); 
    fAntiLambdaBckgDCANegDaug->GetYaxis()->SetTitle("#it{p}_{T} V0"); 
    fOutputQA->Add(fAntiLambdaBckgDCANegDaug);

    // Decay vertex reconstruction
    //     --- signal ---
    fAntiLambdaDecayPos  =  
      new TH3F("fAntiLambdaDecayPos","#bar{#Lambda}: Position of Dec. Vtx",200,-100.,100.,200,-100.,100.,nbins,pMin,pMax);
    fAntiLambdaDecayPos->GetXaxis()->SetTitle("Pos. X"); 
    fAntiLambdaDecayPos->GetYaxis()->SetTitle("Pos. Y"); 
    fAntiLambdaDecayPos->GetZaxis()->SetTitle("#it{p}_{T} V0"); 
    fOutputQA->Add(fAntiLambdaDecayPos);

    fAntiLambdaDecayVertex  =  
      new TH2F("fAntiLambdaDecayVertex","#bar{#Lambda}: decay length",100,0.,100.,nbins,pMin,pMax);
    fAntiLambdaDecayVertex->GetXaxis()->SetTitle("l_{T}"); 
    fAntiLambdaDecayVertex->GetYaxis()->SetTitle("#it{p}_{T} V0"); 
    fOutputQA->Add(fAntiLambdaDecayVertex);

    //     --- background ---
    fAntiLambdaBckgDecayPos  =  
      new TH3F("fAntiLambdaBckgDecayPos","#bar{#Lambda} Bckg: Position of Dec. Vtx",200,-100.,100.,200,-100.,100.,nbins,pMin,pMax);
    fAntiLambdaBckgDecayPos->GetXaxis()->SetTitle("Pos. X"); 
    fAntiLambdaBckgDecayPos->GetYaxis()->SetTitle("Pos. Y"); 
    fAntiLambdaBckgDecayPos->GetZaxis()->SetTitle("#it{p}_{T} V0"); 
    fOutputQA->Add(fAntiLambdaBckgDecayPos);

    fAntiLambdaBckgDecayVertex  =  
      new TH2F("fAntiLambdaBckgDecayVertex","#bar{#Lambda} Bckg: decay length",100,0.,100.,nbins,pMin,pMax);
    fAntiLambdaBckgDecayVertex->GetXaxis()->SetTitle("l_{T}"); 
    fAntiLambdaBckgDecayVertex->GetYaxis()->SetTitle("#it{p}_{T} V0"); 
    fOutputQA->Add(fAntiLambdaBckgDecayVertex);

    // Cosine of the Pointing Angle
    //     --- signal ---
    fAntiLambdaCPA  =  
      new TH2F("fAntiLambdaCPA","#bar{#Lambda}: cosine of the pointing angle",100,0.98,1.,nbins,pMin,pMax);
    fAntiLambdaCPA->GetXaxis()->SetTitle("cpa"); 
    fAntiLambdaCPA->GetYaxis()->SetTitle("#it{p}_{T} V0"); 
    fOutputQA->Add(fAntiLambdaCPA);
    //     --- background ---
    fAntiLambdaBckgCPA  =  
      new TH2F("fAntiLambdaBckgCPA","#bar{#Lambda} Bckg: cosine of the pointing angle",100,0.98,1.,nbins,pMin,pMax);
    fAntiLambdaBckgCPA->GetXaxis()->SetTitle("cpa"); 
    fAntiLambdaBckgCPA->GetYaxis()->SetTitle("#it{p}_{T} V0"); 
    fOutputQA->Add(fAntiLambdaBckgCPA);

    // DCA between daughters
    //     --- signal ---
    fAntiLambdaDCAV0Daug  =  
      new TH2F("fAntiLambdaDCAV0Daug","#bar{#Lambda}: DCA daughters",60,0,1.2,nbins,pMin,pMax);
    fAntiLambdaDCAV0Daug->GetXaxis()->SetTitle("dca between daughters"); 
    fAntiLambdaDCAV0Daug->GetYaxis()->SetTitle("#it{p}_{T} V0"); 
    fOutputQA->Add(fAntiLambdaDCAV0Daug);
    //     --- background ---
    fAntiLambdaBckgDCAV0Daug  =  
      new TH2F("fAntiLambdaBckgDCAV0Daug","#bar{#Lambda} Bckg: DCA daughters",60,0,1.2,nbins,pMin,pMax);
    fAntiLambdaBckgDCAV0Daug->GetXaxis()->SetTitle("dca between daughters"); 
    fAntiLambdaBckgDCAV0Daug->GetYaxis()->SetTitle("#it{p}_{T} V0"); 
    fOutputQA->Add(fAntiLambdaBckgDCAV0Daug);

    // Number of TPC clusters:
    //     --- signal ---
    fAntiLambdaNClustersTPC = 
      new TH3F("fAntiLambdaNClustersTPC","#bar{#Lambda};#varphi;Num. TPC Clusters;#it{p}_{T} (GeV/#it{c})",nbinsPhi,0.,2.*TMath::Pi(),181,0.5,180.5,nbins,-pMax,pMax); 
    fOutputQA->Add(fAntiLambdaNClustersTPC);
    //     --- background ---
    fAntiLambdaBckgNClustersTPC = 
      new TH3F("fAntiLambdaBckgNClustersTPC","#bar{#Lambda} Bckg;#varphi;Num. TPC Clusters;#it{p}_{T} (GeV/#it{c})",nbinsPhi,0.,2.*TMath::Pi(),181,0.5,180.5,nbins,-pMax,pMax); 
    fOutputQA->Add(fAntiLambdaBckgNClustersTPC);
 
    // Number of ITS clusters:
    //     --- signal ---
    fAntiLambdaNClustersITSPos = 
      new TH3F("fAntiLambdaNClustersITSPos","#bar{#Lambda}: Pos. Daug;#varphi;Num. ITS Clusters;#it{p}_{T} (GeV/#it{c})",nbinsPhi,0.,2.*TMath::Pi(),7,-0.5,6.5,nbins,pMin,pMax); 
    fOutputQA->Add(fAntiLambdaNClustersITSPos);

    fAntiLambdaNClustersITSNeg = 
      new TH3F("fAntiLambdaNClustersITSNeg","#bar{#Lambda}: Neg. Daug;#varphi;Num. ITS Clusters;#it{p}_{T} (GeV/#it{c})",nbinsPhi,0.,2.*TMath::Pi(),7,-0.5,6.5,nbins,pMin,pMax); 
    fOutputQA->Add(fAntiLambdaNClustersITSNeg);
    //     --- background ---
    fAntiLambdaBckgNClustersITSPos = 
      new TH3F("fAntiLambdaBckgNClustersITSPos","#bar{#Lambda} Bckg: Pos. Daug;#varphi;Num. ITS Clusters;;#it{p}_{T} (GeV/#it{c})",nbinsPhi,0.,2.*TMath::Pi(),7,-0.5,6.5,nbins,pMin,pMax); 
    fOutputQA->Add(fAntiLambdaBckgNClustersITSPos);

    fAntiLambdaBckgNClustersITSNeg = 
      new TH3F("fAntiLambdaBckgNClustersITSNeg","#bar{#Lambda} Bckg: Neg. Daug;#varphi;Num. ITS Clusters;;#it{p}_{T} (GeV/#it{c})",nbinsPhi,0.,2.*TMath::Pi(),7,-0.5,6.5,nbins,pMin,pMax); 
    fOutputQA->Add(fAntiLambdaBckgNClustersITSNeg);

  // ctau:
    //     --- signal ---
    fAntiLambdaCTau = 
      new TH2F("fAntiLambdaCTau","#bar{#Lambda}: #it{c}#tau;c#tau (cm);#it{p}_{T} (GeV/#it{c})",100,0.,50.,nbins,pMin,pMax); 
    fOutputQA->Add(fAntiLambdaCTau);
    //     --- background ---
    fAntiLambdaBckgCTau = 
      new TH2F("fAntiLambdaBckgCTau","#bar{#Lambda} Bckg: #it{c}#tau;c#tau (cm);#it{p}_{T} (GeV/#it{c})",100,0.,50.,nbins,pMin,pMax); 
    fOutputQA->Add(fAntiLambdaBckgCTau);


  }

  // ============================================================= //
  
  PostData(1, fOutput);
  PostData(2, fOutputME);
  PostData(3, fOutputQA);
  
}

//___________________________________________________________________________________________

static Int_t VtxBin(Double_t vtx)
{
  // Bin in vertez position Z
  Int_t bin = -1;
  for(Int_t i=0;i<kNVtxZ;i++)
    if ( (vtx>=kBinVtxZ[i]) && (vtx<kBinVtxZ[i+1]) )
      bin = i;

  return bin;

}

//___________________________________________________________________________________________

static Int_t PtBin(Double_t pt)
{
  // Bin in pt
  Int_t bin = -1;
  for(Int_t i=0;i<kN1;i++)
    if ( (pt>=kPtBinV0[i]) && (pt<kPtBinV0[i+1]) )
      bin = i;

  return bin;

}

//___________________________________________________________________________________________

static Int_t CentBin(Double_t cent)
{
  // Bin in pt
  Int_t bin = -1;
  for(Int_t i=0;i<kNCent;i++)
    if ( (cent>=kBinCent[i]) && (cent<kBinCent[i+1]) )
      bin = i;

  return bin;

}

//___________________________________________________________________________________________

Bool_t AliAnalysisTaskLambdaOverK0sJets::AcceptTrack(const AliAODTrack *t) 
{
  // Track criteria for primaries particles 
  //if(fTriggerFB!=128 && fTriggerFB!=272) return kFALSE; 
  
  if (TMath::Abs(t->Eta())>0.8 )  return kFALSE; 
  //if (!(t->TestFilterMask(1<<7))) return kFALSE; 
  //if( !(t->TestFilterBit(272)) )  return kFALSE;
  if( !(t->TestFilterBit(fTriggerFB)) )  return kFALSE;
 
  Float_t nCrossedRowsTPC = t->GetTPCClusterInfo(2,1); 
  if (nCrossedRowsTPC < fTriggerNCls) return kFALSE;
  
   // Point in the SPD
  Int_t SPDHits = t->HasPointOnITSLayer(0) + t->HasPointOnITSLayer(1);
  if( SPDHits )
    fTriggerWiSPDHit->Fill(1.5);
  
  return kTRUE;   
}

//___________________________________________________________________________________________

Bool_t AliAnalysisTaskLambdaOverK0sJets::AcceptTrackV0(const AliAODTrack *t) 
{ 
  // Track criteria for daughter particles of V0 candidate 
  if (!t->IsOn(AliAODTrack::kTPCrefit)) return kFALSE;
  Float_t nCrossedRowsTPC = t->GetTPCClusterInfo(2,1); 
  if (nCrossedRowsTPC<fDaugNClsTPC) return kFALSE;

  return kTRUE;   
}

//___________________________________________________________________________________________

Bool_t AliAnalysisTaskLambdaOverK0sJets::AcceptV0(AliAODVertex *vtx, const AliAODv0 *v1) 
{ 
  // Selection for accepting V0 candidates 

  if (v1->GetOnFlyStatus()) return kFALSE;
  
  //if (v1->Pt() < pMin) return kFALSE; ***
  
  const AliAODTrack *ntrack1=(AliAODTrack *)v1->GetDaughter(1);
  const AliAODTrack *ptrack1=(AliAODTrack *)v1->GetDaughter(0);
    
  if( !ntrack1 || !ptrack1 ) return kFALSE;
  if( !AcceptTrackV0(ntrack1) ) return kFALSE;
  if( !AcceptTrackV0(ptrack1) ) return kFALSE;
  
  if( ntrack1->Charge() == ptrack1->Charge()) 
    return kFALSE;

  // Daughters: pseudo-rapidity cut
  if ( TMath::Abs(ntrack1->Eta()) > fMaxEtaDaughter  ||
       TMath::Abs(ptrack1->Eta()) > fMaxEtaDaughter  )
    return kFALSE;

  // Daughters: transverse momentum cut
  if ( ( ntrack1->Pt() < fMinPtDaughter ) || 
       ( ptrack1->Pt() < fMinPtDaughter )  ) 
    return kFALSE;
  
  // Daughters: Impact parameter of daughter to prim vtx
  Float_t xy = v1->DcaNegToPrimVertex();
  if (TMath::Abs(xy)<fDCAToPrimVtx) return kFALSE;
  xy = v1->DcaPosToPrimVertex();
  if (TMath::Abs(xy)<fDCAToPrimVtx) return kFALSE;

  // Daughters: DCA
  Float_t dca = v1->DcaV0Daughters();
  if (dca>fMaxDCADaughter) return kFALSE;

  // V0: Cosine of the pointing angle
  Float_t cpa=v1->CosPointingAngle(vtx);
  if (cpa<fMinCPA) return kFALSE;

  // V0: Fiducial volume
  Double_t xyz[3]; v1->GetSecondaryVtx(xyz);
  Float_t r2=xyz[0]*xyz[0] + xyz[1]*xyz[1];
  if (r2<5.*5.) return kFALSE;
  if (r2>lMax*lMax) return kFALSE;

  return kTRUE;
}

//___________________________________________________________________________________________

Float_t AliAnalysisTaskLambdaOverK0sJets::GetMultiplicity(){
  
  Float_t mult=0;
  Int_t nTrk= fAOD->GetNumberOfTracks();
  for (Int_t i=0; i<nTrk; i++) {
    const AliAODTrack *t = dynamic_cast<const AliAODTrack*>(fAOD->GetTrack(i));
    if(!t) AliFatal("Not a standard AOD");
    if(!AcceptTrack(t)) continue;
    mult++;
  }

  return mult;
}

//___________________________________________________________________________________________

static Float_t dPHI(Float_t phi1, Float_t phi2) 
{ 
  // Calculate the phi difference between two tracks  
  Float_t deltaPhi = phi1 - phi2;
  
  if (deltaPhi<-TMath::PiOver2())    deltaPhi = deltaPhi + 2*(TMath::Pi());
  if (deltaPhi>(3*TMath::PiOver2()))  deltaPhi = deltaPhi - 2*(TMath::Pi());
  return deltaPhi;
}

//___________________________________________________________________________________________

static Float_t MyRapidity(Float_t rE, Float_t rPz)
{ 
  // Local method for rapidity
  return 0.5*TMath::Log((rE+rPz)/(rE-rPz+1.e-13));
} 

//___________________________________________________________________________________________

static Int_t SameTrack(const AliAODTrack *trig, const AliAODTrack *daug)
{ 
  // Local method to compaire the momentum between two tracks

  //double const kEpsilon = 0.01;
  Int_t    isSamePt = 0;

  /*
    Float_t p[3];     trk->GetPxPyPz(p);
    Float_t pNegTrk[3]; nTrk->GetPxPyPz(pNegTrk);
    Float_t pPosTrk[3]; pTrk->GetPxPyPz(pPosTrk);
  
    if( (  fabs(p[0]-pNegTrk[0])<kEpsilon && 
    fabs(p[1]-pNegTrk[1])<kEpsilon && 
    fabs(p[2]-pNegTrk[2])<kEpsilon ) 
    isSamePt = 1;
  */
    
  if(  (TMath::Abs(daug->GetID())+1)==(TMath::Abs(trig->GetID()))  )
    isSamePt = 1;
  

  return isSamePt;

}

//___________________________________________________________________________________________

static Int_t SameLabel(const AliAODTrack *trig, const AliAODTrack *daug)
{ 
  // Compaire the the label value that points back to the Monte Carlo production
  //cout << trig->GetLabel() << "         " << daug->GetLabel() << endl;

  if(  TMath::Abs(trig->GetLabel() ) == 
       TMath::Abs(daug->GetLabel() )  )
    return 1.0;
  
  return 0.;

}

//___________________________________________________________________________________________

static Float_t SpatialResolution(Float_t p1x,Float_t p1y,Float_t p2x,Float_t p2y,Float_t dist)
{
  // Obtains the spacial resolution between trigger and V0
  // within a distance in (deltaPhi,deltaEta) < 0.1

  Float_t res = -100.;

  res = TMath::Sqrt( p1x*p1x + p1y*p1y )*TMath::Sqrt( p2x*p2x + p2y*p2y );
  res = (p1x*p2x + p1y*p2y)/res;

  res = TMath::ACos(res);
  
  return res = TMath::Sin(res)*dist;
 
}

//___________________________________________________________________________________________

static Float_t GetDPhiStar(Float_t phi1, Float_t pt1, Float_t charge1, Float_t phi2, Float_t pt2, Float_t charge2, Float_t radius, Float_t bSign)
{
  //
  // calculates dphistar
  //

  Float_t dphistar = phi1 - phi2 - charge1 * bSign * TMath::ASin(0.075 * radius / pt1) + charge2 * bSign * TMath::ASin(0.075 * radius / pt2);
  static const Double_t kPi = TMath::Pi();

  // circularity
  if (dphistar > kPi)
    dphistar = kPi * 2 - dphistar;
  if (dphistar < -kPi)
    dphistar = -kPi * 2 - dphistar;
  if (dphistar > kPi) // might look funny but is needed
    dphistar = kPi * 2 - dphistar;

  return dphistar;

}


//___________________________________________________________________________________________

static Float_t TwoTrackEfficiencyCut(Float_t phi1, Float_t eta1, Float_t pt1, Float_t charge1, Float_t phi2, Float_t eta2, Float_t pt2, Float_t charge2, Float_t bSign)
{
  // Code taken from the HBT analysis to reject the track splitting
  // It was modified to provide only the value of kDphiStarMax
  // and a possible rejection in the kDphiStarMean

  Float_t kRadiousDphiStarMax = -0.0005;
  Float_t deta = eta1 - eta2;
  Float_t twoTrackEfficiencyCutValue = 0.02;

  // optimization
  if (TMath::Abs(deta) < twoTrackEfficiencyCutValue * 2.5 * 3) {

    // check first boundaries to see if is worth to loop and find the minimum
    Float_t dphistar1 = GetDPhiStar(phi1, pt1, charge1, phi2, pt2, charge2, 0.8, bSign);
    Float_t dphistar2 = GetDPhiStar(phi1, pt1, charge1, phi2, pt2, charge2, 2.5, bSign);

    const Float_t kLimit = twoTrackEfficiencyCutValue * 3;

    //Float_t dphistarminabs = 1e5;
    //Float_t dphistarmin = 1e5;

    if (TMath::Abs(dphistar1) < kLimit || TMath::Abs(dphistar2) < kLimit || dphistar1 * dphistar2 < 0){
  
      kRadiousDphiStarMax = 0;
      //kDphiStarMean = 0; 
      //Int_t i=0;

      for (Double_t rad=0.8; rad<2.51; rad+=0.01){

	if ( TMath::Abs(0.075 * rad / pt2)>1 ) break;

	Float_t dphistar = GetDPhiStar(phi1, pt1, charge1, phi2, pt2, charge2, rad, bSign);
	Float_t dphistarabs = TMath::Abs(dphistar);

	if( ( (dphistarabs*rad) > kRadiousDphiStarMax) && ( TMath::Abs(deta) < twoTrackEfficiencyCutValue ) ){
	  kRadiousDphiStarMax = dphistarabs*rad;
	}

	//kDphiStarMean += dphistarabs;
	//i++;

      }
      
      //kDphiStarMean = kDphiStarMean/i;
      /*if (TMath::Abs(deta) < twoTrackEfficiencyCutValue && kDphiStarMean < twoTrackEfficiencyCutValue ){
	return kFALSE;
	}*/
      
    } // End selection in dphistar
    
  } // End dEta value


  return kRadiousDphiStarMax;

}

//___________________________________________________________________________________________
/*
static Bool_t GoodTPCSharedMap(const AliAODTrack *track)
{
  // Rejects tracks with shared clusters after filling a control histogram
  // This overload is used for primaries
 
  // Get the shared maps
  const TBits sharedMap = track->GetTPCSharedMap();
  // Fill a control histogram
  //fPriHistShare->Fill(sharedMap.CountBits());
  // Reject shared clusters
  if((sharedMap.CountBits()) >= 1){
    // Bad track, has too many shared clusters!
    return kFALSE;
  }
  return kTRUE;
}
*/
//___________________________________________________________________________________________
   
static Float_t GetFractionTPCSharedCls(const AliAODTrack *track)
{
  // Rejects tracks with shared clusters after filling a control histogram
  // This overload is used for primaries
 
  // Get the shared maps
  const TBits sharedMap = track->GetTPCSharedMap();

  return 1.*sharedMap.CountBits()/track->GetTPCNclsF();
  
}

//___________________________________________________________________________________________

Double_t AliAnalysisTaskLambdaOverK0sJets::ThetaS(TString part)
{
  // LINES OBTAINED FROM THE FEMTOSCOPY ANALYSIS:
  // private communication with Hans Beck

  // Returns the longitudinal angle of the particles propagated
  // position at R=1.25m. See
  // https://edms.cern.ch/file/406391/2/ALICE-INT-2003-038.pdf
  // for the ALICE coordinate system. Theta is zero at positive z,
  // pi/2 at z = 0 aka the xy plane and pi at negative z 

  // R^    ^  
  //  |   /
  //  |θ'/
  //  | / θ
  //  |/____>z
  // 
  // Let's compute θ' and θ = π/2 - θ'
  // where θ' can even be and should 
  // sometimes be negative
  // tan(θ') = z/R
  // θ' = arctan(z/R)
  // θ = π/2 - θ'
  //   = π/2 - arctan(z/R)
  // Note that in the doc above theta
  // is calculated as arccos(z/sqrt(x^2+y^2+z^2))

  // Array of positions is 85,105,125,..cm,
  // we take the z position at R=1.25m
  // return TMath::Pi()/2. - TMath::ATan(fXshifted[2][2]/125.);
  /*
    if( part.EqualTo("Trigger") ) 
    return TMath::Pi()/2. - TMath::ATan(fTrigSftR125[2]/125.);
    else if( part.EqualTo("Daughter") )  
    return TMath::Pi()/2. - TMath::ATan(fDaugSftR125[2]/125.);  
  */
  
  Double_t thetaS = -100.;

  if( part.EqualTo("Trigger") ) 
    thetaS = TMath::Pi()/2. - TMath::ATan(fTrigSftR125[2]/fTPCRadius);
  if( part.EqualTo("Daughter") )  
    thetaS = TMath::Pi()/2. - TMath::ATan(fDaugSftR125[2]/fTPCRadius);  

  return thetaS;

}

//___________________________________________________________________________________________

Double_t AliAnalysisTaskLambdaOverK0sJets::EtaS(TString part)
{
  // LINES OBTAINED FROM THE FEMTOSCOPY ANALYSIS:
  // private communication with Hans Beck

  // Returns the corresponding eta of a pri. part. 
  // with this particles pos at R=1.25m

  // http://en.wikipedia.org/wiki/Pseudorapidity
  // η = -ln[ tan(θ/2)]
  // printf("z: %+04.0f, thetaS %+03.2f etaS %+1.2f\n"
  // 	 ,fXshifted[2][2],ThetaS(),-TMath::Log( TMath::Tan(ThetaS()/2.) ));

  return -TMath::Log( TMath::Tan(ThetaS(part)/2.) );
}

//___________________________________________________________________________________________

Float_t AliAnalysisTaskLambdaOverK0sJets::dEtaS()
{
  // LINES OBTAINED FROM THE FEMTOSCOPY ANALYSIS:
  // private communication with Hans Beck

  // Returns the pseudorapidity star difference

  // It is important to keep the calculations easy and separated.
  // The calculation of EtaS is straight forward, one just has to
  // do it step by step to not get confused.
  return EtaS("Trigger") - EtaS("Daughter");
}

//___________________________________________________________________________________________

Float_t AliAnalysisTaskLambdaOverK0sJets::dPhiSAtR125()
{
  // LINES OBTAINED FROM THE FEMTOSCOPY ANALYSIS:
  // private communication with Hans Beck

  // returns delta phi star at R=1.2m
  // position at R=1.2m is stored as second radius
  // const Float_t distSft= TMath::Sqrt(TMath::Power(track1.fXshifted[2][0] - track2.fXshifted[2][0],2)
  // 				     +TMath::Power(track1.fXshifted[2][1] - track2.fXshifted[2][1],2));
  const Float_t distSft= TMath::Sqrt( TMath::Power(fTrigSftR125[0] - fDaugSftR125[0],2) +
				      TMath::Power(fTrigSftR125[1] - fDaugSftR125[1],2) );
  //return 2.0 * TMath::ATan(distSft/2./(125.));
  return 2.0 * TMath::ATan(distSft/2./(fTPCRadius));
}


//___________________________________________________________________________________________

void AliAnalysisTaskLambdaOverK0sJets::SetSftPosR125(const AliAODTrack *track,const Float_t bfield,const Float_t priVtx[3], TString part)
{
  // LINES OBTAINED FROM THE FEMTOSCOPY ANALYSIS:
  // private communication with Hans Beck

  // Sets the spatial position of the track at the radius R=1.25m in the shifted coordinate system
  
  // Initialize the array to something indicating there was no propagation
  if(part.EqualTo("Trigger")){  
    fTrigSftR125[0] = -9999.;
    fTrigSftR125[1] = -9999.;
    fTrigSftR125[2] = -9999.;
  }
  if(part.EqualTo("Daughter")){
    fDaugSftR125[0] = -9999.;
    fDaugSftR125[1] = -9999.;
    fDaugSftR125[2] = -9999.;
  }

  // Make a copy of the track to not change parameters of the track
  AliExternalTrackParam etp;
  etp.CopyFromVTrack(track);
  
  // The global position of the the track
  Double_t xyz[3]={-9999.,-9999.,-9999.};  

  // The radius we want to propagate to, squared
  //const Float_t RSquaredWanted(125.*125.);
  const Float_t RSquaredWanted(fTPCRadius*fTPCRadius);

  // Propagation is done in local x of the track
  for (Float_t x = 58.; x < 247.; x+=1.){
    // Starts at 83 / Sqrt(2) and goes outwards. 85/Sqrt(2) is the smallest local x
    // for global radius 85 cm. x = 245 is the outer radial limit of the TPC when
    // the track is straight, i.e. has inifinite pt and doesn't get bent. 
    // If the track's momentum is smaller than infinite, it will develop a y-component,
    // which adds to the global radius
    // We don't change the propagation steps to not mess up things!

    // Stop if the propagation was not succesful. This can happen for low pt tracks
    // that don't reach outer radii
    if(!etp.PropagateTo(x,bfield)) break;
    etp.GetXYZ(xyz); // GetXYZ returns global coordinates

    // Calculate the shifted radius we are at, squared. 
    // Compare squared radii for faster code
    Float_t shiftedRadiusSquared = (xyz[0]-priVtx[0])*(xyz[0]-priVtx[0])
      + (xyz[1]-priVtx[1])*(xyz[1]-priVtx[1]);

    // Roughly reached the radius we want
    if(shiftedRadiusSquared > RSquaredWanted){
      
      // Bigger loop has bad precision, we're nearly one centimeter too far, 
      // go back in small steps.
      while (shiftedRadiusSquared>RSquaredWanted){
	// Propagate a mm inwards
	x-=.1;
	if(!etp.PropagateTo(x,bfield)){
	  // Propagation failed but we're already with a
	  // cm precision at R=1.25m so we only break the 
	  // inner loop
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
      if(part.EqualTo("Trigger")){
	fTrigSftR125[0] = xyz[0] - priVtx[0];
	fTrigSftR125[1] = xyz[1] - priVtx[1];
	fTrigSftR125[2] = xyz[2] - priVtx[2];

	/*cout << endl
	  << xyz[0] << "   " << xyz[1] << "   " << xyz[2] << endl;
	  cout << fTrigSftR125[0] << "   " << fTrigSftR125[1] << "   " <<fTrigSftR125[2] << endl;*/
      }
      if(part.EqualTo("Daughter")){
	fDaugSftR125[0] = xyz[0] - priVtx[0];
	fDaugSftR125[1] = xyz[1] - priVtx[1];
	fDaugSftR125[2] = xyz[2] - priVtx[2];

	/*cout << endl 
	  << xyz[0] << "   " << xyz[1] << "   " << xyz[2] << endl
	  << fDaugSftR125[0] << "   " << fDaugSftR125[1] << "   " <<fDaugSftR125[2] << endl;*/
      }
 
      // Done
      return;

    } // End of if roughly reached radius
 
  } // End of coarse propagation loop

}

//___________________________________________________________________________________________

void AliAnalysisTaskLambdaOverK0sJets::RecCascade(const AliAODTrack *trk1,const AliAODTrack *trk2,const AliAODTrack *trkBch,TString histo)
{
  // Local method to reconstruct cascades candidates from the combinations of three tracks
  // The input tracks correspond to the trigger particle and the daughter tracks of the V0 candidate (correlation step)
  // The trigger particle track will be always consider as a possible daughter of the V0 which coming from the Cascade decay.
  // The daughters of the V0 candidates are switched to be the bachelor track for the Cascade reconstruction.

  Float_t lMassBach=0., lPtot2Bach=0., lEBach=0.;
  Float_t lMassLambda=0., lPtot2Lambda=0., lELambda = 0.; 
  Float_t pLambda[3] = {0.,0.,0.};
  Float_t pCascade[3] = {0.,0.,0.};
  Float_t lMassCascade = 0., lPtot2Cascade=0.;

  // Two loops are done to consider the posibility to reconstruct a Xi or an Omega
  for(Int_t i=0;i<2;i++){

    // 0. Check the charge for both tracks: trk1 & trk2. 
    //    Usefull in the Lambda step.
    if( trk1->Charge() == trk2->Charge() ) 
      continue;
   
    // 1. Bachelor: Allocation for the track
    if(i==0) // Xi 
      lMassBach = TDatabasePDG::Instance()->GetParticle(kPiMinus)->Mass();
    else if(i==1) //Omega
      lMassBach = TDatabasePDG::Instance()->GetParticle(kKMinus)->Mass();

    lPtot2Bach = TMath::Power(trkBch->P(),2);

    lEBach = TMath::Sqrt(lPtot2Bach + lMassBach*lMassBach);

    // 2. Lambda: Kinematical properties
    lMassLambda = TDatabasePDG::Instance()->GetParticle(kLambda0)->Mass();
      
    pLambda[0] = trk1->Px() + trk2->Px();
    pLambda[1] = trk1->Py() + trk2->Py();
    pLambda[2] = trk1->Pz() + trk2->Pz();

    lPtot2Lambda = pLambda[0]*pLambda[0] +  pLambda[1]*pLambda[1] +  pLambda[2]*pLambda[2];

    lELambda = TMath::Sqrt(lPtot2Lambda + lMassLambda*lMassLambda);

    // 3. Cascade: Reconstruction
    pCascade[0] = pLambda[0] + trkBch->Px();
    pCascade[1] = pLambda[1] + trkBch->Py();
    pCascade[2] = pLambda[2] + trkBch->Pz();

    lPtot2Cascade = pCascade[0]*pCascade[0] + pCascade[1]*pCascade[1] + pCascade[2]*pCascade[2];

    lMassCascade = TMath::Sqrt( TMath::Power(lEBach+lELambda,2) - lPtot2Cascade );
   
    // 4. Filling histograms
    if( histo.Contains("K0s") ) {
      if(i==0) // Xi 
	fV0MassCascade->Fill(lMassCascade,1);
      else if(i==1) //Omega
	fV0MassCascade->Fill(lMassCascade,3);
    }
    else if( histo.Contains("AntiLambda") ) {
      if(i==0) // Xi 
	fV0MassCascade->Fill(lMassCascade,9);
      else if(i==1) //Omega
	fV0MassCascade->Fill(lMassCascade,11);
    }
    else if( histo.Contains("Lambda") ) {
      if(i==0) // Xi 
	fV0MassCascade->Fill(lMassCascade,5);
      else if(i==1) //Omega
	fV0MassCascade->Fill(lMassCascade,7);
    }

  }
  
}

//___________________________________________________________________________________________
 
void AliAnalysisTaskLambdaOverK0sJets::V0Loop(V0LoopStep_t step, Bool_t isTriggered, Int_t iArray, Int_t idTrig) 
{ 
  // Three options for the 'step' variable:
  // 1) TriggerCheck
  // 2) Reconstruction

  const AliAODTrack *trkTrig = 0x0;
  Float_t  ptTrig  = -100.;
  Float_t  phiTrig = -100.;
  Float_t  etaTrig = -100.; 
  Double_t pTrig[3]; 

  if( (step==kTriggerCheck || isTriggered) && idTrig>=0 ){
    trkTrig = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(idTrig));
    if(!trkTrig) AliFatal("Not a standard AOD"); 
    ptTrig  = trkTrig->Pt();
    phiTrig = trkTrig->Phi();
    etaTrig = trkTrig->Eta();
    trkTrig->GetPxPyPz(pTrig); 
  }
  
  // *************************************************
  // Centrality selection
  AliCentrality *cent = fAOD->GetCentrality();
  Float_t centrality = cent->GetCentralityPercentile("V0M");
  Int_t curCentBin = CentBin(centrality);

  // *************************************************
  // MC Event
  TClonesArray *stackMC = 0x0;
  Float_t mcXv=0., mcYv=0., mcZv=0.;
   
  if(fIsMC){
    TList *lst = fAOD->GetList();
    stackMC = (TClonesArray*)lst->FindObject(AliAODMCParticle::StdBranchName());
    if (!stackMC) {
      Printf("ERROR: stack not available");
    }

    AliAODMCHeader *mcHdr = 
      (AliAODMCHeader*)lst->FindObject(AliAODMCHeader::StdBranchName());
    
    mcXv=mcHdr->GetVtxX(); mcYv=mcHdr->GetVtxY(); mcZv=mcHdr->GetVtxZ();
  }
  
  // *************************************************
  // V0 loop - AOD
  const AliAODVertex *vtx=fAOD->GetPrimaryVertex();
  Float_t xv=vtx->GetX(), yv=vtx->GetY(), zv=vtx->GetZ();
  Int_t nV0sTot = fAOD->GetNumberOfV0s();

  for (Int_t iV0 = 0; iV0 < nV0sTot; iV0++) {
    
    AliAODv0 *v0=fAOD->GetV0(iV0);
    if (!v0) continue;
    if (!AcceptV0(fAOD->GetPrimaryVertex(),v0)) continue;
    
    const AliAODTrack *ntrack=(AliAODTrack *)v0->GetDaughter(1);
    const AliAODTrack *ptrack=(AliAODTrack *)v0->GetDaughter(0);

    // Decay vertex
    Double_t xyz[3]; v0->GetSecondaryVtx(xyz);
    Float_t dx=xyz[0]-xv, dy=xyz[1]-yv;//, dz=xyz[2]-zv;
   
    // Momentum: 2D & 3D
    Float_t pt=TMath::Sqrt(v0->Pt2V0());
    //Float_t p=v0->P();

    // Decay length: 2D & 3D 
    Float_t lt=TMath::Sqrt(dx*dx + dy*dy); 
    //Float_t dl=TMath::Sqrt(dx*dx + dy*dy + dz*dz);  
    
    Float_t dlK = 0.4977*lt/pt;
    Float_t dlL = 1.1157*lt/pt; 
    
    //Float_t dlK  = v0->MassK0Short()*dl/p;
    //Float_t dlL  = v0->MassLambda()*dl/p;
    //Float_t dlAL = v0->MassAntiLambda()*dl/p;

    // ctau
    Bool_t ctK=kTRUE;  if (dlK > fMaxCtau*2.68 || dlK < fMinCtau*2.68) ctK=kFALSE; 
    Bool_t ctL=kTRUE;  if (dlL > fMaxCtau*7.89 || dlL < fMinCtau*7.89) ctL=kFALSE; 
    Bool_t ctAL=kTRUE; if (dlL > fMaxCtau*7.89 || dlL < fMinCtau*7.89) ctAL=kFALSE;    

    //  ---- Daughter tracks properties:
    // Pt
    Float_t lPtNeg = ntrack->Pt();
    Float_t lPtPos = ptrack->Pt();  
    // Momentum
    Double_t pNegDaug[3];  ntrack->GetPxPyPz(pNegDaug);                  
    Double_t pPosDaug[3];  ptrack->GetPxPyPz(pPosDaug);
    // Phi
    Float_t phiNeg = ntrack->Phi();
    Float_t phiPos = ptrack->Phi();
    // Eta
    Float_t etaNeg = ntrack->Eta();
    Float_t etaPos = ptrack->Eta();
    //  Number of TPC Clusters 
    Float_t nClsTPCPos = ptrack->GetTPCClusterInfo(2,1);
    Float_t nClsTPCNeg = ntrack->GetTPCClusterInfo(2,1); 
    // Number of clusters of ITS
    Double_t posITSNcls = ptrack->GetITSNcls();   
    Double_t negITSNcls = ntrack->GetITSNcls();

    //  ---- V0 candidate properties:
    // Armenteros variables:
    Float_t lAlphaV0      =  v0->AlphaV0();
    Float_t lPtArmV0      =  v0->PtArmV0();
    // dca to primary vertex
    Float_t dcaNeg = v0->DcaNegToPrimVertex();
    Float_t dcaPos = v0->DcaPosToPrimVertex();
    // dca between daughters
    Float_t dca   = v0->DcaV0Daughters();
    // cpa
    Float_t cpa   = v0->CosPointingAngle(fAOD->GetPrimaryVertex());
    // eta
    Float_t lEta  = v0->PseudoRapV0();
    // phi
    Float_t lPhi  = v0->Phi();
    //lPhi  = ( (lPhi < 0) ? lPhi + 2*TMath::Pi() : lPhi );    

    // Fraction of TPC Shared Cluster 
    Float_t fracPosDaugTPCSharedMap = GetFractionTPCSharedCls(ptrack);
    Float_t fracNegDaugTPCSharedMap = GetFractionTPCSharedCls(ntrack);

    // rapidity
    Float_t rapK0s = v0->RapK0Short();
    Float_t rapLambda = v0->RapLambda();

    if(fUseEtaCut){
           rapK0s = lEta;
           rapLambda = lEta;
    }
   
    // **********************************
    // PID - tracks  
    Float_t pPos = -100.;
    Float_t pNeg = -100.;
    
    //Float_t dedxPos = -1000.;
    //Float_t dedxNeg = -1000.;
    //Float_t nsigPosPion   = 0.;
    //Float_t nsigNegPion   = 0.;
    Float_t nsigPosProton = 0.;
    Float_t nsigNegProton = 0.;

      /*
    if(fUsePID && !fIsMC) {     
      const AliAODPid *pidNeg = ntrack->GetDetPid();
      const AliAODPid *pidPos = ptrack->GetDetPid();
      
      if (pidNeg && pidPos) {
	pPos = pidPos->GetTPCmomentum();
	pNeg = pidNeg->GetTPCmomentum();
	//dedxPos = pidPos->GetTPCsignal()/47.; 
	//dedxNeg = pidNeg->GetTPCsignal()/47.; 
	
	
	if(pPos<1.){
	  //nsigPosPion   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(ptrack,AliPID::kPion));
	  nsigPosProton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(ptrack,AliPID::kProton));
	}
	if(pNeg<1.){
	  //nsigNegPion   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(ntrack,AliPID::kPion));
	  nsigNegProton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(ntrack,AliPID::kProton));
	}
	
      }
      
    }*/

    // **********************************
    // Disentangle the V0 candidate
    Float_t massK0s = 0., mK0s = 0., sK0s = 0.;
    Float_t massLambda = 0., mLambda = 0., sL = 0.;
    Float_t massAntiLambda = 0., sAL = 0.;

    Bool_t isCandidate2K0s = kFALSE;
    massK0s = v0->MassK0Short();
    mK0s = TDatabasePDG::Instance()->GetParticle(kK0Short)->Mass();
    if( fCollision.Contains("PbPb2010") )
      sK0s = kCteK0s2010[curCentBin] + kLinearK0s2010[curCentBin]*pt;
    else if( fCollision.Contains("PbPb2011") ) 
      sK0s = kCteK0s2011[curCentBin] + kLinearK0s2011[curCentBin]*pt;
    if ( TMath::Abs(mK0s-massK0s) < 3*sK0s )  isCandidate2K0s = kTRUE;     
    
    Bool_t isCandidate2Lambda = kFALSE;
    massLambda = v0->MassLambda();
    mLambda = TDatabasePDG::Instance()->GetParticle(kLambda0)->Mass();
    if( fCollision.Contains("PbPb2010") )
      sL = kCteLambda2010[curCentBin] + kLinearLambda2010[curCentBin]*pt;
    else if( fCollision.Contains("PbPb2011") ) 
      sL = kCteLambda2011[curCentBin] + kLinearLambda2011[curCentBin]*pt;
    if (TMath::Abs(mLambda-massLambda) < 3*sL)  isCandidate2Lambda = kTRUE;  
    
    Bool_t isCandidate2LambdaBar = kFALSE;
    massAntiLambda = v0->MassAntiLambda();
    if( fCollision.Contains("PbPb2010") )
      sAL = kCteAntiLambda2010[curCentBin] + kLinearAntiLambda2010[curCentBin]*pt;
    else if( fCollision.Contains("PbPb2011") ) 
      sAL = kCteAntiLambda2011[curCentBin] + kLinearAntiLambda2011[curCentBin]*pt;
    if (TMath::Abs(mLambda-massAntiLambda) < 3*sAL)  isCandidate2LambdaBar = kTRUE; 

    // **********************************
    // MC Association:
    Bool_t lComeFromSigma     = kTRUE; 
    Bool_t lCheckMcK0Short    = kTRUE;
    Bool_t lCheckMcLambda     = kTRUE;
    Bool_t lCheckMcAntiLambda = kTRUE;
    Bool_t lComeFromXi        = kTRUE; 
	
    Int_t lMCAssocNegDaug = -100;
    Int_t lMCAssocPosDaug = -100;  
    
    // ********* MC - Association *********
    // In case of injected-MC, the correlations might be done with only natural particles 
    Bool_t isNaturalPart = kTRUE;
    if(step==kReconstruction){
      
      if(fIsMC){        
      	if(!stackMC) goto noas;

	isNaturalPart = kFALSE;

	lComeFromSigma     = kFALSE; 
	lCheckMcK0Short    = kFALSE;
	lCheckMcLambda     = kFALSE;
	lCheckMcAntiLambda = kFALSE;
	lComeFromXi        = kFALSE;

	Int_t ntrkMC=stackMC->GetEntriesFast();
	
	Int_t nlab = TMath::Abs(ntrack->GetLabel()); // ** UInt_t
	Int_t plab = TMath::Abs(ptrack->GetLabel());
  
	// To avoid futher problems 
	if ( (nlab<0 || plab<0) ||
	     (nlab>=ntrkMC || plab>=ntrkMC) )
	  goto noas;      

	AliAODMCParticle *nPart=(AliAODMCParticle*)stackMC->UncheckedAt(nlab);
	AliAODMCParticle *pPart=(AliAODMCParticle*)stackMC->UncheckedAt(plab);

	if(!nPart || !pPart)   goto noas;

	// MC origin of daughters: Primaries?
	if( nPart->IsPhysicalPrimary() ) lMCAssocNegDaug = 1;
	if( pPart->IsPhysicalPrimary() ) lMCAssocPosDaug = 1;
	
	/*
	  if ( TMath::Abs(nPart->Eta()) > fMaxEtaDaughter ||
	  TMath::Abs(pPart->Eta()) > fMaxEtaDaughter )
	  goto noas;
	*/
	/*
	// Daughter momentum cut
	if ( ( nPart->Pt() < fMinPtDaughter ) || 
	( pPart->Pt() < fMinPtDaughter )  ) 
	goto noas;
	*/

	// ----------------------------------------
	
	Int_t lPDGCodeNegDaughter = nPart->GetPdgCode();
	Int_t lPDGCodePosDaughter = pPart->GetPdgCode();
	
	Int_t ipMother = pPart->GetMother();
	Int_t inMother = nPart->GetMother();
	
	if(inMother<0 || inMother>=ntrkMC) lMCAssocNegDaug = 6;
	if(ipMother<0 || ipMother>=ntrkMC) lMCAssocPosDaug = 6;

	if(inMother<0 || inMother>=ntrkMC) {  goto noas;}
	if(inMother != ipMother) { // did the negative daughter decay ?
	  AliAODMCParticle *negMotherOfMotherPart = (AliAODMCParticle*)stackMC->UncheckedAt(inMother);
	  if (negMotherOfMotherPart->GetMother() != ipMother) 
	    goto noas;
	}
	
	if (ipMother<0 || ipMother>=ntrkMC)
	  goto noas;     
	
	AliAODMCParticle *p0=(AliAODMCParticle*)stackMC->UncheckedAt(ipMother);
	if(!p0) 
	  goto noas; 

	// ----------------------------------------
	
	if ( (ipMother>=fEndOfHijingEvent) && 
	     (fEndOfHijingEvent!=-1)     && 
	     (p0->GetMother()<0) ) 
	  isNaturalPart = kFALSE; 
	else  isNaturalPart = kTRUE; 

	// ----------------------------------------
	
	if(fSeparateInjPart && !isNaturalPart) goto noas;     
	
	Int_t lPDGCodeV0 = p0->GetPdgCode();
	
	// MC origin of daughters:
	//Decay from Weak Decay?
	if( (TMath::Abs(lPDGCodeV0) == kK0Short) || (TMath::Abs(lPDGCodeV0) == kLambda0) || 
	    (TMath::Abs(lPDGCodeV0) == kSigmaMinus) || (TMath::Abs(lPDGCodeV0) == kSigmaPlus) ||
	    (TMath::Abs(lPDGCodeV0) == kSigma0) )
	  { lMCAssocNegDaug = 2; 	  lMCAssocPosDaug = 2; }
	// Cascade Gamma conversion
	if( (TMath::Abs(lPDGCodeV0) == kXiMinus) ||
	    (TMath::Abs(lPDGCodeV0) == kOmegaMinus) )
	  { lMCAssocNegDaug = 3; 	  lMCAssocPosDaug = 3; }
	// Gamma conversion
	else if( TMath::Abs(lPDGCodeV0) == kGamma )
	  { lMCAssocNegDaug = 4; 	  lMCAssocPosDaug = 4; }
	// Unidentied mother:
	else 
	  { lMCAssocNegDaug = 5; 	  lMCAssocPosDaug = 5; }


	Int_t lIndexMotherOfMother   = p0->GetMother();
	Int_t lPdgcodeMotherOfMother = 0;
	if (lIndexMotherOfMother != -1) {
	  AliAODMCParticle *lMCAODMotherOfMother=(AliAODMCParticle*)stackMC->UncheckedAt(lIndexMotherOfMother);
	  if (lMCAODMotherOfMother) {lPdgcodeMotherOfMother = lMCAODMotherOfMother->GetPdgCode();}
	}
	
	/*
3	// Daughter momentum cut: ! FIX it in case of AOD ! //MC or REc
	if ( (nPart->Pt()  < fMinPtDaughter ) ||
	(pPart->Pt()  < fMinPtDaughter ) )
	goto noas;
	*/

	if( (lPDGCodeV0 != kK0Short) &&
	    (lPDGCodeV0 != kLambda0) &&
	    (lPDGCodeV0 != kLambda0Bar) ) 
	  goto noas;
	
	     
	// ----------------------------------------
      
	// K0s
	if( (lPDGCodePosDaughter==+211) && (lPDGCodeNegDaughter==-211) &&
	    (inMother==ipMother) && (lPDGCodeV0==310) ) {
	  
	  if ( ((AliAODMCParticle*)stackMC->UncheckedAt(ipMother))->IsPrimary()  )
	    lCheckMcK0Short  = kTRUE;
	
	}
	// Lambda
	else if( (lPDGCodePosDaughter==+2212) && (lPDGCodeNegDaughter==-211)  &&
		 (inMother==ipMother) && (lPDGCodeV0==3122)  ){
	  
	  if ( ( TMath::Abs(lPdgcodeMotherOfMother) == 3212) /*||
							       ( TMath::Abs(lPdgcodeMotherOfMother) == 3224) ||
							       ( TMath::Abs(lPdgcodeMotherOfMother) == 3214) ||
							       ( TMath::Abs(lPdgcodeMotherOfMother) == 3114)*/
	       ) lComeFromSigma = kTRUE;
	  else lComeFromSigma = kFALSE; 
	  
	  if ( ((AliAODMCParticle*)stackMC->UncheckedAt(ipMother))->IsPrimary() || 
	       ( !(((AliAODMCParticle*)stackMC->UncheckedAt(ipMother))->IsPrimary() ) 
		 && (lComeFromSigma==kTRUE) )
	       ) lCheckMcLambda  = kTRUE; 
	  
	  
	  if ( TMath::Abs(lPdgcodeMotherOfMother) == 3312) 
	    lComeFromXi = kTRUE;
	  
	}
	// AntiLambda
	else if( (lPDGCodePosDaughter==211) && (lPDGCodeNegDaughter==-2212) &&
		 (inMother==ipMother) && (lPDGCodeV0==-3122) ) {
	  
	  
	  if ( ( TMath::Abs(lPdgcodeMotherOfMother) == 3212) /*||
							       ( TMath::Abs(lPdgcodeMotherOfMother) == 3224) ||
							       ( TMath::Abs(lPdgcodeMotherOfMother) == 3214) ||
							       ( TMath::Abs(lPdgcodeMotherOfMother) == 3114)*/
	       ) lComeFromSigma = kTRUE;
	  else lComeFromSigma = kFALSE;  
	  
	  if ( ((AliAODMCParticle*)stackMC->UncheckedAt(ipMother))->IsPrimary() || 
	       ( (!((AliAODMCParticle*)stackMC->UncheckedAt(ipMother))->IsPrimary()) 
		 && (lComeFromSigma==kTRUE) )
	       ) lCheckMcAntiLambda  = kTRUE;
	  
	  if ( TMath::Abs(lPdgcodeMotherOfMother) == 3312 || TMath::Abs(lPdgcodeMotherOfMother) == 3322 ) 
	    lComeFromXi = kTRUE;
	  
	}
	
	//  ----------------------------------------
	
	if ((p0->Pt())<pMin) goto noas;
	if (TMath::Abs(p0->Y())>fYMax ) goto noas;
	
	Float_t dxAs = mcXv - p0->Xv(),  dyAs = mcYv - p0->Yv(),  dzAs = mcZv - p0->Zv();
	Float_t l = TMath::Sqrt(dxAs*dxAs + dyAs*dyAs + dzAs*dzAs);
	
	dxAs = mcXv - pPart->Xv(); dyAs = mcYv - pPart->Yv();
	//Float_t ltAs = TMath::Sqrt(dxAs*dxAs + dyAs*dyAs);
	Float_t ptAs = p0->Pt();
	Float_t rapAs = p0->Y();
	Float_t etaAs = p0->Eta();

	if(fUseEtaCut){
	      rapAs = etaAs;
	}

	// phi resolution for V0-reconstruction and daughter tracks
	Float_t resEta = p0->Eta() - v0->Eta();	
	Float_t resPhi = p0->Phi() - v0->Phi();	
	Float_t resPt  = p0->Pt() - v0->Pt();	
	
	Float_t resEtaPosDaug = pPart->Eta() - ptrack->Eta();	
	Float_t resPhiPosDaug = pPart->Phi() - ptrack->Phi();	
	Float_t resPtPosDaug  = pPart->Pt() - ptrack->Pt();

	Float_t resEtaNegDaug = nPart->Eta() - ntrack->Eta();	
	Float_t resPhiNegDaug = nPart->Phi() - ntrack->Phi();	
	Float_t resPtNegDaug  = nPart->Pt() - ntrack->Pt();

	if ( (l < 0.01)  &&  (ptAs<10.) ) { // Primary V0
	  
	  // K0s:
	  if(ctK && lCheckMcK0Short){ 
	    
	    // Natural particles
	    if(isNaturalPart){

	      if( (dcaPos>fDCAToPrimVtx) && (dcaNeg>fDCAToPrimVtx) && (nClsTPCPos>fDaugNClsTPC) && (nClsTPCNeg>fDaugNClsTPC) && (dca<fMaxDCADaughter) ){

		fK0sAssocPt->Fill(ptAs);
		fK0sAssocPtRap->Fill(ptAs,rapAs,centrality);
		fK0sAssocPtPhiEta[curCentBin]->Fill(p0->Phi(),etaAs,ptAs);
	      
		// Armenteros Pod.  and rapidity cut
		if( (lPtArmV0 > TMath::Abs(0.2*lAlphaV0) ) && TMath::Abs(rapAs)<fYMax ){ 
		 		
		  Double_t effK0sArm[3] = {v0->MassK0Short(),ptAs,rapAs};
		  /*Double_t effK0sVtx[4] = {v0->MassK0Short(),ptAs,rapAs,zv};
		  Double_t effK0sDCA[4] = {v0->MassK0Short(),ptAs,rapAs,dca};
		  Double_t effK0sCPA[4] = {v0->MassK0Short(),ptAs,rapAs,cpa};
		  Double_t effK0sShTPCcls[5] = {v0->MassK0Short(),ptAs,rapAs,fracPosDaugTPCSharedMap,fracNegDaugTPCSharedMap};
		  Double_t effK0sDaugPt[5] = {v0->MassK0Short(),ptAs,rapAs,lPtPos,lPtNeg};
		  Double_t effK0sCtau[4]   = {v0->MassK0Short(),ptAs,rapAs,dlK};
		  Double_t effK0sFidVol[4] = {v0->MassK0Short(),ptAs,rapAs,lt};*/

		  // Distributions for the efficiency (systematics chechks)
		  fK0sAssocPtMassArm[curCentBin]->Fill(effK0sArm);
		  /*fK0sAssocMassPtVtx[curCentBin]->Fill(effK0sVtx);
		  fK0sAssocMassPtDCADaug[curCentBin]->Fill(effK0sDCA);
		  fK0sAssocMassPtCPA[curCentBin]->Fill(effK0sCPA);
		  fK0sAssocMassPtShTPCcls[curCentBin]->Fill(effK0sShTPCcls);
		  fK0sAssocMassPtDaugPt[curCentBin]->Fill(effK0sDaugPt);
		  fK0sAssocMassPtCtau[curCentBin]->Fill(effK0sCtau);
		  fK0sAssocMassPtFidVolume[curCentBin]->Fill(effK0sFidVol);*/
		
		  fK0sMCResEta->Fill(resEta,pt,centrality);
		  fK0sMCResPhi->Fill(resPhi,pt,centrality);
		  fK0sMCResPt->Fill(resPt,pt,centrality);
		
		  fK0sPosMCResEta->Fill(resEtaPosDaug,pt,centrality);
		  fK0sPosMCResPhi->Fill(resPhiPosDaug,pt,centrality);
		  fK0sPosMCResPt->Fill(resPtPosDaug,pt,centrality);

		  fK0sNegMCResEta->Fill(resEtaNegDaug,pt,centrality);
		  fK0sNegMCResPhi->Fill(resPhiNegDaug,pt,centrality);
		  fK0sNegMCResPt->Fill(resPtNegDaug,pt,centrality);

		}

	      } // End selection in the dca to prim. vtx and the number of clusters

	      // Distributions for the efficiency (Systematic checks)
	      /*
	      if( (lPtArmV0 > TMath::Abs(0.2*lAlphaV0) ) && TMath::Abs(rapAs)<fYMax ){ 
	      
		//  Cut in the DCA ToPrim Vtx
		if( (nClsTPCPos>fDaugNClsTPC) && (nClsTPCNeg>fDaugNClsTPC) ){

		  Double_t effK0sdcaPV[5] = {v0->MassK0Short(),ptAs,rapAs,dcaPos,dcaNeg};
		  fK0sAssocMassPtDCAPV[curCentBin]->Fill(effK0sdcaPV);
		}		  

		// cut in the number of tpc clusters
		if( (dcaPos>fDCAToPrimVtx) && (dcaNeg>fDCAToPrimVtx) && TMath::Abs(rapAs)<fYMax ){
		  
		  Double_t effK0sTPCcls[5]  = {v0->MassK0Short(),ptAs,rapAs,nClsTPCPos,nClsTPCNeg};
		  fK0sAssocMassPtDaugNClsTPC[curCentBin]->Fill(effK0sTPCcls);
		  
		}

	      } // End selection for systematics
	      */
	      
	    } // End natural particle selection
	    // Embeded particles
	    if(!isNaturalPart){ 

	      if( (dcaPos>fDCAToPrimVtx) && (dcaNeg>fDCAToPrimVtx) && (nClsTPCPos>fDaugNClsTPC) && (nClsTPCNeg>fDaugNClsTPC) && (dca<fMaxDCADaughter) ){

		fK0sAssocPtRapEmbeded->Fill(ptAs,rapAs,centrality);

		if( (lPtArmV0 > TMath::Abs(0.2*lAlphaV0)) && TMath::Abs(rapAs)<fYMax ){
		  
		  Double_t effK0sArm[3] = {v0->MassK0Short(),ptAs,rapAs};
		  /*Double_t effK0sVtx[4] = {v0->MassK0Short(),ptAs,rapAs,zv};
		  Double_t effK0sDCA[4] = {v0->MassK0Short(),ptAs,rapAs,dca};
		  Double_t effK0sCPA[4] = {v0->MassK0Short(),ptAs,rapAs,cpa};
		  Double_t effK0sShTPCcls[5] = {v0->MassK0Short(),ptAs,rapAs,fracPosDaugTPCSharedMap,fracNegDaugTPCSharedMap};
		  Double_t effK0sDaugPt[5] = {v0->MassK0Short(),ptAs,rapAs,lPtPos,lPtPos};
		  Double_t effK0sCtau[4]   = {v0->MassK0Short(),ptAs,rapAs,dlK};
		  Double_t effK0sFidVol[4] = {v0->MassK0Short(),ptAs,rapAs,lt};*/

		  // Distributions for the efficiency (systematics chechks)
		  fK0sAssocPtMassArmEmbeded[curCentBin]->Fill(effK0sArm);	
		  /*fK0sAssocMassPtVtxEmbeded[curCentBin]->Fill(effK0sVtx);
		  fK0sAssocMassPtDCADaugEmbeded[curCentBin]->Fill(effK0sDCA);
		  fK0sAssocMassPtCPAEmbeded[curCentBin]->Fill(effK0sCPA);
		  fK0sAssocMassPtShTPCclsEmbeded[curCentBin]->Fill(effK0sShTPCcls);
		  fK0sAssocMassPtDaugPtEmbeded[curCentBin]->Fill(effK0sDaugPt);
		  fK0sAssocMassPtCtauEmbeded[curCentBin]->Fill(effK0sCtau);
		  fK0sAssocMassPtFidVolumeEmbeded[curCentBin]->Fill(effK0sFidVol);*/
		}

	      } // End selection in the dca to prim. vtx and the number of clusters

	      // Distributions for the efficiency (Systematic checks)
	      /*
	      if( (lPtArmV0 > TMath::Abs(0.2*lAlphaV0) ) && TMath::Abs(rapAs)<fYMax ){ 

		//  Cut in the DCA ToPrim Vtx
		if( (nClsTPCPos>fDaugNClsTPC) && (nClsTPCNeg>fDaugNClsTPC) ){

		  Double_t effK0sdcaPV[5] = {v0->MassK0Short(),ptAs,rapAs,dcaPos,dcaNeg};
    		  fK0sAssocMassPtDCAPVEmbeded[curCentBin]->Fill(effK0sdcaPV);
		}		  

		// cut in the number of tpc clusters
		if( (dcaPos>fDCAToPrimVtx) && (dcaNeg>fDCAToPrimVtx) && TMath::Abs(rapAs)<fYMax ){

		  Double_t effK0sTPCcls[5]  = {v0->MassK0Short(),ptAs,rapAs,nClsTPCPos,nClsTPCNeg};
		  fK0sAssocMassPtDaugNClsTPCEmbeded[curCentBin]->Fill(effK0sTPCcls);
		}

	      } // End selection for systematics
	      */
	      
	    } // End embeded particle selection

	  }  // End K0s selection
	  
	  // Lambda:
	  if(ctL && lCheckMcLambda && (TMath::Abs(nsigPosProton)<fNSigma) ) {  
	    
	    // Natural particles
	    if(isNaturalPart){

	      if( (dcaPos>fDCAToPrimVtx) && (dcaNeg>fDCAToPrimVtx) && (nClsTPCPos>fDaugNClsTPC) && (nClsTPCNeg>fDaugNClsTPC) && (dca<fMaxDCADaughter) ){

		fLambdaAssocPt->Fill(ptAs);
		fLambdaAssocPtRap->Fill(ptAs,rapAs,centrality);
		fLambdaAssocPtPhiEta[curCentBin]->Fill(p0->Phi(),etaAs,ptAs);

		// Rapidity cut
		if(TMath::Abs(rapAs)<fYMax)  {

		  Double_t effLambda[3] = {v0->MassLambda(),ptAs,rapAs};/*
		  Double_t effLambdaVtx[4] = {v0->MassLambda(),ptAs,rapAs,zv};
		  Double_t effLambdaDCA[4] = {v0->MassLambda(),ptAs,rapAs,dca};
		  Double_t effLambdaCPA[4] = {v0->MassLambda(),ptAs,rapAs,cpa};
		  Double_t effLambdaShTPCcls[5] = {v0->MassLambda(),ptAs,rapAs,fracPosDaugTPCSharedMap,fracNegDaugTPCSharedMap};
		  Double_t effLambdaDaugPt[5] = {v0->MassLambda(),ptAs,rapAs,lPtPos,lPtNeg};
		  Double_t effLambdaCtau[4]   = {v0->MassLambda(),ptAs,rapAs,dlL};
		  Double_t effLambdaFidVol[4] = {v0->MassLambda(),ptAs,rapAs,lt};*/

		  // Distributions for the efficiency (systematics chechks)
		  fLambdaAssocMassPtRap[curCentBin]->Fill(effLambda);
		  /*fLambdaAssocMassPtVtx[curCentBin]->Fill(effLambdaVtx);
		  fLambdaAssocMassPtDCADaug[curCentBin]->Fill(effLambdaDCA);
		  fLambdaAssocMassPtCPA[curCentBin]->Fill(effLambdaCPA);
		  fLambdaAssocMassPtShTPCcls[curCentBin]->Fill(effLambdaShTPCcls);
		  fLambdaAssocMassPtDaugPt[curCentBin]->Fill(effLambdaDaugPt);
		  fLambdaAssocMassPtCtau[curCentBin]->Fill(effLambdaCtau);
		  fLambdaAssocMassPtFidVolume[curCentBin]->Fill(effLambdaFidVol);*/

		  if( !isCandidate2K0s && !isCandidate2LambdaBar)
		    fLambdaAssocMassPtRap2[curCentBin]->Fill(effLambda);

		  fLambdaMCResEta->Fill(resEta,pt,centrality);
		  fLambdaMCResPhi->Fill(resPhi,pt,centrality);
		  fLambdaMCResPt->Fill(resPt,pt,centrality);
		
		  fLambdaPosMCResEta->Fill(resEtaPosDaug,pt,centrality);
		  fLambdaPosMCResPhi->Fill(resPhiPosDaug,pt,centrality);
		  fLambdaPosMCResPt->Fill(resPtPosDaug,pt,centrality);

		  fLambdaNegMCResEta->Fill(resEtaNegDaug,pt,centrality);
		  fLambdaNegMCResPhi->Fill(resPhiNegDaug,pt,centrality);
		  fLambdaNegMCResPt->Fill(resPtNegDaug,pt,centrality);

		}

	      } // End selection in the dca to prim. vtx and the number of clusters
	      
	      // Distributions for the efficiency (Systematic checks)
	      /*
	      if( TMath::Abs(rapAs)<fYMax ){ 
		
		//  Cut in the DCA ToPrim Vtx
		if( (nClsTPCPos>fDaugNClsTPC) && (nClsTPCNeg>fDaugNClsTPC) ){

		  Double_t effLambdadcaPV[5] = {v0->MassLambda(),ptAs,rapAs,dcaPos,dcaNeg};
		  fLambdaAssocMassPtDCAPV[curCentBin]->Fill(effLambdadcaPV);
		}		  

		// cut in the number of tpc clusters
		if( (dcaPos>fDCAToPrimVtx) && (dcaNeg>fDCAToPrimVtx) && TMath::Abs(rapAs)<fYMax){
		 
		  Double_t effLambdaTPCcls[5]  = {v0->MassLambda(),ptAs,rapAs,nClsTPCPos,nClsTPCNeg};
		  fLambdaAssocMassPtDaugNClsTPC[curCentBin]->Fill(effLambdaTPCcls);
		}

	      } // End selection for systematics
	      */
	    } // End natural particle selection
	    // Embeded particles
	    if(!isNaturalPart){

	      if( (dcaPos>fDCAToPrimVtx) && (dcaNeg>fDCAToPrimVtx) && (nClsTPCPos>fDaugNClsTPC) && (nClsTPCNeg>fDaugNClsTPC) && (dca<fMaxDCADaughter) ){
	      
		if( TMath::Abs(rapAs)<fYMax ){

		  Double_t effLambda[3] = {v0->MassLambda(),ptAs,rapAs};
		  /*Double_t effLambdaVtx[4] = {v0->MassLambda(),ptAs,rapAs,zv};
		  Double_t effLambdaDCA[4] = {v0->MassLambda(),ptAs,rapAs,dca};
		  Double_t effLambdaCPA[4] = {v0->MassLambda(),ptAs,rapAs,cpa};
		  Double_t effLambdaShTPCcls[5] = {v0->MassLambda(),ptAs,rapAs,fracPosDaugTPCSharedMap,fracNegDaugTPCSharedMap};
		  Double_t effLambdaDaugPt[5] = {v0->MassLambda(),ptAs,rapAs,lPtPos,lPtNeg};
		  Double_t effLambdaCtau[4]   = {v0->MassLambda(),ptAs,rapAs,dlL};
		  Double_t effLambdaFidVol[4] = {v0->MassLambda(),ptAs,rapAs,lt};*/

		  // Distributions for the efficiency (systematics chechks)
		  fLambdaAssocMassPtRapEmbeded[curCentBin]->Fill(effLambda);
		  /*fLambdaAssocMassPtVtxEmbeded[curCentBin]->Fill(effLambdaVtx);
		  fLambdaAssocMassPtDCADaugEmbeded[curCentBin]->Fill(effLambdaDCA);
		  fLambdaAssocMassPtCPAEmbeded[curCentBin]->Fill(effLambdaCPA);
		  fLambdaAssocMassPtShTPCclsEmbeded[curCentBin]->Fill(effLambdaShTPCcls);
		  fLambdaAssocMassPtDaugPtEmbeded[curCentBin]->Fill(effLambdaDaugPt);
		  fLambdaAssocMassPtCtauEmbeded[curCentBin]->Fill(effLambdaCtau);
		  fLambdaAssocMassPtFidVolumeEmbeded[curCentBin]->Fill(effLambdaFidVol);*/

		  if( !isCandidate2K0s && !isCandidate2LambdaBar)
		    fLambdaAssocMassPtRapEmbeded2[curCentBin]->Fill(effLambda);
		}

	      } // End selection in the dca to prim. vtx and the number of clusters

	      // Distributions for the efficiency (Systematic checks)
	      /*
	      if( TMath::Abs(rapAs)<fYMax ){ 

		//  Cut in the DCA ToPrim Vtx
		if( (nClsTPCPos>fDaugNClsTPC) && (nClsTPCNeg>fDaugNClsTPC) ){
		  Double_t effLambdadcaPV[5] = {v0->MassLambda(),ptAs,rapAs,dcaPos,dcaNeg};
		  fLambdaAssocMassPtDCAPVEmbeded[curCentBin]->Fill(effLambdadcaPV);
		}		  

		// cut in the number of tpc clusters
		if( (dcaPos>fDCAToPrimVtx) && (dcaNeg>fDCAToPrimVtx) ){
		  Double_t effLambdaTPCcls[5]  = {v0->MassLambda(),ptAs,rapAs,nClsTPCPos,nClsTPCNeg};
		  fLambdaAssocMassPtDaugNClsTPCEmbeded[curCentBin]->Fill(effLambdaTPCcls);
		}

	      } // End selection for systematics
	      */
	    }  // End embeded particle selection
	    
	  } // End Lambda selection
	  
	  // AntiLambda:
	  if (ctAL && lCheckMcAntiLambda  && (TMath::Abs(nsigNegProton)<fNSigma) ){
	    
	    if(isNaturalPart){

	      if( (dcaPos>fDCAToPrimVtx) && (dcaNeg>fDCAToPrimVtx) && (nClsTPCPos>fDaugNClsTPC) && (nClsTPCNeg>fDaugNClsTPC) && (dca<fMaxDCADaughter) ){

		fAntiLambdaAssocPt->Fill(ptAs);
		fAntiLambdaAssocPtRap->Fill(ptAs,rapAs,centrality);
		fAntiLambdaAssocPtPhiEta[curCentBin]->Fill(p0->Phi(),etaAs,ptAs);
  
		// Rapidity cut
		if(TMath::Abs(rapAs)<fYMax)  {
  
		  Double_t effAntiLambda[3] = {v0->MassAntiLambda(),ptAs,rapAs};
		  /*Double_t effAntiLambdaVtx[4] = {v0->MassAntiLambda(),ptAs,rapAs,zv};
		  Double_t effAntiLambdaDCA[4] = {v0->MassAntiLambda(),ptAs,rapAs,dca};
		  Double_t effAntiLambdaCPA[4] = {v0->MassAntiLambda(),ptAs,rapAs,cpa};
		  Double_t effAntiLambdaShTPCcls[5] = {v0->MassAntiLambda(),ptAs,rapAs,fracPosDaugTPCSharedMap,fracNegDaugTPCSharedMap};
		  Double_t effAntiLambdaDaugPt[5] = {v0->MassAntiLambda(),ptAs,rapAs,lPtPos,lPtNeg};
		  Double_t effAntiLambdaCtau[4]   = {v0->MassAntiLambda(),ptAs,rapAs,dlL};
		  Double_t effAntiLambdaFidVol[4] = {v0->MassAntiLambda(),ptAs,rapAs,lt};*/

		  // Distributions for the efficiency (systematics chechks)
		  fAntiLambdaAssocMassPtRap[curCentBin]->Fill(effAntiLambda);
		  /*fAntiLambdaAssocMassPtVtx[curCentBin]->Fill(effAntiLambdaVtx);
		  fAntiLambdaAssocMassPtDCADaug[curCentBin]->Fill(effAntiLambdaDCA);
		  fAntiLambdaAssocMassPtCPA[curCentBin]->Fill(effAntiLambdaCPA);
		  fAntiLambdaAssocMassPtShTPCcls[curCentBin]->Fill(effAntiLambdaShTPCcls);
		  fAntiLambdaAssocMassPtDaugPt[curCentBin]->Fill(effAntiLambdaDaugPt);
		  fAntiLambdaAssocMassPtCtau[curCentBin]->Fill(effAntiLambdaCtau);
		  fAntiLambdaAssocMassPtFidVolume[curCentBin]->Fill(effAntiLambdaFidVol);*/


		  if( !isCandidate2K0s && !isCandidate2Lambda )
		    fAntiLambdaAssocMassPtRap2[curCentBin]->Fill(effAntiLambda);
		

		  fAntiLambdaMCResEta->Fill(resEta,pt,centrality);
		  fAntiLambdaMCResPhi->Fill(resPhi,pt,centrality);
		  fAntiLambdaMCResPt->Fill(resPt,pt,centrality);
		
		  fAntiLambdaPosMCResEta->Fill(resEtaPosDaug,pt,centrality);
		  fAntiLambdaPosMCResPhi->Fill(resPhiPosDaug,pt,centrality);
		  fAntiLambdaPosMCResPt->Fill(resPtPosDaug,pt,centrality);

		  fAntiLambdaNegMCResEta->Fill(resEtaNegDaug,pt,centrality);
		  fAntiLambdaNegMCResPhi->Fill(resPhiNegDaug,pt,centrality);
		  fAntiLambdaNegMCResPt->Fill(resPtNegDaug,pt,centrality);

		}

	      } // End selection in the dca to prim. vtx and the number of clusters

	      // Distributions for the efficiency (Systematic checks)
	      /*
	      if( TMath::Abs(rapAs)<fYMax ){ 
		
		//  Cut in the DCA ToPrim Vtx
		if( (nClsTPCPos>fDaugNClsTPC) && (nClsTPCNeg>fDaugNClsTPC) ){
		  
		  Double_t effAntiLambdadcaPV[5] = {v0->MassAntiLambda(),ptAs,rapAs,dcaPos,dcaNeg};
		  fAntiLambdaAssocMassPtDCAPV[curCentBin]->Fill(effAntiLambdadcaPV);
		}		  

		// cut in the number of tpc clusters
		if( (dcaPos>fDCAToPrimVtx) && (dcaNeg>fDCAToPrimVtx) && TMath::Abs(rapAs)<fYMax){
		  Double_t effAntiLambdaTPCcls[5]  = {v0->MassAntiLambda(),ptAs,rapAs,nClsTPCPos,nClsTPCNeg};
		  fAntiLambdaAssocMassPtDaugNClsTPC[curCentBin]->Fill(effAntiLambdaTPCcls);
		}

	      } // End selection for systematics
	      */
	      
	    }  // End natural particle selection
	    // Embeded particles
	    if(!isNaturalPart){

	      if( (dcaPos>fDCAToPrimVtx) && (dcaNeg>fDCAToPrimVtx) && (nClsTPCPos>fDaugNClsTPC) && (nClsTPCNeg>fDaugNClsTPC) && (dca<fMaxDCADaughter) ){

		if( TMath::Abs(rapAs)<fYMax ){

		  Double_t effAntiLambda[3] = {v0->MassAntiLambda(),ptAs,rapAs};
		  /* Double_t effAntiLambdaVtx[4] = {v0->MassAntiLambda(),ptAs,rapAs,zv};
		  Double_t effAntiLambdaDCA[4] = {v0->MassAntiLambda(),ptAs,rapAs,dca};
		  Double_t effAntiLambdaCPA[4] = {v0->MassAntiLambda(),ptAs,rapAs,cpa};
		  Double_t effAntiLambdaShTPCcls[5] = {v0->MassAntiLambda(),ptAs,rapAs,fracPosDaugTPCSharedMap,fracNegDaugTPCSharedMap};
		  Double_t effAntiLambdaDaugPt[5] = {v0->MassAntiLambda(),ptAs,rapAs,lPtPos,lPtNeg};
		  Double_t effAntiLambdaCtau[4]   = {v0->MassAntiLambda(),ptAs,rapAs,dlL};
		  Double_t effAntiLambdaFidVol[4] = {v0->MassAntiLambda(),ptAs,rapAs,lt}; */

		  // Distributions for the efficiency (systematics chechks)
		  fAntiLambdaAssocMassPtRapEmbeded[curCentBin]->Fill(effAntiLambda);
		  /* fAntiLambdaAssocMassPtVtxEmbeded[curCentBin]->Fill(effAntiLambdaVtx);
		  fAntiLambdaAssocMassPtDCADaugEmbeded[curCentBin]->Fill(effAntiLambdaDCA);
		  fAntiLambdaAssocMassPtCPAEmbeded[curCentBin]->Fill(effAntiLambdaCPA);
		  fAntiLambdaAssocMassPtShTPCclsEmbeded[curCentBin]->Fill(effAntiLambdaShTPCcls);
		  fAntiLambdaAssocMassPtDaugPtEmbeded[curCentBin]->Fill(effAntiLambdaDaugPt);
		  fAntiLambdaAssocMassPtCtauEmbeded[curCentBin]->Fill(effAntiLambdaCtau);
		  fAntiLambdaAssocMassPtFidVolumeEmbeded[curCentBin]->Fill(effAntiLambdaFidVol);*/

		  if( !isCandidate2K0s && !isCandidate2Lambda )
		    fAntiLambdaAssocMassPtRapEmbeded2[curCentBin]->Fill(effAntiLambda);
		}

	      } // End selection in the dca to prim. vtx and the number of clusters


	      // Distributions for the efficiency (Systematic checks)
	      /*
	      if( TMath::Abs(rapAs)<fYMax ){ 

		//  Cut in the DCA ToPrim Vtx
		if( (nClsTPCPos>fDaugNClsTPC) && (nClsTPCNeg>fDaugNClsTPC) ){

		  Double_t effAntiLambdadcaPV[5] = {v0->MassAntiLambda(),ptAs,rapAs,dcaPos,dcaNeg};
		  fAntiLambdaAssocMassPtDCAPVEmbeded[curCentBin]->Fill(effAntiLambdadcaPV);
		}		  

		// cut in the number of tpc ckusters
		if( (dcaPos>fDCAToPrimVtx) && (dcaNeg>fDCAToPrimVtx) ){

		  Double_t effAntiLambdaTPCcls[5]  = {v0->MassAntiLambda(),ptAs,rapAs,nClsTPCPos,nClsTPCNeg};
		  fAntiLambdaAssocMassPtDaugNClsTPCEmbeded[curCentBin]->Fill(effAntiLambdaTPCcls);
		}

	      } // End selection for systematics
	      */
	    }  // End embeded particle selection
	  
	  } // End AntiLambda
	  // Xi decay:
	  if( lComeFromXi && isNaturalPart ){
	    if(lPDGCodeV0==3122) { fLambdaAssocFromXi->Fill(ptAs,centrality); }
	    else if(lPDGCodeV0==-3122) { fAntiLambdaAssocFromXi->Fill(ptAs,centrality); }
	  }

	} // End Primary V0 selection
       
	// After the kinematical selection of K0s and Lambdas
	// it might be that the daugthers are not identified through MC Association
	if(lMCAssocNegDaug==0)
	  lMCAssocNegDaug = 6;
	if(lMCAssocPosDaug==0)
	  lMCAssocPosDaug = 6;
		
      } // End MC-Association 
      
    }// End Correlation Step
   
    // ************************************
    
  noas:    

    // Comparing the pt of the trigger particle wrt the v0-candidate's daughter:
    // It is used as well for the side-band subtraction
    Int_t isSameTrkPosDaug = -1;
    Int_t isSameTrkNegDaug = -1;
    if( step==kTriggerCheck ){
      isSameTrkPosDaug = SameTrack(trkTrig,ptrack);
      isSameTrkNegDaug = SameTrack(trkTrig,ntrack);
    }

    // *******************
    //   K0s selection
    // *******************
    if ( ctK && (TMath::Abs(rapK0s)<fYMax) && ( lPtArmV0 > TMath::Abs(0.2*lAlphaV0) ) && (dcaPos>fDCAToPrimVtx) && (dca<fMaxDCADaughter) &&
	 (dcaNeg>fDCAToPrimVtx) && (nClsTPCPos>fDaugNClsTPC) && (nClsTPCNeg>fDaugNClsTPC)  && (massK0s > 0.3979) &&
	 (massK0s < 0.5981) && lCheckMcK0Short ) {
      
      switch(step) {
      case kTriggerCheck: 

	if (isCandidate2K0s){ // Invariant mass cut in +- 3sigma(pt)

	  if(pt>ptTrig){
	    fIsV0LP = 1; 
	    fPtV0LP = pt;
	  }       
	  
	  if( isSameTrkPosDaug==1 || isSameTrkNegDaug==1){
	    Printf("  The LP has the same momentum in X and Y as one of the K0s daughters *** iV0 %d",iV0); 
	    
	    if(fCheckIDTrig){  // Compare properties of daughters nad 
	      Float_t difNegP[3];
	      difNegP[0] = (pTrig[0]-pNegDaug[0])/pTrig[0];  difNegP[1] = (pTrig[1]-pNegDaug[1])/pTrig[1]; difNegP[2] = (pTrig[2]-pNegDaug[2])/pTrig[2];
	      Float_t difPosP[3]; 
	      difPosP[0] = (pTrig[0]-pPosDaug[0])/pTrig[0];  difPosP[1] = (pTrig[1]-pPosDaug[1])/pTrig[1]; difPosP[2] = (pTrig[2]-pPosDaug[2])/pTrig[2];
	      Float_t posDeltaPhi =  phiTrig - phiPos, negDeltaPhi =  phiTrig - phiNeg;
	      Float_t posDeltaEta =  etaTrig - etaPos, negDeltaEta =  etaTrig - etaNeg;
	      
	      /*
		Printf("  The LP has the same momentum in X and Y as one of the K0s daughters *** iV0 %d \n\t\t %d %d %d \n\t\t %lf %lf %lf \n\t\t %lf %lf %lf \n\t\t %lf %lf \n\t\t %lf %lf ",
		iV0, TMath::Abs( trkTrig->GetID() ), ntrack->GetID() ,  ptrack->GetID() ,
		TMath::Abs( difNegP[1] ), TMath::Abs( difNegP[2] ), TMath::Abs( difNegP[0] ),
		TMath::Abs( difPosP[1] ), TMath::Abs( difPosP[2] ), TMath::Abs( difPosP[0] ),
		TMath::Abs( negDeltaPhi ), TMath::Abs( negDeltaEta ),
		TMath::Abs( posDeltaPhi ), TMath::Abs( posDeltaEta )
		);
	      */
	      
	      // Positive daughter
	      if( isSameTrkPosDaug==1 ){
		for(Int_t i=0;i<3;i++)
		  fCheckIDTrigPtK0s->Fill(difPosP[i],i,pt); 
		fCheckIDTrigPhiK0s->Fill(posDeltaPhi,0.,pt);
		fCheckIDTrigEtaK0s->Fill(posDeltaEta,0.,pt);
		fCheckIDTrigNclsK0s->Fill(nClsTPCPos,0.,pt);
	      }
	      // Negative daughter
	      if( isSameTrkNegDaug==1 ){ 
		for(Int_t i=0;i<3;i++)
		  fCheckIDTrigPtK0s->Fill(difNegP[i],i+3,pt); 
		fCheckIDTrigPhiK0s->Fill(negDeltaPhi,2.,pt);
		fCheckIDTrigEtaK0s->Fill(negDeltaEta,2.,pt);
		fCheckIDTrigNclsK0s->Fill(nClsTPCNeg,2.,pt);
	      }
	      
	    } // End check ID
	    
	    
	    fTriggerParticles->RemoveAt(iArray);
	    fTriggerParticles->AddAt( new AliMiniParticle(centrality, zv, idTrig, ptTrig, phiTrig, etaTrig, 0, 0, 0), iArray);
	    

	  } // Close isTrigFromV0daug
	  
	}// End K0s Mass cut
	
	break; // End K0s selection for TriggerCheck
      case kReconstruction:
	
	if(pt<10.){
	  
	  if(isNaturalPart) fK0sMass->Fill(massK0s,pt,centrality);
	  else fK0sMassEmbeded->Fill(massK0s,pt,centrality);
	  
	  fK0sMassPtEta->Fill(massK0s,pt,lEta);
	  fK0sMassPtRap[curCentBin]->Fill(massK0s,pt,rapK0s);
	  fK0sMassPtPhi->Fill(massK0s,pt,lPhi);

	  
	  if( (pt>kPtBinV0[0]) && (pt<kPtBinV0[kN1]) && isNaturalPart )
	    fAssocParticles->Add( new AliMiniParticle(centrality, zv, iV0, pt, lPhi, lEta, lMCAssocNegDaug, lMCAssocPosDaug, 3) );
	  
	  fK0sPosDaugFracShTPCcls->Fill(massK0s,pt,fracPosDaugTPCSharedMap);
	  fK0sNegDaugFracShTPCcls->Fill(massK0s,pt,fracNegDaugTPCSharedMap);

	}

	if( fDoQA && lCheckMcK0Short && isNaturalPart && (pt<10.) ){ // Quality Assurance

	  // Invariant Mass cut
	  if (TMath::Abs(mK0s-massK0s) < 3*sK0s) {
	  
	    fK0sDCAPosDaug->Fill(dcaPos,pt);
	    fK0sDCANegDaug->Fill(dcaNeg,pt);

	    fK0sPtPosDaug->Fill(pt,lPtPos);
	    fK0sPtNegDaug->Fill(pt,lPtNeg);

	    fK0sPhiEtaPosDaug->Fill(phiPos,etaPos,pt);
	    fK0sPhiEtaNegDaug->Fill(phiNeg,etaNeg,pt);
	    
	    fK0sDecayPos->Fill(dx,dy,pt);
	    fK0sDecayVertex->Fill(lt,pt);
	    
	    fK0sCPA->Fill(cpa,pt); 
	    fK0sDCAV0Daug->Fill(dca,pt); 

	    fK0sNClustersITSPos->Fill(phiPos,posITSNcls,pt);
	    fK0sNClustersITSNeg->Fill(phiNeg,negITSNcls,pt);

	    fK0sCTau->Fill(dlK,pt); 
	      

	    fK0sNClustersTPC->Fill(phiPos,nClsTPCPos,pt);
	    fK0sNClustersTPC->Fill(phiNeg,nClsTPCNeg,-pt);
	    

	  } // End selection in mass

	  if( TMath::Abs(mK0s-massK0s + 6.5*sK0s) < 1.5*sK0s ||
	      TMath::Abs(mK0s-massK0s - 6.5*sK0s) < 1.5*sK0s  ) {

	    fK0sBckgDCAPosDaug->Fill(dcaPos,pt);
	    fK0sBckgDCANegDaug->Fill(dcaNeg,pt);
	    	    
	    fK0sBckgPtPosDaug->Fill(pt,lPtPos);
	    fK0sBckgPtNegDaug->Fill(pt,lPtNeg);
	      
	    fK0sBckgPhiEtaPosDaug->Fill(phiPos,etaPos,pt);
	    fK0sBckgPhiEtaNegDaug->Fill(phiNeg,etaNeg,pt);
	      
	    fK0sBckgDecayPos->Fill(dx,dy,pt);
	    fK0sBckgDecayVertex->Fill(lt,pt);
	      
	    fK0sBckgCPA->Fill(cpa,pt); 
	    fK0sBckgDCAV0Daug->Fill(dca,pt); 
	      
	    fK0sBckgNClustersITSPos->Fill(phiPos,posITSNcls,pt);
	    fK0sBckgNClustersITSNeg->Fill(phiNeg,negITSNcls,pt);

	    fK0sBckgCTau->Fill(dlK,pt); 
	      
	    fK0sBckgNClustersTPC->Fill(phiPos,nClsTPCPos,pt);
	    fK0sBckgNClustersTPC->Fill(phiNeg,nClsTPCNeg,-pt);

	  }// End selection in outside the mass cut

	} // End QA
	
	break; // End K0s selection for Corrleation
      default:
	Printf( " Selection of 'step' is not set properly");
	break;
	
      }// End switch

    } // End K0s selection
      
    // *******************
    // Lambda selection
    // *******************
    if ( ctL && (TMath::Abs(rapLambda)<fYMax) && (dcaPos>fDCAToPrimVtx) && (dca<fMaxDCADaughter) &&
	 (dcaNeg>fDCAToPrimVtx) && (nClsTPCPos>fDaugNClsTPC) && (nClsTPCNeg>fDaugNClsTPC) &&
	 (massLambda > 1.0649) && (massLambda < 1.1651) && (TMath::Abs(nsigPosProton)<fNSigma) && lCheckMcLambda ){

      switch(step) {
      case kTriggerCheck: 
	
	if (isCandidate2Lambda && !isCandidate2K0s && !isCandidate2LambdaBar ){

	  if(pt>ptTrig) {
	    fIsV0LP = 1;
	    fPtV0LP = pt;
	  }

	  if( isSameTrkPosDaug==1 || isSameTrkNegDaug==1 ){
	    Printf("  The LP has the same momentum in X and Y as one of the Lambda daughters *** iV0 %d",iV0); 

	    if(fCheckIDTrig){  // Compare properties of daughters nad 
	      Float_t difNegP[3];
	      difNegP[0] = (pTrig[0]-pNegDaug[0])/pTrig[0];  difNegP[1] = (pTrig[1]-pNegDaug[1])/pTrig[1]; difNegP[2] = (pTrig[2]-pNegDaug[2])/pTrig[2];
	      Float_t difPosP[3]; 
	      difPosP[0] = (pTrig[0]-pPosDaug[0])/pTrig[0];  difPosP[1] = (pTrig[1]-pPosDaug[1])/pTrig[1]; difPosP[2] = (pTrig[2]-pPosDaug[2])/pTrig[2];
	      Float_t posDeltaPhi =  phiTrig - phiPos, negDeltaPhi =  phiTrig - phiNeg;
	      Float_t posDeltaEta =  etaTrig - etaPos, negDeltaEta =  etaTrig - etaNeg;
    
	      /*
		Printf("  The LP has the same momentum in X and Y as one of the Lambda daughters *** iV0 %d \n\t\t %d %d %d \n\t\t %lf %lf %lf \n\t\t %lf %lf %lf \n\t\t %lf %lf \n\t\t %lf %lf ",
		iV0, TMath::Abs( trkTrig->GetID() ), ntrack->GetID() ,  ptrack->GetID() ,
		TMath::Abs( difNegP[1] ), TMath::Abs( difNegP[2] ), TMath::Abs( difNegP[0] ),
		TMath::Abs( difPosP[1] ), TMath::Abs( difPosP[2] ), TMath::Abs( difPosP[0] ),
		TMath::Abs( negDeltaPhi ), TMath::Abs( negDeltaEta ),
		TMath::Abs( posDeltaPhi ), TMath::Abs( posDeltaEta )
		);
	      */

	      // Positive daughter
	      if( isSameTrkPosDaug==1 ){
		for(Int_t i=0;i<3;i++)
		  fCheckIDTrigPtLambda->Fill(difPosP[i],i,pt); 
		fCheckIDTrigPhiLambda->Fill(posDeltaPhi,0.,pt);
		fCheckIDTrigEtaLambda->Fill(posDeltaEta,0.,pt);
		fCheckIDTrigNclsLambda->Fill(nClsTPCPos,0.,pt);
	      }
	      // Negative daughter
	      if( isSameTrkNegDaug==1 ){ 
		for(Int_t i=0;i<3;i++)
		  fCheckIDTrigPtLambda->Fill(difNegP[i],i+3,pt); 
		fCheckIDTrigPhiLambda->Fill(negDeltaPhi,2.,pt);
		fCheckIDTrigEtaLambda->Fill(negDeltaEta,2.,pt);
		fCheckIDTrigNclsLambda->Fill(nClsTPCNeg,2.,pt);
	      }

	    } // End check ID

	    fTriggerParticles->RemoveAt(iArray);
	    fTriggerParticles->AddAt( new AliMiniParticle(centrality, zv, idTrig, ptTrig, phiTrig, etaTrig, 0, 0, 0), iArray);

	  } // Close isTrigFromV0daug

	} // End Lambda Mass cut	
	break; // End Lambda selection for TriggerCheck
      case kReconstruction:
	
	if(pt<10.){

	  if(isNaturalPart) fLambdaMass->Fill(massLambda,pt,centrality);
	  else  fLambdaMassEmbeded->Fill(massLambda,pt,centrality);

	  if( !isCandidate2K0s && !isCandidate2LambdaBar){
	    if(isNaturalPart) fLambdaMass2->Fill(massLambda,pt,centrality);
	    else fLambdaMass2Embeded->Fill(massLambda,pt,centrality);
	  }

	  fLambdaMassPtEta->Fill(massLambda,pt,lEta);
	  fLambdaMassPtRap[curCentBin]->Fill(massLambda,pt,rapLambda);	
	  fLambdaMassPtPhi->Fill(massLambda,pt,lPhi);

	  
	  if( (pt>kPtBinV0[0]) && (pt<kPtBinV0[kN1]) && isNaturalPart )
	    fAssocParticles->Add( new AliMiniParticle(centrality, zv, iV0, pt, lPhi, lEta, lMCAssocNegDaug, lMCAssocPosDaug, 4) );
	  
	  
	  fLambdaPosDaugFracShTPCcls->Fill(massLambda,pt,fracPosDaugTPCSharedMap);
	  fLambdaNegDaugFracShTPCcls->Fill(massLambda,pt,fracNegDaugTPCSharedMap);

	}
	
	// Invariant Mass cut
	if(fDoQA && lCheckMcLambda && isNaturalPart && (pt<10.)){ // Quality Assurance
          
	  // Invariant Mass cut
	  if (TMath::Abs(mLambda-massLambda) < 3*sL) {
	   
	    fLambdaDCAPosDaug->Fill(dcaPos,pt);
	    fLambdaDCANegDaug->Fill(dcaNeg,pt);
	    
	    fLambdaPtPosDaug->Fill(pt,lPtPos);
	    fLambdaPtNegDaug->Fill(pt,lPtNeg);

	    fLambdaPhiEtaPosDaug->Fill(phiPos,etaPos,pt);
	    fLambdaPhiEtaNegDaug->Fill(phiNeg,etaNeg,pt);

	    fLambdaDecayPos->Fill(dx,dy,pt);
	    fLambdaDecayVertex->Fill(lt,pt);

	    fLambdaCPA->Fill(cpa,pt); 
	    fLambdaDCAV0Daug->Fill(dca,pt); 

	    fLambdaNClustersITSPos->Fill(phiPos,posITSNcls,pt);
	    fLambdaNClustersITSNeg->Fill(phiNeg,negITSNcls,pt);

	    fLambdaCTau->Fill(dlL,pt); 
	     
	    fLambdaNClustersTPC->Fill(phiPos,nClsTPCPos,pt);
	    fLambdaNClustersTPC->Fill(phiNeg,nClsTPCNeg,-pt);

	  } // End selection in mass
	
	  if( (TMath::Abs(mLambda-massLambda + 6.5*sL) < 1.5*sL) ||
	      (TMath::Abs(mLambda-massLambda - 6.5*sL) < 1.5*sL) ){

	    fLambdaBckgDCAPosDaug->Fill(dcaPos,pt);
	    fLambdaBckgDCANegDaug->Fill(dcaNeg,pt);
	 
	    fLambdaBckgPtPosDaug->Fill(pt,lPtPos);
	    fLambdaBckgPtNegDaug->Fill(pt,lPtNeg);

	    fLambdaBckgPhiEtaPosDaug->Fill(phiPos,etaPos,pt);
	    fLambdaBckgPhiEtaNegDaug->Fill(phiNeg,etaNeg,pt);
	      
	    fLambdaBckgDecayPos->Fill(dx,dy,pt);
	    fLambdaBckgDecayVertex->Fill(lt,pt);
	      
	    fLambdaBckgCPA->Fill(cpa,pt); 
	    fLambdaBckgDCAV0Daug->Fill(dca,pt); 

	    fLambdaBckgNClustersITSPos->Fill(phiPos,posITSNcls,pt);
	    fLambdaBckgNClustersITSNeg->Fill(phiNeg,negITSNcls,pt);

	    fLambdaBckgCTau->Fill(dlL,pt); 
	     	      
	    fLambdaBckgNClustersTPC->Fill(phiPos,nClsTPCPos,pt);
	    fLambdaBckgNClustersTPC->Fill(phiNeg,nClsTPCNeg,-pt);
	    
	  }// End selection in outside the mass cut  
	  
	} // End QA

	break; // End Lambda selection for Correlation
      default:
	Printf(" Selection of 'step' is not set properly");
	break;
	
      }// End switch
      
    } // End Lambda selection

    // *******************
    // AntiLambda selection
    // *******************
    if ( ctAL && (TMath::Abs(rapLambda)<fYMax)  && (dcaPos>fDCAToPrimVtx) && (dca<fMaxDCADaughter) &&
	 (dcaNeg>fDCAToPrimVtx) && (nClsTPCPos>fDaugNClsTPC) && (nClsTPCNeg>fDaugNClsTPC) &&
	 (massAntiLambda > 1.0649) && ( massAntiLambda < 1.1651 ) && (TMath::Abs(nsigNegProton)<fNSigma) && lCheckMcAntiLambda ) {
      
      switch(step) {
      case kTriggerCheck: 
	
	if( isCandidate2LambdaBar && !isCandidate2K0s && !isCandidate2Lambda ){

	  if(pt>ptTrig) {
	    fIsV0LP = 1;
	    fPtV0LP = pt;
	  }
	  
	  if( isSameTrkPosDaug==1 || isSameTrkNegDaug==1 ){
	    Printf("  The LP has the same momentum in X and Y as one of the AntiLambda daughters *** iV0 %d",iV0); 

	    if(fCheckIDTrig){  // Compare properties of daughters nad 
	      Float_t difNegP[3];
	      difNegP[0] = (pTrig[0]-pNegDaug[0])/pTrig[0];  difNegP[1] = (pTrig[1]-pNegDaug[1])/pTrig[1]; difNegP[2] = (pTrig[2]-pNegDaug[2])/pTrig[2];
	      Float_t difPosP[3]; 
	      difPosP[0] = (pTrig[0]-pPosDaug[0])/pTrig[0];  difPosP[1] = (pTrig[1]-pPosDaug[1])/pTrig[1]; difPosP[2] = (pTrig[2]-pPosDaug[2])/pTrig[2];
	      Float_t posDeltaPhi =  phiTrig - phiPos, negDeltaPhi =  phiTrig - phiNeg;
	      Float_t posDeltaEta =  etaTrig - etaPos, negDeltaEta =  etaTrig - etaNeg;

	      /*
		Printf("  The LP has the same momentum in X and Y as one of the AntiLambda daughters *** iV0 %d \n\t\t %d %d %d \n\t\t %lf %lf %lf \n\t\t %lf %lf %lf \n\t\t %lf %lf \n\t\t %lf %lf ",
		iV0, TMath::Abs( trkTrig->GetID() ), ntrack->GetID() ,  ptrack->GetID() ,
		TMath::Abs( difNegP[1] ), TMath::Abs( difNegP[2] ), TMath::Abs( difNegP[0] ),
		TMath::Abs( difPosP[1] ), TMath::Abs( difPosP[2] ), TMath::Abs( difPosP[0] ),
		TMath::Abs( negDeltaPhi ), TMath::Abs( negDeltaEta ),
		TMath::Abs( posDeltaPhi ), TMath::Abs( posDeltaEta )
		);
	      */

	      // Positive daughter
	      if( isSameTrkPosDaug==1 ){
		for(Int_t i=0;i<3;i++)
		  fCheckIDTrigPtAntiLambda->Fill(difPosP[i],i,pt); 
		fCheckIDTrigPhiAntiLambda->Fill(posDeltaPhi,0.,pt);
		fCheckIDTrigEtaAntiLambda->Fill(posDeltaEta,0.,pt);
		fCheckIDTrigNclsAntiLambda->Fill(nClsTPCPos,0.,pt);
	      }
	      // Negative daughter
	      if( isSameTrkNegDaug==1 ){ 
		for(Int_t i=0;i<3;i++)
		  fCheckIDTrigPtAntiLambda->Fill(difNegP[i],i+3,pt); 
		fCheckIDTrigPhiAntiLambda->Fill(negDeltaPhi,2.,pt);
		fCheckIDTrigEtaAntiLambda->Fill(negDeltaEta,2.,pt);
		fCheckIDTrigNclsAntiLambda->Fill(nClsTPCNeg,2.,pt);
	      }

	    } // End check ID  


	    fTriggerParticles->RemoveAt(iArray);
	    fTriggerParticles->AddAt( new AliMiniParticle(centrality, zv, idTrig, ptTrig, phiTrig, etaTrig, 0, 0, 0), iArray);

	  }// Close isTrigFromV0daug
	  
	}// End AntiLambda Mass cut
	break; // End AntiLambda selection for CheckTrigger
      case kReconstruction: 
	
	if( (pt<10.) ) {

	  if(isNaturalPart)  fAntiLambdaMass->Fill(massAntiLambda,pt,centrality);
	  else fAntiLambdaMassEmbeded->Fill(massAntiLambda,pt,centrality);

	  if( !isCandidate2K0s && !isCandidate2Lambda) {
	    if(isNaturalPart) fAntiLambdaMass2->Fill(massAntiLambda,pt,centrality);
	    else fAntiLambdaMass2Embeded->Fill(massAntiLambda,pt,centrality);
	  }

	  fAntiLambdaMassPtEta->Fill(massAntiLambda,pt,lEta);
	  fAntiLambdaMassPtRap[curCentBin]->Fill(massAntiLambda,pt,rapLambda);	  
	  fAntiLambdaMassPtPhi->Fill(massAntiLambda,pt,lPhi);
	
	  
	  if( (pt>kPtBinV0[0]) && (pt<kPtBinV0[kN1]) && isNaturalPart )
	    fAssocParticles->Add( new AliMiniParticle(centrality, zv, iV0, pt, lPhi, lEta, lMCAssocNegDaug, lMCAssocPosDaug, 5) );
	 
	  fAntiLambdaPosDaugFracShTPCcls->Fill(massAntiLambda,pt,fracPosDaugTPCSharedMap);
	  fAntiLambdaNegDaugFracShTPCcls->Fill(massAntiLambda,pt,fracNegDaugTPCSharedMap);
 
	}
 
	if( fDoQA && lCheckMcAntiLambda && isNaturalPart && (pt<10.) ){ // Quality Assurance

	  // Invariant Mass cut
	  if (TMath::Abs(mLambda-massAntiLambda) < 3*sAL) {

	    fAntiLambdaDCAPosDaug->Fill(dcaPos,pt);
	    fAntiLambdaDCANegDaug->Fill(dcaNeg,pt);

	    fAntiLambdaPtPosDaug->Fill(pt,lPtPos);
	    fAntiLambdaPtNegDaug->Fill(pt,lPtNeg);
		  
	    fAntiLambdaPhiEtaPosDaug->Fill(phiPos,etaPos,pt);
	    fAntiLambdaPhiEtaNegDaug->Fill(phiNeg,etaNeg,pt);
		  
	    fAntiLambdaDecayPos->Fill(dx,dy,pt);
	    fAntiLambdaDecayVertex->Fill(lt,pt);
		  
	    fAntiLambdaCPA->Fill(cpa,pt); 
	    fAntiLambdaDCAV0Daug->Fill(dca,pt); 
		  
	    fAntiLambdaNClustersITSPos->Fill(phiPos,posITSNcls,pt);
	    fAntiLambdaNClustersITSNeg->Fill(phiNeg,negITSNcls,pt);

	    fAntiLambdaCTau->Fill(dlL,pt); 
	      	      
	    fAntiLambdaNClustersTPC->Fill(phiPos,nClsTPCPos,pt);
	    fAntiLambdaNClustersTPC->Fill(phiNeg,nClsTPCNeg,-pt);
	    
	  } // End selection in mass
	
	  if( (TMath::Abs(mLambda-massAntiLambda + 6.5*sAL) < 1.5*sAL) ||
	      (TMath::Abs(mLambda-massAntiLambda - 6.5*sAL) < 1.5*sAL) ){

	    fAntiLambdaBckgDCAPosDaug->Fill(dcaPos,pt);
	    fAntiLambdaBckgDCANegDaug->Fill(dcaNeg,pt);
	    	      
	    fAntiLambdaBckgPtPosDaug->Fill(pt,lPtPos);
	    fAntiLambdaBckgPtNegDaug->Fill(pt,lPtNeg);
	      
	    fAntiLambdaBckgPhiEtaPosDaug->Fill(phiPos,etaPos,pt);
	    fAntiLambdaBckgPhiEtaNegDaug->Fill(phiNeg,etaNeg,pt);
	      
	    fAntiLambdaBckgDecayPos->Fill(dx,dy,pt);
	    fAntiLambdaBckgDecayVertex->Fill(lt,pt);
	      
	    fAntiLambdaBckgCPA->Fill(cpa,pt); 
	    fAntiLambdaBckgDCAV0Daug->Fill(dca,pt); 

	    fAntiLambdaBckgNClustersITSPos->Fill(phiPos,posITSNcls,pt);
	    fAntiLambdaBckgNClustersITSNeg->Fill(phiNeg,negITSNcls,pt);

	    fAntiLambdaBckgCTau->Fill(dlL,pt); 
	   
	    fAntiLambdaBckgNClustersTPC->Fill(phiPos,nClsTPCPos,pt);
	    fAntiLambdaBckgNClustersTPC->Fill(phiNeg,nClsTPCNeg,-pt);

	  }// End selection in outside the mass cut
	  
	} // End QA
	
	break;
      default:
	Printf( " Selection of 'step' is not set properly");
	break;
      }// End switch
      
    } // End AntiLambda selection
    
  } // End V0 loop
  
}

//___________________________________________________________________________________________

void AliAnalysisTaskLambdaOverK0sJets::TriggerParticle() 
{ 
  // Obtain the trigger particles of the event to perform the correlations in phi and eta

  // ----------------------------
  // 1. Trigger particles 
  TClonesArray *stack = 0x0;
  if(fIsMC){  
    TList *lst = fAOD->GetList();
    stack = (TClonesArray*)lst->FindObject(AliAODMCParticle::StdBranchName());
    if (!stack) {
      Printf("ERROR: stack not available");
      return;
    }
  }

  Int_t nTrk= fAOD->GetNumberOfTracks();
  AliCentrality *cent = fAOD->GetCentrality();
  Float_t centrality = cent->GetCentralityPercentile("V0M");
  const AliAODVertex *vtx = fAOD->GetPrimaryVertex();
  Float_t zv=vtx->GetZ();
  Double_t pt  = -1000.;
  Double_t eta = -1000.;
  Double_t phi = -1000.;
  Float_t resPt = -1000.; 
  Float_t resEta = -1000.; 
  Float_t resPhi = -1000.;

  for (Int_t i=0; i<nTrk; i++) {
    const AliAODTrack *t = dynamic_cast<const AliAODTrack*>(fAOD->GetTrack(i));
    if(!t) AliFatal("Not a standard AOD");
    if(!AcceptTrack(t)) continue;
    pt=t->Pt();
    eta=t->Eta();
   
    if( (pt>fTrigPtMin)  && (pt<fTrigPtMax) &&  (TMath::Abs(eta)<fTrigEtaMax) ) {

      phi=t->Phi();
      fTriggerParticles->Add( new AliMiniParticle(centrality, zv, i, pt, phi, eta, 0, 0, 1) );    

      if(fIsMC){    
	Int_t lab = TMath::Abs(t->GetLabel());
	AliAODMCParticle *part=(AliAODMCParticle*)stack->UncheckedAt(lab);

	resPt  = (part->Pt()  - pt)/pt;	
	resEta = part->Eta() - eta;	
	resPhi = part->Phi() - phi;

	fTriggerMCResPt->Fill(resPt,pt,centrality);
	fTriggerMCResEta->Fill(resEta,pt,centrality);
	fTriggerMCResPhi->Fill(resPhi,pt,centrality);
      }

    }
  }

  // ----------------------------
  // 2. Checking if the trigger particle 
  // might be a daughter from the V0-candidate
  
  for (Int_t i=0; i<(fTriggerParticles->GetEntriesFast()); i++){
    AliMiniParticle* trig = (AliMiniParticle*) fTriggerParticles->At(i);
    Int_t id = trig->ID();
    V0Loop(kTriggerCheck,kFALSE,i,id);
  }
  
    
}

//___________________________________________________________________________________________

void AliAnalysisTaskLambdaOverK0sJets::UserExec(Option_t *)
{
  // Main loop for the Analysis

  // Initializing global variables for the correlation studies (mandatory for each event).
  // ---- 1) Trigger Particle: id track
  fIdTrigger  = -1;
  // ---- 2) TriggerCheck: Variables used to crosscheck if trigger particle is a V0 daughter ---- //
  fIsV0LP     = 0;
  fPtV0LP     = -10.;
  fIsSndCheck = 0;

  // Getting AOD Event
  fAOD = (AliAODEvent *)InputEvent();
  fEvents->Fill(0); //event counter  

  if (!fAOD) {
    Printf("ERROR: aod not available");
    return;
  }
  fEvents->Fill(1);
  
  // Physics selection
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *hdr=(AliInputEventHandler*)mgr->GetInputEventHandler();
  UInt_t maskIsSelected = hdr->IsEventSelected();
  Bool_t isSelected = kFALSE;

  Bool_t isSelectedCentral     = kFALSE;
  Bool_t isSelectedSemiCentral = kFALSE;
  Bool_t isSelectedMB          = kFALSE;
  if( fCollision.Contains("PbPb2010") )
    isSelected = (maskIsSelected & AliVEvent::kMB);
  else if( fCollision.Contains("PbPb2011") ){
    isSelectedCentral   =  maskIsSelected  &  AliVEvent::kCentral;
    isSelectedSemiCentral   =  maskIsSelected  &  AliVEvent::kSemiCentral;
    isSelectedMB   =  maskIsSelected   &  AliVEvent::kMB;
    if( isSelectedCentral || isSelectedSemiCentral || isSelectedMB ) 	isSelected = kTRUE;
  }

  if(!isSelected) return;
  fEvents->Fill(2);

  // Centrality selection
  AliCentrality *cent = fAOD->GetCentrality();
  Float_t centrality = cent->GetCentralityPercentile("V0M");
  fCentrality->Fill(centrality);

  if (!cent->IsEventInCentralityClass(fCentMin,fCentMax,"V0M")) return;
  fEvents->Fill(3);

  Int_t curCentBin = CentBin(centrality);
  if(curCentBin!=-1) fEvtPerCent->Fill(0.,curCentBin);

  // Global primary vertex 
  const AliAODVertex *vtx = fAOD->GetPrimaryVertex();
  if (vtx->GetNContributors()<3) return;
  fEvents->Fill(4);
  if(curCentBin!=-1) fEvtPerCent->Fill(1,curCentBin);

  // SPD primary vertex 
  const AliAODVertex *vtxSPD = fAOD->GetPrimaryVertexSPD(); 
  if (vtxSPD->GetNContributors()<3) return;
  fEvents->Fill(5);
  if(curCentBin!=-1) fEvtPerCent->Fill(2,curCentBin);

  // Correlaiton between global Zvtx and SPD Zvtx
  Float_t zv=vtx->GetZ(), zvSPD=vtxSPD->GetZ();
  fPrimayVtxGlobalvsSPD->Fill(zv,zvSPD);
  
  Float_t xv=vtx->GetX(), yv=vtx->GetY();
  const Float_t priVtx[3] = {xv,yv,zv};

  if (TMath::Abs(zv) > 10.) return;   
  fEvents->Fill(6);
  if(curCentBin!=-1) fEvtPerCent->Fill(3,curCentBin);

  if( TMath::Abs( zv - zvSPD ) > 0.5) return;
  fEvents->Fill(7);
  if(curCentBin!=-1) fEvtPerCent->Fill(4,curCentBin);

  fPrimaryVertexX->Fill(xv);
  fPrimaryVertexY->Fill(yv);
  fPrimaryVertexZ->Fill(zv);
 
  fCentrality2->Fill(centrality);

  if(isSelectedCentral) fCentralityTrig->Fill(centrality,1);
  if(isSelectedSemiCentral) fCentralityTrig->Fill(centrality,2);
  if(isSelectedMB) fCentralityTrig->Fill(centrality,3);


  // Protect the code: only interested in events with centrality < 40%
  if(centrality>=40.) return;

  //  Event plane 
  AliEventplane *EventPlane = InputEvent()->GetEventplane();
  Double_t eventPlane = EventPlane->GetEventplane("V0",InputEvent(),2);
 
  eventPlane = ( (eventPlane < 0) ? eventPlane + TMath::Pi() : eventPlane );
  eventPlane = ( ( eventPlane > TMath::Pi() ) ? eventPlane - TMath::Pi() : eventPlane );
 
  // Magnetic field sign
  const Float_t bSign = (fAOD->GetMagneticField() > 0) ? 1 : -1;

  // Getting PID Response
  fPIDResponse = hdr->GetPIDResponse();

  Int_t curVtxBin = VtxBin(zv);
 
  // **********************************************
  // Multiplicity
  Float_t mult = GetMultiplicity();
  fChargedMultiplicity->Fill(mult,curCentBin);

  // **********************************************
  // Triggered Particle -  Trigger Particle
  fTriggerParticles = new TObjArray();
  fTriggerParticles->SetOwner(kTRUE);
  TriggerParticle(); 
 
  // V0-candidate is the highest particle in the event:
  if(fIsV0LP) { fEvents->Fill(8); fTriggerIsV0->Fill(fPtV0LP);}

  Int_t NtrigPerEvnt = 0;
  Float_t phi2 = -100.;
  for (Int_t i=0; i<(fTriggerParticles->GetEntriesFast()); i++){
    AliMiniParticle* trig = (AliMiniParticle*) fTriggerParticles->At(i);

    if(trig->WhichCandidate()==0){
      fTriggerComingFromDaug->Fill(trig->Pt());
      fCheckTriggerFromV0Daug->Fill(1);
      if(fIsV0LP)  fCheckTriggerFromV0Daug->Fill(2);
    }
    else if( trig->WhichCandidate()==1){
      fCheckTriggerFromV0Daug->Fill(0);
      fTriggerPtCentCh->Fill(trig->Pt(),centrality,zv);

      phi2 = ( (trig->Phi() > TMath::Pi()) ? trig->Phi() - TMath::Pi() : trig->Phi() )  ;
      fTriggerEventPlane->Fill(phi2);

      NtrigPerEvnt++;

    }

  }

  if(NtrigPerEvnt>0) fEvents->Fill(11); 
  else fEvents->Fill(12);

  fNTrigPerEvt->Fill(NtrigPerEvnt,centrality);
  
  // ******************************************
  // Start loop over MC particles
  fTriggerPartMC = new TObjArray();
  fTriggerPartMC->SetOwner(kTRUE);
  fAssocPartMC = new TObjArray();
  fAssocPartMC->SetOwner(kTRUE);

  fEndOfHijingEvent = -1;
  TClonesArray *stack = 0x0;
  Float_t mcXv=0., mcYv=0., mcZv=0.;
  
  if(fIsMC) {

    TList *lst = fAOD->GetList();
    stack = (TClonesArray*)lst->FindObject(AliAODMCParticle::StdBranchName());
    if (!stack) {
      Printf("ERROR: stack not available");
      return;
    }
    
    AliAODMCHeader *mcHdr = 
      (AliAODMCHeader*)lst->FindObject(AliAODMCHeader::StdBranchName());
  
    mcXv=mcHdr->GetVtxX(); mcYv=mcHdr->GetVtxY(); mcZv=mcHdr->GetVtxZ();
  
    Int_t nTrkMC = stack->GetEntriesFast();
    // -----------------------------------------
    // --------- Trigger particle --------------
    // -----------------------------------------
    for (Int_t iTrkMC = 0; iTrkMC < nTrkMC; iTrkMC++){      
      AliAODMCParticle *p0 = (AliAODMCParticle*)stack->At(iTrkMC);
      if(!p0) continue;

      // ----------------------------------------   
      // For injected MC: it determines where HIJING event ends 
      if (fEndOfHijingEvent==-1) { 
        if ( ( p0->GetStatus() == 21 ) ||
	     ( (p0->GetPdgCode() == 443) &&
	       (p0->GetMother() == -1)   &&
	       (p0->GetDaughterLabel(0) ==  (iTrkMC+1))) ) {
	  fEndOfHijingEvent = iTrkMC; 
        }
      }

      // ----------------------------------------     
      Int_t isNaturalPart = 1;
      if ( (iTrkMC>=fEndOfHijingEvent) && 
	   (fEndOfHijingEvent!=-1)     && 
	   (p0->GetMother()<0) ) 
	isNaturalPart = 0; 
     
      // ----------------------------------------   
      if(p0->Charge()==0) continue;
      if(isNaturalPart == 0) continue;
      if( !p0->IsPhysicalPrimary() ) continue;

      if(TMath::Abs(p0->Eta())>fTrigEtaMax) continue;
      if( ( p0->Pt() < fTrigPtMCMin )  || ( p0->Pt() > fTrigPtMCMax ) ) continue;

      fTriggerPartMC->Add( new AliMiniParticle(centrality, zv, iTrkMC, p0->Pt(), p0->Phi(), p0->Eta(), 0, 0, 1) ); 
    
    } // End loop over charged particles

    // -----------------------------------------
    // ---------- Strange particles ------------
    // -----------------------------------------
    //fEndOfHijingEvent = -1;
    for (Int_t iTrkMC = 0; iTrkMC < nTrkMC; iTrkMC++){     
      AliAODMCParticle *p0 = (AliAODMCParticle*)stack->At(iTrkMC);
      if(!p0) continue;
    
      // ----------------------------------------      
      Int_t lPdgcodeCurrentPart = p0->GetPdgCode();     
      if ( (lPdgcodeCurrentPart != kK0Short) &&
	   (lPdgcodeCurrentPart != kLambda0) &&
	   (lPdgcodeCurrentPart != kLambda0Bar) &&
	   //Adding Xi- and Xi0 particles 
	   (lPdgcodeCurrentPart != kXiMinus) &&
	   (lPdgcodeCurrentPart != 3322) ) continue;
      
      // ----------------------------------------
      Int_t isNaturalPart = 1;
      if ( (iTrkMC>=fEndOfHijingEvent) && 
	   (fEndOfHijingEvent!=-1)     && 
	   (p0->GetMother()<0) ) 
	isNaturalPart = 0; 
     
      if( lPdgcodeCurrentPart != kXiMinus )
	fInjectedParticles->Fill(isNaturalPart);

      if(fSeparateInjPart && !isNaturalPart) continue;

      // ----------------------------------------
      Float_t lRapCurrentPart = MyRapidity(p0->E(),p0->Pz());      
      Float_t lEtaCurrentPart = p0->Eta();
      Float_t lPhiCurrentPart = p0->Phi();
      Float_t lPtCurrentPart  = p0->Pt();

      if(fUseEtaCut){
	lRapCurrentPart = lEtaCurrentPart;
      }

      Int_t iCurrentMother = p0->GetMother();       
      AliAODMCParticle *pCurrentMother = (AliAODMCParticle *)stack->At(iCurrentMother);
      Int_t lPdgCurrentMother = 0;    
      if (iCurrentMother == -1) { lPdgCurrentMother = 0;}
      else { lPdgCurrentMother = pCurrentMother->GetPdgCode(); }

      Int_t id0  = p0->GetDaughterLabel(0);
      Int_t id1  = p0->GetDaughterLabel(1);
    
      //if ( id0 ==  id1 ) continue;
      if ( (id0 < 0 || id1 < 0) ||
	   (id0 >=nTrkMC  || id1 >= nTrkMC) ) continue;

      AliAODMCParticle *pDaughter0 = (AliAODMCParticle *)stack->UncheckedAt(id0);
      AliAODMCParticle *pDaughter1 = (AliAODMCParticle *)stack->UncheckedAt(id1);
      if (!pDaughter0 || !pDaughter1) continue; 
   
      /*
	if ( TMath::Abs(pDaughter0->Eta()) > fMaxEtaDaughter ||
	TMath::Abs(pDaughter1->Eta()) > fMaxEtaDaughter )
	continue;	
      */
      // Daughter momentum cut: ! FIX it in case of AOD !
      /*
	if ( ( pDaughter0->Pt() < fMinPtDaughter ) || 
	( pDaughter1->Pt() < fMinPtDaughter )  ) 
	continue;
      */
      
      if ((p0->Pt())<pMin || (p0->Pt())>100. ) continue;  
      if (TMath::Abs(lRapCurrentPart) > fYMax)  continue;
    
      Float_t dx = mcXv-p0->Xv(),  dy = mcYv-p0->Yv(),  dz = mcZv-p0->Zv();
      Float_t l = TMath::Sqrt(dx*dx + dy*dy + dz*dz);
      
      //Cut in the 3D-distance of the secondary vertex to primary vertex
      if (l > 0.01) continue; // secondary V0 
     
      //Transverse distance to vertex
      dx = mcXv-pDaughter0->Xv(); dy = mcYv-pDaughter0->Yv();
      //Float_t lt=TMath::Sqrt(dx*dx + dy*dy);

      // Pt Selection
      if((p0->Pt())<10.) { 

	// K0s
	if (lPdgcodeCurrentPart == kK0Short) {

	  fK0sMCPt->Fill(lPtCurrentPart);
	  fK0sMCPtRap->Fill(lPtCurrentPart,lRapCurrentPart,centrality); 

	  if(isNaturalPart){
	    fK0sMCPtRap2->Fill(lPtCurrentPart,lRapCurrentPart,centrality);
	    fK0sMCPtPhiEta[curCentBin]->Fill(lPhiCurrentPart,lEtaCurrentPart,lPtCurrentPart);
	  
	    if(TMath::Abs(lRapCurrentPart)<fYMax)  fK0sMCPtRapVtx[curCentBin]->Fill(lPtCurrentPart,lRapCurrentPart,zv);
	  
	    if( (lPtCurrentPart>kPtBinV0[0]) && (lPtCurrentPart<kPtBinV0[kN1]) && isNaturalPart )
	      fAssocPartMC->Add( new AliMiniParticle(centrality, zv, iTrkMC, lPtCurrentPart, lPhiCurrentPart, lEtaCurrentPart, 0, 0, 3) );
	  }
	  else{ 
	    fK0sMCPtRapEmbeded->Fill(lPtCurrentPart,lRapCurrentPart,centrality); 
	    if(TMath::Abs(lRapCurrentPart)<fYMax)  fK0sMCPtRapVtxEmbeded[curCentBin]->Fill(lPtCurrentPart,lRapCurrentPart,zv);
	  }

	} // End K0s selection
	// Lambda
	if (lPdgcodeCurrentPart == kLambda0) {
	
	  fLambdaMCPt->Fill(lPtCurrentPart);
	  fLambdaMCPtRap->Fill(lPtCurrentPart,lRapCurrentPart,centrality);
      
	  if(isNaturalPart){
	    fLambdaMCPtRap2->Fill(lPtCurrentPart,lRapCurrentPart,centrality);
	    fLambdaMCPtPhiEta[curCentBin]->Fill(lPhiCurrentPart,lEtaCurrentPart,lPtCurrentPart);
	    
	    if(TMath::Abs(lRapCurrentPart)<fYMax) fLambdaMCPtRapVtx[curCentBin]->Fill(lPtCurrentPart,lRapCurrentPart,zv);

	    if( (lPtCurrentPart>kPtBinV0[0]) && (lPtCurrentPart<kPtBinV0[kN1]) && isNaturalPart )
	      fAssocPartMC->Add( new AliMiniParticle(centrality, zv, iTrkMC, lPtCurrentPart, lPhiCurrentPart, lEtaCurrentPart, 0, 0, 4) );
	  }
	  else{ 
	    fLambdaMCPtRapEmbeded->Fill(lPtCurrentPart,lRapCurrentPart,centrality); 
	    fLambdaMCPtRapVtxEmbeded[curCentBin]->Fill(lPtCurrentPart,lRapCurrentPart,zv);
	  }

	  if ( isNaturalPart && TMath::Abs(lPdgCurrentMother) == 3312 ) 
	    fLambdaMCFromXi->Fill(lPtCurrentPart,centrality);
	
	} // End Lambda
	// AntiLambda
	if (lPdgcodeCurrentPart == kLambda0Bar) {

	  fAntiLambdaMCPt->Fill(lPtCurrentPart);
	  fAntiLambdaMCPtRap->Fill(lPtCurrentPart,lRapCurrentPart,centrality);

	  if(isNaturalPart){
	    fAntiLambdaMCPtRap2->Fill(lPtCurrentPart,lRapCurrentPart,centrality);	    
	    fAntiLambdaMCPtPhiEta[curCentBin]->Fill(lPhiCurrentPart,lEtaCurrentPart,lPtCurrentPart);

	    if(TMath::Abs(lRapCurrentPart)<fYMax) fAntiLambdaMCPtRapVtx[curCentBin]->Fill(lPtCurrentPart,lRapCurrentPart,zv);
	    
	    if( (lPtCurrentPart>kPtBinV0[0]) && (lPtCurrentPart<kPtBinV0[kN1]) && isNaturalPart )
	      fAssocPartMC->Add( new AliMiniParticle(centrality, zv, iTrkMC, lPtCurrentPart, lPhiCurrentPart, lEtaCurrentPart, 0, 0, 5) );
	  }
	  else{ 
	    fAntiLambdaMCPtRapEmbeded->Fill(lPtCurrentPart,lRapCurrentPart,centrality); 
	    if(TMath::Abs(lRapCurrentPart)<fYMax) fAntiLambdaMCPtRapVtxEmbeded[curCentBin]->Fill(lPtCurrentPart,lRapCurrentPart,zv);
	  }
	  
	  if ( isNaturalPart && TMath::Abs(lPdgCurrentMother) == 3312 ) 
	    fAntiLambdaMCFromXi->Fill(lPtCurrentPart,centrality);
       
	} // End AntiLambda

      } // End pt selection
      // Xi-
      /*
	if(lPdgcodeCurrentPart == kXiMinus || lPdgcodeCurrentPart == 3322){

	if( isNaturalPart )
	fAssocPartMC->Add( new AliMiniParticle(centrality, zv, iTrkMC, lPtCurrentPart, lPhiCurrentPart, lEtaCurrentPart, 0, 0, 6) );

	} //End Xi
      */

    } // End loop over MC
    
    // -----------------------------------------
    // ---------- MC Correlations --------------
    // -----------------------------------------
    
    Float_t triggerMCPt   = -1000.;
    Float_t triggerMCPhi  = -1000.;
    Float_t triggerMCEta  = -1000.;
    
    Float_t dPhiMC = -100.;
    Float_t dEtaMC = -100.;
 
    for(Int_t ii=0; ii<(fTriggerPartMC->GetEntriesFast()); ii++){
      AliMiniParticle* trigMC = (AliMiniParticle*) fTriggerPartMC->At(ii);
      
      triggerMCPt  = trigMC->Pt();
      triggerMCPhi = trigMC->Phi();
      triggerMCEta = trigMC->Eta();

      fTriggerMCPtCent->Fill(triggerMCPt,centrality);
      
      for(Int_t jj=0; jj<(fAssocPartMC->GetEntriesFast()); jj++){
	
	AliMiniParticle* assocMC = (AliMiniParticle*) fAssocPartMC->At(jj);
	if(assocMC->Pt()>triggerMCPt) continue;
	
	dPhiMC = dPHI(triggerMCPhi,assocMC->Phi());
	dEtaMC = triggerMCEta - assocMC->Eta();
     
	// K0s, Lambdas and AntiLambdas (h-V0 correlations)
	if( (triggerMCPt<fTrigPtMax) && ( (assocMC->WhichCandidate()==3) || (assocMC->WhichCandidate()==4) || (assocMC->WhichCandidate()==5) ) )
	  for(Int_t k=0;k<kN1;k++)   // Pt bin
	    if( (assocMC->Pt()>=kPtBinV0[k]) && (assocMC->Pt()<kPtBinV0[k+1]) ){	      
	      if(assocMC->WhichCandidate()==3)
		fK0sdPhidEtaMC[curCentBin*kN1+k]->Fill(dPhiMC,dEtaMC,zv);    
	      if(assocMC->WhichCandidate()==4)
		fLambdadPhidEtaMC[curCentBin*kN1+k]->Fill(dPhiMC,dEtaMC,zv);
	      if(assocMC->WhichCandidate()==5)
		fAntiLambdadPhidEtaMC[curCentBin*kN1+k]->Fill(dPhiMC,dEtaMC,zv);
	    } // End pt bin

      } // End loop over trigger particles

    } // End loop over trigger particles

  } // End MC condition
  
  // *************************************************
  // V0 loop - AOD
  fAssocParticles = new TObjArray(); 
  fAssocParticles->SetOwner(kTRUE);
  if(NtrigPerEvnt>0)
    V0Loop(kReconstruction,kTRUE,-1,-1);
  else 
    V0Loop(kReconstruction,kFALSE,-1,-1); 
 
  //-------------------------------------------------------------
  // Correlations
  //-------------------------------------------------------------

  // variables for correlations:
  Float_t ptTrig=0., pxTrig=0., pyTrig=0.;
  Float_t massK0s=0., mK0s=0., sK0s=0.;
  Float_t massL=0.,   mL=0.,   sL=0.;
  Float_t massAL=0.; //,  mAL=0.,  sAL=0.;
  Float_t pt=-100., pxAssoc=-1000., pyAssoc=-1000.;
  Float_t lPhi=0., lEta=0.;
  Float_t lAlphaV0=0., lPtArmV0=0, dcaPos=0., dcaNeg=0.;
  Float_t dx=-100., dy=-100., lt=-100., res=-100.;
  Float_t dlK=-100., dlL=-100.;
  Float_t dPhi=-100., dEta=-100., radio=-100.;
  Double_t xDCA[2], cov[3];

  // variables for track splititing checks:
  Float_t  posdPhiS = -9999., posdEtaS = -9999., negdPhiS = -9999., negdEtaS = -9999.; 
  Float_t  fracTrigTPCSharedMap=-1., fracPosDaugTPCSharedMap =-1., fracNegDaugTPCSharedMap =-1.;
  //Bool_t   trigTPCMapOk=kTRUE, posDaugTPCMapOk=kTRUE, negDaugTPCMapOk=kTRUE;  
  Float_t  RdPhiStarMaxPosDaug=-1., RdPhiStarMaxNegDaug=-1., den=1.;
  Double_t trigCov[21], posDaugCov[21], negDaugCov[21];
  Double_t trigPos[6], posDaugPos[6], negDaugPos[6];
  Double_t trigXYZ[3], posDaugXYZ[3], negDaugXYZ[3];
  Double_t devPosDaugTrig[9], devNegDaugTrig[9], splitCont[9],  splitCont2[14];
  Int_t    sameSignPosDaug = -1, sameSignNegDaug = -1;
  Float_t  sameLabelPosDaug = 0., sameLabelNegDaug = 0.;
  Int_t    tlab, nlab, plab;
  Double_t resdEtsSdPhiS[6]; 

  // --------------------------------
  // weight to be used for the correlations due to the steps presenteed in the centrality distribution only for 2011 Pb-Pb data;
  Double_t weight = 1.;
  if( fCollision.Contains("PbPb2011") ){
    if( centrality >= 9.0 && centrality < 10.0 ) weight = 1.0675;
    else if( centrality >= 10.0 && centrality < 11.0 ) weight = 0.39188;
    else if( centrality >= 11.0 && centrality < 12.0 ) weight = 0.68262;
    else weight = 1.;
  }

  // --------------------------------
  // h-V0 correlations
  for (Int_t i=0; i<(fTriggerParticles->GetEntriesFast()); i++){
    AliMiniParticle* trig = (AliMiniParticle*) fTriggerParticles->At(i);
    if( trig->WhichCandidate() == 0 ) continue;

    const AliAODTrack *tTrig = (AliAODTrack*)fAOD->GetTrack(trig->ID());
    ptTrig = tTrig->Pt();  pxTrig = tTrig->Px();  pyTrig = tTrig->Py(); 
 
    Bool_t proptodca = ((AliAODTrack*)fAOD->GetTrack(trig->ID()))->PropagateToDCA(vtx,bSign,100.0,xDCA,cov);
    xDCA[0] = TMath::Abs(xDCA[0]);   xDCA[1] = TMath::Abs(xDCA[1]);

    fTriggerDCA->Fill(xDCA[0],1.);
    fTriggerDCA->Fill(xDCA[1],2.);
    //Printf(" %lf    %lf",xDCA[0],xDCA[1]);

    // ---------------- Fraction of TPC Shared Cluster: 
    fracTrigTPCSharedMap = GetFractionTPCSharedCls(tTrig);
    fTrigFracShTPCcls->Fill(ptTrig,fracTrigTPCSharedMap);

    for(Int_t j=0; j<fAssocParticles->GetEntriesFast(); j++){
      AliMiniParticle* trackAssocME = (AliMiniParticle*) (fAssocParticles->At(j));
      AliAODv0 *tAssoc=fAOD->GetV0(trackAssocME->ID());
      const AliAODTrack *ntrack=(AliAODTrack *)tAssoc->GetDaughter(1);
      const AliAODTrack *ptrack=(AliAODTrack *)tAssoc->GetDaughter(0);

      if( SameTrack(tTrig,ntrack) || SameTrack(tTrig,ptrack) )
	continue;

      if( ptTrig < trackAssocME->Pt() ) continue;

      lPhi = trackAssocME->Phi();
      lEta = trackAssocME->Eta();

      // Correlation in deltaPhi & deltaEta
      dPhi = dPHI(trig->Phi(),lPhi);
      dEta = trig->Eta() - lEta;
      radio    = TMath::Sqrt(dPhi*dPhi + dEta*dEta);
     
      // Armenteros variables: 
      lAlphaV0      =  tAssoc->AlphaV0();
      lPtArmV0      =  tAssoc->PtArmV0();

      // 2D momentum
      pt = trackAssocME->Pt(); pxAssoc = tAssoc->Px(); pyAssoc = tAssoc->Py(); 
      // Decay vertex
      Double_t xyz[3]; tAssoc->GetSecondaryVtx(xyz);
      dx=xyz[0]-xv; dy=xyz[1]-yv; //dz=xyz[2]-zv;
      // Decay length: 2D 
      lt=TMath::Sqrt(dx*dx + dy*dy); 
      // Spatial resolution trigger-V0 point decay
      res = SpatialResolution(pxTrig,pyTrig,pxAssoc,pyAssoc,lt);
      // Ctau
      dlK = 0.4977*lt/pt;
      dlL = 1.1157*lt/pt; 

      Int_t binPtv0 = PtBin( pt );
      if(binPtv0==-1) continue;

      Int_t lMCAssocNegDaug = trackAssocME->NegDaugMCLabel();
      Int_t lMCAssocPosDaug = trackAssocME->PosDaugMCLabel();

      // ----------------------------------------------------------------------------

      // -----------------------------------------------------------------
      //   ****************** Track splitting check ******************
      // -----------------------------------------------------------------

      sameLabelPosDaug = 0.; sameLabelNegDaug = 0.;
      sameSignPosDaug = -1; sameSignNegDaug = -1;
      RdPhiStarMaxPosDaug=-1.; RdPhiStarMaxNegDaug=-1.;
      //trigTPCMapOk=kTRUE; posDaugTPCMapOk=kTRUE; negDaugTPCMapOk=kTRUE;
      fracPosDaugTPCSharedMap=0; fracNegDaugTPCSharedMap=0;

      // ---------------- Fraction of TPC Shared Cluster 
      fracPosDaugTPCSharedMap = GetFractionTPCSharedCls(ptrack);
      fracNegDaugTPCSharedMap = GetFractionTPCSharedCls(ntrack);

      // =========== Classical methods for track-splitting  ============= //
      if( TMath::Abs(dPhi)<0.1 && TMath::Abs(dEta)<0.1 ){   
	
	// --------- Check sign of the trigger and daughter track:
	if(tTrig->Charge()==1) { sameSignPosDaug = 1; sameSignNegDaug = 0; }
	else { sameSignPosDaug = 0; sameSignNegDaug = 1; }

	// -------- Shifting charged tracks to the primary vertex.
	// -------- See HBT anlayses: 

	// Trigger particle: 
	SetSftPosR125(tTrig,bSign,priVtx,"Trigger");

	// Positive daughter: calculating delta(phi)* and delta(eta)*
	SetSftPosR125(ptrack,bSign,priVtx,"Daughter");
	posdPhiS = dPhiSAtR125();
	posdEtaS = dEtaS();

	// Negative daughter: calculating delta(phi)* and delta(eta)*
	SetSftPosR125(ntrack,bSign,priVtx,"Daughter");
	negdPhiS = dPhiSAtR125();
	negdEtaS = dEtaS();

	// ------ Get position:
	tTrig->GetXYZ(trigXYZ);
	ptrack->GetXYZ(posDaugXYZ);
	ntrack->GetXYZ(negDaugXYZ);

	// ------ Covaraince matrix for the tracks:
	tTrig->GetCovarianceXYZPxPyPz(trigCov);
	ptrack->GetCovarianceXYZPxPyPz(posDaugCov);
	ntrack->GetCovarianceXYZPxPyPz(negDaugCov);

	// ------- position and momentum:
	// trigger particle
	trigPos[0] = trigXYZ[0];	trigPos[1] = trigXYZ[1];	trigPos[2] = trigXYZ[2];
	trigPos[3] = tTrig->Px();	trigPos[4] = tTrig->Py();	trigPos[5] = tTrig->Pz();

	// positive daughter
	posDaugPos[0] = posDaugXYZ[0];	posDaugPos[1] = posDaugXYZ[1];	posDaugPos[2] = posDaugXYZ[2];
	posDaugPos[3] = ptrack->Px();	posDaugPos[4] = ptrack->Py();	posDaugPos[5] = ptrack->Pz();

	// negative daughter
	negDaugPos[0] = negDaugXYZ[0];	negDaugPos[1] = negDaugXYZ[1];	negDaugPos[2] = negDaugXYZ[2];
	negDaugPos[3] = ntrack->Px();	negDaugPos[4] = ntrack->Py();	negDaugPos[5] = ntrack->Pz();

	// ------- deviation between the two tracks:
	// positive daughter
	for(Int_t ll=0;ll<6;ll++){
	  den = trigCov[ll*(ll+1)/2+ll] +  posDaugCov[ll*(ll+1)/2+ll] ;
	  devPosDaugTrig[ll] = 0.;
	  
	  if(den>0)  devPosDaugTrig[ll] = TMath::Power( trigPos[ll] - posDaugPos[ll] ,2) / den;
	  
	  if(ll<3) devPosDaugTrig[6] +=  devPosDaugTrig[ll];  // sum in X,Y,Z
	  if(ll>2) devPosDaugTrig[7] +=  devPosDaugTrig[ll];  // sum in momemtum
	  devPosDaugTrig[8] +=  devPosDaugTrig[ll];           // sum in all variables
	}

	// negative daughter
	for(Int_t ll=0;ll<6;ll++){
	  den = trigCov[ll*(ll+1)/2+ll]  +  negDaugCov[ll*(ll+1)/2+ll] ;
	  devNegDaugTrig[ll] = 0;

	  if(den>0)  devNegDaugTrig[ll] = TMath::Power( trigPos[ll] - negDaugPos[ll] ,2) / den;
	  
	  if(ll<3) devNegDaugTrig[6] +=  devNegDaugTrig[ll];  // sum in X,Y,Z
	  if(ll>2) devNegDaugTrig[7] +=  devNegDaugTrig[ll];  // sum in momemtum
	  devNegDaugTrig[8] +=  devNegDaugTrig[ll];           // sum in all variables

	}


	// ---------------- Monte Carlo check for track-splitting 
	if(fIsMC){
	     
	  TList *lst = fAOD->GetList();
	  stack = (TClonesArray*)lst->FindObject(AliAODMCParticle::StdBranchName());
	  if (!stack) {
	    Printf("ERROR: stack not available");
	    return;
	  }

	  sameLabelPosDaug = 1.*SameLabel(tTrig,ptrack);
	  sameLabelNegDaug = 1.*SameLabel(tTrig,ntrack);

	  // Resolution of delta(phi)* and delta(eta)*
	  tlab = TMath::Abs(tTrig->GetLabel());
	  plab = TMath::Abs(ptrack->GetLabel());
	  nlab = TMath::Abs(ntrack->GetLabel());

	  AliAODMCParticle *tPart=(AliAODMCParticle*)stack->UncheckedAt(tlab);
	  AliAODMCParticle *pPart=(AliAODMCParticle*)stack->UncheckedAt(plab);
	  AliAODMCParticle *nPart=(AliAODMCParticle*)stack->UncheckedAt(nlab);

	  resdEtsSdPhiS[0] = pt;

	  //positive daughter
	  resdEtsSdPhiS[2] = ptrack->Pt();
	  resdEtsSdPhiS[3] = (tPart->Phi() - pPart->Phi()) - posdPhiS;
	  resdEtsSdPhiS[4] = (tPart->Eta() - pPart->Eta()) - posdEtaS;
	  resdEtsSdPhiS[5] = sameSignPosDaug;

	  if( trackAssocME->WhichCandidate() == 3 ){
	    resdEtsSdPhiS[1] = massK0s;
	    fK0sPosMCResdEtaSdPhiS[curCentBin]->Fill(resdEtsSdPhiS);
	  }
	  if( trackAssocME->WhichCandidate() == 4 ){
	    resdEtsSdPhiS[1] = massL;
	    fLambdaPosMCResdEtaSdPhiS[curCentBin]->Fill(resdEtsSdPhiS);
	  }
	  if( trackAssocME->WhichCandidate() == 5 ){
	    resdEtsSdPhiS[1] = massAL;
	    fAntiLambdaPosMCResdEtaSdPhiS[curCentBin]->Fill(resdEtsSdPhiS);
	  }

	  // negative daughter
	  resdEtsSdPhiS[2] = ntrack->Pt();
	  resdEtsSdPhiS[3] = (tPart->Phi() - nPart->Phi()) - negdPhiS;
	  resdEtsSdPhiS[4] = (tPart->Eta() - nPart->Eta()) - negdEtaS;
	  resdEtsSdPhiS[5] = sameSignNegDaug;

	  if( trackAssocME->WhichCandidate() == 3 ){
	    resdEtsSdPhiS[1] = massK0s;
	    fK0sNegMCResdEtaSdPhiS[curCentBin]->Fill(resdEtsSdPhiS);
	  }
	  if( trackAssocME->WhichCandidate() == 4 ){
	    resdEtsSdPhiS[1] = massL;
	    fLambdaNegMCResdEtaSdPhiS[curCentBin]->Fill(resdEtsSdPhiS);
	  }
	  if( trackAssocME->WhichCandidate() == 5 ){
	    resdEtsSdPhiS[1] = massAL;
	    fAntiLambdaNegMCResdEtaSdPhiS[curCentBin]->Fill(resdEtsSdPhiS);
	    }
	   
	}

	// ================  Alternative methods for track-splitting  ==================
	if(TMath::Abs(dPhi)<0.02 && TMath::Abs(dEta)<0.02){
	  
	  // --------- Calculate TPCRadius*Delta(phi)Star_Max distance:
	  RdPhiStarMaxPosDaug = TwoTrackEfficiencyCut( tTrig->Phi(), tTrig->Eta(), tTrig->Pt(), tTrig->Charge(), ptrack->Phi(), ptrack->Eta(), ptrack->Pt(), 1, bSign);
	  RdPhiStarMaxNegDaug = TwoTrackEfficiencyCut( tTrig->Phi(), tTrig->Eta(), tTrig->Pt(), tTrig->Charge(), ntrack->Phi(), ntrack->Eta(), ntrack->Pt(), -1, bSign);

	  // -------- Comparison between trigger and daughter tracks:
	  // -------- Filling deviation of matrix elements
	  splitCont[0] = pt; splitCont[5] = fracTrigTPCSharedMap; 

	  // ---------------------------
	  // -------- Positive daughter:
	  splitCont[2] = ptrack->Pt();  splitCont[3] = sameSignPosDaug; 
	  splitCont[4] = RdPhiStarMaxPosDaug;   splitCont[6] = fracPosDaugTPCSharedMap; 
	    
	  // ----K0s
	  if( trackAssocME->WhichCandidate() == 3 ){
	    splitCont[1] = massK0s; 
	    for(Int_t ll=0; ll<=8; ll++){
	      splitCont[7] = devPosDaugTrig[ll]; splitCont[8] = ll; 
	      fK0sPosDaugSplCheckCovMat[curCentBin]->Fill(splitCont);
	    } 

	  }
	  // ----Lambda
	  if( trackAssocME->WhichCandidate() == 4 ){
	    splitCont[1] = massL; 
	    for(Int_t ll=0; ll<=8; ll++){
	      splitCont[7] = devPosDaugTrig[ll]; splitCont[8] = ll; 
	      fLambdaPosDaugSplCheckCovMat[curCentBin]->Fill(splitCont);
	    } 

	  }
	  // ----AntiLambda
	  if( trackAssocME->WhichCandidate() == 5 ){
	    splitCont[1] = massAL; 
	    for(Int_t ll=0; ll<=8; ll++){
	      splitCont[7] = devPosDaugTrig[ll]; splitCont[8] = ll; 
	      fAntiLambdaPosDaugSplCheckCovMat[curCentBin]->Fill(splitCont);
	    }
	      
	  }
	  // End: Positive daughter

	  // ---------------------------
	  // -------- Negative daughter:
	  splitCont[2] = ntrack->Pt(); splitCont[3] = sameSignNegDaug; 
	  splitCont[4] = RdPhiStarMaxNegDaug;   splitCont[6] = fracNegDaugTPCSharedMap; 
	
	  // ----K0s
	  if( trackAssocME->WhichCandidate() == 3 ){
	    splitCont[1] = massK0s;  
	    for(Int_t ll=0; ll<=8; ll++){
	      splitCont[7] = devNegDaugTrig[ll]; splitCont[8] = ll; 
	      fK0sNegDaugSplCheckCovMat[curCentBin]->Fill(splitCont);
	    }

	  }
	  // ----Lambda
	  if( trackAssocME->WhichCandidate() == 4 ){
	    splitCont[1] = massL; 
	    for(Int_t ll=0; ll<=8; ll++){
	      splitCont[7] = devNegDaugTrig[ll]; splitCont[8] = ll; 
	      fLambdaNegDaugSplCheckCovMat[curCentBin]->Fill(splitCont);
	    }
	      
	  }
	  // ----AntiLambda
	  if( trackAssocME->WhichCandidate() == 5 ){
	    splitCont[1] = massAL; 
	    for(Int_t ll=0; ll<=8; ll++){
	      splitCont[7] = devNegDaugTrig[ll]; splitCont[8] = ll; 
	      fAntiLambdaNegDaugSplCheckCovMat[curCentBin]->Fill(splitCont);
	    }

	  }
	  // End: Negative daughter
   
	} // end selection in |delta(eta)| < 0.02, |delta(phi)| < 0.02


	// ================  FILLING THnSparse:  Classical track-splitting method: d(phi)* and d(eta)*
	splitCont2[0] = pt;	    splitCont2[6] = fracTrigTPCSharedMap; 
	splitCont2[11] = xDCA[0];   splitCont2[12] = xDCA[1];
	// --------------------------
	// -------- Positive daughter:
	splitCont2[2] = ptrack->Pt();  splitCont2[3] = sameSignPosDaug;  splitCont2[4] = posdPhiS;  splitCont2[5] = posdEtaS; 
	splitCont2[7] = fracPosDaugTPCSharedMap;   splitCont2[8] = fracTrigTPCSharedMap - fracPosDaugTPCSharedMap;
	splitCont2[9] = devPosDaugTrig[7];  splitCont2[10] = tAssoc->DcaPosToPrimVertex(); splitCont2[13] = sameLabelPosDaug; 
 
	// ---- K0s
	if( trackAssocME->WhichCandidate() == 3 ){
	  splitCont2[1] = massK0s;  
	  // Positive daughter 
	  fK0sPosDaugdPhiSdEtaS[curCentBin]->Fill(splitCont2);	  
       	}
	// ---- Lambda
 	if( trackAssocME->WhichCandidate() == 4 ){
	  splitCont2[1] = massL;  
	  // Positive daughter 
	  fLambdaPosDaugdPhiSdEtaS[curCentBin]->Fill(splitCont2);	  
	}
	// ---- AntiLambda
	if( trackAssocME->WhichCandidate() == 5 ){
	  splitCont2[1] = massAL;  
	  // Positive daughter
	  fAntiLambdaPosDaugdPhiSdEtaS[curCentBin]->Fill(splitCont2);	  
	}
	
	// --------------------------
	// ------- Negative daughter:
	splitCont2[2] = ntrack->Pt();  splitCont2[3] = sameSignNegDaug;  splitCont2[4] = negdPhiS;  splitCont2[5] = negdEtaS; 
	splitCont2[7] = fracNegDaugTPCSharedMap;  splitCont2[8] = fracTrigTPCSharedMap - fracNegDaugTPCSharedMap;
	splitCont2[9] = devNegDaugTrig[7];  splitCont2[10] = tAssoc->DcaNegToPrimVertex();  splitCont2[13] = sameLabelNegDaug;  

	// ---- K0s
	if( trackAssocME->WhichCandidate() == 3 ){
	  splitCont2[1] = massK0s;  
      	  // Negative daughter
 	  fK0sNegDaugdPhiSdEtaS[curCentBin]->Fill(splitCont2);
	}
	// ---- Lambda
 	if( trackAssocME->WhichCandidate() == 4 ){
	  splitCont2[1] = massL;  	  
	  // Negative daughter	
	  fLambdaNegDaugdPhiSdEtaS[curCentBin]->Fill(splitCont2);
	}
	// ---- AntiLambda
	if( trackAssocME->WhichCandidate() == 5 ){
	  splitCont2[1] = massAL;  	  
	  // Negative daughter 
	  fAntiLambdaNegDaugdPhiSdEtaS[curCentBin]->Fill(splitCont2);
	}
      
      } // end selection in |delta(eta)| < 0.1, |delta(phi)| < 0.1

      // ----------------------------------------------------------------
      // Reject the 'fake' correlation due to the TPC shared clusters
      // between trigger particle and one of the daughter tracks 
      //    The rejection will affect more the correlations:
      //         - Trigger track - Positive track (from Lambda with pt above 3 GeV/c)
      //         - Trigger track - Negative track (from AntiLambda with pt above 3 GeV/c)
      /* if( fracTrigTPCSharedMap>0.5 && 
	  ( ( sameSignPosDaug==1 && TMath::Abs(fracTrigTPCSharedMap - fracPosDaugTPCSharedMap) < fDiffTrigDaugFracTPCSharedCls ) ||
	  ( sameSignNegDaug==1 && TMath::Abs(fracTrigTPCSharedMap - fracNegDaugTPCSharedMap) < fDiffTrigDaugFracTPCSharedCls ) ) )*/

      if( (fracTrigTPCSharedMap > fFracTPCcls) || (fracPosDaugTPCSharedMap > fFracTPCcls) || (fracNegDaugTPCSharedMap > fFracTPCcls) )
	continue;

      // ----------------------------------------------------------------------------
        
      // *******************
      //   K0s selection
      // *******************
      if( trackAssocME->WhichCandidate() == 3 ){
              
	massK0s = tAssoc->MassK0Short();
	mK0s = TDatabasePDG::Instance()->GetParticle(kK0Short)->Mass();
	if( fCollision.Contains("PbPb2010") )
	  sK0s = kCteK0s2010[curCentBin] + kLinearK0s2010[curCentBin]*pt;
	else if( fCollision.Contains("PbPb2011") ) 
	  sK0s = kCteK0s2011[curCentBin] + kLinearK0s2011[curCentBin]*pt;
	
	// ==== Correlations K0s invariant mass peak ==== //
	// +++++++++++ Pt bin & centrality
	fK0sdPhidEtaPtL[curCentBin*kN1*kNVtxZ + binPtv0*kNVtxZ + curVtxBin]->Fill(dPhi,dEta,massK0s,weight);

	fK0sPosDaugFracShTPCclsTrig->Fill(massK0s,pt,fracPosDaugTPCSharedMap);
	fK0sNegDaugFracShTPCclsTrig->Fill(massK0s,pt,fracNegDaugTPCSharedMap);

	// ==== Correlations K0s invariant mass peak ==== //
	if (TMath::Abs(mK0s-massK0s) < 3*sK0s) {

	  if(radio<0.02){
	    fK0sSpatialRes->Fill(dPhi,res,lt);
	  }
	  if(radio < 0.4){
	    fHistArmPodBckg->Fill(lAlphaV0,lPtArmV0,0);
	    fK0sDCADaugToPrimVtx->Fill(dcaPos,dcaNeg,ptTrig);	    
	    RecCascade(tTrig,ntrack,ptrack,"K0s");
	    RecCascade(tTrig,ptrack,ntrack,"K0s");	
	  }
			
	}
	// ==== Correlations K0s background ==== //
	if( TMath::Abs(mK0s-massK0s + 6.5*sK0s) < 1.5*sK0s ||
	    TMath::Abs(mK0s-massK0s - 6.5*sK0s) < 1.5*sK0s  ) {
	  
	  // Only fills the histograms when it is a triggered event
	  if(j==0){
	    fHistArmenterosPodolanski->Fill(lAlphaV0,lPtArmV0,1);
	    // MC Association of daughter particles 
	    fK0sBckgDCANegDaugToPrimVtx->Fill(lMCAssocNegDaug,dcaNeg);
	    fK0sBckgDCAPosDaugToPrimVtx->Fill(lMCAssocPosDaug,dcaPos);
	  }
	  
	  // Pt bin & centrality
	  //fK0sdPhidEtaPtLBckg[curCentBin*kN1+binPtv0]->Fill(dPhi,dEta,zv);
	    
	  if(radio < 0.4){ // Under the correlation peak
	    fHistArmPodBckg->Fill(lAlphaV0,lPtArmV0,1);
	    fK0sBckgDecLength->Fill(dlK,ptTrig);
	    fK0sBckgDCADaugToPrimVtx->Fill(dcaPos,dcaNeg,ptTrig);
	    fK0sBckgEtaPhi->Fill(lPhi,lEta);
	    fK0sBckgPhiRadio->Fill(lPhi,lt);

	    //RecCascade(trkTrig,ntrack,ptrack,"K0s");
	    //RecCascade(trkTrig,ptrack,ntrack,"K0s");

	  }// End selection in the correlation peak
		
	} // End background selection
	
      } // End K0s selection

      // *******************
      // Lambda selection
      // *******************
      if( trackAssocME->WhichCandidate() == 4 ){
	massL = tAssoc->MassLambda();
	mL = TDatabasePDG::Instance()->GetParticle(kLambda0)->Mass();
	if( fCollision.Contains("PbPb2010") )
	  sL = kCteLambda2010[curCentBin] + kLinearLambda2010[curCentBin]*pt;
	else if( fCollision.Contains("PbPb2011") ) 
	  sL = kCteLambda2011[curCentBin] + kLinearLambda2011[curCentBin]*pt;

	// ==== Correlations Lambda invariant mass peak ==== //
        // +++++++++++ Pt bin & centrality
	fLambdadPhidEtaPtL[curCentBin*kN1*kNVtxZ + binPtv0*kNVtxZ + curVtxBin]->Fill(dPhi,dEta,massL,weight);
	
	fLambdaPosDaugFracShTPCclsTrig->Fill(massL,pt,fracPosDaugTPCSharedMap);
	fLambdaNegDaugFracShTPCclsTrig->Fill(massL,pt,fracNegDaugTPCSharedMap);

	// ==== Correlations Lambda invariant mass peak ==== //
	if (TMath::Abs(mL-massL) < 3*sL) {			  

	  if(radio<0.02)
	    fLambdaSpatialRes->Fill(dPhi,res,lt);
	  if(radio < 0.4){
	    fHistArmPodBckg->Fill(lAlphaV0,lPtArmV0,2);
	    fLambdaDCADaugToPrimVtx->Fill(dcaPos,dcaNeg,ptTrig);
	    RecCascade(tTrig,ntrack,ptrack,"Lambda");
	    RecCascade(tTrig,ptrack,ntrack,"Lambda");
	  }
	    
	} // End mass peak selection
	// ==== Correlations Lambda background ==== //
	if( TMath::Abs(mL-massL + 6.5*sL) < 1.5*sL ||
	    TMath::Abs(mL-massL - 6.5*sL) < 1.5*sL ) {

	  // Only fills the histograms when it is a triggered event
	  if(j==0){
	    fHistArmenterosPodolanski->Fill(lAlphaV0,lPtArmV0,3);
	    // MC Association of daughter particles 
	    fLambdaBckgDCANegDaugToPrimVtx->Fill(lMCAssocNegDaug,dcaNeg);
	    fLambdaBckgDCAPosDaugToPrimVtx->Fill(lMCAssocPosDaug,dcaPos);
	  }

	  // Pt bin & centrality
	  //fLambdadPhidEtaPtLBckg[curCentBin*kN1+binPtv0]->Fill(dPhi,dEta,zv);
	
	  if(radio < 0.4){ // Under the peak
	    fHistArmPodBckg->Fill(lAlphaV0,lPtArmV0,3);
	    fLambdaBckgDecLength->Fill(dlL,ptTrig);
	    fLambdaBckgDCADaugToPrimVtx->Fill(dcaPos,dcaNeg,ptTrig);
	    fLambdaBckgEtaPhi->Fill(lPhi,lEta);
	    fLambdaBckgPhiRadio->Fill(lPhi,lt);
		  
	    //RecCascade(trkTrig,ntrack,ptrack,"Lambda");
	    //RecCascade(trkTrig,ptrack,ntrack,"Lambda");

	  }// End selection in the correlation peak
		
	} // End background selection
	
      }// End Lambda selection
       
       // *******************
      // AntiLambda selection
      // *******************
      if( trackAssocME->WhichCandidate() == 5 ){
	massAL = tAssoc->MassAntiLambda();
	mL = TDatabasePDG::Instance()->GetParticle(kLambda0)->Mass();
	if( fCollision.Contains("PbPb2010") )
	  sL = kCteAntiLambda2010[curCentBin] + kLinearAntiLambda2010[curCentBin]*pt;
	else if( fCollision.Contains("PbPb2011") ) 
	  sL = kCteAntiLambda2011[curCentBin] + kLinearAntiLambda2011[curCentBin]*pt;
	
	// ==== Correlations Lambda invariant mass peak ==== //
        // +++++++++++ Pt bin & centrality
	fAntiLambdadPhidEtaPtL[curCentBin*kN1*kNVtxZ + binPtv0*kNVtxZ + curVtxBin]->Fill(dPhi,dEta,massAL,weight);

	fAntiLambdaPosDaugFracShTPCclsTrig->Fill(massAL,pt,fracPosDaugTPCSharedMap);
	fAntiLambdaNegDaugFracShTPCclsTrig->Fill(massAL,pt,fracNegDaugTPCSharedMap);

	// ==== Correlations AntiLambda invariant mass peak ==== //
	if (TMath::Abs(mL-massAL) < 3*sL) {

	  if(radio<0.1)
	    fAntiLambdaSpatialRes->Fill(dPhi,res,lt);	      
	  if(radio < 0.4){
	    fHistArmPodBckg->Fill(lAlphaV0,lPtArmV0,4);
	    fAntiLambdaDCADaugToPrimVtx->Fill(dcaPos,dcaNeg,ptTrig);
	    RecCascade(tTrig,ntrack,ptrack,"AntiLambda");
	    RecCascade(tTrig,ptrack,ntrack,"AntiLambda");
	  }
	      
	} // End AntiLambda mass peak
	// ==== Correlations AntiLambda background ==== //
	if( (TMath::Abs(mL-massAL + 6.5*sL) < 1.5*sL) ||
	    (TMath::Abs(mL-massAL - 6.5*sL) < 1.5*sL) ){

	  // ----------------------------------------------

	  // Only fills the histograms when it is a triggered event
	  if(j==0){
	    fHistArmenterosPodolanski->Fill(lAlphaV0,lPtArmV0,5);
	    // MC Association of daughter particles 
	    fAntiLambdaBckgDCANegDaugToPrimVtx->Fill(lMCAssocNegDaug,dcaNeg);
	    fAntiLambdaBckgDCAPosDaugToPrimVtx->Fill(lMCAssocPosDaug,dcaPos);
	  }
	    
	  // Pt bin & centrality
	  //fAntiLambdadPhidEtaPtLBckg[curCentBin*kN1+binPtv0]->Fill(dPhi,dEta,zv);
	
	  if(radio < 0.4){ // Under the peak
	    fHistArmPodBckg->Fill(lAlphaV0,lPtArmV0,5);
	    fAntiLambdaBckgDecLength->Fill(dlL,ptTrig);
	    fAntiLambdaBckgDCADaugToPrimVtx->Fill(dcaPos,dcaNeg,ptTrig);
	    fAntiLambdaBckgEtaPhi->Fill(lPhi,lEta);
	    fAntiLambdaBckgPhiRadio->Fill(lPhi,lt);
		  
	    //RecCascade(trkTrig,ntrack,ptrack,"AntiLambda");
	    //RecCascade(trkTrig,ptrack,ntrack,"AntiLambda");

	  }// End selection in the correlation peak
		
	}// End AntiLambda background

      } // End AntiLambda selection

    } // End loop over associated particles
   

    // Filling information of the trigger particle
    // after the rejection in the cut of shared TPC cls
    fTriggerEtaPhi->Fill(trig->Phi(),trig->Eta());
    fTriggerPtCent->Fill(trig->Pt(),centrality,zv);

  } // End loop over trigger particles
 

  //-------------------------------------------------------------
  // Mixing
  //-------------------------------------------------------------
  
  Double_t phiTrigME=0, etaTrigME=0, phiAssocME=0, etaAssocME=0;
  Double_t deltaPhi=0, deltaEta=0;

  TList *evMixList = fMEList[curCentBin*kNVtxZ+curVtxBin];
  Int_t nMixed = evMixList->GetSize(); 
 
  if( nMixed>0 && fAssocParticles->GetEntriesFast() >= 0 ){
    
    for(Int_t ii=0; ii<nMixed; ii++){     
      
      AliMiniParticle* trackTriggerME = (AliMiniParticle*) (evMixList->At(ii));
      phiTrigME = trackTriggerME->Phi();
      etaTrigME = trackTriggerME->Eta();

      // --- V0 associated particles
      for(Int_t j=0; j<fAssocParticles->GetEntriesFast(); j++){
	
	AliMiniParticle* trackAssocME = (AliMiniParticle*) (fAssocParticles->At(j));
	if( CentBin(trackTriggerME->Centrality()) != CentBin(trackAssocME->Centrality()) ) continue;
	if( VtxBin(trackTriggerME->VtxZ()) != VtxBin(trackAssocME->VtxZ()) ) continue;
	if( trackAssocME->WhichCandidate() ==  2 ) continue;

	AliAODv0 *tAssoc=fAOD->GetV0(trackAssocME->ID());
	const AliAODTrack *ntrack=(AliAODTrack *)tAssoc->GetDaughter(1);
	const AliAODTrack *ptrack=(AliAODTrack *)tAssoc->GetDaughter(0);

	// Fraction of TPC Shared Cluster 
	fracPosDaugTPCSharedMap=0; fracNegDaugTPCSharedMap=0;
	fracPosDaugTPCSharedMap = GetFractionTPCSharedCls(ptrack);
	fracNegDaugTPCSharedMap = GetFractionTPCSharedCls(ntrack);

	if( (fracPosDaugTPCSharedMap > fFracTPCcls) || (fracNegDaugTPCSharedMap > fFracTPCcls) )
	  continue;

	pt = tAssoc->Pt();

	massK0s = tAssoc->MassK0Short();
	massL   = tAssoc->MassLambda();
	massAL  = tAssoc->MassAntiLambda();

	/*
	Bool_t IsSelected = kFALSE;
	// K0s
	if( trackAssocME->WhichCandidate() == 3 ){
	  massK0s = tAssoc->MassK0Short();
	  mK0s = TDatabasePDG::Instance()->GetParticle(kK0Short)->Mass();
	  if( fCollision.Contains("PbPb2010") )
	    sK0s = kCteK0s2010[curCentBin] + kLinearK0s2010[curCentBin]*pt;
	  else if( fCollision.Contains("PbPb2011") ) 
	    sK0s = kCteK0s2011[curCentBin] + kLinearK0s2011[curCentBin]*pt;
	  
	  if (TMath::Abs(mK0s-massK0s) < 3*sK0s) IsSelected = kTRUE;
	}
	// Lambda
	if( trackAssocME->WhichCandidate() == 4 ){
	  massL = tAssoc->MassLambda();
	  mL = TDatabasePDG::Instance()->GetParticle(kLambda0)->Mass();	  
	  if( fCollision.Contains("PbPb2010") )
	    sL = kCteLambda2010[curCentBin] + kLinearLambda2010[curCentBin]*pt;
	  else if( fCollision.Contains("PbPb2011") ) 
	    sL = kCteLambda2011[curCentBin] + kLinearLambda2011[curCentBin]*pt;

	  if (TMath::Abs(mL-massL) < 3*sL) IsSelected = kTRUE;
	}
	// AntiLambda
	if( trackAssocME->WhichCandidate() == 5 ){
	  massAL = tAssoc->MassAntiLambda();
	  mL = TDatabasePDG::Instance()->GetParticle(kLambda0)->Mass();
	  if( fCollision.Contains("PbPb2010") )
	    sL = kCteAntiLambda2010[curCentBin] + kLinearAntiLambda2010[curCentBin]*pt;
	  else if( fCollision.Contains("PbPb2011") ) 
	    sL = kCteAntiLambda2011[curCentBin] + kLinearAntiLambda2011[curCentBin]*pt;
	  
	  if (TMath::Abs(mL-massAL) < 3*sL) IsSelected = kTRUE;
	}

	if(!IsSelected) continue;*/

	phiAssocME = trackAssocME->Phi();
	etaAssocME = trackAssocME->Eta();
	 
	deltaPhi = dPHI(phiTrigME,phiAssocME);
	deltaEta = etaTrigME - etaAssocME;

	Int_t binPtv0 = PtBin( trackAssocME->Pt() );
	if(binPtv0==-1) continue;
    
	if( trackAssocME->WhichCandidate() == 3 ) {
	  fK0sdPhidEtaME[curCentBin*kN1*kNVtxZ + binPtv0*kNVtxZ + curVtxBin]->Fill(deltaPhi,deltaEta,massK0s);}
	else if( trackAssocME->WhichCandidate() == 4 )
	  fLambdadPhidEtaME[curCentBin*kN1*kNVtxZ + binPtv0*kNVtxZ + curVtxBin]->Fill(deltaPhi,deltaEta,massL);
	else if( trackAssocME->WhichCandidate() == 5 )
	  fAntiLambdadPhidEtaME[curCentBin*kN1*kNVtxZ + binPtv0*kNVtxZ + curVtxBin]->Fill(deltaPhi,deltaEta,massAL);
	             
      } // End loop over V0's
       
    }
    
  }
  
  //--------------------------------------------------------
  //Add the current event to the list of events for mixing
  //--------------------------------------------------------  
  
  //Add current  event to buffer and Remove redundant events 
  if(fTriggerParticles->GetEntriesFast()>=0){
    
    for(Int_t ii=0; ii<(fTriggerParticles->GetEntriesFast()); ii++){
      AliMiniParticle* trkTrig = (AliMiniParticle*) fTriggerParticles->At(ii);
      //cout << trkTrig->Pt() << "          " << ii << endl;
    
      // Fraction of TPC Shared Cluster 
      const AliAODTrack *tTrig = (AliAODTrack*)fAOD->GetTrack(trkTrig->ID());
      fracTrigTPCSharedMap = GetFractionTPCSharedCls(tTrig);
      if( (fracTrigTPCSharedMap > fFracTPCcls) ) continue;

      if(evMixList->GetSize() < nMaxEvMix)
	evMixList->AddFirst(trkTrig);
      /*
	  if(evMixList->GetSize() >= nMaxEvMix) {
	    AliMiniParticle *tmp = (AliMiniParticle*) (evMixList->Last()) ;
	    evMixList->RemoveLast();
	    delete tmp;
	  }
      */
      
    }// End loop over fTriggerParticles

  }// End adding trigger particles to buffers
  
}

//___________________________________________________________________________________________

void AliAnalysisTaskLambdaOverK0sJets::Terminate(Option_t *)
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.
  
  fOutput=(TList*)GetOutputData(1);
  fOutputME=(TList*)GetOutputData(2);
  fOutputQA=(TList*)GetOutputData(3);

  if (fOutput || fOutputME || fOutputQA) {

    if(fOutput)
      Printf("\n\t *** DONE: fOutput available *** \n");
    if(fOutputME)
      Printf("\n\t *** DONE: fOutputME available *** \n");
    if (fOutputQA)
      Printf("\n\t *** DONE: fOutputQA available *** \n");
  }
  if (!fOutput || !fOutputME || !fOutputQA) {

    if(!fOutput)
      Printf("\n\t *** ERROR: fOutput not available *** \n");
    if(!fOutputME) 
      Printf("\n\t *** ERROR: fOutputME available *** \n");
    if(!fOutputQA)
      Printf("\n\t *** ERROR: fOutputQA not available  *** \n");  
  }

  
  return;

}
