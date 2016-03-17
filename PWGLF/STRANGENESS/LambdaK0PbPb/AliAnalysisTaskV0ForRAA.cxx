/**************************************************************************
 * Authors : Simone Schuchmann                                            *
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

//-----------------------------------------------------------------
// AliAnalysisTaskV0ForRAA class
// This task is for analysing Lambda and K0s pt spectra in PbPb and
// pp as well as with MC. The flag for pp and MC  must be set
// accordingly, default is PbPb data.
// It works with ESD files only.
//-----------------------------------------------------------------


#define AliAnalysisTaskV0ForRAA_cxx


#include "AliAnalysisTaskV0ForRAA.h"

#include "Riostream.h"
//#include "THn.h"
#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
//#include "TH3.h"//xxx
#include "TLorentzVector.h"

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliESDInputHandler.h"
#include "AliMCEventHandler.h"

#include "AliESDVZERO.h"
#include "AliMultiplicity.h"
#include "AliCentrality.h"

#include "AliKFParticle.h"
#include "AliKFVertex.h"
#include "AliESDtrackCuts.h"
#include "AliESDpid.h"
#include "AliESDv0.h"
#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliGenEventHeader.h"


ClassImp(AliAnalysisTaskV0ForRAA)

//________________________________________________________________________
AliAnalysisTaskV0ForRAA::AliAnalysisTaskV0ForRAA()
:AliAnalysisTaskSE("default_AliAnalysisTaskV0ForRAA"),
  fESD(0),
  fMCev(0),
//other objects
  fESDpid(0),
  fESDTrackCuts(0),
  fESDTrackCutsCharged(0),
  fESDTrackCutsLowPt(0),
  fOutputContainer(0),
//thnf
/*
  fTHnFK0s(0),
  fTHnFL(0),
  fTHnFAL(0),

  fTHnFK0sDauEta(0),
  fTHnFLDauEta(0),
  fTHnFALDauEta(0),
  fTHnFK0sDauPhi(0),
  fTHnFLDauPhi(0),
  fTHnFALDauPhi(0),
*/
//event histos
  fHistITSLayerHits(0),
  fHistOneHitWithSDD(0),
  fHistNEvents(0),
  fHistPrimVtxZESDVSNContributors(0),
  fHistPrimVtxZESDTPCVSNContributors(0),
  fHistPrimVtxZESDSPDVSNContributors(0),
  fHistPrimVtxZESD(0),
  fHistPrimVtxZESDTPC(0),
  fHistPrimVtxZESDSPD(0),
  fHistESDVertexZ(0),
  fHistMuliplicity(0),
  fHistMuliplicityRaw(0),
  fHistCentBinRaw(0),
  fHistCentBin(0),
  fHistMultiplicityPrimary(0),
  fHistNPrim(0),
  fHistPiPiK0sVsLambdaMass(0),
  fHistPiPiK0sVsALambdaMass(0),
  fHistPiPK0sVsLambdaMass(0),
  fHistPiAPK0sVsALambdaMass(0),
  fHistPiPALambdaVsLambdaMass(0),
  fHistPiAPLambdaVsALambdaMass(0),
//-----------K0 histos -------------------//
  fHistPiPiMass(0),
  fHistPiPiMassVSPt(0),
  fHistPiPiMassVSPtMCTruth(0),
  fHistPiPiMassVSPtPosMCTruth(0),
  fHistPiPiMassVSPtNegMCTruth(0),
  fHistPiPiMassVSY(0),
  fHistPiPiPtVSY(0),
// fHistPiPiMassVSAlpha(0),
  fHistPiPiRadiusXY(0),
  fHistPiPiCosPointAng(0),
  fHistPiPiDCADaughterPosToPrimVtxVSMass(0),  
  fHistPiPiDecayLengthVsPt(0),
  fHistPiPiDecayLengthVsMass(0),
  fHistPiPiDecayLengthVsCtau(0),
  fHistPiPiDCADaughters(0), 
//  fHistPiPiPtDaughters(0),
  fHistPiPiDCAVSMass(0),
  fHistPiPiDCAZVSMass(0),
  fHistPiPiDCAZPos(0),
  fHistPiPiDCAZNeg(0),
  fHistPiPiTrackLengthPosVsMass(0),
  fHistPiPiTrackLengthNegVsMass(0),  
  fHistPiPiMonitorCuts(0),
  fHistPiPiMonitorMCCuts(0),
  fHistPiPiDecayLengthResolution(0),
  fHistNclsITSPosK0(0),
  fHistNclsITSNegK0(0),
  fHistNclsTPCPosK0(0),
  fHistNclsTPCNegK0(0),
  fHistChi2PerNclsITSPosK0(0),
  fHistChi2PerNclsITSNegK0(0),
  fHistNCRowsTPCPosK0(0),
  fHistNCRowsTPCNegK0(0),
  fHistRatioFoundOverFinableTPCK0Pos(0),
  fHistRatioFoundOverFinableTPCK0Neg(0),
//  fHistPiPiDistDaughtersTPCEntrVsMass(0),
//Lambda Antilambda
// fHistPiPDistDaughtersTPCEntrVsMass(0),
// fHistPiAPDistDaughtersTPCEntrVsMass(0),
//------------MC only histos-----------
  fHistPrimVtxZESDVSNContributorsMC(0),
  fHistPrimVtxZESDTPCVSNContributorsMC(0),
  fHistPrimVtxZESDSPDVSNContributorsMC(0),
  fHistMCVertexZ(0),
  fHistPiPiPDGCode(0),
  fHistPiPPDGCode(0),
  fHistPiAPPDGCode(0),

//-- BG of K0s
// fHistPiPiGA(0),
// fHistPiPiKch(0),
// fHistPiPiPhi(0),
// fHistPiPiL(0),
// fHistPiPiPi0(0),
// fHistPiPiPich(0),
// fHistPiPiRoh(0),
// fHistPiPiOmega(0),
// fHistPiPiKStar(0),
// fHistPiPiNoMother(0),
// fHistPiPiK0s(0),
// fHistPiPiK0L(0),
// fHistPiPiN(0),
// fHistPiPiSigma(0),
// fHistPiPiXi(0),
// fHistPiPiDelta(0),
// fHistPiPiB(0),
// fHistPiPiD(0),
// fHistPiPiEta(0),
// //-- BG of Lambda
// fHistPiPGA(0),
// fHistPiPKch(0),
// fHistPiPK0s(0),
// fHistPiPPi0(0),
// fHistPiPPich(0),
// fHistPiPKStar(0),
// fHistPiPN(0),
// fHistPiPNoMother(0),
// fHistPiPL(0),

//cosine of pointing angle of Xi vs pt histos
  fHistPiPCosPointAngXiVsPt(0),
  fHistPiAPCosPointAngXiVsPt(0),
  fHistPiPMassVSPtSecXiMCTruth(0),
  fHistPiPMassVSPtSecOmegaMCTruth(0),
  fHistPiAPMassVSPtSecXiMCTruth(0),
  fHistPiAPMassVSPtSecOmegaMCTruth(0),
// fHistUserPtShift(0),
// fHistPiPiPhiPosVsPtPosVsMass(0),//xxx
// fHistPiPPhiPosVsPtPosVsMass(0),//xxx
// fHistPiAPPhiPosVsPtPosVsMass(0),//xxx
//selection booleans and values
  fMCMode(0),
  fMCTruthMode(0),
  fSelectInjected(0),
  fSelectMBMotherMC(0),
  fCheckNegLabelReco(0),
  fOnlyFoundRecoV0(0),
  fUseCentrality(0),
  fUseCentralityBin(0),
  fUseCentralityRange(0),
  fAnapp(0),
  fRejectPileUpSPD(0),
  fSelSDD(0),
  fSelNoSDD(0),
  fOntheFly(0),
  fVertexZCut(0),
  fVtxStatus(0),
  fNcr(0),              
  fChi2cls(0),      
  fTPCrefit(0),
  fITSrefit(0),
  fNcrCh(0),      
  fChi2clsCh(0),         
  fTPCrefitCh(0),  
  fITSrefitCh(0),
  fNcrLpt(0),            
  fChi2clsLpt(0),     
  fTPCrefitLpt(0),
  fITSrefitLpt(0),
  fUsePID(0),
  fUsePIDPion(0),
  fNSigma(0),
  fNSigma2(0),
  fPPIDcut(0),
  fPtTPCCut(0),
  fMoreNclsThanRows(0),
  fMoreNclsThanFindable(0),
  fMoreNclsThanFindableMax(0),
  fRatioFoundOverFindable(0),
  fRatioMaxCRowsOverFindable(0),
  fChi2PerClusterITS(0),
  fDistanceTPCInner(0),
  fMinNCLSITSPos(0),
  fMinNCLSITSNeg(0),
  fMaxNCLSITSPos(0),
  fMaxNCLSITSNeg(0),
  fSwitchCaseITSCls(0),
  fCutMITrackLength(0),
  fCutMICrossedR(0),
  fCutMITPCncls(0),
  fCutMITrackLengthLengthF(0),
  fCutMICrossedRLengthF(0),
  fRapCutV0(0),
  fRap(0),
  fEtaCutMCDaughters(0),
  fEtaCutMCDaughtersVal(0),
  fEtaSignCut(0),  
  fUseXi0(0),
  fUseXiM(0),
  fUseOmega(0),
  fCutRapXi(0),
  fMinPt(0),  
  fAlfaCut(0),
  fQtCut(0),
  fQtCutPt(0),
  fQtCutPtLow(0),
  fArmCutK0(0),      
  fArmCutL(0),
  fArmQtSlope(0),
  fExcludeLambdaFromK0s(0),
  fExcludeK0sFromLambda(0),
  fExcludePhotonsFromK0s(0),
  fExcludePhotonsFromLambda(0),
  fDCAToVertexK0(0),
  fDCAToVertexL(0),
  fDCAXK(0),
  fDCAYK(0),
  fDCAXL(0),
  fDCAYL(0),
  fDCAZ(0),
  fDCADaughtersL(0),
  fDCADaughtersAL(0),
  fDCADaughtersK0(0),
  fDCADaughtersToVtxLarge(0),
  fDCADaughtersToVtxSmall(0),
  fDecayRadXYMin(0),
  fDecayRadXYMax(0),
  fPtDecRadMin(0),
  fCosPointAngL(0),
  fCosPointAngK(0),
  fCPAPtCutK0(0),
  fCPAPtCutL(0),
  fOpengAngleDaughters(0),
  fOpAngPtCut(0),
  fDecayLengthMax(0),
  fDecayLengthMin(0),
  fDecRadCutITSMin(0),
  fDecRadCutITSMax(0),
  fCtauK0s(0),
  fCtauL(0),
  fCtauPtCutK0(0),
  fCtauPtCutL(0),
  fChiCutKf(0),			
  fK0sLowMassCut(0),
  fK0sHighMassCut(0),
  fLLowMassCut(0),
  fLHighMassCut(0),
  fSetFillDetAL(0),
  fSetPtDepHist(0),
  fStopLoop(0),
  fDistDauForCheck(0),
  fShift(0),
  fDeltaInvP(0)
{  // Constructor.

  DefineOutput(1,TList::Class());
  // define defaults for globals
  
  fShift = kFALSE;                       // shift in charge/pt yes/no
  fDeltaInvP = 0.00;                     // shift value
    
   
  fMCMode = kFALSE;
  fMCTruthMode = kFALSE;

  fUseCentrality = 0;
  fUseCentralityBin = 0;
  fUseCentralityRange =0;

  fAnapp = kFALSE;
  fRejectPileUpSPD = kFALSE;
  fSelSDD = kFALSE;
  fSelNoSDD= kFALSE;
   
  fSelectInjected = kFALSE;
  fSelectMBMotherMC = kFALSE;
  fCheckNegLabelReco = kFALSE;
  fOnlyFoundRecoV0= kFALSE;

  fVertexZCut = 100000.0;
  fVtxStatus = kFALSE;

  fOntheFly = kTRUE;

  //----- define defaults for V0 and track cuts ----//
  fNcr = 70;              
  fChi2cls = 4;      
  fTPCrefit = kTRUE;   
  fITSrefit = kFALSE;
  fNcrCh = 70;      
  fChi2clsCh =4;         
  fTPCrefitCh = kTRUE; 
  fITSrefitCh = kFALSE;
  fNcrLpt = 70;            
  fChi2clsLpt = 4;     
  fTPCrefitLpt = kTRUE; 
  fITSrefitLpt = kFALSE; 

  fUsePID = kFALSE;
  fUsePIDPion = kFALSE;
  fMoreNclsThanRows = kFALSE;
  fMoreNclsThanFindable = kFALSE;
  fMoreNclsThanFindableMax = kFALSE;
  fRatioFoundOverFindable = -1.0;
  fRatioMaxCRowsOverFindable = 1000.0;


  fChi2PerClusterITS = 100000.0;
  fDistanceTPCInner = -1.0;
  fMinNCLSITSPos = -1;
  fMaxNCLSITSPos = 1000;
  fMinNCLSITSNeg = -1;
  fMaxNCLSITSNeg = 1000;
  fSwitchCaseITSCls = kFALSE;

  fCutMITrackLength = kFALSE;
  fCutMICrossedR    = kFALSE;
  fCutMITPCncls     = kFALSE;
  fCutMITrackLengthLengthF = 1.0;
  fCutMICrossedRLengthF = 0.85;

  fNSigma   = 100000.0;
  fNSigma2  = 100000.0;
  fPPIDcut  = 100.0;
  fPtTPCCut = -1.0;


  fRapCutV0=kFALSE;
  fRap=1000.0;
  fRap=1000.0;

  fAlfaCut= -100.0;
  fQtCut = -1.0;
  fQtCutPt = 100.0;
  fQtCutPtLow = -1.0;
  fArmCutK0=kFALSE;     
  fArmCutL=kFALSE;  
  fArmQtSlope =0.2;
  fExcludeLambdaFromK0s = -1.0;
  fExcludeK0sFromLambda = -1.0;
  fExcludePhotonsFromK0s = -1.0;
  fExcludePhotonsFromLambda = -1.0;

  fEtaCutMCDaughters = kFALSE;
  fEtaCutMCDaughtersVal = 50.0;
  fEtaSignCut = 0.0;

  fUseXi0= kTRUE;
  fUseXiM = kTRUE;
  fUseOmega =kTRUE;
  fCutRapXi = kFALSE;
  
  fMinPt= -1.0;

  fDCAToVertexK0 = 10000.0;
  fDCAToVertexL = 10000.0;
  fDCAXK=10000.0;
  fDCAYK=10000.0;
  fDCAXL=10000.0;
  fDCAYL=10000.0;
  fDCAZ=10000.0;
   
  fDCADaughtersL=10000.0;
  fDCADaughtersAL=10000.0;
  fDCADaughtersK0=10000.0;

  fDCADaughtersToVtxLarge=-1.0;
  fDCADaughtersToVtxSmall=-1.0;

  fDecayRadXYMin = -100000.0;
  fDecayRadXYMax = 1000000.0;
  fPtDecRadMin = 1000000.0;
  fDecayLengthMax = 100000.0;
  fDecayLengthMin = -1000000.0;
   
  fDecRadCutITSMin = 0.0000;
  fDecRadCutITSMax = 10000.0;

  fCosPointAngL=-1.0;
  fCosPointAngK=-1.0;
  fCPAPtCutK0 = 1000.0;
  fCPAPtCutL = -1000.0;//xxx
  fOpengAngleDaughters = -1.0;
  fOpAngPtCut = -1.0;
      
  fCtauK0s=10e6;
  fCtauL=10e6;
  fCtauPtCutK0=10e6;
  fCtauPtCutL=10e6;

  fChiCutKf=1000000.0;

  fK0sLowMassCut  = 0.25;
  fK0sHighMassCut = 0.75;

  fLLowMassCut  = 1.05;
  fLHighMassCut = 1.25;


  fSetFillDetAL = kFALSE;

  fSetPtDepHist=kFALSE;

  fStopLoop = kFALSE;

  fDistDauForCheck =5.0;

  //---- histograms ----//
  for(Int_t j=0;j<2;j++){
    fHistArmenteros[j]=NULL;
    fHistV0RadiusZ[j] =NULL;
    fHistV0RadiusZVSPt[j] =NULL;
    fHistV0RadiusXY[j] =NULL;
    fHistV0RadiusXYVSY[j] =NULL;
   
    //Lambda
    fHistPiPMass[j]=NULL;
    fHistPiPMassVSPt[j]=NULL;
    fHistPiPMassVSY[j] = NULL;
    fHistPiPMassVSPtMCTruth[j]=NULL;
    fHistPiPMassVSPtPosMCTruth[j]=NULL;
    fHistPiPMassVSPtNegMCTruth[j]=NULL;
    fHistPiPRadiusXY[j]=NULL;
    fHistPiPCosPointAng[j]=NULL;
    fHistPiPDecayLengthVsPt[j]=NULL;
    fHistPiPDecayLengthVsMass[j]=NULL;
    fHistPiPDecayLengthVsCtau[j]=NULL;
    fHistPiPDCADaughterPosToPrimVtxVSMass[j]=NULL;
    fHistPiPDCADaughterNegToPrimVtxVSMass[j]=NULL;
    fHistPiPMassVSPtSecSigma[j]=NULL;
    fHistPiPMassVSPtSecXi[j]=NULL;
    fHistPiPMassVSPtSecOmega[j]=NULL;
    fHistPiPMassVSYSecXi[j]=NULL;
    fHistPiPXi0PtVSLambdaPt[j]=NULL;
    fHistPiPXiMinusPtVSLambdaPt[j]=NULL;
    fHistPiPOmegaPtVSLambdaPt[j]=NULL;
    fHistPiPDCADaughters[j]=NULL;
    //  fHistPiPPtDaughters[j]=NULL;
    fHistPiPPtVSY[j]=NULL;
    fHistPiPDCAVSMass[j]=NULL;
    fHistPiPDCAZVSMass[j]=NULL;
    fHistPiPMonitorCuts[j] =NULL;
    fHistPiPMonitorMCCuts[j] =NULL;
    fHistPiPDecayLengthResolution[j] =NULL;
    fHistPiPDCAZPos[j] =NULL;
    fHistPiPDCAZNeg[j] =NULL;
    fHistPiPTrackLengthPosVsMass[j] = NULL;
    fHistPiPTrackLengthNegVsMass[j] = NULL;

    //ALambda
    fHistPiAPMass[j]=NULL;
    fHistPiAPMassVSPt[j]=NULL;
    fHistPiAPMassVSY[j] = NULL;
    fHistPiAPMassVSPtMCTruth[j]=NULL;
    fHistPiAPMassVSPtPosMCTruth[j]=NULL;
    fHistPiAPMassVSPtNegMCTruth[j]=NULL;
    fHistPiAPRadiusXY[j]=NULL;
    fHistPiAPCosPointAng[j]=NULL;
    fHistPiAPDecayLengthVsPt[j]=NULL;
    fHistPiAPDecayLengthVsMass[j]=NULL;
    fHistPiAPDecayLengthVsCtau[j]=NULL;
    fHistPiAPDCADaughterPosToPrimVtxVSMass[j]=NULL;
    fHistPiAPDCADaughterNegToPrimVtxVSMass[j]=NULL;
    fHistPiAPMassVSPtSecSigma[j]=NULL;
    fHistPiAPMassVSPtSecXi[j]=NULL;
    fHistPiAPMassVSPtSecOmega[j]=NULL;
    fHistPiAPMassVSYSecXi[j]=NULL;
    fHistPiAPXi0PtVSLambdaPt[j]=NULL;
    fHistPiAPXiMinusPtVSLambdaPt[j]=NULL;
    fHistPiAPOmegaPtVSLambdaPt[j] =NULL;
    fHistPiAPDCADaughters[j]=NULL;
    // fHistPiAPPtDaughters[j]=NULL;
    fHistPiAPPtVSY[j]=NULL;
    fHistPiAPDCAVSMass[j]=NULL;
    fHistPiAPDCAZVSMass[j]=NULL;
    fHistPiAPMonitorCuts[j] =NULL;
    fHistPiAPMonitorMCCuts[j] =NULL;
    fHistPiAPDecayLengthResolution[j] =NULL;
    //    fHistPiAPDCAZPos[j] =NULL;
    //fHistPiAPDCAZNeg[j] =NULL;
    fHistPiAPTrackLengthPosVsMass[j] = NULL;
    fHistPiAPTrackLengthNegVsMass[j] = NULL;

    //other 
    fHistDedxSecProt[j]=NULL;
    fHistDedxSecAProt[j]=NULL;
    fHistDedxSecPiMinus[j]=NULL;
    fHistDedxSecPiPlus[j]=NULL;
    fHistDedxProt[j]=NULL;
    fHistDedxAProt[j]=NULL;
    fHistDedxPiMinus[j]=NULL;
    fHistDedxPiPlus[j]=NULL;
    fHistNclsITS[j]=NULL;
    fHistNclsTPC[j]=NULL;
    fHistNclsITSPosL[j]=NULL;
    fHistNclsITSNegL[j]=NULL;
    fHistNclsTPCPosL[j]=NULL;
    fHistNclsTPCNegL[j]=NULL;
    fHistChi2PerNclsITSPosL[j]=NULL;
    fHistChi2PerNclsITSNegL[j]=NULL;
    fHistNCRowsTPCPosL[j]=NULL;
    fHistNCRowsTPCNegL[j]=NULL;
    fHistRatioFoundOverFinableTPCLPos[j]=NULL;
    fHistRatioFoundOverFinableTPCLNeg[j]=NULL;
    fHistPiPiEtaDMC[j] = NULL;
    fHistPiPiEtaDReco[j] = NULL;
    fHistPiPEtaDMC[j] = NULL;
    fHistPiPEtaDReco[j] = NULL;
  }
  /*
    for(Int_t m=0;m<3;m++){

    fHistPiPiDistDaughtersPos[m]=NULL;
    fHistPiPiDistDaughtersNeg[m]=NULL;
    fHistPiPiDCADaughtersPos[m]=NULL;
    fHistPiPiDCADaughtersNeg[m]=NULL;
    fHistPiPiRadAtDCA5cmDaughtersPos[m]=NULL;
    fHistPiPiRadAtDCA5cmDaughtersNeg[m]=NULL;
    
    fHistPiPDistDaughtersPos[m]=NULL;
    fHistPiPDistDaughtersNeg[m]=NULL;
    fHistPiPDCADaughtersPos[m]=NULL;
    fHistPiPDCADaughtersNeg[m]=NULL;
    fHistPiPRadAtDCA5cmDaughtersPos[m]=NULL;
    fHistPiPRadAtDCA5cmDaughtersNeg[m]=NULL;

    fHistPiAPDistDaughtersPos[m]=NULL;
    fHistPiAPDistDaughtersNeg[m]=NULL;
    fHistPiAPDCADaughtersPos[m]=NULL;
    fHistPiAPDCADaughtersNeg[m]=NULL;
    fHistPiAPRadAtDCA5cmDaughtersPos[m]=NULL;
    fHistPiAPRadAtDCA5cmDaughtersNeg[m]=NULL;
    }
  */
}
//________________________________________________________________________
AliAnalysisTaskV0ForRAA::AliAnalysisTaskV0ForRAA(const char *name)
  :AliAnalysisTaskSE(name),
  fESD(0),
  fMCev(0),
//other objects
  fESDpid(0),
  fESDTrackCuts(0),
  fESDTrackCutsCharged(0),
  fESDTrackCutsLowPt(0),
  fOutputContainer(0),
//thnf
/*
  fTHnFK0s(0),
  fTHnFL(0),
  fTHnFAL(0),

  fTHnFK0sDauEta(0),
  fTHnFLDauEta(0),
  fTHnFALDauEta(0),
  fTHnFK0sDauPhi(0),
  fTHnFLDauPhi(0),
  fTHnFALDauPhi(0),
*/
//event histos
  fHistITSLayerHits(0),
  fHistOneHitWithSDD(0),
  fHistNEvents(0),
  fHistPrimVtxZESDVSNContributors(0),
  fHistPrimVtxZESDTPCVSNContributors(0),
  fHistPrimVtxZESDSPDVSNContributors(0),
  fHistPrimVtxZESD(0),
  fHistPrimVtxZESDTPC(0),
  fHistPrimVtxZESDSPD(0),
  fHistESDVertexZ(0),
  fHistMuliplicity(0),
  fHistMuliplicityRaw(0),
  fHistCentBinRaw(0),
  fHistCentBin(0),
  fHistMultiplicityPrimary(0),
  fHistNPrim(0),
  fHistPiPiK0sVsLambdaMass(0),
  fHistPiPiK0sVsALambdaMass(0),
  fHistPiPK0sVsLambdaMass(0),
  fHistPiAPK0sVsALambdaMass(0),
  fHistPiPALambdaVsLambdaMass(0),
  fHistPiAPLambdaVsALambdaMass(0),
//-----------K0 histos -------------------//
  fHistPiPiMass(0),
  fHistPiPiMassVSPt(0),
  fHistPiPiMassVSPtMCTruth(0),
  fHistPiPiMassVSPtPosMCTruth(0),
  fHistPiPiMassVSPtNegMCTruth(0),
  fHistPiPiMassVSY(0),
  fHistPiPiPtVSY(0),
// fHistPiPiMassVSAlpha(0),
  fHistPiPiRadiusXY(0),
  fHistPiPiCosPointAng(0),
  fHistPiPiDCADaughterPosToPrimVtxVSMass(0),  
  fHistPiPiDecayLengthVsPt(0),
  fHistPiPiDecayLengthVsMass(0),
  fHistPiPiDecayLengthVsCtau(0),
  fHistPiPiDCADaughters(0), 
//  fHistPiPiPtDaughters(0),
  fHistPiPiDCAVSMass(0),
  fHistPiPiDCAZVSMass(0),
  fHistPiPiDCAZPos(0),
  fHistPiPiDCAZNeg(0),
  fHistPiPiTrackLengthPosVsMass(0),
  fHistPiPiTrackLengthNegVsMass(0),  
  fHistPiPiMonitorCuts(0),
  fHistPiPiMonitorMCCuts(0),
  fHistPiPiDecayLengthResolution(0),
  fHistNclsITSPosK0(0),
  fHistNclsITSNegK0(0),
  fHistNclsTPCPosK0(0),
  fHistNclsTPCNegK0(0),
  fHistChi2PerNclsITSPosK0(0),
  fHistChi2PerNclsITSNegK0(0),
  fHistNCRowsTPCPosK0(0),
  fHistNCRowsTPCNegK0(0),
  fHistRatioFoundOverFinableTPCK0Pos(0),
  fHistRatioFoundOverFinableTPCK0Neg(0),
   //  fHistPiPiDistDaughtersTPCEntrVsMass(0),
//Lambda Antilambda
   // fHistPiPDistDaughtersTPCEntrVsMass(0),
   // fHistPiAPDistDaughtersTPCEntrVsMass(0),
//------------MC only histos-----------
  fHistPrimVtxZESDVSNContributorsMC(0),
  fHistPrimVtxZESDTPCVSNContributorsMC(0),
  fHistPrimVtxZESDSPDVSNContributorsMC(0),
  fHistMCVertexZ(0),
  fHistPiPiPDGCode(0),
  fHistPiPPDGCode(0),
  fHistPiAPPDGCode(0),

//-- BG of K0s
// fHistPiPiGA(0),
// fHistPiPiKch(0),
// fHistPiPiPhi(0),
// fHistPiPiL(0),
// fHistPiPiPi0(0),
// fHistPiPiPich(0),
// fHistPiPiRoh(0),
// fHistPiPiOmega(0),
// fHistPiPiKStar(0),
// fHistPiPiNoMother(0),
// fHistPiPiK0s(0),
// fHistPiPiK0L(0),
// fHistPiPiN(0),
// fHistPiPiSigma(0),
// fHistPiPiXi(0),
// fHistPiPiDelta(0),
// fHistPiPiB(0),
// fHistPiPiD(0),
// fHistPiPiEta(0),
// //-- BG of Lambda
// fHistPiPGA(0),
// fHistPiPKch(0),
// fHistPiPK0s(0),
// fHistPiPPi0(0),
// fHistPiPPich(0),
// fHistPiPKStar(0),
// fHistPiPN(0),
// fHistPiPNoMother(0),
// fHistPiPL(0),

//cosine of pointing angle of Xi vs pt histos
  fHistPiPCosPointAngXiVsPt(0),
  fHistPiAPCosPointAngXiVsPt(0),
  fHistPiPMassVSPtSecXiMCTruth(0),
  fHistPiPMassVSPtSecOmegaMCTruth(0),
  fHistPiAPMassVSPtSecXiMCTruth(0),
  fHistPiAPMassVSPtSecOmegaMCTruth(0),
// fHistUserPtShift(0),
// fHistPiPiPhiPosVsPtPosVsMass(0),//xxx
// fHistPiPPhiPosVsPtPosVsMass(0),//xxx
// fHistPiAPPhiPosVsPtPosVsMass(0),//xxx
//selection booleans and values
  fMCMode(0),
  fMCTruthMode(0),
  fSelectInjected(0),
  fSelectMBMotherMC(0),
  fCheckNegLabelReco(0),
  fOnlyFoundRecoV0(0),
  fUseCentrality(0),
  fUseCentralityBin(0),
  fUseCentralityRange(0),
  fAnapp(0),
  fRejectPileUpSPD(0),
  fSelSDD(0),
  fSelNoSDD(0),
  fOntheFly(0),
  fVertexZCut(0),
  fVtxStatus(0),
  fNcr(0),              
  fChi2cls(0),      
  fTPCrefit(0),
  fITSrefit(0),
  fNcrCh(0),      
  fChi2clsCh(0),         
  fTPCrefitCh(0),  
  fITSrefitCh(0),
  fNcrLpt(0),            
  fChi2clsLpt(0),     
  fTPCrefitLpt(0),
  fITSrefitLpt(0),
  fUsePID(0),
  fUsePIDPion(0),
  fNSigma(0),
  fNSigma2(0),
  fPPIDcut(0),
  fPtTPCCut(0),
  fMoreNclsThanRows(0),
  fMoreNclsThanFindable(0),
  fMoreNclsThanFindableMax(0),
  fRatioFoundOverFindable(0),
  fRatioMaxCRowsOverFindable(0),
  fChi2PerClusterITS(0),
  fDistanceTPCInner(0),
  fMinNCLSITSPos(0),
  fMinNCLSITSNeg(0),
  fMaxNCLSITSPos(0),
  fMaxNCLSITSNeg(0),
  fSwitchCaseITSCls(0),
  fCutMITrackLength(0),
  fCutMICrossedR(0),
  fCutMITPCncls(0),
  fCutMITrackLengthLengthF(0),
  fCutMICrossedRLengthF(0),
  fRapCutV0(0),
  fRap(0),
  fEtaCutMCDaughters(0),
  fEtaCutMCDaughtersVal(0),
  fEtaSignCut(0),    
  fUseXi0(0),
  fUseXiM(0),
  fUseOmega(0),
  fCutRapXi(0),
  fMinPt(0),  
  fAlfaCut(0),
  fQtCut(0),
  fQtCutPt(0),
  fQtCutPtLow(0),
  fArmCutK0(0),      
  fArmCutL(0),
  fArmQtSlope(0),
  fExcludeLambdaFromK0s(0),
  fExcludeK0sFromLambda(0),
  fExcludePhotonsFromK0s(0),
  fExcludePhotonsFromLambda(0),
  fDCAToVertexK0(0),
  fDCAToVertexL(0),
  fDCAXK(0),
  fDCAYK(0),
  fDCAXL(0),
  fDCAYL(0),
  fDCAZ(0),
  fDCADaughtersL(0),
  fDCADaughtersAL(0),
  fDCADaughtersK0(0),
  fDCADaughtersToVtxLarge(0),
  fDCADaughtersToVtxSmall(0),
  fDecayRadXYMin(0),
  fDecayRadXYMax(0),
  fPtDecRadMin(0),
  fCosPointAngL(0),
  fCosPointAngK(0),
  fCPAPtCutK0(0),
  fCPAPtCutL(0),
  fOpengAngleDaughters(0),
  fOpAngPtCut(0),
  fDecayLengthMax(0),
  fDecayLengthMin(0),
  fDecRadCutITSMin(0),
  fDecRadCutITSMax(0),
  fCtauK0s(0),
  fCtauL(0),
  fCtauPtCutK0(0),
  fCtauPtCutL(0),
  fChiCutKf(0),			
  fK0sLowMassCut(0),
  fK0sHighMassCut(0),
  fLLowMassCut(0),
  fLHighMassCut(0),
  fSetFillDetAL(0),
  fSetPtDepHist(0),
  fStopLoop(0),
  fDistDauForCheck(0),
  fShift(0),
  fDeltaInvP(0)
{  // Constructor.

  DefineOutput(1,TList::Class());
  // define defaults for globals
  
  fShift = kFALSE;                       // shift in charge/pt yes/no
  fDeltaInvP = 0.00;                     // shift value
    
   
  fMCMode = kFALSE;
  fMCTruthMode = kFALSE;

  fUseCentrality = 0;
  fUseCentralityBin = 0;
  fUseCentralityRange =0;

  fAnapp = kFALSE;
  fRejectPileUpSPD = kFALSE;
  fSelSDD = kFALSE;
  fSelNoSDD= kFALSE;
   
  fSelectInjected = kFALSE;
  fSelectMBMotherMC = kFALSE;
  fCheckNegLabelReco = kFALSE;
  fOnlyFoundRecoV0= kFALSE;

  fVertexZCut = 100000.0;
  fVtxStatus = kFALSE;

  fOntheFly = kTRUE;

  //----- define defaults for V0 and track cuts ----//
  fNcr = 70;              
  fChi2cls = 4;      
  fTPCrefit = kTRUE;   
  fITSrefit = kFALSE;
  fNcrCh = 70;      
  fChi2clsCh =4;         
  fTPCrefitCh = kTRUE; 
  fITSrefitCh = kFALSE;
  fNcrLpt = 70;            
  fChi2clsLpt = 4;     
  fTPCrefitLpt = kTRUE; 
  fITSrefitLpt = kFALSE; 

  fUsePID = kFALSE;
  fUsePIDPion = kFALSE;
  fMoreNclsThanRows = kFALSE;
  fMoreNclsThanFindable = kFALSE;
  fMoreNclsThanFindableMax = kFALSE;
  fRatioFoundOverFindable = -1.0;
  fRatioMaxCRowsOverFindable = 1000.0;


  fChi2PerClusterITS = 100000.0;
  fDistanceTPCInner = -1.0;
  fMinNCLSITSPos = -1;
  fMaxNCLSITSPos = 1000;
  fMinNCLSITSNeg = -1;
  fMaxNCLSITSNeg = 1000;
  fSwitchCaseITSCls = kFALSE;

  fCutMITrackLength = kFALSE;
  fCutMICrossedR    = kFALSE;
  fCutMITPCncls     = kFALSE;
  fCutMITrackLengthLengthF = 1.0;
  fCutMICrossedRLengthF = 0.85;

  fNSigma   = 100000.0;
  fNSigma2  = 100000.0;
  fPPIDcut  = 100.0;
  fPtTPCCut = -1.0;


  fRapCutV0=kFALSE;
  fRap=1000.0;
  fRap=1000.0;

  fAlfaCut= -100.0;
  fQtCut = -1.0;
  fQtCutPt = 100.0;
  fQtCutPtLow = -1.0;
  fArmCutK0=kFALSE;     
  fArmCutL=kFALSE;  
  fArmQtSlope =0.2;
  fExcludeLambdaFromK0s = -1.0;
  fExcludeK0sFromLambda = -1.0;
  fExcludePhotonsFromK0s = -1.0;
  fExcludePhotonsFromLambda = -1.0;

  fEtaCutMCDaughters = kFALSE;
  fEtaCutMCDaughtersVal = 50.0;
  fEtaSignCut = 0.0;

  fUseXi0= kTRUE;
  fUseXiM = kTRUE;
  fUseOmega =kTRUE;
  fCutRapXi = kFALSE;
  
  fMinPt= -1.0;

  fDCAToVertexK0 = 10000.0;
  fDCAToVertexL = 10000.0;
  fDCAXK=10000.0;
  fDCAYK=10000.0;
  fDCAXL=10000.0;
  fDCAYL=10000.0;
  fDCAZ=10000.0;
   
  fDCADaughtersL=10000.0;
  fDCADaughtersAL=10000.0;
  fDCADaughtersK0=10000.0;

  fDCADaughtersToVtxLarge=-1.0;
  fDCADaughtersToVtxSmall=-1.0;

  fDecayRadXYMin = -100000.0;
  fDecayRadXYMax = 1000000.0;
  fPtDecRadMin = 1000000.0;
  fDecayLengthMax = 100000.0;
  fDecayLengthMin = -1000000.0;
   
  fDecRadCutITSMin = 0.0000;
  fDecRadCutITSMax = 10000.0;

  fCosPointAngL=-1.0;
  fCosPointAngK=-1.0;
  fCPAPtCutK0 = 1000.0;
  fCPAPtCutL = -1000.0;//xxx
  fOpengAngleDaughters = -1.0;
  fOpAngPtCut = -1.0;
      
  fCtauK0s=10e6;
  fCtauL=10e6;
  fCtauPtCutK0=10e6;
  fCtauPtCutL=10e6;

  fChiCutKf=1000000.0;

  fK0sLowMassCut  = 0.25;
  fK0sHighMassCut = 0.75;

  fLLowMassCut  = 1.05;
  fLHighMassCut = 1.25;


  fSetFillDetAL = kFALSE;

  fSetPtDepHist=kFALSE;

  fStopLoop = kFALSE;

  fDistDauForCheck =5.0;

  //---- histograms ----//
  for(Int_t j=0;j<2;j++){
    fHistArmenteros[j]=NULL;
    fHistV0RadiusZ[j] =NULL;
    fHistV0RadiusZVSPt[j] =NULL;
    fHistV0RadiusXY[j] =NULL;
    fHistV0RadiusXYVSY[j] =NULL;
   
    //Lambda
    fHistPiPMass[j]=NULL;
    fHistPiPMassVSPt[j]=NULL;
    fHistPiPMassVSY[j] = NULL;
    fHistPiPMassVSPtMCTruth[j]=NULL;
    fHistPiPMassVSPtPosMCTruth[j]=NULL;
    fHistPiPMassVSPtNegMCTruth[j]=NULL;
    fHistPiPRadiusXY[j]=NULL;
    fHistPiPCosPointAng[j]=NULL;
    fHistPiPDecayLengthVsPt[j]=NULL;
    fHistPiPDecayLengthVsMass[j]=NULL;
    fHistPiPDecayLengthVsCtau[j]=NULL;
    fHistPiPDCADaughterPosToPrimVtxVSMass[j]=NULL;
    fHistPiPDCADaughterNegToPrimVtxVSMass[j]=NULL;
    fHistPiPMassVSPtSecSigma[j]=NULL;
    fHistPiPMassVSPtSecXi[j]=NULL;
    fHistPiPMassVSPtSecOmega[j]=NULL;
    fHistPiPMassVSYSecXi[j]=NULL;
    fHistPiPXi0PtVSLambdaPt[j]=NULL;
    fHistPiPXiMinusPtVSLambdaPt[j]=NULL;
    fHistPiPOmegaPtVSLambdaPt[j]=NULL;
    fHistPiPDCADaughters[j]=NULL;
    //  fHistPiPPtDaughters[j]=NULL;
    fHistPiPPtVSY[j]=NULL;
    fHistPiPDCAVSMass[j]=NULL;
    fHistPiPDCAZVSMass[j]=NULL;
    fHistPiPMonitorCuts[j] =NULL;
    fHistPiPMonitorMCCuts[j] =NULL;
    fHistPiPDecayLengthResolution[j] =NULL;
    fHistPiPDCAZPos[j] =NULL;
    fHistPiPDCAZNeg[j] =NULL;
    fHistPiPTrackLengthPosVsMass[j] = NULL;
    fHistPiPTrackLengthNegVsMass[j] = NULL;

    //ALambda
    fHistPiAPMass[j]=NULL;
    fHistPiAPMassVSPt[j]=NULL;
    fHistPiAPMassVSY[j] = NULL;
    fHistPiAPMassVSPtMCTruth[j]=NULL;
    fHistPiAPMassVSPtPosMCTruth[j]=NULL;
    fHistPiAPMassVSPtNegMCTruth[j]=NULL;
    fHistPiAPRadiusXY[j]=NULL;
    fHistPiAPCosPointAng[j]=NULL;
    fHistPiAPDecayLengthVsPt[j]=NULL;
    fHistPiAPDecayLengthVsMass[j]=NULL;
    fHistPiAPDecayLengthVsCtau[j]=NULL;
    fHistPiAPDCADaughterPosToPrimVtxVSMass[j]=NULL;
    fHistPiAPDCADaughterNegToPrimVtxVSMass[j]=NULL;
    fHistPiAPMassVSPtSecSigma[j]=NULL;
    fHistPiAPMassVSPtSecXi[j]=NULL;
    fHistPiAPMassVSPtSecOmega[j]=NULL;
    fHistPiAPMassVSYSecXi[j]=NULL;
    fHistPiAPXi0PtVSLambdaPt[j]=NULL;
    fHistPiAPXiMinusPtVSLambdaPt[j]=NULL;
    fHistPiAPOmegaPtVSLambdaPt[j] =NULL;
    fHistPiAPDCADaughters[j]=NULL;
    // fHistPiAPPtDaughters[j]=NULL;
    fHistPiAPPtVSY[j]=NULL;
    fHistPiAPDCAVSMass[j]=NULL;
    fHistPiAPDCAZVSMass[j]=NULL;
    fHistPiAPMonitorCuts[j] =NULL;
    fHistPiAPMonitorMCCuts[j] =NULL;
    fHistPiAPDecayLengthResolution[j] =NULL;
    //    fHistPiAPDCAZPos[j] =NULL;
    //fHistPiAPDCAZNeg[j] =NULL;
    fHistPiAPTrackLengthPosVsMass[j] = NULL;
    fHistPiAPTrackLengthNegVsMass[j] = NULL;

    //other 
    fHistDedxSecProt[j]=NULL;
    fHistDedxSecAProt[j]=NULL;
    fHistDedxSecPiMinus[j]=NULL;
    fHistDedxSecPiPlus[j]=NULL;
    fHistDedxProt[j]=NULL;
    fHistDedxAProt[j]=NULL;
    fHistDedxPiMinus[j]=NULL;
    fHistDedxPiPlus[j]=NULL;
    fHistNclsITS[j]=NULL;
    fHistNclsTPC[j]=NULL;
    fHistNclsITSPosL[j]=NULL;
    fHistNclsITSNegL[j]=NULL;
    fHistNclsTPCPosL[j]=NULL;
    fHistNclsTPCNegL[j]=NULL;
    fHistChi2PerNclsITSPosL[j]=NULL;
    fHistChi2PerNclsITSNegL[j]=NULL;
    fHistNCRowsTPCPosL[j]=NULL;
    fHistNCRowsTPCNegL[j]=NULL;
    fHistRatioFoundOverFinableTPCLPos[j]=NULL;
    fHistRatioFoundOverFinableTPCLNeg[j]=NULL;
    fHistPiPiEtaDMC[j] = NULL;
    fHistPiPiEtaDReco[j] = NULL;
    fHistPiPEtaDMC[j] = NULL;
    fHistPiPEtaDReco[j] = NULL;
  }
  /*
    for(Int_t m=0;m<3;m++){

    fHistPiPiDistDaughtersPos[m]=NULL;
    fHistPiPiDistDaughtersNeg[m]=NULL;
    fHistPiPiDCADaughtersPos[m]=NULL;
    fHistPiPiDCADaughtersNeg[m]=NULL;
    fHistPiPiRadAtDCA5cmDaughtersPos[m]=NULL;
    fHistPiPiRadAtDCA5cmDaughtersNeg[m]=NULL;
    
    fHistPiPDistDaughtersPos[m]=NULL;
    fHistPiPDistDaughtersNeg[m]=NULL;
    fHistPiPDCADaughtersPos[m]=NULL;
    fHistPiPDCADaughtersNeg[m]=NULL;
    fHistPiPRadAtDCA5cmDaughtersPos[m]=NULL;
    fHistPiPRadAtDCA5cmDaughtersNeg[m]=NULL;

    fHistPiAPDistDaughtersPos[m]=NULL;
    fHistPiAPDistDaughtersNeg[m]=NULL;
    fHistPiAPDCADaughtersPos[m]=NULL;
    fHistPiAPDCADaughtersNeg[m]=NULL;
    fHistPiAPRadAtDCA5cmDaughtersPos[m]=NULL;
    fHistPiAPRadAtDCA5cmDaughtersNeg[m]=NULL;
    }
  */ 
}
//_____________________________________________________
AliAnalysisTaskV0ForRAA::~AliAnalysisTaskV0ForRAA()
{
  //---- Remove all pointers ----//
  if(fOutputContainer) delete fOutputContainer;fOutputContainer=0;
  if(fESDTrackCuts) delete fESDTrackCuts;fESDTrackCuts=0;
  if(fESDTrackCutsCharged) delete fESDTrackCutsCharged;fESDTrackCutsCharged=0;
  if(fESDTrackCutsLowPt) delete fESDTrackCutsLowPt; fESDTrackCutsLowPt=0;
}
//________________________________________________________________________
void AliAnalysisTaskV0ForRAA::UserCreateOutputObjects(){

  //--- esd track cuts V0 daughters ---//
  TString cutsname = "esdtrackcuts";
  // esd track cuts for pions high pt
  fESDTrackCuts = new AliESDtrackCuts(cutsname);
  fESDTrackCuts->SetMaxChi2PerClusterTPC(fChi2cls);
  fESDTrackCuts->SetMinNCrossedRowsTPC(fNcr);
  fESDTrackCuts->SetAcceptKinkDaughters(kFALSE);
  fESDTrackCuts->SetRequireTPCRefit(fTPCrefit);
  fESDTrackCuts->SetRequireSigmaToVertex(kFALSE);
  fESDTrackCuts->SetRequireITSRefit(fITSrefit);


  // esd track cuts for protons high pt
  TString cutsnameCh = cutsname;
  cutsnameCh +="_charged";
  fESDTrackCutsCharged = new AliESDtrackCuts(cutsnameCh);
  fESDTrackCutsCharged->SetMaxChi2PerClusterTPC(fChi2clsCh);
  fESDTrackCutsCharged->SetMinNCrossedRowsTPC(fNcrCh);
  fESDTrackCutsCharged->SetAcceptKinkDaughters(kFALSE);
  fESDTrackCutsCharged->SetRequireTPCRefit(fTPCrefitCh);
  fESDTrackCutsCharged->SetRequireSigmaToVertex(kFALSE);
  fESDTrackCutsCharged->SetRequireITSRefit(fITSrefit);
 
  // esd track cuts for all low pt
  TString cutsnameLowPt  = cutsname;
  cutsnameLowPt +="_lowpt";
  fESDTrackCutsLowPt = new AliESDtrackCuts(cutsnameLowPt);
  fESDTrackCutsLowPt->SetMaxChi2PerClusterTPC(fChi2clsLpt);
  fESDTrackCutsLowPt->SetMinNCrossedRowsTPC(fNcrLpt);
  fESDTrackCutsLowPt->SetAcceptKinkDaughters(kFALSE);
  fESDTrackCutsLowPt->SetRequireTPCRefit(fTPCrefitLpt);
  fESDTrackCutsLowPt->SetRequireSigmaToVertex(kFALSE);  
  fESDTrackCutsLowPt->SetRequireITSRefit(fITSrefit);


  //create output objects
  Int_t nbMass=500;
  
  //-----------------  create output container -----------------//

  fOutputContainer = new TList() ;
  fOutputContainer->SetName(GetName()) ;
  fOutputContainer->SetOwner();
  
  //  TH1::SetDefaultSumw2();
  //  TH2::SetDefaultSumw2();

  Int_t mchist = 1;// for Data
  if((fMCMode && fMCTruthMode) || fMCTruthMode) mchist = 2;//for MC to create sec. Lambda histos	

  //------------ create allways -----------------------//
  fHistNEvents = new TH1F("fHistNEvents","no of events before cuts =0, after cuts=1, after process =2",5,0.0,5.0);
  fOutputContainer->Add(fHistNEvents);
      
  fHistMuliplicity =  new TH1F("fHistMuliplicity","V0 multiplicity",3000,0.0,30000);
  fOutputContainer->Add(fHistMuliplicity);
      
  fHistMuliplicityRaw =  new TH1F("fHistMuliplicityRaw","V0 multiplicity before process",3000,0.0,30000);      
  fOutputContainer->Add(fHistMuliplicityRaw);
      
  fHistMultiplicityPrimary = new TH1F("fHistMultiplicityPrimary","number of charged tracks",5000,0.0,20000);
  fOutputContainer->Add(fHistMultiplicityPrimary);
      
  fHistESDVertexZ= new TH1F("fHistESDVertexZ"," z vertex distr in cm",500,-50,50);
  fOutputContainer->Add(fHistESDVertexZ);
   
  fHistPrimVtxZESD = new TH1F("fHistPrimVtxZESD","z vertex pos ESD",250,-50,50);
  fOutputContainer->Add(fHistPrimVtxZESD);
      
  fHistPrimVtxZESDVSNContributors = new TH2F("fHistPrimVtxZESDVSNContributors","prim vtx pos z ESD vs no. of contributers TPC",250,-50,50,500,0.0,500.0);
  fOutputContainer->Add(fHistPrimVtxZESDVSNContributors);
      
  fHistNPrim = new TH1F("fHistNPrim","Number of contributers to vertex",2500,0.0,5000);
  fOutputContainer->Add(fHistNPrim);
 
  //------------------------ pp analysis only -------------------------//
  if(fAnapp){
    fHistITSLayerHits = new TH1F("fHistITSLayerHits","SDD layer -1=0,1=1,2=2 ... 5=5,0=nothing",7,-1.5,5.5);
    fOutputContainer->Add(fHistITSLayerHits);
    fHistOneHitWithSDD = new TH1F("fHistOneHitWithSDD","min one hit in SDD",2,-0.5,1.5);
    fOutputContainer->Add(fHistOneHitWithSDD);
    fHistPrimVtxZESDTPC = new TH1F("fHistPrimVtxZESDTPC","z vertex pos TPC",250,-50,50);
    fOutputContainer->Add(fHistPrimVtxZESDTPC);
    fHistPrimVtxZESDSPD = new TH1F("fHistPrimVtxZESDSPD","z vertex pos SPD",250,-50,50);
    fOutputContainer->Add(fHistPrimVtxZESDSPD);  
    fHistPrimVtxZESDTPCVSNContributors = new TH2F("fHistPrimVtxZESDTPCVSNContributors","prim vtx pos z TPC vs no. of contributers TPC",250,-50,50,500,0.0,500.0);
    fOutputContainer->Add(fHistPrimVtxZESDTPCVSNContributors);
    fHistPrimVtxZESDSPDVSNContributors = new TH2F("fHistPrimVtxZESDSPDVSNContributors","prim vtx pos z SPD vs no. of contributers TPC",250,-50,50,500,0.0,500.0);
    fOutputContainer->Add(fHistPrimVtxZESDSPDVSNContributors);

  }
  else {
    Double_t binsCent[12]={0.0,5.0,10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0};
    fHistCentBinRaw = new TH1F("fHistCentBinRaw","centrality bin before cent selection",11,binsCent);
    fOutputContainer->Add(fHistCentBinRaw);
    fHistCentBin = new TH1F("fHistCentBin","centrality bin",11,binsCent);
    fOutputContainer->Add(fHistCentBin);
      
  }
   
  // ------------------- add always ---------------------------//
  //THnF

  // Double_t piForAx = 2.0*TMath::Pi();
  // Int_t binsTHnV0K0s[4] = {100,40,15,40};
  // Int_t binsTHnV0L[4] = {100,40,15,40};
  /*
    Int_t binsTHnV0DauEtaK0s[4] = {150,100,40,40};
    Int_t binsTHnV0DauEtaL[4] = {100,100,40,40};

    Int_t binsTHnV0DauPhiK0s[5] = {150, 18,18, 7,7};
    Int_t binsTHnV0DauPhiL[5] = {100, 18,18, 7,7};
  */

  // Double_t minK0s[4] = {0.35,0.0,0.0,0.0};
  // Double_t maxK0s[4] = {0.65,20.0,60.0,20.0};
  /*
    Double_t minK0sDauEta[4] = {0.35, 0.0,-0.8,-0.8};
    Double_t maxK0sDauEta[4] = {0.65,50.0, 0.8, 0.8};
    Double_t minK0sDauPhi[5] = {0.35,0.0,0.0,-0.5,-0.5};
    Double_t maxK0sDauPhi[5] = {0.65,piForAx,piForAx,6.5,6.5};
  */

  // Double_t minL[4] = {1.07, 0.0,0.0,0.0};
  // Double_t maxL[4] = {1.17,20.0, 60.0,20.0};
  /*
    Double_t minLDauEta[4] = {1.07, 0.0,-0.8,-0.8};
    Double_t maxLDauEta[4] = {1.17,50.0, 0.8, 0.8};
    Double_t minLDauPhi[5] = {1.07,0.0,0.0,-0.5,-0.5};
    Double_t maxLDauPhi[5] = {1.17,piForAx,piForAx,6.5, 6.5};
  */
  
  // char histTitK0s[255];
  // snprintf(histTitK0s,255,"fTHnFK0s");
  // char histTitL[255];
  // snprintf(histTitL,255,"fTHnFL");
  // char histTitAL[255];
  // snprintf(histTitAL,255,"fTHnFAL");

  /*
    char histTitK0sDauEta[255];
    snprintf(histTitK0sDauEta,255,"fTHnFK0sDauEta");
    char histTitLDauEta[255];
    snprintf(histTitLDauEta,255,"fTHnFLDauEta");
    char histTitALDauEta[255];
    snprintf(histTitALDauEta,255,"fTHnFALDauEta");


    char histTitK0sDauPhi[255];
    snprintf(histTitK0sDauPhi,255,"fTHnFK0sDauPhi");
    char histTitLDauPhi[255];
    snprintf(histTitLDauPhi,255,"fTHnFLDauPhi");
    char histTitALDauPhi[255];
    snprintf(histTitALDauPhi,255,"fTHnFALDauPhi");
  */
  // char axTitK0s[255];
  // snprintf(axTitK0s,255,"K^{0}_{s};m_{inv} (GeV/c^{2});p_{T} (GeV/c);#eta(V0);#phi(V0)");
  // char axTitL[255];
  // snprintf(axTitL,255,"#Lambda;m_{inv} (GeV/c^{2});p_{T} (GeV/c);#eta(V0);#phi(V0)");
  // char axTitAL[255];
  // snprintf(axTitAL,255,"#bar{#Lambda};m_{inv} (GeV/c^{2});p_{T} (GeV/c);#eta(V0);#phi(V0)");
 
  /*
    char axTitK0sDauEta[255];
    snprintf(axTitK0sDauEta,255,"K^{0}_{s} daughter;m_{inv} (GeV/c^{2});p_{T} (Gev/c);#eta_{pos};#eta_{neg}");
    char axTitLDauEta[255];
    snprintf(axTitLDauEta,255,"#Lambda daughter;m_{inv} (GeV/c^{2});p_{T} (GeV/c);#eta_{pos};#eta_{neg}");
    char axTitALDauEta[255];
    snprintf(axTitALDauEta,255,"#bar{#Lambda} daughter;m_{inv} (GeV/c^{2});p_{T} (GeV/c);#eta_{pos};#eta_{neg}");


    char axTitK0sDauPhi[255];
    snprintf(axTitK0sDauPhi,255,"K^{0}_{s} daughter;m_{inv} (GeV/c^{2});#phi_{pos};#phi_{neg};ITS hits (pos);ITS hits (neg)");
    char axTitLDauPhi[255];
    snprintf(axTitLDauPhi,255,"#Lambda daughter;m_{inv} (GeV/c^{2});#phi_{pos};#phi_{neg};ITS hits (pos);ITS hits (neg)");
    char axTitALDauPhi[255];
    snprintf(axTitALDauPhi,255,"#bar{#Lambda} daughter;m_{inv} (GeV/c^{2});#phi_{pos};#phi_{neg};ITS hits (pos);ITS hits (neg)");

  */
  // fTHnFK0s = new 	THnF(histTitK0s,axTitK0s,4,binsTHnV0K0s,minK0s,maxK0s);
  // // fTHnFK0s->Sumw2();
  // fTHnFL   = new 	THnF(histTitL  ,axTitL  ,4,binsTHnV0L,minL  ,maxL);
  // // fTHnFL->Sumw2();
  // fTHnFAL  = new 	THnF(histTitAL ,axTitAL ,4,binsTHnV0L,minL  ,maxL);
  //  fTHnFAL->Sumw2();
  /*

    fTHnFK0sDauEta = new 	THnF(histTitK0sDauEta,axTitK0sDauEta,4,binsTHnV0DauEtaK0s,minK0sDauEta,maxK0sDauEta);
    //  fTHnFK0sDauEta->Sumw2();
    fTHnFLDauEta   = new 	THnF(histTitLDauEta  ,axTitLDauEta  ,4,binsTHnV0DauEtaL,minLDauEta  ,maxLDauEta);
    //  fTHnFLDauEta->Sumw2();
    fTHnFALDauEta  = new 	THnF(histTitALDauEta ,axTitALDauEta ,4,binsTHnV0DauEtaL,minLDauEta  ,maxLDauEta);
    // fTHnFALDauEta->Sumw2();

    fTHnFK0sDauPhi = new 	THnF(histTitK0sDauPhi,axTitK0sDauPhi,5,binsTHnV0DauPhiK0s,minK0sDauPhi,maxK0sDauPhi);
    // fTHnFK0sDauPhi->Sumw2();
    fTHnFLDauPhi   = new 	THnF(histTitLDauPhi  ,axTitLDauPhi  ,5,binsTHnV0DauPhiL,minLDauPhi  ,maxLDauPhi);
    // fTHnFLDauPhi->Sumw2();
    fTHnFALDauPhi  = new 	THnF(histTitALDauPhi ,axTitALDauPhi ,5,binsTHnV0DauPhiL,minLDauPhi  ,maxLDauPhi); 
    //fTHnFALDauPhi->Sumw2();
    */
 
  fHistV0RadiusZ[0]  = new TH2F("fHistV0RadiusZ","z of decay radius vs 2D radius",100,0.0,100.0,250,-125.0,125.0);
  fHistV0RadiusZVSPt[0]  = new TH2F("fHistV0RadiusZVSPt","z of decay radius vs pt radius",200,0.0,20.0,125,0.0,125.0);
  fHistV0RadiusXY[0]  = new TH2F("fHistV0RadiusXY","y vs x decay radius",250,-125.0,125.0,250,-125.0,125.0);
  fHistV0RadiusXYVSY[0]  = new TH2F("fHistV0RadiusXYVSY","2D decay radius vs rap",100,-1,1,100,0.0,100.0);
  fHistArmenteros[0] = new TH2F("fHistArmenteros"," pi+pi- armenteros",nbMass,-1.,1.,500,0.,0.5);


  fHistPiPiK0sVsLambdaMass  = new TH2F("fHistPiPiK0sVsLambdaMass","K0s mass vs Lambda mass for all pt for K0s",250,1.05,1.25,250,0.25,0.75);
  fHistPiPiK0sVsALambdaMass = new TH2F("fHistPiPiK0sVsALambdaMass","K0s mass vs ALambda mass for all pt for K0s",250,1.05,1.25,250,0.25,0.75);

  fHistPiPK0sVsLambdaMass   = new TH2F("fHistPiPK0sVsLambdaMass","K0s mass vs Lambda mass for all pt for Lambda",250,1.05,1.25,250,0.25,0.75);

  fHistPiAPK0sVsALambdaMass = new TH2F("fHistPiAPK0sVsALambdaMass","K0s mass vs ALambda mass for all pt for ALambda",250,1.05,1.25,250,0.25,0.75);

  fHistPiPALambdaVsLambdaMass  = new TH2F("fHistPiPALambdaVsLambdaMass","ALambda mass vs Lambda mass for Lambda",250,1.05,1.25,250,1.05,1.25);
  fHistPiAPLambdaVsALambdaMass = new TH2F("fHistPiAPLambdaVsALambdaMass","Lambda mass vs ALambda mass for ALambda",250,1.05,1.25,250,1.05,1.25);

  //-----K0s---------//
  fHistPiPiMass = new TH1F("fHistPiPiMass"," pi+pi- InvMass distribution",2*nbMass,0.,2.);
  fHistPiPiMassVSPt = new TH2F("fHistPiPiMassVSPt","pi+pi- InvMass distribution",nbMass,0.25,0.75,200,0.0,20.0);
  fHistPiPiMassVSPtMCTruth = new TH2F("fHistPiPiMassVSPtMCTruth","pi+pi- InvMass distribution vs pt MCTruth",nbMass,0.25,0.75,200,0.0,20.0);
  fHistPiPiMassVSPtPosMCTruth = new TH2F("fHistPiPiMassVSPtPosMCTruth","pi+ InvMass distribution vs pt MCTruth",nbMass,0.25,0.75,200,0.0,20.0);
  fHistPiPiMassVSPtNegMCTruth = new TH2F("fHistPiPiMassVSPtNegMCTruth","pi- InvMass distribution vs pt MCTruth",nbMass,0.25,0.75,200,0.0,20.0);
  fHistPiPiMassVSY = new TH2F("fHistPiPiMassVSY","pi+pi- InvMass distribution vs rapidity",nbMass,0.25,0.75,200,-1.0,1.0);
  fHistPiPiPtVSY = new TH2F("fHistPiPiPtVSY","phi vs mass",100,-1,1,100,0.0,20);
  fHistPiPiDecayLengthVsPt = new TH2F("fHistPiPiDecayLengthVsPt","K0 decay length vs pt",200,0.0,20.0,220,0.0,110.0);
  fHistPiPiDecayLengthVsMass = new TH2F("fHistPiPiDecayLengthVsMass","K0s decay length vs mass",nbMass,0.25,0.75,220,0.0,110.0);  
  //  fHistPiPiPhiPosVsPtPosVsMass = new TH3F("fHistPiPiPhiPosVsPtPosVsMass","ctau K0s vs pt vs mass",250,0.25,0.75,120,0.0,60.0,200,0.0,20.0);//4.0);//xxx      
  if(!fSetPtDepHist){
    fHistPiPiDecayLengthVsCtau = new TH2F("fHistPiPiDecayLengthVsCtau","K0s ctau vs mass",nbMass,0.25,0.75,250,0.0,50.0);
  }
  else{
    fHistPiPiDecayLengthVsCtau = new TH2F("fHistPiPiDecayLengthVsCtau","K0s ctau vs pt",200,0,20.0,250,0.0,50.0);
  } 
  
  fHistPiPiMonitorCuts = new TH1F("fHistPiPiMonitorCuts","K0 cut monitor",35,0.5,35.5);
  fHistPiPiMonitorMCCuts = new TH1F("fHistPiPiMonitorMCCuts","K0 cut monitor mc",35,0.5,35.5);
  
  //---------------Lambda--------------//
  fHistPiPMass[0] = new TH1F("fHistPiPMass"," p+pi- InvMass distribution",2*nbMass,0.,2.);
  fHistPiPMassVSPt[0] = new TH2F("fHistPiPMassVSPt","p+pi- InvMass distribution",nbMass,1.05,1.25,200,0.0,20.0);
  fHistPiPMassVSPtMCTruth[0] = new TH2F("fHistPiPMassVSPtMCTruth","p+pi- InvMass distribution vs pt MCTruth",nbMass,1.05,1.25,200,0.0,20.0);
  fHistPiPMassVSPtPosMCTruth[0] = new TH2F("fHistPiPMassVSPtPosMCTruth","p+ InvMass distribution vs pt MCTruth",nbMass,1.05,1.25,200,0.0,20.0);
  fHistPiPMassVSPtNegMCTruth[0] = new TH2F("fHistPiPMassVSPtNegMCTruth","pi- InvMass distribution vs pt MCTruth",nbMass,1.05,1.25,200,0.0,20.0);
  fHistPiPMassVSY[0] = new TH2F("fHistPiPMassVSY","p+pi- InvMass distribution vs rapidity",nbMass,1.05,1.25,200,-1.0,1.0);
  fHistPiPPtVSY[0] = new TH2F("fHistPiPPtVSY","p{t} vs y",100,-1,1,100,0.0,20);
  fHistPiPDecayLengthVsPt[0] = new TH2F("fHistPiPDecayLengthVsPt","#Lambda decay length vs pt",200,0.0,20.0,220,0.0,110.0);
  fHistPiPDecayLengthVsMass[0] = new TH2F("fHistPiPDecayLengthVsMass","#Lambda decay length vs mass",nbMass,1.05,1.25,220,0.0,110.0);
  //  fHistPiPPhiPosVsPtPosVsMass  = new TH3F("fHistPiPPhiPosVsPtPosVsMass","ctau L vs pt vs mass",200,1.05,1.25,120,0.0,60.0,200,0.0,20.0);//4.0);//xxx        
  if(!fSetPtDepHist){
    fHistPiPDecayLengthVsCtau[0] = new TH2F("fHistPiPDecayLengthVsCtau","L ctau vs mass",nbMass,1.05,1.25,120,0.0,60.0);
  }
  else{
    fHistPiPDecayLengthVsCtau[0] = new TH2F("fHistPiPDecayLengthVsCtau","L ctau vs pt",200,0.0,20.0,120,0.0,60.0);
  }
  
  fHistPiPMonitorCuts[0] = new TH1F("fHistPiPMonitorCuts","#Lambda cut monitor",35,0.5,35.5);
  fHistPiPMonitorMCCuts[0] = new TH1F("fHistPiPMonitorMCCuts","#Lambda cut monitor mc ",35,0.5,35.5);
  
  //-------------ALamda-------------//
  fHistPiAPMass[0] = new TH1F("fHistPiAPMass"," ap-pi+ InvMass distribution",2*nbMass,0.,2.);
  fHistPiAPMassVSPt[0] = new TH2F("fHistPiAPMassVSPt","p-pi+ InvMass distribution",nbMass,1.05,1.25,200,0.0,20.0);
  fHistPiAPMassVSPtMCTruth[0] = new TH2F("fHistPiAPMassVSPtMCTruth","p-pi+ InvMass distribution vs pt MCTruth",nbMass,1.05,1.25,200,0.0,20.0);
  fHistPiAPMassVSPtPosMCTruth[0] = new TH2F("fHistPiAPMassVSPtPosMCTruth","pi- InvMass distribution vs pt MCTruth",nbMass,1.05,1.25,200,0.0,20.0);
  fHistPiAPMassVSPtNegMCTruth[0] = new TH2F("fHistPiAPMassVSPtNegMCTruth","p- InvMass distribution vs pt MCTruth",nbMass,1.05,1.25,200,0.0,20.0);
  fHistPiAPMassVSY[0] = new TH2F("fHistPiAPMassVSY","p-pi+ InvMass distribution vs rapidity",nbMass,1.05,1.25,200,-1.0,1.0);
  fHistPiAPPtVSY[0] = new TH2F("fHistPiAPPtVSY","p{t} vs y",100,-1,1,100,0.0,20);
  fHistPiAPDecayLengthVsPt[0] = new TH2F("fHistPiAPDecayLengthVsPt","#bar{#Lambda} decay length vs pt",200,0.0,20.0,220,0.0,110.0);
  fHistPiAPDecayLengthVsMass[0] = new TH2F("fHistPiAPDecayLengthVsMass","#bar{#Lambda} decay length vs mass",nbMass,1.05,1.25,220,0.0,110.0);
  //  if(fMCTruthMode) fHistPiAPPhiPosVsPtPosVsMass  = new TH3F("fHistPiAPPhiPosVsPtPosVsMass","ctau AL vs pt vs mass",200,1.05,1.25,120,0.0,60.0,200,0.0,20.0);//4.0);//xxx   
  if(!fSetPtDepHist){
    fHistPiAPDecayLengthVsCtau[0] = new TH2F("fHistPiAPDecayLengthVsCtau","AL ctau vs mass",nbMass,1.05,1.25,120,0.0,60.0);
  }
  else{
    fHistPiAPDecayLengthVsCtau[0] = new TH2F("fHistPiAPDecayLengthVsCtau","AL ctau vs pt",200,0.0,20.0,120,0.0,60.0);
  }
  
  fHistPiAPMonitorCuts[0] = new TH1F("fHistPiAPMonitorCuts","#bar{#Lambda} cut monitor",35,0.5,35.5);
  fHistPiAPMonitorMCCuts[0] = new TH1F("fHistPiAPMonitorMCCuts","#bar{#Lambda} cut monitor mc",35,0.5,35.5);   

  // ---------------------------------------------for MC reco secondaries -----------------------------------------//
  if(mchist==2){
    fHistV0RadiusZ[1]  = new TH2F("fHistV0RadiusZSec","z of decay radius vs 2D radius",100,0.0,100.0,250,-125.0,125.0);
    fHistV0RadiusZVSPt[1]  = new TH2F("fHistV0RadiusZVSPtSec","z of decay radius vs pt radius",200,0.0,20.0,125,0.0,125.0);
    fHistV0RadiusXY[1]  = new TH2F("fHistV0RadiusXYSec","y vs x decay radius",250,-125.0,125.0,250,-125.0,125.0);
    fHistV0RadiusXYVSY[1]  = new TH2F("fHistV0RadiusXYVSYSec","2D decay radius vs rap",100,-1,1,100,0.0,100.0);
    fHistArmenteros[1] = new TH2F("fHistArmenterosSec"," pi+pi- armenteros",nbMass,-1.,1.,500,0.,0.5);

    //-----------------K0s------------//
    //--------------- Lambda----------//
    fHistPiPMass[1] = new TH1F("fHistPiPMassSec"," p+pi- InvMass distribution",2*nbMass,0.,2.);
    fHistPiPMassVSPt[1] = new TH2F("fHistPiPMassVSPtSec","p+pi- InvMass distribution",nbMass,1.05,1.25,200,0.0,20.0);
    fHistPiPMassVSPtMCTruth[1] = new TH2F("fHistPiPMassVSPtMCTruthSec","p+pi- InvMass distribution vs pt MCTruth",nbMass,1.05,1.25,200,0.0,20.0);
    fHistPiPMassVSPtPosMCTruth[1] = new TH2F("fHistPiPMassVSPtPosMCTruthSec","p+ InvMass distribution vs pt MCTruth",nbMass,1.05,1.25,200,0.0,20.0);
    fHistPiPMassVSPtNegMCTruth[1] = new TH2F("fHistPiPMassVSPtNegMCTruthSec","pi- InvMass distribution vs pt MCTruth",nbMass,1.05,1.25,200,0.0,20.0);
    fHistPiPMassVSY[1] = new TH2F("fHistPiPMassVSYSec","p+pi- InvMass distribution vs rapidity",nbMass,1.05,1.25,200,-1.0,1.0);   
    fHistPiPPtVSY[1] = new TH2F("fHistPiPPtVSYSec","p{t} vs y",100,-1,1,100,0.0,20);
    fHistPiPDecayLengthVsPt[1] = new TH2F("fHistPiPDecayLengthVsPtSec","#Lambda decay length vs pt",200,0.0,20.0,220,0.0,110.0);
    fHistPiPDecayLengthVsMass[1] = new TH2F("fHistPiPDecayLengthVsMassSec","#Lambda decay length vs mass",nbMass,1.05,1.25,220,0.0,110.0);
    if(!fSetPtDepHist){
      fHistPiPDecayLengthVsCtau[1] = new TH2F("fHistPiPDecayLengthVsCtauSec","L ctau vs mass",nbMass,1.05,1.25,120,0.0,60.0);
    }
    else{
      fHistPiPDecayLengthVsCtau[1] = new TH2F("fHistPiPDecayLengthVsCtauSec","L ctau vs pt",200,0.0,20.0,120,0.0,60.0);
    }
    
    fHistPiPMonitorCuts[1] = new TH1F("fHistPiPMonitorCutsSec","#Lambda cut monitor",35,0.5,35.5);
    fHistPiPMonitorMCCuts[1] = new TH1F("fHistPiPMonitorMCCutsSec","#Lambda cut monitor mc",35,0.5,35.5);

    //----------------ALambda---------//
    fHistPiAPMass[1] = new TH1F("fHistPiAPMassSec"," ap-pi+ InvMass distribution",2*nbMass,0.,2.);
    fHistPiAPMassVSPt[1] = new TH2F("fHistPiAPMassVSPtSec","p-pi+ InvMass distribution",nbMass,1.05,1.25,200,0.0,20.0);
    fHistPiAPMassVSPtMCTruth[1] = new TH2F("fHistPiAPMassVSPtMCTruthSec","p-pi+ InvMass distribution vs pt MCTruth",nbMass,1.05,1.25,200,0.0,20.0);
    fHistPiAPMassVSPtPosMCTruth[1] = new TH2F("fHistPiAPMassVSPtPosMCTruthSec","p+ InvMass distribution vs pt MCTruth",nbMass,1.05,1.25,200,0.0,20.0);
    fHistPiAPMassVSPtNegMCTruth[1] = new TH2F("fHistPiAPMassVSPtNegMCTruthSec","pi- InvMass distribution vs pt MCTruth",nbMass,1.05,1.25,200,0.0,20.0);
    fHistPiAPMassVSY[1] = new TH2F("fHistPiAPMassVSYSec","p-pi+ InvMass distribution vs rapidity",nbMass,1.05,1.25,200,-1.0,1.0);
    fHistPiAPPtVSY[1] = new TH2F("fHistPiAPPtVSYSec","p{t} vs y",100,-1,1,100,0.0,20);
    fHistPiAPDecayLengthVsPt[1] = new TH2F("fHistPiAPDecayLengthVsPtSec","#bar{#Lambda} decay length vs pt",200,0.0,20.0,220,0.0,110.0);
    fHistPiAPDecayLengthVsMass[1] = new TH2F("fHistPiAPDecayLengthVsMassSec","#bar{#Lambda} decay length vs mass",nbMass,1.05,1.25,220,0.0,110.0);
    if(!fSetPtDepHist){
      fHistPiAPDecayLengthVsCtau[1] = new TH2F("fHistPiAPDecayLengthVsCtauSec","AL ctau vs mass",nbMass,1.05,1.25,120,0.0,60.0);
    }
    else{
      fHistPiAPDecayLengthVsCtau[1] = new TH2F("fHistPiAPDecayLengthVsCtauSec","AL ctau vs pt",200,0.0,20.0,120,0.0,60.0);
    }

    fHistPiAPMonitorCuts[1] = new TH1F("fHistPiAPMonitorCutsSec","#bar{#Lambda} cut monitor",35,0.5,35.5);
    fHistPiAPMonitorMCCuts[1] = new TH1F("fHistPiAPMonitorMCCutsSec","#bar{#Lambda} cut monitor mc",35,0.5,35.5);
  }

  //add to output container
  //------------ K0s ------------------//
  fOutputContainer->Add(fHistPiPiMass);	 
  fOutputContainer->Add(fHistPiPiMassVSPt);
  fOutputContainer->Add(fHistPiPiMassVSPtMCTruth);
  fOutputContainer->Add(fHistPiPiMassVSPtPosMCTruth);
  fOutputContainer->Add(fHistPiPiMassVSPtNegMCTruth);
  fOutputContainer->Add(fHistPiPiMassVSY);
  fOutputContainer->Add(fHistPiPiPtVSY);
  fOutputContainer->Add(fHistPiPiDecayLengthVsPt);
  fOutputContainer->Add(fHistPiPiDecayLengthVsCtau);
  fOutputContainer->Add(fHistPiPiDecayLengthVsMass);  
  fOutputContainer->Add(fHistPiPiMonitorCuts);
  fOutputContainer->Add(fHistPiPiMonitorMCCuts);
  fOutputContainer->Add(fHistPiPiK0sVsLambdaMass);
  fOutputContainer->Add(fHistPiPiK0sVsALambdaMass);

  //  fOutputContainer->Add(fTHnFK0s);
  /*
    fOutputContainer->Add(fTHnFK0sDauEta);
    fOutputContainer->Add(fTHnFK0sDauPhi);
    //fOutputContainer->Add(fHistPiPiPhiPosVsPtPosVsMass);//xxx      
    */
  // --------------- Lambda ---------------//
  fOutputContainer->Add(fHistPiPK0sVsLambdaMass);
  fOutputContainer->Add(fHistPiPALambdaVsLambdaMass);
  //  fOutputContainer->Add(fHistPiPPhiPosVsPtPosVsMass);//xxx
  //  fOutputContainer->Add(fTHnFL);
  /*
    fOutputContainer->Add(fTHnFLDauEta);
    fOutputContainer->Add(fTHnFLDauPhi);
  */
  // --------------- ALambda ---------------//
  fOutputContainer->Add(fHistPiAPK0sVsALambdaMass);
  fOutputContainer->Add(fHistPiAPLambdaVsALambdaMass);
  //  if(fMCTruthMode)  fOutputContainer->Add(fHistPiAPPhiPosVsPtPosVsMass);//xxx

  // fOutputContainer->Add(fTHnFAL);
  /*
    fOutputContainer->Add(fTHnFALDauEta);
    fOutputContainer->Add(fTHnFALDauPhi);
  */
  
  for(Int_t j=0;j<mchist;j++){
    fOutputContainer->Add(fHistArmenteros[j]);
    fOutputContainer->Add(fHistV0RadiusZ[j]);
    fOutputContainer->Add(fHistV0RadiusZVSPt[j]);
    fOutputContainer->Add(fHistV0RadiusXY[j]);
    fOutputContainer->Add(fHistV0RadiusXYVSY[j]);
    fOutputContainer->Add(fHistPiPMass[j]);
    fOutputContainer->Add(fHistPiAPMass[j]);
    fOutputContainer->Add(fHistPiPMassVSPt[j]);
    fOutputContainer->Add(fHistPiAPMassVSPt[j]);
    fOutputContainer->Add(fHistPiPMassVSPtMCTruth[j]);
    fOutputContainer->Add(fHistPiPMassVSPtPosMCTruth[j]);
    fOutputContainer->Add(fHistPiPMassVSPtNegMCTruth[j]);
    fOutputContainer->Add(fHistPiAPMassVSPtMCTruth[j]);
    fOutputContainer->Add(fHistPiAPMassVSPtPosMCTruth[j]);
    fOutputContainer->Add(fHistPiAPMassVSPtNegMCTruth[j]);
    fOutputContainer->Add(fHistPiPMassVSY[j]);
    fOutputContainer->Add(fHistPiAPMassVSY[j]);      
    fOutputContainer->Add(fHistPiPPtVSY[j]);
    fOutputContainer->Add(fHistPiAPPtVSY[j]);
    fOutputContainer->Add(fHistPiPDecayLengthVsPt[j]);
    fOutputContainer->Add(fHistPiAPDecayLengthVsPt[j]);
    fOutputContainer->Add(fHistPiPDecayLengthVsCtau[j]);
    fOutputContainer->Add(fHistPiAPDecayLengthVsCtau[j]);
    fOutputContainer->Add(fHistPiPDecayLengthVsMass[j]);
    fOutputContainer->Add(fHistPiAPDecayLengthVsMass[j]);
    fOutputContainer->Add(fHistPiPMonitorCuts[j]);
    fOutputContainer->Add(fHistPiAPMonitorCuts[j]);
    fOutputContainer->Add(fHistPiPMonitorMCCuts[j]);
    fOutputContainer->Add(fHistPiAPMonitorMCCuts[j]);
  }
  
  //----------------- for reco or data or mc data like MC reco only -----------------//
  if((fMCMode) || (!fMCTruthMode && !fMCMode)){
    
    fHistPiPiEtaDReco[0] = new TH2F("fHistPiPiEtaDRecoRaw","K0s daughters eta raw",300,-6,6,100,0,20);
    fOutputContainer->Add(fHistPiPiEtaDReco[0]);
    fHistPiPiEtaDReco[1] = new TH2F("fHistPiPiEtaDReco","K0s daughters eta after rap V0 cut pos",300,-3,3,300,-3.00,3.0);
    fOutputContainer->Add(fHistPiPiEtaDReco[1]);	 
    fHistPiPEtaDReco[0] = new TH2F("fHistPiPEtaDRecoRaw","#Lambda daughters eta raw",300,-6,6,100,0,20);
    fOutputContainer->Add(fHistPiPEtaDReco[0]);
    fHistPiPEtaDReco[1] = new TH2F("fHistPiPEtaDReco","#Lambda daughters eta after rap V0 cut neg",300,-3,3,300,-3.00,3.0);
    fOutputContainer->Add(fHistPiPEtaDReco[1]);
	
    //-------------K0---------------//
    /*
    // fHistPiPiMassVSAlpha = new TH2F("fHistPiPiMassVSAlpha"," alpha armenteros vs pi+pi- InvMass distribution",nbMass,0.25,0.75,500,-1.,1.);
    fHistPiPiDistDaughtersPos[0]= new TH2F("fHistPiPiDistDaughtersPosPt0","K0s pos daughters distance to other tracks vs mass for pt 0-2GeV/c",250,0.25,0.75,100,0.0,20.0);
    fHistPiPiDistDaughtersNeg[0]= new TH2F("fHistPiPiDistDaughtersNegPt0","K0s neg daughters distance to other tracks vs mass for pt 0-2GeV/c",250,0.25,0.75,100,0.0,20.0);
    fHistPiPiDCADaughtersPos[0]= new TH2F("fHistPiPiDCADaughtersPosPt0","K0s pos daughters DCA to other tracks vs mass for pt 0-2GeV/c",250,0.25,0.75,100,0.0,20.0);
    fHistPiPiDCADaughtersNeg[0]= new TH2F("fHistPiPiDCADaughtersNegPt0","K0s neg daughters DCA to other tracks vs mass for pt 0-2GeV/c",250,0.25,0.75,100,0.0,20.0);
    fHistPiPiRadAtDCA5cmDaughtersPos[0]= new TH2F("fHistPiPiRadAtDCA5cmDaughtersPosPt0","K0s pos daughter position for DCA < 5cm vs mass for pt 0-2GeV/c",250,0.25,0.75,160,0.0,240.0);//xxx
    fHistPiPiRadAtDCA5cmDaughtersNeg[0]= new TH2F("fHistPiPiRadAtDCA5cmDaughtersNegPt0","K0s pos daughter position for DCA < 5cm vs mass for pt 0-2GeV/c",250,0.25,0.75,160,0.0,240.0);

    fHistPiPiDistDaughtersPos[1]= new TH2F("fHistPiPiDistDaughtersPosPt1","K0s pos daughters distance to other tracks vs mass for pt 2-6GeV/c",250,0.25,0.75,100,0.0,20.0);
    fHistPiPiDistDaughtersNeg[1]= new TH2F("fHistPiPiDistDaughtersNegPt1","K0s neg daughters distance to other tracks vs mass for pt 2-6GeV/c",250,0.25,0.75,100,0.0,20.0);
    fHistPiPiDCADaughtersPos[1]= new TH2F("fHistPiPiDCADaughtersPosPt1","K0s pos daughters DCA to other tracks vs mass for pt 2-6GeV/c",250,0.25,0.75,100,0.0,20.0);
    fHistPiPiDCADaughtersNeg[1]= new TH2F("fHistPiPiDCADaughtersNegPt1","K0s neg daughters DCA to other tracks vs mass for pt 2-6GeV/c",250,0.25,0.75,100,0.0,20.0);
    fHistPiPiRadAtDCA5cmDaughtersPos[1]= new TH2F("fHistPiPiRadAtDCA5cmDaughtersPosPt1","K0s pos daughter position for DCA < 5cm vs mass for pt 2-6GeV/c",250,0.25,0.75,160,0.0,240.0);
    fHistPiPiRadAtDCA5cmDaughtersNeg[1]= new TH2F("fHistPiPiRadAtDCA5cmDaughtersNegPt1","K0s pos daughter position for DCA < 5cm vs mass for pt 2-6GeV/c",250,0.25,0.75,160,0.0,240.0);

    fHistPiPiDistDaughtersPos[2]= new TH2F("fHistPiPiDistDaughtersPosPt2","K0s pos daughters distance to other tracks vs mass for pt > 6 GeV/c",250,0.25,0.75,100,0.0,20.0);
    fHistPiPiDistDaughtersNeg[2]= new TH2F("fHistPiPiDistDaughtersNegPt2","K0s neg daughters distance to other tracks vs mass for pt > 6 GeV/c",250,0.25,0.75,100,0.0,20.0);
    fHistPiPiDCADaughtersPos[2]= new TH2F("fHistPiPiDCADaughtersPosPt2","K0s pos daughters DCA to other tracks vs mass for pt > 6 GeV/c",250,0.25,0.75,100,0.0,20.0);
    fHistPiPiDCADaughtersNeg[2]= new TH2F("fHistPiPiDCADaughtersNegPt2","K0s neg daughters DCA to other tracks vs mass for pt > 6 GeV/c",250,0.25,0.75,100,0.0,20.0);
    fHistPiPiRadAtDCA5cmDaughtersPos[2]= new TH2F("fHistPiPiRadAtDCA5cmDaughtersPosPt2","K0s pos daughter position for DCA < 5cm vs mass for pt > 6 GeV/c",250,0.25,0.75,160,0.0,240.0);
    fHistPiPiRadAtDCA5cmDaughtersNeg[2]= new TH2F("fHistPiPiRadAtDCA5cmDaughtersNegPt2","K0s pos daughter position for DCA < 5cm vs mass for pt > 6 GeV/c",250,0.25,0.75,160,0.0,240.0);
    */   
    // fHistPiPiDistDaughtersTPCEntrVsMass = new TH2F("fHistPiPiDistDaughtersTPCEntrVsMass","K0s distance of daughters at TPC entrance vs mass",nbMass,0.25,0.75,100,0.0,20.0);
   
    if(!fSetPtDepHist){
      fHistPiPiDCADaughters = new TH2F("fHistPiPiDCADaughters","dca of K0 daughters",nbMass,0.25,0.75,250,0.0,2);
      fHistPiPiDCADaughterPosToPrimVtxVSMass = new TH2F("fHistPiPiDCADaughterPosToPrimVtxVSMass","pi+ DCA daughter to prim vtx vsinvmass",nbMass,0.25,0.75,250,0.0,10.0);
      fHistPiPiDCAVSMass = new TH2F("fHistPiPiDCAVSMass","pi+pi- dca  vs pt",nbMass,0.25,0.75,250,0.0,5.0);
      fHistPiPiDCAZVSMass = new TH2F("fHistPiPiDCAZVSMass","pi+pi- dca z vs pt",nbMass,0.25,0.75,200,-20.0,20.0);
      fHistPiPiCosPointAng  = new TH2F("fHistPiPiCosPointAng","K0 cosine of pointing angle vs mass ",nbMass,0.25,0.75,200,0.99,1.00);
      fHistPiPiRadiusXY = new TH2F("fHistPiPiRadiusXY","pi+pi- phi dist vs mass",nbMass,0.25,0.75,200,0.0,4.0);
      // fHistPiPiPtDaughters = new TH2F("fHistPiPiPtDaughters","p_{t} pos vs p_{t} neg of daughters",400,0.0,20.0,400,0,20.0);
      fHistPiPiDCAZPos = new TH2F("fHistPiPiDCAZPos","dca z  of K0 pos daughters",nbMass,0.25,0.75,200,-20.0,20.0);
      fHistPiPiDCAZNeg = new TH2F("fHistPiPiDCAZNeg","dca z  of K0 neg daughters",nbMass,0.25,0.75,200,-20.0,20.0);
      fHistPiPiTrackLengthPosVsMass = new TH2F("fHistPiPiTrackLengthPosVsMass","track lenght of pos K0s daughter in TPC",nbMass,0.25,0.75,250,0.0,250.0);
      fHistPiPiTrackLengthNegVsMass = new TH2F("fHistPiPiTrackLengthNegVsMass","track lenght of neg K0s daughter in TPC",nbMass,0.25,0.75,250,0.0,250.0);
    }
    else{//pt dependence
      fHistPiPiDCADaughters = new TH2F("fHistPiPiDCADaughters","dca of K0 daughters",200,0.0,20.0,250,0.0,2);
      fHistPiPiDCADaughterPosToPrimVtxVSMass = new TH2F("fHistPiPiDCADaughterPosToPrimVtxVSMass","pi+ DCA daughter to prim vtx vsinvmass",200,0.0,20.0,250,0.0,10.0);
      fHistPiPiDCAVSMass = new TH2F("fHistPiPiDCAVSMass","pi+pi- dca  vs pt",200,0.0,20.0,250,0.0,5.0);
      fHistPiPiDCAZVSMass = new TH2F("fHistPiPiDCAZVSMass","pi+pi- dca z vs pt",200,0.0,20.0,200,-20.0,20.0);
      fHistPiPiCosPointAng  = new TH2F("fHistPiPiCosPointAng","K0 cosine of pointing angle vs mass ",200,0.0,20.0,200,0.99,1.00);
      fHistPiPiRadiusXY = new TH2F("fHistPiPiRadiusXY","pi+pi- phi dist vs mass",200,0.0,20.0,200,0.0,4.0);
      fHistPiPiDCAZPos = new TH2F("fHistPiPiDCAZPos","dca z  of K0 pos daughters",200,0.0,20.0,200,-20.0,20.0);
      fHistPiPiDCAZNeg = new TH2F("fHistPiPiDCAZNeg","dca z  of K0 neg daughters",200,0.0,20.0,200,-20.0,20.0);
      fHistPiPiTrackLengthPosVsMass = new TH2F("fHistPiPiTrackLengthPosVsMass","track lenght of pos K0s daughter in TPC",200,0.0,20.0,250,0.0,250.0);
      fHistPiPiTrackLengthNegVsMass = new TH2F("fHistPiPiTrackLengthNegVsMass","track lenght of neg K0s daughter in TPC",200,0.0,20.0,250,0.0,250.0);
    }

    //---------------Lambda-------------//
    /*
      fHistPiPDistDaughtersPos[0]= new TH2F("fHistPiPDistDaughtersPosPt0","Lambda pos daughters distance to other tracks vs mass for pt 0-2GeV/c",250,1.05,1.25,100,0.0,20.0);
      fHistPiPDistDaughtersNeg[0]= new TH2F("fHistPiPDistDaughtersNegPt0","Lambda neg daughters distance to other tracks vs mass for pt 0-2GeV/c",250,1.05,1.25,100,0.0,20.0);
      fHistPiPDCADaughtersPos[0]= new TH2F("fHistPiPDCADaughtersPosPt0","Lambda pos daughters DCA to other tracks vs mass for pt 0-2GeV/c",250,1.05,1.25,100,0.0,20.0);
      fHistPiPDCADaughtersNeg[0]= new TH2F("fHistPiPDCADaughtersNegPt0","Lambda neg daughters DCA to other tracks vs mass for pt 0-2GeV/c",250,1.05,1.25,100,0.0,20.0);
      fHistPiPRadAtDCA5cmDaughtersPos[0]= new TH2F("fHistPiPRadAtDCA5cmDaughtersPosPt0","Lambda pos daughter position for DCA < 5cm vs mass for pt 0-2GeV/c",250,1.05,1.25,160,0.0,240.0);
      fHistPiPRadAtDCA5cmDaughtersNeg[0]= new TH2F("fHistPiPRadAtDCA5cmDaughtersNegPt0","Lambda pos daughter position for DCA < 5cm vs mass for pt 0-2GeV/c",250,1.05,1.25,160,0.0,240.0);

      fHistPiPDistDaughtersPos[1]= new TH2F("fHistPiPDistDaughtersPosPt1","Lambda pos daughters distance to other tracks vs mass for pt 2-6GeV/c",250,1.05,1.25,100,0.0,20.0);
      fHistPiPDistDaughtersNeg[1]= new TH2F("fHistPiPDistDaughtersNegPt1","Lambda neg daughters distance to other tracks vs mass for pt 2-6GeV/c",250,1.05,1.25,100,0.0,20.0);
      fHistPiPDCADaughtersPos[1]= new TH2F("fHistPiPDCADaughtersPosPt1","Lambda pos daughters DCA to other tracks vs mass for pt 2-6GeV/c",250,1.05,1.25,100,0.0,20.0);
      fHistPiPDCADaughtersNeg[1]= new TH2F("fHistPiPDCADaughtersNegPt1","Lambda neg daughters DCA to other tracks vs mass for pt 2-6GeV/c",250,1.05,1.25,100,0.0,20.0);
      fHistPiPRadAtDCA5cmDaughtersPos[1]= new TH2F("fHistPiPRadAtDCA5cmDaughtersPosPt1","Lambda pos daughter position for DCA < 5cm vs mass for pt 2-6GeV/c",250,1.05,1.25,160,0.0,240.0);
      fHistPiPRadAtDCA5cmDaughtersNeg[1]= new TH2F("fHistPiPRadAtDCA5cmDaughtersNegPt1","Lambda pos daughter position for DCA < 5cm vs mass for pt 2-6GeV/c",250,1.05,1.25,160,0.0,240.0);

      fHistPiPDistDaughtersPos[2]= new TH2F("fHistPiPDistDaughtersPosPt2","Lambda pos daughters distance to other tracks vs mass for pt > 6 GeV/c",250,1.05,1.25,100,0.0,20.0);
      fHistPiPDistDaughtersNeg[2]= new TH2F("fHistPiPDistDaughtersNegPt2","Lambda neg daughters distance to other tracks vs mass for pt > 6 GeV/c",250,1.05,1.25,100,0.0,20.0);
      fHistPiPDCADaughtersPos[2]= new TH2F("fHistPiPDCADaughtersPosPt2","Lambda pos daughters DCA to other tracks vs mass for pt > 6 GeV/c",250,1.05,1.25,100,0.0,20.0);
      fHistPiPDCADaughtersNeg[2]= new TH2F("fHistPiPDCADaughtersNegPt2","Lambda neg daughters DCA to other tracks vs mass for pt > 6 GeV/c",250,1.05,1.25,100,0.0,20.0);
      fHistPiPRadAtDCA5cmDaughtersPos[2]= new TH2F("fHistPiPRadAtDCA5cmDaughtersPosPt2","Lambda pos daughter position for DCA < 5cm vs mass for pt > 6 GeV/c",250,1.05,1.25,160,0.0,240.0);
      fHistPiPRadAtDCA5cmDaughtersNeg[2]= new TH2F("fHistPiPRadAtDCA5cmDaughtersNegPt2","Lambda pos daughter position for DCA < 5cm vs mass for pt > 6 GeV/c",250,1.05,1.25,160,0.0,240.0);
    */
    // fHistPiPDistDaughtersTPCEntrVsMass = new TH2F("fHistPiPDistDaughtersTPCEntrVsMass","Lambda distance of daughters at TPC entrance vs mass",nbMass,1.05,1.25,100,0.0,20.0);

    if(!fSetPtDepHist){
      fHistPiPDCADaughters[0] = new TH2F("fHistPiPDCADaughters","dca of #Lambda daughters",nbMass,1.05,1.25,250,0.0,2.0);
      fHistPiPDCADaughterPosToPrimVtxVSMass[0] = new TH2F("fHistPiPDCADaughterPosToPrimVtxVSMass","pos DCA daughter to prim vtx vs invmass",nbMass,1.05,1.25,250,0.0,10.0);
      fHistPiPDCADaughterNegToPrimVtxVSMass[0] = new TH2F("fHistPiPDCADaughterNegToPrimVtxVSMass","neg DCA daughter to prim vtx vs invmass",nbMass,1.05,1.25,250,0.0,10.0);
      fHistPiPDCAVSMass[0] = new TH2F("fHistPiPDCAVSMass","ppi- dca  vs pt",nbMass,1.05,1.25,250,0.0,5.0);
      fHistPiPDCAZVSMass[0] = new TH2F("fHistPiPDCAZVSMass","ppi- dca z vs pt",nbMass,1.05,1.25,200,-20.0,20.0);
      fHistPiPCosPointAng[0]  = new TH2F("fHistPiPCosPointAng","#Lambda cosine of pointing angle vs mass ",nbMass,1.05,1.25,200,0.99,1.00);
      fHistPiPRadiusXY[0] = new TH2F("fHistPiPRadiusXY","pi-p+ phi dist vs mass",nbMass,1.05,1.25,200,0.0,4.0);
      // fHistPiPPtDaughters[0] = new TH2F("fHistPiPPtDaughters","p_{t} pos vs p_{t} neg of daughters",400,0.0,20.0,400,0,20.0);
      fHistPiPDCAZPos[0] = new TH2F("fHistPiPDCAZPos","dca z  of Lambda pos daughters",nbMass,1.05,1.25,200,-20.0,20.0);
      fHistPiPDCAZNeg[0] = new TH2F("fHistPiPDCAZNeg","dca z  of Lambda neg daughters",nbMass,1.05,1.25,200,-20.0,20.0);
      fHistPiPTrackLengthPosVsMass[0] = new TH2F("fHistPiPTrackLengthPosVsMass","track length of pos Lambda daughter in TPC",nbMass,1.05,1.25,250,0.0,250.0);
      fHistPiPTrackLengthNegVsMass[0] = new TH2F("fHistPiPTrackLengthNegVsMass","track length of neg Lambda daughter in TPC",nbMass,1.05,1.25,250,0.0,250.0);
    }
    else{//pt dependence
      fHistPiPDCADaughters[0] = new TH2F("fHistPiPDCADaughters","dca of #Lambda daughters",200,0.0,20.0,250,0.0,2.0);
      fHistPiPDCADaughterPosToPrimVtxVSMass[0] = new TH2F("fHistPiPDCADaughterPosToPrimVtxVSMass","pos DCA daughter to prim vtx vs invmass",200,0.0,20.0,250,0.0,10.0);
      fHistPiPDCADaughterNegToPrimVtxVSMass[0] = new TH2F("fHistPiPDCADaughterNegToPrimVtxVSMass","neg DCA daughter to prim vtx vs invmass",200,0.0,20.0,250,0.0,10.0);
      fHistPiPDCAVSMass[0] = new TH2F("fHistPiPDCAVSMass","ppi- dca  vs pt",200,0.0,20.0,250,0.0,5.0);
      fHistPiPDCAZVSMass[0] = new TH2F("fHistPiPDCAZVSMass","ppi- dca z vs pt",200,0.0,20.0,200,-20.0,20.0);
      fHistPiPCosPointAng[0]  = new TH2F("fHistPiPCosPointAng","#Lambda cosine of pointing angle vs mass ",200,0.0,20.0,200,0.99,1.00);
      fHistPiPRadiusXY[0] = new TH2F("fHistPiPRadiusXY","pi-p+ phi dist vs mass",200,0.0,20.0,200,0.0,4.0);
      fHistPiPDCAZPos[0] = new TH2F("fHistPiPDCAZPos","dca z  of Lambda pos daughters",200,0.0,20.0,200,-20.0,20.0);
      fHistPiPDCAZNeg[0] = new TH2F("fHistPiPDCAZNeg","dca z  of Lambda neg daughters",200,0.0,20.0,200,-20.0,20.0);
      fHistPiPTrackLengthPosVsMass[0] = new TH2F("fHistPiPTrackLengthPosVsMass","track length of pos Lambda daughter in TPC",200,0.0,20.0,250,0.0,250.0);
      fHistPiPTrackLengthNegVsMass[0] = new TH2F("fHistPiPTrackLengthNegVsMass","track length of neg Lambda daughter in TPC",200,0.0,20.0,250,0.0,250.0);
    }

    //-------------------AntiLambda-------------//
    /*
      fHistPiAPDistDaughtersPos[0]= new TH2F("fHistPiAPDistDaughtersPosPt0","ALambda pos daughters distance to other tracks vs mass for pt 0-2GeV/c",250,1.05,1.25,100,0.0,20.0);
      fHistPiAPDistDaughtersNeg[0]= new TH2F("fHistPiAPDistDaughtersNegPt0","ALambda neg daughters distance to other tracks vs mass for pt 0-2GeV/c",250,1.05,1.25,100,0.0,20.0);
      fHistPiAPDCADaughtersPos[0]= new TH2F("fHistPiAPDCADaughtersPosPt0","ALambda pos daughters DCA to other tracks vs mass for pt 0-2GeV/c",250,1.05,1.25,100,0.0,20.0);
      fHistPiAPDCADaughtersNeg[0]= new TH2F("fHistPiAPDCADaughtersNegPt0","ALambda neg daughters DCA to other tracks vs mass for pt 0-2GeV/c",250,1.05,1.25,100,0.0,20.0);
      fHistPiAPRadAtDCA5cmDaughtersPos[0]= new TH2F("fHistPiAPRadAtDCA5cmDaughtersPosPt0","ALambda pos daughter position for DCA < 5cm vs mass for pt 0-2GeV/c",250,1.05,1.25,160,0.0,240.0);
      fHistPiAPRadAtDCA5cmDaughtersNeg[0]= new TH2F("fHistPiAPRadAtDCA5cmDaughtersNegPt0","ALambda pos daughter position for DCA < 5cm vs mass for pt 0-2GeV/c",250,1.05,1.25,160,0.0,240.0);

      fHistPiAPDistDaughtersPos[1]= new TH2F("fHistPiAPDistDaughtersPosPt1","ALambda pos daughters distance to other tracks vs mass for pt 2-6GeV/c",250,1.05,1.25,100,0.0,20.0);
      fHistPiAPDistDaughtersNeg[1]= new TH2F("fHistPiAPDistDaughtersNegPt1","ALambda neg daughters distance to other tracks vs mass for pt 2-6GeV/c",250,1.05,1.25,100,0.0,20.0);
      fHistPiAPDCADaughtersPos[1]= new TH2F("fHistPiAPDCADaughtersPosPt1","ALambda pos daughters DCA to other tracks vs mass for pt 2-6GeV/c",250,1.05,1.25,100,0.0,20.0);
      fHistPiAPDCADaughtersNeg[1]= new TH2F("fHistPiAPDCADaughtersNegPt1","ALambda neg daughters DCA to other tracks vs mass for pt 2-6GeV/c",250,1.05,1.25,100,0.0,20.0);
      fHistPiAPRadAtDCA5cmDaughtersPos[1]= new TH2F("fHistPiAPRadAtDCA5cmDaughtersPosPt1","ALambda pos daughter position for DCA < 5cm vs mass for pt 2-6GeV/c",250,1.05,1.25,160,0.0,240.0);
      fHistPiAPRadAtDCA5cmDaughtersNeg[1]= new TH2F("fHistPiAPRadAtDCA5cmDaughtersNegPt1","ALambda pos daughter position for DCA < 5cm vs mass for pt 2-6GeV/c",250,1.05,1.25,160,0.0,240.0);

      fHistPiAPDistDaughtersPos[2]= new TH2F("fHistPiAPDistDaughtersPosPt2","ALambda pos daughters distance to other tracks vs mass for pt > 6 GeV/c",250,1.05,1.25,100,0.0,20.0);
      fHistPiAPDistDaughtersNeg[2]= new TH2F("fHistPiAPDistDaughtersNegPt2","ALambda neg daughters distance to other tracks vs mass for pt > 6 GeV/c",250,1.05,1.25,100,0.0,20.0);
      fHistPiAPDCADaughtersPos[2]= new TH2F("fHistPiAPDCADaughtersPosPt2","ALambda pos daughters DCA to other tracks vs mass for pt > 6 GeV/c",250,1.05,1.25,100,0.0,20.0);
      fHistPiAPDCADaughtersNeg[2]= new TH2F("fHistPiAPDCADaughtersNegPt2","ALambda neg daughters DCA to other tracks vs mass for pt > 6 GeV/c",250,1.05,1.25,100,0.0,20.0);
      fHistPiAPRadAtDCA5cmDaughtersPos[2]= new TH2F("fHistPiAPRadAtDCA5cmDaughtersPosPt2","ALambda pos daughter position for DCA < 5cm vs mass for pt > 6 GeV/c",250,1.05,1.25,160,0.0,240.0);
      fHistPiAPRadAtDCA5cmDaughtersNeg[2]= new TH2F("fHistPiAPRadAtDCA5cmDaughtersNegPt2","ALambda pos daughter position for DCA < 5cm vs mass for pt > 6 GeV/c",250,1.05,1.25,160,0.0,240.0);
    */
    // fHistPiAPDistDaughtersTPCEntrVsMass = new TH2F("fHistPiAPDistDaughtersTPCEntrVsMass","ALambda distance of daughters at TPC entrance vs mass",nbMass,1.05,1.25,100,0.0,20.0);

    if(!fSetPtDepHist){
      fHistPiAPDCADaughters[0] = new TH2F("fHistPiAPDCADaughters","dca of #bar{#Lambda} daughters",nbMass,1.05,1.25,250,0.0,2.0);
      fHistPiAPDCADaughterPosToPrimVtxVSMass[0] = new TH2F("fHistPiAPDCADaughterPosToPrimVtxVSMass","pos DCA daughter to prim vtx vs invmass",nbMass,1.05,1.25,250,0.0,10.0);
      fHistPiAPDCADaughterNegToPrimVtxVSMass[0] = new TH2F("fHistPiAPDCADaughterNegToPrimVtxVSMass","neg DCA daughter to prim vtx vs invmass",nbMass,1.05,1.25,250,0.0,10.0);
      fHistPiAPDCAVSMass[0] = new TH2F("fHistPiAPDCAVSMass","pi+p- dca  vs pt",nbMass,1.05,1.25,250,0.0,5.0);
      fHistPiAPDCAZVSMass[0] = new TH2F("fHistPiAPDCAZVSMass","pi+p- dca z vs pt",nbMass,1.05,1.25,200,-20.0,20.0);
      fHistPiAPCosPointAng[0] = new TH2F("fHistPiAPCosPointAng","#bar{#Lambda} cosine of pointing angle vs mass",nbMass,1.05,1.25,200,0.99,1.00);
      fHistPiAPRadiusXY[0] = new TH2F("fHistPiAPRadiusXY","pi+p- phi dist vs mass",nbMass,1.05,1.25,200,0.0,4.0);
      // fHistPiAPPtDaughters[0] = new TH2F("fHistPiAPPtDaughters","p_{t} pos vs p_{t} neg of daughters",400,0.0,20.0,400,0,20.0);
      //	fHistPiAPDCAZPos[0] = new TH2F("fHistPiAPDCAZPos","dca z  of ALambda pos daughters",nbMass,1.05,1.25,200,-20.0,20.0);
      //fHistPiAPDCAZNeg[0] = new TH2F("fHistPiAPDCAZNeg","dca z  of ALambda neg daughters",nbMass,1.05,1.25,200,-20.0,20.0);
      fHistPiAPTrackLengthPosVsMass[0] = new TH2F("fHistPiAPTrackLengthPosVsMass","track length of pos ALambda daughter in TPC",nbMass,1.05,1.25,250,0.0,250.0);
      fHistPiAPTrackLengthNegVsMass[0] = new TH2F("fHistPiAPTrackLengthNegVsMass","track length of neg ALambda daughter in TPC",nbMass,1.05,1.25,250,0.0,250.0);
    }
    else{//pt dependence
      fHistPiAPDCADaughters[0] = new TH2F("fHistPiAPDCADaughters","dca of #bar{#Lambda} daughters",200,0.0,20.0,250,0.0,2.0);
      fHistPiAPDCADaughterPosToPrimVtxVSMass[0] = new TH2F("fHistPiAPDCADaughterPosToPrimVtxVSMass","pos DCA daughter to prim vtx vs invmass",200,0.0,20.0,250,0.0,10.0);
      fHistPiAPDCADaughterNegToPrimVtxVSMass[0] = new TH2F("fHistPiAPDCADaughterNegToPrimVtxVSMass","neg DCA daughter to prim vtx vs invmass",200,0.0,20.0,250,0.0,10.0);
      fHistPiAPDCAVSMass[0] = new TH2F("fHistPiAPDCAVSMass","pi+p- dca  vs pt",200,0.0,20.0,250,0.0,5.0);
      fHistPiAPDCAZVSMass[0] = new TH2F("fHistPiAPDCAZVSMass","pi+p- dca z vs pt",200,0.0,20.0,200,-20.0,20.0);
      fHistPiAPCosPointAng[0] = new TH2F("fHistPiAPCosPointAng","#bar{#Lambda} cosine of pointing angle vs mass",200,0.0,20.0,200,0.99,1.00);
      fHistPiAPRadiusXY[0] = new TH2F("fHistPiAPRadiusXY","pi+p- phi dist vs mass",200,0.0,20.0,200,0.0,4.0);
      //	fHistPiAPDCAZPos[0] = new TH2F("fHistPiAPDCAZPos","dca z  of ALambda pos daughters",200,0.0,20.0,200,-20.0,20.0);
      //fHistPiAPDCAZNeg[0] = new TH2F("fHistPiAPDCAZNeg","dca z  of ALambda neg daughters",200,0.0,20.0,200,-20.0,20.0);
      fHistPiAPTrackLengthPosVsMass[0] = new TH2F("fHistPiAPTrackLengthPosVsMass","track length of pos ALambda daughter in TPC",200,0.0,20.0,250,0.0,250.0);
      fHistPiAPTrackLengthNegVsMass[0] = new TH2F("fHistPiAPTrackLengthNegVsMass","track length of neg ALambda daughter in TPC",200,0.0,20.0,250,0.0,250.0);
    }
   
    //------dedx--------//
    fHistDedxSecProt[0] = new TH2F("fHistDedxSecProt","proton", 250,0.0,5.0, 100, 0, 400);
    fHistDedxSecPiPlus[0] = new TH2F("fHistDedxSecPiPlus","pi plus", 250,0.0,5.0, 100, 0, 400);
    fHistDedxSecAProt[0] = new TH2F("fHistDedxSecAProt","antiproton", 250,0.0,5.0, 100, 0, 400);
    fHistDedxSecPiMinus[0] = new TH2F("fHistDedxSecPiMinus","pi minus", 250,0.0,5.0, 100, 0, 400);
    fHistDedxProt[0] = new TH2F("fHistDedxProt","proton", 250,0.0,5.0, 100, 0, 400);
    fHistDedxPiPlus[0] = new TH2F("fHistDedxPiPlus","pi plus", 250,0.0,5.0, 100, 0, 400);
    fHistDedxAProt[0] = new TH2F("fHistDedxAProt","antiproton", 250,0.0,5.0, 100, 0, 400);
    fHistDedxPiMinus[0] = new TH2F("fHistDedxPiMinus","pi minus", 250,0.0,5.0, 100, 0, 400);


    // ------------------------------------------ clusters --------------------------------------------------//
    fHistNclsITS[0] = new TH2F("fHistNclsITS","fHistNclsITS pos vs neg L",10,-0.5,9.5,10,-0.5,9.5);
    fHistNclsTPC[0] = new TH2F("fHistNclsTPC","ncls TPC neg vs crossed rows neg L",200,-0.5,199.5,200,-0.5,199.5);

    fHistNclsITS[1] = new TH2F("fHistNclsITSSec","fHistNclsITS pos vs neg K0",10,-0.5,9.5,10,-0.5,9.5);
    fHistNclsTPC[1] = new TH2F("fHistNclsTPCSec","ncls TPC neg vs crossed rows neg K0",200,-0.5,199.5,200,-0.5,199.5);

    if(!fSetPtDepHist){
      //K0s
      fHistNclsITSPosK0 = new TH2F("fHistNclsITSPosK0","fHistNclsITSPos  vs pt K0 pos",nbMass,0.25,0.75,7,-0.5,6.5);
      fHistNclsITSNegK0 = new TH2F("fHistNclsITSNegK0","fHistNclsITSNeg vs pt K0 neg",nbMass,0.25,0.75,7,-0.5,6.5);
	  
      fHistNclsTPCPosK0 = new TH2F("fHistNclsTPCPosK0","K0 mass vs phi pos",nbMass,0.25,0.75,200,0.0,200.0);
      fHistNclsTPCNegK0 = new TH2F("fHistNclsTPCNegK0","K0 mass vs phi neg",nbMass,0.25,0.75,200,0.0,200.0);
	  
      fHistChi2PerNclsITSPosK0 = new TH2F("fHistChi2PerNclsITSPosK0","chi2 per cluster ITS K0 pos",nbMass,0.25,0.75,250,0.0,25.0);
      fHistChi2PerNclsITSNegK0 = new TH2F("fHistChi2PerNclsITSNegK0","chi2 per cluster ITS K0 neg",nbMass,0.25,0.75,250,0.0,25.0);
	  
      fHistNCRowsTPCPosK0 = new TH2F("fHistNCRowsTPCPosK0","n crossed rows vs K0 pos",nbMass,0.25,0.75,200,0.0,200.0);
      fHistNCRowsTPCNegK0 = new TH2F("fHistNCRowsTPCNegK0","n crossed rows vs K0 neg",nbMass,0.25,0.75,200,0.0,200.0);
   
      fHistRatioFoundOverFinableTPCK0Pos = new TH2F("fHistRatioFoundOverFinableTPCK0Pos","ncls found over findable K0 pos sec",nbMass,0.25,0.75,200,0.0,2.0);
      fHistRatioFoundOverFinableTPCK0Neg = new TH2F("fHistRatioFoundOverFinableTPCK0Neg","ncls found over findable K0 neg sec",nbMass,0.25,0.75,200,0.0,2.0);
      //Lambda
      fHistNclsITSPosL[0] = new TH2F("fHistNclsITSPosL","fHistNclsITSPos  vs pt L pos",nbMass,1.05,1.25,7,-0.5,6.5);
      fHistNclsITSNegL[0] = new TH2F("fHistNclsITSNegL","fHistNclsITSNeg vs pt L neg",nbMass,1.05,1.25,7,-0.5,6.5);
	  
      fHistNclsTPCPosL[0] = new TH2F("fHistNclsTPCPosL","L mass vs phi pos",nbMass,1.05,1.25,200,0.0,200.0);
      fHistNclsTPCNegL[0] = new TH2F("fHistNclsTPCNegL","L mass vs phi neg",nbMass,1.05,1.25,200,0.0,200.0);
	  
      fHistChi2PerNclsITSPosL[0] = new TH2F("fHistChi2PerNclsITSPosL","chi2 per cluster ITS L pos",nbMass,1.05,1.25,250,0.0,25.0);
      fHistChi2PerNclsITSNegL[0] = new TH2F("fHistChi2PerNclsITSNegL","chi2 per cluster ITS L neg",nbMass,1.05,1.25,250,0.0,25.0);
	  
      fHistNCRowsTPCPosL[0] = new TH2F("fHistNCRowsTPCPosL","n crossed rows vs L pos",nbMass,1.05,1.25,200,0.0,200.0);
      fHistNCRowsTPCNegL[0] = new TH2F("fHistNCRowsTPCNegL","n crossed rows vs L neg",nbMass,1.05,1.25,200,0.0,200.0);
   
      fHistRatioFoundOverFinableTPCLPos[0] = new TH2F("fHistRatioFoundOverFinableTPCLPos","ncls found over findable L pos sec",nbMass,1.05,1.25,200,0.0,2.0);
      fHistRatioFoundOverFinableTPCLNeg[0] = new TH2F("fHistRatioFoundOverFinableTPCLNeg","ncls found over findable L neg sec",nbMass,1.05,1.25,200,0.0,2.0);
    }
    else{//pt dependence
      //K0s
      fHistNclsITSPosK0 = new TH2F("fHistNclsITSPosK0","fHistNclsITSPos  vs pt L pos",200,0.0,20.0,7,-0.5,6.5);
      fHistNclsITSNegK0 = new TH2F("fHistNclsITSNegK0","fHistNclsITSNeg vs pt L neg",200,0.0,20.0,7,-0.5,6.5);
	  
      fHistNclsTPCPosK0 = new TH2F("fHistNclsTPCPosK0","L mass vs phi pos",200,0.0,20.0,200,0.0,200.0);
      fHistNclsTPCNegK0 = new TH2F("fHistNclsTPCNegK0","L mass vs phi neg",200,0.0,20.0,200,0.0,200.0);
	  
      fHistChi2PerNclsITSPosK0 = new TH2F("fHistChi2PerNclsITSPosK0","chi2 per cluster ITS L pos",200,0.0,20.0,250,0.0,25.0);
      fHistChi2PerNclsITSNegK0 = new TH2F("fHistChi2PerNclsITSNegK0","chi2 per cluster ITS L neg",200,0.0,20.0,250,0.0,25.0);
	  
      fHistNCRowsTPCPosK0 = new TH2F("fHistNCRowsTPCPosK0","n crossed rows vs L pos",200,0.0,20.0,200,0.0,200.0);
      fHistNCRowsTPCNegK0 = new TH2F("fHistNCRowsTPCNegK0","n crossed rows vs L neg",200,0.0,20.0,200,0.0,200.0);
   
      fHistRatioFoundOverFinableTPCK0Pos = new TH2F("fHistRatioFoundOverFinableTPCK0Pos","ncls found over findable L pos sec",200,0.0,20.0,200,0.0,2.0);
      fHistRatioFoundOverFinableTPCK0Neg = new TH2F("fHistRatioFoundOverFinableTPCK0Neg","ncls found over findable L neg sec",200,0.0,20.0,200,0.0,2.0);
      //Lambda
      fHistNclsITSPosL[0] = new TH2F("fHistNclsITSPosL","fHistNclsITSPos  vs pt L pos",200,0.0,20.0,7,-0.5,6.5);
      fHistNclsITSNegL[0] = new TH2F("fHistNclsITSNegL","fHistNclsITSNeg vs pt L neg",200,0.0,20.0,7,-0.5,6.5);
	  
      fHistNclsTPCPosL[0] = new TH2F("fHistNclsTPCPosL","L mass vs phi pos",200,0.0,20.0,200,0.0,200.0);
      fHistNclsTPCNegL[0] = new TH2F("fHistNclsTPCNegL","L mass vs phi neg",200,0.0,20.0,200,0.0,200.0);
	  
      fHistChi2PerNclsITSPosL[0] = new TH2F("fHistChi2PerNclsITSPosL","chi2 per cluster ITS L pos",200,0.0,20.0,250,0.0,25.0);
      fHistChi2PerNclsITSNegL[0] = new TH2F("fHistChi2PerNclsITSNegL","chi2 per cluster ITS L neg",200,0.0,20.0,250,0.0,25.0);
	  
      fHistNCRowsTPCPosL[0] = new TH2F("fHistNCRowsTPCPosL","n crossed rows vs L pos",200,0.0,20.0,200,0.0,200.0);
      fHistNCRowsTPCNegL[0] = new TH2F("fHistNCRowsTPCNegL","n crossed rows vs L neg",200,0.0,20.0,200,0.0,200.0);
   
      fHistRatioFoundOverFinableTPCLPos[0] = new TH2F("fHistRatioFoundOverFinableTPCLPos","ncls found over findable L pos sec",200,0.0,20.0,200,0.0,2.0);
      fHistRatioFoundOverFinableTPCLNeg[0] = new TH2F("fHistRatioFoundOverFinableTPCLNeg","ncls found over findable L neg sec",200,0.0,20.0,200,0.0,2.0);
    }

    // --------------------------------------------- for MC reco secondaries -----------------------------------------//
    if(mchist==2){// for MC reco

      //-----------------K0s---------------------//
      //----------------Lambda-------------------//
      if(!fSetPtDepHist){
	fHistPiPDCADaughters[1] = new TH2F("fHistPiPDCADaughtersSec","dca of #Lambda daughters",nbMass,1.05,1.25,250,0.0,2.0);
	fHistPiPDCADaughterPosToPrimVtxVSMass[1] = new TH2F("fHistPiPDCADaughterPosToPrimVtxVSMassSec","pos sec DCA daughter to prim vtx vs invmass",nbMass,1.05,1.25,250,0.0,10.0);
	fHistPiPDCADaughterNegToPrimVtxVSMass[1] = new TH2F("fHistPiPDCADaughterNegToPrimVtxVSMassSec","neg sec DCA daughter to prim vtx vs invmass",nbMass,1.05,1.25,250,0.0,10.0);
	fHistPiPDCAVSMass[1] = new TH2F("fHistPiPDCAVSMassSec","ppi- dca  vs pt",nbMass,1.05,1.25,250,0.0,5.0);
	fHistPiPDCAZVSMass[1] = new TH2F("fHistPiPDCAZVSMassSec","ppi- dca z vs pt",nbMass,1.05,1.25,200,-20.0,20.0);
	fHistPiPCosPointAng[1]  = new TH2F("fHistPiPCosPointAngSec","#Lambda cosine of pointing angle vs mass",nbMass,1.05,1.25,200,0.99,1.00);
	//	 fHistPiPDecayLengthVsMass[1] = new TH2F("fHistPiPDecayLengthVsMassSec","#Lambda decay length vs mass",nbMass,1.05,1.25,200,0.0,100.0);
	fHistPiPRadiusXY[1] = new TH2F("fHistPiPRadiusXYSec","pi-p+ phi dist vs mass",nbMass,1.05,1.25,200,0.0,4.0);
	// fHistPiPPtDaughters[0] = new TH2F("fHistPiPPtDaughters","p_{t} pos vs p_{t} neg of daughters",400,0.0,20.0,400,0,20.0);
	fHistPiPDCAZPos[1] = new TH2F("fHistPiPDCAZPosSec","dca z  of Lambda sec pos daughters",nbMass,1.05,1.25,200,-20.0,20.0);
	fHistPiPDCAZNeg[1] = new TH2F("fHistPiPDCAZNegSec","dca z  of Lambda sec neg daughters",nbMass,1.05,1.25,200,-20.0,20.0);
	fHistPiPTrackLengthPosVsMass[1] = new TH2F("fHistPiPTrackLengthPosVsMassSec","track length of pos sec Lambda daughter in TPC",nbMass,1.05,1.25,250,0.0,250.0);
	fHistPiPTrackLengthNegVsMass[1] = new TH2F("fHistPiPTrackLengthNegVsMassSec","track length of neg sec Lambda daughter in TPC",nbMass,1.05,1.25,250,0.0,250.0);
      }
      else{
	fHistPiPDCADaughters[1] = new TH2F("fHistPiPDCADaughtersSec","dca of #Lambda daughters",200,0.0,20.0,250,0.0,2.0);
	fHistPiPDCADaughterPosToPrimVtxVSMass[1] = new TH2F("fHistPiPDCADaughterPosToPrimVtxVSMassSec","pos sec DCA daughter to prim vtx vs invmass",200,0.0,20.0,250,0.0,10.0);
	fHistPiPDCADaughterNegToPrimVtxVSMass[1] = new TH2F("fHistPiPDCADaughterNegToPrimVtxVSMassSec","neg sec DCA daughter to prim vtx vs invmass",200,0.0,20.0,250,0.0,10.0);
	fHistPiPDCAVSMass[1] = new TH2F("fHistPiPDCAVSMassSec","ppi- dca  vs pt",200,0.0,20.0,250,0.0,5.0);
	fHistPiPDCAZVSMass[1] = new TH2F("fHistPiPDCAZVSMassSec","ppi- dca  vs pt",200,0.0,20.0,200,-20.0,20.0);
	fHistPiPCosPointAng[1]  = new TH2F("fHistPiPCosPointAngSec","#Lambda cosine of pointing angle vs mass",200,0.0,20.0,200,0.99,1.00);
	fHistPiPRadiusXY[1] = new TH2F("fHistPiPRadiusXYSec","pi-p+ phi dist vs mass",200,0.0,20.0,200,0.0,4.0);
	fHistPiPDCAZPos[1] = new TH2F("fHistPiPDCAZPosSec","dca z  of Lambda sec pos daughters",200,0.0,20.0,200,-20.0,20.0);
	fHistPiPDCAZNeg[1] = new TH2F("fHistPiPDCAZNegSec","dca z  of Lambda sec neg daughters",200,0.0,20.0,200,-20.0,20.0);
	fHistPiPTrackLengthPosVsMass[1] = new TH2F("fHistPiPTrackLengthPosVsMassSec","track length of pos sec Lambda daughter in TPC",200,0.0,20.0,250,0.0,250.0);
	fHistPiPTrackLengthNegVsMass[1] = new TH2F("fHistPiPTrackLengthNegVsMassSec","track length of neg sec Lambda daughter in TPC",200,0.0,20.0,250,0.0,250.0);
      }
	  
      //--------------------ALambda--------------//
      if(!fSetPtDepHist){
	fHistPiAPDCADaughters[1] = new TH2F("fHistPiAPDCADaughtersSec","dca of #bar{#Lambda} daughters",nbMass,1.05,1.25,250,0.0,2.0);
	fHistPiAPDCADaughterPosToPrimVtxVSMass[1] = new TH2F("fHistPiAPDCADaughterPosToPrimVtxVSMassSec","pos sec DCA daughter to prim vtx vs invmass",nbMass,1.05,1.25,250,0.0,10.0);
	fHistPiAPDCADaughterNegToPrimVtxVSMass[1] = new TH2F("fHistPiAPDCADaughterNegToPrimVtxVSMassSec","neg sec DCA daughter to prim vtx vs invmass",nbMass,1.05,1.25,250,0.0,10.0);
	fHistPiAPDCAVSMass[1]   = new TH2F("fHistPiAPDCAVSMassSec","pi+p- dca  vs pt",nbMass,1.05,1.25,250,0.0,5.0);
	fHistPiAPDCAZVSMass[1]   = new TH2F("fHistPiAPDCAZVSMassSec","pi+p- dca z  vs pt",nbMass,1.05,1.25,200,-20.0,20.0);
	fHistPiAPCosPointAng[1] = new TH2F("fHistPiAPCosPointAngSec","#bar{#Lambda} cosine of pointing angle vs mass",nbMass,1.05,1.25,200,0.99,1.00);
	//	 fHistPiAPDecayLengthVsMass[1] = new TH2F("fHistPiAPDecayLengthVsMassSec","#bar{#Lambda} decay length vs mass",nbMass,1.05,1.25,200,0.0,100.0);
	fHistPiAPRadiusXY[1] = new TH2F("fHistPiAPRadiusXYSec","pi+p- phi dist vs mass",nbMass,1.05,1.25,200,0.0,4.0);
	// fHistPiAPPtDaughters[0] = new TH2F("fHistPiAPPtDaughters","p_{t} pos vs p_{t} neg of daughters",400,0.0,20.0,400,0,20.0);
	//	  fHistPiAPDCAZPos[1] = new TH2F("fHistPiAPDCAZPosSec","dca z  of ALambda sec pos daughters",nbMass,1.05,1.25,200,-20.0,20.0);
	//fHistPiAPDCAZNeg[1] = new TH2F("fHistPiAPDCAZNegSec","dca z  of ALambda sec neg daughters",nbMass,1.05,1.25,200,-20.0,20.0);
	fHistPiAPTrackLengthPosVsMass[1] = new TH2F("fHistPiAPTrackLengthPosVsMassSec","track length of pos sec ALambda daughter in TPC",nbMass,1.05,1.25,250,0.0,250.0);
	fHistPiAPTrackLengthNegVsMass[1] = new TH2F("fHistPiAPTrackLengthNegVsMassSec","track length of neg sec ALambda daughter in TPC",nbMass,1.05,1.25,250,0.0,250.0);
      }
      else{
	fHistPiAPDCADaughters[1] = new TH2F("fHistPiAPDCADaughtersSec","dca of #bar{#Lambda} daughters",200,0.0,20.0,250,0.0,2.0);
	fHistPiAPDCADaughterPosToPrimVtxVSMass[1] = new TH2F("fHistPiAPDCADaughterPosToPrimVtxVSMassSec","pos sec DCA daughter to prim vtx vs invmass",200,0.0,20.0,250,0.0,10.0);
	fHistPiAPDCADaughterNegToPrimVtxVSMass[1] = new TH2F("fHistPiAPDCADaughterNegToPrimVtxVSMassSec","neg sec DCA daughter to prim vtx vs invmass",200,0.0,20.0,250,0.0,10.0);
	fHistPiAPDCAVSMass[1]   = new TH2F("fHistPiAPDCAVSMassSec","pi+p- dca  vs pt",200,0.0,20.0,250,0.0,5.0);
	fHistPiAPDCAZVSMass[1]   = new TH2F("fHistPiAPDCAZVSMassSec","pi+p- dca z vs pt",200,0.0,20.0,200,-20.0,20.0);
	fHistPiAPCosPointAng[1] = new TH2F("fHistPiAPCosPointAngSec","#bar{#Lambda} cosine of pointing angle vs mass",200,0.0,20.0,200,0.99,1.00);
	fHistPiAPRadiusXY[1] = new TH2F("fHistPiAPRadiusXYSec","pi+p- phi dist vs mass",200,0.0,20.0,200,0.0,4.0);
	//	  fHistPiAPDCAZPos[1] = new TH2F("fHistPiAPDCAZPosSec","dca z  of ALambda sec pos daughters",200,0.0,20.0,200,-20.0,20.0);
	//fHistPiAPDCAZNeg[1] = new TH2F("fHistPiAPDCAZNegSec","dca z  of ALambda sec neg daughters",200,0.0,20.0,200,-20.0,20.0);
	fHistPiAPTrackLengthPosVsMass[1] = new TH2F("fHistPiAPTrackLengthPosVsMassSec","track length of pos sec ALambda daughter in TPC",200,0.0,20.0,250,0.0,250.0);
	fHistPiAPTrackLengthNegVsMass[1] = new TH2F("fHistPiAPTrackLengthNegVsMassSec","track length of neg sec ALambda daughter in TPC",200,0.0,20.0,250,0.0,250.0);
      }

      //-------------dedx------------//
      fHistDedxSecProt[1] = new TH2F("fHistDedxSecProtSec","proton", 250,0.0,5.0, 100, 0, 400);
      fHistDedxSecPiPlus[1] = new TH2F("fHistDedxSecPiPlusSec","pi plus", 250,0.0,5.0, 100, 0, 400);
      fHistDedxSecAProt[1] = new TH2F("fHistDedxSecAProtSec","antiproton", 250,0.0,5.0, 100, 0, 400);
      fHistDedxSecPiMinus[1] = new TH2F("fHistDedxSecPiMinusSec","pi minus", 250,0.0,5.0, 100, 0, 400);
      fHistDedxProt[1] = new TH2F("fHistDedxProtSec","proton", 250,0.0,5.0, 100, 0, 400);
      fHistDedxPiPlus[1] = new TH2F("fHistDedxPiPlusSec","pi plus", 250,0.0,5.0, 100, 0, 400);
      fHistDedxAProt[1] = new TH2F("fHistDedxAProtSec","antiproton", 250,0.0,5.0, 100, 0, 400);
      fHistDedxPiMinus[1] = new TH2F("fHistDedxPiMinusSec","pi minus", 250,0.0,5.0, 100, 0, 400);

      // ------------------------------------------ clusters --------------------------------------------------//
      if(!fSetPtDepHist){
	fHistNclsITSPosL[1] = new TH2F("fHistNclsITSPosLSec","fHistNclsITSPos  vs pt L pos",nbMass,1.05,1.25,7,-0.5,6.5);
	fHistNclsITSNegL[1] = new TH2F("fHistNclsITSNegLSec","fHistNclsITSNeg vs pt L neg",nbMass,1.05,1.25,7,-0.5,6.5);
	  
	fHistNclsTPCPosL[1] = new TH2F("fHistNclsTPCPosLSec","L mass vs phi pos",nbMass,1.05,1.25,200,0.0,200.0);
	fHistNclsTPCNegL[1] = new TH2F("fHistNclsTPCNegLSec","L mass vs phi neg",nbMass,1.05,1.25,200,0.0,200.0);
	  
	fHistChi2PerNclsITSPosL[1] = new TH2F("fHistChi2PerNclsITSPosLSec","chi2 per cluster ITS L pos",nbMass,1.05,1.25,250,0.0,25.0);
	fHistChi2PerNclsITSNegL[1] = new TH2F("fHistChi2PerNclsITSNegLSec","chi2 per cluster ITS L neg",nbMass,1.05,1.25,250,0.0,25.0);
	  
	fHistNCRowsTPCPosL[1] = new TH2F("fHistNCRowsTPCPosLSec","n crossed rows vs L pos",nbMass,1.05,1.25,200,0.0,200.0);
	fHistNCRowsTPCNegL[1] = new TH2F("fHistNCRowsTPCNegLSec","n crossed rows vs L neg",nbMass,1.05,1.25,200,0.0,200.0);
   
	fHistRatioFoundOverFinableTPCLPos[1] = new TH2F("fHistRatioFoundOverFinableTPCLPosSec","ncls found over findable L pos sec",nbMass,1.05,1.25,200,0.0,2.0);
	fHistRatioFoundOverFinableTPCLNeg[1] = new TH2F("fHistRatioFoundOverFinableTPCLNegSec","ncls found over findable L neg sec",nbMass,1.05,1.25,200,0.0,2.0);
      }
      else{
	fHistNclsITSPosL[1] = new TH2F("fHistNclsITSPosLSec","fHistNclsITSPos  vs pt L pos",200,0.0,20.0,7,-0.5,6.5);
	fHistNclsITSNegL[1] = new TH2F("fHistNclsITSNegLSec","fHistNclsITSNeg vs pt L neg",200,0.0,20.0,7,-0.5,6.5);
	  
	fHistNclsTPCPosL[1] = new TH2F("fHistNclsTPCPosLSec","L mass vs phi pos",200,0.0,20.0,200,0.0,200.0);
	fHistNclsTPCNegL[1] = new TH2F("fHistNclsTPCNegLSec","L mass vs phi neg",200,0.0,20.0,200,0.0,200.0);
	  
	fHistChi2PerNclsITSPosL[1] = new TH2F("fHistChi2PerNclsITSPosLSec","chi2 per cluster ITS L pos",200,0.0,20.0,250,0.0,25.0);
	fHistChi2PerNclsITSNegL[1] = new TH2F("fHistChi2PerNclsITSNegLSec","chi2 per cluster ITS L neg",200,0.0,20.0,250,0.0,25.0);
	  
	fHistNCRowsTPCPosL[1] = new TH2F("fHistNCRowsTPCPosLSec","n crossed rows vs L pos",200,0.0,20.0,200,0.0,200.0);
	fHistNCRowsTPCNegL[1] = new TH2F("fHistNCRowsTPCNegLSec","n crossed rows vs L neg",200,0.0,20.0,200,0.0,200.0);
   
	fHistRatioFoundOverFinableTPCLPos[1] = new TH2F("fHistRatioFoundOverFinableTPCLPosSec","ncls found over findable L pos sec",200,0.0,20.0,200,0.0,2.0);
	fHistRatioFoundOverFinableTPCLNeg[1] = new TH2F("fHistRatioFoundOverFinableTPCLNegSec","ncls found over findable L neg sec",200,0.0,20.0,200,0.0,2.0);
      }
	
    }

    //------ ITS TPC clusters --------------//
    fOutputContainer->Add(fHistNclsITS[0]) ;
    fOutputContainer->Add(fHistNclsTPC[0]);
    fOutputContainer->Add(fHistNclsITS[1]);
    fOutputContainer->Add(fHistNclsTPC[1]);

    //-----------K0s ------------------//
    fOutputContainer->Add(fHistPiPiDCAZNeg);
    fOutputContainer->Add(fHistPiPiDCAZPos);
    fOutputContainer->Add(fHistPiPiDCADaughters); 
    fOutputContainer->Add(fHistPiPiDCADaughterPosToPrimVtxVSMass);
    fOutputContainer->Add(fHistPiPiDCAVSMass);
    fOutputContainer->Add(fHistPiPiDCAZVSMass);
    fOutputContainer->Add(fHistPiPiCosPointAng);
    fOutputContainer->Add(fHistPiPiTrackLengthPosVsMass);
    fOutputContainer->Add(fHistPiPiTrackLengthNegVsMass);
    fOutputContainer->Add(fHistPiPiRadiusXY);
    //	fOutputContainer->Add( fHistPiPiPtDaughters);
    fOutputContainer->Add(fHistNclsITSPosK0);
    fOutputContainer->Add(fHistNclsITSNegK0);
    fOutputContainer->Add(fHistNclsTPCPosK0);
    fOutputContainer->Add(fHistNclsTPCNegK0);
    fOutputContainer->Add(fHistChi2PerNclsITSPosK0);
    fOutputContainer->Add(fHistChi2PerNclsITSNegK0);
    fOutputContainer->Add(fHistNCRowsTPCPosK0);
    fOutputContainer->Add(fHistNCRowsTPCNegK0);
    fOutputContainer->Add(fHistRatioFoundOverFinableTPCK0Pos);
    fOutputContainer->Add(fHistRatioFoundOverFinableTPCK0Neg);
    /*
      fOutputContainer->Add(fHistPiPiDistDaughtersPos[0]);
      fOutputContainer->Add(fHistPiPiDistDaughtersNeg[0]);
      fOutputContainer->Add(fHistPiPiDCADaughtersPos[0]);
      fOutputContainer->Add(fHistPiPiDCADaughtersNeg[0]);
      fOutputContainer->Add(fHistPiPiRadAtDCA5cmDaughtersPos[0]);
      fOutputContainer->Add(fHistPiPiRadAtDCA5cmDaughtersNeg[0]);
    
      fOutputContainer->Add(fHistPiPiDistDaughtersPos[1]);
      fOutputContainer->Add(fHistPiPiDistDaughtersNeg[1]);
      fOutputContainer->Add(fHistPiPiDCADaughtersPos[1]);
      fOutputContainer->Add(fHistPiPiDCADaughtersNeg[1]);
      fOutputContainer->Add(fHistPiPiRadAtDCA5cmDaughtersPos[1]);
      fOutputContainer->Add(fHistPiPiRadAtDCA5cmDaughtersNeg[1]);

      fOutputContainer->Add(fHistPiPiDistDaughtersPos[2]);
      fOutputContainer->Add(fHistPiPiDistDaughtersNeg[2]);
      fOutputContainer->Add(fHistPiPiDCADaughtersPos[2]);
      fOutputContainer->Add(fHistPiPiDCADaughtersNeg[2]);
      fOutputContainer->Add(fHistPiPiRadAtDCA5cmDaughtersPos[2]);
      fOutputContainer->Add(fHistPiPiRadAtDCA5cmDaughtersNeg[2]);
    */
    // fOutputContainer->Add(fHistPiPiDistDaughtersTPCEntrVsMass);
    //----------- Lambda Antilambda -------------//

    for(Int_t j=0;j<mchist;j++){
      fOutputContainer->Add(fHistPiPDCADaughters[j]); 
      fOutputContainer->Add(fHistPiAPDCADaughters[j]);
      fOutputContainer->Add( fHistPiPDCADaughterPosToPrimVtxVSMass[j]);
      fOutputContainer->Add( fHistPiPDCADaughterNegToPrimVtxVSMass[j]);
      fOutputContainer->Add( fHistPiAPDCADaughterPosToPrimVtxVSMass[j]);
      fOutputContainer->Add( fHistPiAPDCADaughterNegToPrimVtxVSMass[j]);
      //fOutputContainer->Add( fHistPiPPtDaughters[j]);
      //fOutputContainer->Add( fHistPiAPPtDaughters[j]);
      fOutputContainer->Add(fHistPiPDCAVSMass[j]);
      fOutputContainer->Add(fHistPiPDCAZVSMass[j]);
      fOutputContainer->Add(fHistPiAPDCAVSMass[j]);
      fOutputContainer->Add(fHistPiAPDCAZVSMass[j]);
      fOutputContainer->Add(fHistPiPCosPointAng[j]);
      fOutputContainer->Add(fHistPiAPCosPointAng[j]);
      fOutputContainer->Add(fHistPiPDCAZNeg[j]);
      fOutputContainer->Add(fHistPiPDCAZPos[j]);
      //fOutputContainer->Add(fHistPiAPDCAZNeg[j]);
      //fOutputContainer->Add(fHistPiAPDCAZPos[j]);
      fOutputContainer->Add(fHistPiPTrackLengthPosVsMass[j]);
      fOutputContainer->Add(fHistPiPTrackLengthNegVsMass[j]);
      fOutputContainer->Add(fHistPiAPTrackLengthPosVsMass[j]);
      fOutputContainer->Add(fHistPiAPTrackLengthNegVsMass[j]);     
      fOutputContainer->Add(fHistPiPRadiusXY[j]);
      fOutputContainer->Add(fHistPiAPRadiusXY[j]);

      //--------- dEdx --------------------------//
      fOutputContainer->Add(fHistDedxSecProt[j]);
      fOutputContainer->Add(fHistDedxSecAProt[j]);
      fOutputContainer->Add(fHistDedxSecPiPlus[j]);
      fOutputContainer->Add(fHistDedxSecPiMinus[j]);
      fOutputContainer->Add(fHistDedxProt[j]);
      fOutputContainer->Add(fHistDedxAProt[j]);
      fOutputContainer->Add(fHistDedxPiPlus[j]);
      fOutputContainer->Add(fHistDedxPiMinus[j]);

      //--------- TPC Lambda-----------------//
      fOutputContainer->Add(fHistNclsITSPosL[j]);
      fOutputContainer->Add(fHistNclsITSNegL[j]);
      fOutputContainer->Add(fHistNclsTPCPosL[j]);
      fOutputContainer->Add(fHistNclsTPCNegL[j]);
      fOutputContainer->Add(fHistChi2PerNclsITSPosL[j]);
      fOutputContainer->Add(fHistChi2PerNclsITSNegL[j]);
      fOutputContainer->Add(fHistNCRowsTPCPosL[j]);
      fOutputContainer->Add(fHistNCRowsTPCNegL[j]);
      fOutputContainer->Add(fHistRatioFoundOverFinableTPCLPos[j]);
      fOutputContainer->Add(fHistRatioFoundOverFinableTPCLNeg[j]);
    }  
    /*
      fOutputContainer->Add(fHistPiPDistDaughtersPos[0]);
      fOutputContainer->Add(fHistPiPDistDaughtersNeg[0]);
      fOutputContainer->Add(fHistPiPDCADaughtersPos[0]);
      fOutputContainer->Add(fHistPiPDCADaughtersNeg[0]);
      fOutputContainer->Add(fHistPiPRadAtDCA5cmDaughtersPos[0]);
      fOutputContainer->Add(fHistPiPRadAtDCA5cmDaughtersNeg[0]);
    
      fOutputContainer->Add(fHistPiPDistDaughtersPos[1]);
      fOutputContainer->Add(fHistPiPDistDaughtersNeg[1]);
      fOutputContainer->Add(fHistPiPDCADaughtersPos[1]);
      fOutputContainer->Add(fHistPiPDCADaughtersNeg[1]);
      fOutputContainer->Add(fHistPiPRadAtDCA5cmDaughtersPos[1]);
      fOutputContainer->Add(fHistPiPRadAtDCA5cmDaughtersNeg[1]);

      fOutputContainer->Add(fHistPiPDistDaughtersPos[2]);
      fOutputContainer->Add(fHistPiPDistDaughtersNeg[2]);
      fOutputContainer->Add(fHistPiPDCADaughtersPos[2]);
      fOutputContainer->Add(fHistPiPDCADaughtersNeg[2]);
      fOutputContainer->Add(fHistPiPRadAtDCA5cmDaughtersPos[2]);
      fOutputContainer->Add(fHistPiPRadAtDCA5cmDaughtersNeg[2]);
   

      fOutputContainer->Add(fHistPiAPDistDaughtersPos[0]);
      fOutputContainer->Add(fHistPiAPDistDaughtersNeg[0]);
      fOutputContainer->Add(fHistPiAPDCADaughtersPos[0]);
      fOutputContainer->Add(fHistPiAPDCADaughtersNeg[0]);
      fOutputContainer->Add(fHistPiAPRadAtDCA5cmDaughtersPos[0]);
      fOutputContainer->Add(fHistPiAPRadAtDCA5cmDaughtersNeg[0]);
    
      fOutputContainer->Add(fHistPiAPDistDaughtersPos[1]);
      fOutputContainer->Add(fHistPiAPDistDaughtersNeg[1]);
      fOutputContainer->Add(fHistPiAPDCADaughtersPos[1]);
      fOutputContainer->Add(fHistPiAPDCADaughtersNeg[1]);
      fOutputContainer->Add(fHistPiAPRadAtDCA5cmDaughtersPos[1]);
      fOutputContainer->Add(fHistPiAPRadAtDCA5cmDaughtersNeg[1]);

      fOutputContainer->Add(fHistPiAPDistDaughtersPos[2]);
      fOutputContainer->Add(fHistPiAPDistDaughtersNeg[2]);
      fOutputContainer->Add(fHistPiAPDCADaughtersPos[2]);
      fOutputContainer->Add(fHistPiAPDCADaughtersNeg[2]);
      fOutputContainer->Add(fHistPiAPRadAtDCA5cmDaughtersPos[2]);
      fOutputContainer->Add(fHistPiAPRadAtDCA5cmDaughtersNeg[2]);
    */   
    // fOutputContainer->Add(fHistPiPDistDaughtersTPCEntrVsMass);   
    // fOutputContainer->Add(fHistPiAPDistDaughtersTPCEntrVsMass);
   

  }

  //----------------------------- MC reco or MC truth only --------------------------//
  if((fMCMode && fMCTruthMode) || fMCTruthMode){//mc reco truth only
    if(fAnapp){
      fHistPrimVtxZESDVSNContributorsMC = new TH2F("fHistPrimVtxZESDVSNContributorsMC","prim vtx pos z ESD vs no. of contributers MC",250,-50,50,500,0.0,500.0);
      fOutputContainer->Add(fHistPrimVtxZESDVSNContributorsMC);
      fHistPrimVtxZESDTPCVSNContributorsMC = new TH2F("fHistPrimVtxZESDTPCVSNContributorsMC","prim vtx pos z TPC vs no. of contributers MC",250,-50,50,500,0.0,500.0);
      fOutputContainer->Add(fHistPrimVtxZESDTPCVSNContributorsMC);
      fHistPrimVtxZESDSPDVSNContributorsMC = new TH2F("fHistPrimVtxZESDSPDVSNContributorsMC","prim vtx pos z SPD vs no. of contributers MC",250,-50,50,500,0.0,500.0);
      fOutputContainer->Add(fHistPrimVtxZESDSPDVSNContributorsMC);
    }
    fHistMCVertexZ= new TH1F("fHistMCVertexZ"," z vertex distr in cm MC",500,-50,50);
    fOutputContainer->Add(fHistMCVertexZ);
    fHistPiPCosPointAngXiVsPt= new TH2F("fHistPiPCosPointAngXiVsPt","pi-p cos of pointing angle vs pt from xi",200,0.0,20.0,250,0.99,1.00);
    fOutputContainer->Add(fHistPiPCosPointAngXiVsPt);
    fHistPiAPCosPointAngXiVsPt= new TH2F("fHistPiAPCosPointAngXiVsPt","pi+p- cos of pointing angle vs pt from xi",200,0.0,20.0,250,0.99,1.00);	
    fOutputContainer->Add(fHistPiAPCosPointAngXiVsPt);    
    fHistPiPiEtaDMC[0] = new TH2F("fHistPiPiEtaDMCRaw","K0s daughters etaMC raw",300,-6,6,100,0,20);//
    fOutputContainer->Add(fHistPiPiEtaDMC[0]);
    fHistPiPiEtaDMC[1] = new TH2F("fHistPiPiEtaDMC","K0s daughters etaMC after rap V0 cut",300,-6,6,100,0,20);
    fOutputContainer->Add(fHistPiPiEtaDMC[1]); 
    fHistPiPEtaDMC[0] = new TH2F("fHistPiPEtaDMCRaw","#Lambda daughters etaMC raw",300,-6,6,100,0,20);
    fOutputContainer->Add(fHistPiPEtaDMC[0]); 
    fHistPiPEtaDMC[1] = new TH2F("fHistPiPEtaDMC","#Lambda daughters etaMC after rap V0 cut",300,-6,6,100,0,20);
    fOutputContainer->Add(fHistPiPEtaDMC[1]);

    //-------------K0s---------------//
   
    fHistPiPiDecayLengthResolution = new TH2F("fHistPiPiDecayLengthResolution","K0s decay length resolution MC",220,0.0,110.0,220,0.0,110);
	 
    //-------------Lambda------------//
    fHistPiPDecayLengthResolution[0] = new TH2F("fHistPiPDecayLengthResolution","Lambda decay length resolution MC",220,0.0,110.0,220,0.0,110);
    fHistPiPDecayLengthResolution[1] = new TH2F("fHistPiPDecayLengthResolutionSec","Lambda sec decay length resolution MC",220,0.0,110.0,220,0.0,110);

    fHistPiPMassVSPtSecSigma[0] = new TH2F("fHistPiPMassVSPtSecSigmaMC"," pi-p+ InvMass distribution secondaries from sigma MC",nbMass,1.05,1.25,200,0.,20);
    fHistPiPMassVSPtSecSigma[1] = new TH2F("fHistPiPMassVSPtSecSigma"," pi-p+ InvMass distribution secondaries from Sigma reco",nbMass,1.05,1.25,200,0.,20);
   
    fHistPiPMassVSPtSecXi[0] = new TH2F("fHistPiPMassVSPtSecXiMC"," pi-p+ InvMass distribution secondaries from  xi MC",nbMass,1.05,1.25,200,0.,20);
    fHistPiPMassVSPtSecXi[1] = new TH2F("fHistPiPMassVSPtSecXi"," pi-p+ InvMass distribution secondaries from  xi  reco",nbMass,1.05,1.25,200,0.,20);

    fHistPiPMassVSPtSecXiMCTruth = new TH2F("fHistPiPMassVSPtSecXiMCTruth","Lambda mass reco vs pt sec Lambda from xi MC truth pt",nbMass,1.05,1.25,200,0.0,20.0);
 
    fHistPiPMassVSYSecXi[0] = new TH2F("fHistPiPMassVSYSecXiMC"," pi-p+ InvMass distribution secondaries from xi MC",nbMass,1.05,1.25,100,-2.,2);
    fHistPiPMassVSYSecXi[1] = new TH2F("fHistPiPMassVSYSecXi"," pi-p+ InvMass distribution secondaries from xi reco",nbMass,1.05,1.25,100,-2.,2);

    fHistPiPXi0PtVSLambdaPt[0]= new TH2F("fHistPiPXi0PtVSLambdaPtMC"," pt xi 0 vs pt lambda MC truth",200,0.0,20.0,200,0.0,20.0);
    fHistPiPXi0PtVSLambdaPt[1]= new TH2F("fHistPiPXi0PtVSLambdaPt"," pt xi 0 truth vs pt lambda reco",200,0.0,20.0,200,0.0,20.0);

    fHistPiPXiMinusPtVSLambdaPt[0]= new TH2F("fHistPiPXiMinusPtVSLambdaPtMC","pt xi- vs pt lambda MC truth",200,0.0,20.0,200,0.0,20.0);
    fHistPiPXiMinusPtVSLambdaPt[1]= new TH2F("fHistPiPXiMinusPtVSLambdaPt","pt xi- truth vs pt lambda reco",200,0.0,20.0,200,0.0,20.0);

    fHistPiPOmegaPtVSLambdaPt[0] = new TH2F("fHistPiPOmegaPtVSLambdaPtMC","pt omega vs pt lambda MC truth",200,0.0,20.0,200,0.0,20.0);
    fHistPiPOmegaPtVSLambdaPt[1] = new TH2F("fHistPiPOmegaPtVSLambdaPt","pt omega vs pt lambda MC reco",200,0.0,20.0,200,0.0,20.0);

    fHistPiPMassVSPtSecOmega[0] = new TH2F("fHistPiPMassVSPtSecOmegaMC","Lambda mass vs pt omega MCtruth",nbMass,1.05,1.25,200,0.0,20.0);
    fHistPiPMassVSPtSecOmega[1] = new TH2F("fHistPiPMassVSPtSecOmega","Lambda mass vs pt omega MCreco",nbMass,1.05,1.25,200,0.0,20.0);
    fHistPiPMassVSPtSecOmegaMCTruth= new TH2F("fHistPiPMassVSPtSecOmegaMCTruth","Lambda mass vs pt sec Lambda from Omega MC truth pt",nbMass,1.05,1.25,200,0.0,20.0);

    //--------------ALambda-----------------//
    fHistPiAPDecayLengthResolution[0] = new TH2F("fHistPiAPDecayLengthResolution","ALambda decay length resolution MC",220,0.0,110.0,220,0.0,110);
    fHistPiAPDecayLengthResolution[1] = new TH2F("fHistPiAPDecayLengthResolutionSec","ALambda sec decay length resolution MC",220,0.0,110.0,220,0.0,110);

    fHistPiAPMassVSPtSecSigma[0] = new TH2F("fHistPiAPMassVSPtSecSigmaMC"," pi+p- InvMass distribution secondaries from Sigma MC",nbMass,1.05,1.25,200,0.,20);
    fHistPiAPMassVSPtSecSigma[1] = new TH2F("fHistPiAPMassVSPtSecSigma"," pi+p- InvMass distribution secondaries from  Sigma  reco",nbMass,1.05,1.25,200,0.,20);

    fHistPiAPMassVSPtSecXi[0] = new TH2F("fHistPiAPMassVSPtSecXiMC"," pi+p- InvMass distribution secondaries from xi MC",nbMass,1.05,1.25,200,0.,20);
    fHistPiAPMassVSPtSecXi[1] = new TH2F("fHistPiAPMassVSPtSecXi"," pi+p- InvMass distribution secondaries from  Xi reco",nbMass,1.05,1.25,200,0.,20);

    fHistPiAPMassVSPtSecXiMCTruth = new TH2F("fHistPiAPMassVSPtSecXiMCTruth","ALambda mass reco vs pt sec Lambda from xi MC truth pt",nbMass,1.05,1.25,200,0.0,20.0);
      
    fHistPiAPMassVSYSecXi[0] = new TH2F("fHistPiAPMassVSYSecXiMC"," pi+p- InvMass distribution secondaries from  xi MC",nbMass,1.05,1.25,100,-2,2);
    fHistPiAPMassVSYSecXi[1] = new TH2F("fHistPiAPMassVSYSecXi"," pi+p- InvMass distribution secondaries from xi reco",nbMass,1.05,1.25,100,-2.,2);

     
    fHistPiAPXi0PtVSLambdaPt[0]= new TH2F("fHistPiAPXi0PtVSLambdaPtMC"," pt xi 0 vs pt Alambda MC truth",200,0.0,20.0,200,0.0,20.0);
    fHistPiAPXi0PtVSLambdaPt[1]= new TH2F("fHistPiAPXi0PtVSLambdaPt"," pt xi 0 truth vs pt Alambda reco",200,0.0,20.0,200,0.0,20.0);

    fHistPiAPXiMinusPtVSLambdaPt[0]= new TH2F("fHistPiAPXiMinusPtVSLambdaPtMC","pt xi- vs pt Alambda MC truth",200,0.0,20.0,200,0.0,20.0);
    fHistPiAPXiMinusPtVSLambdaPt[1]= new TH2F("fHistPiAPXiMinusPtVSLambdaPt","pt xi- truth vs pt Alambda reco",200,0.0,20.0,200,0.0,20.0);

    fHistPiAPOmegaPtVSLambdaPt[0] = new TH2F("fHistPiAPOmegaPtVSLambdaPtMC","pt omega vs pt alambda MC truth",200,0.0,20.0,200,0.0,20.0);
    fHistPiAPOmegaPtVSLambdaPt[1] = new TH2F("fHistPiAPOmegaPtVSLambdaPt","pt omega vs pt alambda MC reco",200,0.0,20.0,200,0.0,20.0);

    fHistPiAPMassVSPtSecOmega[0] = new TH2F("fHistPiAPMassVSPtSecOmegaMC","ALambda mass vs pt omega MCtruth",nbMass,1.05,1.25,200,0.0,20.0);
    fHistPiAPMassVSPtSecOmega[1] = new TH2F("fHistPiAPMassVSPtSecOmega","ALambda mass vs pt omega MCreco",nbMass,1.05,1.25,200,0.0,20.0);
    fHistPiAPMassVSPtSecOmegaMCTruth= new TH2F("fHistPiAPMassVSPtSecOmegaMCTruth","ALambda mass vs pt sec Lambda from Omega MC truth pt",nbMass,1.05,1.25,200,0.0,20.0);
   
    fOutputContainer->Add(fHistPiPMassVSPtSecXiMCTruth);
    fOutputContainer->Add(fHistPiPMassVSPtSecOmegaMCTruth);

    fOutputContainer->Add(fHistPiAPMassVSPtSecXiMCTruth);
    fOutputContainer->Add(fHistPiAPMassVSPtSecOmegaMCTruth);

    fOutputContainer->Add(fHistPiPiDecayLengthResolution);   
  
  
    for(Int_t j=0;j<2;j++){

      fOutputContainer->Add(fHistPiPDecayLengthResolution[j]);  
      fOutputContainer->Add(fHistPiAPDecayLengthResolution[j]);
      fOutputContainer->Add(fHistPiPMassVSPtSecXi[j]);
      fOutputContainer->Add(fHistPiAPMassVSPtSecXi[j]);
      fOutputContainer->Add(fHistPiPMassVSYSecXi[j]);
      fOutputContainer->Add(fHistPiAPMassVSYSecXi[j]);
      fOutputContainer->Add(fHistPiPXi0PtVSLambdaPt[j]);
      fOutputContainer->Add(fHistPiAPXi0PtVSLambdaPt[j]);
      fOutputContainer->Add(fHistPiPXiMinusPtVSLambdaPt[j]);
      fOutputContainer->Add(fHistPiAPXiMinusPtVSLambdaPt[j]);
      fOutputContainer->Add(fHistPiPMassVSPtSecSigma[j]);
      fOutputContainer->Add(fHistPiAPMassVSPtSecSigma[j]);
      fOutputContainer->Add(fHistPiPOmegaPtVSLambdaPt[j]);
      fOutputContainer->Add(fHistPiAPOmegaPtVSLambdaPt[j]);
      fOutputContainer->Add(fHistPiPMassVSPtSecOmega[j]);
      fOutputContainer->Add(fHistPiAPMassVSPtSecOmega[j]);
    }
  }
  if(fMCMode ||fMCTruthMode ){
    fHistPiPiPDGCode = new TH1F("fHistPiPiPDGCode","PDG code of K0s mothers",3503,-2.5,3500.5);
    fOutputContainer->Add(fHistPiPiPDGCode);
    fHistPiPPDGCode = new TH1F("fHistPiPPDGCode","PDG code of #Lambda  mothers",3503,-2.5,3500.5);
    fOutputContainer->Add(fHistPiPPDGCode);
    fHistPiAPPDGCode = new TH1F("fHistPiAPPDGCode","PDG code of #bar{#Lambda} mothers",3503,-2.5,3500.5);
    fOutputContainer->Add(fHistPiAPPDGCode);
  }
  /*
    if(fMCMode && !fMCTruthMode){
    //K0s
    fHistPiPiGA= new TH2F("fHistPiPiGA","photons BG  vs pt K0 ",nbMass,0.25,0.75,200,0,20.0);
    fOutputContainer->Add(fHistPiPiGA);
    fHistPiPiKch= new TH2F("fHistPiPiKch","ch kaons BG  vs pt K0 ",nbMass,0.25,0.75,200,0,20.0);
    fOutputContainer->Add(fHistPiPiKch);
    fHistPiPiPhi= new TH2F("fHistPiPiPhi","phi BG  vs pt K0 ",nbMass,0.25,0.75,200,0,20.0);
    fOutputContainer->Add(fHistPiPiPhi);
    fHistPiPiL= new TH2F("fHistPiPiL","Lambda BG  vs pt K0 ",nbMass,0.25,0.75,200,0,20.0);
    fOutputContainer->Add(fHistPiPiL);
    fHistPiPiPi0= new TH2F("fHistPiPiPi0","pi0 BG vs pt K0 ",nbMass,0.25,0.75,200,0,20.0);
    fOutputContainer->Add(fHistPiPiPi0);
    fHistPiPiPich= new TH2F("fHistPiPiPich","ch pi BG  vs pt K0 ",nbMass,0.25,0.75,200,0,20.0);
    fOutputContainer->Add(fHistPiPiPich);
    fHistPiPiRoh= new TH2F("fHistPiPiRoh","roh BG  vs pt K0 ",nbMass,0.25,0.75,200,0,20.0);
    fOutputContainer->Add(fHistPiPiRoh);
    fHistPiPiOmega= new TH2F("fHistPiPiOmega","omega BG  vs pt K0 ",nbMass,0.25,0.75,200,0,20.0);
    fOutputContainer->Add(fHistPiPiOmega);
    fHistPiPiKStar= new TH2F("fHistPiPiKStar","Kstar BG  vs pt K0 ",nbMass,0.25,0.75,200,0,20.0);
    fOutputContainer->Add(fHistPiPiKStar);
    fHistPiPiNoMother= new TH2F("fHistPiPiNoMother","combi BG  vs pt K0 ",nbMass,0.25,0.75,200,0,20.0);
    fOutputContainer->Add(fHistPiPiNoMother);

    fHistPiPiK0s= new TH2F("fHistPiPiK0s","K0s BG  vs pt K0 ",nbMass,0.25,0.75,200,0,20.0);
    fOutputContainer->Add(fHistPiPiK0s);
    fHistPiPiK0L= new TH2F("fHistPiPiK0L","K0L BG  vs pt K0 ",nbMass,0.25,0.75,200,0,20.0);
    fOutputContainer->Add(fHistPiPiK0L);
    fHistPiPiN= new TH2F("fHistPiPiN","n BG  vs pt K0 ",nbMass,0.25,0.75,200,0,20.0);
    fOutputContainer->Add(fHistPiPiN);
    fHistPiPiSigma= new TH2F("fHistPiPiSigma","sigma BG  vs pt K0 ",nbMass,0.25,0.75,200,0,20.0);
    fOutputContainer->Add(fHistPiPiSigma);
    fHistPiPiXi= new TH2F("fHistPiPiXi","xi BG  vs pt K0 ",nbMass,0.25,0.75,200,0,20.0);
    fOutputContainer->Add(fHistPiPiXi);
    fHistPiPiDelta= new TH2F("fHistPiPiDelta","delta BG  vs pt K0 ",nbMass,0.25,0.75,200,0,20.0);
    fOutputContainer->Add(fHistPiPiDelta);
    fHistPiPiB= new TH2F("fHistPiPiB","b BG  vs pt K0 ",nbMass,0.25,0.75,200,0,20.0);
    fOutputContainer->Add(fHistPiPiB);
    fHistPiPiD= new TH2F("fHistPiPiD","d BG  vs pt K0 ",nbMass,0.25,0.75,200,0,20.0);
    fOutputContainer->Add(fHistPiPiD);
    fHistPiPiEta= new TH2F("fHistPiPiEta","eta BG  vs pt K0 ",nbMass,0.25,0.75,200,0,20.0);
    fOutputContainer->Add(fHistPiPiEta);



    //Lambda
    fHistPiPGA = new TH2F("fHistPiPGA","photons in L BG",nbMass,1.05,1.25,200,0.0,20.0);
    fOutputContainer->Add(fHistPiPGA);
    fHistPiPKch = new TH2F("fHistPiPKch","ch kaons in L BG",nbMass,1.05,1.25,200,0.0,20.0);
    fOutputContainer->Add(fHistPiPKch);
    fHistPiPK0s = new TH2F("fHistPiPK0s","K0s in L BG",nbMass,1.05,1.25,200,0.0,20.0);
    fOutputContainer->Add(fHistPiPK0s);
    fHistPiPPi0 = new TH2F("fHistPiPPi0","pi0 in L BG",nbMass,1.05,1.25,200,0.0,20.0);
    fOutputContainer->Add(fHistPiPPi0);
    fHistPiPPich = new TH2F("fHistPiPPich","ch pions in L BG",nbMass,1.05,1.25,200,0.0,20.0);
    fOutputContainer->Add(fHistPiPPich);
    fHistPiPKStar = new TH2F("fHistPiPKStar","Kstar in L BG",nbMass,1.05,1.25,200,0.0,20.0);
    fOutputContainer->Add(fHistPiPKStar);
    fHistPiPN = new TH2F("fHistPiPN","neutron in L BG",nbMass,1.05,1.25,200,0.0,20.0y);
    fOutputContainer->Add(fHistPiPN);
    fHistPiPNoMother= new TH2F("fHistPiPNoMother","combi BG  vs pt Lambda ",nbMass,1.05,1.25,200,0,20.0);
    fOutputContainer->Add(fHistPiPNoMother);
    fHistPiPL= new TH2F("fHistPiPL","Lambda BG  vs pt K0 ",nbMass,1.05,1.25,200,0,20.0);
    fOutputContainer->Add(fHistPiPL);
    }
  */

  /*    
  //shift q/pt
  fHistUserPtShift = new TH1F("fHistUserPtShift","user defined shift in 1/pt",100,-0.5,1.5);
  */
   
}

//________________________________________________________________________
void AliAnalysisTaskV0ForRAA::UserExec(Option_t *) {
  //user exec

  //-- esd handler --//
  AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> 
    (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if (!esdH) {
    printf("ERROR: Could not get ESDInputHandler");
    return;
  } 
  fESD = (AliESDEvent*)esdH->GetEvent();
  if(!fESD) {
    printf("ERROR: fESD not available \n");
    return ;
  }

  //-- mc handler --//
  if(fMCMode || fMCTruthMode){
    AliMCEventHandler *mcH = dynamic_cast<AliMCEventHandler*> 
      (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
    if(!mcH) {
      printf("ERROR: Could not get MCInputHandler");
      return;
    }
    fMCev = mcH->MCEvent();
    if (!fMCev) {
      printf("ERROR: fMCev not available \n");
      return ;
    }
  }
    
  //--  AliPIDResponse --//
  fESDpid = esdH->GetPIDResponse();
 
  //-- Count events before cuts --//
  fHistNEvents->Fill(0);

  //-- Check object existence --//
  const AliESDVertex *    vtxESD    = fESD->GetPrimaryVertexTracks();
  const AliESDVertex *    vtxESDTPC = fESD->GetPrimaryVertexTPC();  
  const AliESDVertex *    vtxESDSPD = fESD->GetPrimaryVertexSPD();  
  const AliMultiplicity * multESD   = fESD->GetMultiplicity();  

  if ( !vtxESD ){
    AliError("No Tracks Vertex");
    return;
  }

  if ( !vtxESDTPC ){
    AliError("No TPC Vertex");
    return ;
  }

  if ( !vtxESDSPD ){
    AliError("No SPD Vertex");
    return ;
  }

  if ( !multESD ){
    AliError("No Multiplicity");
    return ;
  }
   

  // ----------- MC vertex -----------------------------------//
 
  Int_t nContr =0;
  
  if(fMCTruthMode){
    Double_t vVertexPrim[3];
    fMCev->GetPrimaryVertex()->GetXYZ(vVertexPrim);
    fHistMCVertexZ->Fill(vVertexPrim[2]);
    
    if(fMCMode && fAnapp){
      if (vtxESD->GetStatus()){
	nContr=vtxESD->GetNContributors();
	fHistPrimVtxZESDVSNContributorsMC->Fill(vVertexPrim[2],nContr);
	fHistPrimVtxZESDTPCVSNContributorsMC->Fill(vVertexPrim[2],nContr);
      }
      else {
	if(vtxESDSPD->GetStatus()){
	  nContr=vtxESDSPD->GetNContributors();
	  fHistPrimVtxZESDTPCVSNContributorsMC->Fill(vVertexPrim[2],nContr);
	  fHistPrimVtxZESDSPDVSNContributorsMC->Fill(vVertexPrim[2],nContr);
	}
	else{
	  fHistPrimVtxZESDVSNContributorsMC->Fill(vVertexPrim[2],nContr);//add for correction ESD and ESDPSD!!!!
	  fHistPrimVtxZESDTPCVSNContributorsMC->Fill(vVertexPrim[2],nContr);
	}
      }
    }
  }
  
     
  
  //-- Check fo centrality --//
  Bool_t process = kTRUE;
  Int_t centBin = -1;
  if(fUseCentrality) {
    centBin = CalculateCentralityBin();
    if(!fUseCentralityRange){
      if(centBin!= fUseCentralityBin) process=kFALSE;
    }
    else if(centBin < fUseCentralityBin || centBin > fUseCentralityBin+fUseCentralityRange)
      process = kFALSE;
  }

  AliESDVZERO* esdV0 = fESD->GetVZEROData();
  Float_t multV0 = esdV0->GetMTotV0A() + esdV0->GetMTotV0C();
  
  if(fAnapp){// pp Analysis
  
    // SDD test for 2.76TeV pp
    // select events with SDD
    //   TString trCl = fESD->GetFiredTriggerClasses();
    //if(!(trCl.Contains("ALLNOTRD")) && fSelSDD) return;
    UInt_t maskSel = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
    if(maskSel& AliVEvent::kFastOnly && fSelSDD) return;
    if(!(maskSel& AliVEvent::kFastOnly) && fSelNoSDD) return;
	 

    //-- Monitor event cuts --//
    fHistNEvents->Fill(1);

    //---ask for pileup from SPD---//
    Bool_t pileUpSPD = fESD->IsPileupFromSPD();
    if(fRejectPileUpSPD && pileUpSPD) return;
    
    Int_t ntracks = fESD->GetNumberOfTracks();
    for(Int_t i=0;i<ntracks;i++){//check sdd event selection
      AliESDtrack *tr=   fESD->GetTrack(i);
      
      Bool_t sdd0 = tr->HasPointOnITSLayer(0);
      Bool_t sdd1 = tr->HasPointOnITSLayer(1);
      Bool_t sdd2 = tr->HasPointOnITSLayer(2);
      Bool_t sdd3 = tr->HasPointOnITSLayer(3);
      Bool_t sdd4 = tr->HasPointOnITSLayer(4);
      Bool_t sdd5 = tr->HasPointOnITSLayer(5);
       
      fHistITSLayerHits->Fill(Int_t(sdd0)*(-1),ntracks);
      fHistITSLayerHits->Fill(Int_t(sdd1)*1,ntracks);
      fHistITSLayerHits->Fill(Int_t(sdd2)*2,ntracks);
      fHistITSLayerHits->Fill(Int_t(sdd3)*3,ntracks);
      fHistITSLayerHits->Fill(Int_t(sdd4)*4,ntracks);
      fHistITSLayerHits->Fill(Int_t(sdd5)*5,ntracks);
    }
      
    //--vertex selection--//
    if (vtxESD->GetStatus()){
      fHistNEvents->Fill(2);
      fHistESDVertexZ->Fill(vtxESD->GetZ());
      if(fabs(vtxESD->GetZ()) < fVertexZCut){
	fHistMuliplicityRaw->Fill(multV0);
	fHistNEvents->Fill(3);
	fHistNPrim->Fill(nContr);
	
	Process();
	
	fHistMuliplicity->Fill(multV0);
	
	nContr = vtxESD->GetNContributors();
	//  if(nContr<501){
	fHistPrimVtxZESDVSNContributors->Fill(vtxESD->GetZ(),nContr);
	fHistPrimVtxZESDTPCVSNContributors->Fill(vtxESDTPC->GetZ(),nContr);
	//fHistPrimVtxZESDSPDVSNContributorsTPC->Fill(vtxESDSPD->GetZ(),nContr);
	//   }
	fHistPrimVtxZESD->Fill(vtxESD->GetZ());
	fHistPrimVtxZESDTPC->Fill(vtxESDTPC->GetZ());
	// fHistPrimVtxZESDSPD->Fill(vtxESDSPD->GetZ());
	// -- count events after processing
	fHistNEvents->Fill(4);
      }
    }
    else{
      if(vtxESDSPD->GetStatus()){
	fHistNEvents->Fill(2);
	
	fHistESDVertexZ->Fill(vtxESDSPD->GetZ());
	if(fabs(vtxESDSPD->GetZ()) < fVertexZCut){
	  
	  fHistMuliplicityRaw->Fill(multV0);
	  fHistNEvents->Fill(3);
	  fHistNPrim->Fill(nContr);
	  
	  Process();
	  
	  fHistMuliplicity->Fill(multV0);
	  
	  nContr = vtxESDSPD->GetNContributors();
	  //  if(nContr<501){
	  //fHistPrimVtxZESDVSNContributors->Fill(vtxESD->GetZ(),nContr);
	  fHistPrimVtxZESDTPCVSNContributors->Fill(vtxESDTPC->GetZ(),nContr);
	  fHistPrimVtxZESDSPDVSNContributors->Fill(vtxESDSPD->GetZ(),nContr);
	  // }
	  // fHistPrimVtxZESD->Fill(vtxESD->GetZ());
	  fHistPrimVtxZESDTPC->Fill(vtxESDTPC->GetZ());
	  fHistPrimVtxZESDSPD->Fill(vtxESDSPD->GetZ());
	  // -- count events after processing
	  fHistNEvents->Fill(4);
	}
      }
      //else return;
    }
  }
  else{// PbPb analysis
    //-- Monitor event cuts --//
    fHistNEvents->Fill(1);

    if(vtxESD->GetStatus()){
      Double_t vtxZ = vtxESD->GetZ();
      fHistESDVertexZ->Fill(vtxZ);
      if(process){
	fHistNEvents->Fill(2);
	if(fabs(vtxZ) < fVertexZCut){
	  nContr = vtxESD->GetNContributors();
	  fHistMuliplicityRaw->Fill(multV0);
	  fHistNEvents->Fill(3);
	  fHistNPrim->Fill(nContr);
	  Process();
	  fHistMuliplicity->Fill(multV0);
	  fHistPrimVtxZESD->Fill(vtxZ);
	  fHistPrimVtxZESDVSNContributors->Fill(vtxZ,nContr);
	  // -- count events after processing --//
	  fHistCentBin->Fill(centBin);
	  fHistNEvents->Fill(4);
	}
      }
      if(fabs(vtxZ) < fVertexZCut) fHistCentBinRaw->Fill(centBin);
    }
  }
  PostData(1,fOutputContainer);
}

//________________________________________________________________________
void AliAnalysisTaskV0ForRAA::Terminate(Option_t *) {
  //terminate
}

//________________________________________________________________________
void AliAnalysisTaskV0ForRAA::Process(){//run the analysis
  
  Int_t ntracks = fESD->GetNumberOfTracks();
  Int_t count = 0;

  //-- count number of tracks --//
   
  if(!(!fMCMode && fMCTruthMode)){
    for(Int_t i=0;i<ntracks;i++){
      AliESDtrack *track = (AliESDtrack*)fESD->GetTrack(i);
      if(!fESDTrackCuts->AcceptTrack(track)) continue;
      if( track->Eta() > fEtaCutMCDaughtersVal) continue;
      count++;
    }
    fHistMultiplicityPrimary->Fill(count);
  }
   
  //-- check number of V0s in case of data or mc data like analysis--//
  Int_t nV0 = fESD->GetNumberOfV0s();
  if(!fMCTruthMode) if(nV0 < 1) return;
   
  //-- run analysis --//
  if(fMCTruthMode)  V0MCTruthLoop();
  else  V0RecoLoop(0,0,0,0,0.0,0,0.0,0.0);

}

//________________________________________________________________________
void AliAnalysisTaskV0ForRAA::V0MCTruthLoop(){
  //loop over MC truth particles

  //-- get MC stack --//
  AliStack *stack = fMCev->Stack();

  /*
  //histo for user defined shift in charge/pt 
  if(fShift){
  fHistUserPtShift->Fill(fDeltaInvP);
  }
  */
  /*
    AliKFVertex primVtxStart(*(fESD->GetPrimaryVertex()));
    Int_t nTracksPrim=primVtxStart.GetNContributors();
    fHistNPrim->Fill(nTracksPrim);
  */
  /*
  // MC
    
  Int_t mcPrimaries = stack->GetNprimary();
  Int_t mcParticles    = stack->GetNtrack();
    
  fHistMultiplicityPrimary->Fill(mcPrimaries);
  fHistMCMultiplicityTracks->Fill(mcParticles);
    
  // number of V0
  fHistNV0->Fill(nV0);
  if(nTracksPrim>0) {
  fHistNV0WithVertex->Fill(nV0);
  }
  */

  //-- MC truht loop for V0s --//
  for (Int_t iMc = 0; iMc < (stack->GetNtrack()); iMc++){//MC truth loop
    Int_t fillMCtruth= int(fMCTruthMode);
    if(fMCTruthMode){
      fHistPiPiMonitorMCCuts->Fill(1*fillMCtruth);
      fHistPiPMonitorMCCuts[0]->Fill(1*fillMCtruth);
      fHistPiAPMonitorMCCuts[0]->Fill(1*fillMCtruth);
    }
    TParticle *p0 = stack->Particle(iMc);
    if(!p0) continue;

    if(fMCTruthMode){
      fHistPiPiMonitorMCCuts->Fill(2*fillMCtruth);
      fHistPiPMonitorMCCuts[0]->Fill(2*fillMCtruth);
      fHistPiAPMonitorMCCuts[0]->Fill(2*fillMCtruth);
    }



    Int_t pdgCode = p0->GetPdgCode();

    //-------------- only K0s and Lambda ----------//
    if( (pdgCode != 310 ) && ( fabs(pdgCode) != 3122 ) ) continue;
    Int_t fillFlagK0 = (3122- fabs(pdgCode))/(3122-310)*fillMCtruth;
    Int_t fillFlagL = (fabs(pdgCode) - 310)/(3122-310)*(pdgCode+3122)/(2*3122)*fillMCtruth;
    Int_t fillFlagAL = (fabs(pdgCode) - 310)/(3122-310)*(pdgCode-3122)/(-2*3122)*fillMCtruth;
    
    fHistPiPiMonitorMCCuts->Fill(3*fillFlagK0);
    fHistPiPMonitorMCCuts[0]->Fill(3*fillFlagL);
    fHistPiAPMonitorMCCuts[0]->Fill(3*fillFlagAL);
      
    if(p0->GetNDaughters() !=2) continue;
    fHistPiPiMonitorMCCuts->Fill(4*fillFlagK0);
    fHistPiPMonitorMCCuts[0]->Fill(4*fillFlagL);
    fHistPiAPMonitorMCCuts[0]->Fill(4*fillFlagAL);
      
    //-------------- unique ID check-------------- //
    Int_t uniqueID =  p0->GetUniqueID();
    if(uniqueID==13) continue;
      
    fHistPiPiMonitorMCCuts->Fill(5*fillFlagK0);
    fHistPiPMonitorMCCuts[0]->Fill(5*fillFlagL);
    fHistPiAPMonitorMCCuts[0]->Fill(5*fillFlagAL);
      
    //-------------- daughters --------------------//
    Int_t id0  = p0->GetDaughter(0);
    Int_t id1  = p0->GetDaughter(1);
    if(id0<0 || id1 <0) continue;
      
    fHistPiPiMonitorMCCuts->Fill(6*fillFlagK0);
    fHistPiPMonitorMCCuts[0]->Fill(6*fillFlagL);
    fHistPiAPMonitorMCCuts[0]->Fill(6*fillFlagAL);
            
    Int_t pdgCodeD0 = stack->Particle(id0)->GetPdgCode();
    Int_t pdgCodeD1 = stack->Particle(id1)->GetPdgCode();
      
    if(pdgCodeD0 == pdgCodeD1) continue;
    if(pdgCodeD0*pdgCodeD1>0) continue;
      
    fHistPiPiMonitorMCCuts->Fill(7*fillFlagK0);
    fHistPiPMonitorMCCuts[0]->Fill(7*fillFlagL);
    fHistPiAPMonitorMCCuts[0]->Fill(7*fillFlagAL);
            
    if((fabs(pdgCodeD0) != 211 ) && ( fabs(pdgCodeD0) != 2212 )) continue;
    if((fabs(pdgCodeD1) != 211 ) && ( fabs(pdgCodeD1) != 2212 )) continue;
      
    fHistPiPiMonitorMCCuts->Fill(8*fillFlagK0);
    fHistPiPMonitorMCCuts[0]->Fill(8*fillFlagL);
    fHistPiAPMonitorMCCuts[0]->Fill(8*fillFlagAL);
      
    TParticle *p00 =stack->Particle(id0);
    TParticle *p01 =stack->Particle(id1);
    Double_t etaMC00   = p00->Eta();
    Double_t etaMC01   = p01->Eta();

    //----------- unique ID check daughters-------- //
    Int_t uniqueIDdaughter0 = p00->GetUniqueID();
    Int_t uniqueIDdaughter1 = p01->GetUniqueID();
    if (uniqueIDdaughter0 !=4 || uniqueIDdaughter1 !=4 ) continue;
      
    fHistPiPiMonitorMCCuts->Fill(9*fillFlagK0);
    fHistPiPMonitorMCCuts[0]->Fill(9*fillFlagL);
    fHistPiAPMonitorMCCuts[0]->Fill(9*fillFlagAL);

    fHistPiPMonitorMCCuts[1]->Fill(9*fillFlagL);
    fHistPiAPMonitorMCCuts[1]->Fill(9*fillFlagAL);
      
    //------------ check label reco -------------------//
    if(fCheckNegLabelReco || fOnlyFoundRecoV0){   // check label reco
      Bool_t found =kFALSE;
      Int_t label0=0,label1=0;      
      AliESDv0 * v0MIsMC=NULL;
      AliESDtrack *tr0 = NULL;
      AliESDtrack *tr1 = NULL;
      for(Int_t recL=0;recL < fESD->GetNumberOfV0s();recL++){
	v0MIsMC = fESD->GetV0(recL);
	if(!v0MIsMC) continue;
	tr0 = fESD->GetTrack(v0MIsMC->GetPindex());
	tr1 = fESD->GetTrack(v0MIsMC->GetNindex());
	if(tr0 && tr1){
	  label0 = tr0->GetLabel();
	  label1 = tr1->GetLabel();
	  if((fabs(label0) == id0 && fabs(label1) == id1) || 
	     (fabs(label0) == id1 && fabs(label1) == id0)){
	    found =kTRUE;
	    break; 
	  }     
	}
      }
      if(fCheckNegLabelReco && !fOnlyFoundRecoV0) {
	if(found && (label0 <0 || label1 < 0)) continue;
      }
      else{
	if(!found) continue;
	if(fCheckNegLabelReco && found && (label0 <0 || label1 < 0)) continue;
      }
      
    }
    //-----------get geometric properties --------------//
    // DCA of mother to prim vertex = production vertex
    
    //-- primary and secondary vetex --//
    Double_t vVertexPrimMC[3];
    fMCev->GetPrimaryVertex()->GetXYZ(vVertexPrimMC);
    // Double_t x0=p0->Vx(),y0=p0->Vy(),z0=p0->Vz();//mother production vertex
    
    Double_t x=p00->Vx(),y=p00->Vy(),z=p00->Vz();//daughter vertex =V0 decay vertex
    Double_t rx = x - vVertexPrimMC[0];
    Double_t ry = y - vVertexPrimMC[1];
    Double_t rz = z - vVertexPrimMC[2];
    Double_t sdeclength = rx*rx+ry*ry;//+rz*rz;//=p00->Rho();  
    Double_t declength =0.0;
    if(sdeclength>0) declength = sqrt(sdeclength);
    Double_t declength3d = sqrt( rx*rx+ry*ry+rz*rz);
    
    //-- decay radii --//
    Double_t rMC2D  = sqrt(x*x+y*y);
    const  Double_t xyzMC[3] = {x,y,z};
    // Double_t rMC = p00->R();
      

    //-------------------- V0 variables ----------------//
    Double_t rapidity = p0->Y();
    Double_t massV0MC = p0->GetMass();
    Double_t ptV0MC =  p0->Pt();
    Double_t pV0MC =  p0->P();

     
    //----------------- mother variables-----------------//
    Int_t indexMother1  = p0->GetMother(0);
    Int_t isSecd=0;
    Int_t pdgMother =0;
    // Int_t goodMother=1;
    Int_t uniqueIDmother=0;
    Double_t ptXiMother=0.0;
    Double_t rapXiMother = 0.0;
 
    //------check mother and fill mother histos---------//
    Bool_t isPrim= stack->IsPhysicalPrimary(iMc);
   
    if(!isPrim){//secondary
      isSecd = -1;// is secondary V0s
      if(indexMother1 >-1){// && !isPrim){//secondary V0s
	// if(fSelectMBMotherMC && !fMCev->IsFromBGEvent(indexMother1)) continue;//xxx only standard hijing particles for sec. lambdas:not needed
	
	//-- check for mother --//
	TParticle *mother = stack->Particle(indexMother1);
	if(!mother) {
	  Printf("no mother pointer!");continue;
	}
	pdgMother = mother->GetPdgCode();
	fHistPiPMonitorMCCuts[1]->Fill(10*fillFlagL);
	fHistPiAPMonitorMCCuts[1]->Fill(10*fillFlagAL);

	//-- check for injejcted --//
	Bool_t notinjectedMother = kTRUE;
	notinjectedMother = fMCev->IsFromBGEvent(indexMother1);
	
	if(fSelectInjected && !notinjectedMother ) continue;
	fHistPiPMonitorMCCuts[1]->Fill(11*fillFlagL);
	fHistPiAPMonitorMCCuts[1]->Fill(11*fillFlagAL);

	Bool_t isPrimMother= stack->IsPhysicalPrimary(indexMother1);
	if(!isPrimMother) continue;
	fHistPiPMonitorMCCuts[1]->Fill(12*fillFlagL);
	fHistPiAPMonitorMCCuts[1]->Fill(12*fillFlagAL);
      
	uniqueIDmother =  mother->GetUniqueID();

	if(uniqueIDmother==13){
	  continue;
	}
	fHistPiPMonitorMCCuts[1]->Fill(13*fillFlagL);
	fHistPiAPMonitorMCCuts[1]->Fill(13*fillFlagAL);


	//-- fill secondary V0s histos and pdg histos --// 
	ptXiMother  = mother->Pt();
	rapXiMother = mother->Y();
	

	//-- K0s --//
	if(pdgCode==310){
	  if(fabs(pdgMother)==311 || fabs(pdgMother)==313 || fabs(pdgMother)==323 ) isSecd=0; // from K0L,  K0 and K* as primary
	  else fHistPiPiPDGCode->Fill(fabs(pdgMother));
	}
	
	//-- Lambda --//
	if(pdgCode==3122){
	  fHistPiPPDGCode->Fill(fabs(pdgMother));
	  if (//sigma family
	      ( TMath::Abs(pdgMother) == 3112) || //sigma minus
	      ( TMath::Abs(pdgMother) == 3222) || //sigma plus
	      ( TMath::Abs(pdgMother) == 3224) || //sigma *plus
	      ( TMath::Abs(pdgMother) == 3114) || //sigma *minus
	      ( TMath::Abs(pdgMother) == 3214) || //sigma *0 counts as primary????
	      ( TMath::Abs(pdgMother) == 3212)    //sigma 0 counts as primary
	      )
	    {
	      isSecd=0;
	    }
	   
	  if( pdgMother == 3322) //xi0
	    {
	      if(!fRapCutV0 || fabs(rapidity)<fRap){
		if(!fEtaCutMCDaughters  ||  (fabs(etaMC00)<fEtaCutMCDaughtersVal|| fabs(etaMC01)<fEtaCutMCDaughtersVal)){
		  fHistPiPMassVSPtSecXi[0]->Fill(massV0MC,ptV0MC);
		  fHistPiPMassVSYSecXi[0]->Fill(massV0MC,rapidity);
  		  if(!fCutRapXi || fabs(rapXiMother)<fRap) {
		    fHistPiPXi0PtVSLambdaPt[0]->Fill(ptV0MC,ptXiMother);
		    if(fUseXi0) isSecd=1;// is secondary V0s
		    // cout<<"xi 0 rap ok  "<<fRapCutV0<<" "<<rapXiMother<<endl;
		    //  rapXiMotherOK  = kTRUE;
		  }
		}
	      }
	    }

	  if(pdgMother == 3312) //xi minus
	    {
	      if(!fRapCutV0 || fabs(rapidity)<fRap){
		if(!fEtaCutMCDaughters  ||  (fabs(etaMC00)<fEtaCutMCDaughtersVal|| fabs(etaMC01)<fEtaCutMCDaughtersVal)){
		  fHistPiPMassVSPtSecXi[0]->Fill(massV0MC,ptV0MC);
		  fHistPiPMassVSYSecXi[0]->Fill(massV0MC,rapidity);
		  if(!fCutRapXi || fabs(rapXiMother)<fRap) {
		    fHistPiPXiMinusPtVSLambdaPt[0]->Fill(ptV0MC,ptXiMother);
		    // rapXiMotherOK =kTRUE;
		    if(fUseXiM) isSecd=1;// is secondary V0s
		    // cout<<"xi M rap ok  "<<fRapCutV0<<" "<<rapXiMother<<endl;
		  }
		}
	      }
	    }
	  
	  if(pdgMother == 3334)//omega-
	    {
	      //  fHistPiPDCAtoPrimVtxOmega[0]->Fill(p0->GetMass(),dcaV0ToPrimVertex);
	      fHistPiPMassVSPtSecOmega[0]->Fill(massV0MC,ptV0MC);	      
	      if(!fCutRapXi || fabs(rapXiMother)<fRap){
		fHistPiPOmegaPtVSLambdaPt[0]->Fill(ptV0MC,ptXiMother);
		if(fUseOmega) isSecd=1;// is secondary V0s
		//	rapXiMotherOK =kTRUE;
		
	      }
	    }
	}
	
	//-- AntiLambda --//
	if(pdgCode==-3122 ){
	  fHistPiAPPDGCode->Fill(fabs(pdgMother));
	  if (//sigma family
	      ( TMath::Abs(pdgMother) == 3112) ||//sigma minus
	      ( TMath::Abs(pdgMother) == 3222) ||//sigma plus
	      ( TMath::Abs(pdgMother) == 3224) ||//sigma *plus
	      ( TMath::Abs(pdgMother) == 3114) ||//sigma *minus
	      ( TMath::Abs(pdgMother) == 3214) || //sigma *0
	      ( TMath::Abs(pdgMother) == 3212)    //sigma 0 counts as primary
	      )
	    {
	      isSecd=0;
	    }
		   
	  if( pdgMother == -3322) //xi0
	    {
	      if(!fRapCutV0 || fabs(rapidity)<fRap){
		if(!fEtaCutMCDaughters  ||  (fabs(etaMC00)<fEtaCutMCDaughtersVal|| fabs(etaMC01)<fEtaCutMCDaughtersVal)){
		  fHistPiAPMassVSPtSecXi[0]->Fill(massV0MC,ptV0MC);
		  fHistPiAPMassVSYSecXi[0]->Fill(massV0MC,rapidity);
		  if(!fCutRapXi || fabs(rapXiMother)<fRap)    {
		    fHistPiAPXi0PtVSLambdaPt[0]->Fill(ptV0MC,ptXiMother);
		    //  rapXiMotherOK =kTRUE;
		    if(fUseXi0) isSecd=1;// is secondary V0s
		  }
		   
		}
	      }
	    }
	  
	  if(pdgMother == -3312) //xi plus
	    {
	      if(!fRapCutV0 || fabs(rapidity)<fRap){
		if(!fEtaCutMCDaughters  ||  (fabs(etaMC00)<fEtaCutMCDaughtersVal|| fabs(etaMC01)<fEtaCutMCDaughtersVal)){
		  fHistPiAPMassVSPtSecXi[0]->Fill(massV0MC,ptV0MC);
		  fHistPiAPMassVSYSecXi[0]->Fill(massV0MC,rapidity);
		  if(!fCutRapXi || fabs(rapXiMother)<fRap) {
		    fHistPiAPXiMinusPtVSLambdaPt[0]->Fill(ptV0MC,ptXiMother);
		    //    rapXiMotherOK =kTRUE;
		    if(fUseXiM)  isSecd=1;// is secondary V0s
		  }
		}
	      }
	    }
	  
	  if(pdgMother == -3334)//omega+
	    {
	      fHistPiAPOmegaPtVSLambdaPt[0]->Fill(ptV0MC,ptXiMother);
	      if(!fCutRapXi || fabs(rapXiMother)<fRap) {
		fHistPiAPMassVSPtSecOmega[0]->Fill(massV0MC,ptV0MC);
		//	rapXiMotherOK =kTRUE;
		if(fUseOmega) isSecd=1;// is secondary V0s
	      }
	      // fHistPiAPDCAtoPrimVtxOmega[0]->Fill(p0->GetMass(),dcaV0ToPrimVertex);
	    }
	}
	


      }	
    }//end secondaries
    else{//primaries
      //-- check for injejcted --//
      Bool_t notinjected = kTRUE;
      notinjected = fMCev->IsFromBGEvent(iMc);
      
      if(fSelectInjected && !notinjected ) continue;
      fHistPiPiMonitorMCCuts->Fill(10*fillFlagK0);
      fHistPiPMonitorMCCuts[0]->Fill(10*fillFlagL);
      fHistPiAPMonitorMCCuts[0]->Fill(10*fillFlagAL);
    }
   
    // if(isSecd == 1 && !rapXiMotherOK &&  fCutRapXi) continue;
    if(isSecd == -1) continue;
    //-------------- MC truth or reco mode -----------------//
    if(fMCTruthMode && !fMCMode){//MC true ana
      fHistPiPiMonitorMCCuts->Fill(14*fillFlagK0);
      fHistPiPMonitorMCCuts[isSecd]->Fill(14*fillFlagL);
      fHistPiAPMonitorMCCuts[isSecd]->Fill(14*fillFlagAL);
      
      //-- DCA daughters --//
      // values of one daugher, should be the same      
      /*
      //to primary vertex
      trackPos->GetImpactParameters(tdcaPosToVertex[0],tdcaPosToVertex[1]);
      trackNeg->GetImpactParameters(tdcaNegToVertex[0],tdcaNegToVertex[1]);
	 
      Double_t dcaPosToVertex = TMath::Sqrt(tdcaPosToVertex[0]*tdcaPosToVertex[0]+tdcaPosToVertex[1]*tdcaPosToVertex[1]);
      Double_t dcaNegToVertex = TMath::Sqrt(tdcaNegToVertex[0]*tdcaNegToVertex[0]+tdcaNegToVertex[1]*tdcaNegToVertex[1]);
      fHistDCADaughtersToPrimVtx[isSecd]->Fill(dcaPosToVertex,dcaNegToVertex);
      */
	 
	 
      //-- armenteros values --//
      TVector3 vecPip;
      TVector3 vecPin;
	
      Double_t ptPlus=0, ptMinus=0;
      Double_t pt00 = p00->Pt();
      Double_t pt01 = p01->Pt();
      Double_t  phiMCPos=0.0;
      Double_t  phiMCNeg=0.0;
      Double_t  etaMCPos =0.0;
      Double_t  etaMCNeg =0.0;
      if(p00->GetPdgCode()<0)
	{
	  vecPip.SetXYZ(p01->Px(),p01->Py(),p01->Pz());
	  vecPin.SetXYZ(p00->Px(),p00->Py(),p00->Pz());
	  ptMinus = pt00;
	  ptPlus = pt01;
	  phiMCPos = p01->Phi();
	  phiMCNeg = p00->Phi();
	  etaMCPos = etaMC01;
	  etaMCNeg = etaMC00;
	}
      else{
	vecPin.SetXYZ(p01->Px(),p01->Py(),p01->Pz());
	vecPip.SetXYZ(p00->Px(),p00->Py(),p00->Pz());
	ptMinus = pt01;
	ptPlus = pt00;
	phiMCPos = p00->Phi();
	phiMCNeg = p01->Phi();
	etaMCPos = etaMC00;
	etaMCNeg = etaMC01;
      }
	    
      TVector3 momTot(p0->Px(),p0->Py(),p0->Pz());
      Double_t lQlNeg = fabs(vecPin.Dot(momTot)/momTot.Mag());
      Double_t lQlPos = fabs(vecPip.Dot(momTot)/momTot.Mag());
      Double_t alfa =0.0;
      Double_t den = lQlPos + lQlNeg;
      if(den>0) alfa = (lQlPos - lQlNeg)/den;
      TVector3 qtvec= vecPin.Cross(momTot);//vecPip.Mag()*sqrt(1-pow(thetapip,2));
      Float_t qt = qtvec.Mag()/momTot.Mag();
      
      //clalc masses for test
      Double_t massPi=0.13957018;
      Double_t massP=0.93827203;

      TLorentzVector pionPTest(vecPip, massPi);
      TLorentzVector pionNTest(vecPin, massPi);
      TLorentzVector k0sTest = pionPTest+pionNTest;

      TLorentzVector protPTest(vecPip, massP);
      TLorentzVector lambdaTest = protPTest+pionNTest;

      TLorentzVector protNTest(vecPin, massP);
      TLorentzVector alambdaTest = protNTest+pionPTest;

      Double_t calcK0smass = fabs(k0sTest.M());
      Double_t calcLambdamass = fabs(lambdaTest.M());
      Double_t calcALambdamass = fabs(alambdaTest.M());

      if(pdgCode == 310) {
	fHistPiPiEtaDMC[isSecd]->Fill(etaMCPos,ptV0MC);
	fHistPiPiEtaDMC[isSecd]->Fill(etaMCNeg,ptV0MC);
      }
      if(fabs(pdgCode) == 3122) {
	fHistPiPEtaDMC[isSecd]->Fill(etaMC00,ptV0MC);
	fHistPiPEtaDMC[isSecd]->Fill(etaMC01,ptV0MC);
      }

      //-- rapidity and eta cut --//      
      if(fRapCutV0 && fabs(rapidity)>fRap) continue;
      fHistPiPiMonitorMCCuts->Fill(15*fillFlagK0);
      fHistPiPMonitorMCCuts[isSecd]->Fill(15*fillFlagL);
      fHistPiAPMonitorMCCuts[isSecd]->Fill(15*fillFlagAL);
	 
      if(fEtaCutMCDaughters) { if(fabs(etaMC00)>fEtaCutMCDaughtersVal || fabs(etaMC01)>fEtaCutMCDaughtersVal ) continue; }
      fHistPiPiMonitorMCCuts->Fill(16*fillFlagK0);
      fHistPiPMonitorMCCuts[isSecd]->Fill(16*fillFlagL);
      fHistPiAPMonitorMCCuts[isSecd]->Fill(16*fillFlagAL);

      /*
	Double_t phiMC =   p0->Phi(); 
	Double_t etaMC =   p0->Eta(); 
      */
      /*
	Double_t valTHnMC[4] = {massV0MC,ptV0MC,etaMC,phiMC};
	Double_t valTHnMCDauEta[4] = {massV0MC,ptV0MC,etaMCPos,etaMCNeg};
	Double_t valTHnMCDauPhi[5] = {massV0MC,phiMCPos,phiMCNeg,0.0,0.0};
      */
      //-- Fill Particle histos --//
      if(pdgCode==310){//K0s
	fHistPiPiMonitorMCCuts->Fill(17);

	fHistPiPiEtaDMC[1]->Fill(etaMC00,ptV0MC);
	fHistPiPiEtaDMC[1]->Fill(etaMC01,ptV0MC);
	  
	fHistPiPiMass->Fill(massV0MC);
	fHistPiPiMassVSPt->Fill(massV0MC,ptV0MC);
	if(rapidity>0.0) fHistPiPiMassVSPtPosMCTruth->Fill(massV0MC,ptV0MC);//check rap dependence of ana
	else fHistPiPiMassVSPtNegMCTruth->Fill(massV0MC,ptV0MC);
	fHistPiPiMassVSY->Fill(massV0MC,rapidity);
	// fHistPiPiPtDaughters->Fill(ptMinus,ptPlus);
	fHistPiPiPtVSY->Fill(rapidity,ptV0MC);
	Double_t ctTK0s=0.0,ctK0s=0.0;
	if(pV0MC>0.0) ctK0s=declength3d*0.497614/pV0MC;
	if(ptV0MC>0.0) ctTK0s=declength*0.497614/ptV0MC;
	fHistPiPiDecayLengthResolution->Fill(declength3d,declength);
	fHistPiPiDecayLengthVsPt->Fill(ptV0MC,declength);//ptV0MC,ctK0s);

	if(!fSetPtDepHist) fHistPiPiDecayLengthVsCtau->Fill(massV0MC,ctTK0s);
	else fHistPiPiDecayLengthVsCtau->Fill(ptV0MC,ctTK0s);
	
	fHistPiPiDecayLengthVsMass->Fill(massV0MC,declength);
	//all V0s histo
	fHistArmenteros[isSecd]->Fill(alfa,qt);
	fHistV0RadiusZ[isSecd]->Fill(rMC2D,xyzMC[2]);
	fHistV0RadiusXY[isSecd]->Fill(xyzMC[0],xyzMC[1]);
	fHistV0RadiusZVSPt[isSecd]->Fill(ptV0MC,xyzMC[2]);
	// fHistPiPiPhiPosVsPtPosVsMass->Fill(massV0MC,ctTK0s,ptV0MC);//,ctK0s);//phiPosMC);//xxx
	fHistPiPiK0sVsLambdaMass->Fill(calcLambdamass,calcK0smass);
	fHistPiPiK0sVsALambdaMass->Fill(calcALambdamass,calcK0smass); 
	/*
	  fTHnFK0s->Fill(valTHnMC);
	  fTHnFK0sDauEta->Fill(valTHnMCDauEta);
	  fTHnFK0sDauPhi->Fill(valTHnMCDauPhi);
	*/
      }
      if (pdgCode==3122){ //Lambda
	fHistPiPMonitorMCCuts[isSecd]->Fill(17);
	
	fHistPiPEtaDMC[1]->Fill(etaMC00,ptV0MC);
	fHistPiPEtaDMC[1]->Fill(etaMC01,ptV0MC);

	fHistPiPMassVSPt[isSecd]->Fill(massV0MC,ptV0MC);
	if(rapidity> 0.0) 	fHistPiPMassVSPtPosMCTruth[isSecd]->Fill(massV0MC,ptV0MC);
	else fHistPiPMassVSPtNegMCTruth[isSecd]->Fill(massV0MC,ptV0MC);
	fHistPiPMass[isSecd]->Fill(massV0MC);  
	fHistPiPMassVSY[isSecd]->Fill(massV0MC,rapidity);
	//  fHistPiPPtDaughters[isSecd]->Fill(ptMinus,ptPlus);
	fHistPiPPtVSY[isSecd]->Fill(rapidity,ptV0MC);
	  
	
	Double_t ctTL=0.0, ctL=0.0;
	if(pV0MC>0.0) ctL=declength3d*1.115683/pV0MC;
	if(ptV0MC>0.0) ctTL=declength*1.115683/ptV0MC;
	fHistPiPDecayLengthResolution[0]->Fill(declength3d,declength);
	fHistPiPDecayLengthVsPt[isSecd]->Fill(ptV0MC,declength);//(ptV0MC,ctL);
	if(!fSetPtDepHist)	fHistPiPDecayLengthVsCtau[isSecd]->Fill(massV0MC,ctTL);
	else 	fHistPiPDecayLengthVsCtau[isSecd]->Fill(ptV0MC,ctTL);
	fHistPiPDecayLengthVsMass[isSecd]->Fill(massV0MC,declength);
	//all V0s hito	
	fHistArmenteros[isSecd]->Fill(alfa,qt);
	fHistV0RadiusZ[isSecd]->Fill(rMC2D,xyzMC[2]);	fHistV0RadiusXY[isSecd]->Fill(xyzMC[0],xyzMC[1]);
	fHistV0RadiusZVSPt[isSecd]->Fill(ptV0MC,xyzMC[2]);
	
	fHistPiPK0sVsLambdaMass->Fill(calcLambdamass,calcK0smass);
	/*
	  fTHnFL->Fill(valTHnMC);
	  fTHnFLDauEta->Fill(valTHnMCDauEta);
	  fTHnFLDauPhi->Fill(valTHnMCDauPhi);
	*/
      }
      if (pdgCode==-3122){ //AntiLambda
	fHistPiAPMonitorMCCuts[isSecd]->Fill(17);
	    
	fHistPiPEtaDMC[1]->Fill(etaMC00,ptV0MC);
	fHistPiPEtaDMC[1]->Fill(etaMC01,ptV0MC);

	fHistPiAPMassVSPt[isSecd]->Fill(massV0MC,ptV0MC);
	if(rapidity>0.0)	fHistPiAPMassVSPtPosMCTruth[isSecd]->Fill(massV0MC,ptV0MC);
	else fHistPiAPMassVSPtNegMCTruth[isSecd]->Fill(massV0MC,ptV0MC);
	fHistPiAPMass[isSecd]->Fill(massV0MC);
	fHistPiPMassVSY[isSecd]->Fill(massV0MC,rapidity);
	//  fHistPiAPPtDaughters[isSecd]->Fill(ptMinus,ptPlus);
	fHistPiAPPtVSY[isSecd]->Fill(rapidity,ptV0MC);
	
	Double_t ctTAL=0.0, ctAL=0.0;
	if(pV0MC>0.0) ctAL=declength3d*1.115683/pV0MC;
	if(ptV0MC>0.0) ctTAL=declength*1.115683/ptV0MC;
	fHistPiAPDecayLengthResolution[0]->Fill(declength3d,declength);
	fHistPiAPDecayLengthVsPt[isSecd]->Fill(ptV0MC,declength);//(ptV0MC,ctAL);
	if(!fSetPtDepHist)	fHistPiAPDecayLengthVsCtau[isSecd]->Fill(massV0MC,ctAL);
	else	fHistPiAPDecayLengthVsCtau[isSecd]->Fill(ptV0MC,ctAL);
	fHistPiAPDecayLengthVsMass[isSecd]->Fill(massV0MC,ctTAL);//declength);
	//all V0s histo	   
	fHistArmenteros[isSecd]->Fill(alfa,qt);
	fHistV0RadiusZ[isSecd]->Fill(rMC2D,xyzMC[2]);
	fHistV0RadiusXY[isSecd]->Fill(xyzMC[0],xyzMC[1]);
	fHistV0RadiusZVSPt[isSecd]->Fill(ptV0MC,xyzMC[2]);
	fHistV0RadiusXYVSY[isSecd]->Fill(rapidity,rMC2D);
	fHistPiAPK0sVsALambdaMass->Fill(calcALambdamass,calcK0smass);
	// if(isSecd <1) fHistPiPPhiPosVsPtPosVsMass->Fill(massV0MC,ctTL,ptV0MC);//,ctK0s);//phiPosMC);//xxx
	//else fHistPiAPPhiPosVsPtPosVsMass->Fill(massV0MC,ctTL,ptV0MC);//,ctK0s);//phiPosMC);//xxx	  
	//fHistPiAPPhiPosVsPtPosVsMass->Fill(massV0MC,ctTL,ptV0MC);//,ctK0s);//phiPosMC);//xxx
	/*
	  fTHnFAL->Fill(valTHnMC);
	  fTHnFALDauEta->Fill(valTHnMCDauEta);
	  fTHnFALDauPhi->Fill(valTHnMCDauPhi);
	*/
      }
    }//MC true ana
    else{// V0 reco ana
      V0RecoLoop(id0,id1,isSecd,pdgCode,ptV0MC,pdgMother,ptXiMother,declength);
    }
      
  }//end MC stack loop

}
//________________________________________________________________________
void AliAnalysisTaskV0ForRAA::V0RecoLoop(Int_t id0,Int_t id1,Int_t isSecd,Int_t what,Double_t ptV0MC, Int_t pdgMother,Double_t ptXiMother,Double_t declengthV0MC){
  //loop over reconstructed particles

   
  //--------------------- define variables -----------------------//
  Double_t pp[3];
  Double_t pm[3];
  Double_t xr[3];
   
  Double_t massPi=0.13957018;
  Double_t massP=0.93827203;
  
  TLorentzVector positivesMIP;
  TLorentzVector negativesMIAP;
  TLorentzVector positivesMIPi;
  TLorentzVector negativesMIPi;

  Double_t magField = fESD->GetMagneticField();
  /*
    AliKFParticle::SetField(fESD->GetMagneticField());
    AliKFVertex primVtx(*(fESD->GetPrimaryVertex()));
    AliKFVertex primVtxImproved = primVtx;

    AliKFParticle* negPiKF=NULL;
    AliKFParticle* posPiKF=NULL;
    AliKFParticle* posPKF=NULL;
    AliKFParticle* negAPKF=NULL;
  */

  AliESDtrack* trackPos=NULL;
  AliESDtrack* trackNeg=NULL;
  AliESDtrack* trackPosTest = NULL;
  AliESDtrack* trackNegTest =NULL;

  Double_t primaryVtxPosition[3];
  primaryVtxPosition[0] = fESD->GetPrimaryVertex()->GetX();
  primaryVtxPosition[1] = fESD->GetPrimaryVertex()->GetY();
  primaryVtxPosition[2] = fESD->GetPrimaryVertex()->GetZ();
   
  Int_t nV0 = fESD->GetNumberOfV0s();
  const Int_t sizenV0 = nV0;
  AliESDv0 * v0MIs=NULL;


  //    Int_t on =0,off=0;
  Bool_t stopLoop = kFALSE;
  Int_t trackID[sizenV0][2];

  //----------- store V0 for daughter track id for position mapping ------------//
  Float_t v0idForDauPositionK0s[sizenV0];
  Float_t v0idForDauPositionL[sizenV0];
  Float_t v0idForDauPositionAL[sizenV0];
 
  //---------------------- for MC mode only ------------------//
  AliStack *stackRec = NULL;
  if(fMCMode && !fMCTruthMode) stackRec = fMCev->Stack();
    
  //------------------------ V0 reco loop --------------------//
  for(Int_t iV0MI = 0; iV0MI < nV0; iV0MI++) {//V0 loop
      
    v0idForDauPositionK0s[iV0MI] = 0.0;
    v0idForDauPositionL[iV0MI] = 0.0;
    v0idForDauPositionAL[iV0MI] = 0.0;

    //-- get V0 info --//
    v0MIs = fESD->GetV0(iV0MI);
    if(!v0MIs ) continue;

    fHistPiPiMonitorCuts->Fill(1);
    fHistPiPMonitorCuts[isSecd]->Fill(1);
    fHistPiAPMonitorCuts[isSecd]->Fill(1);

    if(stopLoop && fStopLoop) break;
    //------------ get references of daughters --------------//
    //-- esd tracks --//
    trackPosTest = fESD->GetTrack(v0MIs->GetPindex());
    trackNegTest = fESD->GetTrack(v0MIs->GetNindex());
     
    if ( trackPosTest->GetSign() == trackNegTest->GetSign()) continue;
	 
    fHistPiPiMonitorCuts->Fill(2);
    fHistPiPMonitorCuts[isSecd]->Fill(2);
    fHistPiAPMonitorCuts[isSecd]->Fill(2);

    //-- onthefly selection --//
    Bool_t onthefly = v0MIs->GetOnFlyStatus();
    if(fOntheFly!=onthefly) continue;
      
    fHistPiPiMonitorCuts->Fill(3);
    fHistPiPMonitorCuts[isSecd]->Fill(3);
    fHistPiAPMonitorCuts[isSecd]->Fill(3);
         

    Int_t indexPos = 0,indexNeg=0;
    indexPos = v0MIs->GetPindex();
    indexNeg = v0MIs->GetNindex();

    //-- for MC mode --//
    if(fMCMode){
      //check MC labels (and find partners for MC truth V0 daughters for fMCTruthMode=kTRUE)
      if(!GetMCTruthPartner(trackPosTest,trackNegTest,id0,id1)) continue;
      else stopLoop = kTRUE;
    }
    else{
      //check if V0 was alread found
      if(fStopLoop){
	if(CheckMultipleV0Candidates(indexPos,indexNeg,iV0MI,trackID)) continue;
      }
    }

    fHistPiPiMonitorCuts->Fill(4);
    fHistPiPMonitorCuts[isSecd]->Fill(4);
    fHistPiAPMonitorCuts[isSecd]->Fill(4);
      
         
    //--  get eta from V0 daughters --//
    Double_t posDaughterEta=0.0;
    Double_t negDaughterEta=0.0;
    Double_t posDaughterPhi=0.0;
    Double_t negDaughterPhi=0.0;
	 
    Double_t eta00 = trackPosTest->Eta();
    Double_t eta01 = trackNegTest->Eta();
            

  
    //---------- check sign assignment for daughters --------//
    Bool_t switchSign = kFALSE;
     
    if( trackPosTest->GetSign() >0){//pos
   
	
      v0MIs->GetPPxPyPz(pp[0],pp[1],pp[2]);
      v0MIs->GetNPxPyPz(pm[0],pm[1],pm[2]);

      posDaughterEta = v0MIs->GetParamP()->Eta();
      negDaughterEta = v0MIs->GetParamN()->Eta();
      posDaughterPhi = v0MIs->GetParamP()->Phi();
      negDaughterPhi = v0MIs->GetParamN()->Phi();
      /*      
	      if (negPiKF) delete negPiKF; negPiKF=NULL;
	      if (posPiKF) delete posPiKF; posPiKF=NULL;
	      if (posPKF) delete posPKF; posPKF=NULL;
	      if (negAPKF) delete negAPKF; negAPKF=NULL;

	      negPiKF = new AliKFParticle( *(v0MIs->GetParamN()) ,-211);
	      posPiKF = new AliKFParticle( *(v0MIs->GetParamP()) ,211);
	      posPKF = new AliKFParticle( *(v0MIs->GetParamP()) ,2212);
	      negAPKF = new AliKFParticle( *(v0MIs->GetParamN()) ,-2212);
      */
	    
    }
    if( trackPosTest->GetSign() <0){//neg

      indexPos = v0MIs->GetNindex();
      indexNeg = v0MIs->GetPindex();
	
      v0MIs->GetNPxPyPz(pp[0],pp[1],pp[2]);
      v0MIs->GetPPxPyPz(pm[0],pm[1],pm[2]);
      
      posDaughterEta = v0MIs->GetParamN()->Eta();
      negDaughterEta = v0MIs->GetParamP()->Eta();
      posDaughterPhi = v0MIs->GetParamN()->Phi();
      negDaughterPhi = v0MIs->GetParamP()->Phi();
      /*
	if (negPiKF) delete negPiKF; negPiKF=NULL;
	if (posPiKF) delete posPiKF; posPiKF=NULL;
	if (posPKF) delete posPKF; posPKF=NULL;
	if (negAPKF) delete negAPKF; negAPKF=NULL;


	negPiKF = new AliKFParticle( *(v0MIs->GetParamP()) ,-211);
	posPiKF = new AliKFParticle( *(v0MIs->GetParamN()) ,211);
	posPKF = new AliKFParticle( *(v0MIs->GetParamN()) ,2212);
	negAPKF = new AliKFParticle( *(v0MIs->GetParamP()) ,-2212);
      */
      
      switchSign = kTRUE;
      eta01 = trackPosTest->Eta();
      eta00 = trackNegTest->Eta();

    }
    
    trackPos =fESD->GetTrack(indexPos);
    trackNeg =fESD->GetTrack(indexNeg);
    
    // ------------- calc masses and 4 vectors -------------- //

    positivesMIPi.SetXYZM(pp[0],pp[1],pp[2],massPi);
    negativesMIPi.SetXYZM(pm[0],pm[1],pm[2],massPi);
    positivesMIP.SetXYZM(pp[0],pp[1],pp[2],massP);
    negativesMIAP.SetXYZM(pm[0],pm[1],pm[2],massP);
    //  cout<<positivesMIPi.Pt()<<"  "<< positivesMIPi.Px()<<"  "<< positivesMIPi.Py()<<"  "<< positivesMIPi.Pz()<<"  "<< positivesMIPi.M()<<endl;
    //  //  cout<<negativesMIPi.Pt()<<"  "<< negativesMIPi.Px()<<"  "<< negativesMIPi.Py()<<"  "<< negativesMIPi.Pz()<<"  "<< negativesMIPi.M()<<endl;
    if(fShift && fabs(fDeltaInvP) >0.0){//shift in pT by 1/deltaPt
      //pos pion
      Double_t getPtPos = positivesMIPi.Pt();
      Double_t shiftPtPos = 0.0;
      if(getPtPos >0.0){
	shiftPtPos = 1.0/getPtPos + fDeltaInvP;
	if(fabs(shiftPtPos)>0.0) getPtPos  = 1.0/shiftPtPos;
      }
      Double_t getPzPos = positivesMIPi.Pz();
      Double_t calcPtotPos = sqrt(pow(getPtPos,2.0)+pow(getPzPos,2.0));
      Double_t calcEtaPos = 0.5*log((calcPtotPos + getPzPos)/(calcPtotPos-getPzPos));
           
      //neg pion
      Double_t getPtNeg = negativesMIPi.Pt();
      Double_t shiftPtNeg = 0.0;
      if(getPtNeg >0.0){
	shiftPtNeg = -1.0/getPtNeg + fDeltaInvP;
	if(fabs(shiftPtNeg)>0.0) getPtNeg  = 1.0/shiftPtNeg;
      }
      Double_t getPzNeg = negativesMIPi.Pz();
      Double_t calcPtotNeg = sqrt(pow(getPtNeg,2.0)+pow(getPzNeg,2.0));
      Double_t calcEtaNeg = 0.5*log((calcPtotNeg + getPzNeg)/(calcPtotNeg-getPzNeg));

      //set new pt
      positivesMIPi.SetPtEtaPhiM(fabs(getPtPos),calcEtaPos,positivesMIPi.Phi(),massPi);
      negativesMIPi.SetPtEtaPhiM(fabs(getPtNeg),calcEtaNeg,negativesMIPi.Phi(),massPi);
      positivesMIP.SetPtEtaPhiM(fabs(getPtPos),calcEtaPos,positivesMIPi.Phi(),massP);
      negativesMIAP.SetPtEtaPhiM(fabs(getPtNeg),calcEtaNeg,negativesMIPi.Phi(),massP);     
      if(getPtPos < 0.0 && getPtNeg > 0.0) {
	positivesMIPi.SetPtEtaPhiM(fabs(getPtNeg),calcEtaNeg,negativesMIPi.Phi(),massPi);
	negativesMIPi.SetPtEtaPhiM(fabs(getPtPos),calcEtaPos,positivesMIPi.Phi(),massPi);
	positivesMIP.SetPtEtaPhiM(fabs(getPtNeg),calcEtaNeg,negativesMIPi.Phi(),massP);
	negativesMIAP.SetPtEtaPhiM(fabs(getPtPos),calcEtaPos,positivesMIPi.Phi(),massP);    
      }
      if((getPtPos < 0.0 && getPtNeg< 0.0) ||(getPtPos > 0.0 && getPtNeg > 0.0)) continue;
      // cout<<"*************"<<endl;
      //  //  cout<<positivesMIPi.Pt()<<"  "<< positivesMIPi.Px()<<"  "<< positivesMIPi.Py()<<"  "<< positivesMIPi.Pz()<<"  "<< positivesMIPi.M()<<endl;
      //  cout<<negativesMIPi.Pt()<<"  "<< negativesMIPi.Px()<<"  "<< negativesMIPi.Py()<<"  "<< negativesMIPi.Pz()<<"  "<< negativesMIPi.M()<<endl;
    }
    // cout<<"------------------------------*************"<<endl;
    //K0
    TLorentzVector v0K0=positivesMIPi+negativesMIPi;
    //Lambda
    TLorentzVector v0Lambda=positivesMIP+negativesMIPi;
    //Anitlambda
    TLorentzVector v0ALambda=positivesMIPi+negativesMIAP;

    //---------------------AliKFParticle ---------------------//
    /*  
	Double_t chi2K0C=0.0;
	Double_t chi2LambdaC=0.0;
	Double_t chi2ALambdaC=0.0;

     
	AliKFParticle v0K0KF;
	v0K0KF +=(*negPiKF);
	v0K0KF +=(*posPiKF);
	//v0K0C.SetVtxGuess(xr[0],xr[1],xr[2]);
	v0K0KF.SetProductionVertex(primVtxImproved);
	  
	AliKFParticle v0LambdaKF;
	v0LambdaKF +=(*negPiKF);
	v0LambdaKF +=(*posPKF);
	//v0LambdaC.SetVtxGuess(xr[0],xr[1],xr[2]);
	v0LambdaKF.SetProductionVertex(primVtxImproved);
	  
	AliKFParticle v0ALambdaKF;
	v0ALambdaKF +=(*negAPKF);
	v0ALambdaKF +=(*posPiKF);
	//v0ALambdaC.SetVtxGuess(xr[0],xr[1],xr[2]);
	v0ALambdaKF.SetProductionVertex(primVtxImproved);
	
	if( v0K0KF.GetNDF() != 0) {
	chi2K0C = v0K0KF.GetChi2()/v0K0KF.GetNDF();
	}

	Double_t chi2LambdaC=100000.;
	if( v0LambdaKF.GetNDF() != 0) {
	chi2LambdaC = v0LambdaKF.GetChi2()/v0LambdaKF.GetNDF();
	}

	Double_t chi2ALambdaC=100000.;
	if( v0ALambdaKF.GetNDF() != 0) {
	chi2ALambdaC = v0ALambdaKF.GetChi2()/v0ALambdaKF.GetNDF();
	}
    */
      
    // ----------------- for MC mode ------------------------ //
    Bool_t fillK0sMC = kTRUE;
    Bool_t fillLambdaMC = kTRUE;
    Bool_t fillALambdaMC = kTRUE;

    if(fMCMode && fMCTruthMode) {
      if(what == 310) {
	fillLambdaMC = kFALSE;
	fillALambdaMC = kFALSE;
      }
      else if(what == 3122){
	fillALambdaMC = kFALSE;
	fillK0sMC = kFALSE;
      }
      else if(what == -3122){
	fillLambdaMC = kFALSE;
	fillK0sMC = kFALSE;
      }
    }
   
    //----------------- prepare for V0 ana ------------------//
    TVector3 ppTrack(pp);
    TVector3 pmTrack(pm);
      
    //-- momenta --//
    Double_t ptK0s = v0K0.Pt();
    Double_t ptLambda = v0Lambda.Pt();
    Double_t ptALambda = v0ALambda.Pt();
      
    Double_t pK0s = v0K0.P();
    Double_t pLambda = v0Lambda.P();
    Double_t pALambda = v0ALambda.P();

    Double_t posDaughterP = ppTrack.Mag();
    Double_t negDaughterP = pmTrack.Mag();

    v0MIs->GetXYZ(xr[0],xr[1],xr[2]);  

    // Double_t posDaughterPt = ppTrack.Pt();
    // Double_t negDaughterPt = pmTrack.Pt();
 
    /*
      Double_t v0sPt=v0MIs->Pt();
      if(what == 310 || what ==0){
      fHistPiPiEtaDReco[0]->Fill(posDaughterPt,v0sPt);
      fHistPiPiEtaDReco[0]->Fill(negDaughterPt,v0sPt);
      }
      if(fabs(what) == 3122 || what == 0){
      fHistPiPEtaDReco[0]->Fill(posDaughterPt,v0sPt);
      fHistPiPEtaDReco[0]->Fill(negDaughterPt,v0sPt);
      }
    */
     
     

    //--------------------------------------------------------- general cuts --------------------------------------------------------------//
    //-- track cuts for daughters --//
    //-- eta cut --//
    if( fabs(posDaughterEta) > fEtaCutMCDaughtersVal || fabs(negDaughterEta) > fEtaCutMCDaughtersVal) continue;
    fHistPiPiMonitorCuts->Fill(5);
    fHistPiPMonitorCuts[isSecd]->Fill(5);
    fHistPiAPMonitorCuts[isSecd]->Fill(5);

    //-- esd track cuts --//
    //K0s
    if( ptK0s > fPtTPCCut){
      if(fESDTrackCuts){
	if(!fESDTrackCuts->AcceptTrack(trackPosTest) || !fESDTrackCuts->AcceptTrack(trackNegTest)) continue;
	else  fHistPiPiMonitorCuts->Fill(6); 
      }
    }
    else{
      if(fESDTrackCutsLowPt){
	if(!fESDTrackCutsLowPt->AcceptTrack(trackPosTest) || !fESDTrackCutsLowPt->AcceptTrack(trackNegTest))  continue;
      }
    }
    //Lambda
    if(ptLambda > fPtTPCCut){
      if(fESDTrackCuts && fESDTrackCutsCharged){
	if(!fESDTrackCutsCharged->AcceptTrack(trackPosTest) || !fESDTrackCuts->AcceptTrack(trackNegTest)) continue;
	else  fHistPiPMonitorCuts[isSecd]->Fill(6); 
      }
    }
    else{
      if(fESDTrackCutsLowPt){
	if(!fESDTrackCutsLowPt->AcceptTrack(trackPosTest) || !fESDTrackCutsLowPt->AcceptTrack(trackNegTest))  continue;
      }
    }
    //ALambda
    if(ptALambda > fPtTPCCut){
      if(fESDTrackCuts && fESDTrackCutsCharged){
	if(!fESDTrackCuts->AcceptTrack(trackPosTest) || !fESDTrackCutsCharged->AcceptTrack(trackNegTest)) continue;
	else  fHistPiAPMonitorCuts[isSecd]->Fill(6); 
      }
    }
    else{
      if(fESDTrackCutsLowPt){
	if(!fESDTrackCutsLowPt->AcceptTrack(trackPosTest)|| !fESDTrackCutsLowPt->AcceptTrack(trackNegTest))  continue;
      }
    }

      
    //------------------------ detector values -------------------------------//
    //-- TPC ITS values pos --//
    Int_t nclsTPCPos =  trackPos->GetNcls(1);
    Int_t nclsTPCFindablePos =  trackPos->GetTPCNclsF();
    Int_t nclsITSPos =  trackPos->GetNcls(0);
    Double_t chi2PerClusterITSPos = -1.0;
    if(nclsITSPos>0) chi2PerClusterITSPos = trackPos->GetITSchi2()/Double_t(nclsITSPos);
    Double_t crossedRowsTPCPos = trackPos->GetTPCCrossedRows();
      
    //-- TPC ITS values neg --//
    Int_t nclsTPCNeg =  trackNeg->GetNcls(1);
    Int_t nclsTPCFindableNeg =  trackNeg->GetTPCNclsF();
    Int_t nclsITSNeg =  trackNeg->GetNcls(0);
    Double_t chi2PerClusterITSNeg = -1.0;
    if(nclsITSNeg>0) chi2PerClusterITSNeg =trackNeg->GetITSchi2()/Double_t(nclsITSNeg);
    Double_t crossedRowsTPCNeg = trackNeg->GetTPCCrossedRows();    

    Double_t ratio = 10.0;
    if(nclsTPCFindableNeg >0.0) ratio =double(crossedRowsTPCNeg)/ double(nclsTPCFindableNeg);
    
    Double_t ratioPos = 10.0;
    if(nclsTPCFindablePos >0.0) ratioPos =double(crossedRowsTPCPos)/ double(nclsTPCFindablePos);
    
    Double_t ratioFoFi = 10.0;
    if(nclsTPCFindableNeg >0.0) ratioFoFi =double(nclsTPCNeg)/ double(nclsTPCFindableNeg);
    
    Double_t ratioFoFiPos = 10.0;
    if(nclsTPCFindablePos >0.0) ratioFoFiPos =double(nclsTPCPos)/ double(nclsTPCFindablePos);

    //track length TPC cut
    Double_t lengthTPCPos = trackPos->GetLengthInActiveZone(0,3,236, magField,0,0);//-5 ,0,0);
    Double_t lengthTPCNeg = trackNeg->GetLengthInActiveZone(0,3,236, magField,0,0);//-5 ,0,0);
    if(fCutMITrackLength && lengthTPCPos <=  fCutMITrackLengthLengthF * (130.0 - 5.0*fabs(trackPos->GetSigned1Pt()))) continue;
    fHistPiPiMonitorCuts->Fill(7);
    fHistPiPMonitorCuts[isSecd]->Fill(7);
    fHistPiAPMonitorCuts[isSecd]->Fill(7); 
    if(fCutMITrackLength && lengthTPCNeg <=  fCutMITrackLengthLengthF * (130.0 - 5.0*fabs(trackNeg->GetSigned1Pt()))) continue;
    fHistPiPiMonitorCuts->Fill(8);
    fHistPiPMonitorCuts[isSecd]->Fill(8);
    fHistPiAPMonitorCuts[isSecd]->Fill(8); 

    //crossed rows TPC cut
    if(fCutMICrossedR && trackPos->GetTPCClusterInfo(3,1) <= fCutMICrossedRLengthF *(130.0 - 5.0*fabs(trackPos->GetSigned1Pt()))) continue;
    fHistPiPiMonitorCuts->Fill(9);
    fHistPiPMonitorCuts[isSecd]->Fill(9);
    fHistPiAPMonitorCuts[isSecd]->Fill(9); 
    if(fCutMICrossedR && trackNeg->GetTPCClusterInfo(3,1) <= fCutMICrossedRLengthF *(130.0 - 5.0*fabs(trackNeg->GetSigned1Pt()))) continue;
    fHistPiPiMonitorCuts->Fill(10);
    fHistPiPMonitorCuts[isSecd]->Fill(10);
    fHistPiAPMonitorCuts[isSecd]->Fill(10);

    // ncls TPC cut
    if(fCutMITPCncls &&  nclsTPCPos <= 0.6*(130.0 - 5.0*fabs(trackPos->GetSigned1Pt()))) continue;
    fHistPiPiMonitorCuts->Fill(11);
    fHistPiPMonitorCuts[isSecd]->Fill(11);
    fHistPiAPMonitorCuts[isSecd]->Fill(11); 
    if(fCutMITPCncls &&   nclsTPCNeg <= 0.6*(130.0 - 5.0*fabs(trackNeg->GetSigned1Pt()))) continue;
    fHistPiPiMonitorCuts->Fill(12);
    fHistPiPMonitorCuts[isSecd]->Fill(12);
    fHistPiAPMonitorCuts[isSecd]->Fill(12); 

      
      
    //found
    if(fMoreNclsThanRows && (crossedRowsTPCPos < nclsTPCPos || crossedRowsTPCNeg < nclsTPCNeg  )) continue;
    fHistPiPiMonitorCuts->Fill(13);
    fHistPiPMonitorCuts[isSecd]->Fill(13);
    fHistPiAPMonitorCuts[isSecd]->Fill(13);
      
    if(fMoreNclsThanFindable && (nclsTPCFindablePos < nclsTPCPos || nclsTPCFindableNeg < nclsTPCNeg  )) continue;
    fHistPiPiMonitorCuts->Fill(14);
    fHistPiPMonitorCuts[isSecd]->Fill(14);
    fHistPiAPMonitorCuts[isSecd]->Fill(14);      

    if(fMoreNclsThanFindableMax && ( nclsTPCPos < (nclsTPCFindablePos -60.0) || nclsTPCNeg < (nclsTPCFindableNeg - 60.0) )) continue;
    // if(chi2PerClusterITSNeg > fChi2PerClusterITS || chi2PerClusterITSPos > fChi2PerClusterITS ) continue;
    fHistPiPiMonitorCuts->Fill(15);
    fHistPiPMonitorCuts[isSecd]->Fill(15);
    fHistPiAPMonitorCuts[isSecd]->Fill(15);  
      
    if(ratio > fRatioMaxCRowsOverFindable || ratioPos > fRatioMaxCRowsOverFindable) continue; 

    if(ratioFoFi < fRatioFoundOverFindable || ratioFoFiPos < fRatioFoundOverFindable) continue; 
    fHistPiPiMonitorCuts->Fill(16);
    fHistPiPMonitorCuts[isSecd]->Fill(16);
    fHistPiAPMonitorCuts[isSecd]->Fill(16); 

    Bool_t cutOKITSNegNeg =kTRUE;
    Bool_t cutOKITSPosPos =kTRUE;

    Bool_t cutOKITSNegPos =kTRUE;
    Bool_t cutOKITSPosNeg =kTRUE;

    if(nclsITSNeg < fMinNCLSITSNeg ||  nclsITSNeg > fMaxNCLSITSNeg){
      if(!fSwitchCaseITSCls) continue;
      else cutOKITSNegNeg = kFALSE;
    }
  
    fHistPiPiMonitorCuts->Fill(17);
    fHistPiPMonitorCuts[isSecd]->Fill(17);
    fHistPiAPMonitorCuts[isSecd]->Fill(17); 

    //2D decay radius of V0
    Double_t dim2V0Radius= sqrt(   pow(xr[0] - primaryVtxPosition[0],2.0)
				   +pow(xr[1] - primaryVtxPosition[1],2.0));//TMath::Sqrt(xr[0]*xr[0]+xr[1]*xr[1]);

    if(nclsITSPos < fMinNCLSITSPos || nclsITSPos > fMaxNCLSITSPos){
      if(dim2V0Radius >= fDecRadCutITSMin  && dim2V0Radius < fDecRadCutITSMax){//only for a certain decay radius 2D in xy
	if(!fSwitchCaseITSCls) continue;
	else cutOKITSPosPos = kFALSE;
      }
    }
    
    if(fSwitchCaseITSCls){
      if(nclsITSPos < fMinNCLSITSNeg || nclsITSPos > fMaxNCLSITSNeg){
	if(dim2V0Radius >= fDecRadCutITSMin  && dim2V0Radius < fDecRadCutITSMax){//only for a certain decay radius 2D in xy
	  cutOKITSPosNeg = kFALSE;
	}
      }
      if(nclsITSNeg < fMinNCLSITSPos || nclsITSNeg > fMaxNCLSITSPos){
	if(dim2V0Radius >= fDecRadCutITSMin  && dim2V0Radius < fDecRadCutITSMax){//only for a certain decay radius 2D in xy	 
	  cutOKITSNegPos = kFALSE;
	}      
      }
      
      if(!(cutOKITSNegPos && cutOKITSPosNeg) && !(cutOKITSNegNeg && cutOKITSPosPos) ) continue;
    }

    fHistPiPiMonitorCuts->Fill(18);
    fHistPiPMonitorCuts[isSecd]->Fill(18);
    fHistPiAPMonitorCuts[isSecd]->Fill(18); 

    
     
    
   
    //--------------------- PID ----------------------------//
    //-- dEdx --//
    Float_t nSigmaTPCtrackPosToPion = 0.0;
    Float_t nSigmaTPCtrackNegToPion = 0.0;
    Float_t nSigmaTPCtrackPosToProton = 0.0;
    Float_t nSigmaTPCtrackNegToProton = 0.0;

	 
    if(fESDpid){
      nSigmaTPCtrackPosToPion = fabs(fESDpid->NumberOfSigmasTPC(trackPos,AliPID::kPion));
      nSigmaTPCtrackNegToPion = fabs(fESDpid->NumberOfSigmasTPC(trackNeg,AliPID::kPion));
      nSigmaTPCtrackPosToProton = fabs(fESDpid->NumberOfSigmasTPC(trackPos,AliPID::kProton));
      nSigmaTPCtrackNegToProton = fabs(fESDpid->NumberOfSigmasTPC(trackNeg,AliPID::kProton));
    }
	 
    Bool_t pipidEdx=kTRUE;
    Bool_t pipdEdx =kTRUE;
    Bool_t piapdEdx=kTRUE;

    Double_t tpcsigPos= trackPos->GetTPCsignal();
    Double_t tpcsigNeg= trackNeg->GetTPCsignal();
     
    /*
      Double_t tpcsigNPos= trackPos->GetTPCsignalN();
      Double_t tpcsigNNeg= trackNeg->GetTPCsignalN();
    */
    //     GetYAt(Double_t x, Double_t b, Double_t &y) or GetY()
    Double_t posY =  trackPos->GetInnerParam()->GetY();
    Double_t posZ =  trackPos->GetInnerParam()->GetZ();
    Double_t negY =  trackNeg->GetInnerParam()->GetY();
    Double_t negZ =  trackNeg->GetInnerParam()->GetZ();
    Double_t distTPCinner  = sqrt(pow((posY-negY),2.0)+pow((posZ-negZ),2.0));
    if(distTPCinner < fDistanceTPCInner) continue;
    fHistPiPiMonitorCuts->Fill(19);
    fHistPiPMonitorCuts[isSecd]->Fill(19);
    fHistPiAPMonitorCuts[isSecd]->Fill(19); 

    //AliExternalTrackParam *extTParPos = (AliExternalTrackParam*)trackPos->GetTPCInnerParam();
    //Double_t tpcMomPos = extTParPos->GetP();
    Double_t tpcMomPos = trackPos->GetInnerParam()->GetP();
    // AliExternalTrackParam *extTParNeg = (AliExternalTrackParam*)trackNeg->GetTPCInnerParam();
    // Double_t tpcMomNeg = extTParNeg->GetP();
    Double_t tpcMomNeg = trackNeg->GetInnerParam()->GetP();
	 
    //-- dedx cut --//
    if(fUsePID){
      if(fabs(posDaughterP)<fPPIDcut && tpcsigPos < 5.0){//no zero dedx values!
	pipidEdx =kFALSE;//k0s
	piapdEdx =kFALSE;//antilambda
      }

      if(fabs(negDaughterP)<fPPIDcut &&  tpcsigNeg < 5.0){//no zero dedx values!
	pipidEdx =kFALSE;//k0s
	pipdEdx =kFALSE;//lambda
      }

      if(fabs(posDaughterP)<fPPIDcut && (nSigmaTPCtrackPosToProton > fNSigma || tpcsigPos < 5.0)) pipdEdx =kFALSE;//lambda
      if(fabs(negDaughterP)<fPPIDcut && (nSigmaTPCtrackNegToProton > fNSigma || tpcsigNeg < 5.0)) piapdEdx =kFALSE;//antilambda
     
      if(fabs(fNSigma-fNSigma2) > 0.001){
	if(fabs(posDaughterP) >= fPPIDcut && (nSigmaTPCtrackPosToProton > fNSigma2 || tpcsigPos < 5.0)) pipdEdx =kFALSE;//lambda
	if(fabs(negDaughterP) >= fPPIDcut && (nSigmaTPCtrackNegToProton > fNSigma2 || tpcsigNeg < 5.0)) piapdEdx =kFALSE;//antilambda

	if(fabs(posDaughterP) >= fPPIDcut && tpcsigPos < 5.0){//no zero dedx values!
	  pipidEdx =kFALSE;//k0s
	  piapdEdx =kFALSE;//antilambda
	}
	
	if(fabs(negDaughterP) >= fPPIDcut &&  tpcsigNeg < 5.0){//no zero dedx values!
	  pipidEdx =kFALSE;//k0s
	  pipdEdx =kFALSE;//lambda
	}
	
      }
      
    }


    if(fUsePIDPion){
      if(fabs(posDaughterP)<fPPIDcut && nSigmaTPCtrackPosToPion > fNSigma ){
	pipidEdx =kFALSE;//k0s
      }
      
      if(fabs(negDaughterP)<fPPIDcut && nSigmaTPCtrackNegToPion > fNSigma ){
	pipidEdx =kFALSE;//k0s
      }
    }



    //------------------- DCA  ---------------------//
      
    //-- between the daughters --//
    Double_t dcaDaughters = v0MIs->GetDcaV0Daughters();  
    
    //-- to primary vertex --//
    /* 
       Float_t bP[2],bN[2];
       Float_t bCovP[3],bCovN[3];
    
       trackPos->GetImpactParameters(bP,bCovP);
       trackNeg->GetImpactParameters(bN,bCovN);
    
       if (bCovP[0]<=0 || bCovP[2]<=0) {
       AliDebug(1, "Estimated b resolution lower or equal zero!");
       bCovP[0]=0; bCovP[2]=0;
       }
       if (bCovN[0]<=0 || bCovN[2]<=0) {
       AliDebug(1, "Estimated b resolution lower or equal zero!");
       bCovN[0]=0; bCovN[2]=0;
       }
    
       Float_t dcaToVertexZPos = bP[1];//Float_t dcaToVertexXY = b[0];
       Float_t dcaToVertexZNeg = bN[1];//Float_t dcaToVertexXY = b[0];    
    */
    /*  
	Float_t dcaToVertexZPos = 0.0, dcaToVertexZNeg = 0.0;
	Float_t bP=0.0,bN=0.0;
	trackPos->GetImpactParameters(bP,dcaToVertexZPos);
	trackNeg->GetImpactParameters(bN,dcaToVertexZNeg);
    */

    Double_t dcaToVertexZPos = 0.0, dcaToVertexZNeg = 0.0;
    AliExternalTrackParam *parPos = NULL;
    AliExternalTrackParam *parNeg = NULL;
    Double_t dcaYZP[2],dcaYZN[2],covar[3];
    if(!switchSign){
      parPos = new AliExternalTrackParam( *v0MIs->GetParamP());
      parNeg = new AliExternalTrackParam( *v0MIs->GetParamN());
    }
    else{
      parPos = new AliExternalTrackParam( *v0MIs->GetParamN());
      parNeg = new AliExternalTrackParam( *v0MIs->GetParamP());
    }
    Bool_t checkProp = parPos->PropagateToDCA(fESD->GetPrimaryVertex(),magField,40.0,dcaYZP,covar);
    dcaToVertexZPos =  dcaYZP[1];
    delete parPos;
    checkProp = parNeg->PropagateToDCA(fESD->GetPrimaryVertex(),magField,40.0,dcaYZN,covar);
    dcaToVertexZNeg =  dcaYZN[1];
    delete parNeg;


    Double_t dcaPosToVertex=0.0,dcaNegToVertex=0.0;
    Double_t dzPos=(primaryVtxPosition[0]-xr[0])*ppTrack.Y() - (primaryVtxPosition[1]-xr[1])*ppTrack.X();
    dcaPosToVertex=TMath::Sqrt(dzPos*dzPos/(pow(ppTrack.X(),2)+pow(ppTrack.Y(),2)));
    Double_t dzNeg=(primaryVtxPosition[0]-xr[0])*pmTrack.Y() - (primaryVtxPosition[1]-xr[1])*pmTrack.X();
    dcaNegToVertex=TMath::Sqrt(dzNeg*dzNeg/(pow(pmTrack.X(),2)+pow(pmTrack.Y(),2)));
     
    // Double_t dcaPosToVertex[3];dcaNegToVertex[3];
    // trackPos->GetImpactParameters(dcaPosToVertex[0],dcaPosToVertex[1]);
    // trackNeg->GetImpactParameters(dcaNegToVertex[0],dcaNegToVertex[1]);
	 
    // dcaPosToVertex = TMath::Sqrt(dcaPosToVertex[0]*dcaPosToVertex[0]+dcaPosToVertex[1]*dcaPosToVertex[1]);
    // dcaNegToVertex = TMath::Sqrt(dcaNegToVertex[0]*dcaNegToVertex[0]+dcaNegToVertex[1]*dcaNegToVertex[1]);
	 
    // dcaPosToVertex  =   posPKF->GetDistanceFromVertexXY(primaryVtxPosition);
    // dcaNegToVertex  =   negPiKF->GetDistanceFromVertexXY(primaryVtxPosition);

    Double_t  dcaV0ToPrimVertex= v0MIs->GetD(primaryVtxPosition[0],primaryVtxPosition[1]);////v0K0KF.GetDistanceFromVertexXY(tPrimaryVtxPosition); 

    Double_t dcaZ=(primaryVtxPosition[0]-xr[0])*v0K0.Y() - (primaryVtxPosition[1]-xr[1])*v0K0.X();
    Double_t dcaZToVertex=TMath::Sqrt(dcaZ*dcaZ/(pow(v0K0.X(),2)+pow(v0K0.Y(),2)));
        
    //------------------- decay length V0 -------------//
      

    Double_t decayLength = sqrt( pow(xr[0] - primaryVtxPosition[0],2.0) 
				 +pow(xr[1] - primaryVtxPosition[1],2.0)
				 +pow(xr[2] - primaryVtxPosition[2],2.0)
				 );
    //2D decay radius already calculated for track length cut

    //-- decay radius xy min cut --//
    if(dim2V0Radius < fDecayRadXYMin && ptK0s < fPtDecRadMin) continue;
    //	    if(fabs(xr[1])<fDecayRadY) continue;
    fHistPiPiMonitorCuts->Fill(20);
    fHistPiPMonitorCuts[isSecd]->Fill(20);
    fHistPiAPMonitorCuts[isSecd]->Fill(20);

    //-- decay radius xy max cut --//
    if(dim2V0Radius > fDecayRadXYMax) continue;
    //	    if(fabs(xr[1])<fDecayRadY) continue;
    fHistPiPiMonitorCuts->Fill(21);
    fHistPiPMonitorCuts[isSecd]->Fill(21);
    fHistPiAPMonitorCuts[isSecd]->Fill(21);
      
    //-- 3D decay length min ->ctau --//
    if(decayLength > fDecayLengthMax) continue;
    fHistPiPiMonitorCuts->Fill(22);
    fHistPiPMonitorCuts[isSecd]->Fill(22);
    fHistPiAPMonitorCuts[isSecd]->Fill(22);
	 
    //-- 3D decay length min cut --//
    if(decayLength < fDecayLengthMin) continue;
    fHistPiPiMonitorCuts->Fill(23);
    fHistPiPMonitorCuts[isSecd]->Fill(23);
    fHistPiAPMonitorCuts[isSecd]->Fill(23);
   
   

    //----------------------- V0 variables --------------------//
    //-- armenteros --//
    TVector3 momTot = ppTrack + pmTrack;
    Double_t lQlNeg = fabs(pmTrack.Dot(momTot)/momTot.Mag());
    Double_t lQlPos = fabs(ppTrack.Dot(momTot)/momTot.Mag());
    //return 1.-2./(1.+lQlNeg/lQlPos);
    Double_t alfa =0.0;
    Double_t den = lQlPos + lQlNeg;
    if(den>0) alfa = (lQlPos - lQlNeg)/den;
    TVector3 qtvec= pmTrack.Cross(momTot);//vecPip.Mag()*sqrt(1-pow(thetapip,2));
    Double_t qt = qtvec.Mag()/momTot.Mag();

      
      
    //-- masses --//
    Double_t massK0s = v0K0.M();
    Double_t massLambda = v0Lambda.M();
    Double_t massALambda = v0ALambda.M();

    Double_t energyE1 = sqrt(ppTrack.Mag2()+pow(0.51099e-03,2.0));
    Double_t energyE2 = sqrt(pmTrack.Mag2()+pow(0.51099e-03,2.0));
    TLorentzVector e1(ppTrack,energyE1);
    TLorentzVector e2(pmTrack,energyE2);
    TLorentzVector photon = e1+e2;
    Double_t massPhoton = photon.M();
     
    //-- rapidity --//
    Double_t rapK0s = v0MIs->Y(310);
    Double_t rapL   = v0MIs->Y(3122);
    Double_t rapAL  = v0MIs->Y(3122);

    //-- other variables --//
    Double_t opAng =   fabs(ppTrack.Angle(pmTrack));
    Double_t cosOPAng = v0MIs->GetV0CosineOfPointingAngle();
      
    //    if( ppTrack.Angle(pmTrack)<0.001) continue;  
    //    if( ppTrack.Angle(pmTrack)<0.004) continue;   
    /*    
	  Double_t px = v0K0.Px();
	  Double_t py = v0K0.Py();
	  Double_t phi  = TMath::Pi()+TMath::ATan2(-py, -px);
	  
    */
    Double_t eta =  v0K0.Eta();
    /*     
    //introduce more histo
    Double_t errOnMassK0s = v0MIs->ChangeMassHypothesis(310);
    Double_t errOnMassLambda = 0.0;
    Double_t errOnMassALambda = 0.0;
    if(!switchSign){
    errOnMassLambda  = v0MIs->ChangeMassHypothesis(3122);
    errOnMassALambda = v0MIs->ChangeMassHypothesis(-3122);
    }
    else{
    errOnMassLambda  = v0MIs->ChangeMassHypothesis(-3122);
    errOnMassALambda = v0MIs->ChangeMassHypothesis(3122);
    }
    */

    //------------------ cut flags for V0 type specific cuts --------------//
    Bool_t cutOKK0s = kTRUE;
    Bool_t cutOKLambda = kTRUE;
    Bool_t cutOKALambda = kTRUE;

    //-------------------------- K0 cuts -----------------------------//

    if(dcaV0ToPrimVertex > fDCAToVertexK0)  cutOKK0s = kFALSE;
    else fHistPiPiMonitorCuts->Fill(24);
      
    if(fabs(dcaToVertexZNeg) > fDCAZ  || fabs(dcaToVertexZPos) > fDCAZ ) cutOKK0s = kFALSE;
    //    if(fabs(xr[2])> fDCAZ) cutOKK0s = kFALSE; //like decay radius z component
    else fHistPiPiMonitorCuts->Fill(25);
      
    Double_t ctK0 = 0.0,ctTK0 = 0.0;
    if(fabs(pK0s)>0.0)  ctK0 = decayLength*0.497614/pK0s;
    if(fabs(ptK0s)>0.0)  ctTK0 = dim2V0Radius*0.497614/ptK0s;
    if(ctK0 > fCtauK0s &&  fabs(ptK0s) <fCtauPtCutK0) cutOKK0s = kFALSE;
    else  fHistPiPiMonitorCuts->Fill(26);
      
    if((cosOPAng < fCosPointAngK && fabs(ptK0s) > fCPAPtCutK0)|| cosOPAng<0.99)///xxx
      cutOKK0s = kFALSE;
    else  fHistPiPiMonitorCuts->Fill(27);

    if(dcaDaughters > fDCADaughtersK0 )cutOKK0s = kFALSE;
    else  fHistPiPiMonitorCuts->Fill(28);
	 
    if(dcaNegToVertex < fDCADaughtersToVtxSmall || dcaPosToVertex < fDCADaughtersToVtxSmall)  cutOKK0s = kFALSE;
    else  fHistPiPiMonitorCuts->Fill(29);

    if(fRapCutV0 && fabs(rapK0s) > fRap) cutOKK0s = kFALSE;
    else  fHistPiPiMonitorCuts->Fill(30);  
    
    // if(chi2K0C > fChiCutKf) cutOKK0s = kFALSE;
    if(opAng < fOpengAngleDaughters && fabs(ptK0s) < fOpAngPtCut )  cutOKK0s = kFALSE;
    else fHistPiPiMonitorCuts->Fill(31);
    
    Bool_t ptbinokK0s=kFALSE;
    if( ptK0s < fQtCutPt &&  ptK0s > fQtCutPtLow ) ptbinokK0s=kTRUE;
    
    Double_t qtval = fArmQtSlope*fabs(alfa);
   
    if(fArmCutK0 && ptbinokK0s && qt < qtval) cutOKK0s = kFALSE;
    else  fHistPiPiMonitorCuts->Fill(32);     
    if(fArmCutK0 && ptbinokK0s && qt < fQtCut) cutOKK0s = kFALSE;
      
    if(fEtaSignCut *eta < 0.0) cutOKK0s = kFALSE; 
     
    //-------------------------- Lambda cuts -------------------------//

    if(dcaV0ToPrimVertex > fDCAToVertexL) cutOKLambda = kFALSE;
    else  fHistPiPMonitorCuts[isSecd]->Fill(24);

    if(fabs(dcaToVertexZNeg) > fDCAZ  || fabs(dcaToVertexZPos) > fDCAZ ) cutOKLambda = kFALSE;
    //    if(fabs(xr[2])>fDCAZ) cutOKLambda = kFALSE; //like decay radius z component
    else  fHistPiPMonitorCuts[isSecd]->Fill(25);
         
    Double_t ctL = 0.0,ctTL=0.0;
    if(fabs(pLambda)>0.0)  ctL  = decayLength*1.115683/fabs(pLambda);
    if(fabs(ptLambda)>0.0) ctTL = dim2V0Radius*1.115683/fabs(ptLambda);
	 
    if(ctL > fCtauL && fabs(ptLambda) <fCtauPtCutL)  cutOKLambda = kFALSE;
    else  fHistPiPMonitorCuts[isSecd]->Fill(26);
      
    if((cosOPAng<fCosPointAngL && fabs(ptLambda) > fCPAPtCutL)|| cosOPAng<0.99)///xxx
      cutOKLambda = kFALSE;
    else fHistPiPMonitorCuts[isSecd]->Fill(27);

    if(dcaDaughters > fDCADaughtersL )cutOKLambda = kFALSE;
    else  fHistPiPMonitorCuts[isSecd]->Fill(28);
 
    if( dcaNegToVertex < fDCADaughtersToVtxSmall || dcaPosToVertex < fDCADaughtersToVtxLarge)  cutOKLambda = kFALSE;
    else  fHistPiPMonitorCuts[isSecd]->Fill(29);

    if(fRapCutV0 && fabs(rapL) > fRap) cutOKLambda = kFALSE;
    else  fHistPiPMonitorCuts[isSecd]->Fill(30);
       
   
    /*	 
	 if(chi2LambdaC > fChiCutKf) cutOKLambda = kFALSE;
	 else  fHistPiPMonitorCuts[isSecd]->Fill(20);
    */

    if(opAng < fOpengAngleDaughters && fabs(ptLambda) < fOpAngPtCut )  cutOKLambda = kFALSE;
    else  fHistPiPMonitorCuts[isSecd]->Fill(31);
    

    if(alfa<fAlfaCut  || (fArmCutL && qt > fQtCut)) cutOKLambda = kFALSE;
    else  fHistPiPMonitorCuts[isSecd]->Fill(32);

      
    if(fEtaSignCut *eta < 0.0) cutOKLambda = kFALSE; 
    //--------------------------- ALambda cuts --------------------------//

    if(dcaV0ToPrimVertex > fDCAToVertexL) cutOKALambda = kFALSE;
    else fHistPiAPMonitorCuts[isSecd]->Fill(24);
 
    //    if(fabs(xr[2])> fDCAZ) cutOKALambda = kFALSE;//continue;//like decay radius z component
    if(fabs(dcaToVertexZNeg) > fDCAZ  || fabs(dcaToVertexZPos) > fDCAZ ) cutOKALambda = kFALSE;
    else fHistPiAPMonitorCuts[isSecd]->Fill(25);

    Double_t ctAL = 0.0,ctTAL=0.0;
    if(fabs(pALambda)>0.0)  ctAL  = decayLength*1.115683/fabs(pALambda);
    if(fabs(ptALambda)>0.0) ctTAL = dim2V0Radius*1.115683/fabs(ptALambda);
    if(ctAL > fCtauL &&  fabs(ptALambda) <fCtauPtCutL)  cutOKALambda = kFALSE;
    else  fHistPiAPMonitorCuts[isSecd]->Fill(26);

    if((cosOPAng<fCosPointAngL && fabs(ptALambda) > fCPAPtCutL)|| cosOPAng<0.99)  cutOKALambda = kFALSE;///xxx
    else fHistPiAPMonitorCuts[isSecd]->Fill(27);
      
    if(dcaDaughters > fDCADaughtersAL )cutOKALambda = kFALSE;
    else  fHistPiAPMonitorCuts[isSecd]->Fill(28);
	 
    if( dcaPosToVertex < fDCADaughtersToVtxSmall || dcaNegToVertex < fDCADaughtersToVtxLarge)  cutOKALambda = kFALSE;
    else fHistPiAPMonitorCuts[isSecd]->Fill(29);
	 
    if(fRapCutV0 && fabs(rapAL) > fRap) cutOKALambda = kFALSE;
    else fHistPiAPMonitorCuts[isSecd]->Fill(30);

    /*
      if(chi2ALambdaC > fChiCutKf) cutOKALambda = kFALSE;
      else  fHistPiAPMonitorCuts[isSecd]->Fill(20);
    */
     
    if(opAng < fOpengAngleDaughters && fabs(ptALambda) < fOpAngPtCut )  cutOKALambda = kFALSE;
    else  fHistPiAPMonitorCuts[isSecd]->Fill(31);
    
      
    if((fArmCutL && qt>qtval) || alfa > -1.0*fAlfaCut) cutOKALambda = kFALSE;
    else  fHistPiAPMonitorCuts[isSecd]->Fill(32);

    if(fEtaSignCut *eta < 0.0) cutOKALambda = kFALSE;
    //---------- check pdg codes of BG --------------------//
    
    Int_t pdgBG = 0;
    if(fMCMode && !fMCTruthMode)  pdgBG = FindPDGCode(stackRec,trackPos,trackNeg);
   
    //----------------------------------------------- V0 ana -----------------------------------------------------------------------//

    //-- cut flags for furhter histos--//
    Bool_t k0sOK=kFALSE;
    Bool_t lambdaOK=kFALSE;
    Bool_t alambdaOK=kFALSE;

      
    //------  Check for K0 ------//
    Bool_t exMass = kFALSE;
    if(fabs(1.115 - massLambda)  < fExcludeLambdaFromK0s){
      cutOKK0s = kFALSE;
      exMass = kTRUE;
    }
    if(fabs(1.115 - massALambda) < fExcludeLambdaFromK0s){
      cutOKK0s = kFALSE;
      exMass = kTRUE;
    }
   
    if(fabs(massPhoton) < fExcludePhotonsFromK0s) {
      cutOKK0s = kFALSE;
      exMass = kTRUE;
    }
      
    if(ptK0s >fMinPt){
      if( cutOKK0s  && fillK0sMC ){
	fHistDedxPiPlus[isSecd]->Fill(tpcMomPos,tpcsigPos);
	fHistDedxPiMinus[isSecd]->Fill(tpcMomNeg,tpcsigNeg);
	fHistPiPiMonitorCuts->Fill(33);
	if(pipidEdx){
	  fHistPiPiMonitorCuts->Fill(34);
	  k0sOK = kTRUE;		    
	  if(!exMass && massK0s > fK0sLowMassCut && massK0s < fK0sHighMassCut ){
	    if(!(fMCMode && fMCTruthMode)){
	      ptV0MC = ptK0s;
	      declengthV0MC = dim2V0Radius;
	    }
	    fHistPiPiMonitorCuts->Fill(35);

	    v0idForDauPositionK0s[iV0MI] = massK0s;
	    
	    fHistPiPiMass->Fill(massK0s);
	    fHistPiPiMassVSPt->Fill(massK0s,ptK0s);
	    fHistPiPiMassVSPtMCTruth->Fill(massK0s,ptV0MC);
	    if(rapK0s >0.0) fHistPiPiMassVSPtPosMCTruth->Fill(massK0s,ptV0MC);//check rap dep
	    else   fHistPiPiMassVSPtNegMCTruth->Fill(massK0s,ptV0MC);
	    fHistPiPiMassVSY->Fill(massK0s,rapK0s);
	    fHistPiPiPtVSY->Fill(rapK0s,ptK0s);
	    fHistPiPiDecayLengthVsMass->Fill(massK0s,dim2V0Radius);//decayLength);
	    //  fHistPiPiDistDaughtersTPCEntrVsMass->Fill(massK0s,distTPCinner);
	    // fHistPiPiPhiPosVsPtPosVsMass->Fill(massK0s,ctTK0,ptV0MC);//,ctK0);//posDaughterPhi);//xxx
	    /*
	      Double_t valTHnK0s[4]= {massK0s,ptV0MC,dim2V0Radius,distTPCinner};
	      fTHnFK0s->Fill(valTHnK0s);	    
	   
	      Double_t valTHnK0sDauEta[4]= {massK0s,ptV0MC,posDaughterEta,negDaughterEta};
	      Double_t valTHnK0sDauPhi[5]= {massK0s,posDaughterPhi,negDaughterPhi,double(nclsITSPos),double(nclsITSNeg)};

	      fTHnFK0sDauEta->Fill(valTHnK0sDauEta);
	      fTHnFK0sDauPhi->Fill(valTHnK0sDauPhi);
	    */
	    /*
	      if(fMCMode && !fMCTruthMode){
	      fHistPiPiPDGCode->Fill(pdgBG);
	      if(pdgBG == -1)   fHistPiPiNoMother->Fill(massK0s,ptV0MC);
	      if(pdgBG == 22)   fHistPiPiGA->Fill(massK0s,ptV0MC);
	      if(pdgBG == 321)  fHistPiPiKch->Fill(massK0s,ptV0MC);
	      if(pdgBG == 333)  fHistPiPiPhi->Fill(massK0s,ptV0MC);
	      if(pdgBG == 3122) fHistPiPiL->Fill(massK0s,ptV0MC);
	      if(pdgBG == 111)  fHistPiPiPi0->Fill(massK0s,ptV0MC);
	      if(pdgBG == 211)  fHistPiPiPich->Fill(massK0s,ptV0MC);
	      if(pdgBG == 113)  fHistPiPiRoh->Fill(massK0s,ptV0MC);
	      if(pdgBG == 223)  fHistPiPiOmega->Fill(massK0s,ptV0MC);
	      if(pdgBG == 313)  fHistPiPiKStar->Fill(massK0s,ptV0MC);
	      if(pdgBG == 310)  fHistPiPiK0s->Fill(massK0s,ptV0MC);
	      if(pdgBG == 130)  fHistPiPiK0L->Fill(massK0s,ptV0MC);
	      if(pdgBG == 2112) fHistPiPiN->Fill(massK0s,ptV0MC);
	      if(pdgBG == 3112 || pdgBG ==3222)  fHistPiPiSigma->Fill(massK0s,ptV0MC);
	      if(pdgBG == 3312 || pdgBG ==3322)  fHistPiPiXi->Fill(massK0s,ptV0MC);
	      if(pdgBG == 2114 || pdgBG ==2224)  fHistPiPiDelta->Fill(massK0s,ptV0MC);
	      if(pdgBG >510 && pdgBG <532)   fHistPiPiB->Fill(massK0s,ptV0MC);
	      if(pdgBG >410  && pdgBG <444)  fHistPiPiD->Fill(massK0s,ptV0MC);
	      if(pdgBG == 331 && pdgBG ==221) fHistPiPiEta->Fill(massK0s,ptV0MC);

	      }
	    */
	    if(massK0s > 0.46 && massK0s < 0.53)  fHistPiPiDecayLengthVsPt->Fill(ptV0MC,dim2V0Radius);//decayLength
	    // fHistPiPiPtDaughters->Fill(posDaughterPt,negDaughterPt);
	    if(!fSetPtDepHist){
	      fHistPiPiRadiusXY->Fill(massK0s,opAng);
	      fHistPiPiCosPointAng->Fill(massK0s,cosOPAng);
	      fHistPiPiDecayLengthVsCtau->Fill(massK0s,ctTK0);
	      fHistPiPiTrackLengthPosVsMass->Fill(massK0s,lengthTPCPos);
	      fHistPiPiTrackLengthNegVsMass->Fill(massK0s,lengthTPCNeg);
	      fHistPiPiDCAZPos->Fill(massK0s,dcaToVertexZPos);
	      fHistPiPiDCAZNeg->Fill(massK0s,dcaToVertexZNeg);
	      fHistPiPiDCADaughters->Fill(massK0s,dcaDaughters);
	      fHistPiPiDCADaughterPosToPrimVtxVSMass->Fill(massK0s,dcaPosToVertex);
	      fHistPiPiDCAVSMass->Fill(massK0s,dcaV0ToPrimVertex);
	      fHistPiPiDCAZVSMass->Fill(massK0s,dcaZToVertex);
	
	    }
	    else{
	      fHistPiPiRadiusXY->Fill(ptV0MC,opAng);
	      fHistPiPiCosPointAng->Fill(ptV0MC,cosOPAng);
	      fHistPiPiDecayLengthVsCtau->Fill(ptV0MC,ctTK0);
	      fHistPiPiTrackLengthPosVsMass->Fill(ptV0MC,lengthTPCPos);
	      fHistPiPiTrackLengthNegVsMass->Fill(ptV0MC,lengthTPCNeg);
	      fHistPiPiDCAZPos->Fill(ptV0MC,dcaToVertexZPos);
	      fHistPiPiDCAZNeg->Fill(ptV0MC,dcaToVertexZNeg);
	      fHistPiPiDCADaughters->Fill(ptV0MC,dcaDaughters);
	      fHistPiPiDCADaughterPosToPrimVtxVSMass->Fill(ptV0MC,dcaPosToVertex);
	      fHistPiPiDCAVSMass->Fill(ptV0MC,dcaV0ToPrimVertex);
	      fHistPiPiDCAZVSMass->Fill(ptV0MC,dcaZToVertex);
	    }

	    if(fMCMode && fMCTruthMode)  fHistPiPiDecayLengthResolution->Fill(declengthV0MC,dim2V0Radius);

	    fHistPiPiK0sVsLambdaMass->Fill(massLambda,massK0s);
	    fHistPiPiK0sVsALambdaMass->Fill(massALambda,massK0s);
	  
	    fHistArmenteros[isSecd]->Fill(alfa,qt);
	    fHistDedxSecPiPlus[isSecd]->Fill(tpcMomPos,tpcsigPos);
	    fHistDedxSecPiMinus[isSecd]->Fill(tpcMomNeg,tpcsigNeg);

	    fHistV0RadiusZ[isSecd]->Fill(dim2V0Radius,xr[2]);
	    fHistV0RadiusXY[isSecd]->Fill(xr[0],xr[1]);
	    fHistV0RadiusZVSPt[isSecd]->Fill(ptK0s,dim2V0Radius);
	    fHistV0RadiusXYVSY[isSecd]->Fill(rapK0s,dim2V0Radius);

	    //-- detector values --/
	    fHistNclsITS[1]->Fill(nclsITSPos,nclsITSNeg);
	    fHistNclsTPC[1]->Fill(crossedRowsTPCNeg,nclsTPCNeg);

	    if(!fSetPtDepHist){
	      fHistNclsITSPosK0->Fill(massK0s,nclsITSPos);
	      fHistNclsITSNegK0->Fill(massK0s,nclsITSNeg);
	      fHistNclsTPCPosK0->Fill(massK0s,nclsTPCPos);
	      fHistNclsTPCNegK0->Fill(massK0s,nclsTPCNeg);
	      fHistChi2PerNclsITSPosK0->Fill(massK0s,chi2PerClusterITSPos);
	      fHistChi2PerNclsITSNegK0->Fill(massK0s,chi2PerClusterITSNeg);
	      fHistNCRowsTPCPosK0->Fill(massK0s,crossedRowsTPCPos);
	      fHistNCRowsTPCNegK0->Fill(massK0s,crossedRowsTPCNeg);
	      fHistRatioFoundOverFinableTPCK0Neg->Fill(massK0s,ratioFoFi);
	      fHistRatioFoundOverFinableTPCK0Pos->Fill(massK0s,ratioFoFiPos);
	    }
	    else{
	      fHistNclsITSPosK0->Fill(ptV0MC,nclsITSPos);
	      fHistNclsITSNegK0->Fill(ptV0MC,nclsITSNeg);
	      fHistNclsTPCPosK0->Fill(ptV0MC,nclsTPCPos);
	      fHistNclsTPCNegK0->Fill(ptV0MC,nclsTPCNeg);
	      fHistChi2PerNclsITSPosK0->Fill(ptV0MC,chi2PerClusterITSPos);
	      fHistChi2PerNclsITSNegK0->Fill(ptV0MC,chi2PerClusterITSNeg);
	      fHistNCRowsTPCPosK0->Fill(ptV0MC,crossedRowsTPCPos);
	      fHistNCRowsTPCNegK0->Fill(ptV0MC,crossedRowsTPCNeg);
	      fHistRatioFoundOverFinableTPCK0Neg->Fill(ptV0MC,ratioFoFi);
	      fHistRatioFoundOverFinableTPCK0Pos->Fill(ptV0MC,ratioFoFiPos);
	    }
	  }
	}
      }
    }


    //------  Check for Lambda -------//
    Bool_t  exMassL =kFALSE;
    if(fabs(0.497 - massK0s) < fExcludeK0sFromLambda){
      cutOKLambda = kFALSE;
      exMassL = kTRUE;
    }
    if(fabs(massPhoton) < fExcludePhotonsFromLambda) {
      cutOKLambda = kFALSE;
      exMassL = kTRUE;
    }
    
    if(ptLambda > fMinPt){
      if(cutOKLambda && fillLambdaMC){
	fHistDedxProt[isSecd]->Fill(tpcMomPos,tpcsigPos);
	fHistDedxPiMinus[isSecd]->Fill(tpcMomNeg,tpcsigNeg);
	fHistPiPMonitorCuts[isSecd]->Fill(33);
	if(pipdEdx){
	  fHistPiPMonitorCuts[isSecd]->Fill(34);
	  lambdaOK = kTRUE;
	  if(!exMassL && massLambda > fLLowMassCut && massLambda < fLHighMassCut){// 1.05 && massLambda < 1.25 ){
	    if(!(fMCMode  && fMCTruthMode)) {
	      ptV0MC = ptLambda;
	      declengthV0MC = dim2V0Radius;
	    }
	    fHistPiPMonitorCuts[isSecd]->Fill(35);

	    v0idForDauPositionL[iV0MI] = massLambda;

	    fHistPiPMass[isSecd]->Fill(massLambda);
	    fHistPiPMassVSPt[isSecd]->Fill(massLambda,ptLambda);
	    fHistPiPMassVSPtMCTruth[isSecd]->Fill(massLambda,ptV0MC);
	    if(rapL >0.0) fHistPiPMassVSPtPosMCTruth[isSecd]->Fill(massLambda,ptV0MC);//check rap dep
	    else  fHistPiPMassVSPtNegMCTruth[isSecd]->Fill(massLambda,ptV0MC);
	    fHistPiPMassVSY[isSecd]->Fill(massLambda,rapL);
	    fHistPiPPtVSY[isSecd]->Fill(rapL,ptLambda);
	    //  fHistPiPDistDaughtersTPCEntrVsMass->Fill(massLambda,distTPCinner);
	    //fHistPiPDecayLengthVsPt[isSecd]->Fill(ptLambda,ctL);
	    /*
	      Double_t valTHnL[4]= {massLambda,ptV0MC,dim2V0Radius,distTPCinner};
	      fTHnFL->Fill(valTHnL);
	      Double_t valTHnL[4]= {massLambda,ptV0MC,eta,phi};
	      Double_t valTHnLDauEta[4]= {massLambda,ptV0MC,posDaughterEta,negDaughterEta};
	      Double_t valTHnLDauPhi[5]= {massLambda,posDaughterPhi,negDaughterPhi,double(nclsITSPos),double(nclsITSNeg)};


	      fTHnFLDauEta->Fill(valTHnLDauEta);
	      fTHnFLDauPhi->Fill(valTHnLDauPhi);
	    */
	    /*	      
		      if(fMCMode && !fMCTruthMode) {
		      fHistPiPPDGCode->Fill(pdgBG);
		      if(pdgBG == 22)  fHistPiPGA->Fill(massLambda,ptV0MC);
		      if(pdgBG == 321) fHistPiPKch->Fill(massLambda,ptV0MC);
		      if(pdgBG == 310) fHistPiPK0s->Fill(massLambda,ptV0MC);
		      if(pdgBG == 111) fHistPiPPi0->Fill(massLambda,ptV0MC);
		      if(pdgBG == 211) fHistPiPPich->Fill(massLambda,ptV0MC);
		      if(pdgBG == 313) fHistPiPKStar->Fill(massLambda,ptV0MC);
		      if(pdgBG == 2112) fHistPiPN->Fill(massLambda,ptV0MC);
		      if(pdgBG == 3122) fHistPiPL->Fill(massLambda,ptV0MC);
		      if(pdgBG == -1)  fHistPiPNoMother->Fill(massLambda,ptV0MC);
		      }
	    */
	    if( massLambda > 1.108 && massLambda < 1.123 ) fHistPiPDecayLengthVsPt[isSecd]->Fill(ptV0MC,dim2V0Radius);//decayLength);
	    fHistPiPDecayLengthVsMass[isSecd]->Fill(massLambda,dim2V0Radius);//decayLength);
	      
	    if(!fSetPtDepHist){
	      fHistPiPRadiusXY[isSecd]->Fill(massLambda,opAng);
	      fHistPiPCosPointAng[isSecd]->Fill(massLambda,cosOPAng);
	      fHistPiPTrackLengthPosVsMass[isSecd]->Fill(massLambda,lengthTPCPos);
	      fHistPiPTrackLengthNegVsMass[isSecd]->Fill(massLambda,lengthTPCNeg);
	      fHistPiPDCAZPos[isSecd]->Fill(massLambda,dcaToVertexZPos);
	      fHistPiPDCAZNeg[isSecd]->Fill(massLambda,dcaToVertexZNeg);
	      fHistPiPDCADaughters[isSecd]->Fill(massLambda,dcaDaughters);
	      fHistPiPDCAVSMass[isSecd]->Fill(massLambda,dcaV0ToPrimVertex);
	      fHistPiPDCAZVSMass[isSecd]->Fill(massLambda,dcaZToVertex);
	      fHistPiPDCADaughterPosToPrimVtxVSMass[isSecd]->Fill(massLambda,dcaPosToVertex);
	      fHistPiPDCADaughterNegToPrimVtxVSMass[isSecd]->Fill(massLambda,dcaNegToVertex);
	      fHistPiPDecayLengthVsCtau[isSecd]->Fill(massLambda,ctTL);
	    }
	    else{
	      fHistPiPRadiusXY[isSecd]->Fill(ptV0MC,opAng);
	      fHistPiPCosPointAng[isSecd]->Fill(ptV0MC,cosOPAng);
	      fHistPiPTrackLengthPosVsMass[isSecd]->Fill(ptV0MC,lengthTPCPos);
	      fHistPiPTrackLengthNegVsMass[isSecd]->Fill(ptV0MC,lengthTPCNeg);
	      fHistPiPDCAZPos[isSecd]->Fill(ptV0MC,dcaToVertexZPos);
	      fHistPiPDCAZNeg[isSecd]->Fill(ptV0MC,dcaToVertexZNeg);
	      fHistPiPDCADaughters[isSecd]->Fill(ptV0MC,dcaDaughters);
	      fHistPiPDCAVSMass[isSecd]->Fill(ptV0MC,dcaV0ToPrimVertex);
	      fHistPiPDCAZVSMass[isSecd]->Fill(ptV0MC,dcaZToVertex);
	      fHistPiPDCADaughterPosToPrimVtxVSMass[isSecd]->Fill(ptV0MC,dcaPosToVertex);
	      fHistPiPDCADaughterNegToPrimVtxVSMass[isSecd]->Fill(ptV0MC,dcaNegToVertex);
	      fHistPiPDecayLengthVsCtau[isSecd]->Fill(ptV0MC,ctTL);
	    }

	    if(fMCMode && fMCTruthMode)  fHistPiPDecayLengthResolution[isSecd]->Fill(declengthV0MC,dim2V0Radius);
	    //   fHistPiPPtDaughters[isSecd]->Fill(posDaughterPt,negDaughterPt);
	    fHistPiPK0sVsLambdaMass->Fill(massLambda,massK0s);
	    fHistPiPALambdaVsLambdaMass->Fill(massLambda,massALambda);
	    	    
	    //-- secondaries --//
	    if(isSecd==1){
	      if(fabs(pdgMother) == 3112 || fabs(pdgMother) == 3114 || fabs(pdgMother) == 3222 || fabs(pdgMother) == 3224 || fabs(pdgMother) == 3214 ){
		fHistPiPMassVSPtSecSigma[1]->Fill(massLambda,ptLambda);
	      }
	      if(pdgMother == 3322 || pdgMother == 3312){//Xi0 and xi minus
		fHistPiPCosPointAngXiVsPt->Fill(ptLambda,cosOPAng);
		fHistPiPMassVSPtSecXi[1]->Fill(massLambda,ptLambda);
		fHistPiPMassVSPtSecXiMCTruth->Fill(massLambda,ptV0MC);
		fHistPiPMassVSYSecXi[1]->Fill(massLambda,rapL);
		fHistPiPXi0PtVSLambdaPt[1]->Fill(ptLambda,ptXiMother);
	      }
	      if(pdgMother == 3334){//Omega
		fHistPiPMassVSPtSecOmega[1]->Fill(massLambda,ptLambda);
		fHistPiPMassVSPtSecOmegaMCTruth->Fill(massLambda,ptV0MC);
		fHistPiPOmegaPtVSLambdaPt[1]->Fill(ptLambda,ptXiMother);
	      }  
	    }

	    if(ptLambda > 0.4) fHistArmenteros[isSecd]->Fill(alfa,qt);
	    fHistV0RadiusZ[isSecd]->Fill(dim2V0Radius,xr[2]);
	    fHistV0RadiusXY[isSecd]->Fill(xr[0],xr[1]);
	    fHistV0RadiusZVSPt[isSecd]->Fill(ptLambda,dim2V0Radius);
	    fHistV0RadiusXYVSY[isSecd]->Fill(rapL,dim2V0Radius);
	    fHistDedxSecProt[isSecd]->Fill(tpcMomPos,tpcsigPos);
	    fHistDedxSecPiMinus[isSecd]->Fill(tpcMomNeg,tpcsigNeg);

	    if(!fSetFillDetAL){
	      //-- detector values --//
	      fHistNclsITS[0]->Fill(nclsITSPos,nclsITSNeg);
	      fHistNclsTPC[0]->Fill(crossedRowsTPCNeg,nclsTPCNeg);
	      if(!fSetPtDepHist){
		fHistNclsITSPosL[isSecd]->Fill(massLambda,nclsITSPos);
		fHistNclsITSNegL[isSecd]->Fill(massLambda,nclsITSNeg);
		fHistNclsTPCPosL[isSecd]->Fill(massLambda,nclsTPCPos);
		fHistNclsTPCNegL[isSecd]->Fill(massLambda,nclsTPCNeg);
		fHistChi2PerNclsITSPosL[isSecd]->Fill(massLambda,chi2PerClusterITSPos);
		fHistChi2PerNclsITSNegL[isSecd]->Fill(massLambda,chi2PerClusterITSNeg);
		fHistNCRowsTPCPosL[isSecd]->Fill(massLambda,crossedRowsTPCPos);
		fHistNCRowsTPCNegL[isSecd]->Fill(massLambda,crossedRowsTPCNeg);
		fHistRatioFoundOverFinableTPCLNeg[isSecd]->Fill(massLambda,ratioFoFi);
		fHistRatioFoundOverFinableTPCLPos[isSecd]->Fill(massLambda,ratioFoFiPos);
	      }
	      else{
		fHistNclsITSPosL[isSecd]->Fill(ptV0MC,nclsITSPos);
		fHistNclsITSNegL[isSecd]->Fill(ptV0MC,nclsITSNeg);
		fHistNclsTPCPosL[isSecd]->Fill(ptV0MC,nclsTPCPos);
		fHistNclsTPCNegL[isSecd]->Fill(ptV0MC,nclsTPCNeg);
		fHistChi2PerNclsITSPosL[isSecd]->Fill(ptV0MC,chi2PerClusterITSPos);
		fHistChi2PerNclsITSNegL[isSecd]->Fill(ptV0MC,chi2PerClusterITSNeg);
		fHistNCRowsTPCPosL[isSecd]->Fill(ptV0MC,crossedRowsTPCPos);
		fHistNCRowsTPCNegL[isSecd]->Fill(ptV0MC,crossedRowsTPCNeg);
		fHistRatioFoundOverFinableTPCLNeg[isSecd]->Fill(ptV0MC,ratioFoFi);
		fHistRatioFoundOverFinableTPCLPos[isSecd]->Fill(ptV0MC,ratioFoFiPos);
	      }
	    }
	  }
	  
	}
      }
    }


    //-- Check for AntiLambda --//    
    Bool_t  exMassAL =kFALSE;
    if(fabs(0.497 - massK0s) < fExcludeK0sFromLambda){
      exMassAL = kTRUE;
    }
    if(fabs(massPhoton) < fExcludePhotonsFromLambda) {
      exMassAL = kTRUE;
    }

    if(ptALambda > fMinPt){
      if(cutOKALambda && fillALambdaMC){
	fHistDedxAProt[isSecd]->Fill(tpcMomNeg,tpcsigNeg);
	fHistDedxPiPlus[isSecd]->Fill(tpcMomPos,tpcsigPos);
	fHistPiAPMonitorCuts[isSecd]->Fill(33);
	if(piapdEdx){
	  fHistPiAPMonitorCuts[isSecd]->Fill(34);
	  alambdaOK = kTRUE;
	  if( !exMassAL && massALambda > fLLowMassCut && massALambda < fLHighMassCut){//1.05 && massALambda < 1.25  ){
	    if(!(fMCMode && fMCTruthMode)) {
	      ptV0MC = ptALambda;
	      declengthV0MC = dim2V0Radius;
	    }
	    fHistPiAPMonitorCuts[isSecd]->Fill(35);

	    v0idForDauPositionAL[iV0MI] = massALambda;

	    fHistPiAPMass[isSecd]->Fill(massALambda);
	    fHistPiAPMassVSPt[isSecd]->Fill(massALambda,ptALambda);
	    fHistPiAPMassVSPtMCTruth[isSecd]->Fill(massALambda,ptV0MC);
	    if(rapAL >0.0) fHistPiAPMassVSPtPosMCTruth[isSecd]->Fill(massALambda,ptV0MC);//check rap dep
	    else  fHistPiAPMassVSPtNegMCTruth[isSecd]->Fill(massALambda,ptV0MC);
	    fHistPiAPMassVSY[isSecd]->Fill(massALambda,rapAL);
	    fHistPiAPPtVSY[isSecd]->Fill(rapAL,ptALambda);
	    //   fHistPiAPDistDaughtersTPCEntrVsMass->Fill(massALambda,distTPCinner);
	    //  fHistPiAPPtDaughters[isSecd]->Fill(posDaughterPt,negDaughterPt);
	    //  fHistPiAPPhiPosVsPtPosVsMass->Fill(massALambda,ctTAL,ptV0MC);//
	    //  if(isSecd < 1) fHistPiPPhiPosVsPtPosVsMass->Fill(massALambda,ctTL,ptV0MC);//
	    // else {
	    //	if(fMCTruthMode) fHistPiAPPhiPosVsPtPosVsMass->Fill(massALambda,ctTL,ptV0MC);
	    // }	      
	    /*
	      Double_t valTHnAL[4]= {massALambda,ptV0MC,dim2V0Radius,distTPCinner};
	      fTHnFAL->Fill(valTHnAL);
	      Double_t valTHnAL[4]= {massALambda,ptV0MC,eta,phi};
	      Double_t valTHnALDauEta[4]={massALambda,ptV0MC,posDaughterEta,negDaughterEta};
	      Double_t valTHnALDauPhi[5]={massALambda,posDaughterPhi,negDaughterPhi,double(nclsITSPos),double(nclsITSNeg)};
	       
	      fTHnFALDauEta->Fill(valTHnALDauEta);
	      fTHnFALDauPhi->Fill(valTHnALDauPhi);
	      fTHnFAL->Fill(valTHnAL);
	    */
	    if(fMCMode && !fMCTruthMode) fHistPiAPPDGCode->Fill(pdgBG);
	    if( massALambda>1.108 && massALambda<1.123 )  fHistPiAPDecayLengthVsPt[isSecd]->Fill(ptV0MC,dim2V0Radius);//decayLength);
	    fHistPiAPDecayLengthVsMass[isSecd]->Fill(massALambda,dim2V0Radius);//decayLength);
	    
	    if(!fSetPtDepHist){
	      fHistPiAPRadiusXY[isSecd]->Fill(massALambda,opAng);
	      fHistPiAPCosPointAng[isSecd]->Fill(massALambda,cosOPAng);	    
	      fHistPiAPTrackLengthPosVsMass[isSecd]->Fill(massALambda,lengthTPCPos);
	      fHistPiAPTrackLengthNegVsMass[isSecd]->Fill(massALambda,lengthTPCNeg);
	      //fHistPiAPDCAZPos[isSecd]->Fill(massALambda,dcaToVertexZPos);
	      //fHistPiAPDCAZNeg[isSecd]->Fill(massALambda,dcaToVertexZNeg);
	      fHistPiAPDCADaughters[isSecd]->Fill(massALambda,dcaDaughters);
	      fHistPiAPDCAVSMass[isSecd]->Fill(massALambda,dcaV0ToPrimVertex);
	      fHistPiAPDCAZVSMass[isSecd]->Fill(massALambda,dcaZToVertex);
	      fHistPiAPDCADaughterPosToPrimVtxVSMass[isSecd]->Fill(massALambda,dcaPosToVertex);
	      fHistPiAPDCADaughterNegToPrimVtxVSMass[isSecd]->Fill(massALambda,dcaNegToVertex);
	      fHistPiAPDecayLengthVsCtau[isSecd]->Fill(massALambda,ctTAL);
	    }
	    else{
	      fHistPiAPRadiusXY[isSecd]->Fill(ptV0MC,opAng);
	      fHistPiAPCosPointAng[isSecd]->Fill(ptV0MC,cosOPAng);
	      fHistPiAPTrackLengthPosVsMass[isSecd]->Fill(ptV0MC,lengthTPCPos);
	      fHistPiAPTrackLengthNegVsMass[isSecd]->Fill(ptV0MC,lengthTPCNeg);
	      //fHistPiAPDCAZPos[isSecd]->Fill(ptV0MC,dcaToVertexZPos);
	      //fHistPiAPDCAZNeg[isSecd]->Fill(ptV0MC,dcaToVertexZNeg);
	      fHistPiAPDCADaughters[isSecd]->Fill(ptV0MC,dcaDaughters);
	      fHistPiAPDCAVSMass[isSecd]->Fill(ptV0MC,dcaV0ToPrimVertex);
	      fHistPiAPDCAZVSMass[isSecd]->Fill(ptV0MC,dcaZToVertex);
	      fHistPiAPDCADaughterPosToPrimVtxVSMass[isSecd]->Fill(ptV0MC,dcaPosToVertex);
	      fHistPiAPDCADaughterNegToPrimVtxVSMass[isSecd]->Fill(ptV0MC,dcaNegToVertex);
	      fHistPiAPDecayLengthVsCtau[isSecd]->Fill(ptV0MC,ctTAL);
	    }

	    if(fMCMode && fMCTruthMode)  fHistPiAPDecayLengthResolution[isSecd]->Fill(declengthV0MC,dim2V0Radius);
	    fHistPiAPK0sVsALambdaMass->Fill(massALambda,massK0s);
	    fHistPiAPLambdaVsALambdaMass->Fill(massALambda,massLambda);

	    //-- secondaries --//
	    if(isSecd == 1){
	      if(fabs(pdgMother) == 3112 || fabs(pdgMother) == 3114 || fabs(pdgMother) == 3222 || fabs(pdgMother) == 3224 || fabs(pdgMother) ==  3214 ){
		fHistPiAPMassVSPtSecSigma[1]->Fill(massLambda,ptLambda);
	      }
	      if(pdgMother == -3322 || pdgMother == -3312){//Xi0 and xiplus
		fHistPiAPCosPointAngXiVsPt->Fill(ptALambda,cosOPAng);
		fHistPiAPMassVSPtSecXi[1]->Fill(massALambda,ptALambda);
		fHistPiAPMassVSPtSecXiMCTruth->Fill(massALambda,ptV0MC);
		fHistPiAPMassVSYSecXi[1]->Fill(massALambda,rapAL);
		fHistPiAPXi0PtVSLambdaPt[1]->Fill(ptALambda,ptXiMother);
	      }
	      if(pdgMother == -3334){//Omega
		fHistPiAPMassVSPtSecOmega[1]->Fill(massALambda,ptALambda);
		fHistPiAPMassVSPtSecOmegaMCTruth->Fill(massALambda,ptV0MC);
		fHistPiAPOmegaPtVSLambdaPt[1]->Fill(ptALambda,ptXiMother);
	      }  
	    }
  
	    if(ptALambda > 0.4) fHistArmenteros[isSecd]->Fill(alfa,qt);
	    fHistDedxSecAProt[isSecd]->Fill(tpcMomNeg,tpcsigNeg);
	    fHistDedxSecPiPlus[isSecd]->Fill(tpcMomPos,tpcsigPos);
	    fHistV0RadiusZ[isSecd]->Fill(dim2V0Radius,xr[2]);
	    fHistV0RadiusXY[isSecd]->Fill(xr[0],xr[1]);
	    fHistV0RadiusZVSPt[isSecd]->Fill(ptALambda,dim2V0Radius);
	    fHistV0RadiusXYVSY[isSecd]->Fill(rapAL,dim2V0Radius);

	    if(fSetFillDetAL){
	      //-- detector values --//
	      fHistNclsITS[0]->Fill(nclsITSPos,nclsITSNeg);
	      fHistNclsTPC[0]->Fill(crossedRowsTPCNeg,nclsTPCNeg);
	      if(!fSetPtDepHist){
		fHistNclsITSPosL[isSecd]->Fill(massALambda,nclsITSPos);
		fHistNclsITSNegL[isSecd]->Fill(massALambda,nclsITSNeg);
		fHistNclsTPCPosL[isSecd]->Fill(massALambda,nclsTPCPos);
		fHistNclsTPCNegL[isSecd]->Fill(massALambda,nclsTPCNeg);
		fHistChi2PerNclsITSPosL[isSecd]->Fill(massALambda,chi2PerClusterITSPos);
		fHistChi2PerNclsITSNegL[isSecd]->Fill(massALambda,chi2PerClusterITSNeg);
		fHistNCRowsTPCPosL[isSecd]->Fill(massALambda,crossedRowsTPCPos);
		fHistNCRowsTPCNegL[isSecd]->Fill(massALambda,crossedRowsTPCNeg);
		fHistRatioFoundOverFinableTPCLNeg[isSecd]->Fill(massALambda,ratioFoFi);
		fHistRatioFoundOverFinableTPCLPos[isSecd]->Fill(massALambda,ratioFoFiPos);
	      }
	      else{
		fHistNclsITSPosL[isSecd]->Fill(ptV0MC,nclsITSPos);
		fHistNclsITSNegL[isSecd]->Fill(ptV0MC,nclsITSNeg);
		fHistNclsTPCPosL[isSecd]->Fill(ptV0MC,nclsTPCPos);
		fHistNclsTPCNegL[isSecd]->Fill(ptV0MC,nclsTPCNeg);
		fHistChi2PerNclsITSPosL[isSecd]->Fill(ptV0MC,chi2PerClusterITSPos);
		fHistChi2PerNclsITSNegL[isSecd]->Fill(ptV0MC,chi2PerClusterITSNeg);
		fHistNCRowsTPCPosL[isSecd]->Fill(ptV0MC,crossedRowsTPCPos);
		fHistNCRowsTPCNegL[isSecd]->Fill(ptV0MC,crossedRowsTPCNeg);
		fHistRatioFoundOverFinableTPCLNeg[isSecd]->Fill(ptV0MC,ratioFoFi);
		fHistRatioFoundOverFinableTPCLPos[isSecd]->Fill(ptV0MC,ratioFoFiPos);
	      }
	    }
	  }
	}
      }
    }

    if(lambdaOK || alambdaOK || k0sOK) {
      trackID[iV0MI][0] = v0MIs->GetPindex();
      trackID[iV0MI][1] = v0MIs->GetNindex();
    }

      
    //-- fill detector histos general --//
    if((lambdaOK && !exMassL) || (alambdaOK && !exMassAL) || (k0sOK && !exMass)){
      fHistPiPiEtaDReco[1]->Fill(posDaughterEta,eta00);
      fHistPiPEtaDReco[1]->Fill(negDaughterEta,eta01);
    }
      

    /*
    //-- AliKFParticle --//
    if (negPiKF) delete negPiKF; negPiKF=NULL;
    if (posPiKF) delete posPiKF; posPiKF=NULL;
    if (posPKF)  delete posPKF; posPKF=NULL;
    if (negAPKF) delete negAPKF; negAPKF=NULL;
    */

  }//---- end V0 reco loop----//
   /*
     CheckDistanceOfDaughters( nV0,v0idForDauPositionK0s,magField,0);
     CheckDistanceOfDaughters( nV0,v0idForDauPositionL,magField,1);
     CheckDistanceOfDaughters( nV0,v0idForDauPositionAL,magField,2);
   */
}
  
//________________________________________________________________________

Int_t AliAnalysisTaskV0ForRAA::CalculateCentralityBin(){
  //find centrality bin for centrality selection

  if (fUseCentrality == 0) return -1;

  AliCentrality *esdCentrality = fESD->GetCentrality();

  Float_t centralityVZERO  = esdCentrality->GetCentralityPercentile("V0M");  
  Float_t centralitySPD    = esdCentrality->GetCentralityPercentile("CL1");

  Int_t centralityVZEROBin = -1;
  Int_t centralitySPDBin   = -1;

  //-- SPD centrality --//
  if ( fUseCentrality == 2 ){
    if      ( centralitySPD >=  0. && centralitySPD <   5.) centralitySPDBin =  0;
    else if ( centralitySPD >=  5. && centralitySPD <  10.) centralitySPDBin =  5;
    else if ( centralitySPD >= 10. && centralitySPD <  20.) centralitySPDBin = 10;
    else if ( centralitySPD >= 20. && centralitySPD <  30.) centralitySPDBin = 20;
    else if ( centralitySPD >= 30. && centralitySPD <  40.) centralitySPDBin = 30;
    else if ( centralitySPD >= 40. && centralitySPD <  50.) centralitySPDBin = 40;
    else if ( centralitySPD >= 50. && centralitySPD <  60.) centralitySPDBin = 50;
    else if ( centralitySPD >= 60. && centralitySPD <  70.) centralitySPDBin = 60;
    else if ( centralitySPD >= 70. && centralitySPD <  80.) centralitySPDBin = 70;
    else if ( centralitySPD >= 80. && centralitySPD <  90.) centralitySPDBin = 80;
    else if ( centralitySPD >= 90. && centralitySPD <  99.) centralitySPDBin = 90;
    else if ( centralitySPD >= 99. ) centralitySPDBin = 100;
    else if ( fabs(centralitySPD)< 0.0001 ) centralitySPDBin = 100;
    return centralitySPDBin;
  }

  //-- V0 centrality --//
  if ( fUseCentrality == 1 ){
    if      ( centralityVZERO >  0. && centralityVZERO <   5.) centralityVZEROBin =  0;
    else if ( centralityVZERO >=  5. && centralityVZERO <  10.) centralityVZEROBin =  5;
    else if ( centralityVZERO >= 10. && centralityVZERO <  20.) centralityVZEROBin = 10;
    else if ( centralityVZERO >= 20. && centralityVZERO <  30.) centralityVZEROBin = 20;
    else if ( centralityVZERO >= 30. && centralityVZERO <  40.) centralityVZEROBin = 30;
    else if ( centralityVZERO >= 40. && centralityVZERO <  50.) centralityVZEROBin = 40;
    else if ( centralityVZERO >= 50. && centralityVZERO <  60.) centralityVZEROBin = 50;
    else if ( centralityVZERO >= 60. && centralityVZERO <  70.) centralityVZEROBin = 60;
    else if ( centralityVZERO >= 70. && centralityVZERO <  80.) centralityVZEROBin = 70;
    else if ( centralityVZERO >= 80. && centralityVZERO <  90.) centralityVZEROBin = 80;
    else if ( centralityVZERO >= 90. && centralityVZERO <  99.) centralityVZEROBin = 90;
    else if ( centralityVZERO >= 99. ) centralityVZEROBin = 100;
    else if ( fabs(centralityVZERO)< 0.0001 ) centralityVZEROBin = 100;
    return centralityVZEROBin;
  }
  return -1;
  
}

//__________________________________________________________________________________________________________
Bool_t  AliAnalysisTaskV0ForRAA::GetMCTruthPartner(AliESDtrack *pos,AliESDtrack *neg,Int_t id0,Int_t id1){
  //-- get daughter label and check it --//
  Int_t labelP = fabs(pos->GetLabel());
  Int_t labelN = fabs(neg->GetLabel());
    
  if (labelN==labelP)  return kFALSE;
      
  if(fMCTruthMode){
    if ((labelP!=id0) && (labelP!=id1))  return kFALSE;
    if ((labelN!=id0) && (labelN!=id1))  return kFALSE;
  }
    
  return kTRUE;
}

//__________________________________________________________________________________________________________
Bool_t  AliAnalysisTaskV0ForRAA::CheckMultipleV0Candidates(Int_t part1,Int_t part2,Int_t iV0MI,Int_t trackID[][2]){
 
  Bool_t multFoundV0=kFALSE;
  for(Int_t i = 0; i < iV0MI;i++){
    if(trackID[i][0] == part1 && trackID[1][i] == part2) multFoundV0 = kTRUE;
    if(trackID[i][1] == part2 && trackID[1][i] == part1) multFoundV0 = kTRUE;
  }  
  return multFoundV0;
}
//__________________________________________________________________________________________________________
  /*
void  AliAnalysisTaskV0ForRAA::CheckDistanceOfDaughters(Int_t iV0MI,Float_t V0ID[],Double_t magF,Int_t particle){
  //particle 0=K0s, 1=Lambda, 2 = AntiLambda

  //fESD->GetMagneticField()
  AliESDv0 *v01 = NULL;
  AliExternalTrackParam *parPos1 = NULL;
  AliExternalTrackParam *parNeg1 = NULL;
  AliESDtrack *trackPosTest1 = NULL;
  AliESDtrack *trackNegTest1 = NULL;
  AliESDv0 *v02 = NULL;
  AliExternalTrackParam *parPos2 = NULL;
  AliExternalTrackParam *parNeg2 = NULL;
  AliESDtrack *trackPosTest2 = NULL;
  AliESDtrack *trackNegTest2 = NULL;
 
  for(Int_t i = 0; i < iV0MI-1;i++){
    if(V0ID[i] < 0.00001) continue;
    v01 = fESD->GetV0(i);  
    //-- esd tracks --//
    trackPosTest1 = fESD->GetTrack(v01->GetPindex());
    trackNegTest1 = fESD->GetTrack(v01->GetNindex());
    if(trackPosTest1->GetSign() <0){
      parPos1 = new AliExternalTrackParam( *v01->GetParamN());
      parNeg1 = new AliExternalTrackParam( *v01->GetParamP()); 
    }
    else{
      parPos1 = new AliExternalTrackParam( *v01->GetParamP());
      parNeg1 = new AliExternalTrackParam( *v01->GetParamN()); 
    }
    Bool_t dPos = kFALSE;
    Bool_t dNeg = kFALSE;
    Double_t distPos[3], distNeg[3];
    Double_t  averagePos = 0.0;
    Double_t  averageNeg = 0.0;
    Double_t pt = v01->Pt();
    Int_t fillPtHist = 0;
    if(pt > 2.0 && pt <6.0) fillPtHist =1;
    if( pt >= 6.0)  fillPtHist =2;

    for(Int_t k=i+1;k<iV0MI;k++){
      if(V0ID[k] < 0.00001 ) continue;
      v02 = fESD->GetV0(k);  
      //-- esd tracks --//
      trackPosTest2 = fESD->GetTrack(v02->GetPindex());
      trackNegTest2 = fESD->GetTrack(v02->GetNindex());
      if(trackPosTest2->GetSign() <0){
	parPos2 = new AliExternalTrackParam( *v02->GetParamN());
	parNeg2 = new AliExternalTrackParam( *v02->GetParamP()); 
      }
      else{
	parPos2 = new AliExternalTrackParam( *v02->GetParamP());
	parNeg2 = new AliExternalTrackParam( *v02->GetParamN()); 
      }
      
      for(Int_t j = 0;j<9;j++){
	Double_t atX = 84.5+ j*20.0; //every 20 cm in TPC
	dPos = parPos1->GetDistance(parPos2,atX,distPos,magF);
	dNeg = parNeg1->GetDistance(parNeg2,atX,distNeg,magF);
	Double_t rP = sqrt(pow(distPos[0],2.0)+pow(distPos[1],2.0)+pow(distPos[2],2.0));
	Double_t rN = sqrt(pow(distNeg[0],2.0)+pow(distNeg[1],2.0)+pow(distNeg[2],2.0));
	averagePos += rP;
	averageNeg += rN;
      }
      averagePos /= 9.0;
      averageNeg /= 9.0;
      Double_t atRPos1 = 0.0,atRPos2=0.0;
      Double_t atRNeg1 = 0.0,atRNeg2=0.0;
       
      Double_t dcaPos = parPos1->GetDCA(parPos2,magF, atRPos1,atRPos2);
      Double_t dcaNeg = parNeg1->GetDCA(parNeg2,magF, atRNeg1,atRNeg2);
     
      switch(particle){
      case 0:
	if( averagePos <=20.0) {
	  fHistPiPiDistDaughtersPos[fillPtHist]->Fill(V0ID[k],averagePos);
	  fHistPiPiDCADaughtersPos[fillPtHist]->Fill(V0ID[k],dcaPos);      
	}	
	if( averageNeg <=20.0){
	  fHistPiPiDistDaughtersNeg[fillPtHist]->Fill(V0ID[k],averageNeg);
	  fHistPiPiDCADaughtersNeg[fillPtHist]->Fill(V0ID[k],dcaNeg);
	}
	//	if(dcaPos < fDistDauForCheck && atRPos1 > 80.0 && atRPos1 < 240.0) {
	fHistPiPiRadAtDCA5cmDaughtersPos[fillPtHist]->Fill(V0ID[k],atRPos1);
	//	}
	//	if(dcaNeg < fDistDauForCheck && atRNeg2 > 80.0 && atRNeg2 < 240.0) {
	fHistPiPiRadAtDCA5cmDaughtersNeg[fillPtHist]->Fill(V0ID[k],atRNeg1);
	//	}
	break;
      case 1:
	if( averagePos <=20.0) {
	  fHistPiPDistDaughtersPos[fillPtHist]->Fill(V0ID[k],averagePos);
	  fHistPiPDCADaughtersPos[fillPtHist]->Fill(V0ID[k],dcaPos);	
	}
	if( averageNeg <=20.0){
	  fHistPiPDistDaughtersNeg[fillPtHist]->Fill(V0ID[k],averageNeg);
	  fHistPiPDCADaughtersNeg[fillPtHist]->Fill(V0ID[k],dcaNeg);
	}
	//	if(dcaPos < fDistDauForCheck  && atRPos1 > 80.0 && atRPos1 < 240.0) {
	fHistPiPRadAtDCA5cmDaughtersPos[fillPtHist]->Fill(V0ID[k],atRPos1);
	//	}
	//if(dcaNeg < fDistDauForCheck  && atRNeg1 > 80.0 && atRNeg1 < 240.0) {
	fHistPiPRadAtDCA5cmDaughtersNeg[fillPtHist]->Fill(V0ID[k],atRNeg1);
	//	}
	break;
      case 2:
     	if( averagePos <=20.0) {
	  fHistPiAPDistDaughtersPos[fillPtHist]->Fill(V0ID[k],averagePos);
	  fHistPiAPDCADaughtersPos[fillPtHist]->Fill(V0ID[k],dcaPos);	
	}
	if( averageNeg <=20.0){
	  fHistPiAPDistDaughtersNeg[fillPtHist]->Fill(V0ID[k],averageNeg);
	  fHistPiAPDCADaughtersNeg[fillPtHist]->Fill(V0ID[k],dcaNeg);
	}
	//	if(dcaPos < fDistDauForCheck  && atRPos1 > 80.0 && atRPos1 < 240.0) {
	fHistPiAPRadAtDCA5cmDaughtersPos[fillPtHist]->Fill(V0ID[k],atRPos1);
	//}
	//	if(dcaNeg < fDistDauForCheck  && atRNeg1 > 80.0 && atRNeg1 < 240.0) {
	fHistPiAPRadAtDCA5cmDaughtersNeg[fillPtHist]->Fill(V0ID[k],atRNeg1);
	//	}
	break;
      }
     
      delete parPos2;
      delete parNeg2;
    }
    delete parPos1;
    delete parNeg1;

  }  

}
 */
//__________________________________________________________________________________________________________
Int_t  AliAnalysisTaskV0ForRAA::FindPDGCode(AliStack *stackRec,AliESDtrack *trackPos,AliESDtrack *trackNeg){
  
  //-- get daughter label --//
  Int_t labelP = fabs(trackPos->GetLabel());
  Int_t labelN = fabs(trackNeg->GetLabel());
  TParticle *p0 = stackRec->Particle(labelP);
  TParticle *p1 = stackRec->Particle(labelN);
  Int_t idmother0  = -1;
  if(p0) idmother0 = p0->GetMother(0);
  Int_t idmother1  = -1;
  if(p1) idmother1 = p1->GetMother(0);
  Int_t pdg = -1;
  if(idmother0 != idmother1) return pdg;
  else{
    if(idmother0 > -1) {
      pdg = stackRec->Particle(idmother0)->GetPdgCode();
      if(fabs(pdg) >21) return fabs(pdg);
    }
    return -1;
  }
  
}
