/***************************************************************          *
 * Authors : Simone Schuchmann 
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

#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"//xxx
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
  fHistPiPiMassVSY(0),
  fHistPiPiPtVSY(0),
// fHistPiPiMassVSAlpha(0),
  fHistPiPiRadiusXY(0),
  fHistPiPiCosPointAng(0),
  fHistPiPiDCADaughterPosToPrimVtxVSMass(0),  
  fHistPiPiDecayLengthVsPt(0),
  fHistPiPiDecayLengthVsMass(0),
  fHistPiPiDecayLengthVsCtau(0),
// fHistPiPiMassVSPtK0L(0),
  fHistPiPiDCADaughters(0), 
//    fHistPiPiPtDaughters(0),
  fHistPiPiDCAVSMass(0),
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
//------------MC only histos-----------
  fHistPrimVtxZESDVSNContributorsMC(0),
  fHistPrimVtxZESDTPCVSNContributorsMC(0),
  fHistPrimVtxZESDSPDVSNContributorsMC(0),
  fHistMCVertexZ(0),
  fHistPiPiPDGCode(0),
  fHistPiPPDGCode(0),
  fHistPiAPPDGCode(0),
//cosine of pointing angle of Xi vs pt histos
  fHistPiPCosPointAngXiVsPt(0),
  fHistPiAPCosPointAngXiVsPt(0),
  fHistPiPMassVSPtSecXiMCTruth(0),
  fHistPiPMassVSPtSecOmegaMCTruth(0),
  fHistPiAPMassVSPtSecXiMCTruth(0),
  fHistPiAPMassVSPtSecOmegaMCTruth(0),
 // fHistUserPtShift(0),
//  fHistPiPiPhiPosVsPtPosVsMass(0),//xxx
//  fHistPiPPhiPosVsPtPosVsMass(0),//xxx
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
  fSetPtDepHist(0)
  //  fShift(0),
  // fDeltaInvP(0)
{  // Constructor.

  DefineOutput(1,TList::Class());

  // define defaults for globals
  /*
    fShift = kFALSE;                       // shift in charge/pt yes/no
    fDeltaInvP = 0.00;                     // shift value
  */
   
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
  fUsePID = kFALSE;
  fUsePIDPion = kFALSE;
  fMoreNclsThanRows = kFALSE;
  fMoreNclsThanFindable = kFALSE;
  fMoreNclsThanFindableMax = kFALSE;
  fRatioFoundOverFindable = -1.0;
  fRatioMaxCRowsOverFindable = 1000.0;


  fChi2PerClusterITS = 100000.0;
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

  fEtaCutMCDaughters = kFALSE;
  fEtaCutMCDaughtersVal = 50.0;

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

  fDecayRadXYMin=-100000.0;
  fDecayRadXYMax=1000000.0;
  fDecayLengthMax=100000.0;
  fDecayLengthMin=-1000000.0;
   
  fDecRadCutITSMin = 0.0000;
  fDecRadCutITSMax = 10000.0;

  fCosPointAngL=-1.0;
  fCosPointAngK=-1.0;
  fCPAPtCutK0 = 1000.0;
  fCPAPtCutL =1000.0;
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
    fHistPiPMonitorCuts[j] =NULL;
    fHistPiPMonitorMCCuts[j] =NULL;
    fHistPiPDecayLengthResolution[j] =NULL;
    //    fHistPiPDCAZPos[j] =NULL;
    //fHistPiPDCAZNeg[j] =NULL;
    fHistPiPTrackLengthPosVsMass[j] = NULL;
    fHistPiPTrackLengthNegVsMass[j] = NULL;

    //ALambda
    fHistPiAPMass[j]=NULL;
    fHistPiAPMassVSPt[j]=NULL;
    fHistPiAPMassVSY[j] = NULL;
    fHistPiAPMassVSPtMCTruth[j]=NULL;
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
  //create output objects

  Int_t nbPt=800;
  Int_t nbMass=500;

   
  //-----------------  create output container -----------------//

  fOutputContainer = new TList() ;
  fOutputContainer->SetName(GetName()) ;
  fOutputContainer->SetOwner();
  
  Int_t mchist = 1;// for Data
  if((fMCMode && fMCTruthMode) || fMCTruthMode) mchist = 2;
	

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
  fHistV0RadiusZ[0]  = new TH2F("fHistV0RadiusZ","z of decay radius vs 2D radius",100,0.0,100.0,250,-125.0,125.0);
  fHistV0RadiusZVSPt[0]  = new TH2F("fHistV0RadiusZVSPt","z of decay radius vs pt radius",200,0.0,20.0,125,0.0,125.0);
  fHistV0RadiusXY[0]  = new TH2F("fHistV0RadiusXY","y vs x decay radius",250,-125.0,125.0,250,-125.0,125.0);
  fHistV0RadiusXYVSY[0]  = new TH2F("fHistV0RadiusXYVSY","2D decay radius vs rap",100,-1,1,100,0.0,100.0);
  fHistArmenteros[0] = new TH2F("fHistArmenteros"," pi+pi- armenteros",nbMass,-1.,1.,500,0.,0.5);

  // fHistPiPiPhiPosVsPtPosVsMass = new TH3F("fHistPiPiPhiPosVsPtPosVsMass","ctau K0s vs pt vs mass",250,0.25,0.75,440,0.0,110.0,200,0.0,20.0);//4.0);//xxx
  //  fHistPiPPhiPosVsPtPosVsMass  = new TH3F("fHistPiPPhiPosVsPtPosVsMass","ctau L vs pt vs mass",250,0.25,0.75,440,0.0,110.0,200,0.0,20.0);//4.0);//xxx
  fHistPiPiK0sVsLambdaMass  = new TH2F("fHistPiPiK0sVsLambdaMass","K0s mass vs Lambda mass for all pt for K0s",250,1.05,1.25,250,0.25,0.75);
  fHistPiPiK0sVsALambdaMass = new TH2F("fHistPiPiK0sVsALambdaMass","K0s mass vs ALambda mass for all pt for K0s",250,1.05,1.25,250,0.25,0.75);

  fHistPiPK0sVsLambdaMass   = new TH2F("fHistPiPK0sVsLambdaMass","K0s mass vs Lambda mass for all pt for Lambda",250,1.05,1.25,250,0.25,0.75);

  fHistPiAPK0sVsALambdaMass = new TH2F("fHistPiAPK0sVsALambdaMass","K0s mass vs ALambda mass for all pt for ALambda",250,1.05,1.25,250,0.25,0.75);

  fHistPiPALambdaVsLambdaMass  = new TH2F("fHistPiPALambdaVsLambdaMass","ALambda mass vs Lambda mass for Lambda",250,1.05,1.25,250,1.05,1.25);
  fHistPiAPLambdaVsALambdaMass = new TH2F("fHistPiAPLambdaVsALambdaMass","Lambda mass vs ALambda mass for ALambda",250,1.05,1.25,250,1.05,1.25);

  //-----K0s
  fHistPiPiMass = new TH1F("fHistPiPiMass"," pi+pi- InvMass distribution",2*nbMass,0.,2.);
  fHistPiPiMassVSPt = new TH2F("fHistPiPiMassVSPt","pi+pi- InvMass distribution",nbMass,0.25,0.75,300,0.0,30.0);
  fHistPiPiMassVSPtMCTruth = new TH2F("fHistPiPiMassVSPtMCTruth","pi+pi- InvMass distribution vs pt MCTruth",nbMass,0.25,0.75,300,0.0,30.0);
  fHistPiPiMassVSY = new TH2F("fHistPiPiMassVSY","pi+pi- InvMass distribution vs rapidity",nbMass,0.25,0.75,200,-1.0,1.0);
  fHistPiPiPtVSY = new TH2F("fHistPiPiPtVSY","phi vs mass",100,-1,1,100,0.0,20);
  fHistPiPiDecayLengthVsPt = new TH2F("fHistPiPiDecayLengthVsPt","K0 decay length vs pt",200,0.0,20.0,220,0.0,110.0);
  fHistPiPiDecayLengthVsMass = new TH2F("fHistPiPiDecayLengthVsMass","K0s decay length vs mass",nbMass,0.25,0.75,220,0.0,110.0);  
  if(!fSetPtDepHist){
    fHistPiPiDecayLengthVsCtau = new TH2F("fHistPiPiDecayLengthVsCtau","K0s ctau vs mass",nbMass,0.25,0.75,250,0.0,50.0);
  }
    else{
      fHistPiPiDecayLengthVsCtau = new TH2F("fHistPiPiDecayLengthVsCtau","K0s ctau vs pt",200,0,20.0,250,0.0,50.0);
    } 
  
    fHistPiPiMonitorCuts = new TH1F("fHistPiPiMonitorCuts","K0 cut monitor",35,0.5,35.5);
    fHistPiPiMonitorMCCuts = new TH1F("fHistPiPiMonitorMCCuts","K0 cut monitor mc",35,0.5,35.5);
  
    //---------------Lambda
    fHistPiPMass[0] = new TH1F("fHistPiPMass"," p+pi- InvMass distribution",2*nbMass,0.,2.);
    fHistPiPMassVSPt[0] = new TH2F("fHistPiPMassVSPt","p+pi- InvMass distribution",nbMass,1.05,1.25,300,0.0,30.0);
    fHistPiPMassVSPtMCTruth[0] = new TH2F("fHistPiPMassVSPtMCTruth","p+pi- InvMass distribution vs pt MCTruth",nbMass,1.05,1.25,300,0.0,30.0);
    fHistPiPMassVSY[0] = new TH2F("fHistPiPMassVSY","p+pi- InvMass distribution vs rapidity",nbMass,1.05,1.25,200,-1.0,1.0);
    fHistPiPPtVSY[0] = new TH2F("fHistPiPPtVSY","p{t} vs y",100,-1,1,100,0.0,20);
    fHistPiPDecayLengthVsPt[0] = new TH2F("fHistPiPDecayLengthVsPt","#Lambda decay length vs pt",200,0.0,20.0,220,0.0,110.0);
    fHistPiPDecayLengthVsMass[0] = new TH2F("fHistPiPDecayLengthVsMass","#Lambda decay length vs mass",nbMass,1.05,1.25,220,0.0,110.0);
    if(!fSetPtDepHist){
      fHistPiPDecayLengthVsCtau[0] = new TH2F("fHistPiPDecayLengthVsCtau","L ctau vs mass",nbMass,1.05,1.25,250,0.0,50.0);
    }
    else{
      fHistPiPDecayLengthVsCtau[0] = new TH2F("fHistPiPDecayLengthVsCtau","L ctau vs pt",200,0.0,20.0,250,0.0,50.0);
    }
    fHistPiPMonitorCuts[0] = new TH1F("fHistPiPMonitorCuts","#Lambda cut monitor",35,0.5,35.5);
    fHistPiPMonitorMCCuts[0] = new TH1F("fHistPiPMonitorMCCuts","#Lambda cut monitor mc ",35,0.5,35.5);
  
    //-------------ALamda
    fHistPiAPMass[0] = new TH1F("fHistPiAPMass"," ap-pi+ InvMass distribution",2*nbMass,0.,2.);
    fHistPiAPMassVSPt[0] = new TH2F("fHistPiAPMassVSPt","p-pi+ InvMass distribution",nbMass,1.05,1.25,300,0.0,30.0);
    fHistPiAPMassVSPtMCTruth[0] = new TH2F("fHistPiAPMassVSPtMCTruth","p-pi+ InvMass distribution vs pt MCTruth",nbMass,1.05,1.25,300,0.0,30.0);
    fHistPiAPMassVSY[0] = new TH2F("fHistPiAPMassVSY","p-pi+ InvMass distribution vs rapidity",nbMass,1.05,1.25,200,-1.0,1.0);
    fHistPiAPPtVSY[0] = new TH2F("fHistPiAPPtVSY","p{t} vs y",100,-1,1,100,0.0,20);
    fHistPiAPDecayLengthVsPt[0] = new TH2F("fHistPiAPDecayLengthVsPt","#bar{#Lambda} decay length vs pt",200,0.0,20.0,220,0.0,110.0);
    fHistPiAPDecayLengthVsMass[0] = new TH2F("fHistPiAPDecayLengthVsMass","#bar{#Lambda} decay length vs mass",nbMass,1.05,1.25,220,0.0,110.0);
    if(!fSetPtDepHist){
      fHistPiAPDecayLengthVsCtau[0] = new TH2F("fHistPiAPDecayLengthVsCtau","AL ctau vs mass",nbMass,1.05,1.25,250,0.0,50.0);
    }
    else{
      fHistPiAPDecayLengthVsCtau[0] = new TH2F("fHistPiAPDecayLengthVsCtau","AL ctau vs pt",200,0.0,20.0,250,0.0,50.0);
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

      //K0s
      //--------------- Lambda
      fHistPiPMass[1] = new TH1F("fHistPiPMassSec"," p+pi- InvMass distribution",2*nbMass,0.,2.);
      fHistPiPMassVSPt[1] = new TH2F("fHistPiPMassVSPtSec","p+pi- InvMass distribution",nbMass,1.05,1.25,300,0.0,30.0);
      fHistPiPMassVSPtMCTruth[1] = new TH2F("fHistPiPMassVSPtMCTruthSec","p+pi- InvMass distribution vs pt MCTruth",nbMass,1.05,1.25,300,0.0,30.0);
      fHistPiPMassVSY[1] = new TH2F("fHistPiPMassVSYSec","p+pi- InvMass distribution vs rapidity",nbMass,1.05,1.25,200,-1.0,1.0);   
      fHistPiPPtVSY[1] = new TH2F("fHistPiPPtVSYSec","p{t} vs y",100,-1,1,100,0.0,20);
      fHistPiPDecayLengthVsPt[1] = new TH2F("fHistPiPDecayLengthVsPtSec","#Lambda decay length vs pt",200,0.0,20.0,220,0.0,110.0);
      fHistPiPDecayLengthVsMass[1] = new TH2F("fHistPiPDecayLengthVsMassSec","#Lambda decay length vs mass",nbMass,1.05,1.25,220,0.0,110.0);
      if(!fSetPtDepHist){
	fHistPiPDecayLengthVsCtau[1] = new TH2F("fHistPiPDecayLengthVsCtauSec","L ctau vs mass",nbMass,1.05,1.25,250,0.0,50.0);
      }
      else{
	fHistPiPDecayLengthVsCtau[1] = new TH2F("fHistPiPDecayLengthVsCtauSec","L ctau vs pt",200,0.0,20.0,250,0.0,50.0);
      }
      fHistPiPMonitorCuts[1] = new TH1F("fHistPiPMonitorCutsSec","#Lambda cut monitor",35,0.5,35.5);
      fHistPiPMonitorMCCuts[1] = new TH1F("fHistPiPMonitorMCCutsSec","#Lambda cut monitor mc",35,0.5,35.5);

      //----------------ALambda
      fHistPiAPMass[1] = new TH1F("fHistPiAPMassSec"," ap-pi+ InvMass distribution",2*nbMass,0.,2.);
      fHistPiAPMassVSPt[1] = new TH2F("fHistPiAPMassVSPtSec","p-pi+ InvMass distribution",nbMass,1.05,1.25,300,0.0,30.0);
      fHistPiAPMassVSPtMCTruth[1] = new TH2F("fHistPiAPMassVSPtMCTruthSec","p-pi+ InvMass distribution vs pt MCTruth",nbMass,1.05,1.25,300,0.0,30.0);
      fHistPiAPMassVSY[1] = new TH2F("fHistPiAPMassVSYSec","p-pi+ InvMass distribution vs rapidity",nbMass,1.05,1.25,200,-1.0,1.0);
      fHistPiAPPtVSY[1] = new TH2F("fHistPiAPPtVSYSec","p{t} vs y",100,-1,1,100,0.0,20);
      fHistPiAPDecayLengthVsPt[1] = new TH2F("fHistPiAPDecayLengthVsPtSec","#bar{#Lambda} decay length vs pt",200,0.0,20.0,220,0.0,110.0);
      fHistPiAPDecayLengthVsMass[1] = new TH2F("fHistPiAPDecayLengthVsMassSec","#bar{#Lambda} decay length vs mass",nbMass,1.05,1.25,220,0.0,110.0);
      if(!fSetPtDepHist){
	fHistPiAPDecayLengthVsCtau[1] = new TH2F("fHistPiAPDecayLengthVsCtauSec","AL ctau vs mass",nbMass,1.05,1.25,250,0.0,50.0);
      }
      else{
	fHistPiAPDecayLengthVsCtau[1] = new TH2F("fHistPiAPDecayLengthVsCtauSec","AL ctau vs pt",200,0.0,20.0,250,0.0,50.0);
      }

      fHistPiAPMonitorCuts[1] = new TH1F("fHistPiAPMonitorCutsSec","#bar{#Lambda} cut monitor",35,0.5,35.5);
      fHistPiAPMonitorMCCuts[1] = new TH1F("fHistPiAPMonitorMCCutsSec","#bar{#Lambda} cut monitor mc",35,0.5,35.5);
    }

    //add
    //------------ K0s ------------------
    fOutputContainer->Add(fHistPiPiMass);	 
    fOutputContainer->Add(fHistPiPiMassVSPt);
    fOutputContainer->Add(fHistPiPiMassVSPtMCTruth);
    fOutputContainer->Add(fHistPiPiMassVSY);
    fOutputContainer->Add(fHistPiPiPtVSY);
    fOutputContainer->Add(fHistPiPiDecayLengthVsPt);
    fOutputContainer->Add(fHistPiPiDecayLengthVsCtau);
    fOutputContainer->Add(fHistPiPiDecayLengthVsMass);  
    fOutputContainer->Add(fHistPiPiMonitorCuts);
    fOutputContainer->Add(fHistPiPiMonitorMCCuts);
    fOutputContainer->Add(fHistPiPiK0sVsLambdaMass);
    fOutputContainer->Add(fHistPiPiK0sVsALambdaMass);
    // fOutputContainer->Add(fHistPiPiPhiPosVsPtPosVsMass);//xxx

    // --------------- Lambda ---------------
    fOutputContainer->Add(fHistPiPK0sVsLambdaMass);
    fOutputContainer->Add(fHistPiPALambdaVsLambdaMass);
    // fOutputContainer->Add(fHistPiPPhiPosVsPtPosVsMass);//xxx

    // --------------- ALambda ---------------
    fOutputContainer->Add(fHistPiAPK0sVsALambdaMass);
    fOutputContainer->Add(fHistPiAPLambdaVsALambdaMass);

  
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
      fOutputContainer->Add(fHistPiAPMassVSPtMCTruth[j]);
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
	
      //K0
  
      // fHistPiPiMassVSAlpha = new TH2F("fHistPiPiMassVSAlpha"," alpha armenteros vs pi+pi- InvMass distribution",nbMass,0.25,0.75,500,-1.,1.);
      if(!fSetPtDepHist){
	fHistPiPiDCADaughters = new TH2F("fHistPiPiDCADaughters","dca of K0 daughters",nbMass,0.25,0.75,250,0.0,2);
	fHistPiPiDCADaughterPosToPrimVtxVSMass = new TH2F("fHistPiPiDCADaughterPosToPrimVtxVSMass","pi+ DCA daughter to prim vtx vsinvmass",nbMass,0.25,0.75,250,0.0,10.0);
	fHistPiPiDCAVSMass = new TH2F("fHistPiPiDCAVSMass","pi+pi- dca  vs pt",nbMass,0.25,0.75,250,0.0,5.0);
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
	fHistPiPiCosPointAng  = new TH2F("fHistPiPiCosPointAng","K0 cosine of pointing angle vs mass ",200,0.0,20.0,200,0.99,1.00);
	fHistPiPiRadiusXY = new TH2F("fHistPiPiRadiusXY","pi+pi- phi dist vs mass",200,0.0,20.0,200,0.0,4.0);
	fHistPiPiDCAZPos = new TH2F("fHistPiPiDCAZPos","dca z  of K0 pos daughters",200,0.0,20.0,200,-20.0,20.0);
	fHistPiPiDCAZNeg = new TH2F("fHistPiPiDCAZNeg","dca z  of K0 neg daughters",200,0.0,20.0,200,-20.0,20.0);
	fHistPiPiTrackLengthPosVsMass = new TH2F("fHistPiPiTrackLengthPosVsMass","track lenght of pos K0s daughter in TPC",200,0.0,20.0,250,0.0,250.0);
	fHistPiPiTrackLengthNegVsMass = new TH2F("fHistPiPiTrackLengthNegVsMass","track lenght of neg K0s daughter in TPC",200,0.0,20.0,250,0.0,250.0);
      }

      //Lambda
      if(!fSetPtDepHist){
	fHistPiPDCADaughters[0] = new TH2F("fHistPiPDCADaughters","dca of #Lambda daughters",nbMass,1.05,1.25,250,0.0,2.0);
	fHistPiPDCADaughterPosToPrimVtxVSMass[0] = new TH2F("fHistPiPDCADaughterPosToPrimVtxVSMass","pos DCA daughter to prim vtx vs invmass",nbMass,1.05,1.25,250,0.0,10.0);
	fHistPiPDCADaughterNegToPrimVtxVSMass[0] = new TH2F("fHistPiPDCADaughterNegToPrimVtxVSMass","neg DCA daughter to prim vtx vs invmass",nbMass,1.05,1.25,250,0.0,10.0);
	fHistPiPDCAVSMass[0] = new TH2F("fHistPiPDCAVSMass","ppi- dca  vs pt",nbMass,1.05,1.25,250,0.0,5.0);
	fHistPiPCosPointAng[0]  = new TH2F("fHistPiPCosPointAng","#Lambda cosine of pointing angle vs mass ",nbMass,1.05,1.25,200,0.99,1.00);
	fHistPiPRadiusXY[0] = new TH2F("fHistPiPRadiusXY","pi-p+ phi dist vs mass",nbMass,1.05,1.25,200,0.0,4.0);
	// fHistPiPPtDaughters[0] = new TH2F("fHistPiPPtDaughters","p_{t} pos vs p_{t} neg of daughters",400,0.0,20.0,400,0,20.0);
	//	fHistPiPDCAZPos[0] = new TH2F("fHistPiPDCAZPos","dca z  of Lambda pos daughters",nbMass,1.05,1.25,200,-20.0,20.0);
	//fHistPiPDCAZNeg[0] = new TH2F("fHistPiPDCAZNeg","dca z  of Lambda neg daughters",nbMass,1.05,1.25,200,-20.0,20.0);
	fHistPiPTrackLengthPosVsMass[0] = new TH2F("fHistPiPTrackLengthPosVsMass","track length of pos Lambda daughter in TPC",nbMass,1.05,1.25,250,0.0,250.0);
	fHistPiPTrackLengthNegVsMass[0] = new TH2F("fHistPiPTrackLengthNegVsMass","track length of neg Lambda daughter in TPC",nbMass,1.05,1.25,250,0.0,250.0);
      }
      else{//pt dependence
	fHistPiPDCADaughters[0] = new TH2F("fHistPiPDCADaughters","dca of #Lambda daughters",200,0.0,20.0,250,0.0,2.0);
	fHistPiPDCADaughterPosToPrimVtxVSMass[0] = new TH2F("fHistPiPDCADaughterPosToPrimVtxVSMass","pos DCA daughter to prim vtx vs invmass",200,0.0,20.0,250,0.0,10.0);
	fHistPiPDCADaughterNegToPrimVtxVSMass[0] = new TH2F("fHistPiPDCADaughterNegToPrimVtxVSMass","neg DCA daughter to prim vtx vs invmass",200,0.0,20.0,250,0.0,10.0);
	fHistPiPDCAVSMass[0] = new TH2F("fHistPiPDCAVSMass","ppi- dca  vs pt",200,0.0,20.0,250,0.0,5.0);
	fHistPiPCosPointAng[0]  = new TH2F("fHistPiPCosPointAng","#Lambda cosine of pointing angle vs mass ",200,0.0,20.0,200,0.99,1.00);
	fHistPiPRadiusXY[0] = new TH2F("fHistPiPRadiusXY","pi-p+ phi dist vs mass",200,0.0,20.0,200,0.0,4.0);
	//fHistPiPDCAZPos[0] = new TH2F("fHistPiPDCAZPos","dca z  of Lambda pos daughters",200,0.0,20.0,200,-20.0,20.0);
	//fHistPiPDCAZNeg[0] = new TH2F("fHistPiPDCAZNeg","dca z  of Lambda neg daughters",200,0.0,20.0,200,-20.0,20.0);
	fHistPiPTrackLengthPosVsMass[0] = new TH2F("fHistPiPTrackLengthPosVsMass","track length of pos Lambda daughter in TPC",200,0.0,20.0,250,0.0,250.0);
	fHistPiPTrackLengthNegVsMass[0] = new TH2F("fHistPiPTrackLengthNegVsMass","track length of neg Lambda daughter in TPC",200,0.0,20.0,250,0.0,250.0);
      }

      //AntiLambda
      if(!fSetPtDepHist){
	fHistPiAPDCADaughters[0] = new TH2F("fHistPiAPDCADaughters","dca of #bar{#Lambda} daughters",nbMass,1.05,1.25,250,0.0,2.0);
	fHistPiAPDCADaughterPosToPrimVtxVSMass[0] = new TH2F("fHistPiAPDCADaughterPosToPrimVtxVSMass","pos DCA daughter to prim vtx vs invmass",nbMass,1.05,1.25,250,0.0,10.0);
	fHistPiAPDCADaughterNegToPrimVtxVSMass[0] = new TH2F("fHistPiAPDCADaughterNegToPrimVtxVSMass","neg DCA daughter to prim vtx vs invmass",nbMass,1.05,1.25,250,0.0,10.0);
	fHistPiAPDCAVSMass[0] = new TH2F("fHistPiAPDCAVSMass","pi+p- dca  vs pt",nbMass,1.05,1.25,250,0.0,5.0);
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
	fHistPiAPCosPointAng[0] = new TH2F("fHistPiAPCosPointAng","#bar{#Lambda} cosine of pointing angle vs mass",200,0.0,20.0,200,0.99,1.00);
	fHistPiAPRadiusXY[0] = new TH2F("fHistPiAPRadiusXY","pi+p- phi dist vs mass",200,0.0,20.0,200,0.0,4.0);
	//	fHistPiAPDCAZPos[0] = new TH2F("fHistPiAPDCAZPos","dca z  of ALambda pos daughters",200,0.0,20.0,200,-20.0,20.0);
	//fHistPiAPDCAZNeg[0] = new TH2F("fHistPiAPDCAZNeg","dca z  of ALambda neg daughters",200,0.0,20.0,200,-20.0,20.0);
	fHistPiAPTrackLengthPosVsMass[0] = new TH2F("fHistPiAPTrackLengthPosVsMass","track length of pos ALambda daughter in TPC",200,0.0,20.0,250,0.0,250.0);
	fHistPiAPTrackLengthNegVsMass[0] = new TH2F("fHistPiAPTrackLengthNegVsMass","track length of neg ALambda daughter in TPC",200,0.0,20.0,250,0.0,250.0);
      }
   
      //dedx
      fHistDedxSecProt[0] = new TH2F("fHistDedxSecProt","proton", nbPt, 0, 20, 100, 0, 400);
      fHistDedxSecPiPlus[0] = new TH2F("fHistDedxSecPiPlus","pi plus", nbPt, 0, 20, 100, 0, 400);
      fHistDedxSecAProt[0] = new TH2F("fHistDedxSecAProt","antiproton", nbPt, 0, 20, 100, 0, 400);
      fHistDedxSecPiMinus[0] = new TH2F("fHistDedxSecPiMinus","pi minus", nbPt, 0, 20, 100, 0, 400);
      fHistDedxProt[0] = new TH2F("fHistDedxProt","proton", nbPt, 0, 20, 100, 0, 400);
      fHistDedxPiPlus[0] = new TH2F("fHistDedxPiPlus","pi plus", nbPt, 0, 20, 100, 0, 400);
      fHistDedxAProt[0] = new TH2F("fHistDedxAProt","antiproton", nbPt, 0, 20, 100, 0, 400);
      fHistDedxPiMinus[0] = new TH2F("fHistDedxPiMinus","pi minus", nbPt, 0, 20, 100, 0, 400);


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

	//K0
	//Lambda
	if(!fSetPtDepHist){
	  fHistPiPDCADaughters[1] = new TH2F("fHistPiPDCADaughtersSec","dca of #Lambda daughters",nbMass,1.05,1.25,250,0.0,2.0);
	  fHistPiPDCADaughterPosToPrimVtxVSMass[1] = new TH2F("fHistPiPDCADaughterPosToPrimVtxVSMassSec","pos sec DCA daughter to prim vtx vs invmass",nbMass,1.05,1.25,250,0.0,10.0);
	  fHistPiPDCADaughterNegToPrimVtxVSMass[1] = new TH2F("fHistPiPDCADaughterNegToPrimVtxVSMassSec","neg sec DCA daughter to prim vtx vs invmass",nbMass,1.05,1.25,250,0.0,10.0);
	  fHistPiPDCAVSMass[1] = new TH2F("fHistPiPDCAVSMassSec","ppi- dca  vs pt",nbMass,1.05,1.25,250,0.0,5.0);
	  fHistPiPCosPointAng[1]  = new TH2F("fHistPiPCosPointAngSec","#Lambda cosine of pointing angle vs mass",nbMass,1.05,1.25,200,0.99,1.00);
	  //	 fHistPiPDecayLengthVsMass[1] = new TH2F("fHistPiPDecayLengthVsMassSec","#Lambda decay length vs mass",nbMass,1.05,1.25,200,0.0,100.0);
	  fHistPiPRadiusXY[1] = new TH2F("fHistPiPRadiusXYSec","pi-p+ phi dist vs mass",nbMass,1.05,1.25,200,0.0,4.0);
	  // fHistPiPPtDaughters[0] = new TH2F("fHistPiPPtDaughters","p_{t} pos vs p_{t} neg of daughters",400,0.0,20.0,400,0,20.0);
	  //	  fHistPiPDCAZPos[1] = new TH2F("fHistPiPDCAZPosSec","dca z  of Lambda sec pos daughters",nbMass,1.05,1.25,200,-20.0,20.0);
	  //	  fHistPiPDCAZNeg[1] = new TH2F("fHistPiPDCAZNegSec","dca z  of Lambda sec neg daughters",nbMass,1.05,1.25,200,-20.0,20.0);
	  fHistPiPTrackLengthPosVsMass[1] = new TH2F("fHistPiPTrackLengthPosVsMassSec","track length of pos sec Lambda daughter in TPC",nbMass,1.05,1.25,250,0.0,250.0);
	  fHistPiPTrackLengthNegVsMass[1] = new TH2F("fHistPiPTrackLengthNegVsMassSec","track length of neg sec Lambda daughter in TPC",nbMass,1.05,1.25,250,0.0,250.0);
	}
	else{
	  fHistPiPDCADaughters[1] = new TH2F("fHistPiPDCADaughtersSec","dca of #Lambda daughters",200,0.0,20.0,250,0.0,2.0);
	  fHistPiPDCADaughterPosToPrimVtxVSMass[1] = new TH2F("fHistPiPDCADaughterPosToPrimVtxVSMassSec","pos sec DCA daughter to prim vtx vs invmass",200,0.0,20.0,250,0.0,10.0);
	  fHistPiPDCADaughterNegToPrimVtxVSMass[1] = new TH2F("fHistPiPDCADaughterNegToPrimVtxVSMassSec","neg sec DCA daughter to prim vtx vs invmass",200,0.0,20.0,250,0.0,10.0);
	  fHistPiPDCAVSMass[1] = new TH2F("fHistPiPDCAVSMassSec","ppi- dca  vs pt",200,0.0,20.0,250,0.0,5.0);
	  fHistPiPCosPointAng[1]  = new TH2F("fHistPiPCosPointAngSec","#Lambda cosine of pointing angle vs mass",200,0.0,20.0,200,0.99,1.00);
	  fHistPiPRadiusXY[1] = new TH2F("fHistPiPRadiusXYSec","pi-p+ phi dist vs mass",200,0.0,20.0,200,0.0,4.0);
	  //	  fHistPiPDCAZPos[1] = new TH2F("fHistPiPDCAZPosSec","dca z  of Lambda sec pos daughters",200,0.0,20.0,200,-20.0,20.0);
	  //fHistPiPDCAZNeg[1] = new TH2F("fHistPiPDCAZNegSec","dca z  of Lambda sec neg daughters",200,0.0,20.0,200,-20.0,20.0);
	  fHistPiPTrackLengthPosVsMass[1] = new TH2F("fHistPiPTrackLengthPosVsMassSec","track length of pos sec Lambda daughter in TPC",200,0.0,20.0,250,0.0,250.0);
	  fHistPiPTrackLengthNegVsMass[1] = new TH2F("fHistPiPTrackLengthNegVsMassSec","track length of neg sec Lambda daughter in TPC",200,0.0,20.0,250,0.0,250.0);
	}
	  
	//ALambda
	if(!fSetPtDepHist){
	  fHistPiAPDCADaughters[1] = new TH2F("fHistPiAPDCADaughtersSec","dca of #bar{#Lambda} daughters",nbMass,1.05,1.25,250,0.0,2.0);
	  fHistPiAPDCADaughterPosToPrimVtxVSMass[1] = new TH2F("fHistPiAPDCADaughterPosToPrimVtxVSMassSec","pos sec DCA daughter to prim vtx vs invmass",nbMass,1.05,1.25,250,0.0,10.0);
	  fHistPiAPDCADaughterNegToPrimVtxVSMass[1] = new TH2F("fHistPiAPDCADaughterNegToPrimVtxVSMassSec","neg sec DCA daughter to prim vtx vs invmass",nbMass,1.05,1.25,250,0.0,10.0);
	  fHistPiAPDCAVSMass[1]   = new TH2F("fHistPiAPDCAVSMassSec","pi+p- dca  vs pt",nbMass,1.05,1.25,250,0.0,5.0);
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
	  fHistPiAPCosPointAng[1] = new TH2F("fHistPiAPCosPointAngSec","#bar{#Lambda} cosine of pointing angle vs mass",200,0.0,20.0,200,0.99,1.00);
	  fHistPiAPRadiusXY[1] = new TH2F("fHistPiAPRadiusXYSec","pi+p- phi dist vs mass",200,0.0,20.0,200,0.0,4.0);
	  //	  fHistPiAPDCAZPos[1] = new TH2F("fHistPiAPDCAZPosSec","dca z  of ALambda sec pos daughters",200,0.0,20.0,200,-20.0,20.0);
	  //fHistPiAPDCAZNeg[1] = new TH2F("fHistPiAPDCAZNegSec","dca z  of ALambda sec neg daughters",200,0.0,20.0,200,-20.0,20.0);
	  fHistPiAPTrackLengthPosVsMass[1] = new TH2F("fHistPiAPTrackLengthPosVsMassSec","track length of pos sec ALambda daughter in TPC",200,0.0,20.0,250,0.0,250.0);
	  fHistPiAPTrackLengthNegVsMass[1] = new TH2F("fHistPiAPTrackLengthNegVsMassSec","track length of neg sec ALambda daughter in TPC",200,0.0,20.0,250,0.0,250.0);
	}

	//dedx
	fHistDedxSecProt[1] = new TH2F("fHistDedxSecProtSec","proton", nbPt, 0, 20, 100, 0, 400);
	fHistDedxSecPiPlus[1] = new TH2F("fHistDedxSecPiPlusSec","pi plus", nbPt, 0, 20, 100, 0, 400);
	fHistDedxSecAProt[1] = new TH2F("fHistDedxSecAProtSec","antiproton", nbPt, 0, 20, 100, 0, 400);
	fHistDedxSecPiMinus[1] = new TH2F("fHistDedxSecPiMinusSec","pi minus", nbPt, 0, 20, 100, 0, 400);
	fHistDedxProt[1] = new TH2F("fHistDedxProtSec","proton", nbPt, 0, 20, 100, 0, 400);
	fHistDedxPiPlus[1] = new TH2F("fHistDedxPiPlusSec","pi plus", nbPt, 0, 20, 100, 0, 400);
	fHistDedxAProt[1] = new TH2F("fHistDedxAProtSec","antiproton", nbPt, 0, 20, 100, 0, 400);
	fHistDedxPiMinus[1] = new TH2F("fHistDedxPiMinusSec","pi minus", nbPt, 0, 20, 100, 0, 400);

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

      //------ ITS TPC clusters --------------
      fOutputContainer->Add(fHistNclsITS[0]) ;
      fOutputContainer->Add(fHistNclsTPC[0]);
      fOutputContainer->Add(fHistNclsITS[1]);
      fOutputContainer->Add(fHistNclsTPC[1]);

      //-----------K0s ------------------
      fOutputContainer->Add(fHistPiPiDCAZNeg);
      fOutputContainer->Add(fHistPiPiDCAZPos);
      fOutputContainer->Add(fHistPiPiDCADaughters); 
      fOutputContainer->Add(fHistPiPiDCADaughterPosToPrimVtxVSMass);
      fOutputContainer->Add(fHistPiPiDCAVSMass);
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

      //----------- Lambda Antilambda -------------
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
	fOutputContainer->Add(fHistPiAPDCAVSMass[j]);
	fOutputContainer->Add(fHistPiPCosPointAng[j]);
	fOutputContainer->Add(fHistPiAPCosPointAng[j]);
	//fOutputContainer->Add(fHistPiPDCAZNeg[j]);
	//fOutputContainer->Add(fHistPiPDCAZPos[j]);
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

	//--------- TPC Lambda-----------------
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
      fHistPiPiPDGCode = new TH1F("fHistPiPiPDGCode","PDG code of K0s mothers",3500,0,3500);
      fOutputContainer->Add(fHistPiPiPDGCode);
      fHistPiPPDGCode = new TH1F("fHistPiPPDGCode","PDG code of #Lambda  mothers",3500,0,3500);
      fOutputContainer->Add(fHistPiPPDGCode);
      fHistPiAPPDGCode = new TH1F("fHistPiAPPDGCode","PDG code of #bar{#Lambda} mothers",3500,0,3500);
      fOutputContainer->Add(fHistPiAPPDGCode);
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

      //-------------K0s
   
      fHistPiPiDecayLengthResolution = new TH2F("fHistPiPiDecayLengthResolution","K0s decay length resolution MC",220,0.0,110.0,220,0.0,110);
	 
      //-------------Lambda
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

      //--------------ALambda
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
    fESD = esdH->GetEvent();
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
  
     
  
    // -- Check fo centrality
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

      //ask for pileup from SPD
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
      
      //vertex selection
      if (vtxESD->GetStatus()){
	fHistNEvents->Fill(2);
	fHistESDVertexZ->Fill(vtxESD->GetZv());
	if(fabs(vtxESD->GetZv()) < fVertexZCut){
	  fHistMuliplicityRaw->Fill(multV0);
	  fHistNEvents->Fill(3);
	  fHistNPrim->Fill(nContr);
	
	  Process();
	
	  fHistMuliplicity->Fill(multV0);
	
	  nContr = vtxESD->GetNContributors();
	  //  if(nContr<501){
	  fHistPrimVtxZESDVSNContributors->Fill(vtxESD->GetZv(),nContr);
	  fHistPrimVtxZESDTPCVSNContributors->Fill(vtxESDTPC->GetZv(),nContr);
	  //fHistPrimVtxZESDSPDVSNContributorsTPC->Fill(vtxESDSPD->GetZv(),nContr);
	  //   }
	  fHistPrimVtxZESD->Fill(vtxESD->GetZv());
	  fHistPrimVtxZESDTPC->Fill(vtxESDTPC->GetZv());
	  // fHistPrimVtxZESDSPD->Fill(vtxESDSPD->GetZv());
	  // -- count events after processing
	  fHistNEvents->Fill(4);
	}
      }
      else{
	if(vtxESDSPD->GetStatus()){
	  fHistNEvents->Fill(2);
	
	  fHistESDVertexZ->Fill(vtxESDSPD->GetZv());
	  if(fabs(vtxESDSPD->GetZv()) < fVertexZCut){
	  
	    fHistMuliplicityRaw->Fill(multV0);
	    fHistNEvents->Fill(3);
	    fHistNPrim->Fill(nContr);
	  
	    Process();
	  
	    fHistMuliplicity->Fill(multV0);
	  
	    nContr = vtxESDSPD->GetNContributors();
	    //  if(nContr<501){
	    //fHistPrimVtxZESDVSNContributors->Fill(vtxESD->GetZv(),nContr);
	    fHistPrimVtxZESDTPCVSNContributors->Fill(vtxESDTPC->GetZv(),nContr);
	    fHistPrimVtxZESDSPDVSNContributors->Fill(vtxESDSPD->GetZv(),nContr);
	    // }
	    // fHistPrimVtxZESD->Fill(vtxESD->GetZv());
	    fHistPrimVtxZESDTPC->Fill(vtxESDTPC->GetZv());
	    fHistPrimVtxZESDSPD->Fill(vtxESDSPD->GetZv());
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
	Double_t vtxZ = vtxESD->GetZv();
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
  void AliAnalysisTaskV0ForRAA::Process(){
    //run the analysis

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
    //histo fo user defined shift in charge/pt 
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
      
      //-- phi --//
      Double_t pi = TMath::Pi();
      Double_t phi = p0->Phi();
      if(phi>pi) phi -=2*pi;

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
   
      if(!isPrim){//prim or sec
	isSecd=1;// is secondary V0s
	//    cout<<notinjected<<endl;
	if(indexMother1 >-1){// && !isPrim){//secondary V0s
	  //     isSecd=1;// is secondary V0s
	  //   cout<<fMCev->IsFromBGEvent(indexMother1)<<endl;
	  // if(fSelectMBMotherMC && !fMCev->IsFromBGEvent(indexMother1)) continue;//xxx only standard hijing particles for sec. lambdas:not needed
	  //-- check for mother --//
	  TParticle *mother = stack->Particle(indexMother1);
	  if(!mother) {
	    Printf("no mother pointer!");continue;
	  }
	  pdgMother = mother->GetPdgCode();
	  fHistPiPMonitorMCCuts[1]->Fill(11*fillFlagL);
	  fHistPiAPMonitorMCCuts[1]->Fill(11*fillFlagAL);

	  //-- check for injejcted --//
	  Bool_t notinjectedMother = kTRUE;
	  notinjectedMother = fMCev->IsFromBGEvent(indexMother1);
	
	  if(fSelectInjected && !notinjectedMother ) continue;
	  fHistPiPMonitorMCCuts[1]->Fill(10*fillMCtruth);
	  fHistPiAPMonitorMCCuts[1]->Fill(10*fillMCtruth);

	  Bool_t isPrimMother= stack->IsPhysicalPrimary(indexMother1);
	  if(!isPrimMother) continue;
	  fHistPiPMonitorMCCuts[1]->Fill(12*fillFlagL);
	  fHistPiAPMonitorMCCuts[1]->Fill(12*fillFlagAL);
      
	  uniqueIDmother =  mother->GetUniqueID();

	  if(uniqueIDmother==13){
	    //	goodMother =0;
	    continue;
	  }
	  fHistPiPMonitorMCCuts[1]->Fill(13*fillFlagL);
	  fHistPiAPMonitorMCCuts[1]->Fill(13*fillFlagAL);


	  //-- fill secondary V0s histos and pdg histos --// 
	  ptXiMother=mother->Pt();
	  rapXiMother=mother->Y();
	

	  //-- K0s --//
	  if(pdgCode==310){
	    if(fabs(pdgMother)==311 || fabs(pdgMother)==313 || fabs(pdgMother)==323 ) isSecd=0; // from K0L,  K0 and K* as primary
	    else fHistPiPiPDGCode->Fill(fabs(pdgMother));
	    //goodMother=1;
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
		//goodMother=1;
	      }
	   
	    if( pdgMother == 3322) //xi0
	      {
		if(!fRapCutV0 || fabs(rapidity)<fRap){
		  if(!fEtaCutMCDaughters  ||  (fabs(etaMC00)<fEtaCutMCDaughtersVal|| fabs(etaMC01)<fEtaCutMCDaughtersVal)){
		    fHistPiPMassVSPtSecXi[0]->Fill(massV0MC,ptV0MC);
		    fHistPiPMassVSYSecXi[0]->Fill(massV0MC,rapidity);
		    if(!fRapCutV0 || fabs(rapXiMother)<fRap) fHistPiPXi0PtVSLambdaPt[0]->Fill(ptV0MC,ptXiMother);
		  }
		}
		//goodMother=1;
	      }
	    if(pdgMother == 3312) //xi minus
	      {
		if(!fRapCutV0 || fabs(rapidity)<fRap){
		  if(!fEtaCutMCDaughters  ||  (fabs(etaMC00)<fEtaCutMCDaughtersVal|| fabs(etaMC01)<fEtaCutMCDaughtersVal)){
		    fHistPiPMassVSPtSecXi[0]->Fill(massV0MC,ptV0MC);
		    fHistPiPMassVSYSecXi[0]->Fill(massV0MC,rapidity);
		    if(!fRapCutV0 || fabs(rapXiMother)<fRap)  fHistPiPXiMinusPtVSLambdaPt[0]->Fill(ptV0MC,ptXiMother);
		  }
		}
		//goodMother=1;

	      }	    
	    if(pdgMother == 3334)//omega-
	      {
		//  fHistPiPDCAtoPrimVtxOmega[0]->Fill(p0->GetMass(),dcaV0ToPrimVertex);
		if(!fRapCutV0 || fabs(rapXiMother)<fRap) fHistPiPOmegaPtVSLambdaPt[0]->Fill(ptV0MC,ptXiMother);
		fHistPiPMassVSPtSecOmega[0]->Fill(massV0MC,ptV0MC);
		//goodMother=1;
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
		//	    goodMother=1;
	      }
		   
	    if( pdgMother == -3322) //xi0
	      {
		if(!fRapCutV0 || fabs(rapidity)<fRap){
		  if(!fEtaCutMCDaughters  ||  (fabs(etaMC00)<fEtaCutMCDaughtersVal|| fabs(etaMC01)<fEtaCutMCDaughtersVal)){
		    fHistPiAPMassVSPtSecXi[0]->Fill(massV0MC,ptV0MC);
		    fHistPiAPMassVSYSecXi[0]->Fill(massV0MC,rapidity);
		    if(!fRapCutV0 || fabs(rapXiMother)<fRap) fHistPiAPXi0PtVSLambdaPt[0]->Fill(ptV0MC,ptXiMother);
		  }
		}
		// goodMother=1;
	      }
	    if(pdgMother == -3312) //xi plus
	      {
		if(!fRapCutV0 || fabs(rapidity)<fRap){
		  if(!fEtaCutMCDaughters  ||  (fabs(etaMC00)<fEtaCutMCDaughtersVal|| fabs(etaMC01)<fEtaCutMCDaughtersVal)){
		    fHistPiAPMassVSPtSecXi[0]->Fill(massV0MC,ptV0MC);
		    fHistPiAPMassVSYSecXi[0]->Fill(massV0MC,rapidity);
		    if(!fRapCutV0 || fabs(rapXiMother)<fRap) fHistPiAPXiMinusPtVSLambdaPt[0]->Fill(ptV0MC,ptXiMother);
		  }
		}
		//goodMother=1;
	      }
	    if(pdgMother == -3334)//omega+
	      {
		fHistPiAPOmegaPtVSLambdaPt[0]->Fill(ptV0MC,ptXiMother);
		if(!fRapCutV0 || fabs(rapXiMother)<fRap) fHistPiAPMassVSPtSecOmega[0]->Fill(massV0MC,ptV0MC);
		// fHistPiAPDCAtoPrimVtxOmega[0]->Fill(p0->GetMass(),dcaV0ToPrimVertex);
		// goodMother=1;
	      }
	  }
	}	
      }//end secondaries
      else{
	//-- check for injejcted --//
	Bool_t notinjected = kTRUE;
	notinjected = fMCev->IsFromBGEvent(iMc);
      
	if(fSelectInjected && !notinjected ) continue;
	fHistPiPiMonitorMCCuts->Fill(10*fillMCtruth);
	fHistPiPMonitorMCCuts[0]->Fill(10*fillMCtruth);
	fHistPiAPMonitorMCCuts[0]->Fill(10*fillMCtruth);
      }
      //else goodMother=1;
      // if(!goodMother) continue;// for (A)Lambda important

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
	Double_t  phiPosMC=0.0;
	    
	if(p00->GetPdgCode()<0)
	  {
	    vecPip.SetXYZ(p01->Px(),p01->Py(),p01->Pz());
	    vecPin.SetXYZ(p00->Px(),p00->Py(),p00->Pz());
	    ptMinus = pt00;
	    ptPlus = pt01;
	    phiPosMC = p01->Phi();
	  }
	else{
	  vecPin.SetXYZ(p01->Px(),p01->Py(),p01->Pz());
	  vecPip.SetXYZ(p00->Px(),p00->Py(),p00->Pz());
	  ptMinus = pt01;
	  ptPlus = pt00;
	  phiPosMC = p00->Phi();
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
	  fHistPiPiEtaDMC[isSecd]->Fill(etaMC00,ptV0MC);
	  fHistPiPiEtaDMC[isSecd]->Fill(etaMC01,ptV0MC);
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

    
	//      Double_t phiMC =   p0->Phi(); 
     

	//-- Fill Particle histos --//
	if (pdgCode==310){//K0s
	  fHistPiPiMonitorMCCuts->Fill(17);

	  fHistPiPiEtaDMC[1]->Fill(etaMC00,ptV0MC);
	  fHistPiPiEtaDMC[1]->Fill(etaMC01,ptV0MC);
	  
	  fHistPiPiMass->Fill(massV0MC);
	  fHistPiPiMassVSPt->Fill(massV0MC,ptV0MC);
	  fHistPiPiMassVSY->Fill(massV0MC,rapidity);
	  // fHistPiPiPtDaughters->Fill(ptMinus,ptPlus);
	  fHistPiPiPtVSY->Fill(rapidity,ptV0MC);
	   
	  Double_t ctTK0s=0.0,ctK0s=0.0;
	  if(pV0MC>0.0) ctK0s=declength3d*0.497614/pV0MC;
	  if(ptV0MC>0.0) ctTK0s=declength*0.497614/ptV0MC;
	  fHistPiPiDecayLengthResolution->Fill(declength3d,declength);
	  fHistPiPiDecayLengthVsPt->Fill(ptV0MC,declength);//ptV0MC,ctK0s);
	  fHistPiPiDecayLengthVsCtau->Fill(massV0MC,ctTK0s);
	  fHistPiPiDecayLengthVsMass->Fill(massV0MC,declength);
	  //all V0s histo
	  fHistArmenteros[isSecd]->Fill(alfa,qt);
	  fHistV0RadiusZ[isSecd]->Fill(rMC2D,xyzMC[2]);
	  fHistV0RadiusXY[isSecd]->Fill(xyzMC[0],xyzMC[1]);
	  fHistV0RadiusZVSPt[isSecd]->Fill(ptV0MC,xyzMC[2]);
	  //	fHistPiPiPhiPosVsPtPosVsMass->Fill(massV0MC,declength,ptV0MC);//,ctK0s);//phiPosMC);//xxx
	  fHistPiPiK0sVsLambdaMass->Fill(calcLambdamass,calcK0smass);
	  fHistPiPiK0sVsALambdaMass->Fill(calcALambdamass,calcK0smass);     
	}
	if (pdgCode==3122){ //Lambda
	  fHistPiPMonitorMCCuts[isSecd]->Fill(17);
	
	  fHistPiPEtaDMC[1]->Fill(etaMC00,ptV0MC);
	  fHistPiPEtaDMC[1]->Fill(etaMC01,ptV0MC);

	  fHistPiPMassVSPt[isSecd]->Fill(massV0MC,ptV0MC);
	  fHistPiPMass[isSecd]->Fill(massV0MC);  
	  fHistPiPMassVSY[isSecd]->Fill(massV0MC,rapidity);
	  //  fHistPiPPtDaughters[isSecd]->Fill(ptMinus,ptPlus);
	  fHistPiPPtVSY[isSecd]->Fill(rapidity,ptV0MC);
	   
	  Double_t ctTL=0.0, ctL=0.0;
	  if(pV0MC>0.0) ctL=declength3d*1.115683/pV0MC;
	  if(ptV0MC>0.0) ctTL=declength*1.115683/ptV0MC;
	  fHistPiPDecayLengthResolution[0]->Fill(declength3d,declength);
	  fHistPiPDecayLengthVsPt[isSecd]->Fill(ptV0MC,declength);//(ptV0MC,ctL);
	  fHistPiPDecayLengthVsCtau[isSecd]->Fill(massV0MC,ctTL);
	  fHistPiPDecayLengthVsMass[isSecd]->Fill(massV0MC,declength);
	  //all V0s hito	
	  fHistArmenteros[isSecd]->Fill(alfa,qt);
	  fHistV0RadiusZ[isSecd]->Fill(rMC2D,xyzMC[2]);
	  fHistV0RadiusXY[isSecd]->Fill(xyzMC[0],xyzMC[1]);
	  fHistV0RadiusZVSPt[isSecd]->Fill(ptV0MC,xyzMC[2]);
	  //	fHistPiPPhiPosVsPtPosVsMass->Fill(massV0MC,declength,ptV0MC);//,ctK0s);//phiPosMC);//xxx
	  fHistPiPK0sVsLambdaMass->Fill(calcLambdamass,calcK0smass);
	}
	if (pdgCode==-3122){ //AntiLambda
	  fHistPiAPMonitorMCCuts[isSecd]->Fill(17);
	    
	  fHistPiPEtaDMC[1]->Fill(etaMC00,ptV0MC);
	  fHistPiPEtaDMC[1]->Fill(etaMC01,ptV0MC);

	  fHistPiAPMassVSPt[isSecd]->Fill(massV0MC,ptV0MC);
	  fHistPiAPMass[isSecd]->Fill(massV0MC);
	  fHistPiPMassVSY[isSecd]->Fill(massV0MC,rapidity);
	  //  fHistPiAPPtDaughters[isSecd]->Fill(ptMinus,ptPlus);
	  fHistPiAPPtVSY[isSecd]->Fill(rapidity,ptV0MC);
	  
	  Double_t ctTL=0.0, ctL=0.0;
	  if(pV0MC>0.0) ctL=declength3d*1.115683/pV0MC;
	  if(ptV0MC>0.0) ctTL=declength*1.115683/ptV0MC;
	  fHistPiAPDecayLengthResolution[0]->Fill(declength3d,declength);
	  fHistPiAPDecayLengthVsPt[isSecd]->Fill(ptV0MC,declength);//(ptV0MC,ctAL);
	  fHistPiAPDecayLengthVsCtau[isSecd]->Fill(massV0MC,ctL);
	  fHistPiAPDecayLengthVsMass[isSecd]->Fill(massV0MC,ctTL);//declength);
	  //all V0s histo	   
	  fHistArmenteros[isSecd]->Fill(alfa,qt);
	  fHistV0RadiusZ[isSecd]->Fill(rMC2D,xyzMC[2]);
	  fHistV0RadiusXY[isSecd]->Fill(xyzMC[0],xyzMC[1]);
	  fHistV0RadiusZVSPt[isSecd]->Fill(ptV0MC,xyzMC[2]);
	  fHistV0RadiusXYVSY[isSecd]->Fill(rapidity,rMC2D);
	  fHistPiAPK0sVsALambdaMass->Fill(calcALambdamass,calcK0smass);
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
    primaryVtxPosition[0] = fESD->GetPrimaryVertex()->GetXv();
    primaryVtxPosition[1] = fESD->GetPrimaryVertex()->GetYv();
    primaryVtxPosition[2] = fESD->GetPrimaryVertex()->GetZv();
   
    Int_t nV0 = fESD->GetNumberOfV0s();
    AliESDv0 * v0MIs=NULL;


    //    Int_t on =0,off=0;
    Bool_t stopLoop = kFALSE;
    //------------------------ V0 reco loop --------------------//
    for(Int_t iV0MI = 0; iV0MI < nV0; iV0MI++) {//V0 loop
      //-- get V0 info --//
      // cout<<iV0MI<<endl;
      if(stopLoop) break;
	 
      v0MIs = fESD->GetV0(iV0MI);
      if(!v0MIs ) continue;

      // cout<<"after: "<<iV0MI<<endl;
      fHistPiPiMonitorCuts->Fill(1);
      fHistPiPMonitorCuts[isSecd]->Fill(1);
      fHistPiAPMonitorCuts[isSecd]->Fill(1);

      //------------ get references of daughters --------------//
      //-- esd tracks --//
      trackPosTest = fESD->GetTrack(v0MIs->GetPindex());
      trackNegTest = fESD->GetTrack(v0MIs->GetNindex());
     
      if ( trackPosTest->GetSign() == trackNegTest->GetSign()) continue;
	 
      fHistPiPiMonitorCuts->Fill(2);
      fHistPiPMonitorCuts[isSecd]->Fill(2);
      fHistPiAPMonitorCuts[isSecd]->Fill(2);

      Bool_t onthefly = v0MIs->GetOnFlyStatus();
      if(fOntheFly!=onthefly) continue;
      //  cout<<"onthefly "<<onthefly<<endl;
      // if(onthefly) on++;
      // else off++;
	  
      fHistPiPiMonitorCuts->Fill(3);
      fHistPiPMonitorCuts[isSecd]->Fill(3);
      fHistPiAPMonitorCuts[isSecd]->Fill(3);
         
      //-- for MC mode --//
      if(fMCMode){
	//check MC labels (and find partners for MC truth V0 daughters for fMCTruthMode=kTRUE)
	if(!GetMCTruthPartner(trackPosTest,trackNegTest,id0,id1)) continue;
	else stopLoop = kTRUE;
      }

      fHistPiPiMonitorCuts->Fill(4);
      fHistPiPMonitorCuts[isSecd]->Fill(4);
      fHistPiAPMonitorCuts[isSecd]->Fill(4);
     
      //-- onthefly selection --//
      //  Bool_t onthefly = v0MIs->GetOnFlyStatus();
      // 	 if(fOntheFly!=onthefly) continue;
	 
   
      // 	 fHistPiPiMonitorCuts[isSecd]->Fill(4);
      // 	 fHistPiPMonitorCuts[isSecd]->Fill(4);
      // 	 fHistPiAPMonitorCuts[isSecd]->Fill(4);
         
      //--  get eta from V0 daughters --//
      Double_t posDaughterEta=0.0;
      Double_t negDaughterEta=0.0;
      Double_t posDaughterPhi=0.0;
      Double_t negDaughterPhi=0.0;
	 
      Double_t eta00 = trackPosTest->Eta();
      Double_t eta01 = trackNegTest->Eta();
            

      Int_t indexPos = 0,indexNeg=0;
      //---------- check sign assignment for daughters --------//
      Bool_t switchSign = kFALSE;
     
      if( trackPosTest->GetSign() >0){//pos
	indexPos = v0MIs->GetPindex();
	indexNeg = v0MIs->GetNindex();

	//  trackPos =fESD->GetTrack(indexPos);
	//  trackNeg =fESD->GetTrack(indexNeg);

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

	// trackPos =fESD->GetTrack(indexPos);
	// trackNeg =fESD->GetTrack(indexNeg);
     
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
      //K0
      positivesMIPi.SetXYZM(pp[0],pp[1],pp[2],massPi);
      negativesMIPi.SetXYZM(pm[0],pm[1],pm[2],massPi);
      TLorentzVector v0K0=positivesMIPi+negativesMIPi;
    
      //Lambda
      positivesMIP.SetXYZM(pp[0],pp[1],pp[2],massP);
      TLorentzVector v0Lambda=positivesMIP+negativesMIPi;
     
      //Anitlambda
      negativesMIAP.SetXYZM(pm[0],pm[1],pm[2],massP);
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
      
      //      Double_t posDaughterPt = ppTrack.Pt();
      //      Double_t negDaughterPt = pmTrack.Pt();


      //	 Double_t phiPos = trackPos->Phi();
      //Double_t phiNeg = trackNeg->Phi();
      
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
     
      //------------------- DCA daughters ---------------------//
      v0MIs->GetXYZ(xr[0],xr[1],xr[2]);   
    
      //-- between the daughters --//
      Double_t dcaDaughters = v0MIs->GetDcaV0Daughters();  
    
      //-- to primary vertex --
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
      Bool_t checkProp = parPos->PropagateToDCA(fESD->GetPrimaryVertex(),fESD->GetMagneticField(),40.0,dcaYZP,covar);
      dcaToVertexZPos =  dcaYZP[1];
      delete parPos;
      checkProp = parNeg->PropagateToDCA(fESD->GetPrimaryVertex(),fESD->GetMagneticField(),40.0,dcaYZN,covar);
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
	 
      //dcaPosToVertex  =   posPKF->GetDistanceFromVertexXY(primaryVtxPosition);
      //dcaNegToVertex  =   negPiKF->GetDistanceFromVertexXY(primaryVtxPosition);
        
      //------------------- dca and decay radius V0 -------------//
   
       
      Double_t dim2V0Radius= sqrt(   pow(xr[0] - primaryVtxPosition[0],2.0)
				     +pow(xr[1] - primaryVtxPosition[1],2.0));//TMath::Sqrt(xr[0]*xr[0]+xr[1]*xr[1]);

      Double_t decayLength = sqrt( pow(xr[0] - primaryVtxPosition[0],2.0) 
				   +pow(xr[1] - primaryVtxPosition[1],2.0)
				   +pow(xr[2] - primaryVtxPosition[2],2.0)
				   );
  
	 
      Double_t  dcaV0ToPrimVertex= v0MIs->GetD(primaryVtxPosition[0],primaryVtxPosition[1]);////v0K0KF.GetDistanceFromVertexXY(tPrimaryVtxPosition);

      //-------------------- general cuts -------------------//
      //-- esd track cuts for daughters --//
      // if(fESDTrackCuts){
      //   if(!fESDTrackCuts->AcceptTrack(trackPosTest)) continue;
      //   if(!fESDTrackCuts->AcceptTrack(trackNegTest)) continue;
      // }
      
      // fHistPiPiMonitorCuts->Fill(5);
      // fHistPiPMonitorCuts[isSecd]->Fill(5);
      // fHistPiAPMonitorCuts[isSecd]->Fill(5);
      
      //-- eta cut --//
      if( fabs(posDaughterEta) > fEtaCutMCDaughtersVal || fabs(negDaughterEta) > fEtaCutMCDaughtersVal) continue;
      fHistPiPiMonitorCuts->Fill(5);
      fHistPiPMonitorCuts[isSecd]->Fill(5);
      fHistPiAPMonitorCuts[isSecd]->Fill(5);


      
      //-- radius xy min cut --//
      if(dim2V0Radius < fDecayRadXYMin) continue;
      //	    if(fabs(xr[1])<fDecayRadY) continue;
      fHistPiPiMonitorCuts->Fill(6);
      fHistPiPMonitorCuts[isSecd]->Fill(6);
      fHistPiAPMonitorCuts[isSecd]->Fill(6);

      //-- radius xy max cut --//
      if(dim2V0Radius > fDecayRadXYMax) continue;
      //	    if(fabs(xr[1])<fDecayRadY) continue;
      fHistPiPiMonitorCuts->Fill(7);
      fHistPiPMonitorCuts[isSecd]->Fill(7);
      fHistPiAPMonitorCuts[isSecd]->Fill(7);
      
      //-- decay length min ->ctau --//
      if(decayLength > fDecayLengthMax) continue;
      fHistPiPiMonitorCuts->Fill(8);
      fHistPiPMonitorCuts[isSecd]->Fill(8);
      fHistPiAPMonitorCuts[isSecd]->Fill(8);
	 
      //-- decay length min cut --//
      if(decayLength < fDecayLengthMin) continue;
      fHistPiPiMonitorCuts->Fill(9);
      fHistPiPMonitorCuts[isSecd]->Fill(9);
      fHistPiAPMonitorCuts[isSecd]->Fill(9);
   
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
      
      

      if(fMoreNclsThanRows && (crossedRowsTPCPos < nclsTPCPos || crossedRowsTPCNeg < nclsTPCNeg  )) continue;
      fHistPiPiMonitorCuts->Fill(10);
      fHistPiPMonitorCuts[isSecd]->Fill(10);
      fHistPiAPMonitorCuts[isSecd]->Fill(10);
      
      if(fMoreNclsThanFindable && (nclsTPCFindablePos < nclsTPCPos || nclsTPCFindableNeg < nclsTPCNeg  )) continue;
      fHistPiPiMonitorCuts->Fill(11);
      fHistPiPMonitorCuts[isSecd]->Fill(11);
      fHistPiAPMonitorCuts[isSecd]->Fill(11);      

      if(fMoreNclsThanFindableMax && ( nclsTPCPos < (nclsTPCFindablePos -60.0) || nclsTPCNeg < (nclsTPCFindableNeg - 60.0) )) continue;
      // if(chi2PerClusterITSNeg > fChi2PerClusterITS || chi2PerClusterITSPos > fChi2PerClusterITS ) continue;
      fHistPiPiMonitorCuts->Fill(12);
      fHistPiPMonitorCuts[isSecd]->Fill(12);
      fHistPiAPMonitorCuts[isSecd]->Fill(12);  
      
      if(ratio > fRatioMaxCRowsOverFindable || ratioPos > fRatioMaxCRowsOverFindable) continue; 

      if(ratioFoFi < fRatioFoundOverFindable || ratioFoFiPos < fRatioFoundOverFindable) continue; 
      fHistPiPiMonitorCuts->Fill(13);
      fHistPiPMonitorCuts[isSecd]->Fill(13);
      fHistPiAPMonitorCuts[isSecd]->Fill(13); 

      Bool_t cutOKITSNegNeg =kTRUE;
      Bool_t cutOKITSPosPos =kTRUE;

      Bool_t cutOKITSNegPos =kTRUE;
      Bool_t cutOKITSPosNeg =kTRUE;

      if(nclsITSNeg < fMinNCLSITSNeg ||  nclsITSNeg > fMaxNCLSITSNeg){
	if(!fSwitchCaseITSCls) continue;
	else cutOKITSNegNeg = kFALSE;
      }
  
      fHistPiPiMonitorCuts->Fill(14);
      fHistPiPMonitorCuts[isSecd]->Fill(14);
      fHistPiAPMonitorCuts[isSecd]->Fill(14); 

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

      fHistPiPiMonitorCuts->Fill(15);
      fHistPiPMonitorCuts[isSecd]->Fill(15);
      fHistPiAPMonitorCuts[isSecd]->Fill(15); 

    
      Double_t lengthTPCPos = trackPos->GetLengthInActiveZone(0,3,236, fESD->GetMagneticField(),0,0);//-5 ,0,0);
      Double_t lengthTPCNeg = trackNeg->GetLengthInActiveZone(0,3,236, fESD->GetMagneticField(),0,0);//-5 ,0,0);
   

      if(fCutMITrackLength && lengthTPCPos <=  fCutMITrackLengthLengthF * (130.0 - 5.0*fabs(trackPos->GetSigned1Pt()))) continue;
      fHistPiPiMonitorCuts->Fill(16);
      fHistPiPMonitorCuts[isSecd]->Fill(16);
      fHistPiAPMonitorCuts[isSecd]->Fill(16); 
      if(fCutMITrackLength && lengthTPCNeg <=  fCutMITrackLengthLengthF * (130.0 - 5.0*fabs(trackNeg->GetSigned1Pt()))) continue;
      fHistPiPiMonitorCuts->Fill(17);
      fHistPiPMonitorCuts[isSecd]->Fill(17);
      fHistPiAPMonitorCuts[isSecd]->Fill(17); 


      if(fCutMICrossedR && trackPos->GetTPCClusterInfo(3,1) <= fCutMICrossedRLengthF *(130.0 - 5.0*fabs(trackPos->GetSigned1Pt()))) continue;
      fHistPiPiMonitorCuts->Fill(18);
      fHistPiPMonitorCuts[isSecd]->Fill(18);
      fHistPiAPMonitorCuts[isSecd]->Fill(18); 
      if(fCutMICrossedR && trackNeg->GetTPCClusterInfo(3,1) <= fCutMICrossedRLengthF *(130.0 - 5.0*fabs(trackNeg->GetSigned1Pt()))) continue;
      fHistPiPiMonitorCuts->Fill(19);
      fHistPiPMonitorCuts[isSecd]->Fill(19);
      fHistPiAPMonitorCuts[isSecd]->Fill(19);

      if(fCutMITPCncls &&  nclsTPCPos <= 0.6*(130.0 - 5.0*fabs(trackPos->GetSigned1Pt()))) continue;
      fHistPiPiMonitorCuts->Fill(20);
      fHistPiPMonitorCuts[isSecd]->Fill(20);
      fHistPiAPMonitorCuts[isSecd]->Fill(20); 
      if(fCutMITPCncls &&   nclsTPCNeg <= 0.6*(130.0 - 5.0*fabs(trackNeg->GetSigned1Pt()))) continue;
      fHistPiPiMonitorCuts->Fill(21);
      fHistPiPMonitorCuts[isSecd]->Fill(21);
      fHistPiAPMonitorCuts[isSecd]->Fill(21); 
    
   
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

      Double_t posDaughterP = ppTrack.Mag();
      Double_t negDaughterP = pmTrack.Mag();

      Double_t tpcsigPos= trackPos->GetTPCsignal();
      Double_t tpcsigNeg= trackNeg->GetTPCsignal();
     
      /*
	Double_t tpcsigNPos= trackPos->GetTPCsignalN();
	Double_t tpcsigNNeg= trackNeg->GetTPCsignalN();
      */
     
      //AliExternalTrackParam *extTParPos = (AliExternalTrackParam*)trackPos->GetTPCInnerParam();
      //Double_t tpcMomPos = extTParPos->GetP();
      Double_t tpcMomPos =trackPos->GetInnerParam()->GetP();
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

      //-- momenta --//
      Double_t ptK0s = v0K0.Pt();
      Double_t ptLambda = v0Lambda.Pt();
      Double_t ptALambda = v0ALambda.Pt();
      
      Double_t pK0s = v0K0.P();
      Double_t pLambda = v0Lambda.P();
      Double_t pALambda = v0ALambda.P();
      
      //-- masses --//
      Double_t massK0s = v0K0.M();
      Double_t massLambda = v0Lambda.M();
      Double_t massALambda = v0ALambda.M();

      TLorentzVector e1(ppTrack,0.51099);
      TLorentzVector e2(pmTrack,0.51099);
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
      Double_t phiK0 =  v0K0.Phi();
      Double_t phiL =  v0Lambda.Phi();
      Double_t phiAL =  v0ALambda.Phi();
      */
            
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

      //------------------ cut flags for V0 specific cuts --------------//
      Bool_t cutOKK0s = kTRUE;
      Bool_t cutOKLambda = kTRUE;
      Bool_t cutOKALambda = kTRUE;

      //-------------------------- K0 cuts -----------------------------//

      if(dcaV0ToPrimVertex > fDCAToVertexK0)  cutOKK0s = kFALSE;//continue;
      else fHistPiPiMonitorCuts->Fill(22);
      
      if(fabs(xr[2])> fDCAZ) cutOKK0s = kFALSE; //like decay radius z component
      else fHistPiPiMonitorCuts->Fill(23);
      
      Double_t ctK0 = 0.0,ctTK0 = 0.0;
      if(fabs(pK0s)>0.0)  ctK0 = decayLength*0.497614/pK0s;
      if(fabs(ptK0s)>0.0)  ctTK0 = dim2V0Radius*0.497614/ptK0s;
      if(ctK0 > fCtauK0s &&  fabs(ptK0s) <fCtauPtCutK0) cutOKK0s = kFALSE;
      else  fHistPiPiMonitorCuts->Fill(24);
      
      if((cosOPAng < fCosPointAngK && fabs(ptK0s) < fCPAPtCutK0)|| cosOPAng<0.99)
	cutOKK0s = kFALSE;
      else  fHistPiPiMonitorCuts->Fill(25);

      if(dcaDaughters > fDCADaughtersK0 )cutOKK0s = kFALSE;
      else  fHistPiPiMonitorCuts->Fill(26);
	 
      if(dcaNegToVertex < fDCADaughtersToVtxSmall || dcaPosToVertex < fDCADaughtersToVtxSmall)  cutOKK0s = kFALSE;
      else  fHistPiPiMonitorCuts->Fill(27);

      if(fRapCutV0 && fabs(rapK0s) > fRap) cutOKK0s = kFALSE;
      else  fHistPiPiMonitorCuts->Fill(28);  
    
      // if(chi2K0C > fChiCutKf) cutOKK0s = kFALSE;
      if(opAng < fOpengAngleDaughters && fabs(ptK0s) < fOpAngPtCut )  cutOKK0s = kFALSE;
      else fHistPiPiMonitorCuts->Fill(29);
    
      Bool_t ptbinokK0s=kFALSE;
      if( ptK0s < fQtCutPt &&  ptK0s > fQtCutPtLow ) ptbinokK0s=kTRUE;
    
      Double_t qtval = fArmQtSlope*fabs(alfa);
   
      if(fArmCutK0 && ptbinokK0s && qt < qtval) cutOKK0s = kFALSE;
      if(fArmCutK0 && ptbinokK0s && qt < fQtCut) cutOKK0s = kFALSE;
      else  fHistPiPiMonitorCuts->Fill(30);
    
      if( ptK0s > fPtTPCCut){
	if(fESDTrackCuts){
	  if(!fESDTrackCuts->AcceptTrack(trackPosTest) || !fESDTrackCuts->AcceptTrack(trackNegTest)) cutOKK0s = kFALSE;
	  else  fHistPiPiMonitorCuts->Fill(31); 
	}
      }
      else{
	if(fESDTrackCutsLowPt){
	  if(!fESDTrackCutsLowPt->AcceptTrack(trackPosTest) || !fESDTrackCutsLowPt->AcceptTrack(trackNegTest))  cutOKK0s = kFALSE; 
	}
      }

      //-------------------------- Lambda cuts -------------------------//
      if(dcaV0ToPrimVertex > fDCAToVertexL) cutOKLambda = kFALSE;//continue;
      else  fHistPiPMonitorCuts[isSecd]->Fill(22);

      if(fabs(xr[2])>fDCAZ) cutOKLambda = kFALSE; //like decay radius z component
      else  fHistPiPMonitorCuts[isSecd]->Fill(23);
         
      Double_t ctL = 0.0,ctTL=0.0;
      if(fabs(pLambda)>0.0)  ctL  = decayLength*1.115683/fabs(pLambda);
      if(fabs(ptLambda)>0.0) ctTL = dim2V0Radius*1.115683/fabs(ptLambda);
	 
      if(ctL > fCtauL && fabs(ptLambda) <fCtauPtCutL)  cutOKLambda = kFALSE;
      else  fHistPiPMonitorCuts[isSecd]->Fill(24);
      
      if((cosOPAng<fCosPointAngL && fabs(ptLambda) < fCPAPtCutL)|| cosOPAng<0.99)
	cutOKLambda = kFALSE;
      else fHistPiPMonitorCuts[isSecd]->Fill(25);

      if(dcaDaughters > fDCADaughtersL )cutOKLambda = kFALSE;
      else  fHistPiPMonitorCuts[isSecd]->Fill(26);
 
      if( dcaNegToVertex < fDCADaughtersToVtxSmall || dcaPosToVertex < fDCADaughtersToVtxLarge)  cutOKLambda = kFALSE;
      else  fHistPiPMonitorCuts[isSecd]->Fill(27);

      if(fRapCutV0 && fabs(rapL) > fRap) cutOKLambda = kFALSE;
      else  fHistPiPMonitorCuts[isSecd]->Fill(28);
       
   
      /*	 
		 if(chi2LambdaC > fChiCutKf) cutOKLambda = kFALSE;
		 else  fHistPiPMonitorCuts[isSecd]->Fill(20);
      */

      if(opAng < fOpengAngleDaughters && fabs(ptLambda) < fOpAngPtCut )  cutOKLambda = kFALSE;
      else  fHistPiPMonitorCuts[isSecd]->Fill(29);
    

      if(alfa<fAlfaCut  || (fArmCutL && qt > fQtCut)) cutOKLambda = kFALSE;
      else  fHistPiPMonitorCuts[isSecd]->Fill(30);

      if(ptLambda > fPtTPCCut){
	if(fESDTrackCuts && fESDTrackCutsCharged){
	  if(!fESDTrackCutsCharged->AcceptTrack(trackPosTest) || !fESDTrackCuts->AcceptTrack(trackNegTest)) cutOKLambda = kFALSE;
	  else  fHistPiPMonitorCuts[isSecd]->Fill(31); 
	}
      }
      else{
	if(fESDTrackCutsLowPt){
	  if(!fESDTrackCutsLowPt->AcceptTrack(trackPosTest) || !fESDTrackCutsLowPt->AcceptTrack(trackNegTest))  cutOKLambda = kFALSE; 
	}
      }

      //--------------------------- ALambda cuts --------------------------//

      if(dcaV0ToPrimVertex > fDCAToVertexL) cutOKALambda = kFALSE;//continue;
      else fHistPiAPMonitorCuts[isSecd]->Fill(22);
 
      if(fabs(xr[2])> fDCAZ) cutOKALambda = kFALSE;//continue;//like decay radius z component
      else fHistPiAPMonitorCuts[isSecd]->Fill(23);

      Double_t ctAL = 0.0,ctTAL=0.0;
      if(fabs(pALambda)>0.0)  ctAL  = decayLength*1.115683/fabs(pALambda);
      if(fabs(ptALambda)>0.0) ctTAL = dim2V0Radius*1.115683/fabs(ptALambda);
      if(ctAL > fCtauL &&  fabs(ptALambda) <fCtauPtCutL)  cutOKALambda = kFALSE;
      else  fHistPiAPMonitorCuts[isSecd]->Fill(24);

      if((cosOPAng<fCosPointAngL && fabs(ptALambda) < fCPAPtCutL)|| cosOPAng<0.99)  cutOKALambda = kFALSE;
      else fHistPiAPMonitorCuts[isSecd]->Fill(25);
      
      if(dcaDaughters > fDCADaughtersAL )cutOKALambda = kFALSE;
      else  fHistPiAPMonitorCuts[isSecd]->Fill(26);
	 
      if( dcaPosToVertex < fDCADaughtersToVtxSmall || dcaNegToVertex < fDCADaughtersToVtxLarge)  cutOKALambda = kFALSE;
      else fHistPiAPMonitorCuts[isSecd]->Fill(27);
	 
      if(fRapCutV0 && fabs(rapAL) > fRap) cutOKALambda = kFALSE;
      else fHistPiAPMonitorCuts[isSecd]->Fill(28);

      /*
	if(chi2ALambdaC > fChiCutKf) cutOKALambda = kFALSE;
	else  fHistPiAPMonitorCuts[isSecd]->Fill(20);
      */
     
      if(opAng < fOpengAngleDaughters && fabs(ptALambda) < fOpAngPtCut )  cutOKALambda = kFALSE;
      else  fHistPiAPMonitorCuts[isSecd]->Fill(29);
    
      
      if((fArmCutL && qt>qtval) || alfa > -1.0*fAlfaCut) cutOKALambda = kFALSE;
      else  fHistPiAPMonitorCuts[isSecd]->Fill(30);

      if(ptALambda > fPtTPCCut){
	if(fESDTrackCuts && fESDTrackCutsCharged){
	  if(!fESDTrackCuts->AcceptTrack(trackPosTest) || !fESDTrackCutsCharged->AcceptTrack(trackNegTest)) cutOKALambda = kFALSE;
	  else  fHistPiAPMonitorCuts[isSecd]->Fill(31); 
	}
      }
      else{
	if(fESDTrackCutsLowPt){
	  if(!fESDTrackCutsLowPt->AcceptTrack(trackPosTest))  cutOKALambda = kFALSE; 
	  if(!fESDTrackCutsLowPt->AcceptTrack(trackNegTest))  cutOKALambda = kFALSE; 
	}
      }
  
      //-------------------- V0 ana -------------------------//

      //-- cut flags for furhter histos--//
      Bool_t k0sOK=kFALSE;
      Bool_t lambdaOK=kFALSE;
      Bool_t alambdaOK=kFALSE;

      //--  Check for K0 --//
      Bool_t exMass = kFALSE;
      if(fabs(1.115 - massLambda)  < fExcludeLambdaFromK0s){
	cutOKK0s = kFALSE;
	exMass = kTRUE;
      }
      if(fabs(1.115 - massALambda) < fExcludeLambdaFromK0s){
	cutOKK0s = kFALSE;
	exMass = kTRUE;
      }
   
      if(fabs(massPhoton) < fExcludeLambdaFromK0s) {
	cutOKK0s = kFALSE;
	exMass = kTRUE;
      }
   
      if(ptK0s >fMinPt){
	if( cutOKK0s  && fillK0sMC ){
	  fHistDedxPiPlus[isSecd]->Fill(tpcMomPos,tpcsigPos);
	  fHistDedxPiMinus[isSecd]->Fill(tpcMomNeg,tpcsigNeg);
	  fHistPiPiMonitorCuts->Fill(32);
	  if(pipidEdx){
	    fHistPiPiMonitorCuts->Fill(33);
	    k0sOK = kTRUE;		    
	    if(!exMass && massK0s > fK0sLowMassCut && massK0s < fK0sHighMassCut ){
	      if(!fMCMode){
		ptV0MC = ptK0s;
		declengthV0MC = dim2V0Radius;
	      }
	      fHistPiPiMonitorCuts->Fill(34);
	      fHistPiPiMass->Fill(massK0s);
	      fHistPiPiMassVSPt->Fill(massK0s,ptK0s);
	      fHistPiPiMassVSPtMCTruth->Fill(massK0s,ptV0MC);
	      fHistPiPiMassVSY->Fill(massK0s,rapK0s);
	      fHistPiPiPtVSY->Fill(rapK0s,ptK0s);
	      fHistPiPiDecayLengthVsMass->Fill(massK0s,dim2V0Radius);//decayLength);
	      if(massK0s > 0.46 && massK0s < 0.53)  fHistPiPiDecayLengthVsPt->Fill(ptV0MC,dim2V0Radius);//decayLength
	      // fHistPiPiPtDaughters->Fill(posDaughterPt,negDaughterPt);
	      if(!fSetPtDepHist){
		fHistPiPiRadiusXY->Fill(massK0s,opAng);//posDaughterPhi);//opAng);
		//	    fHistPiPiPhiPosVsPtPosVsMass->Fill(massK0s,dim2V0Radius,ptV0MC);//,ctK0);//posDaughterPhi);//xxx
		fHistPiPiCosPointAng->Fill(massK0s,cosOPAng);
		fHistPiPiDecayLengthVsCtau->Fill(massK0s,ctTK0);
		fHistPiPiTrackLengthPosVsMass->Fill(massK0s,lengthTPCPos);
		fHistPiPiTrackLengthNegVsMass->Fill(massK0s,lengthTPCNeg);
		fHistPiPiDCAZPos->Fill(massK0s,dcaToVertexZPos);
		fHistPiPiDCAZNeg->Fill(massK0s,dcaToVertexZNeg);
		fHistPiPiDCADaughters->Fill(massK0s,dcaDaughters);
		fHistPiPiDCADaughterPosToPrimVtxVSMass->Fill(massK0s,dcaPosToVertex);
		fHistPiPiDCAVSMass->Fill(massK0s,dcaV0ToPrimVertex);
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
	      }


	      if(fMCMode)  fHistPiPiDecayLengthResolution->Fill(declengthV0MC,dim2V0Radius);

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
      //--  Check for Lambda --//
      Bool_t  exMassL =kFALSE;
      if(fabs(0.497 - massK0s) < fExcludeK0sFromLambda){
	cutOKLambda = kFALSE;
	exMassL = kTRUE;
      }
      if(fabs(massPhoton) < fExcludeK0sFromLambda) {
	cutOKLambda = kFALSE;
	exMassL = kTRUE;
      }
    
      if(ptLambda > fMinPt){
	if(cutOKLambda && fillLambdaMC){
	  fHistDedxProt[isSecd]->Fill(tpcMomPos,tpcsigPos);
	  fHistDedxPiMinus[isSecd]->Fill(tpcMomNeg,tpcsigNeg);
	  fHistPiPMonitorCuts[isSecd]->Fill(32);
	  if(pipdEdx){
	    fHistPiPMonitorCuts[isSecd]->Fill(33);
	    lambdaOK = kTRUE;
	    if(!exMassL && massLambda > fLLowMassCut && massLambda < fLHighMassCut){// 1.05 && massLambda < 1.25 ){
	      if(!fMCMode) {
		ptV0MC = ptLambda;
		declengthV0MC = dim2V0Radius;
	      }
	      fHistPiPMonitorCuts[isSecd]->Fill(34);
	      fHistPiPMass[isSecd]->Fill(massLambda);
	      fHistPiPMassVSPt[isSecd]->Fill(massLambda,ptLambda);
	      fHistPiPMassVSPtMCTruth[isSecd]->Fill(massLambda,ptV0MC);
	      fHistPiPMassVSY[isSecd]->Fill(massLambda,rapL);
	      fHistPiPPtVSY[isSecd]->Fill(rapL,ptLambda);
	      //fHistPiPDecayLengthVsPt[isSecd]->Fill(ptLambda,ctL);
	      if( massLambda > 1.108 && massLambda < 1.123 ) fHistPiPDecayLengthVsPt[isSecd]->Fill(ptV0MC,dim2V0Radius);//decayLength);
	      fHistPiPDecayLengthVsMass[isSecd]->Fill(massLambda,dim2V0Radius);//decayLength);
	      //  fHistPiPPhiPosVsPtPosVsMass->Fill(massLambda,dim2V0Radius,ptV0MC);//

	      if(!fSetPtDepHist){
		fHistPiPRadiusXY[isSecd]->Fill(massLambda,opAng);//,posDaughterPhi);//opAng);
		fHistPiPCosPointAng[isSecd]->Fill(massLambda,cosOPAng);
		fHistPiPTrackLengthPosVsMass[isSecd]->Fill(massLambda,lengthTPCPos);
		fHistPiPTrackLengthNegVsMass[isSecd]->Fill(massLambda,lengthTPCNeg);
		//fHistPiPDCAZPos[isSecd]->Fill(massLambda,dcaToVertexZPos);
		//fHistPiPDCAZNeg[isSecd]->Fill(massLambda,dcaToVertexZNeg);
		fHistPiPDCADaughters[isSecd]->Fill(massLambda,dcaDaughters);
		fHistPiPDCAVSMass[isSecd]->Fill(massLambda,dcaV0ToPrimVertex);
		fHistPiPDCADaughterPosToPrimVtxVSMass[isSecd]->Fill(massLambda,dcaPosToVertex);
		fHistPiPDCADaughterNegToPrimVtxVSMass[isSecd]->Fill(massLambda,dcaNegToVertex);
		fHistPiPDecayLengthVsCtau[isSecd]->Fill(massLambda,ctTL);
	      }
	      else{
		fHistPiPRadiusXY[isSecd]->Fill(ptV0MC,opAng);
		fHistPiPCosPointAng[isSecd]->Fill(ptV0MC,cosOPAng);
		fHistPiPTrackLengthPosVsMass[isSecd]->Fill(ptV0MC,lengthTPCPos);
		fHistPiPTrackLengthNegVsMass[isSecd]->Fill(ptV0MC,lengthTPCNeg);
		//fHistPiPDCAZPos[isSecd]->Fill(ptV0MC,dcaToVertexZPos);
		//fHistPiPDCAZNeg[isSecd]->Fill(ptV0MC,dcaToVertexZNeg);
		fHistPiPDCADaughters[isSecd]->Fill(ptV0MC,dcaDaughters);
		fHistPiPDCAVSMass[isSecd]->Fill(ptV0MC,dcaV0ToPrimVertex);
		fHistPiPDCADaughterPosToPrimVtxVSMass[isSecd]->Fill(ptV0MC,dcaPosToVertex);
		fHistPiPDCADaughterNegToPrimVtxVSMass[isSecd]->Fill(ptV0MC,dcaNegToVertex);
		fHistPiPDecayLengthVsCtau[isSecd]->Fill(ptV0MC,ctTL);
	      }

	      if(fMCMode)  fHistPiPDecayLengthResolution[isSecd]->Fill(declengthV0MC,dim2V0Radius);
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
      if(fabs(massPhoton) < fExcludeK0sFromLambda) {
	exMassAL = kTRUE;
      }

      if(ptALambda > fMinPt){
	if(cutOKALambda && fillALambdaMC){
	  fHistDedxAProt[isSecd]->Fill(tpcMomNeg,tpcsigNeg);
	  fHistDedxPiPlus[isSecd]->Fill(tpcMomPos,tpcsigPos);
	  fHistPiAPMonitorCuts[isSecd]->Fill(32);
	  if(piapdEdx){
	    fHistPiAPMonitorCuts[isSecd]->Fill(33);
	    alambdaOK = kTRUE;
	    if( !exMassAL && massALambda > fLLowMassCut && massALambda < fLHighMassCut){//1.05 && massALambda < 1.25  ){
	      if(!fMCMode) {
		ptV0MC = ptALambda;
		declengthV0MC = dim2V0Radius;
	      }
	      fHistPiAPMonitorCuts[isSecd]->Fill(34);
	      fHistPiAPMass[isSecd]->Fill(massALambda);
	      fHistPiAPMassVSPt[isSecd]->Fill(massALambda,ptALambda);
	      fHistPiAPMassVSPtMCTruth[isSecd]->Fill(massALambda,ptV0MC);
	      fHistPiAPMassVSY[isSecd]->Fill(massALambda,rapAL);
	      fHistPiAPPtVSY[isSecd]->Fill(rapAL,ptALambda);
	      //  fHistPiAPPtDaughters[isSecd]->Fill(posDaughterPt,negDaughterPt);
	      if( massALambda>1.108 && massALambda<1.123 )  fHistPiAPDecayLengthVsPt[isSecd]->Fill(ptV0MC,dim2V0Radius);//decayLength);
	      fHistPiAPDecayLengthVsMass[isSecd]->Fill(massALambda,dim2V0Radius);//decayLength);
	    
	      if(!fSetPtDepHist){
		fHistPiAPRadiusXY[isSecd]->Fill(massALambda,posDaughterPhi);//,opAng);
		fHistPiAPCosPointAng[isSecd]->Fill(massALambda,cosOPAng);	    
		fHistPiAPTrackLengthPosVsMass[isSecd]->Fill(massALambda,lengthTPCPos);
		fHistPiAPTrackLengthNegVsMass[isSecd]->Fill(massALambda,lengthTPCNeg);
		//fHistPiAPDCAZPos[isSecd]->Fill(massALambda,dcaToVertexZPos);
		//fHistPiAPDCAZNeg[isSecd]->Fill(massALambda,dcaToVertexZNeg);
		fHistPiAPDCADaughters[isSecd]->Fill(massALambda,dcaDaughters);
		fHistPiAPDCAVSMass[isSecd]->Fill(massALambda,dcaV0ToPrimVertex);
		fHistPiAPDCADaughterPosToPrimVtxVSMass[isSecd]->Fill(massALambda,dcaPosToVertex);
		fHistPiAPDCADaughterNegToPrimVtxVSMass[isSecd]->Fill(massALambda,dcaNegToVertex);
		fHistPiAPDecayLengthVsCtau[isSecd]->Fill(massALambda,ctTAL);
	   
	      }
	      else{
		fHistPiAPRadiusXY[isSecd]->Fill(ptV0MC,posDaughterPhi);
		fHistPiAPCosPointAng[isSecd]->Fill(ptV0MC,cosOPAng);
		fHistPiAPTrackLengthPosVsMass[isSecd]->Fill(ptV0MC,lengthTPCPos);
		fHistPiAPTrackLengthNegVsMass[isSecd]->Fill(ptV0MC,lengthTPCNeg);
		//fHistPiAPDCAZPos[isSecd]->Fill(ptV0MC,dcaToVertexZPos);
		//fHistPiAPDCAZNeg[isSecd]->Fill(ptV0MC,dcaToVertexZNeg);
		fHistPiAPDCADaughters[isSecd]->Fill(ptV0MC,dcaDaughters);
		fHistPiAPDCAVSMass[isSecd]->Fill(ptV0MC,dcaV0ToPrimVertex);
		fHistPiAPDCADaughterPosToPrimVtxVSMass[isSecd]->Fill(ptV0MC,dcaPosToVertex);
		fHistPiAPDCADaughterNegToPrimVtxVSMass[isSecd]->Fill(ptV0MC,dcaNegToVertex);
		fHistPiAPDecayLengthVsCtau[isSecd]->Fill(ptV0MC,ctTAL);
	      }

	      if(fMCMode)  fHistPiAPDecayLengthResolution[isSecd]->Fill(declengthV0MC,dim2V0Radius);
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
	
      //-- fill detector histos general --//
      if((lambdaOK && !exMassL) || (alambdaOK && !exMassAL) || (k0sOK && !exMass)){
	fHistPiPiEtaDReco[1]->Fill(posDaughterEta,eta00);
	fHistPiPEtaDReco[1]->Fill(negDaughterEta,eta01);
      }
   
    

      /*
      //-- AliKFParticle --//
      if (negPiKF) delete negPiKF; negPiKF=NULL;
      if (posPiKF) delete posPiKF; posPiKF=NULL;
      if (posPKF) delete posPKF; posPKF=NULL;
      if (negAPKF) delete negAPKF; negAPKF=NULL;
      */
	
	
    }//---- end V0 reco loop----//
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

    // cout<<pos->GetLabel()<<"  "<<neg->GetLabel()<<"  id "<<id0<<"  "<<id1<<endl;
      
    if (labelN==labelP)  return kFALSE;
      
    if(fMCTruthMode){
      if ((labelP!=id0) && (labelP!=id1))  return kFALSE;
      if ((labelN!=id0) && (labelN!=id1))  return kFALSE;
    }
    // cout<<"yes"<<endl;
    return kTRUE;
  }

