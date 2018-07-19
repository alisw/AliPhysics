/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, proviyaded that the above copyright notice appears in all *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purapose. It is         *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

///////////////////////////////////////////////////////////////////////////
//                                                                       //
//                                                                       //
//          Analysis for event-by-event particle ratio studies           //
//                                                                       //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include "Riostream.h"
#include "TChain.h"
#include "TTree.h"
#include "TCutG.h"
#include "TH1F.h"
#include "THn.h"
#include "THnSparse.h"
#include "TList.h"
#include "TMath.h"
#include "TMatrixF.h"
#include "TVectorF.h"
#include "TCanvas.h"
#include "TObjArray.h"
#include "TF1.h"
#include "TFile.h"
#include "TDatabasePDG.h"
#include "AliExternalTrackParam.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliTPCdEdxInfo.h"
#include "AliKFVertex.h"
#include "AliKFParticle.h"
#include "AliCollisionGeometry.h"
#include "AliGenEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliGenCocktailEventHeader.h"
#include "AliRun.h"
#include "AliHeader.h"
#include "AliPID.h"
#include "AliESDtrackCuts.h"
#include "AliESDv0Cuts.h"
#include "AliESDv0.h"
#include "AliESDVertex.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrack.h"
#include "AliESDpid.h"
#include "AliCentrality.h"
#include "AliESDUtils.h"
#include "AliMultiplicity.h"
#include "AliPIDResponse.h"
#include "AliTPCPIDResponse.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliLog.h"
#include "AliKFParticle.h"
#include "AliAnalysisTaskFilteredTree.h"
#include "AliAnalysisTaskEbyeIterPID.h"
#include "AliMultSelection.h"
#include "AliRunLoader.h"
#include <fstream>
#include <iostream>
#include <iomanip>
using namespace std;
using std::cout;
using std::setw;



ClassImp(AliAnalysisTaskEbyeIterPID)

#define USE_STREAMER 1

// -----------------------------------------------------------------------
//                            Constructor and Destructor
// -----------------------------------------------------------------------
//________________________________________________________________________
AliAnalysisTaskEbyeIterPID::AliAnalysisTaskEbyeIterPID() 
  : AliAnalysisTaskSE("TaskEbyeRatios"), fPIDResponse(0),fESD(0), fListHist(0), fESDtrackCuts(0),
fESDtrackCutsV0(0),
fESDtrackCutsCleanSamp(0),
fPIDCombined(0x0),
fTPCdEdxInfo(0x0),
fMCStack(0x0),
fTree(0x0),
fIdenTree(0x0),
fIdenTreeMC(0x0),
fArmPodTree(0x0),
fTreeSRedirector(0x0),
fTreeMCrec(0x0),
fTreeMCgen(0x0),
fTreeDnchDeta(0x0),
fTreeMC(0x0),
fTreedEdxCheck(0x0),
fTreeBayes(0x0),
fTreeCuts(0x0),
fTreeMCFullAcc(0x0),
fTreeResonance(0x0),
fTreeMCgenMoms(0x0),
fhEta(0),
fhCent(0),
fhPtot(0),
fhndEdx(0),
fhnExpected(0),
fhnCleanEl(0),
fhnCleanKa(0),
fhnCleanDe(0),
fTrackCutBits(0),
fEtaDown(0),
fEtaUp(0),
fnEtaBins(0),
fPercentageOfEvents(0),
fMCtrue(kFALSE),
fTIdentity(kFALSE),
fWeakAndMaterial(kFALSE),
fEffMatrix(kFALSE),
fdEdxCheck(kFALSE),
fCleanSamplesOnly(kFALSE),
fTightCuts(kFALSE),
fIncludeITS(kTRUE),
fFillBayes(kFALSE),
fFillCuts(kFALSE),
fFillDeDxTree(kTRUE),
fFillArmPodTree(kTRUE),
fRunFastSimulation(kFALSE),
fRunFastHighMomentCal(kFALSE),
fFillDnchDeta(kFALSE),
fIncludeTOF(kFALSE),
fUseThnSparse(kFALSE),
fUseCouts(kFALSE),
fnMomBins(0),
fMomDown(0),                
fMomUp(0), 
fdEdxBinWidth(0),
fdEdxUp(0),                 
fdEdxDown (0),                
fdEdxCleanUp (0),  
fArmPodTPCSignal(0),
fArmPodptot(0),
fArmPodEta(0),
fArmPodCentrality(0),
fQt(0),
fAlfa(0),
fPiNSigmasTOF(0),
fPrNSigmasTOF(0),
fdEdxEl(0),
fdEdxKa(0),
fdEdxPi(0),
fdEdxPr(0),
fdEdxDe(0),
fSigmaEl(0),
fSigmaKa(0),
fSigmaPi(0),
fSigmaPr(0), 
fSigmaDe(0),
fNSigmasElTPC(0),
fNSigmasPiTPC(0),
fNSigmasKaTPC(0),
fNSigmasPrTPC(0),  
fNSigmasDeTPC(0),  
fTPCSignalMC(0),
fptotMC(0),
fptotMCtruth(0),
fpTMC(0),
fEtaMC(0),
fCentralityMC(0),
fSignMC(0),
fPxMC(0),
fPyMC(0),
fPzMC(0),
fElMC(0),
fPiMC(0),
fKaMC(0),
fPrMC(0),
fDeMC(0),
fMuMC(0),
fLaMC(0),
fptotMCgen(0),
fpTMCgen(0),
fEtaMCgen(0),
fCentralityMCgen(0),
fSignMCgen(0),
fMCImpactParameter(0),
fElMCgen(0),
fPiMCgen(0),
fKaMCgen(0),
fPrMCgen(0),
fDeMCgen(0),
fMuMCgen(0), 
fLaMCgen(0),
fPx(0),
fPy(0),
fPz(0),
fptot(0),
fpT(0),
fY(0),
fMultiplicity(0),
fMultiplicityMC(0),
fCentrality(0),
fvZ(0),
fEventGID(0),
fEventGIDMC(0),
fEventCountInFile(0),
fEvent(0),
fEventMC(0),
fEventMCgen(0),
fTPCSignal(0),
myDeDx(0),
signNew(0),
myDeDxMC(0),
signNewMC(0),        
fEta(0),
fNContributors(0),
fTheta(0),
fPhi(0),
fSign(0),
fTPCShared(0),
fNcl(0),
fnResBins(0),
fnEtaWinBinsMC(-100),
fnMomBinsMC(-100),
fnCentBinsMC(-100), 
fnResModeMC(2),
fnCentbinsData(10),
fMissingCl(0.),
fIsITSpixel01(0),          
fnITSclusters(0),
fPrimRestriction(0),
fTPCvZ(0),
fCleanPionsFromK0(0),
fCleanPion0FromK0(0),
fCleanPion1FromK0(0),
fCleanPion0FromLambda(0),
fCleanPion1FromLambda(0),
fCleanProton0FromLambda(0),
fCleanProton1FromLambda(0),
fHasTrack0FirstITSlayer(0),
fHasTrack1FirstITSlayer(0),
fHasV0FirstITSlayer(0),
fPionCutG(0),
fAntiProtonCutG(0),
fProtonCutG(0),
fSystCentEstimatetor(0),
fSystCrossedRows(0),
fSystDCAxy(0),
fSystChi2(0),
fSystVz(0),
fetaDownArr(0),           
fetaUpArr(0),
fcentDownArr(0),
fcentUpArr(0),
fpDownArr(0),
fpUpArr(0),
fxCentBins(0),
fResonances(0),
fHistPosEffMatrixRec(0),
fHistNegEffMatrixRec(0),
fHistPosEffMatrixGen(0),
fHistNegEffMatrixGen(0),
fHistEmptyEvent(0),
fHistCentrality(0),
fHistVertex(0),
fHistdEdxTPC(0),
fHistArmPod(0)
{
  // default Constructor
  /* fast compilation test
  gSystem->AddIncludePath("-I$ALICE_ROOT/include"); 
  gSystem->AddIncludePath("-I$ALICE_PHYSICS/include");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  .L /u/marsland/PHD/macros/marsland_EbyeRatios/AliAnalysisTaskEbyeIterPID.cxx++
  .L /hera/alice/marsland/train/trunk/marsland_EbyeRatios/AliAnalysisTaskEbyeIterPID.cxx++
  */
}

//________________________________________________________________________
AliAnalysisTaskEbyeIterPID::AliAnalysisTaskEbyeIterPID(const char *name) 
  : AliAnalysisTaskSE(name), fPIDResponse(0), fESD(0), fListHist(0), fESDtrackCuts(0),
fESDtrackCutsV0(0),
fESDtrackCutsCleanSamp(0),
fPIDCombined(0x0),
fTPCdEdxInfo(0x0),
fMCStack(0x0),
fTree(0),
fIdenTree(0x0),
fIdenTreeMC(0x0),
fArmPodTree(0x0),
fTreeSRedirector(0x0),
fTreeMCrec(0x0),
fTreeMCgen(0x0),
fTreeDnchDeta(0x0),
fTreeMC(0x0),
fTreedEdxCheck(0x0),
fTreeBayes(0x0),
fTreeCuts(0x0),
fTreeMCFullAcc(0x0),
fTreeResonance(0x0),
fTreeMCgenMoms(0x0),
fhEta(0),
fhCent(0),
fhPtot(0),
fhndEdx(0),
fhnExpected(0),
fhnCleanEl(0),
fhnCleanKa(0),
fhnCleanDe(0),
fTrackCutBits(0),
fEtaDown(0),
fEtaUp(0),
fnEtaBins(0),
fPercentageOfEvents(0),
fMCtrue(kFALSE),
fTIdentity(kFALSE),
fWeakAndMaterial(kFALSE),
fEffMatrix(kFALSE),
fdEdxCheck(kFALSE),
fCleanSamplesOnly(kFALSE),
fTightCuts(kFALSE),
fIncludeITS(kTRUE),
fFillBayes(kFALSE),
fFillCuts(kFALSE),
fFillDeDxTree(kTRUE),
fFillArmPodTree(kTRUE),
fRunFastSimulation(kFALSE),
fRunFastHighMomentCal(kFALSE),
fFillDnchDeta(kFALSE),
fIncludeTOF(kFALSE),
fUseThnSparse(kFALSE),
fUseCouts(kFALSE),
fnMomBins(0),
fMomDown(0),                
fMomUp(0),                  
fdEdxBinWidth(0),
fdEdxUp(0),                 
fdEdxDown(0),                
fdEdxCleanUp(0),             
fArmPodTPCSignal(0),
fArmPodptot(0),
fArmPodEta(0),
fArmPodCentrality(0),
fQt(0),
fAlfa(0),
fPiNSigmasTOF(0),
fPrNSigmasTOF(0),
fdEdxEl(0),
fdEdxKa(0),
fdEdxPi(0),
fdEdxPr(0),
fdEdxDe(0),
fSigmaEl(0),
fSigmaKa(0),
fSigmaPi(0),
fSigmaPr(0), 
fSigmaDe(0),
fNSigmasElTPC(0),
fNSigmasPiTPC(0),
fNSigmasKaTPC(0),
fNSigmasPrTPC(0),   
fNSigmasDeTPC(0),  
fTPCSignalMC(0),
fptotMC(0),
fptotMCtruth(0),
fpTMC(0),
fEtaMC(0),
fCentralityMC(0),
fSignMC(0), 
fPxMC(0),
fPyMC(0),
fPzMC(0),
fElMC(0),
fPiMC(0),
fKaMC(0),
fPrMC(0),
fDeMC(0),
fMuMC(0), 
fLaMC(0),
fptotMCgen(0),
fpTMCgen(0),
fEtaMCgen(0),
fCentralityMCgen(0),
fSignMCgen(0),
fMCImpactParameter(0),
fElMCgen(0),
fPiMCgen(0),
fKaMCgen(0),
fPrMCgen(0),
fDeMCgen(0),
fMuMCgen(0),  
fLaMCgen(0),
fPx(0),
fPy(0),
fPz(0),
fptot(0),
fpT(0),
fY(0),
fMultiplicity(0),
fMultiplicityMC(0),
fCentrality(0),
fvZ(0),
fEventGID(0),
fEventGIDMC(0),
fEventCountInFile(0),
fEvent(0),
fEventMC(0),
fEventMCgen(0),
fTPCSignal(0),
myDeDx(0),
signNew(0),
myDeDxMC(0),
signNewMC(0),
fEta(0),
fNContributors(0),
fTheta(0),
fPhi(0),
fSign(0),
fTPCShared(0),
fNcl(0),
fnResBins(0),
fnEtaWinBinsMC(-100),
fnMomBinsMC(-100),
fnCentBinsMC(-100), 
fnResModeMC(2),
fnCentbinsData(10),
fMissingCl(0.),
fIsITSpixel01(0),          
fnITSclusters(0),
fPrimRestriction(0),
fTPCvZ(0),
fCleanPionsFromK0(0),
fCleanPion0FromK0(0),
fCleanPion1FromK0(0),
fCleanPion0FromLambda(0),
fCleanPion1FromLambda(0),
fCleanProton0FromLambda(0),
fCleanProton1FromLambda(0),
fHasTrack0FirstITSlayer(0),
fHasTrack1FirstITSlayer(0),
fHasV0FirstITSlayer(0),
fPionCutG(0),
fAntiProtonCutG(0),
fProtonCutG(0),
fSystCentEstimatetor(0),
fSystCrossedRows(0),
fSystDCAxy(0),
fSystChi2(0),
fSystVz(0),
fetaDownArr(0),           
fetaUpArr(0),
fcentDownArr(0),
fcentUpArr(0),
fpDownArr(0),
fpUpArr(0),
fxCentBins(0),
fResonances(0),
fHistPosEffMatrixRec(0),
fHistNegEffMatrixRec(0),
fHistPosEffMatrixGen(0),
fHistNegEffMatrixGen(0),
fHistEmptyEvent(0),
fHistCentrality(0),
fHistVertex(0),
fHistdEdxTPC(0),
fHistArmPod(0)
{
  //
  //         standard constructur which should be used
  //
   /* fast compilation test
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->AddIncludePath("-I$ALICE_ROOT/include"); 
  gSystem->AddIncludePath("-I$ALICE_PHYSICS/include");
  .L /home/marsland/Desktop/RUN_ON_GRID/Ebye/code/AliAnalysisTaskEbyeIterPID.cxx++
  .L /u/marsland/PHD/macros/marsland_EbyeRatios/AliAnalysisTaskEbyeIterPID.cxx++
  .L /lustre/nyx/alice/users/marsland/train/trunk/marsland_EbyeRatios/AliAnalysisTaskEbyeIterPID.cxx++
  */
  std::cout << " Info::marsland:===================================================================================="<< std::endl;
  std::cout << " Info::marsland:===================================================================================="<< std::endl;
  std::cout << " Info::marsland:===================================================================================="<< std::endl;
  std::cout << " Info::marsland:***************** CONSTRUCTOR CALLED: AliAnalysisTaskEbyeIterPID  *****************"<< std::endl;
  std::cout << " Info::marsland:===================================================================================="<< std::endl;
  std::cout << " Info::marsland:===================================================================================="<< std::endl;
  std::cout << " Info::marsland:===================================================================================="<< std::endl;
  // ==========================================
  //
  //   TFile cutGFile("/lustre/nyx/alice/users/marsland/pFluct/ArmPodGraphicalCuts_Tight.root");
  //   fPionCutG       = (TCutG*)cutGFile.Get("fPionCutG");
  //   fAntiProtonCutG = (TCutG*)cutGFile.Get("fAntiProtonCutG");
  //   fProtonCutG     = (TCutG*)cutGFile.Get("fProtonCutG");
  //   fAntiProtonCutG -> SetName("fAntiProtonCutG");
  //   fAntiProtonCutG -> SetVarX("alfa");
  //   fAntiProtonCutG -> SetVarY("qt");
  //   fProtonCutG     -> SetName("fProtonCutG");
  //   fProtonCutG     -> SetVarX("alfa");
  //   fProtonCutG     -> SetVarY("qt");
  //   fPionCutG       -> SetName("fPionCutG");
  //   fPionCutG       -> SetVarX("alfa");
  //   fPionCutG       -> SetVarY("qt");
  //   
  // ==========================================
  // Initialize arrays
  for (Int_t ires=0; ires<2; ires++){
      for (Int_t imom=0; imom<4; imom++){
          for (Int_t icent=0; icent<20; icent++){
              for (Int_t ieta=0; ieta<20; ieta++){
                  fPiFirstMoments[ires][imom][icent][ieta]=0.;
                  fKaFirstMoments[ires][imom][icent][ieta]=0.;
                  fPrFirstMoments[ires][imom][icent][ieta]=0.;
                  fLaFirstMoments[ires][imom][icent][ieta]=0.;
                  fChFirstMoments[ires][imom][icent][ieta]=0.;
              }
          }
      }
  }
  for (Int_t ibin=0;ibin<3;ibin++) myBin[ibin] = 0;
  for (Int_t ibin=0;ibin<3;ibin++) myBinMC[ibin] = 0; 
  // ==========================================
  //
  // Define outputs
  if (!fdEdxCheck) Initialize();
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
  DefineOutput(3, TTree::Class());
  DefineOutput(4, TTree::Class());
  DefineOutput(5, TTree::Class());
  DefineOutput(6, TTree::Class());
  DefineOutput(7, TTree::Class());
  DefineOutput(8, TTree::Class());
  DefineOutput(9, TTree::Class());
  DefineOutput(10, TTree::Class());
  DefineOutput(11, TTree::Class());
  DefineOutput(12, TTree::Class());
  DefineOutput(13, TTree::Class());
  DefineOutput(14, TTree::Class());
  DefineOutput(15, TTree::Class());
  // ==========================================

}
//________________________________________________________________________
AliAnalysisTaskEbyeIterPID::~AliAnalysisTaskEbyeIterPID() {
  
  //
  // Destructor
  //
  std::cout << " Info::marsland: ===== In the Destructor ===== " << std::endl;
  if (fHistPosEffMatrixRec) delete fHistPosEffMatrixRec;   
  if (fHistNegEffMatrixRec) delete fHistNegEffMatrixRec;  
  if (fHistPosEffMatrixGen) delete fHistPosEffMatrixGen;   
  if (fHistNegEffMatrixGen) delete fHistNegEffMatrixGen;  
  if (fHistdEdxTPC)         delete fHistdEdxTPC;           
  if (fHistEmptyEvent)      delete fHistEmptyEvent;       
  if (fHistCentrality)      delete fHistCentrality;        
  if (fHistVertex)          delete fHistVertex;           
  if (fHistArmPod)          delete fHistArmPod;           
  if (fhEta)                delete fhEta;                 
  if (fhCent)               delete fhCent;                
  if (fhPtot)               delete fhPtot;                
  if (fhnExpected)          delete fhnExpected;           
  if (fhnCleanEl)           delete fhnCleanEl;           
  if (fhnCleanKa)           delete fhnCleanKa;            
  if (fhnCleanDe)           delete fhnCleanDe;           
  if (fhndEdx)              delete fhndEdx;              

}
// 
// ---------------------------------------------------------------------------------
//                                     Functions
// ---------------------------------------------------------------------------------
// 
//________________________________________________________________________
void AliAnalysisTaskEbyeIterPID::Initialize()
{
  //
  // updating parameters in case of changes (standard cuts and the eta window)
  //
  std::cout << " Info::marsland: ===== In the Initialize ===== " << std::endl;
  if (fRunFastSimulation) { std::cout << " Info::marsland: !!! We are running fast simulation return !!! " << std::endl; return; }
  if (fRunFastHighMomentCal) { std::cout << " Info::marsland: !!! We are running fast high moment calculation return !!! " << std::endl; return; }
  AliInfoClass("Creating track cuts for ITS+TPC (2010 definition).");
  // fESDtrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE,1);
  
  fESDtrackCuts = new AliESDtrackCuts;
  fESDtrackCuts -> SetEtaRange(fEtaDown,fEtaUp);
  fESDtrackCuts -> SetPtRange(0.1,1e10);   
  //
  // ------------------------------------------------
 
  if (fMCtrue) {
  // (-2-)  TPC crossed rows cut  
  if (fSystCrossedRows==-1) fESDtrackCuts->SetMinNCrossedRowsTPC(60);
  if (fSystCrossedRows== 0) fESDtrackCuts->SetMinNCrossedRowsTPC(80);
  if (fSystCrossedRows== 1) fESDtrackCuts->SetMinNCrossedRowsTPC(100);
  // ------------------------------------------------
  //
  // ------------------------------------------------
  // (-3-)  Chi2 cut
  if (fSystChi2==-1) fESDtrackCuts->SetMaxChi2PerClusterTPC(3);             
  if (fSystChi2== 0) fESDtrackCuts->SetMaxChi2PerClusterTPC(4);             
  if (fSystChi2== 1) fESDtrackCuts->SetMaxChi2PerClusterTPC(5);        
  // ------------------------------------------------
  // 
  // ------------------------------------------------
  // (-4-)  DCAxy cut
  if (fSystDCAxy==-1) fESDtrackCuts->SetMaxDCAToVertexXYPtDep("0.0156+0.0300/pt^1.01");
  if (fSystDCAxy== 0) fESDtrackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
  if (fSystDCAxy== 1) fESDtrackCuts->SetMaxDCAToVertexXYPtDep("0.0208+0.0400/pt^1.01");
  // ------------------------------------------------
  // ------------------------------------------------
  // vertex z cut
  fESDtrackCuts->SetMaxDCAToVertexZ(2);       
  // ------------------------------------------------
  //
  // ------------------------------------------------
  // require ITS pixels
  if (fIncludeITS){     // Reason for the empty events
    fESDtrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny); // Reason for the structure in phi
      }
  }
  
  // Fixed cuts which are not considered in the systemati checks
  fESDtrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);  
  fESDtrackCuts->SetMaxChi2TPCConstrainedGlobal(36);
  fESDtrackCuts->SetAcceptKinkDaughters(kFALSE);
  fESDtrackCuts->SetRequireTPCRefit(kTRUE);
  fESDtrackCuts->SetRequireITSRefit(kTRUE);               // always use ITS refit
  fESDtrackCuts->SetMaxFractionSharedTPCClusters(0.4);    // ?? FROM MARIAN 
  // ------------------------------------------------
  //
  // ------------------------------------------------
  // Vertex restrictions
  fESDtrackCuts->SetMaxChi2PerClusterITS(36); 
  fESDtrackCuts->SetDCAToVertex2D(kFALSE);
  fESDtrackCuts->SetRequireSigmaToVertex(kFALSE);
  //
  // ------------------------------------------------
  // ------- track cuts to be used for v0s ----------
  // ------------------------------------------------
  //
  fESDtrackCutsCleanSamp = new AliESDtrackCuts("AliESDtrackCutsV0","AliESDtrackCutsV0");
  fESDtrackCutsCleanSamp -> SetEtaRange(-1.5,1.5);
  fESDtrackCutsCleanSamp -> SetPtRange(0.1,1e10);
  fESDtrackCutsCleanSamp -> SetMinNCrossedRowsTPC(70);
  fESDtrackCutsCleanSamp -> SetRequireTPCRefit(kTRUE);
  fESDtrackCutsCleanSamp -> SetAcceptKinkDaughters(kFALSE);
  fESDtrackCutsCleanSamp -> SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
  fESDtrackCutsCleanSamp -> SetMaxChi2PerClusterITS(36);
  fESDtrackCutsCleanSamp -> SetMaxFractionSharedTPCClusters(0.4);
  // ------------------------------------------------
  // V0 selection
  // ------------------------------------------------
  fESDtrackCutsV0   = new AliESDv0Cuts("AliESDCutsV0","AliESDCutsV0");
  fESDtrackCutsV0   ->SetMaxDcaV0Daughters(1.0);
  // ------------------------------------------------
  //
  std::cout << " Info::marsland: ===================================================== " << std::endl;
  std::cout << " Info::marsland: =============== Summary of Track Cuts =============== " << std::endl;
  std::cout << " Info::marsland: ===================================================== " << std::endl;
  fESDtrackCuts->Dump();
}
//________________________________________________________________________
void AliAnalysisTaskEbyeIterPID::UserCreateOutputObjects() 
{
  //
  // Create output histograms, trees of the analysis (called once)
  //
  std::cout << " Info::marsland: ===== In the UserCreateOutputObjects ===== " << std::endl;
  // ------------  setup PIDCombined  ---------------  
  fPIDCombined=new AliPIDCombined;
  fPIDCombined->SetDefaultTPCPriors();
  fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTPC);
 
  // **********************   Input handler to get the PID object *********************
  if (!(fRunFastSimulation || fRunFastHighMomentCal)) {
    AliAnalysisManager* man = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler* inputHandler = dynamic_cast<AliInputEventHandler*>(man->GetInputEventHandler());
    if (!inputHandler)
      AliFatal("Input handler needed");
    else {
      fPIDResponse = inputHandler->GetPIDResponse();       // PID response object
      if (!fPIDResponse) std::cout << " Info::marsland: ======= PIDResponse object was not created ====== " << std::endl;
    }
  }
  // ************************************************************************ 
  //
  OpenFile(1);  // OpenFile for hists  ============================================== 
  fListHist = new TList();
  fListHist->SetOwner(kTRUE);
  
  // ***************   THnSparseF holding the PIDresponse info   ***********************
  Int_t dEdxnBins      = Int_t((fdEdxUp-fdEdxDown)/fdEdxBinWidth);
  Int_t dEdxnBinsClean = Int_t((fdEdxCleanUp-fdEdxDown)/fdEdxBinWidth);
  // 0 --> assumed particle: 0. electron, 1. pion, 2. kaon, 3. proton, 5. deuteron
  // 1 --> sign
  // 2 --> Centrality
  // 3 --> eta
  // 4 --> momentum
  // 5 --> Expected Sigma of a given track
  // 6 --> Expected Mean of a given track
  const Int_t nExpectedbins = 7;
  //                                         0    1,    2,                 3,           4,        5,       6    
  Int_t   binsExpected[nExpectedbins]  = {   5,   1,  fnCentbinsData,  fnEtaBins,   fnMomBins,   600,   4000};   
  Double_t xminExpected[nExpectedbins] = {  0., -2.,   0.,             fEtaDown,    fMomDown,     1.,   fdEdxDown};
  Double_t xmaxExpected[nExpectedbins] = {  5.,  2.,  80.,             fEtaUp,      fMomUp,      61.,   fdEdxUp};
  TString axisNameExpected[nExpectedbins]   = {"particleType","sign","Centrality"    ,"eta" ,"momentum" ,"ExSigma","ExMean"};
  TString axisTitleExpected[nExpectedbins]  = {"particleType","sign","Centrality [%]","#eta","#it{p} (GeV/#it{c})", "#sigma","#mu"};
  fhnExpected = new THnSparseF("hExpected","Expected hists",7,binsExpected,xminExpected,xmaxExpected);
  fhnExpected->GetAxis(2)->Set(fnCentbinsData-1,fxCentBins);
  for (Int_t iaxis=0; iaxis<nExpectedbins;iaxis++){
    fhnExpected->GetAxis(iaxis)->SetName(axisNameExpected[iaxis]);
    fhnExpected->GetAxis(iaxis)->SetTitle(axisTitleExpected[iaxis]);
  }
  // ************************************************************************ 
  //
  // ************************** dEdx histograms ***************************** 
  // 0 --> sign
  // 1 --> Centrality
  // 2 --> eta
  // 3 --> momentum
  // 4 --> TPC dEdx
  const Int_t nhistbins = 5;
  if (!fEffMatrix){
      //                               0,         1,           2,           3,            4
      Int_t   binsdEdx[nhistbins]  = { 2,  fnCentbinsData, fnEtaBins,   fnMomBins,    dEdxnBins};
      Double_t xmindEdx[nhistbins] = {-2,  0.,             fEtaDown,    fMomDown,     fdEdxDown};
      Double_t xmaxdEdx[nhistbins] = { 2,  80.,            fEtaUp,      fMomUp,       fdEdxUp};
      fhndEdx= new THnSparseF("hdEdx","Inclusive dEdx Spectrum"  ,nhistbins,binsdEdx,xmindEdx,xmaxdEdx);
      fhndEdx->GetAxis(1)->Set(fnCentbinsData-1,fxCentBins);
      TString axisNamedEdx[nhistbins]  = {"sign", "Centrality", "eta"  ,"momentum"  ,"TPC dEdx"};
      TString axisTitledEdx[nhistbins] = {"sign" , "cent (%)",   "#eta" ,"#it{p} (GeV/#it{c})" ,"TPC d#it{E}/d#it{x} Signal (a.u.)"};
      for (Int_t iaxis=0; iaxis<nhistbins;iaxis++) {
          fhndEdx->GetAxis(iaxis)->SetName(axisNamedEdx[iaxis]);
          fhndEdx->GetAxis(iaxis)->SetTitle(axisTitledEdx[iaxis]);
      }
  }
  // ************************************************************************ 
  //
  // *************** Clean Kaon, Electrons and dEdx histograms ***************** 
  // 0 --> sign
  // 1 --> Centrality
  // 2 --> eta
  // 3 --> momentum
  // 4 --> TPC dEdx
  // Clean Electrons and Deuterons
  //                                  0,         1,           2,           3,            4
  if (fUseThnSparse){
      Int_t   binsCleanEl[nhistbins]  = { 2,  fnCentbinsData,  fnEtaBins,   fnMomBins,     dEdxnBinsClean};
      Double_t xminCleanEl[nhistbins] = {-2,   0.,             fEtaDown,    fMomDown,      fdEdxDown};
      Double_t xmaxCleanEl[nhistbins] = { 2,  80.,             fEtaUp,      fMomUp,        fdEdxCleanUp};
      fhnCleanEl = new THnSparseF("hCleanEl","Clean Electrons",nhistbins,binsCleanEl,xminCleanEl,xmaxCleanEl);
      fhnCleanEl->GetAxis(1)->Set(fnCentbinsData-1,fxCentBins);
      Int_t   binsCleanDe[nhistbins]  = { 2,  fnCentbinsData,  fnEtaBins,   fnMomBins,     dEdxnBinsClean};
      Double_t xminCleanDe[nhistbins] = {-2,   0.,             fEtaDown,    fMomDown,      fdEdxDown};
      Double_t xmaxCleanDe[nhistbins] = { 2,  80.,             fEtaUp,      fMomUp,        fdEdxCleanUp};
      fhnCleanDe = new THnSparseF("hCleanDe","Clean Deuterons",nhistbins,binsCleanDe,xminCleanDe,xmaxCleanDe);
      fhnCleanDe->GetAxis(1)->Set(fnCentbinsData-1,fxCentBins);
      
  // Clean Kaons
  Int_t   binsCleanKa[nhistbins]  = { 2,  fnCentbinsData,  fnEtaBins,   fnMomBins,     dEdxnBinsClean};
  Double_t xminCleanKa[nhistbins] = {-2,   0.,             fEtaDown,    fMomDown,      fdEdxDown};
  Double_t xmaxCleanKa[nhistbins] = { 2,  80.,             fEtaUp,      fMomUp,        fdEdxCleanUp};
  fhnCleanKa = new THnSparseF("hCleanKa","Clean Kaons"    ,nhistbins,binsCleanKa,xminCleanKa,xmaxCleanKa);
  fhnCleanKa->GetAxis(1)->Set(fnCentbinsData-1,fxCentBins);
  // Set the branch names
  TString axisNamedEdxClean[nhistbins]   = {"sign" ,"Centrality"     ,"eta"  ,"momentum"  ,"TPC dEdx"};
  TString axisTitledEdxClean[nhistbins]  = {"sign" ,"Centrality [%]" ,"#eta" ,"#it{p} (GeV/#it{c})" ,"TPC d#it{E}/d#it{x} Signal (a.u.)"};
  for (Int_t iaxis=0; iaxis<nhistbins;iaxis++){
        fhnCleanEl->GetAxis(iaxis)->SetName(axisNamedEdxClean[iaxis]);  
        fhnCleanEl->GetAxis(iaxis)->SetTitle(axisTitledEdxClean[iaxis]);
        fhnCleanDe->GetAxis(iaxis)->SetName(axisNamedEdxClean[iaxis]);  
        fhnCleanDe->GetAxis(iaxis)->SetTitle(axisTitledEdxClean[iaxis]);
    fhnCleanKa->GetAxis(iaxis)->SetName(axisNamedEdxClean[iaxis]);  fhnCleanKa->GetAxis(iaxis)->SetTitle(axisTitledEdxClean[iaxis]);
      }
      fListHist->Add(fhndEdx);
      fListHist->Add(fhnCleanEl);
      fListHist->Add(fhnCleanDe);
      fListHist->Add(fHistArmPod);

  }
  // ************************************************************************ 
  //
  // ****************** Efficiency matrix histograms ************************ 
  if(fEffMatrix){
    const Int_t ndim=5;
    const Int_t nEtaBins=(fEtaUp-fEtaDown)*20;
    Int_t nbins0[ndim]  ={3,9,50      ,nEtaBins,50  };
    Double_t xmin0[ndim]={0,0,fMomDown,fEtaDown,0.  };
    Double_t xmax0[ndim]={3,9,fMomUp  ,fEtaUp  ,6.25};
    fHistPosEffMatrixRec  =new THnF("fHistPosEffMatrixRec","fHistPosEffMatrixRec",ndim, nbins0,xmin0,xmax0);
    fHistNegEffMatrixRec  =new THnF("fHistNegEffMatrixRec","fHistNegEffMatrixRec",ndim, nbins0,xmin0,xmax0);
    fHistPosEffMatrixGen  =new THnF("fHistPosEffMatrixGen","fHistPosEffMatrixGen",ndim, nbins0,xmin0,xmax0);
    fHistNegEffMatrixGen  =new THnF("fHistNegEffMatrixGen","fHistNegEffMatrixGen",ndim, nbins0,xmin0,xmax0);
    TString axisNameEff[ndim]  = {"particle"      ,"Centrality"     ,"momentum"      ,"eta"  ,"phi"};
    TString axisTitleEff[ndim] = {"particle type" ,"Centrality (%)" ,"#it{p}_{T} (GeV/#it{c})" ,"#eta" ,"#phi"};
    for (Int_t iEff=0; iEff<ndim;iEff++){
      fHistPosEffMatrixRec->GetAxis(iEff)->SetName(axisNameEff[iEff]);  fHistPosEffMatrixRec->GetAxis(iEff)->SetTitle(axisTitleEff[iEff]);
      fHistNegEffMatrixRec->GetAxis(iEff)->SetName(axisNameEff[iEff]);  fHistNegEffMatrixRec->GetAxis(iEff)->SetTitle(axisTitleEff[iEff]);
      fHistPosEffMatrixGen->GetAxis(iEff)->SetName(axisNameEff[iEff]);  fHistPosEffMatrixGen->GetAxis(iEff)->SetTitle(axisTitleEff[iEff]);
      fHistNegEffMatrixGen->GetAxis(iEff)->SetName(axisNameEff[iEff]);  fHistNegEffMatrixGen->GetAxis(iEff)->SetTitle(axisTitleEff[iEff]);
    }
    fListHist->Add(fHistPosEffMatrixRec);
    fListHist->Add(fHistNegEffMatrixRec);
    fListHist->Add(fHistPosEffMatrixGen);
    fListHist->Add(fHistNegEffMatrixGen);
  }
  // ************************************************************************ 
  //
  // **************************** Event histograms ************************** 
  fHistEmptyEvent     = new TH1F("hEmptyEvent",    "control histogram to count empty events"    , 10, 0., 10.);
  fHistCentrality     = new TH1F("hCentrality",    "control histogram to count number of events", 100, 0., 100.);
  fHistVertex         = new TH1F("hVertex",        "control histogram for vertex Z position"    , 100, -50., 50.);
  fHistArmPod         = new TH2F("hArmPod",        "Armenteros-Podolanski plot"                 , 100,-1.,1., 110,0.,0.22);
  // ************************************************************************ 
  //
  // ************************* make up for hists ****************************
  fHistArmPod->GetXaxis()->SetTitle("#alpha"); 
  fHistArmPod->GetYaxis()->SetTitle("q_{t}"); 
  fHistArmPod->SetMarkerStyle(kFullCircle);
  fHistCentrality    ->GetXaxis()->SetTitle("Centrality (%)"); 
  fHistVertex        ->GetXaxis()->SetTitle("vertexZ (cm)"); 

  //   Add histograms to the TList
  fListHist->Add(fhnExpected);
  fListHist->Add(fhnCleanKa);
  fListHist->Add(fHistEmptyEvent);
  fListHist->Add(fHistCentrality);
  fListHist->Add(fHistVertex);

  //   OpenFile(2);  // OpenFile TIden Tree    ==============================================  
  // TIdentity informations
  fIdenTree = new TTree("fIdenTree",   "Tree for TIdentity analysis");
  fIdenTree -> Branch("sign"   ,&signNew   ,"signNew/I");
  fIdenTree -> Branch("myDeDx" ,&myDeDx    ,"myDeDx/D");
  fIdenTree -> Branch("evtNum" ,&fEventGID);
  fIdenTree -> Branch("myBin"  ,myBin      ,"myBin[3]/I");
  //   fIdenTree -> Branch("px"     ,&fPx );
  //   fIdenTree -> Branch("py"     ,&fPy );
  //   fIdenTree -> Branch("pz"     ,&fPz );

  //   OpenFile(3);  // OpenFile TIdenMC tree   ==============================================  
  // MC TIdentity informations
  fIdenTreeMC = new TTree("fIdenTreeMC", "Tree for MC TIdentity analysis");
  fIdenTreeMC -> Branch("sign"   ,&signNewMC   ,"signNewMC/I");
  fIdenTreeMC -> Branch("myDeDx" ,&myDeDxMC    ,"myDeDxMC/D");
  fIdenTreeMC -> Branch("evtNum" ,&fEventGIDMC);
  fIdenTreeMC -> Branch("myBin"  ,myBinMC      ,"myBinMC[3]/I");
  //   fIdenTreeMC -> Branch("px"     ,&fPxMC );
  //   fIdenTreeMC -> Branch("py"     ,&fPyMC );
  //   fIdenTreeMC -> Branch("pz"     ,&fPzMC );
  
  //   OpenFile(4);  // OpenFile for data tree ==============================================
  // Real Data Tree
  fTree = new TTree("fTree",      "Tree for analysis of PID");
  fTree->Branch("dEdx"       , &fTPCSignal);
  fTree->Branch("ptot"       , &fptot);
  fTree->Branch("eta"        , &fEta);
  fTree->Branch("sign"       , &fSign);
  fTree->Branch("cent"       , &fCentrality);
  fTree->Branch("phi"        , &fPhi);
  fTree->Branch("event"      , &fEventGID);

  //   OpenFile(5);  // OpenFile for armPod tree     ==============================================  
  // Armpod tree to be able to use graphical cut
  fArmPodTree  = new TTree("fArmPodTree",  "Tree for Clean Pion and Proton selection");
  fArmPodTree->Branch("dEdx"          , &fArmPodTPCSignal);
  fArmPodTree->Branch("ptot"          , &fArmPodptot);
  fArmPodTree->Branch("eta"           , &fArmPodEta);
  fArmPodTree->Branch("cent"          , &fArmPodCentrality);
  fArmPodTree->Branch("qt"            , &fQt);
  fArmPodTree->Branch("alfa"          , &fAlfa);
  fArmPodTree->Branch("piTOFnSigma"   , &fPiNSigmasTOF);
  fArmPodTree->Branch("prTOFnSigma"   , &fPrNSigmasTOF);
  if (fFillCuts) {
      fArmPodTree->Branch("piFromK0"      , &fCleanPionsFromK0);
      fArmPodTree->Branch("v0haspixel"    , &fHasV0FirstITSlayer);
      fArmPodTree->Branch("pi0FromK0"     , &fCleanPion0FromK0);
      fArmPodTree->Branch("pi1FromK0"     , &fCleanPion1FromK0);
      fArmPodTree->Branch("pi0FromLambda" , &fCleanPion0FromLambda);
      fArmPodTree->Branch("pi1FromLambda" , &fCleanPion1FromLambda);
      fArmPodTree->Branch("pr0FromLambda" , &fCleanProton0FromLambda);
      fArmPodTree->Branch("pr1FromLambda" , &fCleanProton1FromLambda);
      fArmPodTree->Branch("tr0haspixel"   , &fHasTrack0FirstITSlayer);
      fArmPodTree->Branch("tr1haspixel"   , &fHasTrack1FirstITSlayer);
  } else {
      fArmPodTree->Branch("piFromK0"      , &fCleanPionsFromK0);
      fArmPodTree->Branch("v0haspixel"    , &fHasV0FirstITSlayer);
  }

  //   OpenFile(6);
  fTreeSRedirector = new TTreeSRedirector();
  fTreeMCrec     = ((*fTreeSRedirector)<<"mcRec").GetTree();
  fTreeMCgen     = ((*fTreeSRedirector)<<"mcGen").GetTree();
  fTreeDnchDeta  = ((*fTreeSRedirector)<<"dnchdeta").GetTree();
  fTreeMC        = ((*fTreeSRedirector)<<"fTreeMC").GetTree();
  fTreedEdxCheck = ((*fTreeSRedirector)<<"dEdxCheck").GetTree();
  fTreeBayes     = ((*fTreeSRedirector)<<"bayes").GetTree();
  fTreeCuts      = ((*fTreeSRedirector)<<"cuts").GetTree();
  fTreeMCFullAcc = ((*fTreeSRedirector)<<"fullacc").GetTree();
  fTreeResonance = ((*fTreeSRedirector)<<"resonance").GetTree();
  fTreeMCgenMoms = ((*fTreeSRedirector)<<"mcGenMoms").GetTree();

  //  Send data to container   
  PostData(1, fListHist);
  PostData(2, fIdenTree);
  PostData(3, fIdenTreeMC);
  PostData(4, fTree);
  PostData(5, fArmPodTree);
  PostData(6, fTreeMCrec);
  PostData(7, fTreeMCgen);
  PostData(8, fTreeMC);
  PostData(9, fTreedEdxCheck);
  PostData(10, fTreeBayes);
  PostData(11, fTreeCuts);
  PostData(12, fTreeDnchDeta);
  PostData(13, fTreeMCFullAcc);
  PostData(14, fTreeResonance);
  PostData(15, fTreeMCgenMoms);

  std::cout << " Info::marsland: ===== Out of UserCreateOutputObjects ===== " << std::endl;

}
//________________________________________________________________________
void AliAnalysisTaskEbyeIterPID::UserExec(Option_t *) 
{
  //
  // main event loop
  //
  if (fUseCouts) std::cout << " Info::marsland: ===== In the UserExec ===== " << std::endl;
  // Check Monte Carlo information and other access first:
  AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  if (!eventHandler) fMCtrue = kFALSE;
  //     
  // ======================================================================
  // ========================== See if MC or Real =========================
  // ======================================================================
  // impact parameters to use: 0.0, 3.72, 5.23, 7.31, 8.88, 10.20, 11.38, 12.47, 13.50, 14.5
  // corresponding Centrality:  0     5    10    20    30     40     50    60      70    80
  // AliMCEvent *mcEvent = 0x0;
  if (eventHandler) fMCEvent = eventHandler->MCEvent();
  AliGenEventHeader* genHeader = 0x0;
  if (fMCEvent) genHeader = fMCEvent->GenEventHeader(); 
  if (fMCEvent){
      if(!genHeader){ 
          printf("  Event generator header not available!!!\n"); 
          return; 
      }
  }
  //
  // Get rid of "E-AliESDpid::GetTPCsignalTunedOnData: Tune On Data requested, but MC event not set. Call SetCurrentMCEvent before!" errors
  if (!fPIDResponse && !(fRunFastSimulation || fRunFastHighMomentCal)) fPIDResponse = ((AliESDInputHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->GetESDpid();
  if (fMCEvent && !(fRunFastSimulation || fRunFastHighMomentCal)) fPIDResponse->SetCurrentMCEvent(fMCEvent);  
  //
  if(fMCtrue){
      
    // Get the MC stack 
    fMCStack = fMCEvent->Stack(); 
    if (!fMCStack) { printf("  ERROR: No MC stack available !!!\n"); return;}
    //
    // ========================== MC =========================
    //
    // impact parameters to use: 0.0, 3.72, 5.23, 7.31, 8.88, 10.20, 11.38, 12.47, 13.50, 14.5
    // corresponding Centrality:  0     5    10    20    30     40     50    60      70    80
    //
    Double_t impParArr[10] = {0.0, 3.72, 5.23, 7.31, 8.88, 10.20, 11.38, 12.47, 13.50, 14.5};
    fCentrality = -2;
    AliCentrality    *esdCentrality = 0x0;
    AliMultSelection *MultSelection = 0x0;
    if (fESD) {
        esdCentrality = fESD->GetCentrality();
        MultSelection = (AliMultSelection*) fESD-> FindListObject("MultSelection");
    }
    if (MultSelection) {
        fCentrality = MultSelection->GetMultiplicityPercentile("V0M"); 
    } else if (esdCentrality) {
        fCentrality = esdCentrality->GetCentralityPercentile("V0M");
    } else if (!TMath::IsNaN(((AliGenHijingEventHeader*) genHeader)->ImpactParameter()) && fMCEvent){
        fMCImpactParameter = ((AliGenHijingEventHeader*) genHeader)->ImpactParameter();
        if (fMCImpactParameter>=impParArr[0] && fMCImpactParameter<impParArr[1]) fCentrality=0.25;
        if (fMCImpactParameter>=impParArr[1] && fMCImpactParameter<impParArr[2]) fCentrality=7.5;
        if (fMCImpactParameter>=impParArr[2] && fMCImpactParameter<impParArr[3]) fCentrality=15.;
        if (fMCImpactParameter>=impParArr[3] && fMCImpactParameter<impParArr[4]) fCentrality=25.;
        if (fMCImpactParameter>=impParArr[4] && fMCImpactParameter<impParArr[5]) fCentrality=35.;
        if (fMCImpactParameter>=impParArr[5] && fMCImpactParameter<impParArr[6]) fCentrality=45.;
        if (fMCImpactParameter>=impParArr[6] && fMCImpactParameter<impParArr[7]) fCentrality=55.;
        if (fMCImpactParameter>=impParArr[7] && fMCImpactParameter<impParArr[8]) fCentrality=65.;
        if (fMCImpactParameter>=impParArr[8] && fMCImpactParameter<impParArr[9]) fCentrality=75.;
        if (fMCImpactParameter<impParArr[0]  || fMCImpactParameter>impParArr[9]) fCentrality=-10.;
        if (fUseCouts) std::cout << " Info::marsland: impact parameter = " << fMCImpactParameter << std::endl;
    } else {
        std::cout << " Info::marsland: Error: There is no cent info " << std::endl;
    }
    //
    // Use file name in Hashing to create unique event ID
    //
    TTree *chain = (TChain*)GetInputData(0); 
    if(!chain) { Printf("ERROR: Could not receive input chain"); return; }
    TObjString fileName(chain->GetCurrentFile()->GetName());
    TString amptFileName = fileName.GetString(); 
    fEventGIDMC  = TMath::Abs(Int_t(TMath::Hash(chain->GetCurrentFile()->GetName()))/2);       // uniqe id for file 
    //     fEventGIDMC += TMath::Abs(Int_t(fCentrality)+fEventCountInFile+(1000*fMCImpactParameter));
    fEventGIDMC += TMath::Abs(Int_t(fCentrality)+fEventCountInFile);
    fEventGIDMC  = TMath::Abs(Int_t(TString::Hash(&fEventGIDMC,sizeof(Int_t))));    // uniqe event id for real data 
    fEventGID    = fEventGIDMC;
    if (fUseCouts) {
    std::cout << " Info::marsland: ====================================================================================================== " << std::endl; 
    std::cout << fEventCountInFile << " ----- " << "eventIDMC = " << fEventGIDMC << "   " << chain->GetCurrentFile()->GetName() << std::endl;
    std::cout << " Info::marsland: Centrality = " << fCentrality << " ------ Impact Param = " << fMCImpactParameter << std::endl;
    std::cout << " Info::marsland: ====================================================================================================== " << std::endl; 
    }
    fEventCountInFile++; 
    //
  } 
  
  if (!(fRunFastSimulation || fRunFastHighMomentCal)) {
    //
    // ========================== Real =========================
    //
    fESD = dynamic_cast<AliESDEvent*>( InputEvent() );
    if (!fESD)          { Printf("ERROR: fESD not available"); return; }
    if (!fESDtrackCuts) { Printf("ERROR: fESDtrackCuts not available"); return; }
    //
    // =========== monitor vertex position ============
    //
    Bool_t isVertexOk = kTRUE;
    const AliESDVertex *vertex = fESD->GetPrimaryVertexTracks();
    const AliESDVertex *vertexTPC = fESD->GetPrimaryVertexTracks();
    if(vertex->GetNContributors()<1) {
      vertex = fESD->GetPrimaryVertexSPD();    // SPD vertex
      vertexTPC = fESD->GetPrimaryVertexTPC();    // SPD vertex
      TString vertexType = vertex->GetTitle();    // ??? Put condition Abs(vertex-vertexTPC) as a bool_t into ttree
      if (vertexType.Contains("vertexer: Z") && (vertex->GetDispersion() > 0.04 || vertex->GetZRes() > 0.25))  
	isVertexOk = kFALSE;  
      if (vertex->GetNContributors()<1)  isVertexOk = kFALSE;  
    } 
    fTPCvZ=vertexTPC->GetZ();
    fvZ = vertex->GetZ();
    if (isVertexOk) fHistVertex->Fill(vertex->GetZ());
    Int_t vzCut; // cut on vertex and position before everything else
    if (fSystVz == -1) vzCut = 8;
    if (fSystVz ==  0) vzCut = 10;
    if (fSystVz ==  1) vzCut = 12;
    if (!vertex || !isVertexOk) return;
    else if (TMath::Abs(vertex->GetZ()) > vzCut) return;
    //
    // =========== Centrality definition ============
    //
    TString beamType = fESD->GetBeamType();
    fCentrality = -2;
    if (beamType.CompareTo("A-A") == 0) { // PbPb
      AliCentrality *esdCentrality = fESD->GetCentrality();
      AliMultSelection *MultSelection = (AliMultSelection*) fESD-> FindListObject("MultSelection");
      if (MultSelection) {
          if(MultSelection && fSystCentEstimatetor == -1) fCentrality = MultSelection->GetMultiplicityPercentile("TRK"); 
          if(MultSelection && fSystCentEstimatetor ==  0) fCentrality = MultSelection->GetMultiplicityPercentile("V0M"); 
          if(MultSelection && fSystCentEstimatetor ==  1) fCentrality = MultSelection->GetMultiplicityPercentile("CL1");           
      } else if (esdCentrality) {
          if (esdCentrality && fSystCentEstimatetor == -1) fCentrality = esdCentrality->GetCentralityPercentile("TRK");
          if (esdCentrality && fSystCentEstimatetor ==  0) fCentrality = esdCentrality->GetCentralityPercentile("V0M");
          if (esdCentrality && fSystCentEstimatetor ==  1) fCentrality = esdCentrality->GetCentralityPercentile("CL1");
      } else {
          std::cout << " Info::marsland: Error: There is no cent info " << std::endl;
      }
    }
    fHistCentrality->Fill(fCentrality);  // count events after physics and vertex selection
    //
    // Some event infos
    //
    fEvent           = fESD -> GetEventNumberInFile();
    fMultiplicity    = vertex->GetNContributors();    // fMultiplicity = fESD -> GetNumberOfTracks(); 
    fNContributors   = vertex->GetNContributors();
    fMultiplicityMC  = fMultiplicity;
    //
    // Global id for the event --> which is made with Hashing
    //
    // fEventCountInFile++;                 // uniqe event id within job
    ULong64_t orbitID      = (ULong64_t)fESD->GetOrbitNumber();
    ULong64_t bunchCrossID = (ULong64_t)fESD->GetBunchCrossNumber();
    ULong64_t periodID     = (ULong64_t)fESD->GetPeriodNumber();
    ULong64_t gid          = ((periodID << 36) | (orbitID << 12) | bunchCrossID); 
    fEventGID              = TMath::Abs(Int_t(TString::Hash(&gid,sizeof(Int_t))));    // uniqe event id for real data 
    if (fUseCouts) {
        std::cout << " Info::marsland: =============================================================================================== " << std::endl; 
        std::cout << fEventCountInFile << " ----- " << fCentrality << " === gidreal =  " << gid << " hashed = " << fEventGID << std::endl;
        std::cout << " Info::marsland: =============================================================================================== " << std::endl; 
    }
  }
  //
  // in case small stat is enough
  if (fPercentageOfEvents>0 && (fEventCountInFile%fPercentageOfEvents)==0) return;
  cout << "Event counter = " << fEventCountInFile << endl;
  //
  // ======================================================================
  //   
  // ==========================  Filling part  ============================
  //
  if (fRunFastSimulation && fFillDnchDeta) { FillDnchDeta(); return;}
  if (fRunFastHighMomentCal)               { CalculateFastGenHigherMoments(); return;}
  if (fdEdxCheck)                          { FillTPCdEdxCheck(); return;}
  if (fMCtrue && fEffMatrix )              { FillTPCdEdxMC(); FillTPCdEdxMCEffMatrix();return;}
  if (fMCtrue && fWeakAndMaterial)         { FillTPCdEdxMC(); WeakAndMaterial();return;}  
  if (fMCtrue && fFillDeDxTree)            { FillTPCdEdxMC(); WeakAndMaterial();return;}  
  if (fRunFastSimulation)                  { FastGen(); return;}
  //
  //
  if (fFillDeDxTree){
      FillTPCdEdxReal();
      if (!fMCtrue) { FillCleanElectrons(); FillCleanPions(); } 
    }
  //
  //
  // ======================================================================
  //   
  // ======================== End of Filling part =========================
  
}     
//________________________________________________________________________    
void AliAnalysisTaskEbyeIterPID::FillTPCdEdxReal()
{
  //
  // Fill dEdx information for the TPC and also clean kaon and protons
  //
  if (fUseCouts) std::cout << " Info::marsland: ===== In the FillTPCdEdxReal ===== " << std::endl;
  TTree *chain = (TChain*)GetInputData(0); 
  if(!chain) { Printf("ERROR: Could not receive input chain"); return; }
  TObjString fileName(chain->GetCurrentFile()->GetName());
  TString chunkName = fileName.GetString(); 
  if ((fESD->GetEventNumberInFile())<10) std::cout << fESD->GetEventNumberInFile() << "  chunkName  " <<  chunkName << std::endl;

  // Array to hold bayesian probabilities from combined PID
  Double_t probTPC[AliPID::kSPECIES]={0.};
    
  // -------------------------------------------------------------- 
  // Get the event
  AliVEvent *event=InputEvent();
  // -------------------------------------------------------------- 
  //
  // -------------------------------------------------------------- 
  // check if the event is empty 
  if (CountEmptyEvents(2)<1) return;
  // -------------------------------------------------------------- 
  //
  // -------------------------------------------------------------- 
  // fMultiplicity=event->GetNumberOfTracks();
  // Counter only for the TPC track multiplicity
  Int_t tpcMult = 0;
  for (Int_t itrack=0;itrack<event->GetNumberOfTracks();++itrack){
      AliESDtrack *track = fESD->GetTrack(itrack);
      if (track->IsOn(AliESDtrack::kTPCin)) tpcMult++;
  }
  
  // Main track loop
  Int_t trackCounter = 0;
  for (Int_t itrack=0;itrack<event->GetNumberOfTracks();++itrack) {   // Track loop
    
    fdEdxEl=-10.; fSigmaEl=-10.; 
    fdEdxPi=-10.; fSigmaPi=-10.; 
    fdEdxKa=-10.; fSigmaKa=-10.; 
    fdEdxPr=-10.; fSigmaPr=-10.; 
    fdEdxDe=-10.; fSigmaDe=-10.;
    // -------------------------------------------------------------- 
    // Get the track and check if it is in the TPC
    AliESDtrack *track = fESD->GetTrack(itrack); 
    if (!track->GetInnerParam()) continue;            // check if track in TPC
    // -------------------------------------------------------------- 
    //
    // -------------------------------------------------------------- 
    // Extra cut varibales from marion and jens
    Double_t lengthInActiveZone = track->GetLengthInActiveZone(1,3,230, track->GetBz(),0,0);
    Double_t tpcSignalN         = track->GetTPCsignalN();
    // -------------------------------------------------------------- 
    //
    // -------------------------------------------------------------- 
    // Get some track info
    fTPCSignal = track->GetTPCsignal();
    fptot      = track->GetInnerParam()->GetP();
    fEta       = track->Eta();
    fSign      = track->GetSign();
    fpT        = track->Pt();
    fY         = track->Y();
    fTPCShared = track->GetTPCnclsS();
    fNcl       = track->GetTPCsignalN();
    fPx        = track->Px();
    fPy        = track->Py();
    fPz        = track->Pz(); 
    fTheta     = track->Theta();
    fPhi       = track->Phi();
    fMissingCl = track->GetTPCClusterInfo(3,0,0,159);
    fTPCdEdxInfo = (AliTPCdEdxInfo*)track->GetTPCdEdxInfo();
    Float_t pv[2],cov[3]; 
    track->GetImpactParameters(pv,cov); // p[0]=fD; p[1]=fZ; cov[0]=fCdd; cov[1]=fCdz; cov[2]=fCzz;
    Int_t itsNcls = track->GetITSNcls();
    Double_t tpcCrossedRows        = track->GetTPCCrossedRows();
    Double_t tpcChi2               = track->GetTPCchi2();
    Bool_t   fRequireITSRefit      = fESDtrackCuts->GetRequireITSRefit(); 
    Bool_t   fDCAToVertex2D        = fESDtrackCuts->GetDCAToVertex2D();
    Bool_t   fRequireSigmaToVertex = fESDtrackCuts->GetRequireSigmaToVertex();
    Bool_t   fTrackVeto            = fESDtrackCuts->AcceptTrack(track);
    Bool_t   isFirstITSlayer       = track->HasPointOnITSLayer(0);
    Bool_t   isSecondITSlayer      = track->HasPointOnITSLayer(1);
    Bool_t   fAcceptance = ((fEta>=fEtaDown && fEta<=fEtaUp) 
                       && (fptot>=fMomDown && fptot<=fMomUp) 
		       && (fTPCSignal>=fdEdxDown && fTPCSignal<=fdEdxUp)
		       && (fCentrality>1e-5 && fCentrality<=80));
    static Int_t  runNo = fESD->GetRunNumber();
    Bool_t newITScut = ApplyDCAcutIfNoITSPixel(track);

    // assign each bit of fTrackCutBits as a cut variation
    if (tpcCrossedRows>60)  (fTrackCutBits |= 1 << kNCrossedRowsTPC60);
    if (tpcCrossedRows>80)  (fTrackCutBits |= 1 << kNCrossedRowsTPC80);
    if (tpcCrossedRows>100) (fTrackCutBits |= 1 << kNCrossedRowsTPC100);
    if (fRequireITSRefit)   (fTrackCutBits |= 1 << kRequireITSRefit);
    if (tpcChi2<3)  (fTrackCutBits |= 1 << kMaxChi2PerClusterTPC3);
    if (tpcChi2<4)  (fTrackCutBits |= 1 << kMaxChi2PerClusterTPC4);
    if (tpcChi2<5)  (fTrackCutBits |= 1 << kMaxChi2PerClusterTPC5);
    if (TMath::Abs(pv[0])<0.0156+0.0300/TMath::Power(fpT,1.01)) (fTrackCutBits |= 1 << kMaxDCAToVertexXYPtDepSmall);
    if (TMath::Abs(pv[0])<0.0182+0.0350/TMath::Power(fpT,1.01)) (fTrackCutBits |= 1 << kMaxDCAToVertexXYPtDep);
    if (TMath::Abs(pv[0])<0.0208+0.0400/TMath::Power(fpT,1.01)) (fTrackCutBits |= 1 << kMaxDCAToVertexXYPtDepLarge);
    if (isFirstITSlayer || isSecondITSlayer) (fTrackCutBits |= 1 << kClusterRequirementITS);
    if (pv[1]<1) (fTrackCutBits |= 1 << kVertexZSmall);
    if (pv[1]<2) (fTrackCutBits |= 1 << kVertexZ);
    if (pv[1]<3) (fTrackCutBits |= 1 << kVertexZLarge);
    // -------------------------------------------------------------- 
    //
    // -------------------------------------------------------------- 
    // calculate some variables by hand
    Double_t p[3];
    track->GetPxPyPz(p);
    Double_t momentum = TMath::Sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
    Double_t pt       = TMath::Sqrt(p[0]*p[0] + p[1]*p[1]);
    Double_t mass     = track->GetMass();  // assumed to be pion mass
    Double_t energy   = TMath::Sqrt(mass*mass + momentum*momentum);
    Float_t eta = -100.;
    Float_t y   = -100.;
    if((momentum != TMath::Abs(p[2]))&&(momentum != 0)) eta = 0.5*TMath::Log((momentum + p[2])/(momentum - p[2]));
    if((energy != TMath::Abs(p[2]))&&(energy != 0))     y   = 0.5*TMath::Log((energy + p[2])/(energy - p[2]));
    // -------------------------------------------------------------- 
    //
    // -------------------------------------------------------------- 
    // Default cuts
    if (!fESDtrackCuts->AcceptTrack(track)) continue;    // standard track cuts 
    if (!newITScut)                         continue;    // new ITS cut which applies a DCA restriction to the regions without ITS pixels  
    if (fTightCuts) {
        if (track->GetTPCsignalN()<70) continue;         // from jens to increase pid resolution
        if (track->GetLengthInActiveZone(1,3,230, track->GetBz(),0,0)<120) continue;  // from Marian 
    }                                             
    trackCounter++;
    // Tree for the all cut variables 
    if (fFillCuts) {
      if(!fTreeSRedirector) return;
      (*fTreeSRedirector)<<"cuts"<<
      "dEdx="                 << fTPCSignal            <<         //  dEdx of the track
      "ptot="                 << fptot                 <<         //  total momentum 
      "eta="                  << fEta                  <<         //  eta
      "cent="                 << fCentrality           <<         //  centrality
      "sign="                 << fSign                 <<         //  charge
      "cutBit="               << fTrackCutBits         <<
      "cRows="                << tpcCrossedRows        << 
      "fD="                   << pv[0]                 <<
      "fZ="                   << pv[1]                 <<
      "fCdd="                 << cov[0]                <<
      "fCdz="                 << cov[1]                <<
      "fCzz="                 << cov[2]                <<
      "missCl="               << fMissingCl            <<         //  charge
      "nMult="                << fNContributors        << 
      "tpcMult="              << tpcMult               <<
      "itsNcls="              << itsNcls               <<
      "dEdxInfo.="            << fTPCdEdxInfo          <<
      //       "run="                  << runNo                 <<         //  run Number
      //       "pT="                   << fpT                   <<         //  tranverse momentum
      //       "phi="                  << fPhi                  <<         //  ph
      //       "theta="                << fTheta                <<         //  theta
      //       "Y="                    << fY                    <<         //  rapidity
      //       "px="                   << fPx                   <<         //  x momentum
      //       "py="                   << fPy                   <<         //  y momentum
      //       "pz="                   << fPz                   <<         //  z momentum
      //       "trackVeto="            << fTrackVeto            <<         //  bool: track validity flag
      //       "ncluster="             << fNcl                  <<         //  number of points used for dEdx
      //       "tpcshared="            << fTPCShared            <<         //  shared clusters of the track with adjacent ones
      //       "lactivezone="          << lengthInActiveZone    <<         //  lengthInActiveZone in TPC 
      //       "ITSRefit="             << fRequireITSRefit      <<         //  if ITS refit
      //       "isITSpixel="           << fIsITSpixel01         <<         //  if ITS pixels 
      //       "nITScl="               << fnITSclusters         <<         //  number of its clusters
      //       "prim="                 << fPrimRestriction      <<         //  cut for vertex: should be less than 2 when 
      //       "tpcvz="                << fTPCvZ                <<         //  diff between ITS and TPC standalone vertices 
      //       "spdvz="                << fvZ                   <<         //  diff between ITS and TPC standalone vertices 
      "\n";   
    }  
    fTrackCutBits=0;  // reset the bits for the next track
    // -------------------------------------------------------------- 
    //
    // -------------------------------------------------------------- 
    // Momentum range
    if (fptot<fMomDown || fptot>fMomUp)            continue;        // ptot range to be used 
    //
    // -------------------------------------------------------------- 
    //
    // -------------------------------------------------------------- 
    // Fill the trees
    myDeDx   = fTPCSignal;
    myBin[0] = fhEta ->FindBin(fEta)-1;
    myBin[1] = fhCent->FindBin(fCentrality)-1;
    myBin[2] = fhPtot->FindBin(fptot)-1;
    signNew  = fSign;
    if (!fCleanSamplesOnly && !fMCtrue && fAcceptance && fTIdentity)   fIdenTree -> Fill();
    if (!fCleanSamplesOnly && fFillDeDxTree) fTree -> Fill(); 
    // -------------------------------------------------------------- 
    //
    // -------------------------------------------------------------- 
    // Get the bayesian probabilities from combined PID
    if (fFillBayes){
      UInt_t detUsed = fPIDCombined->ComputeProbabilities(track, fPIDResponse, probTPC);
      Double_t wel = probTPC[AliPID::kElectron];
      Double_t wpi = probTPC[AliPID::kPion];
      Double_t wka = probTPC[AliPID::kKaon];
      Double_t wpr = probTPC[AliPID::kProton];
      if (detUsed != 0) {  // TPC is available --> Fill the tree with bayesian probabilities
	if(!fTreeSRedirector) return;
	(*fTreeSRedirector)<<"bayes"<<
	"dEdx="    << fTPCSignal  <<         //  dEdx of the track
	"phi="     << fPhi        <<         //  phi
	"Y="       << fY          <<         //  rapidity
	"eta="     << fEta        <<         //  eta
	"cent="    << fCentrality <<         //  eta
	"sign="    << fSign       <<         //  charge
	"p="       << fptot       <<         //  x momentum
	"px="      << fPx         <<         //  x momentum
	"py="      << fPy         <<         //  y momentum
	"pz="      << fPz         <<         //  z momentum
	"wel="     << wel <<    
	"wpi="     << wpi <<    
	"wka="     << wka <<    
	"wpr="     << wpr <<    
	"\n";
      } 
    }
    // -------------------------------------------------------------- 
    //
    // --------------------------------------------------------------
    // Fil the dEdx histograms
    Double_t trackdEdx[5] = {Double_t(fSign),fCentrality, fEta,fptot, fTPCSignal};
    if(!fEffMatrix) fhndEdx->Fill(trackdEdx);
    // --------------------------------------------------------------
    //
    // -------------------------------------------------------------- 
    // Get the PID splines
    fNSigmasElTPC = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kElectron);
    fNSigmasPiTPC = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion);
    fNSigmasKaTPC = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon);
    fNSigmasPrTPC = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton);
    fNSigmasDeTPC = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kDeuteron);
    Float_t nSigmasKaTOF = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kKaon);
    Float_t nSigmasDeTOF = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kDeuteron);
    Float_t nSigmasPiTOF = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kPion);
     
    // Electron Expected mean and sigma within 3nsigmaTPC
    if (TMath::Abs(fNSigmasElTPC)<2) {
      fdEdxEl  = fPIDResponse->GetTPCResponse().GetExpectedSignal(track, AliPID::kElectron, AliTPCPIDResponse::kdEdxDefault,fPIDResponse->UseTPCEtaCorrection());
      fSigmaEl = fPIDResponse->GetTPCResponse().GetExpectedSigma(track, AliPID::kElectron, AliTPCPIDResponse::kdEdxDefault,fPIDResponse->UseTPCEtaCorrection());
    }
    
    // Pion Expected mean and sigma within 3nsigmaTPC
    if (TMath::Abs(fNSigmasPiTPC)<2) {
      fdEdxPi  = fPIDResponse->GetTPCResponse().GetExpectedSignal(track, AliPID::kPion,     AliTPCPIDResponse::kdEdxDefault,fPIDResponse->UseTPCEtaCorrection());
      fSigmaPi = fPIDResponse->GetTPCResponse().GetExpectedSigma(track, AliPID::kPion,     AliTPCPIDResponse::kdEdxDefault,fPIDResponse->UseTPCEtaCorrection());
      if (fptot>0.5 && TMath::Abs(nSigmasPiTOF)<2){
        fdEdxPi  = fPIDResponse->GetTPCResponse().GetExpectedSignal(track, AliPID::kPion,     AliTPCPIDResponse::kdEdxDefault,fPIDResponse->UseTPCEtaCorrection());
        fSigmaPi = fPIDResponse->GetTPCResponse().GetExpectedSigma(track, AliPID::kPion,     AliTPCPIDResponse::kdEdxDefault,fPIDResponse->UseTPCEtaCorrection());
      }
    } 
    
    // Kaon Expected mean and sigma within 3nsigmaTPC
    if (TMath::Abs(fNSigmasKaTPC)<2) {
      fdEdxKa  = fPIDResponse->GetTPCResponse().GetExpectedSignal(track, AliPID::kKaon,  AliTPCPIDResponse::kdEdxDefault,fPIDResponse->UseTPCEtaCorrection());
      fSigmaKa = fPIDResponse->GetTPCResponse().GetExpectedSigma(track, AliPID::kKaon,   AliTPCPIDResponse::kdEdxDefault,fPIDResponse->UseTPCEtaCorrection());
      if (fptot>0.5 && TMath::Abs(nSigmasKaTOF)<2){
        fdEdxKa  = fPIDResponse->GetTPCResponse().GetExpectedSignal(track, AliPID::kKaon, AliTPCPIDResponse::kdEdxDefault,fPIDResponse->UseTPCEtaCorrection());
        fSigmaKa = fPIDResponse->GetTPCResponse().GetExpectedSigma(track, AliPID::kKaon,  AliTPCPIDResponse::kdEdxDefault,fPIDResponse->UseTPCEtaCorrection());
      }
    }
      
    // Proton Expected mean and sigma within 3nsigmaTPC
    if (TMath::Abs(fNSigmasPrTPC)<2) {
      fdEdxPr  = fPIDResponse->GetTPCResponse().GetExpectedSignal(track, AliPID::kProton, AliTPCPIDResponse::kdEdxDefault,fPIDResponse->UseTPCEtaCorrection());
      fSigmaPr = fPIDResponse->GetTPCResponse().GetExpectedSigma(track,  AliPID::kProton, AliTPCPIDResponse::kdEdxDefault,fPIDResponse->UseTPCEtaCorrection());
    }
    
     // Deuteron Expected mean and sigma within 3nsigmaTPC
    if (TMath::Abs(fNSigmasDeTPC)<3) {
      fdEdxDe  = fPIDResponse->GetTPCResponse().GetExpectedSignal(track, AliPID::kDeuteron, AliTPCPIDResponse::kdEdxDefault,fPIDResponse->UseTPCEtaCorrection());
      fSigmaDe = fPIDResponse->GetTPCResponse().GetExpectedSigma(track,  AliPID::kDeuteron, AliTPCPIDResponse::kdEdxDefault,fPIDResponse->UseTPCEtaCorrection());
    }
    // -------------------------------------------------------------- 
    //
    // -------------------------------------------------------------- 
    // Fill the THnSparseF for the Expected values form PID response
    Double_t exMean[5]  = {fdEdxEl,  fdEdxPi,  fdEdxKa,  fdEdxPr,  fdEdxDe};
    Double_t exSigma[5] = {fSigmaEl, fSigmaPi, fSigmaKa, fSigmaPr, fSigmaDe};
    if (fdEdxEl>20 || fdEdxPi>20 || fdEdxKa>20 || fdEdxPr>20 || fdEdxDe>20) {
      for (Int_t iPart = 0; iPart< 5; iPart++){
        Double_t weightExpected[7] = {Double_t(iPart),Double_t(fSign),fCentrality,fEta,fptot, exSigma[iPart], exMean[iPart]};
        if(!fEffMatrix) fhnExpected->Fill(weightExpected);
      }
    }
    // -------------------------------------------------------------- 
    //
    // --------------------------------------------------------------       
    // Fill clean kaons
    if ((TMath::Abs(nSigmasKaTOF)<=1.2) && (!fMCtrue)) {
      Double_t nclsTRD      = (Float_t)track->GetTRDncls();
      Double_t TOFSignalDx  = track->GetTOFsignalDx();
      Double_t TOFSignalDz  = track->GetTOFsignalDz();
      if (TOFSignalDz<1.2 && TOFSignalDx<1.2 && nclsTRD>80) {
        Double_t weightCleanKa[5] = {Double_t(fSign),fCentrality,fEta,fptot, fTPCSignal};
        if (!fMCtrue && fUseThnSparse) fhnCleanKa->Fill(weightCleanKa); 
      } 
    } 
    // -------------------------------------------------------------- 
    // --------------------------------------------------------------       
    // Fill clean Deuterons
    if ((TMath::Abs(nSigmasDeTOF)<=3) && TMath::Abs(fNSigmasDeTPC)<3 && (!fMCtrue)) {
      Double_t weightCleanDe[5] = {Double_t(fSign),fCentrality,fEta,fptot, fTPCSignal};
      if (!fMCtrue && fUseThnSparse) fhnCleanDe->Fill(weightCleanDe); 
    } 
    // -------------------------------------------------------------- 
    //
  } // end of track loop

}
//________________________________________________________________________    
void AliAnalysisTaskEbyeIterPID::FillTPCdEdxCheck()
{
  //
  // Fill dEdx information for the TPC and also the clean kaon and protons
  //
  if (fUseCouts) std::cout << " Info::marsland: ===== In the FillTPCdEdxCheck ===== " << std::endl;
  AliVEvent *event=InputEvent();
  //    
  // Main track loop
  //   
  Int_t mult=event->GetNumberOfTracks();
  for (Int_t i=0;i<event->GetNumberOfTracks();++i) {   // Track loop

    // get the track object
    AliESDtrack *track = fESD->GetTrack(i); 
                       
    // get track info
    if (!track->GetInnerParam()) continue;            // check if track in TPC 
    fTPCSignal = track->GetTPCsignal();
    fptot      = track->GetInnerParam()->GetP();
    fEta       = track->Eta();
    
    // Track cuts
    if (fTPCSignal>400)                continue;        
    if (fptot>2.)                      continue;        
    if (fEta<fEtaDown && fEta>fEtaUp)  continue;        
    if (track->GetTPCNcls()<80)        continue;       

    // Fill the tree
    if(!fTreeSRedirector) return;
    (*fTreeSRedirector)<<"dEdxCheck"<<
    "dEdx="     << fTPCSignal <<    // dEdx of mc track
    "ptot="     << fptot <<         // mc momentum
    "eta="      << fEta <<          // mc eta
    "mult="     << mult <<          // multiplicity
    "\n";   
  } 

}
//________________________________________________________________________    
void AliAnalysisTaskEbyeIterPID::FillTPCdEdxMC()
{
  
  //
  // Fill dEdx information for the TPC and also the clean kaon and protons
  //
  if (fUseCouts) std::cout << " Info::marsland: ===== In the FillTPCdEdxMC ===== " << std::endl;
  AliMCEventHandler *eventHandler = dynamic_cast<AliMCEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  TTree *chain = (TChain*)GetInputData(0);
  if(!chain) { Printf("ERROR: Could not receive input chain"); return; }
  TObjString fileName(chain->GetCurrentFile()->GetName());
  Int_t runNumber    = fESD->GetRunNumber();
  Int_t evtNuminFile = fESD->GetEventNumberInFile();
  
  // check if the event is empty 
  if (CountEmptyEvents(7)<1) return;
  
  // ======================================================================
  Int_t tpcMult = 0;
  for (Int_t itrack=0;itrack<fESD->GetNumberOfTracks();++itrack){
      AliESDtrack *track = fESD->GetTrack(itrack);
      if (track->IsOn(AliESDtrack::kTPCin)) tpcMult++;
  }
  // 
  // Fill TIdentity tree for MC + Fill MC closure information
  for(Int_t i = 0; i < fESD->GetNumberOfTracks(); i++) {    
      
      // initialize the dummy particle id
      fElMC =-100.; fPiMC =-100.; fKaMC =-100.; fPrMC =-100.; fDeMC =-100.; fMuMC =-100.; 
      AliESDtrack *trackReal = fESD->GetTrack(i); 
      Int_t lab = TMath::Abs(trackReal->GetLabel());           // avoid from negatif labels, they include some garbage
      Bool_t fTrackVeto = fESDtrackCuts->AcceptTrack(trackReal);
      Bool_t newITScut = ApplyDCAcutIfNoITSPixel(trackReal);
      
      Float_t dcaXYprim, dcaZprim;
      Float_t dcaXYweak, dcaZweak;
      Float_t dcaXYmaterial, dcaZmaterial;
      if(fMCStack->IsSecondaryFromMaterial(lab))  trackReal->GetImpactParameters(dcaXYmaterial, dcaZmaterial);
      if(fMCStack->IsSecondaryFromWeakDecay(lab)) trackReal->GetImpactParameters(dcaXYweak, dcaZweak);
      if(fMCStack->IsPhysicalPrimary(lab))        trackReal->GetImpactParameters(dcaXYprim, dcaZprim);
      
      // MC track cuts
      if (!fWeakAndMaterial) { 
          if (!fMCStack->IsPhysicalPrimary(lab)) continue;  // MC primary track check
          // Track cuts from detector
          if (!newITScut)                                    continue;
          if (!trackReal -> GetInnerParam())                 continue;
          if (!fESDtrackCuts -> AcceptTrack(trackReal))      continue;  // real track cuts 
          // if tight cuts are needed
          if (fTightCuts) if (trackReal->GetTPCsignalN()<70) continue;                                          
          if (fTightCuts) if (trackReal->GetLengthInActiveZone(1,3,230, trackReal->GetBz(),0,0)<120) continue;  
      }
      
      // match the track with mc track
      TParticle *trackMC  = fMCStack->Particle(lab);   
      Int_t pdg           = trackMC->GetPdgCode();
      
      Int_t iPart = -10;
      if (TMath::Abs(pdg) == 11)          iPart = 0; // select el-
      if (TMath::Abs(pdg) == 211)         iPart = 1; // select pi+
      if (TMath::Abs(pdg) == 321)         iPart = 2; // select ka+
      if (TMath::Abs(pdg) == 2212)        iPart = 3; // select pr+
      if (TMath::Abs(pdg) == 1000010020)  iPart = 4; // select de
      if (TMath::Abs(pdg) == 13)          iPart = 5; // select mu-
      if (iPart == -10) continue;
      
      if (iPart == 0 ) fElMC = trackReal->GetTPCsignal();
      if (iPart == 1 ) fPiMC = trackReal->GetTPCsignal();
      if (iPart == 2 ) fKaMC = trackReal->GetTPCsignal();
      if (iPart == 3 ) fPrMC = trackReal->GetTPCsignal();
      if (iPart == 4 ) fDeMC = trackReal->GetTPCsignal();
      if (iPart == 5 ) fMuMC = trackReal->GetTPCsignal();
      
      fEtaMC        = trackReal->Eta();
      fTPCSignalMC  = trackReal->GetTPCsignal();
      fptotMC       = trackReal->GetInnerParam()->GetP();
      fpTMC         = trackReal->Pt();
      fSignMC       = trackReal->GetSign();
      fPxMC         = trackReal->Px();
      fPyMC         = trackReal->Py();
      fPzMC         = trackReal->Pz();
      fMissingCl    = trackReal->GetTPCClusterInfo(3,0,0,159);
      Float_t fPhiMC= trackReal->Phi();
      fTPCdEdxInfo = (AliTPCdEdxInfo*)trackReal->GetTPCdEdxInfo();
      Float_t pv[2],cov[3]; 
      trackReal->GetImpactParameters(pv,cov); // p[0]=fD; p[1]=fZ; cov[0]=fCdd; cov[1]=fCdz; cov[2]=fCzz;
      Int_t itsNcls = trackReal->GetITSNcls();
      
      // Fill the ttree for TIdentity study and fit systematics
      if (fCentrality>=0. && fCentrality<80.){
          myDeDxMC        = fTPCSignalMC;
          myBinMC[0]      = fhEta ->FindBin(fEtaMC)-1;
          myBinMC[1]      = fhCent->FindBin(fCentrality)-1;
          myBinMC[2]      = fhPtot->FindBin(fptotMC)-1;
          signNewMC       = fSignMC;
          if (fTIdentity) fIdenTreeMC->Fill();
          // Fill MC closure tree
          if(!fTreeSRedirector) return;
          if (fWeakAndMaterial){
              (*fTreeSRedirector)<<"fTreeMC"<<
              "dEdx="      << fTPCSignalMC <<    // dEdx of mc track
              "ptot="      << fptotMC <<         // mc momentum
              "pT="        << fpTMC <<         // mc momentum
              "eta="       << fEtaMC <<          // mc eta
              "phi="       << fPhiMC <<          // mc eta
              "cent="      << fCentrality <<     // Centrality
              "sign="      << fSignMC <<         // sign
              "el="        << fElMC <<           // electron dEdx
              "pi="        << fPiMC <<           // pion dEdx
              "ka="        << fKaMC <<           // kaon dEdx
              "pr="        << fPrMC <<           // proton dEdx
              "missCl="    << fMissingCl <<     
              "dEdxInfo.=" << fTPCdEdxInfo << 
              "nMult="     << fNContributors << 
              "tpcMult="   << tpcMult               <<
              "itsNcls="   << itsNcls <<
              "trackVeto=" << fTrackVeto << 
              "dcaXYprim="      << dcaXYprim <<           // electron dEdx
              "dcaXYweak="      << dcaXYweak <<           // electron dEdx
              "dcaXYmaterial="  << dcaXYmaterial <<           // electron dEdx
              "dcaZprim="       << dcaZprim <<           // electron dEdx
              "dcaZweak="       << dcaZweak <<           // electron dEdx
              "dcaZmaterial="   << dcaZmaterial <<           // electron dEdx
              "\n";   
          } else {
              (*fTreeSRedirector)<<"fTreeMC"<<
              "dEdx="      << fTPCSignalMC <<    // dEdx of mc track
              "ptot="      << fptotMC <<         // mc momentum
              "pT="        << fpTMC <<           // mc momentum
              "eta="       << fEtaMC <<          // mc eta
              "phi="       << fPhiMC <<          // mc eta
              "cent="      << fCentrality <<     // Centrality
              "sign="      << fSignMC <<         // sign
              "el="        << fElMC <<           // electron dEdx
              "pi="        << fPiMC <<           // pion dEdx
              "ka="        << fKaMC <<           // kaon dEdx
              "pr="        << fPrMC <<           // proton dEdx
              //           "fZ="        << pv[1] <<
              //           "fD="        << pv[0] <<
              //           "fCdd="      << cov[0] <<
              //           "fCdz="      << cov[1] <<
              //           "fCzz="      << cov[2] <<
              //           "dEdxInfo.=" << fTPCdEdxInfo << 
              //           "missCl="    << fMissingCl <<   
              //           "nMult="     << fNContributors << 
              //           "tpcMult="   << tpcMult <<
              //           "itsNcls="   << itsNcls <<
              "\n";   
          }
      }  
  } // ======= end of track loop =======
  //   
  // ======================================================================
  // ======================================================================
  //   
  // ========= Efficiency Check Eta momentum and Centrality scan ==========
  const Int_t nMoments = 11;
  for (Int_t ieta=0; ieta<fnEtaWinBinsMC; ieta++){
      for (Int_t imom=0; imom<fnMomBinsMC; imom++){
          for (Int_t icent=0; icent<fnCentbinsData-1; icent++){
              
              // -----------------------------------------------------------------------------------------
              // ----------------------------   reconstructed MC particles  ------------------------------
              // -----------------------------------------------------------------------------------------
              // vectors to hold moments
              TVectorF recMoments(nMoments);
              TVectorF recMomentsPos(nMoments);
              TVectorF recMomentsNeg(nMoments);
              TVectorF recMomentsCross(nMoments);
              
              // Moments without resonances
              TVectorF noResRecMoments(nMoments);
              TVectorF noResRecMomentsPos(nMoments);
              TVectorF noResRecMomentsNeg(nMoments);
              TVectorF noResRecMomentsCross(nMoments);
              
              // initialize counters 
              for(Int_t i=0;i<nMoments; i++){  
                  recMoments[i]=0.; 
                  recMomentsPos[i]=0.;  
                  recMomentsNeg[i]=0.; 
                  recMomentsCross[i]=0.;
                  noResRecMoments[i]=0.; 
                  noResRecMomentsPos[i]=0.;  
                  noResRecMomentsNeg[i]=0.; 
                  noResRecMomentsCross[i]=0.;
              }
              
              // loop over tracks
              Double_t centBin = (fcentDownArr[icent]+fcentUpArr[icent])/2.;
              Int_t dataType = 0, sampleNo = 0;
              Int_t nTracks=0, trCount=0;
              Int_t nStackTracks = fESD->GetNumberOfTracks();
              for(Int_t i = 0; i < nStackTracks; i++) {    
                  
                  // initialize the dummy particle id
                  fElMC =-100.; fPiMC =-100.; fKaMC =-100.; fPrMC =-100.; fDeMC =-100.; fMuMC =-100.; fLaMC =-100.; 
                  AliESDtrack *trackReal = fESD->GetTrack(i); 
                  Int_t lab = TMath::Abs(trackReal->GetLabel());           // avoid from negatif labels, they include some garbage
                  fEtaMC    = trackReal->Eta();
                  
                  // MC track cuts
                  if ((fEtaMC<fetaDownArr[ieta]) || (fEtaMC>fetaUpArr[ieta])) continue;  // eta Cut
                  if (!fMCStack->IsPhysicalPrimary(lab)) continue;                         // MC primary track check
                  
                  // Track cuts from detector
                  if (!trackReal -> GetInnerParam()) continue;
                  if (!fESDtrackCuts -> AcceptTrack(trackReal)) continue;  // real track cuts 
                  
                  if (fTightCuts) if (trackReal->GetTPCsignalN()<70) continue;                                          
                  if (fTightCuts) if (trackReal->GetLengthInActiveZone(1,3,230, trackReal->GetBz(),0,0)<120) continue;  
                  
                  // match the track with mc track
                  TParticle *trackMC  = fMCStack->Particle(lab);   
                  Int_t pdg           = trackMC->GetPdgCode();
                  Int_t labMom = trackMC->GetFirstMother();
                  
                  TObjString parName(trackMC->GetName());
                  Int_t pdgMom = 0;
                  TObjString momName="xxx";
                  if ((labMom>=0) && (labMom < nStackTracks)){
                      pdgMom = fMCStack->Particle(labMom)->GetPdgCode();
                      momName = fMCStack->Particle(labMom)->GetName();
                  }
                  
                  // Check if the particle is in the black list of resonances
                  Bool_t acceptRes = kTRUE;
                  for (Int_t ires=0;ires<fnResBins;ires++){
                      if (momName.GetString().Contains(fResonances[ires])) {
                          acceptRes=kFALSE; 
                          break;
                      }
                  }
                  
                  // Check if the mother of the particle is rho or phi
                  //           Bool_t ifFromRho     = pdgMom==113;  // for rho(770)
                  //           Bool_t ifFromPhi     = pdgMom==333;  // for phi(1020)
                  //           Bool_t ifFromDeltaPP = pdgMom==2224; // for delta++
                  //           Bool_t ifFromDeltaP  = pdgMom==2214; // for delta+
                  //           Bool_t ifFromDeltaZ  = pdgMom==2114; // for delta0
                  //           Bool_t ifFromDeltaN  = pdgMom==1114; // for delta-
                  
                  // Identify particle wrt pddg code
                  Int_t iPart = -10;
                  if (TMath::Abs(pdg) == 11)          iPart = 0; // select el-
                  if (TMath::Abs(pdg) == 211)         iPart = 1; // select pi+
                  if (TMath::Abs(pdg) == 321)         iPart = 2; // select ka+
                  if (TMath::Abs(pdg) == 2212)        iPart = 3; // select pr+
                  if (TMath::Abs(pdg) == 1000010020)  iPart = 4; // select de
                  if (TMath::Abs(pdg) == 13)          iPart = 5; // select mu-
                  if (TMath::Abs(pdg) == 3122)        {iPart = 6; fLaMC = iPart;} // select Lambda
                  
                  if (iPart == -10) continue;
                  fptotMC = trackReal->GetInnerParam()->GetP();
                  
                  // additional TOF requirement 
                  if (fIncludeTOF && fptotMC>0.8){
                      Float_t nSigmasPiTOF = fPIDResponse->NumberOfSigmasTOF(trackReal, AliPID::kPion);
                      Float_t nSigmasKaTOF = fPIDResponse->NumberOfSigmasTOF(trackReal, AliPID::kKaon);
                      Float_t nSigmasPrTOF = fPIDResponse->NumberOfSigmasTOF(trackReal, AliPID::kProton);
                      if ( !( 
                          ((TMath::Abs(nSigmasPiTOF)<=3) && iPart==1) ||  
                          ((TMath::Abs(nSigmasKaTOF)<=3) && iPart==2) ||
                          ((TMath::Abs(nSigmasPrTOF)<=3) && iPart==3) 
                      ) ) continue;
                  }
                  
                  if (iPart == 0 ) fElMC = trackReal->GetTPCsignal();
                  if (iPart == 1 ) fPiMC = trackReal->GetTPCsignal();
                  if (iPart == 2 ) fKaMC = trackReal->GetTPCsignal();
                  if (iPart == 3 ) fPrMC = trackReal->GetTPCsignal();
                  if (iPart == 4 ) fDeMC = trackReal->GetTPCsignal();
                  if (iPart == 5 ) fMuMC = trackReal->GetTPCsignal();
                  
                  // count first moments for given Centrality and momentum window
                  if ( (fCentrality>=fcentDownArr[icent]) 
                      && (fCentrality<fcentUpArr[icent]) 
                      && (fptotMC>=fpDownArr[imom]) 
                      && (fptotMC<=fpUpArr[imom]) ) 
                  { 
                      nTracks++;                                // count the particles after track cuts
                      if ( fPiMC>-1 || fKaMC>-1 || fPrMC>-1 || fLaMC>-1) trCount++;
                      if ( fPiMC>-1 ) recMoments[kPi]++;
                      if ( fKaMC>-1 ) recMoments[kKa]++;
                      if ( fPrMC>-1 ) recMoments[kPr]++;
                      if ( fPiMC>-1 && pdg<0) recMomentsNeg[kPi]++;
                      if ( fKaMC>-1 && pdg<0) recMomentsNeg[kKa]++;
                      if ( fPrMC>-1 && pdg<0) recMomentsNeg[kPr]++;
                      if ( fPiMC>-1 && pdg>0) recMomentsPos[kPi]++;
                      if ( fKaMC>-1 && pdg>0) recMomentsPos[kKa]++;
                      if ( fPrMC>-1 && pdg>0) recMomentsPos[kPr]++;
                      
                      // Lambdas for alice
                      if ( fLaMC>-1 ) recMoments[kLa]++;
                      if ( fLaMC>-1 && pdg>0) recMomentsPos[kLa]++;
                      if ( fLaMC>-1 && pdg<0) recMomentsNeg[kLa]++;
                      
                      if ( acceptRes ) {
                          if ( fPiMC>-1 ) noResRecMoments[kPi]++;
                          if ( fKaMC>-1 ) noResRecMoments[kKa]++;
                          if ( fPrMC>-1 ) noResRecMoments[kPr]++;
                          if ( fPiMC>-1 && pdg<0) noResRecMomentsNeg[kPi]++;
                          if ( fKaMC>-1 && pdg<0) noResRecMomentsNeg[kKa]++;
                          if ( fPrMC>-1 && pdg<0) noResRecMomentsNeg[kPr]++;
                          if ( fPiMC>-1 && pdg>0) noResRecMomentsPos[kPi]++;
                          if ( fKaMC>-1 && pdg>0) noResRecMomentsPos[kKa]++;
                          if ( fPrMC>-1 && pdg>0) noResRecMomentsPos[kPr]++;
                          
                          // Lambdas for alice
                          if ( fLaMC>-1 ) noResRecMoments[kLa]++;
                          if ( fLaMC>-1 && pdg>0) noResRecMomentsPos[kLa]++;
                          if ( fLaMC>-1 && pdg<0) noResRecMomentsNeg[kLa]++;
                      }
                  } 
                  
              } // ======= end of track loop =======
              
              // calculate second moments with resonances
              recMoments[kPiPi]=recMoments[kPi]*recMoments[kPi]; 
              recMoments[kKaKa]=recMoments[kKa]*recMoments[kKa]; 
              recMoments[kPrPr]=recMoments[kPr]*recMoments[kPr]; 
              recMoments[kPiKa]=recMoments[kPi]*recMoments[kKa]; 
              recMoments[kPiPr]=recMoments[kPi]*recMoments[kPr]; 
              recMoments[kKaPr]=recMoments[kKa]*recMoments[kPr]; 
              recMomentsNeg[kPiPi]=recMomentsNeg[kPi]*recMomentsNeg[kPi]; 
              recMomentsNeg[kKaKa]=recMomentsNeg[kKa]*recMomentsNeg[kKa]; 
              recMomentsNeg[kPrPr]=recMomentsNeg[kPr]*recMomentsNeg[kPr]; 
              recMomentsNeg[kPiKa]=recMomentsNeg[kPi]*recMomentsNeg[kKa]; 
              recMomentsNeg[kPiPr]=recMomentsNeg[kPi]*recMomentsNeg[kPr]; 
              recMomentsNeg[kKaPr]=recMomentsNeg[kKa]*recMomentsNeg[kPr]; 
              recMomentsPos[kPiPi]=recMomentsPos[kPi]*recMomentsPos[kPi]; 
              recMomentsPos[kKaKa]=recMomentsPos[kKa]*recMomentsPos[kKa]; 
              recMomentsPos[kPrPr]=recMomentsPos[kPr]*recMomentsPos[kPr]; 
              recMomentsPos[kPiKa]=recMomentsPos[kPi]*recMomentsPos[kKa]; 
              recMomentsPos[kPiPr]=recMomentsPos[kPi]*recMomentsPos[kPr]; 
              recMomentsPos[kKaPr]=recMomentsPos[kKa]*recMomentsPos[kPr]; 
              recMomentsCross[kPiPosPiNeg]=recMomentsPos[kPi]*recMomentsNeg[kPi]; 
              recMomentsCross[kPiPosKaNeg]=recMomentsPos[kPi]*recMomentsNeg[kKa]; 
              recMomentsCross[kPiPosPrNeg]=recMomentsPos[kPi]*recMomentsNeg[kPr]; 
              recMomentsCross[kKaPosPiNeg]=recMomentsPos[kKa]*recMomentsNeg[kPi]; 
              recMomentsCross[kKaPosKaNeg]=recMomentsPos[kKa]*recMomentsNeg[kKa]; 
              recMomentsCross[kKaPosPrNeg]=recMomentsPos[kKa]*recMomentsNeg[kPr]; 
              recMomentsCross[kPrPosPiNeg]=recMomentsPos[kPr]*recMomentsNeg[kPi]; 
              recMomentsCross[kPrPosKaNeg]=recMomentsPos[kPr]*recMomentsNeg[kKa]; 
              recMomentsCross[kPrPosPrNeg]=recMomentsPos[kPr]*recMomentsNeg[kPr]; 
              
              // net lambda for Alice
              recMoments[kLaLa]=recMoments[kLa]*recMoments[kLa]; 
              recMomentsNeg[kLaLa]=recMomentsNeg[kLa]*recMomentsNeg[kLa]; 
              recMomentsPos[kLaLa]=recMomentsPos[kLa]*recMomentsPos[kLa]; 
              recMomentsCross[kLaPosLaNeg]=recMomentsPos[kLa]*recMomentsNeg[kLa]; 
              
              // calculate second moments with resonances
              noResRecMoments[kPiPi]=noResRecMoments[kPi]*noResRecMoments[kPi]; 
              noResRecMoments[kKaKa]=noResRecMoments[kKa]*noResRecMoments[kKa]; 
              noResRecMoments[kPrPr]=noResRecMoments[kPr]*noResRecMoments[kPr]; 
              noResRecMoments[kPiKa]=noResRecMoments[kPi]*noResRecMoments[kKa]; 
              noResRecMoments[kPiPr]=noResRecMoments[kPi]*noResRecMoments[kPr]; 
              noResRecMoments[kKaPr]=noResRecMoments[kKa]*noResRecMoments[kPr]; 
              noResRecMomentsNeg[kPiPi]=noResRecMomentsNeg[kPi]*noResRecMomentsNeg[kPi]; 
              noResRecMomentsNeg[kKaKa]=noResRecMomentsNeg[kKa]*noResRecMomentsNeg[kKa]; 
              noResRecMomentsNeg[kPrPr]=noResRecMomentsNeg[kPr]*noResRecMomentsNeg[kPr]; 
              noResRecMomentsNeg[kPiKa]=noResRecMomentsNeg[kPi]*noResRecMomentsNeg[kKa]; 
              noResRecMomentsNeg[kPiPr]=noResRecMomentsNeg[kPi]*noResRecMomentsNeg[kPr]; 
              noResRecMomentsNeg[kKaPr]=noResRecMomentsNeg[kKa]*noResRecMomentsNeg[kPr]; 
              noResRecMomentsPos[kPiPi]=noResRecMomentsPos[kPi]*noResRecMomentsPos[kPi]; 
              noResRecMomentsPos[kKaKa]=noResRecMomentsPos[kKa]*noResRecMomentsPos[kKa]; 
              noResRecMomentsPos[kPrPr]=noResRecMomentsPos[kPr]*noResRecMomentsPos[kPr]; 
              noResRecMomentsPos[kPiKa]=noResRecMomentsPos[kPi]*noResRecMomentsPos[kKa]; 
              noResRecMomentsPos[kPiPr]=noResRecMomentsPos[kPi]*noResRecMomentsPos[kPr]; 
              noResRecMomentsPos[kKaPr]=noResRecMomentsPos[kKa]*noResRecMomentsPos[kPr]; 
              noResRecMomentsCross[kPiPosPiNeg]=noResRecMomentsPos[kPi]*noResRecMomentsNeg[kPi]; 
              noResRecMomentsCross[kPiPosKaNeg]=noResRecMomentsPos[kPi]*noResRecMomentsNeg[kKa]; 
              noResRecMomentsCross[kPiPosPrNeg]=noResRecMomentsPos[kPi]*noResRecMomentsNeg[kPr]; 
              noResRecMomentsCross[kKaPosPiNeg]=noResRecMomentsPos[kKa]*noResRecMomentsNeg[kPi]; 
              noResRecMomentsCross[kKaPosKaNeg]=noResRecMomentsPos[kKa]*noResRecMomentsNeg[kKa]; 
              noResRecMomentsCross[kKaPosPrNeg]=noResRecMomentsPos[kKa]*noResRecMomentsNeg[kPr]; 
              noResRecMomentsCross[kPrPosPiNeg]=noResRecMomentsPos[kPr]*noResRecMomentsNeg[kPi]; 
              noResRecMomentsCross[kPrPosKaNeg]=noResRecMomentsPos[kPr]*noResRecMomentsNeg[kKa]; 
              noResRecMomentsCross[kPrPosPrNeg]=noResRecMomentsPos[kPr]*noResRecMomentsNeg[kPr]; 
              
              // net lambda for Alice
              noResRecMoments[kLaLa]=noResRecMoments[kLa]*noResRecMoments[kLa]; 
              noResRecMomentsNeg[kLaLa]=noResRecMomentsNeg[kLa]*noResRecMomentsNeg[kLa]; 
              noResRecMomentsPos[kLaLa]=noResRecMomentsPos[kLa]*noResRecMomentsPos[kLa]; 
              noResRecMomentsCross[kLaPosLaNeg]=noResRecMomentsPos[kLa]*noResRecMomentsNeg[kLa]; 
              
              
              // fill tree which contains moments
              if(!fTreeSRedirector) return;
              // if there is at least one track in an event fill the tree
              if ( trCount>0 ){
                  (*fTreeSRedirector)<<"mcRec"<<
                  "run="          << runNumber <<               // run number
                  "trCount="      << trCount <<                 // number od identified tracks within the given cent and mom range
                  "isample="      << sampleNo <<                // sample id for subsample method
                  "dataType="     << dataType <<                // data type either MCrec(0) or MCgen(1)
                  "centDown="     << fcentDownArr[icent] <<     // lower edge of cent bin
                  "centUp="       << fcentUpArr[icent] <<       // upper edge of cent bin
                  "centBin="      << centBin <<                 // cent bin
                  "impPar="       << fMCImpactParameter <<      // impact parameter taken from MC event header
                  "pDown="        << fpDownArr[imom] <<         // lower edge of momentum bin
                  "pUp="          << fpUpArr[imom] <<           // upper edge of momentum bin
                  "etaDown="      << fetaDownArr[ieta] <<       // lower edge of eta bin
                  "etaUp="        << fetaUpArr[ieta] <<         // upper edge of eta bin
                  "moment.="      << &recMoments <<             // second moments for particle+antiparticle
                  "momentPos.="   << &recMomentsPos <<          // second moment of positive particles
                  "momentNeg.="   << &recMomentsNeg <<          // second moment of negative particles
                  "momentCross.=" << &recMomentsCross <<        // second moment of unlikesign particles
                  "noResmoment.="      << &noResRecMoments <<             // second moments for particle+antiparticle
                  "noResmomentPos.="   << &noResRecMomentsPos <<          // second moment of positive particles
                  "noResmomentNeg.="   << &noResRecMomentsNeg <<          // second moment of negative particles
                  "noResmomentCross.=" << &noResRecMomentsCross <<        // second moment of unlikesign particles
                  "\n";  
              }
              // -----------------------------------------------------------------------------------------
              // ----------------------------   MC generated pure MC particles  --------------------------
              // -----------------------------------------------------------------------------------------
              // vectors to hold moments 
              TVectorF genMoments(nMoments);
              TVectorF genMomentsPos(nMoments);
              TVectorF genMomentsNeg(nMoments);
              TVectorF genMomentsCross(nMoments);
              
              // Moments without resonances
              TVectorF noResGenMoments(nMoments);
              TVectorF noResGenMomentsPos(nMoments);
              TVectorF noResGenMomentsNeg(nMoments);
              TVectorF noResGenMomentsCross(nMoments);
              
              // initialize counters 
              for(Int_t i=0;i<nMoments; i++){  
                  genMoments[i]=0.; 
                  genMomentsPos[i]=0.;  
                  genMomentsNeg[i]=0.; 
                  genMomentsCross[i]=0.;
                  noResGenMoments[i]=0.; 
                  noResGenMomentsPos[i]=0.;  
                  noResGenMomentsNeg[i]=0.; 
                  noResGenMomentsCross[i]=0.;
              }
              
              // go into event loop
              Int_t nTracksgen=0, trCountgen=0; 
              AliMCParticle *trackMCgen;
              Int_t nMCStackTracks = fMCEvent->GetNumberOfTracks();
              for (Int_t iTrack = 0; iTrack < nMCStackTracks; iTrack++) {    // track loop
                  
                  // initialize the dummy particle id
                  fElMCgen =-100.; fPiMCgen =-100.; fKaMCgen =-100.; fPrMCgen =-100.; fDeMCgen =-100.; fMuMCgen =-100.;
                  trackMCgen = (AliMCParticle *)fMCEvent->GetTrack(iTrack);
                  
                  // apply primary vertex and eta cut
                  if ((trackMCgen->Eta()<fetaDownArr[ieta]) || (trackMCgen->Eta()>fetaUpArr[ieta])) continue;
                  if (!fMCStack->IsPhysicalPrimary(iTrack)) continue;
                  Int_t pdg  = trackMCgen->Particle()->GetPdgCode();
                  Int_t labMom = trackMCgen->GetMother();
                  //           TParticle *trackMC  = fMCStack->Particle(iTrack);   
                  //           Int_t labMom = trackMC->GetFirstMother();
                  TObjString parName(trackMCgen->Particle()->GetName());
                  Int_t pdgMom = 0;
                  TObjString momName="xxx";
                  if ((labMom>=0) && (labMom < nMCStackTracks)){
                      pdgMom = fMCStack->Particle(labMom)->GetPdgCode();
                      momName = fMCStack->Particle(labMom)->GetName();
                  }
          // 
                  // Check if the particle is in the black list of resonances
                  Bool_t acceptRes = kTRUE;
                  for (Int_t ires=0;ires<fnResBins;ires++){
                      if (momName.GetString().Contains(fResonances[ires])) {
                          acceptRes=kFALSE; 
                          break;
                      }
                  }
          //
          // select particle of interest
                  Int_t iPart = -10;
                  if (TMath::Abs(pdg) == 11)          {iPart = 0; fElMCgen = iPart;} // select el-
                  if (TMath::Abs(pdg) == 211)         {iPart = 1; fPiMCgen = iPart;} // select pi+
                  if (TMath::Abs(pdg) == 321)         {iPart = 2; fKaMCgen = iPart;} // select ka+
                  if (TMath::Abs(pdg) == 2212)        {iPart = 3; fPrMCgen = iPart;} // select pr+
                  if (TMath::Abs(pdg) == 1000010020)  {iPart = 4; fDeMCgen = iPart;} // select de
                  if (TMath::Abs(pdg) == 13)          {iPart = 5; fMuMCgen = iPart;} // select mu-
                  if (TMath::Abs(pdg) == 3122)        {iPart = 6; fLaMCgen = iPart;} // select Lambda
                  
                  // dump resonance info
                  if(fEventCountInFile==5) {
                      Bool_t parInterest = (fPiMCgen>-1||fKaMCgen>-1||fPrMCgen>-1||fElMCgen>-1||fLaMCgen>-1) ? kTRUE : kFALSE;
                      if(!fTreeSRedirector) return;
                      (*fTreeSRedirector)<<"resonance"<<
                      "acceptRes="  << acceptRes << 
                      "parInterest=" << parInterest <<          // only pi, ka, and proton
                      "centBin="      << centBin <<                 // cent bin
                      "pDown="    << fpDownArr[imom] <<         // lower edge of momentum bin
                      "etaDown="  << fetaDownArr[ieta] <<       // lower edge of eta bin
                      "pdg="      << pdg      <<         // pdg of prim particle
                      "lab="      << iTrack  <<         // index of prim particle
                      "pdgMom="   << pdgMom   <<         // pdg of mother
                      "labMom="   << labMom   <<         // index of mother
                      "parName.=" << &parName <<         //  full path - file name with ESD
                      "momName.=" << &momName <<         //  full path - file name with ESD
                      "\n"; 
                  }
                  
                  // fill the moments 
                  if (iPart == -10) continue;
                  fptotMCgen = trackMCgen->P();      
                  
                  // count first moments
                  if ((fCentrality>=fcentDownArr[icent])
                      &&(fCentrality<fcentUpArr[icent])
                      &&(fptotMCgen>=fpDownArr[imom])
                      &&(fptotMCgen<=fpUpArr[imom])) 
                  { 
                      nTracksgen++;
                      if ( fPiMCgen>-1 || fKaMCgen>-1 || fPrMCgen>-1 || fLaMCgen>-1) trCountgen++;
                      if ( fPiMCgen>-1 ) genMoments[kPi]++;
                      if ( fKaMCgen>-1 ) genMoments[kKa]++;
                      if ( fPrMCgen>-1 ) genMoments[kPr]++;
                      if ( fPiMCgen>-1 && pdg<0) genMomentsNeg[kPi]++;
                      if ( fKaMCgen>-1 && pdg<0) genMomentsNeg[kKa]++;
                      if ( fPrMCgen>-1 && pdg<0) genMomentsNeg[kPr]++;
                      if ( fPiMCgen>-1 && pdg>0) genMomentsPos[kPi]++;
                      if ( fKaMCgen>-1 && pdg>0) genMomentsPos[kKa]++;
                      if ( fPrMCgen>-1 && pdg>0) genMomentsPos[kPr]++;
                      
                      // Lambdas for alice
                      if ( fLaMCgen>-1 ) genMoments[kLa]++;
                      if ( fLaMCgen>-1 && pdg>0) genMomentsPos[kLa]++;
                      if ( fLaMCgen>-1 && pdg<0) genMomentsNeg[kLa]++;
                      
                      // reject resonances 
                      if ( acceptRes ) {
                          // std::cout << pdg << "   " << parName.GetString() << " ---  " << pdgMom << "  " << momName.GetString() << std::endl;
                          if ( fPiMCgen>-1 ) noResGenMoments[kPi]++;
                          if ( fKaMCgen>-1 ) noResGenMoments[kKa]++;
                          if ( fPrMCgen>-1 ) noResGenMoments[kPr]++;
                          if ( fPiMCgen>-1 && pdg<0) noResGenMomentsNeg[kPi]++;
                          if ( fKaMCgen>-1 && pdg<0) noResGenMomentsNeg[kKa]++;
                          if ( fPrMCgen>-1 && pdg<0) noResGenMomentsNeg[kPr]++;
                          if ( fPiMCgen>-1 && pdg>0) noResGenMomentsPos[kPi]++;
                          if ( fKaMCgen>-1 && pdg>0) noResGenMomentsPos[kKa]++;
                          if ( fPrMCgen>-1 && pdg>0) noResGenMomentsPos[kPr]++;
                          
                          // Lambdas for alice
                          if ( fLaMCgen>-1 ) noResGenMoments[kLa]++;
                          if ( fLaMCgen>-1 && pdg>0) noResGenMomentsPos[kLa]++;
                          if ( fLaMCgen>-1 && pdg<0) noResGenMomentsNeg[kLa]++;
                      }
                  }   
              } // ======= end of track loop ======= 
              
              // calculate second moments
              genMoments[kPiPi]=genMoments[kPi]*genMoments[kPi]; 
              genMoments[kKaKa]=genMoments[kKa]*genMoments[kKa]; 
              genMoments[kPrPr]=genMoments[kPr]*genMoments[kPr]; 
              genMoments[kPiKa]=genMoments[kPi]*genMoments[kKa]; 
              genMoments[kPiPr]=genMoments[kPi]*genMoments[kPr]; 
              genMoments[kKaPr]=genMoments[kKa]*genMoments[kPr]; 
              genMomentsNeg[kPiPi]=genMomentsNeg[kPi]*genMomentsNeg[kPi]; 
              genMomentsNeg[kKaKa]=genMomentsNeg[kKa]*genMomentsNeg[kKa]; 
              genMomentsNeg[kPrPr]=genMomentsNeg[kPr]*genMomentsNeg[kPr]; 
              genMomentsNeg[kPiKa]=genMomentsNeg[kPi]*genMomentsNeg[kKa]; 
              genMomentsNeg[kPiPr]=genMomentsNeg[kPi]*genMomentsNeg[kPr]; 
              genMomentsNeg[kKaPr]=genMomentsNeg[kKa]*genMomentsNeg[kPr]; 
              genMomentsPos[kPiPi]=genMomentsPos[kPi]*genMomentsPos[kPi]; 
              genMomentsPos[kKaKa]=genMomentsPos[kKa]*genMomentsPos[kKa]; 
              genMomentsPos[kPrPr]=genMomentsPos[kPr]*genMomentsPos[kPr]; 
              genMomentsPos[kPiKa]=genMomentsPos[kPi]*genMomentsPos[kKa]; 
              genMomentsPos[kPiPr]=genMomentsPos[kPi]*genMomentsPos[kPr]; 
              genMomentsPos[kKaPr]=genMomentsPos[kKa]*genMomentsPos[kPr]; 
              genMomentsCross[kPiPosPiNeg]=genMomentsPos[kPi]*genMomentsNeg[kPi]; 
              genMomentsCross[kPiPosKaNeg]=genMomentsPos[kPi]*genMomentsNeg[kKa]; 
              genMomentsCross[kPiPosPrNeg]=genMomentsPos[kPi]*genMomentsNeg[kPr]; 
              genMomentsCross[kKaPosPiNeg]=genMomentsPos[kKa]*genMomentsNeg[kPi]; 
              genMomentsCross[kKaPosKaNeg]=genMomentsPos[kKa]*genMomentsNeg[kKa]; 
              genMomentsCross[kKaPosPrNeg]=genMomentsPos[kKa]*genMomentsNeg[kPr]; 
              genMomentsCross[kPrPosPiNeg]=genMomentsPos[kPr]*genMomentsNeg[kPi]; 
              genMomentsCross[kPrPosKaNeg]=genMomentsPos[kPr]*genMomentsNeg[kKa]; 
              genMomentsCross[kPrPosPrNeg]=genMomentsPos[kPr]*genMomentsNeg[kPr]; 
              
              // net lambda for Alice
              genMoments[kLaLa]=genMoments[kLa]*genMoments[kLa]; 
              genMomentsNeg[kLaLa]=genMomentsNeg[kLa]*genMomentsNeg[kLa]; 
              genMomentsPos[kLaLa]=genMomentsPos[kLa]*genMomentsPos[kLa]; 
              genMomentsCross[kLaPosLaNeg]=genMomentsPos[kLa]*genMomentsNeg[kLa]; 
              
              // calculate second moments with resonances
              noResGenMoments[kPiPi]=noResGenMoments[kPi]*noResGenMoments[kPi]; 
              noResGenMoments[kKaKa]=noResGenMoments[kKa]*noResGenMoments[kKa]; 
              noResGenMoments[kPrPr]=noResGenMoments[kPr]*noResGenMoments[kPr]; 
              noResGenMoments[kPiKa]=noResGenMoments[kPi]*noResGenMoments[kKa]; 
              noResGenMoments[kPiPr]=noResGenMoments[kPi]*noResGenMoments[kPr]; 
              noResGenMoments[kKaPr]=noResGenMoments[kKa]*noResGenMoments[kPr]; 
              noResGenMomentsNeg[kPiPi]=noResGenMomentsNeg[kPi]*noResGenMomentsNeg[kPi]; 
              noResGenMomentsNeg[kKaKa]=noResGenMomentsNeg[kKa]*noResGenMomentsNeg[kKa]; 
              noResGenMomentsNeg[kPrPr]=noResGenMomentsNeg[kPr]*noResGenMomentsNeg[kPr]; 
              noResGenMomentsNeg[kPiKa]=noResGenMomentsNeg[kPi]*noResGenMomentsNeg[kKa]; 
              noResGenMomentsNeg[kPiPr]=noResGenMomentsNeg[kPi]*noResGenMomentsNeg[kPr]; 
              noResGenMomentsNeg[kKaPr]=noResGenMomentsNeg[kKa]*noResGenMomentsNeg[kPr]; 
              noResGenMomentsPos[kPiPi]=noResGenMomentsPos[kPi]*noResGenMomentsPos[kPi]; 
              noResGenMomentsPos[kKaKa]=noResGenMomentsPos[kKa]*noResGenMomentsPos[kKa]; 
              noResGenMomentsPos[kPrPr]=noResGenMomentsPos[kPr]*noResGenMomentsPos[kPr]; 
              noResGenMomentsPos[kPiKa]=noResGenMomentsPos[kPi]*noResGenMomentsPos[kKa]; 
              noResGenMomentsPos[kPiPr]=noResGenMomentsPos[kPi]*noResGenMomentsPos[kPr]; 
              noResGenMomentsPos[kKaPr]=noResGenMomentsPos[kKa]*noResGenMomentsPos[kPr]; 
              noResGenMomentsCross[kPiPosPiNeg]=noResGenMomentsPos[kPi]*noResGenMomentsNeg[kPi]; 
              noResGenMomentsCross[kPiPosKaNeg]=noResGenMomentsPos[kPi]*noResGenMomentsNeg[kKa]; 
              noResGenMomentsCross[kPiPosPrNeg]=noResGenMomentsPos[kPi]*noResGenMomentsNeg[kPr]; 
              noResGenMomentsCross[kKaPosPiNeg]=noResGenMomentsPos[kKa]*noResGenMomentsNeg[kPi]; 
              noResGenMomentsCross[kKaPosKaNeg]=noResGenMomentsPos[kKa]*noResGenMomentsNeg[kKa]; 
              noResGenMomentsCross[kKaPosPrNeg]=noResGenMomentsPos[kKa]*noResGenMomentsNeg[kPr]; 
              noResGenMomentsCross[kPrPosPiNeg]=noResGenMomentsPos[kPr]*noResGenMomentsNeg[kPi]; 
              noResGenMomentsCross[kPrPosKaNeg]=noResGenMomentsPos[kPr]*noResGenMomentsNeg[kKa]; 
              noResGenMomentsCross[kPrPosPrNeg]=noResGenMomentsPos[kPr]*noResGenMomentsNeg[kPr]; 
              
              // net lambda for Alice
              noResGenMoments[kLaLa]=noResGenMoments[kLa]*noResGenMoments[kLa]; 
              noResGenMomentsNeg[kLaLa]=noResGenMomentsNeg[kLa]*noResGenMomentsNeg[kLa]; 
              noResGenMomentsPos[kLaLa]=noResGenMomentsPos[kLa]*noResGenMomentsPos[kLa]; 
              noResGenMomentsCross[kLaPosLaNeg]=noResGenMomentsPos[kLa]*noResGenMomentsNeg[kLa]; 
              
              // fill tree which contains moments
              dataType = 1, sampleNo = 0;  // dataType-> 1 for MCgen
              if(!fTreeSRedirector) return;
              // if there is at least one track in an event fill the tree
              if ( trCountgen>0 ){   
                  (*fTreeSRedirector)<<"mcGen"<<
                  "run="          << runNumber <<               // run number
                  "trCount="      << trCountgen <<              // number of identified tracks within the given cent and mom range
                  "isample="      << sampleNo <<                // sample id for subsample method
                  "dataType="     << dataType <<                // data type either MCrec(0) or MCgen(1)
                  "centDown="     << fcentDownArr[icent] <<     // lower edge of cent bin
                  "centUp="       << fcentUpArr[icent] <<       // upper edge of cent bin
                  "centBin="      << centBin <<                 // cent bin
                  "impPar="       << fMCImpactParameter <<      // impact parameter taken from MC event header
                  "pDown="        << fpDownArr[imom] <<         // lower edge of momentum bin
                  "pUp="          << fpUpArr[imom] <<           // upper edge of momentum bin
                  "etaDown="      << fetaDownArr[ieta] <<       // lower edge of eta bin
                  "etaUp="        << fetaUpArr[ieta] <<         // upper edge of eta bin
                  "moment.="      << &genMoments <<             // second moments for particle+antiparticle
                  "momentPos.="   << &genMomentsPos <<          // second moment of positive particles
                  "momentNeg.="   << &genMomentsNeg <<          // second moment of negative particles
                  "momentCross.=" << &genMomentsCross <<        // second moment of unlikesign particles
                  "noResmoment.="      << &noResGenMoments <<             // second moments for particle+antiparticle
                  "noResmomentPos.="   << &noResGenMomentsPos <<          // second moment of positive particles
                  "noResmomentNeg.="   << &noResGenMomentsNeg <<          // second moment of negative particles
                  "noResmomentCross.=" << &noResGenMomentsCross <<        // second moment of unlikesign particles
                  "\n";  
              } // tree filling
              
          } // ======= end of Centrality loop ======= 
      }// ======= end of momentum loop ======= 
  } // ======= end of eta loop ======= 
  // 
  // ======================================================================
  // 
  // ======================================================================
   
}
//________________________________________________________________________
void AliAnalysisTaskEbyeIterPID::FastGen()
{
  
  //
  // Fill dEdx information for the TPC and also the clean kaon and protons
  //
  if (fUseCouts) std::cout << " Info::marsland: ===== In the FastGen ===== " << std::endl;
   
  Int_t dataType = 1, sampleNo = 0;  // dataType-> 1 for MCgen
  Int_t evtNuminFile = fMCEvent -> GetEventNumberInFile();

  // ======================================================================
  // ======================================================================
  //   
  // ========= Efficiency Check Eta momentum and Centrality scan ==========
  const Int_t nMoments = 13;
   for (Int_t ieta=0; ieta<fnEtaWinBinsMC; ieta++){
    for (Int_t imom=0; imom<fnMomBinsMC; imom++){
      for (Int_t icent=0; icent<fnCentbinsData-1; icent++){
        
	// vectors to hold moments
	TVectorF genMoments(nMoments);
	TVectorF genMomentsPos(nMoments);
	TVectorF genMomentsNeg(nMoments);
	TVectorF genMomentsCross(nMoments);
	  
        // Moments without resonances
        TVectorF noResGenMoments(nMoments);
	TVectorF noResGenMomentsPos(nMoments);
	TVectorF noResGenMomentsNeg(nMoments);
	TVectorF noResGenMomentsCross(nMoments);
        
        // initialize counters 
	for(Int_t i=0;i<nMoments; i++){  
            genMoments[i]=0.; 
            genMomentsPos[i]=0.;  
            genMomentsNeg[i]=0.; 
            genMomentsCross[i]=0.;
            noResGenMoments[i]=0.; 
            noResGenMomentsPos[i]=0.;  
            noResGenMomentsNeg[i]=0.; 
            noResGenMomentsCross[i]=0.;
        }
        
        Double_t centBin = (fcentDownArr[icent]+fcentUpArr[icent])/2.;
	Float_t nTracksgen=0, trCountgen=0;
	AliMCParticle *trackMCgen;
        Int_t nStackTracks = fMCEvent->GetNumberOfTracks();
        for (Int_t iTrack = 0; iTrack < nStackTracks; iTrack++) {    // track loop
  
	  // initialize the dummy particle id
          fElMCgen =-100.; fPiMCgen =-100.; fKaMCgen =-100.; fPrMCgen =-100.; fDeMCgen =-100.; fMuMCgen =-100.; fLaMCgen =-100.;
          trackMCgen = (AliMCParticle *)fMCEvent->GetTrack(iTrack);
    
	  // apply primary vertex and eta cut
          if ((trackMCgen->Eta()<fetaDownArr[ieta]) || (trackMCgen->Eta()>fetaUpArr[ieta])) continue;
          // iwith or wihout weak decays
          if (fWeakAndMaterial){
	     if ( !(fMCStack->IsPhysicalPrimary(iTrack) || fMCStack->IsSecondaryFromWeakDecay(iTrack)) ) continue;
          } else if (!fMCStack->IsPhysicalPrimary(iTrack)) continue;
          //
          // get the pdg info for maother and daughter  
          Int_t sign = trackMCgen->Particle()->GetPDG()->Charge();
          Int_t pdg  = trackMCgen->Particle()->GetPdgCode();
          Int_t labMom = trackMCgen->GetMother();
          TObjString parName(trackMCgen->Particle()->GetName());
          Int_t pdgMom = 0;
          TObjString momName="xxx";
          if ((labMom>=0) && (labMom < nStackTracks)){
              pdgMom = fMCStack->Particle(labMom)->GetPdgCode();
              momName = fMCStack->Particle(labMom)->GetName();
          }

          // Check if the particle is in the black list of resonances
          Bool_t acceptRes = kTRUE;
          for (Int_t ires=0;ires<fnResBins;ires++){
              
              if (fResonances[ires].Contains("xxx")){
                  // reject all resonances
                  if (!(momName.GetString().Contains(fResonances[ires]))) {acceptRes=kFALSE; break;} 
              } else {
                  // reject resonances in the array
                  if (momName.GetString().Contains(fResonances[ires])) {acceptRes=kFALSE; break;}
              }
          }
        
          Int_t iPart = -10;
          if (TMath::Abs(pdg) == 11)          {iPart = 0; fElMCgen = iPart;} // select el-
          if (TMath::Abs(pdg) == 211)         {iPart = 1; fPiMCgen = iPart;} // select pi+
          if (TMath::Abs(pdg) == 321)         {iPart = 2; fKaMCgen = iPart;} // select ka+
          if (TMath::Abs(pdg) == 2212)        {iPart = 3; fPrMCgen = iPart;} // select pr+
          if (TMath::Abs(pdg) == 1000010020)  {iPart = 4; fDeMCgen = iPart;} // select de
          if (TMath::Abs(pdg) == 13)          {iPart = 5; fMuMCgen = iPart;} // select mu-
          if (TMath::Abs(pdg) == 3122)        {iPart = 6; fLaMCgen = iPart;} // select Lambda

          // dump resonance info
          if(fEventCountInFile==5) {
              Bool_t parInterest = (fPiMCgen>-1||fKaMCgen>-1||fPrMCgen>-1||fElMCgen>-1||fLaMCgen>-1) ? kTRUE : kFALSE;
              if(!fTreeSRedirector) return;
              (*fTreeSRedirector)<<"resonance"<<
              "acceptRes="  << acceptRes << 
              "parInterest=" << parInterest <<          // only pi, ka, and proton
              "centBin="      << centBin <<                 // cent bin
              "pDown="    << fpDownArr[imom] <<         // lower edge of momentum bin
              "etaDown="  << fetaDownArr[ieta] <<       // lower edge of eta bin
              "pdg="      << pdg      <<         // pdg of prim particle
              "lab="      << iTrack  <<         // index of prim particle
              "pdgMom="   << pdgMom   <<         // pdg of mother
              "labMom="   << labMom   <<         // index of mother
              "parName.=" << &parName <<         //  full path - file name with ESD
              "momName.=" << &momName <<         //  full path - file name with ESD
              "\n"; 
          }
          
          // count first moments
	  fptotMCgen = trackMCgen->P();  
          if ((fCentrality>=fcentDownArr[icent])
               &&(fCentrality<fcentUpArr[icent])
               &&(fptotMCgen>=fpDownArr[imom])
               &&(fptotMCgen<=fpUpArr[imom])) 
          { 
	    nTracksgen++;
            // 
            // count charged particles
            if (sign>0 || sign<0) genMoments[kCh]++;
            if (sign>0) genMomentsPos[kCh]++; 
            if (sign<0) genMomentsNeg[kCh]++; 
            if ( acceptRes ) {
                if (sign>0 || sign<0) noResGenMoments[kCh]++;
                if (sign>0) noResGenMomentsPos[kCh]++;
                if (sign<0) noResGenMomentsNeg[kCh]++;
            }   
            // Count identified particles 
            if ( iPart == -10) continue;
	    if ( fPiMCgen>-1 || fKaMCgen>-1 || fPrMCgen>-1 || fLaMCgen>-1) trCountgen++;
            //
            if ( fPiMCgen>-1 ) genMoments[kPi]++;
            if ( fKaMCgen>-1 ) genMoments[kKa]++;
            if ( fPrMCgen>-1 ) genMoments[kPr]++;
            //
	    if ( fPiMCgen>-1 && pdg<0) genMomentsNeg[kPi]++;
            if ( fKaMCgen>-1 && pdg<0) genMomentsNeg[kKa]++;
            if ( fPrMCgen>-1 && pdg<0) genMomentsNeg[kPr]++;
            //
	    if ( fPiMCgen>-1 && pdg>0) genMomentsPos[kPi]++;
            if ( fKaMCgen>-1 && pdg>0) genMomentsPos[kKa]++;
            if ( fPrMCgen>-1 && pdg>0) genMomentsPos[kPr]++; 
            // Lambdas for alice
            if ( fLaMCgen>-1 ) genMoments[kLa]++;
            if ( fLaMCgen>-1 && pdg>0) genMomentsPos[kLa]++;
            if ( fLaMCgen>-1 && pdg<0) genMomentsNeg[kLa]++;
            // reject resonances
            if ( acceptRes ) {
                if ( fPiMCgen>-1 ) noResGenMoments[kPi]++;
                if ( fKaMCgen>-1 ) noResGenMoments[kKa]++;
                if ( fPrMCgen>-1 ) noResGenMoments[kPr]++;
                
                if ( fPiMCgen>-1 && pdg<0) noResGenMomentsNeg[kPi]++;
                if ( fKaMCgen>-1 && pdg<0) noResGenMomentsNeg[kKa]++;
                if ( fPrMCgen>-1 && pdg<0) noResGenMomentsNeg[kPr]++;
                
                if ( fPiMCgen>-1 && pdg>0) noResGenMomentsPos[kPi]++;
                if ( fKaMCgen>-1 && pdg>0) noResGenMomentsPos[kKa]++;
                if ( fPrMCgen>-1 && pdg>0) noResGenMomentsPos[kPr]++;
                // Lambdas for alice
                if ( fLaMCgen>-1 ) noResGenMoments[kLa]++;
                if ( fLaMCgen>-1 && pdg>0) noResGenMomentsPos[kLa]++;
                if ( fLaMCgen>-1 && pdg<0) noResGenMomentsNeg[kLa]++;
            }
          }   
        } // ======= end of track loop ======= 
      
        // calculate second moments
        genMoments[kPiPi]=genMoments[kPi]*genMoments[kPi]; 
	genMoments[kKaKa]=genMoments[kKa]*genMoments[kKa]; 
        genMoments[kPrPr]=genMoments[kPr]*genMoments[kPr]; 
        genMoments[kPiKa]=genMoments[kPi]*genMoments[kKa]; 
        genMoments[kPiPr]=genMoments[kPi]*genMoments[kPr]; 
        genMoments[kKaPr]=genMoments[kKa]*genMoments[kPr]; 
	genMomentsNeg[kPiPi]=genMomentsNeg[kPi]*genMomentsNeg[kPi]; 
	genMomentsNeg[kKaKa]=genMomentsNeg[kKa]*genMomentsNeg[kKa]; 
        genMomentsNeg[kPrPr]=genMomentsNeg[kPr]*genMomentsNeg[kPr]; 
        genMomentsNeg[kPiKa]=genMomentsNeg[kPi]*genMomentsNeg[kKa]; 
        genMomentsNeg[kPiPr]=genMomentsNeg[kPi]*genMomentsNeg[kPr]; 
        genMomentsNeg[kKaPr]=genMomentsNeg[kKa]*genMomentsNeg[kPr]; 
	genMomentsPos[kPiPi]=genMomentsPos[kPi]*genMomentsPos[kPi]; 
	genMomentsPos[kKaKa]=genMomentsPos[kKa]*genMomentsPos[kKa]; 
        genMomentsPos[kPrPr]=genMomentsPos[kPr]*genMomentsPos[kPr];
        genMomentsPos[kPiKa]=genMomentsPos[kPi]*genMomentsPos[kKa]; 
        genMomentsPos[kPiPr]=genMomentsPos[kPi]*genMomentsPos[kPr]; 
        genMomentsPos[kKaPr]=genMomentsPos[kKa]*genMomentsPos[kPr];
	genMomentsCross[kPiPosPiNeg]=genMomentsPos[kPi]*genMomentsNeg[kPi]; 
	genMomentsCross[kPiPosKaNeg]=genMomentsPos[kPi]*genMomentsNeg[kKa]; 
	genMomentsCross[kPiPosPrNeg]=genMomentsPos[kPi]*genMomentsNeg[kPr]; 
	genMomentsCross[kKaPosPiNeg]=genMomentsPos[kKa]*genMomentsNeg[kPi]; 
	genMomentsCross[kKaPosKaNeg]=genMomentsPos[kKa]*genMomentsNeg[kKa]; 
	genMomentsCross[kKaPosPrNeg]=genMomentsPos[kKa]*genMomentsNeg[kPr]; 
	genMomentsCross[kPrPosPiNeg]=genMomentsPos[kPr]*genMomentsNeg[kPi]; 
	genMomentsCross[kPrPosKaNeg]=genMomentsPos[kPr]*genMomentsNeg[kKa]; 
	genMomentsCross[kPrPosPrNeg]=genMomentsPos[kPr]*genMomentsNeg[kPr];
  
        // net lambda for Alice
        genMoments[kLaLa]=genMoments[kLa]*genMoments[kLa]; 
        genMomentsNeg[kLaLa]=genMomentsNeg[kLa]*genMomentsNeg[kLa]; 
        genMomentsPos[kLaLa]=genMomentsPos[kLa]*genMomentsPos[kLa]; 
        genMomentsCross[kLaPosLaNeg]=genMomentsPos[kLa]*genMomentsNeg[kLa]; 

        genMoments[kChCh]=genMoments[kCh]*genMoments[kCh]; 
        genMomentsNeg[kChCh]=genMomentsNeg[kCh]*genMomentsNeg[kCh]; 
        genMomentsPos[kChCh]=genMomentsPos[kCh]*genMomentsPos[kCh]; 
        genMomentsCross[kChPosChNeg]=genMomentsPos[kCh]*genMomentsNeg[kCh]; 

        // calculate second moments with resonances
        noResGenMoments[kPiPi]=noResGenMoments[kPi]*noResGenMoments[kPi]; 
	noResGenMoments[kKaKa]=noResGenMoments[kKa]*noResGenMoments[kKa]; 
        noResGenMoments[kPrPr]=noResGenMoments[kPr]*noResGenMoments[kPr]; 
        noResGenMoments[kPiKa]=noResGenMoments[kPi]*noResGenMoments[kKa]; 
        noResGenMoments[kPiPr]=noResGenMoments[kPi]*noResGenMoments[kPr]; 
        noResGenMoments[kKaPr]=noResGenMoments[kKa]*noResGenMoments[kPr]; 
	noResGenMomentsNeg[kPiPi]=noResGenMomentsNeg[kPi]*noResGenMomentsNeg[kPi]; 
	noResGenMomentsNeg[kKaKa]=noResGenMomentsNeg[kKa]*noResGenMomentsNeg[kKa]; 
        noResGenMomentsNeg[kPrPr]=noResGenMomentsNeg[kPr]*noResGenMomentsNeg[kPr]; 
        noResGenMomentsNeg[kPiKa]=noResGenMomentsNeg[kPi]*noResGenMomentsNeg[kKa]; 
        noResGenMomentsNeg[kPiPr]=noResGenMomentsNeg[kPi]*noResGenMomentsNeg[kPr]; 
        noResGenMomentsNeg[kKaPr]=noResGenMomentsNeg[kKa]*noResGenMomentsNeg[kPr]; 
	noResGenMomentsPos[kPiPi]=noResGenMomentsPos[kPi]*noResGenMomentsPos[kPi]; 
	noResGenMomentsPos[kKaKa]=noResGenMomentsPos[kKa]*noResGenMomentsPos[kKa]; 
        noResGenMomentsPos[kPrPr]=noResGenMomentsPos[kPr]*noResGenMomentsPos[kPr]; 
        noResGenMomentsPos[kPiKa]=noResGenMomentsPos[kPi]*noResGenMomentsPos[kKa]; 
        noResGenMomentsPos[kPiPr]=noResGenMomentsPos[kPi]*noResGenMomentsPos[kPr]; 
        noResGenMomentsPos[kKaPr]=noResGenMomentsPos[kKa]*noResGenMomentsPos[kPr]; 
	noResGenMomentsCross[kPiPosPiNeg]=noResGenMomentsPos[kPi]*noResGenMomentsNeg[kPi]; 
	noResGenMomentsCross[kPiPosKaNeg]=noResGenMomentsPos[kPi]*noResGenMomentsNeg[kKa]; 
	noResGenMomentsCross[kPiPosPrNeg]=noResGenMomentsPos[kPi]*noResGenMomentsNeg[kPr]; 
	noResGenMomentsCross[kKaPosPiNeg]=noResGenMomentsPos[kKa]*noResGenMomentsNeg[kPi]; 
	noResGenMomentsCross[kKaPosKaNeg]=noResGenMomentsPos[kKa]*noResGenMomentsNeg[kKa]; 
	noResGenMomentsCross[kKaPosPrNeg]=noResGenMomentsPos[kKa]*noResGenMomentsNeg[kPr]; 
	noResGenMomentsCross[kPrPosPiNeg]=noResGenMomentsPos[kPr]*noResGenMomentsNeg[kPi]; 
	noResGenMomentsCross[kPrPosKaNeg]=noResGenMomentsPos[kPr]*noResGenMomentsNeg[kKa]; 
	noResGenMomentsCross[kPrPosPrNeg]=noResGenMomentsPos[kPr]*noResGenMomentsNeg[kPr]; 
        
        // net lambda for Alice
        noResGenMoments[kLaLa]=noResGenMoments[kLa]*noResGenMoments[kLa]; 
        noResGenMomentsNeg[kLaLa]=noResGenMomentsNeg[kLa]*noResGenMomentsNeg[kLa]; 
        noResGenMomentsPos[kLaLa]=noResGenMomentsPos[kLa]*noResGenMomentsPos[kLa]; 
        noResGenMomentsCross[kLaPosLaNeg]=noResGenMomentsPos[kLa]*noResGenMomentsNeg[kLa]; 
  
        noResGenMoments[kChCh]=noResGenMoments[kCh]*noResGenMoments[kCh]; 
        noResGenMomentsNeg[kChCh]=noResGenMomentsNeg[kCh]*noResGenMomentsNeg[kCh]; 
        noResGenMomentsPos[kChCh]=noResGenMomentsPos[kCh]*noResGenMomentsPos[kCh]; 
        noResGenMomentsCross[kChPosChNeg]=noResGenMomentsPos[kCh]*noResGenMomentsNeg[kCh]; 
        
  
        // fill tree which contains moments
        if(!fTreeSRedirector) return;
        // if there is at least one track in an event fill the tree
        if ( trCountgen>0 ){   
          (*fTreeSRedirector)<<"mcGen"<<
          "isample="      << sampleNo <<                // sample id for subsample method
          "dataType="     << dataType <<                // data type either MCrec(0) or MCgen(1)
          "centDown="     << fcentDownArr[icent] <<     // lower edge of cent bin
          "centUp="       << fcentUpArr[icent] <<       // upper edge of cent bin
          "impPar="       << fMCImpactParameter <<      // impact parameter taken from MC event header
          "pDown="        << fpDownArr[imom] <<         // lower edge of momentum bin
          "pUp="          << fpUpArr[imom] <<           // upper edge of momentum bin
          "etaDown="      << fetaDownArr[ieta] <<       // lower edge of eta bin
          "etaUp="        << fetaUpArr[ieta] <<         // upper edge of eta bin
          "moment.="      << &genMoments <<             // second moments for particle+antiparticle
          "momentPos.="   << &genMomentsPos <<          // second moment of positive particles
          "momentNeg.="   << &genMomentsNeg <<          // second moment of negative particles
          "momentCross.=" << &genMomentsCross <<        // second moment of unlikesign particles
          "noResmoment.="      << &noResGenMoments <<             // second moments for particle+antiparticle
          "noResmomentPos.="   << &noResGenMomentsPos <<          // second moment of positive particles
          "noResmomentNeg.="   << &noResGenMomentsNeg <<          // second moment of negative particles
          "noResmomentCross.=" << &noResGenMomentsCross <<        // second moment of unlikesign particles
          "\n";  
	} // tree filling
        
      } // ======= end of Centrality loop ======= 
    }// ======= end of momentum loop ======= 
  } // ======= end of eta loop =======
  // 
  // ======================================================================
  // 
  // ======================================================================

     
}
//________________________________________________________________________
void AliAnalysisTaskEbyeIterPID::CalculateFastGenHigherMoments()
{
  
  //
  // Fill dEdx information for the TPC and also the clean kaon and protons
  //
  if (fUseCouts) std::cout << " Info::marsland: ===== In the CalculateFastGenHigherMoments ===== " << std::endl;
   
  Int_t dataType = 1, sampleNo = 0;  // dataType-> 1 for MCgen
  Int_t evtNuminFile = fMCEvent -> GetEventNumberInFile();

  // ======================================================================
  // ======================================================================
  //   
  // ========= Efficiency Check Eta momentum and Centrality scan ==========
  const Int_t nHighMoments = 5;
  const Int_t nMoments = 13;
   for (Int_t ieta=0; ieta<fnEtaWinBinsMC; ieta++){
    for (Int_t imom=0; imom<fnMomBinsMC; imom++){
      for (Int_t icent=0; icent<fnCentbinsData-1; icent++){
        
        // vectors to hold moments 
	TVectorF netPi(nHighMoments);
	TVectorF netKa(nHighMoments);
	TVectorF netPr(nHighMoments);
	TVectorF netLa(nHighMoments);
        TVectorF netCh(nHighMoments);
        TVectorF noResnetPi(nHighMoments);
	TVectorF noResnetKa(nHighMoments);
	TVectorF noResnetPr(nHighMoments);
	TVectorF noResnetLa(nHighMoments);
        TVectorF noResnetCh(nHighMoments);
        
        // initialize counters 
	for(Int_t i=0;i<nHighMoments; i++){  
            netPi[i]=0.; 
            netKa[i]=0.;  
            netPr[i]=0.; 
            netLa[i]=0.;
            netCh[i]=0.;
            noResnetPi[i]=0.; 
            noResnetKa[i]=0.;  
            noResnetPr[i]=0.; 
            noResnetLa[i]=0.;
            noResnetCh[i]=0.;
        }
        
        // vectors to hold moments 
	TVectorF genMomentsPos(nMoments);
	TVectorF genMomentsNeg(nMoments);
        TVectorF genMomentsCross(nMoments);
	TVectorF noResGenMomentsPos(nMoments);
	TVectorF noResGenMomentsNeg(nMoments);
        TVectorF noResGenMomentsCross(nMoments);
        // initialize counters 
	for(Int_t i=0;i<nMoments; i++){  
            genMomentsPos[i]=0.;  
            genMomentsNeg[i]=0.; 
            genMomentsCross[i]=0.; 
            noResGenMomentsPos[i]=0.;  
            noResGenMomentsNeg[i]=0.; 
            noResGenMomentsCross[i]=0.; 
        }
  
        
        Double_t centBin = (fcentDownArr[icent]+fcentUpArr[icent])/2.;
	Float_t nTracksgen=0, trCountgen=0;
	AliMCParticle *trackMCgen;
        Int_t nStackTracks = fMCEvent->GetNumberOfTracks();
        // TRACK LOOP
        for (Int_t iTrack = 0; iTrack < nStackTracks; iTrack++) {    // track loop
  
	  // initialize the dummy particle id
          fElMCgen =-100.; fPiMCgen =-100.; fKaMCgen =-100.; fPrMCgen =-100.; fDeMCgen =-100.; fMuMCgen =-100.; fLaMCgen =-100.;
          trackMCgen = (AliMCParticle *)fMCEvent->GetTrack(iTrack);
          //
	  // apply primary vertex and eta cut
          if ((trackMCgen->Eta()<fetaDownArr[ieta]) || (trackMCgen->Eta()>fetaUpArr[ieta])) continue;
          //
          // iwith or wihout weak decays
          if (fWeakAndMaterial){
	     if ( !(fMCStack->IsPhysicalPrimary(iTrack) || fMCStack->IsSecondaryFromWeakDecay(iTrack)) ) continue;
          } else if (!fMCStack->IsPhysicalPrimary(iTrack)) continue;
          //
          // get the pdg info for maother and daughter
          Int_t sign = trackMCgen->Particle()->GetPDG()->Charge();
          Int_t pdg  = trackMCgen->Particle()->GetPdgCode();
          Int_t labMom = trackMCgen->GetMother();
          TObjString parName(trackMCgen->Particle()->GetName());
          Int_t pdgMom = 0;
          TObjString momName="xxx";
          if ((labMom>=0) && (labMom < nStackTracks)){
              pdgMom = fMCStack->Particle(labMom)->GetPdgCode();
              momName = fMCStack->Particle(labMom)->GetName();
          }
          //
          // Check if the particle is in the black list of resonances
          Bool_t acceptRes = kTRUE;
          for (Int_t ires=0;ires<fnResBins;ires++){
              
              if (fResonances[ires].Contains("xxx")){
                  // reject all resonances
                  if (!(momName.GetString().Contains(fResonances[ires]))) {acceptRes=kFALSE; break;} 
              } else {
                  // reject resonances in the array
                  if (momName.GetString().Contains(fResonances[ires])) {acceptRes=kFALSE; break;}
              }
          }
          
          Int_t iPart = -10;
          if (TMath::Abs(pdg) == 11)          {iPart = 0; fElMCgen = iPart;} // select el-
          if (TMath::Abs(pdg) == 211)         {iPart = 1; fPiMCgen = iPart;} // select pi+
          if (TMath::Abs(pdg) == 321)         {iPart = 2; fKaMCgen = iPart;} // select ka+
          if (TMath::Abs(pdg) == 2212)        {iPart = 3; fPrMCgen = iPart;} // select pr+
          if (TMath::Abs(pdg) == 1000010020)  {iPart = 4; fDeMCgen = iPart;} // select de
          if (TMath::Abs(pdg) == 13)          {iPart = 5; fMuMCgen = iPart;} // select mu-
          if (TMath::Abs(pdg) == 3122)        {iPart = 6; fLaMCgen = iPart;} // select Lambda

          // dump resonance info
          if(fEventCountInFile==5) {
              Bool_t parInterest = (fPiMCgen>-1||fKaMCgen>-1||fPrMCgen>-1||fElMCgen>-1||fLaMCgen>-1) ? kTRUE : kFALSE;
              if(!fTreeSRedirector) return;
              (*fTreeSRedirector)<<"resonance"<<
              "acceptRes="   << acceptRes << 
              "parInterest=" << parInterest <<          // only pi, ka, and proton
              "centBin="     << centBin <<                 // cent bin
              "pDown="       << fpDownArr[imom] <<         // lower edge of momentum bin
              "etaDown="     << fetaDownArr[ieta] <<       // lower edge of eta bin
              "pdg="         << pdg      <<         // pdg of prim particle
              "lab="         << iTrack   <<         // index of prim particle
              "pdgMom="      << pdgMom   <<         // pdg of mother
              "labMom="      << labMom   <<         // index of mother
              "parName.="    << &parName <<         //  full path - file name with ESD
              "momName.="    << &momName <<         //  full path - file name with ESD
              "\n";
          }
          
          // count first moments
	  fptotMCgen = trackMCgen->P();  
          if ((fCentrality>=fcentDownArr[icent])
               &&(fCentrality<fcentUpArr[icent])
               &&(fptotMCgen>=fpDownArr[imom])
               &&(fptotMCgen<=fpUpArr[imom])) 
          { 
	    nTracksgen++;
            // 
            // count charged particles
            if (sign>0) genMomentsPos[kCh]++; 
            if (sign<0) genMomentsNeg[kCh]++; 
            if ( acceptRes ) {
                if (sign>0) noResGenMomentsPos[kCh]++;
                if (sign<0) noResGenMomentsNeg[kCh]++;
            }   
            // Count identified particles 
            if (iPart == -10) continue;
	    if ( fPiMCgen>-1 || fKaMCgen>-1 || fPrMCgen>-1 || fLaMCgen>-1) trCountgen++;
            //
	    if ( fPiMCgen>-1 && pdg<0) genMomentsNeg[kPi]++;
            if ( fKaMCgen>-1 && pdg<0) genMomentsNeg[kKa]++;
            if ( fPrMCgen>-1 && pdg<0) genMomentsNeg[kPr]++;
            //
	    if ( fPiMCgen>-1 && pdg>0) genMomentsPos[kPi]++;
            if ( fKaMCgen>-1 && pdg>0) genMomentsPos[kKa]++;
            if ( fPrMCgen>-1 && pdg>0) genMomentsPos[kPr]++; 
            // Lambdas for alice
            if ( fLaMCgen>-1 && pdg>0) genMomentsPos[kLa]++;
            if ( fLaMCgen>-1 && pdg<0) genMomentsNeg[kLa]++;
            // reject resonances
            if ( acceptRes ) {
                if ( fPiMCgen>-1 && pdg<0) noResGenMomentsNeg[kPi]++;
                if ( fKaMCgen>-1 && pdg<0) noResGenMomentsNeg[kKa]++;
                if ( fPrMCgen>-1 && pdg<0) noResGenMomentsNeg[kPr]++;
                //
                if ( fPiMCgen>-1 && pdg>0) noResGenMomentsPos[kPi]++;
                if ( fKaMCgen>-1 && pdg>0) noResGenMomentsPos[kKa]++;
                if ( fPrMCgen>-1 && pdg>0) noResGenMomentsPos[kPr]++;
                // Lambdas for alice
                if ( fLaMCgen>-1 && pdg>0) noResGenMomentsPos[kLa]++;
                if ( fLaMCgen>-1 && pdg<0) noResGenMomentsNeg[kLa]++;
            }         
          }     
        } // ======= end of track loop ======= 
        
        // moments from Lookup table
        //         std::cout << " Info::marsland: ====================== " << imom << "  " << icent << "  " << ieta << " ====================== " << std::endl;
        //         std::cout << fPiFirstMoments[0][imom][icent][ieta] << std::endl;
        //         std::cout << fKaFirstMoments[0][imom][icent][ieta] << std::endl;
        //         std::cout << fPrFirstMoments[0][imom][icent][ieta] << std::endl;
        //         std::cout << fLaFirstMoments[0][imom][icent][ieta] << std::endl;
        //         std::cout << fChFirstMoments[0][imom][icent][ieta] << std::endl;
        //         
        //         std::cout << fPiFirstMoments[1][imom][icent][ieta] << std::endl;
        //         std::cout << fKaFirstMoments[1][imom][icent][ieta] << std::endl;
        //         std::cout << fPrFirstMoments[1][imom][icent][ieta] << std::endl;
        //         std::cout << fLaFirstMoments[1][imom][icent][ieta] << std::endl;
        //         std::cout << fChFirstMoments[1][imom][icent][ieta] << std::endl;
        //         std::cout << " Info::marsland: ============================================================================================= " << std::endl;
        
        // net lambda for Alice and 
        genMomentsNeg[kLaLa]=genMomentsNeg[kLa]*genMomentsNeg[kLa]; 
        genMomentsPos[kLaLa]=genMomentsPos[kLa]*genMomentsPos[kLa]; 
        genMomentsCross[kLaPosLaNeg]=genMomentsPos[kLa]*genMomentsNeg[kLa]; 

        genMomentsNeg[kChCh]=genMomentsNeg[kCh]*genMomentsNeg[kCh]; 
        genMomentsPos[kChCh]=genMomentsPos[kCh]*genMomentsPos[kCh]; 
        genMomentsCross[kChPosChNeg]=genMomentsPos[kCh]*genMomentsNeg[kCh]; 
        
        netPi[0]=(genMomentsPos[kPi]-genMomentsNeg[kPi]); 
        netKa[0]=(genMomentsPos[kKa]-genMomentsNeg[kKa]); 
        netPr[0]=(genMomentsPos[kPr]-genMomentsNeg[kPr]); 
        netLa[0]=(genMomentsPos[kLa]-genMomentsNeg[kLa]); 
        netCh[0]=(genMomentsPos[kCh]-genMomentsNeg[kCh]); 
        
        netPi[1]=netPi[0]-fPiFirstMoments[0][imom][icent][ieta]; 
        netKa[1]=netKa[0]-fKaFirstMoments[0][imom][icent][ieta]; 
        netPr[1]=netPr[0]-fPrFirstMoments[0][imom][icent][ieta]; 
        netLa[1]=netLa[0]-fLaFirstMoments[0][imom][icent][ieta]; 
        netCh[1]=netCh[0]-fChFirstMoments[0][imom][icent][ieta]; 
        
        netPi[2]=netPi[1]*netPi[1]; 
        netKa[2]=netKa[1]*netKa[1]; 
        netPr[2]=netPr[1]*netPr[1]; 
        netLa[2]=netLa[1]*netLa[1]; 
        netCh[2]=netCh[1]*netCh[1]; 

        netPi[3]=netPi[1]*netPi[1]*netPi[1]; 
        netKa[3]=netKa[1]*netKa[1]*netKa[1]; 
        netPr[3]=netPr[1]*netPr[1]*netPr[1]; 
        netLa[3]=netLa[1]*netLa[1]*netLa[1]; 
        netCh[3]=netCh[1]*netCh[1]*netCh[1]; 
        
        netPi[4]=netPi[1]*netPi[1]*netPi[1]*netPi[1]; 
        netKa[4]=netKa[1]*netKa[1]*netKa[1]*netKa[1]; 
        netPr[4]=netPr[1]*netPr[1]*netPr[1]*netPr[1]; 
        netLa[4]=netLa[1]*netLa[1]*netLa[1]*netLa[1]; 
        netCh[4]=netCh[1]*netCh[1]*netCh[1]*netCh[1]; 

       
        // Moments without resonances
        noResGenMomentsNeg[kLaLa]=noResGenMomentsNeg[kLa]*noResGenMomentsNeg[kLa]; 
        noResGenMomentsPos[kLaLa]=noResGenMomentsPos[kLa]*noResGenMomentsPos[kLa]; 
        noResGenMomentsCross[kLaPosLaNeg]=noResGenMomentsPos[kLa]*noResGenMomentsNeg[kLa]; 

        noResGenMomentsNeg[kChCh]=noResGenMomentsNeg[kCh]*noResGenMomentsNeg[kCh]; 
        noResGenMomentsPos[kChCh]=noResGenMomentsPos[kCh]*noResGenMomentsPos[kCh]; 
        noResGenMomentsCross[kChPosChNeg]=noResGenMomentsPos[kCh]*noResGenMomentsNeg[kCh]; 

        noResnetPi[0]=(noResGenMomentsPos[kPi]-noResGenMomentsNeg[kPi]); 
        noResnetKa[0]=(noResGenMomentsPos[kKa]-noResGenMomentsNeg[kKa]); 
        noResnetPr[0]=(noResGenMomentsPos[kPr]-noResGenMomentsNeg[kPr]); 
        noResnetLa[0]=(noResGenMomentsPos[kLa]-noResGenMomentsNeg[kLa]); 
        noResnetCh[0]=(noResGenMomentsPos[kCh]-noResGenMomentsNeg[kCh]); 

        noResnetPi[1]=noResnetPi[0]-fPiFirstMoments[1][imom][icent][imom]; 
        noResnetKa[1]=noResnetKa[0]-fKaFirstMoments[1][imom][icent][imom]; 
        noResnetPr[1]=noResnetPr[0]-fPrFirstMoments[1][imom][icent][imom]; 
        noResnetLa[1]=noResnetLa[0]-fLaFirstMoments[1][imom][icent][imom]; 
        noResnetCh[1]=noResnetCh[0]-fChFirstMoments[1][imom][icent][imom]; 

        noResnetPi[2]=noResnetPi[1]*noResnetPi[1]; 
        noResnetKa[2]=noResnetKa[1]*noResnetKa[1]; 
        noResnetPr[2]=noResnetPr[1]*noResnetPr[1]; 
        noResnetLa[2]=noResnetLa[1]*noResnetLa[1]; 
        noResnetCh[2]=noResnetCh[1]*noResnetCh[1]; 
        
        noResnetPi[3]=noResnetPi[1]*noResnetPi[1]*noResnetPi[1]; 
        noResnetKa[3]=noResnetKa[1]*noResnetKa[1]*noResnetKa[1]; 
        noResnetPr[3]=noResnetPr[1]*noResnetPr[1]*noResnetPr[1]; 
        noResnetLa[3]=noResnetLa[1]*noResnetLa[1]*noResnetLa[1]; 
        noResnetCh[3]=noResnetCh[1]*noResnetCh[1]*noResnetCh[1]; 
        
        noResnetPi[4]=noResnetPi[1]*noResnetPi[1]*noResnetPi[1]*noResnetPi[1]; 
        noResnetKa[4]=noResnetKa[1]*noResnetKa[1]*noResnetKa[1]*noResnetKa[1]; 
        noResnetPr[4]=noResnetPr[1]*noResnetPr[1]*noResnetPr[1]*noResnetPr[1]; 
        noResnetLa[4]=noResnetLa[1]*noResnetLa[1]*noResnetLa[1]*noResnetLa[1]; 
        noResnetCh[4]=noResnetCh[1]*noResnetCh[1]*noResnetCh[1]*noResnetCh[1]; 
      
        // fill tree which contains moments
        if(!fTreeSRedirector) return;
        // if there is at least one track in an event fill the tree
        if ( trCountgen>0 ){   
          (*fTreeSRedirector)<<"mcGenMoms"<<
          "trCount="      << nTracksgen <<              // number of identified tracks within the given cent and mom range
          "isample="      << sampleNo <<                // sample id for subsample method
          "dataType="     << dataType <<                // data type either MCrec(0) or MCgen(1)
          "centDown="     << fcentDownArr[icent] <<     // lower edge of cent bin
          "centUp="       << fcentUpArr[icent] <<       // upper edge of cent bin
          "centBin="      << centBin <<                 // cent bin
          "impPar="       << fMCImpactParameter <<      // impact parameter taken from MC event header
          "pDown="        << fpDownArr[imom] <<         // lower edge of momentum bin
          "pUp="          << fpUpArr[imom] <<           // upper edge of momentum bin
          "etaDown="      << fetaDownArr[ieta] <<       // lower edge of eta bin
          "etaUp="        << fetaUpArr[ieta] <<         // upper edge of eta bin
          "momPos.="      << &genMomentsPos <<          // second moment of positive particles
          "momNeg.="      << &genMomentsNeg <<          // second moment of negative particles
          "noResmomPos.=" << &noResGenMomentsPos <<     // second moment of positive particles
          "noResmomNeg.=" << &noResGenMomentsNeg <<     // second moment of negative particles
          "netPi.="       << &netPi <<                  // second moments for particle+antiparticle
          "netKa.="       << &netKa <<                  // second moment of positive particles
          "netPr.="       << &netPr <<                  // second moment of negative particles
          "netLa.="       << &netLa <<                  // second moment of unlikesign particles
          "netCh.="       << &netCh <<                  // second moment of unlikesign particles
          "noResnetPi.="  << &noResnetPi <<             // second moments for particle+antiparticle
          "noResnetKa.="  << &noResnetKa <<             // second moment of positive particles
          "noResnetPr.="  << &noResnetPr <<             // second moment of negative particles
          "noResnetLa.="  << &noResnetLa <<             // second moment of unlikesign particles
          "noResnetCh.="  << &noResnetCh <<             // second moment of unlikesign particles
          "\n";  
	} // tree filling
        
      } // ======= end of Centrality loop ======= 
    }// ======= end of momentum loop ======= 
  } // ======= end of eta loop =======
  // 
  // ======================================================================
  // 
  // ======================================================================

     
}
//________________________________________________________________________
void AliAnalysisTaskEbyeIterPID::WeakAndMaterial()
{
  
  //
  // Fill dEdx information for the TPC and also the clean kaon and protons
  //
  if (fUseCouts) std::cout << " Info::marsland: ===== In the WeakAndMaterial ===== " << std::endl;
  //
  // ======================================================================
  // ======================================================================
  //   
  // ========= Efficiency Check Eta momentum and Centrality scan ==========
  AliMCParticle *trackMCgen;
  for (Int_t iTrack = 0; iTrack < fMCEvent->GetNumberOfTracks(); iTrack++) {    // track loop

    // initialize the dummy particle id
    Int_t fElWeak =-100., fPiWeak =-100., fKaWeak =-100., fPrWeak =-100., fDeWeak =-100., fMuWeak =-100.;
    trackMCgen = (AliMCParticle *)fMCEvent->GetTrack(iTrack);
      
    // apply primary vertex and eta cut
    if ( !fMCStack->IsPhysicalPrimary(iTrack) ) continue;
    
    Int_t pdg   = trackMCgen->Particle()->GetPdgCode();
    Int_t iPart = -10;
    if (TMath::Abs(pdg) == 11)          {iPart = 0; fElWeak = iPart;} // select el-
    if (TMath::Abs(pdg) == 211)         {iPart = 1; fPiWeak = iPart;} // select pi+
    if (TMath::Abs(pdg) == 321)         {iPart = 2; fKaWeak = iPart;} // select ka+
    if (TMath::Abs(pdg) == 2212)        {iPart = 3; fPrWeak = iPart;} // select pr+
    if (TMath::Abs(pdg) == 1000010020)  {iPart = 4; fDeWeak = iPart;} // select de
    if (TMath::Abs(pdg) == 13)          {iPart = 5; fMuWeak = iPart;} // select mu-
    if (iPart == -10) continue;
    Float_t fSignWeak = (pdg<0) ? -1:1;
    Float_t fptotWeak = trackMCgen->P();  
    Float_t fpTWeak   = trackMCgen->Pt();
    Float_t fYWeak    = trackMCgen->Y();
    Float_t fEtaWeak  = trackMCgen->Eta();
    Float_t fPhiWeak  = trackMCgen->Phi();    
    
    if(!fTreeSRedirector) return;
    (*fTreeSRedirector)<<"fullacc"<<
    "ptot="      << fptotWeak <<         // mc momentum
    "pT="        << fpTWeak <<         // mc momentum
    "Y="         << fYWeak <<         // mc momentum
    "eta="       << fEtaWeak <<          // mc eta
    "phi="       << fPhiWeak <<          // mc eta
    "cent="      << fCentrality <<     // Centrality
    "sign="      << fSignWeak <<         // sign
    "el="        << fElWeak <<         // sign
    "pi="        << fPiWeak <<         // sign
    "ka="        << fKaWeak <<         // sign
    "pr="        << fPrWeak <<         // sign
    "de="        << fDeWeak <<         // sign
    "mu="        << fMuWeak <<         // sign
    "\n";   
	
	
  } // ======= end of track loop ======= 
  
}
//________________________________________________________________________
void AliAnalysisTaskEbyeIterPID::FillDnchDeta()
{
  
  //
  // Fill dEdx information for the TPC and also the clean kaon and protons
  //
  if (fUseCouts) std::cout << " Info::marsland: ===== In the FillDnchDeta ===== " << std::endl;
  //
  // ======================================================================
  // ======================================================================
  //   
  // ========= Efficiency Check Eta momentum and Centrality scan ==========
  const Int_t netabins = 20;
  Double_t etaDownArray[netabins] ={0.};
  Double_t etaUpArray[netabins]   ={0.};
  for (Int_t i=0; i<netabins; i++){
      etaUpArray[i]=0.1*(i+1);
      etaDownArray[i]=etaUpArray[i]*-1.;
  }
  Double_t centDownArray[9]={0., 5.,  10., 20., 30., 40., 50., 60., 70.};
  Double_t centUpArray[9]  ={5., 10., 20., 30., 40., 50., 60., 70., 80.};
  for (Int_t ieta=0; ieta<netabins; ieta++){
    for (Int_t icent=0; icent<9; icent++){
      
      AliMCParticle *trackMCgen;
      Int_t trCount=0,    elCount=0,    piCount=0,    kaCount=0,    prCount=0;
      for (Int_t iTrack = 0; iTrack < fMCEvent->GetNumberOfTracks(); iTrack++) {    // track loop
	
	// initialize the dummy particle id
	trackMCgen = (AliMCParticle *)fMCEvent->GetTrack(iTrack);
	
	// apply primary vertex and eta cut
	if (!fMCStack->IsPhysicalPrimary(iTrack)) continue;
        // 	if ((trackMCgen->Eta()<fEtaDown) || (trackMCgen->Eta()>fEtaUp)) continue;
	Int_t pdg  = trackMCgen->Particle()->GetPdgCode();
	Double_t etaGen  = trackMCgen->Eta();
	// skip neutral particles
	if ( TMath::Abs(trackMCgen->Charge()) < 0.0001 ) continue; 
	
	Int_t iPart = -10;
	if (TMath::Abs(pdg) == 11)         iPart = 0;  // select el-
	if (TMath::Abs(pdg) == 211)        iPart = 1;  // select pi+
	if (TMath::Abs(pdg) == 321)        iPart = 2;  // select ka+
	if (TMath::Abs(pdg) == 2212)       iPart = 3;  // select pr+
	
	// count first moments
	if ((fCentrality>=centDownArray[icent])
	  &&(fCentrality<centUpArray[icent])
	  &&(etaGen>=etaDownArray[ieta])
	  &&(etaGen<=etaUpArray[ieta])) 
	{ 
	  trCount++;
	  if ( iPart==0   ) elCount++;
	  if ( iPart==1   ) piCount++;
	  if ( iPart==2   ) kaCount++;
	  if ( iPart==3   ) prCount++;
	}   
	
      } // ======= end of track loop ======= 
      
      // fill tree which contains moments
      Double_t etaBin  = (TMath::Abs(etaDownArray[ieta])+etaUpArray[ieta]);
      Double_t centBin = (centDownArray[icent]+centUpArray[icent])/2.;
      if(!fTreeSRedirector) return;
      if ( trCount>0 ){   
	(*fTreeSRedirector)<<"dnchdeta"<<
	"event="    << fEventCountInFile << 
	"centbin="  << centBin <<                 // cent bin
	"etabin="   << etaBin <<                  // eta bin
	"imppar="   << fMCImpactParameter <<      // impact parameter taken from MC event header
	"cent="     << fCentrality <<             // impact parameter taken from MC event header
	"trcount="  << trCount <<                 // number of identified tracks within the given cent and mom range
	"el="       << elCount <<                 // first moment of pions
	"pi="       << piCount <<                 // first moment of pions
	"ka="       << kaCount <<                 // first moment of kaons                 
	"pr="       << prCount <<                 // first moment of protons
	"\n";  
      } // tree filling
      
    } // ======= end of Centrality loop ======= 
  } // ======= end of eta loop =======
  // 
  // ======================================================================
  // 
  // ======================================================================
  
}
//________________________________________________________________________
void AliAnalysisTaskEbyeIterPID::FillTPCdEdxMCEffMatrix()
{
  
  //
  // Fill dEdx information for the TPC and also the clean kaon and protons
  //
  if (fUseCouts) std::cout << " Info::marsland: ===== In the FillTPCdEdxMCEffMatrix ===== " << std::endl;
  AliMCEventHandler *eventHandler = dynamic_cast<AliMCEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  TTree *chain = (TChain*)GetInputData(0);
  if(!chain) { Printf("ERROR: Could not receive input chain"); return; }
  TObjString fileName(chain->GetCurrentFile()->GetName());
  Int_t runNumber = fESD->GetRunNumber();
  //
  // -----------------------------------------------------------------------------------------
  // -----------------------------------------------------------------------------------------
  // Get the real primary vertex
  //   AliGenEventHeader* genHeader = 0x0;
  //   if (mcEvent) genHeader = mcEvent->GenEventHeader(); 
  //   TList* lsth = ((AliGenCocktailEventHeader*)genHeader)->GetHeaders();
  //   TArrayF vtx(3);
  //   TIter next(lsth);
  //   Int_t ivt=0; 
  //   while (genHeader=(AliGenEventHeader*)next()) 
  //   {
  //       if (ivt>0) return;
  //       genHeader->PrimaryVertex(vtx); 
  //       printf("FillTPCdEdxMCEffMatrix First Three vertices: Vtx:%d %f %f %f\n",ivt++,vtx[0],vtx[1],vtx[2]);
  //   }
  //
  // -----------------------------------------------------------------------------------------
  // -----------------------------------------------------------------------------------------
  // ----------------------------   reconstructed MC particles  ------------------------------
  // -----------------------------------------------------------------------------------------
  // loop over tracks
  for(Int_t i = 0; i < fESD->GetNumberOfTracks(); i++) 
  { // track loop  
      AliESDtrack *trackReal = fESD->GetTrack(i); 
      Int_t lab = TMath::Abs(trackReal->GetLabel()); // avoid from negatif labels, they include some garbage
      fEtaMC    = trackReal->Eta();
      // MC track cuts
      if ((fEtaMC<fEtaDown) || (fEtaMC>fEtaUp))          continue; 
      if (!fMCStack->IsPhysicalPrimary(lab))                continue;  // MC primary track check
      // Track cuts from detector
      if (!trackReal -> GetInnerParam())                 continue;
      if (!fESDtrackCuts -> AcceptTrack(trackReal))      continue;  // real track cuts 
      // Extra TPC cuts for better PID
      if (fTightCuts) if (trackReal->GetTPCsignalN()<70) continue;                                          
      if (fTightCuts) if (trackReal->GetLengthInActiveZone(1,3,230, trackReal->GetBz(),0,0)<120) continue; 
      // get track info
      Float_t fpTRec   = trackReal->Pt();
      Float_t fYRec    = trackReal->Y();
      Float_t fEtaRec  = trackReal->Eta();
      Float_t fPhiRec  = trackReal->Phi(); 
      Int_t fCentRec = fhCent->FindBin(fCentrality)-1;
      Int_t fPartID  = -10;
      
      // Efficiency matices for individual particles
      TParticle *trackMC  = fMCStack->Particle(lab); 
      Float_t vZMC        = trackMC->Vz();
      
      Int_t pdg           = trackMC->GetPdgCode();  
      
      if (TMath::Abs(pdg) == 211)   fPartID=0; // select pi
      if (TMath::Abs(pdg) == 321)   fPartID=1; // select ka
      if (TMath::Abs(pdg) == 2212)  fPartID=2; // select pr
      if (fPartID == -10) continue;
      
      // additional TOF requirement 
      if (fIncludeTOF && trackReal->GetInnerParam()->GetP()>0.8){
          Float_t nSigmasPiTOF = fPIDResponse->NumberOfSigmasTOF(trackReal, AliPID::kPion);
          Float_t nSigmasKaTOF = fPIDResponse->NumberOfSigmasTOF(trackReal, AliPID::kKaon);
          Float_t nSigmasPrTOF = fPIDResponse->NumberOfSigmasTOF(trackReal, AliPID::kProton);
          if ( !( 
              ((TMath::Abs(nSigmasPiTOF)<=3) && fPartID==0) ||  
              ((TMath::Abs(nSigmasKaTOF)<=3) && fPartID==1) ||
              ((TMath::Abs(nSigmasPrTOF)<=3) && fPartID==2) 
          ) ) continue;
      }
      
      Double_t xxxRec[5]={Float_t(fPartID),Float_t(fCentRec),fpTRec,fEtaRec,fPhiRec};
      if (pdg>0) fHistPosEffMatrixRec->Fill(xxxRec);
      if (pdg<0) fHistNegEffMatrixRec->Fill(xxxRec);
      
  } // ======= end of track loop =======
  
  // -----------------------------------------------------------------------------------------
  // ----------------------------   MC generated pure MC particles  --------------------------
  // -----------------------------------------------------------------------------------------
  AliMCParticle *trackMCgen;
  for (Int_t iTrack = 0; iTrack < fMCEvent->GetNumberOfTracks(); iTrack++) 
  { // track loop
      trackMCgen = (AliMCParticle *)fMCEvent->GetTrack(iTrack);
      TParticle *trackMC  = fMCStack->Particle(iTrack);  
      Float_t vZMC        = trackMC->Vz();
      if ((trackMCgen->Eta()<fEtaDown) || (trackMCgen->Eta()>fEtaUp)) continue;
      if (!fMCStack->IsPhysicalPrimary(iTrack)) continue;
      
      // get track info
      Float_t fpTGen   = trackMCgen->Pt();
      Float_t fYGen    = trackMCgen->Y();
      Float_t fEtaGen  = trackMCgen->Eta();
      Float_t fPhiGen  = trackMCgen->Phi(); 
      Int_t fCentGen = fhCent->FindBin(fCentrality)-1;
      Int_t fPartID  = -10;
      
      // Efficiency matices for individual particles
      Int_t pdg  = trackMCgen->Particle()->GetPdgCode();    
      if (TMath::Abs(pdg) == 211)   fPartID=0; // select pi
      if (TMath::Abs(pdg) == 321)   fPartID=1; // select ka
      if (TMath::Abs(pdg) == 2212)  fPartID=2; // select pr
      if (fPartID == -10) continue;
      
      Double_t xxxGen[5]={Float_t(fPartID),Float_t(fCentGen),fpTGen,fEtaGen,fPhiGen};
      if (pdg>0) fHistPosEffMatrixGen->Fill(xxxGen);
      if (pdg<0) fHistNegEffMatrixGen->Fill(xxxGen);
  } // ======= end of track loop ======= 
  
}
//________________________________________________________________________
void AliAnalysisTaskEbyeIterPID::FillCleanElectrons()
{

// Fill Clean Electrons from conversion
  if (fUseCouts) std::cout << " Info::marsland: ===== In the FillCleanElectrons ===== " << std::endl;
  if (fPIDResponse) {
    fPIDResponse->GetTPCResponse().SetBetheBlochParameters(1.28778e+00/50., 3.13539e+01, TMath::Exp(-3.16327e+01), 1.87901e+00, 6.41583e+00);
  }
  
  TObjArray* listCrossV0 = fESDtrackCutsV0->GetAcceptedV0s(fESD);  
  Int_t nGoodV0s         = listCrossV0->GetEntries();
  delete listCrossV0;
 
  AliKFParticle::SetField(fESD->GetMagneticField());	
  AliKFVertex primVtx(*(fESD->GetPrimaryVertex()));
  
  Double_t mn[3] = {0,0,0};
  Double_t mp[3] = {0,0,0}; 
  Double_t mm[3] = {0,0,0};
  const Double_t cProtonMass   = TDatabasePDG::Instance()->GetParticle(2212)->Mass();
  const Double_t cPionMass     = TDatabasePDG::Instance()->GetParticle(211)->Mass();
  const Double_t cElectronMass = TDatabasePDG::Instance()->GetParticle(11)->Mass();
  
  // V0 finder for Clean electrons
  for(Int_t iV0MI = 0; iV0MI < nGoodV0s; iV0MI++) {
    AliESDv0 * fV0s = fESD->GetV0(iV0MI);

    Int_t    lOnFlyStatus = 0;
    lOnFlyStatus = fV0s->GetOnFlyStatus();
    if (!lOnFlyStatus) continue;
    
    AliESDtrack* trackPosTest = fESD->GetTrack(fV0s->GetPindex());
    AliESDtrack* trackNegTest = fESD->GetTrack(fV0s->GetNindex());
    
    // my cuts
    if (!fESDtrackCutsCleanSamp->AcceptTrack(trackPosTest)) continue; // To FIX
    if (!fESDtrackCutsCleanSamp->AcceptTrack(trackNegTest)) continue; // To FIX
    if (!trackPosTest->GetInnerParam()) continue;
    if (!trackNegTest->GetInnerParam()) continue;
    
    if (fTightCuts) if (trackPosTest->GetTPCsignalN()<70) continue;
    if (fTightCuts) if (trackNegTest->GetTPCsignalN()<70) continue;
    if (fTightCuts) if (trackPosTest->GetLengthInActiveZone(1,3,230, trackPosTest->GetBz(),0,0)<120) continue;  // to FIX 
    if (fTightCuts) if (trackNegTest->GetLengthInActiveZone(1,3,230, trackNegTest->GetBz(),0,0)<120) continue;  // to FIX 
    
    
    if( trackPosTest->GetSign() >0 && trackNegTest->GetSign() <0){
      fV0s->GetNPxPyPz(mn[0],mn[1],mn[2]); //reconstructed cartesian momentum components of negative daughter
      fV0s->GetPPxPyPz(mp[0],mp[1],mp[2]); //reconstructed cartesian momentum components of positive daughter
    }

    if( trackPosTest->GetSign() <0 && trackNegTest->GetSign() >0){
      fV0s->GetPPxPyPz(mn[0],mn[1],mn[2]); //reconstructed cartesian momentum components of negative daughter
      fV0s->GetNPxPyPz(mp[0],mp[1],mp[2]); //reconstructed cartesian momentum components of positive daughter
    }
	
    fV0s->GetPxPyPz(mm[0],mm[1],mm[2]); //reconstructed cartesian momentum components of mother
      
    TVector3 vecN(mn[0],mn[1],mn[2]);
    TVector3 vecP(mp[0],mp[1],mp[2]);
    TVector3 vecM(mm[0],mm[1],mm[2]);

    Double_t thetaP=0, thetaN=0, alfaEl=0, qtEl=0.;
    // Some if checks for floating point protection
    if ((vecP.Mag() * vecM.Mag())<0.00001) continue;
    if ((vecN.Mag() * vecM.Mag())<0.00001) continue;
    thetaP = acos((vecP * vecM)/(vecP.Mag() * vecM.Mag()));
    thetaN = acos((vecN * vecM)/(vecN.Mag() * vecM.Mag()));
    if ( ((vecP.Mag())*cos(thetaP)+(vecN.Mag())*cos(thetaN)) <0.00001) continue;
    alfaEl = ((vecP.Mag())*cos(thetaP)-(vecN.Mag())*cos(thetaN))/((vecP.Mag())*cos(thetaP)+(vecN.Mag())*cos(thetaN));
    qtEl   = vecP.Mag()*sin(thetaP);

    // fV0s->ChangeMassHypothesis(22);   // ?????

    // V0 particle kinematics
    TLorentzVector posE, negE, photon, posP, negP, posPi, negPi, lambda, antilambda, kaon, sigma0;
    negE .SetXYZM(mn[0],mn[1],mn[2],cElectronMass);
    posE .SetXYZM(mp[0],mp[1],mp[2],cElectronMass);
    negPi.SetXYZM(mn[0],mn[1],mn[2],cPionMass);
    posP .SetXYZM(mp[0],mp[1],mp[2],cProtonMass);
    negP .SetXYZM(mn[0],mn[1],mn[2],cProtonMass);
    posPi.SetXYZM(mp[0],mp[1],mp[2],cPionMass);
    photon=posE+negE;
    kaon=posPi+negPi;
    lambda=posP+negPi;
    antilambda=posPi+negP;

    // Sum cuts
    if (kaon.M()<0.504   && kaon.M()>0.490 ) continue;
    if (lambda.M()>1.113 && lambda.M()<1.118) continue; 
    if (antilambda.M()>1.113 && antilambda.M()<1.118) continue; 
    if (photon.M()>0.005) continue;
    if (TMath::Abs(alfaEl)>0.5) continue;
    if (qtEl > 0.1) continue;
        
    ///////////////// Electron dEdx part///////////////////////         
    Double_t posNTPCSigmaEl = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trackPosTest, AliPID::kElectron)); 
    Double_t negNTPCSigmaEl = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trackNegTest, AliPID::kElectron)); 
        
    Double_t elTPCSignal=0, elptot=0, elEta=0, elSign=0;
    for (Int_t isign = 0; isign<2; isign++){
      if (isign == 0 && (negNTPCSigmaEl<3)) {
        elTPCSignal  = trackPosTest->GetTPCsignal();
        elptot       = trackPosTest->GetInnerParam()->GetP();
        elEta        = trackPosTest->Eta();
        elSign       = trackPosTest->GetSign();
      } else if (isign == 1 && (posNTPCSigmaEl<3)) {
        elTPCSignal  = trackNegTest->GetTPCsignal();
        elptot       = trackNegTest->GetInnerParam()->GetP();
        elEta        = trackNegTest->Eta();
        elSign       = trackNegTest->GetSign();
      }
      
       // Fill the THnSparseF for the inclusive dEdx spectrum
      Double_t weightCleanEl[5]    = {elSign,fCentrality,elEta,elptot, elTPCSignal};
      if (!fMCtrue && fUseThnSparse) fhnCleanEl->Fill(weightCleanEl); 
    } 
  } // end of V0 loop
}
//________________________________________________________________________
void AliAnalysisTaskEbyeIterPID::FillCleanPions()
{

  // Fill Clean Pions from K0s
  if (fUseCouts) std::cout << " Info::marsland: ===== In the FillCleanPions ===== " << std::endl;
  if (fPIDResponse) {
    fPIDResponse->GetTPCResponse().SetBetheBlochParameters(1.28778e+00/50., 3.13539e+01, TMath::Exp(-3.16327e+01), 1.87901e+00, 6.41583e+00);
  }

  TObjArray* listCrossV0   = fESDtrackCutsV0->GetAcceptedV0s(fESD);  
  Int_t nGoodV0s      = listCrossV0->GetEntries();
  delete listCrossV0;
 
  AliKFParticle::SetField(fESD->GetMagneticField());	
  AliKFVertex primVtx(*(fESD->GetPrimaryVertex()));
  
  Double_t mn[3] = {0,0,0};
  Double_t mp[3] = {0,0,0}; 
  Double_t mm[3] = {0,0,0};
  const Double_t cProtonMass=TDatabasePDG::Instance()->GetParticle(2212)->Mass();
  const Double_t cPionMass=TDatabasePDG::Instance()->GetParticle(211)->Mass();
  const Double_t cElectronMass=TDatabasePDG::Instance()->GetParticle(11)->Mass();
  
  //V0 finder for Clean pions!!!!
  for(Int_t iV0MI = 0; iV0MI < nGoodV0s; iV0MI++) {
    AliESDv0 * fV0MIs = fESD->GetV0(iV0MI);

    Int_t    lOnFlyStatus = 0;
    lOnFlyStatus = fV0MIs->GetOnFlyStatus();
    if (!lOnFlyStatus) continue;
    
    AliESDtrack* trackPosTest = fESD->GetTrack(fV0MIs->GetPindex());
    AliESDtrack* trackNegTest = fESD->GetTrack(fV0MIs->GetNindex());
     
    // my cuts
    if (!fESDtrackCutsCleanSamp->AcceptTrack(trackPosTest)) continue; // To FIX
    if (!fESDtrackCutsCleanSamp->AcceptTrack(trackNegTest)) continue; // To FIX
    if (!trackPosTest->GetInnerParam()) continue;
    if (!trackNegTest->GetInnerParam()) continue;
    if (fCentrality>80) continue;
    
    if (fTightCuts) if (trackPosTest->GetTPCsignalN()<70) continue;
    if (fTightCuts) if (trackNegTest->GetTPCsignalN()<70) continue;
    if (fTightCuts) if (trackPosTest->GetLengthInActiveZone(1,3,230, trackPosTest->GetBz(),0,0)<120) continue;  // to FIX 
    if (fTightCuts) if (trackNegTest->GetLengthInActiveZone(1,3,230, trackNegTest->GetBz(),0,0)<120) continue;  // to FIX 
    
    if( trackPosTest->GetSign() >0 && trackNegTest->GetSign() <0){
      fV0MIs->GetNPxPyPz(mn[0],mn[1],mn[2]); //reconstructed cartesian momentum components of negative daughter
      fV0MIs->GetPPxPyPz(mp[0],mp[1],mp[2]); //reconstructed cartesian momentum components of positive daughter
    }

    if( trackPosTest->GetSign() <0 && trackNegTest->GetSign() >0){
      fV0MIs->GetPPxPyPz(mn[0],mn[1],mn[2]); //reconstructed cartesian momentum components of negative daughter
      fV0MIs->GetNPxPyPz(mp[0],mp[1],mp[2]); //reconstructed cartesian momentum components of positive daughter
    }
	
    fV0MIs->GetPxPyPz(mm[0],mm[1],mm[2]); //reconstructed cartesian momentum components of mother
    
    TVector3 vecN(mn[0],mn[1],mn[2]);
    TVector3 vecP(mp[0],mp[1],mp[2]);
    TVector3 vecM(mm[0],mm[1],mm[2]);

    if ((vecP.Mag() * vecM.Mag())<0.00001) continue;
    if ((vecN.Mag() * vecM.Mag())<0.00001) continue;
    Double_t thetaP  = acos((vecP * vecM)/(vecP.Mag() * vecM.Mag()));
    Double_t thetaN  = acos((vecN * vecM)/(vecN.Mag() * vecM.Mag()));
    if ( ((vecP.Mag())*cos(thetaP)+(vecN.Mag())*cos(thetaN)) <0.00001) continue;
    fAlfa = ((vecP.Mag())*cos(thetaP)-(vecN.Mag())*cos(thetaN))/((vecP.Mag())*cos(thetaP)+(vecN.Mag())*cos(thetaN));
    fQt   = vecP.Mag()*sin(thetaP);
    
    fHistArmPod->Fill(fAlfa,fQt);
    
    // main armentoros podolanki cuts
    if (TMath::Abs(fAlfa)>0.9) continue;
    if (fQt<0.06) continue;
    if (fQt >0.22) continue;
    SelectCleanSamplesFromV0s(fV0MIs,trackPosTest,trackNegTest);
    //     if (fQt>0.11 && fQt<0.16) continue;
    //     if (TMath::Abs(fAlfa)<0.5 && fQt<0.15) continue;
    
    // fV0MIs->ChangeMassHypothesis(310); // ?????
    
    TLorentzVector posE, negE, photon, posP, negP, posPi, negPi, lambda, antiLambda, kaon, posProton, k0sProton;   
    negE.SetXYZM(mn[0],mn[1],mn[2],cElectronMass);
    posE.SetXYZM(mp[0],mp[1],mp[2],cElectronMass);
    negPi.SetXYZM(mn[0],mn[1],mn[2],cPionMass);
    posPi.SetXYZM(mp[0],mp[1],mp[2],cPionMass);
    negP.SetXYZM(mn[0],mn[1],mn[2],cProtonMass);
    posP.SetXYZM(mp[0],mp[1],mp[2],cProtonMass);
    kaon=posPi+negPi; 	
    photon=posE+negE;    
    lambda=posP+negPi;
    antiLambda=posPi+negP;
                
    // Fill the tree
    Double_t posNTPCSigmaPi = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trackPosTest, AliPID::kPion)); 
    Double_t negNTPCSigmaPi = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trackNegTest, AliPID::kPion)); 
    Double_t posNTPCSigmaPr = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trackPosTest, AliPID::kProton)); 
    Double_t negNTPCSigmaPr = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trackNegTest, AliPID::kProton)); 
    for (Int_t isign = 0; isign<2; isign++){
      if (isign == 0 && (negNTPCSigmaPi<3 || negNTPCSigmaPr<3)) {
        fArmPodTPCSignal  = trackPosTest->GetTPCsignal();
        fArmPodptot       = trackPosTest->GetInnerParam()->GetP();
        fArmPodEta        = trackPosTest->Eta();
        fPiNSigmasTOF = fPIDResponse->NumberOfSigmasTOF(trackPosTest, AliPID::kPion,  fPIDResponse->GetTOFResponse().GetTimeZero());
        fPrNSigmasTOF = fPIDResponse->NumberOfSigmasTOF(trackPosTest, AliPID::kProton,fPIDResponse->GetTOFResponse().GetTimeZero());
      } else if ( isign == 1 && (posNTPCSigmaPi<3 || posNTPCSigmaPr<3)){
        fArmPodTPCSignal  = trackNegTest->GetTPCsignal();
        fArmPodptot       = trackNegTest->GetInnerParam()->GetP();
        fArmPodEta        = trackNegTest->Eta();
        fPiNSigmasTOF = fPIDResponse->NumberOfSigmasTOF(trackNegTest, AliPID::kPion,  fPIDResponse->GetTOFResponse().GetTimeZero());
        fPrNSigmasTOF = fPIDResponse->NumberOfSigmasTOF(trackNegTest, AliPID::kProton,fPIDResponse->GetTOFResponse().GetTimeZero());
      }
      fArmPodCentrality = fCentrality;
      if (fArmPodTPCSignal<30 || fArmPodTPCSignal>200) continue;
      if (fArmPodptot<fMomDown || fArmPodptot>fMomUp) continue;
      if(fFillArmPodTree) fArmPodTree->Fill(); 
      //       if (fArmPodTPCSignal<30 || fArmPodTPCSignal>200) continue;
      //       if (fArmPodptot<fMomDown || fArmPodptot>fMomUp) continue;
      //       // Fill the THnSparseF for the inclusive dEdx spectrum
      //       if ( (TMath::Abs(fPiNSigmasTOF)<3) || (TMath::Abs(fPrNSigmasTOF)<3)) {
      //           fArmPodCentrality = fCentrality;
      //           if(fFillArmPodTree) fArmPodTree->Fill(); 
      //       }   
    }
    
  } // end of V0 loop

}
//________________________________________________________________________
void  AliAnalysisTaskEbyeIterPID::SelectCleanSamplesFromV0s(AliESDv0 *v0, AliESDtrack *track0, AliESDtrack *track1)
{
  //
  // SetAliases and Metadata for the V0 trees
  //
  AliKFParticle kfparticle; //
  AliAnalysisTaskFilteredTree filteredV0;
  Int_t type=filteredV0.GetKFParticle(v0,fESD,kfparticle);
  //
  TDatabasePDG *pdg = TDatabasePDG::Instance();
  Double_t massLambda = pdg->GetParticle("Lambda0")->Mass();
  Double_t massK0 = pdg->GetParticle("K0")->Mass();
  Double_t massPion = pdg->GetParticle("pi+")->Mass();
  Double_t massProton = pdg->GetParticle("proton")->Mass();
  const Double_t livetimeK0=2.684341668932;  // livetime in cm (surpisely missing info in PDG - see root forum)
  const Double_t livetimeLambda=7.8875395;  // livetime in cm (missing info in PDG - see root forum)
  fHasTrack0FirstITSlayer = track0->HasPointOnITSLayer(0);
  fHasTrack1FirstITSlayer = track1->HasPointOnITSLayer(0);
  fHasV0FirstITSlayer = (fHasTrack0FirstITSlayer||fHasTrack1FirstITSlayer);
  
  //
  //
  Double_t v0Rr = v0->GetRr();  // rec position of the vertex CKBrev
  Double_t v0P  = v0->P();      // TMath::Sqrt(Px()*Px()+Py()*Py()+Pz()*Pz())
  
  //   tree->SetAlias("livetimeLikeK0",TString::Format("exp(-v0.fRr/(sqrt((v0.P()/%f)^2+1)*%f))",massK0, livetimeK0)); 
  //   tree->SetAlias("livetimeLikeLambda",TString::Format("exp(-v0.fRr/(sqrt((v0.P()/%f)^2+1)*%f))",massLambda,livetimeLambda)); 
  //   tree->SetAlias("livetimeLikeGamma","v0.fRr/80");
  //   tree->SetAlias("livetimeLikeBkg","v0.fRr/80"); 
  Double_t livetimeLikeK0 = TMath::Exp(-v0Rr/(TMath::Sqrt((v0P/massK0)*(v0P/massK0)+1)*livetimeK0));
  Double_t livetimeLikeLambda = TMath::Exp(-v0Rr/(TMath::Sqrt((v0P/massLambda)*(v0P/massLambda)+1)*livetimeLambda));
  Double_t livetimeLikeGamma = v0Rr/80.;
  Double_t livetimeLikeBkg   = v0Rr/80.;

  // delta of mass
  Double_t K0Delta = v0->GetEffMass(2,2)-massK0;        //   tree->SetAlias("K0Delta","(v0.GetEffMass(2,2)-massK0)");
  Double_t LDelta  = v0->GetEffMass(4,2)-massLambda;    //   tree->SetAlias("LDelta","(v0.GetEffMass(4,2)-massLambda)");
  Double_t ALDelta = v0->GetEffMass(2,4)-massLambda;    //   tree->SetAlias("ALDelta","(v0.GetEffMass(2,4)-massLambda)");
  Double_t EDelta  = v0->GetEffMass(0,0);               //   tree->SetAlias("EDelta","(v0.GetEffMass(0,0))");

  // pull of the mass
  if (v0->GetKFInfo(2,2,1)==0. || v0->GetKFInfo(4,2,1)==0. || v0->GetKFInfo(2,4,1)==0. || v0->GetKFInfo(0,0,1)==0.) return;
  Double_t K0Pull = (v0->GetEffMass(2,2)-massK0)/v0->GetKFInfo(2,2,1);        //   tree->SetAlias("K0Pull","(v0.GetEffMass(2,2)-massK0)/v0.GetKFInfo(2,2,1)");
  Double_t LPull  = (v0->GetEffMass(4,2)-massLambda)/v0->GetKFInfo(4,2,1);    //   tree->SetAlias("LPull","(v0.GetEffMass(4,2)-massLambda)/v0.GetKFInfo(4,2,1)");
  Double_t ALPull = (v0->GetEffMass(2,4)-massLambda)/v0->GetKFInfo(2,4,1);    //   tree->SetAlias("ALPull","(v0.GetEffMass(2,4)-massLambda)/v0.GetKFInfo(2,4,1)");
  Double_t EPull  = EDelta/v0->GetKFInfo(0,0,1);                              //   tree->SetAlias("EPull","EDelta/v0.GetKFInfo(0,0,1)");

  // effective pull of the mass - (empirical values from fits)
  //   tree->SetAlias("K0PullEff","K0Delta/sqrt((3.63321e-03)**2+(5.68795e-04*v0.Pt())**2)");
  //   tree->SetAlias("LPullEff","LDelta/sqrt((1.5e-03)**2+(1.8e-04*v0.Pt())**2)");
  //   tree->SetAlias("ALPullEff","ALDelta/sqrt((1.5e-03)**2+(1.8e-04*v0.Pt())**2)");
  //   tree->SetAlias("EPullEff","v0.GetEffMass(0,0)/sqrt((5e-03)**2+(1.e-04*v0.Pt())**2)");
  Double_t K0PullEff = K0Delta/TMath::Sqrt((3.63321e-03)*(3.63321e-03)+(5.68795e-04*v0->Pt())*(5.68795e-04*v0->Pt()));
  Double_t LPullEff  = LDelta/TMath::Sqrt((1.5e-03)*(1.5e-03)+(1.8e-04*v0->Pt())*(1.8e-04*v0->Pt()));
  Double_t ALPullEff = ALDelta/TMath::Sqrt((1.5e-03)*(1.5e-03)+(1.8e-04*v0->Pt())*(1.8e-04*v0->Pt()));
  Double_t EPullEff  = v0->GetEffMass(0,0)/TMath::Sqrt((5e-03)*(5e-03)+(1.e-04*v0->Pt())*(1.e-04*v0->Pt()));

  // 
  //   tree->SetAlias("dEdx0DProton","AliExternalTrackParam::BetheBlochAleph(track0.fIp.P()/massProton)");
  //   tree->SetAlias("dEdx1DProton","AliExternalTrackParam::BetheBlochAleph(track1.fIp.P()/massProton)");
  //   tree->SetAlias("dEdx0DPion","AliExternalTrackParam::BetheBlochAleph(track0.fIp.P()/massPion)");
  //   tree->SetAlias("dEdx1DPion","AliExternalTrackParam::BetheBlochAleph(track1.fIp.P()/massPion)");
  Double_t dEdx0DProton = AliExternalTrackParam::BetheBlochAleph(track0->GetInnerParam()->GetP()/massProton);
  Double_t dEdx1DProton = AliExternalTrackParam::BetheBlochAleph(track1->GetInnerParam()->GetP()/massProton);
  Double_t dEdx0DPion   = AliExternalTrackParam::BetheBlochAleph(track0->GetInnerParam()->GetP()/massPion);
  Double_t dEdx1DPion   = AliExternalTrackParam::BetheBlochAleph(track1->GetInnerParam()->GetP()/massPion);

  //   tree->SetAlias("K0Like0","exp(-K0Pull^2)*livetimeLikeK0");
  //   tree->SetAlias("LLike0","exp(-LPull^2)*livetimeLikeLambda");
  //   tree->SetAlias("ALLike0","exp(-ALPull^2)*livetimeLikeLambda");
  //   tree->SetAlias("ELike0","exp(-abs(EPull)*0.2)*livetimeLikeGamma");
  //   tree->SetAlias("V0Like","exp(-acos(v0.fPointAngle)*v0.fRr/0.36)*exp(-sqrt(kf.GetChi2())/0.5)");
  Double_t K0Like0 = TMath::Exp(-K0Pull*K0Pull)*livetimeLikeK0;
  Double_t LLike0  = TMath::Exp(-LPull*LPull)*livetimeLikeLambda;
  Double_t ALLike0 = TMath::Exp(-ALPull*ALPull)*livetimeLikeLambda;
  Double_t ELike0  = TMath::Exp(-abs(EPull)*0.2)*livetimeLikeGamma;
  Double_t V0Like  = TMath::Exp(-TMath::ACos(v0->GetV0CosineOfPointingAngle())*v0Rr/0.36)*TMath::Exp(-TMath::Sqrt(kfparticle.GetChi2())/0.5);
  
  
  //   tree->SetAlias("BkgLike","0.000005*ntracks");  // backround coeefecint  to be fitted - depends on other cuts 
  Int_t ntracks = fESD->GetNumberOfTracks();
  Double_t BkgLike = 0.000005*ntracks;    // backround coeefecint  to be fitted - depends on other cuts 

  //   tree->SetAlias("ELike","(V0Like*ELike0)/(V0Like*(K0Like0+LLike0+ALLike0+ELike0)+BkgLike)");
  Double_t ELike = (V0Like*ELike0)/(V0Like*(K0Like0+LLike0+ALLike0+ELike0)+BkgLike);
  //   tree->SetAlias("K0Like","K0Like0/(K0Like0+LLike0+ALLike0+ELike0+BkgLike)");
  Double_t K0Like = K0Like0/(K0Like0+LLike0+ALLike0+ELike0+BkgLike);
  //   tree->SetAlias("LLike","LLike0/(K0Like0+LLike0+ALLike0+ELike0+BkgLike)");
  Double_t LLike = LLike0/(K0Like0+LLike0+ALLike0+ELike0+BkgLike);
  //   tree->SetAlias("ALLike","ALLike0/(K0Like0+LLike0+ALLike0+ELike0+BkgLike)");
  Double_t ALLike = ALLike0/(K0Like0+LLike0+ALLike0+ELike0+BkgLike);

  
  Double_t tr0NTPCSigmaPi = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track0, AliPID::kPion)); 
  Double_t tr1NTPCSigmaPi = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track1, AliPID::kPion)); 
  Double_t tr0NTPCSigmaPr = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track0, AliPID::kProton)); 
  Double_t tr1NTPCSigmaPr = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track1, AliPID::kProton)); 
 
  //   chainV0->SetAlias("cleanPion0FromK0","K0Like0>0.05&&K0Like0>(LLike0+ALLike0+ELike0)*3&&abs(K0Delta)<0.006&&V0Like>0.1&&abs(track1.fTPCsignal/dEdx1DPion-50)<8&&v0.PtArmV0()>0.06");
  fCleanPion0FromK0 = (K0Like0>0.05) && (K0Like0>(LLike0+ALLike0+ELike0)*3) && (TMath::Abs(K0Delta)<0.006) && (V0Like>0.1) && tr1NTPCSigmaPi<2 && (v0->PtArmV0()>0.06);
  //   chainV0->SetAlias("cleanPion1FromK0","K0Like0>0.05&&K0Like0>(LLike0+ALLike0+ELike0)*3&&abs(K0Delta)<0.006&&V0Like>0.1&&abs(track0.fTPCsignal/dEdx0DPion-50)<8&&v0.PtArmV0()>0.06");
  fCleanPion1FromK0 = (K0Like0>0.05) && (K0Like0>(LLike0+ALLike0+ELike0)*3) && (TMath::Abs(K0Delta)<0.006) && (V0Like>0.1) && tr0NTPCSigmaPi<2 && (v0->PtArmV0()>0.06);
  //   chainV0->SetAlias("cleanPion0FromLambda","ALLike>0.05&&ALLike0>(K0Like0+LLike0+ELike0)*3&&abs(ALDelta)<0.006&&V0Like>0.1&&abs(track1.fTPCsignal/dEdx1DProton-50)<8");
  fCleanPion0FromLambda = (ALLike>0.05) && (ALLike0>(K0Like0+LLike0+ELike0)*3) && (TMath::Abs(ALDelta)<0.006) && (V0Like>0.1) && tr1NTPCSigmaPr<2;
  //   chainV0->SetAlias("cleanPion1FromLambda","LLike>0.05&&LLike0>(K0Like0+ALLike0+ELike0)*3&&abs(LDelta)<0.006&&V0Like>0.1&&abs(track0.fTPCsignal/dEdx0DProton-50)<8");
  fCleanPion1FromLambda = (LLike>0.05) && (LLike0>(K0Like0+ALLike0+ELike0)*3) && (TMath::Abs(LDelta)<0.006) && (V0Like>0.1) && tr0NTPCSigmaPr<2;
  //   chainV0->SetAlias("cleanProton0FromLambda","LLike>0.05&&LLike0>(K0Like0+ALLike0+ELike0)*3&&abs(LDelta)<0.001&&V0Like>0.1&&abs(track1.fTPCsignal/dEdx1DPion-50)<8");
  fCleanProton0FromLambda = (LLike>0.05) && (LLike0>(K0Like0+ALLike0+ELike0)*3) && (TMath::Abs(LDelta)<0.006) && (V0Like>0.1) && tr1NTPCSigmaPi<2;
  //   chainV0->SetAlias("cleanProton1FromLambda","ALLike>0.05&&ALLike0>(K0Like0+LLike0+ELike0)*3&&abs(ALDelta)<0.001&&V0Like>0.1&&abs(track0.fTPCsignal/dEdx0DPion-50)<8");
  fCleanProton1FromLambda = (ALLike>0.05) && (ALLike0>(K0Like0+LLike0+ELike0)*3) && (TMath::Abs(ALDelta)<0.006) && (V0Like>0.1) && tr0NTPCSigmaPi<2;
  
  fCleanPionsFromK0 =  (fCleanPion0FromK0 || fCleanPion1FromK0);
  
  //    //   chainV0->SetAlias("cleanPion0FromK0","K0Like0>0.05&&K0Like0>(LLike0+ALLike0+ELike0)*3&&abs(K0Delta)<0.006&&V0Like>0.1&&abs(track1.fTPCsignal/dEdx1DPion-50)<8&&v0.PtArmV0()>0.06");
  //   fCleanPion0FromK0 = (K0Like0>0.05) && (K0Like0>(LLike0+ALLike0+ELike0)*3) && (TMath::Abs(K0Delta)<0.006) && (V0Like>0.1) && (TMath::Abs(track1->GetTPCsignal()/dEdx1DPion-50)<8) && (v0->PtArmV0()>0.06);
  //   //   chainV0->SetAlias("cleanPion1FromK0","K0Like0>0.05&&K0Like0>(LLike0+ALLike0+ELike0)*3&&abs(K0Delta)<0.006&&V0Like>0.1&&abs(track0.fTPCsignal/dEdx0DPion-50)<8&&v0.PtArmV0()>0.06");
  //   fCleanPion1FromK0 = (K0Like0>0.05) && (K0Like0>(LLike0+ALLike0+ELike0)*3) && (TMath::Abs(K0Delta)<0.006) && (V0Like>0.1) && (TMath::Abs(track0->GetTPCsignal()/dEdx0DPion-50)<8) && (v0->PtArmV0()>0.06);
  //   //   chainV0->SetAlias("cleanPion0FromLambda","ALLike>0.05&&ALLike0>(K0Like0+LLike0+ELike0)*3&&abs(ALDelta)<0.006&&V0Like>0.1&&abs(track1.fTPCsignal/dEdx1DProton-50)<8");
  //   fCleanPion0FromLambda = (ALLike>0.05) && (ALLike0>(K0Like0+LLike0+ELike0)*3) && (TMath::Abs(ALDelta)<0.006) && (V0Like>0.1) && TMath::Abs(track1->GetTPCsignal()/dEdx1DProton-50)<8;
  //   //   chainV0->SetAlias("cleanPion1FromLambda","LLike>0.05&&LLike0>(K0Like0+ALLike0+ELike0)*3&&abs(LDelta)<0.006&&V0Like>0.1&&abs(track0.fTPCsignal/dEdx0DProton-50)<8");
  //   fCleanPion1FromLambda = (LLike>0.05) && (LLike0>(K0Like0+ALLike0+ELike0)*3) && (TMath::Abs(LDelta)<0.006) && (V0Like>0.1) && TMath::Abs(track0->GetTPCsignal()/dEdx0DProton-50)<8;
  //   //   chainV0->SetAlias("cleanProton0FromLambda","LLike>0.05&&LLike0>(K0Like0+ALLike0+ELike0)*3&&abs(LDelta)<0.001&&V0Like>0.1&&abs(track1.fTPCsignal/dEdx1DPion-50)<8");
  //   fCleanProton0FromLambda = (LLike>0.05) && (LLike0>(K0Like0+ALLike0+ELike0)*3) && (TMath::Abs(LDelta)<0.006) && (V0Like>0.1) && TMath::Abs(track1->GetTPCsignal()/dEdx1DPion-50)<8;
  //   //   chainV0->SetAlias("cleanProton1FromLambda","ALLike>0.05&&ALLike0>(K0Like0+LLike0+ELike0)*3&&abs(ALDelta)<0.001&&V0Like>0.1&&abs(track0.fTPCsignal/dEdx0DPion-50)<8");
  //   fCleanProton1FromLambda = (ALLike>0.05) && (ALLike0>(K0Like0+LLike0+ELike0)*3) && (TMath::Abs(ALDelta)<0.006) && (V0Like>0.1) && TMath::Abs(track0->GetTPCsignal()/dEdx0DPion-50)<8;
  //  
  
  
  //   // V0 - cuts - for PID 
  //   tree->SetAlias("cutDist","sqrt((track0.fIp.fP[0]-track1.fIp.fP[0])**2+(track0.fIp.fP[1]-track1.fIp.fP[1])**2)>3");
  //   tree->SetAlias("cutLong","track0.GetTPCClusterInfo(3,1,0)-5*abs(track0.fP[4])>130&&track1.GetTPCClusterInfo(3,1,0)>130-5*abs(track0.fP[4])");
  //   tree->SetAlias("cutPID","track0.fTPCsignal>0&&track1.fTPCsignal>0");
  //   tree->SetAlias("cutResol","sqrt(track0.fC[14]/track0.fP[4])<0.15&&sqrt(track1.fC[14]/track1.fP[4])<0.15");
  //   tree->SetAlias("cutV0","cutPID&&cutLong&&cutResol"); 
  //   //  
  //   tree->SetAlias("K0PullBkg","min(min(abs(LPull),abs(ALPull)),abs(EPull))+0");
  //   tree->SetAlias("LambdaPullBkg","min(min(abs(K0Pull),abs(ALPull)),abs(EPull)+0)");
  //   tree->SetAlias("ALambdaPullBkg","min(min(abs(K0Pull),abs(LPull)),abs(EPull)+0)");
  //   tree->SetAlias("EPullBkg","min(min(abs(K0Pull),abs(LPull)),abs(ALPull)+0)");
  //   //
  //   tree->SetAlias("K0Selected",      "abs(K0Pull)<3. &&abs(K0PullEff)<3.  && abs(LPull)>3  && abs(ALPull)>3 &&v0.PtArmV0()>0.11"); 
  //   tree->SetAlias("LambdaSelected",  "abs(LPull)<3.  &&abs(LPullEff)<3.   && abs(K0Pull)>3 && abs(EPull)>3  && abs(EDelta)>0.05");  
  //   tree->SetAlias("ALambdaSelected", "abs(ALPull)<3. &&abs(ALPullEff)<3   && abs(K0Pull)>3 && abs(EPull)>3  &&abs(EDelta)>0.05");
  //   tree->SetAlias("GammaSelected", "abs(EPull)<3     && abs(K0Pull)>3 && abs(LPull)>3 && abs(ALPull)>3");
  //   tree->SetAlias("BkgLike","0.000005*ntracks");  // backround coeefecint  to be fitted - depends on other cuts 
  //   //
  //   tree->SetAlias("ELike","(V0Like*ELike0)/(V0Like*(K0Like0+LLike0+ALLike0+ELike0)+BkgLike)");
  //   tree->SetAlias("K0Like","K0Like0/(K0Like0+LLike0+ALLike0+ELike0+BkgLike)");
  //   tree->SetAlias("LLike","LLike0/(K0Like0+LLike0+ALLike0+ELike0+BkgLike)");
  //   tree->SetAlias("ALLike","ALLike0/(K0Like0+LLike0+ALLike0+ELike0+BkgLike)");
  //   //
  //   tree->SetAlias("K0PIDPull","(abs(track0.fTPCsignal/dEdx0DPion-50)+abs(track1.fTPCsignal/dEdx1DPion-50))/5.");
  //   tree->SetAlias("mpt","1/v0.Pt()");                 // 
  //   tree->SetAlias("tglV0","v0.Pz()/v0.Pt()");                 // 
  //   tree->SetAlias("alphaV0","atan2(v0.Py(),v0.Px()+0)");
  //   tree->SetAlias("dalphaV0","alphaV0-((int(36+9*(alphaV0/pi))-36)*pi/9.)");

}
//________________________________________________________________________
Bool_t  AliAnalysisTaskEbyeIterPID::ApplyDCAcutIfNoITSPixel(AliESDtrack *track)
{
    
    //     chain->SetAlias("ITS01","(Tracks[].HasPointOnITSLayer(0)||Tracks[].HasPointOnITSLayer(1))");
    //     chain->SetAlias("isPrimPtDep","abs(Tracks[].fD)<0.0182+0.0350/(Tracks[].Pt()^1.01)");
    //     chain->SetAlias("isPrimPtDep2","abs(Tracks[].fD/2)<0.0182+0.0350/(Tracks[].Pt()^1.01)");
    //     chain->SetAlias("isPrim2","sqrt(Tracks[].fD**2/Tracks[].fCdd+Tracks[].fZ**2/Tracks[].fCzz+0)<2");
    //     chain->SetAlias("isPrim5","sqrt(Tracks[].fD**2/Tracks[].fCdd+Tracks[].fZ**2/Tracks[].fCzz+0)<5");
    //     chain->SetAlias("IsPrimCA","((isPrim2&&Tracks[].fITSncls>2)||(isPrim5&&ITS01))");
    
    Float_t p[2],cov[3]; 
    track->GetImpactParameters(p,cov); // p[0]=fD; p[1]=fZ; cov[0]=fCdd; cov[1]=fCdz; cov[2]=fCzz;
    Bool_t isFirstITSlayer  = track->HasPointOnITSLayer(0);
    Bool_t isSecondITSlayer = track->HasPointOnITSLayer(1);
    
    fIsITSpixel01     = (isFirstITSlayer || isSecondITSlayer);
    fnITSclusters = track->GetNumberOfITSClusters();
    
    if (!cov[0] || !cov[2]) {
        return kFALSE;
    } else {
        fPrimRestriction = TMath::Sqrt((p[0]*p[0])/cov[0] + (p[1]*p[1])/cov[2]);
        return (fPrimRestriction<2 && fnITSclusters>2) || (fPrimRestriction<5 && fIsITSpixel01);
    }
            
}
//________________________________________________________________________
void AliAnalysisTaskEbyeIterPID::BinLogAxis(TH1 *h) 
{
  //
  // Method for the correct logarithmic binning of histograms
  //
  if (fUseCouts) std::cout << " Info::marsland: ===== In the BinLogAxis ===== " << std::endl;
  TAxis *axis       = h->GetXaxis();
  Int_t bins        = axis->GetNbins();

  Double_t from     = axis->GetXmin();
  Double_t to       = axis->GetXmax();
  Double_t *newBins = new Double_t[bins + 1];
   
  newBins[0] = from;
  Double_t factor = pow(to/from, 1./bins);
  
  for (Int_t i = 1; i <= bins; i++) {
    newBins[i] = factor * newBins[i-1];
  }
  axis->Set(bins, newBins);
  delete [] newBins;
  
}
//________________________________________________________________________
Int_t AliAnalysisTaskEbyeIterPID::CountEmptyEvents(Int_t counterBin)
{

  //
  // count Empty Events
  //

  TTree *chain = (TChain*)GetInputData(0);
  if(!chain) { Printf("ERROR: Could not receive input chain"); return 0; }
  
  Int_t emptyCount=0;  
  for (Int_t i=0;i<fESD->GetNumberOfTracks();++i) {   // Track loop
    AliESDtrack *track = fESD->GetTrack(i);
    if (!track->GetInnerParam()) continue; 
    Float_t momtrack = track->GetInnerParam()->GetP();
    if (momtrack<fMomDown || momtrack>fMomUp) continue;
    if (!fESDtrackCuts->AcceptTrack(track)) continue;
    if (track->GetTPCsignalN()<60) continue; 
    if (track->GetTPCsignal()>0) emptyCount++;  
  }
  
  // check if the event is empty
  if (emptyCount<1) { 
    fHistEmptyEvent->Fill(counterBin);
    std::cout << " Info::marsland: Empty event in " << chain->GetCurrentFile()->GetName() << std::endl;
  }
  if (fUseCouts) std::cout << " Info::marsland: ====== EVENT IS COOL GO AHEAD ======= " << std::endl;
  return emptyCount;

}
//________________________________________________________________________
void AliAnalysisTaskEbyeIterPID::Terminate(Option_t *) 
{
  std::cout << " Info::marsland: ===== In the Terminate ===== " << std::endl;
  // Draw result to the screen
  // Called once at the end of the query
  
}
