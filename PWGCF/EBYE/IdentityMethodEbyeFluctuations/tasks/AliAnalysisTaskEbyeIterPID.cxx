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
#include "TH1F.h"
#include "THn.h"
#include "TList.h"
#include "TMath.h"
#include "TMatrixF.h"
#include "TVectorF.h"
#include "TCanvas.h"
#include "TObjArray.h"
#include "TF1.h"
#include "TFile.h"
#include "TDatabasePDG.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliKFVertex.h"
#include "AliKFParticle.h"
#include "AliCollisionGeometry.h"
#include "AliGenEventHeader.h"
#include "AliGenHijingEventHeader.h"
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
#include "AliAnalysisTaskEbyeIterPID.h"

#include "iostream"
using namespace std;


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
fhEta(0),
fhCent(0),
fhPtot(0),
fhndEdx(0),
fhnExpected(0),
fhnCleanEl(0),
fhnCleanKa(0),
fhnCleanDe(0),
fEtaDown(0),
fEtaUp(0),
fnEtaBins(0),
fMCtrue(kFALSE),
fEffMatrix(kFALSE),
fdEdxCheck(kFALSE),
fCleanSamplesOnly(kFALSE),
fTightCuts(kFALSE),
fIncludeITS(kTRUE),
fFillBayes(kFALSE),
fFillCuts(kFALSE),
fFillDeDxTree(kTRUE),
fRunFastSimulation(kFALSE),
fFillArmPodTree(kTRUE),
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
fSign(0),
fTPCShared(0),
fNcl(0),
fnEtaWinBinsMC(-100),
fnMomBinsMC(-100),
fnCentBinsMC(-100), 
fnCentbinsData(10),
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
fhEta(0),
fhCent(0),
fhPtot(0),
fhndEdx(0),
fhnExpected(0),
fhnCleanEl(0),
fhnCleanKa(0),
fhnCleanDe(0),
fEtaDown(0),
fEtaUp(0),
fnEtaBins(0),
fMCtrue(kFALSE),
fEffMatrix(kFALSE),
fdEdxCheck(kFALSE),
fCleanSamplesOnly(kFALSE),
fTightCuts(kFALSE),
fIncludeITS(kTRUE),
fFillBayes(kFALSE),
fFillCuts(kFALSE),
fFillDeDxTree(kTRUE),
fRunFastSimulation(kFALSE),
fFillArmPodTree(kTRUE),
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
fSign(0),
fTPCShared(0),
fNcl(0),
fnEtaWinBinsMC(-100),
fnMomBinsMC(-100),
fnCentBinsMC(-100), 
fnCentbinsData(10),
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
  .L /u/marsland/PHD/macros/marsland_EbyeRatios/AliAnalysisTaskEbyeIterPID.cxx++
  .L /hera/alice/marsland/train/trunk/marsland_EbyeRatios/AliAnalysisTaskEbyeIterPID.cxx++
  */
  cout << "===================================================================================="<< endl;
  cout << "===================================================================================="<< endl;
  cout << "===================================================================================="<< endl;
  cout << "***************** CONSTRUCTOR CALLED: AliAnalysisTaskEbyeItPID  ******************"<< endl;
  cout << "===================================================================================="<< endl;
  cout << "===================================================================================="<< endl;
  cout << "===================================================================================="<< endl;
  // ==========================================
  //   
  // ==========================================
  // Create histograms for binning of eta p and cent
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
  // ==========================================

}
//________________________________________________________________________
AliAnalysisTaskEbyeIterPID::~AliAnalysisTaskEbyeIterPID() {
  
  //
  // Destructor
  //
  cout << " ===== In the Destructor ===== " << endl;
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
  cout << " ===== In the Initialize ===== " << endl;
  if (fRunFastSimulation) { cout << " !!! We are running fast simulation return !!! " << endl; return; }
  AliInfoClass("Creating track cuts for ITS+TPC (2010 definition).");
  // fESDtrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE,1);
  
  fESDtrackCuts = new AliESDtrackCuts;
  fESDtrackCuts -> SetEtaRange(fEtaDown,fEtaUp);
  fESDtrackCuts -> SetPtRange(0.1,1e10);   
  //
  // ------------------------------------------------
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
  //
  // ------------------------------------------------
  // Not used in the systematic checks
  // ------------------------------------------------
  if (fIncludeITS){     // Reason for the empty events
    fESDtrackCuts->SetRequireITSRefit(kTRUE);                                              
    fESDtrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny); // Reason for the structure in phi
  }
  fESDtrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);  
  fESDtrackCuts->SetMaxChi2TPCConstrainedGlobal(36);
  fESDtrackCuts->SetAcceptKinkDaughters(kFALSE);
  fESDtrackCuts->SetRequireTPCRefit(kTRUE);
  fESDtrackCuts->SetMaxFractionSharedTPCClusters(0.4);    // ?? FROM MARIAN 
  // ------------------------------------------------
  //
  // ------------------------------------------------
  // Vertex restrictions
  fESDtrackCuts->SetMaxDCAToVertexZ(2);       
  fESDtrackCuts->SetMaxChi2PerClusterITS(36); 
  fESDtrackCuts->SetDCAToVertex2D(kFALSE);
  fESDtrackCuts->SetRequireSigmaToVertex(kFALSE);
  //
  // ------------------------------------------------
  // ------- track cuts to be used for v0s ----------
  // ------------------------------------------------
  //
  fESDtrackCutsCleanSamp = new AliESDtrackCuts("AliESDtrackCutsV0","AliESDtrackCutsV0");
  fESDtrackCutsCleanSamp -> SetEtaRange(-0.8,0.8);
  fESDtrackCutsCleanSamp -> SetPtRange(0.1,1e10);
  fESDtrackCutsCleanSamp -> SetMinNCrossedRowsTPC(80);
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
  cout << " ===================================================== " << endl;
  cout << " =============== Summary of Track Cuts =============== " << endl;
  cout << " ===================================================== " << endl;
  fESDtrackCuts->Dump();
  cout << " ===================================================== " << endl;

}
//________________________________________________________________________
void AliAnalysisTaskEbyeIterPID::UserCreateOutputObjects() 
{
  //
  // Create output histograms, trees of the analysis (called once)
  //
  cout << " ===== In the UserCreateOutputObjects ===== " << endl;
  // ------------  setup PIDCombined  ---------------  
  fPIDCombined=new AliPIDCombined;
  fPIDCombined->SetDefaultTPCPriors();
  fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTPC);
 
  // **********************   Input handler to get the PID object *********************
  if (!fRunFastSimulation) {
    AliAnalysisManager* man = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler* inputHandler = dynamic_cast<AliInputEventHandler*>(man->GetInputEventHandler());
    if (!inputHandler)
      AliFatal("Input handler needed");
    else {
      fPIDResponse = inputHandler->GetPIDResponse();       // PID response object
      if (!fPIDResponse) cout << " ======= PIDResponse object was not created ====== " << endl;
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
  // ************************************************************************ 
  //
  // *************** Clean Kaon, Electrons and dEdx histograms ***************** 
  // 0 --> sign
  // 1 --> Centrality
  // 2 --> eta
  // 3 --> momentum
  // 4 --> TPC dEdx
  // Clean Electrons
  //                                  0,         1,           2,           3,            4
  Int_t   binsCleanEl[nhistbins]  = { 2,  fnCentbinsData,  fnEtaBins,   fnMomBins,     dEdxnBinsClean};
  Double_t xminCleanEl[nhistbins] = {-2,   0.,             fEtaDown,    fMomDown,      fdEdxDown};
  Double_t xmaxCleanEl[nhistbins] = { 2,  80.,             fEtaUp,      fMomUp,        fdEdxCleanUp};
  fhnCleanEl = new THnSparseF("hCleanEl","Clean Electrons",nhistbins,binsCleanEl,xminCleanEl,xmaxCleanEl);
  fhnCleanEl->GetAxis(1)->Set(fnCentbinsData-1,fxCentBins);
  // Clean Kaons
  Int_t   binsCleanKa[nhistbins]  = { 2,  fnCentbinsData,  fnEtaBins,   fnMomBins,     dEdxnBinsClean};
  Double_t xminCleanKa[nhistbins] = {-2,   0.,             fEtaDown,    fMomDown,      fdEdxDown};
  Double_t xmaxCleanKa[nhistbins] = { 2,  80.,             fEtaUp,      fMomUp,        fdEdxCleanUp};
  fhnCleanKa = new THnSparseF("hCleanKa","Clean Kaons"    ,nhistbins,binsCleanKa,xminCleanKa,xmaxCleanKa);
  fhnCleanKa->GetAxis(1)->Set(fnCentbinsData-1,fxCentBins);
  // Clean Deuterons
  Int_t   binsCleanDe[nhistbins]  = { 2,  fnCentbinsData,  fnEtaBins,   fnMomBins,     dEdxnBinsClean};
  Double_t xminCleanDe[nhistbins] = {-2,   0.,             fEtaDown,    fMomDown,      fdEdxDown};
  Double_t xmaxCleanDe[nhistbins] = { 2,  80.,             fEtaUp,      fMomUp,        fdEdxCleanUp};
  fhnCleanDe = new THnSparseF("hCleanDe","Clean Deuterons",nhistbins,binsCleanDe,xminCleanDe,xmaxCleanDe);
  fhnCleanDe->GetAxis(1)->Set(fnCentbinsData-1,fxCentBins);
  // Set the branch names
  TString axisNamedEdxClean[nhistbins]   = {"sign" ,"Centrality"     ,"eta"  ,"momentum"  ,"TPC dEdx"};
  TString axisTitledEdxClean[nhistbins]  = {"sign" ,"Centrality [%]" ,"#eta" ,"#it{p} (GeV/#it{c})" ,"TPC d#it{E}/d#it{x} Signal (a.u.)"};
  for (Int_t iaxis=0; iaxis<nhistbins;iaxis++){
    fhnCleanEl->GetAxis(iaxis)->SetName(axisNamedEdxClean[iaxis]);  fhnCleanEl->GetAxis(iaxis)->SetTitle(axisTitledEdxClean[iaxis]);
    fhnCleanKa->GetAxis(iaxis)->SetName(axisNamedEdxClean[iaxis]);  fhnCleanKa->GetAxis(iaxis)->SetTitle(axisTitledEdxClean[iaxis]);
    fhnCleanDe->GetAxis(iaxis)->SetName(axisNamedEdxClean[iaxis]);  fhnCleanDe->GetAxis(iaxis)->SetTitle(axisTitledEdxClean[iaxis]);
  }
  // ************************************************************************ 
  //
  // ****************** Efficiency matrix histograms ************************ 
  const Int_t ndim=5;
  Int_t nbins0[ndim]  ={3,9,50 , 32 ,50  };
  Double_t xmin0[ndim]={0,0,0.2,-0.8,0.  };
  Double_t xmax0[ndim]={3,9,3.2, 0.8,6.25};
  if(fEffMatrix){
    fHistPosEffMatrixRec  =new THnF("fHistPosEffMatrixRec","fHistPosEffMatrixRec",ndim, nbins0,xmin0,xmax0);
    fHistNegEffMatrixRec  =new THnF("fHistNegEffMatrixRec","fHistNegEffMatrixRec",ndim, nbins0,xmin0,xmax0);
    fHistPosEffMatrixGen  =new THnF("fHistPosEffMatrixGen","fHistPosEffMatrixGen",ndim, nbins0,xmin0,xmax0);
    fHistNegEffMatrixGen  =new THnF("fHistNegEffMatrixGen","fHistNegEffMatrixGen",ndim, nbins0,xmin0,xmax0);
    fHistPosEffMatrixRec->GetAxis(1)->Set(fnCentbinsData-1,fxCentBins);
    fHistNegEffMatrixRec->GetAxis(1)->Set(fnCentbinsData-1,fxCentBins);
    fHistPosEffMatrixGen->GetAxis(1)->Set(fnCentbinsData-1,fxCentBins);
    fHistNegEffMatrixGen->GetAxis(1)->Set(fnCentbinsData-1,fxCentBins);
    TString axisNameEff[ndim]  = {"particle"      ,"Centrality"     ,"momentum"      ,"eta"  ,"phi"};
    TString axisTitleEff[ndim] = {"particle type" ,"Centrality (%)" ,"#it{p}_{T} (GeV/#it{c})" ,"#eta" ,"#phi"};
    for (Int_t iEff=0; iEff<ndim;iEff++){
      fHistPosEffMatrixRec->GetAxis(iEff)->SetName(axisNameEff[iEff]);  fHistPosEffMatrixRec->GetAxis(iEff)->SetTitle(axisTitleEff[iEff]);
      fHistNegEffMatrixRec->GetAxis(iEff)->SetName(axisNameEff[iEff]);  fHistNegEffMatrixRec->GetAxis(iEff)->SetTitle(axisTitleEff[iEff]);
      fHistPosEffMatrixGen->GetAxis(iEff)->SetName(axisNameEff[iEff]);  fHistPosEffMatrixGen->GetAxis(iEff)->SetTitle(axisTitleEff[iEff]);
      fHistNegEffMatrixGen->GetAxis(iEff)->SetName(axisNameEff[iEff]);  fHistNegEffMatrixGen->GetAxis(iEff)->SetTitle(axisTitleEff[iEff]);
    }
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
  fHistArmPod        ->GetXaxis()->SetTitle("#alpha"); fHistArmPod->GetYaxis()->SetTitle("q_{t}"); fHistArmPod->SetMarkerStyle(kFullCircle);
  fHistCentrality    ->GetXaxis()->SetTitle("Centrality (%)"); 
  fHistVertex        ->GetXaxis()->SetTitle("vertexZ (cm)"); 

  //   Add histograms to the TList
  if(fEffMatrix){
    fListHist->Add(fHistPosEffMatrixRec);
    fListHist->Add(fHistNegEffMatrixRec);
    fListHist->Add(fHistPosEffMatrixGen);
    fListHist->Add(fHistNegEffMatrixGen);
  }
  fListHist->Add(fhndEdx);
  fListHist->Add(fhnExpected);
  fListHist->Add(fhnCleanEl);
  fListHist->Add(fhnCleanKa);
  fListHist->Add(fhnCleanDe);
  fListHist->Add(fHistEmptyEvent);
  fListHist->Add(fHistCentrality);
  fListHist->Add(fHistVertex);
  fListHist->Add(fHistArmPod);

  //   OpenFile(2);  // OpenFile TIden Tree    ==============================================  
  // TIdentity informations
  fIdenTree = new TTree("fIdenTree",   "Tree for TIdentity analysis");
  fIdenTree -> Branch("sign"   ,&signNew   ,"signNew/I");
  fIdenTree -> Branch("myDeDx" ,&myDeDx    ,"myDeDx/D");
  fIdenTree -> Branch("evtNum" ,&fEventGID);
  fIdenTree -> Branch("myBin"  ,myBin      ,"myBin[3]/I");
  fIdenTree -> Branch("px"     ,&fPx );
  fIdenTree -> Branch("py"     ,&fPy );
  fIdenTree -> Branch("pz"     ,&fPz );

  //   OpenFile(3);  // OpenFile TIdenMC tree   ==============================================  
  // MC TIdentity informations
  fIdenTreeMC = new TTree("fIdenTreeMC", "Tree for MC TIdentity analysis");
  fIdenTreeMC -> Branch("sign"   ,&signNewMC   ,"signNewMC/I");
  fIdenTreeMC -> Branch("myDeDx" ,&myDeDxMC    ,"myDeDxMC/D");
  fIdenTreeMC -> Branch("evtNum" ,&fEventGIDMC);
  fIdenTreeMC -> Branch("myBin"  ,myBinMC      ,"myBinMC[3]/I");
  fIdenTreeMC -> Branch("px"     ,&fPxMC );
  fIdenTreeMC -> Branch("py"     ,&fPyMC );
  fIdenTreeMC -> Branch("pz"     ,&fPzMC );
  
  //   OpenFile(4);  // OpenFile for data tree ==============================================
  // Real Data Tree
  fTree = new TTree("fTree",      "Tree for analysis of PID");
  fTree->Branch("dEdx"       , &fTPCSignal);
  fTree->Branch("ptot"       , &fptot);
  fTree->Branch("eta"        , &fEta);
  fTree->Branch("sign"       , &fSign);
  fTree->Branch("cent"       , &fCentrality);
//   fTree->Branch("event"      , &fEventGID);

//   OpenFile(5);  // OpenFile for armPod tree     ==============================================  
// Armpod tree to be able to use graphical cut
  fArmPodTree  = new TTree("fArmPodTree",  "Tree for Clean Pion and Proton selection");
  fArmPodTree->Branch("dEdx"        , &fArmPodTPCSignal);
  fArmPodTree->Branch("ptot"        , &fArmPodptot);
  fArmPodTree->Branch("eta"         , &fArmPodEta);
  fArmPodTree->Branch("cent"        , &fArmPodCentrality);
  fArmPodTree->Branch("qt"          , &fQt);
  fArmPodTree->Branch("alfa"        , &fAlfa);
  fArmPodTree->Branch("piTOFnSigma" , &fPiNSigmasTOF);
  fArmPodTree->Branch("prTOFnSigma" , &fPrNSigmasTOF);
  
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

  cout << " ===== Out of UserCreateOutputObjects ===== " << endl;

}
//________________________________________________________________________
void AliAnalysisTaskEbyeIterPID::UserExec(Option_t *) 
{
  //
  // main event loop
  //
  cout << " ===== In the UserExec ===== " << endl;
  // Check Monte Carlo information and other access first:
  AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  if (!eventHandler) fMCtrue = kFALSE;
  //     
  // ======================================================================
  // ========================== See if MC or Real =========================
  // ======================================================================
  // impact parameters to use: 0.0, 3.72, 5.23, 7.31, 8.88, 10.20, 11.38, 12.47, 13.50, 14.5
  // corresponding Centrality:  0     5    10    20    30     40     50    60      70    80
  AliMCEvent *mcEvent = 0x0;
  if (eventHandler) mcEvent = eventHandler->MCEvent();
  // Get impact parameter for the Centrality information
  if (mcEvent){
    AliGenEventHeader* genHeader = mcEvent->GenEventHeader();
    if(!genHeader){ 
      printf("  Event generator header not available!!!\n"); 
      return; 
    }
    fMCImpactParameter = ((AliGenHijingEventHeader*) genHeader)->ImpactParameter();
  }
  
  if(fMCtrue){
    //
    // ========================== MC =========================
    //
    // impact parameters to use: 0.0, 3.72, 5.23, 7.31, 8.88, 10.20, 11.38, 12.47, 13.50, 14.5
    // corresponding Centrality:  0     5    10    20    30     40     50    60      70    80
    //
    Double_t impParArr[10] = {0.0, 3.72, 5.23, 7.31, 8.88, 10.20, 11.38, 12.47, 13.50, 14.5};
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
    //
    // Use file name in Hashing to create unique event ID
    //
    TTree *chain = (TChain*)GetInputData(0); 
    if(!chain) { Printf("ERROR: Could not receive input chain"); return; }
    TObjString fileName(chain->GetCurrentFile()->GetName());
    TString amptFileName = fileName.GetString();   
    fEventGIDMC  = TMath::Abs(Int_t(TMath::Hash(chain->GetCurrentFile()->GetName()))/2);       // uniqe id for file 
    fEventGIDMC += TMath::Abs(Int_t(fCentrality)+fEventCountInFile+(1000*fMCImpactParameter));
    fEventGID    = fEventGIDMC;
    cout << " ====================================================================================================== " << endl; 
    cout << fEventCountInFile << " ----- " << "eventID = " << fEventGIDMC << "   " << chain->GetCurrentFile()->GetName() << endl;
    cout << " Centrality = " << fCentrality << " ------ Impact Param = " << fMCImpactParameter << endl;
    cout << " ====================================================================================================== " << endl; 
    fEventCountInFile++; 
    //
  } 
  
  if (!fRunFastSimulation) {
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
    if(vertex->GetNContributors()<1) {
      vertex = fESD->GetPrimaryVertexSPD();    // SPD vertex
      TString vertexType = vertex->GetTitle();
      if (vertexType.Contains("vertexer: Z") && (vertex->GetDispersion() > 0.04 || vertex->GetZRes() > 0.25))  
	isVertexOk = kFALSE;  
      if (vertex->GetNContributors()<1)  isVertexOk = kFALSE;  
    } 
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
      if (fSystCentEstimatetor == -1) fCentrality = esdCentrality->GetCentralityPercentile("TRK");
      if (fSystCentEstimatetor ==  0) fCentrality = esdCentrality->GetCentralityPercentile("V0M");
      if (fSystCentEstimatetor ==  1) fCentrality = esdCentrality->GetCentralityPercentile("CL1");
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
    fEventCountInFile++;                 // uniqe event id within job
    ULong64_t orbitID      = (ULong64_t)fESD->GetOrbitNumber();
    ULong64_t bunchCrossID = (ULong64_t)fESD->GetBunchCrossNumber();
    ULong64_t periodID     = (ULong64_t)fESD->GetPeriodNumber();
    ULong64_t gid          = ((periodID << 36) | (orbitID << 12) | bunchCrossID); 
    fEventGID              = TMath::Abs(Int_t(TString::Hash(&gid,sizeof(Int_t))));    // uniqe event id for real data 
    cout << " =============================================================================================== " << endl; 
    cout << fEventCountInFile << " ----- " << fCentrality << " === gid =  " << gid << " hashed = " << fEventGID << endl;
    cout << " ====================================================================================================== " << endl; 
  }
  
  // ======================================================================
  //   
  // ==========================  Filling part  ============================
  if(fRunFastSimulation) {
    FastGen();
    FillDnchDeta();
  } else {
    if (fdEdxCheck){
      FillTPCdEdxCheck();
    } else {
      if (fMCtrue && !fEffMatrix) {FillTPCdEdxMC(); WeakAndMaterial();}           // dEdx fill for the MC   data
      if (fMCtrue && fEffMatrix ) FillTPCdEdxMCEffMatrix();
      FillTPCdEdxReal();                                     // dEdx fill for the real data + fill clean kaons
      if (!fMCtrue) FillCleanElectrons();              // dEdx fill for Clean Electrons
      if (!fMCtrue) FillCleanPions();                  // Fill Clean Pions and Protons vi Armanteros-Podolanski plot 
    }
  } 
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
  cout << " ===== In the FillTPCdEdxReal ===== " << endl;
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
  Int_t trackCounter = 0;
  for (Int_t i=0;i<event->GetNumberOfTracks();++i) {   // Track loop
    
    fdEdxEl=-10.; fSigmaEl=-10.; 
    fdEdxPi=-10.; fSigmaPi=-10.; 
    fdEdxKa=-10.; fSigmaKa=-10.; 
    fdEdxPr=-10.; fSigmaPr=-10.; 
    fdEdxDe=-10.; fSigmaDe=-10.;
   
    // -------------------------------------------------------------- 
    // Get the track and check if it is in the TPC
    AliESDtrack *track = fESD->GetTrack(i); 
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
    Float_t fMaxDCAToVertexXY, fMaxDCAToVertexZ;
    track->GetImpactParameters(fMaxDCAToVertexXY, fMaxDCAToVertexZ);
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
    Double_t fPhi = track->Phi();
    Double_t tpcCrossedRows     = track->GetTPCCrossedRows();
    Bool_t fTrackVeto            = fESDtrackCuts->AcceptTrack(track);
    Bool_t fKinkDaughters        = fESDtrackCuts->GetAcceptKinkDaughters();
    Bool_t fRequireTPCRefit      = fESDtrackCuts->GetRequireTPCRefit();
    Bool_t fRequireITSRefit      = fESDtrackCuts->GetRequireITSRefit(); 
    Bool_t fDCAToVertex2D        = fESDtrackCuts->GetDCAToVertex2D();
    Bool_t fRequireSigmaToVertex = fESDtrackCuts->GetRequireSigmaToVertex();
    Float_t fMaxFractionSharedTPCClusters = fESDtrackCuts->GetMaxFractionSharedTPCClusters();
    Bool_t fAcceptance = ((fEta>=fEtaDown && fEta<=fEtaUp) 
                       && (fptot>=fMomDown && fptot<=fMomUp) 
		       && (fTPCSignal>=fdEdxDown && fTPCSignal<=fdEdxUp)
		       && (fCentrality>1e-5 && fCentrality<=80));
    static Int_t  runNo = fESD->GetRunNumber();
    if (fptot>100) continue;  
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
    // Tree for the all cut variables 
    if (fFillCuts && trackCounter%50==0) {
      if(!fTreeSRedirector) return;
      (*fTreeSRedirector)<<"cuts"<<
      "run="                  << runNo                 <<         //  run Number
      "dEdx="                 << fTPCSignal            <<         //  dEdx of the track
      "theta="                << fTheta                <<         //  theta
      "phi="                  << fPhi                  <<         //  phi
      "Y="                    << fY                    <<         //  rapidity
      "px="                   << fPx                   <<         //  x momentum
      "py="                   << fPy                   <<         //  y momentum
      "pz="                   << fPz                   <<         //  z momentum
      "ptot="                 << fptot                 <<         //  total momentum 
      "pT="                   << fpT                   <<         //  tranverse momentum
      "eta="                  << fEta                  <<         //  eta
      "sign="                 << fSign                 <<         //  charge
      "ncluster="             << fNcl                  <<         //  number of points used for dEdx
      "tpcshared="            << fTPCShared            <<         //  shared clusters of the track with adjacent ones
      "lactivezone="          << lengthInActiveZone    <<         //  lengthInActiveZone in TPC 
      "crossedrows="          << tpcCrossedRows        <<         //  ncrossed rows in TPC 
      "DCAxy="                << fMaxDCAToVertexXY     <<        
      "DCAz="                 << fMaxDCAToVertexZ      <<        
      "trackVeto="            << fTrackVeto            <<         //  bool: track validity flag
//       "requireTPCRefit="      << fRequireTPCRefit      <<        
//       "requireITSRefit="      << fRequireITSRefit      <<         
//       "DCAtoVertex2D="        << fDCAToVertex2D        <<         
//       "requireSigmaToVertex=" << fRequireSigmaToVertex <<   
      "\n";   
    }  
    // -------------------------------------------------------------- 
    //
    // -------------------------------------------------------------- 
    // Apply track cuts
    if (fptot<fMomDown || fptot>fMomUp)            continue;   // ptot range to be used 
    if (!fESDtrackCuts->AcceptTrack(track))        continue;        // standard track cuts 
    if (fTightCuts) if (track->GetTPCsignalN()<70) continue;        // from jens to increase pid resolution
    if (fTightCuts) if (track->GetLengthInActiveZone(1,3,230, track->GetBz(),0,0)<120) continue;  // from Marian 
    trackCounter++;
    // -------------------------------------------------------------- 
    //
    // -------------------------------------------------------------- 
    // Fill the trees
    myDeDx   = fTPCSignal;
    myBin[0] = fhEta ->FindBin(fEta)-1;
    myBin[1] = fhCent->FindBin(fCentrality)-1;
    myBin[2] = fhPtot->FindBin(fptot)-1;
    signNew  = fSign;
    if (!fCleanSamplesOnly && !fMCtrue && fAcceptance)   fIdenTree -> Fill();
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
        if (!fMCtrue) fhnCleanKa->Fill(weightCleanKa); 
      } 
    } 
    // -------------------------------------------------------------- 
    // --------------------------------------------------------------       
    // Fill clean Deuterons
    if ((TMath::Abs(nSigmasDeTOF)<=3) && TMath::Abs(fNSigmasDeTPC)<3 && (!fMCtrue)) {
      Double_t weightCleanDe[5] = {Double_t(fSign),fCentrality,fEta,fptot, fTPCSignal};
      if (!fMCtrue) fhnCleanDe->Fill(weightCleanDe); 
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
  cout << " ===== In the FillTPCdEdxCheck ===== " << endl;
  if (!fPIDResponse) fPIDResponse = ((AliESDInputHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->GetESDpid();
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
    if (fTPCSignal>400)         continue;        
    if (fptot>2.)               continue;        
    if (fEta<-0.8 && fEta>0.8)  continue;        
    if (track->GetTPCNcls()<80) continue;       

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
  cout << " ===== In the FillTPCdEdxMC ===== " << endl;
  AliMCEventHandler *eventHandler = dynamic_cast<AliMCEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  TTree *chain = (TChain*)GetInputData(0);
  if(!chain) { Printf("ERROR: Could not receive input chain"); return; }
  TObjString fileName(chain->GetCurrentFile()->GetName());
  Int_t runNumber = fESD->GetRunNumber();
  
  // check if the event is empty 
  if (CountEmptyEvents(7)<1) return;
  
  // =========== get MC event ===========
  AliMCEvent *mcEvent = 0x0;
  AliStack   *stack   = 0x0;     // All particles in an MC event
  if (eventHandler) mcEvent = eventHandler->MCEvent();
  if (!mcEvent) { Printf("ERROR: Could not retrieve MC event"); if (fMCtrue) return; }
  if (fMCtrue)  { stack = mcEvent->Stack(); if (!stack) return; }
  
  // ======================================================================
  // 
  // Fill TIdentity tree for MC + Fill MC closure information
  for(Int_t i = 0; i < fESD->GetNumberOfTracks(); i++) {    
    
    // initialize the dummy particle id
    fElMC =-100.; fPiMC =-100.; fKaMC =-100.; fPrMC =-100.; fDeMC =-100.; fMuMC =-100.; 
    AliESDtrack *trackReal = fESD->GetTrack(i); 
    Int_t lab = TMath::Abs(trackReal->GetLabel());           // avoid from negatif labels, they include some garbage
    Bool_t fTrackVeto = fESDtrackCuts->AcceptTrack(trackReal);

    Float_t dcaXYprim, dcaZprim;
    Float_t dcaXYweak, dcaZweak;
    Float_t dcaXYmaterial, dcaZmaterial;
    if(stack->IsSecondaryFromMaterial(lab))  trackReal->GetImpactParameters(dcaXYmaterial, dcaZmaterial);
    if(stack->IsSecondaryFromWeakDecay(lab)) trackReal->GetImpactParameters(dcaXYweak, dcaZweak);
    if(stack->IsPhysicalPrimary(lab))        trackReal->GetImpactParameters(dcaXYprim, dcaZprim);
    
    // MC track cuts
    if (!stack->IsPhysicalPrimary(lab)) continue;  // MC primary track check
    
    // Tack cuts from detector
    if (!trackReal -> GetInnerParam())                 continue;
    if (!fESDtrackCuts -> AcceptTrack(trackReal))      continue;  // real track cuts 
    if (fTightCuts) if (trackReal->GetTPCsignalN()<70) continue;                                          
    if (fTightCuts) if (trackReal->GetLengthInActiveZone(1,3,230, trackReal->GetBz(),0,0)<120) continue;  
    
    // match the track with mc track
    TParticle *trackMC  = stack->Particle(lab);   
    Int_t pdg           = trackMC->GetPdgCode();
    
    Int_t iPart = -10;
    if (TMath::Abs(pdg) == 11)          iPart = 0; // select el-
    if (TMath::Abs(pdg) == 211)         iPart = 1; // select pi+
    if (TMath::Abs(pdg) == 321)         iPart = 2; // select ka+
    if (TMath::Abs(pdg) == 2212)        iPart = 3; // select pr+
    //     if (TMath::Abs(pdg) == 1000010020)  iPart = 4; // select de
    //     if (TMath::Abs(pdg) == 13)          iPart = 5; // select mu-
    if (iPart == -10) continue;
    
    if (iPart == 0 ) fElMC = trackReal->GetTPCsignal();
    if (iPart == 1 ) fPiMC = trackReal->GetTPCsignal();
    if (iPart == 2 ) fKaMC = trackReal->GetTPCsignal();
    if (iPart == 3 ) fPrMC = trackReal->GetTPCsignal();
    //     if (iPart == 4 ) fDeMC = trackReal->GetTPCsignal();
    //     if (iPart == 5 ) fMuMC = trackReal->GetTPCsignal();
    
    fEtaMC        = trackReal->Eta();
    fTPCSignalMC  = trackReal->GetTPCsignal();
    fptotMC       = trackReal->GetInnerParam()->GetP();
    fpTMC         = trackReal->Pt();
    fSignMC       = trackReal->GetSign();
    fPxMC         = trackReal->Px();
    fPyMC         = trackReal->Py();
    fPzMC         = trackReal->Pz();

    // Fill the ttree for TIdentity study and fit systematics
    if (fCentrality>=0. && fCentrality<80.){
      myDeDxMC        = fTPCSignalMC;
      myBinMC[0]      = fhEta ->FindBin(fEtaMC)-1;
      myBinMC[1]      = fhCent->FindBin(fCentrality)-1;
      myBinMC[2]      = fhPtot->FindBin(fptotMC)-1;
      signNewMC       = fSignMC;
      fIdenTreeMC->Fill();
      // Fill MC closure tree
      if(!fTreeSRedirector) return;
      (*fTreeSRedirector)<<"fTreeMC"<<
      "dEdx="      << fTPCSignalMC <<    // dEdx of mc track
      "ptot="      << fptotMC <<         // mc momentum
      "pT="        << fpTMC <<         // mc momentum
      "eta="       << fEtaMC <<          // mc eta
      "cent="      << fCentrality <<     // Centrality
      "sign="      << fSignMC <<         // sign
      "event="     << fEventGIDMC <<
      "trackVeto=" << fTrackVeto            <<         //  bool: track validity flag
      "el="        << fElMC <<           // electron dEdx
      "pi="        << fPiMC <<           // pion dEdx
      "ka="        << fKaMC <<           // kaon dEdx
      "pr="        << fPrMC <<           // proton dEdx
      "dcaXYprim="      << dcaXYprim <<           // electron dEdx
      "dcaXYweak="      << dcaXYweak <<           // electron dEdx
      "dcaXYmaterial="  << dcaXYmaterial <<           // electron dEdx
      "dcaZprim="       << dcaZprim <<           // electron dEdx
      "dcaZweak="       << dcaZweak <<           // electron dEdx
      "dcaZmaterial="   << dcaZmaterial <<           // electron dEdx
 
      "\n";   
    
    }  
  } // ======= end of track loop =======
  //   
  // ======================================================================
  // ======================================================================
  //   
  // ========= Efficiency Check Eta momentum and Centrality scan ==========
  const Int_t nMoments = 9;
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
	for(Int_t i=0;i<nMoments; i++){  recMoments[i]=0.; recMomentsPos[i]=0.;  recMomentsNeg[i]=0.;}
	
	// loop over tracks
        Int_t nTracks=0, trCount=0;
        for(Int_t i = 0; i < fESD->GetNumberOfTracks(); i++) {    
    
          // initialize the dummy particle id
          fElMC =-100.; fPiMC =-100.; fKaMC =-100.; fPrMC =-100.; fDeMC =-100.; fMuMC =-100.; 
          AliESDtrack *trackReal = fESD->GetTrack(i); 
          Int_t lab = TMath::Abs(trackReal->GetLabel());           // avoid from negatif labels, they include some garbage
          fEtaMC    = trackReal->Eta();
        
          // MC track cuts
          if ((fEtaMC<fetaDownArr[ieta]) || (fEtaMC>fetaUpArr[ieta])) continue;  // eta Cut
          if (!stack->IsPhysicalPrimary(lab)) continue;                         // MC primary track check
        
          // Tack cuts from detector
          if (!trackReal -> GetInnerParam()) continue;
          if (!fESDtrackCuts -> AcceptTrack(trackReal)) continue;  // real track cuts 
          
          if (fTightCuts) if (trackReal->GetTPCsignalN()<70) continue;                                          
          if (fTightCuts) if (trackReal->GetLengthInActiveZone(1,3,230, trackReal->GetBz(),0,0)<120) continue;  
 
          // match the track with mc track
          TParticle *trackMC  = stack->Particle(lab);   
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
    
          // count first moments for given Centrality and momentum window
	  fptotMC = trackReal->GetInnerParam()->GetP();
          if ( (fCentrality>=fcentDownArr[icent]) 
               && (fCentrality<fcentUpArr[icent]) 
               && (fptotMC>=fpDownArr[imom]) 
               && (fptotMC<=fpUpArr[imom]) ) 
          { 
	    nTracks++;                                // count the particles after track cuts
	    if ( fPiMC>-1 || fKaMC>-1 || fPrMC>-1) trCount++;
            if ( fPiMC>-1 ) recMoments[kPi]++;
            if ( fKaMC>-1 ) recMoments[kKa]++;
            if ( fPrMC>-1 ) recMoments[kPr]++;
	    if ( fPiMC>-1 && pdg<0) recMomentsNeg[kPi]++;
            if ( fKaMC>-1 && pdg<0) recMomentsNeg[kKa]++;
            if ( fPrMC>-1 && pdg<0) recMomentsNeg[kPr]++;
	    if ( fPiMC>-1 && pdg>0) recMomentsPos[kPi]++;
            if ( fKaMC>-1 && pdg>0) recMomentsPos[kKa]++;
            if ( fPrMC>-1 && pdg>0) recMomentsPos[kPr]++;
          } 
          
        } // ======= end of track loop =======
    
        // calculate second moments
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

      
        // fill tree which contains moments
        Int_t dataType = 0, sampleNo = 0;
        Double_t centBin = (fcentDownArr[icent]+fcentUpArr[icent])/2.;
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
	for(Int_t i=0;i<nMoments; i++){  genMoments[i]=0.; genMomentsPos[i]=0.;  genMomentsNeg[i]=0.;}
	
	// go into event loop
        Int_t nTracksgen=0, trCountgen=0; 
	AliMCParticle *trackMCgen;
        for (Int_t iTracks = 0; iTracks < fMCEvent->GetNumberOfTracks(); iTracks++) {    // track loop
  
        // initialize the dummy particle id
          fElMCgen =-100.; fPiMCgen =-100.; fKaMCgen =-100.; fPrMCgen =-100.; fDeMCgen =-100.; fMuMCgen =-100.;
          trackMCgen = (AliMCParticle *)fMCEvent->GetTrack(iTracks);
    
         // apply primary vertex and eta cut
          if ((trackMCgen->Eta()<fetaDownArr[ieta]) || (trackMCgen->Eta()>fetaUpArr[ieta])) continue;
          if (!stack->IsPhysicalPrimary(iTracks)) continue;
          Int_t pdg  = trackMCgen->Particle()->GetPdgCode();
        
          Int_t iPart = -10;
          if (TMath::Abs(pdg) == 11)          {iPart = 0; fElMCgen = iPart;} // select el-
          if (TMath::Abs(pdg) == 211)         {iPart = 1; fPiMCgen = iPart;} // select pi+
          if (TMath::Abs(pdg) == 321)         {iPart = 2; fKaMCgen = iPart;} // select ka+
          if (TMath::Abs(pdg) == 2212)        {iPart = 3; fPrMCgen = iPart;} // select pr+
          if (TMath::Abs(pdg) == 1000010020)  {iPart = 4; fDeMCgen = iPart;} // select de
          if (TMath::Abs(pdg) == 13)          {iPart = 5; fMuMCgen = iPart;} // select mu-
          if (iPart == -10) continue;
          fptotMCgen = trackMCgen->P();      
       
          // count first moments
          if ((fCentrality>=fcentDownArr[icent])
               &&(fCentrality<fcentUpArr[icent])
               &&(fptotMCgen>=fpDownArr[imom])
               &&(fptotMCgen<=fpUpArr[imom])) 
          { 
	    nTracksgen++;
	    if ( fPiMCgen>-1 || fKaMCgen>-1 || fPrMCgen>-1) trCountgen++;
            if ( fPiMCgen>-1 ) genMoments[kPi]++;
            if ( fKaMCgen>-1 ) genMoments[kKa]++;
            if ( fPrMCgen>-1 ) genMoments[kPr]++;
	    if ( fPiMCgen>-1 && pdg<0) genMomentsNeg[kPi]++;
            if ( fKaMCgen>-1 && pdg<0) genMomentsNeg[kKa]++;
            if ( fPrMCgen>-1 && pdg<0) genMomentsNeg[kPr]++;
	    if ( fPiMCgen>-1 && pdg>0) genMomentsPos[kPi]++;
            if ( fKaMCgen>-1 && pdg>0) genMomentsPos[kKa]++;
            if ( fPrMCgen>-1 && pdg>0) genMomentsPos[kPr]++;
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
  cout << " ===== In the FastGen ===== " << endl;
  AliMCEventHandler *eventHandler = dynamic_cast<AliMCEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
   
  // =========== get MC event ===========
  AliMCEvent *mcEvent = 0x0;
  AliStack   *stack   = 0x0;     // All particles in an MC event
  if (eventHandler) mcEvent = eventHandler->MCEvent();
  if (!mcEvent) { Printf("ERROR: Could not retrieve MC event"); if (fMCtrue) return; }
  if (fMCtrue)  { stack = mcEvent->Stack(); if (!stack) return; }
  Int_t dataType = 1, sampleNo = 0;  // dataType-> 1 for MCgen
  // ======================================================================
  // ======================================================================
  //   
  // ========= Efficiency Check Eta momentum and Centrality scan ==========
  const Int_t nMoments = 9;
   for (Int_t ieta=0; ieta<fnEtaWinBinsMC; ieta++){
    for (Int_t imom=0; imom<fnMomBinsMC; imom++){
      for (Int_t icent=0; icent<fnCentbinsData-1; icent++){
        
	// vectors to hold moments
	TVectorF genMoments(nMoments);
	TVectorF genMomentsPos(nMoments);
	TVectorF genMomentsNeg(nMoments);
	TVectorF genMomentsCross(nMoments);
	for(Int_t i=0;i<nMoments; i++){  genMoments[i]=0.; genMomentsPos[i]=0.;  genMomentsNeg[i]=0.;}
	  
	Float_t nTracksgen=0, trCountgen=0;
	AliMCParticle *trackMCgen;
        for (Int_t iTracks = 0; iTracks < fMCEvent->GetNumberOfTracks(); iTracks++) {    // track loop
  
	  // initialize the dummy particle id
          fElMCgen =-100.; fPiMCgen =-100.; fKaMCgen =-100.; fPrMCgen =-100.; fDeMCgen =-100.; fMuMCgen =-100.;
          trackMCgen = (AliMCParticle *)fMCEvent->GetTrack(iTracks);
    
	  // apply primary vertex and eta cut
          if ((trackMCgen->Eta()<fetaDownArr[ieta]) || (trackMCgen->Eta()>fetaUpArr[ieta])) continue;
	  // if ( !(stack->IsPhysicalPrimary(iTracks) || stack->IsSecondaryFromWeakDecay(iTracks)) ) continue;
	  if (!stack->IsPhysicalPrimary(iTracks)) continue;
          Int_t pdg  = trackMCgen->Particle()->GetPdgCode();
        
          Int_t iPart = -10;
          if (TMath::Abs(pdg) == 11)          {iPart = 0; fElMCgen = iPart;} // select el-
          if (TMath::Abs(pdg) == 211)         {iPart = 1; fPiMCgen = iPart;} // select pi+
          if (TMath::Abs(pdg) == 321)         {iPart = 2; fKaMCgen = iPart;} // select ka+
          if (TMath::Abs(pdg) == 2212)        {iPart = 3; fPrMCgen = iPart;} // select pr+
          if (TMath::Abs(pdg) == 1000010020)  {iPart = 4; fDeMCgen = iPart;} // select de
          if (TMath::Abs(pdg) == 13)          {iPart = 5; fMuMCgen = iPart;} // select mu-
          if (iPart == -10) continue;
	  fptotMCgen = trackMCgen->P();  
       
          // count first moments
          if ((fCentrality>=fcentDownArr[icent])
               &&(fCentrality<fcentUpArr[icent])
               &&(fptotMCgen>=fpDownArr[imom])
               &&(fptotMCgen<=fpUpArr[imom])) 
          { 
	    nTracksgen++;
	    if ( fPiMCgen>-1 || fKaMCgen>-1 || fPrMCgen>-1) trCountgen++;
            if ( fPiMCgen>-1 ) genMoments[kPi]++;
            if ( fKaMCgen>-1 ) genMoments[kKa]++;
            if ( fPrMCgen>-1 ) genMoments[kPr]++;
	    if ( fPiMCgen>-1 && pdg<0) genMomentsNeg[kPi]++;
            if ( fKaMCgen>-1 && pdg<0) genMomentsNeg[kKa]++;
            if ( fPrMCgen>-1 && pdg<0) genMomentsNeg[kPr]++;
	    if ( fPiMCgen>-1 && pdg>0) genMomentsPos[kPi]++;
            if ( fKaMCgen>-1 && pdg>0) genMomentsPos[kKa]++;
            if ( fPrMCgen>-1 && pdg>0) genMomentsPos[kPr]++; 
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
  
        // fill tree which contains moments
	Double_t centBin = (fcentDownArr[icent]+fcentUpArr[icent])/2.;
        if(!fTreeSRedirector) return;
        // if there is at least one track in an event fill the tree
        if ( trCountgen>0 ){   
          (*fTreeSRedirector)<<"mcGen"<<
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
          "moment.="      << &genMoments <<             // second moments for particle+antiparticle
          "momentPos.="   << &genMomentsPos <<          // second moment of positive particles
          "momentNeg.="   << &genMomentsNeg <<          // second moment of negative particles
          "momentCross.=" << &genMomentsCross <<        // second moment of unlikesign particles
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
  cout << " ===== In the WeakAndMaterial ===== " << endl;
  AliMCEventHandler *eventHandler = dynamic_cast<AliMCEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
   
  // =========== get MC event ===========
  AliMCEvent *mcEvent = 0x0;
  AliStack   *stack   = 0x0;     // All particles in an MC event
  if (eventHandler) mcEvent = eventHandler->MCEvent();
  if (!mcEvent) { Printf("ERROR: Could not retrieve MC event"); if (fMCtrue) return; }
  if (fMCtrue)  { stack = mcEvent->Stack(); if (!stack) return; }
  // ======================================================================
  // ======================================================================
  //   
  // ========= Efficiency Check Eta momentum and Centrality scan ==========
  
  AliMCParticle *trackMCgen;
  for (Int_t iTracks = 0; iTracks < fMCEvent->GetNumberOfTracks(); iTracks++) {    // track loop
    
    // initialize the dummy particle id
    Int_t fElWeak =-100., fPiWeak =-100., fKaWeak =-100., fPrWeak =-100., fDeWeak =-100., fMuWeak =-100.;
    trackMCgen = (AliMCParticle *)fMCEvent->GetTrack(iTracks);
      
    // apply primary vertex and eta cut
    if ( !stack->IsPhysicalPrimary(iTracks) ) continue;
    
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
  cout << " ===== In the FillDnchDeta ===== " << endl;
  AliMCEventHandler *eventHandler = dynamic_cast<AliMCEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  
  // =========== get MC event ===========
  AliMCEvent *mcEvent = 0x0;
  AliStack   *stack   = 0x0;     // All particles in an MC event
  if (eventHandler) mcEvent = eventHandler->MCEvent();
  if (!mcEvent) { Printf("ERROR: Could not retrieve MC event"); if (fMCtrue) return; }
  if (fMCtrue)  { stack = mcEvent->Stack(); if (!stack) return; }
  // ======================================================================
  // ======================================================================
  //   
  // ========= Efficiency Check Eta momentum and Centrality scan ==========
  Double_t etaDownArray[8] ={-0.1, -0.2, -0.3, -0.4, -0.5, -0.6, -0.7, -0.8};
  Double_t etaUpArray[8]   ={ 0.1,  0.2,  0.3,  0.4,  0.5,  0.6,  0.7,  0.8};
  Double_t centDownArray[9]={0., 5.,  10., 20., 30., 40., 50., 60., 70.};
  Double_t centUpArray[9]  ={5., 10., 20., 30., 40., 50., 60., 70., 80.};
  for (Int_t ieta=0; ieta<8; ieta++){
    for (Int_t icent=0; icent<9; icent++){
      
      AliMCParticle *trackMCgen;
      Int_t trCount=0,    elCount=0,    piCount=0,    kaCount=0,    prCount=0;
      for (Int_t iTracks = 0; iTracks < fMCEvent->GetNumberOfTracks(); iTracks++) {    // track loop
	
	// initialize the dummy particle id
	trackMCgen = (AliMCParticle *)fMCEvent->GetTrack(iTracks);
	
	// apply primary vertex and eta cut
	if (!stack->IsPhysicalPrimary(iTracks)) continue;
	if ((trackMCgen->Eta()<-0.8) || (trackMCgen->Eta()>0.8)) continue;
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
  cout << " ===== In the FillTPCdEdxMCEffMatrix ===== " << endl;
  AliMCEventHandler *eventHandler = dynamic_cast<AliMCEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  TTree *chain = (TChain*)GetInputData(0);
  if(!chain) { Printf("ERROR: Could not receive input chain"); return; }
  TObjString fileName(chain->GetCurrentFile()->GetName());
  Int_t runNumber = fESD->GetRunNumber();
   
  // =========== get MC event ===========
  AliMCEvent *mcEvent = 0x0;
  AliStack   *stack   = 0x0;     // All particles in an MC event
  if (eventHandler) mcEvent = eventHandler->MCEvent();
  if (!mcEvent) { Printf("ERROR: Could not retrieve MC event"); if (fMCtrue) return; }
  if (fMCtrue)  { stack = mcEvent->Stack(); if (!stack) return; }
  
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
    if ((fEtaMC<-0.8) || (fEtaMC>0.8))                 continue; 
    if (!stack->IsPhysicalPrimary(lab))                continue;  // MC primary track check
    // Track cuts from detector
    if (!trackReal -> GetInnerParam())                 continue;
    if (!fESDtrackCuts -> AcceptTrack(trackReal))      continue;  // real track cuts 
    if (fTightCuts) if (trackReal->GetTPCsignalN()<70) continue;                                          
    if (fTightCuts) if (trackReal->GetLengthInActiveZone(1,3,230, trackReal->GetBz(),0,0)<120) continue; 
    // get track info
    Float_t fpTRec   = trackReal->Pt();
    Float_t fYRec    = trackReal->Y();
    Float_t fEtaRec  = trackReal->Eta();
    Float_t fPhiRec  = trackReal->Phi(); 
    Float_t fCentRec = fhCent->FindBin(fCentrality)-1;
    Float_t fPartID  = -100;
    
    // Efficiency matices for individual particles
    TParticle *trackMC  = stack->Particle(lab);   
    Int_t pdg           = trackMC->GetPdgCode();  
    
    if (TMath::Abs(pdg) == 211)   fPartID=0; // select pi
    if (TMath::Abs(pdg) == 321)   fPartID=1; // select ka
    if (TMath::Abs(pdg) == 2212)  fPartID=2; // select pr
   
    Double_t xxxRec[5]={fPartID,fCentRec,fpTRec,fEtaRec,fPhiRec};
    if (pdg>0 && fPartID>-1) fHistPosEffMatrixRec->Fill(xxxRec);
    if (pdg<0 && fPartID>-1) fHistNegEffMatrixRec->Fill(xxxRec);

  } // ======= end of track loop =======
  
  // -----------------------------------------------------------------------------------------
  // ----------------------------   MC generated pure MC particles  --------------------------
  // -----------------------------------------------------------------------------------------
  AliMCParticle *trackMCgen;
  for (Int_t iTracks = 0; iTracks < fMCEvent->GetNumberOfTracks(); iTracks++) 
  { // track loop
    trackMCgen = (AliMCParticle *)fMCEvent->GetTrack(iTracks);
    // apply primary vertex and eta cut
    if ((trackMCgen->Eta()<-0.8) || (trackMCgen->Eta()>0.8)) continue;
    if (!stack->IsPhysicalPrimary(iTracks)) continue;
    // get track info
    Float_t fpTGen   = trackMCgen->Pt();
    Float_t fYGen    = trackMCgen->Y();
    Float_t fEtaGen  = trackMCgen->Eta();
    Float_t fPhiGen  = trackMCgen->Phi(); 
    Float_t fCentGen = fhCent->FindBin(fCentrality)-1;
    Float_t fPartID  = -100;
        
    // Efficiency matices for individual particles
    Int_t pdg  = trackMCgen->Particle()->GetPdgCode();    
    if (TMath::Abs(pdg) == 211)   fPartID=0; // select pi
    if (TMath::Abs(pdg) == 321)   fPartID=1; // select ka
    if (TMath::Abs(pdg) == 2212)  fPartID=2; // select pr
   
    Double_t xxxGen[5]={fPartID,fCentGen,fpTGen,fEtaGen,fPhiGen};
    if (pdg>0 && fPartID>-1) fHistPosEffMatrixGen->Fill(xxxGen);
    if (pdg<0 && fPartID>-1) fHistNegEffMatrixGen->Fill(xxxGen);
  } // ======= end of track loop ======= 
      
}
//________________________________________________________________________
void AliAnalysisTaskEbyeIterPID::FillCleanElectrons()
{

// Fill Clean Electrons from conversion
  cout << " ===== In the FillCleanElectrons ===== " << endl;
  if (fPIDResponse) {
    fPIDResponse->GetTPCResponse().SetBetheBlochParameters(1.28778e+00/50., 3.13539e+01, TMath::Exp(-3.16327e+01), 1.87901e+00, 6.41583e+00);
  }

  TObjArray* listCrossV0 = fESDtrackCutsV0->GetAcceptedV0s(fESD);  
  Int_t nGoodV0s         = listCrossV0->GetEntries();
 
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

    Double_t thetaP = acos((vecP * vecM)/(vecP.Mag() * vecM.Mag()));
    Double_t thetaN = acos((vecN * vecM)/(vecN.Mag() * vecM.Mag()));
    Double_t alfaEl = ((vecP.Mag())*cos(thetaP)-(vecN.Mag())*cos(thetaN))/((vecP.Mag())*cos(thetaP)+(vecN.Mag())*cos(thetaN));
    Double_t qtEl   = vecP.Mag()*sin(thetaP);

    fV0s->ChangeMassHypothesis(22);

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
      if (!fMCtrue) fhnCleanEl->Fill(weightCleanEl); 
    } 
  } // end of V0 loop
}
//________________________________________________________________________
void AliAnalysisTaskEbyeIterPID::FillCleanPions()
{

  // Fill Clean Pions from K0s
  cout << " ===== In the FillCleanPions ===== " << endl;
  if (fPIDResponse) {
    fPIDResponse->GetTPCResponse().SetBetheBlochParameters(1.28778e+00/50., 3.13539e+01, TMath::Exp(-3.16327e+01), 1.87901e+00, 6.41583e+00);
  }

  TObjArray* listCrossV0   = fESDtrackCutsV0->GetAcceptedV0s(fESD);  
  Int_t nGoodV0s      = listCrossV0->GetEntries();
 
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

    Double_t thetaP  = acos((vecP * vecM)/(vecP.Mag() * vecM.Mag()));
    Double_t thetaN  = acos((vecN * vecM)/(vecN.Mag() * vecM.Mag()));
    fAlfa = ((vecP.Mag())*cos(thetaP)-(vecN.Mag())*cos(thetaN))/((vecP.Mag())*cos(thetaP)+(vecN.Mag())*cos(thetaN));
    fQt   = vecP.Mag()*sin(thetaP);
    fHistArmPod->Fill(fAlfa,fQt);
    
    // main armentoros podolanki cuts
    if (fQt >0.22) continue;
    if (TMath::Abs(fAlfa)>0.9) continue;
    if (fQt<0.07) continue;
    if (fQt>0.11 && fQt<0.16) continue;
    if (TMath::Abs(fAlfa)<0.5 && fQt<0.15) continue;
    
    fV0MIs->ChangeMassHypothesis(310);
    
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
      
      if (fArmPodTPCSignal<30 || fArmPodTPCSignal>200) continue;
      if (fArmPodptot<fMomDown || fArmPodptot>fMomUp) continue;
       // Fill the THnSparseF for the inclusive dEdx spectrum
      if ( (TMath::Abs(fPiNSigmasTOF)<3) || (TMath::Abs(fPrNSigmasTOF)<3)) {
        fArmPodCentrality = fCentrality;
        if(fFillArmPodTree) fArmPodTree->Fill();  
      }   
    }
    
  } // end of V0 loop

}
//________________________________________________________________________
void AliAnalysisTaskEbyeIterPID::BinLogAxis(TH1 *h) 
{
  //
  // Method for the correct logarithmic binning of histograms
  //
  cout << " ===== In the BinLogAxis ===== " << endl;
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
Int_t AliAnalysisTaskEbyeIterPID::CountEmptyEvents(Int_t counterBin){

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
    cout << "empty event in " << chain->GetCurrentFile()->GetName() << endl;
  }
  cout << " ====== EVENT IS COOL GO AHEAD ======= " << endl;
  return emptyCount;

}
//________________________________________________________________________
void AliAnalysisTaskEbyeIterPID::Terminate(Option_t *) 
{
  cout << " ===== In the Terminate ===== " << endl;
  // Draw result to the screen
  // Called once at the end of the query
  
}
