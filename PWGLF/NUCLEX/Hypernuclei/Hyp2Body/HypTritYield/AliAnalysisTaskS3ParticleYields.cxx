
//--- Task for the determination of p and Lambda Yields in pp ---//
//---     Author: Janik Ditzel; janik.ditzel@cern.ch     ---//


#include <iostream>
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliMultSelection.h"
#include "AliInputEventHandler.h"
#include "AliESDInputHandler.h"
#include "AliESDtrackCuts.h"
#include "AliESDtrack.h"
#include "AliESDpid.h"
#include "AliESDVertex.h"
#include "AliESDEvent.h"
#include "AliESDv0.h"
#include "AliESDVertex.h"
#include "AliVertexerTracks.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCVertex.h"
#include "AliMCParticle.h"
#include "TPDGCode.h"
#include "AliEventCuts.h"
#include "TList.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskS3ParticleYields.h"

class AliAnalysisTaskS3ParticleYields;
using namespace std;

ClassImp(AliAnalysisTaskS3ParticleYields)

// Default Constructor
AliAnalysisTaskS3ParticleYields::AliAnalysisTaskS3ParticleYields()
:AliAnalysisTaskSE("AliAnalysisTaskS3ParticleYields"),    
  fInputHandler(0),
  fPID(0),
  fESDevent(0),
  fStack(),
  fV0(),  
  mcEvent(0),
  trackCutsV0(),
  fHistdEdx(0),
  fHistdEdxV0(0),
  fHistNumEvents(0),
  fHistTrigger(0), 
  fTree(0),
  hTree(0),
  fTreeGen(0),
  fHistogramList(NULL),
  fPrimaryVertex(),
  fPIDCheckOnly(kFALSE),
  fMCtrue(0),
  fEventCuts(),
  fBetheParamsHe(),
  fBetheParamsT(),
  MB(0),
  HMV0(0),
  HMSPD(0),
  HNU(0),
  HQU(0),
  fNumberV0s(-99),
  fonTheFly(-99),
  frunnumber(-99),
  fTrigMB(-99),			
  fTrigHMV0(-99),
  fTrigHMSPD(-99),
  fTrigHNU(0),
  fTrigHQU(0),
  fNTracks(0),
  fMultV0M(-99),
  fMultOfV0M(-99),
  fMultSPDTracklet(-99),
  fMultSPDCluster(-99),
  fMultRef05(-99),
  fMultRef08(-99),
  fSPDCluster(-99),
  fSPDTracklets(-99),
  fSPDFiredChips0(-99),
  fSPDFiredChips1(-99),
  fV0Multiplicity(-99),
  fYear(0),
  fTRDvalidP(0),
  fTRDtrigHNUP(0),
  fTRDtrigHQUP(0),
  fTRDPidP(0),
  fTRDnTrackletsP(0),
  fTRDPtP(0),
  fTRDLayerMaskP(0),
  fTRDSagittaP(-1),
  fTRDvalidPi(0),
  fTRDtrigHNUPi(0),
  fTRDtrigHQUPi(0),
  fTRDPidPi(0),
  fTRDnTrackletsPi(0),
  fTRDPtPi(0),
  fTRDLayerMaskPi(0),
  fTRDSagittaPi(-1),
  fMCtrueP(-99),
  fisPrimaryP(-99),
  fisWeakP(-99),
  fisMaterialP(-99),
  fMCtrueL(-99),
  fisWeakL(-99),
  fisMaterialL(-99),
  fisPrimaryL(-99),
  fCharge(-99),
  fLambdaM(-99),
  fLambdaP(-99),
  fLambdaPx(-99),
  fLambdaPy(-99),
  fLambdaPz(-99),
  fLambdaPt(-99),
  fLambdaE(-99),
  fLambdaCt(-99),
  fLambdaY(-99),
  fLambdaDca(-99),
  fLambdaCos(-99),
  fpiP(-99),
  fpiPt(-99),
  fpiPx(-99),
  fpiPy(-99),
  fpiPz(-99),
  fpiE(-99),
  fpiy(-99),
  fpiDedx(-99),
  fpiDedxSigma(-99),
  fpichi2(-99),
  fpiNcls(-99),
  fpiNclsITS(-99),  
  fpiDca(-99),
  fpiDcaz(-99),
  fpiDcaSec(-99),
  fGeoLengthPi(-99),  
  fTOFSignalPi(-99),
  fTOFnSigmaPi(-99),
  fEtaPi(-99),
  fPhiPi(-99),
  fKinkPi(-1),
  fTPCrefitPi(-1),
  fITSrefitPi(-1),
  fITSchi2Pi(-99),
  fpiSigmaYX(-99),
  fpiSigmaXYZ(-99),
  fpiSigmaZ(-99),
  fpP(-99),
  fpy(-99),
  fpPt(-99),
  fpPx(-99),
  fpPy(-99),
  fpPz(-99),
  fpE(-99),
  fpDedx(-99),
  fpDedxSigma(-99),
  fpchi2(-99),
  fpNcls(-99),
  fpNclsITS(-99),
  fpDcaSec(-99),
  fpDca(-99),
  fpDcaz(-99),
  fGeoLengthP(-99),
  fTOFSignalP(-99),
  fTOFnSigmaP(-99),
  fEtaP(-99),
  fPhiP(-99),
  fKinkP(-1),
  fTPCrefitP(-1),
  fITSrefitP(-1),
  fITSchi2P(-99),
  fpSigmaYX(-99),
  fpSigmaXYZ(-99),
  fpSigmaZ(-99),
  farmalpha(-99),
  farmpt(-99), 
  fthetaP(-99),
  fthetaN(-99)
{
    
}

// Constructor
AliAnalysisTaskS3ParticleYields::AliAnalysisTaskS3ParticleYields(const char *name)
  :AliAnalysisTaskSE(name),
   fInputHandler(0),
  fPID(0),
  fESDevent(0),
  fStack(),
  fV0(),  
  mcEvent(0),
  trackCutsV0(),
  fHistdEdx(0),
  fHistdEdxV0(0),
  fHistNumEvents(0),
  fHistTrigger(0), 
  fTree(0),
  hTree(0),
  fTreeGen(0),
  fHistogramList(NULL),
  fPrimaryVertex(),
  fPIDCheckOnly(kFALSE),
  fMCtrue(0),
  fEventCuts(),
  fBetheParamsHe(),
  fBetheParamsT(),
  MB(0),
  HMV0(0),
  HMSPD(0),
  HNU(0),
  HQU(0),
  fNumberV0s(-99),
  fonTheFly(-99),
  frunnumber(-99),
  fTrigMB(-99),			
  fTrigHMV0(-99),
  fTrigHMSPD(-99),
  fTrigHNU(0),
  fTrigHQU(0),
  fNTracks(0),
  fMultV0M(-99),
  fMultOfV0M(-99),
  fMultSPDTracklet(-99),
  fMultSPDCluster(-99),
  fMultRef05(-99),
  fMultRef08(-99),
  fSPDCluster(-99),
  fSPDTracklets(-99),
  fSPDFiredChips0(-99),
  fSPDFiredChips1(-99),
  fV0Multiplicity(-99),
  fYear(0),
  fTRDvalidP(0),
  fTRDtrigHNUP(0),
  fTRDtrigHQUP(0),
  fTRDPidP(0),
  fTRDnTrackletsP(0),
  fTRDPtP(0),
  fTRDLayerMaskP(0),
  fTRDSagittaP(-1),
  fTRDvalidPi(0),
  fTRDtrigHNUPi(0),
  fTRDtrigHQUPi(0),
  fTRDPidPi(0),
  fTRDnTrackletsPi(0),
  fTRDPtPi(0),
  fTRDLayerMaskPi(0),
  fTRDSagittaPi(-1),
  fMCtrueP(-99),
  fisPrimaryP(-99),
  fisWeakP(-99),
  fisMaterialP(-99),
  fMCtrueL(-99),
  fisWeakL(-99),
  fisMaterialL(-99),
  fisPrimaryL(-99),
  fCharge(-99),
  fLambdaM(-99),
  fLambdaP(-99),
  fLambdaPx(-99),
  fLambdaPy(-99),
  fLambdaPz(-99),
  fLambdaPt(-99),
  fLambdaE(-99),
  fLambdaCt(-99),
  fLambdaY(-99),
  fLambdaDca(-99),
  fLambdaCos(-99),
  fpiP(-99),
  fpiPt(-99),
  fpiPx(-99),
  fpiPy(-99),
  fpiPz(-99),
  fpiE(-99),
  fpiy(-99),
  fpiDedx(-99),
  fpiDedxSigma(-99),
  fpichi2(-99),
  fpiNcls(-99),
  fpiNclsITS(-99),  
  fpiDca(-99),
  fpiDcaz(-99),
  fpiDcaSec(-99),
  fGeoLengthPi(-99),  
  fTOFSignalPi(-99),
  fTOFnSigmaPi(-99),
  fEtaPi(-99),
  fPhiPi(-99),
  fKinkPi(-1),
  fTPCrefitPi(-1),
  fITSrefitPi(-1),
  fITSchi2Pi(-99),
  fpiSigmaYX(-99),
  fpiSigmaXYZ(-99),
  fpiSigmaZ(-99),
  fpP(-99),
  fpy(-99),
  fpPt(-99),
  fpPx(-99),
  fpPy(-99),
  fpPz(-99),
  fpE(-99),
  fpDedx(-99),
  fpDedxSigma(-99),
  fpchi2(-99),
  fpNcls(-99),
  fpNclsITS(-99),
  fpDcaSec(-99),
  fpDca(-99),
  fpDcaz(-99),
  fGeoLengthP(-99),
  fTOFSignalP(-99),
  fTOFnSigmaP(-99),
  fEtaP(-99),
  fPhiP(-99),
  fKinkP(-1),
  fTPCrefitP(-1),
  fITSrefitP(-1),
  fITSchi2P(-99),
  fpSigmaYX(-99),
  fpSigmaXYZ(-99),
  fpSigmaZ(-99),
  farmalpha(-99),
  farmpt(-99), 
  fthetaP(-99),
  fthetaN(-99)
{
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
  DefineOutput(3, TTree::Class());
  DefineOutput(4, TTree::Class());
}

// Destructor
AliAnalysisTaskS3ParticleYields::~AliAnalysisTaskS3ParticleYields() {
  if(fHistogramList) delete fHistogramList;
  ResetVals("Event");    
}
const Int_t AliAnalysisTaskS3ParticleYields::fgkPdgCode[] = {
  211,                //PionPlus
  -211,               //PionMinus
  2212,               //Proton
  -2212,              //Anti-Proton
  -321,               //KaonMinus
  321,                //KaonPlus
  1000010020,         //Deuteron
  -1000010020,        //Anti-Deuteron
  1000010030,         //Triton
  -1000010030,        //Anti-Triton
  1000020030,         //Helium3
  -1000020030,        //Anti-Helium3
  1000020040,         //Helium4
  -1000020040,        //Anti-Helium4
  3122,               //Lambda
  -3122,              //Anti-Lambda
  1114,               //DeltaMinus
  3334,               //OmegaMinus
  -3334,              //OmegaPlus
  3312,               //XiMinus
  -3312,              //XiPlus
  1060020020,         //OmegaOmega
  -1060020020,        //AntiOmegaOmega
  1010000030,         //LambdaNeutronNeutron
  -1010000030,        //Anti-Lambda-Neutron-Neutron
  1030000020,         //Xi0-proton
  -1030000020,        //Anti-Xi0-proton
  1010000020,         //LambdaN
  -1010000020,        //AntiLambdaN
  1030000020,         //OmegaProton
  -1030000020,        //AntiOmegaProton
  1020000021,         //LambdaLambda
  -1020000021,        //AntiLambdaLambda
  1010020040,         //HyperHelium4
  -1010020040,        //AntiHyperHelium 4
  1010020050,         //HyperHelium5
  -1010020050,        //AntiHyperHelium 5
  1020010040,         //DoubleHyperHydrogen 4
  -1020010040,        //AntiDoubleHyperHydrogen 4
  1010010030,         //HyperHydrogen3
  -1010010030,        //AntiHyperHydrogen3
  1010010040,         //HyperHydrogen4
  -1010010040,        //AntiHyperHydrogen4
};

void AliAnalysisTaskS3ParticleYields::UserCreateOutputObjects() {

  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  if (man) {
    fInputHandler = dynamic_cast<AliESDInputHandler*> (man->GetInputEventHandler());
    if (fInputHandler)   fPID = fInputHandler->GetESDpid();
  }
  
  fHistdEdx = new TH2F("fHistdEdX","dE/dx;#frac{#it{p}}{z} (GeV/#it{c});TPC Signal (a.u.)",1000,-5.0,5.0,1000,0.0,1500);
  fHistdEdxV0 = new TH2F("fHistdEdXV0","dE/dx;#frac{#it{p}}{z} (GeV/#it{c});TPC Signal (a.u.)",1000,-5.0,5.0,1000,0.0,1500);
    
  fHistNumEvents = new TH1F("fHistNumEvents","Number of Events",2,0,2);
  fHistNumEvents->GetXaxis()->SetBinLabel(1,"before PhysSel");
  fHistNumEvents->GetXaxis()->SetBinLabel(2,"after PhysSel");
    
  fHistTrigger = new TH1F("fHistTrigger","Trigger",7,0,7);
  fHistTrigger->GetXaxis()->SetBinLabel(1,"Events");
  fHistTrigger->GetXaxis()->SetBinLabel(2,"kINT7");
  fHistTrigger->GetXaxis()->SetBinLabel(3,"kHighMultV0");
  fHistTrigger->GetXaxis()->SetBinLabel(4,"kHighMultSPD");
  fHistTrigger->GetXaxis()->SetBinLabel(5,"HNU");
  fHistTrigger->GetXaxis()->SetBinLabel(6,"HQU");
  fHistTrigger->GetXaxis()->SetBinLabel(7,"HJT");
     
  fHistogramList = new TList();
  fHistogramList->SetOwner(kTRUE);
  fHistogramList->SetName(GetName());
  fHistogramList->Add(fHistdEdx);
  fHistogramList->Add(fHistdEdxV0);
  fHistogramList->Add(fHistNumEvents);
  fHistogramList->Add(fHistTrigger);
  fEventCuts.AddQAplotsToList(fHistogramList);
  
  fTree = new TTree("treeP", "fTree");
  fTree->Branch("frunnumber",      &frunnumber,      "frunnumber/I");
  fTree->Branch("fNTracks",        &fNTracks,        "fNTracks/I");
  fTree->Branch("fTrigMB",         &fTrigMB,         "fTrigMB/I");
  fTree->Branch("fTrigHMV0",       &fTrigHMV0,       "fTrigHMV0/I");
  fTree->Branch("fTrigHMSPD",      &fTrigHMSPD,      "fTrigHMSPD/I");
  fTree->Branch("fTrigHNU",        &fTrigHNU,        "fTrigHNU/I");
  fTree->Branch("fTrigHQU",        &fTrigHQU,        "fTrigHQU/I");	
  fTree->Branch("fMultV0M",        &fMultV0M,        "fMultV0M/I");
  fTree->Branch("fMultOfV0M",      &fMultOfV0M,      "fMultOfV0M/I");
  fTree->Branch("fMultSPDTracklet",&fMultSPDTracklet,"fMultSPDTracklet/I");
  fTree->Branch("fMultSPDCluster", &fMultSPDCluster, "fMultSPDCluster/I");
  fTree->Branch("fMultRef05",      &fMultRef05,      "fMultRef05/I");
  fTree->Branch("fMultRef08",      &fMultRef08,      "fMultRef08/I");
  fTree->Branch("fSPDCluster",     &fSPDCluster,     "fSPDCluster/I");
  fTree->Branch("fSPDTracklets",   &fSPDTracklets,   "fSPDTracklets/I");
  fTree->Branch("fSPDFiredChips0", &fSPDFiredChips0, "fSPDFiredChips0/I");
  fTree->Branch("fSPDFiredChips1", &fSPDFiredChips1, "fSPDFiredChips1/I");
  fTree->Branch("fV0Multiplicity", &fV0Multiplicity, "fV0Multiplicity/I");
  fTree->Branch("fTRDvalidP",      &fTRDvalidP,      "fTRDvalidP/I");
  fTree->Branch("fTRDtrigHNUP",    &fTRDtrigHNUP,    "fTRDtrigHNUP/I");
  fTree->Branch("fTRDtrigHQUP",    &fTRDtrigHQUP,    "fTRDtrigHQUP/I");
  fTree->Branch("fTRDPidP",        &fTRDPidP,        "fTRDPidP/I");
  fTree->Branch("fTRDnTrackletsP", &fTRDnTrackletsP, "fTRDnTrackletsP/I");
  fTree->Branch("fTRDPtP",         &fTRDPtP,         "fTRDPtP/I");
  fTree->Branch("fTRDLayerMaskP",  &fTRDLayerMaskP,  "fTRDLayerMaskP/I");
  fTree->Branch("fTRDSagittaP",    &fTRDSagittaP,    "fTRDSagittaP/F");
  fTree->Branch("fCharge",         &fCharge,         "fCharge/I"); 
  fTree->Branch("fpP",             &fpP,             "fpP/F");
  fTree->Branch("fpE",             &fpE,             "fpE/F");
  fTree->Branch("fpPt",            &fpPt,            "fpPt/F");
  fTree->Branch("fpPx",            &fpPx,            "fpPx/F");
  fTree->Branch("fpPy",            &fpPy,            "fpPy/F");
  fTree->Branch("fpPz",            &fpPz,            "fpPz/F");
  fTree->Branch("fpy",             &fpy,             "fpy/F");  
  fTree->Branch("fpDca",           &fpDca,           "fpDca/F");
  fTree->Branch("fpDcaz",          &fpDcaz,          "fpDcaz/F");
  fTree->Branch("fpSigmaYX",       &fpSigmaYX,       "fpSigmaYX/F");
  fTree->Branch("fpSigmaXYZ",      &fpSigmaXYZ,      "fpSigmaXYZ/F");
  fTree->Branch("fpSigmaZ",        &fpSigmaZ,        "fpSigmaZ/F");
  fTree->Branch("fpchi2",          &fpchi2,          "fpchi2/F");
  fTree->Branch("fpNcls",          &fpNcls,          "fpNcls/I");
  fTree->Branch("fpNclsITS",       &fpNclsITS,       "fpNclsITS/I");
  fTree->Branch("fpDedx",          &fpDedx,          "fpDedx/F");
  fTree->Branch("fpDedxSigma",     &fpDedxSigma,     "fpDedxSigma/F");  
  fTree->Branch("fEtaP",           &fEtaP,           "fEtaP/F");
  fTree->Branch("fPhiP",           &fPhiP,           "fPhiP/F");
  fTree->Branch("fGeoLengthP",     &fGeoLengthP,     "fGeoLengthP/F");
  fTree->Branch("fTOFSignalP",     &fTOFSignalP,     "fTOFSignalP/F");
  fTree->Branch("fTOFnSigmaP",     &fTOFnSigmaP,     "fTOFnSigmaP/F");
  fTree->Branch("fKinkP",          &fKinkP,          "fKinkP/I");
  fTree->Branch("fTPCrefitP",      &fTPCrefitP,      "fTPCrefitP/I");
  fTree->Branch("fITSrefitP",      &fITSrefitP,      "fITSrefitP/I");
  fTree->Branch("fITSchi2P",       &fITSchi2P,       "fITSchi2P/F");
  fTree->Branch("fMCtrueP",        &fMCtrueP,        "fMCtrueP/I");
  fTree->Branch("fisPrimaryP",     &fisPrimaryP,     "fisPrimaryP/I");
  fTree->Branch("fisWeakP",        &fisWeakP,        "fisWeakP/I");
  fTree->Branch("fisMaterialP",    &fisMaterialP,    "fisMaterialP/I");  
  
  hTree = new TTree("treeL", "hTree");
  hTree->Branch("fonTheFly",       &fonTheFly,       "fonTheFly/I");
  hTree->Branch("frunnumber",      &frunnumber,      "frunnumber/I");
  hTree->Branch("fNumberV0s",      &fNumberV0s,      "fNumberV0s/I");
  hTree->Branch("fTrigMB",         &fTrigMB,         "fTrigMB/I");
  hTree->Branch("fTrigHMV0",       &fTrigHMV0,       "fTrigHMV0/I");
  hTree->Branch("fTrigHMSPD",      &fTrigHMSPD,      "fTrigHMSPD/I");
  hTree->Branch("fTrigHNU",        &fTrigHNU,        "fTrigHNU/I");
  hTree->Branch("fTrigHQU",        &fTrigHQU,        "fTrigHQU/I");  
  hTree->Branch("fMultV0M",        &fMultV0M,        "fMultV0M/I");
  hTree->Branch("fMultOfV0M",      &fMultOfV0M,      "fMultOfV0M/I");
  hTree->Branch("fMultSPDTracklet",&fMultSPDTracklet,"fMultSPDTracklet/I");
  hTree->Branch("fMultSPDCluster", &fMultSPDCluster, "fMultSPDCluster/I");
  hTree->Branch("fMultRef05",      &fMultRef05,      "fMultRef05/I");
  hTree->Branch("fMultRef08",      &fMultRef08,      "fMultRef08/I");
  hTree->Branch("fSPDCluster",     &fSPDCluster,     "fSPDCluster/I");
  hTree->Branch("fSPDTracklets",   &fSPDTracklets,   "fSPDTracklets/I");
  hTree->Branch("fSPDFiredChips0", &fSPDFiredChips0, "fSPDFiredChips0/I");
  hTree->Branch("fSPDFiredChips1", &fSPDFiredChips1, "fSPDFiredChips1/I");
  hTree->Branch("fV0Multiplicity", &fV0Multiplicity, "fV0Multiplicity/I");
  hTree->Branch("fTRDvalidP",      &fTRDvalidP,      "fTRDvalidP/I");
  hTree->Branch("fTRDtrigHNUP",    &fTRDtrigHNUP,    "fTRDtrigHNUP/I");
  hTree->Branch("fTRDtrigHQUP",    &fTRDtrigHQUP,    "fTRDtrigHQUP/I");
  hTree->Branch("fTRDPidP",        &fTRDPidP,        "fTRDPidP/I");
  hTree->Branch("fTRDnTrackletsP", &fTRDnTrackletsP, "fTRDnTrackletsP/I");
  hTree->Branch("fTRDPtP",         &fTRDPtP,         "fTRDPtP/I");
  hTree->Branch("fTRDLayerMaskP",  &fTRDLayerMaskP,  "fTRDLayerMaskP/I");
  hTree->Branch("fTRDSagittaP",    &fTRDSagittaP,    "fTRDSagittaP/F");
  hTree->Branch("fTRDvalidPi",     &fTRDvalidPi,     "fTRDvalidPi/I");
  hTree->Branch("fTRDtrigHNUPi",   &fTRDtrigHNUPi,   "fTRDtrigHNUPi/I");
  hTree->Branch("fTRDtrigHQUPi",   &fTRDtrigHQUPi,   "fTRDtrigHQUPi/I");
  hTree->Branch("fTRDPidPi",       &fTRDPidPi,       "fTRDPidPi/I");
  hTree->Branch("fTRDnTrackletsPi",&fTRDnTrackletsPi,"fTRDnTrackletsPi/I");
  hTree->Branch("fTRDPtPi",        &fTRDPtPi,        "fTRDPtPi/I");
  hTree->Branch("fTRDLayerMaskPi", &fTRDLayerMaskPi, "fTRDLayerMaskPi/I");
  hTree->Branch("fTRDSagittaPi",   &fTRDSagittaPi,   "fTRDSagittaPi/F");
  hTree->Branch("fCharge",         &fCharge,         "fCharge/I");
  hTree->Branch("fLambdaM",        &fLambdaM,        "fLambdaM/F");
  hTree->Branch("fLambdaE",        &fLambdaE,        "fLambdaE/F");
  hTree->Branch("fLambdaP",        &fLambdaP,        "fLambdaP/F");
  hTree->Branch("fLambdaPx",       &fLambdaPx,       "fLambdaPx/F");
  hTree->Branch("fLambdaPy",       &fLambdaPy,       "fLambdaPy/F");
  hTree->Branch("fLambdaPz",       &fLambdaPz,       "fLambdaPz/F");
  hTree->Branch("fLambdaPt",       &fLambdaPt,       "fLambdaPt/F");
  hTree->Branch("fLambdaCt",       &fLambdaCt,       "fLambdaCt/F");
  hTree->Branch("fLambdaDca",      &fLambdaDca,      "fLambdaDca/F");
  hTree->Branch("fLambdaCos",      &fLambdaCos,      "fLambdaCos/F");
  hTree->Branch("fLambdaY",        &fLambdaY,        "fLambdaY/F");
  hTree->Branch("farmalpha",       &farmalpha,       "farmalpha/F");
  hTree->Branch("farmpt",          &farmpt,          "farmpt/F");
  hTree->Branch("fpP",             &fpP,             "fpP/F");
  hTree->Branch("fpPt",            &fpPt,            "fpPt/F");
  hTree->Branch("fpPx",            &fpPx,            "fpPx/F");
  hTree->Branch("fpPy",            &fpPy,            "fpPy/F");
  hTree->Branch("fpPz",            &fpPz,            "fpPz/F");
  hTree->Branch("fpE",             &fpE,             "fpE/F");
  hTree->Branch("fpy",             &fpy,             "fpy/F");
  hTree->Branch("fpSigmaYX",       &fpSigmaYX,       "fpSigmaYX/F");
  hTree->Branch("fpSigmaXYZ",      &fpSigmaXYZ,      "fpSigmaXYZ/F");
  hTree->Branch("fpSigmaZ",        &fpSigmaZ,        "fpSigmaZ/F");
  hTree->Branch("fpDca",           &fpDca,           "fpDca/F");
  hTree->Branch("fpDcaz",          &fpDcaz,          "fpDcaz/F");
  hTree->Branch("fpDcaSec",        &fpDcaSec,        "fpDcaSec/F");
  hTree->Branch("fpchi2",          &fpchi2,          "fpchi2/F");
  hTree->Branch("fpNcls",          &fpNcls,          "fpNcls/I");
  hTree->Branch("fpNclsITS",       &fpNclsITS,       "fpNclsITS/I");
  hTree->Branch("fpDedx",          &fpDedx,          "fpDedx/F");
  hTree->Branch("fpDedxSigma",     &fpDedxSigma,     "fpDedxSigma/F");  
  hTree->Branch("fEtaP",           &fEtaP,           "fEtaP/F");
  hTree->Branch("fPhiP",           &fPhiP,           "fPhiP/F");
  hTree->Branch("fGeoLengthP",     &fGeoLengthP,     "fGeoLengthP/F");
  hTree->Branch("fTOFSignalP",     &fTOFSignalP,     "fTOFSignalP/F");
  hTree->Branch("fTOFnSigmaP",     &fTOFnSigmaP,     "fTOFnSigmaP/F");
  hTree->Branch("fKinkP",          &fKinkP,          "fKinkP/I");
  hTree->Branch("fTPCrefitP",      &fTPCrefitP,      "fTPCrefitP/I");
  hTree->Branch("fITSrefitP",      &fITSrefitP,      "fITSrefitP/I");
  hTree->Branch("fITSchi2P",       &fITSchi2P,       "fITSchi2P/F");
  hTree->Branch("fpiP",            &fpiP,            "fpiP/F");
  hTree->Branch("fpiPt",           &fpiPt,           "fpiPt/F");
  hTree->Branch("fpiPx",           &fpiPx,           "fpiPx/F");
  hTree->Branch("fpiPy",           &fpiPy,           "fpiPy/F");
  hTree->Branch("fpiPz",           &fpiPz,           "fpiPz/F");
  hTree->Branch("fpiE",            &fpiE,            "fpiE/F");
  hTree->Branch("fpiy",            &fpiy,            "fpiy/F");
  hTree->Branch("fpiSigmaYX",      &fpiSigmaYX,      "fpiSigmaYX/F");
  hTree->Branch("fpiSigmaXYZ",     &fpiSigmaXYZ,     "fpiSigmaXYZ/F");
  hTree->Branch("fpiSigmaZ",       &fpiSigmaZ,       "fpiSigmaZ/F");
  hTree->Branch("fpiDca",          &fpiDca,          "fpiDca/F");
  hTree->Branch("fpiDcaz",         &fpiDcaz,         "fpiDcaz/F");
  hTree->Branch("fpiDcaSec",       &fpiDcaSec,       "fpiDcaSec/F");
  hTree->Branch("fpichi2",         &fpichi2,         "fpichi2/F");
  hTree->Branch("fpiNcls",         &fpiNcls,         "fpiNcls/I");
  hTree->Branch("fpiNclsITS",      &fpiNclsITS,      "fpiNclsITS/I");
  hTree->Branch("fpiDedx",         &fpiDedx,         "fpiDedx/F");
  hTree->Branch("fpiDedxSigma",    &fpiDedxSigma,    "fpiDedxSigma/F");  
  hTree->Branch("fEtaPi",          &fEtaPi,          "fEtaPi/F");
  hTree->Branch("fPhiPi",          &fPhiPi,          "fPhiPi/F");
  hTree->Branch("fGeoLengthPi",    &fGeoLengthPi,    "fGeoLengthPi/F");
  hTree->Branch("fTOFSignalPi",    &fTOFSignalPi,    "fTOFSignalPi/F");  
  hTree->Branch("fTOFnSigmaPi",    &fTOFnSigmaPi,    "fTOFnSigmaPi/F");
  hTree->Branch("fKinkPi",         &fKinkPi,         "fKinkPi/I");
  hTree->Branch("fTPCrefitPi",     &fTPCrefitPi,     "fTPCrefitPi/I");
  hTree->Branch("fITSrefitPi",     &fITSrefitPi,     "fITSrefitPi/I");
  hTree->Branch("fITSchi2Pi",      &fITSchi2Pi,      "fITSchi2Pi/F");
  hTree->Branch("fMCtrueL",        &fMCtrueL,        "fMCtrueL/I");  
  hTree->Branch("fisPrimaryL",     &fisPrimaryL,     "fisPrimaryL/I");
  hTree->Branch("fisWeakL",        &fisWeakL,        "fisWeakL/I");
  hTree->Branch("fisMaterialL",    &fisMaterialL,    "fisMaterialL/I");
    
  fTreeGen = new TTree("treeGenLp", "fTreeGen");
  fTreeGen->Branch("frunnumber",       &frunnumber,       "frunnumber/I");
  fTreeGen->Branch("fTrigMB",          &fTrigMB,          "fTrigMB/I");
  fTreeGen->Branch("fTrigHMV0",        &fTrigHMV0,        "fTrigHMV0/I");
  fTreeGen->Branch("fTrigHMSPD",       &fTrigHMSPD,       "fTrigHMSPD/I");
  fTreeGen->Branch("fTrigHNU",         &fTrigHNU,         "fTrigHNU/I");
  fTreeGen->Branch("fTrigHQU",         &fTrigHQU,         "fTrigHQU/I");
  fTreeGen->Branch("fMultV0M",         &fMultV0M,         "fMultV0M/I");
  fTreeGen->Branch("fMultOfV0M",       &fMultOfV0M,       "fMultOfV0M/I");
  fTreeGen->Branch("fMultSPDTracklet", &fMultSPDTracklet, "fMultSPDTracklet/I");
  fTreeGen->Branch("fMultSPDCluster",  &fMultSPDCluster,  "fMultSPDCluster/I");
  fTreeGen->Branch("fMultRef05",       &fMultRef05,       "fMultRef05/I");
  fTreeGen->Branch("fMultRef08",       &fMultRef08,       "fMultRef08/I");
  fTreeGen->Branch("fSPDCluster",      &fSPDCluster,      "fSPDCluster/I");
  fTreeGen->Branch("fSPDTracklets",    &fSPDTracklets,    "fSPDTracklets/I");
  fTreeGen->Branch("fSPDFiredChips0",  &fSPDFiredChips0,  "fSPDFiredChips0/I");
  fTreeGen->Branch("fSPDFiredChips1",  &fSPDFiredChips1,  "fSPDFiredChips1/I");
  fTreeGen->Branch("fV0Multiplicity",  &fV0Multiplicity,  "fV0Multiplicity/I");
  fTreeGen->Branch("fCharge",          &fCharge,          "fCharge/I");
  fTreeGen->Branch("fMCtrueP",         &fMCtrueP,         "fMCtrueP/I");
  fTreeGen->Branch("fpP",              &fpP,              "fpP/F");
  fTreeGen->Branch("fpPt",             &fpPt,             "fpPt/F");  
  fTreeGen->Branch("fpy",              &fpy,              "fpy/F");
  fTreeGen->Branch("fisPrimaryP",      &fisPrimaryP,      "fisPrimaryP/I");
  fTreeGen->Branch("fisWeakP",         &fisWeakP,         "fisWeakP/I");
  fTreeGen->Branch("fisMaterialP",     &fisMaterialP,     "fisMaterialP/I");
  fTreeGen->Branch("fLambdaP",         &fLambdaP,         "fLambdaP/F");
  fTreeGen->Branch("fLambdaPt",        &fLambdaPt,        "fLambdaPt/F");
  fTreeGen->Branch("fLambdaY",         &fLambdaY,         "fLambdaY/F");
  fTreeGen->Branch("fLambdaM",         &fLambdaM,         "fLambdaM/F");
  fTreeGen->Branch("fMCtrueL",         &fMCtrueL,         "fMCtrueL/I");
  fTreeGen->Branch("fisPrimaryL",      &fisPrimaryL,      "fisPrimaryL/I");
  fTreeGen->Branch("fisWeakL",         &fisWeakL,         "fisWeakL/I");
  fTreeGen->Branch("fisMaterialL",     &fisMaterialL,     "fisMaterialL/I");
    
  PostData(1, fHistogramList);
  PostData(2, fTree);
  PostData(3, hTree);
  PostData(4, fTreeGen);    
}

void AliAnalysisTaskS3ParticleYields::UserExec(Option_t *) {
  
  // MC
  fMCtrue = kTRUE;
  AliMCEventHandler* mcEventHandler = dynamic_cast<AliMCEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()); 
  if (!mcEventHandler) fMCtrue = kFALSE;
  if (mcEventHandler) mcEvent = mcEventHandler->MCEvent();
  if (!mcEvent && fMCtrue) return;
  if (fMCtrue) {
    fStack = mcEvent->Stack();
    if (!fStack) return;
  }
  fESDevent = dynamic_cast<AliESDEvent*>(InputEvent());
  fHistNumEvents->Fill(0);
  fEventCuts.OverrideAutomaticTriggerSelection(AliVEvent::kINT7 | AliVEvent::kTRD | AliVEvent::kHighMultV0 | AliVEvent::kHighMultSPD);
  if(!fEventCuts.AcceptEvent(fESDevent)) {
    PostData(1,fHistogramList);
    return;
  }
  
  frunnumber = fESDevent->GetRunNumber();
  SetBetheBlochParams(frunnumber);
  AliCDBManager *cdbMgr = AliCDBManager::Instance();
  if (fMCtrue) {
    cdbMgr->SetDefaultStorage("MC","Full");
  }
  else {
    cdbMgr->SetDefaultStorage (Form("alien://Folder=/alice/data/%d/OCDB", fYear));
  }
  cdbMgr->SetRun(frunnumber);
  AliGeomManager::LoadGeometry();
  TriggerSelection();
  SetMultiplicity();

  fHistNumEvents->Fill(1);
  const AliESDVertex *vertex = fESDevent->GetPrimaryVertexSPD();
  fPrimaryVertex.SetXYZ(vertex->GetX(),vertex->GetY(),vertex->GetZ());

  trackCutsV0 = new AliESDtrackCuts("AlitrackCutsV0", "AlitrackCutsV0");   
  trackCutsV0->SetEtaRange(-0.8,0.8);
  trackCutsV0->SetAcceptKinkDaughters(kTRUE);
  trackCutsV0->SetRequireTPCRefit(kFALSE);
  trackCutsV0->SetMaxChi2PerClusterTPC(6);
  trackCutsV0->SetMinNClustersTPC(60);
    
  if(fPIDCheckOnly){
    dEdxCheck();
  }
  else{
    ProtonTracks();
    LambdaV0s();
    if(fMCtrue) MCGenerated();
  }    
  PostData(1, fHistogramList);
  PostData(2, fTree);
  PostData(3, hTree);
  PostData(4, fTreeGen);
}
void AliAnalysisTaskS3ParticleYields::dEdxCheck(){
  AliESDtrackCuts* trackCutsPid = new AliESDtrackCuts("trackCutsPid", "trackCutsPid");
  trackCutsPid = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
  trackCutsPid->SetEtaRange(-0.9,0.9);
  for (Int_t itrack = 0; itrack < fESDevent->GetNumberOfTracks(); itrack++) {
    AliESDtrack* track = fESDevent->GetTrack(itrack);
    if (!trackCutsPid->AcceptTrack(track)) continue;
    Double_t momentum = track->GetInnerParam()->GetP();
    fHistdEdx->Fill(momentum * track->GetSign(), track->GetTPCsignal());
  }
  delete trackCutsPid;
}
void AliAnalysisTaskS3ParticleYields::ProtonTracks(){  
  for (Int_t ATracks = 0; ATracks < fESDevent->GetNumberOfTracks(); ATracks++) {
    
    AliESDtrack* trackA = dynamic_cast<AliESDtrack*>(fESDevent->GetTrack(ATracks));
        
    if (!trackA->GetInnerParam()) continue;
        
    if (!trackCutsV0->AcceptTrack(trackA)) continue;
        
    Double_t ptotA = trackA->GetInnerParam()->GetP();
    Double_t signA = trackA->GetSign();
        
    fHistdEdx->Fill(ptotA*signA, trackA->GetTPCsignal());
    Float_t xv[2],  yv[3];
    if(TMath::Abs(xv[0]) > 3.0) continue;

    TLorentzVector lorentzV(0.,0.,0.,0.);
    lorentzV.SetXYZM(trackA->Px(), trackA->Py(), trackA->Pz(), AliPID::ParticleMass(AliPID::kProton));
    if (TMath::Abs(fPID->NumberOfSigmasTPC(trackA, AliPID::kProton)) < 4){
      if((lorentzV.Pt() > 0.5 && lorentzV.Pt() < 2.0) && (TOFSignal(*trackA) > 1.1 || TOFSignal(*trackA) < 0.8)) continue;
      
      fTrigMB = MB;		
      fTrigHMV0 = HMV0;
      fTrigHMSPD = HMSPD;
      fTrigHNU = HNU;
      fTrigHQU = HQU;

      fKinkP = trackA->GetKinkIndex(0) > 0;
      fTPCrefitP = (trackA->GetStatus() & AliESDtrack::kTPCrefit) != 0;
      fITSrefitP = (trackA->GetStatus() & AliESDtrack::kITSrefit) != 0;
      fCharge = signA;
           
      fpPx = lorentzV.Px();
      fpPy = lorentzV.Py();
      fpPz = lorentzV.Pz();
      fpE  = lorentzV.E();
      fpPt = lorentzV.Pt();
      fpy  = lorentzV.Rapidity();

      fpP = trackA->GetInnerParam()->GetP();      
      fpNcls = trackA->GetTPCNcls();
      fpNclsITS = trackA->GetNumberOfITSClusters();
      fpchi2 = trackA->GetTPCchi2() / (Float_t) trackA->GetTPCclusters(0);
      fITSchi2P = trackA->GetITSchi2() / (Float_t) fpNclsITS;
      fpDedx = trackA->GetTPCsignal();
      fpDedxSigma = fPID->NumberOfSigmasTPC(trackA, AliPID::kProton);//Bethe(*trackA, AliPID::ParticleMass(AliPID::kProton), 1, fBetheParamsT);      
      fPhiP = trackA->Phi();
      fEtaP = trackA->Eta();
      fTOFSignalP = TOFSignal(*trackA);
      fTOFnSigmaP = fPID->NumberOfSigmasTOF(trackA, AliPID::kProton);
      fGeoLengthP = GeoLength(*trackA);
      
      trackA->GetImpactParameters(xv,yv);
      fpDca = xv[0];
      fpDcaz = xv[1];
      fpSigmaYX = yv[0];
      fpSigmaXYZ = yv[1];
      fpSigmaZ = yv[2];

      TRDtrack(trackA, 1);
            
      //MC
      if (fMCtrue) {

	Int_t label = trackA->GetLabel();
	AliMCParticle *particleMother = (AliMCParticle*)(mcEvent->GetTrack(TMath::Abs(label)));

	fMCtrueP = TMath::Abs(particleMother->PdgCode()) == fgkPdgCode[kPDGProton];
	fisPrimaryP = mcEvent->IsPhysicalPrimary(TMath::Abs(label));
	fisWeakP = mcEvent->IsSecondaryFromWeakDecay(TMath::Abs(label));
	fisMaterialP = mcEvent->IsSecondaryFromMaterial(TMath::Abs(label));
      }
      
      fTree->Fill();
      ResetVals("");
    }
  }
}
void AliAnalysisTaskS3ParticleYields::LambdaV0s(){
  Float_t xv[2],  yv[3];
  fNumberV0s = fESDevent->GetNumberOfV0s();
  for (Int_t ivertex = 0; ivertex < fESDevent->GetNumberOfV0s(); ivertex++) {    

    //Get V0
    fV0 = fESDevent->GetV0(ivertex);

    //check if on the fly
    if(!fV0->GetOnFlyStatus()) continue;

    //ChargeCorrection
    Bool_t v0ChargeCorrect = kTRUE;        
    AliESDtrack* trackN = fESDevent->GetTrack(fV0->GetIndex(0));//pos Track
    AliESDtrack* trackP = fESDevent->GetTrack(fV0->GetIndex(1));//neg Track
    if (trackN->GetSign() > 0 ) {
      trackN = fESDevent->GetTrack(fV0->GetIndex(1));//neg Track
      trackP = fESDevent->GetTrack(fV0->GetIndex(0));//pos Track
      v0ChargeCorrect = kFALSE;
    }

    //TrackCuts
    if (!trackCutsV0->AcceptTrack(trackN)) continue;
    if (!trackCutsV0->AcceptTrack(trackP)) continue;

    //dEdx
    fHistdEdxV0->Fill(trackP->GetInnerParam()->GetP() * trackP->GetSign(), trackP->GetTPCsignal());
    fHistdEdxV0->Fill(trackN->GetInnerParam()->GetP() * trackN->GetSign(), trackN->GetTPCsignal());

    //
    Bool_t pionPositive = kFALSE;
    Bool_t pionNegative = kFALSE;
    Bool_t protonPositive = kFALSE;
    Bool_t protonNegative = kFALSE;
    //pion
    if (TMath::Abs(fPID->NumberOfSigmasTPC(trackP, AliPID::kPion)) < 4) {
      pionPositive = kTRUE;
    }
    if (TMath::Abs(fPID->NumberOfSigmasTPC(trackN, AliPID::kPion)) < 4) {
      pionNegative = kTRUE;
    }
    //proton
    if (TMath::Abs(fPID->NumberOfSigmasTPC(trackP, AliPID::kProton)) < 4) {
      protonPositive = kTRUE;
    }
    if (TMath::Abs(fPID->NumberOfSigmasTPC(trackN, AliPID::kProton)) < 4) {
      protonNegative = kTRUE;
    } 
  
    if(!(pionNegative && protonPositive) && !(pionPositive && protonNegative))continue;

    //Lambda
    if (protonPositive && pionNegative) {// 

      fTrigMB = MB;		
      fTrigHMV0 = HMV0;
      fTrigHMSPD = HMSPD;
      fTrigHNU = HNU;
      fTrigHQU = HQU;

      //Vectors
      TLorentzVector proton(0.,0.,0.,0.);
      proton.SetXYZM(trackP->Px(), trackP->Py(), trackP->Pz(), AliPID::ParticleMass(AliPID::kProton));
      TLorentzVector pion(0.,0.,0.,0.);
      pion.SetXYZM(trackN->Px(), trackN->Py(), trackN->Pz(), AliPID::ParticleMass(AliPID::kPion));      
      TVector3 secVertex(fV0->Xv(), fV0->Yv(), fV0->Zv());            

      //Proton Track
      fKinkP = trackP->GetKinkIndex(0) > 0;
      fTPCrefitP = (trackP->GetStatus() & AliESDtrack::kTPCrefit) != 0;
      fITSrefitP = (trackP->GetStatus() & AliESDtrack::kITSrefit) != 0;
      
      fpPx = proton.Px();
      fpPy = proton.Py();
      fpPz = proton.Pz();
      fpE  = proton.E();
      fpPt = proton.Pt();
      fpy  = proton.Rapidity();

      fpP = trackP->GetInnerParam()->GetP();      
      fpNcls = trackP->GetTPCNcls();
      fpNclsITS = trackP->GetNumberOfITSClusters();
      fpchi2 = trackP->GetTPCchi2() / (Float_t) trackP->GetTPCclusters(0);
      fITSchi2P = trackP->GetITSchi2() / (Float_t) fpNclsITS;
      fpDedx = trackP->GetTPCsignal();
      fpDedxSigma = fPID->NumberOfSigmasTPC(trackP, AliPID::kProton);//Bethe(*trackP, AliPID::ParticleMass(AliPID::kProton), 1, fBetheParamsT);
      fpDcaSec = TMath::Abs(trackP->GetD(secVertex.X(), secVertex.Y(), fESDevent->GetMagneticField()));
      fPhiP = trackP->Phi();
      fEtaP = trackP->Eta();
      fTOFSignalP = TOFSignal(*trackP);
      fTOFnSigmaP = fPID->NumberOfSigmasTOF(trackP, AliPID::kProton);
      fGeoLengthP = GeoLength(*trackP);
     
      trackP->GetImpactParameters(xv,yv);
      fpDca = xv[0];
      fpDcaz = xv[1];
      fpSigmaYX = yv[0];
      fpSigmaXYZ = yv[1];
      fpSigmaZ = yv[2];

      //Pion Track
      fKinkPi = trackN->GetKinkIndex(0) > 0;
      fTPCrefitPi = (trackN->GetStatus() & AliESDtrack::kTPCrefit) != 0;
      fITSrefitPi = (trackN->GetStatus() & AliESDtrack::kITSrefit) != 0;
      
      fpiPx = pion.Px();
      fpiPy = pion.Py();
      fpiPz = pion.Pz();
      fpiE  = pion.E();
      fpiPt = pion.Pt();
      fpiy  = pion.Rapidity();

      fpiP = trackN->GetInnerParam()->GetP();      
      fpiNcls = trackN->GetTPCNcls();
      fpiNclsITS = trackN->GetNumberOfITSClusters();
      fpichi2 = trackN->GetTPCchi2() / (Float_t) trackN->GetTPCclusters(0);
      fITSchi2Pi = trackN->GetITSchi2() / (Float_t) fpiNclsITS;
      fpiDedx = trackN->GetTPCsignal();
      fpiDedxSigma = fPID->NumberOfSigmasTPC(trackN, AliPID::kPion);
      fpiDcaSec = TMath::Abs(trackN->GetD(secVertex.X(), secVertex.Y(), fESDevent->GetMagneticField()));
      fPhiPi = trackN->Phi();
      fEtaPi = trackN->Eta();
      fTOFSignalPi = TOFSignal(*trackN);
      fTOFnSigmaPi = fPID->NumberOfSigmasTOF(trackN, AliPID::kPion);
      fGeoLengthPi = GeoLength(*trackN);
     
      trackN->GetImpactParameters(xv,yv);
      fpiDca = xv[0];
      fpiDcaz = xv[1];
      fpiSigmaYX = yv[0];
      fpiSigmaXYZ = yv[1];
      fpiSigmaZ = yv[2];     

      //particle = 1; anti particle = -1; (charge = 0)
      fCharge = 1;

      TRDtrack(trackP, 1);
      TRDtrack(trackN, 2);

      //V0 information
      fLambdaDca = fV0->GetDcaV0Daughters();
      fLambdaCos = fV0->GetV0CosineOfPointingAngle();
      if(fV0->GetOnFlyStatus()) fonTheFly = 1;

      //Mother information
      TLorentzVector Lambda(0.,0.,0.,0.);
      Lambda = pion + proton;
      fLambdaM = Lambda.M();
      fLambdaE = Lambda.E();
      fLambdaP = Lambda.P();
      fLambdaPx = Lambda.Px();
      fLambdaPy = Lambda.Py();
      fLambdaPz = Lambda.Pz();
      fLambdaPt = Lambda.Pt();
      fLambdaY = Lambda.Rapidity();

      //Lifetime
      secVertex = secVertex - fPrimaryVertex;
      fLambdaCt = secVertex.Mag() * fLambdaM / fLambdaP;

      //ArmenterosPodolanski
      TVector3 vecN = pion.Vect();
      TVector3 vecP = proton.Vect();
      TVector3 vecM = Lambda.Vect();           
      fthetaP = TMath::ACos((vecP * vecM)/(vecP.Mag() * vecM.Mag()));
      fthetaN = TMath::ACos((vecN * vecM)/(vecN.Mag() * vecM.Mag()));            
      farmalpha = ((vecP.Mag())*TMath::Cos(fthetaP)-(vecN.Mag())*TMath::Cos(fthetaN))/((vecP.Mag())*TMath::Cos(fthetaP)+(vecN.Mag())*TMath::Cos(fthetaN));
      farmpt = vecP.Mag()*sin(fthetaP);

      //MC info
      if(fMCtrue){
	Int_t label = trackP->GetLabel();
	Int_t labelMother = mcEvent->GetLabelOfParticleMother(TMath::Abs(label));
	AliMCParticle *particleMother = (AliMCParticle*)(mcEvent->GetTrack(TMath::Abs(labelMother)));

	Int_t label1 = trackN->GetLabel();
	Int_t labelMother1 = mcEvent->GetLabelOfParticleMother(TMath::Abs(label1));
	AliMCParticle *particleMother1 = (AliMCParticle*)(mcEvent->GetTrack(TMath::Abs(labelMother1)));
	
	if(particleMother->PdgCode() == fgkPdgCode[kPDGLambda] && particleMother1->PdgCode() == fgkPdgCode[kPDGLambda] && labelMother == labelMother1){
	  if(TMath::Abs(label) == TMath::Abs(GetLabel(labelMother, fgkPdgCode[kPDGProton]))		   
	     && TMath::Abs(label1) == TMath::Abs(GetLabel(labelMother, fgkPdgCode[kPDGPionMinus]))) {
	    fMCtrueL = kTRUE;
	    fisPrimaryL = mcEvent->IsPhysicalPrimary(TMath::Abs(labelMother1));
	    fisWeakL = mcEvent->IsSecondaryFromWeakDecay(TMath::Abs(labelMother1));
	    fisMaterialL = mcEvent->IsSecondaryFromMaterial(TMath::Abs(labelMother1));
	  }
	}
      }     
      hTree->Fill();
      ResetVals("");
    }//end Lambda       
    else if (protonNegative && pionPositive) {

      fTrigMB = MB;		
      fTrigHMV0 = HMV0;
      fTrigHMSPD = HMSPD;
      fTrigHNU = HNU;
      fTrigHQU = HQU;

      //Vectors
      TLorentzVector proton(0.,0.,0.,0.);
      proton.SetXYZM(trackN->Px(), trackN->Py(), trackN->Pz(), AliPID::ParticleMass(AliPID::kProton));
      TLorentzVector pion(0.,0.,0.,0.);
      pion.SetXYZM(trackP->Px(), trackP->Py(), trackP->Pz(), AliPID::ParticleMass(AliPID::kPion));      
      TVector3 secVertex(fV0->Xv(), fV0->Yv(), fV0->Zv());            

      //Proton Track
      fKinkP = trackN->GetKinkIndex(0) > 0;
      fTPCrefitP = (trackN->GetStatus() & AliESDtrack::kTPCrefit) != 0;
      fITSrefitP = (trackN->GetStatus() & AliESDtrack::kITSrefit) != 0;
      
      fpPx = proton.Px();
      fpPy = proton.Py();
      fpPz = proton.Pz();
      fpE  = proton.E();
      fpPt = proton.Pt();
      fpy  = proton.Rapidity();

      fpP = trackN->GetInnerParam()->GetP();      
      fpNcls = trackN->GetTPCNcls();
      fpNclsITS = trackN->GetNumberOfITSClusters();
      fpchi2 = trackN->GetTPCchi2() / (Float_t) trackN->GetTPCclusters(0);
      fITSchi2P = trackN->GetITSchi2() / (Float_t) fpNclsITS;
      fpDedx = trackN->GetTPCsignal();
      fpDedxSigma = fPID->NumberOfSigmasTPC(trackN, AliPID::kProton);//Bethe(*trackN, AliPID::ParticleMass(AliPID::kProton), 1, fBetheParamsT);
      fpDcaSec = TMath::Abs(trackN->GetD(secVertex.X(), secVertex.Y(), fESDevent->GetMagneticField()));
      fPhiP = trackN->Phi();
      fEtaP = trackN->Eta();
      fTOFSignalP = TOFSignal(*trackN);
      fTOFnSigmaP = fPID->NumberOfSigmasTOF(trackN, AliPID::kProton);
      fGeoLengthP = GeoLength(*trackN);
     
      trackN->GetImpactParameters(xv,yv);
      fpDca = xv[0];
      fpDcaz = xv[1];
      fpSigmaYX = yv[0];
      fpSigmaXYZ = yv[1];
      fpSigmaZ = yv[2];

      //Pion Track
      fKinkPi = trackP->GetKinkIndex(0) > 0;
      fTPCrefitPi = (trackP->GetStatus() & AliESDtrack::kTPCrefit) != 0;
      fITSrefitPi = (trackP->GetStatus() & AliESDtrack::kITSrefit) != 0;
      
      fpiPx = pion.Px();
      fpiPy = pion.Py();
      fpiPz = pion.Pz();
      fpiE  = pion.E();
      fpiPt = pion.Pt();
      fpiy =  pion.Rapidity();

      fpiP = trackP->GetInnerParam()->GetP();      
      fpiNcls = trackP->GetTPCNcls();
      fpiNclsITS = trackP->GetNumberOfITSClusters();
      fpichi2 = trackP->GetTPCchi2() / (Float_t) trackP->GetTPCclusters(0);
      fITSchi2Pi = trackP->GetITSchi2() / (Float_t) fpiNclsITS;
      fpiDedx = trackP->GetTPCsignal();
      fpiDedxSigma = fPID->NumberOfSigmasTPC(trackP, AliPID::kPion);
      fpiDcaSec = TMath::Abs(trackP->GetD(secVertex.X(), secVertex.Y(), fESDevent->GetMagneticField()));
      fPhiPi = trackP->Phi();
      fEtaPi = trackP->Eta();
      fTOFSignalPi = TOFSignal(*trackP);
      fTOFnSigmaPi = fPID->NumberOfSigmasTOF(trackP, AliPID::kPion);
      fGeoLengthPi = GeoLength(*trackP);
     
      trackP->GetImpactParameters(xv,yv);
      fpiDca = xv[0];
      fpiDcaz = xv[1];
      fpiSigmaYX = yv[0];
      fpiSigmaXYZ = yv[1];
      fpiSigmaZ = yv[2];

      TRDtrack(trackP, 2);
      TRDtrack(trackN, 1);

      //particle = 1; anti particle = -1; (charge = 0)
      fCharge = -1;

      //V0 information
      fLambdaDca = fV0->GetDcaV0Daughters();
      fLambdaCos = fV0->GetV0CosineOfPointingAngle();
      if(fV0->GetOnFlyStatus()) fonTheFly = 1;

      //Mother information
      TLorentzVector Lambda(0.,0.,0.,0.);
      Lambda = pion + proton;
      fLambdaM = Lambda.M();
      fLambdaE = Lambda.E();
      fLambdaP = Lambda.P();
      fLambdaPx = Lambda.Px();
      fLambdaPy = Lambda.Py();
      fLambdaPz = Lambda.Pz();
      fLambdaPt = Lambda.Pt();
      fLambdaY = Lambda.Rapidity();

      //Lifetime
      secVertex = secVertex - fPrimaryVertex;
      fLambdaCt = secVertex.Mag() * fLambdaM / fLambdaP;

      //ArmenterosPodolanski
      TVector3 vecP = pion.Vect();
      TVector3 vecN = proton.Vect();
      TVector3 vecM = Lambda.Vect();           
      fthetaP = TMath::ACos((vecP * vecM)/(vecP.Mag() * vecM.Mag()));
      fthetaN = TMath::ACos((vecN * vecM)/(vecN.Mag() * vecM.Mag()));            
      farmalpha = ((vecP.Mag())*TMath::Cos(fthetaP)-(vecN.Mag())*TMath::Cos(fthetaN))/((vecP.Mag())*TMath::Cos(fthetaP)+(vecN.Mag())*TMath::Cos(fthetaN));
      farmpt = vecP.Mag()*sin(fthetaP);

      //MC info
      if(fMCtrue){
	Int_t label = trackN->GetLabel();
	Int_t labelMother = mcEvent->GetLabelOfParticleMother(TMath::Abs(label));
	AliMCParticle *particleMother = (AliMCParticle*)(mcEvent->GetTrack(labelMother));

	Int_t label1 = trackP->GetLabel();
	Int_t labelMother1 = mcEvent->GetLabelOfParticleMother(TMath::Abs(label1));
	AliMCParticle *particleMother1 = (AliMCParticle*)(mcEvent->GetTrack(labelMother1));

	if(particleMother->PdgCode() == fgkPdgCode[kPDGAntiLambda] && particleMother1->PdgCode() == fgkPdgCode[kPDGAntiLambda] && labelMother == labelMother1){
	  if(TMath::Abs(label) == TMath::Abs(GetLabel(labelMother, fgkPdgCode[kPDGAntiProton]))		   
	       && TMath::Abs(label1) == TMath::Abs(GetLabel(labelMother, fgkPdgCode[kPDGPionPlus]))) {
	    fMCtrueL = kTRUE;
	    fisPrimaryL = mcEvent->IsPhysicalPrimary(TMath::Abs(labelMother1));
	    fisWeakL = mcEvent->IsSecondaryFromWeakDecay(TMath::Abs(labelMother1));
	    fisMaterialL = mcEvent->IsSecondaryFromMaterial(TMath::Abs(labelMother1));
	  }
	}
      }     
      hTree->Fill();
      ResetVals("");
    }//end Anti-Lambda
  }
}
void AliAnalysisTaskS3ParticleYields::MCGenerated() {   
  for(Int_t stackN = 0; stackN < mcEvent->GetNumberOfTracks(); stackN++){
    
    AliMCParticle *particleMother = (AliMCParticle*)(mcEvent->GetTrack(stackN));//fStack->Particle(TMath::Abs(stackN));
    TLorentzVector part(0.,0.,0.,0.);
    TLorentzVector partd1(0.,0.,0.,0.);
    TLorentzVector partd2(0.,0.,0.,0.);
    AliMCParticle *tparticleFirstDaughter;
    AliMCParticle *tparticleSecondDaughter;
    Int_t PdgCode = particleMother->PdgCode();
    Int_t charge = PdgCode/TMath::Abs(PdgCode);  

    if(TMath::Abs(PdgCode) != fgkPdgCode[kPDGProton] && TMath::Abs(PdgCode) != fgkPdgCode[kPDGLambda]) continue;
    
    //Proton    
    if(TMath::Abs(PdgCode) == fgkPdgCode[kPDGProton]){

      fTrigMB = MB;		
      fTrigHMV0 = HMV0;
      fTrigHMSPD = HMSPD;
      fTrigHNU = HNU;
      fTrigHQU = HQU;

      fMCtrueP = 1;
      fMCtrueL = 0;

      part.SetXYZM(particleMother->Px(), particleMother->Py(), particleMother->Pz(), AliPID::ParticleMass(AliPID::kProton));

      fpPt = part.Pt();
      fpP = part.P();
      fpy = part.Rapidity();

      fisPrimaryP = mcEvent->IsPhysicalPrimary(TMath::Abs(stackN));
      fisWeakP = mcEvent->IsSecondaryFromWeakDecay(TMath::Abs(stackN));
      fisMaterialP = mcEvent->IsSecondaryFromMaterial(TMath::Abs(stackN));
      fCharge = charge;
      
      fTreeGen->Fill();
      ResetVals("");
    }

    //Lambda
    if(TMath::Abs(PdgCode) == fgkPdgCode[kPDGLambda]) {

      fTrigMB = MB;		
      fTrigHMV0 = HMV0;
      fTrigHMSPD = HMSPD;
      fTrigHNU = HNU;
      fTrigHQU = HQU;

      tparticleFirstDaughter = (AliMCParticle*)(mcEvent->GetTrack(TMath::Abs(GetLabel(stackN, charge*fgkPdgCode[kPDGProton]))));
      tparticleSecondDaughter = (AliMCParticle*)(mcEvent->GetTrack(TMath::Abs(GetLabel(stackN, charge*fgkPdgCode[kPDGPionMinus]))));      

      if(tparticleFirstDaughter->PdgCode() == charge*fgkPdgCode[kPDGProton] && tparticleSecondDaughter->PdgCode() == charge*fgkPdgCode[kPDGPionMinus]){

	partd1.SetXYZM(tparticleFirstDaughter->Px(), tparticleFirstDaughter->Py(), tparticleFirstDaughter->Pz(), AliPID::ParticleMass(AliPID::kProton));
	partd2.SetXYZM(tparticleSecondDaughter->Px(), tparticleSecondDaughter->Py(), tparticleSecondDaughter->Pz(), AliPID::ParticleMass(AliPID::kPion));
	part = partd1 + partd2;

	fMCtrueP = 0;
	fMCtrueL = 1;

	fCharge = charge;
	fLambdaPt = part.Pt();
	fLambdaP = part.P();
	fLambdaY = part.Rapidity();
	fLambdaM = part.M();

	fisPrimaryL = mcEvent->IsPhysicalPrimary(TMath::Abs(stackN));
	fisWeakL = mcEvent->IsSecondaryFromWeakDecay(TMath::Abs(stackN));
	fisMaterialL = mcEvent->IsSecondaryFromMaterial(TMath::Abs(stackN));

	fTreeGen->Fill();
	ResetVals("");
      }
    }    
  }
}
//_____________________________________________________________________________
void AliAnalysisTaskS3ParticleYields::Terminate(const Option_t*) {
  //if (!GetOutputData(0)) return;
}
//_____________________________________________________________________________
/// Calculates number of sigma deviation from expected dE/dx in TPC
/// \param track particle track
/// \param mass mass hypothesis of particle
/// \param charge particle charge hypothesis
/// \param params Parameters of Aleph parametrization of Bethe Energy-loss
Double_t AliAnalysisTaskS3ParticleYields::Bethe(const AliESDtrack& track, Double_t mass, Int_t charge, Double_t* params){
  Double_t expected = charge*charge*AliExternalTrackParam::BetheBlochAleph(charge*track.GetInnerParam()->GetP()/mass,params[0],params[1],params[2],params[3],params[4]);
  Double_t sigma = expected*params[5];
  if (TMath::IsNaN(expected)) return -999;
  return (track.GetTPCsignal() - expected) / sigma;
}
//_____________________________________________________________________________
Double_t AliAnalysisTaskS3ParticleYields::GeoLength(const AliESDtrack& track) {
  Double_t deadZoneWidth = 3.0;
  Double_t lengthInActiveZone = track.GetLengthInActiveZone(1, deadZoneWidth, 220, track.GetESDEvent()->GetMagneticField(),0,0);
  return lengthInActiveZone;
}
//_____________________________________________________________________________
Double_t AliAnalysisTaskS3ParticleYields::TOFSignal(const AliESDtrack& track) {
  Float_t mass = 0, time = -1, beta = 0, gamma = 0, length = 0, time0 = 0;
  length = track.GetIntegratedLength();
  time0 = fPID->GetTOFResponse().GetStartTime(track.P());//fESDpid->GetTOFResponse().GetTimeZero();
  time = track.GetTOFsignal() - time0;
  if (time > 0) {
    beta = length / (2.99792457999999984e-02 * time);
    gamma = 1/TMath::Sqrt(1 - beta*beta);
    mass = (track.GetInnerParam()->GetP())/TMath::Sqrt(gamma*gamma - 1); // using inner TPC mom. as approx.
    return mass;
  }
  return -1;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskS3ParticleYields::TriggerSelection() {
  //******************************
  //*   get trigger information  *
  //******************************  
  MB = 0;
  HMV0 = 0;
  HMSPD = 0;
  HNU = 0;
  HQU = 0;
  if (fInputHandler->IsEventSelected() & AliVEvent::kINT7) MB = kTRUE;
  if (fInputHandler->IsEventSelected() & AliVEvent::kHighMultV0) HMV0 = kTRUE;
  if (fInputHandler->IsEventSelected() & AliVEvent::kHighMultSPD) HMSPD = kTRUE;
	
  Int_t nTrdTracks = fESDevent->GetNumberOfTrdTracks();
  if (!fMCtrue){
    // Data: get TRD trigger information from trigger classes 
    TString classes = fESDevent->GetFiredTriggerClasses();   
    if (classes.Contains("HNU")) HNU = 1;
    if (classes.Contains("HQU")) HQU = 1; 
		
  } else {
    // MC: simulate TRD trigger
    if (nTrdTracks > 0) {
      for (Int_t iTrack = 0; iTrack < nTrdTracks; ++iTrack) {
	AliESDTrdTrack* trdTrack = fESDevent->GetTrdTrack(iTrack);
	if (!trdTrack) continue;					
	// simulate HNU
	if((trdTrack->GetPID() >= 255 && trdTrack->GetNTracklets() == 4) || 
	   (trdTrack->GetPID() >= 235 && trdTrack->GetNTracklets() > 4)) {	
	  HNU = 1;
	}
	// simulate HQU
	if (TMath::Abs(trdTrack->GetPt()) >= 256 &&
	    trdTrack->GetPID() >= 130 && trdTrack->GetNTracklets() >= 5 && (trdTrack->GetLayerMask() & 1) ){	
	  Float_t sag = GetInvPtDevFromBC(trdTrack->GetB(), trdTrack->GetC());
	  if (sag < 0.2 && sag > -0.2) {
	    HQU = 1;
	  }
	}
      }     
    }
  } 
  // fill histogram
  fHistTrigger->Fill(0);
  if (MB) fHistTrigger->Fill(1);
  if (HMV0) fHistTrigger->Fill(2);
  if (HMSPD) fHistTrigger->Fill(3);
  if (HNU) fHistTrigger->Fill(4);
  if (HQU) fHistTrigger->Fill(5);
  Bool_t isTriggered = kFALSE;
  if(MB || HMV0 || HMSPD || HNU || HQU) isTriggered = kTRUE;
  return isTriggered;
}
//_____________________________________________________________________________
Float_t AliAnalysisTaskS3ParticleYields::GetInvPtDevFromBC(Int_t b, Int_t c) {
  //returns d(1/Pt) in c/GeV
  //in case of no gtu simulation -> return maximum 0.5
  if(b==0 && c==0) return 0.5;
  Int_t tmp = (((b & 0xfff) << 12) ^ 0x800000) - 0x800000;
  tmp += (c & 0xfff);
  Float_t invPtDev = tmp * 0.000001;
  return invPtDev;
}
//_____________________________________________________________________________
void AliAnalysisTaskS3ParticleYields::SetMultiplicity() {
  // additional multiplicity information
  fNTracks = fESDevent->GetNumberOfTracks();
	
  AliESDVZERO *vzero = fESDevent->GetVZEROData();
  fV0Multiplicity = 0;
  for (Int_t ii = 0; ii < 64; ii++){
    fV0Multiplicity += vzero->GetMultiplicity(ii);
  }
	
  AliMultiplicity *multSPD = fESDevent->GetMultiplicity();
  fSPDCluster	= multSPD->GetNumberOfSPDClusters();
  fSPDTracklets = multSPD->GetNumberOfTracklets();
  fSPDFiredChips0 = multSPD->GetNumberOfFiredChips(0);
  fSPDFiredChips1 = multSPD->GetNumberOfFiredChips(1);

  AliMultSelection *MultSelection = (AliMultSelection*) fESDevent->FindListObject("MultSelection");
  if (MultSelection) {
    fMultV0M = MultSelection->GetMultiplicityPercentile("V0M");
    fMultOfV0M = MultSelection->GetMultiplicityPercentile("OnlineV0M");
    fMultSPDTracklet = MultSelection->GetMultiplicityPercentile("SPDClusters");
    fMultSPDCluster = MultSelection->GetMultiplicityPercentile("SPDTracklets");
    fMultRef05 = MultSelection->GetMultiplicityPercentile("RefMult05");
    fMultRef08 = MultSelection->GetMultiplicityPercentile("RefMult08");
  }	
}
//________________________________________________________________________
Int_t AliAnalysisTaskS3ParticleYields::GetLabel(Int_t labelFirstMother, Int_t particlePdgCode){
  Int_t labelFirstDaughter = mcEvent->GetLabelOfParticleFirstDaughter(TMath::Abs(labelFirstMother));
  Int_t labelLastDaughter = mcEvent->GetLabelOfParticleLastDaughter(TMath::Abs(labelFirstMother));
  Int_t diff = TMath::Abs(labelLastDaughter - labelFirstDaughter) + 1;
  
  Int_t returnval = -99;
  
  for(Int_t i = 0; i<diff; i++){
      
    Int_t labelDaughter = labelFirstDaughter + i;
    AliMCParticle *Daughter = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelDaughter))->Particle());
      
    if(Daughter->PdgCode()== particlePdgCode){
      returnval = labelDaughter;	 
    }
    if(Daughter) delete Daughter;
  }
    
  return returnval;
}
//_____________________________________________________________________________

Double_t AliAnalysisTaskS3ParticleYields::TRDtrack(AliESDtrack* esdTrack, Int_t particleID) {    

  if(!esdTrack) {
    return 0;
  }

  if(fESDevent->GetNumberOfTrdTracks() == 0) {
    return 0;
  }
    
  AliESDTrdTrack* bestGtuTrack = 0x0;
  AliTRDonlineTrackMatching *matching = new AliTRDonlineTrackMatching();
    
  Double_t esdPt = esdTrack->GetSignedPt();
  Double_t mag = fESDevent->GetMagneticField();
  Double_t currentMatch = 0;
  Double_t bestMatch = 0;

  for (Int_t i = 0; i < fESDevent->GetNumberOfTrdTracks(); i++) {

    AliESDTrdTrack* gtuTrack= fESDevent->GetTrdTrack ( i );
    Double_t gtuPt = gtuTrack->Pt();
    if (mag > 0.) gtuPt = gtuPt * (-1.0);

    Double_t ydist;
    Double_t zdist;

    if (matching->EstimateTrackDistance(esdTrack, gtuTrack, mag, &ydist, &zdist) == 0) {
      currentMatch = matching->RateTrackMatch(ydist, zdist, esdPt, gtuPt);
    }
				
    if (currentMatch > bestMatch) {
      bestMatch = currentMatch;
      bestGtuTrack = gtuTrack;
    }
  }
    
  if (!bestGtuTrack) {
    return 0;
  }

  if(particleID == 1){ // protons

    fTRDvalidP = 0;
    fTRDtrigHNUP = 0;
    fTRDtrigHQUP = 0;
    fTRDPidP = 0;
    fTRDnTrackletsP = 0;
    fTRDPtP = 0;
    fTRDLayerMaskP = 0;
    fTRDSagittaP = -1;
    
    fTRDvalidP = 1;
    fTRDPidP = bestGtuTrack->GetPID();
    fTRDnTrackletsP = bestGtuTrack->GetNTracklets();
    fTRDPtP = (TMath::Abs(bestGtuTrack->GetPt()));
    fTRDLayerMaskP =  bestGtuTrack->GetLayerMask();
    fTRDSagittaP = GetInvPtDevFromBC(bestGtuTrack->GetB(), bestGtuTrack->GetC());	
		
    if((bestGtuTrack->GetPID() >= 255 && bestGtuTrack->GetNTracklets() == 4) || 
       (bestGtuTrack->GetPID() >= 235 && bestGtuTrack->GetNTracklets() > 4)) {
      fTRDtrigHNUP = 1;
    }		

    if (TMath::Abs(bestGtuTrack->GetPt()) >= 256 &&
	bestGtuTrack->GetPID() >= 130 && bestGtuTrack->GetNTracklets() >= 5 && (bestGtuTrack->GetLayerMask() & 1) ){	
      Float_t sag = GetInvPtDevFromBC(bestGtuTrack->GetB(), bestGtuTrack->GetC());
      if (sag < 0.2 && sag > -0.2) {
	fTRDtrigHQUP = 1;
      }
    }	
  }
  if(particleID == 2){ // pions

    fTRDvalidPi = 0;
    fTRDtrigHNUPi = 0;
    fTRDtrigHQUPi = 0;
    fTRDPidPi = 0;
    fTRDnTrackletsPi = 0;
    fTRDPtPi = 0;
    fTRDLayerMaskPi = 0;
    fTRDSagittaPi = -1;
    
    fTRDvalidPi = 1;
    fTRDPidPi = bestGtuTrack->GetPID();
    fTRDnTrackletsPi = bestGtuTrack->GetNTracklets();
    fTRDPtPi = (TMath::Abs(bestGtuTrack->GetPt()));
    fTRDLayerMaskPi =  bestGtuTrack->GetLayerMask();
    fTRDSagittaPi = GetInvPtDevFromBC(bestGtuTrack->GetB(), bestGtuTrack->GetC());	
		
    if((bestGtuTrack->GetPID() >= 255 && bestGtuTrack->GetNTracklets() == 4) || 
       (bestGtuTrack->GetPID() >= 235 && bestGtuTrack->GetNTracklets() > 4)) {
      fTRDtrigHNUPi = 1;
    }		

    if (TMath::Abs(bestGtuTrack->GetPt()) >= 256 &&
	bestGtuTrack->GetPID() >= 130 && bestGtuTrack->GetNTracklets() >= 5 && (bestGtuTrack->GetLayerMask() & 1) ){	
      Float_t sag = GetInvPtDevFromBC(bestGtuTrack->GetB(), bestGtuTrack->GetC());
      if (sag < 0.2 && sag > -0.2) {
	fTRDtrigHQUPi = 1;
      }
    }	
  }
  return bestMatch;

}
//_____________________________________________________________________________
void AliAnalysisTaskS3ParticleYields::SetBetheBlochParams(Int_t runNumber) {
  // set Bethe-Bloch parameter
  if (runNumber >= 252235 && runNumber <= 264347 ) { // 2016 pp
    fYear = 2016;
    if(!fMCtrue) { // Data
      // LHC16 + LHC18
      // He3
      fBetheParamsHe[0] = 1.81085;
      fBetheParamsHe[1] = 29.4656;
      fBetheParamsHe[2] = 0.0458225;
      fBetheParamsHe[3] = 2.08689;
      fBetheParamsHe[4] = 2.28772;
      fBetheParamsHe[5] = 0.06;
      // Triton
      fBetheParamsT[0] = 1.58385;
      fBetheParamsT[1] = 25.8334;
      fBetheParamsT[2] = 0.00908038;
      fBetheParamsT[3] = 2.24769;
      fBetheParamsT[4] = 2.87755;
      fBetheParamsT[5] = 0.06;
    } else { // MC
      if (runNumber >= 262424 || runNumber <= 256418 ) {


	//LHC20l7c (-> LHC16)
	// He3
	fBetheParamsHe[0] = 2.74996;
	fBetheParamsHe[1] = 13.98;
	fBetheParamsHe[2] = 0.0251843;
	fBetheParamsHe[3] = 2.04678;
	fBetheParamsHe[4] = 1.37379;
	fBetheParamsHe[5] = 0.06;
	// Triton
	fBetheParamsT[0] = 1.80227;
	fBetheParamsT[1] = 16.8019;
	fBetheParamsT[2] = 2.22419;
	fBetheParamsT[3] = 2.30938;
	fBetheParamsT[4] = 3.52324;
	fBetheParamsT[5] = 0.06;

	/*//LHC18a2b (->LHC16)
	// He3
	fBetheParamsHe[0] = 3.05245;
	fBetheParamsHe[1] = 15.7252;
	fBetheParamsHe[2] = -0.00453331;
	fBetheParamsHe[3] = 2.17241;
	fBetheParamsHe[4] = 2.88422;
	fBetheParamsHe[5] = 0.0834274;
	// Triton
	fBetheParamsT[0] = 2.74259;
	fBetheParamsT[1] = 18.3295;
	fBetheParamsT[2] = 5.91594;
	fBetheParamsT[3] = 1.93471;
	fBetheParamsT[4] = 0.292147;
	fBetheParamsT[5] = 0.0728241;*/
      }
      if (runNumber >= 256941 && runNumber <= 258537 ) {

	//LHC20l7c (-> LHC16)
	// He3
	fBetheParamsHe[0] = 2.74996;
	fBetheParamsHe[1] = 13.98;
	fBetheParamsHe[2] = 0.0251843;
	fBetheParamsHe[3] = 2.04678;
	fBetheParamsHe[4] = 1.37379;
	fBetheParamsHe[5] = 0.06;
	// Triton
	fBetheParamsT[0] = 1.80227;
	fBetheParamsT[1] = 16.8019;
	fBetheParamsT[2] = 2.22419;
	fBetheParamsT[3] = 2.30938;
	fBetheParamsT[4] = 3.52324;
	fBetheParamsT[5] = 0.06;

	/*// LHC18a2b2 (LHC16k)
	// He3
	fBetheParamsHe[0] = 2.80527;
	fBetheParamsHe[1] = 14.2379;
	fBetheParamsHe[2] = 0.0232811;
	fBetheParamsHe[3] = 2.11464;
	fBetheParamsHe[4] = 1.615;
	fBetheParamsHe[5] = 0.0815227;
	// Triton
	fBetheParamsT[0] = 1.31603;
	fBetheParamsT[1] = 36.1798;
	fBetheParamsT[2] = 493.036;
	fBetheParamsT[3] = 2.10841;
	fBetheParamsT[4] = 7.43391;
	fBetheParamsT[5] = 0.0769041;*/
      }
      if (runNumber >= 258962 && runNumber <= 259888 ) {

	//LHC20l7c (-> LHC16)
	// He3
	fBetheParamsHe[0] = 2.74996;
	fBetheParamsHe[1] = 13.98;
	fBetheParamsHe[2] = 0.0251843;
	fBetheParamsHe[3] = 2.04678;
	fBetheParamsHe[4] = 1.37379;
	fBetheParamsHe[5] = 0.06;
	// Triton
	fBetheParamsT[0] = 1.80227;
	fBetheParamsT[1] = 16.8019;
	fBetheParamsT[2] = 2.22419;
	fBetheParamsT[3] = 2.30938;
	fBetheParamsT[4] = 3.52324;
	fBetheParamsT[5] = 0.06;
	
	/*//LHC18a2b3 (->LHC16l)
	// He3
	fBetheParamsHe[0] = 2.80121;
	fBetheParamsHe[1] = 14.2397;
	fBetheParamsHe[2] = 0.0100894;
	fBetheParamsHe[3] = 2.10396;
	fBetheParamsHe[4] = 1.41608;
	fBetheParamsHe[5] = 0.0817429;
	// Triton
	fBetheParamsT[0] = 4.80597;
	fBetheParamsT[1] = 13.8813;
	fBetheParamsT[2] = 189.651;
	fBetheParamsT[3] = 2.05969;
	fBetheParamsT[4] = 4.38013;
	fBetheParamsT[5] = 0.077593;*/
      }
    }
  }
  if (runNumber >= 270581 && runNumber <= 282704) { // 2017 pp
    fYear = 2017;
    if(!fMCtrue) {
      //LHC17
      // He3
      fBetheParamsHe[0] = 3.20025;
      fBetheParamsHe[1] = 16.4971;
      fBetheParamsHe[2] = -0.0116571;
      fBetheParamsHe[3] = 2.3152;
      fBetheParamsHe[4] = 3.11135;
      fBetheParamsHe[5] = 0.06;
      // Triton
      fBetheParamsT[0] = 1.69461;
      fBetheParamsT[1] = 27.6917;
      fBetheParamsT[2] = 0.372214;
      fBetheParamsT[3] = 2.05305;
      fBetheParamsT[4] = -1.25037;
      fBetheParamsT[5] = 0.06;
    } else {

      //LHC20l7b (-> LHC17)
      // He3
      fBetheParamsHe[0] = 3.14546;
      fBetheParamsHe[1] = 16.2277;
      fBetheParamsHe[2] = -0.000523081;
      fBetheParamsHe[3] = 2.28248;
      fBetheParamsHe[4] = 2.60465;
      fBetheParamsHe[5] = 0.06;
      // Triton
      fBetheParamsT[0] = 2.88676;
      fBetheParamsT[1] = 15.3823;
      fBetheParamsT[2] = 0.580675;
      fBetheParamsT[3] = 2.28551;
      fBetheParamsT[4] = 2.47351;
      fBetheParamsT[5] = 0.06;
      
      /* // LHC18a2a (->LHC17)
      // He3
      fBetheParamsHe[0] = 3.12796;
      fBetheParamsHe[1] = 16.1359;
      fBetheParamsHe[2] = -0.00682978;
      fBetheParamsHe[3] = 2.26624;
      fBetheParamsHe[4] = 2.58652;
      fBetheParamsHe[5] = 0.0847009;
      // Triton
      fBetheParamsT[0] = 2.8303;
      fBetheParamsT[1] = 15.4337;
      fBetheParamsT[2] = 3.18352;
      fBetheParamsT[3] = 2.20975;
      fBetheParamsT[4] = 0.218244;
      fBetheParamsT[5] = 0.0780191;*/
    }
  }
  if (runNumber >= 285009 && runNumber <= 294925) { // 2018 pp
    fYear = 2018;
    if(!fMCtrue) {      
      // LHC16 + LHC18
       // He3
      fBetheParamsHe[0] = 1.81085;
      fBetheParamsHe[1] = 29.4656;
      fBetheParamsHe[2] = 0.0458225;
      fBetheParamsHe[3] = 2.08689;
      fBetheParamsHe[4] = 2.28772;
      fBetheParamsHe[5] = 0.06;
      // Triton
      fBetheParamsT[0] = 1.58385;
      fBetheParamsT[1] = 25.8334;
      fBetheParamsT[2] = 0.00908038;
      fBetheParamsT[3] = 2.24769;
      fBetheParamsT[4] = 2.87755;
      fBetheParamsT[5] = 0.06;
    } else {

      //LHC20l7a (-> LHC18)
      // He3
      fBetheParamsHe[0] = 3.07067;
      fBetheParamsHe[1] = 15.8069;
      fBetheParamsHe[2] = -0.0142383;
      fBetheParamsHe[3] = 2.15513;
      fBetheParamsHe[4] = 2.5192;
      fBetheParamsHe[5] = 0.06;
      // Triton
      fBetheParamsT[0] = 2.95171;
      fBetheParamsT[1] = 17.7223;
      fBetheParamsT[2] = 37.7979;
      fBetheParamsT[3] = 2.03313;
      fBetheParamsT[4] = 0.730268;
      fBetheParamsT[5] = 0.06;
      
      /*//LHC18a2d (->LHC18)
      // He3
      fBetheParamsHe[0] = 3.07104;
      fBetheParamsHe[1] = 15.8085;
      fBetheParamsHe[2] = 0.0150992;
      fBetheParamsHe[3] = 2.13909;
      fBetheParamsHe[4] = 2.59495;
      fBetheParamsHe[5] = 0.0865179;
      // Triton
      fBetheParamsT[0] = 2.54486;
      fBetheParamsT[1] = 17.1203;
      fBetheParamsT[2] = -0.0452007;
      fBetheParamsT[3] = 2.00988;
      fBetheParamsT[4] = 0.849292;
      fBetheParamsT[5] = 0.0768715;*/
    }
  }
}
//_____________________________________________________________________________
void AliAnalysisTaskS3ParticleYields::ResetVals(TString mode){

  if(mode == "Event"){
    fNumberV0s = -99;
    frunnumber = -99;
    MB = 0;
    HMV0 = 0;
    HMSPD = 0;
    HNU = 0;
    HQU = 0;
    fMultV0M = -99;
    fMultOfV0M = -99;
    fMultSPDTracklet = -99;
    fMultSPDCluster = -99;
    fMultRef05 = -99;
    fMultRef08 = -99;
    fSPDCluster = -99;
    fSPDTracklets = -99;
    fSPDFiredChips0 = -99;
    fSPDFiredChips1 = -99;
    fV0Multiplicity = -99;
    fMCtrue = -1;
    fNTracks = 0;
    fYear = -99;
  }
  else{
    fTrigMB = -99;			
    fTrigHMV0 = -99;
    fTrigHMSPD = -99;
    fTrigHNU = 0;
    fTrigHQU = 0; 
    fonTheFly = -99;
    fMCtrueP = -99;
    fisPrimaryP = -99;
    fisWeakP = -99;
    fisMaterialP = -99;
    fMCtrueL = -99;
    fisWeakL = -99;
    fisMaterialL = -99;
    fisPrimaryL = -99;
    fCharge = -99;
    fLambdaM = -99;
    fLambdaP = -99;
    fLambdaPx = -99;
    fLambdaPy = -99;
    fLambdaPz = -99;
    fLambdaPt = -99;
    fLambdaE = -99;
    fLambdaCt = -99;
    fLambdaY = -99;
    fLambdaDca = -99;
    fLambdaCos = -99;
    fpiP = -99;
    fpiPt = -99;
    fpiPx = -99;
    fpiPy = -99;
    fpiPz = -99;
    fpiE = -99;
    fpiy = -99;
    fpiDedx = -99;
    fpiDedxSigma = -99;
    fpichi2 = -99;
    fITSchi2Pi = -99;
    fpiNcls = -99;
    fpiNclsITS = -99;  
    fpiDca = -99;
    fpiDcaz = -99;
    fpiDcaSec = -99;
    fGeoLengthPi = -99;  
    fTOFSignalPi = -99;
    fTOFnSigmaPi = -99;
    fEtaPi = -99;
    fPhiPi = -99;
    fKinkPi = -1;
    fTPCrefitPi = -1;
    fITSrefitPi = -1;
    fpiSigmaYX = -99;
    fpiSigmaXYZ = -99;
    fpiSigmaZ = -99;
    fpP = -99;
    fpy = -99;
    fpPt = -99;
    fpPx = -99;
    fpPy = -99;
    fpPz = -99;
    fpE = -99;
    fpDedx = -99;
    fpDedxSigma = -99;
    fpchi2 = -99;
    fITSchi2P = -99;
    fpNcls = -99;
    fpNclsITS = -99;
    fpDcaSec = -99;
    fpDca = -99;
    fpDcaz = -99;
    fGeoLengthP = -99;
    fTOFSignalP = -99;
    fTOFnSigmaP = -99;
    fEtaP = -99;
    fPhiP = -99;
    fKinkP = -1;
    fTPCrefitP = -1;
    fITSrefitP = -1;
    fpSigmaYX = -99;
    fpSigmaXYZ = -99;
    fpSigmaZ = -99;
    farmalpha = -99;
    farmpt = -99; 
    fthetaP = -99;
    fthetaN = -99;
  }
  return;
}
