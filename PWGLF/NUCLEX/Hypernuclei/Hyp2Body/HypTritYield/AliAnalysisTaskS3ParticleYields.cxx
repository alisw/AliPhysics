
//--- Task for the determination of p and Lambda Yields in pp ---
//---     Author: Janik Ditzel; janik.ditzel@cern.ch     ---


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
  fPIDCheckOnly(kFALSE),
  fInputHandler(0),
  fPID(0),
  fESDevent(0),
  fStack(),
  trackCutsV0(),
  mcEvent(0),
  fV0(),
  fHistdEdx(0),
  fHistdEdxV0(0),
  fHistNumEvents(0),
  fHistTrigger(0),
  fHistV0(0),
  fHistEvents(0),
  fTree(0),
  hTree(0),
  fTreeGen(0),
  fHistogramList(NULL),
  fPrimaryVertex(),
  fMagneticField(),
  fNV0Cand(),
  fMCtrue(0),
  fEventCuts(),
  fPeriod(00),
  fTriggerMask(),
  fBetheSplines(kFALSE),
  fBetheParamsHe(),
  fBetheParamsT(),
  fMultV0M(-99),
  fMultOfV0M(-99),
  fMultSPDTracklet(-99),
  fMultSPDCluster(-99),
  fMultRef05(-99),
  fMultRef08(-99),
  tSPDCluster(-99),
  tSPDTracklets(-99),
  tSPDFiredChips0(-99),
  tSPDFiredChips1(-99),
  tV0Multiplicity(-99),
  fpLDca(-99),
  fpiNcls(-99),
  fhe3Ncls(-99),
  fpNcls(-99),
  fpLNcls(-99),
  fpiNclsITS(-99),
  fhe3NclsITS(-99),
  fpNclsITS(-99),
  fpLNclsITS(-99),
  fpiDedxSigma(-99),
  fhe3DedxSigma(-99),
  fpDedxSigma(-99),
  fpLDedxSigma(-99),
  fpiP(-99),
  fhe3P(-99),
  fpP(-99),
  fpPt(-99),
  fpchi2(-99),
  fpDcaz(-99),
  fpLP(-99),
  fpiDedx(-99),
  fhe3Dedx(-99),
  fpDedx(-99),
  fpLDedx(-99),
  farmalpha(-99),
  farmpt(-99),
  ftrig(-99),
  fz(-99),
  fmc(-99),
  fthetaP(-99),
  fthetaN(-99),
  fonTheFly(-99),
  fVertexPosition(),
  fNumberV0s(-99),
  fCentrality(-99),
  frunnumber(-99),
  fTrigger(),
  fTriggerClasses(),
  fEtaHe3(-99),
  fEtaP(-99),
  fEtaPL(-99),
  fEtaPi(-99),
  fPhiHe3(-99),
  fPhiP(-99),
  fPhiPL(-99),
  fPhiPi(-99),
  fGeoLengthHe3(-99),
  fGeoLengthP(-99),
  fGeoLengthPi(-99),
  fGeoLengthPL(-99),
  fTOFSignalHe3(-99),
  fTOFSignalP(-99),
  fTOFSignalPi(-99),
  fTOFSignalPL(-99),
  fMCtrueHe3(-99),
  fisPrimaryHe3(-99),
  fisWeakHe3(-99),
  fisMaterialHe3(-99),
  fisfromHypertriton(-99),
  fisPrimaryP(-99),
  fisWeakP(-99),
  fisMaterialP(-99),
  fMCtrueP(-99),
  fMCtrueL(-99),
  fisWeakL(-99),
  fisMaterialL(-99),
  fisPrimaryL(-99),
  fisWeakGenL(-99),
  fisMaterialGenL(-99),
  fisPrimaryGenL(-99),
  fpHe3Gen(-99),
  fyHe3Gen(-99),
  fisMaterialGenHe3(-99),
  fisPrimaryGenHe3(-99),
  fisSecondaryGenHe3(-99),
  fHe3Charge(-99),
  fpPGen(-99),
  fyPGen(-99),
  fisPrimaryGenP(-99),
  fisMaterialGenP(-99),
  fisSecondaryGenP(-99),
  fPCharge(-99),
  fpLambdaGen(-99),
  fyLambdaGen(-99),
  fmLambdaGen(-99),
  fLambdaCharge(-99),
  fmLambda(-99),
  fpLambda(-99),
  fptLambda(-99),
  fctLambda(-99),
  fdcaLambda(-99),
  fcosLambda(-99),
  fyLambda(-99),
  fhe3Pt(-99),
  fhe3chi2(-99),
  fhe3Dcaz(-99),
  fhe3Dca(-99),
  fpy(-99),
  fpiy(-99),
  fhe3y(-99),
  fpLy(-99),
  fpDcaSec(-99),
  fpLDcaSec(-99),
  fpiDcaSec(-99),
  fpiDca(-99),
  fpDca(-99)
{
    
}

// Constructor
AliAnalysisTaskS3ParticleYields::AliAnalysisTaskS3ParticleYields(const char *name)
  :AliAnalysisTaskSE(name),
   fPIDCheckOnly(kFALSE),
  fInputHandler(0),
  fPID(0),
  fESDevent(0),
  fStack(),
  trackCutsV0(),
  mcEvent(0),
  fV0(),
  fHistdEdx(0),
  fHistdEdxV0(0),
  fHistNumEvents(0),
  fHistTrigger(0),
  fHistV0(0),
  fHistEvents(0),
  fTree(0),
  hTree(0),
  fTreeGen(0),
  fHistogramList(NULL),
  fPrimaryVertex(),
  fMagneticField(),
  fNV0Cand(),
  fMCtrue(0),
  fEventCuts(),
  fPeriod(00),
  fTriggerMask(),
  fBetheSplines(kFALSE),
  fBetheParamsHe(),
  fBetheParamsT(),
  fMultV0M(-99),
  fMultOfV0M(-99),
  fMultSPDTracklet(-99),
  fMultSPDCluster(-99),
  fMultRef05(-99),
  fMultRef08(-99),
  tSPDCluster(-99),
  tSPDTracklets(-99),
  tSPDFiredChips0(-99),
  tSPDFiredChips1(-99),
  tV0Multiplicity(-99),
  fpLDca(-99),
  fpiNcls(-99),
  fhe3Ncls(-99),
  fpNcls(-99),
  fpLNcls(-99),
  fpiNclsITS(-99),
  fhe3NclsITS(-99),
  fpNclsITS(-99),
  fpLNclsITS(-99),
  fpiDedxSigma(-99),
  fhe3DedxSigma(-99),
  fpDedxSigma(-99),
  fpLDedxSigma(-99),
  fpiP(-99),
  fhe3P(-99),
  fpP(-99),
  fpPt(-99),
  fpchi2(-99),
  fpDcaz(-99),
  fpLP(-99),
  fpiDedx(-99),
  fhe3Dedx(-99),
  fpDedx(-99),
  fpLDedx(-99),
  farmalpha(-99),
  farmpt(-99),
  ftrig(-99),
  fz(-99),
  fmc(-99),
  fthetaP(-99),
  fthetaN(-99),
  fonTheFly(-99),
  fVertexPosition(),
  fNumberV0s(-99),
  fCentrality(-99),
  frunnumber(-99),
  fTrigger(),
  fTriggerClasses(),
  fEtaHe3(-99),
  fEtaP(-99),
  fEtaPL(-99),
  fEtaPi(-99),
  fPhiHe3(-99),
  fPhiP(-99),
  fPhiPL(-99),
  fPhiPi(-99),
  fGeoLengthHe3(-99),
  fGeoLengthP(-99),
  fGeoLengthPi(-99),
  fGeoLengthPL(-99),
  fTOFSignalHe3(-99),
  fTOFSignalP(-99),
  fTOFSignalPi(-99),
  fTOFSignalPL(-99),
  fMCtrueHe3(-99),
  fisPrimaryHe3(-99),
  fisWeakHe3(-99),
  fisMaterialHe3(-99),
  fisfromHypertriton(-99),
  fisPrimaryP(-99),
  fisWeakP(-99),
  fisMaterialP(-99),
  fMCtrueP(-99),
  fMCtrueL(-99),
  fisWeakL(-99),
  fisMaterialL(-99),
  fisPrimaryL(-99),
  fisWeakGenL(-99),
  fisMaterialGenL(-99),
  fisPrimaryGenL(-99),
  fpHe3Gen(-99),
  fyHe3Gen(-99),
  fisMaterialGenHe3(-99),
  fisPrimaryGenHe3(-99),
  fisSecondaryGenHe3(-99),
  fHe3Charge(-99),
  fpPGen(-99),
  fyPGen(-99),
  fisPrimaryGenP(-99),
  fisMaterialGenP(-99),
  fisSecondaryGenP(-99),
  fPCharge(-99),
  fpLambdaGen(-99),
  fyLambdaGen(-99),
  fmLambdaGen(-99),
  fLambdaCharge(-99),
  fmLambda(-99),
  fpLambda(-99),
  fptLambda(-99),
  fctLambda(-99),
  fdcaLambda(-99),
  fcosLambda(-99),
  fyLambda(-99),
  fhe3Pt(-99),
  fhe3chi2(-99),
  fhe3Dcaz(-99),
  fhe3Dca(-99),
  fpy(-99),
  fpiy(-99),
  fhe3y(-99),
  fpLy(-99),
  fpDcaSec(-99),
  fpLDcaSec(-99),
  fpiDcaSec(-99),
  fpiDca(-99),
  fpDca(-99)
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
  -1010020040,       //AntiHyperHelium 4
  1010020050,       //HyperHelium5
  -1010020050,        //AntiHyperHelium 5
  1020010040,       //DoubleHyperHydrogen 4
  -1020010040,        //AntiDoubleHyperHydrogen 4
  1010010030,       //HyperHydrogen3
  -1010010030,        //AntiHyperHydrogen3
  1010010040,       //HyperHydrogen4
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
  fHistTrigger->GetXaxis()->SetBinLabel(1,"other");
  fHistTrigger->GetXaxis()->SetBinLabel(2,"kINT7");
  fHistTrigger->GetXaxis()->SetBinLabel(3,"kHighMultV0");
  fHistTrigger->GetXaxis()->SetBinLabel(4,"kHighMultSPD");
  fHistTrigger->GetXaxis()->SetBinLabel(5,"HNU");
  fHistTrigger->GetXaxis()->SetBinLabel(6,"HQU");
  fHistTrigger->GetXaxis()->SetBinLabel(7,"HJT");
  fHistV0 = new TH1F("fHistV0","Trigger V0s",7,0,7);
  fHistV0->GetXaxis()->SetBinLabel(1,"other");
  fHistV0->GetXaxis()->SetBinLabel(2,"kINT7");
  fHistV0->GetXaxis()->SetBinLabel(3,"kHighMultV0");
  fHistV0->GetXaxis()->SetBinLabel(4,"kHighMultSPD");
  fHistV0->GetXaxis()->SetBinLabel(5,"HNU");
  fHistV0->GetXaxis()->SetBinLabel(6,"HQU");
  fHistV0->GetXaxis()->SetBinLabel(7,"HJT");
    
  fHistEvents = new TH1F("fHistV0","Trigger V0s",7,0,7);
  fHistEvents->GetXaxis()->SetBinLabel(1,"other");
  fHistEvents->GetXaxis()->SetBinLabel(2,"kINT7");
  fHistEvents->GetXaxis()->SetBinLabel(3,"kHighMultV0");
  fHistEvents->GetXaxis()->SetBinLabel(4,"kHighMultSPD");
  fHistEvents->GetXaxis()->SetBinLabel(5,"HNU");
  fHistEvents->GetXaxis()->SetBinLabel(6,"HQU");
  fHistEvents->GetXaxis()->SetBinLabel(7,"HJT");
     
  fHistogramList = new TList();
  fHistogramList->SetOwner(kTRUE);
  fHistogramList->SetName(GetName());
  fHistogramList->Add(fHistdEdx);
  fHistogramList->Add(fHistdEdxV0);
  fHistogramList->Add(fHistNumEvents);
  fHistogramList->Add(fHistTrigger);
  fEventCuts.AddQAplotsToList(fHistogramList);
  
  fTree = new TTree("treeP","fTree");
  fTree->Branch("fTrigger", &fTrigger, "fTrigger/I");
  fTree->Branch("fMultV0M",&fMultV0M,"fMultV0M/I");
  fTree->Branch("fMultOfV0M",&fMultOfV0M,"fMultOfV0M/I");
  fTree->Branch("fMultSPDTracklet",&fMultSPDTracklet,"fMultSPDTracklet/I");
  fTree->Branch("fMultSPDCluster",&fMultSPDCluster,"fMultSPDCluster/I");
  fTree->Branch("fMultRef05",&fMultRef05,"fMultRef05/I");
  fTree->Branch("fMultRef08",&fMultRef08,"fMultRef08/I");
  fTree->Branch("tSPDCluster",&tSPDCluster,"tSPDCluster/I");
  fTree->Branch("tSPDTracklets",&tSPDTracklets,"tSPDTracklets/I");
  fTree->Branch("tSPDFiredChips0",&tSPDFiredChips0,"tSPDFiredChips0/I");
  fTree->Branch("tSPDFiredChips1",&tSPDFiredChips1,"tSPDFiredChips1/I");
  fTree->Branch("tV0Multiplicity",&tV0Multiplicity,"tV0Multiplicity/I");    
  fTree->Branch("fpchi2", &fpchi2, "fpchi2/F");
  fTree->Branch("fpP", &fpP, "fpP/F");
  fTree->Branch("fpPt", &fpPt, "fpPt/F");
  fTree->Branch("fpy", &fpy, "fpy/F");
  fTree->Branch("fpDcaz", &fpDcaz ,"fpDcaz/F");
  fTree->Branch("fpDca", &fpDca ,"fpDca/F");
  fTree->Branch("fpNcls", &fpNcls, "fpNcls/F");
  fTree->Branch("fpNclsITS", &fpNclsITS, "fpNclsITS/F");
  fTree->Branch("fpDedxSigma", &fpDedxSigma, "fpDedxSigma/F");
  fTree->Branch("fpDedx", &fpDedx, "fpDedx/F");
  fTree->Branch("fEtaP", &fEtaP, "fEtaP/F");
  fTree->Branch("fPhiP", &fPhiP, "fPhiP/F");
  fTree->Branch("fGeoLengthP", &fGeoLengthP, "fGeoLengthP/F");
  fTree->Branch("fTOFSignalP", &fTOFSignalP, "fTOFSignalP/F");
  fTree->Branch("fMCtrueP", &fMCtrueP, "fMCtrueP/I");
  fTree->Branch("fisPrimaryP", &fisPrimaryP, "fisPrimaryP/I");
  fTree->Branch("fisWeakP", &fisWeakP, "fisWeakP/I");
  fTree->Branch("fisMaterialP", &fisMaterialP, "fisMaterialP/I");
  fTree->Branch("frunnumber", &frunnumber,"frunnumber/I");
  fTree->Branch("fCentrality", &fCentrality, "fCentrality/I");
  fTree->Branch("fz", &fz, "fz/I");
  fTree->Branch("fmc", &fmc, "fmc/I");
  
  hTree = new TTree("treeL","hTree");
  hTree->Branch("fTrigger", &fTrigger, "fTrigger/I");   
  hTree->Branch("fMultV0M",&fMultV0M,"fMultV0M/I");
  hTree->Branch("fMultOfV0M",&fMultOfV0M,"fMultOfV0M/I");
  hTree->Branch("fMultSPDTracklet",&fMultSPDTracklet,"fMultSPDTracklet/I");
  hTree->Branch("fMultSPDCluster",&fMultSPDCluster,"fMultSPDCluster/I");
  hTree->Branch("fMultRef05",&fMultRef05,"fMultRef05/I");
  hTree->Branch("fMultRef08",&fMultRef08,"fMultRef08/I");
  hTree->Branch("tSPDCluster",&tSPDCluster,"tSPDCluster/I");
  hTree->Branch("tSPDTracklets",&tSPDTracklets,"tSPDTracklets/I");
  hTree->Branch("tSPDFiredChips0",&tSPDFiredChips0,"tSPDFiredChips0/I");
  hTree->Branch("tSPDFiredChips1",&tSPDFiredChips1,"tSPDFiredChips1/I");
  hTree->Branch("tV0Multiplicity",&tV0Multiplicity,"tV0Multiplicity/I");
  hTree->Branch("fmLambda", &fmLambda, "fmLambda/F");
  hTree->Branch("fpLambda", &fpLambda, "fpLambda/F");
  hTree->Branch("fptLambda", &fptLambda, "fptLambda/F");
  hTree->Branch("fctLambda", &fctLambda, "fctLambda/F");
  hTree->Branch("fpLP", &fpLP, "fpLP/F");
  hTree->Branch("fpiP", &fpiP, "fpiP/F");
  hTree->Branch("fdcaLambda", &fdcaLambda, "fdcaLambda/F");
  hTree->Branch("fcosLambda", &fcosLambda, "fcosLambda/F");
  hTree->Branch("fyLambda", &fyLambda, "fyLambda/F");
  hTree->Branch("fpLy", &fpLy, "fpLy/F");
  hTree->Branch("fpiy", &fpiy, "fpiy/F");
  hTree->Branch("fpiDca", &fpiDca, "fpiDca/F");
  hTree->Branch("fpLDca", &fpLDca ,"fpLDca/F");
  hTree->Branch("fpiDcaSec", &fpiDcaSec, "fpiDcaSec/F");
  hTree->Branch("fpLDcaSec", &fpLDcaSec ,"fpLDcaSec/F");
  hTree->Branch("fpiNcls", &fpiNcls, "fpiNcls/F");
  hTree->Branch("fpLNcls", &fpLNcls, "fpLNcls/F");
  hTree->Branch("fpiNclsITS", &fpiNclsITS, "fpiNclsITS/F");
  hTree->Branch("fpLNclsITS", &fpLNclsITS, "fpLNclsITS/F");
  hTree->Branch("fpLDedxSigma", &fpLDedxSigma, "fpLDedxSigma/F");
  hTree->Branch("fpiDedxSigma", &fpiDedxSigma, "fpiDedxSigma/F");
  hTree->Branch("fpiDedx", &fpiDedx, "fpiDedx/F");
  hTree->Branch("fpLDedx", &fpLDedx, "fpLDedx/F");
  hTree->Branch("farmalpha", &farmalpha, "farmalpha/F");
  hTree->Branch("farmpt", &farmpt, "farmpt/F");
  hTree->Branch("fEtaPL", &fEtaPL, "fEtaPL/F");
  hTree->Branch("fEtaPi", &fEtaPi, "fEtaPi/F");
  hTree->Branch("fPhiPL", &fPhiPL, "fPhiPL/F");
  hTree->Branch("fPhiPi", &fPhiPi, "fPhiPi/F");
  hTree->Branch("fGeoLengthPL", &fGeoLengthPL, "fGeoLengthPL/F");
  hTree->Branch("fGeoLengthPi", &fGeoLengthPi, "fGeoLengthPi/F");
  hTree->Branch("fTOFSignalPL", &fTOFSignalPL, "fTOFSignalPL/F");
  hTree->Branch("fTOFSignalPi", &fTOFSignalPi, "fTOFSignalPi/F");
  hTree->Branch("fMCtrueL", &fMCtrueL, "fMCtrueL/I");
  hTree->Branch("fonTheFly", &fonTheFly, "fonTheFly/I");
  hTree->Branch("frunnumber", &frunnumber,"frunnumber/I");
  hTree->Branch("fNumberV0s", &fNumberV0s, "fNumberV0s/I");
  hTree->Branch("fCentrality", &fCentrality, "fCentrality/I");
  hTree->Branch("fisPrimaryL", &fisPrimaryL, "fisPrimaryL/I");
  hTree->Branch("fisWeakL", &fisWeakL, "fisWeakL/I");
  hTree->Branch("fisMaterialL", &fisMaterialL, "fisMaterialL/I");
  hTree->Branch("fz", &fz, "fz/I");
  hTree->Branch("fmc", &fmc, "fmc/I");
    
  fTreeGen = new TTree("treeGenLp","fTreeGen");
  fTreeGen->Branch("fTrigger", &fTrigger, "fTrigger/I");
  fTreeGen->Branch("fMultV0M",&fMultV0M,"fMultV0M/I");
  fTreeGen->Branch("fMultOfV0M",&fMultOfV0M,"fMultOfV0M/I");
  fTreeGen->Branch("fMultSPDTracklet",&fMultSPDTracklet,"fMultSPDTracklet/I");
  fTreeGen->Branch("fMultSPDCluster",&fMultSPDCluster,"fMultSPDCluster/I");
  fTreeGen->Branch("fMultRef05",&fMultRef05,"fMultRef05/I");
  fTreeGen->Branch("fMultRef08",&fMultRef08,"fMultRef08/I");
  fTreeGen->Branch("tSPDCluster",&tSPDCluster,"tSPDCluster/I");
  fTreeGen->Branch("tSPDTracklets",&tSPDTracklets,"tSPDTracklets/I");
  fTreeGen->Branch("tSPDFiredChips0",&tSPDFiredChips0,"tSPDFiredChips0/I");
  fTreeGen->Branch("tSPDFiredChips1",&tSPDFiredChips1,"tSPDFiredChips1/I");
  fTreeGen->Branch("tV0Multiplicity",&tV0Multiplicity,"tV0Multiplicity/I");
  fTreeGen->Branch("fpPGen", &fpPGen, "fpPGen/F");
  fTreeGen->Branch("fyPGen", &fyPGen, "fyPGen/F");
  fTreeGen->Branch("fisPrimaryGenP", &fisPrimaryGenP, "fisPrimaryGenP/I");
  fTreeGen->Branch("fisSecondaryGenP", &fisSecondaryGenP, "fisSecondaryGenP/I");
  fTreeGen->Branch("fisMaterialGenP", &fisMaterialGenP, "fisMaterialGenP/I");
  fTreeGen->Branch("fPCharge", &fPCharge, "fPCharge/I");
  fTreeGen->Branch("fpLambdaGen", &fpLambdaGen, "fpLambdaGen/F");
  fTreeGen->Branch("fyLambdaGen", &fyLambdaGen, "fyLambdaGen/F");
  fTreeGen->Branch("fmLambdaGen", &fmLambdaGen, "fmLambdaGen/F");
  fTreeGen->Branch("fLambdaCharge", &fLambdaCharge, "fLambdaCharge/I");
  fTreeGen->Branch("fisPrimaryGenL", &fisPrimaryGenL, "fisPrimaryGenL/I");
  fTreeGen->Branch("fisWeakGenL", &fisWeakGenL, "fisWeakGenL/I");
  fTreeGen->Branch("fisMaterialGenL", &fisMaterialGenL, "fisMaterialGenL/I");
    
  PostData(1, fHistogramList);
  PostData(2, fTree);
  PostData(3, hTree);
  PostData(4, fTreeGen);
    
  //********--------******** Info ********--------********
  //  EParticleType :                                   //
  //  kElectron = 0, kMuon = 1, kPion = 2, kKaon = 3,   //
  //  kProton = 4, kDeuteron = 5, kTriton = 6, kHe3 = 7 //
  //  kAlpha = 8, kPhoton = 9, kPi0 = 10, kNeutron = 11 //
  //  kKaon0 = 12, kEleCon = 13, kUnknown = 14          //
  //                                                    //
  //  Masses: 4Li: 3.74958 GeV/c^2                      //
  //          4He: 3.727379 GeV/c^2                     //
  //          3He: 2.80923 GeV/c^2                      //
  //            t: 2.80925 GeV/c^2                      //
  //            d: 1.875613 GeV/c^2                     //
  //            p: 0.93827 GeV/c^2                      //
  //           pi: 0.13957 GeV/c^2                      //
  //          3LH: 2.99131 GeV/c^2                      //
  //          4LH: 3.931   GeV/c^2                      //
  //         4LHe: 3.929   GeV/c^2                      //
  //         5LHe: 4.841   GeV/c^2                      //
  //         4LLH: 4.106   GeV/c^2                      //
  //*******-------********--------********--------********
    
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
  if (fPeriod == 2016 || fPeriod == 2017 || fPeriod == 2018) {
    if(!fEventCuts.AcceptEvent(fESDevent)) {
      PostData(1,fHistogramList);
      return;
    }
  }
  SetMultiplicity();
  frunnumber = fESDevent->GetRunNumber();
  SetBetheBlochParams(frunnumber);
  TriggerSelection();

  fHistNumEvents->Fill(1);
  fMagneticField  = fESDevent->GetMagneticField();
  const AliESDVertex *vertex = fESDevent->GetPrimaryVertexSPD();
  fPrimaryVertex.SetXYZ(vertex->GetX(),vertex->GetY(),vertex->GetZ());
  fVertexPosition = fPrimaryVertex;

  trackCutsV0 = new AliESDtrackCuts("AlitrackCutsV0", "AlitrackCutsV0");   
  trackCutsV0->SetEtaRange(-0.8,0.8);
  trackCutsV0->SetAcceptKinkDaughters(kFALSE);
  trackCutsV0->SetRequireTPCRefit(kTRUE);
  trackCutsV0->SetMaxChi2PerClusterTPC(5);
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
    Float_t xv[2];
    Float_t yv[3];
    trackA->GetImpactParameters(xv,yv);

    if (TMath::Abs(Bethe(*trackA, AliPID::ParticleMass(AliPID::kProton), 1, fBetheParamsT)) < 3){
      if(TOFSignal(*trackA) > 1.1 || TOFSignal(*trackA) < 0.8) continue;
      TLorentzVector fd(0.,0.,0.,0.);
      fd.SetXYZM(trackA->Px(), trackA->Py(), trackA->Pz(), AliPID::ParticleMass(AliPID::kProton));
      fpPt = fd.Pt();
      fpP =trackA->GetInnerParam()->GetP();
      fpy = fd.Rapidity();
      fpNcls = trackA->GetTPCNcls();
      fpNclsITS = trackA->GetNumberOfITSClusters();
      fpDedx = trackA->GetTPCsignal();
      fGeoLengthP = GeoLength(*trackA);
      fPhiP = trackA->Phi();
      fEtaP = trackA->Eta();
      fTOFSignalP = TOFSignal(*trackA);
      fpchi2 = trackA->GetTPCchi2()/fpNcls;
      fpDca = xv[0];
      fpDcaz = xv[1];
      if(signA < 0) fz = -1;
      if(signA > 0) fz = 1;
      if (fBetheSplines) {
	fpDedxSigma = fPID->NumberOfSigmasTPC(trackA, AliPID::kProton);
      } else {
	fpDedxSigma = Bethe(*trackA, AliPID::ParticleMass(AliPID::kProton), 1, fBetheParamsT);
      }
      if (fMCtrue) {
	Int_t label = trackA->GetLabel();			
	TParticle *particle = fStack->Particle(TMath::Abs(label));
	Int_t labelMother = particle->GetFirstMother();			
	TParticle *particleMother = fStack->Particle(TMath::Abs(labelMother));
                
	fMCtrueP = TMath::Abs(particle->GetPdgCode()) == fgkPdgCode[kPDGProton];
	fisPrimaryP = fStack->IsPhysicalPrimary(TMath::Abs(label));
	fisWeakP = fStack->IsSecondaryFromWeakDecay(TMath::Abs(label));
	fisMaterialP = fStack->IsSecondaryFromMaterial(TMath::Abs(label));
      }
      fTree->Fill();
      ResetVals("");
    }
  }
}
void AliAnalysisTaskS3ParticleYields::LambdaV0s(){
  //N V0's
  fNV0Cand = 0;
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
    if (TMath::Abs(fPID->NumberOfSigmasTPC(trackP, AliPID::kPion)) < 3) {
      pionPositive = kTRUE;
    }
    if (TMath::Abs(fPID->NumberOfSigmasTPC(trackN, AliPID::kPion)) < 3) {
      pionNegative = kTRUE;
    }
    //proton
    if (TMath::Abs(Bethe(*trackP, AliPID::ParticleMass(AliPID::kProton),  1, fBetheParamsT)) < 3) {
      protonPositive = kTRUE;
    }
    if (TMath::Abs(Bethe(*trackN, AliPID::ParticleMass(AliPID::kProton), 1, fBetheParamsT)) < 3) {
      protonNegative = kTRUE;
    }    
    //
    if(!(pionNegative && protonPositive) && !(pionPositive && protonNegative))continue;
    //Lambda
    if (protonPositive && pionNegative) {
      //Vectors
      TLorentzVector proton(0.,0.,0.,0.);
      proton.SetXYZM(trackP->Px(), trackP->Py(), trackP->Pz(), AliPID::ParticleMass(AliPID::kProton));
      TLorentzVector pion(0.,0.,0.,0.);
      pion.SetXYZM(trackN->Px(), trackN->Py(), trackN->Pz(), AliPID::ParticleMass(AliPID::kPion));
      TLorentzVector Lambda(0.,0.,0.,0.);
      Lambda = pion + proton;
      TVector3 secVertex(fV0->Xv(), fV0->Yv(), fV0->Zv());      
      //particle = 1 anti particle = -1 (charge = 0)
      fz = 1;
      //V0 information
      fdcaLambda = fV0->GetDcaV0Daughters();
      fcosLambda = fV0->GetV0CosineOfPointingAngle();
      if(fV0->GetOnFlyStatus()) fonTheFly = 1;
      //daughter information
      fpLP = trackP->GetInnerParam()->GetP();
      fpiP = trackN->GetInnerParam()->GetP();
      fpLy = proton.Rapidity();
      fpiy = pion.Rapidity();
      fpLDedx = trackP->GetTPCsignal();
      fpiDedx = trackN->GetTPCsignal();
      fEtaPL = trackP->Eta();
      fEtaPi = trackN->Eta();
      fPhiPL = trackP->Phi();
      fPhiPi = trackN->Phi();
      fGeoLengthPL = GeoLength(*trackP);
      fGeoLengthPi = GeoLength(*trackN);
      fTOFSignalPL = TOFSignal(*trackP);
      fTOFSignalPi = TOFSignal(*trackN);     
      fpLDedxSigma = Bethe(*trackP, AliPID::ParticleMass(AliPID::kProton), 1, fBetheParamsT);
      fpiDedxSigma = fPID->NumberOfSigmasTPC(trackN, AliPID::kPion);
      fpLDca = TMath::Abs(trackP->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
      fpiDca = TMath::Abs(trackN->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
      fpLDcaSec = TMath::Abs(trackP->GetD(secVertex.X(), secVertex.Y(), fMagneticField));
      fpiDcaSec = TMath::Abs(trackN->GetD(secVertex.X(), secVertex.Y(), fMagneticField));
      fpLNcls = trackP->GetTPCNcls();
      fpiNcls = trackN->GetTPCNcls();           
      fpLNclsITS = trackP->GetNumberOfITSClusters();
      fpiNclsITS = trackN->GetNumberOfITSClusters();                 
      //mother information    
      fmLambda = Lambda.M();
      fpLambda = Lambda.P();
      fptLambda = Lambda.Pt();
      fyLambda = Lambda.Rapidity();
      //Lifetime
      secVertex = secVertex - fPrimaryVertex;
      fctLambda = secVertex.Mag() * fmLambda / fpLambda;
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
	TParticle *particle = fStack->Particle(TMath::Abs(label));
	Int_t labelMother = particle->GetFirstMother();			
	TParticle *particleMother = fStack->Particle(TMath::Abs(labelMother));
                
	Int_t label1 = trackN->GetLabel();
	TParticle *particle1 = fStack->Particle(TMath::Abs(label1));
	Int_t labelMother1 = particle1->GetFirstMother();			
	TParticle *particleMother1 = fStack->Particle(TMath::Abs(labelMother1));
                
	if(particleMother->GetPdgCode() == fgkPdgCode[kPDGLambda] && particleMother1->GetPdgCode() == fgkPdgCode[kPDGLambda] && labelMother == labelMother1){
	  if(TMath::Abs(label) == TMath::Abs(GetLabel(labelMother, fgkPdgCode[kPDGProton]))		   
	       && TMath::Abs(label1) == TMath::Abs(GetLabel(labelMother, fgkPdgCode[kPDGPionMinus]))) {
	    fMCtrueL = kTRUE;
	    fisPrimaryL = fStack->IsPhysicalPrimary(TMath::Abs(labelMother1));
	    fisWeakL = fStack->IsSecondaryFromWeakDecay(TMath::Abs(labelMother1));
	    fisMaterialL = fStack->IsSecondaryFromMaterial(TMath::Abs(labelMother1));
	  }
	}
      }
      fNV0Cand = fNV0Cand + 1;
      fNumberV0s = fNV0Cand;
      hTree->Fill();
      ResetVals("");
    }//end        
    else if (protonNegative && pionPositive) {
       //Vectors
      TLorentzVector proton(0.,0.,0.,0.);
      proton.SetXYZM(trackN->Px(), trackN->Py(), trackN->Pz(), AliPID::ParticleMass(AliPID::kProton));
      TLorentzVector pion(0.,0.,0.,0.);
      pion.SetXYZM(trackP->Px(), trackP->Py(), trackP->Pz(), AliPID::ParticleMass(AliPID::kPion));
      TLorentzVector Lambda(0.,0.,0.,0.);
      Lambda = pion + proton;
      TVector3 secVertex(fV0->Xv(), fV0->Yv(), fV0->Zv());      
      //particle = 1 anti-particle = -1 (charge = 0)
      fz = -1;
      //V0 information
      fdcaLambda = fV0->GetDcaV0Daughters();
      fcosLambda = fV0->GetV0CosineOfPointingAngle();
      if(fV0->GetOnFlyStatus()) fonTheFly = 1;
      //daughter information
      fpLP = trackN->GetInnerParam()->GetP();
      fpiP = trackP->GetInnerParam()->GetP();
      fpLy = proton.Rapidity();
      fpiy = pion.Rapidity();
      fpLDedx = trackN->GetTPCsignal();
      fpiDedx = trackP->GetTPCsignal();
      fEtaPL = trackN->Eta();
      fEtaPi = trackP->Eta();
      fPhiPL = trackN->Phi();
      fPhiPi = trackP->Phi();
      fGeoLengthPL = GeoLength(*trackN);
      fGeoLengthPi = GeoLength(*trackP);
      fTOFSignalPL = TOFSignal(*trackN);
      fTOFSignalPi = TOFSignal(*trackP);     
      fpLDedxSigma = Bethe(*trackN, AliPID::ParticleMass(AliPID::kProton), 1, fBetheParamsT);
      fpiDedxSigma = fPID->NumberOfSigmasTPC(trackP, AliPID::kPion);
      fpLDca = TMath::Abs(trackN->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
      fpiDca = TMath::Abs(trackP->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
      fpLDcaSec = TMath::Abs(trackN->GetD(secVertex.X(), secVertex.Y(), fMagneticField));
      fpiDcaSec = TMath::Abs(trackP->GetD(secVertex.X(), secVertex.Y(), fMagneticField));
      fpLNcls = trackN->GetTPCNcls();
      fpiNcls = trackP->GetTPCNcls();           
      fpLNclsITS = trackN->GetNumberOfITSClusters();
      fpiNclsITS = trackP->GetNumberOfITSClusters();                 
      //mother information    
      fmLambda = Lambda.M();
      fpLambda = Lambda.P();
      fptLambda = Lambda.Pt();
      fyLambda = Lambda.Rapidity();
      //Lifetime
      secVertex = secVertex - fPrimaryVertex;
      fctLambda = secVertex.Mag() * fmLambda / fpLambda;
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
	TParticle *particle = fStack->Particle(TMath::Abs(label));
	Int_t labelMother = particle->GetFirstMother();			
	TParticle *particleMother = fStack->Particle(TMath::Abs(labelMother));
                
	Int_t label1 = trackP->GetLabel();
	TParticle *particle1 = fStack->Particle(TMath::Abs(label1));
	Int_t labelMother1 = particle1->GetFirstMother();			
	TParticle *particleMother1 = fStack->Particle(TMath::Abs(labelMother1));
                
	if(particleMother->GetPdgCode() == fgkPdgCode[kPDGAntiLambda] && particleMother1->GetPdgCode() == fgkPdgCode[kPDGAntiLambda] && labelMother == labelMother1){
	  if(TMath::Abs(label) == TMath::Abs(GetLabel(labelMother, fgkPdgCode[kPDGAntiProton]))		   
	       && TMath::Abs(label1) == TMath::Abs(GetLabel(labelMother, fgkPdgCode[kPDGPionPlus]))) {
	    fMCtrueL = kTRUE;
	    fisPrimaryL = fStack->IsPhysicalPrimary(TMath::Abs(labelMother1));
	    fisWeakL = fStack->IsSecondaryFromWeakDecay(TMath::Abs(labelMother1));
	    fisMaterialL = fStack->IsSecondaryFromMaterial(TMath::Abs(labelMother1));
	  }
	}
      }
      fNV0Cand = fNV0Cand + 1;
      fNumberV0s = (fNV0Cand);
      hTree->Fill();
      ResetVals("");
    }//end
  }
}
void AliAnalysisTaskS3ParticleYields::MCGenerated() {   
  for(Int_t stackN = 0; stackN < fStack->GetNtrack(); stackN++){
    Int_t count = 0;
    TParticle *particleMother = fStack->Particle(TMath::Abs(stackN));
    TLorentzVector part(0.,0.,0.,0.);
    TLorentzVector partd1(0.,0.,0.,0.);
    TLorentzVector partd2(0.,0.,0.,0.);
    TParticle *tparticleFirstDaughter;
    TParticle *tparticleSecondDaughter;
    //Proton    
    if(TMath::Abs(particleMother->GetPdgCode()) == fgkPdgCode[kPDGProton]){
      part.SetXYZM(particleMother->Px(), particleMother->Py(), particleMother->Pz(), AliPID::ParticleMass(AliPID::kProton));
      fpPGen = part.Pt();
      fyPGen = part.Rapidity();
      fisPrimaryGenP = fStack->IsPhysicalPrimary(TMath::Abs(stackN));
      fisSecondaryGenP = fStack->IsSecondaryFromWeakDecay(TMath::Abs(stackN));
      fisMaterialGenP = fStack->IsSecondaryFromMaterial(TMath::Abs(stackN));
      if(particleMother->GetPdgCode() == fgkPdgCode[kPDGProton]) fPCharge = 1;
      if(particleMother->GetPdgCode() == fgkPdgCode[kPDGAntiProton]) fPCharge = -1;
      count++;
    }
    //Lambda
    if(particleMother->GetPdgCode() == fgkPdgCode[kPDGLambda]) {
      tparticleFirstDaughter = fStack->Particle(TMath::Abs(GetLabel(stackN, fgkPdgCode[kPDGProton])));
      tparticleSecondDaughter = fStack->Particle(TMath::Abs(GetLabel(stackN, fgkPdgCode[kPDGPionMinus])));
      fLambdaCharge = 1;
      if(tparticleFirstDaughter->GetPdgCode() == fgkPdgCode[kPDGProton] && tparticleSecondDaughter->GetPdgCode() == fgkPdgCode[kPDGPionMinus]){
	partd1.SetXYZM(tparticleFirstDaughter->Px(), tparticleFirstDaughter->Py(), tparticleFirstDaughter->Pz(), AliPID::ParticleMass(AliPID::kProton));
	partd2.SetXYZM(tparticleSecondDaughter->Px(), tparticleSecondDaughter->Py(), tparticleSecondDaughter->Pz(), AliPID::ParticleMass(AliPID::kPion));
	part = partd1 + partd2;
	fpLambdaGen = part.Pt();
	fyLambdaGen = part.Rapidity();
	fmLambdaGen = part.M();
	fisPrimaryGenL = fStack->IsPhysicalPrimary(TMath::Abs(stackN));
	fisWeakGenL = fStack->IsSecondaryFromWeakDecay(TMath::Abs(stackN));
	fisMaterialGenL = fStack->IsSecondaryFromMaterial(TMath::Abs(stackN));
	count++;
      }
    }
    //AntiLambda
    if(particleMother->GetPdgCode() == fgkPdgCode[kPDGAntiLambda]) {
      tparticleFirstDaughter = fStack->Particle(TMath::Abs(GetLabel(stackN, fgkPdgCode[kPDGAntiProton])));
      tparticleSecondDaughter = fStack->Particle(TMath::Abs(GetLabel(stackN, fgkPdgCode[kPDGPionPlus])));
      fLambdaCharge = -1;
      if(tparticleFirstDaughter->GetPdgCode() == fgkPdgCode[kPDGAntiProton] && tparticleSecondDaughter->GetPdgCode() == fgkPdgCode[kPDGPionPlus]){
	partd1.SetXYZM(tparticleFirstDaughter->Px(), tparticleFirstDaughter->Py(), tparticleFirstDaughter->Pz(), AliPID::ParticleMass(AliPID::kProton));
	partd2.SetXYZM(tparticleSecondDaughter->Px(), tparticleSecondDaughter->Py(), tparticleSecondDaughter->Pz(), AliPID::ParticleMass(AliPID::kPion));
	part = partd1 + partd2;
	fpLambdaGen = part.Pt();
	fyLambdaGen = part.Rapidity();
	fmLambdaGen = part.M();
	fisPrimaryGenL = fStack->IsPhysicalPrimary(TMath::Abs(stackN));
	fisWeakGenL = fStack->IsSecondaryFromWeakDecay(TMath::Abs(stackN));
	fisMaterialGenL = fStack->IsSecondaryFromMaterial(TMath::Abs(stackN));
	count++;
      }
    }
    if(count>0) {
      fTreeGen->Fill();
      ResetVals("");
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
  length =track.GetIntegratedLength();
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
  
  if (!fMCtrue){
    fTriggerClasses = fESDevent->GetFiredTriggerClasses();
    if ((fInputHandler->IsEventSelected() & AliVEvent::kINT7)) fTrigger = 1;
    if ((fInputHandler->IsEventSelected() & AliVEvent::kHighMultV0)) fTrigger = 2;
    if ((fInputHandler->IsEventSelected() & AliVEvent::kHighMultSPD)) fTrigger = 3;    
    if (fTriggerClasses.Contains("HNU")) fTrigger = 4;
    if (fTriggerClasses.Contains("HQU")) fTrigger = 5;
    if (fTriggerClasses.Contains("HJT")) fTrigger = 6;  
    fHistTrigger->Fill(fTrigger);
  } else {
    fTriggerClasses = fESDevent->GetFiredTriggerClasses();
    if ((fInputHandler->IsEventSelected() & AliVEvent::kINT7)) fTrigger = 1;
    if ((fInputHandler->IsEventSelected() & AliVEvent::kHighMultV0)) fTrigger = 2;
    if ((fInputHandler->IsEventSelected() & AliVEvent::kHighMultSPD)) fTrigger = 3;    
    if (fTriggerClasses.Contains("HNU")) fTrigger = 4;
    if (fTriggerClasses.Contains("HQU")) fTrigger = 5;
    if (fTriggerClasses.Contains("HJT")) fTrigger = 6;  
    // MC: simulate TRD trigger
    Int_t nTrdTracks = fESDevent->GetNumberOfTrdTracks();
    if (nTrdTracks > 0) {
      for (Int_t iTrack = 0; iTrack < nTrdTracks; ++iTrack) {
	AliESDTrdTrack* trdTrack = fESDevent->GetTrdTrack(iTrack);
	if (!trdTrack) continue;
	// simulate HNU
	if((trdTrack->GetPID() >= 255 && trdTrack->GetNTracklets() == 4) ||
	   (trdTrack->GetPID() >= 235 && trdTrack->GetNTracklets() > 4)) {
	  fTrigger = 4;
	}
	// simulate HQU
	if (TMath::Abs(trdTrack->GetPt()) >= 256 &&
	    trdTrack->GetPID() >= 130 &&
	    trdTrack->GetNTracklets() >= 5 &&
	    (trdTrack->GetLayerMask() & 1) ){
	  Float_t sag = GetInvPtDevFromBC(trdTrack->GetB(), trdTrack->GetC());
	  if (sag < 0.2 && sag > -0.2) {
	    fTrigger = 5;
	  }
	}
      }
    }
    fHistTrigger->Fill(fTrigger);
  }
  // additional information for high multiplicity trigger
  AliESDVZERO *vzero = fESDevent->GetVZEROData();
  tV0Multiplicity = 0;
  for (Int_t ii = 0; ii < 64; ii++){
    tV0Multiplicity += vzero->GetMultiplicity(ii);
  }
  AliMultiplicity *multSPD = fESDevent->GetMultiplicity();
  tSPDCluster	= multSPD->GetNumberOfSPDClusters();
  tSPDTracklets = multSPD->GetNumberOfTracklets();
  tSPDFiredChips0 = multSPD->GetNumberOfFiredChips(0);
  tSPDFiredChips1 = multSPD->GetNumberOfFiredChips(1);
    
  Bool_t isTriggered = kTRUE;
  if (fTrigger == 0) isTriggered = kFALSE;
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
void AliAnalysisTaskS3ParticleYields::SetBetheBlochParams(Int_t runNumber) {
  // set Bethe-Bloch parameter
  if (runNumber >= 252235 && runNumber <= 264347 ) { // 2016 pp
    if(!fMCtrue) { // Data
      // LHC16 + LHC18
      // He3
      fBetheParamsT[0] = 0.427978;
      fBetheParamsT[1] = 105.46;
      fBetheParamsT[2] =-7.08642e-07;
      fBetheParamsT[3] = 2.23332;
      fBetheParamsT[4] = 18.8231;
      fBetheParamsT[5] = 0.06;
      // Triton
      fBetheParamsHe[0] = 1.81085;
      fBetheParamsHe[1] = 29.4656;
      fBetheParamsHe[2] = 0.0458225;
      fBetheParamsHe[3] = 2.08689;
      fBetheParamsHe[4] = 2.28772;
      fBetheParamsHe[5] = 0.06;
    } else { // MC
      if (runNumber >= 262424 || runNumber <= 256418 ) {
	//LHC18a2b (->LHC16)
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
	fBetheParamsT[5] = 0.0728241;
      }
      if (runNumber >= 256941 && runNumber <= 258537 ) {
	// LHC18a2b2 (LHC16k)
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
	fBetheParamsT[5] = 0.0769041;
      }
      if (runNumber >= 258962 && runNumber <= 259888 ) {
	//LHC18a2b3 (->LHC16l)
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
	fBetheParamsT[5] = 0.077593;
      }
    }
  }
  if (runNumber >= 270581 && runNumber <= 282704) { // 2017 pp
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
      fBetheParamsT[0] = 0.420434;
      fBetheParamsT[1] = 106.102;
      fBetheParamsT[2] = -3.15587e-07;
      fBetheParamsT[3] = 2.32499;
      fBetheParamsT[4] = 21.3439;
      fBetheParamsT[5] = 0.06;
    } else {
      // LHC18a2a (->LHC17)
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
      fBetheParamsT[5] = 0.0780191;
    }
  }
  if (runNumber >= 285009 && runNumber <= 294925) { // 2018 pp
    if(!fMCtrue) {
      // LHC16 + LHC18
      // He3
      fBetheParamsT[0] = 0.427978;
      fBetheParamsT[1] = 105.46;
      fBetheParamsT[2] =-7.08642e-07;
      fBetheParamsT[3] = 2.23332;
      fBetheParamsT[4] = 18.8231;
      fBetheParamsT[5] = 0.06;
      // Triton
      fBetheParamsHe[0] = 1.81085;
      fBetheParamsHe[1] = 29.4656;
      fBetheParamsHe[2] = 0.0458225;
      fBetheParamsHe[3] = 2.08689;
      fBetheParamsHe[4] = 2.28772;
      fBetheParamsHe[5] = 0.06;
    } else {
      //LHC18a2d (->LHC18)
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
      fBetheParamsT[5] = 0.0768715;
    }
  }
}
//_____________________________________________________________________________
void AliAnalysisTaskS3ParticleYields::ResetVals(TString mode){

  if(mode == "Event"){
    fMultV0M = -99;
    fMultOfV0M = -99;
    fMultSPDTracklet = -99;
    fMultSPDCluster = -99;
    fMultRef05 = -99;
    fMultRef08 = -99;
    tSPDCluster = -99;
    tSPDTracklets = -99;
    tSPDFiredChips0 = -99;
    tSPDFiredChips1 = -99;
    tV0Multiplicity = -99;
    fmc = -99;
  }
  else{
    fpLDca = -99;
    fpiNcls = -99;
    fhe3Ncls = -99;
    fpNcls = -99;
    fpLNcls = -99;
    fpiNclsITS = -99;
    fhe3NclsITS = -99;
    fpNclsITS = -99;
    fpLNclsITS = -99;
    fpiDedxSigma = -99;
    fhe3DedxSigma = -99;
    fpDedxSigma = -99;
    fpLDedxSigma = -99;
    fpiP = -99;
    fhe3P = -99;
    fpP = -99;
    fpPt = -99;
    fpchi2 = -99;
    fpDcaz = -99;
    fpLP = -99;
    fpiDedx = -99;
    fhe3Dedx = -99;
    fpDedx = -99;
    fpLDedx = -99;
    farmalpha = -99;
    farmpt = -99;
    ftrig = -99;
    fz = -99;  
    fthetaP = -99;
    fthetaN = -99;
    fonTheFly = -99;
    fEtaHe3 = -99;
    fEtaP = -99;
    fEtaPL = -99;
    fEtaPi = -99;
    fPhiHe3 = -99;
    fPhiP = -99;
    fPhiPL = -99;
    fPhiPi = -99;
    fGeoLengthHe3 = -99;
    fGeoLengthP = -99;
    fGeoLengthPi = -99;
    fGeoLengthPL = -99;
    fTOFSignalHe3 = -99;
    fTOFSignalP = -99;
    fTOFSignalPi = -99;
    fTOFSignalPL = -99;
    fMCtrueHe3 = -99;
    fisPrimaryHe3 = -99;
    fisWeakHe3 = -99;
    fisMaterialHe3 = -99;
    fisfromHypertriton = -99;
    fisPrimaryP = -99;
    fisWeakP = -99;
    fisMaterialP = -99;
    fMCtrueP = -99;
    fMCtrueL = -99;
    fisWeakL = -99;
    fisMaterialL = -99;
    fisPrimaryL = -99;
    fisWeakGenL = -99;
    fisMaterialGenL = -99;
    fisPrimaryGenL = -99;
    fpHe3Gen = -99;
    fyHe3Gen = -99;
    fisMaterialGenHe3 = -99;
    fisPrimaryGenHe3 = -99;
    fisSecondaryGenHe3 = -99;
    fHe3Charge = -99;
    fpPGen = -99;
    fyPGen = -99;
    fisPrimaryGenP = -99;
    fisMaterialGenP = -99;
    fisSecondaryGenP = -99;
    fPCharge = -99;
    fpLambdaGen = -99;
    fyLambdaGen = -99;
    fmLambdaGen = -99;
    fLambdaCharge = -99;
    fmLambda = -99;
    fpLambda = -99;
    fptLambda = -99;
    fctLambda = -99;
    fdcaLambda = -99;
    fcosLambda = -99;
    fyLambda = -99;
    fhe3Pt = -99;
    fhe3chi2 = -99;
    fhe3Dcaz = -99;
    fhe3Dca = -99;
    fpy = -99;
    fpiy = -99;
    fhe3y = -99;
    fpLy = -99;
    fpDcaSec = -99;
    fpLDcaSec = -99;
    fpiDcaSec = -99;
    fpiDca = -99;
    fpDca = -99;
  }
 
}
