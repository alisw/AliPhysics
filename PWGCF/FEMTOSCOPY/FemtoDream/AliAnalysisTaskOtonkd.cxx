/*
 * AliAnalysisTaskOtonkd.cxx
 *
 *
 *  reCreated on: 20 Oct 2021
 *      Authors: bernhardhohlweger, Bhawani, Oton
 */

#include "AliAnalysisTaskOtonkd.h"
#include "AliFemtoDreamBasePart.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliAODInputHandler.h"
#include "AliMCEvent.h"
#include "AliLog.h"
#include "AliVEvent.h"
#include "AliNanoAODTrack.h"



ClassImp(AliAnalysisTaskOtonkd)
AliAnalysisTaskOtonkd::AliAnalysisTaskOtonkd()
  : AliAnalysisTaskSE(),
    fisLightWeight(false),
    fTrackBufferSize(),
    fIsMC(false),
    fIsMCtruth(false),
    fdoSideband(false),
    fSigmaUp(0.0),
    fSigmaLow(0.0),
    fEvent(nullptr),
    fTrack(nullptr),
    fEventCuts(nullptr),
    fTrackCutsDeuteron(nullptr),
    //fTrackCutsDeuteronMass(nullptr),
    fTrackCutsAntiDeuteron(nullptr),
    //fTrackCutsAntiDeuteronMass(nullptr),
    fTrackCutsKaon(nullptr),
    fTrackCutsAntiKaon(nullptr),
    fTrackCutsProton(nullptr),
    fTrackCutsAntiProton(nullptr),
    fConfig(nullptr),
    fPairCleaner(nullptr),
    fPartColl(nullptr),
    fGTI(nullptr),
    fEvtList(nullptr),
    fKaonList(nullptr),
    fKaonMCList(nullptr),
    fAntiKaonList(nullptr),
    fAntiKaonMCList(nullptr),
    fProtonList(nullptr),
    fProtonMCList(nullptr),
    fAntiProtonList(nullptr),
    fAntiProtonMCList(nullptr),
    fDeuteronList(nullptr),
    fDeuteronMCList(nullptr),
    fAntiDeuteronList(nullptr),
    fAntiDeuteronMCList(nullptr),
    //fDeuteronNoTOFList(nullptr),
    //fDeuteronMCNoTOFList(nullptr),
    //fAntiDeuteronNoTOFList(nullptr),
    //fAntiDeuteronMCNoTOFList(nullptr),
    fResults(nullptr),
    fResultsQA(nullptr),
      fTree(0),
    fDeuteronRestMass(nullptr),
    fAntiDeuteronRestMass(nullptr){
    //fDeuteronRestMassNoTOF(nullptr),
    //fAntiDeuteronRestMassNoTOF(nullptr){
}

AliAnalysisTaskOtonkd::AliAnalysisTaskOtonkd(
  const char *name, bool isMC, bool isMCtruth)
  : AliAnalysisTaskSE(name),
    fisLightWeight(false),
    fTrackBufferSize(2000),
    fIsMC(isMC),
    fIsMCtruth(isMCtruth),
    fdoSideband(false),
    fSigmaUp(0.0),
    fSigmaLow(0.0),
    fEvent(nullptr),
    fTrack(nullptr),
    fEventCuts(nullptr),
    fTrackCutsDeuteron(nullptr),
    fTrackCutsAntiDeuteron(nullptr),
    fTrackCutsKaon(nullptr),
    fTrackCutsAntiKaon(nullptr),
    fTrackCutsProton(nullptr),
    fTrackCutsAntiProton(nullptr),
    fConfig(nullptr),
    fPairCleaner(nullptr),
    fPartColl(nullptr),
    fGTI(nullptr),
    fEvtList(nullptr),
    fKaonList(nullptr),
    fKaonMCList(nullptr),
    fAntiKaonList(nullptr),
    fAntiKaonMCList(nullptr),
    fProtonList(nullptr),
    fProtonMCList(nullptr),
    fAntiProtonList(nullptr),
    fAntiProtonMCList(nullptr),
    fDeuteronList(nullptr),
    fDeuteronMCList(nullptr),
    fAntiDeuteronList(nullptr),
    fAntiDeuteronMCList(nullptr),
    fResults(nullptr),
    fResultsQA(nullptr),
      fTree(0),
    fDeuteronRestMass(nullptr),
    fAntiDeuteronRestMass(nullptr){
  DefineOutput(1, TList::Class());  //Output for the Event Cuts
  DefineOutput(2, TList::Class());  //Output for the Kaon Cuts
  DefineOutput(3, TList::Class());  //Output for the AntiKaon Cuts
  DefineOutput(4, TList::Class());  //Output for the Proton Cuts
  DefineOutput(5, TList::Class());  //Output for the AntiProton Cuts
  DefineOutput(6, TList::Class());  //Output for the Dueteron Cuts
  DefineOutput(7, TList::Class());  //Output for the AntiDeuteron Cuts
  DefineOutput(8, TList::Class());  //Output for the Results
  DefineOutput(9, TList::Class());  //Output for the Results QA
  DefineOutput(10, TTree::Class());  // XiTree (former OmegaTree)

}

AliAnalysisTaskOtonkd::~AliAnalysisTaskOtonkd() {
  delete fEvent;
  delete fTrack;
  delete fTrackCutsDeuteron;
  delete fTrackCutsAntiDeuteron;
  delete fTrackCutsKaon;
  delete fTrackCutsAntiKaon;
  delete fTrackCutsProton;
  delete fTrackCutsAntiProton;
  delete fPairCleaner;
  delete fPartColl;
    delete fTree;
}

Float_t AliAnalysisTaskOtonkd::GetMass2sq(
  AliFemtoDreamTrack *track)const{
  Float_t p = track->GetP();
  Float_t mass2sq = -999;
  Float_t beta = track->GetbetaTOF();
  if (!(beta < 0)) {
    mass2sq = ((1 / (beta * beta)) - 1) * (p * p);
  }
  return mass2sq;
}

float AliAnalysisTaskOtonkd::MeanTOFMassSqdDeuteron(AliFemtoDreamTrack *track) const{
  float pTVal = track->GetPt();
  float par0 =  3.55375e+00;
  float par1 = -1.25749e+00;
  float par2 = -3.60444e-01;
  float par3 = -1.00250e-01;
  float par4 = -1.00782e-02;
  return par0 + TMath::Exp(par1 * pTVal + par2 * pTVal * pTVal+ par3 * pTVal * pTVal * pTVal+ par4 * pTVal* pTVal * pTVal * pTVal);
};
float AliAnalysisTaskOtonkd::SigmaTOFMassSqdDeuteron(AliFemtoDreamTrack *track) const{
  float pTVal = track->GetPt();
  Float_t par0 = 1.19287e-02;
  Float_t par1 = 0.202460e-02;
  Float_t par2 = 1.23058e-02;//par[2];
  Float_t par3 = 30.23644e-04;
  Float_t par4 = 45.80006e-05;
  return 0.088 + 0.1*(par0 * pTVal + par1 * pTVal * pTVal + par2 * pTVal * pTVal* pTVal+ par3 * pTVal * pTVal* pTVal* pTVal+ par4 * pTVal * pTVal* pTVal* pTVal* pTVal);
};

void AliAnalysisTaskOtonkd::UserCreateOutputObjects() {


//AOD
//  fGTI = new AliAODTrack*[fTrackBufferSize];
//NANoAOD
  fGTI = new AliVTrack *[fTrackBufferSize];

  if (!fEventCuts) {
    AliError("No Event cuts \n");
  } else {
    fEventCuts->InitQA();
  }

  if (!fTrackCutsKaon) {
    AliError("No Kaon cuts \n");
  } else {
    fTrackCutsKaon->Init();
    fKaonList = fTrackCutsKaon->GetQAHists();
  }

  if (!fTrackCutsAntiKaon) {
    AliError("No AntiKaon cuts \n");
  } else {
    fTrackCutsAntiKaon->Init();
    fAntiKaonList = fTrackCutsAntiKaon->GetQAHists();
  }

  if (!fTrackCutsProton) {
    AliError("No Proton cuts \n");
  } else {
    fTrackCutsProton->Init();
    fProtonList = fTrackCutsProton->GetQAHists();
  }

  if (!fTrackCutsAntiProton) {
    AliError("No AntiProton cuts \n");
  } else {
    fTrackCutsAntiProton->Init();
    fAntiProtonList = fTrackCutsAntiProton->GetQAHists();
  }

  if (!fTrackCutsDeuteron) {
    AliError("No Deuteron cuts \n");
  } else {
    fTrackCutsDeuteron->Init();
    fDeuteronRestMass = new TH2F("fDeuteronRestMass", "Deuteron", 72, 0.5, 8.05,
                                 800, 0.00, 10.0);
    fDeuteronRestMass->GetXaxis()->SetTitle("pT(GeV)");
    fDeuteronRestMass->GetYaxis()->SetTitle("m^2(Gev)^2");
    fDeuteronList = fTrackCutsDeuteron->GetQAHists();
    fDeuteronList->Add(fDeuteronRestMass);
  }

  if (!fTrackCutsAntiDeuteron) {
    AliError("No Deuteron cuts \n");
  } else {
    fTrackCutsAntiDeuteron->Init();
    fAntiDeuteronRestMass = new TH2F("fAntiDeuteronRestMass", "AntiDeuteron", 72,
                                     0.5, 8.05, 800, 0.00, 10.0);
    fAntiDeuteronRestMass->GetXaxis()->SetTitle("pT(GeV)");
    fAntiDeuteronRestMass->GetYaxis()->SetTitle("m^2(Gev)^2");
    fAntiDeuteronList = fTrackCutsAntiDeuteron->GetQAHists();
    fAntiDeuteronList->Add(fAntiDeuteronRestMass);
  }

  if (!fConfig) {
    AliError("No Correlation Config \n");
  } else {
    fPartColl = new AliFemtoDreamPartCollection(fConfig,
        fConfig->GetMinimalBookingME());
    fPairCleaner = new AliFemtoDreamPairCleaner(2, 0,
        fConfig->GetMinimalBookingME());
  }

  fEvent = new AliFemtoDreamEvent(true, !fisLightWeight,GetCollisionCandidates(), true);
  fEvent->SetMultiplicityEstimator(fConfig->GetMultiplicityEstimator());

  fTrack = new AliFemtoDreamTrack();
  fTrack->SetUseMCInfo(fIsMC);

  if (!fEventCuts->GetMinimalBooking()) {
    fEvtList = fEventCuts->GetHistList();
  } else {
    fEvtList = new TList();
    fEvtList->SetName("EventCuts");
    fEvtList->SetOwner();
  }

  fResultsQA = new TList();
  fResultsQA->SetOwner();
  fResultsQA->SetName("ResultsQA");


  if (fConfig->GetUseEventMixing()) {
    fResults = fPartColl->GetHistList();
    if (!fConfig->GetMinimalBookingME()) {
      fResultsQA->Add(fPartColl->GetQAList());
      fResultsQA->Add(fPairCleaner->GetHistList());
    }
  } else {
    fResults = new TList();
    fResults->SetOwner();
    fResults->SetName("Results");
  }


  fTree = new TTree("oTTree","a simple TTree");
//  fTree->Branch("RunNumber",&fTRunNumber,"fTRunNumber/I");
  fTree->Branch("Vz",&fTVz,"fTVz/F");
  fTree->Branch("Mult",&fTMult,"fTMult/I");
  fTree->Branch("Spher",&fTSpher,"fTSpher/F");
  //Kaons:
  fTree->Branch("nKaon",&fTnKaon,"fTnKaon/I");
  fTree->Branch("KaonPx",&fTKaonPx,"fTKaonPx[fTnKaon]/F");
  fTree->Branch("KaonPy",&fTKaonPy,"fTKaonPy[fTnKaon]/F");
  fTree->Branch("KaonPz",&fTKaonPz,"fTKaonPz[fTnKaon]/F");
//  fTree->Branch("KaonEta",&fTKaonEta,"fTKaonEta[fTnKaon]/F");
  fTree->Branch("KaonCharge",&fTKaonCharge,"fTKaonCharge[fTnKaon]/S");
  fTree->Branch("KaonITSsigma_e",&fTKaonITSsigma_e,"fTKaonITSsigma_e[fTnKaon]/F");
  fTree->Branch("KaonTPCsigma_e",&fTKaonTPCsigma_e,"fTKaonTPCsigma_e[fTnKaon]/F");
  fTree->Branch("KaonTOFsigma_e",&fTKaonTOFsigma_e,"fTKaonTOFsigma_e[fTnKaon]/F");
  fTree->Branch("KaonITSsigma_pi",&fTKaonITSsigma_pi,"fTKaonITSsigma_pi[fTnKaon]/F");
  fTree->Branch("KaonTPCsigma_pi",&fTKaonTPCsigma_pi,"fTKaonTPCsigma_pi[fTnKaon]/F");
  fTree->Branch("KaonTOFsigma_pi",&fTKaonTOFsigma_pi,"fTKaonTOFsigma_pi[fTnKaon]/F");
  fTree->Branch("KaonITSsigma_k",&fTKaonITSsigma_k,"fTKaonITSsigma_k[fTnKaon]/F");
  fTree->Branch("KaonTPCsigma_k",&fTKaonTPCsigma_k,"fTKaonTPCsigma_k[fTnKaon]/F");
  fTree->Branch("KaonTOFsigma_k",&fTKaonTOFsigma_k,"fTKaonTOFsigma_k[fTnKaon]/F");
  fTree->Branch("KaonITSsigma_p",&fTKaonITSsigma_p,"fTKaonITSsigma_p[fTnKaon]/F");
  fTree->Branch("KaonTPCsigma_p",&fTKaonTPCsigma_p,"fTKaonTPCsigma_p[fTnKaon]/F");
  fTree->Branch("KaonTOFsigma_p",&fTKaonTOFsigma_p,"fTKaonTOFsigma_p[fTnKaon]/F");
//  fTree->Branch("KaonITSsigma_d",&fTKaonITSsigma_d,"fTKaonITSsigma_d[fTnKaon]/F");
//  fTree->Branch("KaonTPCsigma_d",&fTKaonTPCsigma_d,"fTKaonTPCsigma_d[fTnKaon]/F");
//  fTree->Branch("KaonTOFsigma_d",&fTKaonTOFsigma_d,"fTKaonTOFsigma_d[fTnKaon]/F");
//  fTree->Branch("KaonNcl",&fTKaonNcl,"fTKaonNcl[fTnKaon]/I");
//  fTree->Branch("KaonPhi",&fTKaonPhi,"fTKaonPhi[fTnKaon]/F");
  fTree->Branch("KaonDCA",&fTKaonDCA,"fTKaonDCA[fTnKaon]/F");
  fTree->Branch("KaonID",&fTKaonID,"fTKaonID[fTnKaon]/I");
  fTree->Branch("KaonITStime",&fTKaonITStime,"fTKaonITStime[fTnKaon]/O");
  fTree->Branch("KaonTOFtime",&fTKaonTOFtime,"fTKaonTOFtime[fTnKaon]/O");
//  fTree->Branch("KaonIs",&fTKaonIs,"fTKaonIs[fTnKaon]/O");
//  fTree->Branch("KaonIsFD",&fTKaonIsFD,"fTKaonIsFD[fTnKaon]/O");
//  fTree->Branch("KaonFilterBit",&fTKaonFilterBit,"fTKaonFilterBit[fTnKaon]/O");
  fTree->Branch("KaonPDG",&fTKaonPDG,"fTKaonPDG[fTnKaon]/I");

  //Deuterons:
  fTree->Branch("nDeuteron",&fTnDeuteron,"fTnDeuteron/I");
  fTree->Branch("DeuteronPx",&fTDeuteronPx,"fTDeuteronPx[fTnDeuteron]/F");
  fTree->Branch("DeuteronPy",&fTDeuteronPy,"fTDeuteronPy[fTnDeuteron]/F");
  fTree->Branch("DeuteronPz",&fTDeuteronPz,"fTDeuteronPz[fTnDeuteron]/F");
//  fTree->Branch("DeuteronEta",&fTDeuteronEta,"fTDeuteronEta[fTnDeuteron]/F");
  fTree->Branch("DeuteronCharge",&fTDeuteronCharge,"fTDeuteronCharge[fTnDeuteron]/S");
  fTree->Branch("DeuteronITSsigma_e",&fTDeuteronITSsigma_e,"fTDeuteronITSsigma_e[fTnDeuteron]/F");
  fTree->Branch("DeuteronTPCsigma_e",&fTDeuteronTPCsigma_e,"fTDeuteronTPCsigma_e[fTnDeuteron]/F");
  fTree->Branch("DeuteronTOFsigma_e",&fTDeuteronTOFsigma_e,"fTDeuteronTOFsigma_e[fTnDeuteron]/F");
  fTree->Branch("DeuteronITSsigma_pi",&fTDeuteronITSsigma_pi,"fTDeuteronITSsigma_pi[fTnDeuteron]/F");
  fTree->Branch("DeuteronTPCsigma_pi",&fTDeuteronTPCsigma_pi,"fTDeuteronTPCsigma_pi[fTnDeuteron]/F");
  fTree->Branch("DeuteronTOFsigma_pi",&fTDeuteronTOFsigma_pi,"fTDeuteronTOFsigma_pi[fTnDeuteron]/F");
  fTree->Branch("DeuteronITSsigma_k",&fTDeuteronITSsigma_k,"fTDeuteronITSsigma_k[fTnDeuteron]/F");
  fTree->Branch("DeuteronTPCsigma_k",&fTDeuteronTPCsigma_k,"fTDeuteronTPCsigma_k[fTnDeuteron]/F");
  fTree->Branch("DeuteronTOFsigma_k",&fTDeuteronTOFsigma_k,"fTDeuteronTOFsigma_k[fTnDeuteron]/F");
  fTree->Branch("DeuteronITSsigma_p",&fTDeuteronITSsigma_p,"fTDeuteronITSsigma_p[fTnDeuteron]/F");
  fTree->Branch("DeuteronTPCsigma_p",&fTDeuteronTPCsigma_p,"fTDeuteronTPCsigma_p[fTnDeuteron]/F");
  fTree->Branch("DeuteronTOFsigma_p",&fTDeuteronTOFsigma_p,"fTDeuteronTOFsigma_p[fTnDeuteron]/F");
  fTree->Branch("DeuteronITSsigma_d",&fTDeuteronITSsigma_d,"fTDeuteronITSsigma_d[fTnDeuteron]/F");
  fTree->Branch("DeuteronTPCsigma_d",&fTDeuteronTPCsigma_d,"fTDeuteronTPCsigma_d[fTnDeuteron]/F");
  fTree->Branch("DeuteronTOFsigma_d",&fTDeuteronTOFsigma_d,"fTDeuteronTOFsigma_d[fTnDeuteron]/F");
//  fTree->Branch("DeuteronNcl",&fTDeuteronNcl,"fTDeuteronNcl[fTnDeuteron]/I");
//  fTree->Branch("DeuteronPhi",&fTDeuteronPhi,"fTDeuteronPhi[fTnDeuteron]/F");
  fTree->Branch("DeuteronDCA",&fTDeuteronDCA,"fTDeuteronDCA[fTnDeuteron]/F");
  fTree->Branch("DeuteronID",&fTDeuteronID,"fTDeuteronID[fTnDeuteron]/I");
//  fTree->Branch("DeuteronTOFbeta",&fTDeuteronTOFbeta,"fTDeuteronTOFbeta[fTnDeuteron]/F");
//  fTree->Branch("DeuteronITStime",&fTDeuteronITStime,"fTDeuteronITStime[fTnDeuteron]/O");
//  fTree->Branch("DeuteronTOFtime",&fTDeuteronTOFtime,"fTDeuteronTOFtime[fTnDeuteron]/O");
  fTree->Branch("DeuteronPDG",&fTDeuteronPDG,"fTDeuteronPDG[fTnDeuteron]/I");

  PostData(1, fEvtList);
  PostData(2, fKaonList);
  PostData(3, fAntiKaonList);
  PostData(4, fProtonList);
  PostData(5, fAntiProtonList);
  PostData(6, fDeuteronList);
  PostData(7, fAntiDeuteronList);
  PostData(8, fResults);
  PostData(9, fResultsQA);
 PostData(10, fTree);
}

void AliAnalysisTaskOtonkd::UserExec(Option_t*) {
//AOD
/*
  AliAODEvent *Event = static_cast<AliAODEvent*>(fInputEvent);

  if (!Event) {
    AliWarning("No Input Event");
  } else {
    fEvent->SetEvent(Event);
    if (fEventCuts->isSelected(fEvent)) {
*/

//NANOAOD
//  AliVEvent *fInputEvent = InputEvent();
  if (!fInputEvent) {
    AliError("No input event");
    return;
  }
  fEvent->SetEvent(fInputEvent);
  if (!fEventCuts->isSelected(fEvent)) {
    return;
  }

  //Event properties for tree
  fTRunNumber = 0.; //For NanoAOD filtering trains <100, no info
  Double_t PrimVtx[3];
  fInputEvent->GetPrimaryVertex()->GetXYZ(PrimVtx);
  fTVz = PrimVtx[2];
  fTMult = fEvent->GetMultiplicity();
  fTSpher = fEvent->GetSpher();

//init tree
  for(int ii=0;ii<MAXKaonS;ii++){
   fTKaonEta[ii]=-100000.;
   fTKaonPx[ii]=-100000.;
   fTKaonPy[ii]=-100000.;
   fTKaonPz[ii]=-100000.;
   fTKaonCharge[ii]=-10;
   fTKaonDCA[ii]=-100000.;
   fTKaonNcl[ii]=-100000;
   fTKaonITStime[ii]=kFALSE;
   fTKaonTOFtime[ii]=kFALSE;
   fTKaonIs[ii]=kFALSE;
   fTKaonIsFD[ii]=kFALSE;
   fTKaonFilterBit[ii]=0;
   fTKaonPhi[ii]=-100000.;
   fTKaonID[ii]=-100000;
   fTKaonITSsigma_e[ii]=-100000.;
   fTKaonTPCsigma_e[ii]=-100000.;
   fTKaonTOFsigma_e[ii]=-100000.;
   fTKaonITSsigma_pi[ii]=-100000.;
   fTKaonTPCsigma_pi[ii]=-100000.;
   fTKaonTOFsigma_pi[ii]=-100000.;
   fTKaonITSsigma_k[ii]=-100000.;
   fTKaonTPCsigma_k[ii]=-100000.;
   fTKaonTOFsigma_k[ii]=-100000.;
   fTKaonITSsigma_p[ii]=-100000.;
   fTKaonTPCsigma_p[ii]=-100000.;
   fTKaonTOFsigma_p[ii]=-100000.;
   fTKaonITSsigma_d[ii]=-100000.;
   fTKaonTPCsigma_d[ii]=-100000.;
   fTKaonTOFsigma_d[ii]=-100000.;
   fTKaonPDG[ii]=0.;//sure Zero?
  }
  fTnKaon=0;

  for(int ii=0;ii<MAXDeuteronS;ii++){
   fTDeuteronEta[ii]=-100000.;
   fTDeuteronPx[ii]=-100000.;
   fTDeuteronPy[ii]=-100000.;
   fTDeuteronPz[ii]=-100000.;
   fTDeuteronCharge[ii]=-10;
   fTDeuteronDCA[ii]=-100000.;
   fTDeuteronNcl[ii]=-100000;
   fTDeuteronITStime[ii]=kFALSE;
   fTDeuteronTOFtime[ii]=kFALSE;
   fTDeuteronPhi[ii]=-100000.;
   fTDeuteronID[ii]=-100000;
   fTDeuteronTOFbeta[ii]=-100000;
   fTDeuteronITSsigma_e[ii]=-100000.;
   fTDeuteronTPCsigma_e[ii]=-100000.;
   fTDeuteronTOFsigma_e[ii]=-100000.;
   fTDeuteronITSsigma_pi[ii]=-100000.;
   fTDeuteronTPCsigma_pi[ii]=-100000.;
   fTDeuteronTOFsigma_pi[ii]=-100000.;
   fTDeuteronITSsigma_k[ii]=-100000.;
   fTDeuteronTPCsigma_k[ii]=-100000.;
   fTDeuteronTOFsigma_k[ii]=-100000.;
   fTDeuteronITSsigma_p[ii]=-100000.;
   fTDeuteronTPCsigma_p[ii]=-100000.;
   fTDeuteronTOFsigma_p[ii]=-100000.;
   fTDeuteronITSsigma_d[ii]=-100000.;
   fTDeuteronTPCsigma_d[ii]=-100000.;
   fTDeuteronTOFsigma_d[ii]=-100000.;
   fTDeuteronPDG[ii]=0.;//sure Zero?
  }
  fTnDeuteron=0;
 
  //define buffers for "cheap toy coalesence model"
  vector <AliAODMCParticle> p_Buff;
  vector <AliAODMCParticle> n_Buff;



/*
//AOD
      ResetGlobalTrackReference();
      for (int iTrack = 0; iTrack < Event->GetNumberOfTracks(); ++iTrack) {
        AliAODTrack *track = static_cast<AliAODTrack*>(Event->GetTrack(iTrack));
        if (!track) {
          AliFatal("No Standard AOD");
          return;
        }
        StoreGlobalTrackReference(track);
      }
*/
//NanoAOD
  // Pre-track loop to resetglobaltrackreference and storeglobaltrackreference
  ResetGlobalTrackReference();
  for (int iTrack = 0; iTrack < fInputEvent->GetNumberOfTracks(); ++iTrack) {
    AliVTrack *track = static_cast<AliVTrack *>(fInputEvent->GetTrack(iTrack));
    if (!track) {
      AliFatal("No Standard AOD");
      return;
    }
    StoreGlobalTrackReference(track);
  }

      fTrack->SetGlobalTrackInfo(fGTI, fTrackBufferSize);
      static std::vector<AliFemtoDreamBasePart> Deuterons;
      Deuterons.clear();
      static std::vector<AliFemtoDreamBasePart> AntiDeuterons;
      AntiDeuterons.clear();
      static std::vector<AliFemtoDreamBasePart> Kaons;
      Kaons.clear();
      static std::vector<AliFemtoDreamBasePart> AntiKaons;
      AntiKaons.clear();
      static std::vector<AliFemtoDreamBasePart> Protons;
      Protons.clear();
      static std::vector<AliFemtoDreamBasePart> AntiProtons;
      AntiProtons.clear();

      //Now we loop over all the tracks in the reconstructed event.
/*
//AOD
      for (int iTrack = 0; iTrack < Event->GetNumberOfTracks(); ++iTrack) {
        AliAODTrack *track = static_cast<AliAODTrack*>(Event->GetTrack(iTrack));
        if (!track) {
          AliFatal("No Standard AOD");
          return;
        }
*/
//NANoAOD
     //printf("gonna loop over tracks AGAIN %i\n", fInputEvent->GetNumberOfTracks());
     for (int iTrack = 0; iTrack < fInputEvent->GetNumberOfTracks(); ++iTrack) {
        AliVTrack   *track = static_cast<AliVTrack *>(fInputEvent->GetTrack(iTrack));
        if (!track) {
          AliFatal("No Standard AOD");
          return;
        }

//AOD:
//        fTrack->SetTrack(track);
//NANoaOD
    fTrack->SetTrack(track, fInputEvent);

    //tree selection:
    Bool_t IsKaon = kFALSE;
    Bool_t IsProton = kFALSE;
    Bool_t IsDeuteron = kFALSE;



        if (fTrackCutsKaon->isSelected(fTrack)) {
          Kaons.push_back(*fTrack);
          IsKaon= kTRUE;
          //fTKaonIs[fTnKaon] = kTRUE;
        }
        if (fTrackCutsAntiKaon->isSelected(fTrack)) {
          AntiKaons.push_back(*fTrack);
          IsKaon= kTRUE;
          //fTKaonIs[fTnKaon] = kTRUE;
        }

        if (fTrackCutsDeuteron->isSelected(fTrack)){
         Deuterons.push_back(*fTrack);
         IsDeuteron = kTRUE;
         //fTKaonIsFD[fTnKaon] = kTRUE; //for test purposes
        }
        if (fTrackCutsAntiDeuteron->isSelected(fTrack)){
         AntiDeuterons.push_back(*fTrack);
         IsDeuteron = kTRUE;
         //fTKaonIsFD[fTnKaon] = kTRUE;//for test purposes
        }

        if (fTrackCutsProton->isSelected(fTrack)){
          Protons.push_back(*fTrack);
          IsProton = kTRUE;
          //fTDeuteronIsFD[fTnDeuteron] = kTRUE;//for test purposes
        }
        if (fTrackCutsAntiProton->isSelected(fTrack)){
          AntiProtons.push_back(*fTrack);
          IsProton = kTRUE;
          //fTDeuteronIsFD[fTnDeuteron] = kTRUE;//for test purposes
        }

        //FILL kaons and deuterons
        if(IsKaon) FillKaon(fTrack);
        if(IsDeuteron) FillDeuteron(fTrack);
      }//end of track loop



//---CHEAP TOY MODEL FOR DEUTERON PRODUCTION VIA COLAESENCE (IsMCtruth) --- //
      // loop over the MC stack to get protons and neutrons for coalesence
      if(fIsMC&&fIsMCtruth) {
       AliAODInputHandler *eventHandler =
           dynamic_cast<AliAODInputHandler*>(AliAnalysisManager::GetAnalysisManager()
           ->GetInputEventHandler());
       AliMCEvent* fMC = eventHandler->MCEvent();
       //loop over mcparticles and put protons and neutrons in buffers
       for(int iPart = 0; iPart < (fMC->GetNumberOfTracks()); iPart++) {
        AliAODMCParticle *mcPart = (AliAODMCParticle*) fMC->GetTrack(iPart);
        //cout<<iPart<<"  pdg="<<mcPart->GetPdgCode()<<" mother="<<mcPart->GetMother()<<endl;
        if(mcPart->IsPhysicalPrimary()) {
         if(abs(mcPart->GetPdgCode()) == 2212
         &&iPart>1
         &&fabs(mcPart->Eta()<.9)){
          p_Buff.push_back(*mcPart);  
         }
         if(abs(mcPart->GetPdgCode()) == 2112
         &&fabs(mcPart->Eta()<.9)){
          n_Buff.push_back(*mcPart);  
         }
        }//primary
       }
       //get the smaller momentum in the center of mass (PCM)  for p-n pairs
       //Maximum PCM = 110 MeV. See . https://cds.cern.ch/record/2285500 and https://arxiv.org/pdf/2011.05898.pdf
       //the proton charge is stored in DeuteronPDG

       // CHANGES 14-Dic-2021:
       //    now change and get every pair with PCM<500MeV, store the PCM in the deuteronDCA variable
       //    also added: store the protonr charge pdg in DeuteronCharge (IMPORTANT CHANGE!!! see avobe)
       //    also added: store the proton mother pdg in DeuteronPDG
       //    also added: store the netron mother pdg in DeuteronID
       //    also added: store PCM in DeuteronDCA
       // When reading the ntuple, to selected deuterons formed by coalesence it is enough
       // to require that (for example) DeuteronTPCsigma_d<-999. (and then select the PCM via DeuteronDCA)

       //double PCMmax = 0.110;
       double PCMmax = 0.5;
       Double_t smallerpcm = 999999.;
       Double_t pd[3];
       Int_t DeuteronCharge = 0;
       for(UInt_t ip=0;ip<p_Buff.size();ip++){ //loop over protons
        int p_sign = (p_Buff.at(ip).GetPdgCode()>0)?1:-1;
        for(UInt_t in=0;in<n_Buff.size();in++){ //loop over neutrons
         if(p_sign == (n_Buff.at(in).GetPdgCode()>0)?1:-1 ){
          Double_t p1x = p_Buff.at(ip).Px();
          Double_t p1y = p_Buff.at(ip).Py();
          Double_t p1z = p_Buff.at(ip).Pz();
          //cout<<" the mother is "<<p_Buff.at(ip).GetMother()<<" the sign is "<<p_sign<<endl;
          Int_t protonmotherpdg = -99999;
          if(p_Buff.at(ip).GetMother()<fMC->GetNumberOfTracks()){
           AliAODMCParticle *mcMother = (AliAODMCParticle*) fMC->GetTrack(p_Buff.at(ip).GetMother());
           protonmotherpdg = mcMother->GetPdgCode();
           mcMother->Delete();
          }
          Double_t m1 = 0.938272088;//proton mass
          Double_t p2x = n_Buff.at(in).Px();
          Double_t p2y = n_Buff.at(in).Py();
          Double_t p2z = n_Buff.at(in).Pz();
          Int_t neutronmotherpdg = -99999;
          if(n_Buff.at(in).GetMother()<fMC->GetNumberOfTracks()){
           AliAODMCParticle *mcMother = (AliAODMCParticle*) fMC->GetTrack(n_Buff.at(in).GetMother());
           neutronmotherpdg = mcMother->GetPdgCode();
           mcMother->Delete();
          }
          Double_t m2 = 0.939565420;//neutron mass
	  Double_t E1 = TMath::Sqrt( p1x*p1x + p1y*p1y + p1z*p1z + m1*m1);
	  Double_t E2 = TMath::Sqrt( p2x*p2x + p2y*p2y + p2z*p2z + m2*m2);
	  Double_t s = (E1+E2)*(E1+E2) - ((p1x+p2x)*(p1x+p2x) + (p1y+p2y)*(p1y+p2y) + (p1z+p2z)*(p1z+p2z));
	  Double_t pcm = TMath::Sqrt( (s-(m1-m2)*(m1-m2))*(s-(m1+m2)*(m1+m2)) )/(2.*TMath::Sqrt(s));
          //if(pcm<PCMmax&&pcm<smallerpcm){
          if(pcm<PCMmax){
           //smallerpcm=pcm;
           DeuteronCharge = p_sign;
           pd[0]=p1x+p2x;
           pd[1]=p1y+p2y;
           pd[2]=p1z+p2z;
           fTDeuteronPx[fTnDeuteron] = pd[0];
           fTDeuteronPy[fTnDeuteron] = pd[1];
           fTDeuteronPz[fTnDeuteron] = pd[2];
           fTDeuteronCharge[fTnDeuteron] = DeuteronCharge;
           fTDeuteronPDG[fTnDeuteron] =  protonmotherpdg;//proton pdg
           fTDeuteronID[fTnDeuteron] =  neutronmotherpdg;//neutron pdg
           fTDeuteronDCA[fTnDeuteron] = pcm;
           fTnDeuteron++;
          }
         }//sign check
        }
       }
//       if(DeuteronCharge!=0){
//        fTDeuteronPx[fTnDeuteron] = pd[0];
//        fTDeuteronPy[fTnDeuteron] = pd[1];
//        fTDeuteronPz[fTnDeuteron] = pd[2];
//        fTDeuteronPDG[fTnDeuteron] = DeuteronCharge;
//        fTnDeuteron++;
//       }
       p_Buff.clear();
       n_Buff.clear();
      }//IsMCtruth
//---coalesence cheap toy model end --- //



      //fill tree:
      //----------
      //if(fTnKaon>0) fTree->Fill(); //for test
      if(fTnKaon>0&&fTnDeuteron>0) fTree->Fill(); //kd 


     bool FemtoDreamPairing = false; // Skip FD pairing/mixing for now (to save computing time)
     if(FemtoDreamPairing){
      fPairCleaner->CleanTrackAndDecay(&Kaons, &Deuterons, 0);///NOT SURE AT ALL ABOUT THIS 0 and 1
      fPairCleaner->CleanTrackAndDecay(&AntiKaons, &AntiDeuterons, 1);///NOT SURE AT ALL ABOUT THIS 0 and 1
      fPairCleaner->CleanTrackAndDecay(&Kaons, &Protons, 0);///NOT SURE AT ALL ABOUT THIS 0 and 1
      fPairCleaner->CleanTrackAndDecay(&AntiKaons, &AntiProtons, 1);///NOT SURE AT ALL ABOUT THIS 0 and 1
      fPairCleaner->CleanTrackAndDecay(&Protons, &Deuterons, 0);
      fPairCleaner->CleanTrackAndDecay(&AntiProtons, &AntiDeuterons, 1);
      fPairCleaner->ResetArray();
      fPairCleaner->StoreParticle(Kaons);
      fPairCleaner->StoreParticle(AntiKaons);
      fPairCleaner->StoreParticle(Protons);
      fPairCleaner->StoreParticle(AntiProtons);
      fPairCleaner->StoreParticle(Deuterons);
      fPairCleaner->StoreParticle(AntiDeuterons);

      fPartColl->SetEvent(fPairCleaner->GetCleanParticles(),
                          fEvent->GetZVertex(), fEvent->GetRefMult08(),
                          fEvent->GetV0MCentrality());
  
//??? will this work with NANOAOD ????: (apparently it does)
      void SetEvent(std::vector<AliFemtoDreamBasePart> &vec1,
                    std::vector<AliFemtoDreamBasePart> &vec2,
                    AliFemtoDreamEvent * evt, const int pdg1, const int pdg2);
     }

  PostData(1, fEvtList);
  PostData(2, fKaonList);
  PostData(3, fAntiKaonList);
  PostData(4, fProtonList);
  PostData(5, fAntiProtonList);
  PostData(6, fDeuteronList);
  PostData(7, fAntiDeuteronList);
  PostData(8, fResults);
  PostData(9, fResultsQA);
  PostData(10, fTree);
}

void AliAnalysisTaskOtonkd::ResetGlobalTrackReference() {

  for (UShort_t i = 0; i < fTrackBufferSize; i++) {
    fGTI[i] = 0;
  }
}

//AOD
/*
void AliAnalysisTaskOtonkd::StoreGlobalTrackReference(
  AliAODTrack *track) {

  const int trackID = track->GetID();
  if (trackID < 0) {
    return;
  }

  if (trackID >= fTrackBufferSize) {
    printf("Warning: track ID too big for buffer: ID: %d, buffer %d\n", trackID,
           fTrackBufferSize);
    return;
  }

  if (fGTI[trackID]) {

    if ((!track->GetFilterMap()) && (!track->GetTPCNcls())) {
      return;
    }

    if (fGTI[trackID]->GetFilterMap() || fGTI[trackID]->GetTPCNcls()) {
      // If we come here, there's a problem
      printf("Warning! global track info already there!");
      printf("         TPCNcls track1 %u track2 %u",
             (fGTI[trackID])->GetTPCNcls(), track->GetTPCNcls());
      printf("         FilterMap track1 %u track2 %u\n",
             (fGTI[trackID])->GetFilterMap(), track->GetFilterMap());
    }
  }
  (fGTI[trackID]) = track;
}
*/
//NANoAOD
//____________________________________________________________________________________________________
void AliAnalysisTaskOtonkd::StoreGlobalTrackReference(AliVTrack *track) {
  // see AliFemtoDreamAnalysis for details
  AliNanoAODTrack *nanoTrack = dynamic_cast<AliNanoAODTrack *>(track);
  const int trackID = track->GetID();
  if (trackID < 0) {
    return;
  }
  if (trackID >= fTrackBufferSize) {
    printf("Warning: track ID too big for buffer: ID: %d, buffer %d\n", trackID,
           fTrackBufferSize);
            return;
  }

  if (fGTI[trackID]) {
    if ((!nanoTrack->GetFilterMap()) && (!track->GetTPCNcls())) {
      return;
    }
    if (dynamic_cast<AliNanoAODTrack *>(fGTI[trackID])->GetFilterMap()
        || fGTI[trackID]->GetTPCNcls()) {
      printf("Warning! global track info already there!");
      printf("         TPCNcls track1 %u track2 %u",
             (fGTI[trackID])->GetTPCNcls(), track->GetTPCNcls());
      printf("         FilterMap track1 %u track2 %u\n",
             dynamic_cast<AliNanoAODTrack *>(fGTI[trackID])->GetFilterMap(),
             nanoTrack->GetFilterMap());
    }
  }
  (fGTI[trackID]) = track;
}

//________________________________________________________________________________
Bool_t AliAnalysisTaskOtonkd::FillKaon(AliFemtoDreamTrack *TheTrack) {
 Bool_t Filled = kFALSE;
 TVector3 mom;
 mom = TheTrack->GetMomentum();
 fTKaonPx[fTnKaon] = mom.X();
 fTKaonPy[fTnKaon] = mom.Y();
 fTKaonPz[fTnKaon] = mom.Z();
 fTKaonEta[fTnKaon] = TheTrack->GetEta().at(0);
 fTKaonCharge[fTnKaon] = TheTrack->GetCharge().at(0);
 fTKaonNcl[fTnKaon] = TheTrack->GetNClsTPC();
 fTKaonPhi[fTnKaon] = (TheTrack->GetPhiAtRaidius().at(0)).at(0);//phi for r=85.cm ???
 fTKaonDCA[fTnKaon] = TheTrack->GetDCAXYProp();
 fTKaonID[fTnKaon] = TheTrack->GetIDTracks().at(0);
 fTKaonITSsigma_e[fTnKaon] = (TheTrack->GetnSigmaITS((int) (AliPID::kElectron)));
 fTKaonTPCsigma_e[fTnKaon] = (TheTrack->GetnSigmaTPC((int) (AliPID::kElectron)));
 fTKaonTOFsigma_e[fTnKaon] = (TheTrack->GetnSigmaTOF((int) (AliPID::kElectron)));
 fTKaonITSsigma_pi[fTnKaon] = (TheTrack->GetnSigmaITS((int) (AliPID::kPion)));
 fTKaonTPCsigma_pi[fTnKaon] = (TheTrack->GetnSigmaTPC((int) (AliPID::kPion)));
 fTKaonTOFsigma_pi[fTnKaon] = (TheTrack->GetnSigmaTOF((int) (AliPID::kPion)));
 fTKaonITSsigma_k[fTnKaon] = (TheTrack->GetnSigmaITS((int) (AliPID::kKaon)));
 fTKaonTPCsigma_k[fTnKaon] = (TheTrack->GetnSigmaTPC((int) (AliPID::kKaon)));
 fTKaonTOFsigma_k[fTnKaon] = (TheTrack->GetnSigmaTOF((int) (AliPID::kKaon)));
 fTKaonITSsigma_p[fTnKaon] = (TheTrack->GetnSigmaITS((int) (AliPID::kProton)));
 fTKaonTPCsigma_p[fTnKaon] = (TheTrack->GetnSigmaTPC((int) (AliPID::kProton)));
 fTKaonTOFsigma_p[fTnKaon] = (TheTrack->GetnSigmaTOF((int) (AliPID::kProton)));
 fTKaonITSsigma_d[fTnKaon] = (TheTrack->GetnSigmaITS((int) (AliPID::kDeuteron)));
 fTKaonTPCsigma_d[fTnKaon] = (TheTrack->GetnSigmaTPC((int) (AliPID::kDeuteron)));
 fTKaonTOFsigma_d[fTnKaon] = (TheTrack->GetnSigmaTOF((int) (AliPID::kDeuteron)));
 fTKaonITStime[fTnKaon] = TheTrack->GetHasITSHit();
 fTKaonTOFtime[fTnKaon] = TheTrack->GetTOFTimingReuqirement();
  fTKaonFilterBit[fTnKaon] = TheTrack->GetFilterMap();
 fTKaonPDG[fTnKaon] = TheTrack->GetMCPDGCode();

 fTnKaon++;
 Filled = kTRUE;
 return Filled;
}



//________________________________________________________________________________
Bool_t AliAnalysisTaskOtonkd::FillDeuteron(AliFemtoDreamTrack *TheTrack) {
 Bool_t Filled = kFALSE;
 TVector3 mom;
 mom = TheTrack->GetMomentum();
 fTDeuteronPx[fTnDeuteron] = mom.X();
 fTDeuteronPy[fTnDeuteron] = mom.Y();
 fTDeuteronPz[fTnDeuteron] = mom.Z();
 fTDeuteronEta[fTnDeuteron] = TheTrack->GetEta().at(0);
 fTDeuteronCharge[fTnDeuteron] = TheTrack->GetCharge().at(0);
 fTDeuteronNcl[fTnDeuteron] = TheTrack->GetNClsTPC();
 fTDeuteronPhi[fTnDeuteron] = (TheTrack->GetPhiAtRaidius().at(0)).at(0);//phi for r=85.cm ???
 fTDeuteronDCA[fTnDeuteron] = TheTrack->GetDCAXYProp();
 fTDeuteronID[fTnDeuteron] = TheTrack->GetIDTracks().at(0);
 fTDeuteronTOFbeta[fTnDeuteron] = TheTrack->GetbetaTOF();
 fTDeuteronITSsigma_e[fTnDeuteron] = (TheTrack->GetnSigmaITS((int) (AliPID::kElectron)));
 fTDeuteronTPCsigma_e[fTnDeuteron] = (TheTrack->GetnSigmaTPC((int) (AliPID::kElectron)));
 fTDeuteronTOFsigma_e[fTnDeuteron] = (TheTrack->GetnSigmaTOF((int) (AliPID::kElectron)));
 fTDeuteronITSsigma_pi[fTnDeuteron] = (TheTrack->GetnSigmaITS((int) (AliPID::kPion)));
 fTDeuteronTPCsigma_pi[fTnDeuteron] = (TheTrack->GetnSigmaTPC((int) (AliPID::kPion)));
 fTDeuteronTOFsigma_pi[fTnDeuteron] = (TheTrack->GetnSigmaTOF((int) (AliPID::kPion)));
 fTDeuteronITSsigma_k[fTnDeuteron] = (TheTrack->GetnSigmaITS((int) (AliPID::kKaon)));
 fTDeuteronTPCsigma_k[fTnDeuteron] = (TheTrack->GetnSigmaTPC((int) (AliPID::kKaon)));
 fTDeuteronTOFsigma_k[fTnDeuteron] = (TheTrack->GetnSigmaTOF((int) (AliPID::kKaon)));
 fTDeuteronITSsigma_p[fTnDeuteron] = (TheTrack->GetnSigmaITS((int) (AliPID::kProton)));
 fTDeuteronTPCsigma_p[fTnDeuteron] = (TheTrack->GetnSigmaTPC((int) (AliPID::kProton)));
 fTDeuteronTOFsigma_p[fTnDeuteron] = (TheTrack->GetnSigmaTOF((int) (AliPID::kProton)));
 fTDeuteronITSsigma_d[fTnDeuteron] = (TheTrack->GetnSigmaITS((int) (AliPID::kDeuteron)));
 fTDeuteronTPCsigma_d[fTnDeuteron] = (TheTrack->GetnSigmaTPC((int) (AliPID::kDeuteron)));
 fTDeuteronTOFsigma_d[fTnDeuteron] = (TheTrack->GetnSigmaTOF((int) (AliPID::kDeuteron)));
 fTDeuteronITStime[fTnDeuteron] = TheTrack->GetHasITSHit();
 fTDeuteronTOFtime[fTnDeuteron] = TheTrack->GetTOFTimingReuqirement();
 fTDeuteronPDG[fTnDeuteron] = TheTrack->GetMCPDGCode();

 fTnDeuteron++;
 Filled = kTRUE;
 return Filled;
}
