/*
 * AliAnalysisTaskOtonXx.cxx
 *
 *
 *  reCreated on: 20 Oct 2021
 *      Authors: bernhardhohlweger, Bhawani, Oton
 */

#include "AliAnalysisTaskOtonXx.h"
#include "AliFemtoDreamBasePart.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliAODInputHandler.h"
#include "AliMCEvent.h"
#include "AliLog.h"
#include "AliVEvent.h"
#include "AliNanoAODTrack.h"



ClassImp(AliAnalysisTaskOtonXx)
AliAnalysisTaskOtonXx::AliAnalysisTaskOtonXx()
  : AliAnalysisTaskSE(),
    fisLightWeight(false),
    fTrackBufferSize(),
    fIsMC(false),
    fIsMCtruth(false),
    fdoFDpairing(false),
    fisOmega(false),
    fisPi(false),
    fOnlyXi(false),
    fEvent(nullptr),
    fTrack(nullptr),
    fEventCuts(nullptr),
    fEvtList(nullptr),
    fTrackCutsKaon(nullptr),
    fKaonList(nullptr),
    fTrackCutsAntiKaon(nullptr),
    fAntiKaonList(nullptr),
    fCascade(nullptr),
    fCutsXi(nullptr),
    fXiList(nullptr),
    fCutsAntiXi(nullptr),
    fAntiXiList(nullptr),
    fConfig(nullptr),
    fPairCleaner(nullptr),
    fPartColl(nullptr),
    fGTI(nullptr),
    fResults(nullptr),
    fResultsQA(nullptr),
      fTree(0) {
}

AliAnalysisTaskOtonXx::AliAnalysisTaskOtonXx(const char *name, bool doFDpairing, bool isMC, bool isMCtruth, bool isOmega, bool isPi, bool OnlyXi)
  : AliAnalysisTaskSE(name),
    fisLightWeight(false),
    fTrackBufferSize(2000),
    fIsMC(isMC),
    fIsMCtruth(isMCtruth),
    fdoFDpairing(doFDpairing),
    fisOmega(isOmega),
    fisPi(isPi),
    fOnlyXi(OnlyXi),
    fEvent(nullptr),
    fTrack(nullptr),
    fEventCuts(nullptr),
    fEvtList(nullptr),
    fTrackCutsKaon(nullptr),
    fKaonList(nullptr),
    fTrackCutsAntiKaon(nullptr),
    fAntiKaonList(nullptr),
    fCascade(nullptr),
    fCutsXi(nullptr),
    fXiList(nullptr),
    fCutsAntiXi(nullptr),
    fAntiXiList(nullptr),
    fConfig(nullptr),
    fPairCleaner(nullptr),
    fPartColl(nullptr),
    fGTI(nullptr),
    fResults(nullptr),
    fResultsQA(nullptr),
      fTree(0){
  DefineOutput(1, TList::Class());  //Output for the Event Cuts
  DefineOutput(2, TList::Class());  //Output for the Kaon Cuts
  DefineOutput(3, TList::Class());  //Output for the AntiKaon Cuts
  DefineOutput(4, TList::Class());  //Output for the Xi Cuts
  DefineOutput(5, TList::Class());  //Output for the AntiXi Cuts
  DefineOutput(6, TList::Class());  //Output for the Results
  DefineOutput(7, TList::Class());  //Output for the Results QA
  DefineOutput(8, TTree::Class());  // XiTree (former OmegaTree)

}

AliAnalysisTaskOtonXx::~AliAnalysisTaskOtonXx() {
  if (fEvent) {
    delete fEvent;
  }
  if (fEventCuts) {
    delete fEventCuts;
  }
  if (fTrack) {
    delete fTrack;
  }
  if (fTrackCutsKaon) {
    delete fTrackCutsKaon;
  }
  if (fTrackCutsAntiKaon) {
    delete fTrackCutsAntiKaon;
  }
  if (fCascade) {
    delete fCascade;
  }
  if (fCutsXi) {
    delete fCutsXi;
  }
  if (fCutsAntiXi) {
    delete fCutsAntiXi;
  }
  if (fPairCleaner) {
    delete fPairCleaner;
  }
  if (fPartColl) {
    delete fPartColl;
  }
  if (fTree) {
    delete fTree;
  }
}

void AliAnalysisTaskOtonXx::UserCreateOutputObjects() {


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
  }
  if (!fTrackCutsAntiKaon) {
    AliError("No AntiKaon cuts \n");
  } else {
    fTrackCutsAntiKaon->Init();
  }
  if (!fCutsXi) {
    AliError("No Xi cuts \n");
  } else {
    fCutsXi->Init();
  }

  if (!fCutsAntiXi) {
    AliError("No Xi cuts \n");
  } else {
    fCutsAntiXi->Init();
  }


// THIS SHOULD BE SET UP AT SOME POINT:
//  pair cleaner config ??????????????
  if (!fConfig) {
    AliError("No Correlation Config \n");
  } else {
    fPartColl = new AliFemtoDreamPartCollection(fConfig,fConfig->GetMinimalBookingME());
    //for now, following the same approach as p-Omega:
    fPairCleaner = new AliFemtoDreamPairCleaner(2, 2,fConfig->GetMinimalBookingME()); //??????????????????
  }

// THIS SHOULD BE SET UP AT SOME POINT:
  fEvent = new AliFemtoDreamEvent(true, !fisLightWeight,GetCollisionCandidates(), true);  // true...true ???????
  fEvent->SetMultiplicityEstimator(fConfig->GetMultiplicityEstimator());

  fTrack = new AliFemtoDreamTrack();
  fTrack->SetUseMCInfo(fIsMC);


  fCascade = new AliFemtoDreamCascade();
  fCascade->SetUseMCInfo(fIsMC); 
  //PDG Codes should be set assuming Xi- to also work for Xi+
  if(!fisOmega){
   fCascade->SetPDGCode(3312);
  }else{
   fCascade->SetPDGCode(3334);
  }
  fCascade->SetPDGDaugPos(2212);            //Proton
  fCascade->SetPDGDaugNeg(-211);             //pi^-
  if(!fisOmega){
   fCascade->SetPDGDaugBach(-211);            //EDIT pi^-
  }else{
   fCascade->SetPDGDaugBach(-321);            //EDIT pi^-
  }
  fCascade->Setv0PDGCode(3122);    // Lambda


  //lists
  if (!fEventCuts->GetMinimalBooking()) {
    fEvtList = fEventCuts->GetHistList();
  } else {
    fEvtList = new TList();
    fEvtList->SetName("EventCuts");
    fEvtList->SetOwner();
  }

  if (!fTrackCutsKaon->GetMinimalBooking()) {
    fKaonList = fTrackCutsKaon->GetQAHists();
  } else {
    fKaonList = new TList();
    fKaonList->SetName("TrackCutsKaon");
    fKaonList->SetOwner();
  }
  if (!fTrackCutsAntiKaon->GetMinimalBooking()) {
    fAntiKaonList = fTrackCutsAntiKaon->GetQAHists();
  } else {
    fAntiKaonList = new TList();
    fAntiKaonList->SetName("TrackCutsAntiKaon");
    fAntiKaonList->SetOwner();
  }

  if (!fCutsXi->GetMinimalBooking()) {
    fXiList = fCutsXi->GetQAHists();
  } else {
    fXiList = new TList();
    fXiList->SetName("XiCuts");
    fXiList->SetOwner();
  }
  if (!fCutsAntiXi->GetMinimalBooking()) {
    fAntiXiList = fCutsAntiXi->GetQAHists();
  } else {
    fAntiXiList = new TList();
    fAntiXiList->SetName("AntiXiCuts");
    fAntiXiList->SetOwner();
  }

  //results
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

  // tree booking
  fTree = new TTree("oTTree","a simple TTree");
//  fTree->Branch("RunNumber",&fTRunNumber,"fTRunNumber/I");
  fTree->Branch("Vz",&fTVz,"fTVz/F");
  fTree->Branch("Mult",&fTMult,"fTMult/I");
  fTree->Branch("Spher",&fTSpher,"fTSpher/F");
  //Kaons:
  fTree->Branch("nKaon",&fTnKaon,"fTnKaon/I");
  if(!fOnlyXi){
   fTree->Branch("KaonPx",&fTKaonPx,"fTKaonPx[fTnKaon]/F");
   fTree->Branch("KaonPy",&fTKaonPy,"fTKaonPy[fTnKaon]/F");
   fTree->Branch("KaonPz",&fTKaonPz,"fTKaonPz[fTnKaon]/F");
//  if(!fIsMCtruth)fTree->Branch("KaonPTPC",&fTKaonPTPC,"fTKaonPTPC[fTnKaon]/F");
   if(!fIsMCtruth)fTree->Branch("KaonCharge",&fTKaonCharge,"fTKaonCharge[fTnKaon]/S");
//  fTree->Branch("KaonITSsigma_e",&fTKaonITSsigma_e,"fTKaonITSsigma_e[fTnKaon]/F");
//  if(!fIsMCtruth)fTree->Branch("KaonTPCsigma_e",&fTKaonTPCsigma_e,"fTKaonTPCsigma_e[fTnKaon]/F");
//  fTree->Branch("KaonTOFsigma_e",&fTKaonTOFsigma_e,"fTKaonTOFsigma_e[fTnKaon]/F");
//  fTree->Branch("KaonITSsigma_pi",&fTKaonITSsigma_pi,"fTKaonITSsigma_pi[fTnKaon]/F");
//  if(!fIsMCtruth)fTree->Branch("KaonTPCsigma_pi",&fTKaonTPCsigma_pi,"fTKaonTPCsigma_pi[fTnKaon]/F");
//  if(!fIsMCtruth)fTree->Branch("KaonTOFsigma_pi",&fTKaonTOFsigma_pi,"fTKaonTOFsigma_pi[fTnKaon]/F");
//  fTree->Branch("KaonITSsigma_k",&fTKaonITSsigma_k,"fTKaonITSsigma_k[fTnKaon]/F");
//  if(!fIsMCtruth)fTree->Branch("KaonTPCsigma_k",&fTKaonTPCsigma_k,"fTKaonTPCsigma_k[fTnKaon]/F");
//  if(!fIsMCtruth)fTree->Branch("KaonTOFsigma_k",&fTKaonTOFsigma_k,"fTKaonTOFsigma_k[fTnKaon]/F");
//  fTree->Branch("KaonITSsigma_p",&fTKaonITSsigma_p,"fTKaonITSsigma_p[fTnKaon]/F");
//  if(!fIsMCtruth)fTree->Branch("KaonTPCsigma_p",&fTKaonTPCsigma_p,"fTKaonTPCsigma_p[fTnKaon]/F");
//  if(!fIsMCtruth)fTree->Branch("KaonTOFsigma_p",&fTKaonTOFsigma_p,"fTKaonTOFsigma_p[fTnKaon]/F");
//  fTree->Branch("KaonITSsigma_d",&fTKaonITSsigma_d,"fTKaonITSsigma_d[fTnKaon]/F");
//  fTree->Branch("KaonTPCsigma_d",&fTKaonTPCsigma_d,"fTKaonTPCsigma_d[fTnKaon]/F");
//  fTree->Branch("KaonTOFsigma_d",&fTKaonTOFsigma_d,"fTKaonTOFsigma_d[fTnKaon]/F");
//  if(!fIsMCtruth)fTree->Branch("KaonNcl",&fTKaonNcl,"fTKaonNcl[fTnKaon]/I");
//  fTree->Branch("KaonPhi",&fTKaonPhi,"fTKaonPhi[fTnKaon]/F");
   if(!fIsMCtruth)fTree->Branch("KaonDCA",&fTKaonDCA,"fTKaonDCA[fTnKaon]/F");
//  if(!fIsMCtruth)fTree->Branch("KaonDCAz",&fTKaonDCAz,"fTKaonDCAz[fTnKaon]/F");
   fTree->Branch("KaonID",&fTKaonID,"fTKaonID[fTnKaon]/I"); //used also in MCtruth
   if(!fIsMCtruth)fTree->Branch("KaonSPDtime",&fTKaonSPDtime,"fTKaonSPDtime[fTnKaon]/O");
   if(!fIsMCtruth)fTree->Branch("KaonITStime",&fTKaonITStime,"fTKaonITStime[fTnKaon]/O");
   if(!fIsMCtruth)fTree->Branch("KaonTOFtime",&fTKaonTOFtime,"fTKaonTOFtime[fTnKaon]/O");
//  fTree->Branch("KaonIs",&fTKaonIs,"fTKaonIs[fTnKaon]/O");
//  fTree->Branch("KaonIsFD",&fTKaonIsFD,"fTKaonIsFD[fTnKaon]/O");
//  fTree->Branch("KaonFilterBit",&fTKaonFilterBit,"fTKaonFilterBit[fTnKaon]/O");
//  if(fIsMC||fIsMCtruth)fTree->Branch("KaonPDG",&fTKaonPDG,"fTKaonPDG[fTnKaon]/I");
//  if(fIsMC||fIsMCtruth)fTree->Branch("KaonMotherWeak",&fTKaonMotherWeak,"fTKaonMotherWeak[fTnKaon]/I");
//  if(fIsMC||fIsMCtruth)fTree->Branch("KaonOrigin",&fTKaonOrigin,"fTKaonOrigin[fTnKaon]/I");
   if(fIsMC||fIsMCtruth)fTree->Branch("KaonMotherID",&fTKaonMotherID,"fTKaonMotherID[fTnKaon]/I");
  }


  //Xis:
  fTree->Branch("nXi",&fTnXi,"fTnXi/I");
  fTree->Branch("XiCharge",&fTXiCharge,"fTXiCharge[fTnXi]/S");
  fTree->Branch("XiDCA",&fTXiDCA,"fTXiDCA[fTnXi]/F");
  fTree->Branch("XiDaughtersDCA",&fTXiDaughtersDCA,"fTXiDaughtersDCA[fTnXi]/F");
  fTree->Branch("XiMass",&fTXiMass,"fTXiMass[fTnXi]/F");
  if(fOnlyXi) fTree->Branch("XiOmegaMass",&fTXiOmegaMass,"fTXiOmegaMass[fTnXi]/F");
  if(fOnlyXi) fTree->Branch("XiVr",&fTXiVr,"fTXiVr[fTnXi]/F");
  fTree->Branch("XiPA",&fTXiPA,"fTXiPA[fTnXi]/F");
  if(fOnlyXi) fTree->Branch("XiLambdaDCA",&fTXiLambdaDCA,"fTXiLambdaDCA[fTnXi]/F");
  if(fOnlyXi) fTree->Branch("XiLambdaDaughtersDCA",&fTXiLambdaDaughtersDCA,"fTXiLambdaDaughtersDCA[fTnXi]/F");
  fTree->Branch("XiLambdaMass",&fTXiLambdaMass,"fTXiLambdaMass[fTnXi]/F");
  if(fOnlyXi) fTree->Branch("XiLambdaK0Mass",&fTXiLambdaK0Mass,"fTXiLambdaK0Mass[fTnXi]/F");
  fTree->Branch("XiLambdaVr",&fTXiLambdaVr,"fTXiLambdaVr[fTnXi]/F");
  if(fOnlyXi) fTree->Branch("XiLambdaPA",&fTXiLambdaPA,"fTXiLambdaPA[fTnXi]/F");
  fTree->Branch("XiTrackCharge",&fTXiTrackCharge,"fTXiTrackCharge[fTnXi][3]/S");
  fTree->Branch("XiTrackPx",&fTXiTrackPx,"fTXiTrackPx[fTnXi][3]/F");
  fTree->Branch("XiTrackPy",&fTXiTrackPy,"fTXiTrackPy[fTnXi][3]/F");
  fTree->Branch("XiTrackPz",&fTXiTrackPz,"fTXiTrackPz[fTnXi][3]/F");
  if(fOnlyXi) fTree->Branch("XiTrackTPCmom",&fTXiTrackTPCmom,"fTXiTrackTPCmom[fTnXi][3]/F");
  if(fOnlyXi) fTree->Branch("XiTrackDCA",&fTXiTrackDCA,"fTXiTrackDCA[fTnXi][3]/F");
  if(fOnlyXi) fTree->Branch("XiTrackTPCsigma",&fTXiTrackTPCsigma,"fTXiTrackTPCsigma[fTnXi][3]/F");
  if(fOnlyXi) fTree->Branch("XiTrackTOFsigma",&fTXiTrackTOFsigma,"fTXiTrackTOFsigma[fTnXi][3]/F");
  if(fOnlyXi) fTree->Branch("XiTrackNcl",&fTXiTrackNcl,"fTXiTrackNcl[fTnXi][3]/I");
  if(fOnlyXi) fTree->Branch("XiTrackCrR",&fTXiTrackCrR,"fTXiTrackCrR[fTnXi][3]/F");
  if(fOnlyXi) fTree->Branch("XiTrackCrF",&fTXiTrackCrF,"fTXiTrackCrF[fTnXi][3]/F");
  fTree->Branch("XiTrackSPDtime",&fTXiTrackSPDtime,"fTXiTrackSPDtime[fTnXi][3]/O");
  fTree->Branch("XiTrackITStime",&fTXiTrackITStime,"fTXiTrackITStime[fTnXi][3]/O");
  fTree->Branch("XiTrackTOFtime",&fTXiTrackTOFtime,"fTXiTrackTOFtime[fTnXi][3]/O");
  fTree->Branch("XiTrackID",&fTXiTrackID,"fTXiTrackID[fTnXi][3]/I");
  if(fIsMC||fIsMCtruth) fTree->Branch("XiMotherID",&fTXiMotherID,"fTXiMotherID[fTnXi]/I");
  if(fIsMC||fIsMCtruth) fTree->Branch("XiPDG",&fTXiPDG,"fTXiPDG[fTnXi]/I");
  if(fIsMC||fIsMCtruth) fTree->Branch("XiMotherPDG",&fTXiMotherPDG,"fTXiMotherPDG[fTnXi]/I");
  if(fIsMC||fIsMCtruth) fTree->Branch("XiMotherWeak",&fTXiMotherWeak,"fTXiMotherWeak[fTnXi]/I");
  if(fIsMC||fIsMCtruth) fTree->Branch("XiOrigin",&fTXiOrigin,"fTXiOrigin[fTnXi]/I");

  PostData(1, fEvtList);
  PostData(2, fKaonList);
  PostData(3, fAntiKaonList);
  PostData(4, fXiList);
  PostData(5, fAntiXiList);
  PostData(6, fResults);
  PostData(7, fResultsQA);
  PostData(8, fTree);
}

void AliAnalysisTaskOtonXx::UserExec(Option_t*) {
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
   fTKaonPx[ii]=-100000.;
   fTKaonPy[ii]=-100000.;
   fTKaonPz[ii]=-100000.;
   fTKaonPTPC[ii]=-100000.;
   fTKaonCharge[ii]=-10;
   fTKaonDCA[ii]=-100000.;
   fTKaonDCAz[ii]=-100000.;
   fTKaonNcl[ii]=-100000;
   fTKaonSPDtime[ii]=kFALSE;
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
   fTKaonMotherWeak[ii]=0.;//sure Zero?
   fTKaonOrigin[ii]=-1;
   fTKaonMotherID[ii]=-1;
  }
  fTnKaon=0;



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
      static std::vector<AliFemtoDreamBasePart> Kaons;
      Kaons.clear();
      static std::vector<AliFemtoDreamBasePart> AntiKaons;
      AntiKaons.clear();

      //boolean for MC for Kaon DCA templates
      //Bool_t IsKaonWeak = kFALSE;

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


        if (fTrackCutsKaon->isSelected(fTrack)) {
          Kaons.push_back(*fTrack);
          IsKaon= kTRUE;
          //if(fTrack->GetMotherWeak()!=0) IsKaonWeak = kTRUE;
        }
        if (fTrackCutsAntiKaon->isSelected(fTrack)) {
          AntiKaons.push_back(*fTrack);
          IsKaon= kTRUE;
          //if(fTrack->GetMotherWeak()!=0) IsKaonWeak = kTRUE;
        }

        if(IsKaon) FillKaon(fTrack);


      }//end of track loop


//init cascade xi tree
  for(int ii=0;ii<MAXXiS;ii++){
   fTXiCharge[ii]=-10.;
   fTXiDCA[ii] =  -10000.;
   fTXiDaughtersDCA[ii] =  -10000.;
   fTXiMass[ii] =  -10000.;
   fTXiXiMass[ii] =  -10000.;
   fTXiOmegaMass[ii] =  -10000.;
   fTXiVr[ii] =  -10000.;
   fTXiPA[ii] =  -10000.;
   fTXiLambdaDCA[ii] =  -10000.;
   fTXiLambdaDaughtersDCA[ii] =  -10000.;
   fTXiLambdaMass[ii] =  -10000.;
   fTXiLambdaK0Mass[ii] =  -10000.;
   fTXiLambdaVr[ii] =  -10000.;
   fTXiLambdaPA[ii] =  -10000.;
   for(int jj=0;jj<3;jj++){
    fTXiTrackCharge[ii][jj]=-10.;
    fTXiTrackPx[ii][jj]=-100000.;
    fTXiTrackPy[ii][jj]=-100000.;
    fTXiTrackPz[ii][jj]=-100000.;
    fTXiTrackID[ii][jj]=-100000;
    fTXiTrackTPCmom[ii][jj]=-100000;
    fTXiTrackDCA[ii][jj]=-100000;
    fTXiTrackTPCsigma[ii][jj]=-100000;
    fTXiTrackTOFsigma[ii][jj]=-100000;
    fTXiTrackNcl[ii][jj]=-1;
    fTXiTrackCrR[ii][jj]=-100000;
    fTXiTrackCrF[ii][jj]=-100000;
    fTXiTrackSPDtime[ii][jj]=false;
    fTXiTrackITStime[ii][jj]=false;
    fTXiTrackTOFtime[ii][jj]=false;
   }
   fTXiMotherID[ii]=-1;
   fTXiPDG[ii]=-1;
   fTXiMotherPDG[ii]=-1;
   fTXiMotherWeak[ii]=-1;
   fTXiOrigin[ii]=-1;
  }
  fTnXi=0;
 
      static std::vector<AliFemtoDreamBasePart> Xis;
      Xis.clear();
      static std::vector<AliFemtoDreamBasePart> AntiXis;
      AntiXis.clear();


  AliAODEvent* aodEvt = dynamic_cast<AliAODEvent*>(fInputEvent);
  for (int iCasc = 0;iCasc< static_cast<TClonesArray *>(aodEvt->GetCascades())->GetEntriesFast();++iCasc) {

    Bool_t IsXi = kFALSE;
    Bool_t IsAntiXi = kFALSE;

    AliAODcascade* casc = aodEvt->GetCascade(iCasc);
    fCascade->SetCascade(fInputEvent, casc);


    if (fCutsXi->isSelected(fCascade)) {
      Xis.push_back(*fCascade);
      IsXi = kTRUE;
    }
    if (fCutsAntiXi->isSelected(fCascade)) {
      AntiXis.push_back(*fCascade);
      IsAntiXi = kTRUE;
    }

    if(IsXi||IsAntiXi) FillXi(fCascade,fisOmega);


  }//cascade loop


      //fill tree:
      //----------
      if(fTnKaon>0&&fTnXi>0) fTree->Fill(); // std kd 


     bool FemtoDreamPairing = false; // Skip FD pairing/mixing for now (to save computing time)
     if(fdoFDpairing) FemtoDreamPairing = true;
     if(FemtoDreamPairing){

      // Femto Dream Pair Cleaner. NOT SURE ABOUT THE histcounter !!!! not sure about the config !!!!
      // for now following the same approach as p-Omega
    
      fPairCleaner->ResetArray();

  fPairCleaner->CleanTrackAndDecay(&Kaons, &AntiXis, 0);
  fPairCleaner->CleanTrackAndDecay(&AntiKaons, &Xis, 1);

  fPairCleaner->CleanDecay(&Xis, 0);
  fPairCleaner->CleanDecay(&AntiXis, 1);


      fPairCleaner->StoreParticle(Kaons);
      fPairCleaner->StoreParticle(AntiKaons);
      fPairCleaner->StoreParticle(Xis);
      fPairCleaner->StoreParticle(AntiXis);

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
  PostData(4, fXiList);
  PostData(5, fAntiXiList);
  PostData(6, fResults);
  PostData(7, fResultsQA);
  PostData(8, fTree);
}

//________________________________________________________________________________________________
void AliAnalysisTaskOtonXx::ResetGlobalTrackReference() {
  for (UShort_t i = 0; i < fTrackBufferSize; i++) {
    fGTI[i] = 0;
  }
}
//________________________________________________________________________________________________
//AOD
/*
void AliAnalysisTaskOtonXx::StoreGlobalTrackReference(
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
void AliAnalysisTaskOtonXx::StoreGlobalTrackReference(AliVTrack *track) {
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
Bool_t AliAnalysisTaskOtonXx::FillKaon(AliFemtoDreamTrack *TheTrack) {
 Bool_t Filled = kFALSE;
 TVector3 mom;
 mom = TheTrack->GetMomentum();
 fTKaonPx[fTnKaon] = mom.X();
 fTKaonPy[fTnKaon] = mom.Y();
 fTKaonPz[fTnKaon] = mom.Z();
 fTKaonPTPC[fTnKaon] = TheTrack->GetMomTPC();
 fTKaonCharge[fTnKaon] = TheTrack->GetCharge().at(0);
 fTKaonNcl[fTnKaon] = TheTrack->GetNClsTPC();
 fTKaonPhi[fTnKaon] = (TheTrack->GetPhiAtRaidius().at(0)).at(0);//phi for r=85.cm ???
 fTKaonDCA[fTnKaon] = TheTrack->GetDCAXYProp(); //difference between DCAXY and DCAXYprop???
 fTKaonDCAz[fTnKaon] = TheTrack->GetDCAZProp(); //difference between DCAZ and DCAZprop???
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
 fTKaonSPDtime[fTnKaon] = TheTrack->GetHasSPDHit();
 fTKaonITStime[fTnKaon] = TheTrack->GetHasITSHit();
 fTKaonTOFtime[fTnKaon] = TheTrack->GetTOFTimingReuqirement();
  fTKaonFilterBit[fTnKaon] = TheTrack->GetFilterMap();
 fTKaonPDG[fTnKaon] = TheTrack->GetMCPDGCode();
 fTKaonMotherWeak[fTnKaon] = TheTrack->GetMotherWeak();




 AliFemtoDreamBasePart::PartOrigin org = TheTrack->GetParticleOrigin();
    fTKaonOrigin[fTnKaon] = -1;
    switch (org) {
      case AliFemtoDreamBasePart::kPhysPrimary:
        fTKaonOrigin[fTnKaon] = 0;
        break;
      case AliFemtoDreamBasePart::kWeak:
        fTKaonOrigin[fTnKaon] = 1;
        break;
      case AliFemtoDreamBasePart::kMaterial:
        fTKaonOrigin[fTnKaon] = 2;
        break;
      case AliFemtoDreamBasePart::kContamination:
        fTKaonOrigin[fTnKaon] = 3;
        break;
     }
 fTKaonMotherID[fTnKaon] = TheTrack->GetMotherID();


 fTnKaon++;
 Filled = kTRUE;
 return Filled;
}



//________________________________________________________________________________
Bool_t AliAnalysisTaskOtonXx::FillXi(AliFemtoDreamCascade *TheCasc, bool isomega) {
 Bool_t Filled = kFALSE;
 fTXiCharge[fTnXi] = TheCasc->GetCharge().at(0);

 fTXiDCA[fTnXi] = TheCasc->GetDCAXiPrimVtx();
 fTXiDaughtersDCA[fTnXi] = TheCasc->GetXiDCADaug();

 fTXiMass[fTnXi] = TheCasc->GetXiMass();
 if(isomega)  fTXiMass[fTnXi] = TheCasc->GetOmegaMass();

 fTXiXiMass[fTnXi] = TheCasc->GetXiMass();
 fTXiOmegaMass[fTnXi] = TheCasc->GetOmegaMass();
 fTXiVr[fTnXi] = TheCasc->GetXiTransverseRadius();
 fTXiPA[fTnXi] = TheCasc->GetCPA();
 fTXiLambdaDCA[fTnXi] = TheCasc->Getv0DCAPrimVtx();
 fTXiLambdaDaughtersDCA[fTnXi] = TheCasc->Getv0DCADaug();
 fTXiLambdaMass[fTnXi] = TheCasc->Getv0Mass();
 fTXiLambdaVr[fTnXi] = TheCasc->Getv0TransverseRadius();
 fTXiLambdaPA[fTnXi] = TheCasc->Getv0CPA();

 AliFemtoDreamTrack* TheTrack = TheCasc->GetPosDaug();
 for(int jj=0;jj<3;jj++){
  if(jj==0){ 

//select the proton
   if(fTXiCharge[fTnXi]==-1) {
    TheTrack = TheCasc->GetPosDaug();
    fTXiTrackDCA[fTnXi][jj]= TheCasc->Getv0PosToPrimVtx();
   }
   if(fTXiCharge[fTnXi]==1) {
    TheTrack = TheCasc->GetNegDaug();
    fTXiTrackDCA[fTnXi][jj]= TheCasc->Getv0NegToPrimVtx();
   }
   fTXiTrackTPCsigma[fTnXi][jj]=(TheTrack->GetnSigmaTPC((int) (AliPID::kProton)));
   fTXiTrackTOFsigma[fTnXi][jj]=(TheTrack->GetnSigmaTOF((int) (AliPID::kProton)));

//select the pion
  }else if(jj==1) { 
   if(fTXiCharge[fTnXi]==-1) {
    TheTrack = TheCasc->GetNegDaug();
    fTXiTrackDCA[fTnXi][jj]= TheCasc->Getv0NegToPrimVtx();
   }
   if(fTXiCharge[fTnXi]==1) {
    TheTrack = TheCasc->GetPosDaug();
    fTXiTrackDCA[fTnXi][jj]= TheCasc->Getv0PosToPrimVtx();
   }
   fTXiTrackTPCsigma[fTnXi][jj]=(TheTrack->GetnSigmaTPC((int) (AliPID::kPion)));
   fTXiTrackTOFsigma[fTnXi][jj]=(TheTrack->GetnSigmaTOF((int) (AliPID::kPion)));

//select the bachelor
  }else if(jj==2) {
   TheTrack = TheCasc->GetBach();
   fTXiTrackDCA[fTnXi][jj]= TheCasc->BachDCAPrimVtx();
   fTXiTrackTPCsigma[fTnXi][jj]=(TheTrack->GetnSigmaTPC((int) (AliPID::kKaon)));
   fTXiTrackTOFsigma[fTnXi][jj]=(TheTrack->GetnSigmaTOF((int) (AliPID::kKaon)));
  }

  TVector3 mom;
  mom = TheTrack->GetMomentum();
  fTXiTrackPx[fTnXi][jj] = mom.X();
  fTXiTrackPy[fTnXi][jj] = mom.Y();
  fTXiTrackPz[fTnXi][jj] = mom.Z();
  fTXiTrackCharge[fTnXi][jj] = TheTrack->GetCharge().at(0);
  fTXiTrackID[fTnXi][jj] = TheTrack->GetIDTracks().at(0);
  fTXiTrackNcl[fTnXi][jj] = TheTrack->GetNClsTPC();
  fTXiTrackCrF[fTnXi][jj] = TheTrack->GetRatioCr();
  fTXiTrackCrR[fTnXi][jj] = TheTrack->GetTPCCrossedRows();
  fTXiTrackTPCmom[fTnXi][jj] = TheTrack->GetMomTPC();
  fTXiTrackSPDtime[fTnXi][jj] = TheTrack->GetHasSPDHit();
  fTXiTrackITStime[fTnXi][jj] = TheTrack->GetHasITSHit();
  fTXiTrackTOFtime[fTnXi][jj] = TheTrack->GetTOFTimingReuqirement();
 }

 fTXiMotherID[fTnXi] = TheCasc->GetMotherID();

  fTXiPDG[fTnXi] = TheCasc->GetMCPDGCode();
  fTXiMotherPDG[fTnXi] = TheCasc->GetMotherPDG();
  fTXiMotherWeak[fTnXi] = TheCasc->GetMotherWeak();
  fTXiOrigin[fTnXi] =  TheCasc->GetParticleOrigin();

 fTnXi++;
 Filled = kTRUE;
 return Filled;
}


