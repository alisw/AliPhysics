/*
 * AliAnalysisTaskOtonOmegaNanoAOD.cxx
 *
 *  Created on: May 13, 2019
 *      Author: schmollweger
 */
#include "AliAnalysisTaskOtonOmegaNanoAOD.h"
#include "AliNanoAODTrack.h"
#include "TRandom3.h"

ClassImp(AliAnalysisTaskOtonOmegaNanoAOD)
AliAnalysisTaskOtonOmegaNanoAOD::AliAnalysisTaskOtonOmegaNanoAOD()
    : AliAnalysisTaskSE(),
      fisLightWeight(false),
      fOmegaTreeFlag(false),
      fEvent(nullptr),
      fEventCuts(nullptr),
      fEvtList(nullptr),
      fTrack(nullptr),
      fProton(nullptr),
      fProtonList(nullptr),
      fAntiProton(nullptr),
      fAntiProtonList(nullptr),
      fCascade(nullptr),
      fXi(nullptr),
      fXiList(nullptr),
      fAntiXi(nullptr),
      fAntiXiList(nullptr),
      fOmega(nullptr),
      fOmegaList(nullptr),
      fAntiOmega(nullptr),
      fAntiOmegaList(nullptr),
      fConfig(nullptr),
      fPairCleaner(nullptr),
      fPartColl(nullptr),
      fResults(nullptr),
      fResultsQA(nullptr),
      fTrackBufferSize(2000),
      fOmegaTree(0),
      fGTI(nullptr) {
}

AliAnalysisTaskOtonOmegaNanoAOD::AliAnalysisTaskOtonOmegaNanoAOD(const char* name, bool OmegaTreeFlag)
    : AliAnalysisTaskSE(name),
      fisLightWeight(false),
      fOmegaTreeFlag(false),
      fEvent(nullptr),
      fEventCuts(nullptr),
      fEvtList(nullptr),
      fTrack(nullptr),
      fProton(nullptr),
      fProtonList(nullptr),
      fAntiProton(nullptr),
      fAntiProtonList(nullptr),
      fCascade(nullptr),
      fXi(nullptr),
      fXiList(nullptr),
      fAntiXi(nullptr),
      fAntiXiList(nullptr),
      fOmega(nullptr),
      fOmegaList(nullptr),
      fAntiOmega(nullptr),
      fAntiOmegaList(nullptr),
      fConfig(nullptr),
      fPairCleaner(nullptr),
      fPartColl(nullptr),
      fResults(nullptr),
      fResultsQA(nullptr),
      fTrackBufferSize(2000),
      fOmegaTree(0),
      fGTI(nullptr) {
  DefineOutput(1, TList::Class());  //Output for the Event Cuts
  DefineOutput(2, TList::Class());  //Output for the Proton Cuts
  DefineOutput(3, TList::Class());  //Output for the AntiProton Cuts
  DefineOutput(4, TList::Class());  //Output for the Xi Cuts
  DefineOutput(5, TList::Class());  //Output for the AntiXi Cuts
  DefineOutput(6, TList::Class());  //Output for the Omega Cuts
  DefineOutput(7, TList::Class());  //Output for the AntiOmega Cuts
  DefineOutput(8, TList::Class());  //Output for the Results
  DefineOutput(9, TList::Class());  //Output for the Results QA
  DefineOutput(10, TTree::Class());  // OmegaTree
}

AliAnalysisTaskOtonOmegaNanoAOD::~AliAnalysisTaskOtonOmegaNanoAOD() {
  if (fEvent) {
    delete fEvent;
  }
  if (fEventCuts) {
    delete fEventCuts;
  }
  if (fTrack) {
    delete fTrack;
  }
  if (fProton) {
    delete fProton;
  }
  if (fAntiProton) {
    delete fAntiProton;
  }
  if (fCascade) {
    delete fCascade;
  }
  if (fXi) {
    delete fXi;
  }
  if (fAntiXi) {
    delete fAntiXi;
  }
  if (fOmega) {
    delete fOmega;
  }
  if (fAntiOmega) {
    delete fAntiOmega;
  }
  if (fPairCleaner) {
    delete fPairCleaner;
  }
  if (fPartColl) {
    delete fPartColl;
  }
  if (fOmegaTree) {
    delete fOmegaTree;//not sure about this
  }
}

void AliAnalysisTaskOtonOmegaNanoAOD::UserCreateOutputObjects() {
  fGTI = new AliVTrack *[fTrackBufferSize];

  if (!fEventCuts) {
    AliError("No Event cuts \n");
  } else {
    fEventCuts->InitQA();
  }
  if (!fProton) {
    AliError("No Proton cuts \n");
  } else {
    fProton->Init();
  }
  if (!fAntiProton) {
    AliError("No AntiProton cuts \n");
  } else {
    fAntiProton->Init();
  }
  if (!fXi) {
    AliError("No Xi cuts \n");
  } else {
    fXi->Init();
  }
  if (!fAntiXi) {
    AliError("No AntiXi cuts \n");
  } else {
    fAntiXi->Init();
  }
  if (!fOmega) {
    AliError("No Omega cuts \n");
  } else {
    fOmega->Init();
  }
  if (!fAntiOmega) {
    AliError("No AntiOmega cuts \n");
  } else {
    fAntiOmega->Init();
  }
  if (!fConfig) {
    AliError("No Correlation Config \n");
  } else {
    fPartColl = new AliFemtoDreamPartCollection(fConfig,
                                                fConfig->GetMinimalBookingME());
    fPairCleaner = new AliFemtoDreamPairCleaner(2, 2,
                                                fConfig->GetMinimalBookingME());
  }
  fEvent = new AliFemtoDreamEvent(true, !fisLightWeight,
                                  GetCollisionCandidates(), false);
  fEvent->SetMultiplicityEstimator(fConfig->GetMultiplicityEstimator());

  fTrack = new AliFemtoDreamTrack();
  fTrack->SetUseMCInfo(false);

  fCascade = new AliFemtoDreamCascade();
  fCascade->SetUseMCInfo(false);
  //PDG Codes should be set assuming Xi- to also work for Xi+
  fCascade->SetPDGCode(3334);
  fCascade->SetPDGDaugPos(2212);
  fCascade->GetPosDaug()->SetUseMCInfo(false);
  fCascade->SetPDGDaugNeg(211);
  fCascade->GetNegDaug()->SetUseMCInfo(false);
  fCascade->SetPDGDaugBach(321);
  fCascade->GetBach()->SetUseMCInfo(false);
  fCascade->Setv0PDGCode(3122);


  if (!fEventCuts->GetMinimalBooking()) {
    fEvtList = fEventCuts->GetHistList();
  } else {
    fEvtList = new TList();
    fEvtList->SetName("EventCuts");
    fEvtList->SetOwner();
  }
  if (!fProton->GetMinimalBooking()) {
    fProtonList = fProton->GetQAHists();
  } else {
    fProtonList = new TList();
    fProtonList->SetName("TrackCuts");
    fProtonList->SetOwner();
  }
  if (!fAntiProton->GetMinimalBooking()) {
    fAntiProtonList = fAntiProton->GetQAHists();
  } else {
    fAntiProtonList = new TList();
    fAntiProtonList->SetName("AntiTrackCuts");
    fAntiProtonList->SetOwner();
  }
  if (!fXi->GetMinimalBooking()) {
    fXiList = fXi->GetQAHists();
  } else {
    fXiList = new TList();
    fXiList->SetName("XiCuts");
    fXiList->SetOwner();
  }
  if (!fAntiXi->GetMinimalBooking()) {
    fAntiXiList = fAntiXi->GetQAHists();
  } else {
    fAntiXiList = new TList();
    fAntiXiList->SetName("AntiXiCuts");
    fAntiXiList->SetOwner();
  }

  if (!fOmega->GetMinimalBooking()) {
    fOmegaList = fOmega->GetQAHists();
  } else {
    fOmegaList = new TList();
    fOmegaList->SetName("OmegaCuts");
    fOmegaList->SetOwner();
  }
  if (!fAntiOmega->GetMinimalBooking()) {
    fAntiOmegaList = fAntiOmega->GetQAHists();
  } else {
    fAntiOmegaList = new TList();
    fAntiOmegaList->SetName("AntiOmegaCuts");
    fAntiOmegaList->SetOwner();
  }

  fResultsQA = new TList();
  fResultsQA->SetOwner();
  fResultsQA->SetName("ResultsQA");
  if (!fConfig->GetMinimalBookingME()) {
    fResults = fPartColl->GetHistList();
    fResultsQA->Add(fPartColl->GetQAList());
    fResultsQA->Add(fPairCleaner->GetHistList());
  } else {
    fResults = new TList();
    fResults->SetOwner();
    fResults->SetName("Results");
  }

  //book tree:
  ////////////
  fOmegaTree = new TTree("oTTree","a simple TTree");
  fOmegaTree->SetDirectory(0); // This is to force a memory-resident Tree, and avoid errors. // ????? is this necessary? does it create memory problems?
  fOmegaTree->Branch("RunNumber",&fTRunNumber,"fTRunNumber/I");
  fOmegaTree->Branch("Vz",&fTVz,"fTVz/F");
  fOmegaTree->Branch("Mult",&fTMult,"fTMult/I");
  //protons:
  fOmegaTree->Branch("nProton",&fTnProton,"fTnProton/I");
  fOmegaTree->Branch("ProtonPx",&fTProtonPx,"fTProtonPx[fTnProton]/F");
  fOmegaTree->Branch("ProtonPy",&fTProtonPy,"fTProtonPy[fTnProton]/F");
  fOmegaTree->Branch("ProtonPz",&fTProtonPz,"fTProtonPz[fTnProton]/F");
  fOmegaTree->Branch("ProtonEta",&fTProtonEta,"fTProtonEta[fTnProton]/F");
  fOmegaTree->Branch("ProtonCharge",&fTProtonCharge,"fTProtonCharge[fTnProton]/S");
  fOmegaTree->Branch("ProtonTPCsigma",&fTProtonTPCsigma,"fTProtonTPCsigma[fTnProton]/F");
  fOmegaTree->Branch("ProtonTOFsigma",&fTProtonTOFsigma,"fTProtonTOFsigma[fTnProton]/F");
  fOmegaTree->Branch("ProtonNcl",&fTProtonNcl,"fTProtonNcl[fTnProton]/I");
  fOmegaTree->Branch("ProtonPhi",&fTProtonPhi,"fTProtonPhi[fTnProton]/F");
  fOmegaTree->Branch("ProtonDCA",&fTProtonDCA,"fTProtonDCA[fTnProton]/F");
  fOmegaTree->Branch("ProtonID",&fTProtonID,"fTProtonID[fTnProton]/I");
 //omegas:
 fOmegaTree->Branch("nCascade",&fTnCascade,"fTnCascade/I");
 fOmegaTree->Branch("CascadeCharge",&fTCascadeCharge,"fTCascadeCharge[fTnCascade]/S");
 fOmegaTree->Branch("CascadeDCA",&fTCascadeDCA,"fTCascadeDCA[fTnCascade]/F");
 fOmegaTree->Branch("CascadeDaughtersDCA",&fTCascadeDaughtersDCA,"fTCascadeDaughtersDCA[fTnCascade]/F");
 fOmegaTree->Branch("CascadeXiMass",&fTCascadeXiMass,"fTCascadeXiMass[fTnCascade]/F");
 fOmegaTree->Branch("CascadeOmegaMass",&fTCascadeOmegaMass,"fTCascadeOmegaMass[fTnCascade]/F");
 fOmegaTree->Branch("CascadeVr",&fTCascadeVr,"fTCascadeVr[fTnCascade]/F");
 fOmegaTree->Branch("CascadePA",&fTCascadePA,"fTCascadePA[fTnCascade]/F");
 fOmegaTree->Branch("LambdaDCA",&fTLambdaDCA,"fTLambdaDCA[fTnCascade]/F");
 fOmegaTree->Branch("LambdaDaughtersDCA",&fTLambdaDaughtersDCA,"fTLambdaDaughtersDCA[fTnCascade]/F");
 fOmegaTree->Branch("LambdaMass",&fTLambdaMass,"fTLambdaMass[fTnCascade]/F");
 fOmegaTree->Branch("LambdaVr",&fTLambdaVr,"fTLambdaVr[fTnCascade]/F");
 fOmegaTree->Branch("LambdaPA",&fTLambdaPA,"fTLambdaPA[fTnCascade]/F");
 fOmegaTree->Branch("TrackPx",&fTTrackPx,"fTTrackPx[fTnCascade][3]/F");
 fOmegaTree->Branch("TrackPy",&fTTrackPy,"fTTrackPy[fTnCascade][3]/F");
 fOmegaTree->Branch("TrackPz",&fTTrackPz,"fTTrackPz[fTnCascade][3]/F");
 fOmegaTree->Branch("TrackEta",&fTTrackEta,"fTTrackEta[fTnCascade][3]/F");
 fOmegaTree->Branch("TrackCharge",&fTTrackCharge,"fTTrackCharge[fTnCascade][3]/S");
 fOmegaTree->Branch("TrackDCA",&fTTrackDCA,"fTTrackDCA[fTnCascade][3]/F");
 fOmegaTree->Branch("TrackTPCsigma",&fTTrackTPCsigma,"fTTrackTPCsigma[fTnCascade][3]/F");
 fOmegaTree->Branch("TrackTOFsigma",&fTTrackTOFsigma,"fTTrackTOFsigma[fTnCascade][3]/F");
 fOmegaTree->Branch("TrackNcl",&fTTrackNcl,"fTTrackNcl[fTnCascade][3]/I");
 fOmegaTree->Branch("TrackCrR",&fTTrackCrR,"fTTrackCrR[fTnCascade][3]/F");
 fOmegaTree->Branch("TrackCrF",&fTTrackCrF,"fTTrackCrF[fTnCascade][3]/F");
 fOmegaTree->Branch("TrackITStime",&fTTrackITStime,"fTTrackITStime[fTnCascade][3]/O");
 fOmegaTree->Branch("TrackTOFtime",&fTTrackTOFtime,"fTTrackTOFtime[fTnCascade][3]/O");
 fOmegaTree->Branch("TrackFilterBit",&fTTrackFilterBit,"fTTrackFilterBit[fTnCascade][3]/O");
 fOmegaTree->Branch("TrackPhi",&fTTrackPhi,"fTTrackPhi[fTnCascade][3]/F");
 fOmegaTree->Branch("TrackID",&fTTrackID,"fTTrackID[fTnCascade][3]/I");

  PostData(1, fEvtList);
  PostData(2, fProtonList);
  PostData(3, fAntiProtonList);
  PostData(4, fXiList);
  PostData(5, fAntiXiList);
  PostData(6, fOmegaList);
  PostData(7, fAntiOmegaList);
  PostData(8, fResults);
  PostData(9, fResultsQA);
  PostData(10, fOmegaTree);
}

void AliAnalysisTaskOtonOmegaNanoAOD::UserExec(Option_t *option) {
//  AliVEvent *fInputEvent = InputEvent();
  if (!fInputEvent) {
    AliError("No input event");
    return;
  }
  fEvent->SetEvent(fInputEvent);
  if (!fEventCuts->isSelected(fEvent)) {
    return;
  }

  //fTRunNumber = fInputEvent->GetRunNumber(); // Old method for ESD/AOD (not working for NanoAOD)

  //// For NanoAOD LHC16/17/18 starting with filtering train #100:
  //AliNanoAODHeader* nanoHeader = dynamic_cast<AliNanoAODHeader*>(fInputEvent->GetHeader());
  //fTRunNumber = nanoHeader->GetVarInt(nanoHeader->GetRunNumberIndex());

  fTRunNumber = 0.; //For NanoAOD filtering trains <100, no info

  Double_t PrimVtx[3];
  fInputEvent->GetPrimaryVertex()->GetXYZ(PrimVtx);
  fTVz = PrimVtx[2];
  fTMult = fEvent->GetMultiplicity();

//init tree
  for(int ii=0;ii<MAXPROTONS;ii++){
   fTProtonP[ii]=-100000.;
   fTProtonPt[ii]=-100000.;
   fTProtonmT[ii]=-100000.;
   fTProtonEta[ii]=-100000.;
   fTProtonPx[ii]=-100000.;
   fTProtonPy[ii]=-100000.;
   fTProtonPz[ii]=-100000.;
   fTProtonVPx[ii]=-100000.;
   fTProtonVPy[ii]=-100000.;
   fTProtonVPz[ii]=-100000.;
   fTProtonTPCmom[ii]=-100000.;
   fTProtonCharge[ii]=-10;
   fTProtonDCA[ii]=-100000.;
   fTProtonTPCsigma[ii]=-100000.;
   fTProtonTOFsigma[ii]=-100000.;
   fTProtonNcl[ii]=-100000;
   fTProtonShared[ii]=-100000;
   fTProtonTPCchi2[ii]=-100000.;
   fTProtonITStime[ii]=kFALSE;
   fTProtonTOFtime[ii]=kFALSE;
   fTProtonTPConly[ii]=kFALSE;
   fTProtonITScomplementary[ii]=kFALSE;
   fTProtonITSpure[ii]=kFALSE;
   fTProtonGLOBAL[ii]=kFALSE;
   fTProtonPhi[ii]=-100000.;
   fTProtonID[ii]=-100000;
  }
  fTnProton=0;

  // PROTON SELECTION  (proton loop)
  ResetGlobalTrackReference();
  for (int iTrack = 0; iTrack < fInputEvent->GetNumberOfTracks(); ++iTrack) {
    AliVTrack *track = static_cast<AliVTrack *>(fInputEvent->GetTrack(iTrack));
    if (!track) {
      AliFatal("No Standard AOD");
      return;
    }
    StoreGlobalTrackReference(track);
  }
  std::vector<AliFemtoDreamBasePart> Protons;
  std::vector<AliFemtoDreamBasePart> AntiProtons;
  const int multiplicity = fEvent->GetMultiplicity();
  fTrack->SetGlobalTrackInfo(fGTI, fTrackBufferSize);
  for (int iTrack = 0; iTrack < fInputEvent->GetNumberOfTracks(); ++iTrack) {
    Bool_t IsProton = kFALSE;
    Bool_t IsAntiProton = kFALSE;
    AliVTrack *track = static_cast<AliVTrack *>(fInputEvent->GetTrack(iTrack));
    fTrack->SetTrack(track, fInputEvent, multiplicity);
    if (fProton->isSelected(fTrack)) {
      Protons.push_back(*fTrack);
      IsProton = kTRUE;
    }
    if (fAntiProton->isSelected(fTrack)) {
      AntiProtons.push_back(*fTrack);
      IsAntiProton = kTRUE;
    }

    if(IsProton||IsAntiProton) FillProton(fTrack);
  }

  //init tree
  for(int ii=0;ii<MAXCASCADES;ii++){
   fTCascadePx[ii]=-100000.;
   fTCascadePy[ii]=-100000.;
   fTCascadePz[ii]=-100000.;
   fTCascadeP[ii]=-100000.;
   fTCascadePt[ii]=-100000.;
   fTCascademT[ii]=-100000.;
   fTCascadeCharge[ii]=-10;
   fTCascadeDCA[ii]=-100000.;
   fTCascadeDaughtersDCA[ii]=-100000.;
   fTCascadeXiMass[ii]=-100000.;
   fTCascadeOmegaMass[ii]=-100000.;
   fTCascadeVr[ii]=-100000.;
   fTCascadeVx[ii]=-100000.;
   fTCascadeVy[ii]=-100000.;
   fTCascadeVz[ii]=-100000.;
   fTCascadePA[ii]=-100000.;
   fTLambdaPt[ii]=-100000.;
   fTLambdaPx[ii]=-100000.;
   fTLambdaPy[ii]=-100000.;
   fTLambdaPz[ii]=-100000.;
   fTLambdaDCA[ii]=-100000.;
   fTLambdaDaughtersDCA[ii]=-100000.;
   fTLambdaMass[ii]=-100000.;
   fTLambdaK0Mass[ii]=-100000.;
   fTLambdaVr[ii]=-100000.;
   fTLambdaVx[ii]=-100000.;
   fTLambdaVy[ii]=-100000.;
   fTLambdaVz[ii]=-100000.;
   fTLambdaPA[ii]=-100000.;
   for(int jj=0;jj<3;jj++){
    fTTrackP[ii][jj]=-100000.;
    fTTrackPx[ii][jj]=-100000.;
    fTTrackPy[ii][jj]=-100000.;
    fTTrackPz[ii][jj]=-100000.;
    fTTrackTPCmom[ii][jj]=-100000.;
    fTTrackEta[ii][jj]=-100000.;
    fTTrackCharge[ii][jj]=-10;
    fTTrackDCA[ii][jj]=-100000.;
    fTTrackTPCsigma[ii][jj]=-100000.;
    fTTrackTOFsigma[ii][jj]=-100000.;
    fTTrackNcl[ii][jj]=-100000;
    fTTrackCrR[ii][jj]=-100000.;
    fTTrackCrF[ii][jj]=-100000.;
    fTTrackShared[ii][jj]=-100000;
    fTTrackTPCchi2[ii][jj]=-100000.;
    fTTrackITStime[ii][jj]=kFALSE;
    fTTrackTOFtime[ii][jj]=kFALSE;
    fTTrackTPConly[ii][jj]=kFALSE;
    fTTrackITScomplementary[ii][jj]=kFALSE;
    fTTrackITSpure[ii][jj]=kFALSE;
    fTTrackGLOBAL[ii][jj]=kFALSE;
    fTTrackFilterBit[ii][jj]=0;
    fTTrackPhi[ii][jj]=-100000.;
    fTTrackID[ii][jj]=-100000;
   }
  }
  fTnCascade=0;

  //Xi (bkg) and omegas (cascade loop)
  std::vector<AliFemtoDreamBasePart> Xis;
  std::vector<AliFemtoDreamBasePart> AntiXis;
  std::vector<AliFemtoDreamBasePart> Omegas;
  std::vector<AliFemtoDreamBasePart> AntiOmegas;
  AliAODEvent* aodEvt = dynamic_cast<AliAODEvent*>(fInputEvent);
  for (int iCasc = 0;iCasc< static_cast<TClonesArray *>(aodEvt->GetCascades())->GetEntriesFast();++iCasc) {
    Bool_t IsOmegaBkg = kFALSE;
    Bool_t IsAntiOmegaBkg = kFALSE;
    Bool_t IsOmega = kFALSE;
    Bool_t IsAntiOmega = kFALSE;

    AliAODcascade* casc = aodEvt->GetCascade(iCasc);
    fCascade->SetCascade(fInputEvent, casc);

    if (fXi->isSelected(fCascade)) {
      Xis.push_back(*fCascade);
      IsOmegaBkg = kTRUE;
    }
    if (fAntiXi->isSelected(fCascade)) {
      AntiXis.push_back(*fCascade);
      IsAntiOmegaBkg = kTRUE;
    }

    if (fOmega->isSelected(fCascade)) {
      Omegas.push_back(*fCascade);
      IsOmega = kTRUE;
    }
    if (fAntiOmega->isSelected(fCascade)) {
      AntiOmegas.push_back(*fCascade);
      IsAntiOmega = kTRUE;
    }

    if(IsOmegaBkg||IsAntiOmegaBkg||IsOmega||IsAntiOmega) FillCascade(fCascade);

  }

  //fill Tree
  if(fTnProton>0&&fTnCascade>0) fOmegaTree->Fill(); //Fill when at least 1 proton AND 1 cascade
  //if(fTnProton>0||fTnCascade>0) fOmegaTree->Fill(); //Fill when at least 1 proton OR 1 cascade

  //pair cleaner
  fPairCleaner->ResetArray();

  fPairCleaner->CleanTrackAndDecay(&Protons, &Xis, 0);
  fPairCleaner->CleanTrackAndDecay(&AntiProtons, &AntiXis, 1);
  fPairCleaner->CleanTrackAndDecay(&Protons, &Omegas, 0); //this is apparently working
  fPairCleaner->CleanTrackAndDecay(&AntiProtons, &AntiOmegas, 1); //this is apparently working

  fPairCleaner->CleanDecay(&Xis, 0);
  fPairCleaner->CleanDecay(&AntiXis, 1);
  fPairCleaner->CleanDecay(&Omegas, 0); 
  fPairCleaner->CleanDecay(&AntiOmegas, 1); 

  fPairCleaner->StoreParticle(Protons);
  fPairCleaner->StoreParticle(AntiProtons);
  fPairCleaner->StoreParticle(Xis);
  fPairCleaner->StoreParticle(AntiXis);
  fPairCleaner->StoreParticle(Omegas);
  fPairCleaner->StoreParticle(AntiOmegas);


  fPartColl->SetEvent(fPairCleaner->GetCleanParticles(), fEvent->GetZVertex(),
                      fEvent->GetMultiplicity(), fEvent->GetV0MCentrality());
  PostData(1, fEvtList);
  PostData(2, fProtonList);
  PostData(3, fAntiProtonList);
  PostData(4, fXiList);
  PostData(5, fAntiXiList);
  PostData(6, fOmegaList);
  PostData(7, fAntiOmegaList);
  PostData(8, fResults);
  PostData(9, fResultsQA);
  PostData(10, fOmegaTree);
}

//____________________________________________________________________________________________________
void AliAnalysisTaskOtonOmegaNanoAOD::ResetGlobalTrackReference() {
  // see AliFemtoDreamAnalysis for details
  for (int i = 0; i < fTrackBufferSize; i++) {
    fGTI[i] = 0;
  }
}

//____________________________________________________________________________________________________
void AliAnalysisTaskOtonOmegaNanoAOD::StoreGlobalTrackReference(AliVTrack *track) {
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
Bool_t AliAnalysisTaskOtonOmegaNanoAOD::FillCascade(AliFemtoDreamCascade *TheCasc) {
 Bool_t Filled = kFALSE;

 fTCascadeCharge[fTnCascade] = TheCasc->GetCharge().at(0);
 fTCascadeDCA[fTnCascade] = TheCasc->GetDCAXiPrimVtx();
 fTCascadeDaughtersDCA[fTnCascade] = TheCasc->GetXiDCADaug();
 fTCascadeXiMass[fTnCascade] = TheCasc->GetXiMass();
 fTCascadeOmegaMass[fTnCascade] = TheCasc->GetOmegaMass();
 fTCascadeVr[fTnCascade] = TheCasc->GetXiTransverseRadius();
 fTCascadePA[fTnCascade] = TheCasc->GetCPA();

 fTLambdaDCA[fTnCascade] = TheCasc->Getv0DCAPrimVtx();
 fTLambdaDaughtersDCA[fTnCascade] = TheCasc->Getv0DCADaug();
 fTLambdaMass[fTnCascade] = TheCasc->Getv0Mass();
 fTLambdaVr[fTnCascade] = TheCasc->Getv0TransverseRadius();
 fTLambdaPA[fTnCascade] = TheCasc->Getv0CPA();

 AliFemtoDreamTrack* TheTrack = TheCasc->GetPosDaug();
 for(int jj=0;jj<3;jj++){
  if(jj==0){ //select the proton
   if(fTCascadeCharge[fTnCascade]==-1) {TheTrack = TheCasc->GetPosDaug();fTTrackDCA[fTnCascade][jj]= TheCasc->Getv0PosToPrimVtx();}
   if(fTCascadeCharge[fTnCascade]==1) {TheTrack = TheCasc->GetNegDaug();fTTrackDCA[fTnCascade][jj]= TheCasc->Getv0NegToPrimVtx();}
   fTTrackTPCsigma[fTnCascade][jj]=(TheTrack->GetnSigmaTPC((int) (AliPID::kProton)));
   fTTrackTOFsigma[fTnCascade][jj]=(TheTrack->GetnSigmaTOF((int) (AliPID::kProton)));
  }else if(jj==1) { //select the pion
   if(fTCascadeCharge[fTnCascade]==-1) {TheTrack = TheCasc->GetNegDaug();fTTrackDCA[fTnCascade][jj]= TheCasc->Getv0NegToPrimVtx();}
   if(fTCascadeCharge[fTnCascade]==1) {TheTrack = TheCasc->GetPosDaug();fTTrackDCA[fTnCascade][jj]= TheCasc->Getv0PosToPrimVtx();}
   fTTrackTPCsigma[fTnCascade][jj]=(TheTrack->GetnSigmaTPC((int) (AliPID::kPion)));
   fTTrackTOFsigma[fTnCascade][jj]=(TheTrack->GetnSigmaTOF((int) (AliPID::kPion)));
  }else if(jj==2) { //select the bachelor
   TheTrack = TheCasc->GetBach();
   fTTrackDCA[fTnCascade][jj]= TheCasc->BachDCAPrimVtx();
   fTTrackTPCsigma[fTnCascade][jj]=(TheTrack->GetnSigmaTPC((int) (AliPID::kKaon)));
   fTTrackTOFsigma[fTnCascade][jj]=(TheTrack->GetnSigmaTOF((int) (AliPID::kKaon)));
  }

  TVector3 mom;
  mom = TheTrack->GetMomentum();
  fTTrackPx[fTnCascade][jj] = mom.X();
  fTTrackPy[fTnCascade][jj] = mom.Y();
  fTTrackPz[fTnCascade][jj] = mom.Z();
  fTTrackEta[fTnCascade][jj] = TheTrack->GetEta().at(0);
  fTTrackCharge[fTnCascade][jj] = TheTrack->GetCharge().at(0);
  fTTrackNcl[fTnCascade][jj] = TheTrack->GetNClsTPC();
  fTTrackCrF[fTnCascade][jj] = TheTrack->GetRatioCr();
  fTTrackCrR[fTnCascade][jj] = TheTrack->GetTPCCrossedRows();
  fTTrackITStime[fTnCascade][jj] = TheTrack->GetHasITSHit();
  fTTrackTOFtime[fTnCascade][jj] = TheTrack->GetTOFTimingReuqirement();
  fTTrackFilterBit[fTnCascade][jj] = TheTrack->GetFilterMap();
  fTTrackPhi[fTnCascade][jj] = (TheTrack->GetPhiAtRaidius().at(0)).at(0);//phi for r=85.cm ???
  fTTrackID[fTnCascade][jj] = TheTrack->GetIDTracks().at(0);
 }

 fTnCascade++;
 Filled = kTRUE;
 return Filled;
} 


//________________________________________________________________________________
Bool_t AliAnalysisTaskOtonOmegaNanoAOD::FillProton(AliFemtoDreamTrack *TheTrack) {
 Bool_t Filled = kFALSE;

 TVector3 mom;
 mom = TheTrack->GetMomentum();
 fTProtonPx[fTnProton] = mom.X();
 fTProtonPy[fTnProton] = mom.Y();
 fTProtonPz[fTnProton] = mom.Z();
 fTProtonEta[fTnProton] = TheTrack->GetEta().at(0);
 fTProtonCharge[fTnProton] = TheTrack->GetCharge().at(0);
 fTProtonTPCsigma[fTnProton] = (TheTrack->GetnSigmaTPC((int) (AliPID::kProton)));
 fTProtonTOFsigma[fTnProton] = (TheTrack->GetnSigmaTOF((int) (AliPID::kProton)));
 fTProtonNcl[fTnProton] = TheTrack->GetNClsTPC();
 fTProtonPhi[fTnProton] = (TheTrack->GetPhiAtRaidius().at(0)).at(0);//phi for r=85.cm ???
 fTProtonDCA[fTnProton] = TheTrack->GetDCAXYProp();
 fTProtonID[fTnProton] = TheTrack->GetIDTracks().at(0);

 fTnProton++;

 Filled = kTRUE;
 return Filled;
} 


