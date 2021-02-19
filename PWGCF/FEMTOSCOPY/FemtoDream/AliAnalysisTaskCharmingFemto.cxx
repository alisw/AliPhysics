#include "AliAnalysisTaskCharmingFemto.h"

#include "yaml-cpp/yaml.h"

#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliMultSelection.h"
#include "AliPIDResponse.h"
#include "AliRDHFCutsDplustoKpipi.h"
#include "AliAODRecoDecayHF.h"
#include "AliAODRecoDecayHF3Prong.h"
#include "AliHFMLResponse.h"
#include "AliAODHandler.h"
#include "AliHFMLResponseDplustoKpipi.h"
#include "AliAnalysisTaskSECharmHadronMLSelector.h"
#include "AliMLModelHandler.h"
#include "TDatabasePDG.h"
#include "TRandom.h"

ClassImp(AliAnalysisTaskCharmingFemto)

    //____________________________________________________________________________________________________
    AliAnalysisTaskCharmingFemto::AliAnalysisTaskCharmingFemto()
    : AliAnalysisTaskSE("AliAnalysisTaskCharmingFemto"),
      fInputEvent(nullptr),
      fEvent(nullptr),
      fEvtCuts(nullptr),
      fProtonTrack(nullptr),
      fTrackCutsPartProton(nullptr),
      fTrackCutsPartAntiProton(nullptr),
      fConfig(nullptr),
      fPairCleaner(nullptr),
      fPartColl(nullptr),
      fIsMC(false),
      fIsLightweight(false),
      fTrigger(AliVEvent::kINT7),
      fSystem(kpp13TeV),
      fTrackBufferSize(2500),
      fDmesonPDGs{},
      fGTI(nullptr),
      fQA(nullptr),
      fEvtHistList(nullptr),
      fTrackCutHistList(nullptr),
      fTrackCutHistMCList(nullptr),
      fAntiTrackCutHistList(nullptr),
      fAntiTrackCutHistMCList(nullptr),
      fDChargedHistList(nullptr),
      fResultList(nullptr),
      fResultQAList(nullptr),
      fHistDplusInvMassPt(nullptr),
      fHistDplusInvMassPtSel(nullptr),
      fHistDplusEta(nullptr),
      fHistDplusPhi(nullptr),
      fHistDplusChildPt(),
      fHistDplusChildEta(),
      fHistDplusChildPhi(),
      fHistDplusMCPDGPt(nullptr),
      fHistDplusMCPtRes(nullptr),
      fHistDplusMCPhiRes(nullptr),
      fHistDplusMCThetaRes(nullptr),
      fHistDplusMCOrigin(nullptr),
      fHistDminusInvMassPt(nullptr),
      fHistDminusInvMassPtSel(nullptr),
      fHistDminusEta(nullptr),
      fHistDminusPhi(nullptr),
      fHistDminusChildPt(),
      fHistDminusChildEta(),
      fHistDminusChildPhi(),
      fHistDminusMCPDGPt(nullptr),
      fHistDminusMCPtRes(nullptr),
      fHistDminusMCPhiRes(nullptr),
      fHistDminusMCThetaRes(nullptr),
      fHistDminusMCOrigin(nullptr),
      fDecChannel(kDplustoKpipi),
      fRDHFCuts(nullptr),
      fAODProtection(0),
      fMassSelectionType(kSignal),
      fNSigmaMass(2.),
      fNSigmaOffsetSideband(5.),
      fLowerMassSelection(0.),
      fUpperMassSelection(999.),
      fSidebandWidth(0.2),
      fLowerDstarRemoval(1.992),
      fUpperDstarRemoval(2.028),
      fMCBeautyRejection(false),
      fMCBeautyScalingFactor(1.),
      fUseTrueDOnly(false),
      fApplyML(false),
      fConfigPath(""),
      fMLResponse(nullptr),
      fDependOnMLSelector(false),
      fPtLimsML{},
      fMLScoreCuts{},
      fMLOptScoreCuts{} {}

//____________________________________________________________________________________________________
AliAnalysisTaskCharmingFemto::AliAnalysisTaskCharmingFemto(const char *name,
                                                           const bool isMC)
    : AliAnalysisTaskSE(name),
      fInputEvent(nullptr),
      fEvent(nullptr),
      fEvtCuts(nullptr),
      fProtonTrack(nullptr),
      fTrackCutsPartProton(nullptr),
      fTrackCutsPartAntiProton(nullptr),
      fConfig(nullptr),
      fPairCleaner(nullptr),
      fPartColl(nullptr),
      fIsMC(isMC),
      fIsLightweight(false),
      fTrigger(AliVEvent::kINT7),
      fSystem(kpp13TeV),
      fTrackBufferSize(2500),
      fDmesonPDGs{},
      fGTI(nullptr),
      fQA(nullptr),
      fEvtHistList(nullptr),
      fTrackCutHistList(nullptr),
      fTrackCutHistMCList(nullptr),
      fAntiTrackCutHistList(nullptr),
      fAntiTrackCutHistMCList(nullptr),
      fDChargedHistList(nullptr),
      fResultList(nullptr),
      fResultQAList(nullptr),
      fHistDplusInvMassPt(nullptr),
      fHistDplusInvMassPtSel(nullptr),
      fHistDplusEta(nullptr),
      fHistDplusPhi(nullptr),
      fHistDplusChildPt(),
      fHistDplusChildEta(),
      fHistDplusChildPhi(),
      fHistDplusMCPDGPt(nullptr),
      fHistDminusInvMassPt(nullptr),
      fHistDminusInvMassPtSel(nullptr),
      fHistDplusMCPtRes(nullptr),
      fHistDplusMCPhiRes(nullptr),
      fHistDplusMCThetaRes(nullptr),
      fHistDplusMCOrigin(nullptr),
      fHistDminusEta(nullptr),
      fHistDminusPhi(nullptr),
      fHistDminusChildPt(),
      fHistDminusChildEta(),
      fHistDminusChildPhi(),
      fHistDminusMCPDGPt(nullptr),
      fHistDminusMCPtRes(nullptr),
      fHistDminusMCPhiRes(nullptr),
      fHistDminusMCThetaRes(nullptr),
      fHistDminusMCOrigin(nullptr),
      fDecChannel(kDplustoKpipi),
      fRDHFCuts(nullptr),
      fAODProtection(0),
      fMassSelectionType(kSignal),
      fNSigmaMass(2.),
      fNSigmaOffsetSideband(5.),
      fLowerMassSelection(0.),
      fUpperMassSelection(999.),
      fSidebandWidth(0.2),
      fLowerDstarRemoval(1.992),
      fUpperDstarRemoval(2.028),
      fMCBeautyRejection(false),
      fMCBeautyScalingFactor(1.),
      fUseTrueDOnly(false),
      fApplyML(false),
      fConfigPath(""),
      fMLResponse(nullptr),
      fDependOnMLSelector(false),
      fPtLimsML{},
      fMLScoreCuts{},
      fMLOptScoreCuts{} {
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
  DefineOutput(3, TList::Class());
  DefineOutput(4, TList::Class());
  DefineOutput(5, TList::Class());
  DefineOutput(6, TList::Class());
  DefineOutput(7, TList::Class());
  switch(fDecChannel){ // save cut object for HF particle
    case kDplustoKpipi:
    {
      DefineOutput(8, AliRDHFCutsDplustoKpipi::Class());
      break;
    }
    default:
    {
      AliFatal("Invalid HF channel.");
      break;
    }
  }

  if (fIsMC) {
    DefineOutput(9, TList::Class());
    DefineOutput(10, TList::Class());
  }
}

//____________________________________________________________________________________________________
AliAnalysisTaskCharmingFemto::~AliAnalysisTaskCharmingFemto() {
  delete fPartColl;
  delete fPairCleaner;
  delete fProtonTrack;
  delete fRDHFCuts;

  if(fApplyML && fMLResponse) {
    delete fMLResponse;
  }
}

//____________________________________________________________________________________________________
void AliAnalysisTaskCharmingFemto::LocalInit()
{
    // Initialization
    switch(fDecChannel) {
      case kDplustoKpipi:
        AliRDHFCutsDplustoKpipi* copyCut = new AliRDHFCutsDplustoKpipi(*(static_cast<AliRDHFCutsDplustoKpipi*>(fRDHFCuts)));
        PostData(8, copyCut);
      break;
    }
  return;
}

//____________________________________________________________________________________________________
void AliAnalysisTaskCharmingFemto::UserExec(Option_t * /*option*/) {
  fInputEvent = static_cast<AliAODEvent*>(InputEvent());
  if (fIsMC)
    fMCEvent = MCEvent();

  // PREAMBLE - CHECK EVERYTHING IS THERE
  if (!fInputEvent) {
    AliError("No Input event");
    return;
  }

  // Protection against the mismatch of candidate TRefs:
  // Check if AOD and corresponding deltaAOD files contain the same number of events.
  // In case of discrepancy the event is rejected.
  if(fAODProtection >= 0) {
    // Protection against different number of events in the AOD and deltaAOD
    // In case of discrepancy the event is rejected.
    int matchingAODdeltaAODlevel = AliRDHFCuts::CheckMatchingAODdeltaAODevents();
    if (matchingAODdeltaAODlevel < 0 || (matchingAODdeltaAODlevel == 0 && fAODProtection == 1)) {
      // AOD/deltaAOD trees have different number of entries || TProcessID do not match while it was required
      return;
    }
  }

  // GET HF CANDIDATE ARRAY
  TClonesArray *arrayHF = nullptr;
  int absPdgMom = 0;
  TString mesonName = "";
  if(!fInputEvent && AODEvent() && IsStandardAOD()) {
    // In case there is an AOD handler writing a standard AOD, use the AOD
    // event in memory rather than the input (ESD) event.
    fInputEvent = dynamic_cast<AliAODEvent*>(AODEvent());
    // in this case the braches in the deltaAOD (AliAOD.VertexingHF.root)
    // have to taken from the AOD event hold by the AliAODExtension
    AliAODHandler* aodHandler = (AliAODHandler*)((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
    if(aodHandler->GetExtensions()) {
      AliAODExtension *ext = (AliAODExtension*)aodHandler->GetExtensions()->FindObject("AliAOD.VertexingHF.root");
      AliAODEvent *aodFromExt = ext->GetAOD();
      switch(fDecChannel) {
        case kDplustoKpipi:
          absPdgMom = 411;
          mesonName = "Dplus";
          arrayHF = dynamic_cast<TClonesArray*>(aodFromExt->GetList()->FindObject("Charm3Prong"));
          break;
      }
    }
  }
  else if(fInputEvent){
    switch(fDecChannel) {
      case kDplustoKpipi:
        absPdgMom = 411;
        mesonName = "Dplus";
        arrayHF = dynamic_cast<TClonesArray*>(fInputEvent->GetList()->FindObject("Charm3Prong"));
        break;
    }
  }

  if(!arrayHF) {
    AliError("Branch not found!\n");
    return;
  }

  if (fIsMC && !fMCEvent) {
    AliError("No MC event");
    return;
  }

  if (!fEvtCuts) {
    AliError("Event Cuts missing");
    return;
  }

  if (!fTrackCutsPartProton || !fTrackCutsPartAntiProton) {
    AliError("Proton Cuts missing");
    return;
  }

  // Reset the pair cleaner
  fPairCleaner->ResetArray();

  // EVENT SELECTION
  fEvent->SetEvent(fInputEvent);
  if (!fEvtCuts->isSelected(fEvent))
    return;

  // PROTON SELECTION
  ResetGlobalTrackReference();
  for (int iTrack = 0; iTrack < fInputEvent->GetNumberOfTracks(); ++iTrack) {
    AliAODTrack *track = static_cast<AliAODTrack *>(fInputEvent->GetTrack(
        iTrack));
    if (!track) {
      AliFatal("No Standard AOD");
      return;
    }
    StoreGlobalTrackReference(track);
  }
  static std::vector<AliFemtoDreamBasePart> protons;
  static std::vector<AliFemtoDreamBasePart> antiprotons;
  protons.clear();
  antiprotons.clear();
  const int multiplicity = fEvent->GetMultiplicity();
  fProtonTrack->SetGlobalTrackInfo(fGTI, fTrackBufferSize);
  for (int iTrack = 0; iTrack < fInputEvent->GetNumberOfTracks(); ++iTrack) {
    AliAODTrack *track = static_cast<AliAODTrack *>(fInputEvent->GetTrack(
        iTrack));
    fProtonTrack->SetTrack(track);
    fProtonTrack->SetInvMass(0.938);
    if (fTrackCutsPartProton->isSelected(fProtonTrack)) {
      protons.push_back(*fProtonTrack);
    }
    if (fTrackCutsPartAntiProton->isSelected(fProtonTrack)) {
      antiprotons.push_back(*fProtonTrack);
    }
  }

  // D MESON SELECTION
  int nCand = arrayHF->GetEntriesFast();

  // check if the train includes the common ML selector for the given charm-hadron species
  AliAnalysisTaskSECharmHadronMLSelector *taskMLSelect = nullptr;
  std::vector<int> chHadIdx{};
  std::vector<std::vector<double> > scoresFromMLSelector{};
  if(fDependOnMLSelector) {
    taskMLSelect = dynamic_cast<AliAnalysisTaskSECharmHadronMLSelector*>(AliAnalysisManager::GetAnalysisManager()->GetTask(Form("MLSelector%s", mesonName.Data())));
    if(!taskMLSelect) {
      AliFatal("ML Selector not present in train and ML models not compiled!");
      return;
    }

    chHadIdx = taskMLSelect->GetSelectedCandidates();
    scoresFromMLSelector = taskMLSelect->GetMLSCores();
  }
  else {
    for (int iCand = 0; iCand < nCand; iCand++) {
      chHadIdx.push_back(iCand);
      scoresFromMLSelector.push_back({});
    }
  }

  static std::vector<AliFemtoDreamBasePart> dplus;
  static std::vector<AliFemtoDreamBasePart> dminus;
  dplus.clear();
  dminus.clear();

  // needed to initialise PID response
  fRDHFCuts->IsEventSelected(fInputEvent);

  AliAODRecoDecayHF *dMeson = nullptr;
  for (size_t iCand = 0; iCand < chHadIdx.size(); iCand++) {

    dMeson = dynamic_cast<AliAODRecoDecayHF*>(arrayHF->UncheckedAt(chHadIdx[iCand]));

    bool unsetVtx = false;
    bool recVtx = false;
    AliAODVertex *origOwnVtx = nullptr;

    int isSelected = IsCandidateSelected(dMeson, absPdgMom, unsetVtx, recVtx, origOwnVtx, scoresFromMLSelector[iCand]);
    if(!isSelected) {
      if (unsetVtx) {
        dMeson->UnsetOwnPrimaryVtx();
      }
      if (recVtx) {
        fRDHFCuts->CleanOwnPrimaryVtx(dMeson, fInputEvent, origOwnVtx);
      }
      continue;
    }

    const double mass = dMeson->InvMass(fDmesonPDGs.size(), &fDmesonPDGs[0]);
    if (dMeson->Charge() > 0) {
      fHistDplusInvMassPt->Fill(dMeson->Pt(), mass);
    } else {
      fHistDminusInvMassPt->Fill(dMeson->Pt(), mass);
    }


    if( MassSelection(mass, dMeson->Pt(), absPdgMom) ) {
      if (dMeson->Charge() > 0) {
        AliFemtoDreamBasePart dplusCand(dMeson, fInputEvent, absPdgMom, fDmesonPDGs);
        if (fIsMC && fMCBeautyRejection
            && dplusCand.GetParticleOrigin()
                == AliFemtoDreamBasePart::kBeauty) {
          if (gRandom->Uniform() > fMCBeautyScalingFactor)
            continue;
        }
        if (fIsMC && fUseTrueDOnly
            && std::abs(dplusCand.GetMCPDGCode()) != absPdgMom) {
          continue;
        }
        dplus.push_back(dplusCand);
        if (!fIsLightweight) {
          fHistDplusInvMassPtSel->Fill(dMeson->Pt(), mass);
          fHistDplusEta->Fill(dMeson->Eta());
          fHistDplusPhi->Fill(dMeson->Phi());
          if (fIsMC) {
            fHistDplusMCPDGPt->Fill(dplusCand.GetPt(),
                                    dplusCand.GetMotherPDG());
            fHistDplusMCOrigin->Fill(dplusCand.GetPt(),
                                     dplusCand.GetParticleOrigin());
            if (dplusCand.GetMCPDGCode() != 0) {
              fHistDplusMCPtRes->Fill(dplusCand.GetPt() - dplusCand.GetMCPt(),
                                      dplusCand.GetPt());
              fHistDplusMCPhiRes->Fill(
                  dplusCand.GetPhi().at(0) - dplusCand.GetMCPhi().at(0),
                  dplusCand.GetPt());
              fHistDplusMCThetaRes->Fill(
                  dplusCand.GetTheta().at(0) - dplusCand.GetMCTheta().at(0),
                  dplusCand.GetPt());
            }
          }
          for (unsigned int iChild = 0; iChild < fDmesonPDGs.size(); iChild++) {
            AliAODTrack *track = (AliAODTrack *) dMeson->GetDaughter(iChild);
            fHistDplusChildPt[iChild]->Fill(track->Pt());
            fHistDplusChildEta[iChild]->Fill(track->Eta());
            fHistDplusChildPhi[iChild]->Fill(track->Phi());
          }
        }
      } else {
        AliFemtoDreamBasePart dminusCand(dMeson, fInputEvent, absPdgMom,
                                         fDmesonPDGs);
        if (fIsMC && fMCBeautyRejection
            && dminusCand.GetParticleOrigin()
                == AliFemtoDreamBasePart::kBeauty) {
          if (gRandom->Uniform() > fMCBeautyScalingFactor)
            continue;
        }
        if (fIsMC && fUseTrueDOnly
            && std::abs(dminusCand.GetMCPDGCode()) != absPdgMom) {
          continue;
        }
        dminus.push_back(dminusCand);
        if (!fIsLightweight) {
          fHistDminusInvMassPtSel->Fill(dMeson->Pt(), mass);
          fHistDminusEta->Fill(dMeson->Eta());
          fHistDminusPhi->Fill(dMeson->Phi());
          if (fIsMC) {
            fHistDminusMCPDGPt->Fill(dminusCand.GetPt(),
                                     dminusCand.GetMotherPDG());
            fHistDminusMCOrigin->Fill(dminusCand.GetPt(),
                                      dminusCand.GetParticleOrigin());
            if (dminusCand.GetMCPDGCode() != 0) {
              fHistDminusMCPtRes->Fill(
                  dminusCand.GetPt() - dminusCand.GetMCPt(),
                  dminusCand.GetPt());
              fHistDminusMCPhiRes->Fill(
                  dminusCand.GetPhi().at(0) - dminusCand.GetMCPhi().at(0),
                  dminusCand.GetPt());
              fHistDminusMCThetaRes->Fill(
                  dminusCand.GetTheta().at(0) - dminusCand.GetMCTheta().at(0),
                  dminusCand.GetPt());
            }
          }
          for (unsigned int iChild = 0; iChild < fDmesonPDGs.size(); iChild++) {
            AliAODTrack *track = (AliAODTrack *) dMeson->GetDaughter(iChild);
            fHistDminusChildPt[iChild]->Fill(track->Pt());
            fHistDminusChildEta[iChild]->Fill(track->Eta());
            fHistDminusChildPhi[iChild]->Fill(track->Phi());
          }
        }
      }
    }

    if (unsetVtx) {
      dMeson->UnsetOwnPrimaryVtx();
    }
    if (recVtx) {
      fRDHFCuts->CleanOwnPrimaryVtx(dMeson, fInputEvent, origOwnVtx);
    }
  }

  // PAIR CLEANING AND FEMTO
  fPairCleaner->CleanTrackAndDecay(&protons, &dplus, 0);
  fPairCleaner->CleanTrackAndDecay(&protons, &dminus, 1);
  fPairCleaner->CleanTrackAndDecay(&antiprotons, &dplus, 2);
  fPairCleaner->CleanTrackAndDecay(&antiprotons, &dminus, 3);

  fPairCleaner->StoreParticle(protons);
  fPairCleaner->StoreParticle(antiprotons);
  fPairCleaner->StoreParticle(dplus);
  fPairCleaner->StoreParticle(dminus);

  fPartColl->SetEvent(fPairCleaner->GetCleanParticles(), fEvent);

  // flush the data
  PostData(1, fQA);
  PostData(2, fEvtHistList);
  PostData(3, fTrackCutHistList);
  PostData(4, fAntiTrackCutHistList);
  PostData(5, fDChargedHistList);
  PostData(6, fResultList);
  PostData(7, fResultQAList);
  if (fIsMC) {
    PostData(9, fTrackCutHistMCList);
    PostData(10, fAntiTrackCutHistMCList);
  }
}

//____________________________________________________________________________________________________
void AliAnalysisTaskCharmingFemto::ResetGlobalTrackReference() {
  // see AliFemtoDreamAnalysis for details
  for (int i = 0; i < fTrackBufferSize; i++) {
    fGTI[i] = 0;
  }
}

//____________________________________________________________________________________________________
void AliAnalysisTaskCharmingFemto::StoreGlobalTrackReference(
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
      printf("Warning! global track info already there!");
      printf("         TPCNcls track1 %u track2 %u",
             (fGTI[trackID])->GetTPCNcls(), track->GetTPCNcls());
      printf("         FilterMap track1 %u track2 %u\n",
             (fGTI[trackID])->GetFilterMap(), track->GetFilterMap());
    }
  }
  (fGTI[trackID]) = track;
}

//____________________________________________________________________________________________________
void AliAnalysisTaskCharmingFemto::UserCreateOutputObjects() {
  fGTI = new AliAODTrack *[fTrackBufferSize];

  fEvent = new AliFemtoDreamEvent(true, !fIsLightweight, fTrigger);

  fProtonTrack = new AliFemtoDreamTrack();
  fProtonTrack->SetUseMCInfo(fIsMC);

  fPairCleaner = new AliFemtoDreamPairCleaner(4, 0,
                                              fConfig->GetMinimalBookingME());
  fPartColl = new AliFemtoDreamPartCollection(fConfig,
                                              fConfig->GetMinimalBookingME());

  fQA = new TList();
  fQA->SetName("QA");
  fQA->SetOwner(kTRUE);

  if (fEvtCuts) {
    fEvtCuts->InitQA();
    if (fEvent->GetEvtCutList() && !fIsLightweight) {
      fQA->Add(fEvent->GetEvtCutList());
    }
    if (fEvtCuts->GetHistList() && !fIsLightweight) {
      fEvtHistList = fEvtCuts->GetHistList();
    } else {
      fEvtHistList = new TList();
      fEvtHistList->SetName("EvtCuts");
      fEvtHistList->SetOwner(true);
    }
  } else {
    AliWarning("Event cuts are missing! \n");
  }

  if (!fConfig->GetMinimalBookingME() && fPairCleaner
      && fPairCleaner->GetHistList()) {
    fQA->Add(fPairCleaner->GetHistList());
  }

  fTrackCutsPartProton->Init("Proton");
  // necessary for the non-min booking case
  fTrackCutsPartProton->SetName("Proton");

  if (fTrackCutsPartProton && fTrackCutsPartProton->GetQAHists()) {
    fTrackCutHistList = fTrackCutsPartProton->GetQAHists();
    if (fIsMC && fTrackCutsPartProton->GetMCQAHists()
        && fTrackCutsPartProton->GetIsMonteCarlo()) {
      fTrackCutHistMCList = fTrackCutsPartProton->GetMCQAHists();
    }
  }

  fTrackCutsPartAntiProton->Init("Anti-proton");
  // necessary for the non-min booking case
  fTrackCutsPartAntiProton->SetName("Anti-proton");

  if (fTrackCutsPartAntiProton && fTrackCutsPartAntiProton->GetQAHists()) {
    fAntiTrackCutHistList = fTrackCutsPartAntiProton->GetQAHists();
    if (fIsMC && fTrackCutsPartAntiProton->GetMCQAHists()
        && fTrackCutsPartAntiProton->GetIsMonteCarlo()) {
      fAntiTrackCutHistMCList = fTrackCutsPartAntiProton->GetMCQAHists();
    }
  }

  // Eventually we might put this in a separate class but for the moment let's just leave it floating around here
  fDChargedHistList  = new TList();
  fDChargedHistList->SetName("DChargedQA");
  fDChargedHistList->SetOwner(true);

  fHistDplusInvMassPt = new TH2F(
      "fHistDplusInvMassPt",
      "; #it{p}_{T} (GeV/#it{c}); #it{M}_{K#pi#pi} (GeV/#it{c}^{2})", 100, 0,
      10, 100, 1.77, 1.97);
  fDChargedHistList->Add(fHistDplusInvMassPt);
  if (!fIsLightweight) {
    fHistDplusInvMassPtSel = new TH2F(
        "fHistDplusInvMassPtSel",
        "; #it{p}_{T} (GeV/#it{c}); #it{M}_{K#pi#pi} (GeV/#it{c}^{2})", 100, 0,
        10, 100, 1.57, 2.17);
    fDChargedHistList->Add(fHistDplusInvMassPtSel);
    fHistDplusEta = new TH1F("fHistDplusEta", ";#eta; Entries", 100, -1, 1);
    fDChargedHistList->Add(fHistDplusEta);
    fHistDplusPhi = new TH1F("fHistDplusPhi", ";#phi; Entries", 100, 0.,
                             2. * TMath::Pi());
    fDChargedHistList->Add(fHistDplusPhi);

    if (fIsMC) {
      fHistDplusMCPDGPt = new TH2F("fHistDplusMCPDGPt",
                                   "; #it{p}_{T} (GeV/#it{c}); PDG Code mother",
                                   250, 0, 25, 5000, 0, 5000);
      fHistDplusMCPtRes =
          new TH2F(
              "fHistDplusMCPtRes",
              "; #it{p}_{T, rec} - #it{p}_{T, gen} (GeV/#it{c}); #it{p}_{T} (GeV/#it{c})",
              101, -0.5, 0.5, 100, 0, 10);
      fHistDplusMCPhiRes = new TH2F(
          "fHistDplusMCPhiRes",
          "; #phi_{rec} - #phi_{gen}; #it{p}_{T} (GeV/#it{c})", 101, -0.025, 0.025,
          100, 0, 10);
      fHistDplusMCThetaRes = new TH2F(
          "fHistDplusMCThetaRes",
          "; #theta_{rec} - #theta_{gen}; #it{p}_{T} (GeV/#it{c})", 101, -0.025,
          0.025, 100, 0, 10);
      fHistDplusMCOrigin = new TH2F("fHistDplusMCOrigin",
                                    "; #it{p}_{T} (GeV/#it{c}); Origin", 100, 0,
                                    10, 7, 0, 7);
      fHistDplusMCOrigin->GetYaxis()->SetBinLabel(1, "kPhysPrimary");
      fHistDplusMCOrigin->GetYaxis()->SetBinLabel(2, "kWeak");
      fHistDplusMCOrigin->GetYaxis()->SetBinLabel(3, "kMaterial");
      fHistDplusMCOrigin->GetYaxis()->SetBinLabel(4, "kFake");
      fHistDplusMCOrigin->GetYaxis()->SetBinLabel(5, "kContamination");
      fHistDplusMCOrigin->GetYaxis()->SetBinLabel(6, "kUnknown");
      fHistDplusMCOrigin->GetYaxis()->SetBinLabel(7, "kBeauty");
      fDChargedHistList->Add(fHistDplusMCPDGPt);
      fDChargedHistList->Add(fHistDplusMCPtRes);
      fDChargedHistList->Add(fHistDplusMCPhiRes);
      fDChargedHistList->Add(fHistDplusMCThetaRes);
      fDChargedHistList->Add(fHistDplusMCOrigin);
    }
  }

  fHistDminusInvMassPt = new TH2F(
      "fHistDminusInvMassPt",
      "; #it{p}_{T} (GeV/#it{c}); #it{M}_{K#pi#pi} (GeV/#it{c}^{2})", 100, 0,
      10, 100, 1.57, 2.17);
  fDChargedHistList->Add(fHistDminusInvMassPt);
  if (!fIsLightweight) {
    fHistDminusInvMassPtSel = new TH2F(
        "fHistDminusInvMassPtSel",
        "; #it{p}_{T} (GeV/#it{c}); #it{M}_{K#pi#pi} (GeV/#it{c}^{2})", 100, 0,
        10, 100, 1.77, 1.97);
    fDChargedHistList->Add(fHistDminusInvMassPtSel);
    fHistDminusEta = new TH1F("fHistDminusEta", ";#eta; Entries", 100, -1, 1);
    fDChargedHistList->Add(fHistDminusEta);
    fHistDminusPhi = new TH1F("fHistDminusPhi", ";#phi; Entries", 100, 0.,
                              2. * TMath::Pi());
    fDChargedHistList->Add(fHistDminusPhi);

    if (fIsMC) {
      fHistDminusMCPDGPt = new TH2F(
          "fHistDminusMCPDGPt", "; #it{p}_{T} (GeV/#it{c}); PDG Code mother",
          250, 0, 25, 5000, 0, 5000);
      fHistDminusMCPtRes =
          new TH2F(
              "fHistDminusMCPtRes",
              "; #it{p}_{T, rec} - #it{p}_{T, gen} (GeV/#it{c}); #it{p}_{T} (GeV/#it{c})",
              101, -0.5, 0.5, 100, 0, 10);
      fHistDminusMCPhiRes = new TH2F(
          "fHistDminusMCPhiRes",
          "; #phi_{rec} - #phi_{gen}; #it{p}_{T} (GeV/#it{c})", 101, -0.025, 0.025,
          100, 0, 10);
      fHistDminusMCThetaRes = new TH2F(
          "fHistDminusMCThetaRes",
          "; #theta_{rec} - #theta_{gen}; #it{p}_{T} (GeV/#it{c})", 101, -0.025,
          0.025, 100, 0, 10);
      fHistDminusMCOrigin = new TH2F("fHistDminusMCOrigin",
                                     "; #it{p}_{T} (GeV/#it{c}); Origin", 100,
                                     0, 10, 7, 0, 7);
      fHistDminusMCOrigin->GetYaxis()->SetBinLabel(1, "kPhysPrimary");
      fHistDminusMCOrigin->GetYaxis()->SetBinLabel(2, "kWeak");
      fHistDminusMCOrigin->GetYaxis()->SetBinLabel(3, "kMaterial");
      fHistDminusMCOrigin->GetYaxis()->SetBinLabel(4, "kFake");
      fHistDminusMCOrigin->GetYaxis()->SetBinLabel(5, "kContamination");
      fHistDminusMCOrigin->GetYaxis()->SetBinLabel(6, "kUnknown");
      fHistDminusMCOrigin->GetYaxis()->SetBinLabel(7, "kBeauty");
      fDChargedHistList->Add(fHistDminusMCPDGPt);
      fDChargedHistList->Add(fHistDminusMCPtRes);
      fDChargedHistList->Add(fHistDminusMCPhiRes);
      fDChargedHistList->Add(fHistDminusMCThetaRes);
      fDChargedHistList->Add(fHistDminusMCOrigin);
    }
  }


  if (!fIsLightweight) {
    std::vector<TString> nameVec;
    if (fDecChannel == kDplustoKpipi) {
      nameVec = { {"K", "Pi1", "Pi2"}};
    }
    for (unsigned int iChild = 0; iChild < fDmesonPDGs.size(); ++iChild) {
      fHistDplusChildPt[iChild] = new TH1F(
          TString::Format("fHistDplusChildPt_%s", nameVec.at(iChild).Data()),
          "; #it{p}_{T} (GeV/#it{c}); Entries", 250, 0, 25);
      fHistDplusChildEta[iChild] = new TH1F(
          TString::Format("fHistDplusChildEta_%s", nameVec.at(iChild).Data()),
          "; #eta; Entries", 100, -1, 1);
      fHistDplusChildPhi[iChild] = new TH1F(
          TString::Format("fHistDplusChildPhi_%s", nameVec.at(iChild).Data()),
          "; #phi; Entries", 100, 0, 2. * TMath::Pi());
      fDChargedHistList->Add(fHistDplusChildPt[iChild]);
      fDChargedHistList->Add(fHistDplusChildEta[iChild]);
      fDChargedHistList->Add(fHistDplusChildPhi[iChild]);
      fHistDminusChildPt[iChild] = new TH1F(
          TString::Format("fHistDminusChildPt_%s", nameVec.at(iChild).Data()),
          "; #it{p}_{T} (GeV/#it{c}); Entries", 250, 0, 25);
      fHistDminusChildEta[iChild] = new TH1F(
          TString::Format("fHistDminusChildEta_%s", nameVec.at(iChild).Data()),
          "; #eta; Entries", 100, -1, 1);
      fHistDminusChildPhi[iChild] = new TH1F(
          TString::Format("fHistDminusChildPhi_%s", nameVec.at(iChild).Data()),
          "; #phi; Entries", 100, 0, 2. * TMath::Pi());
      fDChargedHistList->Add(fHistDminusChildPt[iChild]);
      fDChargedHistList->Add(fHistDminusChildEta[iChild]);
      fDChargedHistList->Add(fHistDminusChildPhi[iChild]);
    }
  }
  

  if (fPartColl && fPartColl->GetHistList()) {
    fResultList = fPartColl->GetHistList();
  }
  if (!fConfig->GetMinimalBookingME() && fPartColl && fPartColl->GetQAList()) {
    fResultQAList = fPartColl->GetQAList();
  } else {
    fResultQAList = new TList();
    fResultQAList->SetName("ResultsQA");
    fResultQAList->SetOwner(true);
  }

  //ML model
  if(fApplyML) {
    if(!fDependOnMLSelector) {
      switch(fDecChannel) {
        case kDplustoKpipi:
            fMLResponse = new AliHFMLResponseDplustoKpipi("DplustoKpipiMLResponse", "DplustoKpipiMLResponse", fConfigPath.Data());
            fMLResponse->MLResponseInit();
        break;
      }
    }
    else {
      std::string configLocalPath = AliMLModelHandler::ImportFile(fConfigPath.Data());
      YAML::Node nodeList;
      try {
        nodeList = YAML::LoadFile(configLocalPath);
      } catch (std::exception &e) {
        AliFatal(Form("Yaml-ccp error: %s! Exit", e.what()));
      }
      fPtLimsML = nodeList["BINS"].as<vector<float> >();

      for (const auto &model : nodeList["MODELS"]) {
        fMLScoreCuts.push_back(model["cut"].as<std::vector<double> >());
        fMLOptScoreCuts.push_back(model["cut_opt"].as<std::vector<std::string> >());
      }
    }
  }

  PostData(1, fQA);
  PostData(2, fEvtHistList);
  PostData(3, fTrackCutHistList);
  PostData(4, fAntiTrackCutHistList);
  PostData(5, fDChargedHistList);
  PostData(6, fResultList);
  PostData(7, fResultQAList);
  if (fIsMC) {
    PostData(9, fTrackCutHistMCList);
    PostData(10, fAntiTrackCutHistMCList);
  }
}

//________________________________________________________________________
int AliAnalysisTaskCharmingFemto::IsCandidateSelected(AliAODRecoDecayHF *&dMeson, int absPdgMom, bool &unsetVtx, bool &recVtx, AliAODVertex *&origOwnVtx, std::vector<double> scores) {

  if(!dMeson) {
    return 0;
  }

  bool isSelBit = true;
  switch(fDecChannel) {
    case kDplustoKpipi:
      isSelBit = dMeson->HasSelectionBit(AliRDHFCuts::kDplusCuts);
      if(!isSelBit) {
        return 0;
      }
      break;
  }

  unsetVtx = false;
  if (!dMeson->GetOwnPrimaryVtx())
  {
    dMeson->SetOwnPrimaryVtx(dynamic_cast<AliAODVertex *>(fInputEvent->GetPrimaryVertex()));
    unsetVtx = true;
    // NOTE: the own primary vertex should be unset, otherwise there is a memory leak
    // Pay attention if you use continue inside this loop!!!
  }

  double ptD = dMeson->Pt();
  double yD = dMeson->Y(absPdgMom);
  int ptbin = fRDHFCuts->PtBin(ptD);
  if(ptbin < 0) {
    if (unsetVtx) {
      dMeson->UnsetOwnPrimaryVtx();
    }
    return 0;
  }

  bool isFidAcc = fRDHFCuts->IsInFiducialAcceptance(ptD, yD);
  if(!isFidAcc) {
    if (unsetVtx) {
      dMeson->UnsetOwnPrimaryVtx();
    }
    return 0;
  }

  int isSelected = fRDHFCuts->IsSelected(dMeson, AliRDHFCuts::kAll, fInputEvent);
  if(!isSelected) {
    if (unsetVtx) {
      dMeson->UnsetOwnPrimaryVtx();
    }
    return 0;
  }

  // ML application
  if(fApplyML) {
    if(!fDependOnMLSelector) { //direct application
      AliAODPidHF* pidHF = fRDHFCuts->GetPidHF();
      bool isMLsel = true;
      std::vector<double> modelPred{};
      switch(fDecChannel) {
        case kDplustoKpipi:
          isMLsel = fMLResponse->IsSelectedMultiClass(modelPred, dMeson, fInputEvent->GetMagneticField(), pidHF);
          if(!isMLsel) {
            isSelected = 0;
          }
          break;
      }
    }
    else { // read result from common task
      std::vector<float>::iterator low = std::lower_bound(fPtLimsML.begin(), fPtLimsML.end(), ptD);
      int bin = low - fPtLimsML.begin() - 1;
      if(bin < 0)
        bin = 0;
      else if(bin > fPtLimsML.size()-2)
        bin = fPtLimsML.size()-2;
      for(size_t iScore = 0; iScore < scores.size(); iScore++) {
        if((fMLOptScoreCuts[bin][iScore] == "upper" && scores[iScore] > fMLScoreCuts[bin][iScore]) || (fMLOptScoreCuts[bin][iScore] == "lower" && scores[iScore] < fMLScoreCuts[bin][iScore])){
          isSelected = 0;
          break;
        }
      }
    }
  }

  recVtx = false;
  origOwnVtx = nullptr;

  if (fRDHFCuts->GetIsPrimaryWithoutDaughters())
  {
    if (dMeson->GetOwnPrimaryVtx()) {
      origOwnVtx = new AliAODVertex(*dMeson->GetOwnPrimaryVtx());
    }
    if (fRDHFCuts->RecalcOwnPrimaryVtx(dMeson, fInputEvent)) {
      recVtx = true;
    }
    else {
      fRDHFCuts->CleanOwnPrimaryVtx(dMeson, fInputEvent, origOwnVtx);
    }
  }
  
  return isSelected;
}

//____________________________________________________________________________________________________
bool AliAnalysisTaskCharmingFemto::MassSelection(const double mass,
                                                 const double pt,
                                                 const int pdg) {
  // simple parametrisation from D+ in 5.02 TeV
  double massMean = TDatabasePDG::Instance()->GetParticle(pdg)->Mass() + 0.0025;  // mass shift observed in all Run2 data samples for all
                                                                                  // D-meson species
  double massWidth = 0.;
  switch (fDecChannel) {
    case kDplustoKpipi:
      if (fSystem == kpp5TeV) {
        massWidth = 0.0057 + pt * 0.00066;
      } else if (fSystem == kpp13TeV) {
        massWidth = 0.006758 + pt * 0.0005124;
      }
      break;
  }

  // select D mesons mass window
  if (fMassSelectionType == kSignal) {
    fLowerMassSelection = massMean - fNSigmaMass * massWidth;
    fUpperMassSelection = massMean + fNSigmaMass * massWidth;
  } else if ( fMassSelectionType == kSidebandLeft) {
    fLowerMassSelection = massMean - fNSigmaOffsetSideband * massWidth - fSidebandWidth;
    fUpperMassSelection = massMean - fNSigmaOffsetSideband * massWidth;
  } else if ( fMassSelectionType == kSidebandRight) {
    fLowerMassSelection = massMean + fNSigmaOffsetSideband * massWidth;
    fUpperMassSelection = massMean + fNSigmaOffsetSideband * massWidth + fSidebandWidth;

    // additional removal of D*
    if ( mass > fLowerDstarRemoval && mass < fUpperDstarRemoval) {
      return false;
    }
  }


  if (mass > fLowerMassSelection && mass < fUpperMassSelection) {
    return true;
  }

  return false;
}
