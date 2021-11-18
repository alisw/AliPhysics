#include "AliAnalysisTaskValeNanoTreeLPhi.h"
#include "AliFemtoDreamBasePart.h"
#include "AliLog.h"
#include "AliNanoAODTrack.h"
#include "AliVEvent.h"
#include "TH1F.h"
#include "TList.h"
#include "AliAODMCParticle.h"

ClassImp(AliAnalysisTaskValeNanoTreeLPhi)
    AliAnalysisTaskValeNanoTreeLPhi::AliAnalysisTaskValeNanoTreeLPhi()
    : AliAnalysisTaskSE(),
      fIsMC(false),
      fUseOMixing(false),
      fInvMassCutSBdown(0.0),
      fInvMassCutSBup(0.0),
      fTrigger(AliVEvent::kINT7),
      fQA(nullptr),
      fEvtList(nullptr),
      fLambdaList(nullptr),
      fAntiLambdaList(nullptr),
      fKaonPlusList(nullptr),
      fKaonMinusList(nullptr),
      fPhiList(nullptr),
      fResults(nullptr),
      fResultsQA(nullptr),
      fTree(0),
      fInputEvent(nullptr),
      fEvent(nullptr),
      fTrack(nullptr),
      fLambda(nullptr),
      fPhiParticle(nullptr),
      fEventCuts(nullptr),
      fPosKaonCuts(nullptr),
      fNegKaonCuts(nullptr),
      fPhiCuts(nullptr),
      fLambdaCuts(nullptr),
      fAntiLambdaCuts(nullptr),
      fConfig(nullptr),
      fPairCleaner(nullptr),
      fPartColl(nullptr),
      fSample(nullptr),
      fGTI(nullptr),
      fTrackBufferSize()
{
}

AliAnalysisTaskValeNanoTreeLPhi::AliAnalysisTaskValeNanoTreeLPhi(
    const char *name, bool isMC)
    : AliAnalysisTaskSE(name),
      fIsMC(isMC),
      fUseOMixing(false),
      fInvMassCutSBdown(0.0),
      fInvMassCutSBup(0.0),
      fTrigger(AliVEvent::kINT7),
      fQA(nullptr),
      fEvtList(nullptr),
      fLambdaList(nullptr),
      fAntiLambdaList(nullptr),
      fKaonPlusList(nullptr),
      fKaonMinusList(nullptr),
      fPhiList(nullptr),
      fResults(nullptr),
      fResultsQA(nullptr),
      fTree(0),
      fInputEvent(nullptr),
      fEvent(nullptr),
      fTrack(nullptr),
      fLambda(nullptr),
      fPhiParticle(nullptr),
      fEventCuts(nullptr),
      fPosKaonCuts(nullptr),
      fNegKaonCuts(nullptr),
      fPhiCuts(nullptr),
      fLambdaCuts(nullptr),
      fAntiLambdaCuts(nullptr),
      fConfig(nullptr),
      fPairCleaner(nullptr),
      fPartColl(nullptr),
      fSample(nullptr),
      fGTI(nullptr),
      fTrackBufferSize(2000)
{
  DefineOutput(1, TList::Class());  //Output for the Event Class and Pair Cleaner
  DefineOutput(2, TList::Class());  //Output for the Event Cuts
  DefineOutput(3, TList::Class());  //Output for the Lambda Cuts
  DefineOutput(4, TList::Class());  //Output for the AntiLambda Cuts
  DefineOutput(5, TList::Class());  //Output for the KaonPlus Cuts
  DefineOutput(6, TList::Class());  //Output for the KaonMinus Cuts
  DefineOutput(7, TList::Class());  //Output for the Phi Cuts
  DefineOutput(8, TList::Class());  //Output for the Results
  DefineOutput(9, TList::Class());  //Output for the Results QA
  DefineOutput(10, TList::Class()); //Output for Tree for Lambda
}

AliAnalysisTaskValeNanoTreeLPhi::~AliAnalysisTaskValeNanoTreeLPhi()
{
  delete fEvent;
  delete fTrack;
  delete fLambda;
  delete fPosKaonCuts;
  delete fNegKaonCuts;
  delete fLambdaCuts;
  delete fAntiLambdaCuts;
  delete fPairCleaner;
  delete fPartColl;
  delete fTree;
}

void AliAnalysisTaskValeNanoTreeLPhi::UserCreateOutputObjects()
{

  fEvent = new AliFemtoDreamEvent(true, true, fTrigger);
  fEvent->SetMultiplicityEstimator(fConfig->GetMultiplicityEstimator());

  fTrack = new AliFemtoDreamTrack();
  fTrack->SetUseMCInfo(
      fPosKaonCuts->GetIsMonteCarlo() || fNegKaonCuts->GetIsMonteCarlo());

  fLambda = new AliFemtoDreamv0();
  fLambda->SetUseMCInfo(
      fLambdaCuts->GetIsMonteCarlo() || fAntiLambdaCuts->GetIsMonteCarlo());
  fLambda->SetPDGCode(3122);
  fLambda->SetPDGDaughterPos(2212);
  fLambda->GetPosDaughter()->SetUseMCInfo(
      fLambdaCuts->GetIsMonteCarlo() || fAntiLambdaCuts->GetIsMonteCarlo());
  fLambda->SetPDGDaughterNeg(211);
  fLambda->GetNegDaughter()->SetUseMCInfo(
      fLambdaCuts->GetIsMonteCarlo() || fAntiLambdaCuts->GetIsMonteCarlo());

  fPhiParticle = new AliFemtoDreamv0();
  fPhiParticle->SetPDGCode(fPhiCuts->GetPDGv0());
  fPhiParticle->SetUseMCInfo(fIsMC);
  fPhiParticle->SetPDGDaughterPos(
      fPhiCuts->GetPDGPosDaug()); // order +sign doesnt play a role
  fPhiParticle->GetPosDaughter()->SetUseMCInfo(fIsMC);
  fPhiParticle->SetPDGDaughterNeg(
      fPhiCuts->GetPDGNegDaug()); // only used for MC Matching
  fPhiParticle->GetNegDaughter()->SetUseMCInfo(fIsMC);

  fGTI = new AliVTrack *[fTrackBufferSize];

  if (!fEventCuts)
  {
    AliFatal("Event Cuts not set!");
  }
  fEventCuts->InitQA();

  fQA = new TList();
  fQA->SetOwner();
  fQA->SetName("QA");
  fQA->Add(fEvent->GetEvtCutList());

  if (!fEventCuts->GetMinimalBooking())
  {
    fEvtList = fEventCuts->GetHistList();
  }
  else
  {
    fEvtList = new TList();
    fEvtList->SetName("EventCuts");
    fEvtList->SetOwner();
  }

  if (!fPosKaonCuts)
  {
    AliFatal("Track Cuts for Particle One not set!");
  }
  fPosKaonCuts->Init();
  fPosKaonCuts->SetName("KaonPlus");

  if (fPosKaonCuts->GetIsMonteCarlo())
  {
    fPosKaonCuts->SetMCName("MCParticle1");
    fKaonPlusList->Add(fPosKaonCuts->GetMCQAHists());
  }

  if (!fNegKaonCuts)
  {
    AliFatal("Track Cuts for Particle One not set!");
  }
  fNegKaonCuts->Init();
  fNegKaonCuts->SetName("KaonMinus");
  if (fNegKaonCuts->GetIsMonteCarlo())
  {
    fNegKaonCuts->SetMCName("MCParticle2");
    fKaonMinusList->Add(fNegKaonCuts->GetMCQAHists());
  }

  if (!fPhiCuts)
  {
    AliFatal("Cuts for the phi not set!");
  }
  fPhiCuts->Init();
  fPhiCuts->SetName("Phi");
  if (fPhiCuts->GetIsMonteCarlo())
  {
    fPhiCuts->SetMCName("MCPhi");
    fPhiList->Add(fPhiCuts->GetMCQAHists());
  }

  if (!fLambdaCuts)
  {
    AliFatal("Track Cuts for Lambda not set!");
  }
  fLambdaCuts->Init();
  fLambdaCuts->SetName("Lambda");
  if (fLambdaCuts->GetIsMonteCarlo())
  {
    fLambdaCuts->SetMCName("MCLambda");
    fLambdaList->Add(fLambdaCuts->GetMCQAHists());
  }

  if (!fAntiLambdaCuts)
  {
    AliFatal("Track Cuts for AntiLambda not set!");
  }
  fAntiLambdaCuts->Init();
  fAntiLambdaCuts->SetName("AntiLambda");
  if (fAntiLambdaCuts->GetIsMonteCarlo())
  {
    fAntiLambdaCuts->SetMCName("MCAntiProton");
    fAntiLambdaList->Add(fAntiLambdaCuts->GetMCQAHists());
  }

  fKaonPlusList = fPosKaonCuts->GetQAHists();
  fKaonMinusList = fNegKaonCuts->GetQAHists();
  fLambdaList = fLambdaCuts->GetQAHists();
  fAntiLambdaList = fAntiLambdaCuts->GetQAHists();
  fPhiList = fPhiCuts->GetQAHists();

  fPairCleaner =
      new AliFemtoDreamPairCleaner(0, 2, fConfig->GetMinimalBookingME());
  fPartColl =
      new AliFemtoDreamPartCollection(fConfig, fConfig->GetMinimalBookingME());

  fResultsQA = new TList();
  fResultsQA->SetOwner();
  fResultsQA->SetName("ResultsQA");
  if (fConfig->GetUseEventMixing())
  {
    fResults = fPartColl->GetHistList();
    if (!fConfig->GetMinimalBookingME())
    {
      fResultsQA->Add(fPartColl->GetQAList());
      fResultsQA->Add(fPairCleaner->GetHistList());
    }
  }
  else
  {
    fResults = new TList();
    fResults->SetOwner();
    fResults->SetName("Results");
  }

  ///Creating the entries for the Tree
  fTree = new TTree("LPhiTree", "LPhiTree");
  fTree->Branch("Vz", &fTVz, "fTVz/F");
  fTree->Branch("Mult", &fTMult, "fTMult/I");
  fTree->Branch("sT", &fSpher, "fSpher/F");
  ///Lambda
  fTree->Branch("nLambda", &fNumLambda, "fNumLambda/I");
  fTree->Branch("v0_pT", &v0_pT, "v0_pt[fNumLambda]/F");
  fTree->Branch("v0_decvtx", &v0_decvtx, "v0_decvtx[fNumLambda]/F");
  fTree->Branch("v0_tranrad", &v0_tranrad, "v0_tranrad[fNumLambda]/F");
  fTree->Branch("v0Daugh_dcaDecvtx", &v0Daugh_dcaDecvtx, "v0Daugh_dcaDecvtx[fNumLambda]/F");
  fTree->Branch("v0_cpa", &v0_cpa, "v0_cpa[fNumLambda]/F");
  ///Daughters
  fTree->Branch("Daugh_eta", &Daugh_eta, "Daugh_eta[fNumLambda]/F");
  fTree->Branch("Daugh_nTpcCls", &Daugh_nTpcCls, "Daugh_nTpcCls[fNumLambda]/I");
  fTree->Branch("Daugh_dca", &Daugh_dca, "Daugh_dca[fNumLambda]/F");
  fTree->Branch("Daugh_nSigma", &Daugh_nSigma, "Daugh_nSigma[fNumLambda]/");
  ///Phi
  fTree->Branch("nPhi", &fNumPhi, "fNumPhi/I");

  PostData(1, fQA);
  PostData(2, fEvtList);
  PostData(3, fLambdaList);
  PostData(4, fAntiLambdaList);
  PostData(5, fKaonPlusList);
  PostData(6, fKaonMinusList);
  PostData(7, fPhiList);
  PostData(8, fResults);
  PostData(9, fResultsQA);
  PostData(10, fTree);
}

void AliAnalysisTaskValeNanoTreeLPhi::UserExec(Option_t *)
{
  AliVEvent *fInputEvent = InputEvent();

  // PREAMBLE - CHECK EVERYTHING IS THERE
  if (!fInputEvent)
  {
    AliError("No Input event");
    return;
  }

  // EVENT SELECTION
  fEvent->SetEvent(fInputEvent);
  if (!fEventCuts->isSelected(fEvent))
  {
    return;
  }
  /// Event properties for tree
  fTRunNumber = 0.; //For NanoAOD filtering trains <100, no info
  Double_t PrimVtx[3];
  fInputEvent->GetPrimaryVertex()->GetXYZ(PrimVtx);
  fTVz = PrimVtx[2];
  fTMult = fEvent->GetMultiplicity();
  fSpher = fEvent->GetSpher();

  /// initialization tree
  fNumLambda = 0;
  fNumPhi = 0;
  for (int i = 0; i < nMaxLambda; i++)
  {
    v0_pT[i] = -1000000.;
    v0_decvtx[i] = -1000000.;
    v0_tranrad[i] = -1000000.;
    v0Daugh_dcaDecvtx[i] = -1000000.;
    v0_cpa[i] = -1000000.;
    Daugh_eta[i] = -1000000.;
    Daugh_nTpcCls[i] = -1000000;
    Daugh_dca[i] = -1000000.;
    Daugh_nSigma[i] = -1000000.;
  }

  ResetGlobalTrackReference();
  for (int iTrack = 0; iTrack < fInputEvent->GetNumberOfTracks(); ++iTrack)
  {
    AliVTrack *track = static_cast<AliVTrack *>(fInputEvent->GetTrack(iTrack));
    if (!track)
    {
      AliFatal("No Standard AOD");
      return;
    }
    StoreGlobalTrackReference(track);
  }
  fTrack->SetGlobalTrackInfo(fGTI, fTrackBufferSize);

  std::vector<AliFemtoDreamBasePart> KaonPlus;
  std::vector<AliFemtoDreamBasePart> KaonMinus;
  std::vector<AliFemtoDreamBasePart> PhiParticles;
  std::vector<AliFemtoDreamBasePart> Lambdas;
  std::vector<AliFemtoDreamBasePart> AntiLambdas;

  AliAODEvent *aodEvt = dynamic_cast<AliAODEvent *>(fInputEvent);
  fLambda->SetGlobalTrackInfo(fGTI, fTrackBufferSize);
  for (int iv0 = 0;
       iv0 < static_cast<TClonesArray *>(aodEvt->GetV0s())->GetEntriesFast();
       ++iv0)
  {
    Bool_t IsLambda = kFALSE;
    Bool_t IsAntiLambda = kFALSE;

    AliAODv0 *casc = aodEvt->GetV0(iv0);
    fLambda->Setv0(fInputEvent, casc);
    if (fLambdaCuts->isSelected(fLambda))
    {
      IsLambda = kTRUE;
      Lambdas.push_back(*fLambda);
    }
    if (fAntiLambdaCuts->isSelected(fLambda))
    {
      IsAntiLambda = kTRUE;
      AntiLambdas.push_back(*fLambda);
    }

    if (IsLambda || IsAntiLambda)
      FillLambda(fLambda);
  }

  static float massKaon =
      TDatabasePDG::Instance()->GetParticle(fPosKaonCuts->GetPDGCode())->Mass();

  for (int iTrack = 0; iTrack < fInputEvent->GetNumberOfTracks(); ++iTrack)
  {
    AliVTrack *track = static_cast<AliVTrack *>(fInputEvent->GetTrack(iTrack));

    if (!track)
      continue;
    fTrack->SetTrack(track, fInputEvent);
    fTrack->SetInvMass(massKaon);
    if (fPosKaonCuts->isSelected(fTrack))
    {
      KaonPlus.push_back(*fTrack);
    }
    if (fNegKaonCuts->isSelected(fTrack))
    {
      KaonMinus.push_back(*fTrack);
    }
  }

  fPhiParticle->SetGlobalTrackInfo(fGTI, fTrackBufferSize);
  for (const auto &posK : KaonPlus)
  {
    for (const auto &negK : KaonMinus)
    {
      Bool_t IsPhi = kFALSE;
      fPhiParticle->Setv0(posK, negK, fInputEvent, false, false, true);
      fPhiParticle->SetParticleOrigin(AliFemtoDreamBasePart::kPhysPrimary);

      if (fPhiCuts->isSelected(fPhiParticle))
      {
        IsPhi = kTRUE;
        fPhiParticle->SetCPA(
            gRandom->Uniform()); // cpacode needed for CleanDecay v0;
        PhiParticles.push_back(*fPhiParticle);
      }
      if (IsPhi)
        FillPhi(fPhiParticle);
    }
  }

  if (fNumLambda > 0 && fNumPhi > 0)
    fTree->Fill();

  fPairCleaner->CleanDecayAndDecay(&Lambdas, &PhiParticles, 0);
  fPairCleaner->CleanDecayAndDecay(&AntiLambdas, &PhiParticles, 1);
  fPairCleaner->ResetArray();
  fPairCleaner->StoreParticle(Lambdas);      // 0
  fPairCleaner->StoreParticle(AntiLambdas);  // 1
  fPairCleaner->StoreParticle(PhiParticles); // 2
  fPairCleaner->StoreParticle(KaonPlus);     // 3
  fPairCleaner->StoreParticle(KaonMinus);    // 4

  if (fPairCleaner->GetCounter() > 0)
  {
    fPartColl->SetEvent(fPairCleaner->GetCleanParticles(),
                        fEvent->GetZVertex(), fEvent->GetRefMult08(),
                        fEvent->GetV0MCentrality());
  }

  PostData(1, fQA);
  PostData(2, fEvtList);
  PostData(3, fLambdaList);
  PostData(4, fAntiLambdaList);
  PostData(5, fKaonPlusList);
  PostData(6, fKaonMinusList);
  PostData(7, fPhiList);
  PostData(8, fResults);
  PostData(9, fResultsQA);
  PostData(10, fTree);
}

void AliAnalysisTaskValeNanoTreeLPhi::ResetGlobalTrackReference()
{
  for (UShort_t i = 0; i < fTrackBufferSize; i++)
  {
    fGTI[i] = 0;
  }
}

void AliAnalysisTaskValeNanoTreeLPhi::StoreGlobalTrackReference(
    AliVTrack *track)
{
  // see AliFemtoDreamAnalysis for details
  AliNanoAODTrack *nanoTrack = dynamic_cast<AliNanoAODTrack *>(track);
  const int trackID = track->GetID();
  if (trackID < 0)
  {
    return;
  }
  if (trackID >= fTrackBufferSize)
  {
    printf("Warning: track ID too big for buffer: ID: %d, buffer %d\n", trackID,
           fTrackBufferSize);
    return;
  }

  if (fGTI[trackID])
  {
    if ((!nanoTrack->GetFilterMap()) && (!track->GetTPCNcls()))
    {
      return;
    }
    if (dynamic_cast<AliNanoAODTrack *>(fGTI[trackID])->GetFilterMap() ||
        fGTI[trackID]->GetTPCNcls())
    {
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

Bool_t AliAnalysisTaskValeNanoTreeLPhi::FillLambda(AliFemtoDreamv0 *TheV0)
{
  Bool_t Filled = kFALSE;
  Float_t decvtx = TMath::Sqrt(pow(TheV0->GetDCAv0Vtx(0), 2) + pow(TheV0->GetDCAv0Vtx(1), 2) + pow(TheV0->GetDCAv0Vtx(2), 2));
  v0_pT[fNumLambda] = TheV0->GetPt();
  v0_decvtx[fNumLambda] = decvtx;
  v0_tranrad[fNumLambda] = TheV0->GetTransverseRadius();
  v0Daugh_dcaDecvtx[fNumLambda] = TheV0->GetDCADaugPosVtx();
  v0Daugh_dcaDecvtx[fNumLambda] = TheV0->GetDCADaugNegVtx();
  v0_cpa[fNumLambda] = TheV0->GetCPA();
  AliFemtoDreamTrack *PosDaugh = TheV0->GetPosDaughter();
  AliFemtoDreamTrack *NegDaugh = TheV0->GetNegDaughter();

  Daugh_eta[fNumLambda] = PosDaugh->GetEta().at(1);
  Daugh_eta[fNumLambda] = NegDaugh->GetEta().at(2);
  Daugh_dca[fNumLambda] = TheV0->GetDaugDCA();
  Daugh_nSigma[fNumLambda] = PosDaugh->GetnSigmaTPC((int)(AliPID::kProton));
  Daugh_nSigma[fNumLambda] = NegDaugh->GetnSigmaTPC((int)(AliPID::kPion));
  Daugh_nSigma[fNumLambda] = NegDaugh->GetnSigmaTPC((int)(AliPID::kProton));
  Daugh_nSigma[fNumLambda] = PosDaugh->GetnSigmaTPC((int)(AliPID::kPion));
  fNumLambda++;
  Filled = kTRUE;
  return Filled;
}

Bool_t AliAnalysisTaskValeNanoTreeLPhi::FillPhi(AliFemtoDreamv0 *TheV0) {
  Bool_t Filled = kFALSE;
  fNumPhi++;
  Filled = kTRUE;
  return Filled;
}
