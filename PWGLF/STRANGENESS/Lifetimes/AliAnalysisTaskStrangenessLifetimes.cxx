#include "AliAnalysisTaskStrangenessLifetimes.h"

#include <Riostream.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TList.h>
#include <TMath.h>

#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliExternalTrackParam.h"
#include "AliInputEventHandler.h"
#include "AliLightV0vertexer.h"
#include "AliPIDResponse.h"
#include "AliV0vertexer.h"

#include "AliAnalysisUtils.h"
#include "AliEventCuts.h"

using Lifetimes::MiniV0;
using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskStrangenessLifetimes);

namespace {
constexpr double Sq(double x) { return x * x; }
}  // namespace

AliAnalysisTaskStrangenessLifetimes::AliAnalysisTaskStrangenessLifetimes(
    std::string name)
    : AliAnalysisTaskSE(name.data()),
      fEventCuts{},
      fListHist{nullptr},
      fTreeV0{nullptr},
      fPIDResponse{nullptr},
      fESDtrackCuts{nullptr},
      fDoV0Refit{true},
      fUseLightVertexer{true},
      fV0VertexerSels{33., 0.02, 0.02, 2.0, 0.95, 1.0, 200.},
      fLambdaMassMean{1.116, 0., 0., 0., 0.},
      fLambdaMassSigma{0.002, 0., 0., 0.},
      fMinPtToSave{0.55},
      fMaxPtToSave{100},
      fV0vector{},
      fMultiplicity{} {
  // Standard output
  DefineOutput(1, TList::Class());  // Basic Histograms
  DefineOutput(2, TTree::Class());  // V0 Tree output
}

AliAnalysisTaskStrangenessLifetimes::~AliAnalysisTaskStrangenessLifetimes() {
  if (fListHist) {
    delete fListHist;
    fListHist = 0x0;
  }

  if (fTreeV0) {
    delete fTreeV0;
    fTreeV0 = 0x0;
  }
}

void AliAnalysisTaskStrangenessLifetimes::UserCreateOutputObjects() {
  fTreeV0 = new TTree("fTreeV0", "V0 Candidates");
  fTreeV0->Branch("fMultiplicity", &fMultiplicity, "fMultiplicity/F");
  fTreeV0->Branch("V0s", &fV0vector);

  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *inputHandler =
      (AliInputEventHandler *)(man->GetInputEventHandler());
  fPIDResponse = inputHandler->GetPIDResponse();
  inputHandler->SetNeedField();

  if (!fESDtrackCuts) {
    fESDtrackCuts =
        AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(true, false);
    fESDtrackCuts->SetPtRange(0.15);  // adding pt cut
    fESDtrackCuts->SetEtaRange(-1.0, 1.0);
  }

  fListHist = new TList();
  fListHist->SetOwner();
  fEventCuts.AddQAplotsToList(fListHist);

  PostData(1, fListHist);
  PostData(2, fTreeV0);

}  // end UserCreateOutputObjects

void AliAnalysisTaskStrangenessLifetimes::UserExec(Option_t *) {
  AliESDEvent *esdEvent = dynamic_cast<AliESDEvent *>(InputEvent());
  if (!esdEvent) {
    ::Fatal("AliAnalysisTaskStrangenessLifetimes::UserExec",
            "AliESDEvent not found.");
    return;
  }

  double magneticField = esdEvent->GetMagneticField();

  if (!fEventCuts.AcceptEvent(esdEvent)) {
    PostData(1, fListHist);
    PostData(2, fTreeV0);
    return;
  }

  double primaryVertex[3];
  fMultiplicity = fEventCuts.GetCentrality();
  fEventCuts.GetPrimaryVertex()->GetXYZ(primaryVertex);

  // Variable definition
  double dcaPosToPrimVertex = 0, dcaNegToPrimVertex = 0;
  double lV0CosineOfPointingAngle = 0;
  double v0Pt = 0;

  double fMinV0Pt = 0;
  double fMaxV0Pt = 100;

  // Only reset if not using on-the-fly (or else nothing passes)
  esdEvent->ResetV0s();

  // Decide between regular and light vertexer (default: light)
  if (!fUseLightVertexer) {
    // Instantiate vertexer object
    AliV0vertexer lV0vtxer;
    // Set Cuts
    lV0vtxer.SetDefaultCuts(fV0VertexerSels);
    lV0vtxer.SetCuts(fV0VertexerSels);
    // Redo
    lV0vtxer.Tracks2V0vertices(esdEvent);
  } else {
    // Instantiate vertexer object
    AliLightV0vertexer lV0vtxer;
    // Set do or don't do V0 refit for improved precision
    lV0vtxer.SetDoRefit(false);
    if (fDoV0Refit) lV0vtxer.SetDoRefit(true);
    // Set Cuts
    lV0vtxer.SetDefaultCuts(fV0VertexerSels);
    lV0vtxer.SetCuts(fV0VertexerSels);
    // Redo
    lV0vtxer.Tracks2V0vertices(esdEvent);
  }

  fV0vector.clear();
  for (int iV0 = 0; iV0 < esdEvent->GetNumberOfV0s();
       iV0++) {  // This is the begining of the V0 loop (we analyse only offline
                 // V0s)
    AliESDv0 *v0 = ((AliESDEvent *)esdEvent)->GetV0(iV0);
    if (!v0) continue;
    if (v0->GetOnFlyStatus() != 0) continue;

    // Remove like-sign (will not affect offline V0 candidates!)
    if (v0->GetParamN()->Charge() * v0->GetParamP()->Charge() > 0) continue;

    double decayVtx[3];
    v0->GetXYZ(decayVtx[0], decayVtx[1], decayVtx[2]);

    double tV0mom[3];
    v0->GetPxPyPz(tV0mom[0], tV0mom[1], tV0mom[2]);
    double lV0TotalMomentum = std::sqrt(
        tV0mom[0] * tV0mom[0] + tV0mom[1] * tV0mom[1] + tV0mom[2] * tV0mom[2]);

    double v0Radius = std::hypot(decayVtx[0], decayVtx[1]);

    v0Pt = v0->Pt();
    if ((v0Pt < fMinV0Pt) || (fMaxV0Pt < v0Pt)) continue;

    int lKeyPos = std::abs(v0->GetPindex());
    int lKeyNeg = std::abs(v0->GetNindex());

    double momPos[3];
    v0->GetPPxPyPz(momPos[0], momPos[1], momPos[2]);
    double momNeg[3];
    v0->GetNPxPyPz(momNeg[0], momNeg[1], momNeg[2]);

    // Calculate the sign of the vec prod with momenta projected to xy plane
    // It is unnecessary to to the full calculation like done in the original
    // task
    double lVecProd = momPos[0] * momNeg[1] - momPos[1] * momNeg[0];
    bool isCowboy = lVecProd * magneticField < 0;

    AliESDtrack *pTrack = esdEvent->GetTrack(lKeyPos);
    AliESDtrack *nTrack = esdEvent->GetTrack(lKeyNeg);

    if (!pTrack || !nTrack) {
      ::Fatal("AliAnalysisTaskStrangenessLifetimes::UserExec",
              "Could not retreive one of the daughter track");
      continue;
    }

    /// TODO: check if this extra cleanup is required
    if (std::abs(nTrack->Eta()) > 0.8 || std::abs(pTrack->Eta()) > 0.8)
      continue;
    if (std::abs(v0->RapK0Short()) > 0.5 && std::abs(v0->RapLambda()) > 0.5)
      continue;

    // Filter like-sign V0 (next: add counter and distribution)
    if (pTrack->GetSign() == nTrack->GetSign()) {
      continue;
    }

    // Track quality cuts
    unsigned char posXedRows = pTrack->GetTPCClusterInfo(2, 1);
    unsigned char negXedRows = nTrack->GetTPCClusterInfo(2, 1);

    // TPC refit condition (done during reconstruction for Offline but not for
    // On-the-fly)
    if (!(pTrack->GetStatus() & AliESDtrack::kTPCrefit)) continue;
    if (!(nTrack->GetStatus() & AliESDtrack::kTPCrefit)) continue;

    float negB[2], posB[2], bCov[3];
    pTrack->GetImpactParameters(posB, bCov);
    nTrack->GetImpactParameters(negB, bCov);

    // GetKinkIndex condition
    if (pTrack->GetKinkIndex(0) > 0 || nTrack->GetKinkIndex(0) > 0) continue;

    // Findable cluster s > 0 condition
    if (pTrack->GetTPCNclsF() <= 0 || nTrack->GetTPCNclsF() <= 0) continue;

    float posTrackXedRowsOverFindable =
        float(posXedRows) / pTrack->GetTPCNclsF();
    float negTrackXedRowsOverFindable =
        float(negXedRows) / nTrack->GetTPCNclsF();

    float posChi2PerCluster =
        pTrack->GetTPCchi2() / (pTrack->GetTPCNcls() + 1.e-16);
    float negChi2PerCluster =
        nTrack->GetTPCchi2() / (nTrack->GetTPCNcls() + 1.e-16);

    // Extra track quality: min track length
    float posTrackLength = -1;
    float negTrackLength = -1;
    if (pTrack->GetInnerParam())
      posTrackLength = pTrack->GetLengthInActiveZone(
          1, 2.0, 220.0, esdEvent->GetMagneticField());
    if (nTrack->GetInnerParam())
      negTrackLength = nTrack->GetLengthInActiveZone(
          1, 2.0, 220.0, esdEvent->GetMagneticField());

    float smallestTrackLength =
        (posTrackLength < negTrackLength) ? posTrackLength : negTrackLength;
    if ((((pTrack->GetTPCClusterInfo(2, 1)) < 70) ||
         ((nTrack->GetTPCClusterInfo(2, 1)) < 70)) &&
        smallestTrackLength < 80)
      continue;

    dcaPosToPrimVertex = std::abs(
        pTrack->GetD(primaryVertex[0], primaryVertex[1], magneticField));

    dcaNegToPrimVertex = std::abs(
        nTrack->GetD(primaryVertex[0], primaryVertex[1], magneticField));

    lV0CosineOfPointingAngle = v0->GetV0CosineOfPointingAngle(
        primaryVertex[0], primaryVertex[1], primaryVertex[2]);
    if (lV0CosineOfPointingAngle < MiniV0::fgkV0cosPA_f) continue;

    // Getting invariant mass infos directly from ESD
    v0->ChangeMassHypothesis(310);
    double invMassK0s = v0->GetEffMass();
    v0->ChangeMassHypothesis(3122);
    double invMassLambda = v0->GetEffMass();

    // Official means of acquiring N-sigmas
    float nSigmasPosProton =
        fPIDResponse->NumberOfSigmasTPC(pTrack, AliPID::kProton);
    float nSigmasPosPion =
        fPIDResponse->NumberOfSigmasTPC(pTrack, AliPID::kPion);
    float nSigmasNegProton =
        fPIDResponse->NumberOfSigmasTPC(nTrack, AliPID::kProton);
    float nSigmasNegPion =
        fPIDResponse->NumberOfSigmasTPC(nTrack, AliPID::kPion);

    MiniV0 miniV0;
    miniV0.SetV0pt(v0Pt);
    miniV0.SetV0eta(v0->Eta());
    miniV0.SetLeastNumberOfXedRows(posXedRows < negXedRows ? posXedRows
                                                           : negXedRows);
    miniV0.SetDistOverP(std::sqrt(Sq(decayVtx[0] - primaryVertex[0]) +
                                  Sq(decayVtx[1] - primaryVertex[1]) +
                                  Sq(decayVtx[2] - primaryVertex[2])) /
                        (lV0TotalMomentum + 1e-16));  // avoid division by zero
    miniV0.SetInvMasses(invMassK0s, invMassLambda);
    miniV0.SetArmenterosVariables(v0->AlphaV0(), v0->PtArmV0());
    miniV0.SetV0CosPA(lV0CosineOfPointingAngle);
    miniV0.SetV0Chi2andCowBoy(v0->GetChi2V0(), isCowboy);
    miniV0.SetProngsDCA(v0->GetDcaV0Daughters());
    miniV0.SetProngsPvDCA(dcaPosToPrimVertex, dcaNegToPrimVertex);
    miniV0.SetV0radius(v0Radius);
    miniV0.SetLeastXedRowsOverFindable(posTrackXedRowsOverFindable <
                                               negTrackXedRowsOverFindable
                                           ? posTrackXedRowsOverFindable
                                           : negTrackXedRowsOverFindable);
    miniV0.SetMaxChi2perCluster(posChi2PerCluster > negChi2PerCluster
                                    ? posChi2PerCluster
                                    : negChi2PerCluster);

    // Rugh 20-sigma selection band, parametric.
    // K0Short: Enough to parametrize peak broadening with linear function.
    double lUpperLimitK0Short = (5.63707e-01) + (1.14979e-02) * v0Pt;
    double lLowerLimitK0Short = (4.30006e-01) - (1.10029e-02) * v0Pt;
    // Lambda: Linear (for higher pt) plus exponential (for low-pt broadening)
    //[0]+[1]*x+[2]*TMath::Exp(-[3]*x)
    double lUpperLimitLambda =
        (1.13688e+00) + (5.27838e-03) * v0Pt +
        (8.42220e-02) * TMath::Exp(-(3.80595e+00) * v0Pt);
    double lLowerLimitLambda =
        (1.09501e+00) - (5.23272e-03) * v0Pt -
        (7.52690e-02) * TMath::Exp(-(3.46339e+00) * v0Pt);
    // Do Selection
    if (
        // Case 1: Lambda Selection
        (miniV0.GetLambdaInvMass() < lUpperLimitLambda &&
         miniV0.GetLambdaInvMass() > lLowerLimitLambda &&
         ((std::abs(nSigmasPosProton) < 6.0 &&
           std::abs(nSigmasNegPion) < 6.0) ||
          (std::abs(nSigmasNegProton) < 6.0 &&
           std::abs(nSigmasPosPion) < 6.0))) ||
        // Case 2: K0Short Selection
        (miniV0.GetK0sInvMass() < lUpperLimitK0Short &&
         miniV0.GetK0sInvMass() > lLowerLimitK0Short &&
         (std::abs(nSigmasNegPion) < 6.0 && std::abs(nSigmasPosPion) < 6.0))) {
      // pT window
      if (v0Pt < fMinPtToSave || v0Pt > fMaxPtToSave) continue;
      miniV0.SetProngsTPCnsigmas(nSigmasPosPion, nSigmasPosProton,
                                 nSigmasNegPion, nSigmasNegProton);
      fV0vector.push_back(miniV0);
    }
  }

  if (fV0vector.size()) fTreeV0->Fill();

  PostData(1, fListHist);
  PostData(2, fTreeV0);
}

void AliAnalysisTaskStrangenessLifetimes::Terminate(Option_t *) {}

void AliAnalysisTaskStrangenessLifetimes::SetupStandardVertexing()
// Meant to store standard re-vertexing configuration
{
  // Tell the task to re-run vertexers
  SetDoV0Refit(true);

  // V0-Related topological selections
  SetV0VertexerDCAFirstToPV(0.05);
  SetV0VertexerDCASecondtoPV(0.05);
  SetV0VertexerDCAV0Daughters(1.20);
  SetV0VertexerCosinePA(0.98);
  SetV0VertexerMinRadius(0.9);
  SetV0VertexerMaxRadius(200);
}

void AliAnalysisTaskStrangenessLifetimes::SetupLooseVertexing()
// Meant to store standard re-vertexing configuration
{
  // Tell the task to re-run vertexers
  SetDoV0Refit(true);

  // V0-Related topological selections
  SetV0VertexerDCAFirstToPV(0.1);
  SetV0VertexerDCASecondtoPV(0.1);
  SetV0VertexerDCAV0Daughters(1.40);
  SetV0VertexerCosinePA(0.95);
  SetV0VertexerMinRadius(0.9);
  SetV0VertexerMaxRadius(200);
}
