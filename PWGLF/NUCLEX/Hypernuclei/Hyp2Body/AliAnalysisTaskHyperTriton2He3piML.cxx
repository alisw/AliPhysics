#include "AliAnalysisTaskHyperTriton2He3piML.h"

#include <array>
#include <unordered_map>
#include <TRandom3.h>
#include <Riostream.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TList.h>
#include <TMath.h>
#include <stdio.h>
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliVTrack.h"
#include "AliExternalTrackParam.h"
#include "AliInputEventHandler.h"
#include "AliMCEvent.h"
#include "AliMultiplicity.h"
#include "AliPID.h"
#include "AliPDG.h"
#include "AliPIDResponse.h"
#include "AliVVertex.h"
#include "AliVertexerHyperTriton2Body.h"
#include "AliAnalysisDataContainer.h"
#include <TFile.h>
#include <TSystem.h>
#include <assert.h>
#include "AliNanoAODTrack.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskHyperTriton2He3piML);

namespace
{

bool ComputeMother(AliMCEvent *mcEvent, const AliVTrack *one, const AliVTrack *two, int &label)
{
  int labOne = std::abs(one->GetLabel());
  int labTwo = std::abs(two->GetLabel());

  if (mcEvent->IsPhysicalPrimary(labOne) ||
      mcEvent->IsPhysicalPrimary(labTwo) ||
      mcEvent->IsSecondaryFromMaterial(labOne) ||
      mcEvent->IsSecondaryFromMaterial(labTwo))
    return false;
  else
  {
    AliVParticle *partOne = mcEvent->GetTrack(labOne);
    AliVParticle *partTwo = mcEvent->GetTrack(labTwo);
    if (partOne->GetMother() != partTwo->GetMother())
      return false;
    else
    {
      int coef = one->GetLabel() < 0 || two->GetLabel() < 0 ? -1 : 1;
      label = coef * partTwo->GetMother();
      return true;
    }
  }
}

std::string ImportSplinePath(std::string path)
{
  std::string modelname = path.substr(path.find_last_of("/") + 1);

  std::string newpath = gSystem->pwd() + std::string("/") + modelname.data();
  std::string oldpath = gDirectory->GetPath();

  bool cpStatus = TFile::Cp(path.data(), newpath.data());
  assert(cpStatus == true);

  gDirectory->Cd(oldpath.data());
  return newpath;
}

bool HasTOF(AliVTrack *track)
{
  const bool hasTOFout = track->GetStatus() & AliVTrack::kTOFout;
  const bool hasTOFtime = track->GetStatus() & AliVTrack::kTIME;
  const float len = track->GetIntegratedLength();
  return hasTOFout && hasTOFtime && (len > 350.);
}

} // namespace

AliAnalysisTaskHyperTriton2He3piML::AliAnalysisTaskHyperTriton2He3piML(
    bool mc, std::string name)
    : AliAnalysisTaskSE(name.data()),
      fEventCuts{},
      fCentralityEstimator{0},
      fFillGenericV0s{true},
      fFillGenericTracklets{false},
      fFillTracklet{true},
      fStoreAllEvents{true},
      fSaveFileNames{false},
      fPropagetToPV{true},
      fV0Vertexer{},
      fLambda{false},
      fListHist{nullptr},
      fTreeV0{nullptr},
      fInputHandler{nullptr},
      fPIDResponse{nullptr},
      fCVMFSPath{""},
      fMC{false},
      fUseOnTheFly{false},
      fUseNanoAODs{false},
      fUseCustomBethe{false},
      fCustomBethe{0.f, 0.f, 0.f, 0.f, 0.f},
      fCustomResolution{1.f},
      fHistNsigmaHe3{nullptr},
      fHistNsigmaPi{nullptr},
      fHistInvMass{nullptr},
      fHistTPCdEdx{nullptr},
      fHistTrackletThetaPhi{nullptr},
      fHistTrackletDThetaDPhi{nullptr},
      fHistTrackletCosP{nullptr},
      fMinPtToSave{0.1},
      fMaxPtToSave{100},
      fMaxTPCpiSigma{10.},
      fMaxTPChe3Sigma{10.},
      fMinHe3pt{0.},
      fMinTPCclusters{50},
      fMinPIDclusters{30},
      fMaxDeltaPhi{0.12},
      fMaxDeltaTheta{0.12},
      fMinTrackletCosP{0.8},
      fEnableLikeSign{false},
      fCurrentFileName{""},
      fSHyperTriton{},
      fSGenericV0{},
      fRHyperTriton{},
      fRTracklets{},
      fSGenericTracklets{},
      fRCollision{},
      fFatParticle{AliPID::kHe3},
      fHyperPDG{1010010030}
{

  // Standard output
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class()); // Basic Histograms
  DefineOutput(2, TTree::Class()); // V0 Tree output
  DefineOutput(3, TTree::Class()); // FileName Tree output
}

AliAnalysisTaskHyperTriton2He3piML::~AliAnalysisTaskHyperTriton2He3piML()
{
  if (fListHist)
  {
    delete fListHist;
    fListHist = nullptr;
  }

  if (fTreeV0)
  {
    delete fTreeV0;
    fTreeV0 = nullptr;
  }

}

void AliAnalysisTaskHyperTriton2He3piML::UserCreateOutputObjects()
{
  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  fInputHandler = (AliInputEventHandler *)(man->GetInputEventHandler());
  fInputHandler->SetNeedField();


  OpenFile(2);
  fTreeV0 = new TTree("fTreeV0", "V0 Candidates");
  fTreeV0->Branch("RCollision", &fRCollision);
  fTreeV0->Branch("RHyperTriton", &fRHyperTriton);
  if (fSaveFileNames) {
    fTreeV0->Branch("Filename", &fCurrentFileName);
    fStoreAllEvents = false;
  }
  if (fFillTracklet)
    fTreeV0->Branch("RTracklets", &fRTracklets);

  if (man->GetMCtruthEventHandler())
  {
    fTreeV0->Branch("SHyperTriton", &fSHyperTriton);
    if (fFillGenericV0s)
      fTreeV0->Branch("SGenericV0", &fSGenericV0);
    if (fFillGenericTracklets)
      fTreeV0->Branch("SGenericTracklets", &fSGenericTracklets);
  }

  if (fCVMFSPath != "")
  {
    std::string LocalPath = ImportSplinePath(fCVMFSPath);
    fV0Vertexer.SetSpline(LocalPath);
  }

  fListHist = new TList();
  fListHist->SetOwner();
  fEventCuts.AddQAplotsToList(fListHist);

  fHistNsigmaPi =
      new TH2D("fHistNsigmaPi", ";#it{p}_{T} (GeV/#it{c});n_{#sigma} TPC Pion; Counts", 100, 0, 10, 80, -5, 5);
  fHistNsigmaHe3 =
      new TH2D("fHistNsigmaHe3", ";#it{p}_{T} (GeV/#it{c});n_{#sigma} TPC ^{3}He; Counts", 100, 0, 10, 80, -5, 5);
  fHistInvMass =
      new TH2D("fHistInvMass", ";#it{p}_{T} (GeV/#it{c});Invariant Mass(GeV/#it{c^2}); Counts", 100, 0, 10, 30, 2.96, 3.05);

  fHistTPCdEdx[0] = new TH2D("fHistTPCdEdxPos", ";#it{p} (GeV/#it{c}); dE/dx; Counts", 256, 0, 10.24, 4096, 0, 2048);
  fHistTPCdEdx[1] = new TH2D("fHistTPCdEdxNeg", ";#it{p} (GeV/#it{c}); dE/dx; Counts", 256, 0, 10.24, 4096, 0, 2048);

  fHistTrackletThetaPhi = new TH2D("fHistTrackletThetaPhi", "; #theta_{trkl}; #phi_{trkl} (rad); Counts", 200, -1., 1., 315., 0.0, 630.);
  fHistTrackletDThetaDPhi = new TH2D("fHistTrackletDThetaDPhi", "; #Delta#theta_{trkl}; #Delta#phi_{trkl} (rad); Counts", 60, -0.3, 0.3, 60, -0.3, 0.3);
  fHistTrackletCosP = new TH1D("fHistTrackletCosP", "; cos #theta_{P,trkl}; Counts", 100, 0.0, 1.0);

  fListHist->Add(fHistNsigmaPi);
  fListHist->Add(fHistNsigmaHe3);
  fListHist->Add(fHistInvMass);
  fListHist->Add(fHistTPCdEdx[0]);
  fListHist->Add(fHistTPCdEdx[1]);

  PostData(1, fListHist);
  PostData(2, fTreeV0);

  AliPDG::AddParticlesToPdgDataBase();

  fFatParticle = fLambda ? AliPID::kProton : AliPID::kHe3;
  fHyperPDG = fLambda ? 3122 : 1010010030;
} // end UserCreateOutputObjects

void AliAnalysisTaskHyperTriton2He3piML::UserExec(Option_t *)
{
  AliVEvent *vEvent = InputEvent();
  AliESDEvent *esdEvent = dynamic_cast<AliESDEvent *>(vEvent);

  if (!vEvent)
  {
    ::Fatal("AliAnalysisTaskHyperTriton2He3piML::UserExec",
            "AliVEvent not found.");
    return;
  }

  
  AliMCEvent *mcEvent = MCEvent();
  if (!mcEvent && fMC)
  {
    ::Fatal("AliAnalysisTaskHyperTriton2He3piML::UserExec", "Could not retrieve MC event");
    return;
  }

  
  if (!fEventCuts.AcceptEvent(vEvent))
  {
    PostData(1, fListHist);
    PostData(2, fTreeV0);
  }

  
  if (fSaveFileNames)
  {
    if (fCurrentFileName.String() != CurrentFileName())
    {
      fCurrentFileName = CurrentFileName();
    }
  }
  
  double primaryVertex[3];
  fRCollision.fCent = fEventCuts.GetCentrality(fCentralityEstimator);
  fEventCuts.GetPrimaryVertex()->GetXYZ(primaryVertex);

  fRCollision.fX = primaryVertex[0];
  fRCollision.fY = primaryVertex[1];
  fRCollision.fZ = primaryVertex[2];

  unsigned char tgr = 0x0;

  if (fInputHandler->IsEventSelected() & AliVEvent::kINT7)
    tgr |= kINT7;
  if (fInputHandler->IsEventSelected() & AliVEvent::kCentral)
    tgr |= kCentral;
  if (fInputHandler->IsEventSelected() & AliVEvent::kSemiCentral)
    tgr |= kSemiCentral;
  int magField = vEvent->GetMagneticField() > 0 ? kPositiveB : 0;

  
  fPIDResponse = fInputHandler->GetPIDResponse();
  fRCollision.fTrigger = tgr + magField;  

  std::unordered_map<int, int> mcMap;
  if (fMC)
  {
    
    fSHyperTriton.clear();
    fSGenericV0.clear();
    fSGenericTracklets.clear();
    for (int ilab = 0; ilab < mcEvent->GetNumberOfTracks(); ilab++)
    { // This is the begining of the loop on tracks
      AliVParticle *part = mcEvent->GetTrack(ilab);
      if (!part)
      {
        ::Warning("AliAnalysisTaskHyperTriton2He3piML::UserExec", "Generated loop %d - MC TParticle pointer to current stack particle = 0x0 ! Skipping.", ilab);
        continue;
      }

      int currentPDG = part->PdgCode();
      if (std::abs(currentPDG) == fHyperPDG)
      {
        if (std::abs(part->Y()) > 1.)
          continue;

        double sVtx[3]{0.0, 0.0, 0.0};

        AliVParticle *he3{nullptr}, *pi{nullptr};
        for (int iD = part->GetDaughterFirst(); iD <= part->GetDaughterLast(); ++iD)
        {
          AliVParticle *dau = mcEvent->GetTrack(iD);

          if (mcEvent->IsSecondaryFromWeakDecay(iD) && dau)
          {
            sVtx[0] = dau->Xv();
            sVtx[1] = dau->Yv();
            sVtx[2] = dau->Zv();

            if (std::abs(dau->PdgCode()) == AliPID::ParticleCode(fFatParticle))
              he3 = dau;
            if (std::abs(dau->PdgCode()) == AliPID::ParticleCode(AliPID::kPion))
              pi = dau;
          }
        }
        if (he3 == nullptr)
          continue;

        SHyperTritonHe3pi v0part;
        v0part.fPdgCode = currentPDG;
        v0part.fDecayX = sVtx[0] - part->Xv();
        v0part.fDecayY = sVtx[1] - part->Yv();
        v0part.fDecayZ = sVtx[2] - part->Zv();
        v0part.fPxHe3 = he3->Px();
        v0part.fPyHe3 = he3->Py();
        v0part.fPzHe3 = he3->Pz();
        v0part.fPxPi = pi->Px();
        v0part.fPyPi = pi->Py();
        v0part.fPzPi = pi->Pz();
        v0part.fFake = true;
        v0part.fRecoIndex = -1;
        v0part.fRecoTracklet = -1;
        v0part.fNegativeLabels = true;
        mcMap[ilab] = fSHyperTriton.size();
        fSHyperTriton.push_back(v0part);
      }
    }
  }

  fRHyperTriton.clear();
  std::vector<int> he3TrackIndices;

  std::vector<AliESDv0> V0Vector;
  if (!fUseOnTheFly && !fUseNanoAODs)
  { 
    if(!fLambda)
      esdEvent->ResetV0s();
    V0Vector = fV0Vertexer.Tracks2V0vertices(esdEvent, fPIDResponse, mcEvent, fLambda);
  }

  int nV0s = (fUseOnTheFly || fUseNanoAODs) ? esdEvent->GetNumberOfV0s() : V0Vector.size();
  
  for (int iV0 = 0; iV0 < nV0s; iV0++)
  { // This is the begining of the V0 loop (we analyse only offline
    // V0s)
    RHyperTritonHe3pi v0part;
    int he3index;
    int lKeyPos;
    int lKeyNeg;

    if (!fUseNanoAODs)
    {
      AliESDv0 *v0 = fUseOnTheFly ? esdEvent->GetV0(iV0) : &V0Vector[iV0];
      lKeyPos = std::abs(v0->GetPindex());
      lKeyNeg = std::abs(v0->GetNindex());
      double pP[3]{0.0}, nP[3]{0.0};
      v0->GetPPxPyPz(pP[0], pP[1], pP[2]);
      v0->GetNPxPyPz(nP[0], nP[1], nP[2]);

      Bool_t isFilled = FillHyperCandidate(v0, vEvent, mcEvent, mcMap, pP, nP, lKeyPos, lKeyNeg, v0part, he3index);
      if (!isFilled)
        continue;

      double x{0.}, y{0.}, z{0.};
      v0->GetXYZ(x, y, z);
      v0part.fDecayX = x - fRCollision.fX;
      v0part.fDecayY = y - fRCollision.fY;
      v0part.fDecayZ = z - fRCollision.fZ;
      v0part.fDcaV0daughters = v0->GetDcaV0Daughters();
      v0part.fChi2V0 = v0->GetChi2V0();
    }
    else
    {
      AliAODv0 *v0 = dynamic_cast<AliAODEvent *>(vEvent)->GetV0(iV0);
      double pP[3] = {v0->MomPosX(), v0->MomPosY(), v0->MomPosZ()};
      double nP[3] = {v0->MomNegX(), v0->MomNegY(), v0->MomNegZ()};
      lKeyPos = std::abs(v0->GetPosID());
      lKeyNeg = std::abs(v0->GetNegID());

      Bool_t isFilled = FillHyperCandidate(v0, vEvent, mcEvent, mcMap, pP, nP, lKeyPos, lKeyNeg, v0part, he3index);
      if (!isFilled)
        continue;
      
      double x{v0->AliAODv0::DecayVertexV0X()}, y{v0->AliAODv0::DecayVertexV0Y()}, z{v0->AliAODv0::DecayVertexV0Z()};
      v0part.fDecayX = x - fRCollision.fX;
      v0part.fDecayY = y - fRCollision.fY;
      v0part.fDecayZ = z - fRCollision.fZ;
      v0part.fDcaV0daughters = v0->AliAODv0::DcaV0Daughters();
      v0part.fChi2V0 = v0->Chi2V0();
    }

    fRHyperTriton.push_back(v0part);
    he3TrackIndices.push_back(he3index);
  }
  // loop on tracklets to match them with mother MClabel1m
  fRTracklets.clear();
  if (!fUseNanoAODs)
  {
    AliMultiplicity *tracklets = esdEvent->GetMultiplicity();
    int nTracklets = tracklets->GetNumberOfTracklets();

    for (size_t iHyper{0}; iHyper < fRHyperTriton.size(); ++iHyper)
    {
      const auto &v0 = fRHyperTriton[iHyper];
      for (int iTracklet = 0; iTracklet < nTracklets; iTracklet++)
      {
        double theta = tracklets->GetTheta(iTracklet);
        double phi = tracklets->GetPhi(iTracklet);
        double deltaTheta = tracklets->GetDeltaTheta(iTracklet);
        double deltaPhi = tracklets->GetDeltaPhi(iTracklet);
        fHistTrackletThetaPhi->Fill(theta, phi);
        fHistTrackletDThetaDPhi->Fill(deltaTheta, deltaPhi);

        int id1{-1}, id2{-1};
        tracklets->GetTrackletTrackIDs(iTracklet, 0, id1, id2); // references for eventual Global/ITS_SA tracks

        if (id1 >= 0 && id2 >= 0 && id1 != he3TrackIndices[iHyper] && id2 != he3TrackIndices[iHyper]) /// Both points are used in a track that is not the candidate He3
          continue;

        if (std::abs(deltaPhi) > fMaxDeltaPhi)
          continue;
        if (std::abs(deltaTheta) > fMaxDeltaTheta)
          continue;

        double cx = std::cos(phi) * std::sin(theta);
        double cy = std::sin(phi) * std::sin(theta);
        double cz = std::cos(theta);

        const double cosp = (v0.fDecayX * cx + v0.fDecayY * cy + v0.fDecayZ * cz) / std::sqrt(v0.fDecayX * v0.fDecayX + v0.fDecayY * v0.fDecayY + v0.fDecayZ * v0.fDecayZ);
        fHistTrackletCosP->Fill(cosp);
        if (cosp > fMinTrackletCosP)
        {
          RTracklet trkl;
          trkl.fDeltaPhi = deltaPhi;
          trkl.fDeltaTheta = deltaTheta;
          trkl.fPhi = phi;
          trkl.fTheta = theta;
          trkl.fSharedCluster = (id1 >= 0) || (id2 >= 0);
          if (tracklets->GetLabel(iTracklet, 0) == tracklets->GetLabel(iTracklet, 1) && fMC && tracklets->GetLabel(iTracklet, 0) >= 0)
          {
            int ilab = tracklets->GetLabel(iTracklet, 0);
            AliVParticle *part = mcEvent->GetTrack(ilab);
            if (std::abs(part->PdgCode()) == fHyperPDG)
              fSHyperTriton[mcMap[ilab]].fRecoTracklet = fRTracklets.size();
            else
            {
              AliVParticle *part = mcEvent->GetTrack(ilab);
              SGenericTracklet gen;
              gen.fPdgCode = part->PdgCode();
              gen.fPx = part->Px();
              gen.fPy = part->Py();
              gen.fPz = part->Pz();
              gen.fRecoIndex = fRTracklets.size();
              fSGenericTracklets.emplace_back(gen);
            }
          }
          fRTracklets.push_back(trkl);
          break;
        }
      }
    }
  }

  if (fRHyperTriton.size() != 0 || fStoreAllEvents)
    fTreeV0->Fill();

  PostData(1, fListHist);
  PostData(2, fTreeV0);
}

void AliAnalysisTaskHyperTriton2He3piML::Terminate(Option_t *) {}

void AliAnalysisTaskHyperTriton2He3piML::SetCustomBetheBloch(float res, const float *bethe)
{
  fUseCustomBethe = true;
  fCustomResolution = res;
  std::copy(bethe, bethe + 5, fCustomBethe);
}

double AliAnalysisTaskHyperTriton2He3piML::customNsigma(double mom, double sig)
{
  const float bg = mom / AliPID::ParticleMass(fFatParticle);
  const float *p = fCustomBethe;
  const float expS = AliExternalTrackParam::BetheBlochAleph(bg, p[0], p[1], p[2], p[3], p[4]);
  return (sig - expS) / (fCustomResolution * expS);
}
template <class T, class M>

Bool_t AliAnalysisTaskHyperTriton2He3piML::FillHyperCandidate(T *v0, AliVEvent *event, AliMCEvent *mcEvent, M mcMap,
                                                              double *pP, double *nP, int lKeyPos, int lKeyNeg, RHyperTritonHe3pi &v0part, int &he3index)
{
  if (!v0)
    return false;
  if (v0->GetOnFlyStatus() != 0 && !fUseOnTheFly)
    return false;
  if (fUseOnTheFly && v0->GetOnFlyStatus() == 0)
    return false;

  AliVTrack *pTrack = dynamic_cast<AliVTrack *>(event->GetTrack(lKeyPos));
  AliVTrack *nTrack = dynamic_cast<AliVTrack *>(event->GetTrack(lKeyNeg));

  // Remove like-sign (will not affect offline V0 candidates!)
  if (fEnableLikeSign)
  {
    if (pTrack->GetSign() * nTrack->GetSign() < 0)
      return false;
  }
  else
  {
    if (pTrack->GetSign() * nTrack->GetSign() > 0)
      return false;
  }

  if (!pTrack || !nTrack)
    ::Fatal("AliAnalysisTaskHyperTriton2He3piML::UserExec",
            "Could not retreive one of the daughter track");

  if (std::abs(nTrack->Eta()) > 0.8 || std::abs(pTrack->Eta()) > 0.8)
    return false;

  // TPC refit condition (done during reconstruction for Offline but not for
  // On-the-fly)
  if (!(pTrack->GetStatus() & AliVTrack::kTPCrefit))
    return false;
  if (!(nTrack->GetStatus() & AliVTrack::kTPCrefit))
    return false;

  // GetKinkIndex condition
  if (pTrack->GetKinkIndex(0) > 0 || nTrack->GetKinkIndex(0) > 0)
    return false;

  // Findable cluster s > 0 condition
  if (pTrack->GetTPCNclsF() <= 0 || nTrack->GetTPCNclsF() <= 0)
    return false;

  if ((pTrack->GetTPCClusterInfo(2, 1) < fMinTPCclusters) ||
      (nTrack->GetTPCClusterInfo(2, 1) < fMinTPCclusters))
    return false;

  if ((pTrack->GetTPCsignalN() < fMinPIDclusters) ||
      (nTrack->GetTPCsignalN() < fMinPIDclusters))
    return false;

  // Official means of acquiring N-sigmas
  float nSigmaPosPi = fPIDResponse->NumberOfSigmasTPC(pTrack, AliPID::kPion);
  float nSigmaPosHe3 = fPIDResponse->NumberOfSigmasTPC(pTrack, fFatParticle);
  float nSigmaNegPi = fPIDResponse->NumberOfSigmasTPC(nTrack, AliPID::kPion);
  float nSigmaNegHe3 = fPIDResponse->NumberOfSigmasTPC(nTrack, fFatParticle);

  if (fUseCustomBethe)
  {
    nSigmaNegHe3 = customNsigma(nTrack->GetTPCmomentum(), nTrack->GetTPCsignal());
    nSigmaPosHe3 = customNsigma(pTrack->GetTPCmomentum(), pTrack->GetTPCsignal());
  }

  const float nSigmaNegAbsHe3 = std::abs(nSigmaNegHe3);
  const float nSigmaPosAbsHe3 = std::abs(nSigmaPosHe3);
  const float nSigmaNegAbsPi = std::abs(nSigmaNegPi);
  const float nSigmaPosAbsPi = std::abs(nSigmaPosPi);

  /// Skip V0 not involving (anti-)He3 candidates

  fHistTPCdEdx[0]->Fill(pTrack->GetTPCmomentum(), pTrack->GetTPCsignal());
  fHistTPCdEdx[1]->Fill(nTrack->GetTPCmomentum(), nTrack->GetTPCsignal());



  AliVTrack *he3Track;
  AliVTrack *piTrack;

  bool mHyperTriton = nSigmaPosAbsHe3 < fMaxTPChe3Sigma && nSigmaNegAbsPi < fMaxTPCpiSigma;
  bool aHyperTriton = nSigmaNegAbsHe3 < fMaxTPChe3Sigma && nSigmaPosAbsPi < fMaxTPCpiSigma;

  if (!mHyperTriton && !aHyperTriton)
    return false;
  if(fLambda){
    mHyperTriton = v0->AlphaV0()>0;
    aHyperTriton = v0->AlphaV0()<0;
    he3Track = aHyperTriton ? nTrack : pTrack;
    piTrack = he3Track == nTrack ? pTrack : nTrack;
  }
  else{

    he3Track = aHyperTriton ? nTrack : pTrack;
    piTrack = he3Track == nTrack ? pTrack : nTrack;
  }

  const double charge = fLambda ? 1. : 2.;
  if (he3Track->Pt() * charge < fMinHe3pt)
    return false;

  const double *he3P = (he3Track == pTrack) ? pP : nP;
  const double *piP = (piTrack == pTrack) ? pP : nP;

  LVector_t he3Vector, piVector, hyperVector;
  he3Vector.SetCoordinates(charge * he3P[0], charge * he3P[1], charge * he3P[2], AliPID::ParticleMass(fFatParticle));
  piVector.SetCoordinates(piP[0], piP[1], piP[2], AliPID::ParticleMass(AliPID::kPion));
  hyperVector = piVector + he3Vector;

  float v0Pt = hyperVector.Pt();
  if ((v0Pt < fMinPtToSave) || (fMaxPtToSave < v0Pt))
    return false;

  float he3B[2], piB[2], bCov[3];
  // if (fPropagetToPV)
  // {
  //   piTrack->PropagateToDCA(fEventCuts.GetPrimaryVertex(), event->GetMagneticField(), 25.);
  //   he3Track->PropagateToDCA(fEventCuts.GetPrimaryVertex(), event->GetMagneticField(), 25.);
  // }
  piTrack->GetImpactParameters(piB, bCov);
  he3Track->GetImpactParameters(he3B, bCov);
  const float he3DCA = std::hypot(he3B[0], he3B[1]);
  const float piDCA = std::hypot(piB[0], piB[1]);

  unsigned char posXedRows = pTrack->GetTPCClusterInfo(2, 1);
  unsigned char negXedRows = nTrack->GetTPCClusterInfo(2, 1);
  float posChi2PerCluster =
      pTrack->GetTPCchi2() / (pTrack->GetTPCNcls() + 1.e-16);
  float negChi2PerCluster =
      nTrack->GetTPCchi2() / (nTrack->GetTPCNcls() + 1.e-16);
  float posXedRowsOverFindable = float(posXedRows) / pTrack->GetTPCNclsF();
  float negXedRowsOverFindable = float(negXedRows) / nTrack->GetTPCNclsF();
  float minXedRowsOverFindable =
      posXedRowsOverFindable < negXedRowsOverFindable
          ? posXedRowsOverFindable
          : negXedRowsOverFindable;
  float maxChi2PerCluster = posChi2PerCluster > negChi2PerCluster
                                ? posChi2PerCluster
                                : negChi2PerCluster;

  // Filling the V0 vector
  int ilab = -1;
  bool isFake = true;
  if (fMC && !fUseNanoAODs)
  {
    int label = 0;
    if (ComputeMother(mcEvent, he3Track, piTrack, label))
    {
      ilab = std::abs(label);
      AliVParticle *part = mcEvent->GetTrack(ilab);
      if (part)
      {
        if (std::abs(part->PdgCode()) == fHyperPDG)
        {
          fSHyperTriton[mcMap[ilab]].fRecoIndex = (fRHyperTriton.size());
          fSHyperTriton[mcMap[ilab]].fFake = false;
          fSHyperTriton[mcMap[ilab]].fNegativeLabels = (label < 0);
          isFake = false;
        }
        else
        {
          SGenericV0 genV0;
          genV0.fRecoIndex = fRHyperTriton.size();
          genV0.fPdgCode = part->PdgCode();
          double sVtx[3]{0.0, 0.0, 0.0};
          for (int iD = part->GetDaughterFirst(); iD <= part->GetDaughterLast(); ++iD)
          {
            AliVParticle *dau = mcEvent->GetTrack(iD);
            if (mcEvent->IsSecondaryFromWeakDecay(iD) && dau)
            {
              sVtx[0] = dau->Xv();
              sVtx[1] = dau->Yv();
              sVtx[2] = dau->Zv();
              break;
            }
          }
          genV0.fDecayX = sVtx[0] - part->Xv();
          genV0.fDecayY = sVtx[1] - part->Yv();
          genV0.fDecayZ = sVtx[2] - part->Zv();
          genV0.fPx = part->Px();
          genV0.fPy = part->Py();
          genV0.fPz = part->Pz();
          fSGenericV0.push_back(genV0);
        }
      }
    }
  }
  if (isFake && fV0Vertexer.GetMonteCarloStatus())
    return false;

  v0part.fPxHe3 = he3Vector.Px();
  v0part.fPyHe3 = he3Vector.Py();
  v0part.fPzHe3 = he3Vector.Pz();
  v0part.fPxPi = piVector.Px();
  v0part.fPyPi = piVector.Py();
  v0part.fPzPi = piVector.Pz();
  v0part.fTPCmomHe3 = he3Track->GetTPCmomentum();
  v0part.fTPCmomPi = piTrack->GetTPCmomentum();
  v0part.fDcaHe32PrimaryVertexXY = std::abs(he3B[0]);
  v0part.fDcaPi2PrimaryVertexXY =
   std::abs(piB[0]);
  v0part.fDcaHe32PrimaryVertex = he3DCA;
  v0part.fDcaPi2PrimaryVertex = piDCA;
  v0part.fLeastXedOverFindable = minXedRowsOverFindable;
  v0part.fMaxChi2PerCluster = maxChi2PerCluster;
  v0part.fTPCnSigmaHe3 = (pTrack == he3Track) ? nSigmaPosHe3 : nSigmaNegHe3;
  v0part.fTPCnSigmaPi = (pTrack == piTrack) ? nSigmaPosPi : nSigmaNegPi;
  v0part.fTOFnSigmaHe3 = fPIDResponse->NumberOfSigmasTOF(he3Track, fFatParticle);
  v0part.fTOFnSigmaPi = fPIDResponse->NumberOfSigmasTOF(piTrack, AliPID::kPion);
  v0part.fNpidClustersHe3 = he3Track->GetTPCsignalN();
  v0part.fNpidClustersPi = piTrack->GetTPCsignalN();
  v0part.fTPCsignalHe3 = he3Track->GetTPCsignal();
  v0part.fTPCsignalPi = piTrack->GetTPCsignal();
  v0part.fITSrefitHe3 = he3Track->GetStatus() & AliVTrack::kITSrefit;
  v0part.fITSrefitPi = piTrack->GetStatus() & AliVTrack::kITSrefit;
  v0part.fITSclusHe3 = he3Track->GetITSClusterMap();
  v0part.fITSclusPi = piTrack->GetITSClusterMap();
  v0part.fTOFmatchHe3 = HasTOF(he3Track);
  v0part.fTOFmatchPi = HasTOF(piTrack);
  v0part.fMatter = (pTrack == he3Track);

  fHistNsigmaPi->Fill(piTrack->Pt(), v0part.fTPCnSigmaPi);
  fHistNsigmaHe3->Fill(he3Vector.Pt(), v0part.fTPCnSigmaHe3);
  fHistInvMass->Fill(hyperVector.Pt(), hyperVector.M());

  he3index = aHyperTriton ? lKeyNeg : lKeyPos;
  return true;
}

AliAnalysisTaskHyperTriton2He3piML *AliAnalysisTaskHyperTriton2He3piML::AddTask(bool isMC, TString suffix)
{
  // Get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskHyperTriton2BodyML", "No analysis manager found.");
    return nullptr;
  }

  // Check the analysis type using the event handlers connected to the analysis
  // manager.
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskHyperTriton2BodyML", "This task requires an input event handler");
    return nullptr;
  }

  TString tskname = "AliAnalysisTaskHyperTriton2He3piML";
  tskname.Append(suffix.Data());
  AliAnalysisTaskHyperTriton2He3piML *task = new AliAnalysisTaskHyperTriton2He3piML(isMC, tskname.Data());

  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(
      Form("%s_summary", tskname.Data()), TList::Class(),
      AliAnalysisManager::kOutputContainer, "AnalysisResults.root");

  AliAnalysisDataContainer *coutput2 =
      mgr->CreateContainer(Form("HyperTritonTree%s", suffix.Data()), TTree::Class(),
                           AliAnalysisManager::kOutputContainer, Form("HyperTritonTree.root:%s", suffix.Data()));
  coutput2->SetSpecialOutput();

  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, coutput1);
  mgr->ConnectOutput(task, 2, coutput2);
  return task;
}
