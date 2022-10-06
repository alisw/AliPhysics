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
#include "AliDataFile.h"
#include "AliESDEvent.h"
#include "AliESDtrackCuts.h"
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

#include "AliOADBContainer.h"

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
      fMaxInfo{false},
      fCentralityEstimator{0},
      fFillGenericV0s{true},
      fFillGenericTracklets{false},
      fFillTracklet{true},
      fStoreAllEvents{true},
      fSaveFileNames{false},
      fPropagetToPV{true},
      fV0Vertexer{},
      fLambda{false},
      fUseTPCmomentum{false},
      fNHarm{2},
      fNucleus{AliPID::kHe3},
      fV0CalibrationFile{""},
      fListHist{nullptr},
      fTreeV0{nullptr},
      fInputHandler{nullptr},
      fPIDResponse{nullptr},
      fCVMFSPath{""},
      fMC{mc},
      fUseOnTheFly{false},
      fUseNanoAODs{false},
      fUseCustomBethe{false},
      fCustomBethe{0.f, 0.f, 0.f, 0.f, 0.f},
      fCustomResolution{1.f},
      fHistCentTrigger{nullptr},
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
      fEnableEventMixing{false},
      fCurrentFileName{""},
      fCurrentEventNumber{-1},
      fSHyperTriton{},
      fSGenericV0{},
      fRHyperTriton{},
      fRHyperTritonFull{},
      fRTracklets{},
      fSGenericTracklets{},
      fRCollision{},
      fRPVcovariance{},
      fFatParticle{AliPID::kHe3},
      fHyperPDG{1010010030},
      fEMdepth{10},
      fHe3mixed{},
      fMultV0{nullptr},
      fQxnmV0A{nullptr},
      fQynmV0A{nullptr},
      fQxnsV0A{nullptr},
      fQynsV0A{nullptr},
      fQxnmV0C{nullptr},
      fQynmV0C{nullptr},
      fQxnsV0C{nullptr},
      fQynsV0C{nullptr},
      EPVzAvsCentrality{nullptr},
      EPVzCvsCentrality{nullptr},
      hQVzAQVzCvsCentrality{nullptr},
      hQVzAQTPCvsCentrality{nullptr},
      hQVzCQTPCvsCentrality{nullptr},
      hQxVzAvsCentrality{nullptr},
      hQyVzAvsCentrality{nullptr},
      hQxVzCvsCentrality{nullptr},
      hQyVzCvsCentrality{nullptr},
      hCos2DeltaTPCVzAvsCentrality{nullptr},
      hCos2DeltaTPCVzCvsCentrality{nullptr},
      hCos2DeltaVzAVzCvsCentrality{nullptr},
      hCos2DeltaVzATPCvsCentrality{nullptr},
      hCos2DeltaVzCTPCvsCentrality{nullptr},
      hCos2DeltaVzCVzAvsCentrality{nullptr},
      fESDtrackCutsEP{nullptr}
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
  if (fSaveFileNames) {
    fTreeV0->Branch("Filename", &fCurrentFileName);
    fTreeV0->Branch("EventNumber", &fCurrentEventNumber);
    fStoreAllEvents = false;
  }
  if (fMaxInfo) {
    fTreeV0->Branch("RPVcovariance", &fRPVcovariance, "RPVcovariance[6]/F");
    fTreeV0->Branch("RHyperTriton", &fRHyperTritonFull);
  } else {
    fTreeV0->Branch("RHyperTriton", &fRHyperTriton);
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

  fHistCentTrigger =
      new TH2D("fHistCentTrigger", ";Centrality;trigger", 1000, 0, 100, 32, 0, 32);

  fHistNsigmaPi =
      new TH2D("fHistNsigmaPi", ";#it{p}_{T} (GeV/#it{c});n_{#sigma} TPC Pion; Counts", 100, 0, 10, 80, -5, 5);
  fHistNsigmaHe3 =
      new TH2D("fHistNsigmaHe3", ";#it{p}_{T} (GeV/#it{c});n_{#sigma} TPC ^{3}He; Counts", 100, 0, 10, 80, -5, 5);
  fHistInvMass =
      new TH2D("fHistInvMass", ";#it{p}_{T} (GeV/#it{c});Invariant Mass(GeV/#it{c}^{2}); Counts", 100, 0, 10, 90, 2.96, 3.05);

  fHistTPCdEdx[0] = new TH2D("fHistTPCdEdxPos", ";#it{p} (GeV/#it{c}); dE/dx; Counts", 256, 0, 10.24, 4096, 0, 2048);
  fHistTPCdEdx[1] = new TH2D("fHistTPCdEdxNeg", ";#it{p} (GeV/#it{c}); dE/dx; Counts", 256, 0, 10.24, 4096, 0, 2048);

  fHistTrackletThetaPhi = new TH2D("fHistTrackletThetaPhi", "; #theta_{trkl}; #phi_{trkl} (rad); Counts", 200, -1., 1., 315., 0.0, 630.);
  fHistTrackletDThetaDPhi = new TH2D("fHistTrackletDThetaDPhi", "; #Delta#theta_{trkl}; #Delta#phi_{trkl} (rad); Counts", 60, -0.3, 0.3, 60, -0.3, 0.3);
  fHistTrackletCosP = new TH1D("fHistTrackletCosP", "; cos #theta_{P,trkl}; Counts", 100, 0.0, 1.0);

  fListHist->Add(fHistCentTrigger);
  fListHist->Add(fHistNsigmaPi);
  fListHist->Add(fHistNsigmaHe3);
  fListHist->Add(fHistInvMass);
  fListHist->Add(fHistTPCdEdx[0]);
  fListHist->Add(fHistTPCdEdx[1]);

  EPVzAvsCentrality = new TH2D("EPVzAvsCentrality", "EPVzAvsCentrality", 80, -TMath::Pi(), TMath::Pi(), 105, 0, 105);
  EPVzCvsCentrality = new TH2D("EPVzCvsCentrality", "EPVzCvsCentrality", 80, -TMath::Pi(), TMath::Pi(), 105, 0, 105);

  fListHist->Add(EPVzAvsCentrality);
  fListHist->Add(EPVzCvsCentrality);

  hQVzAQVzCvsCentrality = new TH2D("hQVzAQVzCvsCentrality", "hQVzAQVzCvsCentrality", 1000, -100, 100, 105, 0, 105);
  hQVzAQTPCvsCentrality = new TH2D("hQVzAQTPCvsCentrality", "hQVzAQTPCvsCentrality", 1000, -100, 100, 105, 0, 105);
  hQVzCQTPCvsCentrality = new TH2D("hQVzCQTPCvsCentrality", "hQVzCQTPCvsCentrality", 1000, -100, 100, 105, 0, 105);

  fListHist->Add(hQVzAQVzCvsCentrality);
  fListHist->Add(hQVzAQTPCvsCentrality);
  fListHist->Add(hQVzCQTPCvsCentrality);

  hQxVzAvsCentrality = new TH2D("hQxVzAvsCentrality", "hQxVzAvsCentrality", 100, -20, 20, 105, 0, 105);
  hQyVzAvsCentrality = new TH2D("hQyVzAvsCentrality", "hQyVzAvsCentrality", 100, -20, 20, 105, 0, 105);
  hQxVzCvsCentrality = new TH2D("hQxVzCvsCentrality", "hQxVzCvsCentrality", 100, -20, 20, 105, 0, 105);
  hQyVzCvsCentrality = new TH2D("hQyVzCvsCentrality", "hQyVzCvsCentrality", 100, -20, 20, 105, 0, 105);

  fListHist->Add(hQxVzAvsCentrality);
  fListHist->Add(hQyVzAvsCentrality);
  fListHist->Add(hQxVzCvsCentrality);
  fListHist->Add(hQyVzCvsCentrality);

  hCos2DeltaTPCVzAvsCentrality = new TH2D("hCos2DeltaTPCVzAvsCentrality", "hCos2DeltaTPCVzAvsCentrality", 100, -1.1, 1.1, 105, 0, 105);
  hCos2DeltaTPCVzCvsCentrality = new TH2D("hCos2DeltaTPCVzCvsCentrality", "hCos2DeltaTPCVzCvsCentrality", 100, -1.1, 1.1, 105, 0, 105);
  hCos2DeltaVzAVzCvsCentrality = new TH2D("hCos2DeltaVzAVzCvsCentrality", "hCos2DeltaVzAVzCvsCentrality", 100, -1.1, 1.1, 105, 0, 105);
  hCos2DeltaVzATPCvsCentrality = new TH2D("hCos2DeltaVzATPCvsCentrality", "hCos2DeltaVzATPCvsCentrality", 100, -1.1, 1.1, 105, 0, 105);
  hCos2DeltaVzCTPCvsCentrality = new TH2D("hCos2DeltaVzCTPCvsCentrality", "hCos2DeltaVzCTPCvsCentrality", 100, -1.1, 1.1, 105, 0, 105);
  hCos2DeltaVzCVzAvsCentrality = new TH2D("hCos2DeltaVzCVzAvsCentrality", "hCos2DeltaVzCVzAvsCentrality", 100, -1.1, 1.1, 105, 0, 105);

  fListHist->Add(hCos2DeltaTPCVzAvsCentrality);
  fListHist->Add(hCos2DeltaTPCVzCvsCentrality);
  fListHist->Add(hCos2DeltaVzAVzCvsCentrality);
  fListHist->Add(hCos2DeltaVzATPCvsCentrality);
  fListHist->Add(hCos2DeltaVzCTPCvsCentrality);
  fListHist->Add(hCos2DeltaVzCVzAvsCentrality);

  fESDtrackCutsEP = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();

  PostData(1, fListHist);
  PostData(2, fTreeV0);

  AliPDG::AddParticlesToPdgDataBase();

  fFatParticle = fLambda ? AliPID::kProton : fNucleus;
  fHyperPDG = fLambda ? 3122 : (fNucleus == AliPID::kHe3 ? 1010010030 : 1010010040);
  fV0Vertexer.fNucleus = fNucleus;
} // end UserCreateOutputObjects

void AliAnalysisTaskHyperTriton2He3piML::UserExec(Option_t *)
{
  AliVEvent *vEvent = InputEvent();
  AliESDEvent *esdEvent = dynamic_cast<AliESDEvent *>(vEvent);

  if (fSaveFileNames)
    fCurrentEventNumber = esdEvent->GetHeader()->GetEventNumberInFile();

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
    return;
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

  auto pv = fEventCuts.GetPrimaryVertex();
  std::array<double,6> pvCov;
  pv->GetCovarianceMatrix(pvCov.data());
  std::copy(pvCov.begin(), pvCov.end(), fRPVcovariance);

  unsigned char tgr = 0x0;

  if (fInputHandler->IsEventSelected() & AliVEvent::kINT7)
    tgr |= kINT7;
  if (fInputHandler->IsEventSelected() & AliVEvent::kCentral)
    tgr |= kCentral;
  if (fInputHandler->IsEventSelected() & AliVEvent::kSemiCentral)
    tgr |= kSemiCentral;
  if (fInputHandler->IsEventSelected() & AliVEvent::kHighMultV0)
    tgr |= kHighMultV0;
  int magField = vEvent->GetMagneticField() > 0 ? kPositiveB : 0;


  fPIDResponse = fInputHandler->GetPIDResponse();
  fRCollision.fTrigger = tgr + magField;

  fHistCentTrigger->Fill(fRCollision.fCent, fRCollision.fTrigger);

  /// For the event mixing
  int centBin = std::floor((fRCollision.fCent - 1.e4)/ 10);
  int zBin = std::floor((fRCollision.fZ + 10 - 1.e4) / 2);
  if (fEnableEventMixing && (zBin < 0 || zBin > 9 || centBin < 0 || centBin > 9))
    AliFatal("Event mixing mode cannot work with events with z vertices outside -10, 10 and centralities not in 0,100");


  if (!fUseNanoAODs)
  {

    AliESDVZERO *esdV0 = (AliESDVZERO *)esdEvent->GetVZEROData();

    int iCen = fRCollision.fCent;

    //V0 info
    double Qxan = 0, Qyan = 0;
    double Qxcn = 0, Qycn = 0;
    double sumMa = 0, sumMc = 0;
    double evPlAngV0A = 0, evPlAngV0C = 0;
    double QxanCor = 0, QyanCor = 0;
    double QxcnCor = 0, QycnCor = 0;
    if (!fV0CalibrationFile.empty()) {

      static bool openCalibs{false};
      if (!openCalibs) {
        OpenInfoCalibration(esdEvent->GetRunNumber());
        openCalibs = true;
      }

      for (int iV0 = 0; iV0 < 64; iV0++)
      {

        double phiV0 = TMath::PiOver4() * (0.5 + iV0 % 8);
        float multv0 = esdV0->GetMultiplicity(iV0);

        if (iV0 < 32)
        {

          double multCorC = -10;

          if (iV0 < 8)
            multCorC = multv0 / fMultV0->GetBinContent(iV0 + 1) * fMultV0->GetBinContent(1);
          else if (iV0 >= 8 && iV0 < 16)
            multCorC = multv0 / fMultV0->GetBinContent(iV0 + 1) * fMultV0->GetBinContent(9);
          else if (iV0 >= 16 && iV0 < 24)
            multCorC = multv0 / fMultV0->GetBinContent(iV0 + 1) * fMultV0->GetBinContent(17);
          else if (iV0 >= 24 && iV0 < 32)
            multCorC = multv0 / fMultV0->GetBinContent(iV0 + 1) * fMultV0->GetBinContent(25);

          if (multCorC < 0)
          {
            cout << "Problem with multiplicity in V0C" << endl;
            continue;
          }

          Qxcn += TMath::Cos(fNHarm * phiV0) * multCorC;
          Qycn += TMath::Sin(fNHarm * phiV0) * multCorC;

          sumMc = sumMc + multCorC;
        }
        else
        {

          double multCorA = -10;

          if (iV0 >= 32 && iV0 < 40)
            multCorA = multv0 / fMultV0->GetBinContent(iV0 + 1) * fMultV0->GetBinContent(33);
          else if (iV0 >= 40 && iV0 < 48)
            multCorA = multv0 / fMultV0->GetBinContent(iV0 + 1) * fMultV0->GetBinContent(41);
          else if (iV0 >= 48 && iV0 < 56)
            multCorA = multv0 / fMultV0->GetBinContent(iV0 + 1) * fMultV0->GetBinContent(49);
          else if (iV0 >= 56 && iV0 < 64)
            multCorA = multv0 / fMultV0->GetBinContent(iV0 + 1) * fMultV0->GetBinContent(57);

          if (multCorA < 0)
          {
            cout << "Problem with multiplicity in V0A" << endl;
            continue;
          }

          Qxan += TMath::Cos(fNHarm * phiV0) * multCorA;
          Qyan += TMath::Sin(fNHarm * phiV0) * multCorA;

          sumMa = sumMa + multCorA;
        }
      }
      //if (sumMa < 0 || sumMc < 0)
      //return;
      QxanCor = Qxan;
      QyanCor = (Qyan - fQynmV0A->GetBinContent(iCen + 1)) / fQynsV0A->GetBinContent(iCen + 1);
      QxcnCor = Qxcn;
      QycnCor = (Qycn - fQynmV0C->GetBinContent(iCen + 1)) / fQynsV0C->GetBinContent(iCen + 1);

      if (fNHarm != 4)
      {
        QxanCor = (Qxan - fQxnmV0A->GetBinContent(iCen + 1)) / fQxnsV0A->GetBinContent(iCen + 1);
        QxcnCor = (Qxcn - fQxnmV0C->GetBinContent(iCen + 1)) / fQxnsV0C->GetBinContent(iCen + 1);
      }

      evPlAngV0A = TMath::ATan2(QyanCor, QxanCor) / fNHarm;
      evPlAngV0C = TMath::ATan2(QycnCor, QxcnCor) / fNHarm;

      EPVzAvsCentrality->Fill(evPlAngV0A, iCen);
      EPVzCvsCentrality->Fill(evPlAngV0C, iCen);
    }

    const int nTracks = esdEvent->GetNumberOfTracks();
    double Qxtn = 0, Qytn = 0;

    for (int it1 = 0; it1 < nTracks; it1++)
    {
      AliESDtrack *esdTrk1 = (AliESDtrack *)esdEvent->GetTrack(it1);

      if (!esdTrk1)
        continue;

      if (!fESDtrackCutsEP->AcceptTrack((AliESDtrack *)esdTrk1))
        continue; //tpc only track

      if (TMath::Abs(esdTrk1->Eta()) < 0.8 && esdTrk1->GetTPCNcls() >= 70 && esdTrk1->Pt() >= 0.2 && esdTrk1->Pt() < 3.)
      {
        Qxtn += TMath::Cos(fNHarm * esdTrk1->Phi());
        Qytn += TMath::Sin(fNHarm * esdTrk1->Phi());
      }
    }

    double evPlAngTPC = TMath::ATan2(Qytn, Qxtn) / fNHarm;

    hCos2DeltaTPCVzAvsCentrality->Fill(TMath::Cos(fNHarm * (evPlAngTPC - evPlAngV0A)), iCen);
    hCos2DeltaTPCVzCvsCentrality->Fill(TMath::Cos(fNHarm * (evPlAngTPC - evPlAngV0C)), iCen);
    hCos2DeltaVzAVzCvsCentrality->Fill(TMath::Cos(fNHarm * (evPlAngV0A - evPlAngV0C)), iCen);
    hCos2DeltaVzATPCvsCentrality->Fill(TMath::Cos(fNHarm * (evPlAngV0A - evPlAngTPC)), iCen);
    hCos2DeltaVzCTPCvsCentrality->Fill(TMath::Cos(fNHarm * (evPlAngV0C - evPlAngTPC)), iCen);
    hCos2DeltaVzCVzAvsCentrality->Fill(TMath::Cos(fNHarm * (evPlAngV0C - evPlAngV0A)), iCen);

    //Scalar Product -- Resolutions

    double corV0AV0Cvn = QxanCor * QxcnCor + QyanCor * QycnCor;
    double corV0ATPCvn = QxanCor * Qxtn + QyanCor * Qytn;
    double corV0CTPCvn = QxcnCor * Qxtn + QycnCor * Qytn;

    hQVzAQVzCvsCentrality->Fill(corV0AV0Cvn, iCen);
    hQVzAQTPCvsCentrality->Fill(corV0ATPCvn, iCen);
    hQVzCQTPCvsCentrality->Fill(corV0CTPCvn, iCen);

    //NUA correction

    hQxVzAvsCentrality->Fill(QxanCor, iCen);
    hQyVzAvsCentrality->Fill(QyanCor, iCen);
    hQxVzCvsCentrality->Fill(QxcnCor, iCen);
    hQyVzCvsCentrality->Fill(QycnCor, iCen);

    fRCollision.fEPangleV0A = evPlAngV0A;
    fRCollision.fEPangleV0C = evPlAngV0C;
    fRCollision.fQA[0] = QxanCor;
    fRCollision.fQA[1] = QyanCor;
    fRCollision.fQC[0] = QxcnCor;
    fRCollision.fQC[1] = QycnCor;
  }

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
  fRHyperTritonFull.clear();
  std::vector<int> he3TrackIndices;

  std::vector<AliESDv0> V0Vector;
  if (!fUseOnTheFly && !fUseNanoAODs)
  {
    if(!fLambda)
      esdEvent->ResetV0s();
    if (!fEnableEventMixing)
      V0Vector = fV0Vertexer.Tracks2V0vertices(esdEvent, fPIDResponse, mcEvent, fLambda);
    else {
      std::vector<AliESDtrack*> he3v, antihe3v;
      for (auto& he3 : fHe3mixed[0][centBin][zBin])
        he3v.push_back(&he3);
      for (auto& ahe3 : fHe3mixed[1][centBin][zBin])
        antihe3v.push_back(&ahe3);
      V0Vector = fV0Vertexer.Tracks2V0verticesEM(esdEvent, fPIDResponse, he3v, antihe3v);
    }
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

      if (fUseTPCmomentum) {
        AliExternalTrackParam extP(*(esdEvent->GetTrack(lKeyPos)->GetInnerParam()));
        AliExternalTrackParam extN(*(esdEvent->GetTrack(lKeyNeg)->GetInnerParam()));
        auto vtx = v0->GetVertex();
        extP.PropagateToDCA(&vtx, esdEvent->GetMagneticField(), 25);
        extN.PropagateToDCA(&vtx, esdEvent->GetMagneticField(), 25);
        extP.GetPxPyPz(pP);
        extN.GetPxPyPz(nP);
      }

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
      if (fMaxInfo) {
        RHyperTritonHe3piFull fullV0 = v0part;
        AliESDtrack* he3t = esdEvent->GetTrack((he3index == lKeyPos) ? lKeyPos : lKeyNeg);
        AliESDtrack* pit = esdEvent->GetTrack((he3index == lKeyPos) ? lKeyNeg : lKeyPos);
        fullV0.fRHe3Track = *he3t;
        fullV0.fRPiTrack = *pit;
        fullV0.fRHe3pidHypo = he3t->GetPIDForTracking();
        fRHyperTritonFull.push_back(fullV0);
      } else
        fRHyperTriton.push_back(v0part);
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
      fRHyperTriton.push_back(v0part);
    }

    he3TrackIndices.push_back(he3index);
  }

  if (fEnableEventMixing) {
    for (int idx : he3TrackIndices) {
      AliESDtrack* he3 = esdEvent->GetTrack(idx);
      fHe3mixed[he3->GetSign() < 0][centBin][zBin].push_back(*he3);
    }
    for (int i = 0; i < 2; ++i) {
      while (fHe3mixed[i][centBin][zBin].size() > fEMdepth) {
        fHe3mixed[i][centBin][zBin].pop_front();
      }
    }
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

  if (fRHyperTriton.size() != 0 || fRHyperTritonFull.size() != 0 || fStoreAllEvents || fMC)
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
bool AliAnalysisTaskHyperTriton2He3piML::FillHyperCandidate(T *v0, AliVEvent *event, AliMCEvent *mcEvent, M mcMap,
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
          fSHyperTriton[mcMap[ilab]].fRecoIndex = fMaxInfo ? (fRHyperTritonFull.size()) : (fRHyperTriton.size());
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

//_____________________________________________________________________________
void AliAnalysisTaskHyperTriton2He3piML::OpenInfoCalibration(int run)
{

  if (fV0CalibrationFile.empty())
    return;
  TFile* foadb = TFile::Open((fV0CalibrationFile + (run < 296624 ? "q" : "r") + ".root").data());

  if(!foadb){
    printf("OADB V0 calibration file cannot be opened\n");
    return;
  }

  AliOADBContainer* cont = (AliOADBContainer*) foadb->Get("hMultV0BefCorPfpx");
  if(!cont){
    printf("OADB object hMultV0BefCorr is not available in the file\n");
    return;
  }
  if(!(cont->GetObject(run))){
    printf("OADB object hMultV0BefCorPfpx is not available for run %i\n", run);
    return;
  }
  fMultV0 = ((TH1D*) cont->GetObject(run));

  AliOADBContainer* contQxnam = 0;
  if (fNHarm == 2)
    contQxnam = (AliOADBContainer*) foadb->Get("fqxa2m");
  else
    contQxnam = (AliOADBContainer*) foadb->Get("fqxa3m");

  if(!contQxnam){
    printf("OADB object fqxanm is not available in the file\n");
    return;
  }
  if(!(contQxnam->GetObject(run))){
    printf("OADB object fqxanm is not available for run %i\n", run);
    return;
  }
  fQxnmV0A = ((TH1D*) contQxnam->GetObject(run));

  AliOADBContainer* contQynam = 0;
  if (fNHarm == 2)
    contQynam = (AliOADBContainer*) foadb->Get("fqya2m");
  else if (fNHarm == 3)
    contQynam = (AliOADBContainer*) foadb->Get("fqya3m");
  else if (fNHarm == 4)
    contQynam = (AliOADBContainer*) foadb->Get("fqya4m");

  if(!contQynam){
    printf("OADB object fqyanm is not available in the file\n");
    return;
  }
  if(!(contQynam->GetObject(run))){
    printf("OADB object fqyanm is not available for run %i\n", run);
    return;
  }
  fQynmV0A = ((TH1D*) contQynam->GetObject(run));

  AliOADBContainer* contQxnas = 0;
  if (fNHarm == 2)
    contQxnas = (AliOADBContainer*) foadb->Get("fqxa2s");
  else
    contQxnas = (AliOADBContainer*) foadb->Get("fqxa3s");

  if(!contQxnas){
    printf("OADB object fqxans is not available in the file\n");
    return;
  }
  if(!(contQxnas->GetObject(run))){
    printf("OADB object fqxans is not available for run %i\n", run);
    return;
  }
  fQxnsV0A = ((TH1D*) contQxnas->GetObject(run));

  AliOADBContainer* contQynas = 0;
  if (fNHarm == 2)
    contQynas = (AliOADBContainer*) foadb->Get("fqya2s");
  else if (fNHarm == 3)
    contQynas = (AliOADBContainer*) foadb->Get("fqya3s");
  else if (fNHarm == 4)
    contQynas = (AliOADBContainer*) foadb->Get("fqya4s");

  if(!contQynas){
    printf("OADB object fqyans is not available in the file\n");
    return;
  }
  if(!(contQynas->GetObject(run))){
    printf("OADB object fqyans is not available for run %i\n", run);
    return;
  }
  fQynsV0A = ((TH1D*) contQynas->GetObject(run));

  AliOADBContainer* contQxncm = 0;
  if (fNHarm == 2)
    contQxncm = (AliOADBContainer*) foadb->Get("fqxc2m");
  else
    contQxncm = (AliOADBContainer*) foadb->Get("fqxc3m");

  if(!contQxncm){
    printf("OADB object fqxcnm is not available in the file\n");
    return;
  }
  if(!(contQxncm->GetObject(run))){
    printf("OADB object fqxcnm is not available for run %i\n", run);
    return;
  }
  fQxnmV0C = ((TH1D*) contQxncm->GetObject(run));

  AliOADBContainer* contQyncm = 0;
  if (fNHarm == 2)
    contQyncm = (AliOADBContainer*) foadb->Get("fqyc2m");
  else if (fNHarm == 3)
    contQyncm = (AliOADBContainer*) foadb->Get("fqyc3m");
  else if (fNHarm == 4)
    contQyncm = (AliOADBContainer*) foadb->Get("fqyc4m");

  if(!contQyncm){
    printf("OADB object fqyc2m is not available in the file\n");
    return;
  }
  if(!(contQyncm->GetObject(run))){
    printf("OADB object fqyc2m is not available for run %i\n", run);
    return;
  }
  fQynmV0C = ((TH1D*) contQyncm->GetObject(run));

  AliOADBContainer* contQxncs = 0;
  if (fNHarm == 2)
    contQxncs = (AliOADBContainer*) foadb->Get("fqxc2s");
  else
    contQxncs = (AliOADBContainer*) foadb->Get("fqxc3s");

  if(!contQxncs){
    printf("OADB object fqxc2s is not available in the file\n");
    return;
  }
  if(!(contQxncs->GetObject(run))){
    printf("OADB object fqxc2s is not available for run %i\n", run);
    return;
  }
  fQxnsV0C = ((TH1D*) contQxncs->GetObject(run));

  AliOADBContainer* contQyncs = 0;
  if (fNHarm == 2)
    contQyncs = (AliOADBContainer*) foadb->Get("fqyc2s");
  else if (fNHarm == 3)
    contQyncs = (AliOADBContainer*) foadb->Get("fqyc3s");
  else if (fNHarm == 4)
    contQyncs = (AliOADBContainer*) foadb->Get("fqyc4s");

  if(!contQyncs){
    printf("OADB object fqycnm is not available in the file\n");
    return;
  }
  if(!(contQyncs->GetObject(run))){
    printf("OADB object fqycns is not available for run %i\n", run);
    return;
  }
  fQynsV0C = ((TH1D*) contQyncs->GetObject(run));
}
