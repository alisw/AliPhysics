#include "AliPHOSEmbedggHBT.h"

#include "AliAODCaloCells.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODVertex.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliCaloPhoton.h"
#include "AliConvEventCuts.h"
#include "AliConversionPhotonCuts.h"
#include "AliEMCALGeometry.h"
#include "AliEPFlattener.h"
#include "AliEventplane.h"
#include "AliFemtoPair.h"
#include "AliFemtoThreeVector.h"
#include "AliFemtoTrack.h"
#include "AliLog.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliMagF.h"
#include "AliMultSelection.h"
#include "AliOADBContainer.h"
#include "AliPHOSCalibData.h"
#include "AliPHOSEsdCluster.h"
#include "AliPHOSGeometry.h"
#include "AliV0ReaderV1.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TF1.h"
#include "TFile.h"
#include "TGeoGlobalMagField.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH2I.h"
#include "TH3F.h"
#include "THashList.h"
#include "TObjArray.h"
#include "TParticle.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TRandom.h"
#include "TStyle.h"
#include "TTree.h"

// Analysis task to fill histograms with PHOS ESD clusters and cells
// Authors: Dmitri Peressounko
// Date   : 28.05.2011

ClassImp(AliPHOSEmbedggHBT)

  //________________________________________________________________________
  AliPHOSEmbedggHBT::AliPHOSEmbedggHBT(const char* name)
  : AliAnalysisTaskEtaPhigg(name), fSignalEvent(nullptr), fSignalEvents(nullptr)
{
  // Constructor
}
//________________________________________________________________________
AliPHOSEmbedggHBT::~AliPHOSEmbedggHBT()
{
  // Note that histograms are stored in fOutputContainer and should not be deleted explicitely
  if (fMCEvents) {
    delete fMCEvents;
    fMCEvents = nullptr;
  }
  if (fMCEvents) {
    delete fMCEvents;
    fMCEvents = nullptr;
  }
  if (fSignalEvent) {
    delete fSignalEvent;
    fSignalEvent = nullptr;
  }
  if (fSignalEvents) {
    delete fSignalEvents;
    fSignalEvents = nullptr;
  }
}
//________________________________________________________________________
void AliPHOSEmbedggHBT::UserCreateOutputObjects()
{
  // Create histograms
  // Called once

  if (fOutputContainer != NULL) {
    delete fOutputContainer;
  }
  fOutputContainer = new THashList();
  fOutputContainer->SetOwner(kTRUE);

  //========QA histograms=======

  // Event selection
  fOutputContainer->Add(new TH1F("hTotSelEvents", "Event selection", 10, 0., 10.));

  fOutputContainer->Add(new TH2F("phiRP", "Event plane", 100, 0., TMath::Pi(), 100, 0., 100.));
  fOutputContainer->Add(new TH2F("phiRPV0A", "Event plane", 100, 0., TMath::Pi(), 100, 0., 100.));
  fOutputContainer->Add(new TH2F("phiRPV0Aflat", "Event plane", 100, 0., TMath::Pi(), 100, 0., 100.));
  fOutputContainer->Add(new TH2F("phiRPV0C", "Event plane", 100, 0., TMath::Pi(), 100, 0., 100.));
  fOutputContainer->Add(new TH2F("phiRPV0Cflat", "Event plane", 100, 0., TMath::Pi(), 100, 0., 100.));
  fOutputContainer->Add(new TH2F("phiRPV0AC", "Event plane", 100, 0., TMath::Pi(), 100, 0., TMath::Pi()));
  fOutputContainer->Add(new TH2F("phiRPV0ACflat", "Event plane", 100, 0., TMath::Pi(), 100, 0., TMath::Pi()));
  fOutputContainer->Add(new TH2F("phiRPvsV0A", "Event plane", 100, 0., TMath::Pi(), 100, 0., TMath::Pi()));
  fOutputContainer->Add(new TH2F("phiRPvsV0C", "Event plane", 100, 0., TMath::Pi(), 100, 0., TMath::Pi()));

  // vertex distribution
  fOutputContainer->Add(new TH1F("hZvertex", "Z vertex position", 50, -25., 25.));

  // Centrality
  fOutputContainer->Add(new TH1F("hCentrality", "Event centrality", 100, 0., 100.));

  fOutputContainer->Add(new TH2F("hSignalTofM1", "Time vs E", 100, 0., 10., 200, -400.e-9, 400.e-9));
  fOutputContainer->Add(new TH2F("hSignalTofM2", "Time vs E", 100, 0., 10., 200, -400.e-9, 400.e-9));
  fOutputContainer->Add(new TH2F("hSignalTofM3", "Time vs E", 100, 0., 10., 200, -400.e-9, 400.e-9));
  fOutputContainer->Add(new TH2F("hSignalTofM4", "Time vs E", 100, 0., 10., 200, -400.e-9, 400.e-9));

  fOutputContainer->Add(new TH2F("hTofM1", "Time vs E", 100, 0., 10., 200, -400.e-9, 400.e-9));
  fOutputContainer->Add(new TH2F("hTofM2", "Time vs E", 100, 0., 10., 200, -400.e-9, 400.e-9));
  fOutputContainer->Add(new TH2F("hTofM3", "Time vs E", 100, 0., 10., 200, -400.e-9, 400.e-9));
  fOutputContainer->Add(new TH2F("hTofM4", "Time vs E", 100, 0., 10., 200, -400.e-9, 400.e-9));
  fOutputContainer->Add(new TH2F("hAllSp", "Sp vs cen", 100, 0., 10., 5, 0., 100.));
  fOutputContainer->Add(new TH2F("hCluM1", "Cell (X,Z), M3", 64, 0.5, 64.5, 56, 0.5, 56.5));
  fOutputContainer->Add(new TH2F("hCluM2", "Cell (X,Z), M3", 64, 0.5, 64.5, 56, 0.5, 56.5));
  fOutputContainer->Add(new TH2F("hCluM3", "Cell (X,Z), M3", 64, 0.5, 64.5, 56, 0.5, 56.5));
  fOutputContainer->Add(new TH2F("hCluM4", "Cell (X,Z), M3", 64, 0.5, 64.5, 56, 0.5, 56.5));

  SetCutNames();

  // MC
  fhReQinvMCprim = new TH2F("hReQinvMCprim", "Qinv distribution", 100, 0., 0.25, 40, 0., 2.);
  fOutputContainer->Add(fhReQinvMCprim);
  fhReqMCprim = new TH2F("hReqMCprim", "q distribution", 100, 0., 0.25, 40, 0., 2.);
  fOutputContainer->Add(fhReqMCprim);
  fhReQinvMC = new TH2F("hReQinvMC", "Qinv distribution", 100, 0., 0.25, 40, 0., 2.);
  fOutputContainer->Add(fhReQinvMC);
  fhReqMC = new TH2F("hReqMC", "q distribution", 100, 0., 0.25, 40, 0., 2.);
  fOutputContainer->Add(fhReqMC);
  fhMiQinvMCprim = new TH2F("hMiQinvMCprim", "Qinv distribution", 100, 0., 0.25, 40, 0., 2.);
  fOutputContainer->Add(fhMiQinvMCprim);
  fhMiQinvMC = new TH2F("hMiQinvMC", "Qinv distribution", 100, 0., 0.25, 40, 0., 2.);
  fOutputContainer->Add(fhMiQinvMC);
  fhMiqMCprim = new TH2F("hMiqMCprim", "q distribution", 100, 0., 0.25, 40, 0., 2.);
  fOutputContainer->Add(fhMiqMCprim);
  fhMiqMC = new TH2F("hMiqMC", "q distribution", 100, 0., 0.25, 40, 0., 2.);
  fOutputContainer->Add(fhMiqMC);
  // Signal
  for (Int_t iCut = 0; iCut < kCuts; iCut++) {
    fhReQinvSignal[iCut] =
      new TH2F(Form("hReQinv_signal_%s", fcut[iCut]), "Qinv distribution", 100, 0., 0.25, 40, 0., 2.);
    fOutputContainer->Add(fhReQinvSignal[iCut]);
    fhReQinvSignalConv[iCut] =
      new TH2F(Form("hReQinv_signalConv_%s", fcut[iCut]), "Qinv distribution, conv", 100, 0., 0.25, 40, 0., 2.);
    fOutputContainer->Add(fhReQinvSignalConv[iCut]);
    fhMiQinvSignal[iCut] =
      new TH2F(Form("hMiQinv_signal_%s", fcut[iCut]), "Qinv distribution", 100, 0., 0.25, 40, 0., 2.);
    fOutputContainer->Add(fhMiQinvSignal[iCut]);
    fhReqSignal[iCut] = new TH2F(Form("hReq_signal_%s", fcut[iCut]), "q distribution", 100, 0., 0.25, 40, 0., 2.);
    fOutputContainer->Add(fhReqSignal[iCut]);
    fhReqSignalConv[iCut] =
      new TH2F(Form("hReq_signalConv_%s", fcut[iCut]), "q distribution", 100, 0., 0.25, 40, 0., 2.);
    fOutputContainer->Add(fhReqSignalConv[iCut]);
    fhMiqSignal[iCut] = new TH2F(Form("hMiq_signal_%s", fcut[iCut]), "q distribution", 100, 0., 0.25, 40, 0., 2.);
    fOutputContainer->Add(fhMiqSignal[iCut]);
  }

  // Embedding
  for (Int_t cen = 0; cen < kCentBins; cen++) {
    for (Int_t iCut = 0; iCut < kCuts; iCut++) {
      fhReQinv[cen][iCut] =
        new TH2F(Form("hReQinv_%s_cen%d", fcut[iCut], cen), "Qinv distribution", 100, 0., 0.25, 40, 0., 2.);
      fOutputContainer->Add(fhReQinv[cen][iCut]);
      fhReQinvConv[cen][iCut] =
        new TH2F(Form("hReQinvConv_%s_cen%d", fcut[iCut], cen), "Qinv distribution", 100, 0., 0.25, 40, 0., 2.);
      fOutputContainer->Add(fhReQinvConv[cen][iCut]);
      fhMiQinv[cen][iCut] =
        new TH2F(Form("hMiQinv_%s_cen%d", fcut[iCut], cen), "Qinv distribution", 100, 0., 0.25, 40, 0., 2.);
      fOutputContainer->Add(fhMiQinv[cen][iCut]);
      fhReq[cen][iCut] = new TH2F(Form("hReq_%s_cen%d", fcut[iCut], cen), "q distribution", 100, 0., 0.25, 40, 0., 2.);
      fOutputContainer->Add(fhReq[cen][iCut]);
      fhReqConv[cen][iCut] =
        new TH2F(Form("hReqConv_%s_cen%d", fcut[iCut], cen), "q distribution", 100, 0., 0.25, 40, 0., 2.);
      fOutputContainer->Add(fhReqConv[cen][iCut]);
      fhMiq[cen][iCut] = new TH2F(Form("hMiq_%s_cen%d", fcut[iCut], cen), "q distribution", 100, 0., 0.25, 40, 0., 2.);
      fOutputContainer->Add(fhMiq[cen][iCut]);
    }
  }

  fSignalEvents = new TList();
  fMCEvents = new TList();

  PostData(1, fOutputContainer);
}
//________________________________________________________________________
void AliPHOSEmbedggHBT::UserExec(Option_t*)
{
  printf(" %s: Start event %d\n", GetName(), fEventCounter);
  // Main loop, called for each event
  FillHistogram("hTotSelEvents", 0.5);

  fEvent = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!fEvent) {
    Printf("ERROR: AliPHOSEmbedggHBT::UserExec Could not retrieve event");
    return;
  }
  if (!fEvent->FindListObject("SignalCaloClusters")) {
    Printf("ERROR: can not find embedded and signal event");
    return;
  }

  if (!fPHOSGeo)
    fPHOSGeo = AliPHOSGeometry::GetInstance();

  fMF = fEvent->GetMagneticField();

  // Checks if we have a primary vertex
  // Get primary vertices form AOD
  const AliAODVertex* esdVertex5 = fEvent->GetPrimaryVertex();

  double vtx5[3] = { 0., 0., 0. };

  vtx5[0] = esdVertex5->GetX();
  vtx5[1] = esdVertex5->GetY();
  vtx5[2] = esdVertex5->GetZ();

  FillHistogram("hZvertex", esdVertex5->GetZ());
  if (TMath::Abs(esdVertex5->GetZ()) > 10.) {
    PostData(1, fOutputContainer);
    return;
  }
  FillHistogram("hTotSelEvents", 2.5);

  // Vtx class z-bin
  Int_t zvtx = (Int_t)((vtx5[2] + 10.) / 1.);
  if (zvtx < 0)
    zvtx = 0;
  if (zvtx >= kVtxBins)
    zvtx = kVtxBins - 1;
  FillHistogram("hTotSelEvents", 3.5);

  fCentrality = 300;
  AliMultSelection* MultSelection = (AliMultSelection*)fEvent->FindListObject("MultSelection");
  if (!MultSelection) {
    AliWarning("AliMultSelection object not found!");
  } else {
    fCentrality = MultSelection->GetMultiplicityPercentile("V0M");
  }
  FillHistogram("hCentrality", fCentrality);

  if (fCentrality < 0. || fCentrality > 80.) {
    PostData(1, fOutputContainer);
    return;
  }
  FillHistogram("hTotSelEvents", 4.5);

  if (fCentrality < 5.)
    fCenBin = 0;
  else if (fCentrality < 10.)
    fCenBin = 1;
  else if (fCentrality < 20.)
    fCenBin = 2;
  else if (fCentrality < 30.)
    fCenBin = 3;
  else if (fCentrality < 40.)
    fCenBin = 4;
  else if (fCentrality < 50.)
    fCenBin = 5;
  else if (fCentrality < 80.)
    fCenBin = 6;

  // reaction plane
  AliEventplane* eventPlane = fEvent->GetEventplane();
  if (!eventPlane) { // Event has no event plane
    PostData(1, fOutputContainer);
    return;
  }
  // V0A
  const Int_t harmonics = 2;
  double qx = 0., qy = 0.;
  // Whole V0
  fRP = eventPlane->CalculateVZEROEventPlane(fEvent, 10, harmonics, qx, qy);
  while (fRP < 0)
    fRP += TMath::TwoPi() / harmonics;
  while (fRP > TMath::TwoPi() / harmonics)
    fRP -= TMath::TwoPi() / harmonics;

  // Reaction plane is defined in the range (0;pi)
  // We have 10 bins
  Int_t irp = Int_t(kPRBins * (fRP) / TMath::Pi());
  if (irp < 0)
    irp = 0;
  if (irp >= kPRBins)
    irp = kPRBins - 1;

  if (!fPHOSEvents[zvtx][fCenBin][irp])
    fPHOSEvents[zvtx][fCenBin][irp] = new TList();
  TList* prevPHOS = fPHOSEvents[zvtx][fCenBin][irp];

  if (fSignalEvent)
    fSignalEvent->Clear();
  else
    fSignalEvent = new TClonesArray("AliCaloPhoton", 100);

  if (fPHOSEvent)
    fPHOSEvent->Clear();
  else
    fPHOSEvent = new TClonesArray("AliCaloPhoton", 100);

  if (fCPVEvent)
    fCPVEvent->Clear();
  else
    fCPVEvent = new TClonesArray("AliCaloPhoton", 100);

  TClonesArray* embedded = static_cast<TClonesArray*>(fEvent->FindListObject("EmbeddedCaloClusters"));
  int multClust = embedded->GetEntriesFast();
  // CPV clusters
  Float_t position[3];
  int inCPV = 0;
  const double dxCPV = 0.9;
  const double dyCPV = -2.;           // V3: 5. V2: -2;  V1:18
  const double dzCPV = 4.9;           // V3: 4.6
  const double slopeZCPV = -0.034745; // tilt of module

  for (Int_t j = 0; j < multClust; j++) {
    AliAODCaloCluster* cluCPV = static_cast<AliAODCaloCluster*>(embedded->At(j));
    if (cluCPV->GetType() != AliVCluster::kPHOSCharged)
      continue;
    if (cluCPV->E() < 100.)
      continue;
    cluCPV->GetPosition(position);
    position[0] -= dxCPV;
    position[1] -= dyCPV;
    position[2] -= dzCPV + slopeZCPV * position[2];
    TVector3 globalCPV(position);
    Float_t dxMin, dzMin;
    AliCaloPhoton* p =
      new ((*fCPVEvent)[inCPV]) AliCaloPhoton(globalCPV.Px(), globalCPV.Py(), globalCPV.Pz(), globalCPV.Mag());
    inCPV++;
    int itr = FindTrackMatching(0, globalCPV, 3, dxMin, dzMin);
    p->SetLambdas(dxMin, dzMin);
    p->SetDistToBad(itr);
  }

  Int_t inPHOS = 0, inSignal = 0;
  TVector3 localPos;

  double vtx0[3] = { 0., 0., 0. };

  TClonesArray* signal = static_cast<TClonesArray*>(fEvent->FindListObject("SignalCaloClusters"));
  int multSignal = signal->GetEntriesFast();
  for (Int_t i = 0; i < multSignal; i++) {
    AliAODCaloCluster* clu = static_cast<AliAODCaloCluster*>(signal->At(i));
    if (clu->GetType() != AliVCluster::kPHOSNeutral)
      continue;
    if (clu->E() < 0.1)
      continue;

    Float_t position[3];
    clu->GetPosition(position);
    TVector3 global(position);
    Int_t relId[4];
    fPHOSGeo->GlobalPos2RelId(global, relId);
    Int_t mod = relId[0];
    Int_t cellX = relId[2];
    Int_t cellZ = relId[3];
    TVector3 local;
    fPHOSGeo->Global2Local(local, global, mod);

    FillHistogram(Form("hSignalTofM%d", mod), clu->E(), clu->GetTOF());
    if ((clu->GetTOF() > 30.e-9) || (clu->GetTOF() < -30.e-7))
      continue;

    if (clu->GetNCells() < 2)
      continue;
    if (clu->GetM02() < 0.2)
      continue;

    TLorentzVector pv1;
    clu->GetMomentum(pv1, vtx0);

    if (inSignal >= fSignalEvent->GetSize()) {
      fSignalEvent->Expand(inSignal + 20);
    }
    AliCaloPhoton* ph = new ((*fSignalEvent)[inSignal++]) AliCaloPhoton(pv1.X(), pv1.Py(), pv1.Z(), pv1.E());
    ph->SetModule(mod);
    pv1 *= clu->GetCoreEnergy() / pv1.E();
    ph->SetMomV2(&pv1);
    ph->SetNCells(clu->GetNCells());
    ph->SetDispBit(clu->Chi2() < 2.5 * 2.5);

    float dxMin, dzMin;
    int itr = FindTrackMatching(1, global, mod, dxMin, dzMin);

    // Set veto bits, true: neutral
    int cpvBits = TestCPV(mod, clu->E(), local.X(), local.Z(), dxMin, dzMin, itr);

    ph->SetTagInfo(cpvBits);

    ph->SetEMCx(local.X());
    ph->SetEMCz(local.Z());
    ph->SetLambdas(clu->GetM20(), clu->GetM02());
    ph->SetUnfolded(clu->GetNExMax() < 2); // Remember, if it is unfolded
    ph->SetDistToBad(cellX);
    ph->SetPrimary(clu->GetLabelAt(0));
  }

  for (Int_t i = 0; i < multClust; i++) {
    AliAODCaloCluster* clu = static_cast<AliAODCaloCluster*>(embedded->At(i));
    if (clu->GetType() != AliVCluster::kPHOSNeutral)
      continue;
    if (clu->E() < 0.1)
      continue;
    if (clu->GetNLabels() == 0)
      continue;

    Float_t position[3];
    clu->GetPosition(position);
    TVector3 global(position);
    Int_t relId[4];
    fPHOSGeo->GlobalPos2RelId(global, relId);
    Int_t mod = relId[0];
    Int_t cellX = relId[2];
    Int_t cellZ = relId[3];

    TVector3 local;
    fPHOSGeo->Global2Local(local, global, mod);

    FillHistogram(Form("hTofM%d", mod), clu->E(), clu->GetTOF());
    if ((clu->GetTOF() > 150.e-9) || (clu->GetTOF() < -150.e-7))
      continue;

    if (clu->GetNCells() < 2)
      continue;
    if (clu->GetM02() < 0.2)
      continue;

    TLorentzVector pv1;
    clu->GetMomentum(pv1, vtx5);

    FillHistogram(Form("hCluM%d", mod), cellX, cellZ, 1.);

    FillHistogram("hAllSp", clu->E(), fCentrality);

    if (inPHOS >= fPHOSEvent->GetSize()) {
      fPHOSEvent->Expand(inPHOS + 20);
    }
    AliCaloPhoton* ph = new ((*fPHOSEvent)[inPHOS]) AliCaloPhoton(pv1.X(), pv1.Py(), pv1.Z(), pv1.E());
    inPHOS++;
    ph->SetModule(mod);
    pv1 *= clu->GetCoreEnergy() / pv1.E();
    ph->SetMomV2(&pv1);
    ph->SetNCells(clu->GetNCells());
    bool disp = clu->Chi2() < 2.5 * 2.5;
    ph->SetDispBit(disp);

    float dxMin, dzMin;
    int itr = FindTrackMatching(1, global, mod, dxMin, dzMin);

    // Set veto bits, true: neutral
    int cpvBits = TestCPV(mod, clu->E(), local.X(), local.Z(), dxMin, dzMin, itr);

    ph->SetTagInfo(cpvBits);

    ph->SetEMCx(local.X());
    ph->SetEMCz(local.Z());
    ph->SetLambdas(clu->GetM20(), clu->GetM02());
    ph->SetUnfolded(clu->GetNExMax() < 2); // Remember, if it is unfolded
    ph->SetDistToBad(cellX);
    ph->SetPrimary(clu->GetLabelAt(0));
  }

  // Real
  TClonesArray* mc = static_cast<TClonesArray*>(fEvent->FindListObject(AliAODMCParticle::StdBranchName()));
  TLorentzVector ph1, ph2;
  TLorentzVector rec1, rec2;
  for (Int_t i1 = 0; i1 < mc->GetEntriesFast() - 1; i1++) {
    AliAODMCParticle* p1 = (AliAODMCParticle*)mc->At(i1);
    p1->Momentum(ph1);

    // Module acceptance
    double vtx[3] = { p1->Xv(), p1->Yv(), p1->Zv() };
    int mod1;
    double x1, z1;
    if (!fPHOSGeo->ImpactOnEmc(vtx, p1->Theta(), p1->Phi(), mod1, x1, z1))
      continue;

    TVector3 globaPos;
    fPHOSGeo->Local2Global(mod1, x1, z1, globaPos);
    globaPos.SetMag(p1->E());
    rec1.SetPxPyPzE(globaPos.X(), globaPos.Y(), globaPos.Z(), p1->E());
    for (Int_t i2 = i1 + 1; i2 < mc->GetEntriesFast(); i2++) {
      AliAODMCParticle* p2 = (AliAODMCParticle*)mc->At(i2);

      p2->Momentum(ph2);
      // Module acceptance
      double vtx2[3] = { p2->Xv(), p2->Yv(), p2->Zv() };
      int mod2;
      double x2, z2;
      if (!fPHOSGeo->ImpactOnEmc(vtx2, p2->Theta(), p2->Phi(), mod2, x2, z2)) {
        continue;
      }

      if (mod1 != mod2) {
        continue;
      }
      fPHOSGeo->Local2Global(mod2, x2, z2, globaPos);
      globaPos.SetMag(p2->E());
      rec2.SetPxPyPzE(globaPos.X(), globaPos.Y(), globaPos.Z(), p2->E());

      TLorentzVector sum = ph1 + ph2;
      TVector3 gammaBeta(sum.BoostVector());
      gammaBeta.SetXYZ(0, 0, gammaBeta.Z());
      TLorentzVector gammaCMq(ph1 - ph2);
      gammaCMq.Boost(-gammaBeta);
      double q = gammaCMq.Vect().Mag();

      double qinv = sum.M();
      double kT = 0.5 * sum.Pt();

      TLorentzVector sumRec = rec1 + rec2;
      gammaBeta = sumRec.BoostVector();
      gammaBeta.SetXYZ(0, 0, gammaBeta.Z());
      gammaCMq = rec1 - rec2;
      gammaCMq.Boost(-gammaBeta);
      double qRec = gammaCMq.Vect().Mag();

      double qinvRec = sumRec.M();
      double kTRec = 0.5 * sumRec.Pt();

      // primary
      if (p1->GetMother() == -1 && p2->GetMother() == -1) {
        fhReQinvMCprim->Fill(qinv, kT);
        fhReqMCprim->Fill(q, kT);
      }
      // all
      fhReQinvMC->Fill(qinvRec, kTRec);
      fhReqMC->Fill(qRec, kTRec);
    }
  }

  // Signal
  for (Int_t i1 = 0; i1 < inSignal - 1; i1++) {
    AliCaloPhoton* ph1 = (AliCaloPhoton*)fSignalEvent->At(i1);

    for (Int_t i2 = i1 + 1; i2 < inSignal; i2++) {
      AliCaloPhoton* ph2 = (AliCaloPhoton*)fSignalEvent->At(i2);

      TLorentzVector sum = *ph1 + *ph2;
      TVector3 gammaBeta(sum.BoostVector());
      gammaBeta.SetXYZ(0, 0, gammaBeta.Z());
      TLorentzVector gammaCMq(*ph1 - *ph2);
      gammaCMq.Boost(-gammaBeta);
      double q = gammaCMq.Vect().Mag();

      double qinv = sum.M();
      double kT = 0.5 * sum.Pt();

      bool commonParent = (CommonParent(ph1, ph2) != -1);
      for (Int_t iCut = 0; iCut < kCuts; iCut++) {
        if (!PairCut(ph1, ph2, iCut)) {
          continue;
        }
        fhReQinvSignal[iCut]->Fill(qinv, kT);
        fhReqSignal[iCut]->Fill(q, kT);
        if (commonParent) {
          fhReQinvSignalConv[iCut]->Fill(qinv, kT);
          fhReqSignalConv[iCut]->Fill(q, kT);
        }
      }
    }
  }
  for (Int_t i1 = 0; i1 < inPHOS - 1; i1++) {
    AliCaloPhoton* ph1 = (AliCaloPhoton*)fPHOSEvent->At(i1);
    for (Int_t i2 = i1 + 1; i2 < inPHOS; i2++) {
      AliCaloPhoton* ph2 = (AliCaloPhoton*)fPHOSEvent->At(i2);
      TLorentzVector sum(*ph1 + *ph2);
      double qinv = sum.M();
      double kT = 0.5 * sum.Pt();
      TVector3 gammaBeta(sum.BoostVector());
      gammaBeta.SetXYZ(0, 0, gammaBeta.Z());
      TLorentzVector gammaCMq(*ph1 - *ph2);
      gammaCMq.Boost(-gammaBeta);
      double q = gammaCMq.Vect().Mag();

      bool commonParent = (CommonParent(ph1, ph2) != -1);
      for (Int_t iCut = 0; iCut < kCuts; iCut++) {
        if (!PairCut(ph1, ph2, iCut)) {
          continue;
        }
        fhReQinv[fCenBin][iCut]->Fill(qinv, kT);
        fhReq[fCenBin][iCut]->Fill(q, kT);
        if (commonParent) {
          fhReQinvConv[fCenBin][iCut]->Fill(qinv, kT);
          fhReqConv[fCenBin][iCut]->Fill(q, kT);
        }
      }
    }
  }

  // now mixed
  for (Int_t i1 = 0; i1 < mc->GetEntriesFast(); i1++) {
    AliAODMCParticle* p1 = (AliAODMCParticle*)mc->At(i1);
    p1->Momentum(ph1);
    // Module acceptance
    double vtx[3] = { p1->Xv(), p1->Yv(), p1->Zv() };
    int mod1;
    double x1, z1;
    if (!fPHOSGeo->ImpactOnEmc(vtx, p1->Theta(), p1->Phi(), mod1, x1, z1))
      continue;

    TVector3 globaPos;
    fPHOSGeo->Local2Global(mod1, x1, z1, globaPos);
    globaPos.SetMag(p1->E());
    rec1.SetPxPyPzE(globaPos.X(), globaPos.Y(), globaPos.Z(), p1->E());

    for (Int_t ev = 0; ev < fMCEvents->GetSize(); ev++) {
      TClonesArray* mixMC = static_cast<TClonesArray*>(fMCEvents->At(ev));
      for (Int_t i2 = 0; i2 < mixMC->GetEntriesFast(); i2++) {
        AliAODMCParticle* p2 = (AliAODMCParticle*)mixMC->At(i2);
        p2->Momentum(ph2);
        // Module acceptance
        double vtx2[3] = { p2->Xv(), p2->Yv(), p2->Zv() };
        int mod2;
        double x2, z2;
        if (!fPHOSGeo->ImpactOnEmc(vtx2, p2->Theta(), p2->Phi(), mod2, x2, z2)) {
          continue;
        }

        if (mod1 != mod2) {
          continue;
        }

        fPHOSGeo->Local2Global(mod2, x2, z2, globaPos);
        globaPos.SetMag(p2->E());
        rec2.SetPxPyPzE(globaPos.X(), globaPos.Y(), globaPos.Z(), p2->E());

        TLorentzVector sum = ph1 + ph2;
        TVector3 gammaBeta(sum.BoostVector());
        gammaBeta.SetXYZ(0, 0, gammaBeta.Z());
        TLorentzVector gammaCMq(ph1 - ph2);
        gammaCMq.Boost(-gammaBeta);
        double q = gammaCMq.Vect().Mag();

        double qinv = sum.M();
        double kT = 0.5 * sum.Pt();

        TLorentzVector sumRec = rec1 + rec2;
        gammaBeta = sumRec.BoostVector();
        gammaBeta.SetXYZ(0, 0, gammaBeta.Z());
        gammaCMq = rec1 - rec2;
        gammaCMq.Boost(-gammaBeta);
        double qRec = gammaCMq.Vect().Mag();

        double qinvRec = sumRec.M();
        double kTRec = 0.5 * sumRec.Pt();

        // primary
        if (p1->GetMother() == -1 && p2->GetMother() == -1) {
          fhMiQinvMCprim->Fill(qinv, kT);
          fhMiqMCprim->Fill(q, kT);
        }
        fhMiQinvMC->Fill(qinvRec, kTRec);
        fhMiqMC->Fill(qRec, kTRec);
      }
    }
  }

  // mixed-Signal
  for (Int_t i1 = 0; i1 < inSignal; i1++) {
    AliCaloPhoton* ph1 = (AliCaloPhoton*)fSignalEvent->At(i1);

    for (Int_t ev = 0; ev < fSignalEvents->GetSize(); ev++) {
      TClonesArray* mixSignal = static_cast<TClonesArray*>(fSignalEvents->At(ev));
      for (Int_t i2 = 0; i2 < mixSignal->GetEntriesFast(); i2++) {
        AliCaloPhoton* ph2 = (AliCaloPhoton*)mixSignal->At(i2);

        TLorentzVector sum = *ph1 + *ph2;
        TVector3 gammaBeta(sum.BoostVector());
        gammaBeta.SetXYZ(0, 0, gammaBeta.Z());
        TLorentzVector gammaCMq(*ph1 - *ph2);
        gammaCMq.Boost(-gammaBeta);
        double q = gammaCMq.Vect().Mag();

        double qinv = sum.M();
        double kT = 0.5 * sum.Pt();

        for (Int_t iCut = 0; iCut < kCuts; iCut++) {
          if (!PairCut(ph1, ph2, iCut)) {
            continue;
          }
          fhMiQinvSignal[iCut]->Fill(qinv, kT);
          fhMiqSignal[iCut]->Fill(q, kT);
        }
      }
    }
  }
  // Embedded
  for (Int_t i1 = 0; i1 < inPHOS; i1++) {
    AliCaloPhoton* ph1 = (AliCaloPhoton*)fPHOSEvent->At(i1);
    for (Int_t ev = 0; ev < prevPHOS->GetSize(); ev++) {
      TClonesArray* mixPHOS = static_cast<TClonesArray*>(prevPHOS->At(ev));
      for (Int_t i2 = 0; i2 < mixPHOS->GetEntriesFast(); i2++) {
        AliCaloPhoton* ph2 = (AliCaloPhoton*)mixPHOS->At(i2);
        TLorentzVector sum(*ph1 + *ph2);
        double qinv = sum.M();
        double kT = 0.5 * sum.Pt();
        TVector3 gammaBeta(sum.BoostVector());
        gammaBeta.SetXYZ(0, 0, gammaBeta.Z());
        TLorentzVector gammaCMq(*ph1 - *ph2);
        gammaCMq.Boost(-gammaBeta);
        double q = gammaCMq.Vect().Mag();

        for (Int_t iCut = 0; iCut < kCuts; iCut++) {
          if (!PairCut(ph1, ph2, iCut)) {
            continue;
          }
          fhMiQinv[fCenBin][iCut]->Fill(qinv, kT);
          fhMiq[fCenBin][iCut]->Fill(q, kT);
        }
      }
    }
  }

  if (mc->GetEntriesFast() > 0) {
    TClonesArray* tmpMC = new TClonesArray("AliAODMCParticle", mc->GetEntriesFast());
    for (int i = 0; i < mc->GetEntriesFast(); i++) {
      AliAODMCParticle* a = (AliAODMCParticle*)mc->At(i);
      new ((*tmpMC)[i]) AliAODMCParticle(*a);
    }
    fMCEvents->AddFirst(tmpMC);
    if (fMCEvents->GetSize() > 5) { // Remove redundant events
      TClonesArray* tmp = static_cast<TClonesArray*>(fMCEvents->Last());
      fMCEvents->RemoveLast();
      delete tmp;
    }
  }

  // Now we either add current events to stack or remove
  // If no photons in current event - no need to add it to mixed
  const Int_t kMixEvents[kCentBins] = { 20, 20, 20, 30, 30, 40, 40 };
  if (fPHOSEvent->GetEntriesFast() > 0) {
    fPHOSEvent->Expand(fPHOSEvent->GetEntriesFast());
    prevPHOS->AddFirst(fPHOSEvent);
    fPHOSEvent = 0;
    if (prevPHOS->GetSize() > kMixEvents[fCenBin]) { // Remove redundant events
      TClonesArray* tmp = static_cast<TClonesArray*>(prevPHOS->Last());
      prevPHOS->RemoveLast();
      delete tmp;
    }
  }

  if (fSignalEvent->GetEntriesFast() > 0) {
    fSignalEvent->Expand(fSignalEvent->GetEntriesFast());
    fSignalEvents->AddFirst(fSignalEvent);
    fSignalEvent = 0;
    if (fSignalEvents->GetSize() > kMixEvents[0]) { // Remove redundant events
      TClonesArray* tmp = static_cast<TClonesArray*>(fSignalEvents->Last());
      fSignalEvents->RemoveLast();
      delete tmp;
    }
  }

  // Post output data.
  PostData(1, fOutputContainer);
  fEventCounter++;
}
//___________________________________________________________________________
Int_t AliPHOSEmbedggHBT::JetRejection(Int_t module) const
{
  // We reject events with hard tracks in the vicinity of center of PHOS module
  // track pT thresholds
  const double cutPt1 = 10.; // Maximal Track pT
  const double cutPt2 = 5.;  // Track pT
  const double cutPt3 = 3.;  // Track pT

  // Cone threshold
  const double cutR = 0.5 * 0.5; // Cone radius squared

  Int_t result = 0;

  // azimuthal angle of PHOS module
  double phiPHOS = TMath::DegToRad() * (270. + fPHOSGeo->GetPHOSAngle(module)); // (40,20,0,-20,-40) degrees

  Int_t nt = fEvent->GetNumberOfTracks();
  for (Int_t i = 0; i < nt; i++) {
    AliAODTrack* aodTrack = static_cast<AliAODTrack*>(fEvent->GetTrack(i));

    if (aodTrack->GetTPCncls() < 70) {
      continue;
    }
    if (aodTrack->GetTPCchi2() / Float_t(aodTrack->GetTPCncls()) < 4) {
      continue;
    }
    UInt_t status = aodTrack->GetStatus();
    if ((status & AliAODTrack::kTPCrefit) == 0) {
      continue;
    }
    if (aodTrack->GetKinkIndex(0) > 0)
      continue;
    if (Float_t(aodTrack->GetTPCnclsS()) / Float_t(aodTrack->GetTPCncls()) > 0.4)
      continue;
    // ITS
    if ((status & AliAODTrack::kITSrefit) == 0)
      continue;

    if (aodTrack->Pt() < cutPt3)
      continue;
    double dphi = aodTrack->Phi() - phiPHOS;
    while (dphi < -TMath::Pi())
      dphi += TMath::TwoPi();
    while (dphi > TMath::Pi())
      dphi -= TMath::TwoPi();
    double deta = aodTrack->Eta();
    double r = dphi * dphi + deta * deta;
    if (r < cutR) {
      result = result | 1 << 0;
      if (aodTrack->Pt() > cutPt2) {
        result = result | 1 << 1;
      }
      if (aodTrack->Pt() > cutPt1) {
        result = result | 1 << 2;
      }
      return result;
    }
  }
  return result;
}
//_____________________________________________________________________________
void AliPHOSEmbedggHBT::ReclusterizeCPV()
{
  // Jet-finder like algorithm:
  // find the highest cell, add m*n cells around it,
  // next highest cell, etc. while there will be no cells above the threshold
  typedef std::pair<double, int> pairs;

  const double kSeedAmp = 10.; // Threshold for cluster seed

  const double logWeight = 4.5;

  // CPV geometry
  const Int_t nCPVPadsZ = 60;
  const Int_t nCPVPadsX = 128;

  // CPV cluster parameters
  const Int_t nCluX = 5; // cluster size
  const Int_t nCluZ = 3; // cluster size

  double cpvAmp[nCPVPadsX + 1][nCPVPadsZ + 1] = { 0 };

  Int_t nCpvClu = 0;

  // Copy CPV digits to 2D matrix
  AliAODEvent* event = dynamic_cast<AliAODEvent*>(InputEvent());
  AliAODCaloCells* cells = event->GetPHOSCells();
  const Int_t nCells = cells->GetNumberOfCells();
  Int_t relId[4];
  pairs arr[nCells];
  Int_t nCPV = 0;
  for (Int_t iCell = 0; iCell < nCells; iCell++) {
    Int_t cellAbsId = cells->GetCellNumber(iCell);
    if (cellAbsId > 0)
      continue; // Select CPV cells
    cellAbsId = -cellAbsId + 56 * 64 * 5;

    fPHOSGeo->AbsToRelNumbering(cellAbsId, relId);
    double amp = cells->GetAmplitude(iCell);
    cpvAmp[relId[2]][relId[3]] = amp;
    if (amp > kSeedAmp) {
      arr[nCPV].first = amp;
      arr[nCPV].second = cellAbsId;
      nCPV++;
    }
  }

  // Sort in increasing energy
  std::sort(arr, arr + nCPV);

  // Find cell with highest amplitude and if it is above threshold, combine cluster n*m and remove cells from the matrix
  for (Int_t iCPV = nCPV - 1; iCPV >= 0; iCPV--) {
    fPHOSGeo->AbsToRelNumbering(arr[iCPV].second, relId);
    if (cpvAmp[relId[2]][relId[3]] > 0) {
      // Create CPV cluster with seed in this seed
      Int_t xMin = TMath::Max(relId[2] - (nCluX - 1) / 2, 1);
      Int_t xMax = TMath::Min(relId[2] + (nCluX - 1) / 2, nCPVPadsX);
      Int_t zMin = TMath::Max(relId[3] - (nCluZ - 1) / 2, 1);
      Int_t zMax = TMath::Min(relId[3] + (nCluZ - 1) / 2, nCPVPadsZ);
      double wtot = 0.;
      double x = 0., z = 0.;
      double totE = 0.;
      for (Int_t ix = xMin; ix <= xMax; ix++) {
        for (Int_t iz = zMin; iz <= zMax; iz++) {
          totE += cpvAmp[ix][iz];
        }
      }
      for (Int_t ix = xMin; ix <= xMax; ix++) {
        for (Int_t iz = zMin; iz <= zMax; iz++) {
          Float_t xi = 0., zi = 0.;
          relId[0] = 3;
          relId[1] = 1;
          relId[2] = ix;
          relId[3] = iz;
          fPHOSGeo->RelPosInModule(relId, xi, zi);

          if (cpvAmp[ix][iz] > 0) {
            double w = TMath::Max(0., logWeight + TMath::Log(cpvAmp[ix][iz] / totE));
            x += xi * w;
            z += zi * w;
            wtot += w;

            cpvAmp[ix][iz] = 0; // remove cells
          }
        }
      }
      if (wtot != 0) {
        x /= wtot;
        z /= wtot;
      } else {
        x = 999;
        z = 999;
      }
      // Account mis-alignment
      TVector3 globaPos;
      fPHOSGeo->Local2Global(3, x, z, globaPos);

      double glZ = globaPos.Z() - 3.25 - 5.21 + 0.067 * globaPos.Z();
      double glX = globaPos.X() - 0.6 - 0.48 + 0.067 * globaPos.X();

      if (fCPVEvent->GetSize() <= nCpvClu)
        fCPVEvent->Expand(nCpvClu * 2);
      new ((*fCPVEvent)[nCpvClu++]) TVector3(glX, arr[iCPV].first, glZ);
    }
  }
}
//_________________________________________________________________________
bool AliPHOSEmbedggHBT::IsGoodChannel(Int_t cellX, Int_t cellZ)
{
  if (fBadMap)
    if (fBadMap->GetBinContent(cellX, cellZ) == 1)
      return kFALSE;

  return kTRUE;
}
//_________________________________________________________________________
int AliPHOSEmbedggHBT::CommonParent(const AliCaloPhoton* p1, const AliCaloPhoton* p2) const
{
  // Looks through parents and finds if there was commont pi0 among ancestors
  TClonesArray* fStack = static_cast<TClonesArray*>(fEvent->FindListObject(AliAODMCParticle::StdBranchName()));

  if (!fStack) {
    return -1; // can not say anything
  }
  Int_t prim1 = p1->GetPrimary();
  while (prim1 != -1) {
    Int_t prim2 = p2->GetPrimary();

    while (prim2 != -1) {
      if (prim1 == prim2) {
        return prim1;
      }
      prim2 = ((AliAODMCParticle*)fStack->At(prim2))->GetMother();
    }
    prim1 = ((AliAODMCParticle*)fStack->At(prim1))->GetMother();
  }
  return -1;
}
