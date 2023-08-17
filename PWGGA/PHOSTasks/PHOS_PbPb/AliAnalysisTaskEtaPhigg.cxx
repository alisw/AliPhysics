#include "AliAnalysisTaskEtaPhigg.h"

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

// Task for ggHBT analysis in PHOS
// Authors: Dmitri Peressounko
// Date   : 01.01.2023

ClassImp(AliAnalysisTaskEtaPhigg)

  //________________________________________________________________________
  AliAnalysisTaskEtaPhigg::AliAnalysisTaskEtaPhigg(const char* name)
  : AliAnalysisTaskSE(name),
    // fStack(0x0),
    fOutputContainer(nullptr),
    fEvent(nullptr),
    fPHOSEvent(nullptr),
    fCPVEvent(nullptr),
    fV0AFlat(nullptr),
    fV0CFlat(nullptr),
    fRP(0.),
    fMF(0),
    fRunNumber(0),
    fCentrality(0.),
    fCenBin(0),
    fPHOSGeo(nullptr),
    fEventCounter(0),
    fBadMap(nullptr)
{
  // Constructor
  for (Int_t i = 0; i < kVtxBins; i++) {
    for (Int_t j = 0; j < kCentBins; j++) {
      for (Int_t k = 0; k < kPRBins; k++) { // no RP bins
        fPHOSEvents[i][j][k] = nullptr;
      }
    }
  }

  // Output slots #0 write into a TH1 container
  DefineOutput(1, TList::Class());
}
//________________________________________________________________________
void AliAnalysisTaskEtaPhigg::UserCreateOutputObjects()
{
  // Create histograms
  // Called once
  const Int_t nRuns = 200;

  // ESD histograms
  if (fOutputContainer != NULL) {
    delete fOutputContainer;
  }
  fOutputContainer = new THashList();
  fOutputContainer->SetOwner(kTRUE);

  //========QA histograms=======

  // Event selection
  fOutputContainer->Add(new TH2F("hSelEvents", "Event selection", 10, 0., 10., nRuns, 0., float(nRuns)));
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
  fOutputContainer->Add(new TH2F("hZvertex", "Z vertex position", 50, -25., 25., nRuns, 0., float(nRuns)));

  // Centrality
  fOutputContainer->Add(new TH2F("hCentrality", "Event centrality", 100, 0., 100., nRuns, 0., float(nRuns)));
  fOutputContainer->Add(new TH2F("hCenPHOS", "Centrality vs PHOSclusters", 100, 0., 100., 200, 0., 200.));
  fOutputContainer->Add(new TH2F("hCenPHOSm1", "Centrality vs PHOSclusters in mod1 ", 100, 0., 100., 100, 0., 100.));
  fOutputContainer->Add(new TH2F("hCenPHOSm2", "Centrality vs PHOSclusters in mod2 ", 100, 0., 100., 100, 0., 100.));
  fOutputContainer->Add(new TH2F("hCenPHOSm3", "Centrality vs PHOSclusters in mod3 ", 100, 0., 100., 100, 0., 100.));
  fOutputContainer->Add(new TH2F("hCenPHOSm4", "Centrality vs PHOSclusters in mod4 ", 100, 0., 100., 100, 0., 100.));
  fOutputContainer->Add(new TH2F("hCenPHOSCells", "Centrality vs PHOS cells", 100, 0., 100., 100, 0., 1000.));
  fOutputContainer->Add(new TH1F("hPHOSBadMod", "Number of bad event due to module", 10, 0., 10.));

  fOutputContainer->Add(new TH2F("hCenCPV", "Centrality vs CPVclusters", 100, 0., 100., 2000, 0., 2000.));

  fOutputContainer->Add(new TH2F("hCenTrack", "Centrality vs tracks", 100, 0., 100., 100, 0., 15000.));
  fOutputContainer->Add(new TH2F("hCenTOF", "Centrality vs PHOS TOF", 100, 0., 100., 200, -400.e-9, 400.e-9));

  //   fOutputContainer->Add(new TH2F("hCenEMCAL","Centrality vs EMCAL photons", 100,0.,100.,200,0.,1000.)) ;
  //
  //   fOutputContainer->Add(new TH2F("hPHOSvsEMCAL","PHOS vs EMCAL photons", 100,0.,100.,200,0.,1000.)) ;

  // PHOS QA

  // Bad Map
  fOutputContainer->Add(new TH2F("hCluLowM1", "Cell (X,Z), M1", 64, 0.5, 64.5, 56, 0.5, 56.5));
  fOutputContainer->Add(new TH2F("hCluLowM2", "Cell (X,Z), M2", 64, 0.5, 64.5, 56, 0.5, 56.5));
  fOutputContainer->Add(new TH2F("hCluLowM3", "Cell (X,Z), M3", 64, 0.5, 64.5, 56, 0.5, 56.5));
  fOutputContainer->Add(new TH2F("hCluLowM4", "Cell (X,Z), M4", 64, 0.5, 64.5, 56, 0.5, 56.5));
  fOutputContainer->Add(new TH2F("hDisp", "Cell (X,Z), M3", 64, 0.5, 64.5, 56, 0.5, 56.5));
  fOutputContainer->Add(new TH2F("hTrackVeto", "Cell (X,Z), M3", 64, 0.5, 64.5, 56, 0.5, 56.5));
  fOutputContainer->Add(new TH2F("hDispTrackVeto", "Cell (X,Z), M3", 64, 0.5, 64.5, 56, 0.5, 56.5));
  fOutputContainer->Add(new TH2F("hCPVVeto", "Cell (X,Z), M3", 64, 0.5, 64.5, 56, 0.5, 56.5));
  fOutputContainer->Add(new TH2F("hDispCPVVeto", "Cell (X,Z), M3", 64, 0.5, 64.5, 56, 0.5, 56.5));
  fOutputContainer->Add(new TH2F("hBothVeto", "Cell (X,Z) ", 64, 0.5, 64.5, 56, 0.5, 56.5));
  fOutputContainer->Add(new TH2F("hDispBothVeto", "Cell (X,Z)", 64, 0.5, 64.5, 56, 0.5, 56.5));
  fOutputContainer->Add(new TH2F("hNotCPVVeto", "Cell (X,Z)", 64, 0.5, 64.5, 56, 0.5, 56.5));
  fOutputContainer->Add(new TH2F("hDispNotCPVVeto", "Cell (X,Z)", 64, 0.5, 64.5, 56, 0.5, 56.5));
  fOutputContainer->Add(new TH2F("hNotTrackVeto", "Cell (X,Z)", 64, 0.5, 64.5, 56, 0.5, 56.5));
  fOutputContainer->Add(new TH2F("hDispNotTrackVeto", "Cell (X,Z)", 64, 0.5, 64.5, 56, 0.5, 56.5));
  fOutputContainer->Add(new TH2F("hNotBothVeto", "Cell (X,Z)", 64, 0.5, 64.5, 56, 0.5, 56.5));
  fOutputContainer->Add(new TH2F("hDispNotBothVeto", "Cell (X,Z)", 64, 0.5, 64.5, 56, 0.5, 56.5));

  fOutputContainer->Add(new TH2F("hTofM1", "TOF in M3", 100, 0., 20., 400, -4.e-6, 4.e-6));
  fOutputContainer->Add(new TH2F("hTofM2", "TOF in M3", 100, 0., 20., 400, -4.e-6, 4.e-6));
  fOutputContainer->Add(new TH2F("hTofM3", "TOF in M3", 100, 0., 20., 400, -4.e-6, 4.e-6));
  fOutputContainer->Add(new TH2F("hTofM4", "TOF in M3", 100, 0., 20., 400, -4.e-6, 4.e-6));

  fOutputContainer->Add(new TH2F("hAllSp", "Spectrum in PHOS", 100, 0., 10., 20, 0., 100.));
  fOutputContainer->Add(new TH2F("hMod3Sp", "Spectrum in PHOS", 100, 0., 10., 20, 0., 100.));
  fOutputContainer->Add(new TH2F("hDispSp", "Spectrum in PHOS", 100, 0., 10., 20, 0., 100.));
  fOutputContainer->Add(new TH2F("hTrackVetoSp", "Spectrum in PHOS", 100, 0., 10., 20, 0., 100.));
  fOutputContainer->Add(new TH2F("hDispTrackVetoSp", "Spectrum in PHOS", 100, 0., 10., 20, 0., 100.));
  fOutputContainer->Add(new TH2F("hCPVVetoSp", "Spectrum in PHOS", 100, 0., 10., 20, 0., 100.));
  fOutputContainer->Add(new TH2F("hDispCPVVetoSp", "Spectrum in PHOS", 100, 0., 10., 20, 0., 100.));
  fOutputContainer->Add(new TH2F("hBothVetoSp", "Spectrum in PHOS", 100, 0., 10., 20, 0., 100.));
  fOutputContainer->Add(new TH2F("hDispBothVetoSp", "Spectrum in PHOS", 100, 0., 10., 20, 0., 100.));
  fOutputContainer->Add(new TH2F("hNotCPVVetoSp", "Spectrum in PHOS", 100, 0., 10., 20, 0., 100.));
  fOutputContainer->Add(new TH2F("hDispNotCPVVetoSp", "Spectrum in PHOS", 100, 0., 10., 20, 0., 100.));
  fOutputContainer->Add(new TH2F("hNotTrackVetoSp", "Spectrum in PHOS", 100, 0., 10., 20, 0., 100.));
  fOutputContainer->Add(new TH2F("hDispNotTrackVetoSp", "Spectrum in PHOS", 100, 0., 10., 20, 0., 100.));
  fOutputContainer->Add(new TH2F("hNotBothVetoSp", "Spectrum in PHOS", 100, 0., 10., 20, 0., 100.));
  fOutputContainer->Add(new TH2F("hDispNotBothVetoSp", "Spectrum in PHOS", 100, 0., 10., 20, 0., 100.));

  SetCutNames();
  // HBT part
  for (Int_t cen = 0; cen < kCentBins; cen++) {
    for (Int_t iCut = 0; iCut < kCuts; iCut++) {
      fhReQinv[cen][iCut] =
        new TH2F(Form("hReQinv_%s_cen%d", fcut[iCut], cen), "Qinv distribution", 100, 0., 0.25, 40, 0., 2.);
      fOutputContainer->Add(fhReQinv[cen][iCut]);
      fhMiQinv[cen][iCut] =
        new TH2F(Form("hMiQinv_%s_cen%d", fcut[iCut], cen), "Qinv distribution", 100, 0., 0.25, 40, 0., 2.);
      fOutputContainer->Add(fhMiQinv[cen][iCut]);
      fhReq[cen][iCut] =
        new TH2F(Form("hReq_%s_cen%d", fcut[iCut], cen), "Qinv distribution", 100, 0., 0.25, 40, 0., 2.);
      fOutputContainer->Add(fhReq[cen][iCut]);
      fhMiq[cen][iCut] =
        new TH2F(Form("hMiq_%s_cen%d", fcut[iCut], cen), "Qinv distribution", 100, 0., 0.25, 40, 0., 2.);
      fOutputContainer->Add(fhMiq[cen][iCut]);
    }
  }

  // 3D part
  Int_t nQ = 120;
  Double_t qMax = 0.3;
  for (Int_t cen = 0; cen < kCentBins; cen++) {
    for (Int_t ikT = 0; ikT < kKtbins; ikT++) {
      fhReOSL[cen][ikT] = new TH3F(Form("hReOSL_DzCr_Kt%3.1f-%3.1f_cen%d", fKtBins[ikT], fKtBins[ikT + 1], cen),
                                   "real Out-Side-Long", nQ, -qMax, qMax, nQ, -qMax, qMax, nQ, -qMax, qMax);
      fOutputContainer->Add(fhReOSL[cen][ikT]);
      fhMiOSL[cen][ikT] = new TH3F(Form("hMiOSL_DzCr_Kt%3.1f-%3.1f_cen%d", fKtBins[ikT], fKtBins[ikT + 1], cen),
                                   "real Out-Side-Long", nQ, -qMax, qMax, nQ, -qMax, qMax, nQ, -qMax, qMax);
      fOutputContainer->Add(fhMiOSL[cen][ikT]);
      fhReOSLCTS[cen][ikT] = new TH3F(Form("hReOSL_DzCTS_Kt%3.1f-%3.1f_cen%d", fKtBins[ikT], fKtBins[ikT + 1], cen),
                                      "real Out-Side-Long", nQ, -qMax, qMax, nQ, -qMax, qMax, nQ, -qMax, qMax);
      fOutputContainer->Add(fhReOSLCTS[cen][ikT]);
      fhMiOSLCTS[cen][ikT] = new TH3F(Form("hMiOSL_DzCTS_Kt%3.1f-%3.1f_cen%d", fKtBins[ikT], fKtBins[ikT + 1], cen),
                                      "real Out-Side-Long", nQ, -qMax, qMax, nQ, -qMax, qMax, nQ, -qMax, qMax);
      fOutputContainer->Add(fhMiOSLCTS[cen][ikT]);
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskEtaPhigg::UserExec(Option_t*)
{
  // Main loop, called for each event
  FillHistogram("hTotSelEvents", 0.5);

  fEvent = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!fEvent) {
    Printf("ERROR: AliAnalysisTaskEtaPhigg::UserExec Could not retrieve event");
    return;
  }

  fRunNumber = ConvertRunNumber(fEvent->GetRunNumber());
  FillHistogram("hSelEvents", 1.5, fRunNumber - 0.5);
  FillHistogram("hTotSelEvents", 1.5);

  if (!fPHOSGeo) {
    fPHOSGeo = AliPHOSGeometry::GetInstance();
    fMF = fEvent->GetMagneticField();
  }

  // Checks if we have a primary vertex
  // Get primary vertices form AOD
  const AliAODVertex* esdVertex5 = fEvent->GetPrimaryVertex();

  double vtx5[3] = { 0., 0., 0. };

  vtx5[0] = esdVertex5->GetX();
  vtx5[1] = esdVertex5->GetY();
  vtx5[2] = esdVertex5->GetZ();

  FillHistogram("hZvertex", esdVertex5->GetZ(), fRunNumber - 0.5);
  if (TMath::Abs(esdVertex5->GetZ()) > 10.) {
    PostData(1, fOutputContainer);
    return;
  }
  FillHistogram("hSelEvents", 2.5, fRunNumber - 0.5);
  FillHistogram("hTotSelEvents", 2.5);

  if (fEvent->IsPileupFromSPD()) {
    PostData(1, fOutputContainer);
    return;
  }
  // Vtx class z-bin
  Int_t zvtx = (Int_t)((vtx5[2] + 10.) / 1.);
  if (zvtx < 0)
    zvtx = 0;
  if (zvtx >= kVtxBins)
    zvtx = kVtxBins - 1;

  FillHistogram("hSelEvents", 3.5, fRunNumber - 0.5);
  FillHistogram("hTotSelEvents", 3.5);

  fCentrality = 300;
  AliMultSelection* MultSelection = (AliMultSelection*)fEvent->FindListObject("MultSelection");
  if (!MultSelection) {
    AliWarning("AliMultSelection object not found!");
  } else {
    fCentrality = MultSelection->GetMultiplicityPercentile("V0M");
  }

  if (fCentrality < 0. || fCentrality > 80.) {
    PostData(1, fOutputContainer);
    return;
  }

  FillHistogram("hSelEvents", 4.5, fRunNumber - 0.5);
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

  FillHistogram("hSelEvents", 5.5, fRunNumber - 0.5);
  FillHistogram("hTotSelEvents", 5.5);

  FillHistogram("hSelEvents", 6.5, fRunNumber - 0.5);
  FillHistogram("hTotSelEvents", 6.5);

  // reaction plane
  AliEventplane* eventPlane = fEvent->GetEventplane();
  if (!eventPlane) { // Event has no event plane
    PostData(1, fOutputContainer);
    return;
  }
  // V0A
  const Int_t harmonics = 2;
  double qx = 0., qy = 0.;
  double rpV0A = eventPlane->CalculateVZEROEventPlane(fEvent, 8, harmonics, qx, qy);
  // V0C
  double rpV0C = eventPlane->CalculateVZEROEventPlane(fEvent, 9, harmonics, qx, qy);

  // Whole V0
  fRP = eventPlane->CalculateVZEROEventPlane(fEvent, 10, harmonics, qx, qy);

  while (rpV0A < 0)
    rpV0A += TMath::TwoPi() / harmonics;
  while (rpV0A > TMath::TwoPi() / harmonics)
    rpV0A -= TMath::TwoPi() / harmonics;

  while (rpV0C < 0)
    rpV0C += TMath::TwoPi() / harmonics;
  while (rpV0C > TMath::TwoPi() / harmonics)
    rpV0C -= TMath::TwoPi() / harmonics;

  while (fRP < 0)
    fRP += TMath::TwoPi() / harmonics;
  while (fRP > TMath::TwoPi() / harmonics)
    fRP -= TMath::TwoPi() / harmonics;

  FillHistogram("phiRP", fRP, fCentrality);
  FillHistogram("phiRPV0A", rpV0A, fCentrality);
  FillHistogram("phiRPV0C", rpV0C, fCentrality);
  FillHistogram("phiRPV0AC", rpV0A, rpV0C);

  //   rpV0A = fV0AFlat->MakeFlat(rpV0A,fCentrality) ;
  //   rpV0C = fV0CFlat->MakeFlat(rpV0C,fCentrality) ;
  //   FillHistogram("phiRPV0Aflat",rpV0A,fCentrality) ;
  //   FillHistogram("phiRPV0Cflat",rpV0C,fCentrality) ;
  //   FillHistogram("phiRPV0ACflat",rpV0A,rpV0C) ;

  FillHistogram("phiRPvsV0A", fRP, rpV0A);
  FillHistogram("phiRPvsV0C", fRP, rpV0C);

  FillHistogram("hSelEvents", 7.5, fRunNumber - 0.5);
  FillHistogram("hTotSelEvents", 7.5);
  // All event selections done
  FillHistogram("hCentrality", fCentrality, fRunNumber - 0.5);
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

  if (fPHOSEvent)
    fPHOSEvent->Clear();
  else
    fPHOSEvent = new TClonesArray("AliCaloPhoton", 100);

  if (fCPVEvent)
    fCPVEvent->Clear();
  else
    fCPVEvent = new TClonesArray("AliCaloPhoton", 100);

  TVector3 vertex(vtx5);

  Int_t multClust = fEvent->GetNumberOfCaloClusters();
  Int_t inPHOS = 0;

  AliAODCaloCells* cells = fEvent->GetPHOSCells();
  FillHistogram("hCenPHOSCells", fCentrality, cells->GetNumberOfCells());
  FillHistogram("hCenTrack", fCentrality, fEvent->GetNumberOfTracks());

  // CPV clusters
  Float_t position[3];
  int inCPV = 0;
  const double dxCPV = 0.9;
  const double dyCPV = -2.;           // V3: 5. V2: -2;  V1:18
  const double dzCPV = 4.9;           // V3: 4.6
  const double slopeZCPV = -0.034745; // tilt of module
  for (Int_t j = 0; j < multClust; j++) {
    AliAODCaloCluster* cluCPV = fEvent->GetCaloCluster(j);
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

  TVector3 localPos;
  for (Int_t i = 0; i < multClust; i++) {
    AliAODCaloCluster* clu = fEvent->GetCaloCluster(i);
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

    FillHistogram(Form("hTofM%d", mod), clu->E(), clu->GetTOF());
    if (clu->E() > 1.)
      FillHistogram("hCenTOF", fCentrality, clu->GetTOF());
    if ((clu->GetTOF() > 150.e-9) || (clu->GetTOF() < -150.e-7))
      continue;
    // if((clu->GetTOF()>100.e-9) || (clu->GetTOF() <-100.e-7) )
    //   continue ;

    if (clu->GetNCells() < 2)
      continue;
    if (clu->GetM02() < 0.2)
      continue;

    TLorentzVector pv1;
    clu->GetMomentum(pv1, vtx5);

    FillHistogram(Form("hCluLowM%d", mod), cellX, cellZ, 1.);

    FillHistogram("hAllSp", clu->E(), fCentrality);
    bool disp = clu->Chi2() < 2.5 * 2.5;
    if (disp) {
      FillHistogram("hDisp", cellX, cellZ, 1.);
      FillHistogram("hDispSp", clu->E(), fCentrality);
    }

    float dxMin, dzMin;
    int itr = FindTrackMatching(1, global, mod, dxMin, dzMin);

    // Set veto bits, true: neutral
    int cpvBits = TestCPV(mod, clu->E(), local.X(), local.Z(), dxMin, dzMin, itr);

    if (inPHOS >= fPHOSEvent->GetSize()) {
      fPHOSEvent->Expand(inPHOS + 20);
    }
    AliCaloPhoton* ph = new ((*fPHOSEvent)[inPHOS]) AliCaloPhoton(pv1.X(), pv1.Py(), pv1.Z(), pv1.E());
    ph->SetModule(mod);
    pv1 *= clu->GetCoreEnergy() / pv1.E();
    ph->SetMomV2(&pv1);
    ph->SetNCells(clu->GetNCells());
    ph->SetDispBit(disp);

    ph->SetTagInfo(cpvBits);

    ph->SetEMCx(local.X());
    ph->SetEMCz(local.Z());
    ph->SetDistToBad(cellX);
    ph->SetLambdas(clu->GetM20(), clu->GetM02());
    ph->SetUnfolded(clu->GetNExMax() < 2); // Remember, if it is unfolded
    inPHOS++;
  }

  // Real
  // PHOS-PHOS
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

      for (Int_t iCut = 0; iCut < kCuts; iCut++) {
        if (!PairCut(ph1, ph2, iCut)) {
          continue;
        }
        fhReQinv[fCenBin][iCut]->Fill(qinv, kT);
        fhReq[fCenBin][iCut]->Fill(q, kT);
      }
      if (kT > fKtBins[0] && kT < fKtBins[kKtbins]) {
        int iKt = 0;
        while (kT > fKtBins[iKt + 1]) {
          iKt++;
        }
        double qo = 0.5 * (sum.Px() * gammaCMq.Px() + sum.Py() * gammaCMq.Py()) / kT;
        double qs = (ph1->Px() * ph2->Py() - ph2->Px() * ph1->Py()) / kT;
        double ql = gammaCMq.Pz();
        if (gRandom->Uniform() < 0.5) { // remove ordering during reconstruction
          qo = -qo;
          qs = -qs;
          ql = -ql;
        }

        if (PairCut(ph1, ph2, 170)) {
          fhReOSLCTS[fCenBin][iKt]->Fill(qo, qs, ql);
        }
        if (PairCut(ph1, ph2, 175)) {
          fhReOSL[fCenBin][iKt]->Fill(qo, qs, ql);
        }
      }
    }
  }

  // now mixed
  // mixed-PHOS-PHOS
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
        if (kT > fKtBins[0] && kT < fKtBins[kKtbins]) {
          int iKt = 0;
          while (kT > fKtBins[iKt + 1]) {
            iKt++;
          }

          double qo = 0.5 * (sum.Px() * gammaCMq.Px() + sum.Py() * gammaCMq.Py()) / kT;
          double qs = (ph1->Px() * ph2->Py() - ph2->Px() * ph1->Py()) / kT;
          double ql = gammaCMq.Pz();

          if (PairCut(ph1, ph2, 170)) {
            fhMiOSLCTS[fCenBin][iKt]->Fill(qo, qs, ql);
          }
          if (PairCut(ph1, ph2, 175)) {
            fhMiOSL[fCenBin][iKt]->Fill(qo, qs, ql);
          }
        }
      }
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

  // Post output data.
  PostData(1, fOutputContainer);
  fEventCounter++;
}

//________________________________________________________________________
void AliAnalysisTaskEtaPhigg::Terminate(Option_t*) {}
//_____________________________________________________________________________
void AliAnalysisTaskEtaPhigg::FillHistogram(const char* key, double x) const
{
  // FillHistogram
  TObject* tmp = fOutputContainer->FindObject(key);
  if (!tmp) {
    AliInfo(Form("can not find histogram <%s> ", key));
    return;
  }
  if (tmp->IsA() == TClass::GetClass("TH1I")) {
    ((TH1I*)tmp)->Fill(x);
    return;
  }
  if (tmp->IsA() == TClass::GetClass("TH1F")) {
    ((TH1F*)tmp)->Fill(x);
    return;
  }
  if (tmp->IsA() == TClass::GetClass("TH1D")) {
    ((TH1D*)tmp)->Fill(x);
    return;
  }
  AliInfo(Form("can not find 1D histogram <%s> ", key));
}
//_____________________________________________________________________________
void AliAnalysisTaskEtaPhigg::FillHistogram(const char* key, double x, double y) const
{
  // FillHistogram
  TObject* tmp = fOutputContainer->FindObject(key);
  if (!tmp) {
    AliInfo(Form("can not find histogram <%s> ", key));
    return;
  }
  if (tmp->IsA() == TClass::GetClass("TH1F")) {
    ((TH1F*)tmp)->Fill(x, y);
    return;
  }
  if (tmp->IsA() == TClass::GetClass("TH2F")) {
    ((TH2F*)tmp)->Fill(x, y);
    return;
  }
  AliError(Form("Calling FillHistogram with 2 parameters for histo <%s> of type %s", key, tmp->IsA()->GetName()));
}

//_____________________________________________________________________________
void AliAnalysisTaskEtaPhigg::FillHistogram(const char* key, double x, double y, double z) const
{
  // Fills 1D histograms with Form(
  TObject* tmp = fOutputContainer->FindObject(key);
  if (!tmp) {
    AliInfo(Form("can not find histogram <%s> ", key));
    return;
  }
  if (tmp->IsA() == TClass::GetClass("TH2F")) {
    ((TH2F*)tmp)->Fill(x, y, z);
    return;
  }
  if (tmp->IsA() == TClass::GetClass("TH3F")) {
    ((TH3F*)tmp)->Fill(x, y, z);
    return;
  }
}
//_____________________________________________________________________________
void AliAnalysisTaskEtaPhigg::FillHistogram(const char* key, double x, double y, double z, double w) const
{
  // Fills 1D histograms with Form(
  TObject* tmp = fOutputContainer->FindObject(key);
  if (!tmp) {
    AliInfo(Form("can not find histogram <%s> ", key));
    return;
  }
  if (tmp->IsA() == TClass::GetClass("TH3F")) {
    ((TH3F*)tmp)->Fill(x, y, z, w);
    return;
  }
}

//___________________________________________________________________________
Int_t AliAnalysisTaskEtaPhigg::ConvertRunNumber(Int_t run)
{
  switch (run) {
    default:
      return 199;
  }
}

//___________________________________________________________________________
bool AliAnalysisTaskEtaPhigg::PairCut(const AliCaloPhoton* ph1, const AliCaloPhoton* ph2, Int_t cut) const
{
  // First distance cut based on non-overlapping CPV cuts
  //  if(ph1->Module()!=ph2->Module()){
  //    return false ;
  //  }
  double r = ph1->Angle(ph2->Vect());
  const double r1cut = TMath::ATan2(5., 460.);
  const double r2cut = TMath::ATan2(10., 460.);
  const double r3cut = TMath::ATan2(15., 460.);
  const double r4cut = TMath::ATan2(20., 460.);
  const double r5cut = TMath::ATan2(25., 460.);

  int cpvBits1 = ph1->GetTagInfo();
  bool track1CPV = cpvBits1 & (1 << 1);
  bool track1CPV2 = cpvBits1 & (1 << 2);
  bool cpv1CPV = cpvBits1 & (1 << 3);
  bool cpv1CPV2 = cpvBits1 & (1 << 4);
  bool cpvAndTrack1 = cpvBits1 & (1 << 5);
  int cpvBits2 = ph2->GetTagInfo();
  bool track2CPV = cpvBits2 & (1 << 1);
  bool track2CPV2 = cpvBits2 & (1 << 2);
  bool cpv2CPV = cpvBits2 & (1 << 3);
  bool cpv2CPV2 = cpvBits2 & (1 << 4);
  bool cpvAndTrack2 = cpvBits2 & (1 << 5);
  // Cuts without cross-talks
  bool badPair = false;
  if (ph1->Module() == 1 && ph2->Module() == 1) { // noisy regions 49,50
    badPair =
      ((ph1->DistToBad() == 49) || (ph1->DistToBad() == 50)) && ((ph2->DistToBad() == 49) || (ph2->DistToBad() == 50));
  }
  if (ph1->Module() == 2 && ph2->Module() == 2) { // noisy regions 6,7,  18, 34,37,38, 63
    badPair =
      (((ph1->DistToBad() == 6) || (ph1->DistToBad() == 7)) && ((ph2->DistToBad() == 6) || (ph2->DistToBad() == 7))) ||
      (((ph1->DistToBad() == 34) || (ph1->DistToBad() == 37) || (ph1->DistToBad() == 38)) &&
       ((ph2->DistToBad() == 34) || (ph2->DistToBad() == 37) || (ph2->DistToBad() == 38)));
  }
  if (ph1->Module() == 3 && ph2->Module() == 3) { // noisy regions 6,  18,21,22,  33,34,37,38,  49,50,63
    badPair =
      (((ph1->DistToBad() == 18) || (ph1->DistToBad() == 21) || (ph1->DistToBad() == 22)) &&
       ((ph2->DistToBad() == 18) || (ph2->DistToBad() == 21) || (ph2->DistToBad() == 22))) ||
      (((ph1->DistToBad() == 33) || (ph1->DistToBad() == 34) || (ph1->DistToBad() == 37) || (ph1->DistToBad() == 38)) &&
       ((ph2->DistToBad() == 33) || (ph2->DistToBad() == 34) || (ph2->DistToBad() == 37) ||
        (ph2->DistToBad() == 38))) ||
      (((ph1->DistToBad() == 49) || (ph1->DistToBad() == 50) || (ph1->DistToBad() == 63)) &&
       ((ph2->DistToBad() == 49) || (ph2->DistToBad() == 50) || (ph2->DistToBad() == 63)));
  }
  if (ph1->Module() == 4 && ph2->Module() == 4) { // noisy regions 34, 49,50
    badPair = (((ph1->DistToBad() == 49) || (ph1->DistToBad() == 50)) &&
               ((ph2->DistToBad() == 49) || (ph2->DistToBad() == 50)));
  }

  switch (cut) {
    case 0:
      return kTRUE;
    case 1:
      return r > r1cut;
    case 2:
      return r > r1cut && ph1->IsDispOK() && ph2->IsDispOK();
    case 3:
      return r > r1cut && track1CPV && track2CPV;
    case 4:
      return r > r1cut && track1CPV && ph1->IsDispOK() && track2CPV && ph2->IsDispOK();
    case 5:
      return r > r1cut && ph1->Module() == 3 && ph2->Module() == 3;
    case 6:
      return r > r1cut && ph1->Module() == 3 && ph2->Module() == 3 && cpv1CPV && cpv2CPV;
    case 7:
      return r > r1cut && ph1->Module() == 3 && ph2->Module() == 3 && cpv1CPV && cpv2CPV && ph1->IsDispOK() &&
             ph2->IsDispOK();
    case 8:
      return r > r1cut && ph1->Module() == 3 && ph2->Module() == 3 && cpv1CPV2 && cpv2CPV2;
    case 9:
      return r > r1cut && ph1->Module() == 3 && ph2->Module() == 3 && cpv1CPV2 && cpv2CPV2 && ph1->IsDispOK() &&
             ph2->IsDispOK();
    case 10:
      return r > r1cut && ph1->Module() == 3 && ph2->Module() == 3 && cpvAndTrack1 && cpvAndTrack2;
    // D10
    case 11:
      return r > r2cut && !badPair;
    case 12:
      return r > r2cut && !badPair && ph1->IsDispOK() && ph2->IsDispOK();
    case 13:
      return r > r2cut && !badPair && track1CPV && track2CPV;
    case 14:
      return r > r2cut && !badPair && track1CPV2 && track2CPV2;
    case 15:
      return r > r2cut && !badPair && track1CPV && track2CPV && ph1->IsDispOK() && ph2->IsDispOK();
    case 16:
      return r > r2cut && !badPair && track1CPV2 && track2CPV2 && ph1->IsDispOK() && ph2->IsDispOK();
    case 17:
      return r > r2cut && !badPair && ph1->Module() == 3 && ph2->Module() == 3;
    case 18:
      return r > r2cut && !badPair && ph1->Module() == 3 && ph2->Module() == 3 && cpv1CPV && cpv2CPV;
    case 19:
      return r > r2cut && !badPair && ph1->Module() == 3 && ph2->Module() == 3 && cpv1CPV2 && cpv2CPV2;
    case 20:
      return r > r2cut && !badPair && ph1->Module() == 3 && ph2->Module() == 3 && cpv1CPV && cpv2CPV &&
             ph1->IsDispOK() && ph2->IsDispOK();
    case 21:
      return r > r2cut && !badPair && ph1->Module() == 3 && ph2->Module() == 3 && cpv1CPV2 && cpv2CPV2 &&
             ph1->IsDispOK() && ph2->IsDispOK();
    case 22:
      return r > r2cut && !badPair && ph1->Module() == 3 && ph2->Module() == 3 && cpvAndTrack1 && cpvAndTrack2;
    // D10 Emin 0.2
    case 23:
      return r > r2cut && !badPair && ph1->E() > 0.2 && ph2->E() > 0.2;
    case 24:
      return r > r2cut && !badPair && ph1->E() > 0.2 && ph2->E() > 0.2 && ph1->IsDispOK() && ph2->IsDispOK();
    case 25:
      return r > r2cut && !badPair && ph1->E() > 0.2 && ph2->E() > 0.2 && track1CPV && track2CPV;
    case 26:
      return r > r2cut && !badPair && ph1->E() > 0.2 && ph2->E() > 0.2 && track1CPV2 && track2CPV2;
    case 27:
      return r > r2cut && !badPair && ph1->E() > 0.2 && ph2->E() > 0.2 && track1CPV && track2CPV && ph1->IsDispOK() &&
             ph2->IsDispOK();
    case 28:
      return r > r2cut && !badPair && ph1->E() > 0.2 && ph2->E() > 0.2 && track1CPV2 && track2CPV2 && ph1->IsDispOK() &&
             ph2->IsDispOK();
    case 29:
      return r > r2cut && !badPair && ph1->E() > 0.2 && ph2->E() > 0.2 && ph1->Module() == 3 && ph2->Module() == 3;
    case 30:
      return r > r2cut && !badPair && ph1->E() > 0.2 && ph2->E() > 0.2 && ph1->Module() == 3 && ph2->Module() == 3 &&
             cpv1CPV && cpv2CPV;
    case 31:
      return r > r2cut && !badPair && ph1->E() > 0.2 && ph2->E() > 0.2 && ph1->Module() == 3 && ph2->Module() == 3 &&
             cpv1CPV2 && cpv2CPV2;
    case 32:
      return r > r2cut && !badPair && ph1->E() > 0.2 && ph2->E() > 0.2 && ph1->Module() == 3 && ph2->Module() == 3 &&
             cpv1CPV && cpv2CPV && ph1->IsDispOK() && ph2->IsDispOK();
    case 33:
      return r > r2cut && !badPair && ph1->E() > 0.2 && ph2->E() > 0.2 && ph1->Module() == 3 && ph2->Module() == 3 &&
             cpv1CPV2 && cpv2CPV2 && ph1->IsDispOK() && ph2->IsDispOK();
    case 34:
      return r > r2cut && !badPair && ph1->E() > 0.2 && ph2->E() > 0.2 && ph1->Module() == 3 && ph2->Module() == 3 &&
             cpvAndTrack1 && cpvAndTrack2;
    // D10 Emin 0.3
    case 35:
      return r > r2cut && !badPair && ph1->E() > 0.3 && ph2->E() > 0.3;
    case 36:
      return r > r2cut && !badPair && ph1->E() > 0.3 && ph2->E() > 0.3 && ph1->IsDispOK() && ph2->IsDispOK();
    case 37:
      return r > r2cut && !badPair && ph1->E() > 0.3 && ph2->E() > 0.3 && track1CPV && track2CPV;
    case 38:
      return r > r2cut && !badPair && ph1->E() > 0.3 && ph2->E() > 0.3 && track1CPV2 && track2CPV2;
    case 39:
      return r > r2cut && !badPair && ph1->E() > 0.3 && ph2->E() > 0.3 && track1CPV && track2CPV && ph1->IsDispOK() &&
             ph2->IsDispOK();
    case 40:
      return r > r2cut && !badPair && ph1->E() > 0.3 && ph2->E() > 0.3 && track1CPV2 && track2CPV2 && ph1->IsDispOK() &&
             ph2->IsDispOK();
    case 41:
      return r > r2cut && !badPair && ph1->E() > 0.3 && ph2->E() > 0.3 && ph1->Module() == 3 && ph2->Module() == 3;
    case 42:
      return r > r2cut && !badPair && ph1->E() > 0.3 && ph2->E() > 0.3 && ph1->Module() == 3 && ph2->Module() == 3 &&
             cpv1CPV && cpv2CPV;
    case 43:
      return r > r2cut && !badPair && ph1->E() > 0.3 && ph2->E() > 0.3 && ph1->Module() == 3 && ph2->Module() == 3 &&
             cpv1CPV2 && cpv2CPV2;
    case 44:
      return r > r2cut && !badPair && ph1->E() > 0.3 && ph2->E() > 0.3 && ph1->Module() == 3 && ph2->Module() == 3 &&
             cpv1CPV && cpv2CPV && ph1->IsDispOK() && ph2->IsDispOK();
    case 45:
      return r > r2cut && !badPair && ph1->E() > 0.3 && ph2->E() > 0.3 && ph1->Module() == 3 && ph2->Module() == 3 &&
             cpv1CPV2 && cpv2CPV2 && ph1->IsDispOK() && ph2->IsDispOK();
    case 46:
      return r > r2cut && !badPair && ph1->E() > 0.3 && ph2->E() > 0.3 && ph1->Module() == 3 && ph2->Module() == 3 &&
             cpvAndTrack1 && cpvAndTrack2;
    // D15
    case 47:
      return r > r3cut;
    case 48:
      return r > r3cut && ph1->IsDispOK() && ph2->IsDispOK();
    case 49:
      return r > r3cut && track1CPV && track2CPV;
    case 50:
      return r > r3cut && track1CPV2 && track2CPV2;
    case 51:
      return r > r3cut && track1CPV && track2CPV && ph1->IsDispOK() && ph2->IsDispOK();
    case 52:
      return r > r3cut && track1CPV2 && track2CPV2 && ph1->IsDispOK() && ph2->IsDispOK();
    case 53:
      return r > r3cut && ph1->Module() == 3 && ph2->Module() == 3;
    case 54:
      return r > r3cut && ph1->Module() == 3 && ph2->Module() == 3 && cpv1CPV && cpv2CPV;
    case 55:
      return r > r3cut && ph1->Module() == 3 && ph2->Module() == 3 && cpv1CPV2 && cpv2CPV2;
    case 56:
      return r > r3cut && ph1->Module() == 3 && ph2->Module() == 3 && cpv1CPV && cpv2CPV && ph1->IsDispOK() &&
             ph2->IsDispOK();
    case 57:
      return r > r3cut && ph1->Module() == 3 && ph2->Module() == 3 && cpv1CPV2 && cpv2CPV2 && ph1->IsDispOK() &&
             ph2->IsDispOK();
    case 58:
      return r > r3cut && ph1->Module() == 3 && ph2->Module() == 3 && cpvAndTrack1 && cpvAndTrack2;
    // D15 Emin 0.2
    case 59:
      return r > r3cut && ph1->E() > 0.2 && ph2->E() > 0.2;
    case 60:
      return r > r3cut && ph1->E() > 0.2 && ph2->E() > 0.2 && ph1->IsDispOK() && ph2->IsDispOK();
    case 61:
      return r > r3cut && ph1->E() > 0.2 && ph2->E() > 0.2 && track1CPV && track2CPV;
    case 62:
      return r > r3cut && ph1->E() > 0.2 && ph2->E() > 0.2 && track1CPV2 && track2CPV2;
    case 63:
      return r > r3cut && ph1->E() > 0.2 && ph2->E() > 0.2 && track1CPV && track2CPV && ph1->IsDispOK() &&
             ph2->IsDispOK();
    case 64:
      return r > r3cut && ph1->E() > 0.2 && ph2->E() > 0.2 && track1CPV2 && track2CPV2 && ph1->IsDispOK() &&
             ph2->IsDispOK();
    case 65:
      return r > r3cut && ph1->E() > 0.2 && ph2->E() > 0.2 && ph1->Module() == 3 && ph2->Module() == 3;
    case 66:
      return r > r3cut && ph1->E() > 0.2 && ph2->E() > 0.2 && ph1->Module() == 3 && ph2->Module() == 3 && cpv1CPV &&
             cpv2CPV;
    case 67:
      return r > r3cut && ph1->E() > 0.2 && ph2->E() > 0.2 && ph1->Module() == 3 && ph2->Module() == 3 && cpv1CPV2 &&
             cpv2CPV2;
    case 68:
      return r > r3cut && ph1->E() > 0.2 && ph2->E() > 0.2 && ph1->Module() == 3 && ph2->Module() == 3 && cpv1CPV &&
             cpv2CPV && ph1->IsDispOK() && ph2->IsDispOK();
    case 69:
      return r > r3cut && ph1->E() > 0.2 && ph2->E() > 0.2 && ph1->Module() == 3 && ph2->Module() == 3 && cpv1CPV2 &&
             cpv2CPV2 && ph1->IsDispOK() && ph2->IsDispOK();
    case 70:
      return r > r3cut && ph1->E() > 0.2 && ph2->E() > 0.2 && ph1->Module() == 3 && ph2->Module() == 3 &&
             cpvAndTrack1 && cpvAndTrack2;
    // D15 Emin 0.3
    case 71:
      return r > r3cut && ph1->E() > 0.3 && ph2->E() > 0.3;
    case 72:
      return r > r3cut && ph1->E() > 0.3 && ph2->E() > 0.3 && ph1->IsDispOK() && ph2->IsDispOK();
    case 73:
      return r > r3cut && ph1->E() > 0.3 && ph2->E() > 0.3 && track1CPV && track2CPV;
    case 74:
      return r > r3cut && ph1->E() > 0.3 && ph2->E() > 0.3 && track1CPV2 && track2CPV2;
    case 75:
      return r > r3cut && ph1->E() > 0.3 && ph2->E() > 0.3 && track1CPV && track2CPV && ph1->IsDispOK() &&
             ph2->IsDispOK();
    case 76:
      return r > r3cut && ph1->E() > 0.3 && ph2->E() > 0.3 && track1CPV2 && track2CPV2 && ph1->IsDispOK() &&
             ph2->IsDispOK();
    case 77:
      return r > r3cut && ph1->E() > 0.3 && ph2->E() > 0.3 && ph1->Module() == 3 && ph2->Module() == 3;
    case 78:
      return r > r3cut && ph1->E() > 0.3 && ph2->E() > 0.3 && ph1->Module() == 3 && ph2->Module() == 3 && cpv1CPV &&
             cpv2CPV;
    case 79:
      return r > r3cut && ph1->E() > 0.3 && ph2->E() > 0.3 && ph1->Module() == 3 && ph2->Module() == 3 && cpv1CPV2 &&
             cpv2CPV2;
    case 80:
      return r > r3cut && ph1->E() > 0.3 && ph2->E() > 0.3 && ph1->Module() == 3 && ph2->Module() == 3 && cpv1CPV &&
             cpv2CPV && ph1->IsDispOK() && ph2->IsDispOK();
    case 81:
      return r > r3cut && ph1->E() > 0.3 && ph2->E() > 0.3 && ph1->Module() == 3 && ph2->Module() == 3 && cpv1CPV2 &&
             cpv2CPV2 && ph1->IsDispOK() && ph2->IsDispOK();
    case 82:
      return r > r3cut && ph1->E() > 0.3 && ph2->E() > 0.3 && ph1->Module() == 3 && ph2->Module() == 3 &&
             cpvAndTrack1 && cpvAndTrack2;
    // D20
    case 83:
      return r > r4cut;
    case 84:
      return r > r4cut && ph1->IsDispOK() && ph2->IsDispOK();
    case 85:
      return r > r4cut && track1CPV && track2CPV;
    case 86:
      return r > r4cut && track1CPV2 && track2CPV2;
    case 87:
      return r > r4cut && track1CPV && track2CPV && ph1->IsDispOK() && ph2->IsDispOK();
    case 88:
      return r > r4cut && track1CPV2 && track2CPV2 && ph1->IsDispOK() && ph2->IsDispOK();
    case 89:
      return r > r4cut && ph1->Module() == 3 && ph2->Module() == 3;
    case 90:
      return r > r4cut && ph1->Module() == 3 && ph2->Module() == 3 && cpv1CPV && cpv2CPV;
    case 91:
      return r > r4cut && ph1->Module() == 3 && ph2->Module() == 3 && cpv1CPV2 && cpv2CPV2;
    case 92:
      return r > r4cut && ph1->Module() == 3 && ph2->Module() == 3 && cpv1CPV && cpv2CPV && ph1->IsDispOK() &&
             ph2->IsDispOK();
    case 93:
      return r > r4cut && ph1->Module() == 3 && ph2->Module() == 3 && cpv1CPV2 && cpv2CPV2 && ph1->IsDispOK() &&
             ph2->IsDispOK();
    case 94:
      return r > r4cut && ph1->Module() == 3 && ph2->Module() == 3 && cpvAndTrack1 && cpvAndTrack2;
    // Dz Emin 0.2
    case 95:
      return abs(ph1->Theta() - ph2->Theta()) > TMath::ATan2(6., 460.) && ph1->E() > 0.2 && ph2->E() > 0.2;
    case 96:
      return abs(ph1->Theta() - ph2->Theta()) > TMath::ATan2(6., 460.) && ph1->E() > 0.2 && ph2->E() > 0.2 &&
             ph1->IsDispOK() && ph2->IsDispOK();
    case 97:
      return abs(ph1->Theta() - ph2->Theta()) > TMath::ATan2(6., 460.) && ph1->E() > 0.2 && ph2->E() > 0.2 &&
             track1CPV && track2CPV;
    case 98:
      return abs(ph1->Theta() - ph2->Theta()) > TMath::ATan2(6., 460.) && ph1->E() > 0.2 && ph2->E() > 0.2 &&
             track1CPV2 && track2CPV2;
    case 99:
      return abs(ph1->Theta() - ph2->Theta()) > TMath::ATan2(6., 460.) && ph1->E() > 0.2 && ph2->E() > 0.2 &&
             track1CPV && track2CPV && ph1->IsDispOK() && ph2->IsDispOK();
    case 100:
      return abs(ph1->Theta() - ph2->Theta()) > TMath::ATan2(6., 460.) && ph1->E() > 0.2 && ph2->E() > 0.2 &&
             track1CPV2 && track2CPV2 && ph1->IsDispOK() && ph2->IsDispOK();
    case 101:
      return abs(ph1->Theta() - ph2->Theta()) > TMath::ATan2(6., 460.) && ph1->E() > 0.2 && ph2->E() > 0.2 &&
             ph1->Module() == 3 && ph2->Module() == 3;
    case 102:
      return abs(ph1->Theta() - ph2->Theta()) > TMath::ATan2(6., 460.) && ph1->E() > 0.2 && ph2->E() > 0.2 &&
             ph1->Module() == 3 && ph2->Module() == 3 && cpv1CPV && cpv2CPV;
    case 103:
      return abs(ph1->Theta() - ph2->Theta()) > TMath::ATan2(6., 460.) && ph1->E() > 0.2 && ph2->E() > 0.2 &&
             ph1->Module() == 3 && ph2->Module() == 3 && cpv1CPV2 && cpv2CPV2;
    case 104:
      return abs(ph1->Theta() - ph2->Theta()) > TMath::ATan2(6., 460.) && ph1->E() > 0.2 && ph2->E() > 0.2 &&
             ph1->Module() == 3 && ph2->Module() == 3 && cpv1CPV && cpv2CPV && ph1->IsDispOK() && ph2->IsDispOK();
    case 105:
      return abs(ph1->Theta() - ph2->Theta()) > TMath::ATan2(6., 460.) && ph1->E() > 0.2 && ph2->E() > 0.2 &&
             ph1->Module() == 3 && ph2->Module() == 3 && cpv1CPV2 && cpv2CPV2 && ph1->IsDispOK() && ph2->IsDispOK();
    case 106:
      return abs(ph1->Theta() - ph2->Theta()) > TMath::ATan2(6., 460.) && ph1->E() > 0.2 && ph2->E() > 0.2 &&
             ph1->Module() == 3 && ph2->Module() == 3 && cpvAndTrack1 && cpvAndTrack2;
    // Dz Emin 0.3
    case 107:
      return abs(ph1->Theta() - ph2->Theta()) > TMath::ATan2(6., 460.) && ph1->E() > 0.3 && ph2->E() > 0.3;
    case 108:
      return abs(ph1->Theta() - ph2->Theta()) > TMath::ATan2(6., 460.) && ph1->E() > 0.3 && ph2->E() > 0.3 &&
             ph1->IsDispOK() && ph2->IsDispOK();
    case 109:
      return abs(ph1->Theta() - ph2->Theta()) > TMath::ATan2(6., 460.) && ph1->E() > 0.3 && ph2->E() > 0.3 &&
             track1CPV && track2CPV;
    case 110:
      return abs(ph1->Theta() - ph2->Theta()) > TMath::ATan2(6., 460.) && ph1->E() > 0.3 && ph2->E() > 0.3 &&
             track1CPV2 && track2CPV2;
    case 111:
      return abs(ph1->Theta() - ph2->Theta()) > TMath::ATan2(6., 460.) && ph1->E() > 0.3 && ph2->E() > 0.3 &&
             track1CPV && track2CPV && ph1->IsDispOK() && ph2->IsDispOK();
    case 112:
      return abs(ph1->Theta() - ph2->Theta()) > TMath::ATan2(6., 460.) && ph1->E() > 0.3 && ph2->E() > 0.3 &&
             track1CPV2 && track2CPV2 && ph1->IsDispOK() && ph2->IsDispOK();
    case 113:
      return abs(ph1->Theta() - ph2->Theta()) > TMath::ATan2(6., 460.) && ph1->E() > 0.3 && ph2->E() > 0.3 &&
             ph1->Module() == 3 && ph2->Module() == 3;
    case 114:
      return abs(ph1->Theta() - ph2->Theta()) > TMath::ATan2(6., 460.) && ph1->E() > 0.3 && ph2->E() > 0.3 &&
             ph1->Module() == 3 && ph2->Module() == 3 && cpv1CPV && cpv2CPV;
    case 115:
      return abs(ph1->Theta() - ph2->Theta()) > TMath::ATan2(6., 460.) && ph1->E() > 0.3 && ph2->E() > 0.3 &&
             ph1->Module() == 3 && ph2->Module() == 3 && cpv1CPV2 && cpv2CPV2;
    case 116:
      return abs(ph1->Theta() - ph2->Theta()) > TMath::ATan2(6., 460.) && ph1->E() > 0.3 && ph2->E() > 0.3 &&
             ph1->Module() == 3 && ph2->Module() == 3 && cpv1CPV && cpv2CPV && ph1->IsDispOK() && ph2->IsDispOK();
    case 117:
      return abs(ph1->Theta() - ph2->Theta()) > TMath::ATan2(6., 460.) && ph1->E() > 0.3 && ph2->E() > 0.3 &&
             ph1->Module() == 3 && ph2->Module() == 3 && cpv1CPV2 && cpv2CPV2 && ph1->IsDispOK() && ph2->IsDispOK();
    case 118:
      return abs(ph1->Theta() - ph2->Theta()) > TMath::ATan2(6., 460.) && ph1->E() > 0.3 && ph2->E() > 0.3 &&
             ph1->Module() == 3 && ph2->Module() == 3 && cpvAndTrack1 && cpvAndTrack2;
    // D25
    case 119:
      return r > r5cut && ph1->E() > 0.3 && ph2->E() > 0.3;
    case 120:
      return r > r5cut && ph1->E() > 0.3 && ph2->E() > 0.3 && ph1->IsDispOK() && ph2->IsDispOK();
    case 121:
      return r > r5cut && ph1->E() > 0.3 && ph2->E() > 0.3 && track1CPV && track2CPV;
    case 122:
      return r > r5cut && ph1->E() > 0.3 && ph2->E() > 0.3 && track1CPV2 && track2CPV2;
    case 123:
      return r > r5cut && ph1->E() > 0.3 && ph2->E() > 0.3 && track1CPV && track2CPV && ph1->IsDispOK() &&
             ph2->IsDispOK();
    case 124:
      return r > r5cut && ph1->E() > 0.3 && ph2->E() > 0.3 && track1CPV2 && track2CPV2 && ph1->IsDispOK() &&
             ph2->IsDispOK();
    case 125:
      return r > r5cut && ph1->E() > 0.3 && ph2->E() > 0.3 && ph1->Module() == 3 && ph2->Module() == 3;
    case 126:
      return r > r5cut && ph1->E() > 0.3 && ph2->E() > 0.3 && ph1->Module() == 3 && ph2->Module() == 3 && cpv1CPV &&
             cpv2CPV;
    case 127:
      return r > r5cut && ph1->E() > 0.3 && ph2->E() > 0.3 && ph1->Module() == 3 && ph2->Module() == 3 && cpv1CPV2 &&
             cpv2CPV2;
    case 128:
      return r > r5cut && ph1->E() > 0.3 && ph2->E() > 0.3 && ph1->Module() == 3 && ph2->Module() == 3 && cpv1CPV &&
             cpv2CPV && ph1->IsDispOK() && ph2->IsDispOK();
    case 129:
      return r > r5cut && ph1->E() > 0.3 && ph2->E() > 0.3 && ph1->Module() == 3 && ph2->Module() == 3 && cpv1CPV2 &&
             cpv2CPV2 && ph1->IsDispOK() && ph2->IsDispOK();
    // NExMax
    // D10
    case 130:
      return r > r2cut && ph1->IsntUnfolded() && ph2->IsntUnfolded();
    case 131:
      return r > r2cut && ph1->IsntUnfolded() && ph2->IsntUnfolded() && ph1->IsDispOK() && ph2->IsDispOK();
    case 132:
      return r > r2cut && ph1->IsntUnfolded() && ph2->IsntUnfolded() && track1CPV && track2CPV;
    case 133:
      return r > r2cut && ph1->IsntUnfolded() && ph2->IsntUnfolded() && track1CPV2 && track2CPV2;
    case 134:
      return r > r2cut && ph1->IsntUnfolded() && ph2->IsntUnfolded() && track1CPV && track2CPV && ph1->IsDispOK() &&
             ph2->IsDispOK();
    case 135:
      return r > r2cut && ph1->IsntUnfolded() && ph2->IsntUnfolded() && track1CPV2 && track2CPV2 && ph1->IsDispOK() &&
             ph2->IsDispOK();
    case 136:
      return r > r2cut && ph1->IsntUnfolded() && ph2->IsntUnfolded() && ph1->Module() == 3 && ph2->Module() == 3;
    case 137:
      return r > r2cut && ph1->IsntUnfolded() && ph2->IsntUnfolded() && ph1->Module() == 3 && ph2->Module() == 3 &&
             cpv1CPV && cpv2CPV;
    case 138:
      return r > r2cut && ph1->IsntUnfolded() && ph2->IsntUnfolded() && ph1->Module() == 3 && ph2->Module() == 3 &&
             cpv1CPV2 && cpv2CPV2;
    case 139:
      return r > r2cut && ph1->IsntUnfolded() && ph2->IsntUnfolded() && ph1->Module() == 3 && ph2->Module() == 3 &&
             cpv1CPV && cpv2CPV && ph1->IsDispOK() && ph2->IsDispOK();
    case 140:
      return r > r2cut && ph1->IsntUnfolded() && ph2->IsntUnfolded() && ph1->Module() == 3 && ph2->Module() == 3 &&
             cpv1CPV2 && cpv2CPV2 && ph1->IsDispOK() && ph2->IsDispOK();
    case 141:
      return r > r2cut && ph1->IsntUnfolded() && ph2->IsntUnfolded() && ph1->Module() == 3 && ph2->Module() == 3 &&
             cpvAndTrack1 && cpvAndTrack2;
    // DZ
    case 142:
      return abs(ph1->Theta() - ph2->Theta()) > TMath::ATan2(6., 460.);
    case 143:
      return abs(ph1->Theta() - ph2->Theta()) > TMath::ATan2(6., 460.) && ph1->IsDispOK() && ph2->IsDispOK();
    case 144:
      return abs(ph1->Theta() - ph2->Theta()) > TMath::ATan2(6., 460.) && track1CPV && track2CPV;
    case 145:
      return abs(ph1->Theta() - ph2->Theta()) > TMath::ATan2(6., 460.) && track1CPV2 && track2CPV2;
    case 146:
      return abs(ph1->Theta() - ph2->Theta()) > TMath::ATan2(6., 460.) && track1CPV && track2CPV && ph1->IsDispOK() &&
             ph2->IsDispOK();
    case 147:
      return abs(ph1->Theta() - ph2->Theta()) > TMath::ATan2(6., 460.) && track1CPV2 && track2CPV2 && ph1->IsDispOK() &&
             ph2->IsDispOK();
    case 148:
      return abs(ph1->Theta() - ph2->Theta()) > TMath::ATan2(6., 460.) && ph1->Module() == 3 && ph2->Module() == 3;
    case 149:
      return abs(ph1->Theta() - ph2->Theta()) > TMath::ATan2(6., 460.) && ph1->Module() == 3 && ph2->Module() == 3 &&
             cpv1CPV && cpv2CPV;
    case 150:
      return abs(ph1->Theta() - ph2->Theta()) > TMath::ATan2(6., 460.) && ph1->Module() == 3 && ph2->Module() == 3 &&
             cpv1CPV2 && cpv2CPV2;
    case 151:
      return abs(ph1->Theta() - ph2->Theta()) > TMath::ATan2(6., 460.) && ph1->Module() == 3 && ph2->Module() == 3 &&
             cpv1CPV && cpv2CPV && ph1->IsDispOK() && ph2->IsDispOK();
    case 152:
      return abs(ph1->Theta() - ph2->Theta()) > TMath::ATan2(6., 460.) && ph1->Module() == 3 && ph2->Module() == 3 &&
             cpv1CPV2 && cpv2CPV2 && ph1->IsDispOK() && ph2->IsDispOK();
    case 153:
      return abs(ph1->Theta() - ph2->Theta()) > TMath::ATan2(6., 460.) && ph1->Module() == 3 && ph2->Module() == 3 &&
             cpvAndTrack1 && cpvAndTrack2;
    // DZ2
    case 154:
      return abs(ph1->Theta() - ph2->Theta()) > TMath::ATan2(9., 460.);
    case 155:
      return abs(ph1->Theta() - ph2->Theta()) > TMath::ATan2(9., 460.) && ph1->IsDispOK() && ph2->IsDispOK();
    case 156:
      return abs(ph1->Theta() - ph2->Theta()) > TMath::ATan2(9., 460.) && track1CPV && track2CPV;
    case 157:
      return abs(ph1->Theta() - ph2->Theta()) > TMath::ATan2(9., 460.) && track1CPV2 && track2CPV2;
    case 158:
      return abs(ph1->Theta() - ph2->Theta()) > TMath::ATan2(9., 460.) && track1CPV && track2CPV && ph1->IsDispOK() &&
             ph2->IsDispOK();
    case 159:
      return abs(ph1->Theta() - ph2->Theta()) > TMath::ATan2(9., 460.) && track1CPV2 && track2CPV2 && ph1->IsDispOK() &&
             ph2->IsDispOK();
    case 160:
      return abs(ph1->Theta() - ph2->Theta()) > TMath::ATan2(9., 460.) && ph1->Module() == 3 && ph2->Module() == 3;
    case 161:
      return abs(ph1->Theta() - ph2->Theta()) > TMath::ATan2(9., 460.) && ph1->Module() == 3 && ph2->Module() == 3 &&
             cpv1CPV && cpv2CPV;
    case 162:
      return abs(ph1->Theta() - ph2->Theta()) > TMath::ATan2(9., 460.) && ph1->Module() == 3 && ph2->Module() == 3 &&
             cpv1CPV2 && cpv2CPV2;
    case 163:
      return abs(ph1->Theta() - ph2->Theta()) > TMath::ATan2(9., 460.) && ph1->Module() == 3 && ph2->Module() == 3 &&
             cpv1CPV && cpv2CPV && ph1->IsDispOK() && ph2->IsDispOK();
    case 164:
      return abs(ph1->Theta() - ph2->Theta()) > TMath::ATan2(9., 460.) && ph1->Module() == 3 && ph2->Module() == 3 &&
             cpv1CPV2 && cpv2CPV2 && ph1->IsDispOK() && ph2->IsDispOK();
    case 165:
      return abs(ph1->Theta() - ph2->Theta()) > TMath::ATan2(9., 460.) && ph1->Module() == 3 && ph2->Module() == 3 &&
             cpvAndTrack1 && cpvAndTrack2;
    // CrDZ
    case 166:
      return abs(ph1->Theta() - ph2->Theta()) > TMath::ATan2(6., 460.) && !badPair;
    case 167:
      return abs(ph1->Theta() - ph2->Theta()) > TMath::ATan2(6., 460.) && !badPair && ph1->IsDispOK() &&
             ph2->IsDispOK();
    case 168:
      return abs(ph1->Theta() - ph2->Theta()) > TMath::ATan2(6., 460.) && !badPair && track1CPV && track2CPV;
    case 169:
      return abs(ph1->Theta() - ph2->Theta()) > TMath::ATan2(6., 460.) && !badPair && track1CPV2 && track2CPV2;
    case 170:
      return abs(ph1->Theta() - ph2->Theta()) > TMath::ATan2(6., 460.) && !badPair && track1CPV && track2CPV &&
             ph1->IsDispOK() && ph2->IsDispOK();
    case 171:
      return abs(ph1->Theta() - ph2->Theta()) > TMath::ATan2(6., 460.) && !badPair && track1CPV2 && track2CPV2 &&
             ph1->IsDispOK() && ph2->IsDispOK();
    case 172:
      return abs(ph1->Theta() - ph2->Theta()) > TMath::ATan2(6., 460.) && !badPair && ph1->Module() == 3 &&
             ph2->Module() == 3;
    case 173:
      return abs(ph1->Theta() - ph2->Theta()) > TMath::ATan2(6., 460.) && !badPair && ph1->Module() == 3 &&
             ph2->Module() == 3 && cpv1CPV && cpv2CPV;
    case 174:
      return abs(ph1->Theta() - ph2->Theta()) > TMath::ATan2(6., 460.) && !badPair && ph1->Module() == 3 &&
             ph2->Module() == 3 && cpv1CPV2 && cpv2CPV2;
    case 175:
      return abs(ph1->Theta() - ph2->Theta()) > TMath::ATan2(6., 460.) && !badPair && ph1->Module() == 3 &&
             ph2->Module() == 3 && cpv1CPV && cpv2CPV && ph1->IsDispOK() && ph2->IsDispOK();
    case 176:
      return abs(ph1->Theta() - ph2->Theta()) > TMath::ATan2(6., 460.) && !badPair && ph1->Module() == 3 &&
             ph2->Module() == 3 && cpv1CPV2 && cpv2CPV2 && ph1->IsDispOK() && ph2->IsDispOK();
    case 177:
      return abs(ph1->Theta() - ph2->Theta()) > TMath::ATan2(6., 460.) && !badPair && ph1->Module() == 3 &&
             ph2->Module() == 3 && cpvAndTrack1 && cpvAndTrack2;
    // CrDZ2
    case 178:
      return abs(ph1->Theta() - ph2->Theta()) > TMath::ATan2(9., 460.) && !badPair;
    case 179:
      return abs(ph1->Theta() - ph2->Theta()) > TMath::ATan2(9., 460.) && !badPair && ph1->IsDispOK() &&
             ph2->IsDispOK();
    case 180:
      return abs(ph1->Theta() - ph2->Theta()) > TMath::ATan2(9., 460.) && !badPair && track1CPV && track2CPV;
    case 181:
      return abs(ph1->Theta() - ph2->Theta()) > TMath::ATan2(9., 460.) && !badPair && track1CPV2 && track2CPV2;
    case 182:
      return abs(ph1->Theta() - ph2->Theta()) > TMath::ATan2(9., 460.) && !badPair && track1CPV && track2CPV &&
             ph1->IsDispOK() && ph2->IsDispOK();
    case 183:
      return abs(ph1->Theta() - ph2->Theta()) > TMath::ATan2(9., 460.) && !badPair && track1CPV2 && track2CPV2 &&
             ph1->IsDispOK() && ph2->IsDispOK();
    case 184:
      return abs(ph1->Theta() - ph2->Theta()) > TMath::ATan2(9., 460.) && !badPair && ph1->Module() == 3 &&
             ph2->Module() == 3;
    case 185:
      return abs(ph1->Theta() - ph2->Theta()) > TMath::ATan2(9., 460.) && !badPair && ph1->Module() == 3 &&
             ph2->Module() == 3 && cpv1CPV && cpv2CPV;
    case 186:
      return abs(ph1->Theta() - ph2->Theta()) > TMath::ATan2(9., 460.) && !badPair && ph1->Module() == 3 &&
             ph2->Module() == 3 && cpv1CPV2 && cpv2CPV2;
    case 187:
      return abs(ph1->Theta() - ph2->Theta()) > TMath::ATan2(9., 460.) && !badPair && ph1->Module() == 3 &&
             ph2->Module() == 3 && cpv1CPV && cpv2CPV && ph1->IsDispOK() && ph2->IsDispOK();
    case 188:
      return abs(ph1->Theta() - ph2->Theta()) > TMath::ATan2(9., 460.) && !badPair && ph1->Module() == 3 &&
             ph2->Module() == 3 && cpv1CPV2 && cpv2CPV2 && ph1->IsDispOK() && ph2->IsDispOK();
    case 189:
      return abs(ph1->Theta() - ph2->Theta()) > TMath::ATan2(9., 460.) && !badPair && ph1->Module() == 3 &&
             ph2->Module() == 3 && cpvAndTrack1 && cpvAndTrack2;

    default:
      return false;
  }
}
//_________________________________________________________________________
void AliAnalysisTaskEtaPhigg::SetCutNames()
{
  // Set names corresponding to cut index
  snprintf(fcut[0], 20, "All");
  snprintf(fcut[1], 20, "D5All");
  snprintf(fcut[2], 20, "D5Disp");
  snprintf(fcut[3], 20, "D5CTS");
  snprintf(fcut[4], 20, "D5Both");
  snprintf(fcut[5], 20, "D5Mod3");
  snprintf(fcut[6], 20, "D5CPV");
  snprintf(fcut[7], 20, "D5CPVDisp");
  snprintf(fcut[8], 20, "D5CPV2");
  snprintf(fcut[9], 20, "D5CPV2Disp");
  snprintf(fcut[10], 20, "D5CPVCTS");
  // D10
  snprintf(fcut[11], 20, "D10All");
  snprintf(fcut[12], 20, "D10Disp");
  snprintf(fcut[13], 20, "D10CTS");
  snprintf(fcut[14], 20, "D10CTS2");
  snprintf(fcut[15], 20, "D10CTSDisp");
  snprintf(fcut[16], 20, "D10CTS2Disp");
  snprintf(fcut[17], 20, "D10Mod3");
  snprintf(fcut[18], 20, "D10CPV");
  snprintf(fcut[19], 20, "D10CPV2");
  snprintf(fcut[20], 20, "D10CPVDisp");
  snprintf(fcut[21], 20, "D10CPV2Disp");
  snprintf(fcut[22], 20, "D10CPVCTS");
  // D10 Emin 0.2
  snprintf(fcut[23], 20, "D10E2All");
  snprintf(fcut[24], 20, "D10E2Disp");
  snprintf(fcut[25], 20, "D10E2CTS");
  snprintf(fcut[26], 20, "D10E2CTS2");
  snprintf(fcut[27], 20, "D10E2CTSDisp");
  snprintf(fcut[28], 20, "D10E2CTS2Disp");
  snprintf(fcut[29], 20, "D10E2Mod3");
  snprintf(fcut[30], 20, "D10E2CPV");
  snprintf(fcut[31], 20, "D10E2CPV2");
  snprintf(fcut[32], 20, "D10E2CPVDisp");
  snprintf(fcut[33], 20, "D10E2CPV2Disp");
  snprintf(fcut[34], 20, "D10E2CPVCTS");
  // D10 Emin 0.3
  snprintf(fcut[35], 20, "D10E3All");
  snprintf(fcut[36], 20, "D10E3Disp");
  snprintf(fcut[37], 20, "D10E3CTS");
  snprintf(fcut[38], 20, "D10E3CTS2");
  snprintf(fcut[39], 20, "D10E3CTSDisp");
  snprintf(fcut[40], 20, "D10E3CTS2Disp");
  snprintf(fcut[41], 20, "D10E3Mod3");
  snprintf(fcut[42], 20, "D10E3CPV");
  snprintf(fcut[43], 20, "D10E3CPV2");
  snprintf(fcut[44], 20, "D10E3CPVDisp");
  snprintf(fcut[45], 20, "D10E3CPV2Disp");
  snprintf(fcut[46], 20, "D10E3CPVCTS");
  // D15
  snprintf(fcut[47], 20, "D15All");
  snprintf(fcut[48], 20, "D15Disp");
  snprintf(fcut[49], 20, "D15CTS");
  snprintf(fcut[50], 20, "D15CTS2");
  snprintf(fcut[51], 20, "D15CTSDisp");
  snprintf(fcut[52], 20, "D15CTS2Disp");
  snprintf(fcut[53], 20, "D15Mod3");
  snprintf(fcut[54], 20, "D15CPV");
  snprintf(fcut[55], 20, "D15CPV2");
  snprintf(fcut[56], 20, "D15CPVDisp");
  snprintf(fcut[57], 20, "D15CPV2Disp");
  snprintf(fcut[58], 20, "D15CPVCTS");
  // D15 Emin 0.2
  snprintf(fcut[59], 20, "D15E2All");
  snprintf(fcut[60], 20, "D15E2Disp");
  snprintf(fcut[61], 20, "D15E2CTS");
  snprintf(fcut[62], 20, "D15E2CTS2");
  snprintf(fcut[63], 20, "D15E2CTSDisp");
  snprintf(fcut[64], 20, "D15E2CTS2Disp");
  snprintf(fcut[65], 20, "D15E2Mod3");
  snprintf(fcut[66], 20, "D15E2CPV");
  snprintf(fcut[67], 20, "D15E2CPV2");
  snprintf(fcut[68], 20, "D15E2CPVDisp");
  snprintf(fcut[69], 20, "D15E2CPV2Disp");
  snprintf(fcut[70], 20, "D15E2CPVCTS");
  // D15 Emin 0.3
  snprintf(fcut[71], 20, "D15E3All");
  snprintf(fcut[72], 20, "D15E3Disp");
  snprintf(fcut[73], 20, "D15E3CTS");
  snprintf(fcut[74], 20, "D15E3CTS2");
  snprintf(fcut[75], 20, "D15E3CTSDisp");
  snprintf(fcut[76], 20, "D15E3CTS2Disp");
  snprintf(fcut[77], 20, "D15E3Mod3");
  snprintf(fcut[78], 20, "D15E3CPV");
  snprintf(fcut[79], 20, "D15E3CPV2");
  snprintf(fcut[80], 20, "D15E3CPVDisp");
  snprintf(fcut[81], 20, "D15E3CPV2Disp");
  snprintf(fcut[82], 20, "D15E3CPVCTS");
  // D20
  snprintf(fcut[83], 20, "D20All");
  snprintf(fcut[84], 20, "D20Disp");
  snprintf(fcut[85], 20, "D20CTS");
  snprintf(fcut[86], 20, "D20CTS2");
  snprintf(fcut[87], 20, "D20CTSDisp");
  snprintf(fcut[88], 20, "D20CTS2Disp");
  snprintf(fcut[89], 20, "D20Mod3");
  snprintf(fcut[90], 20, "D20CPV");
  snprintf(fcut[91], 20, "D20CPV2");
  snprintf(fcut[92], 20, "D20CPVDisp");
  snprintf(fcut[93], 20, "D20CPV2Disp");
  snprintf(fcut[94], 20, "D20CPVCTS");
  // D20 Emin 0.2
  snprintf(fcut[95], 20, "DzE2All");
  snprintf(fcut[96], 20, "DzE2Disp");
  snprintf(fcut[97], 20, "DzE2CTS");
  snprintf(fcut[98], 20, "DzE2CTS2");
  snprintf(fcut[99], 20, "DzE2CTSDisp");
  snprintf(fcut[100], 20, "DzE2CTS2Disp");
  snprintf(fcut[101], 20, "DzE2Mod3");
  snprintf(fcut[102], 20, "DzE2CPV");
  snprintf(fcut[103], 20, "DzE2CPV2");
  snprintf(fcut[104], 20, "DzE2CPVDisp");
  snprintf(fcut[105], 20, "DzE2CPV2Disp");
  snprintf(fcut[106], 20, "DzE2CPVCTS");
  // D20 Emin 0.3
  snprintf(fcut[107], 20, "DzE3All");
  snprintf(fcut[108], 20, "DzE3Disp");
  snprintf(fcut[109], 20, "DzE3CTS");
  snprintf(fcut[110], 20, "DzE3CTS2");
  snprintf(fcut[111], 20, "DzE3CTSDisp");
  snprintf(fcut[112], 20, "DzE3CTS2Disp");
  snprintf(fcut[113], 20, "DzE3Mod3");
  snprintf(fcut[114], 20, "DzE3CPV");
  snprintf(fcut[115], 20, "DzE3CPV2");
  snprintf(fcut[116], 20, "DzE3CPVDisp");
  snprintf(fcut[117], 20, "DzE3CPV2Disp");
  snprintf(fcut[118], 20, "DzE3CPVCTS");
  // D25
  snprintf(fcut[119], 20, "D25E3All");
  snprintf(fcut[120], 20, "D25E3Disp");
  snprintf(fcut[121], 20, "D25E3CTS");
  snprintf(fcut[122], 20, "D25E3CTS2");
  snprintf(fcut[123], 20, "D25E3CTSDisp");
  snprintf(fcut[124], 20, "D25E3CTS2Disp");
  snprintf(fcut[125], 20, "D25E3Mod3");
  snprintf(fcut[126], 20, "D25E3CPV");
  snprintf(fcut[127], 20, "D25E3CPV2");
  snprintf(fcut[128], 20, "D25E3CPVDisp");
  snprintf(fcut[129], 20, "D25E3CPV2Disp");
  // ExMax
  snprintf(fcut[130], 20, "D10NoUnfAll");
  snprintf(fcut[131], 20, "D10NoUnfDisp");
  snprintf(fcut[132], 20, "D10NoUnfCTS");
  snprintf(fcut[133], 20, "D10NoUnfCTS2");
  snprintf(fcut[134], 20, "D10NoUnfCTSDisp");
  snprintf(fcut[135], 20, "D10NoUnfCTS2Disp");
  snprintf(fcut[136], 20, "D10NoUnfMod3");
  snprintf(fcut[137], 20, "D10NoUnfCPV");
  snprintf(fcut[138], 20, "D10NoUnfCPV2");
  snprintf(fcut[139], 20, "D10NoUnfCPVDisp");
  snprintf(fcut[140], 20, "D10NoUnfCPV2Disp");
  snprintf(fcut[141], 20, "D10NoUnfCPVCTS");
  // Dz
  snprintf(fcut[142], 20, "DzAll");
  snprintf(fcut[143], 20, "DzDisp");
  snprintf(fcut[144], 20, "DzCTS");
  snprintf(fcut[145], 20, "DzCTS2");
  snprintf(fcut[146], 20, "DzCTSDisp");
  snprintf(fcut[147], 20, "DzCTS2Disp");
  snprintf(fcut[148], 20, "DzMod3");
  snprintf(fcut[149], 20, "DzCPV");
  snprintf(fcut[150], 20, "DzCPV2");
  snprintf(fcut[151], 20, "DzCPVDisp");
  snprintf(fcut[152], 20, "DzCPV2Disp");
  snprintf(fcut[153], 20, "DzCPVCTS");
  // Dz2
  snprintf(fcut[154], 20, "Dz2All");
  snprintf(fcut[155], 20, "Dz2Disp");
  snprintf(fcut[156], 20, "Dz2CTS");
  snprintf(fcut[157], 20, "Dz2CTS2");
  snprintf(fcut[158], 20, "Dz2CTSDisp");
  snprintf(fcut[159], 20, "Dz2CTS2Disp");
  snprintf(fcut[160], 20, "Dz2Mod3");
  snprintf(fcut[161], 20, "Dz2CPV");
  snprintf(fcut[162], 20, "Dz2CPV2");
  snprintf(fcut[163], 20, "Dz2CPVDisp");
  snprintf(fcut[164], 20, "Dz2CPV2Disp");
  snprintf(fcut[165], 20, "Dz2CPVCTS");
  // CrDZ
  snprintf(fcut[166], 20, "CrDzAll");
  snprintf(fcut[167], 20, "CrDzDisp");
  snprintf(fcut[168], 20, "CrDzCTS");
  snprintf(fcut[169], 20, "CrDzCTS2");
  snprintf(fcut[170], 20, "CrDzCTSDisp");
  snprintf(fcut[171], 20, "CrDzCTS2Disp");
  snprintf(fcut[172], 20, "CrDzMod3");
  snprintf(fcut[173], 20, "CrDzCPV");
  snprintf(fcut[174], 20, "CrDzCPV2");
  snprintf(fcut[175], 20, "CrDzCPVDisp");
  snprintf(fcut[176], 20, "CrDzCPV2Disp");
  snprintf(fcut[177], 20, "CrDzCPVCTS");
  // CrDZ2
  snprintf(fcut[178], 20, "CrDz2All");
  snprintf(fcut[179], 20, "CrDz2Disp");
  snprintf(fcut[180], 20, "CrDz2CTS");
  snprintf(fcut[181], 20, "CrDz2CTS2");
  snprintf(fcut[182], 20, "CrDz2CTSDisp");
  snprintf(fcut[183], 20, "CrDz2CTS2Disp");
  snprintf(fcut[184], 20, "CrDz2Mod3");
  snprintf(fcut[185], 20, "CrDz2CPV");
  snprintf(fcut[186], 20, "CrDz2CPV2");
  snprintf(fcut[187], 20, "CrDz2CPVDisp");
  snprintf(fcut[188], 20, "CrDz2CPV2Disp");
  snprintf(fcut[189], 20, "CrDz2CPVCTS");
}

//_________________________________________________________________________
int AliAnalysisTaskEtaPhigg::TestCPV(int mod, double e, double xPHOS, double zPHOS, double dxPHOS, double dzPHOS,
                                     int itr)
{
  // Return true if neutral
  // CTS
  bool trackCPV = true, trackCPV2 = true;
  AliAODTrack* tr = nullptr;
  double pt = 0;
  if (itr >= 0) {
    tr = (AliAODTrack*)fEvent->GetTrack(itr);
    pt = tr->Pt();
    double sx = TMath::Min(5.2, 1.111 + 0.56 * TMath::Exp(-0.031 * pt * pt) + 4.8 / TMath::Power(pt + 0.61, 3));
    double sz = TMath::Min(3.3, 1.12 + 0.35 * TMath::Exp(-0.032 * pt * pt) + 0.75 / TMath::Power(pt + 0.24, 3));
    double meanX = 0, meanZ = 0.102;
    if (fMF < 0.) { // field --
      if (tr->Charge() > 0)
        meanX = TMath::Min(5.8, 0.42 + 0.70 * TMath::Exp(-0.015 * pt * pt) + 35.8 / TMath::Power(pt + 1.41, 3));
      else
        meanX = -TMath::Min(5.8, 0.17 + 0.64 * TMath::Exp(-0.019 * pt * pt) + 26.1 / TMath::Power(pt + 1.21, 3));
    } else { // Field ++
      if (tr->Charge() > 0)
        meanX = -TMath::Min(5.8, 0.58 + 0.68 * TMath::Exp(-0.027 * pt * pt) + 28.0 / TMath::Power(pt + 1.28, 3));
      else
        meanX = TMath::Min(5.8, 0.11 + 0.67 * TMath::Exp(-0.015 * pt * pt) + 29.9 / TMath::Power(pt + 1.29, 3));
    }
    double r2 = pow((dzPHOS - meanZ) / sz, 2) + pow((dxPHOS - meanX) / sx, 2);
    trackCPV = (r2 > 1.);
    trackCPV2 = (r2 > 4.);
  }

  bool cpvCPV = true, cpvCPV2 = true, cpvAndTrack = true;
  if (mod == 3) {
    // test CPV cuts //find closest CPV cluster
    double dMin = 999999.;
    double dxCPV = 99999., dzCPV = 9999.;
    int icpvTr = -2;
    for (Int_t j = 0; j < fCPVEvent->GetEntriesFast(); j++) {
      AliCaloPhoton* cpv = (AliCaloPhoton*)fCPVEvent->At(j);
      double ax = cpv->X() - xPHOS; // huck only for mod 3
      double az = cpv->Z() - zPHOS;
      double d = ax * ax + az * az;
      if (dMin > d) {
        dMin = d;
        dxCPV = ax;
        dzCPV = az;
        icpvTr = cpv->DistToBad(); // track ID stored in this field
      }
    }
    // CPV bit
    if (dxCPV * dxCPV + dzCPV * dzCPV < 25.) {
      double dxMax = 3.36783 / e - 11.5189 / TMath::Sqrt(e) + 3.08283;
      double dxMin = -4.01590 / e + 13.2841 / TMath::Sqrt(e) - 4.31161;
      dxMax = dxCPV - dxMax;
      dxMin = dxCPV - dxMin;
      double sigmaX = TMath::Max(1.5, -2.07124 / e + 6.69554 / TMath::Sqrt(e) - 1.70062);
      double sigmaZ = 1.24859122035152259;
      if (e > 0.4)
        sigmaZ = 5.58984e-01 * TMath::Exp(-e * e / 2.20543 / 2.20543) + 7.07696e-01;
      double r1 = dxMax * dxMax / sigmaX / sigmaX + dzCPV * dzCPV / sigmaZ / sigmaZ;
      double r2 = dxMin * dxMin / sigmaX / sigmaX + dzCPV * dzCPV / sigmaZ / sigmaZ;
      cpvCPV = (r1 > 2.5 * 2.5) && (r2 > 2.5 * 2.5);
    }

    // account conversion
    double mX = -1.07550e+01 * TMath::Exp(-e / 1.75354) + 4.16349e-01;
    double wX = 2.94202e+00 * TMath::Exp(-e / 5.52975) + 2.18075e+00;
    double mZ = 0.2;
    double wZ = 6.32082e-01 * TMath::Exp(-e / 1.04274e+00) + 1.71270e+00;

    cpvCPV2 = ((dxCPV - mX) * (dxCPV - mX) / (2. * wX * wX) + (dzCPV - mZ) * (dzCPV - mZ) / (2. * wZ * wZ) > 9.) &&
              ((dxCPV + mX) * (dxCPV + mX) / (2. * wX * wX) + (dzCPV - mZ) * (dzCPV - mZ) / (2. * wZ * wZ) > 9.);

    // correlated CPV-PHOS vs track
    if (itr == icpvTr) {
      // true if neutral
      float sx = 1.5;
      float sz = 1.2;
      cpvAndTrack =
        (dxPHOS - dxCPV) * (dxPHOS - dxCPV) / (sx * sx) + (dzPHOS - dzCPV) * (dzPHOS - dzCPV) / (sz * sz) > 9.;
    }
    cpvAndTrack &= cpvCPV2; // no conversion and strong track deviations
  }

  // CPV cuts
  int cpvBits = 0;
  cpvBits |= (trackCPV << 1);    // Bit 1: CTS 1 sigma
  cpvBits |= (trackCPV2 << 2);   // Bit 2: CTS 2 sigma
  cpvBits |= (cpvCPV << 3);      // Bit 3: CTS 2 sigma
  cpvBits |= (cpvCPV2 << 4);     // Bit 4: CTS 2 sigma + conversion
  cpvBits |= (cpvAndTrack << 5); // Bit 5: correlated CTS + cvp c

  return cpvBits;
}
//_________________________________________________________________________
bool AliAnalysisTaskEtaPhigg::IsGoodChannel(Int_t cellX, Int_t cellZ)
{
  if (fBadMap)
    if (fBadMap->GetBinContent(cellX, cellZ) == 1)
      return kFALSE;

  return kTRUE;
}
//___________________________________________________________________________________________________
Int_t AliAnalysisTaskEtaPhigg::FindTrackMatching(int det, TVector3& locPHOS, int mod, float& dxMin, float& dzMin)
{
  // Find closest CPV
  if (det == 0) {
    mod = 3;
  }

  dxMin = 999.;
  dzMin = 999.;
  float dMin = 999.;

  // Find track with closest extrapolation to cluster
  Double_t magF = fEvent->GetMagneticField();

  Double_t magSign = 1.0;
  if (magF < 0)
    magSign = -1.0;

  if (!TGeoGlobalMagField::Instance()->GetField()) {
    AliError("Margnetic filed was not initialized, use default");
    AliMagF* field = new AliMagF("Maps", "Maps", magSign, magSign, AliMagF::k5kG);
    TGeoGlobalMagField::Instance()->SetField(field);
  }

  // *** Start the matching
  Int_t nt = fEvent->GetNumberOfTracks();

  // Calculate actual distance to PHOS module
  const Double_t r = TMath::Abs(locPHOS.Y());         // Distance to center of  PHOS/CPV module
  const Double_t kYmax = 72. + 30.;                   // Size of the module (with some reserve) in phi direction
  const Double_t kZmax = 64. + 30.;                   // Size of the module (with some reserve) in z direction
  const Double_t kAlpha0 = 330. / 180. * TMath::Pi(); // First PHOS module angular direction
  const Double_t kAlpha = 20. / 180. * TMath::Pi();   // PHOS module angular size

  Double_t gposTrack[3];

  Double_t bz = ((AliMagF*)TGeoGlobalMagField::Instance()->GetField())->SolenoidField();
  bz = TMath::Sign(0.5 * kAlmost0Field, bz) + bz;

  Double_t b[3];
  Int_t itr = -1;
  AliAODTrack* aodTrack = 0x0;
  Double_t xyz[3] = { 0 }, pxpypz[3] = { 0 }, cv[21] = { 0 };

  //  TPC only tracks: filter bit 128 (1<<7)  - all other bits are reset
  //  Global hybrid tracks: filter bit 256 (1<<8)   - all other bits are reset
  //  Complementary hybrid tracks: filter bit 512 (1<<9) - all other bits are reset
  int filterbits[3] = { 128, 256, 512 };

  for (Int_t i = 0; i < nt; i++) {
    aodTrack = (AliAODTrack*)fEvent->GetTrack(i);
    if (!aodTrack->TestFilterBit(filterbits[0]) && !aodTrack->TestFilterBit(filterbits[1]) &&
        !aodTrack->TestFilterBit(filterbits[2]))
      continue;

    // Continue extrapolation from TPC outer surface
    AliExternalTrackParam outerParam;
    aodTrack->GetPxPyPz(pxpypz);
    aodTrack->GetXYZ(xyz);
    aodTrack->GetCovarianceXYZPxPyPz(cv);

    outerParam.Set(xyz, pxpypz, cv, aodTrack->Charge());

    Double_t z;
    if (!outerParam.GetZAt(r, bz, z)) {
      continue;
    }
    if (TMath::Abs(z) > kZmax) {
      continue; // Some tracks miss the PHOS in Z
    }
    // Direction to the current PHOS module
    Double_t phiMod = kAlpha0 - kAlpha * mod;
    if (!outerParam.RotateParamOnly(phiMod))
      continue; // RS use faster rotation if errors are not needed

    Double_t y; // Some tracks do not reach the PHOS
    if (!outerParam.GetYAt(r, bz, y))
      continue; //    because of the bending

    if (TMath::Abs(y) < kYmax) {
      outerParam.GetBxByBz(b);
      outerParam.PropagateToBxByBz(r, b); // Propagate to the matching module
      // outerParam.CorrectForMaterial(...); // Correct for the TOF material, if needed
      outerParam.GetXYZ(gposTrack);
      TVector3 globalPositionTr(gposTrack);
      TVector3 localPositionTr;
      fPHOSGeo->Global2Local(localPositionTr, globalPositionTr, mod);
      Double_t ddx = locPHOS.X() - localPositionTr.X();
      Double_t ddz = locPHOS.Z() - localPositionTr.Z();

      float x = TMath::Max(0.6, aodTrack->Pt());
      float mean = 0;
      if (aodTrack->Charge() > 0) { // Pos
        if (det == 1) {             // PHOS
          mean = 2.85224e-02 - exp(-(x * x - 4.57197 * x - 3.63253e+01) / (7.48345 + x * 8.14244));
        } else { // cpv
          mean = -2.93946e-01 - exp(-(x * x + 1.15256e+07 * x - 1.77233e+07) / (6.18025e+06 + x * 1.45154e+06));
        }
      } else { // Neg
        if (det == 1) {
          mean = 2.77034e+00 + exp(-(x * x + 2.83814e+07 * x - 1.03977e+08) / (2.80232e+07 + x * 3.80026e+06));
        } else { // cpv
          mean = 4.38793e-02 - exp(-(x * x + 2.34803e+01 * x - 4.88497e+01) / (1.24926e+01 + x * 1.23229e+01));
        }
      }
      ddx -= mean;
      float d = ddx * ddx + ddz * ddz;
      if (d < dMin) {
        dxMin = ddx;
        dzMin = ddz;
        dMin = d;
        itr = i;
      }
    }
  } // Scanned all tracks
  return itr;
}
