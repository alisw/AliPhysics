/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
// Class for bad channels & bad runs identification
//
// This class is designed for the following tasks:
//   1. find bad (dead/noisy/strange) cells on a per run basis;
//   2. find the status/quality of analysed data (e.g. missing RCUs, run quality);
//   3. find the extent of problems related to bad cells: required for
//      systematics estimation for a physical quantity related to a cluster.
//
// See also AliAnalysisTaskCellsQA.
// Works for EMCAL and PHOS. Requires initialized geometry.
//
// Usage example for EMCAL:
//
// a) In LocalInit():
//
//    // to check the quality of analysed data and detect most of the problems:
//    AliEMCALGeometry::GetInstance("EMCAL_COMPLETEV1");
//    cellsQA = new AliCaloCellsQA(10, AliCaloCellsQA::kEMCAL); // 10 supermodules
//
//    // add this line for a full analysis (fills additional histograms):
//    cellsQA->ActivateFullAnalysis();
//
// b) In UserCreateOutputObjects():
//
//    // not required, but suggested
//    cellsQA->InitSummaryHistograms();
//
//
// c) In UserExec():
//
//    AliVEvent *event = InputEvent();
//    AliVCaloCells* cells = event->GetEMCALCells();
//
//    AliVVertex *vertex = (AliVVertex*) event->GetPrimaryVertex();
//    Double_t vertexXYZ[3];
//    vertex->GetXYZ(vertexXYZ);
//
//    // clusters
//    static TObjArray *clusArray = new TObjArray;
//    clusArray->Clear();
//    for (Int_t i = 0; i < event->GetNumberOfCaloClusters(); i++) {
//      AliVCluster *clus = event->GetCaloCluster(i);
//
//      // filter clusters here, if necessary
//      // if (...) continue;
//
//      clusArray->Add(clus);
//    }
//
//    // apply unfolding afterburner, etc.
//    // ...
//
//    cellsQA->Fill(event->GetRunNumber(), clusArray, cells, vertexXYZ);
//
// d) Do not forget to post data, where necessary:
//    PostData(1,cellsQA->GetListOfHistos());
//
//
// TODO: add DCAL (up to 6 supermodules?)
//
//----
//  Author: Olga Driga (SUBATECH)

// --- ROOT system ---
#include <TLorentzVector.h>

// --- AliRoot header files ---
#include <AliEMCALGeometry.h>
#include <AliPHOSGeometry.h>
#include <AliCaloCellsQA.h>

// ClassImp(AliCaloCellsQA)

//_________________________________________________________________________
AliCaloCellsQA::AliCaloCellsQA() :
  fDetector(0),
  fNMods(0),
  fClusElowMin(0),
  fClusEhighMin(0),
  fPi0EClusMin(0),
  fkFullAnalysis(0),
  fNBins_hECells(0),
  fNBins_hPi0Mass(0),
  fNBinsX_hNCellsInCluster(0),
  fNBinsY_hNCellsInCluster(0),
  fXMax_hECells(0),
  fXMax_hPi0Mass(0),
  fXMax_hNCellsInCluster(0),
  fAbsIdMin(0),
  fAbsIdMax(0),
  fListOfHistos(0),
  fhNEventsProcessedPerRun(0),
  fhCellLocMaxNTimesInClusterElow(),
  fhCellLocMaxNTimesInClusterEhigh(),
  fhCellLocMaxETotalClusterElow(),
  fhCellLocMaxETotalClusterEhigh(),
  fhCellNonLocMaxNTimesInClusterElow(),
  fhCellNonLocMaxNTimesInClusterEhigh(),
  fhCellNonLocMaxETotalClusterElow(),
  fhCellNonLocMaxETotalClusterEhigh(),
  fhECells(),
  fhPi0Mass(),
  fhNCellsInCluster(),
  fhCellAmplitude(0),
  fhCellAmplitudeEhigh(0),
  fhCellAmplitudeNonLocMax(0),
  fhCellAmplitudeEhighNonLocMax(0),
  fhCellTime(0)
{
  // Use this constructor if unsure.
  // Defaults: EMCAL, 10 supermodules, fClusElowMin = 0.3 GeV, fClusEhighMin = 1.0 GeV,
  //           fPi0EClusMin = 0.5 GeV, fkFullAnalysis = false.

  Init(10, kEMCAL, 100000, 200000);
}

//_________________________________________________________________________
AliCaloCellsQA::AliCaloCellsQA(Int_t nmods, Int_t det, Int_t startRunNumber, Int_t endRunNumber) :
  fDetector(0),
  fNMods(0),
  fClusElowMin(0),
  fClusEhighMin(0),
  fPi0EClusMin(0),
  fkFullAnalysis(0),
  fNBins_hECells(0),
  fNBins_hPi0Mass(0),
  fNBinsX_hNCellsInCluster(0),
  fNBinsY_hNCellsInCluster(0),
  fXMax_hECells(0),
  fXMax_hPi0Mass(0),
  fXMax_hNCellsInCluster(0),
  fAbsIdMin(0),
  fAbsIdMax(0),
  fListOfHistos(0),
  fhNEventsProcessedPerRun(0),
  fhCellLocMaxNTimesInClusterElow(),
  fhCellLocMaxNTimesInClusterEhigh(),
  fhCellLocMaxETotalClusterElow(),
  fhCellLocMaxETotalClusterEhigh(),
  fhCellNonLocMaxNTimesInClusterElow(),
  fhCellNonLocMaxNTimesInClusterEhigh(),
  fhCellNonLocMaxETotalClusterElow(),
  fhCellNonLocMaxETotalClusterEhigh(),
  fhECells(),
  fhPi0Mass(),
  fhNCellsInCluster(),
  fhCellAmplitude(0),
  fhCellAmplitudeEhigh(0),
  fhCellAmplitudeNonLocMax(0),
  fhCellAmplitudeEhighNonLocMax(0),
  fhCellTime(0)
{
  // Allows to set main class parameters.
  //
  // nmods -- maximum supermodule number + 1:
  //   use 4 for EMCAL <= 2010;
  //   use 4 for PHOS (PHOS numbers start from 1, not from zero);
  //   use 10 for EMCAL >= 2011.
  //
  // det -- detector: AliCaloCellsQA::kEMCAL or AliCaloCellsQA::kPHOS.
  //   Note: DCAL not yet implemented.
  //
  // startRunNumber, endRunNumber -- run range the analysis will run on
  //   (to allow later merging); wide (but reasonable) range can be provided.

  Init(nmods, det, startRunNumber, endRunNumber);
}

//_________________________________________________________________________
AliCaloCellsQA::~AliCaloCellsQA()
{
  delete fListOfHistos;
}

//_________________________________________________________________________
void AliCaloCellsQA::ActivateFullAnalysis()
{
  // Activate initialization and filling of all the designed histograms.
  // Must be called before InitSummaryHistograms() and Fill() methods, if necessary.

  fkFullAnalysis = kTRUE;
}

//_________________________________________________________________________
void AliCaloCellsQA::InitSummaryHistograms(Int_t nbins, Double_t emax,
                                           Int_t nbinsh, Double_t emaxh,
                                           Int_t nbinst, Double_t tmin, Double_t tmax)
{
  // Initialize summary (not per run) histograms.
  // Must be called before calling Fill() method, if necessary.
  //
  // nbins, emax -- number of bins and maximum amplitude for fhCellAmplitude and fhCellAmplitudeNonLocMax;
  // nbinsh, emaxh -- number of bins and maximum amplitude for fhCellAmplitudeEhigh and fhCellAmplitudeEhighNonLocMax;
  // nbinst, tmin, tmax -- number of bins and minimum/maximum time for fhCellTime.

  // do not add histograms to the current directory
  Bool_t ads = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  fhCellAmplitude = new TH2F("hCellAmplitude",
     "Amplitude distribution per cell", fAbsIdMax-fAbsIdMin,fAbsIdMin,fAbsIdMax, nbins,0,emax);
  fhCellAmplitude->SetXTitle("AbsId");
  fhCellAmplitude->SetYTitle("Amplitude, GeV");
  fhCellAmplitude->SetZTitle("Counts");

  fhCellAmplitudeEhigh = new TH2F("hCellAmplitudeEhigh",
     "Amplitude distribution per cell, high energies", fAbsIdMax-fAbsIdMin,fAbsIdMin,fAbsIdMax, nbinsh,0,emaxh);
  fhCellAmplitudeEhigh->SetXTitle("AbsId");
  fhCellAmplitudeEhigh->SetYTitle("Amplitude, GeV");
  fhCellAmplitudeEhigh->SetZTitle("Counts");

  fListOfHistos->Add(fhCellAmplitude);
  fListOfHistos->Add(fhCellAmplitudeEhigh);

  if (fkFullAnalysis) {
    fhCellAmplitudeNonLocMax = new TH2F("hCellAmplitudeNonLocMax",
      "Amplitude distribution per cell which is not a local maximum",
        fAbsIdMax-fAbsIdMin,fAbsIdMin,fAbsIdMax, nbins,0,emax);
    fhCellAmplitudeNonLocMax->SetXTitle("AbsId");
    fhCellAmplitudeNonLocMax->SetYTitle("Amplitude, GeV");
    fhCellAmplitudeNonLocMax->SetZTitle("Counts");

    fhCellAmplitudeEhighNonLocMax = new TH2F("hCellAmplitudeEhighNonLocMax",
      "Amplitude distribution per cell which is not a local maximum, high energies",
        fAbsIdMax-fAbsIdMin,fAbsIdMin,fAbsIdMax, nbinsh,0,emaxh);
    fhCellAmplitudeEhighNonLocMax->SetXTitle("AbsId");
    fhCellAmplitudeEhighNonLocMax->SetYTitle("Amplitude, GeV");
    fhCellAmplitudeEhighNonLocMax->SetZTitle("Counts");

    fhCellTime = new TH2F("hCellTime", "Time distribution per cell",
                          fAbsIdMax-fAbsIdMin,fAbsIdMin,fAbsIdMax, nbinst,tmin,tmax);
    fhCellTime->SetXTitle("AbsId");
    fhCellTime->SetYTitle("Time, microseconds");
    fhCellTime->SetZTitle("Counts");

    fListOfHistos->Add(fhCellAmplitudeNonLocMax);
    fListOfHistos->Add(fhCellAmplitudeEhighNonLocMax);
    fListOfHistos->Add(fhCellTime);
  }

  // return to the previous add directory status
  TH1::AddDirectory(ads);
}

//_________________________________________________________________________
void AliCaloCellsQA::Fill(Int_t runNumber, TObjArray *clusArray, AliVCaloCells *cells, Double_t vertexXYZ[3])
{
  // Does the job: fills histograms for current event.
  //
  // runNumber -- current run number;
  // clusArray -- array of clusters (TObjArray or TClonesArray);
  // cells -- EMCAL or PHOS cells;
  // vertexXYZ -- primary vertex position.

  fhNEventsProcessedPerRun->Fill(runNumber);
  Int_t ri = FindCurrentRunIndex(runNumber);

  FillCellsInCluster(ri, clusArray, cells);
  FillJustCells(ri, cells);
  FillPi0Mass(ri, clusArray, vertexXYZ);
}

//_________________________________________________________________________
void AliCaloCellsQA::SetClusterEnergyCuts(Double_t pi0EClusMin, Double_t elowMin, Double_t ehighMin)
{
  // Set cuts for minimum cluster energies.
  // Must be called before calling Fill() method, if necessary.

  fClusElowMin = elowMin;
  fClusEhighMin = ehighMin;
  fPi0EClusMin = pi0EClusMin;
}

//_________________________________________________________________________
void AliCaloCellsQA::SetBinningParameters(Int_t nbins1, Int_t nbins2, Int_t nbins3x, Int_t nbins3y,
                                          Double_t xmax1, Double_t xmax2, Double_t xmax3)
{
  // Set binning parameters for histograms hECells, hNCellsInCluster and hPi0Mass.
  // Must be called before Fill() method, if necessary.
  //
  // nbins1, xmax1 -- number of bins and maximum X axis value for hECells;
  // nbins2, xmax2 -- number of bins and maximum X axis value for hPi0Mass;
  // nbins3x, nbins3y, xmax3 -- number of bins in X and Y axes and maximum X axis value for hNCellsInCluster.

  fNBins_hECells = nbins1;
  fNBins_hPi0Mass = nbins2;
  fNBinsX_hNCellsInCluster = nbins3x;
  fNBinsY_hNCellsInCluster = nbins3y;
  fXMax_hECells = xmax1;
  fXMax_hPi0Mass = xmax2;
  fXMax_hNCellsInCluster = xmax3;
}

//_________________________________________________________________________
void AliCaloCellsQA::Init(Int_t nmods, Int_t det, Int_t startRunNumber, Int_t endRunNumber)
{
  // Class initialization.
  // Defaults: fClusElowMin = 0.3GeV, fClusEhighMin = 1GeV, fPi0EClusMin = 0.5GeV, fkFullAnalysis = false.

  // check input (for design limitations only)
  if (det != kEMCAL && det != kPHOS) {
    Error("AliCaloCellsQA::Init", "Wrong detector provided");
    Info("AliCaloCellsQA::Init", "I will use EMCAL");
    det = kEMCAL;
  }
  if (nmods < 1 || nmods > 10) {
    Error("AliCaloCellsQA::Init", "Wrong last supermodule number + 1 provided");
    Info("AliCaloCellsQA::Init", "I will use nmods = 10");
    nmods = 10;
  }

  fDetector = det;
  fNMods = nmods;
  fkFullAnalysis = kFALSE;
  SetClusterEnergyCuts();
  SetBinningParameters();

  // minimum/maximum cell absId;
  // straightforward solution avoids complications
  fAbsIdMin = 0;
  fAbsIdMax = fNMods * 1152;
  if (fDetector == kPHOS) {
    fAbsIdMin = 1;
    fAbsIdMax = 1 + (fNMods-1)*3584;
  }

  fListOfHistos = new TObjArray;
  fListOfHistos->SetOwner(kTRUE);

  fhNEventsProcessedPerRun = new TH1D("hNEventsProcessedPerRun",
        "Number of processed events vs run number", endRunNumber - startRunNumber, startRunNumber, endRunNumber);
  fhNEventsProcessedPerRun->SetDirectory(0);
  fhNEventsProcessedPerRun->SetXTitle("Run number");
  fhNEventsProcessedPerRun->SetYTitle("Events");
  fListOfHistos->Add(fhNEventsProcessedPerRun);

  // will indicate whether InitSummaryHistograms() was called
  fhCellAmplitude = NULL;

  // fill with NULLs (will indicate which supermodules can be filled)
  for (Int_t ri = 0; ri < 1000; ri++)
    for (Int_t sm = 0; sm < 10; sm++) {
      fhECells[ri][sm] = NULL;
      fhNCellsInCluster[ri][sm] = NULL;

      for (Int_t sm2 = 0; sm2 < 10; sm2++) fhPi0Mass[ri][sm][sm2] = NULL;
    }
}

//_________________________________________________________________________
Int_t AliCaloCellsQA::FindCurrentRunIndex(Int_t runNumber)
{
  // Return current run index; add a new run if necessary.

  static Int_t ri = -1;           // run index
  static Int_t runNumbers[1000];  // already encountered runs ...
  static Int_t nruns = 0;         // ... and their number

  // try previous value ...
  if (ri >= 0 && runNumbers[ri] == runNumber) return ri;

  // ... or find current run index ...
  for (ri = 0; ri < nruns; ri++)
    if (runNumbers[ri] == runNumber) return ri;

  // ... or add a new run
  if (nruns >= 1000)
    Fatal("AliCaloCellsQA::FindCurrentRunIndex", "Too many runs, how is this possible?");

  // ri = nruns
  runNumbers[ri] = runNumber;
  InitHistosForRun(runNumber, ri);
  nruns++;

  return ri;
}

//_________________________________________________________________________
void AliCaloCellsQA::InitHistosForRun(Int_t run, Int_t ri)
{
  // Initialize per run histograms for a new run number;
  // run -- run number; ri -- run index.

  // do not add histograms to the current directory
  Bool_t ads = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  fhCellLocMaxNTimesInClusterElow[ri] = new TH1F(Form("run%i_hCellLocMaxNTimesInClusterElow",run),
      "Number of times cell was local maximum in a low energy cluster", fAbsIdMax-fAbsIdMin,fAbsIdMin,fAbsIdMax);
  fhCellLocMaxNTimesInClusterElow[ri]->SetXTitle("AbsId");
  fhCellLocMaxNTimesInClusterElow[ri]->SetYTitle("Counts");

  fhCellLocMaxNTimesInClusterEhigh[ri] = new TH1F(Form("run%i_hCellLocMaxNTimesInClusterEhigh",run),
      "Number of times cell was local maximum in a high energy cluster", fAbsIdMax-fAbsIdMin,fAbsIdMin,fAbsIdMax);
  fhCellLocMaxNTimesInClusterEhigh[ri]->SetXTitle("AbsId");
  fhCellLocMaxNTimesInClusterEhigh[ri]->SetYTitle("Counts");

  fhCellLocMaxETotalClusterElow[ri] = new TH1F(Form("run%i_hCellLocMaxETotalClusterElow",run),
      "Total cluster energy for local maximum cell, low energy", fAbsIdMax-fAbsIdMin,fAbsIdMin,fAbsIdMax);
  fhCellLocMaxETotalClusterElow[ri]->SetXTitle("AbsId");
  fhCellLocMaxETotalClusterElow[ri]->SetYTitle("Energy");

  fhCellLocMaxETotalClusterEhigh[ri] = new TH1F(Form("run%i_hCellLocMaxETotalClusterEhigh",run),
      "Total cluster energy for local maximum cell, high energy", fAbsIdMax-fAbsIdMin,fAbsIdMin,fAbsIdMax);
  fhCellLocMaxETotalClusterEhigh[ri]->SetXTitle("AbsId");
  fhCellLocMaxETotalClusterEhigh[ri]->SetYTitle("Energy");

  fListOfHistos->Add(fhCellLocMaxNTimesInClusterElow[ri]);
  fListOfHistos->Add(fhCellLocMaxNTimesInClusterEhigh[ri]);
  fListOfHistos->Add(fhCellLocMaxETotalClusterElow[ri]);
  fListOfHistos->Add(fhCellLocMaxETotalClusterEhigh[ri]);


  if (fkFullAnalysis) {
    fhCellNonLocMaxNTimesInClusterElow[ri] = new TH1F(Form("run%i_hCellNonLocMaxNTimesInClusterElow",run),
        "Number of times cell wasn't local maximum in a low energy cluster", fAbsIdMax-fAbsIdMin,fAbsIdMin,fAbsIdMax);
    fhCellNonLocMaxNTimesInClusterElow[ri]->SetXTitle("AbsId");
    fhCellNonLocMaxNTimesInClusterElow[ri]->SetYTitle("Counts");

    fhCellNonLocMaxNTimesInClusterEhigh[ri] = new TH1F(Form("run%i_hCellNonLocMaxNTimesInClusterEhigh",run),
        "Number of times cell wasn't local maximum in a high energy cluster", fAbsIdMax-fAbsIdMin,fAbsIdMin,fAbsIdMax);
    fhCellNonLocMaxNTimesInClusterEhigh[ri]->SetXTitle("AbsId");
    fhCellNonLocMaxNTimesInClusterEhigh[ri]->SetYTitle("Counts");

    fhCellNonLocMaxETotalClusterElow[ri] = new TH1F(Form("run%i_hCellNonLocMaxETotalClusterElow",run),
        "Total cluster energy for not local maximum cell, low energy", fAbsIdMax-fAbsIdMin,fAbsIdMin,fAbsIdMax);
    fhCellNonLocMaxETotalClusterElow[ri]->SetXTitle("AbsId");
    fhCellNonLocMaxETotalClusterElow[ri]->SetYTitle("Energy");

    fhCellNonLocMaxETotalClusterEhigh[ri] = new TH1F(Form("run%i_hCellNonLocMaxETotalClusterEhigh",run),
        "Total cluster energy for not local maximum cell, high energy", fAbsIdMax-fAbsIdMin,fAbsIdMin,fAbsIdMax);
    fhCellNonLocMaxETotalClusterEhigh[ri]->SetXTitle("AbsId");
    fhCellNonLocMaxETotalClusterEhigh[ri]->SetYTitle("Energy");

    fListOfHistos->Add(fhCellNonLocMaxNTimesInClusterElow[ri]);
    fListOfHistos->Add(fhCellNonLocMaxNTimesInClusterEhigh[ri]);
    fListOfHistos->Add(fhCellNonLocMaxETotalClusterElow[ri]);
    fListOfHistos->Add(fhCellNonLocMaxETotalClusterEhigh[ri]);
  }


  Int_t minsm = 0;
  if (fDetector == kPHOS) minsm = 1;

  // per supermodule histograms
  for (Int_t sm = minsm; sm < fNMods; sm++) {
    fhECells[ri][sm] = new TH1F(Form("run%i_hECellsSM%i",run,sm),
      "Cell amplitude distribution", fNBins_hECells,0,fXMax_hECells);
    fhECells[ri][sm]->SetXTitle("Amplitude, GeV");
    fhECells[ri][sm]->SetYTitle("Number of cells");

    fhNCellsInCluster[ri][sm] = new TH2F(Form("run%i_hNCellsInClusterSM%i",run,sm),
      "Distrubution of number of cells in cluster vs cluster energy",
        fNBinsX_hNCellsInCluster,0,fXMax_hNCellsInCluster, fNBinsY_hNCellsInCluster,0,fNBinsY_hNCellsInCluster);
    fhNCellsInCluster[ri][sm]->SetXTitle("Energy, GeV");
    fhNCellsInCluster[ri][sm]->SetYTitle("Number of cells");
    fhNCellsInCluster[ri][sm]->SetZTitle("Counts");

    fListOfHistos->Add(fhECells[ri][sm]);
    fListOfHistos->Add(fhNCellsInCluster[ri][sm]);
  }

  // pi0 mass spectrum
  for (Int_t sm = minsm; sm < fNMods; sm++)
    for (Int_t sm2 = sm; sm2 < fNMods; sm2++) {
      fhPi0Mass[ri][sm][sm2] = new TH1F(Form("run%i_hPi0MassSM%iSM%i",run,sm,sm2),
                                        "#pi^{0} mass spectrum", fNBins_hPi0Mass,0,fXMax_hPi0Mass);
      fhPi0Mass[ri][sm][sm2]->SetXTitle("M_{#gamma#gamma}, GeV");
      fhPi0Mass[ri][sm][sm2]->SetYTitle("Counts");

      fListOfHistos->Add(fhPi0Mass[ri][sm][sm2]);
    }

  // return to the previous add directory status
  TH1::AddDirectory(ads);
}

//_________________________________________________________________________
void AliCaloCellsQA::FillCellsInCluster(Int_t ri, TObjArray *clusArray, AliVCaloCells *cells)
{
  // Fill histograms related to a cluster;
  // ri -- run index.

  Int_t sm;

  for (Int_t i = 0; i < clusArray->GetEntriesFast(); i++)
  {
    AliVCluster *clus = (AliVCluster*) clusArray->At(i);
    if ((sm = CheckClusterGetSM(clus)) < 0) continue;

    if (fhNCellsInCluster[ri][sm])
      fhNCellsInCluster[ri][sm]->Fill(clus->E(), clus->GetNCells());

    if (clus->E() >= fClusElowMin)
      for (Int_t c = 0; c < clus->GetNCells(); c++) {
        Int_t absId = clus->GetCellAbsId(c);

        if (IsCellLocalMaximum(c, clus, cells)) {// local maximum
          if (clus->E() < fClusEhighMin) {
            fhCellLocMaxNTimesInClusterElow[ri]->Fill(absId);
            fhCellLocMaxETotalClusterElow[ri]->Fill(absId, clus->E());
          } else {
            fhCellLocMaxNTimesInClusterEhigh[ri]->Fill(absId);
            fhCellLocMaxETotalClusterEhigh[ri]->Fill(absId, clus->E());
          }
        }
        else if (fkFullAnalysis) {// not a local maximum
          if (clus->E() < fClusEhighMin) {
            fhCellNonLocMaxNTimesInClusterElow[ri]->Fill(absId);
            fhCellNonLocMaxETotalClusterElow[ri]->Fill(absId, clus->E());
          } else {
            fhCellNonLocMaxNTimesInClusterEhigh[ri]->Fill(absId);
            fhCellNonLocMaxETotalClusterEhigh[ri]->Fill(absId, clus->E());
          }
        }
      } // cells loop

  } // cluster loop

}

//_________________________________________________________________________
void AliCaloCellsQA::FillJustCells(Int_t ri, AliVCaloCells *cells)
{
  // Fill cell histograms not related with a cluster;
  // ri -- run index.

  Short_t absId;
  Double_t amp, time;
  Int_t sm;

  for (Short_t c = 0; c < cells->GetNumberOfCells(); c++) {
    cells->GetCell(c, absId, amp, time);
    if ((sm = GetSM(absId)) < 0) continue;

    if (fhECells[ri][sm]) fhECells[ri][sm]->Fill(amp);

    if (fhCellAmplitude) {// in case InitSummaryHistograms() was not called
      fhCellAmplitude->Fill(absId, amp);
      fhCellAmplitudeEhigh->Fill(absId, amp);
      if (fkFullAnalysis) fhCellTime->Fill(absId, time);

      // fill not a local maximum distributions
      if (fkFullAnalysis && !IsCellLocalMaximum(absId, cells)) {
        fhCellAmplitudeNonLocMax->Fill(absId, amp);
        fhCellAmplitudeEhighNonLocMax->Fill(absId, amp);
      }
    }
  } // cell loop

}

//_________________________________________________________________________
void AliCaloCellsQA::FillPi0Mass(Int_t ri, TObjArray *clusArray, Double_t vertexXYZ[3])
{
  // Fill gamma+gamma invariant mass histograms.
  // ri -- run index.

  Int_t sm1, sm2;
  TLorentzVector p1, p2, psum;

  // cluster loop
  for (Int_t i = 0; i < clusArray->GetEntriesFast(); i++)
  {
    AliVCluster *clus = (AliVCluster*) clusArray->At(i);
    if (clus->E() < fPi0EClusMin) continue;
    if ((sm1 = CheckClusterGetSM(clus)) < 0) continue;

    clus->GetMomentum(p1, vertexXYZ);

    // second cluster loop
    for (Int_t j = i+1; j < clusArray->GetEntriesFast(); j++)
    {
      AliVCluster *clus2 = (AliVCluster*) clusArray->At(j);
      if (clus2->E() < fPi0EClusMin) continue;
      if ((sm2 = CheckClusterGetSM(clus2)) < 0) continue;

      clus2->GetMomentum(p2, vertexXYZ);

      psum = p1 + p2;
      if (psum.M2() < 0) continue;

      // s1 <= s2
      Int_t s1 = (sm1 <= sm2) ? sm1 : sm2;
      Int_t s2 = (sm1 <= sm2) ? sm2 : sm1;

      if (fhPi0Mass[ri][s1][s2])
        fhPi0Mass[ri][s1][s2]->Fill(psum.M());
    } // second cluster loop
  } // cluster loop
}

//_________________________________________________________________________
Int_t AliCaloCellsQA::CheckClusterGetSM(AliVCluster* clus)
{
  // Apply common cluster cuts and return supermodule number on success.
  // Return -1 if cuts not passed or an error occured.

  // check detector and number of cells in cluster
  if (fDetector == kEMCAL && !clus->IsEMCAL()) return -1;
  if (fDetector == kPHOS && !clus->IsPHOS()) return -1;
  if (clus->GetNCells() < 1) return -1;

  return GetSM(clus->GetCellAbsId(0));
}

//_________________________________________________________________________
Int_t AliCaloCellsQA::GetSM(Int_t absId)
{
  // Convert cell absId --> supermodule number. Return -1 in case of error.
  // Note: we use simple and straightforward way to find supermodule number.
  // This allows to avoid unnecessary external dependencies (on geometry).

  Int_t sm = -1;

  if (fDetector == kEMCAL) sm = absId/1152;
  if (fDetector == kPHOS) sm = 1 + (absId-1)/3584;

  // check for data corruption to avoid segfaults
  if (sm < 0 || sm > 9) {
    Error("AliCaloCellsQA::GetSM", "Data corrupted");
    return -1;
  }

  return sm;
}

//____________________________________________________________
Bool_t AliCaloCellsQA::IsCellLocalMaximum(Int_t c, AliVCluster* clus, AliVCaloCells* cells)
{
  // Returns true if cell is a local maximum in cluster (clusterizer dependent).
  // Cell fractions are taken into account.
  //
  // c -- cell number in the cluster.

  Int_t sm, eta, phi, sm2, eta2, phi2;

  Int_t absId = clus->GetCellAbsId(c);
  Double_t amp = cells->GetCellAmplitude(absId);
  if (clus->GetCellAmplitudeFraction(c) > 1e-5) amp *= clus->GetCellAmplitudeFraction(c);

  AbsIdToSMEtaPhi(absId, sm, eta, phi);

  // try to find a neighbour which has bigger amplitude
  for (Int_t c2 = 0; c2 < clus->GetNCells(); c2++) {
    if (c == c2) continue;

    Int_t absId2 = clus->GetCellAbsId(c2);
    Double_t amp2 = cells->GetCellAmplitude(absId2);
    if (clus->GetCellAmplitudeFraction(c2) > 1e-5) amp2 *= clus->GetCellAmplitudeFraction(c2);

    AbsIdToSMEtaPhi(absId2, sm2, eta2, phi2);
    if (sm != sm2) continue;

    Int_t deta = TMath::Abs(eta-eta2);
    Int_t dphi = TMath::Abs(phi-phi2);

    if ((deta == 0 && dphi == 1) || (deta == 1 && dphi == 0) || (deta == 1 && dphi == 1))
      if (amp < amp2) return kFALSE;
  } // cell loop

  return kTRUE;
}

//____________________________________________________________
Bool_t AliCaloCellsQA::IsCellLocalMaximum(Int_t absId, AliVCaloCells* cells)
{
  // Returns true if cell is a local maximum among cells (clusterizer independent).

  Int_t sm, eta, phi, sm2, eta2, phi2;

  Double_t amp = cells->GetCellAmplitude(absId);
  AbsIdToSMEtaPhi(absId, sm, eta, phi);

  // try to find a neighbour which has bigger amplitude
  for (Short_t c2 = 0; c2 < cells->GetNumberOfCells(); c2++) {
    Short_t absId2 = cells->GetCellNumber(c2);
    if (absId2 == absId) continue;

    AbsIdToSMEtaPhi(absId2, sm2, eta2, phi2);
    if (sm != sm2) continue;

    Int_t deta = TMath::Abs(eta-eta2);
    Int_t dphi = TMath::Abs(phi-phi2);

    if ((deta == 0 && dphi == 1) || (deta == 1 && dphi == 0) || (deta == 1 && dphi == 1))
      if (amp < cells->GetAmplitude(c2))
        return kFALSE;
  } // cell loop

  return kTRUE;
}

//____________________________________________________________
void AliCaloCellsQA::AbsIdToSMEtaPhi(Int_t absId, Int_t &sm, Int_t &eta, Int_t &phi)
{
  // Converts absId --> (sm, eta, phi) for a cell.
  // Works both for EMCAL and for PHOS.
  // Geometry must be already initialized.

  // EMCAL
  if (fDetector == kEMCAL) {
    static AliEMCALGeometry *geomEMCAL = AliEMCALGeometry::GetInstance();
    if (!geomEMCAL)
      Fatal("AliCaloCellsQA::AbsIdToSMEtaPhi", "EMCAL geometry is not initialized");

    Int_t nModule, nIphi, nIeta;
    geomEMCAL->GetCellIndex(absId, sm, nModule, nIphi, nIeta);
    geomEMCAL->GetCellPhiEtaIndexInSModule(sm, nModule, nIphi, nIeta, phi, eta);
    return;
  }

  // PHOS
  if (fDetector == kPHOS) {
    static AliPHOSGeometry *geomPHOS = AliPHOSGeometry::GetInstance();
    if (!geomPHOS)
      Fatal("AliCaloCellsQA::AbsIdToSMEtaPhi", "PHOS geometry is not initialized");

    Int_t relid[4];
    geomPHOS->AbsToRelNumbering(absId, relid);
    sm = relid[0];
    eta = relid[2];
    phi = relid[3];
  }

  // DCAL
  // not implemented
}
