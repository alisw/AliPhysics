#include "AliNuclexEventCuts.h"

#include <TH1D.h>
#include <TH1I.h>
#include <TH2D.h>

#include <AliAnalysisManager.h>
#include <AliCentrality.h>
#include <AliInputEventHandler.h>
#include <AliMultSelection.h>
#include <AliVEventHandler.h>
#include <AliVMultiplicity.h>

ClassImp(AliNuclexEventCuts);

AliNuclexEventCuts::AliNuclexEventCuts() : TNamed("AliNuclexEventCuts","AliNuclexEventCuts"),
  fRequireTrackVertex{true},
  fMinVtz{-10.f},
  fMaxVtz{10.f},
  fMaxDeltaSpdTrackAbsolute{0.2f},
  fMaxDeltaSpdTrackNsigmaSPD{10.f},
  fMaxDeltaSpdTrackNsigmaTrack{20.f},
  fRejectDAQincomplete{true},
  fRejectPileupSPD{5},
  fCentralityFramework{2},
  fMinCentrality{0.f},
  fMaxCentrality{90.f},
  fCentEstimators{"V0M","CL0"},
  fCentPercentiles{-1.f},
  fMaxDeltaEstimators{7.5f},
  fTriggerMask{AliVEvent::kINT7},
  fCutStats{new TH1I("fCutStats",";;Number of selected events",7,-.5,6.5)},
  fVtz{0x0},
  fDeltaTrackSPDvtz{0x0},
  fCentrality{0x0},
  fEstimCorrelation{0x0},
  fMultCentCorrelation{0x0}
{
  string labels[2] = {"raw","selected"};
  string titles[2] = {"before event cuts","after event cuts"};
  for (int iS = 0; iS < 2; ++iS) {
    fVtz[iS] = new TH1D(Form("Vtz_%s",labels[iS].data()),Form("Vertex z %s; #it{v_{z}} (cm); Events",titles[iS].data()),400,-20.,20.);
    fDeltaTrackSPDvtz[iS] = new TH1D(Form("DeltaVtz_%s",labels[iS].data()),Form("Vertex tracks - Vertex SPD %s; #Delta#it{v_{z}} (cm); Events",titles[iS].data()),400,-2.,2.);
    fCentrality[iS] = new TH1D(Form("Centrality_%s",labels[iS].data()),Form("Centrality percentile %s; Centrality (%); Events",titles[iS].data()),100,0.,100.);
    fEstimCorrelation[iS] = new TH2D(Form("EstimCorrelation_%s",labels[iS].data()),Form("Correlation estimators %s",titles[iS].data()),100,0.,100.,100,0.,100.);
    fMultCentCorrelation[iS] = new TH2D(Form("MultCentCorrelation_%s",labels[iS].data()),Form("Correlation multiplicity-centrality %s;Centrality (%); Number of tracklets",titles[iS].data()),100,0.,100.,2000,0.,10000.);
  }

  string bin_labels[7] = {"No cuts","DAQ Incomplete","Trigger selection","Vertex selection","Pile-up","Centrality selection","All cuts"};
  for (int iB = 1; iB <= 7; ++iB) fCutStats->GetXaxis()->SetBinLabel(iB,bin_labels[iB-1].data());
}

AliNuclexEventCuts::AliNuclexEventCuts(const AliNuclexEventCuts& ev) : TNamed("AliNuclexEventCuts","AliNuclexEventCuts"),
  fRequireTrackVertex{ev.fRequireTrackVertex},
  fMinVtz{ev.fMinVtz},
  fMaxVtz{ev.fMaxVtz},
  fMaxDeltaSpdTrackAbsolute{ev.fMaxDeltaSpdTrackAbsolute},
  fMaxDeltaSpdTrackNsigmaSPD{ev.fMaxDeltaSpdTrackNsigmaSPD},
  fMaxDeltaSpdTrackNsigmaTrack{ev.fMaxDeltaSpdTrackNsigmaTrack},
  fRejectDAQincomplete{ev.fRejectDAQincomplete},
  fRejectPileupSPD{ev.fRejectPileupSPD},
  fCentralityFramework{ev.fCentralityFramework},
  fMinCentrality{ev.fMinCentrality},
  fMaxCentrality{ev.fMaxCentrality},
  fCentEstimators{ev.fCentEstimators[0],ev.fCentEstimators[1]},
  fCentPercentiles{ev.fCentPercentiles[0],ev.fCentPercentiles[1]},
  fMaxDeltaEstimators{ev.fMaxDeltaEstimators},
  fTriggerMask{ev.fTriggerMask},
  fCutStats{new TH1I(*ev.fCutStats)},
  fVtz{0x0},
  fDeltaTrackSPDvtz{0x0},
  fCentrality{0x0},
  fEstimCorrelation{0x0},
  fMultCentCorrelation{0x0}
{
  for (int iS = 0; iS < 2; ++iS) {
    fVtz[iS] = new TH1D(*ev.fVtz[iS]);
    fDeltaTrackSPDvtz[iS] = new TH1D(*ev.fDeltaTrackSPDvtz[iS]);
    fCentrality[iS] = new TH1D(*ev.fCentrality[iS]);
    fEstimCorrelation[iS] = new TH2D(*ev.fEstimCorrelation[iS]);
  }
}

AliNuclexEventCuts& AliNuclexEventCuts::operator=(const AliNuclexEventCuts& ev) {
  if (this != &ev) {
    fRequireTrackVertex = ev.fRequireTrackVertex;
    fMinVtz = ev.fMinVtz;
    fMaxVtz = ev.fMaxVtz;
    fMaxDeltaSpdTrackAbsolute = ev.fMaxDeltaSpdTrackAbsolute;
    fMaxDeltaSpdTrackNsigmaSPD = ev.fMaxDeltaSpdTrackNsigmaSPD;
    fMaxDeltaSpdTrackNsigmaTrack = ev.fMaxDeltaSpdTrackNsigmaTrack;
    fRejectDAQincomplete = ev.fRejectDAQincomplete;
    fRejectPileupSPD = ev.fRejectPileupSPD;
    fCentralityFramework = ev.fCentralityFramework;
    fMinCentrality = ev.fMinCentrality;
    fMaxCentrality = ev.fMaxCentrality;
    fCentEstimators[0] = ev.fCentEstimators[0];
    fCentEstimators[1] = ev.fCentEstimators[1];
    fCentPercentiles[0] = ev.fCentPercentiles[0];
    fCentPercentiles[1] = ev.fCentPercentiles[1];
    fMaxDeltaEstimators = ev.fMaxDeltaEstimators;
    fTriggerMask = ev.fTriggerMask;
    fCutStats = new TH1I(*ev.fCutStats);
    for (int iS = 0; iS < 2; ++iS) {
      fVtz[iS] = new TH1D(*ev.fVtz[iS]);
      fDeltaTrackSPDvtz[iS] = new TH1D(*ev.fDeltaTrackSPDvtz[iS]);
      fCentrality[iS] = new TH1D(*ev.fCentrality[iS]);
      fEstimCorrelation[iS] = new TH2D(*ev.fEstimCorrelation[iS]);
    }
  }
  return *this;
}

AliNuclexEventCuts::~AliNuclexEventCuts() {
  if (fCutStats) delete fCutStats;
  for (int iS = 0; iS < 2; ++iS) {
    if (fVtz[iS]) delete fVtz[iS];
    if (fDeltaTrackSPDvtz[iS]) delete fDeltaTrackSPDvtz[iS];
    if (fCentrality[iS]) delete fCentrality[iS];
    if (fEstimCorrelation[iS]) delete fEstimCorrelation[iS];      
  }
}

void AliNuclexEventCuts::SetupLHC15o() {
  fRequireTrackVertex = true;
  fMinVtz = -10.f;
  fMaxVtz = 10.f;
  fMaxDeltaSpdTrackAbsolute = 0.2f;
  fMaxDeltaSpdTrackNsigmaSPD = 10.f;
  fMaxDeltaSpdTrackNsigmaTrack = 20.f;

  fRejectDAQincomplete = true;

  fRejectPileupSPD = 5;

  fCentralityFramework = 2;
  fCentEstimators[0] = "V0M";
  fCentEstimators[1] = "CL0";
  fMinCentrality = 0.f;
  fMaxCentrality = 90.f;
  fMaxDeltaEstimators = 7.5f;

  fTriggerMask = AliVEvent::kINT7;

}

bool AliNuclexEventCuts::AcceptEvent(AliVEvent *ev) {
  /// The first bin of the cut stats corresponds to the total number of events
  fCutStats->Fill(0);

  /// Event selection flag, as soon as the event does not pass one cut this becomes false.
  int pass = true;

  /// Rejection of the DAQ incomplete events
  if (ev->IsIncompleteDAQ()) pass = false;
  else fCutStats->Fill(1);

  /// Trigger mask
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* handl = (AliInputEventHandler*)mgr->GetInputEventHandler();
  if (handl->IsEventSelected() & fTriggerMask != fTriggerMask) pass = false;
  else fCutStats->Fill(2);

  /// Vertex selection
  const AliVVertex* vtTrc = ev->GetPrimaryVertex();
  const AliVVertex* vtSPD = ev->GetPrimaryVertexSPD();
  const AliVVertex* &vtx = (vtTrc->GetNContributors() < 2) ? vtSPD : vtTrc;
  double covTrc[6],covSPD[6];
  vtTrc->GetCovarianceMatrix(covTrc);
  vtSPD->GetCovarianceMatrix(covSPD);
  double dz = vtTrc->GetZ() - vtSPD->GetZ();
  double errTot = TMath::Sqrt(covTrc[5]+covSPD[5]);
  double errTrc = TMath::Sqrt(covTrc[5]);
  double nsigTot = TMath::Abs(dz) / errTot, nsigTrc = TMath::Abs(dz) / errTrc;
  if ((vtTrc->GetNContributors() < 2 || (vtSPD->GetNContributors() < 1 && fRequireTrackVertex)) || // Check if SPD vertex is there and (if required) check if Track vertex is present.
      (TMath::Abs(dz) > fMaxDeltaSpdTrackAbsolute || nsigTot > fMaxDeltaSpdTrackNsigmaSPD || nsigTrc > fMaxDeltaSpdTrackNsigmaTrack) || // discrepancy track-SPD vertex
      (vtx->GetZ() < fMinVtz || vtx->GetZ() > fMaxVtz)) // min-max limits of the vertex z
    pass = false;
  else fCutStats->Fill(3);

  /// SPD pile-up rejection
  if (ev->IsPileupFromSPD(fRejectPileupSPD,0.8,3.,2.,5.)) pass = false; // pile-up
  else fCutStats->Fill(4);

  /// Centrality cuts:
  /// * Check for min and max centrality
  /// * Cross check correlation between two centrality estimators
  AliVMultiplicity* mult = ev->GetMultiplicity();
  const int ntrkl = mult->GetNumberOfTracklets();
  if (fCentralityFramework) {
    if (fCentralityFramework == 1) {
      AliCentrality* cent = ev->GetCentrality();
      fCentPercentiles[0] = cent->GetCentralityPercentile(fCentEstimators[0].data());
      fCentPercentiles[1] = cent->GetCentralityPercentile(fCentEstimators[1].data());
    } else {
      AliMultSelection* cent = (AliMultSelection*)ev->FindListObject("MultSelection");
      fCentPercentiles[0] = cent->GetMultiplicityPercentile(fCentEstimators[0].data(), true);
      fCentPercentiles[1] = cent->GetMultiplicityPercentile(fCentEstimators[1].data(), true);
    }
    fCentrality[0]->Fill(fCentPercentiles[0]);
    fEstimCorrelation[0]->Fill(fCentPercentiles[0],fCentPercentiles[1]);
    if (TMath::Abs(fCentPercentiles[1] - fCentPercentiles[0]) > fMaxDeltaEstimators 
        || fCentPercentiles[0] < fMinCentrality || fCentPercentiles[0] > fMaxCentrality) pass = false;
    else fCutStats->Fill(5);
  } else fCutStats->Fill(5);

  /// Filling the monitoring histograms (first iteration always filled, second iteration only for selected events.
  for (int befaft = 0; befaft < 2; ++befaft) {
    fCentrality[befaft]->Fill(fCentPercentiles[0]);
    fEstimCorrelation[befaft]->Fill(fCentPercentiles[0],fCentPercentiles[1]);
    fMultCentCorrelation[befaft]->Fill(fCentPercentiles[0],ntrkl);
    fVtz[befaft]->Fill(vtx->GetZ());
    fDeltaTrackSPDvtz[befaft]->Fill(dz);
    if (!pass) return pass; /// Do not fill the "after" histograms if the event does not pass the cuts.
  }

  /// Last bin of the cut stats contains the number of selected events
  fCutStats->Fill(6);

  return pass;
}
