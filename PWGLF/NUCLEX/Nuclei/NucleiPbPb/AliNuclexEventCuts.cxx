#include "AliNuclexEventCuts.h"

#include <TH1D.h>
#include <TH1I.h>
#include <TH2D.h>
#include <TList.h>

#include <AliAnalysisManager.h>
#include <AliCentrality.h>
#include <AliInputEventHandler.h>
#include <AliMultSelection.h>
#include <AliVEventHandler.h>
#include <AliVMultiplicity.h>

const string AliNuclexEventCuts::fgkLabels[2] = {"raw","selected"};

ClassImp(AliNuclexEventCuts);

/// Standard constructor with null selection
///
AliNuclexEventCuts::AliNuclexEventCuts() : TNamed("AliNuclexEventCuts","AliNuclexEventCuts"),
  fRequireTrackVertex{false},
  fMinVtz{-1000.f},
  fMaxVtz{1000.f},
  fMaxDeltaSpdTrackAbsolute{1000.f},
  fMaxDeltaSpdTrackNsigmaSPD{1000.f},
  fMaxDeltaSpdTrackNsigmaTrack{20000.f},
  //fMaxDispersionSPDvertex{0.03f},
  fMaxResolutionSPDvertex{1000.f},
  fRejectDAQincomplete{false},
  fSPDpileupMinContributors{1000},
  fSPDpileupMinZdist{-1.f},
  fSPDpileupNsigmaZdist{-1.f},
  fSPDpileupNsigmaDiamXY{-1.f},
  fSPDpileupNsigmaDiamZ{-1.f},
  fCentralityFramework{0},
  fMinCentrality{-1000.f},
  fMaxCentrality{1000.f},
  fMaxDeltaEstimators{1000.f},
  fRequireExactTriggerMask{false},
  fTriggerMask{AliVEvent::kAny},
  fCentEstimators{"V0M","CL0"},
  fCentPercentiles{-1.f},
  fPrimaryVertex{nullptr}
{
}

void AliNuclexEventCuts::SetupLHC15o() {
  fRequireTrackVertex = true;
  fMinVtz = -10.f;
  fMaxVtz = 10.f;
  fMaxDeltaSpdTrackAbsolute = 0.2f;
  fMaxDeltaSpdTrackNsigmaSPD = 10.f;
  fMaxDeltaSpdTrackNsigmaTrack = 20.f;
  //fMaxDispersionSPDvertex = 0.03f;
  fMaxResolutionSPDvertex = 0.25f;

  fRejectDAQincomplete = true;

  fSPDpileupMinContributors = 5;
  fSPDpileupMinZdist = 0.8;
  fSPDpileupNsigmaZdist = 3.;
  fSPDpileupNsigmaDiamXY = 2.;
  fSPDpileupNsigmaDiamZ = 5.;

  fCentralityFramework = 2;
  fCentEstimators[0] = "V0M";
  fCentEstimators[1] = "CL0";
  fMinCentrality = 0.f;
  fMaxCentrality = 90.f;
  fMaxDeltaEstimators = 7.5f;

  fTriggerMask = AliVEvent::kINT7;

}

bool AliNuclexEventCuts::AcceptEvent(AliVEvent *ev, TList *qaList) {
  TH1I* fCutStats = nullptr;
  TH1D* fVtz[2] = {nullptr};
  TH1D* fDeltaTrackSPDvtz[2] = {nullptr};
  TH1D* fCentrality[2] = {nullptr};
  TH2D* fEstimCorrelation[2] = {nullptr};
  TH2D* fMultCentCorrelation[2] = {nullptr};
  if (qaList) {
    fCutStats = (TH1I*)qaList->FindObject("fCutStats");
    for (int iS = 0; iS < 2; ++iS) {
      fVtz[iS] = (TH1D*)qaList->FindObject(Form("Vtz_%s",fgkLabels[iS].data()));
      fDeltaTrackSPDvtz[iS] = (TH1D*)qaList->FindObject(Form("DeltaVtz_%s",fgkLabels[iS].data()));
      fCentrality[iS] = (TH1D*)qaList->FindObject(Form("Centrality_%s",fgkLabels[iS].data()));
      fEstimCorrelation[iS] = (TH2D*)qaList->FindObject(Form("EstimCorrelation_%s",fgkLabels[iS].data()));
      fMultCentCorrelation[iS] = (TH2D*)qaList->FindObject(Form("MultCentCorrelation_%s",fgkLabels[iS].data()));
    }
  }

  /// The first bin of the cut stats corresponds to the total number of events
  if (fCutStats != nullptr) fCutStats->Fill(0);

  /// Event selection flag, as soon as the event does not pass one cut this becomes false.
  int pass = true;

  /// Rejection of the DAQ incomplete events
  if (fRejectDAQincomplete && ev->IsIncompleteDAQ()) pass = false;
  else if (fCutStats != nullptr) fCutStats->Fill(1);

  /// Trigger mask
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* handl = (AliInputEventHandler*)mgr->GetInputEventHandler();
  auto selected_trigger = handl->IsEventSelected() & fTriggerMask;
  if ((selected_trigger != fTriggerMask && fRequireExactTriggerMask) || !selected_trigger) pass = false;
  else if (fCutStats != nullptr) fCutStats->Fill(2);

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
  if (((vtTrc->GetNContributors() < 2 && fRequireTrackVertex) || vtSPD->GetNContributors() < 1) || // Check if SPD vertex is there and (if required) check if Track vertex is present.
      (TMath::Abs(dz) > fMaxDeltaSpdTrackAbsolute || nsigTot > fMaxDeltaSpdTrackNsigmaSPD || nsigTrc > fMaxDeltaSpdTrackNsigmaTrack) || // discrepancy track-SPD vertex
      (vtx->GetZ() < fMinVtz || vtx->GetZ() > fMaxVtz) || // min-max limits of the vertex z
      (vtSPD->IsFromVertexerZ() && TMath::Sqrt(covSPD[5]) > fMaxResolutionSPDvertex)) // quality cut on vertexer SPD z
    pass = false;
  else if (fCutStats != nullptr) fCutStats->Fill(3);
  fPrimaryVertex = const_cast<AliVVertex*>(vtx);

  /// SPD pile-up rejection
  if (ev->IsPileupFromSPD(fSPDpileupMinContributors,fSPDpileupMinZdist,fSPDpileupNsigmaZdist,fSPDpileupNsigmaDiamXY,fSPDpileupNsigmaDiamZ)) pass = false; // pile-up
  else if (fCutStats != nullptr) fCutStats->Fill(4);

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
    else if (fCutStats != nullptr) fCutStats->Fill(5);
  } else if (fCutStats != nullptr) fCutStats->Fill(5);

  /// Filling the monitoring histograms (first iteration always filled, second iteration only for selected events.
  for (int befaft = 0; befaft < 2; ++befaft) {
    if (fCentrality[befaft] != nullptr) fCentrality[befaft]->Fill(fCentPercentiles[0]);
    if (fEstimCorrelation[befaft] != nullptr) fEstimCorrelation[befaft]->Fill(fCentPercentiles[0],fCentPercentiles[1]);
    if (fMultCentCorrelation[befaft] != nullptr) fMultCentCorrelation[befaft]->Fill(fCentPercentiles[0],ntrkl);
    if (fVtz[befaft] != nullptr) fVtz[befaft]->Fill(vtx->GetZ());
    if (fDeltaTrackSPDvtz[befaft] != nullptr) fDeltaTrackSPDvtz[befaft]->Fill(dz);
    if (!pass) return pass; /// Do not fill the "after" histograms if the event does not pass the cuts.
  }

  /// Last bin of the cut stats contains the number of selected events
  if (fCutStats != nullptr) fCutStats->Fill(6);

  return pass;
}

void AliNuclexEventCuts::AddQAplotsToList(TList *qaList) {
  TH1I* fCutStats = new TH1I("fCutStats",";;Number of selected events",7,-.5,6.5);
  string bin_labels[7] = {"No cuts","DAQ Incomplete","Trigger selection","Vertex selection","Pile-up","Centrality selection","All cuts"};
  for (int iB = 1; iB <= 7; ++iB) fCutStats->GetXaxis()->SetBinLabel(iB,bin_labels[iB-1].data());
  qaList->Add(fCutStats);

  TH1D* fVtz[2];                 //<! Vertex z distribution
  TH1D* fDeltaTrackSPDvtz[2];    //<! Difference between the vertex computed using SPD and the track vertex
  TH1D* fCentrality[2];          //<! Centrality percentile distribution
  TH2D* fEstimCorrelation[2];    //<! Correlation between centrality estimators
  TH2D* fMultCentCorrelation[2]; //<! Correlation between main centrality estimator and multiplicity
  string titles[2] = {"before event cuts","after event cuts"};
  for (int iS = 0; iS < 2; ++iS) {
    fVtz[iS] = new TH1D(Form("Vtz_%s",fgkLabels[iS].data()),Form("Vertex z %s; #it{v_{z}} (cm); Events",titles[iS].data()),400,-20.,20.);
    fDeltaTrackSPDvtz[iS] = new TH1D(Form("DeltaVtz_%s",fgkLabels[iS].data()),Form("Vertex tracks - Vertex SPD %s; #Delta#it{v_{z}} (cm); Events",titles[iS].data()),400,-2.,2.);
    fCentrality[iS] = new TH1D(Form("Centrality_%s",fgkLabels[iS].data()),Form("Centrality percentile %s; Centrality (%%); Events",titles[iS].data()),100,0.,100.);
    fEstimCorrelation[iS] = new TH2D(Form("EstimCorrelation_%s",fgkLabels[iS].data()),Form("Correlation estimators %s;%s;%s",titles[iS].data(),fCentEstimators[0].data(),fCentEstimators[1].data()),100,0.,100.,100,0.,100.);
    fMultCentCorrelation[iS] = new TH2D(Form("MultCentCorrelation_%s",fgkLabels[iS].data()),Form("Correlation multiplicity-centrality %s;Percentile of %s; Number of tracklets",titles[iS].data(),fCentEstimators[0].data()),100,0.,100.,2000,0.,10000.);
    
    qaList->Add(fVtz[iS]);
    qaList->Add(fDeltaTrackSPDvtz[iS]);
    qaList->Add(fCentrality[iS]);
    qaList->Add(fEstimCorrelation[iS]);
    qaList->Add(fMultCentCorrelation[iS]);
  }

}

float AliNuclexEventCuts::GetCentrality (unsigned int estimator) const { 
  if (estimator > 1) {
    /// Escalate to Fatal
    ::Fatal("AliNuclexEventCuts::GetCentrality","You asked for the centrality estimator with index %i, but you should choose between index 0 and 1.", estimator);
    return -1.f;
  } else 
    return fCentPercentiles[estimator];
}

string AliNuclexEventCuts::GetCentralityEstimator (unsigned int estimator) const { 
  if (estimator > 1) {
    /// Escalate to Fatal
    ::Fatal("AliNuclexEventCuts::GetCentralityEstimator","You asked for the centrality estimator with index %i, but you should choose between index 0 and 1.", estimator);
    return "";
  } else 
    return fCentEstimators[estimator];
}
