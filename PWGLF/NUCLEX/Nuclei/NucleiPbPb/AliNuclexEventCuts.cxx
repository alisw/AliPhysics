#include "AliNuclexEventCuts.h"

#include <vector>
using std::vector;

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

/// Standard constructor with null selection
///
AliNuclexEventCuts::AliNuclexEventCuts(bool saveplots) : TList(),
  fRequireTrackVertex{false},
  fMinVtz{-1000.f},
  fMaxVtz{1000.f},
  fMaxDeltaSpdTrackAbsolute{1000.f},
  fMaxDeltaSpdTrackNsigmaSPD{1000.f},
  fMaxDeltaSpdTrackNsigmaTrack{20000.f},
  fMaxResolutionSPDvertex{1000.f},
  fRejectDAQincomplete{false},
  fRequiredSolenoidPolarity{0},
  fSPDpileupMinContributors{1000},
  fSPDpileupMinZdist{-1.f},
  fSPDpileupNsigmaZdist{-1.f},
  fSPDpileupNsigmaDiamXY{-1.f},
  fSPDpileupNsigmaDiamZ{-1.f},
  fTrackletBGcut{false},
  fCentralityFramework{0},
  fMinCentrality{-1000.f},
  fMaxCentrality{1000.f},
  fUseEstimatorsCorrelationCut{false},
  fEstimatorsCorrelationCoef{0.,1.},
  fEstimatorsSigmaPars{10000.,0.,0.,0.},
  fDeltaEstimatorNsigma{1.,1.},
  fRequireExactTriggerMask{false},
  fTriggerMask{AliVEvent::kAny},
  fkLabels{"raw","selected"},
  fManualMode{false},
  fSavePlots{saveplots},
  fCurrentRun{-1},
  fFlag{BIT(kNoCuts)},
  fCentEstimators{"V0M","CL0"},
  fCentPercentiles{-1.f},
  fPrimaryVertex{nullptr},
  fCutStats{nullptr},
  fVtz{nullptr},
  fDeltaTrackSPDvtz{nullptr},
  fCentrality{nullptr},
  fEstimCorrelation{nullptr},   
  fMultCentCorrelation{nullptr}
{
  SetName("AliNuclexEventCuts");
  SetOwner(true);
}

bool AliNuclexEventCuts::AcceptEvent(AliVEvent *ev) {
  /// If not specified the cuts are set according to the run period 
  const int current_run = ev->GetRunNumber();
  if (!fManualMode && current_run != fCurrentRun) {
    ::Info("AliNuclexEventCuts::AcceptEvent","Current run (%i) is different from the previous (%i): setting automatically the corresponding event cuts.",current_run,fCurrentRun);
    fCurrentRun = current_run;
    AutomaticSetup();
  }

  if (fSavePlots && !this->Last()) {
    AddQAplotsToList();
  }

  /// Event selection flag, as soon as the event does not pass one cut this becomes false.
  fFlag = BIT(kNoCuts);

  /// Rejection of the DAQ incomplete events
  if (!fRejectDAQincomplete || !ev->IsIncompleteDAQ()) fFlag |= BIT(kDAQincomplete); 

  /// Magnetic field selection
  float bField = ev->GetMagneticField();
  if (fRequiredSolenoidPolarity == 0 || fRequiredSolenoidPolarity * bField > 0.) fFlag |= BIT(kBfield);

  /// Trigger mask
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* handl = (AliInputEventHandler*)mgr->GetInputEventHandler();
  auto selected_trigger = handl->IsEventSelected() & fTriggerMask;
  if ((selected_trigger == fTriggerMask && fRequireExactTriggerMask) || (selected_trigger && !fRequireExactTriggerMask))
    fFlag |= BIT(kTrigger);

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

  /// Vertex position cut
  if (vtx->GetZ() >= fMinVtz && vtx->GetZ() <= fMaxVtz) fFlag |= BIT(kVertexPosition);

  /// Vertex quality cuts
  if (((vtTrc->GetNContributors() >= 2 ||  !fRequireTrackVertex) && vtSPD->GetNContributors() >= 1) && // Check if SPD vertex is there and (if required) check if Track vertex is present.
      (TMath::Abs(dz) <= fMaxDeltaSpdTrackAbsolute && nsigTot <= fMaxDeltaSpdTrackNsigmaSPD && nsigTrc <= fMaxDeltaSpdTrackNsigmaTrack) && // discrepancy track-SPD vertex
      (!vtSPD->IsFromVertexerZ() || TMath::Sqrt(covSPD[5]) <= fMaxResolutionSPDvertex)
     ) // quality cut on vertexer SPD z
    fFlag |= BIT(kVertexQuality);
  fPrimaryVertex = const_cast<AliVVertex*>(vtx);

  /// SPD pile-up rejection
  if (!ev->IsPileupFromSPD(fSPDpileupMinContributors,fSPDpileupMinZdist,fSPDpileupNsigmaZdist,fSPDpileupNsigmaDiamXY,fSPDpileupNsigmaDiamZ) && 
      (!fTrackletBGcut || !fUtils.IsSPDClusterVsTrackletBG(ev)))
    fFlag |= BIT(kPileUp);

  /// Centrality cuts:
  /// * Check for min and max centrality
  /// * Cross check correlation between two centrality estimators
  AliVMultiplicity* mult = ev->GetMultiplicity();
  const int ntrkl = mult->GetNumberOfTracklets();
  if (fCentralityFramework) {
    if (fCentralityFramework == 2) {
      AliCentrality* cent = ev->GetCentrality();
      fCentPercentiles[0] = cent->GetCentralityPercentile(fCentEstimators[0].data());
      fCentPercentiles[1] = cent->GetCentralityPercentile(fCentEstimators[1].data());
    } else {
      AliMultSelection* cent = (AliMultSelection*)ev->FindListObject("MultSelection");
      fCentPercentiles[0] = cent->GetMultiplicityPercentile(fCentEstimators[0].data(), true);
      fCentPercentiles[1] = cent->GetMultiplicityPercentile(fCentEstimators[1].data(), true);
    }
    const auto& x = fCentPercentiles[1];
    const double center = x * fEstimatorsCorrelationCoef[1] + fEstimatorsCorrelationCoef[0];
    const double sigma = fEstimatorsSigmaPars[0] + fEstimatorsSigmaPars[1] * x + fEstimatorsSigmaPars[2] * x * x + fEstimatorsSigmaPars[3] * x * x * x;
    fCentrality[0]->Fill(fCentPercentiles[0]);
    fEstimCorrelation[0]->Fill(fCentPercentiles[1],fCentPercentiles[0]);
    if ((!fUseEstimatorsCorrelationCut || 
          (fCentPercentiles[0] >= center - fDeltaEstimatorNsigma[0] * sigma
            && fCentPercentiles[0] <= center + fDeltaEstimatorNsigma[1] * sigma))
        && fCentPercentiles[0] >= fMinCentrality 
        && fCentPercentiles[0] <= fMaxCentrality) fFlag |= BIT(kMultiplicity);
  } else fFlag |= BIT(kMultiplicity);

  bool allcuts = (fFlag == (BIT(kAllCuts) - 1));
  if (allcuts) fFlag |= BIT(kAllCuts);
  for (int iCut = kNoCuts; iCut <= kAllCuts; ++iCut) {
    if (TESTBIT(fFlag,iCut)) 
      fCutStats->Fill(iCut);
  }

  /// Filling the monitoring histograms (first iteration always filled, second iteration only for selected events.
  for (int befaft = 0; befaft < 2; ++befaft) {
    if (fCentrality[befaft] != nullptr) fCentrality[befaft]->Fill(fCentPercentiles[0]);
    if (fEstimCorrelation[befaft] != nullptr) fEstimCorrelation[befaft]->Fill(fCentPercentiles[1],fCentPercentiles[0]);
    if (fMultCentCorrelation[befaft] != nullptr) fMultCentCorrelation[befaft]->Fill(fCentPercentiles[0],ntrkl);
    if (fVtz[befaft] != nullptr) fVtz[befaft]->Fill(vtx->GetZ());
    if (fDeltaTrackSPDvtz[befaft] != nullptr) fDeltaTrackSPDvtz[befaft]->Fill(dz);
    if (!allcuts) return false; /// Do not fill the "after" histograms if the event does not pass the cuts.
  }

  return true;
}

void AliNuclexEventCuts::AddQAplotsToList(TList *qaList) {
  if (!qaList) {
    if (fSavePlots)
      qaList = static_cast<TList*>(this);
    else
      ::Fatal("AliNuclexEventCuts::AddQAplotsToList","Invalid TList pointer provided");
  }
  static_cast<TList*>(this)->Delete();

  vector<string> bin_labels = {"No cuts","DAQ Incomplete","Magnetic field choice","Trigger selection","Vertex position","Vertex quality","Pile-up","Centrality selection","All cuts"};
  fCutStats = new TH1I("fCutStats",";;Number of selected events",bin_labels.size(),-.5,bin_labels.size() - 0.5);
  for (int iB = 1; iB <= bin_labels.size(); ++iB) fCutStats->GetXaxis()->SetBinLabel(iB,bin_labels[iB-1].data());
  qaList->Add(fCutStats);

  string titles[2] = {"before event cuts","after event cuts"};
  for (int iS = 0; iS < 2; ++iS) {
    fVtz[iS] = new TH1D(Form("Vtz_%s",fkLabels[iS].data()),Form("Vertex z %s; #it{v_{z}} (cm); Events",titles[iS].data()),400,-20.,20.);
    fDeltaTrackSPDvtz[iS] = new TH1D(Form("DeltaVtz_%s",fkLabels[iS].data()),Form("Vertex tracks - Vertex SPD %s; #Delta#it{v_{z}} (cm); Events",titles[iS].data()),400,-2.,2.);
    fCentrality[iS] = new TH1D(Form("Centrality_%s",fkLabels[iS].data()),Form("Centrality percentile %s; Centrality (%%); Events",titles[iS].data()),100,0.,100.);
    fEstimCorrelation[iS] = new TH2D(Form("EstimCorrelation_%s",fkLabels[iS].data()),Form("Correlation estimators %s;%s;%s",titles[iS].data(),fCentEstimators[1].data(),fCentEstimators[0].data()),100,0.,100.,100,0.,100.);
    fMultCentCorrelation[iS] = new TH2D(Form("MultCentCorrelation_%s",fkLabels[iS].data()),Form("Correlation multiplicity-centrality %s;Percentile of %s; Number of tracklets",titles[iS].data(),fCentEstimators[0].data()),100,0.,100.,2000,0.,10000.);

    qaList->Add(fVtz[iS]);
    qaList->Add(fDeltaTrackSPDvtz[iS]);
    qaList->Add(fCentrality[iS]);
    qaList->Add(fEstimCorrelation[iS]);
    qaList->Add(fMultCentCorrelation[iS]);
  }

}

void AliNuclexEventCuts::AutomaticSetup() {
  if ((fCurrentRun >= 225000 && fCurrentRun <= 244628) ||
      (fCurrentRun >= 256146 && fCurrentRun <= 260187)) {
    SetupRun2pp();
    return;
  }

  /// Run 2 Pb-Pb
  if ( fCurrentRun >= 244917 && fCurrentRun <= 256145 ) {
    SetupLHC15o();
    return;
  }

  ::Fatal("AliNuclexEventCuts::AutomaticSetup","Automatic period detection failed, please use the manual mode.");
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

void AliNuclexEventCuts::SetupRun2pp() {
  ::Info("AliNuclexEventCuts::SetupRun2pp","Setup event cuts for the Run2 pp periods.");
  SetName("StandardRun2ppEventCuts");

  fRequireTrackVertex = true;
  fMinVtz = -10.f;
  fMaxVtz = 10.f;
  fMaxDeltaSpdTrackAbsolute = 0.5f;
  fMaxDeltaSpdTrackNsigmaSPD = 2000.f;
  fMaxDeltaSpdTrackNsigmaTrack = 2000.f;
  fMaxResolutionSPDvertex = 0.25f;

  fRejectDAQincomplete = true;

  fRequiredSolenoidPolarity = 0;

  fSPDpileupMinContributors = 3;
  fSPDpileupMinZdist = 0.8;
  fSPDpileupNsigmaZdist = 3.;
  fSPDpileupNsigmaDiamXY = 2.;
  fSPDpileupNsigmaDiamZ = 5.;
  fTrackletBGcut = true;

  if (fCentralityFramework > 1)
    ::Fatal("AliNuclexEventCutsSetupRun2pp","You cannot use the legacy centrality framework in pp. Please set the fCentralityFramework to 0 to disable the multiplicity selection or to 1 to use AliMultSelection.");
  fCentEstimators[0] = "V0M";
  fCentEstimators[1] = "CL0";
  fMinCentrality = 0.f;
  fMaxCentrality = 100.f;

  fTriggerMask = AliVEvent::kINT7;

}

void AliNuclexEventCuts::SetupLHC15o() {
  ::Info("AliNuclexEventCuts::SetupLHC15o","Setup event cuts for the LHC15o period.");
  SetName("StandardLHC15oEventCuts");

  fRequireTrackVertex = true;
  fMinVtz = -10.f;
  fMaxVtz = 10.f;
  fMaxDeltaSpdTrackAbsolute = 0.2f;
  fMaxDeltaSpdTrackNsigmaSPD = 10.f;
  fMaxDeltaSpdTrackNsigmaTrack = 20.f;
  fMaxResolutionSPDvertex = 0.25f;

  fRejectDAQincomplete = true;

  fRequiredSolenoidPolarity = 0;

  fSPDpileupMinContributors = 5;
  fSPDpileupMinZdist = 0.8;
  fSPDpileupNsigmaZdist = 3.;
  fSPDpileupNsigmaDiamXY = 2.;
  fSPDpileupNsigmaDiamZ = 5.;
  fTrackletBGcut = false;

  fCentralityFramework = 1;
  fCentEstimators[0] = "V0M";
  fCentEstimators[1] = "CL0";
  fMinCentrality = 0.f;
  fMaxCentrality = 90.f;

  fEstimatorsCorrelationCoef[0] = 0.0157497;
  fEstimatorsCorrelationCoef[1] = 0.973488;
  fEstimatorsSigmaPars[0] = 0.673612;
  fEstimatorsSigmaPars[1] = 0.0290718;
  fEstimatorsSigmaPars[2] = -0.000546728;
  fEstimatorsSigmaPars[3] = 5.82749e-06;
  fDeltaEstimatorNsigma[0] = 5.;
  fDeltaEstimatorNsigma[1] = 5.5;

  fTriggerMask = AliVEvent::kINT7;

}


