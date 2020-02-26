#include "AliEventCuts.h"

#include <algorithm>
#include <array>
using std::array;
#include <memory>
using std::string;
using std::vector;

#include <TClonesArray.h>
#include <TH1D.h>
#include <TH1I.h>
#include <TH2D.h>
#include <TH2F.h>
#include <TObjArray.h>
#include <TObjString.h>

#include <AliAnalysisManager.h>
#include <AliAODMCParticle.h>
#include <AliAODMCHeader.h>
#include <AliCentrality.h>
#include <AliESDtrackCuts.h>
#include <AliInputEventHandler.h>
#include <AliMCEvent.h>
#include <AliMCEventHandler.h>
#include <AliMultSelection.h>
#include <AliMultSelectionTask.h>
#include <AliVEventHandler.h>
#include <AliVMultiplicity.h>
#include <AliVVZERO.h>

#include <AliAODEvent.h>
#include <AliESDEvent.h>

#include <AliLog.h>

ClassImp(AliEventCutsContainer);
ClassImp(AliEventCuts);



/// Standard constructor with null selection
///
AliEventCuts::AliEventCuts(bool saveplots) : TList(),
  fUtils{},
  fGreenLight{false},
  fMC{false},
  fRequireTrackVertex{false},
  fMinVtz{-1000.f},
  fMaxVtz{1000.f},
  fMaxDeltaSpdTrackAbsolute{1.e14},
  fMaxDeltaSpdTrackNsigmaSPD{1.e14},
  fMaxDeltaSpdTrackNsigmaTrack{1.e14},
  fMaxResolutionSPDvertex{1.e14f},
  fMaxDispersionSPDvertex{1.e14f},
  fCheckAODvertex{true},
  fRejectDAQincomplete{false},
  fRequiredSolenoidPolarity{0},
  fUseCombinedMVSPDcut{false},
  fUseMultiplicityDependentPileUpCuts{false},
  fUseSPDpileUpCut{false},
  fSPDpileupMinContributors{1000},
  fSPDpileupMinZdist{-1.f},
  fSPDpileupNsigmaZdist{-1.f},
  fSPDpileupNsigmaDiamXY{-1.f},
  fSPDpileupNsigmaDiamZ{-1.f},
  fTrackletBGcut{false},
  fPileUpCutMV{false},
  fCentralityFramework{0},
  fMinCentrality{-1000.f},
  fMaxCentrality{1000.f},
  fUseVariablesCorrelationCuts{false},
  fUseEstimatorsCorrelationCut{false},
  fUseStrongVarCorrelationCut{false},
  fEstimatorsCorrelationCoef{0.,1.},
  fEstimatorsSigmaPars{10000.,0.,0.,0.},
  fDeltaEstimatorNsigma{1.,1.},
  fTOFvsFB32correlationPars{0.},
  fTOFvsFB32sigmaPars{1.e8,0.,0.,0.,0.,0.},
  fTOFvsFB32nSigmaCut{1.,1.},
  fESDvsTPConlyLinearCut{1.e8,0.},
  fMultiplicityV0McorrCut{nullptr},
  fFB128vsTrklLinearCut{1.e8,0.},
  fVZEROvsTPCoutPolCut{1.e8,0.,0.,0.,0.},
  fRequireExactTriggerMask{false},
  fTriggerMask{AliVEvent::kAny},
  fTriggerClasses{},
  fContainer{},
  fkLabels{"raw","selected"},
  fManualMode{false},
  fSavePlots{saveplots},
  fCurrentRun{-0xBADCAFE},
  fFlag{BIT(kNoCuts)},
  fCentEstimators{"V0M","CL0"},
  fCentPercentiles{-1.f},
  fPrimaryVertex{nullptr},
  fNewEvent{true},
  fOverrideAutoTriggerMask{false},
  fOverrideAutoPileUpCuts{false},
  fMultSelectionEvCuts{false},  
  fUseTimeRangeCut{false},
  fSelectInelGt0{false},
  fOverrideInelGt0{false},
  fOverrideCentralityFramework{false},
  fTimeRangeCut{},
  fCutStats{nullptr},
  fCutStatsAfterTrigger{nullptr},
  fCutStatsAfterMultSelection{nullptr},
  fNormalisationHist{nullptr},
  fVtz{nullptr},
  fDeltaTrackSPDvtz{nullptr},
  fCentrality{nullptr},
  fEstimCorrelation{nullptr},
  fMultCentCorrelation{nullptr},
  fTOFvsFB32{nullptr},
  fTPCvsAll{nullptr},
  fMultvsV0M{nullptr},
  fTPCvsTrkl{nullptr},
  fVZEROvsTPCout{nullptr},
  fFB32trackCuts{nullptr},
  fTPConlyCuts{nullptr}
{
  SetName("AliEventCuts");
  SetOwner(true);
}

AliEventCuts::~AliEventCuts() { 
  delete fMultiplicityV0McorrCut; 
  delete fFB32trackCuts;
  delete fFB32trackCuts;
}

bool AliEventCuts::AcceptEvent(AliVEvent *ev) {
  if (fGreenLight) return true; /// Bypass all the selections

  /// If not specified the cuts are set according to the run period
  const int current_run = ev->GetRunNumber();
  if (current_run != fCurrentRun) {
    fCurrentRun = current_run;
    if (!fManualMode) {
      ::Info("AliEventCuts::AcceptEvent","Current run (%i) is different from the previous (%i): setting automatically the corresponding event cuts.",current_run,fCurrentRun);
      AutomaticSetup(ev);
    }
    if (fUseTimeRangeCut) {
      fTimeRangeCut.InitFromRunNumber(fCurrentRun);
    }
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
  unsigned int selected_trigger = handl->IsEventSelected() & fTriggerMask;
  if ((selected_trigger == fTriggerMask && fRequireExactTriggerMask) || (selected_trigger && !fRequireExactTriggerMask))
    fFlag |= BIT(kTrigger);

  /// Use of trigger classes overrides the trigger mask
  /// (i.e. if trigger mask is not fired but we see the trigger class we want we enable the trigger bit)
  /// A special bit is set in this case
  TString classes = ev->GetFiredTriggerClasses();
  if (fTriggerClasses.empty())
    fFlag |= BIT(kTriggerClasses);
  for (const std::string& myClass : fTriggerClasses) {
    if (classes.Contains(myClass.data()) && !myClass.empty()) {
      fFlag |= BIT(kTrigger);
      fFlag |= BIT(kTriggerClasses);
      break;
    }
  }

  /// Vertex existance
  const AliVVertex* vtTrc = ev->GetPrimaryVertex();
  bool isTrackV = true;
  if(vtTrc->IsFromVertexer3D() || vtTrc->IsFromVertexerZ()) isTrackV=false;
  const AliVVertex* vtSPD = ev->GetPrimaryVertexSPD();
  /// On current AODs primary vertex could be from TPC or invalid SPD vertex
  /// The following check should be applied only on AOD.
  bool goodAODvtx = (dynamic_cast<AliAODEvent*>(ev) ? GoodPrimaryAODVertex(ev) : true) || !fCheckAODvertex;

  if (vtSPD->GetNContributors() > 0) fFlag |= BIT(kVertexSPD);
  if (vtTrc->GetNContributors() > 1 && isTrackV && goodAODvtx) fFlag |= BIT(kVertexTracks);
  if (((fFlag & BIT(kVertexTracks)) ||  !fRequireTrackVertex) && (fFlag & BIT(kVertexSPD))) fFlag |= BIT(kVertex);
  const AliVVertex* &vtx = bool(fFlag & BIT(kVertexTracks)) ? vtTrc : vtSPD;
  fPrimaryVertex = const_cast<AliVVertex*>(vtx);

  /// Vertex position cut
  if (vtSPD->GetZ() >= fMinVtz && vtSPD->GetZ() <= fMaxVtz) fFlag |= BIT(kVertexPositionSPD);
  if (vtTrc->GetZ() >= fMinVtz && vtTrc->GetZ() <= fMaxVtz) fFlag |= BIT(kVertexPositionTracks);
  if (vtx->GetZ()   >= fMinVtz && vtx->GetZ()   <= fMaxVtz) fFlag |= BIT(kVertexPosition);

  /// Vertex quality cuts
  double covTrc[6],covSPD[6];
  vtTrc->GetCovarianceMatrix(covTrc);
  vtSPD->GetCovarianceMatrix(covSPD);
  double dz = bool(fFlag & kVertexSPD) && bool(fFlag & kVertexTracks) ? vtTrc->GetZ() - vtSPD->GetZ() : 0.; /// If one of the two vertices is not available this cut is always passed.
  double errTot = TMath::Sqrt(covTrc[5]+covSPD[5]);
  double errTrc = bool(fFlag & kVertexTracks) ? TMath::Sqrt(covTrc[5]) : 1.;
  double nsigTot = TMath::Abs(dz) / errTot, nsigTrc = TMath::Abs(dz) / errTrc;
  /// vertex dispersion for run1, only for ESD, AOD code to be added here
  const AliESDVertex* vtSPDESD = dynamic_cast<const AliESDVertex*>(vtSPD);
  double vtSPDdispersion = vtSPDESD ? vtSPDESD->GetDispersion() : 0;
  if (
      (TMath::Abs(dz) <= fMaxDeltaSpdTrackAbsolute && nsigTot <= fMaxDeltaSpdTrackNsigmaSPD && nsigTrc <= fMaxDeltaSpdTrackNsigmaTrack) && // discrepancy track-SPD vertex
      (!vtSPD->IsFromVertexerZ() || TMath::Sqrt(covSPD[5]) <= fMaxResolutionSPDvertex) &&
      (!vtSPD->IsFromVertexerZ() || vtSPDdispersion <= fMaxDispersionSPDvertex) /// vertex dispersion cut for run1, only for ESD
     ) // quality cut on vertexer SPD z
    fFlag |= BIT(kVertexQuality);  

  /// Pile-up rejection
  bool usePileUpMV = (fUseCombinedMVSPDcut && vtx != vtSPD) || fPileUpCutMV;
  bool usePileUpSPD = (fUseCombinedMVSPDcut && vtx == vtSPD) || fUseSPDpileUpCut;
  AliVMultiplicity* mult = ev->GetMultiplicity();
  const int ntrkl = mult->GetNumberOfTracklets();
  if (fUseMultiplicityDependentPileUpCuts) {
    if (ntrkl < 20) fSPDpileupMinContributors = 3;
    else if (ntrkl < 50) fSPDpileupMinContributors = 4;
    else fSPDpileupMinContributors = 5;
  }
  if ((!usePileUpSPD || !ev->IsPileupFromSPD(fSPDpileupMinContributors,fSPDpileupMinZdist,fSPDpileupNsigmaZdist,fSPDpileupNsigmaDiamXY,fSPDpileupNsigmaDiamZ)) &&
      (!fTrackletBGcut || !fUtils.IsSPDClusterVsTrackletBG(ev)) &&
      (!usePileUpMV || !fUtils.IsPileUpMV(ev)))
    fFlag |= BIT(kPileUp);

  /// Centrality cuts:
  /// * Check for min and max centrality
  /// * Cross check correlation between two centrality estimators
  if (fCentralityFramework) {
    if (fCentralityFramework == 2) {
      AliCentrality* cent = ev->GetCentrality();
      if (!cent) {
        AliFatal("The legacy centrality framework has been request but no AliCentrality object was found attached to the Event."
                 " Did you run the Centrality Framework?");
      }
      fCentPercentiles[0] = cent->GetCentralityPercentile(fCentEstimators[0].data());
      fCentPercentiles[1] = cent->GetCentralityPercentile(fCentEstimators[1].data());
    } else {
      AliMultSelection* cent = (AliMultSelection*)ev->FindListObject("MultSelection");
      if (!cent) {
        AliFatal("The multiplicity selection framework has been request but no AliMultSelection object was found attached to the Event."
                 " Did you run the AliMultSelectionTask?");
      }
      fCentPercentiles[0] = cent->GetMultiplicityPercentile(fCentEstimators[0].data(), fMultSelectionEvCuts);
      fCentPercentiles[1] = cent->GetMultiplicityPercentile(fCentEstimators[1].data(), fMultSelectionEvCuts);
    }
    const auto& x = fCentPercentiles[1];
    const double center = x * fEstimatorsCorrelationCoef[1] + fEstimatorsCorrelationCoef[0];
    const double sigma = fEstimatorsSigmaPars[0] + fEstimatorsSigmaPars[1] * x + fEstimatorsSigmaPars[2] * x * x + fEstimatorsSigmaPars[3] * x * x * x;
    if ((!fUseEstimatorsCorrelationCut || fMC ||
          (fCentPercentiles[0] >= center - fDeltaEstimatorNsigma[0] * sigma && fCentPercentiles[0] <= center + fDeltaEstimatorNsigma[1] * sigma))
        && fCentPercentiles[0] >= fMinCentrality
        && fCentPercentiles[0] <= fMaxCentrality) {
          fFlag |= BIT(kMultiplicity);
    }
  } else
    fFlag |= BIT(kMultiplicity);

  if (AliMultSelectionTask::IsINELgtZERO(ev) || !fSelectInelGt0) {
    fFlag |= BIT(kINELgt0);
  } else if (!fSelectInelGt0) { // for simple INEL events multiplicity is not defined.
    fCentPercentiles[0] = -0.5;
    fCentPercentiles[1] = -0.5;
  }

  /// If the correlation plots are defined, we should fill them
  if (fUseVariablesCorrelationCuts || fTOFvsFB32[0]) {
    ComputeTrackMultiplicity(ev);
    const double fb32 = fContainer.fMultTrkFB32;
    const double fb32acc = fContainer.fMultTrkFB32Acc;
    const double fb32tof = fContainer.fMultTrkFB32TOF;
    const double fb128 = fContainer.fMultTrkTPC;
    const double esd = fContainer.fMultESD;

    const double mu32tof = PolN(fb32,fTOFvsFB32correlationPars,3);
    const double sigma32tof = PolN(fb32,fTOFvsFB32sigmaPars, 5);
    const double vzero_tpcout_limit = PolN(double(fContainer.fMultTrkTPCout),fVZEROvsTPCoutPolCut,4);

    const bool multV0Mcut = (fMultiplicityV0McorrCut) ? fb32acc > fMultiplicityV0McorrCut->Eval(fCentPercentiles[0]) : true;

    if (((fb32tof <= mu32tof + fTOFvsFB32nSigmaCut[0] * sigma32tof && fb32tof >= mu32tof - fTOFvsFB32nSigmaCut[1] * sigma32tof) &&
        (esd < fESDvsTPConlyLinearCut[0] + fESDvsTPConlyLinearCut[1] * fb128) &&
        multV0Mcut &&
        (fb128 < fFB128vsTrklLinearCut[0] + fFB128vsTrklLinearCut[1] * ntrkl) &&
        (!fUseStrongVarCorrelationCut || fContainer.fMultVZERO > vzero_tpcout_limit))
        || fMC || !fUseVariablesCorrelationCuts)
      fFlag |= BIT(kCorrelations);
  } else fFlag |= BIT(kCorrelations);

  /// Time Range masking
  if (fUseTimeRangeCut) {
    if ( fTimeRangeCut.CutEvent(ev) == kFALSE ) {
      // good event: should be accepted
      fFlag |= BIT(kTimeRangeCut);
    }
  } else {
    fFlag |= BIT(kTimeRangeCut);
  }

  /// Ignore SPD/tracks vertex position and reconstruction individual flags
  bool allcuts = CheckNormalisationMask(kPassesAllCuts);
  if (allcuts) {
    fFlag |= BIT(kAllCuts);
  }
  if (fCutStats) {
    for (int iCut = kNoCuts; iCut <= kAllCuts; ++iCut) {
      if (TESTBIT(fFlag,iCut)) {
        fCutStats->Fill(iCut);
        if (TESTBIT(fFlag,kTrigger)) {
          fCutStatsAfterTrigger->Fill(iCut);
        }
        if (TESTBIT(fFlag,kMultiplicity)) {
          fCutStatsAfterMultSelection->Fill(iCut);
        }
      }
    }
  }

  /// Filling normalisation histogram
  array <NormMask,5> norm_masks {
    kAnyEvent,
    kTriggeredEvent,
    kPassesNonVertexRelatedSelections,
    kHasReconstructedVertex,
    kPassesAllCuts
  };
  for (int iC = 0; iC < 5; ++iC) {
    if (CheckNormalisationMask(norm_masks[iC])) {
      if (fNormalisationHist) {
        fNormalisationHist->Fill(iC);
      }
    }
  }

  /// Filling the monitoring histograms (first iteration always filled, second iteration only for selected events.
  for (int befaft = 0; befaft < 2; ++befaft) {
    if (fCentrality[befaft]) fCentrality[befaft]->Fill(fCentPercentiles[0]);
    if (fEstimCorrelation[befaft]) fEstimCorrelation[befaft]->Fill(fCentPercentiles[1],fCentPercentiles[0]);
    if (fMultCentCorrelation[befaft]) fMultCentCorrelation[befaft]->Fill(fCentPercentiles[0],ntrkl);
    if (fVtz[befaft]) fVtz[befaft]->Fill(vtx->GetZ());
    if (fDeltaTrackSPDvtz[befaft]) fDeltaTrackSPDvtz[befaft]->Fill(dz);
    if (fTOFvsFB32[befaft]) fTOFvsFB32[befaft]->Fill(fContainer.fMultTrkFB32,fContainer.fMultTrkFB32TOF);
    if (fTPCvsAll[befaft])  fTPCvsAll[befaft]->Fill(fContainer.fMultTrkTPC,float(fContainer.fMultESD) - fESDvsTPConlyLinearCut[1] * fContainer.fMultTrkTPC);
    if (fMultvsV0M[befaft]) fMultvsV0M[befaft]->Fill(GetCentrality(),fContainer.fMultTrkFB32Acc);
    if (fTPCvsTrkl[befaft]) fTPCvsTrkl[befaft]->Fill(ntrkl,fContainer.fMultTrkTPC);
    if (fVZEROvsTPCout[befaft]) fVZEROvsTPCout[befaft]->Fill(fContainer.fMultTrkTPCout,fContainer.fMultVZERO);
    if (!allcuts) return false; /// Do not fill the "after" histograms if the event does not pass the cuts.
  }

  return true;
}

void AliEventCuts::AddQAplotsToList(TList *qaList, bool addCorrelationPlots) {
  if (!qaList) {
    if (fSavePlots)
      qaList = static_cast<TList*>(this);
    else
      ::Fatal("AliEventCuts::AddQAplotsToList","Invalid TList pointer provided");
  }
  static_cast<TList*>(this)->Delete();

  vector<string> bin_labels = {
    "No cuts",
    "DAQ Incomplete",
    "Magnetic field choice",
    "Trigger selection",
    "Trigger classes",
    "SPD vertex reconstruction",
    "Track vertex reconstruction",
    "Vertex reconstruction",
    "SPD vertex position",
    "Track vertex position",
    "Vertex position",
    "Vertex quality",
    "Pile-up",
    "Multiplicity selection",
    "INEL > 0",
    "Correlations",
    "TimeRangeCut",
    "All cuts"
  };

  vector<string> norm_labels = {
    "No cuts",
    "Trigger selection",
    "Event selection",
    "Vertex reconstruction and quality",
    "Vertex position"
  };

  fCutStats = new TH1D("fCutStats",";;Number of selected events",bin_labels.size(),-.5,bin_labels.size() - 0.5);
  fCutStatsAfterTrigger = new TH1D("fCutStatsAfterTrigger",";;Number of selected events",bin_labels.size(),-.5,bin_labels.size() - 0.5);
  fCutStatsAfterMultSelection = new TH1D("fCutStatsAfterMultSelection",";;Number of selected events",bin_labels.size(),-.5,bin_labels.size() - 0.5);
  fNormalisationHist = new TH1D("fNormalisationHist",";;Number of selected events",norm_labels.size(),-.5,norm_labels.size() - 0.5);
  for (size_t iB = 1; iB <= bin_labels.size(); ++iB) {
    fCutStats->GetXaxis()->SetBinLabel(iB,bin_labels[iB-1].data());
    fCutStatsAfterTrigger->GetXaxis()->SetBinLabel(iB,bin_labels[iB-1].data());
    fCutStatsAfterMultSelection->GetXaxis()->SetBinLabel(iB,bin_labels[iB-1].data());
  }
  for (size_t iB = 1; iB <= norm_labels.size(); ++iB) fNormalisationHist->GetXaxis()->SetBinLabel(iB,norm_labels[iB-1].data());
  qaList->Add(fCutStats);
  qaList->Add(fCutStatsAfterTrigger);
  qaList->Add(fCutStatsAfterMultSelection);
  qaList->Add(fNormalisationHist);

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

    if (addCorrelationPlots) {
      fTOFvsFB32[iS] = new TH2F(Form("fTOFvsFB32_%s",fkLabels[iS].data()),Form("%s;Multiplicity FB32;Multiplicity FB32+TOF",titles[iS].data()),500,0.,3000.,500,0.,2000.);
      fTPCvsAll[iS] = new TH2F(Form("fTPCvsAll_%s",fkLabels[iS].data()),Form("%s;Multiplicity TPC;Multiplicity ESD - a_{0} #times Multiplicity TPC",titles[iS].data()),300,0.,3000.,3000,-200.,28800.);
      fMultvsV0M[iS] = new TH2F(Form("fMultvsV0M_%s",fkLabels[iS].data()),Form("%s;Centrality (V0M);FB32",titles[iS].data()),100,0.,100.,500,0.,3000.);
      fTPCvsTrkl[iS] = new TH2F(Form("fTPCvsTrkl_%s",fkLabels[iS].data()),Form("%s;SPD tracklets;FB128",titles[iS].data()),2500,0.,5000.,2500,0.,5000.);
      fVZEROvsTPCout[iS] = new TH2F(Form("fVZEROvsTPCout_%s",fkLabels[iS].data()),Form("%s;Tracks with kTPCout on;VZERO multiplicity",titles[iS].data()),2000,0.,20000.,2000,0.,40000.);
      qaList->Add(fTOFvsFB32[iS]);
      qaList->Add(fTPCvsAll[iS]);
      qaList->Add(fMultvsV0M[iS]);
      qaList->Add(fTPCvsTrkl[iS]);
      qaList->Add(fVZEROvsTPCout[iS]);
    }
  }

}

void AliEventCuts::AutomaticSetup(AliVEvent *ev) {
  if (dynamic_cast<AliAODEvent*>(ev)) {
    TClonesArray *stack = (TClonesArray*)ev->GetList()->FindObject(AliAODMCParticle::StdBranchName());
    fMC = (stack) ? true : false;
  } else if (!dynamic_cast<AliESDEvent*>(ev))
    AliFatal("I don't find the AOD event nor the ESD one, aborting.");
  else {
    AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
    TClonesArray* aodMC = (TClonesArray*)ev->GetList()->FindObject(AliAODMCParticle::StdBranchName());
    fMC = (eventHandler || aodMC);
  }
  
  if (fCurrentRun == -1 && fMC) {
    ::Info("AliEventCuts::AutomaticSetup","MCGEN train / Kinematics only production detected, disabling all the cuts.");
    fGreenLight = true;
    return;
  }

  if (fCurrentRun == 280234 || fCurrentRun == 280235) {
    SetupLHC17n();
    return;
  }

  if ((fCurrentRun >= 166423 && fCurrentRun <= 170593) ||
      (fCurrentRun >= 136833 && fCurrentRun <= 139517)) {
    SetupRun1PbPb();
    return;
  }

  /// Run 2 Pb-Pb: LHC15o
  if (fCurrentRun >= 244917 && fCurrentRun <= 246994) {
    SetupRun2PbPb();
    return;
  }

  if (fCurrentRun >= 295369 && fCurrentRun <= 297624) {
    SetupPbPb2018();
    return;
  }

  /// Run 2 p-Pb
  if ((fCurrentRun >= 264896 && fCurrentRun <= 265533) || /// LHC16q: p-Pb 5 TeV
      (fCurrentRun >= 265534 && fCurrentRun <= 266329) || /// LHC16r: p-Pb 8 TeV
      (fCurrentRun >= 267139 && fCurrentRun <= 267166)) { /// LHC16t: p-Pb 5 TeV
    SetupRun2pA(0);
    return;
  }
  if (fCurrentRun >= 266330 && fCurrentRun <= 267138) {   /// LHC16s: Pb-p 5 TeV
    SetupRun2pA(1);
    return;
  }

  /// Run 1 p-Pb
  if (fCurrentRun >= 195344 && fCurrentRun <= 196311) {  /// LHC13bcde: p-Pb 5 TeV
    SetupRun1pA(0);
    return;
  }
  if (fCurrentRun >= 196433 && fCurrentRun <= 197388) {   /// LHC13f: Pb-p 5 TeV
    SetupRun1pA(1);
    return;
  }

  if ((fCurrentRun >= 225000 && fCurrentRun <= 244628) || // 2015 5+13 TeV sample
      (fCurrentRun >= 252235 && fCurrentRun <= 264347) || // 2016 13 TeV sample
      (fCurrentRun >= 270531 && fCurrentRun <= 294960)) { // 2017 13 and 5 TeV + 2018 13 TeV samples
    SetupRun2pp();
    return;
  }
  
  /// Run 1 pp
  if ((fCurrentRun >= 140164 && fCurrentRun <= 146860) ||  // pp 2.76 TeV (LHC11a)
      (fCurrentRun >= 121692 && fCurrentRun <= 126437) ||  // pp 7 TeV (LHC10d)
      (fCurrentRun >= 172439 && fCurrentRun <= 193766)) {  // pp 8 TeV (LHC12a-i)
    SetupRun1pp();
    return;
  }

  ::Fatal("AliEventCuts::AutomaticSetup","Automatic period detection failed, please use the manual mode. If you need this period in AliEventCuts send an email to the DPG event selection mailing list.");
}

float AliEventCuts::GetCentrality (unsigned int estimator) const {
  if (estimator > 1) {
    /// Escalate to Fatal
    ::Fatal("AliEventCuts::GetCentrality","You asked for the centrality estimator with index %i, but you should choose between index 0 and 1.", estimator);
    return -1.f;
  } else
    return fCentPercentiles[estimator];
}

string AliEventCuts::GetCentralityEstimator (unsigned int estimator) const {
  if (estimator > 1) {
    /// Escalate to Fatal
    ::Fatal("AliEventCuts::GetCentralityEstimator","You asked for the centrality estimator with index %i, but you should choose between index 0 and 1.", estimator);
    return "";
  } else
    return fCentEstimators[estimator];
}


void AliEventCuts::ComputeTrackMultiplicity(AliVEvent *ev) {
  AliEventCutsContainer* tmp_cont = static_cast<AliEventCutsContainer*>(ev->FindListObject("AliEventCutsContainer"));
  if (tmp_cont) {
    unsigned long evid = ((unsigned long)(ev->GetBunchCrossNumber()) << 32) + ev->GetTimeStamp();
    fNewEvent = (tmp_cont->fEventId != evid);
    tmp_cont->fEventId = evid;

    if (!fNewEvent) {
      fContainer = *tmp_cont;
      return;
    }
  } else {
    tmp_cont = new AliEventCutsContainer;
    ev->AddObject(tmp_cont);
  }

  bool isAOD = false;
  if (dynamic_cast<AliAODEvent*>(ev))
    isAOD = true;
  else if (!dynamic_cast<AliESDEvent*>(ev))
    AliFatal("I don't find the AOD event nor the ESD one, aborting.");

  if (!fFB32trackCuts) fFB32trackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011();
  if (!fTPConlyCuts) fTPConlyCuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();

  const int nTracks = ev->GetNumberOfTracks();
  tmp_cont->fMultESD = (isAOD) ? ((AliAODHeader*)ev->GetHeader())->GetNumberOfESDTracks() : dynamic_cast<AliESDEvent*>(ev)->GetNumberOfTracks();
  tmp_cont->fMultTrkFB32 = 0;
  tmp_cont->fMultTrkFB32Acc = 0;
  tmp_cont->fMultTrkFB32TOF = 0;
  tmp_cont->fMultTrkTPC = 0;
  tmp_cont->fMultTrkTPCout = 0;
  for (int it = 0; it < nTracks; it++) {
    if (isAOD) {
      AliAODTrack* trk = (AliAODTrack*)ev->GetTrack(it);
      if (!trk) continue;
      if ((trk->GetStatus() & AliESDtrack::kTPCout) &&
          trk->GetID() > 0) tmp_cont->fMultTrkTPCout++;
      if (trk->TestFilterBit(32)) {
        tmp_cont->fMultTrkFB32++;
        if ( TMath::Abs(trk->GetTOFsignalDz()) <= 10. && trk->GetTOFsignal() >= 12000. && trk->GetTOFsignal() <= 25000.)
          tmp_cont->fMultTrkFB32TOF++;
        if ((fabs(trk->Eta()) < 0.8) && (trk->GetTPCNcls() >= 70) && (trk->Pt() >= 0.2) && (trk->Pt() < 50))
          tmp_cont->fMultTrkFB32Acc++;
      }
      if (trk->TestFilterBit(128))
        tmp_cont->fMultTrkTPC++;
    } else {
      AliESDtrack* esdTrack = (AliESDtrack*)ev->GetTrack(it);
      if (!esdTrack) continue;

      if (esdTrack->GetStatus() & AliESDtrack::kTPCout) tmp_cont->fMultTrkTPCout++;

      if (fFB32trackCuts->AcceptTrack(esdTrack)) {
        tmp_cont->fMultTrkFB32++;
        if (TMath::Abs(esdTrack->GetTOFsignalDz()) <= 10 && esdTrack->GetTOFsignal() >= 12000 && esdTrack->GetTOFsignal() <= 25000)
          tmp_cont->fMultTrkFB32TOF++;

        if ((TMath::Abs(esdTrack->Eta()) < 0.8) && (esdTrack->GetTPCNcls() > 70) && (esdTrack->Pt() > 0.2) && (esdTrack->Pt() < 50))
          tmp_cont->fMultTrkFB32Acc++;
      }

      /// TPC only tracks, with the same cuts of the filter bit 128
      AliESDtrack tpcParam;
      if (!esdTrack->FillTPCOnlyTrack(tpcParam)) continue;
      if (!fTPConlyCuts->AcceptTrack(&tpcParam)) continue;
      if (tpcParam.Pt() > 0.) {
        // only constrain tracks above threshold
        AliExternalTrackParam exParam;
        // take the B-field from the ESD, no 3D fieldMap available at this point
        bool relate = false;
        relate = tpcParam.RelateToVertexTPC((AliESDVertex*)ev->GetPrimaryVertexSPD(),ev->GetMagneticField(),kVeryBig, &exParam);
        if(!relate) continue;
      }
      tmp_cont->fMultTrkTPC++;
    }
  }
  AliVVZERO *vzero = (AliVVZERO*)ev->GetVZEROData();
  if(vzero) {
    tmp_cont->fMultVZERO = 0.;
    for(int ich=0; ich < 64; ich++)
      tmp_cont->fMultVZERO += vzero->GetMultiplicity(ich);
  }
  fContainer = *tmp_cont;
}

void AliEventCuts::SetupRun1pp() {
  ::Info("AliEventCuts::SetupRun1pp","EXPERIMENTAL: Setup event cuts for the Run1 pp periods. Currently only for ESD");
  SetName("StandardRun1ppEventCuts");

  fMinVtz = -10.f;
  fMaxVtz = 10.f;
  fMaxDeltaSpdTrackAbsolute = 0.5f;
  fMaxDeltaSpdTrackNsigmaSPD = 1.e14f;
  fMaxDeltaSpdTrackNsigmaTrack = 1.e14;
  fMaxResolutionSPDvertex = 0.25f;
  fMaxDispersionSPDvertex = 0.04f;

  fRejectDAQincomplete = true;
  
  fTrackletBGcut = true;

  fRequiredSolenoidPolarity = 0;
  /*
  //remove pileup cuts, to be confirmed
  if (!fOverrideAutoPileUpCuts) {
    fUseCombinedMVSPDcut = true;
    fSPDpileupMinZdist = 0.8;
    fSPDpileupNsigmaZdist = 3.;
    fSPDpileupNsigmaDiamXY = 2.;
    fSPDpileupNsigmaDiamZ = 5.;
    fTrackletBGcut = true;
    fUseMultiplicityDependentPileUpCuts = true;
    fUtils.SetMaxPlpChi2MV(5);
    fUtils.SetMinWDistMV(15);
    fUtils.SetCheckPlpFromDifferentBCMV(kFALSE);
  }
  */

  if (!fOverrideAutoTriggerMask) fTriggerMask = AliVEvent::kINT7 |  AliVEvent::kCINT5 |  AliVEvent::kMB;

}

void AliEventCuts::SetupRun2pp() {
  ::Info("AliEventCuts::SetupRun2pp","Setup event cuts for the Run2 pp periods.");
  SetName("StandardRun2ppEventCuts");

  fMinVtz = -10.f;
  fMaxVtz = 10.f;
  fMaxDeltaSpdTrackAbsolute = 0.5f;
  fMaxDeltaSpdTrackNsigmaSPD = 1.e14f;
  fMaxDeltaSpdTrackNsigmaTrack = 1.e14;
  fMaxResolutionSPDvertex = 0.25f;

  fRejectDAQincomplete = true;

  fRequiredSolenoidPolarity = 0;

  if (!fOverrideAutoPileUpCuts) {
    fUseCombinedMVSPDcut = true;
    fSPDpileupMinZdist = 0.8;
    fSPDpileupNsigmaZdist = 3.;
    fSPDpileupNsigmaDiamXY = 2.;
    fSPDpileupNsigmaDiamZ = 5.;
    fTrackletBGcut = true;
    fUseMultiplicityDependentPileUpCuts = true;
    fUtils.SetMaxPlpChi2MV(5);
    fUtils.SetMinWDistMV(15);
    fUtils.SetCheckPlpFromDifferentBCMV(kFALSE);
  }

  if (fCentralityFramework > 1)
    ::Fatal("AliEventCutsSetupRun2pp","You cannot use the legacy centrality framework in pp. Please set the fCentralityFramework to 0 to disable the multiplicity selection or to 1 to use AliMultSelection.");
  else if (fCentralityFramework == 1) {
    fCentEstimators[0] = "V0M";
    fCentEstimators[1] = "CL0";
    fSelectInelGt0 = fOverrideInelGt0 ? fSelectInelGt0 : true;
  }

  fFB128vsTrklLinearCut[0] = 32.077;
  fFB128vsTrklLinearCut[1] = 0.932;

  if (!fOverrideAutoTriggerMask) fTriggerMask = AliVEvent::kINT7;

}

void AliEventCuts::SetupPbPb2018() {
  ::Info("AliEventCuts::SetupPbPb2018","EXPERIMENTAL: Setup event cuts for the 2018 Pb-Pb periods.");
  SetName("StandardSetupPbPb2018");

  fRequireTrackVertex = true;
  fMinVtz = -10.f;
  fMaxVtz = 10.f;
  fMaxDeltaSpdTrackAbsolute = 0.2f;
  fMaxDeltaSpdTrackNsigmaSPD = 10.f;
  fMaxDeltaSpdTrackNsigmaTrack = 20.f;
  fMaxResolutionSPDvertex = 0.25f;

  fRejectDAQincomplete = true;

  fRequiredSolenoidPolarity = 0;

  if (!fOverrideAutoPileUpCuts) {
    fUseMultiplicityDependentPileUpCuts = true; // If user specify a value it is not overwritten
    fSPDpileupMinZdist = 0.8;
    fSPDpileupNsigmaZdist = 3.;
    fSPDpileupNsigmaDiamXY = 2.;
    fSPDpileupNsigmaDiamZ = 5.;
    fTrackletBGcut = false;
  }

  if (!fOverrideCentralityFramework) {
    fCentralityFramework = 1;
    fCentEstimators[0] = "V0M";
    fCentEstimators[1] = "CL0";
    fMinCentrality = 0.f;
    fMaxCentrality = 90.f;
  }

  fUseEstimatorsCorrelationCut = false;
  fEstimatorsCorrelationCoef[0] = -0.669108;
  fEstimatorsCorrelationCoef[1] = 1.04489;
  fEstimatorsSigmaPars[0] = 0.933321;
  fEstimatorsSigmaPars[1] = 0.0416976;
  fEstimatorsSigmaPars[2] = -0.000936344;
  fEstimatorsSigmaPars[3] = 8.92179e-06;
  fDeltaEstimatorNsigma[0] = 5.5;
  fDeltaEstimatorNsigma[1] = 5.;

  array<double,4> tof_fb32_corr = {-1.0178, 0.333132, 9.10282e-05, -1.61861e-08};
  std::copy(tof_fb32_corr.begin(),tof_fb32_corr.end(),fTOFvsFB32correlationPars);
  array<double,6> tof_fb32_sigma = {1.47848, 0.0385923, -5.06153e-05, 4.37641e-08, -1.69082e-11, 2.35085e-15};
  std::copy(tof_fb32_sigma.begin(),tof_fb32_sigma.end(),fTOFvsFB32sigmaPars);
  fTOFvsFB32nSigmaCut[0] = 4.;
  fTOFvsFB32nSigmaCut[1] = 4.;

  fESDvsTPConlyLinearCut[0] = 15000.;
  fESDvsTPConlyLinearCut[1] = 3.38;

  array<double,5> vzero_tpcout_polcut = {-2000.,2.1,3.5e-5,0.,0.};
  std::copy(vzero_tpcout_polcut.begin(),vzero_tpcout_polcut.end(),fVZEROvsTPCoutPolCut);

  if (fCentralityFramework != 0) {
    if(!fMultiplicityV0McorrCut) fMultiplicityV0McorrCut = new TF1("fMultiplicityV0McorrCut","[0]+[1]*x+[2]*exp([3]-[4]*x) - 5.*([5]+[6]*exp([7]-[8]*x))",0,100);
    fMultiplicityV0McorrCut->SetParameters(-6.15980e+02, 4.89828e+00, 4.84776e+03, -5.22988e-01, 3.04363e-02, -1.21144e+01, 2.95321e+02, -9.20062e-01, 2.17372e-02);
  }

  if (!fOverrideAutoTriggerMask) {
    fTriggerMask = AliVEvent::kINT7 | AliVEvent::kCentral | AliVEvent::kSemiCentral;
  }

  fUseTimeRangeCut = true;

}

void AliEventCuts::SetupRun2PbPb() {
  ::Info("AliEventCuts::SetupRun2PbPb","Setup event cuts for the Run2 Pb-Pb periods.");
  SetName("StandardRun2PbPbEventCuts");

  fRequireTrackVertex = true;
  fMinVtz = -10.f;
  fMaxVtz = 10.f;
  fMaxDeltaSpdTrackAbsolute = 0.2f;
  fMaxDeltaSpdTrackNsigmaSPD = 10.f;
  fMaxDeltaSpdTrackNsigmaTrack = 20.f;
  fMaxResolutionSPDvertex = 0.25f;

  fRejectDAQincomplete = true;

  fRequiredSolenoidPolarity = 0;

  if (!fOverrideAutoPileUpCuts) {
    fUseMultiplicityDependentPileUpCuts = true; // If user specify a value it is not overwritten
    fSPDpileupMinZdist = 0.8;
    fSPDpileupNsigmaZdist = 3.;
    fSPDpileupNsigmaDiamXY = 2.;
    fSPDpileupNsigmaDiamZ = 5.;
    fTrackletBGcut = false;
  }

  if (!fOverrideCentralityFramework) {
    fCentralityFramework = 1;
    fCentEstimators[0] = "V0M";
    fCentEstimators[1] = "CL0";
    fMinCentrality = 0.f;
    fMaxCentrality = 90.f;
  }

  fUseEstimatorsCorrelationCut = fCentralityFramework != 0;
  fEstimatorsCorrelationCoef[0] = 0.0157497;
  fEstimatorsCorrelationCoef[1] = 0.973488;
  fEstimatorsSigmaPars[0] = 0.673612;
  fEstimatorsSigmaPars[1] = 0.0290718;
  fEstimatorsSigmaPars[2] = -0.000546728;
  fEstimatorsSigmaPars[3] = 5.82749e-06;
  fDeltaEstimatorNsigma[0] = 5.;
  fDeltaEstimatorNsigma[1] = 5.5;

  array<double,4> tof_fb32_corr = {-1.0178, 0.333132, 9.10282e-05, -1.61861e-08};
  std::copy(tof_fb32_corr.begin(),tof_fb32_corr.end(),fTOFvsFB32correlationPars);
  array<double,6> tof_fb32_sigma = {1.47848, 0.0385923, -5.06153e-05, 4.37641e-08, -1.69082e-11, 2.35085e-15};
  std::copy(tof_fb32_sigma.begin(),tof_fb32_sigma.end(),fTOFvsFB32sigmaPars);
  fTOFvsFB32nSigmaCut[0] = 4.;
  fTOFvsFB32nSigmaCut[1] = 4.;

  fESDvsTPConlyLinearCut[0] = 15000.;
  fESDvsTPConlyLinearCut[1] = 3.38;

  array<double,5> vzero_tpcout_polcut = {-2000.,2.1,3.5e-5,0.,0.};
  std::copy(vzero_tpcout_polcut.begin(),vzero_tpcout_polcut.end(),fVZEROvsTPCoutPolCut);

  if (fCentralityFramework != 0) {
    if(!fMultiplicityV0McorrCut) fMultiplicityV0McorrCut = new TF1("fMultiplicityV0McorrCut","[0]+[1]*x+[2]*exp([3]-[4]*x) - 5.*([5]+[6]*exp([7]-[8]*x))",0,100);
    fMultiplicityV0McorrCut->SetParameters(-6.15980e+02, 4.89828e+00, 4.84776e+03, -5.22988e-01, 3.04363e-02, -1.21144e+01, 2.95321e+02, -9.20062e-01, 2.17372e-02);
  }
  if (!fOverrideAutoTriggerMask) {
    fTriggerMask = AliVEvent::kINT7;
    if (fCurrentRun >= 295369)
      fTriggerMask |= AliVEvent::kCentral | AliVEvent::kSemiCentral;
  }

}

void AliEventCuts::SetupRun1PbPb() {
  fRequireTrackVertex = true;
  fMinVtz = -10.f;
  fMaxVtz = 10.f;

  fCentralityFramework = 2;
  fCentEstimators[0] = "V0M";
  fCentEstimators[1] = "CL0";
  fMinCentrality = 0.f;
  fMaxCentrality = 90.f;

  fUseEstimatorsCorrelationCut = false;
  array<double,5> vzero_tpcout_polcut{-1391.46,1.44779,4.78352e-05,-5.96914e-10,0.0};
  std::copy(vzero_tpcout_polcut.begin(),vzero_tpcout_polcut.end(),fVZEROvsTPCoutPolCut);

  if (!fOverrideAutoTriggerMask) fTriggerMask = AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral;
}

void AliEventCuts::SetupRun1pA(int iPeriod) {
  ::Info("AliEventCuts::SetupRun1pA","Event cuts for pA are being setup.");
  SetName("StandardRun1pAEventCuts");

  /// p--Pb requires nsigma cuts on primary vertex
  fMaxDeltaSpdTrackNsigmaSPD = 20.f;
  fMaxDeltaSpdTrackNsigmaTrack = 40.f;

  if (!fOverrideAutoPileUpCuts) {
    fUseCombinedMVSPDcut = true; // If user specify a value it is not overwritten
    fSPDpileupMinZdist = 0.8;
    fSPDpileupNsigmaZdist = 3.;
    fSPDpileupNsigmaDiamXY = 2.;
    fSPDpileupNsigmaDiamZ = 5.;
    fTrackletBGcut = false;
    fUseMultiplicityDependentPileUpCuts = true;
    fUtils.SetMaxPlpChi2MV(5);
    fUtils.SetMinWDistMV(15);
    fUtils.SetCheckPlpFromDifferentBCMV(kFALSE);
  }

  fMinVtz = -10.f;
  fMaxVtz = 10.f;
  fMaxDeltaSpdTrackAbsolute = 0.5f;
  fMaxDeltaSpdTrackNsigmaSPD = 1.e14f;
  fMaxDeltaSpdTrackNsigmaTrack = 1.e14;
  fMaxResolutionSPDvertex = 0.25f;

  fRejectDAQincomplete = true;

  if (!fOverrideAutoTriggerMask) fTriggerMask = AliVEvent::kINT7;

  /// p-Pb MC do not have the ZDC, the following line avoid any crashes.
  if (fMC) {
    ::Warning("AliEventCuts::SetupRun1pA","Be careful, in p-Pb MC ZDC is not available, if you require the ZN centrality you may encounter a crash.");
  }
  fUseEstimatorsCorrelationCut = false;

  if (fCentralityFramework == 1)
    ::Fatal("AliEventCuts::SetupRun1pA","You cannot use the new centrality framework in Run 1 p-Pb. Please set the fCentralityFramework to 0 to disable the multiplicity selection or to 2 to use AliCentrality.");
  else if (fCentralityFramework == 2) {
    fCentEstimators[0] = iPeriod ? "ZNC" : "ZNA";
    fCentEstimators[1] = iPeriod ? "V0C" : "V0A";
  }
}

void AliEventCuts::SetupRun2pA(int iPeriod) {
  /// iPeriod: 0 p-Pb 5&8 TeV, 1 Pb-p 8 TeV
  ::Info("AliEventCuts::SetupRun2pA","Event cuts for pA are being set on top of Run2 pp standard selections.");
  SetName("StandardRun2pAEventCuts");  

  /// p--Pb requires nsigma cuts on primary vertex
  fMaxDeltaSpdTrackNsigmaSPD = 20.f;
  fMaxDeltaSpdTrackNsigmaTrack = 40.f;

  if (!fOverrideAutoPileUpCuts) {
    fUseCombinedMVSPDcut = true; // If user specify a value it is not overwritten
    fSPDpileupMinZdist = 0.8;
    fSPDpileupNsigmaZdist = 3.;
    fSPDpileupNsigmaDiamXY = 2.;
    fSPDpileupNsigmaDiamZ = 5.;
    fTrackletBGcut = false;
    fUseMultiplicityDependentPileUpCuts = true;
    fUtils.SetMaxPlpChi2MV(5);
    fUtils.SetMinWDistMV(15);
    fUtils.SetCheckPlpFromDifferentBCMV(kFALSE);
  }

  fMinVtz = -10.f;
  fMaxVtz = 10.f;
  fMaxDeltaSpdTrackAbsolute = 0.5f;
  fMaxDeltaSpdTrackNsigmaSPD = 1.e14f;
  fMaxDeltaSpdTrackNsigmaTrack = 1.e14;
  fMaxResolutionSPDvertex = 0.25f;

  fRejectDAQincomplete = true;

  if (!fOverrideAutoTriggerMask) fTriggerMask = AliVEvent::kINT7;

  /// p-Pb MC do not have the ZDC, the following line avoid any crashes.
  if (fMC) {
    ::Warning("AliEventCuts::SetupRun2pA","Be careful, in p-Pb MC ZDC is not available, if you require the ZN centrality you may encounter a crash.");
  }
  fUseEstimatorsCorrelationCut = false;

  if (fCentralityFramework > 1)
    ::Fatal("AliEventCuts::SetupRun2pA","You cannot use the legacy centrality framework in pp. Please set the fCentralityFramework to 0 to disable the multiplicity selection or to 1 to use AliMultSelection.");
  else if (fCentralityFramework == 1) {
    fCentEstimators[0] = iPeriod ? "ZNC" : "ZNA";
    fCentEstimators[1] = iPeriod ? "V0C" : "V0A";
  }

}

void  AliEventCuts::OverridePileUpCuts(int minContrib, float minZdist, float nSigmaZdist, float nSigmaDiamXY, float nSigmaDiamZ, bool ov) {
  fSPDpileupMinContributors = minContrib;
  fSPDpileupMinZdist = minZdist;
  fSPDpileupNsigmaZdist = nSigmaZdist;
  fSPDpileupNsigmaDiamXY = nSigmaDiamXY;
  fSPDpileupNsigmaDiamZ = nSigmaDiamZ;
  fOverrideAutoPileUpCuts = ov;
}

bool AliEventCuts::GoodPrimaryAODVertex(AliVEvent* ev) {
  AliAODEvent* aodEv = dynamic_cast<AliAODEvent*>(ev);
  if (!aodEv) {
    ::Fatal("AliEventCuts::GoodPrimaryAODVertex","Passed argument is not an AOD event.");
  }
  if (!aodEv->GetPrimaryVertex()) return kFALSE;
  if (aodEv->GetPrimaryVertex()->GetType()!=AliAODVertex::kPrimary) return kFALSE;
  const AliAODVertex *vtPrim = aodEv->GetPrimaryVertex();
  const AliAODVertex *vtTPC  = aodEv->GetPrimaryVertexTPC();
  if (std::abs(vtPrim->GetZ()-vtTPC->GetZ())<1e-6 && 
      std::abs(vtPrim->GetChi2perNDF()-vtTPC->GetChi2perNDF())<1e-6) {
      ::Warning("AliEventCuts::GoodPrimaryAODVertex","TPC vertex used as primary");
      return false;
  }
  return true;
}

void AliEventCuts::UseMultSelectionEventSelection(bool useIt) {
  fMultSelectionEvCuts = useIt;
  if (useIt) {
    SetCentralityRange(0.,100.);
  }
}

void AliEventCuts::SetupLHC17n() {
  ::Info("AliEventCuts::AutomaticSetup","Xe-Xe runs found: we will setup the same LHC15o cuts with somewhat different correlation cuts settings (switched off by default).");
  SetName("StandardLHC17nEventCuts");
  SetupRun2PbPb();
  fUseEstimatorsCorrelationCut = false;

  array<double,5> vzero_tpcout_polcut = {-1000.0, 2.8, 1.2e-5,0.,0.};
  std::copy(vzero_tpcout_polcut.begin(),vzero_tpcout_polcut.end(),fVZEROvsTPCoutPolCut);

}

bool AliEventCuts::IsTrueINELgtZero(AliVEvent *ev, bool chkGenVtxZ) {
  if (!ev) return false;
  if (dynamic_cast<AliAODEvent*>(ev)) {
    TClonesArray *stack = (TClonesArray*)ev->GetList()->FindObject(AliAODMCParticle::StdBranchName());
    if (!stack) {
      AliWarning("MC particles array not found...");
      return false;
    }
    AliAODMCHeader* mcHeader = (AliAODMCHeader*)ev->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    if (!mcHeader) {
      AliWarning("AliAODHeader not found...");
      return false;
    }
    if (chkGenVtxZ && !(mcHeader->GetVtxZ() >= fMinVtz &&  mcHeader->GetVtxZ() <= fMaxVtz)) return false;
    for (int i_mcTrk = 0; i_mcTrk < stack->GetEntriesFast(); ++i_mcTrk) {
      AliAODMCParticle* aliMCPart = dynamic_cast<AliAODMCParticle*>(stack->At(i_mcTrk));
      if (!aliMCPart) continue;
      if (aliMCPart->Charge() == -99) continue;                      // dummy value for particles w/o TParticlePDG informations
      if ((double)TMath::Abs(aliMCPart->Charge()/3.) < 0.001) continue; // consider only charged tracks
      if (!aliMCPart->IsPhysicalPrimary()) continue;                 // consider only physical primary particle
      if (TMath::Abs(aliMCPart->Eta()) > 1) continue;                // out acceptances
      return true;
    }
  } else if (!dynamic_cast<AliESDEvent*>(ev))
    AliFatal("I don't find the AOD event nor the ESD one, aborting.");
  else {
    AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
    AliMCEvent* mcEvent = eventHandler->MCEvent();
    const AliVVertex* gVtx = mcEvent->GetPrimaryVertex();
    if (chkGenVtxZ && !(gVtx->GetZ() >= fMinVtz && gVtx->GetZ() <= fMaxVtz)) return false;
    for (int i_mcTrk = 0; i_mcTrk < mcEvent->GetNumberOfTracks(); ++i_mcTrk) {
      AliMCParticle *aliMCPart = (AliMCParticle *)mcEvent->GetTrack(i_mcTrk);
      if (!aliMCPart) continue;
      if (aliMCPart->Charge() == -99) continue;                         // dummy value for particles w/o TParticlePDG informations
      if ((double)TMath::Abs(aliMCPart->Charge()/3.) < 0.001) continue; // consider only charged particle
      if (!mcEvent->IsPhysicalPrimary(i_mcTrk)) continue;               // consider only physical primary particle
      if (TMath::Abs(aliMCPart->Eta()) > 1) continue;                   // out acceptances
      return true;
    }
  }
  return false;
}

/// Different trigger classes are comma separated
void AliEventCuts::SetAcceptedTriggerClasses(TString classes) {
  TObjArray* classArray = classes.Tokenize(",");
  for (int iClass = 0; iClass < classArray->GetSize(); ++iClass) {
    TObjString* thisClass = static_cast<TObjString*>(classArray->At(iClass));
    fTriggerClasses.push_back(thisClass->GetString().Data());
  }
}

