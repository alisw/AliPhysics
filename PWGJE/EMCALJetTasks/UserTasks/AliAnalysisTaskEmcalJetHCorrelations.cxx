//////////
//Measure Jet-hadron correlations
//Does event Mixing using AliEventPoolManager
/////////

#include "AliAnalysisTaskEmcalJetHCorrelations.h"

#include <bitset>

#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THnSparse.h>
#include <TVector3.h>
#include <TFile.h>
#include <TGrid.h>
#include <TRandom3.h>

#include <AliAnalysisManager.h>
#include <AliInputEventHandler.h>
#include <AliEventPoolManager.h>
#include <AliLog.h>
#include <AliVAODHeader.h>
#include <AliVTrack.h>

#include "AliEmcalJet.h"
#include "AliTLorentzVector.h"
#include "AliBasicParticle.h"
#include "AliEmcalContainerUtils.h"
#include "AliClusterContainer.h"
#include "AliTrackContainer.h"
#include "AliJetContainer.h"
#include "AliAnalysisTaskEmcalEmbeddingHelper.h"
#include "AliAnalysisTaskEmcalJetHUtils.h"

ClassImp(PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalJetHCorrelations);

namespace PWGJE {
namespace EMCALJetTasks {

// 0-10% centrality: Semi-Good Runs
Double_t AliAnalysisTaskEmcalJetHCorrelations::p0_10SG[17] = {0.906767, 0.0754127, 1.11638, -0.0233078, 0.795454, 0.00935385, -0.000327857, 1.08903, 0.0107272, 0.443252, -0.143411, 0.965822, 0.359156, -0.581221, 1.0739, 0.00632828, 0.706356};
// 10-30% centrality: Semi-Good Runs
Double_t AliAnalysisTaskEmcalJetHCorrelations::p10_30SG[17] = {0.908011, 0.0769254, 1.11912, -0.0249449, 0.741488, 0.0361252, -0.00367954, 1.10424, 0.011472, 0.452059, -0.133282, 0.980633, 0.358222, -0.620256, 1.06871, 0.00564449, 0.753168};
// 30-50% centrality: Semi-Good Runs
Double_t AliAnalysisTaskEmcalJetHCorrelations::p30_50SG[17] = {0.958708, 0.0799197, 1.10817, -0.0357678, 0.75051, 0.0607808, -0.00929713, 0.998801, 0.00692244, 0.615452, -0.0480328, 0.968431, 0.321634, -0.619066, 1.03412, 0.00656201, 0.798666};
// 50-90% centrality: Semi-Good Runs
Double_t AliAnalysisTaskEmcalJetHCorrelations::p50_90SG[17] = {0.944565, 0.0807258, 1.12709, -0.0324746, 0.666452, 0.0842476, -0.00963837, 1.02829, 0.00666852, 0.549625, -0.0603107, 0.981374, 0.309374, -0.619181, 1.05367, 0.005925, 0.744887};

// 0-10% centrality: Good Runs
Double_t AliAnalysisTaskEmcalJetHCorrelations::p0_10G[17] = {0.971679, 0.0767571, 1.13355, -0.0274484, 0.856652, 0.00536795, 3.90795e-05, 1.06889, 0.011007, 0.447046, -0.146626, 0.919777, 0.192601, -0.268515, 1.00243, 0.00620849, 0.709477};
// 10-30% centrality: Good Runs
Double_t AliAnalysisTaskEmcalJetHCorrelations::p10_30G[17] = {0.97929, 0.0776039, 1.12213, -0.0300645, 0.844722, 0.0134788, -0.0012333, 1.07955, 0.0116835, 0.456608, -0.132743, 0.930964, 0.174175, -0.267154, 0.993118, 0.00574892, 0.765256};
// 30-50% centrality: Good Runs
Double_t AliAnalysisTaskEmcalJetHCorrelations::p30_50G[17] = {0.997696, 0.0816769, 1.14341, -0.0353734, 0.752151, 0.0744259, -0.0102926, 1.01561, 0.00713274, 0.57203, -0.0640248, 0.947747, 0.102007, -0.194698, 0.999164, 0.00568476, 0.7237};
// 50-90% centrality: Good Runs
Double_t AliAnalysisTaskEmcalJetHCorrelations::p50_90G[17] = {0.97041, 0.0813559, 1.12151, -0.0368797, 0.709327, 0.0701501, -0.00784043, 1.06276, 0.00676173, 0.53607, -0.0703117, 0.982534, 0.0947881, -0.18073, 1.03229, 0.00580109, 0.737801};

/**
 * Default constructor.
 */
AliAnalysisTaskEmcalJetHCorrelations::AliAnalysisTaskEmcalJetHCorrelations() :
  AliAnalysisTaskEmcalJet("AliAnalysisTaskEmcalJetHCorrelations", kFALSE),
  fTrackBias(5),
  fClusterBias(5),
  fDoEventMixing(kFALSE),
  fNMixingTracks(50000), fMinNTracksMixedEvents(5000), fMinNEventsMixedEvents(5), fNCentBinsMixedEvent(10),
  fPoolMgr(nullptr),
  fTriggerType(AliVEvent::kEMCEJE), fMixingEventType(AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral),
  fDisableFastPartition(kFALSE),
  fSingleTrackEfficiencyCorrectionType(AliAnalysisTaskEmcalJetHCorrelations::kEffDisable),
  fArtificialTrackInefficiency(1.0),
  fNoMixedEventJESCorrection(kFALSE),
  fJESCorrectionHist(nullptr),
  fDoLessSparseAxes(kFALSE), fDoWiderTrackBin(kFALSE),
  fRequireMatchedJetWhenEmbedding(kTRUE),
  fMinSharedMomentumFraction(0.),
  fRequireMatchedPartLevelJet(false),
  fMaxMatchedJetDistance(-1),
  fHistManager(),
  fHistTrackPt(nullptr),
  fHistJetEtaPhi(nullptr),
  fHistJetHEtaPhi(nullptr),
  fhnMixedEvents(nullptr),
  fhnJH(nullptr),
  fhnTrigger(nullptr)
{
  // Default Constructor
  InitializeArraysToZero();
}

/**
 * Standard constructor
 */
AliAnalysisTaskEmcalJetHCorrelations::AliAnalysisTaskEmcalJetHCorrelations(const char *name) :
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fTrackBias(5),
  fClusterBias(5),
  fDoEventMixing(kFALSE),
  fNMixingTracks(50000), fMinNTracksMixedEvents(5000), fMinNEventsMixedEvents(5), fNCentBinsMixedEvent(10),
  fPoolMgr(nullptr),
  fTriggerType(AliVEvent::kEMCEJE), fMixingEventType(AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral),
  fDisableFastPartition(kFALSE),
  fSingleTrackEfficiencyCorrectionType(AliAnalysisTaskEmcalJetHCorrelations::kEffDisable),
  fArtificialTrackInefficiency(1.0),
  fNoMixedEventJESCorrection(kFALSE),
  fJESCorrectionHist(nullptr),
  fDoLessSparseAxes(kFALSE), fDoWiderTrackBin(kFALSE),
  fRequireMatchedJetWhenEmbedding(kTRUE),
  fMinSharedMomentumFraction(0.),
  fRequireMatchedPartLevelJet(false),
  fMaxMatchedJetDistance(-1),
  fHistManager(name),
  fHistTrackPt(nullptr),
  fHistJetEtaPhi(nullptr),
  fHistJetHEtaPhi(nullptr),
  fhnMixedEvents(nullptr),
  fhnJH(nullptr),
  fhnTrigger(nullptr)
{
  // Constructor
  InitializeArraysToZero();
  // Ensure that additional general histograms are created
  SetMakeGeneralHistograms(kTRUE);
}

/**
 * Initialize all member arrays to nullptr. Used as a convenience function in the constructors.
 */
void AliAnalysisTaskEmcalJetHCorrelations::InitializeArraysToZero()
{
  for(Int_t trackPtBin = 0; trackPtBin < kMaxTrackPtBins; trackPtBin++){
    fHistTrackEtaPhi[trackPtBin] = nullptr;
  }
  for(Int_t centralityBin = 0; centralityBin < kMaxCentralityBins; ++centralityBin){
    fHistJetPt[centralityBin] = nullptr;
    fHistJetPtBias[centralityBin] = nullptr;
  }
}

/**
 * Perform run independent initializations, such as histograms and the event pool.
 */
void AliAnalysisTaskEmcalJetHCorrelations::UserCreateOutputObjects() {
  // Called once 
  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  // Create histograms
  fHistTrackPt = new TH1F("fHistTrackPt", "P_{T} distribution", 1000, 0.0, 100.0);
  fHistJetEtaPhi = new TH2F("fHistJetEtaPhi","Jet eta-phi",900,-1.8,1.8,720,-3.2,3.2);
  fHistJetHEtaPhi = new TH2F("fHistJetHEtaPhi","Jet-Hadron deta-dphi",900,-1.8,1.8,720,-1.6,4.8);

  fOutput->Add(fHistTrackPt);
  fOutput->Add(fHistJetEtaPhi);
  fOutput->Add(fHistJetHEtaPhi);

  TString name;
  for(Int_t trackPtBin = 0; trackPtBin < kMaxTrackPtBins; ++trackPtBin){
    name = Form("fHistTrackEtaPhi_%i", trackPtBin);
    fHistTrackEtaPhi[trackPtBin] = new TH2F(name,name,400,-1,1,720,0.0,2.0*TMath::Pi());
    fOutput->Add(fHistTrackEtaPhi[trackPtBin]);
  }

  for(Int_t centralityBin = 0; centralityBin < kMaxCentralityBins; ++centralityBin){
    name = Form("fHistJetPt_%i",centralityBin);
    fHistJetPt[centralityBin] = new TH1F(name,name,200,0,200);
    fOutput->Add(fHistJetPt[centralityBin]);

    name = Form("fHistJetPtBias_%i",centralityBin);
    fHistJetPtBias[centralityBin] = new TH1F(name,name,200,0,200);
    fOutput->Add(fHistJetPtBias[centralityBin]);
  }

  // Jet matching cuts
  std::vector<std::string> binLabels = {"noMatch", "matchedJet", "sharedMomentumFraction", "partLevelMatchedJet", "jetDistance", "passedAllCuts"};
  for (auto histName : std::vector<std::string>({"SameEvent", "MixedEvent"})) {
    name = std::string("fHistJetMatching") + histName.c_str() + "Cuts";
    std::string title = std::string("Jets which passed matching jet cuts for ") + histName;
    auto histMatchedJetCuts = fHistManager.CreateTH1(name, title.c_str(), binLabels.size(), 0, binLabels.size());
    // Set label names
    for (unsigned int i = 1; i <= binLabels.size(); i++) {
      histMatchedJetCuts->GetXaxis()->SetBinLabel(i, binLabels.at(i-1).c_str());
    }
    histMatchedJetCuts->GetYaxis()->SetTitle("Number of jets");
  }

  UInt_t cifras = 0; // bit coded, see GetDimParams() below 
  if(fDoLessSparseAxes) {
    cifras = 1<<0 | 1<<1 | 1<<2 | 1<<3 | 1<<4 | 1<<5 | 1<<9;
  } else {
    cifras = 1<<0 | 1<<1 | 1<<2 | 1<<3 | 1<<4 | 1<<5 | 1<<7 | 1<<9;
  }
  fhnJH = NewTHnSparseF("fhnJH", cifras);
  fhnJH->Sumw2();
  fOutput->Add(fhnJH);

  if(fDoEventMixing){    
    // The event plane angle does not need to be included because the semi-central determined that the EP angle didn't change
    // significantly for any of the EP orientations. However, it will be included so this can be demonstrated for the central
    // analysis if so desired.
    if(fDoLessSparseAxes) { 
      cifras = 1<<0 | 1<<1 | 1<<2 | 1<<3 | 1<<4 | 1<<5 | 1<<9;
    } else {
      cifras = 1<<0 | 1<<1 | 1<<2 | 1<<3 | 1<<4 | 1<<5 | 1<<7 | 1<<9;
    }
    fhnMixedEvents = NewTHnSparseF("fhnMixedEvents", cifras);
    fhnMixedEvents->Sumw2();
    fOutput->Add(fhnMixedEvents);
  }

  // Trigger THnSparse
  cifras = 1<<0 | 1<<1 | 1<<9;
  fhnTrigger = NewTHnSparseF("fhnTrigger", cifras);
  fhnTrigger->Sumw2();
  fOutput->Add(fhnTrigger);

  // Store hist manager output in the output list
  TIter next(fHistManager.GetListOfHistograms());
  TObject* obj = 0;
  while ((obj = next())) {
    fOutput->Add(obj);
  }
  
  PostData(1, fOutput);

  // Event Mixing
  Int_t poolSize = -1;  // Maximum number of events. Set to -1 to avoid limits on number of events
  // ZVertex
  Int_t nZVertexBins = 10;
  Double_t* zVertexBins = GenerateFixedBinArray(nZVertexBins, -10, 10);
  // Event activity (centrality of multiplicity)
  Int_t nEventActivityBins = 8;
  Double_t* eventActivityBins = 0;
  // +1 to accomodate the fact that we define bins rather than array entries.
  Double_t multiplicityBins[kMixedEventMultiplicityBins+1] = {0., 4., 9., 15., 25., 35., 55., 100., 500.};

  // Cannot use GetBeamType() since it is not available until UserExec()
  if (fForceBeamType != AliAnalysisTaskEmcal::kpp ) {   //all besides pp
    // Event Activity is centrality in AA, pA
    nEventActivityBins = fNCentBinsMixedEvent;
    eventActivityBins = GenerateFixedBinArray(nEventActivityBins, 0, 100);
  }
  else if (fForceBeamType == AliAnalysisTaskEmcal::kpp) { //for pp only
    // Event Activity is multiplicity in pp
    eventActivityBins = multiplicityBins;
  }

  fPoolMgr = new AliEventPoolManager(poolSize, fNMixingTracks, nEventActivityBins, eventActivityBins, nZVertexBins, zVertexBins);
}

void AliAnalysisTaskEmcalJetHCorrelations::ExecOnce()
{
  if (fArtificialTrackInefficiency < 1.) {
    if (gRandom) delete gRandom;
    gRandom = new TRandom3(0);
  }

  // Call the base class
  AliAnalysisTaskEmcalJet::ExecOnce();
}

/**
 * Get the proper bin based on the track pt value.
 *
 * @param[in] pt Track pt value to be binned.
 * @return Bin corresponding to the input value.
 */
Int_t AliAnalysisTaskEmcalJetHCorrelations::GetTrackPtBin(Double_t pt) const
{
  Int_t ptBin = -1;
  if      (pt <  0.5) ptBin = 0;
  else if (pt <  1  ) ptBin = 1;
  else if (pt <  2  ) ptBin = 2;
  else if (pt <  3  ) ptBin = 3;
  else if (pt <  5  ) ptBin = 4;
  else if (pt <  8  ) ptBin = 5;
  else if (pt < 20  ) ptBin = 6;

  return ptBin;
}

/**
 * Retrieve the trigger mask from the input event or embedding helper as approapriate.
 *
 * @return The event trigger mask.
 */
UInt_t AliAnalysisTaskEmcalJetHCorrelations::RetrieveTriggerMask() const
{
  UInt_t eventTrigger = 0;
  if (fIsEmbedded) {
    auto embeddingHelper = AliAnalysisTaskEmcalEmbeddingHelper::GetInstance();
    if (embeddingHelper) {
      auto aodHeader = dynamic_cast<AliVAODHeader *>(embeddingHelper->GetEventHeader());
      if (aodHeader) {
        AliDebugStream(5) << "Retrieving trigger mask from embedded event\n";
        eventTrigger = aodHeader->GetOfflineTrigger();
      }
      else {
        AliErrorStream() << "Failed to retrieve requested AOD header from embedding helper\n";
      }
    }
    else {
      AliErrorStream() << "Failed to retrieve requested embedding helper\n";
    }
  }
  else {
    AliDebugStream(5) << "Retrieving trigger mask from internal event\n";
    eventTrigger = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  }

  return eventTrigger;
}

/**
 * Main loop called for each event by AliAnalysisTaskEmcal.
 */
Bool_t AliAnalysisTaskEmcalJetHCorrelations::Run()
{
  // NOTE: Clusters are never used directly in the task, so the container is neither created not retrieved
  // Retrieve tracks
  AliTrackContainer * tracks = static_cast<AliTrackContainer * >(GetParticleContainer("tracksForCorrelations"));
  if (!tracks) {
    AliError(Form("%s: Unable to retrieve tracks!", GetName()));
    return kFALSE;
  }

  // Retrieve jets
  AliJetContainer * jets = GetJetContainer(0);
  if (!jets) {
    AliError(Form("%s: Unable to retrieve jets!", GetName()));
    return kFALSE;
  }

  // Keep track of the tracks which are rejected with an aritificial track inefficiency
  std::vector<unsigned int> rejectedTrackIndices;
  bool useListOfRejectedIndices = false;

  // Get z vertex
  Double_t zVertex=fVertex[2];
  // Flags
  Bool_t biasedJet = kFALSE;
  Bool_t leadJet = kFALSE;
  // Relative angles and distances
  Double_t deltaPhi = 0;
  Double_t deltaEta = 0;
  Double_t deltaR = 0;
  Double_t epAngle = 0;
  // Event activity (centrality or multipilicity)
  Double_t eventActivity = 0;
  // Efficiency correction
  Double_t efficiency = -999;
  // For comparison to the current jet
  AliEmcalJet * leadingJet = jets->GetLeadingJet();
  // For getting the proper properties of tracks
  AliTLorentzVector track;

  // Determine the trigger for the current event
  UInt_t eventTrigger = RetrieveTriggerMask();

  AliDebugStream(5) << "Beginning main processing. Number of jets: " << jets->GetNJets() << ", accepted jets: " << jets->GetNAcceptedJets() << "\n";

  // Handle fast partition if selected
  if ((eventTrigger & AliVEvent::kFastOnly) && fDisableFastPartition) {
    AliDebugStream(4) << GetName() << ": Fast partition disabled\n";
    if (fGeneralHistograms) {
      fHistEventRejection->Fill("Fast Partition", 1);
    }
    return kFALSE;
  }

  for (auto jet : jets->accepted()) {
    // Selects only events that we are interested in (ie triggered)
    if (!(eventTrigger & fTriggerType)) {
      // The two bitwise and to zero yet are still equal when both are 0, so we allow for that possibility
      if (eventTrigger == fTriggerType && eventTrigger == 0) {
        AliDebugStream(5) << "Event accepted because the physics selection is \"0\".\n";
      }
      else {
        AliDebugStream(5) << "Rejected jets due to physics selection. Phys sel: " << std::bitset<32>(eventTrigger) << ", requested triggers: " << std::bitset<32>(fTriggerType) << " \n";
        // We can break here - the physics selection is not going to change within an event.
        break;
      }
    }

    AliDebugStream(5) << "Jet passed event selection!\nJet: " << jet->toString().Data() << "\n";

    // Require the found jet to be matched
    // This match should be between detector and particle level MC
    if (fIsEmbedded && fRequireMatchedJetWhenEmbedding) {
      bool foundMatchedJet = CheckForMatchedJet(jets, jet, "fHistJetMatchingSameEventCuts");
      if (foundMatchedJet == false) {
        continue;
      }
    }

    // Determine event activity
    if (fBeamType == kAA || fBeamType == kpA) {
      eventActivity = fCent;
    }
    else if (fBeamType == kpp) {
      eventActivity = static_cast<Double_t>(tracks->GetNTracks());
    }

    // Jet properties
    // Determine if we have the lead jet
    leadJet = kFALSE;
    if (jet == leadingJet) leadJet = kTRUE;
    biasedJet = BiasedJet(jet);
    epAngle = PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalJetHUtils::RelativeEPAngle(jet->Phi(), fEPV0);

    // Fill jet properties
    fHistJetEtaPhi->Fill(jet->Eta(), jet->Phi());
    FillHist(fHistJetPt[fCentBin], jet->Pt());
    if (biasedJet == kTRUE) {
      FillHist(fHistJetPtBias[fCentBin], jet->Pt());

      const double triggerInfo[] = {eventActivity, jet->Pt(), epAngle};
      fhnTrigger->Fill(triggerInfo);
    }

    // Cut on jet pt of 15 to reduce the size of the sparses
    if (jet->Pt() > 15) {

      AliDebugStream(4) << "Passed min jet pt cut of 15. Jet: " << jet->toString().Data() << "\n";
      AliDebugStream(4) << "N accepted tracks: " << tracks->GetNAcceptedTracks() << "\n";
      auto tracksIter = tracks->accepted_momentum();
      for (auto trackIter = tracksIter.begin(); trackIter != tracksIter.end(); trackIter++ ) {
        // Get proper track properties
        track.Clear();
        track = trackIter->first;

        // Artificial inefficiency
        // Note that we already randomly rejected tracks so that the same tracks will be rejected for the mixed events
        bool rejectParticle = CheckArtificialTrackEfficiency(trackIter.current_index(), rejectedTrackIndices, useListOfRejectedIndices);
        if (rejectParticle) {
          AliDebugStream(4) << "Track rejected in signal correlation loop.\n";
          continue;
        }

        // Determine relative angles and distances and set the respective variables
        GetDeltaEtaDeltaPhiDeltaR(track, jet, deltaEta, deltaPhi, deltaR);

        // Fill track properties
        fHistTrackPt->Fill(track.Pt());
        fHistJetHEtaPhi->Fill(deltaEta, deltaPhi);

        // Calculate single particle tracking efficiency for correlations
        efficiency = EffCorrection(track.Eta(), track.Pt());
        AliDebugStream(6) << "track eta: " << track.Eta() << ", track pt: " << track.Pt() << ", efficiency: " << efficiency << "\n";

        if (biasedJet == kTRUE) {
          if(fDoLessSparseAxes) { // check if we want all dimensions
            double triggerEntries[] = {eventActivity, jet->Pt(), track.Pt(), deltaEta, deltaPhi, static_cast<Double_t>(leadJet), epAngle};
            FillHist(fhnJH, triggerEntries, 1.0/efficiency);
          } else { 
            double triggerEntries[] = {eventActivity, jet->Pt(), track.Pt(), deltaEta, deltaPhi, static_cast<Double_t>(leadJet), deltaR, epAngle};
            FillHist(fhnJH, triggerEntries, 1.0/efficiency);
          }
        }

      } //track loop

      // After one jet (and looping over whatever tracks are available in this event), we want to use the list of rejected indices,
      // both for the next possible signal jet in the event and for the mixed events
      AliDebugStream(4) << "Switching to list of rejected track indices. Number of indices: " << rejectedTrackIndices.size() << "\n";
      useListOfRejectedIndices = true;

    }//jet pt cut
  }//jet loop

  //Prepare to do event mixing

  // create a list of reduced objects. This speeds up processing and reduces memory consumption for the event pool
  TObjArray* tracksClone = 0;

  if(fDoEventMixing == kTRUE){

    // event mixing

    // 1. First get an event pool corresponding in mult (cent) and
    //    zvertex to the current event. Once initialized, the pool
    //    should contain nMix (reduced) events. This routine does not
    //    pre-scan the chain. The first several events of every chain
    //    will be skipped until the needed pools are filled to the
    //    specified depth. If the pool categories are not too rare, this
    //    should not be a problem. If they are rare, you could lose
    //    statistics.

    // 2. Collect the whole pool's content of tracks into one TObjArray
    //    (bgTracks), which is effectively a single background super-event.

    // 3. The reduced and bgTracks arrays must both be passed into
    //    FillCorrelations(). Also nMix should be passed in, so a weight
    //    of 1./nMix can be applied.

    AliEventPool *pool = 0;
    if (fBeamType == kAA || fBeamType == kpA) {//everything but pp
      pool = fPoolMgr->GetEventPool(fCent, zVertex);
    }
    else if (fBeamType == kpp) {//pp only
      pool = fPoolMgr->GetEventPool(static_cast<Double_t>(tracks->GetNTracks()), zVertex);
    }

    if (!pool){
      if (fBeamType == kAA || fBeamType == kpA) AliFatal(Form("No pool found for centrality = %f, zVertex = %f", fCent, zVertex));
      else if (fBeamType == kpp) AliFatal(Form("No pool found for ntracks_pp = %d, zVertex = %f", tracks->GetNTracks(), zVertex));
      return kTRUE;
    }

    // The number of events in the pool
    Int_t nMix = pool->GetCurrentNEvents();

    // The two bitwise and to zero yet are still equal when both are 0, so we allow for that possibility
    if((eventTrigger & fTriggerType) || eventTrigger == fTriggerType) {
      // check for a trigger jet
      if (pool->IsReady() || pool->NTracksInPool() >= fMinNTracksMixedEvents || nMix >= fMinNEventsMixedEvents) {

        for (auto jet : jets->accepted()) {
          // Require the found jet to be matched
          // This match should be between detector and particle level MC
          if (fIsEmbedded && fRequireMatchedJetWhenEmbedding) {
            bool foundMatchedJet = CheckForMatchedJet(jets, jet, "fHistJetMatchingMixedEventCuts");
            if (foundMatchedJet == false) {
              continue;
            }
          }

          if (fBeamType == kAA || fBeamType == kpA) { //pA and AA
            eventActivity = fCent;
          }
          else if (fBeamType == kpp) {
            eventActivity = static_cast<Double_t>(tracks->GetNTracks());
          }

          // Jet properties
          // Determine if we have the lead jet
          leadJet = kFALSE;
          if (jet == leadingJet) { leadJet = kTRUE; }
          biasedJet = BiasedJet(jet);
          epAngle = PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalJetHUtils::RelativeEPAngle(jet->Phi(), fEPV0);

          // Make sure event contains a biased jet above our threshold (reduce stats of sparse)
          if (jet->Pt() < 15 || biasedJet == kFALSE) continue;

          // Fill mixed-event histos here
          for (Int_t jMix=0; jMix < nMix; jMix++) {
            TObjArray* bgTracks = pool->GetEvent(jMix);

            for (Int_t ibg=0; ibg < bgTracks->GetEntries(); ibg++){
              AliBasicParticle *bgTrack = static_cast<AliBasicParticle*>(bgTracks->At(ibg));
              if(!bgTrack) {
                AliError(Form("%s:Failed to retrieve tracks from mixed events", GetName()));
              }

              // NOTE: We don't need to apply the artificial track inefficiency here because we already applied
              //       it when will filling into the event pool (in CloneAndReduceTrackList()).

              // Fill into TLorentzVector for use with functions below
              track.Clear();
              track.SetPtEtaPhiE(bgTrack->Pt(), bgTrack->Eta(), bgTrack->Phi(), 0);

              // Calculate single particle tracking efficiency of mixed events for correlations
              efficiency = EffCorrection(track.Eta(), track.Pt());

              // Phi is [-0.5*TMath::Pi(), 3*TMath::Pi()/2.]
              GetDeltaEtaDeltaPhiDeltaR(track, jet, deltaEta, deltaPhi, deltaR);

              if (fDoLessSparseAxes) {  // check if we want all the axis filled
                double triggerEntries[] = {eventActivity, jet->Pt(), track.Pt(), deltaEta, deltaPhi, static_cast<Double_t>(leadJet), epAngle};
                FillHist(fhnMixedEvents, triggerEntries, 1./(nMix*efficiency), fNoMixedEventJESCorrection);
              } else {
                double triggerEntries[] = {eventActivity, jet->Pt(), track.Pt(), deltaEta, deltaPhi, static_cast<Double_t>(leadJet), deltaR, epAngle};
                FillHist(fhnMixedEvents, triggerEntries, 1./(nMix*efficiency), fNoMixedEventJESCorrection);
              }
            }
          }
        }
      }
    }

    // The two bitwise and to zero yet are still equal when both are 0, so we allow for that possibility
    if ((eventTrigger & fMixingEventType) || eventTrigger == fMixingEventType) {
      tracksClone = CloneAndReduceTrackList(rejectedTrackIndices, useListOfRejectedIndices);

      //update pool if jet in event or not
      pool->UpdatePool(tracksClone);
    }

  } // end of event mixing

  return kTRUE;
}      

/**
 * Determine if a jet passes the track or cluster bias and is therefore a "biased" jet.
 *
 * @param[in] jet The jet to be checked.
 * @return true if the jet is biased.
 */
Bool_t AliAnalysisTaskEmcalJetHCorrelations::BiasedJet(AliEmcalJet * jet)
{
  if ((jet->MaxTrackPt() > fTrackBias) || (jet->MaxClusterPt() > fClusterBias))
  {
    return kTRUE;
  }
  return kFALSE;
}

/**
 * Get \f$\Delta\phi\f$, \f$\Delta\eta\f$, and \f$\Delta R\f$ between two given particles.
 *
 * @param[in] particleOne The first particle.
 * @param[in] particleTwo The second particle.
 * @param[out] deltaEta The calculated \f$\Delta\eta\f$ value.
 * @param[out] deltaPhi The calculated \f$\Delta\phi\f$ value.
 * @param[out] deltaR The calculated \f$\Delta R\f$ value.
 */
void AliAnalysisTaskEmcalJetHCorrelations::GetDeltaEtaDeltaPhiDeltaR(AliTLorentzVector & particleOne, AliVParticle * particleTwo, Double_t & deltaEta, Double_t & deltaPhi, Double_t & deltaR)
{
  // Define dPhi = jet.Phi() - particle.Phi() and similarly for dEta
  // Returns deltaPhi in symmetric range so that we can calculate DeltaR.
  deltaPhi = DeltaPhi(particleOne.Phi(), particleTwo->Phi(), -1.0*TMath::Pi(), TMath::Pi());
  deltaEta = particleTwo->Eta() - particleOne.Eta();
  deltaR = TMath::Sqrt(deltaPhi*deltaPhi + deltaEta*deltaEta);

  // Adjust to the normal range after the DeltaR caluclation
  deltaPhi = DeltaPhi(particleTwo->Phi(), particleOne.Phi(), -0.5*TMath::Pi(), 3*TMath::Pi()/2.);
}

/**
 * Check whether a track should be rejected due to artificial track inefficiency.
 *
 * @param[in] trackIndex Index of the current track.
 * @param[in] rejectedTrackIndices Vector of track indices which have been already rejected. Can be used to store the track index if rejected or reject a track if matched (depending on useRejectedList).
 * @param[in] useRejectedList If true, check if the current track index is in the rejected index list and reject if so. If false, will randomly reject the track according to the track inefficiency.
 * @return True if particle should be rejected
 */
bool AliAnalysisTaskEmcalJetHCorrelations::CheckArtificialTrackEfficiency(unsigned int trackIndex, std::vector<unsigned int> & rejectedTrackIndices, bool useRejectedList)
{
  bool returnValue = false;
  if (fArtificialTrackInefficiency < 1.) {
    if (useRejectedList) {
      if (std::find(rejectedTrackIndices.begin(), rejectedTrackIndices.end(), trackIndex) != rejectedTrackIndices.end()) {
        AliDebugStream(4) << "Track " << trackIndex << " rejected due to artificial tracking inefficiency (from list)\n";
        returnValue = true;
      }
    }
    else {
      // Rejet randomly
      Double_t rnd = gRandom->Rndm();
      if (fArtificialTrackInefficiency < rnd) {
        // Store index so we can reject it again if it is also filled for mixed events
        rejectedTrackIndices.push_back(trackIndex);
        AliDebugStream(4) << "Track " << trackIndex << " rejected due to artificial tracking inefficiency (from random)\n";
        returnValue = true;
      }
    }
  }

  return returnValue;
}

/**
 * Check for whether a matched jet should be accepted based on:
 * - Jet being identified as matched to another jet
 * - The shared momentum fraction being larger than some minimum value
 * - Their matched distance being below the max matching distance
 * - A particle level jet being matched to the detector level jet
 *
 * All of the above options are configurable and off by default (except for requiring a basic match).
 *
 * NOTE: AliEmcalJet::ClosestJet() is called instead of AliEmcalJet::MatchedJet() because ClosestJet() will work
 * with both the EMCal Jet Tagger and the Response Maker, while MatchedJet() will only work with the Response Maker
 * due to the design of the classes.
 *
 * @param[in] jets Jet container corresponding to the matched jet
 * @param[in] jet Jet to be checked
 * @param[in] histName Name of the hist in the hist manager where QA information will be filled
 * @return true if the jet passes the criteria. false otherwise.
 */
bool AliAnalysisTaskEmcalJetHCorrelations::CheckForMatchedJet(AliJetContainer * jets, AliEmcalJet * jet, const std::string & histName)
{
  bool returnValue = false;
  if (jet->ClosestJet()) {
    fHistManager.FillTH1(histName.c_str(), "matchedJet");
    returnValue = true;
    // TODO: Can it be merged with the function in JetHPerformance?
    AliDebugStream(4) << "Jet is matched!\nJet: " << jet->toString() << "\n";
    // Check shared momentum fraction
    // We explicitly want to use indices instead of geometric matching
    double sharedFraction = jets->GetFractionSharedPt(jet, nullptr);
    if (sharedFraction < fMinSharedMomentumFraction) {
      AliDebugStream(4) << "Jet rejected due to shared momentum fraction of " << sharedFraction << ", which is smaller than the min momentum fraction of " << fMinSharedMomentumFraction << "\n";
      returnValue = false;
    }
    else {
      AliDebugStream(4) << "Passed shared momentum fraction with value of " << sharedFraction << "\n";
      fHistManager.FillTH1(histName.c_str(), "sharedMomentumFraction");
    }

    if (fRequireMatchedPartLevelJet) {
      AliEmcalJet * detLevelJet = jet->ClosestJet();
      AliEmcalJet * partLevelJet = detLevelJet->ClosestJet();
      if (!partLevelJet) {
        AliDebugStream(4) << "Jet rejected due to no matching part level jet.\n";
        returnValue = false;
      }
      else {
        AliDebugStream(4) << "Det level jet has a required match to a part level jet.\n" << "Part level jet: " << partLevelJet->toString() << "\n";
        fHistManager.FillTH1(histName.c_str(), "partLevelMatchedJet");
      }
    }

    // Only check matched jet distance if a value has been set
    if (fMaxMatchedJetDistance > 0) {
      double matchedJetDistance = jet->ClosestJetDistance();
      if (matchedJetDistance > fMaxMatchedJetDistance) {
        AliDebugStream(4) << "Jet rejected due to matching distance of " << matchedJetDistance << ", which is larger than the max distance of " << fMaxMatchedJetDistance << "\n";
        returnValue = false;
      }
      else {
        AliDebugStream(4) << "Jet passed distance cut with distance of " << matchedJetDistance << "\n";
        fHistManager.FillTH1(histName.c_str(), "jetDistance");
      }
    }

    // Record all cuts passed
    if (returnValue == true) {
      fHistManager.FillTH1(histName.c_str(), "passedAllCuts");
    }
  }
  else {
    AliDebugStream(5) << "Rejected jet because it was not matched to a external event jet.\n";
    fHistManager.FillTH1(histName.c_str(), "noMatch");
    returnValue = false;
  }

  return returnValue;
}

/**
 * Creates a new THnSparseF based on the desired number of entries. Note that the axes are defined in
 * GetDimParams().
 *
 * @param[in] name Name of the new THnSparseF
 * @param[in] entries Bits corresponding to entires to include in the THnSparseF, with the axis deteremined in GetDimParams().
 */
THnSparse* AliAnalysisTaskEmcalJetHCorrelations::NewTHnSparseF(const char* name, UInt_t entries)
{
  Int_t count = 0;
  UInt_t tmp = entries;
  while(tmp!=0){
    count++;
    tmp = tmp &~ -tmp;  // clear lowest bit
  }

  TString hnTitle(name);
  const Int_t dim = count;
  Int_t nbins[dim];
  Double_t xmin[dim];
  Double_t xmax[dim];

  Int_t i=0;
  Int_t c=0;
  while(c<dim && i<32){
    if(entries&(1<<i)){

      TString label("");
      GetDimParams(i, label, nbins[c], xmin[c], xmax[c]);
      hnTitle += Form(";%s",label.Data());
      c++;
    }

    i++;
  }
  hnTitle += ";";

  return new THnSparseF(name, hnTitle.Data(), dim, nbins, xmin, xmax);
}

/**
 * Stores labels and binning of axes for creating a new THnSparseF.
 *
 * @param[in] iEntry Which axis to pick out of the cases.
 * @param[out] label Label of the axis.
 * @param[out] nbins Number of bins for the axis.
 * @param[out] xmin Minimum value of the axis.
 * @param[out] xmax Maximum value of the axis.
 */
void AliAnalysisTaskEmcalJetHCorrelations::GetDimParams(Int_t iEntry, TString &label, Int_t &nbins, Double_t &xmin, Double_t &xmax)
{
  const Double_t pi = TMath::Pi();

  switch(iEntry){

    case 0:
      label = "V0 centrality (%)";
      nbins = 10;
      xmin = 0.;
      xmax = 100.;
      // Adjust for pp, since we are retrieving multiplicity instead
      if (fForceBeamType == AliAnalysisTaskEmcal::kpp) {
        label = "Multiplicity";
        xmax = 200.;
      }
      break;

    case 1:
      label = "Jet p_{T}";
      nbins = 20;
      xmin = 0.;
      xmax = 200.;
      break;

    case 2:
      if(fDoWiderTrackBin) {
        label = "Track p_{T}";
        nbins = 40;
        xmin = 0.;
        xmax = 10.;
      } else {
        label = "Track p_{T}";
        nbins = 100;
        xmin = 0.;
        xmax = 10;
      }
      break;

    case 3:
      label = "#Delta#eta";
      nbins = 24;
      xmin = -1.2;
      xmax = 1.2;
      break;

    case 4:
      label = "#Delta#phi";
      nbins = 72;
      xmin = -0.5*pi;
      xmax = 1.5*pi;
      break;         

    case 5:
      label = "Leading Jet";
      nbins = 3;
      xmin = -0.5;
      xmax = 2.5;
      break;

    case 6:
      label = "Trigger track";
      nbins = 10;
      xmin = 0;
      xmax = 50;
      break;

    case 7:
      label = "deltaR";
      nbins = 10;
      xmin = 0.;
      xmax = 5.0;
      break;

    case 8:
      label = "Leading track";
      nbins = 20;
      xmin = 0;
      xmax = 50;
      break;

    case 9:
      label = "Event plane angle";
      nbins = 3;
      xmin = 0;
      xmax = TMath::Pi()/2.;
      break;
  }
}

/**
 * Clone tracks into a lighter object (AliBasicParticle) to store in the event pool. By using
 * lighter objects, it reduces the event pool size in memory. Adapted from CF event mixing
 * code PhiCorrelations.
 *
 * @return Array containing the lighter track objects.
 */
TObjArray* AliAnalysisTaskEmcalJetHCorrelations::CloneAndReduceTrackList(std::vector<unsigned int> & rejectedTrackIndices, const bool useRejectedList)
{
  // clones a track list by using AliBasicTrack which uses much less memory (used for event mixing)
  TObjArray* tracksClone = new TObjArray;
  tracksClone->SetOwner(kTRUE);

  // Loop over all tracks
  AliVParticle * particle = 0;
  AliBasicParticle * clone = 0;
  AliTrackContainer * tracks = GetTrackContainer("tracksForCorrelations");

  auto particlesIter = tracks->accepted_momentum();
  for (auto particleIter = particlesIter.begin(); particleIter != particlesIter.end(); particleIter++)
  {
    // Retrieve the particle
    particle = particleIter->second;

    // Artificial inefficiency
    bool rejectParticle = CheckArtificialTrackEfficiency(particleIter.current_index(), rejectedTrackIndices, useRejectedList);
    if (rejectParticle) {
      AliDebugStream(4) << "Track rejected in CloneAndReduceTrackList()\n";
      continue;
    }

    // Fill some QA information about the tracks
    Int_t trackPtBin = GetTrackPtBin(particle->Pt());
    if(trackPtBin > -1) fHistTrackEtaPhi[trackPtBin]->Fill(particle->Eta(),particle->Phi());

    // Create new particle
    clone = new AliBasicParticle(particle->Eta(), particle->Phi(), particle->Pt(), particle->Charge());
    // Set so that we can do comparisons using the IsEqual() function.
    clone->SetUniqueID(particle->GetUniqueID());

    tracksClone->Add(clone);
  }

  return tracksClone;
}

/**
 * Utility function to apply the efficiency correction. This function always uses fBeamType, which
 * is preferred when speed is desired (ie for analysis). The function below is used for external
 * testing of the efficiency correction.
 *
 * @param trackETA Eta of the track
 * @param trackPT pT of the track
 *
 * @return Track efficiency of the track (the entry should be weighted as 1/(return value))
 */
Double_t AliAnalysisTaskEmcalJetHCorrelations::EffCorrection(Double_t trackETA, Double_t trackPT) const {
  return EffCorrection(trackETA, trackPT, fBeamType);
}

/**
 * Determine the efficiency correction for a given track pT and eta. fDoEffCorrection determines
 * the mode of the correction:
 * - 0 disables the correction.
 * - 1 enables the correction. In Pb-Pb, this will automatically select the proper efficiency based
 *   on run list (Good vs Semi-good) and centrality.
 * - 2-9 Explicitly select an Pb-Pb efficiency correction function. It will not be automatically
 *   selected later!
 *
 * @param trackETA Eta of the track.
 * @param trackPT pT of the track.
 * @param beamType Type of collision system.
 *
 * @return Track efficiency of the track (the entry should be weighted as 1/(return value))
 */
Double_t AliAnalysisTaskEmcalJetHCorrelations::EffCorrection(Double_t trackETA, Double_t trackPT, AliAnalysisTaskEmcal::BeamType beamType) const
{
  // default (current) parameters
  // x-variable = track pt, y-variable = track eta
  Double_t x = trackPT;
  Double_t y = trackETA;
  double etaaxis = 0;
  double ptaxis = 0;

  // Efficiency paramaters
  double TRefficiency = -1;
  int effSwitch = 1;

  if (fSingleTrackEfficiencyCorrectionType == AliAnalysisTaskEmcalJetHCorrelations::kEffAutomaticConfiguration) {
    if (beamType == AliAnalysisTaskEmcal::kAA) {
      // Setup for Pb--Pb
      int runQuality = -1;
      // Semi-Good OROC C08 Runlists
      if (fCurrentRunNumber == 169975 || fCurrentRunNumber == 169981 || fCurrentRunNumber == 170038 || fCurrentRunNumber == 170040 || fCurrentRunNumber == 170083 || fCurrentRunNumber == 170084 || fCurrentRunNumber == 170085 || fCurrentRunNumber == 170088 || fCurrentRunNumber == 170089 || fCurrentRunNumber == 170091 || fCurrentRunNumber == 170152 || fCurrentRunNumber == 170155 || fCurrentRunNumber == 170159 || fCurrentRunNumber == 170163 || fCurrentRunNumber == 170193 || fCurrentRunNumber == 170195 || fCurrentRunNumber == 170203 || fCurrentRunNumber == 170204 || fCurrentRunNumber == 170228 || fCurrentRunNumber == 170230 || fCurrentRunNumber == 170268 || fCurrentRunNumber == 170269 || fCurrentRunNumber == 170270 || fCurrentRunNumber == 170306 || fCurrentRunNumber == 170308 || fCurrentRunNumber == 170309) runQuality = 0;

      // Good Runlists
      if (fCurrentRunNumber == 167902 || fCurrentRunNumber == 167903 || fCurrentRunNumber == 167915 || fCurrentRunNumber == 167920 || fCurrentRunNumber == 167987 || fCurrentRunNumber == 167988 || fCurrentRunNumber == 168066 || fCurrentRunNumber == 168068 || fCurrentRunNumber == 168069 || fCurrentRunNumber == 168076 || fCurrentRunNumber == 168104 || fCurrentRunNumber == 168107 || fCurrentRunNumber == 168108 || fCurrentRunNumber == 168115 || fCurrentRunNumber == 168212 || fCurrentRunNumber == 168310 || fCurrentRunNumber == 168311 || fCurrentRunNumber == 168322 || fCurrentRunNumber == 168325 || fCurrentRunNumber == 168341 || fCurrentRunNumber == 168342 || fCurrentRunNumber == 168361 || fCurrentRunNumber == 168362 || fCurrentRunNumber == 168458 || fCurrentRunNumber == 168460 || fCurrentRunNumber == 168461 || fCurrentRunNumber == 168464 || fCurrentRunNumber == 168467 || fCurrentRunNumber == 168511 || fCurrentRunNumber == 168512 || fCurrentRunNumber == 168777 || fCurrentRunNumber == 168826 || fCurrentRunNumber == 168984 || fCurrentRunNumber == 168988 || fCurrentRunNumber == 168992 || fCurrentRunNumber == 169035 || fCurrentRunNumber == 169091 || fCurrentRunNumber == 169094 || fCurrentRunNumber == 169138 || fCurrentRunNumber == 169143 || fCurrentRunNumber == 169144 || fCurrentRunNumber == 169145 || fCurrentRunNumber == 169148 || fCurrentRunNumber == 169156 || fCurrentRunNumber == 169160 || fCurrentRunNumber == 169167 || fCurrentRunNumber == 169238 || fCurrentRunNumber == 169411 || fCurrentRunNumber == 169415 || fCurrentRunNumber == 169417 || fCurrentRunNumber == 169835 || fCurrentRunNumber == 169837 || fCurrentRunNumber == 169838 || fCurrentRunNumber == 169846 || fCurrentRunNumber == 169855 || fCurrentRunNumber == 169858 || fCurrentRunNumber == 169859 || fCurrentRunNumber == 169923 || fCurrentRunNumber == 169956 || fCurrentRunNumber == 170027 || fCurrentRunNumber == 170036 || fCurrentRunNumber == 170081) runQuality = 1;

      // Determine which efficiency to use.
      // This is just a way to map all possible values of the cent bin and runQuality to a unique flag.
      // 4 is the number of cent bins, and we want to index the effSwitch starting at 2.
      if (runQuality != -1) {
        effSwitch = 2 + runQuality*4 + fCentBin;
      }
    }
    else if (beamType == AliAnalysisTaskEmcal::kpA) {
      // pA
      AliErrorStream() << "Single track efficiency for pA is not available.\n";
      return 0;
    }
    else if (beamType == AliAnalysisTaskEmcal::kpp) {
      // Setup for pp
      effSwitch = 10;
    }
    else {
      AliErrorStream() << "Beam type " << fBeamType << " is not defined\n";
      return 0;
    }
  }
  else if (fSingleTrackEfficiencyCorrectionType == AliAnalysisTaskEmcalJetHCorrelations::kEffPP) {
    // Manually set to pp (for example, during embedding)
    effSwitch = 10;
  }
  else if (fSingleTrackEfficiencyCorrectionType == AliAnalysisTaskEmcalJetHCorrelations::kEffDisable) {
    // Use the default value, which is to have a constant efficiency of 1
    effSwitch = 1;
  }
  else {
    AliErrorStream() << "Single track efficiency correction type " << fSingleTrackEfficiencyCorrectionType << " is not recongized!\n";
  }

  AliDebugStream(5) << "Using efficiency switch value of " << effSwitch << "\n";

  switch(effSwitch) {
    case 1 :
      // first switch value - TRefficiency not used so = 1
      // In this case, the run number isn't in any run list, so efficiency = 1
      TRefficiency = 1.0;
      break;

    case 2 :
      // Parameter values for Semi-GOOD TPC (LHC11h) runs (0-10%):
      ptaxis = (x<2.9)*(p0_10SG[0]*exp(-pow(p0_10SG[1]/x,p0_10SG[2])) + p0_10SG[3]*x) + (x>=2.9)*(p0_10SG[4] + p0_10SG[5]*x + p0_10SG[6]*x*x);
      etaaxis = (y<-0.07)*(p0_10SG[7]*exp(-pow(p0_10SG[8]/TMath::Abs(y+0.91),p0_10SG[9])) + p0_10SG[10]*y) + (y>=-0.07 && y<=0.4)*(p0_10SG[11] + p0_10SG[12]*y + p0_10SG[13]*y*y) + (y>0.4)*(p0_10SG[14]*exp(-pow(p0_10SG[15]/TMath::Abs(-y+0.91),p0_10SG[16])));
      TRefficiency = ptaxis*etaaxis;
      break;

    case 3 :
      // Parameter values for Semi-GOOD TPC (LHC11h) runs (10-30%):
      ptaxis = (x<2.9)*(p10_30SG[0]*exp(-pow(p10_30SG[1]/x,p10_30SG[2])) + p10_30SG[3]*x) + (x>=2.9)*(p10_30SG[4] + p10_30SG[5]*x + p10_30SG[6]*x*x);
      etaaxis = (y<-0.07)*(p10_30SG[7]*exp(-pow(p10_30SG[8]/TMath::Abs(y+0.91),p10_30SG[9])) + p10_30SG[10]*y) + (y>=-0.07 && y<=0.4)*(p10_30SG[11] + p10_30SG[12]*y + p10_30SG[13]*y*y) + (y>0.4)*(p10_30SG[14]*exp(-pow(p10_30SG[15]/TMath::Abs(-y+0.91),p10_30SG[16])));
      TRefficiency = ptaxis*etaaxis;
      break;

    case 4 :
      // Parameter values for Semi-GOOD TPC (LHC11h) runs (30-50%):
      ptaxis = (x<2.9)*(p30_50SG[0]*exp(-pow(p30_50SG[1]/x,p30_50SG[2])) + p30_50SG[3]*x) + (x>=2.9)*(p30_50SG[4] + p30_50SG[5]*x + p30_50SG[6]*x*x);
      etaaxis = (y<-0.07)*(p30_50SG[7]*exp(-pow(p30_50SG[8]/TMath::Abs(y+0.91),p30_50SG[9])) + p30_50SG[10]*y) + (y>=-0.07 && y<=0.4)*(p30_50SG[11] + p30_50SG[12]*y + p30_50SG[13]*y*y) + (y>0.4)*(p30_50SG[14]*exp(-pow(p30_50SG[15]/TMath::Abs(-y+0.91),p30_50SG[16])));
      TRefficiency = ptaxis*etaaxis;
      break;

    case 5 :
      // Parameter values for Semi-GOOD TPC (LHC11h) runs (50-90%):
      ptaxis = (x<2.9)*(p50_90SG[0]*exp(-pow(p50_90SG[1]/x,p50_90SG[2])) + p50_90SG[3]*x) + (x>=2.9)*(p50_90SG[4] + p50_90SG[5]*x + p50_90SG[6]*x*x);
      etaaxis = (y<-0.07)*(p50_90SG[7]*exp(-pow(p50_90SG[8]/TMath::Abs(y+0.91),p50_90SG[9])) + p50_90SG[10]*y) + (y>=-0.07 && y<=0.4)*(p50_90SG[11] + p50_90SG[12]*y + p50_90SG[13]*y*y) + (y>0.4)*(p50_90SG[14]*exp(-pow(p50_90SG[15]/TMath::Abs(-y+0.91),p50_90SG[16])));
      TRefficiency = ptaxis*etaaxis;
      break;

    case 6 :
      // Parameter values for GOOD TPC (LHC11h) runs (0-10%):
      ptaxis = (x<2.9)*(p0_10G[0]*exp(-pow(p0_10G[1]/x,p0_10G[2])) + p0_10G[3]*x) + (x>=2.9)*(p0_10G[4] + p0_10G[5]*x + p0_10G[6]*x*x);
      etaaxis = (y<0.0)*(p0_10G[7]*exp(-pow(p0_10G[8]/TMath::Abs(y+0.91),p0_10G[9])) + p0_10G[10]*y) + (y>=0.0 && y<=0.4)*(p0_10G[11] + p0_10G[12]*y + p0_10G[13]*y*y) + (y>0.4)*(p0_10G[14]*exp(-pow(p0_10G[15]/TMath::Abs(-y+0.91),p0_10G[16])));
      TRefficiency = ptaxis*etaaxis;
      break;

    case 7 :
      // Parameter values for GOOD TPC (LHC11h) runs (10-30%):
      ptaxis = (x<2.9)*(p10_30G[0]*exp(-pow(p10_30G[1]/x,p10_30G[2])) + p10_30G[3]*x) + (x>=2.9)*(p10_30G[4] + p10_30G[5]*x + p10_30G[6]*x*x);
      etaaxis = (y<0.0)*(p10_30G[7]*exp(-pow(p10_30G[8]/TMath::Abs(y+0.91),p10_30G[9])) + p10_30G[10]*y) + (y>=0.0 && y<=0.4)*(p10_30G[11] + p10_30G[12]*y + p10_30G[13]*y*y) + (y>0.4)*(p10_30G[14]*exp(-pow(p10_30G[15]/TMath::Abs(-y+0.91),p10_30G[16])));
      TRefficiency = ptaxis*etaaxis;
      break;

    case 8 :
      // Parameter values for GOOD TPC (LHC11h) runs (30-50%):
      ptaxis = (x<2.9)*(p30_50G[0]*exp(-pow(p30_50G[1]/x,p30_50G[2])) + p30_50G[3]*x) + (x>=2.9)*(p30_50G[4] + p30_50G[5]*x + p30_50G[6]*x*x);
      etaaxis = (y<0.0)*(p30_50G[7]*exp(-pow(p30_50G[8]/TMath::Abs(y+0.91),p30_50G[9])) + p30_50G[10]*y) + (y>=0.0 && y<=0.4)*(p30_50G[11] + p30_50G[12]*y + p30_50G[13]*y*y) + (y>0.4)*(p30_50G[14]*exp(-pow(p30_50G[15]/TMath::Abs(-y+0.91),p30_50G[16])));
      TRefficiency = ptaxis*etaaxis;
      break;

    case 9 :
      // Parameter values for GOOD TPC (LHC11h) runs (50-90%):
      ptaxis = (x<2.9)*(p50_90G[0]*exp(-pow(p50_90G[1]/x,p50_90G[2])) + p50_90G[3]*x) + (x>=2.9)*(p50_90G[4] + p50_90G[5]*x + p50_90G[6]*x*x);
      etaaxis = (y<0.0)*(p50_90G[7]*exp(-pow(p50_90G[8]/TMath::Abs(y+0.91),p50_90G[9])) + p50_90G[10]*y) + (y>=0.0 && y<=0.4)*(p50_90G[11] + p50_90G[12]*y + p50_90G[13]*y*y) + (y>0.4)*(p50_90G[14]*exp(-pow(p50_90G[15]/TMath::Abs(-y+0.91),p50_90G[16])));
      TRefficiency = ptaxis*etaaxis;
      break;

    case 10 :
      {
        // Track efficiency for pp
        // Calculated using LHC12f1a. See analysis note for more details!
        // If the trackPt > 6 GeV, then all we need is this coefficient
        Double_t coefficient = 0.898052;                                                // p6
        if (trackPT < 6) {
          coefficient =  (1 + -0.442232 * trackPT                                     // p0
                   +  0.501831 * std::pow(trackPT, 2)                        // p1
                   + -0.252024 * std::pow(trackPT, 3)                        // p2
                   +  0.062964 * std::pow(trackPT, 4)                        // p3
                   + -0.007681 * std::pow(trackPT, 5)                        // p4
                   +  0.000365 * std::pow(trackPT, 6));                      // p5
        }

        // Calculate track eff
        TRefficiency = coefficient * (1 +  0.402825 * std::abs(trackETA)                // p7
                        + -2.213152 * std::pow(trackETA, 2)             // p8
                        +  4.311098 * std::abs(std::pow(trackETA, 3))   // p9
                        + -2.778200 * std::pow(trackETA, 4));           // p10
        break;
      }

    default :
      // no Efficiency Switch option selected.. therefore don't correct, and set eff = 1
      // ie. The efficiency correction is disabled.
      AliErrorStream() << "No single track efficiency setting selected! Please select one.\n";
      TRefficiency = 0.;
  }

  return TRefficiency;
}

/**
 * Utility function to fill a histogram with a given weight. If requested, it will also apply the JES correction derived weight
 * to the hist.
 *
 * @param[in] hist Histogram to be filled.
 * @param[in] fillValue The value to be filled into the hist.
 * @param[in] weight The weight to be applied when filling the hist.
 * @param[in] noCorrection True if the JES correction _should not_ be applied.
 */
void AliAnalysisTaskEmcalJetHCorrelations::FillHist(TH1 * hist, Double_t fillValue, Double_t weight, Bool_t noCorrection)
{
  if (fJESCorrectionHist == 0 || noCorrection == kTRUE)
  {
    AliDebugStream(3) << GetName() << ":" << hist->GetName() << ": " << std::boolalpha << "Using normal weights: JESHist: " << (fJESCorrectionHist ? fJESCorrectionHist->GetName() : "Null") << ", noCorrection: " << noCorrection << std::endl;
    hist->Fill(fillValue, weight);
  }
  else
  {
    // Determine where to get the values in the correction hist
    Int_t xBin = fJESCorrectionHist->GetXaxis()->FindBin(fillValue);

    std::vector <Double_t> yBinsContent;
    AliDebug(3, TString::Format("%s: Attempt to access weights from JES correction hist %s with jet pt %f!", GetName(), hist->GetName(), fillValue));
    AccessSetOfYBinValues(fJESCorrectionHist, xBin, yBinsContent);
    AliDebug(3, TString::Format("weights size: %zd", yBinsContent.size()));

    // Loop over all possible bins to contribute.
    // If content is 0 then calling Fill won't make a difference
    for (Int_t index = 1; index <= fJESCorrectionHist->GetYaxis()->GetNbins(); index++)
    {
      // Don't bother trying to fill in the weight is 0
      if (yBinsContent.at(index-1) > 0) {
        // Determine the value to fill based on the center of the bins.
        // This in principle allows the binning between the correction and hist to be different
        Double_t fillLocation = fJESCorrectionHist->GetYaxis()->GetBinCenter(index); 
        AliDebug(4, TString::Format("fillLocation: %f, weight: %f", fillLocation, yBinsContent.at(index-1)));
        // minus 1 since loop starts at 1
        hist->Fill(fillLocation, weight*yBinsContent.at(index-1));
      }
    }

    //TEMP
    //hist->Draw();
    //END TEMP
  }
}

/**
 * Utility function to fill a histogram with a given weight. If requested, it will also apply the JES correction derived weight
 * to the hist. It assumes that the corrected jet pt value is located in the second element of the fillValue (as defined in
 * GetDimParams()). If that is changed or the first entry is not included, then this function must be updated!
 *
 * @param[in] hist Histogram to be filled.
 * @param[in] fillValue Array of values to be filled into the hist.
 * @param[in] weight The weight to be applied when filling the hist.
 * @param[in] noCorrection True if the JES correction _should not_ be applied.
 */
void AliAnalysisTaskEmcalJetHCorrelations::FillHist(THnSparse * hist, Double_t *fillValue, Double_t weight, Bool_t noCorrection)
{
  if (fJESCorrectionHist == 0 || noCorrection == kTRUE)
  {
    AliDebugStream(3) << GetName() << ":" << hist->GetName() << ": " << std::boolalpha << "Using normal weights: JESHist: " << (fJESCorrectionHist ? fJESCorrectionHist->GetName() : "Null") << ", noCorrection: " << noCorrection << std::endl;
    hist->Fill(fillValue, weight);
  }
  else
  {
    // Jet pt is always located in the second position
    Double_t jetPt = fillValue[1];

    // Determine where to get the values in the correction hist
    Int_t xBin = fJESCorrectionHist->GetXaxis()->FindBin(jetPt);

    std::vector <Double_t> yBinsContent;
    AliDebug(3, TString::Format("%s: Attempt to access weights from JES correction hist %s with jet pt %f!", GetName(), hist->GetName(), jetPt));
    AccessSetOfYBinValues(fJESCorrectionHist, xBin, yBinsContent);
    AliDebug(3, TString::Format("weights size: %zd", yBinsContent.size()));

    // Loop over all possible bins to contribute.
    // If content is 0 then calling Fill won't make a difference
    for (Int_t index = 1; index <= fJESCorrectionHist->GetYaxis()->GetNbins(); index++)
    {
      // Don't bother trying to fill in the weight is 0
      if (yBinsContent.at(index-1) > 0) {
        // Determine the value to fill based on the center of the bins.
        // This in principle allows the binning between the correction and hist to be different
        fillValue[1] = fJESCorrectionHist->GetYaxis()->GetBinCenter(index); 
        AliDebug(4,TString::Format("fillValue[1]: %f, weight: %f", fillValue[1], yBinsContent.at(index-1)));
        // minus 1 since loop starts at 1
        hist->Fill(fillValue, weight*yBinsContent.at(index-1));
      }
    }
  }
}

/**
 * Access the array of y bin values for a given x bin. If the scaleFactor is greater than 0, then the values in yBinsContent will be set
 * in the given histogram with the values scaled by the scale factor. This scaling can be used to normalize the histogram.
 *
 * @param[in] hist Histogram from which the bins should be extracted.
 * @param[in] xBin X bin where the y bins should be extracted.
 * @param[in,out] yBinsContent Array containing the y bins contents for the given x bin.
 * @param[in] scaleFactor Scale factor to be applied to the
 */
void AliAnalysisTaskEmcalJetHCorrelations::AccessSetOfYBinValues(TH2D * hist, Int_t xBin, std::vector <Double_t> & yBinsContent, Double_t scaleFactor)
{
  for (Int_t index = 1; index <= hist->GetYaxis()->GetNbins(); index++)
  {
    //yBinsContent[index-1] = hist->GetBinContent(hist->GetBin(xBin,index));
    yBinsContent.push_back(hist->GetBinContent(hist->GetBin(xBin,index)));

    if (scaleFactor >= 0)
    {
      // -1 since index starts at 1
      hist->SetBinContent(hist->GetBin(xBin,index), yBinsContent.at(index-1)/scaleFactor);
    }
  }
}

/**
 * Attempt to retrieve and initialize the jet energy scale correction histogram with a given name.
 * If successfully initialized, "_JESCorr" is added to the task name before the last underscore
 * (which is usually the suffix).
 *
 * @param filename Name of file which contais the JES correction histogram.
 * @param histName Name of the JES correction histogram in the file.
 * @param trackBias The track bias used to create the JES correction histogram. Usually should match
 *     the track bias set in the task. The string ("_track%.2f", trackBias) will be appended to the hist name
 *     if this value is not set to be disabled.
 * @param clusterBias The cluster bias used to create the JES correction histogram. Usually should match
 *     the track bias set in the task. The string ("_clus%.2f", clusterBias) will be appended to the hist name
 *     if this value is not set to be disabled.
 *
 * @return kTRUE if the histogram is successfully retrieve and initialized.
 */
Bool_t AliAnalysisTaskEmcalJetHCorrelations::RetrieveAndInitializeJESCorrectionHist(TString filename, TString histName, Double_t trackBias, Double_t clusterBias)
{
  // Initialize grid connection if necessary
  if (filename.Contains("alien://") && !gGrid) {
    TGrid::Connect("alien://");
  }

  // Setup hist name if a track or cluster bias was defined.
  // NOTE: This can always be disabled by setting kDisableBias.
  //       We arbitrarily add 0.1 to test since the values are doubles and cannot be
  //       tested directly for equality. If we are still less than disable bins, then
  //       it has been set and we should format it.
  // NOTE: To ensure we can disable, we don't just take the member values!
  // NOTE: The histBaseName will be attempted if the formatted name cannot be found.
  TString histBaseName = histName;
  if (trackBias + 0.1 < AliAnalysisTaskEmcalJetHCorrelations::kDisableBias) {
    histName = TString::Format("%s_Track%.2f", histName.Data(), trackBias);
  }
  if (clusterBias + 0.1 < AliAnalysisTaskEmcalJetHCorrelations::kDisableBias) {
    histName = TString::Format("%s_Clus%.2f", histName.Data(), clusterBias);
  }

  // Open file containing the correction
  TFile * jesCorrectionFile = TFile::Open(filename);
  if (!jesCorrectionFile || jesCorrectionFile->IsZombie()) {
    AliError(TString::Format("%s: Could not open JES correction file %s", GetName(), filename.Data()));
    return kFALSE;
  }

  // Retrieve the histogram containing the correction and safely add it to the task.
  TH2D * JESCorrectionHist = dynamic_cast<TH2D*>(jesCorrectionFile->Get(histName.Data()));
  if (JESCorrectionHist) {
    AliInfo(TString::Format("%s: JES correction hist name \"%s\" loaded from file %s.", GetName(), histName.Data(), filename.Data()));
  }
  else {
    AliError(TString::Format("%s: JES correction hist name \"%s\" not found in file %s.", GetName(), histName.Data(), filename.Data()));

    // Attempt the base name instead of the formatted hist name
    JESCorrectionHist = dynamic_cast<TH2D*>(jesCorrectionFile->Get(histBaseName.Data()));
    if (JESCorrectionHist) {
      AliInfo(TString::Format("%s: JES correction hist name \"%s\" loaded from file %s.", GetName(), histBaseName.Data(), filename.Data()));
      histName = histBaseName;
    }
    else
    {
      AliError(TString::Format("%s: JES correction with base hist name %s not found in file %s.", GetName(), histBaseName.Data(), filename.Data()));
      return kFALSE;
    }
  }

  // Clone to ensure that the hist is available
  TH2D * tempHist = static_cast<TH2D *>(JESCorrectionHist->Clone());
  tempHist->SetDirectory(0);
  SetJESCorrectionHist(tempHist);

  // Close file
  jesCorrectionFile->Close();

  // Append to task name for clarity
  // Unfortunately, this doesn't change the name of the output list (it would need to be
  // changed in the AnalysisManager output container), so the suffix is still important
  // if this correction is manually configured!
  TString tempName = GetName();
  TString tag = "_JESCorr";
  // Append the tag if it isn't already included
  if (tempName.Index(tag) == -1) {
    // Insert before the suffix
    Ssiz_t suffixLocation = tempName.Last('_');
    tempName.Insert(suffixLocation, tag.Data());

    // Set the new name
    AliDebug(3, TString::Format("%s: Setting task name to %s", GetName(), tempName.Data()));
    SetName(tempName.Data());
  }

  // Successful
  return kTRUE;
}

/**
 * AddTask for the jet-hadron task. We benefit for actually having compiled code, as opposed to
 * struggling with CINT.
 */
AliAnalysisTaskEmcalJetHCorrelations * AliAnalysisTaskEmcalJetHCorrelations::AddTaskEmcalJetHCorrelations(
   const char *nTracks,
   const char *nCaloClusters,
   // Jet options
   const Double_t trackBias,
   const Double_t clusterBias,
   // Mixed event options
   const Int_t nTracksMixedEvent,  // Additionally acts as a switch for enabling mixed events
   const Int_t minNTracksMixedEvent,
   const Int_t minNEventsMixedEvent,
   const UInt_t nCentBinsMixedEvent,
   // Triggers
   UInt_t trigEvent,
   UInt_t mixEvent,
   // Options
   const Bool_t lessSparseAxes,
   const Bool_t widerTrackBin,
   // Corrections
   const AliAnalysisTaskEmcalJetHCorrelations::ESingleTrackEfficiency_t singleTrackEfficiency,
   const Bool_t JESCorrection,
   const char * JESCorrectionFilename,
   const char * JESCorrectionHistName,
   const char *suffix
   )
{
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    AliErrorClass("No analysis manager to connect to.");
    return nullptr;
  }

  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  // Determine cluster and track names
  TString trackName(nTracks);
  TString clusName(nCaloClusters);

  if (trackName == "usedefault") {
    trackName = AliEmcalContainerUtils::DetermineUseDefaultName(AliEmcalContainerUtils::kTrack);
  }

  if (clusName == "usedefault") {
    clusName = AliEmcalContainerUtils::DetermineUseDefaultName(AliEmcalContainerUtils::kCluster);
  }

  TString name("AliAnalysisTaskJetH");
  if (!trackName.IsNull()) {
    name += TString::Format("_%s", trackName.Data());
  }
  if (!clusName.IsNull()) {
    name += TString::Format("_%s", clusName.Data());
  }
  if (strcmp(suffix, "") != 0) {
    name += TString::Format("_%s", suffix);
  }

  AliAnalysisTaskEmcalJetHCorrelations *correlationTask = new AliAnalysisTaskEmcalJetHCorrelations(name);
  // Set jet bias
  correlationTask->SetTrackBias(trackBias);
  correlationTask->SetClusterBias(clusterBias);
  // Mixed events
  correlationTask->SetEventMixing(static_cast<Bool_t>(nTracksMixedEvent));
  correlationTask->SetNumberOfMixingTracks(nTracksMixedEvent);
  correlationTask->SetMinNTracksForMixedEvents(minNTracksMixedEvent);
  correlationTask->SetMinNEventsForMixedEvents(minNEventsMixedEvent);
  correlationTask->SetNCentBinsMixedEvent(nCentBinsMixedEvent);
  // Triggers
  correlationTask->SetTriggerType(trigEvent);
  correlationTask->SetMixedEventTriggerType(mixEvent);
  // Options
  correlationTask->SetNCentBins(5);
  correlationTask->SetVzRange(-10,10);
  correlationTask->SetDoLessSparseAxes(lessSparseAxes);
  correlationTask->SetDoWiderTrackBin(widerTrackBin);
  // Corrections
  correlationTask->SetSingleTrackEfficiencyType(singleTrackEfficiency);
  if (JESCorrection == kTRUE)
  {
    Bool_t result = correlationTask->RetrieveAndInitializeJESCorrectionHist(JESCorrectionFilename, JESCorrectionHistName, correlationTask->GetTrackBias(), correlationTask->GetClusterBias());
    if (!result) {
      AliErrorClass("Failed to successfully retrieve and initialize the JES correction! Task initialization continuing without JES correction (can be set manually later).");
    }
  }

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------

  mgr->AddTask(correlationTask);

  // Create containers for input/output
  mgr->ConnectInput (correlationTask, 0, mgr->GetCommonInputContainer() );
  AliAnalysisDataContainer * cojeth = mgr->CreateContainer(correlationTask->GetName(),
                TList::Class(),
                AliAnalysisManager::kOutputContainer,
							    Form("%s", AliAnalysisManager::GetCommonFileName()));
  mgr->ConnectOutput(correlationTask, 1, cojeth);

  return correlationTask;
}

/**
 * Call after the AddTask to setup the standard Jet-h configuration.
 *
 * @return true if sucessfully configured for the standard configuration
 */
bool AliAnalysisTaskEmcalJetHCorrelations::ConfigureForStandardAnalysis(std::string trackName,
    std::string clusName,
    const double jetConstituentPtCut,
    const double trackEta,
    const double jetRadius)
{
  bool returnValue = false;
  AliInfoStream() << "Configuring Jet-H Correlations task for a standard analysis.\n";

  // Add Containers
  // Clusters
  if (clusName == "usedefault") {
    clusName = AliEmcalContainerUtils::DetermineUseDefaultName(AliEmcalContainerUtils::kCluster);
  }
  // For jet finding
  AliClusterContainer * clustersForJets = new AliClusterContainer(clusName.c_str());
  clustersForJets->SetName("clustersForJets");
  clustersForJets->SetMinE(jetConstituentPtCut);

  // Tracks
  // For jet finding
  if (trackName == "usedefault") {
    trackName = AliEmcalContainerUtils::DetermineUseDefaultName(AliEmcalContainerUtils::kTrack);
  }
  AliParticleContainer * particlesForJets = CreateParticleOrTrackContainer(trackName.c_str());
  particlesForJets->SetName("particlesForJets");
  particlesForJets->SetMinPt(jetConstituentPtCut);
  particlesForJets->SetEtaLimits(-1.0*trackEta, trackEta);
  // Don't need to adopt the container - we'll just use it to find the right jet collection
  // For correlations
  AliParticleContainer * particlesForCorrelations = CreateParticleOrTrackContainer(trackName.c_str());
  if (particlesForCorrelations)
  {
    particlesForCorrelations->SetName("tracksForCorrelations");
    particlesForCorrelations->SetMinPt(0.15);
    particlesForCorrelations->SetEtaLimits(-1.0*trackEta, trackEta);
    // Adopt the container
    this->AdoptParticleContainer(particlesForCorrelations);
  }
  else {
    AliWarningStream() << "No particle container was successfully created!\n";
  }

  // Jets
  AliJetContainer * jetContainer = this->AddJetContainer(AliJetContainer::kFullJet,
                              AliJetContainer::antikt_algorithm,
                              AliJetContainer::pt_scheme,
                              jetRadius,
                              AliEmcalJet::kEMCALfid,
                              particlesForJets,
                              clustersForJets);
  // 0.6 * jet area
  jetContainer->SetJetAreaCut(jetRadius * jetRadius * TMath::Pi() * 0.6);
  jetContainer->SetMaxTrackPt(100);
  jetContainer->SetJetPtCut(0.1);

  // Successfully configured
  returnValue = true;

  return returnValue;
}

/**
 * Call after the AddTask to setup the embedded Jet-h configuration.
 *
 * @return true if successfully configured for embedding
 */
bool AliAnalysisTaskEmcalJetHCorrelations::ConfigureForEmbeddingAnalysis(std::string trackName,
    std::string clusName,
    const double jetConstituentPtCut,
    const double trackEta,
    const double jetRadius,
    const std::string & jetTag,
    const std::string & correlationsTracksCutsPeriod)
{
  bool returnValue = false;
  AliInfoStream() << "Configuring Jet-H Correlations task for an embedding analysis.\n";

  // Set the task to know it that is embedded
  this->SetIsEmbedded(true);

  // Add Containers
  // Clusters
  if (clusName == "usedefault") {
    clusName = AliEmcalContainerUtils::DetermineUseDefaultName(AliEmcalContainerUtils::kCluster);
  }
  // For jet finding
  AliClusterContainer * clustersForJets = new AliClusterContainer(clusName.c_str());
  clustersForJets->SetName("clustersForJets");
  clustersForJets->SetMinE(jetConstituentPtCut);
  // We need the combined clusters, which should be available in the internal event.
  // However, we don't need to adopt the container - we'll just use it to find the right jet collection
  // For correlations
  /*AliClusterContainer * clustersforCorrelations = new AliClusterContainer("usedefault");
  clustersForCorrelations->SetName("clustersForCorrelations");
  clustersForCorrelations->SetMinE(0.30);
  clustersForCorrelations->SetIsEmbedding(true);
  this->AdoptClusterContainer(clustersForCorrelations);*/

  // Tracks
  // For jet finding
  if (trackName == "usedefault") {
    trackName = AliEmcalContainerUtils::DetermineUseDefaultName(AliEmcalContainerUtils::kTrack);
  }
  AliParticleContainer * particlesForJets = CreateParticleOrTrackContainer(trackName.c_str());
  particlesForJets->SetName("particlesForJets");
  particlesForJets->SetMinPt(jetConstituentPtCut);
  particlesForJets->SetEtaLimits(-1.0*trackEta, trackEta);
  // Don't need to adopt the container - we'll just use it to find the right jet collection
  // For correlations
  AliParticleContainer * particlesForCorrelations = CreateParticleOrTrackContainer(trackName.c_str());
  // Ensure that we don't operate on a null pointer
  if (particlesForCorrelations)
  {
    particlesForCorrelations->SetName("tracksForCorrelations");
    particlesForCorrelations->SetMinPt(0.15);
    particlesForCorrelations->SetEtaLimits(-1.0*trackEta, trackEta);
    particlesForCorrelations->SetIsEmbedding(true);
    AliTrackContainer * trackCont = dynamic_cast<AliTrackContainer *>(particlesForCorrelations);
    if (trackCont) {
      // This option only exists for track containers
      trackCont->SetTrackCutsPeriod(correlationsTracksCutsPeriod.c_str());
    }
    // Adopt the container
    this->AdoptParticleContainer(particlesForCorrelations);
  }
  else {
    AliWarningStream() << "No particle container was successfully created!\n";
  }

  // Jets
  // The tag "hybridLevelJets" is defined in the jet finder
  AliJetContainer * jetContainer = this->AddJetContainer(AliJetContainer::kFullJet,
                              AliJetContainer::antikt_algorithm,
                              AliJetContainer::pt_scheme,
                              jetRadius,
                              AliEmcalJet::kEMCALfid,
                              particlesForJets,
                              clustersForJets,
                              jetTag);
  // 0.6 * jet area
  jetContainer->SetJetAreaCut(jetRadius * jetRadius * TMath::Pi() * 0.6);
  jetContainer->SetMaxTrackPt(100);
  jetContainer->SetJetPtCut(0.1);

  // Successfully configured
  returnValue = true;

  return returnValue;
}

/**
 * Utility function to create a particle or track container given the collection name of the desired container.
 *
 * @param[in] collectionName Name of the particle or track collection name.
 *
 * @return A newly created particle or track container.
 */
AliParticleContainer * AliAnalysisTaskEmcalJetHCorrelations::CreateParticleOrTrackContainer(const std::string & collectionName) const
{
  AliParticleContainer * partCont = 0;
  if (collectionName == AliEmcalContainerUtils::DetermineUseDefaultName(AliEmcalContainerUtils::kTrack)) {
    AliTrackContainer * trackCont = new AliTrackContainer(collectionName.c_str());
    partCont = trackCont;
  }
  else if (collectionName != "") {
    partCont = new AliParticleContainer(collectionName.c_str());
  }

  return partCont;
}

} /* namespace EMCALJetTasks */
} /* namespace PWGJE */
