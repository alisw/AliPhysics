//////////
// Measure Jet-hadron correlations
// Does event Mixing using AliEventPoolManager
/////////

#include "AliAnalysisTaskEmcalJetHdEdxCorrelations.h"

#include <bitset>
#include <vector>
#include <math.h>

#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THnSparse.h>
#include <TVector3.h>
#include <TFile.h>
#include <TGrid.h>
#include <TList.h>
#include <TMap.h>

#include <AliAnalysisManager.h>
#include <AliInputEventHandler.h>
#include <AliEventPoolManager.h>
#include <AliLog.h>
#include <AliVAODHeader.h>
#include <AliVTrack.h>

#include "AliEmcalJet.h"
#include "AliTLorentzVector.h"
#include "AliAODTrack.h"
#include "AliBasicParticle.h"
#include "AliEmcalContainerUtils.h"
#include "AliClusterContainer.h"
#include "AliTrackContainer.h"
#include "AliJetContainer.h"
#include "AliAnalysisTaskEmcalEmbeddingHelper.h"
#include "AliAnalysisTaskEmcalJetHUtils.h"
#include "AliPIDResponse.h"

ClassImp(PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalJetHdEdxCorrelations);
using namespace std;
namespace PWGJE
{
  namespace EMCALJetTasks
  {

    /**
     * Default constructor.
     */
    AliAnalysisTaskEmcalJetHdEdxCorrelations::AliAnalysisTaskEmcalJetHdEdxCorrelations() : AliAnalysisTaskEmcalJet("AliAnalysisTaskEmcalJetHdEdxCorrelations", kFALSE),
											   fqnVectorReader(0),
                                                                                           fYAMLConfig(),
                                                                                           fConfigurationInitialized(false),
                                                                                           fTrackBias(5),
                                                                                           fClusterBias(5),
                                                                                           fDoEventMixing(kFALSE),
                                                                                           fNMixingTracks(50000), fMinNTracksMixedEvents(5000), fMinNEventsMixedEvents(5), fNCentBinsMixedEvent(10),
                                                                                           fPoolMgr(nullptr),
                                                                                           fTriggerType(AliVEvent::kEMCEJE), fMixingEventType(AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral),
                                                                                           fDisableFastPartition(kFALSE),
                                                                                           fRandom(0),
                                                                                           fEfficiencyPeriodIdentifier(AliAnalysisTaskEmcalJetHUtils::kLHC18qr),
                                                                                           fArtificialTrackInefficiency(1.0),
                                                                                           fNoMixedEventJESCorrection(kFALSE),
                                                                                           fJESCorrectionHist(nullptr),
                                                                                           fDoLessSparseAxes(kFALSE), fDoWiderTrackBin(kFALSE),
                                                                                           fRequireMatchedJetWhenEmbedding(kTRUE),
                                                                                           fMinSharedMomentumFraction(0.),
                                                                                           fRequireMatchedPartLevelJet(false),
                                                                                           fMaxMatchedJetDistance(-1),
                                                                                           fHistManager(),
                                                                                           fHistJetHTrackPt(nullptr),
                                                                                           fHistEPAngle(nullptr),
                                                                                           fHistJetEtaPhi(nullptr),
                                                                                           fHistJetHEtaPhi(nullptr),
                                                                                           fhnMixedEvents(nullptr),
                                                                                           fhnJH(nullptr),
                                                                                           fhnTrigger(nullptr),
											   hShiftCoeff(nullptr), hShiftEvents(nullptr),
											   previousrunnumber(-1),
											   fMakeShiftCorrectionCalibHists(kFALSE),
											   fDoShiftCorrection(kFALSE),
											   fEPshiftHistoList(nullptr)
    {
      // Default Constructor
      InitializeArraysToZero();
    }

    /**
     * Standard constructor
     */
    AliAnalysisTaskEmcalJetHdEdxCorrelations::AliAnalysisTaskEmcalJetHdEdxCorrelations(const char *name) : AliAnalysisTaskEmcalJet(name, kTRUE),
													   fqnVectorReader(0),
                                                                                                           fYAMLConfig(),
                                                                                                           fConfigurationInitialized(false),
                                                                                                           fTrackBias(5),
                                                                                                           fClusterBias(5),
                                                                                                           fDoEventMixing(kFALSE),
                                                                                                           fNMixingTracks(50000), fMinNTracksMixedEvents(5000), fMinNEventsMixedEvents(5), fNCentBinsMixedEvent(10),
                                                                                                           fPoolMgr(nullptr),
                                                                                                           fTriggerType(AliVEvent::kEMCEJE), fMixingEventType(AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral),
                                                                                                           fDisableFastPartition(kFALSE),
                                                                                                           fRandom(0),
                                                                                                           fEfficiencyPeriodIdentifier(AliAnalysisTaskEmcalJetHUtils::kLHC18qr),
                                                                                                           fArtificialTrackInefficiency(1.0),
                                                                                                           fNoMixedEventJESCorrection(kFALSE),
                                                                                                           fJESCorrectionHist(nullptr),
                                                                                                           fDoLessSparseAxes(kFALSE), fDoWiderTrackBin(kFALSE),
                                                                                                           fRequireMatchedJetWhenEmbedding(kTRUE),
                                                                                                           fMinSharedMomentumFraction(0.),
                                                                                                           fRequireMatchedPartLevelJet(false),
                                                                                                           fMaxMatchedJetDistance(-1),
                                                                                                           fHistManager(name),
                                                                                                           fHistJetHTrackPt(nullptr),
                                                                                                           fHistEPAngle(nullptr),
                                                                                                           fHistJetEtaPhi(nullptr),
                                                                                                           fHistJetHEtaPhi(nullptr),
                                                                                                           fhnMixedEvents(nullptr),
                                                                                                           fhnJH(nullptr),
                                                                                                           fhnTrigger(nullptr),
													   hShiftCoeff(nullptr), hShiftEvents(nullptr),
													   previousrunnumber(-1),
													   fMakeShiftCorrectionCalibHists(kFALSE),
													   fDoShiftCorrection(kFALSE),
													   fEPshiftHistoList(nullptr)
    {
      // Constructor
      InitializeArraysToZero();
      // Ensure that additional general histograms are created
      SetMakeGeneralHistograms(kTRUE);
    }

    /**
     * Initialize all member arrays to nullptr. Used as a convenience function in the constructors.
     */
    void AliAnalysisTaskEmcalJetHdEdxCorrelations::InitializeArraysToZero()
    {
      for (Int_t trackPtBin = 0; trackPtBin < kMaxTrackPtBins; trackPtBin++)
      {
        fHistTrackEtaPhi[trackPtBin] = nullptr;
      }
      for (Int_t centralityBin = 0; centralityBin < kMaxCentralityBins; ++centralityBin)
      {
        fHistJetPt[centralityBin] = nullptr;
        fHistJetPtBias[centralityBin] = nullptr;
      }
    }

    /**
     * Initialize task.
     */
    bool AliAnalysisTaskEmcalJetHdEdxCorrelations::Initialize()
    {
      fConfigurationInitialized = false;

      // Ensure that we have at least one configuration in the YAML config.
      if (fYAMLConfig.DoesConfigurationExist(0) == false)
      {
        std::cerr << "No configurations exist in the YAML configuration. Returning immediately.\n";
        // No configurations exist. Return immediately.
        return fConfigurationInitialized;
      }

      // Always initialize for streaming purposes
      fYAMLConfig.Initialize();

      // Setup task based on the properties defined in the YAML config
      AliDebugStream(2) << "Configuring task from the YAML configuration.\n";
      RetrieveAndSetTaskPropertiesFromYAMLConfig();
      AliDebugStream(2) << "Finished configuring via the YAML configuration.\n";

      // Print the results of the initialization
      // Print outside of the ALICE Log system to ensure that it is always available!

      fConfigurationInitialized = true;
      return fConfigurationInitialized;
    }

    /**
     * Perform task configuration via YAML.
     */
    void AliAnalysisTaskEmcalJetHdEdxCorrelations::RetrieveAndSetTaskPropertiesFromYAMLConfig()
    {
      // Base class options
      // Task physics (trigger) selection.
      std::string baseName = "eventCuts";
      std::vector<std::string> physicsSelection;
      bool res = fYAMLConfig.GetProperty(std::vector<std::string>({"eventCuts", "physicsSelection"}), physicsSelection, false);
      if (res)
      {
        fOfflineTriggerMask = AliEmcalContainerUtils::DeterminePhysicsSelectionFromYAML(physicsSelection);
      }

      // Event cuts
      // If event cuts are enabled (which they exceptionally are by default), then we want to configure them here.
      // If the event cuts are explicitly disabled, then we invert that value to enable the AliAnylsisTaskEmcal
      // builtin event selection.
      bool tempBool;
      fYAMLConfig.GetProperty({baseName, "enabled"}, tempBool, false);
      fUseBuiltinEventSelection = !tempBool;
      if (fUseBuiltinEventSelection == false)
      {
        // Need to include the namespace so that AliDebug will work properly...
        std::string taskName = "PWGJE::EMCALJetTasks::";
        taskName += GetName();
        AliAnalysisTaskEmcalJetHUtils::ConfigureEventCuts(fAliEventCuts, fYAMLConfig, fOfflineTriggerMask, baseName, taskName);
      }

      // General task options
      baseName = "general";
      fYAMLConfig.GetProperty({baseName, "nCentBins"}, fNcentBins, false);

      // Efficiency
      std::string tempStr = "";
      baseName = "efficiency";
      res = fYAMLConfig.GetProperty({baseName, "periodIdentifier"}, tempStr, false);
      if (res)
      {
        fEfficiencyPeriodIdentifier = AliAnalysisTaskEmcalJetHUtils::fgkEfficiencyPeriodIdentifier.at(tempStr);
      }
    }

    void AliAnalysisTaskEmcalJetHdEdxCorrelations::CreateHistograms(){

      fHistJetHTrackPt = new TH1F("fHistJetHTrackPt", "P_{T} distribution", 1000, 0.0, 100.0);
      fHistEPAngle = new TH1F("fHistEPAngle", "#Psi_{2} distribution for 2018 calib", 50, 0, TMath::Pi());
      fHistJetEtaPhi = new TH2F("fHistJetEtaPhi", "Jet eta-phi", 900, -1.8, 1.8, 720, -3.2, 3.2);
      fHistJetHEtaPhi = new TH2F("fHistJetHEtaPhi", "Jet-Hadron deta-dphi", 900, -1.8, 1.8, 720, -1.6, 4.8);

      fOutput->Add(fHistJetHTrackPt);
      fOutput->Add(fHistEPAngle);
      fOutput->Add(fHistJetEtaPhi);
      fOutput->Add(fHistJetHEtaPhi);

      TString name;
      for (Int_t trackPtBin = 0; trackPtBin < kMaxTrackPtBins; ++trackPtBin)
      {
        name = Form("fHistTrackEtaPhi_%i", trackPtBin);
        fHistTrackEtaPhi[trackPtBin] = new TH2F(name, name, 400, -1, 1, 720, 0.0, 2.0 * TMath::Pi());
        fOutput->Add(fHistTrackEtaPhi[trackPtBin]);
      }

      for (Int_t centralityBin = 0; centralityBin < kMaxCentralityBins; ++centralityBin)
      {
        name = Form("fHistJetPt_%i", centralityBin);
        fHistJetPt[centralityBin] = new TH1F(name, name, 240, -40, 200);
        fOutput->Add(fHistJetPt[centralityBin]);

        name = Form("fHistJetPtBias_%i", centralityBin);
        fHistJetPtBias[centralityBin] = new TH1F(name, name, 240, -40, 200);
        fOutput->Add(fHistJetPtBias[centralityBin]);
      }

      // Jet matching cuts
      // Only need if we actually jet matching
      if (fRequireMatchedJetWhenEmbedding)
      {
        std::vector<std::string> binLabels = {"noMatch", "matchedJet", "sharedMomentumFraction", "partLevelMatchedJet", "jetDistance", "passedAllCuts"};
        for (auto histName : std::vector<std::string>({"SameEvent", "MixedEvent"}))
        {
          name = std::string("fHistJetMatching") + histName.c_str() + "Cuts";
          std::string title = std::string("Jets which passed matching jet cuts for ") + histName;
          auto histMatchedJetCuts = fHistManager.CreateTH1(name, title.c_str(), binLabels.size(), 0, binLabels.size());
          // Set label names
          for (unsigned int i = 1; i <= binLabels.size(); i++)
          {
            histMatchedJetCuts->GetXaxis()->SetBinLabel(i, binLabels.at(i - 1).c_str());
          }
          histMatchedJetCuts->GetYaxis()->SetTitle("Number of jets");
        }
      }

      // Store hist manager output in the output list
      TIter next(fHistManager.GetListOfHistograms());
      TObject *obj = 0;
      while ((obj = next()))
      {
        fOutput->Add(obj);
      }
    }

    enum SparseAxes{
      centrality = 0,
      trigger_pT = 1,
      associated_pT = 2,
      delta_eta = 3,
      delta_phi = 4,
      is_leading_jet = 5,
      is_trigger_track = 6, // I have never used this, nor do I understand it. It is not used in the analysis.
      event_plane_angle = 7,
      z_vertex = 8,
      delta_R = 9,
      is_leading_track = 10,
      track_eta = 11,
      pion_TPC_nSigma = 12,
      pion_TOF_nSigma = 13,
      proton_TOF_nSigma = 14,
      kaon_TOF_nSigma = 15,
      electron_TOF_nSigma = 16,
      has_TOF_hit = 17,
    };

    void AliAnalysisTaskEmcalJetHdEdxCorrelations::CreateSparses(){
      // NOTE: The bit encoding doesn't preserve the order defined here. It's
      //       just using the bit values.
      UInt_t jetHadronAxesBitCode = 0; // bit coded, see GetDimParams() below
      if (fForceBeamType == AliAnalysisTaskEmcal::kpp)
      {
        vector<SparseAxes> jetHadronAxes = {
            centrality,
            trigger_pT,
            associated_pT,
            delta_eta,
            delta_phi,
            z_vertex,
            track_eta,
            pion_TPC_nSigma,
            pion_TOF_nSigma,
            proton_TOF_nSigma,
            kaon_TOF_nSigma,
            has_TOF_hit
            };
        for (auto axis : jetHadronAxes)
        {
          jetHadronAxesBitCode |= 1 << axis;
        }
      }
      else
      {
        vector<SparseAxes> jetHadronAxes = {
            centrality,
            trigger_pT,
            associated_pT,
            delta_eta,
            delta_phi,
            event_plane_angle,
            z_vertex,
            track_eta,
            pion_TPC_nSigma,
            pion_TOF_nSigma,
            proton_TOF_nSigma,
            kaon_TOF_nSigma,
            has_TOF_hit
            };
        for (auto axis : jetHadronAxes)
        {
          jetHadronAxesBitCode |= 1 << axis;
        }
      }
      fhnJH = NewTHnSparseF("fhnJH", jetHadronAxesBitCode);
      fhnJH->Sumw2();
      fOutput->Add(fhnJH);

      if (fDoEventMixing)
      {
        // The event plane angle does not need to be included because the semi-central determined that the EP angle didn't change
        // significantly for any of the EP orientations. However, it will be included so this can be demonstrated for the central
        // analysis if so desired.
        UInt_t mixedEventAxesBitCode = 0;
        if (fForceBeamType == AliAnalysisTaskEmcal::kpp)
        {
          vector<SparseAxes> mixedEventAxes = {
              centrality,
              trigger_pT,
              associated_pT,
              delta_eta,
              delta_phi,
              z_vertex,
              has_TOF_hit
              };
          for (auto axis : mixedEventAxes)
          {
            mixedEventAxesBitCode |= 1 << axis;
          }
        }
        else
        {
          vector<SparseAxes> mixedEventAxes = {
              centrality,
              trigger_pT,
              associated_pT,
              delta_eta,
              delta_phi,
              event_plane_angle,
              z_vertex,
              has_TOF_hit
              };
          for (auto axis : mixedEventAxes)
          {
            mixedEventAxesBitCode |= 1 << axis;
          }
        }
        fhnMixedEvents = NewTHnSparseF("fhnMixedEvents", mixedEventAxesBitCode);
        fhnMixedEvents->Sumw2();
        fOutput->Add(fhnMixedEvents);
      }

      // Trigger THnSparse
      UInt_t triggerAxesBitCode = 0;
      if (fForceBeamType == AliAnalysisTaskEmcal::kpp)
      {
        vector<SparseAxes> triggerAxes = {
            centrality,
            trigger_pT,
            z_vertex,
            };
        for(auto axis : triggerAxes){
          triggerAxesBitCode |= 1 << axis;
        }
      }
      else{
        vector<SparseAxes> triggerAxes = {
            centrality,
            trigger_pT,
            event_plane_angle,
            z_vertex,
            };
        for (auto axis : triggerAxes)
        {
          triggerAxesBitCode |= 1 << axis;
        }
      }
      fhnTrigger = NewTHnSparseF("fhnTrigger", triggerAxesBitCode);
      fhnTrigger->Sumw2();
      fOutput->Add(fhnTrigger);
    }

    void AliAnalysisTaskEmcalJetHdEdxCorrelations::CreateEventPool(){
      Int_t poolSize = -1; // Maximum number of events. Set to -1 to avoid limits on number of events
      // ZVertex
      Int_t nZVertexBins = 10;
      Double_t *zVertexBins = GenerateFixedBinArray(nZVertexBins, -10, 10);
      // Event activity (centrality of multiplicity)
      Int_t nEventActivityBins = 8;
      Double_t *eventActivityBins = 0;
      // +1 to accomodate the fact that we define bins rather than array entries.
      Double_t multiplicityBins[kMixedEventMultiplicityBins + 1] = {0., 4., 9., 15., 25., 35., 55., 100., 700.};

      // Cannot use GetBeamType() since it is not available until UserExec()
      if (fForceBeamType != AliAnalysisTaskEmcal::kpp)
      { // all besides pp
        // Event Activity is centrality in AA, pA
        nEventActivityBins = fNCentBinsMixedEvent;
        eventActivityBins = GenerateFixedBinArray(nEventActivityBins, 0, 100);
      }
      else if (fForceBeamType == AliAnalysisTaskEmcal::kpp)
      { // for pp only
        // Event Activity is multiplicity in pp
        eventActivityBins = multiplicityBins;
      }

      fPoolMgr = new AliEventPoolManager(poolSize, fNMixingTracks, nEventActivityBins, eventActivityBins, nZVertexBins, zVertexBins);

      // Print pool properties
      fPoolMgr->Validate();
    }

    /**
     * Perform run independent initializations, such as histograms and the event pool.
     */
    void AliAnalysisTaskEmcalJetHdEdxCorrelations::UserCreateOutputObjects()
    {
      // Call the base class implementation to get the base output.
      AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

      // Check that the task was initialized
      if (!fConfigurationInitialized)
      {
        AliFatal("Task was not initialized. Please ensure that Initialize() was called!");
      }
      // Reinitialize the YAML configuration
      fYAMLConfig.Reinitialize();

      CreateHistograms();

      CreateSparses();

      PostData(1, fOutput);

      // Event Mixing
      CreateEventPool();
    }

    /**
     * User specific initializations to perform before the first event.
     */
    void AliAnalysisTaskEmcalJetHdEdxCorrelations::UserExecOnce()
    {
      // Base class.
      AliAnalysisTaskEmcalJet::UserExecOnce();

      // Ensure that the random number generator is seeded in each job.
      fRandom.SetSeed(0);
    }

    /**
     * Get the proper bin based on the track pt value.
     *
     * @param[in] pt Track pt value to be binned.
     * @return Bin corresponding to the input value.
     */
    Int_t AliAnalysisTaskEmcalJetHdEdxCorrelations::GetTrackPtBin(Double_t pt) const
    {
      Int_t ptBin = -1;
      if (pt < 0.5)
        ptBin = 0;
      else if (pt < 1)
        ptBin = 1;
      else if (pt < 2)
        ptBin = 2;
      else if (pt < 3)
        ptBin = 3;
      else if (pt < 5)
        ptBin = 4;
      else if (pt < 8)
        ptBin = 5;
      else if (pt < 20)
        ptBin = 6;

      return ptBin;
    }

    /**
     * Retrieve the trigger mask from the input event or embedding helper as approapriate.
     *
     * @return The event trigger mask.
     */
    UInt_t AliAnalysisTaskEmcalJetHdEdxCorrelations::RetrieveTriggerMask() const
    {
      UInt_t eventTrigger = 0;
      if (fIsEmbedded)
      {
        auto embeddingHelper = AliAnalysisTaskEmcalEmbeddingHelper::GetInstance();
        if (embeddingHelper)
        {
          auto aodHeader = dynamic_cast<AliVAODHeader *>(embeddingHelper->GetEventHeader());
          if (aodHeader)
          {
            AliDebugStream(5) << "Retrieving trigger mask from embedded event\n";
            eventTrigger = aodHeader->GetOfflineTrigger();
          }
          else
          {
            AliErrorStream() << "Failed to retrieve requested AOD header from embedding helper\n";
          }
        }
        else
        {
          AliErrorStream() << "Failed to retrieve requested embedding helper\n";
        }
      }
      else
      {
        AliDebugStream(5) << "Retrieving trigger mask from internal event\n";
        eventTrigger = ((AliInputEventHandler *)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
      }

      return eventTrigger;
    }

    void AliAnalysisTaskEmcalJetHdEdxCorrelations::MixEvents(AliJetContainer *jets, AliPIDResponse *pidResponse, std::vector<unsigned int> rejectedTrackIndices, bool useListOfRejectedIndices, Int_t current_event_multiplicity, Double_t zVertex, UInt_t eventTrigger, Double_t flattened_EP_angle)
    {
        
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

      Double_t deltaPhi = 0;
      Double_t deltaEta = 0;
      Double_t deltaR = 0;
      Double_t epAngle = 0;
      Double_t eventActivity = 0;
      Double_t jetPt = 0;
      Float_t efficiency = -999;
      Bool_t leadJet = kFALSE;
      Bool_t isBiasedJet = kFALSE;
      AliPIDResponse::EDetPidStatus TOFPIDstatus;

      AliEventPool *pool = 0;
      if (fBeamType == kAA || fBeamType == kpA)
      { // everything but pp
        pool = fPoolMgr->GetEventPool(fCent, zVertex);
      }
      else if (fBeamType == kpp)
      { // pp only
        pool = fPoolMgr->GetEventPool(static_cast<Double_t>(current_event_multiplicity), zVertex);
      }

      if (!pool)
      {
        if (fBeamType == kAA || fBeamType == kpA)
          AliFatal(Form("No pool found for centrality = %f, zVertex = %f", fCent, zVertex));
        else if (fBeamType == kpp)
          AliFatal(Form("No pool found for ntracks_pp = %d, zVertex = %f", current_event_multiplicity, zVertex));
          return;
      }

      // The number of events in the pool
      Int_t nMix = pool->GetCurrentNEvents();
      Double_t rhoVal = jets->GetRhoVal();
      AliEmcalJet *leadingJet = jets->GetLeadingJet();
      // The two bitwise and to zero yet are still equal when both are 0, so we allow for that possibility
      if ((eventTrigger & fTriggerType) || eventTrigger == fTriggerType)
      {
        // check for a trigger jet
        if (pool->IsReady() || pool->NTracksInPool() >= fMinNTracksMixedEvents || nMix >= fMinNEventsMixedEvents)
        {

          for (auto jet : jets->accepted())
          {
            // Require the found jet to be matched
            // This match should be between detector and particle level MC
            if (fIsEmbedded && fRequireMatchedJetWhenEmbedding)
            {
              bool foundMatchedJet = CheckForMatchedJet(jets, jet, "fHistJetMatchingMixedEventCuts");
              if (foundMatchedJet == false)
              {
                continue;
              }
            }

            if (fBeamType == kAA || fBeamType == kpA)
            { // pA and AA
              eventActivity = fCent;
            }
            else if (fBeamType == kpp)
            {
              eventActivity = static_cast<Double_t>(current_event_multiplicity);
            }

            // Jet properties
            jetPt = AliAnalysisTaskEmcalJetHUtils::GetJetPt(jet, rhoVal);
            // Determine if we have the lead jet
            leadJet = kFALSE;
            if (jet == leadingJet)
            {
              leadJet = kTRUE;
            }
            isBiasedJet = BiasedJet(jet);
            // epAngle = PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalJetHUtils::RelativeEPAngle(jet->Phi(), fEPV0);
            // new way of getting qnvectors
            if (fBeamType != kpp)
            {
              epAngle = PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalJetHUtils::RelativeEPAngle(jet->Phi(), flattened_EP_angle);
            }

            // Make sure event contains a biased jet above our threshold (reduce stats of sparse)
            if (jetPt < 15 || isBiasedJet == kFALSE){
              continue;
            }

            // Fill mixed-event histos here
            for (Int_t jMix = 0; jMix < nMix; jMix++)
            {
              TObjArray *bgTracks = pool->GetEvent(jMix);

              for (Int_t ibg = 0; ibg < bgTracks->GetEntries(); ibg++)
              {
                
                const AliVTrack* bgTrack = dynamic_cast<const AliVTrack*>(bgTracks->At(ibg));
                if(!bgTrack){
                  AliError(Form("%s: Could not receive track %d in mixed event %d", GetName(), ibg, jMix));
                  continue;
                }

                // NOTE: We don't need to apply the artificial track inefficiency here because we already applied
                //       it when will filling into the event pool (in CloneAndReduceTrackList()).
                Bool_t hasTOFhit = kFALSE;
                AliPIDResponse* pidResponse = fInputHandler->GetPIDResponse();
                if(pidResponse){
                  TOFPIDstatus = pidResponse->CheckPIDStatus(AliPIDResponse::kTOF, bgTrack);
                  if(TOFPIDstatus == AliPIDResponse::kDetPidOk){
                    hasTOFhit = kTRUE;
                  }
                }


                // Calculate single particle tracking efficiency of mixed events for correlations
                efficiency = EffCorrection(bgTrack->Eta(), bgTrack->Pt());
                if(isinf(1.0/efficiency) || isnan(1.0/efficiency)){
                  cout<<efficiency<<" <- Efficiency for eta="<<bgTrack->Eta()<<" and pt="<<bgTrack->Pt()<<" Skipping"<<endl;
                  continue;
                }

                AliTLorentzVector temp_track;
                temp_track.SetPtEtaPhiM(bgTrack->Pt(), bgTrack->Eta(), bgTrack->Phi(), 0);
                // Phi is [-0.5*TMath::Pi(), 3*TMath::Pi()/2.]
                GetDeltaEtaDeltaPhiDeltaR(temp_track, jet, deltaEta, deltaPhi, deltaR);
                if (fBeamType != AliAnalysisTaskEmcal::kpp)
                {
                  if(nMix*efficiency==0){
                    cout<<nMix<<" <- nMix and "<<efficiency<<" <- Efficiency for eta="<<bgTrack->Eta()<<" and pt="<<bgTrack->Pt()<<" Skipping"<<endl;
                    continue;
                  }
		  double triggerEntries[] = {eventActivity, jetPt, bgTrack->Pt(), deltaEta, deltaPhi, epAngle, zVertex, (double)hasTOFhit};
                    FillHist(fhnMixedEvents, triggerEntries, 1. / (nMix * efficiency), fNoMixedEventJESCorrection);
                }
                else
                {
                  double triggerEntries[] = {eventActivity, jetPt, bgTrack->Pt(), deltaEta, deltaPhi, zVertex, (double)hasTOFhit};
                  FillHist(fhnMixedEvents, triggerEntries, 1. / (nMix * efficiency), fNoMixedEventJESCorrection);
                }
              }
            }
          }
        }
      }
      
      
      // The two bitwise and to zero yet are still equal when both are 0, so we allow for that possibility
      if ((eventTrigger & fMixingEventType) || eventTrigger == fMixingEventType)
      {
        AliTrackContainer *tracks = dynamic_cast<AliTrackContainer*>(GetParticleContainer("tracksForCorrelations"));

        TObjArray* tracksClone = new TObjArray;
        tracksClone->SetOwner(kTRUE);

        auto trackIter = tracks->accepted_momentum();
        for (auto particleIter = trackIter.begin(); particleIter != trackIter.end(); particleIter++)
        {
          AliVTrack *particle = (AliVTrack*)(dynamic_cast<const AliVTrack *>(particleIter->second))->Clone();
          if (!particle)
          {
            AliError(Form("%s: Could not receive track in updating track pool", GetName()));
            continue;
          }

          tracksClone->Add(particle);
        }

          // update pool if jet in event or not
          pool->UpdatePool(tracksClone);
        }
      return;
    }

    /**
     * Main loop called for each event by AliAnalysisTaskEmcal.
     */
    Bool_t AliAnalysisTaskEmcalJetHdEdxCorrelations::Run()
    {
      // NOTE: Clusters are never used directly in the task, so the container is neither created not retrieved
      // Retrieve tracks
      AliTrackContainer *tracks = static_cast<AliTrackContainer *>(GetParticleContainer("tracksForCorrelations"));
      if (!tracks)
      {
        AliError(Form("%s: Unable to retrieve tracks!", GetName()));
        return kFALSE;
      }

      // Retrieve jets
      AliJetContainer *jets = GetJetContainer(0);
      if (!jets)
      {
        AliError(Form("%s: Unable to retrieve jets!", GetName()));
        return kFALSE;
      }
    

      // Keep track of the tracks which are rejected with an aritificial track inefficiency
      std::vector<unsigned int> rejectedTrackIndices;
      bool useListOfRejectedIndices = false;

      // Get z vertex
      Double_t zVertex = fVertex[2];
      if(zVertex < -10 || zVertex > 10)
      {
        AliError(Form("%s: Rejecting event because zVertex = %f is out of range!", GetName(), zVertex));
        return kFALSE;
      }
      // Flags
      Bool_t isBiasedJet = kFALSE;
      Bool_t leadJet = kFALSE;
      // Jet pt
      Double_t jetPt = 0;
      // Rho
      // NOTE: Defaults to 0 if rho was not specified.
      Double_t rhoVal = jets->GetRhoVal();
      // Relative angles and distances
      Double_t deltaPhi = 0;
      Double_t deltaEta = 0;
      Double_t deltaR = 0;
      Double_t epAngle = 0;
      // Event plane angle from V0C
      Double_t EP_angle_from_calib = 0;
      Double_t flattenedEPangle = 0;
          // Event activity (centrality or multipilicity)
          Double_t eventActivity = 0;
      // Efficiency correction
      Double_t efficiency = -999;
      // For comparison to the current jet
      AliEmcalJet *leadingJet = jets->GetLeadingJet();
      // For getting the proper properties of tracks
      AliTLorentzVector track;
      AliPIDResponse::EDetPidStatus TOFPIDstatus;

      // Get PID Response
      AliPIDResponse *pidResponse = fInputHandler->GetPIDResponse();
      if (!pidResponse)
      {
        AliErrorStream() << "PID Response not available\n";
        return kFALSE;
      }



      // Determine the trigger for the current event
      UInt_t eventTrigger = RetrieveTriggerMask();

      // new way of getting qnvectors
      if (fBeamType != kpp)
      {
        fqnVectorReader = (AliAnalysisTaskJetQnVectors *)AliAnalysisManager::GetAnalysisManager()->GetTask("AliAnalysisTaskJetQnVectors");
        if (!fqnVectorReader)
        {
          printf("Error: No AliAnalysisTaskJetQnVectors");
          return 0;
        } // GetQnVectorReader
        EP_angle_from_calib = fqnVectorReader->GetEPangleV0C();
        // At this point we need to apply the flattening to the EP angle
	
	Double_t flattenedEPangle = GetFlattenedEPAngle(EP_angle_from_calib);	
        fHistEPAngle->Fill(flattenedEPangle);

	if(fMakeShiftCorrectionCalibHists)
	  return kTRUE; // do not continue to process the event, just go to the next one
      }

      AliDebugStream(5) << "Beginning main processing. Number of jets: " << jets->GetNJets() << ", accepted jets: " << jets->GetNAcceptedJets() << "\n";

      // Handle fast partition if selected
      if ((eventTrigger & AliVEvent::kFastOnly) && fDisableFastPartition)
      {
        AliDebugStream(4) << GetName() << ": Fast partition disabled\n";
        if (fGeneralHistograms)
        {
          fHistEventRejection->Fill("Fast Partition", 1);
        }
        return kFALSE;
      }

      for (auto jet : jets->accepted())
      {
        // Selects only events that we are interested in (ie triggered)
        if (!(eventTrigger & fTriggerType))
        {
          // The two bitwise and to zero yet are still equal when both are 0, so we allow for that possibility
          if (eventTrigger == fTriggerType && eventTrigger == 0)
          {
            AliDebugStream(5) << "Event accepted because the physics selection is \"0\".\n";
          }
          else
          {
            AliDebugStream(5) << "Rejected jets due to physics selection. Phys sel: " << std::bitset<32>(eventTrigger) << ", requested triggers: " << std::bitset<32>(fTriggerType) << " \n";
            // We can break here - the physics selection is not going to change within an event.
            break;
          }
        }

        AliDebugStream(5) << "Jet passed event selection!\nJet: " << jet->toString().Data() << "\n";

        // Require the found jet to be matched
        // This match should be between detector and particle level MC
        if (fIsEmbedded && fRequireMatchedJetWhenEmbedding)
        {
          bool foundMatchedJet = CheckForMatchedJet(jets, jet, "fHistJetMatchingSameEventCuts");
          if (foundMatchedJet == false)
          {
            continue;
          }
        }

        // Determine event activity
        if (fBeamType == kAA || fBeamType == kpA)
        {
          eventActivity = fCent;
        }
        else if (fBeamType == kpp)
        {
          eventActivity = static_cast<Double_t>(tracks->GetNTracks());
        }

        // Jet properties
        jetPt = AliAnalysisTaskEmcalJetHUtils::GetJetPt(jet, rhoVal);
        // Determine if we have the lead jet
        leadJet = kFALSE;
        if (jet == leadingJet)
          leadJet = kTRUE;
        isBiasedJet = BiasedJet(jet);
        // epAngle = PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalJetHUtils::RelativeEPAngle(jet->Phi(), fEPV0);
        // new way of getting qnvectors
        if (fBeamType != kpp)
        {
          epAngle = PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalJetHUtils::RelativeEPAngle(jet->Phi(), flattenedEPangle);
        }
        // Fill jet properties
        fHistJetEtaPhi->Fill(jet->Eta(), jet->Phi());
        FillHist(fHistJetPt[fCentBin], jetPt);
        if (isBiasedJet == kTRUE)
        {
          FillHist(fHistJetPtBias[fCentBin], jetPt);

          if (fBeamType != kpp)
          {
            const double triggerInfo[] = {eventActivity, jetPt, epAngle, zVertex};
            fhnTrigger->Fill(triggerInfo);
          }
          else
          {
            const double triggerInfo[] = {eventActivity, jetPt, zVertex};
            fhnTrigger->Fill(triggerInfo);
          }
        }

        // Cut on jet pt of 15 to reduce the size of the sparses
        if (jetPt > 15)
        {

          AliDebugStream(4) << "Passed min jet pt cut of 15. Jet: " << jet->toString().Data() << "\n";
          AliDebugStream(4) << "N accepted tracks: " << tracks->GetNAcceptedTracks() << "\n";
          auto tracksIter = tracks->accepted_momentum();
          for (auto trackIter = tracksIter.begin(); trackIter != tracksIter.end(); trackIter++)
          {
            // Get proper track properties
            track.Clear();
            track = trackIter->first;
            const AliVTrack *vTrack = dynamic_cast<const AliVTrack *>(trackIter->second);
            if (!vTrack)
            {
              AliErrorStream() << "Could not retrieve associated track from trackIter, skipping track.\n";
              continue;
            }

            // Artificial inefficiency
            // Note that we already randomly rejected tracks so that the same tracks will be rejected for the mixed events
            bool rejectParticle = CheckArtificialTrackEfficiency(trackIter.current_index(), rejectedTrackIndices, useListOfRejectedIndices);
            if (rejectParticle)
            {
              AliDebugStream(4) << "Track rejected in signal correlation loop.\n";
              continue;
            }

              GetDeltaEtaDeltaPhiDeltaR(track, jet, deltaEta, deltaPhi, deltaR);

              Double_t pionTPCnSigma;
              Double_t pionTOFnSigma;
              Double_t protonTOFnSigma;
              Double_t kaonTOFnSigma;

              Double_t hasTOFhit;

              pionTPCnSigma = pidResponse->NumberOfSigmasTPC(vTrack, (AliPID::EParticleType)2);

              TOFPIDstatus = pidResponse->CheckPIDStatus(AliPIDResponse::kTOF, vTrack);
              if (TOFPIDstatus == AliPIDResponse::kDetPidOk){
                hasTOFhit = 1;
              } else {
                hasTOFhit = 0;
              }
            
              pionTOFnSigma = pidResponse->NumberOfSigmasTOF(vTrack, (AliPID::EParticleType)2);
              protonTOFnSigma = pidResponse->NumberOfSigmasTOF(vTrack, (AliPID::EParticleType)4);
              kaonTOFnSigma = pidResponse->NumberOfSigmasTOF(vTrack, (AliPID::EParticleType)3);


              // Double_t TPCsignal = aodTrack->GetTPCsignal();

              // Fill track properties
              fHistJetHTrackPt->Fill(track.Pt());
              fHistJetHEtaPhi->Fill(deltaEta, deltaPhi);

              // Calculate single particle tracking efficiency for correlations
              efficiency = EffCorrection(track.Eta(), track.Pt());
              AliDebugStream(6) << "track eta: " << track.Eta() << ", track pt: " << track.Pt() << ", efficiency: " << efficiency << ", TPCsignal: "
                                << "\n";

              // get track eta
              Double_t trackEta = track.Eta();

              if (isBiasedJet == kTRUE)
              {
                if (fBeamType != kpp)
                {
                  if (fDoLessSparseAxes)
                  { // check if we want all dimensions
                    double triggerEntries[] = {eventActivity, jetPt, track.Pt(), deltaEta, deltaPhi, /*static_cast<Double_t>(leadJet),*/ epAngle, zVertex, trackEta, pionTPCnSigma, pionTOFnSigma, protonTOFnSigma, kaonTOFnSigma, hasTOFhit};
                    FillHist(fhnJH, triggerEntries, 1.0 / efficiency);
                  }
                  else
                  {
                    double triggerEntries[] = {eventActivity, jetPt, track.Pt(), deltaEta, deltaPhi, /*static_cast<Double_t>(leadJet),*/ epAngle, zVertex, deltaR, trackEta, pionTPCnSigma, pionTOFnSigma, protonTOFnSigma, kaonTOFnSigma, hasTOFhit};
                    FillHist(fhnJH, triggerEntries, 1.0 / efficiency);
                  }
                }
                else
                {
                  double triggerEntries[] = {eventActivity, jetPt, track.Pt(), deltaEta, deltaPhi, /*static_cast<Double_t>(leadJet),*/ zVertex, trackEta, pionTPCnSigma, pionTOFnSigma, protonTOFnSigma, kaonTOFnSigma};
                  FillHist(fhnJH, triggerEntries, 1.0 / efficiency);
                }
              }

            }// track loop

            // After one jet (and looping over whatever tracks are available in this event), we want to use the list of rejected indices,
            // both for the next possible signal jet in the event and for the mixed events
            AliDebugStream(4) << "Switching to list of rejected track indices. Number of indices: " << rejectedTrackIndices.size() << "\n";
            useListOfRejectedIndices = true;

          } // jet pt cut
        }   // jet loop
    
        // Prepare to do event mixing


        if (fDoEventMixing == kTRUE)
        {
          MixEvents(jets, pidResponse, rejectedTrackIndices, useListOfRejectedIndices, tracks->GetNTracks(), zVertex, eventTrigger, flattenedEPangle);
        } // end of event mixing

        return kTRUE;
      }

    /**
     * Determine if a jet passes the track or cluster bias and is therefore a "biased" jet.
     *
     * @param[in] jet The jet to be checked.
     * @return true if the jet is biased.
     */
    Bool_t AliAnalysisTaskEmcalJetHdEdxCorrelations::BiasedJet(AliEmcalJet *jet)
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
    void AliAnalysisTaskEmcalJetHdEdxCorrelations::GetDeltaEtaDeltaPhiDeltaR(AliTLorentzVector &particleOne, AliVParticle *particleTwo, Double_t &deltaEta, Double_t &deltaPhi, Double_t &deltaR)
    {
      // Define dPhi = jet.Phi() - particle.Phi() and similarly for dEta
      // Returns deltaPhi in symmetric range so that we can calculate DeltaR.
      deltaPhi = DeltaPhi(particleOne.Phi(), particleTwo->Phi(), -1.0 * TMath::Pi(), TMath::Pi());
      deltaEta = particleTwo->Eta() - particleOne.Eta();
      deltaR = TMath::Sqrt(deltaPhi * deltaPhi + deltaEta * deltaEta);

      // Adjust to the normal range after the DeltaR caluclation
      deltaPhi = DeltaPhi(particleTwo->Phi(), particleOne.Phi(), -0.5 * TMath::Pi(), 3 * TMath::Pi() / 2.);
    }

    /**
     * Check whether a track should be rejected due to artificial track inefficiency.
     *
     * @param[in] trackIndex Index of the current track.
     * @param[in] rejectedTrackIndices Vector of track indices which have been already rejected. Can be used to store the track index if rejected or reject a track if matched (depending on useRejectedList).
     * @param[in] useRejectedList If true, check if the current track index is in the rejected index list and reject if so. If false, will randomly reject the track according to the track inefficiency.
     * @return True if particle should be rejected
     */
    bool AliAnalysisTaskEmcalJetHdEdxCorrelations::CheckArtificialTrackEfficiency(unsigned int trackIndex, std::vector<unsigned int> &rejectedTrackIndices, bool useRejectedList)
    {
      bool returnValue = false;
      if (fArtificialTrackInefficiency < 1.)
      {
        if (useRejectedList)
        {
          if (std::find(rejectedTrackIndices.begin(), rejectedTrackIndices.end(), trackIndex) != rejectedTrackIndices.end())
          {
            AliDebugStream(4) << "Track " << trackIndex << " rejected due to artificial tracking inefficiency (from list)\n";
            returnValue = true;
          }
        }
        else
        {
          // Rejet randomly
          Double_t rnd = fRandom.Rndm();
          if (fArtificialTrackInefficiency < rnd)
          {
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
     * @param[in] jets Jet container corresponding to the jet to be checked
     * @param[in] jet Jet to be checked
     * @param[in] histName Name of the hist in the hist manager where QA information will be filled
     * @return true if the jet passes the criteria. false otherwise.
     */
    bool AliAnalysisTaskEmcalJetHdEdxCorrelations::CheckForMatchedJet(AliJetContainer *jets, AliEmcalJet *jet, const std::string &histName)
    {
      bool returnValue = false;
      if (jet->ClosestJet())
      {
        fHistManager.FillTH1(histName.c_str(), "matchedJet");
        returnValue = true;
        // TODO: Can it be merged with the function in JetHPerformance?
        AliDebugStream(4) << "Jet is matched!\nJet: " << jet->toString() << "\n";
        // Check shared momentum fraction
        // We explicitly want to use indices instead of geometric matching
        double sharedFraction = jets->GetFractionSharedPt(jet, nullptr);
        if (sharedFraction < fMinSharedMomentumFraction)
        {
          AliDebugStream(4) << "Jet rejected due to shared momentum fraction of " << sharedFraction << ", which is smaller than the min momentum fraction of " << fMinSharedMomentumFraction << "\n";
          returnValue = false;
        }
        else
        {
          AliDebugStream(4) << "Passed shared momentum fraction with value of " << sharedFraction << "\n";
          fHistManager.FillTH1(histName.c_str(), "sharedMomentumFraction");
        }

        if (fRequireMatchedPartLevelJet)
        {
          AliEmcalJet *detLevelJet = jet->ClosestJet();
          AliEmcalJet *partLevelJet = detLevelJet->ClosestJet();
          if (!partLevelJet)
          {
            AliDebugStream(4) << "Jet rejected due to no matching part level jet.\n";
            returnValue = false;
          }
          else
          {
            AliDebugStream(4) << "Det level jet has a required match to a part level jet.\n"
                              << "Part level jet: " << partLevelJet->toString() << "\n";
            fHistManager.FillTH1(histName.c_str(), "partLevelMatchedJet");
          }
        }

        // Only check matched jet distance if a value has been set
        if (fMaxMatchedJetDistance > 0)
        {
          double matchedJetDistance = jet->ClosestJetDistance();
          if (matchedJetDistance > fMaxMatchedJetDistance)
          {
            AliDebugStream(4) << "Jet rejected due to matching distance of " << matchedJetDistance << ", which is larger than the max distance of " << fMaxMatchedJetDistance << "\n";
            returnValue = false;
          }
          else
          {
            AliDebugStream(4) << "Jet passed distance cut with distance of " << matchedJetDistance << "\n";
            fHistManager.FillTH1(histName.c_str(), "jetDistance");
          }
        }

        // Record all cuts passed
        if (returnValue == true)
        {
          fHistManager.FillTH1(histName.c_str(), "passedAllCuts");
        }
      }
      else
      {
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
    THnSparse *AliAnalysisTaskEmcalJetHdEdxCorrelations::NewTHnSparseF(const char *name, UInt_t entries)
    {
      Int_t count = 0;
      UInt_t tmp = entries;
      while (tmp != 0)
      {
        count++;
        tmp = tmp & ~-tmp; // clear lowest bit
      }

      TString hnTitle(name);
      const Int_t dim = count;
      Int_t nbins[dim];
      Double_t xmin[dim];
      Double_t xmax[dim];

      Int_t i = 0;
      Int_t c = 0;
      while (c < dim && i < 32)
      {
        if (entries & (1 << i))
        {
          TString label("");
          GetDimParams(i, label, nbins[c], xmin[c], xmax[c]);
          hnTitle += Form(";%s", label.Data());
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
    void AliAnalysisTaskEmcalJetHdEdxCorrelations::GetDimParams(Int_t iEntry, TString &label, Int_t &nbins, Double_t &xmin, Double_t &xmax)
    {
      const Double_t pi = TMath::Pi();

      switch (iEntry)
      {

      case 0:
        label = "V0 centrality (%)";
        nbins = 10;
        xmin = 0.;
        xmax = 100.;
        // Adjust for pp, since we are retrieving multiplicity instead
        if (fForceBeamType == AliAnalysisTaskEmcal::kpp)
        {
          label = "Multiplicity";
          xmax = 200.;
        }
        break;

      case 1:
        label = "Jet p_{T}";
        nbins = 24;
        xmin = -40.;
        xmax = 200.;
        break;

      case 2:
        if (fDoWiderTrackBin)
        {
          label = "Track p_{T}";
          nbins = 40;
          xmin = 0.;
          xmax = 10.;
        }
        else
        {
          label = "Track p_{T}";
          nbins = 100;
          xmin = 0.;
          xmax = 10;
        }
        break;

      case 3:
        label = "#Delta#eta";
        nbins = 28;
        xmin = -1.4;
        xmax = 1.4;
        break;

      case 4:
        label = "#Delta#phi";
        nbins = 72;
        xmin = -0.5 * pi;
        xmax = 1.5 * pi;
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
        label = "Event plane angle";
        nbins = 3;
        xmin = 0;
        xmax = TMath::Pi() / 2.;
        break;

      case 8:
        label = "Z vertex (cm)";
        nbins = 10;
        xmin = -10;
        xmax = 10;
        break;

      case 9:
        label = "deltaR";
        nbins = 10;
        xmin = 0.;
        xmax = 5.0;
        break;

      case 10:
        label = "Leading track";
        nbins = 20;
        xmin = 0;
        xmax = 50;
        break;

      case 11:
        label = "Track eta";
        nbins = 50;
        xmin = -1.2;
        xmax = 1.2;
        break;

      case 12:
        label = "Pion TPC nSigma";
        nbins = 100;
        xmin = -10;
        xmax = 10;
        break;

      case 13:
        label = "Pion TOF nSigma";
        nbins = 60;
        xmin = -5;
        xmax = 5;
        break;

      case 14:
        label = "Proton TOF nSigma";
        nbins = 60;
        xmin = -5;
        xmax = 5;
        break;

      case 15:
        label = "Kaon TOF nSigma";
        nbins = 60;
        xmin = -5;
        xmax = 5;
        break;

      case 16:
        label = "Electron TOF nSigma";
        nbins = 60;
        xmin = -5;
        xmax = 5;
        break;

      case 17:
        label = "Has TOF hit";
        nbins = 2;
        xmin = 0;
        xmax = 2;
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
    TObjArray *AliAnalysisTaskEmcalJetHdEdxCorrelations::CloneAndReduceTrackList(std::vector<unsigned int> &rejectedTrackIndices, const bool useRejectedList)
    {
      // clones a track list by using AliBasicTrack which uses much less memory (used for event mixing)
      TObjArray *tracksClone = new TObjArray;
      tracksClone->SetOwner(kTRUE);

      // Loop over all tracks
      // TPair* particle;
      // AliBasicParticle *clone = 0;
      AliTrackContainer *tracks = static_cast<AliTrackContainer*>(GetTrackContainer("tracksForCorrelations"));

      auto particlesIter = tracks->accepted_momentum();
      for (auto particleIter = particlesIter.begin(); particleIter != particlesIter.end(); particleIter++)
      {
        // Get the aliVtrack and Get the TOF PID status
        AliVTrack *bgTrackSecond = dynamic_cast<AliVTrack *>(particleIter->second);
        if (!bgTrackSecond)
        {
          AliErrorStream() << "Could not retrieve associated track from particleIter in CloneAndReduceTrackList, skipping track.\n";
          continue;
        }

        TParameter<Bool_t> hasTOFhit = TParameter("hasTOFhit", kFALSE);
        AliPIDResponse *pidResponse = fInputHandler->GetPIDResponse();
        AliPIDResponse::EDetPidStatus TOFPIDstatus = pidResponse->CheckPIDStatus(AliPIDResponse::kTOF, bgTrackSecond);
        if (TOFPIDstatus == AliPIDResponse::kDetPidOk)
        {
          hasTOFhit = TParameter("hasTOFhit", kTRUE);
        }

        // particle = TPair(const_cast<const TObject*>(&(particleIter->first)), dynamic_cast<TObject*>(&hasTOFhit));

        // Artificial inefficiency
        bool rejectParticle = kFALSE;//CheckArtificialTrackEfficiency(particleIter->current_index(), rejectedTrackIndices, useRejectedList);
        if (rejectParticle)
        {
          AliDebugStream(4) << "Track rejected in CloneAndReduceTrackList()\n";
          continue;
        }

        // Fill some QA information about the tracks
        // Int_t trackPtBin = GetTrackPtBin(const_cast<const AliTLorentzVector*>(particle->Key())->Pt());
        // if (trackPtBin > -1)
        //   fHistTrackEtaPhi[trackPtBin]->Fill(const_cast<const AliTLorentzVector*>(particle->Key())->Eta(), const_cast<const AliTLorentzVector*>(particle->Key())->Phi());

        // // Create new particle
        // clone = new AliBasicParticle(particle->Eta(), particle->Phi(), particle->Pt(), particle->Charge());
        // // Set so that we can do comparisons using the IsEqual() function.
        // clone->SetUniqueID(particle->GetUniqueID());
        
        tracksClone->Add(bgTrackSecond);
      }

      return tracksClone;
    }

    /**
     * Utility function to apply the determine the single track efficiency.
     *
     * @param trackEta Eta of the track
     * @param trackPt pT of the track
     *
     * @return Track efficiency of the track (the entry in a histogram should be weighted as 1/(return value))
     */
    Double_t AliAnalysisTaskEmcalJetHdEdxCorrelations::EffCorrection(Double_t trackEta, Double_t trackPt) const
    {
      return AliAnalysisTaskEmcalJetHUtils::DetermineTrackingEfficiency(
          trackPt, trackEta, fCentBin, fEfficiencyPeriodIdentifier,
          "PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalJetHdEdxCorrelations");
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
    void AliAnalysisTaskEmcalJetHdEdxCorrelations::FillHist(TH1 *hist, Double_t fillValue, Double_t weight, Bool_t noCorrection)
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

        std::vector<Double_t> yBinsContent;
        AliDebug(3, TString::Format("%s: Attempt to access weights from JES correction hist %s with jet pt %f!", GetName(), hist->GetName(), fillValue));
        AccessSetOfYBinValues(fJESCorrectionHist, xBin, yBinsContent);
        AliDebug(3, TString::Format("weights size: %zd", yBinsContent.size()));

        // Loop over all possible bins to contribute.
        // If content is 0 then calling Fill won't make a difference
        for (Int_t index = 1; index <= fJESCorrectionHist->GetYaxis()->GetNbins(); index++)
        {
          // Don't bother trying to fill in the weight is 0
          if (yBinsContent.at(index - 1) > 0)
          {
            // Determine the value to fill based on the center of the bins.
            // This in principle allows the binning between the correction and hist to be different
            Double_t fillLocation = fJESCorrectionHist->GetYaxis()->GetBinCenter(index);
            AliDebug(4, TString::Format("fillLocation: %f, weight: %f", fillLocation, yBinsContent.at(index - 1)));
            // minus 1 since loop starts at 1
            hist->Fill(fillLocation, weight * yBinsContent.at(index - 1));
          }
        }

        // TEMP
        // hist->Draw();
        // END TEMP
      }
    }

    /**
     * Utility function to fill a histogram with a given weight. If requested, it will also apply the JES correction
     * derived weight to the hist. It assumes that the corrected jet pt value is located in the second element of the
     * fillValue (as defined in GetDimParams()). If that is changed or the first entry is not included, then this
     * function must be updated!
     *
     * @param[in] hist Histogram to be filled.
     * @param[in] fillValue Array of values to be filled into the hist.
     * @param[in] weight The weight to be applied when filling the hist.
     * @param[in] noCorrection True if the JES correction _should not_ be applied.
     */
    void AliAnalysisTaskEmcalJetHdEdxCorrelations::FillHist(THnSparse *hist, Double_t *fillValue, Double_t weight, Bool_t noCorrection)
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

        std::vector<Double_t> yBinsContent;
        AliDebug(3, TString::Format("%s: Attempt to access weights from JES correction hist %s with jet pt %f!", GetName(), hist->GetName(), jetPt));
        AccessSetOfYBinValues(fJESCorrectionHist, xBin, yBinsContent);
        AliDebug(3, TString::Format("weights size: %zd", yBinsContent.size()));

        // Loop over all possible bins to contribute.
        // If content is 0 then calling Fill won't make a difference
        for (Int_t index = 1; index <= fJESCorrectionHist->GetYaxis()->GetNbins(); index++)
        {
          // Don't bother trying to fill in the weight is 0
          if (yBinsContent.at(index - 1) > 0)
          {
            // Determine the value to fill based on the center of the bins.
            // This in principle allows the binning between the correction and hist to be different
            fillValue[1] = fJESCorrectionHist->GetYaxis()->GetBinCenter(index);
            AliDebug(4, TString::Format("fillValue[1]: %f, weight: %f", fillValue[1], yBinsContent.at(index - 1)));
            // minus 1 since loop starts at 1
            hist->Fill(fillValue, weight * yBinsContent.at(index - 1));
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
    void AliAnalysisTaskEmcalJetHdEdxCorrelations::AccessSetOfYBinValues(TH2D *hist, Int_t xBin, std::vector<Double_t> &yBinsContent, Double_t scaleFactor)
    {
      for (Int_t index = 1; index <= hist->GetYaxis()->GetNbins(); index++)
      {
        // yBinsContent[index-1] = hist->GetBinContent(hist->GetBin(xBin,index));
        yBinsContent.push_back(hist->GetBinContent(hist->GetBin(xBin, index)));

        if (scaleFactor >= 0)
        {
          // -1 since index starts at 1
          hist->SetBinContent(hist->GetBin(xBin, index), yBinsContent.at(index - 1) / scaleFactor);
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
    Bool_t AliAnalysisTaskEmcalJetHdEdxCorrelations::RetrieveAndInitializeJESCorrectionHist(TString filename, TString histName, Double_t trackBias, Double_t clusterBias)
    {
      // Initialize grid connection if necessary
      if (filename.Contains("alien://") && !gGrid)
      {
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
      if (trackBias + 0.1 < AliAnalysisTaskEmcalJetHdEdxCorrelations::kDisableBias)
      {
        histName = TString::Format("%s_Track%.2f", histName.Data(), trackBias);
      }
      if (clusterBias + 0.1 < AliAnalysisTaskEmcalJetHdEdxCorrelations::kDisableBias)
      {
        histName = TString::Format("%s_Clus%.2f", histName.Data(), clusterBias);
      }

      // Open file containing the correction
      TFile *jesCorrectionFile = TFile::Open(filename);
      if (!jesCorrectionFile || jesCorrectionFile->IsZombie())
      {
        AliError(TString::Format("%s: Could not open JES correction file %s", GetName(), filename.Data()));
        return kFALSE;
      }

      // Retrieve the histogram containing the correction and safely add it to the task.
      TH2D *JESCorrectionHist = dynamic_cast<TH2D *>(jesCorrectionFile->Get(histName.Data()));
      if (JESCorrectionHist)
      {
        AliInfo(TString::Format("%s: JES correction hist name \"%s\" loaded from file %s.", GetName(), histName.Data(), filename.Data()));
      }
      else
      {
        AliError(TString::Format("%s: JES correction hist name \"%s\" not found in file %s.", GetName(), histName.Data(), filename.Data()));

        // Attempt the base name instead of the formatted hist name
        JESCorrectionHist = dynamic_cast<TH2D *>(jesCorrectionFile->Get(histBaseName.Data()));
        if (JESCorrectionHist)
        {
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
      TH2D *tempHist = static_cast<TH2D *>(JESCorrectionHist->Clone());
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
      if (tempName.Index(tag) == -1)
      {
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

    Double_t AliAnalysisTaskEmcalJetHdEdxCorrelations::GetFlattenedEPAngle(Double_t uncorrectedAngle){
      if(fMakeShiftCorrectionCalibHists)
	{
	  if(fRunNumber != previousrunnumber)
	    {
	      hShiftCoeff = (TH2D*)fOutput->FindObject(Form("hShiftCoeff_run%i",fRunNumber));
	      hShiftEvents = (TH1I*)fOutput->FindObject(Form("hShiftEvents_run%i",fRunNumber));
	      if((!hShiftCoeff) || (!hShiftEvents))
		{
		  hShiftCoeff = new TH2D(Form("hShiftCoeff_run%i",fRunNumber),Form("hShiftCoeff_run%i",fRunNumber),10,0,100,6,0.5,6.5);
		  fOutput->Add(hShiftCoeff);
		  hShiftEvents = new TH1I(Form("hShiftEvents_run%i",fRunNumber),Form("hShiftEvents_run%i",fRunNumber),10,0,100);
		  fOutput->Add(hShiftEvents);
		}
	    }
	  hShiftCoeff->Fill(fCent,1,TMath::Cos(2.*uncorrectedAngle));
	  hShiftCoeff->Fill(fCent,2,TMath::Sin(2.*uncorrectedAngle));
	  hShiftCoeff->Fill(fCent,3,TMath::Cos(4.*uncorrectedAngle));
	  hShiftCoeff->Fill(fCent,4,TMath::Sin(4.*uncorrectedAngle));
	  hShiftCoeff->Fill(fCent,5,TMath::Cos(6.*uncorrectedAngle));
	  hShiftCoeff->Fill(fCent,6,TMath::Sin(6.*uncorrectedAngle));
	  hShiftEvents->Fill(fCent);
	  previousrunnumber = fRunNumber;
	  return uncorrectedAngle;
	}

      if(!fDoShiftCorrection) return uncorrectedAngle;
      if(!fEPshiftHistoList) return uncorrectedAngle;

      if(fRunNumber != previousrunnumber)
	{
	  hShiftCoeff = (TH2D*)fEPshiftHistoList->FindObject(Form("hShiftCoeff_run%i",fRunNumber));
	  hShiftEvents = (TH1I*)fEPshiftHistoList->FindObject(Form("hShiftEvents_run%i",fRunNumber));
	}
      
      if((!hShiftCoeff) || (!hShiftEvents)){
	AliWarningStream() << "Warning! Cannot find event plane shift calibration for run " << fRunNumber << "!\n";
	return uncorrectedAngle;
      }

      Double_t avCos2Psi = 0;
      Double_t avSin2Psi = 0;
      Double_t avCos4Psi = 0;
      Double_t avSin4Psi = 0;
      Double_t avCos6Psi = 0;
      Double_t avSin6Psi = 0;
      
      Int_t centbin = hShiftEvents->GetXaxis()->FindBin(fCent);
      avCos2Psi = hShiftCoeff->GetBinContent(centbin,1)/hShiftEvents->GetBinContent(centbin);
      avSin2Psi = hShiftCoeff->GetBinContent(centbin,2)/hShiftEvents->GetBinContent(centbin);
      avCos4Psi = hShiftCoeff->GetBinContent(centbin,3)/hShiftEvents->GetBinContent(centbin);
      avSin4Psi = hShiftCoeff->GetBinContent(centbin,4)/hShiftEvents->GetBinContent(centbin);
      avCos6Psi = hShiftCoeff->GetBinContent(centbin,5)/hShiftEvents->GetBinContent(centbin);
      avSin6Psi = hShiftCoeff->GetBinContent(centbin,6)/hShiftEvents->GetBinContent(centbin);

      previousrunnumber = fRunNumber;
      
      Double_t correctedAngle = uncorrectedAngle
	                        + avCos2Psi*TMath::Sin(2.*uncorrectedAngle) - avSin2Psi*TMath::Cos(2.*uncorrectedAngle)
	                        + (avCos4Psi*TMath::Sin(4.*uncorrectedAngle) - avSin4Psi*TMath::Cos(4.*uncorrectedAngle))/2.
	                        + (avCos6Psi*TMath::Sin(6.*uncorrectedAngle) - avSin6Psi*TMath::Cos(6.*uncorrectedAngle))/3.;

      if(correctedAngle < 0) correctedAngle += TMath::Pi();
      if(correctedAngle > TMath::Pi()) correctedAngle -= TMath::Pi();

      return correctedAngle;
    }

    // EPCorrectionsFileSetter
    void AliAnalysisTaskEmcalJetHdEdxCorrelations::OpenEPcorrectionsFile(const char* filename) {
      //Check if filename contains alien:// and connect to TGrid if it does
      TString *filenameStr = new TString(filename);
      if (filenameStr->Contains("alien://")) {
	TGrid::Connect("alien://");
      }
      TFile *fEPcorrectionFile = TFile::Open(filename);
      if (!fEPcorrectionFile) 
	AliError(Form("Could not open file %s", filename));
      fEPshiftHistoList = new TList();
      fEPshiftHistoList->SetOwner(kTRUE);
      TIter next(fEPcorrectionFile->GetListOfKeys());
      TKey *key;
      while ( (key = (TKey*)next()) )
	{
	  TString histoname = key->GetName();
	  if(histoname.BeginsWith("hShiftCoeff"))
	    {
	      TH2D* hShiftCoeffTemp = (TH2D*)((TH2D*)fEPcorrectionFile->Get(histoname.Data()))->Clone();
	      hShiftCoeffTemp->SetDirectory(0);
	      fEPshiftHistoList->Add(hShiftCoeffTemp);
	    }
	  if(histoname.BeginsWith("hShiftEvents"))
	    {
	      TH1I* hShiftEventsTemp = (TH1I*)((TH1I*)fEPcorrectionFile->Get(histoname.Data()))->Clone();
	      hShiftEventsTemp->SetDirectory(0);
	      fEPshiftHistoList->Add(hShiftEventsTemp);
	    }
	}
      fEPcorrectionFile->Close();
    }
    
    /**
     * AddTask for the jet-hadron task. We benefit for actually having compiled code, as opposed to
     * struggling with CINT.
     */
    AliAnalysisTaskEmcalJetHdEdxCorrelations *AliAnalysisTaskEmcalJetHdEdxCorrelations::AddTaskEmcalJetHdEdxCorrelations(
        const char *nTracks,
        const char *nCaloClusters,
        // Jet options
        const Double_t trackBias,
        const Double_t clusterBias,
        // Mixed event options
        const Int_t nTracksMixedEvent, // Additionally acts as a switch for enabling mixed events
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
        const Bool_t JESCorrection,
        const char *JESCorrectionFilename,
        const char *JESCorrectionHistName,
	Bool_t makeEPcorrectionsFile,
	Bool_t doEPshiftCorrection,
        const char *epCorrectionsFilename,
        const char *suffix)
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

      if (trackName == "usedefault")
      {
        trackName = AliEmcalContainerUtils::DetermineUseDefaultName(AliEmcalContainerUtils::kTrack);
      }

      if (clusName == "usedefault")
      {
        clusName = AliEmcalContainerUtils::DetermineUseDefaultName(AliEmcalContainerUtils::kCluster);
      }

      TString name("AliAnalysisTaskJetH");
      if (!trackName.IsNull())
      {
        name += TString::Format("_%s", trackName.Data());
      }
      if (!clusName.IsNull())
      {
        name += TString::Format("_%s", clusName.Data());
      }
      if (strcmp(suffix, "") != 0)
      {
        name += TString::Format("_%s", suffix);
      }

      AliAnalysisTaskEmcalJetHdEdxCorrelations *correlationTask = new AliAnalysisTaskEmcalJetHdEdxCorrelations(name);
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
      correlationTask->SetDoLessSparseAxes(lessSparseAxes);
      correlationTask->SetDoWiderTrackBin(widerTrackBin);
      // Corrections
      if (JESCorrection == kTRUE)
      {
        Bool_t result = correlationTask->RetrieveAndInitializeJESCorrectionHist(JESCorrectionFilename, JESCorrectionHistName, correlationTask->GetTrackBias(), correlationTask->GetClusterBias());
        if (!result)
        {
          AliErrorClass("Failed to successfully retrieve and initialize the JES correction! Task initialization continuing without JES correction (can be set manually later).");
        }
      }
    
      // shift correction calibration
      correlationTask->SetMakeShiftCorrectionCalibHists(makeEPcorrectionsFile);
      correlationTask->SetDoEPShiftCorrection(doEPshiftCorrection);
      if(doEPshiftCorrection)
	correlationTask->OpenEPcorrectionsFile(epCorrectionsFilename);
      
      //-------------------------------------------------------
      // Final settings, pass to manager and set the containers
      //-------------------------------------------------------

      mgr->AddTask(correlationTask);

      // Create containers for input/output
      mgr->ConnectInput(correlationTask, 0, mgr->GetCommonInputContainer());
      AliAnalysisDataContainer *cojeth =
          mgr->CreateContainer(correlationTask->GetName(), TList::Class(), AliAnalysisManager::kOutputContainer,
                               Form("%s", AliAnalysisManager::GetCommonFileName()));
      mgr->ConnectOutput(correlationTask, 1, cojeth);

      return correlationTask;
    }
    
    /**
     * Call after the AddTask to setup the standard Jet-h configuration.
     *
     * @return true if sucessfully configured for the standard configuration
     */
    bool AliAnalysisTaskEmcalJetHdEdxCorrelations::ConfigureForStandardAnalysis(std::string trackName,
                                                                                std::string clusName,
                                                                                const double jetConstituentPtCut,
                                                                                const double trackEta,
                                                                                const double jetRadius)
    {
      bool returnValue = false;
      AliInfoStream() << "Configuring Jet-H Correlations task for a standard analysis.\n";

      // Add Containers
      // Clusters
      if (clusName == "usedefault")
      {
        clusName = AliEmcalContainerUtils::DetermineUseDefaultName(AliEmcalContainerUtils::kCluster);
      }
      // For jet finding
      AliClusterContainer *clustersForJets = new AliClusterContainer(clusName.c_str());
      clustersForJets->SetName("clustersForJets");
      clustersForJets->SetMinE(jetConstituentPtCut);

      // Tracks
      // For jet finding
      if (trackName == "usedefault")
      {
        trackName = AliEmcalContainerUtils::DetermineUseDefaultName(AliEmcalContainerUtils::kTrack);
      }
      AliParticleContainer *particlesForJets = AliAnalysisTaskEmcalJetHUtils::CreateParticleOrTrackContainer(trackName.c_str());
      particlesForJets->SetName("particlesForJets");
      particlesForJets->SetMinPt(jetConstituentPtCut);
      particlesForJets->SetEtaLimits(-1.0 * trackEta, trackEta);

      // Don't need to adopt the container - we'll just use it to find the right jet collection
      // For correlations
      AliParticleContainer *particlesForCorrelations = AliAnalysisTaskEmcalJetHUtils::CreateParticleOrTrackContainer(trackName.c_str());
      if (particlesForCorrelations)
      {

        particlesForCorrelations->SetName("tracksForCorrelations");
        particlesForCorrelations->SetMinPt(0.15);
        particlesForCorrelations->SetEtaLimits(-1.0 * trackEta, trackEta);

        // Adopt the container
        this->AdoptParticleContainer(particlesForCorrelations);
      }
      else
      {
        AliWarningStream() << "No particle container was successfully created!\n";
      }

      // Jets
      AliJetContainer *jetContainer = this->AddJetContainer(AliJetContainer::kFullJet,
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
    bool AliAnalysisTaskEmcalJetHdEdxCorrelations::ConfigureForEmbeddingAnalysis(std::string trackName,
                                                                                 std::string clusName,
                                                                                 const double jetConstituentPtCut,
                                                                                 const double trackEta,
                                                                                 const double jetRadius,
                                                                                 const std::string &jetTag,
                                                                                 const std::string &correlationsTracksCutsPeriod)
    {
      bool returnValue = false;
      AliInfoStream() << "Configuring Jet-H Correlations task for an embedding analysis.\n";

      // Set the task to know it that is embedded
      this->SetIsEmbedded(true);

      // Add Containers
      // Clusters
      if (clusName == "usedefault")
      {
        clusName = AliEmcalContainerUtils::DetermineUseDefaultName(AliEmcalContainerUtils::kCluster);
      }
      // For jet finding
      AliClusterContainer *clustersForJets = new AliClusterContainer(clusName.c_str());
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
      if (trackName == "usedefault")
      {
        trackName = AliEmcalContainerUtils::DetermineUseDefaultName(AliEmcalContainerUtils::kTrack);
      }
      AliParticleContainer *particlesForJets = AliAnalysisTaskEmcalJetHUtils::CreateParticleOrTrackContainer(trackName.c_str());
      particlesForJets->SetName("particlesForJets");
      particlesForJets->SetMinPt(jetConstituentPtCut);
      particlesForJets->SetEtaLimits(-1.0 * trackEta, trackEta);
      // Don't need to adopt the container - we'll just use it to find the right jet collection
      // For correlations
      AliParticleContainer *particlesForCorrelations = AliAnalysisTaskEmcalJetHUtils::CreateParticleOrTrackContainer(trackName.c_str());
      // Ensure that we don't operate on a null pointer
      if (particlesForCorrelations)
      {
        particlesForCorrelations->SetName("tracksForCorrelations");
        particlesForCorrelations->SetMinPt(0.15);
        particlesForCorrelations->SetEtaLimits(-1.0 * trackEta, trackEta);
        particlesForCorrelations->SetIsEmbedding(true);
        AliTrackContainer *trackCont = dynamic_cast<AliTrackContainer *>(particlesForCorrelations);
        if (trackCont)
        {
          // This option only exists for track containers
          trackCont->SetTrackCutsPeriod(correlationsTracksCutsPeriod.c_str());
        }
        // Adopt the container
        this->AdoptParticleContainer(particlesForCorrelations);
      }
      else
      {
        AliWarningStream() << "No particle container was successfully created!\n";
      }

      // Jets
      // The tag "hybridLevelJets" is defined in the jet finder
      AliJetContainer *jetContainer = this->AddJetContainer(AliJetContainer::kFullJet,
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
     * Prints information about the jet-hadron correlations task.
     *
     * @return std::string containing information about the task.
     */
    std::string AliAnalysisTaskEmcalJetHdEdxCorrelations::toString() const
    {
      std::stringstream tempSS;
      tempSS << std::boolalpha;
      tempSS << "Jet collections:\n";
      TIter next(&fJetCollArray);
      AliJetContainer *jetCont;
      while ((jetCont = static_cast<AliJetContainer *>(next())))
      {
        tempSS << "\t" << jetCont->GetName() << ": " << jetCont->GetArrayName() << "\n";
      }
      tempSS << "Event selection\n";
      tempSS << "\tUse AliEventCuts: " << !fUseBuiltinEventSelection << "\n";
      tempSS << "\tTrigger event selection: " << std::bitset<32>(fTriggerType) << "\n";
      tempSS << "\tMixed event selection: " << std::bitset<32>(fMixingEventType) << "\n";
      tempSS << "\tEnabled only for non-fast partition: " << fDisableFastPartition << "\n";
      tempSS << "Jet settings:\n";
      tempSS << "\tTrack bias: " << fTrackBias << "\n";
      tempSS << "\tCluster bias: " << fClusterBias << "\n";
      tempSS << "Event mixing:\n";
      tempSS << "\tEnabled: " << fDoEventMixing << "\n";
      tempSS << "\tN mixed tracks: " << fNMixingTracks << "\n";
      tempSS << "\tMin number of tracks for mixing: " << fMinNTracksMixedEvents << "\n";
      tempSS << "\tMin number of events for mixing: " << fMinNEventsMixedEvents << "\n";
      tempSS << "\tNumber of centrality bins for mixing: " << fNCentBinsMixedEvent << "\n";
      tempSS << "Histogramming options:\n";
      tempSS << "\tLess sparse axes: " << fDoLessSparseAxes << "\n";
      tempSS << "\tWider associated track pt bins: " << fDoWiderTrackBin << "\n";
      tempSS << "Jet matching options for embedding:\n";
      tempSS << "\tRequire the jet to match to a embedded jets: " << fRequireMatchedJetWhenEmbedding << "\n";
      tempSS << "\tRequire an additional match to a part level jet: " << fRequireMatchedPartLevelJet << "\n";
      tempSS << "\tMinimum shared momentum fraction: " << fMinSharedMomentumFraction << "\n";
      tempSS << "\tMax matched jet distance: " << fMaxMatchedJetDistance << "\n";
      tempSS << "Efficiency\n";
      tempSS << "\tSingle track efficiency identifier: " << fEfficiencyPeriodIdentifier << "\n";
      tempSS << "\tArtifical track inefficiency: " << fArtificialTrackInefficiency << "\n";

      return tempSS.str();
    }

    /**
     * Print jet-hadron correlations task information on an output stream using the string representation provided by
     * AliAnalysisTaskEmcalJetHdEdxCorrelations::toString. Used by operator<<
     * @param in output stream stream
     * @return reference to the output stream
     */
    std::ostream &AliAnalysisTaskEmcalJetHdEdxCorrelations::Print(std::ostream &in) const
    {
      in << toString();
      return in;
    }

    /**
     * Print basic jet-hadron correlations task information using the string representation provided by
     * AliAnalysisTaskEmcalJetHdEdxCorrelations::toString
     *
     * @param[in] opt Unused
     */
    void AliAnalysisTaskEmcalJetHdEdxCorrelations::Print(Option_t *opt) const
    {
      Printf("%s", toString().c_str());
    }

  } /* namespace EMCALJetTasks */
} /* namespace PWGJE */

/**
 * Implementation of the output stream operator for AliAnalysisTaskEmcalJetHdEdxCorrelations. Printing
 * basic jet-hadron correlations task information provided by function toString
 * @param in output stream
 * @param myTask Task which will be printed
 * @return Reference to the output stream
 */
std::ostream &operator<<(std::ostream &in, const PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalJetHdEdxCorrelations &myTask)
{
  std::ostream &result = myTask.Print(in);
  return result;
}
