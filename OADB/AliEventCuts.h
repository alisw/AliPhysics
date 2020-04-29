#ifndef _AliEventCuts_h_
#define _AliEventCuts_h_

#include <TF1.h>
#include <TList.h>
#include <TNamed.h>
#include <TString.h> /// required to have easy access to tokenize
#include <cmath>
#include <string>
#include <vector>

#include "AliVEvent.h"
#include "AliAnalysisCuts.h"
#include "AliAnalysisUtils.h"
#include "AliTimeRangeMasking.h"
#include "AliTimeRangeCut.h"
#include "AliEMCALLEDEventsCut.h"

class AliESDtrackCuts;
class TList;
class TH1D;
class TH1I;
class TH2D;
class TH2F;

class AliEventCutsContainer : public TNamed {
  public:
    AliEventCutsContainer() : TNamed("AliEventCutsContainer","AliEventCutsContainer"),
    fEventId(0u),
    fMultESD(-1),
    fMultTrkFB32(-1),
    fMultTrkFB32Acc(-1),
    fMultTrkFB32TOF(-1),
    fMultTrkTPC(-1),
    fMultTrkTPCout(-1),
    fMultVZERO(-1.) {}

    unsigned long fEventId;
    int fMultESD;
    int fMultTrkFB32;
    int fMultTrkFB32Acc;
    int fMultTrkFB32TOF;
    int fMultTrkTPC;
    int fMultTrkTPCout;
    double fMultVZERO;
  ClassDef(AliEventCutsContainer,2)
};

class AliEventCuts : public TList {
  public:
    AliEventCuts(bool savePlots = false);
    virtual ~AliEventCuts();
    enum CutsBin {
      kNoCuts = 0,
      kDAQincomplete,
      kBfield,
      kTrigger,
      kTriggerClasses,
      kVertexSPD,
      kVertexTracks,
      kVertex,
      kVertexPositionSPD,
      kVertexPositionTracks,
      kVertexPosition,
      kVertexQuality,
      kPileUp,
      kMultiplicity,
      kINELgt0,
      kCorrelations,
      kTimeRangeCut,
      kEMCALEDCut,
      kAllCuts
    };

    enum NormMask {
      kAnyEvent = BIT(kNoCuts),
      kTriggeredEvent = BIT(kTrigger),
      kPassesAllCuts = (BIT(kAllCuts) - 1) ^ (BIT(kVertexPositionSPD) | BIT(kVertexPositionTracks) | BIT(kVertexSPD) | BIT(kVertexTracks) | BIT(kTriggerClasses)),
      kPassesNonVertexRelatedSelections = kPassesAllCuts ^ (BIT(kVertex) | BIT(kVertexPosition) | BIT(kVertexQuality)),
      kHasReconstructedVertex = kPassesAllCuts ^ BIT(kVertexPosition)
    };


    bool   AcceptEvent (AliVEvent *ev);
    bool   PassedCut (AliEventCuts::CutsBin cut) { return fFlag & BIT(cut); }
    bool   CheckNormalisationMask (AliEventCuts::NormMask mask) { return (fFlag & mask) == mask; }
    bool   IsTrueINELgtZero (AliVEvent *ev, bool chkGenVtxZ = false);
    void   AddQAplotsToList(TList *qaList = 0x0, bool addCorrelationPlots = false);
    void   OverrideAutomaticTriggerSelection(unsigned long tr, bool ov = true) { fTriggerMask = tr; fOverrideAutoTriggerMask = ov; }
    void   OverridePileUpCuts(int minContrib, float minZdist, float nSigmaZdist, float nSigmaDiamXY, float nSigmaDiamZ, bool ov = true);
    void   OverrideCentralityFramework(int centFramework = 0) { fOverrideCentralityFramework = true; fCentralityFramework = centFramework; }
    void   SetManualMode (bool man = true) { fManualMode = man; }
    void   SetupRun1PbPb();
    void   SetupLHC15o() { SetupRun2PbPb(); }
    void   SetupPbPb2018();
    void   SetupRun2PbPb();
    void   SetupLHC17n();
    void   SetupRun2pp();
    void   SetupRun1pp();
    void   SetupRun1pA(int iPeriod);
    void   SetupRun2pA(int iPeriod);
    void   UseMultSelectionEventSelection(bool useIt = true);
    void   SetAcceptedTriggerClasses(TString classes);

    static bool GoodPrimaryAODVertex(AliVEvent *ev);

    /// set up the usage of the time range cut
    void UseTimeRangeCut() { fUseTimeRangeCut = true;}

    ///
    const AliTimeRangeCut& GetTimeRangeCut() const { return fTimeRangeCut; }

    /// set up the usage of the time range cut
    void UseEMCALLEDEventsCut() { fUseEMCALLEDEventsCut = true;}

    ///
    AliEMCALLEDEventsCut& GetEMCALLEDEventsCut() { return fEMCALLEDEventsCut; }
  
    /// While the general philosophy here is to avoid setters and getters
    /// for some variables (like the max z vertex position) standard the cuts usually follow some patterns
    /// (e.g. the max vertex z position is always symmetric wrt the nominal beam collision point) thus
    /// is convenient having special setters in this case.
    float             GetCentrality (unsigned int estimator = 0) const;
    std::string       GetCentralityEstimator (unsigned int estimator = 0) const;
    const AliVVertex* GetPrimaryVertex() const { return fPrimaryVertex; }

    void          SetCentralityEstimators (std::string first = "V0M", std::string second = "CL0") { fCentEstimators[0] = first; fCentEstimators[1] = second; }
    void          SetCentralityRange (float min, float max) { fMinCentrality = min; fMaxCentrality = max; }
    void          SetMaxVertexZposition (float max) { fMinVtz = -fabs(max); fMaxVtz = fabs(max); }
    void          SelectOnlyInelGt0(bool toogle) { fOverrideInelGt0 = true; fSelectInelGt0 = toogle; }

    AliAnalysisUtils fUtils;                      ///< Analysis utils for the pileup rejection

    bool          fGreenLight;                    ///< If true it will bypass all the selections.
    bool          fMC;                            ///< Set to true by the automatic setup when analysing MC (in manual mode *you* are responsible for it). In MC the correlations cuts are disabled.

    bool          fRequireTrackVertex;            ///< if true all the events with only the SPD vertex are rejected
    float         fMinVtz;                        ///< Min z position for the primary vertex
    float         fMaxVtz;                        ///< Max z position for the primary vertex
    float         fMaxDeltaSpdTrackAbsolute;      ///<
    float         fMaxDeltaSpdTrackNsigmaSPD;     ///<
    float         fMaxDeltaSpdTrackNsigmaTrack;   ///<
    float         fMaxResolutionSPDvertex;        ///<
    float         fMaxDispersionSPDvertex;        ///<
    bool          fCheckAODvertex;                ///< if true it rejects the AOD primary vertices coming from TPC ESD vertices or SPD placeholder vertices

    bool          fRejectDAQincomplete;           ///< Reject events that have incomplete information

    int           fRequiredSolenoidPolarity;      ///< 0: does not require any particular polarity. Positive numbers -> positive B field, negative numbers -> negative B field

    bool          fUseCombinedMVSPDcut;           ///< If true, MV cut is used when a track vertex is present and if not, it used the SPDvsMultCut

    bool          fUseMultiplicityDependentPileUpCuts; ///< If true fSPDpileupMinContributors is set according to the event multiplicity (automatic setup set it true only if the user does not specify any custom value)
    bool          fUseSPDpileUpCut;               ///< Enable the SPD pileup cut
    int           fSPDpileupMinContributors;      ///< Reject all the events with SPD pile-up vertices with more than fRejectPileupSPD contributors
    double        fSPDpileupMinZdist;             ///<
    double        fSPDpileupNsigmaZdist;          ///<
    double        fSPDpileupNsigmaDiamXY;         ///<
    double        fSPDpileupNsigmaDiamZ;          ///<
    bool          fTrackletBGcut;                 ///<

    bool          fPileUpCutMV;                   ///< Reject pile-up based on multi-vertexer

    unsigned int  fCentralityFramework;           ///< 0: skip centrality checks, 1: multiplicity framework, 2: legacy centrality framework
    float         fMinCentrality;                 ///< Minimum centrality to be analised
    float         fMaxCentrality;                 ///< Maximum centrality to be analised

    bool          fUseVariablesCorrelationCuts;   ///< Switch on/off the cuts on the correlation between event variables
    bool          fUseEstimatorsCorrelationCut;   ///< Switch on/off the cut on the correlation between centrality estimators
    bool          fUseStrongVarCorrelationCut;    ///< Switch on/off the strong cuts on the correlation between event variables
    bool          fUseITSTPCCluCorrelationCut;    ///< Switch on/off the strong cuts on the correlation between event variables
    double        fEstimatorsCorrelationCoef[2];  ///< fCentEstimators[0] = [0] + [1] * fCentEstimators[1]
    double        fEstimatorsSigmaPars[4];        ///< Sigma parametrisation fCentEstimators[1] vs fCentEstimators[0]
    double        fDeltaEstimatorNsigma[2];       ///< Number of sigma to cut on fCentEstimators[1] vs fCentEstimators[0]
    double        fTOFvsFB32correlationPars[4];   ///< Polynomial parametrisation of the correlation between FB32+TOF tracks and FB32 tracks
    double        fTOFvsFB32sigmaPars[6];         ///< Polynomial parametrisation of the rms of the correlation between FB32+TOF tracks and FB32 tracks
    double        fTOFvsFB32nSigmaCut[2];         ///< N sigma cut in the FB32+TOF vs FB32 plane
    double        fESDvsTPConlyLinearCut[2];      ///< Linear cut in the ESD track vs TPC only track plane
    TF1          *fMultiplicityV0McorrCut;        //!<! Cut on the FB128 vs V0M plane
    double        fFB128vsTrklLinearCut[2];       ///< Cut on the FB128 vs Tracklet plane
    double        fVZEROvsTPCoutPolCut[5];        ///< Cut on VZERO multipliciy vs the number of tracks with kTPCout on
    double        fITSvsTPCcluPolCut[3];          ///< Cut on SDD+SSD clusters vs TPC clusters
    
    bool          fRequireExactTriggerMask;       ///< If true the event selection mask is required to be equal to fTriggerMask
    unsigned long fTriggerMask;                   ///< Trigger mask
    std::vector<std::string> fTriggerClasses;     ///< Trigger classes
  
    AliEventCutsContainer fContainer;       //!<! Local copy of the event cuts container (safe against user changes)
    const std::string  fkLabels[2];                    ///< Histograms labels (raw/selected)

  private:
    AliEventCuts(const AliEventCuts& copy);
    AliEventCuts operator=(const AliEventCuts& copy);
    void          AutomaticSetup (AliVEvent *ev);
    void          ComputeTrackMultiplicity(AliVEvent *ev);
    template<typename F> F PolN(F x, F* coef, int n);

    bool          fManualMode;                    ///< if true the cuts are not loaded automatically looking at the run number
    bool          fSavePlots;                     ///< if true the plots are automatically added to this object
    int           fCurrentRun;                    ///<
    unsigned long fFlag;                          ///< Flag of the passed cuts

    std::string   fCentEstimators[2];             ///< Centrality estimators: the first is used as main estimators, that is correlated with the second to monitor spurious events.
    float         fCentPercentiles[2];            ///< Centrality percentiles
    AliVVertex   *fPrimaryVertex;                 //!<! Primary vertex pointer

    ///
    bool          fNewEvent;                      ///<  True if the AliVEvent identifier in the AcceptEvent and fIdentifier are different
    /// Overrides
    bool          fOverrideAutoTriggerMask;       ///<  If true the trigger mask chosen by the user is not overridden by the Automatic Setup
    bool          fOverrideAutoPileUpCuts;        ///<  If true the pile-up cuts are defined by the user.
    bool          fMultSelectionEvCuts;           ///< Enable/Disable the event selection applied in the AliMultSelection framework
    bool          fUseTimeRangeCut;               ///< If to use the time range cut
    bool          fUseEMCALLEDEventsCut;          ///< If EMCal LED events are to be removed

    bool          fSelectInelGt0;                 ///< Select only INEL > 0 events
    bool          fOverrideInelGt0;               ///< If the user ask for a configuration, let's not touch it
    bool          fOverrideCentralityFramework;   ///< If the user ask (not) to run a centrality framework this should be onored by AliEventCuts 

    AliTimeRangeCut fTimeRangeCut;       ///< Time Range cut
  
    AliEMCALLEDEventsCut fEMCALLEDEventsCut;      ///< EMCAL LED events cut

    /// The following pointers are used to avoid the intense usage of FindObject. The objects pointed are owned by (TList*)this.
    TH1D* fCutStats;               //!<! Cuts statistics: every column keeps track of how many times a cut is passed independently from the other cuts.
    TH1D* fCutStatsAfterTrigger;   //!<! Cuts statistics: every column keeps track of how many times a cut is passed for events that pass the trigger selection
    TH1D* fCutStatsAfterMultSelection; //!<! Cuts statistics: every column keeps track of how many times a cut is passed for events that pass the multiplicity selection
    TH1D* fNormalisationHist;      //!<! Cuts statistics: every column keeps track of how many times a cut is passed once that all the other are passed.
    TH1D* fVtz[2];                 //!<! Vertex z distribution
    TH1D* fDeltaTrackSPDvtz[2];    //!<! Difference between the vertex computed using SPD and the track vertex
    TH1D* fCentrality[2];          //!<! Centrality percentile distribution
    TH2D* fEstimCorrelation[2];    //!<! Correlation between centrality estimators
    TH2D* fMultCentCorrelation[2]; //!<! Correlation between main centrality estimator and multiplicity

    TH2F* fTOFvsFB32[2];           //!<!
    TH2F* fTPCvsAll[2];            //!<!
    TH2F* fMultvsV0M[2];           //!<!
    TH2F* fTPCvsTrkl[2];           //!<!
    TH2F* fVZEROvsTPCout[2];       //!<!

    AliESDtrackCuts* fFB32trackCuts; //!<! Cuts corresponding to FB32 in the ESD (used only for correlations cuts in ESDs)
    AliESDtrackCuts* fTPConlyCuts;   //!<! Cuts corresponding to the standalone TPC cuts in the ESDs (used only for correlations cuts in ESDs)

    ClassDef(AliEventCuts, 14)
};

template<typename F> F AliEventCuts::PolN(F x,F* coef, int n) {
  if (n < 1) ::Fatal("AliEventCuts::PolN","PolN should be used only for n>=1.");
  F ret = coef[0] + coef[1] * x;
  for (int i = 2; i <= n; ++i)
    ret += coef[i] * std::pow(x,i);
  return ret;
}

#endif

