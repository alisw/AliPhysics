#ifndef ALICSEVENTCUTS_H
#define ALICSEVENTCUTS_H

/// \file AliCSEventCuts.h
/// \brief Event cuts support for the correlations studies tasks
///
/// Based on ideas taken from Gamma Conversion analysis under PWGGA and from AliESDTrackCuts

#include "AliAnalysisUtils.h"
#include "AliCSAnalysisCutsBase.h"

class TH1F;
class TH2F;
class TF1;
class TFormula;
class AliESDtrackCuts;


/// \class AliCSEventCuts
/// \brief Class which implements event cuts
///
/// Provides support for configuring event selection cuts.
/// It is derived from #AliCSAnalysisCutsBase so it incorporates
/// the cuts string in all output histograms.
///
/// Based on Gamma Conversion analysis implementation within PWGGA
///
/// \author Víctor González <victor.gonzalez@cern.ch>, UCM
/// \date Oct 17, 2016

class AliCSEventCuts : public AliCSAnalysisCutsBase {
public:

  /// \enum cutsParametersIds
  /// \brief The ids of the different event cuts parameters supported
  enum cutsParametersIds {
    kSystem,            ///< system type
    kCentralityType,    ///< centrality estimator
    kCentralityMin,     ///< centrality cut minimum value
    kCentralityMax,     ///< centrality cut maximum value
    kActiveTrigger,     ///< trigger selection
    kActiveSubTrigger,  ///< sub-trigger selection
    kRemovePileUp,      ///< pile up removal
    kVertex,            ///< vertex cut
    kRemove2015PileUp,  ///< 2015 PbPb additional pile up removal
    kNCutsParameters    ///< the number of event cuts parameters
  };

  /// \enum cutsIds
  /// \brief The ids of the different cuts currently supported
  enum cutsIds {
    kMCdataQuality,             ///< MC analysis without proper MC data
    kDAQIncompleteCut,          ///< Incomplete data event cut
    kNoTracks,                  ///< Event with no tracks
    kOfflineTriggerCut,         ///< Offline trigger cut
    kVertexContributorsCut,     ///< Vertex contributor cut
    kVertexQualityCut,          ///< Vertex quality cut
    kSPDTrackVtxDistance,       ///< Distance between SPD and tracks vertex cut
    kVertexCut,                 ///< Vertex z cut
    kPileUpCut,                 ///< Pile up cut
    k2015PileUpCut,             ///< 2015 PbPb additional pile up cut
    kSPDClsVsTrkaletsCut,       ///< SPD clusters vs tracklets cut
    kCentralityCut,             ///< Centrality cut
    kNCuts                      ///< The number of supported cuts
  };

  /// \enum SystemType
  /// \brief The type of the system under analysis
  enum SystemType {
    kNoSystem = 0,    ///< no system defined
    kpp,          ///< **p-p** system
    kpPb,         ///< **p-Pb** system
    kPbPb,        ///< **Pb-Pb** system
    kXeXe,        ///< **Xe-Xe** system
    kPbp,         ///< **Pb-p** system
    knSystems     ///< number of handled systems
  };

private:
  /// \enum baseSystems
  /// \brief The ids of the systems used as base for the filter track counting cuts
  enum baseSystems {
    kUnknownBase,                    ///< no system base
    k2010based,                      ///< 2010 based system
    k2011based                       ///< 2011 based system
  };

public:

                      AliCSEventCuts();
                      AliCSEventCuts(const char *name, const char * title);
  virtual            ~AliCSEventCuts();

  virtual void        InitCuts(const char *name = NULL);
  virtual void        NotifyRun();
  virtual void        NotifyEvent();
  virtual Bool_t      IsEventAccepted(AliVEvent *fInputEvent);

                      /// Get the system type
                      /// \return the system type being analysed
  SystemType          GetSystem() const { return fSystem; }
                      /// Get the accepted centrality range
                      /// \param minmax storage for the return values: index 0 minimum, index 1 maximum
  void                GetCentralityRange(Float_t minmax[2]) const { minmax[0] = fCentralityMin; minmax[1] = fCentralityMax; }
                      /// Get the vertex \f$z\f$ coordinate
                      /// \return the vertex \f$z\f$ coordinate
  Double_t            GetVertexZ() { return fVertexZ; }
                      /// Gets the centrality of the current event
                      /// \return event centrality in percentage
  Double_t            GetCentrality() { return fCentrality; }

private:
  /* we set them private to force cuts string consistency */
  Bool_t              SetSystemType(SystemType system);
  Bool_t              SetCentralityType(Int_t code);
  Bool_t              SetCentralityMin(Int_t code);
  Bool_t              SetCentralityMax(Int_t code);
  Bool_t              SetSelectActiveTrigger(Int_t code);
  Bool_t              SetSelectActiveSubTrigger(Int_t code);
  Bool_t              SetRemovePileUp(Int_t code);
  Bool_t              SetVertexCut(Int_t code);
  Bool_t              SetRemove2015PileUp(Int_t code);


private:
  virtual Bool_t      SetCutAndParams(Int_t paramID, Int_t value);
  virtual void        PrintCutWithParams(Int_t paramID) const;
  virtual void        PrintTrigger(UInt_t &printed, UInt_t trigger, const char *name) const;

  Int_t               GetNumberOfVertexContributors(AliVEvent *event) const;
  Bool_t              PassVertexResolutionAndDispersionTh(AliVEvent *event) const;
  Bool_t              AcceptSPDTracksVtxDist(AliVEvent *event) const;
  Bool_t              Is2015PileUpEvent() const;
  Double_t            GetVertexZ(AliVEvent *event) const;
  Bool_t              UseNewMultiplicityFramework() const;
  Float_t             GetEventCentrality(AliVEvent *event) const;
  Float_t             GetEventAltCentrality(AliVEvent *event) const;
  Bool_t              StoreEventCentralities(AliVEvent *event);
  Bool_t              StoreEventMultiplicities(AliVEvent *event);
  void                GetCentralityEstimatorNames(const char *&sel, const char *&alt) const;
  void                SetActualActiveTrigger();
  void                SetActualSystemType();
  void                SetActualVertexQuality();
  void                SetActual2015PileUpRemoval();
  void                SetActualFilterTracksCuts();

  void                DefineHistograms();


private:
  static const char  *fgkCutsNames[kNCuts];                         ///< the names of the different event cuts
  static Float_t      fgkVertexResolutionThreshold;                 ///< the vertex resolution threshold default value
  static Float_t      fgkVertexResolutionThreshold_pPb;             ///< the vertex resolution threshold default value for pPb
  static Float_t      fgkVertexDispersionThreshold;                 ///< the vertex dispersion threshold default value (for ESD only)
  static Float_t      fgkSPDTracksVtxDistanceThreshold;             ///< the threshold for the distance between SPD and tracks vertex
  static Float_t      fgkSPDTracksVtxDistanceThreshold_pPb;         ///< the threshold for the distance between SPD and tracks vertex for pPb
  static Float_t      fgkSPDTracksVtxDistanceSigmas;                ///< number of tolerated sigmas for the total distance between SPD and tracks vertex
  static Float_t      fgkSPDTracksVtxDistanceSigmas_pPb;            ///< number of tolerated sigmas for the total distance between SPD and tracks vertex for pPb
  static Float_t      fgkTrackVertexSigmas;                         ///< number of tolerated track vertex sigmas for the distance between SPD and tracks vertex
  static Float_t      fgkTrackVertexSigmas_pPb;                     ///< number of tolerated track vertex sigmas for the distance between SPD and tracks vertex for pPb

  SystemType          fSystem;                ///< the type of system being analyzed
  Double_t            fVertexZ;               ///< the vertex \f$z\f$ coordinate
  Double_t            fCentrality;            ///< the event centrality
  Double_t            fAltCentrality;         ///< the event centrality from the alternate detector
  Int_t               fCentralityDetector;    ///< the detector to estimate the centrality
  Int_t               fCentralityModifier;    ///< the modifier of the centrality cut value
  Float_t             fCentralityMin;         ///< the minimum value for centrality cut
  Float_t             fCentralityMax;         ///< the maximum value for centrality cut
  UInt_t              fOfflineTriggerMask;    ///< the selected trigger mask
  Float_t             fMaxVertexZ;            ///< the maximum value of the absolute z vertex coordinate
  Bool_t              fUseSPDTracksVtxDist;   ///< check the distance between SPD and tracks vertex
  Float_t             fVertexResolutionTh;    ///< vertex resolution threshold
  Float_t             fVertexDispersionTh;    ///< vertex dispersion threshold (for ESD only)
  Float_t             fSPDTrkVtxDistTh;       ///< SPD tracks vertexes distance threshold
  Float_t             fSPDTrkVtxDistSigmas;   ///< n total sigmas for the SPD tracks vertexes distance
  Float_t             fTrkVtxDistSigmas;      ///< track vertex n sigmas for the SPD tracks vertexes distance
  Bool_t              fUseNewMultFramework;   ///< kTRUE if the new multiplicity framework for centrality estimation must be used
  TFormula           *fRun2V0MBasedPileUpCorrelation;    ///< formula to evaluate Run2 additional pileup cut
  TF1*                fCentOutLowCut;         ///< cut low for centrality outliers
  TF1*                fCentOutHighCut;        ///< cut high for centrality outliers
  TF1*                fTOFMultOutLowCut;      ///< cut low for TOF multiplicity outliers
  TF1*                fTOFMultOutHighCut;     ///< cut high for TOF multiplicity outliers
  TF1*                fMultCentOutLowCut;     ///< cut low for multiplicity centrality outliers
  Float_t             fV0MCentrality;         ///< the event V0M centrality
  Float_t             fV0ACentrality;         ///< the event V0A centrality
  Float_t             fV0CCentrality;         ///< the event V0C centrality
  Float_t             fCL0Centrality;         ///< the event CL0 centrality
  Float_t             fCL1Centrality;         ///< the event CL1 centrality
  Int_t               fReferenceMultiplicity; ///< event reference multiplicity
  Int_t               fV0Multiplicity;        ///< the event V0 multiplicity
  Int_t               fNoOfAODTracks;         ///< the number of AOD tracks
  Int_t               fNoOfESDTracks;         ///< the number of ESD tracks
  Int_t               fNoOfFB32Tracks;        ///< the number of globals tracks with tight DCA
  Int_t               fNoOfFB128Tracks;       ///< the number of TPC only able to be constrained to SPD vertex tracks
  Int_t               fNoOfFB32AccTracks;     ///< the number of global tracks with tight DCA accepted by standard cuts
  Int_t               fNoOfFB32TOFTracks;     ///< the number of global tracks with tight DCA acceptable by TOF
  Int_t               fNoOfTPCoutTracks;      ///< the number of tracks with TPCout flag on
  Int_t               fNoOfInitialTPCoutTracks;      ///< the number of tracks with TPCout flag on, initial track counting method
  Int_t               fNoOfTotalTPCClusters;  ///< the total number of TPC clusters for the event


  AliAnalysisUtils    fAnalysisUtils;         ///< analysis utilities for pile up detection
  AliESDtrackCuts    *fESDFB32;               ///< ESD tracks cuts to mirror AOD FB32
  AliESDtrackCuts    *fESDFB128;              ///< ESD tracks cuts to mirror AOD FB128

  /* histograms with two indices: before (0) / after (1) applying the cut */
  TH1F               *fhCutsStatistics;                ///< the cuts statistics
  TH1F               *fhUniqueCutsStatistics;          ///< the unique cuts statistics
  TH2F               *fhCutsCorrelation;               ///< cuts correlation
  TH1F               *fhCentrality[2];                 ///< the event centrality histogram (b/a)
  TH1F               *fhVertexZ[2];                    ///< the event vertex z histograms (b/a)
  TH2F               *fhSPDClustersVsTracklets[2];     ///< SPD clusters vs number of tracklets histogram (b/a)
  TH2F               *fhV0MvsTracksTPCout[2];          ///< V0 multiplicity vs number of TPCout tracks histogram (b/a)
  TH2F               *fhV0MvsTracksInitialTPCout[2];   ///< V0 multiplicity vs number of initial method TPCout tracks histogram (b/a)
  TH2F               *fhV0MvsTotalTPCClusters[2];      ///< V0 multiplicity vs number of total TPC clusters (b/a)
  TH2F               *fhCentralityAltVsSel[2];         ///< Centrality correlation alternative vs selected detector
  TH2F               *fhCL0vsV0MCentrality[2];         ///< Centrality correlation CL0 vs V0M
  TH2F               *fhESDvsTPConlyMultiplicity[2];   ///< Multiplicity ESD tracks vs TPC only tracks
  TH2F               *fhTOFvsGlobalMultiplicity[2];    ///< Multiplicity TOF global tracks vs global tracks
  TH2F               *fhAccTrkvsV0MCentrality[2];      ///< Multiplicity accepted global tracks vs V0M centraliy
  TH1F               *fhTriggerClass[2];               ///< the fired triggers (b/a)

  /// Copy constructor
  /// Not allowed. Forced private.
  AliCSEventCuts(const AliCSEventCuts&);
  /// Assignment operator
  /// Not allowed. Forced private.
  /// \return l-value reference object
  AliCSEventCuts& operator=(const AliCSEventCuts&);

  /// \cond CLASSIMP
  ClassDef(AliCSEventCuts,11);
  /// \endcond
};

#endif /* ALICSEVENTCUTS_H */
