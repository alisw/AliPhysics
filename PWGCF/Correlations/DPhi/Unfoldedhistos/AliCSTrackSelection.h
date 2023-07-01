#ifndef ALICSTRACKSELECTION_H
#define ALICSTRACKSELECTION_H

/// \file AliCSTrackSelection.h
/// \brief Track selection cuts support for the correlations studies tasks
///

#include <TNamed.h>
#include <TObjArray.h>
#include "AliCSTrackCutsBase.h"
#include "AliPID.h"

class TBits;
class TH1F;
class TH2F;
class TH3F;
class AliVTrack;
class AliVParticle;
class AliPIDResponse;

/// \class AliCSTrackSelection
/// \brief Class which implements the whole set of track selection cuts
///
/// The track selection is based in two set of cuts one of them is the
/// inclusive set and the other is the exclusive one. For a potential track
/// to be accepted it has to pass at least one of the inclusive cuts and
/// none of the exclusive ones. Otherwise the track will be rejected.
///
/// \author Víctor González <victor.gonzalez@cern.ch>, UCM
/// \date Feb 12, 2017

class AliCSTrackSelection : public TNamed {
 public:
  /// \enum poiIds
  /// \brief The internal ids of the particles/tracks of interest
  enum poiIds {
    kWrongPOIid = -1, ///< wrong id
    kPOIh = 0,        ///< hadrons
    kPOIpi,           ///< pions
    kPOIka,           ///< kaons
    kPOIpr,           ///< protons
    kNPOIids          ///< the number of POD ids
  };

 public:
  AliCSTrackSelection();
  AliCSTrackSelection(const char* name, const char* title);
  virtual                           ~AliCSTrackSelection();

  Bool_t                             InitializeFromString(const char *confstring);
                                     /// Sets the desired level for the QA histograms output
                                     /// \param level the desired QA histograms output level
  void                               SetQALevelOutput(AliCSAnalysisCutsBase::QALevel level) { fQALevel = level; }
                                     /// Get the histograms list
                                     /// \return the histograms list
  virtual TList                     *GetHistogramsList() { return fHistogramsList; }


  virtual void                       InitCuts(const char *name = NULL);
  virtual void                       NotifyRun();
  virtual void                       NotifyEvent();
  virtual Bool_t                     IsTrackAccepted(AliVTrack *trk);
  virtual poiIds GetAcceptedTrackPOIId() { return fTrackId; };
  virtual Bool_t                     IsTrueTrackAccepted(AliVTrack *trk);
  virtual Bool_t                     IsTrueTrackAccepted(Int_t itrk);
  virtual poiIds GetAcceptedParticlePOIId() { return fParticleId; };

  static Bool_t IsTruePrimary(AliVTrack* trk);
  Bool_t                             IsFromMCInjectedSignal(Int_t itrk);
  static poiIds poiid(AliPID::EParticleType sp);
  std::vector<std::string> getPoiNames() { return fPoiNames; }

 private:
  void DefineHistograms();

 private:
  std::vector<std::string> fPoiNames;             ///< the names of the particles of interest, i.e. {"Ha"."Pi","Ka","Pr"}
  poiIds fTrackId;                                ///< the track PID encoded as POI Id
  poiIds fParticleId;                             ///< the particle PID encoded as POI Id
  AliCSAnalysisCutsBase::QALevel fQALevel;        ///< the level of the QA histograms output
  TObjArray fInclusiveTrackCuts;                  ///< the set of inclusive track cuts
  TObjArray fInclusivePIDCuts;                    ///< the set of inclusive PID track cuts
  TObjArray fInclusiveCutsStrings;                ///< the strings to configure the inclusive track cuts
  TObjArray fInclusivePidCutsStrings;             ///< the strings to configure the PID inclusive cuts
  TObjArray fExclusiveCuts;                       ///< the set of exclusive track cuts
  TObjArray fExclusiveCutsStrings;                ///< the strings to configure the exclusive track cuts
  TObjArray fExclusivePidCutsStrings;             ///< the strings to configure the PID exclusive cuts
  TBits* fCutsActivatedMask;                      ///< the mask of cut activated for the ongoing event
  TString fCutsString;                            ///< the string that configures the set of cuts

  AliCSAnalysisCutsBase::ProdPeriods fDataPeriod; ///< the current period under analysis
  Int_t fHighestPrimaryLabel;                     ///< the highest primary label accepted when MC signals are injected
  AliPIDResponse* fPIDResponse;                   ///< the PID response instance

  /* histograms with two indices: before (0) / after (1) applying the cut */
  TH1F* fhCutsStatistics;             ///< the cuts statistics
  TH2F* fhCutsCorrelation;            ///< cuts correlation
  TH2F* fhPtVsDCAxy[2];               ///< \f$ p_{T} \f$ vs \f$ \mbox{DCA}_{XY} \f$ (b/a)
  TH2F* fhPtVsDCAz[2];                ///< \f$ p_{T} \f$ vs \f$ \mbox{DCA}_{Z} \f$ (b/a)
  TH2F* fhPtVsTPCCls[2];              ///< \f$ p_{T} \f$ vs TPC clusters (b/a)
  TH2F* fhPtVsTPCRows[2];             ///< \f$ p_{T} \f$ vs TPC crossed rows (b/a)
  TH2F* fhPtVsTPCRowOverFindCls[2];   ///< \f$ p_{T} \f$ vs TPC crossed rows findable clusters ratio (b/a)
  TH2F* fhEtaVsPhi[2];                ///< \f$ p_{T} \f$ vs \f$ \eta \f$ (b/a)
  TH2F* fhPtVsEta[2];                 ///< \f$ p_{T} \f$ vs \f$ \eta \f$ (b/a)
  TH2F* fhITSdEdxSignalVsP[2];        ///< ITS dE/dx signal vs P (b/a)
  TH2F* fhTPCdEdxSignalVsP[2];        ///< TPC dE/dx signal vs P(b/a)
  TH2F* fhTPCdEdxSelSignalVsP[3];     ///< TPC dE/dx signal vs P for the different selected species
  TH2F* fhTOFSignalVsP[2];            ///< TOF signal vs P(b/a)
  TH2F* fhTOFSelSignalVsP[3];         ///< TOF signal vs P for the different selected species
  TH3F* fhTPCTOFSigmaVsP[3][2];       ///< TOF n sigmas vs TPC n sigmas vs P  to the pi, k, p lines (b/a)
  TH2F* fhTPCdEdxSignalDiffVsP[3][2]; ///< TPC dE/dx signal difference to the pi, k, p lines vs P (b/a)
  TH2F* fhTOFSignalDiffVsP[3][2];     ///< TOF signal difference to the pi, k, p lines vs P (b/a)
  TH2F* fhTPCdEdxSelSignalDiffVsP
    [3][3];                           ///< TPC dE/dx signal difference to the pi, k, p lines vs P  for the different selected species
  TH2F* fhTOFSelSignalDiffVsP
    [3][3];                           ///< TOF signal difference to the pi, k, p lines vs P  for the different selected species
  TH2F* fhPvsTOFMassSq[2];            ///< momentum vs TOF estimated squared mass (b/a)
  TH2F* fhSelPvsTOFMassSq[3];         ///< momentum vs TOF estimated squared mass for the different selected species

  TList* fHistogramsList;             ///< the list of histograms used

  /// Copy constructor
  /// Not allowed. Forced private.
  AliCSTrackSelection(const AliCSTrackSelection&);
  /// Assignment operator
  /// Not allowed. Forced private.
  /// \return l-value reference object
  AliCSTrackSelection& operator=(const AliCSTrackSelection&);

  /// \cond CLASSIMP
  ClassDef(AliCSTrackSelection, 1);
  /// \endcond
};



#endif /* ALICSTRACKSELECTION_H */
