#ifndef ALICSTRACKCUTS_H
#define ALICSTRACKCUTS_H

/// \file AliCSTrackCuts.h
/// \brief Track selection cuts support for the correlations studies tasks
///
/// Based on ideas taken from Gamma Conversion analysis under PWGGA and from AliESDTrackCuts

#include "AliAnalysisUtils.h"
#include "AliCSTrackCutsBase.h"

class TH1F;
class TH2F;
class AliESDtrackCuts;
class AliVTrack;

/// \class AliCSTrackCuts
/// \brief Class which implements track selection cuts
///
/// Provides support for configuring track selection cuts.
/// It is derived from #AliCSAnalysisCutsBase so it incorporates
/// the cuts string in all output histograms.
///
/// Based on Gamma Conversion analysis implementation within PWGGA
///
/// \author Víctor González <victor.gonzalez@cern.ch>, UCM
/// \date Oct 17, 2016

class AliCSTrackCuts : public AliCSTrackCutsBase {
public:

  /// \enum cutsParametersIds
  /// \brief The ids of the different track cuts parameters supported
  enum cutsParametersIds {
    kTrackTypeCutParam,              ///< The type of track cut parameters
    kEtaCutParam,                    ///< \f$ \eta \f$ cut parameters
    kITSClsCutParam,                 ///< ITS clusters cut parameters
    kTPCClsCutParam,                 ///< TPC clusters cut parameters
    kDCACutParam,                    ///< DCA cut parameters
    kPtCutParam,                     ///< \f$ p_{T} \f$ cut parameters
    kNCutsParameters                 ///< the number of track cuts parameters
  };

  /// \enum cutsIds
  /// \brief The ids of the different track cuts currently supported
  enum cutsIds {
    kNotConstrainedCut,              ///< not constrained information cut
    kEtaCut,                         ///< \f$ \eta \f$ cut
    kTrackTypeCuts,                  ///< ESD track cuts machinery. Includes ITS clusters, TPC clusters and DCA cuts
    kPtCut,                          ///< \f$ p_{T} \f$ cut
    kNCuts                           ///< The number of supported cuts
  };

private:
  /// \enum baseSystems
  /// \brief The ids of the systems used as base for the track cuts
  enum baseSystems {
    kUnknownBase,                    ///< no system base
    k2010based,                      ///< 2010 based system
    k2011based                       ///< 2011 based system
  };

public:

                                     AliCSTrackCuts();
                                     AliCSTrackCuts(const char *name, const char * title);
  virtual                           ~AliCSTrackCuts();

  virtual void                       InitCuts(const char *name = NULL);
  virtual void                       NotifyRun();
  virtual void                       NotifyEvent();
  virtual Bool_t                     IsTrackAccepted(AliVTrack *trk);
  virtual Bool_t                     IsTrueTrackAccepted(AliVTrack *trk);
  virtual Bool_t                     IsTrueTrackAccepted(Int_t itrk);
                                     /// Does the current track need to be constrained
                                     /// \return kTRUE if the track needs to be constrained kFALSE otherwise
  virtual Bool_t                     GetConstrain() { return fConstrain; }

  static Bool_t                      IsTruePrimary(AliVTrack *trk);
  static Bool_t                      IsPhysicalPrimary(Int_t itrk);

private:
  virtual Bool_t                     AcceptTrackType(AliVTrack *trk, AliVTrack *&ttrk);
  void                               SetActualTypeOfTrackCuts();
  void                               SetActualITSClustersCut();
  void                               SetActualTPCClustersCut();
  void                               SetActualDCACut();

  void                               DefineHistograms();

private:
  /* we set them private to force cuts string consistency */
  Bool_t                             SetTypeOfTrackCut(Int_t ttype);
  Bool_t                             SetEtaCut(Int_t etacode);
  Bool_t                             SetITSClustersCut(Int_t ITSclscode);
  Bool_t                             SetTPCClustersCut(Int_t TPCclscode);
  Bool_t                             SetDCACut(Int_t dcacode);
  Bool_t                             SetPtCut(Int_t ptcode);

private:
  virtual Bool_t                     SetCutAndParams(Int_t paramID, Int_t value);
  virtual void                       PrintCutWithParams(Int_t paramID) const;
  void                               PrintTypeOfTrackCut() const;
  void                               PrintITSClustersCut() const;
  void                               PrintTPCClustersCut() const;
  void                               PrintDCACut() const;

private:
  static const char                 *fgkCutsParametersNames[kNCutsParameters];     ///< the names of the different event cuts parameters
  static const char                 *fgkCutsNames[kNCuts];                         ///< the names of the different event cuts

  Bool_t                             fConstrain;                         ///< the current track has been constrained so, it needs to be constrained

  baseSystems                        fBaseSystem;                        ///< the id of the system used as base for the cuts in the current analysis
  Double_t                           fEtaCut;                            ///< the maximum value of \f$ | \eta | \f$
  Double_t                           fMinPt;                             ///< the minimum \f$ p_{T} \f$ value
  Double_t                           fMaxPt;                             ///< the maximum \f$ p_{T} \f$ value
  Double_t                           fEtaShift;                          ///< shift in pseudo-rapidity when asymmetric system (eg. pPb)

  AliESDtrackCuts                   *fESDTrackCuts;                      ///< the ESD track cuts instance
  UInt_t                             fAODFilterBits;                     ///< the AOD track filter bits selected


  /* histograms with two indices: before (0) / after (1) applying the cut */
  TH1F                              *fhCutsStatistics;                   ///< the cuts statistics
  TH2F                              *fhCutsCorrelation;                  ///< cuts correlation
  TH2F                              *fhPtVsDCAxy[2];                     ///< \f$ p_{T} \f$ vs \f$ \mbox{DCA}_{XY} \f$ (b/a)
  TH2F                              *fhPtVsDCAz[2];                      ///< \f$ p_{T} \f$ vs \f$ \mbox{DCA}_{Z} \f$ (b/a)
  TH2F                              *fhPtVsTPCCls[2];                    ///< \f$ p_{T} \f$ vs TPC clusters (b/a)
  TH2F                              *fhPtVsTPCRowOverFindCls[2];         ///< \f$ p_{T} \f$ vs TPC crossed rows findable clusters ratio (b/a)
  TH2F                              *fhEtaVsPhi[2];                      ///< \f$ p_{T} \f$ vs \f$ \eta \f$ (b/a)
  TH2F                              *fhPtVsEta[2];                       ///< \f$ p_{T} \f$ vs \f$ \eta \f$ (b/a)

  /// Copy constructor
  /// Not allowed. Forced private.
  AliCSTrackCuts(const AliCSTrackCuts&);
  /// Assignment operator
  /// Not allowed. Forced private.
  /// \return l-value reference object
  AliCSTrackCuts& operator=(const AliCSTrackCuts&);

  /// \cond CLASSIMP
  ClassDef(AliCSTrackCuts,3);
  /// \endcond
};



#endif /* ALICSTRACKCUTS_H */
