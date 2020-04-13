#ifndef ALICSPIDCUTS_H
#define ALICSPIDCUTS_H

/// \file AliCSPIDCuts.h
/// \brief PID track selection cuts support for the correlations studies tasks
///
/// Based on ideas taken from Gamma Conversion analysis under PWGGA and from AliESDTrackCuts

#include <TBits.h>
#include "AliPID.h"
#include "AliCSTrackCutsBase.h"

class TH1F;
class TH2F;
class AliVTrack;
class AliVParticle;
class AliPIDResponse;

/// \class AliCSPIDCuts
/// \brief Class which implements PID track selection cuts
///
/// Provides support for configuring PID track selection cuts.
/// It is derived from #AliCSAnalysisCutsBase so it incorporates
/// the cuts string in all output histograms.
///
/// There are two basic use cases
/// - identify a track by its proximity to the selected family line in ITS and/or TPC and/or TOF
/// - the same as above but additionally requiring a minimum separation to the other families lines
///
/// Both use cases can be configured on the same object instance. The parameters of the selected
/// family are interpreted as proximity cuts while the parameter of the other families are interpreted
/// as separation cuts. All above applicable within a configurable momentum range
///
/// Based on Gamma Conversion analysis implementation within PWGGA
///
/// \author Víctor González <victor.gonzalez@cern.ch>, UCM
/// \date Feb 09, 2017

class AliCSPIDCuts : public AliCSTrackCutsBase {
public:

  /// \enum cutsParametersIds
  /// \brief The ids of the different PID track cuts parameters supported
  enum cutsParametersIds {
    kPRangeCutParam,                 ///< momentum range cut parameters
    kITSdEdxSigmaCutParam_e,         ///< ITS dEdx \f$ n \sigma \f$ cut parameters for electron
    kITSdEdxSigmaCutParam_mu,        ///< ITS dEdx \f$ n \sigma \f$ cut parameters for muon
    kITSdEdxSigmaCutParam_pi,        ///< ITS dEdx \f$ n \sigma \f$ cut parameters for pion
    kITSdEdxSigmaCutParam_k,         ///< ITS dEdx \f$ n \sigma \f$ cut parameters for kaon
    kITSdEdxSigmaCutParam_p,         ///< ITS dEdx \f$ n \sigma \f$ cut parameters for proton
    kTPCdEdxSigmaCutParam_e,         ///< TPC dEdx \f$ n \sigma \f$ cut parameters for electron
    kTPCdEdxSigmaCutParam_mu,        ///< TPC dEdx \f$ n \sigma \f$ cut parameters for muon
    kTPCdEdxSigmaCutParam_pi,        ///< TPC dEdx \f$ n \sigma \f$ cut parameters for pion
    kTPCdEdxSigmaCutParam_k,         ///< TPC dEdx \f$ n \sigma \f$ cut parameters for kaon
    kTPCdEdxSigmaCutParam_p,         ///< TPC dEdx \f$ n \sigma \f$ cut parameters for proton
    kTOFSigmaCutParam_e,             ///< TOF \f$ n \sigma \f$ cut parameters for electron
    kTOFSigmaCutParam_mu,            ///< TOF \f$ n \sigma \f$ cut parameters for muon
    kTOFSigmaCutParam_pi,            ///< TOF \f$ n \sigma \f$ cut parameters for pion
    kTOFSigmaCutParam_k,             ///< TOF \f$ n \sigma \f$ cut parameters for kaon
    kTOFSigmaCutParam_p,             ///< TOF \f$ n \sigma \f$ cut parameters for proton
    kNCutsParameters                 ///< the number of track cuts parameters
  };

  /// \enum cutsIds
  /// \brief The ids of the different track cuts currently supported
  enum cutsIds {
    kITSdEdxSigmaCut,                ///< ITS dEdx \f$ n \sigma \f$ cut
    kTPCdEdxSigmaCut,                ///< TPC dEdx \f$ n \sigma \f$ cut
    kTOFSigmaCut,                    ///< TOF \f$ n \sigma \f$ cut
    kNCuts                           ///< The number of supported cuts
  };

public:

                                     AliCSPIDCuts();
                                     AliCSPIDCuts(const char *name, const char * title, AliPID::EParticleType target);
  virtual                           ~AliCSPIDCuts();

  virtual void                       InitCuts(const char *name = NULL);
  virtual void                       NotifyRun();
  virtual void                       NotifyEvent();
  virtual Bool_t                     IsTrackAccepted(AliVTrack *trk);
  virtual Bool_t                     IsTrueTrackAccepted(AliVTrack *trk);
  virtual Bool_t                     IsTrueTrackAccepted(Int_t itrk);

  static AliPID::EParticleType       GetTrueSpecies(AliVTrack *trk);
  static AliPID::EParticleType       GetTrueSpecies(AliVParticle *par);

private:
  void                               DefineHistograms();

private:
  /* we set them private to force cuts string consistency */
  Bool_t                             SetPRange(Int_t pcode);
  Bool_t                             SetITSdEdxSigmaCut(AliPID::EParticleType id, Int_t dEdxCode);
  Bool_t                             SetTPCdEdxSigmaCut(AliPID::EParticleType id, Int_t dEdxCode);
  Bool_t                             SetTOFSigmaCut(AliPID::EParticleType id, Int_t tofcode);
  void                               PrintITSdEdxSigmaCut(AliPID::EParticleType id) const;
  void                               PrintTPCdEdxSigmaCut(AliPID::EParticleType id) const;
  void                               PrintTOFSigmaCut(AliPID::EParticleType id) const;

private:
  virtual Bool_t                     SetCutAndParams(Int_t paramID, Int_t value);
  virtual void                       PrintCutWithParams(Int_t paramID) const;

private:
  static const char                 *fgkCutsParametersNames[kNCutsParameters];     ///< the names of the different event cuts parameters
  static const char                 *fgkCutsNames[kNCuts];                         ///< the names of the different event cuts

  Double_t                           fMinP;                                        ///< the minimum momentum value for applying the cut
  Double_t                           fMaxP;                                        ///< the maximum momentum value for applying the cut
  Double_t                           fITSnSigmaAbove[AliPID::kSPECIESC];           ///< \f$ n \sigma \f$ above a species ITS line cut
  Double_t                           fITSnSigmaBelow[AliPID::kSPECIESC];           ///< \f$ n \sigma \f$ below a species ITS line cut
  Double_t                           fTPCnSigmaAbove[AliPID::kSPECIESC];           ///< \f$ n \sigma \f$ above a species TPC line cut
  Double_t                           fTPCnSigmaBelow[AliPID::kSPECIESC];           ///< \f$ n \sigma \f$ below a species TPC line cut
  Bool_t                             fTOFRequired[AliPID::kSPECIESC];              ///< track presence in TOF required
  Double_t                           fTOFnSigmaAbove[AliPID::kSPECIESC];           ///< \f$ n \sigma \f$ above a species TOF line cut
  Double_t                           fTOFnSigmaBelow[AliPID::kSPECIESC];           ///< \f$ n \sigma \f$ below a species TOF line cut

  TBits                              fITSEnabledSpeciesMask;                       ///< mask for the species to consider within the ITS
  TBits                              fTPCEnabledSpeciesMask;                       ///< mask for the spdcies to consider within the TPC
  TBits                              fTOFEnabledSpeciesMask;                       ///< mask for the species to consider within the TOF

  AliPIDResponse                    *fPIDResponse;                                 ///< the PID response instance
  AliPID::EParticleType              fTargetSpecies;                               ///< the species to address with the PID cuts


  /* histograms with two indices: before (0) / after (1) applying the cut */
  TH1F                              *fhCutsStatistics;                             ///< the cuts statistics
  TH2F                              *fhCutsCorrelation;                            ///< cuts correlation
  TH2F                              *fhITSdEdxSigmaVsP[2];                         ///< ITS dE/dx \f$ \mbox{n} \sigma \f$ vs P (b/a)
  TH2F                              *fhITSdEdxSignalVsP[2];                        ///< ITS dE/dx signal vs P (b/a)
  TH2F                              *fhTPCdEdxSigmaVsP[2];                         ///< TPC dE/dx \f$ \mbox{n} \sigma \f$ vs P(b/a)
  TH2F                              *fhTPCdEdxSignalVsP[2];                        ///< TPC dE/dx signal vs P(b/a)
  TH2F                              *fhTOFSigmaVsP[2];                             ///< TOF \f$ \mbox{n} \sigma \f$ vs P(b/a)
  TH2F                              *fhTOFSignalVsP[2];                            ///< TOF signal vs P(b/a)

  /// Copy constructor
  /// Not allowed. Forced private.
  AliCSPIDCuts(const AliCSPIDCuts&);
  /// Assignment operator
  /// Not allowed. Forced private.
  /// \return l-value reference object
  AliCSPIDCuts& operator=(const AliCSPIDCuts&);

  /// \cond CLASSIMP
  ClassDef(AliCSPIDCuts,1);
  /// \endcond
};

#endif /* ALICSPIDCUTS_H */

