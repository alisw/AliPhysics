#ifndef ALICSPIDCUTS_H
#define ALICSPIDCUTS_H

/// \file AliCSPIDCuts.h
/// \brief PID track selection cuts support for the correlations studies tasks
///
/// Based on ideas taken from Gamma Conversion analysis under PWGGA and from AliESDTrackCuts

#include <TBits.h>
#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliPID.h"
#include "AliCSTrackMaps.h"
#include "AliCSTrackCutsBase.h"
#include "AliPIDResponse.h"
#include "AliVTrack.h"

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
    kPMinCutParam,            ///< minimum momentum cut parameters
    kPMaxCutParam,            ///< maximum momentum cut parameters
    kITSdEdxSigmaCutParam_e,  ///< ITS dEdx \f$ n \sigma \f$ cut parameters for electron
    kITSdEdxSigmaCutParam_mu, ///< ITS dEdx \f$ n \sigma \f$ cut parameters for muon
    kITSdEdxSigmaCutParam_pi, ///< ITS dEdx \f$ n \sigma \f$ cut parameters for pion
    kITSdEdxSigmaCutParam_k,  ///< ITS dEdx \f$ n \sigma \f$ cut parameters for kaon
    kITSdEdxSigmaCutParam_p,  ///< ITS dEdx \f$ n \sigma \f$ cut parameters for proton
    kTPCdEdxSigmaCutParam_e,  ///< TPC dEdx \f$ n \sigma \f$ cut parameters for electron
    kTPCdEdxSigmaCutParam_mu, ///< TPC dEdx \f$ n \sigma \f$ cut parameters for muon
    kTPCdEdxSigmaCutParam_pi, ///< TPC dEdx \f$ n \sigma \f$ cut parameters for pion
    kTPCdEdxSigmaCutParam_k,  ///< TPC dEdx \f$ n \sigma \f$ cut parameters for kaon
    kTPCdEdxSigmaCutParam_p,  ///< TPC dEdx \f$ n \sigma \f$ cut parameters for proton
    kTOFSigmaCutParam_e,      ///< TOF \f$ n \sigma \f$ cut parameters for electron
    kTOFSigmaCutParam_mu,     ///< TOF \f$ n \sigma \f$ cut parameters for muon
    kTOFSigmaCutParam_pi,     ///< TOF \f$ n \sigma \f$ cut parameters for pion
    kTOFSigmaCutParam_k,      ///< TOF \f$ n \sigma \f$ cut parameters for kaon
    kTOFSigmaCutParam_p,      ///< TOF \f$ n \sigma \f$ cut parameters for proton
    kTPCTOFCutParam,          ///< TOF required and TPC+TOF 2D cut
    kNCutsParameters          ///< the number of track cuts parameters
  };

  /// \enum cutsIds
  /// \brief The ids of the different track cuts currently supported
  enum cutsIds {
    kMom,              ///< momentum cut
    kITSdEdxSigmaCut,  ///< ITS dEdx \f$ n \sigma \f$ cut
    kTPCdEdxSigmaCut,  ///< TPC dEdx \f$ n \sigma \f$ cut
    kTOFSigmaCut,      ///< TOF \f$ n \sigma \f$ cut
    kTPCTOF2DSigmaCut, ///< 2D TPC+TOF \f$ n \sigma \f$ cut
    kNCuts             ///< The number of supported cuts
  };

  /// \typedef \enum detectors
  /// \brief The ids of the detectors used for PID
  typedef enum {
    kITS,
    kTPC,
    kTOF,
    kTPCTOF2D
  } detectors;

 public:
  AliCSPIDCuts();
  AliCSPIDCuts(const char* name, const char* title, AliPID::EParticleType target, int ncut);
  virtual ~AliCSPIDCuts();

  virtual void InitCuts(const char* name = NULL);
  virtual void NotifyRun();
  virtual void NotifyEvent();
  AliPID::EParticleType GetTargetSpecies() { return fTargetSpecies; }
  virtual Bool_t IsTrackAccepted(AliVTrack* trk, float*);
  virtual Bool_t IsTrueTrackAccepted(AliVTrack* trk);
  virtual Bool_t IsTrueTrackAccepted(Int_t itrk);

  static AliPID::EParticleType GetTrueSpecies(AliVTrack* trk);
  static AliPID::EParticleType GetTrueSpecies(AliVParticle* par);

 private:
  void DefineHistograms();

 public:
  template <detectors det>
  float nsigmas(AliPID::EParticleType s, AliVTrack* t);

 protected:
  template <detectors det>
  bool nearerIn(AliPID::EParticleType s, AliVTrack* t);
  template <detectors det>
  bool fartherIn(AliVTrack* t);
  template <detectors det>
  bool check(AliVTrack* trk);
  bool accept(AliVTrack* ttrk);

 private:
  /* we set them private to force cuts string consistency */
  bool SetPMin(Int_t pcode);
  bool SetPMax(Int_t pcode);
  Bool_t SetITSdEdxSigmaCut(AliPID::EParticleType id, Int_t dEdxCode);
  Bool_t SetTPCdEdxSigmaCut(AliPID::EParticleType id, Int_t dEdxCode);
  Bool_t SetTOFSigmaCut(AliPID::EParticleType id, Int_t tofcode);
  bool SetTPCTOFCut(int tpctofcode);
  void PrintITSdEdxSigmaCut(AliPID::EParticleType id) const;
  void PrintTPCdEdxSigmaCut(AliPID::EParticleType id) const;
  void PrintTOFSigmaCut(AliPID::EParticleType id) const;

 private:
  virtual Bool_t SetCutAndParams(Int_t paramID, Int_t value);
  virtual void PrintCutWithParams(Int_t paramID) const;

 private:
  static const char* fgkCutsParametersNames[kNCutsParameters]; ///< the names of the different event cuts parameters
  static const char* fgkCutsNames[kNCuts];                     ///< the names of the different event cuts

 private:
  float fMinP;                              ///< the minimum momentum value for applying the cut
  float fMaxP;                              ///< the maximum momentum value for applying the cut
  bool fTOFRequired;                        ///< track presence in TOF required
  bool fTPCTOF2Dcut;                        ///< bidimensional TPC+TOF cut
  float fITSnSigmaAbove[AliPID::kSPECIESC]; ///< \f$ n \sigma \f$ above a species ITS line cut
  float fITSnSigmaBelow[AliPID::kSPECIESC]; ///< \f$ n \sigma \f$ below a species ITS line cut
  float fTPCnSigmaAbove[AliPID::kSPECIESC]; ///< \f$ n \sigma \f$ above a species TPC line cut
  float fTPCnSigmaBelow[AliPID::kSPECIESC]; ///< \f$ n \sigma \f$ below a species TPC line cut
  float fTOFnSigmaAbove[AliPID::kSPECIESC]; ///< \f$ n \sigma \f$ above a species TOF line cut
  float fTOFnSigmaBelow[AliPID::kSPECIESC]; ///< \f$ n \sigma \f$ below a species TOF line cut

  TBits fCutsEnabledMask;   ///< the mask of enabled cuts
  TBits fCutsActivatedMask; ///< the mask of cut activated for the ongoing event

  TBits fITSEnabledSpeciesMask;      ///< mask for the species to consider within the ITS
  TBits fTPCEnabledSpeciesMask;      ///< mask for the spdcies to consider within the TPC
  TBits fTOFEnabledSpeciesMask;      ///< mask for the species to consider within the TOF
  TBits fTPCTOF2DEnabledSpeciesMask; ///< mask for the species to consider within the 2D TPC+TOF

  AliPIDResponse* fPIDResponse;         ///< the PID response instance
  AliPID::EParticleType fTargetSpecies; ///< the species to address with the PID cuts
  int fCutNumber;                       ///< the sequential cut number for this species

  /* histograms with two indices: before (0) / after (1) applying the cut */
  TH1F* fhCutsStatistics;      ///< the cuts statistics
  TH2F* fhCutsCorrelation;     ///< cuts correlation
  TH2F* fhITSdEdxSigmaVsP[2];  ///< ITS dE/dx \f$ \mbox{n} \sigma \f$ vs P (b/a)
  TH2F* fhITSdEdxSignalVsP[2]; ///< ITS dE/dx signal vs P (b/a)
  TH2F* fhTPCdEdxSigmaVsP[2];  ///< TPC dE/dx \f$ \mbox{n} \sigma \f$ vs P(b/a)
  TH2F* fhTPCdEdxSignalVsP[2]; ///< TPC dE/dx signal vs P(b/a)
  TH2F* fhTOFSigmaVsP[2];      ///< TOF \f$ \mbox{n} \sigma \f$ vs P(b/a)
  TH2F* fhTOFSignalVsP[2];     ///< TOF signal vs P(b/a)
  TH2F* fhTPCTOFSigma[2];      ///< TPC+TOF 2D \f$ \mbox{n} \sigma \f$

  /// Copy constructor
  /// Not allowed. Forced private.
  AliCSPIDCuts(const AliCSPIDCuts&);
  /// Assignment operator
  /// Not allowed. Forced private.
  /// \return l-value reference object
  AliCSPIDCuts& operator=(const AliCSPIDCuts&);

  /// \cond CLASSIMP
  ClassDef(AliCSPIDCuts, 2);
  /// \endcond
};

template <AliCSPIDCuts::detectors det>
inline float AliCSPIDCuts::nsigmas(AliPID::EParticleType s, AliVTrack* t)
{
  if constexpr (det == kITS) {
    return fPIDResponse->NumberOfSigmasITS(t, s);
  }
  if constexpr (det == kTPC) {
    return fPIDResponse->NumberOfSigmasTPC(t, s);
  }
  if constexpr (det == kTOF) {
    return fPIDResponse->NumberOfSigmasTOF(t, s);
  }
  ::Fatal("AliCSPIDCuts::nsigmas", "Unrecoginzed detector id %d ", det);
  return -999.9;
}

template <AliCSPIDCuts::detectors det>
inline bool AliCSPIDCuts::nearerIn(AliPID::EParticleType s, AliVTrack* t)
{
  float val;
  float lowlim;
  float uplim;

  auto nearer = [](float value, float below, float above) {
    if (value < below || above < value) {
      return false;
    } else {
      return true;
    }
  };

  if constexpr (det == kITS) {
    val = nsigmas<det>(s, t);
    lowlim = fITSnSigmaBelow[s];
    uplim = fITSnSigmaAbove[s];
  }
  if constexpr (det == kTPC) {
    val = nsigmas<det>(s, t);
    lowlim = fTPCnSigmaBelow[s];
    uplim = fTPCnSigmaAbove[s];
  }
  if constexpr (det == kTOF) {
    val = nsigmas<det>(s, t);
    lowlim = fTOFnSigmaBelow[s];
    uplim = fTOFnSigmaAbove[s];
  }
  if constexpr (det == kTPCTOF2D) {
    /* special handling for 2D TPC+TOF */
    float vtpc = nsigmas<kTPC>(s, t);
    float vtof = nsigmas<kTOF>(s, t);
    val = vtpc * vtpc / fTPCnSigmaAbove[s] / fTPCnSigmaAbove[s] + vtof * vtof / fTOFnSigmaAbove[s] / fTOFnSigmaAbove[s];
    return nearer(val, 0.0, 1.0);
  }
  if constexpr (det != kITS && det != kTPC && det != kTOF && det != kTPCTOF2D) {
    ::Fatal("AliCSPIDCuts", "Unrecoginzed detector id %d", det);
    return false;
  }
  return nearer(val, lowlim, uplim);
}

template <AliCSPIDCuts::detectors det>
inline bool AliCSPIDCuts::fartherIn(AliVTrack* t)
{
  const auto& enabled = [&]() {
    if constexpr (det == kITS) {
      return fITSEnabledSpeciesMask;
    }
    if constexpr (det == kTPC) {
      return fTPCEnabledSpeciesMask;
    }
    if constexpr (det == kTOF) {
      return fTOFEnabledSpeciesMask;
    }
    if constexpr (det == kTPCTOF2D) {
      return fTPCTOF2DEnabledSpeciesMask;
    }
    ::Fatal("AliCSPIDCuts", "Unrecoginzed detector id %d", det);
    return fTPCEnabledSpeciesMask;
  };

  const TBits& mask = enabled();
  if (mask.CountBits() > 1) {
    AliPID::EParticleType spid = AliPID::EParticleType(mask.FirstSetBit());
    while (spid < AliPID::kDeuteron) {
      if (spid != fTargetSpecies) {
        if (nearerIn<det>(spid, t)) {
          return false;
        }
      }
      spid = AliPID::EParticleType(mask.FirstSetBit(spid + 1));
    }
  }
  return true;
}

template <AliCSPIDCuts::detectors det>
inline bool AliCSPIDCuts::check(AliVTrack* trk)
{
  auto id = []() {
    if constexpr (det == kITS) {
      return kITSdEdxSigmaCut;
    }
    if constexpr (det == kTPC) {
      return kTPCdEdxSigmaCut;
    }
    if constexpr (det == kTOF) {
      return kTOFSigmaCut;
    }
    if constexpr (det == kTPCTOF2D) {
      return kTPCTOF2DSigmaCut;
    }
    ::Fatal("AliCSPIDCuts", "Unrecoginzed detector id %d", det);
    return kTPCdEdxSigmaCut;
  };
  AliCSPIDCuts::cutsIds cutid = id();
  if (fCutsEnabledMask.TestBitNumber(cutid)) {
    if (!nearerIn<det>(fTargetSpecies, trk)) {
      fCutsActivatedMask.SetBitNumber(cutid);
      return false;
    } else {
      /* proximity accepted, now check the separation if required */
      if (!fartherIn<det>(trk)) {
        fCutsActivatedMask.SetBitNumber(cutid);
        return false;
      }
    }
  }
  return true;
}

inline bool AliCSPIDCuts::accept(AliVTrack* ttrk)
{
  /* sanity check */
  if (ttrk == nullptr) {
    return false;
  }

  /* initialize the mask of activated cuts */
  fCutsActivatedMask.ResetAllBits();

  /* if not in the momentum range it is not accepted */
  if (ttrk->P() < fMinP || fMaxP < ttrk->P()) {
    fCutsActivatedMask.SetBitNumber(kMom);
    return false;
  }

  /* for the time being */
  bool accepted = true;

  /* consider a potential constrained track */
  AliVTrack* trk = ttrk;
  if (trk->GetID() < 0) {
    /* let's switch to the original one which has the PID information */
    trk = AliCSTrackMaps::GetOriginalTrack(dynamic_cast<AliAODTrack*>(ttrk));
  }

  bool bits = check<kITS>(trk);
  bool btpc = check<kTPC>(trk);
  /* special handling for TOF */
  bool btof = true;
  if ((trk->GetStatus() & AliESDtrack::kTOFin) && (!(trk->GetStatus() & AliESDtrack::kTOFmismatch))) {
    if (fTPCTOF2Dcut) {
      btof = btof && check<kTPCTOF2D>(trk);
      btpc = btof;
    } else {
      btof = btof && check<kTOF>(trk);
    }
  } else if (fTOFRequired) {
    /* TOF is required but is not there */
    btof = false;
  }

  accepted = bits && btpc && btof;
  return accepted;
}

#endif /* ALICSPIDCUTS_H */
