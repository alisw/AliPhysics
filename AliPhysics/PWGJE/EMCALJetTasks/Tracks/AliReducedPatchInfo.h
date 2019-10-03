/**
 * \file AliReducedPatchInfo.h
 * \brief Declaration of class AlliReducedPatchInfo
 *
 * In this header file, the class AliReducedPatchInfo is declared. This object as part of the high-\f$ p_{t} \f$ event
 * structure contains basic information about EMCAL trigger patches.
 *
 * \author Markus Fasel <markus.fasel\/gcern.ch>, Lawrence Berkeley National Laboratory
 * \date Apr 14, 2015
 */
#ifndef ALIREDUCEDPATCHINFO_H
#define ALIREDUCEDPATCHINFO_H
/* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TObject.h>

/**
 * \namespace HighPtTracks
 * \brief Namespace for classes creating trees of events with jets
 *
 * This namespace contains classes descibing reduced events with high
 * jets. A jet event consists of the following classes.
 *    - AliReducedJetEvent
 *    - AliReducedJetParticle
 *    - AliReducedJetConstituent
 *    - AliReducedMatchedTrack
 *  The task AliHighPtReconstructionEfficiency produces the reduced jet events
 *  and is included in the namespace as well. Also the helper class AliParticleMap
 *  is part of this namespace.
 */
namespace HighPtTracks {

/**
 * \class AliReducedPatchInfo
 * \brief Reduced information of EMCAL trigger patches
 *
 * This class stores reduced information about EMCAL trigger patches. Information available are:
 * - Offline energy
 * - Online amplitude
 * - patch center position
 * This class is part of the reduced high-\f$ p_{t} \f$ event structure
 */
class AliReducedPatchInfo : public TObject {
public:
  AliReducedPatchInfo();
  AliReducedPatchInfo(Float_t energy, Float_t amplitude, Float_t eta, Float_t phi);
  /**
   * Destructor
   */
  virtual ~AliReducedPatchInfo() {}

  Float_t GetEnergy() const { return fEnergy; }
  Int_t GetAmplitude() const { return fAmplitude; }
  Float_t GetEta() const { return fEta; }
  Float_t GetPhi() const { return fPhi; }

  /**
   * Set Patch information
   * \param energy Patch energy
   * \param amplitude Patch amplitude
   * \param eta Patch \f$ \eta \f$ (center)
   * \param phi Patch \f$ \phi \f$ (center)
   */
  void Set(Float_t energy, Float_t amplitude, Float_t eta, Float_t phi){
    fEnergy = energy;
    fAmplitude = amplitude;
    fEta = eta;
    fPhi = phi;
  }
  /**
   * Set the patch energy
   * \param energy Patch energy
   */
  void SetEnergy(Float_t energy){ fEnergy = energy; }
  /**
   * Set the patch \f$ \eta \f$ position at the patch center
   * \param eta Patch \f$ \eta \f$
   */
  void SetEta(Float_t eta) { fEta = eta; }
  /**
   * Set the patch \f$ \phi \f$ at the patch center
   * \param phi Patch \f$ \eta \f$
   */
  void SetPhi(Float_t phi){ fPhi = phi; }
  /**
   * Set the patch amplitude
   * \param amplitude Patch amplitude
   */
  void SetAmplitude(Int_t amplitude) { fAmplitude = amplitude; }
protected:
  Float_t             fEnergy;              ///< Patch energy
  Int_t               fAmplitude;           ///< Patch amplitude
  Float_t             fEta;                 ///< Patch \f$ \eta \f$ (center)
  Float_t             fPhi;                 ///< Patch \f$ \phi \f$ (center)

  /// \cond CLASSIMP
  ClassDef(AliReducedPatchInfo, 1);
  /// \endcond
};

} /* namespace HighPtTracks */

#endif /* ALIREDUCEDPATCHINFO_H */
