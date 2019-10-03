/**
 * \file AliReducedGeneratedParticle.h
 * \brief Declaration of class AliReducedGeneratedParticle
 *
 * In this header file a structure for reduced particle information at generator level
 * is declared.
 *
 * \author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 * \date Apr 14, 2015
 */
#ifndef ALIREDUCEDGENERATEDPARTICLE_H
#define ALIREDUCEDGENERATEDPARTICLE_H
/* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TObject.h>

class TVector3;

/**
 * \namespace HighPtTracks
 * \brief Namespace for classes creating trees of events with jets
 *
 * This namespace contains classes describing reduced events with high
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
 * \class AliReducedGeneratedParticle
 * \brief Structure for reduced particle information at generator level
 *
 * This structure stores basic information at generator level
 * - momentum vector of the particle
 * - energy
 * - PDG code
 * - Unique ID for matching with reconstructed track
 * This class is part of the reduced event for the high-\f$ p_{t} \f$ analysis
 */
class AliReducedGeneratedParticle: public TObject {
public:
  AliReducedGeneratedParticle();
  AliReducedGeneratedParticle(Int_t id, Int_t pdg, Double_t px, Double_t py, Double_t pz, Double_t energy);
  virtual ~AliReducedGeneratedParticle();

  /**
   * Get the unique ID of the particle
   * \return The unique ID of the particle
   */
  Int_t GetID() const { return fParticleID; }
  /**
   * Get the PDG code of the particle
   * \return The pdg code of the particle
   */
  Int_t GetPdgCode() const { return fPDGCode; }
  /**
   * Get the energy of the particle
   * \return The energy of the particle
   */
  Double_t GetParticleEnergy() const { return fEnergy; }

  Double_t Pt() const;
  Double_t Eta() const;
  Double_t Phi() const;

  void FillMomentumVector(TVector3 &) const;

  /**
   * Access to components of the  3-momentum vector.
   * \param px x-component
   * \param py y-component
   * \param pz z-component
   */
  void GetMomentumVector(Double_t &px, Double_t &py, Double_t &pz) const {
    px = fPVec[0];
    py = fPVec[1];
    pz = fPVec[2];
  }

  /**
   * Set the basic paricle parameters
   * \param id ID of the particle
   * \param pdgcode PDG code of the particle
   * \param px x-component of the momentum vector
   * \param py y-component of the momentum vector
   * \param pz z-component of the momentum vector
   * \param energy Particle energy
   */
  void Set(Int_t id, Int_t pdgcode, Double_t px, Double_t py, Double_t pz, Double_t energy){
    fParticleID = id;
    fPDGCode = pdgcode;
    fPVec[0] = px;
    fPVec[1] = py;
    fPVec[2] = px;
    fEnergy = energy;
  }
  /**
   * Set the unique ID of the particle
   * \param id Unique ID of the particle
   */
  void SetID(Int_t id) { fParticleID = id; }
  /**
   * Set the particle PDG code
   * \param pdgcode PDG code of the particle
   */
  void SetPDGCode(Int_t pdgcode) { fPDGCode = pdgcode; }
  /**
   * Set the particle 3-momentum vector
   * \param px x-component of the momentum vector
   * \param py y-component of the momentum vector
   * \param pz z-component of the momentum vector
   */
  void SetPVector(Double_t px, Double_t py, Double_t pz){
    fPVec[0] = px;
    fPVec[1] = py;
    fPVec[2] = px;
  }
  /**
   * Set the particle energy
   * \param energy The particle energy
   */
  void SetEnergy(Double_t energy) { fEnergy = energy; }

protected:
  Int_t                       fParticleID;          ///< Unique ID of the particle
  Int_t                       fPDGCode;             ///< PDG code
  Double_t                    fPVec[3];             ///< Particle momentum vector
  Double_t                    fEnergy;              ///< Particle energy

  /// \cond CLASSIMP
  ClassDef(AliReducedGeneratedParticle, 1);
  /// \endcond
};

} /* namespace HighPtTracks */

#endif /* ALIREDUCEDGENERATEDPARTICLE_H */
