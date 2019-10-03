/**
 * \file AliReducedEmcalCluster.h
 * \brief Reduced information about reconstructed EMCAL clusters
 *
 * Declaration of class AliReducedEmcalCluster, a structure for reduced information of reconstructed EMCAL clusters
 *
 * \author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 * \date Apr 14, 2015
 */
#ifndef ALIREDUCEDEMCALCLUSTER_H
#define ALIREDUCEDEMCALCLUSTER_H
/* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TObject.h>

class TArrayD;
class TLorentzVector;
class TObjArray;

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
 * \class AliReducedClusterParticle
 * \brief MC true contributor to a reconstructed EMCAL cluster
 *
 * This class is part of the reduced EMCAL cluster and stores true information of particles
 * matched to the EMCAL cluster based on the label
 * - 4-momentum vector
 * - PDG code
 */
class AliReducedClusterParticle : public TObject {
public:
  AliReducedClusterParticle();
  AliReducedClusterParticle(Int_t pdg, Double_t px, Double_t py, Double_t pz, Double_t energy);
  virtual ~AliReducedClusterParticle();

  void FillLorentzVector(TLorentzVector &target) const;
  /**
   * Get the PDG code of the particle
   * \return PDG code of the particle
   */
  Int_t GetPdgCode() const { return fPdg; }

  /**
   * Set the particle momentum
   * \param px x-component of the particle momentum
   * \param py y-component of the particle momentum
   * \param pz z-component of the particle momentum
   */
  void SetMomentum(Double_t px, Double_t py, Double_t pz){
    fPvec[0]= px;
    fPvec[1]= py;
    fPvec[2]= pz;
  }
  /**
   * Set the particle energy
   * \param energy The particle energy
   */
  void SetEnergy(Double_t energy) { fEnergy = energy; }
  /**
   * Set the particle PDG code
   * \param The PDG code of the particle
   */
  void SetPdgCode(Int_t pdg) { fPdg = pdg; }

protected:
  Int_t               fPdg;               ///< Particle PDG code
  Double_t            fPvec[3];           ///< Particle momentum vector
  Double_t            fEnergy;            ///< Particle energy

  /// \cond CLASSIMP
  ClassDef(AliReducedClusterParticle, 1);
  /// \endcond
};

/**
 * \class AliReducedEmcalCluster
 * \brief Reduced EMCAL cluster information
 *
 * This object stores a reduced set of information of a reconstructed EMCAL cluster
 * - unique ID of a cluster in the event (for cluster-track matching)
 * - energy of the cluster
 * - position of the cluster in (\f$ \eta \f$, \f$ \phi \f$)
 * - shower shape parameters
 * This object is part of the reduced event structure for the high-\f$ p_{t} \f$ analysis
 */
class AliReducedEmcalCluster: public TObject {
public:
  AliReducedEmcalCluster();
  AliReducedEmcalCluster(Int_t id, Float_t energy, Float_t eta, Float_t phi, Float_t m02, Float_t m20);
  AliReducedEmcalCluster(const AliReducedEmcalCluster &ref);
  AliReducedEmcalCluster &operator=(const AliReducedEmcalCluster &ref);
  virtual ~AliReducedEmcalCluster();
  void Copy(TObject &target) const;

  /**
   * Get the ID of the cluster inside the event
   * \return ID of the cluster
   */
  Int_t           GetClusterID() const        { return fClusterID; }
  /**
   * Get the cluster energy
   * \return Energy of the cluster
   */
  Float_t         GetClusterEnergy() const    { return fEnergy; }
  /**
   * Get the \f$ \eta \f$ position of the cluster
   * \return \f$ \eta \f$ position of the cluster relative to the primary vertex
   */
  Float_t         GetEta() const              { return fEta; }
  /**
   * Get the \f$ \phi \f$ position of the cluster
   * \return \f$ \phi \f$ position of the cluster relative to the primary vertex
   */
  Float_t         GetPhi() const              { return fPhi; }
  /**
   * Get the M02 shower shape parameter
   * \return the M02 shower shape parameter
   */
  Float_t         GetM02() const              { return fM02; }
  /**
   * Get the M20 shower shape parameter
   * \return the M20 shower shape parameter
   */
  Float_t         GetM20() const              { return fM20; }
  void            FillCellEnergies(TArrayD &target);
  /**
   * Get MC-true contributing particles
   * \return Array of contributing particles (NULL if not available)
   */
  TObjArray      *GetClusterContributors() const { return fContributors; }

  /**
   * Set the cluster ID
   * \param id Cluster ID
   */
  void            SetClusterID(Int_t id)                                          { fClusterID = id; }
  /**
   * Set the cluster energy
   * \param energy
   */
  void            SetClusterEnergy(Float_t energy)                                { fEnergy = energy; }
  /**
   * Set the cluster position in \f$ \eta \f$ and \f$ \phi \f$
   * \param eta \f$ \eta \f$ position of the cluster
   * \param phi \f$ \phi \f$ position of the cluster
   */
  void            SetClusterPosition(Float_t eta, Float_t phi)                    { fEta = eta; fPhi = phi; }
  /**
   * Set the cluster shower shape parameters
   * \param m02 The m02 shower shape parameter
   * \param m20 The m20 shower shape parameter
   */
  void            SetShowerShapeParameters(Float_t m02, Float_t m20)              { fM02 = m02; fM20 = m20; }
  /**
   * Set the leading cell energies
   * \param e1 max cell energy
   * \param e2 2nd highest cell energy
   * \param e3 3rd highest cell energy
   */
  void            SetLeadingCellEnergies(Double_t e1, Double_t e2, Double_t e3)   {
    fCellEnergies[0] = e1;
    fCellEnergies[1] = e2;
    fCellEnergies[2] = e3;
  }
  /**
   * Set cluster parameters
   * \param id ID of the cluster
   * \param energy Energy of the cluster
   * \param eta \f$ \eta \f$ position of the cluster
   * \param phi \f$ \phi \f$ position of the cluster
   * \param m02 The m02 shower shape parameter
   * \param m20 The m20 shower shape parameter
   */
  void  Set(Int_t id, Float_t energy, Float_t eta, Float_t phi, Float_t m02, Float_t m20){
    fClusterID = id;
    fEnergy = energy;
    fEta = eta;
    fPhi = phi;
    fM02 = m02;
    fM20 = m20;
  }
  void AddTrueContributor(Int_t pdg, Double_t px, Double_t py, Double_t pz, Double_t energy);
protected:
  Int_t                     fClusterID;             ///< ID of the cluster
  Float_t                   fEnergy;                ///< Energy of the cluster
  Float_t                   fEta;                   ///< Cluster position in \f$ \eta \f$ relative to the primary vertex
  Float_t                   fPhi;                   ///< Cluster position in \f$ \phi \f$ relative to the primary vertex
  Float_t                   fM02;                   ///< M02 shower shape parameter
  Float_t                   fM20;                   ///< M20 shower shape parameter
  Float_t                   fCellEnergies[3];       ///< Leading cell energies
  TObjArray                 *fContributors;         ///< True particles contributing to the cluster

  /// \cond CLASSIMP
  ClassDef(AliReducedEmcalCluster, 1);
  /// \endcond
};

} /* namespace HighPtTracks */

#endif /* ALIREDUCEDEMCALCLUSTER_H */
