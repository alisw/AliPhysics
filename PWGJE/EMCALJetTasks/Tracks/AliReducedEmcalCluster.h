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
  virtual ~AliReducedEmcalCluster();

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
protected:
  Int_t                     fClusterID;             ///< ID of the cluster
  Float_t                   fEnergy;                ///< Energy of the cluster
  Float_t                   fEta;                   ///< Cluster position in \f$ \eta \f$ relative to the primary vertex
  Float_t                   fPhi;                   ///< Cluster position in \f$ \phi \f$ relative to the primary vertex
  Float_t                   fM02;                   ///< M02 shower shape parameter
  Float_t                   fM20;                   ///< M20 shower shape parameter

  /// \cond CLASSIMP
  ClassDef(AliReducedEmcalCluster, 1);
  /// \endcond
};

} /* namespace HighPtTracks */

#endif /* ALIREDUCEDEMCALCLUSTER_H */
