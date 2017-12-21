/************************************************************************************
 * Copyright (C) 2017, Copyright Holders of the ALICE Collaboration                 *
 * All rights reserved.                                                             *
 *                                                                                  *
 * Redistribution and use in source and binary forms, with or without               *
 * modification, are permitted provided that the following conditions are met:      *
 *     * Redistributions of source code must retain the above copyright             *
 *       notice, this list of conditions and the following disclaimer.              *
 *     * Redistributions in binary form must reproduce the above copyright          *
 *       notice, this list of conditions and the following disclaimer in the        *
 *       documentation and/or other materials provided with the distribution.       *
 *     * Neither the name of the <organization> nor the                             *
 *       names of its contributors may be used to endorse or promote products       *
 *       derived from this software without specific prior written permission.      *
 *                                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND  *
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED    *
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE           *
 * DISCLAIMED. IN NO EVENT SHALL ALICE COLLABORATION BE LIABLE FOR ANY              *
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES       *
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;     *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND      *
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT       *
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS    *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                     *
 ************************************************************************************/
#ifndef ALIEMCALCLUSTERJETCONSTITUENT_H_
#define ALIEMCALCLUSTERJETCONSTITUENT_H_

#include "AliEmcalJetConstituent.h"
#include "AliVCluster.h"

/**
 * @namespace PWG
 * @brief Basic namespace for general framework objects
 */
namespace PWG {

/**
 * @namespace JETFW
 * @brief Namespace for objects belonging to the ALICE jet framework
 * @ingroup JETFW
 */
namespace JETFW {

/**
 * @class AliEmcalClusterJetConstituent
 * @brief Implementation of a jet constituent for constituent clusters
 * @ingroup JETFW
 * @author Markus Fasel <markus.fasel@cern.ch>, Oak Ridge National Laboratory
 * @since Sept 22, 2017
 *
 * Cluster constituents provide access to the kinematic quantities as used in the
 * jetfinder. Furthermore it handles settings used in the cluster selection in the
 * jetfinder like default cluster energy, and provides access to the underlying cluster
 * itself.
 */
class AliEmcalClusterJetConstituent : public AliEmcalJetConstituent {
public:
  /**
   * @brief Default constructor
   */
  AliEmcalClusterJetConstituent();

  /**
   * @brief Main constructor
   *
   * Fully initializing the cluster consituent with the EMCAL cluster and the momentum vector used in the jetfinder
   * @param[in] clust EMCAL cluster
   * @param[in] energydef Cluster
   * @param[in] pvec Cluster momentum vector
   */
  AliEmcalClusterJetConstituent(const AliVCluster *const clust, AliVCluster::VCluUserDefEnergy_t energydef, Double_t *pvec);

  /**
   * @brief Cluster-only constructor
   *
   * Only initializing cluster pointer. This constructer is intended for the usage in the comparator
   * @param[in] clust EMCAL cluster
   */
  AliEmcalClusterJetConstituent(const AliVCluster *const clust);

  /**
   * @brief Copy constructor
   * Constituent object is non-owning
   * @param other Object from which to copy
   */
  AliEmcalClusterJetConstituent(const AliEmcalClusterJetConstituent &other);

  /**
   * @brief Assignment operator
   * Constituent object is non-owning, therefore shallow copy is sufficient
   * @param[in] other
   * @return Cluster constituent object after assignment
   */
  AliEmcalClusterJetConstituent &operator=(const AliEmcalClusterJetConstituent &other);

  /**
   * @brief Check for equality
   * Check done on the cluster itself
   * @param[in] rhs Object to be checked against
   * @return True if underlying clusters are equal, false otherwise
   */
  bool operator==(const AliEmcalClusterJetConstituent &rhs) const;

  /**
   * @brief Destructor
   */
  virtual ~AliEmcalClusterJetConstituent();

  /**
   * @brief Access to x-component of the momentum vector
   * @return x-component of the momentum vector
   */
  virtual double Px() const;

  /**
   * @brief Access to y-component of the momentum vector
   * @return y-component of the momentum vector
   */
  virtual double Py() const;

  /**
   * @brief Access to z-component of the momentum vector
   * @return z-component of the momentum vector
   */
  virtual double Pz() const;

  /**
   * @brief Access to transverse momentum
   * @return Constituent transverse momentum
   */
  virtual double Pt() const;

  /**
   * @brief Access to constituent energy
   * @return Constituent energy
   */
  virtual double E() const;

  /**
   * @brief Access to pseudorapidity
   * @return Constituent pseudorapidity
   */
  virtual double Eta() const;

  /**
   * @brief Access to azimuthal angle
   * @return Constituent azimuthal angle
   */
  virtual double Phi() const;

  /**
   * @brief Get the underlying cluster
   * @return
   */
  const AliVCluster *GetCluster() const { return fkCaloCluster; }

  /**
   * @brief Get the energy definition used to calculate the cluster energy in the jetfinder
   * @return Cluster energy definition
   */
  AliVCluster::VCluUserDefEnergy_t GetDefaultEnergyType() const { return fDefaultEnergyDefinition; }

private:
  const AliVCluster                     *fkCaloCluster;              ///< Underlying calorimeter cluster
  AliVCluster::VCluUserDefEnergy_t      fDefaultEnergyDefinition;   ///< Default energy definition used in jetfinder
  TVector3                              fPVec;                      ///< Momentum vector

  /// \cond CLASSIMP
  ClassDef(AliEmcalClusterJetConstituent, 1);
  /// \endcond
};

}

}

#endif
