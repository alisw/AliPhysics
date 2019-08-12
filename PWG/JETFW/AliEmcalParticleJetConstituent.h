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
#ifndef ALIEMCALPARTICLEJETCONSTITUENT_H
#define ALIEMCALPARTICLEJETCONSTITUENT_H

#include "AliEmcalJetConstituent.h"

class AliVParticle;

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

class AliEmcalParticleJetConstituent : public AliEmcalJetConstituent {
public:
  AliEmcalParticleJetConstituent();

  /**
   * @brief Main constructor
   *
   * Fully initializing the particle constituent with the particle
   * @param[in] part Particle trajectory
   */
  AliEmcalParticleJetConstituent(const AliVParticle *const part);

  /**
   * @brief Copy constructor
   * Constituent object is non-owning
   * @param other Object from which to copy
   */
  AliEmcalParticleJetConstituent(const AliEmcalParticleJetConstituent &other);

  /**
   * @brief Assignment operator
   * Constituent object is non-owning, therefore shallow copy is sufficient
   * @param[in] other
   * @return Particle constituent object after assignment
   */
  AliEmcalParticleJetConstituent &operator=(const AliEmcalParticleJetConstituent &other);

  /**
   * @brief Check for equality
   * Check done on the cluster itself
   * @param[in] rhs Object to be checked against
   * @return True if underlying clusters are equal, false otherwise
   */
  bool operator==(const AliEmcalParticleJetConstituent &rhs) const;

  /**
   * @brief Comparison operator (implemented for pT)
   * @param[in] rhs Object to be checked against
   * @return True if constituent pT is larger, false otherwise
   */
  bool operator>(const AliEmcalParticleJetConstituent &rhs) const;

  /**
   * @brief Comparison operator (implemented for pT)
   * @param[in] rhs Object to be checked against
   * @return True if constituent pT is smaller, false otherwise
   */
  bool operator<(const AliEmcalParticleJetConstituent &rhs) const;

  /**
   * @brief Destructor
   */
  virtual ~AliEmcalParticleJetConstituent();

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

  const AliVParticle *GetParticle() const { return fkParticle;}

private:
  const AliVParticle *fkParticle;         ///< Underlying particle


  /// \cond CLASSIMP
  ClassDef(AliEmcalParticleJetConstituent, 0);
  /// \endcond
};

}
}

#endif
