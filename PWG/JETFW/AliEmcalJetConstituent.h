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
#ifndef ALIEMCALJETCONSTITUENT_H
#define ALIEMCALJETCONSTITUENT_H

#include <TObject.h>

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
 * @class AliEmcalJetConstituent
 * @brief Interface class for constituent objects (clusters / particles) in an ALICE jet
 * @ingroup JETFW
 * @author Markus Fasel <markus.fasel@cern.ch>, Oak Ridge National Laboratory
 * @since Sept 22, 2017
 *
 * The interface object defines access functions to kinematic quantities (4-momentum vector).
 * Implementation is done within the classes implementing the constituents for particles and
 * clusters. Here also the access to the underlying object is provided.
 */
class AliEmcalJetConstituent : public TObject {
public:
  /**
   * @brief Constructor
   */
  AliEmcalJetConstituent();

  /**
   * @brief Destructor
   */
  virtual ~AliEmcalJetConstituent();

  /**
   * @brief Access to x-component of the momentum vector
   * @return x-component of the momentum vector
   */
  virtual double Px() const = 0;

  /**
   * @brief Access to y-component of the momentum vector
   * @return y-component of the momentum vector
   */
  virtual double Py() const = 0;

  /**
   * @brief Access to z-component of the momentum vector
   * @return z-component of the momentum vector
   */
  virtual double Pz() const = 0;

  /**
   * @brief Access to transverse momentum
   * @return Constituent transverse momentum
   */
  virtual double Pt() const = 0;

  /**
   * @brief Access to constituent energy
   * @return Constituent energy
   */
  virtual double E() const = 0;

  /**
   * @brief Access to pseudorapidity
   * @return Constituent pseudorapidity
   */
  virtual double Eta() const = 0;

  /**
   * @brief Access to azimuthal angle
   * @return Constituent azimuthal angle
   */
  virtual double Phi() const = 0;

  /**
   * @brief Checks whether the constituent is from an embedded event
   * @return True if the constituent is from an embedded event
   */
  Bool_t IsFromEmbeddedEvent() const { return fIsFromEmbeddedEvent; }

  /**
   * @brief Get the index of the constituent in the global index map
   * @return Index of the constituent in the global index map
   */
  ULong_t GetGlobalIndex() const  { return fGlobalIndex; }

  /**
   * @brief Specify whether constituent is from embedded event
   * @param[in] isEmbedded Flag of embedded status (true - constituent is from embedded event)
   */
  void SetIsFromEmbeddedEvent(Bool_t isEmbedded) { fIsFromEmbeddedEvent = isEmbedded; }

  /**
   * @brief Set the index in the globl
   * @param[in] index Index of the constituent in the global index map
   */
  void SetGlobalIndex(ULong_t index) { fGlobalIndex = index; }

protected:
  Bool_t 			fIsFromEmbeddedEvent;	///< Flag whether constituent is from embedded event
  ULong_t		  fGlobalIndex;			    ///< Index of the constituent in the global index map

  /// \cond CLASSIMP
  ClassDef(AliEmcalJetConstituent, 0);
  /// \endcond
};

}

}

#endif /* ALIEMCALJETCONSTITUENT_H */
