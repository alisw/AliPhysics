/************************************************************************************
 * Copyright (C) 2013, Copyright Holders of the ALICE Collaboration                 *
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
#ifndef ALIMCPARTICLECONTAINER_H
#define ALIMCPARTICLECONTAINER_H

class AliVEvent;
class AliVParticle;
class AliTLorentzVector;

#include <TArrayC.h>

#include "AliAODMCParticle.h"
#include "AliParticleContainer.h"

#if !(defined(__CINT__) || defined(__MAKECINT__))
typedef EMCALIterableContainer::AliEmcalIterableContainerT<AliAODMCParticle, EMCALIterableContainer::operator_star_object<AliAODMCParticle> > AliMCParticleIterableContainer;
typedef EMCALIterableContainer::AliEmcalIterableContainerT<AliAODMCParticle, EMCALIterableContainer::operator_star_pair<AliAODMCParticle> > AliMCParticleIterableMomentumContainer;
#endif

/**
 * @class AliMCParticleContainer
 * @brief Container for MC-true particles within the EMCAL framework
 * @ingroup EMCALCOREFW
 * @author Salvatore Aiola <salvatore.aiola@cern.ch>, Yale University
 * 
 * Container with name, TClonesArray and cuts for particles
 */
class AliMCParticleContainer : public AliParticleContainer {
 public:

  /**
   * @brief Default constructor.
   */
  AliMCParticleContainer();

  /**
   * @brief Standard constructor.
   * @param[in] name Name of the container (= name of the array operated on)
   */
  AliMCParticleContainer(const char *name);

  /**
   * @brief Destructor
   */
  virtual ~AliMCParticleContainer(){;}

  /**
   * @brief Apply MC particle cuts, e.g. primary particle selection.
   * @param[in] vp Particle to be checked
   * @param[in] rejectionReason Bitmap encoding the reason why the
   * particle was rejected. Note: The variable is not set to NULL
   * inside this function before changing its value.
   * @return True if the particle is accepted, false otherwise
   */
  virtual Bool_t              ApplyMCParticleCuts(const AliAODMCParticle* vp, UInt_t &rejectionReason) const;

  /**
   * @brief Check whether particle at a given index is accepted
   * @param i Index of the particle in the container
   * @param[out] rejectionReason Bitmap containing selections criteria which were not passed
   * @return True if the particle is accepted, false otherwise 
   * 
   * The function is used by AliEmcalContainer in order to select the particle at a given position 
   */
  virtual Bool_t              AcceptObject(Int_t i, UInt_t &rejectionReason) const { return AcceptMCParticle(i, rejectionReason);}

  /**
   * @brief Check whether a given particle is accepted
   * @param obj Particle to be checked
   * @param[out] rejectionReason Bitmap containing selections criteria which were not passed
   * @return True if the particle is accepted, false otherwise 
   * 
   * The function is used by AliEmcalContainer in order to select a given particle
   */
  virtual Bool_t              AcceptObject(const TObject* obj, UInt_t &rejectionReason) const { return AcceptMCParticle(dynamic_cast<const AliAODMCParticle*>(obj), rejectionReason);}

  /**
   * @brief Check whether a particle is selected among the selection cut defined for this container
   * @param vp Particle for which to perform the selection
   * @param rejectionReason Bitmap with the reason why the particle was accepted.
   * Note: The value is not set to 0 in the function in order to combine the information
   * with other selection steps.
   * @return True if the particle is accepted, false otherwise.
   */
  virtual Bool_t              AcceptParticle(Int_t i, UInt_t &rejectionReason) const { return AcceptMCParticle(i, rejectionReason);}

  /**
   * @brief Check whether the particle at a given index in the container is selected among the selection cut defined for this container 
   * @param[in] i Index of the particle to select.
   * @param[out] rejectionReason Bitmap with the reason why the particle was accepted.
   * Note: The value is not set to 0 in the function in order to combine the information
   * with other selection steps.
   * @return True if the particle was accepted, false otherwise
   */
  virtual Bool_t              AcceptParticle(const AliVParticle* vp, UInt_t &rejectionReason) const { return AcceptMCParticle(dynamic_cast<const AliAODMCParticle*>(vp), rejectionReason);}

  /**
   * @brief Check whether the MC particle is accepted
   * @param[in] vp Particle to be checked
   * @param[in] rejectionReason Bitmap encoding the reason why the
   * particle was rejected. Note: The variable is not set to NULL
     * inside this function before changing its value.
   * @return True if the particle is accepted, false otherwise
   * 
   * Perform full MC particle selection for the particle vp, consisting
   * of kinematical particle selection and MC-specific cuts
   */
  virtual Bool_t              AcceptMCParticle(const AliAODMCParticle* vp, UInt_t &rejectionReason) const;

  /**
   * @brief Check whether the MC particle at a given index in the container is accepted
   * @param[in] i Index of the particle to check
   * @param[in] rejectionReason Bitmap encoding the reason why the
   * particle was rejected. Note: The variable is not set to NULL
   * inside this function before changing its value.
   * @return True if the particle is accepted, false otherwise
   * 
   * Perform full MC particle selection for the particle vp, consisting
   * of kinematical particle selection and MC-specific cuts
   */
  virtual Bool_t              AcceptMCParticle(Int_t i, UInt_t &rejectionReason) const;

  /**
   * @brief Get MC particle using the MC label
   * @param[in] lab Label of the particle
   * @return pointer to particle if particle is found, NULL otherwise
   */
  virtual AliAODMCParticle   *GetMCParticleWithLabel(Int_t lab)         const;

  /**
   * @brief Get MC particle using the MC label in case the particle was accepted
   * @param[in] lab Label of the particle
   * @return pointer to particle if particle is found and accepted, NULL otherwise
   */
  virtual AliAODMCParticle   *GetAcceptMCParticleWithLabel(Int_t lab)        ;

  /**
   * @brief Get the highest energetic MC particle 
   * @param opt Option among which to select the leading partice
   * @return Leading particle in the container
   */
  virtual AliAODMCParticle   *GetLeadingMCParticle(const char* opt="")        { return static_cast<AliAODMCParticle*>(GetLeadingParticle(opt)); }

  /**
   * @brief Get particle at index in the container
   * @param[in] i Index of the particle in the container
   * @return pointer to particle if particle is found, NULL otherwise
   */
  virtual AliAODMCParticle   *GetMCParticle(Int_t i=-1)                 const;

  /**
   * @brief Get track at index in the container
   * @param[in] i Index of the particle in the container
   * @return pointer to particle if particle is accepted, NULL otherwise
   */
  virtual AliAODMCParticle   *GetAcceptMCParticle(Int_t i=-1)           const;

  /**
   * @brief Get next accepted particle in the container selected using the track cuts provided.
   * @return Next accepted particle (NULL if the end of the array is reached)
   * @deprecated Old style iterator - for compatibility reasons, use AliParticleContainer::accept_iterator instead
   */
  virtual AliAODMCParticle   *GetNextAcceptMCParticle()                      ;

  /**
   * @brief Get next particle in the container
   * @return Next track in the container (NULL if end of the container is reached)
   * @deprecated Old style iterator - for compatibility reasons, use AliParticleContainer::all_iterator instead
   */
  virtual AliAODMCParticle   *GetNextMCParticle()                            ;

  /**
   * @brief Get the particle at a given index in the container 
   * @param i Index in the container
   * @return Particle at the index (nullptr if the index is larger than then number of particles in the container)
   */
  virtual AliVParticle       *GetParticle(Int_t i=-1)                   const { return GetMCParticle(i)           ; }

  /**
   * @brief Get the particle at a given index in the container in case it was accepted
   * @param i Index in the container
   * @return Particle at the index in case it was accepted (nullptr if the particle was not accepted or he index is larger than then number of particles in the container)
   */
  virtual AliVParticle       *GetAcceptParticle(Int_t i=-1)             const { return GetAcceptMCParticle(i)     ; }

  /**
   * @brief Get the next accepted particle in the container
   * @return Next accepted particle in the container (nullptr if no more acceped particles can be found)
   * @deprecated Old style iterator - for compatibility reasons, use AliParticleContainer::all_iterator instead
   */
  virtual AliVParticle       *GetNextAcceptParticle()                         { return GetNextAcceptMCParticle()  ; }

  /**
   * @brief Get the next particle in the container
   * @return Next particle in the container (nullptr if no more particles can be found)
   * @deprecated Old style iterator - for compatibility reasons, use AliParticleContainer::accepted_iterator instead
   */
  virtual AliVParticle       *GetNextParticle()                               { return GetNextMCParticle()        ; }

  void                        SetMCFlag(UInt_t m)                             { fMCFlag          = m ; }

  /**
   * @brief Set min. pt cut applied on charged particles
   * @param minpt Min. pt (in GeV/c)
   */
  void                        SetMinPtCharged(Double_t minpt)                 { fMinPtCharged = minpt; }

  /**
   * @brief Set min. pt cut applied on neutral particles
   * @param minpt Min. pt (in GeV/c)
   */
  void                        SetMinPtNeutral(Double_t minpt)                 { fMinPtNeutral = minpt; }

  /**
   * @brief Set min. energy cut applied on charged particles
   * @param minE Min. energy (in GeV)
   */
  void                        SetMinECharged(Double_t minE)                   { fMinECharged = minE; }

  /**
   * @brief Set min. energy cut applied on neutral particles
   * @param minE Min. energy (in GeV)
   */
  void                        SetMinENeutral(Double_t minE)                   { fMinENeutral = minE; }

  /**
   * @brief Require particle to be a physical primary particle
   * @param s If true only physical primary particles are selected 
   */
  void                        SelectPhysicalPrimaries(Bool_t s)               { if (s) fMCFlag |=  AliAODMCParticle::kPhysicalPrim ;   }

  /**
   * @brief Get the title of the MC particle container
   * @return Title of the container
   * 
   * Build title of the container consisting of the container name
   * and a string encoding the minimum \f$ p_{t} \f$ cut applied
   * in the kinematic particle selection.
   */
  const char*                 GetTitle() const;

#if !(defined(__CINT__) || defined(__MAKECINT__))

  /**
   * @brief Create an iterable container interface over all objects in the
   * EMCAL container.
   * @return iterable container over all objects in the EMCAL container
   */
  const AliMCParticleIterableContainer      all() const;

  /**
   * @brief Create an iterable container interface over accepted objects in the
   * EMCAL container.
   * @return iterable container over accepted objects in the EMCAL container
   */
  const AliMCParticleIterableContainer      accepted() const;

  /**
   * @brief Create an iterable container interface over all objects in the
   * EMCAL container.
   * @return iterable container over all objects in the EMCAL container
   */
  const AliMCParticleIterableMomentumContainer      all_momentum() const;

  /**
   * @brief Create an iterable container interface over accepted objects in the
   * EMCAL container.
   * @return iterable container over accepted objects in the EMCAL container
   */
  const AliMCParticleIterableMomentumContainer      accepted_momentum() const;
#endif

 protected:
  virtual TString             GetDefaultArrayName(const AliVEvent * const ev) const { return "mcparticles"; }

  UInt_t                      fMCFlag;                        ///< select MC particles with flags
  Double_t                    fMinPtCharged;                  ///< min. charged particle pt
  Double_t                    fMinPtNeutral;                  ///< min. neutral particle pt
  Double_t                    fMinECharged;                   ///< min. charged particle energy
  Double_t                    fMinENeutral;                   ///< min. neutral particle energy

 private:
  AliMCParticleContainer(const AliMCParticleContainer& obj); // copy constructor
  AliMCParticleContainer& operator=(const AliMCParticleContainer& other); // assignment

  /// \cond CLASSIMP
  ClassDef(AliMCParticleContainer,2);
  /// \endcond
};

#endif

