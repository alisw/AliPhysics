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
#ifndef ALIPARTICLECONTAINER_H
#define ALIPARTICLECONTAINER_H

#include <TArrayI.h>
#include <AliVParticle.h>

class AliVEvent;
class AliTLorentzVector;

#include "AliEmcalContainer.h"
#if !(defined(__CINT__) || defined(__MAKECINT__))
#include "AliEmcalContainerIndexMap.h"
#endif

#if !(defined(__CINT__) || defined(__MAKECINT__))
typedef EMCALIterableContainer::AliEmcalIterableContainerT<AliVParticle, EMCALIterableContainer::operator_star_object<AliVParticle> > AliParticleIterableContainer;
typedef EMCALIterableContainer::AliEmcalIterableContainerT<AliVParticle, EMCALIterableContainer::operator_star_pair<AliVParticle> > AliParticleIterableMomentumContainer;
#endif

/**
 * @class AliParticleContainer
 * @brief Container for particles within the EMCAL framework
 * @ingroup EMCALCOREFW
 * @author Marta Verweij
 * @author Salvatore Aiola <salvatore.aiola@cern.ch>, Yale University
 *
 * Container with name, TClonesArray and cuts for particles
 */
class AliParticleContainer : public AliEmcalContainer {
 public:
  enum EChargeCut_t { kNoChargeCut, kCharged, kNeutral, kPositiveCharge, kNegativeCharge };

  /**
   * @brief Default constructor.
   */
  AliParticleContainer();

  /**
   * @brief Standard constructor.
   * @param name Name of the particle branch (TClonesArray)
   */
  AliParticleContainer(const char *name);

  /**
   * @brief Destructor
   */
  virtual ~AliParticleContainer(){;}

  /**
   * @brief Index operator 
   * @param[in] index Index of the particle in the array
   * @return Particle at the given index
   * 
   * Providing access to track in the container with the given index.
   */
  virtual TObject *operator[](int index) const { return GetParticle(index); }

  /**
   * @brief Select particle based on particle level cuts
   * @param[in] vp Particle for which the selection is performed
   * @param[out] rejectionReason Bitmap with the reason why the particle was accepted.
   * Note: The value is not set to 0 in the function in order to combine the information
   * with other selection steps.
   * @return True if the particle is accepted, False otherwise
   * 
   * Apply cuts on the particle properties. Implemented are
   * - Monte-Carlo label
   * - Charge
   * - Generator index
   */
  virtual Bool_t              ApplyParticleCuts(const AliVParticle* vp, UInt_t &rejectionReason) const;

  /**
   * @brief Select particle based on kinematic cuts
   * @param[in] mom Momentum vector for which
   * @param[out] rejectionReason  Bitmap with the reason why the particle was accepted.
   * Note: The value is not set to 0 in the function in order to combine the information
   * with other selection steps.
   * @return True if the momentum vector was selected under the given cuts, false otherwise
   * 
   * Apply kinematical cuts to the momentum vector provided. In addition
   * to the standard kinematical cuts in \f$ p_{t} \f$, \f$ \eta \f$ and
   * \f$ \phi \f$, implemented in the AliEmcalContainer::ApplyKinematicCuts,
   * also the distance to the TPC sector boundary is checked.
   */
  virtual Bool_t              ApplyKinematicCuts(const AliTLorentzVector& mom, UInt_t &rejectionReason) const;

  /**
   * @brief Check whether particle at a given index is accepted
   * @param i Index of the particle in the container
   * @param[out] rejectionReason Bitmap containing selections criteria which were not passed
   * @return True if the particle is accepted, false otherwise 
   * 
   * The function is used by AliEmcalContainer in order to select the particle at a given position 
   */
  virtual Bool_t              AcceptObject(Int_t i, UInt_t &rejectionReason) const              { return AcceptParticle(i, rejectionReason);}

  /**
   * @brief Check whether a given particle is accepted
   * @param obj Particle to be checked
   * @param[out] rejectionReason Bitmap containing selections criteria which were not passed
   * @return True if the particle is accepted, false otherwise 
   * 
   * The function is used by AliEmcalContainer in order to select a given particle
   */
  virtual Bool_t              AcceptObject(const TObject* obj, UInt_t &rejectionReason) const   { return AcceptParticle(dynamic_cast<const AliVParticle*>(obj), rejectionReason);}

  /**
   * @brief Check whether a particle is selected among the selection cut defined for this container
   * @param vp Particle for which to perform the selection
   * @param rejectionReason Bitmap with the reason why the particle was accepted.
   * Note: The value is not set to 0 in the function in order to combine the information
   * with other selection steps.
   * @return True if the particle is accepted, false otherwise.
   * 
   * Perform full particle selection consisting of kinematical and particle property
   * selection on the particle provided. In case the particle is rejected, the
   * reason for the rejection is encoded in the bitmap rejectionReason.
   */
  virtual Bool_t              AcceptParticle(const AliVParticle* vp, UInt_t &rejectionReason) const        ;

  /**
   * @brief Check whether the particle at a given index in the container is selected among the selection cut defined for this container 
   * @param[in] i Index of the particle to select.
   * @param[out] rejectionReason Bitmap with the reason why the particle was accepted.
   * Note: The value is not set to 0 in the function in order to combine the information
   * with other selection steps.
   * @return True if the particle was accepted, false otherwise
   * 
   * Perform full particle selection consisting of kinematical and particle property
   * selection on the \f$ i^{th} \f$ particle in the array. In case the particle is
   * rejected, the reason for the rejection is encoded in the bitmap rejectionReason.
   */
  virtual Bool_t              AcceptParticle(Int_t i, UInt_t &rejectionReason) const                       ;

  Double_t                    GetParticlePtCut()                        const   { return GetMinPt()     ; }
  Double_t                    GetParticleEtaMin()                       const   { return GetMinEta()    ; }
  Double_t                    GetParticleEtaMax()                       const   { return GetMaxEta()    ; }
  Double_t                    GetParticlePhiMin()                       const   { return GetMinPhi()    ; }
  Double_t                    GetParticlePhiMax()                       const   { return GetMaxPhi()    ; }
  void                        SetParticlePtCut(Double_t cut)                    { SetMinPt(cut)         ; }
  void                        SetParticleEtaLimits(Double_t min, Double_t max)  { SetEtaLimits(min, max); }
  void                        SetParticlePhiLimits(Double_t min, Double_t max)  { SetPhiLimits(min, max); }

  /**
   * @brief Get the leading particle in the container. 
   * @param[in] opt Options for the selection of the leading particle
   * @return Leading particle in the container
   * 
   * If "p" is contained in the parameter opt, 
   * then the absolute momentum is use instead of the transverse momentum.
   */
  virtual AliVParticle       *GetLeadingParticle(const char* opt="")         ;

  /**
   * @brief Get \f$ i^{th} \f$ particle in the container.
   * @param[in] i Index of the particle to access
   * @return Parrticle at the given index (NULL if the index is out of range)
   */
  virtual AliVParticle       *GetParticle(Int_t i=-1)                   const;

  /**
   * @brief Get \f$ i^{th} \f$ particle in the container if it is accepted.
   * @param[in] i Index of the particle
   * @return Particle at the index if it is accepted, NULL otherwise
   */
  virtual AliVParticle       *GetAcceptParticle(Int_t i=-1)             const;

  /**
   * @brief Get the next accepted particle in the array. 
   * @return Next accepted particle in the array (NULL if the end is reached)
   * @deprecated Only for backward compatibility - use accept_iterator instead
   */
  virtual AliVParticle       *GetNextAcceptParticle()                        ;

  /**
   * @brief Get the next particle in the array.
   * @return Next particle in the array (NULL if the end is reached)
   * @deprecated Only for backward compatibility - use all_iterator instead.
   */
  virtual AliVParticle       *GetNextParticle()                              ;

  /**
   * @brief Get particle momentum for a given particle under a user-defined mass hypothesis
   * @param[out] mom Momentum vector to be filled
   * @param[in] part Particle from which the momentum information is obtained.
   * @param[in] mass (Optional) Mass hypothesis
   * @return True if the particle is valid, false if it is invalid
   * 
   * Retrieve momentum information of a particle (part) and fill a TLorentzVector
   * with it. In case the optional parameter mass is provided, it is used as mass
   * hypothesis, otherwise the mass hypothesis from the particle itself is used.
   */
  virtual Bool_t              GetMomentumFromParticle(TLorentzVector &mom, const AliVParticle* part, Double_t mass) const;

  /**
   * @brief Get the momentum vector for a given particle using the default mass hypothesis
   * @param[out] mom Momentum vector of the particle provided
   * @param[in] part Particle from which to obtain the momentum information
   * @return Always true
   * 
   * Fills a TLorentzVector with the momentum information of the particle provided
   * under a global mass hypothesis.
   */
  virtual Bool_t              GetMomentumFromParticle(TLorentzVector &mom, const AliVParticle* part) const;

  /**
   * Get the momentum for the particle at a given index in the container
   * @param[out] mom Momentum vector of the \f$ i^{th} \f$ particle in the array
   * @param[in] i Index of th particle to check
   * @return True if the request was successfull, false otherwise
   * 
   * Fills a TLorentzVector with the monentum infomation of the
   * \f$ i^{th} \f$ particle in the container, using a global
   * mass hypothesis. In case the provided index is out of
   * range, false is returned as return value.
   */
  virtual Bool_t              GetMomentum(TLorentzVector &mom, Int_t i) const;

  /**
   * @brief Get the momentum vector for the particle at a given position in the container if it is accepted
   * @param[out] mom Momentum vector of the accepted particle
   * @param[in] i Index to check
   * @return True if the request was successfull, false otherwise
   * 
   * Fills a TLorentzVector with the monentum infomation of the
   * \f$ i^{th} \f$ accepted particle in the container, using a
   * global mass hypothesis. In case the provided index is out of
   * range, or the particle under the index is not accepted, false
   * is returned as return value.
   */
  virtual Bool_t              GetAcceptMomentum(TLorentzVector &mom, Int_t i) const;

  /**
   * @brief Get momentum vector for the next particle in the container
   * @param[out] mom Momentum vector of the next particle
   * @return True if the request was successfull, false otherwise
   * @deprecated Old style iterator - use all_iterator instead
   * 
   * Fills a TLorentzVector with the monentum infomation of the
   * next particle in the container, using a global mass hypothesis.
   * In case the iterator reached the end of the array, false
   * is returned as return value.
   */
  virtual Bool_t              GetNextMomentum(TLorentzVector &mom);

  /**
   * @brief Get the momentum vector of the next accepted particle in the container
   * @param[out] mom Momentum vector of the next particle in the array
   * @return True if the request was successfull, false (no more entries) otherwise
   * @deprecated Old style iterator - use accept_iterator instead
   * 
   * Fills a TLorentzVector with the monentum infomation of the
   * next accepted particle in the container, using a global
   * mass hypothesis. In case the iteration reached the end of
   * the array, false is returned as return value.
   */
  virtual Bool_t              GetNextAcceptMomentum(TLorentzVector &mom);

  /**
   * @brief Get numebr of particles in this container
   * @return Number of particles
   */
  Int_t                       GetNParticles()                           const   {return GetNEntries();}

  /**
   * @brief Get number of accepted particles handled by this container 
   * @return Number of selected particles under the given particle selection
   * 
   * In order to get this number, the selection has to be applied to each 
   * particle within this function.
   */
  Int_t                       GetNAcceptedParticles()                   const;

  void                        SetMinDistanceTPCSectorEdge(Double_t min)         { fMinDistanceTPCSectorEdge = min; }
  void                        SetCharge(EChargeCut_t c)                         { fChargeCut = c       ; }
  void                        SelectHIJING(Bool_t s)                            { if (s) fGeneratorIndex = 0; else fGeneratorIndex = -1; }
  void                        SetGeneratorIndex(Short_t i)                      { fGeneratorIndex = i  ; }

  /**
   * @brief Connect container to input event
   * @param event Input event containing the array with content.
   * 
   * Connect the container to the array with content stored inside the virtual event.
   * The object name in the event must match the name given in the constructor.
   *
   * Additionally register the array into the index map.
   */
  void                        SetArray(const AliVEvent * event);

  /**
   * @brief Get title of the container
   * @return Title of the container
   * 
   * Make a title of the container name based on the min \f$ p_{t} \f$ used
   * in the particle selection process.
   */
  const char*                 GetTitle() const;

#if !(defined(__CINT__) || defined(__MAKECINT__))
  /// Get the EMCal container utils associated with particle containers
  static const AliEmcalContainerIndexMap <TClonesArray, AliVParticle>& GetEmcalContainerIndexMap() { return fgEmcalContainerIndexMap; }

  /**
   * @brief Create an iterable container interface over all objects in the
   * EMCAL container.
   * @return iterable container over all objects in the EMCAL container
   */
  const AliParticleIterableContainer      all() const;

  /**
   * @brief Create an iterable container interface over accepted objects in the
   * EMCAL container.
   * @return iterable container over accepted objects in the EMCAL container
   */
  const AliParticleIterableContainer      accepted() const;

  /**
   * @brief Create an iterable container interface over all objects in the
   * EMCAL container.
   * @return iterable container over all objects in the EMCAL container
   */
  const AliParticleIterableMomentumContainer      all_momentum() const;

  /**
   * @brief Create an iterable container interface over accepted objects in the
   * EMCAL container.
   * @return iterable container over accepted objects in the EMCAL container
   */
  const AliParticleIterableMomentumContainer      accepted_momentum() const;
#endif

 protected:

#if !(defined(__CINT__) || defined(__MAKECINT__))
  static AliEmcalContainerIndexMap <TClonesArray, AliVParticle> fgEmcalContainerIndexMap; //!<! Mapping from containers to indices
#endif

  Double_t                    fMinDistanceTPCSectorEdge;      ///< require minimum distance to edge of TPC sector edge
  EChargeCut_t                fChargeCut;                     ///< select particles according to their charge
  Short_t                     fGeneratorIndex;                ///< select MC particles with generator index (default = -1 = switch off selection)

 private:
  AliParticleContainer(const AliParticleContainer& obj); // copy constructor
  AliParticleContainer& operator=(const AliParticleContainer& other); // assignment

  ClassDef(AliParticleContainer,11);
};

/**
 * Unit test for the iterators. Comparing iterators against for-loop of particles.
 * All particles selected in the for-loop must be found in order to pass the test.
 * @ingroup EMCALCOREFW
 * @param cont Particle container used for the test.
 * @param iteratorType type of the iterator (0 = accept_iterator, 1 = all_iterator)
 * @param verbose Switch on verbosity in case of true
 * @return Result of the unit test (0 - passed, 1 - particles missing, 2 - excess particles)
 */
int TestParticleContainerIterator(const AliParticleContainer *const cont, int iteratorType = 0, bool verbose = false);

#endif

