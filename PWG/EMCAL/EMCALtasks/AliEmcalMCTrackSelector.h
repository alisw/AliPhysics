/**************************************************************************************
 * Copyright (C) 2012, Copyright Holders of the ALICE Collaboration                   *
 * All rights reserved.                                                               *
 *                                                                                    *
 * Redistribution and use in source and binary forms, with or without                 *
 * modification, are permitted provided that the following conditions are met:        *
 *     * Redistributions of source code must retain the above copyright               *
 *       notice, this list of conditions and the following disclaimer.                *
 *     * Redistributions in binary form must reproduce the above copyright            *
 *       notice, this list of conditions and the following disclaimer in the          *
 *       documentation and/or other materials provided with the distribution.         *
 *     * Neither the name of the <organization> nor the                               *
 *       names of its contributors may be used to endorse or promote products         *
 *       derived from this software without specific prior written permission.        *
 *                                                                                    *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND    *
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED      *
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE             *
 * DISCLAIMED. IN NO EVENT SHALL ALICE COLLABORATION BE LIABLE FOR ANY                *
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES         *
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;       *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND        *
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT         *
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS      *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                       *
 **************************************************************************************/
#ifndef ALIEMCALMCTRAKCSELECTOR_H
#define ALIEMCALMCTRAKCSELECTOR_H

class TClonesArray;
class TString;
class AliVEvent;
class AliMCEvent;
class AliNamedArrayI;
class AliAODMCParticle;

#include "AliAnalysisTaskSE.h"

/**
 * @class AliEmcalMCTrackSelector
 * @brief Class to select particles in MC events.
 * @author Salvatore Aiola, Yale Univeristy
 * @ingroup EMCALFWTASKS
 * @since Aug 5, 2012
 */
class AliEmcalMCTrackSelector : public AliAnalysisTaskSE {
 public:

  /**
   * @brief Dummy constructor
   */
  AliEmcalMCTrackSelector();

  /**
   * @brief Main constructor
   *
   * @param name Name of the task
   */
  AliEmcalMCTrackSelector(const char *name);

  /**
   * @brief Destructor
   */
  virtual ~AliEmcalMCTrackSelector() {}

  /**
   * @brief Select only physical primary particles
   * @param s If true only physical primary particles are used
   */
  void SetOnlyPhysPrim(Bool_t s)                        { fOnlyPhysPrim     = s    ; }

  /**
   * @brief Select only charged particles
   * @param c If true only charged particles are selected
   */
  void SetChargedMC(Bool_t c = kTRUE)                   { fChargedMC        = c    ; }

  /**
   * @brief Set the eta acceptance
   * @param e Maximum eta acceptance
   */
  void SetEtaMax(Double_t e)                            { fEtaMax           = e    ; }

  /**
   * @brief Reject neutrons and K0long particles
   * @param r If true neutrons and K0long particles are rejected
   */
  void SetRejectNK(Bool_t r = kTRUE)                    { fRejectNK         = r    ; }

  /**
   * @brief Reject photon in case it is the mother of another photon
   *
   * In order to mimic processes PYTHIA8 puts mothers and daugthers of
   * the process on the stack, which can lead in case of photons to
   * double counting. Photon mothers that are duaghters of photon mothers
   * need to be rejected.
   *
   * @param doReject If true photons are rejected if they are mothers of other photons
   */
  void SetRejectPhotonMother(bool doReject)             { fRejectPhotonMothers = doReject; }

  /**
   * @brief Make Pi0 stable
   *
   * Select Pi0 as stable (physical primary) particle, and reject particles where
   * an ancestor in the decay chain is a pi0. Ancestor is usually the mother, however
   * in order to drop conversions of photons from pi0 to electrons any particle in the
   * decay chain is checked.
   *
   * @param doSelect If true Pi0s are selected as primary particles and their daughters are rejected
   */
  void SelectStablePi0(bool doSelect)                   { fSelectStablePi0 = doSelect; }

  void SetOnlyHIJING(Bool_t s)                          { fOnlyHIJING       = s    ; }

  /**
   * @brief Set the name of the output container
   *
   * This container is attached to the input event with the corresponding name.
   * This name has to be used in the user tasks to connect the MC particle container
   * to the particles selected by this instance of the task.
   *
   * @param name Name of the output container attached to the input event
   */
  void SetParticlesOutName(const char *name)            { fParticlesOutName = name ; }

  /**
   * @brief Set the name of the input track container
   *
   * If there is a track container attached to the event other than the default,
   * this is where to specify the name. For example, if AliEmcalRejectMCBackground
   * is run first, it will attach a new track container to the event. Otherwise,
   * leave this as default.
   *
   * @param nTrackCont Name of the track container attached to the input event
   */
  void SetTrackContainerName(const char *nTrackCont)    {fTrackContainerName = nTrackCont ; }

  void SetIsContainerSpecified(bool isContSpecified)    {fIsContainerSpecified = isContSpecified ; }


  /**
   * @brief Create new AliEmcalMCTrackSelector task and add it to the analysis manager
   *
   * @param outname name of the output contaienr
   * @param nk Reject neutrons and K0long
   * @param ch Select only charged particles
   * @param etamax  Max eta acceptance
   * @param physPrim Require physical primary particles
   * @return AliEmcalMCTrackSelector*
   */
  static AliEmcalMCTrackSelector* AddTaskMCTrackSelector(TString outname = "mcparticles", TString nTrackCont = "usedefault", Bool_t nk = kFALSE, Bool_t ch = kFALSE, Double_t etamax = 1, Bool_t physPrim = kTRUE);

 protected:

  /**
   * @brief Creating user output
   *
   * Not used in this task
   */
  void UserCreateOutputObjects() {}

  /**
   * @brief Main event loop
   *
   * Run selection of particles and convert them to AliAODMCParticles and copy them
   * to the output container. Set AliEmcalMCTrackSelector::AccpetParticle for the
   * definition of selected particles.
   *
   * @param option Not used
   */
  void UserExec(Option_t *option);

  /**
   * @brief Check whether paricle is selected
   *
   * Acceptance criteria:
   * - Physical primary
   * - Charged / neutral
   * - Is neutron or K0long
   * - Eta range
   * - Generator index (for HIJING prodctions)
   *
   * @param part Particle to be checked
   * @return True if the particle is accepted, false otherwise
   */
  virtual Bool_t            AcceptParticle(AliAODMCParticle* part) const;

  /**
   * @brief Convert MC particles in MC AOD articles (for ESD analysis).
   * @param mcEvent Input event
   * @param partOut Output particle container with selected particles
   * @param partMap
   */
  void                      ConvertMCParticles(AliMCEvent* mcEvent, TClonesArray* partOut, AliNamedArrayI* partMap=0);

  /**
   * @brief Convert standard MC AOD particles in a new array, and filter if requested (for AOD analysis).
   * @param partIn Input particle container
   * @param partOut Output particle container with selected particles
   * @param partMap Index map between particles in input and output container
   */
  void                      CopyMCParticles(TClonesArray* partIn, TClonesArray* partOut, AliNamedArrayI* partMap=0);

  /**
   * @brief Recursive check if the mother or another ancestor is a Pi0
   *
   * @param part Particle to be checked
   * @return True if any ancestor is a pi0, false otherwise
   */
  bool                      IsFromPi0Mother(const AliVParticle &part) const;


  TString                   fParticlesOutName;     ///< name of output particle array
  TString                   fTrackContainerName;   ///< name of input track container
  bool                      fIsContainerSpecified; ///< true = track container other than default
  Bool_t                    fOnlyPhysPrim;         ///< true = only physical primary particles
  Bool_t                    fRejectNK;             ///< true = reject K_0^L and neutrons
  Bool_t                    fChargedMC;            ///< true = only charged particles
  Bool_t                    fOnlyHIJING;           ///< true = only HIJING particles
  Bool_t                    fRejectPhotonMothers;  ///< Reject photons that are mothers of other photons
  Bool_t                    fSelectStablePi0;      ///< Select Pi0 as stable particle, and reject daughters of Pi0 in case they are physical primary
  Double_t                  fEtaMax;               ///< maximum eta to accept particles
  TString                   fParticlesMapName;     //!<! name of the particle map
  Bool_t                    fInit;                 //!<! true = task initialized
  TClonesArray             *fParticlesIn;          //!<! particle array in (AOD)
  TClonesArray             *fParticlesOut;         //!<! particle array out
  AliNamedArrayI           *fParticlesMap;         //!<! particle index/label
  AliVEvent                *fEvent;                //!<! event
  AliMCEvent               *fMC;                   //!<! MC event (ESD)
  Bool_t                    fIsESD;                //!<! ESD or AOD analysis
  Bool_t                    fDisabled;             //!<! Disable task if a problem occurs at initialization

 private:
  AliEmcalMCTrackSelector(const AliEmcalMCTrackSelector&);            // not implemented
  AliEmcalMCTrackSelector &operator=(const AliEmcalMCTrackSelector&); // not implemented

  ClassDef(AliEmcalMCTrackSelector, 6);
};
#endif
