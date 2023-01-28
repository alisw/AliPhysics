/**************************************************************************************
 * Copyright (C) 2022, Copyright Holders of the ALICE Collaboration                   *
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
#ifndef AliEmcalFilterAODMCHeaders_H
#define AliEmcalFilterAODMCHeaders_H

class TClonesArray;
class TString;
class AliVEvent;
class AliMCEvent;
class AliNamedArrayI;
class AliAODMCParticle;
class AliAODMCHeader;

#include "AliAnalysisTaskSE.h"

/**
 * @class AliEmcalFilterAODMCHeaders
 * @brief Class to select particles in MC events.
 * @author Salvatore Aiola, Yale Univeristy
 * @ingroup EMCALFWTASKS
 * @since Aug 5, 2012
 */
class AliEmcalFilterAODMCHeaders : public AliAnalysisTaskSE {
 public:

  /**
   * @brief Dummy constructor
   */
  AliEmcalFilterAODMCHeaders();

  /**
   * @brief Main constructor
   *
   * @param name Name of the task
   */
  AliEmcalFilterAODMCHeaders(const char *name);

  /**
   * @brief Destructor
   */
  virtual ~AliEmcalFilterAODMCHeaders();

  /**
   * @brief
   * @param
   */
  Bool_t IsParticleFromBGEvent(Int_t index);

  /**
   * @brief Set the name of the particle output container
   *
   * This container is attached to the input event with the corresponding name.
   * This name has to be used in the user tasks to connect the MC particle container
   * to the particles selected by this instance of the task.
   *
   * @param name Name of the output container attached to the input event
   */
  void SetParticlesOutName(const char *namePart)            { fParticlesOutName = namePart ; }

  /**
   * @brief Set the name of the track output container
   *
   * This container is attached to the input event with the corresponding name.
   * This name has to be used in the user tasks to connect the track container
   * to the tracks selected by this instance of the task.
   *
   * @param name Name of the output container attached to the input event
   */
  void SetTracksOutName(const char *nameTrack)            { fTracksOutName = nameTrack ; }

  /**
   * @brief Set the name of the cluster output container
   *
   * This container is attached to the input event with the corresponding name.
   * This name has to be used in the user tasks to connect the cluster container
   * to the clusters selected by this instance of the task.
   *
   * @param name Name of the output container attached to the input event
   */
  void SetClustersOutName(const char *nameClus)            { fClustersOutName = nameClus ; }

  /**
   * @brief Set the debug level
   * @param fDebugLevel -- 0 for no output, 3 for max output
   */
  void SetDebugLevel(Int_t debug)                       { fDebugLevel = debug      ; }

  /**
   * @brief Create new AliEmcalFilterAODMCHeaders task and add it to the analysis manager
   *
   * @param outname name of the output contaienr
   * @param nk Reject neutrons and K0long
   * @param ch Select only charged particles
   * @param etamax  Max eta acceptance
   * @param physPrim Require physical primary particles
   * @return AliEmcalFilterAODMCHeaders*
   */
  static AliEmcalFilterAODMCHeaders* AddTaskFilterAODMCHeaders(TString nParticlesOut = "MCParticlesNotRejected", TString nTracksOut = "MCTracksNotRejected", TString nClustersOut = "MCClustersNotRejected", Int_t debug = 0);

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
   * to the output container. Set AliEmcalFilterAODMCHeaders::AccpetParticle for the
   * definition of selected particles.
   *
   * @param option Not used
   */
  void UserExec(Option_t *option);

  /**
   * @brief Convert standard MC AOD particles in a new array, and filter if requested (for AOD analysis).
   * @param partIn Input particle container
   * @param partOut Output particle container with selected particles
   * @param partMap Index map between particles in input and output container
   */
  void ProcessParticles();

  /**
   * @brief
   * @param
   */
  void ProcessClusters();

  /**
   * @brief
   * @param
   */
  void ProcessTracks();

  /**
   * @brief
   * @param
   */
  void LinkMothers();

  /**
   * @brief
   * @param
   */
  TString GetNameGenerator(Int_t label, AliAODMCHeader* header);

  TString                   fParticlesOutName;     ///< name of output particle array
  TString                   fParticlesMapName;     //!<! name of the particle map
  TString                   fTracksOutName;        ///< name of output track array
  TString                   fTracksInName;         ///< name of input track array
  TString                   fClustersOutName;      ///< name of output cluster array
  TString                   fClustersInName;       ///< name of input cluster array
  Bool_t                    fInit;                 //!<! true = task initialized
  TClonesArray             *fParticlesIn;          //!<! particle array in (AOD)
  TClonesArray             *fParticlesOut;         //!<! particle array out
  AliNamedArrayI           *fParticlesMap;         //!<! particle index/label
  TClonesArray             *fTracksIn;             //!<! track array in
  TClonesArray             *fTracksOut;            //!<! track array out
  TClonesArray             *fClustersIn;           //!<! cluster array in
  TClonesArray             *fClustersOut;          //!<! cluster array out
  AliVEvent                *fEvent;                //!<! event
  AliMCEvent               *fMC;                   //!<! MC event (ESD)
  Bool_t                    fDisabled;             //!<! Disable task if a problem occurs at initialization
  Int_t                     fDebugLevel;           ///< debug level for interactive debugging
  TClonesArray             *fAODMCTrackArray;      ///< pointer to track array
  AliAODMCHeader           *fMCHeader;             //!

 private:
  AliEmcalFilterAODMCHeaders(const AliEmcalFilterAODMCHeaders&);            // not implemented
  AliEmcalFilterAODMCHeaders &operator=(const AliEmcalFilterAODMCHeaders&); // not implemented

  ClassDef(AliEmcalFilterAODMCHeaders, 1);
};
#endif
