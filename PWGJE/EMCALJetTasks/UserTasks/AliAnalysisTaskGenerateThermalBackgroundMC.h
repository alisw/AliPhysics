#ifndef AliAnalysisTaskGenerateThermalBackgroundMC_H
#define AliAnalysisTaskGenerateThermalBackgroundMC_H

/**********************************************************************************
 * Copyright (C) 2016, Copyright Holders of the ALICE Collaboration                *
 * All rights reserved.                                                            *
 *                                                                                 *
 * Redistribution and use in source and binary forms, with or without              *
 * modification, are permitted provided that the following conditions are met:     *
 *   * Redistributions of source code must retain the above copyright              *
 *     notice, this list of conditions and the following disclaimer.               *
 *   * Redistributions in binary form must reproduce the above copyright           *
 *     notice, this list of conditions and the following disclaimer in the         *
 *     documentation and/or other materials provided with the distribution.        *
 *   * Neither the name of the <organization> nor the                              *
 *     names of its contributors may be used to endorse or promote products        *
 *     derived from this software without specific prior written permission.       *
 *                                                                                 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND *
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED   *
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE          *
 * DISCLAIMED. IN NO EVENT SHALL ALICE COLLABORATION BE LIABLE FOR ANY             *
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES      *
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;    *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND     *
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT      *
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS   *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                    *
 * *********************************************************************************/

/**
 * \file AliAnalysisTaskGenerateThermalBackgroundMC.h
 * \brief This task uses a thermal model of background to generate a collection of
 * truth-level background particles, and attaches it to the event.
 *
 * \author James Mulligan <james.mulligan@yale.edu>, Yale University
 * \date Aug 15, 2018
 */

class AliVEvent;

#include <TRandom3.h>

#include <AliAnalysisTaskSE.h>
#include <THistManager.h>

class AliAnalysisTaskGenerateThermalBackgroundMC : public AliAnalysisTaskSE {
 public:

  AliAnalysisTaskGenerateThermalBackgroundMC()                                          ;
  AliAnalysisTaskGenerateThermalBackgroundMC(const char *name)                          ;
  virtual ~AliAnalysisTaskGenerateThermalBackgroundMC()                                 ;
  
  static AliAnalysisTaskGenerateThermalBackgroundMC* AddTaskGenerateThermalBackgroundMC(const char *outputName = "thermalparticles",
                                                                                        const Double_t beta    = 0.3,
                                                                                        const char *suffix     = "");
  // Methods from AliAnalysisTaskSE
  void UserCreateOutputObjects();
  void UserExec(Option_t* option);
  
  // Additional steering functions
  void ExecOnce();
  void Run();
  void FillHistograms();
  
  // Analysis functions
  void                        CreateNewObjectBranch();
  Double_t                    GetDeltaR(Double_t eta, Double_t phi, Double_t etaRef, Double_t phiRef);
  
  // Thermal model setters
  void SetOutputCollectionName(std::string name)            { fOutputCollectionName = name; }
  void SetChargedParticleFraction(Double_t d)               { fChargedParticleFraction = d; }
  
  void SetAlpha(Double_t d)                                 { fAlpha = d; }
  void SetBeta(Double_t d)                                  { fBeta = d; }
  void SetMinChargedPt(Int_t d)                             { fMinChargedPt = d; }
  void SetMinNeutralPt(Int_t d)                             { fMinNeutralPt = d; }
  void SetMaxPt(Int_t d)                                    { fMaxPt = d; }
  
  void SetUseGaussianForN(bool b)                           { fUseGaussianForN = b; }
  void SetNGaussianMean(Double_t d)                         { fNGaussianMean = d; }
  void SetNGaussianSigma(Double_t d)                        { fNGaussianSigma = d; }
  void SetMinN(Int_t d)                                     { fMinN = d; }
  void SetMaxN(Int_t d)                                     { fMaxN = d; }
  
  void SetMinEta(Int_t d)                                   { fMinEta = d; }
  void SetMaxEta(Int_t d)                                   { fMaxEta = d; }

 protected:
  // Event parameters
  bool                        fEventInitialized;                ///< If the event is initialized properly
  AliVEvent*                  fEvent;                           //!<! Pointer to the current event
  
  /////////////////////////////////////////////////////////////////////////////////
  // Thermal model parameters
  
  // General
  std::string                 fOutputCollectionName;            ///< Name of TClonesArray output the thermal particles to the event
  TClonesArray*               fThermalParticlesArray;           //!<! Thermal particle collection
  Double_t                    fChargedParticleFraction;         ///< Fraction of thermal particles set to have nonzero charge
  TRandom3                    fRandom;                          //!<! Random number generator
  
  // pT, based on a gamma distribution. Default choice is a=2, b=<pT>/2. Choosing a=1 gives exponential.
  TF1*                        fPtDistribution;                  //!<! Distribution to sample particle pT from
  Double_t                    fAlpha;                           ///< Value of a in the gamma distribution f_gamma ~ x^(a-1) * e^(-x/b)
  Double_t                    fBeta;                            ///< Value of b in the gamma distribution f_gamma ~ x^(a-1) * e^(-x/b)
  Double_t                    fMinChargedPt;                    ///< Min pT of charged thermal particles
  Double_t                    fMinNeutralPt;                    ///< Min pT of neutral thermal particles
  Double_t                    fMaxPt;                           ///< Max pT of thermal particles
  
  // N particles
  bool                        fUseGaussianForN;                 ///< Use a Gaussian for number of particles per event. Otherwise, use flat distribution.
  Double_t                    fNGaussianMean;                   ///< Mean of Gaussian (if enabled)
  Double_t                    fNGaussianSigma;                  ///< Sigma of Gaussian (if enabled)
  Int_t                       fMinN;                            ///< Min number of particles in thermal model (for flat distribution)
  Int_t                       fMaxN;                            ///< Max number of particles in thermal model (for flat distribution)
  
  // Other
  Double_t                    fMinEta;                          ///< Min eta for thermal particles
  Double_t                    fMaxEta;                          ///< Max eta for thermal particles
  
  /////////////////////////////////////////////////////////////////////////////////
  
  // Histogram output
  TList *                     fOutput;                          //!<! Output for histograms
  THistManager                fHistManager;                     ///< Histogram manager

 private:
  AliAnalysisTaskGenerateThermalBackgroundMC(const AliAnalysisTaskGenerateThermalBackgroundMC&)           ; // not implemented
  AliAnalysisTaskGenerateThermalBackgroundMC &operator=(const AliAnalysisTaskGenerateThermalBackgroundMC&); // not implemented

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskGenerateThermalBackgroundMC, 1);
  /// \endcond
};
#endif
