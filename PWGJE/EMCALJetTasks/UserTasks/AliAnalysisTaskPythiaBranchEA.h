#ifndef AliAnalysisTaskPythiaBranchEA_H
#define AliAnalysisTaskPythiaBranchEA_H

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
 * \file AliAnalysisTaskPythiaBranchEA.h
 * \brief This task uses a thermal model of background to generate a collection of
 * truth-level background particles, and attaches it to the event.
 *
 * \author James Mulligan <james.mulligan@yale.edu>, Yale University
 * \date Aug 15, 2018
 */

class AliVEvent;

#include <TRandom3.h>
#include <TString.h>

#include <AliAnalysisTaskSE.h>
#include <THistManager.h>

namespace PWGJE {

namespace EMCALJetTasks {

class AliAnalysisTaskPythiaBranchEA : public AliAnalysisTaskSE {
 public:

  AliAnalysisTaskPythiaBranchEA()                                          ;
  AliAnalysisTaskPythiaBranchEA(const char *name)                          ;
  virtual ~AliAnalysisTaskPythiaBranchEA()                                 ;
  
  static AliAnalysisTaskPythiaBranchEA* AddTaskPythiaBranchEA(const char *outputName = "pyparticles",
                                                                                        const char *pyfilepath ="",
                                                                                        const char *pyfilemask ="",
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
  
  // Thermal model setters
  void SetOutputCollectionName(std::string name)            { fOutputCollectionName = name; }
  
  
  void SetMinEta(Int_t d)                                   { fMinEta = d; }
  void SetMaxEta(Int_t d)                                   { fMaxEta = d; }

  void SetPythiaFilePath(const char* pyfile)                { fPyFilePath = pyfile;}
  void SetPythiaFileMask(const char* pyfilemask)            { fPyFileMask = pyfilemask;}
  void SetNFiles(Int_t n)                                   { fNfiles = n;}


 protected:
  // Event parameters
  bool                        fEventInitialized;                ///< If the event is initialized properly
  AliVEvent*                  fEvent;                           //!<! Pointer to the current event
  TRandom3                    fRandom;                          //!<! Random number generator 
  /////////////////////////////////////////////////////////////////////////////////
  // Thermal model parameters
  
  // General
  std::string                 fOutputCollectionName;            ///< Name of TClonesArray output the thermal particles to the event
  TClonesArray*               fThermalParticlesArray;           //!<! Thermal particle collection
  
  
  Double_t                    fMinEta;                          ///< Min eta for thermal particles
  Double_t                    fMaxEta;                          ///< Max eta for thermal particles

  TString                     fPyFilePath;                      ///< Path to PYTHIA text file
  TString                     fPyFileMask;                      ///< Mask of PYTHIA text file
  TString                     fPyFile;                          ///< Mask of PYTHIA text file
  Int_t                       fNumber;                          //!<! Random number
  Int_t                       fNfiles;                          ///< Number of files   
  /////////////////////////////////////////////////////////////////////////////////
  
  // Histogram output
  TList *                     fOutput;                          //!<! Output for histograms
  THistManager                fHistManager;                     ///< Histogram manager
  FILE*                       fInput;                           //!<! PYTHIA txt file     

 private:
  AliAnalysisTaskPythiaBranchEA(const AliAnalysisTaskPythiaBranchEA&)           ; // not implemented
  AliAnalysisTaskPythiaBranchEA &operator=(const AliAnalysisTaskPythiaBranchEA&); // not implemented

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskPythiaBranchEA, 2);
  /// \endcond
};
}
}
#endif
