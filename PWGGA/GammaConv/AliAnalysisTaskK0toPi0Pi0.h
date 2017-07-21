/**************************************************************************************
 * Copyright (C) 2017, Copyright Holders of the ALICE Collaboration                   *
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
#ifndef ALIANALYSISTASKK0TOPI0PI0_H
#define ALIANALYSISTASKK0TOPI0PI0_H

#include <vector>
#include <TString.h>

#include "AliAODConversionPhoton.h"
#include "AliAnalysisTaskSE.h"

class THistManager;
class TList;
class AliAODConversionPhoton;
class AliAODConversionMother;
class AliCaloPhotonCuts;
class AliConvEventCuts;
class AliConversionMesonCuts;
class AliConversionPhotonCuts;
class AliClusterContainer;
class AliV0ReaderV1;
class AliGammaConversionAODBGHandler; 

class AliAnalysisTaskK0toPi0Pi0 : public AliAnalysisTaskSE {
public:
  enum PhotonType_t {
    kPCMPhoton = 0,
    kEMCALPhoton = 1,
    kUndefined = -1
  };

  enum MesonType_t {
    kPi0,
    kK0
  };

  AliAnalysisTaskK0toPi0Pi0();
  AliAnalysisTaskK0toPi0Pi0(const char *name);
  virtual ~AliAnalysisTaskK0toPi0Pi0();

  AliClusterContainer *AddClusterContainer(const char *name);
  void SetNameV0Reader(const char *name) { fV0ReaderName = name; }
  void SetEventCuts(AliConvEventCuts *cuts) { fEventCuts = cuts; }
  void SetConversionPhotonCuts(AliConversionPhotonCuts *cuts) { fConvPhotonCuts = cuts; }
  void SetCaloPhotonCuts(AliCaloPhotonCuts *cuts) { fCaloPhotonCuts = cuts; }
  void SetPi0CutsConvConv(AliConversionMesonCuts *cuts) { fPi0CutsConvConv = cuts; }
  void SetPi0CutsCaloCalo(AliConversionMesonCuts *cuts) { fPi0CutsCaloCalo = cuts; }
  void SetPi0CutsConvCalo(AliConversionMesonCuts *cuts) { fPi0CutsConvCalo = cuts; }
  void SetK0Cuts(AliConversionMesonCuts *cuts){ fK0Cuts = cuts;}

protected:
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *);
  virtual bool UserNotify() { fNewFile = kTRUE; return kTRUE; } // File changed logics in FileChanged function, called when the first event from the file is fully constructed
  virtual void ExecOnce();
  virtual void RunChanged()   {}
  virtual void FileChanged()  {}

  // Selectors
  std::vector<AliAODConversionPhoton> MakeCaloPhotonCandidates(const AliClusterContainer &inputcont, AliCaloPhotonCuts &cuts);
  std::vector<AliAODConversionPhoton> MakeConversionPhotonCandidates(const AliV0ReaderV1 &reader, AliConversionPhotonCuts &cuts);
  std::vector<AliAODConversionMother> SelectMeson(std::vector<AliAODConversionMother> &candidates, AliConversionMesonCuts &cuts, MesonType_t meson, const char *reccase);
  std::vector<AliAODConversionMother> MakePi0Candidates(const std::vector<AliAODConversionPhoton> *primaryLeg, const std::vector<AliAODConversionPhoton> *secondaryLeg, AliConversionMesonCuts &cuts);
  std::vector<AliAODConversionMother> MakeK0ShortCandidates(const std::vector<AliAODConversionMother> *primaryLeg, const std::vector<AliAODConversionMother> *secondaryLeg, AliConversionMesonCuts &cuts);

  void MakePhotonQACalo(const std::vector<AliAODConversionPhoton> &photons, AliConvEventCuts &cuts);
  void MakePhotonQAConv(const std::vector<AliAODConversionPhoton> &photons, AliConvEventCuts &cuts);
  void MakePi0QA(const std::vector<AliAODConversionMother> &pi0s, const char *reccase,TString selectionStatus);
  void MakeK0ShortQA(const std::vector<AliAODConversionMother> &k0s, const char *reccase,TString selectionStatus);

private:
  Bool_t                                       fLocalInitialized;     ///< Check whether the task was initialized (triggers ExecOnce)
  Int_t                                        fCurrentRun;           ///< Current run number (triggers RunChanged)
  Bool_t                                       fNewFile;              ///< New file loaded (triggers fileChanged)
  AliV0ReaderV1                               *fV0Reader;             //!<! V0 reader
  TString                                      fV0ReaderName;         ///< Name of the V0 reader
  AliClusterContainer                         *fClusterContainer;     ///< Cluster container
  Bool_t                                       fIsMC;                 ///< Switch whether we run over data or MC
  Double_t                                     fWeightJetJetMC;       ///< Weight of the jet-jet event
  Double_t                                     fEventPlaneAngle;      ///< Event Plane Angle

  AliConvEventCuts                            *fEventCuts;            ///< Event cuts
  AliConversionPhotonCuts                     *fConvPhotonCuts;       ///< Cuts on conversion photons
  AliCaloPhotonCuts                           *fCaloPhotonCuts;       ///< Calo photon cuts
  AliConversionMesonCuts                      *fPi0CutsConvConv;      ///< Cuts on the pi0 for the conv conv case
  AliConversionMesonCuts                      *fPi0CutsCaloCalo;      ///< Cuts on the pi0 for the calo calo case
  AliConversionMesonCuts                      *fPi0CutsConvCalo;      ///< Cuts on the pi0 for the conv calo case
  AliConversionMesonCuts                      *fK0Cuts;               ///< Cuts on the K0

  AliGammaConversionAODBGHandler              *fBGHandler;            //!<!   Background Handler
  THistManager                                *fHistos;               ///< Container for Histograms
  TList                                       *fOutput;               ///< Global output container


  AliAnalysisTaskK0toPi0Pi0(const AliAnalysisTaskK0toPi0Pi0 &);
  AliAnalysisTaskK0toPi0Pi0 &operator=(const AliAnalysisTaskK0toPi0Pi0 &);

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskK0toPi0Pi0, 1);
  /// \endcond
};

#endif /* ALIANALYSISTASKK0TOPI0PI0_H */
