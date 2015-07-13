#ifndef ALIANALYSISTASKCHARGEDPARTICLESREFMC_H
#define ALIANALYSISTASKCHARGEDPARTICLESREFMC_H
/* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliAnalysisTaskSE.h"
#include <TString.h>

class TArrayD;
class TClonesArray;

class AliAnalysisUtil;
class AliAODTrack;
class AliESDtrack;
class AliESDtrackCuts;
class AliGenPythiaEventHeader;
class AliVParticle;
class AliMCEvent;

namespace EMCalTriggerPtAnalysis {

class AliEMCalHistoContainer;
/**
 * \class AliAnalysisTaskChargedParticlesRefMC
 * \brief Unit test class for charged particle distributions (MC case)
 */
class AliAnalysisTaskChargedParticlesRefMC: public AliAnalysisTaskSE {
public:
  AliAnalysisTaskChargedParticlesRefMC();
  AliAnalysisTaskChargedParticlesRefMC(const char *name);
  virtual ~AliAnalysisTaskChargedParticlesRefMC();

  void                        UserCreateOutputObjects();
  void                        UserExec(Option_t *);
  Bool_t                      UserNotify();
  void                        Terminate(Option_t *) {}

protected:
  void                        CreateOldPtBinning(TArrayD &binning) const;
  void                        CreateNewPtBinning(TArrayD &binning) const;

  Bool_t                      TrackSelectionESD(AliESDtrack *track);
  Bool_t                      TrackSelectionAOD(AliAODTrack *track);
  TString                     GetFiredTriggerClasses(const TClonesArray * triggerpatches);
  Bool_t                      PythiaInfoFromFile(const char* currFile, Float_t &fXsec, Float_t &fTrials, Int_t &pthard) const;
  AliGenPythiaEventHeader    *GetPythiaHeader() const;
  Bool_t                      IsPhysicalPrimary(const AliVParticle *const part, AliMCEvent *const mcevent);


  AliESDtrackCuts                 *fTrackCuts;                ///< Standard track selection
  AliAnalysisUtils                *fAnalysisUtil;             ///< Event selection
  AliEMCalHistoContainer          *fHistos;                   ///< Histogram manager

  // Monte-Carlo specific information
  Double_t                        fPtHard;                    ///< event pt hard
  Int_t                           fPtHardBin;                 ///< event pt hard bin
  Int_t                           fNTrials;                   ///< event trials
  Float_t                         fXsection;                  ///< x-section from pythia header

private:
  AliAnalysisTaskChargedParticlesRefMC(const AliAnalysisTaskChargedParticlesRefMC &);
  AliAnalysisTaskChargedParticlesRefMC &operator=(const AliAnalysisTaskChargedParticlesRefMC &);

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskChargedParticlesRefMC, 1);
  /// \endcond
};

} /* namespace EMCalTriggerPtAnalysis */

#endif /* ALIANALYSISTASKCHARGEDPARTICLESREFMC_H */
