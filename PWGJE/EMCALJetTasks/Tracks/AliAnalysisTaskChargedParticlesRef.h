#ifndef ALIANALYSISTASKCHARGEDPARTICLESREF_H
#define ALIANALYSISTASKCHARGEDPARTICLESREF_H
/* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliAnalysisTaskSE.h"

class TArrayD;
class AliAnalysisUtils;
class AliAODTrack;
class AliESDtrackCuts;
class AliESDtrack;

namespace EMCalTriggerPtAnalysis {

class AliEMCalHistoContainer;

/**
 * \class AliAnalysisTaskChargedParticlesRef
 * \brief Unit test class for charged particle distributions
 */
class AliAnalysisTaskChargedParticlesRef : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskChargedParticlesRef();
  AliAnalysisTaskChargedParticlesRef(const char *name);
  virtual ~AliAnalysisTaskChargedParticlesRef();

  void UserCreateOutputObjects();
  void UserExec(Option_t *);
  void Terminate(Option_t *) {}

protected:
  void CreateOldPtBinning(TArrayD &binning) const;
  void CreateNewPtBinning(TArrayD &binning) const;

  Bool_t TrackSelectionESD(AliESDtrack *track);
  Bool_t TrackSelectionAOD(AliAODTrack *track);

  AliESDtrackCuts                 *fTrackCuts;                ///< Standard track selection
  AliAnalysisUtils                *fAnalysisUtil;             ///< Event selection
  AliEMCalHistoContainer          *fHistos;                   ///< Histogram manager

private:
  AliAnalysisTaskChargedParticlesRef(const AliAnalysisTaskChargedParticlesRef &);
  AliAnalysisTaskChargedParticlesRef &operator=(const AliAnalysisTaskChargedParticlesRef &);

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskChargedParticlesRef, 1);
  /// \endcond
};

} /* namespace EMCalTriggerPtAnalysis */

#endif /* ALIANALYSISTASKCHARGEDPARTICLESREF_H */
