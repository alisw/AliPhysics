#ifndef ALIANALYSISTASKTRACKDCA_H
#define ALIANALYSISTASKTRACKDCA_H
/* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliAnalysisTaskSE.h"

class TArray;
class THistManager;
class AliAnalysisUtils;
class AliESDtrackCuts;

namespace PWGJE {
  
namespace EMCALJetTasks {

/**
 * \class AliAnalysisTaskTrackDCA
 * \brief Simple analysis task monitoring the track DCA distribution
 *
 * This is a simple analysis task monitoring the track DCA distribution to the primary vertex
 */
class AliAnalysisTaskTrackDCA: public AliAnalysisTaskSE {
public:
  AliAnalysisTaskTrackDCA();
  AliAnalysisTaskTrackDCA(const char *name);
  virtual ~AliAnalysisTaskTrackDCA();

  void UserCreateOutputObjects();
  void UserExec(Option_t  *opt);
  /**
   * Terminate function - does nothing
   * \param opt Option
   */
  void Terminate(Option_t *opt) {}

  /**
   * Set the track selection. The user has to define them in the Add macro.
   * \param trackcuts The track selection cuts to be applied
   */
  void SetTrackCuts(AliESDtrackCuts *trackcuts) { fStandardCuts = trackcuts; }

protected:
  void CreatePtBinning(TArrayD& binning) const;
  void CreateLinearBinning(TArrayD& binning, int nbins, double min, double max) const;

private:
  AliAnalysisTaskTrackDCA(const AliAnalysisTaskTrackDCA &);
  AliAnalysisTaskTrackDCA &operator=(const AliAnalysisTaskTrackDCA &);

  AliAnalysisUtils                    *fAnalysisUtils;            ///< Analysis utils (helper)
  AliESDtrackCuts                     *fStandardCuts;             ///< Standard track cuts
  THistManager                        *fHistos;                   ///< Histogram container

  ClassDef(AliAnalysisTaskTrackDCA, 1)
};

} /* namespace EMCALTriggerJets */

} /* namespace PWGJE */

#endif /* ALIANALYSISTASKTRACKDCA_H_ */
