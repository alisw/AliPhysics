/*
 * AliAnalysisTaskTrackDensityData.h
 *
 *  Created on: Mar 11, 2016
 *      Author: markus
 */

#ifndef ALIANALYSISTASKTRACKDENSITYDATA_H
#define ALIANALYSISTASKTRACKDENSITYDATA_H

#include "AliAnalysisTaskEmcalJet.h"
#include <TString.h>

class AliEmcalJet;
class AliParticleContainer;
class THistManager;

class AliEmcalTrackSelection;

namespace PWGJE {

namespace EMCALJetTasks {

class AliEMCalTriggerBinningComponent;

class AliAnalysisTaskTrackDensityData : public AliAnalysisTaskEmcalJet {
public:
  AliAnalysisTaskTrackDensityData();
  AliAnalysisTaskTrackDensityData(const char *name);
  virtual ~AliAnalysisTaskTrackDensityData();

  void SetEmcalTrackSelection(AliEmcalTrackSelection *sel) { fTrackSelection = sel; }
  void SetNameJetContainer(TString name) { fNameJetContainer = name; }
  void SetNameTrackContainer(TString name) { fNameTrackContainer = name; }

  AliEMCalTriggerBinningComponent *GetBinningHandler() const { return fBinHandler; }

protected:
  virtual void UserCreateOutputObjects();
  virtual bool Run();

  int GetParticleMultiplicity(const AliEmcalJet &jet, const AliParticleContainer &partcont, double ptmin, double ptmax, double rmin, double rmax) const;
  void FindJetPtBin(const AliEmcalJet *const jet, double &ptmin, double &ptmax) const;

private:
  THistManager                          *fHistos;                 //!<!
  AliEmcalTrackSelection    *fTrackSelection;         /// EMCAL track selection
  AliEMCalTriggerBinningComponent       *fBinHandler;             /// Binning handler

  TString                               fNameJetContainer;        /// name of the jet container
  TString                               fNameTrackContainer;      /// name of the track container

  ClassDef(AliAnalysisTaskTrackDensityData, 1);
};

} /* namespace EMCALJetTasks */

} /* namespace PWGJE */

#endif /* ALIANALYSISTASKTRACKDENSITYDATA_H */
