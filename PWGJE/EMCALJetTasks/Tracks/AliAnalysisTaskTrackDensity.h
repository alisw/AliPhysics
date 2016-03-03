/*
 * AliAnalysisTaskTrackDensity.h
 *
 *  Created on: Mar 2, 2016
 *      Author: markus
 */

#ifndef ALIANALYSISTASKTRACKDENSITY_H
#define ALIANALYSISTASKTRACKDENSITY_H

#include "AliAnalysisTaskEmcalJet.h"

#include <TArrayD.h>
#include <TString.h>

class AliEmcalJet;
class AliParticleContainer;

namespace EMCalTriggerPtAnalysis {

class AliAnalysisTaskTrackDensity : public AliAnalysisTaskEmcalJet {
public:
  AliAnalysisTaskTrackDensity();
  AliAnalysisTaskTrackDensity(const char *name);
  virtual ~AliAnalysisTaskTrackDensity();

  void SetJetRadiusBinning(TArrayD binning) { fJetRadii = binning; }
  void SetJtPtBinning(TArrayD binning) { fJetPtBins = binning; }
  void SetParticlePtSteps(TArrayD binning) { fPtMinSteps = binning; }
  void SetParticlePtBinning(TArrayD binning) { fParticlePtBinning = binning; }

  void SetMCJetContainer(TString contname) { fMCJetContainerName = contname; }
  void SetMCParticleContainer(TString contname) { fMCParticleContainerName = contname; }

protected:

  virtual void UserCreateOutputObjects();
  virtual bool Run();

  int GetParticleMultiplicity(const AliEmcalJet &jet, const AliParticleContainer &partcont, double ptmin, double ptmax, double rmin, double rmax) const;
  void FindJetPtBin(const AliEmcalJet *const jet, double &ptmin, double &ptmax) const;

private:
  THistManager                *fHistos;                     //!<! Histogram manager

  TString                     fMCJetContainerName;          /// Name of the MC jet container
  TString                     fMCParticleContainerName;     /// Name of the MC particle container

  TArrayD                     fJetRadii;                    ///
  TArrayD                     fJetPtBins;                   ///
  TArrayD                     fPtMinSteps;                  ///
  TArrayD                     fParticlePtBinning;           ///

  ClassDef(AliAnalysisTaskTrackDensity, 1);
};

} /* namespace EMCalTriggerPtAnalysis */

#endif /* ALIANALYSISTASKTRACKDENSITY_H */
