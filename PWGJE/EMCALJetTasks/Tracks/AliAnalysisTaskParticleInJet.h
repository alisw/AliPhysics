/*
 * AliAnalysisTaskParticleInJet.h
 *
 *  Created on: Feb 17, 2016
 *      Author: markus
 */

#ifndef ALIANALYSISTASKPARTICLEINJET_H
#define ALIANALYSISTASKPARTICLEINJET_H

#include "AliAnalysisTaskEmcalJet.h"

#include <vector>

class TArrayD;
class THistManager;

class AliVParticle;
class AliVTrack;
class AliEmcalTrackSelection;

class AliAnalysisTaskParticleInJet: public AliAnalysisTaskEmcalJet {
public:
  AliAnalysisTaskParticleInJet();
  AliAnalysisTaskParticleInJet(const char *name);
  virtual ~AliAnalysisTaskParticleInJet();

  void SetTrackSelection(AliEmcalTrackSelection *sel)     { fTrackSelection = sel; }
  void SetParticleContainerNameRec(TString name)          { fParticleContainerNameRec = name; }
  void SetParticleContainerNameMC(TString name)           { fParticleContainerNameMC = name;}
  void SetJetContainerNameRec(TString name)               { fJetContainerNameRec = name; }
  void SetJetContainerNameMC(TString name)                { fJetContainerNameMC = name; }

protected:
  void UserCreateOutputObjects();
  Bool_t Run();

private:
  std::vector<const AliVParticle *> GetSelectedParticles(AliParticleContainer *const cont) const;
  Bool_t AcceptParticle(const AliVParticle * const part) const;
  Bool_t AcceptTrack(AliVTrack * const track) const;
  Bool_t IsPhysicalPrimary(const AliVParticle * const part) const;

  void CreatePtBinning(TArrayD &binning) const;
  void CreateLinearBinning(TArrayD& binning, int nbins, double min, double max) const;

  THistManager                    *fHistMgr;
  AliEmcalTrackSelection          *fTrackSelection;

  // Container names
  TString                         fParticleContainerNameRec;
  TString                         fParticleContainerNameMC;
  TString                         fJetContainerNameRec;
  TString                         fJetContainerNameMC;

  ClassDef(AliAnalysisTaskParticleInJet, 1);

};

#endif /* ALIANALYSISTASKPARTICLEINJET_H */
