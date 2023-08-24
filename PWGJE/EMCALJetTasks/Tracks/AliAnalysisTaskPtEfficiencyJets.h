#ifndef ALIANALYSISTASKPTEFFICIENCYJETS_H
#define ALIANALYSISTASKPTEFFICIENCYJETS_H
/* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <AliAnalysisTaskEmcalJet.h>
#include <TString.h>

class AliAnalysisUtils;
class AliEmcalJet;
class AliJetContainer;
class AliVParticle;
class AliVTrack;
class TNtuple;

class AliEmcalTrackSelection;

namespace PWGJE {
  
namespace EMCALJetTasks {

class AliAnalysisTaskPtEfficiencyJets: public AliAnalysisTaskEmcalJet {
public:
  AliAnalysisTaskPtEfficiencyJets();
  AliAnalysisTaskPtEfficiencyJets(const char *name);
  virtual ~AliAnalysisTaskPtEfficiencyJets();

  virtual void UserCreateOutputObjects();
  virtual Bool_t Run();

  void SetMCJetContainer(const char *name) { fMCJetContainer = name; }
  void SetTrackCuts(AliEmcalTrackSelection *cuts) { fTrackCuts = cuts; }

protected:
  AliVTrack *FindAssociatedTrack(AliVParticle *trueParticle);
  AliEmcalJet *FindAssociatedJet(AliVParticle *trueParticle, AliJetContainer *jets);
  bool SelectTrueParticle(AliVParticle *part);

private:
  AliAnalysisTaskPtEfficiencyJets(const AliAnalysisTaskPtEfficiencyJets &);
  AliAnalysisTaskPtEfficiencyJets &operator=(const AliAnalysisTaskPtEfficiencyJets &);

  AliAnalysisUtils                    *fAnalysisUtils;
  TString                             fMCJetContainer;
  AliEmcalTrackSelection  *fTrackCuts;
  TNtuple                             *fTrackNtuple;

  ClassDef(AliAnalysisTaskPtEfficiencyJets, 1);
};

} /* namespace EMCALTriggerJets */

} /* namespace PWGJE */

#endif /* ALIANALYSISTASKPTEFFICIENCYJETS_H */
