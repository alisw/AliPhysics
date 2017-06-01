#ifndef ALIANALYSISTASKEMCALTRIGGERJETS_H
#define ALIANALYSISTASKEMCALTRIGGERJETS_H
// Copyright (C) 2017, Copyright Holders of the ALICE Collaboration
// All rights reserved.

#include <AliAnalysisTaskEmcalJet.h>

class AliJetContainer;
class THistManager;
class TString;

namespace EmcalTriggerJets {

class AliAnalysisTaskEmcalTriggerJets: public AliAnalysisTaskEmcalJet {
public:
  AliAnalysisTaskEmcalTriggerJets();
  AliAnalysisTaskEmcalTriggerJets(const char *name);
  virtual ~AliAnalysisTaskEmcalTriggerJets();

  static AliAnalysisTaskEmcalTriggerJets *AddTaskEmcalTriggerJets(const char *name);

protected:

  virtual void UserCreateOutputObjects();
  virtual bool Run();

private:
  AliAnalysisTaskEmcalTriggerJets(const AliAnalysisTaskEmcalTriggerJets &);
  AliAnalysisTaskEmcalTriggerJets &operator=(const AliAnalysisTaskEmcalTriggerJets &);

  THistManager          *fHistos;             //!<! Histogram handler

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskEmcalTriggerJets, 1)
  /// \endcond
};

} /* namespace EmcalTriggerJets */

#endif /* ALIANALYSISTASKEMCALTRIGGERJETS_H */
