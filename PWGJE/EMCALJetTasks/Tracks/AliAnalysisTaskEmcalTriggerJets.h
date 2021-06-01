#ifndef ALIANALYSISTASKEMCALTRIGGERJETS_H
#define ALIANALYSISTASKEMCALTRIGGERJETS_H
// Copyright (C) 2017, Copyright Holders of the ALICE Collaboration
// All rights reserved.

#include <AliAnalysisTaskEmcalJet.h>

class AliJetContainer;
class AliPIDResponse;
class THistManager;
class TString;

namespace PWGJE {
  
namespace EMCALJetTasks {

class AliAnalysisTaskEmcalTriggerJets: public AliAnalysisTaskEmcalJet {
public:
  AliAnalysisTaskEmcalTriggerJets();
  AliAnalysisTaskEmcalTriggerJets(const char *name);
  virtual ~AliAnalysisTaskEmcalTriggerJets();

  static AliAnalysisTaskEmcalTriggerJets *AddTaskEmcalTriggerJets(const char *name);

protected:

  virtual void UserCreateOutputObjects();
  virtual void UserExecOnce();
  virtual bool Run();

  void FillJetPIDPlots(const AliEmcalJet *jet, double radius, const char *trigger, const char *detector);
  void FillJetPIDPlotsLeading(const AliVTrack *leading, double ptjet, double radius, const char *trigger, const char *detector);

private:
  AliAnalysisTaskEmcalTriggerJets(const AliAnalysisTaskEmcalTriggerJets &);
  AliAnalysisTaskEmcalTriggerJets &operator=(const AliAnalysisTaskEmcalTriggerJets &);

  AliPIDResponse        *fPIDResponse;        //!<! PID Response handler
  THistManager          *fHistos;             //!<! Histogram handler

  ClassDef(AliAnalysisTaskEmcalTriggerJets, 1)
};

} /* namespace EMCALJetTasks */

} /* namespace PWGJE */

#endif /* ALIANALYSISTASKEMCALTRIGGERJETS_H */
