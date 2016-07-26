#ifndef ALIEMCALTRIGGERDECISIONCONTAINER_H
#define ALIEMCALTRIGGERDECISIONCONTAINER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Author: Markus Fasel

#include <TList.h>
#include <TNamed.h>

class AliEmcalTriggerDecision;

class AliEmcalTriggerDecisionContainer: public TNamed {
public:
  AliEmcalTriggerDecisionContainer();
  AliEmcalTriggerDecisionContainer(const char *name);
  virtual ~AliEmcalTriggerDecisionContainer() {}

  void Reset();

  void AddTriggerDecision(AliEmcalTriggerDecision * const decision);
  const AliEmcalTriggerDecision *FindTriggerDecision(const char *name) const;

protected:
  TList     fContainer;         // List of trigger decisions

  ClassDef(AliEmcalTriggerDecisionContainer, 1);    // Container for trigger decisions
};

#endif /* ALIEMCALTRIGGERDECISIONCONTAINER_H */
