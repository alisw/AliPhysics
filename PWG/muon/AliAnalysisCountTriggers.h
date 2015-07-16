#ifndef ALIANALYSISCOUNTRIGGERS_H
#define ALIANALYSISCOUNTRIGGERS_H

/**
 *
 * \class AliAnalysisCountTriggers
 *
 * \brief Simple task to count trigger classes
 *
 *  \author: L. Aphecetche (Subatech)
 */

#ifndef ALIANALYSISTASKSE_H
#  include "AliAnalysisTaskSE.h"
#endif

class AliCounterCollection;
class TH1;

class AliAnalysisCountTriggers : public AliAnalysisTaskSE
{
public:
  AliAnalysisCountTriggers();
  virtual ~AliAnalysisCountTriggers();

  virtual void Terminate(Option_t *);
  
  void UserCreateOutputObjects();

  virtual void UserExec(Option_t* opt);
  
private:
  
  AliAnalysisCountTriggers(const AliAnalysisCountTriggers&); // not implemented (on purpose)
  AliAnalysisCountTriggers& operator=(const AliAnalysisCountTriggers&); // not implemented (on purpose)

  AliVEvent* Event() const;

private:
  
  AliCounterCollection* fEventCounters; //! event counters
  TH1* fHTriggerMask; //! trigger masks
  
  ClassDef(AliAnalysisCountTriggers,1) // a class to count trigger classes
};

#endif

