#ifndef ALIANALYSISMUMUGLOBAL_H
#define ALIANALYSISMUMUGLOBAL_H

/**
 * \class AliAnalysisMuMuGlobal
 * \brief Basic histogramming of global event properties (vertex, pile-up, background)
 * \author L. Aphecetche (Subatech)
 */

#include "AliAnalysisMuMuBase.h"

class AliMergeableCollectionProxy;

class AliAnalysisMuMuGlobal : public AliAnalysisMuMuBase
{
public:
  AliAnalysisMuMuGlobal();
  virtual ~AliAnalysisMuMuGlobal() {}
  
  void FillHistosForEvent(const char* eventSelection, const char* triggerClassName,
                          const char* centrality);

  void FillHistosForMCEvent(const char* eventSelection, const char* triggerClassName,
                            const char* centrality);

  virtual void DefineHistogramCollection(const char* eventSelection, const char* triggerClassName,
                                         const char* centrality, Bool_t mix =kFALSE);

  Bool_t SelectAnyTriggerClass(const TString& firedTriggerClasses, TString& acceptedTriggerClasses) const;
  

private:
  
  void FillHistosForEvent(AliMergeableCollectionProxy& p);
  void FillHistosForMCEvent(AliMergeableCollectionProxy& p);
  
  ClassDef(AliAnalysisMuMuGlobal,1) // implementation of AliAnalysisMuMuBase for global event properties
};

#endif
