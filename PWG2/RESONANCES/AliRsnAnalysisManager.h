//
// Class AliRsnAnalysisManager
//
// This is the uppermost level of analysis objects collection.
// It contains a list of pair managers, which all will process
// a pool of events passed to this object, and fill their histograms.
//
// The utility of this object is to define a unique implementation
// of the whole processing, which can then be included in the different
// designs of AnalysisTask provided for SE and ME analysis.
//
// The base architecture is still AliRsnVManager, but in this case
// all the objects in the list will be AliRsnPairManager's.
//
// author     : M. Vala       [martin.vala@cern.ch]
// revised by : A. Pulvirenti [alberto.pulvirenti@ct.infn.it]
//

#ifndef AliRsnAnalysisManager_H
#define AliRsnAnalysisManager_H

#include <TROOT.h>

#include "AliRsnVManager.h"
#include "AliRsnPairManager.h"

class AliRsnAnalysisManager : public AliRsnVManager
{
  public:

    AliRsnAnalysisManager(const char*name = "defaultAnalysisMgr");

    virtual void   Add(AliRsnPairManager *pair);
    virtual void   AddConfig(TString config, TString prefix, TString functionName = "");
    virtual void   PrintArray() const;
    virtual void   Print(Option_t *option = "") const;

    TList*         InitAllPairMgrs();
    void           ProcessAllPairMgrs(AliRsnPIDIndex *pidIndexes1, AliRsnEvent *ev1, AliRsnPIDIndex *pidIndexes2 = 0, AliRsnEvent *ev2 = 0);

  private:

    ClassDef(AliRsnAnalysisManager, 1)
};

#endif
