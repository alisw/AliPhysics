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

#ifndef ALIRSNANALYSISMANAGER_H
#define ALIRSNANALYSISMANAGER_H

#include <TObjArray.h>

#include "AliRsnCutSet.h"

class AliRsnEvent;
class AliRsnPair;

class AliRsnAnalysisManager : public TNamed
{
  public:

    AliRsnAnalysisManager(const char*name = "RSN");
    AliRsnAnalysisManager(const AliRsnAnalysisManager& copy);
    AliRsnAnalysisManager& operator=(const AliRsnAnalysisManager& copy);
    virtual ~AliRsnAnalysisManager() { }

    virtual void   Add(AliRsnPair *pair);
    virtual void   PrintArray() const;
    virtual void   Print(Option_t *option = "") const;

    void           InitAllPairs(TList*list);
    void           ProcessAllPairs();
    void           ProcessAllPairsMC();
    AliRsnCutSet*  GetGlobalTrackCuts() {return &fGlobalTrackCuts;}

  private:
  
    TList        *fList;             // container for output histograms (external object)
    TObjArray     fPairs;            // collection of pair objects for the different outputs
    AliRsnCutSet  fGlobalTrackCuts;  // a set of cuts which are applied to all tracks for all analysis

    ClassDef(AliRsnAnalysisManager, 1)
};

#endif
