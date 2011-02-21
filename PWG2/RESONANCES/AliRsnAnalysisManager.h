//
// Class AliRsnAnalysisManager
//
// This is the uppermost level in the analysis task and is buil as a separate object
// since up to this level, the execution of the analysis is not directly related to
// the real analysis task structure.
//
// This object collects all computations whic must be done and processes all of them
// by means of a unique single- or double-loop on all tracks in each event or pair
// of events (in case of mixing), by passing all candidate pairs to all computation
// objects, where each of the latters will process all pairs of tracks which satisfy
// their requirements.
//
// authors : M. Vala       [martin.vala@cern.ch]
//         : A. Pulvirenti [alberto.pulvirenti@ct.infn.it]
//

#ifndef ALIRSNANALYSISMANAGER_H
#define ALIRSNANALYSISMANAGER_H

#include <TObjArray.h>

#include "AliRsnCutSet.h"

class AliRsnPair;
class AliRsnMonitor;

class AliRsnAnalysisManager : public TNamed {
public:

   AliRsnAnalysisManager(const char*name = "RSN");
   AliRsnAnalysisManager(const AliRsnAnalysisManager& copy);
   AliRsnAnalysisManager& operator=(const AliRsnAnalysisManager& copy);
   virtual ~AliRsnAnalysisManager() { }

   virtual void   Add(AliRsnPair *pair);
   virtual void   Add(AliRsnMonitor *monitor);
   virtual void   PrintArray() const;
   virtual void   Print(Option_t *option = "") const;

   void           InitAllPairs(TList*list);
   void           ProcessAllPairs();
   void           ProcessAllPairsMC();
   AliRsnCutSet*  GetGlobalTrackCuts() {return &fGlobalTrackCuts;}

private:

   Bool_t        fAddUsageHist;     // flag to switch on the production of usage histograms
   TList        *fList;             // container for output histograms (external object)
   TObjArray     fPairs;            // collection of pair objects for the different outputs
   TObjArray     fMonitors;         // collection of monitor objects for the different outputs
   AliRsnCutSet  fGlobalTrackCuts;  // a set of cuts which are applied to all tracks for all analysis

   ClassDef(AliRsnAnalysisManager, 1)
};

#endif
