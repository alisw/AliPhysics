#ifndef ALIRSNANALYSISMANAGER_H
#define ALIRSNANALYSISMANAGER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
////////////////////////////////////////////////////////////////////////////////
//
//  Manager for resonance analysis.
//
////////////////////////////////////////////////////////////////////////////////

#include <TObjArray.h>

#include "AliRsnCutSet.h"

class TList;
class AliRsnPair;
class AliRsnMonitor;

class AliRsnAnalysisManager : public TNamed {
public:

   AliRsnAnalysisManager(const char *name = "RSN");
   AliRsnAnalysisManager(const AliRsnAnalysisManager& copy);
   AliRsnAnalysisManager& operator=(const AliRsnAnalysisManager& copy);
   virtual ~AliRsnAnalysisManager() { }

   virtual void   Add(AliRsnPair *pair);
   virtual void   Add(AliRsnMonitor *monitor);
   virtual void   PrintArray() const;
   virtual void   Print(Option_t *option = "") const;

   void           InitAllPairs(TList *list);
   void           ProcessAll(Bool_t pureMC = kFALSE);
   AliRsnCutSet*  GetGlobalTrackCuts() {return &fGlobalTrackCuts;}

private:

   Bool_t        fAddUsageHist;     //  flag to switch on the production of usage histograms
   TList        *fList;             //! container for output histograms (external object)
   TObjArray     fPairs;            //  collection of pair objects for the different outputs
   TObjArray     fMonitors;         //  collection of monitor objects for the different outputs
   AliRsnCutSet  fGlobalTrackCuts;  //  a set of cuts which are applied to all tracks for all analysis

   ClassDef(AliRsnAnalysisManager, 1)
};

#endif
