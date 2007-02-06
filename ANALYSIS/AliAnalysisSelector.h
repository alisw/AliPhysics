#ifndef ALIANALYSISSELECTOR_H
#define ALIANALYSISSELECTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Author: Andrei Gheata, 31/05/2006

//==============================================================================
//   AliAnalysisSelector - Transparent selector class instantiated by an
// analysis manager object.
//==============================================================================

#ifndef ROOT_TSelector
#include "TSelector.h"
#endif

class AliAnalysisManager;

class AliAnalysisSelector : public TSelector {

protected:
   Bool_t              fInitialized; // Flag that initialization was done
   AliAnalysisManager *fAnalysis;    // Analysis manager to be processed
   
private:
   AliAnalysisSelector(const AliAnalysisSelector&);            // not implemented
   AliAnalysisSelector& operator=(const AliAnalysisSelector&); // not implemented
   void                RestoreAnalysisManager();

public:
   AliAnalysisSelector() : TSelector(), fInitialized(kFALSE), fAnalysis(NULL) {}
   AliAnalysisSelector(AliAnalysisManager *mgr);
   virtual ~AliAnalysisSelector();
   
   virtual int         Version() const {return 1;}
   virtual void        Init(TTree *tree);
   virtual void        Begin(TTree *);
   virtual void        SlaveBegin(TTree *tree);
   virtual Bool_t      Notify() {return kTRUE;}   
   virtual Bool_t      Process(Long64_t entry);
   virtual void        SlaveTerminate();
   virtual void        Terminate();

   ClassDef(AliAnalysisSelector,0)  //A class for processing jobs using AliAnalysisManager
};

#endif
