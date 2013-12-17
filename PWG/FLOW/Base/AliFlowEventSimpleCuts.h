/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

// AliFlowEventSimpleCuts:
// An event cut base class
// origin: Mikolaj Krzewicki (mikolaj.krzewicki@cern.ch)

#ifndef ALIFLOWEVENTSIMPLECUTS_H
#define ALIFLOWEVENTSIMPLECUTS_H

#include "TNamed.h"

class AliFlowEventSimpleCuts : public TNamed {

 public:
  AliFlowEventSimpleCuts();
  AliFlowEventSimpleCuts(const char* name, const char* title = "AliFlowEventSimpleCuts");
  AliFlowEventSimpleCuts(const AliFlowEventSimpleCuts& someCuts);
  AliFlowEventSimpleCuts& operator=(const AliFlowEventSimpleCuts& someCuts);
  virtual  ~AliFlowEventSimpleCuts();
  
  virtual Bool_t IsSelected(TObject* obj, TObject* mcobj);
  virtual void SetCentralityPercentileRange(Float_t min, Float_t max){ fCentralityPercentileMin=min;
                                                               fCentralityPercentileMax=max;
                                                               fCutCentralityPercentile=kTRUE; }

 protected:
  Bool_t fCutCentralityPercentile; //cut on centrality perc. from AliESDCentrality
  Float_t fCentralityPercentileMax; // max centr. perc
  Float_t fCentralityPercentileMin; // min centr. perc

  ClassDef(AliFlowEventSimpleCuts,1)
};

#endif


