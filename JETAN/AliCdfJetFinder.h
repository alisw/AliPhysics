#ifndef ALICDFJETFINDER_H
#define ALICDFJETFINDER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


//---------------------------------------------------------------------
//
//
//---------------------------------------------------------------------

#include "AliJetFinder.h"

class AliCdfJetHeader;

class AliCdfJetFinder : public AliJetFinder
  {
  public:

    AliCdfJetFinder();
    virtual ~AliCdfJetFinder();

    void           CreateOutputObjects(TList *histos);
    void           FindJets();
    virtual void   FinishRun();

  protected:
    AliCdfJetFinder ( const AliCdfJetFinder& jf );
    AliCdfJetFinder& operator = ( const AliCdfJetFinder& jf );

    TList         *fHistos;    // List of histograms
    ClassDef ( AliCdfJetFinder, 1 )
  };//
#endif
