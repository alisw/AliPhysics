#ifndef ALIRAWREADERCHAIN_H
#define ALIRAWREADERCHAIN_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////////////////////////////////
///
/// This is a class for reading raw data from a root chain.
///
///////////////////////////////////////////////////////////////////////////////

#include "AliRawReaderRoot.h"

class TChain;
class TFileCollection;

class AliRawReaderChain: public AliRawReaderRoot {
  public :
    AliRawReaderChain();
    AliRawReaderChain(const char* listFileName);
    AliRawReaderChain(TFileCollection *collection);
    AliRawReaderChain(const AliRawReaderChain& rawReader);
    AliRawReaderChain& operator = (const AliRawReaderChain& rawReader);
    virtual ~AliRawReaderChain();

    virtual Bool_t   NextEvent();
    virtual Bool_t   RewindEvents();

  protected :
    TChain*          fChain;        // root chain with raw events

    ClassDef(AliRawReaderChain, 0) // class for reading raw digits from a root file
};

#endif
