#ifndef ALIITSRAWSTREAMSSDV1_H
#define ALIITSRAWSTREAMSSDV1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
///
/// This class provides access to ITS SSD digits in raw data.
///
///////////////////////////////////////////////////////////////////////////////

#include "AliITSRawStreamSSD.h"

class AliRawReader;


class AliITSRawStreamSSDv1: public AliITSRawStreamSSD {
  public :
    AliITSRawStreamSSDv1(AliRawReader* rawReader);
    virtual ~AliITSRawStreamSSDv1() {};

    virtual Bool_t   Next();

    Int_t            GetADC() const {return fADC;}
    Int_t            GetADModule() const {return fADModule;}

  protected :

  Int_t  fADModule;               // FEROM module
  Int_t  fADC;                    // ADC within the ADModule
         

    ClassDef(AliITSRawStreamSSDv1, 0) // class for reading beam test ITS digits
};

#endif

