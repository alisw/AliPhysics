#ifndef ALIITSRAWSTREAMSDDBEAMTESTNOV04_H
#define ALIITSRAWSTREAMSDDBEAMTESTNOV04_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: */

///////////////////////////////////////////////////////////////////////////////
///
/// This class provides access to ITS SDD digits in test beam raw data.
///
///////////////////////////////////////////////////////////////////////////////

#include "AliITSRawStreamSDDBeamTest.h"

class AliRawReader;


class AliITSRawStreamSDDBeamTestNov04: public AliITSRawStreamSDDBeamTest {
  public :
    AliITSRawStreamSDDBeamTestNov04(AliRawReader* rawReader);
    virtual ~AliITSRawStreamSDDBeamTestNov04() {};

    virtual Bool_t   Next();
    virtual Int_t    ReadJitter();
  private :


    ClassDef(AliITSRawStreamSDDBeamTestNov04, 1) // class for reading ITS SDD raw digits
};

#endif
