#ifndef ALIITSRAWSTREAMSDDV3_H
#define ALIITSRAWSTREAMSDDV3_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
///
/// This class provides access to ITS SDD digits in test beam raw data.
///
///////////////////////////////////////////////////////////////////////////////

#include "AliITSRawStreamSDD.h"

class AliRawReader;


class AliITSRawStreamSDDv3: public AliITSRawStreamSDD {
  public :
    AliITSRawStreamSDDv3(AliRawReader* rawReader);
    virtual ~AliITSRawStreamSDDv3() {};

    virtual Bool_t   Next();
    virtual Int_t    GetJitter();
  private :


    ClassDef(AliITSRawStreamSDDv3, 1) // class for reading ITS SDD raw digits
};

#endif
