#ifndef ALIITSRAWSTREAMSPD_H
#define ALIITSRAWSTREAMSPD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliITSRawStream.h"


class AliITSRawStreamSPD: public AliITSRawStream {
  public :
    AliITSRawStreamSPD(AliRawReader* rawReader);
    virtual ~AliITSRawStreamSPD() {};

    virtual Bool_t   Next();

    Int_t            GetRow() const {return fCoord1;};
    Int_t            GetColumn() const {return fCoord2;};

  private :
    UShort_t         fData;         // data read for file
    UInt_t           fOffset;       // offset for cell column
    UInt_t           fHitCount;     // counter of hits

    ClassDef(AliITSRawStreamSPD, 0) // class for reading ITS SPD raw digits
};

#endif
