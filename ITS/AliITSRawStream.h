#ifndef ALIITSRAWSTREAM_H
#define ALIITSRAWSTREAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TObject.h>
#include "AliRawReader.h"


class AliITSRawStream: public TObject {
  public :
    AliITSRawStream(AliRawReader* rawReader);

    virtual Bool_t   Next() = 0;

    inline Int_t     GetModuleID() const {return fModuleID;};
    inline Int_t     GetPrevModuleID() const {return fPrevModuleID;};
    inline Bool_t    IsNewModule() const {return fModuleID != fPrevModuleID;};
    inline Int_t     GetCoord1() const {return fCoord1;};
    inline Int_t     GetCoord2() const {return fCoord2;};
    inline Int_t     GetSignal() const {return fSignal;};

  protected :
    AliRawReader*    fRawReader;    // object for reading the raw data

    Int_t            fModuleID;     // index of current module
    Int_t            fPrevModuleID; // index of previous module
    Int_t            fCoord1;       // current 1st coordinate
                                    //  SPD: column cell number (z)
                                    //  SDD: anode cell number (z)
                                    //  SSD: N/P, flag for side
    Int_t            fCoord2;       // current 2nd coordinate
                                    //  SPD: row cell number (y)
                                    //  SDD: time bin number (y)
                                    //  SSD: strip number
    Int_t            fSignal;       // signal in ADC counts

    ClassDef(AliITSRawStream, 0) // base class for reading ITS raw digits
};

#endif
