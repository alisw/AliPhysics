#ifndef ALIITSRAWSTREAMSDDV2_H
#define ALIITSRAWSTREAMSDDV2_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliITSRawStream.h"

class AliRawReader;


class AliITSRawStreamSDDv2: public AliITSRawStream {
  public :
    AliITSRawStreamSDDv2(AliRawReader* rawReader);
    virtual ~AliITSRawStreamSDDv2() {};

    virtual Bool_t   Next();

    Int_t            GetAnode() const {return fCoord1;};
    Int_t            GetTimeBin() const {return fCoord2;};

  private :
    UInt_t           ReadBits();
    Int_t            DecompAmbra(Int_t value) const;

    UInt_t           fData;         // data read for file

    static const UInt_t fgkCodeLength[8];  // length of coded data word
    Int_t            fSkip;            // number of skipped words
    Int_t            fEventId;         // event ID from the header
    Int_t            fCarlosId;        // carlos ID from the header
    Int_t            fChannel;         // current channel
    ULong64_t        fChannelData[2];  // packed data for the 2 channels
    UInt_t           fLastBit[2];      // last filled bit in fChannelData
    UInt_t           fChannelCode[2];  // current channel code
    Bool_t           fReadCode[2];     // next bits are code or data
    UInt_t           fReadBits[2];     // number of bits to read
    Int_t            fTimeBin[2];      // current time bin
    Int_t            fAnode[2];        // current anode number
    Int_t            fLowThreshold[2]; // low Carlos threshold

    ClassDef(AliITSRawStreamSDDv2, 0) // class for reading ITS SDD raw digits
};

#endif
