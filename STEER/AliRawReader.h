#ifndef ALIRAWREADER_H
#define ALIRAWREADER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TObject.h>
#ifdef __CINT__
class fstream;
#else
#include <Riostream.h>
#endif

struct AliMiniHeader {
  UInt_t    fSize;
  UChar_t   fDetectorID;
  UChar_t   fMagicWord[3];
  UChar_t   fVersion;
  UChar_t   fCompressionFlag;
  UShort_t  fDDLID;
};

class AliRawReader: public TObject {
  public :
    AliRawReader(const char* fileName, Bool_t addNumber = kTRUE);
    virtual ~AliRawReader();

    inline Int_t     GetDetectorID() const {return fMiniHeader.fDetectorID;};
    inline Int_t     GetDDLID() const {return fMiniHeader.fDDLID;};
    inline Int_t     GetVersion() const {return fMiniHeader.fVersion;};
    inline Bool_t    IsCompressed() const {return fMiniHeader.fCompressionFlag != 0;};

    Bool_t           ReadNextInt(UInt_t& data);
    Bool_t           ReadNextShort(UShort_t& data);
    Bool_t           ReadNextChar(UChar_t& data);

  protected :
    Bool_t           OpenNextFile();

    Bool_t           ReadMiniHeader();

    const char*      fFileName;    // name of input files
    Int_t            fFileNumber;  // number of current input file
    fstream*         fStream;      // stream of raw digits
    AliMiniHeader    fMiniHeader;  // current mini header
    Int_t            fCount;       // counter of bytes to be read for current DDL

    ClassDef(AliRawReader, 0) // base class for reading raw digits
};

#endif
