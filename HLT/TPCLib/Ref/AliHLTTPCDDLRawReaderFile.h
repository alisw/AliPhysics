// @(#) $Id$

#ifndef ALIHLTTPCDDLRAWREADERFILE_H
#define ALIHLTTPCDDLRAWREADERFILE_H

#include "AliHLTTPCDDLRawReader.h"

class AliHLTTPCDDLRawReaderFile: public AliHLTTPCDDLRawReader 
{
  public :
    AliHLTTPCDDLRawReaderFile(const Char_t* name, Bool_t addnum = kTRUE);
    virtual ~AliHLTTPCDDLRawReaderFile();

    virtual Bool_t   ReadMiniHeader();
    virtual Bool_t   ReadNextData(UChar_t*& data);

    virtual Bool_t   Reset();

  protected :
    Bool_t           OpenNextFile();

    virtual Bool_t   ReadNext(UChar_t* data, Int_t size);

    Char_t*          fFileName;    //! name of input files
    Int_t            fFileNumber;  //  number of current input file
    fstream*         fStream;      //! stream of raw digits
    UChar_t*         fBuffer;      //!  buffer for payload
    Int_t            fBufferSize;  //  size of fBuffer in bytes

    ClassDef(AliHLTTPCDDLRawReaderFile, 1) //AliHLTTPCDDLRawReaderFile
};

#endif
