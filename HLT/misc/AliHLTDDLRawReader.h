// @(#) $Id$

#ifndef ALIL3DDLRAWREADER_H
#define ALIL3DDLRAWREADER_H

#include "AliHLTRootTypes.h"

// see description in upcoming ALICE note
// by D.Favretto and A.K.Mohanty
struct AliHLTDDLMiniHeader 
{
  UInt_t    fSize;           // size
  UChar_t   fDetectorID;     // detector ID
  UChar_t   fMagicWord[3];   // magic word
  UChar_t   fVersion;        // version
  UChar_t   fCompressionFlag;// compression flag
  UShort_t  fDDLID;          // DDL ID
};

class AliHLTDDLRawReader 
{
  public :
    AliHLTDDLRawReader();
    virtual ~AliHLTDDLRawReader();

    void Select(Int_t detectorID, Int_t minDDLID = -1, Int_t maxDDLID = -1);

    Int_t     GetDataSize()   const {return fMiniHeader->fSize;};
    Int_t     GetDetectorID() const {return fMiniHeader->fDetectorID;};
    Int_t     GetDDLID()      const {return fMiniHeader->fDDLID;};
    Int_t     GetVersion()    const {return fMiniHeader->fVersion;};
    Bool_t    IsCompressed()  const {return fMiniHeader->fCompressionFlag != 0;};

    virtual Bool_t   ReadMiniHeader() = 0;
    virtual Bool_t   ReadNextData(UChar_t*& data) = 0;
    virtual Bool_t   ReadNextInt(UInt_t& data);
    virtual Bool_t   ReadNextShort(UShort_t& data);
    virtual Bool_t   ReadNextChar(UChar_t& data);

    virtual Bool_t   Reset() = 0;

  protected :
    Bool_t           IsSelected() const;

    Bool_t           CheckMiniHeader() const;
    virtual Bool_t   ReadNext(UChar_t* data, Int_t size) = 0;

    AliHLTDDLMiniHeader*  fMiniHeader;  // current mini header
    Int_t                fCount;       // counter of bytes to be read for current DDL

    Int_t            fSelectDetectorID;  // id of selected detector (<0 = no selection)
    Int_t            fSelectMinDDLID;    // minimal index of selected DDLs (<0 = no selection)
    Int_t            fSelectMaxDDLID;    // maximal index of selected DDLs (<0 = no selection)

    ClassDef(AliHLTDDLRawReader,1) //AliHLTDDLRawReader
};

typedef AliHLTDDLRawReader AliL3DDLRawReader; // for backward compatibility

#endif
