// @(#) $Id$

#ifndef ALIHLTTPCDDLRAWREADER_H
#define ALIHLTTPCDDLRAWREADER_H

// see description in upcoming ALICE note
// by D.Favretto and A.K.Mohanty
struct AliHLTTPCDDLMiniHeader 
{
  UInt_t    fSize;
  UChar_t   fDetectorID;
  UChar_t   fMagicWord[3];
  UChar_t   fVersion;
  UChar_t   fCompressionFlag;
  UShort_t  fDDLID;
};

class AliHLTTPCDDLRawReader 
{
  public :
    AliHLTTPCDDLRawReader();
    virtual ~AliHLTTPCDDLRawReader();

    void Select(Int_t detectorID, Int_t minDDLID = -1, Int_t maxDDLID = -1);

    inline Int_t     GetDataSize()   const {return fMiniHeader->fSize;};
    inline Int_t     GetDetectorID() const {return fMiniHeader->fDetectorID;};
    inline Int_t     GetDDLID()      const {return fMiniHeader->fDDLID;};
    inline Int_t     GetVersion()    const {return fMiniHeader->fVersion;};
    inline Bool_t    IsCompressed()  const {return fMiniHeader->fCompressionFlag != 0;};

    virtual Bool_t   ReadMiniHeader() = 0;
    virtual Bool_t   ReadNextData(UChar_t*& data) = 0;
    virtual Bool_t   ReadNextInt(UInt_t& data);
    virtual Bool_t   ReadNextShort(UShort_t& data);
    virtual Bool_t   ReadNextChar(UChar_t& data);

    virtual Bool_t   Reset() = 0;

  protected :
    Bool_t           IsSelected();

    Bool_t           CheckMiniHeader();
    virtual Bool_t   ReadNext(UChar_t* data, Int_t size) = 0;

    AliHLTTPCDDLMiniHeader*  fMiniHeader;  // current mini header
    Int_t                fCount;       // counter of bytes to be read for current DDL

    Int_t            fSelectDetectorID;  // id of selected detector (<0 = no selection)
    Int_t            fSelectMinDDLID;    // minimal index of selected DDLs (<0 = no selection)
    Int_t            fSelectMaxDDLID;    // maximal index of selected DDLs (<0 = no selection)

    ClassDef(AliHLTTPCDDLRawReader,1) //AliHLTTPCDDLRawReader
};

#endif
