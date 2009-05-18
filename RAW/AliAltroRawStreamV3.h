#ifndef ALIALTRORAWSTREAMV3_H
#define ALIALTRORAWSTREAMV3_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////////////////////////////////
///
/// This is a base class for reading raw data digits in Altro format.
/// The class is able to read the RCU v3 and above formats.
/// The main difference between the format V3 and older ones is in
/// the coding of the 10-bit Altro payload words. In V3 3 10-bit words
/// are coded in one 32-bit word. The bits 30 and 31 are used to identify
/// the payload, altro header and RCU trailer contents.
///
///
/// cvetan.cheshkov@cern.ch 1/04/2009
///////////////////////////////////////////////////////////////////////////////

#include <TObject.h>

class AliRawReader;

class AliAltroRawStreamV3: public TObject {
  public :
    AliAltroRawStreamV3(AliRawReader* rawReader);
    virtual ~AliAltroRawStreamV3();

    AliAltroRawStreamV3& operator = (const AliAltroRawStreamV3& stream);
    AliAltroRawStreamV3(const AliAltroRawStreamV3& stream);

    void SelectRawData(Int_t detId);                           // Select raw data for specific detector id
    void SelectRawData(const char *detName);                   // Select raw data for specific detector name

    virtual void   Reset();                                    // Reset the raw-stream object

    virtual Bool_t NextDDL();                                  // Iterate over DDLs/RCUs
    virtual Bool_t NextChannel();                              // Iterate over altro channels
    virtual Bool_t NextBunch();                                // Iterate over altro bunches

    Int_t  GetDDLNumber()      const { return fDDLNumber; }    // Provide current DDL number
    Int_t  GetHWAddress()      const { return fHWAddress; }    // Provide current hardware address
    Int_t  GetRCUId()          const { return fRCUId; }        // Provide current RCU identifier
 
    UInt_t GetStartTimeBin()   const { return fStartTimeBin; } // Provide the index if the first time-bin in current bunch
    UInt_t GetEndTimeBin()     const { return fStartTimeBin-fBunchLength+1; } // Provide the index of the last time-bin in current bunch
    Int_t  GetBunchLength()    const { return fBunchLength; }  // Provide the current bunch length
    const UShort_t* GetSignals() const { return fBunchDataPointer; }  // Provide access to altro data itself
    Bool_t IsChannelBad()      const { return fBadChannel; }   // Is the channel data bad or not


    Int_t GetBranch()     const; // Provide the branch index for the current hardware address
    Int_t GetFEC()        const; // Provide the front-end card index for the current hardware address
    Int_t GetAltro()      const; // Provide the altro chip index for the current hardware address
    Int_t GetChannel()    const; // Provide the channel index for the current hardware address

    Bool_t  GetRCUTrailerData(UChar_t*& data) const;              // Provide a pointer to RCU trailer
    Int_t   GetRCUTrailerSize() const { return fRCUTrailerSize; } // Provide size of RCU trailer

    // RCU trailer related getters
    UInt_t  GetFECERRA() const { return fFECERRA; }
    UInt_t  GetFECERRB() const { return fFECERRB; }
    UShort_t GetERRREG2() const { return fERRREG2; }
    UShort_t GetNChAddrMismatch() const { return (fERRREG3 & 0xFFF); }
    UShort_t GetNChLengthMismatch() const { return ((fERRREG3 >> 12) & 0x1FFF); }

    UShort_t GetActiveFECsA() const { return fActiveFECsA; }
    UShort_t GetActiveFECsB() const { return fActiveFECsB; }

    UInt_t  GetAltroCFG1() const { return fAltroCFG1; }
    UChar_t GetBaselineCorr() const { return fAltroCFG1 & 0xF; }
    Bool_t  GetPolarity() const { return (fAltroCFG1 >> 4) & 0x1; }
    UChar_t GetNPresamples() const  { return (fAltroCFG1 >> 5) & 0x3; }
    UChar_t GetNPostsamples() const { return (fAltroCFG1 >> 7) & 0xF; }
    Bool_t  GetSecondBaselineCorr() const { return (fAltroCFG1 >> 11) & 0x1; }
    UChar_t GetGlitchFilter() const { return (fAltroCFG1 >> 12) & 0x3; }
    UChar_t GetNNonZSPostsamples() const { return (fAltroCFG1 >> 14) & 0x7; }
    UChar_t GetNNonZSPresamples() const  { return (fAltroCFG1 >> 17) & 0x3; }
    Bool_t  GetZeroSupp() const          { return (fAltroCFG1 >> 19) & 0x1; }
    
    UInt_t   GetAltroCFG2() const { return fAltroCFG2; }
    Bool_t   GetNAltroBuffers() const     { return (fAltroCFG2 >> 24) & 0x1; }
    UChar_t  GetNPretriggerSamples() const{ return (fAltroCFG2 >> 20) & 0xF; }
    UShort_t GetNSamplesPerCh() const     { return (fAltroCFG2 >> 10) & 0x3FF; }
    Bool_t   GetSparseRO() const          { return (fAltroCFG2 >> 9) & 0x1; }
    Double_t GetTSample() const;
    Double_t GetL1Phase() const;
    void     PrintRCUTrailer() const;
 
    void  SetShortDataHeader(Bool_t flag) { fIsShortDataHeader = flag; } // Specify whenever to assume or not a short CDH format

    enum EAltroRawStreamV3Error {
      kRCUTrailerErr = 1,
      kRCUVerErr = 2,
      kRCUTrailerSizeErr = 3,
      kAltroBunchHeadErr = 4,
      kBunchLengthErr = 5,
      kAltroPayloadErr = 6,
      kBadAltroMapping = 7
    };

    enum {kMaxNTimeBins = 1024};

  protected:

    void             AddMappingErrorLog(const char *message = NULL);

    Bool_t           fIsShortDataHeader; // flag used to select between normal and short CDH format

  private:

    UInt_t           Get32bitWord(Int_t index) const;
    Bool_t           ReadRCUTrailer(UChar_t rcuVer);

    Int_t            fDDLNumber;    // index of current DDL number
    Int_t            fRCUId;        // current RCU identifier
    Short_t          fHWAddress;    // current hardware address

    AliRawReader*    fRawReader;    // object for reading the raw data

    UChar_t*         fData;         // raw data

    Int_t            fPosition;     // current position (32-bit words) in fData
    Int_t            fCount;        // 
    Int_t            fStartTimeBin; //
    Int_t            fBunchLength;  //

    Bool_t           fBadChannel;   //
    Int_t            fPayloadSize;  //

    UShort_t         fBunchData[kMaxNTimeBins];    // cache for the decoded altro payload
    UShort_t*        fBunchDataPointer;            // pointer to the current bunch samples
    Int_t            fBunchDataIndex;              // current position in the payload

    UChar_t*         fRCUTrailerData; // pointer to RCU trailer data
    Int_t            fRCUTrailerSize; // size of RCU trailer data in bytes

    // RCU trailer contents
    UInt_t           fFECERRA;      // contains errors related to ALTROBUS transactions
    UInt_t           fFECERRB;      // contains errors related to ALTROBUS transactions
    UShort_t         fERRREG2;      // contains errors related to ALTROBUS transactions or trailer of ALTRO channel block
    UInt_t           fERRREG3;      // contains number of altro channels skipped due to an address mismatch 
    UShort_t         fActiveFECsA;  // bit pattern of active FECs in branch A
    UShort_t         fActiveFECsB;  // bit pattern of active FECs in branch B
    UInt_t           fAltroCFG1;    // ALTROCFG1 register
    UInt_t           fAltroCFG2;    // ALTROCFG2 and ALTROIF registers

    ClassDef(AliAltroRawStreamV3, 0)  // base class for reading Altro raw digits
};

#endif
