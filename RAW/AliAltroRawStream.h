#ifndef ALIALTRORAWSTREAM_H
#define ALIALTRORAWSTREAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
///
/// This is a base class for reading raw data digits in Altro format
/// The class is able to read both old and new RCU trailer formats
/// Switch between formats is done automatically using the last payload word.
/// In case the Common Data Header is 7 32-bit words long, one
/// can use the fIsShortDataHeader flag.
///
///////////////////////////////////////////////////////////////////////////////

#include <TObject.h>

class AliRawReader;

class AliAltroRawStream: public TObject {
  public :
    AliAltroRawStream(AliRawReader* rawReader);
    virtual ~AliAltroRawStream();

    virtual void             Reset();
    virtual Bool_t           Next();

    virtual Bool_t NextDDL(UChar_t* data = NULL);              // Iterate over DDLs/RCUs
    virtual Bool_t NextChannel();                              // Iterate over altro channels
    virtual Bool_t NextBunch(UShort_t *bunchData,
			     Int_t &bunchLength,
			     Int_t &startTimeBin);             // Iterate over altro bunches

    Int_t GetDDLNumber()  const { return fDDLNumber; }  // Provide current DDL number
    Int_t GetPrevDDLNumber() const { return fPrevDDLNumber; }// Provide previous DDL number
    Bool_t  IsNewDDLNumber() const {return (fDDLNumber != fPrevDDLNumber);};
    Int_t GetRCUId()      const { return fRCUId; }      // Provide current RCU identifier
    Int_t GetPrevRCUId()  const { return fPrevRCUId; }  // Provide previous RCU identifier
    Bool_t  IsNewRCUId() const {return (fRCUId != fPrevRCUId);};
    Int_t GetHWAddress()  const { return fHWAddress; }  // Provide current hardware address
    Int_t GetPrevHWAddress() const { return fPrevHWAddress; }  // Provide previous hardware address
    Bool_t  IsNewHWAddress() const {return (fHWAddress != fPrevHWAddress) || IsNewDDLNumber();};
    Int_t GetTime()       const { return fTime; }       // Provide index of current time bin
    Int_t GetPrevTime()   const { return fPrevTime; }   // Provide index of previous time bin
    Bool_t  IsNewTime()   const {return (fTime != fPrevTime) || IsNewHWAddress();};
    Int_t GetSignal()     const { return fSignal; }     // Provide signal in ADC counts
    Int_t GetTimeLength() const { return fTimeBunch; }  // Provide total length of current time bunch

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
    UShort_t GetNChAddrMismatch() const { return fERRREG3; }
    UShort_t GetNChLengthMismatch() const { return fERRREG4; }

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
 
    void SelectRawData(Int_t detId);                           // Select raw data for specific detector id
    void SelectRawData(const char *detName);                   // Select raw data for specific detector name

    void  SetShortDataHeader(Bool_t flag) { fIsShortDataHeader = flag; } // Specify whenever to assume or not a short CDH format

    void PrintDebug() const; // Print debug information in case of decoding errors
    void AddMappingErrorLog(const char *message = NULL);

    enum EAltroRawStreamError {
      kRCUTrailerSizeErr = 1,
      kAltroTrailerErr = 2,
      kBunchLengthReadErr = 3,
      kTimeBinReadErr = 4,
      kAmplitudeReadErr = 5,
      k32bitWordReadErr = 6,
      kBadAltroMapping = 7,
      kRCUTrailerErr = 8
    };

    AliAltroRawStream& operator = (const AliAltroRawStream& stream);
    AliAltroRawStream(const AliAltroRawStream& stream);

  protected:

    Bool_t           fIsShortDataHeader; // flag used to select between normal and short CDH format

  private :

    UShort_t         GetNextWord();
    Bool_t           ReadTrailer();
    void             ReadBunch();
    void             ReadAmplitude();
    Int_t            GetPosition();
    UInt_t           Get32bitWord(Int_t &index);
    Int_t            ReadRCUTrailer(Int_t &index, Int_t trailerSize);

    Int_t            fDDLNumber;    // index of current DDL number
    Int_t            fPrevDDLNumber;// index of previous DDL number
    Int_t            fRCUId;        // current RCU identifier
    Int_t            fPrevRCUId;    // previous RCU identifier
    Short_t          fHWAddress;    // current hardware address
    Short_t          fPrevHWAddress;// previous hardware address
    Int_t            fTime;         // index of current time bin
    Int_t            fPrevTime;     // index of previous time bin
    Int_t            fSignal;       // signal in ADC counts
    Int_t            fTimeBunch;    // total length of the current time bunch

    AliRawReader*    fRawReader;    // object for reading the raw data


    UChar_t*         fData;         // raw data

    Int_t            fPosition;     // current (10 bit) position in fData
    Int_t            fCount;        // counter of words to be read for current trailer
    Int_t            fBunchLength;  // remaining number of signal bins in the current bunch

    UChar_t*         fRCUTrailerData; // pointer to RCU trailer data
    Int_t            fRCUTrailerSize; // size of RCU trailer data in bytes

    // RCU trailer contents
    UInt_t           fFECERRA;      // contains errors related to ALTROBUS transactions
    UInt_t           fFECERRB;      // contains errors related to ALTROBUS transactions
    UShort_t         fERRREG2;      // contains errors related to ALTROBUS transactions or trailer of ALTRO channel block
    UShort_t         fERRREG3;      // contains number of altro channels skipped due to an address mismatch 
    UShort_t         fERRREG4;      // contains number of altro channels skipped due to a block length mismatch 
    UShort_t         fActiveFECsA;  // bit pattern of active FECs in branch A
    UShort_t         fActiveFECsB;  // bit pattern of active FECs in branch B
    UInt_t           fAltroCFG1;    // ALTROCFG1 register
    UInt_t           fAltroCFG2;    // ALTROCFG2 and ALTROIF registers

    ClassDef(AliAltroRawStream, 0)  // base class for reading Altro raw digits
};

#endif
