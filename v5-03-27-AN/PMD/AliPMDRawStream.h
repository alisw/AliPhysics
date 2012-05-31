#ifndef ALIPMDRAWSTREAM_H
#define ALIPMDRAWSTREAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
///
/// This class provides access to PMD digits in raw data.
///
///////////////////////////////////////////////////////////////////////////////

#include <TObject.h>

class AliRawReader;
class AliCDBManager;
class AliCDBStorage;
class AliCDBEntry;
class AliPMDMappingData;


class AliPMDRawStream: public TObject {
  public :
    AliPMDRawStream(AliRawReader* rawReader);
    virtual ~AliPMDRawStream();

    Int_t DdlData(TObjArray *pmdddlcont);

    enum {kDDLOffset = 0xC00};      // offset for DDL numbers

    enum EPMDRawStreamError {
      kDDLIndexMismatch = 1,
      kNoMappingFile = 2,
      kParityError = 3
    };

  private :
    AliPMDRawStream(const AliPMDRawStream& stream);
    AliPMDRawStream& operator = (const AliPMDRawStream& stream);

    void             GetRowCol(Int_t imodule, Int_t pbusid, 
			       UInt_t mcmno, UInt_t chno,
			       Int_t startRowBus[], Int_t endRowBus[],
			       Int_t startColBus[], Int_t endColBus[],
			       Int_t &row, Int_t &col) const;
    void             ConvertDDL2SMN(Int_t iddl, Int_t imodule,
				    Int_t &smn, Int_t &detector) const;
    void             TransformH2S(Int_t smn, Int_t &row, Int_t &col) const;
    Int_t            ComputeParity(UInt_t data1);
    UInt_t           GetNextWord();
    void             Ddl0Mapping(Int_t moduleNo[],    Int_t mcmperBus[],
				 Int_t startRowBus[], Int_t endRowBus[],
				 Int_t startColBus[], Int_t endColBus[]);
    void             Ddl1Mapping(Int_t moduleNo[],    Int_t mcmperBus[],
				 Int_t startRowBus[], Int_t endRowBus[],
				 Int_t startColBus[], Int_t endColBus[]);
    void             Ddl2Mapping(Int_t moduleNo[],    Int_t mcmperBus[],
				 Int_t startRowBus[], Int_t endRowBus[],
				 Int_t startColBus[], Int_t endColBus[]);
    void             Ddl3Mapping(Int_t moduleNo[],    Int_t mcmperBus[],
				 Int_t startRowBus[], Int_t endRowBus[],
				 Int_t startColBus[], Int_t endColBus[]);
    void             Ddl4Mapping(Int_t moduleNo[],    Int_t mcmperBus[],
				 Int_t startRowBus[], Int_t endRowBus[],
				 Int_t startColBus[], Int_t endColBus[]);
    void             Ddl5Mapping(Int_t moduleNo[],    Int_t mcmperBus[],
				 Int_t startRowBus[], Int_t endRowBus[],
				 Int_t startColBus[], Int_t endColBus[]);

    AliPMDMappingData *GetMappingData() const;

    AliRawReader*    fRawReader;    // object for reading the raw data
    UChar_t*         fData;         // pointer to the data
    Int_t            fPosition;

    AliPMDMappingData  *fMapData;   //! Mapping data

    ClassDef(AliPMDRawStream, 8)    // class for reading PMD raw digits
};

#endif
