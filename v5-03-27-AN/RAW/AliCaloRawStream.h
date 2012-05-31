#ifndef ALICALORAWSTREAM_H
#define ALICALORAWSTREAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
///
/// This class provides access to Calo digits in raw data.
///
///Modification: Class exported from PHOS to be used by EMCAL and PHOS
///November 2006 Gustavo Conesa Balbastre 
///////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---
#include "TString.h"

// --- AliRoot header files ---
#include "AliAltroRawStream.h"
class AliRawReader;
class AliAltroMapping;

class AliCaloRawStream: public AliAltroRawStream {

public :
  AliCaloRawStream(AliRawReader* rawReader,  TString calo, AliAltroMapping **mapping = NULL);
  virtual ~AliCaloRawStream();
 
  virtual void             Reset();
  virtual Bool_t           Next();
  
  Int_t            GetModule()     const {return fModule;}
  Int_t            GetRow()        const {return fRow;}
  Int_t            GetColumn()     const {return fColumn;}
  Int_t            GetPrevModule() const {return fPrevModule;}
  Int_t            GetPrevRow()    const {return fPrevRow;}
  Int_t            GetPrevColumn() const {return fPrevColumn;}
  Bool_t           IsNewModule()   const {return GetModule() != GetPrevModule();}
  Bool_t           IsNewRow()      const {return (GetRow() != GetPrevRow()) || IsNewModule();}
  Bool_t           IsNewColumn()   const {return (GetColumn() != GetPrevColumn()) || IsNewRow();}
  Int_t            GetNRCU() const {return fNRCU;}
  Int_t            GetNSides() const {return fNSides;}
  TString          GetCalorimeter() const {return fCalo;}
  enum EAliCaloFlag { kLowGain=0, kHighGain=1, kTRUData=2, kLEDMonData=3 };
  Bool_t           IsLowGain()     const {return (fCaloFlag == kLowGain);}
  Bool_t           IsHighGain()    const {return (fCaloFlag == kHighGain);}
  Bool_t           IsTRUData()     const {return (fCaloFlag == kTRUData);}
  Bool_t           IsLEDMonData()  const {return (fCaloFlag == kLEDMonData);} 

  Int_t GetCaloFlag() const { return fCaloFlag; } 
  Int_t GetFilter() const { return fFilter; } 

  void SkipData(EAliCaloFlag caloFlag=kLEDMonData) 
    { fFilter |= (1<<caloFlag); }

protected:

  AliCaloRawStream(const AliCaloRawStream& stream);
  AliCaloRawStream& operator = (const AliCaloRawStream& stream);

  virtual void ApplyAltroMapping();

  Int_t            fModule;       // index of current module
  Int_t            fPrevModule;   // index of previous module
  Int_t            fRow;          // index of current row
  Int_t            fPrevRow;      // index of previous row
  Int_t            fColumn;       // index of current column
  Int_t            fPrevColumn;   // index of previous column
  Int_t            fCaloFlag;     // low (0) or (1) high gain; see enum EAliCaloFlag above
  Int_t            fFilter; // default 0 = let everything through
  Int_t            fNRCU;   // number of RCU per (super)module
  Int_t            fNSides;   // Division of EMCal in "A" "C" sides
  TString          fCalo; // Calorimeter name
  Bool_t           fExternalMapping;   // use external mapping or create a default one
  AliAltroMapping *fMapping[4];   // pointers to ALTRO mapping

  ClassDef(AliCaloRawStream, 1)   // class for reading PHOS/EMCAL raw digits

};

#endif

