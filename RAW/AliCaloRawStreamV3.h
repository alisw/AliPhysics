#ifndef ALICALORAWSTREAMV3_H
#define ALICALORAWSTREAMV3_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: $ */

///////////////////////////////////////////////////////////////////////////////
///
/// This class provides access to Calo digits in raw data.
///
/// Yuri Kharlov. 23 June 2009
/// 
///////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---
#include "TString.h"

// --- AliRoot header files ---
#include "AliAltroRawStreamV3.h"
class AliRawReader;
class AliAltroMapping;

class AliCaloRawStreamV3: public AliAltroRawStreamV3 {

public :
  AliCaloRawStreamV3(AliRawReader* rawReader, TString calo, AliAltroMapping **mapping = NULL);
  virtual ~AliCaloRawStreamV3();
  
  virtual void             Reset();
  virtual Bool_t           NextChannel();
  
  Int_t   GetModule()      const {return fModule;}
  Int_t   GetRow()         const {return fRow   ;} // EMCAL notation
  Int_t   GetColumn()      const {return fColumn;} // EMCAL notation
  Int_t   GetCellX()       const {return fRow   ;} // PHOS notation
  Int_t   GetCellZ()       const {return fColumn;} // PHOS notation
  Int_t   GetNRCU()        const {return fNRCU  ;}
  Int_t   GetNSides()      const {return fNSides;}
  TString GetCalorimeter() const {return fCalo  ;}

  enum EAliCaloFlag { kLowGain=0, kHighGain=1, kTRUData=2, kLEDMonData=3 };
  Bool_t  IsLowGain()      const {return (fCaloFlag == kLowGain)   ;}
  Bool_t  IsHighGain()     const {return (fCaloFlag == kHighGain)  ;}
  Bool_t  IsTRUData()      const {return (fCaloFlag == kTRUData)   ;}
  Bool_t  IsLEDMonData()   const {return (fCaloFlag == kLEDMonData);} 

  Int_t   GetCaloFlag() const { return fCaloFlag; } 
  Int_t   GetFilter() const { return fFilter; } 

  void SkipData(EAliCaloFlag caloFlag=kLEDMonData) 
    { fFilter |= (1<<caloFlag); }

protected:

  AliCaloRawStreamV3& operator = (const AliCaloRawStreamV3& stream);
  AliCaloRawStreamV3(const AliCaloRawStreamV3& stream);

  virtual void ApplyAltroMapping();

  Int_t            fModule;   // index of current module
  Int_t            fRow;      // index of current row
  Int_t            fColumn;   // index of current column
  Int_t            fCaloFlag; // low (0) or (1) high gain; see enum EAliCaloFlag above
  Int_t            fFilter;   // default 0 = let everything through
  Int_t            fNModules; // number of (super)modules
  Int_t            fNRCU;     // number of RCU per (super)module
  Int_t            fNSides;   // Division of EMCal in "A" "C" sides
  TString          fCalo;     // Calorimeter name
  Bool_t           fExternalMapping; // use external mapping or create a default one
  AliAltroMapping *fMapping[20];     // pointers to ALTRO mapping

  ClassDef(AliCaloRawStreamV3, 2)   // class for reading PHOS/EMCAL raw digits

};

#endif
