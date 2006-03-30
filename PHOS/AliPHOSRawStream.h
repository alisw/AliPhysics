#ifndef ALIPHOSRAWSTREAM_H
#define ALIPHOSRAWSTREAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
///
/// This class provides access to PHOS digits in raw data.
///
///////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---

// --- AliRoot header files ---
#include "AliAltroRawStream.h"
class AliRawReader;


class AliPHOSRawStream: public AliAltroRawStream {

public :
  AliPHOSRawStream(AliRawReader* rawReader);
  virtual ~AliPHOSRawStream();
 
  virtual void             Reset();
  virtual Bool_t           Next();
  
  Int_t            GetColumn() const {return fColumn;}
  Int_t            GetModule() const {return fModule;}
  Int_t            GetPrevColumn() const {return fPrevColumn;}
  Int_t            GetPrevModule() const {return fPrevModule;}
  Int_t            GetPrevRow() const {return fPrevRow;}
  Int_t            GetRow() const {return fRow;}
  Bool_t           IsNewColumn() const {return (GetColumn() != GetPrevColumn()) || IsNewRow();}
  Bool_t           IsNewModule() const {return GetModule() != GetPrevModule();}
  Bool_t           IsNewRow() const {return (GetRow() != GetPrevRow()) || IsNewModule();}

protected:
    AliPHOSRawStream(const AliPHOSRawStream& stream);
    AliPHOSRawStream& operator = (const AliPHOSRawStream& stream);

    virtual void ApplyAltroMapping();

    Int_t            fModule;       // index of current module
    Int_t            fPrevModule;   // index of previous module
    Int_t            fRow;          // index of current row
    Int_t            fPrevRow;      // index of previous row
    Int_t            fColumn;       // index of current column
    Int_t            fPrevColumn;   // index of previous column
  
  ClassDef(AliPHOSRawStream, 0)   // class for reading PHOS raw digits
    };

#endif

