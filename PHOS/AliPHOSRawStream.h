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
  
  Int_t            GetColumn() const {return fPad;}
  Int_t            GetModule() const {return fSector;}
  Int_t            GetPrevColumn() const {return fPrevPad;}
  Int_t            GetPrevModule() const {return fPrevSector;}
  Int_t            GetPrevRow() const {return fPrevRow;}
  Int_t            GetRow() const {return fRow;}
  Int_t            GetSignal() const {return fSignal;}
  Int_t            GetTime() const {return fTime;}
  Bool_t           IsNewColumn() const {return (GetColumn() != GetPrevColumn()) || IsNewRow();}
  Bool_t           IsNewModule() const {return GetModule() != GetPrevModule();}
  Bool_t           IsNewRow() const {return (GetRow() != GetPrevRow()) || IsNewModule();}
  
  ClassDef(AliPHOSRawStream, 0)   // class for reading PHOS raw digits
    };

#endif

