#ifndef ALIEMCALRAWSTREAM_H
#define ALIEMCALRAWSTREAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
///
/// This class provides access to EMCAL digits in raw data.
///
///////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---

// --- AliRoot header files ---
#include "AliAltroRawStream.h"
class AliRawReader;


class AliEMCALRawStream: public AliAltroRawStream {

public :
  
  AliEMCALRawStream(AliRawReader* rawReader);
  
  Int_t            GetId() const {return fPad;};
  Int_t            GetModule() const {return fSector;}
  Int_t            GetPrevId() const {return fPrevPad;};
  Int_t            GetSignal() const {return fSignal;};
  Int_t            GetTime() const {return fTime;};
  Bool_t           IsNewId() const {return (GetId() != GetPrevId());};
  
  ClassDef(AliEMCALRawStream, 0)   // class for reading EMCAL raw digits
    };

#endif
