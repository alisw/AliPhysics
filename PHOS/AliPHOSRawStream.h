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

#include "AliAltroRawStream.h"

class AliRawReader;


class AliPHOSRawStream: public AliAltroRawStream {
  public :
    AliPHOSRawStream(AliRawReader* rawReader);

    Int_t            GetModule() const {return fSector;};
    Int_t            GetPrevModule() const {return fPrevSector;};
    Bool_t           IsNewModule() const {return fSector != fPrevSector;};
    Int_t            GetRow() const {return fRow;};
    Int_t            GetPrevRow() const {return fPrevRow;};
    Bool_t           IsNewRow() const {return (fRow != fPrevRow) || IsNewModule();};
    Int_t            GetColumn() const {return fPad;};
    Int_t            GetPrevColumn() const {return fPrevPad;};
    Bool_t           IsNewColumn() const {return (fPad != fPrevPad) || IsNewRow();};
    Int_t            GetTime() const {return fTime;};
    Int_t            GetSignal() const {return fSignal;};

    ClassDef(AliPHOSRawStream, 0)   // class for reading PHOS raw digits
};

#endif
