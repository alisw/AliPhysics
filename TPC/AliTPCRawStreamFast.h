#ifndef ALITPCRAWSTREAMFAST_H
#define ALITPCRAWSTREAMFAST_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////////////////////////////////
///
/// This class provides access to TPC digits in raw data.
///
///////////////////////////////////////////////////////////////////////////////

#include "AliAltroRawStreamFast.h"

class AliRawReader;
class AliAltroMapping;

class AliTPCRawStreamFast: public AliAltroRawStreamFast {
  public :
    AliTPCRawStreamFast(AliRawReader* rawReader);
    virtual ~AliTPCRawStreamFast();

    virtual void             Reset();
    virtual Bool_t           NextChannel();

    inline Int_t GetSector()     const { return fSector; }     // Provide index of current sector
    inline Int_t GetPrevSector() const { return fPrevSector; } // Provide index of previous sector
    inline Bool_t  IsNewSector() const {return fSector != fPrevSector;};
    inline Int_t GetRow()        const { return fRow; }        // Provide index of current row
    inline Int_t GetPrevRow()    const { return fPrevRow; }    // Provide index of previous row
    inline Bool_t  IsNewRow()    const {return (fRow != fPrevRow) || IsNewSector();};
    inline Int_t GetPad()        const { return fPad; }        // Provide index of current pad
    inline Int_t GetPrevPad()    const { return fPrevPad; }    // Provide index of previous pad
    inline Bool_t  IsNewPad()    const {return (fPad != fPrevPad) || IsNewRow();};


  protected :
    AliTPCRawStreamFast& operator = (const AliTPCRawStreamFast& stream);
    AliTPCRawStreamFast(const AliTPCRawStreamFast& stream);

    virtual void ApplyAltroMapping();

    Int_t            fSector;       // index of current sector
    Int_t            fPrevSector;   // index of previous sector
    Int_t            fRow;          // index of current row
    Int_t            fPrevRow;      // index of previous row
    Int_t            fPad;          // index of current pad
    Int_t            fPrevPad;      // index of previous pad

    AliAltroMapping *fMapping[6];   // Pointers to ALTRO mapping
    Bool_t           fIsMapOwner;   // does object own its mappings?

    ClassDef(AliTPCRawStreamFast, 0)    // base class for reading TPC raw digits using the fast algorithm
};

#endif
