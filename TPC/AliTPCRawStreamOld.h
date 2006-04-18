#ifndef ALITPCRAWSTREAMOLD_H
#define ALITPCRAWSTREAMOLD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////////////////////////////////
///
/// This class provides access to TPC digits in raw data.
///
///////////////////////////////////////////////////////////////////////////////

#include "AliAltroRawStreamOld.h"

class AliRawReader;
class AliAltroMapping;

class AliTPCRawStreamOld: public AliAltroRawStreamOld {
  public :
    AliTPCRawStreamOld(AliRawReader* rawReader);
    virtual ~AliTPCRawStreamOld();

    virtual void             Reset();
    virtual Bool_t           Next();

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
    AliTPCRawStreamOld(const AliTPCRawStreamOld& stream);
    AliTPCRawStreamOld& operator = (const AliTPCRawStreamOld& stream);

    virtual void ApplyAltroMapping();

    Int_t            fSector;       // index of current sector
    Int_t            fPrevSector;   // index of previous sector
    Int_t            fRow;          // index of current row
    Int_t            fPrevRow;      // index of previous row
    Int_t            fPad;          // index of current pad
    Int_t            fPrevPad;      // index of previous pad

    AliAltroMapping *fMapping[6];   // Pointers to ALTRO mapping

    ClassDef(AliTPCRawStreamOld, 0)    // base class for reading TPC raw digits
};

#endif
