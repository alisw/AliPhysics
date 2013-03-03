#ifndef ALITPCRAWSTREAMV3_H
#define ALITPCRAWSTREAMV3_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////////////////////////////////
///
/// This class provides access to TPC digits in raw data.
///
///////////////////////////////////////////////////////////////////////////////

#include "AliAltroRawStreamV3.h"

class AliRawReader;
class AliAltroMapping;

class AliTPCRawStreamV3: public AliAltroRawStreamV3 {
  public :
    AliTPCRawStreamV3(AliRawReader* rawReader, AliAltroMapping **mapping = NULL);
    virtual ~AliTPCRawStreamV3();

    virtual void             Reset();
    virtual Bool_t           NextChannel();
    virtual Bool_t           NextDDL();
  
    inline Int_t GetSector()     const { return fSector; }     // Provide index of current sector
    inline Int_t GetPrevSector() const { return fPrevSector; } // Provide index of previous sector
    inline Bool_t  IsNewSector() const {return fSector != fPrevSector;};
    inline Int_t GetRow()        const { return fRow; }        // Provide index of current row
    inline Int_t GetPrevRow()    const { return fPrevRow; }    // Provide index of previous row
    inline Bool_t  IsNewRow()    const {return (fRow != fPrevRow) || IsNewSector();};
    inline Int_t GetPad()        const { return fPad; }        // Provide index of current pad
    inline Int_t GetPrevPad()    const { return fPrevPad; }    // Provide index of previous pad
    inline Bool_t  IsNewPad()    const {return (fPad != fPrevPad) || IsNewRow();};
    inline Int_t GetPatchIndex() const { return fPatchIndex; }        // Provide index of current patch


  protected :
    AliTPCRawStreamV3& operator = (const AliTPCRawStreamV3& stream);
    AliTPCRawStreamV3(const AliTPCRawStreamV3& stream);

    virtual void ApplyAltroMapping();

    Int_t            fSector;       // index of current sector
    Int_t            fPrevSector;   // index of previous sector
    Int_t            fRow;          // index of current row
    Int_t            fPrevRow;      // index of previous row
    Int_t            fPad;          // index of current pad
    Int_t            fPrevPad;      // index of previous pad
    Int_t            fPatchIndex;   // current patch

    AliAltroMapping *fMapping[6];   // Pointers to ALTRO mapping
    Bool_t           fIsMapOwner;   // does object own its mappings?

    ClassDef(AliTPCRawStreamV3, 0)    // base class for reading TPC raw digits using the fast algorithm
};

#endif
