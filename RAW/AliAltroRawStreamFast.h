#ifndef ALIALTRORAWSTREAMFAST_H
#define ALIALTRORAWSTREAMFAST_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////////////////////////////////
///
/// This is a base class for reading raw data in Altro format.
/// It makes use of the fast altro decoder written by Per Thomas.
/// More information can be found in AliAltroDecoder class.
///
///////////////////////////////////////////////////////////////////////////////

#include <TObject.h>
#include <Riostream.h>

#include "AliRawReader.h"
#include "AliAltroDecoder.h"
#include "AliAltroData.h"
#include "AliAltroBunch.h"

class AliAltroRawStreamFast: public TObject {
  public :
    AliAltroRawStreamFast(AliRawReader* rawReader);
    virtual ~AliAltroRawStreamFast();

    void SelectRawData(Int_t detId);                           // Select raw data for specific detector id
    void SelectRawData(const char *detName);                   // Select raw data for specific detector name

    virtual void   Reset();
    virtual Bool_t NextDDL();
    virtual Bool_t NextChannel();
    virtual Bool_t NextBunch();

    Int_t  GetDDLNumberOld()   const { return fRawReader->GetDDLID();         } // Provide current DDL number
    Int_t  GetDDLNumber()      const { return fDDLNumber;                     } // Provide current DDL number
    Int_t  GetPrevDDLNumber()  const { return fPrevDDLNumber;                 } // Provide previous DDL number
    Bool_t IsNewDDLNumber()    const { return (fDDLNumber != fPrevDDLNumber); }

    Int_t  GetHWAddress()      const { return fData.GetHadd();                } // Provide current hardware address
    Int_t  GetPrevHWAddress()  const { return fData.GetPrevHadd();            } // Provide previous hardware address
    Bool_t IsNewHWAddress()    const { return fData.IsNewHadd() || IsNewDDLNumber(); }

    UInt_t GetStartTimeBin()   const { return fBunch.GetStartTimeBin();       } // Provide the index if the first time-bin in current bunch
    UInt_t GetEndTimeBin()     const { return fBunch.GetEndTimeBin();         } // Provide the index if the last time-bin in current bunch
    const UInt_t* GetSignals() const { return fBunch.GetData();               } // Provide access to altro data itself


  private :
    AliAltroRawStreamFast& operator = (const AliAltroRawStreamFast& stream);
    AliAltroRawStreamFast(const AliAltroRawStreamFast& stream);


    Int_t            fDDLNumber;    // index of current DDL number
    Int_t            fPrevDDLNumber;// index of previous DDL number

    AliAltroDecoder  fDecoder;      // decoder for altro payload
    AliAltroData     fData;         // container for altro payload
    AliAltroBunch    fBunch;        // container for altro bunches

    AliRawReader*    fRawReader;    // object for reading the raw data

    ClassDef(AliAltroRawStreamFast, 0)  // base class for fast reading of Altro raw data
};

#endif
