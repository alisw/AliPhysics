#ifndef ALIMUONRAWSTREAMTRIGGER_H
#define ALIMUONRAWSTREAMTRIGGER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*$Id$*/

/// \ingroup base
/// \class AliMUONRawStreamTrigger
/// \brief Class for reading MUON raw digits
///
/// \author Christian Finck
///
///////////////////////////////////////////////////////////////////////////////
///
/// This class provides access to MUON digits in raw data.
///
///////////////////////////////////////////////////////////////////////////////

#include <TObject.h>
#include <TClonesArray.h>
#include "AliMUONPayloadTrigger.h"

class AliRawReader;
class AliMUONDDLTrigger;


class AliMUONRawStreamTrigger: public TObject {
  public :
    AliMUONRawStreamTrigger();
    AliMUONRawStreamTrigger(AliRawReader* rawReader);
    virtual ~AliMUONRawStreamTrigger();

    virtual Bool_t   Next();
    virtual Bool_t   NextDDL();

    Int_t GetMaxDDL() const {return fMaxDDL;}
    Int_t GetMaxReg() const {return fPayload->GetMaxReg();}
    Int_t GetMaxLoc() const {return fPayload->GetMaxLoc();}


    void SetMaxDDL(Int_t ddl);
    void SetMaxReg(Int_t reg);
    void SetMaxLoc(Int_t loc);


    void SetReader(AliRawReader* rawReader) {fRawReader = rawReader;}

    AliMUONDDLTrigger* GetDDLTrigger() const {return fPayload->GetDDLTrigger();}
    Int_t              GetDDL()        const {return fDDL - 1;}

  private :

    AliRawReader*    fRawReader;     ///< object for reading the raw data
    AliMUONPayloadTrigger* fPayload; ///< pointer to payload decoder

    Int_t  fDDL;          ///< number of DDL
    Int_t  fSubEntries;   ///< entries of buspatch structure
    Bool_t fNextDDL;      ///< flag for next DDL to be read
    Int_t  fMaxDDL;       ///< maximum number of DDL in DATE file

    AliMUONRawStreamTrigger(const AliMUONRawStreamTrigger& stream);
    AliMUONRawStreamTrigger& operator = (const AliMUONRawStreamTrigger& stream);

    ClassDef(AliMUONRawStreamTrigger, 3)    // base class for reading MUON trigger rawdata
};

#endif
