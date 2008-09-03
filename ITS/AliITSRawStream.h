#ifndef ALIITSRAWSTREAM_H
#define ALIITSRAWSTREAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
///
/// This is a base class for providing access to ITS digits in raw data.
///
///////////////////////////////////////////////////////////////////////////////

#include <TObject.h>
#include "AliITSDDLModuleMapSDD.h"
#include "AliLog.h"

class AliRawReader;


class AliITSRawStream: public TObject {
  public :
    AliITSRawStream(AliRawReader* rawReader);
    AliITSRawStream(const AliITSRawStream& stream);
    AliITSRawStream& operator = (const AliITSRawStream& stream);
    virtual ~AliITSRawStream() {};

    virtual Bool_t   Next() = 0;

    Int_t            GetModuleID() const {return fModuleID;};
    Int_t            GetPrevModuleID() const {return fPrevModuleID;};
    Bool_t           IsNewModule() const {return fModuleID != fPrevModuleID;};
    Int_t            GetCoord1() const {return fCoord1;};
    Int_t            GetCoord2() const {return fCoord2;};
    Int_t            GetSignal() const {return fSignal;};
    virtual Bool_t   IsCompletedModule() const {return fCompletedModule;}; // to be implemented in derived class
    virtual void     SetDDLModuleMap(AliITSDDLModuleMapSDD* /*ddlsdd*/){
      AliError("This method must be implemented in a derived class");
    };
    virtual void     SetZeroSuppLowThreshold(Int_t /*iMod*/, Int_t /*iSid*/, Int_t /*th*/) {
      AliError("This method must be implemented in a derived class");
    };
    virtual Int_t     GetCarlosId() const {
      AliError("This method must be implemented in a derived class");
      return -1;
    };
    virtual Int_t     GetChannel() const {
      AliError("This method must be implemented in a derived class");
      return -1;
    };


  protected :
    AliRawReader*    fRawReader;    // object for reading the raw data

    Int_t            fModuleID;     // index of current module
    Int_t            fPrevModuleID; // index of previous module
    Int_t            fCoord1;       // current 1st coordinate
                                    //  SPD: column cell number (z)
                                    //  SDD: anode cell number (z)
                                    //  SSD: N/P, flag for side
    Int_t            fCoord2;       // current 2nd coordinate
                                    //  SPD: row cell number (y)
                                    //  SDD: time bin number (y)
                                    //  SSD: strip number
    Int_t            fSignal;       // signal in ADC counts
    Bool_t           fCompletedModule; // set to kTRUE when all data from a module (SDD) are read

    ClassDef(AliITSRawStream, 1) // base class for reading ITS raw digits
};

#endif
