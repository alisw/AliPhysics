#ifndef ALITRDTRACKLETBASE_H
#define ALITRDTRACKLETBASE_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDtrackletBase.h 26327 2008-06-02 15:36:18Z cblume $ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
// TRD tracklet                                                           //
// abstract base class for TRD tracklets                                  //
//                                                                        //
// Authors                                                                //
//  Alex Bercuci (A.Bercuci@gsi.de)                                       //
//  Jochen Klein (jochen.klein@cern.ch)                                   //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#ifndef ROOT_TObject
#include "TObject.h"
#endif

class AliTRDtrackletBase : public TObject {

 public:

    AliTRDtrackletBase() : TObject() {}
    AliTRDtrackletBase(const AliTRDtrackletBase &o) : TObject(o) {}
    virtual ~AliTRDtrackletBase() {}

    virtual Bool_t   CookPID() = 0;

    virtual Int_t    GetDetector() const = 0 ;

    virtual Float_t  GetX() const  = 0;
    virtual Float_t  GetY() const  = 0;
    virtual Float_t  GetZ() const  = 0;
    virtual Float_t  GetdYdX() const = 0;
    virtual Float_t  GetdZdX() const { return 0; }

    virtual Int_t    GetdY() const = 0;     // in units of 140um
    virtual Int_t    GetYbin() const  = 0;  // in units of 160um
    virtual Int_t    GetZbin() const  = 0;  // in pad length units

    virtual Double_t GetPID(Int_t is=-1) const = 0;

    virtual void     LocalToGlobal(Float_t&, Float_t&, Float_t&, Float_t&) {}

    virtual void     Print(Option_t * /*option=""*/) const {}

    virtual UInt_t   GetTrackletWord() const = 0;

    virtual void     SetDetector(Int_t id) = 0;

 protected:

    ClassDef(AliTRDtrackletBase, 1);        // Base class for TRD on- and offline tracklets

};

#endif
