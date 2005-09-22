#ifndef ALIMUONDIGIT_H
#define ALIMUONDIGIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Revision of includes 07/05/2004

#include <TObject.h>

static const Int_t kMAXTRACKS=10;

class AliMUONDigit : public TObject {

 public:
    AliMUONDigit();
    AliMUONDigit(const AliMUONDigit& rhs);
    AliMUONDigit(Int_t *digits);
    AliMUONDigit(Int_t *tracks, Int_t *charges, Int_t *digits);
    virtual ~AliMUONDigit();

    AliMUONDigit& operator=(const AliMUONDigit& rhs);
    
    virtual Bool_t IsSortable() const {return kTRUE;}
    virtual int Compare(const TObject *obj) const;

    virtual Int_t    DetElemId()const     {return fDetElemId;}
    virtual Int_t    PadX() const         {return fPadX;}
    virtual Int_t    PadY() const         {return fPadY;}
    virtual Int_t    Signal() const       {return fSignal;}
    virtual Int_t    Physics() const      {return fPhysics;}
    virtual Int_t    Hit() const          {return fHit;}    
    virtual Int_t    Cathode() const      {return fCathode;}
    virtual Int_t    Track(Int_t i) const {return fTracks[i];}
    virtual Int_t    TrackCharge(Int_t i) const {return fTcharges[i];} 

    virtual void     SetDetElemId(Int_t id)    {fDetElemId = id;}
    virtual void     SetPadX(Int_t pad)        {fPadX = pad;}
    virtual void     SetPadY(Int_t pad)        {fPadY = pad;}
    virtual void     SetSignal(Int_t q)        {fSignal = q;}
    virtual void     AddSignal(Int_t q)        {fSignal += q;}
    virtual void     AddPhysicsSignal(Int_t q) {fPhysics += q;}
    virtual void     SetHit(Int_t n)           {fHit = n;}    
    virtual void     SetCathode(Int_t c)       {fCathode = c;}
  
	    
 private:
    Int_t     fPadX;          // Pad number along x
    Int_t     fPadY;          // Pad number along y
    Int_t     fCathode;       // Cathode number
    
    Int_t     fSignal;        // Signal amplitude
    Int_t     fTcharges[kMAXTRACKS];  // charge per track making this digit (up to 10)
    Int_t     fTracks[kMAXTRACKS];    // primary tracks making this digit (up to 10)
    Int_t     fPhysics;       // physics contribution to signal 
    Int_t     fHit;           // hit number - temporary solution
    Int_t     fDetElemId;     // Detection element ID

    ClassDef(AliMUONDigit,2)  //Digits for MUON
};
#endif
