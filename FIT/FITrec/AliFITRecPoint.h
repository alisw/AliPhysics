#ifndef ALIFITRECPOINT_H
#define ALIFITRECPOINT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TObject.h>


//___________________________________________
class AliFITRecPoint: public TObject  {
////////////////////////////////////////////////////////////////////////
 public:
    AliFITRecPoint();
    AliFITRecPoint(const AliFITRecPoint &o);
    AliFITRecPoint& operator= (const AliFITRecPoint &) { return *this;}
    virtual ~AliFITRecPoint() {}

    
     
    void    SetTime (Int_t ipmt, Float_t time) { fTime[ipmt] = time;}
    Float_t GetTime (Int_t ipmt)const { return fTime[ipmt];}
    void    SetAmp (Int_t ipmt, Float_t adc) { fADCQTC[ipmt] = adc;}
    Float_t GetAmp (Int_t ipmt) const{ return fADCQTC[ipmt];}
    
 
 
 
  private: 
 
    Float_t fTime[300];    // array's TDC
    Float_t fADCQTC[300];    // array's amplitude

 
    ClassDef(AliFITRecPoint,1)  // RecPoints (Header) object for set:T0
};

#endif



