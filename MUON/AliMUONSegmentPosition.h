#ifndef ALIMUONSEGMENTPOSITION_H
#define ALIMUONSEGMENTPOSITION_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//===================================================================
//  Segment element position in local coordinates of the detection element   
//        Gines MARTINEZ, SUBATECH July 04                
//  This class is one of the basic component of 
//  AliMUONSegmentationDetectionElement and contains al the 
//  info about a segment (pad or strip):
//          Id-indetectionelement,  x_local, y_local 
//  Detailed information in Alice Technical Note xxxxxxxx (2004)
//====================================================================

#include <TNamed.h>

class AliMUONSegmentPosition : public TNamed
{
 public:
    AliMUONSegmentPosition();
    AliMUONSegmentPosition(const Int_t channelId, const Float_t x, const  Float_t y, const Int_t cathode);
    virtual ~AliMUONSegmentPosition();
      
    Int_t   Compare(const TObject *obj) const;
    Int_t   GetChannelId()const {return fChannelId;}
    Float_t GetXlocal()   const {return fX;}
    Float_t GetYlocal()   const {return fY;}
    Int_t   GetCathode()  const {return fCathode;}

    void Print() const;
 private:
    Int_t fChannelId;   // Id of the channel within the detection element
    Float_t fX;
    Float_t fY;
    Int_t fCathode;
     
    ClassDef(AliMUONSegmentPosition,1) // Loal positions of segments
	
};
#endif






