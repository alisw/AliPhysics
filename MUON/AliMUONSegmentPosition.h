#ifndef ALIMUONSEGMENTPOSITION_H
#define ALIMUONSEGMENTPOSITION_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/// \ingroup base
/// \class AliMUONSegmentPosition
/// \brief Local positions of segments

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
#include <TString.h>

class AliMUONSegmentPosition : public TNamed
{
 public:
    AliMUONSegmentPosition();
    AliMUONSegmentPosition(Int_t channelId, Float_t x, Float_t y, Int_t cathode);
    virtual ~AliMUONSegmentPosition();
      
    Int_t   Compare(const TObject *obj) const;
    Float_t Distance(Float_t x, Float_t y);
    Int_t   GetChannelId()const {return fChannelId;}
    Float_t GetXlocal()   const {return fX;}
    Float_t GetYlocal()   const {return fY;}
    Int_t   GetCathode()  const {return fCathode;}

    static  Float_t GetUnit()            {return fUnit;} 
    static  TString Name(Float_t x, Float_t y, Int_t cathode) ;

    void    Print(const char* opt="") const;

 private:
    Int_t   fChannelId;   // Id of the channel within the detection element
    Float_t fX;           // Position X of the center of the segment (pad, strip, etc...)
    Float_t fY;           // Position Y of the center of the segment (pad, strip, etc...)
    Int_t   fCathode;     // Cathode Side Bending 1  or non bending 0 
    Float_t fPadSizeX;
    Float_t fPadSizeY;

    static Float_t fUnit;  // Unit for generation of the name 3mm has been choses     

    ClassDef(AliMUONSegmentPosition,1) // Local positions of segments
	
};
#endif






