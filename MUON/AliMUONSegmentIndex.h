#ifndef ALIMUONSEGMENTINDEX_H
#define ALIMUONSEGMENTINDEX_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//===================================================================
//  Segment element indexing in a detection element    
//        Gines MARTINEZ, SUBATECH July 04                
//  This class is the basic component of 
//  AliMUONSegmentationDetectionElement and contains al the 
//  info about a segment (pad or strip):
//          Id-indetectionelement,  ix ,iy 
//  Detailed information in Alice Technical Note xxxxxxxx (2004)
//====================================================================

#include <TNamed.h>

class AliMUONSegmentIndex : public TNamed {
 public:
  AliMUONSegmentIndex();
  AliMUONSegmentIndex(const Int_t channelId, const Int_t padX, const Int_t padY, const Int_t cathode);
  virtual ~AliMUONSegmentIndex();
  

  Int_t Compare(const TObject *obj) const;
  Int_t GetChannelId() const {return fChannelId;}
  Int_t GetPadX()      const {return fPadX;} 
  Int_t GetPadY()      const {return fPadX;} 
  Int_t GetCathode()   const {return fCathode;} 
  
  void Print() const;

 private:
  Int_t fChannelId; // Id of the channel within the detection element
  Int_t fPadX;
  Int_t fPadY;
  Int_t fCathode;
  
  ClassDef(AliMUONSegmentIndex,1) // Segmenation for MUON detection elements	
};
#endif






