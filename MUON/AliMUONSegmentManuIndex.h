#ifndef ALIMUONSEGMENTMANUINDEX_H
#define ALIMUONSEGMENTMANUINDEX_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//===================================================================
//  Segment element indexing in a detection element for electronics   
//        Gines MARTINEZ, SUBATECH July 04                
//  This class is the basic component of 
//  AliMUONSegmentationDetectionElement and contains al the 
//  info about a segment (pad or strip):
//          Id-indetectionelement,  #manu, #manuchannel 
//  Detailed information in Alice Technical Note xxxxxxxx (2004)
//====================================================================

#include <TNamed.h>

class TString;

class AliMUONSegmentManuIndex : public TNamed {
 public:
  AliMUONSegmentManuIndex();
  AliMUONSegmentManuIndex(Int_t channelId, Int_t manuId, Int_t busPatchId, Int_t manuChannelId);

  virtual ~AliMUONSegmentManuIndex();

  Int_t Compare(const TObject *obj) const;

  Int_t GetChannelId()     const{return fChannelId;}
  Int_t GetManuId()        const{return fManuId;}
  Int_t GetBusPatchId()    const{return fBusPatchId;}
  Int_t GetManuChannelId() const{return fManuChannelId;}

  static TString Name(Int_t manuId, Int_t manuchannel);
  
  void   Print() const;

 private:
  Int_t fChannelId; // Id of the channel within the detection element
  Int_t fManuId; // Manu id in the detection element
  Int_t fBusPatchId; // BusPatchId in the detection element up to 4 for slats
  Int_t fManuChannelId; // ChannelId in the manu card 1-64
  
  ClassDef(AliMUONSegmentManuIndex,1) // Segmenation for MUON detection elements
    
};

#endif






