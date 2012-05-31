#ifndef ALISEGMENTID_H
#define ALISEGMENTID_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Manager class generaol Alice segment 
//  segment is for example one pad row in TPC //
////////////////////////////////////////////////

#include "TObject.h"

class AliSegmentID: public TObject{
public:
  AliSegmentID() : fSegmentID(0) {}
  AliSegmentID(Int_t index) : fSegmentID(index) {}
  Int_t GetID() {return fSegmentID;}
  void  SetID(Int_t index){fSegmentID = index;} 
protected:
  Int_t fSegmentID;   //identification number of Segment
  ClassDef(AliSegmentID,1) 
};
   
#endif //ALISEGMENTID_H

