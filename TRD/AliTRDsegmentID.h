#ifndef ALITRDSEGMENTID_H
#define ALITRDSEGMENTID_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDsegmentID.h,v */

////////////////////////////////////////////////
//  Manager class generaol Alice segment 
//  segment is for example one pad row in TPC //
////////////////////////////////////////////////

#include "TObject.h"

class AliTRDsegmentID: public TObject{
public:
  AliTRDsegmentID();
  AliTRDsegmentID(Int_t index){fSegmentID = index;}
  Int_t GetID() {return fSegmentID;}
  void  SetID(Int_t index){fSegmentID = index;} 
  virtual Int_t  GetSize() {return sizeof(*this);} //function which return size of object 
protected:
  Int_t fSegmentID;   //identification number of Segment
  ClassDef(AliTRDsegmentID,1) 
};
   
#endif //ALISEGMENTID_H

