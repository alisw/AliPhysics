#ifndef ALITRDSEGMENTID_H
#define ALITRDSEGMENTID_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDsegmentID.h,v */

////////////////////////////////////////////////
//     Base class for a detector segment      // 
////////////////////////////////////////////////

#include <TObject.h>

class AliTRDsegmentID : public TObject {

 public:

  AliTRDsegmentID();
  AliTRDsegmentID(Int_t index);
  virtual ~AliTRDsegmentID();

          Int_t  GetID() const      { return fSegmentID;    }
  virtual Int_t  GetSize()          { return sizeof(*this); }

          void   SetID(Int_t index) { fSegmentID = index;} 

 protected:

  Int_t fSegmentID;           // Identification number of a segment

  ClassDef(AliTRDsegmentID,1) // Detector segment base class

};
   
#endif 

