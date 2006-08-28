#ifndef ALITRDTIMEBIN_H
#define ALITRDTIMEBIN_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$*/

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Hit compression class                                                 //
//  Provides tools to address clusters which lie within one time bin      //
//  Adapted from AliTPCTimeBin by Marian                                  //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include <TObject.h>

class AliTRDcluster;

class AliTRDtimeBin : public TObject {

 public: 

  AliTRDtimeBin();
  virtual ~AliTRDtimeBin() { };

                   operator Int_t() const                         { return fN;        }
  AliTRDcluster   *operator[](Int_t i);

          void     InsertCluster(AliTRDcluster *c, UInt_t index);
          UInt_t   GetIndex(Int_t i) const                        { return fIndex[i]; } 
          Int_t    Find(Double_t y) const; 

 protected:

   enum { kMaxClusterPerTimeBin = 3500 };
 
          UInt_t   fN;                                 // ????
   AliTRDcluster  *fClusters[kMaxClusterPerTimeBin];   // ????
          UInt_t   fIndex[kMaxClusterPerTimeBin];      // ????

  ClassDef(AliTRDtimeBin,1)                            // Provides tools to address clusters within one time bin

}; 

#endif 

