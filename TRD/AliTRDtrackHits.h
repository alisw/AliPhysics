#ifndef ALITRDTRACKHITS_H
#define ALITRDTRACKHITS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////
//                                            //
//  Manager class for TRD   hits              //
//                                            //
////////////////////////////////////////////////

#include "../TPC/AliTPCTrackHitsV2.h"
class AliTRDhit;

class AliTRDtrackHits : public AliTPCTrackHitsV2 {
public:
  void AddHitTRD(Int_t volumeID, Int_t trackID, Double_t x, 
		    Double_t y, Double_t z,Int_t q, Bool_t inDrift);
  Bool_t First(); //set current hit to first hit 
  Bool_t Next();  //set current hit to next
public:
  ClassDef(AliTRDtrackHits,1) 
};


#endif //ALITRDTRACKHITS_H
