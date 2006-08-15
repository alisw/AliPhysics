#ifndef ALITRDTRACKHITS_H
#define ALITRDTRACKHITS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Manager class for TRD hits                                            //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "../TPC/AliTPCTrackHitsV2.h"

class AliTRDhit;

class AliTRDtrackHits : public AliTPCTrackHitsV2 {

 public:

  AliTRDtrackHits()          { };
  virtual ~AliTRDtrackHits() { };

          void     AddHitTRD(Int_t volumeID, Int_t trackID, Double_t x
	                   , Double_t y, Double_t z,Int_t q, Bool_t inDrift);
          Bool_t   First();
          Bool_t   Next();

 public:

  ClassDef(AliTRDtrackHits,1) // Manager class for TRD hits

};

#endif
