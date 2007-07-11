#ifndef ALITRDPIDESD_H
#define ALITRDPIDESD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Assigns the PID probabilities based on TRD information to the ESDs    //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include <Rtypes.h>

#include <TObject.h>

class AliESDEvent;
class AliESDtrack;
class AliExternalTrackParam;
class AliTRDpidESD : public TObject {

 public:

  AliTRDpidESD();
  AliTRDpidESD(const AliTRDpidESD &p);
  virtual ~AliTRDpidESD();
  AliTRDpidESD &operator=(const AliTRDpidESD &p);

  virtual void    Copy(TObject &p) const;
  static  Bool_t  CheckTrack(AliESDtrack *t);
          Int_t   MakePID(AliESDEvent *event);

          void    SetCheckTrackStatus(Bool_t status = kTRUE) { fCheckTrackStatus = status; };
          void    SetCheckKinkStatus(Bool_t status = kTRUE)  { fCheckKinkStatus  = status; };
          void    SetMinPlane(Int_t plane)                   { fMinPlane         = plane;  };

	  Bool_t  GetCheckTrackStatus()                      { return fCheckTrackStatus;   };      
	  Bool_t  GetCheckKinkStatus()                       { return fCheckKinkStatus;    };      
          Int_t   GetMinPlane()                              { return fMinPlane;           };

private:
  Bool_t  RecalculateTrackSegmentKine(AliESDtrack *t, Int_t plan, Float_t &mom, Float_t &length);

private:

  static  Bool_t  fCheckTrackStatus;    // Enable check on ESD track status
  static  Bool_t  fCheckKinkStatus;     // Enable check on ESD kink track
  static  Int_t   fMinPlane;            // Minimum number of planes

	AliExternalTrackParam *fTrack;				//! Memory holder for Track segment calculations
	
  ClassDef(AliTRDpidESD,2)              // TRD PID class

};

#endif


