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

class AliESD;

class AliTRDpidESD : public TObject {

 public:

  AliTRDpidESD();
  AliTRDpidESD(const AliTRDpidESD &p);
  virtual ~AliTRDpidESD() {}
  AliTRDpidESD &operator=(const AliTRDpidESD &p);

  virtual void    Copy(TObject &p) const;

  static  Int_t   MakePID(AliESD *event);

          void    SetCheckTrackStatus(Bool_t status = kTRUE) { fCheckTrackStatus = status; };
          void    SetCheckKinkStatus(Bool_t status = kTRUE)  { fCheckKinkStatus  = status; };
          void    SetMinPlane(Int_t plane)                   { fMinPlane         = plane;  };

	  Bool_t  GetCheckTrackStatus()                      { return fCheckTrackStatus;   };      
	  Bool_t  GetCheckKinkStatus()                       { return fCheckKinkStatus;    };      
          Int_t   GetMinPlane()                              { return fMinPlane;           };

 private:

  static  Bool_t  fCheckTrackStatus;    // Enable check on ESD track status
  static  Bool_t  fCheckKinkStatus;     // Enable check on ESD kink track
  static  Int_t   fMinPlane;            // Minimum number of planes

  ClassDef(AliTRDpidESD,2)              // TRD PID class

};

#endif


