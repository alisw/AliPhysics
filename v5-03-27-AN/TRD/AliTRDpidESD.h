#ifndef ALITRDPIDESD_H
#define ALITRDPIDESD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
// Assigns the PID probabilities based on TRD information to the ESDs     //
//                                                                        //
// Authors :                                                              //
//   Prashant Shukla <shukla@pi0.physi.uni-heidelberg.de> (orig. version) //
//   Alex Bercuci (a.bercuci@gsi.de)                                      //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

//#include <Rtypes.h>

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
  static  Bool_t  CheckTrack(AliESDtrack * const t);
          Int_t   MakePID(AliESDEvent * const event);

          void    SetCheckTrackStatus(Bool_t status = kTRUE) { fgCheckTrackStatus = status; };
          void    SetCheckKinkStatus(Bool_t status = kTRUE)  { fgCheckKinkStatus  = status; };
          void    SetMinPlane(Int_t plane)                   { fgMinPlane         = plane;  };

	  Bool_t  GetCheckTrackStatus() const                { return fgCheckTrackStatus;   };      
	  Bool_t  GetCheckKinkStatus() const                 { return fgCheckKinkStatus;    };      
          Int_t   GetMinPlane() const                        { return fgMinPlane;           };

private:

          Bool_t  RecalculateTrackSegmentKine(AliESDtrack * const t
                                            , Int_t plan
                                            , Float_t &mom
                                            , Float_t &length);

  static  Bool_t  fgCheckTrackStatus;           //  Enable check on ESD track status
  static  Bool_t  fgCheckKinkStatus;            //  Enable check on ESD kink track
  static  Int_t   fgMinPlane;                   //  Minimum number of planes

  AliExternalTrackParam *fTrack;                //! Memory holder for Track segment calculations
	
  ClassDef(AliTRDpidESD,2)                      //  TRD PID class

};

#endif


