#ifndef ALITRDMODULE_H
#define ALITRDMODULE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDmodule.h 19198 2007-06-19 13:50:47Z cblume $ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD module class                                                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TObject.h>

class AliTRDgtuTrack;
class AliTRDltuTracklet;

class AliTRDmodule : public TObject {

 public:

  enum { kNplan = 6, kNmaxZchan = 100, kNsubZchan = 16, kNmaxTrk = 12 };

  AliTRDmodule();
  AliTRDmodule(const AliTRDmodule &m);
  virtual         ~AliTRDmodule();
  AliTRDmodule    &operator=(const AliTRDmodule &m);

  virtual void     Copy(TObject &m) const;

          Int_t    GetNtracklets() const;
          Int_t    GetNtracks() const;

          void     Reset();
          void     AddTracklet(Int_t det, Int_t row, Float_t rowz, Float_t slope, Float_t offset 
  		             , Float_t time, Int_t ncl, Int_t label, Float_t q);
          void     AddTrack();
 
          void     ResetTracklets();
          void     ResetTracks();
          void     SortTracklets();
          void     SortTracks();
          void     RemoveMultipleTracklets();
          void     RemoveMultipleTracks();
          void     RemoveTracklet(Int_t pos);
          void     RemoveTrack(Int_t pos);
          void     SortZ(Int_t cha);
          void     InitZLUT();
          void     FindTracks();
          void     FindTracksCombi(Int_t zchan);

          TObjArray         *Tracklets(); 
          TObjArray         *Tracks();
          AliTRDltuTracklet *GetTracklet(Int_t pos) const;
          AliTRDgtuTrack    *GetTrack(Int_t pos) const;

 protected:

          Float_t            fXprojPlane;                              //! X (time) coordinate of the projection plane
          Float_t            fField;                                   //! Magnetic field
          TObjArray         *fTracklets;                               //! Array of LTU tracklets
          TObjArray         *fTracks;                                  //! Array of GTU tracks

          Int_t              fZnchan[kNplan][kNsubZchan];              //! Number of LTU tracklets in each subchannel
          Int_t              fZtrkid[kNplan][kNmaxZchan][kNsubZchan];  //! List of LTU tracklet id's for each subchannel

          Float_t            fDeltaY;                                  //  Y (offset) matching window in the GTU
          Float_t            fDeltaS;                                  //  Slope matching window in the GTU

          AliTRDltuTracklet *fLTUtrk;                                  //! Current LTU tracklet
          AliTRDgtuTrack    *fGTUtrk;                                  //! Current GTU track

  ClassDef(AliTRDmodule,2)                                             //  TRD module class

};

#endif
