// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel & Bogdan Vulpescu: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/
#ifndef AliEveMUONTrack_H
#define AliEveMUONTrack_H

#include <TEveTrack.h>

class AliMUONTrack;
class AliMUONTriggerTrack;
class AliMagF;
class AliESDMuonTrack;
class AliTrackReference;

class TParticle;

class TEveTrackPropagator;


class AliEveMUONTrack: public TEveTrack
{
 public:

  AliEveMUONTrack(TEveRecTrack* t, TEveTrackPropagator* rs);
  virtual ~AliEveMUONTrack();

  virtual void MakeTrack(Bool_t /*recurse*/=kFALSE) {}

  void  MakeMUONTrack(AliMUONTrack *mtrack);
  void  MakeMUONTriggerTrack(AliMUONTriggerTrack *mtrack);
  void  MakeESDTrack(AliESDMuonTrack *mtrack);
  void  MakeMCTrack(TParticle *part);
  void  MakeRefTrack(AliMUONTrack *mtrack);
  void  Propagate(Float_t *xr, Float_t *yr, Float_t *zr, Int_t i1, Int_t i2);
  void  OneStepRungekutta(Double_t charge, Double_t step,
			  Double_t* vect, Double_t* vout);
  Int_t ColorIndex(Float_t val);

  Bool_t IsMUONTrack()        const { return fIsMUONTrack; }
  Bool_t IsMUONTriggerTrack() const { return fIsMUONTrack; }
  Bool_t IsESDTrack()         const { return fIsESDTrack;  }
  Bool_t IsMCTrack()          const { return fIsMCTrack;   }
  Bool_t IsRefTrack()         const { return fIsRefTrack;  }

  void PrintMCTrackInfo();
  void PrintMUONTrackInfo();
  void PrintMUONTriggerTrackInfo();
  void PrintESDTrackInfo();

  void  MUONTrackInfo();          // *MENU*
  void  MUONTriggerInfo();        // *MENU*

 private:

  AliEveMUONTrack(const AliEveMUONTrack&);            // Not implemented
  AliEveMUONTrack& operator=(const AliEveMUONTrack&); // Not implemented

  AliMUONTrack *fTrack;              // pointer to the MUON track
  TParticle    *fPart;               // pointer to the MC particle
  Int_t         fCount;              // track points counter
  Bool_t        fIsMUONTrack;        // track from MUON.Tracks.root
  Bool_t        fIsMUONTriggerTrack; // trigger track from MUON.Tracks.root
  Bool_t        fIsESDTrack;         // track from AliESDs.root
  Bool_t        fIsMCTrack;          // track from Kinematics.root
  Bool_t        fIsRefTrack;         // track from TrackRefs.root

  ClassDef(AliEveMUONTrack, 0);    // Produce TEveUtil:TEveTrack from AliMUONTrack
};

class AliEveMUONTrackList : public TEveTrackList
{
public:
  AliEveMUONTrackList(TEveTrackPropagator* rs=0) : TEveTrackList(rs) {}
  AliEveMUONTrackList(const Text_t* name, TEveTrackPropagator* rs=0) : TEveTrackList(name, rs) {}
  virtual ~AliEveMUONTrackList() {}

  void HackMomentumLimits(Bool_t recurse=kTRUE);

  ClassDef(AliEveMUONTrackList, 0);    // Temporary workaround for deficiency in TEveTrackList
};

#endif

