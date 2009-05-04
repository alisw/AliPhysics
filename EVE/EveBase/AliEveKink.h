// $Id$
// Main authors: Paraskevi Ganoti: 2009

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/
#ifndef AliEveKink_H
#define AliEveKink_H

#include <TEvePointSet.h>

class TEveTrack;
class TEveTrackPropagator;

class TH1F;
class TH2F;

class AliEveKinkList;

//------------------------------------------------------------------------------
// AliEveKink
//
// Graphical representation of a kink.
//------------------------------------------------------------------------------

class AliEveKink : public TEvePointSet
{
  friend class AliEveKinkList;
  friend class AliEveKinkEditor;

public:
  AliEveKink();
  AliEveKink(TEveRecTrack* tMoth, TEveRecTrack* tDaug, TEveRecKink* kink, TEveTrackPropagator* rsMoth, TEveTrackPropagator* rsDaugh);
  virtual ~AliEveKink();

  void MakeKink();

  virtual void SetMainColor(Color_t col)
  {
    TEvePointSet::SetMainColor(col);
  }

  void SetRnrStyleMother(TEveTrackPropagator* rsMoth)    { fRnrStyleMoth  = rsMoth;  }
  void SetRnrStyleDaughter(TEveTrackPropagator* rsDaugh) { fRnrStyleDaugh = rsDaugh; }
  Double_t GetKinkAngle(Int_t i) const { return fKinkAngle[i]; }
  void SetKinkAngle(Int_t i, Double_t anglekink) { fKinkAngle[i] = anglekink; }

  Float_t GetKinkRadius() const { return fRecKinkPosition.Perp(); }

  Float_t GetKinkPMother()     const { return fMotherMomentum.Mag(); }
  Float_t GetKinkPMotherPerp() const { return fMotherMomentum.Perp(); }
  Float_t GetKinkPDaughter()   const { return fDaughterMomentum.Mag(); }

  Float_t GetInvMass(Int_t dPdgCode) const;
  Float_t GetQt() const;

  void SetMaxProbPdgPid(Int_t rPdg, Float_t rPid);
  Int_t   GetDaugMaxProbPdg() const { return fDaugMaxProbPdg; }
  Float_t GetDaugMaxProbPid() const { return fDaugMaxProbPid; }

  Int_t GetESDKinkIndex() const { return fESDKinkIndex; }
  void  SetESDKinkIndex(Int_t ind) { fESDKinkIndex = ind;}

  virtual const Text_t* GetName()  const { return Form("ESDkink_%i",fESDKinkIndex); }
  virtual const Text_t* GetTitle() const { return Form("ESDkink_%i",fESDKinkIndex); }

  TEveTrack* GetMotherTrack()   const { return fMotherTrack; }
  TEveTrack* GetDaughterTrack() const { return fDaughterTrack; }

protected:
  TEveVector       fRecKinkPosition;  // Reconstructed position of kink.
  TEveVector       fMotherMomentum;   // Momentum of mother track.
  TEveVector       fDaughterMomentum; // Momentum of daugter track.

  TEveTrack        *fMotherTrack;     // Graphical representation of mother track.
  TEveTrack        *fDaughterTrack;   // Graphical representation of daughter track.

  TEveTrackPropagator  *fRnrStyleMoth;  // Track-propagator for mother track.
  TEveTrackPropagator  *fRnrStyleDaugh; // Track-propagator for daughter track.

  Int_t             fESDKinkIndex;    // Index in ESD Kink array.
  Double_t          fKinkAngle[3];    // TODO

  Int_t             fDaugMaxProbPdg;  // Maximum PDG probability for the daughter
  Float_t           fDaugMaxProbPid;  // Maximum PID probability for the daughter

private:
  AliEveKink(const AliEveKink&);            // Not implemented
  AliEveKink& operator=(const AliEveKink&); // Not implemented

  ClassDef(AliEveKink, 0); // Visual representation of a AliEveKink.
};


//------------------------------------------------------------------------------
// AliEveKinkList
//
// Container for AliEveKink objects
// Provides managmenet methods for setting cuts and common visualization
// parameters.
//------------------------------------------------------------------------------

class AliEveKinkList : public TEveElementList
{
  friend class AliEveKinkListEditor;

public:
  AliEveKinkList();
  AliEveKinkList(TEveTrackPropagator* rsMoth, TEveTrackPropagator* rsDaugh);
  AliEveKinkList(const Text_t* name, TEveTrackPropagator* rsMoth=0, TEveTrackPropagator* rsDaugh=0);
  virtual ~AliEveKinkList() {}

  virtual void SetTracksColor(Color_t cMoth, Color_t cDaug)
  { fMothColor = cMoth; fDaugColor = cDaug; }

  virtual Bool_t CanEditMainColor() const { return kTRUE; }

  void  SetRnrStyleMoth(TEveTrackPropagator* rstMoth) { fRnrStyleMoth = rstMoth; }
  TEveTrackPropagator* GetPropagatorMoth() const      { return fRnrStyleMoth; }

  void  SetRnrStyleDaugh(TEveTrackPropagator* rstDaugh) { fRnrStyleDaugh = rstDaugh; }
  TEveTrackPropagator* GetPropagatorDaugh() const       { return fRnrStyleDaugh; }

  Bool_t GetRnrKinkvtx()      const { return fRnrKinkvtx; }
  Bool_t GetRnrKinkDaughter() const { return fRnrKinkDaughter; }   //not yet be sure about this!!!

  void   MakeKinks();

  void   FilterByRadius(Float_t minR, Float_t maxR);
  void   FilterByKinkAngle(Float_t minKinkAngle, Float_t maxKinkAngle);
  void   FilterByPt(Float_t minPt, Float_t maxPt);
  void   FilterByInvariantMass(Float_t minPt, Float_t maxPt, Int_t dPdgCode);

  void   FilterByCheckedPidMinProb(Int_t rFlag, Int_t rPid, Float_t rProb);
  void   SetDaugCheckedPid(Int_t dcpid) { fDaugCheckedPid = dcpid; }
  Int_t  GetDaugCheckedPid()      const { return fDaugCheckedPid; }

  void    SetDaugCheckedProb(Float_t dcprob) { fDaugCheckedProb = dcprob; }
  Float_t GetDaugCheckedProb()         const { return fDaugCheckedProb; }


protected:
  TEveTrackPropagator *fRnrStyleMoth;     // Default track-propagator for mother tracks.
  TEveTrackPropagator *fRnrStyleDaugh;    // Default track-propagator for daughter tracks.

  Bool_t               fRnrKinkDaughter;  // Flag - show daughter tracks.
  Bool_t               fRnrKinkvtx;       // Flag - show kink vertex.

  Color_t              fMothColor;        // Color of mother tracks.
  Color_t              fDaugColor;        // Color of daughter tracks.

  Float_t              fMinRCut;          // Cut - minimum kink radius.
  Float_t              fMaxRCut;          // Cut - maximum kink radius.

  Float_t              fMinKinkAngle;     // Cut - minimum kink angle.
  Float_t              fMaxKinkAngle;     // Cut - maximum kink angle.

  Float_t              fMinPt;            // Cut - minimum pT of mother track.
  Float_t              fMaxPt;            // Cut - maximum pT of mother track.

  Float_t              fMinInvariantMass; // Cut - minimum invariant mass.
  Float_t              fMaxInvariantMass; // Cut - maximum invariant mass.

  Int_t                fDaugCheckedPid;   // Cut - PID of daughter track.
  Float_t              fDaugCheckedProb;  // Cut - min PID probability of the daughter track.

private:
  void Init();

  AliEveKinkList(const AliEveKinkList&);            // Not implemented
  AliEveKinkList& operator=(const AliEveKinkList&); // Not implemented

  ClassDef(AliEveKinkList, 0); // A list of AliEveKink objecs.
};


#endif
