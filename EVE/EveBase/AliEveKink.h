// $Id$
// Main authors: Paraskevi Ganoti: 2009

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/
#ifndef AliEveKink_H
#define AliEveKink_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               
***************************************************************************/

#include <TEveVSDStructs.h>
#include <TEveElement.h>
#include <TEveTrack.h>
#include <TEveTrackPropagator.h>

#include <TPolyMarker3D.h>
#include <TPolyLine3D.h>

#include <TPDGCode.h>

class TH1F;
class TH2F;

class AliEveKinkList;

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

  void SetRnrStyleMother(TEveTrackPropagator* rsMoth)  {fRnrStyleMoth=rsMoth;}
  void SetRnrStyleDaughter(TEveTrackPropagator* rsDaugh)  {fRnrStyleDaugh=rsDaugh;}   
  Double_t GetKinkAngle(Int_t i) const { return fKinkAngle[i]; }
  void SetKinkAngle(Int_t i, Double_t anglekink) {fKinkAngle[i] = anglekink; }
  
  Float_t GetKinkRadius() const { return fRecKinkPosition.Perp(); }
  
  Float_t GetKinkPMother() const { return fMotherMomentum.Mag(); }
  Float_t GetKinkPMotherPerp() const { return fMotherMomentum.Perp(); } 
  Float_t GetKinkPDaughter() const { return fDaughterMomentum.Mag(); }
  
  Float_t GetInvMass(Int_t dPdgCode) const;
  Float_t GetQt() const;
  
  void SetMaxProbPdgPid(Int_t rPdg, Float_t rPid);
  Int_t   GetDaugMaxProbPdg() const { return fDaugMaxProbPdg; }
  Float_t GetDaugMaxProbPid() const { return fDaugMaxProbPid; }

  Int_t GetESDKinkIndex() const { return fESDKinkIndex; }
  void  SetESDKinkIndex(Int_t ind) { fESDKinkIndex = ind;}

  virtual const Text_t* GetName() const    { return Form("ESDkink_%i",fESDKinkIndex); }
  virtual const Text_t* GetTitle() const   { return Form("ESDkink_%i",fESDKinkIndex); }

  TEveTrack* GetMotherTrack() { return fMotherTrack; }
  TEveTrack* GetDaughterTrack() { return fDaughterTrack; }

protected:
  TEveVector       fRecKinkPosition;
  TEveVector       fMotherMomentum;
  TEveVector       fDaughterMomentum;

  TEveTrack        *fMotherTrack;
  TEveTrack        *fDaughterTrack;

  TEveTrackPropagator  *fRnrStyleMoth;
  TEveTrackPropagator  *fRnrStyleDaugh;
  
  Int_t             fESDKinkIndex;    // Index in ESD Kink array.
  Double_t          fKinkAngle[3]; //
  
  Int_t             fDaugMaxProbPdg; // Maximum PDG probability for the daughter
  Float_t           fDaugMaxProbPid; // Maximum PID probability for the daughter 

private:
  AliEveKink(const AliEveKink&);            // Not implemented
  AliEveKink& operator=(const AliEveKink&); // Not implemented

  ClassDef(AliEveKink, 0); // Visual representation of a AliEveKink.
};


/******************************************************************************/
// AliEveKinkList
/******************************************************************************/

class AliEveKinkList : public TEveElementList
{
  friend class AliEveKinkListEditor;

public:
  AliEveKinkList();
  AliEveKinkList(TEveTrackPropagator* rsMoth, TEveTrackPropagator* rsDaugh);
  AliEveKinkList(const Text_t* name, TEveTrackPropagator* rsMoth=0, TEveTrackPropagator* rsDaugh=0);
  virtual ~AliEveKinkList() {}

  virtual const Text_t* GetTitle() const { return fTitle; }
  virtual void SetTitle(const Text_t* t) { fTitle = t; }
  virtual void SetTracksColor(Color_t cMoth, Color_t cDaug) {
  fMothColor = cMoth; fDaugColor = cDaug;}

  virtual Bool_t CanEditMainColor() const { return kTRUE; }

  void  SetRnrStyleMoth(TEveTrackPropagator* rstMoth) { fRnrStyleMoth = rstMoth; }
  TEveTrackPropagator* GetPropagatorMoth()        { return fRnrStyleMoth; }
   
  void  SetRnrStyleDaugh(TEveTrackPropagator* rstDaugh) { fRnrStyleDaugh = rstDaugh; }
  TEveTrackPropagator* GetPropagatorDaugh()        { return fRnrStyleDaugh; }

  Bool_t GetRnrKinkvtx()     const { return fRnrKinkvtx; }
  Bool_t GetRnrKinkDaughter() const { return fRnrKinkDaughter; }   //not yet be sure about this!!!

  void   MakeKinks();

  void   FilterByRadius(Float_t minR, Float_t maxR);
  void   FilterByKinkAngle(Float_t minKinkAngle, Float_t maxKinkAngle);
  void   FilterByPt(Float_t minPt, Float_t maxPt);
  void   FilterByInvariantMass(Float_t minPt, Float_t maxPt, Int_t dPdgCode);
  
  void   FilterByCheckedPidMinProb(Int_t rFlag, Int_t rPid, Float_t rProb);
  void   SetDaugCheckedPid(Int_t rDaugCheckedPid) {fDaugCheckedPid = rDaugCheckedPid;}
  Int_t  GetDaugCheckedPid() {return fDaugCheckedPid;}

  void   SetDaugCheckedProb(Float_t rDaugCheckedProb) {fDaugCheckedProb = rDaugCheckedProb;}
  Float_t  GetDaugCheckedProb() {return fDaugCheckedProb;}


protected:
  TString              fTitle;

  TEveTrackPropagator *fRnrStyleMoth;
  TEveTrackPropagator *fRnrStyleDaugh; 

  Bool_t               fRnrKinkDaughter;
  Bool_t               fRnrKinkvtx;

  Color_t              fMothColor;
  Color_t              fDaugColor;

  Float_t              fMinRCut;
  Float_t              fMaxRCut;

  Float_t              fMinKinkAngle;
  Float_t              fMaxKinkAngle;

  Float_t              fMinPt;
  Float_t              fMaxPt;

  Float_t              fMinInvariantMass;
  Float_t              fMaxInvariantMass;
  
  Int_t                fDaugCheckedPid;
  Float_t              fDaugCheckedProb;  

private:
  void Init();

  AliEveKinkList(const AliEveKinkList&);            // Not implemented
  AliEveKinkList& operator=(const AliEveKinkList&); // Not implemented

  ClassDef(AliEveKinkList, 0); // A list of AliEveKink objecs.
};


#endif
