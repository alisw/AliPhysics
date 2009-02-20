// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/
#ifndef ALIEVECASCADE_H
#define ALIEVECASCADE_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/***********************************************************************
* This code defines the reconstructed cascade visualized with EVE
*
* Boris Hippolyte, IPHC (hippolyt@in2p3.fr)
************************************************************************/

#include "AliEveTrack.h"
#include <TEveVSDStructs.h>

class AliEveCascadeList;

class TH1F;
class TH2F;

class AliEveCascade : public TEvePointSet
{
  friend class AliEveCascadeList;

public:
  AliEveCascade();
  AliEveCascade(TEveRecTrack* tBac, TEveRecV0* v0, TEveRecCascade* cascade, TEveTrackPropagator* rs);
  virtual ~AliEveCascade();

  void MakeCascade();

  virtual void  SetMainColor(Color_t col)
  {
    TEvePointSet::SetMainColor(col);
    fPointingCurve->SetLineColor(fMarkerColor);
    fV0Path->SetLineColor(fMarkerColor);
  }

  void SetRnrStyle(TEveTrackPropagator* rs) { fRnrStyle = rs; }

  Float_t GetDaughterDCA() const { return fDaughterDCA; }
  void SetDaughterDCA(Float_t dca) { fDaughterDCA = dca; }

  Float_t GetRadius() const { return fRecDecayV.Perp(); }
  Float_t GetPt()     const { return fRecDecayP.Perp(); }

  Int_t GetESDIndex() const { return fESDIndex; }
  void  SetESDIndex(Int_t ind) { fESDIndex = ind;}

  virtual const Text_t* GetName() const    { return Form("ESDcascade_%i",fESDIndex); }
  virtual const Text_t* GetTitle() const   { return Form("ESDcascade_%i",fESDIndex); }

  TEveTrackPropagator* GetPropagator() const  { return fRnrStyle; }

  AliEveTrack* GetBacTrack() { return fBacTrack; }

  TEveLine*  GetPointingCurve() { return fPointingCurve; }
  TEveLine*  GetV0Path() { return fV0Path; }

protected:
  TEveVector fRecBirthV;    // Assumed birth point of cascade
  TEveVector fRecDecayV;    // Point of closest approach
  TEveVector fRecDecayP;    // Reconstructed momentum at the decay
  TEveVector fRecDecayV0;   // Reconstructed birth point of neutral daughter

  AliEveTrack        *fBacTrack;

  TEveTrackPropagator *fRnrStyle;

  TEveLine         *fPointingCurve;
  TEveLine         *fV0Path;

  Int_t             fESDIndex;    // Index in ESD V0 array.
  Float_t           fDaughterDCA; // Distance at the point of closest approach. 
  Float_t           fChi2Cascade; // Some Chi-square.

private:
  AliEveCascade(const AliEveCascade&);            // Not implemented
  AliEveCascade& operator=(const AliEveCascade&); // Not implemented

  ClassDef(AliEveCascade, 0); // Visual representation of a AliEveCascade.
};


/******************************************************************************/
// AliEveCascadeList
/******************************************************************************/

class AliEveCascadeList : public TEveElementList
{
  friend class AliEveCascadeListEditor;

public:
  AliEveCascadeList();
  AliEveCascadeList(TEveTrackPropagator* rs);
  AliEveCascadeList(const Text_t* name, TEveTrackPropagator* rs=0);
  virtual ~AliEveCascadeList() {}

  virtual const Text_t* GetTitle() const { return fTitle; }
  virtual void SetTitle(const Text_t* t) { fTitle = t; }
  virtual void SetTracksColor(Color_t cBac) { fBacColor = cBac;}

  virtual Bool_t CanEditMainColor() const { return kTRUE; }

  void  SetRnrStyle(TEveTrackPropagator* rst) { fRnrStyle = rst; }
  TEveTrackPropagator* GetPropagator()        { return fRnrStyle; }

  Bool_t GetRnrCascadevtx()     const { return fRnrCascadevtx; }
  Bool_t GetRnrCascadepath()    const { return fRnrCascadepath; }
  Bool_t GetRnrDaughters()      const { return fRnrDaughters; }

  void   MakeCascades();

  void   FilterByRadius(Float_t minR, Float_t maxR);
  void   FilterByDaughterDCA(Float_t minDaughterDCA, Float_t maxDaughterDCA);
  void   FilterByPt(Float_t minPt, Float_t maxPt);

protected:
  TString              fTitle;

  TEveTrackPropagator *fRnrStyle;

  Bool_t               fRnrDaughters;
  Bool_t               fRnrCascadevtx;
  Bool_t               fRnrCascadepath;

  Color_t              fBacColor;

  Float_t              fMinRCut;
  Float_t              fMaxRCut;

  Float_t              fMinDaughterDCA;
  Float_t              fMaxDaughterDCA;

  Float_t              fMinPt;
  Float_t              fMaxPt;

private:
  void Init();

  AliEveCascadeList(const AliEveCascadeList&);            // Not implemented
  AliEveCascadeList& operator=(const AliEveCascadeList&); // Not implemented

  ClassDef(AliEveCascadeList, 0); // A list of AliEveCascade objecs.
};


#endif
