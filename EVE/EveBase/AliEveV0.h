// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/
#ifndef AliEveV0_H
#define AliEveV0_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/***********************************************************************
* This code defines the reconstructed v0 visualized with EVE
*
* Ludovic Gaudichet (gaudichet@to.infn.it)
************************************************************************/

#include <TEveVSDStructs.h>
#include <TEveElement.h>
#include <TEveTrack.h>

#include <TPolyMarker3D.h>
#include <TPolyLine3D.h>

class TH1F;
class TH2F;


class AliEveV0List;

class AliEveV0 : public TEvePointSet
{
  friend class AliEveV0List;
  friend class AliEveV0Editor;

public:
  AliEveV0();
  AliEveV0(TEveRecTrack* tNeg, TEveRecTrack* tPos, TEveRecV0* v0,
     TEveTrackPropagator* rs);
  virtual ~AliEveV0();

  void MakeV0();

  virtual void  SetMainColor(Color_t col)
  {
    TEvePointSet::SetMainColor(col);
    fPointingLine->SetLineColor(fMarkerColor);
  }

  void SetRnrStyle(TEveTrackPropagator* rs) { fRnrStyle = rs; }

  Float_t GetDaughterDCA() const { return fDaughterDCA; }
  void SetDaughterDCA(Float_t dca) { fDaughterDCA = dca; }

  Float_t GetRadius() const { return fRecDecayV.Perp(); }
  Float_t GetPt()     const { return fRecDecayP.Perp(); }

  Bool_t GetOnFlyStatus()    const { return fOnFlyStatus; }
  void   SetOnFlyStatus(Bool_t fs) { fOnFlyStatus = fs; }

  Int_t GetESDIndex() const { return fESDIndex; }
  void  SetESDIndex(Int_t ind) { fESDIndex = ind;}

  virtual const Text_t* GetName() const    { return Form("ESDv0_%i",fESDIndex); }
  virtual const Text_t* GetTitle() const   { return Form("ESDv0_%i",fESDIndex); }

  TEveTrackPropagator* GetPropagator() const  { return fRnrStyle; }

  TEveTrack* GetNegTrack() { return fNegTrack; }
  TEveTrack* GetPosTrack() { return fPosTrack; }

  TEveLine*  GetPointingLine() { return fPointingLine; }

protected:
  TEveVector fRecBirthV;    // Reconstucted birth point of neutral particle
  TEveVector fRecDecayV;    // Point of closest approach
  TEveVector fRecDecayP;

  TEveTrack        *fNegTrack;
  TEveTrack        *fPosTrack;

  TEveTrackPropagator *fRnrStyle;

  TEveLine         *fPointingLine;

  Int_t             fESDIndex;    // Index in ESD V0 array.
  Bool_t            fOnFlyStatus; // Reconstructed during tracking.
  Float_t           fDaughterDCA; // Distance at the point of closest approach. 
  Float_t           fChi2V0;      // Some Chi-square.

private:
  AliEveV0(const AliEveV0&);            // Not implemented
  AliEveV0& operator=(const AliEveV0&); // Not implemented

  ClassDef(AliEveV0, 0); // Visual representation of a AliEveV0.
};


/******************************************************************************/
// AliEveV0List
/******************************************************************************/

class AliEveV0List : public TEveElementList
{
  friend class AliEveV0ListEditor;

public:
  AliEveV0List();
  AliEveV0List(TEveTrackPropagator* rs);
  AliEveV0List(const Text_t* name, TEveTrackPropagator* rs=0);
  virtual ~AliEveV0List() {}

  virtual const Text_t* GetTitle() const { return fTitle; }
  virtual void SetTitle(const Text_t* t) { fTitle = t; }
  virtual void SetTracksColor(Color_t cNeg, Color_t cPos) {
    fNegColor = cNeg; fPosColor = cPos;}

  virtual Bool_t CanEditMainColor() const { return kTRUE; }

  void  SetRnrStyle(TEveTrackPropagator* rst) { fRnrStyle = rst; }
  TEveTrackPropagator* GetPropagator()        { return fRnrStyle; }

  Bool_t GetRnrV0vtx()     const { return fRnrV0vtx; }
  Bool_t GetRnrV0path()    const { return fRnrV0path; }
  Bool_t GetRnrDaughters() const { return fRnrDaughters; }

  void   MakeV0s();

  void   FilterByRadius(Float_t minR, Float_t maxR);
  void   FilterByDaughterDCA(Float_t minDaughterDCA, Float_t maxDaughterDCA);
  void   FilterByPt(Float_t minPt, Float_t maxPt);

protected:
  TString              fTitle;

  TEveTrackPropagator *fRnrStyle;

  Bool_t               fRnrDaughters;
  Bool_t               fRnrV0vtx;
  Bool_t               fRnrV0path;

  Color_t              fNegColor;
  Color_t              fPosColor;

  Float_t              fMinRCut;
  Float_t              fMaxRCut;

  Float_t              fMinDaughterDCA;
  Float_t              fMaxDaughterDCA;

  Float_t              fMinPt;
  Float_t              fMaxPt;

private:
  void Init();

  AliEveV0List(const AliEveV0List&);            // Not implemented
  AliEveV0List& operator=(const AliEveV0List&); // Not implemented

  ClassDef(AliEveV0List, 0); // A list of AliEveV0 objecs.
};


#endif
