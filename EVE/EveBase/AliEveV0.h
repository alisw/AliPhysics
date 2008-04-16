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

public:
  AliEveV0();
  AliEveV0(TEveRecTrack* tNeg, TEveRecTrack* tPos, TEveRecV0* v0,
     TEveTrackPropagator* rs);
  virtual ~AliEveV0();

  //void  AddPathMarkPos(TEvePathMark* pm) { fPathMarksPos.push_back(pm); }
  //void  AddPathMarkNeg(TEvePathMark* pm) { fPathMarksNeg.push_back(pm); }

  void Reset(TPolyLine3D* polyLine);

  void MakeV0path();
  void MakeV0();

  virtual void PaintDaughters(Option_t* option="")
  {
    if (fRnrSelf) { fNegTrack->Paint(option); fPosTrack->Paint(option);}
  }

  virtual void Paint(Option_t* option="")
  {
    if (fRnrSelf) TEvePointSet::Paint(option);
  }

  virtual void PaintPath(Option_t* option="")
  {
    if (fRnrSelf) fPolyLineV0.Paint(option);
  }

  virtual void  SetMainColor(Color_t col)
  {
    fMarkerColor = col; fMainColorPtr = &fMarkerColor;
    fPolyLineV0.SetLineColor(fMarkerColor);
  }

  void SetRnrStyle(TEveTrackPropagator* rs) { fRnrStyle = rs; }

  Float_t GetDaughterDCA() const { return fDaughterDCA; }
  void SetDaughterDCA(Float_t dca) { fDaughterDCA = dca; }

  Float_t GetRadius() const { return fRecDecayV.Perp(); }

  Int_t GetESDIndex() const { return fESDIndex; }
  void  SetESDIndex(Int_t ind) { fESDIndex = ind;}

  virtual const Text_t* GetName() const    { return Form("ESDv0_%i",fESDIndex); }
  virtual const Text_t* GetTitle() const   { return Form("ESDv0_%i",fESDIndex); }

  //Int_t          GetLabelPos() const { return fLabel_pos; }
  //Int_t          GetLabelNeg() const { return fLabel_neg; }
  TEveTrackPropagator* GetPropagator() const  { return fRnrStyle; }

  TEveTrack*   GetNegTrack() { return fNegTrack; }
  TEveTrack*   GetPosTrack() { return fPosTrack; }

  TPolyLine3D*  GetPolyLineV0() { return &fPolyLineV0; }

protected:
  TEveVector fRecBirthV;    // Reconstucted birth point of neutral particle
  TEveVector fRecDecayV;    // Point of closest approach
  TEveVector fRecDecayP;

  TEveTrack        *fNegTrack;
  TEveTrack        *fPosTrack;

  TEveTrackPropagator *fRnrStyle;

  TPolyLine3D       fPolyLineV0;

  Int_t             fESDIndex;
  Float_t           fDaughterDCA;
  Float_t           fChi2V0;

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
public:
  AliEveV0List();
  AliEveV0List(TEveTrackPropagator* rs);
  AliEveV0List(const Text_t* name, TEveTrackPropagator* rs=0);
  virtual ~AliEveV0List();

  virtual const Text_t* GetTitle() const { return fTitle; }
  virtual void SetTitle(const Text_t* t) { fTitle = t; }
  virtual void SetTracksColor(Color_t cNeg, Color_t cPos) {
    fNegColor = cNeg; fPosColor = cPos;}

  virtual Bool_t CanEditMainColor() const { return kTRUE; }

  virtual void Paint(Option_t* option="");

  void  SetRnrStyle(TEveTrackPropagator* rst) { fRnrStyle = rst; }
  TEveTrackPropagator* GetPropagator()        { return fRnrStyle; }

  Bool_t GetRnrV0vtx()     const { return fRnrV0vtx; }
  Bool_t GetRnrV0path()    const { return fRnrV0path; }
  Bool_t GetRnrDaughters() const { return fRnrDaughters; }

  void   SetRnrV0vtx(Bool_t);
  void   SetRnrV0path(Bool_t);
  void   SetRnrDaughters(Bool_t);

  void   MakeV0s();
  void   MakeMarkers();

protected:
  TString              fTitle;

  TEveTrackPropagator *fRnrStyle;

  Bool_t               fRnrDaughters;
  Bool_t               fRnrV0vtx;
  Bool_t               fRnrV0path;

  Color_t              fNegColor;
  Color_t              fPosColor;

private:
  void  Init();

  AliEveV0List(const AliEveV0List&);            // Not implemented
  AliEveV0List& operator=(const AliEveV0List&); // Not implemented

  ClassDef(AliEveV0List, 0); // A list of AliEveV0 objecs.
};


#endif
