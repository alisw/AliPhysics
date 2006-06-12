// $Header$

#ifndef REVE_Track_H
#define REVE_Track_H

#include <Reve/PODs.h>
#include <Reve/RenderElement.h>

#include <TPolyMarker3D.h>
#include <TPolyLine3D.h>

namespace Reve {

class TrackRnrStyle;
class TrackList;

class Track : public TPolyLine3D, public RenderElement
{
  friend class TrackList;

private:
  void                   Init();

protected:
  Reve::Vector           fV;
  Reve::Vector           fP;
  Double_t               fBeta;
  Int_t                  fCharge;
  Int_t                  fLabel;

  TrackRnrStyle* 	 fRnrStyle;
    
  TString                fName; 
  TString                fTitle; 

public: 
  Track();
  Track(Reve::MCTrack*  t, TrackRnrStyle* rs);
  Track(Reve::RecTrack* t, TrackRnrStyle* rs);
  std::vector<Reve::PathMark*> fPathMarks;
  virtual ~Track();

  void Reset(Int_t n_points=0);

  virtual void SetLineColor(Color_t col)
  { SetMainColor(col); }

  virtual void Paint(Option_t* option="")
  { if(fRnrElement) TPolyLine3D::Paint(option); }

  void MakeTrack();

  TrackRnrStyle* GetRnrStyle() const  { return fRnrStyle; }
  void SetRnrStyle(TrackRnrStyle* rs) { fRnrStyle = rs; }

  virtual const Text_t* GetName() const    { return fName; }
  virtual void SetName(const Text_t* name) { fName = name; }

  virtual const Text_t* GetTitle() const     { return fTitle; }
  virtual void SetTitle(const Text_t* title) { fTitle = title; }

  Int_t GetLabel() const { return fLabel; }

  //--------------------------------

  void ImportHits();     // *MENU*
  void ImportClusters(); // *MENU*

  ClassDef(Track, 1);
}; // endclass Track


/**************************************************************************/
// TrackRnrStyle
/**************************************************************************/

// This is decoupled from Track/TrackList to allow sharing of the
// RnrStyle among several instances. The interface is half cooked and
// there is no good way to set RnrStyle after the object has been
// created (shouldn't be too hard to fix).
//
// TrackList has Get/Set methods for RnrStlye and
// TrackListEditor provides editor access to them.

class TrackRnrStyle : public TObject 
{
private:
  void         Init();

public:
  Color_t                  fColor;
  Float_t                  fMagField;  
  // track limits
  Float_t                  fMaxR;       
  Float_t                  fMaxZ;       
  // helix limits 
  Float_t                  fMaxOrbs; // Maximal angular path of tracks' orbits (1 ~ 2Pi).
  Float_t                  fMinAng;  // Minimal angular step between two helix points.
  Float_t                  fDelta;   // Maximal error at the mid-point of the line connecting to helix points.

  Bool_t                   fFitDaughters;   
  Bool_t                   fFitDecay;   

  TrackRnrStyle() { Init(); }

  void    SetColor(Color_t c) { fColor = c; }
  Color_t GetColor() const    { return fColor; }

  Float_t GetMagField() const     { return fMagField; }
  void    SetMagField(Float_t mf) { fMagField = mf; }

  static Float_t       fgDefMagField;
  static const Float_t fgkB2C;
  static TrackRnrStyle fgDefStyle;

  ClassDef(TrackRnrStyle, 1);
}; // endclass TrackRnrStyle


/**************************************************************************/
// TrackList
/**************************************************************************/

class TrackList : public TPolyMarker3D, public RenderElementListBase
{
private:
  void  Init();

protected:
  TString              fTitle;

  TrackRnrStyle*       mRnrStyle;

  Bool_t               fRnrMarkers;
  Bool_t               fRnrTracks;

public:
  TrackList(Int_t n_tracks=0);
  TrackList(const Text_t* name, Int_t n_tracks=0);

  void Reset(Int_t n_tracks=0);

  virtual const Text_t* GetTile() const  { return fTitle; }
  virtual void SetTitle(const Text_t* t) { fTitle = t; }

  virtual Bool_t CanEditMainColor()  { return kTRUE; }

  virtual void Paint(Option_t* option="");

  virtual void AddElement(RenderElement* el);

  void  SetRnrStyle(TrackRnrStyle* rst) { mRnrStyle= rst; }
  TrackRnrStyle* GetRnrStyle()          { return mRnrStyle; } 

  Bool_t GetRnrTracks() const { return fRnrTracks; }
  void   SetRnrTracks(Bool_t);

  Bool_t GetRnrMarkers() const { return fRnrMarkers; }
  void   SetRnrMarkers(Bool_t);

  void   MakeTracks();
  void   MakeMarkers();

  Float_t GetMaxR()         const { return mRnrStyle->fMaxZ; }
  Float_t GetMaxZ()         const { return mRnrStyle->fMaxR; }
  Float_t GetMaxOrbs()      const { return mRnrStyle->fMaxOrbs; }
  Float_t GetMinAng()       const { return mRnrStyle->fMinAng; }
  Float_t GetDelta()        const { return mRnrStyle->fDelta; }
  Bool_t  GetFitDaughters() const { return mRnrStyle->fFitDaughters; }
  Bool_t  GetFitDecay()     const { return mRnrStyle->fFitDecay; }

  void SetMaxR(Float_t x);
  void SetMaxZ(Float_t x);
  void SetMaxOrbs(Float_t x);
  void SetMinAng(Float_t x);
  void SetDelta(Float_t x);
  void SetFitDaughters(Bool_t x);
  void SetFitDecay(Bool_t x);


  // void  UpdateBounds();
  Int_t   GetNTracks() { return fN; }

  void SelectByPt(Float_t min_pt=0.2, Float_t max_pt=10); // *MENU*

  void MakePtScrollbar(); // *MENU*
  void HandlePtScrollEvent();

  ClassDef(TrackList, 1);
};


} // namespace Reve

#endif
