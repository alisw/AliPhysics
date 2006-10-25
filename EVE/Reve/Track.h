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

class Track : public RenderElement,
	      public TPolyLine3D
{
  friend class TrackList;

  Track(const Track&);            // Not implemented
  Track& operator=(const Track&); // Not implemented

protected:
  typedef std::vector<Reve::PathMark*>           vpPathMark_t;
  typedef std::vector<Reve::PathMark*>::iterator vpPathMark_i;

  Reve::Vector      fV;
  Reve::Vector      fP;
  Double_t          fBeta;
  Int_t             fCharge;
  Int_t             fLabel;
  vpPathMark_t      fPathMarks;

  TrackRnrStyle*    fRnrStyle;
    
  TString           fName;
  TString           fTitle;

public: 
  Track();
  Track(TParticle* t, Int_t label, TrackRnrStyle* rs);
  Track(Reve::MCTrack*  t, TrackRnrStyle* rs);
  Track(Reve::RecTrack* t, TrackRnrStyle* rs);
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
  void  AddPathMark(Reve::PathMark* pm) { fPathMarks.push_back(pm); }

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
public:
  Color_t                  fColor;
  Width_t                  fWidth;
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

  TrackRnrStyle();

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

class TrackList : public RenderElementListBase,
		  public TPolyMarker3D
{
  TrackList(const TrackList&);            // Not implemented
  TrackList& operator=(const TrackList&); // Not implemented

private:
  void  Init();

protected:
  TString              fTitle;

  TrackRnrStyle*       fRnrStyle;

  Bool_t               fRnrMarkers;
  Bool_t               fRnrTracks;

public:
  TrackList(Int_t n_tracks=0, TrackRnrStyle* rs=0);
  TrackList(const Text_t* name, Int_t n_tracks=0, TrackRnrStyle* rs=0);

  void Reset(Int_t n_tracks=0);

  virtual const Text_t* GetTitle() const { return fTitle; }
  virtual void SetTitle(const Text_t* t) { fTitle = t; }

  virtual Bool_t CanEditMainColor()  { return kTRUE; }

  virtual void Paint(Option_t* option="");

  virtual void AddElement(RenderElement* el);

  void  SetRnrStyle(TrackRnrStyle* rst) { fRnrStyle= rst; }
  TrackRnrStyle* GetRnrStyle()          { return fRnrStyle; } 

  Bool_t GetRnrTracks() const { return fRnrTracks; }
  void   SetRnrTracks(Bool_t);

  Bool_t GetRnrMarkers() const { return fRnrMarkers; }
  void   SetRnrMarkers(Bool_t);

  void   MakeTracks();
  void   MakeMarkers();

  Width_t GetWidth() const { return fRnrStyle->fWidth; }
  void SetWidth(Width_t w);

  Float_t GetMaxR()         const { return fRnrStyle->fMaxZ; }
  Float_t GetMaxZ()         const { return fRnrStyle->fMaxR; }
  Float_t GetMaxOrbs()      const { return fRnrStyle->fMaxOrbs; }
  Float_t GetMinAng()       const { return fRnrStyle->fMinAng; }
  Float_t GetDelta()        const { return fRnrStyle->fDelta; }
  Bool_t  GetFitDaughters() const { return fRnrStyle->fFitDaughters; }
  Bool_t  GetFitDecay()     const { return fRnrStyle->fFitDecay; }

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

  //--------------------------------

  void ImportHits();     // *MENU*
  void ImportClusters(); // *MENU*

  ClassDef(TrackList, 1);
};


} // namespace Reve

#endif
