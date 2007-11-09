#ifndef REVE_Track_H
#define REVE_Track_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <vector>

#include <Reve/PODs.h>
#include <Reve/RenderElement.h>
#include <Reve/Line.h>

#include <TPolyMarker3D.h>
#include <TMarker.h>
#include <TQObject.h>

namespace Reve {

class TrackRnrStyle;
class TrackList;

class Track : public Line, 
              public TQObject
{
  friend class TrackRnrStyle;
  friend class TrackList;
  friend class TrackCounter;
  friend class TrackGL;

public:
  typedef std::vector<Reve::PathMark*>           vpPathMark_t;
  typedef std::vector<Reve::PathMark*>::iterator vpPathMark_i;

protected:
  Reve::Vector      fV;          // Starting vertex
  Reve::Vector      fP;          // Starting momentum
  Double_t          fBeta;       // Relativistic beta factor
  Int_t             fPdg;        // PDG code
  Int_t             fCharge;     // Charge in units of e0
  Int_t             fLabel;      // Simulation label
  Int_t             fIndex;      // Reconstruction index
  vpPathMark_t      fPathMarks;  // Vector of known points along the track

  TrackRnrStyle*    fRnrStyle;   // Pointer to shared render-style

public:
  Track();
  Track(TParticle* t, Int_t label, TrackRnrStyle* rs);
  Track(Reve::MCTrack*  t, TrackRnrStyle* rs);
  Track(Reve::RecTrack* t, TrackRnrStyle* rs);
  Track(const Track& t);
  virtual ~Track();

  virtual void SetStdTitle();

  virtual void SetTrackParams(const Track& t);
  virtual void SetPathMarks  (const Track& t);

  virtual void MakeTrack(Bool_t recurse=kTRUE);

  TrackRnrStyle* GetRnrStyle() const  { return fRnrStyle; }
  void SetRnrStyle(TrackRnrStyle* rs);
  void SetAttLineAttMarker(TrackList* tl);

  Int_t GetPdg()    const   { return fPdg;   }
  void SetPdg(Int_t pdg)    { fPdg = pdg;    }
  Int_t GetCharge() const   { return fCharge; }
  void SetCharge(Int_t chg) { fCharge = chg; }
  Int_t GetLabel()  const   { return fLabel; }
  void  SetLabel(Int_t lbl) { fLabel = lbl;  }
  Int_t GetIndex()  const   { return fIndex; }
  void  SetIndex(Int_t idx) { fIndex = idx;  }

  void          AddPathMark(Reve::PathMark* pm) { fPathMarks.push_back(pm); }
  vpPathMark_t& GetPathMarksRef()               { return fPathMarks; }
  const vpPathMark_t& GetPathMarksRef() const   { return fPathMarks; }
  void          SortPathMarksByTime();

  //--------------------------------

  void ImportHits();              // *MENU*
  void ImportClusters();          // *MENU*
  void ImportClustersFromIndex(); // *MENU*
  void ImportKine();              // *MENU*
  void ImportKineWithArgs(Bool_t importMother=kTRUE, Bool_t impDaugters=kTRUE,
			  Bool_t colorPdg    =kTRUE, Bool_t recurse    =kTRUE); // *MENU*
  void PrintKineStack();          // *MENU*
  void PrintPathMarks();          // *MENU*

  //--------------------------------

  virtual void CtrlClicked(Reve::Track*); // *SIGNAL*
  virtual void SetLineStyle(Style_t lstyle);

  virtual const TGPicture* GetListTreeIcon() { return fgListTreeIcons[4]; }; 

  virtual TClass* ProjectedClass() const;

  ClassDef(Track, 1); // Visual representation of a track.
}; // endclass Track


/**************************************************************************/
// TrackRnrStyle
/**************************************************************************/

class TrackRnrStyle : public TObject,
                      public ReferenceBackPtr
{
private:
  void                     RebuildTracks();

public:
  Float_t                  fMagField;      // Constant magnetic field along z.

  // Track limits
  Float_t                  fMaxR;          // Max radius for track extrapolation
  Float_t                  fMaxZ;          // Max z-coordinate for track extrapolation.
  // Helix limits
  Float_t                  fMaxOrbs;       // Maximal angular path of tracks' orbits (1 ~ 2Pi).
  Float_t                  fMinAng;        // Minimal angular step between two helix points.
  Float_t                  fDelta;         // Maximal error at the mid-point of the line connecting to helix points.

  // Path-mark control
  Bool_t                   fEditPathMarks; // Show widgets for path-mark control in GUI editor.
  TMarker                  fPMAtt;         // Marker attributes for rendering of path-marks.

  Bool_t                   fFitDaughters;  // Pass through daughter creation points when extrapolating a track.
  Bool_t                   fFitReferences; // Pass through given track-references when extrapolating a track.
  Bool_t                   fFitDecay;      // Pass through decay point when extrapolating a track.

  Bool_t                   fRnrDaughters;  // Render daughter path-marks.
  Bool_t                   fRnrReferences; // Render track-reference path-marks.
  Bool_t                   fRnrDecay;      // Render decay path-marks.

  // First vertex control
  Bool_t                   fRnrFV;         // Render first vertex.
  TMarker                  fFVAtt;         // Marker attributes for fits vertex.

  TrackRnrStyle();

  // callbacks
  void   SetEditPathMarks(Bool_t x) { fEditPathMarks = x; }
  void   SetRnrDaughters(Bool_t x);
  void   SetRnrReferences(Bool_t x);
  void   SetRnrDecay(Bool_t x);

  void   SetRnrFV(Bool_t x){  fRnrFV = x;}

  void   SetFitDaughters(Bool_t x);
  void   SetFitReferences(Bool_t x);
  void   SetFitDecay(Bool_t x);

  void   SetMaxR(Float_t x);
  void   SetMaxZ(Float_t x);
  void   SetMaxOrbs(Float_t x);
  void   SetMinAng(Float_t x);
  void   SetDelta(Float_t x);

  Float_t GetMagField() const     { return fMagField; }
  void    SetMagField(Float_t mf) { fMagField = mf; }

  static Float_t       fgDefMagField; // Default value for constant solenoid magnetic field.
  static const Float_t fgkB2C;        // Constant for conversion of momentum to curvature.
  static TrackRnrStyle fgDefStyle;    // Default track render-style.

  ClassDef(TrackRnrStyle, 1); // Rendering parameters for tracks.
}; // endclass TrackRnrStyle


/**************************************************************************/
// TrackList
/**************************************************************************/

class TrackList : public RenderElementList,
                  public NLTProjectable,
                  public TAttMarker,
                  public TAttLine
{
  friend class TrackListEditor;

private:
  TrackList(const TrackList&);            // Not implemented
  TrackList& operator=(const TrackList&); // Not implemented

  Bool_t               fRecurse;    // Recurse when propagating marker/line attributes to tracks.

protected:
  TrackRnrStyle*       fRnrStyle;   // Basic track rendering parameters, not enforced to elements.

  Bool_t               fRnrLine;    // Render track as line.
  Bool_t               fRnrPoints;  // Render track as points.

  Float_t              fMinPt;      // Minimum track pT for display selection.
  Float_t              fMaxPt;      // Maximum track pT for display selection.
  Float_t              fLimPt;      // Highest track pT in the container.
  Float_t              fMinP;       // Minimum track p for display selection.
  Float_t              fMaxP;       // Maximum track p for display selection.
  Float_t              fLimP;       // Highest track p in the container.

  Float_t RoundMomentumLimit(Float_t x);

public:
  TrackList(TrackRnrStyle* rs=0);
  TrackList(const Text_t* name, TrackRnrStyle* rs=0);
  virtual ~TrackList();

  void  MakeTracks(Bool_t recurse=kTRUE);
  void  FindMomentumLimits(RenderElement* el, Bool_t recurse);

  void  SetRnrStyle(TrackRnrStyle* rs);
  TrackRnrStyle*  GetRnrStyle(){return fRnrStyle;}

  //--------------------------------

  virtual void   SetMainColor(Color_t c);
  virtual void   SetLineColor(Color_t c){SetMainColor(c);}
  virtual void   SetLineColor(Color_t c, RenderElement* el);
  virtual void   SetLineWidth(Width_t w);
  virtual void   SetLineWidth(Width_t w, RenderElement* el);
  virtual void   SetLineStyle(Style_t s);
  virtual void   SetLineStyle(Style_t s, RenderElement* el);

  virtual void   SetMarkerColor(Color_t c);
  virtual void   SetMarkerColor(Color_t c, RenderElement* el);
  virtual void   SetMarkerSize(Size_t s);
  virtual void   SetMarkerSize(Size_t s, RenderElement* el);
  virtual void   SetMarkerStyle(Style_t s);
  virtual void   SetMarkerStyle(Style_t s, RenderElement* el);

  void SetRnrLine(Bool_t rnr);
  void SetRnrLine(Bool_t rnr, RenderElement* el);
  Bool_t GetRnrLine(){return fRnrLine;}

  void SetRnrPoints(Bool_t r);
  void SetRnrPoints(Bool_t r, RenderElement* el);
  Bool_t GetRnrPoints(){return fRnrPoints;}

  void SelectByPt(Float_t min_pt, Float_t max_pt);
  void SelectByPt(Float_t min_pt, Float_t max_pt, RenderElement* el);
  void SelectByP (Float_t min_p,  Float_t max_p);
  void SelectByP (Float_t min_p,  Float_t max_p,  RenderElement* el);

  //--------------------------------

  Track* FindTrackByLabel(Int_t label); // *MENU*
  Track* FindTrackByIndex(Int_t index); // *MENU*

  void ImportHits();     // *MENU*
  void ImportClusters(); // *MENU*

  virtual TClass* ProjectedClass() const;

  ClassDef(TrackList, 1); // A list of tracks.
};


/**************************************************************************/
// TrackCounter
/**************************************************************************/

class TrackCounter : public RenderElement, public TNamed
{
  friend class TrackCounterEditor;

public:
  enum ClickAction_e { CA_PrintTrackInfo, CA_ToggleTrack };

private:
  TrackCounter(const TrackCounter&);            // Not implemented
  TrackCounter& operator=(const TrackCounter&); // Not implemented

protected:
  Int_t fBadLineStyle;
  Int_t fClickAction;

  Int_t fEventId;

  Int_t fAllTracks;
  Int_t fGoodTracks;

  TList fTrackLists;

public:
  TrackCounter(const Text_t* name="TrackCounter", const Text_t* title="");
  virtual ~TrackCounter();

  Int_t GetEventId() const { return fEventId; }
  void  SetEventId(Int_t id) { fEventId = id; }

  void Reset();

  void RegisterTracks(TrackList* tlist, Bool_t goodTracks);

  void DoTrackAction(Track* track);

  Int_t GetClickAction() const  { return fClickAction; }
  void  SetClickAction(Int_t a) { fClickAction = a; }

  void OutputEventTracks(FILE* out=0);

  static TrackCounter* fgInstance;

  ClassDef(TrackCounter, 1);
}; // endclass TrackCounter


} // namespace Reve

#endif
