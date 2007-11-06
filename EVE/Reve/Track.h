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
  Reve::Vector      fV;
  Reve::Vector      fP;
  Double_t          fBeta;
  Int_t             fCharge;
  Int_t             fLabel;
  Int_t             fIndex;
  vpPathMark_t      fPathMarks;

  TrackRnrStyle*    fRnrStyle; 

public:
  Track();
  Track(TParticle* t, Int_t label, TrackRnrStyle* rs);
  Track(Reve::MCTrack*  t, TrackRnrStyle* rs);
  Track(Reve::RecTrack* t, TrackRnrStyle* rs);
  virtual ~Track();

  Track(const Track& t);
  virtual void SetTrackParams(const Track& t);
  virtual void SetPathMarks  (const Track& t);

  virtual void MakeTrack(Bool_t recurse=kTRUE);

  TrackRnrStyle* GetRnrStyle() const  { return fRnrStyle; }
  void SetRnrStyle(TrackRnrStyle* rs);
  void SetAttLineAttMarker(TrackList* tl);

  Int_t GetLabel() const    { return fLabel; }
  void  SetLabel(Int_t lbl) { fLabel = lbl;  }
  Int_t GetIndex() const    { return fIndex; }
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
  void ImportKineWithArgs(Bool_t importMother=kTRUE,
			  Bool_t impDaugters =kTRUE); // *MENU*
  void PrintKineStack();          // *MENU*
  void PrintPathMarks();          // *MENU*

  //--------------------------------

  virtual void CtrlClicked(Reve::Track*); // *SIGNAL*
  virtual void SetLineStyle(Style_t lstyle);

  virtual const TGPicture* GetListTreeIcon() { return fgListTreeIcons[4]; }; 

  virtual TClass* ProjectedClass() const;

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

class TrackRnrStyle : public TObject,
                      public ReferenceBackPtr
{
private:
  void                     RebuildTracks();

public:
  Float_t                  fMagField;
  // track limits
  Float_t                  fMaxR;
  Float_t                  fMaxZ;
  // helix limits
  Float_t                  fMaxOrbs; // Maximal angular path of tracks' orbits (1 ~ 2Pi).
  Float_t                  fMinAng;  // Minimal angular step between two helix points.
  Float_t                  fDelta;   // Maximal error at the mid-point of the line connecting to helix points.

  Bool_t                   fEditPathMarks;
  TMarker                  fPMAtt;

  Bool_t                   fFitDaughters;
  Bool_t                   fFitReferences;
  Bool_t                   fFitDecay;

  Bool_t                   fRnrDaughters;
  Bool_t                   fRnrReferences;
  Bool_t                   fRnrDecay;
 
  Bool_t                   fRnrFV; // first vertex
  TMarker                  fFVAtt;

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

  static Float_t       fgDefMagField;
  static const Float_t fgkB2C;
  static TrackRnrStyle fgDefStyle;

  ClassDef(TrackRnrStyle, 1);
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
  Bool_t               fRecurse;

protected:
  TrackRnrStyle*       fRnrStyle;

  Bool_t               fRnrLine;
  Bool_t               fRnrPoints;

  Float_t              fMinPt;
  Float_t              fMaxPt;
  Float_t              fLimPt;
  Float_t              fMinP;
  Float_t              fMaxP;
  Float_t              fLimP;

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

  void ImportHits();     // *MENU*
  void ImportClusters(); // *MENU*

  virtual TClass* ProjectedClass() const;

  ClassDef(TrackList, 1);
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
