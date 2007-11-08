// $Header$

#include "Track.h"
#include "MCHelixLine.hi"
#include "PointSet.h"

#include <TPolyLine3D.h>
#include <TMarker.h>
#include <TPolyMarker3D.h>
#include <TColor.h>

// Updates
#include <Reve/ReveManager.h>
#include <Reve/RGBrowser.h>
#include <Reve/NLTTrack.h>
#include <TCanvas.h>

#include <vector>
#include <algorithm>
#include <functional>
#include <iostream>

using namespace Reve;

//______________________________________________________________________________
// Track
//
// Visual representation of a track.
//

ClassImp(Reve::Track)

//______________________________________________________________________________
Track::Track() :
  Line(),

  fV(),
  fP(),
  fBeta(0),
  fCharge(0),
  fLabel(kMinInt),
  fIndex(kMinInt),
  fPathMarks(),

  fRnrStyle(0)
{}

//______________________________________________________________________________
Track::Track(TParticle* t, Int_t label, TrackRnrStyle* rs):
  Line(),

  fV(t->Vx(), t->Vy(), t->Vz()),
  fP(t->Px(), t->Py(), t->Pz()),
  fBeta(t->P()/t->Energy()),
  fCharge(0),
  fLabel(label),
  fIndex(kMinInt),
  fPathMarks(),

  fRnrStyle(0)
{
  SetRnrStyle(rs);
  fMainColorPtr = &fLineColor;

  TParticlePDG* pdgp = t->GetPDG();
  if (pdgp)
    fCharge = (Int_t) TMath::Nint(pdgp->Charge()/3);

  SetName(t->GetName());
}

//______________________________________________________________________________
Track::Track(Reve::MCTrack* t, TrackRnrStyle* rs):
  Line(),

  fV(t->Vx(), t->Vy(), t->Vz()),
  fP(t->Px(), t->Py(), t->Pz()),
  fBeta(t->P()/t->Energy()),
  fCharge(0),
  fLabel(t->label),
  fIndex(t->index),
  fPathMarks(),

  fRnrStyle(0)
{
  SetRnrStyle(rs);
  fMainColorPtr = &fLineColor;

  TParticlePDG* pdgp = t->GetPDG();
  if(pdgp == 0) {
    t->ResetPdgCode(); pdgp = t->GetPDG();
  }
  fCharge = (Int_t) TMath::Nint(pdgp->Charge()/3);

  SetName(t->GetName());
}

//______________________________________________________________________________
Track::Track(Reve::RecTrack* t, TrackRnrStyle* rs) :
  Line(),

  fV(t->V),
  fP(t->P),
  fBeta(t->beta),
  fCharge(t->sign),
  fLabel(t->label),
  fIndex(t->index),
  fPathMarks(),

  fRnrStyle(0)
{
  SetRnrStyle(rs);
  fMainColorPtr = &fLineColor;

  SetName(t->GetName());
}

//______________________________________________________________________________
Track::~Track()
{
  SetRnrStyle(0);
  for (vpPathMark_i i=fPathMarks.begin(); i!=fPathMarks.end(); ++i)
    delete *i;
}

//______________________________________________________________________________
Track::Track(const Track& t) :
  Line(),
  TQObject(),
  fV(t.fV),
  fP(t.fP),
  fBeta(t.fBeta),
  fCharge(t.fCharge),
  fLabel(t.fLabel),
  fIndex(t.fIndex),
  fPathMarks(),
  fRnrStyle(0)
{
  // Copy constructor.

  SetMainColor(t.GetMainColor());
  // Line
  fRnrLine   = t.fRnrLine;
  fRnrPoints = t.fRnrPoints;
  // TLineAttrib
  fLineColor = t.fLineColor;
  fLineStyle = t.fLineStyle;
  fLineWidth = t.fLineWidth;
  SetPathMarks(t);
  SetRnrStyle (t.fRnrStyle);
}

//______________________________________________________________________________
void Track::SetTrackParams(const Track& t)
{
  // Copy track parameters from t.
  // PathMarks are cleared.

  fV         = t.fV;
  fP         = t.fP;
  fBeta      = t.fBeta;
  fCharge    = t.fCharge;
  fLabel     = t.fLabel;
  fIndex     = t.fIndex;

  SetMainColor(t.GetMainColor());
  // Line
  fRnrLine   = t.fRnrLine;
  fRnrPoints = t.fRnrPoints;
  // TLineAttrib
  fLineColor = t.fLineColor;
  fLineStyle = t.fLineStyle;
  fLineWidth = t.fLineWidth;
  fPathMarks.clear();
  SetRnrStyle(t.fRnrStyle);
}

//______________________________________________________________________________
void Track::SetPathMarks(const Track& t)
{
  // Copy path-marks from t.

  const std::vector<PathMark*>& refs = t.GetPathMarksRef();
  for(std::vector<PathMark*>::const_iterator i=refs.begin(); i!=refs.end(); ++i)
  {
    fPathMarks.push_back(new PathMark(**i));
  }
}

/******************************************************************************/

//______________________________________________________________________________
void Track::SetRnrStyle(TrackRnrStyle* rs)
{
 if (fRnrStyle == rs) return;
  if (fRnrStyle) fRnrStyle->DecRefCount(this);
  fRnrStyle = rs;
  if (fRnrStyle) rs->IncRefCount(this);
}

/******************************************************************************/

//______________________________________________________________________________
void Track::SetAttLineAttMarker(TrackList* tl)
{
  SetLineColor(tl->GetLineColor());
  SetLineStyle(tl->GetLineStyle());
  SetLineWidth(tl->GetLineWidth());

  SetMarkerColor(tl->GetMarkerColor());
  SetMarkerStyle(tl->GetMarkerStyle());
  SetMarkerSize(tl->GetMarkerSize());
}

/******************************************************************************/

//______________________________________________________________________________
void Track::MakeTrack(Bool_t recurse)
{
  TrackRnrStyle& RS((fRnrStyle != 0) ? *fRnrStyle : TrackRnrStyle::fgDefStyle);

  Float_t px = fP.x, py = fP.y, pz = fP.z;  

  MCVertex  mc_v0;
  mc_v0.x = fV.x;
  mc_v0.y = fV.y; 
  mc_v0.z = fV.z; 
  mc_v0.t = 0;

  std::vector<MCVertex> track_points;
  Bool_t decay = kFALSE;

  if ((TMath::Abs(fV.z) > RS.fMaxZ) || (fV.x*fV.x + fV.y*fV.y > RS.fMaxR*RS.fMaxR)) 
    goto make_polyline;
  
  if (fCharge != 0 && TMath::Abs(RS.fMagField) > 1e-5 && fP.Perp2() > 1e-6)
  {
    // Charged particle in magnetic field with non-zero pT.

    Float_t a = RS.fgkB2C * RS.fMagField * fCharge;
   
    MCHelix helix(fRnrStyle, &mc_v0, TMath::C()*fBeta, &track_points, a); //m->cm
    helix.Init(TMath::Sqrt(px*px+py*py), pz);
    // Set max number of points for loop-to-vertex.
    // loop-to-bounds (last step) does this separately.
    helix.NMax = 4096;
   
    if (!fPathMarks.empty())
    {
      for(std::vector<Reve::PathMark*>::iterator i=fPathMarks.begin(); i!=fPathMarks.end(); ++i)
      {
	Reve::PathMark* pm = *i;
        
	if (RS.fFitReferences && pm->type == Reve::PathMark::Reference)
	{
	  if(TMath::Abs(pm->V.z) > RS.fMaxZ ||
	     TMath::Sqrt(pm->V.x*pm->V.x + pm->V.y*pm->V.y) > RS.fMaxR)
	    goto helix_bounds;

	  // printf("%s fit reference  \n", fName.Data()); 
	  helix.LoopToVertex(px, py, pz, pm->V.x, pm->V.y, pm->V.z);
	  px = pm->P.x;
	  py = pm->P.y;
	  pz = pm->P.z;
	}
	else if(RS.fFitDaughters &&  pm->type == Reve::PathMark::Daughter)
	{
	  if(TMath::Abs(pm->V.z) > RS.fMaxZ ||
	     TMath::Sqrt(pm->V.x*pm->V.x + pm->V.y*pm->V.y) > RS.fMaxR)
	    goto helix_bounds;

          // printf("%s fit daughter  \n", fName.Data()); 
	  helix.LoopToVertex(px, py, pz, pm->V.x, pm->V.y, pm->V.z);
	  px -= pm->P.x;
	  py -= pm->P.y;
	  pz -= pm->P.z;
	}
	else if(RS.fFitDecay &&  pm->type == Reve::PathMark::Decay)
	{
	  if(TMath::Abs(pm->V.z) > RS.fMaxZ ||
	     TMath::Sqrt(pm->V.x*pm->V.x + pm->V.y*pm->V.y) > RS.fMaxR)
	    goto helix_bounds;
	  helix.LoopToVertex(px, py, pz, pm->V.x, pm->V.y, pm->V.z);
          decay = true;
          break;
	}
	if (track_points.size() > 4096)
	{
	  Warning("Track::MakeTrack", "exceeding 4k points (%u) for '%s'; aborting extrapolation.",
		  track_points.size(), GetName());
	  goto make_polyline;
	}
      }
    }
  helix_bounds:
    // go to bounds
    if(!decay || RS.fFitDecay == kFALSE)
    {
      helix.LoopToBounds(px,py,pz);
      // printf("%s loop to bounds  \n",fName.Data() );
    }

  } else {

    // Neutral particle or no field

    MCLine line(fRnrStyle, &mc_v0, TMath::C()*fBeta, &track_points);
   
    if(!fPathMarks.empty())
    {
      for(std::vector<Reve::PathMark*>::iterator i=fPathMarks.begin(); i!=fPathMarks.end(); ++i)
      {
	Reve::PathMark* pm = *i;

	if(RS.fFitDaughters &&  pm->type == Reve::PathMark::Daughter)
        {
          if(TMath::Abs(pm->V.z) > RS.fMaxZ ||
	     TMath::Sqrt(pm->V.x*pm->V.x + pm->V.y*pm->V.y) > RS.fMaxR)
          {
	    goto line_bounds;
          }
	  line.GotoVertex(pm->V.x, pm->V.y, pm->V.z);
	  fP.x -= pm->P.x;
	  fP.y -= pm->P.y;
	  fP.z -= pm->P.z;
	}

	if(RS.fFitDecay &&  pm->type == Reve::PathMark::Decay)
        {
	  if(TMath::Abs(pm->V.z) > RS.fMaxZ ||
	     TMath::Sqrt(pm->V.x*pm->V.x + pm->V.y*pm->V.y) > RS.fMaxR)
          {
	    goto line_bounds;
          }
	  line.GotoVertex(pm->V.x, pm->V.y, pm->V.z);
          decay = true;
	  break;
	}
      }
    }

  line_bounds:
    if(!decay || RS.fFitDecay == kFALSE)
      line.GotoBounds(px,py,pz);

  }
make_polyline:
  {
    Int_t size = TMath::Min(4096, (Int_t) track_points.size());
    // printf("track '%s'   N = %u\n", GetName(), track_points.size());
    Reset(size);
    for(Int_t i=0; i<size; ++i)
    {
      const MCVertex& v = track_points[i];
      SetNextPoint(v.x, v.y, v.z);
    }
  }

  if(recurse)
  {
    for(List_i i=fChildren.begin(); i!=fChildren.end(); ++i)
    {
      Track* t = dynamic_cast<Track*>(*i);
      if(t) t->MakeTrack(recurse); 
    }
  }
}

/******************************************************************************/

//______________________________________________________________________________
TClass* Track::ProjectedClass() const
{
  return NLTTrack::Class();
}

/******************************************************************************/

namespace {

struct cmp_pathmark
{
  bool operator()(PathMark* const & a, PathMark* const & b)
  { return a->time < b->time; }
};

}

//______________________________________________________________________________
void Track::SortPathMarksByTime()
{
  std::sort(fPathMarks.begin(), fPathMarks.end(), cmp_pathmark());
}

/******************************************************************************/

//______________________________________________________________________________
void Track::ImportHits()
{
  Reve::LoadMacro("hits_from_label.C");
  gROOT->ProcessLine(Form("hits_from_label(%d, (Reve::RenderElement*)%p);", 
			  fLabel, this));
}

//______________________________________________________________________________
void Track::ImportClusters()
{
  Reve::LoadMacro("clusters_from_label.C");
  gROOT->ProcessLine(Form("clusters_from_label(%d, (Reve::RenderElement*)%p);", 
			  fLabel, this));
}

//______________________________________________________________________________
void Track::ImportClustersFromIndex()
{
  static const Exc_t eH("Track::ImportClustersFromIndex ");

  if (fIndex == kMinInt)
    throw(eH + "index not set.");

  Reve::LoadMacro("clusters_from_index.C");
  gROOT->ProcessLine(Form("clusters_from_index(%d, (Reve::RenderElement*)%p);", 
			  fIndex, this));
}

/******************************************************************************/

//______________________________________________________________________________
void Track::ImportKine()
{
  static const Exc_t eH("Track::ImportKine ");

  if (fLabel == kMinInt)
    throw(eH + "label not set.");

  Int_t label;
  if (fLabel < 0) {
    Warning(eH, "label negative, taking absolute value.");
    label = -fLabel;
  } else {
    label = fLabel;
  }

  Reve::LoadMacro("kine_tracks.C");
  gROOT->ProcessLine(Form("kine_track(%d, kFALSE, kTRUE, (Reve::RenderElement*)%p);", 
			  label, this));

}

//______________________________________________________________________________
void Track::ImportKineWithArgs(Bool_t importMother, Bool_t importDaugters)
{
  static const Exc_t eH("Track::ImportKineWithArgs ");

  if (fLabel == kMinInt)
    throw(eH + "label not set.");

  Int_t label;
  if (fLabel < 0) {
    Warning(eH, "label negative, taking absolute value.");
    label = -fLabel;
  } else {
    label = fLabel;
  }

  Reve::LoadMacro("kine_tracks.C");
  gROOT->ProcessLine(Form("kine_track(%d, %d, %d, (Reve::RenderElement*)%p);", 
			   label, importMother, importDaugters, this));
}

/******************************************************************************/

//______________________________________________________________________________
void Track::PrintKineStack()
{
  Reve::LoadMacro("print_kine_from_label.C");
  gROOT->ProcessLine(Form("print_kine_from_label(%d);", fLabel));
}

//______________________________________________________________________________
void Track::PrintPathMarks()
{
  static const Exc_t eH("Track::PrintPathMarks ");

  printf("Track '%s', number of path marks %d, label %d\n",
	 GetName(), fPathMarks.size(), fLabel);

  PathMark* pm;
  for(vpPathMark_i i=fPathMarks.begin(); i!=fPathMarks.end(); i++) 
  {
    pm = *i;
    printf("  %-9s  p: %8f %8f %8f Vertex: %8e %8e %8e %g \n",
	   pm->type_name(),
	   pm->P.x,  pm->P.y, pm->P.z,
	   pm->V.x,  pm->V.y, pm->V.z,
	   pm->time);
  }
}

/******************************************************************************/

//______________________________________________________________________________
void Track::CtrlClicked(Reve::Track* track)
{
  Emit("CtrlClicked(Reve::Track*)", (Long_t)track);
}

//______________________________________________________________________________
void Track::SetLineStyle(Style_t lstyle)
{
  TAttLine::SetLineStyle(lstyle);
  std::list<NLTProjected*>::iterator pi = fProjectedList.begin();
  while (pi != fProjectedList.end())
  {
    Track* pt = dynamic_cast<Track*>(*pi);
    if (pt)
    {
      pt->SetLineStyle(lstyle);
      pt->ElementChanged();
    }
    ++pi;
  }
}


/******************************************************************************/
/******************************************************************************/

//______________________________________________________________________________
// TrackRnrStyle
//
// Holding structure for a number of track rendering parameters.
//
// This is decoupled from Track/TrackList to allow sharing of the
// RnrStyle among several instances. Back references are kept so the
// tracks can be recreated when the parameters change.
//
// TrackList has Get/Set methods for RnrStlye. TrackEditor and
// TrackListEditor provide editor access.

ClassImp(Reve::TrackRnrStyle)

Float_t       TrackRnrStyle::fgDefMagField = 5;
const Float_t TrackRnrStyle::fgkB2C        = 0.299792458e-3;
TrackRnrStyle TrackRnrStyle::fgDefStyle;

//______________________________________________________________________________
TrackRnrStyle::TrackRnrStyle() :
  TObject(),
  ReferenceBackPtr(),

  fMagField(fgDefMagField),

  fMaxR  (350),
  fMaxZ  (450),

  fMaxOrbs (0.5),
  fMinAng  (45),
  fDelta   (0.1),

  fEditPathMarks(kFALSE),
  fPMAtt(),

  fFitDaughters  (kTRUE),
  fFitReferences (kTRUE),
  fFitDecay      (kTRUE),

  fRnrDaughters  (kTRUE),
  fRnrReferences (kTRUE),
  fRnrDecay      (kTRUE),

  fRnrFV(kFALSE),
  fFVAtt()
{
  // Default constructor.

  fPMAtt.SetMarkerColor(4);
  fPMAtt.SetMarkerStyle(2);

  fFVAtt.SetMarkerSize(0.6);
  fFVAtt.SetMarkerColor(4);
  fFVAtt.SetMarkerStyle(2);
}

/**************************************************************************/

//______________________________________________________________________________
void TrackRnrStyle::RebuildTracks()
{
  // Rebuild all tracks using this render-style.

  Track* track;
  std::list<RenderElement*>::iterator i = fBackRefs.begin();  
  while (i != fBackRefs.end())
  {
    track = dynamic_cast<Track*>(*i);
    track->MakeTrack();
    ++i;
  }
}

/******************************************************************************/

//______________________________________________________________________________
void TrackRnrStyle::SetMaxR(Float_t x)
{
  // Set maximum radius and rebuild tracks.

  fMaxR = x;
  RebuildTracks();
}

//______________________________________________________________________________
void TrackRnrStyle::SetMaxZ(Float_t x)
{
  // Set maximum z and rebuild tracks.

  fMaxZ = x;
  RebuildTracks();
}

//______________________________________________________________________________
void TrackRnrStyle::SetMaxOrbs(Float_t x)
{
  // Set maximum number of orbits and rebuild tracks.

  fMaxOrbs = x;
  RebuildTracks();
}

//______________________________________________________________________________
void TrackRnrStyle::SetMinAng(Float_t x)
{
  // Set minimum step angle and rebuild tracks.

  fMinAng = x;
  RebuildTracks();
}

//______________________________________________________________________________
void TrackRnrStyle::SetDelta(Float_t x)
{
  // Set maximum error and rebuild tracks.

  fDelta = x;
  RebuildTracks();
}

//______________________________________________________________________________
void TrackRnrStyle::SetFitDaughters(Bool_t x)
{
  // Set daughter creation point fitting and rebuild tracks.

  fFitDaughters = x;
  RebuildTracks();
}

//______________________________________________________________________________
void TrackRnrStyle::SetFitReferences(Bool_t x)
{
  // Set track-reference fitting and rebuild tracks.

  fFitReferences = x;
  RebuildTracks();
}

//______________________________________________________________________________
void TrackRnrStyle::SetFitDecay(Bool_t x)
{
  // Set decay fitting and rebuild tracks.

  fFitDecay = x;
  RebuildTracks();
}

//______________________________________________________________________________
void TrackRnrStyle::SetRnrDecay(Bool_t rnr)
{
  // Set decay rendering and rebuild tracks.

  fRnrDecay = rnr;
  RebuildTracks();
}

//______________________________________________________________________________
void TrackRnrStyle::SetRnrDaughters(Bool_t rnr)
{
  // Set daughter rendering and rebuild tracks.

  fRnrDaughters = rnr;
  RebuildTracks();
}

//______________________________________________________________________________
void TrackRnrStyle::SetRnrReferences(Bool_t rnr)
{
  // Set track-reference rendering and rebuild tracks.

  fRnrReferences = rnr;
  RebuildTracks();
}
 

/**************************************************************************/
/**************************************************************************/

//______________________________________________________________________________
// TrackList
//

ClassImp(Reve::TrackList)

//______________________________________________________________________________
TrackList::TrackList(TrackRnrStyle* rs) :
  RenderElementList(),
  TAttMarker(1, 20, 1),
  TAttLine(1,1,1),

  fRecurse(kTRUE),
  fRnrStyle(0),
  fRnrLine(kTRUE),
  fRnrPoints(kFALSE),

  fMinPt (0), fMaxPt (0), fLimPt (0),
  fMinP  (0), fMaxP  (0), fLimP  (0)
{
  // Constructor. If TrackRenderStyle argument is 0, a new default
  // render-style is created.

  fChildClass = Track::Class(); // override member from base RenderElementList

  fMainColorPtr = &fLineColor;
  if (fRnrStyle== 0) rs = new TrackRnrStyle;
  SetRnrStyle(rs);
}

//______________________________________________________________________________
TrackList::TrackList(const Text_t* name, TrackRnrStyle* rs) :
  RenderElementList(name),
  TAttMarker(1, 20, 1),
  TAttLine(1,1,1),

  fRecurse(kTRUE),
  fRnrStyle      (0),
  fRnrLine(kTRUE),
  fRnrPoints(kFALSE),

  fMinPt (0), fMaxPt (0), fLimPt (0),
  fMinP  (0), fMaxP  (0), fLimP  (0)
{
  // Constructor. If TrackRenderStyle argument is 0, a new default
  // render-style is created.

  fChildClass = Track::Class(); // override member from base RenderElementList

  fMainColorPtr = &fLineColor;
  if (fRnrStyle== 0) rs = new TrackRnrStyle;
  SetRnrStyle(rs);
}

//______________________________________________________________________________
TrackList::~TrackList()
{
  // Destructor.

  SetRnrStyle(0);
}

/******************************************************************************/

//______________________________________________________________________________
void TrackList::SetRnrStyle(TrackRnrStyle* rs)
{
  // Set default render-style for tracks.
  // This is not enforced onto the tracks themselves but this is the
  // render-style that is show in the TrackListEditor.

  if (fRnrStyle == rs) return;
  if (fRnrStyle) fRnrStyle->DecRefCount();
  fRnrStyle = rs;
  if (fRnrStyle) rs->IncRefCount();
}

/**************************************************************************/

//______________________________________________________________________________
void TrackList::MakeTracks(Bool_t recurse)
{
  // Regenerate the visual representations of tracks.
  // The momentum limits are rescanned during the same traversal.

  fLimPt = fLimP = 0;

  for (List_i i=fChildren.begin(); i!=fChildren.end(); ++i)
  {
    Track* track = (Track*)(*i);
    track->MakeTrack(recurse);

    fLimPt = TMath::Max(fLimPt, track->fP.Perp());
    fLimP  = TMath::Max(fLimP,  track->fP.Mag());
    if (recurse)
      FindMomentumLimits(*i, recurse);
  }

  fLimPt = RoundMomentumLimit(fLimPt);
  fLimP  = RoundMomentumLimit(fLimP);
  if (fMaxPt == 0) fMaxPt = fLimPt;
  if (fMaxP  == 0) fMaxP  = fLimP;

  gReve->Redraw3D();
}

//______________________________________________________________________________
void TrackList::FindMomentumLimits(RenderElement* el, Bool_t recurse)
{
  // Loop over track elements of argument el and find highest pT and p.
  // These are stored in members fLimPt and fLimP.

  for (List_i i=el->BeginChildren(); i!=el->EndChildren(); ++i)
  {
    Track* track = dynamic_cast<Track*>(*i);
    if (track)
    {
      fLimPt = TMath::Max(fLimPt, track->fP.Perp());
      fLimP  = TMath::Max(fLimP,  track->fP.Mag());
      if (recurse)
        FindMomentumLimits(*i, recurse);
    }
  }
}

//______________________________________________________________________________
Float_t TrackList::RoundMomentumLimit(Float_t x)
{
  // Round the momentum limit up to a nice value.

  using namespace TMath;
  Double_t fac = Power(10, 1 - Floor(Log10(x)));
  return Ceil(fac*x) / fac;
}

/**************************************************************************/

//______________________________________________________________________________
void TrackList::SetRnrLine(Bool_t rnr)
{
  // Set rendering of track as line for the list and the elements.

  for (List_i i=BeginChildren(); i!=EndChildren(); ++i)
  {
    Track* track = (Track*)(*i);
    if (track->GetRnrLine() == fRnrLine)
      track->SetRnrLine(rnr);
    if (fRecurse)
      SetRnrLine(rnr, *i);
  }
  fRnrLine = rnr;
}

//______________________________________________________________________________
void TrackList::SetRnrLine(Bool_t rnr, RenderElement* el)
{
  // Set rendering of track as line for children of el.

  Track* track;
  for (List_i i=el->BeginChildren(); i!=el->EndChildren(); ++i)
  {
    track = dynamic_cast<Track*>(*i);
    if (track && (track->GetRnrLine() == fRnrLine))
      track->SetRnrLine(rnr);
    if (fRecurse)
      SetRnrLine(rnr, *i);
  }
}

/******************************************************************************/

//______________________________________________________________________________
void TrackList::SetRnrPoints(Bool_t rnr)
{
  // Set rendering of track as points for the list and the elements.

  for (List_i i=BeginChildren(); i!=EndChildren(); ++i)
  {
    Track* track = (Track*)(*i);
    if (track->GetRnrPoints() == fRnrPoints)
      track->SetRnrPoints(rnr);
    if (fRecurse)
      SetRnrPoints(rnr, *i);
  }
  fRnrPoints = rnr;
}

//______________________________________________________________________________
void TrackList::SetRnrPoints(Bool_t rnr, RenderElement* el)
{
  // Set rendering of track as points for children of el.

  Track* track;
  for (List_i i=el->BeginChildren(); i!=el->EndChildren(); ++i)
  {
    track = dynamic_cast<Track*>(*i);
    if (track) 
      if (track->GetRnrPoints() == fRnrPoints)
	track->SetRnrPoints(rnr);
    if (fRecurse)
      SetRnrPoints(rnr, *i);
  }
}

/******************************************************************************/

//______________________________________________________________________________
void TrackList::SetMainColor(Color_t col)
{
  // Set main (line) color for the list and the elements.

  for (List_i i=BeginChildren(); i!=EndChildren(); ++i)
  {
    Track* track = (Track*)(*i);
    if (track->GetLineColor() == fLineColor)
      track->SetLineColor(col);
    if (fRecurse)
      SetLineColor(col, *i);
  }
  RenderElement::SetMainColor(col);
}

//______________________________________________________________________________
void TrackList::SetLineColor(Color_t col, RenderElement* el)
{
  // Set line color for children of el.

  Track* track;
  for (List_i i=el->BeginChildren(); i!=el->EndChildren(); ++i)
  {
    track = dynamic_cast<Track*>(*i);
    if (track && track->GetLineColor() == fLineColor)
      track->SetLineColor(col);
    if (fRecurse)
      SetLineColor(col, *i);
  }
}

/******************************************************************************/

//______________________________________________________________________________
void TrackList::SetLineWidth(Width_t width)
{
  // Set line width for the list and the elements.

  for (List_i i=BeginChildren(); i!=EndChildren(); ++i)
  {
    Track* track = (Track*)(*i);
    if (track->GetLineWidth() == fLineWidth)
      track->SetLineWidth(width);
    if (fRecurse)
      SetLineWidth(width, *i);
  }
  fLineWidth=width; 
}

//______________________________________________________________________________
void TrackList::SetLineWidth(Width_t width, RenderElement* el)
{
  // Set line width for children of el.

  Track* track;
  for (List_i i=el->BeginChildren(); i!=el->EndChildren(); ++i)
  {
    track = dynamic_cast<Track*>(*i);
    if (track && track->GetLineWidth() == fLineWidth)
      track->SetLineWidth(width);
    if (fRecurse)
      SetLineWidth(width, *i);
  }
}

/******************************************************************************/

//______________________________________________________________________________
void TrackList::SetLineStyle(Style_t style)
{
  // Set line style for the list and the elements.

  for (List_i i=BeginChildren(); i!=EndChildren(); ++i)
  {
    Track* track = (Track*)(*i);
    if (track->GetLineStyle() == fLineStyle)
      track->SetLineStyle(style);
    if (fRecurse)
      SetLineStyle(style, *i);
  }
  fLineStyle=style; 
}

//______________________________________________________________________________
void TrackList::SetLineStyle(Style_t style, RenderElement* el)
{
  // Set line style for children of el.

  Track* track;
  for (List_i i=el->BeginChildren(); i!=el->EndChildren(); ++i)
  {
    track = dynamic_cast<Track*>(*i);
    if (track && track->GetLineStyle() == fLineStyle)
      track->SetLineStyle(style);
    if (fRecurse)
      SetLineStyle(style, *i);
  }
}

/******************************************************************************/

//______________________________________________________________________________
void TrackList::SetMarkerStyle(Style_t style)
{
  // Set marker style for the list and the elements.

  for (List_i i=BeginChildren(); i!=EndChildren(); ++i)
  {
    Track* track = (Track*)(*i);
    if (track->GetMarkerStyle() == fMarkerStyle)
      track->SetMarkerStyle(style);
    if (fRecurse)
      SetMarkerStyle(style, *i);
  }
  fMarkerStyle=style; 
}

//______________________________________________________________________________
void TrackList::SetMarkerStyle(Style_t style, RenderElement* el)
{
  // Set marker style for children of el.

  Track* track;
  for (List_i i=el->BeginChildren(); i!=el->EndChildren(); ++i)
  {
    track = dynamic_cast<Track*>(*i);
    if (track && track->GetMarkerStyle() == fMarkerStyle)
      track->SetMarkerStyle(style);
    if(fRecurse)
      SetMarkerStyle(style, *i);
  }
}

/******************************************************************************/

//______________________________________________________________________________
void TrackList::SetMarkerColor(Color_t col)
{
  // Set marker color for the list and the elements.

  for (List_i i=BeginChildren(); i!=EndChildren(); ++i)
  {
    Track* track = (Track*)(*i);
    if (track->GetMarkerColor() == fMarkerColor)
      track->SetMarkerColor(col);
    if (fRecurse)
      SetMarkerColor(col, *i);
  }
  fMarkerColor=col; 
}

//______________________________________________________________________________
void TrackList::SetMarkerColor(Color_t col, RenderElement* el)
{
  // Set marker color for children of el.

  Track* track;
  for (List_i i=el->BeginChildren(); i!=el->EndChildren(); ++i)
  {
    track = dynamic_cast<Track*>(*i);
    if (track && track->GetMarkerColor() == fMarkerColor)
      track->SetMarkerColor(col);
    if (fRecurse)
      SetMarkerColor(col, *i);
  }
}

/******************************************************************************/

//______________________________________________________________________________
void TrackList::SetMarkerSize(Size_t size)
{
  // Set marker size for the list and the elements.

  for (List_i i=BeginChildren(); i!=EndChildren(); ++i)
  {
    Track* track = (Track*)(*i);
    if (track->GetMarkerSize() == fMarkerSize)
      track->SetMarkerSize(size);
    if (fRecurse)
      SetMarkerSize(size, *i);
  }
  fMarkerSize=size; 
}

//______________________________________________________________________________
void TrackList::SetMarkerSize(Size_t size, RenderElement* el)
{
  // Set marker size for children of el.

  Track* track;
  for (List_i i=el->BeginChildren(); i!=el->EndChildren(); ++i)
  {
    track = dynamic_cast<Track*>(*i);
    if (track && track->GetMarkerSize() == fMarkerSize)
      track->SetMarkerSize(size);
    if (fRecurse)
      SetMarkerSize(size, *i);
  }
}

/******************************************************************************/

//______________________________________________________________________________
void TrackList::SelectByPt(Float_t min_pt, Float_t max_pt)
{
  fMinPt = min_pt;
  fMaxPt = max_pt;

  const Float_t minptsq = min_pt*min_pt;
  const Float_t maxptsq = max_pt*max_pt;

  for (List_i i=BeginChildren(); i!=EndChildren(); ++i)
  {
    const Float_t ptsq = ((Track*)(*i))->fP.Perp2();
    Bool_t on = ptsq >= minptsq && ptsq <= maxptsq;
    (*i)->SetRnrState(on);
    if (on && fRecurse)
      SelectByPt(min_pt, max_pt, *i);
  }
}

//______________________________________________________________________________
void TrackList::SelectByPt(Float_t min_pt, Float_t max_pt, RenderElement* el)
{
  const Float_t minptsq = min_pt*min_pt;
  const Float_t maxptsq = max_pt*max_pt;

  for (List_i i=el->BeginChildren(); i!=el->EndChildren(); ++i)
  {
    Track* track = dynamic_cast<Track*>(*i);
    if (track)
    {
      const Float_t ptsq = track->fP.Perp2();
      Bool_t on = ptsq >= minptsq && ptsq <= maxptsq;
      track->SetRnrState(on);
      if (on && fRecurse)
        SelectByPt(min_pt, max_pt, *i);
    }
  }
}

//______________________________________________________________________________
void TrackList::SelectByP(Float_t min_p, Float_t max_p)
{
  fMinP = min_p;
  fMaxP = max_p;

  const Float_t minpsq = min_p*min_p;
  const Float_t maxpsq = max_p*max_p;

  for (List_i i=BeginChildren(); i!=EndChildren(); ++i)
  {
    const Float_t psq  = ((Track*)(*i))->fP.Mag2();
    Bool_t on = psq >= minpsq && psq <= maxpsq;
    (*i)->SetRnrState(psq >= minpsq && psq <= maxpsq);
    if (on && fRecurse)
      SelectByP(min_p, max_p, *i);
  }
}

//______________________________________________________________________________
void TrackList::SelectByP(Float_t min_p, Float_t max_p, RenderElement* el)
{
  const Float_t minpsq = min_p*min_p;
  const Float_t maxpsq = max_p*max_p;

  for (List_i i=el->BeginChildren(); i!=el->EndChildren(); ++i)
  {
    Track* track = dynamic_cast<Track*>(*i);
    if (track)
    {
      const Float_t psq  = ((Track*)(*i))->fP.Mag2();
      Bool_t on = psq >= minpsq && psq <= maxpsq;
      track->SetRnrState(on);
      if (on && fRecurse)
        SelectByP(min_p, max_p, *i);
    }
  }
}

/******************************************************************************/

//______________________________________________________________________________
Track* TrackList::FindTrackByLabel(Int_t label)
{
  // Find track by label, select it and display it in the editor.

  for (List_i i=fChildren.begin(); i!=fChildren.end(); ++i) {
    if (((Track*)(*i))->GetLabel() == label) {
      TGListTree     *lt   = gReve->GetLTEFrame()->GetListTree();
      TGListTreeItem *mlti = lt->GetSelected();
      if (mlti->GetUserData() != this)
	mlti = FindListTreeItem(lt);
      TGListTreeItem *tlti = (*i)->FindListTreeItem(lt, mlti);
      lt->HighlightItem(tlti);
      lt->SetSelected(tlti);
      gReve->EditRenderElement(*i);
      return (Track*) *i;
    }
  }
  return 0;
}

//______________________________________________________________________________
Track* TrackList::FindTrackByIndex(Int_t index)
{
  // Find track by index, select it and display it in the editor.

  for (List_i i=fChildren.begin(); i!=fChildren.end(); ++i) {
    if (((Track*)(*i))->GetIndex() == index) {
      TGListTree     *lt   = gReve->GetLTEFrame()->GetListTree();
      TGListTreeItem *mlti = lt->GetSelected();
      if (mlti->GetUserData() != this)
	mlti = FindListTreeItem(lt);
      TGListTreeItem *tlti = (*i)->FindListTreeItem(lt, mlti);
      lt->HighlightItem(tlti);
      lt->SetSelected(tlti);
      gReve->EditRenderElement(*i);
      return (Track*) *i;
    }
  }
  return 0;
}

//______________________________________________________________________________
void TrackList::ImportHits()
{
  for (List_i i=fChildren.begin(); i!=fChildren.end(); ++i) {
    ((Track*)(*i))->ImportHits();
  }
}

//______________________________________________________________________________
void TrackList::ImportClusters()
{
  for (List_i i=fChildren.begin(); i!=fChildren.end(); ++i) {
    ((Track*)(*i))->ImportClusters();
  }
}

/******************************************************************************/

//______________________________________________________________________________
TClass* TrackList::ProjectedClass() const
{
  return NLTTrackList::Class();
}


/******************************************************************************/
/******************************************************************************/

#include "RGEditor.h"

//______________________________________________________________________________
// TrackCounter
//

ClassImp(TrackCounter)

TrackCounter* TrackCounter::fgInstance = 0;

TrackCounter::TrackCounter(const Text_t* name, const Text_t* title) :
  RenderElement(),
  TNamed(name, title),

  fBadLineStyle (6),
  fClickAction  (CA_ToggleTrack),
  fAllTracks    (0),
  fGoodTracks   (0),
  fTrackLists   ()
{
  if (fgInstance == 0) fgInstance = this;
  TQObject::Connect("Reve::Track", "CtrlClicked(Reve::Track*)",
		    "Reve::TrackCounter", this, "DoTrackAction(Reve::Track*)");
}

TrackCounter::~TrackCounter()
{
  TQObject::Disconnect("Reve::Track", "DoTrackAction(Reve::Track*)");
  if (fgInstance == this) fgInstance = 0;
}

/**************************************************************************/

void TrackCounter::Reset()
{
  printf("TrackCounter::Reset()\n");
  fAllTracks  = 0;
  fGoodTracks = 0;
  TIter next(&fTrackLists);
  TrackList* tlist;
  while ((tlist = dynamic_cast<TrackList*>(next())))
    tlist->DecDenyDestroy();
  fTrackLists.Clear("nodelete");
}

void TrackCounter::RegisterTracks(TrackList* tlist, Bool_t goodTracks)
{
  // printf("TrackCounter::RegisterTracks '%s', %s\n",
  //   tlist->GetObject()->GetName(), goodTracks ? "good" : "bad");

  tlist->IncDenyDestroy();
  fTrackLists.Add(tlist);

  List_i i = tlist->BeginChildren();
  while (i != tlist->EndChildren())
  {
    Track* t = dynamic_cast<Track*>(*i);
    if (t != 0)
    {
      if (goodTracks)
      {
	++fGoodTracks;
      } else {
	t->SetLineStyle(fBadLineStyle);
      }
      ++fAllTracks;
    }
    ++i;
  }
}

void TrackCounter::DoTrackAction(Track* track)
{
  // !!!! No check done if ok.
  // !!!! Should also override RemoveElementLocal
  // !!!! But then ... should also store local information if track is ok.

  switch (fClickAction)
  {

    case CA_PrintTrackInfo:
    {
      printf("Track '%s'\n", track->GetObject()->GetName());
      Vector &v = track->fV, &p = track->fP;
      printf("  Vx=%f, Vy=%f, Vz=%f; Pt=%f, Pz=%f, phi=%f)\n",
	     v.x, v.y, v.z, p.Perp(), p.z, TMath::RadToDeg()*p.Phi());
      printf("  <other information should be printed ... full AliESDtrack>\n");
      break;
    }

    case CA_ToggleTrack:
    {
      if (track->GetLineStyle() == 1)
      {
	track->SetLineStyle(fBadLineStyle);
	--fGoodTracks;
      } else {
	track->SetLineStyle(1);
	++fGoodTracks;
      }
      track->ElementChanged();
      gReve->Redraw3D();

      printf("TrackCounter::CountTrack All=%d, Good=%d, Bad=%d\n",
	     fAllTracks, fGoodTracks, fAllTracks-fGoodTracks);

      if (gReve->GetEditor()->GetModel() == GetObject())
	gReve->EditRenderElement(this);

      break;
    }

  } // end switch fClickAction
}

/**************************************************************************/

void TrackCounter::OutputEventTracks(FILE* out)
{
  if (out == 0)
  {
    out = stdout;
    fprintf(out, "TrackCounter::FinalizeEvent()\n");
  }

  fprintf(out, "Event = %d  Ntracks = %d\n", fEventId, fGoodTracks);

  TIter tlists(&fTrackLists);
  TrackList* tlist;
  Int_t cnt = 0;
  while ((tlist = (TrackList*) tlists()) != 0)
  {
    List_i i = tlist->BeginChildren();
    while (i != tlist->EndChildren())
    {
      Track* t = dynamic_cast<Track*>(*i);
      if (t != 0 && t->GetLineStyle() == 1)
      {
	++cnt;
	fprintf(out, " %2d: chg=%+2d  pt=%8.5f  eta=%+8.5f\n",
	       cnt, t->fCharge, t->fP.Perp(), t->fP.Eta());
      }
      ++i;
    }
  }
}
