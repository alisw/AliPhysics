// $Header$

#include "Track.h"
#include "MCHelixLine.hi"

#include <TPolyLine3D.h>
#include <TPolyMarker3D.h>
#include <TColor.h>

// Updates
#include <Reve/RGTopFrame.h>
#include <TCanvas.h>

#include <vector>

using namespace Reve;

//______________________________________________________________________
// Track
//

ClassImp(Reve::Track)

Track::Track() :
  Line(),

  fV(),
  fP(),
  fBeta(0),
  fCharge(0),
  fLabel(0),
  fPathMarks(),

  fRnrStyle(0)
{}

Track::Track(TParticle* t, Int_t label, TrackRnrStyle* rs):
  Line(),

  fV(t->Vx(), t->Vy(), t->Vz()),
  fP(t->Px(), t->Py(), t->Pz()),
  fBeta(t->P()/t->Energy()),
  fCharge(0),
  fLabel(label),
  fPathMarks(),

  fRnrStyle(rs)
{
  fLineColor = fRnrStyle->GetColor();
  fMainColorPtr = &fLineColor;

  TParticlePDG* pdgp = t->GetPDG();
  if (pdgp)
    fCharge = (Int_t) TMath::Nint(pdgp->Charge()/3);

  SetName(t->GetName());
}

Track::Track(Reve::MCTrack* t, TrackRnrStyle* rs):
  Line(),

  fV(t->Vx(), t->Vy(), t->Vz()),
  fP(t->Px(), t->Py(), t->Pz()),
  fBeta(t->P()/t->Energy()),
  fCharge(0),
  fLabel(t->label),
  fIndex(t->index),
  fPathMarks(),

  fRnrStyle(rs)
{
  fLineColor = fRnrStyle->GetColor();
  fMainColorPtr = &fLineColor;

  TParticlePDG* pdgp = t->GetPDG();
  if(pdgp == 0) {
    t->ResetPdgCode(); pdgp = t->GetPDG();
  }
  fCharge = (Int_t) TMath::Nint(pdgp->Charge()/3);

  SetName(t->GetName());
}

Track::Track(Reve::RecTrack* t, TrackRnrStyle* rs) :
  Line(),

  fV(t->V),
  fP(t->P),
  fBeta(t->beta),
  fCharge(t->sign),
  fLabel(t->label),
  fIndex(t->index),
  fPathMarks(),

  fRnrStyle(rs)
{
  fLineColor = fRnrStyle->GetColor();
  fMainColorPtr = &fLineColor;

  SetName(t->GetName());
}

Track::~Track()
{
  for (vpPathMark_i i=fPathMarks.begin(); i!=fPathMarks.end(); ++i)
    delete *i;
}

/*
void Track::Reset(Int_t n_points)
{
  delete [] TPolyLine3D::fP; TPolyLine3D::fP = 0;
  fN = n_points;
  if(fN) TPolyLine3D::fP = new Float_t [3*fN];
  memset(TPolyLine3D::fP, 0, 3*fN*sizeof(Float_t));
  fLastPoint = -1;
}
*/

/**************************************************************************/

void Track::MakeTrack()
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
  
  if (fCharge != 0 && TMath::Abs(RS.fMagField) > 1e-5) {

    // Charged particle in magnetic field

    Float_t a = RS.fgkB2C * RS.fMagField * fCharge;
   
    MCHelix helix(fRnrStyle, &mc_v0, TMath::C()*fBeta, &track_points, a); //m->cm
    helix.Init(TMath::Sqrt(px*px+py*py), pz);
   
    if(!fPathMarks.empty()){
      for(std::vector<Reve::PathMark*>::iterator i=fPathMarks.begin(); i!=fPathMarks.end(); ++i) {
	Reve::PathMark* pm = *i;
        
	if(RS.fFitDaughters &&  pm->type == Reve::PathMark::Daughter){
	  if(TMath::Abs(pm->V.z) > RS.fMaxZ 
	     || TMath::Sqrt(pm->V.x*pm->V.x + pm->V.y*pm->V.y) > RS.fMaxR )
	    goto helix_bounds;

          //printf("%s fit daughter  \n", fName.Data()); 
	  helix.LoopToVertex(fP.x, fP.y, fP.z, pm->V.x, pm->V.y, pm->V.z);
	  fP.x -=  pm->P.x;
	  fP.y -=  pm->P.y;
	  fP.z -=  pm->P.z;
	}
	if(RS.fFitDecay &&  pm->type == Reve::PathMark::Decay){
	  
	  if(TMath::Abs(pm->V.z) > RS.fMaxZ 
	     || TMath::Sqrt(pm->V.x*pm->V.x + pm->V.y*pm->V.y) > RS.fMaxR )
	    goto helix_bounds;
	  helix.LoopToVertex(fP.x, fP.y, fP.z, pm->V.x, pm->V.y, pm->V.z);
          decay = true;
          break;
	}
      }
    }
  helix_bounds:
    //go to bounds
    if(!decay || RS.fFitDecay == kFALSE){
      helix.LoopToBounds(px,py,pz);
      // printf("%s loop to bounds  \n",fName.Data() );
    }

  } else {

    // Neutral particle or no field

    MCLine line(fRnrStyle, &mc_v0, TMath::C()*fBeta, &track_points);
   
    if(!fPathMarks.empty()){
      for(std::vector<Reve::PathMark*>::iterator i=fPathMarks.begin(); i!=fPathMarks.end(); ++i) {
	Reve::PathMark* pm = *i;

	if(RS.fFitDaughters &&  pm->type == Reve::PathMark::Daughter){
          if(TMath::Abs(pm->V.z) > RS.fMaxZ 
	     || TMath::Sqrt(pm->V.x*pm->V.x + pm->V.y*pm->V.y) > RS.fMaxR )
	    goto line_bounds;
	  line.GotoVertex(pm->V.x, pm->V.y, pm->V.z);
	  fP.x -=  pm->P.x;
	  fP.y -=  pm->P.y;
	  fP.z -=  pm->P.z;
	}

	if(RS.fFitDecay &&  pm->type == Reve::PathMark::Decay){
	  if(TMath::Abs(pm->V.z) > RS.fMaxZ 
	     || TMath::Sqrt(pm->V.x*pm->V.x + pm->V.y*pm->V.y) > RS.fMaxR )
	    goto line_bounds;
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
  Reset(track_points.size());
  for(std::vector<MCVertex>::iterator i=track_points.begin(); i!=track_points.end(); ++i)
    SetNextPoint(i->x, i->y, i->z);
}

/**************************************************************************/

void Track::ImportHits()
{
  Reve::LoadMacro("hits_from_label.C");
  gROOT->ProcessLine(Form("hits_from_label(%d);", fLabel));
}

void Track::ImportClusters()
{
  Reve::LoadMacro("clusters_from_label.C");
  gROOT->ProcessLine(Form("clusters_from_label(%d);", fLabel));
}

void Track::ImportClustersFromIndex()
{
  Reve::LoadMacro("clusters_from_index.C");
  gROOT->ProcessLine(Form("clusters_from_index(%d);", fIndex));
}

void Track::PrintKineStack()
{
  Reve::LoadMacro("print_kine_from_label.C");
  gROOT->ProcessLine(Form("print_kine_from_label(%d);", fLabel));
}

/**************************************************************************/

void Track::CtrlClicked(Reve::Track* track)
{
  Emit("CtrlClicked(Reve::Track*)", (Long_t)track);
}


/**************************************************************************/
/**************************************************************************/

//______________________________________________________________________
// TrackRnrStyle
//

ClassImp(Reve::TrackRnrStyle)

Float_t       TrackRnrStyle::fgDefMagField = 5;
const Float_t TrackRnrStyle::fgkB2C        = 0.299792458e-3;
TrackRnrStyle TrackRnrStyle::fgDefStyle;

TrackRnrStyle::TrackRnrStyle() :
  TObject(),

  fColor(1),
  fWidth(1),
  fMagField(fgDefMagField),

  fMaxR  (350),
  fMaxZ  (450),

  fMaxOrbs (0.5),
  fMinAng  (45),
  fDelta   (0.1),

  fFitDaughters(kTRUE),
  fFitDecay    (kTRUE)
{}

/**************************************************************************/
/**************************************************************************/

//______________________________________________________________________
// TrackList
//

ClassImp(Reve::TrackList)

void TrackList::Init()
{
  fMarkerStyle = 6;
  fMarkerColor = 5;
  // fMarker->SetMarkerSize(0.05);

  if (fRnrStyle== 0) fRnrStyle = new TrackRnrStyle;
  SetMainColorPtr(&fRnrStyle->fColor);
}

TrackList::TrackList(Int_t n_tracks, TrackRnrStyle* rs) :
  RenderElementListBase(),
  TPolyMarker3D(n_tracks),

  fTitle(),

  fRnrStyle   (rs),
  fRnrMarkers (kTRUE),
  fRnrTracks  (kTRUE)
{
  Init();
}

TrackList::TrackList(const Text_t* name, Int_t n_tracks, TrackRnrStyle* rs) :
  RenderElementListBase(),
  TPolyMarker3D(n_tracks),
  
  fTitle(),

  fRnrStyle   (rs),
  fRnrMarkers (kTRUE),
  fRnrTracks  (kTRUE)
{
  Init();
  SetName(name);
}

void TrackList::Reset(Int_t n_tracks)
{
  delete [] fP; fP = 0;
  fN = n_tracks;
  if(fN) fP = new Float_t [3*fN];
  memset(fP, 0, 3*fN*sizeof(Float_t));
  fLastPoint = -1;
}

/**************************************************************************/

void TrackList::Paint(Option_t* option)
{
  if(fRnrElement) {
    if(fRnrMarkers) {
      TPolyMarker3D::Paint(option);
    }
    if(fRnrTracks) {
      for(List_i i=fChildren.begin(); i!=fChildren.end(); ++i) {
	if((*i)->GetRnrElement())
	  (*i)->GetObject()->Paint(option);
      }
    }
  }
}

/**************************************************************************/

void TrackList::AddElement(RenderElement* el)
{
  static const Exc_t eH("TrackList::AddElement ");
  if (dynamic_cast<Track*>(el)  == 0)
    throw(eH + "new element not a Track.");
  RenderElementListBase::AddElement(el);
}

/**************************************************************************/

void TrackList::SetRnrMarkers(Bool_t rnr)
{
  fRnrMarkers = rnr;
  gReve->Redraw3D();
}

void TrackList::SetRnrTracks(Bool_t rnr)
{

  fRnrTracks = rnr;
  gReve->Redraw3D();
}

/**************************************************************************/

void TrackList::MakeTracks()
{
  for(List_i i=fChildren.begin(); i!=fChildren.end(); ++i) {
    ((Track*)(*i))->MakeTrack();
  }
  gReve->Redraw3D();
}


void TrackList::MakeMarkers()
{
  Reset(fChildren.size());
  for(List_i i=fChildren.begin(); i!=fChildren.end(); ++i) {
    Track& t = *((Track*)(*i));
    if(t.GetN() > 0)
      SetNextPoint(t.fV.x, t.fV.y, t.fV.z);
  }
  gReve->Redraw3D();
}

/**************************************************************************/
/*************************************************************************/

void TrackList::SetWidth(Width_t w)
{
  Width_t oldw = fRnrStyle->fWidth;
  fRnrStyle->fWidth = w;
  for (List_i i=fChildren.begin(); i!=fChildren.end(); ++i) {
    Track& t = *((Track*)(*i));
    if (t.GetLineWidth() == oldw)
      t.SetLineWidth(w);
  }
}

void TrackList::SetMaxR(Float_t x)
{
  fRnrStyle->fMaxR = x;
  MakeTracks();
  MakeMarkers();
}

void TrackList::SetMaxZ(Float_t x)
{
  fRnrStyle->fMaxZ = x;
  MakeTracks();
  MakeMarkers();
}

void TrackList::SetMaxOrbs(Float_t x)
{
  fRnrStyle->fMaxOrbs = x;
  MakeTracks();
}

void TrackList::SetMinAng(Float_t x)
{
  fRnrStyle->fMinAng = x;
  MakeTracks();
}

void TrackList::SetDelta(Float_t x)
{
  fRnrStyle->fDelta = x;
  MakeTracks();
}

void TrackList::SetFitDaughters(Bool_t x)
{
  fRnrStyle->fFitDaughters = x;
  MakeTracks();
}

void TrackList::SetFitDecay(Bool_t x)
{
  fRnrStyle->fFitDecay = x;
  MakeTracks();
}

/**************************************************************************/
/**************************************************************************/

void TrackList::SelectByPt(Float_t min_pt, Float_t max_pt)
{
  Float_t minptsq = min_pt*min_pt;
  Float_t maxptsq = max_pt*max_pt;
  Float_t ptsq;

  for(List_i i=fChildren.begin(); i!=fChildren.end(); ++i) {
    ptsq = ((Track*)(*i))->fP.Perp2();
    (*i)->SetRnrElement(ptsq >= minptsq && ptsq <= maxptsq);
  }
}

/**************************************************************************/

void TrackList::ImportHits()
{
  for(List_i i=fChildren.begin(); i!=fChildren.end(); ++i) {
    ((Track*)(*i))->ImportHits();
  }
}

void TrackList::ImportClusters()
{
  for(List_i i=fChildren.begin(); i!=fChildren.end(); ++i) {
    ((Track*)(*i))->ImportClusters();
  }
}

/**************************************************************************/
/**************************************************************************/
/**************************************************************************/

#include "RGEditor.h"

//______________________________________________________________________
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
    tlist->RemoveParent(this);
  fTrackLists.Clear();
}

void TrackCounter::RegisterTracks(TrackList* tlist, Bool_t goodTracks)
{
  // printf("TrackCounter::RegisterTracks '%s', %s\n",
  //   tlist->GetObject()->GetName(), goodTracks ? "good" : "bad");

  tlist->AddParent(this);
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
  // !!!! But then ... should also sore local information if track is ok.

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
      gReve->Redraw3D();

      printf("TrackCounter::CountTrack All=%d, Good=%d, Bad=%d\n",
	     fAllTracks, fGoodTracks, fAllTracks-fGoodTracks);

      if (gReve->GetEditor()->GetModel() == GetObject())
	gReve->EditRenderElement(this);

      break;
    }

  } // end switch fClickAction
}
