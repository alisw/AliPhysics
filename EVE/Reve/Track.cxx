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
  RenderElement(),
  TPolyLine3D(),

  fV(),
  fP(),
  fBeta(0),
  fCharge(0),
  fLabel(0),
  fPathMarks(),

  fRnrStyle(0),

  fName(),
  fTitle()
{}

Track::Track(Reve::MCTrack* t, TrackRnrStyle* rs):
  RenderElement(),
  TPolyLine3D(),

  fV(t->Vx(), t->Vy(), t->Vz()),
  fP(t->Px(), t->Py(), t->Pz()),
  fBeta(t->P()/t->Energy()),
  fCharge(0),
  fLabel(t->label),
  fPathMarks(),

  fRnrStyle(rs),

  fName(t->GetName()),
  fTitle()
{
  fLineColor = fRnrStyle->GetColor();
  fMainColorPtr = &fLineColor;

  TParticlePDG* pdgp = t->GetPDG();
  if(pdgp == 0) {
    t->ResetPdgCode(); pdgp = t->GetPDG();
  }
  fCharge = (Int_t) TMath::Nint(pdgp->Charge()/3);
}

Track::Track(Reve::RecTrack* t, TrackRnrStyle* rs) :
  RenderElement(),
  TPolyLine3D(),

  fV(t->V),
  fP(t->P),
  fBeta(t->beta),
  fCharge(t->sign),
  fLabel(t->label),
  fPathMarks(),

  fRnrStyle(rs),

  fName(t->GetName()),
  fTitle()
{
  fLineColor = fRnrStyle->GetColor();
  fMainColorPtr = &fLineColor;
}

Track::~Track()
{
  for (vpPathMark_i i=fPathMarks.begin(); i!=fPathMarks.end(); ++i)
    delete *i;
}

void Track::Reset(Int_t n_points)
{
  delete [] TPolyLine3D::fP; TPolyLine3D::fP = 0;
  fN = n_points;
  if(fN) TPolyLine3D::fP = new Float_t [3*fN];
  memset(TPolyLine3D::fP, 0, 3*fN*sizeof(Float_t));
  fLastPoint = -1;
}

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
      for(lpRE_i i=fChildren.begin(); i!=fChildren.end(); ++i) {
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
  for(lpRE_i i=fChildren.begin(); i!=fChildren.end(); ++i) {
    ((Track*)(*i))->MakeTrack();
  }
  gReve->Redraw3D();
}


void TrackList::MakeMarkers()
{
  Reset(fChildren.size());
  for(lpRE_i i=fChildren.begin(); i!=fChildren.end(); ++i) {
    Track& t = *((Track*)(*i));
    if(t.GetN() > 0)
      SetNextPoint(t.fV.x, t.fV.y, t.fV.z);
  }
  gReve->Redraw3D();
}

/**************************************************************************/
/*************************************************************************/

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

  for(lpRE_i i=fChildren.begin(); i!=fChildren.end(); ++i) {
    ptsq = ((Track*)(*i))->fP.Perp2();
    (*i)->SetRnrElement(ptsq >= minptsq && ptsq <= maxptsq);
  }
}

/**************************************************************************/

void TrackList::ImportHits()
{
  for(lpRE_i i=fChildren.begin(); i!=fChildren.end(); ++i) {
    ((Track*)(*i))->ImportHits();
  }
}

void TrackList::ImportClusters()
{
  for(lpRE_i i=fChildren.begin(); i!=fChildren.end(); ++i) {
    ((Track*)(*i))->ImportClusters();
  }
}
