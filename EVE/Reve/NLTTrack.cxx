// $Header$

#include "NLTTrack.h"
#include <Reve/NLTProjector.h>
#include <Reve/PODs.h>

using namespace Reve;

//______________________________________________________________________
// NLTTrack
//

ClassImp(NLTTrack)

NLTTrack::NLTTrack() :
  Track     (),
  fOrigPnts(0),
  fProjection(0)
{}

/**************************************************************************/
void NLTTrack::SetProjection(NLTProjector* proj, NLTProjectable* model)
{
  NLTProjected::SetProjection(proj, model);
  Track* origTrack = dynamic_cast<Track*>(fProjectable);

  SetTrackParams(*origTrack);

  // unfortunately fPathMarks is a vector of PathMark pointers
  PathMark* pm;
  std::vector<PathMark*>& refs = origTrack->GetPathMarksRef();
  for(std::vector<PathMark*>::iterator i=refs.begin(); i!=refs.end(); ++i)
  {
   pm = new PathMark();
   pm->V = (*i)->V;
   pm->P = (*i)->P;
   pm->type = (*i)->type;
   pm->time = (*i)->time;
   fPathMarks.push_back(pm);
  }
}

/**************************************************************************/
void NLTTrack::UpdateProjection()
{
  fProjection = fProjector->GetProjection();
  MakeTrack(kFALSE); //NLTProjector makes recursive calls
}

//______________________________________________________________________________
void  NLTTrack::GetBreakPoint(Int_t idx, Bool_t back,  Float_t& x, Float_t& y, Float_t& z)
{
  Vector vL = fOrigPnts[idx];
  Vector vR = fOrigPnts[idx+1];
  // printf("out of tolerance:%d (%f, %f, %f)(%f, %f, %f) \n",
  // 	 idx, vL.x, vL.y, vL.z, vR.x, vR.y, vR.z );

  Vector vM, vLP, vRP, vMP;
  while((vL-vR).Mag() > 0.1)
  { 
    vM.Mult(vL+vR, 0.5f);
    vLP.Set(vL); fProjection->ProjectPoint(vLP.x, vLP.y, vLP.z);
    vMP.Set(vM); fProjection->ProjectPoint(vMP.x, vMP.y, vMP.z);
    // if(fProjection->AcceptSegment(vL, vM, GetRnrStyle()->fDelta*0.1))
    if(fProjection->AcceptSegment(vL, vM, 0.0f))
    {
      vL.Set(vM);
    }
    else 
    {
      vR.Set(vM);
    }
    //printf("new interval Mag %f (%f, %f, %f)(%f, %f, %f) \n",(vL-vR).Mag(), vL.x, vL.y, vL.z, vR.x, vR.y, vR.z);
  }

  if(back)
  {
    x = vL.x; y = vL.y; z = vL.z;
  }
  else
  {
    x = vR.x; y = vR.y; z = vR.z;
  }
  fProjection->ProjectPoint(x, y, z);
  // printf("NLTTrack::GetBreakPoint %d (%f, %f, %f) \n", idx, x, y, z);
}

//______________________________________________________________________________
Int_t  NLTTrack::GetBreakPointIdx(Int_t start)
{
  Int_t val = fLastPoint;

  Vector v1; Vector v2;
  if( Size() > 1 )
  {
    Bool_t broken = kFALSE;
    Int_t i = start;
    while(i < fLastPoint)
    {
      GetPoint(i,   v1.x, v1.y, v1.z);
      GetPoint(i+1, v2.x, v2.y, v2.z);
      if(fProjection->AcceptSegment(v1, v2, fRnrStyle->fDelta) == kFALSE)
      {
	broken = kTRUE;
        break;
      }
      i++;
    }
    if(broken) val = i;
  }
  // printf("BreakPoint IDX start:%d, BREAK %d,  total:%d \n", start, val, Size());
  return val;
}

/**************************************************************************/

void NLTTrack::MakeTrack(Bool_t recurse)
{
  Track::MakeTrack(recurse);

  fBreakPoints.clear();
  if(Size() == 0) return; // it is possible to be outside the limits of MaxR, MaxZ ...

  // poject line points
  Float_t *p = GetP();
  fOrigPnts = new Vector[Size()];
  for(Int_t i = 0; i < Size(); ++i, p+=3)
  {
    fOrigPnts[i].Set(p);
    fProjection->ProjectPoint(p[0], p[1], p[2]);
    p[2] = fDepth;
  } 
  Float_t x, y, z;
  std::vector<Vector> vvec;
  Int_t bL = 0, bR = GetBreakPointIdx(0); 
  while (1)
  {
    for(Int_t i=bL; i<=bR; i++)
    {
      GetPoint(i, x, y,z);
      vvec.push_back(Vector(x, y, z));
    }
    if (bR == fLastPoint)
      break;

    GetBreakPoint(bR, kTRUE,  x, y, z); vvec.push_back(Vector(x, y, z));
    fBreakPoints.push_back(vvec.size());
    GetBreakPoint(bR, kFALSE, x, y, z); vvec.push_back(Vector(x, y, z));

    bL = bR + 1;
    bR = GetBreakPointIdx(bL);
  }
  fBreakPoints.push_back(fLastPoint+1); // enforce drawing to end of line
  Reset(vvec.size());
  for (std::vector<Reve::Vector>::iterator i=vvec.begin(); i!=vvec.end(); ++i)
    SetNextPoint((*i).x, (*i).y, (*i).z); 
  delete [] fOrigPnts;
}

/**************************************************************************/
void NLTTrack::PrintLineSegments()
{
  printf("%s LineSegments:\n", GetName());
  Int_t start = 0;
  Int_t segment = 0;
  Vector S;
  Vector E;
  for (std::vector<Int_t>::iterator bpi = fBreakPoints.begin();
       bpi != fBreakPoints.end(); ++bpi)
  {
    Int_t size = *bpi - start;

    GetPoint(start, S.x, S.y, S.z);
    GetPoint((*bpi)-1, E.x, E.y, E.z);
    printf("seg %d size %d start %d ::(%f, %f, %f) (%f, %f, %f)\n", 
	   segment, size, start, S.x, S.y, S.z, E.x, E.y, E.z);
    start   += size;
    segment ++;
  }
}

/**************************************************************************/

void NLTTrack::CtrlClicked(Reve::Track* /*track*/)
{
  Track* t = dynamic_cast<Track*>(fProjectable);
  if (t)
    t->CtrlClicked(t);
}

//______________________________________________________________________
// NLTTrackList
//

ClassImp(NLTTrackList)

NLTTrackList::NLTTrackList() :
  TrackList    (),
  NLTProjected ()
{
}

/**************************************************************************/
void NLTTrackList::SetProjection(NLTProjector* proj, NLTProjectable* model)
{
  NLTProjected::SetProjection(proj, model);
 
  TrackList& tl   = * dynamic_cast<TrackList*>(model);
  SetLineColor(tl.GetLineColor());
  SetLineStyle(tl.GetLineStyle());
  SetLineWidth(tl.GetLineWidth());
  SetMarkerColor(tl.GetMarkerColor());
  SetMarkerStyle(tl.GetMarkerStyle());
  SetMarkerSize(tl.GetMarkerSize());
  SetRnrLine(tl.GetRnrLine());
  SetRnrPoints(tl.GetRnrPoints());
 
  SetRnrStyle(tl.GetRnrStyle());
}
