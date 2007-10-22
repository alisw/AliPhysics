// $Header$

#include "TrackGL.h"
#include <Reve/Track.h>
#include <Reve/GLUtilNS.h>

#include <TGLSelectRecord.h>

#include <TGLIncludes.h>

using namespace Reve;

//______________________________________________________________________
// TrackGL
//

ClassImp(TrackGL)

TrackGL::TrackGL() : LineGL()
{
  // fDLCache = false; // Disable display list.
}

TrackGL::~TrackGL()
{}

/**************************************************************************/

Bool_t TrackGL::SetModel(TObject* obj, const Option_t* /*opt*/)
{
  if(LineGL::SetModel(obj) == kFALSE) return kFALSE;
  if(SetModelCheckClass(obj, Track::Class())) {
    fTrack = dynamic_cast<Track*>(obj);
    return kTRUE;
  }
  return kFALSE;
}
/**************************************************************************/

void TrackGL::ProcessSelection(TGLRnrCtx & /*rnrCtx*/, TGLSelectRecord & rec)
{
  // Processes secondary selection from TGLViewer.
  // Calls TPointSet3D::PointSelected(Int_t) with index of selected
  // point as an argument.

  printf("TrackGL::ProcessSelection %d names on the stack (z1=%g, z2=%g).\n",
	 rec.GetN(), rec.GetMinZ(), rec.GetMaxZ());
  printf("  Names: ");
  for (Int_t j=0; j<rec.GetN(); ++j) printf ("%d ", rec.GetItem(j));
  printf("\n");

  ((Track*)fM)->CtrlClicked((Track*)fM);
}

/**************************************************************************/
void TrackGL::DirectDraw(TGLRnrCtx & rnrCtx) const
{
  LineGL::DirectDraw(rnrCtx);  

  // path-marks
  std::vector<PathMark*>& pm = fTrack->fPathMarks;
  TrackRnrStyle& RS = *fTrack->GetRnrStyle();
  if(pm.size())
  {
    Float_t* pnts = new Float_t[3*pm.size()]; // maximum
    Int_t N = 0;  
    Bool_t accept;
    for(std::vector<PathMark*>::iterator i=pm.begin(); i!=pm.end(); ++i) 
    {
      accept = kFALSE;
      switch((*i)->type)
      {
	case(PathMark::Daughter):
	  if(RS.fRnrDaughters) accept = kTRUE;
	  break;
	case(PathMark::Reference):
	  if(RS.fRnrReferences) accept = kTRUE;
	  break;
	case(PathMark::Decay):
	  if(RS.fRnrDecay) accept = kTRUE;
	  break;
      } 
      if(accept)
      {
	if((TMath::Abs((*i)->V.z) < RS.fMaxZ) && ((*i)->V.Perp() < RS.fMaxR))
	{
	  pnts[3*N  ] =(*i)->V.x; 
	  pnts[3*N+1] =(*i)->V.y; 
	  pnts[3*N+2] =(*i)->V.z;
	  N++;
	}
      }
    } 
    GLUtilNS::RenderPolyMarkers(RS.fPMAtt, pnts, N);
    delete [] pnts;
  }

  // fist vertex
 if(RS.fRnrFV && fTrack->GetLastPoint())
    GLUtilNS::RenderPolyMarkers(RS.fFVAtt, fTrack->GetP(), 1);
}
