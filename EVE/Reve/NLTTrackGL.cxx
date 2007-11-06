// $Header$

#include "NLTTrackGL.h"
#include <Reve/NLTTrack.h>
#include <Reve/NLTProjector.h>
#include <Reve/GLUtilNS.h>

#include <TGLRnrCtx.h>
#include <TGLIncludes.h>

using namespace Reve;

//______________________________________________________________________
// NLTTrackGL
//

ClassImp(NLTTrackGL)

NLTTrackGL::NLTTrackGL() : TrackGL(), fM(0)
{
  // fDLCache = kFALSE; // Disable display list.
}

NLTTrackGL::~NLTTrackGL()
{}

/**************************************************************************/

Bool_t NLTTrackGL::SetModel(TObject* obj, const Option_t* /*opt*/)
{
  if(TrackGL::SetModel(obj) == kFALSE) return kFALSE;
  if(SetModelCheckClass(obj, NLTTrack::Class())) {
    fM = dynamic_cast<NLTTrack*>(obj);
    return kTRUE;
  }
  return kFALSE;
}

/**************************************************************************/

void NLTTrackGL::DirectDraw(TGLRnrCtx & rnrCtx) const
{
  // printf("NLTTrackGL::DirectDraw Style %d, LOD %d\n", flags.Style(), flags.LOD());
  if (rnrCtx.DrawPass() == TGLRnrCtx::kPassOutlineLine || fM->Size() == 0)
    return;

  // lines
  Int_t start = 0;
  Float_t* p = fM->GetP();
  for (std::vector<Int_t>::iterator bpi = fM->fBreakPoints.begin();
       bpi != fM->fBreakPoints.end(); ++bpi)
  {
    Int_t size = *bpi - start;
    if (fM->fRnrLine)   GLUtilNS::RenderLine(*fM, p, size);
    if (fM->fRnrPoints) GLUtilNS::RenderPolyMarkers(*fM, p, size);
    p     += 3*size;
    start +=   size;
  }

  // path-marks
  std::vector<PathMark*>& pm = fM->fPathMarks;
  TrackRnrStyle& RS = *fM->GetRnrStyle();
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
          fM->fProjection->ProjectPoint(pnts[3*N  ], pnts[3*N+1], pnts[3*N+2]);
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
