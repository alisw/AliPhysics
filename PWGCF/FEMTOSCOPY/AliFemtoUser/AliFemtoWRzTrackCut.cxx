#include "AliFemtoWRzTrackCut.h"


AliFemtoWRzTrackCut::AliFemtoWRzTrackCut():
  AliFemtoMJTrackCut()
  {
     fDCA_ptdependence[0]=2.4; fDCA_ptdependence[1]=3.2;
     fPt_dca[0]=0.5; fPt_dca[1]=4;
  }

AliFemtoWRzTrackCut::AliFemtoWRzTrackCut(const AliFemtoWRzTrackCut &aCut) : AliFemtoMJTrackCut(aCut)
{
    //copy constructor 
}

AliFemtoWRzTrackCut::~AliFemtoWRzTrackCut()
{
  // Destructor
}


AliFemtoWRzTrackCut& AliFemtoWRzTrackCut::operator=(const AliFemtoWRzTrackCut& aCut)
{
  // assignment operator
  if (this == &aCut)
    return *this;

  AliFemtoMJTrackCut::operator=(aCut);
 
  return *this;
}

bool AliFemtoWRzTrackCut::Pass( const AliFemtoTrack* track)
{
  AliFemtoMJTrackCut::Pass(track);

  float tPt = ::sqrt((track->P().x())*(track->P().x())+(track->P().y())*(track->P().y()));

  if (tPt>fPt_dca[0] && tPt <fPt_dca[1])
	if( (TMath::Abs(track->ImpactD()) > fDCA_ptdependence[0]) || (TMath::Abs(track->ImpactZ()) > fDCA_ptdependence[1]))
             return false;

}
