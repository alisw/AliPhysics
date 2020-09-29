#include "AliFemtoWRzTrackCut.h"


AliFemtoWRzTrackCut::AliFemtoWRzTrackCut():
  AliFemtoMJTrackCut()
  {
     fDCA_ptdependence[0]=2.4; fDCA_ptdependence[4]=3.2;
     fDCA_ptdependence[1]=2.4; fDCA_ptdependence[5]=3.2;
     fDCA_ptdependence[2]=2.4; fDCA_ptdependence[6]=3.2;
     fDCA_ptdependence[3]=2.4; fDCA_ptdependence[7]=3.2;
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


void AliFemtoWRzTrackCut::SetDCAcutPtDependence(const float &DCAxy1,const float &DCAxy2, const float &DCAxy3, const float &DCAxy4, const float &DCAz1, const float &DCAz2, const float &DCAz3, const float &DCAz4)
{
  fDCA_ptdependence[0]=DCAxy1;
  fDCA_ptdependence[1]=DCAxy2;
  fDCA_ptdependence[2]=DCAxy3;
  fDCA_ptdependence[3]=DCAxy4;
  fDCA_ptdependence[4]=DCAz1;
  fDCA_ptdependence[5]=DCAz2;
  fDCA_ptdependence[6]=DCAz3;
  fDCA_ptdependence[7]=DCAz4;
}


bool AliFemtoWRzTrackCut::Pass( const AliFemtoTrack* track)
{
    AliFemtoMJTrackCut::Pass(track);

    float tPt = ::sqrt((track->P().x())*(track->P().x())+(track->P().y())*(track->P().y()));


    if (tPt > 0.5 && tPt < 1.0){
       if( (TMath::Abs(track->ImpactD()) > fDCA_ptdependence[0]) || (TMath::Abs(track->ImpactZ()) > fDCA_ptdependence[4]))
          return false;
    }
    else if (tPt >= 1.0 && tPt < 1.5){
       if( (TMath::Abs(track->ImpactD()) > fDCA_ptdependence[1]) || (TMath::Abs(track->ImpactZ()) > fDCA_ptdependence[5]))
          return false;
    }
    else if (tPt >= 1.5 && tPt < 2.0){
       if( (TMath::Abs(track->ImpactD()) > fDCA_ptdependence[2]) || (TMath::Abs(track->ImpactZ()) > fDCA_ptdependence[6]))
          return false;
    }
    else if (tPt >= 2.0 && tPt < 4.0){
       if( (TMath::Abs(track->ImpactD()) > fDCA_ptdependence[3]) || (TMath::Abs(track->ImpactZ()) > fDCA_ptdependence[7]))
          return false;
    }


}
