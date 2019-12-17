///////////////////////////////////////////////////////////////////////////
//                                                                       //
// AliFemtoWRzTrackCut: An extension to AliFemtoMJTrackCut               //
// The method allows you to set different cut on DCA for different pt ranges //
///////////////////////////////////////////////////////////////////////////

#ifndef AliFemtoWRzTrackCut_hh
#define AliFemtoWRzTrackCut_hh


#include "AliFemtoMJTrackCut.h"

class AliFemtoWRzTrackCut : public AliFemtoMJTrackCut{

public:
   AliFemtoWRzTrackCut();
   AliFemtoWRzTrackCut(const AliFemtoWRzTrackCut &aCut);
   virtual ~AliFemtoWRzTrackCut();
   AliFemtoWRzTrackCut& operator =(const AliFemtoWRzTrackCut &aCut);
   virtual bool Pass(const AliFemtoTrack* aTrack);

   void SetDCAcutPtDependence(const float &DCAxy,const float &DCAz,const float& Pt_min, const float& Pt_max);
   
private:
   float fPt_dca[2]; //bounds for transverse momentum to DCA settings
   float fDCA_ptdependence[2]; //DCA values (cm) in XY and Z plane for given ranges of pt 
};

inline void AliFemtoWRzTrackCut::SetDCAcutPtDependence(const float &DCAxy, const float &DCAz, const float& Pt_min,const float& Pt_max){ 
fPt_dca[0]=Pt_min, fPt_dca[1]=Pt_max,fDCA_ptdependence[0]=DCAxy,fDCA_ptdependence[1]=DCAz;
}


#endif

