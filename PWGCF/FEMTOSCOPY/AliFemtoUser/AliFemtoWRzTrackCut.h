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

   void SetDCAcutPtDependence(const float &DCAxy1,const float &DCAxy2,const float &DCAxy3,const float &DCAxy4, const float &DCAz1, const float &DCAz2, const float &DCAz3, const float &DCAz4);
   
private:
   float fDCA_ptdependence[8]; //DCA values (cm) in XY and Z plane for given ranges of pt 
};



#endif

