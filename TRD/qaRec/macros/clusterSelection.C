// Example macro to select TRD tracks
// The macro selects for displaying only the tracks 
// with less than 60 clusters
// 
// Authors:
// B.Hess <Hess@Stud.Uni-Heidelberg.de>
// A.Bercuci <A.Bercuci@gsi.de>
//

#ifndef __CINT__
#include <AliTRDtrackV1.h>
#endif

Bool_t clusterSelection(const AliTRDtrackV1* track=0x0)
{
  if(!track) return kFALSE;
  if (track->GetNumberOfClusters() > 60) return kFALSE;
  return kTRUE;
}
