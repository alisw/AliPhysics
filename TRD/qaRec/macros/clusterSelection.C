// Example macro to select TRD tracks
// The macro selects for displaying only the tracks 
// with less than 60 clusters
// 
// Authors:
// B.Hess <Hess@Stud.Uni-Heidelberg.de>
// A.Bercuci <A.Bercuci@gsi.de>
//

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TRD/AliTRDtrackV1.h>
#endif

Bool_t clusterSelection(const AliTRDtrackV1 *track)
{
  if(!track) return kFALSE;
  if (track->GetNumberOfClusters() > 60) return kFALSE;
  return kTRUE;
}
