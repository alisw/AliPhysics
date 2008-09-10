// Example macro to select TRD tracks
// The macro selects for displaying only the tracks 
// with less than 60 clusters
// 
// Authors:
// B.Hess <Hess@Stud.Uni-Heidelberg.de>
// A.Bercuci <A.Bercuci@gsi.de>
//
Bool_t clusterSelection(AliTRDtrackV1* track)
{
  if (track->GetNumberOfClusters() > 60) return kFALSE;
  return kTRUE;
}
