#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TRD/AliTRDcluster.h>
#include <TRD/AliTRDseedV1.h>
#include <TRD/AliTRDtrackV1.h>
#endif

void nclusters(const AliTRDtrackV1* track, Double_t* res, Int_t& n)
{
  if (!track){
    Error("nclusters()", "Missing track.");
    return;
  }
  n = 3;
  res[0] = track->GetNumberOfClusters();
  res[1] = track->GetNumberOfTracklets();
  res[2] = track->Pt();
  //printf("-> ncl = %3d ntrklt[%d] pt[%f]\n", Int_t(res[0]), Int_t(res[1]), res[2]);
  return;
}

