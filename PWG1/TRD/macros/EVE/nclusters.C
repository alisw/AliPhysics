#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TRD/AliTRDcluster.h>
#include <TRD/AliTRDseedV1.h>
#include <TRD/AliTRDtrackV1.h>
#endif

void nclusters(const AliTRDtrackV1* track, Double_t* &res, Int_t& n)
{
  if (!track)  return;

  n = 1;
  res = new Double_t[n];
  res[0] = track->GetNumberOfClusters();
  return;
}

