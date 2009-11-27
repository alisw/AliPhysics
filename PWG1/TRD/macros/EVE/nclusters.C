#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TRD/AliTRDcluster.h>
#include <TRD/AliTRDseedV1.h>
#include <TRD/AliTRDtrackV1.h>
#endif

void nclusters(const TObject* object, Double_t* &res, Int_t& n)
{
  if (!object) return;
  if (object->IsA() != AliTRDtrackV1::Class()) return;

  const AliTRDtrackV1* track = dynamic_cast<const AliTRDtrackV1*>(object); 
  if (!track)  return;

  n = 1;
  res = new Double_t[n];
  res[0] = track->GetNumberOfClusters();
  return;
}

