#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TRD/AliTRDgeometry.h>
#include <TRD/AliTRDcluster.h>
#include <TRD/AliTRDseedV1.h>
#include <TRD/AliTRDtrackV1.h>
#endif

void clusterResiduals(const AliTRDtrackV1* track, Double_t* &res, Int_t& n)
{
  if (!track)  return;

  Int_t in = 0;
  n = track->GetNumberOfClusters();
  res = new Double_t[n];
  memset(res, 0, n*sizeof(Double_t));
  
  AliTRDseedV1 *fTracklet = 0x0;
  for(Int_t ily=0; ily<AliTRDgeometry::kNlayer; ily++){
    if(!(fTracklet = track->GetTracklet(ily))) continue;
    if(!fTracklet->IsOK()) continue;
    if(!fTracklet->Fit(kFALSE)) continue;
      
    AliTRDcluster *c = 0x0;
    for(Int_t ic=AliTRDseed::knTimebins-1; ic>=0; ic--){
      if(!(c = fTracklet->GetClusters(ic))) continue;
      
      res[in++] = fTracklet->GetYat(c->GetX()) - c->GetY();
    }
  }
  return;
}
