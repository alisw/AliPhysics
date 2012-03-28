#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TRD/AliTRDcluster.h>
#include <TRD/AliTRDseedV1.h>
#include <TRD/AliTRDtrackV1.h>
#endif

void chargeDistr(const AliTRDtrackV1* track, Double_t* results, Int_t& nResults)
{
  if (!track){
    Error("chargeDistr()", "Missing track.");
    return;
  }
  Int_t Nt = track->GetNumberOfTracklets();
  AliTRDcluster* cls(NULL);
  AliTRDseedV1 *tracklet(NULL);

  // Count clusters
  Int_t nCls = 0;
  for (Int_t ily = 0; ily < 6; ily++) {
    if(!(tracklet = track->GetTracklet(ily))){
      //Warning("chargeDistr()", "Missing tracklet in ly[%d]", ily);
      continue;
    }
    if(!tracklet->IsOK()){
      //Warning("chargeDistr()", "Bad tracklet in ly[%d]", ily);
      continue;
    }
    for (Int_t icl = 0; icl < AliTRDseedV1::kNclusters; icl++) {
      if(!(cls = tracklet->GetClusters(icl))){
        //Warning("chargeDistr()", "Missing cls[%2d] tracklet in ly[%d]", icl, ily);
        continue;
      }
      nCls++;
    }
  }
  
  // Nothing to do?
  if(!nCls){
    Warning("chargeDistr()", "Missing clusters.");
    return;
  }
  //Info("chargeDistr()", "Found %3d clusters.", nCls);

  // Allocate memory for the results (AliEveTRDAnalyseObjectList will clean this memory automatically)
  //results = new Double_t[nResults];
  nResults = 2;
  for (Int_t i = 0; i < nResults; i++)  results[i] = 0.;
  Int_t currentIndex = 0;

  for (Int_t trackletInd = 0; trackletInd < Nt && currentIndex < nResults; trackletInd++){
    if(!(tracklet = track->GetTracklet(trackletInd))) continue;
    if(!tracklet->IsOK()) continue;
    for (Int_t clusterInd = 0; clusterInd < AliTRDseedV1::kNclusters; clusterInd++){
      if(!(cls = tracklet->GetClusters(clusterInd))) continue;
      results[1] += cls->GetQ();
    }
  }
  results[0] =nCls;
  results[1]/=nCls;
}
