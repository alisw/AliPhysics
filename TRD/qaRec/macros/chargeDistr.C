#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TRD/AliTRDcluster.h>
#include <TRD/AliTRDseedV1.h>
#include <TRD/AliTRDtrackV1.h>
#endif

void chargeDistr(const AliTRDtrackV1* track, Double_t* &results, Int_t& nResults)
{
  if (!track)  return;
  
  Int_t Nt = track->GetNumberOfTracklets();
  AliTRDcluster* cls = 0;
  AliTRDseedV1 *tracklet = 0x0;

  // Count clusters
  nResults = 0;
  for (Int_t trackletInd = 0; trackletInd < Nt; trackletInd++)
  {
    if(!(tracklet = track->GetTracklet(trackletInd))) continue;
    if(!tracklet->IsOK()) continue;
    
    for (Int_t clusterInd = 0; clusterInd < AliTRDseed::knTimebins; clusterInd++) 
    {
      if(!(cls = tracklet->GetClusters(clusterInd))) continue;
            
	    nResults++;
	  }
  }
  
  // Nothing to do?
  if (nResults == 0)  return;

  // Allocate memory for the results (AliEveTRDTrackList will clean this memory automatically)
  results = new Double_t[nResults];
  for (Int_t i = 0; i < nResults; i++)  results[i] = -100;
  Int_t currentIndex = 0;

  for (Int_t trackletInd = 0; trackletInd < Nt && currentIndex < nResults; trackletInd++)
  {
    if(!(tracklet = track->GetTracklet(trackletInd))) continue;
    if(!tracklet->IsOK()) continue;
    
    for (Int_t clusterInd = 0; clusterInd < AliTRDseed::knTimebins; clusterInd++) 
    {
      if(!(cls = tracklet->GetClusters(clusterInd))) continue;
            
	    results[currentIndex++] = cls->GetQ();
	  }
  }
}
