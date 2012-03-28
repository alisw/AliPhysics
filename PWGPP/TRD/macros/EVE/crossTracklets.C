#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TRD/AliTRDseedV1.h>
#include <TRD/AliTRDtrackV1.h>
#include <TRD/AliTRDgeometry.h>
#endif


Bool_t crossTracklets(const AliTRDtrackV1* track)
{
  if (!track){
    Error("crossTracklets()", "Missing track.");
    return kFALSE;
  }
  AliTRDseedV1 *tracklet(NULL);
  Bool_t cross = kFALSE;
  for(Int_t ily=0; ily<AliTRDgeometry::kNlayer; ily++){
    if(!(tracklet = track->GetTracklet(ily))) continue;
    if(!tracklet->IsOK()) continue;
    if(tracklet->IsRowCross()){
      cross = kTRUE;
      break;
    }
  }

  return cross;
}
