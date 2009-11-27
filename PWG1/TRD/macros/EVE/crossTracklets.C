#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TRD/AliTRDseedV1.h>
#include <TRD/AliTRDtrackV1.h>
#include <TRD/AliTRDgeometry.h>
#endif


Bool_t crossTracklets(const TObject* object)
{
  if (!object) return kFALSE;
  if (object->IsA() != AliTRDtrackV1::Class()) return kFALSE;

  const AliTRDtrackV1* track = dynamic_cast<const AliTRDtrackV1*>(object); 
  if(!track) return kFALSE;
  if(track->GetNumberOfTracklets() != AliTRDgeometry::kNlayer) return kFALSE;

  AliTRDseedV1 *tracklet = 0x0;
  Bool_t cross = kFALSE;
  for(Int_t ily=0; ily<AliTRDgeometry::kNlayer; ily++){
    if(!(tracklet = track->GetTracklet(ily))) continue;
    if(!tracklet->IsOK()) continue;
   // if(!tracklet->GetNChange()) continue;
    cross = kTRUE;
    break;
  }

  return cross;
}
