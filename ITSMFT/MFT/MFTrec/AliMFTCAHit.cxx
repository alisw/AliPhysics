#include "AliMFTCAHit.h"

ClassImp(AliMFTCAHit)


//___________________________________________________________________________
AliMFTCAHit::AliMFTCAHit() :
TObject(),
fTrackGID(-1),
fLayer(-1),
fID(-1),
fPos(),
fIsUsed(kFALSE),
fDetElemID(-1),
fNRoads(0),
fInRoad(),
fIsFace(0),
fNInL(-1)
{
  
  for (Int_t i = 0; i < 100; i++) fInRoad[i] = -1;
  
}

//___________________________________________________________________________
AliMFTCAHit::AliMFTCAHit(const AliMFTCAHit &hit) :
TObject(hit),
fTrackGID(hit.fTrackGID),
fLayer(hit.fLayer),
fID(hit.fID),
fIsUsed(hit.fIsUsed),
fDetElemID(hit.fDetElemID),
fNRoads(hit.fNRoads),
fIsFace(hit.fIsFace),
fNInL(hit.fNInL)
{
  
  // copy constructor
  
  fPos[0] = hit.fPos[0]; fPos[1] = hit.fPos[1]; fPos[2] = hit.fPos[2];
  
  for (Int_t i = 0; i < hit.fNRoads; i++) {
    fInRoad[i] = hit.fInRoad[i];
  }
  
}
