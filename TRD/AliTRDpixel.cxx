///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Contains the information for one TRD pixel                               //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliTRDpixel.h"

ClassImp(AliTRDpixel)

//_____________________________________________________________________________
AliTRDpixel::AliTRDpixel():TObject()
{
  //
  // Create a TRD pixel
  // 

  fSignal   = 0;
  fTrack[0] = 0;
  fTrack[1] = 0;
  fTrack[2] = 0;

}
