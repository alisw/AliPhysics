/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

////////////////////////////////////////////////////////////////////////////
//                                                                        
//       === Class for properties specific to pad regions ===
//
//    Each segment of the TPC (i.e. IROC and corresponding OROC) consists
//    of three different pad sizes (short, medium, long). This class
//    is useful for scenarios, where it is appropriate to describe
//    some behaviour per pad size region. It provides an easy interface
//    for getting and setting arbitrary objects for each pad size region.
//    There is no need that this object is of the same type for each
//    pad size region (though it probably will be in most of the cases),
//    nor that it is set at all (e.g. when no data for this region
//    exists and such an object is not needed).
//
//    An example that makes usage of this class is the AliTPCFitPad class
//    which stores TLinearFitter objects for each pad region.
//
////////////////////////////////////////////////////////////////////////////

#include "AliTPCCalPadRegion.h"
#include "AliTPCROC.h"

ClassImp(AliTPCCalPadRegion)

AliTPCCalPadRegion::AliTPCCalPadRegion():
   TNamed(),
   fObjects(0)
{
   //
   // Default constructor.
   //
}

AliTPCCalPadRegion::AliTPCCalPadRegion(const char *name, const char *title) :
   TNamed(name, title),
   fObjects(0)
{
   //
   // Constructor.
   //

   fObjects = new TObjArray(fgkNSegments * fgkNPadTypes);
   fObjects->SetOwner(kTRUE);
}

AliTPCCalPadRegion::AliTPCCalPadRegion(const AliTPCCalPadRegion& obj) :
  TNamed(obj),
  fObjects(0)
{
   //
   // Copy constructor.
   //

   fObjects = new TObjArray(*(obj.fObjects));
   fObjects->SetOwner(kTRUE);
}

AliTPCCalPadRegion& AliTPCCalPadRegion::operator=(const AliTPCCalPadRegion& rhs) {
   //
   // Assignment operator.
   //

   if (this != &rhs) {
      TNamed::operator=(rhs);
      fObjects = new TObjArray(*(rhs.fObjects));
   }
   return *this;
}


void       AliTPCCalPadRegion::SetObject(TObject* obj, UInt_t segment, UInt_t padType)
{
  //
  // Set the object for given segment
  //
  if (!fObjects) {
    fObjects = new TObjArray(fgkNSegments * fgkNPadTypes);
    fObjects->SetOwner(kTRUE);
  }
  if (fObjects->GetEntriesFast()<Int_t(fgkNSegments * fgkNPadTypes)){
    fObjects->Expand(fgkNSegments * fgkNPadTypes);
  }
  if (BoundsOk("SetObject", segment, padType)){ 
    if (segment+fgkNSegments*padType>static_cast<UInt_t>(fObjects->GetEntriesFast())) fObjects->Expand(fgkNSegments * fgkNPadTypes);
    fObjects->AddAt(obj, segment+fgkNSegments*padType); 
  }
}

TObject*   AliTPCCalPadRegion::GetObject(UInt_t segment, UInt_t padType){  
  //
  //
  //
  if (fObjects->GetEntriesFast()<Int_t(fgkNSegments * fgkNPadTypes)){
    fObjects->Expand(fgkNSegments * fgkNPadTypes);
  }
  return fObjects->At(segment+fgkNSegments*padType); 
}



void AliTPCCalPadRegion::GetPadRegionCenterLocal(UInt_t padType, Double_t* xy) {
   //
   // Return the center of the pad size region in local
   // coordinates as an Double_t array xy of length 2.
   //
   
   Float_t centerPad[3] = {0};
   AliTPCROC* tpcROC = AliTPCROC::Instance();

   Int_t IOROC = (padType == 0) ? 0 : tpcROC->GetNInnerSector();
   //tpcROC->GetPositionLocal(IOROC, tpcROC->GetNRows(IOROC)/2, tpcROC->GetNPads(IOROC, tpcROC->GetNRows(IOROC)/2)/2, centerPad);  // use this instead of the switch statement if you want to calculate the center of the ROC and not the center of the regions with the same pad size
   switch (padType) {
      case 0:  // short pads
         tpcROC->GetPositionLocal(IOROC, tpcROC->GetNRows(IOROC)/2, tpcROC->GetNPads(IOROC, tpcROC->GetNRows(IOROC)/2)/2, centerPad);
         break;
      case 1:  // medium pads
         tpcROC->GetPositionLocal(IOROC, 64/2, tpcROC->GetNPads(IOROC, 64/2)/2, centerPad);
         break;
      case 2:  // long pads
         tpcROC->GetPositionLocal(IOROC, 64+32/2, tpcROC->GetNPads(IOROC, 64+32/2)/2, centerPad);
         break;
   }

   xy[0] = centerPad[0];
   xy[1] = centerPad[1];
}

/*UInt_t AliTPCCalPadRegion::GetStartRow(UInt_t padType) {
   //
   // Returns the index of the 
   //
}

UInt_t AliTPCCalPadRegion::GetEndRow(UInt_t padType) {

}*/
