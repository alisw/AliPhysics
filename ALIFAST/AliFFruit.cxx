
//////////////////////////////////////////////////////////////////////////
//                                                                      //
// AliFFruit                                                            //
//                                                                      //
// Utility class to draw Electrons, photons, Jets, Clusters,etc         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <TROOT.h>
#include <TMath.h>

#include "AliFDisplay.h"
#include "AliFFruit.h"
#include "AliFast.h"

ClassImp(AliFFruit)


//_____________________________________________________________________________
AliFFruit::AliFFruit(TObject *obj, Float_t eta, Float_t phi, Float_t pt, Int_t type)
          : TPolyLine3D(2)
{
  // Create a fruit object.
  // Current implementation uses a 3-d polyline to visualize this fruit

   fFruit = obj;
   SetBit(kCanDelete);

   const Int_t color[7] = {0,7,3,2,6,4,0};
   const Int_t width[7] = {8,8,8,8,8,8,8};
   Int_t lwidth = width[type];
   AliFDisplay *display = (AliFDisplay*)gAliFast->Display();
   if (display->AllViews()) lwidth /= 2;
   const Float_t PTMAX = 100;
   if (pt <= 0) return;
   Float_t rin    = display->Rin();
   Float_t rout   = display->Rout();
   Float_t theta  = 2*TMath::ATan(TMath::Exp(-eta));
   Float_t tantet = TMath::Tan(theta);
   Float_t cosphi = TMath::Cos(phi);
   Float_t sinphi = TMath::Sin(phi);
   Float_t zz = pt/PTMAX;
   if (zz > 3) zz = 3;
   Float_t rex = rin + 3*zz*(rout - rin);
   Float_t z1,z2;
   if (eta != 0) {
      z1 = rin/tantet;
      z2 = rex/tantet;
   } else {
      z1 = z2 = 0;
   }
   SetPoint(0, rin*cosphi,rin*sinphi, z1);
   SetPoint(1, rex*cosphi,rex*sinphi, z2);
   SetLineColor(color[type]);
   SetLineWidth(width[type]);
}

//_____________________________________________________________________________
void AliFFruit::Delete(Option_t *)
{
//    Dummy

}

//______________________________________________________________________________
char *AliFFruit::GetObjectInfo(Int_t px, Int_t py)
{
   return fFruit->GetObjectInfo(px, py);
}























