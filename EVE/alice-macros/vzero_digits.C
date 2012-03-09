// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/


#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TClonesArray.h>
#include <TGeoManager.h>
#include <TMath.h>
#include <TStyle.h>
#include <TTree.h>
#include <TEveManager.h>
#include <TEveElement.h>
#include <TEvePointSet.h>
#include <TEveQuadSet.h>
#include <TEveTrans.h>

#include <AliRunLoader.h>
#include <AliVZEROdigit.h>
#include <AliEveEventManager.h>
#endif


//   fV0CHeight1  =  2.5; // height of cell 1
//   fV0CHeight2  =  4.4; // height of cell 2
//   fV0CHeight3  =  7.4; // height of cell 3
//   fV0CHeight4  = 12.5; // height of cell 4
//   fV0CRMin     =  4.6; // inner radius of box
//
//   fV0AR0       =  4.2; // Radius of hole
//   fV0AR1       =  7.6; // Maximun radius of 1st cell
//   fV0AR2       = 13.8; // Maximun radius of 2nd cell
//   fV0AR3       = 22.7; // Maximun radius of 3rd cell
//   fV0AR4       = 41.3; // Maximun radius of 4th cell

void vzero_digits()
{
  static const Float_t RadC[] = { 4.6, 7.1, 11.5, 18.9, 31.4 };
  static const Float_t RadA[] = { 4.2, 7.6, 13.8, 22.7, 41.4 };
  static const Float_t RadEps = 0.4;
  static const Float_t PhiEps = 0.025;
  static const Float_t PhiStp = TMath::TwoPi()/8.0;

  gStyle->SetPalette(1, 0);

  AliRunLoader* rl =  AliEveEventManager::AssertRunLoader();
  rl->LoadDigits("VZERO");

  TTree* dt = rl->GetTreeD("VZERO", false);
  TClonesArray* dca = 0;
  dt->SetBranchAddress("VZERODigit", &dca);
  dt->GetEntry(0);

  Float_t v[12];

  TEveQuadSet* modR = new TEveQuadSet("V0R");
  modR->Reset(TEveQuadSet::kQT_FreeQuad, kFALSE, 32);

  TEveQuadSet* modL = new TEveQuadSet("V0L");
  modL->Reset(TEveQuadSet::kQT_FreeQuad, kFALSE, 32);

  Int_t numEntr = dca->GetEntriesFast();
  for (Int_t entr=0; entr<numEntr; ++entr)
  {
    AliVZEROdigit* d = (AliVZEROdigit*) dca->UncheckedAt(entr);
    Int_t i = d->PMNumber();

    if (i < 32)   // AliEveV0 Right
    {
      TEveQuadSet* module = modR;
      Int_t ri = i / 8;
      Int_t pi = i % 8;
      Float_t minR = RadC[ri]  + RadEps, maxR = RadC[ri+1] - RadEps;
      Float_t minP = pi*PhiStp + PhiEps, maxP = (pi+1)*PhiStp - PhiEps;

      v[ 0] = minR*TMath::Cos(minP); v[ 1] = minR*TMath::Sin(minP); v[ 2] = 0;
      v[ 3] = maxR*TMath::Cos(minP); v[ 4] = maxR*TMath::Sin(minP); v[ 5] = 0;
      v[ 6] = maxR*TMath::Cos(maxP); v[ 7] = maxR*TMath::Sin(maxP); v[ 8] = 0;
      v[ 9] = minR*TMath::Cos(maxP); v[10] = minR*TMath::Sin(maxP); v[11] = 0;

      module->AddQuad(v);
      module->QuadValue(d->ADC());
      module->QuadId(d);
    }
    else          // AliEveV0 Left
    {
      TEveQuadSet* module = modL;
      Int_t ri = (i-32) / 8;
      Int_t pi = i % 8;
      Float_t minR = RadA[ri]  + RadEps, maxR = RadA[ri+1] - RadEps;
      Float_t minP = pi*PhiStp + PhiEps, maxP = (pi+1)*PhiStp - PhiEps;

      v[ 0] = minR*TMath::Cos(minP); v[ 1] = minR*TMath::Sin(minP); v[ 2] = 0;
      v[ 3] = maxR*TMath::Cos(minP); v[ 4] = maxR*TMath::Sin(minP); v[ 5] = 0;
      v[ 6] = maxR*TMath::Cos(maxP); v[ 7] = maxR*TMath::Sin(maxP); v[ 8] = 0;
      v[ 9] = minR*TMath::Cos(maxP); v[10] = minR*TMath::Sin(maxP); v[11] = 0;

      module->AddQuad(v);
      module->QuadValue(d->ADC());
      module->QuadId(d);
    }
  }

  modL->RefMainTrans().SetPos(0, 0, 324);
  modR->RefMainTrans().SetPos(0, 0, -84);

  gEve->AddElement(modL);
  gEve->AddElement(modR);

  gEve->Redraw3D();
}
