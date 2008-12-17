// @(#)root/eve:$Id$
// Author: Matevz Tadel 2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TEveJetCone.h"
#include "TEveTrans.h"

#include "TBuffer3D.h"
#include "TBuffer3DTypes.h"
#include "TVirtualPad.h"
#include "TVirtualViewer3D.h"

//______________________________________________________________________________
// Description of TEveJetCone
//
// The base must have at least three points.

ClassImp(TEveJetCone)

//______________________________________________________________________________
TEveJetCone::TEveJetCone(const Text_t* n, const Text_t* t) :
   TEveElementList(n, t, kTRUE),
   TAttBBox(),
   fApex(),
   fBasePoints()
{
   // Constructor.

   fColor = kYellow;
}


/******************************************************************************/

//______________________________________________________________________________
void TEveJetCone::ComputeBBox()
{
   // Compute bounding-box of the data.

   BBoxInit();
   BBoxCheckPoint(fApex);
   for (vTEveVector_ci i = fBasePoints.begin(); i != fBasePoints.end(); ++i)
   {
      BBoxCheckPoint(*i);
   }
}

//______________________________________________________________________________
void TEveJetCone::Paint(Option_t*)
{
   // Paint object.
   // This is for direct rendering (using TEveJetConeGL class).

   static const TEveException eh("TEveJetCone::Paint ");

   if (fRnrSelf == kFALSE) return;

   TBuffer3D buff(TBuffer3DTypes::kGeneric);

   // Section kCore
   buff.fID           = this;
   buff.fColor        = GetMainColor();
   buff.fTransparency = GetMainTransparency();
   if (HasMainTrans())
      RefMainTrans().SetBuffer3D(buff);
   buff.SetSectionsValid(TBuffer3D::kCore);

   Int_t reqSections = gPad->GetViewer3D()->AddObject(buff);
   if (reqSections != TBuffer3D::kNone)
      Error(eh, "only direct GL rendering supported.");
}
