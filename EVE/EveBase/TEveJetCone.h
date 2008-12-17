// @(#)root/eve:$Id$
// Author: Matevz Tadel 2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TEveJetCone
#define ROOT_TEveJetCone

#include "TEveElement.h"
#include "TEveVSDStructs.h"
#include "TAttBBox.h"

class TEveJetCone : public TEveElementList,
                    public TAttBBox
{
   friend class TEveJetConeGL;

private:
   TEveJetCone(const TEveJetCone&);            // Not implemented
   TEveJetCone& operator=(const TEveJetCone&); // Not implemented

protected:
   typedef std::vector<TEveVector>        vTEveVector_t;
   typedef vTEveVector_t::iterator        vTEveVector_i;
   typedef vTEveVector_t::const_iterator  vTEveVector_ci;

   TEveVector      fApex;
   vTEveVector_t   fBasePoints;

public:
   TEveJetCone(const Text_t* n="TEveJetCone", const Text_t* t="");
   virtual ~TEveJetCone() {}

   void SetApex(const TEveVector& a)      { fApex = a; }
   void AddBasePoint(const TEveVector& p) { fBasePoints.push_back(p); }

   // void SetBaseFromEtaPhi(radius, eta, phi, deta, dphi);

   virtual Bool_t  CanEditMainTransparency() const { return kTRUE; }

   // For TAttBBox:
   virtual void ComputeBBox();
   // If painting is needed:
   virtual void Paint(Option_t* option="");

   ClassDef(TEveJetCone, 0); // Short description.
};

#endif
