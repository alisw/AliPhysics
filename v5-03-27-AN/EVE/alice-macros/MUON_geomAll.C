// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

/// \ingroup evemacros
/// \file MUON_geomAll.C
///
/// \author B. Vulpescu, LPC

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TGeoManager.h>
#include <TGeoNode.h>
#include <TEveManager.h>
#include <TEveGeoNode.h>
#endif

void MUON_geomAll()
{
  gGeoManager = gEve->GetGeometry("geometry.root");

  TEveGeoTopNode* topn_re = new TEveGeoTopNode
    (gGeoManager, gGeoManager->GetTopNode());

  gEve->AddGlobalElement(topn_re);

  gEve->Redraw3D(kTRUE);

}
