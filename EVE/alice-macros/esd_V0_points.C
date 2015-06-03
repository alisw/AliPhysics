// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TEveManager.h>
#include <TEvePointSet.h>

#include <AliESDEvent.h>
#include <AliESDv0.h>
#include <AliEveEventManager.h>
#endif

void esd_VO_fill_pointset(TEvePointSet* ps, Bool_t onFly)
{
  AliESDEvent* esd = AliEveEventManager::GetMaster()->AssertESD();

  Int_t NV0s = esd->GetNumberOfV0s();

  Double_t x, y, z;
  for (Int_t n = 0; n < NV0s; ++n)
  {
    AliESDv0* av = esd->GetV0(n);
    if (av->GetOnFlyStatus() == onFly)
    {
      av->GetXYZ(x, y, z);
      ps->SetNextPoint(x, y, z);
      ps->SetPointId(av);
    }
  }
}

TEvePointSet* esd_V0_points_offline()
{
  TEvePointSet* points = new TEvePointSet("V0 offline vertex locations");

  esd_VO_fill_pointset(points, kFALSE);

  points->SetTitle(Form("N=%d", points->Size()));
  points->SetMarkerStyle(4);
  points->SetMarkerSize(1.5);
  points->SetMarkerColor(kOrange+8);

  gEve->AddElement(points);
  gEve->Redraw3D();

  return points;
}

TEvePointSet* esd_V0_points_onfly()
{
  TEvePointSet* points = new TEvePointSet("V0 on-the-fly vertex locations");

  esd_VO_fill_pointset(points, kTRUE);

  points->SetTitle(Form("N=%d", points->Size()));
  points->SetMarkerStyle(4);
  points->SetMarkerSize(1.5);
  points->SetMarkerColor(kPink+10);

  gEve->AddElement(points);
  gEve->Redraw3D();

  return points;
}


void esd_V0_points()
{
  esd_V0_points_offline();
  esd_V0_points_onfly();
}
