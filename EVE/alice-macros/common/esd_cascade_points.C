// $Id$

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/
#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TEveManager.h>
#include <TEvePointSet.h>

#include <AliESDEvent.h>
#include <AliESDcascade.h>
#include <AliRunLoader.h>
#include <AliEveEventManager.h>
#endif

void esd_cascade_fill_pointset(TEvePointSet* ps)
{
  AliESDEvent* esd = AliEveEventManager::GetMaster()->AssertESD();

  Int_t NCascades = esd->GetNumberOfCascades();

  Double_t x, y, z;
  for (Int_t n = 0; n < NCascades; ++n)
  {
    AliESDcascade* av = esd->GetCascade(n);
    av->GetXYZcascade(x, y, z);
    ps->SetNextPoint(x, y, z);
    ps->SetPointId(av);
  }
}

TEvePointSet* esd_cascade_points()
{
  TEvePointSet* points = new TEvePointSet("Cascade vertex locations");

  esd_cascade_fill_pointset(points);

  points->SetTitle(Form("N=%d", points->Size()));
  points->SetMarkerStyle(4);
  points->SetMarkerSize(1.5);
  points->SetMarkerColor(kMagenta-9);

  gEve->AddElement(points);
  gEve->Redraw3D();

  return points;
}
