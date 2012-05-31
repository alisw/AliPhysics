// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/
#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TPolyMarker3D.h>
#include <TEveManager.h>
#include <TEveElement.h>
#include <TEveUtil.h>
#endif


void region_marker(Float_t a=10, Float_t b=10, Float_t c=20,
		   Float_t x=0, Float_t y=0, Float_t z=0)
{
  TPolyMarker3D* mark = new TPolyMarker3D(8);
  mark->SetName("Origin marker");
  mark->SetMarkerStyle(3);
  mark->SetMarkerColor(6);
  mark->SetPoint(0, x+a, y+b, z+c);
  mark->SetPoint(1, x+a, y-b, z+c);
  mark->SetPoint(2, x-a, y-b, z+c);
  mark->SetPoint(3, x-a, y+b, z+c);

  mark->SetPoint(4, x+a, y+b, z-c);
  mark->SetPoint(5, x+a, y-b, z-c);
  mark->SetPoint(6, x-a, y+b, z-c);
  mark->SetPoint(7, x-a, y-b, z-c);
  Color_t* colp = TEveUtil::FindColorVar(mark, "fMarkerColor");
  TEveElementObjectPtr* rnrEl = new TEveElementObjectPtr(mark, *colp);
  gEve->AddGlobalElement(rnrEl);
  gEve->Redraw3D();
}
