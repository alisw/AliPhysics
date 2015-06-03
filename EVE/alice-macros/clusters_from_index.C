// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/
#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TEveManager.h>
#include <TEveElement.h>
#include <TEvePointSet.h>

#include <AliESDEvent.h>
#include <AliTrackPointArray.h>
#include <AliEveEventManager.h>
#include <AliEveMultiView.h>
#endif

TEvePointSet* clusters_from_index(Int_t index=0, TEveElement* cont=0)
{
  AliESDEvent* esd = AliEveEventManager::GetMaster()->AssertESD();

  if (index < 0) {
    Warning("clusters_from_index", "index not set.");
    return 0;
  }

  if (index >= esd->GetNumberOfTracks()) {
    Warning("clusters_from_index", "index out of range");
    return 0;
  }

  TEvePointSet* clusters = new TEvePointSet(64);
  clusters->SetOwnIds(kTRUE);

  AliESDtrack* at = esd->GetTrack(index);
  const AliTrackPointArray* pArr = at->GetTrackPointArray();
  if (pArr == 0) {
    Warning("clusters_from_index", "TrackPointArray not stored with ESD track.");
  }
  
  Int_t np =  pArr->GetNPoints();
  const Float_t* x = pArr->GetX();
  const Float_t* y = pArr->GetY();
  const Float_t* z = pArr->GetZ();
  for (Int_t i=0; i<np; ++i) {
    clusters->SetNextPoint(x[i], y[i], z[i]);
    AliTrackPoint *atp = new AliTrackPoint;
    pArr->GetPoint(*atp, i);
    clusters->SetPointId(atp);    }


  if(clusters->Size() == 0 && gEve->GetKeepEmptyCont() == kFALSE) {
    Warning("clusters_from_index", "No clusters for index '%d'", index);
    delete clusters;
    return 0;
  }

  clusters->SetMarkerStyle(2);
  clusters->SetMarkerSize(2);
  clusters->SetMarkerColor(4);

  clusters->SetName(Form("Clusters idx=%d", index));
  clusters->SetTitle(Form("N=%d", clusters->Size()));

  gEve->AddElement(clusters);

  if (AliEveMultiView::Instance())
  {
    AliEveMultiView::Instance()->ImportEventRPhi(clusters);
    AliEveMultiView::Instance()->ImportEventRhoZ(clusters);
  }

  gEve->Redraw3D();

  return clusters;
}
