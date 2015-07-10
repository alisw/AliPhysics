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
#include <AliESDtrack.h>
#include <AliTrackPointArray.h>
#include <AliEveEventManager.h>
#include <AliEveMultiView.h>
#endif

TEvePointSet* clusters_from_label(Int_t label=0, TEveElement* cont=0)
{
  AliESDEvent* esd = AliEveEventManager::GetMaster()->AssertESD();
  TEvePointSet* clusters = new TEvePointSet(64);
  clusters->SetOwnIds(kTRUE);

  for (Int_t n=0; n<esd->GetNumberOfTracks(); n++)
  {
    AliESDtrack* at = esd->GetTrack(n);
    if (at->GetLabel() == label) {
      const AliTrackPointArray* pArr = at->GetTrackPointArray();
      if (pArr == 0) {
	Warning("clusters_from_label", "TrackPointArray not stored with ESD track.");
	continue;
      }
      Int_t np =  pArr->GetNPoints();
      const Float_t* x = pArr->GetX();
      const Float_t* y = pArr->GetY();
      const Float_t* z = pArr->GetZ();
      for (Int_t i=0; i<np; ++i) {
	clusters->SetNextPoint(x[i], y[i], z[i]);
	AliTrackPoint *atp = new AliTrackPoint;
	pArr->GetPoint(*atp, i);
	clusters->SetPointId(atp);
      }
    }
  }

  if(clusters->Size() == 0 && gEve->GetKeepEmptyCont() == kFALSE) {
    Warning("clusters_from_label", "No clusters match label '%d'", label);
    delete clusters;
    return 0;
  }

  clusters->SetMarkerStyle(2);
  clusters->SetMarkerSize(0.5);
  clusters->SetMarkerColor(4);
  //PH The line below is replaced waiting for a fix in Root
  //PH which permits to use variable siza arguments in CINT
  //PH on some platforms (alphalinuxgcc, solariscc5, etc.)
  //PH  clusters->SetName(Form("Clusters lab=%d", label));
  char form[1000];
  sprintf(form,"Clusters lab=%d", label);
  clusters->SetName(form);

  char tip[1000];
  sprintf(tip,"N=%d", clusters->Size());
  clusters->SetTitle(tip);
  gEve->AddElement(clusters, cont);
  gEve->Redraw3D();

  return clusters;
}
