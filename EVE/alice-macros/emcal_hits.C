
/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

//************************************************************************
///
/// \file emcal_hits.C
/// \brief Visualize EMCAL digits
///
/// A macro to read and visualize EMCAL hits. Standalone.
///
/// \author Magali Estienne <magali.estienne@cern.ch>, SUBATECH. EMCal implementation, June 2008
//************************************************************************


#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TTree.h>
#include <TString.h>
#include <TEveManager.h>
#include <TEveElement.h>
#include <TEvePointSet.h>
#include <TEveTreeTools.h>

#include <AliRunLoader.h>
#include <AliEveEventManager.h>
#endif

TEvePointSet*
emcal_hits(const char *varexp    = "fX:fY:fZ",
           const char *selection = "",
           TEveElement* cont = 0)
{
  AliRunLoader* rl =  AliEveEventManager::AssertRunLoader();
  rl->LoadHits("EMCAL");

  TTree* ht = rl->GetTreeH("EMCAL", false);

  TEvePointSet* points = new TEvePointSet(Form("EMCAL Hits '%s'", selection));

  TEvePointSelector ps(ht, points, varexp, selection);
  ps.Select();

  if (points->Size() == 0 && gEve->GetKeepEmptyCont() == kFALSE) 
  {
    Warning("emcal_hits", "No hits match '%s'",selection);
    delete points;
    return 0;
  }

  points->SetName(Form("EMCAL Hits"));
  const TString viz_tag("SIM Hits EMCAL");
  points->ApplyVizTag(viz_tag, "Hits");

  points->SetTitle(Form("N=%d", points->Size()));
  points->SetMarkerSize(.5);
  points->SetMarkerColor(2);

  gEve->AddElement(points, cont);
  gEve->Redraw3D();

  return points;
}
