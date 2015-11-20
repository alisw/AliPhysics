/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

//
// Select points from TEvePointSet(clusters, hits, etc.) with Alt+mouse-left
// click action.
//
// In AliEvetrackFitEditor press "Fit" button to make track fit on the
// selected points. To fit new track, press "Reset".

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TEveManager.h>
#include <TEveSelection.h>

#include <AliEveTrackFitter.h>
#include <AliEveCosmicRayFitter.h>
#endif

void track_fitter(Int_t mode = 1)
{
  gEve->GetSelection()->SetPickToSelect(1);
  gEve->GetHighlight()->SetPickToSelect(0);

  if (mode == 0)
  {
    // helix fit
    AliEveTrackFitter* t = new AliEveTrackFitter();
    gEve->AddElement(t);
    t->Start();
  }
  else
  {
    // straight line fit
    AliEveCosmicRayFitter* t = new AliEveCosmicRayFitter();
    gEve->AddElement(t);
    t->Start();
  }
}
