//
// Select points from TEvePointSet(clusters, hits, etc.) with Alt+mouse-left
// click action.
//
// In AliEvetrackFitEditor press "Fit" button to make track fit on the
// selected points. To fit new track, press "Reset".

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
