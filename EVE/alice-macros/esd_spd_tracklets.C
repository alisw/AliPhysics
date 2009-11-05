// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

// Lines commented with //x should be reactivated when we
// move to a newer ROOT (after 5.25.2).
// Corresponding lines that are now replacing them should be removed.
//
// The problem was the TEveTrackProjected did not support projection
// of tracks with locked points -- and we do that for tracklets.
//
// Maybe it would be even better to have AliEveTracklet : public TEveLine
// and support both in AliEveTrackCounter.
// Or have trackelt counter -- as not all histograms collected for tracks
// are relevant for tracklets.

TEveTrackList* esd_spd_tracklets(Float_t radius=8, Width_t line_width=3,
				 Float_t d_theta=0.025, Float_t d_phi=0.08)
{
  // radius - cylindrical radius to which the tracklets should be extrapolated

  AliESDEvent     *esd = AliEveEventManager::AssertESD();
  AliESDVertex    *pv  = esd->GetPrimaryVertexSPD();
  AliMultiplicity *mul = esd->GetMultiplicity();

  AliMagF *field = AliEveEventManager::AssertMagField();

  TEveTrackList *cont = new TEveTrackList("SPD Tracklets");
  cont->SetTitle(Form("N=%d", mul->GetNumberOfTracklets()));
  cont->SetMainColor(7);
  cont->SetLineWidth(line_width);

  TEveTrackPropagator* prop = cont->GetPropagator();
  prop->SetMaxR(radius);
  gEve->AddElement(cont);

  for (Int_t i = 0; i < mul->GetNumberOfTracklets(); ++i)
  {
    Float_t theta = mul->GetTheta(i);
    Float_t phi   = mul->GetPhi(i);

    AliEveTracklet* t = new AliEveTracklet(i, pv, theta, phi, prop);
    t->SetAttLineAttMarker(cont);
    t->SetElementName(Form("Tracklet %d", i));
    t->SetElementTitle(Form("id=%d: eta=%.3f, theta=%.3f, phi=%.3f",
			    i, mul->GetEta(i), theta, phi));
    // if some condition
    mul->SetLabel(i, 0, 3);
    // else mul->SetLabel(i, 0, 0);

    cont->AddElement(t);
  }

  cont->MakeTracks();

  gEve->Redraw3D();

  return cont;
}
