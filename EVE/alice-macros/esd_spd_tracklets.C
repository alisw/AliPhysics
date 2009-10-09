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

TEveElementList* esd_spd_tracklets(Float_t radius=8, Width_t line_width=3)
//x TEveTrackList* esd_spd_tracklets(Float_t rad=8)
{
  // radius - cylindrical radius to which the tracklets should be extrapolated

  AliESDEvent     *esd = AliEveEventManager::AssertESD();
  AliESDVertex    *pv  = esd->GetPrimaryVertexSPD();
  AliMultiplicity *mul = esd->GetMultiplicity();

  Double_t pvx[3], pve[3];
  pv->GetXYZ(pvx);
  pv->GetSigmaXYZ(pve);

  TEveCompound *cont = new TEveCompound("SPD Tracklets");
  cont->OpenCompound();
  //x TEveTrackList *cont = new TEveTrackList("SPD Tracklets");
  cont->SetTitle(Form("N=%d", mul->GetNumberOfTracklets()));
  cont->SetMainColor(7);
  //x cont->SetLineWidth(line_width);

  gEve->AddElement(cont);

  for (Int_t i=0; i<mul->GetNumberOfTracklets(); ++i)
  {
    using namespace TMath;
    Float_t theta = mul->GetTheta(i);
    Float_t phi   = mul->GetPhi(i);
    Float_t dr[3];
    dr[0] = radius*Cos(phi);
    dr[1] = radius*Sin(phi);
    dr[2] = radius/Tan(theta);

    TEveLine* track = new TEveLine;
    track->SetMainColor(7);
    track->SetLineWidth(line_width);
    //x AliEveTrack* track = new AliEveTrack;
    //x track->SetPropagator(cont->GetPropagator());
    //x track->SetAttLineAttMarker(cont);
    track->SetElementName(Form("Tracklet %d", i));
    track->SetElementTitle(Form("id=%d: theta=%.3f, phi=%.3f", i, theta, phi));

    track->SetPoint(0, pvx[0], pvx[1], pvx[2]);
    track->SetPoint(1,pvx[0]+dr[0], pvx[1]+dr[1], pvx[2]+dr[2]);

    //x track->SetLockPoints(kTRUE);

    cont->AddElement(track);
  }

  gEve->Redraw3D();

  return cont;
}
