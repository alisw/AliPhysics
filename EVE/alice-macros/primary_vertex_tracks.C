// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

TEveTrackList* primary_vertex_tracks()
{
  TEveUtil::LoadMacro("esd_tracks.C");
  AliESDEvent   *esd = AliEveEventManager::AssertESD();
  AliESDVertex *pv  = esd->GetPrimaryVertex();

  TEveTrackList* cont = new TEveTrackList("Tracks for Primary Vertex");
  cont->SetMainColor(7);
  TEveTrackPropagator* rnrStyle = cont->GetPropagator();
  rnrStyle->SetMagField( 0.1*esd->GetMagneticField() );
  rnrStyle->fRnrFV = kTRUE;
  rnrStyle->fFVAtt->SetMarkerColor(2);
  gEve->AddElement(cont);

  for (Int_t n=0; n<pv->GetNIndices(); n++)
  {
    AliESDtrack* at = esd->GetTrack(pv->GetIndices()[n]);
    AliEveTrack* track = esd_make_track(rnrStyle, n, at, at);
    track->SetLineWidth(4);
    track->SetLineColor(cont->GetMainColor());
    track->SetLineStyle(7);
    gEve->AddElement(track, cont);
  }

  cont->MakeTracks();
  gEve->Redraw3D();

  return cont;
}
