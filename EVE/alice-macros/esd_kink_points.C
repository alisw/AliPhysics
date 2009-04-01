// $Id$
// Main authors: Paraskevi Ganoti: 2009

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

void esd_kink_fill_pointset(TEvePointSet* ps)
{
  AliESDEvent* esd = AliEveEventManager::AssertESD();

  for (Int_t n=0; n<esd->GetNumberOfTracks(); ++n)
  { 
      AliESDtrack* track = esd->GetTrack(n);
      if(track->GetKinkIndex(0)<0){
    
          AliESDkink *kink = esd->GetKink(TMath::Abs(track->GetKinkIndex(0))-1);
	  const TVector3 Position(kink->GetPosition());
	  ps->SetNextPoint(Position.X(), Position.Y(), Position.Z());
          ps->SetPointId(kink);
      }
  }

}

TEvePointSet* esd_kink_points()
{
  TEvePointSet* points = new TEvePointSet("Kink vertex locations");

  esd_kink_fill_pointset(points);

  points->SetTitle(Form("N=%d", points->Size()));
  points->SetMarkerStyle(4);
  points->SetMarkerSize(1.5);
  points->SetMarkerColor(kOrange+8);

  gEve->AddElement(points);
  gEve->Redraw3D();

  return points;
}
