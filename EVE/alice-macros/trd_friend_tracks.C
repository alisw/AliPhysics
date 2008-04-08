// $Id: esd_tracks.C 24485 2008-03-13 15:27:38Z mtadel $
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

TEvePointSet* trd_friend_tracks(TEveElement *cont = 0)
{
  if(!AliEveEventManager::AssertESDfriend())
  {
    printf("AliESDfriend not found\n");
    return 0;
  }
  AliESDEvent* esd = AliEveEventManager::AssertESD();
  AliESDfriend *eventESDfriend = AliEveEventManager::AssertESDfriend();
  AliEveEventManager::AssertGeometry();
  printf("eventESDfriend = %p\n", eventESDfriend);
  //	esd->SetESDfriend(eventESDfriend);
  //	if(!eventESDfriend) return 0;

  Int_t nTracks = esd->GetNumberOfTracks();
  printf("Number of tracks: %d\n", nTracks);
  Int_t kMaxClusters = nTracks * 6 * 40 * 3;
  printf("kMaxClusters = %d\n", kMaxClusters);
  TEvePointSet *friendClusters = new TEvePointSet(kMaxClusters);

  Int_t kTotalClusters = 0;
  for (Int_t n=0; n< nTracks; n++)
  {
    printf("Taking track %d\n", n);
    AliESDtrack* at = esd->GetTrack(n);
    AliESDfriendTrack *friendTrack = eventESDfriend->GetTrack(n);

    TObject *cal = 0x0;
    Int_t ical = 0;
    while (cal = friendTrack->GetCalibObject(ical++))
    {
      if (strcmp(cal->IsA()->GetName(), "AliTRDtrackV1") != 0) continue;
      printf("Got the track\n");
      AliTRDtrackV1 *trdTrack = dynamic_cast<AliTRDtrackV1 *>(cal);
      printf("trdTrack = %p\n", trdTrack);

      for (Int_t iseed = 0; iseed < 6; iseed++)
      {
        printf("processing tracklet %d\n", iseed);
        AliTRDseedV1 *tracklet = trdTrack->GetTracklet(iseed);
        if (!tracklet->IsOK()) continue;
        for (Int_t itb = 0; itb < 30; itb++)
        {
          if (!tracklet->IsUsable(itb)) continue;
          AliTRDcluster *cl = tracklet->GetClusters(itb);
          if (!cl) continue;
          printf("cluster position: x = %f, y = %f, z = %f\n", cl->GetX(), cl->GetY(), cl->GetZ());
          Float_t globalCoords[3]; //global coordinates
          cl->GetGlobalXYZ(globalCoords);
          friendClusters->SetPoint(kTotalClusters, globalCoords[0], globalCoords[1], globalCoords[2]);
          AliCluster *atp = new AliCluster(*cl);
          friendClusters->SetPointId(atp);
          kTotalClusters++;
        }
      }
    }
  }

  printf("kTotalClusters = %d\n", kTotalClusters);

  if (friendClusters->Size() == 0 && gEve->GetKeepEmptyCont() == kFALSE)
  {
    Warning("trd_clusters", "No TRD clusters");
    delete friendClusters;
    return 0;
  }

  friendClusters->SetMarkerStyle(2);
  friendClusters->SetMarkerSize(0.2);
  friendClusters->SetMarkerColor(4);

  char form[1000];
  sprintf(form,"TRD Clusters");
  friendClusters->SetName(form);

  char tip[1000];
  sprintf(tip,"N=%d", friendClusters->Size());
  friendClusters->SetTitle(tip);
  gEve->AddElement(friendClusters, cont);
  gEve->Redraw3D();
  return friendClusters;
}
