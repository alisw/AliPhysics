TEveElementList* trd_friend_tracks(TEveElement *cont = 0)
{

  // Link data containers
  AliESDfriend *eventESDfriend = 0x0;
  if(!(eventESDfriend = AliEveEventManager::AssertESDfriend())){
    Warning("trd_friend_tracks", "AliESDfriend not found");
    return 0x0;
  }

  AliESDEvent* esd = AliEveEventManager::AssertESD();

  AliEveEventManager::AssertGeometry();

  // TRD related objects
  AliTRDgeometry geo;

  Int_t nTracks = esd->GetNumberOfTracks();
  TEveElementList *tracks = new TEveElementList("TRD Tracks");
  TEveElementList *track = 0x0;
  TEveLine *tracklet = 0x0;
  TEvePointSet *clusters = 0x0;
  for (Int_t n=0; n<nTracks; n++){
    AliESDtrack* esdTrack = esd->GetTrack(n);
    AliESDfriendTrack *friendTrack = eventESDfriend->GetTrack(n);
  
    TObject *cal = 0x0;
    Int_t ical = 0;
    while(cal = friendTrack->GetCalibObject(ical++)){
      if(strcmp(cal->IsA()->GetName(), "AliTRDtrackV1") != 0) continue;
      AliTRDtrackV1 *trackObj = dynamic_cast<AliTRDtrackV1 *>(cal);
      
      // define track
      Int_t ncl = 0;
      track = new TEveElementList();
      track->SetName(Form("Track%d", esdTrack->GetLabel()));
//       track->SetLineColor(kWhite);
//       track->SetOwnIds(kTRUE);

      for(Int_t iseed = 0; iseed < 6; iseed++){
        AliTRDseedV1 *trackletObj = trackObj->GetTracklet(iseed);
        if(!trackletObj->IsOK()) continue;
       
        // define clusters
        clusters = new TEvePointSet();
        clusters->SetOwnIds(kTRUE);
        clusters->SetMarkerStyle(4);
        clusters->SetMarkerSize(0.5);
        clusters->SetMarkerColor(kGreen);
        clusters->SetName("Clusters");

        AliTRDcluster *cl = 0x0; Int_t sec, det = -1;
        for(Int_t itb = 0; itb < 30; itb++){
          if(!trackletObj->IsUsable(itb)) continue;
          if(!(cl = trackletObj->GetClusters(itb))) continue;

          if(det<0){
            det = cl->GetDetector();
            sec = geo.GetSector(det);
          }
          Float_t globalCoords[3]; //global coordinates
          cl->GetGlobalXYZ(globalCoords);
          Int_t id = clusters->SetNextPoint( globalCoords[0], globalCoords[1], globalCoords[2]);
          clusters->SetPointId(id, new AliTRDcluster(*cl));
        }
        ncl += clusters->Size();
        clusters->SetTitle(Form("Clusters %d", clusters->Size()));

        // define tracklet
        tracklet = new TEveLine(Form("Tracklet%d", iseed)); 
        tracklet->SetTitle(Form("P = %7.3f [GeV/c]", trackletObj->GetMomentum()));        tracklet->SetLineColor(kYellow);
        tracklet->SetOwnIds(kTRUE);
        Double_t alpha = AliTRDgeometry::GetAlpha() * (sec<9 ? sec + .5 : sec - 17.5); 
        Double_t x0 = trackletObj->GetX0(), 
          y0f = trackletObj->GetYfit(0), 
          ysf = trackletObj->GetYfit(1),
          z0r = trackletObj->GetZref(0), 
          zsr = trackletObj->GetZref(1);
        Double_t xg =  x0 * TMath::Cos(alpha) - y0f * TMath::Sin(alpha); 
        Double_t yg = x0 * TMath::Sin(alpha) + y0f * TMath::Cos(alpha);
        tracklet->SetPoint(0, xg, yg, z0r);
        tracklet->SetPointId(0, new AliTRDseedV1(*trackletObj));
        Double_t x1 = x0-3.5, 
          y1f = y0f - ysf*3.5,
          z1r = z0r - zsr*3.5; 
        xg =  x1 * TMath::Cos(alpha) - y1f * TMath::Sin(alpha); 
        yg = x1 * TMath::Sin(alpha) + y1f * TMath::Cos(alpha);
        tracklet->SetPoint(1, xg, yg, z1r);


        tracklet->AddElement(clusters);
        track->AddElement(tracklet);
      }
      track->SetTitle(Form("Clusters %d", ncl));       tracks->AddElement(track);
    }
  }
	
  tracks->SetTitle(Form("Tracks %d", tracks->GetNChildren()));       	
  gEve->AddElement(tracks, cont);

  gEve->Redraw3D();
  
  return tracks;
}
