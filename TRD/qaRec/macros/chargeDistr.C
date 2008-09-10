void chargeDistr(AliTRDtrackV1* track, Double_t* &results, Int_t& nResults)
{
  if (track == 0)  return;
  
  Int_t Nt = track->AliKalmanTrack::GetNumberOfTracklets();;
  AliTRDcluster* cls = 0;
  
  // Count cluster
  nResults = 0;
  for (Int_t trackletInd = 0; trackletInd < Nt; trackletInd++)
  {
    if(!(tracklet = track->GetTracklet(trackletInd))) continue;
    if(!tracklet->IsOK()) continue;
    
    for (Int_t clusterInd = 0; clusterInd < 34; clusterInd++) 
    {
      if(!(cls = tracklet->GetClusters(clusterInd))) continue;
            
	    nResults++;
	  }
  }
  
  // Allocate memory for the results (AliEveTRDTrackList will clean this memory automatically)
  results = new Double_t[nResults];
  for (Int_t i = 0; i < nResults; i++)  results[i] = -100;
  Int_t currentIndex = 0;

  for (Int_t trackletInd = 0; trackletInd < Nt; trackletInd++)
  {
    if(!(tracklet = track->GetTracklet(trackletInd))) continue;
    if(!tracklet->IsOK()) continue;
    
    for (Int_t clusterInd = 0; clusterInd < 34; clusterInd++) 
    {
      if(!(cls = tracklet->GetClusters(clusterInd))) continue;
            
	    results[currentIndex++] = cls->GetQ();
	  }
  }
}

/*
Double_t chargeDistr(AliTRDtrackV1* track)
{
  Double_t returnValue = 0;
  AliTRDseedV1 *tracklet = 0x0;
  Int_t Nt = track->AliKalmanTrack::GetNumberOfTracklets();

  for (Int_t trackletInd = 0; trackletInd < Nt; trackletInd++)
  {
    if(!(tracklet = track->GetTracklet(trackletInd))) continue;
    if(!tracklet->IsOK()) continue;
    
    for (Int_t clusterInd = 0; clusterInd < 34; clusterInd++) 
    {
      if(!(cls = tracklet->GetClusters(clusterInd))) continue;
            
	    returnValue += cls->GetQ();
	  }
  }

  return returnValue;
}

void chargeDistr(AliEveTRDTrackList* trackList)
{
  if (trackList == 0)  return;
  
  // hist
  TH1D *mQ = new TH1D("Q", ";Charge", 200, 0, 200);
  
  Int_t Nt = 0;

  AliTRDcluster* cls = 0;
  AliEveTRDTrack* track = 0;

  for (TEveElement::List_i iterator = trackList->BeginChildren(); iterator != trackList->EndChildren(); ++iterator)
  {
    track = dynamic_cast<AliEveTRDTrack*>(*iterator);

    if (!track) continue;
    if (!track->GetRnrState()) continue;
      
    AliTRDtrackV1 *trackv1 = 0x0;
    AliTRDseedV1 *tracklet = 0x0;
    trackv1 = (AliTRDtrackV1*)track->GetUserData();

    Nt = trackv1->AliKalmanTrack::GetNumberOfTracklets();

    for (Int_t trackletInd = 0; trackletInd < Nt; trackletInd++)
    {
      if(!(tracklet = trackv1->GetTracklet(trackletInd))) continue;
      if(!tracklet->IsOK()) continue;
    
      for (Int_t clusterInd = 0; clusterInd < 34; clusterInd++) 
      {
        if(!(cls = tracklet->GetClusters(clusterInd))) continue;
            
	      mQ->Fill(cls->GetQ());
	    }
    }
  }
  
  gEve->AddCanvasTab("Charge distribution");
  mQ->Draw();
}
*/
