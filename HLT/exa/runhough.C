void runhough(Int_t slice,Char_t *path,Int_t n_eta_segments)
{
  

  hough = new AliL3Hough();
  hough->Init(path,kTRUE,n_eta_segments,kTRUE);
  
  hough->ReadData(slice);

  hough->Transform();

  hough->SetPeakThreshold(1);
  hough->AddAllHistograms();
  hough->FindTrackCandidates();
  
  hough->WriteTracks(slice);
  
  //hough->Evaluate();
  tracks = (AliL3TrackArray*)hough->GetTracks(0);
  for(int i=0; i<tracks->GetNTracks(); i++)
    {
      track = (AliL3HoughTrack*)tracks->GetCheckedTrack(i);
      if(!track) continue;
      cout<<"pt "<<track->GetPt()<<" psi "<<track->GetPsi()<<" etaindex "<<track->GetEtaIndex()<<" weight "<<track->GetWeight()<<endl;
    }
  

  display(hough,0);
  
}

void display(AliL3Hough *hough,Int_t eta_index)
{
  //Display the data/tracks in eta_index
  
  hough->InitEvaluate();
  digitd = new AliL3Histogram("Digits display","",250,0,250,250,-125,125);
  trackd = new AliL3Histogram("Found tracks display","",250,0,250,250,-125,125);
  for(int i=0; i<6; i++)
    hough->GetEval(i)->DisplayEtaSlice(eta_index,digitd);
  
  tracks = (AliL3TrackArray*)hough->GetTracks(0);
  float xyz[3];
  for(int i=0; i<tracks->GetNTracks(); i++)
    {
      AliL3HoughTrack *track = (AliL3HoughTrack*)tracks->GetCheckedTrack(i);
      if(!track) continue;
      if(track->GetEtaIndex() != eta_index) continue;

      for(int j=0; j<176; j++)
	{
	  track->GetCrossingPoint(j,xyz);
	  trackd->Fill(xyz[0],xyz[1],1);
	}
    }
  
  //Draw the parameter space
  c1 = new TCanvas("c1","",2);
  hough->GetTransformer(0)->GetHistogram(eta_index)->Draw("lego");
  return;
  //Draw the tracks
  c2 = new TCanvas("c2","",2);
  digitd->Draw();
  trackd->Draw("same");
  ((TH1F*)trackd->GetRootHisto())->SetMarkerColor(2);
}
