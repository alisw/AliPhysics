void testPF(char *rootfile,int patch)
{
  gStyle->SetOptStat(0);
  
  Int_t xbin = 70;
  Int_t ybin = 70;
  Float_t xrange[2] = {-0.006 , 0.006}; //Pt 0.2->
  //Float_t yrange[2] = {-0.17 , 0.17}; //slice 2 0.55->0.88
  Float_t yrange[2] = {-0.26 , 0.26};// -15 - 15
  //Float_t yrange[2] = {0.55 , 0.88}; //slice 2 0.55->0.88
  
  TH2F *hist = new TH2F("hist","Parameter space",xbin,xrange[0],xrange[1],ybin,yrange[0],yrange[1]);
  hist->GetXaxis()->SetTitle("#kappa [cm^{-1}]");
  hist->GetYaxis()->SetTitle("#Phi [rad]");
  SetTH1Options(hist);

  Int_t xr[2] = {0,250};
  Int_t yr[2] = {-125,125};
 
  TH2F *raw = new TH2F("raw","",250,xr[0],xr[1],250,yr[0],yr[1]);
  TH2F *road = new TH2F("road","",250,0,250,250,yr[0],yr[1]);
  TH2F *peaks = new TH2F("peaks","",xbin,xrange[0],xrange[1],ybin,yrange[0],yrange[1]);
  peaks->SetMarkerStyle(3);
  peaks->SetMarkerColor(2);
  road->SetMarkerStyle(6);
  road->SetMarkerColor(2);
  
  real_peaks = new TH2F("real_peaks","",xbin,xrange[0],xrange[1],ybin,yrange[0],yrange[1]);
  real_peaks->SetMarkerStyle(3);
  real_peaks->SetMarkerColor(4);

  int slice = 2,n_phi_segments=30;
  //double eta[2] = {0.3,0.4};
  double eta[2] = {0.04,0.05};
  
  a = new AliL3HoughTransformer(slice,0,eta,1);
  a->GetPixels(rootfile,raw);
  a->InitTemplates(hist);
  a->Transform2Circle(hist,3);
  
  b = new AliL3HoughMaxFinder("KappaPhi");
  c = new AliL3HoughEval(a);
  //  tracks = (AliL3TrackArray*)b->FindMaxima(hist);
  /*
  int n_iter=1;
  AliL3TrackArray **track_cand = new AliL3TrackArray*[n_iter];
  for(int k=0; k<n_iter; k++)
    {
      hist->Reset();
      a->Transform2Circle(hist);
      track_cand[k] = (AliL3TrackArray*)b->FindPeak(hist,3,0.95,5);
      c->LookInsideRoad(track_cand[k],true);
    }
  
  //tracks = (AliL3TrackArray*)b->LookInWindows(hist,4,4,0.95,5);
  
  //tracks = (AliL3TrackArray*)b->FindMaxima(hist);
  //cout << "Found "<<tracks->GetNTracks()<<" peaks"<<endl;
  
  for(int j=0; j<n_iter; j++)
    {
      track = (AliL3HoughTrack*)track_cand[j]->GetCheckedTrack(0);
      if(!track) continue;
      printf("Pt %f phi0 %f weight %d\n",track->GetPt(),track->GetPhi0(),track->GetWeight());
      for(int padrow=0; padrow<174; padrow++)
	{
	  Float_t xyz[3];
	  track->GetCrossingPoint(padrow,xyz);
	  road->Fill(xyz[0],xyz[1],1);
	}
      peaks->Fill(track->GetKappa(),track->GetPhi0(),1);
    }
  */
  
  printf("Looking for big maxima\n");
  tracks = (AliL3TrackArray*)b->FindPeak(hist,3,0.95,5);
  cout<<"NTracks after checking "<<tracks->GetNTracks()<<endl;
  for(int i=0; i<tracks->GetNTracks(); i++)
    {
      track = (AliL3HoughTrack*)tracks->GetCheckedTrack(i);
      peaks->Fill(track->GetKappa(),track->GetPhi0(),1);
      if(!track) continue;
      printf("Pt %f phi0 %f weight %d\n",track->GetPt(),track->GetPhi0(),track->GetWeight());
    }  
  
  
  c2 = new TCanvas("c2","",900,900);
  c2->Divide(2,2);
  c2->cd(1);
  raw->Draw("");
  road->Draw("same");
  c2->cd(2);
  raw->DrawCopy();
  c2->cd(3);
  hist->Draw("box");
  peaks->Draw("same");
  //real_peaks->Draw("same");
  c2->cd(4);
  hist->DrawCopy("lego1");

  delete a;
  delete b;
  delete c;

}
