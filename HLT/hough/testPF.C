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
  
  real_peaks = new TH2F("real_peaks","",xbin,xrange[0],xrange[1],ybin,yrange[0],yrange[1]);
  real_peaks->SetMarkerStyle(3);
  real_peaks->SetMarkerColor(4);

  int slice = 2,n_phi_segments=30;
  float eta[2] = {0.3,0.4};
  //float eta[2] = {0.03,0.04};
  
  a = new AliL3HoughTransformer(slice,patch,eta,n_phi_segments);
  a->GetPixels(rootfile,raw);
  a->InitTemplates(hist);
  a->Transform2Circle(hist);
    
  b = new AliL3HoughMaxFinder("KappaPhi");
  tracks = (AliL3TrackArray*)b->FindPeak(hist,3,0.95,5);

  c = new AliL3HoughEval(a);
  c->LookInsideRoad(tracks);

  for(int i=0; i<tracks->GetNTracks(); i++)
    {
      track = (AliL3HoughTrack*)tracks->GetCheckedTrack(i);
      peaks->Fill(track->GetKappa(),track->GetPhi0(),1);
      if(!track) continue;
      printf("Pt %f phi0 %f weight %d\n",track->GetPt(),track->GetPhi0(),track->GetWeight());
    }  
  
  c2 = new TCanvas("c2","",800,800);
  //  c2->Divide(2,2);
  //c2->cd(1);
  hist->Draw("box");
  peaks->Draw("same");
  real_peaks->Draw("same");
  
  
  /*
  AliL3TrackArray *tracks = b->FindMaxima(hist);
  
  track = (AliL3HoughTrack*)tracks->GetCheckedTrack(0);
  peaks->Fill(track->GetKappa(),track->GetPhi0(),track->GetWeight());
  int xbin = hist->GetXaxis()->FindBin(track->GetKappa());
  int ybin = hist->GetYaxis()->FindBin(track->GetPhi0());
  
  TH1D *proj_kappa = hist->ProjectionX("proj_kappa",ybin,ybin);
  TH1D *proj_phi = hist->ProjectionY("proj_phi",xbin,xbin);
  
  TCanvas *c1 = new TCanvas("c1","",800,800);
  c1->Divide(2,2);
  c1->cd(1);
  hist->Draw("box");
  peaks->Draw("same");
  
  c1->cd(3);
  proj_kappa->Draw();
  c1->cd(4);
  proj_phi->Draw();
  return;
  */
  delete a;
  delete b;
  delete c;

}
