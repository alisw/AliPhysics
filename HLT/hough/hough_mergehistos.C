void hough_mergehistos(char *rootfile)
{
  gStyle->SetOptStat(0);
  
  Int_t xbin = 60;
  Int_t ybin = 60;
  Float_t xrange[2] = {-0.006 , 0.006}; //Pt 0.1->4GeV 0.006-0.006
  //Float_t yrange[2] = {-0.17 , 0.17}; //slice 2 0.55->0.88
  Float_t yrange[2] = {-0.26 , 0.26}; //slice 2 0.55->0.88
  //Float_t yrange[2] = {0.55 , 0.88}; //slice 2 0.55->0.88

  Int_t xr[2] = {0,250};
  Int_t yr[2] = {-125,125};
  TH1F *ntracks = new TH1F("ntracks","",100,0,200);
  TH2F *raw = new TH2F("raw","",250,xr[0],xr[1],250,yr[0],yr[1]);
  TH2F *road = new TH2F("road","",250,0,250,250,yr[0],yr[1]);
  TH2F *fake = new TH2F("fake","",300,0,300,250,-125,125);
  TH2F *peaks = new TH2F("peaks","",xbin,xrange[0],xrange[1],ybin,yrange[0],yrange[1]);
  peaks->SetMarkerStyle(3);
  peaks->SetMarkerColor(2);
  
  int slice = 2,patch=0,n_phi_segments=40;
  float eta[2] = {0.3,0.4};
  // float eta[2] = {0.04,0.05};
  
  AliL3HoughTransformer *a;
  AliL3HoughMaxFinder *b = new AliL3HoughMaxFinder("KappaPsi");
  
  TH2F **hist = new TH2F*[5];
  
  for(int pat=0; pat<5; pat++)
    {
      a = new AliL3HoughTransformer(slice,pat,eta,n_phi_segments);
      hist[pat] = new TH2F("hist","Parameter space",xbin,xrange[0],xrange[1],ybin,yrange[0],yrange[1]);
      a->GetPixels(rootfile,raw);
      a->InitTemplates(hist[pat]);
      a->Transform2Circle(hist[pat]);
      
      //AliL3TrackArray *tracks = b->FindMaxima(hist[pat]);
      
      delete a;
      
    }
  
  
  TCanvas *c1 = new TCanvas("c1","",800,800);
  c1->Divide(2,2);
  c1->cd(1);
  hist[0]->DrawCopy("box");
  c1->cd(2);
  hist[4]->DrawCopy("box");
  c1->cd(3);
  hist[2]->DrawCopy("box");
  
  c1->cd(4);
  hist[4]->Add(hist[0],hist[2]);
  hist[4]->Draw("box");

  delete b;
  
  
  delete [] hist;
}
