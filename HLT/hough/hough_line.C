void hough_line(char *rootfile,int patch,Bool_t MC=false)
{
  gStyle->SetOptStat(0);
  
  
  Double_t pi = TMath::Pi();
  Int_t xbin = 60;
  Int_t ybin = 60;
  Float_t xrange[2] = {-pi/2,pi/2};
  Float_t yrange[2] = {-40,40};
  
  TH2F *hist = new TH2F("hist","Parameter space",xbin,xrange[0],xrange[1],ybin,yrange[0],yrange[1]);
  hist->GetXaxis()->SetTitle("#Psi");
  hist->GetYaxis()->SetTitle("D");
  SetTH1Options(hist);

  Int_t xr[2] = {0,250};
  Int_t yr[2] = {-125,125};
  TH1F *ntracks = new TH1F("ntracks","",100,0,200);
  TH2F *raw = new TH2F("raw","",250,xr[0],xr[1],250,yr[0],yr[1]);
  TH2F *road = new TH2F("road","",250,0,250,250,yr[0],yr[1]);
  TH2F *container = new TH2F("container","",250,0,250,250,yr[0],yr[1]);
  TH2F *peaks = new TH2F("peaks","",xbin,xrange[0],xrange[1],ybin,yrange[0],yrange[1]);
  peaks->SetMarkerStyle(3);
  peaks->SetMarkerColor(2);
  
  int slice = 2,n_phi_segments=20;
  float eta[2] = {0.3,0.4};
  //float eta[2] = {0.05,0.06};
  
  AliL3HoughTransformer *a = new AliL3HoughTransformer(slice,patch,eta,n_phi_segments);

  const Int_t NRows[5][2] = {{ 0, 45},{46,77},{78,109},{110,141},{142,173}};
  const int ref_row[5]={22,60,93,125,157};
  int rowrange[2] = {6,12};
  double phirange[2] = {-10.,10.};
  a->GetPixels(rootfile,raw);
  a->Transform2Line(hist,9,rowrange,phirange);

  AliL3HoughMaxFinder *b = new AliL3HoughMaxFinder("DPsi");
  AliL3TrackArray *tracks = b->FindMaxima(hist,rowrange,5);

  TH2F *cross = new TH2F("cross","",250,xr[0],xr[1],250,yr[0],yr[1]);
  for(int i=0; i<2; i++)
  //for(int i=0; i<tracks->GetNTracks(); i++)
    {
      AliL3HoughTrack *track = (AliL3HoughTrack*)tracks->GetCheckedTrack(i);
      if(!track) continue;
      Double_t xy[2];
      for(int j=0; j<174; j++)
	{
	  track->GetLineCrossingPoint(j,xy);
	  cross->Fill(xy[0],xy[1],1);
	}
      printf("Segment %f %f weight %d\n",track->GetPsiLine(),track->GetDLine(),track->GetWeight());
    }
  cross->SetMarkerStyle(6);
  cross->SetMarkerColor(2);
  
  
  TCanvas *c1 = new TCanvas("c1","",900,900);
  SetCanvasOptions(c1);
  c1->Divide(2,2);
  c1->cd(1);
  hist->Draw("box");
  peaks->SetMarkerStyle(3);
  peaks->SetMarkerColor(2);
  //peaks->Draw("same");

  c1->cd(2);
  container->Draw();

  c1->cd(3);
  raw->Draw("");
  cross->Draw("same");
  return;
  //TCanvas *c3 = new TCanvas("c3","",2);
  //road->SetMarkerStyle(5);
  //road->SetMarkerColor(3);
  c1->cd(4);
  road->SetMarkerStyle(9);
  road->Draw();
  /*
  TCanvas *c4 = new TCanvas("c4","",2);
  fake->Draw("cont1");

  TCanvas *c5 = new TCanvas("c5","",2);
  ntracks->Draw();
  */
  delete a;
  delete b;
  delete c;

  //----------------------------------------------------------------------------
  if(MC)
    {
      TFile *file = new TFile(rootfile);
      file->cd();
      
      gAlice = (AliRun*)file->Get("gAlice");
      gAlice->GetEvent(0);
      
      TClonesArray *particles=gAlice->Particles();
      Int_t n=particles->GetEntriesFast();
      Float_t torad=TMath::Pi()/180;
      Float_t phi_min = slice*20 - 10;
      Float_t phi_max = slice*20 + 10;
      
      for (Int_t j=0; j<n; j++) {
	TParticle *p=(TParticle*)particles->UncheckedAt(j);
	if (p->GetFirstMother()>=0) continue;  //secondary particle
	if(p->Eta() < eta[0] || p->Eta() > eta[1]) continue;
	Double_t ptg=p->Pt(),pxg=p->Px(),pyg=p->Py(),pzg=p->Pz();
	Double_t phi_part = TMath::ATan2(pyg,pxg);
	if (phi_part < 0) phi_part += 2*TMath::Pi();
	
	if(phi_part < phi_min*torad || phi_part > phi_max*torad) {continue;}
	if(ptg<0.100) continue;
		
	printf("Generated particle %d, with pt %f phi %f\n",p->GetPdgCode(),ptg,phi_part);
      }
      file->Close();
      delete file;
    }
  
}
