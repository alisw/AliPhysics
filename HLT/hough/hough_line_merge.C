void hough_line_merge(char *rootfile,int patch,Bool_t MC=false)
{
  gStyle->SetOptStat(0);
  
  
  Double_t pi = TMath::Pi();
  Int_t xbin = 60;
  Int_t ybin = 60;
  Float_t xrange[2] = {-pi/2,pi/2};
  Float_t yrange[2] = {-40,40};
  
  TH2F *hist;
  
  Int_t xr[2] = {0,250};
  Int_t yr[2] = {-125,125};
  TH1F *ntracks = new TH1F("ntracks","",100,0,200);
  TH2F *raw = new TH2F("raw","",250,xr[0],xr[1],250,yr[0],yr[1]);
  TH2F *road = new TH2F("road","",250,0,250,250,yr[0],yr[1]);
  TH2F *fake = new TH2F("fake","",250,0,250,250,-125,125);
  TH2F *peaks = new TH2F("peaks","",xbin,xrange[0],xrange[1],ybin,yrange[0],yrange[1]);
  peaks->SetMarkerStyle(3);
  peaks->SetMarkerColor(2);
  
  int slice = 2,n_phi_segments=20;
  float eta[2] = {0.3,0.4};
  //float eta[2] = {0.05,0.06};
 
  AliL3HoughTransformer *a = new AliL3HoughTransformer(slice,patch,eta,n_phi_segments);
  a->GetPixels(rootfile,raw);
  
  AliL3HoughMaxFinder *b = new AliL3HoughMaxFinder("DPsi");
  
  const Int_t NRows[5][2] = {{ 0, 45},{46,77},{78,109},{110,141},{142,173}};
  //const int ref_row[4]={5,15,25,35};
  //int rowrange[4][2]={{0,10},{10,20},{20,30},{30,45}};
  int ref_row[7] = {3,9,15,21,27,33,39};
  int rowrange[7][2] = {{0,6},{6,12},{12,18},{18,24},{24,30},{30,36},{36,42}};
  
  double phirange[7][2] = {{-10.,-5.},{-7.5,-2.5},{-5.,0.},{-2.5,2.5},{0.,5.},{2.5,7.5},{5.,10.}};
  int npatches = 7;
  int count=0;

  //AliL3TrackArray **track_pieces = new AliL3TrackArray*[npatches];
  AliL3TrackArray *track_pieces = new AliL3TrackArray("AliL3HoughTrack");

  double pran[2] = {-10,10};
  
  // int pa=3;
  for(int row_pat=0; row_pat < 7; row_pat++)
    {
      
      for(int pat=0; pat<7; pat++)
	{
	  
	  TH2F *histo = new TH2F("histo","Parameter space",xbin,xrange[0],xrange[1],ybin,yrange[0],yrange[1]);
	  
	  printf("transforming ref_row %d rowrange %d %d phirange %f %f\n",ref_row[row_pat],rowrange[row_pat][0],rowrange[row_pat][1],phirange[pat][0],phirange[pat][1]);
	  
	  a->Transform2Line(histo,ref_row[row_pat],rowrange[row_pat],phirange[pat]);
	  //a->Transform2Line(histo,ref_row[row_pat],rowrange[row_pat],pran);

	  AliL3TrackArray *tmp = (AliL3TrackArray*)b->FindMaxima(histo,rowrange[row_pat],ref_row[row_pat]);
	  
	  track_pieces->AddTracks(tmp);
	  
	  //AliL3HoughTrack *tra = (AliL3HoughTrack*)b->GetTracks()->GetCheckedTrack(0);
	  //printf("psi %f rowrange %d %d\n",tra->GetPsiLine(),tra->GetFirstRow(),tra->GetLastRow());
	  
	  delete histo;
	  delete tmp;
	}
      
    }
  
  
  printf("Transformation complete, genereted %d tracks\n",track_pieces->GetNTracks());

    
  //track_pieces->Sort();
  TH2F *cross = new TH2F("cross","",250,xr[0],xr[1],250,yr[0],yr[1]);  
  Double_t xy[2];
  //for(int i=0; i<npatches; i++)
  for(int i=0; i<4; i++)
    {
      if(track_pieces->GetNTracks()==0) continue;
      AliL3HoughTrack *track = (AliL3HoughTrack*)track_pieces->GetCheckedTrack(i);
      if(!track) {printf("NO TRACK\n"); continue;}
      //for(int j=rowrange[pa][0]; j<=rowrange[pa][1]; j++)
      /*
      for(int j=0; j<174; j++)
	{
	  track->GetLineCrossingPoint(j,xy);
	  cross->Fill(xy[0],xy[1],1);
	}
      */
      //printf("weight %d psi %f D %f rowrange %d %d\n",track->GetWeight(),track->GetPsiLine(),track->GetDLine(),track->GetFirstRow(),track->GetLastRow());
      
    }
  cross->SetMarkerStyle(6);
  cross->SetMarkerColor(2);
  //cross->Draw();
  
  

  AliL3HoughTransformer *a2 = new AliL3HoughTransformer(slice,pat,eta,n_phi_segments);
  
  xbin = 60 ;
  ybin = 60;
  xrange[0] = -0.006;
  xrange[1] = 0.006;
  yrange[0] = -0.26;
  yrange[1] = 0.26;
  hist = new TH2F("hist","Parameter space",xbin,xrange[0],xrange[1],ybin,yrange[0],yrange[1]);
  
  TH2F *cross_helix = new TH2F("cross_helix","",250,xr[0],xr[1],250,yr[0],yr[1]);  
  
  a2->TransformLines2Circle(hist,(AliL3TrackArray*)track_pieces);
      
  AliL3HoughMaxFinder *b2 = new AliL3HoughMaxFinder("KappaPhi");
  AliL3TrackArray *final_tracks = b2->FindMaxima(hist);
  
  printf("Number of track candidates before checking %d\n",final_tracks->GetNTracks());
  AliL3HoughEval *eval = new AliL3HoughEval(slice);
  eval->SetTransformer(a);
  eval->LookInsideRoad(final_tracks);
  
  printf("Total number of tracks after checking %d\n",final_tracks->GetNTracks());
  for(int i=0; i<3; i++)//final_tracks->GetNTracks(); i++)
    {
      
      AliL3HoughTrack *track = (AliL3HoughTrack*)final_tracks->GetCheckedTrack(i);
      if(!track) continue;
      
      for(int r=0; r<174; r++)
	{
	  Float_t xyz_helix[3];
	  track->GetCrossingPoint(r,xyz_helix);
	  cross_helix->Fill(xyz_helix[0],xyz_helix[1],1);
	}
      
      peaks->Fill(track->GetKappa(),track->GetPhi0(),1);
      printf("Found track, pt %f phi0 %f\n",track->GetPt(),track->GetPhi0());
    }
  cross_helix->SetMarkerStyle(6);
  cross_helix->SetMarkerColor(2);
  TCanvas *c1 = new TCanvas("c1","",900,900);
  SetCanvasOptions(c1);
  c1->Divide(2,2);
  c1->cd(1);
  hist->Draw("box");

  peaks->SetMarkerStyle(3);
  peaks->SetMarkerColor(2);
  //peaks->Draw("same");

  c1->cd(2);
  cross_helix->Draw();

  c1->cd(3);
  raw->Draw("");
  cross_helix->DrawCopy("same");
  //cross->Draw("same");

  c1->cd(4);
  //hist->Draw("box");
  raw->DrawCopy();
  
  delete track_pieces;
  delete a;
  delete b;
  delete a2;
  delete eval;
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
