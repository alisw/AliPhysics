void hough(char *rootfile,int patch,Bool_t MC=false)
{
  gStyle->SetOptStat(0);
  
  Int_t xbin = 60;
  Int_t ybin = 60;
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
  TH1F *ntracks = new TH1F("ntracks","",100,0,200);
  TH2F *raw = new TH2F("raw","",250,xr[0],xr[1],250,yr[0],yr[1]);
  TH2F *road = new TH2F("road","",250,0,250,250,yr[0],yr[1]);
  TH2F *fake = new TH2F("fake","",250,0,250,250,-125,125);
  TH2F *peaks = new TH2F("peaks","",xbin,xrange[0],xrange[1],ybin,yrange[0],yrange[1]);
  peaks->SetMarkerStyle(3);
  peaks->SetMarkerColor(2);
  
  int slice = 2,n_phi_segments=30;
  //float eta[2] = {0.3,0.4};
  float eta[2] = {0.03,0.04};
  
  AliL3HoughTransformer *a = new AliL3HoughTransformer(slice,patch,eta,n_phi_segments);
  a->GetPixels(rootfile,raw);
  a->InitTemplates(hist);
  
  //a->Transform2Circle(hist,78);
  TStopwatch sw;
  a->Transform2Circle(hist);
  sw.Stop();
  printf("Transformation done in %f ms\n",sw.CpuTime()*1000);
  
  AliL3HoughMaxFinder *b = new AliL3HoughMaxFinder("KappaPhi");
  //b->SetThreshold(10000);
  AliL3TrackArray *tracks = b->FindMaxima(hist);
    
  AliL3HoughEval *c = new AliL3HoughEval(a);
  printf("number of peaks before %d\n",tracks->GetNTracks());
  
  /*
    for(int i=0; i<tracks->GetNTracks(); i++)
    {
      AliL3HoughTrack *track = (AliL3HoughTrack*)tracks->GetCheckedTrack(i);
      if(!track) printf("no track\n");
      printf("Pt %f Phi %f kappa %f weigth %d charge %d\n",track->GetPt(),track->GetPhi0(),track->GetKappa(),track->GetNHits(),track->GetCharge());
    }
  */
  //c->LookInsideRoad(a,tracks,road,fake,ntracks);
  c->LookInsideRoad(tracks);
  printf("Number of trackcandidates %d\n",tracks->GetNTracks());

  AliL3Transform *transform = new AliL3Transform();
  Int_t NRows[5][2] = {{ 0, 45},{46,77},{78,109},{110,141},{142,173}};
  for(Int_t i=0; i<tracks->GetNTracks(); i++)
    {
      AliL3HoughTrack *track = (AliL3HoughTrack*)tracks->GetCheckedTrack(i);
      if(!track) {printf("NO TRACK\n"); continue;}
      peaks->Fill(track->GetKappa(),track->GetPhi0(),1);
      //for(int row=0; row<174; row++)
      for(int row=NRows[patch][0]; row<=NRows[patch][1]; row++)
	{
	  Float_t xyz[3];
	  track->GetCrossingPoint(row,xyz);
	  //track->GetCrossingPoint(slice,row,xyz);
	    //printf("Track does not cross line at row %d\n",row);
	  //transform->Local2Global(xyz,slice);
	  road->Fill(xyz[0],xyz[1],1);
	}
      //float an[2] = {track->GetPhi0(),0};
      //transform->Local2GlobalAngle(an,slice);
      printf("pt %f phi0 %f weight %d\n",track->GetPt(),track->GetPhi0(),track->GetWeight());
    }
  
  
  TCanvas *c1 = new TCanvas("c1","",900,900);
  SetCanvasOptions(c1);
  c1->Divide(2,2);
  c1->cd(1);
  hist->Draw("lego1");
  //  peaks->Draw("same");

  //  TCanvas *c2 = new TCanvas("c2","",2);
  c1->cd(3);
  raw->Draw("");

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
  delete tracks;
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
