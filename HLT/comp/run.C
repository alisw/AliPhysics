void run()
{
  
  a = new AliL3Modeller();
  //a->Init(1,0,"/prog/alice/data/Rawdata/6_patch/1track_s1/");
  a->Init(1,0,"/prog/alice/data/Rawdata/1_patch/hg_1000_s1-3/");
  a->Process();
  
  a->WriteRemaining("test.raw");
  return;
  
  tracks = a->GetTracks();
  
  //plot(tracks);
  
  for(int i=0; i<tracks->GetNTracks(); i++)
    {
      track = (AliL3ModelTrack*)tracks->GetCheckedTrack(i);
      if(!track) continue;
      track->FillModel();
    }
  
  c = new AliL3Compress();
  c->Write2File(tracks,"mtracks.raw");
    
  //c->ReadFile("mtracks.raw");
  c->CompressFile("mtracks.raw","ctracks.raw");
  //c->ExpandFile();
  delete c;
  
  //delete a;
}

void plot(AliL3TrackArray *tracks)
{
  hist = new TH1F("hist","",256,0,255);

  for(int i=0; i<tracks->GetNTracks(); i++)
    {
      track = (AliL3ModelTrack*)tracks->GetCheckedTrack(i);
      if(!track) continue;
      if(track->GetNHits()<150) break;
      
      track->Print();
      
      for(int j=0; j<30; j++)
	{
	  Float_t res;
	  if(track->GetPadResidual(j,res))
	    hist->Fill(res);
	}
      
      
    }
  return;
  c1 = new TCanvas("c1","",2);
  hist->Draw();
}

void getcharge()
{
  
  TNtuple *ntuppel = new TNtuple("ntuppel","","charge");

  int patch=0;
  file = new AliL3MemHandler();
  for(int event=0; event<25; event++)
    {
      for(int slice=0; slice<35; slice++)
	{
	  char fname[100];
	  sprintf(fname,"/prog/alice/data/Rawdata/1_patch/pp/recon_%d/points_%d_%d.raw",event,slice,patch);
	  file->SetBinaryInput(fname);
	  
	  UInt_t npoints;
	  AliL3SpacePointData *points = (AliL3SpacePointData *) file->Allocate();
	  file->Binary2Memory(npoints,points);
	  file->CloseBinaryInput();
	  
	  for(int i=0; i<npoints; i++)
	    {
	      //cout<<""<<points[i].fX<<" "<<points[i].fY<<" "<<points[i].fZ<<endl;
	      //cout<<"Charge "<<points[i].fCharge<<endl;
	      Float_t charge[1] = {(float)points[i].fCharge};
	      ntuppel->Fill(charge);
	    }
	  file->Free();
	}
    }
  delete file;
  rootfile = TFile::Open("average_charge.root","RECREATE");
  ntuppel->Write();
  rootfile->Close();

}

void plotcharge()
{
  gStyle->SetOptFit(1);
  file = TFile::Open("average_charge.root");
  
  hist = new TH1F("hist","",100,0,2000);
  ntuppel->Draw("charge>>hist","","goff");
  
  f1 = new TF1("f1","landau",0,2000);
  hist->Draw();
  hist->Fit(f1,"R");
}
