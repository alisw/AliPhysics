// $Id$

/**
   Shows how to use the hough code. Stores tracks parameters
   in files.
*/

void runhough(Int_t slice,Char_t *path,Int_t n_eta_segments,Int_t vseg=-1)
{
  
  //AliL3Transform::Init("/prog/alice/data/new/fixed-slice0/");
  AliL3Transform::Init(path);  

  hough = new AliL3Hough();
  Bool_t binary = kTRUE;
  Bool_t bit8 = kTRUE;
  hough->Init(path,binary,n_eta_segments,bit8);
  
  hough->GetMaxFinder()->SetThreshold(14000);
  
  hough->ReadData(slice);
  hough->Transform();
  hough->AddAllHistograms();
  hough->FindTrackCandidates();
  
  //hough->Evaluate();
  tracks = (AliL3TrackArray*)hough->GetTracks(0);
  for(int i=0; i<tracks->GetNTracks(); i++)
    {
      track = (AliL3HoughTrack*)tracks->GetCheckedTrack(i);
      if(!track) continue;
      cout<<"pt "<<track->GetPt()<<" psi "<<track->GetPsi()<<" eta "<<track->GetEta()<<" etaindex "<<track->GetEtaIndex()<<" weight "<<track->GetWeight()<<endl;
      if(vseg==-1) vseg=track->GetEtaIndex();
    }
  
  hough->WriteTracks(slice);
  cout<<"Found in slice " << slice << " total "<<tracks->GetNTracks()<<" tracks"<<endl;
  display(hough,vseg);
  
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
  hough->GetTransformer(0)->GetHistogram(eta_index)->Draw("box");
  
  //Draw the tracks
  c2 = new TCanvas("c2","",2);
  digitd->Draw();
  trackd->Draw("same");
  ((TH1F*)trackd->GetRootHisto())->SetMarkerColor(2);
}

struct GoodTrack
{
  Int_t event;
  Int_t label;
  Double_t eta;
  Int_t code;
  Double_t px,py,pz;
  Double_t pt;
  Int_t nhits;
};

void geteff(char *fname)
{
  GoodTrack gt[15000];
  int counter=0;
  ifstream in(fname);
  if(!in)
    {
      cerr<<"Could not open "<<fname<<endl;
      return;
    }
  while(in>>gt[counter].event>>gt[counter].label>>gt[counter].code
	>>gt[counter].px>>gt[counter].py>>gt[counter].pz>>gt[counter].pt>>gt[counter].eta>>gt[counter].nhits)
    counter++;
  
  char filename[100];
  file = new AliL3MemHandler();
  
  
}
