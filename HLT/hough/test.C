void test(char *file="/prog/alice/data/Rawdata/1_patch/pp/event_0/",bool bin=true)
{
  AliL3Logger l;
//  l.UnSet(AliL3Logger::kDebug);
//  l.UnSet(AliL3Logger::kAll);
  l.Set(AliL3Logger::kAll);
  //l.UseStdout();
  l.UseStream();
  
  Int_t fNEtaSegments = 1;
  
  Char_t histname[50];
  Int_t i;
  
  int slice=1,patch=0;
  
  fHoughTransformer = new AliL3HoughTransformer(slice,patch,fNEtaSegments);
  fMemHandler = new AliL3FileHandler();
  fMaxFinder = new AliL3HoughMaxFinder("KappaPhi");
  
  fTransform = new AliL3Transform();
  
  fHoughTransformer->CreateHistograms(64,-0.012,0.012,64,-0.35,0.35);
  fHoughTransformer->SetThreshold(10);
  fTracks = new AliL3TrackArray("AliL3HoughTrack");
    
  TH2F *road = new TH2F("road","",250,0,250,250,-125,125);
  TH2F *peaks = new TH2F("peaks","",64,-0.006,0.006,64,-0.26,0.26);
  peaks->SetMarkerStyle(3);
  peaks->SetMarkerColor(2);
  road->SetMarkerStyle(5);
  road->SetMarkerColor(2);
  
  Char_t name[256];
  
  UInt_t ndigits=0;
  AliL3DigitRowData *digits =0;
  
  if(bin) //read data from binary files
    {
      Char_t fPath[256];
      strcpy(fPath,file);
      fMemHandler->Free();
      sprintf(name,"%sdigits_%d_%d.raw",fPath,slice,patch);
      fMemHandler->SetBinaryInput(name);
      digits = (AliL3DigitRowData *)fMemHandler->CompBinary2Memory(ndigits);
      fMemHandler->CloseBinaryInput();
    }
  else //read data from root file
    {
      const Int_t NRows[6][2] = {{0,31},{32,63},{64,91},{92,119},{120,143},{144,175}};
      const int nr[2] = {0,175};
      fMemHandler->Free();
      fMemHandler->SetAliInput(file);
      fMemHandler->Init(slice,patch,NRows[patch]);
      fMemHandler->Init(fTransform);
      digits=(AliL3DigitRowData *)fMemHandler->AliDigits2Memory(ndigits); 
      fMemHandler->CloseAliInput();
    }
    
  fHoughTransformer->SetInputData(ndigits,digits);
  fEval = new AliL3HoughEval();
  fEval->InitTransformer(fHoughTransformer);
  //fEval->SetNumOfRowsToMiss(0);
  //fEval->SetNumOfPadsToLook(2);
  //fEval->RemoveFoundTracks();
  AliL3HoughTrack *track;
  AliL3Histogram *histo;
  Int_t good_count;
  int eind=0;
  while(1)
    {
      printf("Transforming\n");
      double init = AliL3Benchmark::GetCpuTime();
      fHoughTransformer->TransformCircle();
      double final = AliL3Benchmark::GetCpuTime();
      printf("done in %f ms\n",(final-init)*1000);
      good_count=0;
      for(Int_t e=0; e<fNEtaSegments; e++)
	{
	  histo = (AliL3Histogram*)fHoughTransformer->GetHistogram(e);
	  if(!histo) continue;
	  
	  fMaxFinder->SetHistogram(histo);
	  Int_t n=10;
	  Int_t weight[10];
	  Float_t x[10];
	  Float_t y[10];
	  fMaxFinder->FindPeak1(x,y,weight,n);
	  track = new AliL3HoughTrack();
	  track->SetTrackParameters(x[0],y[0],1);
	  
	  //	  if(!fEval->LookInsideRoad(track,e))
	  // continue;
	  
	  
	  if(e==eind)
	    {
	      for(int t=0; t<176; t++)
		{
		  if(t%10) continue;
		  float xyz_tr[3];
		  track->GetCrossingPoint(t,xyz_tr);
		  road->Fill(xyz_tr[0],xyz_tr[1],1);
		}
	    }
	  delete track;
	  
	}
	  
      break;
      if(good_count==0)
	break;
    }
  

  image = new AliL3Histogram("image","",250,0,250,250,-125,125);
  fEval->DisplayEtaSlice(eind,image);
  
  c1 = new TCanvas("c1","",1000,500);
  c1->Divide(2);
  c1->cd(1);
  histo = (AliL3Histogram*)fHoughTransformer->GetHistogram(eind);
  if(!histo)
    {printf("No histogram\n"); return;}
  histo->Draw("box");

  peaks->Draw("same");
  c1->cd(2);
  image->Draw("");
  road->Draw("same");
} 

void process(char *digitfile="/prog/alice/data/Rawdata/1_patch/pp/event_0/",bool bin=true,char *trackfile="good_tracks")
{
  AliL3Logger l;
  //  l.UnSet(AliL3Logger::kDebug);
  //  l.UnSet(AliL3Logger::kAll);
  //l.Set(AliL3Logger::kError);
  //l.UseStdout();
  l.UseStream();
  
  double torad = 3.1415/180;
  int slice=1;
  a = new AliL3Hough(digitfile,bin,100);
 
  a->ReadData(slice);
  a->Transform();
  a->AddAllHistograms();
  a->FindTrackCandidates();
  // a->Evaluate(2);
  a->WriteTracks();
  return;
  //a->MergePatches();

  //a->MergeInternally();
  
  //tracks = (AliL3TrackArray*)(a->GetInterMerger())->GetOutTracks();
  //tracks = (AliL3TrackArray*)(a->GetMerger())->GetOutTracks();
  tracks = (AliL3TrackArray*)a->GetTracks(0);
  //a->GetEval(0)->CompareMC(tracks,"good_tracks_hg4000_s1");
 
  for(int i=0; i<tracks->GetNTracks(); i++)
    {
      AliL3HoughTrack *track = (AliL3HoughTrack*)tracks->GetCheckedTrack(i);
      if(!track) continue;
      printf("found pt %f phi0 %f eta %f weight %d rowange %d %d\n",track->GetPt(),track->GetPhi0(),track->GetEta(),track->GetWeight(),track->GetFirstRow(),track->GetLastRow());
    }
  return;

}

void cf()
{
    AliL3Logger l;
  //  l.UnSet(AliL3Logger::kDebug);
  //  l.UnSet(AliL3Logger::kAll);
  //l.Set(AliL3Logger::kError);
  //l.UseStdout();
  l.UseStream();
  
  AliL3FileHandler *file = new AliL3FileHandler();
  file->SetBinaryInput("tracks.raw");
  fTracks = new AliL3TrackArray("AliL3HoughTrack");
  file->Binary2TrackArray(fTracks);
  cout<<"Ntracks "<<fTracks->GetNTracks()<<endl;
  file->CloseBinaryInput();
  delete file;
  

  return;
  /*
  a = new AliL3ClusterFinder();
  a->LoadData();
  a->Process();
  */
}


void HistogramParticles(char *trackfile)
{
  struct GoodTrack 
  {
    
    Int_t label;
    Double_t eta;
    Int_t code;
    Double_t px,py,pz;
    Double_t pt;
    Int_t nhits;
  };
  struct GoodTrack goodtracks[15000];
  Int_t nt=0;
  ifstream in(trackfile);
  if(in)
    {
      printf("Reading good tracks from file %s\n",trackfile);
      while (in>>goodtracks[nt].label>>goodtracks[nt].code>>
	     goodtracks[nt].px>>goodtracks[nt].py>>goodtracks[nt].pz>>
	     goodtracks[nt].pt>>goodtracks[nt].eta>>goodtracks[nt].nhits) 
	{
	  nt++;
	  if (nt==15000) 
	    {
	      cerr<<"Too many good tracks"<<endl;
	      break;
	    }
	}
    }
  
  int fNEtaSegments = 150;
  Double_t etaslice = (0.9 - 0)/fNEtaSegments;
  Int_t *particles = new Int_t[fNEtaSegments];
  for(Int_t i=0; i<fNEtaSegments; i++)
    particles[i]=0;
  
  for(Int_t i=0; i<nt; i++)
    {
      //if(goodtracks[i].nhits < 150) continue;
      if(goodtracks[i].pt < 0.5) continue;
      Int_t particleindex = (Int_t)(goodtracks[i].eta/etaslice);
      if(particleindex < 0 || particleindex >= fNEtaSegments) continue;
      particles[particleindex]++;
    }
  
  hist = new TH1F("hist","",21,0,20);
  for(int i=0; i<fNEtaSegments; i++)
    hist->Fill(particles[i]);
  
  
  c1 = new TCanvas("c1","",2);
  hist->Draw();
  
}
