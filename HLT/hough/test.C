void test(char *file="/prog/alice/data/Rawdata/6_patch/hg_42105_s1-3/",bool bin=true)
{
  
  Int_t fNEtaSegments = 1 ;
  
  Char_t histname[50];
  Int_t i;
  
  int slice=1,patch=0;
  
  fHoughTransformer = new AliL3HoughTransformer(slice,patch,fNEtaSegments);
  fMemHandler = new AliL3FileHandler();
  fMaxFinder = new AliL3HoughMaxFinder("KappaPhi");
  
  fTransform = new AliL3Transform();
  
  fHoughTransformer->CreateHistograms();//64,0,3.1415,128,-120,120);
  fHoughTransformer->SetThreshold(10);
  fTracks = new AliL3TrackArray("AliL3HoughTrack");
  
  
  TH2F *road = new TH2F("road","",250,0,250,250,-125,125);
  TH2F *peaks = new TH2F("peaks","",64,-0.006,0.006,64,-0.26,0.26);
  peaks->SetMarkerStyle(3);
  peaks->SetMarkerColor(2);
  road->SetMarkerStyle(6);
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
      fMemHandler->Free();
      TFile *rootfile = new TFile(file);
      fMemHandler->SetAliInput(rootfile);
      fMemHandler->Init(slice,patch,NRows[patch]);
      fMemHandler->Init(fTransform);
      digits=(AliL3DigitRowData *)fMemHandler->AliDigits2Memory(ndigits); 
      rootfile->Close();
      delete rootfile;
    }
  fHoughTransformer->SetInputData(ndigits,digits);
  fEval = new AliL3HoughEval(fHoughTransformer);
  fEval->SetNumOfRowsToMiss(2);
  //fEval->SetNumOfPadsToLook(2);
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
      /*
      for(Int_t e=0; e<fNEtaSegments; e++)
	{
	  histo = (AliL3Histogram*)fHoughTransformer->GetHistogram(e);
	  if(!histo) continue;
	  
	  fMaxFinder->SetHistogram(histo);
	  Int_t n=10;
	  Float_t x[10];
	  Float_t y[10];
	  fMaxFinder->FindPeak1(x,y,n);
	  track = new AliL3HoughTrack();
	  track->SetTrackParameters(x[0],y[0],1);
	  
	  //if(!fEval->LookInsideRoad(track,e))
	  //continue;
	  for(int t=0; t<176; t++)
	    {
	      float xyz_tr[3];
	      track->GetCrossingPoint(t,xyz_tr);
	      road->Fill(xyz_tr[0],xyz_tr[1],1);
	    }
	  
	}
      */
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
  fHoughTransformer->GetHistogram(eind)->Draw("box");
  //  peaks->Draw("same");
  c1->cd(2);
  image->Draw("");
  //road->Draw("same");
} 

void process(char *path="/prog/alice/data/Rawdata/6_patch/hg_42105_s1-3/",bool bin=true)
{
  double torad = 3.1415/180;
  a = new AliL3Hough(path,bin,1);
  a->TransformSlice(1);
  
  hist = (AliL3Histogram*)a->AddHistograms();
  
  //hist->SetThreshold(10000);
  
  b = new AliL3HoughMaxFinder("KappaPhi");
  b->SetHistogram(hist);
  
  Int_t xbin,ybin;
  
  Int_t n=10;
  Float_t x[10];
  Float_t y[10];
  b->FindPeak1(x,y,n);
  printf("peak at pt %f phi0 %f\n",0.2*0.003/x[0],y[0]/torad);
  
  track = new AliL3HoughTrack();
  track->SetTrackParameters(x[0],y[0],1);
  
  image = new AliL3Histogram("image","",250,0,250,250,-125,125);
  a->Evaluate(image);
  TH2F *road = new TH2F("road","",250,0,250,250,-125,125);
  road->SetMarkerStyle(5);
  road->SetMarkerColor(2);
  
  float xyz[3];
  for(int i=0; i<176; i++)
    {
      if(i%10) continue;
      track->GetCrossingPoint(i,xyz);
      road->Fill(xyz[0],xyz[1],1);
    }
  c1 = new TCanvas("c1","",1000,500);
  c1->Divide(2);
  c1->cd(1);
  hist->Draw("box");
  c1->cd(2);
  image->Draw();
  road->Draw("same");
}
