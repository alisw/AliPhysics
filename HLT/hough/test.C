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

  fHoughTransformer->CreateHistograms();
  fHoughTransformer->SetThreshold(3);
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
      fHoughTransformer->Transform();
      double final = AliL3Benchmark::GetCpuTime();
      printf("done in %f ms\n",(final-init)*1000);
      
      good_count=0;
      
      for(Int_t e=0; e<fNEtaSegments; e++)
	{
	  histo = (AliL3Histogram*)fHoughTransformer->GetHistogram(e);
	  if(!histo) continue;
	  
	  fMaxFinder->SetHistogram(histo);
	  track = (AliL3HoughTrack*)fMaxFinder->FindPeak1();
	  peaks->Fill(track->GetKappa(),track->GetPhi0(),1);
	  
	  if(!fEval->LookInsideRoad(track,e))
	    continue;
	  for(int t=0; t<175; t++)
	    {
	      float xyz_tr[3];
	      track->GetCrossingPoint(t,xyz_tr);
	      road->Fill(xyz_tr[0],xyz_tr[1],1);
	    }
	  
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
  fHoughTransformer->GetHistogram(eind)->Draw("box");
  peaks->Draw("same");
  c1->cd(2);
  image->Draw("colz");
  road->Draw("same");
}
