void test(int minslice,int maxslice,int nevent=1,char *clufile="")
{
  AliL3Logger l;
  //  l.UnSet(AliL3Logger::kDebug);
  //  l.UnSet(AliL3Logger::kAll);
  l.Set(AliL3Logger::kAll);
  //l.UseStdout();
  l.UseStream();
  
  TFile * in = new TFile(clufile,"read");
  TFile * out= new TFile("dummy.root","recreate");
  
  //a = new AliLevel3(in,out); 
  
  
  for(Int_t ev=0; ev<nevent; ev++)
    {
      a = new AliLevel3(in,out);
      // a->UseBinaryInput("/heim/franken/binary/onepatch/");
      //a->UseBinaryInput("/prog/alice/data/Rawdata/6_patch/20tracks_s2/");
      //a->UseBinaryInput("/prog/alice/data/Rawdata/6_patch/hg_84210_s1-3/");
      //a->UseBinaryInput("/prog/alice/data/Rawdata/1_patch/pp/event_0/");
      /*
	char fname[256];
	sprintf(fname,"/prog/alice/data/Rawdata/1_patch/pp/event_%d/",ev);
	a->UseBinaryInput(fname);
	sprintf(fname,"/prog/alice/data/Rawdata/1_patch/pp/recon_%d/",ev);
	a->WriteFiles(fname);
      */
      //  a->DoRoi();
      //a->DoMc();
      
      Int_t phi_segments,eta_segments,trackletlength,tracklength;
      Int_t rowscopetracklet,rowscopetrack;
      Double_t min_pt_fit,maxangle,goodDist,hitChi2Cut;
      Double_t goodHitChi2,trackChi2Cut,maxphi,maxeta;
      
      phi_segments = 50;//50;
      eta_segments = 100;//100;
      trackletlength = 5;
      tracklength = 10;
      rowscopetracklet = 2;
      rowscopetrack = 2;
      min_pt_fit = 0;
      maxangle = 1.31;
      goodDist = 5;
      maxphi=100;
      maxeta=100;
      hitChi2Cut = 200;
      goodHitChi2 = 50;
      trackChi2Cut = 100;
      /*
      //main vertex tracking parameters:
      a->SetTrackerParam(phi_segments,eta_segments,trackletlength,tracklength,
      rowscopetracklet,rowscopetrack,
      min_pt_fit,maxangle,goodDist,hitChi2Cut,
      goodHitChi2,trackChi2Cut,50,maxphi,maxeta,kTRUE);
      */
      a->SetTrackerParam();
      
      //non vertex tracking parameters:
      trackletlength = 2;
      tracklength = 10;
      rowscopetracklet = 2;
      rowscopetrack = 2;
      min_pt_fit = 0;
      maxangle = 1.31;
      goodDist = 0;
      hitChi2Cut = 100;
      goodHitChi2 = 5;
      trackChi2Cut = 80;
      
      a->SetTrackerParam(phi_segments,eta_segments,trackletlength,tracklength,
			 rowscopetracklet,rowscopetrack,
			 min_pt_fit,maxangle,goodDist,hitChi2Cut,
			 goodHitChi2,trackChi2Cut,50,maxphi,maxeta,kFALSE);
      
      // a->SetTrackerParam();
      //a->DoNonVertexTracking();
      a->WriteFiles();
      //a->ProcessEvent(0,2);
      a->ProcessEvent(minslice,maxslice);
      //a->DoBench("benchmark_0");
      
      delete a;
    }
}
