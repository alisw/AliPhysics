/**
   Run this macro for cluster finder and track follower (see steering class
   AliLevel3. 
*/

void runtracker(int minslice,int maxslice,char* path=0,char *rootfile=0,int nevent=1)
{
  AliL3Logger l;
  l.Set(AliL3Logger::kAll);
  l.UseStdout();
  //l.UseStream();

  char path_to_use[1024];
  if(!path) 
    strcpy(path_to_use,"/usr/local/anders/data/hg_42105_s1-3/");
  //strcpy(path_to_use,"/prog/alice/data/Rawdata/6_patch/hg_1000_s1-3/");
  //strcpy(path_to_use,"/prog/alice/data/Rawdata/6_patch/1track_s1/");
  //strcpy(path_to_use,"/prog/alice/data/Rawdata/6_patch/pp/event_0/");
  //strcpy(path_to_use,"/prog/alice/data/Rawdata/1_patch/pp/pileups/event_1/");
  else strcpy(path_to_use,path);

  //AliL3Transform::Init(path_to_use);  

  Int_t phi_segments,eta_segments,trackletlength,tracklength;
  Int_t rowscopetracklet,rowscopetrack;
  Double_t min_pt_fit,maxangle,goodDist,hitChi2Cut;
  Double_t goodHitChi2,trackChi2Cut,maxphi,maxeta;

  for(Int_t ev=0; ev<nevent; ev++)
    {
      if(!rootfile)
	a = new AliLevel3();
      else
	a = new AliLevel3(rootfile);

      a->Init(path_to_use);

      phi_segments = 50;
      eta_segments = 100;
      trackletlength = 3;
      tracklength = 5;
      rowscopetracklet = 2;
      rowscopetrack = 2;
      min_pt_fit = 0;
      maxangle = 1.31;
      goodDist = 5;
      maxphi=100;
      maxeta=100;
      hitChi2Cut = 100;
      goodHitChi2 = 20;
      trackChi2Cut = 50;
       
      //main vertex tracking parameters:
      a->SetTrackerParam(phi_segments,eta_segments,trackletlength,tracklength,
			 rowscopetracklet,rowscopetrack,
			 min_pt_fit,maxangle,goodDist,hitChi2Cut,
			 goodHitChi2,trackChi2Cut,50,maxphi,maxeta,kTRUE);

      //a->DoRoi();    /*do region of interest*/
      //a->DoMc();     /*do monte carlo identification*/
      a->WriteFiles(); /*enable output*/

      a->ProcessEvent(minslice,maxslice);
      //a->DoBench("benchmark_0");

      delete a;
    }
}


