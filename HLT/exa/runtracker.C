// $Id$

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
    //strcpy(path_to_use,"/usr/local/anders/data/hg_42105_s1-3/");
    strcpy(path_to_use,"/tmp/data/new/hijing/bfact1/1000/rawdata/");
  else strcpy(path_to_use,path);
  
  AliL3Transform::Init(path_to_use);  

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

      phi_segments = 50;   //devide the space into phi_segments and eta_segments
      eta_segments = 100;  //to access the search of points to that area!
      trackletlength = 3;  //number of hits a tracklet has to have
      tracklength = 5;     //number of hits a track has to have 
      rowscopetracklet = 2;//search range of rows for a tracklet
      rowscopetrack = 2;   //search range of rows for a track
      min_pt_fit = 0;      
      maxangle = 1.31;     //maximum angle for the three point look ahead
      goodDist = 5;        //threshold distance between two hits when building tracklets
      maxphi=100;          //maximum phi difference for neighboring hits
      maxeta=100;          //maximum eta difference for neighboring hits
      hitChi2Cut = 100;    //maximum chi2 of added hit to track
      goodHitChi2 = 20;    //stop looking for next hit to add if chi2 is less than goodHitChi2
      trackChi2Cut = 50;   //maximum chi2 for track after final fit
       
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


