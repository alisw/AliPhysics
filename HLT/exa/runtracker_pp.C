// $Id$

/**
   Run this macro for cluster finder and track 
   follower on pp pileup events.  
*/

void runtracker_pp(Int_t trigev,Char_t* path=0,Char_t *rootfile=0,Char_t *wpath="./")
{
  AliHLTLogger l;
  l.Set(AliHLTLogger::kAll);
  l.UseStdout();
  //l.UseStream();

  Char_t path_to_use[1024];
  if(!path){
    strcpy(path_to_use,"/prog/alice/data/V3-08-Rev-01/hijing/bfact1/rawdata/100/");
    //strcpy(path_to_use,"/usr/local/anders/data/hg_42105_s1-3/");
    //strcpy(path_to_use,"/tmp/data/new/hijing/bfact1/1000/rawdata/");
  } else strcpy(path_to_use,path);
  
  if(getenv("TRANSFORMCONFIGPATH")){
    AliHLTTransform::Init(getenv("TRANSFORMCONFIGPATH"));
  } else AliHLTTransform::Init(path_to_use);  

  Int_t phi_segments,eta_segments,trackletlength,tracklength;
  Int_t rowscopetracklet,rowscopetrack;
  Double_t min_pt_fit,maxangle,goodDist,hitChi2Cut;
  Double_t goodHitChi2,trackChi2Cut,maxphi,maxeta;

  Bool_t binary;

  if(!rootfile){
    a = new AliLevel3();
    binary = kTRUE;
  }
  else {
    a = new AliLevel3(rootfile);
    binary = kFALSE;
  }


  a->Init(path_to_use,binary,1); //set one patch (not 6) for the pileup!

  a->DoPileup(); // this is the difference to non pileup events!!!
                 // because binary files with MC info are not compressed!

  gROOT->LoadMacro("SetTrackingParameters_pp.C");
  //gROOT->LoadMacro("SetTrackingParameters_1000bf04.C");
  SetTrackingParameters(a);

#if 0
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
#endif


  //a->DoRoi();    /*do region of interest*/
  //a->DoMc();     /*do monte carlo identification*/

  Char_t Carry[1024];
  sprintf(Carry,"mkdir -p %s/",wpath);
  gSystem->Exec(Carry);
  a->WriteFiles(wpath); /*enable output*/

  a->ProcessEvent(0,35,trigev);
  //a->DoBench("benchmark_0");

  delete a;
}
