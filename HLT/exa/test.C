void test(char *clufile="")
{
//  AliL3Logger l;
//  l.UnSet(AliL3Logger::kDebug);
//  l.UnSet(AliL3Logger::kAll);
//  l.Set(AliL3Logger::kInformational);
//  l.UseStdout();
//  l.UseStream();

//  TFile * in = new TFile(clufile,"read");
//  TFile * out= new TFile("dummy.root","recreate");

//  a = new AliLevel3(in,out); 
  a = new AliLevel3();
  a->UseBinaryInput("/prog/alice/data/Rawdata/6_patch/hg_84210_s1-3/");
//  a->UseBinaryInput("/prog/alice/data/Rawdata/6_patch/20tracks/");
//  a->DoRoi(0.1,0.2);


  Int_t phi_segments,eta_segments,trackletlength,tracklength;
  Int_t rowscopetracklet,rowscopetrack;
  Double_t min_pt_fit,maxangle,goodDist,hitChi2Cut;
  Double_t goodHitChi2,trackChi2Cut;


  phi_segments = 50;
  eta_segments = 100;
  trackletlength = 2;
  tracklength = 3;
  rowscopetracklet = 2;
  rowscopetrack = 2;
  min_pt_fit = 0;
  maxangle = 1.31;
  goodDist = 5;
  hitChi2Cut = 5;
  goodHitChi2 = 10;
  trackChi2Cut = 50;

  a->SetTrackerParam(phi_segments,eta_segments,trackletlength,tracklength,
                     rowscopetracklet,rowscopetrack,
                     min_pt_fit,maxangle,goodDist,hitChi2Cut,
                     goodHitChi2,trackChi2Cut);

//  a->SetTrackerParam();

  a->WriteFiles();
  a->ProcessEvent(1,3);
  a->DoBench("benchmark_0");

  delete a;
}
