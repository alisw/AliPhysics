int MakeTrend(char *infile,int run) {

  //char input_file[300];
  //sprintf(input_file,"Run_%d/AnalysisResults.root",run);
  
  TStopwatch timer;
  timer.Start();
  gSystem->SetIncludePath("-I$ROOTSYS/include -I$ALICE_PHYSICS/include -I$ALICE_PHYSICS/ITS -I$ALICE_PHYSICS -I$ALICE_PHYSICS/TRD");
  //gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/TPC/macros/ConfigOCDB.C");
  
  gSystem->Load("libSTAT");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libANALYSIScalib");
  gSystem->Load("libCORRFW");
  gSystem->Load("libTPCcalib");
  gSystem->Load("libTRDcalib");
  gSystem->Load("libT0calib");
  gSystem->Load("libTOFcalib");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libANALYSIScalib");
  gSystem->Load("libTender");
  gSystem->Load("libPWGPP");
                    
  
  /*
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libTender");
  gSystem->Load("libCORRFW");
  gSystem->Load("libPWG0base");
  gSystem->Load("libPWG0dep");
  gSystem->Load("libPWG0selectors");
  gSystem->Load("libPWG1");
    gSystem->Load("libPWG2");
  gSystem->Load("libPWG4base");
  gSystem->Load("libPWG3muon");
  gSystem->Load("libPWG3muondep");
  gSystem->Load("libPWG2forward");
   gSystem->Load("libPWG4PartCorrBase");
   gSystem->Load("libPWG4PartCorrDep");
 
  */

  // config OCDB
  //
  
  // ConfigOCDB(run,"local:///lustre/alice/alien/alice/data/2011/OCDB");
  
  char *outfile = "trending.root";
  
  if (!infile) return -1;
  if (!outfile) return -1;
  TFile *f =0;
  f=TFile::Open(infile,"read");
  if (!f) {
    printf("File %s not available\n", infile);
    return -1;
  }

  cout<< "test breakpoint -------------------------------------"<<endl; 
  TList* list = 0;
  list = dynamic_cast<TList*>(f->Get("TPC")); 
  cout<< "test0"<< endl;
  if (!list) { list = dynamic_cast<TList*>(f->Get("TPCQA")); }
  cout<< "test1"<< endl;
  if (!list) { list = dynamic_cast<TList*>(f->Get("TPCQA_v0_c0")); }
  if (!list) { list = dynamic_cast<TList*>(f->Get("TPCQA_v0_c30")); }
  if (!list) { list = dynamic_cast<TList*>(f->Get("TPCQA_v0_c70")); }
  if (!list) { list = dynamic_cast<TList*>(f->Get("TPC_PerformanceQA/TPCQA")); }
  if (!list) { list = dynamic_cast<TList*>(f->Get("TPC_PerformanceQA")); }
  if (!list) { list = dynamic_cast<TList*>(f->Get("ITSTPCMatch")); }
  cout<< "test2"<< endl;
  if (!list) {
    printf("QA %s not available\n", infile);
    return -1;
  } 
  cout<< "test3"<< endl;
  
  AliPerformanceTPC* pTPC = 0;
  AliPerformanceDEdx* pTPCgain = 0; 
  AliPerformanceMatch* pTPCmatch = 0; 
  AliPerformanceMatch* pTPCPull = 0; 
  AliPerformanceMatch* pConstrain = 0; 
  if (list) {  pTPC = dynamic_cast<AliPerformanceTPC*>(list->FindObject("AliPerformanceTPC")); }
  if (list) {  pTPCgain = dynamic_cast<AliPerformanceDEdx*>(list->FindObject("AliPerformanceDEdxTPCInner")); }
  if (list) {  pTPCmatch = dynamic_cast<AliPerformanceMatch*>(list->FindObject("AliPerformanceMatchTPCITS")); }
  if (list) {  pTPCPull = dynamic_cast<AliPerformanceMatch*>(list->FindObject("AliPerformanceMatchITSTPC")); }
  if (list) {  pConstrain = dynamic_cast<AliPerformanceMatch*>(list->FindObject("AliPerformanceMatchTPCConstrain")); }

  
  cout<< "ss "<<pTPC<<pTPCgain<<pTPCmatch <<pTPCPull<<pConstrain<<endl;
  /*
  cout<< ((TH2D*)(pTPC->GetTPCTrackHisto()->Projection(4,5)))->GetEntries()<<endl;
  
  
  TH1* his1D=0;
  //TH2* his2D=0;
  TH3* his3D_0=0;
  
  
  his3D_0 = dynamic_cast<TH3*>(pTPC->GetHistos()->FindObject("h_tpc_track_all_recvertex_0_5_7"));
  his3D_0->GetYaxis()->SetRangeUser(-1,1);
  his3D_0->GetZaxis()->SetRangeUser(0.25,10);
  
  his1D = his3D_0->Project3D("x");
  Double_t meanTPCncl= his1D->GetMean();
  
  cout<< "ss "<< meanTPCncl<<endl;
  */

  Int_t returncode = 0;
  AliTPCPerformanceSummary::WriteToFile(pTPC, pTPCgain, pTPCmatch ,pTPCPull,pConstrain, outfile, run);
  if (f) { delete f; f=0; }
  
  ofstream fout("trend_exist");
  fout<<9<<endl;
  fout.close();
  
  return 9;

}
