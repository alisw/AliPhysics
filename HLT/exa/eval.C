// $Id$

void eval(char *inFile)
{
  // Connect the Root Galice file containing Geometry, Kine and Hits

  file = new TFile(inFile);
  
  gAlice = (AliRun*)file->Get("gAlice");
  if (gAlice) printf("AliRun object found on file\n");
  if (!gAlice) gAlice = new AliRun("gAlice","Alice test program");
  
  gAlice->GetEvent(0);
  AliTPC *TPC = (AliTPC*)gAlice->GetDetector("TPC");      
  
  AliTPCParam *param = (AliTPCParam*)file->Get("75x40_100x60");
  
  TFile *nfile = new TFile("/prog/alice/data/exact_clusters.root","RECREATE");
  nfile->cd();

  //setup AliTPCClustersArray
  AliTPCClustersArray * arr=new AliTPCClustersArray;
  arr->SetClusterType("AliComplexCluster");
  arr->Setup(param);
  TPC->SetParam(param);
  arr->MakeTree();

  TPC->SetClustersArray(arr); 
  TPC->Hits2ExactClustersSector(1);
  TPC->Hits2ExactClustersSector(37);
    

  //write results
  char treeName[100];
  sprintf(treeName,"TreeCExact_%s",param->GetTitle());
  TPC->GetClustersArray()->GetTree()->Write(treeName);
  param->Write(param->GetTitle());
  file->Close();
  nfile->Close();
  
  return;
}

void calc(char *rootfile,char *digitsfile,char *cfile)
{

  int slice[2] = {1,1};
  a = new AliL3Evaluate(rootfile,digitsfile,slice);
  a->SetupSlow("tracks.raw","./");
  ntuppel = (TNtuple*)a->EvaluatePoints(cfile);

  f = new TFile("results.root","RECREATE");
  f->cd();
  ntuppel->Write();
  f->Close();

}

void plot(char *file)
{
  gStyle->SetOptFit(0110);
  f = new TFile(file);
  
  TCanvas *c1 = new TCanvas("c1","",2);
  TH1F *hist = new TH1F("hist","",100,-2,2);
  
  ntuppel->Draw("resy>>hist","pt>1");

  float rms = hist->GetRMS();
  hist->SetXTitle("#delta_{T} [cm]");

  TF1 *f1 = new TF1("f1","gaus",-rms,rms);
  hist->Fit("f1","R");

}
