// $Id$

void Eval(char *rootfile="")
{
  AliL3Logger l;
//  l.UnSet(AliL3Logger::kDebug);
//  l.UnSet(AliL3Logger::kAll);
  l.Set(AliL3Logger::kError);
  l.UseStdout();
  //l.UseStream();
    
  int slice[2] = {0,35};
  e = new AliL3Evaluate(rootfile,slice);
  e->SetupSlow("tracks.raw",".");
  //e->SetupFast("tracks.raw","/nfs/david/subatom/alice/data/V3.04/fast/clusters/hg_8k_v0_s1-3_e0_cl.root",".");
  
  TNtuple *ntuppel = (TNtuple*)e->EvaluatePoints();
  file = new TFile("CFeval_nodeconv.root","RECREATE");
  file->cd();
  ntuppel->Write();
  file->Close();
  delete file;
}

void plotPt(char *rootfile)
{
  gStyle->SetStatColor(10);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1100);  
  
  c = new TCanvas("c","",2);
  SetCanvasOptions(c);
  
  f1 = new TFile(rootfile);

  hist = new TH1F("hist","",50,-10,10);
  SetTH1Options(hist);
  fNtuppel->Draw("(pt_found-pt_gen)/pt_gen*100>>hist","nHits>30");
  hist->GetXaxis()->SetTitle("%");
  hist->GetYaxis()->SetTitle("Counts");

  float rms = hist->GetRMS();
  printf("Rms value : %f\n",rms);
  
  TF1 *f = new TF1("f","gaus",-rms,rms);
  hist->Fit("f","R");
}

void plot(char *rootfile)
{
  gStyle->SetStatColor(10);
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1);
  
  file = new TFile(rootfile);
  hist = new TH1F("hist","",100,-0.6,0.6);
  SetTH1Options(hist);
  
  can = new TCanvas("can","Residuals",900,600);
  can->Divide(2);
  SetCanvasOptions(can);
  can->cd(1);
  //ntuppel->Draw("residual_trans>>hist","zHit < 50 && padrow > 55");//beta < 10*4.1515/180");
  ntuppel->Draw("resy>>hist","ptgen > 1.0");
  
  float rms = hist->GetRMS();
  printf("Rms value : %f\n",rms);
  
  TF1 *f = new TF1("f","gaus",-rms,rms);
  hist->Fit("f","R");
    
  hist->GetXaxis()->SetTitle("#delta_{T} [cm]");
  hist->GetYaxis()->SetTitle("Counts");
  
  f2 = new TFile("results_fast_oldparams.root");
  hist2 = new TH1F("hist2","",100,-0.6,0.6);
  SetTH1Options(hist2);
  can->cd(2);
  ntuppel_fast->Draw("residual_trans>>hist2","nHits>100 && pt>1.0 && padrow > 0 && zHit < 50");//dipangle < 20*3.1415/180");
  
  can->Update();
}


