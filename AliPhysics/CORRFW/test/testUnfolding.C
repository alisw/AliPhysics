
void testUnfolding() {
  TBenchmark b;

  b.Start("init");

  Load();

  AliLog::SetGlobalDebugLevel(0);

  // get the essential
  AliCFDataGrid *measured    = (AliCFDataGrid*) GetMeasuredSpectrum();
  AliCFDataGrid *generated   = (AliCFDataGrid*) GetGeneratedSpectrum();
  AliCFEffGrid  *efficiency  = (AliCFEffGrid*)  GetEfficiency();
  THnSparse     *response    = (THnSparse*)     GetResponseMatrix();
  
  // create a guessed "a priori" distribution using binning of MC
  THnSparse* guessed = CreateGuessed(((AliCFGridSparse*)generated->GetData())->GetGrid()) ;
  //----

  TF2* fit2D = new TF2("fit2D","[0]*([2]+[3]*y+[4]*y*y)*x*TMath::Exp(-x/[1])",0.1,0.8,-1.5,1.5);
  fit2D->SetParameter(0,1);
  fit2D->SetParameter(1,1);
  fit2D->SetParameter(2,1);
  fit2D->SetParameter(3,1);
  fit2D->SetParameter(4,1);
  fit2D->SetParLimits(1,0,1);

  AliCFUnfolding unfolding("unfolding","",2,response,efficiency->GetGrid(),((AliCFGridSparse*)measured->GetData())->GetGrid()/*,guessed*/);
  unfolding.SetMaxNumberOfIterations(100);
  unfolding.SetMaxChi2PerDOF(1.e-07);
  unfolding.UseSmoothing(fit2D,"ren");
  //unfolding.UseSmoothing();
  
  b.Stop("init");
  b.Start("unfolding");
  unfolding.Unfold();
  b.Stop("unfolding");
  b.Start("finish");

  canvas->cd();

  TH2D* h_gen = generated->Project(0,1);
  h_gen->SetTitle("generated");
  h_gen->Draw("lego2");
  canvas->SaveAs("/tmp/generated.eps");
  canvas->SaveAs("/tmp/generated.gif");

  TH2D* h_meas = measured->Project(0,1);
  h_meas->SetTitle("measured");
  h_meas->Draw("lego2");
  canvas->SaveAs("/tmp/measured.eps");
  canvas->SaveAs("/tmp/measured.gif");

  TH2D* h_guessed = unfolding.GetOriginalPrior()->Projection(1,0);
  h_guessed->SetTitle("a priori");
  h_guessed->Draw("lego2");
  canvas->SaveAs("/tmp/apriori.eps");
  canvas->SaveAs("/tmp/apriori.gif");

  TH2D* h_eff = efficiency->Project(0,1);
  h_eff->SetTitle("efficiency");
  h_eff->Draw("e");

  TH2D* h_unf = unfolding.GetUnfolded()->Projection(1,0);
  h_unf->SetTitle("unfolded");
  h_unf->Draw("lego2");
  fit2D->Draw("surf same");
  return;
  canvas->SaveAs("/tmp/unfolded.eps");
  canvas->SaveAs("/tmp/unfolded.gif");


  TH2D* ratio = (TH2D*)h_unf->Clone();
  ratio->SetTitle("unfolded / generated");
  ratio->Divide(h_unf,h_gen,1,1);
  ratio->GetZaxis()->SetRangeUser(0.5,1.5);
  //ratio->Draw("cont4z");
  ratio->Draw("lego2");
  //ratio->Draw("e");
//   ratio->ProjectionY()->Draw();
  canvas->SaveAs("/tmp/ratio.eps");
  canvas->SaveAs("/tmp/ratio.gif");
   

  canvas->cd(7);
  TH2* orig = unfolding.GetOriginalPrior()->Projection(1,0);
  orig->SetTitle("original prior");
  orig->Draw("lego2");

  canvas->cd(8);
  TH2D* h_est = (TH2D*) unfolding.GetEstMeasured()->Projection(1,0);
  h_est->SetTitle("est. measured");
  h_est->Draw("e");

  canvas->cd(9);
  unfolding.GetUnfolded()->Projection(0)->Draw();

  canvas->cd(10);
  unfolding.GetUnfolded()->Projection(1)->Draw();

  canvas->cd(11);
  h_gen->ProjectionX()->Draw();

  canvas->cd(12);
  h_gen->ProjectionY()->Draw();

  b.Stop("finish");
  Float_t x,y;
  b.Summary(x,y);
}

// ====================================================================


void Load(){
  gSystem->Load("libANALYSIS");
  gSystem->Load("libCORRFW");
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  TCanvas * canvas = new TCanvas("canvas","",600,400);
  canvas->SetFillColor(10);
  canvas->SetFrameFillColor(10);
}

TObject* GetMeasuredSpectrum() {
  TFile * f = TFile::Open("test/output.root","read");
  AliCFContainer* c = (AliCFContainer*)f->Get("container");
  AliCFDataGrid* data = new AliCFDataGrid("data","",*c);
  data->SetMeasured(1);
  return data;
}
TObject* GetGeneratedSpectrum() {
  TFile * f = TFile::Open("test/output.root","read");
  AliCFContainer* c = (AliCFContainer*)f->Get("container");
  AliCFDataGrid* data = new AliCFDataGrid("data","",*c);
  data->SetMeasured(0);
  return data;
}
TObject* GetEfficiency() {
  TFile * f = TFile::Open("test/output.root","read");
  AliCFContainer* c = (AliCFContainer*)f->Get("container");
  AliCFEffGrid* eff = new AliCFEffGrid("eff","",*c);
  eff->CalculateEfficiency(2,0);
  return eff;
}
THnSparse* GetResponseMatrix() {
  TFile * f = TFile::Open("test/output.root","read");
  THnSparse* h = (THnSparse*)f->Get("correlation");
  return h;
}
THnSparse* CreateGuessed(const THnSparse* h) {
  THnSparse* guessed = (THnSparse*) h->Clone();
  return guessed ;
}
