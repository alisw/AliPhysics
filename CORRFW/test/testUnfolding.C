
void testUnfolding() {
  Load();

  // get the essential
  AliCFDataGrid *measured    = (AliCFDataGrid*) GetMeasuredSpectrum();
  AliCFDataGrid *generated   = (AliCFDataGrid*) GetGeneratedSpectrum();
  AliCFEffGrid  *efficiency  = (AliCFEffGrid*)  GetEfficiency();
  THnSparse     *response    = (THnSparse*)     GetResponseMatrix();
  
  // create a guessed "a priori" distribution using binning of MC
  THnSparse* guessed = CreateGuessed(((AliCFGridSparse*)generated->GetData())->GetGrid()) ;
  //----
  AliCFUnfolding unfolding("unfolding","",2,response,efficiency->GetGrid(),((AliCFGridSparse*)measured->GetData())->GetGrid(),guessed);
  unfolding.SetMaxNumberOfIterations(100);
  unfolding.SetMaxChi2PerDOF(1.e-07);
  //unfolding.UseSmoothing();
  unfolding.Unfold();

  THnSparse* result = unfolding.GetUnfolded();
  //----
  
  TCanvas * canvas = new TCanvas("canvas","",1000,700);
  canvas->Divide(3,3);

  canvas->cd(1);
  TH2D* h_gen = generated->Project(0,1);
  h_gen->SetTitle("generated");
  h_gen->Draw("lego2");

  canvas->cd(2);
  TH2D* h_meas = measured->Project(0,1);
  h_meas->SetTitle("measured");
  h_meas->Draw("lego2");
  
  canvas->cd(3);
  TH2D* h_guessed = guessed->Projection(1,0);
  h_guessed->SetTitle("a priori");
  h_guessed->Draw("lego2");

  canvas->cd(4);
  TH2D* h_eff = efficiency->Project(0,1);
  h_eff->SetTitle("efficiency");
  h_eff->Draw("lego2");

  canvas->cd(5);
  TH2D* h_unf = result->Projection(1,0);
  h_unf->SetTitle("unfolded");
  h_unf->Draw("lego2");

  canvas->cd(6);
  TH2D* ratio = (TH2D*)h_unf->Clone();
  ratio->SetTitle("unfolded / generated");
  ratio->Divide(h_unf,h_gen,1,1);
//   ratio->Draw("cont4z");
//   ratio->Draw("surf2");
  ratio->Draw("lego2");

  canvas->cd(7);
  TH2* orig = unfolding.GetOriginalPrior()->Projection(1,0);
  orig->SetTitle("original prior");
  orig->Draw("lego2");

  canvas->cd(8);
  AliCFDataGrid* corrected = (AliCFDataGrid*)measured->Clone();
  corrected->ApplyEffCorrection(*efficiency);
  TH2D* corr = corrected->Project(0,1);
  corr->SetTitle("simple correction");
  corr->Draw("lego2");

  canvas->cd(9);
  TH2D* ratio2 = (TH2D*) corr->Clone();
  ratio2->Divide(corr,h_gen,1,1);
  ratio2->SetTitle("simple correction / generated");
  ratio2->Draw("cont4z");

  return;
}

// ====================================================================


void Load(){
  gSystem->Load("libANALYSIS");
  gSystem->Load("libCORRFW");
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
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
