void Load();
TList *GetResults(const Char_t *testfile);
TObject* GetMeasuredSpectrum(TList *fContainer);
TObject* GetGeneratedSpectrum(TList *fContainer);
TObject* GetEfficiency(TList *fContainer);
THnSparse* GetResponseMatrix(TList *fContainer);
THnSparse* CreateGuessed(const THnSparse* h);

void testUnfolding(const char *testfile = "HFEtask.root") {
  Load();
  TList *results = GetResults(testfile);
  if(!results){
    Error("No output objects: Calculation will terminate here");
    return;
  }
  // get the essential
  AliCFDataGrid *measured    = (AliCFDataGrid*) GetMeasuredSpectrum(testfile);
  AliCFDataGrid *generated   = (AliCFDataGrid*) GetGeneratedSpectrum(testfile);
  AliCFEffGrid  *efficiency  = (AliCFEffGrid*)  GetEfficiency(testfile);
  THnSparse     *response    = (THnSparse*)     GetResponseMatrix(testfile);
  
  // create a guessed "a priori" distribution using binning of MC
  THnSparse* guessed = CreateGuessed(((AliCFGridSparse*)generated->GetData())->GetGrid()) ;
  //----
  AliCFUnfolding unfolding("unfolding","",3,response,efficiency->GetGrid(),((AliCFGridSparse*)measured->GetData())->GetGrid(),guessed);
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
  //printf("c1\n");
  h_gen->Draw("lego2");

  canvas->cd(2);
  TH2D* h_meas = measured->Project(0,1);
  h_meas->SetTitle("measured");
  //printf("c2\n");
  h_meas->Draw("lego2");
  
  canvas->cd(3);
  TH2D* h_guessed = guessed->Projection(1,0);
  h_guessed->SetTitle("a priori");
  //printf("c3\n");
  h_guessed->Draw("lego2");

  canvas->cd(4);
  TH2D* h_eff = efficiency->Project(0,1);
  h_eff->SetTitle("efficiency");
  //printf("c4\n");
  h_eff->Draw("lego2");

  canvas->cd(5);
  TH2D* h_unf = result->Projection(1,0);
  h_unf->SetTitle("unfolded");
  //printf("c5\n");
  h_unf->Draw("lego2");

  canvas->cd(6);
  TH2D* ratio = (TH2D*)h_unf->Clone();
  ratio->SetTitle("unfolded / generated");
  ratio->Divide(h_unf,h_gen,1,1);
//   ratio->Draw("cont4z");
//   ratio->Draw("surf2");
  //printf("c6\n");
  ratio->Draw("lego2");

  canvas->cd(7);
  TH2* orig = unfolding.GetOriginalPrior()->Projection(1,0);
  orig->SetTitle("original prior");
  //printf("c7\n");
  orig->Draw("lego2");

  canvas->cd(8);
  AliCFDataGrid* corrected = (AliCFDataGrid*)measured->Clone();
  //corrected->ApplyEffCorrection(*efficiency);
  TH2D* corr = corrected->Project(0,1);
  corr->SetTitle("simple correction");
  //printf("c8\n");
  corr->Draw("lego2");

  canvas->cd(9);
  TH2D* ratio2 = (TH2D*) corr->Clone();
  ratio2->Divide(corr,h_gen,1,1);
  ratio2->SetTitle("simple correction / generated");
  //printf("c9\n");
  ratio2->Draw("lego2");

  return;
}

// ====================================================================


void Load(){
  gSystem->Load("libANALYSIS");
  gSystem->Load("libCORRFW");
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
}

TList *GetResults(const Char_t *testfile){
  //
  // read output
  //
  TFile *f = TFile::Open(testfile);
  if(!f || f->IsZombie()){
    Error("File not readable");
    return 0x0;
  }
  TList *l = dynamic_cast<TList *>(f->Get("Results"));
  if(l){
    Error("Output container not found");
    f->Close(); delete f;
    return 0x0;
  } 
  TList *returnlist = dynamic_cast<TList *>(l->Clone());
  f->Close(); delete f;
  return returnlist;
}

TObject* GetMeasuredSpectrum(TList *fContainer) {
  AliCFContainer* c = (AliCFContainer*)fContainer->FindObject("container");
  AliCFDataGrid* data = new AliCFDataGrid("data","",*c);
  //data->SetMeasured(AliHFEcuts::kStepHFEcuts + 1);
  data->SetMeasured(5);
  return data;
}

TObject* GetGeneratedSpectrum(TList *fContainer) {
  AliCFContainer* c = (AliCFContainer*)fContainer->FindObject("container");
  AliCFDataGrid* data = new AliCFDataGrid("data","",*c);
  //data->SetMeasured(AliHFEcuts::kStepMCGenerated);
  data->SetMeasured(0);
  return data;
}

TObject* GetEfficiency(TList *fContainer) {
  AliCFContainer* c = (AliCFContainer*)fContainer->FindObject("container");
  AliCFEffGrid* eff = new AliCFEffGrid("eff","",*c);
  //eff->CalculateEfficiency(AliHFEcuts::kStepHFEcuts + 2,AliHFEcuts::kStepMCGenerated);
  eff->CalculateEfficiency(6,0);
  return eff;
}

THnSparse* GetResponseMatrix(TList *fContainer) {
  THnSparse* h = (THnSparse*)f->Get("correlation");
  return h;
}

THnSparse* CreateGuessed(const THnSparse* h) {
  THnSparse* guessed = (THnSparse*) h->Clone();
  return guessed ;
}
