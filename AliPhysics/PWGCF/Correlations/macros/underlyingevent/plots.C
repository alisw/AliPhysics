void loadlibs()
{
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCORRFW");
  gSystem->Load("libJETAN");
  gSystem->Load("libPWG4JetTasks");
}

void PlotTrackingEfficiencies()
{
  TGrid::Connect("alien://");

  loadlibs();

  ifstream in;
  in.open("list");

  Int_t i=0;

  legend = new TLegend(0.5, 0.1, 0.9, 0.5);
  legend->SetNColumns(2);

  TString line;
  while (in.good())
  {
    in >> line;

    if (line.Length() == 0)
      continue;

    TString fileName;
    //fileName.Form("%s/%s/PWG4_JetTasksOutput.root", "alien:///alice/cern.ch/user/k/kleinb/analysisESD/LHC10d4/output_pwg4train_LHC10d4_101001a", line.Data());
    fileName.Form("%s/%s/PWG4_JetTasksOutput.root", "maps", line.Data());
    Printf("%s", fileName.Data());
    file = TFile::Open(fileName);
    if (!file)
      continue;
    list = (TList*) file->Get("PWG4_LeadingTrackUE/histosLeadingTrackUE");
    AliUEHistograms* h = (AliUEHistograms*) list->FindObject("AliUEHistograms");
    h->SetEtaRange(-0.79, 0.79);

    TH1::AddDirectory(0);
    eff = (TH1*) h->GetNumberDensitypT()->GetTrackingEfficiency(1)->Clone(Form("trackingefficiency%s", line.Data()));

    delete h;
    //delete list;
    file->Close();

    eff->SetLineColor(i+1);
    eff->GetXaxis()->SetRangeUser(0, 10);
    eff->GetYaxis()->SetRangeUser(0.6, 0.8);
    eff->Draw((i == 0) ? "" : "SAME");
 
    legend->AddEntry(eff, line);

    i++;
    //if (i > 1)
    //  break;
  }

  legend->Draw();

  in.close();
}

void TrackingEfficiencySpecies(const char* fileName = "PWG4_JetTasksOutput.root")
{
  loadlibs();

  file = TFile::Open(fileName);
  list = (TList*) file->Get("PWG4_LeadingTrackUE/histosLeadingTrackUE");
  AliUEHistograms* h = (AliUEHistograms*) list->FindObject("AliUEHistograms");
  h->SetEtaRange(-0.79, 0.79);

  eff = (TH2*) h->GetNumberDensitypT()->GetTrackEfficiency(AliUEHist::kCFStepAnaTopology, AliUEHist::kCFStepTrackedOnlyPrim, 1, 2);
  eff->Draw("colz");

  for (Int_t i=0; i<4; i++)
  {
    proj = eff->ProjectionX(Form("p%d", i), i+1, i+1);
    proj->SetLineColor(i+1);
    proj->Draw((i==0) ? "" : "SAME");
  }
}

void CheckTrackingEfficiency(const char* reference = 0, const char* fileName = "PWG4_JetTasksOutput.root", Bool_t all = kFALSE)
{
  loadlibs();

  file = TFile::Open(fileName);
  list = (TList*) file->Get("PWG4_LeadingTrackUE/histosLeadingTrackUE");
  tree = (TTree*) list->FindObject("UEAnalysisSettings");
  
  if (reference)
  {
    file2 = TFile::Open(reference);
    refHist = (TH1D*) file2->Get("trackingefficiency");
  }
  else
  {
    AliUEHistograms* h = (AliUEHistograms*) list->FindObject("AliUEHistograms");
    h->SetEtaRange(-0.79, 0.79);
    refHist = h->GetNumberDensitypT()->GetTrackingEfficiency(1);
  }

  refHist->SetLineWidth(3);
  refHist->GetListOfFunctions()->Clear();
  refHist->Draw();

  TH1D* hist = 0;
  
  tree->SetBranchAddress("fkTrackingEfficiency", &hist);
  
  for (Int_t i=0; i<tree->GetEntries(); i++)
  {
    tree->GetEntry(i);
    hist->GetListOfFunctions()->Clear();
    hist->SetLineColor(i+2);
    hist->DrawCopy("SAME");
    if (!all)
      break;
    Printf("%d", i);
  }
}

void PlotSingleTrackingEfficiency(const char* fileName, Int_t what = 0)
{
  loadlibs();

  file = TFile::Open(fileName);
  list = (TList*) file->Get("PWG4_LeadingTrackUE/histosLeadingTrackUE");
  AliUEHistograms* h = (AliUEHistograms*) list->FindObject("AliUEHistograms");
  h->SetEtaRange(-0.79, 0.79);

  if (what == 0)
    eff = (TH2*) h->GetNumberDensitypT()->GetTrackEfficiency(AliUEHist::kCFStepAnaTopology, AliUEHist::kCFStepTrackedOnlyPrim, 0, 1);
  else
    eff = (TH2*) h->GetNumberDensitypT()->GetTrackingContamination();

  eff->Draw("colz");
}

void ExtendTrackingEfficiency(const char* fileName)
{
  loadlibs();

  file = TFile::Open(fileName);
  list = (TList*) file->Get("PWG4_LeadingTrackUE/histosLeadingTrackUE");
  AliUEHistograms* h = (AliUEHistograms*) list->FindObject("AliUEHistograms");
  h->SetEtaRange(-0.79, 0.79);

  h->GetUEHist(0)->ExtendTrackingEfficiency(1);
}

void PlotSystUncertainties()
{
  //
}
 
