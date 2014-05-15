Float_t gpTMin_T1 = 6.0;
Float_t gpTMax_T1 = 14.0;
Float_t gpTMin_T2 = 4.0;
Float_t gpTMax_T2 = 10.0;
Float_t gpTMin_assoc = 1.0;
Float_t gpTMax_assoc = 6.0;
Float_t gZVtxRange = -1;
Float_t gAxis = 4;

void loadlibs()
{
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCORRFW");
  gSystem->Load("libPWGTools");
  gSystem->Load("libPWGCFCorrelationsBase");
}

void* events = 0;

void* GetTwoPlusOne(const char* fileName, TList** listRef = 0, Bool_t mixed = kFALSE, const char* tag = "")
{
    file = TFile::Open(fileName);
    if (!file)
      return 0;
 
    list = (TList*) gFile->Get("PWGCF_TwoPlusOne/histosTwoPlusOne");
      
    if (!list)
      return 0;
      
    if (listRef)
      *listRef = list;
      
    return events = list->FindObject("AliTwoPlusOneContainer");
}

void PlotQA(const char* fileName, const char* tag = "")
{
  loadlibs();

  TFile::Open(fileName);

  AliTwoPlusOneContainer* h = (AliTwoPlusOneContainer*) GetTwoPlusOne(fileName, 0, kFALSE, tag);

  if (h->GetData()->GetTrackHist(0)->GetGrid(6)->GetGrid()->GetNbins() == 0)
  {
    Printf("We have %d axes", ((AliTHn*) h->GetData()->GetTrackHist(0)->GetNVar()));
    
    ((AliTHn*) h->GetData()->GetTrackHist(0))->FillParent();
    ((AliTHn*) h->GetData()->GetTrackHist(0))->DeleteContainers();
  }

  TCanvas* c1 = new TCanvas("can1", "can1", 1200, 800);
  c1->Divide(2, 1);

  c1->cd(1);
  AliCFGridSparse* near_plot = h->GetData()->GetTrackHist(0)->GetGrid((AliUEHist::CFStep) AliTwoPlusOneContainer::kSameNS);
  //near_plot->SetRangeUser(2, gpTMin_T1, gpTMax_T1);
  //near_plot->SetRangeUser(6, gpTMin_T2, gpTMax_T2);
  //near_plot->SetRangeUser(1, gpTMin_assoc, gpTMax_assoc);
  TH1D* tracks_near = near_plot->Project(gAxis);

  tracks_near->DrawCopy();
  c1->cd(2);
 
  Printf("\n output_TwoPlusOne: entries %f \n", h->GetData()->GetTrackHist(AliUEHist::kToward)->GetGrid((AliUEHist::CFStep) AliTwoPlusOneContainer::kSameNS)->Project(1)->GetEntries());

  AliCFGridSparse* away_plot = h->GetData()->GetTrackHist(0)->GetGrid((AliUEHist::CFStep) AliTwoPlusOneContainer::kMixedNS);
  //away_plot->SetRangeUser(2, gpTMin_T1, gpTMax_T1);
  //away_plot->SetRangeUser(6, gpTMin_T2, gpTMax_T2);
  //away_plot->SetRangeUser(1, gpTMin_assoc, gpTMax_assoc);
  TH1D* tracks_away = away_plot->Project(gAxis);
  tracks_away->DrawCopy();

  AliCFGridSparse* trigger = (AliCFGridSparse*) h->GetData()->GetEventHist()->GetGrid((AliUEHist::CFStep) AliTwoPlusOneContainer::kSameNS);
  //trigger->SetRangeUser(0, gpTMin_T1, gpTMax_T1);

  TCanvas* c2 = new TCanvas("can2", "can2", 900, 600);
  TH1D* trigger_hist = trigger->Project(0);
  trigger_hist->DrawCopy();


  TCanvas* c3 = new TCanvas("can3", "can3", 900, 600);
  TH2* h2_near = h->GetData()->GetSumOfRatios2(h->GetData(), (AliUEHist::CFStep) AliTwoPlusOneContainer::kSameNS, 0, 4.0, 14.0, 0, 100, kFALSE, (AliUEHist::CFStep) AliTwoPlusOneContainer::kMixedNS);

  h2_near->DrawCopy("surf1");
}




