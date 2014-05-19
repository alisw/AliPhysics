Float_t gpTMin_T1 = 6.0;
Float_t gpTMax_T1 = 14.0;
Float_t gpTMin_T2 = 4.0;
Float_t gpTMax_T2 = 10.0;
Float_t gpTMin_assoc = 3.0;
Float_t gpTMax_assoc = 6.0;
Float_t gZVtxRange = -1;
Float_t gAxis = 0;
Float_t gBasisSize = 350;

void loadlibs()
{
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCORRFW");
  gSystem->Load("libPWGTools");
  gSystem->Load("libPWGCFCorrelationsBase");
}

void* getList(const char* fileName, const char* folder)
{
    file = TFile::Open(fileName);
    if (!file)
      return 0;
 
    list = (TList*) gFile->Get(folder);
    //list = (TList*) gFile->Get("PWGCF_TwoPlusOne/histosTwoPlusOne");
      
    if (!list)
      return 0;

    return list;
}

void* GetTwoPlusOne(const char* fileName)
{
  list = (TList*) getList(fileName, "PWGCF_TwoPlusOne/histosTwoPlusOne");
  
  AliTwoPlusOneContainer* container = (AliTwoPlusOneContainer*) list->FindObject("AliTwoPlusOneContainer");

  if (container->GetData()->GetTrackHist(0)->GetGrid(6)->GetGrid()->GetNbins() == 0)
  {
    Printf("We have %d axes", ((AliTHn*) container->GetData()->GetTrackHist(0)->GetNVar()));
    
    ((AliTHn*) container->GetData()->GetTrackHist(0))->FillParent();
    ((AliTHn*) container->GetData()->GetTrackHist(0))->DeleteContainers();
  }

  return container;
}

void* GetPhiCorrelations(const char* fileName)
{
  list = (TList*) getList(fileName, "PWG4_PhiCorrelations/histosPhiCorrelations");
  

  AliUEHistograms* container = list->FindObject("AliUEHistogramsSame");

  if (container->GetUEHist(2)->GetTrackHist(0)->GetGrid(6)->GetGrid()->GetNbins() == 0)
  {
    Printf("We have %d axes", ((AliTHn*) container->GetUEHist(2)->GetTrackHist(0)->GetNVar()));
    
    ((AliTHn*) container->GetUEHist(2)->GetTrackHist(0))->FillParent();
    ((AliTHn*) container->GetUEHist(2)->GetTrackHist(0))->DeleteContainers();
  }

  return container;
}

void* showEvents(const char* fileName, int twoPlusOne = 1)
{
  if(twoPlusOne)
    list = (TList*) getList(fileName, "PWGCF_TwoPlusOne/histosTwoPlusOne");
  else
    list = (TList*) getList(fileName, "PWG4_PhiCorrelations/histosPhiCorrelations");

  TH1F* event_stat = (TH1F*) list->FindObject("eventStat");
  TCanvas* can = new TCanvas();
  event_stat->DrawCopy();
}


//returns maximum of the flow
void subtractFlow(TH2* etaPhi){
  int firstbin = etaPhi->GetYaxis()->FindBin(1.01);
  int lastbin = etaPhi->GetYaxis()->FindBin(1.59);

  TH1D* etaPhi_proj = etaPhi->ProjectionX("_px", firstbin, lastbin);
  int usedBins = lastbin - firstbin + 1;

  firstbin = etaPhi->GetYaxis()->FindBin(-1.59);
  lastbin = etaPhi->GetYaxis()->FindBin(-1.01);
  TH1D* signal2 = etaPhi->ProjectionX("_px2", firstbin, lastbin);
  usedBins += lastbin - firstbin + 1;
  etaPhi_proj->Add(signal2, 1.0);

  etaPhi_proj->Scale(1.0/usedBins);

  etaPhi->GetYaxis()->SetRangeUser(-1.0, 1.0);

  for(int i=0;i<=etaPhi->GetNbinsX();i++){
    double subtract = etaPhi_proj->GetBinContent(i);

    for(int j=0;j<etaPhi->GetNbinsY();j++){
      double content = etaPhi->GetBinContent(i,j)-subtract;
      etaPhi->SetBinContent(i,j,content);
    }
  }

}


void compareTrigger1(const char* fileName){
  loadlibs();

  double pt1Min = 6.0;
  double pt1Max = 8.0;

  double ptAssocMin = 6.0;
  double ptAssocMax = 14.0;

  AliTwoPlusOneContainer* twoPlusOne = (AliTwoPlusOneContainer*) GetTwoPlusOne(fileName);
  
  AliUEHistograms* phiCorr = (AliUEHistograms*) GetPhiCorrelations(fileName);

  THnBase* trackSameAll_1 = 0;
  THnBase* trackSameAll_2 = 0;

  TH2* eventSameAll_1 = 0;
  TH2* eventSameAll_2 = 0;

  //twoPlusOne->GetData()->SetPtRange(ptAssocMin, ptAssocMax);
  //phiCorr->GetUEHist(2)->SetPtRange(ptAssocMin, ptAssocMax);

  twoPlusOne->GetData()->GetHistsZVtxMult((AliUEHist::CFStep) AliTwoPlusOneContainer::kSameNS, AliUEHist::kToward, pt1Min, pt1Max, &trackSameAll_1, &eventSameAll_1);

  phiCorr->GetUEHist(2)->GetHistsZVtxMult(6, AliUEHist::kToward, pt1Min, pt1Max, &trackSameAll_2, &eventSameAll_2);
  
  int projection = 1;
  
  AliCFGridSparse* trigger = (AliCFGridSparse*) twoPlusOne->GetData()->GetTrackHist(0)->GetGrid((AliUEHist::CFStep) AliTwoPlusOneContainer::kSameNS);
  TH1D* trigger_hist = trigger->Project(projection);
  //trigger_hist->Rebin(6);

  TCanvas* c1 =  new TCanvas();
  cout<<"entries twoplusone "<<trigger_hist->GetEntries()<<endl;
  trigger_hist->DrawCopy();

  AliCFGridSparse* trigger_2 = (AliCFGridSparse*) phiCorr->GetUEHist(2)->GetTrackHist(0)->GetGrid(6);
  TH1D* trigger_hist_2 = trigger_2->Project(projection);
  //trigger_hist_2->Rebin(6);

  //TCanvas* c2 = new TCanvas("can2", "can2", 900, 600);
  trigger_hist_2->SetLineColor(kRed);
  cout<<"entries phicorr "<<trigger_hist_2->GetEntries()<<endl;
  trigger_hist_2->DrawCopy("same");
  

  /*
  AliCFGridSparse* trigger = (AliCFGridSparse*) twoPlusOne->GetData()->GetEventHist()->GetGrid((AliUEHist::CFStep) AliTwoPlusOneContainer::kSameNS);
  //trigger->SetRangeUser(0, gpTMin_T1, gpTMax_T1);

  TCanvas* c1 = new TCanvas("can1", "can1", 900, 600);
  TH1D* trigger_hist = trigger->Project(0);
  trigger_hist->DrawCopy();

  AliCFGridSparse* trigger_2 = (AliCFGridSparse*) phiCorr->GetUEHist(2)->GetEventHist()->GetGrid(6);
  //trigger->SetRangeUser(0, gpTMin_T1, gpTMax_T1);

  TCanvas* c2 = new TCanvas("can2", "can2", 900, 600);
  TH1D* trigger_hist_2 = trigger_2->Project(0);
  trigger_hist_2->DrawCopy();
  */
  //TCanvas* can2 = new TCanvas("can2", "can2");
  //tracksSame_2->DrawCopy("colz");
}


void showAllTriggerPt(const char* fileName, Int_t multBinBegin = 0, Int_t multBinEnd = 16, Int_t side = 0){
  
  const Int_t pt1_bins_length = 3;
  //Int_t pt1_bins[pt1_bins_length+1] = {6.0, 8.0, 10.0, 12.0, 14.0};
  Int_t pt1_bins[pt1_bins_length+1] = {6.0, 8.0, 10.0, 12.0};
  const Int_t pt2_bins_length = 2;
  Int_t pt2_bins[pt2_bins_length] = {4.0, 6.0};

  for(Int_t i=0; i<pt1_bins_length; i++){
    Double_t pt1Minimum = pt1_bins[i];
    Double_t pt1Maximum = pt1_bins[i+1];
    for(Int_t j=0; j<pt2_bins_length; j++){
      Double_t pt2Minimum = pt2_bins[j];
      getAnalysis(fileName, pt1Minimum, pt1Maximum, pt2Minimum, gpTMin_assoc, gpTMax_assoc, multBinBegin, multBinEnd, side, i, j);
    }
  }

}

void getAnalysis(const char* fileName, double pt1Min = 4.0, double pt1Max = 14.0, Double_t pt2Min = 2.0, Double_t ptAssocMin = 0.5, Double_t ptAssocMax = 8.0, Int_t multBinBegin = 1, Int_t multBinEnd = 5, Int_t step = 0, Int_t posX = 1, Int_t posY = 1)
{
  loadlibs();
  TFile::Open(fileName);
  
  AliTwoPlusOneContainer* h = (AliTwoPlusOneContainer*) GetTwoPlusOne(fileName);

  AliUEHist::CFStep step_same = step;//(AliUEHist::CFStep) AliTwoPlusOneContainer::kSameNS;
  AliUEHist::CFStep step_mixed = step+2;//(AliUEHist::CFStep) AliTwoPlusOneContainer::kMixedNS;

  //GetSumOfRatios2(mixed, step, ptLeadMin, ptLeadMax, multBinBegin, multBinEnd, normalizePerTrigger, stepForMixed)
  h->GetData()->SetPtRange(ptAssocMin, ptAssocMax);
  h->GetData()->SetPt2Min(pt2Min);
  TH2* h2_etaPhi = h->GetData()->GetSumOfRatios2(h->GetData(), step_same, 0, pt1Min, pt1Max, multBinBegin, multBinEnd, kTRUE, step_mixed);

  subtractFlow(h2_etaPhi);
  //h2_etaPhi->GetXaxis()->SetRangeUser(-TMath::Pi()/2, TMath::Pi()/2);
  //h2_etaPhi->GetYaxis()->SetRangeUser(-1.6, 1.6);//corresponds to the really used era range

  TCanvas* c1 = new TCanvas(Form("can %i %i", posX, posY), Form("can %i %i", posX, posY), posX*gBasisSize+50, posY*gBasisSize+50, gBasisSize, gBasisSize);

  //if(h2_etaPhi)
  //h2_etaPhi->DrawCopy("surf1");

  int firstbin = h2_etaPhi->GetYaxis()->FindBin(-0.19);
  int lastbin = h2_etaPhi->GetYaxis()->FindBin(0.19);

  //int firstbin = -1;
  //int lastbin = -1;

  TH1D* h1_phi = h2_etaPhi->ProjectionX(Form("px_%i_%i", posX, posY), firstbin, lastbin);
  h1_phi->SetYTitle("1/N_{trig} \\  dN_{assoc}/d \\varphi");
  h1_phi->SetTitle("");
  h1_phi->SetStats(kFALSE);
  h1_phi->GetXaxis()->SetRangeUser(-TMath::Pi()/2, TMath::Pi()/2);
  h1_phi->DrawCopy();
}

//do correction which is done in GetSumOfRatios2 by hand
void getRawAnalysis(const char* fileName, double pt1Min = 4.0, double pt1Max = 14.0, Double_t ptAssocMin = 0.5, Double_t ptAssocMax = 8.0, Int_t posX = 1, Int_t posY = 1)
{
  loadlibs();

  AliTwoPlusOneContainer* twoPlusOne = (AliTwoPlusOneContainer*) GetTwoPlusOne(fileName);
  
  TH2* totalTracks = 0;
  TH2* totalMixed =0;

  THnBase* trackSameAll = 0;
  THnBase* trackMixedAll = 0;

  TH2* eventSameAll = 0;
  TH2* eventMixedAll = 0;

  twoPlusOne->GetData()->SetPtRange(ptAssocMin, ptAssocMax);

  twoPlusOne->GetData()->GetHistsZVtxMult((AliUEHist::CFStep) AliTwoPlusOneContainer::kSameAS, AliUEHist::kToward, pt1Min, pt1Max, &trackSameAll, &eventSameAll);
  twoPlusOne->GetData()->GetHistsZVtxMult((AliUEHist::CFStep) AliTwoPlusOneContainer::kMixedAS, AliUEHist::kToward, pt1Min, pt1Max, &trackMixedAll, &eventMixedAll);

  //const Int_t num_mult_bins = 4;
  //Int_t mult_bins[num_mult_bins+1] = {0, 5, 6, 10, 15};
  //0, 1, 2, 3, 4, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100.1 //standard mult bins
  
  //const Int_t num_vtx_bins = 3;
  //Int_t vtx_bins[num_vtx_bins+1] = {0, 3, 4, 7};
  //-7, -5, -3, -1, 1, 3, 5, 7 //standard vertex bins

  const Int_t num_mult_bins = 1;
  Int_t mult_bins[num_mult_bins+1] = {1, 16};
  const Int_t num_vtx_bins = 1;
  Int_t vtx_bins[num_vtx_bins+1] = {1, 7};

  Int_t trigger_same = 0;
  Int_t trigger_mixed = 0;

  for(Int_t multBin=0; multBin<num_mult_bins; multBin++){
    trackSameAll->GetAxis(3)->SetRange(mult_bins[multBin], mult_bins[multBin+1]);
    trackMixedAll->GetAxis(3)->SetRange(mult_bins[multBin], mult_bins[multBin+1]);

    Printf("mult Bin start %i, bin center %f", mult_bins[multBin], trackSameAll->GetAxis(3)->GetBinCenter(mult_bins[multBin]));
    Printf("mult Bin end %i, bin center %f", mult_bins[multBin+1], trackSameAll->GetAxis(3)->GetBinCenter(mult_bins[multBin+1]));

    for(Int_t vertexBin = 0; vertexBin<num_vtx_bins; vertexBin++){
      trackSameAll->GetAxis(2)->SetRange(vtx_bins[vertexBin], vtx_bins[vertexBin+1]);
      trackMixedAll->GetAxis(2)->SetRange(vtx_bins[vertexBin], vtx_bins[vertexBin+1]);

      Printf("vertex Bin start %i, bin center %f", vtx_bins[vertexBin], trackSameAll->GetAxis(2)->GetBinCenter(vtx_bins[vertexBin]));
      Printf("vertex Bin end %i, bin center %f", vtx_bins[vertexBin+1], trackSameAll->GetAxis(2)->GetBinCenter(vtx_bins[vertexBin+1]));

      trigger_same += eventSameAll->Integral(vtx_bins[vertexBin], vtx_bins[vertexBin+1], mult_bins[multBin], mult_bins[multBin+1]);
      trigger_mixed += eventMixedAll->Integral(vtx_bins[vertexBin], vtx_bins[vertexBin+1], mult_bins[multBin], mult_bins[multBin+1]);

      TH2* tracksSame = trackSameAll->Projection(1, 0, "E");
      TH2* tracksMixed = trackMixedAll->Projection(1, 0, "E");

      tracksSame->Rebin2D(2, 2);
      tracksMixed->Rebin2D(2, 2);

      //tracksSame->Divide(tracksMixed);

      if (!totalTracks)
	totalTracks = (TH2*) tracksSame->Clone("totalTracks");
      else
	totalTracks->Add(tracksSame);

      if (!totalMixed)
	totalMixed = (TH2*) tracksMixed->Clone("totalTracks");
      else
	totalMixed->Add(tracksMixed);

      delete tracksSame;
      delete tracksMixed;
    }
  }

  //totalTracks->Divide(totalMixed);

  TCanvas* c1 = new TCanvas(Form("can %i %i", posX, posY), Form("can %i %i", posX, posY), posX*gBasisSize+50, posY*gBasisSize+50, gBasisSize, gBasisSize);
  //if(totalTracks)
  //totalTracks->DrawCopy("surf1");
  TH1D* etaTracks = totalTracks->ProjectionY();
  etaTracks->Scale(1.0/etaTracks->Integral());
  etaTracks->DrawCopy();

  Printf(Form("trigger same side %i",trigger_same));
  Printf(Form("trigger mixed side %i",trigger_mixed));

  posX++;

  //TCanvas* c2 = new TCanvas(Form("can %i %i", posX, posY), Form("can %i %i", posX, posY), posX*gBasisSize+50, posY*gBasisSize+50, gBasisSize, gBasisSize);
  //if(totalMixed)
  //totalMixed->DrawCopy("surf1");

    TH1D* etaMixed = totalMixed->ProjectionY();
    if(etaMixed->Integral()>0)
      etaMixed->Scale(1.0/etaMixed->Integral());
    etaMixed->SetLineColor(kRed);
    etaMixed->DrawCopy("same");

}



//macro to see if the analysis worked
//the near side of the same event and the near side of the mixed event are shown
//the number of trigger particles of the same event (near side) are shown
//GetSumOfRatios2 is tested, it shows the near side of the same event divided by the same side of the mixed event. This is done per multiplicity and z vtx bin
void PlotQA(const char* fileName)
{
  loadlibs();

  TFile::Open(fileName);

  AliTwoPlusOneContainer* h = (AliTwoPlusOneContainer*) GetTwoPlusOne(fileName);
  
  TCanvas* c1 = new TCanvas("can1", "can1", 1200, 800);
  c1->Divide(3, 1);
  
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
  
  c1->cd(3);

  AliCFGridSparse* mixed_comb_plot = h->GetData()->GetTrackHist(0)->GetGrid((AliUEHist::CFStep) AliTwoPlusOneContainer::kMixedCombNS);
  //away_plot->SetRangeUser(2, gpTMin_T1, gpTMax_T1);
  //away_plot->SetRangeUser(6, gpTMin_T2, gpTMax_T2);
  //away_plot->SetRangeUser(1, gpTMin_assoc, gpTMax_assoc);
  TH1D* tracks_mixed_comb = mixed_comb_plot->Project(gAxis);
  tracks_mixed_comb->DrawCopy();

  
  AliCFGridSparse* trigger = (AliCFGridSparse*) h->GetData()->GetEventHist()->GetGrid((AliUEHist::CFStep) AliTwoPlusOneContainer::kSameNS);
  //trigger->SetRangeUser(0, gpTMin_T1, gpTMax_T1);

  TCanvas* c2 = new TCanvas("can2", "can2", 900, 600);
  TH1D* trigger_hist = trigger->Project(0);
  trigger_hist->DrawCopy();
  
  /*
  TCanvas* c3 = new TCanvas("can3", "can3", 900, 600);
  TH2* h2_near = h->GetData()->GetSumOfRatios2(h->GetData(), (AliUEHist::CFStep) AliTwoPlusOneContainer::kSameNS, 0, 4.0, 14.0, 0, 100, kFALSE, (AliUEHist::CFStep) AliTwoPlusOneContainer::kMixedNS);
  h2_near->GetYaxis()->SetRangeUser(-1.8, 1.8);`

  if(h2_near)
    h2_near->DrawCopy("surf1");
  */
}

