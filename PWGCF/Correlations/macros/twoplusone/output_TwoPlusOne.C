Float_t gpTMin_T1 = 6.0;
Float_t gpTMax_T1 = 14.0;
Float_t gpTMin_T2 = 4.0;
Float_t gpTMax_T2 = 10.0;
Float_t gpTMin_assoc = 1.0;
Float_t gpTMax_assoc = 4.0;
Float_t gVertex = 6.9;
Float_t gZVtxRange = -1;
Float_t gAxis = 0;
Float_t gBasisSize = 350;
Float_t g_phi_bin = 0.34;//0.85;//0.68;//0.32;//bins of 0.174532  //maximum is pi/4, because max(delta phi) + 2 alpha < pi/2 (because at Delta phi = pi/2 the number of background triggers are measured)
Float_t g_eta_bin = 0.39;//0.79;//0.59;//bins of 0.2

Float_t gSubtractFlow_start = 1.01;//1.01;//variable in eta
Float_t gSubtractFlow_end = 1.39;//1.39;

//Float_t gSubtractBaseline_start = 2.0;//0.699;//variable in phi
//Float_t gSubtractBaseline_end = 4.3;//1.396;
Float_t gSubtractBaseline_start = 1.0;//0.8;//1.0;//0.88;//0.699;//variable in phiss
Float_t gSubtractBaseline_end = 1.396;//1.5;//1.396;//1.396;

Double_t ptMaxAdd = 2; //pt2Max = pt2Min+4;

//char* path = "PWGCF_TwoPlusOne/histosTwoPlusOne";//before 1091
//char* path = "PWGCF_TwoPlusOne/addedEvents_";
//char* path = "PWGCF_TwoPlusOne/histosTwoPlusOne_lowPt";//1095
//char* path = "PWGCF_TwoPlusOne/histosTwoPlusOne_bgSame";//1181
//char* path = "PWGCF_TwoPlusOne/histosTwoPlusOne_leadingPt";//1269
//char* path = "PWGCF_TwoPlusOne/histosTwoPlusOne_triggerTwo";//1269 + 1278
//char* path = "PWGCF_TwoPlusOne/histosTwoPlusOne_highPt";//1326
//char* path = "PWGCF_TwoPlusOne/histosTwoPlusOne_vertexBin";//1614
//char* path = "PWGCF_TwoPlusOne/histosTwoPlusOne_globalCuts";//1662
//char* path = "PWGCF_TwoPlusOne/histosTwoPlusOne_STAR";//2076
//char* path = "PWGCF_TwoPlusOne/histosTwoPlusOne_highMixedStats";//2129
char* path = "PWGCF_TwoPlusOne/histosTwoPlusOne";//AR

//static const int pt_assoc_bins_number = 8;
//Double_t pt_assoc_bins[pt_assoc_bins_number] = {0.5, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0};
//Double_t pt_assoc_bin_center[pt_assoc_bins_number-1] = {0.75, 1.5, 2.5, 3.5, 5.0, 7.0, 9.0};
//Double_t pt_assoc_bin_error[pt_assoc_bins_number-1] = {0.25, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0};

//static const int pt_assoc_bins_number = 6;
//Double_t pt_assoc_bins[pt_assoc_bins_number] = {0.5, 1.0, 2.0, 3.0, 4.0, 6.0};
//Double_t pt_assoc_bin_center[pt_assoc_bins_number-1] = {0.75, 1.5, 2.5, 3.5, 5.0};
//Double_t pt_assoc_bin_error[pt_assoc_bins_number-1] = {0.25, 0.5, 0.5, 0.5, 1.0};

//static const int pt_assoc_bins_number = 8;
//Double_t pt_assoc_bins[pt_assoc_bins_number] = {0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0};
//Double_t pt_assoc_bin_center[pt_assoc_bins_number-1] = {0.75, 1.5, 2.5, 3.5, 4.5, 5.5, 7.0};
//Double_t pt_assoc_bin_error[pt_assoc_bins_number-1] = {0.25, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0};

static const int pt_assoc_bins_number = 7;
Double_t pt_assoc_bins[pt_assoc_bins_number] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0};
Double_t pt_assoc_bin_center[pt_assoc_bins_number-1] = {1.5, 2.5, 3.5, 4.5, 5.5, 7.0};
Double_t pt_assoc_bin_error[pt_assoc_bins_number-1] = {0.5, 0.5, 0.5, 0.5, 0.5, 1.0};

//if there are no bins otherwise, so that the code runs
//static const int pt_assoc_bins_number = 5;
//Double_t pt_assoc_bins[pt_assoc_bins_number] = {0.5, 1.0, 2.0, 3.0, 4.0};
//Double_t pt_assoc_bin_center[pt_assoc_bins_number-1] = {0.75, 1.5, 2.5, 3.5};
//Double_t pt_assoc_bin_error[pt_assoc_bins_number-1] = {0.25, 0.5, 0.5, 0.5};

//static const int pt_assoc_bins_number = 7;
//Double_t pt_assoc_bins[pt_assoc_bins_number] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0};
//Double_t pt_assoc_bin_center[pt_assoc_bins_number-1] = {1.5, 2.5, 3.5, 4.5, 5.5, 7.0};
//Double_t pt_assoc_bin_error[pt_assoc_bins_number-1] = {0.5, 0.5, 0.5, 0.5, 0.5, 1.0};

//static const int pt_assoc_bins_number = 6;
//Double_t pt_assoc_bins[pt_assoc_bins_number] = {2.0, 3.0, 4.0, 6.0, 8.0, 10};
//Double_t pt_assoc_bin_center[pt_assoc_bins_number-1] = {2.5, 3.5, 5.0, 7.0, 9.0};
//Double_t pt_assoc_bin_error[pt_assoc_bins_number-1] = {0.5, 0.5, 1.0, 1.0, 1.0};

//static const int pt_assoc_bins_number = 12;
//Double_t pt_assoc_bins[pt_assoc_bins_number] = {0.5, 0.75, 1.0, 1.25, 1.50, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0};
//Double_t pt_assoc_bin_center[pt_assoc_bins_number-1] = {0.625, 0.875, 1.125, 1.375, 1.75, 2.25, 2.75, 3.25, 3.75, 4.5, 5.5};
//Double_t pt_assoc_bin_error[pt_assoc_bins_number-1] = {0.125, 0.125, 0.125, 0.125, 0.25, 0.25, 0.25, 0.25, 0.25, 0.5, 0.5};

//static const int pt_assoc_bins_number = 10;
//Double_t pt_assoc_bins[pt_assoc_bins_number] = {0.5, 1.0, 1.50, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0};
//Double_t pt_assoc_bin_center[pt_assoc_bins_number-1] = {0.75, 1.25, 1.75, 2.25, 2.75, 3.25, 3.75, 4.5, 5.5};
//Double_t pt_assoc_bin_error[pt_assoc_bins_number-1] = {0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.5, 0.5};


//static const int pt_assoc_bins_number = 8;
//Double_t pt_assoc_bins[pt_assoc_bins_number] = {0.5, 1.0, 1.50, 2.0, 2.5, 3.0, 4.0, 6.0};
//Double_t pt_assoc_bin_center[pt_assoc_bins_number-1] = {0.75, 1.25, 1.75, 2.25, 2.75, 3.5, 5.0};
//Double_t pt_assoc_bin_error[pt_assoc_bins_number-1] = {0.25, 0.25, 0.25, 0.25, 0.25, 0.5, 1.0};


//static const int pt_assoc_bins_number = 5;
//Double_t pt_assoc_bins[pt_assoc_bins_number] = {0.5, 1.50, 2.5, 4.0, 6.0};
//Double_t pt_assoc_bin_center[pt_assoc_bins_number-1] = {1.0, 2.0, 3.25, 5.0};
//Double_t pt_assoc_bin_error[pt_assoc_bins_number-1] = {0.5, 0.5, 0.75, 1.0};


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
      
    if (!list)
      return 0;

    return list;
}

const char* lastFileName = 0;
const char* lastFileName2 = 0;
void* cacheEvent = 0;
void* cacheEvent2 = 0;


void* GetTwoPlusOne(const char* fileName)
{
  if (lastFileName && strcmp(lastFileName, fileName) == 0){
    Printf("GetTwoPlusOne --> Using cache for %s", fileName);
    return cacheEvent;
  }
  if (lastFileName2 && strcmp(lastFileName2, fileName) == 0){
    Printf("GetTwoPlusOne --> Using cache for %s", fileName);
    return cacheEvent2;
  }
  if(lastFileName)
    lastFileName2 = fileName;
  else
    lastFileName = fileName;

  //list = (TList*) getList(fileName, "PWGCF_TwoPlusOne/histosTwoPlusOne");
  list = (TList*) getList(fileName, path);
  //path = "PWGCF_TwoPlusOne/histosTwoPlusOne_highPt";
  
  AliTwoPlusOneContainer* container = (AliTwoPlusOneContainer*) list->FindObject("AliTwoPlusOneContainer");
  
  if (container->GetData()->GetTrackHist(0)->GetGrid(6)->GetGrid()->GetNbins() == 0)
  {
    list->Remove(container);
    Printf("We have %d axes", ((AliTHn*) container->GetData()->GetTrackHist(0)->GetNVar()));
    
    ((AliTHn*) container->GetData()->GetTrackHist(0))->FillParent();
    ((AliTHn*) container->GetData()->GetTrackHist(0))->DeleteContainers();

    list->Add(container);

    TString newFileName(fileName);
    newFileName.ReplaceAll("AnalysisResults_", "");
    newFileName = "AR_"+newFileName;
 
    file = TFile::Open(newFileName, "RECREATE");
    file->mkdir("PWGCF_TwoPlusOne");
    file->cd("PWGCF_TwoPlusOne");
    list->Write("histosTwoPlusOne", TObject::kSingleKey);
    file->Close();
  }

  
  if(lastFileName2 && strcmp(lastFileName2, fileName) == 0){
    cacheEvent2 = container;
    return cacheEvent2;
  }

  cacheEvent = container;

  return cacheEvent;
}

void* showEvents(const char* fileName, int twoPlusOne = 1)
{
  if(twoPlusOne)
    //list = (TList*) getList(fileName, "PWGCF_TwoPlusOne/histosTwoPlusOne");
    list = (TList*) getList(fileName, path);
  else
    list = (TList*) getList(fileName, "PWG4_PhiCorrelations/histosPhiCorrelations");

  TH1F* event_stat = (TH1F*) list->FindObject("eventStat");
  TCanvas* can = new TCanvas();
  event_stat->DrawCopy();
}

void* showEventsAtStage(const char* fileName)
{
  list = (TList*) getList(fileName, path);

  TH2F* event_stat = (TH2F*) list->FindObject("eventStatCent");
  TCanvas* can = new TCanvas();
  event_stat->DrawCopy();

  TCanvas* can2 = new TCanvas("can2", "can2");
  TLegend* leg = getLegend();

  for(int i=0; i<=2; i++){
    int bin = event_stat->GetXaxis()->FindBin(i);
    TH1F* event_stat_proj = (TH1F*) event_stat->ProjectionY("_px", bin, bin, "e")->Clone("_clone");

    if(i==0){
      event_stat_proj->SetLineColor(kRed);
      event_stat_proj->GetYaxis()->SetTitle("events");
      leg->AddEntry(event_stat_proj->Clone(), "all events", "l");
    }else if(i==1){
      event_stat_proj->SetLineColor(kBlue);
      leg->AddEntry(event_stat_proj->Clone(), "after centrality and vertex cut", "l");
    }else if(i==2){
      event_stat_proj->SetLineColor(kGreen);
      leg->AddEntry(event_stat_proj->Clone(), "event pool is ready", "l");
    }

    event_stat_proj->SetLineWidth(3.0);

    if(i==0){
      event_stat_proj->SetStats(kFALSE);
      event_stat_proj->DrawCopy();
    }else
      event_stat_proj->DrawCopy("same");

    //leg->AddEntry(event_stat_proj->Clone(), Form("stage %i", i), "l");
  }
      
  leg->Draw("same");
}







//returns maximum of the flow
//flow from 1.0 > |eta| > 1.6 is subtracted from the rest of the histogram
//the flow in the TH2 is subtracted
TH1D* subtractFlow(TH2* etaPhi){

  return subtractFlow(etaPhi, gSubtractFlow_start, gSubtractFlow_end);
}

//returns flow distribution
//flow from 1.0 > |eta| > 1.6 is subtracted from the rest of the histogram
//the flow in the TH2 is subtracted
TH1D* subtractFlow(TH2* etaPhi, Double_t start, Double_t end){
  start += 0.01;
  end -= 0.01;
  int firstbin = etaPhi->GetYaxis()->FindBin(start);
  int lastbin = etaPhi->GetYaxis()->FindBin(end);

  TH1D* etaPhi_proj = etaPhi->ProjectionX("_px", firstbin, lastbin, "e")->Clone("_clone");
  int usedBins = lastbin - firstbin + 1;

  firstbin = etaPhi->GetYaxis()->FindBin(-1.0*end);
  lastbin = etaPhi->GetYaxis()->FindBin(-1.0*start);
  TH1D* signal2 = etaPhi->ProjectionX("_px2", firstbin, lastbin, "e");
  usedBins += lastbin - firstbin + 1;
  etaPhi_proj->Add(signal2, 1.0);

  etaPhi_proj->Scale(1.0/usedBins);

  
  for(int i=0;i<=etaPhi->GetNbinsX();i++){
    double subtract = etaPhi_proj->GetBinContent(i);
    double subtract_err = etaPhi_proj->GetBinError(i);

    for(int j=1;j<=etaPhi->GetNbinsY();j++){
      double content = etaPhi->GetBinContent(i,j)-subtract;
      double error = TMath::Sqrt(TMath::Power(subtract_err, 2) + TMath::Power(etaPhi->GetBinError(i,j), 2));
      etaPhi->SetBinContent(i,j,content);
      etaPhi->SetBinError(i, j, error);
    }
  }
  
  return etaPhi_proj;
}

//returns constant yield at high Delta eta
//flow from 1.0 > |eta| > 1.6 is subtracted from the rest of the histogram via a subtraction of a constant yield
//the flow in the TH2 is subtracted, this makes it impossible to make a baseline subtraction because the baseline is not on 0
Double_t subtractConstBkg(TH2* etaPhi, Double_t start, Double_t end, Double_t ptAssocMin){
  start += 0.01;
  end -= 0.01;
  int firstbin = etaPhi->GetYaxis()->FindBin(start);
  int lastbin = etaPhi->GetYaxis()->FindBin(end);

  TH1D* etaPhi_proj = etaPhi->ProjectionX("_px", firstbin, lastbin, "e")->Clone("_clone");
  int usedBins = lastbin - firstbin + 1;

  firstbin = etaPhi->GetYaxis()->FindBin(-1.0*end);
  lastbin = etaPhi->GetYaxis()->FindBin(-1.0*start);
  TH1D* signal2 = etaPhi->ProjectionX("_px2", firstbin, lastbin, "e");
  usedBins += lastbin - firstbin + 1;
  etaPhi_proj->Add(signal2, 1.0);

  etaPhi_proj->Scale(1.0/usedBins);

  //Double_t phi_bin = 0.5;
  Double_t phi_bin = getPhiBinsForAnalysis(ptAssocMin);


  TF1* fit = new TF1("fit", "[0]", -1.0*phi_bin, phi_bin);
  etaPhi_proj->Fit("fit", "MN", "", -1.0*phi_bin, phi_bin);
  Double_t par0 = fit->GetParameter(0);
  Double_t par0Err = fit->GetParError(0);
  
  Double_t phi_bin_0 = etaPhi->GetXaxis()->FindBin(-0.01);
  Double_t eta_bin_0 = etaPhi->GetYaxis()->FindBin(-0.01);

  for(int i=0;i<=etaPhi->GetNbinsX();i++){
    for(int j=0;j<etaPhi->GetNbinsY();j++){
      double content = etaPhi->GetBinContent(i,j)-par0;
      double error = TMath::Sqrt(TMath::Power(par0Err, 2) + TMath::Power(etaPhi->GetBinError(i,j), 2));
      etaPhi->SetBinContent(i,j,content);
      etaPhi->SetBinError(i, j, error);
    }
  }

  return par0;
}


//returns flow distribution
//flow from 1.0 > |eta| > 1.6 is subtracted from the rest of the histogram
//the flow in the TH2 is subtracted
TH1D* subtractFlow(TH2* etaPhi, TH2* etaPhi_1plus1, Double_t start, Double_t end){
  start += 0.01;
  end -= 0.01;
  int firstbin = etaPhi_1plus1->GetYaxis()->FindBin(start);
  int lastbin = etaPhi_1plus1->GetYaxis()->FindBin(end);

  TH1D* etaPhi_proj = etaPhi_1plus1->ProjectionX("_px", firstbin, lastbin, "e")->Clone("_clone");
  int usedBins = lastbin - firstbin + 1;

  firstbin = etaPhi_1plus1->GetYaxis()->FindBin(-1.0*end);
  lastbin = etaPhi_1plus1->GetYaxis()->FindBin(-1.0*start);
  TH1D* signal2 = etaPhi_1plus1->ProjectionX("_px2", firstbin, lastbin, "e");
  usedBins += lastbin - firstbin + 1;
  etaPhi_proj->Add(signal2, 1.0);

  etaPhi_proj->Scale(1.0/usedBins);

  for(int i=0;i<=etaPhi->GetNbinsX();i++){
    double subtract = etaPhi_proj->GetBinContent(i);
    double subtract_err = etaPhi_proj->GetBinError(i);
    //Printf("bin center of x axis: %f ", etaPhi_1plus1->GetXaxis()->GetBinCenter(i));
    //Printf("subtract %i, %f ", i, subtract);

    for(int j=0;j<etaPhi->GetNbinsY();j++){
      //Printf("content of %i (=%f) is %f ", j, etaPhi->GetYaxis()->GetBinCenter(j), etaPhi->GetBinContent(i,j));
      double content = etaPhi->GetBinContent(i,j)-subtract;
      double error = TMath::Sqrt(TMath::Power(subtract_err, 2) + TMath::Power(etaPhi->GetBinError(i,j), 2));
      etaPhi->SetBinContent(i,j,content);
      etaPhi->SetBinError(i, j, error);
    }
  }

  return etaPhi_proj;
}


//implemented if needed
//in the normal analysis this screws up the whole pT spectrum
void removeWing(TH2D* h2_phi){

  Float_t width = 1.5;
  //TH1* proj = h2_phi->ProjectionY(Form("projx", h2_phi->GetName()), h2_phi->GetXaxis()->FindBin(TMath::Pi() - width), h2_phi->GetXaxis()->FindBin(TMath::Pi()+ width), "e");
  TH1* proj = h2_phi->ProjectionY(Form("projx", h2_phi->GetName()), h2_phi->GetXaxis()->FindBin(3.0), h2_phi->GetXaxis()->FindBin(4.2), "e");
  proj->Fit("pol0", "0");
  //new TCanvas; proj->DrawCopy();
  proj->Divide(proj->GetFunction("pol0"));
  //new TCanvas; proj->DrawCopy();
  //new TCanvas; hist->DrawCopy("SURF1");
  for (Int_t x=1; x<=h2_phi->GetNbinsX(); x++)
    for (Int_t y=1; y<=h2_phi->GetNbinsY(); y++)
      {
	if (proj->GetBinContent(y) <= 0)
	  continue;
	h2_phi->SetBinContent(x, y, h2_phi->GetBinContent(x, y) / proj->GetBinContent(y));
	h2_phi->SetBinError(x, y, h2_phi->GetBinError(x, y) / proj->GetBinContent(y));
      }
  
}

  //use fit method to subtract the baseline if the subtractFlow method doesn't work because of poor statistics
Double_t subtractBaseline(TH1D* h1_phi, Int_t execute = 1){
  
  //do not use information inside the cone to subtract the background
  if(g_phi_bin > gSubtractBaseline_start) 
    gSubtractBaseline_start = g_phi_bin;
  
  TF1* fit = new TF1("fit", "[0]", -1.0*gSubtractBaseline_end, -1.0*gSubtractBaseline_start);
  
  //h1_phi->Fit("fit", "N", "", -TMath::Pi()/2+0.01, -1.0);
  h1_phi->Fit("fit", "MN", "", -1*gSubtractBaseline_end, -1*gSubtractBaseline_start);
  Double_t par0 = fit->GetParameter(0);
  Double_t par0Err = fit->GetParError(0);
  
  TF1* fit2 = new TF1("fit2", "[0]", gSubtractBaseline_start, gSubtractBaseline_end);

  h1_phi->Fit("fit2", "MN", "", gSubtractBaseline_start, gSubtractBaseline_end);
  Double_t par1 = fit2->GetParameter(0);
  Double_t par1Err = fit2->GetParError(0);

  //this is done so that the fit is not drawn
  //h1_phi->Fit("fit", "0", "", 1.0, TMath::Pi()/2-0.01);
 
  Double_t subtract = (par0+par1)/2;//h1_phi->GetMinimum();
  Double_t subtract_err = 0.5*TMath::Sqrt(TMath::Power(par0Err, 2) + TMath::Power(par1Err, 2));
  //Double_t subtract = par1;
  //Double_t subtract_err = par1Err;
  Printf("subtract %f +/- %f ", subtract, subtract_err);
  
  //subtract =  h1_phi->GetBinContent(h1_phi->GetMinimumBin());
  //subtract_err = h1_phi->GetBinError(h1_phi->GetMinimumBin());

  if(execute>0)
    for(int i=0; i<=h1_phi->GetNbinsX(); i++){
      h1_phi->SetBinContent(i, h1_phi->GetBinContent(i) - subtract);
      Double_t bin_error = TMath::Sqrt(TMath::Power(subtract_err, 2) + TMath::Power(h1_phi->GetBinError(i), 2));
      h1_phi->SetBinError(i, bin_error);
    }

  return subtract;
}


Double_t getPhiBinsForAnalysis(Double_t p_t_assoc){

  //if(true) return 0.51;//0.34;

  //2 bins in the first 2 pT bins is showing that 2+1 analysis is already difficult with too many bins, because fluctuations rise up
  if(p_t_assoc<1.0){
    return 0.51;//0.68;
  }else if(p_t_assoc<2.0){
    return 0.51;//0.68;
  }else if(p_t_assoc<3.0){
    return 0.51;
  }else if(p_t_assoc<4.0){
    return 0.34;//0.51;//0.34;
  }else if(p_t_assoc<6.0){
    return 0.34;
  }else{
    return 0.34;
  }
  
}

Double_t getEtaBinsForAnalysis(Double_t p_t_assoc){

  if(true) return 0.39;

  //2 bins in the first 2 pT bins is showing that 2+1 analysis is already difficult with too many bins, because fluctuations rise up
  if(p_t_assoc<1.0){
    return 0.79;
  }else if(p_t_assoc<2.0){
    return 0.59;//0.79;
  }else if(p_t_assoc<3.0){
    return 0.59;//0.59;
  }else if(p_t_assoc<4.0){
    return 0.39;//0.59;
  }else if(p_t_assoc<6.0){
    return 0.39;
  }else{
    return 0.39;
  }
}


void showAllTriggerPt(const char* fileName, Int_t multBinBegin = 5, Int_t multBinEnd = 7, Int_t side = 0){
  
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
      getAnalysis(fileName, pt1Minimum, pt1Maximum, pt2Minimum, gpTMin_assoc, gpTMax_assoc, gVertex, multBinBegin, multBinEnd, side, i, j, 1, 0);
    }
  }
}

void showAllTrigger1Pt(const char* fileName, Int_t multBinBegin = 5, Int_t multBinEnd = 7, Int_t pt2Minimum = 4.0, Double_t ptAssocMin = 0.5, Double_t ptAssocMax = 8.0){

  const Int_t pt1_bins_length = 3;
  Int_t pt1_bins[pt1_bins_length+1] = {6.0, 8.0, 10.0, 12.0};

  for(Int_t i=0; i<pt1_bins_length; i++){
    TH1D* near = getAnalysis(fileName, pt1_bins[i], pt1_bins[i+1], pt2Minimum, ptAssocMin, ptAssocMax, gVertex, multBinBegin, multBinEnd, 0, i, 0, 1, 0);

    TH1D* away = getAnalysis(fileName, pt1_bins[i], pt1_bins[i+1], pt2Minimum, ptAssocMin, ptAssocMax, gVertex, multBinBegin, multBinEnd, 1, i, 1, 1, 0);
    /*
    near->Divide(near, away);
    TCanvas* c1 = new TCanvas(Form("can %i %i", i, 0), Form("can %i %i", i, 0), i*gBasisSize+50, 0, gBasisSize, gBasisSize);
    near->DrawCopy();
    */
  }

}

void showAllAssocPt(const char* fileName, Int_t multBinBegin = 0, Int_t multBinEnd = 16, Int_t pt1Minimum = 4.0, Int_t pt1Maximum = 14.0, Int_t pt2Minimum = 2.0){

  //const Int_t ptAssoc_bins_length = 5;
  //Int_t ptAssoc_bins[ptAssoc_bins_length+1] = {1.0, 2.0, 3.0, 4.0, 6.0, 8.0};
  const Int_t ptAssoc_bins_length = 3;
  Int_t ptAssoc_bins[ptAssoc_bins_length+1] = {1.0, 3.0, 6.0, 8.0};

  for(Int_t i=0; i<ptAssoc_bins_length; i++){
    Int_t draw = 0;

    TH1D* near = getAnalysis(fileName, pt1Minimum, pt1Maximum, pt2Minimum, ptAssoc_bins[i], ptAssoc_bins[i+1], gVertex, multBinBegin, multBinEnd, 0, i, 0, draw, 0)->Clone();

    TH1D* away = getAnalysis(fileName, pt1Minimum, pt1Maximum, pt2Minimum, ptAssoc_bins[i], ptAssoc_bins[i+1], gVertex, multBinBegin, multBinEnd, 1, i, 1, draw, 0)->Clone();

    if(!draw){
      near->Divide(near, away);
      
      TCanvas* c1 = new TCanvas(Form("can %i %i", i, 0), Form("can %i %i", i, 0), i*gBasisSize+50, 0, gBasisSize, gBasisSize);
      near->DrawCopy();
      //TCanvas* c2 = new TCanvas(Form("can2 %i %i", i, 0), Form("can2 %i %i", i, 0), i*gBasisSize+50, gBasisSize+50, gBasisSize, gBasisSize);
      //away->DrawCopy();
    }
  }
}

void showAllMult(const char* fileName, Int_t pt1Minimum = 4.0, Int_t pt1Maximum = 14.0, Int_t pt2Minimum = 2.0, Double_t ptAssocMin = 0.5, Double_t ptAssocMax = 8.0){

  //in multiplicity this are the bins
  //standard binning: 0, 1, 2, 3, 4, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100.1
  const Int_t mult_bins_length = 5;
  //Int_t mult_bins[mult_bins_length+1] = {1, 5, 6, 8, 11, 15};
  Int_t mult_bins[mult_bins_length+1] = {1, 1, 5, 7, 9, 10};

    for(Int_t i=0; i<mult_bins_length; i++){
      TH1D* near = getAnalysis(fileName, pt1Minimum, pt1Maximum, pt2Minimum, ptAssocMin, ptAssocMax, gVertex, mult_bins[i], mult_bins[i+1], 0, i, 0, 1, 0);

      TH1D* away = getAnalysis(fileName, pt1Minimum, pt1Maximum, pt2Minimum, ptAssocMin, ptAssocMax, gVertex, mult_bins[i], mult_bins[i+1], 1, i, 1, 1, 0);
      /*
      near->Divide(near, away);

      TCanvas* c1 = new TCanvas(Form("can %i %i", i, 0), Form("can %i %i", i, 0), i*gBasisSize+50, 0, gBasisSize, gBasisSize);
      near->DrawCopy();
      */
    }

}

void showMultCompare(const char* filename, Int_t pt1Minimum = 4.0, Int_t pt1Maximum = 14.0, Int_t pt2Minimum = 2.0, Double_t ptAssocMin = 0.5, Double_t ptAssocMax = 8.0, Int_t side = 0, Double_t vtx = 8.9, Int_t subtractMixedComb = 1, Int_t subtractFlow = 1, Int_t showPlot = 1){

  Int_t show2Dplots = 2;

  TLegend *leg  = new TLegend(0.65,0.7,0.85,0.9);
  leg->SetFillColor(10);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextFont(gStyle->GetTextFont());
  leg->SetTextSize(gStyle->GetTextSize()*1.0);

  TH1D* central = (TH1D*)(getAnalysis(filename, pt1Minimum, pt1Maximum, pt2Minimum, ptAssocMin, ptAssocMax, vtx, 1, 6, side, showPlot, 0, show2Dplots, subtractMixedComb, subtractFlow))->Clone("_central");

  central->SetLineColor(kBlue);
  leg->AddEntry(central->Clone(),"central","l");

  TH1D* semiCentral = getAnalysis(filename, pt1Minimum, pt1Maximum, pt2Minimum, ptAssocMin, ptAssocMax, vtx, 10, 11, side, showPlot+1, 0, show2Dplots, subtractMixedComb, subtractFlow);

  TCanvas* multCompare = new TCanvas(Form("multCompare  %i", showPlot), Form("multCompare %i", showPlot), showPlot*gBasisSize+100, showPlot*gBasisSize+50, 2*gBasisSize, 2*gBasisSize);
  central->DrawCopy();

  semiCentral->SetLineColor(kRed);
  semiCentral->DrawCopy("same");  
  leg->AddEntry(semiCentral->Clone(),"semiCentral","l");


  leg->Draw("same");

  //show same and mixed event of one vertex and multiplicity bin
  //getRawAnalysis(filename, pt1Minimum, pt1Maximum, pt2Minimum, ptAssocMin, ptAssocMax, 3, 3);
}

void showAnalysis_all(const char* fileName, double pt1Min = 4.0, double pt1Max = 14.0, double pt2Min = 2.0, Double_t ptAssocMin = 0.5, Double_t ptAssocMax = 8.0, Int_t subtractMixedComb = 0, Int_t subtractFlow = 0){

  showAnalysis(fileName, pt1Min, pt1Max, pt2Min, ptAssocMin, ptAssocMax, 1, 6, subtractMixedComb, subtractFlow, 0);
  showAnalysis(fileName, pt1Min, pt1Max, pt2Min, ptAssocMin, ptAssocMax, 10, 11, subtractMixedComb, subtractFlow, 1);

}


void showAnalysis(const char* fileName, double pt1Min = 4.0, double pt1Max = 14.0, double pt2Min = 2.0, Double_t ptAssocMin = 0.5, Double_t ptAssocMax = 8.0, Int_t multBinBegin = 1, Int_t multBinEnd = 5, Int_t subtractMixedComb = 0, Int_t subtractFlow = 0, Int_t subtractBaseline = 0){

  Double_t pt2Max = pt2Min+ptMaxAdd;

  Int_t showPlots = 0;
  Double_t vertex = 7.9;
  /*
  TH1* near = (TH1*)getAnalysis(fileName, pt1Min, pt1Max, pt2Min, ptAssocMin, ptAssocMax, vertex, multBinBegin, multBinEnd, 0, 1, 1, showPlots, subtractMixedComb, subtractFlow, NULL, NULL, 0)->Clone();
  TH1* near_pure = getAnalysis(fileName, pt1Min, pt1Max, pt2Min, ptAssocMin, ptAssocMax, vertex, multBinBegin, multBinEnd, 6, 2, 1, showPlots, subtractMixedComb, subtractFlow, NULL, NULL, 0)->Clone();
  */
  //near - away side
  for(i=0; i<=1; i++){
    TH2D* h2_near = getAnalysis_2D(fileName, pt1Min, pt1Max, pt2Min, ptAssocMin, ptAssocMax, vertex, multBinBegin, multBinEnd, i, 1, 1, NULL, NULL)->Clone("h2_etaphi_clone");
    TH2D* h2_near_pure;
    if(i==0)
      h2_near_pure = (TH2D*)getAnalysis_2D(fileName, pt1Min, pt1Max, pt2Min, ptAssocMin, ptAssocMax, vertex, multBinBegin, multBinEnd, 6, 1, 1, NULL, NULL)->Clone("h2_etaphi_pure_clone");
    else
      h2_near_pure = (TH2D*)getAnalysis_2D(fileName, pt1Min, pt2Min, pt2Max, ptAssocMin, ptAssocMax, vertex, multBinBegin, multBinEnd, 6, 1, 1, NULL, NULL)->Clone("h2_etaphi_pure_clone");

    if(subtractFlow){
      subtractFlow(h2_near, 1.0, 1.4);
      subtractFlow(h2_near_pure, 1.0, 1.4);
    }
    
    Double_t eta_bin = getEtaBinsForAnalysis(ptAssocMin);
    int firstbin = h2_near->GetYaxis()->FindBin(-1*eta_bin);
    int lastbin = h2_near->GetYaxis()->FindBin(eta_bin);
    
    Double_t phi_bin = getPhiBinsForAnalysis(ptAssocMin);
    int firstbin_phi = h2_near->GetYaxis()->FindBin(-1*eta_bin);
    int lastbin_phi = h2_near->GetYaxis()->FindBin(eta_bin);
    
    TH1D* near;
    TH1D* near_pure;
    TCanvas* c1;

    
    //phi - eta projection
    for(int j=0; j<=1; j++){
      if(j==0){
	near = (TH1D*)projectToTH1D(h2_near, "fully corrected", firstbin, lastbin, 1.0, subtractBaseline)->Clone("_near");
	near_pure = (TH1D*)projectToTH1D(h2_near_pure, "fully corrected", firstbin, lastbin, 1.0, subtractBaseline)->Clone("_near_pure");
	c1 = new TCanvas(Form("side %i, phi projection", i), Form("side %i, phi projection", i), (1+i)*gBasisSize+50, (1+j)*gBasisSize+50, gBasisSize, gBasisSize);
      }else if (j==1){
	near = (TH1D*)projectToTH1D_eta(h2_near, "fully corrected", firstbin_phi, lastbin_phi, 1.0)->Clone("_near");
	near_pure = (TH1D*)projectToTH1D_eta(h2_near_pure, "fully corrected", firstbin_phi, lastbin_phi, 1.0)->Clone("_near_pure");
	c1 = new TCanvas(Form("side %i, eta projection", i), Form("side %i, eta projection", i), (1+i)*gBasisSize+50, (1+j)*gBasisSize+50, gBasisSize, gBasisSize);
      }
    
      near->DrawCopy();  

      near_pure->SetLineColor(kRed);
      near_pure->DrawCopy("same");

      TLegend* leg = getLegend();
      
      if(i==0){
	leg->AddEntry(near->Clone(), "near side", "l");
	leg->AddEntry(near_pure, "near pure BG", "l");
      }else if(i==1){
	leg->AddEntry(near->Clone(), "away side", "l");
	leg->AddEntry(near_pure, "away pure BG", "l");
      }
      
      leg->Draw("same");
    }
  }




  /*
    Double_t phi_bin = getPhiBinsForAnalysis(ptAssocMin);
    Int_t bin_start = near->FindBin(-1*phi_bin);
    Int_t bin_end = near->FindBin(phi_bin);

  Double_t error[4];

  Double_t base1 = subtractBaseline((TH1D*)near, 0);
  Double_t base2 = subtractBaseline((TH1D*)near_pure, 0);
  Double_t base3 = subtractBaseline((TH1D*)away, 0);
  Double_t base4 = subtractBaseline((TH1D*)away_pure, 0);

    Double_t error;
    Double_t content = near->IntegralAndError(bin_start, bin_end, error[0]);
    Printf("near %f +/- %f, subtraction %f ", content, error[0], base1);

    Double_t error2;
    Double_t content2 = near_pure->IntegralAndError(bin_start, bin_end, error[1]);
    Printf("near pure BG %f +/- %f, subtraction %f ", content2, error[1], base2);

    Double_t error3;
    Double_t content3 = away->IntegralAndError(bin_start, bin_end, error[2]);
    Printf("away %f +/- %f, subtraction %f ", content3, error[2], base3);

    Double_t error4;
    Double_t content4 = away_pure->IntegralAndError(bin_start, bin_end, error[3]);
    Printf("away pure BG %f +/- %f, subtraction %f ", content4, error[3], base4);

  */
}



void showNearAwaySide(const char* fileName, double pt1Min = 4.0, double pt1Max = 14.0, double pt2Min = 2.0, Double_t ptAssocMin = 0.5, Double_t ptAssocMax = 8.0, Int_t multBinBegin = 1, Int_t multBinEnd = 5, Int_t subtractMixedComb = 0, Int_t subtractFlow = 0){

  Double_t pt2Max = pt2Min+ptMaxAdd;

  Int_t showPlots = 0;
  Double_t vertex = 7.9;

  TH2D* h2_near = getAnalysis_2D(fileName, pt1Min, pt1Max, pt2Min, ptAssocMin, ptAssocMax, vertex, multBinBegin, multBinEnd, 0, subtractMixedComb, subtractFlow, NULL, NULL)->Clone("h2_etaphi_near");

  TH2D* h2_away = getAnalysis_2D(fileName, pt1Min, pt1Max, pt2Min, ptAssocMin, ptAssocMax, vertex, multBinBegin, multBinEnd, 1, subtractMixedComb, subtractFlow, NULL, NULL)->Clone("h2_etaphi_away");

  if(subtractFlow){
    subtractFlow(h2_near, 1.0, 1.4);
    subtractFlow(h2_away, 1.0, 1.4);
  }
    
  h2_near->GetXaxis()->SetRangeUser(-TMath::Pi()/2+0.01, TMath::Pi()/2-0.01);
  h2_away->GetXaxis()->SetRangeUser(-TMath::Pi()/2+0.01, TMath::Pi()/2-0.01);

  Double_t eta_bin = getEtaBinsForAnalysis(ptAssocMin);
  int firstbin = h2_near->GetYaxis()->FindBin(-1*eta_bin);
  int lastbin = h2_near->GetYaxis()->FindBin(eta_bin);

  TH1D* near = (TH1D*)projectToTH1D(h2_near, "fully corrected", firstbin, lastbin, 1.0, 0)->Clone("_near");
  TH1D* away = (TH1D*)projectToTH1D(h2_away, "fully corrected", firstbin, lastbin, 1.0, 0)->Clone("_away");

  TCanvas *c1 = new TCanvas("phi projection", "phi projection", gBasisSize+50, gBasisSize+50, gBasisSize, gBasisSize);

  near->DrawCopy();
  away->SetLineColor(kRed);
  away->DrawCopy("same");

  TLegend* leg = getLegend();
  leg->AddEntry(near->Clone(), "T1 associated", "l");
  leg->AddEntry(away, "T2 associated", "l");
      
  leg->Draw("same");
}





/*

void peakDifference_mult(const char* fileName, double pt1Min = 4.0, double pt1Max = 14.0, double pt2Min = 2.0, Int_t subtractMixedComb = 1, Int_t subtractFlow = 1){

  Int_t side = 0;//no effect because here both sides for one multiplicity are shown
  Int_t mode = 0;
  peakDifference(fileName, pt1Min, pt1Max, pt2Min, 1, 5, side, 0, mode, subtractMixedComb, subtractFlow);
  peakDifference(fileName, pt1Min, pt1Max, pt2Min, 9, 10, side, 1, mode, subtractMixedComb, subtractFlow);

}

void peakDifference_side(const char* fileName, double pt1Min = 4.0, double pt1Max = 14.0, double pt2Min = 2.0, Int_t subtractMixedComb = 1, Int_t subtractFlow = 1){

  Int_t mult = 0;//no effect because here both multiplicities for one side are shown
  Int_t mode = 1;
  peakDifference(fileName, pt1Min, pt1Max, pt2Min, mult, mult, 0, 0, mode, subtractMixedComb, subtractFlow);
  peakDifference(fileName, pt1Min, pt1Max, pt2Min, mult, mult, 1, 1, mode, subtractMixedComb, subtractFlow);

}

//mode 0 near and away side (two different multiplicities in the two rows)
//mode 1 central and semi central (near side in top and away side in bottom row)
void peakDifference(const char* fileName, double pt1Min = 4.0, double pt1Max = 14.0, double pt2Min = 2.0, Int_t multBinBegin = 1, Int_t multBinEnd = 5, Int_t side = 0, Int_t yPos = 0, Int_t mode = 0,  Int_t subtractMixedComb = 1, Int_t subtractFlow = 1, Int_t draw = 0){

  Double_t z_vtx = 8.9;

  Double_t near_content[pt_assoc_bins_number-1];
  Double_t away_content[pt_assoc_bins_number-1];
  Double_t onePlusOne_content[pt_assoc_bins_number-1];
  Double_t near_error[pt_assoc_bins_number-1];
  Double_t away_error[pt_assoc_bins_number-1];
  Double_t onePlusOne_error[pt_assoc_bins_number-1];

  for(int i=0; i<pt_assoc_bins_number-1; i++){
    TH1* near;
    TH1* away;
    TH1* onePlusOne;
    TCanvas* can;

    if(mode==0){
      near = (TH1*)getAnalysis(fileName, pt1Min, pt1Max, pt2Min, pt_assoc_bins[i], pt_assoc_bins[i+1], z_vtx, multBinBegin, multBinEnd, 0, 1, 1, 0, subtractMixedComb, subtractFlow)->Clone();
      away = (TH1*)getAnalysis(fileName, pt1Min, pt1Max, pt2Min, pt_assoc_bins[i], pt_assoc_bins[i+1], z_vtx, multBinBegin, multBinEnd, 1, 2, 1, 0, subtractMixedComb, subtractFlow)->Clone();
      onePlusOne = (TH1*)getAnalysis(fileName, pt1Min, pt1Max, pt2Min, pt_assoc_bins[i], pt_assoc_bins[i+1], 6.9, multBinBegin, multBinEnd, 6, 2, 1, 0, subtractMixedComb, subtractFlow);
	//onePlusOne = (TH1*)getAnalysis(fileName, pt2Min, pt1Min, pt2Min, pt_assoc_bins[i], pt_assoc_bins[i+1], 6.9, multBinBegin, multBinEnd, 6, 2, 1, 0, subtractMixedComb, subtractFlow);

      can = new TCanvas(Form(" %i, mult %i-%i ", i, multBinBegin, multBinEnd), Form(" %i, mult %i-%i ", i, multBinBegin, multBinEnd), i*gBasisSize+100, yPos*(gBasisSize+50), gBasisSize, gBasisSize);
      
    }else if(mode==1){
      near = (TH1*)getAnalysis(fileName, pt1Min, pt1Max, pt2Min, pt_assoc_bins[i], pt_assoc_bins[i+1], z_vtx, 1, 5, side, 1, 1, 0, subtractMixedComb, subtractFlow)->Clone();
      away = (TH1*)getAnalysis(fileName, pt1Min, pt1Max, pt2Min, pt_assoc_bins[i], pt_assoc_bins[i+1], z_vtx, 9, 10, side, 2, 1, 0, subtractMixedComb, subtractFlow);

      can = new TCanvas(Form(" %i, side %i ", i, side), Form(" %i, side %i ", i, side), i*gBasisSize+50, yPos*(gBasisSize+50), gBasisSize, gBasisSize);
    }else{
      Printf("mode does not exist!");
      return;
    }
    away->SetLineColor(kRed);
    if(mode==0)
      onePlusOne->SetLineColor(kGreen);
    
    Double_t phi_bin = getPhiBinsForAnalysis(pt_assoc_bins[i]);

    Int_t bin_start = near->FindBin(-1*phi_bin);
    Int_t bin_end = near->FindBin(phi_bin);

    near_content[i] = near->IntegralAndError(bin_start, bin_end, near_error[i]);
    away_content[i] = away->IntegralAndError(bin_start, bin_end, away_error[i]);
    if(mode==0)
      onePlusOne_content[i] = onePlusOne->IntegralAndError(bin_start, bin_end, onePlusOne_error[i]);

    ((TH1*)(near->Clone()))->DrawCopy();
    ((TH1*)(away->Clone()))->DrawCopy("same");
    if(mode==0)
      ((TH1*)(onePlusOne->Clone()))->DrawCopy("same");

    if(i==4){
      TLegend* leg = getLegend();
      if(mode==0){
	leg->AddEntry(near,"near side","l");
	leg->AddEntry(away,"away side","l");
	leg->AddEntry(onePlusOne, "1plus1", "l");
      }else if(mode==1){
	leg->AddEntry(near,"central","l");
	leg->AddEntry(away,"semi central","l");
      }
      leg->Draw("same");
    }
  }

  double sum_near = 0;
  double sum_away = 0;
  double sum_1plus1 = 0;
  double err_near = 0;
  double err_away = 0;
  double err_1plus1 = 0;

  for(int i=0; i<pt_assoc_bins_number-1; i++){
    Double_t diff = near_content[i] - away_content[i];
    Double_t err = TMath::Sqrt(TMath::Power(near_error[i], 2) + TMath::Power(away_error[i], 2));

    if(mode==0){
      Printf("near side (%.1f - %.1f): %f +/- %f, away side %f +/- %f, difference: %f +/- %f ", pt_assoc_bins[i], pt_assoc_bins[i+1], near_content[i], near_error[i], away_content[i], away_error[i], diff, err);
      sum_near += near_content[i]*(pt_assoc_bins[i] + pt_assoc_bins[i+1])/2;
      sum_away += away_content[i]*(pt_assoc_bins[i] + pt_assoc_bins[i+1])/2;
      sum_1plus1 += onePlusOne_content[i]*(pt_assoc_bins[i] + pt_assoc_bins[i+1])/2;
      err_near += near_error[i]*(pt_assoc_bins[i] + pt_assoc_bins[i+1])/2;
      err_away += away_error[i]*(pt_assoc_bins[i] + pt_assoc_bins[i+1])/2;
      err_1plus1 += onePlusOne_error[i]*(pt_assoc_bins[i] + pt_assoc_bins[i+1])/2;
    }else
      Printf("central (%.1f - %.1f): %f +/- %f, semi central %f +/- %f, difference: %f +/- %f ", pt_assoc_bins[i], pt_assoc_bins[i+1], near_content[i], near_error[i], away_content[i], away_error[i], diff, err);
  }

  if(mode==0){
    Printf("Near side pT: %f +/- %f", sum_near, err_near);
    Printf("Away side pT: %f +/- %f", sum_away, err_away);
    Printf("1plus1 pT: %f +/- %f", sum_1plus1, err_1plus1);
  }

  //scale content with the bin width
  for(int i=0; i<pt_assoc_bins_number-1; i++){
    near_content[i] /= 2*pt_assoc_bin_error[i];
    near_error[i] /= 2*pt_assoc_bin_error[i];
    away_content[i] /= 2*pt_assoc_bin_error[i];
    away_error[i] /= 2*pt_assoc_bin_error[i];
    if(mode==0){
      onePlusOne_content[i] /= 2*pt_assoc_bin_error[i];
      onePlusOne_error[i] /= 2*pt_assoc_bin_error[i];
    }
  }
  / *
  //divide content by 1+1 content
  for(int i=0; i<pt_assoc_bins_number-1; i++){
    near_content[i] /= onePlusOne_content[i];
    near_error[i] /= onePlusOne_content[i];
    away_content[i] /= onePlusOne_content[i];
    away_error[i] /= onePlusOne_content[i];
  }
  * /
  TGraphErrors* graph_near = new TGraphErrors(pt_assoc_bins_number, pt_assoc_bin_center, near_content, pt_assoc_bin_error, near_error);
  TGraphErrors* graph_away = new TGraphErrors(pt_assoc_bins_number, pt_assoc_bin_center, away_content, pt_assoc_bin_error, away_error);
  TGraphErrors* graph_onePlusone = new TGraphErrors(pt_assoc_bins_number, pt_assoc_bin_center, onePlusOne_content, pt_assoc_bin_error, onePlusOne_error);

  TCanvas* can_graph = new TCanvas(Form("result %i", yPos), Form("result %i", yPos), pt_assoc_bins_number*gBasisSize+50, yPos*(gBasisSize+50), gBasisSize, gBasisSize);

  graph_near->GetXaxis()->SetTitle("p_{T,assoc} (GeV/c)");
  graph_near->GetYaxis()->SetTitle("1/N dN/dpT");
  graph_near->SetMarkerSize(2);
  graph_near->SetLineWidth(2);
  graph_near->SetMarkerColor(kBlue);
  graph_near->SetLineColor(kBlue);
  graph_near->SetMarkerStyle(20);
  graph_away->GetXaxis()->SetTitle("p_{T,assoc} (GeV/c)");
  graph_away->GetYaxis()->SetTitle("1/N dN/dpT");
  graph_away->SetMarkerSize(2);
  graph_away->SetLineWidth(2);
  graph_away->SetMarkerColor(kRed);
  graph_away->SetLineColor(kRed);
  graph_away->SetMarkerStyle(22);
  graph_near->Draw("AP");
  graph_away->Draw("P");

  if(false && mode==0){
    graph_onePlusone->GetXaxis()->SetTitle("p_{T,assoc} (GeV/c)");
    graph_onePlusone->GetYaxis()->SetTitle("1/N dN/dpT");
    graph_onePlusone->SetMarkerSize(2);
    graph_onePlusone->SetLineWidth(2);
    graph_onePlusone->SetMarkerColor(kGreen);
    graph_onePlusone->SetLineColor(kGreen);
    graph_onePlusone->SetMarkerStyle(22);
    graph_onePlusone->Draw("P");
  }
  gPad->SetLogy();

  TLegend* leg2 = getLegend();
  if(mode==0){
    leg2->AddEntry(graph_near,"near side","l");
    leg2->AddEntry(graph_away,"away side","l");
    leg2->AddEntry(graph_onePlusone,"1plus1","l");
  }else if(mode==1){
    leg2->AddEntry(graph_near,"central","l");
    leg2->AddEntry(graph_away,"semi central","l");
  }
  leg2->Draw("same");
  
}
*/






void create_peakDifference_pictures(const char* fileName){

  Int_t subtractMixedComb = 1;
  Int_t subtractFlow = 1;//if this is 2, turn of subtractBaseline

  gROOT->SetBatch(kTRUE);
  //useful for big axis labels on the y axis
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadRightMargin(0.07);
  
  for(int mode=1; mode<=5; mode++){
    if(mode!=3) continue;
    //peakDifference_draw(fileName, 4, 6, 4, 1, subtractMixedComb, subtractFlow, mode);
    //peakDifference_draw(fileName, 8, 16, 4, 1, subtractMixedComb, subtractFlow, mode);
  
    peakDifference_draw(fileName, 6, 8, 4, 1, subtractMixedComb, subtractFlow, mode);
    //peakDifference_draw(fileName, 6, 8, 6, 1, subtractMixedComb, subtractFlow, mode);
    //peakDifference_draw(fileName, 8, 12, 4, 1, subtractMixedComb, subtractFlow, mode);
    //peakDifference_draw(fileName, 8, 12, 6, 1, subtractMixedComb, subtractFlow, mode);

    //peakDifference_draw(fileName, 10, 15, 4, 1, subtractMixedComb, subtractFlow, mode);
    //peakDifference_draw(fileName, 5, 10, 4, 1, subtractMixedComb, subtractFlow, mode);

    //peakDifference_draw(fileName, 8, 10, 4, 1, subtractMixedComb, subtractFlow, mode);

    //peakDifference_draw(fileName, 12, 16, 4, 1, subtractMixedComb, subtractFlow, mode);
    //peakDifference_draw(fileName, 12, 16, 6, 1, subtractMixedComb, subtractFlow, mode);
    //peakDifference_draw(fileName, 12, 16, 8, 1, subtractMixedComb, subtractFlow, mode);
    //peakDifference_draw(fileName, 8, 12, 8, 1, subtractMixedComb, subtractFlow, mode);
    //peakDifference_draw(fileName, 8, 12, 10, 1, subtractMixedComb, subtractFlow, mode);

    //peakDifference_draw(fileName, 8, 12, 6, 1, subtractMixedComb, subtractFlow, mode);

    //peakDifference_draw(fileName, 16, 20, 4, 1, subtractMixedComb, subtractFlow, mode);

    //peakDifference_draw(fileName, 12, 20, 8, 1, subtractMixedComb, subtractFlow, mode);

    //peakDifference_draw(fileName, 8, 10, 4, 1, subtractMixedComb, subtractFlow, mode);
    //peakDifference_draw(fileName, 10, 12, 4, 1, subtractMixedComb, subtractFlow, mode);
    //peakDifference_draw(fileName, 12, 14, 4, 1, subtractMixedComb, subtractFlow, mode);

    /*
      Printf("peakDifference_draw(fileName, 6, 8, 4, 1, 1, 1");
      peakDifference_draw(fileName, 6, 8, 4, 1, subtractMixedComb, subtractFlow); 
      Printf("peakDifference_draw(fileName, 6, 8, 6, 1, 1, 1");
      peakDifference_draw(fileName, 6, 8, 6, 1, subtractMixedComb, subtractFlow);
      Printf("peakDifference_draw(fileName, 8, 10, 4, 1, 1, 1");
      peakDifference_draw(fileName, 8, 10, 4, 1, subtractMixedComb, subtractFlow);
      Printf("peakDifference_draw(fileName, 8, 10, 6, 1, 1, 1");
      peakDifference_draw(fileName, 8, 10, 6, 1, subtractMixedComb, subtractFlow);
      Printf("peakDifference_draw(fileName, 10, 12, 4, 1, 1, 1");
      peakDifference_draw(fileName, 10, 12, 4, 1, subtractMixedComb, subtractFlow);
      Printf("peakDifference_draw(fileName, 10, 12, 6, 1, 1, 1");
      peakDifference_draw(fileName, 10, 12, 6, 1, subtractMixedComb, subtractFlow);
    */


    gDirectory->GetList()->Delete();
  }

  TCanvas* can = new TCanvas("can_filename", "can_filename");
  TH1F* mixedDist = getMixedDist(fileName);
  mixedDist->SetTitle(fileName);
  mixedDist->Draw("colz");

  can->SaveAs("pt_spectra/filename.eps");
  
  TCanvas* info_can = new TCanvas("run info", "run info");
  TLatex * name = new TLatex(0.35,0.945,fileName);
  name->SetTextFont(62);
  name->SetTextSize(0.05);
  name->Draw();
  

  TLatex * phi = new TLatex(0.65,0.645,Form("|\\Delta \\varphi| < %.2f", g_phi_bin));
  //TLatex * phi = new TLatex(0.65,0.645,"use phi areas dependent on the associated pt");
  phi->SetTextFont(62);
  phi->SetTextSize(0.05);
  phi->Draw();

  TLatex * eta = new TLatex(0.65,0.545,Form("|\\Delta \\eta| < %.2f", g_eta_bin));
  //TLatex * eta = new TLatex(0.65,0.545,"use eta areas dependent on the associated pt");
  eta->SetTextFont(62);
  eta->SetTextSize(0.05);
  eta->Draw("same");

  info_can->SaveAs("pt_spectra/run_info.eps");
}











void compareSpectra(const char* fileName, const char* fileName2, double pt1Min = 4.0, double pt1Max = 14.0, double pt2Min = 2.0,  Int_t subtractMixedComb = 1, Int_t subtractFlow = 1){

  Double_t pt2Max = pt2Min+ptMaxAdd;

  Double_t saveGraph_min = 0.;
  Double_t saveGraph_max = 3.;

  Int_t trigger_central_near = 0;
  Int_t trigger_semi_near = 0;
  Int_t trigger_pp = 0;

  Int_t trigger_bg_central_near = 0;
  Int_t trigger_bg_semi_near = 0;
  Int_t trigger_bg_pp = 0;

  Int_t mult_cent_min = 1;
  Int_t mult_cent_max = 6;//9;

  Int_t mult_semi_min = 10;
  Int_t mult_semi_max = 11;

  Int_t mult_pp_min = 1;
  Int_t mult_pp_max = 1;

  Int_t yPos = 0;
  /*
  TGraphErrors* graph_central_near = (TGraphErrors*)(peakDifference_graph(fileName, pt1Min, pt1Max, pt2Min, mult_cent_min, mult_cent_max, 6, yPos, subtractMixedComb, subtractFlow, &trigger_central_near, &trigger_bg_central_near)->Clone("test1"));
  TGraphErrors* graph_semi_near = (TGraphErrors*)(peakDifference_graph(fileName, pt1Min, pt1Max, pt2Min, mult_semi_min, mult_semi_max, 6, yPos, subtractMixedComb, subtractFlow, &trigger_semi_near, &trigger_bg_semi_near)->Clone("test2"));
  TGraphErrors* graph_pp_near = (TGraphErrors*)(peakDifference_graph(fileName2, pt1Min, pt1Max, pt2Min, mult_pp_min, mult_pp_max, 6, yPos, subtractMixedComb, subtractFlow, &trigger_pp, &trigger_bg_pp)->Clone());

  save_graph_ratio(graph_central_near, graph_semi_near, graph_central_near, graph_pp_near,  pt_assoc_bins_number, Form("ICP and IAA near side for %.1f < p_{T,1} < %.1f and %.1f < p_{T,2} < %.1f", pt1Min, pt1Max, pt2Min, pt2Max), Form("pt_spectra/ICP_IAA_1plus1_near_%.0f_%.0f_%.0f.eps", pt1Min, pt1Max, pt2Min), "I\_{CP}", "I\_{AA}", saveGraph_min, saveGraph_max, pt2Min);
//save_graph_ratio(graph_central_near, graph_pp_near, graph_semi_near, graph_pp_near,  pt_assoc_bins_number, Form("ICP and IAA near side for %.1f < p_{T,1} < %.1f and %.1f < p_{T,2} < %.1f", pt1Min, pt1Max, pt2Min, pt2Max), Form("pt_spectra/ICP_IAA_1plus1_near_%.0f_%.0f_%.0f.eps", pt1Min, pt1Max, pt2Min), "I\_{CP}", "I\_{AA}", saveGraph_min, saveGraph_max, pt2Min);
  

  TGraphErrors* graph_central_away = (TGraphErrors*)(peakDifference_graph(fileName, pt2Min, pt2Max, pt2Min, mult_cent_min, mult_cent_max, 6, yPos, subtractMixedComb, subtractFlow, &trigger_central_near, &trigger_bg_central_near)->Clone("test1"));
  TGraphErrors* graph_semi_away = (TGraphErrors*)(peakDifference_graph(fileName, pt2Min, pt2Max, pt2Min, mult_semi_min, mult_semi_max, 6, yPos, subtractMixedComb, subtractFlow, &trigger_semi_near, &trigger_bg_semi_near)->Clone("test2"));
  TGraphErrors* graph_pp_away = (TGraphErrors*)(peakDifference_graph(fileName2, pt2Min, pt2Max, pt2Min, mult_pp_min, mult_pp_max, 6, yPos, subtractMixedComb, subtractFlow, &trigger_pp, &trigger_bg_pp)->Clone());

  save_graph_ratio(graph_central_away, graph_semi_away, graph_central_away, graph_pp_away,  pt_assoc_bins_number, Form("ICP and IAA near side for %.1f < p_{T,1} < %.1f and %.1f < p_{T,2} < %.1f", pt1Min, pt1Max, pt2Min, pt2Max), Form("pt_spectra/ICP_IAA_1plus1_away_%.0f_%.0f_%.0f.eps", pt1Min, pt1Max, pt2Min), "I\_{CP}", "I\_{AA}", saveGraph_min, saveGraph_max, pt2Min);
 //save_graph_ratio(graph_central_away, graph_pp_away, graph_semi_away, graph_pp_away,  pt_assoc_bins_number, Form("ICP and IAA near side for %.1f < p_{T,1} < %.1f and %.1f < p_{T,2} < %.1f", pt1Min, pt1Max, pt2Min, pt2Max), Form("pt_spectra/ICP_IAA_1plus1_away_%.0f_%.0f_%.0f.eps", pt1Min, pt1Max, pt2Min), "I\_{CP}", "I\_{AA}", saveGraph_min, saveGraph_max, pt2Min);
  
  */
  
  TGraphErrors* graph_central_near = (TGraphErrors*)(peakDifference_graph(fileName, pt1Min, pt1Max, pt2Min, mult_cent_min, mult_cent_max, 0, yPos, subtractMixedComb, subtractFlow, &trigger_central_near, &trigger_bg_central_near)->Clone());
  TGraphErrors* graph_semi_near = (TGraphErrors*)(peakDifference_graph(fileName, pt1Min, pt1Max, pt2Min, mult_semi_min, mult_semi_max, 0, yPos, subtractMixedComb, subtractFlow, &trigger_semi_near, &trigger_bg_semi_near)->Clone());
  TGraphErrors* graph_pp_near = (TGraphErrors*)(peakDifference_graph(fileName2, pt1Min, pt1Max, pt2Min, mult_pp_min, mult_pp_max, 0, yPos, subtractMixedComb, subtractFlow, &trigger_pp, &trigger_bg_pp)->Clone());

  save_graph_ratio(graph_central_near, graph_semi_near, graph_central_near, graph_pp_near,  pt_assoc_bins_number, Form("ICP and IAA near side for %.1f < p_{T,1} < %.1f and %.1f < p_{T,2} < %.1f", pt1Min, pt1Max, pt2Min, pt2Max), Form("pt_spectra/ICP_IAA_near_%.0f_%.0f_%.0f.eps", pt1Min, pt1Max, pt2Min), "I\_{CP}", "I\_{AA}", saveGraph_min, saveGraph_max, pt2Min);
  //save_graph_ratio(graph_central_near, graph_pp_near, graph_semi_near, graph_pp_near,  pt_assoc_bins_number, Form("ICP and IAA near side for %.1f < p_{T,1} < %.1f and %.1f < p_{T,2} < %.1f", pt1Min, pt1Max, pt2Min, pt2Max), Form("pt_spectra/ICP_IAA_near_%.0f_%.0f_%.0f.eps", pt1Min, pt1Max, pt2Min), "I\_{CP}", "I\_{AA}", saveGraph_min, saveGraph_max, pt2Min);
  
 
  TGraphErrors* graph_central_away = (TGraphErrors*)(peakDifference_graph(fileName, pt1Min, pt1Max, pt2Min, mult_cent_min, mult_cent_max, 1, yPos, subtractMixedComb, subtractFlow, &trigger_central_near, &trigger_bg_central_near)->Clone());
  TGraphErrors* graph_semi_away = (TGraphErrors*)(peakDifference_graph(fileName, pt1Min, pt1Max, pt2Min, mult_semi_min, mult_semi_max, 1, yPos, subtractMixedComb, subtractFlow, &trigger_semi_near, &trigger_bg_semi_near)->Clone());
  TGraphErrors* graph_pp_away = (TGraphErrors*)(peakDifference_graph(fileName2, pt1Min, pt1Max, pt2Min, mult_pp_min, mult_pp_max, 1, yPos, subtractMixedComb, subtractFlow, &trigger_pp, &trigger_bg_pp)->Clone());

  save_graph_ratio(graph_central_away, graph_semi_away, graph_central_away, graph_pp_away,  pt_assoc_bins_number, Form("ICP and IAA away side for %.1f < p_{T,1} < %.1f and %.1f < p_{T,2} < %.1f", pt1Min, pt1Max, pt2Min, pt2Max), Form("pt_spectra/ICP_IAA_away_%.0f_%.0f_%.0f.eps", pt1Min, pt1Max, pt2Min), "I\_{CP}", "I\_{AA}", saveGraph_min, saveGraph_max, pt2Min);
  //save_graph_ratio(graph_central_away, graph_pp_away, graph_semi_away, graph_pp_away,  pt_assoc_bins_number, Form("ICP and IAA away side for %.1f < p_{T,1} < %.1f and %.1f < p_{T,2} < %.1f", pt1Min, pt1Max, pt2Min, pt2Max), Form("pt_spectra/ICP_IAA_away_%.0f_%.0f_%.0f.eps", pt1Min, pt1Max, pt2Min), "I\_{CP}", "I\_{AA}", saveGraph_min, saveGraph_max, pt2Min);
  
}

void save_IAA_all(){

  gROOT->SetBatch(kTRUE);

  Int_t subtractMixedComb = 1;
  Int_t subtractFlow = 1;
  
  Double_t pt2Min = 4;

  for(int i=0; i<pt_assoc_bins_number-1; i++){
    //if(i!=3) continue;
    show_IAA_cent("AR_1509.root", "AR_pp_687.root", "I_AA_Npart_near", 8, 12, pt2Min, pt_assoc_bins[i], pt_assoc_bins[i+1], 0, 7.9, subtractMixedComb, subtractFlow);
    gDirectory->GetList()->Delete();

    show_IAA_cent("AR_1509.root", "AR_pp_687.root", "I_AA_Npart_away", 8, 12, pt2Min, pt_assoc_bins[i], pt_assoc_bins[i+1], 1, 7.9, subtractMixedComb, subtractFlow);
    gDirectory->GetList()->Delete();

    show_IAA_cent("AR_1509.root", "AR_pp_687.root", "I_AA_Npart_1plus1", 8, 12, pt2Min, pt_assoc_bins[i], pt_assoc_bins[i+1], 6, 7.9, subtractMixedComb, subtractFlow);
    gDirectory->GetList()->Delete();
    
    show_IAA_cent("AR_1509.root", "AR_pp_687.root", "I_AA_Npart_1plus1_away", pt2Min, pt2Min+ptMaxAdd, pt2Min, pt_assoc_bins[i], pt_assoc_bins[i+1], 6, 7.9, subtractMixedComb, subtractFlow);
    gDirectory->GetList()->Delete();

  }

}

void show_IAA_cent(const char* fileName, const char* fileName2, char* name, double pt1Min = 4.0, double pt1Max = 14.0, double pt2Min = 2.0,  double ptAssoc_min = 1.0, double ptAssoc_max = 2.0, Int_t side, Double_t vtx, Int_t subtractMixedComb = 1, Int_t subtractFlow = 1){

  Double_t pt2Max = pt2Min+ptMaxAdd;

  Double_t saveGraph_min = 0.;
  Double_t saveGraph_max = 3.;

  static const Int_t bins = 6;

  Double_t content[bins];
  Double_t error[bins];

  Double_t multBin[bins+1] = {0, 5, 7, 8, 9, 10, 11};
  Double_t multBin_center[bins] = {2.5, 7.5, 15, 25, 35, 45};
  Double_t multBin_error[bins] = {2.5, 2.5, 5, 5, 5, 5};

  TH1D* pp_plot = (TH1D*)getAnalysis(fileName2, pt1Min, pt1Max, pt2Min, ptAssoc_min, ptAssoc_max, vtx, 1, 1, side, 1, 1, 0, subtractMixedComb, subtractFlow, NULL, NULL)->Clone("pp");

  TCanvas* can_graph_ratio = new TCanvas("can", "can", gBasisSize+50, gBasisSize+50, gBasisSize, gBasisSize);
  pp_plot->DrawCopy();

  Double_t phi_bin = getPhiBinsForAnalysis(ptAssoc_min);

  //Int_t bin_start = near->FindBin(3.1415-1*phi_bin);
  //Int_t bin_end = near->FindBin(3.1415+phi_bin);
  Int_t bin_start = pp_plot->FindBin(-1*phi_bin);
  Int_t bin_end = pp_plot->FindBin(+phi_bin);
 
  Double_t error_pp; 
  Double_t content_pp = pp->IntegralAndError(bin_start, bin_end, error_pp);

  pp->~TH1D();//without this destructor the results are wrong if show_IAA_cent is called multiple times

  for(int i=0; i<=bins; i++){
    TH1D* pbpb = (TH1D*)getAnalysis(fileName, pt1Min, pt1Max, pt2Min, ptAssoc_min, ptAssoc_max, vtx, multBin[i]+1, multBin[i+1], side, 1, 1, 0, subtractMixedComb, subtractFlow, NULL, NULL)->Clone("PbPb");
    if(i==0){
      pbpb->SetLineColor(kRed);
      pbpb->DrawCopy("same");
    }
    content[i] = pbpb->IntegralAndError(bin_start, bin_end, error[i]);

    error[i] = TMath::Sqrt(TMath::Power(error[i]/content_pp, 2) + TMath::Power(content[i]/(content_pp*content_pp)*error_pp, 2));
    content[i] /= content_pp;

    pbpb->~TH1D();
  }

  //TGraphErrors* graph = new TGraphErrors(pt_assoc_bins_number-1, pt_assoc_bin_center, content, pt_assoc_bin_error, error);
  TGraphErrors* graph = new TGraphErrors(bins, multBin_center, content, multBin_error, error);
 
  graph->GetXaxis()->SetTitle("centrality");
  graph->GetYaxis()->SetTitle("I_{AA}");
  graph->SetMarkerSize(1);
  graph->SetLineWidth(3);
  graph->SetMarkerColor(kBlue);
  graph->SetLineColor(kBlue);
  graph->SetMarkerStyle(20);

  TCanvas* can_graph_ratio = new TCanvas("can_saveGraph2", "can_saveGraph2", gBasisSize+50, gBasisSize+50, gBasisSize, gBasisSize);
  graph->Draw("AP");

  can_graph_ratio->SetGridx();
  can_graph_ratio->SetGridy();
  can_graph_ratio->SaveAs(Form("pt_spectra/I_AA/%s_%.0f_%.0f.eps", name, ptAssoc_min, ptAssoc_max));
  can_graph_ratio->SaveAs(Form("pt_spectra/I_AA/%s_%.0f_%.0f.C", name, ptAssoc_min, ptAssoc_max));
}



void peakDifference_draw(const char* fileName, double pt1Min = 4.0, double pt1Max = 14.0, double pt2Min = 2.0, Int_t yPos = 0,  Int_t subtractMixedComb = 1, Int_t subtractFlow = 1, Int_t mode = 3){

  Double_t pt2Max = pt2Min+ptMaxAdd;
  char* foldername = "H";

  Double_t saveGraph_min = 0.;
  Double_t saveGraph_max = 3.;

  Int_t trigger_central_near = 0;
  Int_t trigger_central_away = 0;
  Int_t trigger_semi_near = 0;
  Int_t trigger_semi_away = 0;

  Int_t trigger_bg_central_near = 0;
  Int_t trigger_bg_central_away = 0;
  Int_t trigger_bg_semi_near = 0;
  Int_t trigger_bg_semi_away = 0;

  Int_t mult_cent_min = 1;
  Int_t mult_cent_max = 6;//9;

  Int_t mult_semi_min = 8;//10;
  Int_t mult_semi_max = 9;//11;

  const char* fileName2 = "AR_pp_768.root";//"AR_pp_789.root";//"AR_pp_768.root";//"AR_pp_758.root";//"AR_pp_757.root";//"AR_pp_687.root";//"AR_pp_671.root";

  Double_t vtx = 7.9;

  //standard+pure background central
  if(mode==1){
    foldername = "D1_central";
    saveGraph_min = 0.;
    saveGraph_max = 10.;
    TGraphErrors* graph_central_near = (TGraphErrors*)(peakDifference_graph(fileName, pt1Min, pt1Max, pt2Min, mult_cent_min, mult_cent_max, 0, yPos, subtractMixedComb, subtractFlow, &trigger_central_near, &trigger_bg_central_near, vtx)->Clone());
    TGraphErrors* graph_central_away = (TGraphErrors*)(peakDifference_graph(fileName, pt1Min, pt1Max, pt2Min, mult_cent_min, mult_cent_max, 1, yPos, subtractMixedComb, subtractFlow, &trigger_central_away, &trigger_bg_central_away, vtx)->Clone());
    TGraphErrors* graph_semi_near = (TGraphErrors*)(peakDifference_graph(fileName, pt1Min, pt1Max, pt2Min, mult_cent_min, mult_cent_max, 6, yPos, subtractMixedComb, subtractFlow, &trigger_semi_near, &trigger_bg_semi_near, vtx)->Clone());
    TGraphErrors* graph_semi_away = (TGraphErrors*)(peakDifference_graph(fileName, pt2Min, pt2Max, pt2Min, mult_cent_min, mult_cent_max, 6, yPos, subtractMixedComb, subtractFlow, &trigger_semi_away, &trigger_bg_semi_away, vtx)->Clone());
  }else if(mode==2){
    //standard+pure background central
    saveGraph_min = 0.;
    saveGraph_max = 10.;
    foldername = "D1_semi";
    TGraphErrors* graph_central_near = (TGraphErrors*)(peakDifference_graph(fileName, pt1Min, pt1Max, pt2Min, mult_semi_min, mult_semi_max, 0, yPos, subtractMixedComb, subtractFlow, &trigger_central_near, &trigger_bg_central_near, vtx)->Clone());
    TGraphErrors* graph_central_away = (TGraphErrors*)(peakDifference_graph(fileName, pt1Min, pt1Max, pt2Min, mult_semi_min, mult_semi_max, 1, yPos, subtractMixedComb, subtractFlow, &trigger_central_away, &trigger_bg_central_away, vtx)->Clone());
    TGraphErrors* graph_semi_near = (TGraphErrors*)(peakDifference_graph(fileName, pt1Min, pt1Max, pt2Min, mult_semi_min, mult_semi_max, 6, yPos, subtractMixedComb, subtractFlow, &trigger_semi_near, &trigger_bg_semi_near, vtx)->Clone());
    TGraphErrors* graph_semi_away = (TGraphErrors*)(peakDifference_graph(fileName, pt2Min, pt2Max, pt2Min, mult_semi_min, mult_semi_max, 6, yPos, subtractMixedComb, subtractFlow, &trigger_semi_away, &trigger_bg_semi_away, vtx)->Clone());
  }else if(mode==3){
    //standard
    foldername = "D2";
    saveGraph_min = 0.0;
    saveGraph_max = 3.0;//5.0;
    TGraphErrors* graph_central_near = (TGraphErrors*)(peakDifference_graph(fileName, pt1Min, pt1Max, pt2Min, mult_cent_min, mult_cent_max, 0, yPos, subtractMixedComb, subtractFlow, &trigger_central_near, &trigger_bg_central_near, vtx, 2)->Clone());
    TGraphErrors* graph_central_away = (TGraphErrors*)(peakDifference_graph(fileName, pt1Min, pt1Max, pt2Min, mult_cent_min, mult_cent_max, 1, yPos, subtractMixedComb, subtractFlow, &trigger_central_away, &trigger_bg_central_away, vtx, 2)->Clone());
    TGraphErrors* graph_semi_near = (TGraphErrors*)(peakDifference_graph(fileName, pt1Min, pt1Max, pt2Min, mult_semi_min, mult_semi_max, 0, yPos, subtractMixedComb, subtractFlow, &trigger_semi_near, &trigger_bg_semi_near, vtx, 2)->Clone());
    TGraphErrors* graph_semi_away = (TGraphErrors*)(peakDifference_graph(fileName, pt1Min, pt1Max, pt2Min, mult_semi_min, mult_semi_max, 1, yPos, subtractMixedComb, subtractFlow, &trigger_semi_away, &trigger_bg_semi_away, vtx, 2)->Clone());
  }else if(mode==4){
    //pure background
    foldername = "D3";
    saveGraph_min = 1.0;
    saveGraph_max = 3.0;//1.3;
    TGraphErrors* graph_central_near = (TGraphErrors*)(peakDifference_graph(fileName, pt1Min, pt1Max, pt2Min, mult_cent_min, mult_cent_max, 6, yPos, subtractMixedComb, subtractFlow, &trigger_central_near, &trigger_bg_central_near, vtx, 2)->Clone());
    TGraphErrors* graph_central_away = (TGraphErrors*)(peakDifference_graph(fileName, pt2Min, pt2Max, pt2Min, mult_cent_min, mult_cent_max, 6, yPos, subtractMixedComb, subtractFlow, &trigger_central_away, &trigger_bg_central_away, vtx, 2)->Clone());
    TGraphErrors* graph_semi_near = (TGraphErrors*)(peakDifference_graph(fileName, pt1Min, pt1Max, pt2Min, mult_semi_min, mult_semi_max, 6, yPos, subtractMixedComb, subtractFlow, &trigger_semi_near, &trigger_bg_semi_near, vtx, 2)->Clone());
    TGraphErrors* graph_semi_away = (TGraphErrors*)(peakDifference_graph(fileName, pt2Min, pt2Max, pt2Min, mult_semi_min, mult_semi_max, 6, yPos, subtractMixedComb, subtractFlow, &trigger_semi_away, &trigger_bg_semi_away, vtx, 2)->Clone());
  }else if(mode==5){
    //pure background
    foldername = "D4";
    saveGraph_min = 0.;//0.8;
    saveGraph_max = 3.0;//1.2;
    TGraphErrors* graph_central_near = (TGraphErrors*)(peakDifference_graph(fileName, pt1Min, pt1Max, pt2Min, mult_cent_min, mult_cent_max, 1, yPos, subtractMixedComb, subtractFlow, &trigger_central_near, &trigger_bg_central_near, vtx, 1)->Clone());
    TGraphErrors* graph_central_away = (TGraphErrors*)(peakDifference_graph(fileName, pt1Min, pt1Max, pt2Min, mult_cent_min, mult_cent_max, 1, yPos, subtractMixedComb, subtractFlow, &trigger_central_away, &trigger_bg_central_away, vtx, 2)->Clone());
    TGraphErrors* graph_semi_near = (TGraphErrors*)(peakDifference_graph(fileName, pt1Min, pt1Max, pt2Min, mult_semi_min, mult_semi_max, 1, yPos, subtractMixedComb, subtractFlow, &trigger_semi_near, &trigger_bg_semi_near, vtx, 1)->Clone());
    TGraphErrors* graph_semi_away = (TGraphErrors*)(peakDifference_graph(fileName, pt1Min, pt1Max, pt2Min, mult_semi_min, mult_semi_max, 1, yPos, subtractMixedComb, subtractFlow, &trigger_semi_away, &trigger_bg_semi_away, vtx, 2)->Clone());
  }
  

  /*
  TGraphErrors* graph_central_near = (TGraphErrors*)(peakDifference_graph(fileName, pt1Min, pt1Max, pt2Min, 1, 5, 8, yPos, subtractMixedComb, subtractFlow, &trigger_central_near, &trigger_bg_central_near)->Clone());
  TGraphErrors* graph_central_away = (TGraphErrors*)(peakDifference_graph(fileName, pt2Min, pt1Min, pt2Min, 1, 5, 6, yPos, subtractMixedComb, subtractFlow, &trigger_central_near, &trigger_bg_central_near)->Clone());
  TGraphErrors* graph_semi_near = (TGraphErrors*)(peakDifference_graph(fileName, pt1Min, pt1Max, pt2Min, 1, 5, 1, yPos, subtractMixedComb, subtractFlow, &trigger_central_near, &trigger_bg_central_near)->Clone());
  TGraphErrors* graph_semi_away = (TGraphErrors*)(peakDifference_graph(fileName, pt2Min, pt1Max, pt2Min, 1, 5, 6, yPos, subtractMixedComb, subtractFlow, &trigger_central_near, &trigger_bg_central_near)->Clone());
  */
  
  /*
  //near side over pure background
  TGraphErrors* graph_central_near = (TGraphErrors*)(peakDifference_graph(fileName, pt1Min, pt1Max, pt2Min, 1, 5, 0, yPos, subtractMixedComb, subtractFlow, &trigger_central_near, &trigger_bg_central_near)->Clone());
  TGraphErrors* graph_central_away = (TGraphErrors*)(peakDifference_graph(fileName, pt1Min, pt1Max, pt2Min, 1, 5, 6, yPos, subtractMixedComb, subtractFlow, &trigger_central_away, &trigger_bg_central_away)->Clone());

  TGraphErrors* graph_semi_near = (TGraphErrors*)(peakDifference_graph(fileName, pt1Min, pt1Max, pt2Min, 9, 10, 0, yPos, subtractMixedComb, subtractFlow, &trigger_semi_near, &trigger_bg_semi_near)->Clone());
  TGraphErrors* graph_semi_away = (TGraphErrors*)(peakDifference_graph(fileName, pt1Min, pt1Max, pt2Min, 9, 10, 6, yPos, subtractMixedComb, subtractFlow, &trigger_semi_away, &trigger_bg_semi_away)->Clone());
  */
  /*
  //away side over pure background
  TGraphErrors* graph_central_near = (TGraphErrors*)(peakDifference_graph(fileName, pt1Min, pt1Max, pt2Min, 1, 5, 1, yPos, subtractMixedComb, subtractFlow, &trigger_central_near, &trigger_bg_central_near)->Clone());
  //TGraphErrors* graph_central_away = (TGraphErrors*)(peakDifference_graph(fileName, pt2Min, pt1Min, pt2Min, 1, 5, 6, yPos, subtractMixedComb, subtractFlow, &trigger_central_away, &trigger_bg_central_away)->Clone());
  TGraphErrors* graph_central_away = (TGraphErrors*)(peakDifference_graph(fileName, pt1Min, pt1Max, pt2Min, 1, 5, 6, yPos, subtractMixedComb, subtractFlow, &trigger_central_away, &trigger_bg_central_away)->Clone());

  TGraphErrors* graph_semi_near = (TGraphErrors*)(peakDifference_graph(fileName, pt1Min, pt1Max, pt2Min, 9, 10, 1, yPos, subtractMixedComb, subtractFlow, &trigger_semi_near, &trigger_bg_semi_near)->Clone());
  //TGraphErrors* graph_semi_away = (TGraphErrors*)(peakDifference_graph(fileName, pt2Min, pt1Min, pt2Min, 9, 10, 6, yPos, subtractMixedComb, subtractFlow, &trigger_semi_away, &trigger_bg_semi_away)->Clone());
  TGraphErrors* graph_semi_away = (TGraphErrors*)(peakDifference_graph(fileName, pt1Min, pt1Max, pt2Min, 9, 10, 6, yPos, subtractMixedComb, subtractFlow, &trigger_semi_away, &trigger_bg_semi_away)->Clone());
  */
  /*
  //mixed comb as signal
  TGraphErrors* graph_central_near = (TGraphErrors*)(peakDifference_graph(fileName, pt1Min, pt1Max, pt2Min, 1, 5, 4, yPos, subtractMixedComb, subtractFlow, &trigger_central_near, &trigger_bg_central_near)->Clone());
  TGraphErrors* graph_central_away = (TGraphErrors*)(peakDifference_graph(fileName, pt1Min, pt1Max, pt2Min, 1, 5, 5, yPos, subtractMixedComb, subtractFlow, &trigger_central_away, &trigger_bg_central_away)->Clone());
  TGraphErrors* graph_semi_near = (TGraphErrors*)(peakDifference_graph(fileName, pt1Min, pt1Max, pt2Min, 9, 10, 4, yPos, subtractMixedComb, subtractFlow, &trigger_semi_near, &trigger_bg_semi_near)->Clone());
  TGraphErrors* graph_semi_away = (TGraphErrors*)(peakDifference_graph(fileName, pt1Min, pt1Max, pt2Min, 9, 10, 5, yPos, subtractMixedComb, subtractFlow, &trigger_semi_away, &trigger_bg_semi_away)->Clone());
  */
  /*
  //yield from background same
  TGraphErrors* graph_central_near = (TGraphErrors*)(peakDifference_graph(fileName, pt1Min, pt1Max, pt2Min, 1, 5, 7, yPos, subtractMixedComb, subtractFlow, &trigger_central_near, &trigger_bg_central_near)->Clone());
  TGraphErrors* graph_central_away = (TGraphErrors*)(peakDifference_graph(fileName, pt1Min, pt1Max, pt2Min, 1, 5, 8, yPos, subtractMixedComb, subtractFlow, &trigger_central_away, &trigger_bg_central_away)->Clone());

  TGraphErrors* graph_semi_near = (TGraphErrors*)(peakDifference_graph(fileName, pt1Min, pt1Max, pt2Min, 9, 10, 7, yPos, subtractMixedComb, subtractFlow, &trigger_semi_near, &trigger_bg_semi_near)->Clone());
  TGraphErrors* graph_semi_away = (TGraphErrors*)(peakDifference_graph(fileName, pt1Min, pt1Max, pt2Min, 9, 10, 8, yPos, subtractMixedComb, subtractFlow, &trigger_semi_away, &trigger_bg_semi_away)->Clone());
  */

  /*
    //mixed combinatorics pure background 
  TGraphErrors* graph_central_near = (TGraphErrors*)(peakDifference_graph(fileName, pt1Min, pt1Max, pt2Min, 1, 5, 4, yPos, subtractMixedComb, subtractFlow, &trigger_central_near, &trigger_bg_central_near)->Clone());
  TGraphErrors* graph_central_away = (TGraphErrors*)(peakDifference_graph(fileName, pt1Min, pt1Max, pt2Min, 1, 5, 5, yPos, subtractMixedComb, subtractFlow, &trigger_central_away, &trigger_bg_central_away)->Clone());

  TGraphErrors* graph_semi_near = (TGraphErrors*)(peakDifference_graph(fileName, pt1Min, pt1Max, pt2Min, 9, 10, 4, yPos, subtractMixedComb, subtractFlow, &trigger_semi_near, &trigger_bg_semi_near)->Clone());
  TGraphErrors* graph_semi_away = (TGraphErrors*)(peakDifference_graph(fileName, pt1Min, pt1Max, pt2Min, 9, 10, 5, yPos, subtractMixedComb, subtractFlow, &trigger_semi_away, &trigger_bg_semi_away)->Clone());
  */
  


  //write out trigger info
  TCanvas* trigger_info = new TCanvas(Form("trigger info for %.1f > p\\_{T,1} > %.1f,  p\\_{T,12} > %.1f ", pt1Min, pt1Max, pt2Min), Form("trigger info for %.1f > p\\_{T,1} > %.1f,  p\\_{T,12} > %.1f ", pt1Min, pt1Max, pt2Min));
  Double_t rate_central_near = ((Double_t)trigger_central_near)/(trigger_central_near+trigger_bg_central_near);
  TLatex * trigger_value_central_near = new TLatex(0.05,0.90, Form("central near: %i, background: %i, trigger rate %.2f ", trigger_central_near, trigger_bg_central_near, rate_central_near));
  trigger_value_central_near->SetTextFont(62);
  trigger_value_central_near->SetTextSize(0.05);
  trigger_value_central_near->Draw();
  Double_t rate_central_away = ((Double_t)trigger_central_away)/(trigger_central_away+trigger_bg_central_away);
  TLatex * trigger_value_central_away = new TLatex(0.05,0.80, Form("central away: %i, background: %i, trigger rate %.2f  ", trigger_central_away, trigger_bg_central_away, rate_central_away));
  trigger_value_central_away->SetTextFont(62);
  trigger_value_central_away->SetTextSize(0.05);
  trigger_value_central_away->Draw();
  Double_t rate_semi_near = ((Double_t)trigger_semi_near)/(trigger_semi_near+trigger_bg_semi_near);
  TLatex * trigger_value_semi_near = new TLatex(0.05,0.60, Form("semi central near: %i, background: %i, trigger rate %.2f  ", trigger_semi_near, trigger_bg_semi_near, rate_semi_near));
  trigger_value_semi_near->SetTextFont(62);
  trigger_value_semi_near->SetTextSize(0.05);
  trigger_value_semi_near->Draw();
  Double_t rate_semi_away = ((Double_t)trigger_semi_away)/(trigger_semi_away+trigger_bg_semi_away);
  TLatex * trigger_value_semi_away = new TLatex(0.05,0.50, Form("semi central away: %i, background: %i, trigger rate %.2f  ", trigger_semi_away, trigger_bg_semi_away, rate_semi_away));
  trigger_value_semi_away->SetTextFont(62);
  trigger_value_semi_away->SetTextSize(0.05);
  trigger_value_semi_away->Draw();
  trigger_info->SaveAs(Form("pt_spectra/%s/trigger_info_%.0f_%.0f_%.0f.eps", foldername, pt1Min, pt1Max, pt2Min));
  //trigger info is saved


  TCanvas* can_graph = new TCanvas(Form("result %i", yPos), Form("result %i", yPos), gBasisSize+50, yPos*(gBasisSize+50), gBasisSize, gBasisSize);

  Int_t elements = pt_assoc_bins_number;
  
  Double_t maximum = 0;
  Double_t max_hard = 1000000000;
  Double_t max  = TMath::MaxElement(elements, graph_central_near->GetY());
  if(max>maximum && max<max_hard) maximum = max;
  Printf("maximum of central near is %f ", max);
  max  = TMath::MaxElement(elements, graph_central_away->GetY());
  Printf("maximum of central away is %f ", max);
  if(max>maximum && max<max_hard) maximum = max;
  max = TMath::MaxElement(elements, graph_semi_near->GetY());
  Printf("maximum of semi near is %f ", max);
  if(max>maximum && max<max_hard) maximum = max;
  max = TMath::MaxElement(elements, graph_semi_away->GetY());
  Printf("maximum of semi away is %f ", max);
  if(max>maximum && max<max_hard) maximum = max;
  graph_central_near->SetMaximum(maximum*1.25);
  
  Double_t minimum = maximum;

  Double_t min = TMath::MinElement(elements, graph_central_near->GetY());
  Printf("minimum of central near is %f ", minimum);
  if(min>0 && min<minimum) minimum = min;
  min = TMath::MinElement(elements, graph_central_away->GetY());
  Printf("minimum of central away is %f ", min);
  if(min>0 && min<minimum) minimum = min;
  min = TMath::MinElement(elements, graph_semi_near->GetY());
  Printf("minimum of semi near is %f ", min);
  if(min>0 && min<minimum) minimum = min;
  min = TMath::MinElement(elements, graph_semi_away->GetY());
  Printf("minimum of semi away is %f ", min);
  if(min>0 && min<minimum) minimum = min;

  if(minimum<0.00000000001)
    minimum = 0.01;//0.0002;
 
  graph_central_near->SetMinimum((double)minimum/1.30);
  
    
  //graph_central_near->SetTitle(Form("p_{T,assoc} spectrum for %.1f < p_{T,1} < %.1f and %.1f < p_{T,2} < %.1f", pt1Min, pt1Max, pt2Min, pt2Max));
  graph_central_near->SetTitle("");
  graph_central_near->Draw("AP");

  graph_central_away->SetMarkerColor(kRed);
  graph_central_away->SetLineColor(kRed);
  graph_central_away->Draw("P");
  
  graph_semi_near->SetMarkerColor(kCyan);
  graph_semi_near->SetLineColor(kCyan);
  graph_semi_near->Draw("P");

  graph_semi_away->SetMarkerColor(kOrange);
  graph_semi_away->SetLineColor(kOrange);
  graph_semi_away->Draw("P");

  gPad->SetLogy();

  TLegend* leg = getLegend(0.30,0.73,0.85,0.88);
  
  leg->AddEntry(graph_central_near,"central T1 assoc","l");
  leg->AddEntry(graph_central_away,"central T2 assoc","l");
  leg->AddEntry(graph_semi_near,"semi central T1 assoc","l");
  leg->AddEntry(graph_semi_away,"semi central T2 assoc","l");
  
  /*
  leg->AddEntry(graph_central_near,"cent bin counting","l");
  leg->AddEntry(graph_central_away,"cent fit function","l");
  leg->AddEntry(graph_semi_near,"semi bin counting","l");
  leg->AddEntry(graph_semi_away,"semi fit function","l");
  */
  /*
  leg->AddEntry(graph_central_near,"cent bin counting","l");
  leg->AddEntry(graph_central_away,"cent fit function","l");
  leg->AddEntry(graph_semi_near,"semi bin counting","l");
  leg->AddEntry(graph_semi_away,"semi fit function","l");
  */
  /*
  leg->AddEntry(graph_central_near,"T1 bin counting","l");
  leg->AddEntry(graph_central_away,"T1 fit function","l");
  leg->AddEntry(graph_semi_near,"T2 bin counting","l");
  leg->AddEntry(graph_semi_away,"T2 fit function","l");
  */
  //leg->AddEntry(graph_semi_near,"pp T1 assoc","l");
  //leg->AddEntry(graph_semi_away,"pp T2 assoc","l");
  //leg->AddEntry(graph_semi_near,"pure BG central near","l");
  //leg->AddEntry(graph_semi_away,"pure BG central away","l");

  /*
  leg->AddEntry(graph_central_near,"background same","l");
  leg->AddEntry(graph_central_away,"pure Background","l");
  leg->AddEntry(graph_semi_near,"2+1 no BG Correction","l");
  leg->AddEntry(graph_semi_away,"pure Background near side","l");
  */
  leg->Draw("same");
  
  writePoints(graph_central_near, Form("pt_spectra/%s/central_near_%.0f_%.0f_%.0f", foldername, pt1Min, pt1Max, pt2Min));
  writePoints(graph_central_away, Form("pt_spectra/%s/central_away_%.0f_%.0f_%.0f", foldername, pt1Min, pt1Max, pt2Min));
  writePoints(graph_semi_near, Form("pt_spectra/%s/semi_near_%.0f_%.0f_%.0f", foldername, pt1Min, pt1Max, pt2Min));
  writePoints(graph_semi_away, Form("pt_spectra/%s/semi_away_%.0f_%.0f_%.0f", foldername, pt1Min, pt1Max, pt2Min));

  can_graph->SaveAs(Form("pt_spectra/%s/pt_spectrum_%.0f_%.0f_%.0f.eps", foldername, pt1Min, pt1Max, pt2Min));
  can_graph->SaveAs(Form("pt_spectra/%s/pt_spectrum_%.0f_%.0f_%.0f.C", foldername, pt1Min, pt1Max, pt2Min));
  
  save_graph_ratio(graph_central_near, graph_central_away, graph_semi_near, graph_semi_away, elements, "\\frac{T1 assoc}{T2 assoc}", Form("pt_spectra/%s/centralANDsemi_near_away_ratio_%.0f_%.0f_%.0f", foldername, pt1Min, pt1Max, pt2Min), "central", "semi central", saveGraph_min, saveGraph_max, pt2Min);
  //save_graph_ratio(graph_central_near, graph_central_away, graph_semi_near, graph_semi_away, elements, Form("central and semi near over away p_{T,assoc} spectrum for %.1f < p_{T,1} < %.1f and %.1f < p_{T,2} < %.1f", pt1Min, pt1Max, pt2Min, pt2Max), Form("pt_spectra/%s/centralANDsemi_near_away_ratio_%.0f_%.0f_%.0f", foldername, pt1Min, pt1Max, pt2Min), "T1 assoc", "T2 assoc", saveGraph_min, saveGraph_max, pt2Min);
  // //save_graph_ratio(graph_central_near, graph_central_away, graph_semi_near, graph_semi_away, elements, Form("central and semi near over away p_{T,assoc} spectrum for %.1f < p_{T,1} < %.1f and %.1f < p_{T,2} < %.1f", pt1Min, pt1Max, pt2Min, pt2Max), Form("pt_spectra/%s/centralANDsemi_near_away_ratio_%.0f_%.0f_%.0f", foldername, pt1Min, pt1Max, pt2Min), "PbPb", "pp", saveGraph_min, saveGraph_max, pt2Min);

 save_graph_ratio(graph_central_near, graph_semi_near, graph_central_away, graph_semi_away, elements, "I_{CP}", Form("pt_spectra/%s/nearANDaway_central_semi_ratio_%.0f_%.0f_%.0f", foldername, pt1Min, pt1Max, pt2Min), "T1 associated", "T2 associated", saveGraph_min, saveGraph_max, pt2Min);

  //mode 5
  //save_graph_ratio(graph_central_near, graph_semi_near, graph_central_away, graph_semi_away, elements, Form("near and semi central over semi p_{T,assoc} spectrum for %.1f < p_{T,1} < %.1f and %.1f < p_{T,2} < %.1f", pt1Min, pt1Max, pt2Min, pt2Max), Form("pt_spectra/%s/nearANDaway_central_semi_ratio_%.0f_%.0f_%.0f", foldername, pt1Min, pt1Max, pt2Min), "bin counting", "fit function", saveGraph_min, saveGraph_max, pt2Min);

  //save_graph_difference(graph_central_near, graph_central_away, graph_semi_near, graph_semi_away, elements, Form("central and semi near minus away p_{T,assoc} spectrum for %.1f < p_{T,1} < %.1f and %.1f < p_{T,2} < %.1f", pt1Min, pt1Max, pt2Min, pt2Max), Form("pt_spectra/%s/centralANDsemi_near_away_diff_%.0f_%.0f_%.0f.eps", foldername, pt1Min, pt1Max, pt2Min));
  
}


void writePoints(TGraphErrors* graph, char* name){
   //write out trigger info

  FILE *points; 
  points = fopen(name, "w");

  for(int i=0; i<pt_assoc_bins_number-1; i++){
    Double_t content_first;
    Double_t content_second;
    Double_t error_first;
    Double_t error_second;
    graph->GetPoint(i, content_first, content_second);
    error_first = graph->GetErrorX(i);
    error_second = graph->GetErrorY(i);

    fprintf(points, Form("%f +/- %f, %f +/- %f \n", content_first, error_first, content_second, error_second));
  }
 
  fclose(points);


}

TGraphErrors* save_graph_ratio(TGraphErrors* first, TGraphErrors* second, TGraphErrors* third, TGraphErrors* forth, const Int_t bins, char* title, char* name, char* legend1, char* legend2, double min, double max, Double_t pt2Min){

  Double_t start_x = 0.45;
  if(legend1=="near")
    start_x = 0.70;

  TLegend *leg  = getLegend(start_x,0.75,0.85,0.9);

  TGraphErrors* graph = (TGraphErrors*)(save_graph_compute(first, second, bins, title, 0, NULL, NULL)->Clone());
  if(third!=NULL && forth!=NULL){
    TGraphErrors* graph2 = (TGraphErrors*)(save_graph_compute(third, forth, bins, title, 0, NULL, NULL)->Clone());

    leg->AddEntry(graph, legend1, "p");
    leg->AddEntry(graph2, legend2, "p");
    
    save_graph(graph, graph2, name, 0, leg, min, max, pt2Min);
  }else{
    save_graph(graph, NULL, name, 0, leg, min, max, pt2Min);
  }

  return graph;
  
}


void save_graph_difference(TGraphErrors* first, TGraphErrors* second, TGraphErrors* third, TGraphErrors* forth, const Int_t bins, char* title, char* name){

  Double_t diff1 = 0;
  Double_t diff1_err = 0;
  Double_t diff2 = 0;
  Double_t diff2_err = 0;

  TGraphErrors* graph = (TGraphErrors*)(save_graph_compute(first, second, bins, title, 1, &diff1, &diff1_err)->Clone());
  TGraphErrors* graph2 = (TGraphErrors*)(save_graph_compute(third, forth, bins, title, 1, &diff2, &diff2_err)->Clone());

  TLegend *leg  = getLegend(0.45,0.15,0.85,0.3);
  leg->AddEntry(graph, "central", "p");
  leg->AddEntry(graph2, "semi central", "p");
  leg->AddEntry(graph,  Form("\\Delta E = %.1f \\pm %.1f", diff1, diff1_err),"");
  leg->AddEntry(graph2,  Form("\\Delta E = %.1f \\pm %.1f", diff2, diff2_err),"");

  save_graph(graph, graph2, name, 1, leg, 0, 0);

}


//mode: 0 ratio of the graphs
//mode: 1 difference of the graphs
//returns difference of the graphs
void save_graph(TGraphErrors* first, TGraphErrors* second, char* name, int mode, TLegend *leg, double min, double max, Double_t pt2Min){
  TCanvas* can_graph_ratio = new TCanvas("can_saveGraph", "can_saveGraph", gBasisSize+50, gBasisSize+50, gBasisSize, gBasisSize);
  if(mode==0){
    //first->SetMinimum(0.4);
    //first->SetMaximum(1.6);

    first->SetMinimum(min);
    first->SetMaximum(max);

    //first->SetMinimum(0.2);
    //first->SetMaximum(1.8);
  }else if(mode==1){

    //first->SetMinimum(-0.9);
    //first->SetMaximum(0.5);

    first->SetMinimum(-0.4);
    first->SetMaximum(0.3);

  }
  first->SetMarkerSize(1);
  first->SetLineWidth(3);
  first->SetMarkerColor(kBlue);
  first->SetLineColor(kBlue);
  first->SetMarkerStyle(20);
  first->Draw("AP");
  writePoints(first, name);

  Double_t maximum = 0.067;
  TBox *box = new TBox(0.5,1.0-maximum,8,1.0+maximum);
  //box->SetLineColor(kRed);
  box->SetFillColor(kGreen);
  box->Draw();
  first->Draw("P");

  
  TLine *line = new TLine(pt2Min,min,pt2Min,max);
  line->SetLineStyle(2);
  line->SetLineColor(kRed);
  line->Draw();

  if(mode==0){
    /*
    TLatex * cent = new TLatex(5.5,0.35, "central events 0-7.5%");
    cent->SetTextFont(62);
    cent->SetTextSize(0.04);
    cent->Draw();
    TLatex * semi = new TLatex(5.5,0.10, "semi events 30-50%");
    semi->SetTextFont(62);
    semi->SetTextSize(0.04);
    semi->Draw();
    */
    TLatex * semi = new TLatex(5.5,0.15, "\\frac{events \\ at \\ 0-7.5}{events \\ at \\ 30-50}");//TODO finish this information, the percentage sign is not shown yet
    semi->SetTextFont(62);
    semi->SetTextSize(0.05);
    //semi->Draw();
  }
  

  if(second!=NULL){
    second->SetMarkerSize(1);
    second->SetLineWidth(3);
    second->SetMarkerColor(kRed);
    second->SetLineColor(kRed);
    second->SetMarkerStyle(20);
    second->Draw("P");
    writePoints(second, Form("%s_second", name));
  }

  leg->Draw("same");
 
  /*
  TLine *line = new TLine(0,maximum,8,maximum);
  line->SetLineColor(kRed);
  line->Draw();
  */

  can_graph_ratio->SetGridx();
  can_graph_ratio->SetGridy();
  can_graph_ratio->SaveAs(Form("%s.eps", name));
  can_graph_ratio->SaveAs(Form("%s.C", name));
}

TGraphErrors* save_graph_compute(TGraphErrors* first, TGraphErrors* second, Int_t bins, char* title, int mode, Double_t* diff, Double_t* diff_err){

  Int_t used_bins = bins;

  Printf("bins %i ", bins);

  //for not identifiable reasons the use of bins doesn't work and produced an array out of index error

  //Double_t content_x[bins];
  //Double_t content_y[bins];
  //Double_t x_error[bins];
  //Double_t y_error[bins];

  Double_t content_x[pt_assoc_bins_number];
  Double_t content_y[pt_assoc_bins_number];
  Double_t x_error[pt_assoc_bins_number];
  Double_t y_error[pt_assoc_bins_number];

  for(int i=0; i<bins; i++){
    content_x[i] = 0;
    content_y[i] = 0;
    x_error[i] = 0;
    y_error[i] = 0;

    Double_t error_first = 0;
    Double_t error_second = 0;
    first->GetPoint(i, content_x[i], content_y[i]);
    x_error[i] = first->GetErrorX(i);
    error_first = first->GetErrorY(i);
    Double_t content_second = 0;
    second->GetPoint(i, content_x[i], content_second);
    error_second = second->GetErrorY(i);
    if(content_second>0 && error_second>0){
      if(mode==0){
	content_y[i] /= content_second;
	y_error[i] = TMath::Sqrt(TMath::Power(error_first/content_second, 2) + TMath::Power(error_second/content_second*content_y[i], 2));
      }else if(mode==1){
	content_y[i] -= content_second;
	y_error[i] = TMath::Sqrt(TMath::Power(error_first, 2) + TMath::Power(error_second, 2));
      }else{
	Printf("calls save graph with no correct mode");
	return;
      }
    }else{
      Printf("content (%f) or content error (%f) is smaller or equal to 0", content_second, error_second); 
      content_second = 0.; 
      error_second = 0.;
      used_bins--;
    }
  }
   //if this is mode 1 (difference of the graphs) this calculates the difference of the plots
  //Double_t diff = 0;
  //Double_t diff_err = 0;
  if(mode==1){
    for(int i=0; i<bins; i++){
      *diff += content_y[i] * content_x[i];
      *diff_err += TMath::Power(y_error[i]*content_x[i], 2);
    }
    *diff_err = TMath::Sqrt(*diff_err);

  }
  TGraphErrors* graph = new TGraphErrors(used_bins, content_x, content_y, x_error, y_error);
  //graph->SetTitle(title);
  graph->SetTitle("");
  graph->GetYaxis()->SetTitle(title);
  graph->GetXaxis()->SetTitle("p_{T,assoc} (GeV/c)");

  return graph;
}

void printPeakYield(const char* fileName, double pt1Min = 4.0, double pt1Max = 14.0, double pt2Min = 2.0, double ptAssocMin = 1.0, double ptAssocMax = 2.0, Int_t multBinBegin = 1, Int_t multBinEnd = 5, Int_t side = 0, Int_t subtractMixedComb = 1, Int_t subtractFlow = 1, Int_t* trigger = NULL, Int_t* triggerBackground = NULL, double vtx = 9.9, Int_t method = 2){

  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadRightMargin(0.07);

  Double_t error;

  Double_t content =     peakYield(fileName, pt1Min, pt1Max, pt2Min, ptAssocMin, ptAssocMax, multBinBegin, multBinEnd, side, subtractMixedComb, subtractFlow, trigger, triggerBackground, vtx, &error, method);

  printf("content: %f +/- %f\n", content, error);
}

Double_t peakYield(const char* fileName, double pt1Min = 4.0, double pt1Max = 14.0, double pt2Min = 2.0, double ptAssocMin = 1.0, double ptAssocMax = 2.0, Int_t multBinBegin = 1, Int_t multBinEnd = 5, Int_t side = 0, Int_t subtractMixedComb = 1, Int_t subtractFlow = 1, Int_t* trigger = NULL, Int_t* triggerBackground = NULL, double vtx = 9.9, Double_t* error = NULL, Int_t method = 2){

  if(side==1 && pt2Min+ptMaxAdd<=ptAssocMin){
    *error = 0;
    return 0;
  }
  if(side==6 && pt1Max<=ptAssocMin){
    *error = 0;
    return 0;
  }
  if(side==0 && ptAssocMin>=pt1Min){
    *error = 0;
    return 0;
  }
  if(side==1 && ptAssocMin>=pt2Min){
    *error = 0;
    return 0;
  }

  /*
  if(side==1 && pt2Min+0.5<=ptAssocMin){
    *error = 0;
    return 0;
  }
  */
  Printf("analyze next: pT assoc : %f,  for side %i ", ptAssocMin, side);

  TH1D* near = (TH1D*)getAnalysis(fileName, pt1Min, pt1Max, pt2Min, ptAssocMin, ptAssocMax, vtx, multBinBegin, multBinEnd, side, 1, 1, 0, subtractMixedComb, subtractFlow, trigger, triggerBackground, method)->Clone();
    
  if(near==NULL){
    *error = 0;
    return 0;
  }

  Int_t zero_bin = near->FindBin(-0.01);
  Printf("#############################################################");
  Printf("pT assoc : %f,  yield %f ", ptAssocMin, near->GetBinContent(zero_bin));
  Printf("#############################################################");

  TCanvas* can = new TCanvas("can_filename", "can_filename");
  near->SetTitle("");
  
  //show more from the top of the plot to be able to see more from the fit and to have space for the fit information
  double max = near->GetMaximum();
  double min = near->GetMinimum();
  near->SetMaximum(max*1.25-0.25*min);
  near->DrawCopy();

  if(method==1){
    Double_t phi_bin = getPhiBinsForAnalysis(ptAssocMin);

    //Int_t bin_start = near->FindBin(3.1415-1*phi_bin);
    //Int_t bin_end = near->FindBin(3.1415+phi_bin);
    Int_t bin_start = near->FindBin(-1*phi_bin);
    Int_t bin_end = near->FindBin(phi_bin);

    Double_t error2[2];
    /*
    Double_t content = near->IntegralAndError(bin_start, bin_end, error2[0])*near->GetBinWidth(bin_start);
    if(error)
      *error = error2[0]*near->GetBinWidth(bin_start);
      */
    Double_t content = near->IntegralAndError(bin_start, bin_end, error2[0]);
    if(error)
      *error = error2[0];


    near->~TH1D();
    
    TLatex* latex_par0;
    if(*error>3.0)
      latex_par0 = new TLatex(0.7,1.17*max-0.17*min,Form("b = %.0f \\pm %.0f", content, *error));
    else if(*error>0.3)
      latex_par0 = new TLatex(0.7,1.17*max-0.17*min,Form("b = %.1f \\pm %.1f", content, *error));
    else
      latex_par0 = new TLatex(0.7,1.17*max-0.17*min,Form("b = %.2f \\pm %.2f", content, *error));


    latex_par0->SetTextFont(62);
    latex_par0->SetTextSize(0.05);
    latex_par0->Draw();
    //can->SaveAs(Form("pt_spectra/D2/plots/plot_%.1f_%.1f_%.1f_%.1f_%.1f_side_%i_mult_%i_%i.eps", pt1Min, pt1Max, pt2Min, ptAssocMin, ptAssocMax, side, multBinBegin, multBinEnd));

    return content;
  }

  Double_t fit_limit = 1.4;//1.0;//0.8;

  //TF1* fit = new TF1("fit", "[0]+[1]*exp(-0.5*(x/[2])**2)", -1.0*fit_limit, 1.0*fit_limit);
  TF1* fit = new TF1("fit", "[0]+[1]/([2]*2.506628275)*exp(-0.5*(x/[2])**2)", -1.0*fit_limit, 1.0*fit_limit);
  Bool_t fit_accepted = kFALSE;
  Double_t chiSquare = 0;

    Double_t par0 = 0;
    Double_t par0Err = 0;
    Double_t par1 = 0;
    Double_t par1Err = 0;
    Double_t par2 = 0;
    Double_t par2Err = 0;

    Double_t result = 0;

    fit->SetParameters(max, 0.1, 0.3);
    Double_t start_sigma = 0.2;
    fit->SetParLimits(2, 0.07, 0.5);

    Int_t corrections = assignNotZeroError(near);

    Double_t peakhight = (max-min)/(start_sigma*2.506628275);
    fit->SetParameters(near->GetMinimum(), peakhight, start_sigma);
    TFitResultPtr res = near->Fit("fit", "MSI", "", -1.0*fit_limit, 1.0*fit_limit);
    Double_t binWidth = near->GetBinWidth(near->FindBin(0));
    par0 = fit->GetParameter(0);
    par0Err = fit->GetParError(0);
    par1 = fit->GetParameter(1);
    par1Err = fit->GetParError(1);
    par2 = fit->GetParameter(2);
    par2Err = fit->GetParError(2);
    
    result = par1/binWidth;


    /*
    while(!fit_accepted){
      //fit->SetParameters(near->GetMinimum(), near->GetMaximum()-near->GetMinimum(), start_sigma);
      Double_t peakhight = (max-min)/(start_sigma*2.506628275);
      fit->SetParameters(near->GetMinimum(), peakhight, start_sigma);
    
      //h1_phi->Fit("fit", "N", "", -TMath::Pi()/2+0.01, -1.0);
      TFitResultPtr res = near->Fit("fit", "MSI", "", -1.0*fit_limit, 1.0*fit_limit);
      Int_t fitstatus = res;
      Double_t newChiSquare = fit->GetChisquare();
      Printf("fitstatus: %i, chi square: %f", fitstatus, chiSquare);
      //near->Fit("fit", "MI", "", -1.0*fit_limit, 1.0*fit_limit);
      if(chiSquare==0||newChiSquare<chiSquare){
	par0 = fit->GetParameter(0);
	par0Err = fit->GetParError(0);
	par1 = fit->GetParameter(1);
	par1Err = fit->GetParError(1);
	par2 = fit->GetParameter(2);
	par2Err = fit->GetParError(2);
	chiSquare = newChiSquare;
      }
      
      //result = par1*par2*TMath::Sqrt(2*TMath::Pi());//sqrt 2 pi is the integral of the normal distribution
      result = par1;
      / *
      if(par2==0.5||par2==0.07||(result<max-min)/10){
	start_sigma += 0.1;
      }else
	fit_accepted = kTRUE;
      * /
      start_sigma+= 0.1;
      if(start_sigma>0.45)
	fit_accepted = kTRUE;
    }
    */
  
    TLatex* latex_par0;
    TLatex* latex_par1;
    
    if(par0Err>3.0)
      latex_par0 = new TLatex(0.7,1.17*max-0.17*min,Form("a = %.0f \\pm %.0f", par0, par0Err));
    else if(par0Err>0.3)
      latex_par0 = new TLatex(0.7,1.17*max-0.17*min,Form("a = %.1f \\pm %.1f", par0, par0Err));
    else
      latex_par0 = new TLatex(0.7,1.17*max-0.17*min,Form("a = %.2f \\pm %.2f", par0, par0Err));

    if(par1Err>3.0)
      latex_par1 = new TLatex(0.7,1.09*max-0.09*min,Form("b = %.0f \\pm %.0f", par1, par1Err));
    else if(par1Err>0.3)
      latex_par1 = new TLatex(0.7,1.09*max-0.09*min,Form("b = %.1f \\pm %.1f", par1, par1Err));
    else
      latex_par1 = new TLatex(0.7,1.09*max-0.09*min,Form("b = %.2f \\pm %.2f", par1, par1Err));
    TLatex* latex_par2 = new TLatex(0.7,1.01*max-0.01*min,Form("\\sigma = %.2f \\pm %.2f", par2, par2Err));

    TLatex* latex_par3 = new TLatex(0.7,0.93*max+0.07*min,Form("%i corrections", corrections));

    TLatex* latex_par4;
    if(par1Err/binWidth>3.0)
      latex_par4 = new TLatex(0.7,0.85*max+0.15*min,Form("I = %.0f \\pm %.0f", par1/binWidth, par1Err/binWidth));
    else if(par1Err/binWidth>0.3)
      latex_par4 = new TLatex(0.7,0.85*max+0.15*min,Form("I = %.1f \\pm %.1f", par1/binWidth, par1Err/binWidth));
    else
      latex_par4 = new TLatex(0.7,0.85*max+0.15*min,Form("I = %.2f \\pm %.2f", par1/binWidth, par1Err/binWidth));

    latex_par0->SetTextFont(62);
    latex_par0->SetTextSize(0.05);
    latex_par0->Draw();
    latex_par1->SetTextFont(62);
    latex_par1->SetTextSize(0.05);
    latex_par1->Draw();
    latex_par2->SetTextFont(62);
    latex_par2->SetTextSize(0.05);
    latex_par2->Draw();
    latex_par4->SetTextFont(62);
    latex_par4->SetTextSize(0.05);
    latex_par4->Draw();

    latex_par3->SetTextFont(62);
    latex_par3->SetTextSize(0.05);
    if(corrections>0)latex_par3->Draw();

    //can->SaveAs(Form("pt_spectra/D2/plots/plot_%.1f_%.1f_%.1f_%.1f_%.1f_side_%i_mult_%i_%i.eps", pt1Min, pt1Max, pt2Min, ptAssocMin, ptAssocMax, side, multBinBegin, multBinEnd));

  if(par2==0.5||par2==0.07){
    if(error) *error = near->GetMaximum();
    return 0;
  }

  fit->Draw("same");
  if(error)
    *error = par1Err/binWidth;//TMath::Sqrt(2*TMath::Pi()*TMath::Power(par1Err*par2, 2) + TMath::Power(par1*par2Err, 2));
  Printf("par1 %f +/- %f", par1, par1Err);
 
  return result;
}


Int_t assignNotZeroError(TH1D* plot){
  
  Double_t err = 0;
  Int_t bins = 0;
  Int_t didCorrection = 0;

  for(int i=0; i<plot->GetNbinsX(); i++){
    if(plot->GetBinContent(i)>0 && (plot->GetBinCenter(i)>0.4||plot->GetBinCenter(i)<-0.4)){
      err+= plot->GetBinError(i);
      bins++;
    }
  }
  if(bins>0)
    err/=bins;
  else
    return;

  for(int i=0; i<plot->GetNbinsX(); i++){
    if(plot->GetBinCenter(i)<-1.45||plot->GetBinCenter(i)>1.45) continue;
    if(plot->GetBinContent(i)<=0){
      plot->SetBinError(i, err);
      didCorrection++;
    }
  }

  return didCorrection;
}


TGraphErrors* peakDifference_graph(const char* fileName, double pt1Min = 4.0, double pt1Max = 14.0, double pt2Min = 2.0, Int_t multBinBegin = 1, Int_t multBinEnd = 5, Int_t side = 0, Int_t yPos = 0, Int_t subtractMixedComb = 1, Int_t subtractFlow = 1, Int_t* trigger = NULL, Int_t* triggerBackground = NULL, double vtx = 9.9, Int_t method = 1){
  //TGraphErrors* peakDifference_graph(const char* fileName, double pt1Min = 4.0, double pt1Max = 14.0, double pt2Min = 2.0, Int_t multBinBegin = 1, Int_t multBinEnd = 5, Int_t side = 0, Int_t yPos = 0, Int_t subtractMixedComb = 1, Int_t subtractFlow = 1, Int_t* trigger = NULL, Int_t* triggerBackground = NULL){
  Double_t content[pt_assoc_bins_number-1];
  Double_t error[pt_assoc_bins_number-1];
  Int_t bins = 0;

  for(int i=0; i<pt_assoc_bins_number-1; i++){
    content[i] = peakYield(fileName, pt1Min, pt1Max, pt2Min, pt_assoc_bins[i], pt_assoc_bins[i+1], multBinBegin, multBinEnd, side, subtractMixedComb, subtractFlow, trigger, triggerBackground, vtx, &error[i], method);
    if(content>0) bins++;
  }

  //scale content with the bin width
  for(int i=0; i<pt_assoc_bins_number-1; i++){
    content[i] /= 2*pt_assoc_bin_error[i];
    error[i] /= 2*pt_assoc_bin_error[i];
  }

  TGraphErrors* graph = new TGraphErrors(bins, pt_assoc_bin_center, content, pt_assoc_bin_error, error);

  graph->GetXaxis()->SetTitle("p_{T,assoc} (GeV/c)");
  graph->GetYaxis()->SetTitle("1/N dN/dpT");
  graph->SetMarkerSize(1);
  graph->SetLineWidth(3);
  graph->SetMarkerColor(kBlue);
  graph->SetLineColor(kBlue);
  graph->SetMarkerStyle(20);

  return graph;
  
}

void compare2D(const char* fileName, const char* fileName2, double pt1Min = 4.0, double pt1Max = 14.0, Double_t pt2Min = 2.0, Double_t ptAssocMin = 0.5, Double_t ptAssocMax = 8.0, Double_t setVertex = 7, Int_t multBinBegin = 1, Int_t multBinEnd = 5, Int_t step = 0, Int_t subtractMixedComb = 0, Int_t subtractFlow =0){

  Int_t trigger1 = 0;
  Int_t background1 = 0;
  Int_t trigger2 = 0;
  Int_t background2 = 0;

  TCanvas* c1 = new TCanvas("file1", "file1", gBasisSize+50, gBasisSize+50, gBasisSize, gBasisSize);
  TH2D* file1 = getAnalysis_2D(fileName, pt1Min, pt1Max, pt2Min, ptAssocMin, ptAssocMax, setVertex, multBinBegin, multBinEnd, step, subtractMixedComb, subtractFlow, &trigger1, &background1);
  file1->DrawCopy("surf1");

  TCanvas* c2 = new TCanvas("file2", "file2", 2*gBasisSize+50, gBasisSize+50, gBasisSize, gBasisSize);
  TH2D* file2 = getAnalysis_2D(fileName2, pt1Min, pt1Max, pt2Min, ptAssocMin, ptAssocMax, setVertex, multBinBegin, multBinEnd, step, subtractMixedComb, subtractFlow, &trigger2, &background2);
  file2->DrawCopy("surf1");

  file1->Divide(file2);
  TCanvas* c3 = new TCanvas("mixed-ratio", "mixed-ratio", 3*gBasisSize+50, gBasisSize+50, gBasisSize, gBasisSize);
  file1->DrawCopy("surf1");
}


//uses different files
//have to hack the path so that two different paths can be loaded in the files if they are different
void comparePtSpectrum(const char* fileName, const char* fileName2, double pt1Min = 4.0, double pt1Max = 14.0, double pt2Min = 2.0, Int_t multBinBegin = 1, Int_t multBinEnd = 5, Int_t side = 0, Int_t yPos = 0, Int_t subtractMixedComb = 1, Int_t subtractFlow = 1){

  Int_t trigger_natural = 0;
  Int_t trigger_natural_background = 0;

  TGraphErrors* graph_cent_natural = (TGraphErrors*)(peakDifference_graph(fileName, pt1Min, pt1Max, pt2Min, multBinBegin, multBinEnd, side, yPos, subtractMixedComb, subtractFlow, &trigger_natural, &trigger_natural_background)->Clone());

  //TGraphErrors* graph_semi_natural = (TGraphErrors*)(peakDifference_graph(fileName, pt1Min, pt1Max, pt2Min, 10, 11, side, yPos, subtractMixedComb, subtractFlow, &trigger_natural, &trigger_natural_background)->Clone());

  TGraphErrors* graph_cent_efficiency = (TGraphErrors*)(peakDifference_graph(fileName2, pt1Min, pt1Max, pt2Min, multBinBegin, multBinEnd, side, yPos, subtractMixedComb, subtractFlow, &trigger_natural, &trigger_natural_background)->Clone());

  //TGraphErrors* graph_semi_efficiency = (TGraphErrors*)(peakDifference_graph(fileName2, pt1Min, pt1Max, pt2Min, 10, 11, side, yPos, subtractMixedComb, subtractFlow, &trigger_natural, &trigger_natural_background)->Clone());

  Double_t saveGraph_min = 0;
  Double_t saveGraph_max = 3;

  Double_t pt2Max = pt2Min+ptMaxAdd;
  
  save_graph_ratio(graph_cent_natural, graph_cent_efficiency, NULL, NULL, pt_assoc_bins_number, Form("no efficiency over efficiency pt spectrum for %.1f < p_{T,1} < %.1f and %.1f < p_{T,2} < %.1f", pt1Min, pt1Max, pt2Min, pt2Max), Form("pt_spectra/efficiency_ratio_%.0f_%.0f_%.0f.eps", pt1Min, pt1Max, pt2Min), "eff", "no eff", saveGraph_min, saveGraph_max, pt2Min);
  
  /*
  TGraphErrors* ratio_natural = (TGraphErrors*)(save_graph_ratio(graph_cent_natural, graph_semi_natural, NULL, NULL, pt_assoc_bins_number, Form("no efficiency central over semi central pt spectrum for %.1f < p_{T,1} < %.1f and %.1f < p_{T,2} < %.1f", pt1Min, pt1Max, pt2Min, pt2Max), Form("pt_spectra/no_efficiency_ratio_%.0f_%.0f_%.0f.eps", pt1Min, pt1Max, pt2Min), "eff", "no eff", saveGraph_min, saveGraph_max, pt2Min)->Clone());

  TGraphErrors* ratio_efficiency = (TGraphErrors*)(save_graph_ratio(graph_cent_efficiency, graph_semi_efficiency, NULL, NULL, pt_assoc_bins_number, Form("efficiency central over semi central pt spectrum for %.1f < p_{T,1} < %.1f and %.1f < p_{T,2} < %.1f", pt1Min, pt1Max, pt2Min, pt2Max), Form("pt_spectra/efficiency_ratio_%.0f_%.0f_%.0f.eps", pt1Min, pt1Max, pt2Min), "eff", "no eff", saveGraph_min, saveGraph_max, pt2Min)->Clone());

save_graph_ratio(ratio_natural, ratio_efficiency, NULL, NULL, pt_assoc_bins_number, Form("efficiency double ratio over efficiency pt spectrum for %.1f < p_{T,1} < %.1f and %.1f < p_{T,2} < %.1f", pt1Min, pt1Max, pt2Min, pt2Max), Form("pt_spectra/efficiency_double_ratio_%.0f_%.0f_%.0f.eps", pt1Min, pt1Max, pt2Min), "eff", "no eff", saveGraph_min, saveGraph_max, pt2Min);
  */

}



void showAll1plus1(const char* fileName, double pt1Min = 4.0, double pt1Max = 14.0, double pt2Min = 2.0, Int_t subtractMixedComb = 1, Int_t subtractFlow = 1){

  //show1plus1("AnalysisResults_1045.root", pt1Min, pt1Max, 1, 5, 0, subtractFlow);
  //show1plus1("AnalysisResults_1045.root", pt1Min, pt1Max, 9, 10, 1, subtractFlow);

  show1plus1(fileName, pt1Min, pt1Max, 1, 5, 0, subtractFlow);

  show1plus1(fileName, pt1Min, pt1Max, 9, 10, 1, subtractFlow);
}

void show1plus1(const char* fileName, double pt1Min = 4.0, double pt1Max = 14.0, Int_t multBinBegin = 1, Int_t multBinEnd = 5, Int_t yPos = 0, Int_t subtractFlow = 1){

  //for 1+1 there is no mixed comb
  Int_t subtractMixedComb = 0;

  const int pt1_bins = (pt1Max-pt1Min)/2;
  static const int pt1_bins_max = 20;

  if(pt1_bins>pt1_bins_max){
    printf("error not enough pt1_bins");
    return;
  }

  Double_t content[pt1_bins_max][pt_assoc_bins_number-1];
  Double_t error[pt1_bins_max][pt_assoc_bins_number-1];

  for(int i=0; i<pt_assoc_bins_number-1; i++){
    TCanvas* can = new TCanvas(Form("1plus1 %i, pT %i", yPos, i), Form("1plus1 %i, pT %i", yPos, i), i*gBasisSize+50, yPos*(gBasisSize+50), gBasisSize, gBasisSize);
     
    TLegend *leg  = getLegend();

    for(int j=0; j<pt1_bins; j++){
      double pt1Minimum = pt1Min+j*2.0;
      TH1D* onePlusOne = (TH1*)getAnalysis(fileName, pt1Minimum, pt1Minimum + 2.0, 0, pt_assoc_bins[i], pt_assoc_bins[i+1], 8.9, multBinBegin, multBinEnd, 6, 2, 1, 0, subtractMixedComb, subtractFlow);
      
      Double_t phi_bin = getPhiBinsForAnalysis(pt_assoc_bins[i]);

      Int_t bin_start = onePlusOne->FindBin(-1*phi_bin);
      Int_t bin_end = onePlusOne->FindBin(phi_bin);

      content[j][i] = onePlusOne->IntegralAndError(bin_start, bin_end, error[j][i]);

      if(j==1)
	onePlusOne->SetLineColor(kRed);
      else if(j==2)
	onePlusOne->SetLineColor(kGreen);
      else if(j==3)
	onePlusOne->SetLineColor(kMagenta);

      if(j==0){
	((TH1*)(onePlusOne->Clone()))->DrawCopy();
      }else{
	((TH1*)(onePlusOne->Clone()))->DrawCopy("same");
      }

      leg->AddEntry((TH1*)(onePlusOne->Clone()), Form("pT %.0f", pt1Minimum) ,"l");
      if(j==pt1_bins-1)
	leg->Draw("same");
    }
  }


  TCanvas* can_graph = new TCanvas(Form("result %i", yPos), Form("result %i", yPos), pt_assoc_bins_number*gBasisSize+50, yPos*(gBasisSize+50), gBasisSize, gBasisSize);

  TLegend *leg  = getLegend();

  for(int j=0; j<pt1_bins; j++){
    TGraphErrors* graph = new TGraphErrors(pt_assoc_bins_number, pt_assoc_bin_center, content[j], pt_assoc_bin_error, error[j]);

    graph->GetXaxis()->SetTitle("p_{T,assoc} (GeV/c)");
    graph->GetYaxis()->SetTitle("1/N dN/dpT");
    graph->SetMarkerSize(2);
    graph->SetLineWidth(2);
    graph->SetMarkerStyle(20);

    if(j==0){
      //graph->SetMarkerColor(kCyan);
      //graph->SetLineColor(kCyan);
      graph->SetMarkerColor(kBlue);
      graph->SetLineColor(kBlue);
      (TGraphErrors*)(graph->Clone())->Draw("AP");
      //(TGraphErrors*)(graph->Clone())->Draw("P");
    }else if(j==1){
      //graph->SetMarkerColor(kOrange);
      //graph->SetLineColor(kOrange);
      graph->SetMarkerColor(kRed);
      graph->SetLineColor(kRed);
      (TGraphErrors*)(graph->Clone())->Draw("P");
    }else if(j==2){
      //graph->SetMarkerColor(kMagenta);
      //graph->SetLineColor(kMagenta);
      graph->SetMarkerColor(kGreen);
      graph->SetLineColor(kGreen);
      (TGraphErrors*)(graph->Clone())->Draw("P");
    }

    double pt1Minimum = pt1Min+j*2.0;
    leg->AddEntry((TGraphErrors*)(graph->Clone()),Form("pT %.0f", pt1Minimum),"l");
    if(j==pt1_bins-1)
      leg->Draw("same");
  }
  

}


//does the reading out of the results
//divides the same event by the mixed events
//subtracts the mixed combinatorics
TH2D* getAnalysis_2D(const char* fileName, double pt1Min = 4.0, double pt1Max = 14.0, Double_t pt2Min = 2.0, Double_t ptAssocMin = 0.5, Double_t ptAssocMax = 8.0, Double_t setVertex = 7, Int_t multBinBegin = 1, Int_t multBinEnd = 5, Int_t step = 0, Int_t subtractMixedComb = 0, Int_t subtractFlow =0, Int_t* trigger = NULL, Int_t* triggerBackgroundSame = NULL)
{
  loadlibs();
  //to guarantee to pick only the bins from the choosen pt on
  pt1Min += 0.01;
  pt1Max -= 0.01;
  pt2Min += 0.01;
  ptAssocMin += 0.01;
  ptAssocMax -= 0.01;

  Double_t pt2Max = pt2Min+ptMaxAdd-0.02;
  
  AliTwoPlusOneContainer* h = (AliTwoPlusOneContainer*) GetTwoPlusOne(fileName);

  AliUEHist::CFStep step_same = step;//(AliUEHist::CFStep) AliTwoPlusOneContainer::kSameNS;
  AliUEHist::CFStep step_mixed = (step%2)+2;//(AliUEHist::CFStep) AliTwoPlusOneContainer::kMixedNS;
  if(step>=7) step_mixed = ((step+1)%2)+2;
  AliUEHist::CFStep step_mixedComb = (step%2)+4;//(AliUEHist::CStep) AliTwoPlusOneContainer::kMixedCombNS;
  AliUEHist::CFStep step_backgroundSame = (step%2)+7;//(AliUEHist::CStep) AliTwoPlusOneContainer::kBackgroundSameNS;

  AliUEHist::CFStep step_1plus1 = (AliUEHist::CFStep) AliTwoPlusOneContainer::k1plus1;
  AliUEHist::CFStep step_1plus1_mixed = (AliUEHist::CFStep) AliTwoPlusOneContainer::kMixed1plus1;

  h->GetData()->SetPtRange(ptAssocMin, ptAssocMax);
  h->GetData()->SetPt2Min(pt2Min);
  h->GetData()->SetPt2Max(pt2Max);
  h->GetData()->SetZVtxRange(-1*setVertex, setVertex);
  //  h->GetData()->SetZVtxRange(8.1, 9.9);//run this with getAnalysis("AR_1502.root", 4, 8, 4, 1, 2, 7.9, 1, 6, 6, 1, 1, 2, 0, 0, NULL, NULL, 0); it shows that the last vertex bin produces a ridge on Delta eta=0;
  
  //GetSumOfRatios2(mixed, step, region, ptLeadMin, ptLeadMax, multBinBegin, multBinEnd, normalizePerTrigger, stepForMixed)
  //TH2D* h2_etaPhi;
  TH2D* h2_etaPhi_mixedComb;
  TH2D* h2_etaPhi_backgroundSame;

  //TH2D* h2_etaPhi2;
  //TH2D* h2_etaPhi_mixedComb2;

  if(step>=2||!subtractMixedComb){
    Int_t trigger_etaPhi;
    if(step!=6)
      h2_etaPhi = (TH2D*) h->GetData()->GetSumOfRatios2(h->GetData(), step_same, 0, pt1Min, pt1Max, multBinBegin, multBinEnd, kTRUE, step_mixed, &trigger_etaPhi);
    else  if(step==6){
      h->GetData()->SetPt2Max(-1);
      h->GetData()->SetPt2Min(-1);
      h2_etaPhi = (TH2D*) h->GetData()->GetSumOfRatios2(h->GetData(), step_1plus1, 0, pt1Min, pt1Max, multBinBegin, multBinEnd, kTRUE, step_1plus1_mixed, &trigger_etaPhi);
      h->GetData()->SetPt2Max(pt2Max);
      h->GetData()->SetPt2Min(pt2Min);
    }

    Printf("########################################");
    Printf("trigger 1+1: %i", trigger_etaPhi);
    Printf("########################################");

    if(trigger!=NULL)
      *trigger = (Int_t)trigger_etaPhi;
  }else if(step<2 && subtractMixedComb){
  
    Int_t trigger_same;
    Int_t trigger_mixed_comb;
    Int_t trigger_background;
    Int_t badSMratio = 0;
	  
    h2_etaPhi = (TH2D*) h->GetData()->GetSumOfRatios2(h->GetData(), step_same, 0, pt1Min, pt1Max, multBinBegin, multBinEnd, kFALSE, step_mixed, &trigger_same);
    //h2_etaPhi = (TH2D*) h->GetData()->GetSumOfRatios2(h->GetData(), step_mixedComb, 0, pt1Min, pt1Max, multBinBegin, multBinEnd, kTRUE, step_mixed, &trigger_same);
    //h2_etaPhi->Scale(trigger_same);

    /*
    AliTwoPlusOneContainer* h2;
    if(fileName=="AR_1899.root")
      h2 = (AliTwoPlusOneContainer*) GetTwoPlusOne("AR_1916.root");
    else
      h2 = (AliTwoPlusOneContainer*) GetTwoPlusOne(fileName);
    h2->GetData()->SetPtRange(ptAssocMin, ptAssocMax);
    h2->GetData()->SetPt2Min(pt2Min);
    h2->GetData()->SetPt2Max(pt2Max);
    h2->GetData()->SetZVtxRange(-1*setVertex, setVertex);
    h2->GetData()->SetPt2Max(-1);//1+1 spectrum has no real information about pt2
    h2->GetData()->SetPt2Min(-1);
    */
    int bin_x = h2_etaPhi->GetXaxis()->FindBin(0.01);
    int bin_y = h2_etaPhi->GetYaxis()->FindBin(0.01);
    int global_bin = h2_etaPhi->GetBin(bin_x, bin_y);
    Printf("##################################################################################");
    Printf("loot at same bin %f +/- %f", h2_etaPhi->GetBinContent(global_bin), h2_etaPhi->GetBinError(global_bin));

    //don't need to use getMixedComb_scaled_backgroundSame, see compareScaledMixedComb which shows that both methods are identical but here getSumOfRatios is used which makes it easier (and safer against errors)
    //the 1+1 analysis can only be used if the full analysis was done within the same pt bins
    h->GetData()->SetPt2Max(-1);//1+1 spectrum has no real information about pt2
    h->GetData()->SetPt2Min(-1);
    if(step==0){
      Printf("USE normal method, step %i", step);
      h2_etaPhi_mixedComb = (TH2D*) h->GetData()->GetSumOfRatios2(h->GetData(), step_1plus1, 0, pt1Min, pt1Max, multBinBegin, multBinEnd, kFALSE, step_1plus1_mixed, &trigger_mixed_comb);
    }else if(step==1){
      Printf("USE normal METHOD, step %i", step);
      h2_etaPhi_mixedComb = (TH2D*) h->GetData()->GetSumOfRatios2(h->GetData(), step_1plus1, 0, pt2Min, pt2Max, multBinBegin, multBinEnd, kFALSE, step_1plus1_mixed, &trigger_mixed_comb);
    }else{ //this is never used, it's the old method of getting the background
      Printf("#########################");
      Printf("USE OLD BACKGROUND METHOD, step %i", step);
      Printf("#########################");
      h2_etaPhi_mixedComb = (TH2D*) h->GetData()->GetSumOfRatios2(h->GetData(), step_mixedComb, 0, pt1Min, pt1Max, multBinBegin, multBinEnd, kFALSE, step_mixed, &trigger_mixed_comb);
    }
    Printf("##################################################################################");
    Printf("loot at 1+1 bin %f +/- %f", h2_etaPhi_mixedComb->GetBinContent(global_bin), h2_etaPhi_mixedComb->GetBinError(global_bin));

    /*
    h->GetData()->SetPtRange(pt2Min, pt2Max);
    Int_t trig2_1plus1_trigger = 0;
    TH2D* trig2_1plus1 = (TH2D*) h->GetData()->GetSumOfRatios2(h->GetData(), step_1plus1, 0, pt1Min, pt1Max, multBinBegin, multBinEnd, kFALSE, step_1plus1_mixed, &trig2_1plus1_trigger);
    TH1D* h1_trig2_1plus1 = trig2_1plus1->ProjectionX("px_trigger estmation", -1, -1);

    Int_t guessed_trigger = h1_trig2_1plus1->GetMinimum();
    guessed_trigger *= TMath::Pi()/4;//scale by the alpha aceptance
    h->GetData()->SetPtRange(ptAssocMin, ptAssocMax);

    TCanvas* can = new TCanvas("can", "can", gBasisSize+100, 50, gBasisSize, gBasisSize);
    h1_trig2_1plus1->DrawCopy();
    */

    h->GetData()->SetPt2Max(pt2Max);
    h->GetData()->SetPt2Min(pt2Min);
    //h2_etaPhi_mixedComb->Scale(trigger_mixed_comb);

    h2_etaPhi_backgroundSame = (TH2D*) h->GetData()->GetSumOfRatios2(h->GetData(), step_backgroundSame, 0, pt1Min, pt1Max, multBinBegin, multBinEnd, kFALSE, step_mixed, &trigger_background);
    //h->GetData()->GetSumOfRatios2(h->GetData(), step_mixedComb, 0, pt1Min, pt1Max, multBinBegin, multBinEnd, kFALSE, step_mixed, &trigger_background);

    //trigger_background -= (trigger_same-trigger_background);//assumption that the background is too big

    //Printf("trigger same: %i, trigger 1plus1 minimum: %i ", trigger_same, guessed_trigger);
    //trigger_background_same *= 2;//error in the plots for runs 968 and 973
    //trigger_background = trigger_same*0.9;
    if(false){
      double trigger_ratio = (double)trigger_background;
      if(trigger_mixed_comb>0)
	trigger_ratio /= trigger_mixed_comb;
      h2_etaPhi_mixedComb->Scale(trigger_ratio);
      Printf("trigger 1+1: %i", trigger_mixed_comb);
      Printf("########################################");
      Printf("trigger same: %i, trigger background %i, trigger ratio ", trigger_same, trigger_background, trigger_ratio);
      Printf("########################################");     
      
      h2_etaPhi->Add(h2_etaPhi_mixedComb, -1);
      //h2_etaPhi->Add(h2_etaPhi_backgroundSame, -1);
    
      h2_etaPhi_mixedComb->Scale(1.0/(double)trigger_background);

      Printf("##################################################################################");
      Printf("loot at same bin after the correction %f +/- %f", h2_etaPhi->GetBinContent(global_bin), h2_etaPhi->GetBinError(global_bin));
      
      if((trigger_same-trigger_background)>0)
	h2_etaPhi->Scale(1.0/(trigger_same-trigger_background));

    }else{
      Printf("########################################");
      Printf("trigger same: %i, trigger background %i, trigger 1+1 %i", trigger_same, trigger_background, trigger_mixed_comb);
      Printf("########################################");     
 

      //correct handling of all the errors
      subtract_uncorrBG(h2_etaPhi, h2_etaPhi_mixedComb, trigger_same, trigger_mixed_comb, trigger_background);
      //subtract_uncorrBG(h2_etaPhi, h2_etaPhi_backgroundSame, trigger_same, trigger_mixed_comb, trigger_background);
    }
      Printf("loot at same bin after the correction and scaling %f +/- %f", h2_etaPhi->GetBinContent(global_bin), h2_etaPhi->GetBinError(global_bin));


    if(trigger!=NULL)
      *trigger = trigger_same-trigger_background;
    if(triggerBackgroundSame!=NULL)
      *triggerBackgroundSame = trigger_background;



  }else{
    Printf("cannot subtract mixed combinatorics from step %i ", step);
    return 0;
  }

  //if(subtractFlow)
  //subtractFlow(h2_etaPhi, 1.0, 1.4);

  if(h2_etaPhi==NULL)
    return NULL;
  
  h2_etaPhi->GetYaxis()->SetRangeUser(-1.6, 1.6);//corresponds to the really used eta range
 
  //undo normalization
  Float_t normalization = h2_etaPhi->GetXaxis()->GetBinWidth(1);
  h2_etaPhi->Scale(normalization);
  //h2_etaPhi->Rebin2D(2,2);

  //do eta normalization
  //Float_t normalization_eta = h2_etaPhi->GetYaxis()->GetBinWidth(1);
  //h2_etaPhi->Scale(1/normalization_eta);

  return h2_etaPhi;

}



void* subtract_uncorrBG(TH2D* h2_sameEvent, TH2D* h2_bgEvent, Int_t sameEvent_trigger, Int_t onePlusOne_trigger, Int_t bg_trigger){

  TH2D* h2_sameEvent2 = h2_sameEvent->Clone("_clone1");
  Int_t realCorr_trigger = sameEvent_trigger-bg_trigger;

  double trigger_ratio = (double)bg_trigger;
  if(onePlusOne_trigger>0)
    trigger_ratio /= onePlusOne_trigger;

  //comment this out for BG same systematic error
  h2_bgEvent->Scale(trigger_ratio);

  h2_sameEvent->Add(h2_bgEvent, -1);
  if((sameEvent_trigger-bg_trigger)>0)
    h2_sameEvent->Scale(1.0/realCorr_trigger);

  Double_t err_mult = TMath::Sqrt((double)sameEvent_trigger)/realCorr_trigger;

  TH2D* h2_bgEvent2 = h2_bgEvent->Clone("_clone2");
  h2_bgEvent2->Scale((double)sameEvent_trigger/bg_trigger);
  h2_sameEvent2->Add(h2_bgEvent2, -1);
  Double_t err_scaling = TMath::Sqrt(bg_trigger)/TMath::Power(realCorr_trigger, 2);
  h2_sameEvent2->Scale(err_scaling);

  TCanvas* can = new TCanvas();
  h2_bgEvent->DrawCopy("surf1");


    int bin_x = h2_sameEvent->GetXaxis()->FindBin(0.01);
    int bin_y = h2_sameEvent->GetYaxis()->FindBin(0.01);
    int global_bin = h2_sameEvent->GetBin(bin_x, bin_y);
    Printf("##################################################################################");
    Printf("look at same bin %f +/- %f", h2_sameEvent->GetBinContent(global_bin), h2_sameEvent->GetBinError(global_bin));

  for(int i=0; i<h2_sameEvent->GetNbinsX(); i++){
    for(int j=0; j<h2_sameEvent->GetNbinsX(); j++){
      Double_t cont = h2_sameEvent->GetBinContent(i, j);
      Double_t err = h2_sameEvent->GetBinError(i, j);
      Double_t err_sameTrigger = err_mult*cont;
      Double_t err_bgTrigger = h2_sameEvent2->GetBinContent(i, j);
      Double_t final_err = TMath::Sqrt(TMath::Power(err, 2) + TMath::Power(err_sameTrigger, 2) + TMath::Power(err_bgTrigger, 2));
      h2_sameEvent->SetBinError(i, j, final_err);

      if(i==bin_x && j==bin_y){
	Printf("normal stat error: %f, error same trigger: %f, error bgTrigger %f", err, err_sameTrigger, err_bgTrigger);
      }

    }
  }
}






TH1D* getFlowDist(const char* fileName, double pt1Min = 4.0, double pt1Max = 14.0, Double_t pt2Min = 2.0, Double_t ptAssocMin = 0.5, Double_t ptAssocMax = 8.0, Double_t setVertex = 7, Int_t multBinBegin = 1, Int_t multBinEnd = 5, Int_t step = 0, Int_t subtractMixedComb = 0, Double_t flow_start, Double_t flow_end){

  TH2D* h2_etaPhi = getAnalysis_2D(fileName, pt1Min, pt1Max, pt2Min, ptAssocMin, ptAssocMax, setVertex, multBinBegin, multBinEnd, step, subtractMixedComb, 0, NULL, NULL);

  return subtractFlow(h2_etaPhi, flow_start, flow_end);

}


Double_t getConstBkg(const char* fileName, double pt1Min = 4.0, double pt1Max = 14.0, Double_t pt2Min = 2.0, Double_t ptAssocMin = 0.5, Double_t ptAssocMax = 8.0, Double_t setVertex = 7, Int_t multBinBegin = 1, Int_t multBinEnd = 5, Int_t step = 0, Int_t subtractMixedComb = 0, Double_t bkg_start, Double_t bkg_end){

  TH2D* h2_etaPhi = getAnalysis_2D(fileName, pt1Min, pt1Max, pt2Min, ptAssocMin, ptAssocMax, setVertex, multBinBegin, multBinEnd, step, subtractMixedComb, 0, NULL, NULL);

  return subtractConstBkg(h2_etaPhi, bkg_start, bkg_end, ptAssocMin);

}

//this method is more or less identical to the subtractFlow. Just the statistical error is smaller
void compareBkgSubtractions(const char* fileName, double pt1Min = 4.0, double pt1Max = 14.0, Double_t pt2Min = 2.0, Double_t ptAssocMin = 0.5, Double_t ptAssocMax = 8.0, Double_t setVertex = 7, Int_t multBinBegin = 1, Int_t multBinEnd = 5, Int_t step = 0, Int_t subtractMixedComb = 0, Double_t bkg_start, Double_t bkg_end){

  TH1D* flow = getFlowDist(fileName, pt1Min, pt1Max, pt2Min, ptAssocMin, ptAssocMax, setVertex, multBinBegin, multBinEnd, step, subtractMixedComb, bkg_start, bkg_end);

  Double_t const_bkg = getConstBkg(fileName, pt1Min, pt1Max, pt2Min, ptAssocMin, ptAssocMax, setVertex, multBinBegin, multBinEnd, step, subtractMixedComb, bkg_start, bkg_end);

  TCanvas* can = new TCanvas();
  flow->DrawCopy();

  Printf("constant Bkg %f", const_bkg);

  Double_t phi_bin = getPhiBinsForAnalysis(ptAssocMin);
  Int_t start_bin = flow->FindBin(-1.0*phi_bin);
  Int_t end_bin = flow->FindBin(phi_bin);

  Double_t avg = 0;
  for(int i = start_bin; i<=end_bin; i++){
    Printf("bin %i, content: %f", i, flow->GetBinContent(i));
    avg += flow->GetBinContent(i);
  }

  avg /= end_bin-start_bin+1;

  Printf("Average: %f", avg);

}


Int_t getNumberOfTrigger(const char* fileName, Double_t pt1Min, Double_t pt1Max, Double_t pt2Min, Double_t ptAssocMin, Double_t ptAssocMax, Double_t setVertex, Int_t multBinBegin, Int_t multBinEnd)
{

  AliTwoPlusOneContainer* twoPlusOne = (AliTwoPlusOneContainer*) GetTwoPlusOne(fileName);

  twoPlusOne->GetData()->SetPtRange(ptAssocMin, ptAssocMax);
  twoPlusOne->GetData()->SetPt2Min(pt2Min);
  twoPlusOne->GetData()->SetZVtxRange(-1*setVertex, setVertex);

  AliUEHist::CFStep step_same = (AliUEHist::CFStep) AliTwoPlusOneContainer::kSameAS;
  AliUEHist::CFStep step_mixed = (AliUEHist::CFStep) AliTwoPlusOneContainer::kMixedAS;
  AliUEHist::CFStep step_backgroundSame = (AliUEHist::CFStep) AliTwoPlusOneContainer::kBackgroundSameAS;

  Int_t trigger_same = 0;
  Int_t trigger_background_same = 0;

  twoPlusOne->GetData()->GetSumOfRatios2(twoPlusOne->GetData(), step_same, 0, pt1Min, pt1Max, multBinBegin, multBinEnd, kFALSE, step_mixed, &trigger_same);
  twoPlusOne->GetData()->GetSumOfRatios2(twoPlusOne->GetData(), step_backgroundSame, 0, pt1Min, pt1Max, multBinBegin, multBinEnd, kFALSE, step_mixed, &trigger_background_same);

  return trigger_same-trigger_background_same;

}



//does the reading out of the results
//divides the same event by the mixed events
//subtracts the mixed combinatorics
TH1D* getAnalysis(const char* fileName, double pt1Min = 4.0, double pt1Max = 14.0, Double_t pt2Min = 2.0, Double_t ptAssocMin = 0.5, Double_t ptAssocMax = 8.0, Double_t setVertex = 7, Int_t multBinBegin = 1, Int_t multBinEnd = 5, Int_t step = 0, Int_t posX = 1, Int_t posY = 1, Int_t showPlots = 0, Int_t subtractMixedComb = 0, Int_t subtractFlow =0, Int_t* trigger = NULL, Int_t* triggerBackground = NULL, Int_t plotAxis = 1)
{
 
  TH2D* h2_etaPhi = getAnalysis_2D(fileName, pt1Min, pt1Max, pt2Min, ptAssocMin, ptAssocMax, setVertex, multBinBegin, multBinEnd, step, subtractMixedComb, subtractFlow, trigger, triggerBackground)->Clone("h2_etaphi_clone");
  //TH2D* h2_etaPhi2 = getAnalysis_2D(fileName, pt1Min, pt1Max, pt2Min, ptAssocMin, ptAssocMax, -1*setVertex, multBinBegin, multBinEnd, step, subtractMixedComb, subtractFlow, trigger, triggerBackground)->Clone("h2_etaphi_clone");
  //h2_etaPhi->Add(h2_etaPhi2, 1.0);

  if(h2_etaPhi==NULL)
    return NULL;

  //removeWing(h2_etaPhi);

  Int_t trigger2;
  Int_t triggerBackground2;
  TH2D* h1_etaPhi;
  /*
  if(step==0)
    h1_etaPhi = getAnalysis_2D(fileName, pt1Min, pt1Max, -1, ptAssocMin, ptAssocMax, setVertex, multBinBegin, multBinEnd, 6, subtractMixedComb, subtractFlow, &trigger2, &triggerBackground2);
  else if(pt1Min>pt2Min)
    h1_etaPhi = getAnalysis_2D(fileName, pt2Min, pt1Min, -1, ptAssocMin, ptAssocMax, setVertex, multBinBegin, multBinEnd, 6, subtractMixedComb, subtractFlow, &trigger2, &triggerBackground2);
  else
    h1_etaPhi = getAnalysis_2D(fileName, pt2Min, pt1Max, -1, ptAssocMin, ptAssocMax, setVertex, multBinBegin, multBinEnd, 6, subtractMixedComb, subtractFlow, &trigger2, &triggerBackground2);
  */
  //h1_etaPhi = (TH2D*)getAnalysis_2D(fileName, pt1Min, pt1Max, -1, ptAssocMin, ptAssocMax, setVertex, multBinBegin, multBinEnd, 6, subtractMixedComb, subtractFlow, &trigger2, &triggerBackground2)->Clone("h1_etaphi_clone");

  if(subtractFlow==1 && plotAxis !=2){
    //subtractFlow(h2_etaPhi, h1_etaPhi, 1.0, 1.4);
    subtractFlow(h2_etaPhi, 1.0, 1.4);
  }else if(subtractFlow==2){
    subtractConstBkg(h2_etaPhi, 1.0, 1.4, ptAssocMin);//make sure subtractBaseline is deactivated
   }

  Double_t eta_bin = getEtaBinsForAnalysis(ptAssocMin);
  Printf("eta_bin first %f", eta_bin);

  int firstbin = h2_etaPhi->GetYaxis()->FindBin(-1*eta_bin);
  int lastbin = h2_etaPhi->GetYaxis()->FindBin(eta_bin);
  //Printf("got bins from %i to %i", firstbin, lastbin);
  //Printf("bin center : %f and %f", h2_etaPhi->GetYaxis()->GetBinCenter(firstbin), h2_etaPhi->GetYaxis()->GetBinCenter(lastbin));

  TH1D* h1_phi = projectToTH1D(h2_etaPhi, "fully corrected", firstbin, lastbin, 1.0);
  //TH1D* h1_phi = projectToTH1D(h2_etaPhi_mixedComb, "fully corrected", firstbin, lastbin, 1.0);//shows background distributions
  
  Double_t phi_bin = 0.51;//getPhiBinsForAnalysis(ptAssocMin);//0.51;

  int firstbin_phi = h2_etaPhi->GetXaxis()->FindBin(-1*phi_bin);
  int lastbin_phi = h2_etaPhi->GetXaxis()->FindBin(phi_bin);
  TH1D* h1_eta = projectToTH1D_eta(h2_etaPhi, "fully corrected", firstbin_phi, lastbin_phi, 1.0);
  //Printf("bin center found at x axis %f and %f", h2_etaPhi->GetXaxis()->GetBinCenter(firstbin_phi), h2_etaPhi->GetXaxis()->GetBinCenter(lastbin_phi));
  //Printf("bin center h1_eta: %f and %f", h1_eta->GetBinCenter(firstbin), h1_eta->GetBinCenter(lastbin));
  


  if(showPlots>1){ 
    TCanvas* can_4 = new TCanvas(Form("dim2_fully_corrected_%i_%i", posX, posY), Form("dim2_fully_corrected_%i_%i", posX, posY), posX*gBasisSize+100, 50, gBasisSize, gBasisSize);
    //TCanvas* can_4 = new TCanvas("2dim_fully_corrected_1_1", "2dim_fully_corrected_1_1", posX*gBasisSize+100, 50, gBasisSize, gBasisSize);
    //h2_etaPhi->GetYaxis()->SetRangeUser(-1.6, 1.6);
    h2_etaPhi->SetTitle("");
    h2_etaPhi->DrawCopy("surf1");
  }
  if(showPlots>2){

    TCanvas* can_4b = new TCanvas(Form("1dim phi fully corrected, %i, %i ", posX, posY), Form("1dim phi fully corrected, %i, %i ", posX, posY), posX*gBasisSize+100, gBasisSize+50, gBasisSize, gBasisSize);
    //h2_etaPhi_mixedComb->DrawCopy("surf1");
    h1_phi->SetTitle("");
    h1_phi->DrawCopy();
   
    Double_t error;
    Double_t content = h1_phi->IntegralAndError(firstbin_phi, lastbin_phi, error);
    //Printf("bin center h1_phi: %f and %f", h1_phi->GetBinCenter(firstbin_phi), h1_phi->GetBinCenter(lastbin_phi));
    //Printf("bin content h1_phi: %f and %f", h1_phi->GetBinContent(firstbin_phi), h1_phi->GetBinContent(lastbin_phi));
    Printf("found Integral phi: %f +/- %f", content, error);


    TCanvas* can_4c = new TCanvas(Form("1dim eta fully corrected, %i, %i ", posX, posY), Form("1dim eta fully corrected, %i, %i ", posX, posY), (posX+1)*gBasisSize+100, gBasisSize+50, gBasisSize, gBasisSize);
    h1_eta->SetTitle("");
    h1_eta->DrawCopy();

    Double_t error_eta;
    Double_t content_eta = h1_eta->IntegralAndError(firstbin-1, lastbin-1, error_eta);//the eta bins are moved compared to the 2d version before
    //Printf("bin center h1_eta: %f and %f", h1_eta->GetBinCenter(firstbin-1), h1_eta->GetBinCenter(lastbin-1));
    //Printf("bin content h1_eta: %f and %f", h1_eta->GetBinContent(firstbin), h1_eta->GetBinContent(lastbin));
    Printf("found Integral eta: %f +/- %f", content_eta, error_eta);

    /*
    TH1D* test_phi = h2_etaPhi->ProjectionX("abc px", firstbin, lastbin);
    Double_t test_phi_err;
    Printf("phi bins: %i from %f to %f ", test_phi->GetNbinsX(), test_phi->GetBinCenter(bin_start), test_phi->GetBinCenter(bin_end));
    Double_t test_phi_cont = test_phi->IntegralAndError(bin_start, bin_end, test_phi_err);
    Printf("found phi Integral: %f +/- %f", test_phi_cont, test_phi_err);

    Printf("2dim bins per axis: %i, %i ", h2_etaPhi->GetNbinsX(), h2_etaPhi->GetNbinsY());
    TH1D* test_eta = h2_etaPhi->ProjectionY("abc py", bin_start, bin_end);
    Double_t test_eta_err;
    Printf("eta bins: %i from %f to %f ", test_eta->GetNbinsX(), test_eta->GetBinCenter(firstbin), test_eta->GetBinCenter(lastbin));
    Double_t test_eta_cont = test_eta->IntegralAndError(firstbin, lastbin, test_eta_err);
    Printf("found eta Integral: %f +/- %f", test_eta_cont, test_eta_err);
    */
  }

  if(showPlots>5){//==2
    TCanvas* can_5 = new TCanvas(Form("compare projection, %i, %i ", posX, posY), Form("compare projection, %i, %i ", posX, posY), posX*gBasisSize+100, posY*gBasisSize+50, gBasisSize, gBasisSize);
    TLegend *leg2  = getLegend();

    h1_phi->DrawCopy();
    leg2->AddEntry(h1_phi,"chosen correction","l");
    TH1D* h1_phi_raw = projectToTH1D(h2_etaPhi_clone, "raw same event", firstbin, lastbin, 1.0);
    h1_phi_raw->SetLineColor(kRed);
    h1_phi_raw->DrawCopy("same");
    leg2->AddEntry(h1_phi_raw,"raw same event","l");
    TH1D* h1_phi_mixedComb = projectToTH1D(h2_etaPhi_mixedComb, "mixed comb", firstbin, lastbin, 1.0);
    h1_phi_mixedComb->SetLineColor(kGreen);
    h1_phi_mixedComb->DrawCopy("same");
    leg2->AddEntry(h1_phi_mixedComb,"mixed comb","l");
    leg2->Draw("same");
  }

  Printf("pt1min %f, pt1max %f, pt2min %f, ptassocmin %f,  ptassocmax %f ", pt1Min, pt1Max, pt2Min, ptAssocMin, ptAssocMax);
  /*
    TF1* fit = new TF1("fit", "[0]+[1]*exp(-1.0*[2]*x**2)", -1.0*gSubtractBaseline_end, gSubtractBaseline_end);
    //fit with gaussian distribution
    h1_phi->Fit("fit", "N", "", -1.0*gSubtractBaseline_end, gSubtractBaseline_end);
    Double_t par1 = fit->GetParameter(1);
    Double_t par1Err = fit->GetParError(1);

    Double_t par2 = fit->GetParameter(2);
    Double_t par2Err = fit->GetParError(2);

    //content[i] = par1*TMath::Sqrt(TMath::Pi()/par2);
    */

  if(plotAxis==1)
    return h1_phi;
  else
    return h1_eta;
}













void show_steps_of_getAnalysis(const char* fileName, double pt1Min = 4.0, double pt1Max = 14.0, Double_t pt2Min = 2.0, Double_t ptAssocMin = 0.5, Double_t ptAssocMax = 8.0, Double_t setVertex = 7, Int_t multBinBegin = 1, Int_t multBinEnd = 5, Int_t step = 0, Int_t subtractMixedComb = 0, Int_t subtractFlow =0)
{
  if(step>1){
    Printf("choose step between 0 and 1!");
  }

  TCanvas* can1 = new TCanvas("sameEvent near side", "sameEvent near side", 100, 50, gBasisSize, gBasisSize);
  Int_t trigger = 0;
  TH2D* h2_etaPhi = (TH2D*)(getAnalysis_2D(fileName, pt1Min, pt1Max, pt2Min, ptAssocMin, ptAssocMax, setVertex, multBinBegin, multBinEnd, step, 0, subtractFlow, &trigger))->Clone();
  //h2_etaPhi->DrawCopy("surf1");
  h2_etaPhi->DrawCopy("colz");

  Double_t eta_bin = getEtaBinsForAnalysis(ptAssocMin);

  int firstbin = h2_etaPhi->GetYaxis()->FindBin(-1*eta_bin);
  int lastbin = h2_etaPhi->GetYaxis()->FindBin(eta_bin);

  TH1D* h1_phi = (projectToTH1D(h2_etaPhi, "same near", firstbin, lastbin, trigger))->Clone();
  h1_phi->SetLineColor(kRed);
  TCanvas* can_1 = new TCanvas("sameEvent near side projection", "sameEvent near side projection", 100, gBasisSize+50, gBasisSize, gBasisSize);
  h1_phi->DrawCopy();

  TCanvas* can2 = new TCanvas("sameEvent mixed comb", "sameEvent mixed comb", gBasisSize+100, 50, gBasisSize, gBasisSize);
  Int_t trigger_mixed_comb = 0;
  TH2D* h2_etaPhi_mixedComb = (TH2D*)(getAnalysis_2D(fileName, pt1Min, pt1Max, pt2Min, ptAssocMin, ptAssocMax, setVertex, multBinBegin, multBinEnd, step+4, 0, subtractFlow, &trigger_mixed_comb))->Clone();

  Printf("#####################################################################");
  Printf("trigger same: %i, trigger mixed comb: %i", trigger, trigger_mixed_comb);
  Printf("#####################################################################");

  h2_etaPhi_mixedComb->DrawCopy("surf1");
  TH1D* h1_phi_mixedComb = (projectToTH1D(h2_etaPhi_mixedComb, "mixedComb near", firstbin, lastbin, trigger_mixed_comb))->Clone();
  Printf(Form("trigger mixed_comb", trigger_mixed_comb));
  TCanvas* can_2 = new TCanvas("sameEvent mixed comb projection", "sameEvent mixed comb projection", gBasisSize+100, gBasisSize+50, gBasisSize, gBasisSize);
  h1_phi_mixedComb->DrawCopy();

  TCanvas* can2b = new TCanvas("1+1", "1+1", 2*gBasisSize+100, 50, gBasisSize, gBasisSize);
  Int_t trigger_1plus1 = 0;
  TH2D* h2_etaPhi_1plus1 = getAnalysis_2D(fileName, pt1Min, pt1Max, pt2Min, ptAssocMin, ptAssocMax, setVertex, multBinBegin, multBinEnd, 6, 1, subtractFlow, &trigger_1plus1);
  h2_etaPhi_1plus1->DrawCopy("surf1");
  TH1D* h1_phi_1plus1 = (projectToTH1D(h2_etaPhi_1plus1, "1plus1 phi", firstbin, lastbin, trigger_1plus1))->Clone();
  TCanvas* can_2b = new TCanvas("1+1 projection", "1+1 projection", 2*gBasisSize+100, gBasisSize+50, gBasisSize, gBasisSize);
  h1_phi_1plus1->DrawCopy();



  TCanvas* can3 = new TCanvas("mixed comb subtracted", "mixed comb subtracted", 3*gBasisSize+100, 50, gBasisSize, gBasisSize);
  Int_t trigger2 = 0;
  TH2D* h2_etaPhi_subtracted = getAnalysis_2D(fileName, pt1Min, pt1Max, pt2Min, ptAssocMin, ptAssocMax, setVertex, multBinBegin, multBinEnd, step, 1, subtractFlow, &trigger2);
  h2_etaPhi_subtracted->DrawCopy("surf1");
  TH1D* h1_phi_subtracted = (projectToTH1D(h2_etaPhi_subtracted, "mixedComb near", firstbin, lastbin, 1.0))->Clone();
  TCanvas* can_3 = new TCanvas("mixed comb subtracted projection", "mixed comb subtracted projection", 3*gBasisSize+100, gBasisSize+50, gBasisSize, gBasisSize);
  h1_phi_subtracted->DrawCopy();

  TCanvas* can_3_2 = new TCanvas("mixed comb subtracted 2 projection", "mixed comb subtracted 2 projection", 3*gBasisSize+100, 2*gBasisSize+50, gBasisSize, gBasisSize);
  h1_phi->Add(h1_phi_1plus1, -1);
  h1_phi->DrawCopy();

  h2_etaPhi->Add(h2_etaPhi_mixedComb, -1);
  if(trigger-trigger_mixed_comb>0)
    h2_etaPhi->Scale(1.0/(trigger-trigger_mixed_comb));
  TCanvas* can4 = new TCanvas("mixed comb subtracted reproduced", "mixed comb subtracted reproduced", 4*gBasisSize+100, 50, gBasisSize, gBasisSize);
  h2_etaPhi->DrawCopy("surf1");
  TH1D* h1_phi_reproduced = (projectToTH1D(h2_etaPhi, "same near reproduced", firstbin, lastbin, 1.0))->Clone();
  TCanvas* can_4 = new TCanvas("mixed comb subtracted reproduced projection", "mixed comb subtracted reproduced projection", 4*gBasisSize+100, gBasisSize+50, gBasisSize, gBasisSize);
  h1_phi_reproduced->DrawCopy();


  Printf("trigger same: %i, trigger mixed comb: %i", trigger, trigger_mixed_comb);

}


void show_steps_of_getAnalysis_projection(const char* fileName, double pt1Min = 4.0, double pt1Max = 14.0, Double_t pt2Min = 2.0, Double_t ptAssocMin = 0.5, Double_t ptAssocMax = 8.0, Double_t setVertex = 7, Int_t multBinBegin = 1, Int_t multBinEnd = 5, Int_t step = 0, Int_t subtractMixedComb = 0, Int_t subtractFlow =0)
{
  if(step>1){
    Printf("choose step between 0 and 1!");
  }

  Int_t trigger = 0;
  TH2D* h2_etaPhi = (TH2D*)(getAnalysis_2D(fileName, pt1Min, pt1Max, pt2Min, ptAssocMin, ptAssocMax, setVertex, multBinBegin, multBinEnd, step, 0, subtractFlow, &trigger))->Clone();

  Double_t eta_bin = getEtaBinsForAnalysis(ptAssocMin);

  int firstbin = h2_etaPhi->GetYaxis()->FindBin(-1*eta_bin);
  int lastbin = h2_etaPhi->GetYaxis()->FindBin(eta_bin);

  TH1D* h1_phi = (projectToTH1D(h2_etaPhi, "same near", firstbin, lastbin, trigger))->Clone();
  h1_phi->SetLineColor(kRed);

  Int_t trigger_mixed_comb = 0;
  TH2D* h2_etaPhi_mixedComb = (TH2D*)(getAnalysis_2D(fileName, pt1Min, pt1Max, pt2Min, ptAssocMin, ptAssocMax, setVertex, multBinBegin, multBinEnd, step+4, 0, subtractFlow, &trigger_mixed_comb))->Clone();
  TH1D* h1_phi_mixedComb = (projectToTH1D(h2_etaPhi_mixedComb, "mixedComb near", firstbin, lastbin, trigger_mixed_comb))->Clone();
  h1_phi_mixedComb->SetLineColor(kOrange);

  Printf("#####################################################################");
  Printf("trigger same: %i, trigger mixed comb: %i", trigger, trigger_mixed_comb);
  Printf("#####################################################################");

  Int_t trigger_1plus1 = 0;
  TH2D* h2_etaPhi_1plus1;
  if(step==0)
    h2_etaPhi_1plus1 = getAnalysis_2D(fileName, pt1Min, pt1Max, pt2Min, ptAssocMin, ptAssocMax, setVertex, multBinBegin, multBinEnd, 6, 1, subtractFlow, &trigger_1plus1);
  else
    h2_etaPhi_1plus1 = getAnalysis_2D(fileName, pt2Min, pt1Min, pt2Min, ptAssocMin, ptAssocMax, setVertex, multBinBegin, multBinEnd, 6, 1, subtractFlow, &trigger_1plus1); 
  TH1D* h1_phi_1plus1 = (projectToTH1D(h2_etaPhi_1plus1, "1plus1 phi", firstbin, lastbin, trigger_1plus1))->Clone();
  h1_phi_1plus1->SetLineColor(kGreen);
  //h1_phi_1plus1->SetLineWidth(4.0);

  Int_t trigger2 = 0;
  TH2D* h2_etaPhi_subtracted = getAnalysis_2D(fileName, pt1Min, pt1Max, pt2Min, ptAssocMin, ptAssocMax, setVertex, multBinBegin, multBinEnd, step, 1, subtractFlow, &trigger2);
  TH1D* h1_phi_subtracted = (projectToTH1D(h2_etaPhi_subtracted, "mixedComb near", firstbin, lastbin, 1.0))->Clone();
  h1_phi_subtracted->SetLineColor(kBlue);

  h2_etaPhi->Add(h2_etaPhi_1plus1, -1.0);
  TH1D* h1_phi_2 = (projectToTH1D(h2_etaPhi, "same near", firstbin, lastbin, 1.0))->Clone();
  h1_phi_2->SetLineColor(kCyan);

  TCanvas* can = new TCanvas("can", "can", gBasisSize+100, gBasisSize+50, 2*gBasisSize, 2*gBasisSize);
  h1_phi->DrawCopy();
  //h1_phi_mixedComb->DrawCopy("same");
  h1_phi_1plus1->DrawCopy("same");
  h1_phi_subtracted->DrawCopy("same");
  //h1_phi_2->DrawCopy("same");

  TLegend *leg  = getLegend();
  leg->SetTextSize(gStyle->GetTextSize()*0.8);
  leg->AddEntry(h1_phi,"same event raw","l");
  leg->AddEntry(h1_phi_1plus1,"1+1 correlations","l");
  leg->AddEntry(h1_phi_subtracted,"same event corrected","l");

  leg->Draw("same");

  Printf("trigger same: %i, trigger mixed comb: %i", trigger, trigger_mixed_comb);

}










TH1D* projectToTH1D_eta(TH2D* etaPhi, char* name, int firstbin, int lastbin, int trigger){
  //  Printf(Form("name %s", name));

  TH1D* h1_phi_1 = etaPhi->ProjectionY(Form("px_%s", name), firstbin, lastbin, "e");
  h1_phi_1->SetLineWidth(3);
  h1_phi_1->SetStats(kFALSE);
  h1_phi_1->Scale(1.0/trigger);
  //h1_phi_1->Scale(1.0/(double)(lastbin-firstbin+1));
  //h1_phi_1->Scale(etaPhi->GetXaxis()->GetBinWidth(firstbin));//scale with phi bin width 
  //h1_phi_1->GetXaxis()->SetRangeUser(-TMath::Pi()/2+0.01, TMath::Pi()/2-0.01);
  //h1_phi_1->GetXaxis()->SetRangeUser(-TMath::Pi()/2+0.01, 3*TMath::Pi()/2-0.01);
  //h1_phi_1->SetYTitle("1/N \\ dN/(d \\Delta \\eta)");
  h1_phi_1->GetYaxis()->SetTitleOffset(1.2);
  h1_phi_1->GetYaxis()->SetTitleSize(0.05);
  h1_phi_1->GetXaxis()->SetTitleSize(0.05);
  h1_phi_1->GetYaxis()->SetLabelSize(0.05);
  h1_phi_1->GetXaxis()->SetLabelSize(0.05);
  h1_phi_1->SetYTitle("\\frac{1}{N_{trig}} \\frac{d^2N_{assoc}}{d\\Delta\\varphi d\\Delta\\eta}");

  //symmetrize(h1_phi_1);
  //subtractBaseline(h1_phi_1);

  return h1_phi_1;
}


TH1D* projectToTH1D(TH2D* etaPhi, char* name, int firstbin, int lastbin, int trigger){
  //  Printf(Form("name %s", name));

  //firstbin = -1;
  //lastbin = -1;

  TH1D* h1_phi_1 = etaPhi->ProjectionX(Form("py_%s", name), firstbin, lastbin, "e");
  h1_phi_1->SetLineWidth(3);
  h1_phi_1->SetStats(kFALSE);
  h1_phi_1->Scale(1.0/trigger);
  //h1_phi_1->Scale(1.0/(double)(lastbin-firstbin+1));
  //h1_phi_1->Scale(etaPhi->GetYaxis()->GetBinWidth(firstbin));//scale with eta bin width 
  //h1_phi_1->GetXaxis()->SetRangeUser(-TMath::Pi()/2+0.01, TMath::Pi()/2-0.01);
  //h1_phi_1->GetXaxis()->SetRangeUser(-TMath::Pi()/2+0.01, 3*TMath::Pi()/2-0.01);
  h1_phi_1->SetYTitle("1/N \\ dN/(d \\Delta \\varphi)");

  //symmetrize(h1_phi_1);
  subtractBaseline(h1_phi_1);

  return h1_phi_1;
}

TLegend* getLegend(double x_start = 0.65, double y_start = 0.7, double x_end = 0.85, double y_end = 0.9){
  TLegend *leg  = new TLegend(x_start, y_start, x_end, y_end);
  leg->SetFillColor(10);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextFont(gStyle->GetTextFont());
  leg->SetTextSize(gStyle->GetTextSize()*1.0);
  return leg;
}

void compareAllBackground(const char* fileName, double pt1Min = 4.0, double pt1Max = 14.0, Double_t pt2Min = 2.0, Int_t subtractMixedComb = 0, Int_t subtractFlow = 0){

  //static const int pt_assoc_bins_number = 8;
  //Double_t pt_assoc_bins[pt_assoc_bins_number] = {0.5, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0};

  for(int i=0; i<pt_assoc_bins_number-1; i++){
    compareBackground(fileName, pt1Min, pt1Max, pt2Min, pt_assoc_bins[i], pt_assoc_bins[i+1], 1, 5,  i, 0, subtractMixedComb, subtractFlow);
    compareBackground(fileName, pt1Min, pt1Max, pt2Min, pt_assoc_bins[i], pt_assoc_bins[i+1], 9, 10,  i, 1, subtractMixedComb, subtractFlow);
  }

}


//compares the background of the mixed comb analysis with the 1+1 background
//statistics for background same is the worst because in each event at 2 places is searched for T2
//for mixed comb a lot of events are searched for the second trigger particle
//for 1+1 no second trigger particle needs to be found at all so the statistics is the highest
void compareBackground(const char* fileName, double pt1Min = 4.0, double pt1Max = 14.0, Double_t pt2Min = 2.0, Double_t ptAssocMin = 0.5, Double_t ptAssocMax = 8.0, Int_t multBinBegin = 1, Int_t multBinEnd = 5, Int_t posX = 1, Int_t posY = 1, Int_t subtractMixedComb = 0, Int_t subtractFlow = 0){

  Int_t trigger_mixed_comb=0;
  Int_t trigger_1plus1=0;
  Int_t trigger_background_same=0;

  TH1D* mixedComb = (TH1*)getAnalysis(fileName, pt1Min, pt1Max, pt2Min, ptAssocMin, ptAssocMax, 8.9, multBinBegin, multBinEnd, 0, 1, 1, 0, subtractMixedComb, subtractFlow, &trigger_mixed_comb)->Clone();

  TH1D* backgroundSame = (TH1*)getAnalysis(fileName, pt1Min, pt1Max, pt2Min, ptAssocMin, ptAssocMax, 8.9, multBinBegin, multBinEnd, 8, 1, 1, 0, subtractMixedComb, subtractFlow, &trigger_background_same)->Clone();

  //near side
  //TH1D* onePlusOne = (TH1*)getAnalysis(fileName, pt1Min, pt1Max, pt2Min, ptAssocMin, ptAssocMax, 8.9, multBinBegin, multBinEnd, 6, 1, 1, 0, subtractMixedComb, subtractFlow, &trigger_1plus1)->Clone();
  //away side
  TH1D* onePlusOne = (TH1*)getAnalysis(fileName, pt2Min, pt1Min, pt2Min, ptAssocMin, ptAssocMax, 8.9, multBinBegin, multBinEnd, 6, 1, 1, 0, subtractMixedComb, subtractFlow, &trigger_1plus1)->Clone();

  onePlusOne->SetLineColor(kRed);
  //backgroundSame->SetLineColor(kGreen);
  backgroundSame->SetLineColor(kBlue);

  Printf("found trigger: mixed comb %i, 1plus1 %i, background_same %i", trigger_mixed_comb, trigger_1plus1, trigger_background_same);
  //do not need to scale these distributions because both are already divided by the number of triggers

  TCanvas* can = new TCanvas(Form("compare background, %i, %i ", posX, posY), Form("compare background, %i, %i ", posX, posY), posX*gBasisSize+100, posY*gBasisSize+50, gBasisSize, gBasisSize);

  backgroundSame->SetLineWidth(3);
  mixedComb->SetLineWidth(3);
  onePlusOne->SetLineWidth(3);

  ((TH1*)(backgroundSame->Clone()))->DrawCopy();
  //((TH1*)(mixedComb->Clone()))->DrawCopy("same");
  ((TH1*)(onePlusOne->Clone()))->DrawCopy("same");

  TLegend *leg  = getLegend();
  //leg->AddEntry(mixedComb, "mixed combinatorics","l");
  //leg->AddEntry(onePlusOne, "1plus1 background","l");
  //leg->AddEntry(backgroundSame, "background same","l");

  //leg->AddEntry(mixedComb, "2+1 near side","l");
  leg->AddEntry(onePlusOne, "1+1 scaled","l");
  leg->AddEntry(backgroundSame, "background  at pi/2","l");
  leg->Draw("same");
}


void compareDeltaEtaGap(const char* fileName, double pt1Min = 4.0, double pt1Max = 14.0, Double_t pt2Min = 2.0, Double_t ptAssocMin = 0.5, Double_t ptAssocMax = 8.0, Int_t multBinBegin = 1, Int_t multBinEnd = 5, Int_t step = 0, Int_t subtractMixedComb = 0){
  /*
  TH1D* first = (TH1*)getFlowDist(fileName, pt1Min, pt1Max, pt2Min, ptAssocMin, ptAssocMax, 8.9, multBinBegin, multBinEnd, step, subtractMixedComb, 1.0, 1.4)->Clone();

  TH1D* second = (TH1*)getFlowDist(fileName, pt1Min, pt1Max, pt2Min, ptAssocMin, ptAssocMax, 8.9, multBinBegin, multBinEnd, step, subtractMixedComb, 1.2, 1.6)->Clone();

  TH1D* third = (TH1*)getFlowDist(fileName, pt1Min, pt1Max, pt2Min, ptAssocMin, ptAssocMax, 8.9, multBinBegin, multBinEnd, step, subtractMixedComb, 1.0, 1.6)->Clone();
  */
  /*
  TH1D* first = (TH1*)getFlowDist(fileName, pt1Min, pt1Max, pt2Min, ptAssocMin, ptAssocMax, 7.9, multBinBegin, multBinEnd, 0, subtractMixedComb, 1.0, 1.4)->Clone();

  TH1D* second = (TH1*)getFlowDist(fileName, pt1Min, pt1Max, pt2Min, ptAssocMin, ptAssocMax, 7.9, multBinBegin, multBinEnd, 7, subtractMixedComb, 1.0, 1.4)->Clone();

  TH1D* third = (TH1*)getFlowDist(fileName, pt1Min, pt1Max, pt2Min, ptAssocMin, ptAssocMax, 7.9, multBinBegin, multBinEnd, 6, subtractMixedComb, 1.0, 1.4)->Clone();
  */
  TH1D* first = (TH1*)getFlowDist(fileName, pt1Min, pt1Max, pt2Min, ptAssocMin, ptAssocMax, 7.9, 1, 6, 0, subtractMixedComb, 1.0, 1.4)->Clone();
  TH1D* second = (TH1*)getFlowDist(fileName, pt1Min, pt1Max, pt2Min, ptAssocMin, ptAssocMax, 7.9, 10, 11, 0, subtractMixedComb, 1.0, 1.4)->Clone();

  Double_t subtract =  first->GetBinContent(first->GetMinimumBin());
  Double_t subtract_err = first->GetBinError(first->GetMinimumBin());

    for(int i=0; i<=first->GetNbinsX(); i++){
      first->SetBinContent(i, first->GetBinContent(i) - subtract);
      Double_t bin_error = TMath::Sqrt(TMath::Power(subtract_err, 2) + TMath::Power(first->GetBinError(i), 2));
      first->SetBinError(i, bin_error);
    }
  subtract =  second->GetBinContent(second->GetMinimumBin());
  subtract_err = second->GetBinError(second->GetMinimumBin());

    for(int i=0; i<=second->GetNbinsX(); i++){
      second->SetBinContent(i, second->GetBinContent(i) - subtract);
      Double_t bin_error = TMath::Sqrt(TMath::Power(subtract_err, 2) + TMath::Power(second->GetBinError(i), 2));
      second->SetBinError(i, bin_error);
    }



  second->SetLineColor(kRed);
  //third->SetLineColor(kGreen);

  TCanvas* can = new TCanvas("compare Delta eta", "compare Delta eta", gBasisSize+100, gBasisSize+50, gBasisSize, gBasisSize);

  first->SetLineWidth(5);
  second->SetLineWidth(4);
  //third->SetLineWidth(5);

  first->DrawCopy();
  second->DrawCopy("same");
  //third->DrawCopy("same");

  TLegend *leg  = getLegend();
  /*
  leg->AddEntry(first, "1.0 - 1.6","l");
  leg->AddEntry(second, "1.0 - 1.4","l");
  leg->AddEntry(third, "1.2 - 1.6","l");
  */
  //leg->AddEntry(first, "2+1","l");
  //leg->AddEntry(second, "background same","l");
  //leg->AddEntry(third, "1+1","l");

leg->AddEntry(first, "central","l");
leg->AddEntry(second, "semi central","l");

  leg->Draw("same");
}




void compareAllMixedComb(const char* fileName, double pt1Min = 4.0, double pt1Max = 14.0, Double_t pt2Min = 2.0, Int_t subtractMixedComb = 0, Int_t subtractFlow = 0){

  //static const int pt_assoc_bins_number = 8;
  //Double_t pt_assoc_bins[pt_assoc_bins_number] = {0.5, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0};

  for(int i=0; i<pt_assoc_bins_number-1; i++){
    compareMixedComb_sides(fileName, pt1Min, pt1Max, pt2Min, pt_assoc_bins[i], pt_assoc_bins[i+1], 1, 5,  i, 0, subtractMixedComb, subtractFlow);
    compareMixedComb_sides(fileName, pt1Min, pt1Max, pt2Min, pt_assoc_bins[i], pt_assoc_bins[i+1], 9, 10,  i, 1, subtractMixedComb, subtractFlow);
  }

}


void compareMixedComb_sides(const char* fileName, double pt1Min = 4.0, double pt1Max = 14.0, Double_t pt2Min = 2.0, Double_t ptAssocMin = 0.5, Double_t ptAssocMax = 8.0, Int_t multBinBegin = 1, Int_t multBinEnd = 5, Int_t posX = 1, Int_t posY = 1, Int_t subtractMixedComb = 0, Int_t subtractFlow = 0){

  Int_t trigger_mixed_comb=0;

  TH1D* mixedComb_near = (TH1*)getAnalysis(fileName, pt1Min, pt1Max, pt2Min, ptAssocMin, ptAssocMax, 6.9, multBinBegin, multBinEnd, 4, 1, 1, 0, subtractMixedComb, subtractFlow, &trigger_mixed_comb)->Clone();

  TH1D* mixedComb_away = (TH1*)getAnalysis(fileName, pt1Min, pt1Max, pt2Min, ptAssocMin, ptAssocMax, 6.9, multBinBegin, multBinEnd, 5, 1, 1, 0, subtractMixedComb, subtractFlow, &trigger_mixed_comb)->Clone();

 TCanvas* can = new TCanvas(Form("compare mixedComb, %i, %i ", posX, posY), Form("compare mixedComb, %i, %i ", posX, posY), posX*gBasisSize+100, posY*gBasisSize+50, gBasisSize, gBasisSize);

 mixedComb_away->SetLineColor(kRed);

  ((TH1*)(mixedComb_near->Clone()))->DrawCopy();
  ((TH1*)(mixedComb_away->Clone()))->DrawCopy("same");

  TLegend *leg  = getLegend();
  leg->AddEntry(mixedComb_near, "near side","l");
  leg->AddEntry(mixedComb_away, "away side","l");
  //leg->Draw("same");
}




//shows alpha dependent significance of the signal, signal/background for different centralities
void get1plus1(const char* fileName, Double_t pt1Min, Double_t pt1Max, Double_t ptAssocMin, Double_t ptAssocMax){
  loadlibs();
  TFile::Open(fileName);

  //to guarantee to pick only the bins from the choosen pt on
  pt1Min += 0.01;
  pt1Max -= 0.01;
  ptAssocMin += 0.01;
  ptAssocMax -= 0.01;

  Int_t alpha_bins = 5;//number of alpha bins maximal used around the two bins close to pi, maximum usable would be 5 because than already the full away side is searched for a trigger 2
  Int_t showAlpha = 1;
  Double_t alpha_bin_width = 0.174532;
  Int_t draw = 1;

  Int_t multBinBegin = 9;
  Int_t multBinEnd = 10;

  AliTwoPlusOneContainer* h = (AliTwoPlusOneContainer*) GetTwoPlusOne(fileName);

  AliUEHist::CFStep step_1plus1 = (AliUEHist::CFStep) AliTwoPlusOneContainer::k1plus1;
  AliUEHist::CFStep step_mixed = (AliUEHist::CFStep) AliTwoPlusOneContainer::kMixedNS;
  
  TH1D* alpha_plot_significance = new TH1D("alpha_plot", "", alpha_bins, 0, alpha_bins*alpha_bin_width);
  alpha_plot_significance->SetXTitle("\\alpha");
  alpha_plot_significance->SetYTitle("significance");

  TH1D* alpha_plot_sigBack = new TH1D("alpha_plot", "", alpha_bins, 0, alpha_bins*alpha_bin_width);
  alpha_plot_sigBack->SetXTitle("\\alpha");
  alpha_plot_sigBack->SetYTitle("signal/background");

  for(Int_t cent = 0; cent<4; cent++){
    char* cent_string = "";
    if(cent==0){
      multBinBegin = 1;
      multBinEnd = 5;
      cent_string = "0-5%";
    }else if(cent==1){
      multBinBegin = 6;
      multBinEnd = 6;
      cent_string = "5-10%";
    }else if(cent==2){
      multBinBegin = 7;
      multBinEnd = 8;
      cent_string = "10-30%";
    }else if(cent==3){
      multBinBegin = 9;
      multBinEnd = 10;
      cent_string = "30-50%";
    }

    for(Int_t alpha = 0; alpha <alpha_bins; alpha++){
      if(alpha==showAlpha) draw = 1; else draw = 0;
      
      h->GetData()->SetPtRange(ptAssocMin, ptAssocMax);
      TH2* h2_etaPhi = h->GetData()->GetSumOfRatios2(h->GetData(), step_1plus1, 0, pt1Min, pt1Max, multBinBegin, multBinEnd, kFALSE, step_mixed);
      Double_t arr1[4];
      showResult(h2_etaPhi, cent, 0, alpha, arr1, draw);
      alpha_plot_significance->SetBinContent(alpha+1, arr1[0]);
      alpha_plot_significance->SetBinError(alpha+1, arr1[2]);
      alpha_plot_sigBack->SetBinContent(alpha+1, arr1[1]);
      alpha_plot_sigBack->SetBinError(alpha+1, arr1[3]);
    }

    TCanvas* c1 = new TCanvas(Form("alpha_plot_significance %s ", cent_string), Form("alpha_plot_significance %s ", cent_string), cent*gBasisSize+50, 1*gBasisSize+50, gBasisSize, gBasisSize);
    alpha_plot_significance->DrawCopy("E");

    TCanvas* c2 = new TCanvas(Form("alpha_plot_signalBackground %s ", cent_string), Form("alpha_plot_signalBackground %s ", cent_string), cent*gBasisSize+50, 2*gBasisSize+50, gBasisSize, gBasisSize);
    alpha_plot_sigBack->DrawCopy("E");
  }

}



//method used just in get1plus1 to calculate the significance and the signal/background
void showResult(TH2* h2_etaPhi, Int_t posX, Int_t posY, Int_t alpha_bins, Double_t* arr, Int_t draw){

  Printf("#################################################");
  Printf("g_eta_bin is used instead of getEtaBinsForAnalysis!!! ");//change the finding below to getEtaBinsForAnalysis
  Printf("#################################################");

  int firstbin = h2_etaPhi->GetYaxis()->FindBin(-1*g_eta_bin);
  int lastbin = h2_etaPhi->GetYaxis()->FindBin(g_eta_bin);

  Printf("firstbin %i, lastbin %i, \n", firstbin, lastbin);

  TH1D* h1_phi = h2_etaPhi->ProjectionX(Form("px_%i_%i", posX, posY), firstbin, lastbin, "e");
  
  double alpha = 0.15;// use 0.15 to get only the next two bins which are a bit smaller than 0.2//0.2;
  //Pi is exactely the bin edge
  int bin_1 = 27 - alpha_bins;//h1_phi->FindBin(TMath::Pi()-alpha+0.01);
  int bin_2 = 28 + alpha_bins;//h1_phi->FindBin(TMath::Pi()+alpha-0.01);

  Int_t bin_start = h1_phi->FindBin(TMath::Pi()/3+0.01);
  Int_t bin_end = h1_phi->FindBin(2*TMath::Pi()/3-0.01);
  Double_t par0 = 0;
  Double_t err0 = 0;

  for(int i=bin_start; i<=bin_end; i++){
    par0 += h1_phi->GetBinContent(i);
    err0 += h1_phi->GetBinError(i);
  }
  par0 /= bin_end - bin_start + 1;  
  err0 /= bin_end - bin_start + 1;


  double sum = 0;
  double sum_err = 0;
  for(int bin = bin_1; bin<=bin_2; bin++){
    sum += h1_phi->GetBinContent(bin);
    sum_err += h1_phi->GetBinError(bin);
  }
  sum_err /= (bin_2-bin_1+1);
  double background = (bin_2 - bin_1 + 1)*par0;
  double background_err = (bin_2 - bin_1 + 1)*err0;
  double signal2 = sum - background;
  double signal2_err = TMath::Sqrt(sum_err*sum_err+background_err*background_err);

  Double_t significance = signal2/h1_phi->GetBinError((bin_1+bin_2)/2);
  Double_t significance_err = signal2_err/(TMath::Sqrt(sum));
  //Double_t significance = signal2/(TMath::Sqrt(sum));
  //Double_t significance_err = TMath::Sqrt(TMath::Power(signal2_err/(TMath::Sqrt(sum)), 2) + TMath::Power(signal2/(TMath::Sqrt(sum*sum*sum)*sum_err), 2));

  Double_t significance_rel = signal2/background;//signal over background
  Double_t significance_rel_err = TMath::Sqrt(TMath::Power(signal2_err/background,2)+TMath::Power(signal2/(background*background)*background_err,2));//signal over background

  Printf("X %i, Y %i, significance: %f, signal over background %f ", posX, posY, significance, significance_rel);

  h1_phi->SetYTitle("1/N_{trig} \\  dN_{assoc}/d \\varphi");
  h1_phi->SetTitle("");
  h1_phi->SetStats(kFALSE);
  //h1_phi->Scale(1/(lastbin - firstbin + 1));
  if(draw>0){
    TCanvas* c1 = new TCanvas(Form("can %i %i", posX, posY), Form("can %i %i", posX, posY), posX*gBasisSize+50, posY*gBasisSize+50, gBasisSize, gBasisSize);
    //h1_phi->DrawCopy();
    h2_etaPhi->DrawCopy("surf1");
  }

  arr[0] = significance;
  arr[1] = significance_rel;
  arr[2] = significance_err;
  arr[3] = significance_rel_err;
}


void symmetrize(TH1D* tracks)
{
  int bins = tracks->FindBin(TMath::Pi()/2-0.01);
  for(int i=1; i<bins/2; i++){
    double tmp =  (tracks->GetBinContent(i)+tracks->GetBinContent(bins-i+1))/2;
    tracks->SetBinContent(i, tmp);
    tracks->SetBinContent(bins-i+1, tmp);    
  }
}


void compareMixedEvent(const char* fileName, double pt1Min = 4.0, double pt1Max = 14.0, double pt2Min = 2.0, Double_t ptAssocMin = 0.5, Double_t ptAssocMax = 8.0, Int_t multBinBegin = 1, Int_t multBinEnd = 5, Int_t vertexBinBegin = 1, Int_t vertexBinEnd = 5, Int_t step = 1, Int_t symmetrize = 0){

  pt1Min += 0.01;
  pt1Max -= 0.01;
  pt2Min += 0.01;
  ptAssocMin += 0.01;
  ptAssocMax -= 0.01;

  TH2D* twoPlusOne = getMixedEvent(fileName, pt1Min, pt1Max ,pt2Min, ptAssocMin, ptAssocMax, multBinBegin, multBinEnd, vertexBinBegin, vertexBinEnd, 3, symmetrize)->Clone("_clone");

  //TH2D* onePlusOne = getMixedEvent(fileName, pt1Min, pt1Max ,pt2Min, ptAssocMin, ptAssocMax, multBinBegin, multBinEnd, vertexBinBegin, vertexBinEnd, 9, symmetrize)->Clone("_clone2");
  TH2D* onePlusOne = getMixedEvent(fileName, pt2Min, pt2Min+ptMaxAdd-0.02 ,pt2Min, ptAssocMin, ptAssocMax, multBinBegin, multBinEnd, vertexBinBegin, vertexBinEnd, 9, symmetrize)->Clone("_clone2");

 TCanvas* c1a = new TCanvas("mixedEvent1a", "mixedEvent1a", 50, 1*gBasisSize+50, gBasisSize, gBasisSize);
 twoPlusOne->DrawCopy("surf1");

 TCanvas* c1b = new TCanvas("mixedEvent1b", "mixedEvent1b", 50, 2*gBasisSize+50, gBasisSize, gBasisSize);
 onePlusOne->DrawCopy("surf1");


  twoPlusOne->Divide(onePlusOne);

  TCanvas* c2 = new TCanvas("mixedEvent", "mixedEvent", gBasisSize+50, 2*gBasisSize+50, gBasisSize, gBasisSize);
  twoPlusOne->DrawCopy("surf1");

}




TH2D* getMixedEvent(const char* fileName, double pt1Min = 4.0, double pt1Max = 14.0, double pt2Min = 2.0, Double_t ptAssocMin = 0.5, Double_t ptAssocMax = 8.0, Int_t multBinBegin = 1, Int_t multBinEnd = 5, Int_t vertexBinBegin = 1, Int_t vertexBinEnd = 5, Int_t step = 1, Int_t symmetrize = 0)
{

  TH2D* tracksMixed = getCorrelationStep(fileName, pt1Min, pt1Max, pt2Min, ptAssocMin, ptAssocMax, multBinBegin, multBinEnd, vertexBinBegin, vertexBinEnd, step);

  double b1 = tracksMixed->GetBinContent(tracksMixed->FindBin(0.001,0.001));
  double b2 = tracksMixed->GetBinContent(tracksMixed->FindBin(0.001,-0.001));
  double b3 = tracksMixed->GetBinContent(tracksMixed->FindBin(-0.001,0.001));
  double b4 = tracksMixed->GetBinContent(tracksMixed->FindBin(-0.001,-0.001));
  double divide = (b1 + b2 + b3 + b4)/4;//=value at zero
  double error = TMath::Power(b1/divide-1, 2)+TMath::Power(b2/divide-1, 2)+TMath::Power(b3/divide-1, 2)+TMath::Power(b4/divide-1, 2);//not further used
  Printf(Form(" filename: %s, four bins: %f, %f, %f, %f, divide: %f", fileName, b1/divide, b2/divide, b3/divide, b4/divide, divide));

  /*  
  for(int i=0; i<10000000; i++){
    Double_t tmp = tracksMixed->GetBinContent(tracksMixed->GetNbinsX()/2, tracksMixed->GetNbinsY()/2);
    tracksMixed->SetBinContent(tracksMixed->GetNbinsX()/2, tracksMixed->GetNbinsY()/2, tmp+1);
  }
  */

  //symmetrize
  if(symmetrize>0)
    for(i=0; i<=; i++){
      int y_bins = tracksMixed->GetNbinsY();
      for(j=0; j<=y_bins/2; j++){
	double tmp =  (tracksMixed->GetBinContent(i, j)+tracksMixed->GetBinContent(i, y_bins-j+1))/2;
	tracksMixed->SetBinContent(i, j, tmp);
	tracksMixed->SetBinContent(i, y_bins-j+1, tmp);
      }
    }
  
  //finite bin correction
  //ftrackEtaCut = 0.9; 
  Double_t finiteBinCorrection = -1.0 / (2*0.9) * tracksMixed->GetYaxis()->GetBinWidth(1) / 2 + 1;
  Printf("Finite bin correction: %f", finiteBinCorrection);
  divide /= finiteBinCorrection;
  error /= finiteBinCorrection;//not further used

  //tracksMixed->Scale(1.0/divide);

  TH1D* proj = tracksMixed->ProjectionY("_norm");//, "e"
  Int_t bins_X = tracksMixed->GetNbinsX();

  double norm1 = proj->GetBinContent(proj->FindBin(-0.001));
  double norm2 = proj->GetBinContent(proj->FindBin(0.001));
  double norm = (norm1+norm2)/2.0;

  //tracksMixed->Scale(1.0/(bins_X*norm));
  Printf("norm1 %f, norm2 %f, norm %f", norm1/norm, norm2/norm, norm);

  TCanvas* c2 = new TCanvas("mixedEvent", "mixedEvent", gBasisSize+50, 2*gBasisSize+50, gBasisSize, gBasisSize);
  tracksMixed->DrawCopy("surf1");

  return tracksMixed;
}

void getParticleDist(const char* fileName, Double_t ptAssocMin = 0.5, Double_t ptAssocMax = 8.0, Int_t multBinBegin = 1, Int_t multBinEnd = 5, Int_t vertexBinBegin = 1, Int_t vertexBinEnd = 5, Int_t step = 1){

  TH2D* particleDist = getCorrelationStep(fileName, 1, 20, 4, ptAssocMin, ptAssocMax, multBinBegin, multBinEnd, vertexBinBegin, vertexBinEnd, step);

  TCanvas* can = new TCanvas();
  particleDist->DrawCopy("surf1");
}



//old method with a double cache
//used to compare results but very difficult to read
const char* lastFileName_corStep = 0;
double lastPt1Min = 0;
double lastPt1Max = 0;
double lastPt2Min = 0;
double lastPtAssocMin = 0;
double lastPtAssocMax = 0;
int laststep = 0;
void* cacheHistsZVTxMult = 0;
void* cacheEvtZVTxMult = 0;

const char* lastFileName2_corStep = 0;
double lastPt1Min2 = 0;
double lastPt1Max2 = 0;
double lastPt2Min2 = 0;
double lastPtAssocMin2 = 0;
double lastPtAssocMax2 = 0;
int laststep2 = 0;
void* cacheHistsZVTxMult2 = 0;
void* cacheEvtZVTxMult2 = 0;

int next = 1;

TH2D* getCorrelationStep(const char* fileName, double pt1Min = 4.0, double pt1Max = 14.0, double pt2Min = 2.0, Double_t ptAssocMin = 0.5, Double_t ptAssocMax = 8.0, Int_t multBinBegin = 1, Int_t multBinEnd = 5, Int_t vertexBinBegin = 1, Int_t vertexBinEnd = 5, Int_t step = 1, Int_t* trigger = NULL){

    //to guarantee to pick only the bins from the choosen pt on
    pt1Min += 0.01;
    pt1Max -= 0.01;
    pt2Min += 0.01;
    ptAssocMin += 0.01;
    ptAssocMax -= 0.01;
  

  THnBase* trackMixedAll = 0;
  TH2* eventMixedAll = 0;

  //look if the THnBase is already in the cache
  if(lastFileName_corStep && strcmp(lastFileName_corStep, fileName) == 0 && lastPt1Min==pt1Min && lastPt1Max==pt1Max && lastPt2Min==pt2Min && lastPtAssocMin==ptAssocMin && lastPtAssocMax==ptAssocMax && laststep==step){
    Printf("use cache 1");

      trackMixedAll = cacheHistsZVTxMult;
      eventMixedAll = cacheEvtZVTxMult;

  }else if(lastFileName2_corStep && strcmp(lastFileName2_corStep, fileName) == 0 && lastPt1Min2==pt1Min && lastPt1Max2==pt1Max && lastPt2Min2==pt2Min && lastPtAssocMin2==ptAssocMin && lastPtAssocMax2==ptAssocMax && laststep2==step){
    Printf("use cache 2");

      trackMixedAll = cacheHistsZVTxMult2;
      eventMixedAll = cacheEvtZVTxMult2;

      //if it is in no cache continue getting it (by always recomputing it a memory leak is introduced)
  }else{
    loadlibs();

    AliTwoPlusOneContainer* twoPlusOne = (AliTwoPlusOneContainer*) GetTwoPlusOne(fileName);

    twoPlusOne->GetData()->SetPtRange(ptAssocMin, ptAssocMax);
    twoPlusOne->GetData()->SetPt2Min(pt2Min);
    
    twoPlusOne->GetData()->GetHistsZVtxMult((AliUEHist::CFStep)step, AliUEHist::kToward, pt1Min, pt1Max, &trackMixedAll, &eventMixedAll);

    if(next == 1){
      cacheHistsZVTxMult = trackMixedAll;
      cacheEvtZVTxMult = eventMixedAll;

      lastFileName_corStep = fileName;
      lastPt1Min = pt1Min;
      lastPt1Max = pt1Max;
      lastPt2Min = pt2Min;
      lastPtAssocMin = ptAssocMin;
      lastPtAssocMax = ptAssocMax;
      laststep = step;

      next = 2;

    }else if(next ==2){
      cacheHistsZVTxMult2 = trackMixedAll;
      cacheEvtZVTxMult2 = eventMixedAll;

      lastFileName2_corStep = fileName;
      lastPt1Min2 = pt1Min;
      lastPt1Max2 = pt1Max;
      lastPt2Min2 = pt2Min;
      lastPtAssocMin2 = ptAssocMin;
      lastPtAssocMax2 = ptAssocMax;
      laststep2 = step;
     
      next = 1;
    }
  }

  trackMixedAll->GetAxis(3)->SetRange(multBinBegin, multBinEnd);
  trackMixedAll->GetAxis(2)->SetRange(vertexBinBegin, vertexBinEnd);

  TH2D* tracksMixed = trackMixedAll->Projection(1, 0, "E");

  if(trigger!=NULL)
    *trigger = eventMixedAll->Integral(vertexBinBegin, vertexBinEnd, multBinBegin, multBinEnd);

  //TCanvas* c1 = new TCanvas("sameEvent", "sameEvent", gBasisSize+50, gBasisSize+50, gBasisSize, gBasisSize);
  //eventMixedAll->Draw("colz");

 Printf("trigger %i ", (int)eventMixedAll->Integral(vertexBinBegin, vertexBinEnd, multBinBegin, multBinEnd));

  return tracksMixed;
}

void trivial1plus1Division(const char* fileName, double pt1Min = 4.0, double pt1Max = 14.0, double pt2Min = 2.0, Double_t ptAssocMin = 0.5, Double_t ptAssocMax = 8.0, Int_t multBinBegin = 1, Int_t multBinEnd = 5, Int_t vertexBinBegin = 1, Int_t vertexBinEnd = 5, Int_t step = 1, Int_t* trigger = NULL){

  TH2D* tracksSame = getCorrelationStep(fileName, pt1Min, pt1Max, pt2Min, ptAssocMin, ptAssocMax, multBinBegin, multBinEnd, vertexBinBegin, vertexBinEnd, 6, trigger);

  TH2D* tracksMixed = getCorrelationStep(fileName, pt1Min, pt1Max, pt2Min, ptAssocMin, ptAssocMax, multBinBegin, multBinEnd, vertexBinBegin, vertexBinEnd, 9, trigger);

  TCanvas* c1 = new TCanvas("sameEvent", "sameEvent", gBasisSize+50, gBasisSize+50, gBasisSize, gBasisSize);
  tracksSame->DrawCopy("surf1");

  TCanvas* c1 = new TCanvas("mixedEvent", "mixedEvent", 2*gBasisSize+50, gBasisSize+50, gBasisSize, gBasisSize);
  tracksMixed->DrawCopy("surf1");

  TCanvas* c1 = new TCanvas("dividedEvent", "dividedEvent", 3*gBasisSize+50, gBasisSize+50, gBasisSize, gBasisSize);
  tracksSame->Divide(tracksMixed);
  tracksSame->DrawCopy("surf1");

}




void showMixedDist(const char* fileName){
  loadlibs();

  TH1F* mixedDist = getMixedDist(fileName);
  TCanvas* can = new TCanvas();
  mixedDist->DrawCopy("colz");
}

TH1F* getMixedDist(const char* fileName){
  loadlibs();

  //list = (TList*) getList(fileName, "PWGCF_TwoPlusOne/histosTwoPlusOne");
  //list = (TList*) getList(fileName, "PWGCF_TwoPlusOne/addedEvents_");
  list = (TList*) getList(fileName, path);

  TH1F* mixedDist = (TH1F*) list->FindObject("mixedDist");
  return mixedDist;
}



void Plot_oneAxis(const char* fileName, int project)
{
  //kSameNS = 0, kSameAS, kMixedNS, kMixedAS, kMixedCombNS, kMixedCombAS, k1plus1, kBackgroundSameNS, kBackgroundSameAS, kMixed1plus1
  loadlibs();

  TFile::Open(fileName);

  AliTwoPlusOneContainer* h = (AliTwoPlusOneContainer*) GetTwoPlusOne(fileName);

  TCanvas* c1 = new TCanvas("can1", "can1", 1200, 800);
  c1->Divide(2, 1);

  AliCFGridSparse* near_plot = h->GetData()->GetTrackHist(0)->GetGrid((AliUEHist::CFStep) AliTwoPlusOneContainer::k1plus1);
  //near_plot->SetRangeUser(2, gpTMin_T1, gpTMax_T1);
  //near_plot->SetRangeUser(6, gpTMin_T2, gpTMax_T2);
  //near_plot->SetRangeUser(1, gpTMin_assoc, gpTMax_assoc);
  TH1D* tracks_near = near_plot->Project(project);

  c1->cd(1);
  tracks_near->DrawCopy();

  AliCFGridSparse* away_plot = h->GetData()->GetTrackHist(0)->GetGrid((AliUEHist::CFStep) AliTwoPlusOneContainer::kSameNS);
  TH1D* tracks_away = away_plot->Project(project);

  c1->cd(2);
  //TCanvas* c2 = new TCanvas("can2", "can2", 1200, 800);
 tracks_away->SetLineWidth(3);
  tracks_away->DrawCopy();
}  

void test(const char* fileName, double pt1Min, double pt1Max, double pt2Min, double ptAssocMin, double ptAssocMax)
{
  loadlibs();
  Int_t multBinBegin = 1;
  Int_t multBinEnd = 5;
  Int_t trigger_mixed_comb;

  TFile::Open(fileName);

  AliTwoPlusOneContainer* h = (AliTwoPlusOneContainer*) GetTwoPlusOne(fileName);

  h->GetData()->SetPtRange(ptAssocMin, ptAssocMax);
  h->GetData()->SetPt2Min(pt2Min);
  h->GetData()->SetZVtxRange(-1*6.9, 6.9);

  AliUEHist::CFStep step_same = (AliUEHist::CFStep) AliTwoPlusOneContainer::kSameNS;
  AliUEHist::CFStep step_mixed = (AliUEHist::CFStep) AliTwoPlusOneContainer::kMixedNS;

  AliUEHist::CFStep step_1plus1 = (AliUEHist::CFStep) AliTwoPlusOneContainer::k1plus1;

  //GetSumOfRatios2(mixed, step, region, ptLeadMin, ptLeadMax, multBinBegin, multBinEnd, normalizePerTrigger, stepForMixed)
  TH2D* h2_etaPhi;
  TH2D* h2_etaPhi_mixedComb;

    h2_etaPhi_mixedComb = (TH2D*) h->GetData()->GetSumOfRatios2(h->GetData(), step_1plus1, 0, pt1Min, pt1Max, multBinBegin, multBinEnd, kFALSE, step_mixed, &trigger_mixed_comb);

 TCanvas* c1 = new TCanvas("can1", "can1", 1200, 800);
 h2_etaPhi_mixedComb->DrawCopy("surf1");
}

void showEvents_2D(const char* fileName, Double_t ptLeadMin, Double_t ptLeadMax, Double_t pt2Min, Double_t pt2Max){
  loadlibs();
  
  AliUEHist::CFStep step = (AliUEHist::CFStep) AliTwoPlusOneContainer::kMixedCombAS;

  AliTwoPlusOneContainer* twoPlusOne = (AliTwoPlusOneContainer*) GetTwoPlusOne(fileName);

  THnBase* trackAll = 0;
  TH2* eventAll = 0;


  twoPlusOne->GetData()->SetPt2Min(pt2Min);
  twoPlusOne->GetData()->SetPt2Max(pt2Max);
  twoPlusOne->GetData()->GetHistsZVtxMult((AliUEHist::CFStep)step, AliUEHist::kToward, ptLeadMin, ptLeadMax, &trackAll, &eventAll);

  TCanvas* c1 = new TCanvas("can1", "can1", 1200, 800);
  eventAll->Draw("surf1");
}

void showAsymmetry(const char* fileName, Int_t mixed){
  loadlibs();

  TFile::Open(fileName);
  AliTwoPlusOneContainer* h = (AliTwoPlusOneContainer*) GetTwoPlusOne(fileName);

  TH1F* asymm_same = h->GetAsymmetry();
  TH1F* asymm_mixed = h->GetAsymmetryMixed();

  TCanvas* c1 = new TCanvas("can1", "can1", 1200, 800);
  if(mixed==0)
    asymm_same->DrawCopy();
  else if(mixed==1)
    asymm_mixed->DrawCopy();
  else if(mixed==2){
    asymm_same->DrawCopy();
    asymm_mixed->SetLineColor(kRed);
    asymm_mixed->Scale(1.0/80);
    asymm_mixed->DrawCopy("same");
  }
}

void showTriggerPt(const char* fileName){
  loadlibs();

  TFile::Open(fileName);
  AliTwoPlusOneContainer* h = (AliTwoPlusOneContainer*) GetTwoPlusOne(fileName);

  TCanvas* c1 = new TCanvas("can1", "can1", 1200, 800);
  h->GetTriggerPt()->DrawCopy("colz");
}

