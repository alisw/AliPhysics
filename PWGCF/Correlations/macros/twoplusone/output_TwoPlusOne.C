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
Float_t g_phi_bin = 0.32;//bins of 0.174532
Float_t g_eta_bin = 0.39;//bins of 0.2

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

  list = (TList*) getList(fileName, "PWGCF_TwoPlusOne/histosTwoPlusOne");
  
  AliTwoPlusOneContainer* container = (AliTwoPlusOneContainer*) list->FindObject("AliTwoPlusOneContainer");

  if (container->GetData()->GetTrackHist(0)->GetGrid(6)->GetGrid()->GetNbins() == 0)
  {
    Printf("We have %d axes", ((AliTHn*) container->GetData()->GetTrackHist(0)->GetNVar()));
    
    ((AliTHn*) container->GetData()->GetTrackHist(0))->FillParent();
    ((AliTHn*) container->GetData()->GetTrackHist(0))->DeleteContainers();
  }

  if(lastFileName2 && strcmp(lastFileName2, fileName) == 0){
    cacheEvent2 = container;
    return cacheEvent2;
  }

  cacheEvent = container;

  return cacheEvent;
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
//flow from 1.0 > |eta| > 1.6 is subtracted from the rest of the histogram
//the flow in the TH2 is subtracted
void subtractFlow(TH2* etaPhi){
  int firstbin = etaPhi->GetYaxis()->FindBin(1.01);
  int lastbin = etaPhi->GetYaxis()->FindBin(1.39);

  TH1D* etaPhi_proj = etaPhi->ProjectionX("_px", firstbin, lastbin);
  int usedBins = lastbin - firstbin + 1;

  firstbin = etaPhi->GetYaxis()->FindBin(-1.39);
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

  //use fit method to subtract the baseline if the subtractFlow method doesn't work because of poor statistics
Double_t subtractBaseline(TH1D* h1_phi){
  
  TF1* fit = new TF1("fit", "[0]", -TMath::Pi()/2+0.01, TMath::Pi()/2-0.01);

  h1_phi->Fit("fit", "N", "", -TMath::Pi()/2+0.01, -1.0);
  Double_t par0 = fit->GetParameter(0);
 
  h1_phi->Fit("fit", "N", "", 1.0, TMath::Pi()/2-0.01);
  Double_t par1 = fit->GetParameter(0);

  //this is done so that the fit is not drawn
  //h1_phi->Fit("fit", "0", "", 1.0, TMath::Pi()/2-0.01);
 
  Double_t subtract = (par0+par1)/2;//h1_phi->GetMinimum();

  
  for(int i=0; i<h1_phi->GetNbinsX(); i++)
    h1_phi->SetBinContent(i, h1_phi->GetBinContent(i) - subtract);

  return subtract;
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

void showMultCompare(const char* filename, Int_t pt1Minimum = 4.0, Int_t pt1Maximum = 14.0, Int_t pt2Minimum = 2.0, Double_t ptAssocMin = 0.5, Double_t ptAssocMax = 8.0, Int_t side = 0, Int_t subtractMixedComb = 1, Int_t subtractFlow = 1){

  TLegend *leg  = new TLegend(0.65,0.7,0.85,0.9);
  leg->SetFillColor(10);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextFont(gStyle->GetTextFont());
  leg->SetTextSize(gStyle->GetTextSize()*1.0);

  TH1D* central = getAnalysis(filename, pt1Minimum, pt1Maximum, pt2Minimum, ptAssocMin, ptAssocMax, 6.9, 1, 5, side, 1, 1, 0, subtractMixedComb, subtractFlow);

  TCanvas* multCompare = new TCanvas("multCompare", "multCompare", gBasisSize+100, gBasisSize+50, 2*gBasisSize, 2*gBasisSize);
  central->SetLineColor(kBlue);
  leg->AddEntry(central->Clone(),"central","l");
  subtractBaseline(central);
  central->DrawCopy();

  TH1D* semiCentral = getAnalysis(filename, pt1Minimum, pt1Maximum, pt2Minimum, ptAssocMin, ptAssocMax, 6.9, 9, 10, side, 1, 1, 0, subtractMixedComb, subtractFlow);

  semiCentral->SetLineColor(kRed);
  subtractBaseline(semiCentral);
  semiCentral->DrawCopy("same");  
  leg->AddEntry(semiCentral->Clone(),"semiCentral","l");
  /*
  TH1D* peripheral = getAnalysis(filename, pt1Minimum, pt1Maximum, pt2Minimum, ptAssocMin, ptAssocMax, 6.9, 11, 14, side, 1, 1, 0, subtractMixedComb);

  peripheral->SetLineColor(kGreen);
   leg->AddEntry(peripheral,"peripheral","l");
  peripheral->DrawCopy("same");
  */

  leg->Draw("same");

  //show same and mixed event of one vertex and multiplicity bin
  //getRawAnalysis(filename, pt1Minimum, pt1Maximum, pt2Minimum, ptAssocMin, ptAssocMax, 3, 3);
}

void showAnalysis(const char* fileName, double pt1Min = 4.0, double pt1Max = 14.0, double pt2Min = 2.0, Double_t ptAssocMin = 0.5, Double_t ptAssocMax = 8.0, Int_t multBinBegin = 1, Int_t multBinEnd = 5, Int_t subtractMixedComb = 0, Int_t subtractFlow = 0){

  TLegend* leg = getLegend();
  Int_t showPlots = 1;

  TH1* near = (TH1*)getAnalysis(fileName, pt1Min, pt1Max, pt2Min, ptAssocMin, ptAssocMax, 6.9, multBinBegin, multBinEnd, 0, 1, 1, showPlots, subtractMixedComb, subtractFlow)->Clone();
  /*
  subtractFlow(h2_etaPhi);
  h2_etaPhi->GetYaxis()->SetRangeUser(-1.6, 1.6);//corresponds to the really used eta range
 
  int firstbin = h2_etaPhi->GetYaxis()->FindBin(-0.59);
  int lastbin = h2_etaPhi->GetYaxis()->FindBin(0.59);

  TH1D* h1_phi = projectToTH1D(h2_etaPhi, "fully corrected", firstbin, lastbin, 1.0);//last number is triggers
  */
  leg->AddEntry(near->Clone(), "near side", "l");
  if(true) return;
  TH1* away = getAnalysis(fileName, pt1Min, pt1Max, pt2Min, ptAssocMin, ptAssocMax, 6.9, multBinBegin, multBinEnd, 1, 2, 1, showPlots, subtractMixedComb, subtractFlow);
  leg->AddEntry(away, "away side", "l");



  TCanvas* c1 = new TCanvas("can", "can", gBasisSize+50, gBasisSize+50, gBasisSize, gBasisSize);
  near->DrawCopy();

  away->SetLineColor(kRed);
  away->DrawCopy("same");
  leg->Draw("same");
}

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
void peakDifference(const char* fileName, double pt1Min = 4.0, double pt1Max = 14.0, double pt2Min = 2.0, Int_t multBinBegin = 1, Int_t multBinEnd = 5, Int_t side = 0, Int_t yPos = 0, Int_t mode = 0,  Int_t subtractMixedComb = 1, Int_t subtractFlow = 1){

  static const int pt_assoc_bins_number = 8;
  Double_t pt_assoc_bins[pt_assoc_bins_number] = {0.5, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0};
  Double_t near_content[pt_assoc_bins_number-1];
  Double_t away_content[pt_assoc_bins_number-1];
  Double_t onePlusOne_content[pt_assoc_bins_number-1];
  Double_t near_error[pt_assoc_bins_number-1];
  Double_t away_error[pt_assoc_bins_number-1];
  Double_t onePlusOne_error[pt_assoc_bins_number-1];

  Double_t pt_assoc_bin_center[pt_assoc_bins_number-1] = {0.75, 1.5, 2.5, 3.5, 5.0, 7.0, 9.0};
  Double_t pt_assoc_bin_error[pt_assoc_bins_number-1] = {0.25, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0};

  for(int i=0; i<pt_assoc_bins_number-1; i++){
    TH1* near;
    TH1* away;
    TH1* onePlusOne;
    TCanvas* can;

    if(mode==0){
      near = (TH1*)getAnalysis(fileName, pt1Min, pt1Max, pt2Min, pt_assoc_bins[i], pt_assoc_bins[i+1], 6.9, multBinBegin, multBinEnd, 0, 1, 1, 0, subtractMixedComb, subtractFlow)->Clone();
      away = (TH1*)getAnalysis(fileName, pt1Min, pt1Max, pt2Min, pt_assoc_bins[i], pt_assoc_bins[i+1], 6.9, multBinBegin, multBinEnd, 1, 2, 1, 0, subtractMixedComb, subtractFlow)->Clone();
      onePlusOne = (TH1*)getAnalysis(fileName, pt1Min, pt1Max, pt2Min, pt_assoc_bins[i], pt_assoc_bins[i+1], 6.9, multBinBegin, multBinEnd, 6, 2, 1, 0, subtractMixedComb, subtractFlow);
      

      can = new TCanvas(Form(" %i, mult %i-%i ", i, multBinBegin, multBinEnd), Form(" %i, mult %i-%i ", i, multBinBegin, multBinEnd), i*gBasisSize+100, yPos*(gBasisSize+50), gBasisSize, gBasisSize);
    }else if(mode==1){
      near = (TH1*)getAnalysis(fileName, pt1Min, pt1Max, pt2Min, pt_assoc_bins[i], pt_assoc_bins[i+1], 6.9, 1, 5, side, 1, 1, 0, subtractMixedComb, subtractFlow)->Clone();
      away = (TH1*)getAnalysis(fileName, pt1Min, pt1Max, pt2Min, pt_assoc_bins[i], pt_assoc_bins[i+1], 6.9, 9, 10, side, 2, 1, 0, subtractMixedComb, subtractFlow);

      can = new TCanvas(Form(" %i, side %i ", i, side), Form(" %i, side %i ", i, side), i*gBasisSize+50, yPos*(gBasisSize+50), gBasisSize, gBasisSize);
    }else{
      Printf("mode does not exist!");
      return;
    }
    away->SetLineColor(kRed);
    
    Int_t bin_start = near->FindBin(-1*g_phi_bin);
    Int_t bin_end = near->FindBin(g_phi_bin);

    near_content[i] = near->IntegralAndError(bin_start, bin_end, near_error[i]);
    away_content[i] = away->IntegralAndError(bin_start, bin_end, away_error[i]);
    if(mode==0)
      onePlusOne_content[i] = onePlusOne->IntegralAndError(bin_start, bin_end, onePlusOne_error[i]);

    ((TH1*)(near->Clone()))->DrawCopy();
    ((TH1*)(away->Clone()))->DrawCopy("same");

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
    onePlusOne_content[i] /= 2*pt_assoc_bin_error[i];
    onePlusOne_error[i] /= 2*pt_assoc_bin_error[i];
  }

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
  graph_onePlusone->GetXaxis()->SetTitle("p_{T,assoc} (GeV/c)");
  graph_onePlusone->GetYaxis()->SetTitle("1/N dN/dpT");
  graph_onePlusone->SetMarkerSize(2);
  graph_onePlusone->SetLineWidth(2);
  graph_onePlusone->SetMarkerColor(kGreen);
  graph_onePlusone->SetLineColor(kGreen);
  graph_onePlusone->SetMarkerStyle(22);
  graph_near->Draw("AP");
  graph_away->Draw("P");
  graph_onePlusone->Draw("P");

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


//does the reading out of the results
//divides the same event by the mixed events
//subtracts the mixed combinatorics
TH1D* getAnalysis(const char* fileName, double pt1Min = 4.0, double pt1Max = 14.0, Double_t pt2Min = 2.0, Double_t ptAssocMin = 0.5, Double_t ptAssocMax = 8.0, Double_t setVertex = 7, Int_t multBinBegin = 1, Int_t multBinEnd = 5, Int_t step = 0, Int_t posX = 1, Int_t posY = 1, Int_t showPlots = 0, Int_t subtractMixedComb = 0, Int_t subtractFlow =0)
{
  loadlibs();
  //to guarantee to pick only the bins from the choosen pt on
  pt1Min += 0.01;
  pt1Max -= 0.01;
  pt2Min += 0.01;
  ptAssocMin += 0.01;
  ptAssocMax -= 0.01;
  
  AliTwoPlusOneContainer* h = (AliTwoPlusOneContainer*) GetTwoPlusOne(fileName);

  AliUEHist::CFStep step_same = step;//(AliUEHist::CFStep) AliTwoPlusOneContainer::kSameNS;
  AliUEHist::CFStep step_mixed = step+2;//(AliUEHist::CFStep) AliTwoPlusOneContainer::kMixedNS;
  AliUEHist::CFStep step_mixedComb = step+4;//(AliUEHist::CStep) AliTwoPlusOneContainer::kMixedCombNS;
  AliUEHist::CFStep step_backgroundSame = step+7;//(AliUEHist::CStep) AliTwoPlusOneContainer::kBackgroundSameNS;

  AliUEHist::CFStep step_1plus1 = (AliUEHist::CFStep) AliTwoPlusOneContainer::k1plus1;
  AliUEHist::CFStep step_1plus1_mixed = (AliUEHist::CFStep) AliTwoPlusOneContainer::kMixed1plus1;

  h->GetData()->SetPtRange(ptAssocMin, ptAssocMax);
  h->GetData()->SetPt2Min(pt2Min);
  h->GetData()->SetZVtxRange(-1*setVertex, setVertex);
 
  //GetSumOfRatios2(mixed, step, region, ptLeadMin, ptLeadMax, multBinBegin, multBinEnd, normalizePerTrigger, stepForMixed)
  TH2D* h2_etaPhi;
  TH2D* h2_etaPhi_mixedComb;
  TH2D* h2_etaPhi_backgroundSame;

  //TH2D* h2_etaPhi2;
  //TH2D* h2_etaPhi_mixedComb2;

  TH1D* h1_phi_cloneProject;

  if(step>=2||!subtractMixedComb){
    Int_t trigger;
    if(step!=6)
      h2_etaPhi = (TH2D*) h->GetData()->GetSumOfRatios2(h->GetData(), step_same, 0, pt1Min, pt1Max, multBinBegin, multBinEnd, kTRUE, step_mixed, &trigger);
    else  if(step==6)
      h2_etaPhi = (TH2D*) h->GetData()->GetSumOfRatios2(h->GetData(), step_1plus1, 0, pt1Min, pt1Max, multBinBegin, multBinEnd, kTRUE, step_1plus1_mixed, &trigger);
  }else if(step<2 && subtractMixedComb){
    Int_t trigger_same;
    Int_t trigger_mixed_comb;
    Int_t trigger_background_same;

    h2_etaPhi = (TH2D*) h->GetData()->GetSumOfRatios2(h->GetData(), step_same, 0, pt1Min, pt1Max, multBinBegin, multBinEnd, kFALSE, step_mixed, &trigger_same);
    TH2D* h2_etaPhi_clone = h2_etaPhi->Clone();
    if(subtractFlow)
      subtractFlow(h2_etaPhi_clone);
    h2_etaPhi_clone->Scale(1.0/trigger_same);

    //don't need to use getMixedComb_scaled_backgroundSame, see compareScaledMixedComb which shows that both methods are identical but here getSumOfRatios is used which makes it easier (and safer against errors)
   h2_etaPhi_mixedComb = (TH2D*) h->GetData()->GetSumOfRatios2(h->GetData(), step_mixedComb, 0, pt1Min, pt1Max, multBinBegin, multBinEnd, kFALSE, step_mixed, &trigger_mixed_comb);

    h2_etaPhi_backgroundSame = (TH2D*) h->GetData()->GetSumOfRatios2(h->GetData(), step_backgroundSame, 0, pt1Min, pt1Max, multBinBegin, multBinEnd, kFALSE, step_mixed, &trigger_background_same);

    //trigger_background_same *= 2;//error in the plots for runs 968 and 973
    double trigger_ratio = (double)trigger_background_same;
    if(trigger_mixed_comb>0)
      trigger_ratio /= trigger_mixed_comb;
    Printf("trigger mixed comb: %i, trigger background same %i, ratio %f", trigger_mixed_comb, trigger_background_same, trigger_ratio);
    h2_etaPhi_mixedComb->Scale(trigger_ratio);
    trigger_mixed_comb = trigger_background_same;
    Printf("trigger same: %i, trigger mixed comb %i ", trigger_same, trigger_mixed_comb);

    int firstbin = h2_etaPhi->GetYaxis()->FindBin(-1*g_eta_bin);
    int lastbin = h2_etaPhi->GetYaxis()->FindBin(g_eta_bin);

    if(showPlots==1){
      TCanvas* can_1 = new TCanvas("sameEvent near side", "sameEvent near side", 100, 50, gBasisSize, gBasisSize);
      h2_etaPhi->DrawCopy("surf1");
    
      //TH1D* h1_phi_1 = projectToTH1D(h2_etaPhi, "same near", firstbin, lastbin, trigger_same);
      TH1D* h1_phi_1 = projectToTH1D(h2_etaPhi, "same near", firstbin, lastbin, 1.0);
      h1_phi_1->SetLineColor(kRed);
      //h1_phi_1->DrawCopy();
      
      TLegend *leg  = getLegend();
      leg->AddEntry(h1_phi_1,"same event","l");
    
      TCanvas* can_2 = new TCanvas("sameEvent mixed comb", "sameEvent mixed comb", gBasisSize+100, 50, gBasisSize, gBasisSize);
      h2_etaPhi_mixedComb->DrawCopy("surf1");
      //TH1D* h1_phi_2 = projectToTH1D(h2_etaPhi_mixedComb, "mixedComb near", firstbin, lastbin, trigger_mixed_comb);
      TH1D* h1_phi_2 = projectToTH1D(h2_etaPhi_mixedComb, "mixedComb near", firstbin, lastbin, 1.0);
      Printf(Form("trigger mixed_comb", trigger_mixed_comb));
      //h1_phi_2->DrawCopy();
      leg->AddEntry(h1_phi_2->Clone(),"mixed comb","l");
      //leg->Draw("same");
      //if(true) return;
      TCanvas* can_2 = new TCanvas("sameEvent flow subtracted", "sameEvent flow subtracted", 100, 2*(gBasisSize+50), gBasisSize, gBasisSize);
      //TH2D* h2_etaPhi_clone = h2_etaPhi->Clone(); // is already cloned earlier
      //subtractFlow(h2_etaPhi_clone);
      
      h1_phi_cloneProject = projectToTH1D(h2_etaPhi_clone, "same flow subtracted", firstbin, lastbin, trigger_same);
      h1_phi_cloneProject->DrawCopy();
    
      TCanvas* can_3 = new TCanvas("mixed comb subtracted", "mixed comb subtracted", 2*gBasisSize+100, 2*(gBasisSize+50), gBasisSize, gBasisSize);
      h2_etaPhi->Add(h2_etaPhi_mixedComb, -1);
      //h2_etaPhi2->Add(h2_etaPhi_mixedComb2, -1);
    
      if(trigger_same-trigger_mixed_comb>0)
	h2_etaPhi->Scale(1.0/(trigger_same-trigger_mixed_comb));
      //h2_etaPhi2->Scale(1.0/(trigger_same2-trigger_mixed_comb2));

      //h2_etaPhi->DrawCopy("surf1");
      h1_phi_3 = projectToTH1D(h2_etaPhi, "same minus mixedComb", firstbin, lastbin, 1.0);
      //h1_phi_3_b = projectToTH1D(h2_etaPhi2, "same minus mixedComb away", firstbin, lastbin, 1.0);
      h1_phi_3->SetLineColor(kRed);
      h1_phi_3->DrawCopy();
      //h1_phi_3_b->DrawCopy("same");
    
      TLegend* leg2 = getLegend();
      leg2->AddEntry(h1_phi_3, "near side", "l");
      //leg2->AddEntry(h1_phi_3_b, "away side", "l");
      leg2->Draw("same");
    }else{
      
       h2_etaPhi->Add(h2_etaPhi_mixedComb, -1);
    
      if(trigger_same-trigger_mixed_comb>0)
	h2_etaPhi->Scale(1.0/(trigger_same-trigger_mixed_comb));
      //using backgroundSame as correction, can only be used if the used phi and eta bins of the associated particles are very close to the trigger particles
      /*
      h2_etaPhi->Add(h2_etaPhi_backgroundSame, -1);
      h2_etaPhi->Scale(1.0/(trigger_same-trigger_background_same));
      */

      if(subtractFlow)
	subtractFlow(h2_etaPhi_mixedComb);
      if(trigger_mixed_comb>0)
	h2_etaPhi_mixedComb->Scale(1.0/trigger_mixed_comb);
    }
    
  }else{
    Printf("cannot subtract mixed combinatorics from step %i ", step);
    return 0;
  }

  if(subtractFlow)
    subtractFlow(h2_etaPhi);
  h2_etaPhi->GetYaxis()->SetRangeUser(-1.6, 1.6);//corresponds to the really used eta range
 
  int firstbin = h2_etaPhi->GetYaxis()->FindBin(-1*g_eta_bin);
  int lastbin = h2_etaPhi->GetYaxis()->FindBin(g_eta_bin);

  TH1D* h1_phi = projectToTH1D(h2_etaPhi, "fully corrected", firstbin, lastbin, 1.0);
  //TH1D* h1_phi = projectToTH1D(h2_etaPhi_mixedComb, "fully corrected", firstbin, lastbin, 1.0);//shows background distributions

  if(showPlots>0){
    TLegend *leg  = getLegend();
    
    TCanvas* can_4 = new TCanvas(Form("fully corrected, %i, %i ", posX, posY), Form("fully corrected, %i, %i ", posX, posY), posX*gBasisSize+100, 50, gBasisSize, gBasisSize);
    h2_etaPhi->GetYaxis()->SetRangeUser(-1.4, 1.4);
    h2_etaPhi->DrawCopy("surf1");
    //h2_etaPhi_mixedComb->DrawCopy("surf1");
    //h1_phi->DrawCopy();
    //h1_phi_cloneProject->SetLineColor(kRed);
    //h1_phi_cloneProject->DrawCopy("same");

    leg->AddEntry(h1_phi,"with mixed comb","l");
    leg->AddEntry(h1_phi_cloneProject,"without mixed comb","l");

    //leg->Draw("same");
  }

  if(showPlots==2){
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

  return h1_phi;
}


TH1D* projectToTH1D(TH2D* etaPhi, char* name, int firstbin, int lastbin, int trigger){
  //  Printf(Form("name %s", name));

  TH1D* h1_phi_1 = etaPhi->ProjectionX(Form("px_%s", name), firstbin, lastbin);
  h1_phi_1->SetLineWidth(3);
  h1_phi_1->SetStats(kFALSE);
  h1_phi_1->Scale(1.0/trigger);
  h1_phi_1->Scale(1.0/(double)(lastbin-firstbin+1));
  h1_phi_1->GetXaxis()->SetRangeUser(-TMath::Pi()/2+0.01, TMath::Pi()/2-0.01);
  //h1_phi_1->GetXaxis()->SetRangeUser(-TMath::Pi()/2+0.01, 3*TMath::Pi()/2-0.01);
  h1_phi_1->SetYTitle("1/N \\ dN/(d \\Delta \\varphi)");

  //symmetrize(h1_phi_1);

  return h1_phi_1;
}

TLegend* getLegend(){
  TLegend *leg  = new TLegend(0.65,0.7,0.85,0.9);
  leg->SetFillColor(10);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextFont(gStyle->GetTextFont());
  leg->SetTextSize(gStyle->GetTextSize()*1.0);
  return leg;
}


//calculates the mixed comb by scaling it up to the amount of trigger in same background
//this method can probably not be used because the statistics is too smal
//within one vertex and multiplicity bin there are 12 to 30 trigger combinations -> the error on the correction factor is too high.
TH2D* getMixedComb_scaled_backgroundSame(const char* fileName, double pt1Min = 4.0, double pt1Max = 14.0, double pt2Min = 4.0, double ptAssocMin = 2.0, double ptAssocMax = 4.0, Int_t multBinBegin = 1, Int_t multBinEnd = 5, Double_t setVertex = 7, Int_t step = 4, Int_t* trigger = NULL){

 loadlibs();

  printf("#############################");
  printf(Form("step %i", step));
  printf("#############################");
  AliTwoPlusOneContainer* h = (AliTwoPlusOneContainer*) GetTwoPlusOne(fileName);
  h->GetData()->SetPtRange(ptAssocMin, ptAssocMax);
  h->GetData()->SetPt2Min(pt2Min);

  TH2D* return_mixed_comb = 0;
  
  if(step!=4 && step!=5){
    Printf("#############################");
    Printf("ERROR wrong step in getMixedComb_scaled_backgroundSame");
    Printf("#############################");
    return return_mixed_comb;
  }

  THnBase* trackMixedCombAll = 0;
  TH2* eventMixedCombAll = 0;
  h->GetData()->GetHistsZVtxMult((AliUEHist::CFStep)step, AliUEHist::kToward, pt1Min, pt1Max, &trackMixedCombAll, &eventMixedCombAll);
  
  THnBase* trackBackgroundSameAll = 0;
  TH2* eventBackgroundSameAll = 0;
  h->GetData()->GetHistsZVtxMult((AliUEHist::CFStep)(step+3), AliUEHist::kToward, pt1Min, pt1Max, &trackBackgroundSameAll, &eventBackgroundSameAll);
  
  TAxis* vertexAxis = trackMixedCombAll->GetAxis(2);
  Int_t vertexBinBegin = vertexAxis->FindBin(-1*setVertex);
  Int_t vertexBinEnd = vertexAxis->FindBin(setVertex);
  
  Printf(Form("vertex bin begin %i, end %i ", vertexBinBegin, vertexBinEnd));

  for(Int_t multBin=multBinBegin; multBin<=multBinEnd; multBin++){
    for(Int_t vertexBin=vertexBinBegin; vertexBin<=vertexBinEnd; vertexBin++){
      Printf(" multBin %i, vertexBin %i ", multBin, vertexBin);

      trackMixedCombAll->GetAxis(3)->SetRange(multBin, multBin);
      trackMixedCombAll->GetAxis(2)->SetRange(vertexBin, vertexBin);

      trackBackgroundSameAll->GetAxis(3)->SetRange(multBin, multBin);
      trackBackgroundSameAll->GetAxis(2)->SetRange(vertexBin, vertexBin);

      TH2D* tracks_mixedComb = trackMixedCombAll->Projection(1, 0, "E");
      Int_t trigger_mixed_comb = eventMixedCombAll->Integral(vertexBin, vertexBin, multBin, multBin);

      TH2D* tracks_background_same = trackBackgroundSameAll->Projection(1, 0, "E");
      Int_t trigger_background_same = eventBackgroundSameAll->Integral(vertexBin, vertexBin, multBin, multBin);
 
      //this has to be done because of a bug in the AliRoot code
      //because of an absolute cast in the dphi computation only half of this background is considered
      //wrong runs are 968 and 973
      //trigger_background_same *= 2;

      Double_t trigger_ratio = 0;
      if(trigger_background_same>0)
	trigger_ratio = (Double_t)trigger_background_same/trigger_mixed_comb;
      Printf("(getMixedComb_scaled_backgroundSame) Triggers: mixed comb %i, background same %i, trigger ratio %f ", trigger_mixed_comb, trigger_background_same, trigger_ratio);
      
      double b1 = tracks_mixedComb->GetBinContent(tracks_mixedComb->FindBin(0.001,0.001));
      double b2 = tracks_mixedComb->GetBinContent(tracks_mixedComb->FindBin(0.001,-0.001));
      double b3 = tracks_mixedComb->GetBinContent(tracks_mixedComb->FindBin(-0.001,0.001));
      double b4 = tracks_mixedComb->GetBinContent(tracks_mixedComb->FindBin(-0.001,-0.001));

      double b1_bsame = tracks_background_same->GetBinContent(tracks_background_same->FindBin(0.001,0.001));
      double b2_bsame = tracks_background_same->GetBinContent(tracks_background_same->FindBin(0.001,-0.001));
      double b3_bsame = tracks_background_same->GetBinContent(tracks_background_same->FindBin(-0.001,0.001));
      double b4_bsame = tracks_background_same->GetBinContent(tracks_background_same->FindBin(-0.001,-0.001));

      Double_t content_ratio = 0;
      if(b1 + b2 + b3 + b4 > 0)
	content_ratio = (b1_bsame + b2_bsame + b3_bsame + b4_bsame)/(b1 + b2 + b3 + b4);
      Printf(" trigger_ratio %f, content_ratio %f", trigger_ratio, content_ratio);
      
      //Printf("Maximum: mixed comb %f, background same %f, normalization ratio %f ", content_mixed_comb, content_background_same, content_ratio);
      
      tracks_mixedComb->Scale(trigger_ratio);
      //tracks_mixedComb->Scale(content_ratio);

      //correct event with mixed event
      TH2D* mixedEvent =  getMixedEvent(fileName, pt1Min, pt1Max, pt2Min, ptAssocMin, ptAssocMax, multBin, multBin, vertexBin, vertexBin, (AliUEHist::CFStep)(step-2), 0);
      //TH2D* mixedEvent =  getMixedEvent(fileName, pt1Min, pt1Max, pt2Min, ptAssocMin, ptAssocMax, multBin, multBin, vertexBinBegin, vertexBinEnd, (AliUEHist::CFStep)(step-2), 1);
      tracks_mixedComb->Divide(mixedEvent);

      //TCanvas* can_1 = new TCanvas("mixed comb", "mixed comb", 100, 50, gBasisSize, gBasisSize);
      //tracks_mixedComb->DrawCopy("surf1");
      //TCanvas* can_2 = new TCanvas("mixed", "mixed", gBasisSize+100, 50, gBasisSize, gBasisSize);
      //mixedEvent->DrawCopy("surf1");

      if(!return_mixed_comb)
	return_mixed_comb = (TH2D*) tracks_mixedComb->Clone("total tracks");
      else
	return_mixed_comb->Add(tracks_mixedComb);

      if(trigger!=NULL)
	*trigger += trigger_background_same ;
    }
  }

  // normalizate to dphi widht
  Float_t normalization = return_mixed_comb->GetXaxis()->GetBinWidth(1);
  return_mixed_comb->Scale(1.0 / normalization);

  return return_mixed_comb;
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

  int firstbin = h2_etaPhi->GetYaxis()->FindBin(-1*g_eta_bin);
  int lastbin = h2_etaPhi->GetYaxis()->FindBin(g_eta_bin);

  Printf("firstbin %i, lastbin %i, \n", firstbin, lastbin);

  TH1D* h1_phi = h2_etaPhi->ProjectionX(Form("px_%i_%i", posX, posY), firstbin, lastbin);
  
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

  //Double_t significance = signal2/h1_phi->GetBinError((bin_1+bin_2)/2);
  //Double_t significance_err = signal2_err/(TMath::Sqrt(sum));
  Double_t significance = signal2/(TMath::Sqrt(sum));
  Double_t significance_err = TMath::Sqrt(TMath::Power(signal2_err/(TMath::Sqrt(sum)), 2) + TMath::Power(signal2/(TMath::Sqrt(sum*sum*sum)*sum_err), 2));

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


//this method is not further maintained after the copy into this file
//it was used to compare the different kind of background distributions
//this method compares the expected amount of background trigger combinations for 3 methods:
//1. get background from the arithmetical average between delta phi = pi/3 and 2pi/3
//2. mixed combinatorics (they are lacking some combinations because the buffers have to be filled first)
//3. background same, search at delta phi = pi/2 for trigger particles the same way as it is done in mixed comb
void printBackground2plus1(const char* fileName, double pt1Min = 4.0, double pt1Max = 14.0, Double_t pt2Min = 2.0, Double_t ptAssocMin = 2.0, Double_t ptAssocMax = 4.0, Double_t setVertex = 7, Int_t multBinBegin = 1, Int_t multBinEnd = 5, Int_t step = 0){

  loadlibs();
  TFile::Open(fileName);

  //to guarantee to pick only the bins from the choosen pt on
  pt1Min += 0.01;
  pt1Max -= 0.01;
  pt2Min += 0.01;
  ptAssocMin += 0.01;
  ptAssocMax -= 0.01;

  Double_t alpha_radius = 0.393;// =pi/8 //0.2;//use radius which was used by the mixed comb and background same method in the loaded file

  AliTwoPlusOneContainer* h = (AliTwoPlusOneContainer*) GetTwoPlusOne(fileName);
  AliUEHist::CFStep step_1plus1 = (AliUEHist::CFStep) AliTwoPlusOneContainer::k1plus1;
  AliUEHist::CFStep step_mixed = (AliUEHist::CFStep) AliTwoPlusOneContainer::kMixedAS;
  AliUEHist::CFStep step_mixedComb = (AliUEHist::CFStep) AliTwoPlusOneContainer::kMixedCombAS;
  AliUEHist::CFStep step_backgroundSame = (AliUEHist::CFStep) AliTwoPlusOneContainer::kBackgroundSameAS;

  //get error from the 1+1 correlation
  
  THnBase* trackSameAll = 0;
  TH2* eventSameAll = 0;
  //h->GetData()->SetPtRange(pt2Min, pt1Max);
  h->GetData()->SetPtRange(ptAssocMin, ptAssocMax);
  //h->GetData()->SetPtRange(6, 10);
  h->GetData()->SetPt2Min(pt2Min);
  h->GetData()->SetZVtxRange(-1*setVertex, setVertex);
  /*
  h->GetData()->GetHistsZVtxMult(step_1plus1, AliUEHist::kToward, pt1Min, pt1Max, &trackSameAll, &eventSameAll);
  trackSameAll->GetAxis(3)->SetRange(multBinBegin, multBinEnd);

  TAxis* vertexAxis = trackSameAll->GetAxis(2);
  Int_t vertexBinBegin = vertexAxis->FindBin(-1*setVertex);
  Int_t vertexBinEnd = vertexAxis->FindBin(setVertex);

  trackSameAll->GetAxis(2)->SetRange(vertexBinBegin, vertexBinEnd);

  TH2* tracksSame = trackSameAll->Projection(1, 0, "E");
  TH1D* h1_phi = tracksSame->ProjectionX();
  h1_phi->Scale(1/h1_phi->GetXaxis()->GetBinWidth(1));
  
  */
  
  TH2* h2_etaPhi = h->GetData()->GetSumOfRatios2(h->GetData(), step_1plus1, 0, pt1Min, pt1Max, multBinBegin, multBinEnd, kFALSE, step_mixed);

  int firstbin = -1;//h2_etaPhi->GetYaxis()->FindBin(-0.99);
  int lastbin = -1;//h2_etaPhi->GetYaxis()->FindBin(0.99);
  TH1D* h1_phi = h2_etaPhi->ProjectionX("px", firstbin, lastbin);
  

  //undo phi normalization
  //h1_phi->Scale(h1_phi->GetXaxis()->GetBinWidth(1));
  /*
  TF1* fit = new TF1("fit", "[0]", 1.01, 1.99);
  h1_phi->Fit("fit", "N", "", 1.01, 1.99);
  Double_t par0 = fit->GetParameter(0);
  Double_t err0 = fit->GetParError(0);
  */
  Int_t bin_start = h1_phi->FindBin(TMath::Pi()/3+0.01);
  Int_t bin_end = h1_phi->FindBin(2*TMath::Pi()/3-0.01);
  //Int_t bin_start = h1_phi->FindBin(1.37);
  //Int_t bin_end = h1_phi->FindBin(1.57);
  Double_t sum = 0;
  Double_t err = 0;

  for(int i=bin_start; i<=bin_end; i++){
    sum += h1_phi->GetBinContent(i);
    err += h1_phi->GetBinError(i);
  }
  Printf("bin start %i, bin end %i, sum %f ", bin_start, bin_end, sum);
  sum /= (bin_end - bin_start + 1);
  err /= (bin_end - bin_start + 1);
  Printf(" sum %f ", sum);  

  Printf("sum %f, bin width %f ", sum, h1_phi->GetBinWidth(bin_start));

  double background = 2*alpha_radius*sum;//the plot is already normalized for the bin_width by getting it from getSumOfRatios2
  double background_err = 2*alpha_radius*err;

  //TCanvas* can = new TCanvas();
  //h1_phi->DrawCopy();

  //get error from mixed combinatorics
  Int_t trigger_mixed_comb;




  //TH2* h2_etaPhi_mixedComb = (TH2D*) h->GetData()->GetSumOfRatios2(h->GetData(), step_mixedComb, 0, pt1Min, pt1Max, multBinBegin, multBinEnd, kFALSE, step_mixed, &trigger_mixed_comb);
  h2_etaPhi_mixedComb = getMixedComb_scaled_backgroundSame(fileName, pt1Min, pt1Max, pt2Min, ptAssocMin, ptAssocMax, multBinBegin, multBinEnd, setVertex, step_mixedComb, &trigger_mixed_comb);

  int firstbin_eta = h2_etaPhi_mixedComb->GetYaxis()->FindBin(-0.39);
  int lastbin_eta = h2_etaPhi_mixedComb->GetYaxis()->FindBin(0.39);
  TH1* h1_phi_mixedComb = h2_etaPhi_mixedComb->ProjectionX("mixedComb_px", firstbin_eta, lastbin_eta);

  //get error from mixed combinatorics
  Int_t trigger_backgroundSame;
  
  TH2* h2_etaPhi_backgroundSame = (TH2D*) h->GetData()->GetSumOfRatios2(h->GetData(), step_backgroundSame, 0, pt1Min, pt1Max, multBinBegin, multBinEnd, kFALSE, step_mixed, &trigger_backgroundSame);
  TH1* h1_phi_backgroundSame = h2_etaPhi_backgroundSame->ProjectionX("backgroundSame_px", firstbin_eta, lastbin_eta);

  
  TCanvas* can = new TCanvas("can", "can",  gBasisSize+50, gBasisSize+50, 2*gBasisSize, 2*gBasisSize);
  h1_phi_mixedComb->SetLineWidth(3);
  h1_phi_mixedComb->SetStats(kFALSE);
  h1_phi_mixedComb->DrawCopy();
  h1_phi_backgroundSame->SetLineWidth(3);
  h1_phi_backgroundSame->SetLineColor(kRed);
  h1_phi_backgroundSame->DrawCopy("same");

  TLegend *leg  = getLegend();
  leg->AddEntry(h1_phi_mixedComb, "mixed combinatorics","l");
  leg->AddEntry(h1_phi_backgroundSame, "background same","l");
  leg->Draw("same");
  
  Printf("expect from 1+1 %f +/- %f", background, background_err);
  Printf("expect from trigger_mixed_comb %i ", trigger_mixed_comb);
  Printf("expect from trigger_backgroundSame %i ", trigger_backgroundSame); 
}











void symmetrize(TH1D* tracks)
{
  int bins = tracks->FindBin(TMath::Pi()/2-0.01);
  for(int i=0; i<=bins/2; i++){
    double tmp =  (tracks->GetBinContent(i)+tracks->GetBinContent(bins-i+1))/2;
    tracks->SetBinContent(i, tmp);
    tracks->SetBinContent(bins-i+1, tmp);
  }
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
  
  //symmetrize
  if(symmetrize>0)
    for(i=0; i<=tracksMixed->GetNbinsX(); i++){
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

  tracksMixed->Scale(1.0/divide);

  TH1D* proj = tracksMixed->ProjectionY("_norm");
  Int_t bins_X = tracksMixed->GetNbinsX();

  double norm1 = proj->GetBinContent(proj->FindBin(-0.001));
  double norm2 = proj->GetBinContent(proj->FindBin(0.001));
  double norm = (norm1+norm2)/2.0;

  //tracksMixed->Scale(1.0/(bins_X*norm));
  Printf("norm1 %f, norm2 %f, norm %f", norm1/norm, norm2/norm, norm);

  //TCanvas* c2 = new TCanvas("mixedEvent", "mixedEvent", gBasisSize+50, 2*gBasisSize+50, gBasisSize, gBasisSize);
  //tracksMixed->DrawCopy("surf1");

  return tracksMixed;
}




//old method with an double cache
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

  return tracksMixed;
}


void showMixedDist(const char* fileName){
  loadlibs();

  list = (TList*) getList(fileName, "PWGCF_TwoPlusOne/histosTwoPlusOne");

  TH1F* mixedDist = (TH1F*) list->FindObject("mixedDist");
  TCanvas* can = new TCanvas();
  mixedDist->DrawCopy("colz");

}



void compareScaledMixedComb(const char* fileName, double pt1Min = 4.0, double pt1Max = 14.0, Double_t pt2Min = 2.0, Double_t ptAssocMin = 0.5, Double_t ptAssocMax = 8.0, Double_t setVertex = 7, Int_t multBinBegin = 1, Int_t multBinEnd = 5, Int_t step = 0){

  loadlibs();

  //to guarantee to pick only the bins from the choosen pt on
  pt1Min += 0.01;
  pt1Max -= 0.01;
  pt2Min += 0.01;
  ptAssocMin += 0.01;
  ptAssocMax -= 0.01;
  
  AliTwoPlusOneContainer* h = (AliTwoPlusOneContainer*) GetTwoPlusOne(fileName);

  AliUEHist::CFStep step_same = step;//(AliUEHist::CFStep) AliTwoPlusOneContainer::kSameNS;
  AliUEHist::CFStep step_mixed = step+2;//(AliUEHist::CFStep) AliTwoPlusOneContainer::kMixedNS;
  AliUEHist::CFStep step_mixedComb = step+4;//(AliUEHist::CStep) AliTwoPlusOneContainer::kMixedCombNS;
  AliUEHist::CFStep step_backgroundSame = step+7;//(AliUEHist::CStep) AliTwoPlusOneContainer::kBackgroundSameNS;

  h->GetData()->SetPtRange(ptAssocMin, ptAssocMax);
  h->GetData()->SetPt2Min(pt2Min);
  h->GetData()->SetZVtxRange(-1*setVertex, setVertex);
 
  //GetSumOfRatios2(mixed, step, region, ptLeadMin, ptLeadMax, multBinBegin, multBinEnd, normalizePerTrigger, stepForMixed)
  TH2D* h2_etaPhi;
  TH2D* h2_etaPhi_mixedComb;
  TH2D* h2_etaPhi_backgroundSame;


    Int_t trigger_same;
    Int_t trigger_mixed_comb;
    Int_t trigger_mixed_comb_method;
    Int_t trigger_background_same;


    h2_etaPhi_mixedComb_sumOfRatios = (TH2D*) h->GetData()->GetSumOfRatios2(h->GetData(), step_mixedComb, 0, pt1Min, pt1Max, multBinBegin, multBinEnd, kFALSE, step_mixed, &trigger_mixed_comb);
    h2_etaPhi_backgroundSame_sumOfRatios = (TH2D*) h->GetData()->GetSumOfRatios2(h->GetData(), step_backgroundSame, 0, pt1Min, pt1Max, multBinBegin, multBinEnd, kFALSE, step_mixed, &trigger_background_same);
    //trigger_background_same *= 2;//error in the plots for run 968 and 973
    double trigger_ratio = (double)trigger_background_same/trigger_mixed_comb;
    Printf("trigger mixed comb: %i, trigger background same %i, ratio %f", trigger_mixed_comb, trigger_background_same, trigger_ratio);
    h2_etaPhi_mixedComb_sumOfRatios->Scale(trigger_ratio);
    trigger_mixed_comb = trigger_background_same;
      
    h2_etaPhi_mixedComb_method = getMixedComb_scaled_backgroundSame(fileName, pt1Min, pt1Max, pt2Min, ptAssocMin, ptAssocMax, multBinBegin, multBinEnd, setVertex, step_mixedComb, &trigger_mixed_comb_method);

    
    TCanvas* can = new TCanvas("can", "can",  gBasisSize+50, gBasisSize+50, gBasisSize, gBasisSize);
    h2_etaPhi_mixedComb_sumOfRatios->DrawCopy("surf1");

    TCanvas* can2 = new TCanvas("can2", "can2",  gBasisSize+50, 2*(gBasisSize+50), gBasisSize, gBasisSize);
    h2_etaPhi_mixedComb_method->DrawCopy("surf1");

   TCanvas* can3 = new TCanvas("can3", "can3",  gBasisSize+50, 3*(gBasisSize+50), gBasisSize, gBasisSize);
   h2_etaPhi_mixedComb_sumOfRatios->Add(h2_etaPhi_mixedComb_method, -1);
   h2_etaPhi_mixedComb_sumOfRatios->Divide(h2_etaPhi_mixedComb_method);
   h2_etaPhi_mixedComb_sumOfRatios->GetYaxis()->SetRangeUser(-1.5, 1.5);
   h2_etaPhi_mixedComb_sumOfRatios->DrawCopy("surf1");

    Printf("mixed comb sumofratios %i, method %i", trigger_mixed_comb, h2_etaPhi_mixedComb_sumOfRatios);


}


void Plot_oneAxis(const char* fileName, int project)
{
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
