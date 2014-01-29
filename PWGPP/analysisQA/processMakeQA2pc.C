const char * tag = "_Hadrons_QA_PbPb";
Float_t gpTMin = 0.51;
Float_t gpTMax = 49.99;
const char* lastFileName = 0;
void* cacheSameEvent = 0;
void* cacheMixedEvent = 0;

void processMakeQA2pc(const char * filename="AnalysisResults",  //file name without ".root"
		      TString suffix = "eps",
		      const char *outfile="MakeQA2pc_output.root"){
  
  loadlibs();
  if (!gGrid && TString(filename).BeginsWith("alien://"))
    TGrid::Connect("alien://");
  
  //fill CFContainers
  FillParentTHnSparse(Form("%s.root",filename),kFALSE,tag);
  
  Double_t ptTmin[]={0.5,0.5,2.0,2.0};
  Double_t ptTmax[]={1.0,1.0,4.0,4.0};
  Double_t ptAmin[]={0.5,0.5,2.0,2.0};
  Double_t ptAmax[]={1.0,1.0,4.0,4.0};
  Int_t   centBegin[]={  0, 60,  0, 60};
  Int_t     centEnd[]={ 10, 90, 10, 90};
  


  for(Int_t i=0;i<4;i++){
    TCanvas *c=new TCanvas(Form("c%d",i),Form("c%d",i),1200,1200);
    DrawSameMixedSV(Form("%s_zvtx.root",filename),ptTmin[i],ptTmax[i],ptAmin[i],ptAmax[i],centBegin[i],centEnd[i],8,c,i,outfile);
    c->SaveAs(Form("fig_cf_c%d.%s",i,suffix.Data()));
  }
  
  TCanvas *ccorr=new TCanvas("ccorr","ccorr",1200,1200);
  DrawCorrelation(Form("%s_zvtx.root",filename),ccorr,outfile);


  ccorr->SaveAs(Form("fig_cf_ccorr.%s",suffix.Data()));
  
}

void SetupRanges(void* obj)
{
  if (!obj)
    return;
  ((AliUEHistograms*) obj)->SetEtaRange(0, 0);
  ((AliUEHistograms*) obj)->SetPtRange(gpTMin, gpTMax);
  ((AliUEHistograms*) obj)->SetCombineMinMax(kTRUE);
}

void DrawCorrelation(const char* fileName,TCanvas *c,const char *outfile){
  
  AliUEHistograms* h = (AliUEHistograms*) GetUEHistogram(fileName);
  c->Divide(1,3);
  c->cd(1);
  gPad->SetLogz();
  h->((TH2F*)h->GetCorrelationpT())->DrawCopy("colz")->GetYaxis()->SetRangeUser(0.,20.);
  c->cd(2);
  h->((TH2F*)h->GetCorrelationEta())->DrawCopy("colz")->GetYaxis()->SetRangeUser(-1.5,1.5);
  c->cd(3);
  h->((TH2F*)h->GetCorrelationPhi())->DrawCopy("colz");

  //Added by sjena
  TFile *fout = TFile::Open(outfile,"UPDATE");
  fout->ls();

  TDirectoryFile *cdd = NULL;
  cdd = (TDirectoryFile*)fout->Get("CF");
  if(!cdd) {
    Printf("Warning: CF <dir> doesn't exist, creating a new one");
    cdd = (TDirectoryFile*)fout->mkdir("CF");
  }
  cdd->cd();
  cdd->ls();
  

  TH2F *h1 = (TH2F*)h->((TH2F*)h->GetCorrelationpT())->Clone();
  TH2F *h2 = (TH2F*)h->((TH2F*)h->GetCorrelationEta())->Clone();
  TH2F *h3 = (TH2F*)h->((TH2F*)h->GetCorrelationPhi())->Clone();
  h1->Write();
  h2->Write();
  h3->Write();
  fout->Close();

}


void DrawSameMixed(const char* fileName,Double_t ptTmin,Double_t ptTmax,Double_t ptAmin,Double_t ptAmax,Int_t centBegin,Int_t centEnd,Int_t step = 8,TCanvas *c)
{
  c->Divide(3, 3);
  TPaveText* paveText = new TPaveText(0.2, 0.9, 1., 1., "BRNDC");
  paveText->SetTextSize(0.04);
  paveText->SetFillColor(0);
  paveText->SetShadowColor(0);
  paveText->SetBorderSize(0);
  paveText->SetFillStyle(0);
  paveText->AddText(Form("%.1f<p_{T,trig}<%.1f, %.1f<p_{T,ass}<%.1f, %d-%d %%",ptTmin,ptTmax,ptAmin,ptAmax,centBegin,centEnd));
  
  loadlibs();
  
  AliUEHistograms* h = (AliUEHistograms*) GetUEHistogram(fileName);
  hMixed = (AliUEHistograms*) GetUEHistogram(fileName, 0, kTRUE);
  
  gpTMin=ptAmin;
  gpTMax=ptAmax;
  
  SetupRanges(h);
  SetupRanges(hMixed);
  
  TH1* hist1 = 0;
  TH1* hist2 = 0;
  GetDistAndFlow(h, &hist1, step, centBegin,  centEnd, ptTmin, ptTmax); 
  
  NormalizeToBinWidth(hist1);
  
  c->cd(1);
  gPad->SetLeftMargin(0.15);
  hist1->SetTitle("");
  hist1->GetYaxis()->SetRangeUser(-1.79, 1.79);
  hist1->GetZaxis()->SetTitleOffset(1.8);
  hist1->GetXaxis()->SetTitleOffset(1.5);
  hist1->GetYaxis()->SetTitleOffset(2);
  hist1->GetZaxis()->SetTitle("same event pairs (a.u.)");
  hist1->SetStats(kFALSE);
  hist1->DrawCopy("SURF1");
  paveText->Draw();
  c->cd(4);
  ((TH2*) hist1)->ProjectionX()->DrawCopy();
  c->cd(7);
  ((TH2*) hist1)->ProjectionY()->DrawCopy();
  
  hist2 = hist1;
  
  GetDistAndFlow(hMixed, &hist1, step, centBegin,  centEnd, ptTmin, ptTmax); 
  NormalizeToBinWidth(hist1);
  
  c->cd(2);
  gPad->SetLeftMargin(0.15);
  hist1->SetTitle("");
  hist1->GetYaxis()->SetRangeUser(-1.79, 1.79);
  hist1->GetZaxis()->SetTitleOffset(1.8);
  hist1->GetXaxis()->SetTitleOffset(1.5);
  hist1->GetYaxis()->SetTitleOffset(2);
  hist1->GetZaxis()->SetTitle("mixed event pairs (a.u.)");
  hist1->SetStats(kFALSE);
  hist1->DrawCopy("SURF1");
  c->cd(5);
  ((TH2*) hist1)->ProjectionX()->DrawCopy();
  c->cd(8);
  ((TH2*) hist1)->ProjectionY()->DrawCopy();
  
  c->cd(3);
  gPad->SetLeftMargin(0.15);
  hist2->GetZaxis()->SetTitle("same/mixed event pairs (a.u.)");
  hist2->Divide(hist1);
  hist2->DrawCopy("SURF1");
  c->cd(6);
  ((TH2*) hist2)->ProjectionX()->DrawCopy();
  c->cd(9);
  ((TH2*) hist2)->ProjectionY()->DrawCopy();
  
}


void DrawSameMixedSV(const char* fileName,Double_t ptTmin,Double_t ptTmax,Double_t ptAmin,Double_t ptAmax,Int_t centBegin,Int_t centEnd,Int_t step = 8,TCanvas *c, Int_t key, const char *outfile)
{
  c->Divide(3, 3);
  TPaveText* paveText = new TPaveText(0.2, 0.9, 1., 1., "BRNDC");
  paveText->SetTextSize(0.04);
  paveText->SetFillColor(0);
  paveText->SetShadowColor(0);
  paveText->SetBorderSize(0);
  paveText->SetFillStyle(0);
  paveText->AddText(Form("%.1f<p_{T,trig}<%.1f, %.1f<p_{T,ass}<%.1f, %d-%d %%",ptTmin,ptTmax,ptAmin,ptAmax,centBegin,centEnd));
  
  loadlibs();
  
  AliUEHistograms* h = (AliUEHistograms*) GetUEHistogram(fileName);
  hMixed = (AliUEHistograms*) GetUEHistogram(fileName, 0, kTRUE);
  
  gpTMin=ptAmin;
  gpTMax=ptAmax;
  
  SetupRanges(h);
  SetupRanges(hMixed);
  
  TH1* hist1 = 0;
  TH1* hist2 = 0;
  GetDistAndFlow(h, &hist1, step, centBegin,  centEnd, ptTmin, ptTmax); 
  
  NormalizeToBinWidth(hist1);
  
  c->cd(1);
  gPad->SetLeftMargin(0.15);
  hist1->SetTitle("");
  hist1->GetYaxis()->SetRangeUser(-1.79, 1.79);
  hist1->GetZaxis()->SetTitleOffset(1.8);
  hist1->GetXaxis()->SetTitleOffset(1.5);
  hist1->GetYaxis()->SetTitleOffset(2);
  hist1->GetZaxis()->SetTitle("same event pairs (a.u.)");
  hist1->SetStats(kFALSE);
  hist1->DrawCopy("SURF1");
  
  //added by sjena
  TFile *fout = TFile::Open(outfile,"UPDATE");
  fout->ls();
  
  TDirectoryFile *cdd = NULL;
  cdd = (TDirectoryFile*)fout->Get("CF");
  if(!cdd) {
    Printf("Warning: CF <dir> doesn't exist, creating a new one");
    cdd = (TDirectoryFile*)fout->mkdir("CF");
  }
  cdd->cd();
  cdd->ls();
  

  TH1D *h1 = (TH1D*)hist1->Clone();
  h1->SetName(Form("fig_cf_cr1_%d",key));
  h1->SetTitle(hist1->GetName());
  h1->Write();
  

  paveText->Draw();
  c->cd(4);
  ((TH2*) hist1)->ProjectionX()->DrawCopy();
 
  TH2D *h2 =  (TH2D*)((TH2*) hist1)->ProjectionX();
  h2->SetTitle(h2->GetName());
  h2->SetName(Form("fig_cf_cr2_%d", key));
  h2->Write();

  c->cd(7);
  ((TH2*) hist1)->ProjectionY()->DrawCopy();
  //ofile->cd();
  TH2D *h3 =  (TH2D*)((TH2*) hist1)->ProjectionY();
  h3->SetTitle(h3->GetName());
  h3->SetName(Form("fig_cf_cr3_%d",key));
  h3->Write();




  hist2 = hist1;
  
  GetDistAndFlow(hMixed, &hist1, step, centBegin,  centEnd, ptTmin, ptTmax); 
  NormalizeToBinWidth(hist1);
  
  c->cd(2);
  gPad->SetLeftMargin(0.15);
  hist1->SetTitle("");
  hist1->GetYaxis()->SetRangeUser(-1.79, 1.79);
  hist1->GetZaxis()->SetTitleOffset(1.8);
  hist1->GetXaxis()->SetTitleOffset(1.5);
  hist1->GetYaxis()->SetTitleOffset(2);
  hist1->GetZaxis()->SetTitle("mixed event pairs (a.u.)");
  hist1->SetStats(kFALSE);
  hist1->DrawCopy("SURF1");
  //ofile->cd();
   TH1D *h4 =  (TH1D*)hist1->Clone();
  h4->SetTitle(hist1->GetName());
  h4->SetName(Form("fig_cf_cr4_%d",key));
  h4->Write();





  c->cd(5);
  ((TH2*) hist1)->ProjectionX()->DrawCopy();
  //ofile->cd();
  TH2D *h5 =  (TH2D*)((TH2*) hist1)->ProjectionX();
  h5->SetTitle(h5->GetName());
  h5->SetName(Form("fig_cf_cr5_%d",key));
  h5->Write();




  c->cd(8);
  ((TH2*) hist1)->ProjectionY()->DrawCopy();
  //ofile->cd();
  TH2D *h6 =  (TH2D*)((TH2*) hist1)->ProjectionY();
  h6->SetTitle(h6->GetName());
   h6->SetName(Form("fig_cf_cr6_%d",key));
   h6->Write();  


  c->cd(3);
  gPad->SetLeftMargin(0.15);
  hist2->GetZaxis()->SetTitle("same/mixed event pairs (a.u.)");
  hist2->Divide(hist1);
  hist2->DrawCopy("SURF1");
  //ofile->cd();
 
  TH1D *h7 =  (TH1D*)hist2->Clone();
  h7->SetTitle(hist2->GetName());
  h7->SetName(Form("fig_cf_cr7_%d",key));
  h7->Write();


  c->cd(6);
  ((TH2*) hist2)->ProjectionX()->DrawCopy();
  //ofile->cd();
  TH2D *h8 =  (TH2D*)((TH2*) hist2)->ProjectionX();
  h8->SetTitle(h8->GetName());
  h8->SetName(Form("fig_cf_cr8_%d",key));
  h8->Write();  
  

  c->cd(9);
  ((TH2*) hist2)->ProjectionY()->DrawCopy();
  //ofile->cd();
  TH2D *h9 =  (TH2D*)((TH2*) hist2)->ProjectionY();
  h9->SetTitle(h9->GetName());
  h9->SetName(Form("fig_cf_cr9_%d",key));
  h9->Write();  

  fout->Close();
}


void loadlibs()
{
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCORRFW");
  gSystem->Load("libPWGTools");
  gSystem->Load("libPWGCFCorrelationsBase");
}

void FillParentTHnSparse(const char* fileName, Bool_t reduce = kFALSE, const char* tag = "")
{
  if (TString(fileName).BeginsWith("alien:"))
    TGrid::Connect("alien:");
  
  loadlibs();

  TList* list = 0;
  
  AliUEHistograms* h = (AliUEHistograms*) GetUEHistogram(fileName, &list, kFALSE, tag);

  
  Printf("We have %d axes", ((AliTHn*) h->GetUEHist(2)->GetTrackHist(0)->GetNVar()));
  
  if (reduce)
    ((AliTHn*) h->GetUEHist(2)->GetTrackHist(0))->ReduceAxis();
  ((AliTHn*) h->GetUEHist(2)->GetTrackHist(0))->FillParent();
  ((AliTHn*) h->GetUEHist(2)->GetTrackHist(0))->DeleteContainers();
  
  AliUEHistograms* hMixed = (AliUEHistograms*) GetUEHistogram(fileName, 0, kTRUE, tag);
  if (reduce)
    ((AliTHn*) hMixed->GetUEHist(2)->GetTrackHist(0))->ReduceAxis();
  ((AliTHn*) hMixed->GetUEHist(2)->GetTrackHist(0))->FillParent();
  ((AliTHn*) hMixed->GetUEHist(2)->GetTrackHist(0))->DeleteContainers();
  
  TString newFileName(fileName);

  if (TString(fileName).BeginsWith("alien:"))
    newFileName = gSystem->BaseName(newFileName);
  
  newFileName.ReplaceAll(".root", "");
  if (reduce)
    newFileName += "_.root";
  else
    newFileName += "_zvtx.root";

  file3 = TFile::Open(newFileName, "RECREATE");
  file3->mkdir("PWG4_PhiCorrelations");
  file3->cd("PWG4_PhiCorrelations");
  list->Write("histosPhiCorrelations", TObject::kSingleKey);
  file3->Close();
}

void* GetUEHistogram(const char* fileName, TList** listRef = 0, Bool_t mixed = kFALSE, const char* tag = "")
{
  if (!lastFileName || strcmp(lastFileName, fileName) != 0)
    {
      lastFileName = fileName;
      file = TFile::Open(fileName);
      if (!file)
	return 0;
      
      list = (TList*) gFile->Get("PWG4_LeadingTrackUE/histosLeadingTrackUE");
      if (!list)
	list = (TList*) gFile->Get(Form("PWG4_PhiCorrelations/histosPhiCorrelations%s", tag));
      if (!list)
	list = (TList*) gFile->Get("PWG4_PhiCorrelations/histosPhiCorrelations_Syst");
      
      if (!list)
	return 0;
      
      if (listRef)
	*listRef = list;
      
      cacheMixedEvent = list->FindObject("AliUEHistogramsMixed");
      cacheSameEvent = list->FindObject("AliUEHistogramsSame");

      if (mixed)
	return cacheMixedEvent;
    
      if (list->FindObject("AliUEHistograms"))
	return list->FindObject("AliUEHistograms");
      
      return cacheSameEvent;
    }
  else
    {
      Printf("GetUEHistogram --> Using cache for %s", fileName);
    
      if (mixed)
	return cacheMixedEvent;
      else
	return cacheSameEvent;
    }
}


Int_t gHistCount = 0;

void GetDistAndFlow(void* hVoid, TH1** hist, Int_t step, Int_t centralityBegin, Int_t centralityEnd, Float_t ptBegin, Float_t ptEnd)
{
  h = (AliUEHistograms*) hVoid;
  
  Int_t centralityBeginBin = 0;
  Int_t centralityEndBin = -1;
  
  if (centralityEnd >= centralityBegin)
    {
      centralityBeginBin = h->GetUEHist(2)->GetEventHist()->GetGrid(step)->GetGrid()->GetAxis(1)->FindBin(0.01 + centralityBegin);
      centralityEndBin = h->GetUEHist(2)->GetEventHist()->GetGrid(step)->GetGrid()->GetAxis(1)->FindBin(-0.01 + centralityEnd);
    }
  
  // 2d same and mixed event
  TH2* sameTwoD  = h->GetUEHist(2)->GetUEHist(step, 0, ptBegin, ptEnd, centralityBeginBin, centralityEndBin, 1, kFALSE);
  
  TString histName;
  histName.Form("GetDistAndFlow%d", gHistCount++);
  
  *hist = sameTwoD;
    
  TString str;
  str.Form("%.1f < p_{T,trig} < %.1f", ptBegin - 0.01, ptEnd + 0.01);
  
  TString str2;
  str2.Form("%.2f < p_{T,assoc} < %.2f", gpTMin - 0.01, gpTMax + 0.01);
  
  TString newTitle;
  newTitle.Form("%s - %s - %d-%d%%", str.Data(), str2.Data(), centralityBegin, centralityEnd);
  (*hist)->SetTitle(newTitle);
  
}

void NormalizeToBinWidth(TH1* hist)
{
  //
  // normalizes a 1-d histogram to its bin width
  //

  if (hist->GetDimension() == 1)
    {
      for (Int_t i=1; i<=hist->GetNbinsX(); ++i)
	{
	  hist->SetBinContent(i, hist->GetBinContent(i) / hist->GetBinWidth(i));
	  hist->SetBinError(i, hist->GetBinError(i) / hist->GetBinWidth(i));
	}
    }
  else if (hist->GetDimension() == 2)
    {
      for (Int_t i=1; i<=hist->GetNbinsX(); ++i)
	{
	  for (Int_t j=1; j<=hist->GetNbinsY(); ++j)
	    {
	      hist->SetBinContent(i, j, hist->GetBinContent(i, j) / hist->GetXaxis()->GetBinWidth(i) / hist->GetYaxis()->GetBinWidth(j));
	      hist->SetBinError(i, j, hist->GetBinError(i) / hist->GetXaxis()->GetBinWidth(i) / hist->GetYaxis()->GetBinWidth(j));
	    }
	}
    }
}
