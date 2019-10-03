{
  TFile* ef = TFile::Open("../AnalysisResults.root");
  
  TList *list_Event[2];
  list_Event[0] = (TList*)ef->Get("list_kINT7_Event");
  
  TH1F* fHistAnalyzedEvents[2];
  fHistAnalyzedEvents[0] = (TH1F*)list_Event[0]->FindObject("fHistAnalyzedEvents");

  TH1F* fHistAnalyzed0PH0Events[2];
  fHistAnalyzed0PH0Events[0] = (TH1F*)list_Event[0]->FindObject("fHistAnalyzed0PH0Events");

  TList *list_Cluster[2];
  list_Cluster[0] = (TList*)ef->Get("list_kINT7_Cluster");

  TH1F* fHistClustEne[2];
  fHistClustEne[0] = (TH1F*)list_Cluster[0]->FindObject("fHistGoodTRUClustEneNcellCutTOFcut1");
  fHistClustEne[0]->Sumw2();
  
  TH1F* fHist0PH0ClustEne[2];
  fHist0PH0ClustEne[0] = (TH1F*)list_Cluster[0]->FindObject("fHistGoodTRU0PH0ClustEneNcellCutTOFcut1");
  fHist0PH0ClustEne[0]->Sumw2();
  
  Double_t bins[] = {0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.5,4.0,4.5,5.0,6.0,7.0,8.0,9.0,10,12,14,16,20,25,30,40};
  const Int_t binnum = sizeof(bins)/sizeof(Double_t) - 1;

  TH1F* fHistClustEne_Rebin[2];
  TH1F* fHist0PH0ClustEne_Rebin[2];
  fHistClustEne_Rebin[0]     = new TH1F("fHistClustEneINT7","",binnum,bins);
  fHist0PH0ClustEne_Rebin[0] = new TH1F("fHist0PH0ClustEneINT7","",binnum,bins);

  for(Int_t i=0; i<binnum; ++i){
    
    Int_t min_bin = fHistClustEne[0]->GetXaxis()->FindBin(bins[i]);
    Int_t max_bin = fHistClustEne[0]->GetXaxis()->FindBin(bins[i+1]);

    Double_t width = bins[i+1]-bins[i];

    Double_t nINT7   = 0;
    Double_t e_nINT7 = 0;
    Double_t nPHI7   = 0;
    Double_t e_nPHI7 = 0;

    Double_t n0PH0INT7   = 0;
    Double_t e_n0PH0INT7 = 0;
    Double_t n0PH0PHI7   = 0;
    Double_t e_n0PH0PHI7 = 0;

    for(Int_t j=min_bin; j<max_bin; ++j){

      nINT7   += fHistClustEne[0]->GetBinContent(j);
      e_nINT7 += pow(fHistClustEne[0]->GetBinError(j),2);

      n0PH0INT7   += fHist0PH0ClustEne[0]->GetBinContent(j);
      e_n0PH0INT7 += pow(fHist0PH0ClustEne[0]->GetBinError(j),2);
    }

    fHistClustEne_Rebin[0]->SetBinContent(i+1,nINT7/width);
    fHistClustEne_Rebin[0]->SetBinError(i+1,sqrt(e_nINT7)/width);
    fHist0PH0ClustEne_Rebin[0]->SetBinContent(i+1,n0PH0INT7/width);
    fHist0PH0ClustEne_Rebin[0]->SetBinError(i+1,sqrt(e_n0PH0INT7)/width);
  }
  
  TGraphAsymmErrors* graphTrigEff = new TGraphAsymmErrors();
  graphTrigEff->BayesDivide(fHist0PH0ClustEne_Rebin[0],fHistClustEne_Rebin[0]);
  graphTrigEff->SetName("graphTrigEff");
  
  TFile *tfout = new TFile("TrigEff.root","recreate");
  tfout->WriteTObject(graphTrigEff);
  
  
}
