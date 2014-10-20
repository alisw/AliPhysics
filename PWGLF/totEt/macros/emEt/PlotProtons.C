TH1* bayneseffdiv(TH1* numerator, TH1* denominator,Char_t* name);
void SetStyles(TH1 *histo,int marker, int color){
  histo->Sumw2();
  histo->SetMarkerStyle(marker);
  histo->SetMarkerColor(color);
  histo->SetLineColor(color);
  //histo->GetXaxis()->SetTitle(xtitle);
  //histo->GetYaxis()->SetTitle(ytitle);
}
void Rescale(TH1 *histo,Int_t rescale){
  histo->Rebin(rescale);
  histo->Scale(1.0/((Float_t)rescale));
}
void SetStyles(TH1 *histo,int marker, int color,char *name){
  SetStyles(histo,marker,color);
  histo->SetName(name);
}
void PrintInfo(TH1 *histo){
  cout<<histo->GetName()<<":"<<endl;
  cout<<"x range "<<histo->GetBinLowEdge(1)<<" - "<<histo->GetBinLowEdge(histo->GetNbinsX()+1)<<endl;
  cout<<"y range "<<histo->GetMinimum()<<" - "<<histo->GetMaximum()<<endl;
}
Int_t colors[] = {0,TColor::kRed, TColor::kOrange, TColor::kGreen+3, TColor::kBlue, TColor::kBlack, 
		    TColor::kRed, TColor::kOrange, TColor::kGreen+3, TColor::kBlue, TColor::kBlack, 
		    TColor::kRed, TColor::kOrange, TColor::kGreen+3, TColor::kBlue, TColor::kBlack, 
		    TColor::kRed, TColor::kOrange, TColor::kGreen+3, TColor::kBlue, TColor::kBlack};
Int_t markers[] = {20,21,22,23,33, 24,25,26,32,27, 20,21,22,23,33, 24,25,26,32,27};


Float_t nNeutronsSim[] = {0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0};
Float_t nNeutronsData[] = {0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0};
Float_t nNeutronsDataErr[] = {0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0};
Float_t nNeutronsShortSim[] = {0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0};
Float_t nNeutronsShortData[] = {0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0};
Float_t nNeutronsShortDataErr[] = {0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0};
Float_t nNeutronsSimCl[] = {0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0};
Float_t nNeutronsDataCl[] = {0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0};
Float_t nNeutronsShortSimCl[] = {0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0};
Float_t nNeutronsShortDataCl[] = {0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0};

void WriteLatex();
Float_t neutronCorrEmcal[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t neutronCorrPhos[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t neutronErrorEmcal[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t neutronErrorPhos[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
TH1D* Divide(TH1D* numerator, TH1D* denominator,Char_t* name){
  //TH1D *output = numerator->Clone(name);
  TH1D *output = new TH1D(name,name,numerator->GetNbinsX(),numerator->GetBinLowEdge(1),numerator->GetBinLowEdge(numerator->GetNbinsX()+1));
  output->GetXaxis()->SetTitle(numerator->GetXaxis()->GetTitle());
  output->GetYaxis()->SetTitle(numerator->GetYaxis()->GetTitle());
  Int_t nbinsNum = numerator->GetNbinsX();
  Int_t nbinsDen = denominator->GetNbinsX();
  for(Int_t i=1;i<=output->GetNbinsX();i++){
    if(nbinsNum!=nbinsDen){
      if(denominator->GetBinContent(i)>0){
	output->SetBinContent(i,numerator->GetBinContent(i)/denominator->GetBinContent(i));
      }
    }
    else{
      output->SetBinContent(i,numerator->GetBinContent(i));
    }
    //output->SetBinError((Int_t)i,0,0,(Double_t)0.0);
    //denominator->SetBinError((Int_t)i,0,0,(Double_t)0.0);
    //cout<<" "<<output->GetBinError(i)<<" "<<denominator->GetBinError(i)<<endl;
  }
  if(nbinsNum==nbinsDen){
    output->Divide(denominator);
  }
  //cout<<"Fixing "<<output->GetName()<<endl;
  Double_t newerror = 1e-5;
  for(Int_t i=1;i<=output->GetNbinsX();i++){
    //cout<<"Bin "<<i<<" setting bin error "<<output->GetBinError(i)<<" to 0.  New bin error: ";
    output->SetBinError((Int_t)i,(Double_t)newerror);
    //output->SetBinError((Int_t)i,0,0,(Double_t)0.0);
    //cout<<" "<<output->GetBinError(i)<<endl;
  }
  return output;
}

void PlotProtons(TString filename = "rootFiles/LHC11a10a_bis/Et.ESD.simPbPb.EMCal.LHC11a10a_bis.Run139465.root", TString filenameData = "rootFiles/LHC10hPass2/Et.ESD.realPbPb.EMCal.LHC10hPass2.Run139465.root", Bool_t effCorr = kTRUE){
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  Bool_t isPhos = kTRUE;
  TString detector = "Phos";
  if(filename.Contains("EMC")){
    detector = "Emcal";
    isPhos = kFALSE;
  }

  TString tag = "";
  if(!effCorr) tag = "NoEffCorr";

  ofstream myfile;
  TString textfilename = "Neutrons"+detector+tag+".dat";
  myfile.open (textfilename.Data());
  ofstream myfile2;
  TString textfilename2 = "Neutrons"+detector+tag+"Short.dat";
  myfile2.open (textfilename2.Data());

  //Reading histograms
  TFile *f = TFile::Open(filename, "READ");
  TList *l = dynamic_cast<TList*>(f->Get("out1"));
  TH2F *fHistPIDProtonsTrackMatchedDepositedVsNch;
  TH2F *fHistPIDAntiProtonsTrackMatchedDepositedVsNch;
  TH2F *fHistPiKPTrackMatchedDepositedVsNch;
  TH2F *fHistNeutronsDepositedVsNch;
  TH2F *fHistAntiNeutronsDepositedVsNch;
  TH2F *fHistProtonsDepositedVsNch;
  TH2F *fHistAntiProtonsDepositedVsNch;
  TH2F *fHistPiKPDepositedVsNch;
  TH2F *fHistProtonsNotTrackMatchedDepositedVsNch;
  TH2F *fHistAntiProtonsNotTrackMatchedDepositedVsNch;
  TH2F *fHistPIDProtonsTrackMatchedDepositedVsNcl;
  TH2F *fHistPIDAntiProtonsTrackMatchedDepositedVsNcl;
  TH2F *fHistPiKPTrackMatchedDepositedVsNcl;
  TH2F *fHistNeutronsDepositedVsNcl;
  TH2F *fHistAntiNeutronsDepositedVsNcl;
  TH2F *fHistProtonsDepositedVsNcl;
  TH2F *fHistAntiProtonsDepositedVsNcl;
  TH2F *fHistPiKPDepositedVsNcl;
  TH2F *fHistProtonsNotTrackMatchedDepositedVsNcl;
  TH2F *fHistAntiProtonsNotTrackMatchedDepositedVsNcl;

  if(effCorr){
    fHistPIDProtonsTrackMatchedDepositedVsNch = (TH2F *)l->FindObject("fHistPIDProtonsTrackMatchedDepositedVsNch");
    fHistPIDAntiProtonsTrackMatchedDepositedVsNch = (TH2F *)l->FindObject("fHistPIDAntiProtonsTrackMatchedDepositedVsNch");
    fHistPiKPTrackMatchedDepositedVsNch = (TH2F *)l->FindObject("fHistPiKPTrackMatchedDepositedVsNch");
    fHistNeutronsDepositedVsNch = (TH2F *)l->FindObject("fHistNeutronsDepositedVsNch");
    fHistAntiNeutronsDepositedVsNch = (TH2F *)l->FindObject("fHistAntiNeutronsDepositedVsNch");
    fHistProtonsDepositedVsNch = (TH2F *)l->FindObject("fHistProtonsDepositedVsNch");
    fHistAntiProtonsDepositedVsNch = (TH2F *)l->FindObject("fHistAntiProtonsDepositedVsNch");
    fHistPiKPDepositedVsNch = (TH2F *)l->FindObject("fHistPiKPDepositedVsNch");
    fHistProtonsNotTrackMatchedDepositedVsNch = (TH2F *)l->FindObject("fHistProtonsNotTrackMatchedDepositedVsNch");
    fHistAntiProtonsNotTrackMatchedDepositedVsNch = (TH2F *)l->FindObject("fHistAntiProtonsNotTrackMatchedDepositedVsNch");

    fHistPIDProtonsTrackMatchedDepositedVsNcl = (TH2F *)l->FindObject("fHistPIDProtonsTrackMatchedDepositedVsNcl");
    fHistPIDAntiProtonsTrackMatchedDepositedVsNcl = (TH2F *)l->FindObject("fHistPIDAntiProtonsTrackMatchedDepositedVsNcl");
    fHistPiKPTrackMatchedDepositedVsNcl = (TH2F *)l->FindObject("fHistPiKPTrackMatchedDepositedVsNcl");
    fHistNeutronsDepositedVsNcl = (TH2F *)l->FindObject("fHistNeutronsDepositedVsNcl");
    fHistAntiNeutronsDepositedVsNcl = (TH2F *)l->FindObject("fHistAntiNeutronsDepositedVsNcl");
    fHistProtonsDepositedVsNcl = (TH2F *)l->FindObject("fHistProtonsDepositedVsNcl");
    fHistAntiProtonsDepositedVsNcl = (TH2F *)l->FindObject("fHistAntiProtonsDepositedVsNcl");
    fHistPiKPDepositedVsNcl = (TH2F *)l->FindObject("fHistPiKPDepositedVsNcl");
    fHistProtonsNotTrackMatchedDepositedVsNcl = (TH2F *)l->FindObject("fHistProtonsNotTrackMatchedDepositedVsNcl");
    fHistAntiProtonsNotTrackMatchedDepositedVsNcl = (TH2F *)l->FindObject("fHistAntiProtonsNotTrackMatchedDepositedVsNcl");
  }
  else{
    fHistPIDProtonsTrackMatchedDepositedVsNch = (TH2F *)l->FindObject("fHistPIDProtonsTrackMatchedDepositedVsNchNoEff");
    fHistPIDAntiProtonsTrackMatchedDepositedVsNch = (TH2F *)l->FindObject("fHistPIDAntiProtonsTrackMatchedDepositedVsNchNoEff");
    fHistPiKPTrackMatchedDepositedVsNch = (TH2F *)l->FindObject("fHistPiKPTrackMatchedDepositedVsNchNoEff");
    fHistNeutronsDepositedVsNch = (TH2F *)l->FindObject("fHistNeutronsDepositedVsNchNoEffCorr");
    fHistAntiNeutronsDepositedVsNch = (TH2F *)l->FindObject("fHistAntiNeutronsDepositedVsNchNoEffCorr");
    fHistProtonsDepositedVsNch = (TH2F *)l->FindObject("fHistProtonsDepositedVsNchNoEffCorr");
    fHistAntiProtonsDepositedVsNch = (TH2F *)l->FindObject("fHistAntiProtonsDepositedVsNchNoEffCorr");
    fHistPiKPDepositedVsNch = (TH2F *)l->FindObject("fHistPiKPDepositedVsNchNoEffCorr");
    fHistProtonsNotTrackMatchedDepositedVsNch = (TH2F *)l->FindObject("fHistProtonsNotTrackMatchedDepositedVsNchNoEffCorr");
    fHistAntiProtonsNotTrackMatchedDepositedVsNch = (TH2F *)l->FindObject("fHistAntiProtonsNotTrackMatchedDepositedVsNchNoEffCorr");



    fHistPIDProtonsTrackMatchedDepositedVsNcl = (TH2F *)l->FindObject("fHistPIDProtonsTrackMatchedDepositedVsNclNoEff");
    fHistPIDAntiProtonsTrackMatchedDepositedVsNcl = (TH2F *)l->FindObject("fHistPIDAntiProtonsTrackMatchedDepositedVsNclNoEff");
    fHistPiKPTrackMatchedDepositedVsNcl = (TH2F *)l->FindObject("fHistPiKPTrackMatchedDepositedVsNclNoEff");
    fHistNeutronsDepositedVsNcl = (TH2F *)l->FindObject("fHistNeutronsDepositedVsNclNoEffCorr");
    fHistAntiNeutronsDepositedVsNcl = (TH2F *)l->FindObject("fHistAntiNeutronsDepositedVsNclNoEffCorr");
    fHistProtonsDepositedVsNcl = (TH2F *)l->FindObject("fHistProtonsDepositedVsNclNoEffCorr");
    fHistAntiProtonsDepositedVsNcl = (TH2F *)l->FindObject("fHistAntiProtonsDepositedVsNclNoEffCorr");
    fHistPiKPDepositedVsNcl = (TH2F *)l->FindObject("fHistPiKPDepositedVsNclNoEffCorr");
    fHistProtonsNotTrackMatchedDepositedVsNcl = (TH2F *)l->FindObject("fHistProtonsNotTrackMatchedDepositedVsNclNoEffCorr");
    fHistAntiProtonsNotTrackMatchedDepositedVsNcl = (TH2F *)l->FindObject("fHistAntiProtonsNotTrackMatchedDepositedVsNclNoEffCorr");
  }

  TH3F *fHistCentVsNchVsNcl = l->FindObject("fHistCentVsNchVsNclReco");
  fHistCentVsNchVsNcl->GetXaxis()->SetTitle("cent");
  fHistCentVsNchVsNcl->GetYaxis()->SetTitle("N_{Ch}");
  fHistCentVsNchVsNcl->GetZaxis()->SetTitle("N_{Cl}");
  //fHistCentVsNchVsNcl->SetName("fHistCentVsNchVsNclSim");

  fHistPIDProtonsTrackMatchedDepositedVsNch->SetName(Form("%sSim","fHistPIDProtonsTrackMatchedDepositedVsNch"));
  fHistPIDAntiProtonsTrackMatchedDepositedVsNch->SetName(Form("%sSim","fHistPIDAntiProtonsTrackMatchedDepositedVsNch"));
  fHistPiKPTrackMatchedDepositedVsNch->SetName(Form("%sSim","fHistPiKPTrackMatchedDepositedVsNch"));

  TFile *fData = TFile::Open(filenameData, "READ");
  TList *lData = dynamic_cast<TList*>(fData->Get("out1"));
  TH2F *fHistPIDProtonsTrackMatchedDepositedVsNchData;
  TH2F *fHistPIDAntiProtonsTrackMatchedDepositedVsNchData;
  TH2F *fHistPiKPTrackMatchedDepositedVsNchData;
  TH2F *fHistPIDProtonsTrackMatchedDepositedVsNclData;
  TH2F *fHistPIDAntiProtonsTrackMatchedDepositedVsNclData;
  TH2F *fHistPiKPTrackMatchedDepositedVsNclData;
  if(effCorr){
    fHistPIDProtonsTrackMatchedDepositedVsNchData =(TH2F *) lData->FindObject("fHistPIDProtonsTrackMatchedDepositedVsNch");
    fHistPIDAntiProtonsTrackMatchedDepositedVsNchData = (TH2F *)lData->FindObject("fHistPIDAntiProtonsTrackMatchedDepositedVsNch");
    fHistPiKPTrackMatchedDepositedVsNchData = (TH2F *)lData->FindObject("fHistPiKPTrackMatchedDepositedVsNch");
    fHistPIDProtonsTrackMatchedDepositedVsNclData =(TH2F *) lData->FindObject("fHistPIDProtonsTrackMatchedDepositedVsNcl");
    fHistPIDAntiProtonsTrackMatchedDepositedVsNclData = (TH2F *)lData->FindObject("fHistPIDAntiProtonsTrackMatchedDepositedVsNcl");
    fHistPiKPTrackMatchedDepositedVsNclData = (TH2F *)lData->FindObject("fHistPiKPTrackMatchedDepositedVsNcl");
  }
  else{
    fHistPIDProtonsTrackMatchedDepositedVsNchData =(TH2F *) lData->FindObject("fHistPIDProtonsTrackMatchedDepositedVsNchNoEff");
    fHistPIDAntiProtonsTrackMatchedDepositedVsNchData = (TH2F *)lData->FindObject("fHistPIDAntiProtonsTrackMatchedDepositedVsNchNoEff");
    fHistPiKPTrackMatchedDepositedVsNchData = (TH2F *)lData->FindObject("fHistPiKPTrackMatchedDepositedVsNchNoEff");
    fHistPIDProtonsTrackMatchedDepositedVsNclData =(TH2F *) lData->FindObject("fHistPIDProtonsTrackMatchedDepositedVsNclNoEff");
    fHistPIDAntiProtonsTrackMatchedDepositedVsNclData = (TH2F *)lData->FindObject("fHistPIDAntiProtonsTrackMatchedDepositedVsNclNoEff");
    fHistPiKPTrackMatchedDepositedVsNclData = (TH2F *)lData->FindObject("fHistPiKPTrackMatchedDepositedVsNclNoEff");
  }
  TH3F *fHistCentVsNchVsNclData = lData->FindObject("fHistCentVsNchVsNclReco");
  fHistPIDProtonsTrackMatchedDepositedVsNchData->SetName(Form("%sData","fHistPIDProtonsTrackMatchedDepositedVsNch"));
  fHistPIDAntiProtonsTrackMatchedDepositedVsNchData->SetName(Form("%sData","fHistPIDAntiProtonsTrackMatchedDepositedVsNch"));
  fHistPiKPTrackMatchedDepositedVsNchData->SetName(Form("%sData","fHistPiKPTrackMatchedDepositedVsNch"));
  fHistPIDProtonsTrackMatchedDepositedVsNclData->SetName(Form("%sData","fHistPIDProtonsTrackMatchedDepositedVsNcl"));
  fHistPIDAntiProtonsTrackMatchedDepositedVsNclData->SetName(Form("%sData","fHistPIDAntiProtonsTrackMatchedDepositedVsNcl"));
  //fHistPiKPTrackMatchedDepositedVsNclData->SetName(Form("%sData","fHistPiKPTrackMatchedDepositedVsNcl"));
  fHistCentVsNchVsNclData->SetName("fHistCentVsNchVsNclRecoData");
  //TH3F *fHistCentVsNchVsNclData = lData->FindObject("fHistCentVsNchVsNcl");

  TCanvas *c1 = new TCanvas("c1","Simulation",600,400);
  c1->SetTopMargin(0.02);
  c1->SetRightMargin(0.149329);
  c1->SetBorderSize(0);
  c1->SetFillColor(0);
  c1->SetFillColor(0);
  c1->SetBorderMode(0);
  c1->SetFrameFillColor(0);
  c1->SetFrameBorderMode(0);
  //c1->SetLogz();
  //fHistPIDProtonsTrackMatchedDepositedVsNch->Draw("colz");
  TH1D *fHistPiKPTrackMatchedDepositedVsNchProf = fHistPiKPTrackMatchedDepositedVsNch->ProfileY();
  fHistPiKPTrackMatchedDepositedVsNchProf->Scale(1.0/30.0);
  fHistPiKPTrackMatchedDepositedVsNchProf->GetXaxis()->SetTitle("N_{ch}");
  fHistPiKPTrackMatchedDepositedVsNchProf->GetYaxis()->SetTitle("<E_{T}>");
  TH1D *fHistPiKPTrackMatchedDepositedVsNchProfCopy = fHistPiKPTrackMatchedDepositedVsNchProf->Clone("fHistPiKPTrackMatchedDepositedVsNchProfCopy");
  fHistPiKPTrackMatchedDepositedVsNchProfCopy->Draw();
  TH1D *fHistPIDAntiProtonsTrackMatchedDepositedVsNchProf = fHistPIDAntiProtonsTrackMatchedDepositedVsNch->ProfileY("fHistPIDAntiProtonsTrackMatchedDepositedVsNchProf");
  TH1D *fHistPIDAntiProtonsTrackMatchedDepositedVsNclProf = fHistPIDAntiProtonsTrackMatchedDepositedVsNcl->ProfileY("fHistPIDAntiProtonsTrackMatchedDepositedVsNclProf");
  fHistPIDAntiProtonsTrackMatchedDepositedVsNchProf->Draw("same");
  TH1D *fHistPIDProtonsTrackMatchedDepositedVsNchProf = fHistPIDProtonsTrackMatchedDepositedVsNch->ProfileY("fHistPIDProtonsTrackMatchedDepositedVsNchProf");
  TH1D *fHistPIDProtonsTrackMatchedDepositedVsNclProf = fHistPIDProtonsTrackMatchedDepositedVsNcl->ProfileY("fHistPIDProtonsTrackMatchedDepositedVsNclProf");
  fHistPIDProtonsTrackMatchedDepositedVsNchProf->Draw("same");
  SetStyles(fHistPiKPTrackMatchedDepositedVsNchProf,29,1);
  SetStyles(fHistPIDAntiProtonsTrackMatchedDepositedVsNchProf,20,TColor::kRed);
  SetStyles(fHistPIDProtonsTrackMatchedDepositedVsNchProf,24,TColor::kBlue);
  TLegend *leg = new TLegend(0.157718,0.709677,0.278523,0.935484);
  leg->SetFillStyle(0);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);
  leg->SetTextSize(0.038682);
  leg->AddEntry(fHistPiKPTrackMatchedDepositedVsNchProf,"Track matched E_{T}/30");
  leg->AddEntry(fHistPIDAntiProtonsTrackMatchedDepositedVsNchProf,"Identified #bar{p} E_{T}");
  leg->AddEntry(fHistPIDProtonsTrackMatchedDepositedVsNchProf,"Identified p E_{T}");
  leg->Draw();
  TString name1 = "/tmp/Sim"+detector+".png";
  c1->SaveAs(name1.Data());
  //cerr<<"176"<<endl;

  TCanvas *c2 = new TCanvas("c2","Data",600,400);
  c2->SetTopMargin(0.02);
  c2->SetRightMargin(0.149329);
  c2->SetBorderSize(0);
  c2->SetFillColor(0);
  c2->SetFillColor(0);
  c2->SetBorderMode(0);
  c2->SetFrameFillColor(0);
  c2->SetFrameBorderMode(0);
  //c2->SetLogz();
  //fHistPIDProtonsTrackMatchedDepositedVsNchData->Draw("colz");
  TH1D *fHistPiKPTrackMatchedDepositedVsNchDataProf = fHistPiKPTrackMatchedDepositedVsNchData->ProfileY("fHistPiKPTrackMatchedDepositedVsNchDataProf");
  fHistPiKPTrackMatchedDepositedVsNchDataProf->Scale(1.0/30.0);
  fHistPiKPTrackMatchedDepositedVsNchDataProf->GetXaxis()->SetTitle("N_{ch}");
  fHistPiKPTrackMatchedDepositedVsNchDataProf->GetYaxis()->SetTitle("<E_{T}>");
  TH1D *fHistPiKPTrackMatchedDepositedVsNchDataProfCopy = fHistPiKPTrackMatchedDepositedVsNchDataProf->Clone("fHistPiKPTrackMatchedDepositedVsNchDataProfCopy");
  fHistPiKPTrackMatchedDepositedVsNchDataProfCopy->Draw();
  TH1D *fHistPIDAntiProtonsTrackMatchedDepositedVsNchDataProf = fHistPIDAntiProtonsTrackMatchedDepositedVsNchData->ProfileY("fHistPIDAntiProtonsTrackMatchedDepositedVsNchDataProf");
  TH1D *fHistPIDAntiProtonsTrackMatchedDepositedVsNclDataProf = fHistPIDAntiProtonsTrackMatchedDepositedVsNclData->ProfileY("fHistPIDAntiProtonsTrackMatchedDepositedVsNclDataProf");
  fHistPIDAntiProtonsTrackMatchedDepositedVsNchDataProf->Draw("same");
  TH1D *fHistPIDProtonsTrackMatchedDepositedVsNchDataProf = fHistPIDProtonsTrackMatchedDepositedVsNchData->ProfileY("fHistPIDProtonsTrackMatchedDepositedVsNchDataProf");
  TH1D *fHistPIDProtonsTrackMatchedDepositedVsNclDataProf = fHistPIDProtonsTrackMatchedDepositedVsNclData->ProfileY("fHistPIDProtonsTrackMatchedDepositedVsNclDataProf");
  fHistPIDProtonsTrackMatchedDepositedVsNchDataProf->Draw("same");
  SetStyles(fHistPiKPTrackMatchedDepositedVsNchDataProf,29,1);
  SetStyles(fHistPIDAntiProtonsTrackMatchedDepositedVsNchDataProf,20,TColor::kRed);
  SetStyles(fHistPIDProtonsTrackMatchedDepositedVsNchDataProf,24,TColor::kBlue);
  TLegend *legData = new TLegend(0.157718,0.709677,0.278523,0.935484);
  legData->SetFillStyle(0);
  legData->SetFillColor(0);
  legData->SetBorderSize(0);
  legData->SetTextSize(0.03);
  legData->SetTextSize(0.038682);
  legData->AddEntry(fHistPiKPTrackMatchedDepositedVsNchDataProf,"Track matched E_{T}/30");
  legData->AddEntry(fHistPIDAntiProtonsTrackMatchedDepositedVsNchDataProf,"Identified #bar{p} E_{T}");
  legData->AddEntry(fHistPIDProtonsTrackMatchedDepositedVsNchDataProf,"Identified p E_{T}");
  legData->Draw();
  //cerr<<"212"<<endl;

  TString name2 = "/tmp/Data"+detector+".png";
  c2->SaveAs(name2.Data());
  fHistPiKPTrackMatchedDepositedVsNchDataProf->Scale(30.0);
  fHistPiKPTrackMatchedDepositedVsNchProf->Scale(30.0);

  TCanvas *c3 = new TCanvas("c3","Ratios at low pT",600,400);
  c3->SetTopMargin(0.02);
  c3->SetRightMargin(0.149329);
  c3->SetBorderSize(0);
  c3->SetFillColor(0);
  c3->SetFillColor(0);
  c3->SetBorderMode(0);
  c3->SetFrameFillColor(0);
  c3->SetFrameBorderMode(0);
  TH1D *hAntiProtonSim = Divide(fHistPIDAntiProtonsTrackMatchedDepositedVsNchProf,fHistPiKPTrackMatchedDepositedVsNchProf,"hAntiProtonSim");
  TH1D *hProtonSim = Divide(fHistPIDProtonsTrackMatchedDepositedVsNchProf,fHistPiKPTrackMatchedDepositedVsNchProf,"hProtonSim");
  TH1D *htemp =  fHistPIDProtonsTrackMatchedDepositedVsNchProf->Clone("temp");
  //cout<<htemp->GetEntries();
  htemp->Add(fHistPIDAntiProtonsTrackMatchedDepositedVsNchProf);
  htemp->Scale(2);//because the profile averages when you add two profiles
  //cout<<" "<<htemp->GetEntries()<<endl;
  //cout<<"entries "<<fHistPiKPTrackMatchedDepositedVsNchProf->GetEntries()<<endl;
  TH1D *hAllProtonSim = Divide(htemp,fHistPiKPTrackMatchedDepositedVsNchProf,"hAllProtonSim");
  //delete htemp;


  //cerr<<"240"<<endl;

  TH1D *hAntiProtonData = Divide(fHistPIDAntiProtonsTrackMatchedDepositedVsNchDataProf,fHistPiKPTrackMatchedDepositedVsNchDataProf,"hAntiProtonData");
  TH1D *hProtonData = Divide(fHistPIDProtonsTrackMatchedDepositedVsNchDataProf,fHistPiKPTrackMatchedDepositedVsNchDataProf,"hProtonData");
  TH1D *htemp2 =  fHistPIDProtonsTrackMatchedDepositedVsNchDataProf->Clone("temp2");

  //cout<<htemp2->GetEntries();
  htemp2->Add(fHistPIDAntiProtonsTrackMatchedDepositedVsNchDataProf);
  htemp2->Scale(2);//because the profile averages when you add two profiles
  //cout<<" "<<htemp2->GetEntries()<<endl;
  //cout<<"entries "<<fHistPiKPTrackMatchedDepositedVsNchDataProf->GetEntries()<<endl;
  TH1D *hAllProtonData = Divide(htemp2,fHistPiKPTrackMatchedDepositedVsNchDataProf,"hAllProtonData");
  //delete htemp2;
  SetStyles(hAntiProtonSim,20,TColor::kRed);
  SetStyles(hProtonSim,24,TColor::kRed);
  SetStyles(hAllProtonSim,28,1);
  SetStyles(hAntiProtonData,22,TColor::kRed);
  SetStyles(hProtonData,26,TColor::kRed);
  SetStyles(hAllProtonData,33,1);
  hAllProtonSim->GetXaxis()->SetTitle("N_ch");
  hAllProtonSim->GetYaxis()->SetTitle("E_{T}^{p,#bar{p}}/E_{T}^{#pi,K,p TM}");
  TF1 *funcData = new TF1("funcData","[0]*exp(-x/[1])+[2]-[3]*x",0,2500);
  funcData->SetParameter(0,0.12);
  funcData->SetParameter(1,10);
  funcData->SetParLimits(1,1e-9,1e9);
  funcData->SetParameter(2,0.09);
  funcData->SetParameter(3,2e-5);
  funcData->SetLineColor(hAllProtonData->GetLineColor());
  hAllProtonData->Fit(funcData,"","",50,2700);
  TF1 *funcSim = new TF1("funcSim","[0]*exp(-x/[1])+[2]-[3]*x",0,2500);
  funcSim->SetParameter(0,0.06);
  funcSim->SetParLimits(0,0,0.12);
  funcSim->SetParameter(1,10);
  funcSim->SetParLimits(1,1e-3,1e2);
  funcSim->SetParameter(2,0.09);
  funcSim->SetParameter(3,2e-5);
  funcSim->SetLineColor(hAllProtonSim->GetLineColor());
  hAllProtonSim->Fit(funcSim,"","",0,2700);
  hAllProtonSim->Draw();
  hAllProtonData->Draw("same");
  hAntiProtonData->Draw("same");
  hProtonData->Draw("same");
  hAntiProtonSim->Draw("same");
  hProtonSim->Draw("same");
  TLegend *legRatio = new TLegend(0.157718,0.709677,0.278523,0.935484);
  legRatio->SetFillStyle(0);
  legRatio->SetFillColor(0);
  legRatio->SetBorderSize(0);
  legRatio->SetTextSize(0.03);
  legRatio->SetTextSize(0.038682);
  legRatio->AddEntry(hAllProtonSim,"p+#bar{p} Sim");
  legRatio->AddEntry(hAllProtonData,"p+#bar{p} Data");
  legRatio->AddEntry(hProtonSim,"p Sim");
  legRatio->AddEntry(hProtonData,"p Data");
  legRatio->AddEntry(hAntiProtonSim,"#bar{p} Sim");
  legRatio->AddEntry(hAntiProtonData,"#bar{p} Data");
  legRatio->Draw();
  //     c1->cd();
  //     SetStyles(htemp,21,1);
  //     htemp->Draw("same");
  //     c2->cd();
  //     SetStyles(htemp2,21,1);
  //     htemp2->Draw("same");

  TString name3 = "/tmp/ProtonRatio"+detector+".eps";
  c3->SaveAs(name3.Data());

  cerr<<"306"<<endl;


  TCanvas *c4 = new TCanvas("c4","Sim: ratios of particles to low pT protons",600,400);
  c4->SetTopMargin(0.02);
  c4->SetRightMargin(0.149329);
  c4->SetBorderSize(0);
  c4->SetFillColor(0);
  c4->SetFillColor(0);
  c4->SetBorderMode(0);
  c4->SetFrameFillColor(0);
  c4->SetFrameBorderMode(0);
  TH1D *fHistNeutronsDepositedVsNchProf = fHistNeutronsDepositedVsNch->ProfileY("fHistNeutronsDepositedVsNchProf");
  TH1D *fHistAntiNeutronsDepositedVsNchProf = fHistAntiNeutronsDepositedVsNch->ProfileY("fHistAntiNeutronsDepositedVsNchProf");
  TH1D *fHistProtonsDepositedVsNchProf = fHistProtonsDepositedVsNch->ProfileY("fHistProtonsDepositedVsNchProf");
  TH1D *fHistAntiProtonsDepositedVsNchProf = fHistAntiProtonsDepositedVsNch->ProfileY("fHistAntiProtonsDepositedVsNchProf");
  TH1D *hNeutronsOverLowPtProtons = Divide(fHistNeutronsDepositedVsNchProf,fHistPIDProtonsTrackMatchedDepositedVsNchProf,"hAllNeutronsOverLowPtProtonsName");
  TH1D *hAntiNeutronsOverLowPtAntiProtons = Divide(fHistAntiNeutronsDepositedVsNchProf,fHistPIDAntiProtonsTrackMatchedDepositedVsNchProf,"hAllNeutronsOverLowPtAntiProtonsName");
  TH1D *hProtonsOverLowPtProtons = Divide(fHistProtonsDepositedVsNchProf,fHistPIDProtonsTrackMatchedDepositedVsNchProf,"hAllProtonsOverLowPtProtonsName");
  TH1D *hAntiProtonsOverLowPtAntiProtons = Divide(fHistAntiProtonsDepositedVsNchProf,fHistPIDAntiProtonsTrackMatchedDepositedVsNchProf,"hAllProtonsOverLowPtAntiProtonsName");

  TH1D *fHistNeutronsDepositedVsNclProf = fHistNeutronsDepositedVsNcl->ProfileY("fHistNeutronsDepositedVsNclProf");
  TH1D *fHistAntiNeutronsDepositedVsNclProf = fHistAntiNeutronsDepositedVsNcl->ProfileY("fHistAntiNeutronsDepositedVsNclProf");
  TH1D *fHistProtonsDepositedVsNclProf = fHistProtonsDepositedVsNcl->ProfileY("fHistProtonsDepositedVsNclProf");
  TH1D *fHistAntiProtonsDepositedVsNclProf = fHistAntiProtonsDepositedVsNcl->ProfileY("fHistAntiProtonsDepositedVsNclProf");
  cout<<"Dividing cluster histos"<<endl;
  TH1D *hNeutronsOverLowPtProtonsCl = Divide(fHistNeutronsDepositedVsNclProf,fHistPIDProtonsTrackMatchedDepositedVsNclProf,"hAllNeutronsOverLowPtProtonsClName");
  TH1D *hAntiNeutronsOverLowPtAntiProtonsCl = Divide(fHistAntiNeutronsDepositedVsNclProf,fHistPIDAntiProtonsTrackMatchedDepositedVsNclProf,"hAllNeutronsOverLowPtAntiProtonsClName");
  TH1D *hProtonsOverLowPtProtonsCl = Divide(fHistProtonsDepositedVsNclProf,fHistPIDProtonsTrackMatchedDepositedVsNclProf,"hAllProtonsOverLowPtProtonsClName");
  TH1D *hAntiProtonsOverLowPtAntiProtonsCl = Divide(fHistAntiProtonsDepositedVsNclProf,fHistPIDAntiProtonsTrackMatchedDepositedVsNclProf,"hAllProtonsOverLowPtAntiProtonsClName");
  cout<<"Done dividing cluster histos"<<endl;

  TH1D *hTmpAllLowPtProtons = fHistPIDProtonsTrackMatchedDepositedVsNchProf->Clone("hTmpAllLowPtProtons");
  hTmpAllLowPtProtons->Add(fHistPIDAntiProtonsTrackMatchedDepositedVsNchProf);
  TH1D *hTmpAllLowPtProtonsData = fHistPIDProtonsTrackMatchedDepositedVsNchDataProf->Clone("hTmpAllLowPtProtonsData");
  hTmpAllLowPtProtonsData->Add(fHistPIDAntiProtonsTrackMatchedDepositedVsNchDataProf);
  TH1D *hTmpAllProtons = fHistProtonsDepositedVsNchProf->Clone("hTmpAllProtons");
  hTmpAllProtons->Add(fHistAntiProtonsDepositedVsNchProf);
  TH1D *hTmpAllNeutrons = fHistNeutronsDepositedVsNchProf->Clone("hTmpAllNeutrons");
  hTmpAllNeutrons->Add(fHistAntiNeutronsDepositedVsNchProf);

  TH1D *hTmpAllLowPtProtonsCl = fHistPIDProtonsTrackMatchedDepositedVsNclProf->Clone("hTmpAllLowPtProtonsCl");
  hTmpAllLowPtProtonsCl->Add(fHistPIDAntiProtonsTrackMatchedDepositedVsNclProf);
  TH1D *hTmpAllLowPtProtonsDataCl = fHistPIDProtonsTrackMatchedDepositedVsNclDataProf->Clone("hTmpAllLowPtProtonsDataCl");
  hTmpAllLowPtProtonsDataCl->Add(fHistPIDAntiProtonsTrackMatchedDepositedVsNclDataProf);
  TH1D *hTmpAllProtonsCl = fHistProtonsDepositedVsNclProf->Clone("hTmpAllProtonsCl");
  hTmpAllProtonsCl->Add(fHistAntiProtonsDepositedVsNclProf);
  TH1D *hTmpAllNeutronsCl = fHistNeutronsDepositedVsNclProf->Clone("hTmpAllNeutronsCl");
  hTmpAllNeutronsCl->Add(fHistAntiNeutronsDepositedVsNclProf);

  TH1D *hAllNeutronsOverLowPtProtons = Divide(hTmpAllNeutrons,hTmpAllLowPtProtons,"hAllNeutronsOverLowPtProtons");
  TH1D *hAllNeutronsOverLowPtProtonsData = Divide(hTmpAllNeutrons,hTmpAllLowPtProtonsData,"hAllNeutronsOverLowPtProtonsData");
  TH1D *hAllProtonsOverLowPtProtons = Divide(hTmpAllProtons,hTmpAllLowPtProtons,"hAllProtonsOverLowPtProtons");
  TH1D *hAllProtonsOverLowPtProtonsData = Divide(hTmpAllProtons,hTmpAllLowPtProtonsData,"hAllProtonsOverLowPtProtonsData");
  SetStyles(hAntiNeutronsOverLowPtAntiProtons,20,TColor::kBlue);
  SetStyles(hNeutronsOverLowPtProtons,24,TColor::kBlue);
  SetStyles(hAntiProtonsOverLowPtAntiProtons,22,TColor::kRed);
  SetStyles(hProtonsOverLowPtProtons,26,TColor::kRed);
  SetStyles(hAllNeutronsOverLowPtProtonsData,29,TColor::kGreen+3);
  SetStyles(hAllNeutronsOverLowPtProtons,30,TColor::kGreen+3);
  SetStyles(hAllProtonsOverLowPtProtons,24,TColor::kRed);
  SetStyles(hAllProtonsOverLowPtProtonsData,20,TColor::kRed);

  TH1D *hAllNeutronsOverLowPtProtonsCl = Divide(hTmpAllNeutronsCl,hTmpAllLowPtProtonsCl,"hAllNeutronsOverLowPtProtonsCl");
  TH1D *hAllNeutronsOverLowPtProtonsDataCl = Divide(hTmpAllNeutronsCl,hTmpAllLowPtProtonsDataCl,"hAllNeutronsOverLowPtProtonsDataCl");
  TH1D *hAllProtonsOverLowPtProtonsCl = Divide(hTmpAllProtonsCl,hTmpAllLowPtProtonsCl,"hAllProtonsOverLowPtProtonsCl");
  TH1D *hAllProtonsOverLowPtProtonsDataCl = Divide(hTmpAllProtonsCl,hTmpAllLowPtProtonsDataCl,"hAllProtonsOverLowPtProtonsDataCl");
  SetStyles(hAntiNeutronsOverLowPtAntiProtonsCl,20,TColor::kBlue);
  SetStyles(hNeutronsOverLowPtProtonsCl,24,TColor::kBlue);
  SetStyles(hAntiProtonsOverLowPtAntiProtonsCl,22,TColor::kRed);
  SetStyles(hProtonsOverLowPtProtonsCl,26,TColor::kRed);
  SetStyles(hAllNeutronsOverLowPtProtonsDataCl,29,TColor::kGreen+3);
  SetStyles(hAllNeutronsOverLowPtProtonsCl,30,TColor::kGreen+3);
  SetStyles(hAllProtonsOverLowPtProtonsCl,24,TColor::kRed);
  SetStyles(hAllProtonsOverLowPtProtonsDataCl,20,TColor::kRed);

  //cerr<<"349"<<endl;

  Rescale(hAntiNeutronsOverLowPtAntiProtons,2);
  Rescale(hNeutronsOverLowPtProtons,2);
  Rescale(hAntiProtonsOverLowPtAntiProtons,2);
  Rescale(hProtonsOverLowPtProtons,2);
  Rescale(hAllNeutronsOverLowPtProtonsData,2);
  Rescale(hAllNeutronsOverLowPtProtons,2);
  Rescale(hAllProtonsOverLowPtProtons,2);
  Rescale(hAllProtonsOverLowPtProtonsData,2);

  Rescale(hAntiNeutronsOverLowPtAntiProtonsCl,2);
  Rescale(hNeutronsOverLowPtProtonsCl,2);
  Rescale(hAntiProtonsOverLowPtAntiProtonsCl,2);
  Rescale(hProtonsOverLowPtProtonsCl,2);
  Rescale(hAllNeutronsOverLowPtProtonsDataCl,2);
  Rescale(hAllNeutronsOverLowPtProtonsCl,2);
  Rescale(hAllProtonsOverLowPtProtonsCl,2);
  Rescale(hAllProtonsOverLowPtProtonsDataCl,2);


  //SetStyles(hAllProtonSim,28,1);
  //SetStyles(hAntiProtonData,22,TColor::kRed);
  //SetStyles(hProtonData,26,TColor::kRed);
  //SetStyles(hAllProtonData,33,1);
  hAntiNeutronsOverLowPtAntiProtons->SetMinimum(0.0);
  hAntiNeutronsOverLowPtAntiProtons->SetMaximum(8.0);
  hAntiNeutronsOverLowPtAntiProtons->GetXaxis()->SetTitle("N_{ch}");
  hAntiNeutronsOverLowPtAntiProtons->GetYaxis()->SetTitle("Ratio of energy deposited in sim");
  hAntiNeutronsOverLowPtAntiProtons->Draw("");
  hNeutronsOverLowPtProtons->Draw("same");
  hAntiProtonsOverLowPtAntiProtons->Draw("same");
  hProtonsOverLowPtProtons->Draw("same");
  hAllNeutronsOverLowPtProtons->Draw("same");
  hAllProtonsOverLowPtProtons->Draw("same");
  TLegend *legRatio2 = new TLegend(0.157718,0.709677,0.278523,0.935484);
  legRatio2->SetFillStyle(0);
  legRatio2->SetFillColor(0);
  legRatio2->SetBorderSize(0);
  legRatio2->SetTextSize(0.03);
  legRatio2->SetTextSize(0.038682);
  legRatio2->AddEntry(hAntiNeutronsOverLowPtAntiProtons,"#bar{n}/PID'd #bar{p}");
  legRatio2->AddEntry(hNeutronsOverLowPtProtons,"n/PID'd p");
  legRatio2->AddEntry(hAntiProtonsOverLowPtAntiProtons,"#bar{p}/PID'd #bar{p}");
  legRatio2->AddEntry(hProtonsOverLowPtProtons,"p/PID'd p");
  legRatio2->AddEntry(hAllNeutronsOverLowPtProtons,"(#bar{n}+n)/(PID'd #bar{p}+p)");
  legRatio2->AddEntry(hAllProtonsOverLowPtProtons,"(#bar{p}+p)/(PID'd #bar{p}+p)");
  legRatio2->Draw();

  //cerr<<"387"<<endl;
  TString name4 = "/tmp/RatiosToLowPt"+detector+".eps";
  c4->SaveAs(name4.Data());

  TCanvas *c5 = new TCanvas("c5","Neutron ratios used for error bounds",600,400);
  c5->SetTopMargin(0.02);
  c5->SetRightMargin(0.149329);
  c5->SetBorderSize(0);
  c5->SetFillColor(0);
  c5->SetFillColor(0);
  c5->SetBorderMode(0);
  c5->SetFrameFillColor(0);
  c5->SetFrameBorderMode(0);

  TF1 *funcNominal = new TF1("funcNominal","([0]*exp(-x/[1])+[2]-[3]*x)*[4]",0,3700);
  funcNominal->SetParameter(0,3.38452e+00);
  funcNominal->SetParameter(1,1.32723e+02);
  funcNominal->SetParameter(2,2.51044e+00);
  funcNominal->SetParameter(3,4.78729e-04);
  funcNominal->SetParLimits(1,1,1e3);
  funcNominal->FixParameter(4,1);
  //cout<<"funcnominal"<<endl;
  funcNominal->SetLineColor(hAllNeutronsOverLowPtProtonsData->GetLineColor());
  hAllNeutronsOverLowPtProtonsData->Fit(funcNominal,"","",50,2700);
  TF1 *funcLow = funcNominal->Clone("funcLow");
  funcLow->FixParameter(4,1.15);
  funcLow->SetLineStyle(2);
  TF1 *funcHigh = funcNominal->Clone("funcHigh");
  funcHigh->FixParameter(4,0.85);
  funcHigh->SetLineStyle(2);


  hAllNeutronsOverLowPtProtons->SetMinimum(0.0);
  hAllNeutronsOverLowPtProtons->SetMaximum(8.0);
  hAllNeutronsOverLowPtProtons->GetXaxis()->SetTitle("N_{ch}");
  hAllNeutronsOverLowPtProtons->GetYaxis()->SetTitle("Ratio of energy deposited in sim");
  hAllNeutronsOverLowPtProtons->Draw();
  hAllNeutronsOverLowPtProtonsData->Draw("same");
  hAllProtonsOverLowPtProtons->Draw("same");
  hAllProtonsOverLowPtProtonsData->Draw("same");
  TLegend *legRatioForData = new TLegend(0.157718,0.709677,0.278523,0.935484);
  legRatioForData->SetFillStyle(0);
  legRatioForData->SetFillColor(0);
  legRatioForData->SetBorderSize(0);
  legRatioForData->SetTextSize(0.03);
  legRatioForData->SetTextSize(0.038682);
  legRatioForData->AddEntry(hAllNeutronsOverLowPtProtons,"(#bar{n}+n)/(PID'd #bar{p}+p in Sim)");
  legRatioForData->AddEntry(hAllNeutronsOverLowPtProtonsData,"(#bar{n}+n)/(PID'd #bar{p}+p in Data)");
  legRatioForData->AddEntry(hAllProtonsOverLowPtProtons,"(#bar{p}+p)/(PID'd #bar{p}+p in Sim)");
  legRatioForData->AddEntry(hAllProtonsOverLowPtProtonsData,"(#bar{p}+p)/(PID'd #bar{p}+p in Data)");
  legRatioForData->Draw();
  funcHigh->Draw("same");
  funcLow->Draw("same");

  TString name5 = "/tmp/RatiosForErrors"+detector+".eps";
  c5->SaveAs(name5.Data());

  TCanvas *c5a = new TCanvas("c5a","Neutron ratios used for error bounds",600,400);
  c5a->SetTopMargin(0.02);
  c5a->SetRightMargin(0.149329);
  c5a->SetBorderSize(0);
  c5a->SetFillColor(0);
  c5a->SetFillColor(0);
  c5a->SetBorderMode(0);
  c5a->SetFrameFillColor(0);
  c5a->SetFrameBorderMode(0);

  TF1 *funcANominal = new TF1("funcANominal","([0]*exp(-x/[1])+[2]-[3]*x)*[4]",0,3700);
  funcANominal->SetParameter(0,3.09123e+00);
  funcANominal->SetParameter(1,2.69186e+01);
  funcANominal->SetParameter(2,3.07388e+00);
  funcANominal->SetParameter(3,2.88916e-03);
  funcANominal->SetParLimits(1,0,1e3);
  //funcANominal->SetParLimits(2,0,3);
  //funcANominal->SetParLimits(3,1e-3,1e3);
  funcANominal->SetParLimits(0,1,1e2);
  funcANominal->FixParameter(4,1);
  //cout<<"funcAnominal"<<endl;
  funcANominal->SetLineColor(hAllNeutronsOverLowPtProtonsDataCl->GetLineColor());
  Float_t fitlim = 280;
  if(isPhos) fitlim = 160;
  hAllNeutronsOverLowPtProtonsDataCl->Fit(funcANominal,"","",20,fitlim);
  TF1 *funcALow = funcANominal->Clone("funcALow");
  funcALow->FixParameter(4,1.15);
  funcALow->SetLineStyle(2);
  TF1 *funcAHigh = funcANominal->Clone("funcAHigh");
  funcAHigh->FixParameter(4,0.85);
  funcAHigh->SetLineStyle(2);


  hAllNeutronsOverLowPtProtonsCl->SetMinimum(0.0);
  hAllNeutronsOverLowPtProtonsCl->SetMaximum(8.0);
  hAllNeutronsOverLowPtProtonsCl->GetXaxis()->SetTitle("N_{cl}");
  hAllNeutronsOverLowPtProtonsCl->GetYaxis()->SetTitle("Ratio of energy deposited in sim");
  hAllNeutronsOverLowPtProtonsCl->GetXaxis()->SetRange(1,hAllNeutronsOverLowPtProtonsCl->FindBin(fitlim));
  hAllNeutronsOverLowPtProtonsCl->Draw();
  hAllNeutronsOverLowPtProtonsDataCl->Draw("same");
  hAllProtonsOverLowPtProtonsCl->Draw("same");
  hAllProtonsOverLowPtProtonsDataCl->Draw("same");
  TLegend *legRatioForDataCl = new TLegend(0.157718,0.709677,0.278523,0.935484);
  legRatioForDataCl->SetFillStyle(0);
  legRatioForDataCl->SetFillColor(0);
  legRatioForDataCl->SetBorderSize(0);
  legRatioForDataCl->SetTextSize(0.03);
  legRatioForDataCl->SetTextSize(0.038682);
  legRatioForDataCl->AddEntry(hAllNeutronsOverLowPtProtonsCl,"(#bar{n}+n)/(PID'd #bar{p}+p in Sim)");
  legRatioForDataCl->AddEntry(hAllNeutronsOverLowPtProtonsDataCl,"(#bar{n}+n)/(PID'd #bar{p}+p in Data)");
  legRatioForDataCl->AddEntry(hAllProtonsOverLowPtProtonsCl,"(#bar{p}+p)/(PID'd #bar{p}+p in Sim)");
  legRatioForDataCl->AddEntry(hAllProtonsOverLowPtProtonsDataCl,"(#bar{p}+p)/(PID'd #bar{p}+p in Data)");
  legRatioForDataCl->Draw();
  funcAHigh->Draw("same");
  funcALow->Draw("same");

  TString name5a = "/tmp/RatiosForErrorsCl"+detector+".eps";
  c5a->SaveAs(name5a.Data());



//     TCanvas *c5b = new TCanvas("c5b","Neutron ratios used for error bounds",600,400);
//     c5b->SetTopMargin(0.02);
//     c5b->SetRightMargin(0.149329);
//     c5b->SetBorderSize(0);
//     c5b->SetFillColor(0);
//     c5b->SetFillColor(0);
//     c5b->SetBorderMode(0);
//     c5b->SetFrameFillColor(0);
//     c5b->SetFrameBorderMode(0);
//     fHistNeutronsDepositedVsNclProf->Draw();
//     fHistAntiNeutronsDepositedVsNclProf->Draw("same");
//     fHistPIDProtonsTrackMatchedDepositedVsNclProf->Draw("same");
//     fHistPIDAntiProtonsTrackMatchedDepositedVsNclProf->Draw("same");
//     //fHistPIDProtonsTrackMatchedDepositedVsNclProf->Draw("same");
//     //fHistPIDAntiProtonsTrackMatchedDepositedVsNclProf->Draw("same");

//     SetStyles(fHistNeutronsDepositedVsNclProf,20, TColor::kBlue);
//     SetStyles(fHistAntiNeutronsDepositedVsNclProf,24, TColor::kBlue);
//     SetStyles(fHistPIDProtonsTrackMatchedDepositedVsNclProf,21,TColor::kRed);
//     SetStyles(fHistPIDAntiProtonsTrackMatchedDepositedVsNclProf,25,TColor::kRed);
//     PrintInfo(fHistNeutronsDepositedVsNclProf);
//     PrintInfo(fHistAntiNeutronsDepositedVsNclProf);
//     PrintInfo(fHistPIDProtonsTrackMatchedDepositedVsNclProf);
//     PrintInfo(fHistPIDAntiProtonsTrackMatchedDepositedVsNclProf);

    //cerr<<"442"<<endl;

    TObjArray trackmultiplicity(20);
    TObjArray trackmultiplicityData(20);
    TObjArray clustermultiplicity(20);
    TObjArray clustermultiplicityData(20);
    int nbinsChMult = fHistCentVsNchVsNcl->GetYaxis()->GetNbins();
    int nbinsClMult = fHistCentVsNchVsNcl->GetZaxis()->GetNbins();
    fHistCentVsNchVsNcl->GetXaxis()->SetTitle("Cent Bin");
    fHistCentVsNchVsNcl->GetYaxis()->SetTitle("N_{tr}");
    fHistCentVsNchVsNcl->GetZaxis()->SetTitle("N_{cl}");
    for(int cb=0;cb<20;cb++){
      //x axis is centrality
      //y axis is charged track multiplicity
      //z axis is cluster multiplicity
      fHistCentVsNchVsNcl->GetXaxis()->SetRange(cb+1,cb+1);
      trackmultiplicity[cb] = fHistCentVsNchVsNcl->Project3D("y");
      SetStyles((TH1*)trackmultiplicity[cb],markers[cb],colors[cb],Form("tr%i",cb));
      clustermultiplicity[cb] = fHistCentVsNchVsNcl->Project3D("z");
      SetStyles((TH1*)clustermultiplicity[cb],markers[cb],colors[cb],Form("cl%i",cb));
      fHistCentVsNchVsNclData->GetXaxis()->SetRange(cb+1,cb+1);
      trackmultiplicityData[cb] = fHistCentVsNchVsNclData->Project3D("y");
      SetStyles((TH1*)trackmultiplicityData[cb],markers[cb],colors[cb],Form("tr%i",cb));
      clustermultiplicityData[cb] = fHistCentVsNchVsNclData->Project3D("z");
      SetStyles((TH1*)clustermultiplicityData[cb],markers[cb],colors[cb],Form("cl%i",cb));
    }

    int currentShortCB = 0;
    //cerr<<"464"<<endl;
    Float_t neventsShortData[] = {0,0,0,0,0,0,0,0,0,0,0};
    Float_t neventsShortSim[] = {0,0,0,0,0,0,0,0,0,0,0};
    Float_t neventsShortDataCl[] = {0,0,0,0,0,0,0,0,0,0,0};
    Float_t neventsShortSimCl[] = {0,0,0,0,0,0,0,0,0,0,0};
    for(int cb=0;cb<19;cb++){
      //cout<<"cb "<<cb<<" current short cb "<<currentShortCB<<endl;
      cout<<"cb "<<cb<<" mean mult "<< ((TH1*)trackmultiplicity[cb])->GetMean()<<endl;
      int nbins = ((TH1*)trackmultiplicity[cb])->GetNbinsX();
      Float_t nevents = 0.0;
      for(int binn = 1;binn<=nbins;binn++){
      	float neventsBinn = ((TH1*)trackmultiplicityData[cb])->GetBinContent(binn);
	float mult= ((TH1*)trackmultiplicityData[cb])->GetBinCenter(binn);
	//mean et deposited in this event class =  (et of n+nbar / measured et of p+pbar) * (measured et of p+pbar)
	float meanetppbar = fHistPIDAntiProtonsTrackMatchedDepositedVsNchDataProf->GetBinContent(fHistPIDAntiProtonsTrackMatchedDepositedVsNchDataProf->FindBin(mult)) +  fHistPIDProtonsTrackMatchedDepositedVsNchDataProf->GetBinContent(fHistPIDProtonsTrackMatchedDepositedVsNchDataProf->FindBin(mult));
	float meanet = funcNominal->Eval(mult) * meanetppbar;
	nNeutronsData[cb] += neventsBinn*meanet;
	nevents +=neventsBinn;
	nNeutronsShortData[currentShortCB] += neventsBinn*meanet;
	neventsShortData[currentShortCB] +=neventsBinn;
      }
      if(nevents>0) nNeutronsData[cb] = nNeutronsData[cb]/nevents;
      nevents = 0.0;
      for(int binn = 1;binn<=nbins;binn++){
      	float neventsBinn = ((TH1*)trackmultiplicity[cb])->GetBinContent(binn);
	float mult= ((TH1*)trackmultiplicity[cb])->GetBinCenter(binn);
	//mean et deposited in this event class =  (et of n+nbar / measured et of p+pbar) * (measured et of p+pbar)
	float meanetppbar = fHistPIDAntiProtonsTrackMatchedDepositedVsNchProf->GetBinContent(fHistPIDAntiProtonsTrackMatchedDepositedVsNchProf->FindBin(mult)) +  fHistPIDProtonsTrackMatchedDepositedVsNchProf->GetBinContent(fHistPIDProtonsTrackMatchedDepositedVsNchProf->FindBin(mult));
	float meanet = funcNominal->Eval(mult) * meanetppbar;
	nNeutronsSim[cb] += neventsBinn*meanet;
	nevents +=neventsBinn;
	nNeutronsShortSim[currentShortCB] += neventsBinn*meanet;
	neventsShortSim[currentShortCB] +=neventsBinn;
      }
      if(nevents>0) nNeutronsSim[cb] = nNeutronsSim[cb]/nevents;


      nbins = ((TH1*)clustermultiplicity[cb])->GetNbinsX();
      nevents = 0.0;
      for(int binn = 1;binn<=nbins;binn++){
      	float neventsBinn = ((TH1*)clustermultiplicityData[cb])->GetBinContent(binn);
	float mult= ((TH1*)clustermultiplicityData[cb])->GetBinCenter(binn);
	//mean et deposited in this event class =  (et of n+nbar / measured et of p+pbar) * (measured et of p+pbar)
	float meanetppbar = fHistPIDAntiProtonsTrackMatchedDepositedVsNclDataProf->GetBinContent(fHistPIDAntiProtonsTrackMatchedDepositedVsNclDataProf->FindBin(mult)) +  fHistPIDProtonsTrackMatchedDepositedVsNclDataProf->GetBinContent(fHistPIDProtonsTrackMatchedDepositedVsNclDataProf->FindBin(mult));
	float meanet = funcANominal->Eval(mult) * meanetppbar;
	nNeutronsDataCl[cb] += neventsBinn*meanet;
	nevents +=neventsBinn;
	//cout<<"et "<< neventsBinn<<"*"<<meanet<<" events "<<neventsBinn<<endl;
	nNeutronsShortDataCl[currentShortCB] += neventsBinn*meanet;
	neventsShortDataCl[currentShortCB] +=neventsBinn;
      }
      if(nevents>0) nNeutronsDataCl[cb] = nNeutronsDataCl[cb]/nevents;
      nevents = 0.0;
      for(int binn = 1;binn<=nbins;binn++){
      	float neventsBinn = ((TH1*)clustermultiplicity[cb])->GetBinContent(binn);
	float mult= ((TH1*)clustermultiplicity[cb])->GetBinCenter(binn);
	//mean et deposited in this event class =  (et of n+nbar / measured et of p+pbar) * (measured et of p+pbar)
	float meanetppbar = fHistPIDAntiProtonsTrackMatchedDepositedVsNclProf->GetBinContent(fHistPIDAntiProtonsTrackMatchedDepositedVsNclProf->FindBin(mult)) +  fHistPIDProtonsTrackMatchedDepositedVsNclProf->GetBinContent(fHistPIDProtonsTrackMatchedDepositedVsNclProf->FindBin(mult));
	float meanet = funcANominal->Eval(mult) * meanetppbar;
	nNeutronsSimCl[cb] += neventsBinn*meanet;
	nevents +=neventsBinn;
	nNeutronsShortSimCl[currentShortCB] += neventsBinn*meanet;
	neventsShortSimCl[currentShortCB] +=neventsBinn;
      }
      if(nevents>0) nNeutronsSimCl[cb] = nNeutronsSimCl[cb]/nevents;


      //cout<<"cb "<<cb;
      //cout<<" data "<<  nNeutronsData[cb];// <<" +/- "<< 0.15*nNeutronsData[cb];
      //cout<<" sim "<<  nNeutronsSim[cb];// <<" +/- "<< 0.15*nNeutronsSim[cb]<<endl;
      //cout<<" data cl "<<  nNeutronsDataCl[cb];// <<" +/- "<< 0.15*nNeutronsDataCl[cb];
      //cout<<" sim cl "<<  nNeutronsSimCl[cb];// <<" +/- "<< 0.15*nNeutronsSimCl[cb]<<endl;
      //cout<<" w/err data ";
      float val = (nNeutronsData[cb]+nNeutronsDataCl[cb])/2.0;
      float err = TMath::Abs(nNeutronsData[cb]-nNeutronsDataCl[cb])/2.0;
      //cout<< val<<" +/- "<<err;
      //if(val>0)cout<<" ("<<err/val<<")";
      float valsim = (nNeutronsSim[cb]+nNeutronsSimCl[cb])/2.0;
      float errsim = TMath::Abs(nNeutronsSim[cb]-nNeutronsSimCl[cb])/2.0;
      //cout<<" sim "<< valsim <<" +/- "<<errsim;
      //if(valsim>0)cout<<" ("<<errsim/valsim<<")";
      //cout<<endl;
      nNeutronsData[cb] = val;
      nNeutronsDataErr[cb] = err;
      myfile<<Form("%2.3f %2.3f",nNeutronsData[cb],nNeutronsDataErr[cb])<<endl;
      if(cb<2 || cb%2==1){//normalize
	if(neventsShortSimCl[currentShortCB]>0) nNeutronsShortSimCl[currentShortCB] = nNeutronsShortSimCl[currentShortCB]/neventsShortSimCl[currentShortCB];
	if(neventsShortDataCl[currentShortCB]>0) nNeutronsShortDataCl[currentShortCB] = nNeutronsShortDataCl[currentShortCB]/neventsShortDataCl[currentShortCB];
	if(neventsShortSim[currentShortCB]>0) nNeutronsShortSim[currentShortCB] = nNeutronsShortSim[currentShortCB]/neventsShortSim[currentShortCB];
	if(neventsShortData[currentShortCB]>0) nNeutronsShortData[currentShortCB] = nNeutronsShortData[currentShortCB]/neventsShortData[currentShortCB];
	//cout<<"cbShort "<<currentShortCB<<" data "<<  nNeutronsShortData[currentShortCB] <<" +/- "<< 0.15*nNeutronsShortData[currentShortCB];
	//cout<<" sim "<<  nNeutronsShortSim[currentShortCB] <<" +/- "<< 0.15*nNeutronsShortSim[currentShortCB];//<<endl;

	val = (nNeutronsShortData[currentShortCB]+nNeutronsShortDataCl[currentShortCB])/2.0;
	err = TMath::Abs(nNeutronsShortData[currentShortCB]-nNeutronsShortDataCl[currentShortCB])/2.0;
	//cout<<" val "<<val<<" err "<<err<<endl;
	//cout<<nNeutronsShortData[currentShortCB]<<" "<<nNeutronsShortDataCl[currentShortCB]<<endl<<endl;
	nNeutronsShortData[currentShortCB] = val;
	nNeutronsShortDataErr[currentShortCB] = err;


	myfile2<<Form("%2.3f %2.3f",nNeutronsShortData[currentShortCB],0.15*nNeutronsShortData[currentShortCB])<<endl;
	if(cb<2) currentShortCB++;
      }
      if(cb>2 && cb%2==1){//increment the counter
	//cout<<"Here!"<<endl;
	currentShortCB++;
      }

    }
    myfile.close();
    myfile2.close();
    WriteLatex();
}


void WriteLatex(){
  TString detector = "Emcal";
    string inline;
    //NeutronsEmcal.dat
    TString neutronInfileName = "Neutrons"+detector+".dat";
    ifstream myneutronfile3 (neutronInfileName.Data());
    Float_t value = 0;
    Float_t error = 0;
    Int_t i=0;
    if (myneutronfile3.is_open()){
      while ( myneutronfile3.good() )
	{
	  getline (myneutronfile3,inline);
	  istringstream tmp(inline);
	  tmp >> value;
	  tmp >> error;
	  if(i<20){
	    neutronCorrEmcal[i] = value;
	    neutronErrorEmcal[i] = error;
	  }
	  i++;
	}
        myneutronfile3.close();
    }

    detector = "Phos";
    neutronInfileName = "Neutrons"+detector+".dat";
    ifstream myneutronfile4 (neutronInfileName.Data());
    Float_t value = 0;
    Float_t error = 0;
    Int_t i=0;
    if (myneutronfile4.is_open()){
      while ( myneutronfile4.good() )
	{
	  getline (myneutronfile4,inline);
	  istringstream tmp(inline);
	  tmp >> value;
	  tmp >> error;
	  if(i<20){
	    neutronCorrPhos[i] = value;
	    neutronErrorPhos[i] = error;
	  }
	  i++;
	}
        myneutronfile4.close();
    }
    ofstream myfile3;
    myfile3.open ("Protons.tex");
    for(int i=0;i<20;i++){
      TString line = Form("%i-%i & %2.3f $\\pm$ %2.3f & %2.3f $\\pm$ %2.3f \\\\",i*5,(i+1)*5,neutronCorrPhos[i],neutronErrorPhos[i],neutronCorrEmcal[i],neutronErrorEmcal[i]);
      myfile3<<line.Data()<<endl;

    }

    myfile3.close();

}


TH1* bayneseffdiv(TH1* numerator, TH1* denominator,Char_t* name) 
{
    if(!numerator){
      cerr<<"Error:  numerator does not exist!"<<endl;
      return NULL;
    }
    if(!denominator){
      cerr<<"Error:  denominator does not exist!"<<endl;
      return NULL;
    }
    TH1* result = (TH1*)numerator->Clone(name);
    Int_t nbins = numerator->GetNbinsX();
    for (Int_t ibin=0; ibin<= nbins+1; ++ibin) {
      Double_t numeratorVal = numerator->GetBinContent(ibin);
      Double_t denominatorVal = denominator->GetBinContent(ibin);
      // Check if the errors are right or the thing is scaled
      Double_t numeratorValErr = numerator->GetBinError(ibin);
      if (!(numeratorValErr==0. || numeratorVal ==0.) ) {
	Double_t rescale = numeratorValErr*numeratorValErr/numeratorVal;
	numeratorVal /= rescale;
      }
      Double_t denominatorValErr = denominator->GetBinError(ibin);
      if (!(denominatorValErr==0. || denominatorVal==0. )) {
	Double_t rescale = denominatorValErr*denominatorValErr/denominatorVal;
	denominatorVal /= rescale;
      }
      Double_t quotient = 0.;
      if (denominatorVal!=0.) {
	quotient = numeratorVal/denominatorVal;
      }
      Double_t quotientErr=0;
      quotientErr = TMath::Sqrt(
				(numeratorVal+1.0)/(denominatorVal+2.0)*
				((numeratorVal+2.0)/(denominatorVal+3.0)-(numeratorVal+1.0)/(denominatorVal+2.0)));
      result->SetBinContent(ibin,quotient);
      result->SetBinError(ibin,quotientErr);
      //cout<<"Setting bin "<<ibin<<" to "<<quotient<<" "<<numeratorVal<<"/"<<denominatorVal<<endl;
    }
    return result;
}
