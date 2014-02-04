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
Int_t colors[] = {0,TColor::kRed, TColor::kOrange, TColor::kGreen+3, TColor::kBlue, TColor::kBlack, 
		    TColor::kRed, TColor::kOrange, TColor::kGreen+3, TColor::kBlue, TColor::kBlack, 
		    TColor::kRed, TColor::kOrange, TColor::kGreen+3, TColor::kBlue, TColor::kBlack, 
		    TColor::kRed, TColor::kOrange, TColor::kGreen+3, TColor::kBlue, TColor::kBlack};
Int_t markers[] = {20,21,22,23,33, 24,25,26,32,27, 20,21,22,23,33, 24,25,26,32,27};


Float_t nNeutronsSim[] = {0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0};
Float_t nNeutronsData[] = {0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0};
Float_t nNeutronsShortSim[] = {0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0};
Float_t nNeutronsShortData[] = {0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0};

TH1D* Divide(TH1D* numerator, TH1D* denominator,Char_t* name){
  //TH1D *output = numerator->Clone(name);
  TH1D *output = new TH1D(name,name,numerator->GetNbinsX(),numerator->GetBinLowEdge(1),numerator->GetBinLowEdge(numerator->GetNbinsX()+1));
  output->GetXaxis()->SetTitle(numerator->GetXaxis()->GetTitle());
  output->GetYaxis()->SetTitle(numerator->GetYaxis()->GetTitle());
  for(Int_t i=1;i<=output->GetNbinsX();i++){
    output->SetBinContent(i,numerator->GetBinContent(i));
    //output->SetBinError((Int_t)i,0,0,(Double_t)0.0);
    //denominator->SetBinError((Int_t)i,0,0,(Double_t)0.0);
    //cout<<" "<<output->GetBinError(i)<<" "<<denominator->GetBinError(i)<<endl;
  }
  output->Divide(denominator);
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

void PlotProtons(TString filename = "rootFiles/LHC11a10a_bis/Et.ESD.simPbPb.EMCal.LHC11a10a_bis.Run139465.root", TString filenameData = "rootFiles/LHC10hPass2/Et.ESD.realPbPb.EMCal.LHC10hPass2.Run139465.root"){
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  TString detector = "Phos";
  if(filename.Contains("EMC")){
    detector = "Emcal";
  }
  ofstream myfile;
  TString textfilename = "Neutrons"+detector+".dat";
  myfile.open (textfilename.Data());
  ofstream myfile2;
  TString textfilename2 = "Neutrons"+detector+"Short.dat";
  myfile2.open (textfilename2.Data());


    TFile *f = TFile::Open(filename, "READ");
    TList *l = dynamic_cast<TList*>(f->Get("out1"));
    TH2F *fHistPIDProtonsTrackMatchedDepositedVsNch = l->FindObject("fHistPIDProtonsTrackMatchedDepositedVsNch");
    TH2F *fHistPIDAntiProtonsTrackMatchedDepositedVsNch = l->FindObject("fHistPIDAntiProtonsTrackMatchedDepositedVsNch");
    TH2F *fHistPiKPTrackMatchedDepositedVsNch = l->FindObject("fHistPiKPTrackMatchedDepositedVsNch");
    fHistPIDProtonsTrackMatchedDepositedVsNch->SetName(Form("%sSim","fHistPIDProtonsTrackMatchedDepositedVsNch"));
    fHistPIDAntiProtonsTrackMatchedDepositedVsNch->SetName(Form("%sSim","fHistPIDAntiProtonsTrackMatchedDepositedVsNch"));
    fHistPiKPTrackMatchedDepositedVsNch->SetName(Form("%sSim","fHistPiKPTrackMatchedDepositedVsNch"));

    TH2F *fHistNeutronsDepositedVsNch = l->FindObject("fHistNeutronsDepositedVsNch");
    TH2F *fHistAntiNeutronsDepositedVsNch = l->FindObject("fHistAntiNeutronsDepositedVsNch");
    TH2F *fHistProtonsDepositedVsNch = l->FindObject("fHistProtonsDepositedVsNch");
    TH2F *fHistAntiProtonsDepositedVsNch = l->FindObject("fHistAntiProtonsDepositedVsNch");
    TH2F *fHistPiKPDepositedVsNch = l->FindObject("fHistPiKPDepositedVsNch");

    TH2F *fHistProtonsNotTrackMatchedDepositedVsNch = l->FindObject("fHistProtonsNotTrackMatchedDepositedVsNch");
    TH2F *fHistAntiProtonsNotTrackMatchedDepositedVsNch = l->FindObject("fHistAntiProtonsNotTrackMatchedDepositedVsNch");

    TH3F *fHistCentVsNchVsNcl = l->FindObject("fHistCentVsNchVsNclReco");
    //fHistCentVsNchVsNcl->SetName("fHistCentVsNchVsNclSim");

    TFile *fData = TFile::Open(filenameData, "READ");
    TList *lData = dynamic_cast<TList*>(fData->Get("out1"));
    TH2F *fHistPIDProtonsTrackMatchedDepositedVsNchData = lData->FindObject("fHistPIDProtonsTrackMatchedDepositedVsNch");
    TH2F *fHistPIDAntiProtonsTrackMatchedDepositedVsNchData = lData->FindObject("fHistPIDAntiProtonsTrackMatchedDepositedVsNch");
    TH2F *fHistPiKPTrackMatchedDepositedVsNchData = lData->FindObject("fHistPiKPTrackMatchedDepositedVsNch");
    fHistPIDProtonsTrackMatchedDepositedVsNchData->SetName(Form("%sData","fHistPIDProtonsTrackMatchedDepositedVsNch"));
    fHistPIDAntiProtonsTrackMatchedDepositedVsNchData->SetName(Form("%sData","fHistPIDAntiProtonsTrackMatchedDepositedVsNch"));
    fHistPiKPTrackMatchedDepositedVsNchData->SetName(Form("%sData","fHistPiKPTrackMatchedDepositedVsNch"));
    TH3F *fHistCentVsNchVsNclData = lData->FindObject("fHistCentVsNchVsNclReco");
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
    TH1D *fHistPIDAntiProtonsTrackMatchedDepositedVsNchProf = fHistPIDAntiProtonsTrackMatchedDepositedVsNch->ProfileY();
    fHistPIDAntiProtonsTrackMatchedDepositedVsNchProf->Draw("same");
    TH1D *fHistPIDProtonsTrackMatchedDepositedVsNchProf = fHistPIDProtonsTrackMatchedDepositedVsNch->ProfileY();
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
    TH1D *fHistPiKPTrackMatchedDepositedVsNchDataProf = fHistPiKPTrackMatchedDepositedVsNchData->ProfileY();
    fHistPiKPTrackMatchedDepositedVsNchDataProf->Scale(1.0/30.0);
    fHistPiKPTrackMatchedDepositedVsNchDataProf->GetXaxis()->SetTitle("N_{ch}");
    fHistPiKPTrackMatchedDepositedVsNchDataProf->GetYaxis()->SetTitle("<E_{T}>");
    TH1D *fHistPiKPTrackMatchedDepositedVsNchDataProfCopy = fHistPiKPTrackMatchedDepositedVsNchDataProf->Clone("fHistPiKPTrackMatchedDepositedVsNchDataProfCopy");
    fHistPiKPTrackMatchedDepositedVsNchDataProfCopy->Draw();
    TH1D *fHistPIDAntiProtonsTrackMatchedDepositedVsNchDataProf = fHistPIDAntiProtonsTrackMatchedDepositedVsNchData->ProfileY();
    fHistPIDAntiProtonsTrackMatchedDepositedVsNchDataProf->Draw("same");
    TH1D *fHistPIDProtonsTrackMatchedDepositedVsNchDataProf = fHistPIDProtonsTrackMatchedDepositedVsNchData->ProfileY();
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
    funcData->SetParameter(2,0.09);
    funcData->SetParameter(3,2e-5);
    funcData->SetLineColor(hAllProtonData->GetLineColor());
    hAllProtonData->Fit(funcData,"","",50,2700);
    TF1 *funcSim = new TF1("funcSim","[0]*exp(-x/[1])+[2]-[3]*x",0,2500);
    funcSim->SetParameter(0,0.06);
    funcSim->SetParLimits(0,0,0.12);
    funcSim->SetParameter(1,10);
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

    TString name3 = "/tmp/Ratio"+detector+".png";
    c3->SaveAs(name3.Data());

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
    TH1D *hNeutronsOverLowPtProtons = Divide(fHistNeutronsDepositedVsNchProf,fHistPIDProtonsTrackMatchedDepositedVsNchProf,"hAllNeutronsOverLowPtProtons");
    TH1D *hAntiNeutronsOverLowPtAntiProtons = Divide(fHistAntiNeutronsDepositedVsNchProf,fHistPIDAntiProtonsTrackMatchedDepositedVsNchProf,"hAllNeutronsOverLowPtAntiProtons");
    TH1D *hProtonsOverLowPtProtons = Divide(fHistProtonsDepositedVsNchProf,fHistPIDProtonsTrackMatchedDepositedVsNchProf,"hAllProtonsOverLowPtProtons");
    TH1D *hAntiProtonsOverLowPtAntiProtons = Divide(fHistAntiProtonsDepositedVsNchProf,fHistPIDAntiProtonsTrackMatchedDepositedVsNchProf,"hAllProtonsOverLowPtAntiProtons");

    TH1D *hTmpAllLowPtProtons = fHistPIDProtonsTrackMatchedDepositedVsNchProf->Clone("hTmpAllLowPtProtons");
    hTmpAllLowPtProtons->Add(fHistPIDAntiProtonsTrackMatchedDepositedVsNchProf);
    TH1D *hTmpAllLowPtProtonsData = fHistPIDProtonsTrackMatchedDepositedVsNchDataProf->Clone("hTmpAllLowPtProtonsData");
    hTmpAllLowPtProtonsData->Add(fHistPIDAntiProtonsTrackMatchedDepositedVsNchDataProf);
    TH1D *hTmpAllProtons = fHistProtonsDepositedVsNchProf->Clone("hTmpAllProtons");
    hTmpAllProtons->Add(fHistAntiProtonsDepositedVsNchProf);
    TH1D *hTmpAllNeutrons = fHistNeutronsDepositedVsNchProf->Clone("hTmpAllNeutrons");
    hTmpAllNeutrons->Add(fHistAntiNeutronsDepositedVsNchProf);

    TH1D *hAllNeutronsOverLowPtProtons = Divide(hTmpAllNeutrons,hTmpAllLowPtProtons,"hAllNeutronsOverLowPtProtons");
    TH1D *hAllNeutronsOverLowPtProtonsData = Divide(hTmpAllNeutrons,hTmpAllLowPtProtonsData,"hAllNeutronsOverLowPtProtonsData");
    TH1D *hAllProtonsOverLowPtProtons = Divide(hTmpAllProtons,hTmpAllLowPtProtons,"hAllProtonsOverLowPtProtons");
    TH1D *hAllProtonsOverLowPtProtonsData = Divide(hTmpAllProtons,hTmpAllLowPtProtonsData,"hAllProtonsOverLowPtProtonsData");
    SetStyles(hAntiNeutronsOverLowPtAntiProtons,20,TColor::kBlue);
    SetStyles(hNeutronsOverLowPtProtons,24,TColor::kBlue);
    SetStyles(hAntiProtonsOverLowPtAntiProtons,22,TColor::kRed);
    SetStyles(hProtonsOverLowPtProtons,26,TColor::kRed);
    SetStyles(hAllNeutronsOverLowPtProtonsData,22,TColor::kBlue+3);
    SetStyles(hAllNeutronsOverLowPtProtons,28,TColor::kGreen+3);
    SetStyles(hAllProtonsOverLowPtProtons,33,TColor::kGreen+3);
    SetStyles(hAllProtonsOverLowPtProtonsData,24,TColor::kBlue+3);

    Rescale(hAntiNeutronsOverLowPtAntiProtons,2);
    Rescale(hNeutronsOverLowPtProtons,2);
    Rescale(hAntiProtonsOverLowPtAntiProtons,2);
    Rescale(hProtonsOverLowPtProtons,2);
    Rescale(hAllNeutronsOverLowPtProtonsData,2);
    Rescale(hAllNeutronsOverLowPtProtons,2);
    Rescale(hAllProtonsOverLowPtProtons,2);
    Rescale(hAllProtonsOverLowPtProtonsData,2);
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

    TString name4 = "/tmp/RatiosToLowPt"+detector+".png";
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
    funcNominal->SetParameter(0,0.12);
    funcNominal->SetParameter(1,10);
    funcNominal->SetParameter(2,0.09);
    funcNominal->SetParameter(3,2e-5);
    funcNominal->FixParameter(4,1);
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

    TString name5 = "/tmp/RatiosForErrors"+detector+".png";
    c5->SaveAs(name5.Data());


    TObjArray trackmultiplicity(20);
    TObjArray trackmultiplicityData(20);
    int nbinsChMult = fHistCentVsNchVsNcl->GetYaxis()->GetNbins();
    int nbinsClMult = fHistCentVsNchVsNcl->GetZaxis()->GetNbins();
    fHistCentVsNchVsNcl->GetXaxis()->SetTitle("Cent Bin");
    fHistCentVsNchVsNcl->GetYaxis()->SetTitle("N_{tr}");
    fHistCentVsNchVsNcl->GetZaxis()->SetTitle("N_{cl}");
    int currentShortCB = 0;
    for(int cb=0;cb<20;cb++){
      //x axis is centrality
      //y axis is charged track multiplicity
      //z axis is cluster multiplicity
      fHistCentVsNchVsNcl->GetXaxis()->SetRange(cb+1,cb+1);
      trackmultiplicity[cb] = fHistCentVsNchVsNcl->Project3D("y");
      SetStyles((TH1*)trackmultiplicity[cb],markers[cb],colors[cb],Form("tr%i",cb));
      fHistCentVsNchVsNclData->GetXaxis()->SetRange(cb+1,cb+1);
      trackmultiplicityData[cb] = fHistCentVsNchVsNclData->Project3D("y");
      SetStyles((TH1*)trackmultiplicityData[cb],markers[cb],colors[cb],Form("tr%i",cb));
    }

    Float_t neventsShortData[] = {0,0,0,0,0,0,0,0,0,0,0};
    Float_t neventsShortSim[] = {0,0,0,0,0,0,0,0,0,0,0};
    for(int cb=0;cb<19;cb++){
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
      nNeutronsData[cb] = nNeutronsData[cb]/nevents;
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
      nNeutronsSim[cb] = nNeutronsSim[cb]/nevents;
      cout<<"cb "<<cb<<" data "<<  nNeutronsData[cb] <<" +/- "<< 0.15*nNeutronsData[cb];
      cout<<" sim "<<  nNeutronsSim[cb] <<" +/- "<< 0.15*nNeutronsSim[cb]<<endl;
      myfile<<Form("%2.3f %2.3f",nNeutronsData[cb],0.15*nNeutronsData[cb])<<endl;
      if(cb<2 || cb%2==1){//normalize
	nNeutronsShortSim[currentShortCB] = nNeutronsShortSim[currentShortCB]/neventsShortSim[currentShortCB];
	nNeutronsShortData[currentShortCB] = nNeutronsShortData[currentShortCB]/neventsShortData[currentShortCB];
	cout<<"cbShort "<<currentShortCB<<" data "<<  nNeutronsShortData[currentShortCB] <<" +/- "<< 0.15*nNeutronsShortData[currentShortCB];
	cout<<" sim "<<  nNeutronsShortSim[currentShortCB] <<" +/- "<< 0.15*nNeutronsShortSim[currentShortCB]<<endl;

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
