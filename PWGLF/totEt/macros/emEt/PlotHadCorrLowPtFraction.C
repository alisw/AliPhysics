void SetStyles(TH1 *histo,int marker, int color){
  histo->Sumw2();
  histo->SetMarkerStyle(marker);
  histo->SetMarkerColor(color);
  histo->SetLineColor(color);
  //histo->GetXaxis()->SetTitle(xtitle);
  //histo->GetYaxis()->SetTitle(ytitle);
}
void PlotHadCorrLowPtFraction(Bool_t isPhos = kFALSE, Int_t mycase = 0){
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  TString filename, detname;
  if(mycase==0){
    if(isPhos){
      detname = "PHOS";
      //filename = "rootFiles/LHC11a10a_bis/Et.ESD.simPbPb.PHOS.LHC11a10a_bis.Run139465.root";
      filename = "rootFiles/LHC11a10a_bis/Et.ESD.simPbPb.PHOSOldHadMethod.LHC11a10a_bis.root";
    }
    else{
      //filename = "rootFiles/LHC11a10a_bis/Et.ESD.simPbPb.EMCal.LHC11a10a_bis.Run139465.root";
      filename = "rootFiles/LHC11a10a_bis/Et.ESD.simPbPb.EMCalOldHadMethod.LHC11a10a_bis.root";
      detname = "EMCal";
    }
  }
  else{
    if(isPhos){
      detname = "PHOS";
	filename = "rootFiles/LHC13e1abcCombined/Et.ESD.simPbPb.PHOS.LHC13e1abc.root";
    }
    else{
      filename = "rootFiles/LHC13e1abcCombined/Et.ESD.simPbPb.EMCal.LHC13e1abc.root";
      detname = "EMCal";
    }
  }
  
  TFile *f = TFile::Open(filename, "READ");
  TList *l = dynamic_cast<TList*>(f->Get("out1"));
  TH2F *fHistChargedTrackDepositsAcceptedVsPt = l->FindObject("fHistChargedTrackDepositsAcceptedVsPt");
  TH2F *fHistChargedTrackDepositsAllVsPt = l->FindObject("fHistChargedTrackDepositsAllVsPt");
    TCanvas *c1 = new TCanvas("c1","Simulation",600,400);
    c1->SetTopMargin(0.02);
    c1->SetRightMargin(0.03);
    c1->SetLeftMargin(0.11745);
    c1->SetBottomMargin(0.11745);
    c1->SetBorderSize(0);
    c1->SetFillColor(0);
    c1->SetFillColor(0);
    c1->SetBorderMode(0);
    c1->SetFrameFillColor(0);
    c1->SetFrameBorderMode(0);
    fHistChargedTrackDepositsAcceptedVsPt->Draw("colz");
    TH1D *proffHistChargedTrackDepositsAcceptedVsPt = fHistChargedTrackDepositsAcceptedVsPt->ProfileX("acceptedprofx");
    proffHistChargedTrackDepositsAcceptedVsPt->Draw("same");
    TCanvas *c2 = new TCanvas("c2","Simulation",600,400);
    c2->SetTopMargin(0.02);
    c2->SetRightMargin(0.03);
    c2->SetLeftMargin(0.11745);
    c2->SetBottomMargin(0.11745);
    c2->SetBorderSize(0);
    c2->SetFillColor(0);
    c2->SetFillColor(0);
    c2->SetBorderMode(0);
    c2->SetFrameFillColor(0);
    c2->SetFrameBorderMode(0);
    fHistChargedTrackDepositsAllVsPt->Draw("colz");
    TH1D *proffHistChargedTrackDepositsAllVsPt = fHistChargedTrackDepositsAllVsPt->ProfileX("allprofx");
    proffHistChargedTrackDepositsAllVsPt->Draw("same");
    TCanvas *c3 = new TCanvas("c3","Simulation",600,400);
    c3->SetTopMargin(0.02);
    c3->SetRightMargin(0.03);
    c3->SetLeftMargin(0.11745);
    c3->SetBottomMargin(0.11745);
    c3->SetBorderSize(0);
    c3->SetFillColor(0);
    c3->SetFillColor(0);
    c3->SetBorderMode(0);
    c3->SetFrameFillColor(0);
    c3->SetFrameBorderMode(0);
    TH1D *profyfHistChargedTrackDepositsAcceptedVsPt = fHistChargedTrackDepositsAcceptedVsPt->ProfileY("acceptedprofy");
    TH1D *profyfHistChargedTrackDepositsAllVsPt = fHistChargedTrackDepositsAllVsPt->ProfileY("allprofy");
    //profyfHistChargedTrackDepositsAllVsPt->Sumw2();
    //profyfHistChargedTrackDepositsAcceptedVsPt->Sumw2();
    TH1D *profyfHistChargedTrackDepositsAllVsPtClone = profyfHistChargedTrackDepositsAllVsPt->Clone("clone");
    profyfHistChargedTrackDepositsAllVsPt->Divide(profyfHistChargedTrackDepositsAcceptedVsPt);
    profyfHistChargedTrackDepositsAllVsPt->Draw();
    profyfHistChargedTrackDepositsAllVsPtClone->Add(profyfHistChargedTrackDepositsAcceptedVsPt,-1);
    profyfHistChargedTrackDepositsAllVsPtClone->Divide(profyfHistChargedTrackDepositsAcceptedVsPt);
    profyfHistChargedTrackDepositsAllVsPtClone->Draw("same");
    SetStyles(profyfHistChargedTrackDepositsAllVsPt,20,TColor::kBlue);
    SetStyles(profyfHistChargedTrackDepositsAllVsPtClone,29,TColor::kRed);
    Float_t low = 1.0;
    Float_t high = 0.0;
    for(int bin = 1;bin<17;bin++){
      if(profyfHistChargedTrackDepositsAllVsPtClone->GetBinContent(bin)<low) low = profyfHistChargedTrackDepositsAllVsPtClone->GetBinContent(bin);
      if(profyfHistChargedTrackDepositsAllVsPtClone->GetBinContent(bin)>high) high = profyfHistChargedTrackDepositsAllVsPtClone->GetBinContent(bin);
    }
    //cout<<"low "<<low<<" high "<<high<<" mean "<<(low+high)/2.0<<" +/- "<<(high-low)/2.0<<endl;
    cout<<"corrfac = "<<(low+high)/2.0<<";"<<endl;
    cout<<"corrfacerr = "<<(high-low)/2.0<<";"<<endl;
    cout<<Form("%2.3f $\\pm$ %2.3f",(low+high)/2.0,(high-low)/2.0)<<endl;
}
