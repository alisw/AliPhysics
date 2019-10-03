AliCFContainer *MakeSlice(AliCFContainer *cont, Double_t pmin = 2., Int_t charge = 0){
  TArrayI steps(cont->GetNStep()); for(Int_t istep = 0; istep < cont->GetNStep(); istep++) steps[istep] = istep;
  Int_t vars[2] = {1,2};
  TArrayD min(cont->GetNVar()), max(cont->GetNVar());
  for(Int_t ivar = 0; ivar < cont->GetNVar(); ivar++){
    if(ivar == 0){
      min[ivar] = pmin;
      max[ivar] = 10.;
    } else if(ivar == 3){
      if(charge == -1){
        min[ivar] = max[ivar] = -1;
      } else if(charge == 1){
        min[ivar] = max[ivar] = 1;
      } else {
        min[ivar] = 1;
        max[ivar] = -1;
      }
    } else {
      min[ivar] = cont->GetAxis(ivar, 0)->GetXmin();
      max[ivar] = cont->GetAxis(ivar, 0)->GetXmax();
    }
  }
  return cont->MakeSlice(cont->GetNStep(), steps.GetArray(), 2, vars, min.GetArray(), max.GetArray());
}

void DrawEtaPhi(Double_t pmin = 2.){
  TFile *in = TFile::Open("HFEtaskTRD2.root");
  TList *res = (TList *)in->Get("TRD_HFE_Results2");
  // Container with tracks
  AliHFEcontainer *hfecont = (AliHFEcontainer *)res->FindObject("trackContainer");
  AliCFContainer *cont = hfecont->GetCFContainer("recTrackContReco");
  // Container with number of events for normalization
  AliCFContainer *evc = (AliCFContainer *)res->FindObject("eventContainer");
  TH1 *h1 = (TH1 *)evc->Project(4, 0);
  Int_t nEvents = h1->GetEntries();
  cout << "Number of events: " << nEvents << endl;
  // For plotting
  Double_t rm = 0.2;
  Double_t ts = 0.045;

  Int_t steps = cont->GetNStep() - 1;
  TH2 *htpr = 0;
  TCanvas *plotall = new TCanvas("plot", "Eta-phi maps for different cut steps (all charge)", 1000, 600);
  plotall->Divide(steps/2, 2);
  AliCFContainer *scall= MakeSlice(cont, pmin, 0);
  for(Int_t istep = 1; istep < cont->GetNStep(); istep++){
    plotall->cd(istep);
    gPad->SetRightMargin(rm);
    hptr = (TH2 *)scall->Project(istep, 0,1);
    hptr->Scale(1e5/static_cast<Double_t>(nEvents));
    hptr->SetStats(0);
    hptr->SetName(Form("%sall", cont->GetStepTitle(istep)));
    hptr->SetTitle(Form("Step %s", cont->GetStepTitle(istep)));
    hptr->GetXaxis()->SetTitle("#eta");
    hptr->GetYaxis()->SetTitle("#phi");
    hptr->GetZaxis()->SetTitle("Tracks / 100000 Events");
    hptr->GetZaxis()->SetTitleOffset(1.2);
    hptr->GetXaxis()->SetTitleSize(ts);
    hptr->GetYaxis()->SetTitleSize(ts);
    hptr->GetZaxis()->SetTitleSize(ts);
    hptr->GetZaxis()->SetRangeUser(0, 2.2);
    hptr->Draw("colz");
  }
  plotall->cd();

  TCanvas *plotpos = new TCanvas("plotpos", "Eta-phi maps for different cut steps (pos charge)", 1000, 600);
  plotpos->Divide(steps/2, 2);
  AliCFContainer *scpos = MakeSlice(cont, pmin, 1);
  for(Int_t istep = 1; istep < cont->GetNStep(); istep++){
    plotpos->cd(istep);
    gPad->SetRightMargin(rm);
    hptr = (TH2 *)scpos->Project(istep, 0,1);
    hptr->Scale(1e5/static_cast<Double_t>(nEvents));
    hptr->SetStats(0);
    hptr->SetTitle(Form("Step %s", cont->GetStepTitle(istep)));
    hptr->SetName(Form("%spos", cont->GetStepTitle(istep)));
    hptr->GetXaxis()->SetTitle("#eta");
    hptr->GetYaxis()->SetTitle("#phi");
    hptr->GetZaxis()->SetTitle("Tracks / 100000 Events");
    hptr->GetZaxis()->SetTitleOffset(1.2);
    hptr->GetXaxis()->SetTitleSize(ts);
    hptr->GetYaxis()->SetTitleSize(ts);
    hptr->GetZaxis()->SetTitleSize(ts);
    hptr->GetZaxis()->SetRangeUser(0, 1.1);
    hptr->Draw("colz");
  }
  plotpos->cd();

  TCanvas *plotneg = new TCanvas("plotneg", "Eta-phi maps for different cut steps (neg charge)", 1000, 600);
  plotneg->Divide(steps/2, 2);
  AliCFContainer *scneg = MakeSlice(cont, pmin, -1);
  for(Int_t istep = 1; istep < cont->GetNStep(); istep++){
    plotneg->cd(istep);
    gPad->SetRightMargin(rm);
    hptr = (TH2 *)scneg->Project(istep, 0,1);
    hptr->Scale(1e5/static_cast<Double_t>(nEvents));
    hptr->SetStats(0);
    hptr->SetName(Form("%sneg", cont->GetStepTitle(istep)));
    hptr->SetTitle(Form("Step %s", cont->GetStepTitle(istep)));
    hptr->GetXaxis()->SetTitle("#eta");
    hptr->GetYaxis()->SetTitle("#phi");
    hptr->GetZaxis()->SetTitle("Tracks / 100000 Events");
    hptr->GetZaxis()->SetTitleOffset(1.2);
    hptr->GetXaxis()->SetTitleSize(ts);
    hptr->GetYaxis()->SetTitleSize(ts);
    hptr->GetZaxis()->SetTitleSize(ts);
    hptr->GetZaxis()->SetRangeUser(0, 1.1);
    hptr->Draw("colz");
  }
  plotneg->cd();

}
