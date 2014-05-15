// Author: Ionut-Cristian Arsene
// Contacts: iarsene@cern.ch
// Make QA pictures for the J/psi->ee analysis



// modified by sjena 
// -- to add it in automatic script
// -- changes Naming conventions
// -- saving objects instead of canvas


//______________________________________________________________________________
void processJpsi2eeQAplots(const Char_t* filename="jpsi_Default.root", 
			   TString suffix = "eps",
			   const Char_t* outfile="Jpsi2eeQAplots_output.root") {
  //
 
  //  
  TFile* file = TFile::Open(filename);
  
  // event wise histograms 
  TH1F* zDistrib = (TH1F*)GetHistogram(file, "default", "Event", "Z");
  
  // electron candidate histograms
  TH2F* tpcDedx = (TH2F*)GetHistogram(file, "default", "Track_ev1+", "P_InnerParam_TPC_signal");
  tpcDedx->Add((TH2F*)GetHistogram(file, "default", "Track_ev1-", "P_InnerParam_TPC_signal"));
  TH2F* tpcNsigmaEle = (TH2F*)GetHistogram(file, "default", "Track_ev1+", "P_InnerParam_TPC_nSigma_Electrons");
  tpcNsigmaEle->Add((TH2F*)GetHistogram(file, "default", "Track_ev1-", "P_InnerParam_TPC_nSigma_Electrons"));
  
  // pair candidate histograms
  TH1F* invmass_pp = (TH1F*)GetHistogram(file, "default", "Pair_ev1+_ev1+", "pM");
  TH1F* invmass_pm = (TH1F*)GetHistogram(file, "default", "Pair_ev1+_ev1-", "pM");
  TH1F* invmass_mm = (TH1F*)GetHistogram(file, "default", "Pair_ev1-_ev1-", "pM");
    
  // draw stuff
  TLatex* latex=new TLatex();
  latex->SetNDC();
  latex->SetTextSize(0.04);
  latex->SetTextFont(42);
  TCanvas* c1=new TCanvas("fig_dq_tpcDedx", "");
  if(tpcDedx) {
    tpcDedx->SetStats(kFALSE);
    tpcDedx->GetXaxis()->SetRangeUser(0.0,10.0);
    tpcDedx->GetYaxis()->SetRangeUser(40.0,120.0);
    tpcDedx->SetTitle("");
    tpcDedx->Draw("colz");
    latex->DrawLatex(0.12, 0.83, "J/#psi electron candidates");
  }
  
  TCanvas* c2=new TCanvas("fig_dq_tpcNsigmaElectron", "");
  if(tpcNsigmaEle) {
    tpcNsigmaEle->SetStats(kFALSE);
    tpcNsigmaEle->GetYaxis()->SetRangeUser(-5.0,5.0);
    tpcNsigmaEle->GetXaxis()->SetRangeUser(0.0,10.0);
    tpcNsigmaEle->SetTitle("");
    tpcNsigmaEle->Draw("colz");
    latex->DrawLatex(0.12, 0.83, "J/#psi electron candidates");
  }
  
  TCanvas* c3=new TCanvas("fig_dq_eeInvMass", "");
  if(invmass_pm) {
    invmass_pm->SetStats(kFALSE);
    invmass_pm->SetTitle("");
    invmass_pm->SetLineColor(1);
    invmass_pm->GetYaxis()->SetTitle(Form("Pairs per %.0f MeV/c^{2}", 1000.0*invmass_pm->GetXaxis()->GetBinWidth(1)));
    invmass_pm->GetXaxis()->SetTitle("m_{e^{+}e^{-}} (GeV/c^{2})");
    invmass_pm->Draw();
  }
  if(invmass_pp) {
    invmass_pp->SetLineColor(2);
    invmass_pp->Draw((invmass_pm ? "same" : ""));
  }
  if(invmass_mm) {
    invmass_mm->SetLineColor(4);
    invmass_mm->Draw((invmass_mm ? "same" : ""));
  }
  if(invmass_pm || invmass_mm || invmass_pp)
    latex->DrawLatex(0.12, 0.85, "J/#psi candidates");
  if(invmass_pm && zDistrib && zDistrib->Integral()>1.);
    latex->DrawLatex(0.12, 0.80, Form("candidates / event = %.3f", invmass_pm->Integral() / zDistrib->Integral()));
  TLegend* legend=new TLegend(0.7,0.7,0.89,0.89);
  legend->SetTextFont(42);
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  if(invmass_pm) legend->AddEntry(invmass_pm, "+- pairs", "l");
  if(invmass_pp) legend->AddEntry(invmass_pp, "++ pairs", "l");
  if(invmass_mm) legend->AddEntry(invmass_mm, "-- pairs", "l");
  if(invmass_mm || invmass_pm || invmass_pp) legend->Draw();
  
  c1->SaveAs(Form("fig_dq_tpcDedx.%s",suffix.Data()));
  c2->SaveAs(Form("fig_dq_tpcNsigmaElectron.%s",suffix.Data()));
  c3->SaveAs(Form("fig_dq_eeInvMass.%s",suffix.Data()));


  //  TFile* save=new TFile(outputName, "RECREATE");
  // c1->Write();
  // c2->Write();
  // c3->Write();
  // save->Close();


  // Added by jsatya

  TFile *fout = TFile::Open(outfile,"UPDATE");
  fout->ls();
  
  TDirectoryFile *cdd = NULL;
  cdd = (TDirectoryFile*)fout->Get("DQ");
  if(!cdd) {
    Printf("Warning: DQ <dir> doesn't exist, creating a new one");
    cdd = (TDirectoryFile*)fout->mkdir("DQ");
  }
  cdd->cd();
  cdd->ls();



  if (invmass_pp){
    invmass_pp->SetName(Form("fig_dq_%s_pp", invmass_pp->GetName()));
    invmass_pp->Write();
  }
  
  if (invmass_pm){
    invmass_pm->SetName(Form("fig_dq_%s_pm", invmass_pm->GetName()));
    invmass_pm->Write();
  }
  if (invmass_mm){
    invmass_mm->SetName(Form("fig_dq_%s_mm", invmass_mm->GetName()));
    invmass_mm->Write();
  }
  
  tpcNsigmaEle->SetName(Form("fig_dq_%s", tpcNsigmaEle->GetName()));
  tpcNsigmaEle->Write();
  
  tpcDedx->SetName(Form("fig_dq_%s", tpcDedx->GetName()));
  tpcDedx->Write();
  
  fout->cd();
  fout->Close();
  
}


//______________________________________________________________________________
TObject* GetHistogram(TFile* file, 
		      const Char_t* cutClass, 
		      const Char_t* histList, 
		      const Char_t* histName) {
  //
  // get a histogram from the J/psi output
  //
  TKey* qaKey = (TKey*)file->GetListOfKeys()->At(0);
  TList* qaList = (TList*)qaKey->ReadObj();
  
  THashList* qaCutClass = (THashList*)qaList->FindObject(cutClass);
  THashList* qaHistList = (THashList*)qaCutClass->FindObject(histList);
  return qaHistList->FindObject(histName);
}
