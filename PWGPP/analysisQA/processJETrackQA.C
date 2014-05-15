/********************************************************* 
   processJETrackQA: 
   Post processing macro for the JET Track QA


 *********************************************************/


const Float_t ptmin =  0.0 ;
void processJETrackQA(TString strFileIn   = "AnalysisResults.root", 
		      TString suffix      = "eps",
		      Int_t cent          = 10, 
		      Int_t trig          = 1, 
		      Bool_t bESD         = kFALSE, 
		      Int_t run           = 0, 
		      const char *outfile ="JETrackQA_output.root") {

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TString prefix = "fig_je_TrackQA_";

  TString strTrigger = "";
  if(trig==1) strTrigger = "kCentral";
  if(trig==2) strTrigger = "kSemiCentral";
  if(trig==3) strTrigger = "kMBkCentralkSemiCentral";
  if(trig==4) strTrigger = "kEMCEJE";
  if(trig==5) strTrigger = "kINT7";
  if(trig==6) strTrigger = "kMB";

  TString strTrigger2 = "";
  if(trig==1) strTrigger2 = "kCentral";
  if(trig==2) strTrigger2 = "kSemiCentral";
  if(trig==3) strTrigger2 = "kMBkCentralSemiCentral";
  if(trig==4) strTrigger2 = "kEMCEJE";
  if(trig==5) strTrigger2 = "kINT7";
  if(trig==6) strTrigger2 = "kMB";

  Int_t globStTrackType = 0;
  Int_t globStCuts = 5;
  Int_t globCnoSPDTrackType = 7;
  Int_t globCnoSPDCuts = 5;

  TFile * f1 = TFile::Open(strFileIn.Data());

  //Load histograms
  TList *histsGlobSt = 0x0;
  histsGlobSt = (TList*)f1->Get(Form("PWG4_HighPtTrackQACent%dTrackType%dCuts%d%s/qa_histsQAtrackCent%dType%dcuts%d%s",cent,globStTrackType,globStCuts,strTrigger.Data(),cent,globStTrackType,globStCuts,strTrigger2.Data()));

  if(!histsGlobSt) {
    Printf(">>>> ERROR: histsGlobSt is not Available");
    return;
  }

  TH1F *fNEventSelGlobSt = histsGlobSt->FindObject("fNEventSel");
  float NEventsGlobSt = fNEventSelGlobSt->GetEntries();

  TList *histsGlobCnoSPD = 0x0;
  histsGlobCnoSPD = (TList*)f1->Get(Form("PWG4_HighPtTrackQACent%dTrackType%dCuts%d%s/qa_histsQAtrackCent%dType%dcuts%d%s",cent,globCnoSPDTrackType,globCnoSPDCuts,strTrigger.Data(),cent,globCnoSPDTrackType,globCnoSPDCuts,strTrigger2.Data()));
  //histsGlobCnoSPD->Print();
  TH1F *fNEventSelGlobCnoSPD = histsGlobCnoSPD->FindObject("fNEventSel");
  float NEventsGlobCnoSPD = fNEventSelGlobCnoSPD->GetEntries();

  if(NEventsGlobSt==0)
    NEventsGlobSt=0.1;
  if(NEventsGlobCnoSPD==0)
    NEventsGlobCnoSPD=0.1;

  //---------------------------------------------------------------------------------------------------
  //                       phi distribution of hybrid tracks
  //---------------------------------------------------------------------------------------------------

  TH2F *fPtPhiGlobSt = histsGlobSt->FindObject("fPtPhi");
  TH2F *fPtPhiGlobCnoSPD = histsGlobCnoSPD->FindObject("fPtPhi");
  fPtPhiGlobSt->SetXTitle("p_{T,track} [GeV]");
  fPtPhiGlobSt->SetYTitle("#varphi");
  fPtPhiGlobCnoSPD->SetXTitle("p_{T,track} [GeV]");
  fPtPhiGlobCnoSPD->SetYTitle("#varphi");

  TCanvas *c2 =new TCanvas("c2","c2: Phi",600,450);
  Int_t binMin = 1;
  if(ptmin>0.) binMin = fPtPhiGlobSt->GetXaxis()->FindBin(ptmin+0.00001);
  TH1D *fPhiGlobSt = fPtPhiGlobSt->ProjectionY("fPhiGlobSt",binMin,fPtPhiGlobSt->GetNbinsX());
  TH1D *fPhiGlobCnoSPD = fPtPhiGlobCnoSPD->ProjectionY("fPhiGlobCnoSPD",binMin,fPtPhiGlobSt->GetNbinsX());

  fPhiGlobSt->SetLineColor(2);
  fPhiGlobSt->SetLineWidth(3);

  fPhiGlobCnoSPD->SetLineStyle(1);
  fPhiGlobCnoSPD->SetLineColor(4);
  fPhiGlobCnoSPD->SetLineWidth(3);

  TH1D *fPhiGlobSum = fPhiGlobSt->Clone();
  fPhiGlobSum->SetTitle("fPhiGlobSum");
  fPhiGlobSum->SetName("fPhiGlobSum");
  fPhiGlobSum->Add(fPhiGlobCnoSPD);

  fPhiGlobSt->Scale(1./NEventsGlobSt,"width");
  fPhiGlobCnoSPD->Scale(1./NEventsGlobSt,"width");
  fPhiGlobSum->Scale(1./NEventsGlobSt,"width");

  gPad->SetLeftMargin(0.2);
  gPad->SetRightMargin(0.02);
  gPad->SetBottomMargin(0.12);

  TH1F *frame2 = gPad->DrawFrame(0.,0.,2.*TMath::Pi(),fPhiGlobSum->GetBinContent(fPhiGlobSum->GetMaximumBin())*1.5);
  frame2->SetXTitle("#varphi");
  frame2->SetYTitle("#frac{1}{N_{evts}} #frac{dN}{d#varphi}");
  frame2->GetXaxis()->SetTitleSize(0.06);
  frame2->GetYaxis()->SetTitleSize(0.06);
  frame2->GetXaxis()->SetTitleOffset(0.75);
  frame2->GetYaxis()->SetTitleOffset(1.4);

  fPhiGlobSt->DrawCopy("same");
  fPhiGlobCnoSPD->DrawCopy("same");
 
  fPhiGlobSum->SetLineColor(1);
  fPhiGlobSum->SetMarkerColor(1);
  fPhiGlobSum->DrawCopy("same");

  TLegend *leg2 = NULL;
  if(run>0) leg2 = new TLegend(0.22,0.65,0.88,0.88,Form("Hybrid tracks. run:%d",run));
  else leg2 = new TLegend(0.22,0.65,0.88,0.88,"Hybrid tracks");
  leg2->SetFillColor(10);
  leg2->SetBorderSize(0);
  leg2->SetFillStyle(0);
  leg2->SetTextSize(0.06);
  leg2->AddEntry(fPhiGlobSt,"w/ SPD & ITSrefit","l");
  leg2->AddEntry(fPhiGlobCnoSPD,"w/o SPD & w/ ITSrefit","l");
  leg2->AddEntry(fPhiGlobSum,"sum","l");
  leg2->Draw();

  TLatex textNEvents;
  textNEvents.SetNDC();
  textNEvents.DrawLatex(0.52,0.42,Form("#it{N}_{events} = %.0f",NEventsGlobSt));

  c2->SaveAs(Form("%sPhiCent%d%sRun%d.%s",prefix.Data(),cent,strTrigger.Data(),run,suffix.Data()));

  //---------------------------------------------------------------------------------------------------
  //                       pt distribution of hybrid tracks
  //---------------------------------------------------------------------------------------------------

  TCanvas *c1 =new TCanvas("c1","c1: Phi",600,450);
  Int_t binMin = 1;
  if(ptmin>0.) binMin = fPtPhiGlobSt->GetXaxis()->FindBin(ptmin+0.00001);
  TH1D *fPtGlobSt = fPtPhiGlobSt->ProjectionX("fPtGlobSt");
  TH1D *fPtGlobCnoSPD = fPtPhiGlobCnoSPD->ProjectionX("fPtGlobCnoSPD");

  fPtGlobSt->SetLineColor(2);
  fPtGlobSt->SetLineWidth(3);

  fPtGlobCnoSPD->SetLineStyle(1);
  fPtGlobCnoSPD->SetLineColor(4);
  fPtGlobCnoSPD->SetLineWidth(3);

  TH1D *fPtGlobSum = fPtGlobSt->Clone();
  fPtGlobSum->SetTitle("fPtGlobSum");
  fPtGlobSum->SetName("fPtGlobSum");
  fPtGlobSum->Add(fPtGlobCnoSPD);

  fPtGlobSt->Scale(1./NEventsGlobSt,"width");
  fPtGlobCnoSPD->Scale(1./NEventsGlobSt,"width");
  fPtGlobSum->Scale(1./NEventsGlobSt,"width");

  gPad->SetLeftMargin(0.2);
  gPad->SetRightMargin(0.02);
  gPad->SetBottomMargin(0.12);

  TH1F *frame1 = gPad->DrawFrame(0.,1e-7,100.,fPtGlobSum->GetBinContent(fPtGlobSum->GetMaximumBin())*1.5);
  frame1->SetXTitle("p_{T} [GeV]");
  frame1->SetYTitle("#frac{1}{N_{evts}}#frac{dN}{dp_{T}} [GeV^{-1}]");
  frame1->GetXaxis()->SetTitleSize(0.06);
  frame1->GetYaxis()->SetTitleSize(0.06);
  frame1->GetXaxis()->SetTitleOffset(0.75);
  frame1->GetYaxis()->SetTitleOffset(1.4);

  gPad->SetLogy();

  fPtGlobSt->DrawCopy("same");
  fPtGlobCnoSPD->DrawCopy("same");
 
  fPtGlobSum->SetLineColor(1);
  fPtGlobSum->SetMarkerColor(1);
  fPtGlobSum->DrawCopy("same");

  TLegend *leg1 = NULL;
  if(run>0) leg1 = new TLegend(0.35,0.65,0.88,0.88,Form("Hybrid tracks. run:%d",run));
  else leg1 = new TLegend(0.35,0.65,0.88,0.88,"Hybrid tracks");
  leg1->SetFillColor(10);
  leg1->SetBorderSize(0);
  leg1->SetFillStyle(0);
  leg1->SetTextSize(0.06);
  leg1->AddEntry(fPtGlobSt,"w/ SPD & ITSrefit","l");
  leg1->AddEntry(fPtGlobCnoSPD,"w/o SPD & w/ ITSrefit","l");
  leg1->AddEntry(fPtGlobSum,"sum","l");
  leg1->Draw();

  TLatex textNEvents;
  textNEvents.SetNDC();
  textNEvents.DrawLatex(0.55,0.45,Form("#it{N}_{events} = %.0f",NEventsGlobSt));

  c1->SaveAs(Form("%sPtCent%d%sRun%d.%s",prefix.Data(),cent,strTrigger.Data(),run,suffix.Data()));

  //---------------------------------------------------------------------------------------------------
  //                       pT resolution
  //---------------------------------------------------------------------------------------------------


  TProfile *fProfPtPtSigma1PtGlobSt = histsGlobSt->FindObject("fProfPtPtSigma1Pt");
  TProfile *fProfPtPtSigma1PtGlobCnoSPD = histsGlobCnoSPD->FindObject("fProfPtPtSigma1Pt");

  fProfPtPtSigma1PtGlobSt->SetTitle("fProfPtPtSigma1PtGlobSt");
  fProfPtPtSigma1PtGlobSt->SetName("fProfPtPtSigma1PtGlobSt");
  fProfPtPtSigma1PtGlobCnoSPD->SetTitle("fProfPtPtSigma1PtGlobCnoSPD");
  fProfPtPtSigma1PtGlobCnoSPD->SetName("fProfPtPtSigma1PtGlobCnoSPD");

  for(Int_t i =1 ; i<fProfPtPtSigma1PtGlobSt->GetNbinsX(); i++) {
    if(fProfPtPtSigma1PtGlobSt->GetBinEffectiveEntries(i)<10.) {
      fProfPtPtSigma1PtGlobSt->SetBinContent(i,0);
      fProfPtPtSigma1PtGlobSt->SetBinError(i,0);
    }
    if(fProfPtPtSigma1PtGlobCnoSPD->GetBinEffectiveEntries(i)<10.) {
      fProfPtPtSigma1PtGlobCnoSPD->SetBinContent(i,0);
      fProfPtPtSigma1PtGlobCnoSPD->SetBinError(i,0);
    }
  }

  TCanvas *c3 =new TCanvas("c3","c3: pT resolution",600,450);
  gPad->SetLeftMargin(0.12);
  gPad->SetRightMargin(0.02);
  gPad->SetBottomMargin(0.12);
  TH1F *frame3 = gPad->DrawFrame(0.,0.,100.,0.3);
  frame3->SetXTitle("#it{p}_{T,track} [GeV/c]");
  frame3->SetYTitle("#it{p}_{T,track}#sigma(1/#it{p}_{T,track})");
  frame3->GetXaxis()->SetTitleSize(0.06);
  frame3->GetYaxis()->SetTitleSize(0.06);
  frame3->GetXaxis()->SetTitleOffset(0.8);
  frame3->GetYaxis()->SetTitleOffset(0.8);

  fProfPtPtSigma1PtGlobSt->SetLineColor(2);
  fProfPtPtSigma1PtGlobSt->SetLineWidth(3);
  fProfPtPtSigma1PtGlobSt->SetMarkerStyle(21);
  fProfPtPtSigma1PtGlobSt->SetMarkerColor(2);

  fProfPtPtSigma1PtGlobCnoSPD->SetLineColor(4);
  fProfPtPtSigma1PtGlobCnoSPD->SetLineWidth(3);
  fProfPtPtSigma1PtGlobCnoSPD->SetMarkerStyle(24);
  fProfPtPtSigma1PtGlobCnoSPD->SetMarkerColor(4);

  fProfPtPtSigma1PtGlobSt->DrawCopy("same");
  fProfPtPtSigma1PtGlobCnoSPD->DrawCopy("same");

  TLegend *leg3 = NULL;
  if(run>0) leg3 = new TLegend(0.16,0.6,0.88,0.88,Form("Hybrid tracks. run:%d",run));
  else leg3 = new TLegend(0.16,0.6,0.88,0.88,"Hybrid tracks");
  leg3->SetFillColor(10);
  leg3->SetBorderSize(0);
  leg3->SetFillStyle(0);
  leg3->SetTextSize(0.06);
  leg3->AddEntry(fProfPtPtSigma1PtGlobSt,"w/ SPD & ITSrefit","lp");
  leg3->AddEntry(fProfPtPtSigma1PtGlobCnoSPD,"w/o SPD & w/ ITSrefit","lp");
  leg3->Draw();

  c3->SaveAs(Form("%sPtResolutionCent%d%sRun%d.%s",prefix.Data(),cent,strTrigger.Data(),run,suffix.Data()));

  //---------------------------------------------------------------------------------------------------
  //                       WRITE OUTPUT TO ROOT FILE
  //---------------------------------------------------------------------------------------------------

  // Modifed by satya to have unique name of objects and same name of file
  //
  //  TFile *histOut = new TFile(Form("%sHybridCent%d%sRun%d.root",prefix.Data(),cent,strTrigger.Data(),run),"RECREATE");

 

  // Added by sjena
  TFile *fout = TFile::Open(outfile,"UPDATE");
  fout->ls();
  
  TDirectoryFile *cdd = NULL;
  cdd = (TDirectoryFile*)fout->Get("JE");
  if(!cdd) {
    Printf("Warning: JE <dir> doesn't exist, creating a new one");
    cdd = (TDirectoryFile*)fout->mkdir("JE");
  }
  cdd->cd();
  cdd->ls();

  fPtGlobSt->SetName(Form("%s%s",prefix.Data(), fPtGlobSt->GetName()));
  fPtGlobSt->Write();
  fPtGlobCnoSPD->SetName(Form("%s%s",prefix.Data(), fPtGlobCnoSPD->GetName()));
  fPtGlobCnoSPD->Write();
  fPtGlobSum->SetName(Form("%s%s",prefix.Data(), fPtGlobSum->GetName()));
  fPtGlobSum->Write();
  fPhiGlobSt->SetName(Form("%s%s",prefix.Data(), fPhiGlobSt->GetName()));
  fPhiGlobSt->Write();
  fPhiGlobCnoSPD->SetName(Form("%s%s",prefix.Data(), fPhiGlobCnoSPD->GetName()));
  fPhiGlobCnoSPD->Write();
  fPhiGlobSum->SetName(Form("%s%s",prefix.Data(), fPhiGlobSum->GetName()));
  fPhiGlobSum->Write();
  fProfPtPtSigma1PtGlobSt->SetName(Form("%s%s",prefix.Data(), fProfPtPtSigma1PtGlobSt->GetName()));
  fProfPtPtSigma1PtGlobSt->Write();
  fProfPtPtSigma1PtGlobCnoSPD->SetName(Form("%s%s",prefix.Data(), fProfPtPtSigma1PtGlobCnoSPD->GetName()));
  fProfPtPtSigma1PtGlobCnoSPD->Write();

  fout->Close();


}
