/***************************************************************
 processHFEQAtask:

 Post Processing of HF QA task in Analysis QA train 
 
 Modification Done // sjena
 To save file and putting an unique naming convention

***************************************************************/

void processHFEQAtask(const char *fnamedata = "AnalysisResults.root", 
		      TString suffix = "eps",
		      const char * outfile = "HFEQAtask_output.root") {

  const char *fnamedir  = "t131073f16s0p1TPC110r60p80ITS4Pi2DCAr100z200TOF30TPCe50V0D2er8i0t-20t50";
  const char *fnamelist = "list_t131073f16s0p1TPC110r60p80ITS4Pi2DCAr100z200TOF30TPCe50V0D2er8i0t-20t50";


  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libCORRFW.so");
  gSystem->Load("libTENDER.so");
  gSystem->Load("libPWGPP.so");
  gSystem->Load("libPWGflowBase.so");
  gSystem->Load("libPWGflowTasks.so");
  gSystem->Load("libPWGHFhfe.so");

  
  gStyle->SetPalette(1);
  gStyle->SetOptStat(1111);
  gStyle->SetPadBorderMode(0);
  gStyle->SetCanvasColor(10);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadRightMargin(0.13);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadTopMargin(0.13);

  
  //////////////////////
  // Data
  //////////////////////
  TFile *indata = TFile::Open(fnamedata);
  indata->ls();
  TDirectoryFile *d = (TDirectoryFile *) indata->Get(fnamedir);

  if(!d) {
    Printf("FATAL: %s :  Not found", fnamedir);
    Printf("Exiting <<<<<<<<<<<<<<<< processHFEQAtask >>>>>>>>>>>>>");
    return;
  }

  TList *qadata = (TList *) d->Get(fnamelist);

  
  if(!qadata) {
    Printf("FATAL: %s :  Not found", fnamelist);
    Printf("Exiting <<<<<<<<<<<<<<<< processHFEQAtask >>>>>>>>>>>>>");
    return;
  }



  TList *qaadata = (TList *)qadata->FindObject("HFEpidQA");
 
  if(!qaadata) {
    Printf("FATAL: HFEpidQA   Not found");
    Printf("Exiting <<<<<<<<<<<<<<<< processHFEQAtask >>>>>>>>>>>>>");
    return;
  }


  // Make Plots for TPC
  AliHFEtpcPIDqa *tpcqadata = (AliHFEtpcPIDqa *)qaadata->FindObject("TPCQA");
  AliHFEcollection *collectiontpc = tpcqadata->GetHistograms();
  THnSparseF *tpcsparseF = (THnSparseF *) collectiontpc->Get("tpcnSigma");
 
  // Make Plots for TOF
  AliHFEtofPIDqa *tofqadata = (AliHFEtofPIDqa *)qaadata->FindObject("TOFQA");
  AliHFEcollection *collectiontof = tofqadata->GetHistoCollection();
  THnSparseF *tofsparseF = (THnSparseF *) collectiontof->Get("tofnSigma");
  THnSparseF *toftpcsparseF = (THnSparseF *) collectiontof->Get("tofMonitorTPC");
  

  ////////////////
  // Projection
  /////////////////

  TAxis *pidaxistpc = tpcsparseF->GetAxis(0);
  TAxis *pidaxistof = tofsparseF->GetAxis(0);
  TAxis *pidaxistoftpc = toftpcsparseF->GetAxis(0);

  TAxis *stepaxistpc = tpcsparseF->GetAxis(3);
  TAxis *stepaxistof = tofsparseF->GetAxis(3);
  TAxis *stepaxistoftpc = toftpcsparseF->GetAxis(3);

  TAxis *centralityaxistpc = tpcsparseF->GetAxis(4);
  TAxis *centralityaxistof = tofsparseF->GetAxis(4);
  TAxis *centralityaxistoftpc = toftpcsparseF->GetAxis(4);

  stepaxistpc->SetRange(1,1);
  stepaxistof->SetRange(1,1);
  stepaxistoftpc->SetRange(2,2);

  Int_t nbinsc = centralityaxistpc->GetNbins();
  printf("There are %d centrality bins \n",nbinsc);
  

  //
  centralityaxistpc->SetRange(8,9);
  centralityaxistof->SetRange(8,9);
  centralityaxistoftpc->SetRange(8,9);
  stepaxistoftpc->SetRange(1,1);  
  TH2D *toftpc2Dsumper = toftpcsparseF->Projection(2,1);
  toftpc2Dsumper->SetName("toftpc2Dsumper_21"); 
  TH2D *tof2Dsumper = tofsparseF->Projection(2,1);
  tof2Dsumper->SetName("tof2Dsumper_21");

  stepaxistof->SetRange(2,2);
  TH2D *tof2Dsumbper = tofsparseF->Projection(2,1);
  tof2Dsumbper->SetName("tof2Dsumbper_21");

  TH2D *tpc2Dsumper = tpcsparseF->Projection(2,1);
  tpc2Dsumper->SetName("toftpc2Dsumper_21");

  stepaxistpc->SetRange(2,2);
  TH2D *tpc2Dsumbper = tpcsparseF->Projection(2,1);
  tpc2Dsumbper->SetName("tpc2Dsumbper_21");

  //
  centralityaxistpc->SetRange(0,1);
  centralityaxistof->SetRange(0,1);
  centralityaxistoftpc->SetRange(0,1);
  stepaxistpc->SetRange(1,1);
  stepaxistof->SetRange(1,1);
  stepaxistoftpc->SetRange(1,1);
  TH2D *toftpc2Dsumc = toftpcsparseF->Projection(2,1);
  toftpc2Dsumc->SetName("toftpc2Dsumc_21");

  TH2D *tof2Dsumc = tofsparseF->Projection(2,1);
  tof2Dsumc->SetName("tof2Dsumc_21");

  stepaxistof->SetRange(2,2);

  TH2D *tof2Dsumbc = tofsparseF->Projection(2,1);
  tof2Dsumbc->SetName("tof2Dsumbc_21");

  TH2D *tpc2Dsumc = tpcsparseF->Projection(2,1);
  tpc2Dsumc->SetName("tpc2Dsumc_21");

  stepaxistpc->SetRange(2,2);
  TH2D *tpc2Dsumbc = tpcsparseF->Projection(2,1);
  tpc2Dsumbc->SetName("tpc2Dsumbc_21");

  
  ///////////////
  // Plots
  ///////////////

  //added by sjena 
  TFile *fout = TFile::Open(outfile,"UPDATE");
  fout->ls();
  
  TDirectoryFile *cdd = NULL;
  cdd = (TDirectoryFile*)fout->Get("HF");
  if(!cdd) {
    Printf("Warning: HF <dir> doesn't exist, creating a new one");
    cdd = (TDirectoryFile*)fout->mkdir("HF");
  }
  cdd->cd();
  cdd->ls();



  // Draw Plots
  TCanvas *cElectronsTPCper = new TCanvas("ElectronTPCPID_70-90%", "ElectronTPCPID_70-90%");
  cElectronsTPCper->Divide(3,1);
  cElectronsTPCper->cd(1);
  gPad->SetLogz();
  gPad->SetLogx();
  gPad->SetTicks();
  tof2Dsumper->SetStats(0);
  tof2Dsumper->SetTitle("TOF 70-90% Pb-Pb");
  tof2Dsumper->GetXaxis()->SetTitle("p_{t} (GeV/c)");
  tof2Dsumper->GetYaxis()->SetTitle("TOF time_{e} - expected time|_{e}  (#sigma)");
  tof2Dsumper->GetXaxis()->SetRangeUser(0.3,10.0);
  tof2Dsumper->GetXaxis()->SetTitleSize(0.05);
  tof2Dsumper->GetYaxis()->SetTitleSize(0.05);
  tof2Dsumper->Draw("colz");
  
  //Added By Satya
  tof2Dsumper->SetName(Form("fig_hf_1_%s", tof2Dsumper->GetName()));
  tof2Dsumper->Write();

  cElectronsTPCper->cd(2);
  gPad->SetLogz();
  gPad->SetLogx();
  gPad->SetTicks();
  tof2Dsumbper->SetStats(0);
  tof2Dsumbper->SetTitle("TOF cut 70-90% Pb-Pb");
  tof2Dsumbper->GetXaxis()->SetTitle("p_{t} (GeV/c)");
  tof2Dsumbper->GetYaxis()->SetTitle("TOF time_{e} - expected time|_{e}  (#sigma)");
  tof2Dsumbper->GetXaxis()->SetRangeUser(0.3,10.0);
  tof2Dsumbper->GetXaxis()->SetTitleSize(0.05);
  tof2Dsumbper->GetYaxis()->SetTitleSize(0.05);
  tof2Dsumbper->Draw("colz");


  //Added By Satya
  tof2Dsumbper->SetName(Form("fig_hf_2_%s", tof2Dsumbper->GetName()));
  tof2Dsumbper->Write();

  cElectronsTPCper->cd(3);
  gPad->SetLogz();
  gPad->SetLogx();
  gPad->SetTicks();
  tpc2Dsumper->SetStats(0);
  tpc2Dsumper->SetTitle("TPC 70-90% Pb-Pb");
  tpc2Dsumper->GetXaxis()->SetTitle("p_{t} (GeV/c)");
  tpc2Dsumper->GetYaxis()->SetTitle("TPC dE/dx - expected dE/dx|_{e}  (#sigma)");
  tpc2Dsumper->GetXaxis()->SetRangeUser(0.3,10.0);
  tpc2Dsumper->GetXaxis()->SetTitleSize(0.05);
  tpc2Dsumper->GetYaxis()->SetTitleSize(0.05);
  tpc2Dsumper->Draw("colz");
  //Added by satya

  
  tpc2Dsumper->SetName(Form("fig_hf_3_%s", tpc2Dsumper->GetName()));
  tpc2Dsumper->Write();

  cElectronsTPCper->SaveAs(Form("fig_hf_ElectronTPCPID_70-90.%s",suffix.Data()));



  TCanvas *cElectronsTPCc = new TCanvas("ElectronTPCPID_0-10%", "ElectronTPCPID_0-10%");
  cElectronsTPCc->Divide(3,1);
  cElectronsTPCc->cd(1);
  gPad->SetLogz();
  gPad->SetLogx();
  gPad->SetTicks();
  tof2Dsumc->SetStats(0);
  tof2Dsumc->SetTitle("TOF 0-10% Pb-Pb");
  tof2Dsumc->GetXaxis()->SetTitle("p_{t} (GeV/c)");
  tof2Dsumc->GetYaxis()->SetTitle("TOF time_{e} - expected time|_{e}  (#sigma)");
  tof2Dsumc->GetXaxis()->SetRangeUser(0.3,10.0);
  tof2Dsumc->GetXaxis()->SetTitleSize(0.05);
  tof2Dsumc->GetYaxis()->SetTitleSize(0.05);
  tof2Dsumc->Draw("colz");

  tof2Dsumc->SetName(Form("fig_hf_4_%s", tof2Dsumc->GetName()));
  tof2Dsumc->Write();

  cElectronsTPCc->cd(2);
  gPad->SetLogz();
  gPad->SetLogx();
  gPad->SetTicks();
  tof2Dsumbc->SetStats(0);
  tof2Dsumbc->SetTitle("TOF cut 0-10% Pb-Pb");
  tof2Dsumbc->GetXaxis()->SetTitle("p_{t} (GeV/c)");
  tof2Dsumbc->GetYaxis()->SetTitle("TOF time_{e} - expected time|_{e}  (#sigma)");
  tof2Dsumbc->GetXaxis()->SetRangeUser(0.3,10.0);
  tof2Dsumbc->GetXaxis()->SetTitleSize(0.05);
  tof2Dsumbc->GetYaxis()->SetTitleSize(0.05);
  tof2Dsumbc->Draw("colz");

  tof2Dsumbc->SetName(Form("fig_hf_5_%s", tof2Dsumbc->GetName()));
  tof2Dsumbc->Write();

  cElectronsTPCc->cd(3);
  gPad->SetLogz();
  gPad->SetLogx();
  gPad->SetTicks();
  tpc2Dsumc->SetStats(0);
  tpc2Dsumc->SetTitle("TPC 0-10% Pb-Pb");
  tpc2Dsumc->GetXaxis()->SetTitle("p_{t} (GeV/c)");
  tpc2Dsumc->GetYaxis()->SetTitle("TPC dE/dx - expected dE/dx|_{e}  (#sigma)");
  tpc2Dsumc->GetXaxis()->SetRangeUser(0.3,10.0);
  tpc2Dsumc->GetXaxis()->SetTitleSize(0.05);
  tpc2Dsumc->GetYaxis()->SetTitleSize(0.05);
  tpc2Dsumc->Draw("colz");
  // Added by satya

  tpc2Dsumc->SetName(Form("fig_hf_6_%s", tpc2Dsumc->GetName()));
  tpc2Dsumc->Write();

  cElectronsTPCc->SaveAs(Form("fig_hf_ElectronTPCPID_0-10.%s",suffix.Data()));

  fout->Close();


}
