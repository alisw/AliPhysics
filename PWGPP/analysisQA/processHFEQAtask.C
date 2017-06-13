/***************************************************************
 processHFEQAtask:

 Post Processing of HF QA task in Analysis QA train 
 
 Modification  by Noor Alam on 4.5.2017

***************************************************************/

void processHFEQAtask(TString fnamedata = "AnalysisResults.root", 
		      TString suffix = "eps",
		      int type=1) {


  TString fnamedir  = "PWGHF_hfeHFEemc";
  TString QAINT8=fnamedir+"QAINT8_EMC";
  TString QAINT7=fnamedir+"QAINT7_EMC";
  TString QATrigGAEG1=fnamedir+"QATrigGAEG1_EMC";
  TString QATrigGAEG2=fnamedir+"QATrigGAEG2_EMC";
  TString QATrigGADG1=fnamedir+"QATrigGADG1_EMC";
  TString QATrigGADG2=fnamedir+"QATrigGADG2_EMC";


  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCORRFW");
  gSystem->Load("libTender");
  gSystem->Load("libPWGPP");
  gSystem->Load("libPWGflowBase");
  gSystem->Load("libPWGflowTasks");
  gSystem->Load("libPWGHFhfe");

  
  gStyle->SetPalette(1);
  gStyle->SetOptStat(1111);
  gStyle->SetPadBorderMode(0);
  gStyle->SetCanvasColor(10);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadRightMargin(0.13);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadTopMargin(0.13);

  
  TFile *qaEvent = TFile::Open(fnamedata.Data());
//  qaEvent->ls();


TString dirName;

if(type == 1) dirName=QAINT7;
else if(type == 2) dirName=QATrigGAEG1;
else if(type == 3) dirName=QATrigGAEG2;

  TDirectoryFile *d = (TDirectoryFile *) qaEvent->Get(Form("%s",dirName.Data()));

  if(!d) {
    Printf("FATAL: %s :  Not found", fnamedata.Data());
    Printf("Exiting <<<<<<<<<<<<<<<< processHFEQAtask >>>>>>>>>>>>>");
    return;
  }

  cout<<" dir name : "<<dirName<<'\t'<<"Tlist Name : "<<d->GetListOfKeys()->At(0)->GetName()<<endl;

  TList *qadata = (TList *) (d->Get(d->GetListOfKeys()->At(0)->GetName()));

qadata->ls();

// now it is time to play with histogram 

  TCanvas *c1 = new TCanvas();
  c1->Divide(3,2);
  c1->cd(1);
  gPad->SetLogy();
  TH1F *ClusE = qadata->FindObject("fHistClustE");
  ClusE->Draw();

  c1->cd(2);
  TH2F *EMCClsEtaPhi= qadata->FindObject("fEMCClsEtaPhi");
  EMCClsEtaPhi->Draw("colz");

  c1->cd(3);
  gPad->SetLogz();
  TH2F *EMCDeltaR = qadata->FindObject("fEMCTrkMatch");
  EMCDeltaR->Draw("colz");


  c1->cd(4);
  gPad->SetLogz();
  TH2F *NsigEovP = qadata->FindObject("fHistNsigEop");
  NsigEovP->Draw("colz");

  c1->cd(5);
  TH1D* EovP = NsigEovP->ProjectionX("EovP",110,114);
  EovP->SetTitle("E/p;E/p;counts"); 
  EovP->Draw(); 

  c1->cd(6);
  gPad->SetLogy();
  TH1F *fULS = qadata->FindObject("fInvmassULS");
  TH1F *fLS = qadata->FindObject("fInvmassLS");
  
  fULS->SetAxisRange(0,0.3);
  fULS->Draw();
  fLS->Draw("same");

}
