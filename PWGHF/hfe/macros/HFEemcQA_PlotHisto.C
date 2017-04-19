#include "TLegend.h"

void SetTrigQA(char *name0, char *name1);

////////////////////////////////////////////////////////////
// Macro for plotting histos from AliAnalysisTaskHFEemcQA //
//            Author: Deepa Thomas, UT Austin             //
////////////////////////////////////////////////////////////


void HFEemcQA_PlotHisto()
{
  gROOT->Reset();
  gStyle->SetPalette(1);
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetCanvasBorderMode(0);     
  gStyle->SetPadBorderMode(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetOptStat(0);

  
  char dQA0[100];
  char dQA1[100];
  Bool_t iMC = kFALSE;
  SetTrigQA(dQA0,dQA1,iMC);
  cout << "MC = " << iMC << endl;

  TString file;
  file.Form("%s.pdf",dQA1);
  //file.Form("EMC.pdf");

  TFile *f1 = new TFile("AnalysisResults.root");
  //PWGHF_hfeHFEemcQAINT7 - QA1
  //PWGHF_hfeHFEemcQAINT8 - HFEemcQAINT8_woTender
  //PWGHF_hfeHFEemcQAEMC7 - HFEemcQAEMC7_woTender
  //PWGHF_hfeHFEemcQAEMC8 - HFEemcQAEMC8_woTender
  //PWGHF_hfeHFEemcQATrigGA - HFEemcQATrigGA_woTender
  //PWGHF_hfeHFEemcQATrigJE - HFEemcQATrigJE_woTender

   TDirectory *QA0 = (TDirectory*)f1->Get(dQA0);
   TList *QA1 = (TList*)QA0->Get(dQA1);

  TH1F *NEvents = QA1->FindObject("fNevents");
  TH1F *VtxZ = QA1->FindObject("fVtxZ");
  TH2F *TrigMulti = QA1->FindObject("fTrigMulti");
  TH1F *ClusE = QA1->FindObject("fHistClustE");
  TH2F *ClusEcent = QA1->FindObject("fHistClustEcent");
  TH2F *EMCClsEtaPhi= QA1->FindObject("fEMCClsEtaPhi");
  TH1F *NCls = QA1->FindObject("fHistoNCls");
  TH1F *NClsE1 = QA1->FindObject("fHistoNClsE1");
  TH1F *NClsE2 = QA1->FindObject("fHistoNClsE2");
  TH1F *NClsE3 = QA1->FindObject("fHistoNClsE3");
  TH2F *NCellPerCls = QA1->FindObject("fHistoNCells");
  TH2F *CellE = QA1->FindObject("fHistoCalCell");
  TH1F *TrkPt = QA1->FindObject("fTrkPt");
  TH1F *Trketa= QA1->FindObject("fTrketa");
  TH1F *Trkphi= QA1->FindObject("fTrkphi");
  TH2F *dEdx = QA1->FindObject("fdEdx");
  TH2F *TPCNsig = QA1->FindObject("fTPCnsig");
  TH2F *TPCNsigEta0 = QA1->FindObject("fTPCnsigEta0");
  TH2F *TPCNsigEta1 = QA1->FindObject("fTPCnsigEta1");
  TH2F *TPCNsigEta2 = QA1->FindObject("fTPCnsigEta2");
  TH1F *EMCTrkMatchPt = QA1->FindObject("fHistPtMatch");
  TH2F *EMCDeltaR = QA1->FindObject("fEMCTrkMatch");
  //TH2F *EMCDeltaR = QA1->FindObject("fEMCTrkMatch2");
  TH1F *EMCTrkPt = QA1->FindObject("fEMCTrkPt");
  TH1F *EMCTrketa = QA1->FindObject("fEMCTrketa");
  TH1F *EMCTrkphi = QA1->FindObject("fEMCTrkphi");
  TH2F *ClsEAftTrkMatch = QA1->FindObject("fClsEtaPhiAftMatch");
  TH2F *EMCTrkTPCNsig = QA1->FindObject("fEMCTPCnsig");
  //TH2F *EovP = QA1->FindObject("fHistEop");
  //TH2F *NsigEovP = QA1->FindObject("fHistdEdxEop");
  TH2F *NsigEovP = QA1->FindObject("fHistNsigEop");
  TH2F *M20 = QA1->FindObject("fM20");
  TH2F *M02 = QA1->FindObject("fM02");
  TH2F *M02EovP = QA1->FindObject("fM02EovP");
  TH2F *M20EovP = QA1->FindObject("fM20EovP");
  THnSparse *electron = QA1->FindObject("Electron");  
  TH1F *cent = QA1->FindObject("fCent");
  TH1F *fULS = QA1->FindObject("fInvmassULS");
  TH1F *fLS = QA1->FindObject("fInvmassLS");
  TH2F *fEPV0 = QA1->FindObject("fEvPlaneV0");
  TH2F *fEPV0A = QA1->FindObject("fEvPlaneV0A");
  TH2F *fEPV0C = QA1->FindObject("fEvPlaneV0C");
  TH2F *fEPTPC = QA1->FindObject("fEvPlaneTPC");
  TH2F *fTrMult = QA1->FindObject("fMult"); 
  TH2F *fMomInfo = QA1->FindObject("fMCcheckMother");
  TH2F *fULSmc = QA1->FindObject("fInvmassULS_MCtrue");

  cout << "No of events with NTrk > 2, ZVtx < 10 cm = " << NEvents->GetBinContent(3) <<endl;

  TCanvas *NE = new TCanvas("NoEvents", "Number of events",50,50,700,500);
  ProcessHisto(NEvents);
  NEvents->Draw();
  NE->Print(file + "(","pdf");

  TCanvas *c1 = new TCanvas("Centrality", "Centrality distribution",50,50,700,500);
  ProcessHisto(cent);
  //VtxZ->GetXaxis()->SetRangeUser(-25,25);
  gPad->SetLogy();
  cent->Draw();
  c1->Print(file,"pdf");

  TCanvas *c2 = new TCanvas("VtxZ", "Z vertex postion",50,50,700,500);
  ProcessHisto(VtxZ);
  VtxZ->GetXaxis()->SetRangeUser(-25,25);
  VtxZ->Draw();
  c2->Print(file,"pdf");

  TCanvas *c3 = new TCanvas("TrigMultiplicity", "Trigger multiplicity",50,50,700,500);
  ProcessHisto2D(TrigMulti);
  TrigMulti->GetXaxis()->SetTitle("");
  TrigMulti->GetXaxis()->SetBinLabel(1,"All"); 
  TrigMulti->GetXaxis()->SetBinLabel(2,"kAny"); 
  TrigMulti->GetXaxis()->SetBinLabel(3,"kMB"); 
  TrigMulti->GetXaxis()->SetBinLabel(4,"kINT7"); 
  TrigMulti->GetXaxis()->SetBinLabel(5,"kINT8"); 
  TrigMulti->GetXaxis()->SetBinLabel(6,"kEMC1"); 
  TrigMulti->GetXaxis()->SetBinLabel(7,"kEMC7"); 
  TrigMulti->GetXaxis()->SetBinLabel(8,"kEMC8"); 
  TrigMulti->GetXaxis()->SetBinLabel(9,"kEMCEJE"); 
  TrigMulti->GetXaxis()->SetBinLabel(10,"kEMCEGA"); 
  TrigMulti->GetXaxis()->SetBinLabel(11,"kEMCEGA & EG2"); 
  TrigMulti->Draw("COLZ");
  gPad->SetLogz();
  //c3->Print(file,"pdf");

  /*
  TCanvas *c4 = new TCanvas();
  c4->Divide(2,2);
  c4->cd(1);
  TH1D *EPV0 = fEPV0->ProjectionY("EPV0",31,50);
  EPV0->SetMinimum(1);
  EPV0->Draw();
  c4->cd(2);
  TH1D *EPV0A = fEPV0A->ProjectionY("EPV0A",31,50);
  EPV0A->SetMinimum(1);
  EPV0A->Draw();
  c4->cd(3);
  TH1D *EPV0C = fEPV0C->ProjectionY("EPV0C",31,50);
  EPV0C->SetMinimum(1);
  EPV0C->Draw();
  c4->cd(4);
  TH1D *EPTPC = fEPTPC->ProjectionY("EPTPC",31,50);
  EPTPC->SetMinimum(1);
  EPTPC->Draw();
  c4->Print(file,"pdf");
  */

  TCanvas *c5 = new TCanvas();
  gPad->SetLogz();
  fTrMult->Draw("colz");
  c5->Print(file,"pdf");


  TCanvas *c10 = new TCanvas("EMCClusEtaPhi", "EMCAL cluster eta and phi distribution",50,50,700,500);
  ProcessHisto2D(EMCClsEtaPhi);
  EMCClsEtaPhi->GetXaxis()->SetRangeUser(-0.7,0.7);
  gPad->SetLogz();
  EMCClsEtaPhi->Draw("COLZ");
  c10->Print(file,"pdf");


  TCanvas *ncel = new TCanvas("NcellPerCls", "N cell per cluster",50,50,700,500);
  ProcessHisto2D(NCellPerCls);
  gPad->SetLogz();
  NCellPerCls->Draw("COLZ");
  ncel->Print(file,"pdf");

  TCanvas *cE = new TCanvas("EnergyCell", "Energy of Cell",50,50,700,500);
  ProcessHisto2D(CellE);
  gPad->SetLogz();
  CellE->Draw("COLZ");
  cE->Print(file,"pdf");

  TCanvas *c9 = new TCanvas("EMCClusE", "EMCAL cluster energy",50,50,700,500);
  ClusE->Sumw2();
  ProcessHisto(ClusE);
  ClusE->Rebin(10);
  gPad->SetLogy();
  ClusE->Draw();
  //c9->Print(file + ")","pdf");
  c9->Print(file,"pdf");

  TCanvas *ncls = new TCanvas("EMC_NCls", "Number of clusters in event",50,50,700,500);
  ProcessHisto(NCls,1.4,1,20);
  ProcessHisto(NClsE1,1.4,2,20);
  ProcessHisto(NClsE2,1.4,3,20);
  ProcessHisto(NClsE3,1.4,4,20);
  gPad->SetLogy();
  NClsE3->SetTitle("No. of EMCAL clusters in the event");
  NClsE3->GetXaxis()->SetRangeUser(0,100);
  NClsE3->Draw();
  NClsE1->Draw("same");
  NClsE2->Draw("same");
  NCls->Draw("same");
  TLegend *leg = new TLegend(0.33,0.72,0.58,0.89);
  ProcessLegend(leg);
  leg->AddEntry(NCls,"All","pl");
  leg->AddEntry(NClsE1,"ClsE > 0.1 GeV","pl");
  leg->AddEntry(NClsE2,"ClsE > 0.2 GeV","pl");
  leg->AddEntry(NClsE3,"ClsE > 0.5 GeV","pl");
  ncls->Print(file,"pdf");

  TCanvas *c11 = new TCanvas("TrackPt", "TrackPt",50,50,700,500);
  TrkPt->Sumw2();
  EMCTrkPt->Sumw2();
  EMCTrkMatchPt->Sumw2();
  TrkPt->Rebin(4);
  EMCTrkPt->Rebin(4);
  EMCTrkMatchPt->Rebin(4);
  TrkPt->GetYaxis()->SetRangeUser(1,2*TrkPt->GetBinContent(TrkPt->GetMaximumBin()));
  TrkPt->SetTitle("P_{T} ditribution");
  TrkPt->GetXaxis()->SetRangeUser(0,25);
  gPad->SetLogy();
  ProcessHisto(TrkPt);
  ProcessHisto(EMCTrkPt,1.4,2,20);
  ProcessHisto(EMCTrkMatchPt,1.4,4,20);
  TrkPt->Draw();
  EMCTrkPt->Draw("same");
  EMCTrkMatchPt->Draw("same");
  TLegend *leg = new TLegend(0.33,0.72,0.58,0.89);
  ProcessLegend(leg);
  leg->AddEntry(TrkPt,"All track","pl");
  leg->AddEntry(EMCTrkMatchPt,"EMC matched track","pl");
  leg->AddEntry(EMCTrkPt,"EMCAL selected track","pl");
  c11->Print(file,"pdf");

  TCanvas *c12 = new TCanvas("TrackEta", "Track eta distribution",50,50,700,500);
  Trketa->SetTitle("Track #eta ditribution");
  gPad->SetLogy();
  ProcessHisto(Trketa);
  ProcessHisto(EMCTrketa,1.4,2);
  Trketa->Draw();
  EMCTrketa->Draw("same");
  TLegend *leg = new TLegend(0.353,0.249,0.603,0.418);
  ProcessLegend(leg);
  leg->AddEntry(Trketa,"All track","l");
  leg->AddEntry(EMCTrketa,"EMCAL selected track","l");
  c12->Print(file,"pdf");

  TCanvas *c14 = new TCanvas("TrackPhi", "Track Phi distribution",50,50,700,500);
  Trkphi->SetTitle("Track #phi ditribution");
  Trkphi->GetYaxis()->SetRangeUser(1,10*Trkphi->GetBinContent(Trkphi->GetMaximumBin()));
  gPad->SetLogy();
  ProcessHisto(Trkphi);
  ProcessHisto(EMCTrkphi,1.4,2);
  Trkphi->Draw();
  EMCTrkphi->Draw("same");
  TLegend *leg = new TLegend(0.205,0.524,0.455,0.693);
  ProcessLegend(leg);
  leg->AddEntry(Trkphi,"All track","l");
  leg->AddEntry(EMCTrkphi,"EMCAL selected track","l");
  c14->Print(file,"pdf");

  TCanvas *c16 = new TCanvas("dedx","dEdx distribution",50,50,700,500);
  ProcessHisto2D(dEdx);
  dEdx->GetXaxis()->SetRangeUser(0,15);
  gPad->SetLogz();
  dEdx->Draw("colz");
  c16->Print(file,"pdf");

  TCanvas *c23 = new TCanvas("TPNnSig","TPC nsigma distribution",50,50,700,500);
  ProcessHisto2D(TPCNsig);
  TPCNsig->GetXaxis()->SetRangeUser(0,15);
  gPad->SetLogz();
  TPCNsig->Draw("colz");
  c23->Print(file,"pdf");

  TCanvas *nSigE0 = new TCanvas("TPNnSig_Eta_Pt0","TPC nsigma vs #eta for pT>2",50,50,700,500);
  gPad->SetLogz();
  TPCNsigEta0->SetTitle("TPC Nsigma vs #eta, p_{T} > 2 GeV/c");
  ProcessHisto2D(TPCNsigEta0);
  TPCNsigEta0->Draw("colz");
  nSigE0->Print(file,"pdf");
  /*
  TCanvas *nSigE1 = new TCanvas("TPNnSig_Eta_Pt1","TPC nsigma vs #eta for pT>3",50,50,700,500);
  gPad->SetLogz();
  TPCNsigEta1->SetTitle("TPC Nsigma vs #eta, p_{T} > 3 GeV/c");
  ProcessHisto2D(TPCNsigEta1);
  TPCNsigEta1->Draw("colz");
  nSigE1->Print(file,"pdf");

  TCanvas *nSigE2 = new TCanvas("TPNnSig_Eta_Pt2","TPC nsigma vs #eta for pT>5",50,50,700,500);
  gPad->SetLogz();
  TPCNsigEta2->SetTitle("TPC Nsigma vs #eta, p_{T} > 5 GeV/c");
  ProcessHisto2D(TPCNsigEta2);
  TPCNsigEta2->Draw("colz");
  nSigE2->Print(file,"pdf");
  */

  TCanvas *c18 = new TCanvas("EMCdeltaR","Distance of EMC cluster to its closest track",50,50,700,500);
  //EMCDeltaR->GetXaxis()->SetRangeUser(-0.08,0.08);
  //EMCDeltaR->GetYaxis()->SetRangeUser(-0.04,0.04);
  EMCDeltaR->GetXaxis()->SetTitle("#Delta#phi");
  EMCDeltaR->GetYaxis()->SetTitle("#Delta#eta");
  gPad->SetLogz();
  //ProcessHisto2D(EMCDeltaR);
  //EMCDeltaR->Draw("colz");
  c18->Divide(2,1);
  c18->cd(1);
  TH1D *hDeltaRX = EMCDeltaR->ProjectionX("hDeltaR");
  hDeltaRX->Draw();
  c18->cd(2);
  TH1D *hDeltaRY = EMCDeltaR->ProjectionY("hDeltaRY");
  hDeltaRY->Draw();
  c18->Print(file,"pdf");

  //TCanvas *c18_1 = new TCanvas("EMCdeltaR","Distance of EMC cluster to its closest track",50,50,700,500);
  //EMCDeltaR->Draw("colz");


  TCanvas *ClsPE = new TCanvas("ClusEtaPhi_AftMat", "EMCAL cluster eta and phi distribution aft trk matching",50,50,700,500);
  ProcessHisto2D(ClsEAftTrkMatch);
  ClsEAftTrkMatch->GetXaxis()->SetRangeUser(-0.7,0.7);
  gPad->SetLogz();
  ClsEAftTrkMatch->Draw("COLZ");
  ClsPE->Print(file,"pdf");
 

  /*
  TCanvas *c25 = new TCanvas("emcTrkTPCnsig","EMC matched track nsig distribution",50,50,700,500);
  EMCTrkTPCNsig->SetTitle("#phi,#eta distribution of clusters matched to tracks");
  EMCTrkTPCNsig->GetXaxis()->SetRangeUser(0,15);
  gPad->SetLogz();
  ProcessHisto2D(EMCTrkTPCNsig);
  EMCTrkTPCNsig->Draw("colz");
  c25->Print(file,"pdf");
  */

  TCanvas *c22 = new TCanvas("nsigEovP","TPCnsig vs E/p distribution",50,50,700,500);
  gPad->SetLogz();
  NsigEovP->SetTitle("TPCnsigma vs E/p, p_{T} > 3 GeV/c");
  NsigEovP->GetXaxis()->SetTitle("E/p");
  NsigEovP->GetYaxis()->SetTitle("#sigma_{TPC-dE/dx}");
  ProcessHisto2D(NsigEovP);
  NsigEovP->Draw("colz");

  TLine *LNsigEop[4];
  //LNsigEop[0] = new TLine(0.,-1,3.,-1);
  //LNsigEop[1] = new TLine(0.,3,3.,3);
  //LNsigEop[2] = new TLine(0.9,-10,0.9,10);
  //LNsigEop[3] = new TLine(1.3,-10,1.3,10);
  LNsigEop[0] = new TLine(0.9,-1,1.3,-1);
  LNsigEop[1] = new TLine(0.9, 3,1.3,3);
  LNsigEop[2] = new TLine(0.9,-1,0.9,3);
  LNsigEop[3] = new TLine(1.3,-1,1.3,3);
  for(int i=0; i<4; i++)LNsigEop[i]->SetLineStyle(2);
  for(int i=0; i<4; i++)LNsigEop[i]->Draw();

  //NsigEovP->ProjectionY()->Draw();
  //NsigEovP->ProjectionX("a",535,700);
  //a->Rebin(2);
  //a->Draw("e");
  c22->Print(file,"pdf");

  
  TCanvas *c22_1 = new TCanvas("EovP","E/p distribution",50,50,700,500);
  //ProcessHisto2D(EovP);
  TH1D* EovP = NsigEovP->ProjectionX("EovP",110,114);
  EovP->SetTitle("E/p with -1<n#sigma<3");
  EovP->GetXaxis()->SetRangeUser(0,15);
  gPad->SetLogz();
  EovP->Draw();
  TLine *LEop[2];
  LEop[0] = new TLine(0.85,0,0.85,EovP->GetBinContent(EovP->GetMaximumBin()));
  LEop[1] = new TLine(1.3,0,1.3,EovP->GetBinContent(EovP->GetMaximumBin()));
  for(int i=0; i<2; i++)LEop[i]->SetLineStyle(2);
  for(int i=0; i<2; i++)LEop[i]->Draw();
  //EovP->ProjectionY()->Draw();
  c22_1->Print(file,"pdf");

  TCanvas *c30 = new TCanvas("Invmass","Invariant mass",50,50,700,500);
  gPad->SetLogy();
  fULS->SetAxisRange(0,0.2);
  //fULS->SetMinimum(1);
  fULS->Draw();
  fLS->Draw("same");
  c30->Print(file,"pdf");


  TCanvas *c24 = new TCanvas("M20","M20 distribution",50,50,700,500);
  c24->Divide(2,1);
  c24->cd(1);
  M20->GetXaxis()->SetRangeUser(0,15);
  M20->SetTitle("");
  M20->GetXaxis()->SetTitle("p_{T} GeV/c");
  M20->GetYaxis()->SetTitle("M20");
  ProcessHisto2D(M20);
  gPad->SetLogz();
  M20->Draw("colz");
  c24->cd(2);
  M20->ProjectionY("fm20",20,100);
  fm20->GetXaxis()->SetRangeUser(0.02,1.0);
  fm20->Draw();
  c24->Print(file,"pdf");


  TCanvas *c25 = new TCanvas("M02","M02 distribution",50,50,700,500);
  c25->Divide(2,1);
  c25->cd(1);
  M02->GetXaxis()->SetRangeUser(0,15);
  M02->SetTitle("");
  M02->GetXaxis()->SetTitle("p_{T} GeV/c");
  M02->GetYaxis()->SetTitle("M02");
  ProcessHisto2D(M02);
  gPad->SetLogz();
  M02->Draw("colz");
  c25->cd(2);
  M02->ProjectionY("fm02",20,100);
  fm02->GetXaxis()->SetRangeUser(0.02,1);
  fm02->Draw();
  //c25->Print(file,"pdf");

  if(!iMC)
    {
     c25->Print(file + ")","pdf");
    }
  else
    {

     TCanvas *c99 = new TCanvas();
     gPad->SetLogy();
     TH1D *fPDG = fMomInfo->ProjectionX();
     fPDG->SetTitle("PDG;PDG;Counts");
     fPDG->Draw();
     c99->Print(file,"pdf");

     TCanvas *c99_0 = new TCanvas();
     c99_0->Divide(2,2);
     c99_0->cd(1);
     gPad->SetLogy();
     TH1D *fenPi = fMomInfo->ProjectionY("fenPi",108,114);
     fenPi->SetMinimum(1);
     fenPi->SetMaximum(1e+6);
     fenPi->SetTitle("enhanced #pi^{0};p_{T}(GeV/c);Counts");
     fenPi->Draw();
     c99_0->cd(2);
     gPad->SetLogy();
     TH1D *fenEta = fMomInfo->ProjectionY("fenEta",200,230);
     fenEta->SetMinimum(1);
     fenEta->SetMaximum(1e+6);
     fenEta->SetTitle("enhanced #eta;p_{T}(GeV/c);Counts");
     fenEta->Draw();
     c99_0->cd(3);
     gPad->SetLogy();
     TH1D *fenD = fMomInfo->ProjectionY("fenD",400,500);
     fenD->SetMinimum(1);
     fenD->SetTitle("enhanced D mesons;p_{T}(GeV/c);Counts");
     fenD->Draw();
     c99_0->cd(4);
     gPad->SetLogy();
     TH1D *fenB = fMomInfo->ProjectionY("fenB",500,600);
     fenB->SetMinimum(1);
     fenB->SetTitle("enhanced B mesons;p_{T}(GeV/c);Counts");
     fenB->Draw();

     //c99_0->Print(file + ")","pdf");
     c99_0->Print(file,"pdf");

     TCanvas *c31 = new TCanvas();
     c31->Divide(2,1);
     c31->cd(1);
     gPad->SetLogy();
     TH1D *fMassDal0 = fULSmc->ProjectionY("fMassDal0",2,2);
     fMassDal0->SetTitle("e+e- from #pi^{0} Dalitz;mass;Counts");
     fMassDal0->Draw();
     c31->cd(2);
     gPad->SetLogy();
     TH1D *fMassDal1 = fULSmc->ProjectionY("fMassDal1",3,3);
     fMassDal1->SetTitle("e+e- from #eta Dalitz;mass;Counts");
     fMassDal1->Draw();
     c31->Print(file + ")","pdf");

    }

  /*
  TCanvas *c26 = new TCanvas("M02EovP","M02 vs E/p distribution",50,50,700,500);
  M02EovP->SetTitle("M02 vs E/p, p_{T} > 1 GeV/c");
  M02EovP->GetXaxis()->SetTitle("E/p");
  M02EovP->GetYaxis()->SetTitle("M02");
  gPad->SetLogz();
  ProcessHisto2D(M02EovP);
  M02EovP->Draw("colz");

  TCanvas *c26 = new TCanvas("M20EovP","M20 vs E/p distribution",50,50,700,500);
  M20EovP->SetTitle("M20 vs E/p, p_{T} > 1 GeV/c");
  M20EovP->GetXaxis()->SetTitle("E/p");
  M20EovP->GetYaxis()->SetTitle("M20");
  gPad->SetLogz();
  ProcessHisto2D(M20EovP);
  M20EovP->Draw("colz");
  c26->Print(file + ")","pdf");
  */

  TFile *fout = new TFile("ElectronsMB.root","recreate");
  electron->Write("eArray");
  ClusEcent->Write("clE");
  cent->Write("Cent");
}
//----------------------------------------
void ProcessHisto(TH1 *h, Double_t size=1.4, Int_t col=1, Int_t style=20)
{
  gPad->SetTickx();
  gPad->SetTicky();
  h->SetMarkerSize(size);
  h->SetMarkerColor(col);
  h->SetLineColor(col);
  h->SetMarkerStyle(style);
  h->GetYaxis()->SetTitleOffset(1.0);
  h->GetYaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetLabelFont(42);
  h->GetXaxis()->SetLabelFont(42);
  h->GetYaxis()->SetTitleFont(42);
  h->GetXaxis()->SetTitleFont(42);
  h->GetXaxis()->SetTitleOffset(1.0);
  h->GetXaxis()->SetTitleSize(0.045);
  h->GetXaxis()->SetLabelSize(0.045);
}
void ProcessHisto2D(TH2 *h)
{
  //  h->SetLogz();
  h->GetYaxis()->SetTitleOffset(0.9);
  h->GetYaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetLabelSize(0.045);
  h->GetYaxis()->SetLabelFont(42);
  h->GetXaxis()->SetLabelFont(42);
  h->GetYaxis()->SetTitleFont(42);
  h->GetXaxis()->SetTitleFont(42);
  h->GetZaxis()->SetLabelFont(42);
  h->GetZaxis()->SetLabelSize(0.04);
  h->GetXaxis()->SetTitleOffset(0.9);
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetXaxis()->SetLabelSize(0.045);
}
void ProcessLegend(TLegend *leg)
{
  leg->SetFillColor(kWhite);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.04);
  leg->Draw();
}

void SetTrigQA(char *name0, char *name1, Bool_t &iMC)
{

  int iTender;
  cout << "Not Used Tender : 0" << endl;
  cout << "Used Tender : 1" << endl;
  cin >> iTender;  


  int iEMC;
  cout << "select EMCal or DCal" << endl;
  cout << "EMCal:0" << endl;
  cout << "DCal:1" << endl;
  cin >> iEMC;

  int itrig;
  cout << "select event selection" << endl;
  cout << "INT7 0" << endl;  //change for different triggers
  cout << "INT8 1" << endl;   //HFEemcQAINT8_woTender
  cout << "EMC7 2" << endl;   //HFEemcQAEMC7_woTender
  cout << "EMC8 3" << endl;   //HFEemcQAEMC8_woTender
  cout << "TrigGA 4" << endl;   //HFEemcQATrigGA_woTender
  cout << "TrigJE 5" << endl;   //HFEemcQATrigJE_woTender
  cout << "MC 6" << endl;
  cin >> itrig;
  if(itrig==6)iMC = kTRUE;


  if(iEMC==0)
    {
    if(itrig==0 || itrig==6)sprintf(name0,"PWGHF_hfeHFEemcQAINT7_EMC"); //change for different triggers
    if(itrig==1)sprintf(name0,"PWGHF_hfeHFEemcQAINT8_EMC");  //HFEemcQAINT8_woTender
    if(itrig==2)sprintf(name0,"PWGHF_hfeHFEemcQAEMC7_EMC");  //HFEemcQAEMC7_woTender
    if(itrig==3)sprintf(name0,"PWGHF_hfeHFEemcQAEMC8_EMC");  //HFEemcQAEMC8_woTender
    if(itrig==4)sprintf(name0,"PWGHF_hfeHFEemcQATrigGAEG1_EMC");  //HFEemcQATrigGA_woTender
    if(itrig==5)sprintf(name0,"PWGHF_hfeHFEemcQATrigJE_EMC");  //HFEemcQATrigJE_woTender
  
    if(iTender==0)
      {
       if(itrig==0 || itrig==6)sprintf(name1,"HFEemcQAINT7_woTender_EMC"); //change for different triggers
       if(itrig==1)sprintf(name1,"HFEemcQAINT8_woTender_EMC");
       if(itrig==2)sprintf(name1,"HFEemcQAEMC7_woTender_EMC");
       if(itrig==3)sprintf(name1,"HFEemcQAEMC8_woTender_EMC");
       if(itrig==4)sprintf(name1,"HFEemcQATrigGAEG1_woTender_EMC");
       if(itrig==5)sprintf(name1,"HFEemcQATrigJE_woTender_EMC");
     }
    else
     {
       if(itrig==0 || itrig==6)sprintf(name1,"HFEemcQAINT7_wTender_EMC"); //change for different triggers
       if(itrig==1)sprintf(name1,"HFEemcQAINT8_wTender_EMC");
       if(itrig==2)sprintf(name1,"HFEemcQAEMC7_wTender_EMC");
       if(itrig==3)sprintf(name1,"HFEemcQAEMC8_wTender_EMC");
       if(itrig==4)sprintf(name1,"HFEemcQATrigGAEG1_wTender_EMC");
       if(itrig==5)sprintf(name1,"HFEemcQATrigJE_wTender_EMC");
     }
   }
 else
   {
    if(itrig==0 || itrig==6)sprintf(name0,"PWGHF_hfeHFEemcQAINT7_DCAL"); //change for different triggers
    if(itrig==1)sprintf(name0,"PWGHF_hfeHFEemcQAINT8_DCAL");  //HFEemcQAINT8_woTender
    if(itrig==2)sprintf(name0,"PWGHF_hfeHFEemcQAEMC7_DCAL");  //HFEemcQAEMC7_woTender
    if(itrig==3)sprintf(name0,"PWGHF_hfeHFEemcQAEMC8_DCAL");  //HFEemcQAEMC8_woTender
    if(itrig==4)sprintf(name0,"PWGHF_hfeHFEemcQATrigGAEG1_DCAL");  //HFEemcQATrigGA_woTender
    if(itrig==5)sprintf(name0,"PWGHF_hfeHFEemcQATrigJE_DCAL");  //HFEemcQATrigJE_woTender

    if(iTender==0)
      {
       if(itrig==0)sprintf(name1,"HFEemcQAINT7_woTender_DCAL"); //change for different triggers
       if(itrig==1)sprintf(name1,"HFEemcQAINT8_woTender_DCAL");
       if(itrig==2)sprintf(name1,"HFEemcQAEMC7_woTender_DCAL");
       if(itrig==3)sprintf(name1,"HFEemcQAEMC8_woTender_DCAL");
       if(itrig==4)sprintf(name1,"HFEemcQATrigGAEG1_woTender_DCAL");
       if(itrig==5)sprintf(name1,"HFEemcQATrigJE_woTender_DCAL");
      }
      else
      {
       if(itrig==0)sprintf(name1,"HFEemcQAINT7_wTender_DCAL"); //change for different triggers
       if(itrig==1)sprintf(name1,"HFEemcQAINT8_wTender_DCAL");
       if(itrig==2)sprintf(name1,"HFEemcQAEMC7_wTender_DCAL");
       if(itrig==3)sprintf(name1,"HFEemcQAEMC8_wTender_DCAL");
       if(itrig==4)sprintf(name1,"HFEemcQATrigGAEG1_wTender_DCAL");
       if(itrig==5)sprintf(name1,"HFEemcQATrigJE_wTender_DCAL");
      }
   }
}
