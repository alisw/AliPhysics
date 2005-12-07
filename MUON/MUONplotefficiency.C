/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

// Macro to compare/plot the Ntuple stored in MUONefficiency.root
// comparison is done between the generated and reconstructed upsilons
// allowing to determine several important quantities 
// reconstruction efficiency (versus pt,y), invariant mass peak variations (vs pt,t)
// reconstructed tracks and trigger tracks  matching efficiency

// Christophe Suire, IPN Orsay

void MUONplotefficiency(Int_t fittype = 1, Int_t select = 0){
  
  //gStyle->Reset();
 
  Int_t SAVE=1;
  Int_t WRITE=1;

  Int_t fitfunc = fittype;
  // 0 gaussian fit +/- 450 MeV/c2
  // 1 gaussian fit +/- 250 MeV/c2
  // 2 approximated landau fit reverted 
  // 3 landau fit reverted convoluated with a gaussian 
  // 4 reverted landau fitlanUpsilon

  Float_t gaussianFitWidth = 0.450 ; // in GeV/c2
  if (fitfunc==1) gaussianFitWidth = 0.250 ; 

  // fit range : very important for LandauGauss Fit
  Float_t FitLow = 8.5 ;
  Float_t FitHigh = 10.2 ;
  
  // two ways to set the stat box 
  // for all histos
  gStyle->SetOptStat(1110);
  
  // place the stat box
  gStyle->SetStatX(0.5);
  gStyle->SetStatY(0.550);
  //gStyle->SetStatW(0.2);
  gStyle->SetStatH(0.2);

  
  // choose histos with stat
  // gStyle->SetOptStat(1111);
  //hist->SetStats(kTRUE); or hist->SetStats(kFALSE);
  
  //gStyle->SetOptFit(1111);
  gStyle->SetOptFit(1);
  //gStyle->SetOptFit(0);

  gStyle->SetOptTitle(0);
  gStyle->SetPalette(1,0);
  //gStyle->SetOptLogy(1);
  //gStyle->SetOptLogx(0);
  
  gStyle->SetFrameFillColor(10);
  gStyle->SetFillColor(10);
  gStyle->SetCanvasColor(10);


  TFile *effFile = new TFile("MUONefficiency.root","READ");

  TNtuple *Ktuple = (TNtuple*)effFile->Get("Ktuple");
  TNtuple *ESDtuple = (TNtuple*)effFile->Get("ESDtuple");
  
  
  Double_t MUON_MASS = 0.105658369;
  Double_t UPSILON_MASS = 9.4603 ;
  Double_t PI = 3.14159265358979312; 
  
  /*********************************/
  // Histograms limits and binning
  Float_t ptmin = 0.0;      Float_t ptmax =  20.0;     Int_t ptbins = 10; 
  Float_t ymin = -4.5;       Float_t ymax =  -2.;       Int_t ybins = 10;
  Float_t thetamin = 165.;  Float_t thetamax = 180.;   Int_t thetabins = 100;  
  
  Float_t etacutmin = -4.04813 ;
  Float_t etacutmax = -2.54209 ;
  
  /*********************************/
  // Values used to define the acceptance
  //
  Float_t thetacutmin = 171.0;  
  Float_t thetacutmax = 178.;
  
  Float_t ptcutmin = 0.;
  Float_t ptcutmax = 20.;

  Float_t masssigma  = 0.1; // 100 MeV/c2 is a correct estimation 
  Float_t masscutmin = 0; Float_t masscutmax = 0 ;    
  if (masssigma){
    masscutmin = UPSILON_MASS - 3.0*masssigma ; 
    masscutmax = UPSILON_MASS + 3.0*masssigma ;
  }
  else {
    masscutmin = UPSILON_MASS - 1.0 ; 
    masscutmax = UPSILON_MASS + 1.0 ;
  }


  Int_t REALISTIC_BACKGROUND = 0 ; 


  // here no cut on theta is necesary since during the simulation 
  // gener->SetChildThetaRange(171.0,178.0);
  // thus the upsilon are generated in 4 PI but only the ones for which the decay muons 
  // are in the dimuon arm acceptance
  // the pt cut is also "critical", in the generation part, upsilon are generated within  the range 
  // (ptcutmin,ptcutmax). During the reconstruction, the upsilon could have a pt greater than ptcutmax !!
  // Should these upsilons be cut or not.
  // probably not but they will be for now (~ 1/2000)

  // Another acceptance cut to add is a range for the  recontructed  invariant mass, it is 
  // obvious that an upsilon reconstructed with a mass  of 5 GeV/c2 is not correct. Thus
  Char_t UpsilonAccCutMC[200];  
  sprintf(UpsilonAccCutMC,"pt>= %.2f && pt<= %.2f",ptcutmin,ptcutmax);

  Char_t UpsilonAccCutESD[200];
  sprintf(UpsilonAccCutESD,"pt >=  %.2f && pt<= %.2f && minv>=  %.2f && minv<= %.2f",ptcutmin,ptcutmax,masscutmin,masscutmax);

  /*********************************/
  // Cut conditions (Id,Background...)

  Char_t IdcutUpsilonMC[100];            sprintf(IdcutUpsilonMC,"id==553");
  Char_t IdcutMuminusMC[100];            sprintf(IdcutMuminusMC,"id==13 && idmo==553");
  Char_t IdcutMuplusMC[100];             sprintf(IdcutMuplusMC,"id==-13 && idmo==553");  
  //means no cuts since we don't have the trackid propagated yet pt>0
  // now we have it 10/05/2005 but it's still not being used
  
  
  // Background is calculated in the MUONmass_ESD.C macro
  // Background calculation is meaningful only when Upsilon have been generated with realistic background
  // when it's a pure Upsilon file, the background lies artificially at the Upsilon mass 
  
  //no realistic background
  //Char_t BckgdCutUpsilonESD[100];           sprintf(BckgdCutUpsilonESD,"pt>0 && minv>7 && pt1>1 && pt2>1"); 
  //Char_t BckgdCutUpsilonESD[100];           sprintf(BckgdCutUpsilonESD,"pt>0 && minv>0 && pt1>0 && pt2>0"); 
  Char_t BckgdCutUpsilonESD[100];           sprintf(BckgdCutUpsilonESD,"pt>0 && minv>%f && minv<%f && pt1>0 && pt2>0",masscutmin,masscutmax); 

  // realistic background 
  // same cut + substract the background from hInvMassBg
  
  //sprintf(IdcutMuplusRec,"id==-13 && ch==0 && Ptrec!=0 && Yrec!=0");



  Char_t txt[50];
  TLatex *tex;
  TLegend *leg;

  /*******************************************************************************************************************/

  
  /*******************************/
  /*      Monte Carlo Part       */
  /*******************************/

  //----------------------------
  // Pt-rapidity distributions from Kinematics
  //----------------------------
  TH2F *hptyUpsilonMC = new TH2F("hptyUpsilonMC", " Monte Carlo #Upsilon",ptbins,ptmin,ptmax,ybins,ymin,ymax);
  hptyUpsilonMC->SetLineColor(1);
  hptyUpsilonMC->SetLineStyle(1);
  hptyUpsilonMC->SetLineWidth(2);
  Ktuple->Project("hptyUpsilonMC","y:pt",IdcutUpsilonMC);
  hptyUpsilonMC->GetXaxis()->SetTitle("P_{#perp}   [GeV/c]");
  hptyUpsilonMC->GetYaxis()->SetTitle("Rapidity");
 
  cout << "  " << endl;
  cout << "********************" << endl;
  cout << " Monte Carlo Tracks  "<< endl;
  cout << " " << hptyUpsilonMC->GetEntries() << " Upsilon in simulation " << endl;
  

  //******** Add acceptance cuts  - Theta and Pt
  TH2F *hptyUpsilonMCAcc = new TH2F("hptyUpsilonMCAcc", " Monte Carlo #Upsilon in acceptance",ptbins,ptmin,ptmax,ybins,ymin,ymax);
  hptyUpsilonMCAcc->SetLineColor(1);
  hptyUpsilonMCAcc->SetLineStyle(1);
  hptyUpsilonMCAcc->SetLineWidth(2);

  TString m1MC(IdcutUpsilonMC);
  TString m2MC(UpsilonAccCutMC);
  TString m3MC = m1MC + " && " + m2MC ;

  Ktuple->Project("hptyUpsilonMCAcc","y:pt",m3MC.Data());
  hptyUpsilonMCAcc->GetYaxis()->SetTitle("Rapidity");
  hptyUpsilonMCAcc->GetXaxis()->SetTitle("P_{#perp}   [GeV/c]");
  hptyUpsilonMCAcc->SetStats(1);
  hptyUpsilonMCAcc->Sumw2();
  
  TCanvas *c1 = new TCanvas("c1", "#Upsilon  Monte Carlo:  Pt vs Y",30,30,700,500);
  c1->Range(-1.69394,-0.648855,15.3928,2.77143);
  c1->SetBorderSize(2);
  c1->SetRightMargin(0.0229885);
  c1->SetTopMargin(0.0275424);
  c1->SetFrameFillColor(0);
  c1->cd();

  hptyUpsilonMCAcc->SetStats(0);
  hptyUpsilonMCAcc->Draw("LEGO2ZFB"); 
  //TLatex *tex = new TLatex(2,6,"#mu^{-}");
  //tex->SetTextSize(0.06);
  //tex->SetLineWidth(2);  
  //tex->Draw();
 
  sprintf(txt,"#Upsilon : %d entries",hptyUpsilonMCAcc->GetEntries());
  tex = new TLatex(-0.854829,0.794436,txt);
  tex->SetLineWidth(2);
  tex->Draw();

  c1->Modified();
  c1->Update();
  if (SAVE){
    c1->SaveAs("ptyUpsilonMCAcc.gif");
    c1->SaveAs("ptyUpsilonMCAcc.eps");
  }


  TH1F *hptUpsilonMCAcc = new TH1F("hptUpsilonMCAcc", " Monte Carlo #Upsilon in acceptance",ptbins,ptmin,ptmax);
  hptUpsilonMCAcc->SetLineColor(1);
  hptUpsilonMCAcc->SetLineStyle(1);
  hptUpsilonMCAcc->SetLineWidth(2);
  Ktuple->Project("hptUpsilonMCAcc","pt",m3MC.Data());
  hptUpsilonMCAcc->GetXaxis()->SetTitle("P_{#perp}   [GeV/c]");
  hptUpsilonMCAcc->SetStats(1);
  hptUpsilonMCAcc->Sumw2();
  
  TH1F *hyUpsilonMCAcc = new TH1F("hyUpsilonMCAcc", " Monte Carlo #Upsilon in acceptance",ybins,ymin,ymax);
  hyUpsilonMCAcc->SetLineColor(1);
  hyUpsilonMCAcc->SetLineStyle(1);
  hyUpsilonMCAcc->SetLineWidth(2);
  Ktuple->Project("hyUpsilonMCAcc","y",m3MC.Data());
  hyUpsilonMCAcc->GetXaxis()->SetTitle("Rapidity");
  hyUpsilonMCAcc->SetStats(1);
  hyUpsilonMCAcc->Sumw2();
  
  
  //   TH2F *hKAccptyCutsMuplus = new TH2F("hKAccptyCutsMuplus", "Monte Carlo #mu^{+} ",ptbins,ptmin,ptmax,ybins,ymin,ymax);
  //   hKAccptyCutsMuplus->SetLineColor(1);
  //   hKAccptyCutsMuplus->SetLineStyle(1);
  //   hKAccptyCutsMuplus->SetLineWidth(2);
  
  //   TString p1MC(IdcutMuplusMC);
  //   TString p2MC(AccCutMC);
  //   TString p3MC = p1MC + " && " + p2MC ;
  
  //   Ktuple->Project("hKAccptyCutsMuplus","y:pt",p3MC.Data());
  //   hKAccptyCutsMuplus->GetYaxis()->SetTitle("Rapidity");
  //   hKAccptyCutsMuplus->GetXaxis()->SetTitle("P_{#perp}   [GeV/c]");
  //   hKAccptyCutsMuplus->SetStats(1);
  //   hKAccptyCutsMuplus->Sumw2();
  
  
  cout << "  " << endl;
  cout << "********************" << endl;
  cout << " Monte Carlo Tracks  "<< endl;
  cout << " " << hptyUpsilonMCAcc->GetEntries() << " Upsilon in acceptance cuts " << endl;
  

  /*******************************************************************************************************************/

  /*******************************/
  /*  Reconstructed Tracks Study */
  /*******************************/

  //----------------------------
  // Pt-rapidity distributions from ESD : reconstructed tracks/particle
  //----------------------------
  TH2F *hptyUpsilonESD = new TH2F("hptyUpsilonESD", " Reconstucted #Upsilon",ptbins,ptmin,ptmax,ybins,ymin,ymax);
  hptyUpsilonESD->SetLineColor(1);
  hptyUpsilonESD->SetLineStyle(1);
  hptyUpsilonESD->SetLineWidth(2);
  ESDtuple->Project("hptyUpsilonESD","y:pt",BckgdCutUpsilonESD);
  hptyUpsilonESD->GetXaxis()->SetTitle("P_{#perp}   [GeV/c]");
  hptyUpsilonESD->GetYaxis()->SetTitle("Rapidity");

  cout << "  " << endl;
  cout << "********************" << endl;
  cout << " Reconstructed Tracks" << endl ;
  cout << " " << hptyUpsilonESD->GetEntries() << " Upsilon reconstructed " << endl;

  if (REALISTIC_BACKGROUND){
    TH2F *hptyUpsilonESDBck = new TH2F("hptyUpsilonESDBck", "Background",ptbins,ptmin,ptmax,ybins,ymin,ymax);
    hptyUpsilonESDBck->SetLineColor(1);
    hptyUpsilonESDBck->SetLineStyle(1);
    hptyUpsilonESDBck->SetLineWidth(2);
    ESDtupleBck->Project("hptyUpsilonESDBck","y:pt",BckgdCutUpsilonESD);
    hptyUpsilonESDBck->GetXaxis()->SetTitle("P_{#perp}   [GeV/c]");
    hptyUpsilonESDBck->GetYaxis()->SetTitle("Rapidity");
    cout << " with " << hptyUpsilonESDBck->GetEntries() << " Upsilons from Background (random mixing) " << endl;
  }

  // if something is wrong 
  if ( hptyUpsilonESD->GetEntries()==0) {cout << " No entries in hptyUpsilonESD " << endl ; break ;}

  //with Acc cuts - Theta and Pt
  TH2F *hptyUpsilonESDAcc = new TH2F("hptyUpsilonESDAcc", "Reconstructed #Upsilon",ptbins,ptmin,ptmax,ybins,ymin,ymax);
  hptyUpsilonESDAcc->SetLineColor(1);
  hptyUpsilonESDAcc->SetLineStyle(1);
  hptyUpsilonESDAcc->SetLineWidth(2);
  
  TString m1Rec(BckgdCutUpsilonESD);
  TString m2Rec(UpsilonAccCutESD);
  TString m3Rec = m1Rec + " && " + m2Rec ;

  ESDtuple->Project("hptyUpsilonESDAcc","y:pt",m3Rec.Data());
  hptyUpsilonESDAcc->GetYaxis()->SetTitle("Rapidity");
  hptyUpsilonESDAcc->GetXaxis()->SetTitle("P_{#perp}   [GeV/c]");
  hptyUpsilonESDAcc->SetStats(1);
  hptyUpsilonESDAcc->Sumw2();

  cout << "  " << endl;
  cout << "********************" << endl;
  cout << " Reconstructed Tracks" << endl ;
  cout << " " << hptyUpsilonESDAcc->GetEntries() << " Upsilon in acceptance cuts " << endl;

  if (REALISTIC_BACKGROUND){
    //with Acc cuts - Theta and Pt
    TH2F *hptyUpsilonESDBckAcc = new TH2F("hptyUpsilonESDBckAcc", "Reconstructed #Upsilon",ptbins,ptmin,ptmax,ybins,ymin,ymax);
    hptyUpsilonESDBckAcc->SetLineColor(1);
    hptyUpsilonESDBckAcc->SetLineStyle(1);
    hptyUpsilonESDBckAcc->SetLineWidth(2);
    
    TString m1Rec(BckgdCutUpsilonESD);
    TString m2Rec(UpsilonAccCutESD);
    TString m3Rec = m1Rec + " && " + m2Rec ;
    
    ESDtupleBck->Project("hptyUpsilonESDBckAcc","y:pt",m3Rec.Data());
    hptyUpsilonESDBckAcc->GetYaxis()->SetTitle("Rapidity");
    hptyUpsilonESDBckAcc->GetXaxis()->SetTitle("P_{#perp}   [GeV/c]");
    hptyUpsilonESDBckAcc->SetStats(1);
    hptyUpsilonESDBckAcc->Sumw2();
    cout << " with " << hptyUpsilonESDBckAcc->GetEntries() << " Upsilons from Background (random mixing) " << endl;
  }


  TCanvas *c100 = new TCanvas("c100", "#Upsilon Reconstructed in Acceptance: Pt vs Y",210,30,700,500);
  c100->Range(-1.69394,-0.648855,15.3928,2.77143);
  c100->SetBorderSize(2);
  c100->SetRightMargin(0.0229885);
  c100->SetTopMargin(0.0275424);
  c100->SetFrameFillColor(0);
  
  c100->cd();
  hptyUpsilonESDAcc->SetStats(0);
  hptyUpsilonESDAcc->Draw("LEGO2ZFB"); 
  sprintf(txt,"#Upsilon : %d entries",hptyUpsilonESDAcc->GetEntries());
  tex = new TLatex(-0.854829,0.794436,txt);
  tex->SetLineWidth(2);
  tex->Draw();

  c100->Update();
  if (SAVE){
    c100->SaveAs("ptyUpsilonESDAcc.gif");
    c100->SaveAs("ptyUpsilonESDAcc.eps");
  }


  if (REALISTIC_BACKGROUND){
    TCanvas *c110 = new TCanvas("c110", "#Upsilon Background Reconstructed in Acceptance: Pt vs Y",215,35,700,500);
    c110->Range(-1.69394,-0.648855,15.3928,2.77143);
    c110->SetBorderSize(2);
    c110->SetRightMargin(0.0229885);
    c110->SetTopMargin(0.0275424);
    c110->SetFrameFillColor(0);
    
    c110->cd();
    hptyUpsilonESDBckAcc->SetStats(0);
    hptyUpsilonESDBckAcc->Draw("LEGO2ZFB"); 
    sprintf(txt,"#Upsilon backround : %d entries",hptyUpsilonESDBckAcc->GetEntries());
    tex = new TLatex(-0.854829,0.794436,txt);
    tex->SetLineWidth(2);
    tex->Draw();
    
    
    c110->Update();
    if (SAVE){
      c110->SaveAs("ptyUpsilonESDBckAcc.gif");
      c110->SaveAs("ptyUpsilonESDBckAcc.eps");
    }
  }
  
  TH1F *hptUpsilonESDAcc = new TH1F("hptUpsilonESDAcc", " Monte Carlo #Upsilon in acceptance",ptbins,ptmin,ptmax);
  hptUpsilonESDAcc->SetLineColor(1);
  hptUpsilonESDAcc->SetLineStyle(1);
  hptUpsilonESDAcc->SetLineWidth(2);
  ESDtuple->Project("hptUpsilonESDAcc","pt",m3Rec.Data());
  hptUpsilonESDAcc->GetXaxis()->SetTitle("P_{#perp}   [GeV/c]");
  hptUpsilonESDAcc->SetStats(1);
  hptUpsilonESDAcc->Sumw2();
  
  TH1F *hyUpsilonESDAcc = new TH1F("hyUpsilonESDAcc", " Monte Carlo #Upsilon in acceptance",ybins,ymin,ymax);
  hyUpsilonESDAcc->SetLineColor(1);
  hyUpsilonESDAcc->SetLineStyle(1);
  hyUpsilonESDAcc->SetLineWidth(2);
  ESDtuple->Project("hyUpsilonESDAcc","y",m3Rec.Data());
  hyUpsilonESDAcc->GetXaxis()->SetTitle("Rapidity");
  hyUpsilonESDAcc->SetStats(1);
  hyUpsilonESDAcc->Sumw2();
  

  /*******************************************************************************************************************/

  /*******************************/
  /* Efficiencies calculations   */
  /*******************************/
  

  cout << "  " << endl;
  cout << "***********************" << endl;
  cout << " Integrated efficiency" << endl ;
  
  if (REALISTIC_BACKGROUND)
    cout << " UpsilonESDAcc/UpsilonMCAcc = " << (hptyUpsilonESDAcc->GetEntries()-hptyUpsilonESDBckAcc->GetEntries())/hptyUpsilonMCAcc->GetEntries()  << endl;
  else 
    cout << " UpsilonESDAcc/UpsilonMCAcc = " << hptyUpsilonESDAcc->GetEntries()/hptyUpsilonMCAcc->GetEntries()  << endl;

  TH2F *hEffptyUpsilon = new TH2F("hEffptyUpsilon", " #Upsilon Efficiency",ptbins,ptmin,ptmax,ybins,ymin,ymax);

  if (REALISTIC_BACKGROUND){
    TH2F *hptyUpsilonESDAccBckSubstracted = new TH2F("hptyUpsilonESDAccBckSubstracted","hptyUpsilonESDAccBckSubstracted",ptbins,ptmin,ptmax,ybins,ymin,ymax);
    hptyUpsilonESDAccBckSubstracted->Add(hptyUpsilonESDAcc,hptyUpsilonESDBckAcc,1,-1);
    hEffptyUpsilon->Divide(hptyUpsilonESDAccBckSubstracted,hptyUpsilonMCAcc,1,1);
  }
  else 
  hEffptyUpsilon->Divide(hptyUpsilonESDAcc,hptyUpsilonMCAcc,1,1);

  TCanvas *c1000 = new TCanvas("c1000", "#Upsilon efficiency : Pt vs Y",390,30,700,500);
  c1000->Range(-1.69394,-0.648855,15.3928,2.77143);
  c1000->SetBorderSize(2);
  c1000->SetRightMargin(0.0229885);
  c1000->SetTopMargin(0.0275424);
  c1000->SetFrameFillColor(0);

  c1000->cd();
  hEffptyUpsilon->SetStats(0);
  hEffptyUpsilon->Draw("LEGO2fz");



  TH1F *hEffptUpsilon = new TH1F("hEffptUpsilon", "#Upsilon Efficiency vs pt",ptbins,ptmin,ptmax);
  hEffptUpsilon->Divide(hptUpsilonESDAcc,hptUpsilonMCAcc,1,1);
  hEffptUpsilon->SetLineWidth(2);
  hEffptUpsilon->SetMinimum(0);
  hEffptUpsilon->SetStats(1);


  TCanvas *c1100 = new TCanvas("c1100", "#Upsilon efficiency : Pt",410,50,700,500);
  c1100->Range(-1.69394,-0.648855,15.3928,2.77143);
  c1100->SetBorderSize(2);
  c1100->SetRightMargin(0.0229885);
  c1100->SetTopMargin(0.0275424);
  c1100->SetFrameFillColor(0);

  c1100->cd();
  hEffptUpsilon->SetStats(0);
  hEffptUpsilon->GetXaxis()->SetTitle("P_{#perp}  [GeV/c]");
  hEffptUpsilon->GetYaxis()->SetTitle("Efficiency ");
  hEffptUpsilon->Draw("E");
   
  
  
  TH1F *hEffyUpsilon = new TH1F("hEffyUpsilon", "#Upsilon Efficiency vs y",ybins,ymin,ymax);
  hEffyUpsilon->Divide(hyUpsilonESDAcc,hyUpsilonMCAcc,1,1);
  hEffyUpsilon->SetLineWidth(2);
  hEffyUpsilon->SetStats(1);
  
  TCanvas *c1200 = new TCanvas("c1200", "#Upsilon efficiency : Y",430,70,700,500);
  c1200->Range(-1.69394,-0.648855,15.3928,2.77143);
  c1200->SetBorderSize(2);
  c1200->SetRightMargin(0.0229885);
  c1200->SetTopMargin(0.0275424);
  c1200->SetFrameFillColor(0);

  c1200->cd();
  hEffyUpsilon->SetStats(0);
  hEffyUpsilon->GetXaxis()->SetTitle("Rapidity");
  hEffyUpsilon->GetYaxis()->SetTitle("Efficiency ");
  hEffyUpsilon->Draw("E");

  c1000->Update();
  c1100->Update();
  c1200->Update();
  if (SAVE){
    c1000->SaveAs("EffptyUpsilon.gif");
    c1000->SaveAs("EffptyUpsilon.eps");
    c1100->SaveAs("EffptUpsilon.gif");
    c1100->SaveAs("EffptUpsilon.eps");
    c1200->SaveAs("EffyUpsilon.gif");
    c1200->SaveAs("EffyUpsilon.eps");
  }

 /*******************************************************************************************************************/

  /*******************************/
  /* Trigger matching            */
  /*******************************/


  Float_t triggerChi2Min = 0.; 
  Float_t triggerChi2Max = 7.5;

  Float_t invMassMin = 7.0 ;
  Float_t invMassMax = 11.0 ;
  Int_t   invMassBins = 100 ;
  
  //TString m1Rec(BckgdCutUpsilonESD);
  //TString m2Rec(UpsilonAccCutESD);
  TString  m1TG = m1Rec + " && " + m2Rec ;
  //cout << m1TG.Data() << endl;
  
  TH1F *hInvMassNoTriggerCut= new TH1F("hInvMassNoTriggerCut","hInvMassNoTriggerCut",invMassBins,invMassMin,invMassMax);
  hInvMassNoTriggerCut->GetXaxis()->SetTitle("Mass [GeV/c^{2}]");
  hInvMassNoTriggerCut->GetYaxis()->SetTitle("Counts");
  ESDtuple->Project("hInvMassNoTriggerCut","minv",m1TG.Data());

  //Add trigger UnlikePairAllPt
  TString  m2TG = m1TG + " && (tw & 0x800) == 2048" ;
  TH1F *hInvMassUpsilonTrigger= new TH1F("hInvMassUpsilonTrigger","hInvMassUpsilonTrigger",invMassBins,invMassMin,invMassMax);
  hInvMassUpsilonTrigger->GetXaxis()->SetTitle("Mass [GeV/c^{2}]");
  hInvMassUpsilonTrigger->GetYaxis()->SetTitle("Counts");
  ESDtuple->Project("hInvMassUpsilonTrigger","minv",m2TG.Data());


  cout << "  " << endl;
  cout << "********************" << endl;
  cout << " TriggerMatching" << endl ;
  cout << " " << hInvMassUpsilonTrigger->GetEntries() << " Upsilon with trigger UnlikePairAllPt " << endl;

  TString  m2TGNo = m1TG + " && (tw & 0x800) != 2048" ;
  TH1F *hUpsilonTriggerNo= new TH1F("hUpsilonTriggerNo","hUpsilonTriggerNo",32769,0,32768);
  hUpsilonTriggerNo->GetXaxis()->SetTitle("Mass [GeV/c^{2}]");
  hUpsilonTriggerNo->GetYaxis()->SetTitle("Counts");
  ESDtuple->Project("hUpsilonTriggerNo","tw",m2TGNo.Data());

  


  //Add matching rec/trig for 2 tracks 
  TString m3TG = m2TG + " && trig1 > 0 && trig2 > 0" ;  
  TH1F *hInvMassUpsilonTriggerTwoMatch = new TH1F("hInvMassUpsilonTriggerTwoMatch","hInvMassUpsilonTrigger with 2 matching tracks",invMassBins,invMassMin,invMassMax);
  hInvMassUpsilonTriggerTwoMatch->GetXaxis()->SetTitle("Mass [GeV/c^{2}]");
  hInvMassUpsilonTriggerTwoMatch->GetYaxis()->SetTitle("Counts");
  ESDtuple->Project("hInvMassUpsilonTriggerTwoMatch","minv",m3TG.Data());


  //Add matching rec/trig for 1  tracks 
  TString m4TG = m2TG + " && (trig1 > 0 || trig2 > 0)" ;  
  TH1F *hInvMassUpsilonTriggerOneMatch= new TH1F("hInvMassUpsilonTriggerOneMatch","hInvMassUpsilonTrigger with 1 matching tracks",invMassBins,invMassMin,invMassMax);
  hInvMassUpsilonTriggerOneMatch->GetXaxis()->SetTitle("Mass [GeV/c^{2}]");
  hInvMassUpsilonTriggerOneMatch->GetYaxis()->SetTitle("Counts");
  ESDtuple->Project("hInvMassUpsilonTriggerOneMatch","minv",m4TG.Data());

  TString m5TG = m2TG + " && (trig1 == 0 && trig2 == 0)" ;  
  TH1F *hInvMassUpsilonTriggerNoMatch= new TH1F("hInvMassUpsilonTriggerNoMatch","hInvMassUpsilonTrigger with no matching tracks",invMassBins,invMassMin,invMassMax);
  hInvMassUpsilonTriggerNoMatch->GetXaxis()->SetTitle("Mass [GeV/c^{2}]");
  hInvMassUpsilonTriggerNoMatch->GetYaxis()->SetTitle("Counts");
  ESDtuple->Project("hInvMassUpsilonTriggerNoMatch","minv",m5TG.Data());

  TCanvas *c2 = new TCanvas("c2", "#Upsilon Trigger efficiency",30,90,700,500);
  c2->Range(-1.69394,-0.648855,15.3928,2.77143);
  c2->SetBorderSize(2);
  c2->SetRightMargin(0.0229885);
  c2->SetTopMargin(0.0275424);
  c2->SetFrameFillColor(0);
  c2->SetLogy(1);

  c2->cd();
  hInvMassNoTriggerCut->SetStats(0);
  hInvMassNoTriggerCut->SetLineWidth(2);
  hInvMassNoTriggerCut->SetLineStyle(3);
  hInvMassNoTriggerCut->SetLineColor(1);
  hInvMassNoTriggerCut->GetYaxis()->SetTitleOffset(1.2);
  hInvMassNoTriggerCut->Draw();
  hInvMassUpsilonTrigger->SetLineWidth(2);
  hInvMassUpsilonTrigger->SetLineStyle(0);
  hInvMassUpsilonTrigger->SetLineColor(2);
  hInvMassUpsilonTrigger->Draw("same");
  hInvMassUpsilonTriggerTwoMatch->SetLineWidth(2);
  hInvMassUpsilonTriggerTwoMatch->SetLineStyle(1);
  hInvMassUpsilonTriggerTwoMatch->SetLineColor(4);
  hInvMassUpsilonTriggerTwoMatch->Draw("same");
  hInvMassUpsilonTriggerOneMatch->SetLineWidth(2);
  hInvMassUpsilonTriggerOneMatch->SetLineStyle(2);
  hInvMassUpsilonTriggerOneMatch->SetLineColor(51);
  hInvMassUpsilonTriggerOneMatch->Draw("same");
  hInvMassUpsilonTriggerNoMatch->SetLineWidth(2);
  hInvMassUpsilonTriggerNoMatch->SetLineStyle(2);
  hInvMassUpsilonTriggerNoMatch->SetLineColor(46);
  hInvMassUpsilonTriggerNoMatch->Draw("same");

  TLegend *leg = new TLegend(0.12,0.6,0.50,0.89);
  leg->SetHeader("Reconstructed #Upsilon Invariant Mass");

  leg->AddEntry(hInvMassNoTriggerCut,Form("All  (%.0f cnts)",hInvMassNoTriggerCut->GetEntries()),"l");
  leg->AddEntry(hInvMassUpsilonTrigger,Form("UnlikePairAllPt Trigger (%.0f cnts)",hInvMassUpsilonTrigger->GetEntries()),"l");
  leg->AddEntry(hInvMassUpsilonTriggerTwoMatch,Form("UPAllPt Trig. and 2 tracks matches (%.0f cnts)",hInvMassUpsilonTriggerTwoMatch->GetEntries()),"l");  
  leg->AddEntry(hInvMassUpsilonTriggerOneMatch,Form("UPAllPt Trig. and 1 track match (%.0f cnts)",hInvMassUpsilonTriggerOneMatch->GetEntries()),"l");  
  leg->AddEntry(hInvMassUpsilonTriggerNoMatch,Form("UPAllPt Trig. and no matching track (%.0f cnts)",hInvMassUpsilonTriggerNoMatch->GetEntries()),"l");  
  leg->Draw();

  c2->Update();
  if(SAVE){
    c2->SaveAs("TriggerMatching.gif");
    c2->SaveAs("TriggerMatching.eps");
  }
  
 /*******************************************************************************************************************/

  /*******************************/
  /* Mass resolution             */
  /*******************************/
  
  const Int_t nofMassHistograms = ptbins ; 
  invMassMin = 7.0 ;
  invMassMax = 11.0 ;
  invMassBins = 100 ;
  
 //Float_t ptmin = 0.0;      
  //Float_t ptmax =  16.0;  Int_t ptbins = 8;    
  Float_t ptbinssize = ptmax/ptbins; 


  Float_t ptbinslimits[nofMassHistograms+1];
  Float_t ptbinscenters[nofMassHistograms];
  
  for(Int_t as = 0 ; as < nofMassHistograms+1 ; as++){
    ptbinslimits[as] = as*ptbinssize ;
    //cout << "ptbinslimits["<<as<<"] = " << ptbinslimits[as] << endl ;
  }
  for(Int_t as = 0 ; as < nofMassHistograms ; as++) {
    ptbinscenters[as] = ptbinslimits[as] + (ptbinslimits[as+1] - ptbinslimits[as])/2 ;
    //cout << "ptbinscenters[as] = " << ptbinscenters[as] << endl ;
  }
  

  TH1F **hInvMassInPtBins = new TH1F* [nofMassHistograms] ;
  TString histExt;
  Char_t  theExt[20];
  TString histCut;
  Char_t  theCut[100];
  
  
  // Load  the fit functions 
  //gROOT->ProcessLine(".L /home/suire/AliRoot/macros/FIT/FitFunction.C");
  
  TF1 *fitFunc ;  
  
  if (fitfunc==3){
    fitFunc = new TF1("fitlangaussUpsilon",fitlangaussUpsilon,FitLow,FitHigh,6);
    fitFunc->SetParNames("Constant","Mean value","Width","SigmaGauss","LowFitVal","HighFitVal");
    fitFunc->SetLineColor(2);
    fitFunc->SetLineWidth(2);
    fitFunc->SetLineStyle(2);
    
    Float_t  *tParams = new Float_t[6] ;
    tParams[0] = 5000;
    tParams[1] = 9.47;   
    tParams[2] = 0.05;   //0.5
    tParams[3] = 0.05;   // 1.
    
    tParams[4] = FitLow;
    tParams[5] = FitHigh;

    fitFunc->SetParameters(tParams[0],tParams[1],tParams[2],tParams[3],tParams[4],tParams[5]);
    fitFunc->FixParameter(4,FitLow);
    fitFunc->FixParameter(5,FitHigh);
    
    fitFunc->SetParLimits(1, 9.1 ,9.7);
    fitFunc->SetParLimits(2, 0.001 ,0.1);
    fitFunc->SetParLimits(3, 0.001 ,0.5);

    // for special case one may use 
    //fitFunc->SetParameter(1,9.35);

    TString FitFuncName = "fitlangaussUpsilon"; 
  }
  
  if (fitfunc==2){
    fitFunc = new TF1("fitlan_a",fitlan_a,FitLow,FitHigh,3);
    fitFunc->SetParNames("Constant","Mean value","Width");
    fitFunc->SetLineColor(2);
    fitFunc->SetLineWidth(2);
    fitFunc->SetLineStyle(2);
    
    Float_t  *tParams = new Float_t[3] ;
    tParams[0] = 1000;
    tParams[1] = 9.47;
    tParams[2] = 0.05;
    fitFunc->SetParameters(tParams[0],tParams[1],tParams[2]);
    
    TString FitFuncName = "fitlan_a";  
  }
  
  if (fitfunc==0 || fitfunc==1){
    //mass upsilon 9.4603 GeV/c2
    //fitFunc = new TF1("gaussian",gaussian,9,10,3);
    fitFunc = new TF1("gaussian","[0]*TMath::Exp(-0.5*(x-[1])*(x-[1])/[2]/[2])",9,10);
    //TF1 *fb2 = new TF1("fa3","TMath::Landau(x,[0],[1],0)",-5,10);
    //if (sigma == 0) return 1.e30;
    // Double_t arg = (x-mean)/sigma; (x-[1])*(x-[1])/[2]/[2]
    //Double_t res = TMath::Exp(-0.5*arg*arg);
    //if (!norm) return res;
    //return res/(2.50662827463100024*sigma); //sqrt(2*Pi)=2.50662827463100024

   fitFunc->SetParNames("Constant","Mean value","Width");
    fitFunc->SetLineColor(2);
    fitFunc->SetLineWidth(2);
    fitFunc->SetLineStyle(2);
    
    Float_t  *tParams = new Float_t[3] ;
    tParams[0] = 500;
    tParams[1] = 9.47;
    tParams[2] = 0.1;
    fitFunc->SetParameters(tParams[0],tParams[1],tParams[2]);
    
    TString FitFuncName = "gaussian";  
  }
  
  if (fitfunc==4){
    fitFunc = new TF1("fitlanUpsilon",fitlanUpsilon,FitLow,FitHigh,3);
    fitFunc->SetParNames("Constant","Mean value","Width");
    fitFunc->SetLineColor(2);
    fitFunc->SetLineWidth(2);
    fitFunc->SetLineStyle(2);
    
    Float_t  *tParams = new Float_t[3] ;
    tParams[0] = 1000;
    tParams[1] = 9.47;
    tParams[2] = 0.05;
    fitFunc->SetParameters(tParams[0],tParams[1],tParams[2]);
    
    TString FitFuncName = "fitlanUpsilon";  
  }


  Float_t fitResultsWidth[nofMassHistograms];
  Float_t fitResultsWidthErr[nofMassHistograms];
  
  Float_t fitResultsMean[nofMassHistograms];
  Float_t fitResultsMeanErr[nofMassHistograms];
  
  Float_t errx[nofMassHistograms];



  TCanvas *call = new TCanvas("call", "#Upsilon : invariant mass spectra",30,330,700,500);
  call->Range(-1.69394,-0.648855,15.3928,2.77143);
  call->SetBorderSize(2);
  //call->SetRightMargin(0.2229885);
  call->SetRightMargin(0.0229885);
  call->SetTopMargin(0.0275424);
  call->SetFrameFillColor(0);
  call->cd();
  TH1F* hInvMass = new TH1F("hInvMass","Inv. Mass Pt,Y integrated",30,0,12);
  ESDtuple->Project("hInvMass","minv");
  hInvMass->SetLineWidth(2);
  hInvMass->Draw("HE");
  if(REALISTIC_BACKGROUND){  
    TH1F* hInvMassBck = new TH1F("hInvMassBck","Background Pt,Y integrated",30,0,12);
    ESDtupleBck->Project("hInvMassBck","minv");
    hInvMassBck->SetLineWidth(2);
    hInvMassBck->SetLineStyle(2);
    hInvMassBck->SetLineColor(2);    
    hInvMassBck->Draw("same");
  }
  call->Modified();
  call->Update();


  TCanvas *cfit = new TCanvas("cfit", "#Upsilon : invariant mass fit",30,330,700,500);
  cfit->Range(-1.69394,-0.648855,15.3928,2.77143);
  cfit->SetBorderSize(2);
  //cfit->SetRightMargin(0.2229885);
  cfit->SetRightMargin(0.0229885);
  cfit->SetTopMargin(0.0275424);
  cfit->SetFrameFillColor(0);
  
  //TH1F *hInvMassAllClone = (TH1F*) gROOT->FindObject("hInvMassAll");
  
  TH1F* hInvMassAll = new TH1F("hInvMassAll","Inv. Mass Pt,Y integrated",invMassBins,invMassMin,invMassMax);
  sprintf(theCut,"pt> %.2f && pt<= %.2f",ptbinslimits[0],ptbinslimits[nofMassHistograms]);
  //cout << "theCut" << theCut << endl ;
  ESDtuple->Project("hInvMassAll","minv",theCut);
  hInvMassAll->Draw();
  
  cfit->Update();
  if (SAVE){
    cfit->SaveAs("UpsilonMass.gif");
    cfit->SaveAs("UpsilonMass.eps");
  }

  if(REALISTIC_BACKGROUND){  
    TH1F* hInvMassAllBck = new TH1F("hInvMassAllBck","Background Pt,Y integrated",invMassBins,invMassMin,invMassMax);
    ESDtupleBck->Project("hInvMassAllBck","minv",theCut);
    hInvMassAllBck->SetLineWidth(2);
    hInvMassAllBck->SetLineStyle(2);
    hInvMassAllBck->SetLineColor(2);    
    hInvMassAllBck->Draw("same");
  }
  cfit->Modified();
  cfit->Update();

#if 1
  //Fit also the pt-integrated mass histogram 
   cfit->cd();
   hInvMassAll->Fit(FitFuncName.Data(),"rv");
   if (FitFuncName == "gaussian")
     fitFunc->SetRange(fitFunc->GetParameter(1)-gaussianFitWidth,fitFunc->GetParameter(1)+gaussianFitWidth);
   hInvMassAll->Fit(FitFuncName.Data(),"rv");
   cfit->Modified();
   cfit->Update();
   
   cout << " " << endl;
   cout << "********************************************" << endl;
   cout << " Fit ("<<FitFuncName.Data()<<" results of the pt-integrated mass peak " << endl ;
   cout << "   Mean value : " << fitFunc->GetParameter(1) << " +/- " << fitFunc->GetParError(1) << endl ;
   cout << "   Width value : " << fitFunc->GetParameter(2) << " +/- " << fitFunc->GetParError(2) << endl ;
   cout << "   ChiSquare of the fit : " << fitFunc->GetChisquare()<< " / " << fitFunc->GetNDF() << endl ;
   cout << "********************************************" << endl;
   cout << " " << endl;
#endif
   
   TCanvas *cptfits = new TCanvas("cptfits", "#Upsilon : invariant mass fit",30,330,700,500);
   cptfits->Range(-1.69394,-0.648855,15.3928,2.77143);
   cptfits->SetBorderSize(2);
   //cptfits->SetRightMargin(0.2229885);
   cptfits->SetRightMargin(0.0229885);
   cptfits->SetTopMargin(0.0275424);
   cptfits->SetFrameFillColor(0);

   Float_t chi2sum = 0.;
  
  for (Int_t qw = 0 ; qw < nofMassHistograms ; qw++){
    sprintf(theExt,"_%d_%d",ptbinslimits[qw],ptbinslimits[qw+1]);
    histExt= theExt;
    hInvMassInPtBins[qw] = new TH1F("hInvMassInPtBins"+histExt,"hInvMassInPtBins"+histExt,invMassBins,invMassMin,invMassMax);
    hInvMassInPtBins[qw]->GetXaxis()->SetTitle("Mass [GeV/c^{2}]");
    hInvMassInPtBins[qw]->GetYaxis()->SetTitle("Counts");
     
    sprintf(theCut,"pt> %.2f && pt<= %.2f",ptbinslimits[qw],ptbinslimits[qw+1]);
    histCut = theCut; 
    ESDtuple->Project("hInvMassInPtBins"+histExt,"minv",histCut.Data());
    
    if (hInvMassInPtBins[qw]->GetEntries() > 80){
      cptfits->cd();

      if (FitFuncName == "gaussian"){
	//prefit with a gaussian to get the maximum set a correct range
	hInvMassInPtBins[qw]->Fit(FitFuncName.Data(),"rv");
	fitFunc->SetRange(fitFunc->GetParameter(1)-gaussianFitWidth,fitFunc->GetParameter(1)+gaussianFitWidth);
      }
      
      hInvMassInPtBins[qw]->SetStats(1);
      hInvMassInPtBins[qw]->Draw("HE");
      hInvMassInPtBins[qw]->Fit(FitFuncName.Data(),"rv");

      
      cout << "==============" << endl;
      cout << "ChiSquare of the fit : " << fitFunc->GetChisquare()<< " / " << fitFunc->GetNDF() 
	   << " = " << fitFunc->GetChisquare()/fitFunc->GetNDF() <<endl;
      cout << "==============" << endl;
      chi2sum += fitFunc->GetChisquare()/fitFunc->GetNDF();
 
      cptfits->Modified();
      cptfits->Update();
      
      fitResultsWidth[qw] = TMath::Abs(fitFunc->GetParameter(2));
      fitResultsWidthErr[qw] =  fitFunc->GetParError(2);

      fitResultsMean[qw] = TMath::Abs(fitFunc->GetParameter(1));
      fitResultsMeanErr[qw] =  fitFunc->GetParError(1);
      
      if(FitFuncName == "fitlangaussUpsilon"){	
	// width = gauss_width + width landau
	fitResultsWidth[qw] = TMath::Abs(fitFunc->GetParameter(2)) +  TMath::Abs(fitFunc->GetParameter(3)); 
	fitResultsWidthErr[qw] =  TMath::Abs(fitFunc->GetParError(2)) +  TMath::Abs(fitFunc->GetParError(3));
	
	// problem with the mean of the fit (parameter1)  which doesn't give the peak position 
	fitResultsMean[qw] = TMath::Abs(fitFunc->GetMaximumX());
	fitResultsMeanErr[qw] =  fitFunc->GetParError(1); 
	//this could work but it's not very nice....
      }

    }
    else {
      fitResultsWidth[qw] = 0.;
      fitResultsWidthErr[qw] = 0.;
      fitResultsMean[qw] = 0.;
      fitResultsMeanErr[qw] = 0.;
    }

    errx[qw] = 0.;
    
  }
  
  
  //From fit results extract a Mean value for the Mass and width
  cout << " " << endl;
  cout << "********************************************" << endl;
  cout << " Fit Function : " << FitFuncName  << " " ; 
  if (FitFuncName == "gaussian" ) cout << "+/- " << gaussianFitWidth << " MeV/c2" ;  
  cout <<  endl ;  
  cout << " Fit results :  mean values  " << endl ;
  Float_t meanMass = 0 ;
  Float_t meanMassError = 0 ;
  Float_t meanWidth = 0 ;
  Float_t meanWidthError = 0 ;
  Int_t cnt = 0 ;
  for (Int_t qw = 0 ; qw < nofMassHistograms ; qw++){
    if ( fitResultsMean[qw] > 9. &&   fitResultsMean[qw] <10 &&  fitResultsWidth[qw] > 0 &&  fitResultsWidth[qw] < 1){
      meanWidth += fitResultsWidth[qw] ;
      meanWidthError += fitResultsWidthErr[qw];
      meanMass += fitResultsMean[qw] ;
      meanMassError += fitResultsMeanErr[qw];
      cnt++;
    }
  }
  if (cnt==0) {
    cout << "Fitting procedure didn't work" << endl;  
    break ;
  }

  cout << " Mass : " << meanMass/cnt << " +/- " << meanMassError/cnt << endl ;
  cout << " Width : " << meanWidth/cnt << " +/- " << meanWidthError/cnt << endl ;
  cout << " Mean Chi2 of the fits : " << chi2sum/nofMassHistograms << endl ;
  cout << "********************************************" << endl;
  cout << " " << endl;
  

  cout << "  " << endl;
  cout << "***********************" << endl;
  cout << " Integrated efficiency (+/- "<< 3*masssigma << " GeV around Upsilon mass cut)" << endl ;
  cout << " UpsilonESDAcc/UpsilonMCAcc = " << hptyUpsilonESDAcc->GetEntries() << "/" << hptyUpsilonMCAcc->GetEntries()  << " = "  << hptyUpsilonESDAcc->GetEntries()/hptyUpsilonMCAcc->GetEntries()  <<endl;
  cout << "***********************" << endl;
  cout << "  " << endl;


  cout << "  " << endl;
  cout << "***********************" << endl;
  cout << " Trigger  efficiency" << endl ;
  cout << " Two muons matching = " << hInvMassUpsilonTriggerTwoMatch->GetEntries() << "/" << hInvMassUpsilonTrigger->GetEntries() << " = " <<  hInvMassUpsilonTriggerTwoMatch->GetEntries()/hInvMassUpsilonTrigger->GetEntries() << endl;
  cout << " Single muon matching  = " << hInvMassUpsilonTriggerOneMatch->GetEntries() << "/" << hInvMassUpsilonTrigger->GetEntries() << " = " <<  hInvMassUpsilonTriggerOneMatch->GetEntries()/hInvMassUpsilonTrigger->GetEntries() << endl;
  cout << " No matching  = " << hInvMassUpsilonTriggerNoMatch->GetEntries() << "/" << hInvMassUpsilonTrigger->GetEntries() << " = " <<  hInvMassUpsilonTriggerNoMatch->GetEntries()/hInvMassUpsilonTrigger->GetEntries() <<endl;
  cout << "***********************" << endl;
  cout << "  " << endl;

  
  

  TH2F *hFrame = new TH2F("hFrame", "A hFrame",ptbins,ptmin,ptmax,100,fitResultsWidth[1]-fitResultsWidth[1]*2,fitResultsWidth[1]+fitResultsWidth[1]*2);
  hptyUpsilonMCAcc->SetLineColor(1);
  hptyUpsilonMCAcc->SetLineStyle(1);
  hptyUpsilonMCAcc->SetLineWidth(2);
  
  
  TCanvas *cA = new TCanvas("cA", "#Upsilon : Width from fit",30,330,700,500);
  cA->Range(-1.69394,-0.648855,15.3928,2.77143);
  cA->SetBorderSize(2);
  cA->SetRightMargin(0.0229885);
  cA->SetTopMargin(0.0275424);
  cA->SetFrameFillColor(0);
  cA->cd();
  
  //TH2F *hFrameA = new TH2F("hFrameA", "A hFrameA",ptbins,ptmin,ptmax,100,fitResultsWidth[1]-fitResultsWidth[1]*2,fitResultsWidth[1]+fitResultsWidth[1]*2);
  TH2F *hFrameA = new TH2F("hFrameA", "A hFrameA",ptbins,ptmin,ptmax,100,0.,2*fitResultsWidth[1]);
  hFrameA->GetYaxis()->SetTitle("Peak Width [GeV/c^{2}]");
  hFrameA->GetYaxis()->SetTitleOffset(1.2);
  hFrameA->GetXaxis()->SetTitle("P_{#perp}   [GeV/c]");
  hFrameA->SetStats(0);
  hFrameA->Draw(); 
  
  
  TGraphErrors *grWidth = new TGraphErrors(nofMassHistograms,ptbinscenters,fitResultsWidth,errx,fitResultsWidthErr);
  grWidth->SetTitle("Upsilon Mass Width");
  grWidth->SetMarkerColor(4);
  grWidth->SetMarkerStyle(21);
  grWidth->Draw("P");
  
  
  
  TCanvas *cB = new TCanvas("cB", "#Upsilon : Mean from fit",50,350,700,500);
  cB->Range(-1.69394,-0.648855,15.3928,2.77143);
  cB->SetBorderSize(2);
  cB->SetRightMargin(0.0229885);
  cB->SetTopMargin(0.0275424);
  cB->SetFrameFillColor(0);
  cB->cd();
  
  //TH2F *hFrameB = new TH2F("hFrameB", "A hFrameB",ptbins,ptmin,ptmax,100,fitResultsMean[1]-fitResultsMean[1]*2,fitResultsMean[1]+fitResultsMean[1]*2);
  TH2F *hFrameB = new TH2F("hFrameB", "A hFrameB",ptbins,ptmin,ptmax,100,8.0,10.);
  hFrameB->GetYaxis()->SetTitle("Peak Position [GeV/c^{2}]");
  hFrameB->GetXaxis()->SetTitle("P_{#perp}   [GeV/c]");
  hFrameB->SetStats(0);
  hFrameB->Draw(); 
  
  
  
  
  //  TGaxis *Mean = new TGaxis(16,fitResultsWidth[1]-fitResultsWidth[1]*2,16,fitResultsWidth[1]+fitResultsWidth[1]*2,7,12,510,"+");
  //  Mean->SetLabelOffset(0.050);
  //   Mean->SetLabelSize(0.04);
  //   Mean->SetTickSize(0.03);
  //   Mean->SetGridLength(0);
  //   Mean->SetTitleOffset(0.1);
  //   Mean->SetTitleOffset(-1.1);
  //   Mean->SetTitleSize(0.05);
  //   Mean->SetTitle("Mean [Gev/c^2]");
  //   Mean->Draw();
  
  
  TGraphErrors *grMean = new TGraphErrors(nofMassHistograms,ptbinscenters,fitResultsMean,errx,fitResultsMeanErr);
  grMean->SetTitle("Upsilon Mass Mean");
  grMean->SetMarkerColor(4);
  grMean->SetMarkerStyle(21);
  grMean->Draw("P");
  
  cA->Update();
  cB->Update();
  if (SAVE){
    cB->SaveAs("UpsilonMeanVsPt.gif");
    cB->SaveAs("UpsilonMeanVsPt.eps");
    cA->SaveAs("UpsilonWidthVsPt.gif");
    cA->SaveAs("UpsilonWidthVsPt.eps");
  }
  
  if (WRITE){
    TFile *myFile = new TFile("MUONplotefficiency.root", "RECREATE");
    hptyUpsilonMCAcc->Write();
    hptyUpsilonMC->Write();
    hptUpsilonMCAcc ->Write();
    hyUpsilonMCAcc->Write();

    hptyUpsilonESDAcc->Write();
    hptyUpsilonESD->Write();
    hptUpsilonESDAcc ->Write();
    hyUpsilonESDAcc->Write();

    hEffptyUpsilon->Write();
    hEffyUpsilon->Write();
    hEffptUpsilon->Write();

    hInvMass->Write();
    hInvMassUpsilonTrigger->Write();
    hUpsilonTriggerNo->Write();
    hInvMassUpsilonTriggerOneMatch->Write();
    hInvMassUpsilonTriggerTwoMatch->Write();
    hInvMassUpsilonTriggerNoMatch->Write();

    for (Int_t qw = 0 ; qw < nofMassHistograms ; qw++) 
      hInvMassInPtBins[qw]->Write();
    
    hInvMassAll->Write();
    
    grWidth->Write();
    grMean->Write();      
    
    myFile->Close();
  }
  
} // macro ends
