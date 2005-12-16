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
// comparison is done between the generated and reconstructed resonance
// Default is Upsilon but Jpsi can be studied if the ResType argument is changed
// This allows to determine several important quantities 
// reconstruction efficiency (versus pt,y), invariant mass peak variations (vs pt,y)
// reconstructed tracks and trigger tracks  matching efficiency

// Christophe Suire, IPN Orsay

void MUONplotefficiency(Int_t ResType = 553, Int_t fittype = 1){
 
  Int_t SAVE=1;
  Int_t WRITE=1;
  
  Double_t MUON_MASS = 0.105658369;
  Double_t UPSILON_MASS = 9.4603 ;
  Double_t JPSI_MASS = 3.096916;
  Double_t PI = 3.14159265358979312; 
  Double_t RESONANCE_MASS; 

  if (ResType==553)  
     RESONANCE_MASS = UPSILON_MASS;
  if (ResType==443)  
     RESONANCE_MASS = JPSI_MASS ;
    
  // 553 is the pdg code for Upsilon 
  // 443 is the pdg code for Jpsi 

  Int_t fitfunc = fittype;
  // 0 gaussian fit +/- 450 MeV/c2
  // 1 gaussian fit +/- 250 MeV/c2
  // the following are not implemented in this version
  // 2 approximated landau fit reverted 
  // 3 landau fit reverted convoluated with a gaussian 
  // 4 reverted landau fitlanUpsilon

  
  Float_t gaussianFitWidth = 0.450 ; // in GeV/c2
  if (fitfunc==1) gaussianFitWidth = 0.250 ; 
  
  // fit range : very important for LandauGauss Fit
  Float_t FitLow = 0.0  ;  Float_t FitHigh = 100. ;
  if (ResType==443) FitLow = 2.2 ; FitHigh = 3.9 ;
  if (ResType==553) FitLow = 8.5 ; FitHigh = 10.2 ;

  
  //gStyle->Reset();
  // two ways to set the stat box for all histos
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
  //("Ktuple","Kinematics NTuple","ev:npart:id:idmo:idgdmo:p:pt:y:theta:pseudorap:vx:vy:vz");
  TNtuple *ESDtuple = (TNtuple*)effFile->Get("ESDtuple");
  //("ESDtuple","Reconstructed Mu+Mu- pairs","ev:pt:y:theta:minv:pt1:y1:theta1:q1:pt2:y2:theta2:q2");


  /*********************************/
  // Histograms limits and binning
  Float_t ptmin = 0.0;      Float_t ptmax =  20.0;     Int_t ptbins = 10; 
  Float_t ymin = -4.5;       Float_t ymax =  -2.;       Int_t ybins = 10;
  Float_t thetamin = 165.;  Float_t thetamax = 180.;   Int_t thetabins = 100;  
  
  Float_t etacutmin = -4.04813 ;
  Float_t etacutmax = -2.54209 ;

  Float_t invMassMin = .0 ;   Float_t invMassMax = 15.0 ;   Int_t   invMassBins = 100 ;
  Float_t ptmin = 0.0;      Float_t ptmax =  20.0;     Int_t ptbins = 10; 
 if (ResType==443){
    invMassMin = 1.0 ;  invMassMax = 4.5 ;  
    ptmin = 0.0;     ptmax =  10.0; 
  }
  if (ResType==553){
   invMassMin =  7.0 ;  invMassMax = 11.0 ;  
   ptmin = 0.0;   ptmax =  20.0; 
  }
  
  /*********************************/
  // Values used to define the acceptance
  //
  Float_t thetacutmin = 171.0;  
  Float_t thetacutmax = 178.;
  
  Float_t ptcutmin = 0.;
  Float_t ptcutmax = 20.;

  if (ResType==443){
    ptcutmin = 0.; ptcutmax = 10.;
  }

  if (ResType==553){
    ptcutmin = 0.; ptcutmax = 20.;
  }
  

  Float_t masssigma  = 0.1; // 100 MeV/c2 is a correct estimation 
  Float_t masscutmin = 0; Float_t masscutmax = 0 ;    
  if (masssigma){
    masscutmin = RESONANCE_MASS - 3.0*masssigma ; 
    masscutmax = RESONANCE_MASS + 3.0*masssigma ;
  }
  else {
    masscutmin = RESONANCE_MASS - 1.0 ; 
    masscutmax = RESONANCE_MASS + 1.0 ;
  }


  // here no cut on theta is necesary since during the simulation 
  // gener->SetChildThetaRange(171.0,178.0);
  // thus the resonance are generated in 4 PI but only the ones for which the decay muons 
  // are in the dimuon arm acceptance
  // the pt cut is also "critical", in the generation part, resonance are generated within  the range 
  // (ptcutmin,ptcutmax). During the reconstruction, the resonance could have a pt greater than ptcutmax !!
  // Should these resonances be cut or not ?
  // probably not but they will be for now (~ 1/2000)

  // Another acceptance cut to add is a range for the  recontructed  invariant mass, it is 
  // obvious that an upsilon reconstructed with a mass  of 5 GeV/c2 is not correct. Thus
  Char_t ResonanceAccCutMC[200];  
  sprintf(ResonanceAccCutMC,"pt>= %.2f && pt<= %.2f",ptcutmin,ptcutmax);

  Char_t ResonanceAccCutESD[200];
  sprintf(ResonanceAccCutESD,"pt >=  %.2f && pt<= %.2f && minv>=  %.2f && minv<= %.2f",ptcutmin,ptcutmax,masscutmin,masscutmax);



  /*********************************/
  // Cut conditions (Id,Background...)

  Char_t IdcutResonanceMC[100];      
  Char_t IdcutMuminusMC[100];   
  Char_t IdcutMuplusMC[100];      
 
  if (ResType==553){
    sprintf(IdcutResonanceMC,"id==553");
    sprintf(IdcutMuminusMC,"id==13 && idmo==553");
    sprintf(IdcutMuplusMC,"id==-13 && idmo==553");  
  }
  if (ResType==443){
    sprintf(IdcutResonanceMC,"id==443");
    sprintf(IdcutMuminusMC,"id==13 && idmo==443");
    sprintf(IdcutMuplusMC,"id==-13 && idmo==443");  
  }
  
  //means no cuts since we don't have the trackid propagated yet pt>0
  // now we have it 10/05/2005 but it's still not being used
  
  
  // Background is calculated in the MUONefficiency.C macro 
  // Background calculation is meaningful only when Resonance have been generated with realistic background
  // when it's a pure Resonance file, the background lies artificially at the Resonance mass 
  
  //no realistic background
  Int_t REALISTIC_BACKGROUND = 0 ; 

  // to take background into account, i.e. substraction, if necessary
  // this is not ready yet
  Char_t BckgdCutResonanceESD[100];      
  sprintf(BckgdCutResonanceESD,"pt>0 && minv>%f && minv<%f && pt1>0 && pt2>0",masscutmin,masscutmax); 
  
  // if realistic background 
  // same cut + substract the background from hInvMassBg
  


  /*******************************************************************************************************************/
  Char_t txt[50];
  TLatex *tex;
  TLegend *leg;
  
  /*******************************/
  /*      Monte Carlo Part       */
  /*******************************/

  //----------------------------
  // Pt-rapidity distributions from Kinematics
  //----------------------------
  TH2F *hptyResonanceMC = new TH2F("hptyResonanceMC", " Monte Carlo Resonance",ptbins,ptmin,ptmax,ybins,ymin,ymax);
  hptyResonanceMC->SetLineColor(1);
  hptyResonanceMC->SetLineStyle(1);
  hptyResonanceMC->SetLineWidth(2);
  Ktuple->Project("hptyResonanceMC","y:pt",IdcutResonanceMC);
  hptyResonanceMC->GetXaxis()->SetTitle("P_{#perp}   [GeV/c]");
  hptyResonanceMC->GetYaxis()->SetTitle("Rapidity");
 
  cout << "  " << endl;
  cout << "********************" << endl;
  cout << " Monte Carlo Tracks  "<< endl;
  cout << " " << hptyResonanceMC->GetEntries() << " Resonance in simulation " << endl;
  

  //******** Add acceptance cuts  - Theta and Pt
  TH2F *hptyResonanceMCAcc = new TH2F("hptyResonanceMCAcc", " Monte Carlo Resonance in acceptance",ptbins,ptmin,ptmax,ybins,ymin,ymax);
  hptyResonanceMCAcc->SetLineColor(1);
  hptyResonanceMCAcc->SetLineStyle(1);
  hptyResonanceMCAcc->SetLineWidth(2);

  TString m1MC(IdcutResonanceMC);
  TString m2MC(ResonanceAccCutMC);
  TString m3MC = m1MC + " && " + m2MC ;

  Ktuple->Project("hptyResonanceMCAcc","y:pt",m3MC.Data());
  hptyResonanceMCAcc->GetYaxis()->SetTitle("Rapidity");
  hptyResonanceMCAcc->GetXaxis()->SetTitle("P_{#perp}   [GeV/c]");
  hptyResonanceMCAcc->SetStats(1);
  hptyResonanceMCAcc->Sumw2();
  
  TCanvas *c1 = new TCanvas("c1", "Resonance  Monte Carlo:  Pt vs Y",30,30,700,500);
  c1->Range(-1.69394,-0.648855,15.3928,2.77143);
  c1->SetBorderSize(2);
  c1->SetRightMargin(0.0229885);
  c1->SetTopMargin(0.0275424);
  c1->SetFrameFillColor(0);
  c1->cd();

  hptyResonanceMCAcc->SetStats(0);
  hptyResonanceMCAcc->Draw("LEGO2ZFB"); 
  //TLatex *tex = new TLatex(2,6,"#mu^{-}");
  //tex->SetTextSize(0.06);
  //tex->SetLineWidth(2);  
  //tex->Draw();
 
  sprintf(txt,"Resonance : %d entries",hptyResonanceMCAcc->GetEntries());
  tex = new TLatex(-0.854829,0.794436,txt);
  tex->SetLineWidth(2);
  tex->Draw();

  c1->Update();
  if (SAVE){
    c1->SaveAs("ptyResonanceMCAcc.gif");
    c1->SaveAs("ptyResonanceMCAcc.eps");
  }


  TH1F *hptResonanceMCAcc = new TH1F("hptResonanceMCAcc", " Monte Carlo Resonance in acceptance",ptbins,ptmin,ptmax);
  hptResonanceMCAcc->SetLineColor(1);
  hptResonanceMCAcc->SetLineStyle(1);
  hptResonanceMCAcc->SetLineWidth(2);
  Ktuple->Project("hptResonanceMCAcc","pt",m3MC.Data());
  hptResonanceMCAcc->GetXaxis()->SetTitle("P_{#perp}   [GeV/c]");
  hptResonanceMCAcc->SetStats(1);
  hptResonanceMCAcc->Sumw2();
  
  TH1F *hyResonanceMCAcc = new TH1F("hyResonanceMCAcc", " Monte Carlo Resonance in acceptance",ybins,ymin,ymax);
  hyResonanceMCAcc->SetLineColor(1);
  hyResonanceMCAcc->SetLineStyle(1);
  hyResonanceMCAcc->SetLineWidth(2);
  Ktuple->Project("hyResonanceMCAcc","y",m3MC.Data());
  hyResonanceMCAcc->GetXaxis()->SetTitle("Rapidity");
  hyResonanceMCAcc->SetStats(1);
  hyResonanceMCAcc->Sumw2();
  
  
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
  cout << " " << hptyResonanceMCAcc->GetEntries() << " Resonance in acceptance cuts " << endl;
  

  /*******************************************************************************************************************/

  /*******************************/
  /*  Reconstructed Tracks Study */
  /*******************************/

  //----------------------------
  // Pt-rapidity distributions from ESD : reconstructed tracks/particle
  //----------------------------
  TH2F *hptyResonanceESD = new TH2F("hptyResonanceESD", " Reconstucted Resonances",ptbins,ptmin,ptmax,ybins,ymin,ymax);
  hptyResonanceESD->SetLineColor(1);
  hptyResonanceESD->SetLineStyle(1);
  hptyResonanceESD->SetLineWidth(2);
  ESDtuple->Project("hptyResonanceESD","y:pt",BckgdCutResonanceESD);
  hptyResonanceESD->GetXaxis()->SetTitle("P_{#perp}   [GeV/c]");
  hptyResonanceESD->GetYaxis()->SetTitle("Rapidity");

  cout << "  " << endl;
  cout << "********************" << endl;
  cout << " Reconstructed Tracks" << endl ;
  cout << " " << hptyResonanceESD->GetEntries() << " Resonance reconstructed " << endl;

  if (REALISTIC_BACKGROUND){
    TH2F *hptyResonanceESDBck = new TH2F("hptyResonanceESDBck", "Background",ptbins,ptmin,ptmax,ybins,ymin,ymax);
    hptyResonanceESDBck->SetLineColor(1);
    hptyResonanceESDBck->SetLineStyle(1);
    hptyResonanceESDBck->SetLineWidth(2);
    ESDtupleBck->Project("hptyResonanceESDBck","y:pt",BckgdCutResonanceESD);
    hptyResonanceESDBck->GetXaxis()->SetTitle("P_{#perp}   [GeV/c]");
    hptyResonanceESDBck->GetYaxis()->SetTitle("Rapidity");
    cout << " with " << hptyResonanceESDBck->GetEntries() << " Resonances from Background (random mixing) " << endl;
  }

  // if something is wrong 
  if ( hptyResonanceESD->GetEntries()==0) {cout << " No entries in hptyResonanceESD " << endl ; break ;}

  //with Acc cuts - Theta and Pt
  TH2F *hptyResonanceESDAcc = new TH2F("hptyResonanceESDAcc", "Reconstructed Resonances",ptbins,ptmin,ptmax,ybins,ymin,ymax);
  hptyResonanceESDAcc->SetLineColor(1);
  hptyResonanceESDAcc->SetLineStyle(1);
  hptyResonanceESDAcc->SetLineWidth(2);
  
  TString m1Rec(BckgdCutResonanceESD);
  TString m2Rec(ResonanceAccCutESD);
  TString m3Rec = m1Rec + " && " + m2Rec ;

  ESDtuple->Project("hptyResonanceESDAcc","y:pt",m3Rec.Data());
  hptyResonanceESDAcc->GetYaxis()->SetTitle("Rapidity");
  hptyResonanceESDAcc->GetXaxis()->SetTitle("P_{#perp}   [GeV/c]");
  hptyResonanceESDAcc->SetStats(1);
  hptyResonanceESDAcc->Sumw2();

  cout << "  " << endl;
  cout << "********************" << endl;
  cout << " Reconstructed Tracks" << endl ;
  cout << " " << hptyResonanceESDAcc->GetEntries() << " Resonance in acceptance cuts " << endl;

  if (REALISTIC_BACKGROUND){
    //with Acc cuts - Theta and Pt
    TH2F *hptyResonanceESDBckAcc = new TH2F("hptyResonanceESDBckAcc", "Reconstructed Resonances",ptbins,ptmin,ptmax,ybins,ymin,ymax);
    hptyResonanceESDBckAcc->SetLineColor(1);
    hptyResonanceESDBckAcc->SetLineStyle(1);
    hptyResonanceESDBckAcc->SetLineWidth(2);
    
    TString m1Rec(BckgdCutResonanceESD);
    TString m2Rec(ResonanceAccCutESD);
    TString m3Rec = m1Rec + " && " + m2Rec ;
    
    ESDtupleBck->Project("hptyResonanceESDBckAcc","y:pt",m3Rec.Data());
    hptyResonanceESDBckAcc->GetYaxis()->SetTitle("Rapidity");
    hptyResonanceESDBckAcc->GetXaxis()->SetTitle("P_{#perp}   [GeV/c]");
    hptyResonanceESDBckAcc->SetStats(1);
    hptyResonanceESDBckAcc->Sumw2();
    cout << " with " << hptyResonanceESDBckAcc->GetEntries() << " Resonances from Background (random mixing) " << endl;
  }


  TCanvas *c100 = new TCanvas("c100", "Resonance Reconstructed in Acceptance: Pt vs Y",210,30,700,500);
  c100->Range(-1.69394,-0.648855,15.3928,2.77143);
  c100->SetBorderSize(2);
  c100->SetRightMargin(0.0229885);
  c100->SetTopMargin(0.0275424);
  c100->SetFrameFillColor(0);
  
  c100->cd();
  hptyResonanceESDAcc->SetStats(0);
  hptyResonanceESDAcc->Draw("LEGO2ZFB"); 
  sprintf(txt,"Resonance : %d entries",hptyResonanceESDAcc->GetEntries());
  tex = new TLatex(-0.854829,0.794436,txt);
  tex->SetLineWidth(2);
  tex->Draw();

  c100->Update();
  if (SAVE){
    c100->SaveAs("ptyResonanceESDAcc.gif");
    c100->SaveAs("ptyResonanceESDAcc.eps");
  }

  if (REALISTIC_BACKGROUND){
    TCanvas *c110 = new TCanvas("c110", "Resonance Background Reconstructed in Acceptance: Pt vs Y",215,35,700,500);
    c110->Range(-1.69394,-0.648855,15.3928,2.77143);
    c110->SetBorderSize(2);
    c110->SetRightMargin(0.0229885);
    c110->SetTopMargin(0.0275424);
    c110->SetFrameFillColor(0);
    
    c110->cd();
    hptyResonanceESDBckAcc->SetStats(0);
    hptyResonanceESDBckAcc->Draw("LEGO2ZFB"); 
    sprintf(txt,"Resonance backround : %d entries",hptyResonanceESDBckAcc->GetEntries());
    tex = new TLatex(-0.854829,0.794436,txt);
    tex->SetLineWidth(2);
    tex->Draw();
    
    
    c110->Update();
    if (SAVE){
      c110->SaveAs("ptyResonanceESDBckAcc.gif");
      c110->SaveAs("ptyResonanceESDBckAcc.eps");
    }
  }
  
  TH1F *hptResonanceESDAcc = new TH1F("hptResonanceESDAcc", " Monte Carlo Resonance in acceptance",ptbins,ptmin,ptmax);
  hptResonanceESDAcc->SetLineColor(1);
  hptResonanceESDAcc->SetLineStyle(1);
  hptResonanceESDAcc->SetLineWidth(2);
  ESDtuple->Project("hptResonanceESDAcc","pt",m3Rec.Data());
  hptResonanceESDAcc->GetXaxis()->SetTitle("P_{#perp}   [GeV/c]");
  hptResonanceESDAcc->SetStats(1);
  hptResonanceESDAcc->Sumw2();
  
  TH1F *hyResonanceESDAcc = new TH1F("hyResonanceESDAcc", " Monte Carlo Resonance in acceptance",ybins,ymin,ymax);
  hyResonanceESDAcc->SetLineColor(1);
  hyResonanceESDAcc->SetLineStyle(1);
  hyResonanceESDAcc->SetLineWidth(2);
  ESDtuple->Project("hyResonanceESDAcc","y",m3Rec.Data());
  hyResonanceESDAcc->GetXaxis()->SetTitle("Rapidity");
  hyResonanceESDAcc->SetStats(1);
  hyResonanceESDAcc->Sumw2();
  

  /*******************************************************************************************************************/

  /*******************************/
  /* Efficiencies calculations   */
  /*******************************/
  

  cout << "  " << endl;
  cout << "***********************" << endl;
  cout << " Integrated efficiency" << endl ;
  
  if (REALISTIC_BACKGROUND)
    cout << " ResonanceESDAcc/ResonanceMCAcc = " << (hptyResonanceESDAcc->GetEntries()-hptyResonanceESDBckAcc->GetEntries())/hptyResonanceMCAcc->GetEntries()  << endl;
  else 
    cout << " ResonanceESDAcc/ResonanceMCAcc = " << hptyResonanceESDAcc->GetEntries()/hptyResonanceMCAcc->GetEntries()  << endl;

  TH2F *hEffptyResonance = new TH2F("hEffptyResonance", " Resonance Efficiency",ptbins,ptmin,ptmax,ybins,ymin,ymax);

  if (REALISTIC_BACKGROUND){
    TH2F *hptyResonanceESDAccBckSubstracted = new TH2F("hptyResonanceESDAccBckSubstracted","hptyResonanceESDAccBckSubstracted",ptbins,ptmin,ptmax,ybins,ymin,ymax);
    hptyResonanceESDAccBckSubstracted->Add(hptyResonanceESDAcc,hptyResonanceESDBckAcc,1,-1);
    hEffptyResonance->Divide(hptyResonanceESDAccBckSubstracted,hptyResonanceMCAcc,1,1);
  }
  else 
  hEffptyResonance->Divide(hptyResonanceESDAcc,hptyResonanceMCAcc,1,1);

  TCanvas *c1000 = new TCanvas("c1000", "Resonance efficiency : Pt vs Y",390,30,700,500);
  c1000->Range(-1.69394,-0.648855,15.3928,2.77143);
  c1000->SetBorderSize(2);
  c1000->SetRightMargin(0.0229885);
  c1000->SetTopMargin(0.0275424);
  c1000->SetFrameFillColor(0);

  c1000->cd();
  hEffptyResonance->SetStats(0);
  hEffptyResonance->Draw("LEGO2fz");



  TH1F *hEffptResonance = new TH1F("hEffptResonance", "Resonance Efficiency vs pt",ptbins,ptmin,ptmax);
  hEffptResonance->Divide(hptResonanceESDAcc,hptResonanceMCAcc,1,1);
  hEffptResonance->SetLineWidth(2);
  hEffptResonance->SetMinimum(0);
  hEffptResonance->SetStats(1);


  TCanvas *c1100 = new TCanvas("c1100", "Resonance efficiency : Pt",410,50,700,500);
  c1100->Range(-1.69394,-0.648855,15.3928,2.77143);
  c1100->SetBorderSize(2);
  c1100->SetRightMargin(0.0229885);
  c1100->SetTopMargin(0.0275424);
  c1100->SetFrameFillColor(0);

  c1100->cd();
  hEffptResonance->SetStats(0);
  hEffptResonance->GetXaxis()->SetTitle("P_{#perp}  [GeV/c]");
  hEffptResonance->GetYaxis()->SetTitle("Efficiency ");
  hEffptResonance->Draw("E");
   
  
  
  TH1F *hEffyResonance = new TH1F("hEffyResonance", "Resonance Efficiency vs y",ybins,ymin,ymax);
  hEffyResonance->Divide(hyResonanceESDAcc,hyResonanceMCAcc,1,1);
  hEffyResonance->SetLineWidth(2);
  hEffyResonance->SetStats(1);
  
  TCanvas *c1200 = new TCanvas("c1200", "Resonance efficiency : Y",430,70,700,500);
  c1200->Range(-1.69394,-0.648855,15.3928,2.77143);
  c1200->SetBorderSize(2);
  c1200->SetRightMargin(0.0229885);
  c1200->SetTopMargin(0.0275424);
  c1200->SetFrameFillColor(0);

  c1200->cd();
  hEffyResonance->SetStats(0);
  hEffyResonance->GetXaxis()->SetTitle("Rapidity");
  hEffyResonance->GetYaxis()->SetTitle("Efficiency ");
  hEffyResonance->Draw("E");

  c1000->Update();
  c1100->Update();
  c1200->Update();
  if (SAVE){
    c1000->SaveAs("EffptyResonance.gif");
    c1000->SaveAs("EffptyResonance.eps");
    c1100->SaveAs("EffptResonance.gif");
    c1100->SaveAs("EffptResonance.eps");
    c1200->SaveAs("EffyResonance.gif");
    c1200->SaveAs("EffyResonance.eps");
  }

 /*******************************************************************************************************************/

  /*******************************/
  /* Trigger matching            */
  /*******************************/


  Float_t triggerChi2Min = 0.; 
  Float_t triggerChi2Max = 7.5;

  //TString m1Rec(BckgdCutResonanceESD);
  //TString m2Rec(ResonanceAccCutESD);
  TString  m1TG = m1Rec + " && " + m2Rec ;
  //cout << m1TG.Data() << endl;
  
  TH1F *hInvMassNoTriggerCut= new TH1F("hInvMassNoTriggerCut","hInvMassNoTriggerCut",invMassBins,invMassMin,invMassMax);
  hInvMassNoTriggerCut->GetXaxis()->SetTitle("Mass [GeV/c^{2}]");
  hInvMassNoTriggerCut->GetYaxis()->SetTitle("Counts");
  ESDtuple->Project("hInvMassNoTriggerCut","minv",m1TG.Data());

  //Add trigger UnlikePairAllPt
  TString  m2TG = m1TG + " && (tw & 0x800) == 2048" ;
  TH1F *hInvMassResonanceTrigger= new TH1F("hInvMassResonanceTrigger","hInvMassResonanceTrigger",invMassBins,invMassMin,invMassMax);
  hInvMassResonanceTrigger->GetXaxis()->SetTitle("Mass [GeV/c^{2}]");
  hInvMassResonanceTrigger->GetYaxis()->SetTitle("Counts");
  ESDtuple->Project("hInvMassResonanceTrigger","minv",m2TG.Data());


  cout << "  " << endl;
  cout << "********************" << endl;
  cout << " TriggerMatching" << endl ;
  cout << " " << hInvMassResonanceTrigger->GetEntries() << " Resonance with trigger UnlikePairAllPt " << endl;

  TString  m2TGNo = m1TG + " && (tw & 0x800) != 2048" ;
  TH1F *hResonanceTriggerNo= new TH1F("hResonanceTriggerNo","hResonanceTriggerNo",32769,0,32768);
  hResonanceTriggerNo->GetXaxis()->SetTitle("Mass [GeV/c^{2}]");
  hResonanceTriggerNo->GetYaxis()->SetTitle("Counts");
  ESDtuple->Project("hResonanceTriggerNo","tw",m2TGNo.Data());

  


  //Add matching rec/trig for 2 tracks 
  TString m3TG = m2TG + " && trig1 > 0 && trig2 > 0" ;  
  TH1F *hInvMassResonanceTriggerTwoMatch = new TH1F("hInvMassResonanceTriggerTwoMatch","hInvMassResonanceTrigger with 2 matching tracks",invMassBins,invMassMin,invMassMax);
  hInvMassResonanceTriggerTwoMatch->GetXaxis()->SetTitle("Mass [GeV/c^{2}]");
  hInvMassResonanceTriggerTwoMatch->GetYaxis()->SetTitle("Counts");
  ESDtuple->Project("hInvMassResonanceTriggerTwoMatch","minv",m3TG.Data());


  //Add matching rec/trig for 1  tracks 
  TString m4TG = m2TG + " && (trig1 > 0 || trig2 > 0)" ;  
  TH1F *hInvMassResonanceTriggerOneMatch= new TH1F("hInvMassResonanceTriggerOneMatch","hInvMassResonanceTrigger with 1 matching tracks",invMassBins,invMassMin,invMassMax);
  hInvMassResonanceTriggerOneMatch->GetXaxis()->SetTitle("Mass [GeV/c^{2}]");
  hInvMassResonanceTriggerOneMatch->GetYaxis()->SetTitle("Counts");
  ESDtuple->Project("hInvMassResonanceTriggerOneMatch","minv",m4TG.Data());

  TString m5TG = m2TG + " && (trig1 == 0 && trig2 == 0)" ;  
  TH1F *hInvMassResonanceTriggerNoMatch= new TH1F("hInvMassResonanceTriggerNoMatch","hInvMassResonanceTrigger with no matching tracks",invMassBins,invMassMin,invMassMax);
  hInvMassResonanceTriggerNoMatch->GetXaxis()->SetTitle("Mass [GeV/c^{2}]");
  hInvMassResonanceTriggerNoMatch->GetYaxis()->SetTitle("Counts");
  ESDtuple->Project("hInvMassResonanceTriggerNoMatch","minv",m5TG.Data());

  TCanvas *c2 = new TCanvas("c2", "Resonance Trigger efficiency",30,90,700,500);
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
  hInvMassResonanceTrigger->SetLineWidth(2);
  hInvMassResonanceTrigger->SetLineStyle(0);
  hInvMassResonanceTrigger->SetLineColor(2);
  hInvMassResonanceTrigger->Draw("same");
  hInvMassResonanceTriggerTwoMatch->SetLineWidth(2);
  hInvMassResonanceTriggerTwoMatch->SetLineStyle(1);
  hInvMassResonanceTriggerTwoMatch->SetLineColor(4);
  hInvMassResonanceTriggerTwoMatch->Draw("same");
  hInvMassResonanceTriggerOneMatch->SetLineWidth(2);
  hInvMassResonanceTriggerOneMatch->SetLineStyle(2);
  hInvMassResonanceTriggerOneMatch->SetLineColor(51);
  hInvMassResonanceTriggerOneMatch->Draw("same");
  hInvMassResonanceTriggerNoMatch->SetLineWidth(2);
  hInvMassResonanceTriggerNoMatch->SetLineStyle(2);
  hInvMassResonanceTriggerNoMatch->SetLineColor(46);
  hInvMassResonanceTriggerNoMatch->Draw("same");

  TLegend *leg = new TLegend(0.12,0.6,0.50,0.89);
  leg->SetHeader("Reconstructed Resonance Invariant Mass");

  leg->AddEntry(hInvMassNoTriggerCut,Form("All  (%.0f cnts)",hInvMassNoTriggerCut->GetEntries()),"l");
  leg->AddEntry(hInvMassResonanceTrigger,Form("UnlikePairAllPt Trigger (%.0f cnts)",hInvMassResonanceTrigger->GetEntries()),"l");
  leg->AddEntry(hInvMassResonanceTriggerTwoMatch,Form("UPAllPt Trig. and 2 tracks matches (%.0f cnts)",hInvMassResonanceTriggerTwoMatch->GetEntries()),"l");  
  leg->AddEntry(hInvMassResonanceTriggerOneMatch,Form("UPAllPt Trig. and 1 track match (%.0f cnts)",hInvMassResonanceTriggerOneMatch->GetEntries()),"l");  
  leg->AddEntry(hInvMassResonanceTriggerNoMatch,Form("UPAllPt Trig. and no matching track (%.0f cnts)",hInvMassResonanceTriggerNoMatch->GetEntries()),"l");  
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
  // cs invMassMin = 7.0 ;
  // cs invMassMax = 11.0 ;
  // cs invMassBins = 100 ;
  
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
  
  if (fitfunc==0 || fitfunc==1){
    fitFunc = new TF1("gaussian","gaus(0)",FitLow,FitHigh);
    if (!fitFunc)
      fitFunc = new TF1("gaussian","[0]*TMath::Exp(-0.5*(x-[1])*(x-[1])/[2]/[2])",FitLow,FitHigh);
    fitFunc->SetParNames("Constant","Mean value","Width");
    fitFunc->SetLineColor(2);
    fitFunc->SetLineWidth(2);
    fitFunc->SetLineStyle(2);
    
    Float_t  *tParams = new Float_t[3] ;
    if (ResType==553){
      tParams[0] = 500;
      tParams[1] = 9.47;
      tParams[2] = 0.1;
    }
    if (ResType==443){
      tParams[0] = 500;
      tParams[1] = 3.1;
      tParams[2] = 0.1;
    }

    fitFunc->SetParameters(tParams[0],tParams[1],tParams[2]);
    
    TString FitFuncName = "gaussian";  
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



  TCanvas *call = new TCanvas("call", "Resonance : invariant mass spectra",30,330,700,500);
  call->Range(-1.69394,-0.648855,15.3928,2.77143);
  call->SetBorderSize(2);
  //call->SetRightMargin(0.2229885);
  call->SetRightMargin(0.0229885);
  call->SetTopMargin(0.0275424);
  call->SetFrameFillColor(0);
  call->cd();

  TH1F* hInvMass = new TH1F("hInvMass","Inv. Mass Pt,Y integrated",30,0,3+RESONANCE_MASS);
  ESDtuple->Project("hInvMass","minv");
  hInvMass->SetLineWidth(2);
  hInvMass->Draw("HE");

  if(REALISTIC_BACKGROUND){  
    TH1F* hInvMassBck = new TH1F("hInvMassBck","Background Pt,Y integrated",30,0,3+RESONANCE_MASS);
    ESDtupleBck->Project("hInvMassBck","minv");
    hInvMassBck->SetLineWidth(2);
    hInvMassBck->SetLineStyle(2);
    hInvMassBck->SetLineColor(2);    
    hInvMassBck->Draw("same");
  }
  call->Modified();
  call->Update();


  TCanvas *cfit = new TCanvas("cfit", "Resonance : invariant mass fit",30,330,700,500);
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
    cfit->SaveAs("ResonanceMass.gif");
    cfit->SaveAs("ResonanceMass.eps");
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
   
   TCanvas *cptfits = new TCanvas("cptfits", "Resonance : invariant mass fit",30,330,700,500);
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
    if ( fitResultsMean[qw] > invMassMin &&   fitResultsMean[qw] < invMassMax &&  fitResultsWidth[qw] > 0 &&  fitResultsWidth[qw] < 1){
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
  cout << " Integrated efficiency (+/- "<< 3*masssigma << " GeV around Resonance mass cut)" << endl ;
  cout << " ResonanceESDAcc/ResonanceMCAcc = " << hptyResonanceESDAcc->GetEntries() << "/" << hptyResonanceMCAcc->GetEntries()  << " = "  << hptyResonanceESDAcc->GetEntries()/hptyResonanceMCAcc->GetEntries()  <<endl;
  cout << "***********************" << endl;
  cout << "  " << endl;


  cout << "  " << endl;
  cout << "***********************" << endl;
  cout << " Trigger  efficiency" << endl ;
  cout << " Two muons matching = " << hInvMassResonanceTriggerTwoMatch->GetEntries() << "/" << hInvMassResonanceTrigger->GetEntries() << " = " <<  hInvMassResonanceTriggerTwoMatch->GetEntries()/hInvMassResonanceTrigger->GetEntries() << endl;
  cout << " Single muon matching  = " << hInvMassResonanceTriggerOneMatch->GetEntries() << "/" << hInvMassResonanceTrigger->GetEntries() << " = " <<  hInvMassResonanceTriggerOneMatch->GetEntries()/hInvMassResonanceTrigger->GetEntries() << endl;
  cout << " No matching  = " << hInvMassResonanceTriggerNoMatch->GetEntries() << "/" << hInvMassResonanceTrigger->GetEntries() << " = " <<  hInvMassResonanceTriggerNoMatch->GetEntries()/hInvMassResonanceTrigger->GetEntries() <<endl;
  cout << "***********************" << endl;
  cout << "  " << endl;

  
  

  TH2F *hFrame = new TH2F("hFrame", "A hFrame",ptbins,ptmin,ptmax,100,fitResultsWidth[1]-fitResultsWidth[1]*2,fitResultsWidth[1]+fitResultsWidth[1]*2);
  hptyResonanceMCAcc->SetLineColor(1);
  hptyResonanceMCAcc->SetLineStyle(1);
  hptyResonanceMCAcc->SetLineWidth(2);
  
  
  TCanvas *cA = new TCanvas("cA", "Resonance : Width from fit",30,330,700,500);
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
  grWidth->SetTitle("Resonance Mass Width");
  grWidth->SetMarkerColor(4);
  grWidth->SetMarkerStyle(21);
  grWidth->Draw("P");
  
  
  
  TCanvas *cB = new TCanvas("cB", "Resonance : Mean from fit",50,350,700,500);
  cB->Range(-1.69394,-0.648855,15.3928,2.77143);
  cB->SetBorderSize(2);
  cB->SetRightMargin(0.0229885);
  cB->SetTopMargin(0.0275424);
  cB->SetFrameFillColor(0);
  cB->cd();
  
  //TH2F *hFrameB = new TH2F("hFrameB", "A hFrameB",ptbins,ptmin,ptmax,100,fitResultsMean[1]-fitResultsMean[1]*2,fitResultsMean[1]+fitResultsMean[1]*2);
  TH2F *hFrameB = new TH2F("hFrameB", "A hFrameB",ptbins,ptmin,ptmax,100,invMassMin,invMassMax);
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
  grMean->SetTitle("Resonance Mass Mean");
  grMean->SetMarkerColor(4);
  grMean->SetMarkerStyle(21);
  grMean->Draw("P");
  
  cA->Update();
  cB->Update();
  if (SAVE){
    cB->SaveAs("ResonanceMeanVsPt.gif");
    cB->SaveAs("ResonanceMeanVsPt.eps");
    cA->SaveAs("ResonanceWidthVsPt.gif");
    cA->SaveAs("ResonanceWidthVsPt.eps");
  }
  
  if (WRITE){
    TFile *myFile = new TFile("MUONplotefficiency.root", "RECREATE");
    hptyResonanceMCAcc->Write();
    hptyResonanceMC->Write();
    hptResonanceMCAcc ->Write();
    hyResonanceMCAcc->Write();

    hptyResonanceESDAcc->Write();
    hptyResonanceESD->Write();
    hptResonanceESDAcc ->Write();
    hyResonanceESDAcc->Write();

    hEffptyResonance->Write();
    hEffyResonance->Write();
    hEffptResonance->Write();

    hInvMass->Write();
    hInvMassResonanceTrigger->Write();
    hResonanceTriggerNo->Write();
    hInvMassResonanceTriggerOneMatch->Write();
    hInvMassResonanceTriggerTwoMatch->Write();
    hInvMassResonanceTriggerNoMatch->Write();

    for (Int_t qw = 0 ; qw < nofMassHistograms ; qw++) 
      hInvMassInPtBins[qw]->Write();
    
    hInvMassAll->Write();
    
    myFile->Close();
  }
  
} // macro ends
