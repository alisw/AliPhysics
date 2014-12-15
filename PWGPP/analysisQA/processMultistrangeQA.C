//////////////////////////////////////////////////
//
//  This macro was written by Domenico Colella (domenico.colella@cern.ch).
//  12 November 2013
//
//   ------ Arguments
//   --  icasType          =  0) Xi- 1) Xi+ 2) Omega- 3) Omega+
//   --  collidingsystem   =  0) PbPb  1) pp
//   --  fileDir           =  "Input file directory"
//   --  filein            =  "Input file name"
//
//   ------ QATask output content
//   The output produced by the QATask is a CFContainer with 4 steps and 21 variables.
//   The meaning of each variable within the container are listed here:
//   --  0   = Max DCA Cascade Daughters                 pp: 2.0     PbPb: 0.3
//   --  1   = Min DCA Bach To PV                        pp: 0.01    PbPb: 0.03
//   --  2   = Min Cascade Cosine Of PA                  pp: 0.98    PbPb: 0.999
//   --  3   = Min Cascade Radius Fid. Vol.              pp: 0.2     PbPb: 0.9
//   --  4   = Window Invariant Mass Lambda              pp: 0.008   PbPb: 0.0008
//   --  5   = Max DCA V0 Daughters                      pp: 1.5     PbPb: 1.0
//   --  6   = Min V0 Cosine Of PA To PV                 pp: pT dep. PbPb: 0.98
//   --  7   = Min V0 Radius Fid. Vol.                   pp: 0.2     PbPb: 0.9
//   --  8   = Min DCA V0 To PV                          pp: 0.01    PbPb: 0.05
//   --  9   = Min DCA Pos To PV                         pp: 0.05    PbPb: 0.1
//   --  10  = Min DCA Neg To PV                         pp: 0.05    PbPb: 0.1
//   --  11  = Invariant Mass distribution for Xi
//   --  12  = Invariant Mass distribution for Omega
//   --  13  = Transverse Momentum distribution
//   --  14  = Rapidity distribution for Xi
//   --  15  = Rapidity distribution for Omega
//   --  16  = Proper length distribution for the cascade
//   --  17  = Proper length distribution for the V0
//   --  18  = Min V0 Cosine Of PA To Xi Vertex         pp: pT dep. PbPb: pT dep.
//   --  19  = Centrality
//   --  20  = ESD track multiplicity
//   The last two variables are empty in the case proton-proton collisions.
//
//   ------ Present Macro Check
//   With this macro are produced the plots of the distributions for the topological
//   variables used during the reconstruction of the cascades and defined in the two
//   classes AliCascadeVertexer.cxx and AliV0vertexer.cxx contained in /STEER/ESD/ 
//   folder in Aliroot.
//
//   -- First Canvas: DCA cascade daughters, Bachelor IP to PV, Cascade cosine of PA
//                    Cascade radius of fiducial volume, Invariant mass Lambda,
//                    DCA V0 daughters.
//   -- Second Canvas: V0 cosine of PA to PV, Min V0 Radius Fid. Vol., Min DCA V0 To PV
//                     Min DCA Pos To PV, Min DCA Neg To PV, V0 cosine of PA to XiV
//
//   In this plots, in correspondence to the min/max cut value adopted in the reconstruction
//   a line is drawn, to check if there is same problem in the cuts definition.
//
//   Also, other specific distribution for the selected cascades are ploted as: the
//   invariant mass distribution, the transverse momentum distribution, the rapidity
//   distribution, proper length distribution for the cascade and the v0.
//
//   -- Third Canvas: InvMass, Transverse momentum, Cascade proper length, V0 proper length
//
//   In the end, only for thr PbPb collision the centrality distribution and the
//   track multiplicity distribution are sored.
//
//   -- Fourth Canvas: Centrality, track multiplicity
//
//
//////////////////////////////////////////////////////




class AliCFContainer;

void processMultistrangeQA(Int_t   icasType        = 0,                        // 0) Xi- 1) Xi+ 2) Omega- 3) Omega+
			   Int_t   collidingsystem = 0,                        // 0) PbPb  1) pp
			   Char_t *fileDir         = "./",                     // Input file directory
			   Char_t *filein          = "AnalysisResults.root",   // Input file name
			   TString suffix          = "eps",
			   const char *outfile     = "output.root"             // output
			   ) {


      /////////////
      gStyle->SetOptStat(1110);
      gStyle->SetOptStat(kFALSE);
      gStyle->SetOptTitle(kFALSE);
      gStyle->SetFrameLineWidth(2.5);
      gStyle->SetCanvasColor(0);
      gStyle->SetPadColor(0);
      gStyle->SetHistLineWidth(2.5);
      gStyle->SetLabelSize(0.05, "x");
      gStyle->SetLabelSize(0.05, "y");
      gStyle->SetTitleSize(0.05, "x");
      gStyle->SetTitleSize(0.05, "y");
      gStyle->SetTitleOffset(1.1, "x");
      gStyle->SetPadBottomMargin(0.14);
      gSystem->Load("libANALYSIS");
      gSystem->Load("libANALYSISalice");
      gSystem->Load("libCORRFW");

 

     TFile *f1 = new TFile(Form("%s/%s",fileDir,filein));
     AliCFContainer *cf = (AliCFContainer*) (f1->Get("PWGLFStrangeness.outputCheckCascade/fCFContCascadeCuts"));
     if(!cf) {
       Printf("FATAL: PWGLFStrangeness.outputCheckCascade/fCFContCascadeCuts ===> Not Available");
       Printf("Exiting processMultistrangeQA");
       return;
     }


     //DEEFINE TEXT
     TLatex* t1 = new TLatex(0.6,0.55,"#color[3]{OK!!}");
     t1->SetTextSize(0.1);
     t1->SetNDC();
     TLatex* t2 = new TLatex(0.6,0.55,"#color[2]{NOT OK!!}");
     t2->SetTextSize(0.1);
     t2->SetNDC();
     t2->SetTextColor(2);
     TLatex* tcasc;
     if      (icasType == 0) tcasc = new TLatex(0.8,0.7,"#color[1]{#Xi^{-}}");
     else if (icasType == 1) tcasc = new TLatex(0.8,0.7,"#color[1]{#Xi^{+}}");
     else if (icasType == 2) tcasc = new TLatex(0.8,0.7,"#color[1]{#Omega^{-}}");
     else if (icasType == 3) tcasc = new TLatex(0.8,0.7,"#color[1]{#Omega^{+}}");
     tcasc->SetTextSize(0.2);
     tcasc->SetNDC();
     tcasc->SetTextColor(2);
     TLatex* tpdgmass;
     if      (icasType == 0) tpdgmass = new TLatex(0.55,0.7,"#color[1]{PDG mass: 1.321 GeV/c^{2}}");
     else if (icasType == 1) tpdgmass = new TLatex(0.55,0.7,"#color[1]{PDG mass: 1.321 GeV/c^{2}}");
     else if (icasType == 2) tpdgmass = new TLatex(0.55,0.7,"#color[1]{PDG mass: 1.672 GeV/c^{2}}");
     else if (icasType == 3) tpdgmass = new TLatex(0.55,0.7,"#color[1]{PDG mass: 1.672 GeV/c^{2}}");
     tpdgmass->SetTextSize(0.07);
     tpdgmass->SetNDC();
     tpdgmass->SetTextColor(2);
 
     // Added by sjena
     TFile *fout = TFile::Open(outfile,"UPDATE");
     fout->ls();
     
     TDirectoryFile *cdd = NULL;
     cdd = (TDirectoryFile*)fout->Get("LF");
     if(!cdd) {
       Printf("Warning: LF <dir> doesn't exist, creating a new one");
       cdd = (TDirectoryFile*)fout->mkdir("LF");
     }
     cdd->cd();
     cdd->ls();
     


     // TFile *f = TFile::Open("OutPut.root","UPDATE");


     //DEFINE 1st CANVAS AND DRAW PLOTS
     TCanvas *c1 = new TCanvas("c1","",1200,800);
     c1->Divide(2,3); 
       //Pad 1: DCA cascade daughters
       c1->cd(1);
       gPad->SetLogy();
       TH1D *hvar0 = cf->ShowProjection(0,icasType);
       hvar0->Draw("histo");

       //   hvar0->SetName(Form("fig_lf_multistrange_%s", hvar0->GetName()));
       hvar0->SetName(Form("fig_lf_multistrange_0"));
       hvar0->Write();



       Double_t x0;
       if      (collidingsystem == 0) x0 = 0.3;
       else if (collidingsystem == 1) x0 = 2.0;
       TLine *line0 = new TLine(x0,0.,x0,hvar0->GetBinContent(hvar0->GetMaximumBin()));
       line0->SetLineColor(kRed);
       line0->SetLineStyle(9);
       line0->SetLineWidth(2.0);
       line0->Draw("same");
          Bool_t check_0 = checkOverTheLimit(hvar0, x0);
          if (check_0) { cout<<"The cut is OK!!"<<endl; t1->Draw(); }
          else         { cout<<"The cut is NOT OK!!"<<endl; t2->Draw(); }
       tcasc->Draw();
       //Pad 2: Bachelor IP to PV
       c1->cd(2);
       gPad->SetLogy();
       TH1D *hvar1 = cf->ShowProjection(1,icasType);
       hvar1->GetXaxis()->SetRangeUser(0.,0.24);
       hvar1->Draw("histo");

       //   hvar1->SetName(Form("fig_lf_multistrange_%s", hvar1->GetName()));
       hvar1->SetName(Form("fig_lf_multistrange_1"));
       hvar1->Write();


       Double_t x1;
       if      (collidingsystem == 0) x1 = 0.03;
       else if (collidingsystem == 1) x1 = 0.01;
       TLine *line1 = new TLine(x1,0.,x1,hvar1->GetBinContent(hvar1->GetMaximumBin()));
       line1->SetLineColor(kRed);
       line1->SetLineStyle(9);
       line1->SetLineWidth(2.0);
       line1->Draw("same");
          Bool_t check_1 = checkUnderTheLimit(hvar1, x1);
          if (check_1) { cout<<"The cut is OK!!"<<endl; t1->Draw(); }
          else         { cout<<"The cut is NOT OK!!"<<endl; t2->Draw(); }
       //Pad 3: Cascade cosine of Pointing Angle
       c1->cd(3);
       gPad->SetLogy();
       TH1D *hvar2 = cf->ShowProjection(2,icasType);
       Double_t max2 = hvar2->GetBinContent(hvar2->GetMaximumBin());
       hvar2->GetYaxis()->SetRangeUser(0.01,max2*1.5);
       hvar2->Draw("histo");

       //     hvar2->SetName(Form("fig_lf_multistrange_%s", hvar2->GetName()));
       hvar2->SetName(Form("fig_lf_multistrange_2"));
       hvar2->Write();
       

       Double_t x2;
       if      (collidingsystem == 0) x2 = 0.999;
       else if (collidingsystem == 1) x2 = 0.98;
       TLine *line2 = new TLine(x2,0.,x2,hvar2->GetBinContent(hvar2->GetMaximumBin()));
       line2->SetLineColor(kRed);
       line2->SetLineStyle(9);
       line2->SetLineWidth(2.0);
       line2->Draw("same");
       line1->Draw("same");
          Bool_t check_2 = checkUnderTheLimit(hvar2, x2);
          if (check_2) { cout<<"The cut is OK!!"<<endl; t1->Draw(); }
          else         { cout<<"The cut is NOT OK!!"<<endl; t2->Draw(); }
       //Pad 4: Cascade radius of fiducial volume
       c1->cd(4);
       gPad->SetLogy();
       TH1D *hvar3 = cf->ShowProjection(3,icasType);
       hvar3->GetXaxis()->SetRangeUser(0.,3.8);
       hvar3->Draw("histo");

       //   hvar3->SetName(Form("fig_lf_multistrange_%s", hvar3->GetName()));
       hvar3->SetName(Form("fig_lf_multistrange_3"));
       hvar3->Write();

       Double_t x3;
       if      (collidingsystem == 0) x3 = 0.9;
       else if (collidingsystem == 1) x3 = 0.2;
       TLine *line3 = new TLine(x3,0.,x3,hvar3->GetBinContent(hvar3->GetMaximumBin()));
       line3->SetLineColor(kRed);
       line3->SetLineStyle(9);
       line3->SetLineWidth(2.0);
       line3->Draw("same");
          Bool_t check_3 = checkUnderTheLimit(hvar3, x3);
          if (check_3) { cout<<"The cut is OK!!"<<endl; t1->Draw(); }
          else         { cout<<"The cut is NOT OK!!"<<endl; t2->Draw(); }
       //Pad 5: Invariant mass Lambda
       c1->cd(5);
       TH1D *hvar4 = cf->ShowProjection(4,icasType);
       hvar4->Draw("histo");

       hvar4->SetName(Form("fig_lf_multistrange_4"));
       hvar4->Write();

       Double_t x41 = 1.116 + 0.008;
       TLine *line41 = new TLine(x41,0.,x41,hvar4->GetBinContent(hvar4->GetMaximumBin()));
       line41->SetLineColor(kRed);
       line41->SetLineStyle(9);
       line41->SetLineWidth(2.0);
       line41->Draw("same");
       Double_t x42 = 1.115 - 0.008;
       TLine *line42 = new TLine(x42,0.,x42,hvar4->GetBinContent(hvar4->GetMaximumBin()));
       line42->SetLineColor(kRed);
       line42->SetLineStyle(9);
       line42->SetLineWidth(2.0);
       line42->Draw("same");
          Bool_t check_4_1 = checkUnderTheLimit(hvar3, x3);
          Bool_t check_4_2 = checkOverTheLimit(hvar0, x0);
          if (check_4_1 && check_4_2) { cout<<"The cut is OK!!"<<endl; t1->Draw(); }
          else                        { cout<<"The cut is NOT OK!!"<<endl; t2->Draw(); }
       //Pad 6: DCA V0 daughters
       c1->cd(6);
       gPad->SetLogy();
       TH1D *hvar5 = cf->ShowProjection(5,icasType);
       hvar5->Draw("histo");

       //   hvar5->SetName(Form("fig_lf_multistrange_%s", hvar5->GetName()));
       hvar5->SetName(Form("fig_lf_multistrange_5"));
       hvar5->Write();

       Double_t x5;
       if      (collidingsystem == 0) x5 = 1.0;
       else if (collidingsystem == 1) x5 = 1.5;
       TLine *line5 = new TLine(x5,0.,x5,hvar5->GetBinContent(hvar5->GetMaximumBin()));
       line5->SetLineColor(kRed);
       line5->SetLineStyle(9);
       line5->SetLineWidth(2.0);
       line5->Draw("same");
          Bool_t check_5 = checkOverTheLimit(hvar5, x5);
          if (check_5) { cout<<"The cut is OK!!"<<endl; t1->Draw(); }
          else         { cout<<"The cut is NOT OK!!"<<endl; t2->Draw(); }

	  c1->SaveAs(Form("fig_lf_MultistrangeQA_1.%s",suffix.Data()));
	  
    

     //DEFINE 2st CANVAS AND DRAW PLOTS
     TCanvas *c2 = new TCanvas("c2","",1200,800);
     c2->Divide(2,3);
       //Pad 1: V0 cosine of Pointing Angle to PV
       c2->cd(1);
       gPad->SetLogy();
       TH1D *hvar6 = cf->ShowProjection(6,icasType);
       Double_t max6 = hvar6->GetBinContent(hvar6->GetMaximumBin());
       hvar6->GetYaxis()->SetRangeUser(0.01,max6*1.5);
       hvar6->Draw("histo");

       //  hvar6->SetName(Form("fig_lf_multistrange_%s", hvar6->GetName()));
       hvar6->SetName(Form("fig_lf_multistrange_6"));
       hvar6->Write();

       //Pad 2: Min V0 Radius Fid. Vol.  
       c2->cd(2);
       gPad->SetLogy();
       TH1D *hvar7 = cf->ShowProjection(7,icasType);
       hvar7->GetXaxis()->SetRangeUser(0.,3.0);
       hvar7->Draw("histo");

       hvar7->SetName(Form("fig_lf_multistrange_7"));
       hvar7->Write();

       Double_t x7;
       if      (collidingsystem == 0) x7 = 0.9;
       else if (collidingsystem == 1) x7 = 0.2;
       TLine *line7 = new TLine(x7,0.,x7,hvar7->GetBinContent(hvar7->GetMaximumBin()));
       line7->SetLineColor(kRed);
       line7->SetLineStyle(9);
       line7->SetLineWidth(2.0);
       line7->Draw("same");
          Bool_t check_7 = checkUnderTheLimit(hvar7, x7);
          if (check_7) { cout<<"The cut is OK!!"<<endl; t1->Draw(); }
          else         { cout<<"The cut is NOT OK!!"<<endl; t2->Draw(); }
       //Pad3: Min DCA V0 To PV
       c2->cd(3);
       gPad->SetLogy();
       TH1D *hvar8 = cf->ShowProjection(8,icasType);
       hvar8->GetXaxis()->SetRangeUser(0.,0.3);
       hvar8->Draw("histo");

       //      hvar8->SetName(Form("fig_lf_multistrange_%s", hvar8->GetName()));
       hvar8->SetName(Form("fig_lf_multistrange_8"));
       hvar8->Write();

       Double_t x8;
       if      (collidingsystem == 0) x8 = 0.05;
       else if (collidingsystem == 1) x8 = 0.01;
       TLine *line8 = new TLine(x8,0.,x8,hvar8->GetBinContent(hvar8->GetMaximumBin()));
       line8->SetLineColor(kRed);
       line8->SetLineStyle(9);
       line8->SetLineWidth(2.0);
       line8->Draw("same");
          Bool_t check_8 = checkUnderTheLimit(hvar8, x8);
          if (check_8) { cout<<"The cut is OK!!"<<endl; t1->Draw(); }
          else         { cout<<"The cut is NOT OK!!"<<endl; t2->Draw(); }
       //Pad 4: Min DCA Pos To PV
       c2->cd(4);
       gPad->SetLogy();
       TH1D *hvar9 = cf->ShowProjection(9,icasType);
       hvar9->GetXaxis()->SetRangeUser(0.,0.2);
       hvar9->Draw("histo");

       //   hvar9->SetName(Form("fig_lf_multistrange_%s", hvar9->GetName()));
       hvar9->SetName(Form("fig_lf_multistrange_9"));
       hvar9->Write();

       Double_t x9;
       if      (collidingsystem == 0) x9 = 0.1;
       else if (collidingsystem == 1) x9 = 0.05;
       TLine *line9 = new TLine(x9,0.,x9,hvar9->GetBinContent(hvar9->GetMaximumBin()));
       line9->SetLineColor(kRed);
       line9->SetLineStyle(9);
       line9->SetLineWidth(2.0);
       line9->Draw("same");

          Bool_t check_9 = checkUnderTheLimit(hvar9, x9);
          if (check_9) { cout<<"The cut is OK!!"<<endl; t1->Draw(); }
          else         { cout<<"The cut is NOT OK!!"<<endl; t2->Draw(); }
       //Pad 5: Min DCA Neg To PV
       c2->cd(5);
       gPad->SetLogy();
       TH1D *hvar10 = cf->ShowProjection(10,icasType);
       hvar10->GetXaxis()->SetRangeUser(0.,0.2);
       hvar10->Draw("histo");

       //    hvar10->SetName(Form("fig_lf_multistrange_%s", hvar10->GetName()));
       hvar10->SetName(Form("fig_lf_multistrange_10"));
       hvar10->Write();


       Double_t x10;
       if      (collidingsystem == 0) x10 = 0.1;
       else if (collidingsystem == 1) x10 = 0.05;
       TLine *line10 = new TLine(x10,0.,x10,hvar10->GetBinContent(hvar10->GetMaximumBin()));
       line10->SetLineColor(kRed);
       line10->SetLineStyle(9);
       line10->SetLineWidth(2.0);
       line10->Draw("same");
          Bool_t check_10 = checkUnderTheLimit(hvar10, x10);
          if (check_10) { cout<<"The cut is OK!!"<<endl; t1->Draw(); }
          else         { cout<<"The cut is NOT OK!!"<<endl; t2->Draw(); }
       //Pad 6: V0 cosine of Pointing Angle to XiV
       c2->cd(6);
       gPad->SetLogy();
       TH1D *hvar20 = cf->ShowProjection(18,icasType);
       Double_t max20 = hvar20->GetBinContent(hvar20->GetMaximumBin());
       hvar20->GetYaxis()->SetRangeUser(0.01,max20*1.5);
       hvar20->Draw("histo");

       //  hvar20->SetName(Form("fig_lf_multistrange_%s", hvar20->GetName()));
       hvar20->SetName(Form("fig_lf_multistrange_20"));
       hvar20->Write();


       c2->SaveAs(Form("fig_lf_MultistrangeQA_2.%s",suffix.Data()));
       

     //DEFINE 3st CANVAS AND DRAW PLOTS
     TCanvas *c3 = new TCanvas("c3","",1200,800);
     c3->Divide(2,3);
       //Pad 1: InvMass
       c3->cd(1);
       TH1D *hvar12 = cf->ShowProjection(11+icasType/2,icasType);
       hvar12->Draw("histo");

       hvar12->SetName(Form("fig_lf_multistrange_12", hvar12->GetName()));
       hvar12->Write();

       tpdgmass->Draw(); 
       TLine *linemass;
       if      (icasType == 0) linemass = new TLine(1.32171,0.,1.32171,0.5*hvar12->GetBinContent(hvar12->GetMaximumBin()));
       else if (icasType == 1) linemass = new TLine(1.32171,0.,1.32171,0.5*hvar12->GetBinContent(hvar12->GetMaximumBin()));
       else if (icasType == 2) linemass = new TLine(1.67245,0.,1.67245,0.5*hvar12->GetBinContent(hvar12->GetMaximumBin()));
       else if (icasType == 3) linemass = new TLine(1.67245,0.,1.67245,0.5*hvar12->GetBinContent(hvar12->GetMaximumBin()));
       linemass->SetLineColor(kRed);
       linemass->SetLineStyle(1);
       linemass->SetLineWidth(2.0);
       linemass->Draw("same");
       //Pad 2: Transverse momentum
       c3->cd(2);
       TH1D *hvar13 = cf->ShowProjection(13,icasType);
       hvar13->Draw("histo");

       //  hvar13->SetName(Form("fig_lf_multistrange_%s", hvar13->GetName()));
       hvar13->SetName(Form("fig_lf_multistrange_13"));
       hvar13->Write();
       //Pad 3: Y
       c3->cd(3);
       TH1D *hvar14 = cf->ShowProjection(14+icasType/2,icasType);
       hvar14->Draw("histo");

       //  hvar14->SetName(Form("fig_lf_multistrange_%s", hvar14->GetName()));
       hvar14->SetName(Form("fig_lf_multistrange_14"));
       hvar14->Write();

       //Pad 4: Cascade proper length
       c3->cd(4);
       TH1D *hvar18;
       hvar18 = cf->ShowProjection(16,icasType);
       hvar18->GetXaxis()->SetRangeUser(0.,90.);
       hvar18->Draw("histo");

       //   hvar18->SetName(Form("fig_lf_multistrange_%s", hvar18->GetName()));
       hvar18->SetName(Form("fig_lf_multistrange_18", hvar18->GetName()));
       hvar18->Write();

       //Pad 5: V0 proper length 
       c3->cd(5);
       TH1D *hvar19;
       hvar19 = cf->ShowProjection(17,icasType);
       hvar19->GetXaxis()->SetRangeUser(0.,90.);
       hvar19->Draw("histo");
       
       hvar19->SetName(Form("fig_lf_multistrange_19", hvar19->GetName()));
       hvar19->Write();
       

      //Pad 6
      // empty 
       if      (collidingsystem == 1) { c3->SaveAs(Form("fig_lf_MultistrangeQA_3.%s",suffix.Data())); }
       else if (collidingsystem == 0) { c3->SaveAs(Form("fig_lf_MultistrangeQA_3.%s",suffix.Data())); }

    
     //DEFINE 4st CANVAS AND DRAW PLOTS
    TCanvas *c4 = new TCanvas("c4","",600,400);
    c4->Divide(2,1);
      //Pad1: invariant mass fit
      c4->cd(1);
      TH1D *hvar18 = cf->ShowProjection(11+icasType/2,icasType);
      hvar18->Draw("histo");

      hvar18->SetName(Form("fig_lf_multistrange_18_1", hvar18->GetName()));
      hvar18->Write();

       // - SOME PARAMETER VALUE
       Bool_t kfitgauss = kFALSE;
       Bool_t kfitleft  = kFALSE;
       Bool_t kfitright = kFALSE;
       Int_t  ptbinNarrowY = 0;
       if (icasType < 2) ptbinNarrowY = 10;   // 6;
       else              ptbinNarrowY =  3;   // 2;
       // - SOME DEFINITIONS
       Float_t lowlimmass;
       Float_t uplimmass;
       Float_t lowgausslim;
       Float_t upgausslim;
       if (icasType==0||icasType==1) {
           lowlimmass=1.30;
           uplimmass=1.34;
           lowgausslim=1.312;
           upgausslim=1.332;
       } else {
           lowlimmass=1.645;
           uplimmass=1.70;
           lowgausslim=1.668;
           upgausslim=1.678;
       }
       TF1*  fitinvmass = new TF1("fitinvmass","gaus(0)+pol2(3)",lowlimmass,uplimmass);
       fitinvmass->SetParName(0, "cnstntG");
       fitinvmass->SetParName(1, "meanG");
       fitinvmass->SetParName(2, "sigmaG");
       fitinvmass->SetParLimits(0,0.,500000.);
       if (icasType==0||icasType==1) {
           fitinvmass->SetParameter(1, 1.32171);
           fitinvmass->SetParLimits(1, 1.31,1.33);
           fitinvmass->SetParLimits(2,0.001,0.005);
       } else {
           fitinvmass->SetParameter(1, 1.67245);
           fitinvmass->SetParLimits(1, 1.664,1.68);
           fitinvmass->SetParLimits(2,0.0008,0.006);
       }
       hvar18->Fit("fitinvmass","rimeN");
       fitinvmass->SetLineColor(kRed);
       fitinvmass->Draw("same");
       Float_t meanGauss   = fitinvmass->GetParameter(1);
       Float_t sigmaGauss  = fitinvmass->GetParameter(2);
       cout<<"Mean: "<<meanGauss<<endl;
       cout<<"Sigma: "<<sigmaGauss<<endl;
     //Pad2: Text
      c4->cd(2);
       Float_t refwidth = 0.002;
      TPaveText *pave1 = new TPaveText(0.05,0.3,0.95,0.5);
      pave1->SetFillColor(0);
      pave1->SetTextSize(0.04);
      pave1->SetTextAlign(12);
      if (icasType < 2) pave1->AddText("PDG mass: 1.32171 GeV/c^{2}");
      else              pave1->AddText("PDG mass: 1.67245 GeV/c^{2}");
      pave1->AddText(Form("#color[1]{Mass form Fit: %.5f #pm %.5f GeV/c^{2}}",meanGauss,sigmaGauss));
      if (sigmaGauss > refwidth - 0.0003 && sigmaGauss < refwidth + 0.0003) pave1->AddText("#color[3]{OK!! The width is compatible with standard.}");
      else                                                                  pave1->AddText("#color[2]{NOT OK!! Problem.}");
      pave1->Draw();
      cout<<"   "<<refwidth - 0.0003<<"<"<<sigmaGauss<<"<"<<refwidth + 0.0003<<endl;
    
     //DEFINE 5st CANVAS AND DRAW PLOTS
     if (collidingsystem == 0) {
       TCanvas *c5 = new TCanvas("c5","" );
         c5->Divide(2,1);
            //Pad 1: centrality
            c5->cd(1);
            TH1D *hvar16 = cf->ShowProjection(19,icasType);
            hvar16->Draw("histo");

	    hvar16->SetName(Form("fig_lf_multistrange_16", hvar16->GetName()));
	    hvar16->Write();

            //Pad 2: track multiplicity
            c5->cd(2);
            TH1D *hvar17 = cf->ShowProjection(20,icasType);
            hvar17->Draw("histo");
	    
	    hvar17->SetName(Form("fig_lf_multistrange_17", hvar17->GetName()));
	    hvar17->Write();
	    
	    
	    c5->SaveAs(Form("fig_lf_MultistrangeQA_5.%s",suffix.Data()));
     }
     
     fout->Close();
}




Bool_t checkUnderTheLimit(TH1D *lHist, Double_t limit) {

      Int_t binlimit = lHist->FindBin(limit);
      Bool_t checkOk = kTRUE;
      for (Int_t i = 1; i < binlimit; i++) {
           Int_t content = 0;
           content = lHist->GetBinContent(i);
           if (content != 0) checkOk = kFALSE;
           //cout<<"Content bin "<<i<<": "<<content<<endl;
      }
      return checkOk;

}


Bool_t checkOverTheLimit(TH1D *lHist, Double_t limit) {

      Int_t binlimit = lHist->FindBin(limit);
      Int_t lastbin  = lHist->GetNbinsX();
      Bool_t checkOk = kTRUE;
      for (Int_t i = binlimit; i < lastbin+1; i++) {
           Int_t content = 0;
           content = lHist->GetBinContent(i);
           if (content != 0) checkOk = kFALSE;
           //cout<<"Content bin "<<i<<": "<<content<<endl;
      }
      return checkOk;

}


