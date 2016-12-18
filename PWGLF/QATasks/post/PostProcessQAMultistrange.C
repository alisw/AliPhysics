//////////////////////////////////////////////////
//
//  This macro was written by Domenico Colella (domenico.colella@cern.ch).
//  12 November 2013
//   - Updated for Run2 (24 October 2016)
//
//   ------------------------
//   ------ Arguments -------
//   ------------------------
//   --  collidingsystem   =  0) PbPb  1) pp  2) pPb
//   --  isMC              =  kTRUE if running on MC production 
//   --  fileDir           =  "Input file directory"
//   --  filein            =  "Input file name"
//
//
//   -------------------------------------
//   ------ QATask output content --------
//   -------------------------------------
//   The output produced by the QATask is a CFContainer with 4 steps and 21 variables.
//   The meaning of each variable within the container are listed here:
//   --  0   = Max DCA Cascade Daughters                 pp: 2.0     PbPb: 0.3     pPb: 2.0
//   --  1   = Min DCA Bach To PV                        pp: 0.01    PbPb: 0.03    pPb: 0.03
//   --  2   = Min Cascade Cosine Of PA                  pp: 0.98    PbPb: 0.999   pPb: 0.95
//   --  3   = Min Cascade Radius Fid. Vol.              pp: 0.2     PbPb: 0.9     pPb: 0.4
//   --  4   = Window Invariant Mass Lambda              pp: 0.008   PbPb: 0.0008  pPb: 0.010
//   --  5   = Max DCA V0 Daughters                      pp: 1.5     PbPb: 1.0     pPb: 2.0
//   --  6   = Min V0 Cosine Of PA To PV                 pp: pT dep. PbPb: 0.98    pPb: 0.95
//   --  7   = Min V0 Radius Fid. Vol.                   pp: 0.2     PbPb: 0.9     pPb: 1.0
//   --  8   = Min DCA V0 To PV                          pp: 0.01    PbPb: 0.05    pPb: 0.05
//   --  9   = Min DCA Pos To PV                         pp: 0.05    PbPb: 0.1     pPb: 0.02
//   --  10  = Min DCA Neg To PV                         pp: 0.05    PbPb: 0.1     pPb: 0.02
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
//   In case of MC production one more CFContainer is produced, containing infos on the 
//   generated particles. As the previous container, this one is composed by 4 steps, one
//   for each cascade and 7 variables:
//    -- 0   = Total momentum
//    -- 1   = Transverse momentum
//    -- 2   = Rapidity
//    -- 3   = Pseudo-rapidity
//    -- 4   = Theta angle
//    -- 5   = Phi angle
//    -- 6   = Centrality
//    The previous container is still produced with the informations from the reconstructed
//    particles.
//
//
//   -----------------------------------
//   ------ Present Macro Checks -------
//   -----------------------------------
//   Using this macro many checks on the cascade topological reconstruction procedure
//   can be performed. In particular, the shape and the limit for the topological 
//   variable distributions as well as other kinematical variable distributions. The
//   reconstruction of the cascades are performed using two classes AliCascadeVertexer.cxx 
//   and AliV0vertexer.cxx contained in /STEER/ESD/ folder in Aliroot.
//   In the following are listed the contents of each page of the produced pdf:
//   
//   -- [Page 1] Distributions for the variables: 
//                DCA cascade daughters,  Bachelor IP to PV, 
//                Cascade cosine of PA,   Cascade radius of fiducial volume, 
//                Invariant mass Lambda,  DCA V0 daughters.
//   -- [Page 2] Distributions for the variables:
//                V0 cosine of PA to PV,  Min V0 Radius fiducial volume, 
//                Min DCA V0 To PV,       Min DCA positive To PV, 
//                Min DCA negative To PV, V0 cosine of PA to XiV
//   -- [Page 3] Distributions for the variables;
//                Transverse momentum,    Rapidity
//                Cascade proper length,  V0 proper length.
//   -- [Page 4] Check on the invariant mass distribution fit
//   -- [Page 5] Centrality percentile for PbPb and pPb and multiplicity percentile for pp.
//   -- [Page 6] Generated cascade multiplicity distribution
//   -- [Page 7] Only in case of MC production, distributions for the MC generated
//               particles, of the variables:
//                Total momentum,         Transverse momentum,
//                Rapidity,               Pseudo-rapidity,
//                Theta angle,            Phi angle,
//
//////////////////////////////////////////////////////




class AliCFContainer;



//=================
// - MAIN FUNCTION
//=================
void PostProcessQAMultistrange(Int_t   collidingsystem = 0,                             // 0) PbPb  1) pp 2) pPb
                               Bool_t  isMC            = kTRUE,                         // kTRUE-->MC and kFALSE-->Exp.
                               Char_t *fileDir         = ".",                           // Input file directory
                               Char_t *filein          = "AnalysisResults.root",        // Input file name
                               Char_t *output          = "pdf"                          // "eps", "png" or "pdf"
                              ) {




     //___________________
     //DEFINE DRAW OPTIONS
      myOptions();
      gROOT->ForceStyle();

     //_______________________
     //SOURCE USEFUL LIBRARIES
     gSystem->Load("libANALYSIS.so");
     gSystem->Load("libANALYSISalice.so");
     gSystem->Load("libCORRFW.so");

     //_________________________________
     //SOURCE THE FILE AND THE CONTAINER
     TFile *f = new TFile(Form("%s/%s",fileDir,filein));
     AliCFContainer *cf = (AliCFContainer*) (f->Get("PWGLFStrangeness.outputCheckCascade/fCFContCascadeCuts"));  

     //___________
     //DEFINE TEXT
     TLatex* t1 = new TLatex(0.6,0.7,"#color[3]{OK!!}");         myLatexMakeUp(t1,42,0.2,1);
     TLatex* t2 = new TLatex(0.6,0.7,"#color[2]{NOT OK!!}");     myLatexMakeUp(t2,42,0.2,2);
     Char_t *pname[4] = {"#color[1]{#Xi^{-}}", "#color[1]{#bar{#Xi}^{+}}", "#color[1]{#Omega^{-}}", "#color[1]{#bar{#Omega}^{+}}"};
     TLatex** tcas = new TLatex*[4];
     for (Int_t icas = 0; icas < 4; icas++) { tcas[icas] = new TLatex(0.2,0.7,Form("%s",pname[icas]));  myLatexMakeUp(tcas[icas],42,0.25,2); }
     const Int_t color[4] = {kRed+1,kOrange+1,kAzure+2,kViolet-4};
     Char_t *pdgmass[4] = {"PDG mass: 1.32171 GeV/c^{2}","PDG mass: 1.32171 GeV/c^{2}","PDG mass: 1.67245 GeV/c^{2}","PDG mass: 1.67245 GeV/c^{2}"}; 

     //_________________________________ 
     // - PAGE 1 : topological variables
     cout<<"--- Build the first page: topological variable distributions for reco. candidate ---"<<endl;
     TCanvas** c1 = new TCanvas*[4];
     TH1D** hvar0 = new TH1D*[4];  
     TH1D** hvar1 = new TH1D*[4];  
     TH1D** hvar2 = new TH1D*[4];  
     TH1D** hvar3 = new TH1D*[4];  
     TH1D** hvar4 = new TH1D*[4];  
     TH1D** hvar5 = new TH1D*[4];  
     Double_t x0 = 0, x1 = 0, x2 = 0, x3 = 0, x41 = 0, x42 = 0, x5 = 0;
     if      (collidingsystem == 0) { x0 = 0.3; x1 = 0.03; x2 = 0.999; x3 = 0.9; x41 = 1.116 + 0.008; x42 = 1.115 - 0.008; x5 = 1.0; }
     else if (collidingsystem == 1) { x0 = 2.0; x1 = 0.01; x2 = 0.980; x3 = 0.2; x41 = 1.116 + 0.008; x42 = 1.115 - 0.008; x5 = 1.5; }
     else if (collidingsystem == 2) { x0 = 2.0; x1 = 0.03; x2 = 0.950; x3 = 0.4; x41 = 1.116 + 0.010; x42 = 1.115 - 0.010; x5 = 2.0; }
     TLine** line0  = new TLine*[4];
     TLine** line1  = new TLine*[4];
     TLine** line2  = new TLine*[4];      
     TLine** line3  = new TLine*[4];
     TLine** line41 = new TLine*[4];
     TLine** line42 = new TLine*[4];
     TLine** line5  = new TLine*[4];
     for (Int_t icas = 0; icas < 4; icas++)  {
           c1[icas] = new TCanvas(Form("c1_%i",icas),"",1200,800);
           c1[icas]->Divide(2,3); 
           // -- Pad 1: DCA cascade daughters
           c1[icas]->cd(1);
            myPadSetUp(gPad); 
            hvar0[icas] = cf->ShowProjection(0,icas);  
            hvar0[icas]->SetLineColor(color[icas]);
            hvar0[icas]->Draw("histo");
            line0[icas] = new TLine(x0,0.,x0,hvar0[icas]->GetBinContent(hvar0[icas]->GetMaximumBin()));
            myLineMakeUp(line0[icas],2,9,2.0);
            line0[icas]->Draw("same");
            if (checkOverTheLimit(hvar0[icas],x0)) { cout<<"The cut is OK!!"<<endl;     t1->Draw(); }
            else                                   { cout<<"The cut is NOT OK!!"<<endl; t2->Draw(); }
            tcas[icas]->Draw();
           // -- Pad 2: Bachelor IP to PV
           c1[icas]->cd(2);
            myPadSetUp(gPad); 
            hvar1[icas] = cf->ShowProjection(1,icas);
            hvar1[icas]->SetLineColor(color[icas]);
            hvar1[icas]->GetXaxis()->SetRangeUser(0.,0.24);
            hvar1[icas]->Draw("histo");
            line1[icas] = new TLine(x1,0.,x1,hvar1[icas]->GetBinContent(hvar1[icas]->GetMaximumBin()));
            myLineMakeUp(line1[icas],2,9,2.0);
            line1[icas]->Draw("same");
            if (checkUnderTheLimit(hvar1[icas], x1)) { cout<<"The cut is OK!!"<<endl;     t1->Draw(); }
            else                                     { cout<<"The cut is NOT OK!!"<<endl; t2->Draw(); }
           // -- Pad 3: Cascade cosine of Pointing Angle
           c1[icas]->cd(3);
            myPadSetUp(gPad); 
            hvar2[icas] = cf->ShowProjection(2,icas);
            hvar2[icas]->SetLineColor(color[icas]);
            hvar2[icas]->GetYaxis()->SetRangeUser(0.01,(hvar2[icas]->GetMaximumBin())*1.5);
            hvar2[icas]->Draw("histo");
            line2[icas] = new TLine(x2,0.,x2,hvar2[icas]->GetBinContent(hvar2[icas]->GetMaximumBin()));
            myLineMakeUp(line2[icas],2,9,2.0);
            line2[icas]->Draw("same");
            if (checkUnderTheLimit(hvar2[icas], x2)) { cout<<"The cut is OK!!"<<endl;     t1->Draw(); }
            else                                     { cout<<"The cut is NOT OK!!"<<endl; t2->Draw(); }
           // -- Pad 4: Cascade radius of fiducial volume
           c1[icas]->cd(4);
            myPadSetUp(gPad); 
            hvar3[icas] = cf->ShowProjection(3,icas);
            hvar3[icas]->SetLineColor(color[icas]);
            hvar3[icas]->GetXaxis()->SetRangeUser(0.,3.8);
            hvar3[icas]->Draw("histo");
            line3[icas] = new TLine(x3,0.,x3,hvar3[icas]->GetBinContent(hvar3[icas]->GetMaximumBin()));
            myLineMakeUp(line3[icas],2,9,2.0);
            line3[icas]->Draw("same");
            if (checkUnderTheLimit(hvar3[icas], x3)) { cout<<"The cut is OK!!"<<endl;     t1->Draw(); }
            else                                     { cout<<"The cut is NOT OK!!"<<endl; t2->Draw(); }
           // -- Pad 5: Invariant mass Lambda
           c1[icas]->cd(5);
            myPadSetUp(gPad,kFALSE); 
            hvar4[icas] = cf->ShowProjection(4,icas);
            hvar4[icas]->SetLineColor(color[icas]);
            hvar4[icas]->Draw("histo");
            line41[icas] = new TLine(x41,0.,x41,hvar4[icas]->GetBinContent(hvar4[icas]->GetMaximumBin()));
            myLineMakeUp(line41[icas],2,9,2.0);
            line41[icas]->Draw("same");
            line42[icas] = new TLine(x42,0.,x42,hvar4[icas]->GetBinContent(hvar4[icas]->GetMaximumBin()));
            myLineMakeUp(line42[icas],2,9,2.0);
            line42[icas]->Draw("same");
            if (checkUnderTheLimit(hvar4[icas], x42) && checkOverTheLimit(hvar4[icas], x41)) { cout<<"The cut is OK!!"<<endl;     t1->Draw(); }
            else                                                                             { cout<<"The cut is NOT OK!!"<<endl; t2->Draw(); }
           // -- Pad 6: DCA V0 daughters
           c1[icas]->cd(6);
            myPadSetUp(gPad); 
            hvar5[icas] = cf->ShowProjection(5,icas);
            hvar5[icas]->SetLineColor(color[icas]);
            hvar5[icas]->Draw("histo");
            line5[icas] = new TLine(x5,0.,x5,hvar5[icas]->GetBinContent(hvar5[icas]->GetMaximumBin()));
            myLineMakeUp(line5[icas],2,9,2.0); 
            line5[icas]->Draw("same");
            if (checkOverTheLimit(hvar5[icas], x5)) { cout<<"The cut is OK!!"<<endl;     t1->Draw(); }
            else                                    { cout<<"The cut is NOT OK!!"<<endl; t2->Draw(); }

           if      (output == "png") c1[icas]->SaveAs(Form("fig_lf_Multistrange_page1_%i.png",icas));
           else if (output == "eps") c1[icas]->SaveAs(Form("fig_lf_Multistrange_page1_%i.eps",icas));
           else if (output == "pdf") c1[icas]->SaveAs("fig_lf_Multistrange.pdf(");
     }


     //_________________________________
     // - PAGE 2 : topological variables
     cout<<"--- Build the second page: topological variable distributions for reco. candidate ---"<<endl;
     TCanvas** c2 = new TCanvas*[4];
     TH1D** hvar6  = new TH1D*[4];
     TH1D** hvar7  = new TH1D*[4];
     TH1D** hvar8  = new TH1D*[4];
     TH1D** hvar9  = new TH1D*[4];
     TH1D** hvar10 = new TH1D*[4];
     TH1D** hvar11 = new TH1D*[4];
     Double_t x7 = 0, x8 = 0, x9 = 0, x10 = 0;
     if      (collidingsystem == 0) { x7 = 0.9; x8 = 0.05; x9 = 0.10; x10 = 0.10; }
     else if (collidingsystem == 1) { x7 = 0.2; x8 = 0.01; x9 = 0.05; x10 = 0.05; }
     else if (collidingsystem == 2) { x7 = 0.4; x8 = 0.05; x9 = 0.02; x10 = 0.02; }
     TLine** line7  = new TLine*[4];
     TLine** line8  = new TLine*[4];
     TLine** line9  = new TLine*[4];
     TLine** line10 = new TLine*[4];

     for (Int_t icas = 0; icas < 4; icas++)  {
           c2[icas] = new TCanvas(Form("c2_%i",icas),"",1200,800);
           c2[icas]->Divide(2,3);
           // -- Pad 1: V0 cosine of Pointing Angle to PV
           c2[icas]->cd(1);
            myPadSetUp(gPad);
            hvar6[icas] = cf->ShowProjection(6,icas);
            hvar6[icas]->SetLineColor(color[icas]);
            hvar6[icas]->GetYaxis()->SetRangeUser(0.01,(hvar6[icas]->GetBinContent(hvar6[icas]->GetMaximumBin()))*1.5);
            hvar6[icas]->Draw("histo");
            tcas[icas]->Draw();
           // -- Pad 2: Min V0 Radius Fid. Vol.  
           c2[icas]->cd(2);
            myPadSetUp(gPad); 
            hvar7[icas] = cf->ShowProjection(7,icas);
            hvar7[icas]->SetLineColor(color[icas]);
            hvar7[icas]->GetXaxis()->SetRangeUser(0.,3.0);
            hvar7[icas]->Draw("histo");
            line7[icas] = new TLine(x7,0.,x7,hvar7[icas]->GetBinContent(hvar7[icas]->GetMaximumBin()));
            myLineMakeUp(line7[icas],2,9,2.0);
            line7[icas]->Draw("same");
            if (checkUnderTheLimit(hvar7[icas], x7)) { cout<<"The cut is OK!!"<<endl;     t1->Draw(); }
            else                                     { cout<<"The cut is NOT OK!!"<<endl; t2->Draw(); }
           // -- Pad3: Min DCA V0 To PV
           c2[icas]->cd(3);
            myPadSetUp(gPad); 
            hvar8[icas] = cf->ShowProjection(8,icas);
            hvar8[icas]->SetLineColor(color[icas]);
            hvar8[icas]->GetXaxis()->SetRangeUser(0.,0.3);
            hvar8[icas]->Draw("histo");
            line8[icas] = new TLine(x8,0.,x8,hvar8[icas]->GetBinContent(hvar8[icas]->GetMaximumBin()));
            myLineMakeUp(line8[icas],2,9,2.0);
            line8[icas]->Draw("same");
            if (checkUnderTheLimit(hvar8[icas], x8)) { cout<<"The cut is OK!!"<<endl;     t1->Draw(); }
            else                                     { cout<<"The cut is NOT OK!!"<<endl; t2->Draw(); }
           // -- Pad 4: Min DCA Pos To PV
           c2[icas]->cd(4);
            myPadSetUp(gPad); 
            hvar9[icas] = cf->ShowProjection(9,icas);
            hvar9[icas]->SetLineColor(color[icas]);
            hvar9[icas]->GetXaxis()->SetRangeUser(0.,0.2);
            hvar9[icas]->Draw("histo");
            line9[icas] = new TLine(x9,0.,x9,hvar9[icas]->GetBinContent(hvar9[icas]->GetMaximumBin()));
            myLineMakeUp(line9[icas],2,9,2.0);
            line9[icas]->Draw("same");
            if (checkUnderTheLimit(hvar9[icas], x9)) { cout<<"The cut is OK!!"<<endl;     t1->Draw(); }
            else                                     { cout<<"The cut is NOT OK!!"<<endl; t2->Draw(); }
           // -- Pad 5: Min DCA Neg To PV
           c2[icas]->cd(5);
            myPadSetUp(gPad); 
            hvar10[icas] = cf->ShowProjection(10,icas);
            hvar10[icas]->SetLineColor(color[icas]);
            hvar10[icas]->GetXaxis()->SetRangeUser(0.,0.2);
            hvar10[icas]->Draw("histo");
            line10[icas] = new TLine(x10,0.,x10,hvar10[icas]->GetBinContent(hvar10[icas]->GetMaximumBin()));
            myLineMakeUp(line10[icas],2,9,2.0);
            line10[icas]->Draw("same");
            if (checkUnderTheLimit(hvar10[icas], x10)) { cout<<"The cut is OK!!"<<endl;     t1->Draw(); }
            else                                       { cout<<"The cut is NOT OK!!"<<endl; t2->Draw(); }
           // -- Pad 6: V0 cosine of Pointing Angle to Xi vtx
           c2[icas]->cd(6);
            myPadSetUp(gPad); 
            hvar11[icas] = cf->ShowProjection(18,icas);
            hvar11[icas]->SetLineColor(color[icas]);
            hvar11[icas]->GetYaxis()->SetRangeUser(0.01,(hvar11[icas]->GetBinContent(hvar11[icas]->GetMaximumBin()))*1.5);
            hvar11[icas]->Draw("histo");
     
           if      (output == "png") c2[icas]->SaveAs(Form("fig_lf_Multistrange_page2_%i.png",icas));
           else if (output == "eps") c2[icas]->SaveAs(Form("fig_lf_Multistrange_page2_%i.eps",icas));
           else if (output == "pdf") c2[icas]->SaveAs("fig_lf_Multistrange.pdf");
     }


     //_____________________________
     // - PAGE 3 : further variables
     cout<<"--- Build the third page: further variables for reco. candidate ---"<<endl;
     TCanvas** c3 = new TCanvas*[4];
     TH1D** hvar12 = new TH1D*[4];
     TH1D** hvar13 = new TH1D*[4];
     TH1D** hvar14 = new TH1D*[4];
     TH1D** hvar15 = new TH1D*[4];

     for (Int_t icas = 0; icas < 4; icas++)  {
           c3[icas] = new TCanvas(Form("c3_%i",icas),"",1200,800);
           c3[icas]->Divide(2,3);
           // -- Pad 1: Transverse momentum
           c3[icas]->cd(1);
            myPadSetUp(gPad,kFALSE);
            hvar12[icas] = cf->ShowProjection(13,icas);
            hvar12[icas]->SetLineColor(color[icas]);
            hvar12[icas]->Draw("histo");
            tcas[icas]->Draw();
           // -- Pad 2: Y
           c3[icas]->cd(2);
            myPadSetUp(gPad,kFALSE);
            hvar13[icas] = cf->ShowProjection(14+icas/2,icas);
            hvar13[icas]->SetLineColor(color[icas]);
            hvar13[icas]->Draw("histo");
           // -- Pad 3: Cascade proper length
           c3[icas]->cd(3);
            myPadSetUp(gPad,kFALSE);
            hvar14[icas] = cf->ShowProjection(16,icas);
            hvar14[icas]->GetXaxis()->SetRangeUser(0.,90.);
            hvar14[icas]->SetLineColor(color[icas]);
            hvar14[icas]->Draw("histo");
           // -- Pad 4: V0 proper length 
           c3[icas]->cd(4);
            myPadSetUp(gPad,kFALSE);
            hvar15[icas] = cf->ShowProjection(17,icas);
            hvar15[icas]->SetLineColor(color[icas]);
            hvar15[icas]->GetXaxis()->SetRangeUser(0.,90.);
            hvar15[icas]->Draw("histo");
           // -- Pad 5 & 6
           // empty 
           if      (output == "png") c3[icas]->SaveAs(Form("fig_lf_Multistrange_page3_%i.png",icas));
           else if (output == "eps") c3[icas]->SaveAs(Form("fig_lf_Multistrange_page3_%i.eps",icas));
           else if (output == "pdf") c3[icas]->SaveAs("fig_lf_Multistrange.pdf");
     }


     //____________________________________ 
     // - PAGE 4 : invariant mass fit check
     cout<<"--- Build the fourth page: invariant mass fit for reco. candidate ---"<<endl;
     TCanvas *c4 = new TCanvas("c4","",1200,800);
     c4->Divide(2,2);
     Float_t lowlimmass[4] = {1.300, 1.300, 1.645, 1.645};
     Float_t uplimmass[4]  = {1.340, 1.340, 1.700, 1.700};
     Double_t param1[4]    = {1.32171, 1.32171, 1.67245, 1.67245};
     Double_t param2min[4] = {1.310, 1.310, 1.664, 1.664};
     Double_t param2max[4] = {1.330, 1.330, 1.680, 1.680};
     Double_t param3min[4] = {0.001, 0.001, 0.0008, 0.0008};
     Double_t param3max[4] = {0.005, 0.005, 0.0060, 0.0060};     
     Float_t refwidth[4]   = {0.002, 0.002, 0.0025, 0.0025};
     TH1D** hvar16 = new TH1D*[4];
     TF1** fitinvmass = new TF1*[4];
     TPaveText** pave = new TPaveText[4];
     Double_t meanGauss[4], sigmaGauss[4];

     for (Int_t icas = 0; icas < 4; icas++)  {
           c4->cd(icas+1); 
            myPadSetUp(gPad,kFALSE);
            hvar16[icas] = cf->ShowProjection(11+icas/2,icas);
            hvar16[icas]->Draw("histo");
            fitinvmass[icas] = new TF1(Form("fitinvmass[%i]",icas),"gaus(0)+pol2(3)",lowlimmass[icas],uplimmass[icas]);
            fitinvmass[icas]->SetParName(0, "cnstntG");
            fitinvmass[icas]->SetParName(1, "meanG");
            fitinvmass[icas]->SetParName(2, "sigmaG");
            fitinvmass[icas]->SetParLimits(0,0.,500000.);
            fitinvmass[icas]->SetParameter(1, param1[icas]);
            fitinvmass[icas]->SetParLimits(1, param2min[icas],param2max[icas]);
            fitinvmass[icas]->SetParLimits(2, param3min[icas],param3max[icas]);
            hvar16[icas]->Fit(Form("fitinvmass[%i]",icas),"rimeN");
            fitinvmass[icas]->SetLineColor(kRed);
            fitinvmass[icas]->Draw("same");
            meanGauss[icas]  = fitinvmass[icas]->GetParameter(1);
            sigmaGauss[icas] = fitinvmass[icas]->GetParameter(2);
            pave[icas] = new TPaveText(0.55,0.75,0.95,0.95,"NDC");
            myPaveMakeUp(pave[icas],Form("%s",pdgmass[icas]),0,22,0.04,0,0);
            pave[icas]->AddText(Form("#color[1]{Mass form Fit: %.5f #pm %.5f GeV/c^{2}}",meanGauss[icas],sigmaGauss[icas]));
            if (sigmaGauss[icas] > refwidth[icas] - 0.001 && sigmaGauss[icas] < refwidth[icas] + 0.001) pave[icas]->AddText("#color[3]{OK!! The width is compatible with standard.}");
            else                                                                                        pave[icas]->AddText("#color[2]{NOT OK!! Problem.}"); 
            pave[icas]->Draw("same");
            tcas[icas]->Draw();
     }
     if      (output == "png") c4->SaveAs("fig_lf_Multistrange_page4.png");
     else if (output == "eps") c4->SaveAs("fig_lf_Multistrange_page4.eps");
     else if (output == "pdf") c4->SaveAs("fig_lf_Multistrange.pdf");   


     //________________________________________________ 
     // - PAGE 5 : centrality/multiplicity distribution 
     cout<<"--- Build the fifth page: invariant mass fit for reco. candidate ---"<<endl;
     TCanvas *c5 = new TCanvas("c5","",600,400);
     myPadSetUp(gPad,kFALSE,0.11,0.06,0.04,0.15);
     TH1D *hvar17 = cf->ShowProjection(19,0);
     hvar17->Draw("histo");
     if      (output == "png")          c5->SaveAs("fig_lf_Multistrange_page5.png");
     else if (output == "eps")          c5->SaveAs("fig_lf_Multistrange_page5.eps");
     else if (output == "pdf" && isMC)  c5->SaveAs("fig_lf_Multistrange.pdf");
     else if (output == "pdf" && !isMc) c5->SaveAs("fig_lf_Multistrange.pdf)");


     //_____________________________________________
     // - PAGE 6 : cascade multiplicity distribution
     cout<<"--- Build the sixth page: generated cascade multiplicity distribution ---"<<endl;
     if (isMC) {
        TList *l = (TList*) f->Get("PWGLFStrangeness.outputCheckCascade/fListHistMultistrangeQA");
        Char_t *histoname[4] = {"fHistCascadeMultiplicityXiMinus","fHistCascadeMultiplicityXiPlus","fHistCascadeMultiplicityOmegaMinus","fHistCascadeMultiplicityOmegaPlus"}; 
        TCanvas *c6 = new TCanvas("c6","",1200,800);
        c6->Divide(2,2);
        TH1D** hvar18 = new TH1D*[4];
        Double_t mean[4] = {};
        Double_t integral[4] = {};
        Int_t    entries[4] = {};
        TPaveText** pave1 = new TPaveText[4];

        for (Int_t icas = 0; icas < 4; icas++)  {
              c6->cd(icas+1);
               myPadSetUp(gPad,kFALSE); 
               hvar18[icas] = (TH1D*) l->FindObject(Form("%s",histoname[icas]));
               hvar18[icas]->Draw();
               for (Int_t j = 0; j < hvar18[icas]->GetNbinsX(); j++) {
                     mean[icas]     = mean[icas] + hvar18[icas]->GetBinCenter(j) * hvar18[icas]->GetBinWidth(j) * hvar18[icas]->GetBinContent(j);
                     integral[icas] = integral[icas] + hvar18[icas]->GetBinContent(j) * hvar18[icas]->GetBinWidth(j);
                     entries[icas] = entries[icas] + hvar18[icas]->GetBinContent(j);
               }
               mean[icas] = mean[icas] / integral[icas];
               tcas[icas]->Draw();  
               pave1[icas] = new TPaveText(0.55,0.75,0.95,0.95,"NDC");
               myPaveMakeUp(pave1[icas],Form("Mean : %.2f",mean[icas]),0,22,0.1,0,0);
               pave1[icas]->Draw("same");
        }
        if      (output == "png") c6->SaveAs("fig_lf_Multistrange_page5.png");
        else if (output == "eps") c6->SaveAs("fig_lf_Multistrange_page5.eps");
        else if (output == "pdf") c6->SaveAs("fig_lf_Multistrange.pdf");
     }


     //_______________________________________________
     // - PAGE 7 : MC generated particles check
     cout<<"--- Build the seventh page: general variables for MC generated particles ---"<<endl;
     if (isMC) { 
          AliCFContainer *cfMC = (AliCFContainer*) (f->Get("PWGLFStrangeness.outputCheckCascade/fCFContCascadeMCgen"));
          TCanvas** c7 = new TCanvas*[4];
          TH1D** hvar19 = new TH1D*[4];
          TH1D** hvar20 = new TH1D*[4];
          TH1D** hvar21 = new TH1D*[4]; 
          TH1D** hvar22 = new TH1D*[4];
          TH1D** hvar23 = new TH1D*[4];
          TH1D** hvar24 = new TH1D*[4];

          for (Int_t icas = 0; icas < 4; icas++)  {
                c7[icas] = new TCanvas(Form("c7_%i",icas),"",1200,800);
                c7[icas]->Divide(2,3);
                // -- Pad 1: Total Momentum
                c7[icas]->cd(1);
                 myPadSetUp(gPad,kFALSE);
                 hvar19[icas] = cfMC->ShowProjection(0,icas);
                 hvar19[icas]->SetLineColor(color[icas]);
                 hvar19[icas]->Draw("histo");
                 tcas[icas]->Draw();
                // -- Pad 2: Transverse Momentum
                c7[icas]->cd(2);
                 myPadSetUp(gPad,kFALSE);
                 hvar20[icas] = cfMC->ShowProjection(1,icas);
                 hvar20[icas]->SetLineColor(color[icas]);
                 hvar20[icas]->Draw("histo");
                // -- Pad 3: Rapidity (y)
                c7[icas]->cd(3);
                 myPadSetUp(gPad,kFALSE);
                 hvar21[icas] = cfMC->ShowProjection(2,icas);
                 hvar21[icas]->SetLineColor(color[icas]);
                 hvar21[icas]->Draw("histo");
                // -- Pad 4: Pseudo-rapidity (eta)
                c7[icas]->cd(4);
                 myPadSetUp(gPad,kFALSE);
                 hvar22[icas] = cfMC->ShowProjection(3,icas);
                 hvar22[icas]->SetLineColor(color[icas]);
                 hvar22[icas]->Draw("histo");
                // -- Pad 5: Theta
                c7[icas]->cd(5);
                 myPadSetUp(gPad,kFALSE);
                 hvar23[icas] = cfMC->ShowProjection(4,icas);
                 hvar23[icas]->SetLineColor(color[icas]);
                 hvar23[icas]->Draw("histo");
                // -- Pad 6: Phi
                c7[icas]->cd(6);
                 myPadSetUp(gPad,kFALSE);
                 hvar24[icas] = cfMC->ShowProjection(5,icas);
                 hvar24[icas]->SetLineColor(color[icas]);
                 hvar24[icas]->Draw("histo");
		if      (output == "png") c7[icas]->SaveAs(Form("fig_lf_Multistrange_page6_%i.png",icas));
                else if (output == "eps") c7[icas]->SaveAs(Form("fig_lf_Multistrange_page6_%i.eps",icas));
                else if (output == "pdf") {
                    if (icas < 3) c7[icas]->SaveAs("fig_lf_Multistrange.pdf");
                    else          c7[icas]->SaveAs("fig_lf_Multistrange.pdf)");   
                }
          }
     }

}




//====================
// - USEFUL FUNCTIONS
//====================
//______________________
Bool_t checkUnderTheLimit(TH1D *lHist, Double_t limit) {
         Int_t binlimit = lHist->FindBin(limit);
         Bool_t checkOk = kTRUE;
         for (Int_t i = 1; i < binlimit; i++) {
              Int_t content = 0;
              content = lHist->GetBinContent(i);
              if (content != 0) checkOk = kFALSE;
         }
         return checkOk;
}
//______________________
Bool_t checkOverTheLimit(TH1D *lHist, Double_t limit) {
         Int_t binlimit = lHist->FindBin(limit);
         Int_t lastbin  = lHist->GetNbinsX();
         Bool_t checkOk = kTRUE;
         for (Int_t i = binlimit; i < lastbin+1; i++) {
              Int_t content = 0;
              content = lHist->GetBinContent(i);
              if (content != 0) checkOk = kFALSE;
         }
         return checkOk;
}
//-----------------------------
void myOptions(Int_t lStat=0) {
        cout << "Set my personal style!" << endl;
        gStyle->SetOptStat("ie");
        gStyle->SetOptStat(kFALSE);
        gStyle->SetOptTitle(kFALSE);
        gStyle->SetOptLogy(kFALSE);
        gStyle->SetFrameLineWidth(2.5);
        gStyle->SetHistLineWidth(2.5);
        // Set gStyle
        int font = 42;
        // From plain
        gStyle->SetFrameBorderMode(0);
        gStyle->SetFrameFillColor(0);
        gStyle->SetCanvasBorderMode(0);
        gStyle->SetPadBorderMode(0);
        gStyle->SetPadColor(0);
        gStyle->SetCanvasColor(0);
        gStyle->SetTitleFillColor(0);
        gStyle->SetTitleBorderSize(1);
        gStyle->SetStatColor(0);
        gStyle->SetStatBorderSize(1);
        gStyle->SetLegendBorderSize(1);
        gStyle->SetDrawBorder(0);
        gStyle->SetTextFont(font);
        gStyle->SetStatFont(font);
        gStyle->SetStatFontSize(0.05);
        gStyle->SetStatX(0.97);
        gStyle->SetStatY(0.98);
        gStyle->SetStatH(0.03);
        gStyle->SetStatW(0.3);
        gStyle->SetTickLength(0.02,"y");
        gStyle->SetEndErrorSize(3);
        gStyle->SetLabelSize(0.05,"xyz");
        gStyle->SetLabelFont(font,"xyz");
        gStyle->SetLabelOffset(0.01,"xyz");
        gStyle->SetTitleFont(font,"xyz");
        gStyle->SetTitleOffset(1.0,"xyz");
        gStyle->SetTitleSize(0.06,"xyz");
        gStyle->SetMarkerSize(1);
        gStyle->SetGridColor(kGray);
        gStyle->SetPalette(1,0);
        if (lStat){
              gStyle->SetOptTitle(1);
              gStyle->SetOptStat(1111);
              gStyle->SetOptFit(1111);
        } else {
              gStyle->SetOptTitle(0);
              gStyle->SetOptStat(0);
              gStyle->SetOptFit(0);
        }
}
//----------------------------------------------------------
void myLatexMakeUp(TLatex *currentLatex, Int_t textfont, Double_t textsize, Int_t textcolor){
       currentLatex->SetNDC();
       currentLatex->SetTextFont(textfont);
       currentLatex->SetTextSize(textsize);
       currentLatex->SetTextColor(textcolor);
       return;
}
//----------------------------------------------------------
void myLineMakeUp(TLine *currentLine, Int_t linecolor, Int_t linestyle, Double_t linewidth){
       currentLine->SetLineColor(linecolor);
       currentLine->SetLineStyle(linestyle);
       currentLine->SetLineWidth(linewidth);
       return;
}
//-----------------------------------------------------------
void myPadSetUp(TPad *currentPad, Bool_t setlogy=kTRUE, Float_t currentLeft=0.11, Float_t currentTop=0.04, Float_t currentRight=0.04, Float_t currentBottom=0.15){
        if (setlogy) currentPad->SetLogy();
        currentPad->SetLeftMargin(currentLeft);
        currentPad->SetTopMargin(currentTop);
        currentPad->SetRightMargin(currentRight);
        currentPad->SetBottomMargin(currentBottom);
        return;
}
//----------------------------------------------------------
void myPaveMakeUp(TPaveText *currentPave, Char_t *texture, Double_t textangle, Int_t textalign, Double_t textsize, Int_t pavefillcolor, Double_t pavebordersize){
       TText *text = currentPave->AddText(Form("%s",texture));
       text->SetTextAngle(textangle);
       text->SetTextAlign(textalign);
       text->SetTextSize(textsize);
       currentPave->SetFillColor(pavefillcolor);
       currentPave->SetFillStyle(0);
       currentPave->SetBorderSize(pavebordersize);
       return;
}
