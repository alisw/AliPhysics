// $Id$
/*
 * Plotting macro for comparing offline- and HLT- ESD trees from
 * HLT-OFFLINE-GLOBAL-comparison.root produced using $ALICE_ROOT/HLT/QA/tasks/AliAnalysisTaskHLT.*
 *
 * Usage: aliroot drawGlobalESDHistograms.C'("HLT-OFFLINE-GLOBAL-comparison.root")'
 *
 * It saves the canvas with the output histograms in a png file.
 *
 * @ingroup alihlt_qa
 * @author Camilla.Stokkevag@student.uib.no, Kalliopi.Kanaki@ift.uib.no 
 */

void drawGlobalESDHistograms(const char* filename="HLT-OFFLINE-GLOBAL-comparison.root"){
 
 gROOT->SetStyle("Plain");
 gStyle->SetPalette(1);
 gStyle->SetOptStat(10);
 gStyle->SetTitleX(gStyle->GetPadLeftMargin());

 TFile *f1 = TFile::Open(filename); 
 if(!f1 || f1->IsZombie()) {
    printf("file %s does not exist or there is an error opening it\n", filename);
    return;
 }

 TList *l1 = (TList*)f1->Get("global_histograms");
 if(!l1){
    printf("No list %s contained in your input file\n", l1->GetName()); 
    return; 
 }
 

//=================================================================//
//--------------- READ HISTOGRAMS ---------------------------------//
//=================================================================//

 TH1F *h1 = l1->FindObject("fNcluster_hlt");
 TH1F *h2 = l1->FindObject("fNcluster_off"); 
 h1->SetTitle("cluster distribution");
 h1->GetXaxis()->SetTitle("TPC clusters per track");
 h2->SetLineColor(2);
 if(h1->GetMaximum() >= h2->GetMaximum())  h2->SetMaximum(1.1*h1->GetMaximum());
 else h1->SetMaximum(1.1*h2->GetMaximum());
  


 TH1F *h3  = l1->FindObject("fDCA_hlt");
 TH1F *h4  = l1->FindObject("fDCA_off"); 
 TH1F *hSG = l1->FindObject("fDCA_hltSG");
 h3->SetTitle("DCA between track and vertex on XY plane");
 h3->GetXaxis()->SetTitle("DCAr (cm)");
 h4->SetLineColor(2);
 hSG->SetLineColor(kBlue);
 if(h3->GetMaximum() >= h4->GetMaximum()) h4->SetMaximum(1.1*h3->GetMaximum());
 else h3->SetMaximum(1.1*h4->GetMaximum());
 

 TH1F *h5 = l1->FindObject("fMult_hlt");
 TH1F *h6 = l1->FindObject("fMult_off"); 
 h5->SetTitle("track multiplicity");
 h6->SetLineColor(2);
 if(h5->GetMaximum() > h6->GetMaximum()) h6->SetMaximum(1.1*h5->GetMaximum());
 else h5->SetMaximum(1.1*h6->GetMaximum());


 TH1F *h7 = l1->FindObject("fCharge_hlt");
 TH1F *h8 = l1->FindObject("fCharge_off");
 h7->GetXaxis()->SetTitle("polarity"); 
 h7->SetTitle("charge distribution");
 h8->SetLineColor(2);
 if(h7->GetMaximum() > h8->GetMaximum()) h8->SetMaximum(1.1*h7->GetMaximum());
 else h7->SetMaximum(1.1*h8->GetMaximum());
 

 TH1F *h9  = l1->FindObject("fMomentum_hlt");
 TH1F *h10 = l1->FindObject("fMomentum_off"); 
 h9->GetXaxis()->SetTitle("p_{t} (GeV/c)"); 
 h9->SetTitle("transverse momentum");
 h10->SetLineColor(2);
 if(h9->GetMaximum() > h10->GetMaximum()) h10->SetMaximum(1.1*h9->GetMaximum());
 else h9->SetMaximum(1.1*h10->GetMaximum());


 TH1F *h11 = l1->FindObject("fEta_hlt");
 TH1F *h12 = l1->FindObject("fEta_off"); 
 h11->SetTitle("pseudorapidity");
 h11->GetXaxis()->SetTitle("#eta");
 h12->SetLineColor(2);
 if(h11->GetMaximum() > h12->GetMaximum()) h12->SetMaximum(1.1*h11->GetMaximum());
 else h11->SetMaximum(1.1*h12->GetMaximum());


 TH1F *h13 = l1->FindObject("fXvertex_hlt");
 TH1F *h14 = l1->FindObject("fXvertex_off"); 
 h13->GetXaxis()->SetTitle("x (cm)");
 h13->SetTitle("x of primary vertex");
 h14->SetLineColor(2);
 if(h13->GetMaximum() > h14->GetMaximum()) h14->SetMaximum(1.1*h13->GetMaximum());
 else h13->SetMaximum(1.1*h14->GetMaximum());

 TH1F *h15 = l1->FindObject("fYvertex_hlt");
 TH1F *h16 = l1->FindObject("fYvertex_off"); 
 h15->GetXaxis()->SetTitle("y (cm)");
 h15->SetTitle("y of primary vertex");
 h16->SetLineColor(2);
 if(h15->GetMaximum() > h16->GetMaximum()) h16->SetMaximum(1.1*h15->GetMaximum());
 else h15->SetMaximum(1.1*h16->GetMaximum());

 TH1F *h17 = l1->FindObject("fZvertex_hlt");
 TH1F *h18 = l1->FindObject("fZvertex_off"); 
 h17->GetXaxis()->SetTitle("z (cm)");
 h17->SetTitle("z of primary vertex");
 h18->SetLineColor(2);
 if(h17->GetMaximum() > h18->GetMaximum()) h18->SetMaximum(1.1*h17->GetMaximum());
 else h17->SetMaximum(1.1*h18->GetMaximum());


//  TH2F *h15 = l1->FindObject("fXYvertex_off");
//  h15->GetXaxis()->SetTitle("X (cm)");
//  h15->GetYaxis()->SetTitle("Y (cm)");
//  h15->SetTitle("XY primary vertex Offline");
// 
//  TH2F *h16 = l1->FindObject("fXYvertex_hlt"); 
//  h16->GetXaxis()->SetTitle("X (cm)");
//  h16->GetYaxis()->SetTitle("Y (cm)");
//  h16->SetTitle("XY of primary vertex HLT");

 TLegend *leg1 = new TLegend(0.6,0.6,0.8,0.8);
 leg1->SetFillColor(10);
 leg1->SetLineColor(10);
 leg1->AddEntry(h1,"HLT", "l");
 leg1->AddEntry(h2,"OFF", "l");



//=================================================================//
//--------------------- DRAW HISTOGRAMS ---------------------------//
//=================================================================//



 TCanvas *c1 = new TCanvas("c1","HLT vs offline",1300,800);
 c1->Divide(3,3);

 c1->cd(1);
 h1->Draw();
 h2->Draw("sames");
 leg1->Draw("same");
 
 gPad->Update();
 TPaveStats *st1 = (TPaveStats*)h1->FindObject("stats");
 st1->SetLineColor(0);

 gPad->Update();
 TPaveStats *st2 = (TPaveStats*)h2->FindObject("stats");
 st2->SetY1NDC(st1->GetY1NDC()-0.05);
 st2->SetY2NDC(st1->GetY2NDC()-0.05);
 st2->SetLineColor(0);
 st2->SetTextColor(h2->GetLineColor());
 st2->SetFillStyle(0);
 st2->Draw();

 //======================

 c1->cd(2)->SetLogy();
 h3->Draw();
 h4->Draw("sames");
 hSG->Draw("sames");
 leg1->Draw("same");

 gPad->Update();
 TPaveStats *st3 = (TPaveStats*)h3->FindObject("stats");
 st3->SetLineColor(0);

 gPad->Update();
 TPaveStats *st4 = (TPaveStats*)h4->FindObject("stats");
 st4->SetY1NDC(st3->GetY1NDC()-0.05);
 st4->SetY2NDC(st3->GetY2NDC()-0.05);
 st4->SetLineColor(0);
 st4->SetTextColor(h4->GetLineColor());
 st4->SetFillStyle(0);
 st4->Draw();

 gPad->Update();
 TPaveStats *stSG = (TPaveStats*)hSG->FindObject("stats");
 stSG->SetY1NDC(st4->GetY1NDC()-0.05);
 stSG->SetY2NDC(st4->GetY2NDC()-0.05);
 stSG->SetLineColor(0);
 stSG->SetTextColor(hSG->GetLineColor());
 stSG->SetFillStyle(0);
 stSG->Draw();
 
//======================

 c1->cd(3);
 h5->Draw();
 h6->Draw("sames");
 leg1->Draw("same");
 
 gPad->Update();
 TPaveStats *st5 = (TPaveStats*)h5->FindObject("stats");
 st5->SetLineColor(0);

 gPad->Update();
 TPaveStats *st6 = (TPaveStats*)h6->FindObject("stats");
 st6->SetY1NDC(st5->GetY1NDC()-0.05);
 st6->SetY2NDC(st5->GetY2NDC()-0.05);
 st6->SetLineColor(0);
 st6->SetTextColor(h6->GetLineColor());
 st6->SetFillStyle(0);
 st6->Draw();
 
//======================

 c1->cd(4);
 h7->Draw();
 h8->Draw("sames");
 leg1->Draw("same");
 
 gPad->Update();
 TPaveStats *st7 = (TPaveStats*)h7->FindObject("stats");
 st7->SetLineColor(0);

 gPad->Update();
 TPaveStats *st8 = (TPaveStats*)h8->FindObject("stats");
 st8->SetY1NDC(st7->GetY1NDC()-0.05);
 st8->SetY2NDC(st7->GetY2NDC()-0.05);
 st8->SetLineColor(0);
 st8->SetTextColor(h8->GetLineColor());
 st8->SetFillStyle(0);
 st8->Draw();
 
//======================
 
 c1->cd(5)->SetLogy();
 h9->Draw();
 h10->Draw("sames");
 leg1->Draw("same");
 
 gPad->Update();
 TPaveStats *st9 = (TPaveStats*)h9->FindObject("stats");
 st9->SetLineColor(0);

 gPad->Update();
 TPaveStats *st10 = (TPaveStats*)h10->FindObject("stats");
 st10->SetY1NDC(st9->GetY1NDC()-0.05);
 st10->SetY2NDC(st9->GetY2NDC()-0.05);
 st10->SetLineColor(0);
 st10->SetTextColor(h10->GetLineColor());
 st10->SetFillStyle(0);
 st10->Draw();
 
//======================
 
 c1->cd(6);
 h11->Draw();
 h12->Draw("sames");
 leg1->Draw("same");
 
 gPad->Update();
 TPaveStats *st11 = (TPaveStats*)h11->FindObject("stats");
 st11->SetLineColor(0);

 gPad->Update();
 TPaveStats *st12 = (TPaveStats*)h12->FindObject("stats");
 st12->SetY1NDC(st11->GetY1NDC()-0.05);
 st12->SetY2NDC(st11->GetY2NDC()-0.05);
 st12->SetLineColor(0);
 st12->SetTextColor(h12->GetLineColor());
 st12->SetFillStyle(0);
 st12->Draw();
 
//======================

 c1->cd(7);
 h13->Draw();
 h14->Draw("sames");
 leg1->Draw("same");
 
 gPad->Update();
 TPaveStats *st13 = (TPaveStats*)h13->FindObject("stats");
 st13->SetLineColor(0);

 gPad->Update();
 TPaveStats *st14 = (TPaveStats*)h14->FindObject("stats");
 st14->SetY1NDC(st13->GetY1NDC()-0.05);
 st14->SetY2NDC(st13->GetY2NDC()-0.05);
 st14->SetLineColor(0);
 st14->SetTextColor(h14->GetLineColor());
 st14->SetFillStyle(0);
 st14->Draw();

//======================

 c1->cd(8);
 h15->Draw();
 h16->Draw("sames");
 leg1->Draw("same");
 
 gPad->Update();
 TPaveStats *st15 = (TPaveStats*)h15->FindObject("stats");
 st15->SetLineColor(0);

 gPad->Update();
 TPaveStats *st16 = (TPaveStats*)h16->FindObject("stats");
 st16->SetY1NDC(st15->GetY1NDC()-0.05);
 st16->SetY2NDC(st15->GetY2NDC()-0.05);
 st16->SetLineColor(0);
 st16->SetTextColor(h16->GetLineColor());
 st16->SetFillStyle(0);
 st16->Draw();

//======================

 c1->cd(9);
 h17->Draw();
 h18->Draw("sames");
 leg1->Draw("same");
 
 gPad->Update();
 TPaveStats *st17 = (TPaveStats*)h17->FindObject("stats");
 st17->SetLineColor(0);

 gPad->Update();
 TPaveStats *st18 = (TPaveStats*)h18->FindObject("stats");
 st18->SetY1NDC(st17->GetY1NDC()-0.05);
 st18->SetY2NDC(st17->GetY2NDC()-0.05);
 st18->SetLineColor(0);
 st18->SetTextColor(h18->GetLineColor());
 st18->SetFillStyle(0);
 st18->Draw();

//======================

 c1->SaveAs("HLT-offline.png");  
}
