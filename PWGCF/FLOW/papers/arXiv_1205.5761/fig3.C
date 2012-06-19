{
//=========Macro generated from canvas: cInt/cInt
//=========  (Wed Jun 13 23:32:46 2012) by ROOT version5.33/03
   TCanvas *cInt = new TCanvas("cInt", "cInt",666,160,700,500);
   gStyle->SetOptStat(0);
   cInt->Range(-11.09756,-0.127561,74.26829,0.2382927);
   cInt->SetFillColor(0);
   cInt->SetBorderMode(0);
   cInt->SetBorderSize(2);
   cInt->SetLeftMargin(0.13);
   cInt->SetRightMargin(0.05);
   cInt->SetTopMargin(0.05);
   cInt->SetBottomMargin(0.13);
   cInt->SetFrameBorderMode(0);
   cInt->SetFrameBorderMode(0);
   Double_t xAxis1[9] = {0, 5, 10, 20, 30, 40, 50, 60, 70}; 
   
   TH1D *vnInt = new TH1D("vnInt","",8, xAxis1);
   vnInt->SetMinimum(-0.08);
   vnInt->SetMaximum(0.22);
   vnInt->SetStats(0);
   vnInt->SetLineStyle(2);
   vnInt->GetXaxis()->SetTitle(" centrality percentile");
   vnInt->GetXaxis()->SetTitleSize(0.05);
   vnInt->GetYaxis()->SetTitle(" v_{n}");
   vnInt->GetYaxis()->SetNdivisions(505);
   vnInt->GetYaxis()->SetTitleSize(0.06);
   vnInt->Draw("");
   
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(5);
   grae->SetName("v2_WHDG");
   grae->SetTitle("Graph");

   Int_t ci;   // for color index setting
   ci = TColor::GetColor("#99cccc");
   grae->SetFillColor(ci);

   ci = TColor::GetColor("#0000ff");
   grae->SetLineColor(ci);
   grae->SetLineStyle(2);

   ci = TColor::GetColor("#99cccc");
   grae->SetMarkerColor(ci);
   grae->SetPoint(0,25,0.0559156);
   grae->SetPointError(0,0,0,0.003094881,0.001750282);
   grae->SetPoint(1,35,0.0657101);
   grae->SetPointError(1,0,0,0.005230785,0.001979145);
   grae->SetPoint(2,45,0.06971244);
   grae->SetPointError(2,0,0,0.007281904,0.002649773);
   grae->SetPoint(3,55,0.06053538);
   grae->SetPointError(3,0,0,0.006458776,0.002522587);
   grae->SetPoint(4,65,0.04610242);
   grae->SetPointError(4,0,0,0.006477253,0.003075159);
   
   TH1F *Graph_v2_WHDG1 = new TH1F("Graph_v2_WHDG1","Graph",100,21,69);
   Graph_v2_WHDG1->SetMinimum(0.03635146);
   Graph_v2_WHDG1->SetMaximum(0.07563592);
   Graph_v2_WHDG1->SetDirectory(0);
   Graph_v2_WHDG1->SetStats(0);
   Graph_v2_WHDG1->GetYaxis()->SetNdivisions(505);
   grae->SetHistogram(Graph_v2_WHDG1);
   
   grae->Draw("3");
   
   TGraph *graph = new TGraph(5);
   graph->SetName("v2b_WHDG");
   graph->SetTitle("Graph");
   graph->SetFillColor(1);

   ci = TColor::GetColor("#0000ff");
   graph->SetLineColor(ci);
   graph->SetLineStyle(2);
   graph->SetPoint(0,25,0.05591559579);
   graph->SetPoint(1,35,0.06571010211);
   graph->SetPoint(2,45,0.06971243881);
   graph->SetPoint(3,55,0.06053538029);
   graph->SetPoint(4,65,0.04610241589);
   
   TH1F *Graph_v2b_WHDG1 = new TH1F("Graph_v2b_WHDG1","Graph",100,21,69);
   Graph_v2b_WHDG1->SetMinimum(0.04374141);
   Graph_v2b_WHDG1->SetMaximum(0.07207344);
   Graph_v2b_WHDG1->SetDirectory(0);
   Graph_v2b_WHDG1->SetStats(0);
   Graph_v2b_WHDG1->GetYaxis()->SetNdivisions(505);
   graph->SetHistogram(Graph_v2b_WHDG1);
   
   graph->Draw(" l");
   
   grae = new TGraphAsymmErrors(8);
   grae->SetName("v2_syst_quad");
   grae->SetTitle("Graph");
   grae->SetFillColor(1);
   grae->SetMarkerStyle(20);
   grae->SetMarkerSize(1.2);
   grae->SetPoint(0,2.5,0.02326507);
   grae->SetPointError(0,0,0,0.01288447,0.007621022);
   grae->SetPoint(1,7.5,0.03821626);
   grae->SetPointError(1,0,0,0.01278374,0.005990243);
   grae->SetPoint(2,15,0.06364278);
   grae->SetPointError(2,0,0,0.01384455,0.004617843);
   grae->SetPoint(3,25,0.07246827);
   grae->SetPointError(3,0,0,0.01676203,0.005172315);
   grae->SetPoint(4,35,0.0695591);
   grae->SetPointError(4,0,0,0.02074983,0.006119452);
   grae->SetPoint(5,45,0.09679133);
   grae->SetPointError(5,0,0,0.02681535,0.008728469);
   grae->SetPoint(6,55,0.09471041);
   grae->SetPointError(6,0,0,0.03631868,0.01365573);
   grae->SetPoint(7,65,0.08335175);
   grae->SetPointError(7,0,0,0.05478074,0.02957645);
   
   TH1F *Graph_v2_syst_quad2 = new TH1F("Graph_v2_syst_quad2","Graph",100,0,71.25);
   Graph_v2_syst_quad2->SetMinimum(0.0001258374);
   Graph_v2_syst_quad2->SetMaximum(0.123183);
   Graph_v2_syst_quad2->SetDirectory(0);
   Graph_v2_syst_quad2->SetStats(0);
   Graph_v2_syst_quad2->GetYaxis()->SetNdivisions(505);
   grae->SetHistogram(Graph_v2_syst_quad2);
   
   grae->Draw("p");
   
   TGraphErrors *gre = new TGraphErrors(6);
   gre->SetName("v3_syst_quad");
   gre->SetTitle("Graph");
   gre->SetFillColor(1);
   gre->SetLineColor(2);
   gre->SetMarkerColor(2);
   gre->SetMarkerStyle(21);
   gre->SetMarkerSize(1.2);
   gre->SetPoint(0,3.5,0.02606927);
   gre->SetPointError(0,0,0.01366336);
   gre->SetPoint(1,8.5,0.02960644);
   gre->SetPointError(1,0,0.01454145);
   gre->SetPoint(2,16,0.03202531);
   gre->SetPointError(2,0,0.01348101);
   gre->SetPoint(3,26,0.00545194);
   gre->SetPointError(3,0,0.01686044);
   gre->SetPoint(4,36,-0.03126624);
   gre->SetPointError(4,0,0.02280925);
   gre->SetPoint(5,46,-0.02060315);
   gre->SetPointError(5,0,0.03378564);
   
   TH1F *Graph_v3_syst_quad1 = new TH1F("Graph_v3_syst_quad1","Graph",100,0,50.25);
   Graph_v3_syst_quad1->SetMinimum(-0.0643783);
   Graph_v3_syst_quad1->SetMaximum(0.05549583);
   Graph_v3_syst_quad1->SetDirectory(0);
   Graph_v3_syst_quad1->SetStats(0);
   Graph_v3_syst_quad1->GetYaxis()->SetNdivisions(505);
   gre->SetHistogram(Graph_v3_syst_quad1);
   
   gre->Draw("p");
   
   gre = new TGraphErrors(5);
   gre->SetName("v4_syst_quad");
   gre->SetTitle("Graph");
   gre->SetFillColor(1);
   gre->SetLineColor(4);
   gre->SetMarkerColor(4);
   gre->SetMarkerStyle(22);
   gre->SetMarkerSize(1.2);
   gre->SetPoint(0,2,-0.02976129);
   gre->SetPointError(0,0,0.03563585);
   gre->SetPoint(1,7,-0.0115773);
   gre->SetPointError(1,0,0.0384893);
   gre->SetPoint(2,14.5,0.02238646);
   gre->SetPointError(2,0,0.03119325);
   gre->SetPoint(3,24.5,0.03452682);
   gre->SetPointError(3,0,0.0397024);
   gre->SetPoint(4,34.5,0.0009197578);
   gre->SetPointError(4,0,0.05208576);
   
   TH1F *Graph_v4_syst_quad2 = new TH1F("Graph_v4_syst_quad2","Graph",100,0,37.75);
   Graph_v4_syst_quad2->SetMinimum(-0.07935977);
   Graph_v4_syst_quad2->SetMaximum(0.08819186);
   Graph_v4_syst_quad2->SetDirectory(0);
   Graph_v4_syst_quad2->SetStats(0);
   Graph_v4_syst_quad2->GetYaxis()->SetNdivisions(505);
   gre->SetHistogram(Graph_v4_syst_quad2);
   
   gre->Draw("p");
   
   gre = new TGraphErrors(6);
   gre->SetName("v42_syst_quad");
   gre->SetTitle("Graph");
   gre->SetFillColor(1);

   ci = TColor::GetColor("#009900");
   gre->SetLineColor(ci);

   ci = TColor::GetColor("#009900");
   gre->SetMarkerColor(ci);
   gre->SetMarkerStyle(23);
   gre->SetMarkerSize(1.2);
   gre->SetPoint(0,1.5,0.03777945);
   gre->SetPointError(0,0,0.02229117);
   gre->SetPoint(1,6.5,0.01429825);
   gre->SetPointError(1,0,0.01500416);
   gre->SetPoint(2,14,0.007410989);
   gre->SetPointError(2,0,0.01239348);
   gre->SetPoint(3,24,0.0001361442);
   gre->SetPointError(3,0,0.01457756);
   gre->SetPoint(4,34,0.02065072);
   gre->SetPointError(4,0,0.0185912);
   gre->SetPoint(5,44,-0.01572954);
   gre->SetPointError(5,0,0.02623625);
   
   TH1F *Graph_v42_syst_quad3 = new TH1F("Graph_v42_syst_quad3","Graph",100,0,48.25);
   Graph_v42_syst_quad3->SetMinimum(-0.05216943);
   Graph_v42_syst_quad3->SetMaximum(0.07027425);
   Graph_v42_syst_quad3->SetDirectory(0);
   Graph_v42_syst_quad3->SetStats(0);
   Graph_v42_syst_quad3->GetYaxis()->SetNdivisions(505);
   gre->SetHistogram(Graph_v42_syst_quad3);
   
   gre->Draw("p");
   
   gre = new TGraphErrors(6);
   gre->SetName("Graph_v24");
   gre->SetTitle("Graph_v24");
   gre->SetFillColor(1);
   gre->SetMarkerStyle(24);
   gre->SetMarkerSize(1.2);
   gre->SetPoint(0,8,0.02865856);
   gre->SetPointError(0,0,0.01892823);
   gre->SetPoint(1,15.5,0.04904932);
   gre->SetPointError(1,0,0.00890763);
   gre->SetPoint(2,25.5,0.07272248);
   gre->SetPointError(2,0,0.00884885);
   gre->SetPoint(3,35.5,0.09656275);
   gre->SetPointError(3,0,0.01160921);
   gre->SetPoint(4,45.5,0.1316306);
   gre->SetPointError(4,0,0.01961132);
   gre->SetPoint(5,55.5,0.08108739);
   gre->SetPointError(5,0,0.04998008);
   
   TH1F *Graph_Graph4 = new TH1F("Graph_v24","Graph_v24",100,3.25,60.25);
   Graph_Graph4->SetMinimum(0);
   Graph_Graph4->SetMaximum(0.1653931);
   Graph_Graph4->SetDirectory(0);
   Graph_Graph4->SetStats(0);
   Graph_Graph4->GetYaxis()->SetNdivisions(505);
   gre->SetHistogram(Graph_Graph4);
   
   gre->Draw("p");
   
   TLegend *leg = new TLegend(0.17,0.585,0.47,0.93,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextFont(132);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   TLegendEntry *entry=leg->AddEntry("v2_syst_quad","v_{2}{EP, |#Delta#eta|>2.0}","P");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1.2);
   entry=leg->AddEntry("Graph","v_{2}{4}","P");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(24);
   entry->SetMarkerSize(1.2);
   entry=leg->AddEntry("v3_syst_quad","v_{3}{EP, |#Delta#eta|>2.0}","P");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(2);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1.2);
   entry=leg->AddEntry("v4_syst_quad","v_{4/#Psi_{4}}{EP, |#Delta#eta|>2.0}","P");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(4);
   entry->SetMarkerStyle(22);
   entry->SetMarkerSize(1.2);
   entry=leg->AddEntry("v42_syst_quad","v_{4/#Psi_{2}}{EP, |#Delta#eta|>2.0}","P");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#009900");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(23);
   entry->SetMarkerSize(1.2);
   entry=leg->AddEntry("v2_WHDG","#pi^{0} v_{2} WHDG LHC","FL");

   ci = TColor::GetColor("#99cccc");
   entry->SetFillColor(ci);
   entry->SetFillStyle(1001);

   ci = TColor::GetColor("#0000ff");
   entry->SetLineColor(ci);
   entry->SetLineStyle(2);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("v2_WHDG","Extrapolation","");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   leg->Draw();
   TLatex *   tex = new TLatex(0.63,0.82,"10 < p_{t} < 20 GeV/c");
tex->SetNDC();
   tex->SetTextFont(132);
   tex->SetTextSize(0.048);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.56,0.89,"ALICE Pb-Pb #sqrt{s_{NN}}=2.76 TeV");
tex->SetNDC();
   tex->SetTextFont(132);
   tex->SetTextSize(0.048);
   tex->SetLineWidth(2);
   tex->Draw();
   cInt->Modified();
   cInt->cd();
   cInt->SetSelected(cInt);
}
