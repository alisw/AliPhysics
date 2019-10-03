{
//=========Macro generated from canvas: cInt/cInt
//=========  (Fri Dec 28 05:23:36 2012) by ROOT version5.34/03
   TCanvas *cInt = new TCanvas("cInt", "cInt",66,52,700,500);
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
   
   TGraphErrors *gre = new TGraphErrors(6);
   gre->SetName("Graph");
   gre->SetTitle("Graph");

   ci = TColor::GetColor("#999999");
   gre->SetFillColor(ci);
   gre->SetLineColor(0);
   gre->SetLineWidth(-1);
   gre->SetMarkerStyle(24);
   gre->SetMarkerSize(1.1);
   gre->SetPoint(0,8,0.02865856);
   gre->SetPointError(0,1,0.000598431);
   gre->SetPoint(1,15.5,0.04904932);
   gre->SetPointError(1,1,0.000339083);
   gre->SetPoint(2,25.5,0.07272248);
   gre->SetPointError(2,1,0.000404397);
   gre->SetPoint(3,35.5,0.09656275);
   gre->SetPointError(3,1,0.000598475);
   gre->SetPoint(4,45.5,0.1316306);
   gre->SetPointError(4,1,0.00084029);
   gre->SetPoint(5,55.5,0.08108739);
   gre->SetPointError(5,1,0.00140364);
   
   TH1F *Graph_Graph1 = new TH1F("Graph_Graph1","Graph",100,2.05,61.45);
   Graph_Graph1->SetMinimum(0.01761905);
   Graph_Graph1->SetMaximum(0.142912);
   Graph_Graph1->SetDirectory(0);
   Graph_Graph1->SetStats(0);
   Graph_Graph1->GetYaxis()->SetNdivisions(505);
   gre->SetHistogram(Graph_Graph1);
   
   gre->Draw(" 2");
   
   grae = new TGraphAsymmErrors(8);
   grae->SetName("v2_sys");
   grae->SetTitle("Graph");

   ci = TColor::GetColor("#999999");
   grae->SetFillColor(ci);
   grae->SetLineColor(0);
   grae->SetMarkerStyle(20);
   grae->SetMarkerSize(1.1);
   grae->SetPoint(0,2.5,0.02326507);
   grae->SetPointError(0,0.5,0.5,0.01043052,0.0009306028);
   grae->SetPoint(1,7.5,0.03821626);
   grae->SetPointError(1,0.5,0.5,0.01139639,0.00152865);
   grae->SetPoint(2,15,0.06364278);
   grae->SetPointError(2,0.5,0.5,0.01329766,0.002545711);
   grae->SetPoint(3,25,0.07246827);
   grae->SetPointError(3,0.5,0.5,0.01620542,0.002898731);
   grae->SetPoint(4,35,0.0695591);
   grae->SetPointError(4,0.5,0.5,0.02002122,0.002782364);
   grae->SetPoint(5,45,0.09679133);
   grae->SetPointError(5,0.5,0.5,0.02564891,0.003871653);
   grae->SetPoint(6,55,0.09471041);
   grae->SetPointError(6,0.5,0.5,0.0338662,0.003788416);
   grae->SetPoint(7,65,0.08335175);
   grae->SetPointError(7,0.5,0.5,0.04623071,0.00333407);
   
   TH1F *Graph_v2_sys2 = new TH1F("Graph_v2_sys2","Graph",100,0,71.85);
   Graph_v2_sys2->SetMinimum(0.004051711);
   Graph_v2_sys2->SetMaximum(0.1094458);
   Graph_v2_sys2->SetDirectory(0);
   Graph_v2_sys2->SetStats(0);
   Graph_v2_sys2->GetYaxis()->SetNdivisions(505);
   grae->SetHistogram(Graph_v2_sys2);
   
   grae->Draw(" 2");
   
   gre = new TGraphErrors(6);
   gre->SetName("v3_sys");
   gre->SetTitle("Graph");

   ci = TColor::GetColor("#ffcccc");
   gre->SetFillColor(ci);
   gre->SetLineColor(0);
   gre->SetLineWidth(-1);
   gre->SetMarkerColor(2);
   gre->SetMarkerStyle(25);
   gre->SetMarkerSize(1.1);
   gre->SetPoint(0,3.5,0.02606927);
   gre->SetPointError(0,0.5,0.007164576);
   gre->SetPoint(1,8.5,0.02960644);
   gre->SetPointError(1,0.5,0.007797881);
   gre->SetPoint(2,16,0.03202531);
   gre->SetPointError(2,0.5,0.00899564);
   gre->SetPoint(3,26,0.00545194);
   gre->SetPointError(3,0.5,0.01084998);
   gre->SetPoint(4,36,-0.03126624);
   gre->SetPointError(4,0.5,0.01356202);
   gre->SetPoint(5,46,-0.02060315);
   gre->SetPointError(5,0.5,0.01727465);
   
   TH1F *Graph_v3_sys2 = new TH1F("Graph_v3_sys2","Graph",100,0,50.85);
   Graph_v3_sys2->SetMinimum(-0.05341319);
   Graph_v3_sys2->SetMaximum(0.04960587);
   Graph_v3_sys2->SetDirectory(0);
   Graph_v3_sys2->SetStats(0);
   Graph_v3_sys2->GetYaxis()->SetNdivisions(505);
   gre->SetHistogram(Graph_v3_sys2);
   
   gre->Draw(" 2");
   
   gre = new TGraphErrors(5);
   gre->SetName("v4_sys");
   gre->SetTitle("Graph");

   ci = TColor::GetColor("#ccccff");
   gre->SetFillColor(ci);
   gre->SetLineColor(0);
   gre->SetLineWidth(-1);
   gre->SetMarkerColor(4);
   gre->SetMarkerStyle(26);
   gre->SetMarkerSize(1.1);
   gre->SetPoint(0,2,-0.02976129);
   gre->SetPointError(0,0.5,0.007481394);
   gre->SetPoint(1,7,-0.0115773);
   gre->SetPointError(1,0.5,0.007735906);
   gre->SetPoint(2,14.5,0.02238646);
   gre->SetPointError(2,0.5,0.00906379);
   gre->SetPoint(3,24.5,0.03452682);
   gre->SetPointError(3,0.5,0.01120995);
   gre->SetPoint(4,34.5,0.0009197578);
   gre->SetPointError(4,0.5,0.01347619);
   
   TH1F *Graph_v4_sys3 = new TH1F("Graph_v4_sys3","Graph",100,0,38.35);
   Graph_v4_sys3->SetMinimum(-0.04554063);
   Graph_v4_sys3->SetMaximum(0.05403472);
   Graph_v4_sys3->SetDirectory(0);
   Graph_v4_sys3->SetStats(0);
   Graph_v4_sys3->GetYaxis()->SetNdivisions(505);
   gre->SetHistogram(Graph_v4_sys3);
   
   gre->Draw(" 2");
   
   gre = new TGraphErrors(6);
   gre->SetName("v42_sys");
   gre->SetTitle("Graph");

   ci = TColor::GetColor("#99ff99");
   gre->SetFillColor(ci);
   gre->SetLineColor(0);
   gre->SetLineWidth(-1);

   ci = TColor::GetColor("#009900");
   gre->SetMarkerColor(ci);
   gre->SetMarkerStyle(32);
   gre->SetMarkerSize(1.1);
   gre->SetPoint(0,1.5,0.03777945);
   gre->SetPointError(0,0.5,0.008698125);
   gre->SetPoint(1,6.5,0.01429825);
   gre->SetPointError(1,0.5,0.008898175);
   gre->SetPoint(2,14,0.007410989);
   gre->SetPointError(2,0.5,0.01021012);
   gre->SetPoint(3,24,0.0001361442);
   gre->SetPointError(3,0.5,0.01245007);
   gre->SetPoint(4,34,0.02065072);
   gre->SetPointError(4,0.5,0.01557679);
   gre->SetPoint(5,44,-0.01572954);
   gre->SetPointError(5,0.5,0.01996606);
   
   TH1F *Graph_v42_sys4 = new TH1F("Graph_v42_sys4","Graph",100,0,48.85);
   Graph_v42_sys4->SetMinimum(-0.04391292);
   Graph_v42_sys4->SetMaximum(0.05469489);
   Graph_v42_sys4->SetDirectory(0);
   Graph_v42_sys4->SetStats(0);
   Graph_v42_sys4->GetYaxis()->SetNdivisions(505);
   gre->SetHistogram(Graph_v42_sys4);
   
   gre->Draw(" 2");
   
   gre = new TGraphErrors(6);
   gre->SetName("Graph");
   gre->SetTitle("Graph");
   gre->SetFillColor(1);
   gre->SetMarkerStyle(24);
   gre->SetMarkerSize(1.1);
   gre->SetPoint(0,8,0.02865856);
   gre->SetPointError(0,0,0.0189188);
   gre->SetPoint(1,15.5,0.04904932);
   gre->SetPointError(1,0,0.00890117);
   gre->SetPoint(2,25.5,0.07272248);
   gre->SetPointError(2,0,0.0088396);
   gre->SetPoint(3,35.5,0.09656275);
   gre->SetPointError(3,0,0.0115938);
   gre->SetPoint(4,45.5,0.1316306);
   gre->SetPointError(4,0,0.0195933);
   gre->SetPoint(5,55.5,0.08108739);
   gre->SetPointError(5,0,0.0499604);
   
   TH1F *Graph_Graph5 = new TH1F("Graph_Graph5","Graph",100,3.25,60.25);
   Graph_Graph5->SetMinimum(0);
   Graph_Graph5->SetMaximum(0.1653723);
   Graph_Graph5->SetDirectory(0);
   Graph_Graph5->SetStats(0);
   Graph_Graph5->GetYaxis()->SetNdivisions(505);
   gre->SetHistogram(Graph_Graph5);
   
   gre->Draw("p");
   
   gre = new TGraphErrors(8);
   gre->SetName("v2_stat");
   gre->SetTitle("Graph");
   gre->SetFillColor(1);
   gre->SetMarkerStyle(20);
   gre->SetPoint(0,2.5,0.02326507);
   gre->SetPointError(0,0,0.007563991);
   gre->SetPoint(1,7.5,0.03821626);
   gre->SetPointError(1,0,0.005791912);
   gre->SetPoint(2,15,0.06364278);
   gre->SetPointError(2,0,0.003852769);
   gre->SetPoint(3,25,0.07246827);
   gre->SetPointError(3,0,0.004283713);
   gre->SetPoint(4,35,0.0695591);
   gre->SetPointError(4,0,0.005450335);
   gre->SetPoint(5,45,0.09679133);
   gre->SetPointError(5,0,0.007822818);
   gre->SetPoint(6,55,0.09471041);
   gre->SetPointError(6,0,0.01311971);
   gre->SetPoint(7,65,0.08335175);
   gre->SetPointError(7,0,0.02938793);
   
   TH1F *Graph_v2_stat6 = new TH1F("Graph_v2_stat6","Graph",100,0,71.25);
   Graph_v2_stat6->SetMinimum(0.005997218);
   Graph_v2_stat6->SetMaximum(0.1224435);
   Graph_v2_stat6->SetDirectory(0);
   Graph_v2_stat6->SetStats(0);
   Graph_v2_stat6->GetYaxis()->SetNdivisions(505);
   gre->SetHistogram(Graph_v2_stat6);
   
   gre->Draw("p");
   
   gre = new TGraphErrors(6);
   gre->SetName("v3_stat");
   gre->SetTitle("Graph");
   gre->SetFillColor(1);
   gre->SetLineColor(2);
   gre->SetMarkerColor(2);
   gre->SetMarkerStyle(25);
   gre->SetMarkerSize(1.1);
   gre->SetPoint(0,3.5,0.02606927);
   gre->SetPointError(0,0,0.01163427);
   gre->SetPoint(1,8.5,0.02960644);
   gre->SetPointError(1,0,0.01227383);
   gre->SetPoint(2,16,0.03202531);
   gre->SetPointError(2,0,0.01004072);
   gre->SetPoint(3,26,0.00545194);
   gre->SetPointError(3,0,0.01290552);
   gre->SetPoint(4,36,-0.03126624);
   gre->SetPointError(4,0,0.0183394);
   gre->SetPoint(5,46,-0.02060315);
   gre->SetPointError(5,0,0.02903542);
   
   TH1F *Graph_v3_stat7 = new TH1F("Graph_v3_stat7","Graph",100,0,50.25);
   Graph_v3_stat7->SetMinimum(-0.05880904);
   Graph_v3_stat7->SetMaximum(0.0512365);
   Graph_v3_stat7->SetDirectory(0);
   Graph_v3_stat7->SetStats(0);
   Graph_v3_stat7->GetYaxis()->SetNdivisions(505);
   gre->SetHistogram(Graph_v3_stat7);
   
   gre->Draw("p");
   
   gre = new TGraphErrors(5);
   gre->SetName("v4_stat");
   gre->SetTitle("Graph");
   gre->SetFillColor(1);
   gre->SetLineColor(4);
   gre->SetMarkerColor(4);
   gre->SetMarkerStyle(26);
   gre->SetMarkerSize(1.1);
   gre->SetPoint(0,2,-0.02976129);
   gre->SetPointError(0,0,0.03484168);
   gre->SetPoint(1,7,-0.0115773);
   gre->SetPointError(1,0,0.03770387);
   gre->SetPoint(2,14.5,0.02238646);
   gre->SetPointError(2,0,0.02984738);
   gre->SetPoint(3,24.5,0.03452682);
   gre->SetPointError(3,0,0.03808698);
   gre->SetPoint(4,34.5,0.0009197578);
   gre->SetPointError(4,0,0.05031221);
   
   TH1F *Graph_v4_stat8 = new TH1F("Graph_v4_stat8","Graph",100,0,37.75);
   Graph_v4_stat8->SetMinimum(-0.07832464);
   Graph_v4_stat8->SetMaximum(0.08633548);
   Graph_v4_stat8->SetDirectory(0);
   Graph_v4_stat8->SetStats(0);
   Graph_v4_stat8->GetYaxis()->SetNdivisions(505);
   gre->SetHistogram(Graph_v4_stat8);
   
   gre->Draw("p");
   
   gre = new TGraphErrors(6);
   gre->SetName("v42_stat");
   gre->SetTitle("Graph");
   gre->SetFillColor(1);

   ci = TColor::GetColor("#009900");
   gre->SetLineColor(ci);

   ci = TColor::GetColor("#009900");
   gre->SetMarkerColor(ci);
   gre->SetMarkerStyle(32);
   gre->SetMarkerSize(1.1);
   gre->SetPoint(0,1.5,0.03777945);
   gre->SetPointError(0,0,0.0205241);
   gre->SetPoint(1,6.5,0.01429825);
   gre->SetPointError(1,0,0.01208086);
   gre->SetPoint(2,14,0.007410989);
   gre->SetPointError(2,0,0.007025086);
   gre->SetPoint(3,24,0.0001361442);
   gre->SetPointError(3,0,0.007582955);
   gre->SetPoint(4,34,0.02065072);
   gre->SetPointError(4,0,0.01014872);
   gre->SetPoint(5,44,-0.01572954);
   gre->SetPointError(5,0,0.01702049);
   
   TH1F *Graph_v42_stat9 = new TH1F("Graph_v42_stat9","Graph",100,0,48.25);
   Graph_v42_stat9->SetMinimum(-0.04185539);
   Graph_v42_stat9->SetMaximum(0.06740891);
   Graph_v42_stat9->SetDirectory(0);
   Graph_v42_stat9->SetStats(0);
   Graph_v42_stat9->GetYaxis()->SetNdivisions(505);
   gre->SetHistogram(Graph_v42_stat9);
   
   gre->Draw("p");
   
   TLegend *leg = new TLegend(0.17,0.585,0.46,0.93,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextFont(132);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   TLegendEntry *entry=leg->AddEntry("v2_sys","v_{2}{EP, |#Delta#eta|>2.0}","PF");

   ci = TColor::GetColor("#999999");
   entry->SetFillColor(ci);
   entry->SetFillStyle(1001);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1.1);
   entry=leg->AddEntry("Graph","v_{2}{4}","PF");

   ci = TColor::GetColor("#999999");
   entry->SetFillColor(ci);
   entry->SetFillStyle(1001);
   entry->SetLineStyle(1);
   entry->SetLineWidth(-1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(24);
   entry->SetMarkerSize(1.1);
   entry=leg->AddEntry("v3_sys","v_{3}{EP, |#Delta#eta|>2.0}","PF");

   ci = TColor::GetColor("#ffcccc");
   entry->SetFillColor(ci);
   entry->SetFillStyle(1001);
   entry->SetLineStyle(1);
   entry->SetLineWidth(-1);
   entry->SetMarkerColor(2);
   entry->SetMarkerStyle(25);
   entry->SetMarkerSize(1.1);
   entry=leg->AddEntry("v4_sys","v_{4/#Psi_{4}}{EP, |#Delta#eta|>2.0}","PF");

   ci = TColor::GetColor("#ccccff");
   entry->SetFillColor(ci);
   entry->SetFillStyle(1001);
   entry->SetLineStyle(1);
   entry->SetLineWidth(-1);
   entry->SetMarkerColor(4);
   entry->SetMarkerStyle(26);
   entry->SetMarkerSize(1.1);
   entry=leg->AddEntry("v42_sys","v_{4/#Psi_{2}}{EP, |#Delta#eta|>2.0}","PF");

   ci = TColor::GetColor("#99ff99");
   entry->SetFillColor(ci);
   entry->SetFillStyle(1001);
   entry->SetLineStyle(1);
   entry->SetLineWidth(-1);

   ci = TColor::GetColor("#009900");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(32);
   entry->SetMarkerSize(1.1);
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
   TLatex *   tex = new TLatex(0.63,0.82,"10 < p_{T} < 20 GeV/c");
tex->SetNDC();
   tex->SetTextFont(132);
   tex->SetTextSize(0.048);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.56,0.89,"ALICE Pb-Pb #sqrt{s_{NN}} = 2.76 TeV");
tex->SetNDC();
   tex->SetTextFont(132);
   tex->SetTextSize(0.048);
   tex->SetLineWidth(2);
   tex->Draw();
   cInt->Modified();
   cInt->cd();
   cInt->SetSelected(cInt);
}
