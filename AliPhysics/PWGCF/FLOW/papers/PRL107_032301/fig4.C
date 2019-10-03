void fig4()
{
    
    
    myOptions();
    gROOT->ForceStyle();   
    TDatime now;
    int iDate = now.GetDate();
    int iYear=iDate/10000;
    int iMonth=(iDate%10000)/100;
    int iDay=iDate%100;
    char* cMonth[12]={"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"};
    char cStamp1[25],cStamp2[25];
    sprintf(cStamp1,"%i %s %i",iDay,cMonth[iMonth-1],iYear);
    sprintf(cStamp2,"%i/%.2d/%i",iDay,iMonth,iYear);
    
    
    
    //=========Macro generated from canvas: c/c
    //=========  (Wed Apr 20 21:53:35 2011) by ROOT version5.27/05
    TCanvas *c = new TCanvas("c", "correlation function",4,23,800,700);
    gStyle->SetOptStat(0);
    
    // ------------>Primitives in pad: c_1
    TPad *c_1 = new TPad("c_1", "c_1",0,0.,1,1);
    c_1->Draw();
    c_1->cd();
    c_1->Range(-3.246312,0.9897024,5.131268,1.019464);
    c_1->SetFillColor(0);
    c_1->SetTickx(1);
    c_1->SetTicky(1);
    c_1->SetLeftMargin(0.2);
    c_1->SetRightMargin(0.05);
    c_1->SetTopMargin(0.02);
    c_1->SetBottomMargin(0.15);
    
    
    Double_t xAxis5[37] = {-1.570796, -1.396263, -1.22173, -1.047198, -0.8726646, -0.6981317, -0.5235988, -0.3490659, -0.1745329, 0, 0.1745329, 0.3490659, 0.5235988, 0.6981317, 0.8726646, 1.047198, 1.22173, 1.396263, 1.570796, 1.745329, 1.919862, 2.094395, 2.268928, 2.443461, 2.617994, 2.792527, 2.96706, 3.141593, 3.316126, 3.490659, 3.665191, 3.839724, 4.014257, 4.18879, 4.363323, 4.537856, 4.712389}; 
    
    TH1D *hist = new TH1D("hist","",36, xAxis5);
    hist->SetBinContent(1,0.9969217);
    hist->SetBinContent(2,0.994212);
    hist->SetBinContent(3,0.9924679);
    hist->SetBinContent(4,0.9927132);
    hist->SetBinContent(5,0.9941483);
    hist->SetBinContent(6,0.99839);
    hist->SetBinContent(7,1.004135);
    hist->SetBinContent(8,1.009233);
    hist->SetBinContent(9,1.012218);
    hist->SetBinContent(10,1.012218);
    hist->SetBinContent(11,1.009233);
    hist->SetBinContent(12,1.004135);
    hist->SetBinContent(13,0.99839);
    hist->SetBinContent(14,0.9941483);
    hist->SetBinContent(15,0.9927132);
    hist->SetBinContent(16,0.9924679);
    hist->SetBinContent(17,0.994212);
    hist->SetBinContent(18,0.9969217);
    hist->SetBinContent(19,0.9993932);
    hist->SetBinContent(20,1.000563);
    hist->SetBinContent(21,1.001791);
    hist->SetBinContent(22,1.001702);
    hist->SetBinContent(23,1.001143);
    hist->SetBinContent(24,1.000711);
    hist->SetBinContent(25,1.000244);
    hist->SetBinContent(26,0.9996513);
    hist->SetBinContent(27,0.9998416);
    hist->SetBinContent(28,0.9998416);
    hist->SetBinContent(29,0.9996513);
    hist->SetBinContent(30,1.000244);
    hist->SetBinContent(31,1.000711);
    hist->SetBinContent(32,1.001143);
    hist->SetBinContent(33,1.001702);
    hist->SetBinContent(34,1.001791);
    hist->SetBinContent(35,1.000563);
    hist->SetBinContent(36,0.9993932);
    hist->SetBinError(1,0.0001929295);
    hist->SetBinError(2,0.0001921493);
    hist->SetBinError(3,0.0001918335);
    hist->SetBinError(4,0.000192209);
    hist->SetBinError(5,0.0001920719);
    hist->SetBinError(6,0.0001927733);
    hist->SetBinError(7,0.0001938162);
    hist->SetBinError(8,0.0001943988);
    hist->SetBinError(9,0.0001947595);
    hist->SetBinError(10,0.0001947595);
    hist->SetBinError(11,0.0001943988);
    hist->SetBinError(12,0.0001938162);
    hist->SetBinError(13,0.0001927733);
    hist->SetBinError(14,0.0001920719);
    hist->SetBinError(15,0.000192209);
    hist->SetBinError(16,0.0001918335);
    hist->SetBinError(17,0.0001921493);
    hist->SetBinError(18,0.0001929295);
    hist->SetBinError(19,0.0001931266);
    hist->SetBinError(20,0.00019325);
    hist->SetBinError(21,0.0001934904);
    hist->SetBinError(22,0.0001933326);
    hist->SetBinError(23,0.0001931764);
    hist->SetBinError(24,0.0001931073);
    hist->SetBinError(25,0.0001930693);
    hist->SetBinError(26,0.0001929318);
    hist->SetBinError(27,0.0001934329);
    hist->SetBinError(28,0.0001934329);
    hist->SetBinError(29,0.0001929318);
    hist->SetBinError(30,0.0001930693);
    hist->SetBinError(31,0.0001931073);
    hist->SetBinError(32,0.0001931764);
    hist->SetBinError(33,0.0001933326);
    hist->SetBinError(34,0.0001934904);
    hist->SetBinError(35,0.00019325);
    hist->SetBinError(36,0.0001931266);
    hist->SetMinimum(0.99);
    hist->SetMaximum(1.015);
    hist->SetEntries(4.826814e+08);
    hist->SetStats(0);
    hist->SetLineWidth(2);
    hist->GetXaxis()->SetTitle("#Delta#phi (rad.)");
    hist->GetXaxis()->SetLabelSize(0.05);
    hist->GetXaxis()->SetTitleSize(0.06);
    hist->GetYaxis()->SetTitle("C(#Delta#phi)");
    hist->GetYaxis()->SetLabelOffset(0.01);
    hist->GetYaxis()->SetLabelSize(0.05);
    hist->GetYaxis()->SetTitleSize(0.06);
    hist->GetYaxis()->SetTitleOffset(1.7);
    hist->GetZaxis()->SetLabelSize(0.05);
    hist->GetZaxis()->SetTitleSize(0.06);
    
    hist->SetMarkerStyle(kOpenCircle);
    //hist->Draw("");
    /*
     
     Double_t xAxis6[37] = {-1.570796, -1.396263, -1.22173, -1.047198, -0.8726646, -0.6981317, -0.5235988, -0.3490659, -0.1745329, 0, 0.1745329, 0.3490659, 0.5235988, 0.6981317, 0.8726646, 1.047198, 1.22173, 1.396263, 1.570796, 1.745329, 1.919862, 2.094395, 2.268928, 2.443461, 2.617994, 2.792527, 2.96706, 3.141593, 3.316126, 3.490659, 3.665191, 3.839724, 4.014257, 4.18879, 4.363323, 4.537856, 4.712389}; 
     
     TH1D *GetDistAndFlow_1__2 = new TH1D("GetDistAndFlow_1__2","",36, xAxis6);
     GetDistAndFlow_1__2->SetBinContent(1,0.9974817);
     GetDistAndFlow_1__2->SetBinContent(2,0.9945214);
     GetDistAndFlow_1__2->SetBinContent(3,0.9927582);
     GetDistAndFlow_1__2->SetBinContent(4,0.993248);
     GetDistAndFlow_1__2->SetBinContent(5,0.9938889);
     GetDistAndFlow_1__2->SetBinContent(6,0.9976373);
     GetDistAndFlow_1__2->SetBinContent(7,1.003408);
     GetDistAndFlow_1__2->SetBinContent(8,1.007627);
     GetDistAndFlow_1__2->SetBinContent(9,1.010439);
     GetDistAndFlow_1__2->SetBinContent(10,1.010439);
     GetDistAndFlow_1__2->SetBinContent(11,1.007627);
     GetDistAndFlow_1__2->SetBinContent(12,1.003408);
     GetDistAndFlow_1__2->SetBinContent(13,0.9976373);
     GetDistAndFlow_1__2->SetBinContent(14,0.9938889);
     GetDistAndFlow_1__2->SetBinContent(15,0.993248);
     GetDistAndFlow_1__2->SetBinContent(16,0.9927582);
     GetDistAndFlow_1__2->SetBinContent(17,0.9945214);
     GetDistAndFlow_1__2->SetBinContent(18,0.9974817);
     GetDistAndFlow_1__2->SetBinContent(19,1.000223);
     GetDistAndFlow_1__2->SetBinContent(20,1.000993);
     GetDistAndFlow_1__2->SetBinContent(21,1.002592);
     GetDistAndFlow_1__2->SetBinContent(22,1.002206);
     GetDistAndFlow_1__2->SetBinContent(23,1.001343);
     GetDistAndFlow_1__2->SetBinContent(24,1.000875);
     GetDistAndFlow_1__2->SetBinContent(25,1.000312);
     GetDistAndFlow_1__2->SetBinContent(26,0.9996157);
     GetDistAndFlow_1__2->SetBinContent(27,1.000308);
     GetDistAndFlow_1__2->SetBinContent(28,1.000308);
     GetDistAndFlow_1__2->SetBinContent(29,0.9996157);
     GetDistAndFlow_1__2->SetBinContent(30,1.000312);
     GetDistAndFlow_1__2->SetBinContent(31,1.000875);
     GetDistAndFlow_1__2->SetBinContent(32,1.001343);
     GetDistAndFlow_1__2->SetBinContent(33,1.002206);
     GetDistAndFlow_1__2->SetBinContent(34,1.002592);
     GetDistAndFlow_1__2->SetBinContent(35,1.000993);
     GetDistAndFlow_1__2->SetBinContent(36,1.000223);
     GetDistAndFlow_1__2->SetBinError(1,0.0003562986);
     GetDistAndFlow_1__2->SetBinError(2,0.000354732);
     GetDistAndFlow_1__2->SetBinError(3,0.0003541372);
     GetDistAndFlow_1__2->SetBinError(4,0.0003549467);
     GetDistAndFlow_1__2->SetBinError(5,0.000354547);
     GetDistAndFlow_1__2->SetBinError(6,0.0003558383);
     GetDistAndFlow_1__2->SetBinError(7,0.0003578351);
     GetDistAndFlow_1__2->SetBinError(8,0.0003588045);
     GetDistAndFlow_1__2->SetBinError(9,0.0003594313);
     GetDistAndFlow_1__2->SetBinError(10,0.0003594313);
     GetDistAndFlow_1__2->SetBinError(11,0.0003588045);
     GetDistAndFlow_1__2->SetBinError(12,0.0003578351);
     GetDistAndFlow_1__2->SetBinError(13,0.0003558383);
     GetDistAndFlow_1__2->SetBinError(14,0.000354547);
     GetDistAndFlow_1__2->SetBinError(15,0.0003549467);
     GetDistAndFlow_1__2->SetBinError(16,0.0003541372);
     GetDistAndFlow_1__2->SetBinError(17,0.000354732);
     GetDistAndFlow_1__2->SetBinError(18,0.0003562986);
     GetDistAndFlow_1__2->SetBinError(19,0.000356631);
     GetDistAndFlow_1__2->SetBinError(20,0.0003568194);
     GetDistAndFlow_1__2->SetBinError(21,0.0003573007);
     GetDistAndFlow_1__2->SetBinError(22,0.0003569495);
     GetDistAndFlow_1__2->SetBinError(23,0.0003566265);
     GetDistAndFlow_1__2->SetBinError(24,0.0003564958);
     GetDistAndFlow_1__2->SetBinError(25,0.0003564227);
     GetDistAndFlow_1__2->SetBinError(26,0.0003561458);
     GetDistAndFlow_1__2->SetBinError(27,0.0003572533);
     GetDistAndFlow_1__2->SetBinError(28,0.0003572533);
     GetDistAndFlow_1__2->SetBinError(29,0.0003561458);
     GetDistAndFlow_1__2->SetBinError(30,0.0003564227);
     GetDistAndFlow_1__2->SetBinError(31,0.0003564958);
     GetDistAndFlow_1__2->SetBinError(32,0.0003566265);
     GetDistAndFlow_1__2->SetBinError(33,0.0003569495);
     GetDistAndFlow_1__2->SetBinError(34,0.0003573007);
     GetDistAndFlow_1__2->SetBinError(35,0.0003568194);
     GetDistAndFlow_1__2->SetBinError(36,0.000356631);
     GetDistAndFlow_1__2->SetEntries(1.416101e+08);
     GetDistAndFlow_1__2->SetDirectory(0);
     GetDistAndFlow_1__2->SetStats(0);
     GetDistAndFlow_1__2->SetLineColor(4);
     GetDistAndFlow_1__2->SetLineWidth(2);
     GetDistAndFlow_1__2->GetXaxis()->SetTitle("#Delta#phi (rad.)");
     GetDistAndFlow_1__2->GetXaxis()->SetLabelSize(0.05);
     GetDistAndFlow_1__2->GetXaxis()->SetTitleSize(0.06);
     GetDistAndFlow_1__2->GetYaxis()->SetTitle("C(#Delta#phi)");
     GetDistAndFlow_1__2->GetYaxis()->SetLabelOffset(0.01);
     GetDistAndFlow_1__2->GetYaxis()->SetLabelSize(0.05);
     GetDistAndFlow_1__2->GetYaxis()->SetTitleSize(0.06);
     GetDistAndFlow_1__2->GetYaxis()->SetTitleOffset(1.3);
     GetDistAndFlow_1__2->GetZaxis()->SetLabelSize(0.05);
     GetDistAndFlow_1__2->GetZaxis()->SetTitleSize(0.06);
     GetDistAndFlow_1__2->SetMarkerStyle(kFullCircle);
     
     GetDistAndFlow_1__2->Draw("");
     */
    
    
    Double_t xAxis2[37] = {-1.570796, -1.396263, -1.22173, -1.047198, -0.8726646, -0.6981317, -0.5235988, -0.3490659, -0.1745329, 0, 0.1745329, 0.3490659, 0.5235988, 0.6981317, 0.8726646, 1.047198, 1.22173, 1.396263, 1.570796, 1.745329, 1.919862, 2.094395, 2.268928, 2.443461, 2.617994, 2.792527, 2.96706, 3.141593, 3.316126, 3.490659, 3.665191, 3.839724, 4.014257, 4.18879, 4.363323, 4.537856, 4.712389}; 
    
    TH1D *GetDistAndFlow_1__1 = new TH1D("GetDistAndFlow_1__1","",36, xAxis2);
    GetDistAndFlow_1__1->SetBinContent(1,0.9978028);
    GetDistAndFlow_1__1->SetBinContent(2,0.9945293);
    GetDistAndFlow_1__1->SetBinContent(3,0.9926859);
    GetDistAndFlow_1__1->SetBinContent(4,0.99349);
    GetDistAndFlow_1__1->SetBinContent(5,0.9934892);
    GetDistAndFlow_1__1->SetBinContent(6,0.9974816);
    GetDistAndFlow_1__1->SetBinContent(7,1.003301);
    GetDistAndFlow_1__1->SetBinContent(8,1.007196);
    GetDistAndFlow_1__1->SetBinContent(9,1.010189);
    GetDistAndFlow_1__1->SetBinContent(10,1.010189);
    GetDistAndFlow_1__1->SetBinContent(11,1.007196);
    GetDistAndFlow_1__1->SetBinContent(12,1.003301);
    GetDistAndFlow_1__1->SetBinContent(13,0.9974816);
    GetDistAndFlow_1__1->SetBinContent(14,0.9934892);
    GetDistAndFlow_1__1->SetBinContent(15,0.99349);
    GetDistAndFlow_1__1->SetBinContent(16,0.9926859);
    GetDistAndFlow_1__1->SetBinContent(17,0.9945293);
    GetDistAndFlow_1__1->SetBinContent(18,0.9978028);
    GetDistAndFlow_1__1->SetBinContent(19,1.000258);
    GetDistAndFlow_1__1->SetBinContent(20,1.000897);
    GetDistAndFlow_1__1->SetBinContent(21,1.002466);
    GetDistAndFlow_1__1->SetBinContent(22,1.001989);
    GetDistAndFlow_1__1->SetBinContent(23,1.001441);
    GetDistAndFlow_1__1->SetBinContent(24,1.000879);
    GetDistAndFlow_1__1->SetBinContent(25,1.000151);
    GetDistAndFlow_1__1->SetBinContent(26,0.99953);
    GetDistAndFlow_1__1->SetBinContent(27,1.000522);
    GetDistAndFlow_1__1->SetBinContent(28,1.000522);
    GetDistAndFlow_1__1->SetBinContent(29,0.99953);
    GetDistAndFlow_1__1->SetBinContent(30,1.000151);
    GetDistAndFlow_1__1->SetBinContent(31,1.000879);
    GetDistAndFlow_1__1->SetBinContent(32,1.001441);
    GetDistAndFlow_1__1->SetBinContent(33,1.001989);
    GetDistAndFlow_1__1->SetBinContent(34,1.002466);
    GetDistAndFlow_1__1->SetBinContent(35,1.000897);
    GetDistAndFlow_1__1->SetBinContent(36,1.000258);
    GetDistAndFlow_1__1->SetBinError(1,0.0004580698);
    GetDistAndFlow_1__1->SetBinError(2,0.0004559522);
    GetDistAndFlow_1__1->SetBinError(3,0.0004551655);
    GetDistAndFlow_1__1->SetBinError(4,0.0004563055);
    GetDistAndFlow_1__1->SetBinError(5,0.0004556771);
    GetDistAndFlow_1__1->SetBinError(6,0.0004573828);
    GetDistAndFlow_1__1->SetBinError(7,0.0004599905);
    GetDistAndFlow_1__1->SetBinError(8,0.0004611871);
    GetDistAndFlow_1__1->SetBinError(9,0.0004620077);
    GetDistAndFlow_1__1->SetBinError(10,0.0004620077);
    GetDistAndFlow_1__1->SetBinError(11,0.0004611871);
    GetDistAndFlow_1__1->SetBinError(12,0.0004599905);
    GetDistAndFlow_1__1->SetBinError(13,0.0004573828);
    GetDistAndFlow_1__1->SetBinError(14,0.0004556771);
    GetDistAndFlow_1__1->SetBinError(15,0.0004563055);
    GetDistAndFlow_1__1->SetBinError(16,0.0004551655);
    GetDistAndFlow_1__1->SetBinError(17,0.0004559522);
    GetDistAndFlow_1__1->SetBinError(18,0.0004580698);
    GetDistAndFlow_1__1->SetBinError(19,0.0004584189);
    GetDistAndFlow_1__1->SetBinError(20,0.0004586533);
    GetDistAndFlow_1__1->SetBinError(21,0.0004592669);
    GetDistAndFlow_1__1->SetBinError(22,0.0004587856);
    GetDistAndFlow_1__1->SetBinError(23,0.0004584035);
    GetDistAndFlow_1__1->SetBinError(24,0.0004582231);
    GetDistAndFlow_1__1->SetBinError(25,0.000458117);
    GetDistAndFlow_1__1->SetBinError(26,0.0004577663);
    GetDistAndFlow_1__1->SetBinError(27,0.000459302);
    GetDistAndFlow_1__1->SetBinError(28,0.000459302);
    GetDistAndFlow_1__1->SetBinError(29,0.0004577663);
    GetDistAndFlow_1__1->SetBinError(30,0.000458117);
    GetDistAndFlow_1__1->SetBinError(31,0.0004582231);
    GetDistAndFlow_1__1->SetBinError(32,0.0004584035);
    GetDistAndFlow_1__1->SetBinError(33,0.0004587856);
    GetDistAndFlow_1__1->SetBinError(34,0.0004592669);
    GetDistAndFlow_1__1->SetBinError(35,0.0004586533);
    GetDistAndFlow_1__1->SetBinError(36,0.0004584189);
    //GetDistAndFlow_1__1->SetMinimum(0.9430516);
    //GetDistAndFlow_1__1->SetMaximum(1.060699);
    GetDistAndFlow_1__1->SetEntries(8.569633e+07);
    GetDistAndFlow_1__1->SetDirectory(0);
    GetDistAndFlow_1__1->SetStats(0);
    GetDistAndFlow_1__1->SetLineColor(4);
    GetDistAndFlow_1__1->SetLineWidth(2);
    GetDistAndFlow_1__1->GetXaxis()->SetTitle("#Delta#phi (rad.)");
    GetDistAndFlow_1__1->GetXaxis()->SetLabelSize(0.05);
    GetDistAndFlow_1__1->GetXaxis()->SetTitleSize(0.06);
    GetDistAndFlow_1__1->GetYaxis()->SetTitle("C(#Delta#phi)");
    GetDistAndFlow_1__1->GetYaxis()->SetLabelOffset(0.01);
    GetDistAndFlow_1__1->GetYaxis()->SetLabelSize(0.05);
    GetDistAndFlow_1__1->GetYaxis()->SetTitleSize(0.06);
    GetDistAndFlow_1__1->GetYaxis()->SetTitleOffset(1.3);
    GetDistAndFlow_1__1->GetZaxis()->SetLabelSize(0.05);
    GetDistAndFlow_1__1->GetZaxis()->SetTitleSize(0.06);
    GetDistAndFlow_1__1->SetMarkerStyle(kFullCircle);
    GetDistAndFlow_1__1->SetMarkerSize(1.2)
    GetDistAndFlow_1__1->SetMarkerColor(kBlue);
    
    
    
    GetDistAndFlow_1__1->Draw(""); 
    
    TF1 *flowFunc = new TF1("flowFunc","[0]*(1+2*[1]*cos(2*x)+2*[2]*cos(3*x)+2*[3]*cos(4*x)+2*[4]*cos(5*x))",-1.570796,4.712389);
    flowFunc->SetFillColor(19);
    flowFunc->SetFillStyle(0);
    flowFunc->SetLineColor(2);
    flowFunc->SetLineWidth(3);
    flowFunc->SetLineStyle(1);
    flowFunc->SetChisquare(62.12704);
    flowFunc->SetNDF(35);
    flowFunc->GetXaxis()->SetLabelSize(0.05);
    flowFunc->GetXaxis()->SetTitleSize(0.06);
    flowFunc->GetYaxis()->SetLabelOffset(0.01);
    flowFunc->GetYaxis()->SetLabelSize(0.05);
    flowFunc->GetYaxis()->SetTitleSize(0.06);
    flowFunc->GetYaxis()->SetTitleOffset(1.3);
    flowFunc->SetParameter(0,0.999971);
    flowFunc->SetParError(0,5.941878e-05);
    flowFunc->SetParLimits(0,0,0);
    flowFunc->SetParameter(1,0.001639691);
    flowFunc->SetParError(1,0);
    flowFunc->SetParLimits(1,0.001639691,0.001639691);
    flowFunc->SetParameter(2,0.002702756);
    flowFunc->SetParError(2,0);
    flowFunc->SetParLimits(2,0.002702756,0.002702756);
    flowFunc->SetParameter(3,0.001165077);
    flowFunc->SetParError(3,0);
    flowFunc->SetParLimits(3,0.001165077,0.001165077);
    flowFunc->SetParameter(4,0.0001550361);
    flowFunc->SetParError(4,0);
    flowFunc->SetParLimits(4,0.0001550361,0.0001550361);
    flowFunc->Draw("SAME");
    
    TF1 *flowFuncPart = new TF1("flowFuncPart","[0]*(1+2*[1]*cos([2]*x))",-1.570796,4.712389);
    flowFuncPart->SetFillColor(19);
    flowFuncPart->SetFillStyle(0);
    flowFuncPart->SetLineWidth(1);
    flowFuncPart->GetXaxis()->SetLabelSize(0.05);
    flowFuncPart->GetXaxis()->SetTitleSize(0.06);
    flowFuncPart->GetYaxis()->SetLabelOffset(0.01);
    flowFuncPart->GetYaxis()->SetLabelSize(0.05);
    flowFuncPart->GetYaxis()->SetTitleSize(0.06);
    flowFuncPart->GetYaxis()->SetTitleOffset(1.3);
    flowFuncPart->SetParameter(0,0.999971);
    flowFuncPart->SetParError(0,0);
    flowFuncPart->SetParLimits(0,0,0);
    flowFuncPart->SetParameter(1,0);
    flowFuncPart->SetParError(1,0);
    flowFuncPart->SetParLimits(1,0,0);
    flowFuncPart->SetParameter(2,1);
    flowFuncPart->SetParError(2,0);
    flowFuncPart->SetParLimits(2,0,0);
    flowFuncPart->Draw("SAME");
    
    TF1 *flowFuncPart = new TF1("flowFuncPart","[0]*(1+2*[1]*cos([2]*x))",-1.570796,4.712389);
    flowFuncPart->SetFillColor(19);
    flowFuncPart->SetFillStyle(0);
    flowFuncPart->SetLineWidth(1);
    flowFuncPart->SetLineStyle(2);
    flowFuncPart->GetXaxis()->SetLabelSize(0.05);
    flowFuncPart->GetXaxis()->SetTitleSize(0.06);
    flowFuncPart->GetYaxis()->SetLabelOffset(0.01);
    flowFuncPart->GetYaxis()->SetLabelSize(0.05);
    flowFuncPart->GetYaxis()->SetTitleSize(0.06);
    flowFuncPart->GetYaxis()->SetTitleOffset(1.3);
    flowFuncPart->SetParameter(0,0.999971);
    flowFuncPart->SetParError(0,0);
    flowFuncPart->SetParLimits(0,0,0);
    flowFuncPart->SetParameter(1,0.001639691);
    flowFuncPart->SetParError(1,0);
    flowFuncPart->SetParLimits(1,0,0);
    flowFuncPart->SetParameter(2,2);
    flowFuncPart->SetParError(2,0);
    flowFuncPart->SetParLimits(2,0,0);
    flowFuncPart->Draw("SAME");
    
    TF1 *flowFuncPart = new TF1("flowFuncPart","[0]*(1+2*[1]*cos([2]*x))",-1.570796,4.712389);
    flowFuncPart->SetFillColor(19);
    flowFuncPart->SetFillStyle(0);
    flowFuncPart->SetLineWidth(1);
    flowFuncPart->SetLineStyle(3);
    flowFuncPart->GetXaxis()->SetLabelSize(0.05);
    flowFuncPart->GetXaxis()->SetTitleSize(0.06);
    flowFuncPart->GetYaxis()->SetLabelOffset(0.01);
    flowFuncPart->GetYaxis()->SetLabelSize(0.05);
    flowFuncPart->GetYaxis()->SetTitleSize(0.06);
    flowFuncPart->GetYaxis()->SetTitleOffset(1.3);
    flowFuncPart->SetParameter(0,0.999971);
    flowFuncPart->SetParError(0,0);
    flowFuncPart->SetParLimits(0,0,0);
    flowFuncPart->SetParameter(1,0.002702756);
    flowFuncPart->SetParError(1,0);
    flowFuncPart->SetParLimits(1,0,0);
    flowFuncPart->SetParameter(2,3);
    flowFuncPart->SetParError(2,0);
    flowFuncPart->SetParLimits(2,0,0);
    flowFuncPart->Draw("SAME");
    
    TF1 *flowFuncPart = new TF1("flowFuncPart","[0]*(1+2*[1]*cos([2]*x))",-1.570796,4.712389);
    flowFuncPart->SetFillColor(19);
    flowFuncPart->SetFillStyle(0);
    flowFuncPart->SetLineWidth(1);
    flowFuncPart->SetLineStyle(4);
    flowFuncPart->GetXaxis()->SetLabelSize(0.05);
    flowFuncPart->GetXaxis()->SetTitleSize(0.06);
    flowFuncPart->GetYaxis()->SetLabelOffset(0.01);
    flowFuncPart->GetYaxis()->SetLabelSize(0.05);
    flowFuncPart->GetYaxis()->SetTitleSize(0.06);
    flowFuncPart->GetYaxis()->SetTitleOffset(1.3);
    flowFuncPart->SetParameter(0,0.999971);
    flowFuncPart->SetParError(0,0);
    flowFuncPart->SetParLimits(0,0,0);
    flowFuncPart->SetParameter(1,0.001165077);
    flowFuncPart->SetParError(1,0);
    flowFuncPart->SetParLimits(1,0,0);
    flowFuncPart->SetParameter(2,4);
    flowFuncPart->SetParError(2,0);
    flowFuncPart->SetParLimits(2,0,0);
    flowFuncPart->Draw("SAME");
    
    TF1 *flowFuncPart = new TF1("flowFuncPart","[0]*(1+2*[1]*cos([2]*x))",-1.570796,4.712389);
    flowFuncPart->SetFillColor(19);
    flowFuncPart->SetFillStyle(0);
    flowFuncPart->SetLineWidth(1);
    flowFuncPart->SetLineStyle(5);
    flowFuncPart->GetXaxis()->SetLabelSize(0.05);
    flowFuncPart->GetXaxis()->SetTitleSize(0.06);
    flowFuncPart->GetYaxis()->SetLabelOffset(0.01);
    flowFuncPart->GetYaxis()->SetLabelSize(0.05);
    flowFuncPart->GetYaxis()->SetTitleSize(0.06);
    flowFuncPart->GetYaxis()->SetTitleOffset(1.3);
    flowFuncPart->SetParameter(0,0.999971);
    flowFuncPart->SetParError(0,0);
    flowFuncPart->SetParLimits(0,0,0);
    flowFuncPart->SetParameter(1,0.0001550361);
    flowFuncPart->SetParError(1,0);
    flowFuncPart->SetParLimits(1,0,0);
    flowFuncPart->SetParameter(2,5);
    flowFuncPart->SetParError(2,0);
    flowFuncPart->SetParLimits(2,0,0);
    flowFuncPart->Draw("SAME");
    
    //hist->SetMarkerStyle(kOpenCircle);
    //hist->Draw("SAME");
    
    
    GetDistAndFlow_1__1->Draw("SAME");
    TLatex *   tex = new TLatex(0.6,0.27,"2.0 < p_{t,trig} < 3.0");
    tex->SetNDC();
    tex->SetTextSize(0.04);
    tex->SetLineWidth(2);
    tex->Draw();
    tex = new TLatex(0.6,0.22,"1.0 < p_{t,assoc} < 2.0");
    tex->SetNDC();
    tex->SetTextSize(0.04);
    tex->SetLineWidth(2);
    tex->Draw();
    
    
    //TLegend *leg = new TLegend(0.48,0.82,0.91,0.95,NULL,"brNDC");
    TLegend *leg = new TLegend(0.48,0.78,0.91,0.95,"Centrality 0-1%, |#eta| < 0.8");
    myLegendSetUp(leg,0.04); 
    //leg->SetTextAlign(22);
    //leg->SetTextFont(22);
    //leg->SetTextSize(0.034);
    //leg->SetFillColor(10);
    //leg->SetFillStyle(0);
    //leg->AddEntry("hist","All #Delta#eta","p");
    leg->AddEntry("GetDistAndFlow_1__1","|#Delta#eta| > 1","p");
    leg->AddEntry("flowFuncPart","v_{2,3,4,5}{2, |#Delta#eta| > 1}","l");
    leg->Draw();
}


void myLegendSetUp(TLegend *currentLegend=0,float currentTextSize=0.07){
    currentLegend->SetTextFont(42);
    currentLegend->SetBorderSize(0);
    currentLegend->SetFillStyle(0);
    currentLegend->SetFillColor(0);
    currentLegend->SetMargin(0.25);
    currentLegend->SetTextSize(currentTextSize);
    currentLegend->SetEntrySeparation(0.5);
    return;
}

void myPadSetUp(TPad *currentPad, float currentLeft=0.11, float currentTop=0.04, float currentRight=0.04, float currentBottom=0.15){
    currentPad->SetLeftMargin(currentLeft);
    currentPad->SetTopMargin(currentTop);
    currentPad->SetRightMargin(currentRight);
    currentPad->SetBottomMargin(currentBottom);
    return;
}

void myGraphSetUp(TGraphErrors *currentGraph=0, Float_t currentMarkerSize = 1.0,
                  int currentMarkerStyle=20, int currentMarkerColor=0,
                  int currentLineStyle=1, int currentLineColor=0, int currentFillColor=0, int currentFillStyle=0)
{
    currentGraph->SetMarkerSize(currentMarkerSize);
    currentGraph->SetMarkerStyle(currentMarkerStyle);
    currentGraph->SetMarkerColor(currentMarkerColor);
    currentGraph->SetLineStyle(currentLineStyle);
    currentGraph->SetLineWidth(2);
    currentGraph->SetLineColor(currentLineColor);
    currentGraph->SetFillColor(currentFillColor);
    currentGraph->SetFillStyle(currentFillStyle);
    return;
}

void ShiftAlongXaxis(TGraphErrors *ge, Double_t shift)
{
    // Shift original TGraphErrors along x-axis by Double_t shift.
    
    if(!ge){cout<<"!ge"<<endl;exit(0);}
    
    Int_t nPoints = ge->GetN();
    Double_t x = 0.;
    Double_t y = 0.;
    for(Int_t p=0;p<nPoints;p++)
    { 
        ge->GetPoint(p,x,y);
        x+=shift;
        ge->SetPoint(p,x,y);
    } // end of for(Int_t p=0;p<nPoints;p++)
    
} // end of void ShiftAlongXaxis(TGraphErrors *ge, Double_t shift)

void myOptions(Int_t lStat=0){
    // Set gStyle
    int font = 42;
    // From plain
    gStyle->SetFrameBorderMode(0);
    gStyle->SetFrameFillColor(0);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetPadBorderMode(0);
    gStyle->SetPadColor(10);
    gStyle->SetCanvasColor(10);
    gStyle->SetTitleFillColor(10);
    gStyle->SetTitleBorderSize(1);
    gStyle->SetStatColor(10);
    gStyle->SetStatBorderSize(1);
    gStyle->SetLegendBorderSize(1);
    //
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
    gStyle->SetPalette(1,0); 
    if (lStat){
        gStyle->SetOptTitle(1);
        gStyle->SetOptStat(1111);
        gStyle->SetOptFit(1111);
    }
    else {
        gStyle->SetOptTitle(0);
        gStyle->SetOptStat(0);
        gStyle->SetOptFit(0);
    }
}

