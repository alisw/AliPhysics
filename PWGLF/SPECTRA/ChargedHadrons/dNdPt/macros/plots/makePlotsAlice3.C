//------------------------------------------------------------------------------
// makePlotsAlice3.C
//
// create figure with ALICE INEL data and ratios to mod Hagedorn and powerlaw fit
// fit ratios in seperate panels
//------------------------------------------------------------------------------


void makePlotsAlice3() 
{

TGraphErrors *graphStatNsdAlice = new TGraphErrors(binsNsdAlice,ptNsdAlice,nsdAlice,0,statNsdAlice);
graphStatNsdAlice->SetMarkerColor(colorAliceFit);
graphStatNsdAlice->SetLineColor(colorAliceFit);
graphStatNsdAlice->SetLineWidth(1);
graphStatNsdAlice->SetMarkerStyle(20);
graphStatNsdAlice->SetTitle("");

TGraphErrors *graphStatInelAlice = new TGraphErrors(binsInelAlice,ptInelAlice,inelAlice,0,statInelAlice);
graphStatInelAlice->SetMarkerColor(colorAliceFit);
graphStatInelAlice->SetLineColor(colorAliceFit);
graphStatInelAlice->SetLineWidth(1);
graphStatInelAlice->SetMarkerStyle(20);
graphStatInelAlice->SetTitle("");

TGraphErrors *graphSystNsdAlice = new TGraphErrors(binsNsdAlice,centerPtNsdAlice,nsdAlice,0,systNsdAlice);
graphSystNsdAlice->SetMarkerColor(colorAlice);
graphSystNsdAlice->SetFillColor(colorAliceErrors);
graphSystNsdAlice->SetLineColor(20);
graphSystNsdAlice->SetTitle("");
graphSystNsdAlice->SetMarkerStyle(20);

TGraphErrors *graphSystInelAlice =  new TGraphErrors(binsInelAlice,centerPtInelAlice,inelAlice,0,systInelAlice);
graphSystInelAlice->SetMarkerColor(colorAlice);
graphSystInelAlice->SetFillColor(colorAliceErrors);
graphSystInelAlice->SetLineColor(20);
graphSystInelAlice->SetTitle("");
graphSystInelAlice->SetMarkerStyle(20);





TCanvas *can23 = new TCanvas("can23","PlotsAlice3_INEL",520,876);

TPad *pad23_1 = new TPad("pad23_1","pad23_1",0.0,425.0/876.0,1.0,1.0);
setAttrib(pad23_1);
pad23_1->SetBottomMargin(0);
//pad23_1->SetTopMargin(0);
pad23_1->SetBorderSize(0);



TPad *pad23_2 = new TPad("pad23_2","pad23_2",0.0,(425.0-175.0)/876.0,1.0,425.0/876.0);
setAttrib(pad23_2);.
pad23_2->SetTopMargin(0.0);
pad23_2->SetBottomMargin(0.0);
pad23_2->SetBorderSize(0);

TPad *pad23_3 = new TPad("pad23_2","pad23_2",0.0,0.0,1.0,(425.0-175.0)/876.0);
setAttrib(pad23_3);
pad23_3->SetTopMargin(0.0);
//pad23_3->SetBottomMargin(0);
pad23_3->SetBottomMargin(0.30);
pad23_3->SetBorderSize(0);


can23->cd();

pad23_1->Draw();
pad23_1->cd();
pad23_1->SetLogx();
pad23_1->SetLogy();

graphStatInelAlice->Draw("APZ");
//graphSystInelAlice->Draw("AE3");
graphStatInelAlice->GetXaxis()->SetLimits(minPt,maxPt);
graphStatInelAlice->GetXaxis()->SetTitle("p_{T} (GeV/c)");
graphStatInelAlice->GetXaxis()->SetTitleOffset(1.6);
graphStatInelAlice->GetYaxis()->SetTitleOffset(1.6);
graphStatInelAlice->GetYaxis()->SetTitle("1/N_{evt} 1/(2#pi p_{T}) (d^{2}N_{ch})/(d#eta dp_{T}) (GeV/c)^{-2}");
graphStatInelAlice->GetYaxis()->SetLimits(1.01e-7,10);
graphStatInelAlice->GetYaxis()->SetRangeUser(1.01e-7,10);
setAttrib(graphStatInelAlice);
graphStatInelAlice->GetYaxis()->SetTitleOffset(2.8);

//ptfunction2 = new TF1("ptfunction2","[0]*(1+(sqrt(0.14*0.14+x*x)-0.14)/([1]*[2]))^ ((-1) * [1])",0,12.0); 
ptfunction2 = new TF1("ptfunction2","[0]*(x/sqrt(0.14*0.14+x*x))*(1+(x)/([2]))^((-1) * [1])",0.0,10.0); 

powerlaw4 = new TF1("powerlaw4","[0]*x^[1]",0.0,10.0);

ptfunction2->SetParameter(0,54.931915);
ptfunction2->SetParameter(1,7.931081);
ptfunction2->SetParameter(2,0.148710);
ptfunction2->SetLineWidth(2);

graphStatInelAlice->Fit("ptfunction2","MN","",0.0,10.0);
graphStatInelAlice->Fit("powerlaw4","MN","",3.0,10.0);

ptfunction2->SetLineColor(colorCms);
powerlaw4->SetLineWidth(2);
powerlaw4->SetLineColor(colorAtlas);

powerlaw4->SetLineStyle(7); //was 7 for dashed
powerlaw4->SetRange(1,3);
powerlaw4->DrawCopy("SAME");

ptfunction2->SetLineStyle(1); //was 7 for dashed
ptfunction2->SetRange(4,10);
ptfunction2->DrawCopy("SAME");

powerlaw4->SetLineStyle(7);
powerlaw4->SetRange(3,10);
powerlaw4->DrawCopy("SAME");

ptfunction2->SetLineStyle(1);
ptfunction2->SetRange(0,4);
ptfunction2->DrawCopy("SAME");



//TLegend *leg23_1 = new TLegend(0.2,0.145,0.65,0.376);
TLegend *leg23_1 = new TLegend(0.2,0.03,0.65,0.261);
leg23_1->SetHeader("pp, INEL, #sqrt{s} = 900 GeV, | #eta | < 0.8");
//leg23_1->AddEntry("","ALICE pp, INEL","");
//leg23_1->AddEntry("","#sqrt{s} = 900 GeV, | #eta | < 0.8","");
//leg23_1->AddEntry("","#LT p_{T} #GT = (469 #pm 1 #pm 12) MeV/c","");
leg23_1->AddEntry(graphStatInelAlice,"ALICE data","LP");
leg23_1->AddEntry(ptfunction2,"mod. Hagedorn fit","L");
leg23_1->AddEntry(powerlaw4,"power law fit, p_{T} > 3 GeV/c","L");
//leg23_1->AddEntry("","n (3-10 GeV/c) = 6.63 #pm 0.12 #pm 0.01","");

leg23_1->SetFillColor(0);
leg23_1->SetLineColor(0);
leg23_1->SetTextSize(legendTextSize);
leg23_1->Draw();
/*
TLegend *leg23_11 = new TLegend(0.2,0.03,0.3,0.088);
leg23_11->SetFillColor(0);
leg23_11->SetLineColor(0);
leg23_11->SetTextSize(legendTextSize);
leg23_11->AddEntry("","A*p_{T}^{-n}; n = 6.63 #pm 0.12","");
leg23_11->Draw();
*/

can23->cd();
/*
pad23_11->Draw();
pad23_11->cd();
*/
pad23_2->SetTopMargin(0.0);
pad23_2->SetBottomMargin(0.0);
pad23_2->Draw();
pad23_2->cd();
pad23_2->SetLogx();


Double_t one[binsInelAlice];  
for (Int_t i=0; i < binsInelAlice; i++) { 
    one[i] = 1.0;
}

TGraphErrors *systInelAlice = new TGraphErrors(binsInelAlice*2);
for (int i=0; i < binsInelAlice; i++) {
   systInelAlice->SetPoint(i*2,lowPtInelAlice[i],1);
   systInelAlice->SetPointError(i*2,0,relSystInelAlice[i]);
   systInelAlice->SetPoint(i*2+1,highPtInelAlice[i],1);   
   systInelAlice->SetPointError(i*2+1,0,relSystInelAlice[i]);
}
//TGraphErrors *ratioAliceAlice = new TGraphErrors(binsInelAlice,centerPtInelAlice,one,0,relErr2InelAlice);
   
systInelAlice->Draw("AE3");
systInelAlice->SetFillColor(colorAliceErrors);
systInelAlice->SetLineColor(colorAliceErrors);
systInelAlice->SetTitle("");
systInelAlice->SetLineWidth(0);
systInelAlice->GetXaxis()->SetLimits(minPt,maxPt);
systInelAlice->GetXaxis()->SetTitle("p_{T} (GeV/c)");
systInelAlice->GetXaxis()->SetTitleOffset(1.4);
systInelAlice->GetYaxis()->SetTitleOffset(0.9);
systInelAlice->GetYaxis()->SetTitle("fit / data");
systInelAlice->GetYaxis()->CenterTitle();
systInelAlice->GetYaxis()->SetLabelSize(0.07);
systInelAlice->GetXaxis()->SetLabelSize(0.07);
systInelAlice->GetXaxis()->SetTitleSize(0.07);
systInelAlice->GetYaxis()->SetTitleSize(0.07);
systInelAlice->GetYaxis()->SetRangeUser(0.1,1.9);
systInelAlice->GetYaxis()->SetLimits(0.1,1.9);
setAttrib(systInelAlice);
systInelAlice->GetYaxis()->SetNdivisions(505);
systInelAlice->SetFillColor(colorAliceErrors);
systInelAlice->SetLineColor(20);
//systInelAlice->Clone()->Draw("L");
systInelAlice->GetYaxis()->SetTitleOffset(2.8);

TGraphErrors *ratioFitInelAlice = new TGraphErrors(binsInelAlice);
for (int i=0; i < binsInelAlice; i++) {
   if (ptInelAlice[i] > 4.0) break;
   ratioFitInelAlice->SetPoint(i,ptInelAlice[i],ptfunction2->Eval(ptInelAlice[i])/inelAlice[i]);
   ratioFitInelAlice->SetPointError(i,0,statInelAlice[i]*ptfunction2->Eval(ptInelAlice[i])/(inelAlice[i]*inelAlice[i]));
}

TGraphErrors *ratioPFitInelAlice = new TGraphErrors(binsInelAlice);
for (int i=0; i < binsInelAlice; i++) {
   ratioPFitInelAlice->SetPoint(i,ptInelAlice[i],powerlaw4->Eval(ptInelAlice[i])/inelAlice[i]);
   ratioPFitInelAlice->SetPointError(i,0,statInelAlice[i]*powerlaw4->Eval(ptInelAlice[i])/(inelAlice[i]*inelAlice[i]));
}

TGraphErrors *ratioFitInelAliceNF = new TGraphErrors(binsInelAlice);
for (int i=0; i < binsInelAlice; i++) {
   ratioFitInelAliceNF->SetPoint(i,ptInelAlice[i],ptfunction2->Eval(ptInelAlice[i])/inelAlice[i]);
   ratioFitInelAliceNF->SetPointError(i,0,statInelAlice[i]*ptfunction2->Eval(ptInelAlice[i])/(inelAlice[i]*inelAlice[i]));
}

TGraphErrors *ratioPFitInelAliceNF = new TGraphErrors(binsInelAlice);
for (int i=0; i < binsInelAlice; i++) {
   if (ptInelAlice[i] > 3.0) break;
   ratioPFitInelAliceNF->SetPoint(i,ptInelAlice[i],powerlaw4->Eval(ptInelAlice[i])/inelAlice[i]);
   ratioPFitInelAliceNF->SetPointError(i,0,statInelAlice[i]*powerlaw4->Eval(ptInelAlice[i])/(inelAlice[i]*inelAlice[i]));
}


systInelAlice->GetYaxis()->SetLabelOffset(0.023);

setAttrib(ratioFitInelAliceNF);
ratioFitInelAliceNF->SetLineWidth(1);
ratioFitInelAliceNF->SetMarkerStyle(20);
ratioFitInelAliceNF->SetMarkerColor(colorCms);
ratioFitInelAliceNF->SetLineColor(colorCms);
ratioFitInelAliceNF->Draw("PZ");

setAttrib(ratioFitInelAlice);
ratioFitInelAlice->SetLineWidth(1);
ratioFitInelAlice->SetMarkerStyle(20);
ratioFitInelAlice->SetMarkerColor(colorCms);
ratioFitInelAlice->SetLineColor(colorCms);
ratioFitInelAlice->Draw("PZ");


setAttrib(ratioPFitInelAlice);
ratioPFitInelAlice->SetLineWidth(1);
ratioPFitInelAlice->SetMarkerStyle(20);
ratioPFitInelAlice->SetMarkerColor(colorAtlas);
ratioPFitInelAlice->SetLineColor(colorAtlas);
//ratioPFitInelAlice->Draw("PZ");





TF1 *fOne = new TF1("one","1",0.1,20);
fOne->SetRange(0.1,20);
fOne->SetLineWidth(1);
fOne->Draw("SAME");

TLegend *l23_2 = new TLegend(0.200,0.053,0.65,0.353);
l23_2->AddEntry(systInelAlice,"ALICE systematic uncertainties","F");
l23_2->AddEntry(ratioFitInelAlice,"mod. Hagedorn","LP");
l23_2->SetFillColor(0);
l23_2->SetLineColor(0);
l23_2->SetTextSize(legendTextSize);
l23_2->Draw();

//ratioPFitInelAlice->Draw("PZ");



can23->cd();
// pad23_11->cd();

pad23_3->Draw();
pad23_3->cd();
pad23_3->SetLogx();

systInelAlice->Draw("AE3");
systInelAlice->SetFillColor(colorAliceErrors);
systInelAlice->SetLineColor(colorAliceErrors);
systInelAlice->SetTitle("");
systInelAlice->SetLineWidth(0);
systInelAlice->GetXaxis()->SetLimits(minPt,maxPt);
systInelAlice->GetXaxis()->SetTitle("p_{T} (GeV/c)");
systInelAlice->GetXaxis()->SetTitleOffset(3.6);
systInelAlice->GetYaxis()->SetTitleOffset(0.9);
systInelAlice->GetYaxis()->SetTitle("fit / data");
systInelAlice->GetYaxis()->CenterTitle();
systInelAlice->GetYaxis()->SetLabelSize(0.07);
systInelAlice->GetXaxis()->SetLabelSize(0.07);
systInelAlice->GetXaxis()->SetTitleSize(0.07);
systInelAlice->GetYaxis()->SetTitleSize(0.07);
systInelAlice->GetYaxis()->SetRangeUser(0.1,1.9);
systInelAlice->GetYaxis()->SetLimits(0.1,1.9);
systInelAlice->GetYaxis()->SetNdivisions(505);
setAttrib(systInelAlice);
systInelAlice->SetFillColor(colorAliceErrors);
systInelAlice->SetLineColor(20);
systInelAlice->GetYaxis()->SetTitleOffset(2.8);
systInelAlice->GetXaxis()->SetTitleOffset(4.6);

setAttrib(ratioPFitInelAlice);
ratioPFitInelAlice->SetLineWidth(1);
ratioPFitInelAlice->SetMarkerStyle(20);
ratioPFitInelAlice->SetMarkerColor(colorAtlas);
ratioPFitInelAlice->SetLineColor(colorAtlas);
ratioPFitInelAlice->Draw("PZ");

setAttrib(ratioPFitInelAliceNF);
ratioPFitInelAliceNF->SetLineWidth(1);
ratioPFitInelAliceNF->SetMarkerStyle(20);
ratioPFitInelAliceNF->SetMarkerColor(kWhite);
ratioPFitInelAliceNF->SetLineColor(kWhite);
ratioPFitInelAliceNF->DrawClone("PZ");

setAttrib(ratioPFitInelAliceNF);
ratioPFitInelAliceNF->SetLineWidth(1);
ratioPFitInelAliceNF->SetMarkerStyle(24);
ratioPFitInelAliceNF->SetMarkerColor(colorAtlas);
ratioPFitInelAliceNF->SetLineColor(colorAtlas);
ratioPFitInelAliceNF->Draw("PZ");

fOne->Draw("SAME");


TLegend *l23_3 = new TLegend(0.200,0.323,0.65,0.537);
l23_3->AddEntry(systInelAlice,"ALICE systematic uncertainties","F");
l23_3->AddEntry(ratioPFitInelAlice,"power law","LP");
l23_3->SetFillColor(0);
l23_3->SetLineColor(0);
l23_3->SetTextSize(legendTextSize);
l23_3->Draw();

logoPrelim(can23);
}
