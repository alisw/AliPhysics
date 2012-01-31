//------------------------------------------------------------------------------
// makeCompNSD.C
//
// create figure with ALICE comparison to CMS and ATLAS
//------------------------------------------------------------------------------



void makeCompNSD()
{

gROOT->LoadMacro("divide.C");
 
TGraphErrors *graphCms = new TGraphErrors(binsCms,ptCms,nsdCms,errPtCms,err2NsdCms);
graphCms->SetMarkerColor(colorCms);
graphCms->SetLineColor(colorCms);
graphCms->SetMarkerStyle(26);
 
TGraphAsymmErrors *graphAtlas = new TGraphAsymmErrors(binsAtlas,ptAtlas,nsdAtlas,lowErrPtAtlas,highErrPtAtlas,err2NsdAtlas,err2NsdAtlas);
graphAtlas->SetMarkerColor(colorAtlas);
graphAtlas->SetLineColor(colorAtlas);
graphAtlas->SetMarkerStyle(25);

TGraphErrors *graphAlice     = new TGraphErrors(binsNsdAlice,ptNsdAlice,nsdAlice,0,err2NsdAlice);
graphAlice->SetMarkerColor(colorAliceFit);
graphAlice->SetLineColor(colorAliceFit);
graphAlice->SetMarkerStyle(20);

TMultiGraph *all = new TMultiGraph("all","");


all->Add(graphAtlas);
all->Add(graphCms);
//all->Add(graphAlice);
all->Add(graphAlice);
//all->Draw("AP");
//setAttrib(all);


TGraphErrors *ratioAtlasAlice = divide(graphAtlas,graphAlice);
ratioAtlasAlice->SetFillColor(colorAtlas);
ratioAtlasAlice->SetLineColor(colorAtlas);
ratioAtlasAlice->SetFillStyle(3354);

TGraphErrors *ratioCmsAlice = divide(graphCms,graphAlice);
ratioCmsAlice->SetFillColor(colorCms);
ratioCmsAlice->SetLineColor(colorCms);
ratioCmsAlice->SetFillStyle(3345);

Double_t one[binsNsdAlice];  
for (Int_t i=0; i < binsNsdAlice; i++) { 
    one[i] = 1.0;
}

TGraphErrors *ratioAliceAlice = new TGraphErrors(binsNsdAlice,ptNsdAlice,one,0,relErr2NsdAlice);
ratioAliceAlice->SetFillColor(colorAliceErrors);
ratioAliceAlice->SetLineColor(20);

TCanvas *can3 = new TCanvas("can3","CompNSD",520,700);

TPad *pad3_1 = new TPad("pad3_1","pad3_1",0.0,0.35,1.0,1.0);
setAttrib(pad3_1);

TPad *pad3_2 = new TPad("pad3_2","pad3_2",0.0,0.0,1.0,0.35);
setAttrib(pad3_2);
//
can3->cd();

pad3_1->Draw();
pad3_1->cd();
pad3_1->SetLogx();
pad3_1->SetLogy();
//all->Draw("AE4");
all->Draw("APZ");
all->GetXaxis()->SetTitle("p_{T} [GeV/c]");
all->GetXaxis()->SetTitleOffset(1.6);
all->GetYaxis()->SetTitleOffset(1.6);
all->GetYaxis()->SetTitle("1/N_{evt} 1/(2#pi p_{T}) (d^{2}N_{ch})/(d#eta dp_{T}) (GeV/c)^{-2}");
all->GetXaxis()->SetLimits(minPt,maxPt);
setAttrib(all);
/*
all->SetTitle("");
all->UseCurrentStyle();
all->GetXaxis()->SetRangeUser(minPt,maxPt);
all->GetXaxis()->SetLabelOffset(0.0);
all->GetXaxis()->SetTitleOffset(1.0);
all->GetYaxis()->SetTitleOffset(1.5);
all->GetYaxis()->SetTitleSize(0.06);
all->GetYaxis()->SetLabelSize(0.06);
all->GetXaxis()->SetNdivisions(405);
all->GetXaxis()->SetLabelSize(0.0);
*/



graphAlice->SetFillColor(2);
TLegend *l3_1 = new TLegend(0.2,0.03,0.65,0.261);
l3_1->SetHeader("pp, NSD, #sqrt{s} = 900 GeV");
l3_1->AddEntry(graphAlice,"ALICE | #eta | < 0.8","LP");
l3_1->AddEntry(graphAtlas,"ATLAS | #eta | < 2.5","LP");
l3_1->AddEntry(graphCms,"CMS | #eta | < 2.4","LP");
//la->AddEntry(graphUA1,"UA1","p");
//la->AddEntry(graphCMS,"CMS","p");
//la->AddEntry(graphCMSred,"CMS, #eta<0.8","p");
//la->AddEntry(graphATLAS,"ATLAS","p");
l3_1->SetTextSize(legendTextSize);
l3_1->SetFillColor(0);
l3_1->SetLineColor(0);
l3_1->Draw();


//histo1->GetXaxis()->SetNdivisions(405);

can3->cd();
pad3_2->Draw();
pad3_2->cd();
pad3_2->SetLogx();

TMultiGraph *ratios = new TMultiGraph("ratios","");


ratios->Add(ratioAliceAlice);
ratios->Add(ratioAtlasAlice);
ratios->Add(ratioCmsAlice);

;
//ratio->SetTitle("ratio ALICE/ATLAS (different #eta)");
ratios->Draw("AE3");
ratios->GetXaxis()->SetLimits(minPt,maxPt);
ratios->GetXaxis()->SetTitle("p_{T} (GeV/c)");
ratios->GetXaxis()->SetTitleOffset(1.4);
ratios->GetYaxis()->SetTitleOffset(0.9);
ratios->GetYaxis()->SetTitle("ratio");
ratios->GetYaxis()->SetLabelSize(0.07);
ratios->GetXaxis()->SetLabelSize(0.07);
ratios->GetXaxis()->SetTitleSize(0.07);
ratios->GetYaxis()->SetTitleSize(0.07);
ratios->GetYaxis()->CenterTitle();
setAttrib(ratios);

ratios->GetYaxis()->SetLabelOffset(0.023);

//ratios->GetYaxis()->SetRangeUser(0.2,1.7);
//ratios->GetYaxis()->SetLimits(0.2,1.7);

/*
ratios->UseCurrentStyle();
ratios->GetXaxis()->SetTitleOffset(1.2);
ratios->GetYaxis()->SetTitleOffset(0.6);
ratios->GetXaxis()->SetTitleSize(0.135);
ratios->GetYaxis()->SetTitleSize(0.135);
ratios->GetXaxis()->SetLabelSize(0.135);
ratios->GetYaxis()->SetLabelSize(0.135);
ratios->GetXaxis()->SetLabelOffset(0.025);
ratios->GetYaxis()->SetLabelOffset(0.06);
ratios->GetXaxis()->SetTickLength(0.09);
ratios->GetYaxis()->SetLabelOffset(0.025);
ratios->GetXaxis()->SetNdivisions(405);
ratios->GetYaxis()->SetNdivisions(402);
ratios->GetXaxis()->SetRangeUser(minPt,maxPt);
*/

graphAlice->SetFillColor(2);
TLegend *l3_2 = new TLegend(0.2,0.323,0.65,0.637);
l3_2->AddEntry(ratioAliceAlice,"ALICE uncertainties","F");
l3_2->AddEntry(ratioAtlasAlice,"ATLAS / ALICE","F");
l3_2->AddEntry(ratioCmsAlice,"CMS / ALICE","F");
l3_2->SetFillColor(0);
l3_2->SetLineColor(0);
l3_2->SetTextSize(legendTextSize);
l3_2->Draw();


//aliceData->Draw("AE4");
TF1 *fOne = new TF1("one","1",0.1,20);
fOne->SetRange(0.1,20);
fOne->SetLineWidth(1);
fOne->Draw("SAME");
//ratio->GetYaxis()->SetTitle("ratio ALICE/ATLAS");
/*
histo1c->Draw();
histo1c->GetXaxis()->SetLabelSize(0.08);
histo1c->GetYaxis()->SetLabelSize(0.08);
histo1c->GetXaxis()->SetNdivisions(405);
histo1c->GetYaxis()->SetNdivisions(405);
*/

ratios->Draw("E3");
fOne->Draw("SAME");
logoPrelim(can3);
}
