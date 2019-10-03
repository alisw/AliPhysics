//------------------------------------------------------------------------------
// makeCompYield.C
//
// create figure with ALICE vs. UA1 data
//------------------------------------------------------------------------------


void makeCompYield()
{
using namespace std;
 
TGraphErrors *graphYieldUa1 = new TGraphErrors(binsUa1,ptUa1,yieldUa1,errPtUa1,err2YieldUa1);
graphYieldUa1->SetMarkerColor(colorUa1);
graphYieldUa1->SetLineColor(colorUa1);
graphYieldUa1->SetMarkerStyle(26);

TGraphErrors *graphYieldAlice = new TGraphErrors(binsYieldAlice,ptYieldAlice,yieldAlice,0,err2YieldAlice);
graphYieldAlice->SetMarkerColor(colorAliceFit);
graphYieldAlice->SetLineColor(colorAliceFit);
graphYieldAlice->SetLineWidth(1);
graphYieldAlice->SetMarkerStyle(20);


TMultiGraph *all = new TMultiGraph("all","");


all->Add(graphYieldUa1);
//all->Add(graphCms);
//all->Add(graphAlice);
all->Add(graphYieldAlice);
//all->Draw("AP");

Double_t one[binsYieldAlice];  
for (Int_t i=0; i < binsYieldAlice; i++) { 
    one[i] = 1.0;
}

TGraphErrors *ratioUa1Alice = divide(graphYieldUa1,graphYieldAlice);
ratioUa1Alice->SetFillColor(colorUa1);
ratioUa1Alice->SetLineColor(colorUa1);
ratioUa1Alice->SetFillStyle(3345);


TGraphErrors *ratioAliceAlice = new TGraphErrors(binsYieldAlice,ptYieldAlice,one,0,relErr2YieldAlice);
ratioAliceAlice->SetFillColor(colorAliceErrors);
ratioAliceAlice->SetLineColor(20);

TCanvas *can4 = new TCanvas("can4","CompYield",520,700);

TPad *pad4_1 = new TPad("pad4_1","pad4_1",0.0,0.35,1.0,1.0);
setAttrib(pad4_1);

TPad *pad4_2 = new TPad("pad4_2","pad4_2",0.0,0.0,1.0,0.35);
setAttrib(pad4_2);

//
can4->cd();

pad4_1->Draw();
pad4_1->cd();
pad4_1->SetLogx();
pad4_1->SetLogy();
//all->Draw("AE4");
all->Draw("APZ");
all->GetXaxis()->SetTitle("p_{T} (GeV/c)");
all->GetXaxis()->SetTitleOffset(1.6);
all->GetYaxis()->SetTitleOffset(1.6);
all->GetYaxis()->SetTitle("1/N_{evt} 1/(2#pi p_{T}) (d^{2}N_{ch})/(dy dp_{T}) (GeV/c)^{-2}");
all->GetXaxis()->SetLimits(minPt,maxPt);
setAttrib(all);


graphYieldAlice->SetFillColor(2);

TLegend *l4_1 = new TLegend(0.200,0.03,0.65,0.203);
//TLegend *l4_1 = new TLegend(0.200,0.031,0.302,0.137);
l4_1->SetHeader("#sqrt{s} = 900 GeV");
l4_1->AddEntry(graphYieldAlice,"ALICE pp, NSD, | #eta | < 0.8","LP");
l4_1->AddEntry(graphYieldUa1,"UA1 p#bar{p}, NSD, | #eta | < 2.5","LP");
//l1->AddEntry(graphCms,"CMS |#eta|<2.4","LP");
//la->AddEntry(graphUA1,"UA1","p");
//la->AddEntry(graphCMS,"CMS","p");
//la->AddEntry(graphCMSred,"CMS, #eta<0.8","p");
//la->AddEntry(graphATLAS,"ATLAS","p");
l4_1->SetFillColor(0);
l4_1->SetLineColor(0);
l4_1->SetTextSize(legendTextSize);
l4_1->Draw();


//histo1->GetXaxis()->SetNdivisions(405);

can4->cd();
pad4_2->Draw();
pad4_2->cd();
pad4_2->SetLogx();

TMultiGraph *ratios = new TMultiGraph("ratios","");


ratios->Add(ratioAliceAlice);
//ratios->Add(ratioAtlasAlice);
ratios->Add(ratioUa1Alice);


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
//ratios->GetYaxis()->SetRangeUser(0.2,1.5);
//ratios->GetYaxis()->SetLimits(0.2,1.5);

ratios->GetYaxis()->SetLabelOffset(0.023);

graphYieldAlice->SetFillColor(2);

TLegend *l4_2 = new TLegend(0.200,0.323,0.65,0.537);
l4_2->AddEntry(ratioAliceAlice,"ALICE uncertainties","F");
l4_2->AddEntry(ratioUa1Alice,"UA1 / ALICE","F");
//l2->AddEntry(ratioCmsAlice,"CMS / ALICE","F");
l4_2->SetFillColor(0);
l4_2->SetLineColor(0);
l4_2->SetTextSize(legendTextSize);
l4_2->Draw();


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

logoPrelim(can4);
}
