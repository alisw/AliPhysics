//------------------------------------------------------------------------------
// makeCompInel.C
//
// create figure with MC model comparisons
//------------------------------------------------------------------------------




void makeCompInel()
{
using namespace std;

TGraphErrors *graphInelAliceFit = new TGraphErrors(binsInelAlice,ptInelAlice,inelAlice,0,err2InelAlice);
graphInelAliceFit->SetMarkerColor(colorAliceFit);
graphInelAliceFit->SetLineColor(colorAliceFit);
graphInelAliceFit->SetLineWidth(1);
graphInelAliceFit->SetMarkerStyle(20);

TGraphErrors *graphInelAlice     = new TGraphErrors(binsInelAlice,centerPtInelAlice,inelAlice,0,0);
graphInelAlice->SetMarkerColor(colorAlice);
graphInelAlice->SetLineColor(colorAlice);
graphInelAlice->SetMarkerStyle(20);

TGraphErrors *graphPhojet = new TGraphErrors(binsPhojet,ptPhojet,inelPhojet,errPtPhojet,errInelPhojet);
graphPhojet->SetMarkerColor(colorPhojet);
graphPhojet->SetLineColor(colorPhojet);
graphPhojet->SetMarkerStyle(10);
graphPhojet->SetLineWidth(2);
graphPhojet->SetLineStyle(1);

TGraphErrors *graphPythia109 = new TGraphErrors(binsPythia109,ptPythia109,inelPythia109,errPtPythia109,errInelPythia109);
graphPythia109->SetMarkerColor(colorPythia109);
graphPythia109->SetLineColor(colorPythia109);
graphPythia109->SetMarkerStyle(10);
graphPythia109->SetLineWidth(2);
graphPythia109->SetLineStyle(7);

TGraphErrors *graphPythia306 = new TGraphErrors(binsPythia306,ptPythia306,inelPythia306,errPtPythia306,errInelPythia306);
graphPythia306->SetMarkerColor(colorPythia306);
graphPythia306->SetLineColor(colorPythia306);
graphPythia306->SetMarkerStyle(10);
graphPythia306->SetLineWidth(2);
graphPythia306->SetLineStyle(3);

TGraphErrors *graphPythia320 = new TGraphErrors(binsPythia320,ptPythia320,inelPythia320,errPtPythia320,errInelPythia320);
graphPythia320->SetMarkerColor(colorPythia320);
graphPythia320->SetLineColor(colorPythia320);
graphPythia320->SetMarkerStyle(10);
graphPythia320->SetLineWidth(2);
graphPythia320->SetLineStyle(9);

Double_t one[binsInelAlice];  
for (Int_t i=0; i < binsInelAlice; i++) { 
    one[i] = 1.0;
}

TGraphErrors *ratioAliceAlice = new TGraphErrors(binsInelAlice*2);
for (int i=0; i < binsInelAlice; i++) {
   ratioAliceAlice->SetPoint(i*2,lowPtInelAlice[i],1);
   ratioAliceAlice->SetPointError(i*2,0,relErr2InelAlice[i]);
   ratioAliceAlice->SetPoint(i*2+1,highPtInelAlice[i],1);   
   ratioAliceAlice->SetPointError(i*2+1,0,relErr2InelAlice[i]);
}
//TGraphErrors *ratioAliceAlice = new TGraphErrors(binsInelAlice,centerPtInelAlice,one,0,relErr2InelAlice);
ratioAliceAlice->SetFillColor(colorAliceErrors);
ratioAliceAlice->SetLineColor(20);

TGraph *ratioPhojetAlice = new TGraph(binsInelAlice*2);
for (int i=0; i < binsInelAlice; i++) {
   ratioPhojetAlice->SetPoint(i*2,lowPtInelAlice[i],inelPhojet[i+3]/inelAlice[i]);
   ratioPhojetAlice->SetPoint(i*2+1,highPtInelAlice[i],inelPhojet[i+3]/inelAlice[i]);
   cout << "inelAlice[" << i << "]           = " << inelAlice[i] <<endl;
   cout << "centerPtInelAlice[" << i << "]   = " << centerPtInelAlice[i] <<endl;
   cout << "inelPhojet[" << i+3 << "]          = " << inelPhojet[i+3] <<endl;
   cout << "ptPhojet[" << i+3 << "]            = " << ptPhojet[i+3] <<endl;
}
//TGraphErrors *ratioPhojetAlice = divide(graphPhojet,graphInelAlice);
ratioPhojetAlice->SetFillColor(colorPhojet);
ratioPhojetAlice->SetLineColor(colorPhojet);
ratioPhojetAlice->SetLineWidth(2);
ratioPhojetAlice->SetLineStyle(1);
ratioPhojetAlice->SetMarkerStyle(26);
ratioPhojetAlice->SetFillStyle(3354);


TGraph *ratioPythia109Alice = new TGraph(binsInelAlice*2);
for (int i=0; i < binsInelAlice; i++) {
   ratioPythia109Alice->SetPoint(i*2,lowPtInelAlice[i],inelPythia109[i+3]/inelAlice[i]);
   ratioPythia109Alice->SetPoint(i*2+1,highPtInelAlice[i],inelPythia109[i+3]/inelAlice[i]);
   cout << "inelAlice[" << i << "]           = " << inelAlice[i] <<endl;
   cout << "centerPtInelAlice[" << i << "]   = " << centerPtInelAlice[i] <<endl;
   cout << "inelPythia109[" << i+3 << "]          = " << inelPythia109[i+3] <<endl;
   cout << "ptPythia109[" << i+3 << "]            = " << ptPythia109[i+3] <<endl;
}
//TGraphErrors *ratioPythia109Alice = divide(graphPythia109,graphInelAlice);
ratioPythia109Alice->SetFillColor(colorPythia109);
ratioPythia109Alice->SetLineColor(colorPythia109);
ratioPythia109Alice->SetLineWidth(2);
ratioPythia109Alice->SetLineStyle(7);
ratioPythia109Alice->SetMarkerStyle(26);
ratioPythia109Alice->SetFillStyle(3354);

TGraph *ratioPythia306Alice = new TGraph(binsInelAlice*2);
for (int i=0; i < binsInelAlice; i++) {
   ratioPythia306Alice->SetPoint(i*2,lowPtInelAlice[i],inelPythia306[i+3]/inelAlice[i]);
   ratioPythia306Alice->SetPoint(i*2+1,highPtInelAlice[i],inelPythia306[i+3]/inelAlice[i]);
   cout << "inelAlice[" << i << "]           = " << inelAlice[i] <<endl;
   cout << "centerPtInelAlice[" << i << "]   = " << centerPtInelAlice[i] <<endl;
   cout << "inelPythia306[" << i+3 << "]          = " << inelPythia306[i+3] <<endl;
   cout << "ptPythia306[" << i+3 << "]            = " << ptPythia306[i+3] <<endl;
}
//TGraphErrors *ratioPythia306Alice = divide(graphPythia306,graphInelAlice);
ratioPythia306Alice->SetFillColor(colorPythia306);
ratioPythia306Alice->SetLineColor(colorPythia306);
ratioPythia306Alice->SetLineWidth(2);
ratioPythia306Alice->SetLineStyle(3);
ratioPythia306Alice->SetMarkerStyle(26);
ratioPythia306Alice->SetFillStyle(3354);

TGraph *ratioPythia320Alice = new TGraph(binsInelAlice*2);
for (int i=0; i < binsInelAlice; i++) {
   ratioPythia320Alice->SetPoint(i*2,lowPtInelAlice[i],inelPythia320[i+3]/inelAlice[i]);
   ratioPythia320Alice->SetPoint(i*2+1,highPtInelAlice[i],inelPythia320[i+3]/inelAlice[i]);
   cout << "inelAlice[" << i << "]           = " << inelAlice[i] <<endl;
   cout << "centerPtInelAlice[" << i << "]   = " << centerPtInelAlice[i] <<endl;
   cout << "inelPythia320[" << i+3 << "]          = " << inelPythia320[i+3] <<endl;
   cout << "ptPythia320[" << i+3 << "]            = " << ptPythia320[i+3] <<endl;
}
//TGraphErrors *ratioPythia320Alice = divide(graphPythia320,graphInelAlice);
ratioPythia320Alice->SetFillColor(colorPythia320);
ratioPythia320Alice->SetLineColor(colorPythia320);
ratioPythia320Alice->SetLineWidth(2);
ratioPythia320Alice->SetLineStyle(9);
ratioPythia320Alice->SetMarkerStyle(26);
ratioPythia320Alice->SetFillStyle(3354);


TCanvas *can5 = new TCanvas("can5","CompInel",520,700);

TPad *pad5_1 = new TPad("pad5_1","pad5_1",0.0,0.35,1.0,1.0);
setAttrib(pad5_1);

TPad *pad5_2 = new TPad("pad5_2","pad5_2",0.0,0.0,1.0,0.35);
setAttrib(pad5_2);


//
can5->cd();

pad5_1->Draw();
pad5_1->cd();
pad5_1->SetLogx();
pad5_1->SetLogy();
//all->Draw("AE4");

TMultiGraph *all = new TMultiGraph("all","");

all->Add(graphPhojet);

all->Add(graphPythia109);
all->Add(graphPythia306);
all->Add(graphPythia320);
all->Draw("AXL");
all->GetXaxis()->SetLimits(minPt,maxPt);
all->GetXaxis()->SetTitle("p_{T} (GeV/c)");
all->GetXaxis()->SetTitleOffset(1.6);
all->GetYaxis()->SetTitleOffset(1.6);
all->GetYaxis()->SetTitle("1/N_{evt} 1/(2#pi p_{T}) (d^{2}N_{ch})/(d#eta dp_{T}) (GeV/c)^{-2}");
all->GetXaxis()->SetLimits(minPt,maxPt);
all->GetYaxis()->SetRangeUser(5e-8,10);
all->GetYaxis()->SetLimits(5e-8,10);
setAttrib(all);
graphInelAliceFit->Draw("PZ");

TLegend *l5_1 = new TLegend(0.200,0.03,0.65,0.376);
l5_1->SetHeader("pp, INEL, #sqrt{s} = 900 GeV, | #eta | < 0.8");
l5_1->AddEntry(graphInelAliceFit,"ALICE data","PL");
l5_1->AddEntry(graphPhojet,"PHOJET","L");
l5_1->AddEntry(graphPythia109,"PYTHIA D6T (109)","L");
l5_1->AddEntry(graphPythia306,"PYTHIA ATLAS-CSC (306)","L");
l5_1->AddEntry(graphPythia320,"PYTHIA Perugia0 (320)","L");

//l1->AddEntry(graphCms,"CMS |#eta|<2.4","LP");
//la->AddEntry(graphUA1,"UA1","p");
//la->AddEntry(graphCMS,"CMS","p");
//la->AddEntry(graphCMSred,"CMS, #eta<0.8","p");
//la->AddEntry(graphATLAS,"ATLAS","p");
l5_1->SetFillColor(0);
l5_1->SetLineColor(0);
l5_1->SetTextSize(legendTextSize);
l5_1->Draw();


//histo1->GetXaxis()->SetNdivisions(405);

can5->cd();
pad5_2->Draw();
pad5_2->cd();
pad5_2->SetLogx();

TMultiGraph *ratios = new TMultiGraph("ratios","");


//ratios->Add(ratioAliceAlice);
//ratios->Add(ratioAtlasAlice);
ratios->Add(ratioPhojetAlice);
ratios->Add(ratioPythia109Alice);
ratios->Add(ratioPythia306Alice);
ratios->Add(ratioPythia320Alice);

ratioAliceAlice->Draw("AE3");
ratioAliceAlice->SetTitle("");
ratioAliceAlice->GetXaxis()->SetLimits(minPt,maxPt);
ratioAliceAlice->GetXaxis()->SetTitle("p_{T} (GeV/c)");
ratioAliceAlice->GetXaxis()->SetTitleOffset(1.4);
ratioAliceAlice->GetYaxis()->SetTitleOffset(0.9);
ratioAliceAlice->GetYaxis()->SetTitle("ratio");
ratioAliceAlice->GetYaxis()->SetLabelSize(0.07);
ratioAliceAlice->GetXaxis()->SetLabelSize(0.07);
ratioAliceAlice->GetXaxis()->SetTitleSize(0.07);
ratioAliceAlice->GetYaxis()->SetTitleSize(0.07);
ratioAliceAlice->GetYaxis()->SetRangeUser(0.3,1.7);
ratioAliceAlice->GetYaxis()->SetLimits(0.3,1.7);
ratioAliceAlice->GetYaxis()->CenterTitle();
setAttrib(ratioAliceAlice);
ratioAliceAlice->GetYaxis()->SetLabelOffset(0.023); 


//ratio->SetTitle("ratio ALICE/ATLAS (different #eta)");
//ratios->Draw("XL");
ratios->Draw("XL");
//ratios->GetYaxis()->SetRangeUser(0.41,1.49);
//ratios->GetYaxis()->SetLimits(0.41,1.49);


graphInelAlice->SetFillColor(2);

TGraph* tmpG = new TGraph();
tmpG->SetMarkerColor(colorAliceFit);
tmpG->SetLineColor(colorAliceFit);
tmpG->SetMarkerStyle(10);
tmpG->SetLineWidth(2);
tmpG->SetLineStyle(1);

TLegend *l5_2 = new TLegend(0.200,0.759+0.005,0.65,0.973);
l5_2->AddEntry(ratioAliceAlice,"ALICE data uncertainties","F");
l5_2->AddEntry(tmpG,"MC / data","L");
//l2->AddEntry(ratioCmsAlice,"CMS / ALICE","F");
l5_2->SetFillColor(0);
l5_2->SetLineColor(0);
l5_2->SetTextSize(legendTextSize);
l5_2->Draw();



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

logoPrelim(can5);
}
