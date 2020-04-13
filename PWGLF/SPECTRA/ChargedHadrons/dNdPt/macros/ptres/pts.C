#include "TMath.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1D.h"
#include "TRandom3.h"
#include "TPad.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TObjString.h"


void pts(Double_t power = 0., const char* infile = "scaledcov.root", const char* outfile = "ptrescorr.root", Double_t ntracks=1e8, Double_t resfactor = 1.0, Double_t minpt=7., Double_t extrasmear = 0, Double_t meanshifta = 0., Double_t sigmashifta = 0., Double_t meanshiftc = 0., Double_t sigmashiftc = 0., TF1* ptdist = 0, const char* tag="") 
{

if (!ptdist) {
    ptdist = new TF1("ptdist","x^[0]",minpt,100);
    ptdist->SetParameter(0,-1.*TMath::Abs(power));
} else {
    power=0;
}

Bool_t shifts = kFALSE;
if (meanshifta != 0) shifts = kTRUE;
if (sigmashifta != 0) shifts = kTRUE;
if (meanshiftc != 0) shifts = kTRUE;
if (sigmashiftc != 0) shifts = kTRUE;
    
delete gRandom;
gRandom = new TRandom3(0);
gStyle->SetOptStat(kFALSE);

    TDatime d;
    TString ids = TString("file generated ");
    ids += d.AsString();
    ids += " inputfile=";
    ids += infile;
    ids += " outputfile=";
    ids += outfile;
    ids += " power=";
    ids += power;
    ids += " ntracks=";
    ids += ntracks;
    ids += " resfactor=";
    ids += resfactor;
    ids += " minpt=";
    ids += minpt;
    ids += " extrasmear=";
    ids += extrasmear;    
    if (shifts) {
        ids += " meanshifta=";
        ids += meanshifta;
        ids += " sigmashifta=";
        ids += sigmashifta;
        ids += " meanshiftc=";
        ids += meanshiftc;
        ids += " sigmashiftc=";
        ids += sigmashiftc;
    }    
    ids += " tag=";
    ids += tag;

    TObjString oids(ids);    

TFile* fc = TFile::Open(infile,"READ");
TH1F* sc = (TH1F*) fc->Get("scaledcov");
TF1*  nr = (TF1*)  fc->Get("normres");
const Int_t nbins = 41;
Double_t xbins[42] = {0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.6, 4.0, 5.0, 6.0, 8.0,  10.0, 13.0, 20.0, 30.0, 50.0, 80.0, 100.0};

TH1D* ptgen = new TH1D("ptgen","ptgen",41,xbins);
TH1D* ptrec = new TH1D("ptrec","ptrec",41,xbins);
TH1D* ptgenb = new TH1D("ptgenb","ptgenb",100,0,100);
TH1D* ptrecb = new TH1D("ptrecb","ptrecb",100,0,100);


ptgen->Sumw2();
ptrec->Sumw2();
ptgenb->Sumw2();
ptrecb->Sumw2();

for (ULong64_t i=0; i<ntracks; i++) {
    if (i%1000000 == 0) cout<<i/1000000<<endl;
    Double_t pt  = ptdist->GetRandom();
    Double_t pt1 = 1./pt; 
    Double_t sspt1 = sc->GetRandom();    
    Double_t spt1 = sspt1 + nr->Eval(pt1);
    spt1 *= resfactor;
    if (extrasmear!=0) spt1 = TMath::Sqrt(spt1*spt1+extrasmear*extrasmear);
    if (shifts) {
        int side = gRandom->Integer(2);
        int charge = gRandom->Integer(2);
        Double_t dpt1 = 0;        
        if (side)    { dpt1 = gRandom->Gaus(meanshifta,sigmashifta); }
                else { dpt1 = gRandom->Gaus(meanshiftc,sigmashiftc); }            
        if (charge) { pt1 += dpt1; }
               else { pt1 -= dpt1; }
    }
    Double_t pt1s = gRandom->Gaus(pt1,spt1);
    Double_t pts = 1./pt1s;
    ptgen->Fill(pt);
    ptgenb->Fill(pt);
    ptrec->Fill(pts);
    ptrecb->Fill(pts);
}

    TCanvas* c1 = new TCanvas();
    c1->cd();

ptgen->SetLineWidth(2);
ptgen->SetLineColor(kBlack);
ptgenb->SetLineColor(kRed);

TH1D* corr = (TH1D*) ptgen->Clone("corr");
TH1D* corrb = (TH1D*) ptgenb->Clone("corrb");

corr->Divide(ptgen,ptrec,1,1,"");
corrb->Divide(ptgenb,ptrecb,1,1,"");
corrb->Draw("");
corr->Draw("SAME h");

corrb->GetXaxis()->SetRangeUser(minpt,50);
corrb->GetYaxis()->SetRangeUser(0.9,1.);

gPad->SetGridy();

TFile* fout = TFile::Open(outfile,"RECREATE");
ptgen->Write("ptgen");
ptgenb->Write("ptgenb");
ptrec->Write("ptrec");
ptrecb->Write("ptrecb");
corr->Write("corr");
corrb->Write("corrb");

c1->Write("c1_CORRFACTOR");

sc->Write("scaledcov");
nr->Write("normres");
ptdist->Write("ptdist");
oids.Write("INFO");

fout->Close();
cout<<ids<<endl;

}
