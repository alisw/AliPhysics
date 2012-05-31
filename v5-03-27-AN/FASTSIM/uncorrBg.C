#if !defined(__CINT__) || defined(__MAKECINT__)
#include <stdio.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TF1.h>
#include <TDatabasePDG.h>
#include <TRandom.h>
#include <TStopwatch.h>
#include "AliFastGlauber.h"
#include "AliFastMuonTrackingRes.h"
#include "AliFastMuonTrackingAcc.h"
#include "AliFastMuonTrackingEff.h"
#include "AliFastMuonTriggerEff.h"
#endif

void uncorrBg(Int_t nev = 1000000, Double_t bmin = 0., Double_t bmax = 5.)
{
//
//
//
    AliFastGlauber*  glauber = new AliFastGlauber();
    glauber->Init(1);
    
    Float_t lumi   = 5.e26;   // cm^-2 s^-1
    Float_t time   = 1.e6;    // s
    Float_t rate   = lumi/1.e27 * glauber->CrossSection(bmin, bmax); // Hz
    Float_t fhard  = glauber->FractionOfHardCrossSection(bmin, bmax);
    Float_t fgeo   = glauber->CrossSection(bmin, bmax) / glauber->CrossSection(0, 100.);
    Float_t events = rate * time;

    printf("-------------------------------------------------------------------\n");
    printf("Impact parameter range: %10.3f - %10.3f fm \n", bmin, bmax);
    printf("Luminosity: %10.3e cm^-2 s^-1 \n", lumi);
    printf("Rate: %10.3f Hz\n", rate);
    printf("Fraction of hard  cross-section: %10.3f \n", fhard * 100.);
    printf("Fraction of geom. cross-section: %10.3f \n", fgeo  * 100.);
    printf("Events in 10^6 s: %10.3e \n", events);
    printf("-------------------------------------------------------------------\n");

//
//    
    Float_t ptMinCut  = 3.;     // GeV
    Float_t etamin    = 2.543;
    Float_t etar      = 1.457;
    Float_t ptUp      = 20.;    // GeV
    Float_t dpt       = 0.01;   // GeV
//    
//  For b = 0
//  (factor 1.28 to scale from 10% most central to b=0)
//
    Float_t aGlauber = 208. * 208. * 7.03e-4;    // 1/mb
    
    Float_t scaleC0 = aGlauber * ptUp / dpt;
    Float_t scaleB0 = aGlauber * ptUp / dpt;
    Float_t scaleD0 = 1.28 * etar * ptUp / 1.35; // scaled by 1.35 to match ALICE-INT-2002-6
    
//
//
//  Fast response
//    
    AliFastMuonTriggerEff  *trigeff = new AliFastMuonTriggerEff();
    AliFastMuonTrackingAcc *acc     = new AliFastMuonTrackingAcc();
    AliFastMuonTrackingEff *eff     = new AliFastMuonTrackingEff();
    AliFastMuonTrackingRes *res     = new AliFastMuonTrackingRes();
    acc->SetBackground(1.);
    eff->SetBackground(1.);
    res->SetBackground(1.);  

    acc    ->Init(); 
    eff    ->Init(); 
    res    ->Init(); 
    trigeff->Init();
    
/*
//  To be improved
//
    AliFastDetector* tracker = new AliFastDetector("Tracker", "Muon Tracker");
    tracker->AddResponse(acc);
    tracker->AddResponse(eff);
    tracker->AddResponse(res);

    AliFastDetector* trigger = new AliFastDetector("Trigger", "Muon Trigger");
    trigger->AddResponse(trigeff);

    AliFastDetector* spectro = new AliFastDetector("Spectro", "Muon Spectrometer");
    spectro->AddSubdetector(tracker, "");
    spectro->AddSubdetector(trigger, "");    
    spectro->Init();
*/
	    
//
//  Heavy Flavors
//
    
    TF1*  ptBBLf = new TF1("ptBBLf", "[0] * x / (1. + (x/[1])**2)**[2]", 0., 3.);
    ptBBLf->SetParameter(0, 4.390e-05);
    ptBBLf->SetParameter(1, 1.8706);
    ptBBLf->SetParameter(2, 2.6623);

    TF1*  ptBBHf = new TF1("ptBBHf", "[0] * x / (1. + (x/[1])**2)**[2]", 3., ptUp);
    ptBBHf->SetParameter(0, 2.5329e-05);
    ptBBHf->SetParameter(1, 2.6067);
    ptBBHf->SetParameter(2, 3.3821);

    TF1*  ptCCHf = new TF1("ptCCHf", "[0] * x / (1. + (x/[1])**2)**([2] + [3] * x)", 1.5, ptUp);
    ptCCHf->SetParameter(0, 4.8234e-03);
    ptCCHf->SetParameter(1, 7.5656e-01);
    ptCCHf->SetParameter(2, 2.7707e+00);
    ptCCHf->SetParameter(3, 2.3362e-02);

    TF1*  ptCCLf = new TF1("ptCCLf", "[0] * x / (1. + (x/[1])**2)**([2] + [3] * x)", 0., 1.5);
    ptCCLf->SetParameter(0, 1.190e-02);
    ptCCLf->SetParameter(1, 3.6343e-01);
    ptCCLf->SetParameter(2, 1.4689e+00);
    ptCCLf->SetParameter(3, 2.5772e-01);
    
    TF1*  ptBf = new TF1("ptBf", "[0] * x / (1. + (x/[1])**2)**[2]", 0., ptUp);
    ptBf->SetParameter(0, 1.e5 * 0.7 * 1.125);
    ptBf->SetParameter(1, 6.27);
    ptBf->SetParameter(2, 3.2);
//
//  pi/K -> mu
//
    TFile* f = new TFile("$(ALICE_ROOT)/FASTSIM/data/pikmu.root");
    TH2F*  etaptPiK = (TH2F*) f->Get("etaptH");
    TAxis* etaAxis  = etaptPiK->GetXaxis();
    TAxis* ptAxis   = etaptPiK->GetYaxis();    
//
//
// Book histograms
    TH1F* massBBH = new TH1F("massBBH", "Mass Spectrum: b-b        ", 150, 0., 15.);
    TH1F* massCCH = new TH1F("massCCH", "Mass Spectrum: c-c        ", 150, 0., 15.);
    TH1F* massBCH = new TH1F("massBCH", "Mass Spectrum: b-c        ", 150, 0., 15.);
    TH1F* massDDH = new TH1F("massDDH", "Mass Spectrum: decay-decay", 150, 0., 15.);
    TH1F* massBDH = new TH1F("massBDH", "Mass Spectrum: decay-b    ", 150, 0., 15.);
    TH1F* massCDH = new TH1F("massCDH", "Mass Spectrum: decay-c    ", 150, 0., 15.);    
    TH1F* ptCH    = new TH1F("ptCH", "pt Spectrum mu from c",          20, 0., 10.);    
    TH1F* ptBH    = new TH1F("ptBH", "pt Spectrum mu from b",          20, 0., 10.);    
    TH1F* ptDH    = new TH1F("ptDH", "pt Spectrum mu from pi/K",       20, 0., 10.);    
    TH1F* ptBH2   = new TH1F("ptBH2", "pt Spectrum mu from b",         20, 0., 10.);    
    TH1F* costBBH = new TH1F("costBBH", "cos theta*    ", 20, -1, 1.);    
//
// Event Loop
//
    Int_t iev;
    for (iev = 0; iev < nev; iev++) {
//
//  Collision geometry
//
	Float_t b;
	b = glauber->GetRandomImpactParameter(bmin, bmax);
	Double_t nbinary = glauber->Binaries(b);
	Float_t  scaleC  = scaleC0 * nbinary;
	Float_t  scaleB  = scaleB0 * nbinary;
	Float_t  scaleD  = scaleD0 * nbinary;
//
// pT
	Float_t pT1 = ptUp * gRandom->Rndm();
	Float_t pT2 = ptUp * gRandom->Rndm();
//
// phi
	Float_t phi1 = 2. * TMath::Pi() * gRandom->Rndm() - TMath::Pi();
	Float_t phi2 = 2. * TMath::Pi() * gRandom->Rndm() - TMath::Pi();
	Float_t dphi = phi1 - phi2;
//
// eta
	Float_t eta1 = etar * gRandom->Rndm() + etamin;
	Float_t eta2 = etar * gRandom->Rndm() + etamin;	
	Float_t deta = eta1 - eta2;
//
// invariant mass
	Float_t m2   = 2. * pT1 * pT2 * (TMath::CosH(deta) - TMath::Cos(dphi));
	Float_t m    = TMath::Sqrt(m2);
	Float_t mt2  = pT1 * pT1 + pT2 * pT2 + 2. * pT1 * pT2 * TMath::CosH(deta);
	Float_t mt   = TMath::Sqrt(mt2);
	Float_t cost = 2. * pT1 * pT2 * TMath::SinH(deta) / m / mt;  

//
// Smearing (to be improved)
	Float_t dm = m * 0.01;
	m += gRandom->Gaus(0., dm);	
//
// Weights
//
//      Heavy Flavour
//
	Int_t ibin;
	Float_t wgtB1, wgtB2;
	Float_t wgtC1, wgtC2;

	if (pT1 > 1.5) {
	    wgtC1 = ptCCHf->Eval(pT1) * scaleC;
	} else {
	    wgtC1 = ptCCLf->Eval(pT1) * scaleC;
	}
	if (pT2 > 1.5) {
	    wgtC2 = ptCCHf->Eval(pT2) * scaleC;
	} else {
	    wgtC2 = ptCCLf->Eval(pT2) * scaleC;
	}


	if (pT1 > 3.) {
	    wgtB1 = ptBBHf->Eval(pT1) * scaleB;
	} else {
	    wgtB1 = ptBBLf->Eval(pT1) * scaleB;
	}
	if (pT2 > 3.) {
	    wgtB2 = ptBBHf->Eval(pT2) * scaleB;
	} else {
	    wgtB2 = ptBBLf->Eval(pT2) * scaleB;
	}


//
//      Weight  for decays
//
	Int_t etaBin, ptBin;
	Float_t wgtD1, wgtD2;
	
	etaBin = etaAxis->FindBin(eta1);
	ptBin  = ptAxis ->FindBin(pT1);	
	wgtD1  = etaptPiK->GetBinContent(etaBin, ptBin) * scaleD;
	
	etaBin = etaAxis->FindBin(eta2);
	ptBin  = ptAxis ->FindBin(pT2);	
	wgtD2  = etaptPiK->GetBinContent(etaBin, ptBin) * scaleD;
	
//
//      Efficiencies
//	
	Float_t theta1 = 2. * TMath::ATan(TMath::Exp(-eta1)) * 180./TMath::Pi();
	Float_t theta2 = 2. * TMath::ATan(TMath::Exp(-eta2)) * 180./TMath::Pi();
	Float_t phid1  = phi1 * 180./TMath::Pi();
	Float_t phid2  = phi2 * 180./TMath::Pi();
	Float_t p1     = pT1/TMath::Sin(theta1 * TMath::Pi()/180.);
	Float_t p2     = pT2/TMath::Sin(theta2 * TMath::Pi()/180.);
	
	res->SetCharge(1);
	eff->SetCharge(1);
	acc->SetCharge(1);
	Float_t eff1  = eff->Evaluate(pT1, theta1, phid1);
	Float_t acc1  = acc->Evaluate(pT1, theta1, phid1);
	Float_t tri1  = trigeff->Evaluate(1, pT1, theta1, phid1);
	res->SetCharge(-1);
	eff->SetCharge(-1);
	acc->SetCharge(-1);
	Float_t eff2  = eff->Evaluate(pT2, theta2, phid2);
	Float_t acc2  = acc->Evaluate(pT2, theta2, phid2);
	Float_t tri2  = trigeff->Evaluate(-1, pT2, theta2, phid2);

	Float_t effA   = eff1 * eff2 * acc1 * acc2 * tri1 * tri2;

	Float_t ptMax = pT1;
	Float_t ptMin = pT2;
	if (pT2 > pT1) {
	    ptMax = pT2;
	    ptMin = pT1;
	}
	
	if (ptMin > ptMinCut && p1 > 4. && p2 > 4.) {
	    massBBH->Fill(m, wgtB1 * wgtB2 / 4. * effA);
	    massCCH->Fill(m, wgtC1 * wgtC2 / 4. * effA);
	    massBCH->Fill(m, wgtC1 * wgtB2 / 4. * effA);
	    massBCH->Fill(m, wgtC2 * wgtB1 / 4. * effA);
	    massDDH->Fill(m, wgtD1 * wgtD2 / 4. * effA);
	    massBDH->Fill(m, wgtB1 * wgtD2 / 4. * effA);
	    massBDH->Fill(m, wgtB2 * wgtD1 / 4. * effA);
	    massCDH->Fill(m, wgtC1 * wgtD2 / 4. * effA);
	    massCDH->Fill(m, wgtC2 * wgtD1 / 4. * effA);

	    costBBH->Fill(cost, wgtB1 * wgtB2 / 4. * effA);
	}
	//
	// pT - Spectra
	//
	for (Int_t ipt = 0; ipt < 20; ipt++)
	{
	    Float_t pt = 0.5 * ipt;
	    ptBH2->Fill(pT1, wgtB1);		
	    if (pT1 > pt) {
		ptCH->Fill(pt, wgtC1);
		ptBH->Fill(pt, wgtB1);
		ptDH->Fill(pt, wgtD1);
	    }
	} // bins
    } // event loop

    Float_t evtWgt = events / Float_t(nev);
    
    massBBH->Scale(evtWgt);
    massCCH->Scale(evtWgt);
    massBCH->Scale(evtWgt);
    massDDH->Scale(evtWgt);
    massBDH->Scale(evtWgt);
    massCDH->Scale(evtWgt);
    
    TH1F * massALH = new TH1F(*massCDH);
    massALH->Add(massBDH);
    massALH->Add(massDDH);
    massALH->Add(massBCH);    
    massALH->Add(massCCH);
    massALH->Add(massBBH);     

    TCanvas *c0 = new TCanvas("c0","Canvas 1",400,10,600,700);
    massCCH->SetLineColor(4);
    massCCH->SetMinimum(1.);
    massCCH->SetMaximum(1.e4);
    massCCH->SetXTitle("m_{#mu#mu} [GeV]");
    massCCH->SetYTitle("Entries/100 MeV /10^{6} s");
    massCCH->Draw("");
    massALH->SetLineColor(3);
    massALH->Draw("same");
    massBBH->SetLineColor(6);
    massBBH->Draw("same");

    TCanvas *c2 = new TCanvas("c2","Canvas 3",400,10,600,700);
    massDDH->SetLineColor(2);
    massDDH->SetMinimum(1.e2);
    massDDH->SetMaximum(1.e6);
    massDDH->SetXTitle("m_{#mu#mu} [GeV]");
    massDDH->SetYTitle("Entries/100 MeV /10^{6} s");
    massDDH->Draw("");
    massALH->SetLineColor(3);
    massALH->Draw("same");
    massCCH->SetLineColor(4);
    massCCH->Draw("same");
    massBBH->SetLineColor(6);
    massBBH->Draw("same");

    TCanvas *c3 = new TCanvas("c3","Canvas 4",400,10,600,700);
    ptCH->Scale(1./float(nev));
    ptBH->Scale(1./float(nev));    
    ptDH->Scale(1./float(nev));    
    ptCH->SetLineColor(4);
    ptBH->SetLineColor(6);
    ptDH->SetLineColor(2);
    ptCH->SetXTitle("p_{T}^{min} [GeV]");
    ptCH->SetYTitle("<n>_{#mu}/event");
    
    ptDH->Draw();
    ptBH->Draw("same");
    ptCH->Draw("same");
    TCanvas *c4 = new TCanvas("c4","Canvas 5",400,10,600,700);
    ptBH2->Draw();

    TCanvas *c5 = new TCanvas("c5","Canvas 6",400,10,600,700);
    costBBH->Draw();
 
}
