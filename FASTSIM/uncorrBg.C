void uncorrBg(Int_t nev = 1000000)
{
// b: 0.00261     (0.0375 per unit of rapidity per event)
// c: 1.029       (1.35 per unit of rapidity per event)    
//
// rate of central collisions  400 Hz
// central events 4 10^8
// Get the single muon pT spectra
//
// 
    Float_t scaleC = 200.;
    Float_t scaleB = 200.;
    Float_t scaleD = 0.34 * 1.362e-3 * 1.e3;
 
//
    Float_t bbx[15] = 
	{
	    0.21, 0.55, 0.65, 0.65,  0.5, 
	    0.40, 0.26, 0.19, 0.11,  0.1, 
	    0.075, 0.045, 0.035, 0.02, 0.017 
	};
//
//  Fast response
//    
    AliFastMuonTriggerEff *trigeff = new AliFastMuonTriggerEff();
    AliFastMuonTrackingAcc *acc = new AliFastMuonTrackingAcc();
    AliFastMuonTrackingEff *eff = new AliFastMuonTrackingEff();
    AliFastMuonTrackingRes *res = new AliFastMuonTrackingRes();
    acc->SetBackground(1.);
    eff->SetBackground(1.);
    res->SetBackground(1.);  
    acc->Init(); 
    eff->Init(); 
    res->Init(); 
    AliFastDetector* tracker = new AliFastDetector("Tracker", "Muon Tracker");
    tracker->AddResponse(trigeff);
    tracker->AddResponse(acc);
    tracker->AddResponse(eff);
    tracker->AddResponse(res);
    tracker->Init();
    tracker->Dump();

	    
//
//  Heavy Flavors
//
    f = new TFile("HVQinc.root");
    TH1F* ptBB = (TH1F*) f->Get("hPtCorra");
    TH1F* ptCC = (TH1F*) f->Get("hpta");   
    TCanvas *c5 = new TCanvas("c5","Canvas 6",400,10,600,700);

    TF1*  ptBBLf = new TF1("ptBBLf", "[0] * x / (1. + (x/[1])**2)**[2]", 0., 3.);
    ptBBLf->SetParameter(0, 4.46695e-03);
    ptBBLf->SetParameter(1, 1.60242e+00);
    ptBBLf->SetParameter(2, 2.24948e+00);

    TF1*  ptBBHf = new TF1("ptBBHf", "[0] * x / (1. + (x/[1])**2)**[2]", 3., 20.);
    ptBBHf->SetParameter(0, 2.59961e-03);
    ptBBHf->SetParameter(1, 2.41);
    ptBBHf->SetParameter(2, 3.075);

    TF1*  ptCCHf = new TF1("ptCCHf", "[0] * x / (1. + (x/[1])**2)**([2] + [3] * x)", 1.5, 20.);
    ptCCHf->SetParameter(0, 6.72360e-01);
    ptCCHf->SetParameter(1, 7.06043e-01);
    ptCCHf->SetParameter(2, 2.74240e+00);
    ptCCHf->SetParameter(3, 8.45018e-03);
//    ptCCHf->Draw("ac");

    TF1*  ptCCLf = new TF1("ptCCLf", "[0] * x / (1. + (x/[1])**2)**([2] + [3] * x)", 0., 1.5);
    ptCCLf->SetParameter(0, 1.40260e+00);
    ptCCLf->SetParameter(1, 3.75762e-01);
    ptCCLf->SetParameter(2, 1.54345e+00);
    ptCCLf->SetParameter(3, 2.49806e-01);
//    ptCCLf->Draw("ac");
    /*    
    
    TF1*  ptCCHf = new TF1("ptCCHf", "[0] * x / (1. + (x/[1])**2)**[2]", 0., 20.);
    ptCCHf->SetParameter(0, 0.249);
    ptCCHf->SetParameter(1, 1.15);
    ptCCHf->SetParameter(2, 3.33);
//    ptCCHf->Draw("ac");

    TF1*  ptCCLf = new TF1("ptCCLf", "[0] * x / (1. + (x/[1])**2)**[2]", 0., 20.);
    ptCCLf->SetParameter(0, 1.125);
    ptCCLf->SetParameter(1, 0.525);
    ptCCLf->SetParameter(2, 2.42);
//    ptCCLf->Draw("ac");
*/

    
    TF1*  ptBf = new TF1("ptBf", "[0] * x / (1. + (x/[1])**2)**[2]", 0., 20.);
    ptBf->SetParameter(0, 1.e5 * 0.7 * 1.125);
    ptBf->SetParameter(1, 6.27);
    ptBf->SetParameter(2, 3.2);
//    ptBf->Draw("ac");
//
//  pi/K -> mu
//
    f->Close();
    f = new TFile("pikmu.root");
    TH2F* etaptPiK = (TH2F*) f->Get("etaptH");
    TAxis* etaAxis = etaptPiK->GetXaxis();
    TAxis* ptAxis  = etaptPiK->GetYaxis();    
    TH1F* ptPiK    = (TH1F*) f->Get("ptH3");
//    ptAxis = ptPiK->GetXaxis();
    
//
// Book histograms
    TH1F* massBBH = new TH1F("massBBH", "Mass Spectrum: b-b        ", 150, 0., 15.);
    TH1F* massCCH = new TH1F("massCCH", "Mass Spectrum: c-c        ", 150, 0., 15.);
    TH1F* massBCH = new TH1F("massBCH", "Mass Spectrum: b-c        ", 150, 0., 15.);
    TH1F* massDDH = new TH1F("massDDH", "Mass Spectrum: decay-decay", 150, 0., 15.);
    TH1F* massBDH = new TH1F("massBDH", "Mass Spectrum: decay-b    ", 150, 0., 15.);
    TH1F* massCDH = new TH1F("massCDH", "Mass Spectrum: decay-c    ", 150, 0., 15.);    
    TH1F* ptCH    = new TH1F("ptCH", "pt Spectrum mu from c", 20., 0., 10.);    
    TH1F* ptBH    = new TH1F("ptBH", "pt Spectrum mu from b", 20., 0., 10.);    
    TH1F* ptDH    = new TH1F("ptDH", "pt Spectrum mu from pi/K", 20., 0., 10.);    
    TH1F* ptBH2    = new TH1F("ptBH2", "pt Spectrum mu from b", 20., 0., 10.);    
    TH1F* ptBH3    = new TH1F("ptBH3", "pt Spectrum mu from b", 15., 0., 15.);
    for (Int_t i = 0; i < 15; i++)
    {
	ptBH3->SetBinContent(i+1, bbx[i]);
    }
    ptBH3->Draw();
    
    
//
// Event Loop
//
    Int_t iev;
    for (iev = 0; iev < nev; iev++) {
//
// pT
	Float_t pT1 = 20. * gRandom->Rndm();
	Float_t pT2 = 20. * gRandom->Rndm();
//
// phi
	Float_t phi1 = 2. * TMath::Pi() * gRandom->Rndm() - TMath::Pi();
	Float_t phi2 = 2. * TMath::Pi() * gRandom->Rndm() - TMath::Pi();
	Float_t dphi = phi1 - phi2;
//
// eta
	Float_t eta1 = 1.457 * gRandom->Rndm() + 2.543;
	Float_t eta2 = 1.457 * gRandom->Rndm() + 2.543;	
	Float_t deta = eta1 - eta2;
//
// invariant mass
	Float_t m2 = 2. * pT1 * pT2 * (TMath::CosH(deta) - TMath::Cos(dphi));
	Float_t m  = TMath::Sqrt(m2);

//
// Smearing
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

	/*
	ptBin  = ptAxis ->FindBin(pT1);	
	wgtD1  = ptPiK->GetBinContent(ptBin) * scaleD;
	ptBin  = ptAxis ->FindBin(pT2);	
	wgtD2  = ptPiK->GetBinContent(ptBin) * scaleD;
	*/

//
//   efficiencies
//	
	Float_t theta1 = 2. * TMath::ATan(TMath::Exp(-eta1)) * 180./TMath::Pi();
	Float_t theta2 = 2. * TMath::ATan(TMath::Exp(-eta2)) * 180./TMath::Pi();
	Float_t phid1  = phi1 * 180./TMath::Pi();
	Float_t phid2  = phi2 * 180./TMath::Pi();

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
	Float_t tri2  = trigeff->Evaluate(1, pT2, theta2, phid2);

	Float_t effA   = eff1 * eff2 * acc1 * acc2 * tri1 * tri2;
	
	Float_t ptMax = pT1;
	Float_t ptMin = pT2;
	if (pT2 > pT1) {
	    Float_t ptMax = pT2;
	    Float_t ptMin = pT1;
	}
	

//	if (
//	    (ptMax > 6. && ptMin > 3.) ||
//	    (ptMax < 6. && ptMin > (6. - 0.5 * ptMax))
//	    ) 
//	{
	if (ptMin > 3.) {
	    massBBH->Fill(m, wgtB1 * wgtB2 / 4. * effA);
	    massCCH->Fill(m, wgtC1 * wgtC2 / 4. * effA);
	    massBCH->Fill(m, wgtC1 * wgtB2 / 4. * effA);
	    massBCH->Fill(m, wgtC2 * wgtB1 / 4. * effA);
	    massDDH->Fill(m, wgtD1 * wgtD2 / 4. * effA);
	    massBDH->Fill(m, wgtB1 * wgtD2 / 4. * effA);
	    massBDH->Fill(m, wgtB2 * wgtD1 / 4. * effA);
	    massCDH->Fill(m, wgtC1 * wgtD2 / 4. * effA);
	    massCDH->Fill(m, wgtC2 * wgtD1 / 4. * effA);
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
	}

    } // event loop
//    Float_t eff    = 0.6 * 1.0;
    Float_t evtWgt = 1. / Float_t(nev) * 4.e8;
    
    
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
 
}
