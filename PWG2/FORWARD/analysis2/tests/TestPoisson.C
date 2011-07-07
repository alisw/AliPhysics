/** 
 * 
 * 
 * @param nBins 
 * @param min 
 * @param max 
 *
 * @ingroup pwg2_forward_scripts_tests
 */
void
MakeIntegerAxis(Int_t& nBins, Double_t& min, Double_t& max)
{
  if (nBins <= 0) return;
  if (max <= 0)   max = nBins;
  
  Double_t dBin = max / nBins;
  min   = - dBin / 2;
  max   = max + dBin / 2;
  nBins = nBins + 1;
}

/** 
 * 
 * 
 * @param o 
 * @param useWeights 
 *
 * @ingroup pwg2_forward_scripts_tests
 */
void
TestPoisson(Double_t o=.3, bool useWeights=false)
{
  const char* load = "$ALICE_ROOT/PWG2/FORWARD/analysis2/scripts/LoadLibs.C";
  gROOT->Macro(load);
  gROOT->GetInterpreter()->UnloadFile(gSystem->ExpandPathName(load));
  
  // --- Parameters of this script -----------------------------------
  Int_t      nBin =  5;  // Our detector matrix size 
  Int_t      nMax = TMath::Max(Int_t(nBin * nBin * o + .5)+nBin/2,nBin);  
  Int_t      nEv  = 10000; // Number of events
  Double_t   mp   = o;   // The 'hit' probability 


  TH2D* base = new TH2D("base", "Basic histogram", 
			nBin,-.5, nBin-.5, nBin, -.5, nBin-.5);
  base->SetXTitle("#eta");
  base->SetYTitle("#varphi");
  base->SetDirectory(0);
  base->SetOption("colz");

  Int_t tN1=nMax;    Double_t tMin1; Double_t tMax1;
  Int_t tN2=nMax*10; Double_t tMin2; Double_t tMax2=nMax;
  MakeIntegerAxis(tN1, tMin1, tMax1);
  MakeIntegerAxis(tN2, tMin2, tMax2);
  TH2D* corr = new TH2D("comp", "Comparison", 
			tN1, tMin1, tMax1, tN2, tMin2, tMax2);
  corr->SetXTitle("Input");
  corr->SetYTitle("Poisson");
  corr->SetDirectory(0);
  corr->SetOption("colz");
  corr->SetStats(0);
  TLine* lcorr = new TLine(0, 0, tMax2, tMax2);

  Int_t mm = TMath::Max(Int_t(nBin * o + .5),nBin/2);
  tN2=mm*10; tMax2 = mm;
  MakeIntegerAxis(tN2, tMin2, tMax2);
  Info("", "Making mean w/nbins=%d,range=[%f,%f]", tN2, tMin2, tMax2);
  TH2D* mean = new TH2D("mean", "Mean comparison", 
			tN2, tMin2, tMax2, tN2, tMin2, tMax2);
  mean->SetXTitle("Input");
  mean->SetYTitle("Poisson");
  mean->SetDirectory(0);
  mean->SetOption("colz");
  mean->SetStats(0);
  TLine* lmean = new TLine(tMin2, tMin2, tMax2, tMax2);

  TH1D* dist = new TH1D("dist", "Distribution of hits", tN1, tMin1, tMax1);
  dist->SetXTitle("s");
  dist->SetYTitle("P(s)");
  dist->SetFillColor(kRed+1);
  dist->SetFillStyle(3001);
  dist->SetDirectory(0);
			
  AliPoissonCalculator* c = new AliPoissonCalculator("ignored");
  c->Init(nBin ,nBin);

  for (Int_t i = 0; i < nEv; i++) { 
    c->Reset(base);
    base->Reset();

    for (Int_t iEta = 0; iEta < nBin; iEta++) { 
      for (Int_t iPhi = 0; iPhi < nBin; iPhi++) { 
	// Throw a die 
	Int_t m = gRandom->Poisson(mp);
	dist->Fill(m);

	// Fill into our base histogram 
	base->Fill(iEta, iPhi, m);

	// Fill into poisson calculator 
	c->Fill(iEta, iPhi, m > 0, (useWeights ? m : 1));
      }
    }
    // Calculate the result 
    TH2D* res = c->Result();
    
    // Now loop and compare 
    Double_t mBase = 0;
    Double_t mPois = 0;
    for (Int_t iEta = 0; iEta < nBin; iEta++) { 
      for (Int_t iPhi = 0; iPhi < nBin; iPhi++) { 
	Double_t p = res->GetBinContent(iEta, iPhi);
	Double_t t = base->GetBinContent(iEta, iPhi);

	mBase += t;
	mPois += p;
	corr->Fill(t, p);
      }
    }
    Int_t nn = nBin * nBin;
    mean->Fill(mBase / nn, mPois / nn);
  }

  TCanvas* cc = new TCanvas("c", "c", 900, 900);
  cc->SetFillColor(0);
  cc->SetFillStyle(0);
  cc->SetBorderMode(0);
  cc->SetRightMargin(0.02);
  cc->SetTopMargin(0.02);
  cc->Divide(2,2);
  
  TVirtualPad* pp = cc->cd(1);
  pp->SetFillColor(0);
  pp->SetFillStyle(0);
  pp->SetBorderMode(0);
  pp->SetRightMargin(0.15);
  pp->SetTopMargin(0.02);
  pp->SetLogz();
  pp->SetGridx();
  pp->SetGridy();
  corr->Draw();
  lcorr->Draw();

  pp = cc->cd(2);
  pp->SetFillColor(0);
  pp->SetFillStyle(0);
  pp->SetBorderMode(0);
  pp->SetRightMargin(0.02);
  pp->SetTopMargin(0.02);
#if 1
  c->GetMean()->Draw();
#elif 1
  c->GetOccupancy()->Draw();
#else
  pp->SetLogy();
  dist->SetStats(0);
  dist->Scale(1. / dist->Integral());
  dist->Draw();
  TH1D* m1 = c->GetMean();
  m1->Scale(1. / m1->Integral());
  m1->Draw("same");
  Double_t eI;
  Double_t ii = 100 * dist->Integral(2, 0);
  TLatex* ll = new TLatex(.97, .85, 
			  Form("Input #bar{m}: %5.3f", mp));
  ll->SetNDC();
  ll->SetTextFont(132);
  ll->SetTextAlign(31);
  ll->Draw();
  ll->DrawLatex(.97, .75, Form("Result #bar{m}: %5.3f", dist->GetMean()));
  ll->DrawLatex(.97, .65, Form("Occupancy: #int_{1}^{#infty}P(s)ds = %6.2f%%",
			       ii));
			 
#endif

  pp = cc->cd(3);
  pp->SetFillColor(0);
  pp->SetFillStyle(0);
  pp->SetBorderMode(0);
  pp->SetRightMargin(0.15);
  pp->SetTopMargin(0.02);
  pp->SetGridx();
  pp->SetGridy();
  c->GetCorrection()->Draw();

  pp = cc->cd(4);
  pp->SetFillColor(0);
  pp->SetFillStyle(0);
  pp->SetBorderMode(0);
  pp->SetRightMargin(0.15);
  pp->SetTopMargin(0.02);
  pp->SetLogz();
  pp->SetGridx();
  pp->SetGridy();
  mean->Draw();
  lmean->Draw();

  cc->cd();
}

//
// EOF
//


