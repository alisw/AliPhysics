//____________________________________________________________________
//
// $Id$
//
// Script that contains a class to draw hits, using the
// AliFMDInputHits class in the util library. 
//
// It draws the energy loss versus the p/(mq^2).  It can be overlayed
// with the Bethe-Bloc curve to show how the simulation behaves
// relative to the expected. 
//
// Use the script `Compile.C' to compile this class using ACLic. 
//
#include <TH2D.h>
#include <AliFMDHit.h>
#include <AliFMDInput.h>
#include <AliStack.h>
#include <iostream>
#include <TStyle.h>
#include <TArrayF.h>
#include <TParticle.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TDatabasePDG.h>
#include <TParticlePDG.h>
#include <TLegend.h>
#include <TArrow.h>
#include <TLatex.h>
#include <TF1.h>

/** @class DrawHits
    @brief Draw hit energy loss
    @code 
    Root> .L Compile.C
    Root> Compile("DrawHits.C")
    Root> DrawHits c
    Root> c.Run();
    @endcode
    @ingroup FMD_script
 */
class DrawHits : public AliFMDInput
{
private:
  TH2D* fElossVsPMQ; // Histogram 
  TH1D* fEloss;
  TH1D* fBetaGamma;
  TParticlePDG* fPdg;
  const Double_t fRho;
  const Double_t fBetaGammaMip;
public:
  //__________________________________________________________________
  DrawHits(const char* pdgName="pi+",
	   Int_t m=500, Double_t emin=1, Double_t emax=1000, 
	   Int_t n=900, Double_t tmin=1e-2, Double_t tmax=1e3) 
    : AliFMDInput("galice.root"),
      fElossVsPMQ(0),
      fEloss(0),
      fBetaGamma(0),
      fPdg(0), 
      fRho(2.33), // 2.33),
      fBetaGammaMip(3.4601)
  { 
    AddLoad(kKinematics);
    AddLoad(kHits);
    TDatabasePDG* pdgDB = TDatabasePDG::Instance();
    fPdg                = pdgDB->GetParticle(pdgName);
    if (!fPdg) Warning("DrawHits", "Particle %s not found", pdgName);

    TArrayF tkine(MakeLogScale(n,    tmin, tmax));
    TArrayF betag(MakeLogScale(n/10, tmin, tmax));
    TArrayF eloss(MakeLogScale(m,    emin, emax));
    TString name("elossVsPMQ");
    TString title(Form("#Delta E/#Delta x / q^{2} vs. p/m, %s", 
		       (pdgName ? pdgName : "")));
    fElossVsPMQ = new TH2D(name.Data(), title.Data(), 
			   tkine.fN-1, tkine.fArray, 
			   eloss.fN-1, eloss.fArray);
    fElossVsPMQ->SetXTitle("p/(mq^{2})=#beta#gamma/q^{2}");
    fElossVsPMQ->SetYTitle("#Delta E/#Delta x / q^{2} [MeV/cm]");
    fElossVsPMQ->Sumw2();
    fEloss = new TH1D("eloss", "#Delta E/#Delta x / q^{2}", 
		      eloss.fN-1, eloss.fArray);
    fEloss->SetFillColor(kRed);
    fEloss->SetFillStyle(3001);
    fEloss->SetXTitle("#Delta E/#Delta x / q^{2} [MeV/cm]");
    fEloss->Sumw2();

    fBetaGamma = new TH1D("betaGamma", "#beta#gamma of particles", 
			  betag.fN-1, betag.fArray);
    fBetaGamma->SetXTitle("#beta#gamma");
    fBetaGamma->SetFillColor(kBlue+1);
    fBetaGamma->SetFillStyle(3001);
  }
  //__________________________________________________________________
  Bool_t ProcessHit(AliFMDHit* hit, TParticle* p) 
  {
    if (!hit) {
      std::cout << "No hit" << std::endl;
      return kFALSE;
    }

    if (!p) {
      std::cout << "No track" << std::endl;
      return kFALSE;
    }
    // if (!p->IsPrimary()) return kTRUE;
    if (hit->IsStop()) return kTRUE;
    if (hit->Length() == 0) { 
      Warning("ProcessHit", "Hit in %s has 0 length", hit->GetName());
      return kTRUE;
    }
    
    Float_t q = hit->Q() / 3.;
    Float_t m = hit->M();
    if (m == 0 || q == 0) return kTRUE;

    TLorentzVector pp;
    p->Momentum(pp);
    Double_t betagamma = 0;
    Info("ProcessHit", "%s (%s) beta=%f", p->GetPDG()->GetName(), 
	 fStack->IsPhysicalPrimary(hit->Track()) ? "primary" : "secondary", 
	 pp.Beta());
    if (pp.Beta() <= 1 && pp.Beta() >= 0) 
      betagamma = pp.Beta() * pp.Gamma();
    fBetaGamma->Fill(betagamma);
#if 0
    if (betagamma < 10) { 
      Info("ProcessHit", "%s (%s) beta=%f gamma=%f beta*gamma=%f", 
	   p->GetPDG()->GetName(), 
	   fStack->IsPhysicalPrimary(hit->Track()) ? "primary" : 
	   "secondary", 
	   pp.Beta(), pp.Gamma(), betagamma);
      return kTRUE;
    }
#endif

    Float_t x = hit->P();
    Float_t y = hit->Edep()/hit->Length();

    x /= hit->M();
    // y /= q * q;
    fElossVsPMQ->Fill(x, y);
    fEloss->Fill(y);
    // fNHits++;
    return kTRUE;
  }
  //__________________________________________________________________
  void ShowFit(Double_t x1, Double_t y1, const char* title, 
	       TF1* f, Double_t dx=0, Double_t dy=0.05)
  {
    Double_t x = x1, y = y1;
    TLatex* latex = new TLatex(x, y, title);
    latex->SetTextFont(132);
    latex->SetTextSize(0.8*dy);
    latex->SetNDC();
    latex->Draw();
    x -= dx;
    y -= dy;
    const Double_t eqDx=0.1;
    Double_t chi2 = f->GetChisquare();
    Int_t    ndf  = f->GetNDF();
    Double_t prob = f->GetProb();
    latex->DrawLatex(x, y, "#chi^{2}/NDF");
    latex->DrawLatex(x+eqDx, y, Form("= %7.4f/%3d=%5.2f (%3d%%)", 
				     chi2, ndf, chi2/ndf, int(100*prob)));
    Int_t     n = f->GetNpar();
    Double_t* p = f->GetParameters();
    Double_t* e = f->GetParErrors();
    for (int i = 0; i < n; i++) { 
      x -= dx;
      y -= dy;
      latex->DrawLatex(x, y, f->GetParName(i));
      latex->DrawLatex(x+eqDx, y, Form("= %7.4f", p[i]));
      latex->DrawLatex(x+2*eqDx, y, Form("#pm %7.4f", e[i]));
    }
  }      
  //__________________________________________________________________
  Bool_t Finish()
  {
    gStyle->SetPalette(1);
    // gStyle->SetOptTitle(0);
    gStyle->SetTitleBorderSize(1);
    gStyle->SetTitleFillColor(0);
    gStyle->SetCanvasColor(0);
    gStyle->SetCanvasColor(0);
    gStyle->SetCanvasBorderSize(0);
    gStyle->SetPadColor(0);
    gStyle->SetPadBorderSize(0);
    gStyle->SetOptStat(0);
    TCanvas* c = new TCanvas("elossVsP", "Energy loss versus momentum", 
			     1200, 800);
    c->SetLogy();
    c->SetLogx();

    TString title(Form("%s, %d events", fElossVsPMQ->GetTitle(), fEventCount));
    fElossVsPMQ->SetTitle(title.Data());
    fElossVsPMQ->SetStats(kFALSE);
    fElossVsPMQ->Draw("AXIS");
    fElossVsPMQ->Draw("ACOLZ same");
    TGraph*  mate    = FromGFMATE();
    TGraph*  bb      = FromRPPFull();
    TGraph*  nodelta = FromRPPNoDelta();
    TGraph*  norad   = FromRPPNoRad();
    TGraph*  mean    = FromRPPMean();
    // fElossVsPMQ->Draw("ACOLZ same");
    TLegend* l       = new TLegend(.5, .6, .89, .89);
    // l->AddEntry(fElossVsPMQ, fElossVsPMQ->GetTitle(), "pf");
    l->SetFillColor(0);
    l->AddEntry(mate,        mate->GetTitle(),    "lf");
    l->AddEntry(bb,          bb->GetTitle(),      "l");
    l->AddEntry(nodelta,     nodelta->GetTitle(), "l");
    l->AddEntry(norad,       norad->GetTitle(),   "l");
    l->AddEntry(mean,        mean->GetTitle(),    "l");
    l->Draw("same");
    Double_t min = fElossVsPMQ->GetYaxis()->GetFirst();
    TArrow* a = new TArrow(fBetaGammaMip, min, fBetaGammaMip, 35*min, 03, "<|");
    a->SetAngle(30);
    a->Draw("same");
    TLatex* t = new TLatex(fBetaGammaMip, 40*min, "Minimum Ionising");
    t->SetTextSize(0.04);
    t->SetTextAlign(21);
    t->Draw("same");
    c->Modified();
    c->Update();
    c->cd();
    c->SaveAs("eloss_bethe.png");

    c = new TCanvas("cEloss", "Energy loss per unit material",
		    1200, 800);
    c->SetRightMargin(0.05);
    c->SetTopMargin(0.05);
    c->SetLeftMargin(0.05);
    fEloss->SetStats(kFALSE);
    // c->SetLogx();
    TF1* land     = new TF1("land", "landau");
    land->SetParNames("A", "MPV", "width");
    land->SetLineWidth(2);
    land->SetLineColor(kGreen+1);

    TF1* landgaus = new TF1("landgaus", "gaus(0)+landau(3)");
    landgaus->SetParNames("A", "#mu", "#sigma", "B", "MPV", "width");
    landgaus->SetLineWidth(2);
    landgaus->SetLineColor(kMagenta+1);
    TGraph*  corr  = GetCorr();
    TGraph*  resp  = GetResp();
    if (fEloss->GetEntries() != 0) { 
      c->SetLogy();
      fEloss->Scale(1. / fEloss->GetEntries());
      fEloss->GetXaxis()->SetRangeUser(1, 10);
      fEloss->Fit(land, "+", "", 2, 10);
      landgaus->SetParameters(land->GetParameter(0) / 100, 
			      land->GetParameter(1), 
			      land->GetParameter(2), 
			      land->GetParameter(0),
			      land->GetParameter(1), 
			      land->GetParameter(2));
      fEloss->Fit(landgaus, "+", "", 1, 10);
      fEloss->DrawCopy("HIST SAME");
    }
    
    // fEloss->DrawCopy("E SAME");
    // land->Draw("same");
    // landgaus->Draw("same");
    Double_t max = fEloss->GetMaximum();
    Double_t* x   = resp->GetX();
    Double_t* y   = resp->GetY();
    TGraph*   g   = new TGraph(resp->GetN());
    g->SetName(Form("%sCorr", resp->GetName()));
    g->SetTitle(resp->GetTitle());
    g->SetLineStyle(resp->GetLineStyle());
    g->SetLineColor(resp->GetLineColor());
    g->SetLineWidth(resp->GetLineWidth());
    Double_t  xs2 = corr->Eval(fBetaGammaMip);
    Double_t xss   = 1.1;
    Double_t  xs  = fRho * xss;
    std::cout << "Correction factor: " << xs2 << std::endl;
    for (Int_t i = 0; i < g->GetN(); i++) 
      g->SetPoint(i, x[i] * xs, y[i] * max);
    g->Draw("C same");
    
    l = new TLegend(.05, .6, .4, .95);
    l->SetFillColor(0);
    l->SetBorderSize(1);
    l->AddEntry(fEloss, fEloss->GetTitle(), "lf");
    if (land) 
      l->AddEntry(land,   Form("Landau fit\t- #chi^{2}/NDF=%7.5f", 
			       land->GetChisquare()/land->GetNDF()), "l");
    if (landgaus) 
      l->AddEntry(landgaus,   
		  Form("Landau+Gauss fit\t- #chi^{2}/NDF=%7.5f", 
		       landgaus->GetChisquare()/landgaus->GetNDF()), "l");
    l->AddEntry(resp,   Form("f(%s#Delta/x) 320#mum Si [RPP fig 27.7]",
			     xss != 1 ? Form("%4.2f#times", xss) : ""), 
		"l");
    l->Draw("same");
    
    fEloss->GetYaxis()->SetRangeUser(1e-4, 100);
    ShowFit(0.45,.9, "Landau+Gaus", landgaus, 0, 0.02);
    ShowFit(0.70,.9, "Landau", land, 0, 0.02);

    c->Modified();
    c->Update();
    c->cd();
    c->SaveAs("eloss_landau.png");


    c = new TCanvas("cBetaGamma", "beta gamma", 1200, 800);
    c->SetLogx();
    fBetaGamma->Draw();
    Int_t mipbin = fBetaGamma->FindBin(fBetaGammaMip) + 1;
    Int_t maxbin = fBetaGamma->GetNbinsX();
    Int_t total  = fBetaGamma->Integral();
    Int_t over   = fBetaGamma->Integral(mipbin,maxbin);
    TH1*  res    = (TH1*)fBetaGamma->Clone("overMip");
    res->SetFillColor(kRed+1);
    for (int i = 0; i < mipbin; i++) 
      res->SetBinContent(i, 0);
    res->Draw("same");
    std::cout << "Percentage over MIP : " << float(over) / total << std::endl;

    Double_t yy = fBetaGamma->GetBinContent(mipbin) * 1.1; 
    a = new TArrow(fBetaGammaMip, 0, fBetaGammaMip, yy, 3, "<|");
    a->SetAngle(30);
    a->Draw("same");
    t = new TLatex(fBetaGammaMip, yy, "Minimum Ionising");
    t->SetTextSize(0.04);
    t->SetTextAlign(21);
    t->Draw("same");
    c->Modified();
    c->Update();
    c->cd();
    c->SaveAs("eloss_betagamma.png");
    
    return kTRUE;
  }

  /** Scale a graph by density (multiply) and mass (divide). 
      @param graph Graph to scale 
      @param density If @c true, scale by the Si density
      ($\rho=2.33$/cm^3$).  The y axis is assumed to have units of
      $MeVg^{-1}cm^2$. 
      @param mass Mass to scale with. The x axis is assumed to be the
      kinetic energy of the particles in units of $GeV$.  */
  void ScaleGraph(TGraph* graph, bool density=true, double mass=1) 
  {
      Double_t*      x   = graph->GetX();
      Double_t*      y   = graph->GetY();
      const Double_t rho = (density ? fRho : 1);
      for (Int_t i = 0; i < graph->GetN(); i++) 
	graph->SetPoint(i, x[i] / mass, y[i] * rho); 
  }    
  /** Draw pure Bethe-Bloc from Review of Particle Physics, fig. 27.1 
      @return TGraph object */ 
  TGraph* FromRPPFull() 
  {
    static TGraph* graph = 0;
    if (!graph) { 
      graph = new TGraph(20);
      graph->GetHistogram()->SetXTitle("#beta#gamma");
      graph->GetHistogram()->SetYTitle("#Delta E/#Delta x [MeV/cm]");
      graph->SetFillColor(0);
      graph->SetLineColor(kRed+1);
      graph->SetLineStyle(2);
      graph->SetLineWidth(2);
      graph->SetName("full_stop");
      graph->SetTitle("Stopping (MeVcm^{2}/g) [RPP fig 27.1]");
      graph->SetPoint(0,0.001461622,40.17542);
      graph->SetPoint(1,0.003775053,91.28429);
      graph->SetPoint(2,0.01178769,202.7359);
      graph->SetPoint(3,0.01722915,212.1938);
      graph->SetPoint(4,0.03162278,172.8318);
      graph->SetPoint(5,0.06028646,91.28429);
      graph->SetPoint(6,0.09506529,51.62633);
      graph->SetPoint(7,0.433873,5.281682);
      graph->SetPoint(8,1.255744,1.808947);
      graph->SetPoint(9,2.393982,1.440177);
      graph->SetPoint(10,3.499097,1.407715);
      graph->SetPoint(11,10.92601,1.542122);
      graph->SetPoint(12,60.28646,1.85066);
      graph->SetPoint(13,236.3885,2.121938);
      graph->SetPoint(14,468.0903,2.324538);
      graph->SetPoint(15,1208.976,2.987085);
      graph->SetPoint(16,6670.768,7.961412);
      graph->SetPoint(17,23341.67,24.3298);
      graph->SetPoint(18,110651.2,104.6651);
      graph->SetPoint(19,264896.9,260.5203);
      ScaleGraph(graph);
    }
    graph->Draw("C same");
    return graph;
  }

  /** Draw pure Bethe-Bloc from Review of Particle Physics, fig. 27.1,
      but without delta electrons  
      @return TGraph object */ 
  TGraph* FromRPPNoDelta() 
  {
    static TGraph* graph = 0;
    if (!graph) { 
      graph = new TGraph(20);
      graph->SetName("stop_nodelta");
      graph->SetTitle("Stopping w/o #delta's [RPP fig 27.1]");
      graph->GetHistogram()->SetYTitle("(MeVcm^{2}/g)");
      graph->GetHistogram()->SetXTitle("#beta#gamma");
      graph->SetFillColor(0);
      graph->SetLineColor(kGreen+1);
      graph->SetLineStyle(2);
      graph->SetLineWidth(2);
      graph->SetPoint(0,0.001461622,40.17542);
      graph->SetPoint(1,0.003775053,91.28429);
      graph->SetPoint(2,0.01178769,202.7359);
      graph->SetPoint(3,0.01722915,212.1938);
      graph->SetPoint(4,0.03162278,172.8318);
      graph->SetPoint(5,0.06028646,91.28429);
      graph->SetPoint(6,0.09506529,51.62633);
      graph->SetPoint(7,0.433873,5.281682);
      graph->SetPoint(8,1.255744,1.808947);
      graph->SetPoint(9,2.304822,1.473387);
      graph->SetPoint(10,3.921088,1.473387);
      graph->SetPoint(11,8.064796,1.614064);
      graph->SetPoint(12,26.15667,1.936996);
      graph->SetPoint(13,264.8969,2.489084);
      graph->SetPoint(14,544.8334,2.665278);
      graph->SetPoint(15,1163.949,2.853945);
      graph->SetPoint(16,5312.204,3.19853);
      graph->SetPoint(17,15374.93,3.424944);
      graph->SetPoint(18,49865.73,3.667384);
      graph->SetPoint(19,634158.5,4.110185);
      ScaleGraph(graph);
    }
    graph->Draw("C same");
    return graph;
  }

  /** Draw pure Bethe-Bloc from Review of Particle Physics, fig. 27.1,
      but without delta electrons  
      @return TGraph object */ 
  TGraph* FromRPPNoRad() 
  {
    static TGraph* graph = 0;
    if (!graph) { 
      graph = new TGraph(18);
      graph->SetName("norad_stop");
      graph->SetTitle("Stopping w/o radiative loss [RPP fig. 27.1]");
      graph->GetHistogram()->SetYTitle("(MeVcm^{2}/g)");
      graph->GetHistogram()->SetXTitle("#beta#gamma");
      graph->SetFillColor(0);
      graph->SetLineColor(kBlue+1);
      graph->SetLineWidth(2);
      graph->SetLineStyle(2);
      graph->SetPoint(0,0.001,24.3298);
      graph->SetPoint(1,0.003117649,74.35105);
      graph->SetPoint(2,0.008675042,172.8318);
      graph->SetPoint(3,0.01782497,212.1938);
      graph->SetPoint(4,0.02704573,189.3336);
      graph->SetPoint(5,0.07481082,70.29816);
      graph->SetPoint(6,0.3300035,8.524974);
      graph->SetPoint(7,0.819559,2.489084);
      graph->SetPoint(8,1.447084,1.651284);
      graph->SetPoint(9,2.555097,1.440177);
      graph->SetPoint(10,4.026598,1.407715);
      graph->SetPoint(11,32.38084,1.728318);
      graph->SetPoint(12,97.19733,1.893336);
      graph->SetPoint(13,1732.539,2.170869);
      graph->SetPoint(14,11098.58,2.324538);
      graph->SetPoint(15,32075.46,2.378141);
      graph->SetPoint(16,221655.8,2.546482);
      graph->SetPoint(17,593830.6,2.605203);
      ScaleGraph(graph);
    }
    graph->Draw("C same");
    return graph;
  }

  /** Draw pure Bethe-Bloc from Review of Particle Physics, fig. 27.6 
      @return TGraph object */ 
  TGraph* FromRPPMean() 
  {
    static TGraph* graph = 0;
    if (!graph) { 
      graph = new TGraph(12);
      graph->SetName("mean_eloss");
      graph->SetTitle("Mean #Delta E/#Delta x - "
		      "electronic only  [RPP fig. 27.6]");
      graph->GetHistogram()->SetYTitle("(MeVcm^{2}/g)");
      graph->GetHistogram()->SetXTitle("#mu E_{kin} (GeV)");
      graph->SetFillColor(0);
      graph->SetLineStyle(2);
      graph->SetLineWidth(2);
      graph->SetLineColor(kMagenta+1);
      graph->SetMarkerStyle(21);
      graph->SetMarkerSize(0.6);
      graph->SetPoint(0,0.1,1.346561);
      graph->SetPoint(1,0.1435819,1.230159);
      graph->SetPoint(2,0.2061576,1.156085);
      graph->SetPoint(3,0.3698076,1.124339);
      graph->SetPoint(4,0.4620113,1.124339);
      graph->SetPoint(5,0.8521452,1.145503);
      graph->SetPoint(6,1.909707,1.177249);
      graph->SetPoint(7,4.048096,1.198413);
      graph->SetPoint(8,12.66832,1.219577);
      graph->SetPoint(9,48.17031,1.230159);
      graph->SetPoint(10,285.8863,1.230159);
      graph->SetPoint(11,894.6674,1.230159);
      const Double_t m   = 0.10566; // Muon 
      ScaleGraph(graph, true, m);
    }
    graph->Draw("C same");
    return graph;
  }

  /** Draw energy loss as obtained from GEANT 3.21 GFMATE. 
      @return TGraph object */
  TGraph* FromGFMATE() 
  {
    static TGraphErrors* gre = 0;
    if (!gre) {
      gre = new TGraphErrors(91);
      gre->SetName("ELOSS");
      gre->SetTitle("Energy loss 300#mu Si [GFMATE]");
      gre->GetHistogram()->SetXTitle("#beta#gamma");
      gre->GetHistogram()->SetYTitle("#Delta E/#Delta x [MeV/cm]");
      gre->SetFillColor(kGray);
      gre->SetFillStyle(3001); // 0 /* 3002 */);
      gre->SetLineColor(kGray+1);
      gre->SetLineStyle(1);
      gre->SetLineWidth(2);
      gre->SetPoint(0,7.16486e-05,1218.84);
      gre->SetPoint(1,9.25378e-05,1221.38);
      gre->SetPoint(2,0.000119517,1180.12);
      gre->SetPoint(3,0.000154362,1100.31);
      gre->SetPoint(4,0.000199367,996.621);
      gre->SetPoint(5,0.000257492,886.005);
      gre->SetPoint(6,0.000332563,780.483);
      gre->SetPoint(7,0.000429522,684.927);
      gre->SetPoint(8,0.000554749,599.407);
      gre->SetPoint(9,0.000716486,522.375);
      gre->SetPoint(10,0.000925378,452.497);
      gre->SetPoint(11,0.00119517,389.101);
      gre->SetPoint(12,0.00154362,331.974);
      gre->SetPoint(13,0.00199367,280.969);
      gre->SetPoint(14,0.00257492,235.689);
      gre->SetPoint(15,0.00332564,196.156);
      gre->SetPoint(16,0.00429522,162.402);
      gre->SetPoint(17,0.00554749,133.87);
      gre->SetPoint(18,0.00716486,109.959);
      gre->SetPoint(19,0.00925378,90.2035);
      gre->SetPoint(20,0.0119517,74.1317);
      gre->SetPoint(21,0.0154362,60.8988);
      gre->SetPoint(22,0.0199367,49.9915);
      gre->SetPoint(23,0.0257492,40.9812);
      gre->SetPoint(24,0.0332564,33.5739);
      gre->SetPoint(25,0.0429522,27.5127);
      gre->SetPoint(26,0.0554749,22.5744);
      gre->SetPoint(27,0.0716486,18.5674);
      gre->SetPoint(28,0.0925378,15.3292);
      gre->SetPoint(29,0.119517,12.7231);
      gre->SetPoint(30,0.154362,10.6352);
      gre->SetPoint(31,0.199367,8.97115);
      gre->SetPoint(32,0.257492,7.65358);
      gre->SetPoint(33,0.332564,6.61909);
      gre->SetPoint(34,0.429522,5.79837);
      gre->SetPoint(35,0.554749,5.148);
      gre->SetPoint(36,0.716486,4.65024);
      gre->SetPoint(37,0.925378,4.27671);
      gre->SetPoint(38,1.19517,3.99831);
      gre->SetPoint(39,1.54362,3.79877);
      gre->SetPoint(40,1.99367,3.6629);
      gre->SetPoint(41,2.57492,3.57594);
      gre->SetPoint(42,3.32564,3.52565);
      gre->SetPoint(43,4.29522,3.50206);
      gre->SetPoint(44,5.54749,3.49715);
      gre->SetPoint(45,7.16486,3.50467);
      gre->SetPoint(46,9.25378,3.51988);
      gre->SetPoint(47,11.9517,3.53932);
      gre->SetPoint(48,15.4362,3.56054);
      gre->SetPoint(49,19.9367,3.58189);
      gre->SetPoint(50,25.7492,3.60231);
      gre->SetPoint(51,33.2564,3.62113);
      gre->SetPoint(52,42.9522,3.638);
      gre->SetPoint(53,55.4749,3.65275);
      gre->SetPoint(54,71.6486,3.66537);
      gre->SetPoint(55,92.5378,3.67586);
      gre->SetPoint(56,119.517,3.68433);
      gre->SetPoint(57,154.362,3.69105);
      gre->SetPoint(58,199.367,3.6962);
      gre->SetPoint(59,257.492,3.69997);
      gre->SetPoint(60,332.564,3.70257);
      gre->SetPoint(61,429.522,3.70421);
      gre->SetPoint(62,554.749,3.70511);
      gre->SetPoint(63,716.486,3.7055);
      gre->SetPoint(64,925.378,3.70559);
      gre->SetPoint(65,1195.17,3.70558);
      gre->SetPoint(66,1543.62,3.70557);
      gre->SetPoint(67,1993.67,3.70555);
      gre->SetPoint(68,2574.92,3.70553);
      gre->SetPoint(69,3325.64,3.70552);
      gre->SetPoint(70,4295.22,3.7055);
      gre->SetPoint(71,5547.49,3.70548);
      gre->SetPoint(72,7164.86,3.70547);
      gre->SetPoint(73,9253.78,3.70545);
      gre->SetPoint(74,11951.7,3.70544);
      gre->SetPoint(75,15436.2,3.70544);
      gre->SetPoint(76,19936.7,3.70544);
      gre->SetPoint(77,25749.2,3.70544);
      gre->SetPoint(78,33256.4,3.70544);
      gre->SetPoint(79,42952.2,3.70544);
      gre->SetPoint(80,55474.9,3.70544);
      gre->SetPoint(81,71648.6,3.70544);
      gre->SetPoint(82,92537.8,3.70544);
      gre->SetPoint(83,119517,3.70544);
      gre->SetPoint(84,154362,3.70544);
      gre->SetPoint(85,199367,3.70544);
      gre->SetPoint(86,257492,3.70544);
      gre->SetPoint(87,332563,3.70544);
      gre->SetPoint(88,429522,3.70544);
      gre->SetPoint(89,554749,3.70544);
      gre->SetPoint(90,716486,3.70544);
      // Double_t* x = gre->GetX();
      Double_t* y = gre->GetY();
      for (Int_t i = 0; i < gre->GetN(); i++) 
	gre->SetPointError(i, 0, 2 * 0.1 * y[i]); // ! 1 sigma
    }
    gre->Draw("c3 same");
    return gre;
  }

  /** Get the response functin @f$ f(\Delta_p/x)@f$ from Review of
      Particle Physics (fig. 27.7).  It is scaled to the value at
      MPV. */ 
  TGraph* GetResp()
  {
    static TGraph*  graph = 0;
    if (!graph) {
      graph = new TGraph;
      graph->SetName("si_resp");
      graph->SetTitle("f(#Delta/x) scaled to the MPV value ");
      graph->GetHistogram()->SetXTitle("#Delta/x (MeVcm^{2}/g)");
      graph->GetHistogram()->SetYTitle("f(#Delta/x)");
      graph->SetLineColor(kBlue+1);
      graph->SetLineWidth(2);
      graph->SetFillColor(kBlue+1);
      graph->SetMarkerStyle(21);
      graph->SetMarkerSize(0.6);
#if 0
      // Figure 27.7 or Review of Particle physics - Straggeling function in 
      // silicon of 500 MeV Pions, normalised to unity at the most probable 
      // value.   
      graph->SetPoint(0,0.808094,0.00377358);
      graph->SetPoint(1,0.860313,0.0566038);
      graph->SetPoint(2,0.891645,0.116981);
      graph->SetPoint(3,0.912533,0.181132);
      graph->SetPoint(4,0.928198,0.260377);
      graph->SetPoint(5,0.938642,0.320755);
      graph->SetPoint(6,0.954308,0.377358);
      graph->SetPoint(7,0.964752,0.433962);
      graph->SetPoint(8,0.975196,0.490566);
      graph->SetPoint(9,0.98564,0.550943);
      graph->SetPoint(10,0.996084,0.611321);
      graph->SetPoint(11,1.00653,0.667925);
      graph->SetPoint(12,1.02219,0.732075);
      graph->SetPoint(13,1.03264,0.784906);
      graph->SetPoint(14,1.0483,0.845283);
      graph->SetPoint(15,1.06397,0.901887);
      graph->SetPoint(16,1.09008,0.958491);
      graph->SetPoint(17,1.10574,0.984906);
      graph->SetPoint(18,1.13708,1);
      graph->SetPoint(19,1.13708,1);
      graph->SetPoint(20,1.15796,0.988679);
      graph->SetPoint(21,1.17363,0.966038);
      graph->SetPoint(22,1.19974,0.916981);
      graph->SetPoint(23,1.2154,0.89434);
      graph->SetPoint(24,1.23629,0.837736);
      graph->SetPoint(25,1.2624,0.784906);
      graph->SetPoint(26,1.28329,0.724528);
      graph->SetPoint(27,1.3094,0.664151);
      graph->SetPoint(28,1.32507,0.611321);
      graph->SetPoint(29,1.3564,0.550943);
      graph->SetPoint(30,1.41384,0.445283);
      graph->SetPoint(31,1.44517,0.392453);
      graph->SetPoint(32,1.48695,0.335849);
      graph->SetPoint(33,1.52872,0.286792);
      graph->SetPoint(34,1.58094,0.237736);
      graph->SetPoint(35,1.63838,0.196226);
      graph->SetPoint(36,1.68016,0.169811);
      graph->SetPoint(37,1.75326,0.135849);
      graph->SetPoint(38,1.81593,0.113208);
      graph->SetPoint(39,1.89426,0.0981132);
      graph->SetPoint(40,1.96214,0.0830189);
      graph->SetPoint(41,2.0718,0.0641509);
      graph->SetPoint(42,2.19191,0.0490566);
      graph->SetPoint(43,2.31723,0.0415094);
      graph->SetPoint(44,2.453,0.0301887);
      graph->SetPoint(45,2.53133,0.0264151);
      graph->SetPoint(46,2.57833,0.0264151);
#else
      graph->SetPoint(0,0.8115124,0.009771987);
      graph->SetPoint(1,0.9198646,0.228013);
      graph->SetPoint(2,0.996614,0.5895765);
      graph->SetPoint(3,1.041761,0.8241042);
      graph->SetPoint(4,1.059819,0.8794788);
      graph->SetPoint(5,1.077878,0.9348534);
      graph->SetPoint(6,1.100451,0.980456);
      graph->SetPoint(7,1.141084,0.9967427);
      graph->SetPoint(8,1.204289,0.9153094);
      graph->SetPoint(9,1.276524,0.742671);
      graph->SetPoint(10,1.402935,0.465798);
      graph->SetPoint(11,1.515801,0.3029316);
      graph->SetPoint(12,1.73702,0.1465798);
      graph->SetPoint(13,1.985327,0.08143322);
      graph->SetPoint(14,2.301354,0.04234528);
      graph->SetPoint(15,2.56772,0.02931596);
#endif
    }
    return graph;
  }

  /** Get the correction to Bethe-Bloc from Review of Particle Physics
      (fig 27.8). 
  */
  TGraph* GetCorr() 
  {
    static TGraph* graph = 0;
    if (!graph) {
      graph = new TGraph(14);
      graph->SetName("graph");
      graph->SetTitle("(#Delta_{p}/x)/(dE/dx)|_{mip} for 320#mu Si");
      graph->GetHistogram()->SetXTitle("#beta#gamma = p/m");
      graph->SetFillColor(1);
      graph->SetLineColor(7);
      graph->SetMarkerStyle(21);
      graph->SetMarkerSize(0.6);
      graph->SetPoint(0,1.196058,0.9944915);
      graph->SetPoint(1,1.28502,0.9411017);
      graph->SetPoint(2,1.484334,0.8559322);
      graph->SetPoint(3,1.984617,0.7491525);
      graph->SetPoint(4,2.658367,0.6983051);
      graph->SetPoint(5,3.780227,0.6779661);
      graph->SetPoint(6,4.997358,0.6741525);
      graph->SetPoint(7,8.611026,0.684322);
      graph->SetPoint(8,15.28296,0.6995763);
      graph->SetPoint(9,41.54516,0.7186441);
      graph->SetPoint(10,98.91461,0.7288136);
      graph->SetPoint(11,203.2734,0.7326271);
      graph->SetPoint(12,505.6421,0.7338983);
      graph->SetPoint(13,896.973,0.7338983);
    }
    return graph;
  }
  
  ClassDef(DrawHits,0);
};

//____________________________________________________________________
//
// EOF
//
