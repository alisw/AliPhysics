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
#include <iostream>
#include <TStyle.h>
#include <TArrayF.h>
#include <TParticle.h>
#include <TCanvas.h>
#include <TGraphErrors.h>

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
  Int_t fPdg;
public:
  //__________________________________________________________________
  TArrayF MakeLogScale(Int_t n, Double_t min, Double_t max) 
  {
    TArrayF bins(n+1);
    bins[0]      = min;
    if (n <= 20) {
      for (Int_t i = 1; i < n+1; i++) bins[i] = bins[i-1] + (max-min)/n;
      return bins;
    }
    Float_t dp   = n / TMath::Log10(max / min);
    Float_t pmin = TMath::Log10(min);
    for (Int_t i = 1; i < n+1; i++) {
      Float_t p = pmin + i / dp;
      bins[i]   = TMath::Power(10, p);
    }
    return bins;
  }
  //__________________________________________________________________
  DrawHits(Int_t m=1000, Double_t emin=1, Double_t emax=1000, 
	   Int_t n=900, Double_t tmin=1e-2, Double_t tmax=1e3, 
	   Int_t pdg=211) 
    : fPdg(pdg)
  { 
    AddLoad(kKinematics);
    AddLoad(kHits);
    TArrayF tkine(MakeLogScale(n, tmin, tmax));
    TArrayF eloss(MakeLogScale(m, emin, emax));
    TString name(Form("elossVsP%s", (fPdg == 0 ? "MQ" : "")));
    TString title(Form("#Delta E/#Delta x vs. p%s", 
		       (fPdg == 0 ? "/(mq^{2})" : "")));
    fElossVsPMQ = new TH2D(name.Data(), title.Data(), 
			   tkine.fN-1, tkine.fArray, 
			   eloss.fN-1, eloss.fArray);
    fElossVsPMQ->SetXTitle(Form("p%s", (fPdg == 0 ? "/(mq^{2})" : " [GeV]")));
    fElossVsPMQ->SetYTitle("#Delta E/#Delta x [MeV/cm]");
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
    if (!p->IsPrimary()) return kTRUE;
    if (hit->IsStop()) return kTRUE;
    Float_t x = hit->P();
    Float_t y = hit->Edep()/hit->Length();
    if (fPdg != 0) {
      if (TMath::Abs(hit->Pdg()) != fPdg) return kTRUE;
    }
    else {
      Float_t q = hit->Q() / 3.;
      if (hit->M() == 0 || q == 0) return kTRUE;
      x /= hit->M();
      y /= q * q;
    }
    fElossVsPMQ->Fill(x, y);
    return kTRUE;
  }
  //__________________________________________________________________
  Bool_t Finish()
  {
    TCanvas* c = new TCanvas("c", "C");
    c->SetLogy();
    c->SetLogx();
    gStyle->SetPalette(1);
    gStyle->SetOptTitle(0);
    gStyle->SetCanvasColor(0);
    gStyle->SetCanvasBorderSize(0);
    gStyle->SetPadColor(0);
    gStyle->SetPadBorderSize(0);
    fElossVsPMQ->SetStats(kFALSE);
    fElossVsPMQ->Draw("COLZ box");
    c->Modified();
    c->Update();
    c->cd();
    return kTRUE;
  }
  void SuperImposeBetheBloc() 
  {
    // This is for pi+, made with MakeXsection.C and DrawXsection.C
    TGraphErrors *gre = new TGraphErrors(91);
    gre->SetName("BetheBlocPi");
    gre->SetTitle("BetheBlocPi");
    gre->SetFillColor(6);
    gre->SetFillStyle(3001);
    gre->SetLineWidth(2);
    gre->SetPoint(0,7.16486e-05,1218.84);
    gre->SetPointError(0,0,609.418);
    gre->SetPoint(1,9.25378e-05,1221.38);
    gre->SetPointError(1,0,610.689);
    gre->SetPoint(2,0.000119517,1180.12);
    gre->SetPointError(2,0,590.058);
    gre->SetPoint(3,0.000154362,1100.31);
    gre->SetPointError(3,0,550.156);
    gre->SetPoint(4,0.000199367,996.621);
    gre->SetPointError(4,0,498.31);
    gre->SetPoint(5,0.000257492,886.005);
    gre->SetPointError(5,0,443.003);
    gre->SetPoint(6,0.000332563,780.483);
    gre->SetPointError(6,0,390.241);
    gre->SetPoint(7,0.000429522,684.927);
    gre->SetPointError(7,0,342.463);
    gre->SetPoint(8,0.000554749,599.407);
    gre->SetPointError(8,0,299.703);
    gre->SetPoint(9,0.000716486,522.375);
    gre->SetPointError(9,0,261.187);
    gre->SetPoint(10,0.000925378,452.497);
    gre->SetPointError(10,0,226.249);
    gre->SetPoint(11,0.00119517,389.101);
    gre->SetPointError(11,0,194.551);
    gre->SetPoint(12,0.00154362,331.974);
    gre->SetPointError(12,0,165.987);
    gre->SetPoint(13,0.00199367,280.969);
    gre->SetPointError(13,0,140.485);
    gre->SetPoint(14,0.00257492,235.689);
    gre->SetPointError(14,0,117.844);
    gre->SetPoint(15,0.00332564,196.156);
    gre->SetPointError(15,0,98.078);
    gre->SetPoint(16,0.00429522,162.402);
    gre->SetPointError(16,0,81.2012);
    gre->SetPoint(17,0.00554749,133.87);
    gre->SetPointError(17,0,66.9351);
    gre->SetPoint(18,0.00716486,109.959);
    gre->SetPointError(18,0,54.9797);
    gre->SetPoint(19,0.00925378,90.2035);
    gre->SetPointError(19,0,45.1017);
    gre->SetPoint(20,0.0119517,74.1317);
    gre->SetPointError(20,0,37.0658);
    gre->SetPoint(21,0.0154362,60.8988);
    gre->SetPointError(21,0,30.4494);
    gre->SetPoint(22,0.0199367,49.9915);
    gre->SetPointError(22,0,24.9957);
    gre->SetPoint(23,0.0257492,40.9812);
    gre->SetPointError(23,0,20.4906);
    gre->SetPoint(24,0.0332564,33.5739);
    gre->SetPointError(24,0,16.7869);
    gre->SetPoint(25,0.0429522,27.5127);
    gre->SetPointError(25,0,13.7563);
    gre->SetPoint(26,0.0554749,22.5744);
    gre->SetPointError(26,0,11.2872);
    gre->SetPoint(27,0.0716486,18.5674);
    gre->SetPointError(27,0,9.28372);
    gre->SetPoint(28,0.0925378,15.3292);
    gre->SetPointError(28,0,7.66462);
    gre->SetPoint(29,0.119517,12.7231);
    gre->SetPointError(29,0,6.36156);
    gre->SetPoint(30,0.154362,10.6352);
    gre->SetPointError(30,0,5.31759);
    gre->SetPoint(31,0.199367,8.97115);
    gre->SetPointError(31,0,4.48558);
    gre->SetPoint(32,0.257492,7.65358);
    gre->SetPointError(32,0,3.82679);
    gre->SetPoint(33,0.332564,6.61909);
    gre->SetPointError(33,0,3.30955);
    gre->SetPoint(34,0.429522,5.81614);
    gre->SetPointError(34,0,2.90807);
    gre->SetPoint(35,0.554749,5.20286);
    gre->SetPointError(35,0,2.60143);
    gre->SetPoint(36,0.716486,4.74533);
    gre->SetPointError(36,0,2.37267);
    gre->SetPoint(37,0.925378,4.40987);
    gre->SetPointError(37,0,2.20494);
    gre->SetPoint(38,1.19517,4.17077);
    gre->SetPointError(38,0,2.08538);
    gre->SetPoint(39,1.54362,4.014);
    gre->SetPointError(39,0,2.007);
    gre->SetPoint(40,1.99367,3.92577);
    gre->SetPointError(40,0,1.96288);
    gre->SetPoint(41,2.57492,3.89199);
    gre->SetPointError(41,0,1.946);
    gre->SetPoint(42,3.32564,3.90063);
    gre->SetPointError(42,0,1.95032);
    gre->SetPoint(43,4.29522,3.94146);
    gre->SetPointError(43,0,1.97073);
    gre->SetPoint(44,5.54749,4.00597);
    gre->SetPointError(44,0,2.00299);
    gre->SetPoint(45,7.16486,4.08725);
    gre->SetPointError(45,0,2.04362);
    gre->SetPoint(46,9.25378,4.17985);
    gre->SetPointError(46,0,2.08993);
    gre->SetPoint(47,11.9517,4.27962);
    gre->SetPointError(47,0,2.13981);
    gre->SetPoint(48,15.4362,4.38347);
    gre->SetPointError(48,0,2.19174);
    gre->SetPoint(49,19.9367,4.48919);
    gre->SetPointError(49,0,2.2446);
    gre->SetPoint(50,25.7492,4.59523);
    gre->SetPointError(50,0,2.29762);
    gre->SetPoint(51,33.2564,4.70052);
    gre->SetPointError(51,0,2.35026);
    gre->SetPoint(52,42.9522,4.80435);
    gre->SetPointError(52,0,2.40218);
    gre->SetPoint(53,55.4749,4.90625);
    gre->SetPointError(53,0,2.45312);
    gre->SetPoint(54,71.6486,5.00589);
    gre->SetPointError(54,0,2.50295);
    gre->SetPoint(55,92.5378,5.10279);
    gre->SetPointError(55,0,2.55139);
    gre->SetPoint(56,119.517,5.19654);
    gre->SetPointError(56,0,2.59827);
    gre->SetPoint(57,154.362,5.28758);
    gre->SetPointError(57,0,2.64379);
    gre->SetPoint(58,199.367,5.37581);
    gre->SetPointError(58,0,2.6879);
    gre->SetPoint(59,257.492,5.46109);
    gre->SetPointError(59,0,2.73055);
    gre->SetPoint(60,332.564,5.54335);
    gre->SetPointError(60,0,2.77167);
    gre->SetPoint(61,429.522,5.62248);
    gre->SetPointError(61,0,2.81124);
    gre->SetPoint(62,554.749,5.69843);
    gre->SetPointError(62,0,2.84922);
    gre->SetPoint(63,716.486,5.77122);
    gre->SetPointError(63,0,2.88561);
    gre->SetPoint(64,925.378,5.84093);
    gre->SetPointError(64,0,2.92046);
    gre->SetPoint(65,1195.17,5.9077);
    gre->SetPointError(65,0,2.95385);
    gre->SetPoint(66,1543.62,5.97165);
    gre->SetPointError(66,0,2.98582);
    gre->SetPoint(67,1993.67,6.03292);
    gre->SetPointError(67,0,3.01646);
    gre->SetPoint(68,2574.92,6.09171);
    gre->SetPointError(68,0,3.04586);
    gre->SetPoint(69,3325.64,6.14827);
    gre->SetPointError(69,0,3.07413);
    gre->SetPoint(70,4295.22,6.20286);
    gre->SetPointError(70,0,3.10143);
    gre->SetPoint(71,5547.49,6.25577);
    gre->SetPointError(71,0,3.12788);
    gre->SetPoint(72,7164.86,6.30725);
    gre->SetPointError(72,0,3.15363);
    gre->SetPoint(73,9253.78,6.35757);
    gre->SetPointError(73,0,3.17878);
    gre->SetPoint(74,11951.7,6.38446);
    gre->SetPointError(74,0,3.19223);
    gre->SetPoint(75,15436.2,6.38446);
    gre->SetPointError(75,0,3.19223);
    gre->SetPoint(76,19936.7,6.38446);
    gre->SetPointError(76,0,3.19223);
    gre->SetPoint(77,25749.2,6.38446);
    gre->SetPointError(77,0,3.19223);
    gre->SetPoint(78,33256.4,6.38446);
    gre->SetPointError(78,0,3.19223);
    gre->SetPoint(79,42952.2,6.38446);
    gre->SetPointError(79,0,3.19223);
    gre->SetPoint(80,55474.9,6.38446);
    gre->SetPointError(80,0,3.19223);
    gre->SetPoint(81,71648.6,6.38446);
    gre->SetPointError(81,0,3.19223);
    gre->SetPoint(82,92537.8,6.38446);
    gre->SetPointError(82,0,3.19223);
    gre->SetPoint(83,119517,6.38446);
    gre->SetPointError(83,0,3.19223);
    gre->SetPoint(84,154362,6.38446);
    gre->SetPointError(84,0,3.19223);
    gre->SetPoint(85,199367,6.38446);
    gre->SetPointError(85,0,3.19223);
    gre->SetPoint(86,257492,6.38446);
    gre->SetPointError(86,0,3.19223);
    gre->SetPoint(87,332563,6.38446);
    gre->SetPointError(87,0,3.19223);
    gre->SetPoint(88,429522,6.38446);
    gre->SetPointError(88,0,3.19223);
    gre->SetPoint(89,554749,6.38446);
    gre->SetPointError(89,0,3.19223);
    gre->SetPoint(90,716486,6.38446);
    gre->SetPointError(90,0,3.19223);
    gre->Draw("l same");
    gre->DrawClone("l3 same");
    gPad->Modified();
    gPad->Update();
    gPad->cd();
  }
  
  ClassDef(DrawHits,0);
};

//____________________________________________________________________
//
// EOF
//
