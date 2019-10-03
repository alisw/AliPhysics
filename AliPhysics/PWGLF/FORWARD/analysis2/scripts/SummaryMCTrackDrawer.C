#include "SummaryDrawer.C"
#ifndef __CINT__
# include <TGraph.h>
# include <TGraphErrors.h>
# include <TF1.h>
# include <TArrow.h>
#else
class TGraph;
class TFile;
class TGraphErrors;
#endif

/**
 * Class to draw a summary of the AOD production
 *
 * @par Input: 
 * - The merged <tt>forward.root</tt> file.
 *   If the file isn't merged, it should still work. 
 *
 * @par Output:
 * - A PDF file named after the input, but with <tt>.root</tt>
 *   replaced with <tt>pdf</tt>
 * 
 */
class SummaryMCTrackDrawer : public SummaryDrawer
{
public:
  enum EFlags { 
    kEventInspector    = 0x001, 
    kTrackDensity      = 0x002, 
    kELossFits         = 0x004,
    kNormal            = 0x107
  };
  SummaryMCTrackDrawer() 
    : SummaryDrawer(),
      fSums(0),
      fResults(0)
  {}
  //__________________________________________________________________
  TFile* Init(const char* fname)
  {
    // --- Open the file ---------------------------------------------
    TString filename(fname);
    TFile*  file = TFile::Open(filename, "READ");
    if (!file) { 
      Error("Run", "Failed to open \"%s\"", filename.Data());
      return 0;
    }

    // --- Get top-level collection ----------------------------------
    fSums = GetCollection(file, "ForwardTracksSums");
    if (!fSums) {
      if (!fSums) return 0;
    }

    // --- Do the results ----------------------------------------------
    fResults = GetCollection(file, "ForwardTracksResults");
    if (!fResults) fResults = fSums; // Old-style

    return file;
  }
  //__________________________________________________________________
  /** 
   * 
   * 
   * @param fname 
   * @param what 
   */
  void Run(const char* fname, UShort_t what=kNormal)
  {
    // --- Initialize ------------------------------------------------
    TFile* file = 0;
    if (!(file = Init(fname))) return;

    // --- Make our canvas -------------------------------------------
    TString pdfName(fname);
    pdfName.ReplaceAll(".root", ".pdf");
    CreateCanvas(pdfName, what & kLandscape);
    DrawTitlePage();

    // --- Set pause flag --------------------------------------------
    fPause = what & kPause;

    // --- Do each sub-algorithm -------------------------------------
    if (what & kEventInspector)    DrawEventInspector(fSums);
    if (what & kTrackDensity)      DrawTrackDensity(fSums);
    if (what & kELossFits)         DrawELossFits(fResults);
    
    CloseCanvas();
  }
protected:
  //____________________________________________________________________
  void DrawBetaGammadEdx(TCollection* dir, 
			 const char* prefix,
			 UShort_t extra=0x1F)
  {
    Info("DrawBetaGammadEdx", "Drawing beta*gamma & dE/dx from %s", 
	 dir->GetName());
    fBody->Divide(2,2);
    
    TH2* betaGammadEdx = GetH2(dir, "betaGammadEdx");
    TH2* betaGammaEta  = GetH2(dir, "betaGammaEta");
    TH2* dEdxEta       = GetH2(dir, "dEdxEta");

    if (betaGammadEdx) {
      Int_t n = betaGammadEdx->GetEntries();
      betaGammadEdx->Scale(1./n, "width");
      UInt_t opt = kLogx|kLogy|kLogz|kGridx;
      DrawInPad(fBody, 1, betaGammadEdx, "colz", opt);
      TGraph* gs[] = { 0, 0, 0, 0, 0 }; 

      if (extra & 0x01) DrawInPad(fBody, 1, gs[0]=FromGFMATE(),     "C3 SAME");
      if (extra & 0x02) DrawInPad(fBody, 1, gs[1]=FromRPPFull(),    "C SAME");
      if (extra & 0x04) DrawInPad(fBody, 1, gs[2]=FromRPPNoDelta(), "C SAME");
      if (extra & 0x08) DrawInPad(fBody, 1, gs[3]=FromRPPNoRad(),   "C SAME");
      if (extra & 0x10) DrawInPad(fBody, 1, gs[4]=FromRPPMean(),    "C SAME");

      TArrow* mip = new TArrow(3.5, 1.1*betaGammadEdx->GetYaxis()->GetXmin(), 
			       3.5, .6, 0.02, "<|");
      mip->SetAngle(30);
      mip->SetLineWidth(2);
      mip->SetFillColor(mip->GetLineColor());
      fBody->cd(1);
      mip->Draw();

      TLegend* l = new TLegend(0.05, 0.05, 0.95, 0.95);
      l->SetBorderSize(0);
      l->SetFillColor(0);
      l->SetFillStyle(0);
      // l->SetNColumns(2);
      // l->SetNDC();
      for (Int_t i = 0; i < 5; i++) { 
	if (!gs[i]) continue; 
	l->AddEntry(gs[i], gs[i]->GetTitle(), "l");
      }
      fBody->cd(4);
      l->Draw();
    }
    if (betaGammaEta) {
      Int_t n = betaGammaEta->GetEntries();
      betaGammaEta->Scale(1./n, "width");      
      DrawInPad(fBody, 2, betaGammaEta, "colz", kLogy|kLogz|kGridx);
    }
    if (dEdxEta) {
      Int_t n = dEdxEta->GetEntries();
      dEdxEta->Scale(1. / n, "width");      
      DrawInPad(fBody, 3, dEdxEta,      "colz", kLogy|kLogz|kGridy);
    }

    TString tit(prefix);
    if (!tit.IsNull()) tit.Append(" - ");
    tit.Append("BetaGamma and dE/dx");
    PrintCanvas(tit);
  }
  //____________________________________________________________________
  void DrawTrackDensity(TCollection* parent)
  {
    const char* folderName = "mcTrackELoss";
    SummaryDrawer::DrawTrackDensity(parent, folderName);

    TCollection* mc = GetCollection(parent, folderName, false);
    if (!mc) return; 
    DrawBetaGammadEdx(mc, "");
  }
  //____________________________________________________________________
  void DrawELossFits(TCollection* parent)
  {
    Info("DrawELossFits", "Drawing energy loss fits");
    TCollection* ef = GetCollection(parent, "fmdEnergyFitter");
    if (!ef) return ; 

    
    const Char_t** ringNames = GetRingNames();
    const Char_t** ringPtr   = ringNames;
    while ((*ringPtr)) { 
      TCollection* ringCol = GetCollection(ef, *ringPtr);
      if (!ringCol) { 
	ringPtr++;
	continue; 
      }
      
      DrawBetaGammadEdx(ringCol, *ringPtr);
      Info("DrawELossFits", "Drawing energy loss fits for %s", *ringPtr);

      TCollection* all  = GetCollection(ringCol, "elossDists");
      TCollection* prim = GetCollection(ringCol, "primaryDists");
      TCollection* sec  = GetCollection(ringCol, "secondaryDists");
      if (!all || !prim || !sec) {
	ringPtr++;
	continue; 
      }
      
      TIter next(all);
      TH1*  allHist = 0;
      Int_t iPad    = 0;
      Int_t nPad    = 2;
      Int_t nTotal  = 0;
      while ((allHist = static_cast<TH1*>(next()))) {
	TH1* primHist = static_cast<TH1*>(prim->FindObject(allHist->GetName()));
	TH1* secHist  = static_cast<TH1*>(sec ->FindObject(allHist->GetName()));
	if (!primHist  || !secHist) continue; 

	allHist ->GetXaxis()->SetRangeUser(0,2);
	primHist->SetMarkerStyle(24);
	secHist ->SetMarkerStyle(25);
	
	if (iPad == 0) fBody->Divide(2,1);
	iPad++;

	DrawInPad(fBody, iPad, allHist,  "",     kLogy);
	DrawInPad(fBody, iPad, primHist, "same", kLogy);
	DrawInPad(fBody, iPad, secHist,  "same", kLogy|kLegend);
	nTotal++;

	if (iPad == nPad) { 
	  PrintCanvas(Form("%s Delta fits - page %d", 
			   *ringPtr, iPad/nPad));
	  
	  iPad = 0;
	}
      }
      Info("DrawELossFits", "Drew %d energy loss fits for %s", 
	   nTotal, *ringPtr);
      
      ringPtr++;
    }
  }
  //____________________________________________________________________
  void DrawTitlePage()
  {
    fBody->cd();

    TLatex* ltx = new TLatex(.5, .7, "MC track inspection");
    ltx->SetNDC();
    ltx->SetTextSize(0.07);
    ltx->SetTextAlign(22);
    ltx->Draw();

    TCollection* fwd = fSums; // GetCollection(f, "ForwardSums");
    Double_t y = .6;
    
    Double_t save = fParName->GetTextSize();
    fParName->SetTextSize(0.03);
    fParVal->SetTextSize(0.03);

    TCollection* ei = GetCollection(fwd, "fmdEventInspector");
    if (ei) { 

      UShort_t sys=0, sNN=0;
      Int_t field=0;
      ULong_t runNo=0;
      Bool_t mc=false;
      GetParameter(ei, "sys", sys);
      GetParameter(ei, "sNN", sNN);
      GetParameter(ei, "field", field);
      GetParameter(ei, "runNo", runNo);
      if (!GetParameter(ei, "mc", mc, false)) mc = false;
      
      TString sysString;    SysString(sys, sysString);
      TString sNNString;    SNNString(sNN, sNNString);
      
      DrawParameter(y, "System", sysString);
      DrawParameter(y, "#sqrt{s_{NN}}", sNNString);
      DrawParameter(y, "L3 B field", Form("%+2dkG", field));
      DrawParameter(y, "Run #", Form("%6lu", runNo));
      DrawParameter(y, "Simulation", (mc ? "yes" : "no"));	
    }
    PrintCanvas("Title page");
    fParName->SetTextSize(save);
    fParVal->SetTextSize(save);
  }
  //__________________________________________________________________
  /** 
   * @{ 
   * @name Bethe-bloc 
   */
  /** Scale a graph by density (multiply) and mass (divide). 
   * 
   * @param graph Graph to scale 
   *
   * @param density If @c true, scale by the Si density
   * (@f$ \rho=2.33/\mbox{cm}^3@f$).  The y axis is assumed to have units of
   * @f$ MeVg^{-1}cm^2@f$.
   * 
   * @param mass Mass to scale with. The x axis is assumed to be the
   * kinetic energy of the particles in units of GeV.  
  */
  void ScaleGraph(TGraph* graph, bool density=true, double mass=1) 
  {
      Double_t*      x   = graph->GetX();
      Double_t*      y   = graph->GetY();
      const Double_t rho = (density ? 2.33 : 1);
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
    return graph;
  }

  /**
   * Draw energy loss as obtained from GEANT 3.21 GFMATE. 
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
    // gre->Draw("c3 same");
    return gre;
  }

  /** 
   * Get the response functin @f$ f(\Delta_p/x)@f$ from Review of
   * Particle Physics (fig. 27.7).  It is scaled to the value at
   * MPV. 
   */ 
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
    }
    return graph;
  }

  /** 
   * Get the correction to Bethe-Bloc from Review of Particle Physics
   * (fig 27.8). 
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
  /* @} */

  TCollection* fSums;
  TCollection* fResults;
};

// #endif
