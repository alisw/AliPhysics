// ELossTimeTrain.C
#ifndef NO_TRAIN
#ifndef __CINT__
#include <AliAnalysisManager.h>
#include <fstream>
#else
class AliAnalysisManager;
#endif
#include "TrainSetup.C"
#include "ParUtilities.C"

/** 
 * Train to record time of each event 
 * 
 * @ingroup pwglf_forward_eventtime
 */        
class ELossTimeTrain : public TrainSetup
{
public:
  /** 
   * Constructor 
   * 
   * @param name The name 
   */
  ELossTimeTrain(const char* name="eventTime") : TrainSetup(name)
  { 
    fOptions.Add("map", "FILE", "File containg map", "map.root");
    fOptions.Set("type", "ESD");
  }
  /** 
   * Create our tasks 
   */
  void CreateTasks(AliAnalysisManager*)
  {
    if (!fRailway->LoadLibrary("PWGLFforward2")) 
      Fatal("CreateTasks", "Failed to load PWGLFforward2");
    
    if (!ParUtilities::MakeScriptPAR(fRailway->Mode() == Railway::kLocal,
				     "EventTimeTask.C",
				     // Gui because of CDB - sigh!
				     // XMLParser because of CDB 
				     "Gui,XMLParser,"
				     "STEERBase,CDB,ESD,AOD,ANALYSIS,OADB,"
				     "ANALYSISalice",
				     fRailway)) 
      Fatal("","Failed to make support PAR");
    if (!fRailway->LoadLibrary("EventTimeTask")) 
      Fatal("CreateTasks", "Failed to load EventTimeTask");

    if (!ParUtilities::MakeScriptPAR(fRailway->Mode() == Railway::kLocal,
				     "ELossTimeTask.C",
				     "Gui,STEERBase,CDB,ESD,AOD,ANALYSIS,OADB,"
				     "ANALYSISalice,PWGLFforward2,"
				     "EventTimeTask",
				     fRailway)) 
      Fatal("","Failed to make PAR");
    if (!fRailway->LoadLibrary("ELossTimeTask")) 
      Fatal("CreateTasks", "Failed to load ELossTimeTask");

    TString mapfile = fOptions.Get("map");
    gROOT->ProcessLine(Form("ELossTimeTask::Create(\"%s\")", mapfile.Data()));

    fRailway->LoadAux(mapfile.Data(), true);
  }
  /** 
   * Do not create a physics selection
   */
  // void CreatePhysicsSelection(Bool_t, AliAnalysisManager*) {}
  /** 
   * Do not create a centrality selection
   */
  void CreateCentralitySelection(Bool_t) {}
  /** 
   * Do not create an output handler
   */
  AliVEventHandler* CreateOutputHandler(UShort_t) { return 0; }
  /** 
   * The train class name 
   * 
   * @return Train class name
   */
  const char* ClassName() const { return "ELossTimeTrain"; }
  /** 
   * Overloaded to create new dNdeta.C and dndeta.sh in the output 
   * directory
   * 
   * @param asShellScript 
   */
  void SaveSetup(Bool_t asShellScript)
  {
    TrainSetup::SaveSetup(asShellScript);
    
    SaveDraw();
  }
  void SaveDraw()
  {
    std::ofstream o("Draw.C");
    o << "// Written by " << ClassName() << "\n"
      << "void Draw(const char* fileName=\"AnalysisResults.root\")\n"
      << "{\n"
      << "  gSystem->AddIncludePath(\"-DNO_TRAIN -DSUMMARY\");\n"
      << "  const char* fwd = \"$ALICE_PHYSICS/PWGLF/FORWARD/analysis2\";\n"
      << "  gSystem->AddIncludePath(Form(\"-I%s/scripts\", fwd));\n"
      << "  gROOT->SetMacroPath(Form(\"%s/trains:%s\", fwd,\n"
      << "                           gROOT->GetMacroPath()));\n"
      << "  gROOT->LoadMacro(\"ELossTimeTrain.C+\");\n"
      << "  ELossTimeSummary s;\n"
      << "  s.Run(fileName);\n"
      << "}\n"
      << std::endl;
    o.close();
  }
  void PostShellCode(std::ostream& f)
  {
    f << "  echo \"=== Draw results ...\"\n"
      << "  aliroot -l -b -q ${prefix}Draw.C\\(\\\"AnalysisResults.root\\\"\\)\n"
      << std::endl;
  }
};
#endif
#ifdef SUMMARY
# include <SummaryDrawer.C>
# include <TColor.h>

/**
 * Draw summary of the above train
 * 
 * @ingroup pwglf_forward_eventtime
 */
struct ELossTimeSummary : public SummaryDrawer 
{
  enum EFlags { 
    kEventInspector    = 0x001, 
  };
  /** 
   * Run the class 
   * 
   * @param fname Filename 
   * @param flags Flags
   */
  void Run(const char* fname, UShort_t flags=0x01)
  {
    // --- Open the file -----------------------------------------------
    TString filename(fname);
    TFile* file = TFile::Open(filename.Data(), "READ");
    if (!file) { 
      Error("Run", "Failed to open \"%s\"", filename.Data());
      return;
    }
    fPause         = flags & kPause;
    
    // --- Make our canvas ---------------------------------------------
    TString pdfName(filename);
    pdfName.ReplaceAll(".root", ".pdf");
    CreateCanvas(pdfName, flags & kLandscape);

    // --- Make title page -------------------------------------------
    TCollection* c = GetCollection(file, "elossTimeSums");
    DrawTitlePage(c);

    if (flags & kEventInspector)    DrawEventInspector(c);
    
    TH1* dt = GetH1(c, "dt");
    DrawInPad(fBody, 0, dt, "", kLogy);
    PrintCanvas("#Deltat");

    const char* rings[] = { "FMD1i", "FMD2i", "FMD2o", "FMD3o", "FMD3i", 0 };
    const char** pring   = rings;

    while (*pring) {
      DrawRing(c, *pring, dt);
      pring++;
    }

    CloseCanvas();
  }
  /** 
   * Draw the title page 
   * 
   * @param c Parent collection
   */
  void DrawTitlePage(TCollection* c)
  {
    fBody->cd();

    Double_t y = .7;
    TLatex* ltx = new TLatex(.5, y, "#Deltat vs #Delta/#Delta_{mip}");
    ltx->SetTextSize(0.07);
    ltx->SetTextFont(62);
    ltx->SetTextAlign(22);
    ltx->SetNDC();
    ltx->Draw();

    Double_t save = fParName->GetTextSize();
    fParName->SetTextSize(0.03);
    fParVal->SetTextSize(0.03);
    y = .6;
    
    TCollection* fc = GetCollection(c, "fmdEventInspector");
    UShort_t sys;
    UShort_t sNN;
    ULong_t  runNo;
    GetParameter(fc, "sys",     sys);
    GetParameter(fc, "sNN",     sNN);
    GetParameter(fc, "runNo",   runNo);
    
    DrawParameter(y, "Run #", Form("%lu", runNo));
    TString tS; SysString(sys, tS);      DrawParameter(y, "System", tS);
    TString tE; SNNString(sNN, tE);      DrawParameter(y, "#sqrt{s_{NN}}", tE);

    PrintCanvas("Title page");
    fParName->SetTextSize(save);
    fParVal->SetTextSize(save);
  }
  /** 
   * Draw a single ring
   * 
   * @param c     Parent collection
   * @param ring  Ring name 
   * @param dt    Histogram of delta time 
   */
  void DrawRing(TCollection* c, const char* ring, TH1* dt)
  {
    TCollection* lring = GetCollection(c, ring);
    if (!lring) return;
    
    TH2* h2 = GetH2(lring, "dtVsELoss");
    if (!h2) return;

    THStack* stack = new THStack(ring, ring);
    // stack->SetTitle(ring);

    THStack* ratios = new THStack(Form("Ratios for %s",ring), ring);
    // stack->SetTitle(ring);

    Printf(ring);
    Int_t j = 2;
    TH1* first = 0;
    Double_t lfirst = 0;
    Double_t rmax = 0;
    Double_t max = 3;
    for (Int_t i = 1; i <= h2->GetNbinsY(); i++) { 
      TH1*     h     = h2->ProjectionX(Form("%s_%03d", ring, i), i,i);
      Double_t logDt = h2->GetYaxis()->GetBinCenter(i);

      Int_t    nFill = h->GetEntries();
      if (nFill <= 1000) continue;
      Double_t norm = dt->GetBinContent(i);
      if (norm <= 1e-6) {
	Warning("", "Normalization=%f<1e-6 but got "
		"%d>1000 entries for log10(dt)=%5.3f", norm, nFill, logDt);
	continue;
      }
      if (!first && logDt > TMath::Log10(25.)) {
	lfirst = logDt;
	first = h;
      }
      // Info("", "Normalization is %f", norm);
      h->Sumw2();
      h->Scale(1. / norm);
      h->SetTitle(Form("log_{10}(#Deltat)=%5.3f", logDt));

      Float_t r, g, b;
      TColor::HSV2RGB((j-1)*45, 1, .8, r, g, b);
      Int_t col = TColor::GetColor(r, g, b);
      j++;
      h->SetLineColor(col);
      h->SetLineStyle(j % 3+1);
      h->SetLineWidth(2);
      // h->SetFillColor(col);
      // h->SetFillStyle(3002);
      stack->Add(h);

      if (h == first) continue;
      TH1* rh = static_cast<TH1*>(h->Clone(Form("ratio%s", h->GetName())));
      // rh->SetTitle(Form("log_{10}(#Deltat)=%5.3f", logDt));
      rh->Divide(first);
      for (Int_t k = 1; k <= rh->GetNbinsX(); k++) {
	if (rh->GetXaxis()->GetBinCenter(k) > max) break;
	rmax = TMath::Max(rmax, rh->GetBinContent(k));
      }

      ratios->Add(rh);
    }
    Double_t savX = fParVal->GetX();
    Double_t savY = fParVal->GetY();
    fParVal->SetX(0.12);
    fParVal->SetY(0.12);
    fBody->Divide(1,2,0,0);
    DrawInPad(fBody,1,stack,"nostack hist", kLogy|kLegend);
    stack->GetXaxis()->SetRangeUser(-.1,max);
    stack->GetYaxis()->SetTitle("1/N_{ev} dN/d(#Delta/#Delta_{mip})");

    fParVal->SetX(0.6);
    fParVal->SetY(0.4);
    DrawInPad(fBody,2,ratios,"nostack hist", kLegend);
    ratios->GetXaxis()->SetRangeUser(-.1, max);
    ratios->GetXaxis()->SetTitle("#Delta/#Delta_{mip}");
    ratios->GetYaxis()
      ->SetTitle(Form("X/(1/N_{ev}dN/d(#Delta/#Delta_{mip}))|_{%5.3f}",lfirst));
    Printf("Max: %f (%f)", ratios->GetMaximum(), rmax);
    ratios->SetMaximum(rmax*1.2);


    fParVal->SetX(savX);
    fParVal->SetY(savY);
    PrintCanvas(ring);
  }
};

#endif 
// EOF
