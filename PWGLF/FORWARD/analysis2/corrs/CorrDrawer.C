#ifndef __CINT__
# include "SummaryDrawer.C"
# include "AliFMDCorrAcceptance.h"
# include "AliFMDCorrSecondaryMap.h"
# include "AliFMDCorrELossFit.h"
# include "AliFMDCorrNoiseGain.h"
# include "AliForwardUtil.h"
# include "AliForwardCorrectionManager.h"
# include "AliLog.h"
# include <TString.h>
# include <TError.h>
#else
class SummaryDrawer;
class TObject;
// class TAxis;
class AliFMDCorrAcceptance;
class AliFMDCorrSecondaryMap;
class AliFMDCorrELossFit;
class AliFMDCorrNoiseGain;
#include <TString.h>
#endif

class CorrDrawer : public SummaryDrawer
{
public:
  TString  fELossExtra;
  UShort_t fMinQuality;
  /** 
   * Constructor 
   */
  CorrDrawer() 
  {
    fELossExtra = ""; // forward_eloss.root";
    fMinQuality = AliFMDCorrELossFit::kDefaultQuality;
  }
  /** 
   * Destructor.  Closes the PDF 
   */
  ~CorrDrawer() 
  {
    CloseCanvas();
  }
  /** 
   * Create output file name 
   * 
   * @param out    Output file name on return 
   * @param prefix Prefix of the file name 
   */
  static void MakeFileName(TString&        out,
			   const TString&  prefix)
  {
    out = TString::Format("forward_%s.pdf", prefix.Data());
  }

  /** 
   * Run the correction drawer, fetching information from extra file
   * 
   * @param what     What to draw
   * @param extra    Extra file 
   * @param options  Options
   * @param local    Local DB
   */
  void Run(const Char_t* what, 
	   const Char_t* extra, 
	   Option_t*     options="",
	   const Char_t* local="") 
  { 
    Run(AliForwardCorrectionManager::ParseFields(what), 
	extra, options, local);
  }
  /** 
   * Run the correction drawer, fetching information from extra file
   * 
   * @param what     What to draw
   * @param extra    Extra file 
   * @param options  Options
   * @param local    Local DB
   */
  void Run(UShort_t      what, 
	   const Char_t* extra, 
	   Option_t*     options="",
	   const Char_t* local="") 
  { 
    fELossExtra    = extra;
    ULong_t  runNo = 0;
    UShort_t sys   = 0;
    UShort_t sNN   = 0;
    Short_t  fld   = 0;
    Bool_t   mc    = false;
    Bool_t   sat   = false;
    if (!GetInformation(runNo, sys, sNN, fld, mc, sat)) return;
    
    Run(what, runNo, sys, sNN, fld, mc, sat, options, local);
  }

  /** 
   * Draw corrections using the correction manager to get them 
   *  
   * @param what    What to draw 
   * @param runNo   Run Number
   * @param sys     Collision system 
   * @param sNN     Center of mass energy 
   * @param field   L3 Field 
   * @param mc      Simulations or not
   * @param sat     Satellite interactions or not 
   * @param options Options
   * @param local   Local database file 
   */
  void Run(const Char_t* what, 
	   ULong_t       runNo, 
	   const Char_t* sys, 
	   UShort_t      sNN, 
	   UShort_t      field,
	   Bool_t        mc=false, 
	   Bool_t        sat=false,
	   Option_t*     options="",
	   const char*   local="")	   
  {
    Run(AliForwardCorrectionManager::ParseFields(what), 
	runNo, AliForwardUtil::ParseCollisionSystem(sys), 
	sNN, field, mc, sat, options, local);
  }
  /** 
   * Draw corrections using the correction manager to get them 
   *  
   * @param what    What to draw 
   * @param runNo   Run Number
   * @param sys     Collision system 
   * @param sNN     Center of mass energy 
   * @param field   L3 Field 
   * @param mc      Simulations or not
   * @param sat     Satellite interactions or not 
   * @param options Options
   * @param local   Local database file 
   */
  void Run(UShort_t    what, 
	   ULong_t     runNo, 
	   UShort_t    sys, 
	   UShort_t    sNN, 
	   UShort_t    field,
	   Bool_t      mc=false, 
	   Bool_t      sat=false,
	   Option_t*   options="",
	   const char* local="")
  {
    AliForwardCorrectionManager& mgr = AliForwardCorrectionManager::Instance();
    mgr.SetDebug(true);

    // Set local prefix 
    if (local) mgr.SetPrefix(gSystem->DirName(local));

    // Get output file name 
    TString  name;
    if (what & AliForwardCorrectionManager::kSecondaryMap)
      AppendName(name, AliForwardCorrectionManager::kSecondaryMap);
    if (what & AliForwardCorrectionManager::kAcceptance) 
      AppendName(name, AliForwardCorrectionManager::kAcceptance);
    if (what & AliForwardCorrectionManager::kELossFits) 
      AppendName(name, AliForwardCorrectionManager::kELossFits);
    if (what & AliForwardCorrectionManager::kNoiseGain) 
      AppendName(name, AliForwardCorrectionManager::kNoiseGain);
    if (what & AliForwardCorrectionManager::kVertexBias) 
      Warning("CorrDrawer","Vertex bias not implemented yet");
    if (what & AliForwardCorrectionManager::kDoubleHit) 
      Warning("CorrDrawer","Double hit not implemented yet");    
    if (what & AliForwardCorrectionManager::kMergingEfficiency) 
      Warning("CorrDrawer","Merging efficiency not implemented yet");

    // Filter the ones we can handle 
    UShort_t flags = what & (AliForwardCorrectionManager::kELossFits|
			     AliForwardCorrectionManager::kAcceptance|
			     AliForwardCorrectionManager::kSecondaryMap|
			     AliForwardCorrectionManager::kNoiseGain);
    if (!mgr.Init(runNo, sys, sNN, field, mc, sat, flags, true)) {
      Error("CorrDrawer", "Failed to initialize for flags=0x%02x"
		"run=%lu, sys=%hu, sNN=%hu, field=%hd, mc=%d, sat=%d",
		flags, runNo, sys, sNN, field, mc, sat);
      return;
    }

    TString out;
    MakeFileName(out, name); // , runNo, sys, sNN, field, mc, sat);

    TString opts(options);
    opts.ToUpper();
    Bool_t landscape = opts.Contains("LANDSCAPE");
    Bool_t few       = opts.Contains("FEW");
    Bool_t details   = !opts.Contains("SINGLE");
    if (opts.Contains("PORTRAIT")) landscape = false;
    CreateCanvas(out, landscape);
    

    fBody->cd();
    if (details) {
      Double_t y = .8;
      DrawParameter(y, "Run #", Form("%lu", runNo));
      DrawParameter(y, "System", AliForwardUtil::CollisionSystemString(sys));
      DrawParameter(y, "#sqrt{s_{NN}}", 
		    AliForwardUtil::CenterOfMassEnergyString(sNN));
      DrawParameter(y, "L3 field", AliForwardUtil::MagneticFieldString(field));
      DrawParameter(y, "Simulation", Form("%s", mc ? "yes" : "no"));
      DrawParameter(y, "Satellite", Form("%s", sat ? "yes" : "no"));
      PrintCanvas("Title");
    }
      
    if (what & AliForwardCorrectionManager::kSecondaryMap) 
      DrawIt(mgr.GetSecondaryMap(), details);
    if (what & AliForwardCorrectionManager::kAcceptance) 
      DrawIt(mgr.GetAcceptance(), details);
    if (what & AliForwardCorrectionManager::kELossFits)
      DrawIt(mgr.GetELossFit(), details, few);
    if (what & AliForwardCorrectionManager::kNoiseGain)
      DrawIt(mgr.GetNoiseGain(), details);

    // Done
    CloseCanvas();
  }
  /** 
   * Fall-back method
   * 
   * @param o Object to draw
   */
  virtual void Draw(const TObject* o) 
  {
    if (!o) return;
    Warning("CorrDrawer", "Don't know how to draw a %s object", 
	    o->ClassName());
  }
  /** 
   * Draw a single plot of the mean acceptance correction
   * 
   * @param acc Acceptance correction
   */
  virtual void Draw(const AliFMDCorrAcceptance* acc) { Summarize(acc, false); }
  /** 
   * Draw a single plot of the mean secondary correction 
   * 
   * @param sec Secondary correction
   */
  virtual void Draw(const AliFMDCorrSecondaryMap* sec) { Summarize(sec, false);}
  /** 
   * Draw a single plot summarizing the energy loss fits
   * 
   * @param fits Energy loss fits
   */
  virtual void Draw(const AliFMDCorrELossFit* fits) { Summarize(fits, false); }
  /** 
   * A generalized entry to the summarization functions
   * 
   * @param what    What to show - only one field
   * @param runNo   Run number
   * @param sys     System 
   * @param sNN     Center of mass energy in GeV
   * @param field   L3 magnetic field
   * @param mc      Simulation flag
   * @param sat     Satellite interaction flag
   * @param options Options 
   * @param local   Local storage
   *
   * @deprecated See Run instead 
   */
  virtual void Summarize(const TString& what, 
			 ULong_t        runNo, 
			 const Char_t*  sys, 
			 UShort_t       sNN, 
			 Short_t        field,
			 Bool_t         mc=false, 
			 Bool_t         sat=false,
			 Option_t*      options="",
			 const char*    local="")
  {
    Summarize(AliForwardCorrectionManager::ParseFields(what), 
	      runNo, AliForwardUtil::ParseCollisionSystem(sys), 
	      sNN, field, mc, sat, options, local);
  }
  /** 
   * A generalized entry to the summarization functions
   * 
   * @param what    What to show - only one field
   * @param runNo   Run number
   * @param sys     System 
   * @param sNN     Center of mass energy in GeV
   * @param field   L3 magnetic field
   * @param mc      Simulation flag
   * @param sat     Satellite interaction flag
   * @param options Options 
   * @param local   Local storage
   *
   * @deprecated See Run instead 
   */
  virtual void Summarize(UShort_t    what, 
			 ULong_t     runNo, 
			 UShort_t    sys, 
			 UShort_t    sNN, 
			 Short_t     field,
			 Bool_t      mc=false, 
			 Bool_t      sat=false,
			 Option_t*   options="",
			 const char* local="")
  {
    Run(what, runNo, sys, sNN, field, mc, sat, options, local);
  }
  /** 
   * Fall-back method
   * 
   * @param o Object to draw
   * @param pdf Not used
   */
  virtual void Summarize(const TObject* o, Bool_t pdf=true) 
  {
    if (!o) return;
    Warning("CorrDrawer", "Don't know how to draw a %s object (PDF: %s)", 
	    o->ClassName(), pdf ? "yes" : "no");
  }
  /** 
   * Draw a single summary plot or multiple plots of the acceptance
   * correction. A new Canvas is created for this.
   * 
   * @param acc Acceptance correction
   * @param pdf If true, do multiple plots. Otherwise a single summary plot
   */
  virtual void Summarize(const AliFMDCorrAcceptance* acc, Bool_t pdf=true) 
  { 
    CreateCanvas(CanvasName("acceptance.pdf"), false, pdf);
    DrawIt(acc, pdf); 
    if (pdf) CloseCanvas();
  }
  /** 
   * Draw a single summary plot multiple plots of the secondary
   * correction. A new canvas is created for this.
   * 
   * @param sec Secondary correction
   * @param pdf If true, do multiple plots. Otherwise a single summary plot
   */
  virtual void Summarize(const AliFMDCorrSecondaryMap* sec, Bool_t pdf=true) 
  { 
    CreateCanvas(CanvasName("secondarymap.pdf"), false, pdf);
    DrawIt(sec, pdf); 
    if (pdf) CloseCanvas();
  }
  /** 
   * Draw a single summary plot multiple plots of the energy loss
   * fits.  A new canvas is created for this.
   * 
   * @param corr Energy loss fits
   * @param pdf  If true, do multiple plots. Otherwise a single summary plot
   */
  virtual void Summarize(const AliFMDCorrNoiseGain* corr, Bool_t pdf=true) 
  { 
    CreateCanvas(CanvasName("noisegain.pdf"), true, pdf);
    DrawIt(corr, pdf); 
    if (pdf) CloseCanvas();
  }
  /** 
   * Draw a single summary plot multiple plots of the energy loss
   * fits.  A new canvas is created for this.
   * 
   * @param fits Energy loss fits
   * @param pdf  If true, do multiple plots. Otherwise a single summary plot
   */
  virtual void Summarize(const AliFMDCorrELossFit* fits, Bool_t pdf=true) 
  { 
    CreateCanvas(CanvasName("elossfits.pdf"), true, pdf);
    DrawIt(fits, pdf); 
    if (pdf) CloseCanvas();
  }
  /** 
   * Draw a single summary plot/multiple plots of the correction.
   * A new canvas is created for this.
   * 
   * @param what     What to plot
   * @param output   Output of correction pass (must exist)
   * @param local    Local storage of correction
   * @param options  Various options
   *
   * @deprecated Use Run instead 
   */
  static void Summarize(const TString& what   = "", 
			Bool_t         /*mc*/ = false,
			const TString& output = "",
			const TString& local  = "fmd_corrections.root",
			Option_t*      options= "")
  {
    CorrDrawer* drawer = new CorrDrawer;
    drawer->Run(AliForwardCorrectionManager::ParseFields(what), 
		output, local, options);
  }
  /** 
   * Draw a single summary plot/multiple plots of the correction.
   * A new canvas is created for this.
   * 
   * @param what     What to plot
   * @param output   Output of correction pass (must exist)
   * @param local    Local storage of correction
   * @param options  Various options
   *
   * @deprecated Use Run instead 
   */
  static void Summarize(UShort_t       what, 
			Bool_t         /*mc*/ = false,
			const TString& output = "",
			const TString& local  = "fmd_corrections.root",
			Option_t*      options= "")
  {
    CorrDrawer* drawer = new CorrDrawer;
    drawer->Run(what, output, options, local);
  }
protected:
  /** 
   * Append a name to output prefix 
   * 
   * @param what  What to append to 
   * @param which Which string to append
   */
  void AppendName(TString& what, UShort_t which)
  {
    if (!what.IsNull()) what.Append("_");
    switch (which) {
    case AliForwardCorrectionManager::kSecondaryMap:
      what.Append("secondary"); break;
    case AliForwardCorrectionManager::kAcceptance:
      what.Append("acceptance"); break;
    case AliForwardCorrectionManager::kELossFits:		  
      what.Append("elossfits"); break;
    case AliForwardCorrectionManager::kNoiseGain:		  
      what.Append("noisegain"); break;
    default:
      what.Append("unknown"); break;
    }
  }

  /** 
   * Get information from auxillary file 
   * 
   * @param runNo  On return, the run number
   * @param sys    On return, the collision system
   * @param sNN    On return, the collision energy 
   * @param fld    On return, the L3 magnetic field   
   * @param mc     On return, true for MC input 
   * @param sat    On return, true for satellite input enabled
   * 
   * @return true on success, false otherwise 
   */
  virtual Bool_t GetInformation(ULong_t&  runNo,
				UShort_t& sys,
				UShort_t& sNN,
				Short_t&  fld, 
				Bool_t&   mc, 
				Bool_t&   sat)
  { 
    TFile* fout = TFile::Open(fELossExtra, "READ");
    if (!fout) { 
      Warning("SummarizeELoss", "Correction task output \"%s\" not found",
	      fELossExtra.Data());
      return false;
    }
    Bool_t ret = false;
    try {
      TCollection* forward = GetCollection(fout, "ForwardELossSums");
      if (!forward) {
	forward = GetCollection(fout, "forwardQAResults");
	if (!forward) throw false;
      }
      
      TCollection* eventInsp = GetCollection(forward, "fmdEventInspector");
      if (!eventInsp) throw false;

      if (!GetParameter(eventInsp, "sys",       sys))   throw false;
      if (!GetParameter(eventInsp, "sNN",       sNN))   throw false;
      if (!GetParameter(eventInsp, "field",     fld))   throw false;
      if (!GetParameter(eventInsp, "satellite", sat))   throw false;
      if (!GetParameter(eventInsp, "runNo",     runNo)) throw false;
      if (!GetParameter(eventInsp, "mc",        mc))    throw false;

      ret = true;
    }
    catch (bool e) {
      ret = e;
    }
    if (fout) fout->Close();
    return ret;
  }
  /** 
   * Get the canvas name.  If the auxillary file has been set, use
   * that as the base of the canvas name.  Otherwise use @a def.
   * 
   * @param def Default value 
   * 
   * @return Canvas name 
   */
  virtual TString CanvasName(const char* def) 
  { 
    TString canName(def);
    if (!fELossExtra.IsNull()) { 
      canName = gSystem->BaseName(fELossExtra.Data());
      canName.ReplaceAll(".root", ".pdf");
    }
    return canName;
  }
  /** 
   * Fall-back method
   * 
   * @param o Object to summarize
   */
  virtual void DrawIt(const TObject* o) 
  {
    if (!o) return;
    Warning("CorrDrawer", "Don't know how to summarize a %s object", 
	    o->ClassName());
  }
  /** 
   * Draw the acceptance correction 
   * 
   * @param corr    Correction
   * @param details If true, make a multipage PDF, otherwise plot the mean. 
   */
  virtual void DrawIt(const AliFMDCorrAcceptance* corr, Bool_t details=true)
  {
    if (!corr) {
      Warning("CorrDrawer","No acceptance available");
      return;
    }

    if (!fCanvas) {
      Warning("CorrDrawer", "No canvas");
      return;
    }

    // --- Get vertex axis ---------------------------------------------
    const TAxis& vtxAxis = corr->GetVertexAxis();
    Int_t        nVtx    = vtxAxis.GetNbins();

    // --- Create stacks for summaries ---------------------------------
    TObjArray* stacks  = CreateVtxStacks(vtxAxis);
    TObjArray* stacks2 = (corr->HasOverflow() && details 
			  ? CreateVtxStacks(vtxAxis) : 0);
    
    //__________________________________________________________________
    // Create a title page 
    if (details) {
      fBody->cd();
      TLatex* ll = new TLatex(.5,.8, fCanvas->GetTitle());
      ll->SetTextAlign(22);
      ll->SetTextSize(0.03);
      ll->SetNDC();
      ll->Draw();
      
      TLatex* l = new TLatex(.5,.8, "");
      l->SetNDC();
      l->SetTextSize(0.03);
      l->SetTextFont(132);
      l->SetTextAlign(12);
      l->DrawLatex(0.2, 0.70, "Acceptance due to dead channels");
      l->SetTextAlign(22);
      l->DrawLatex(0.5, 0.55, "c_{v,r}(#eta,#phi) = #frac{"
		   "#sum active strips #in (#eta,#phi)}{"
		 "#sum strips #in (#eta,#phi)}");
      
      PrintCanvas("Acceptance");
    }
    
    // --- Loop over vertex ------------------------------------------
    for (UShort_t v=1; v <= nVtx; v++) { 
      Double_t vzMin = vtxAxis.GetBinLowEdge(v);
      Double_t vzMax = vtxAxis.GetBinUpEdge(v);

      if (details) DivideForRings(true, true);

      // --- Loop over detectors -------------------------------------
      for (UShort_t d = 1; d <= 3; d++) {
	UShort_t     nQ = (d == 1 ? 1 : 2);
	for (UShort_t q = 0; q < nQ; q++) { 
	  Char_t r = (q == 0 ? 'I' : 'O');
	  
	  TH2* h2 = corr->GetCorrection(d, r, v);
	  if (!h2) { 
	    Warning("DrawCorrAcc", "No correction for FMD%d%c, v=%d", d, r, v);
	    corr->ls();
	    continue;
	  }
	  
	  if (details) DrawInRingPad(d, r, h2, "colz");

	  Int_t nY = h2->GetNbinsY();
	  TH1* hh = h2->ProjectionX(Form("FMD%d%c", d, r), 1, nY);
	  hh->Scale(1. / nY);
	  hh->SetDirectory(0);
	  hh->SetMarkerColor(AliForwardUtil::RingColor(d, r));
	  hh->SetLineColor(AliForwardUtil::RingColor(d, r));
	  hh->SetFillColor(AliForwardUtil::RingColor(d, r));
	  hh->SetFillStyle(3001);
	  
	  THStack* stack = static_cast<THStack*>(stacks->At(v-1));
	  if (!stack) { 
	    Error("", "No stack at v=%d", v-1);
	    continue;
	  }
	  stack->Add(hh);

	  if (!stacks2) {
	    Warning("", "No phi acceptance defined");
	    continue;
	  }
	  stack = static_cast<THStack*>(stacks2->At(v-1));
	  if (!stack) { 
	    Error("", "No stack at v=%d", v-1);
	    continue;
	  }
	  TH1* hp = corr->GetPhiAcceptance(d, r, v);
	  if (!hp) { 
	    Error("", "No phi acceptance at v=%d", v-1);
	    continue;
	  }
	  hp->SetDirectory(0);
	  hp->SetMarkerColor(AliForwardUtil::RingColor(d, r));
	  hp->SetLineColor(AliForwardUtil::RingColor(d, r));
	  hp->SetFillColor(AliForwardUtil::RingColor(d, r));
	  hp->SetFillStyle(3001);
	  // Info("", "Adding phi acceptance plot %d", Int_t(hp->GetEntries()));
	  stack->Add(hp);

	}
      }
      if (details) 
	PrintCanvas(Form("%+5.1fcm<IP_{z}<%+5.1fcm", vzMin, vzMax));
    }
    if (DrawVtxStacks(stacks2, 1.2)) {
      PrintCanvas("#phi acceptance");
    }
    if (DrawVtxStacks(stacks, 1.2)) {
      PrintCanvas("#LTacceptance#GT");
    }
  }
  /** 
   * Draw the secondary correction 
   * 
   * @param corr       Correction
   * @param details If true, make a multipage PDF, otherwise plot the mean. 
   */
  virtual void DrawIt(const AliFMDCorrSecondaryMap* corr, bool details) 
  {
    if (!corr) {
      Warning("CorrDrawer","No secondary map available");
      return;
    }

    if (!fCanvas) {
      Warning("CorrDrawer", "No canvas");
      return;
    }
    
    const TAxis& vtxAxis = corr->GetVertexAxis();
    Int_t        nVtx    = vtxAxis.GetNbins();
    TObjArray*   stacks  = CreateVtxStacks(vtxAxis);
    
    //__________________________________________________________________
    // Create a title page 
    if (details) {
      fBody->cd();
      TLatex* ll = new TLatex(.5,.8, fCanvas->GetTitle());
      ll->SetTextAlign(22);
      ll->SetTextSize(0.03);
      ll->SetNDC();
      ll->Draw();
      
      TLatex* l = new TLatex(.5,.8, "");
      l->SetNDC();
      l->SetTextSize(0.03);
      l->SetTextFont(132);
      l->SetTextAlign(12);
      l->DrawLatex(0.2, 0.70, "Secondary map");
      l->SetTextAlign(22);
      l->DrawLatex(0.5, 0.60, "c_{v,r}(#eta,#phi)=#frac{"
		   "#sum N_{ch,primary,i}(#eta,#phi)}{"
		   "#sum N_{ch,FMD,i}(#eta,#phi)}");
      l->SetTextAlign(12);
      l->DrawLatex(0.2, 0.50, "N: Number of events");
      l->DrawLatex(0.2, 0.45, "N_{ch,primary,i}(#eta,#phi): Number of charged, "
		   "primary particles in (#eta,#phi) bin");
      l->DrawLatex(0.2, 0.40, "N_{ch,primary,i}(#eta,#phi): Number of charged, "
		   "particles that hit the FMD in (#eta,#phi) bin");
      l->DrawLatex(0.2, 0.35, "All quantities determined in MC");
      
      PrintCanvas("Secondary maps");
    }
    
    // --- Loop over vertex ------------------------------------------
    for (UShort_t v=1; v <= nVtx; v++) { 
      Double_t vzMin = vtxAxis.GetBinLowEdge(v);
      Double_t vzMax = vtxAxis.GetBinUpEdge(v);

      if (details) DivideForRings(true, true);

      // --- Loop over detectors -------------------------------------
      for (UShort_t d = 1; d <= 3; d++) {
	UShort_t     nQ = (d == 1 ? 1 : 2);
	for (UShort_t q = 0; q < nQ; q++) { 
	  Char_t r = (q == 0 ? 'I' : 'O');
	  
	  TH2* h2 = corr->GetCorrection(d, r, v);
	  if (!h2) { 
	    Warning("DrawCorrSec", "No correction for FMD%d%c, v=%d", d, r, v);
	    continue;
	  }
	  
	  if (details) DrawInRingPad(d, r, h2, "colz");

	  Int_t nY = h2->GetNbinsY();
	  TH1* hh = h2->ProjectionX(Form("FMD%d%c", d, r), 1, nY);
	  hh->Scale(1. / nY);
	  hh->SetDirectory(0);
	  hh->SetMarkerColor(AliForwardUtil::RingColor(d, r));
	  hh->SetLineColor(AliForwardUtil::RingColor(d, r));
	  hh->SetFillColor(AliForwardUtil::RingColor(d, r));
	  hh->SetFillStyle(3001);
	  
	  THStack* stack = static_cast<THStack*>(stacks->At(v-1));
	  if (!stack) { 
	    Error("", "No stack at v=%d", v-1);
	    continue;
	  }
	  stack->Add(hh);
	}
      }
      if (details) 
	PrintCanvas(Form("%+5.1fcm<IP_{z}<%+5.1fcm", vzMin, vzMax));
    }
    if (DrawVtxStacks(stacks, 3.5)) {
      PrintCanvas("#LTsecondary map#GT");
    }
  }
  virtual void DrawIt(const AliFMDCorrNoiseGain* corr, bool /*details*/) 
  {
    if (!corr) {
      Warning("CorrDrawer","No noise-gain correction available");
      return;
    }

    if (!fCanvas) {
      Warning("CorrDrawer", "No canvas");
      return;
    }

    DivideForRings(false,false);
    
    for (UShort_t d = 1; d <= 3; d++) { 
      UShort_t nQ = d == 1 ? 1 : 2;
      for (UShort_t q = 0; q < nQ; q++) { 
	Char_t   r  = q == 0 ? 'I' : 'O';
	UShort_t nS = q == 0 ?  20 :  40;
	UShort_t nT = q == 0 ? 512 : 256;
	
	TH2* h = new TH2D(Form("fmd%d%c", d, r),
			  Form("FMD%d%c", d, r), 
			  nT, -.5, nT-.5,  nS, -.5, nS-.5);
	h->SetDirectory(0);
	h->SetXTitle("Strip");
	h->SetYTitle("Sector");

	for (UShort_t s = 0; s < nS; s++) { 
	  for (UShort_t t = 0; t < nT; t++) { 
	    Float_t c = corr->Get(d,r,s,t);
	    h->Fill(t,s,c);
	  }
	}
	h->GetZaxis()->SetRangeUser(0,0.05);
	DrawInRingPad(d,r,h,"COLZ");
      }
    }
    PrintCanvas("Noise correction");
  }
  /** 
   * Draw the energy loss fits correction 
   * 
   * @param corr       Correction
   * @param details If true, make a multipage PDF, 
   *                   otherwise plot the parameters.
   * @param few     Only a few
   */
  virtual void DrawIt(const AliFMDCorrELossFit* corr, bool details,
		      bool few=true) 
  {
    if (!corr) {
      Warning("CorrDrawer","No energy loss fits available");
      return;
    }

    if (!fCanvas) {
      Warning("CorrDrawer", "No canvas");
      return;
    }

    AliFMDCorrELossFit* fits = const_cast<AliFMDCorrELossFit*>(corr);
    fits->CacheBins(8);
    fits->Print("C");
    TList* fitter = 0;
    if (details) { 
      TFile* hists = 0;
      TDirectory* savDir = gDirectory;
      if (!gSystem->AccessPathName(fELossExtra.Data())) {
	hists = TFile::Open(fELossExtra, "READ");
	Info("", "Opened %s -> %p", fELossExtra.Data(), hists);
      }
      else 
	Warning("", "Couldn't open %s", fELossExtra.Data());
      if (hists) {
	TList* fr = static_cast<TList*>(hists->Get("ForwardELossResults"));
	if (!fr)
	  fr = static_cast<TList*>(hists->Get("forwardQAResults"));
	// Info("", "Got forward results -> %p", fr);
	if (fr) { 
	  fitter = static_cast<TList*>(fr->FindObject("fmdEnergyFitter"));
	  // Info("", "Got fitter -> %p", fitter);
	  // fr->ls();
	  DrawEventInspector(fr);
	}
	hists->Close();
	savDir->cd();
      }
      fBody->cd();
      TLatex* ll = new TLatex(.5,.9, "ESD #rightarrow #Delta-fits"
			      /* fCanvas->GetTitle() */);
      ll->SetTextAlign(22);
      ll->SetTextSize(0.05);
      ll->SetNDC();
      ll->Draw();

      const Double_t fontSize = 0.03;
#define DL(X,Y,T) do { l->DrawLatex(X,Y,T); Y -= fontSize; } while (false)
      TLatex* l = new TLatex(.5,.8, "");
      l->SetNDC();
      l->SetTextSize(fontSize);
      l->SetTextFont(132);
      l->SetTextAlign(12);
      Double_t y = 0.80;
      Double_t x = 0.20;
      Double_t z = 0.30;
      DL(x,y,"1^{st} page is a summary of fit parameters");
      DL(x,y,"2^{nd} page is a summary of relative errors");
      DL(x,y,"Subsequent pages shows the fitted functions");
      y -= 0.01;
      DL(z,y,"Black line is the full fitted function");
      DL(z,y,"Coloured lines are the individual N-mip comp.");
      DL(x,y,"Each component has the form");
      y -= 0.02;
      DL(z,y,"f_{n}(x; #Delta, #xi, #sigma') = "
	 "#int_{-#infty}^{+#infty}dx' "
	 "landau(x'; #Delta, #xi)gaus(x'; x, #sigma')");
      y -= 0.02;
      DL(x,y,"The full function is given by");
      y -= 0.02;
      DL(z,y,"f_{N}(x; #Delta, #xi, #sigma', #bf{a}) = "
	 "C #sum_{i=1}^{N} a_{i} "
	 "f_{i}(x; #Delta_{i}, #xi_{i}, #sigma_{i}')");
      y -= 0.03;
      DL(z,y,"#Delta_{i} = i (#Delta_{1} + #xi_{1} log(i)) +#delta_{i}");
      DL(z,y,"#xi_{i} = i #xi_{1}");
      DL(z,y,"#sigma_{i} = #sqrt{i} #sigma_{1}");
      DL(z,y,"#sigma_{n} #dot{=} 0");
      DL(z,y,"#sigma_{i}'^{2} = #sigma^{2}_{n} + #sigma_{i}^{2}");
      DL(z,y,"#delta_{i} = c#sigmau/(1+1/i)^{pu#sqrt{u}}");
      DL(z,y,"u = #sigma/#xi");
      DL(z,y,"a_{1} #dot{=} 1");
      y -= 0.02;
      DL(z,y,Form("Least quality: %d", fMinQuality));
      y -= 0.02;
      if (fitter) {
	TObject* refit = fitter->FindObject("refitted");
	if (refit) DL(z,y, "Refitted distributions");//.10
      }
      PrintCanvas("Energy loss fits");
    }
    
    if (details && fitter) { 
      // Draw parameter from the fitter 
      fBody->cd();
      Double_t y = 0.90;
      Double_t s = fParName->GetTextSize();
      Double_t t = fParVal->GetTextSize();
      fParName->SetTextSize(0.04);
      fParVal->SetTextSize(0.04);
      DrawTParameter<double>(y, fitter, "lowCut");
      DrawTParameter<int>   (y, fitter, "nParticles");
      DrawTParameter<int>   (y, fitter, "minEntries");
      DrawTParameter<int>   (y, fitter, "subtractBins");
      DrawTParameter<bool>  (y, fitter, "doFits");
      DrawTParameter<double>(y, fitter, "maxE");
      DrawTParameter<int>   (y, fitter, "nEbins");
      DrawTParameter<bool>  (y, fitter, "increasingBins");
      DrawTParameter<double>(y, fitter, "maxRelPerError");
      DrawTParameter<double>(y, fitter, "maxChi2PerNDF");
      DrawTParameter<double>(y, fitter, "minWeight");
      DrawTParameter<double>(y, fitter, "regCut");
      DrawParameter(y,"Use #delta#Delta(#sigma/#xi)",
		    Form("%s", 
			 fits->TestBit(AliFMDCorrELossFit::kHasShift) 
			 ? "yes" : "no"));
      PrintCanvas("Fitter settings");
      fParName->SetTextSize(s);
      fParVal->SetTextSize(t);
    }

    fBody->cd();
    fits->Draw("error good");
    PrintCanvas("Fit overview");
    if (!details) return;

    //__________________________________________________________________
    // Draw relative parameter errors 
    fBody->cd();
    fits->Draw("relative good");
    PrintCanvas("Relative parameter errors");

    //__________________________________________________________________
    // Draw all fits individually
    for (UShort_t d=1; d<=3; d++) { 
      UShort_t nQ = (d == 1 ? 1 : 2);
      for (UShort_t q = 0; q < nQ; q++) { 
	Char_t r = (q == 0 ? 'I' : 'O');
	TList* dists = 0;
	TList* resis = 0;
	if (fitter) { 
	  // Info("", "Fitter: %s", fitter->GetName());
	  TList* dl = 
	    static_cast<TList*>(fitter->FindObject(Form("FMD%d%c",d,r)));
	  // Info("", "Got detector list -> %p", dl);
	  if (dl) { 
	    // Info("", "Detector list: %s", dl->GetName());
	    dists = static_cast<TList*>(dl->FindObject("elossDists"));
	    // Info("", "Got distributions -> %p", dists);
	    resis = static_cast<TList*>(dl->FindObject("elossResiduals"));
	    // Info("", "Got residuals -> %p", resis);
	  }
	}
	
	printf("FMD%d%c ", d, r);
	ClearCanvas();
	TObjArray*  ra = fits->GetRingArray(d, r);
	if (!ra) continue;
	DrawELossFits(d, r, ra, dists, resis, few);
      }
    }
  }
  /** 
   * CINT does too much when optimizing on a loop, so we take this out
   * to force CINT to not optimize the third nested loop.
   * 
   * @param d     Detector
   * @param r     Ring 
   * @param ra    Ring array 
   * @param dists Distributions (optional)
   * @param resis Residuals (optional)   
   * @param few   A few
   */
  void DrawELossFits(UShort_t d, Char_t r, TObjArray* ra, 
		     TList* dists, TList* resis, bool few)
  {
    Int_t nRow = 3;
    Int_t nCol = (few ? 1 : 2);
    Int_t nPad = nRow * nCol;
    AliFMDCorrELossFit::ELossFit* fit = 0;
    TIter next(ra);
    Int_t i = 0;
    Int_t j = 0;
    DividedPad divided(fBody, fLandscape, nCol, nRow);
    while ((fit = static_cast<AliFMDCorrELossFit::ELossFit*>(next()))) {
      j           = i % nPad;
      Bool_t last = j == nPad-1;
      if (j == 0) divided.Divide(true, true);

      Bool_t        same    = false;
      TVirtualPad*  drawPad = divided.GetPad(j); // fBody->GetPad(j+1);
      Int_t         subPad  = 0;
      // fBody->ls();
      // Info("", "Now in sub-pad %d of %d: %p", j, nPad, drawPad);
      // Info("", "Pad %s", drawPad->GetName());
      if (dists) { 
	// Info("", "Distributions: %s", dists->GetName());
	TString hName(Form("FMD%d%c_etabin%03d", d,r,fit->GetBin()));
	TH1* dist = static_cast<TH1*>(dists->FindObject(hName));
	TH1* resi = 0;
	if (resis) resi = static_cast<TH1*>(resis->FindObject(hName));
	// Info("", "Got histogram -> %p", dist);
	if (resi) { 
	  Bool_t err = resi->GetUniqueID() <= 1;
	  drawPad->SetGridx();
	  if (err) {
	    resi->SetYTitle("#chi^{2}_{bin}=(h-f)^{2}/#delta^{2}h");
	    for (Int_t k=1; k<=resi->GetNbinsX(); k++) { 
	      Double_t c = resi->GetBinContent(k);
	      Double_t e = resi->GetBinError(k);
	      if (e <= 0) continue;
	      c          *= c;
	      c          /= (e*e);
	      resi->SetBinContent(k, c);
	      resi->SetBinError(k, 0);
	    }
	  }
	  drawPad->Divide(1,2,0,0);
	  DrawInPad(drawPad, 2, resi, "HIST", kGridx);
	  subPad = 1;
	  Double_t red = fit->GetNu() > 0 ? fit->GetChi2() / fit->GetNu() : 0;
	  if (red > 0) {
	    drawPad->cd(2);
	    TLine* l = new TLine(resi->GetXaxis()->GetXmin(), red,
				 resi->GetXaxis()->GetXmax(), red);
	    l->SetLineWidth(2);
	    l->SetLineStyle(2);
	    l->Draw();
	    TLatex* cltx = new TLatex(0.5, 0.5,
				      Form("#chi^{2}/#nu=%6.2f", red));
	    cltx->SetNDC();
	    cltx->SetTextAlign(22);
	    cltx->SetTextFont(42);
	    cltx->SetTextSize(0.07);
	    cltx->Draw();
	    cltx->DrawLatex(0.5,0.4,Form("%g", dist->GetEntries()));
	  }
	}
	if (dist) { 
	  // Info("", "Histogram: %s", dist->GetName());
	  dist->SetFillStyle(3001);
	  dist->SetMarkerStyle(0);
	  DrawInPad(drawPad, subPad, dist, "HIST E", (subPad * kGridx) + kLogy);
	  same = true;	  
	}
      }
      // if (same)
      DrawInPad(drawPad, subPad, fit, 
		Form("comp good values legend peak %s", (same ? "same" : "")),
		kLogy);
      if (fit->GetQuality() < fMinQuality) { 
	TLatex* ltx = new TLatex(.2, .2, "NOT USED");
	ltx->SetNDC();
	ltx->SetTextFont(62);
	ltx->SetTextColor(kRed+1);
	ltx->SetTextAngle(30);
	ltx->SetTextSize(0.2);
	DrawInPad(fBody, j+1, ltx, "", 0);
	// ltx->Draw();
      }

      // else 
      // DrawInPad(fBody, j+1, fit, "comp good values legend", kLogy);
      printf(".");
      
      if (last) 
	PrintCanvas(Form("FMD%d%c page %d", d, r, (i/nPad)+1));
      i++;
    }
    j = i % nPad;
    if (j != 0) 
      PrintCanvas(Form("FMD%d%c page %d", d, r, (i/nPad)+1));
    printf(" done\n");
  }

  /** 
   * Create an array of per-vertex bin stacks 
   * 
   * @param vtxAxis Vertex axis 
   * 
   * @return Array of stacks 
   */
  TObjArray* CreateVtxStacks(const TAxis& vtxAxis)
  {
    // --- Create stacks for summaries ---------------------------------
    Int_t      nVtx    = vtxAxis.GetNbins();
    TObjArray* stacks  = new TObjArray(nVtx);
    for (UShort_t v = 1; v <= nVtx; v++) { 
      THStack* stack = new THStack(Form("vtx%02d", v),
				   Form("%+5.1f<v_{z}<%+5.1f",
					vtxAxis.GetBinLowEdge(v),
					vtxAxis.GetBinUpEdge(v)));
      stacks->AddAt(stack, v-1);
    }
    return stacks;
  }
  /** 
   * Draw the vertex stacks in the canvas 
   * 
   * @param stacks Stacks to draw
   * @param max    Possible maximum of the stacks 
   * 
   * @return true on success 
   */
  Bool_t DrawVtxStacks(TObjArray* stacks, Double_t max=-1)
  {
    if (!stacks) return false;
    // --- Make summary page -------------------------------------------
    Int_t nVtx = stacks->GetEntries();

    fBody->Divide(3, (nVtx+2)/3, 0, 0);
    Int_t ipad = 0;
    for (UShort_t v = 1; v <= nVtx; v++) {
      ipad++;
    
      if (nVtx == 10 && (ipad == 1 || ipad == 12)) ipad++;

      THStack*     stack = static_cast<THStack*>(stacks->At(v-1));
      if (!stack) { 
	Error("", "No stack at v=%d", v-1);
	continue;
      }
      TVirtualPad* pad   = fBody->cd(ipad);
      if (!pad) { 
	Error("", "No pad at %d", ipad);
	continue;
      }
      pad->SetFillColor(kWhite);
    
      if (max > 0) stack->SetMaximum(max);
      stack->Draw("nostack hist");
    }
    return true;
  }

};


  
//
// EOF
//
