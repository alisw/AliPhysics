#ifndef __CINT__
# include "SummaryDrawer.C"
# include "AliFMDCorrAcceptance.h"
# include "AliFMDCorrSecondaryMap.h"
# include "AliFMDCorrELossFit.h"
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
    fELossExtra = "forward_eloss.root";
    fMinQuality = 8;
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
  /*
   * @param runNo  Run Number
   * @param sys    Collision system 
   * @param sNN    Center of mass energy 
   * @param field  L3 Field 
   * @param mc     Simulations or not
   * @param sat    Satellite interactions or not 
			   ULong_t         runNo, 
			   UShort_t        sys, 
			   UShort_t        sNN, 
			   Short_t         field,
			   Bool_t          mc=false, 
			   Bool_t          sat=false) */
  {
    out = TString::Format("forward_%s.pdf", prefix.Data());
#if 0
    out = TString::Format("%s_run%09lu_%s_%04dGeV_%c%dkG_%s_%s.pdf",
			  prefix.Data(), runNo, 
			  (sys == 1 ? "pp" : 
			   sys == 2 ? "PbPb" : 
			   sys == 3 ? "pPb" : "unknown"), sNN, 
			  (field >= 0 ? 'p' : 'm'), TMath::Abs(field), 
			  (mc ? "MC" : "real"),
			  (sat ? "satellite" : "nominal"));
#endif
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
    default:
      what.Append("unknown"); break;
    }
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
   * @param local   Local database file 
   */
  void Run(UShort_t    what, 
	   ULong_t     runNo, 
	   UShort_t    sys, 
	   UShort_t    sNN, 
	   UShort_t    field,
	   Bool_t      mc=false, 
	   Bool_t      sat=false,
	   Option_t*   /*options*/="",
	   const char* local="")
  {
    AliForwardCorrectionManager& mgr = AliForwardCorrectionManager::Instance();
    mgr.SetDebug(true);
    UShort_t flags = 0;

    TString name;
    if (what & AliForwardCorrectionManager::kSecondaryMap) {
      flags |= AliForwardCorrectionManager::kSecondaryMap;
      if (local) mgr.SetSecondaryMapPath(local);
      AppendName(name, AliForwardCorrectionManager::kSecondaryMap);
    }
    if (what & AliForwardCorrectionManager::kAcceptance) {
      flags |= AliForwardCorrectionManager::kAcceptance;
      if (local) mgr.SetAcceptancePath(local);
      AppendName(name, AliForwardCorrectionManager::kAcceptance);
    }
    if (what & AliForwardCorrectionManager::kELossFits) {
      flags |= AliForwardCorrectionManager::kELossFits;
      if (local) mgr.SetELossFitsPath(local);
      AppendName(name, AliForwardCorrectionManager::kELossFits);
    }
    if (what & AliForwardCorrectionManager::kVertexBias) 
      Warning("CorrDrawer","Vertex bias not implemented yet");
    if (what & AliForwardCorrectionManager::kDoubleHit) 
      Warning("CorrDrawer","Double hit not implemented yet");    
    if (what & AliForwardCorrectionManager::kMergingEfficiency) 
      Warning("CorrDrawer","Merging efficiency not implemented yet");
    
    if (!mgr.Init(runNo, sys, sNN, field, mc, sat, flags, true)) {
      Error("CorrDrawer", "Failed to initialize for flags=0x%02x"
		"run=%lu, sys=%hu, sNN=%hu, field=%hd, mc=%d, sat=%d",
		flags, runNo, sys, sNN, field, mc, sat);
      return;
    }

    TString out;
    MakeFileName(out, name); // , runNo, sys, sNN, field, mc, sat);
    CreateCanvas(out);

    fBody->cd();
    Double_t y = .8;
    DrawParameter(y, "Run #", Form("%lu", runNo));
    DrawParameter(y, "System", AliForwardUtil::CollisionSystemString(sys));
    DrawParameter(y, "#sqrt{s_{NN}}", 
		  AliForwardUtil::CenterOfMassEnergyString(sNN));
    DrawParameter(y, "L3 field", AliForwardUtil::MagneticFieldString(field));
    DrawParameter(y, "Simulation", Form("%s", mc ? "yes" : "no"));
    DrawParameter(y, "Satellite", Form("%s", sat ? "yes" : "no"));
    PrintCanvas("Title");
      
    if (what & AliForwardCorrectionManager::kSecondaryMap) {
      const AliFMDCorrSecondaryMap* sec = mgr.GetSecondaryMap();
      if (!sec) 
	Warning("CorrDrawer","No secondary map available");
      else 
	DrawIt(sec, true);
    }
    if (what & AliForwardCorrectionManager::kAcceptance) {
      const AliFMDCorrAcceptance* acc = mgr.GetAcceptance();
      if (!acc) 
	Warning("CorrDrawer","No acceptance available");
      else 
	DrawIt(acc, true);
    }
    if (what & AliForwardCorrectionManager::kELossFits) {
      const AliFMDCorrELossFit* fit = mgr.GetELossFit();
      if (!fit) 
	Warning("CorrDrawer","No energy loss fits available");
      else 
	DrawIt(fit, true);
    }
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
   * @param local   Local storage
   */
  virtual void Summarize(UShort_t    what, 
			 ULong_t     runNo, 
			 UShort_t    sys, 
			 UShort_t    sNN, 
			 Short_t     field,
			 Bool_t      mc=false, 
			 Bool_t      sat=false,
			 Option_t*   /*options*/="",
			 const char* local="")
  {
    AliForwardCorrectionManager& mgr = AliForwardCorrectionManager::Instance();
    mgr.SetDebug(true);
    if (local) mgr.SetPrefix(gSystem->DirName(local));
    UShort_t flag = 0;
    
    if (what & AliForwardCorrectionManager::kSecondaryMap) 
      flag = AliForwardCorrectionManager::kSecondaryMap;
    if (what & AliForwardCorrectionManager::kAcceptance) 
      flag = AliForwardCorrectionManager::kAcceptance;
    if (what & AliForwardCorrectionManager::kELossFits) 
      flag = AliForwardCorrectionManager::kELossFits;
    if (what & AliForwardCorrectionManager::kVertexBias) 
      Warning("CorrDrawer","Vertex bias not implemented yet");
    if (what & AliForwardCorrectionManager::kDoubleHit) 
      Warning("CorrDrawer","Double hit not implemented yet");    
    if (what & AliForwardCorrectionManager::kMergingEfficiency) 
      Warning("CorrDrawer","Merging efficiency not implemented yet");
    if (flag == 0) { 
      Warning("CorrDrawer", "Nothing to draw");
      return;
    }
    
    if (!mgr.Init(runNo, sys, sNN, field, mc, sat, flag, true)) {
      Error("CorrDrawer", "Failed to initialize for flags=0x%02x "
	    "run=%lu, sys=%hu, sNN=%hu, field=%hd, mc=%d, sat=%d",
	    flag, runNo, sys, sNN, field, mc, sat);
      return;
    }

    TString prefix;
    if      (flag == AliForwardCorrectionManager::kSecondaryMap) 
      prefix = "secondarymap";
    else if (flag == AliForwardCorrectionManager::kAcceptance)
      prefix = "acceptance";
    else if (flag == AliForwardCorrectionManager::kELossFits) 
      prefix = "elossfits";
    else 
      prefix = "unknown";
    TString out;
    MakeFileName(out, prefix); // , runNo, sys, sNN, field, mc, sat);
    CreateCanvas(out);

    fBody->cd();
    Double_t y = .8;
    DrawParameter(y, "Run #", Form("%lu", runNo));
    DrawParameter(y, "System", AliForwardUtil::CollisionSystemString(sys));
    DrawParameter(y, "#sqrt{s_{NN}}", 
		  AliForwardUtil::CenterOfMassEnergyString(sys));
    DrawParameter(y, "L3 field", AliForwardUtil::MagneticFieldString(field));
    DrawParameter(y, "Simulation", Form("%s", mc ? "yes" : "no"));
    DrawParameter(y, "Satellite", Form("%s", sat ? "yes" : "no"));
    PrintCanvas("Title");
      
    if (flag == AliForwardCorrectionManager::kSecondaryMap) {
      const AliFMDCorrSecondaryMap* sec = mgr.GetSecondaryMap();
      if (!sec) 
	Warning("CorrDrawer","No secondary map available");
      else 
	DrawIt(sec, true);
    }
    else if (flag == AliForwardCorrectionManager::kAcceptance) {
      const AliFMDCorrAcceptance* acc = mgr.GetAcceptance();
      if (!acc) 
	Warning("CorrDrawer","No acceptance available");
      else 
	DrawIt(acc, true);
    }
    if (flag == AliForwardCorrectionManager::kELossFits) {
      const AliFMDCorrELossFit* fit = mgr.GetELossFit();
      if (!fit) 
	Warning("CorrDrawer","No energy loss fits available");
      else 
	DrawIt(fit, true);
    }
    CloseCanvas();
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
    CreateCanvas("acceptance.pdf", false, pdf);
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
    CreateCanvas("secondarymap.pdf", false, pdf);
    DrawIt(sec, pdf); 
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
    CreateCanvas("elossfits.pdf", false, pdf);
    DrawIt(fits, pdf); 
    if (pdf) CloseCanvas();
  }

  static void Summarize(const TString& what   = TString(""), 
			Bool_t         mc     = false,
			const TString& output = TString("forward_eloss.root"), 
			const TString& local  = TString("fmd_corrections.root"),
			Option_t*      options= "")
  {
    Summarize(AliForwardCorrectionManager::ParseFields(what), mc, 
	      output, local, options);
  }
  static void Summarize(UShort_t       what, 
			Bool_t         mc     = false,
			const TString& output = "forward_eloss.root", 
			const TString& local  = "fmd_corrections.root",
			Option_t*      options= "")
  {
    TFile* fout = TFile::Open(output, "READ");
    if (!fout) { 
      Warning("SummarizeELoss", "Energy loss task output %s not found",
	      output.Data());
      return;
    }
    TCollection* forward = GetCollection(fout, "ForwardELossSums");
    if (!forward) return;
    
    TCollection* eventInsp = GetCollection(forward, "fmdEventInspector");
    if (!eventInsp) return;

    UShort_t sys   = 0, sNN = 0;
    Int_t    field = 0;
    ULong_t  runNo = 0; 
    Bool_t satellite;
    if (!GetParameter(eventInsp, "sys",       sys))       return;
    if (!GetParameter(eventInsp, "sNN",       sNN))       return;
    if (!GetParameter(eventInsp, "field",     field))     return;
    if (!GetParameter(eventInsp, "satellite", satellite)) return;
    if (!GetParameter(eventInsp, "runNo",     runNo))     return;
    
    CorrDrawer* drawer = new CorrDrawer;
    drawer->Run(what, runNo, sys, sNN, field, mc, satellite,
		options, local);
  }
protected:
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
    if (!corr || !fCanvas) return;

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
    if (!corr || !fCanvas) return;
    
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
  /** 
   * Draw the energy loss fits correction 
   * 
   * @param corr       Correction
   * @param details If true, make a multipage PDF, 
   *                   otherwise plot the parameters. 
   */
  virtual void DrawIt(const AliFMDCorrELossFit* corr, bool details) 
  {
    if (!corr || !fCanvas) return;

    AliFMDCorrELossFit* fits = const_cast<AliFMDCorrELossFit*>(corr);
    fits->CacheBins(8);
    fits->Print("C");
    TList* fitter = 0;
    if (details) { 
      TFile* hists = 0;
      TDirectory* savDir = gDirectory;
      if (!gSystem->AccessPathName(fELossExtra.Data())) {
	hists = TFile::Open(fELossExtra, "READ");
	Info("", "Opened forward_eloss.root -> %p", hists);
      }
      else 
	Warning("", "Couldn't open %s", fELossExtra.Data());
      if (hists) {
	TList* fr = static_cast<TList*>(hists->Get("ForwardELossResults"));
	// Info("", "Got forward results -> %p", fr);
	if (fr) { 
	  fitter = static_cast<TList*>(fr->FindObject("fmdEnergyFitter"));
	  // Info("", "Got fitter -> %p", fitter);
	}
	hists->Close();
	savDir->cd();
      }
      fBody->cd();
      TLatex* ll = new TLatex(.5,.8, "ESD #rightarrow #Delta-fits"
			      /* fCanvas->GetTitle() */);
      ll->SetTextAlign(22);
      ll->SetTextSize(0.05);
      ll->SetNDC();
      ll->Draw();
      
      TLatex* l = new TLatex(.5,.8, "");
      l->SetNDC();
      l->SetTextSize(0.03);
      l->SetTextFont(132);
      l->SetTextAlign(12);
      l->DrawLatex(0.2, 0.70, "1^{st} page is a summary of fit parameters");
      l->DrawLatex(0.2, 0.67, "2^{nd} page is a summary of relative errors");
      l->DrawLatex(0.2, 0.64, "Subsequent pages shows the fitted functions");
      l->DrawLatex(0.3, 0.60, "Black line is the full fitted function");
      l->DrawLatex(0.3, 0.57, "Coloured lines are the individual N-mip comp.");
      //l->DrawLatex(0.3, 0.54, "Full drawn lines correspond to used components");
      //l->DrawLatex(0.3, 0.51, "Dashed lines correspond to ignored components");
      l->DrawLatex(0.2, 0.54, "Each component has the form");
      l->DrawLatex(0.3, 0.49, "f_{n}(x; #Delta, #xi, #sigma') = "
		   "#int_{-#infty}^{+#infty}d#Delta' "
		   "landau(x; #Delta', #xi)gaus(#Delta'; #Delta, #sigma')");
      l->DrawLatex(0.2, 0.44, "The full function is given by");
      l->DrawLatex(0.3, 0.41, "f_{N}(x; #Delta, #xi, #sigma', #bf{a}) = "
		 "C #sum_{i=1}^{N} a_{i} "
		   "f_{i}(x; #Delta_{i}, #xi_{i}, #sigma_{i}')");
      l->DrawLatex(0.3, 0.35, "#Delta_{i} = i (#Delta_{1} + #xi_{1} log(i))");
      l->DrawLatex(0.3, 0.32, "#xi_{i} = i #xi_{1}");
      l->DrawLatex(0.3, 0.29, "#sigma_{i} = #sqrt{i} #sigma_{1}");
      l->DrawLatex(0.3, 0.26, "#sigma_{n} #dot{=} 0");
      l->DrawLatex(0.3, 0.23, "#sigma_{i}'^{2} = #sigma^{2}_{n} + #sigma_{i}^{2}");
      l->DrawLatex(0.3, 0.20, "a_{1} #dot{=} 1");
      l->DrawLatex(0.3, 0.15, Form("Least quality: %d", fMinQuality));
      if (fitter) {
	TObject* refit = fitter->FindObject("refitted");
	if (refit) l->DrawLatex(0.3, .10, "Refitted distributions");
      }
      PrintCanvas("Energy loss fits");
    }

    fBody->cd();
    fits->Draw("error");
    PrintCanvas("Fit overview");
    if (!details) return;

    //__________________________________________________________________
    // Draw relative parameter errors 
    fBody->cd();
    fits->Draw("relative");
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
	DrawELossFits(d, r, ra, dists, resis);
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
   */
  void DrawELossFits(UShort_t d, Char_t r, TObjArray* ra, 
		     TList* dists, TList* resis)
  {
    Int_t nPad = 6;
    AliFMDCorrELossFit::ELossFit* fit = 0;
    TIter next(ra);
    Int_t i = 0;
    Int_t j = 0;
    while ((fit = static_cast<AliFMDCorrELossFit::ELossFit*>(next()))) {
      j           = i % nPad;
      Bool_t last = j == nPad-1;
      if (j == 0) DivideForRings(true, true);

      Bool_t        same    = false;
      TVirtualPad*  drawPad = fBody->GetPad(j+1);
      Int_t         subPad  = 0;
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
		Form("comp good values legend %s", (same ? "same" : "")),
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
    Int_t nVtx = 10; // stacks->GetEntries();

    fBody->Divide(3, (nVtx+2)/3, 0, 0);
    Int_t ipad = 0;
    for (UShort_t v = 1; v <= nVtx; v++) {
      ipad++;
    
      if (ipad == 1 || ipad == 12) ipad++;

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
