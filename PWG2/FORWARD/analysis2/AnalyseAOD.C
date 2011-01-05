/**
 * @file 
 * 
 * @ingroup pwg2_forward_scripts
 */
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TFile.h>
#include <TTree.h>
#include <TError.h>
#include <TStyle.h>
#include <THStack.h>
#include <TLegend.h>
#include <TMath.h>
#include <TMultiGraph.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TSystem.h>
#include "AliAODForwardMult.h"
#include "OtherData.C"

/** 
 * Draw the data stored in the AOD 
 *
 * To use this, do 
 * 
 * @code 
 *   Root> .L $ALICE_ROOT/PWG2/FORWARD/analysis2/Compile.C
 *   Root> Compile("$ALICE_ROOT/PWG2/FORWARD/analysis2/AnalyseAOD.C")
 *   Root> AnalyseAOD dr
 *   Root> dr.Run("AliAODs.root",-10,10,5,AliAODForwardMult::kInel,900)
 * @endcode 
 * 
 * The output is stored in a ROOT file 
 * 
 * See also the script Pass2.C 
 * 
 * @ingroup pwg2_forward_scripts
 */
struct AnalyseAOD 
{
public: 
  /** AOD tree */
  TTree*             fTree;
  /** AOD object */
  AliAODForwardMult* fAOD;
  /** AOD object */
  AliAODForwardMult* fMCAOD;
  /** Output file */
  TFile*             fOut;
  /** Summed histogram */
  TH2D*              fSum;
  /** Summed histogram */
  TH2D*              fMCSum;
  /** Vertex efficiency */
  Double_t           fVtxEff;
  /** Title to put on the plot */
  TString            fTitle;
  /** Do HHD comparison */
  Bool_t             fDoHHD;

  //__________________________________________________________________
  /** 
   * Constructor 
   * 
   */
  AnalyseAOD()
    : fTree(0), 
      fAOD(0),
      fMCAOD(0),
      fOut(0), 
      fSum(0),
      fMCSum(0),
      fVtxEff(0),
      fTitle(""),
      fDoHHD(kTRUE)
  {}
  //__________________________________________________________________
  /** 
   * Reset internal structures 
   * 
   */
  void Clear(Option_t* )
  {
    if (fTree && fTree->GetCurrentFile()) { 
      fTree->GetCurrentFile()->Close();
      delete fTree;
    }
    if (fOut) {
      fOut->Close();
      delete fOut;
    }
    if (fSum)   delete fSum;
    if (fMCSum) delete fMCSum;
    fTree   = 0;
    fOut    = 0;
    fSum    = 0;
    fMCSum  = 0;
    fVtxEff = 0;
    
  }
  //__________________________________________________________________
  /** 
   * Run 
   * 
   * @param file   Input file with AOD tree
   * @param vzMin  Minimum interaction point z coordinate 
   * @param vzMax  Maximum interaction point z coordinate 
   * @param rebin  How many times to re-bin the @f$ dN_{ch}/d\eta@f$
   * @param mask   Trigger mask 
   * @param energy Collision energy @f$ \sqrt{s}@f$ or @f$\sqrt{s_{NN}}@f$ 
   * @param title  Title to put on the plot
   * @param doHHD  Whether to do HHD comparison
   * @param doComp Whether to do comparisons 
   * 
   * @return True on success, false otherwise 
   */
  Bool_t Run(const char* file="AliAODs.root", 
	     Double_t vzMin=-10, Double_t vzMax=10, Int_t rebin=1,
	     Int_t mask=AliAODForwardMult::kInel, Int_t energy=900,
	     const char* title="", Bool_t doHHD=false, Bool_t doComp=false)
  {
    TString trgName;
    if    (mask & AliAODForwardMult::kInel)    trgName.Append("inel-");
    if    (mask & AliAODForwardMult::kInelGt0) trgName.Append("inelgt0-");
    if    (mask & AliAODForwardMult::kNSD)     trgName.Append("nsd-");
    if (trgName.EndsWith("-")) trgName.Remove(trgName.Length()-1);
    if (trgName.IsNull()) trgName = "unknown";
    TString outName = 
      TString::Format("dndeta_%04dGeV_%c%02d-%c%02dcm_rb%02d_%s.root",
		      energy, 
		      vzMin < 0 ? 'm' : 'p', Int_t(TMath::Abs(vzMin)+.5), 
		      vzMax < 0 ? 'm' : 'p', Int_t(TMath::Abs(vzMax)+.5), 
		      rebin, trgName.Data());
    fTitle = title;
    if (!Open(file, outName)) return kFALSE;
    if (!Process(vzMin,vzMax,mask)) return kFALSE;
    if (!Finish(rebin, mask, energy,doHHD,doComp)) return kFALSE;

    return kTRUE;
  }
  //__________________________________________________________________
  /** 
   * Open the files @a fname containg the tree with AliAODEvent objects, 
   * and also open the output file @a outname for writting. 
   * 
   * @param fname   Input file name 
   * @param outname Output file name 
   * 
   * @return true on success, false otherwise 
   */
  Bool_t Open(const char* fname, const char* outname) 
  {
    Clear("");
    
    TFile* file = TFile::Open(fname, "READ");
    if (!file) { 
      Error("Open", "Cannot open file %s", fname);
      return kFALSE;
    }

    // Get the AOD tree 
    fTree = static_cast<TTree*>(file->Get("aodTree"));
    if (!fTree) {
      Error("Init", "Couldn't get the tree");
      return kFALSE;
    }

    // Set the branch pointer 
    fTree->SetBranchAddress("Forward", &fAOD);

    // Set the branch pointer 
    fTree->SetBranchAddress("ForwardMC", &fMCAOD);

    fOut = TFile::Open(outname, "RECREATE");
    if (!fOut) { 
      Error("Open", "Couldn't open %s", outname);
      return kFALSE;
    }
    return kTRUE;
  }
  //__________________________________________________________________
  /** 
   * Process the events 
   * 
   * @param vzMin  Minimum interaction point z coordinate 
   * @param vzMax  Maximum interaction point z coordinate 
   * @param mask   Trigger mask 
   *
   * @return true on success, false otherwise 
   */
  Bool_t Process(Double_t vzMin, Double_t vzMax, Int_t mask) 
  {
    Int_t nTriggered  = 0;                    // # of triggered ev.
    Int_t nWithVertex = 0;                    // # of ev. w/vertex
    Int_t nAccepted   = 0;                    // # of ev. used
    Int_t nAvailable  = fTree->GetEntries();  // How many entries
 
    for (int i = 0; i < nAvailable; i++) { 
      fTree->GetEntry(i);
      if (((i+1) % 100) == 0) {
	fprintf(stdout,"Event # %9d of %9d, %9d accepted so far\r", 
	       i+1, nAvailable, nAccepted);
	fflush(stdout);
      }

      // Create sum histogram on first event - to match binning to input
      if (!fSum) {
	fSum = static_cast<TH2D*>(fAOD->GetHistogram().Clone("d2ndetadphi"));
	Info("Process", "Created sum histogram (%d,%f,%f)x(%d,%f,%f)", 
	     fSum->GetNbinsX(), 
	     fSum->GetXaxis()->GetXmin(),
	     fSum->GetXaxis()->GetXmax(),
	     fSum->GetNbinsY(), 
	     fSum->GetYaxis()->GetXmin(),
	     fSum->GetYaxis()->GetXmax());
      }
      if (!fMCSum && fTree->GetBranch("ForwardMC")) { 
	fMCSum = 
	  static_cast<TH2D*>(fMCAOD->GetHistogram().Clone("d2ndetadphiMC"));
	Info("Process", "Created MC sum histogram (%d,%f,%f)x(%d,%f,%f)", 
	     fMCSum->GetNbinsX(), 
	     fMCSum->GetXaxis()->GetXmin(),
	     fMCSum->GetXaxis()->GetXmax(),
	     fMCSum->GetNbinsY(), 
	     fMCSum->GetYaxis()->GetXmin(),
	     fMCSum->GetYaxis()->GetXmax());
      }
 
      // fAOD->Print();
      // Other trigger/event requirements could be defined 
      if (!fAOD->IsTriggerBits(mask)) continue; 
      nTriggered++;

      // Check if there's a vertex 
      if (!fAOD->HasIpZ()) continue; 
      nWithVertex++;

      // Select vertex range (in centimeters) 
      if (!fAOD->InRange(vzMin, vzMax)) continue; 
      nAccepted++;
 
      // Add contribution from this event
      fSum->Add(&(fAOD->GetHistogram()));

      // Add contribution from this event
      if (fMCSum) fMCSum->Add(&(fMCAOD->GetHistogram()));
    }
    printf("\n");
    fVtxEff = Double_t(nWithVertex)/nTriggered;

    Info("Process", "Total of %9d events\n"
	 "               %9d has a trigger\n" 
	 "               %9d has a vertex\n" 
	 "               %9d was used\n",
	 nAvailable, nTriggered, nWithVertex, nAccepted);

    return kTRUE;
  }
  //__________________________________________________________________
  /** 
   * Finish the stuff and draw 
   * 
   * @param rebin  How many times to rebin
   * @param energy Collision energy 
   * @param mask   Trigger mask 
   * @param doHHD  Whether to do HHD comparison
   * @param doComp Whether to do comparisons 
   * 
   * @return true on success, false otherwise 
   */
  Bool_t Finish(Int_t rebin, Int_t mask, Int_t energy, 
		Bool_t doHHD, Bool_t doComp)
  {
    fOut->cd();

    // Get acceptance normalisation from underflow bins 
    TH1D* norm   = fSum->ProjectionX("norm", 0, 1, "");
    // Project onto eta axis - _ignoring_underflow_bins_!
    TH1D* dndeta = fSum->ProjectionX("dndeta", 1, -1, "e");
    dndeta->SetTitle("ALICE Forward");
    // Normalize to the acceptance 
    dndeta->Divide(norm);
    // Scale by the vertex efficiency 
    dndeta->Scale(fVtxEff, "width");
    dndeta->SetMarkerColor(kRed+1);
    dndeta->SetMarkerStyle(20);
    dndeta->SetMarkerSize(1);
    dndeta->SetFillStyle(0);
    Rebin(dndeta, rebin);

    TH1D* dndetaMC = 0;
    if (fMCSum) { 
      // Get acceptance normalisation from underflow bins 
      norm = fMCSum->ProjectionX("norm", 0, 1, "");
      // Project onto eta axis - _ignoring_underflow_bins_!
      dndetaMC = fMCSum->ProjectionX("dndetaMC", 1, -1, "e");
      dndetaMC->SetTitle("ALICE Forward (MC)");
      // Normalize to the acceptance 
      dndetaMC->Divide(norm);
      // Scale by the vertex efficiency 
      dndetaMC->Scale(fVtxEff, "width");
      dndetaMC->SetMarkerColor(kRed+3);
      dndetaMC->SetMarkerStyle(21);
      dndetaMC->SetMarkerSize(1);
      dndetaMC->SetFillStyle(0);
      Rebin(dndetaMC, rebin);
    }

    DrawIt(dndeta, dndetaMC, mask, energy, doHHD, doComp, rebin);

    return kTRUE;
  }
  //__________________________________________________________________
  /** 
   */
  void DrawIt(TH1* dndeta, TH1* dndetaMC, Int_t mask, Int_t energy, 
	      Bool_t doHHD, Bool_t doComp, Int_t rebin)
  {
    // --- 1st part - prepare data -----------------------------------
    TH1* sym = Symmetrice(dndeta);
    // Rebin(sym, rebin);

    Double_t max = dndeta->GetMaximum();

    // Make our histogram stack 
    THStack* stack = new THStack("results", "Results");

    // If available, plot the MC truth 
    TH1* mc = 0;
    TH1* mcsym = 0;
    if (!gSystem->AccessPathName("AnalysisResults.root")) { 
      TFile* anares = TFile::Open("AnalysisResults.root", "READ");
      if (anares) { 
	TList* list = static_cast<TList*>(anares->Get("Forward"));
	if (list) { 
	  mc = static_cast<TH1*>(list->FindObject("mcSumEta"));
	  Rebin(mc, rebin);
	  if (mc) {
	    mcsym = Symmetrice(mc);
	    stack->Add(mc);
	    stack->Add(mcsym);
	    max = TMath::Max(mc->GetMaximum(),max);
	  }
	}
	anares->Close();
	fOut->cd();
      }
    }

    // Get the data from HHD's analysis - if available 
    TH1* hhd    = 0;
    TH1* hhdsym = 0;
    Info("DrawIt", "Will %sdraw HHD results", (doHHD ? "" : "not "));
    if (doHHD) {
      TString hhdName(fOut->GetName());
      hhdName.ReplaceAll("dndeta", "hhd");    
      hhd    = GetHHD(hhdName.Data(), mask & AliAODForwardMult::kNSD);
      hhdsym = 0;
      if (hhd) { 
	// Symmetrice and add to stack 
	hhdsym = Symmetrice(hhd);
	stack->Add(hhd);
	stack->Add(hhdsym);
	max = TMath::Max(hhd->GetMaximum(),max);
      }
      else 
	Warning("DrawIt", "No HHD data found");
      fOut->cd();
    }

    // If we have MC analysed data, plot it 
    if (dndetaMC) { 
      Info("DrawIt", "Maximum of MC dndeta: %f, was %f", 
	   dndetaMC->GetMaximum(),dndeta->GetMaximum());
      TH1* mcSym = Symmetrice(dndetaMC);
      stack->Add(dndetaMC);
      stack->Add(mcSym);
      max = TMath::Max(dndetaMC->GetMaximum(),max);
    }

    // Add the analysis results to the list 
    stack->Add(dndeta);
    stack->Add(sym);

    // Get graph of 'other' data - e.g., UA5, CMS, ... - and check if
    // there's any graphs.  Set the pad division based on that.
    Info("DrawIt", "Will %sdraw other results", (doComp ? "" : "not "));
    TMultiGraph* other    = (doComp ? GetData(energy, mask) : 0);
    THStack*     ratios   = MakeRatios(dndeta, sym, hhd, hhdsym, 
				       mc, mcsym, other);

    // Check if we have ratios 

    // --- 2nd part - Draw in canvas ---------------------------------
    // Make a canvas 
    gStyle->SetOptTitle(0);
    gStyle->SetTitleFont(132, "xyz");
    gStyle->SetLabelFont(132, "xyz");
    TCanvas* c = new TCanvas("c", "C", 900, 800);
    c->SetFillColor(0);
    c->SetBorderMode(0);
    c->SetBorderSize(0);
    c->cd();

    Double_t yd = (ratios ? 0.3 : 0);

    // Make a sub-pad for the result itself
    TPad* p1 = new TPad("p1", "p1", 0, yd, 1.0, 1.0, 0, 0);
    p1->SetTopMargin(0.05);
    p1->SetBottomMargin(ratios ? 0.001 : 0.1);
    p1->SetRightMargin(0.05);
    p1->SetGridx();
    p1->SetTicks(1,1);
    p1->Draw();
    p1->cd();

    // Fix the apperance of the stack and redraw. 
    if (other) {
      TGraphAsymmErrors* o      = 0;
      TIter              nextG(other->GetListOfGraphs());
      while ((o = static_cast<TGraphAsymmErrors*>(nextG()))) {
	Double_t gmax = TMath::MaxElement(o->GetN(), o->GetY());
	Info("DrawIt", "Maximum of %s is %f, was %f", o->GetName(), gmax, max);
	max = TMath::Max(max, gmax);
      }
    }
    max *= 1.1;
    Info("DrawIt", "Setting maximum to %f", max);
    stack->SetMinimum(ratios ? -0.1 : 0);
    stack->SetMaximum(max);
    FixAxis(stack, 1/(1-yd)/1.5, "#frac{1}{N} #frac{dN_{ch}}{#eta}");
    p1->Clear();
    stack->DrawClone("nostack e1");

    // Draw other data 
    if (other) {
      TGraphAsymmErrors* o      = 0;
      TIter              nextG(other->GetListOfGraphs());
      while ((o = static_cast<TGraphAsymmErrors*>(nextG())))
	o->Draw("same p");
    }

    // Make a legend in the main result pad 
    TString trigs(AliAODForwardMult::GetTriggerString(mask));
    TLegend* l = p1->BuildLegend(.15,p1->GetBottomMargin()+.01,.90,.35);
    l->SetNColumns(2);
    l->SetFillColor(0);
    l->SetFillStyle(0);
    l->SetBorderSize(0);
    l->SetTextFont(132);
    p1->cd();

    // Put a title on top 
    TLatex* tit = new TLatex(0.10, 0.95, fTitle.Data());
    tit->SetNDC();
    tit->SetTextFont(132);
    tit->SetTextSize(0.05);
    tit->Draw();

    // Put a nice label in the plot 
    TLatex* tt = new TLatex(.93, .93, 
			    Form("#sqrt{s}=%dGeV, %s", energy,
				 AliAODForwardMult::GetTriggerString(mask)));
    tt->SetNDC();
    tt->SetTextFont(132);
    tt->SetTextAlign(33);
    tt->Draw();

    // Mark the plot as preliminary 
    TLatex* pt = new TLatex(.93, .88, "Preliminary");
    pt->SetNDC();
    pt->SetTextFont(22);
    pt->SetTextSize(0.07);
    pt->SetTextColor(kRed+1);
    pt->SetTextAlign(33);
    pt->Draw();
    c->cd();

    // If we do not have the ratios, fix up the display 
    // p1->SetPad(0, 0, 1, 1);
    // p1->SetBottomMargin(0.1);
    // l->SetY1(0.11);
    // stack->SetMinimum(0);
    // FixAxis(stack, (1-yd)/1,  "#frac{1}{N} #frac{dN_{ch}}{#eta}",10,false);
    if (ratios) {
      // If we do have the ratios, then make a new pad and draw the 
      // ratios there 
      c->cd();
      TPad* p2 = new TPad("p2", "p2", 0, 0.0, 1.0, yd, 0, 0);
      p2->SetTopMargin(0.001);
      p2->SetRightMargin(0.05);
      p2->SetBottomMargin(1/yd * 0.07);
      p2->SetGridx();
      p2->SetTicks(1,1);
      p2->Draw();
      p2->cd();

      // Fix up axis 
      FixAxis(ratios, 1/yd/1.5, "Ratios", 5);

      // Fix up y range and redraw 
      ratios->SetMinimum(.58);
      ratios->SetMaximum(1.22);
      p2->Clear();
      ratios->DrawClone("nostack e1");
      
      // Make a legend 
      TLegend* l2 = p2->BuildLegend(.15,p2->GetBottomMargin()+.01,.9,.6);
      l2->SetNColumns(2);
      l2->SetFillColor(0);
      l2->SetFillStyle(0);
      l2->SetBorderSize(0);
      l2->SetTextFont(132);

      // Make a nice band from 0.9 to 1.1 
      TGraphErrors* band = new TGraphErrors(2);
      band->SetPoint(0, sym->GetXaxis()->GetXmin(), 1);
      band->SetPoint(1, dndeta->GetXaxis()->GetXmax(), 1);
      band->SetPointError(0, 0, .1);
      band->SetPointError(1, 0, .1);
      band->SetFillColor(kYellow+2);
      band->SetFillStyle(3002);
      band->Draw("3 same");

      // Replot the ratios on top 
      ratios->DrawClone("nostack e1 same");

      c->cd();
    }
    
    // Plot to disk
    TString imgName(fOut->GetName());
    imgName.ReplaceAll(".root", ".png");
    c->SaveAs(imgName.Data());
    imgName.ReplaceAll(".png", ".C");
    c->SaveAs(imgName.Data());
    
    fOut->cd();
    stack->Write();
    if (other)  other->Write();
    if (ratios) ratios->Write();

    // Close our file 
    fOut->Close();
  }
  //__________________________________________________________________
  /** 
   * Get the result from previous analysis code 
   * 
   * @param fn  File to open 
   * @param nsd Whether this is NSD
   * 
   * @return null or result of previous analysis code 
   */
  TH1* GetHHD(const char* fn="fmd_dNdeta_mult.root", bool nsd=false) 
  {
    TDirectory* savdir = gDirectory;
    if (gSystem->AccessPathName(fn)) { 
      Warning("GetHHD", "Output of HHD analysis (%s) not available", fn);
      return 0;
    }
    TFile* file = TFile::Open(fn, "READ");
    if (!file) { 
      Warning("GetHHD", "couldn't open HHD file %s", fn);
      savdir->cd();
      return 0;
    }
    TString hist(Form("dNdeta_dNdeta%s", (nsd ? "NSD" : "")));
    TH1* h = static_cast<TH1*>(file->Get(hist.Data()));
    if (!h) { 
      Warning("GetHHD", "Couldn't find HHD histogram %s in %s", 
	      hist.Data(), fn);
      file->Close();
      savdir->cd();
      return 0;
    }
    TH1* r = static_cast<TH1*>(h->Clone("dndeta_hhd"));
    r->SetTitle("1/N dN_{ch}/d#eta (HHD)");
    r->SetFillStyle(0);
    r->SetFillColor(0);
    r->SetMarkerStyle(21);
    r->SetMarkerColor(kPink+1);
    r->SetDirectory(savdir);

    file->Close();
    savdir->cd();
    return r;
  }
  //__________________________________________________________________
  /** 
   */ 
  THStack* MakeRatios(const TH1* dndeta, const TH1* sym, 
		      const TH1* hhd,    const TH1* hhdsym, 
		      const TH1* mc,     const TH1* mcsym,
		      TMultiGraph* other) const 
  {
    // If we have 'other' data, then do the ratio of the results to that
    Bool_t hasOther = (other && other->GetListOfGraphs() && 
		       other->GetListOfGraphs()->GetEntries() > 0);
    Bool_t hasHhd   = (hhd && hhdsym);
    if (!hasOther && !hasHhd && !mc && !mcsym) return 0;

    THStack* ratios = new THStack("ratios", "Ratios");
    if (hasOther) {
      TGraphAsymmErrors* o      = 0;
      TIter              nextG(other->GetListOfGraphs());
      while ((o = static_cast<TGraphAsymmErrors*>(nextG()))) {
	ratios->Add(Ratio(dndeta, o));
	ratios->Add(Ratio(sym,    o));
	ratios->Add(Ratio(hhd,    o));
	ratios->Add(Ratio(hhdsym, o));
      }
    }

    // If we have data from HHD's analysis, then do the ratio of 
    // our result to that data. 
    if (hasHhd) { 
      TH1F* t1 = static_cast<TH1F*>(dndeta->Clone(Form("%s_%s", 
						       dndeta->GetName(), 
						       hhd->GetName())));
      TH1F* t2 = static_cast<TH1F*>(sym->Clone(Form("%s_%s", 
						    sym->GetName(), 
						    hhdsym->GetName())));
      t1->SetTitle(Form("%s / %s", dndeta->GetTitle(), hhd->GetTitle()));
      t2->SetTitle(Form("%s / %s", sym->GetTitle(), hhdsym->GetTitle()));
      t1->Divide(hhd);
      t2->Divide(hhdsym);
      t1->SetMarkerColor(hhd->GetMarkerColor());
      t2->SetMarkerColor(hhdsym->GetMarkerColor());
      ratios->Add(t1);
      ratios->Add(t2);
    }

    // Do comparison to MC 
    if (mc) { 
      TH1D* t1 = static_cast<TH1D*>(dndeta->Clone(Form("%s_%s", 
						       dndeta->GetName(), 
						       mc->GetName())));
      TH1D* t2 = static_cast<TH1D*>(sym->Clone(Form("%s_%s", 
						    sym->GetName(), 
						    mcsym->GetName())));
      t1->SetTitle(Form("%s / %s", dndeta->GetTitle(), mc->GetTitle()));
      t2->SetTitle(Form("%s / %s", sym->GetTitle(), mcsym->GetTitle()));
      t1->Divide(mc);
      t2->Divide(mcsym);
      t1->SetMarkerColor(mc->GetMarkerColor());
      t2->SetMarkerColor(mcsym->GetMarkerColor());
      ratios->Add(t1);
      ratios->Add(t2);
    }

    // Check if we have ratios 
    Bool_t   hasRatios = (ratios->GetHists() && 
			  (ratios->GetHists()->GetEntries() > 0));
    Info("MakeRatios", "Got a total of %d ratios", !hasRatios ? 0 :
	 ratios->GetHists()->GetEntries());

    if (!hasRatios) { delete ratios; ratios = 0; }
    return ratios;
  }

  //__________________________________________________________________
  /** 
   * Fix the apperance of the axis in a stack 
   * 
   * @param stack  stack of histogram
   * @param s      Scaling factor 
   * @param ytitle Y axis title 
   * @param force  Whether to draw the stack first or not 
   * @param ynDiv  Divisions on Y axis 
   */
  void FixAxis(THStack* stack, Double_t s, const char* ytitle, 
	       Int_t ynDiv=10, Bool_t force=true) 
  {
    if (force) stack->Draw("nostack e1");

    TH1* h = stack->GetHistogram();
    if (!h) return;

    h->SetXTitle("#eta");
    h->SetYTitle(ytitle);
    TAxis* xa = h->GetXaxis();
    TAxis* ya = h->GetYaxis();
    if (xa) { 
      xa->SetTitle("#eta");
      // xa->SetTicks("+-");
      xa->SetTitleSize(s*xa->GetTitleSize());
      xa->SetLabelSize(s*xa->GetLabelSize());
      xa->SetTickLength(s*xa->GetTickLength());
    }
    if (ya) { 
      ya->SetTitle(ytitle);
      ya->SetDecimals();
      // ya->SetTicks("+-");
      ya->SetNdivisions(ynDiv);
      ya->SetTitleSize(s*ya->GetTitleSize());
      ya->SetLabelSize(s*ya->GetLabelSize());
    }      
  }
  //__________________________________________________________________
  /** 
   * Compute the ratio of @a h to @a g.  @a g is evaluated at the bin
   * centers of @a h
   * 
   * @param h  Numerator 
   * @param g  Divisor 
   * 
   * @return h/g 
   */
  TH1* Ratio(const TH1* h, const TGraph* g) const 
  {
    if (!h || !g) return 0;

    TH1* ret = static_cast<TH1*>(h->Clone("tmp"));
    ret->SetName(Form("%s_over_%s", h->GetName(), g->GetName()));
    ret->SetTitle(Form("%s / %s", h->GetTitle(), g->GetTitle()));
    ret->Reset();
    ret->SetMarkerStyle(h->GetMarkerStyle());
    ret->SetMarkerColor(g->GetMarkerColor());
    ret->SetMarkerSize(0.9*g->GetMarkerSize());
    Double_t xlow  = g->GetX()[0];
    Double_t xhigh = g->GetX()[g->GetN()-1];
    if (xlow > xhigh) { Double_t t = xhigh; xhigh = xlow; xlow = t; }

    for (Int_t i = 1; i <= h->GetNbinsX(); i++) { 
      Double_t c = h->GetBinContent(i);
      if (c <= 0) continue;

      Double_t x = h->GetBinCenter(i);
      if (x < xlow || x > xhigh) continue; 

      Double_t f = g->Eval(x);
      if (f <= 0) continue; 

      ret->SetBinContent(i, c / f);
      ret->SetBinError(i, h->GetBinError(i) / f);
    }
    if (ret->GetEntries() <= 0) { delete ret; ret = 0; }
    return ret;
  }

  //__________________________________________________________________
  /** 
   * Make an extension of @a h to make it symmetric about 0 
   * 
   * @param h Histogram to symmertrice 
   * 
   * @return Symmetric extension of @a h 
   */
  TH1* Symmetrice(const TH1* h) const
  {
    fOut->cd();

    Int_t nBins = h->GetNbinsX();
    TH1*  s     = new TH1D(Form("%s_mirror", h->GetName()),
			   Form("%s (mirrored)", h->GetTitle()), 
			   nBins, 
			   -h->GetXaxis()->GetXmax(), 
			   -h->GetXaxis()->GetXmin());
    s->SetMarkerColor(h->GetMarkerColor());
    s->SetMarkerSize(h->GetMarkerSize());
    s->SetMarkerStyle(h->GetMarkerStyle()+4);
    s->SetFillColor(h->GetFillColor());
    s->SetFillStyle(h->GetFillStyle());
    // s->SetDirectory(0);

    // Find the first and last bin with data 
    Int_t first = nBins+1;
    Int_t last  = 0;
    for (Int_t i = 1; i <= nBins; i++) { 
      if (h->GetBinContent(i) <= 0) continue; 
      first = TMath::Min(first, i);
      last  = TMath::Max(last,  i);
    }
    
    Double_t xfirst = h->GetBinCenter(first-1);
    Int_t    f1     = h->GetXaxis()->FindBin(-xfirst);
    Int_t    l2     = s->GetXaxis()->FindBin(xfirst);
    for (Int_t i = f1, j=l2; i <= last; i++,j--) { 
      s->SetBinContent(j, h->GetBinContent(i));
      s->SetBinError(j, h->GetBinError(i));
    }
    return s;
  }
  //__________________________________________________________________
  /** 
   * Rebin a histogram 
   * 
   * @param h     Histogram to rebin
   * @param rebin Rebinning factor 
   * 
   * @return 
   */
  virtual void Rebin(TH1* h, Int_t rebin) const
  { 
    if (rebin <= 1) return;

    Int_t nBins = h->GetNbinsX();
    if(nBins % rebin != 0) {
      Warning("Rebin", "Rebin factor %d is not a devisor of current number "
	      "of bins %d in the histogram %s", rebin, nBins, h->GetName());
      return;
    }
    
    // Make a copy 
    TH1* tmp = static_cast<TH1*>(h->Clone("tmp"));
    tmp->Rebin(rebin);
    tmp->SetDirectory(0);

    // The new number of bins 
    Int_t nBinsNew = nBins / rebin;
    for(Int_t i = 1;i<= nBinsNew; i++) {
      Double_t content = 0;
      Double_t sumw    = 0;
      Double_t wsum    = 0;
      Int_t    nbins   = 0;
      for(Int_t j = 1; j<=rebin;j++) {
	Int_t bin = (i-1)*rebin + j;
	if(h->GetBinContent(bin) <= 0) continue;
	Double_t c =  h->GetBinContent(bin);
	Double_t w = 1 / TMath::Power(c,2);
	content    += c;
	sumw       += w;
	wsum       += w * c;
	nbins++;
      }
      
      if(content > 0 ) {
	tmp->SetBinContent(i, wsum / sumw);
	tmp->SetBinError(i,TMath::Sqrt(sumw));
      }
    }

    // Finally, rebin the histogram, and set new content
    h->Rebin(rebin);
    for(Int_t i =1;i<=nBinsNew; i++) {
      h->SetBinContent(i,tmp->GetBinContent(i));
      // h->SetBinError(i,tmp->GetBinError(i));
    }
    
    delete tmp;
  }
};

//____________________________________________________________________
//
// EOF
//

