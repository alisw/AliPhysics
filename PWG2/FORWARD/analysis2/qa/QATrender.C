/**
 * @file   QATrender.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Thu Nov 17 12:28:17 2011
 * 
 * @brief  Class to make the QA trending tree
 * 
 * @ingroup pwg2_forward_qa_scripts
 */
#ifndef __CINT__
# include <TFile.h>
# include <TTree.h>
# include <TH1.h>
# include <TH2.h>
# include <TList.h>
# include <TStyle.h>
# include <TCanvas.h>
# include <TPad.h>
# include <TLatex.h>
# include <TMath.h>
# include <THStack.h>
# include <TLegend.h>
# include <TLegendEntry.h>
# include <TLine.h>
# include <TLinearFitter.h>
# include <TF1.h>
# include <TSystem.h>
# include <fstream>
# include <TFile.h>
# include "QAStructs.h"
# include "QARing.h"
# include "QABase.h"
#else
class TList;
class TFile;
class TH1;
class Quantity;
class RingQuantity;
class Merge;
class FitStatus;
class Correlation;
class QARing;
class QABase;
class Global;
class THStack;
class TLatex;
class TVirtualPad;
#endif 


// --- A quantity ----------------------------------------------------
/** 
 * Class to make the QA trending tree
 * 
 * @ingroup pwg2_forward_qa_scripts
 */
struct QATrender : public QABase
{
public:
  /******************************************************************/
  /** 
   * A ring object
   * 
   */
  struct Ring : public QARing
  { 
    /** 
     * Constructor 
     * 
     * @param d Detector 
     * @param r Ring
     */
    Ring(UShort_t d, Char_t r) : QARing(d, r)  {}
    /** 
     * Get the detector specific list 
     * 
     * @param parent Parent list  
     * 
     * @return Found list or null
     */
    TList* GetDetectorList(const TList* parent) const 
    {
      TObject* o = QATrender::GetSubList(parent, Form("FMD%d%c", fD, fR));
      return static_cast<TList*>(o);
    }
    /** 
     * Calculate numbers from read histograms
     * 
     * @param parent Parent list 
     * @param name   Name of histogram
     * @param q      Qauntity to store in
     * 
     * @return true on success 
     */
    Bool_t ExtractYQuantity(const TList* parent, const char* name, 
			    Quantity* q)
    {
      q->mean  = 0;
      q->var   = 0;
      q->min   = 0;
      q->max   = 0;

      TH1* h = QATrender::GetHistogram(parent, name);
      if (!h) {
	// Warning("ExtractYQuantity", "Histogram %s not found", name);
	return false;
      }
      
      Double_t sum   = 0;
      Double_t sumw  = 0;
      Double_t min   = 1e12;
      Double_t max   = 0;
      for (Int_t i = 1; i <= h->GetNbinsX(); i++) {
	Double_t y = h->GetBinContent(i);
	Double_t e = h->GetBinError(i);
	Double_t w = 1;
	if (y <= 1e-12) continue;
	if (e != 0) w = 1 / (e*e);
	sum   += w * y;
	sumw  += w;
	min  =  TMath::Min(min, y);
	max  =  TMath::Max(max, y);
	// Info("", "  %s %3d: y=%f,\t e=%f,\t w=%f", name, i, y, e, w); 
      }
      q->min  = min;
      q->max  = max;
      if (sumw <= 1e-6) {
	// Warning("ExtractYQuantity", 
	//         "Sum of weights for %s too small %g<1e-6",
	//         name, sumw);
	return true;
      }
      q->mean = sum / sumw;
      q->var  = TMath::Sqrt(1/sumw);
      // Info(name, " mean=%f, var=%f", q->mean, q->var);
      return true;
    }
    /** 
     * Process data from the energy loss fits
     * 
     * @param parent Parent list 
     * 
     * @return true on success 
     */
    Bool_t ProcessEnergyLoss(const TList* parent)
    {
      TList* det = GetDetectorList(parent);
      if (!det) return false;

      TList* res = QATrender::GetSubList(det, "FitResults");
      if (!res) return false;
      
      ExtractYQuantity(res, Form("FMD%d%c_chi2",  fD, fR), fChi2);
      ExtractYQuantity(res, Form("FMD%d%c_c",     fD, fR), fC);
      ExtractYQuantity(res, Form("FMD%d%c_delta", fD, fR), fDelta);
      ExtractYQuantity(res, Form("FMD%d%c_xi",    fD, fR), fXi);
      ExtractYQuantity(res, Form("FMD%d%c_sigma", fD, fR), fSigma);

      TH1* status = QATrender::GetHistogram(res, "status");
      if (!status) return false;
      fFitStatus->nLow        = UShort_t(status->GetBinContent(3));
      fFitStatus->nCandidates = UShort_t(status->GetBinContent(4));
      fFitStatus->nFitted     = UShort_t(status->GetBinContent(5));
      
      return true;
    }
    /** 
     * Process data on neighbors 
     * 
     * @param parent Parent list 
     * @param p      (optional) Pad to draw
     * 
     * @return true on success
     */
    Bool_t ProcessNeighbors(const TList* parent, TVirtualPad* p) 
    {
      if (!p) return true;  // Case of no drawing 

      TList* det = GetDetectorList(parent);
      if (!det) return false;

      TH1* before = QATrender::GetHistogram(det, "neighborsBefore");
      TH1* after  = QATrender::GetHistogram(det, "neighborsAfter");
      if (!before || !after) return false;
      
      if (before->GetMaximum() > 0) p->SetLogz();
      p->SetFillColor(0);
      if (fD == 3) p->SetRightMargin(0.15);

      before->SetTitle(Form("FMD%d%c",fD,fR));
      before->Draw("colz");
      after->Draw("same box");

      before->GetXaxis()->SetRangeUser(-.5, 2);
      before->GetYaxis()->SetRangeUser(-.5, 2);

      TLatex* ltx = new TLatex(p->GetLeftMargin()+.01, 
			       p->GetBottomMargin()+.01, 
			       before->GetTitle());
      ltx->SetNDC();
      ltx->SetTextSize(.07);
      ltx->Draw();      

      return true;
    }
    /** 
     * Process data on single, double, and triple hits
     * 
     * @param parent Parent list
     * @param p      (optional) Pad to draw
     * 
     * @return true on success
     */
    Bool_t Process123(const TList* parent, TVirtualPad* p) 
    {
      TList* det = GetDetectorList(parent);
      if (!det) return false;

      fMerge->one     = 0;
      fMerge->two     = 0;
      fMerge->three   = 0;

      TH1* one   = QATrender::GetHistogram(det, "singleEloss");
      TH1* two   = QATrender::GetHistogram(det, "doubleEloss");
      TH1* three = QATrender::GetHistogram(det, "tripleEloss");
      if (!one || !two || !three) return false;

      Int_t    nOne   = one->GetEntries();
      Int_t    nTwo   = two->GetEntries();
      Int_t    nThree = three->GetEntries();
      Int_t    total  = nOne + nTwo + nThree; 
      if (total > 0) { 
	fMerge->one     = Double_t(nOne) / total;
	fMerge->two     = Double_t(nTwo) / total;
	fMerge->three   = Double_t(nThree) / total;
      }

      if (!p) return true;

      one->SetStats(0);
      one->SetTitle(Form("FMD%d%c", fD, fR));
      one->GetXaxis()->SetRangeUser(0, 8);
      
      if (one->GetMaximum() > 0) p->SetLogy();
      p->SetFillColor(0);
      if ((p->GetNumber() % 3) != 1) p->SetLeftMargin(.1);
      
      one->SetLineStyle(1);
      two->SetLineStyle(2);
      three->SetLineStyle(3);

      one->Draw();
      two->Draw("same");
      three->Draw("same");

      // Double_t ts  = 0.03;
      Double_t y = 1-p->GetTopMargin()-1.2*gStyle->GetTitleH(); // +ts;
      Double_t x = p->GetLeftMargin() + .25;
      TLatex* ltx = new TLatex(x, y, "Fraction of ...");
      ltx->SetNDC();
      ltx->SetTextColor(kBlue+3);
      ltx->SetTextAlign(13);
      // ltx->SetTextAlign(33);
      ltx->SetTextSize(.07);
      ltx->Draw();      

      // y -= .1;
      DrawText(ltx, x, y, "Singles:", Form("%5.2f%%", fMerge->one * 100));
      DrawText(ltx, x, y, "Doubles:", Form("%5.2f%%", fMerge->two * 100));
      DrawText(ltx, x, y, "Triples:", Form("%5.2f%%", fMerge->three * 100));

      p->cd();

      return true;
    }
    /** 
     * Process energy loss distributions 
     * 
     * @param p1     First parent list (sharing filter)
     * @param p2     Second parent list (density calculator)
     * @param p      (optional) Pad to draw
     * 
     * @return true on success
     */
    Bool_t ProcessELoss(TList* p1, TList* p2, TVirtualPad* p)
    {
      TList* det1 = GetDetectorList(p1);
      if (!det1) return false;

      TList* det2 = GetDetectorList(p2);
      if (!det2) return false;

      TH1* before    =  QATrender::GetHistogram(det1, "esdEloss");
      TH1* after     =  QATrender::GetHistogram(det1, "anaEloss");
      TH1* presented =  QATrender::GetHistogram(det2, "eloss");
      TH1* used      =  QATrender::GetHistogram(det2, "elossUsed");

      if (!before || !after || !presented || !used) return false;

      const Double_t lowCut = 0.15;

      Int_t low = before->GetXaxis()->FindBin(lowCut);
      Int_t ib  = Int_t(before->Integral(low,before->GetNbinsX()));
      Int_t ia  = Int_t(after->Integral(low,after->GetNbinsX()));
      Int_t ip  = Int_t(presented->Integral(low,presented->GetNbinsX()));
      Int_t iu  = Int_t(used->Integral(low,used->GetNbinsX()));

      Double_t dba = ib > 0 ? (100.*(ia-ib))/ib : 0;
      Double_t dbp = ib > 0 ? (100.*(ip-ib))/ib : 0;
      Double_t dbu = ib > 0 ? (100.*(iu-ib))/ib : 0;
      Double_t dap = ia > 0 ? (100.*(ip-ia))/ia : 0;
      Double_t dau = ia > 0 ? (100.*(iu-ia))/ia : 0;
      Double_t dpu = ip > 0 ? (100.*(iu-ip))/ip : 0;
      
      fDataLoss->merge   = dba;
      fDataLoss->density = dau;
      fDataLoss->full    = dbu;

      if (!p) return true;

      if (before->GetMaximum() > 0) p->SetLogy();
      p->SetFillColor(0);
      before->SetTitle(Form("FMD%d%c",fD,fR));

      before->SetLineStyle(1);
      after->SetLineStyle(2);
      used->SetLineStyle(3);

      before->Draw("");
      after->Draw("same");
      presented->Draw("same");
      used->Draw("same");
  
      Double_t ts  = 0.03;
      Double_t x   = p->GetLeftMargin() + .25;
      Double_t y   = 1-p->GetTopMargin()-gStyle->GetTitleH()+ts;
      TLatex*  ltx = new TLatex(x, y, Form("FMD%d%c", fD, fR));
      ltx->SetNDC();
      ltx->SetTextAlign(13);
      ltx->SetTextSize(ts);
      ltx->SetTextColor(kBlue+3);
      // ltx->Draw();
      // ltx->SetTextSize(.05);
      TString inte(Form("Integral [%4.2f,#infty]", lowCut));
      DrawText(ltx, x, y, Form("%s before:", inte.Data()), Form("%9d", ib));
      DrawText(ltx, x, y, Form("%s after:",  inte.Data()), Form("%9d", ia));
      DrawText(ltx, x, y, Form("%s input:",  inte.Data()), Form("%9d", ip));
      DrawText(ltx, x, y, Form("%s user:",   inte.Data()), Form("%9d", iu));
      TLine* l = new TLine;
      l->SetLineWidth(1);
      l->DrawLineNDC(x, y-0.9*ts, 1-p->GetRightMargin()-0.01, y-0.9*ts);
      if (ib != 0 && ia != 0) {
	DrawText(ltx, x, y, "Change (merging)", Form("%5.1f%%", dba));
	DrawText(ltx, x, y, "Change (input)",   Form("%5.1f%% (%5.1f%%)", 
						     dap, dbp));
	DrawText(ltx, x, y, "Change (use)",     Form("%5.1f%% (%5.1f%%)", 
						     dpu, dbu));
      }
      before->GetXaxis()->SetRangeUser(0, 4);
      p->cd();

      return true;
    }
    /** 
     * Process occupancy information 
     * 
     * @param parent Parent list 
     * @param p      (optional) Pad to draw
     * 
     * @return true on success
     */
    Bool_t ProcessOccupancy(const TList* parent, TVirtualPad* p) 
    {
      TList* det = GetDetectorList(parent);
      if (!det) return false;

      TH1* occ   = QATrender::GetHistogram(det, "occupancy");
      if (!occ) return false;

      fOccupancy->mean = occ->GetMean();
      fOccupancy->var  = occ->GetRMS();
      fOccupancy->min  = 1;
      fOccupancy->max  = 0;

      for (Int_t i = occ->GetNbinsX(); i >= 1; i--) {
	Double_t y = occ->GetBinContent(i);
	if (y < 1e-6) continue;
	fOccupancy->max = occ->GetXaxis()->GetBinUpEdge(i);
	break;
      }
      
      if (!p) return true;
      p->SetGridy();
      p->SetGridx();
      if (occ->GetMaximum() > 0) p->SetLogy();
      p->SetFillColor(0);
      p->SetRightMargin(0.01);
      occ->Rebin(4);
      occ->SetTitle(Form("FMD%d%c", fD, fR));
      occ->Draw("hist");

      // Double_t ts  = 0.03;
      Double_t y = 1-p->GetTopMargin()-1.2*gStyle->GetTitleH(); // +ts;
      Double_t x = p->GetLeftMargin() + .25;
      TLatex* ltx = new TLatex(x, y, "");
      ltx->SetNDC();
      ltx->SetTextColor(kBlue+3);
      ltx->SetTextAlign(13);
      // ltx->SetTextAlign(33);
      ltx->SetTextSize(.07);
      ltx->Draw();      

      // y -= .1;
      DrawText(ltx, x, y, "Mean:",     Form("%5.2f%%", fOccupancy->mean));
      DrawText(ltx, x, y, "Variance:", Form("%5.2f%%", fOccupancy->var));
      DrawText(ltx, x, y, "Max:",      Form("%5.2f%%", fOccupancy->max));

      p->cd();

      return true;
    }
    /** 
     * Process method correlations 
     * 
     * @param parent Parent list
     * @param p      (optional) Pad to draw
     * 
     * @return true on success
     */
    Bool_t ProcessCorrelation(const TList* parent, TVirtualPad* p) 
    {
      TList* det = GetDetectorList(parent);
      if (!det) return false;

      TH1* co   = QATrender::GetHistogram(det, "elossVsPoisson");
      if (!co) return false;
      TH2* corr = static_cast<TH2*>(co);

      fCorrelation->alpha = 0;
      fCorrelation->beta  = 0;
      fCorrelation->a     = 0;
      fCorrelation->ea    = 0;
      fCorrelation->b     = 0;
      fCorrelation->eb    = 0;
      fCorrelation->chi2  = 0;

      Double_t xmin = -1;
      Double_t xmax = corr->GetXaxis()->GetXmax();

      if (corr->GetEntries() > 0) {
	// Calculate the linear regression 
	// Double_t dx    = (xmax-xmin);
	// Double_t rxy   = corr->GetCorrelationFactor();
	Double_t sx    = corr->GetRMS(1);
	Double_t sy    = corr->GetRMS(2);
	Double_t sx2   = sx*sx;
	Double_t sy2   = sy*sy;
	Double_t cxy   = corr->GetCovariance();
	Double_t mx    = corr->GetMean(1);
	Double_t my    = corr->GetMean(2);
	Double_t delta = 1;
	if (TMath::Abs(cxy) > 1e-6) {
	  fCorrelation->beta  = ((sy2 - delta*sx2 + 
				  TMath::Sqrt(TMath::Power(sy2-delta*sx2,2) + 
					      4*delta*cxy*cxy)) / 2 / cxy);
	  fCorrelation->alpha = my - fCorrelation->beta * mx;
	}


#if 0
	TLinearFitter* fitter = new TLinearFitter(1);
	fitter->SetFormula("1 ++ x");
	for (Int_t i = 1; i <= corr->GetNbinsX(); i++) { 
	  Double_t x = corr->GetXaxis()->GetBinCenter(i);
	  if (x < -1 || x > xmax) continue;
	  for (Int_t j = 1; j <= corr->GetNbinsY(); j++) {
	    Double_t y = corr->GetYaxis()->GetBinCenter(j);
	    if (y < -1 || y > xmax) continue;
	    Double_t c = corr->GetBinContent(i,j);
	    if (c < .1) continue;
	    fitter->AddPoint(&x, y, c);
	  }
	}
	fitter->Eval();
	fCorrelation->a    = fitter->GetParameter(0);
	fCorrelation->ea   = fitter->GetParError(0); 
	fCorrelation->b    = fitter->GetParameter(1);
	fCorrelation->eb   = fitter->GetParError(1);
	Double_t chi2 = fitter->GetChisquare();
	Int_t    ndf  = (fitter->GetNpoints() - 
			 fitter->GetNumberFreeParameters() );
	fCorrelation->chi2 = chi2 / ndf;
#endif
      }
      if (!p) return true;

      p->SetGridy();
      p->SetGridx();
      if (corr->GetMaximum() > 0) p->SetLogz();
      // if (fD == 3) p->SetRightMargin(0.15);
      p->SetFillColor(0);
      corr->GetXaxis()->SetRangeUser(-1,xmax);
      corr->GetYaxis()->SetRangeUser(-1,xmax);
      corr->SetTitle(Form("FMD%d%c",fD,fR));
      corr->Draw("colz");

      TF1* f = new TF1("f", "[0]+[1]*x", xmin, xmax);
      f->SetParameters(fCorrelation->alpha, fCorrelation->beta);
      f->SetLineWidth(1);
      f->SetLineStyle(1);
      f->Draw("same");
      
      TLine* l = new TLine(-1,-1,xmax,xmax);
      l->SetLineWidth(1);
      l->SetLineStyle(2);
      l->SetLineColor(kBlack);
      l->Draw();

      Double_t x   = p->GetLeftMargin() + .05;
      Double_t y   = 1-p->GetTopMargin()-gStyle->GetTitleH(); // +ts;
      TLatex* ltx = new TLatex(x, y, "Deming regression: y=#alpha+#beta x");
      ltx->SetNDC();
      ltx->SetTextAlign(13);
      ltx->SetTextSize(0.06);
      ltx->SetTextColor(kBlue+3);
      ltx->Draw();

      DrawText(ltx, x, y, "#alpha:", Form("%5.3f", fCorrelation->alpha));
      DrawText(ltx, x, y, "#beta:",  Form("%5.3f", fCorrelation->beta));

#if 0
      DrawText(ltx, x, y, "A:",      Form("%5.3f#pm%5.3f", 
					  fCorrelation->a,
					  fCorrelation->ea));
      DrawText(ltx, x, y, "B:",      Form("%5.3f#pm%5.3f", 
					  fCorrelation->b,
					  fCorrelation->eb));
      DrawText(ltx, x, y, "#chi^{2}/#nu:", Form("%5.3f", fCorrelation->chi2));
#endif
      p->cd();

      return true;
    }
  };

  /******************************************************************/
  /** 
   * CTOR
   */
  QATrender(Bool_t keep=false, Bool_t single=false) 
    : QABase(single), 
      fCurrentFile(0),
      fSharingFilter(0),
      fEventInspector(0),
      fDensityCalculator(0),
      fEnergyFitter(0),
      fFiles(0), 
      fKeep(keep)
  {
    fFMD1i = new Ring(1, 'I'); 
    fFMD2i = new Ring(2, 'I'); 
    fFMD2o = new Ring(2, 'O'); 
    fFMD3i = new Ring(3, 'I'); 
    fFMD3o = new Ring(3, 'O'); 
  }
  /** 
   * DTOR
   */
  virtual ~QATrender() {} 
  /**
   * Copy CTOR
   * 
   */
  QATrender(const QATrender& o) 
    : QABase(o), 
      fCurrentFile(o.fCurrentFile),
      fSharingFilter(o.fSharingFilter),
      fEventInspector(o.fEventInspector),
      fDensityCalculator(o.fDensityCalculator),
      fEnergyFitter(o.fEnergyFitter),
      fFiles(0), 
      fKeep(o.fKeep)
  {} 
  /**
   * Assignment operator 
   *
   * @return Reference to this 
   */
  QATrender operator=(const QATrender&) { return *this; }


  // --- Interface ---------------------------------------------------
  /** 
   * Add a file to be processed
   * 
   * @param filename Name of file 
   */
  void AddFile(const char* filename)
  {
    fFiles.Add(new TObjString(filename));
  }
  /** 
   * Run the job
   * 
   */
  void Run()
  {
    Init(false);
    TIter next(&fFiles);
    TObject* o = 0;
    while ((o = next())) {
      ProcessOne(o->GetName());
    }
    Finish();
  }
  /** 
   * Finish the job
   * 
   */
  void Finish()
  {
    if (!fOutput) return;
    fOutput->Write();
    fOutput->Close();
    fOutput = 0;
    fTree   = 0;
    gSystem->Exec(Form("chmod g+rw %s", OutputName()));
  }
  /** 
   * Process a single file 
   * 
   * @param filename File to open. 
   * 
   * @return true on success 
   */
  Bool_t ProcessOne(const char* filename)
  {
    if (fCurrentFile) { 
      fCurrentFile->Close();
      fCurrentFile = 0;
    }

    fCurrentFile = TFile::Open(filename, "READ");
    if (!fCurrentFile) { 
      Error("ProcessOne", "Failed to open %s", filename);
      return false;
    }
    
    if (!GetLists()) { 
      Error("ProcessOne", "Failed to get lists from %s", filename);
      return false;
    }
    
    if (!ProcessGlobal()) { 
      Error("ProcessOne", "Failed to get global stuff from %s", filename);
      return false;
    }
    fTeXName = Form("qa_%09d", fGlobal->runNo);
    MakeCanvas(Form("QA plots for run %d", fGlobal->runNo));
    Bool_t eloss = ProcessELossFitter();
    Bool_t merge = ProcessSharingFilter();
    Bool_t dense = ProcessDensityCalculator();
    if (fTree) fTree->Fill();
    CloseCurrent();

    return eloss && merge && dense;
  }
  // --- Processing member functions ---------------------------------
  /** 
   * Get global stuff 
   * 
   * @return true on success
   */
  Bool_t ProcessGlobal()
  {
    TObject* oRun = GetObject(fEventInspector, "runNo");
    if (!oRun) return false;
    
    fGlobal->runNo = oRun->GetUniqueID();
    TH1* oAcc = GetHistogram(fEventInspector,"nEventsAccepted");
    if (!oAcc) return false; 

    fGlobal->nAccepted = oAcc->GetEntries();
    fGlobal->meanVz    = oAcc->GetMean();
    fGlobal->sigmaVz   = oAcc->GetRMS();
    
    return true;
  }
  /** 
   * Clean a stack for histograms that match @a what 
   * 
   * @param stack Stack to clean 
   * @param what  Pattern 
   */
  void CleanStack(THStack* stack, const char* what)
  {
    TList*   l = stack->GetHists();

    // Clean up list of histogram.  Histograms with no entries or 
    // no functions are deleted.  We have to do this using the TObjLink 
    // objects stored in the list since ROOT cannot guaranty the validity 
    // of iterators when removing from a list - tsck.  Should just implement
    // TIter::Remove(). 
    TObjLink* lnk = l->FirstLink();
    while (lnk) {
      TObject* o = lnk->GetObject();
      TString  s(o->GetName());
      if (s.Contains(what)) {
	TObjLink* keep = lnk->Next();
	l->Remove(lnk);
	lnk = keep;
	continue;
      }
      lnk = lnk->Next();
    }
  }
  /** 
   * Process the information from the energy loss fitter 
   * 
   * 
   * @return true on succes
   */
  Bool_t ProcessELossFitter()
  {
    MakeCanvasTitle("Summary of energy loss fits");

    THStack* chi2  = static_cast<THStack*>(GetObject(fEnergyFitter, "chi2"));
    THStack* c     = static_cast<THStack*>(GetObject(fEnergyFitter, "c"));
    THStack* delta = static_cast<THStack*>(GetObject(fEnergyFitter, "delta"));
    THStack* xi    = static_cast<THStack*>(GetObject(fEnergyFitter, "xi"));
    THStack* sigma = static_cast<THStack*>(GetObject(fEnergyFitter, "sigma"));

    if (!chi2)  return false;
    if (!c)     return false;
    if (!delta) return false;
    if (!xi)    return false;
    if (!sigma) return false;

    CleanStack(chi2, "_n");
    CleanStack(c,    "status");
    if (fCanvas) {
      const int nL = 3;
      THStack* stacks[] = { chi2, c, delta, xi, sigma, 0 };
      for (int i = 0; i < 5; i++) { 
	THStack*     stack = stacks[i];
	TVirtualPad* p     = GetPad(i+1);
	// stack->GetHists()->ls();

	p->SetLeftMargin(.6/nL);
	p->SetTopMargin(.01);
	p->SetRightMargin(.01);
	p->SetFillColor(0);
	p->SetFillStyle(0);
	p->SetGridx();
	stack->Draw("nostack");
	stack->GetHistogram()->SetYTitle(stack->GetTitle());
	stack->GetHistogram()->SetXTitle("#eta");
	
	TAxis* yaxis = stack->GetHistogram()->GetYaxis();
	if (i == 0) yaxis->SetRangeUser(0,20); // Chi2
	if (i == 1) stack->SetMaximum(1);      // c
	if (i == 2) stack->SetMaximum(1);      // delta
	if (i == 3) stack->SetMaximum(0.1);   // xi
	if (i == 4) stack->SetMaximum(0.5);    // sigma{,n}
	// if (i == 0) p->SetLogy();
	yaxis->SetTitleSize(0.3/nL);
	yaxis->SetLabelSize(0.08);
	yaxis->SetTitleOffset(3/nL);
	yaxis->SetNdivisions(5);
	yaxis->SetTitleFont(42);
	yaxis->SetLabelFont(42);
	yaxis->SetDecimals();
	
	TAxis* xaxis = stack->GetHistogram()->GetXaxis();
	xaxis->SetTitleSize(0.3/nL);
	xaxis->SetLabelSize(0.08);
	xaxis->SetTitleOffset(2./nL);
	xaxis->SetNdivisions(10);
	xaxis->SetTitleFont(42);
	xaxis->SetLabelFont(42);
	xaxis->SetDecimals();      
	
	stack->Draw("nostack");
	p->cd();
      }
      TVirtualPad* p     = GetPad(6);
      p->SetFillColor(kWhite);
      Double_t x = .3;
      Double_t y = .8;
      TLatex* l = new TLatex(x, y, "Fits to #Delta (energy loss) spectra");
      l->SetTextColor(kBlue+3);
      l->SetNDC();
      l->Draw();
      x = .05;
      y -= 2 * 1.2*l->GetTextSize();
      l->DrawLatex(x, y, "F(#Delta;c,#Delta_{p},#xi,#sigma)="
		   "#frac{c}{#sqrt{2#pi}#sigma}#int_{-#infty}^{#infty}d#Delta'"
		   "L(#Delta;#Delta',#xi) G(#Delta_{p};#Delta',#sigma^{2})");
      y -= 1.2*l->GetTextSize();
      x += .1;
      DrawText(l, x, y, "#chi^{2}/#nu", "Goodness of fit",         .2);
      DrawText(l, x, y, "c",            "Overall constant",        .2); 
      DrawText(l, x, y, "#Delta_{p}",   "Most probable value",     .2);
      DrawText(l, x, y, "#xi",          "'Width' of Landau (L)",   .2);
      DrawText(l, x, y, "#sigma",       "'Width' of Gaussian (G)", .2);
      
      // stack->GetHists()->ls();

      PrintCanvas("fitResults", fGlobal->runNo);
    }
    
    static_cast<Ring*>(fFMD1i)->ProcessEnergyLoss(fEnergyFitter);
    static_cast<Ring*>(fFMD2i)->ProcessEnergyLoss(fEnergyFitter);
    static_cast<Ring*>(fFMD2o)->ProcessEnergyLoss(fEnergyFitter);
    static_cast<Ring*>(fFMD3i)->ProcessEnergyLoss(fEnergyFitter);
    static_cast<Ring*>(fFMD3o)->ProcessEnergyLoss(fEnergyFitter);

    return true;
  }
  /** 
   * Process the information from the sharing filter 
   * 
   * 
   * @return true on success
   */
  Bool_t ProcessSharingFilter()
  {
    // --- Neighbors -------------------------------------------------
    MakeCanvasTitle("Correlation of neighboring strips");
    
    for (Int_t i = 1; i <= 6; i++) { 
      TVirtualPad* p = GetPad(i);

      Ring* r = GetRing(i);
      if (!r) continue;

      r->ProcessNeighbors(fSharingFilter, p);
    }
    TVirtualPad* p = 0;
    if ((p = GetPad(4))) {
      p->SetFillColor(kWhite);
      
      TLatex* l = new TLatex(.2, .7, "Gradient: before merging");
      l->SetNDC();
      l->SetTextColor(kBlue+3);
      l->Draw();
      l->DrawText(.2, .6, "Boxes: after merging");
      
      fCanvas->cd();
      PrintCanvas("neighbors", fGlobal->runNo);
    }

    // --- 123 -------------------------------------------------------
    MakeCanvasTitle("#Delta for singles, doubles, and triples");
    
    for (Int_t i = 1; i <= 6; i++) { 
      p = GetPad(i);

      Ring* r = GetRing(i);
      if (!r) continue;

      r->Process123(fSharingFilter, p);
    }
    if ((p = GetPad(4))) { 
      TLegend* ll = new TLegend(.2, .2, .8, .8);
      ll->SetFillColor(0);
      ll->SetBorderSize(0);
      TLegendEntry* e = ll->AddEntry("dummy", "Singles", "l");
      e->SetLineStyle(1);
      e = ll->AddEntry("dummy", "Doubles", "l");
      e->SetLineStyle(2);
      e = ll->AddEntry("dummy", "Triples", "l");
      e->SetLineStyle(3);
      ll->Draw();
      
      PrintCanvas("123", fGlobal->runNo);
    }
    return true;
  }
  /** 
   * Process the information from the density calculator 
   * 
   * 
   * @return true on success
   */
  Bool_t ProcessDensityCalculator()
  {
    // --- ELoss -----------------------------------------------------
    MakeCanvasTitle("Energy loss from ESD, after merging, used");
    
    for (Int_t i = 1; i <= 6; i++) { 
      TVirtualPad* p = GetPad(i);

      Ring* r = GetRing(i);
      if (!r) continue;

      r->ProcessELoss(fSharingFilter, fDensityCalculator, p);
    }
    TVirtualPad* p = 0;
    if ((p = GetPad(4))) {
      TLegend* ll = new TLegend(.2, .2, .8, .8);
      ll->SetFillColor(0);
      ll->SetBorderSize(0);
      TLegendEntry* e = ll->AddEntry("dummy", "From ESDs", "l");
      e->SetLineStyle(1);
      e = ll->AddEntry("dummy", "After Merging", "l");
      e->SetLineStyle(2);
      e = ll->AddEntry("dummy", "Used", "l");
      e->SetLineStyle(3);
      ll->Draw();
      
      PrintCanvas("recAna", fGlobal->runNo);
    }

    // --- Occupancy -------------------------------------------------
    MakeCanvasTitle("Occupancy");
    
    for (Int_t i = 1; i <= 6; i++) { 
      p = GetPad(i);

      Ring* r = GetRing(i);
      if (!r) continue;

      r->ProcessOccupancy(fDensityCalculator, p);
    }
    if ((p = GetPad(4))) {
      TLatex* ltx = new TLatex(.2, .8, "Calculated assuming Poisson stat.");
      ltx->SetNDC();
      ltx->SetTextColor(kBlue+3);
      ltx->Draw();
      
      TObject* etaL = GetObject(fDensityCalculator, "etaLumping");
      TObject* phiL = GetObject(fDensityCalculator, "phiLumping");
      if (etaL && phiL) 
	ltx->DrawLatex(.2, .7, Form("Regions of %s strips #times %s sectors",
				    etaL->GetTitle(), phiL->GetTitle()));
      
      PrintCanvas("occupancy", fGlobal->runNo);
    }

    // --- Correlation of methods ------------------------------------
    MakeCanvasTitle("Correlation of N_{ch} methods");
    
    for (Int_t i = 1; i <= 6; i++) { 
      p = GetPad(i);

      Ring* r = GetRing(i);
      if (!r) continue;

      r->ProcessCorrelation(fDensityCalculator, p);
    }
    if ((p = GetPad(4))) {
      TLatex* ltx = new TLatex(.2, .8, "Correlation of N_{ch} methods");
      ltx->SetNDC();
      ltx->SetTextColor(kBlue+3);
      ltx->Draw();
      ltx->DrawLatex(.24, .7, "From #DeltaE fits along x");
      ltx->DrawLatex(.24, .6, "From Poisson assumption along y");
      ltx->DrawLatex(.24, .4, "Solid line: regression");
      ltx->DrawLatex(.24, .3, "Dashed line: x=y to guide the eye");
      
      PrintCanvas("elossVsPoisson", fGlobal->runNo);
    }

    return true;
  }

  // --- Utilities ---------------------------------------------------
  /** 
   * Close current file if open, and reset other pointers
   * 
   */
  void CloseCurrent()
  {
    if (fCurrentFile) {
      fCurrentFile->Close();
      fCurrentFile       = 0;
      fSharingFilter     = 0;
      fEventInspector    = 0;
      fDensityCalculator = 0;
      fEnergyFitter      = 0;
    }
    bool keep = fKeep; 
    if (fSingle) keep = true;
    Close(!keep);
  }
  /** 
   * Get the ring corresponding to a pad 
   * 
   * @param padNo Pad number 
   * 
   * @return Pointer to ring
   */
  Ring* GetRing(Int_t padNo) 
  {
    QARing* r = 0;
    switch(padNo) { 
    case 1: r = fFMD1i; break;
    case 2: r = fFMD2i; break;
    case 3: r = fFMD3i; break;
    case 5: r = fFMD2o; break;
    case 6: r = fFMD3o; break;
    }
    if (r) return static_cast<Ring*>(r);
    return 0;
  }
  /** 
   * Get the pad corresponding to a number.  Also clears the pad 
   * 
   * @param padNo Pad number to get
   * 
   * @return Pointer to pad or null
   */
  TVirtualPad* GetPad(Int_t padNo)
  {
    if (!fCanvas) return 0;
    TVirtualPad* p = fCanvas->cd(padNo);
    if (!p) return 0;
    p->Clear();
    p->SetFillColor(kWhite);
    return p;
  }

  /** 
   * Make canvas title. Canvas is divided into 3x2 pads 
   * 
   * @param what Title on canvas
   */
  void MakeCanvasTitle(const char* what)
  {
    if (!fCanvas) return;
    fCanvas->cd();

    CanvasTitle(what);
    fCanvas->Divide(3,2,0,0);
  }
  // --- List utilities ----------------------------------------------
  /** 
   * Get a sub-list from parent list 
   * 
   * @param parent Parent list 
   * @param name   Name of sub-list 
   * 
   * @return Pointer to the list 
   */
  static TList* GetSubList(const TList* parent, const char* name)
  {
    TList* tmp = static_cast<TList*>(parent->FindObject(name));
    if (!tmp) { 
      Error("GetLists", "List %s not found in %s", name, 
	    parent->GetName());
      return 0;
    }
    return tmp;
  }
  /** 
   * Get the lists from the file 
   * 
   * 
   * @return true on success
   */    
  Bool_t GetLists()
  {
    if (!fCurrentFile) return 0;
    
    const char* folder = "ForwardResults";
    TList* forward = static_cast<TList*>(fCurrentFile->Get(folder));
    if (!forward) { 
      Error("GetLists", "List %s not found in %s", folder, 
	    fCurrentFile->GetName());
      return false;
    }

    fSharingFilter     = GetSubList(forward, "fmdSharingFilter");
    fEventInspector    = GetSubList(forward, "fmdEventInspector");
    fDensityCalculator = GetSubList(forward, "fmdDensityCalculator");
    fEnergyFitter      = GetSubList(forward, "fmdEnergyFitter");
    
    if (!fSharingFilter)     return false; 
    if (!fEventInspector)    return false;
    if (!fDensityCalculator) return false;
    if (!fEnergyFitter)     return false;

    return true;
  }
  /** 
   * Find and object in a list
   * 
   * @param list List to look in
   * @param name Name of object 
   * 
   * @return Pointer to object or null
   */
  static TObject* GetObject(const TList* list, const char* name)
  {
    if (!list) return 0;
    TObject* o = list->FindObject(name);
    if (!o) { 
      Error("GetObject", "Failed to find object %s in %s", 
	    name, list->GetName());
      return 0;
    }
    return o;
  }
  /** 
   * Get a histogram from a list 
   * 
   * @param list List 
   * @param name Name of histogram
   * 
   * @return Pointer to object or null
   */
  static TH1* GetHistogram(const TList* list, const char* name)
  {
    return static_cast<TH1*>(GetObject(list, name));
  }
  /** 
   * Draw some text on a pad 
   * 
   * @param l   LaTeX object to use 
   * @param x   x coordinate 
   * @param y   y coordinate (incremented on return)
   * @param c1  First text 
   * @param c2  Second text 
   * @param dx  Distance between c1 start and c2 start 
   */
  static void 
  DrawText(TLatex* l, Double_t x, Double_t& y, const char* c1, const char* c2,
	   Double_t dx=.4)
  {
    y -= 1.2*l->GetTextSize();
    l->DrawLatex(x,    y, c1);
    l->DrawLatex(x+dx, y, c2);
  }
  

  // --- Members -----------------------------------------------------
  TFile*  fCurrentFile;          // Current input file 
  TList*  fSharingFilter;        // Sharing filter list
  TList*  fEventInspector;       // Event inspector list
  TList*  fDensityCalculator;    // Density calculator list 
  TList*  fEnergyFitter;         // Energy fitter list 
  TList   fFiles;                // List of files to process 
  Bool_t  fKeep;                 // Keep PNGs
};
//
// EOF
//



