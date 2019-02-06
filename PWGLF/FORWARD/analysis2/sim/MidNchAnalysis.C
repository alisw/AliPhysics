#ifndef MIDNCHANALYSIS_C
#define MIDNCHANALYSIS_C
#ifndef __CINT__
# include <TH1.h>
# include <TH2.h>
# include <TMath.h>
# include <TParticle.h>
# include <TObjArray.h>
# include <TString.h>
# include <TGraphErrors.h>
# include <TProfile.h>
#else
class TH1;
class TH2;
class TParticle;
class TMap;
class TProfile;
#endif
#include "FastAnalysis.C"
#include "FastCentHelper.C"

namespace MidNch {
  /** 
   * Base class 
   */
  struct Base : public FastAnalysis
  {
    FastCentHelper fHelper;
    TProfile* fNpartVsNch; //!
    TProfile* fNpartB;     //!
    TProfile* fMeanNpart;  //!
    TProfile* fMeanNch;    //!
    
    Int_t     fNch;
    Int_t     fNpart;
    Double_t  fCent;
    /** 
     * Constructor
     * 
     * @param verbose 
     * @param monitor 
     */
    Base(Bool_t verbose=false, Int_t monitor=0)
      : FastAnalysis(verbose, monitor),
	fHelper("B"),
	fNpartVsNch(0),
	fNpartB(0),
	fMeanNpart(0),
	fMeanNch(0),
	fNch(0),
	fNpart(0),
	fCent(-1)
    {
      fCentMethod = fHelper.fCentMeth;
    }
    /** 
     * Set-up centrality estimator
     *
     * @return true on success 
     */
    virtual Bool_t SetupEstimator()
    {
      if (!FastAnalysis::SetupEstimator()) return false;
    
      fHelper.CreateDiagnostics(fOutput, fCentHist);
      return true;
    }    

    /** 
     * Create output object
     * 
     * @return Output object 
     */
    virtual TProfile* CreateOutput() = 0;
    /** 
     * Called on each slave before execution 
     */
    virtual void SlaveBegin(TTree*)
    {
      Info("SlaveBegin", "Creating histogram");
      fNpartVsNch = CreateOutput();
      fNpartVsNch->SetDirectory(0);
      fNpartVsNch->SetName("nPartVsnCh");
      fNpartVsNch->SetXTitle("\\langle N_{\\mathrm{part}}\\rangle"); 
      fNpartVsNch->SetYTitle("\\left\\langle\\frac{\\mathrm{d}N_{\\mathrm{ch}}}"
			     "{\\mathrm{d}\\eta}\\right\\rangle");
      fNpartVsNch->SetMarkerColor(kRed+2);
      fNpartVsNch->SetMarkerStyle(20);
      fNpartVsNch->SetFillStyle(1001);
      fNpartVsNch->SetFillColor(kYellow+1);

      fNpartB = new TProfile("nPartB", "N_{\\mathrm{part}}\\hbox{ vs. }b",
			     80, 0, 20);
      fNpartB->SetXTitle("b\\hbox{ (fm)}");
      fNpartB->SetYTitle("\\langle N_{\\mathrm{part}}\\rangle");
      fNpartB->SetDirectory(0);
      fNpartB->SetMarkerColor(kBlue+2);
      fNpartB->SetMarkerStyle(21);

      fMeanNpart = 0;
      fMeanNch   = 0;
      const TAxis& centAxis = (*fHelper.fCentAxis);
      if (centAxis.GetXbins() && centAxis.GetXbins()->GetArray()) {
	fMeanNpart = new TProfile("meanNpart", "", centAxis.GetNbins(),
				  centAxis.GetXbins()->GetArray());
	fMeanNch   = new TProfile("meanNch", "", centAxis.GetNbins(),
				  centAxis.GetXbins()->GetArray());
      }
      else {
	fMeanNpart = new TProfile("meanNpart","", centAxis.GetNbins(),
				  centAxis.GetXmin(), centAxis.GetXmax());
	fMeanNch   = new TProfile("meanNch","", centAxis.GetNbins(),
				  centAxis.GetXmin(), centAxis.GetXmax());
      }
      fMeanNpart->SetDirectory(0);
      fMeanNch  ->SetDirectory(0);
      fMeanNpart->SetXTitle("Centrality (%)");
      fMeanNch  ->SetXTitle("Centrality (%)");
      fMeanNpart->SetXTitle(fNpartVsNch->GetXaxis()->GetTitle());
      fMeanNch  ->SetYTitle(fNpartVsNch->GetYaxis()->GetTitle());
      fMeanNpart->SetMarkerStyle(20);
      fMeanNch  ->SetMarkerStyle(21);
      fMeanNpart->SetMarkerColor(kRed +1);
      fMeanNch  ->SetMarkerStyle(kBlue+1);
      fMeanNpart->SetLineColor(kRed +1);
      fMeanNch  ->SetLineStyle(kBlue+1);

      fHelper.CreateHistos(fOutput, 0);
      
      fOutput->Add(fNpartVsNch); 
      fOutput->Add(fNpartB);
      fOutput->Add(fMeanNpart);
      fOutput->Add(fMeanNch);
    }
    /** 
     * Clear internal caches.  Called at start of each event.  Can be
     * overloaded to do some more stuff if needed.
     */
    virtual void Clear(Option_t* option="")
    {
      FastAnalysis::Clear(option);
      fNch = 0;
      fNpart = 0;
    }
    /** 
     * Process the header.  Shall return true if the event is accepted,
     * false otherwise.  Must be overloaded by derived class. 
     * 
     * @return True if the event is to be taken. 
     */
    virtual Bool_t ProcessHeader()
    {
      if (!FastAnalysis::ProcessHeader()) return false;
      fNpart = fHeader->fNtgt + fHeader->fNproj;
      fCent  = GetCentrality();
      fNpartB->Fill(fHeader->fB, fNpart);
      return (fHelper.CheckCentrality(fCent,				      
				      fEventMult,
				      fHeader->fB,
				      fNpart,
				      fHeader->fNbin) >= 0);
    }
    /** 
     * Process a single particle.  
     *
     * @param p Pointer to TParticle object 
     * 
     * @return true if the particle was accepted.  
     */      
    virtual Bool_t ProcessParticle(const TParticle* p)
    {
      Double_t   pT    = p->Pt();
      Double_t   pZ    = p->Pz();
      Double_t   theta = TMath::ATan2(pT, pZ);
      Double_t   eta   = (pT < 1e-10 ? 1024 : -TMath::Log(TMath::Tan(theta/2)));
      if (TMath::Abs(eta) > 0.5) return false;
      fNch++;
      return true;
    }
    /** 
     * Called on each event 
     * 
     * @param entry Entry in chain 
     * 
     * @return true on success, false otherwise 
     */
    virtual Bool_t Process(Long64_t entry)
    {
      if (!FastAnalysis::Process(entry)) return false;
      // Printf("Npart=%4d  Nch=%4d", fNpart, fNch);
						 
      fNpartVsNch->Fill(fNpart, fNch);
      fMeanNpart ->Fill(fCent, fNpart);
      fMeanNch   ->Fill(fCent, fNch);
      return true;
    }
    /** 
     * Final processing.  Scales the histogram to the nubmer of events
     * and the bin width.
     */
    virtual void Terminate()
    {
      fOK = GetEventCount();
      if (fOK <= 0) {
	SetStatus(-1);
	Warning("Terminate", "No events selected");
	return;
      }
      TClass* pcl = TProfile::Class();
      fNpartVsNch      = static_cast<TProfile*>(GetOutputObject("nPartVsnCh",
								pcl));
      fMeanNpart       = static_cast<TProfile*>(GetOutputObject("meanNpart",
								pcl));
      fMeanNch         = static_cast<TProfile*>(GetOutputObject("meanNch",pcl));
      if (!fNpartVsNch || !fMeanNpart || !fMeanNch) {
	SetStatus(-1);
	Warning("Terminate", "No profiles found found (%p,%p,%p)",
		fNpartVsNch, fMeanNpart, fMeanNch);
	return;
      }
      Printf("A total of %ld events", fOK);
      
      TH1*      eg     = static_cast<TH1*>(GetOutputObject("eg",TH1::Class()));

      Color_t c = kYellow+1;
      Style_t s = 1;
      if (eg) {
	TString t(eg->GetTitle());
	t.ToUpper();
	if      (t.Contains("HIJING"))   { c =  kRed     - 3; s = 1; }
	else if (t.Contains("AMPT"))     { c =  kBlue    - 3; s = 2; }
	else if (t.Contains("DPMJET"))   { c =  kMagenta - 3; s = 3; }
	else if (t.Contains("EPOS-LHC")) { c =  kGreen   - 3; s = 4; }
	else if (t.Contains("AMPT2"))    { c =  kMagenta - 3; s = 3; }
	else if (t.Contains("KLN"))      { c =  kBlue    - 3; s = 2; }
	else if (t.Contains("HIJING2"))  { c =  kRed     + 3; s = 1; }
	else if (t.Contains("ANGANTYR")) { c =  kOrange  + 3; s = 3; }
      }

#if 0
      TGraphErrors* res = new TGraphErrors(mNch->GetNbinsX());
      res->SetName("result");
      res->SetTitle(Form("%s \\frac{2}{%s}",
			 fNpartVsNch->GetYaxis()->GetTitle(),
			 fNpartVsNch->GetXaxis()->GetTitle()));
      res->SetMarkerColor(c);
      res->SetFillColor(c);
      res->SetFillStyle(1001);
      res->SetLineColor(c);
      res->SetLineStyle(s);
      for (Int_t i = 1; i <= mNch->GetNbinsX(); i++) {
	Double_t nPart = mNch->GetXaxis()->GetBinCenter(i);
	Double_t ePart = mNch->GetXaxis()->GetBinWidth(i);
	Double_t nCh   = mNch->GetBinContent(i);
	Double_t eCh   = mNch->GetBinError  (i);
	res->SetPoint     (i-1, nPart,   nCh/(nPart/2));
	res->SetPointError(i-1, ePart/2, eCh/(nPart/2));
      }
#endif

      TGraphErrors* ret = new TGraphErrors(fMeanNpart->GetNbinsX());
      ret->SetName("resultB");
      ret->SetTitle(Form("%s \\frac{2}{%s}",
			 fNpartVsNch->GetYaxis()->GetTitle(),
			 fNpartVsNch->GetXaxis()->GetTitle()));
      ret->SetMarkerColor(c);
      ret->SetFillColor(c);
      ret->SetFillStyle(1001);
      ret->SetLineColor(c);
      ret->SetLineStyle(s);
      for (Int_t i = 1; i <= fMeanNpart->GetNbinsX(); i++) {
	Double_t nPart = fMeanNpart->GetBinContent(i);
	Double_t ePart = fMeanNpart->GetBinError  (i);
	Double_t nCh   = fMeanNch  ->GetBinContent(i);
	Double_t eCh   = fMeanNch  ->GetBinError  (i);
	Double_t y     = nCh/(nPart/2);
	Double_t ey    = y*TMath::Sqrt(TMath::Power(eCh/nCh,2)
				       +TMath::Power(ePart/nPart,2));
	ret->SetPoint     (i-1, nPart, y);
	ret->SetPointError(i-1, ePart, ey);
      }
      
      ret->Draw("al3");
      ret->GetHistogram()->SetTitle("");
      ret->GetHistogram()->SetXTitle(fMeanNpart->GetYaxis()->GetTitle());
      ret->GetHistogram()->SetYTitle(Form("%s \\frac{2}{%s}",
					  fMeanNch->GetYaxis()->GetTitle(),
					  fMeanNpart->GetYaxis()->GetTitle()));

      // fOutput->Add(res);
      fOutput->Add(ret);

      fHelper.Finalize(fOutput, 0, 0);
    }
    /** 
     * Get the list of monitored objects 
     * 
     * @return The list of monitored objects 
     */
    virtual TList* GetMonitorObjects()
    {
      TObject* m1 = new TNamed("nPartVsnCh", "");
      TObject* m2 = new TNamed("nPartB", "");
      TObject* m3 = new TNamed("meanNpart", "");
      TObject* m4 = new TNamed("meanNch", "");

      TList* ret = new TList;
      ret->Add(m1);
      ret->Add(m2);
      ret->Add(m3);
      ret->Add(m4);
    
      return ret;
    }
    ClassDef(Base,1);
  };
  /**
   * For the Xe-Xe run at sqrt(sNN)=5.44TeV
   * 
   */
  struct XeXe : public Base
  {
    XeXe(Bool_t verbose=false, Int_t monitor=0)
      : Base(verbose, monitor)
    {
      fTrigMask = 0x0;
    }
    virtual TProfile* CreateOutput()
    {
      return new TProfile("", "Xe-Xe", 10, 0, 260);
#if 0
		      
      Double_t mid[21] = {
	245,
	230,
	215,
	200,
	178,
	152,
	129,
	108,
	90.1,
	74.1,
	60.4,
	48.5,
	38.3,
	29.8,
	22.7,
	16.9,
	12.6,
	9.37,
	6.99,
	5.26,
	3.45 };
      TArrayD bins(22);
      bins[21] = mid[0] + (mid[0]-mid[1])/2;
      for (size_t i = 1; i < 21; i++) 
	bins[21-i] = mid[i] + (mid[i-1]-mid[i]) / 2;
      bins[0] = mid[20] - (mid[19]-mid[20]) / 2;

      return new TH2D("", "Xe-Xe", bins.GetSize()-1, bins.GetArray(),
		      700, 0, 1400);
#endif 
    }
    ClassDef(XeXe,1);
  };
}
      
/** 
 * Create analysis objects
 */
struct MidNchMaker : public FastAnalysis::Maker
{
  MidNchMaker() : FastAnalysis::Maker("MidNch") {}
  FastAnalysis* Make(const TString& subtype,
		     Int_t          monitor,
		     Bool_t         verbose,
		     TMap&          uopt)
  {
    FastAnalysis* ret = 0;
    TString t(subtype);
    if (t.EqualTo("XeXe")) ret = new MidNch::XeXe(verbose, monitor);

    if (ret) {
      TPair* tp = static_cast<TPair*>(uopt.FindObject("trig"));
      if( tp) ret->SetTrigger(tp->Value()->GetName());
    }
    else 
      Printf("Error: MidNchMaker::Make: Invalid spec: %s", t.Data());
    return ret;
    
  }
  void List() const
  {
    Printf("XeXe          - For Xe-Xe collisions at sqrt(sNN)=5.44TeV");
  }
  const char* Script() const { return __FILE__; }
};

MidNchMaker* _midNchMaker = new MidNchMaker;

#ifdef __MAKECINT__
#pragma link C++ nestedclasses;
#pragma link C++ namespace MidNch;
#endif
#endif
//
// EOF
//
