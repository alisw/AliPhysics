/**
 * @file   PureMCWeights.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Tue Sep 20 17:18:17 2016
 * 
 * @brief  
 * 
 * 
 * @ingroup pwglf_forward_tracklets
 */
#include "AliTrackletAODUtils.C"
#ifndef __CINT__
#include "AliTrackletWeights.C"
#include <THStack.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#else
class THStack;
class TFile;
class TCanvas;
#endif
/** 
 * A structure to define weights (@f$ p_{T}@f$, centrality, and
 * species dependent) from one simulation to another.
 * 
 * @ingroup pwglf_forward_tracklets
 * @relates AliTrackletBaseWeights
 */
struct PureMCWeights : public AliTrackletAODUtils 
{
  /** 
   * Centrality name 
   * 
   * @param c1 Lower bound 
   * @param c2 Upper bound 
   * 
   * @return The name 
   */
  const char* CentName(Double_t c1, Double_t c2)
  {
    static TString tmp;
    tmp.Form("cent%06.2f_%06.2f", c1, c2);
    tmp.ReplaceAll(".", "d");
    return tmp.Data();
  }
  /** 
   * Run the class 
   * 
   * @param nFile Numerator (truth) file 
   * @param dFile Denominator (correction) file 
   * @param oFile Output file 
   */
  void Run(const char* nFile, const char* dFile, const char* oFile=0) {

    const Color_t cc[] = { kMagenta+2, // 0
			   kBlue+2,    // 1
			   kAzure-1,   // 2 // 10,
			   kCyan+2,    // 3
			   kGreen+1,   // 4 
			   kSpring+5,  // 5 //+10,
			   kYellow+1,  // 6
			   kOrange+5,  // 7 //+10,
			   kRed+1,     // 8
			   kPink+5,    // 9 //+10,
			   kBlack };   // 10
    TString oN(oFile);
    if (oN.IsNull()) oN.Form("%c2%c.root", dFile[0], nFile[0]);
    Printf("********************************************************\n"
	   " Generating pure MC weights\n"
	   "\n"
	   "   Numerator (truth) file:        %s\n"
	   "   Denominator (correction) file: %s\n"
	   "   Output (weight) file:          %s\n",
	   nFile, dFile, oN.Data());
	   

    TFile* nF = TFile::Open(nFile, "READ");
    TFile* dF = TFile::Open(dFile, "READ");
    if (!nF || !dF) return;
    
    Container* nT = GetC(nF, "MidRapidityMCResults");
    Container* dT = GetC(dF, "MidRapidityMCResults");
    if (!nT || !dT) return;

    TH1* nC = GetH1(nT, "cent");
    TH1* dC = GetH1(dT, "cent");
    if (!nC || !dC || !CheckConsistency(nC,dC)) return;

    TList* pdg = new TList;
    TH2D*  pt  = 0;
    pdg->SetOwner();
    
    for (Int_t b = 1; b <= nC->GetNbinsX(); b++) {
      Double_t    c1 = nC->GetXaxis()->GetBinLowEdge(b);
      Double_t    c2 = nC->GetXaxis()->GetBinUpEdge (b);
      const char* cN = CentName(c1, c2);
      Color_t     cC = cc[b%10];
      Container*  nB = GetC(nT, cN);
      Container*  dB = GetC(dT, cN);
      if (!nB || !dB) continue;

      Container*  nG = GetC(nB, "generated");
      Container*  dG = GetC(dB, "generated");
      if (!nG || !dG) continue;

      TH2* nPt2 = GetH2(nG, "etaPt");
      TH2* dPt2 = GetH2(dG, "etaPt");
      if (!nPt2 || !dPt2) continue;

      TH1* nPt  = nPt2->ProjectionY(Form("n%s",cN)); nPt->SetDirectory(0);
      TH1* dPt  = dPt2->ProjectionY(Form("d%s",cN)); dPt->SetDirectory(0);
      TH1* rPt  = static_cast<TH1*>(nPt->Clone(Form("r%s", cN)));
      rPt->Divide(dPt);
      rPt->SetMarkerColor(cC);
      rPt->SetDirectory(0);

      if (!pt) {
	pt = static_cast<TH2D*>(Make2D(0, "pt", "pt weights",
				       kBlack, 20,  *(nC->GetXaxis()),
				       *(nPt->GetXaxis())));
	pt->SetDirectory(0);
      }
      for (Int_t p = 1; p <= rPt->GetNbinsX(); p++) { 
	pt->SetBinContent(b, p, rPt->GetBinContent(p));
	pt->SetBinError  (b, p, rPt->GetBinError  (p));
      }
      // nPt->SetMarkerColor(kRed);  pt->Add(nPt);
      // dPt->SetMarkerColor(kBlue); pt->Add(dPt);
      // pt->Add(nPt2->ProjectionY());
      // pt->Add(dPt2->ProjectionY());
      delete nPt;
      delete dPt;

      TH2* nPdg2 = GetH2(nG, "etaPdg");
      TH2* dPdg2 = GetH2(dG, "etaPdg");
      if (!nPdg2 || !dPdg2) continue;
      if (!CheckConsistency(nPdg2,dPdg2)) continue;

      for (Int_t p = 1; p <= nPdg2->GetNbinsY(); p++) {
	TString sPdg = nPdg2->GetYaxis()->GetBinLabel(p);
	Int_t   iPdg = sPdg.Atoi();
	if (iPdg < 0) continue;

	TH1*    nPdg = nPdg2->ProjectionX("ntmp", p, p);
	TH1*    dPdg = dPdg2->ProjectionX("dtmp", p, p);
	if (!nPdg || !dPdg) continue;

	TH1*    rPdg =
	  static_cast<TH1*>(nPdg->Clone(Form("r%s_%s",cN,sPdg.Data())));
	rPdg->Divide(dPdg);
	rPdg->SetDirectory(0);
	if (rPdg->GetEntries() < 1) {
	  delete nPdg;
	  delete dPdg;
	  delete rPdg;
	  continue;
	}
	
	TF1* fPdg = new TF1("fPdg", "pol0", -5, .5);
	rPdg->Fit(fPdg, "NQR", "", -.5, +.5);

	TH1* hPdg = static_cast<TH1*>(pdg->FindObject(sPdg));
	if (!hPdg) {
	  TString tmp;
	  Color_t c;
	  Style_t s;
	  PdgAttr(iPdg, tmp, c, s);
	  hPdg = Make1D(pdg,sPdg, tmp, c, s, *(nC->GetXaxis()));
	  hPdg->SetBinContent(0, iPdg);
	}
	hPdg->SetBinContent(b, fPdg->GetParameter(0));
	hPdg->SetBinError  (b, fPdg->GetParError (0));

	delete nPdg;
	delete dPdg;
	delete rPdg;
	delete fPdg;
      }

    }
    pt->Draw("colz");
    // pdg->Draw("nostack");

    AliTrackletPtPidStrWeights* w = new AliTrackletPtPidStrWeights("weights");
    w->SetPtWeight(pt);
    TIter next(pdg);
    TH1D*  hPdg = 0;
    while ((hPdg = static_cast<TH1D*>(next()))) {
      Int_t iPdg = hPdg->GetBinContent(0);
      hPdg->SetBinContent(0, 0);
      w->AddAbundanceWeight(iPdg, hPdg);
    }
    TCanvas* c = new TCanvas("c","c");
    w->Draw();

    TFile* out = TFile::Open(oN, "RECREATE");
    w->Write();
    out->Write();
    
  }
};
