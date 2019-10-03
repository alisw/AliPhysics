/**
 * @file   YieldExpectations.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Fri Oct  7 12:01:42 2016
 * 
 * @brief  Do back of the envelope calculation of effect of reweighing
 * 
 * @ingroup pwglf_forward_tracklets
 * 
 */
#ifndef __CINT__
# include "AliTrackletAODUtils.C"
# include <vector>
# include <TFile.h>
# include <TError.h>
# include <TCollection.h>
# include <TString.h>
# include <TH1.h>
# include <TH2.h>
# include <THStack.h>
# include <TClass.h>
# include <TBrowser.h>
# include <TCanvas.h>
# include <TLegend.h>
# include <TLegendEntry.h>
# include <TStyle.h>
#else
class TFile;
class TClass;
class TH1;
class TH2;
class TDirectory;
class THStack;
class TCollection;
class TBrowser;
class TCanvas;
class TVirtualPad;
class TLegend;
class TAxis;
#endif
/**
 * 
 * 
 * @ingroup pwglf_forward_tracklets
 */
struct ClusterCalculations
{
#ifndef __CINT__
  typedef AliTrackletAODUtils U;
#endif
  void FixLabels(TAxis* axis)
  {
    for (Int_t i = 1; i <= axis->GetNbins(); i++) {
      TString lbl = axis->GetBinLabel(i);
      Int_t   pdg = lbl.Atoi();
      Color_t d1; Style_t d2;
      U::PdgAttr(pdg, lbl, d1, d2);
      axis->SetBinLabel(i, lbl);
    }
    axis->LabelsOption("v");
    axis->SetLabelSize(0.05);
  }
  void Run(const char* filename="PbPb_5023_LHC15k1a1_245064_gridAOD_stk/root_archive_245064/AnalysisResults.root")
  {
    TFile*        file  = U::OpenFile(filename);
    U::Container* top   = U::GetC(file, "trackletsSums");
    TH1*          seenP = U::GetH1(top, "seenTrackPdg");
    TH1*          usedP = U::GetH1(top, "usedTrackPdg");
    TH2*          seenC = U::GetH2(top, "seenClusterPdg");
    TH2*          usedC = U::GetH2(top, "usedClusterPdg");
    
    THStack* stackP = new THStack("stackP","");
    stackP->Add(seenP);
    stackP->Add(usedP);

    TH1* ratioP = static_cast<TH1*>(seenP->Clone("ratioP"));
    ratioP->Divide(usedP);

    TH1*          seenC1 = seenC->ProjectionX("seenC1",1,1);
    TH1*          usedC1 = usedC->ProjectionX("usedC1",1,1);
    TH1*          seenC2 = seenC->ProjectionX("seenC2",2,2);
    TH1*          usedC2 = usedC->ProjectionX("usedC2",2,2);
    seenC1->SetMarkerStyle(20); seenC1->SetMarkerColor(kGreen+1);
    usedC1->SetMarkerStyle(24); usedC1->SetMarkerColor(kBlue+1);
    seenC2->SetMarkerStyle(21); seenC2->SetMarkerColor(kGreen+3);
    usedC2->SetMarkerStyle(25);	usedC2->SetMarkerColor(kBlue+3);	 
    TH1*          ratio1 = static_cast<TH1*>(seenC1->Clone("ratio1"));
    TH1*          ratio2 = static_cast<TH1*>(seenC2->Clone("ratio2"));    
    ratio1->Divide(usedC1);
    ratio2->Divide(usedC2);
    FixLabels(ratio1->GetXaxis());
    FixLabels(ratio2->GetXaxis());
    
    THStack* stackC = new THStack("stackC","");
    stackC->Add(seenC1);
    stackC->Add(usedC1);
    stackC->Add(seenC2);
    stackC->Add(usedC2);
    
    THStack* ratioC = new THStack("ratioC","");
    ratioC->Add(ratio1);
    ratioC->Add(ratio2);

    Double_t smin = 0.2;
    Double_t smax = 2e11;
    Double_t rmin = 0;
    Double_t rmax = 3.48;
    Int_t    mode = 1; // 0: square, 1: landscape, 2: portrait
    Int_t    cw   = (mode == 1 ? 1600 : mode == 2 ? 800    : 1000);
    Int_t    ch   = (mode == 1 ? cw/2 : mode == 2 ? 1.5*cw : cw);
    TCanvas* c    = new TCanvas("clusterExpectations","Cluster canvas",
				cw,ch);
    c->SetTopMargin(0.01);
    c->SetRightMargin(0.01);
    c->SetBottomMargin(0.13);
    c->Divide(2,2,0,0);

    TVirtualPad* q = c->cd(1);
    q->SetLogy();
    q->SetTicks();
    q->SetGridx();
    q->SetGridy();
    stackP->SetMaximum(smax);
    stackP->SetMinimum(smin);
    stackP->Draw("nostack");
    stackP->GetHistogram()->SetYTitle("#it{N}_{tracks}");

    TLegend* l1 = new TLegend(q->GetLeftMargin()+.5,
			      1-q->GetTopMargin()-.2,
			      1-q->GetRightMargin(),
			      1-q->GetTopMargin(),
			      "Generated tracks");
    l1->SetBorderSize(0);
    l1->SetFillStyle(0);
    TLegendEntry* e = l1->AddEntry("d1", "Before suppression", "p");
    e->SetMarkerStyle(20);
    e = l1->AddEntry("d1", "After suppression", "p");
    e->SetMarkerStyle(24);
    l1->Draw();
    
    q = c->cd(3);
    q->SetTicks();
    q->SetGridx();
    q->SetGridy();
    ratioP->SetMinimum(rmin);
    ratioP->SetMaximum(rmax);
    ratioP->SetStats(0);
    ratioP->SetTitle("");
    ratioP->Draw();
    ratioP->SetYTitle("Before/After");
    FixLabels(ratioP->GetXaxis());

    
    q = c->cd(2);
    q->SetLogy();
    q->SetTicks();
    q->SetGridx();
    q->SetGridy();
    q->SetRightMargin(0.01);
    stackC->SetMinimum(smin);
    stackC->SetMaximum(smax);
    stackC->Draw("nostack");
    stackC->GetHistogram()->SetYTitle("#it{N}_{cluster}");

    TLegend* l2 = new TLegend(q->GetLeftMargin()+.6,
			      1-q->GetTopMargin()-.2,
			      1-q->GetRightMargin(),
			      1-q->GetTopMargin(),
			      "Associated clusters");
    l2->SetBorderSize(0);
    l2->SetFillStyle(0);
    // l2->SetNColumns(2);
    e = l2->AddEntry("d1","Layer 1", "p");
    e->SetMarkerStyle(20);
    e = l2->AddEntry("d1","Layer 2", "p");
    e->SetMarkerStyle(21);
    l2->Draw();
    
		     
    q = c->cd(4);
    q->SetTicks();
    q->SetRightMargin(0.01);
    q->SetGridx();
    q->SetGridy();
    ratioC->SetMinimum(rmin);
    ratioC->SetMaximum(rmax);
    ratioC->Draw("nostack");    
    ratioC->GetHistogram()->SetYTitle("Before/After");
    ratioC->GetHistogram()->GetXaxis()->LabelsOption("v");
    ratioC->GetHistogram()->GetXaxis()->SetLabelSize(ratioP->GetXaxis()
						     ->GetLabelSize());
    // FixLabels(ratioC->GetHistogram()->GetXaxis());

    c->Modified();
    c->Update();
    c->cd();
    c->SaveAs(Form("%s.png",c->GetName()));

    TH1* eff1 = static_cast<TH1*>(seenC1->Clone("eff1"));
    TH1* eff2 = static_cast<TH1*>(seenC2->Clone("eff2"));
    eff1->Divide(seenP);
    eff2->Divide(seenP);
    FixLabels(eff1->GetXaxis());
    FixLabels(eff2->GetXaxis());
    THStack* eff = new THStack("eff","");
    eff->Add(eff1);
    eff->Add(eff2);
    TCanvas* cc = new TCanvas("clusterEfficiency","EffCanvas",cw,ch);
    cc->SetTopMargin(0.01);
    cc->SetRightMargin(0.01);
    cc->SetTicks();
    eff->Draw("nostack");    
    eff->GetHistogram()->SetYTitle("#it{N}_{cluster}/#it{N}_{tracks}");
    eff->GetHistogram()->GetXaxis()->LabelsOption("v");
    eff->GetHistogram()->GetXaxis()->SetLabelSize(ratioP->GetXaxis()
						     ->GetLabelSize());
    
    
    cc->Modified();
    cc->Update();
    cc->cd();
    cc->SaveAs(Form("%s.png",cc->GetName()));
  }        
};

/** 
 * 
 * 
 * @ingroup pwglf_forward_tracklets
 */
void ClusterExpectations(Double_t c1=0, Double_t c2=0)
{
  if (!gROOT->GetClass("AliTrackletAODUtils")) {
    Printf("Loading utilities");    
    gROOT->LoadMacro("$ANA_SRC/dndeta/tracklets3/AliTrackletAODUtils.C+g");
  }
#if 0
  gSystem->AddIncludePath("-DSELF_COMPILE__");
  if (!gROOT->GetClass("ClusterCalculations")) { 
    gInterpreter->ClearFileBusy();
    // gInterpreter->UnloadFile("YieldExpectations.C");
    Printf("Reload self compiled");  
    gROOT->LoadMacro("$ANA_SRC/dndeta/tracklets3/YieldExpectations.C+g");
  }
#endif
  ClusterCalculations c;
  c.Run();
}


