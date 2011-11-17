/**
 * @file   DrawCuts.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Thu Nov 17 11:20:22 2011
 * 
 * @brief  Draw the cuts used in the analysis
 * 
 * @ingroup pwg2_forward_scripts_qa
 * 
 */
/** 
 * Draw cuts used in analysis
 * 
 * @param filename Input file name 
 *
 * @ingroup pwg2_forward_scripts_qa
 */
void
DrawCuts(const char* filename="forward.root")
{
  gStyle->SetPalette(1);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetTitleW(.4);
  gStyle->SetTitleH(.1);
  gStyle->SetTitleColor(0);
  gStyle->SetTitleStyle(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleX(.6);

  TFile* file = TFile::Open(filename, "READ");
  if (!file) { 
    Error("DrawCuts", "failed to open %s", filename);
    return;
  }

  TList* forward = static_cast<TList*>(file->Get("Forward"));
  if (!forward) { 
    Error("DrawCuts", "List Forward not found in %s", filename);
    return;
  }

  TList* dc = static_cast<TList*>(forward->FindObject("fmdDensityCalculator"));
  if (!dc) { 
    Error("DrawCuts", "List fmdDensityCalculator not found in Forward");
    return;
  }
  TList* sf = static_cast<TList*>(forward->FindObject("fmdSharingFilter"));
  if (!dc) { 
    Error("DrawCuts", "List fmdSharingFilter not found in Forward");
    return;
  }
  TList* hc = static_cast<TList*>(forward->FindObject("fmdHistCollector"));
  if (!hc) { 
    Error("DrawCuts", "List fmdHistCollector not found in Forward");
    return;
  }
  TH2* hC = static_cast<TH2*>(sf->FindObject("highCuts"));
  if (!hC) { 
    Error("DrawCuts", "Histogram highCuts found in %s", sf->GetName());
    return;
  }
  TH2* lC = static_cast<TH2*>(dc->FindObject("lowCuts"));
  if (!lC) { 
    Error("DrawCuts", "Histogram lowCuts found in %s", dc->GetName());
    return;
  }
  TH2* co = static_cast<TH2*>(hc->FindObject("coverage"));
  if (!co) { 
    Error("DrawCuts", "Histogram coverage found in %s", hc->GetName());
    return;
  }
  TCanvas* c = new TCanvas("cuts", "Cuts used in the analysis", 900, 700);
  c->SetFillColor(0);
  c->SetBorderSize(0);
  c->Divide(3,1);
  
  c->cd(1); hC->Draw("colz");
  c->cd(2); lC->Draw("colz");
  c->cd(3); co->Draw("colz");
  c->cd();

}

  
  
 
