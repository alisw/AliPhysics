/** 
 * 
 * 
 * @param fname1 
 * @param fname2 
 *
 * @ingroup pwg2_forward_scripts_corr
 */
void
CompELossFits(const char* fname1, const char* fname2)
{

  TFile* file1 = TFile::Open(fname1, "READ");
  if (!file1) { 
    Error("CompELossFits", "Couldn't open %s", fname1);
    return;
  }

  TFile* file2 = TFile::Open(fname2, "READ");
  if (!file2) { 
    Error("CompELossFits", "Couldn't open %s", fname2);
    return;
  }

  AliFMDCorrELossFit* fit1 = 
    static_cast<AliFMDCorrELossFit*>(file1->Get("elossfits"));
  if (!fit1) { 
    Error("CompELossFits", "Couldn't get elossfits from %s", fname1);
    return;
  }

  AliFMDCorrELossFit* fit2 = 
    static_cast<AliFMDCorrELossFit*>(file2->Get("elossfits"));
  if (!fit2) { 
    Error("CompELossFits", "Couldn't get elossfits from %s", fname2);
    return;
  }


  TList* stacks1 = fit1->GetStacks(true, false, 4);
  TList* stacks2 = fit2->GetStacks(true, false, 4);

  Int_t nStacks = stacks1->GetEntries();

  TCanvas* c = new TCanvas("c", "c", 900, 1200);
  c->SetRightMargin(0.02);
  c->SetTopMargin(0.02);
  c->SetFillColor(0);
  c->SetBorderSize(0);
  c->SetBorderMode(0);

  c->cd();
  TPad* top = new TPad("top", "Top", 0, .95, 1, 1, 0, 0, 0);
  top->Draw();
  top->cd();
  TLatex* l = new TLatex(.5,.5, Form("%s / %s", fname1, fname2));
  l->SetTextSize(0.3);
  l->SetNDC();
  l->SetTextAlign(22);
  l->Draw();

  c->cd();
  TPad* body = new TPad("body", "body", 0, 0, 1, .95, 0, 0, 0);
  body->Draw();
  body->cd();
  body->Divide(2, (nStacks+1)/2, 0, 0);

  Int_t nPad2 = nStacks;
  for (Int_t i = 0; i < nStacks; i++) {
    Int_t iPad = 1 + i/nPad2 + 2 * (i % nPad2);
    TVirtualPad* p = body->cd(i+1);
    p->SetLeftMargin(0.15);
    p->SetRightMargin(0.01);
    THStack* stack1 = static_cast<THStack*>(stacks1->At(i));
    THStack* stack2 = static_cast<THStack*>(stacks2->At(i));

    THStack* ratio  = static_cast<THStack*>(stack1->Clone());
    Int_t    nHists = stack1->GetHists()->GetEntries();
    for (Int_t j = 0; j < nHists; j++) { 
      TH1* h1 = static_cast<TH1*>(stack1->GetHists()->At(j));
      TH1* h2 = static_cast<TH1*>(stack2->GetHists()->At(j));
      TH1* hr = static_cast<TH1*>(ratio->GetHists()->At(j));
      hr->Divide(h1, h2);
    }
    ratio->Draw("nostack");
  }
}

  


  

  
//
// EOF
//
