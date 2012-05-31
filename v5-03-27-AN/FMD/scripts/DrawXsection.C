//____________________________________________________________________
//
// Script to draw a X-section, LOSS, or range made with MakeXsection
//
/** @ingroup xsec_script
    @param scale 
    @param filename 
    @param var 
    @param medName 
    @param thick 
    @param pdgName 
*/
void
DrawXsection(Bool_t scale=kFALSE, 
	     const char* filename="xsec.root", 
	     const char* var="LOSS", 
	     const char* medName="FMD_Si$", 
	     Double_t thick=.03,
	     const char* pdgName="pi+")
{
  TFile*   file = TFile::Open(filename, "READ");
  TTree*   tree = static_cast<TTree*>(file->Get(Form("%s_%s",medName,
						     pdgName)));
  TLeaf* tb   = tree->GetLeaf("T");
  TLeaf* vb   = tree->GetLeaf(var);
  if (!vb) {
    std::cerr << "Leaf " << var << " not found" << std::endl;
    return;
  }
  Float_t tkine, value;
  tb->SetAddress(&tkine);
  vb->SetAddress(&value);
  Int_t n = tree->GetEntries();

  Float_t xscale = 1;
  Float_t yscale = 1;
  if (scale) {
    TDatabasePDG* pdgDb = TDatabasePDG::Instance();
    TParticlePDG* pdgP  = pdgDb->GetParticle(pdgName);
    if (!pdgP) {
      std::cerr << "Couldn't find particle " << pdgName << std::endl;
      return;
    }
    Double_t m = pdgP->Mass();
    Double_t q = pdgP->Charge() / 3;
    if (m == 0 || q == 0) {
      std::cerr  << "Mass is 0" << std::endl;
      return;
    }
    xscale = 1 / m;
    yscale = 1 / (q * q);
  }
  
  TGraphErrors* graph = new TGraphErrors(n);
  for (Int_t i = 0; i < n; i++) {
    tree->GetEntry(i);
    Double_t x = tkine*xscale;
    Double_t y = value*yscale;
    graph->SetPoint(i, x, y); 
    // 5 sigma
    graph->SetPointError(i, 0, 5 * .1 * y);
  }
  TCanvas* c = new TCanvas("c","c");
  c->SetLogx();
  c->SetLogy();
  graph->SetLineWidth(2);
  graph->SetFillStyle(3001);
  graph->SetFillColor(6);
  graph->Draw("L");
  graph->DrawClone("AL3");
  c->Modified();
  c->Update();
  c->cd();
  c->SaveAs("xsec.C");
  
}

//____________________________________________________________________
//
// EOF
//
