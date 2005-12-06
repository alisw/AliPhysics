void
DrawXsection(const char* file="xsec.root", 
	     const char* var="LOSS", 
	     const char* medName="FMD_Si$", 
	     Double_t thick=.03,
	     const char* pdgName="pi+")
{
  TFile*   file = TFile::Open("xsec.root", "READ");
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

  TDatabasePDG* pdgDb = TDatabasePDG::Instance();
  TParticlePDG* pdgP  = pdgDb->GetParticle(pdgName);
  if (!pdgP) {
    std::cerr << "Couldn't find particle " << pdgName << std::endl;
    return;
  }
  Double_t m = pdgP->Mass();
  Double_t q = pdgP->Charge() / 3;
  if (m == 0) {
    std::cerr  << "Mass is 0" << std::endl;
    return;
  }
  
  TGraph* graph = new TGraph(n);
  for (Int_t i = 0; i < n; i++) {
    tree->GetEntry(i);
    graph->SetPoint(i, tkine/m/q/q, value*thick);
  }
  graph->Draw("ALP");
}

//____________________________________________________________________
//
// EOF
//
