//e.g. fileName = AnalysisResults-HERWIG7.root -> New file is then AnalysisResults-HERWIG7_commonStructure.root

void AdaptFileStructure(TString fileName)
{
  TFile* f = TFile::Open(fileName.Data(), "READ");
  TList* cList = new TList();
  cList->SetName("cList");
  cList->SetOwner(kTRUE);
  for (Int_t i = 0; i < gFile->GetNkeys(); i++) {
    TObject* obj = f->GetListOfKeys()->At(i);
    TString name = obj->GetName();
    TObject* obj2 = f->Get(name.Data());
    cList->Add(obj2); 
  }
  
  TString fileNameOutput = fileName;
  fileNameOutput.ReplaceAll(".root", "_commonStructure.root");
  TFile* fOut = TFile::Open(fileNameOutput.Data(), "RECREATE");
  fOut->mkdir("PWGLF_PPVsMultXCheck_MC");
  fOut->cd("PWGLF_PPVsMultXCheck_MC");
  cList->Write("cList", TObject::kSingleKey);
  fOut->Close();
  f->Close();
}
