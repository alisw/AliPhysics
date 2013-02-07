/** 
 * Create a file collection in a ROOT file 
 * 
 * @param dir         Input directory
 * @param tN          Tree name
 * @param pa          File name pattern
 * @param mc          If true, simulations
 * @param recursive   Scan recursively 
 *
 * @ingroup pwglf_forward_trains_helper
 */
void
CreateFileCollection(const TString& dir="/data/alice/data/ppb/LHC12g/pass1/188359/",
		     const TString& tN="esdTree", 
		     const TString& pa="AliESD", 
		     Bool_t mc=false, 
		     Bool_t recursive=false)
)
{
  gROOT->LoadMacro("ChainBuilder.C+");
  
  UShort_t type = ChainBuilder::CheckSource(dir);
  Info("", "type=%d", type);
  TChain* chain = ChainBuilder::Create(dir, tN, pa, mc, recursive);
  if (!chain) { 
    Error("CreateFileCollection", "Failed to make chain");
    return;
  }
  Int_t port;
  TString host;
  { 
    TUrl u(Form("root://%s//foo", gSystem->HostName()));
    port = u.GetPort() * 10;
    host = u.GetHostFQDN();
  }
  

  TFileCollection* fc  = new TFileCollection("files");
  TObjArray*       cEs = chain->GetListOfFiles();
  TChainElement*   cE  = 0;
  TIter            next(cEs);
  while ((cE= static_cast<TChainElement*>(next()))) {
    TString fN(cE->GetTitle());
    TFile* f = TFile::Open(fN, "READ");
    TTree* t = static_cast<TTree*>(f->Get(tN));
    
    fN.Prepend(Form("root://%s:%d/", host.Data(), port));
    TFileInfo* fi = new TFileInfo(Form("%s tree:%s,%d", 
				       fN.Data(), tN.Data(), t->GetEntries()),
				  f->GetSize());
    f->Close();
    fc->Add(fi);
  }
  fc->Print("F");
  
  TFile* files = TFile::Open("files.root", "RECREATE");
  fc->Write();
  files->Close();
  
}
//
// EOF
//
