// Make a summary tree
// Use guiTime with summary file  - it is faster as with chain
// guiTime calibTimeSummary.root 


void MakeChain(const char *prefix){
  gSystem->Exec(Form("find %s |grep calibTree | grep root > calib.list",prefix));
  AliXRDPROOFtoolkit::FilterList("calib.list","* dcs",1);
  AliXRDPROOFtoolkit toolkit;
  TChain * chain = toolkit.MakeChain("calib.list.Good","dcs",0,2000);
  chain->Lookup();
  TTree * tree = chain->CopyTree("1");
  TFile f("calibTimeSummary.root","recreate");
  tree->Write("dcs");
  f.Close();
}
