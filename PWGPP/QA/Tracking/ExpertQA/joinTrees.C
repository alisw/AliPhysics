
void joinTrees( const char* inlist = "qpt.list" ){

  TChain * chain = AliXRDPROOFtoolkit::MakeChain( inlist,"TrendingTree",0,10000,0);
TFile f("TrendingTree.root","recreate");
TTree * tree = chain->CopyTree("1");
tree->Write();
f.Close(); 

}
