void joinTrees( const char* inlist = "trending.list" )
{

  TChain * chain = AliXRDPROOFtoolkit::MakeChain( inlist,"trending",0,10000,0);
  TFile f("trending.root","recreate");
  TTree * tree = chain->CopyTree("1");
  tree->Write();
  f.Close(); 
}
