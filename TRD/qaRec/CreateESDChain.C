TChain* CreateESDChain(const char* filename = "ESDfiles.txt", Int_t nfiles=-1 )
{
  // Create the chain
  TChain* chain = new TChain("esdTree");

  if(!filename){
    chain->Add(Form("%s/AliESDs.root", gSystem->pwd()));
    return chain;
  }


  // read ESD files from the input list.
  ifstream in;
  in.open(filename);
  TString esdfile;
  while(in.good() && (nfiles--) ) {
    in >> esdfile;
    if (!esdfile.Contains("root")) continue; // protection
    chain->Add(esdfile.Data());
  }

  in.close();

  return chain;
}
