TChain* CreateESDChain(const char* filename = "ESDfiles.txt", Int_t nfiles=-1 )
{
  // Create the chain
  TChain* chain = new TChain("esdTree");

  // Open the input stream
  ifstream in;
  in.open(filename);

  // Read the input list of files and add them to the chain
  TString esdfile;
  while(in.good() && (nfiles--) ) {
    in >> esdfile;
    if (!esdfile.Contains("root")) continue; // protection

    chain->Add(esdfile.Data());
  }

  in.close();

  return chain;
}
