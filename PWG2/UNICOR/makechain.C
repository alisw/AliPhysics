// Author: Dariusz Miskowiec <mailto:d.miskowiec@gsi.de> 2008

//=============================================================================
TChain *makechainexp(char *exp, char *path, int nfil=1000000, int nskip=0) {
  char treename[1000];
  if (strcmp(exp,"ceres3c2")==0) sprintf(treename,"T");
  if (strcmp(exp,"aliceesd")==0) sprintf(treename,"esdTree");
  if (strcmp(exp,"cbm")==0) sprintf(treename,"cbmsim");
  return makechain(treename, path, nfil, nskip);
}
//=============================================================================
TChain *makechain(char *name, char *path, int nfil=1000000, int nskip=0) {

  // make chain of root files
  // if path ends with "root" then add all these files; 
  // otherwise interprete path as the list of files and add nfil files 
  // after skipping nskip

  printf("path=%s\n",path);

  TChain *chain = new TChain(name);
  TString str(path);
  if (str.EndsWith("root")) chain->Add(path);
  else {
    fstream ascii_in;
    ascii_in.open(path, ios::in);
    char filnam[1000];
    for (int i=0; i<nfil+nskip; i++) {
      ascii_in >> filnam;
      if (ascii_in.eof()) break;
      if (i>=nskip) chain->Add(filnam);
    }
  }
  chain->Lookup();
  chain->GetListOfFiles()->Print();
  return chain;
}
//=============================================================================
