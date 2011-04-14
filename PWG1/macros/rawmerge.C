//
// Macro to create the "raw" data file with selected events
// Paramaters:
//  eventListFileName - input ascii file with list of chunks and event numbers - two collumns
//  osplit            - split the file after given threhshold value
//
void rawmerge(const char *eventListFileName,
              const char *outputDirectoryURI,
              Long64_t osplit=-1)
{
  Int_t eventNumber;
  FILE *files=fopen(eventListFileName,"r");
  if (!files) {
    fprintf(stderr,"error: could not read event list file \"%s\". Exiting.\n",eventListFileName);
    return;
  }
  char iURI[1000];
  char oURI[1000];
  TFile *ifile=0;
  TFile *ofile=0;
  TTree *itree=0;
  TTree *otree=0;
  Long64_t ievent;
  Long64_t oevent;
  Int_t ofilenumber=0;
  Int_t line=0;
  while (!feof(files)) {
    delete ifile;ifile=0;
    ++line;
    if (fscanf(files,"%s %d\n",iURI,&ievent)!=2) {
      fprintf(stderr,"warning: corrupted event line (%d) in input file, skipping it...\n",line);
      continue;
    }
    printf("> processing \"%s\" event %d...\n",iURI,ievent);
    TFile *ifile=TFile::Open(iURI);
    if (!ifile) {
      fprintf(stderr,"warning: could not open file for event \"%s\", skipping it...\n",iURI);
      continue;
    }
    TTree *itree=dynamic_cast<TTree*>(ifile->Get("RAW"));
    if (!itree) {
      fprintf(stderr,"warning: could not find RAW tree for event \"%s\", skipping it...\n",iURI);
      continue;
    }

    // create (new) output file and tree
    if (!ofile || (osplit>0 && oevent%osplit==0)) {
      delete ofile;
      ++ofilenumber;
      sprintf(oURI,"%s/merged_%d.root",outputDirectoryURI,ofilenumber);
      printf("< creating output file \"%s\"\n",oURI);
      ofile=TFile::Open(oURI,"RECREATE");
      if (!ofile) {
        fprintf(stderr,"error: could not create output file: \"%s\" Exiting.\n",oURI);
        break;
      }
      otree=itree->CloneTree(0);
    }

    // copy event and write to file
    otree->CopyAddresses(itree);
    itree->GetEntry(ievent);
    otree->Fill();
//    otree->CopyEntries(itree,Form("Entry$==%d",ievent),1);
    ofile->Write();
    ++oevent;

    // reset input
    itree->ResetBranchAddresses();
  }

  printf("Merged %d events.\n",oevent);
  delete ifile;
  delete ofile;
}

