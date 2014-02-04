merger(const Char_t *filelist, const Char_t *output, Bool_t alien = kFALSE)
{

  gSystem->Setenv("TMPDIR", gSystem->pwd());
  if (alien) TGrid::Connect("alien");

  ifstream filein(filelist);
  Int_t nfiles = 0;
  Char_t filename[4096];
  TFileMerger m(kFALSE);
  m.OutputFile(output);

  while(1) {
    filein.getline(filename, 4096);
    if (filein.eof()) break;
    if (alien) sprintf(filename, "alien://%s", filename);
    printf("adding file: %s\n", filename);
    if (m.AddFile(filename)) {
      printf("file %s successfully added\n", filename);
      nfiles++;
    }
  }

  printf("start merging %d files\n", nfiles);
  m.Merge();
  printf("%d files merged: %s\n", nfiles, output);

  /* create dummy file to tell we are done */
  gSystem->Exec("touch done");
}
