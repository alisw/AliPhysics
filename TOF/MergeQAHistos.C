/*
 * this macros connects to ALIEN and searches for TOFQA.root
 * files in the provided search dir. the files are merged in
 * a single local TOFQA.root file.
 * the output file can be used both to check QA histos merged 
 * over more files and to use them to provide reference QA
 * histos merged over more runs (the starting reference histos
 * will be histo from MC production/reconstruction process).
 */

MergeQAHistos(const Char_t *searchDir = "/alice/sim/PDC_09/LHC09a4/80050", Int_t maxFiles = kMaxInt)
{

  const Char_t *fileName = "TOFQA.root";

  TGrid *gGrid = TGrid::Connect("alien");
  if (!gGrid || !gGrid->IsConnected()) {
    AliError("cannot connect to ALIEN");
    return;
  }
  
  TGridResult *gr = gGrid->Query(searchDir, fileName);
  if (gr->GetEntries() < 1) {
    printf("less than one input files: abort\n");
    return;
  }
  printf("%d files found\n", gr->GetEntries());
  
  TFileMerger merger;
  merger.OutputFile(fileName);
  
  Int_t mergedFiles = 0;
  Int_t nFiles = gr->GetEntries();
  for (Int_t i = 0; i < nFiles && i < maxFiles; i++) {
    
    if (merger.AddFile(gr->GetKey(i, "turl")))
      mergedFiles++;

  }

  merger.Merge();
  printf("merged %d files\n", mergedFiles);
  printf("output written on %s\n", fileName);
}

