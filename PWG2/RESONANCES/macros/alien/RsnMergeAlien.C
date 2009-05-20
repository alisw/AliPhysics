void RsnMergeAlien()
{
  // evaluate necessary arguments
  char   searchpath[1024];
  char   searchpattern[1024];
  char   outputfile[1024];
  char   username[1024];
  char   masterjob[1024];
  for (int i=0; i< gApplication->Argc();i++)
  {
    if (!(strcmp(gApplication->Argv(i), "--path")))       sprintf(searchpath, "%s", gApplication->Argv(i+1));
    if (!(strcmp(gApplication->Argv(i), "--name")))       sprintf(searchpattern,"%s", gApplication->Argv(i+1));
    if (!(strcmp(gApplication->Argv(i), "--out")))        sprintf(outputfile, "%s", gApplication->Argv(i+1));
    if (!(strcmp(gApplication->Argv(i), "--user")))       sprintf(username, "%s", gApplication->Argv(i+1));
    if (!(strcmp(gApplication->Argv(i), "--masterjob")))  sprintf(masterjob, "%s", gApplication->Argv(i+1));
  }

  // connect to grid
  TGrid::Connect("alien://");

  // query the grid to find all files in a given pattern
  TGridResult* result = gGrid->Query(searchpath, searchpattern);

  // initialize the file merger
  TFileMerger merger;
  merger.OutputFile(outputfile);

  // loop on query result to merge
  Int_t nMerged = 0;
  while (result->GetKey(nMerged, "turl"))
  {
    merger.AddFile(result->GetKey(nMerged, "turl"));
    nMerged++;
  }
  if (nMerged) merger.Merge();
}
