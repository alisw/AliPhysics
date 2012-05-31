void MergeESDTags()
{
  char spath[2048];
  char sglob[1024];
  sprintf(spath,"");
  sprintf(sglob,"");
  for (int i=0; i< gApplication->Argc();i++){
    if (!(strcmp(gApplication->Argv(i),"--path")))
      sprintf(spath, gApplication->Argv(i+1));
    if (!(strcmp(gApplication->Argv(i),"--glob")))
      sprintf(sglob, gApplication->Argv(i+1));
  }

  if (!strcmp(sglob, "")) sprintf(sglob, "ESD.tag.root");

  printf("*** Connect to AliEn ***\n");
  TGrid::Connect("alien://");
  gSystem->Load("libProofPlayer.so");
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  
  // Create A tag creator object 
  AliESDTagCreator *tagCre = new AliESDTagCreator();
  tagCre->SetStorage(0);

  // Find all the event tag files in the GRID directory
  TGridResult* tagResult = gGrid->Query(spath,"ESD.tag.root");

  // Merge the tags	
  tagCre->MergeTagsForRun("ESD",tagResult);

  return;
}
