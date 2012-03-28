//___________________________________________________________________________

LoadLibraries()
{
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  
}

//___________________________________________________________________________

MakeOCDB(const Char_t *filename = "AliESDfriends_v1.root", const Char_t *dbString = "raw://", Int_t runNB)
{
  LoadLibraries();
  AliMeanVertexPreprocessorOffline meanVertexCalib;
  meanVertexCalib.ProcessOutput(filename, dbString, runNB);
}


