//
// Macro to merge root files
//
void merge(Bool_t isOnGrid=kTRUE,
	   Char_t *outputFileName="AnalysisResults.root",
	   Char_t *inputFileName="listOfFiles.txt")
{

  if (isOnGrid) TGrid::Connect("alien://");

  TFileMerger m;
  m.OutputFile(outputFileName);

  ifstream ftxt(inputFileName);
  Char_t fileName[4096];
  while (!ftxt.eof()) {
    ftxt.getline(fileName,4096);
    if (ftxt.eof()) break;
    m.AddFile(fileName);
  }
  ftxt.close();

  m.Merge();

}
