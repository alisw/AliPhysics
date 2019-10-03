#include <fstream>

TChain* MakeAODInputChain(const char* fileListName, int nFiles)
{
  // Create AOD Chain from a list of files.

  TChain *chainAOD = new TChain("aodTree");
  // TChain *chainAODfriend = new TChain("aodTree"); // No support for as of yet.
 
  std::ifstream flstream(fileListName);
  if(!flstream) {
    ::Error("MakeAODInputChain.C", Form("Error opening file %s", fileListName));
    return 0;
  }
  int line = 0;
  while(!flstream.eof()){
    if(nFiles > 0 && line > nFiles)
      break;
    
    const char fileName[256] = "";
    flstream >> fileName;
    Printf("adding %s, ", fileName);
    chainAOD->Add(fileName);
    //chainAODfriend->Add(aodHFname.Data());
  }

  flstream.close();
  return chainAOD;
}
