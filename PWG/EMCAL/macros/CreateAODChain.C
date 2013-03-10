/* $Id$ */

// This helper macros creates a chain of AOD files for you. 
// Source can be either a text file with the file paths or 
// a directory. In the latter case all AOD files in all 
// subdirectories are considered.
//
// Author: Jan.Fiete.Grosse-Oetringhaus@cern.ch

TChain* CreateAODChain(
  const char* aDataDir = "AODfiles.txt", 
  Int_t aRuns          = 20, 
  Int_t offset         = 0, 
  Bool_t addFileName   = kFALSE, 
  const char* check    = 0)
{
  // creates chain of files in a given directory or file containing a list.
  // In case of directory the structure is expected as:
  // <aDataDir>/<dir0>/AliAOD.root
  // <aDataDir>/<dir1>/AliAOD.root
  // ...
  //
  // if addFileName is true the list only needs to contain the directories that contain the AliAODs.root files
  // if check is != 0 the files that work are written back into the textfile with the name check

  if (!aDataDir)
    return 0;

  Long_t id, size, flags, modtime;
  if (gSystem->GetPathInfo(aDataDir, &id, &size, &flags, &modtime))
  {
    printf("%s not found.\n", aDataDir);
    return 0;
  }

  TChain* chain = new TChain("aodTree");
  TChain* chainFriend = 0;

  if (flags & 2)
  {
    TString execDir(gSystem->pwd());
    TSystemDirectory* baseDir = new TSystemDirectory(".", aDataDir);
    TList* dirList            = baseDir->GetListOfFiles();
    Int_t nDirs               = dirList->GetEntries();
    gSystem->cd(execDir);

    Int_t count = 0;

    for (Int_t iDir=0; iDir<nDirs; ++iDir)
    {
      TSystemFile* presentDir = (TSystemFile*) dirList->At(iDir);
      if (!presentDir || !presentDir->IsDirectory() || strcmp(presentDir->GetName(), ".") == 0 || strcmp(presentDir->GetName(), "..") == 0)
        continue;

      if (offset > 0)
      {
        --offset;
        continue;
      }

      if (count++ == aRuns)
        break;

      TString presentDirName(aDataDir);
      presentDirName += "/";
      presentDirName += presentDir->GetName();

      chain->Add(presentDirName + "/AliAOD.root/aodTree");
    }
  }
  else
  {
    // Open the input stream
    ifstream in;
    in.open(aDataDir);

    ofstream outfile;
    if (check)
      outfile.open(check);

    Int_t count = 0;

    // Read the input list of files and add them to the chain
    TString line;
    while (in.good())
    {
      in >> line;

      if (line.Length() == 0)
        continue;

      if (offset > 0)
      {
        offset--;
        continue;
      }

      if (count++ == aRuns)
        break;

      TString aodFile(line);

      if (addFileName)
        aodFile += "/AliAOD.root";

      if (line.EndsWith(".zip"))
	aodFile += "#AliAOD.root";
 
      if (check)
      {
        TFile* file = TFile::Open(aodFile);
        if (!file)
          continue;
        file->Close();
 
        outfile << line.Data() << endl;
        printf("%s\n", line.Data());
      }        
        
        // add aod file
      chain->Add(aodFile);
    }

    in.close();
    
    if (check)
      outfile.close();
  }

  return chain;
}

void ChainToTextFile(TChain* chain, const char* target)
{
  // write a text list of the files in the chain
  
  TObjArray* list = chain->GetListOfFiles();
  TIterator* iter = list->MakeIterator();
  TObject* obj = 0;

  ofstream outfile;
  outfile.open(target);

  while ((obj = iter->Next())) {
    TString fileName(obj->GetTitle());
    
    outfile << fileName.Data() << endl;
  }

  outfile.close();

  delete iter;
} 

TObjArray* Chain2List(TChain* chain)
{
  // returns a TObjArray of TObjStrings of the file names in the chain

  TObjArray* result = new TObjArray;

  for (Int_t i=0; i<chain->GetListOfFiles()->GetEntries(); i++)
    result->Add(new TObjString(chain->GetListOfFiles()->At(i)->GetTitle()));

  return result;
}

void LookupWrite(TChain* chain, const char* target)
{
  // looks up the chain and writes the remaining files to the text file target

  chain->Lookup();

  ChainToTextFile(chain, target);
}

TChain* CreateChain(const char* treeName, const char* aDataDir, Int_t aRuns, Int_t offset = 0)
{
  // creates chain of files in a given directory or file containing a list.

  if (!treeName || !aDataDir)
    return 0;

  TChain* chain = new TChain(treeName);
  
  // Open the input stream
  ifstream in;
  in.open(aDataDir);

  Int_t count = 0;

  // Read the input list of files and add them to the chain
  TString line;
  while(in.good()) 
  {
    in >> line;
      
    if (line.Length() == 0)
      continue;

    if (line.EndsWith(".zip"))
      line += "#AliAOD.root";
    
    if (offset > 0)
    {
      --offset;
      continue;
    }

    if (count++ == aRuns)
      break;

    chain->Add(line);
  }

  in.close();

  return chain;
}
