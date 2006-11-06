/* $Id$ */

/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// This helper macros creates a chain of ESD files for you. Source can be either a text
// file with the file paths or a directory. In the latter case all ESD files in all subdirectories
// are considered.
//
// Author: Jan.Fiete.Grosse-Oetringhaus@cern.ch

TChain* CreateESDChain(const char* aDataDir = "ESDfiles.txt", Int_t aRuns = 20, Int_t offset = 0, Bool_t addFileName = kFALSE, Bool_t addFriend = kFALSE, Bool_t check = kFALSE)
{
  // creates chain of files in a given directory or file containing a list.
  // In case of directory the structure is expected as:
  // <aDataDir>/<dir0>/AliESDs.root
  // <aDataDir>/<dir1>/AliESDs.root
  // ...
  //
  // if addFileName is true the list only needs to contain the directories that contain the AliESDs.root files
  // if addFriend is true a file AliESDfriends.root is expected in the same directory and added to the chain as friend

  if (!aDataDir)
    return 0;

  Long_t id, size, flags, modtime;
  if (gSystem->GetPathInfo(aDataDir, &id, &size, &flags, &modtime))
  {
    printf("%s not found.\n", aDataDir);
    return 0;
  }

  TChain* chain = new TChain("esdTree");
  TChain* chainFriend = 0;
  
  if (addFriend)
    chainFriend = new TChain("esdFriendTree");

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

      chain->Add(presentDirName + "/AliESDs.root/esdTree");
    }
  }
  else
  {
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
      
      if (offset > 0)
      {
        --offset;
        continue;
      }

      if (count++ == aRuns)
        break;

      TString esdFile(line);

      if (addFileName)
        esdFile += "/AliESDs.root";
        
      TString esdFileFriend(esdFile);
      esdFileFriend.ReplaceAll("AliESDs.root", "AliESDfriends.root");
        
      if (check)
      {
        TFile* file = TFile::Open(esdFile);
        if (!file)
          continue;
        file->Close();
        
        if (chainFriend)
        {
          TFile* file = TFile::Open(esdFileFriend);
          if (!file)
            continue;
          file->Close();
        }
        
        printf("%s\n", line.Data());
      }        
        
        // add esd file
      chain->Add(esdFile);

        // add file
      if (chainFriend)
        chainFriend->Add(esdFileFriend);
    }

    in.close();
  }
  
  if (chainFriend)
    chain->AddFriend(chainFriend);

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

  while ((obj = iter->Next()))
    outfile << obj->GetTitle() << endl;

  outfile.close();

  delete iter;
} 

void LookupWrite(TChain* chain, const char* target)
{
  // looks up the chain and writes the remaining files to the text file target

  chain->Lookup();

  ChainToTextFile(chain, target);
}
