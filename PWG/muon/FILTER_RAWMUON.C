#include "TTree.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TString.h"
#include "TGrid.h"
#include "TLeaf.h"
#include "TError.h"
#include "TFile.h"

void DisableBranches(TTree* tree)
{
  TObjArray* list = tree->GetListOfLeaves();
  TIter next(list);
  TLeaf* leaf;
  
  TObjArray branchesToKeep;
  branchesToKeep.SetOwner(kTRUE);
  
  branchesToKeep.Add(new TObjString("rawevent"));
  branchesToKeep.Add(new TObjString("MUON"));
  branchesToKeep.Add(new TObjString("Common"));
  branchesToKeep.Add(new TObjString("TRG"));
  branchesToKeep.Add(new TObjString("T0"));
  branchesToKeep.Add(new TObjString("VZERO"));
  branchesToKeep.Add(new TObjString("ZDC"));
  branchesToKeep.Add(new TObjString("ITSSPD"));
  
  TIter nit(&branchesToKeep);

  tree->SetBranchStatus("*",0);

  Bool_t on(kTRUE);

  while ( ( leaf = static_cast<TLeaf*>(next()) ) )
  {
    TString name(leaf->GetName());
    if (!name.BeginsWith("f") && name != "rawevent" )
    {
      TObjString* str;
      
      nit.Reset();
      
      on = kFALSE;
      
      while ( ( str = static_cast<TObjString*>(nit()) ) )
      {
        if ( name.BeginsWith(str->String()) )
        {
          on = kTRUE;
        }
      }
    }
    if ( on )
    {
      leaf->GetBranch()->SetStatus(1);
    }
  }
}

int FILTER_RAWMUON(const char* from, const char* to)
{
  TString sinputfile(from);
  
  if (sinputfile.BeginsWith("alien://"))
  {
    if (!gGrid)
    {
      TGrid::Connect("alien://");
    }
    if (!gGrid)
    {
      Error("FILTER_RAWMUON","Cannot connect to Grid !");
      return -1;
    }
  }
  
  TFile* rawFile = TFile::Open(from);
  
  TTree *rawTree=(TTree *)rawFile->Get("RAW");
  if(!rawTree)
  {
    Error("FILTER_RAWMUON","Error getting RAW tree from file %s",from);
    return -2;
  }
  
  TString snewfile(to);
  
  DisableBranches(rawTree);
  
  TFile* newfile =  new TFile(snewfile.Data(),"recreate");
  
  TTree* newTree = rawTree->CloneTree();
  
  newTree->Print();
  
  newTree->AutoSave();
  
  newfile->Close();
  
  delete newfile;

  delete rawFile;
  
  return 0;
}