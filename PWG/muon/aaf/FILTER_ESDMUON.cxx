#include "TTree.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TString.h"
#include "TGrid.h"
#include "TLeaf.h"
#include "TError.h"
#include "TFile.h"
#include "Riostream.h"

namespace AAF {

void DisableESDBranches(TTree* tree)
{
  TString branchesToKeep[] = {"MuonTracks","MuonClusters","MuonPads","AliESDRun.","AliESDHeader.","AliMultiplicity.","AliESDFMD.","AliESDVZERO.","AliESDTZERO.","SPDVertex.","PrimaryVertex.","AliESDZDC.","SPDPileupVertices"};

  tree->SetBranchStatus("*",0);

  Int_t nBranches = sizeof(branchesToKeep)/sizeof(branchesToKeep[0]);
  for ( Int_t ibranch=0; ibranch<nBranches; ibranch++ ) {
    tree->SetBranchStatus(Form("*%s*",branchesToKeep[ibranch].Data()), 1);
  }
}

int FILTER_ESDMUON(const char* from, const char* to)
{
  std::cout << "FILTER_ESDMUON(" << from << "," << to << ")" << std::endl;

  TString sinputfile(from);

  if (sinputfile.BeginsWith("alien://"))
  {
    if (!gGrid)
    {
      TGrid::Connect("alien://");
    }
    if (!gGrid)
    {
      Error("FILTER_ESDMUON","Cannot connect to Grid !");
      return -1;
    }
  }

  TFile* esdFile = TFile::Open(from);

  if (!esdFile || !esdFile->IsOpen())
  {
    Error("FILTER_ESDMUON","Cannot open input file %s",from);
    return -3;
  }

  TTree *esdTree=(TTree *)esdFile->Get("esdTree");
  if(!esdTree)
  {
    Error("FILTER_ESDMUON","Error getting ESD tree from file %s",from);
    return -2;
  }

  TString snewfile(to);

  DisableESDBranches(esdTree);

  TFile* newfile =  new TFile(snewfile.Data(),"recreate");

  TTree* newTree = esdTree->CloneTree();

  newTree->Print();

  newTree->AutoSave();

  newfile->Close();

  delete newfile;

  delete esdFile;

  return 0;
}

}
