#include "AliRawVEvent.h"
#include "Riostream.h"
#include "TFile.h"
#include "TString.h"
#include "TTree.h"
#include <cstdlib>
#include <sstream>
#include <memory>

// TTree* tree = (TTree*) fFile->Get("RAW");
// if (!tree) {
//   Error("AliRawReaderRoot", "no raw data tree found");
//   fIsValid = kFALSE;
//   return;
// }
// fBranch = tree->GetBranch("rawevent");
// if (!fBranch) {
//   Error("AliRawReaderRoot", "no raw data branch found");
//   fIsValid = kFALSE;
//   return;
// }
//
// fBranch->SetAddress(&fEvent);

Bool_t CopyEntries(const char* filename, const std::vector<int>& entries, TTree*& newtree, TFile* file)
{
  std::cout << "CopyEntries from " << filename << std::endl;

  std::unique_ptr<TFile> input(TFile::Open(filename));

  TTree* oldtree = static_cast<TTree*>(input->Get("RAW"));

  if (!oldtree) {
    std::cout << "Could not get RAW tree from " << filename << std::endl;
    return kFALSE;
  }
  AliRawVEvent* event{ nullptr };
  TBranch* branch = oldtree->GetBranch("rawevent");

  branch->SetAddress(&event);

  if (!newtree) {
    file->cd();
    newtree = oldtree->CloneTree(0);
  }
  else {
    oldtree->CopyAddresses(newtree);
  }

  for (auto i : entries) {
    oldtree->GetEntry(i);
    newtree->Fill();
  }

  newtree->AutoSave();
  return kTRUE;
}

Int_t ApplyEntryList(const char* entryListFileName, const char* outputFileName)
{
  /// Copy all the entries found in entryListFileName to a new
  /// tree in outputFileName
  /// The format of entryListFileName is :
  /// filename1
  /// n1
  /// id1 id2 ... idn1
  /// filename2
  /// n2
  /// id2_1 id2_2 .... id2_n2
  /// where n is the number of entries to be copied, and the ids of the entries
  /// to be copied are on the next line

  std::unique_ptr<TFile> file(TFile::Open(outputFileName, "RECREATE"));

  TTree* newtree = nullptr;

  std::ifstream in(entryListFileName);
  std::string filename;
  std::string numberofentries;
  std::string entries;

  while (std::getline(in, filename)) {

    std::getline(in, numberofentries);
    std::getline(in, entries);

    std::istringstream s(entries);
    int x;
    std::vector<int> v;
    while (s >> x) {
      v.push_back(x);
    }

    CopyEntries(filename.c_str(), v, newtree, file.get());
  }

  return 0;
}

int main(int argc, char** argv)
{
  if (argc < 2) {
    std::cout << "Usage : " << argv[0] << " --entrylist filename.txt --output destination" << std::endl;
    return 1;
  }

  TString output{ "" };
  TString entryList{ "" };

  for (int i = 0; i < argc; ++i) {
    if (TString(argv[i]) == "--entrylist") {
      entryList = argv[i + 1];
    }
    if (TString(argv[i]) == "--output") {
      output = argv[i + 1];
    }
  }

  return ApplyEntryList(entryList.Data(), output.Data());
}
