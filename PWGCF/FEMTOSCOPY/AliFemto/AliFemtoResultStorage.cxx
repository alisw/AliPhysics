///
/// \file AliFemtoResultStorage.cxx
///

#include "AliFemtoResultStorage.h"

#include <iostream>

#include <TObjArray.h>
#include <TFile.h>

/// \cond CLASSIMP
ClassImp(AliFemtoResultStorage);
/// \endcond

AliFemtoResultStorage::AliFemtoResultStorage()
  : TNamed("FemtoResult", "AliFemtoResultStorage")
  , fObjects(new TObjArray())
{
}

AliFemtoResultStorage::AliFemtoResultStorage(const TString &name)
  : TNamed(name.Data(), "FemtoResultStorage")
  , fObjects(new TObjArray())
{
}


AliFemtoResultStorage::AliFemtoResultStorage(const TString &name, TList *)
  : TNamed(name.Data(), "FemtoResultStorage")
  , fObjects(new TObjArray())
{
}

static Int_t recursive_directory_write(const TObject *objects, TDirectory *dir)
{
  Int_t result = 0;
  if (auto *collection = dynamic_cast<const TCollection*>(objects)) {
    TString path = collection->GetName();
    dir->mkdir(path);
    TDirectory *subdir = dir->GetDirectory(path);

    TIter next_object(collection);
    while (TObject *obj = next_object()) {
      result += recursive_directory_write(obj, subdir);
    }
  } else {
    result += dir->WriteTObject(objects);
  }
  return result;
}

Int_t
AliFemtoResultStorage::Write(const char *name, Int_t option, Int_t bufsize)
{
  std::cout << "\n\n[FemtoResultStorage::Write]\n";
  std::cout << "  gDirectory: " << gDirectory->GetName() << "\n";
  std::cout << "  gFile: " << gFile->GetName() << "\n";
  // return const_cast<const AliFemtoResultStorage*>(this)->Write(name, option, bufsize);

  TString path = name == nullptr ? name : fName.Data();
  // TDirectory *original_directory = gDirectory;

  gDirectory->mkdir(path);
  TDirectory *outdir = gDirectory->GetDirectory(path);

  Int_t result = recursive_directory_write(fObjects, outdir);

  // gDirectory = original_directory;
  return result;
}

Int_t
AliFemtoResultStorage::Write(const char *name, Int_t option, Int_t bufsize) const
{
  std::cout << "\n\n[FemtoResultStorage::Write] CONST\n\n";
  return 0;
}
//*/
