///
/// \file AliFemtoResultStorage.cxx
///

#include "AliFemtoResultStorage.h"
#include "AliFemtoManager.h"

#include <TFile.h>
#include <TObjArray.h>

#include <iostream>


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


AliFemtoResultStorage::AliFemtoResultStorage(const TString &name, TList *list)
  : TNamed(name.Data(), "FemtoResultStorage")
  , fObjects(new TObjArray())
{
  fObjects->Add(list);
}

AliFemtoResultStorage::AliFemtoResultStorage(const TString &name, AliFemtoManager &mgr)
  : AliFemtoResultStorage(name)
{
  for (auto &analysis : *mgr.AnalysisCollection()) {
    TList *output_list = analysis->GetOutputList();
    fObjects->Add(output_list);
  }
}

AliFemtoResultStorage::~AliFemtoResultStorage()
{
  delete fObjects;
}

static Int_t recursive_directory_write(const TObject *objects, TDirectory &dir)
{
  Int_t result = 0;
  if (auto *collection = dynamic_cast<const TCollection*>(objects)) {
    TString path = collection->GetName();
    dir.mkdir(path);
    TDirectory *subdir = dir.GetDirectory(path);

    TIter next_object(collection);
    while (TObject *obj = next_object()) {
      result += recursive_directory_write(obj, *subdir);
    }
  } else {
    result += dir.WriteTObject(objects);
  }
  return result;
}

Int_t
AliFemtoResultStorage::Write(const char *name, Int_t option, Int_t bufsize)
{
  return const_cast<const AliFemtoResultStorage*>(this)->Write(name, option, bufsize);
}

Int_t
AliFemtoResultStorage::Write(const char *name, Int_t option, Int_t bufsize) const
{
  TString path(name);

  // override if 'path' is empty
  if (path.IsWhitespace()) {
    path = fName.Strip(TString::kBoth);
  }

  gDirectory->mkdir(path);
  TDirectory *outdir = gDirectory->GetDirectory(path);
  if (!outdir) {
    Error("AliFemtoResultStorage::Write", "Could not create path %s", path.Data());
    return 0;
  }

  Int_t result = 0;
  for (TObject *obj : *fObjects) {
    auto *output_object_list = dynamic_cast<TList*>(obj);
    if (!output_object_list) {
      AliWarning(Form("Unexpected type '%s' in output list (name: '%s'). Skipping.",
                      obj->ClassName(), obj->GetName()));
      continue;
    }

    for (TObject *output_object : *output_object_list) {
      result += recursive_directory_write(output_object, *outdir);
    }
  }

  return result;
}
