
#include "AliForwardFlowResultStorage.h"
#include <TFile.h>
#include <TObjArray.h>
#include <iostream>


/// \cond CLASSIMP
ClassImp(AliForwardFlowResultStorage);
/// \endcond

AliForwardFlowResultStorage::AliForwardFlowResultStorage()
  : TNamed("ForwardFlowResult", "AliForwardFlowResultStorage")
  , fObjects(new TObjArray())
{
}

AliForwardFlowResultStorage::AliForwardFlowResultStorage(const TString &name)
  : TNamed(name.Data(), "ForwardFlowResultStorage")
  , fObjects(new TObjArray())
{
}


AliForwardFlowResultStorage::AliForwardFlowResultStorage(const TString &name, TList *list)
  : TNamed(name.Data(), "ForwardFlowResultStorage")
  , fObjects(new TObjArray())
{
  fObjects->Add(list);
}


AliForwardFlowResultStorage::~AliForwardFlowResultStorage()
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
AliForwardFlowResultStorage::Write(const char *name, Int_t option, Int_t bufsize)
{
  return const_cast<const AliForwardFlowResultStorage*>(this)->Write(name, option, bufsize);
}

Int_t
AliForwardFlowResultStorage::Write(const char *name, Int_t option, Int_t bufsize) const
{
  TString path(name);

  // override if 'path' is empty
  if (path.IsWhitespace()) {
    path = fName.Strip(TString::kBoth);
  }

  gDirectory->mkdir(path);
  TDirectory *outdir = gDirectory->GetDirectory(path);
  if (!outdir) {
    Error("AliForwardFlowResultStorage::Write", "Could not create path %s", path.Data());
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