#include "AliSelectorRL.h"

#include <AliLog.h>
#include <AliRun.h>
#include <AliRunLoader.h>

ClassImp(AliSelectorRL)

AliSelectorRL::AliSelectorRL() :
  AliSelector(),
  fRunLoader(0)
{
  //
  // Constructor. Initialization of pointers
  //
}

AliSelectorRL::~AliSelectorRL()
{
  //
  // Destructor
  //

  // histograms are in the output list and deleted when the output
  // list is deleted by the TSelector dtor
}

Bool_t AliSelectorRL::Notify()
{
  // Calls base class Notify
  // On top of that run loader is closed, because we change the input file

  if (AliSelector::Notify() == kFALSE)
    return kFALSE;

  DeleteRunLoader();

  return kTRUE;
}

void AliSelectorRL::SlaveTerminate()
{
  // removes runloader

  AliSelector::SlaveTerminate();

  DeleteRunLoader();
}

AliRun* AliSelectorRL::GetAliRun()
{
  // Returns AliRun instance corresponding to current ESD active in fChain
  // Loads galice.root, the file is identified by replacing "AliESDs" to
  // "galice" in the file path of the ESD file. This is a hack, to be changed!

  if (!fRunLoader)
  {
    if (!fChain->GetCurrentFile())
      return 0;

    TString fileName(fChain->GetCurrentFile()->GetName());
    fileName.ReplaceAll("AliESDs", "galice");

    fRunLoader = AliRunLoader::Open(fileName);
    if (!fRunLoader)
      return 0;

    fRunLoader->LoadgAlice();
  }

  return fRunLoader->GetAliRun();
}

void AliSelectorRL::DeleteRunLoader()
{
  //
  // deletes the runloader
  //

  if (fRunLoader)
  {
    fRunLoader->Delete();
    fRunLoader = 0;
  }
}
