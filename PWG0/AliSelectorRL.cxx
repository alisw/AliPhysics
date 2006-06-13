/* $Id$ */

#include "AliSelectorRL.h"

#include <AliLog.h>
#include <AliRunLoader.h>

#include <TTree.h>
#include <TFile.h>

//
// This selector depends on the RunLoader, therefore to use it you have to have the whole AliRoot available
// The benefit is that you can use the RunLoader to access everything in the data structure
// If you only have the ESD library use AliSelector instead
//

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

Bool_t AliSelectorRL::Process(Long64_t entry)
{
  // Call the baseclass Process and set event number in runLoader (if loaded)

  if (AliSelector::Process(entry) == kFALSE)
    return kFALSE;

  if (fRunLoader)
    fRunLoader->GetEvent(entry);

  return kTRUE;
}

void AliSelectorRL::SlaveTerminate()
{
  // removes runloader

  AliSelector::SlaveTerminate();

  DeleteRunLoader();
}

AliRunLoader* AliSelectorRL::GetRunLoader()
{
  // Returns AliRun instance corresponding to current ESD active in fTree
  // Loads galice.root, the file is identified by replacing "AliESDs" to
  // "galice" in the file path of the ESD file. This is a hack, to be changed!

  if (!fRunLoader)
  {
    if (!fTree->GetCurrentFile())
      return 0;

    TString fileName(fTree->GetCurrentFile()->GetName());
    fileName.ReplaceAll("AliESDs", "galice");

    fRunLoader = AliRunLoader::Open(fileName);
    if (!fRunLoader)
      return 0;

    fRunLoader->LoadgAlice();
    fRunLoader->GetEvent(fTree->GetTree()->GetReadEntry());
  }

  return fRunLoader;
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

AliHeader* AliSelectorRL::GetHeader()
{
  // Returns header retrieved from RunLoader

  AliRunLoader* runLoader = GetRunLoader();
  if (!runLoader)
    return 0;

  if (runLoader->GetHeader() == 0)
    runLoader->LoadHeader();

  return runLoader->GetHeader();
}

AliStack* AliSelectorRL::GetStack()
{
  // Returns stack retrieved from RunLoader

  AliRunLoader* runLoader = GetRunLoader();
  if (!runLoader)
    return 0;

  if (runLoader->Stack() == 0)
    runLoader->LoadKinematics();

  return runLoader->Stack();
}
