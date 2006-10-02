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

/* $Id$ */

#include "AliSelectorRL.h"

#include <AliLog.h>
#include <AliRunLoader.h>

#include <TTree.h>
#include <TFile.h>

// Selector base class for analysis based on whole AliRoot
// Please derive your selector-based analysis from this class, if you need to access the
// RunLoader or anything provided by it.
//
// This selector depends on the RunLoader, therefore to use it you have to have the whole AliRoot available
// The benefit is that you can use the RunLoader to access everything in the data structure
//
// The class has functions to return the stack (GetStack()), the RunLoader (GetRunLoader()), and
// the event header (GetHeader()).
//
// Author: Jan.Fiete.Grosse-Oetringhaus@cern.ch

ClassImp(AliSelectorRL)

AliSelectorRL::AliSelectorRL() :
  AliSelector(),
  fRunLoader(0),
  fKinematicsLoaded(kFALSE),
  fHeaderLoaded(kFALSE)
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
  {
    if (fRunLoader->GetEvent(entry) != 0)
      return kFALSE;
  }

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
  // "galice" in the file path of the ESD file. 

  if (!fRunLoader)
  {
    if (!fTree->GetCurrentFile())
      return 0;

    TString fileName(fTree->GetCurrentFile()->GetName());
    fileName.ReplaceAll("AliESDs", "galice");

    // temporary workaround for PROOF bug #18505
    fileName.ReplaceAll("#galice.root#galice.root", "#galice.root");

    fRunLoader = AliRunLoader::Open(fileName);
    if (!fRunLoader)
      return 0;

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

  fKinematicsLoaded = kFALSE;
  fHeaderLoaded = kFALSE;
}

AliHeader* AliSelectorRL::GetHeader()
{
  // Returns header retrieved from RunLoader

  AliRunLoader* runLoader = GetRunLoader();
  if (!runLoader)
    return 0;

  if (fHeaderLoaded == kFALSE)
    if (runLoader->LoadHeader() != 0)
      return 0;

  fHeaderLoaded = kTRUE;

  return runLoader->GetHeader();
}

AliStack* AliSelectorRL::GetStack()
{
  // Returns stack retrieved from RunLoader

  AliRunLoader* runLoader = GetRunLoader();
  if (!runLoader)
    return 0;

  if (fKinematicsLoaded == kFALSE)
    if (runLoader->LoadKinematics() != 0)
      return 0;

  fKinematicsLoaded = kTRUE;

  return runLoader->Stack();
}
