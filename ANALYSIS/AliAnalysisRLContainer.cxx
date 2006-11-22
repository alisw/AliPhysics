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
// Author: Andrei Gheata, 31/05/2006

//==============================================================================
//   AliAnalysysRLContainer - 
//
//
//
//==============================================================================

#include "AliAnalysisRLContainer.h"
#include "TTree.h"
#include "TFile.h"

#include "AliLog.h"
#include "AliRunLoader.h"
#include "AliESD.h"
#include "AliHeader.h"
#include "AliStack.h"
#include "AliAnalysisDataContainer.h"
#include "AliAnalysisDataSlot.h"
#include "AliAnalysisTask.h"

ClassImp(AliAnalysisRLContainer)

//______________________________________________________________________________
AliAnalysisRLContainer::AliAnalysisRLContainer()
                       :AliAnalysisDataContainer(),
                        fRunLoader(NULL),
                        fESD(NULL),
                        fKineFile(NULL),
                        fKinematicsLoaded(kFALSE),
                        fHeaderLoaded(kFALSE)
{
// Dummy ctor.
}

//______________________________________________________________________________
AliAnalysisRLContainer::AliAnalysisRLContainer(const char *name)
                       :AliAnalysisDataContainer(name, TTree::Class()),
                        fRunLoader(NULL),
                        fESD(NULL),
                        fKineFile(NULL),
                        fKinematicsLoaded(kFALSE),
                        fHeaderLoaded(kFALSE)
{
// Normal constructor.
}

//______________________________________________________________________________
AliAnalysisRLContainer::AliAnalysisRLContainer(const AliAnalysisRLContainer &rlc)
                       :AliAnalysisDataContainer(rlc),
                        fRunLoader(rlc.fRunLoader),
                        fESD(rlc.fESD),
                        fKineFile(rlc.fKineFile),
                        fKinematicsLoaded(rlc.fKinematicsLoaded),
                        fHeaderLoaded(rlc.fHeaderLoaded)
{
// Copy constructor.
}   

//______________________________________________________________________________
AliAnalysisRLContainer::~AliAnalysisRLContainer()
{
// Destructor. Deletes data ! (What happens if data is a container ???)
}

//______________________________________________________________________________
AliAnalysisRLContainer &AliAnalysisRLContainer::operator=(const AliAnalysisRLContainer &rlc)
{
// Assignment.
   if (&rlc != this) {
      AliAnalysisDataContainer::operator=(rlc);
      fRunLoader = rlc.fRunLoader;
      fESD = rlc.fESD;
      fKineFile = rlc.fKineFile;
      fKinematicsLoaded = rlc.fKinematicsLoaded;
      fHeaderLoaded = rlc.fHeaderLoaded;
   }
   return *this;   
}

//______________________________________________________________________________
Bool_t AliAnalysisRLContainer::SetData(TObject *data, Option_t */*option*/)
{
// Data must be a tree here.
   fData = data;
   TTree *tree = (TTree *)data;
   // Set branch address
   tree->SetBranchAddress("ESD", &fESD);
   fDataReady = kTRUE;
   return kTRUE;
}

//______________________________________________________________________________
void AliAnalysisRLContainer::GetEntry(Long64_t ientry)
{
// If data is ready and derives from TTree or from TBranch, this will get the
// requested entry in memory if not already loaded.
   if (!fDataReady) return;
   TTree *tree = (TTree*)fData;
   tree->GetTree()->GetEntry(ientry);
   if (fRunLoader) fRunLoader->GetEvent(ientry);
}   

//______________________________________________________________________________
void AliAnalysisRLContainer::NotifyChange(ENotifyMessage type)
{
// Notify container that file has changed.
   AliAnalysisDataContainer::NotifyChange(type);
   if (!type==kFileChange) return;
   DeleteKinematicsFile();
   DeleteRunLoader();
}

//______________________________________________________________________________
void AliAnalysisRLContainer::DeleteRunLoader() 
{
// Deletes the runloader.
   if (fRunLoader) {
      fRunLoader->Delete();
      fRunLoader = 0;
   }
   fKinematicsLoaded = kFALSE;
   fHeaderLoaded = kFALSE;
}

//______________________________________________________________________________
void AliAnalysisRLContainer::DeleteKinematicsFile() 
{
// Closes the kinematics file and deletes the pointer.
   if (fKineFile) {
      fKineFile->Close();
      delete fKineFile;
      fKineFile = 0;
   }
}   

//______________________________________________________________________________
AliRunLoader* AliAnalysisRLContainer::GetRunLoader()
{
// Returns AliRun instance corresponding to current ESD active in fTree
// Loads galice.root, the file is identified by replacing "AliESDs" to
// "galice" in the file path of the ESD file. This is a hack, to be changed!

   if (!fDataReady) return 0;
   TTree *tree = (TTree*)fData;
   if (!fRunLoader) {
      if (!tree->GetCurrentFile()) return 0;
      TString fileName(tree->GetCurrentFile()->GetName());
      fileName.ReplaceAll("AliESDs", "galice");
      fRunLoader = AliRunLoader::Open(fileName);
      if (!fRunLoader) return 0;
      if (fRunLoader->LoadgAlice() != 0) {
         delete fRunLoader;
         fRunLoader = 0;
         return 0;
      }
      fRunLoader->GetEvent(tree->GetTree()->GetReadEntry());
  }
  return fRunLoader;
}

//______________________________________________________________________________
AliHeader* AliAnalysisRLContainer::GetHeader()
{
// Returns header retrieved from RunLoader
   AliRunLoader* runLoader = GetRunLoader();
   if (!runLoader) return 0;
   if (!fHeaderLoaded) 
      if (runLoader->LoadHeader() != 0) return 0;
   fHeaderLoaded = kTRUE;
   return runLoader->GetHeader();
}

//______________________________________________________________________________
TTree* AliAnalysisRLContainer::GetKinematics()
{
// Returns kinematics tree corresponding to current ESD active in fTree
// Loads the kinematics from the kinematics file, the file is identified by replacing "AliESDs" to
// "Kinematics" in the file path of the ESD file. This is a hack, to be changed!

   if (!fDataReady) return 0;
   TTree *tree = (TTree*)fData;
   if (!fKineFile) {
      if (!tree->GetCurrentFile()) return 0;
      TString fileName(tree->GetCurrentFile()->GetName());
      fileName.ReplaceAll("AliESDs", "Kinematics");

      AliDebug(AliLog::kInfo, Form("Opening %s", fileName.Data()));

      fKineFile = TFile::Open(fileName);
      if (!fKineFile) return 0;
   }
   return dynamic_cast<TTree*> (fKineFile->Get(Form("Event%d/TreeK", tree->GetTree()->GetReadEntry())));
}

//______________________________________________________________________________
AliStack* AliAnalysisRLContainer::GetStack()
{
// Returns stack retrieved from RunLoader

   AliRunLoader* runLoader = GetRunLoader();
   if (!runLoader) return 0;
   if (!fKinematicsLoaded)
      if (runLoader->LoadKinematics() != 0) return 0;
   fKinematicsLoaded = kTRUE;
   return runLoader->Stack();
}
