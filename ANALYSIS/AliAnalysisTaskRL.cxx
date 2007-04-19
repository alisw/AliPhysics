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

//-----------------------------------------------------------------
//           AliAnalysisTaskRL class
//     Task that gives access to the run loader
//   Origin: Panos Christakoglou, UOA-CERN, Panos.Christakoglou@cern.ch
//-----------------------------------------------------------------

#include "AliAnalysisTaskRL.h"
#include "AliRunLoader.h"
#include "AliStack.h"
#include "AliHeader.h"


#include <TTree.h>
#include <TChain.h>
#include <TFile.h>

ClassImp(AliAnalysisTaskRL)

//___________________________________________________________________________
AliAnalysisTaskRL::AliAnalysisTaskRL() :
  AliAnalysisTask(),
  fTree(0), fRunLoader(0),
  fKinematicsLoaded(kFALSE), fHeaderLoaded(kFALSE), fTreeNumber(-1) {
  //
  // Constructor. Initialization of pointers
  //
}

//___________________________________________________________________________
AliAnalysisTaskRL::AliAnalysisTaskRL(const char *name, const char *title) :
  AliAnalysisTask(name,title),
  fTree(0), fRunLoader(0),
  fKinematicsLoaded(kFALSE), fHeaderLoaded(kFALSE) {
  // Constructor.
}

//___________________________________________________________________________
AliAnalysisTaskRL::~AliAnalysisTaskRL() {
  //
  // Destructor
  //
  DeleteRunLoader();
}

//___________________________________________________________________________
Bool_t AliAnalysisTaskRL::GetEntry(Long64_t ientry) {
  //returns the entry of the run loader
  if(fRunLoader) {
    if(fRunLoader->GetEvent((Int_t)ientry) != 0)
      return kFALSE;
  }
  return kTRUE;
}

//___________________________________________________________________________
AliRunLoader *AliAnalysisTaskRL::GetRunLoader() {
  // Returns AliRun instance corresponding to current ESD active in fTree
  // Loads galice.root, the file is identified by replacing "AliESDs" to
  // "galice" in the file path of the ESD file. 

  fTree = (TTree *)AliAnalysisTask::GetInputData(0);
  Int_t iTree = ((TChain *)AliAnalysisTask::GetInputData(0))->GetTreeNumber();
  if (iTree != fTreeNumber) {
      DeleteRunLoader();
      fTreeNumber = iTree;
  }
  
      
  if (!fRunLoader) {
    if (!fTree->GetCurrentFile())
      return 0;
    
    TString fileName(fTree->GetCurrentFile()->GetName());
    printf("Current file %s \n", fileName.Data());
    
    fileName.ReplaceAll("AliESDs", "galice");
    
    // temporary workaround for PROOF bug #18505
    fileName.ReplaceAll("#galice.root#galice.root", "#galice.root");
    
    fRunLoader = AliRunLoader::Open(fileName);
    if (!fRunLoader)
      return 0;
    fRunLoader->GetEvent((Int_t)(fTree->GetTree()->GetReadEntry()));
  }
  
  return fRunLoader;
}

//___________________________________________________________________________
void AliAnalysisTaskRL::DeleteRunLoader() {
  //
  // deletes the runloader
  //
  if (fRunLoader) {
    fRunLoader->Delete();
    fRunLoader = 0;
  }
  
  fKinematicsLoaded = kFALSE;
  fHeaderLoaded = kFALSE;
}

//___________________________________________________________________________
AliHeader* AliAnalysisTaskRL::GetHeader() {
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

//___________________________________________________________________________
AliStack* AliAnalysisTaskRL::GetStack() {
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
