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


#include "AliROCRawAnalysisSelector.h"

#include <AliLog.h>

#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>


ClassImp(AliROCRawAnalysisSelector)

AliROCRawAnalysisSelector::AliROCRawAnalysisSelector() :
  AliSelector()
{
  //
  // Constructor. Initialization of pointers
  //
  
  //  for (Int_t i=0; i<kTPCSectors; i++)
  // fClusterHistograms[i] = 0;
}

AliROCRawAnalysisSelector::~AliROCRawAnalysisSelector()
{
  //
  // Destructor
  //
}

void AliROCRawAnalysisSelector::SlaveBegin(TTree* tree)
{
  //
  
  AliSelector::SlaveBegin(tree);
} 

void AliROCRawAnalysisSelector::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses of the tree
  // will be set. It is normaly not necessary to make changes to the
  // generated code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running with PROOF.

  AliSelector::Init(tree);

  // Set branch address
  //if (tree)
  //  tree->SetBranchAddress("rawevent", &fRawEvent);
    
}

Bool_t AliROCRawAnalysisSelector::Process(Long64_t entry)
{
  //
  // Implement your analysis here. Do not forget to call the parent class Process by
  // if (AliSelector::Process(entry) == kFALSE)
  //   return kFALSE;
  //

  if (AliSelector::Process(entry) == kFALSE)
    return kFALSE;

  
   
  return kTRUE;
}

void AliROCRawAnalysisSelector::SlaveTerminate()
{
  //
  
  //for (Int_t i=0; i<kTPCSectors; i++)
  // if (fClusterHistograms[i])
  //    fOutput->Add(fClusterHistograms[i]);
} 

void AliROCRawAnalysisSelector::Terminate()
{
  // TODO read from output list for PROOF
    
  TFile* file = TFile::Open("rocRaw.root", "RECREATE");
  
  //  for (Int_t i=0; i<kTPCSectors; i++)
  //  if (fClusterHistograms[i])
  //    fClusterHistograms[i]->SaveHistograms();

  file->Close();
} 
