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

//_________________________________________________________________________
//  Base class for the clusterization algorithm (pure abstract)
//*--
//*-- Author: Yves Schutz  SUBATECH 
// Modif: 
//  August 2002 Yves Schutz: clone PHOS as closely as possible and intoduction
//                           of new  IO (à la PHOS)
//////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---
#include "TClonesArray.h"
#include "TTree.h"

// --- Standard library ---


// --- AliRoot header files ---
#include "AliEMCALClusterizer.h"
#include "AliLog.h"

ClassImp(AliEMCALClusterizer)

//____________________________________________________________________________
AliEMCALClusterizer::AliEMCALClusterizer():
  fDigitsArr(NULL),
  fTreeR(NULL),
  fRecPoints(NULL)
{
  // ctor
}

//____________________________________________________________________________
AliEMCALClusterizer::~AliEMCALClusterizer()
{
  // dtor
  if (fDigitsArr) {
    fDigitsArr->Delete();
    delete fDigitsArr;
  }
  if (fRecPoints) {
    fRecPoints->Delete();
    delete fRecPoints;
  }
}

//____________________________________________________________________________
void AliEMCALClusterizer::SetInput(TTree *digitsTree)
{
  // Read the digits from the input tree
  TBranch *branch = digitsTree->GetBranch("EMCAL");
  if (!branch) { 
    AliError("can't get the branch with the EMCAL digits !");
    return;
  }
  fDigitsArr = new TClonesArray("AliEMCALDigit",100);
  branch->SetAddress(&fDigitsArr);
  branch->GetEntry(0);
}

//____________________________________________________________________________
void AliEMCALClusterizer::SetOutput(TTree *clustersTree)
{
  // Read the digits from the input tree
  fTreeR = clustersTree;
  
  AliDebug(9, "Making array for EMCAL clusters");
  fRecPoints = new TObjArray(100) ;
  Int_t split = 0;
  Int_t bufsize = 32000;
  fTreeR->Branch("EMCALECARP", "TObjArray", &fRecPoints, bufsize, split);
}
