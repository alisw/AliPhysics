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

//_________________________________________________________________________
//  Base class for the clusterization algorithm (pure abstract)
//*--
//*-- Author: Yves Schutz  SUBATECH 
//////////////////////////////////////////////////////////////////////////////

#include <TClonesArray.h>
#include <TTree.h>

#include "AliPHOSClusterizer.h"
#include "AliPHOSDigit.h"
#include "AliLog.h"

ClassImp(AliPHOSClusterizer)

//____________________________________________________________________________
AliPHOSClusterizer::AliPHOSClusterizer():
  fGeom(NULL),
  fDigitsArr(0),
  fTreeR(0),
  fEMCRecPoints(0),
  fCPVRecPoints(0)
{
  // ctor
}

//____________________________________________________________________________
AliPHOSClusterizer::AliPHOSClusterizer(AliPHOSGeometry *geom):
  fGeom(geom),
  fDigitsArr(0),
  fTreeR(0),
  fEMCRecPoints(0),
  fCPVRecPoints(0)
{
  // ctor
 
}

//____________________________________________________________________________
AliPHOSClusterizer::~AliPHOSClusterizer()
{
  // dtor
  if (fDigitsArr) {
    fDigitsArr->Delete();
    delete fDigitsArr;
  }
  if (fEMCRecPoints) {
    fEMCRecPoints->Delete();
    delete fEMCRecPoints;
  }
  if (fCPVRecPoints) {
    fCPVRecPoints->Delete();
    delete fCPVRecPoints;
  }
}

//____________________________________________________________________________
void AliPHOSClusterizer::SetInput(TTree * digitsTree) 
{
  // Get the tree with digits and sets
  // the input array with digits for PHOS
  TBranch *branch = digitsTree->GetBranch("PHOS");
  if (!branch) { 
    AliError("can't get the branch with the PHOS digits !");
    return;
  }
  fDigitsArr = new TClonesArray("AliPHOSDigit",100);
  branch->SetAddress(&fDigitsArr);
  branch->GetEntry(0);
}

//____________________________________________________________________________
void AliPHOSClusterizer::SetOutput(TTree * clustersTree) 
{
  // Set the output clusters tree,
  // creates the arrays for EMC and CPV,
  // and set the corresponding branch addresses
  fTreeR = clustersTree;

  AliDebug(9, "Making array for EMC clusters");
  fEMCRecPoints = new TObjArray(100) ;
  fEMCRecPoints->SetName("EMCRECPOINTS") ;
  Int_t split = 0;
  Int_t bufsize = 32000;
  fTreeR->Branch("PHOSEmcRP", "TObjArray", &fEMCRecPoints, bufsize, split);

  AliDebug(9, "Making array for CPV clusters");
  fCPVRecPoints = new TObjArray(100) ;
  fCPVRecPoints->SetName("CPVRECPOINTS") ;
  fTreeR->Branch("PHOSCpvRP", "TObjArray", &fCPVRecPoints, bufsize, split);
}
