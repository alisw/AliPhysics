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
//  Algorithm class for the identification of particles detected in PHOS        
//  base  class  of identified particle  
//  Why should I put meaningless comments
//  just to satisfy
//  the code checker                
                         
//                  
//*-- Author: Yves Schutz (SUBATECH) & Dmitri Peressounko


// --- ROOT system ---
#include "TBranch.h"
#include "TClonesArray.h"
#include "TTree.h"

// --- Standard library ---

// --- AliRoot header files ---
#include "AliConfig.h"
#include "AliLog.h"
#include "AliPHOSPID.h"

ClassImp(AliPHOSPID)

//____________________________________________________________________________
AliPHOSPID::AliPHOSPID():
  TObject(),
  fGeom(NULL),
  fESD(0x0),
  fEMCRecPoints(NULL),
  fCPVRecPoints(NULL),
  fTrackSegments(NULL),
  fRecParticles(NULL),
  fEnergyCorrectionOn(kTRUE)
{
  // ctor
}


//____________________________________________________________________________
AliPHOSPID::AliPHOSPID(AliPHOSGeometry *geom):
  TObject(),
  fGeom(geom),
  fESD(0x0),
  fEMCRecPoints(NULL),
  fCPVRecPoints(NULL),
  fTrackSegments(NULL),
  fRecParticles(NULL),
  fEnergyCorrectionOn(kTRUE)
{
  // ctor
  fEMCRecPoints = new TObjArray(100) ;
  fCPVRecPoints = new TObjArray(100) ;
  fRecParticles = new TClonesArray("AliPHOSRecParticle",100) ;
  fRecParticles->SetName("RECPARTICLES");

}

//____________________________________________________________________________
AliPHOSPID::AliPHOSPID(const AliPHOSPID & pid) :
  TObject(pid),
  fGeom(pid.fGeom),
  fESD(pid.fESD), 
  fEMCRecPoints(pid.fEMCRecPoints),
  fCPVRecPoints(pid.fCPVRecPoints),
  fTrackSegments(pid.fTrackSegments),
  fRecParticles(pid.fRecParticles),
  fEnergyCorrectionOn(pid.fEnergyCorrectionOn)
{
  // Copy constructor
}

//____________________________________________________________________________
AliPHOSPID::~AliPHOSPID()
{
  // dtor
  if (fEMCRecPoints) {
    fEMCRecPoints->Delete();
    delete fEMCRecPoints;
  }
  if (fCPVRecPoints) {
    fCPVRecPoints->Delete();
    delete fCPVRecPoints;
  }
  if (fRecParticles) {
    fRecParticles->Delete();
    delete fRecParticles;
  }
}

//____________________________________________________________________________
void AliPHOSPID::SetInput(TTree *clustersTree, TClonesArray *trackSegments)
{
  // Read the clusters tree and creates the
  // arrays with the EMC and CPV
  // clusters.
  // and set the corresponding branch addresses

  fTrackSegments = trackSegments;

  TBranch *emcbranch = clustersTree->GetBranch("PHOSEmcRP");
  if (!emcbranch) { 
    AliError("can't get the branch with the PHOS EMC clusters !");
    return;
  }
  emcbranch->SetAddress(&fEMCRecPoints);
  fEMCRecPoints->Delete();
  emcbranch->GetEntry(0);

  TBranch *cpvbranch = clustersTree->GetBranch("PHOSCpvRP");
  if (!cpvbranch) { 
    AliError("can't get the branch with the PHOS CPV clusters !");
    return;
  }
  cpvbranch->SetAddress(&fCPVRecPoints);
  fCPVRecPoints->Delete();
  cpvbranch->GetEntry(0);
}
