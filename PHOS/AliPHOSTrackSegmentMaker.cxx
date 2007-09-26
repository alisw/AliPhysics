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

/* History of cvs commits:
 *
 * $Log$
 * Revision 1.29  2007/08/28 12:55:08  policheh
 * Loaders removed from the reconstruction code (C.Cheshkov)
 *
 * Revision 1.28  2007/08/07 14:12:03  kharlov
 * Quality assurance added (Yves Schutz)
 *
 * Revision 1.27  2006/08/25 16:56:30  kharlov
 * Compliance with Effective C++
 *
 * Revision 1.26  2006/08/25 16:00:53  kharlov
 * Compliance with Effective C++AliPHOSHit.cxx
 *
 * Revision 1.25  2005/05/28 14:19:05  schutz
 * Compilation warnings fixed by T.P.
 *
 */

//_________________________________________________________________________
// Algorithm Base class to construct PHOS track segments
// Associates EMC and PPSD clusters
// Unfolds the EMC cluster   
//*-- 
//*-- Author: Dmitri Peressounko (RRC Ki & SUBATECH)


// --- ROOT system ---
#include "TTree.h"

// --- Standard library ---

// --- AliRoot header files ---
#include "AliPHOSTrackSegmentMaker.h"
#include "AliLog.h"

ClassImp( AliPHOSTrackSegmentMaker) 


//____________________________________________________________________________
AliPHOSTrackSegmentMaker:: AliPHOSTrackSegmentMaker() : 
  TObject(),
  fESD(0), 
  fGeom(0),
  fEMCRecPoints(0),
  fCPVRecPoints(0)
{
 // ctor
}

//____________________________________________________________________________
AliPHOSTrackSegmentMaker::AliPHOSTrackSegmentMaker(AliPHOSGeometry *geom):
  TObject(),
  fESD(0), 
  fGeom(geom),
  fEMCRecPoints(0),
  fCPVRecPoints(0)
{
  // ctor
}

//____________________________________________________________________________
AliPHOSTrackSegmentMaker::AliPHOSTrackSegmentMaker(const AliPHOSTrackSegmentMaker & tsmaker) :
  TObject(tsmaker),
  fESD(tsmaker.GetESD()), 
  fGeom(tsmaker.fGeom),
  fEMCRecPoints(tsmaker.fEMCRecPoints),
  fCPVRecPoints(tsmaker.fCPVRecPoints)
{
  //Copy constructor
} 

//____________________________________________________________________________
AliPHOSTrackSegmentMaker::~AliPHOSTrackSegmentMaker()
{
 //Remove this from the parental task before destroying
  //  if(AliPHOSGetter::Instance()->PhosLoader())
  //    AliPHOSGetter::Instance()->PhosLoader()->CleanTracker();
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
void AliPHOSTrackSegmentMaker::SetInput(TTree *clustersTree)
{
  // Read the clusters tree and creates the
  // arrays with the EMC and CPV
  // clusters.
  // and set the corresponding branch addresses

  TBranch *emcbranch = clustersTree->GetBranch("PHOSEmcRP");
  if (!emcbranch) { 
    AliError("can't get the branch with the PHOS EMC clusters !");
    return;
  }
  fEMCRecPoints = new TObjArray(100) ;
  emcbranch->SetAddress(&fEMCRecPoints);
  emcbranch->GetEntry(0);

  TBranch *cpvbranch = clustersTree->GetBranch("PHOSCpvRP");
  if (!cpvbranch) { 
    AliError("can't get the branch with the PHOS CPV clusters !");
    return;
  }
  fCPVRecPoints = new TObjArray(100) ;
  cpvbranch->SetAddress(&fCPVRecPoints);
  cpvbranch->GetEntry(0);
}
