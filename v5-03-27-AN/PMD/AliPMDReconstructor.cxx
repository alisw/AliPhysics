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

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// class for PMD reconstruction                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "Riostream.h"
#include "AliPMDReconstructor.h"
#include "AliRun.h"
#include "AliPMDClusterFinder.h"
#include "AliPMDtracker.h"
#include "AliRawReader.h"
#include "AliESDPmdTrack.h"
#include "AliESDEvent.h"
#include "AliLog.h"
#include "AliRunInfo.h"

ClassImp(AliPMDReconstructor)

// ------------------------------------------------------------------------ //

void AliPMDReconstructor::Reconstruct(AliRawReader *rawReader,
				      TTree *clustersTree) const
{
// reconstruct clusters from Raw Data
// Event Species Added By satyajit
  Int_t gRecoMode = 1;
  TString beamType = GetRunInfo()->GetBeamType();
  
  if (((beamType.CompareTo("pp"))==0) || 
      ((beamType.CompareTo("p-p"))==0)||
      ((beamType.CompareTo("PP"))==0) || 
      ((beamType.CompareTo("P-P"))==0)) {
    gRecoMode=1;
  }
  
  else if ((beamType.CompareTo("A-A")) == 0 || 
	   (beamType.CompareTo("AA")) == 0) {
    gRecoMode=2;
  }

    static AliPMDClusterFinder pmdClus;
    pmdClus.Digits2RecPoints(rawReader, clustersTree, gRecoMode);
}

// ------------------------------------------------------------------------ //
void AliPMDReconstructor::Reconstruct(TTree *digitsTree,
				      TTree *clustersTree) const
{
  // reconstruct clusters from Digits
  // Setting reconstruction mode
 
  // Added to Have Sepatrate Event Spcies
  Int_t gRecoMode = 1;
  TString beamType = GetRunInfo()->GetBeamType();
 
  if (((beamType.CompareTo("pp"))==0) || 
     ((beamType.CompareTo("p-p"))==0)||
     ((beamType.CompareTo("PP"))==0) || 
     ((beamType.CompareTo("P-P"))==0)) {
    gRecoMode=1;
  }
  
  else if ((beamType.CompareTo("A-A")) == 0 || 
	   (beamType.CompareTo("AA")) == 0) {
    gRecoMode=2;
  }
  
  static AliPMDClusterFinder pmdClus;
  pmdClus.Digits2RecPoints(digitsTree, clustersTree,gRecoMode);

}

// ------------------------------------------------------------------------ //
void AliPMDReconstructor::FillESD(AliRawReader* /*rawReader*/,
				  TTree* clustersTree, AliESDEvent* esd) const
{
    static AliPMDtracker pmdtracker;
    pmdtracker.LoadClusters(clustersTree);
    pmdtracker.Clusters2Tracks(esd);
}
// ------------------------------------------------------------------------ //
void AliPMDReconstructor::FillESD(TTree * /*digitsTree*/,
				  TTree* clustersTree, AliESDEvent* esd) const
{
    static AliPMDtracker pmdtracker;
    pmdtracker.LoadClusters(clustersTree);
    pmdtracker.Clusters2Tracks(esd);
}


