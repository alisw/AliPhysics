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

#include <Riostream.h>
#include <TMath.h>
#include <AliITSlayerUpgrade.h>
#include <AliITStrackerMI.h>

AliITSlayerUpgrade::AliITSlayerUpgrade():
  fPhiOffset(0.),
  fZOffset(0.),
  fNsel(0),
  fN(0)
{  //--------------------------------------------------------------------
  //default AliITSlayerUpgrade constructor
  //--------------------------------------------------------------------
  for (Int_t i=0; i<AliITSRecoParam::kMaxClusterPerLayer; i++) {
  fIndex[i] = -1;
  fClusters[i]= 0x0;
  }
}
//___________________________________________________________________________
AliITSlayerUpgrade::AliITSlayerUpgrade(Double_t p,Double_t z):
  fPhiOffset(p),
  fZOffset(z),
  fNsel(0),
  fN(0)
{
  //--------------------------------------------------------------------
  //main AliITSlayerUpgrade constructor
  //--------------------------------------------------------------------

  for (Int_t i=0; i<AliITSRecoParam::kMaxClusterPerLayer; i++) {
  fIndex[i] = -1;
  fClusters[i]= 0x0;
  }

}

AliITSlayerUpgrade::~AliITSlayerUpgrade() {
  //--------------------------------------------------------------------
  // AliITSlayerUpgrade destructor
  //--------------------------------------------------------------------
ResetClusters();

}
void AliITSlayerUpgrade::ResetClusters() {
  //--------------------------------------------------------------------
  // This function removes loaded clusters
  //--------------------------------------------------------------------

  for (Int_t i=0; i<fN; i++) {
  fIndex[i]=-1;
  if(fClusters[i]) {
    delete fClusters[i];
    }
  }
  fN=0;
  return;
}


Int_t AliITSlayerUpgrade::InsertCluster(AliITSRecPoint *c) {
//--------------------------------------------------------------------
  // This function inserts a cluster to this layer in increasing
  // order of the cluster's fZ
  //--------------------------------------------------------------------
  fClusters[fN]=c;
  fN++;
  return 0;
}

//-----------------------------------------------------------------------
const AliITSRecPoint *AliITSlayerUpgrade::GetNextCluster(Int_t &ci){
  //--------------------------------------------------------------------
  // This function returns clusters within the "window"
  //--------------------------------------------------------------------
  AliITSRecPoint *c=0;
  ci=-1;
  if (fNsel) {
     fNsel--;
     ci=fIndex[fNsel];
     c=fClusters[ci];
  }
  return c;
}

