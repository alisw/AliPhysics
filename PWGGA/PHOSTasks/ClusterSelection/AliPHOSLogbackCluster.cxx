/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
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

#include "AliPHOSLogbackCluster.h"

ClassImp(AliPHOSLogbackCluster)

//______________________________________________________________________________
AliPHOSLogbackCluster::AliPHOSLogbackCluster(AliVCluster* cluster) 
: fE( cluster->E() )
{
  //TODO: fCoreE

  cluster->GetPosition( fPosition );
}

AliPHOSLogbackCluster::~AliPHOSLogbackCluster()
{
}


TLorentzVector AliPHOSLogbackCluster::GetMomentum(Double_t * vertex)
{
   // Returns TLorentzVector with momentum of the cluster. Only valid for clusters 
  // identified as photons or pi0 (overlapped gamma) produced on the vertex
  // Vertex can be recovered with esd pointer doing:  
  //" Double_t vertex[3] ; esd->GetVertex()->GetXYZ(vertex) ; "

  Double_t pos[3]={ fPosition[0] - vertex[0], 
		    fPosition[1] - vertex[1], 
		    fPosition[2] - vertex[2]  };
  
  Double_t r = TMath::Sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]   ) ; 
    
  return TLorentzVector( fEnergy*pos[0]/r,  fEnergy*pos[1]/r,  fEnergy*pos[2]/r,  fEnergy) ;
}
