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

//-------------------------------------------------------------------------
//     AOD PMD cluster class
//     Author: Markus Oldenburg, CERN
//-------------------------------------------------------------------------

#include "AliAODPmdCluster.h"

ClassImp(AliAODPmdCluster)

//______________________________________________________________________________
AliAODPmdCluster::AliAODPmdCluster() : 
  AliAODCluster(),
  fAssocCluster(NULL)
{
  // default constructor
}

//______________________________________________________________________________
AliAODPmdCluster::AliAODPmdCluster(Int_t id,
				   Int_t nLabel,
				   Int_t *label, 
				   Double_t energy,
				   Double_t x[3],
				   Double_t pid[13],
				   Char_t ttype,
				   UInt_t /*selectInfo*/,
				   AliAODPmdCluster* assoc) :
  AliAODCluster(id, nLabel, label, energy, x, pid, ttype),
  fAssocCluster(assoc)
{
  // constructor
}

//______________________________________________________________________________
AliAODPmdCluster::AliAODPmdCluster(Int_t id,
				   Int_t nLabel,
				   Int_t *label, 
				   Float_t energy,
				   Float_t x[3],
				   Float_t pid[13],
				   Char_t ttype,
				   UInt_t /*selectInfo*/,
				   AliAODPmdCluster* assoc) :
  AliAODCluster(id, nLabel, label, energy, x, pid, ttype),
  fAssocCluster(assoc)
{
  // constructor
}


//______________________________________________________________________________
AliAODPmdCluster::~AliAODPmdCluster() 
{
  // destructor
}


//______________________________________________________________________________
AliAODPmdCluster::AliAODPmdCluster(const AliAODPmdCluster& clus) :
  AliAODCluster(clus),
  fAssocCluster(clus.fAssocCluster)
{
  // Copy constructor
}

//______________________________________________________________________________
AliAODPmdCluster& AliAODPmdCluster::operator=(const AliAODPmdCluster& clus)
{
  // Assignment operator
  if(this!=&clus) {

    AliAODCluster::operator=(clus);
    
    fAssocCluster = clus.fAssocCluster;
  }

  return *this;
}
