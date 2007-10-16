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
//     AOD FMD cluster class
//     Author: Markus Oldenburg, CERN
//-------------------------------------------------------------------------

#include "AliAODFmdCluster.h"

ClassImp(AliAODFmdCluster)
  
//______________________________________________________________________________
AliAODFmdCluster::AliAODFmdCluster() : 
  AliAODCluster(),
  fProdVertex(NULL),
  fPrimTrack(NULL)
{
  // default constructor
}

//______________________________________________________________________________
AliAODFmdCluster::AliAODFmdCluster(Int_t id,
				   Int_t nLabel,
				   Int_t *label, 
				   Double_t energy,
				   Double_t x[3],
				   Double_t pid[9],
				   Char_t ttype,
				   AliAODVertex *prodVertex,
				   AliAODTrack *primTrack) :
  AliAODCluster(id, nLabel, label, energy, x, pid, ttype),
  fProdVertex(prodVertex),
  fPrimTrack(primTrack)
{
  // constructor
}

//______________________________________________________________________________
AliAODFmdCluster::AliAODFmdCluster(Int_t id,
				   Int_t nLabel,
				   Int_t *label, 
				   Float_t energy,
				   Float_t x[3],
				   Float_t pid[9],
				   Char_t ttype,
				   AliAODVertex *prodVertex,
				   AliAODTrack *primTrack) :
  AliAODCluster(id, nLabel, label, energy, x, pid, ttype),
  fProdVertex(prodVertex),
  fPrimTrack(primTrack)
{
  // constructor
}


//______________________________________________________________________________
AliAODFmdCluster::~AliAODFmdCluster() 
{
  // destructor
}


//______________________________________________________________________________
AliAODFmdCluster::AliAODFmdCluster(const AliAODFmdCluster& clus) :
  AliAODCluster(clus),
  fProdVertex(clus.fProdVertex),
  fPrimTrack(clus.fPrimTrack)
{
  // Copy constructor
}

//______________________________________________________________________________
AliAODFmdCluster& AliAODFmdCluster::operator=(const AliAODFmdCluster& clus)
{
  // Assignment operator
  if(this!=&clus) {

    AliAODCluster::operator=(clus);

    fProdVertex = clus.fProdVertex;
    fPrimTrack = clus.fPrimTrack;
  }

  return *this;
}
