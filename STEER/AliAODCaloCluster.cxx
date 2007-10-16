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
//     AOD calorimeter cluster class (for PHOS and EMCAL)
//     Author: Markus Oldenburg, CERN
//-------------------------------------------------------------------------

#include "AliAODCaloCluster.h"

ClassImp(AliAODCaloCluster)

//______________________________________________________________________________
AliAODCaloCluster::AliAODCaloCluster() : 
  AliAODCluster(),
  fDistToBadChannel(-999.),
  fDispersion(-1),
  fM20(0.),
  fM01(0.),
  fM02(0.),
  fM11(0.),
  fEmcCpvDistance(-999.),
  fNExMax(0),
  fTracksMatched(),
  fCellNumber()
{
  // default constructor

}

//______________________________________________________________________________
AliAODCaloCluster::AliAODCaloCluster(Int_t id,
				     UInt_t nLabel,
				     Int_t *label, 
				     Double_t energy,
				     Double_t x[3],
				     Double_t pid[9],
				     Char_t ttype,
				     UInt_t selectInfo) :
  AliAODCluster(id, nLabel, label, energy, x, pid, ttype, selectInfo),
  fDistToBadChannel(-999.),
  fDispersion(-1),
  fM20(0.),
  fM01(0.),
  fM02(0.),
  fM11(0.),
  fEmcCpvDistance(-999.),
  fNExMax(0),
  fTracksMatched(),
  fCellNumber()
{
  // constructor

}

//______________________________________________________________________________
AliAODCaloCluster::AliAODCaloCluster(Int_t id,
				     UInt_t nLabel,
				     Int_t *label, 
				     Float_t energy,
				     Float_t x[3],
				     Float_t pid[9],
				     Char_t ttype,
				     UInt_t selectInfo) :
  AliAODCluster(id, nLabel, label, energy, x, pid, ttype, selectInfo),
  fDistToBadChannel(-999.),
  fDispersion(-1),
  fM20(0.),
  fM01(0.),
  fM02(0.),
  fM11(0.),
  fEmcCpvDistance(-999.),
  fNExMax(0),
  fTracksMatched(),
  fCellNumber()
{
  // constructor
}


//______________________________________________________________________________
AliAODCaloCluster::~AliAODCaloCluster() 
{
  // destructor
}


//______________________________________________________________________________
AliAODCaloCluster::AliAODCaloCluster(const AliAODCaloCluster& clus) :
  AliAODCluster(clus),
  fDistToBadChannel(clus.fDistToBadChannel),
  fDispersion(clus.fDispersion),
  fM20(clus.fM20),
  fM01(clus.fM01),
  fM02(clus.fM02),
  fM11(clus.fM11),
  fEmcCpvDistance(clus.fEmcCpvDistance),
  fNExMax(clus.fNExMax),
  fTracksMatched(clus.fTracksMatched),
  fCellNumber(clus.fCellNumber)
{
  // Copy constructor
}

//______________________________________________________________________________
AliAODCaloCluster& AliAODCaloCluster::operator=(const AliAODCaloCluster& clus)
{
  // Assignment operator
  if(this!=&clus) {

    AliAODCluster::operator=(clus);

    fDistToBadChannel = clus.fDistToBadChannel;
    fDispersion = clus.fDispersion;
    fM20 = clus.fM20;
    fM01 = clus.fM01;
    fM02 = clus.fM02;
    fM11 = clus.fM11;
    fEmcCpvDistance = clus.fEmcCpvDistance;
    fNExMax = clus.fNExMax;
    fTracksMatched = clus.fTracksMatched;
    fCellNumber = clus.fCellNumber;
  }

  return *this;
}

//_______________________________________________________________________
Bool_t AliAODCaloCluster::HasTrackMatched(TObject *trk) const
{
  // Checks if the given track contributed to this cluster.

  TRefArrayIter iter(&fTracksMatched);
  while (TObject *track = iter.Next()) {
    if (trk == track) return kTRUE;
  }
  return kFALSE;
}
