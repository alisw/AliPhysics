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
//     AOD cluster base class
//     Author: Markus Oldenburg, CERN
//-------------------------------------------------------------------------

#include "AliAODCluster.h"

ClassImp(AliAODCluster)

//______________________________________________________________________________
AliAODCluster::AliAODCluster() : 
  fEnergy(0),
  fChi2(-999.),
  fID(-999),
  fNLabel(0),
  fLabel(0x0),
  fFilterMap(0),
  fType(kUndef)
{
  // default constructor

  SetPosition((Float_t*)NULL);
  SetPID((Float_t*)NULL);
}

//______________________________________________________________________________
AliAODCluster::AliAODCluster(Int_t id,
			     UInt_t nLabel,
			     Int_t *label, 
			     Double_t energy,
			     Double_t x[3],
			     Double_t pid[13],
			     Char_t ttype,
			     UInt_t selectInfo) :
  fEnergy(energy),
  fChi2(-999.),
  fID(id),
  fNLabel(0),
  fLabel(0x0),
  fFilterMap(selectInfo),
  fType(ttype)
{
  // constructor
 
  SetPosition(x);
  SetPID(pid);
  SetLabel(label, nLabel);
}

//______________________________________________________________________________
AliAODCluster::AliAODCluster(Int_t id,
			     UInt_t nLabel,
			     Int_t *label, 
			     Float_t energy,
			     Float_t x[3],
			     Float_t pid[13],
			     Char_t ttype,
			     UInt_t selectInfo) :
  fEnergy(energy),
  fChi2(-999.),
  fID(id),
  fNLabel(0),
  fLabel(0x0),
  fFilterMap(selectInfo),
  fType(ttype)
{
  // constructor
 
  SetPosition(x);
  SetPID(pid);
  SetLabel(label, nLabel);
}


//______________________________________________________________________________
AliAODCluster::~AliAODCluster() 
{
  // destructor

  RemoveLabel();
}


//______________________________________________________________________________
AliAODCluster::AliAODCluster(const AliAODCluster& clus) :
  TObject(clus),
  fEnergy(clus.fEnergy),
  fChi2(clus.fChi2),
  fID(clus.fID),
  fNLabel(0),
  fLabel(0x0),
  fFilterMap(clus.fFilterMap),
  fType(clus.fType)
{
  // Copy constructor

  clus.GetPosition(fPosition);
  SetPID(clus.fPID);
  SetLabel(clus.fLabel, clus.fNLabel);
}

//______________________________________________________________________________
AliAODCluster& AliAODCluster::operator=(const AliAODCluster& clus)
{
  // Assignment operator
  if(this!=&clus) {

    clus.GetPosition(fPosition);
    clus.GetPID(fPID);

    fEnergy = clus.fEnergy;
    fChi2 = clus.fChi2;

    fID = clus.fID;
    SetLabel(clus.fLabel, clus.fNLabel);
    fFilterMap = clus.fFilterMap;

    fType = clus.fType;
  }

  return *this;
}

//______________________________________________________________________________
template <class T> void AliAODCluster::SetPosition(const T *x) 
{
  // set the position

  if (x) {
      fPosition[0] = x[0];
      fPosition[1] = x[1];
      fPosition[2] = x[2];
  } else {

    fPosition[0] = -999.;
    fPosition[1] = -999.;
    fPosition[2] = -999.;
  }
}

//______________________________________________________________________________
AliAODCluster::AODCluPID_t AliAODCluster::GetMostProbablePID() const 
{
  // Returns the most probable PID array element.
  
  Int_t nPID = 13;
  if (fPID) {
    AODCluPID_t loc = kUnknown;
    Double_t max = 0.;
    Bool_t allTheSame = kTRUE;
    
    for (Int_t iPID = 0; iPID < nPID; iPID++) {
      if (fPID[iPID] >= max) {
	if (fPID[iPID] > max) {
	  allTheSame = kFALSE;
	  max = fPID[iPID];
	  loc = (AODCluPID_t)iPID;
	} else {
	  allTheSame = kTRUE;
	}
      }
    }
    
    return allTheSame ? kUnknown : loc;
  } else {
    return kUnknown;
  }
}

//______________________________________________________________________________
void AliAODCluster::SetLabel(Int_t *label, UInt_t size) 
{
  if (label && size>0) {
    if (size != (UInt_t)fNLabel) {
      RemoveLabel();
      fNLabel = size;
      fLabel = new Int_t[fNLabel];
    }
    
    for (Int_t i = 0; i < fNLabel; i++) {
      fLabel[i] = label[i];
    }
  } else {
    RemoveLabel();
  }

  return;
}

//______________________________________________________________________________
Int_t AliAODCluster::GetLabel(UInt_t i) const
{
  if (fLabel && i < (UInt_t)fNLabel) {
    return fLabel[i];
  } else {
    return -999;
  }
}

//______________________________________________________________________________
void AliAODCluster::RemoveLabel()
{
  delete[] fLabel;
  fLabel = 0x0;
  fNLabel = 0;

  return;
}

//______________________________________________________________________________
void AliAODCluster::Print(Option_t* /* option */) const
{
  // prints information about AliAODCluster

  printf("Cluster type: %d\n", GetType()); 
  printf("     energy = %f\n", E());
  printf("       chi2 = %f\n", Chi2());
  printf("  PID object: %p\n", static_cast<const void*>(PID()));
}
