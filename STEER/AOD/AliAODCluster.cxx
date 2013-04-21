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
  AliVCluster(),
  fEnergy(0),
  fChi2(-999.),
  fID(-999),
  fNLabel(0),
  fLabel(0x0),
  fFilterMap(0),
  fType(kUndef),
  fMCEnergyFraction(0.)
{
  // default constructor

  SetPosition(NULL);
  SetPID(NULL);
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
  AliVCluster(),
  fEnergy(energy),
  fChi2(-999.),
  fID(id),
  fNLabel(0),
  fLabel(0x0),
  fFilterMap(selectInfo),
  fType(ttype),
  fMCEnergyFraction(0.)
{
  // constructor
  for (Int_t i = 0; i <  3; i++) fPosition[i] = 0.;
  for (Int_t i = 0; i < 13; i++) fPID[i]      = 0;
 
  if(x)   {for (Int_t i = 0; i < 3  ; i++) SetPositionAt(x[i]  ,i);}
  if(pid) {for (Int_t i = 0; i < 13 ; i++) SetPIDAt     (pid[i],i);}
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
  AliVCluster(),
  fEnergy(energy),
  fChi2(-999.),
  fID(id),
  fNLabel(0),
  fLabel(0x0),
  fFilterMap(selectInfo),
  fType(ttype),
  fMCEnergyFraction(0.)
{
  // constructor
  for (Int_t i = 0; i <  3; i++) fPosition[i] = 0.;
  for (Int_t i = 0; i < 13; i++) fPID[i]      = 0;

  if(x)   {for (Int_t i = 0; i < 3  ; i++) SetPositionAt(x[i]  ,i);}
  if(pid) {for (Int_t i = 0; i < 13 ; i++) SetPIDAt     (pid[i],i);}
  SetLabel(label, nLabel);
}


//______________________________________________________________________________
AliAODCluster::~AliAODCluster() 
{
  // destructor

  RemoveLabel();
}

//______________________________________________________________________________
void AliAODCluster::Clear(const Option_t*) 
{
  // Clear
  
  RemoveLabel();
}


//______________________________________________________________________________
AliAODCluster::AliAODCluster(const AliAODCluster& clus) :
  AliVCluster(clus),
  fEnergy(clus.fEnergy),
  fChi2(clus.fChi2),
  fID(clus.fID),
  fNLabel(0),
  fLabel(0x0),
  fFilterMap(clus.fFilterMap),
  fType(clus.fType),
  fMCEnergyFraction(clus.fMCEnergyFraction)
{
  // Copy constructor

  for(Int_t i = 0; i < 3  ; i++) fPosition[i]  = clus.fPosition[i];
  for(Int_t i = 0; i < 13 ; i++)  fPID[i]      = clus.fPID[i];

  SetLabel(clus.fLabel, clus.fNLabel);
}

//______________________________________________________________________________
AliAODCluster& AliAODCluster::operator=(const AliAODCluster& clus)
{
  // Assignment operator
  if(this!=&clus) {
    
	for(Int_t i = 0; i < 3 ;  i++) fPosition[i] = clus.fPosition[i];
	for(Int_t i = 0; i < 13 ; i++) fPID[i]      = clus.fPID[i];
    
    fEnergy = clus.fEnergy;
    fChi2 = clus.fChi2;

    fID = clus.fID;
    SetLabel(clus.fLabel, clus.fNLabel);
    fFilterMap = clus.fFilterMap;

    fType = clus.fType;

    fMCEnergyFraction = clus.fMCEnergyFraction;
  }

  return *this;
}

//______________________________________________________________________________
void AliAODCluster::SetPosition(Float_t *x) 
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
UShort_t AliAODCluster::GetMostProbablePID() const 
{
  // Returns the most probable PID array element.
  
  Int_t nPID = 13;
  UShort_t unknown = AliVCluster::kUnknown;
  
  UShort_t loc = unknown;
  Double_t max = 0.;
  Bool_t allTheSame = kTRUE;
  
  for (Int_t iPID = 0; iPID < nPID; iPID++) {
    if (fPID[iPID] >= max) {
      if (fPID[iPID] > max) {
	allTheSame = kFALSE;
	max = fPID[iPID];
	loc = (UShort_t)iPID;
      } else {
	allTheSame = kTRUE;
      }
    }
  }
  return allTheSame ? unknown : loc;
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
Int_t AliAODCluster::GetLabelAt(UInt_t i) const
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
  const Double_t *pid = GetPID();
  printf("PID weights: photon %0.2f, pi0 %0.2f, electron %0.2f, conversion electron %0.2f\n, hadrons: pion %0.2f, kaon %0.2f, proton %0.2f , neutron %0.2f, kaon %0.2f \n",
	 pid[AliVCluster::kPhoton],   pid[AliVCluster::kPi0],
	 pid[AliVCluster::kElectron], pid[AliVCluster::kEleCon],
	 pid[AliVCluster::kPion],     pid[AliVCluster::kKaon],   pid[AliVCluster::kProton],
	 pid[AliVCluster::kNeutron],  pid[AliVCluster::kKaon0]);
}
