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
/* $Log $ */

//-----------------------------------------------------------------
//           Implementation of the ESD Calorimeter cluster class
//   ESD = Event Summary Data
//   This is the class to deal with during the phisics analysis of data
//
//   J.L. Klay (LLNL)
//-----------------------------------------------------------------

#include <TLorentzVector.h>
#include "AliESDCaloCluster.h"

ClassImp(AliESDCaloCluster)

//_______________________________________________________________________
AliESDCaloCluster::AliESDCaloCluster() : 
  AliVCluster(),
  fTracksMatched(0x0),
  fLabels(0x0),
  fNCells(0),
  fCellsAbsId(0x0),
  fCellsAmpFraction(0x0),
  fEnergy(0),
  fDispersion(0),
  fChi2(0),
  fM20(0),
  fM02(0),
  fEmcCpvDistance(1024),
  fTrackDx(1024),fTrackDz(1024),
  fDistToBadChannel(1024),
  fID(0),
  fNExMax(0),
  fClusterType(kUndef), 
  fTOF(0.),
  fMCEnergyFraction(0.)
{
  //
  // The default ESD constructor 
  //
  fGlobalPos[0] = fGlobalPos[1] = fGlobalPos[2] = 0.;
  for(Int_t i=0; i<AliPID::kSPECIESCN; i++) fPID[i] = 0.;
}

//_______________________________________________________________________
AliESDCaloCluster::AliESDCaloCluster(const AliESDCaloCluster& clus) : 
  AliVCluster(clus),
  fTracksMatched(clus.fTracksMatched?new TArrayI(*clus.fTracksMatched):0x0),
  fLabels(clus.fLabels?new TArrayI(*clus.fLabels):0x0),
  fNCells(clus.fNCells),
  fCellsAbsId(),
  fCellsAmpFraction(),
  fEnergy(clus.fEnergy),
  fDispersion(clus.fDispersion),
  fChi2(clus.fChi2),
  fM20(clus.fM20),
  fM02(clus.fM02),
  fEmcCpvDistance(clus.fEmcCpvDistance),
  fTrackDx(clus.fTrackDx),
  fTrackDz(clus.fTrackDz),
  fDistToBadChannel(clus.fDistToBadChannel),
  fID(clus.fID),
  fNExMax(clus.fNExMax),
  fClusterType(clus.fClusterType),
  fTOF(clus.fTOF),
  fMCEnergyFraction(clus.fMCEnergyFraction)
{
  //
  // The copy constructor 
  //
  fGlobalPos[0] = clus.fGlobalPos[0];
  fGlobalPos[1] = clus.fGlobalPos[1];
  fGlobalPos[2] = clus.fGlobalPos[2];

  for(Int_t i=0; i<AliPID::kSPECIESCN; i++) fPID[i] = clus.fPID[i];

  if (clus.fNCells > 0) {

    if(clus.fCellsAbsId){
      fCellsAbsId = new UShort_t[clus.fNCells];
      for (Int_t i=0; i<clus.fNCells; i++)
	fCellsAbsId[i]=clus.fCellsAbsId[i];
    }
    
    if(clus.fCellsAmpFraction){
      fCellsAmpFraction = new Double32_t[clus.fNCells];
      for (Int_t i=0; i<clus.fNCells; i++)
	fCellsAmpFraction[i]=clus.fCellsAmpFraction[i];
    }
    
  }

}

//_______________________________________________________________________
AliESDCaloCluster &AliESDCaloCluster::operator=(const AliESDCaloCluster& source)
{
  // assignment operator

  if(&source == this) return *this;
  AliVCluster::operator=(source);
  fGlobalPos[0] = source.fGlobalPos[0];
  fGlobalPos[1] = source.fGlobalPos[1];
  fGlobalPos[2] = source.fGlobalPos[2];

  fEnergy = source.fEnergy;
  fDispersion = source.fDispersion;
  fChi2 = source.fChi2;
  fM20 = source.fM20;
  fM02 = source.fM02;
  fEmcCpvDistance = source.fEmcCpvDistance;
  fTrackDx= source.fTrackDx ;
  fTrackDz= source.fTrackDz ;
  fDistToBadChannel = source.fDistToBadChannel ;
  for(Int_t i=0; i<AliPID::kSPECIESCN; i++) fPID[i] = source.fPID[i];
  fID = source.fID;

  fNCells= source.fNCells;

  if (source.fNCells > 0) {
    if(source.fCellsAbsId){
      if(fNCells != source.fNCells||!fCellsAbsId){
	if(fCellsAbsId)delete [] fCellsAbsId;
	fCellsAbsId = new UShort_t[source.fNCells];
      }
      for (Int_t i=0; i<source.fNCells; i++){
	fCellsAbsId[i]=source.fCellsAbsId[i];
      }
    }
    
    if(source.fCellsAmpFraction){
      if(fNCells != source.fNCells||!fCellsAmpFraction){
	if(fCellsAmpFraction) delete [] fCellsAmpFraction;
	fCellsAmpFraction = new Double32_t[source.fNCells];
      }
      for (Int_t i=0; i<source.fNCells; i++)
	fCellsAmpFraction[i]=source.fCellsAmpFraction[i];
    }  
  }

  fNExMax = source.fNExMax;
  fClusterType = source.fClusterType;
  fTOF = source.fTOF;

  //not in use
  if(source.fTracksMatched){
    // assign or copy construct
    if(fTracksMatched){
      *fTracksMatched = *source.fTracksMatched;
    }
    else fTracksMatched = new TArrayI(*source.fTracksMatched);
  }
  else{
    if(fTracksMatched)delete fTracksMatched;
    fTracksMatched = 0;
  }

  if(source.fLabels){
    // assign or copy construct
    if(fLabels){ 
      *fLabels = *source.fLabels;
    }
    else fLabels = new TArrayI(*source.fLabels);
  }
  else{
    if(fLabels)delete fLabels;
    fLabels = 0;
  }

  fMCEnergyFraction = source.fMCEnergyFraction;
  
  return *this;

}

//_______________________________________________________________________
void AliESDCaloCluster::Copy(TObject &obj) const {
  
  // this overwrites the virtual TOBject::Copy()
  // to allow run time copying without casting
  // in AliESDEvent

  if(this==&obj)return;
  AliESDCaloCluster *robj = dynamic_cast<AliESDCaloCluster*>(&obj);
  if(!robj)return; // not an AliESDCluster
  *robj = *this;

}

//_______________________________________________________________________
AliESDCaloCluster::~AliESDCaloCluster(){ 
  //
  // This is destructor according Coding Conventions 
  //
  if(fTracksMatched)delete fTracksMatched;fTracksMatched = 0;
  if(fLabels) delete fLabels; fLabels = 0;
  if(fCellsAmpFraction){ delete[] fCellsAmpFraction; fCellsAmpFraction=0;}
  if(fCellsAbsId){ delete[] fCellsAbsId;  fCellsAbsId = 0;}
}

//_______________________________________________________________________
void AliESDCaloCluster::Clear(const Option_t*){ 
  //
  // This is destructor according Coding Conventions 
  //
  if(fTracksMatched)delete fTracksMatched;fTracksMatched = 0;
  if(fLabels) delete fLabels; fLabels = 0;
  if(fCellsAmpFraction){ delete[] fCellsAmpFraction; fCellsAmpFraction=0;}
  if(fCellsAbsId){ delete[] fCellsAbsId;  fCellsAbsId = 0;}
}


//_______________________________________________________________________
void AliESDCaloCluster::SetPID(const Float_t *p) {
  // Sets the probability of each particle type
  // Copied from AliESDtrack SetPIDValues
  // This function copies "n" PID weights from "scr" to "dest"
  // and normalizes their sum to 1 thus producing conditional
  // probabilities.
  // The negative weights are set to 0.
  // In case all the weights are non-positive they are replaced by
  // uniform probabilities

  Int_t n = AliPID::kSPECIESCN;

  Float_t uniform = 1./(Float_t)n;

  Float_t sum = 0;
  for (Int_t i=0; i<n; i++)
    if (p[i]>=0) {
      sum+=p[i];
      fPID[i] = p[i];
    }
    else {
      fPID[i] = 0;
    }

  if(sum>0)
    for (Int_t i=0; i<n; i++) fPID[i] /= sum;
  else
    for (Int_t i=0; i<n; i++) fPID[i] = uniform;

}

//_______________________________________________________________________
void AliESDCaloCluster::GetMomentum(TLorentzVector& p, Double_t *vertex ) {
  // Returns TLorentzVector with momentum of the cluster. Only valid for clusters 
  // identified as photons or pi0 (overlapped gamma) produced on the vertex
  //Vertex can be recovered with esd pointer doing:  
  //" Double_t vertex[3] ; esd->GetVertex()->GetXYZ(vertex) ; "

  Double32_t pos[3]={ fGlobalPos[0], fGlobalPos[1], fGlobalPos[2]};
  if(vertex){//calculate direction from vertex
    pos[0]-=vertex[0];
    pos[1]-=vertex[1];
    pos[2]-=vertex[2];
  }
  
  Double_t r = TMath::Sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]   ) ; 
  
  p.SetPxPyPzE( fEnergy*pos[0]/r,  fEnergy*pos[1]/r,  fEnergy*pos[2]/r,  fEnergy) ;   
}

//_______________________________________________________________________
void  AliESDCaloCluster::SetCellsAbsId(UShort_t *array)
{
    //  Set the array of cell absId numbers 
  if (fNCells) {
    fCellsAbsId = new  UShort_t[fNCells];
    for (Int_t i = 0; i < fNCells; i++) fCellsAbsId[i] = array[i];
  }
}

//_______________________________________________________________________
void  AliESDCaloCluster::SetCellsAmplitudeFraction(Double32_t *array)
{
  //  Set the array of cell amplitude fraction
  if (fNCells) {
    fCellsAmpFraction = new  Double32_t[fNCells];
    for (Int_t i = 0; i < fNCells; i++) fCellsAmpFraction[i] = array[i];
  }
}

//______________________________________________________________________________
void AliESDCaloCluster::SetPosition(Float_t *x) 
{
  // set the position
  
  if (x) {
    fGlobalPos[0] = x[0];
    fGlobalPos[1] = x[1];
    fGlobalPos[2] = x[2];
  } else {
    
    fGlobalPos[0] = -999.;
    fGlobalPos[1] = -999.;
    fGlobalPos[2] = -999.;
  }
}


