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

//-------------------------------------------------------------------------
//     AOD calorimeter cluster class (for PHOS and EMCAL)
//     Author: Markus Oldenburg, CERN
//-------------------------------------------------------------------------

#include <TLorentzVector.h>
#include "AliLog.h"
#include "AliAODCaloCluster.h"

ClassImp(AliAODCaloCluster)

//______________________________________________________________________________
AliAODCaloCluster::AliAODCaloCluster() : 
  AliAODCluster(),
  fDistToBadChannel(-999.),
  fDispersion(-1),
  fM20(0.),
  fM02(0.),
  fEmcCpvDistance(-999.),
  fTrackDx(-999),
  fTrackDz(-999),
  fNExMax(0), 
  fTOF(0.),
  fCoreEnergy(0.),
  fTracksMatched(),
  fNCells(0),
  fCellsAbsId(0x0),
  fCellsAmpFraction(0x0),
  fMCEnergyFraction(0.),
  fIsExotic(kFALSE)
{
  // default constructor

  for (Int_t i = 0; i <= kLastUserDefEnergy; i++) {
    fUserDefEnergy[i] = 1.;
  }
}

//______________________________________________________________________________
AliAODCaloCluster::AliAODCaloCluster(Int_t id,
				     UInt_t nLabel,
				     Int_t *label, 
				     Double_t energy,
				     Double_t x[3],
				     Double_t pid[13],
				     Char_t ttype,
				     UInt_t selectInfo) :
  AliAODCluster(id, nLabel, label, energy, x, pid, ttype, selectInfo),
  fDistToBadChannel(-999.),
  fDispersion(-1),
  fM20(0.),
  fM02(0.),
  fEmcCpvDistance(-999.),
  fTrackDx(-999),
  fTrackDz(-999),
  fNExMax(0),
  fTOF(0.),
  fCoreEnergy(0.),
  fTracksMatched(),
  fNCells(0),
  fCellsAbsId(0x0),
  fCellsAmpFraction(0x0),
  fMCEnergyFraction(0.),
  fIsExotic(kFALSE)
{
  // constructor

  for (Int_t i = 0; i <= kLastUserDefEnergy; i++) {
    fUserDefEnergy[i] = 1.;
  }
}

//______________________________________________________________________________
AliAODCaloCluster::AliAODCaloCluster(Int_t id,
				     UInt_t nLabel,
				     Int_t *label, 
				     Float_t energy,
				     Float_t x[3],
				     Float_t pid[13],
				     Char_t ttype,
				     UInt_t selectInfo) :
  AliAODCluster(id, nLabel, label, energy, x, pid, ttype, selectInfo),
  fDistToBadChannel(-999.),
  fDispersion(-1),
  fM20(0.),
  fM02(0.),
  fEmcCpvDistance(-999.),
  fTrackDx(-999),
  fTrackDz(-999),
  fNExMax(0),
  fTOF(0.),
  fCoreEnergy(0.),
  fTracksMatched(),
  fNCells(0),
  fCellsAbsId(0x0),
  fCellsAmpFraction(0x0),
  fMCEnergyFraction(0.),
  fIsExotic(kFALSE)
{
  // constructor

  for (Int_t i = 0; i <= kLastUserDefEnergy; i++) {
    fUserDefEnergy[i] = 1.;
  }
}


//______________________________________________________________________________
AliAODCaloCluster::~AliAODCaloCluster() 
{
  // destructor
  if(fCellsAmpFraction) delete[] fCellsAmpFraction; fCellsAmpFraction=0;
  if(fCellsAbsId)       delete[] fCellsAbsId;       fCellsAbsId = 0;
}

//______________________________________________________________________________
void AliAODCaloCluster::Clear(const Option_t*) 
{
  // clear
  RemoveLabel();
  if(fCellsAmpFraction) delete[] fCellsAmpFraction; fCellsAmpFraction=0;
  if(fCellsAbsId)       delete[] fCellsAbsId;       fCellsAbsId = 0;
}

//______________________________________________________________________________
AliAODCaloCluster::AliAODCaloCluster(const AliAODCaloCluster& clus) :
  AliAODCluster(clus),
  fDistToBadChannel(clus.fDistToBadChannel),
  fDispersion(clus.fDispersion),
  fM20(clus.fM20),
  fM02(clus.fM02),
  fEmcCpvDistance(clus.fEmcCpvDistance),
  fTrackDx(clus.fTrackDx),
  fTrackDz(clus.fTrackDz),
  fNExMax(clus.fNExMax),
  fTOF(clus.fTOF),
  fCoreEnergy(clus.fCoreEnergy),
  fTracksMatched(clus.fTracksMatched),
  fNCells(clus.fNCells),
  fCellsAbsId(0x0),
  fCellsAmpFraction(0x0),
  fMCEnergyFraction(clus.fMCEnergyFraction),
  fIsExotic(clus.fIsExotic)
{
  // Copy constructor

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

  for (Int_t i = 0; i <= kLastUserDefEnergy; i++) {
    fUserDefEnergy[i] = clus.fUserDefEnergy[i];
  }
  
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
    fM02 = clus.fM02;
    fEmcCpvDistance = clus.fEmcCpvDistance;
    fTrackDx=clus.fTrackDx;
    fTrackDz=clus.fTrackDz;
    fNExMax = clus.fNExMax;
    fTOF = clus.fTOF;
    fCoreEnergy = clus.fCoreEnergy;
    fTracksMatched = clus.fTracksMatched;

    fNCells= clus. fNCells;
    // delete anyway 
    if(fCellsAbsId)delete [] fCellsAbsId;
    if(fCellsAmpFraction)delete [] fCellsAmpFraction;

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

  fMCEnergyFraction = clus.fMCEnergyFraction;
  fIsExotic = clus.fIsExotic;

  for (Int_t i = 0; i <= kLastUserDefEnergy; i++) {
    fUserDefEnergy[i] = clus.fUserDefEnergy[i];
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

//_______________________________________________________________________
void AliAODCaloCluster::GetMomentum(TLorentzVector& p, Double_t *vertex ) const {
  // Returns TLorentzVector with momentum of the cluster. Only valid for clusters 
  // identified as photons or pi0 (overlapped gamma) produced on the vertex
  //Vertex can be recovered with esd pointer doing:  
  //" Double_t vertex[3] ; esd->GetVertex()->GetXYZ(vertex) ; "

  Double32_t energy = E();
  Float_t    pos[3];
  GetPosition(pos);
  
  if(vertex){//calculate direction from vertex
    pos[0]-=vertex[0];
    pos[1]-=vertex[1];
    pos[2]-=vertex[2];
  }
  
  Double_t r = TMath::Sqrt(pos[0]*pos[0]+
			   pos[1]*pos[1]+
			   pos[2]*pos[2]   ) ; 
     
  if ( r > 0 ) 
    p.SetPxPyPzE( energy*pos[0]/r,  energy*pos[1]/r,  energy*pos[2]/r,  energy) ; 
  else
    AliInfo("Null cluster radius, momentum calculation not possible");
}

//_______________________________________________________________________
void AliAODCaloCluster::GetMomentum(TLorentzVector& p, Double_t *vertex, VCluUserDefEnergy_t t ) const {
  // Returns TLorentzVector with momentum of the cluster. Only valid for clusters 
  // identified as photons or pi0 (overlapped gamma) produced on the vertex
  // Uses the user defined energy t
  //Vertex can be recovered with esd pointer doing:  
  //" Double_t vertex[3] ; esd->GetVertex()->GetXYZ(vertex) ; "

  Double32_t energy = GetUserDefEnergy(t);
  Float_t    pos[3];
  GetPosition(pos);
  
  if(vertex){//calculate direction from vertex
    pos[0]-=vertex[0];
    pos[1]-=vertex[1];
    pos[2]-=vertex[2];
  }
  
  Double_t r = TMath::Sqrt(pos[0]*pos[0]+
			   pos[1]*pos[1]+
			   pos[2]*pos[2]   ) ; 
  
  p.SetPxPyPzE( energy*pos[0]/r,  energy*pos[1]/r,  energy*pos[2]/r,  energy) ; 
  
}


void  AliAODCaloCluster::SetCellsAbsId(UShort_t *array)
{
    //  Set the array of cell absId numbers 
    if (fNCells) {
      if(!fCellsAbsId)fCellsAbsId = new  UShort_t[fNCells];
      for (Int_t i = 0; i < fNCells; i++) fCellsAbsId[i] = array[i];
    }
}

void  AliAODCaloCluster::SetCellsAmplitudeFraction(Double32_t *array)
{
    //  Set the array of cell amplitude fraction
    if (fNCells) {
      if(!fCellsAmpFraction)fCellsAmpFraction = new  Double32_t[fNCells];
      for (Int_t i = 0; i < fNCells; i++) fCellsAmpFraction[i] = array[i];
    }
}
