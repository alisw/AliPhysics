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

#include <TLorentzVector.h>
#include "AliLog.h"
#include "AliAODCaloCluster.h"

/// \cond CLASSIMP
ClassImp(AliAODCaloCluster) ;
/// \endcond

///
/// Default constructor.
///
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
  fIsExotic(kFALSE),
  fCellsMCEdepFractionMap(0x0)
{
  for (Int_t i = 0; i <= kLastUserDefEnergy; i++) {
    fUserDefEnergy[i] = 1.;
  }
}

///
/// Constructor.
///
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
  fIsExotic(kFALSE),
  fCellsMCEdepFractionMap(0x0)
{
  for (Int_t i = 0; i <= kLastUserDefEnergy; i++) {
    fUserDefEnergy[i] = 1.;
  }
}

///
/// Constructor.
///
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
  fIsExotic(kFALSE),
  fCellsMCEdepFractionMap(0x0)
{
  for (Int_t i = 0; i <= kLastUserDefEnergy; i++) {
    fUserDefEnergy[i] = 1.;
  }
}

///
/// Destructor.
///
//______________________________________________________________________________
AliAODCaloCluster::~AliAODCaloCluster() 
{
  if(fCellsAmpFraction)       delete[] fCellsAmpFraction;       fCellsAmpFraction       = 0 ;
  if(fCellsAbsId)             delete[] fCellsAbsId;             fCellsAbsId             = 0 ;
  if(fCellsMCEdepFractionMap) delete[] fCellsMCEdepFractionMap; fCellsMCEdepFractionMap = 0 ;
}

///
/// Clear array pointers.
///
//______________________________________________________________________________
void AliAODCaloCluster::Clear(const Option_t*) 
{
  RemoveLabel();
  if(fCellsAmpFraction)       delete[] fCellsAmpFraction;       fCellsAmpFraction       = 0 ;
  if(fCellsAbsId)             delete[] fCellsAbsId;             fCellsAbsId             = 0 ;
  if(fCellsMCEdepFractionMap) delete[] fCellsMCEdepFractionMap; fCellsMCEdepFractionMap = 0 ;
}

///
/// Copy constructor.
///
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
  fIsExotic(clus.fIsExotic),
  fCellsMCEdepFractionMap(0x0)
{
  if (clus.fNCells > 0) 
  {
    if(clus.fCellsAbsId)
    {
      fCellsAbsId = new UShort_t[clus.fNCells];
      for (Int_t i=0; i<clus.fNCells; i++)
        fCellsAbsId[i]=clus.fCellsAbsId[i];
    }
    
    if(clus.fCellsAmpFraction)
    {
      fCellsAmpFraction = new Double32_t[clus.fNCells];
      for (Int_t i=0; i<clus.fNCells; i++)
        fCellsAmpFraction[i]=clus.fCellsAmpFraction[i];
    }
    
    if(clus.fCellsMCEdepFractionMap)
    {
      fCellsMCEdepFractionMap = new UInt_t[clus.fNCells];
      for (Int_t i=0; i<clus.fNCells; i++) fCellsMCEdepFractionMap[i]=clus.fCellsMCEdepFractionMap[i];
    }
  }
  
  for (Int_t i = 0; i <= kLastUserDefEnergy; i++) 
    fUserDefEnergy[i] = clus.fUserDefEnergy[i];
  
}

///
/// Assignment operator.
///
//______________________________________________________________________________
AliAODCaloCluster& AliAODCaloCluster::operator=(const AliAODCaloCluster& clus)
{
  if(this!=&clus) 
  {
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
    if(fCellsAbsId)             delete [] fCellsAbsId;
    if(fCellsAmpFraction)       delete [] fCellsAmpFraction;
    if(fCellsMCEdepFractionMap) delete [] fCellsMCEdepFractionMap;
    
    if (clus.fNCells > 0) 
    {
      if(clus.fCellsAbsId)
      {
        fCellsAbsId = new UShort_t[clus.fNCells];
        for (Int_t i=0; i<clus.fNCells; i++)
          fCellsAbsId[i]=clus.fCellsAbsId[i];
      }
      
      if(clus.fCellsAmpFraction)
      {
        fCellsAmpFraction = new Double32_t[clus.fNCells];
        for (Int_t i=0; i<clus.fNCells; i++)
          fCellsAmpFraction[i]=clus.fCellsAmpFraction[i];
      }
      
      if(clus.fCellsMCEdepFractionMap)
      {
        fCellsMCEdepFractionMap = new UInt_t[clus.fNCells];
        for (Int_t i=0; i<clus.fNCells; i++) 
          fCellsMCEdepFractionMap[i]=clus.fCellsMCEdepFractionMap[i];
      }
    } // fNCells > 0
    
    fMCEnergyFraction = clus.fMCEnergyFraction;
    fIsExotic = clus.fIsExotic;
    
    for (Int_t i = 0; i <= kLastUserDefEnergy; i++) 
      fUserDefEnergy[i] = clus.fUserDefEnergy[i];
  }
  
  return *this;
}

///
/// Checks if the given track contributed to this cluster.
///
//_______________________________________________________________________
Bool_t AliAODCaloCluster::HasTrackMatched(TObject *trk) const
{
  TRefArrayIter iter(&fTracksMatched);
  while (TObject *track = iter.Next()) {
    if (trk == track) return kTRUE;
  }
  return kFALSE;
}

///
/// Returns TLorentzVector with momentum of the cluster. Only valid for clusters 
/// identified as photons or pi0 (overlapped gamma) produced on the vertex
/// Vertex can be recovered with esd pointer doing:  
///  " Double_t vertex[3] ; esd->GetVertex()->GetXYZ(vertex) ; "
///
//_______________________________________________________________________
void AliAODCaloCluster::GetMomentum(TLorentzVector& p, Double_t *vertex ) const 
{
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

///
/// Returns TLorentzVector with momentum of the cluster. Only valid for clusters 
/// identified as photons or pi0 (overlapped gamma) produced on the vertex
/// Uses the user defined energy t
/// Vertex can be recovered with esd pointer doing:  
///  " Double_t vertex[3] ; esd->GetVertex()->GetXYZ(vertex) ; "
///
//_______________________________________________________________________
void AliAODCaloCluster::GetMomentum(TLorentzVector& p, Double_t *vertex, VCluUserDefEnergy_t t ) const 
{
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

///
///  Set the array of cell absolute Id. numbers. 
///
//______________________________________________________________________________
void  AliAODCaloCluster::SetCellsAbsId(UShort_t *array)
{
    if (fNCells) {
      if(!fCellsAbsId)fCellsAbsId = new  UShort_t[fNCells];
      for (Int_t i = 0; i < fNCells; i++) fCellsAbsId[i] = array[i];
    }
}

///
///  Set the array of cell amplitude fractions. 
///  Cell can be shared between 2 clusters, here the fraction of energy
///  assigned to each cluster is stored. Only in unfolded clusters.
///
//______________________________________________________________________________
void  AliAODCaloCluster::SetCellsAmplitudeFraction(Double32_t *array)
{
    if (fNCells) {
      if(!fCellsAmpFraction)fCellsAmpFraction = new  Double32_t[fNCells];
      for (Int_t i = 0; i < fNCells; i++) fCellsAmpFraction[i] = array[i];
    }
}

///
/// \param cellIndex: position of cell in array fCellsAbsId
/// \param eDep: Filled float array with 4 entries, each is the fraction of deposited 
///              energy by 4 most significant MC particles (GetLabels()) in a cell of the cluster.
/// In this method, the 4 fractions  stored in % values (0 to 100) 
/// in each bit of the integer fCellsMCEdepFractionMap[cellIndex] are unpacked. 
//______________________________________________________________________________
void  AliAODCaloCluster::GetCellMCEdepFractionArray(Int_t cellIndex, Float_t * eDep) const
{ 
  if ( cellIndex >= fNCells || fNCells < 0 || !fCellsMCEdepFractionMap )
  {
    eDep[0] = eDep[1] = eDep[2] = eDep[3] = 0. ;
    return;
  }
  
  eDep[0] =  (fCellsMCEdepFractionMap[cellIndex]&0x000000ff)        / 100.;
  eDep[1] = ((fCellsMCEdepFractionMap[cellIndex]&0x0000ff00) >>  8) / 100.;
  eDep[2] = ((fCellsMCEdepFractionMap[cellIndex]&0x00ff0000) >> 16) / 100.;
  eDep[3] = ((fCellsMCEdepFractionMap[cellIndex]&0xff000000) >> 24) / 100.;  
}

///
/// \param eDep: Float array with 4 entries, each is the fraction of deposited 
///              energy by an MC particle in a cell of the cluster.
/// 
/// The MC particle must correspond one of the 4 first labels in GetLabels(). This method
/// packs the 4 floats into an integer, assigning each bit a value between 0 and 100
//______________________________________________________________________________
UInt_t  AliAODCaloCluster::PackMCEdepFraction(Float_t * eDep) const
{ 
  UInt_t intEDep[4];
  
  for(Int_t i = 0; i < 4; i++)
    intEDep[i] = TMath::Nint(eDep[i]*100) ;
  
  UInt_t map = intEDep[0]|(intEDep[1]<<8)|(intEDep[2]<<16)|(intEDep[3]<<24);
  
  return map;
}

///
/// Set the array with the fraction of deposited energy in a cell belonging to 
/// the cluster by a given primary  particle. Each entry of the array corresponds 
/// to the same entry in fCellsAbsId. Each entry is an integer where a maximum 
/// of 4 energy deposition fractions are encoded, each corresponding to the 
/// first 4 entries in GetLabels()
//______________________________________________________________________________
void  AliAODCaloCluster::SetCellsMCEdepFractionMap(UInt_t *array)
{
  if ( fNCells <= 0 || !array ) return; 
  
  fCellsMCEdepFractionMap = new  UInt_t[fNCells];
  
  for (Int_t i = 0; i < fNCells; i++) 
    fCellsMCEdepFractionMap[i] = array[i];
}



