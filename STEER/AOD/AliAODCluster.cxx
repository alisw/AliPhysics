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

#include "AliAODCluster.h"

/// \cond CLASSIMP
ClassImp(AliAODCluster) ;
/// \endcond

///
/// Default constructor.
///
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
  fMCEnergyFraction(0.),  
  fClusterMCEdepFraction(0x0)
{

  SetPosition(NULL);
  SetPID(NULL);
}

///
/// Constructor.
///
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
  fMCEnergyFraction(0.),  
  fClusterMCEdepFraction(0x0)
{
  for (Int_t i = 0; i <  3; i++) fPosition[i] = 0.;
  for (Int_t i = 0; i < 13; i++) fPID[i]      = 0;
 
  if(x)   {for (Int_t i = 0; i < 3  ; i++) SetPositionAt(x[i]  ,i);}
  if(pid) {for (Int_t i = 0; i < 13 ; i++) SetPIDAt     (pid[i],i);}
  SetLabel(label, nLabel);
}

///
/// Constructor.
///
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
  fMCEnergyFraction(0.),  
  fClusterMCEdepFraction(0x0)
{
  for (Int_t i = 0; i <  3; i++) fPosition[i] = 0.;
  for (Int_t i = 0; i < 13; i++) fPID[i]      = 0;

  if(x)   {for (Int_t i = 0; i < 3  ; i++) SetPositionAt(x[i]  ,i);}
  if(pid) {for (Int_t i = 0; i < 13 ; i++) SetPIDAt     (pid[i],i);}
  SetLabel(label, nLabel);
}

///
/// Destructor
///
//______________________________________________________________________________
AliAODCluster::~AliAODCluster() 
{
  RemoveLabel();
}

///
/// Clear
///
//______________________________________________________________________________
void AliAODCluster::Clear(const Option_t*) 
{  
  RemoveLabel();
}

///
/// Copy constructor.
///
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
  fMCEnergyFraction(clus.fMCEnergyFraction),  
  fClusterMCEdepFraction(0x0)
{
  for(Int_t i = 0; i < 3  ; i++) fPosition[i]  = clus.fPosition[i];
  for(Int_t i = 0; i < 13 ; i++) fPID[i]       = clus.fPID[i];

  SetLabel(clus.fLabel, clus.fNLabel);
  SetClusterMCEdepFraction(clus.fClusterMCEdepFraction);
}

///
/// Assignment operator.
///
//______________________________________________________________________________
AliAODCluster& AliAODCluster::operator=(const AliAODCluster& clus)
{
  if(this!=&clus) {
    
	for(Int_t i = 0; i < 3 ;  i++) fPosition[i] = clus.fPosition[i];
	for(Int_t i = 0; i < 13 ; i++) fPID[i]      = clus.fPID[i];
    
    fEnergy = clus.fEnergy;
    fChi2 = clus.fChi2;

    fID = clus.fID;
    
    SetLabel(clus.fLabel, clus.fNLabel);
    SetClusterMCEdepFraction(clus.fClusterMCEdepFraction);
    
    fFilterMap = clus.fFilterMap;

    fType = clus.fType;

    fMCEnergyFraction = clus.fMCEnergyFraction;
  }

  return *this;
}

///
/// Set cluster global position.
///
//______________________________________________________________________________
void AliAODCluster::SetPosition(Float_t *x) 
{  
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

///
/// \return The most probable PID array element.
///
//______________________________________________________________________________
UShort_t AliAODCluster::GetMostProbablePID() const 
{  
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

///
/// Set the array with MC particle labels and number of labels.
///
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

///
/// Get the label i of the possible fNLabel in the array fLabel
///
//______________________________________________________________________________
Int_t AliAODCluster::GetLabelAt(UInt_t i) const
{
  if (fLabel && i < (UInt_t)fNLabel) {
    return fLabel[i];
  } else {
    return -999;
  }
}

///
/// Delete MC arrays, used in destructor and Clear() and SetLabel() methods
///
//______________________________________________________________________________
void AliAODCluster::RemoveLabel()
{
  delete[] fLabel;
  fLabel = 0x0;
  fNLabel = 0;

  if(fClusterMCEdepFraction)  
    delete[] fClusterMCEdepFraction;  
  fClusterMCEdepFraction  = 0 ;
  
  return;
}

///
/// Prints information about AliAODCluster.
///
//______________________________________________________________________________
void AliAODCluster::Print(Option_t* /* option */) const
{
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

///
/// \return Fraction of deposited energy by the one of the particles in array fLable
/// 
/// \param mcIndex: position of MC particle in array GetLabels()
///
/// The parameter is stored as %, return the corresponding float.
//______________________________________________________________________________
Float_t  AliAODCluster::GetClusterMCEdepFraction(Int_t mcIndex) const
{ 
  if ( mcIndex < 0 ||  mcIndex >= fNLabel ||  !fClusterMCEdepFraction ) return 0. ;

  return  fClusterMCEdepFraction[mcIndex]/100. ; 
}

///
/// Set the array with the fraction of deposited energy in cluster by a given primary 
/// particle. Each entry of the array corresponds to the same entry in GetLabels().
/// Set the fraction in % with respect the cluster energy, store a value between 0 and 100
///
/// Not sure this method is usable for AODs
///
/// \param array: energy deposition array
//______________________________________________________________________________
void  AliAODCluster::SetClusterMCEdepFractionFromEdepArray(Float_t *array)
{
  if ( fNLabel <= 0 || !array ) return ;
  
  fClusterMCEdepFraction = new  UShort_t[fNLabel];
  
  // Get total deposited energy (can be different from reconstructed energy)
  Float_t totalE = 0;
  for (Int_t i = 0; i < fNLabel; i++) totalE+=array[i];

  // Set the fraction of energy per MC contributor in %
  for (Int_t i = 0; i < fNLabel; i++) 
    fClusterMCEdepFraction[i] = TMath::Nint(array[i]/totalE*100.);
}

///
/// Set the array with the fraction of deposited energy in cluster by a given primary 
/// particle. Each entry of the array corresponds to the same entry in GetLabels().
///
/// The fraction must already be in % with respect the cluster energy, store a value between 0 and 100
///
/// Execute after setting of fLable and fNLabel
///
/// \param array: array of fraction of energy deposition / cluster energy 
//______________________________________________________________________________
void  AliAODCluster::SetClusterMCEdepFraction(UShort_t *array)
{
  if ( fNLabel <= 0 || !array ) 
    return ;

  fClusterMCEdepFraction = new  UShort_t[fNLabel];
  
  for (Int_t i = 0; i < fNLabel; i++) 
    fClusterMCEdepFraction[i] = array[i];
}




