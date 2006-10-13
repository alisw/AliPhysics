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

//*-- Author: Renan Cabrera (Creighton U.)

#include "AliEMCALParton.h"
#include "Ecommon.h"
#include <Riostream.h>
ClassImp(AliEMCALParton)   
    
//____________________________________________________________________________
AliEMCALParton::AliEMCALParton()
  : fEnergy(0.), fEta(0.), fPhi(0.),
    fNTracks(0), fPartonCode(0)
{
  // Default constructor
}

AliEMCALParton::AliEMCALParton(Float_t energy, Float_t phi, Float_t eta)
  : fEnergy(energy), fEta(eta), fPhi(phi),
    fNTracks(0), fPartonCode(0)
{
  // Constructor
}

void AliEMCALParton::SetTrackList(Int_t NTracks, Float_t* Energy,  Float_t* Eta, Float_t* Phi, Int_t* PDG)
{
// Set the stored tracklist
  fNTracks     = NTracks;
  for (Int_t i=0;i<NTracks;i++)
  {
    fTrackEnergy[i] = Energy[i] ;
    fTrackEta[i]    = Eta[i];
    fTrackPhi[i]    = Phi[i];
    fTrackPDG[i]    = PDG[i];
  }
}

void AliEMCALParton::GetTrackList(Float_t* Energy,  Float_t* Eta, Float_t* Phi, Int_t* PDG) const
{
// retrieves the stored tracklist	
  for (Int_t i=0;i<fNTracks;i++)
  {
    Energy[i] = fTrackEnergy[i] ;
    Eta[i]    = fTrackEta[i];
    Phi[i]    = fTrackPhi[i];
    PDG[i]    = fTrackPDG[i];
  } 
}


//____________________________________________________________________________

AliEMCALParton::~AliEMCALParton()
{
  // Destructor
  
}
