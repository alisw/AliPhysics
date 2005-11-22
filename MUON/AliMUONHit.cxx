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

// MUON classe for MonteCarlo Hits, inherited from AliHit for the 
// In addition to the ALiHit data member fX, fY, fZ and fTrack, AliMUONHit contains some info about the particle crossing the chamber:
// Impulsion: fPtot, fPx, fPy and fPz
// Reference position at the center of the chamber (wire plane) fXref, fYref and fZref
// Cumulated path along the active volume fTlength for spliting of hits for very inclined tracks 
// Energy loss of the particle inside the gas active volume.
// Incident fTheta and fPhi angle with respect of the wire plane of the chamber.
//

#include <TMath.h>

#include "AliMUONHit.h"
#include "AliMUONGeometryDEIndexing.h"
#include "AliLog.h"

ClassImp(AliMUONHit)
 
//___________________________________________
AliMUONHit::AliMUONHit()
  : AliHit() 
{
// Default constructor
}

//___________________________________________
AliMUONHit::AliMUONHit(Int_t shunt, Int_t track, Int_t *vol, Float_t *hits):
	AliHit(shunt, track)
{
// Constructor
// TBR
   
    fIsDetElemId = kFALSE;
    fDetElemId = vol[0];
    fParticle  = hits[0];
    fX         = hits[1];
    fY         = hits[2];
    fZ         = hits[3];
    fTheta     = hits[4];
    fPhi       = hits[5];
    fTlength   = hits[6];
    fEloss     = hits[7];
    fPHfirst   = (Int_t) hits[8];
    fPHlast    = (Int_t) hits[9];
    fPTot      = hits[10];
    fPx        = hits[11];
    fPy        = hits[12];
    fPz        = hits[13];
    fAge       = hits[14];
    fXref      = 0.;
    fYref      = 0.;
    fZref      = 0.;
}

//___________________________________________
AliMUONHit::AliMUONHit(Int_t shunt, Int_t track, Int_t *vol, Float_t *hits, 
                       Bool_t /*isDetElemId*/) :
	AliHit(shunt, track)
{
// Constructor
   
    fIsDetElemId = kTRUE;
    fDetElemId = vol[0];
    fParticle  = hits[0];
    fX         = hits[1];
    fY         = hits[2];
    fZ         = hits[3];
    fTheta     = hits[4];
    fPhi       = hits[5];
    fTlength   = hits[6];
    fEloss     = hits[7];
    fPHfirst   = (Int_t) hits[8];
    fPHlast    = (Int_t) hits[9];
    fPTot      = hits[10];
    fPx        = hits[11];
    fPy        = hits[12];
    fPz        = hits[13];
    fAge       = hits[14];
    fXref      = 0.;
    fYref      = 0.;
    fZref      = 0.;
}

//___________________________________________
AliMUONHit::AliMUONHit(Int_t shunt, Int_t track, Int_t iChamber, Int_t idpart, 
		       Float_t X, Float_t Y, Float_t Z, Float_t tof, Float_t momentum, 
		       Float_t theta, Float_t phi, Float_t length, Float_t destep):
	AliHit(shunt, track)
{
// Constructor
// TBR

    fIsDetElemId = kFALSE;
    fDetElemId = iChamber;
    fParticle  = idpart;
    fX         = X;
    fY         = Y;
    fZ         = Z;
    fTheta     = theta;
    fPhi       = phi;
    fTlength   = length;
    fEloss     = destep;
    fPHfirst   = 0;
    fPHlast    = 0;
    fPTot      = momentum;
    fPx        = momentum * TMath::Sin(theta) * TMath::Cos(phi);
    fPy        = momentum * TMath::Sin(theta) * TMath::Sin(phi);
    fPz        = momentum * TMath::Cos(theta) ;
    fAge       = tof;
    fXref      = 0.;
    fYref      = 0.;
    fZref      = 0.;
}

//___________________________________________
AliMUONHit::AliMUONHit(Int_t shunt, Int_t track, Int_t detElemId, Int_t idpart, 
		       Float_t X, Float_t Y, Float_t Z, Float_t tof, Float_t momentum, 
		       Float_t theta, Float_t phi, Float_t length, Float_t destep,
		       Bool_t /*isDetElemId*/):
	AliHit(shunt, track)
{
// Constructor
    fIsDetElemId = kTRUE;
    fDetElemId = detElemId;
    fParticle  = idpart;
    fX         = X;
    fY         = Y;
    fZ         = Z;
    fTheta     = theta;
    fPhi       = phi;
    fTlength   = length;
    fEloss     = destep;
    fPHfirst   = 0;
    fPHlast    = 0;
    fPTot      = momentum;
    fPx        = momentum * TMath::Sin(theta) * TMath::Cos(phi);
    fPy        = momentum * TMath::Sin(theta) * TMath::Sin(phi);
    fPz        = momentum * TMath::Cos(theta) ;
    fAge       = tof;
    fXref      = 0.;
    fYref      = 0.;
    fZref      = 0.;
}

//-----------------------------------------------------------------------------------------------
AliMUONHit::AliMUONHit(Int_t shunt, Int_t track, Int_t iChamber, Int_t idpart, 
		       Float_t X, Float_t Y, Float_t Z, Float_t tof, Float_t momentum, 
		       Float_t theta, Float_t phi, Float_t length, Float_t destep,
		       Float_t Xref,Float_t Yref,Float_t Zref):
	AliHit(shunt, track)
{
// Constructor
// TBR

    fIsDetElemId = kFALSE;
    fDetElemId = iChamber;
    fParticle  = idpart;
    fX         = X;
    fY         = Y;
    fZ         = Z;
    fTheta     = theta;
    fPhi       = phi;
    fTlength   = length;
    fEloss     = destep;
    fPHfirst   = 0;
    fPHlast    = 0;
    fPTot      = momentum;
    fPx        = momentum * TMath::Sin(theta) * TMath::Cos(phi);
    fPy        = momentum * TMath::Sin(theta) * TMath::Sin(phi);
    fPz        = momentum * TMath::Cos(theta) ;
    fAge       = tof;
    fXref      = Xref;
    fYref      = Yref;
    fZref      = Zref;
}
//-----------------------------------------------------------------------------------------------
AliMUONHit::AliMUONHit(Int_t shunt, Int_t track, Int_t detElemId, Int_t idpart, 
		       Float_t X, Float_t Y, Float_t Z, Float_t tof, Float_t momentum, 
		       Float_t theta, Float_t phi, Float_t length, Float_t destep,
		       Float_t Xref,Float_t Yref,Float_t Zref,
		       Bool_t /*isDetElemId*/):
	AliHit(shunt, track)
{
// Constructor
    fIsDetElemId = kTRUE;
    fDetElemId = detElemId;
    fParticle  = idpart;
    fX         = X;
    fY         = Y;
    fZ         = Z;
    fTheta     = theta;
    fPhi       = phi;
    fTlength   = length;
    fEloss     = destep;
    fPHfirst   = 0;
    fPHlast    = 0;
    fPTot      = momentum;
    fPx        = momentum * TMath::Sin(theta) * TMath::Cos(phi);
    fPy        = momentum * TMath::Sin(theta) * TMath::Sin(phi);
    fPz        = momentum * TMath::Cos(theta) ;
    fAge       = tof;
    fXref      = Xref;
    fYref      = Yref;
    fZref      = Zref;
}

//-----------------------------------------------------------------------------------------------
Int_t AliMUONHit::DetElemId()const
{
// Return detection element ID

  if (!fIsDetElemId) {
    AliWarning("Detection element Id is not defined.");
    return 0;
  }  
  // end of TBR
  
  return fDetElemId;
}

//-----------------------------------------------------------------------------------------------
Int_t  AliMUONHit::Chamber()  const
{  
// Return chamber ID

  if (!fIsDetElemId) 
    return fDetElemId;
  else  
    return AliMUONGeometryDEIndexing::GetModuleId(fDetElemId)+1;  
}

