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

// MUON class for MonteCarlo Hits, inherited from AliHit for the 
// In addition to the ALiHit data member fX, fY, fZ and fTrack, AliMUONHit contains some info about the particle crossing the chamber:
// Impulsion: fPtot, fPx, fPy and fPz
// Reference position at the center of the chamber (wire plane) fXref, fYref and fZref
// Cumulated path along the active volume fTlength for spliting of hits for very inclined tracks 
// Energy loss of the particle inside the gas active volume.
// Incident fTheta and fPhi angle with respect of the wire plane of the chamber.
//

#include "AliMUONHit.h"
#include "AliMpDEManager.h"

#include "AliLog.h"

#include <Riostream.h>
#include <TMath.h>
#include <TString.h>

/// \cond CLASSIMP
ClassImp(AliMUONHit)
/// \endcond
 
//___________________________________________
AliMUONHit::AliMUONHit()
  : AliHit(), 
    fDetElemId(0),
    fParticle(0),
    fTheta(0),
    fPhi(0),
    fTlength(0),
    fEloss(0),
    fAge(0),
    fPHfirst(0),
    fPHlast(0),
    fPTot(0),
    fPx(0),
    fPy(0),
    fPz(0),
    fXref(0),
    fYref(0),
    fZref(0)
{
/// Default constructor
}

//___________________________________________
AliMUONHit::AliMUONHit(Int_t shunt, Int_t track, Int_t *vol, Float_t *hits)
  : AliHit(shunt, track),
    fDetElemId(vol[0]),
    fParticle(hits[0]),
    fTheta(hits[4]),
    fPhi(hits[5]),
    fTlength(hits[6]),
    fEloss(hits[7]),
    fAge(hits[14]),
    fPHfirst((Int_t)hits[8]),
    fPHlast((Int_t)hits[9]),
    fPTot(hits[10]),
    fPx(hits[11]),
    fPy(hits[12]),
    fPz(hits[13]),
    fXref(0),
    fYref(0),
    fZref(0)
{
/// Constructor
   
    fX         = hits[1];
    fY         = hits[2];
    fZ         = hits[3];
}

//___________________________________________
AliMUONHit::AliMUONHit(Int_t shunt, Int_t track, Int_t detElemId, Int_t idpart, 
		       Float_t X, Float_t Y, Float_t Z, Float_t tof, Float_t momentum, 
		       Float_t theta, Float_t phi, Float_t length, Float_t destep)
  : AliHit(shunt, track),
    fDetElemId(detElemId),
    fParticle(idpart),
    fTheta(theta),
    fPhi(phi),
    fTlength(length),
    fEloss(destep),
    fAge(tof),
    fPHfirst(0),
    fPHlast(0),
    fPTot(momentum),
    fPx(momentum * TMath::Sin(theta) * TMath::Cos(phi)),
    fPy(momentum * TMath::Sin(theta) * TMath::Sin(phi)),
    fPz(momentum * TMath::Cos(theta)),
    fXref(0),
    fYref(0),
    fZref(0)
{
/// Constructor
    fX         = X;
    fY         = Y;
    fZ         = Z;
}

//-----------------------------------------------------------------------------------------------
AliMUONHit::AliMUONHit(Int_t shunt, Int_t track, Int_t detElemId, Int_t idpart, 
		       Float_t X, Float_t Y, Float_t Z, Float_t tof, Float_t momentum, 
		       Float_t theta, Float_t phi, Float_t length, Float_t destep,
		       Float_t Xref,Float_t Yref,Float_t Zref)
  : AliHit(shunt, track),
    fDetElemId(detElemId),
    fParticle(idpart),
    fTheta(theta),
    fPhi(phi),
    fTlength(length),
    fEloss(destep),
    fAge(tof),
    fPHfirst(0),
    fPHlast(0),
    fPTot(momentum),
    fPx(momentum * TMath::Sin(theta) * TMath::Cos(phi)),
    fPy(momentum * TMath::Sin(theta) * TMath::Sin(phi)),
    fPz(momentum * TMath::Cos(theta)),
    fXref(Xref),
    fYref(Yref),
    fZref(Zref)
{
/// Constructor

    fX         = X;
    fY         = Y;
    fZ         = Z;
}

//-----------------------------------------------------------------------------------------------
AliMUONHit::~AliMUONHit()
{
/// Dectructor
}

//-----------------------------------------------------------------------------------------------
Int_t  AliMUONHit::Chamber()  const
{  
/// Return chamber ID

  return AliMpDEManager::GetChamberId(fDetElemId) + 1;  
}

//-----------------------------------------------------------------------------------------------
void AliMUONHit::Print(Option_t* opt) const
{
/// Printing hit information 
/// "full" option for printing all the information about the hit

  TString sopt(opt);
  sopt.ToUpper();
 
  if ( sopt.Contains("FULL") ) { 
    cout <<"<AliMUONHit>: Geant track="   << setw(4)  << Track() <<
      ", DetEle="        << setw(4)  << DetElemId() <<  
      ", (x,y,z)=(" << setw(7) << setprecision(5) << X() << "," << setw(7) << setprecision(5) << Y() <<  "," << setw(7) << setprecision(5) << Z() << 
      " )cm, Delta E=" << setw(8) << setprecision(3) << Eloss() << " GeV" << endl;
  }
  else {
    cout << "<AliMUONHit>: DetEle="        << setw(4)  << DetElemId() << 
      ", (x,y,z)=(" << setw(7) << setprecision(5) << X() << "," << setw(7) << setprecision(5) << Y() <<  "," << setw(7) << setprecision(5) << Z() << 
      " ) cm" <<endl;
  }
    
}
