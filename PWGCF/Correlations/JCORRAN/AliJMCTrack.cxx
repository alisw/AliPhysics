/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notifce   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// $Id: AliJMCTrack.cxx,v 1.2 2008/05/08 13:44:45 djkim Exp $

////////////////////////////////////////////////////
//
//  \file AliJMCTrack.cxx
//  \brief
//  \author J. Rak, D.J.Kim, R.Diaz (University of Jyvaskyla)
//  \email: djkim@jyu.fi
//  \version $Revision: 1.2 $
//  \date $Date: 2008/05/08 13:44:45 $
//
//  class which encapsulates track  monte carlo information
////////////////////////////////////////////////////

#include "AliJMCTrack.h"
#include <TMath.h>
#include <TROOT.h>

ClassImp(AliJMCTrack);

//______________________________________________________________________________
AliJMCTrack::AliJMCTrack() : 
    AliJBaseTrack(),
    fPdgCode(0),
    fVx(0),
    fVy(0),
    fVz(0)
{
  // default constructor
  fMother[0] = fMother[1] = 0;
  fDaughter[0] = fDaughter[1] = 0;


}

//______________________________________________________________________________
AliJMCTrack::AliJMCTrack(const AliJMCTrack& a):
    AliJBaseTrack(a),
    fPdgCode(a.fPdgCode),
    fVx(a.fVx),
    fVy(a.fVy),
    fVz(a.fVz)
{
  //copy constructor
  for(Int_t i=0;i<2;i++){
    fMother[i] = a.fMother[i];
    fDaughter[i] = a.fDaughter[i];
  }
}


//______________________________________________________________________________
void AliJMCTrack::SetPdgCode(Int_t icode) {
  // Set PDG code, charge,  recalculate E
  if( TMath::Abs(icode) > 32767 ) icode = 0; // Short_t
  fPdgCode=icode;
  SetVectM(Vect(), GetPDGData().Mass());
  SetCharge( TMath::Nint(GetPDGData().Charge()) ); // is this right?
}


//______________________________________________________________________________
Bool_t AliJMCTrack::IsHadron() const{
  // Check is hadron 
  int absID = TMath::Abs(GetPdgCode());
  if( absID >= 211 && absID<=533 ) return true; //meson
  if( absID >1000 && absID<6000 ) return true; // barion
  return false;
}

//______________________________________________________________________________
AliJMCTrack& AliJMCTrack::operator=(const AliJMCTrack& trk){
  //operator= 
  if(this != &trk){
    AliJBaseTrack::operator=(trk);
    fPdgCode   =  trk.fPdgCode;
    fLabel     =  trk.fLabel;
    for(Int_t i=0;i<2;i++){
      fMother[i]   = trk.fMother[i];
      fDaughter[i] = trk.fDaughter[i];
    }
    fVx  = trk.fVx;
    fVy  = trk.fVy;
    fVz  = trk.fVz;
  }
  return *this;
}

const TParticlePDG& AliJMCTrack::GetPDGData() const { 
  return *(TDatabasePDG::Instance()->GetParticle( fPdgCode )); 
}





