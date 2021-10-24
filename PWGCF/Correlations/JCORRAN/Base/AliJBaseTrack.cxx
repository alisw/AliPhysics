/**************************************************************************
 * Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
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

// Comment describing what this class does needed!

// $Id: AliJBaseTrack.cxx,v 1.5 2008/05/08 15:19:52 djkim Exp $
////////////////////////////////////////////////////
/*!
  \file AliJBaseTrack.cxx
  \brief
  \author J. Rak, D.J.Kim, R.Diaz (University of Jyvaskyla)
  \email: djkim@jyu.fi
  \version $Revision: 1.5 $
  \date $Date: 2008/05/08 15:19:52 $
  */
////////////////////////////////////////////////////

#include "AliJBaseTrack.h"
#include <TString.h>

//______________________________________________________________________________
AliJBaseTrack::AliJBaseTrack():
    fID(-1),
    fLabel(-9999), 
    fParticleType(-1), 
    fCharge(0), 
    fStatus(0), 
    fFlags(0),
    fTriggID(-1),
    fAssocID(-1),
    fTracEff(1.0),
    fMCIndex(-9999),
    fWeight(1.0)
{
  // constructor
}

//_____________________________________________________________
AliJBaseTrack::AliJBaseTrack(float px,float py, float pz, float e, Int_t id, Short_t ptype, Char_t charge):
    TLorentzVector( px, py, pz, e ), 
    fID(id),
    fLabel(-9999),
    fParticleType(ptype), 
    fCharge(charge), 
    fStatus(0), 
    fFlags(0),
    fTriggID(-1),
    fAssocID(-1),
    fTracEff(1.0),
    fMCIndex(-9999),
    fWeight(1.0)
{
  // constructor
}

//_____________________________________________________________
AliJBaseTrack::AliJBaseTrack(const AliJBaseTrack& a):
    TLorentzVector  (a), 
    fID             (a.fID),
    fLabel          (a.fLabel),
    fParticleType   ( a.fParticleType ), 
    fCharge         ( a.fCharge ), 
    fStatus         ( a.fStatus ), 
    fFlags       ( a.fFlags ),
    fTriggID( a.fTriggID ),
    fAssocID( a.fAssocID ),
    fTracEff( a.fTracEff ),
    fMCIndex( a.fMCIndex ),
    fWeight ( a.fWeight  )
{
  //copy constructor
}

//_____________________________________________________________
AliJBaseTrack::AliJBaseTrack(const TLorentzVector& a):
    TLorentzVector  (a),
    fID             ( -1 ),
    fLabel          ( -9999 ),
    fParticleType   ( -1 ),
    fCharge         ( 0 ),
    fStatus         ( 0 ),
    fFlags       ( 0 ),
    fTriggID(-1),
    fAssocID(-1),
    fTracEff(-1),
    fMCIndex(-9999),
    fWeight(1.0)
{
  //copy constructor
}
//_____________________________________________________________
AliJBaseTrack& AliJBaseTrack::operator=(const AliJBaseTrack& trk){
  //operator =  
  if(this != &trk){
    TLorentzVector::operator=(trk);
    fID           = trk.fID;
    fLabel        = trk.fLabel;
    fParticleType = trk.fParticleType;
    fCharge       = trk.fCharge;
    fStatus       = trk.fStatus;
    fFlags     = trk.fFlags;
    fTriggID   = trk.fTriggID;
    fAssocID   = trk.fAssocID;
    fTracEff   = trk.fTracEff;
    fMCIndex   = trk.fMCIndex;
    fWeight       = trk.fWeight;
  }
  return *this;
}

//_____________________________________________________________
void AliJBaseTrack::Print(Option_t *option) const{
  //object print out
  JUNUSED(option);
  std::cout<<Form("(ID,Type,Charge,Flags)=(%d, %d, %d, %d)" , 
                  fID, fParticleType, fCharge,  fFlags );
  TLorentzVector::Print();
  cout<<"ID ="<<fID <<endl; 
  cout<<"fLabel="<<fLabel <<endl; 
  cout<<"fParticleType="<<fParticleType <<endl; 
  cout<<"fCharge="<<fCharge <<endl; 
  cout<<"fStatus="<<fStatus <<endl; 
  cout<<"fFlags="<<fFlags <<endl; 
  cout<<"fTriggID="<<fTriggID <<endl; 
  cout<<"fAssocID="<<fAssocID <<endl; 
  cout<<"fTracEff="<<fTracEff <<endl; 
  cout<<"fWeight="<<fWeight <<endl;
}

ClassImp(AliJBaseTrack)

