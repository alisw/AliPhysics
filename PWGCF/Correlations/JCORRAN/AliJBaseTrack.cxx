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

#include <iostream>
#include <TLorentzVector.h>
#include "AliJBaseTrack.h"

//______________________________________________________________________________
AliJBaseTrack::AliJBaseTrack():
    fID(-1),
    fLabel(-1), 
    fParticleType(-1), 
    fCharge(0), 
    fStatus(0), 
    fFlags(0)
{
  // constructor
}

//_____________________________________________________________
AliJBaseTrack::AliJBaseTrack(float px,float py, float pz, float e, Int_t id, Short_t ptype, Char_t charge):
    TLorentzVector( px, py, pz, e ), 
    fID(id),
    fLabel(-1), 
    fParticleType(ptype), 
    fCharge(charge), 
    fStatus(0), 
    fFlags(0)
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
    fFlags       ( a.fFlags )
{
  //copy constructor
}

//_____________________________________________________________
AliJBaseTrack::AliJBaseTrack(const TLorentzVector& a):
    TLorentzVector  (a),
    fID             ( -1 ),
    fLabel          ( -1 ), 
    fParticleType   ( -1 ),
    fCharge         ( 0 ),
    fStatus         ( 0 ),
    fFlags       ( 0 )
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
  }
  return *this;
}

//_____________________________________________________________
void AliJBaseTrack::Print(Option_t* option = "" ) const{
  //object print out
  std::cout<<Form("(ID,Type,Charge,Flags)=(%d, %d, %d, %d)" , 
                  fID, fParticleType, fCharge,  fFlags );
  TLorentzVector::Print(option);
}

ClassImp(AliJBaseTrack)

