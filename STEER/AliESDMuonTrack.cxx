#include "AliESDMuonTrack.h"

AliESDMuonTrack::AliESDMuonTrack (const AliESDMuonTrack& MUONTrack):TObject(MUONTrack)
{
  fInverseBendingMomentum = MUONTrack.fInverseBendingMomentum; 
  fThetaX = MUONTrack.fThetaX;           
  fThetaY = MUONTrack.fThetaY ;           
  fZ = MUONTrack.fZ;                
  fBendingCoor = MUONTrack.fBendingCoor;      
  fNonBendingCoor = MUONTrack.fNonBendingCoor;   
  fChi2 = MUONTrack.fChi2;             
  fNHit= MUONTrack.fNHit ; 

  fX11  = MUONTrack.fX11;  
  fY11  = MUONTrack.fY11;
  fThetaX11  = MUONTrack.fThetaX11; 
  fThetaY11  = MUONTrack.fThetaY11;      
}

AliESDMuonTrack& AliESDMuonTrack::operator=(const AliESDMuonTrack& MUONTrack)
{
  if (this == &MUONTrack)
    return *this;

  fInverseBendingMomentum = MUONTrack.fInverseBendingMomentum; 
  fThetaX = MUONTrack.fThetaX;           
  fThetaY = MUONTrack.fThetaY ;           
  fZ = MUONTrack.fZ;                
  fBendingCoor = MUONTrack.fBendingCoor;      
  fNonBendingCoor = MUONTrack.fNonBendingCoor;   
  fChi2 = MUONTrack.fChi2;             
  fNHit= MUONTrack.fNHit ; 

  fX11  = MUONTrack.fX11;  
  fY11  = MUONTrack.fY11;
  fThetaX11  = MUONTrack.fThetaX11; 
  fThetaY11  = MUONTrack.fThetaY11;  

  return *this;
}

ClassImp(AliESDMuonTrack)

