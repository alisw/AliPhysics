#include "AliESDMuonTrack.h"

AliESDMuonTrack::AliESDMuonTrack (const AliESDMuonTrack& MUONTrack):TObject(MUONTrack)
{
  fInverseBendingMomentum = MUONTrack.fInverseBendingMomentum; 
  fThetaX                 = MUONTrack.fThetaX;           
  fThetaY                 = MUONTrack.fThetaY ;           
  fZ                      = MUONTrack.fZ;                
  fBendingCoor            = MUONTrack.fBendingCoor;      
  fNonBendingCoor         = MUONTrack.fNonBendingCoor;   
  fChi2                   = MUONTrack.fChi2;             
  fNHit                   = MUONTrack.fNHit ; 

  fMatchTrigger           = MUONTrack.fMatchTrigger;  
  fChi2MatchTrigger       = MUONTrack.fChi2MatchTrigger; 
}

AliESDMuonTrack& AliESDMuonTrack::operator=(const AliESDMuonTrack& MUONTrack)
{
  if (this == &MUONTrack)
    return *this;

  fInverseBendingMomentum = MUONTrack.fInverseBendingMomentum; 
  fThetaX                 = MUONTrack.fThetaX;           
  fThetaY                 = MUONTrack.fThetaY ;           
  fZ                      = MUONTrack.fZ;                
  fBendingCoor            = MUONTrack.fBendingCoor;      
  fNonBendingCoor         = MUONTrack.fNonBendingCoor;   
  fChi2                   = MUONTrack.fChi2;             
  fNHit                   = MUONTrack.fNHit ; 

  fMatchTrigger           = MUONTrack.fMatchTrigger;  
  fChi2MatchTrigger       = MUONTrack.fChi2MatchTrigger; 
 
  return *this;
}

ClassImp(AliESDMuonTrack)

