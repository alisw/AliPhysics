////////////////////////////////////////////////
//  Digits classes for set:ITS                //
////////////////////////////////////////////////


#include "AliITSdigit.h"

ClassImp(AliITSdigit)
//_____________________________________________________________________________

AliITSdigit::AliITSdigit(Int_t *digits) {
  //
  // Creates a real data digit object
  //
  fCoord1       = digits[0];
  fCoord2       = digits[1];
  fSignal       = digits[2];
}


ClassImp(AliITSdigitSPD)
  
  //_________________________________________________________________________
  AliITSdigitSPD::AliITSdigitSPD(Int_t *digits) {
  //
  // Creates a SPD digit object 
  //
  
  fCoord1        = digits[0];
  fCoord2        = digits[1];
  fSignal        = digits[2];
  
}

//_____________________________________________________________________________
AliITSdigitSPD::AliITSdigitSPD(Int_t *digits,Int_t *tracks,Int_t *hits) {  
  //
  // Creates a simulated SPD digit object 
  //
  
  fCoord1        = digits[0];
  fCoord2        = digits[1];
  fSignal        = digits[2];
  
  for(Int_t i=0; i<3; i++) {
    fTracks[i]    = tracks[i];
    fHits[i]      = hits[i];
  }
}


ClassImp(AliITSdigitSDD)
  //________________________________________________________________________
AliITSdigitSDD::AliITSdigitSDD(Float_t phys,Int_t *digits) {
  //
  // Creates a simulated SDD digit object to be updated
  //
  fCoord1       = digits[0];
  fCoord2       = digits[1];
  fSignal       = digits[2];
  fPhysics      = phys;
}

//_____________________________________________________________________________
AliITSdigitSDD::AliITSdigitSDD(Float_t phys,Int_t *digits,Int_t *tracks,Int_t *hits,Float_t *charges) {
  //
  // Creates a simulated SDD digit object 
  //
  fCoord1        = digits[0];
  fCoord2        = digits[1];
  fSignal        = digits[2];
  fPhysics       = phys;
  
  for(Int_t i=0; i<3; i++) {
    fTcharges[i]  = charges[i];
    fTracks[i]    = tracks[i];
    fHits[i]      = hits[i];
  }
}


ClassImp(AliITSTransientDigit)
  //_______________________________________________________________________
AliITSTransientDigit::AliITSTransientDigit(Float_t phys,Int_t *digits): 
    AliITSdigitSDD(phys,digits) {
  //
  // Creates a digit object in a list of digits to be updated
  //  
  fTrackList   = new TObjArray;  
}

//__________________________________________________________________________
AliITSTransientDigit::AliITSTransientDigit(const AliITSTransientDigit &source){
  //     Copy Constructor 
  if(&source == this) return;
  this->fTrackList = source.fTrackList;
  return;
}

//_________________________________________________________________________
AliITSTransientDigit& 
  AliITSTransientDigit::operator=(const AliITSTransientDigit &source) {
  //    Assignment operator
  if(&source == this) return *this;
  this->fTrackList = source.fTrackList;
  return *this;
}

ClassImp(AliITSdigitSSD)
  //__________________________________________________________________________
AliITSdigitSSD::AliITSdigitSSD(Int_t *digits) {
  //
  // Creates a real SSD digit object 
  //
  
  fCoord1        = digits[0];
  fCoord2        = digits[1];
  fSignal        = digits[2];
  
}

//_____________________________________________________________________________
AliITSdigitSSD::AliITSdigitSSD(Int_t *digits,Int_t *tracks,Int_t *hits) {
  //
  // Creates a simulated SSD digit object 
  // 
  
  fCoord1        = digits[0];
  fCoord2        = digits[1];
  fSignal        = digits[2];
  
  for(Int_t i=0; i<3; i++) {
    fTracks[i]    = tracks[i];
    fHits[i]      = hits[i];
  }
}

