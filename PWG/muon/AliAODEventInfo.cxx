
/* $Id$ */

// AliAODEventInfo: a class for AODs for the MUON Arm of the ALICE Experiment
// Author: P. Cortese, Universita' del Piemonte Orientale in Alessandria and
// INFN of Torino - Italy
//
// This class provides additional information about the AliAODEvent, some
// of this is specific to MUON. The information stored in this class will 
// follow the evolution of the framework (i.e. some data members  may be 
// moved into the header in the future).
// Currently it stores the beam energy and the information needed to decode
// the trigger pattern (which are the bits corresponding to the different
// muon triggers).
//

#include "AliAODEventInfo.h"
#include "AliAODHeader.h"
#define AliAODEventInfo_CXX

// ************************************************************************** //
ClassImp(AliAODEventInfo)

//______________________________________________________________________________
AliAODEventInfo::AliAODEventInfo():fBeamEnergy(0),
	fMuonSingleLPtL0(0),fMuonSingleHPtL0(0),fMuonLikeLPtL0(0),
	fMuonLikeHPtL0(0),fMuonUnlikeLPtL0(0),fMuonUnlikeHPtL0(0),
	fEv(0),fEi(this),fHe(0),fTr(0),fDi(0)
{
  // Default constructor
}

//______________________________________________________________________________
AliAODEventInfo::~AliAODEventInfo()
{
  // Default destructor
}

//______________________________________________________________________________
Bool_t AliAODEventInfo::IsHeaderAccessible(const Char_t *msg) const {
  // Tests if the header pointer is set
  if(fHe!=0){
    return 1;
  }else{
    if(msg)printf("Error! Header is not accessible: %s",msg);
    return 0;
  }
}

//______________________________________________________________________________
void AliAODEventInfo::SelectTriggerBits(UChar_t muonSingleLPtL0,
    UChar_t muonSingleHPtL0,UChar_t muonLikeLPtL0, UChar_t muonLikeHPtL0,
    UChar_t muonUnlikeLPtL0, UChar_t muonUnlikeHPtL0){
  // Define which bit in the trigger pattern corresponds to the given trigger condition 
  fMuonSingleLPtL0=muonSingleLPtL0;
  fMuonSingleHPtL0=muonSingleHPtL0;
  fMuonLikeLPtL0=muonLikeLPtL0;  
  fMuonLikeHPtL0=muonLikeHPtL0;  
  fMuonUnlikeLPtL0=muonUnlikeLPtL0;
  fMuonUnlikeHPtL0=muonUnlikeHPtL0;
}

//______________________________________________________________________________
Bool_t AliAODEventInfo::MuonSingleLPtL0() const {
  // Test if the event fired MUON_Single_LPt_L0
  if(IsHeaderAccessible("MuonSingleLPtL0"))
    return ((((AliAODHeader*)fHe.GetObject())->GetTriggerMask())>>fMuonSingleLPtL0)&0x1;
  else
    return 0;
}

//__________________________________________
Bool_t AliAODEventInfo::MuonSingleHPtL0() const {
  // Test if the event fired MUON_Single_HPt_L0
  if(IsHeaderAccessible("MuonSingleHPtL0"))
    return ((((AliAODHeader*)fHe.GetObject())->GetTriggerMask())>>fMuonSingleHPtL0)&0x1;
  else
    return 0;
}

//__________________________________________
Bool_t AliAODEventInfo::MuonLikeLPtL0() const {
  // Test if the event fired MUON_Like_LPt_L0
  if(IsHeaderAccessible("MuonLikeLPtL0"))
    return ((((AliAODHeader*)fHe.GetObject())->GetTriggerMask())>>fMuonLikeLPtL0)&0x1;
  else
    return 0;
}

//__________________________________________
Bool_t AliAODEventInfo::MuonLikeHPtL0() const {
  // Test if the event fired MUON_Like_HPt_L0
  if(IsHeaderAccessible("MuonLikeHPtL0"))
    return ((((AliAODHeader*)fHe.GetObject())->GetTriggerMask())>>fMuonLikeHPtL0)&0x1;
  else
    return 0;
}

//__________________________________________
Bool_t AliAODEventInfo::MuonUnlikeLPtL0() const {
  // Test if the event fired MUON_Unlike_LPt_L0
  if(IsHeaderAccessible("MuonUnlikeLPtL0"))
    return ((((AliAODHeader*)fHe.GetObject())->GetTriggerMask())>>fMuonUnlikeLPtL0)&0x1;
  else
    return 0;
}

//__________________________________________
Bool_t AliAODEventInfo::MuonUnlikeHPtL0() const {
  // Test if the event fired MUON_Unlike_HPt_L0
  if(IsHeaderAccessible("MuonUnlikeHPtL0"))
    return ((((AliAODHeader*)fHe.GetObject())->GetTriggerMask())>>fMuonUnlikeHPtL0)&0x1;
  else
    return 0;
}
