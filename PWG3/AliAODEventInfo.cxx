/* AliAODEventInfo: a class for AODs for the MUON Arm of the ALICE Experiment
 * Author: P. Cortese, Universita' del Piemonte Orientale in Alessandria and
 * INFN of Torino - Italy
 */

#include "AliAODEventInfo.h"
#define AliAODEventInfo_CXX

// ************************************************************************** //
ClassImp(AliAODEventInfo)

//______________________________________________________________________________
AliAODEventInfo::AliAODEventInfo():fBeamEnergy(0),
	fMUON_Single_LPt_L0(0),fMUON_Single_HPt_L0(0),fMUON_Like_LPt_L0(0),
	fMUON_Like_HPt_L0(0),fMUON_Unlike_LPt_L0(0),fMUON_Unlike_HPt_L0(0),
	ev(0),ei(this),he(0),tr(0),di(0)
{
}

//______________________________________________________________________________
AliAODEventInfo::~AliAODEventInfo()
{
}

//______________________________________________________________________________
Bool_t AliAODEventInfo::IsHeaderAccessible(Char_t *msg){
  if(he!=0){
    return 1;
  }else{
    if(msg)printf("Error! Header is not accessible: %s",msg);
    return 0;
  }
}

//______________________________________________________________________________
Bool_t AliAODEventInfo::MUON_Single_LPt_L0(){
  if(IsHeaderAccessible("MUON_Single_LPt_L0"))
    return ((((AliAODHeader*)he.GetObject())->GetTriggerMask())>>fMUON_Single_LPt_L0)&0x1;
  else
    return 0;
}

//______________________________________________________________________________
Bool_t AliAODEventInfo::MUON_Single_HPt_L0(){
  if(IsHeaderAccessible("MUON_Single_HPt_L0"))
    return ((((AliAODHeader*)he.GetObject())->GetTriggerMask())>>fMUON_Single_HPt_L0)&0x1;
  else
    return 0;
}

//______________________________________________________________________________
Bool_t AliAODEventInfo::MUON_Like_LPt_L0(){
  if(IsHeaderAccessible("MUON_Like_LPt_L0"))
    return ((((AliAODHeader*)he.GetObject())->GetTriggerMask())>>fMUON_Like_LPt_L0)&0x1;
  else
    return 0;
}

//______________________________________________________________________________
Bool_t AliAODEventInfo::MUON_Like_HPt_L0(){
  if(IsHeaderAccessible("MUON_Like_HPt_L0"))
    return ((((AliAODHeader*)he.GetObject())->GetTriggerMask())>>fMUON_Like_HPt_L0)&0x1;
  else
    return 0;
}

//______________________________________________________________________________
Bool_t AliAODEventInfo::MUON_Unlike_LPt_L0(){
  if(IsHeaderAccessible("MUON_Unlike_LPt_L0"))
    return ((((AliAODHeader*)he.GetObject())->GetTriggerMask())>>fMUON_Unlike_LPt_L0)&0x1;
  else
    return 0;
}

//______________________________________________________________________________
Bool_t AliAODEventInfo::MUON_Unlike_HPt_L0(){
  if(IsHeaderAccessible("MUON_Unlike_HPt_L0"))
    return ((((AliAODHeader*)he.GetObject())->GetTriggerMask())>>fMUON_Unlike_HPt_L0)&0x1;
  else
    return 0;
}
