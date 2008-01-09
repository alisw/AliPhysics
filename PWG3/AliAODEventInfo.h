#include "AliAODHeader.h"
#include "TClonesArray.h"

#ifndef AliAODEventInfo_H
#define AliAODEventInfo_H

#include "AliAODEvent.h"


/* AliAODEventInfo: a class for AODs for the MUON Arm of the ALICE Experiment
 * Author: P. Cortese, Universita' del Piemonte Orientale in Alessandria and
 * INFN of Torino - Italy
 */

/* 2007/11/06 v0.00 Initial version */
/* 2007/12/18 v0.01 More compact information for what regards trigger info */

class AliAODEventInfo : public TNamed {
public:
  AliAODEventInfo();
  ~AliAODEventInfo();

  // Missing in AliAODHeader and added here
  Double_t fBeamEnergy; // Add beam energy not present in AliAODHeader
  UChar_t fMUON_Single_LPt_L0; // Decode trigger info
  UChar_t fMUON_Single_HPt_L0; // Decode trigger info
  UChar_t fMUON_Like_LPt_L0;   // Decode trigger info
  UChar_t fMUON_Like_HPt_L0;   // Decode trigger info
  UChar_t fMUON_Unlike_LPt_L0; // Decode trigger info
  UChar_t fMUON_Unlike_HPt_L0; // Decode trigger info

  void SetBeamEnergy(Double_t BeamEnergy){ fBeamEnergy=BeamEnergy; };
  Double_t EBeam(){ return fBeamEnergy; };
  Double_t SqrtS(){ return 2*fBeamEnergy; };

  // Data members and functions that provide automatic access
  TRef ev; // Event
  TRef ei; // EventInfo
  TRef he; // Header
  TRef tr; // Tracks
  TRef di; // Dimuons
  Bool_t IsHeaderAccessible(Char_t *msg=0);

  // -- Trigger information (to be updated in future)
  Bool_t MUON_Single_LPt_L0();
  Bool_t MUON_Single_HPt_L0();
  Bool_t MUON_Like_LPt_L0();
  Bool_t MUON_Like_HPt_L0();
  Bool_t MUON_Unlike_LPt_L0();
  Bool_t MUON_Unlike_HPt_L0();

  Int_t         GetNDimuons()           const { return (di!=0) ? ((TClonesArray*)di.GetObject())->GetSize() : 0;}
  Int_t         NDimu()           { return GetNDimuons();}

  ClassDef(AliAODEventInfo,1)  // Additional header for MUON arm
};

#endif
