#ifndef ALIRSNCUTDAUGHTERD0_H
#define ALIRSNCUTDAUGHTERD0_H

//
// Cuts for selecting good pion candidates for D0 analysis
// with the data samples from pp/PbPb/pPb runs in 2010/2011/2012/2013.
// Applies track quality selection plus PID selection,
// with different tolerance ranges depending on the momentum.
//

#include "AliVTrack.h"
#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliRsnCut.h"
#include "AliRsnCutTrackQuality.h"

class AliRsnCutDaughterD0 : public AliRsnCut {

 public:

  AliRsnCutDaughterD0(const char *name = "", AliPID::EParticleType pid = AliPID::kPion);
  virtual ~AliRsnCutDaughterD0() { }

  void           SetNoPID(Bool_t yn = kTRUE)                  {fNoPID = yn;}
  //void           SetIsCheckOnMother(Bool_t yn = kTRUE)        {fIsCheckOnMother = yn;}
  void           SetPtDependentPIDCut(Bool_t yn = kTRUE)      {fPtDepPIDCut = yn;}
  void           SetPID(AliPID::EParticleType type)           {fPID = type;}
   
  void           SetTPCPionPIDCut(Double_t value)             {fPionTPCPIDCut = value;}
  void           SetTPCKaonPIDCut(Double_t value)             {fKaonTPCPIDCut = value;}
   
  void           SetTOFPionPIDCut(Double_t value)             {fPionTOFPIDCut = value;}
  void           SetTOFKaonPIDCut(Double_t value)             {fKaonTOFPIDCut = value;}
   
  AliRsnCutTrackQuality *CutQuality()                       {return &fCutQuality;}
  Bool_t                 MatchTOF(const AliVTrack *vtrack);
  Bool_t                 MatchTPC(const AliVTrack *vtrack);
  virtual Bool_t         IsSelected(TObject *obj);
   
 private:
  AliRsnCutDaughterD0(const AliRsnCutDaughterD0 &copy); // Not implemented
  AliRsnCutDaughterD0 &operator=(const AliRsnCutDaughterD0 &copy); // Not implemented
  Bool_t                fNoPID;            // flag to switch off PID check
  //Bool_t                fIsCheckOnMother;  // flag to switch off tracks check
  AliPID::EParticleType fPID;              // PID for track
  AliRsnCutTrackQuality fCutQuality;       // track quality cut 

 protected:
  Double_t         fPionTPCPIDCut;                // TPC nsigmas for pions
  Double_t         fKaonTPCPIDCut;                // TPC nsigmas for kaons
  Double_t         fPionTOFPIDCut;                // TOF nsigmas for pions
  Double_t         fKaonTOFPIDCut;                // TOF nsigmas for kaons
  Bool_t           fPtDepPIDCut;                  // flag for setting a pt dependent or independent PID cut

  ClassDef(AliRsnCutDaughterD0, 3)          // cut definitions for D0
    };

//__________________________________________________________________________________________________
inline Bool_t AliRsnCutDaughterD0::MatchTOF(const AliVTrack *vtrack)
{
  //
  // Checks if the track has matched the TOF detector
  //

  if (!vtrack) {
    AliWarning("NULL argument: impossible to check status");
    return kFALSE;
  }
  AliPIDResponse *fPidResponse = fEvent->GetPIDResponse();
  if (!fPidResponse) {
    AliFatal("NULL PID response");
    return kFALSE;
  }
  
  AliPIDResponse::EDetPidStatus status = fPidResponse->CheckPIDStatus(AliPIDResponse::kTOF,vtrack);
  if (status != AliPIDResponse::kDetPidOk) return kFALSE; 
  Float_t probMis = fPidResponse->GetTOFMismatchProbability(vtrack);
  if (probMis > 0.01) return kFALSE;
  if ((vtrack->GetStatus()&AliESDtrack::kTOFpid )==0 && vtrack->GetStatus()&AliESDtrack::kITSrefit )   return kFALSE;
  return kTRUE;
  

  //if (!(vtrack->GetStatus() & AliESDtrack::kTOFout)) return kFALSE;
  //if (!(vtrack->GetStatus() & AliESDtrack::kTIME  )) return kFALSE;

  return kTRUE;
}

//___________________________________________________________________________________________________
inline Bool_t AliRsnCutDaughterD0::MatchTPC(const AliVTrack *vtrack) 
{
  // Check if the track is good for TPC PID
  
  
  if (!vtrack) {
    AliWarning("NULL argument: impossible to check status");
    return kFALSE;
  }
  AliPIDResponse *fPidResponse = fEvent->GetPIDResponse();
  if (!fPidResponse) {
    AliFatal("NULL PID response");
    return kFALSE;
  }
  
  AliPIDResponse::EDetPidStatus status = fPidResponse->CheckPIDStatus(AliPIDResponse::kTPC,vtrack);
  if (status != AliPIDResponse::kDetPidOk) return kFALSE;
  /* UInt_t nclsTPCPID = vtrack->GetTPCsignalN(); */
  /* if(nclsTPCPID<0) return kFALSE; */
  return kTRUE;
}

#endif
