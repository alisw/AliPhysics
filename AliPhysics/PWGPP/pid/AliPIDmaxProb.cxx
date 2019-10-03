#include"AliPIDmaxProb.h"
#include"AliPIDResponse.h"
#include"AliAnalysisManager.h"
#include"AliInputEventHandler.h"

ClassImp(AliPIDmaxProb);

AliPIDmaxProb::AliPIDmaxProb(const char *name):
  AliPIDperfCut(name),
  fPIDCombined(NULL),
  fMaskPID(AliPIDResponse::kDetTPC|AliPIDResponse::kDetTOF),
  fTPCin(kFALSE),
  fTOFin(kFALSE)
{
  fPIDCombined=new AliPIDCombined;
  fPIDCombined->SetDefaultTPCPriors();
}
//--------------------------------------------------------------------
AliPIDmaxProb::AliPIDmaxProb():
  AliPIDperfCut("AliPIDmaxProb"),
  fPIDCombined(NULL),
  fMaskPID(AliPIDResponse::kDetTPC|AliPIDResponse::kDetTOF),
  fTPCin(kFALSE),
  fTOFin(kFALSE)
{
  fPIDCombined=new AliPIDCombined;
  fPIDCombined->SetDefaultTPCPriors();
}
//--------------------------------------------------------------------
Bool_t AliPIDmaxProb::IsSelected(AliVTrack *track,AliPID::EParticleType type) const{
  Int_t nspecies = AliPID::kSPECIES;
  if(type >= nspecies || type < 0) return 0; // type has not a proper value

  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  AliPIDResponse *PIDResponse=inputHandler->GetPIDResponse();

  fPIDCombined->SetDetectorMask(fMaskPID);

  Double_t prob[10] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

  UInt_t detUsed = fPIDCombined->ComputeProbabilities(track, PIDResponse, prob);

  if(fTPCin && !(detUsed & AliPIDResponse::kDetTPC)) return kFALSE;
  if(fTOFin && !(detUsed & AliPIDResponse::kDetTOF)) return kFALSE;

  Bool_t status = kTRUE;
  for(Int_t i=0;i<nspecies;i++){
    if(i != type && prob[type] <= prob[i]){
      status = kFALSE;
      break;
    }
  }

  return status;
}
