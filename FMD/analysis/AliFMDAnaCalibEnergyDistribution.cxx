
#include "AliFMDAnaCalibEnergyDistribution.h"

ClassImp(AliFMDAnaCalibEnergyDistribution)

AliFMDAnaCalibEnergyDistribution::AliFMDAnaCalibEnergyDistribution() : TObject(),
  fArray(), fIsInit(kFALSE){
  
  
  
}
//____________________________________________________________________
void AliFMDAnaCalibEnergyDistribution::Init() {
  
  fArray.SetOwner();
  for(Int_t det = 1; det<=3;det++) {
    TObjArray* detArray = new TObjArray();
    fArray.AddAtAndExpand(detArray,det);
  }
  fIsInit = kTRUE;
}


//____________________________________________________________________
TH1F* AliFMDAnaCalibEnergyDistribution::GetEnergyDistribution(Int_t det, Char_t ring) {

  Int_t ringNumber     = (ring == 'I' ? 0 : 1);
  TObjArray* detArray  = (TObjArray*)fArray.At(det);
  TH1F* hEdist         = (TH1F*)detArray->At(ringNumber);    
  
  return hEdist;
}

//____________________________________________________________________
void AliFMDAnaCalibEnergyDistribution::SetEnergyDistribution(Int_t det, Char_t ring, TH1F* edist) {
  
  if(!fIsInit)
    Init();
  Int_t ringNumber     = (ring == 'I' ? 0 : 1);
  TObjArray* detArray  = (TObjArray*)fArray.At(det);
  detArray->AddAtAndExpand(edist,ringNumber);
  

}

//____________________________________________________________________
//
// EOF
//
