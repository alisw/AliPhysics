
#include "AliFMDAnaCalibEnergyDistribution.h"
#include "TAxis.h"
#include <AliLog.h>
#include <iostream>
ClassImp(AliFMDAnaCalibEnergyDistribution)

AliFMDAnaCalibEnergyDistribution::AliFMDAnaCalibEnergyDistribution() : TObject(),
  fArray(), 
  fIsInit(kFALSE),
  fNetaBins(0),
  fEtaMax(0),
  fEtaMin(0){
   
  
}
//____________________________________________________________________
void AliFMDAnaCalibEnergyDistribution::Init() {
  
  if(fNetaBins == 0)
    AliFatal("Set Eta bins before doing Init or anything else");
  
  fArray.SetOwner();
  
  for(Int_t i = 0; i<=fNetaBins+1; i++) {
    TObjArray* etaArray = new TObjArray();
    fArray.AddAtAndExpand(etaArray,i);
    for(Int_t det = 1; det<=3;det++) {
      TObjArray* detArray = new TObjArray();
      etaArray->AddAtAndExpand(detArray,det);
      
      
    }
  
  }
  fIsInit = kTRUE;
}


//____________________________________________________________________
TH1F* AliFMDAnaCalibEnergyDistribution::GetEnergyDistribution(Int_t det, Char_t ring, Float_t eta) {
  
  TAxis testaxis(fNetaBins,fEtaMin,fEtaMax);
  Int_t binnumber = testaxis.FindBin(eta);
  
  Int_t ringNumber     = (ring == 'I' ? 0 : 1);
  TObjArray* etaArray  = (TObjArray*)fArray.At(binnumber); 
  TObjArray* detArray  = (TObjArray*)etaArray->At(det); 
  TH1F* hEdist         = (TH1F*)detArray->At(ringNumber);    
  
  return hEdist;
}

//____________________________________________________________________
void AliFMDAnaCalibEnergyDistribution::SetEnergyDistribution(Int_t det, Char_t ring, Float_t eta, TH1F* edist ) {
  
  if(!fIsInit)
    Init();
  
  TAxis testaxis(fNetaBins,fEtaMin,fEtaMax);
  Int_t binnumber = testaxis.FindBin(eta);
  
  Int_t ringNumber     = (ring == 'I' ? 0 : 1);
  TObjArray* etaArray  = (TObjArray*)fArray.At(binnumber);
  TObjArray* detArray  = (TObjArray*)etaArray->At(det);
  
  detArray->AddAtAndExpand(edist,ringNumber);
  

}

//____________________________________________________________________
//
// EOF
//
