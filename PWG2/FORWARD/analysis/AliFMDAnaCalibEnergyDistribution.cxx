
#include "AliFMDAnaCalibEnergyDistribution.h"
#include "TAxis.h"
#include <AliLog.h>
#include <iostream>
#include "TH2F.h"
#include "AliFMDAnaParameters.h"
ClassImp(AliFMDAnaCalibEnergyDistribution)

AliFMDAnaCalibEnergyDistribution::AliFMDAnaCalibEnergyDistribution() : TObject(),
  fArray(), 
  fEmptyArray(), 
  fRingArray(), 
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
  fEmptyArray.SetOwner();
  
  for(Int_t i = 0; i<fNetaBins; i++) {
    TObjArray* etaArray = new TObjArray();
    fArray.AddAtAndExpand(etaArray,i);
    for(Int_t det = 1; det<=3;det++) {
      TObjArray* detArray = new TObjArray();
      etaArray->AddAtAndExpand(detArray,det);
            
    }
  
  }
  

  
  for(Int_t det = 1; det<=3;det++) {
    TObjArray* detArray = new TObjArray();
    fEmptyArray.AddAtAndExpand(detArray,det);
  }
    
  for(Int_t det = 1; det<=3;det++) {
    TObjArray* detArray = new TObjArray();
    fRingArray.AddAtAndExpand(detArray,det);
  }
    
    
  fIsInit = kTRUE;
}


//____________________________________________________________________
TH1F* AliFMDAnaCalibEnergyDistribution::GetEnergyDistribution(Int_t det, Char_t ring, Float_t eta) {
  
  //TAxis testaxis(fNetaBins,fEtaMin,fEtaMax);
  //  Int_t binnumber = testaxis.FindBin(eta);
  AliFMDAnaParameters* pars = AliFMDAnaParameters::Instance();
    
  Int_t binnumber = pars->GetEtaBin(eta);
  
  Int_t ringNumber     = (ring == 'I' ? 0 : 1);
  TObjArray* etaArray  = (TObjArray*)fArray.At(binnumber); 
  TObjArray* detArray  = (TObjArray*)etaArray->At(det); 
  TH1F* hEdist         = (TH1F*)detArray->At(ringNumber);    
  
  return hEdist;
}
//____________________________________________________________________
TH1F* AliFMDAnaCalibEnergyDistribution::GetEmptyEnergyDistribution(Int_t det, Char_t ring) {
  
  Int_t ringNumber     = (ring == 'I' ? 0 : 1);
  
  TObjArray* detArray  = (TObjArray*)fEmptyArray.At(det); 
  TH1F* hEdist         = (TH1F*)detArray->At(ringNumber);    
  
  return hEdist;
}
//____________________________________________________________________
TH1F* AliFMDAnaCalibEnergyDistribution::GetRingEnergyDistribution(Int_t det, Char_t ring) {
  
  Int_t ringNumber     = (ring == 'I' ? 0 : 1);
  
  TObjArray* detArray  = (TObjArray*)fRingArray.At(det); 
  TH1F* hEdist         = (TH1F*)detArray->At(ringNumber);    
  
  return hEdist;
}
//____________________________________________________________________
void  AliFMDAnaCalibEnergyDistribution::SetEnergyDistribution(Int_t det, Char_t ring, Float_t eta, TH1F* edist) {
  AliFMDAnaParameters* pars = AliFMDAnaParameters::Instance();
  Int_t binnumber = pars->GetEtaBin(eta);
  //std::cout<<binnumber<<std::endl;
  SetEnergyDistribution(det, ring, binnumber, edist );
}
//____________________________________________________________________
void AliFMDAnaCalibEnergyDistribution::SetEnergyDistribution(Int_t det, Char_t ring, Int_t etabin, TH1F* edist ) {
  
  if(!fIsInit)
    Init();
  
  Int_t ringNumber     = (ring == 'I' ? 0 : 1);
  TObjArray* etaArray  = (TObjArray*)fArray.At(etabin);
  TObjArray* detArray  = (TObjArray*)etaArray->At(det);
  
  detArray->AddAtAndExpand(edist,ringNumber);
  

}

//____________________________________________________________________
void AliFMDAnaCalibEnergyDistribution::SetEmptyEnergyDistribution(Int_t det, Char_t ring, TH1F* edist ) {
  
  if(!fIsInit)
    Init();
    
  Int_t ringNumber     = (ring == 'I' ? 0 : 1);
  TObjArray* detArray  = (TObjArray*)fEmptyArray.At(det);
  
  detArray->AddAtAndExpand(edist,ringNumber);
  

}
//____________________________________________________________________
void AliFMDAnaCalibEnergyDistribution::SetRingEnergyDistribution(Int_t det, Char_t ring, TH1F* edist ) {
  
  if(!fIsInit)
    Init();
    
  Int_t ringNumber     = (ring == 'I' ? 0 : 1);
  TObjArray* detArray  = (TObjArray*)fRingArray.At(det);
  
  detArray->AddAtAndExpand(edist,ringNumber);
  

}
//____________________________________________________________________
//
// EOF
//
