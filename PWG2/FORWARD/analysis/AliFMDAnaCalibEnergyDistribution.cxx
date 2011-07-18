//
// Object to store energy distribution corrections as used in the FMD
// analysis.  
//
//
#include "AliFMDAnaCalibEnergyDistribution.h"
#include "TAxis.h"
#include <AliLog.h>
#include <iostream>
#include "TH2F.h"
#include <TBrowser.h>
#include "AliFMDAnaParameters.h"
ClassImp(AliFMDAnaCalibEnergyDistribution)
#if 0  
; // This is for Emacs 
#endif 

//____________________________________________________________________
AliFMDAnaCalibEnergyDistribution::AliFMDAnaCalibEnergyDistribution() 
: TObject(),
  fArray(), 
  fEmptyArray(), 
  fRingArray(), 
  fIsInit(kFALSE),
  fNetaBins(0),
  fEtaMax(0),
  fEtaMin(0)
{
  
  
}
//____________________________________________________________________
void AliFMDAnaCalibEnergyDistribution::Init() 
{
  //Init object
  if(fNetaBins == 0)
    AliFatal("Set Eta bins before doing Init or anything else");
  
  fArray.SetOwner();
  fArray.SetName("etaBins");
  fEmptyArray.SetOwner();
  fEmptyArray.SetName("empty");
  fRingArray.SetName("rings");

  for(Int_t i = 0; i<fNetaBins; i++) {
    TObjArray* etaArray = new TObjArray();
    etaArray->SetName(Form("etabin_%03d", i+1));
    fArray.AddAtAndExpand(etaArray,i);
    for(Int_t det = 1; det<=3;det++) {
      TObjArray* detArray = new TObjArray();
      detArray->SetName(Form("FMD%d", det));
      etaArray->AddAtAndExpand(detArray,det);
    }
  }
  

  
  for(Int_t det = 1; det<=3;det++) {
    TObjArray* detArray = new TObjArray();
    detArray->SetName(Form("FMD%d", det));
    fEmptyArray.AddAtAndExpand(detArray,det);
  }
    
  for(Int_t det = 1; det<=3;det++) {
    TObjArray* detArray = new TObjArray();
    detArray->SetName(Form("FMD%d", det));
    fRingArray.AddAtAndExpand(detArray,det);
  }
    
    
  fIsInit = kTRUE;
}


//____________________________________________________________________
TH1F* 
AliFMDAnaCalibEnergyDistribution::GetEnergyDistribution(Int_t det, 
							Char_t ring, 
							Float_t eta) {
 
  //Get Energy dist
  //TAxis testaxis(fNetaBins,fEtaMin,fEtaMax);
  //  Int_t binnumber = testaxis.FindBin(eta);
  AliFMDAnaParameters* pars = AliFMDAnaParameters::Instance();
    
  Int_t binnumber       = pars->GetEtaBin(eta);
  
  Int_t      ringNumber = (ring == 'I' ? 0 : 1);
  TObjArray* etaArray   = (TObjArray*)fArray.At(binnumber); 
  TObjArray* detArray   = (TObjArray*)etaArray->At(det); 
  TH1F*      hEdist     = (TH1F*)detArray->At(ringNumber);    
  
  return hEdist;
}
//____________________________________________________________________
TH1F* 
AliFMDAnaCalibEnergyDistribution::GetEmptyEnergyDistribution(Int_t det, 
							     Char_t ring) 
{
  //Get e dist of empty
  Int_t ringNumber     = (ring == 'I' ? 0 : 1);
  
  TObjArray* detArray  = (TObjArray*)fEmptyArray.At(det); 
  TH1F* hEdist         = (TH1F*)detArray->At(ringNumber);    
  
  return hEdist;
}
//____________________________________________________________________
TH1F* 
AliFMDAnaCalibEnergyDistribution::GetRingEnergyDistribution(Int_t det, 
							    Char_t ring) {
  // Get E dist of ring
  Int_t ringNumber     = (ring == 'I' ? 0 : 1);
  
  TObjArray* detArray  = (TObjArray*)fRingArray.At(det); 
  TH1F* hEdist         = (TH1F*)detArray->At(ringNumber);    
  
  return hEdist;
}
//____________________________________________________________________
void  
AliFMDAnaCalibEnergyDistribution::SetEnergyDistributionUser(Int_t det, 
							    Char_t ring, 
							    Float_t eta, 
							    TH1F* edist) 
{
  //Set E dist (user)
  AliFMDAnaParameters* pars = AliFMDAnaParameters::Instance();
  Int_t binnumber = pars->GetEtaBin(eta);
  //std::cout<<binnumber<<std::endl;
  SetEnergyDistribution(det, ring, binnumber, edist );
}
//____________________________________________________________________
void 
AliFMDAnaCalibEnergyDistribution::SetEnergyDistribution(Int_t  det, 
							Char_t ring, 
							Int_t  etabin, 
							TH1F*  edist) 
{  
  //Set E dist
  if(!fIsInit)  Init();
  
  Int_t ringNumber     = (ring == 'I' ? 0 : 1);
  TObjArray* etaArray  = (TObjArray*)fArray.At(etabin);
  TObjArray* detArray  = (TObjArray*)etaArray->At(det);
  
  detArray->AddAtAndExpand(edist,ringNumber);
}

//____________________________________________________________________
void 
AliFMDAnaCalibEnergyDistribution::SetEmptyEnergyDistribution(Int_t  det, 
							     Char_t ring, 
							     TH1F*  edist) 
{  
  //Set the empty dist
  if(!fIsInit)
    Init();
    
  Int_t ringNumber     = (ring == 'I' ? 0 : 1);
  TObjArray* detArray  = (TObjArray*)fEmptyArray.At(det);
  
  detArray->AddAtAndExpand(edist,ringNumber);
  

}
//____________________________________________________________________
void AliFMDAnaCalibEnergyDistribution::SetRingEnergyDistribution(Int_t  det, 
								 Char_t ring, 
								 TH1F*  edist) 
{
  // Set E dist of ring
  if(!fIsInit) Init();
    
  Int_t ringNumber     = (ring == 'I' ? 0 : 1);
  TObjArray* detArray  = (TObjArray*)fRingArray.At(det);
  
  detArray->AddAtAndExpand(edist,ringNumber);
}
//____________________________________________________________________
void AliFMDAnaCalibEnergyDistribution::Browse(TBrowser* b)
{
  //Browse object
  for(Int_t i = 0; i<fNetaBins; i++) {
    TObjArray* etaArray = static_cast<TObjArray*>(fArray.At(i));
    etaArray->SetName(Form("etabin_%03d", i+1));
    for(Int_t det = 1; det<=3;det++) {
      TObjArray* detArray = static_cast<TObjArray*>(etaArray->At(det));
      detArray->SetName(Form("FMD%d", det));
    }
  }
  
  for(Int_t det = 1; det<=3;det++) {
    TObjArray* detArray = static_cast<TObjArray*>(fEmptyArray.At(det));
    detArray->SetName(Form("FMD%d", det));
  }
    
  for(Int_t det = 1; det<=3;det++) {
    TObjArray* detArray = static_cast<TObjArray*>(fRingArray.At(det));
    detArray->SetName(Form("FMD%d", det));
  }

  b->Add(&fArray,     "etabins");
  b->Add(&fEmptyArray,"empty");
  b->Add(&fRingArray, "rings");
}

//____________________________________________________________________
//
// EOF
//
