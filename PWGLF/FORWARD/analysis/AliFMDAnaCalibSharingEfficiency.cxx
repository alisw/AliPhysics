
#include "AliFMDAnaCalibSharingEfficiency.h"
#include <TH2F.h>
#include <TH1F.h>
#include <TBrowser.h>

ClassImp(AliFMDAnaCalibSharingEfficiency)
#if 0
; // For Emacs
#endif 

//____________________________________________________________________
AliFMDAnaCalibSharingEfficiency::AliFMDAnaCalibSharingEfficiency() : TObject(),
								     fArray(),
								     fArrayTrVtx(),
								     fIsInit(kFALSE)
{
  
  
  
}


//____________________________________________________________________
AliFMDAnaCalibSharingEfficiency::AliFMDAnaCalibSharingEfficiency(const AliFMDAnaCalibSharingEfficiency& o)
  : TObject(o), fArray(o.fArray), fArrayTrVtx(o.fArrayTrVtx), fIsInit(o.fIsInit)
{
  // Copy ctor 
}
//____________________________________________________________________
AliFMDAnaCalibSharingEfficiency&
AliFMDAnaCalibSharingEfficiency::operator=(const AliFMDAnaCalibSharingEfficiency& o) 
{
  // Assignment operator 
  
  fArray       = o.fArray;
  fArrayTrVtx  = o.fArrayTrVtx;
  
  return (*this);
}
//____________________________________________________________________
void AliFMDAnaCalibSharingEfficiency::Init() {
  
  fArray.SetOwner();
  for(Int_t det = 1; det<=3;det++) {
    TObjArray* detArray = new TObjArray();
    detArray->SetOwner();
    fArray.AddAtAndExpand(detArray,det);
    Int_t nRings = (det == 1 ? 1 : 2);
    for(Int_t ring = 0;ring<nRings;ring++) {
      TObjArray* ringArray = new TObjArray();
      ringArray->SetOwner();
      detArray->AddAtAndExpand(ringArray,ring);
      
    }
  }
  fArrayTrVtx.SetOwner();
  for(Int_t det = 1; det<=3;det++) {
    TObjArray* detArrayTrVtx = new TObjArray();
    detArrayTrVtx->SetOwner();
    fArrayTrVtx.AddAtAndExpand(detArrayTrVtx,det);
    Int_t nRings = (det == 1 ? 1 : 2);
    for(Int_t ring = 0;ring<nRings;ring++) {
      TObjArray* ringArrayTrVtx = new TObjArray();
      ringArrayTrVtx->SetOwner();
      detArrayTrVtx->AddAtAndExpand(ringArrayTrVtx,ring);
      
    }
  }
  fIsInit = kTRUE;
 
}
//____________________________________________________________________
TObjArray* AliFMDAnaCalibSharingEfficiency::GetRingArrayTrVtx(Int_t det, 
							      Char_t ring) {
  
  Int_t ringNumber      = (ring == 'I' ? 0 : 1);
  TObjArray* detArray  = (TObjArray*)fArrayTrVtx.At(det);
  TObjArray* ringArray = (TObjArray*)detArray->At(ringNumber);
  
  return ringArray;
}
//____________________________________________________________________
void AliFMDAnaCalibSharingEfficiency::SetSharingEffTrVtx(Int_t det, 
							 Char_t ring, 
							 Int_t vtxbin, 
							 TH1F* hCorrection) {
  if(!fIsInit)
    Init();
  
  TObjArray* ringArray = GetRingArrayTrVtx(det,ring);
  ringArray->AddAtAndExpand(hCorrection,vtxbin);
  
}
//____________________________________________________________________
TH1F* AliFMDAnaCalibSharingEfficiency::GetSharingEffTrVtx(Int_t det, 
							  Char_t ring, 
							  Int_t vtxbin) {
  TObjArray* ringArray = GetRingArrayTrVtx(det,ring);
  TH1F* hCorrection    = (TH1F*)ringArray->At(vtxbin);
  return hCorrection;
}

//____________________________________________________________________
TObjArray* AliFMDAnaCalibSharingEfficiency::GetRingArray(Int_t det, 
							 Char_t ring) {
  
  Int_t ringNumber      = (ring == 'I' ? 0 : 1);
  TObjArray* detArray  = (TObjArray*)fArray.At(det);
  TObjArray* ringArray = (TObjArray*)detArray->At(ringNumber);
  
  return ringArray;
}
//____________________________________________________________________
void AliFMDAnaCalibSharingEfficiency::SetSharingEff(Int_t det, 
						    Char_t ring, 
						    Int_t vtxbin, 
						    TH1F* hCorrection) {
  if(!fIsInit)
    Init();
  
  TObjArray* ringArray = GetRingArray(det,ring);
  ringArray->AddAtAndExpand(hCorrection,vtxbin);
  
}
//____________________________________________________________________
TH1F* AliFMDAnaCalibSharingEfficiency::GetSharingEff(Int_t det, 
						     Char_t ring, 
						     Int_t vtxbin) {
  TObjArray* ringArray = GetRingArray(det,ring);
  TH1F* hCorrection    = (TH1F*)ringArray->At(vtxbin);
  return hCorrection;
}

//____________________________________________________________________
void
AliFMDAnaCalibSharingEfficiency::Browse(TBrowser* b)
{
  b->Add(&fArray, "Array of histograms w/sharing eff corrections");
}

//____________________________________________________________________
//
// EOF
//
