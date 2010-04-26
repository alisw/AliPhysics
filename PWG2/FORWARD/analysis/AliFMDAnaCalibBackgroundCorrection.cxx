
#include "AliFMDAnaCalibBackgroundCorrection.h"
#include <TH2F.h>
#include <TH1F.h>
#include <TBrowser.h>

ClassImp(AliFMDAnaCalibBackgroundCorrection)
#if 0
; // For Emacs
#endif 

//____________________________________________________________________
AliFMDAnaCalibBackgroundCorrection::AliFMDAnaCalibBackgroundCorrection() : TObject(),
									   fArray(),
									   fAxis(),
									   fIsInit(kFALSE),
									   fListOfDoubleHitCorrection()
{
  
  
  
}


//____________________________________________________________________
AliFMDAnaCalibBackgroundCorrection::AliFMDAnaCalibBackgroundCorrection(const AliFMDAnaCalibBackgroundCorrection& o)
  : TObject(o), fArray(o.fArray), fAxis(o.fAxis), fIsInit(o.fIsInit), fListOfDoubleHitCorrection()
{
  // Copy ctor 
}
//____________________________________________________________________
AliFMDAnaCalibBackgroundCorrection&
AliFMDAnaCalibBackgroundCorrection::operator=(const AliFMDAnaCalibBackgroundCorrection& o) 
{
  // Assignment operator 
  
  fArray     = o.fArray;
  
  return (*this);
}

//____________________________________________________________________
TH2F* AliFMDAnaCalibBackgroundCorrection::GetBgCorrection(Int_t det, 
							  Char_t ring, 
							  Int_t vtxbin) {
  TObjArray* ringArray = GetRingArray(det,ring);
  TH2F* hCorrection    = (TH2F*)ringArray->At(vtxbin);
  return hCorrection;
}

//____________________________________________________________________
void AliFMDAnaCalibBackgroundCorrection::SetBgCorrection(Int_t det, 
							 Char_t ring, 
							 Int_t vtxbin, 
							 TH2F* hCorrection) {
  if(!fIsInit)
    Init();
  
  TObjArray* ringArray = GetRingArray(det,ring);
  ringArray->AddAtAndExpand(hCorrection,vtxbin);
  
}
//____________________________________________________________________
TH1F* AliFMDAnaCalibBackgroundCorrection::GetDoubleHitCorrection(Int_t  det, 
								 Char_t ring) {
  
  
  TH1F* hCorrection    = (TH1F*)fListOfDoubleHitCorrection.FindObject(Form("hDoubleHitCorrection_FMD%d%c",det,ring));
  return hCorrection;
}

//____________________________________________________________________
void AliFMDAnaCalibBackgroundCorrection::SetDoubleHitCorrection(Int_t det, 
								Char_t ring, 
								TH1F* hCorrection) {
  hCorrection->SetName(Form("hDoubleHitCorrection_FMD%d%c",det,ring));
  fListOfDoubleHitCorrection.Add(hCorrection);    
}
//____________________________________________________________________
TH1F* AliFMDAnaCalibBackgroundCorrection::GetSPDDeadCorrection(Int_t  vtxbin) {
  
  TH1F* hCorrection    = (TH1F*)fListOfDoubleHitCorrection.FindObject(Form("hSPDDeadCorrection_vtx%d",vtxbin));
  return hCorrection;
}

//____________________________________________________________________
void AliFMDAnaCalibBackgroundCorrection::SetSPDDeadCorrection(Int_t vtxbin, 
							      TH1F* hCorrection) {
  hCorrection->SetName(Form("hSPDDeadCorrection_vtx%d",vtxbin));
  fListOfDoubleHitCorrection.Add(hCorrection);    
}
//____________________________________________________________________
void AliFMDAnaCalibBackgroundCorrection::SetRefAxis(TAxis* axis) {
  
 
  fAxis.Set(axis->GetNbins(),axis->GetXmin(),axis->GetXmax());
    
}
//____________________________________________________________________
void AliFMDAnaCalibBackgroundCorrection::Init() {
  
  fArray.SetOwner();
  
  TObjArray* spdArray = new TObjArray();
  spdArray->SetOwner();
  fArray.AddAtAndExpand(spdArray,0);
  
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
  fIsInit = kTRUE;
  
}
//____________________________________________________________________
TObjArray* AliFMDAnaCalibBackgroundCorrection::GetRingArray(Int_t det, 
							    Char_t ring) {
  
  if(det==0 || det == 4) {
    TObjArray* spdArray  = (TObjArray*)fArray.At(det);
    return spdArray;
  }
  Int_t ringNumber      = (ring == 'I' ? 0 : 1);
  TObjArray* detArray  = (TObjArray*)fArray.At(det);
  TObjArray* ringArray = (TObjArray*)detArray->At(ringNumber);
  
  return ringArray;
}
//____________________________________________________________________
Int_t AliFMDAnaCalibBackgroundCorrection::GetNvtxBins() {
  
  return fAxis.GetNbins();
  
}
//____________________________________________________________________
Float_t AliFMDAnaCalibBackgroundCorrection::GetVtxCutZ() {
  
  return fAxis.GetXmax();
  
}

//____________________________________________________________________
void
AliFMDAnaCalibBackgroundCorrection::Browse(TBrowser* b)
{
  b->Add(&fAxis, "Vertex bins");
  b->Add(&fArray, "Array of histograms w/BG corrections");
}

//____________________________________________________________________
//
// EOF
//
