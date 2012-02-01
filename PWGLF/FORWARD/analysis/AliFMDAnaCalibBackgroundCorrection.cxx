//
// Calibration/Correction object that stores the secondary to primary
// correction used in the FMD analysis.
// 
#include "AliFMDAnaCalibBackgroundCorrection.h"
#include <TH2F.h>
#include <TH1F.h>
#include <TBrowser.h>

ClassImp(AliFMDAnaCalibBackgroundCorrection)
#if 0
; // For Emacs
#endif 

//____________________________________________________________________
AliFMDAnaCalibBackgroundCorrection::AliFMDAnaCalibBackgroundCorrection() 
  : TObject(),
    fArray(),
    fAxis(),
    fIsInit(kFALSE),
    fListOfDoubleHitCorrection(),
    fListOfNSDBgMaps()
{
}


//____________________________________________________________________
AliFMDAnaCalibBackgroundCorrection::AliFMDAnaCalibBackgroundCorrection(const AliFMDAnaCalibBackgroundCorrection& o)
  : TObject(o), 
    fArray(o.fArray), 
    fAxis(o.fAxis), 
    fIsInit(o.fIsInit), 
    fListOfDoubleHitCorrection(), 
    fListOfNSDBgMaps()
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
void 
AliFMDAnaCalibBackgroundCorrection::SetBgCorrection(Int_t det, 
						    Char_t ring, 
						    Int_t vtxbin, 
						    TH2F* hCorrection) 
{
  if(!fIsInit) Init();
  
  TObjArray* ringArray = GetRingArray(det,ring);
  if (!ringArray) return;

  ringArray->AddAtAndExpand(hCorrection,vtxbin);
  
}
//____________________________________________________________________
TH2F* 
AliFMDAnaCalibBackgroundCorrection::GetBgCorrection(Int_t det, 
						    Char_t ring, 
						    Int_t vtxbin) const
{
  TObjArray* ringArray   = GetRingArray(det,ring);
  if (!ringArray) return 0;

  TH2F*      hCorrection = static_cast<TH2F*>(ringArray->At(vtxbin));
  return hCorrection;
}
//____________________________________________________________________
void 
AliFMDAnaCalibBackgroundCorrection::SetNSDBgCorrection(Int_t det, 
						       Char_t ring, 
						       Int_t vtxbin, 
						       TH2F* hCorrection) 
{
  if(!fIsInit)  Init();
  hCorrection->SetName(Form("FMDNSD%d%c_vtxbin_%d_correction",det,ring,vtxbin));
  fListOfNSDBgMaps.Add(hCorrection);
    
}
//____________________________________________________________________
TH2F* 
AliFMDAnaCalibBackgroundCorrection::GetNSDBgCorrection(Int_t det, 
						       Char_t ring, 
						       Int_t vtxbin) const
{
  TH2F* hCorrection    = 
    static_cast<TH2F*>(fListOfNSDBgMaps
		       .FindObject(Form("FMDNSD%d%c_vtxbin_%d_correction",
					det,ring,vtxbin)));
  return hCorrection;
}
//____________________________________________________________________
TH1F* 
AliFMDAnaCalibBackgroundCorrection::GetDoubleHitCorrection(Int_t  det, 
							   Char_t ring) const
{  
  TH1F* hCorrection    = 
    static_cast<TH1F*>(fListOfDoubleHitCorrection
		       .FindObject(Form("hDoubleHitCorrection_FMD%d%c",
					det,ring)));
  return hCorrection;
}

//____________________________________________________________________
void 
AliFMDAnaCalibBackgroundCorrection::SetDoubleHitCorrection(Int_t det, 
							   Char_t ring, 
							   TH1F* hCorrection) 
{
  hCorrection->SetName(Form("hDoubleHitCorrection_FMD%d%c",det,ring));
  fListOfDoubleHitCorrection.Add(hCorrection);    
}
//____________________________________________________________________
TH1F* 
AliFMDAnaCalibBackgroundCorrection::GetSPDDeadCorrection(Int_t  vtxbin) const
{
  TH1F* hCorrection    = 
    static_cast<TH1F*>(fListOfDoubleHitCorrection
		       .FindObject(Form("hSPDDeadCorrection_vtx%d",vtxbin)));
  return hCorrection;
}

//____________________________________________________________________
void 
AliFMDAnaCalibBackgroundCorrection::SetSPDDeadCorrection(Int_t vtxbin, 
							 TH1F* hCorrection) 
{
  hCorrection->SetName(Form("hSPDDeadCorrection_vtx%d",vtxbin));
  fListOfDoubleHitCorrection.Add(hCorrection);    
}

//____________________________________________________________________
TH1F* 
AliFMDAnaCalibBackgroundCorrection::GetFMDDeadCorrection(Int_t  vtxbin) 
{
  TH1F* hCorrection    = 
    static_cast<TH1F*>(fListOfDoubleHitCorrection
		       .FindObject(Form("hFMDDeadCorrection_vtx%d",vtxbin)));
  return hCorrection;
}

//____________________________________________________________________
void 
AliFMDAnaCalibBackgroundCorrection::SetFMDDeadCorrection(Int_t vtxbin, 
							 TH1F* hCorrection) 
{
  hCorrection->SetName(Form("hFMDDeadCorrection_vtx%d",vtxbin));
  fListOfDoubleHitCorrection.Add(hCorrection);    
}

//____________________________________________________________________
void 
AliFMDAnaCalibBackgroundCorrection::SetRefAxis(TAxis* axis) 
{
  // Set the reference axis 
  fAxis.Set(axis->GetNbins(),axis->GetXmin(),axis->GetXmax());
}
//____________________________________________________________________
void AliFMDAnaCalibBackgroundCorrection::Init() 
{
  // Initialize 
  fArray.SetOwner();
  
  TObjArray* spdArray = new TObjArray();
  spdArray->SetOwner();
  fArray.AddAtAndExpand(spdArray,0);
  
  for(Int_t det = 1; det<=3;det++) {
    TObjArray* detArray = new TObjArray();
    detArray->SetOwner();
    detArray->SetName(Form("FMD%d", det));
    // detArray->SetTitle(Form("Array of FMD%d corrections", det));
    fArray.AddAtAndExpand(detArray,det);
    Int_t nRings = (det == 1 ? 1 : 2);
    for(Int_t ring = 0;ring<nRings;ring++) {
      TObjArray* ringArray = new TObjArray();
      Char_t r = (ring == 0 ? 'I' : 'O');
      ringArray->SetOwner();
      ringArray->SetName(Form("FMD%d%c", det, r));
      // ringArray->SetTitle(Form("Array of FMD%d%c corrections", det, r));
      detArray->AddAtAndExpand(ringArray,ring);
      
    }
  }
  fIsInit = kTRUE;
  
}
//____________________________________________________________________
TObjArray* AliFMDAnaCalibBackgroundCorrection::GetRingArray(Int_t det, 
							    Char_t ring) const
{
  // Find array corresponding to det, ring.  
  // Note, 0 and 4 refers to the SPD arrays 
  if(det==0 || det == 4) {
    TObjArray* spdArray  = (TObjArray*)fArray.At(det);
    return spdArray;
  }

  if (det < 0 || det >= fArray.GetEntriesFast()) return 0; 

  Int_t      ringNumber = (ring == 'I' ? 0 : 1);
  TObjArray* detArray   = static_cast<TObjArray*>(fArray.At(det));
  TObjArray* ringArray  = static_cast<TObjArray*>(detArray->At(ringNumber));
  
  return ringArray;
}

//____________________________________________________________________
void
AliFMDAnaCalibBackgroundCorrection::Browse(TBrowser* b)
{
  // Browse 
  b->Add(&fAxis, "Vertex bins");
  b->Add(&fArray, "Array of histograms w/BG corrections");
}

//____________________________________________________________________
//
// EOF
//
