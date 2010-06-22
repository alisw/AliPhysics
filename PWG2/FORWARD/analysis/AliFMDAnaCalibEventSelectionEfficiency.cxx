
#include "AliFMDAnaCalibEventSelectionEfficiency.h"

#include <TH1F.h>
#include <TBrowser.h>

ClassImp(AliFMDAnaCalibEventSelectionEfficiency)
#if 0
; // For Emacs
#endif 

//____________________________________________________________________
AliFMDAnaCalibEventSelectionEfficiency::AliFMDAnaCalibEventSelectionEfficiency() : TObject(),
								     fCorrection(),
								     fIsInit(kFALSE)
{
  
  
  
}


//____________________________________________________________________
AliFMDAnaCalibEventSelectionEfficiency::AliFMDAnaCalibEventSelectionEfficiency(const AliFMDAnaCalibEventSelectionEfficiency& o) : TObject(o), 								     fCorrection(o.fCorrection),						     fIsInit(o.fIsInit)
{
  // Copy ctor 
}
//____________________________________________________________________
AliFMDAnaCalibEventSelectionEfficiency&
AliFMDAnaCalibEventSelectionEfficiency::operator=(const AliFMDAnaCalibEventSelectionEfficiency& /*o*/) 
{
  // Assignment operator 
  
  return (*this);
}
//____________________________________________________________________
void AliFMDAnaCalibEventSelectionEfficiency::Init() {

  fCorrection.SetName("EventSelectionEffCorrection");
  
  fIsInit = kTRUE;

}
//____________________________________________________________________
void AliFMDAnaCalibEventSelectionEfficiency::SetCorrection(TH1F* hCorrection) {
  
  fCorrection.SetBins(hCorrection->GetNbinsX(),
		      hCorrection->GetXaxis()->GetXmin(),
		      hCorrection->GetXaxis()->GetXmax());
  for(Int_t i=1; i<=hCorrection->GetNbinsX(); i++) {
    fCorrection.SetBinContent(i,hCorrection->GetBinContent(i));
    fCorrection.SetBinError(i,hCorrection->GetBinError(i));
  }
  

}
//____________________________________________________________________
void AliFMDAnaCalibEventSelectionEfficiency::SetCorrection(Int_t vtxbin, 
							   Char_t ring,
							   TH2F* hCorrection) {
  hCorrection->SetName(Form("correction_%c_%d",ring,vtxbin));
  fCorrectionList.Add(hCorrection);
  
}
//____________________________________________________________________
TH2F* AliFMDAnaCalibEventSelectionEfficiency::GetCorrection(Int_t vtxbin, 
							    Char_t ring) {
  
  TH2F* hCorrection = (TH2F*)fCorrectionList.FindObject(Form("correction_%c_%d",ring,vtxbin));
  return hCorrection;

}

//____________________________________________________________________
Float_t AliFMDAnaCalibEventSelectionEfficiency::GetCorrection(Int_t vtxbin) {

  if( (vtxbin-1) > fCorrection.GetNbinsX() || vtxbin < 0)
    return 0;
  
  Float_t correction = fCorrection.GetBinContent(vtxbin+1);
  
  return correction;

}

//____________________________________________________________________
void AliFMDAnaCalibEventSelectionEfficiency::Browse(TBrowser* /*b*/)
{
  
}

//____________________________________________________________________
//
// EOF
//
