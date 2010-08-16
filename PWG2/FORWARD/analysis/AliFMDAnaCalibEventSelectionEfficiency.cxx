
#include "AliFMDAnaCalibEventSelectionEfficiency.h"

#include <TH1F.h>
#include <TBrowser.h>
#include "AliLog.h"
#include "iostream"

ClassImp(AliFMDAnaCalibEventSelectionEfficiency)
#if 0
; // For Emacs
#endif 

//____________________________________________________________________
AliFMDAnaCalibEventSelectionEfficiency::AliFMDAnaCalibEventSelectionEfficiency() : TObject(),
										   fCorrection(),
										   fIsInit(kFALSE),
										   fCorrectionList(),
										   fVtxEff(1)
{
  
  
  
}


//____________________________________________________________________
AliFMDAnaCalibEventSelectionEfficiency::
AliFMDAnaCalibEventSelectionEfficiency(const AliFMDAnaCalibEventSelectionEfficiency& o) : TObject(o),			     
											  fCorrection(o.fCorrection),
											  fIsInit(o.fIsInit),
											  fCorrectionList(),
											  fVtxEff(o.fVtxEff)
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
void AliFMDAnaCalibEventSelectionEfficiency::SetCorrection(Char_t* trig,
							   Int_t vtxbin, 
							   Char_t ring,
							   TH2F* hCorrection) {
  if(trig != "INEL" && trig != "NSD")
    AliWarning("Please choose NSD or INEL!");
  
  if(trig == "INEL")
    hCorrection->SetName(Form("correction_%c_%d",ring,vtxbin));
  if(trig == "NSD") 
    hCorrection->SetName(Form("correction%s_%c_%d","NSD",ring,vtxbin));
  
  fCorrectionList.Add(hCorrection);
  
}
//____________________________________________________________________
TH2F* AliFMDAnaCalibEventSelectionEfficiency::GetCorrection(Char_t* trig,
							    Int_t vtxbin, 
							    Char_t ring) {
  
  TString name;
  
  if(trig == "INEL")
    name.Form("correction_%c_%d",ring,vtxbin);
  if(trig == "NSD") 
    name.Form("correction%s_%c_%d","NSD",ring,vtxbin);
  
  TH2F* hCorrection = (TH2F*)fCorrectionList.FindObject(name);
   
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
