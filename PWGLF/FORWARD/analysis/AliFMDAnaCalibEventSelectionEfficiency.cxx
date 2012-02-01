
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
										   fCorrectionList(),
										   fIsInit(kFALSE),
										   fVtxEff(1)
{
  
  
  
}


//____________________________________________________________________
AliFMDAnaCalibEventSelectionEfficiency::
AliFMDAnaCalibEventSelectionEfficiency(const AliFMDAnaCalibEventSelectionEfficiency& o) : TObject(o),			     
											  fCorrection(o.fCorrection),
											  fCorrectionList(),
											  fIsInit(o.fIsInit),
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
void AliFMDAnaCalibEventSelectionEfficiency::SetCorrection(TString trig,
							   Int_t vtxbin, 
							   Char_t ring,
							   TH2F* hCorrection) {
  //TString test = trig;
  if(!trig.Contains("INEL") && !trig.Contains("NSD"))
    AliWarning("Please choose NSD or INEL!");
  
  if(trig.Contains("INEL"))
    hCorrection->SetName(Form("correction_%c_%d",ring,vtxbin));
  if(trig.Contains("NSD")) 
    hCorrection->SetName(Form("correction%s_%c_%d","NSD",ring,vtxbin));
  
  fCorrectionList.Add(hCorrection);
  
}
//____________________________________________________________________
TH2F* AliFMDAnaCalibEventSelectionEfficiency::GetCorrection(TString name,
							    Int_t vtxbin, 
							    Char_t ring) {
  
  //TString name = trig;
  
  if(name.Contains("INEL"))
    name.Form("correction_%c_%d",ring,vtxbin);
  if(name.Contains("NSD")) 
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
