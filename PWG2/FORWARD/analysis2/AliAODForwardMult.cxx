#include "AliAODForwardMult.h"
#include <TBrowser.h>
#include <iostream>
#include <TMath.h>
#include <TObjString.h>

ClassImp(AliAODForwardMult)
#if 0 
; // For Emacs 
#endif

//____________________________________________________________________
const Float_t AliAODForwardMult::fgkInvalidIpZ = 1e6;

//____________________________________________________________________
AliAODForwardMult::AliAODForwardMult()
  : fHist(),
    fTriggers(0),
    fIpZ(fgkInvalidIpZ)
{}

//____________________________________________________________________
AliAODForwardMult::AliAODForwardMult(Bool_t) 
  : fHist("forwardMult", "d^{2}N_{ch}/d#etad#varphi in the forward regions", 
	  200, -4, 6, 20, 0, 2*TMath::Pi()),
    fTriggers(0),
    fIpZ(fgkInvalidIpZ)
{
  fHist.SetXTitle("#eta");
  fHist.SetYTitle("#varphi [radians]");
  fHist.SetZTitle("#frac{d^{2}N_{ch}}{d#etad#varphi}");
  fHist.Sumw2();
}

//____________________________________________________________________
void
AliAODForwardMult::Init(const TAxis& etaAxis)
{
  fHist.SetBins(etaAxis.GetNbins(), etaAxis.GetXmin(), etaAxis.GetXmax(), 
		20, 0, 2*TMath::Pi());
}

//____________________________________________________________________
void
AliAODForwardMult::Clear(Option_t* option)
{
  fHist.Reset(option);
  fTriggers = 0;
  fIpZ      = fgkInvalidIpZ;
}
//____________________________________________________________________
Bool_t
AliAODForwardMult::HasIpZ() const
{
  return TMath::Abs(fIpZ - fgkInvalidIpZ) > 1;
}

//____________________________________________________________________
void
AliAODForwardMult::Browse(TBrowser* b)
{
  static TObjString ipz;
  static TObjString trg;
  ipz = Form("ip_z=%fcm", fIpZ);
  trg = GetTriggerString(fTriggers);
  b->Add(&fHist);
  b->Add(&ipz);
  b->Add(&trg);
}

//____________________________________________________________________
const Char_t*
AliAODForwardMult::GetTriggerString(UInt_t mask)
{
  static TString trg;
  trg = "";
  if ((mask & kInel)    != 0x0) trg.Append("INEL ");
  if ((mask & kInelGt0) != 0x0) trg.Append("INEL>0 ");
  if ((mask & kNSD)     != 0x0) trg.Append("NSD ");
  if ((mask & kA)       != 0x0) trg.Append("A ");
  if ((mask & kB)       != 0x0) trg.Append("B ");
  if ((mask & kC)       != 0x0) trg.Append("C ");
  if ((mask & kE)       != 0x0) trg.Append("E ");
  return trg.Data();
}
  
//____________________________________________________________________
void
AliAODForwardMult::Print(Option_t* option) const
{
  fHist.Print(option);
  std::cout << "Ipz:      " << fIpZ << "cm " << (HasIpZ() ? "" : "in") 
	    << "valid\n"
	    << "Triggers: " << GetTriggerString(fTriggers) << std::endl;
}

//____________________________________________________________________
//
// EOF
//
