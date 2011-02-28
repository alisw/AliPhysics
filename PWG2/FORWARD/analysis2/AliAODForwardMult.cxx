//
// Class that contains the forward multiplicity data per event 
//
// This class contains a histogram of 
// @f[
//   \frac{d^2N_{ch}}{d\eta d\phi}\quad,
// @f]
// as well as a trigger mask for each analysed event.  
// 
// The eta acceptance of the event is stored in the underflow bins of
// the histogram.  So to build the final histogram, one needs to
// correct for this acceptance (properly weighted by the events), and
// the vertex efficiency.  This simply boils down to defining a 2D
// histogram and summing the event histograms in that histogram.  One
// should of course also do proper book-keeping of the accepted event.
//
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
  : fIsMC(false),
    fHist(),
    fTriggers(0),
    fIpZ(fgkInvalidIpZ), 
    fCentrality(-1)
{
  // 
  // Constructor 
  // 
}

//____________________________________________________________________
AliAODForwardMult::AliAODForwardMult(Bool_t isMC) 
  : fIsMC(isMC),
    fHist("forwardMult", "d^{2}N_{ch}/d#etad#varphi in the forward regions", 
	  200, -4, 6, 20, 0, 2*TMath::Pi()),
    fTriggers(0),
    fIpZ(fgkInvalidIpZ), 
    fCentrality(-1)
{
  // 
  // Constructor 
  // 
  // Parameters: 
  //  isMC   If set to true this is for MC data (effects branch name)
  // 
  fHist.SetXTitle("#eta");
  fHist.SetYTitle("#varphi [radians]");
  fHist.SetZTitle("#frac{d^{2}N_{ch}}{d#etad#varphi}");
  fHist.SetDirectory(0);
  fHist.Sumw2();
}

//____________________________________________________________________
void
AliAODForwardMult::Init(const TAxis& etaAxis)
{
  // Initialize the histogram with an eta axis 
  // 
  // Parameters: 
  //   etaAxis       Eta axis to use 
  // 
  fHist.SetBins(etaAxis.GetNbins(), etaAxis.GetXmin(), etaAxis.GetXmax(), 
		20, 0, 2*TMath::Pi());
}

//____________________________________________________________________
void
AliAODForwardMult::Clear(Option_t* option)
{
  // Clear (or reset) internal values 
  // 
  // Parameters: 
  //  option   Passed to TH1::Reset 
  // 
  fHist.Reset(option);
  fTriggers = 0;
  fIpZ      = fgkInvalidIpZ;
}
//____________________________________________________________________
void
AliAODForwardMult::SetSNN(UShort_t snn)
{
  // set the center of mass energy per nucleon pair (GeV). 
  // This is stored in bin (0,0) of the histogram 
  // 
  // Parameters: 
  //   sNN   Center of mass energy per nuclean 
  fHist.SetBinContent(0,0,snn);
}
//____________________________________________________________________
void
AliAODForwardMult::SetSystem(UShort_t sys)
{
  // set the center of mass energy per nucleon pair (GeV). 
  // This is stored in bin (N+1,0) of the histogram 
  // 
  // Parameters: 
  //   sys   Collision system number 
  fHist.SetBinContent(fHist.GetNbinsX()+1,0,sys);
}

//____________________________________________________________________
Bool_t
AliAODForwardMult::HasIpZ() const
{
  // Check if we have valid z coordinate of the interaction point 
  // 
  // Return:
  //   true if the z coordinate of the interaction point is valid 
  // 
  return TMath::Abs(fIpZ - fgkInvalidIpZ) > 1;
}

//____________________________________________________________________
UShort_t
AliAODForwardMult::GetSNN() const
{
  // set the center of mass energy per nucleon pair (GeV). 
  // This is stored in bin (0,0) of the histogram 
  // 
  // Parameters: 
  //   sNN   Center of mass energy per nuclean 
  return UShort_t(fHist.GetBinContent(0,0));
}

//____________________________________________________________________
UShort_t
AliAODForwardMult::GetSystem() const
{
  // set the center of mass energy per nucleon pair (GeV). 
  // This is stored in bin (N+1,0) of the histogram 
  // 
  // Parameters: 
  //   sNN   Center of mass energy per nuclean 
  return UShort_t(fHist.GetBinContent(fHist.GetNbinsX()+1,0));
}

//____________________________________________________________________
void
AliAODForwardMult::Browse(TBrowser* b)
{
  // Browse this object 
  // 
  // Parameters: 
  //   b   Browser to use 
  static TObjString ipz;
  static TObjString trg;
  static TObjString cnt;
  ipz = Form("ip_z=%fcm", fIpZ);
  trg = GetTriggerString(fTriggers);
  cnt = Form("%+6.1f%%", fCentrality);
  b->Add(&fHist);
  b->Add(&ipz);
  b->Add(&trg);
  b->Add(&cnt);
}

//____________________________________________________________________
const Char_t*
AliAODForwardMult::GetTriggerString(UInt_t mask)
{
  // Get a string that describes the triggers 
  // 
  // Parameters: 
  //   mask  Bit pattern of triggers 
  // Return: 
  //   Character string representation of mask 
  static TString trg;
  trg = "";
  if ((mask & kInel)        != 0x0) trg.Append("INEL ");
  if ((mask & kInelGt0)     != 0x0) trg.Append("INEL>0 ");
  if ((mask & kNSD)         != 0x0) trg.Append("NSD ");
  if ((mask & kA)           != 0x0) trg.Append("A ");
  if ((mask & kB)           != 0x0) trg.Append("B ");
  if ((mask & kC)           != 0x0) trg.Append("C ");
  if ((mask & kE)           != 0x0) trg.Append("E ");
  if ((mask & kMCNSD)       != 0x0) trg.Append("MCNSD ");
  return trg.Data();
}
  
//____________________________________________________________________
void
AliAODForwardMult::Print(Option_t* option) const
{
  // Print this object 
  // 
  // Parameters: 
  //  option   Passed to TH1::Print 
  fHist.Print(option);
  UShort_t sys = GetSystem();
  TString  str = "unknown";
  switch (sys) { 
  case 1:  str = "pp"; break;
  case 2:  str = "PbPb"; break;
  }
  std::cout << "Ipz:         " << fIpZ << "cm " << (HasIpZ() ? "" : "in") 
	    << "valid\n"
	    << "Triggers:    " << GetTriggerString(fTriggers) 
	    << "sNN:         " << GetSNN() << "GeV\n" 
	    << "System:      " << str 
	    << "Centrality:  " << fCentrality << "%" 
	    << std::endl;
}

//____________________________________________________________________
//
// EOF
//
