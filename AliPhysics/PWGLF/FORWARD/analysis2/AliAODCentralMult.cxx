//
// Class that contains the central multiplicity data per event 
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
#include "AliAODCentralMult.h"
#include <TBrowser.h>
#include <iostream>
#include <TMath.h>
#include <TObjString.h>

ClassImp(AliAODCentralMult)
#if 0 
; // For Emacs 
#endif

//____________________________________________________________________
AliAODCentralMult::AliAODCentralMult()
  : fIsMC(false),
    fHist()
{
  // 
  // Constructor 
  // 
}

//____________________________________________________________________
AliAODCentralMult::AliAODCentralMult(Bool_t isMC) 
  : fIsMC(isMC),
    fHist("centralMult", "d^{2}N_{ch}/d#etad#varphi in the central regions", 
	  200, -4, 6, 20, 0, 2*TMath::Pi())
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
AliAODCentralMult::Clear(Option_t*) {
  
  fHist.Reset();
  
}
//____________________________________________________________________
void
AliAODCentralMult::Init(const TAxis& etaAxis)
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
AliAODCentralMult::Browse(TBrowser* b)
{
  // Browse this object 
  // 
  // Parameters: 
  //   b   Browser to use 

  b->Add(&fHist);

}
//____________________________________________________________________
void
AliAODCentralMult::Print(Option_t* option) const
{
  // Print this object 
  // 
  // Parameters: 
  //  option   Passed to TH1::Print 
  fHist.Print(option);
}

//____________________________________________________________________
//
// EOF
//
