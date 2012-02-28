// 
// Class that contains results from FMD eventplane calculations
//
#include "AliAODForwardEP.h"
#include <TBrowser.h>
#include <iostream>
#include <TMath.h>
#include <TObjString.h>
#include <TObjArray.h>
#include "AliLog.h"
ClassImp(AliAODForwardEP)
#ifdef DOXY_INPUT
; // For Emacs 
#endif

//____________________________________________________________________
AliAODForwardEP::AliAODForwardEP()
  : fIsMC(false),
    fEpT(-1),
    fEpA(-1),
    fEpC(-1),
    fEp1(-1),
    fEp2(-1),
    fHist()
{
  // 
  // Constructor 
  // 
}

//____________________________________________________________________
AliAODForwardEP::AliAODForwardEP(Bool_t isMC) 
  : fIsMC(isMC),
    fEpT(-1),
    fEpA(-1),
    fEpC(-1),
    fEp1(-1),
    fEp2(-1),
    fHist()
{
  // 
  // Constructor 
  // 
  // Parameters: 
  //  isMC   If set to true this is for MC data (effects branch name)
  // 
  fHist.SetXTitle("#eta");
  fHist.SetDirectory(0);
  fHist.Sumw2();
}

//____________________________________________________________________
void
AliAODForwardEP::Init(const TAxis& etaAxis)
{
  // Initialize the histogram with an eta axis 
  // 
  // Parameters: 
  //   etaAxis       Eta axis to use 
  // 
  fHist.SetBins(etaAxis.GetNbins()/10, etaAxis.GetXmin(), etaAxis.GetXmax()); 
}

//____________________________________________________________________
void
AliAODForwardEP::Clear(Option_t* option)
{
  // Clear (or reset) internal values 
  // 
  // Parameters: 
  //  option   Passed to TH1::Reset 
  // 
  fHist.Reset(option);
  fEpT = -1;
  fEpA = -1;
  fEpC = -1;
  fEp1 = -1;
  fEp2 = -1;
}
//____________________________________________________________________
void
AliAODForwardEP::Browse(TBrowser* /*b*/)
{
  // Browse this object 
  // 
  // Parameters: 
  //   b   Browser to use 
  // TODO: Make nice
/*  static TObjString ipz;
  static TObjString trg;
  static TObjString cnt;
  static TObjString ncl;
  ipz = Form("ip_z=%fcm", fIpZ);
  trg = GetTriggerString(fTriggers);
  cnt = Form("%+6.1f%%", fCentrality);
  ncl = Form("%d clusters", fNClusters);
  b->Add(&fHist);
  b->Add(&ipz);
  b->Add(&trg);
  b->Add(&cnt);
  b->Add(&ncl);
*/
}

//____________________________________________________________________
void
AliAODForwardEP::Print(Option_t* option) const
{
  // Print this object 
  // 
  // Parameters: 
  //  option   Not used 
  fHist.Print(option);
  std::cout << "Total FMD EP: \t" << fEpT <<std::endl; 
  std::cout << "FMD1+2 EP: \t" << fEpA <<std::endl; 
  std::cout << "FMD3   EP: \t" << fEpC <<std::endl; 
  std::cout << "Random Subep 1: \t" << fEp1 <<std::endl; 
  std::cout << "Random Subep 2: \t" << fEp2 <<std::endl; 
}

//____________________________________________________________________
//
// EOF
//
