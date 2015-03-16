/**
 * @file   AliPtEtaPidWeights.cxx
 * @author Marek Chojnacki <mchojnac#cern.ch>
 * @date   Wed Feb  4 00:09:32 2015
 * 
 * @brief  Implementation of pt-eta weights
 * 
 * 
 */

#include "AliPtEtaPidWeights.h"
#include "AliForwardUtil.h"
#include <TROOT.h>
#include <TH2.h>
#include <iostream>

//____________________________________________________________________
AliPtEtaPidWeights::AliPtEtaPidWeights()
  : AliBaseMCWeights(),
    fPdgs(10),
    fWeights()
{
  fWeights.SetOwner();
}
//____________________________________________________________________
AliPtEtaPidWeights::AliPtEtaPidWeights(const AliPtEtaPidWeights& o)
  : AliBaseMCWeights(o),
    fPdgs(o.fPdgs),
    fWeights()
{
  fWeights.SetOwner();
  for (int i=0;i<o.fWeights.GetEntries();i++) {
    TObject* obj = o.fWeights.At(i);
    fWeights.Add(obj->Clone());
  }
}

//____________________________________________________________________
AliPtEtaPidWeights&
AliPtEtaPidWeights::operator=(const AliPtEtaPidWeights& o)
{
  if (&o == this) return *this; 

  fWeights.Clear();
  
  AliBaseMCWeights::operator=(o);
  fPdgs = o.fPdgs;
  for (int i = 0; i < o.fWeights.GetEntries(); i++) {
    TObject* obj = o.fWeights.At(i);		
    fWeights.AddAt(obj->Clone(), i);
  }

  return *this;
}
//____________________________________________________________________
void AliPtEtaPidWeights::Init(TList*)
{}

//____________________________________________________________________
Double_t AliPtEtaPidWeights::CalcWeight(Double_t eta,
					Double_t pt,
					Double_t phi,
					Int_t id,
					Double_t phiR,
					Double_t b) const 
{
  for (Int_t i = 0; i < fPdgs.GetSize(); i++)  {
    Int_t v = fPdgs.At(i);
    if (v == 0) 
      break;
    if (v == id) {
      TH2*  hist = static_cast<TH2*>(fWeights.At(i));
      Int_t binpt = hist->GetYaxis()->FindBin(pt);
      Int_t bineta = hist->GetXaxis()->FindBin(eta);
      return hist->GetBinContent(bineta,binpt);
    }
  }
  return 1.;
}
//____________________________________________________________________
void AliPtEtaPidWeights::AddPDGCode(Int_t  pdg, TH2* weight)
{

  Int_t n   = fPdgs.GetSize();
  Int_t pos = 0;
  // Find first empty slot 
  for (pos = 0; pos < n; pos++) if (fPdgs.At(pos) == 0) break;

  // If no empty slot is found, resize to twice the size 
  if (pos >= n) {
    fPdgs.Set(2*n);
  }

  // Now set the slot 
  fPdgs.AddAt(pdg, pos);
  fWeights.AddAt(weight, pos);
}
//____________________________________________________________________
#define PF(N,V,...)					\
  AliForwardUtil::PrintField(N,V, ## __VA_ARGS__)
#define PFB(N,FLAG)							\
  do {									\
    AliForwardUtil::PrintName(N);					\
    std::cout << std::boolalpha << (FLAG) << std::noboolalpha << std::endl; \
  } while(false)
#define PFV(N,VALUE)					\
  do {							\
    AliForwardUtil::PrintName(N);			\
    std::cout << (VALUE) << std::endl; } while(false)
  
//____________________________________________________________________
void AliPtEtaPidWeights::Print(Option_t* option) const
{
  PFV("MC Weights", "PtEta PID based");
#if 0
  gROOT->IncreaseDirLevel();
  for (Int_t i = 0; i < fPdgs.GetSize(); i++)  {
    Int_t v = fPdgs.At(i);
    if (v == 0) break;
    // PFV(Form("%d", v));
    // fWeights.Print("all"); 		
  }
  gROOT->DecreaseDirLevel();
#endif
}

//____________________________________________________________________
//
// EOF
// 
