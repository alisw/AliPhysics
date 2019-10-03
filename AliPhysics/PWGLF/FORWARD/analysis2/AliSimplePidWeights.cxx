/**
 * @file   AliSimplePidWeights.cxx
 * @author Marek Chojnacki <mchojnac#cern.ch>
 * @date   Wed Feb  4 00:09:32 2015
 * 
 * @brief  Implementation of simple weights
 * 
 * 
 */

#include "AliSimplePidWeights.h"
#include "AliForwardUtil.h"
#include <TROOT.h>
#include <iostream>

//____________________________________________________________________
AliSimplePidWeights::AliSimplePidWeights()
  : AliBaseMCWeights(),
    fPdgs(10),
    fWeights(10)
{
}
//____________________________________________________________________
AliSimplePidWeights::AliSimplePidWeights(const AliSimplePidWeights& o)
  : AliBaseMCWeights(o),
    fPdgs(o.fPdgs),
    fWeights(o.fWeights)
{
}

//____________________________________________________________________
AliSimplePidWeights&
AliSimplePidWeights::operator=(const AliSimplePidWeights& o)
{
  if (&o == this) return *this; 

  AliBaseMCWeights::operator=(o);
  fWeights = o.fWeights;
  fPdgs = o.fPdgs;
  
  return *this;
}
//____________________________________________________________________
void AliSimplePidWeights::Init(TList*)
{}

//____________________________________________________________________
Double_t AliSimplePidWeights::CalcWeight(Double_t,
					 Double_t,
					 Double_t,
					 Int_t id,
					 Double_t,
					 Double_t) const 
{
  for (Int_t i = 0; i < fPdgs.GetSize(); i++)  {
    Int_t v = fPdgs.At(i);
    if (v == 0) break;
    if (v == id) return fWeights.At(i);
  }
  return 1.;
}
//____________________________________________________________________
void AliSimplePidWeights::AddPDGCode(Int_t pdg, Double_t factor, Bool_t anti)
{

  Int_t n   = fPdgs.GetSize();
  Int_t pos = 0;

  // Find first empty slot 
  for (pos = 0; pos < n; pos++) if (fPdgs.At(pos) == 0) break;

  // If no empty slot is found, resize to twice the size 
  if (pos >= n) {
    fPdgs.Set(2*n);
    fWeights.Set(2*n);
  }

  // Now set the slot 
  fPdgs.AddAt(pdg, pos);
  fWeights.AddAt(factor, pos);

  if (anti) AddPDGCode(-pdg,factor,false);
}
#define PF(N,V,...)					\
  AliForwardUtil::PrintField(N,V, ## __VA_ARGS__)
#define PFB(N,FLAG)				\
  do {									\
    AliForwardUtil::PrintName(N);					\
    std::cout << std::boolalpha << (FLAG) << std::noboolalpha << std::endl; \
  } while(false)
#define PFV(N,VALUE)					\
  do {							\
    AliForwardUtil::PrintName(N);			\
    std::cout << (VALUE) << std::endl; } while(false)
  
//____________________________________________________________________
void AliSimplePidWeights::Print(Option_t* option) const
{
  PFV("MC Weights", "Simple PID based");
  gROOT->IncreaseDirLevel();
  for (Int_t i = 0; i < fPdgs.GetSize(); i++)  {
    Int_t v = fPdgs.At(i);
    if (v == 0) break;
    PFV(Form("%d", v), fWeights.At(i));
  }
  gROOT->DecreaseDirLevel();
}

//____________________________________________________________________
//
// EOF
// 
