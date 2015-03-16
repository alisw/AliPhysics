/**
 * @file   AliBaseMCWeights.cxx
 * @author Marek Chojnacki <mchojnac#cern.ch>
 * @date   Wed Feb  4 00:30:56 2015
 * 
 * @brief  
 * 
 * 
 */
#include "AliBaseMCWeights.h"
#include "AliMCParticle.h"
#include "AliForwardUtil.h"
#include <iostream>
//____________________________________________________________________
AliBaseMCWeights&
AliBaseMCWeights::operator=(const AliBaseMCWeights& o)
{
  if (&o == this) return *this; 

  TObject::operator=(o);

  return *this;
}
//____________________________________________________________________
void
AliBaseMCWeights::Init(TList*)
{}
#define PFV(N,VALUE)					\
  do {							\
    AliForwardUtil::PrintName(N);			\
    std::cout << (VALUE) << std::endl; } while(false)

//____________________________________________________________________
void
AliBaseMCWeights::Print(Option_t* option) const
{
  PFV("MC weights", "Not specified");
}

//____________________________________________________________________
Double_t
AliBaseMCWeights::CalcWeight(const AliMCParticle* p,
			     Bool_t,
			     Double_t phiR,
			     Double_t b) const
{
  if (!p) return 1;
  return CalcWeight(p->Eta(), p->Pt(), p->Phi(), p->PdgCode(), phiR, b);
}
//____________________________________________________________________
//
// EOF
//

