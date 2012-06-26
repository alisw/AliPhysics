//_________________________________________________________________________
//  Utility Class for transverse energy studies
//  Selector Base class 
//  -  
// implementation file
//
//*-- Authors: Oystein Djuvsland (Bergen)
//_________________________________________________________________________
#include "AliAnalysisEtSelector.h"
#include "AliAnalysisEtCuts.h"
#include "AliESDCaloCluster.h"
#include "AliStack.h"
#include "TParticle.h"
#include "TParticlePDG.h"
#include "AliAnalysisEtCommon.h"
#include <iostream>

ClassImp(AliAnalysisEtSelector);

AliAnalysisEtSelector::AliAnalysisEtSelector(AliAnalysisEtCuts *cuts) : AliAnalysisEtCommon()
,fEvent(0)
,fCuts(cuts)
,fClusterArray(0)
,fRunNumber(0)
{
}

AliAnalysisEtSelector::~AliAnalysisEtSelector()
{ // dtor
  if(fClusterArray)
  {
    delete fClusterArray;
  }
}

Bool_t AliAnalysisEtSelector::CutNeutralMcParticle(Int_t pIdx, AliStack& s, const TParticlePDG& pdg) const
{
  return s.IsPhysicalPrimary(pIdx) &&(TMath::Abs(TMath::Abs(pdg.Charge()) - fCuts->GetMonteCarloNeutralParticle())<1e-3);
}

Bool_t AliAnalysisEtSelector::IsEmEtParticle(const Int_t pdgCode) const
{
  return pdgCode == fgGammaCode || pdgCode == fgPi0Code || pdgCode == fgEtaCode || pdgCode == fgEPlusCode || pdgCode == fgEMinusCode;
}


Bool_t AliAnalysisEtSelector::PrimaryIsEmEtParticle(const Int_t pIdx, AliStack& stack) const
{
  return IsEmEtParticle(stack.Particle(GetPrimary(pIdx, stack))->GetPdgCode());
}
Int_t AliAnalysisEtSelector::GetPrimary(const Int_t partIdx, AliStack& stack) const
{ // get primary
  if(partIdx >= 0) 
  {
    Int_t mothIdx = stack.Particle(partIdx)->GetMother(0);
    if(mothIdx < 0) return -1;
    TParticle *mother = stack.Particle(mothIdx);
    if(mother)
    {
      if(stack.IsPhysicalPrimary(mothIdx)) return mothIdx;
      else return GetPrimary(mothIdx, stack);
    }
    else 
    {
      std::cout << "WAT!" << std::endl;
      return -1;
    }
  }
  return -1;
}
