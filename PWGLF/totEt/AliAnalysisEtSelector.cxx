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
,fClusterArray(0)
,fCuts(cuts)
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
Bool_t AliAnalysisEtSelector::FromSecondaryInteraction(const TParticle& part, AliStack &stack) const
{
//   Bool_t partVtxSecondary = (
// 			      TMath::Sqrt(part.Vx()*part.Vx() + part.Vy()*part.Vy()) > fCuts->GetPrimaryVertexCutXY() 
// 			      || TMath::Abs(part.Vz()) > fCuts->GetPrimaryVertexCutZ()
// 			    )
// 			    && TMath::Sqrt(part.Vx()*part.Vx()+part.Vy()*part.Vy() + part.Vz()*part.Vz())<(fCuts->GetGeometryPhosDetectorRadius()-10);
  
  //Let's find suspect decay (typical for secondary interaction)...
  
  return SuspeciousDecayInChain(211, 111, part, stack);
  
			    
  
  
}

Bool_t AliAnalysisEtSelector::SuspeciousDecayInChain(const UInt_t suspectMotherPdg, const UInt_t suspectDaughterPdg, const TParticle &part, AliStack& stack) const
{
  UInt_t partPdg = TMath::Abs(part.GetPdgCode());
  if(part.GetFirstMother() == -1)
  {
    return kFALSE;
  }
  TParticle *mother = stack.Particle(part.GetFirstMother()); 
  UInt_t motherPdg = TMath::Abs(mother->GetPdgCode());
  if((suspectDaughterPdg==partPdg || 2112 == partPdg )&& suspectMotherPdg == motherPdg)
  {
    return kTRUE;
  }
  return SuspeciousDecayInChain(suspectMotherPdg, suspectDaughterPdg, *mother, stack);
}
