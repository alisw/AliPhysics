/*
 * AliFemtoPPbpbLamZVtxMultContainer.cxx
 *
 *  Created on: Aug 30, 2017
 *      Author: gu74req
 */
//#include "AliLog.h"
#include <iostream>
#include "AliFemtoDreamZVtxMultContainer.h"
#include "TLorentzVector.h"
#include "TDatabasePDG.h"
ClassImp(AliFemtoDreamPartContainer)
AliFemtoDreamZVtxMultContainer::AliFemtoDreamZVtxMultContainer()
:fPartContainer(0),
 fPDGParticleSpecies(0)
{

}

AliFemtoDreamZVtxMultContainer::AliFemtoDreamZVtxMultContainer(
    AliFemtoDreamCollConfig *conf)
:fPartContainer(conf->GetNParticles(),
                AliFemtoDreamPartContainer(conf->GetMixingDepth()))
,fPDGParticleSpecies(conf->GetPDGCodes())
{
  std::cout<<"Number of Particles: "<<conf->GetNParticles()<<
      "MixingDepth: "<<conf->GetMixingDepth()<<std::endl;
}

AliFemtoDreamZVtxMultContainer::~AliFemtoDreamZVtxMultContainer() {
  // TODO Auto-generated destructor stub
}

void AliFemtoDreamZVtxMultContainer::SetEvent(
    std::vector<std::vector<AliFemtoDreamBasePart>> &Particles)
{
  //This method sets the particles of an event only in the case, that
  //more than one particle was identified, to avoid empty events.
  //  if (Particles->size()!=fParticleSpecies){
  //    TString errMessage = Form("Number of Input Particlese (%d) doese not"
  //        "correspond to the Number of particles Set (%d)",Particles->size(),
  //        fParticleSpecies);
  //    AliFatal(errMessage.Data());
  //  } else {
  std::vector<std::vector<AliFemtoDreamBasePart>>::iterator itInput=
      Particles.begin();
  std::vector<AliFemtoDreamPartContainer>::iterator itContainer=
      fPartContainer.begin();
  while (itContainer!=fPartContainer.end()) {
    if (itInput->size()>0) {
      itContainer->SetEvent(*itInput);
    }
    ++itInput;
    ++itContainer;
  }
  //  }
}
void AliFemtoDreamZVtxMultContainer::PairParticlesSE(
    std::vector<std::vector<AliFemtoDreamBasePart>> &Particles,
    AliFemtoDreamCorrHists *ResultsHist,int iMult)
{
  double RelativeK = 0;
  int _intSpec1 = 0;
  int HistCounter=0;
  //First loop over all the different Species
 auto itPDGPar1 = fPDGParticleSpecies.begin();
  for (auto itSpec1=Particles.begin();itSpec1!=Particles.end();++itSpec1) {
    auto itPDGPar2 = fPDGParticleSpecies.begin();
    itPDGPar2+=itSpec1-Particles.begin();
    int _intSpec2 = itSpec1-Particles.begin();
    for (auto itSpec2=itSpec1;itSpec2!=Particles.end();++itSpec2) {
      ResultsHist->FillPartnersSE(HistCounter,itSpec1->size(),itSpec2->size());
      //Now loop over the actual Particles and correlate them
      for (auto itPart1=itSpec1->begin();itPart1!=itSpec1->end();++itPart1) {
        std::vector<AliFemtoDreamBasePart>::iterator itPart2;
        if (itSpec1==itSpec2) {
          itPart2=itPart1+1;
        } else {
          itPart2=itSpec2->begin();
        }
        while (itPart2!=itSpec2->end()) {
          RelativeK=RelativePairMomentum(itPart1->GetMomentum(),*itPDGPar1,
                                         itPart2->GetMomentum(),
                                         *itPDGPar2);
          ResultsHist->FillSameEventDist(HistCounter,RelativeK);
          if (ResultsHist->GetDoMultBinning()) {
            ResultsHist->FillSameEventMultDist(HistCounter,iMult+1,RelativeK);
          }
          ++itPart2;
        }
      }
      ++HistCounter;
      _intSpec2++;
      itPDGPar2++;
    }
    _intSpec1++;
    itPDGPar1++;
  }
}

void AliFemtoDreamZVtxMultContainer::PairParticlesME(
    std::vector<std::vector<AliFemtoDreamBasePart>> &Particles,
    AliFemtoDreamCorrHists *ResultsHist,int iMult)
{
  double RelativeK = 0;
  int _intSpec1 = 0;
  int HistCounter=0;
  auto itPDGPar1 = fPDGParticleSpecies.begin();
  //First loop over all the different Species
  for (auto itSpec1=Particles.begin();itSpec1!=Particles.end();++itSpec1) {
    //We dont want to correlate the particles twice. Mixed Event Dist. of
    //Particle1 + Particle2 == Particle2 + Particle 1
    int SkipPart=itSpec1-Particles.begin();
    int _intSpec2 = SkipPart;
    auto itPDGPar2 = fPDGParticleSpecies.begin()+
        SkipPart;
    for (auto itSpec2=fPartContainer.begin()+SkipPart;
        itSpec2!=fPartContainer.end();++itSpec2) {
      for(int iDepth=0;iDepth<(int)itSpec2->GetMixingDepth();++iDepth){
        std::vector<AliFemtoDreamBasePart> ParticlesOfEvent=
            itSpec2->GetEvent(iDepth);
        ResultsHist->FillPartnersME(
            HistCounter,itSpec1->size(),ParticlesOfEvent.size());
        for (auto itPart1=itSpec1->begin();itPart1!=itSpec1->end();++itPart1) {
          for(auto itPart2=ParticlesOfEvent.begin();
              itPart2!=ParticlesOfEvent.end();++itPart2) {
            RelativeK=RelativePairMomentum(itPart1->GetMomentum(),*itPDGPar1,
                                           itPart2->GetMomentum(),*itPDGPar2);
            ResultsHist->FillMixedEventDist(HistCounter,RelativeK);
            if (ResultsHist->GetDoMultBinning()) {
              ResultsHist->FillMixedEventMultDist(HistCounter,iMult+1,RelativeK);
            }
          }
        }
      }
      ++HistCounter;
      ++_intSpec2;
      ++itPDGPar2;
    }
    ++_intSpec1;
    ++itPDGPar1;
  }
}
double AliFemtoDreamZVtxMultContainer::RelativePairMomentum(TVector3 Part1Momentum,
                                                            int PDGPart1,
                                                            TVector3 Part2Momentum,
                                                            int PDGPart2)
{
    if(PDGPart1 == 0 || PDGPart2== 0){
      AliError("Invalid PDG Code");
    }
  double results = 0.;
  TLorentzVector SPtrack,TPProng,trackSum,SPtrackCMS,TPProngCMS;
  //Even if the Daughter tracks were switched up during PID doesn't play a role here cause we are
  //only looking at the mother mass
  SPtrack.SetXYZM(Part1Momentum.X(), Part1Momentum.Y(),Part1Momentum.Z(),
                  TDatabasePDG::Instance()->GetParticle(PDGPart1)->Mass());
  TPProng.SetXYZM(Part2Momentum.X(), Part2Momentum.Y(),Part2Momentum.Z(),
                  TDatabasePDG::Instance()->GetParticle(PDGPart2)->Mass());
  trackSum = SPtrack + TPProng;

  double beta = trackSum.Beta();
  double betax = beta*cos(trackSum.Phi())*sin(trackSum.Theta());
  double betay = beta*sin(trackSum.Phi())*sin(trackSum.Theta());
  double betaz = beta*cos(trackSum.Theta());

  SPtrackCMS = SPtrack;
  TPProngCMS = TPProng;

  SPtrackCMS.Boost(-betax,-betay,-betaz);
  TPProngCMS.Boost(-betax,-betay,-betaz);

  TLorentzVector trackRelK;

  trackRelK = SPtrackCMS - TPProngCMS;
  results = 0.5*trackRelK.P();
  return results;
}
