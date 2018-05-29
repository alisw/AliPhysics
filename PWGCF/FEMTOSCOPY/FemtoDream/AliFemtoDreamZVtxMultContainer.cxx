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
{}

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
    AliFemtoDreamCorrHists *ResultsHist,int iMult,float cent)
{
  float RelativeK = 0;
  int HistCounter=0;
  //First loop over all the different Species
  auto itPDGPar1 = fPDGParticleSpecies.begin();
  for (auto itSpec1=Particles.begin();itSpec1!=Particles.end();++itSpec1) {
    auto itPDGPar2 = fPDGParticleSpecies.begin();
    itPDGPar2+=itSpec1-Particles.begin();
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
          if (ResultsHist->GetDokTBinning()) {
            ResultsHist->FillSameEventkTDist(
                HistCounter,
                RelativePairkT(
                    itPart1->GetMomentum(),*itPDGPar1,
                    itPart2->GetMomentum(),*itPDGPar2),
                    RelativeK,cent);
          }
          if (ResultsHist->GetDomTBinning()) {
            ResultsHist->FillSameEventmTDist(
                HistCounter,
                RelativePairmT(
                    itPart1->GetMomentum(),*itPDGPar1,
                    itPart2->GetMomentum(),*itPDGPar2),
                    RelativeK);
          }
          if (ResultsHist->GetEtaPhiPlots()) {
            DeltaEtaDeltaPhi(HistCounter,&(*itPart1),&(*itPart2),true,ResultsHist);
          }
          if (ResultsHist->GetDodPhidEtaPlots()) {
            float deta=itPart1->GetEta().at(0)-itPart2->GetEta().at(0);
            float dphi=itPart1->GetPhi().at(0)-itPart2->GetPhi().at(0);
            if (dphi < 0) {
              ResultsHist->FilldPhidEtaSE(HistCounter,dphi+2*TMath::Pi(),deta);
            } else {

              ResultsHist->FilldPhidEtaSE(HistCounter,dphi,deta);
            }
          }
          ++itPart2;
        }
      }
      ++HistCounter;
      itPDGPar2++;
    }
    itPDGPar1++;
  }
}
void AliFemtoDreamZVtxMultContainer::PairMCParticlesSE(
    std::vector<std::vector<AliFemtoDreamBasePart>> &Particles,
    AliFemtoDreamCorrHists *ResultsHist,int iMult)
{
  float RelativeK = 0;
  int HistCounter=0;
  //First loop over all the different Species
  auto itPDGPar1 = fPDGParticleSpecies.begin();
  for (auto itSpec1=Particles.begin();itSpec1!=Particles.end();++itSpec1) {
    auto itPDGPar2 = fPDGParticleSpecies.begin();
    itPDGPar2+=itSpec1-Particles.begin();
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
                                         itPart2->GetMomentum(),*itPDGPar2);
          //If the ancestor is the same fill one hist, if it isnt the other
          //          std::cout << itPart1->GetMotherID() << '\t' << itPart2->GetMotherID() << std::endl;
          if (itPart1->GetMotherID()==itPart2->GetMotherID()) {//common ancestor
            ResultsHist->FillSameEventCommonAncestDist(HistCounter,RelativeK);
          } else {//different ancestor
            ResultsHist->FillSameEventNonCommonAncestDist(
                HistCounter,RelativeK);
          }
          ++itPart2;
        }
      }
      ++HistCounter;
      itPDGPar2++;
    }
    itPDGPar1++;
  }
}

void AliFemtoDreamZVtxMultContainer::PairParticlesME(
    std::vector<std::vector<AliFemtoDreamBasePart>> &Particles,
    AliFemtoDreamCorrHists *ResultsHist,int iMult,float cent)
{
  float RelativeK = 0;
  int HistCounter=0;
  auto itPDGPar1 = fPDGParticleSpecies.begin();
  //First loop over all the different Species
  for (auto itSpec1=Particles.begin();itSpec1!=Particles.end();++itSpec1) {
    //We dont want to correlate the particles twice. Mixed Event Dist. of
    //Particle1 + Particle2 == Particle2 + Particle 1
    int SkipPart=itSpec1-Particles.begin();
    auto itPDGPar2 = fPDGParticleSpecies.begin()+SkipPart;
    for (auto itSpec2=fPartContainer.begin()+SkipPart;
        itSpec2!=fPartContainer.end();++itSpec2) {
      if(itSpec1->size()>0) {
        ResultsHist->FillEffectiveMixingDepth(
            HistCounter,(int)itSpec2->GetMixingDepth());
      }
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
            if (ResultsHist->GetDokTBinning()) {
              ResultsHist->FillMixedEventkTDist(
                  HistCounter,
                  RelativePairkT(
                      itPart1->GetMomentum(),*itPDGPar1,
                      itPart2->GetMomentum(),*itPDGPar2),
                      RelativeK,cent);
            }
            if (ResultsHist->GetDomTBinning()) {
              ResultsHist->FillMixedEventmTDist(
                  HistCounter,
                  RelativePairmT(
                      itPart1->GetMomentum(),*itPDGPar1,
                      itPart2->GetMomentum(),*itPDGPar2),
                      RelativeK);
            }
            if (ResultsHist->GetObtainMomentumResolution()) {
              //It is sufficient to do this in Mixed events, which allows
              //to increase the statistics. The Resolution of the tracks and therefore
              //of the pairs does not change event by event.
              //Now we only want to use the momentum of particles we are after, hence
              //we check the PDG Code!
              if ((*itPDGPar1==TMath::Abs(itPart1->GetMCPDGCode()))&&
                  ((*itPDGPar2==TMath::Abs(itPart2->GetMCPDGCode())))) {
                float RelKTrue=
                    RelativePairMomentum(itPart1->GetMCMomentum(),*itPDGPar1,
                                         itPart2->GetMCMomentum(),*itPDGPar2);
                ResultsHist->FillMomentumResolution(
                    HistCounter,RelKTrue,RelativeK);
              }
            }
            if (ResultsHist->GetEtaPhiPlots()) {
              DeltaEtaDeltaPhi(HistCounter,&(*itPart1),&(*itPart2),false,ResultsHist);
            }
            if (ResultsHist->GetDodPhidEtaPlots()) {
              float deta=itPart1->GetEta().at(0)-itPart2->GetEta().at(0);
              float dphi=itPart1->GetPhi().at(0)-itPart2->GetPhi().at(0);
              if (dphi < 0) {
                ResultsHist->FilldPhidEtaME(HistCounter,dphi+2*TMath::Pi(),deta);
              } else {
                ResultsHist->FilldPhidEtaME(HistCounter,dphi,deta);
              }
            }
          }
        }
      }
      ++HistCounter;
      ++itPDGPar2;
    }
    ++itPDGPar1;
  }
}
float AliFemtoDreamZVtxMultContainer::RelativePairMomentum(TVector3 Part1Momentum,
                                                           int PDGPart1,
                                                           TVector3 Part2Momentum,
                                                           int PDGPart2)
{
  if(PDGPart1 == 0 || PDGPart2== 0){
    AliError("Invalid PDG Code");
  }
  float results = 0.;
  TLorentzVector SPtrack,TPProng,trackSum,SPtrackCMS,TPProngCMS;
  //Even if the Daughter tracks were switched up during PID doesn't play a role here cause we are
  //only looking at the mother mass
  SPtrack.SetXYZM(Part1Momentum.X(), Part1Momentum.Y(),Part1Momentum.Z(),
                  TDatabasePDG::Instance()->GetParticle(PDGPart1)->Mass());
  TPProng.SetXYZM(Part2Momentum.X(), Part2Momentum.Y(),Part2Momentum.Z(),
                  TDatabasePDG::Instance()->GetParticle(PDGPart2)->Mass());
  trackSum = SPtrack + TPProng;

  float beta = trackSum.Beta();
  float betax = beta*cos(trackSum.Phi())*sin(trackSum.Theta());
  float betay = beta*sin(trackSum.Phi())*sin(trackSum.Theta());
  float betaz = beta*cos(trackSum.Theta());

  SPtrackCMS = SPtrack;
  TPProngCMS = TPProng;

  SPtrackCMS.Boost(-betax,-betay,-betaz);
  TPProngCMS.Boost(-betax,-betay,-betaz);

  TLorentzVector trackRelK;

  trackRelK = SPtrackCMS - TPProngCMS;
  results = 0.5*trackRelK.P();
  return results;
}
float AliFemtoDreamZVtxMultContainer::RelativePairkT(TVector3 Part1Momentum,
                                                     int PDGPart1,
                                                     TVector3 Part2Momentum,
                                                     int PDGPart2) {
  if(PDGPart1 == 0 || PDGPart2== 0){
    AliError("Invalid PDG Code");
  }
  float results = 0.;
  TLorentzVector SPtrack,TPProng,trackSum,SPtrackCMS,TPProngCMS;
  //Even if the Daughter tracks were switched up during PID doesn't play a role here cause we are
  //only looking at the mother mass
  SPtrack.SetXYZM(Part1Momentum.X(), Part1Momentum.Y(),Part1Momentum.Z(),
                  TDatabasePDG::Instance()->GetParticle(PDGPart1)->Mass());
  TPProng.SetXYZM(Part2Momentum.X(), Part2Momentum.Y(),Part2Momentum.Z(),
                  TDatabasePDG::Instance()->GetParticle(PDGPart2)->Mass());
  trackSum = SPtrack + TPProng;

  results=0.5*trackSum.Pt();
  return results;
}
float AliFemtoDreamZVtxMultContainer::RelativePairmT(TVector3 Part1Momentum,
                                                     int PDGPart1,
                                                     TVector3 Part2Momentum,
                                                     int PDGPart2) {
  if(PDGPart1 == 0 || PDGPart2== 0){
    AliError("Invalid PDG Code");
  }
  float results = 0.;
  TLorentzVector SPtrack,TPProng,trackSum,SPtrackCMS,TPProngCMS;
  //Even if the Daughter tracks were switched up during PID doesn't play a role here cause we are
  //only looking at the mother mass
  SPtrack.SetXYZM(Part1Momentum.X(), Part1Momentum.Y(),Part1Momentum.Z(),
                  TDatabasePDG::Instance()->GetParticle(PDGPart1)->Mass());
  TPProng.SetXYZM(Part2Momentum.X(), Part2Momentum.Y(),Part2Momentum.Z(),
                  TDatabasePDG::Instance()->GetParticle(PDGPart2)->Mass());
  trackSum = SPtrack + TPProng;
  float pairKT=0.5*trackSum.Pt();
  float averageMass = 0.5*(TDatabasePDG::Instance()->GetParticle(PDGPart1)->Mass()
      + TDatabasePDG::Instance()->GetParticle(PDGPart2)->Mass());
  results= TMath::Sqrt(pow(pairKT,2.) + pow(averageMass,2.));
  return results;
}

void AliFemtoDreamZVtxMultContainer::DeltaEtaDeltaPhi(
    int Hist,AliFemtoDreamBasePart *part1, AliFemtoDreamBasePart *part2,
    bool SEorME, AliFemtoDreamCorrHists *ResultsHist) {
  //used to check for track splitting/merging
  //this function only produces meaningful results for track with x Daughter
  //looking at this quantity makes only sense anyways for Track - Track not
  //for v0 - v0 ...
  float eta1=part1->GetEta().at(0);
  std::vector<float> Phirad1=part1->GetPhiAtRaidius().at(0);

  std::vector<float> eta2=part2->GetEta();
  for (unsigned int iDaug=0;iDaug<part2->GetPhiAtRaidius().size();++iDaug) {
    std::vector<float> phiAtRad2=part2->GetPhiAtRaidius().at(iDaug);
    float etaPar2;
    if (part2->GetPhiAtRaidius().size()==1) {
      etaPar2=eta2.at(0);
    } else {
      etaPar2=eta2.at(iDaug+1);
    }
    float deta=TMath::Abs(eta1-etaPar2);
    for (int iRad=0;iRad<9;++iRad) {
      float dphi=TMath::Abs(Phirad1.at(iRad)-phiAtRad2.at(iRad));
      if (SEorME) {
        ResultsHist->FillEtaPhiAtRadiiSE(Hist,iDaug,iRad,dphi,deta);
      } else {
        ResultsHist->FillEtaPhiAtRadiiME(Hist,iDaug,iRad,dphi,deta);
      }
    }
  }
}
