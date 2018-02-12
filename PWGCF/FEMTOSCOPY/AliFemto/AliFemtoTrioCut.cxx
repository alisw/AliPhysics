///
/// \file    AliFemtoTrioCut.cxx
/// \author  Jeremi Niedziela

#include "AliFemtoTrioCut.h"

#include <TString.h>
#include <TObjString.h>

#ifdef __ROOT__
ClassImp(AliFemtoTrioCut)
#endif

AliFemtoTrioCut::AliFemtoTrioCut():
  fNfailed(0),
  fNpassed(0)
{

}

AliFemtoTrioCut::~AliFemtoTrioCut()
{

}

bool AliFemtoTrioCut::Pass(AliFemtoTrio *trio)
{
  AliFemtoTrio::EPart trioType[3];
  
  trioType[0] = trio->GetTrack1type();
  trioType[1] = trio->GetTrack2type();
  trioType[2] = trio->GetTrack3type();
  
  AliFemtoTrio::EPart pairType[2];
  int trioPart1index=-1,trioPart2index=-1;
  
  AliFemtoParticle *track1=nullptr;
  AliFemtoParticle *track2=nullptr;
  
  double pairMinv, mass, delta;
  
  for(int i=0;i<fExcludedPairsMasses.size();i++){
    pairType[0] = fExcludedPairsType1[i];
    pairType[1] = fExcludedPairsType2[i];
    
    for(int p1=0;p1<3;p1++){
      if(pairType[0] == trioType[p1]) trioPart1index=p1;
      if(pairType[1] == trioType[p1]) trioPart2index=p1;
    }
    if(trioPart1index<0 || trioPart2index<0) continue;
    
    if(trioPart1index == 0) track1 = trio->Track1();
    if(trioPart1index == 1) track1 = trio->Track2();
    if(trioPart1index == 2) track1 = trio->Track3();
    
    if(trioPart2index == 0) track2 = trio->Track1();
    if(trioPart2index == 1) track2 = trio->Track2();
    if(trioPart2index == 2) track2 = trio->Track3();
    
    if(track1 == track2){
      cout<<"AliFemtoTrioCut::Pass() -- Warning -- same track counted twice!"<<endl;
      continue;
    }
    
    pairMinv = GetPairMInv(track1, track2);
    mass = fExcludedPairsMasses[i];
    delta = fExcludedPairsDeltas[i];
    
    if(pairMinv > (mass-delta) && pairMinv < (mass+delta)){
      fNfailed++;
      return false;
    }
  }
  
  fNpassed++;
  return true;
}

void AliFemtoTrioCut::SetExcludePair(double mass, double delta,
                                     AliFemtoTrio::EPart type1, AliFemtoTrio::EPart type2)
{
  fExcludedPairsMasses.push_back(mass);
  fExcludedPairsDeltas.push_back(delta);
  fExcludedPairsType1.push_back(type1);
  fExcludedPairsType2.push_back(type2);
}

double AliFemtoTrioCut::GetPairMInv(AliFemtoParticle *track1,AliFemtoParticle *track2)
{
  if(!track1 || !track2){
    cout<<"W - AliFemtoTrioCut::GetPairMInv - track missing:"<<endl;
    cout<<"track1:"<<track1<<endl;
    cout<<"track2:"<<track2<<endl;
    return -1.0;
  }
  AliFemtoLorentzVector p1 = track1->FourMomentum();
  AliFemtoLorentzVector p2 = track2->FourMomentum();
  
  double E  = p1.e()  + p2.e();
  double px = p1.px() + p2.px();
  double py = p1.py() + p2.py();
  double pz = p1.pz() + p2.pz();
  
  return sqrt(E*E-(px*px+py*py+pz*pz));
}
