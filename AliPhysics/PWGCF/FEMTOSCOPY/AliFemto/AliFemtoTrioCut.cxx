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
  fNpassed(0),
  fIncludeTrioMass(-1.0),
  fIncludeTrioDelta(9999)
{

}

AliFemtoTrioCut::~AliFemtoTrioCut()
{

}

bool AliFemtoTrioCut::Pass(AliFemtoTrio *trio)
{
  if(fIncludeTrioMass > 0){
    if(trio->MInv() < (fIncludeTrioMass-fIncludeTrioDelta) ||
       trio->MInv() > (fIncludeTrioMass+fIncludeTrioDelta)
       ){
      fNfailed++;
      return false;
    }
  }
  AliFemtoTrio::EPart trioType[3];
  
  trioType[0] = trio->GetTrack1type();
  trioType[1] = trio->GetTrack2type();
  trioType[2] = trio->GetTrack3type();
  
  AliFemtoTrio::EPart pairType[2];
  int trioPart1index=-1,trioPart2index=-1;
  
  double pairMinv=0;
  double mass;
  double delta;

  for(unsigned int i=0;i<fExcludedPairsMasses.size();i++){
    pairType[0] = fExcludedPairsType1[i];
    pairType[1] = fExcludedPairsType2[i];
    
    trioPart1index=-1;trioPart2index=-1;
    
    for(int p1=0;p1<3;p1++){
      if(pairType[0] == trioType[p1]) trioPart1index=p1;
      if(pairType[1] == trioType[p1]) trioPart2index=p1;
    }
    if(trioPart1index<0 || trioPart2index<0) continue;
    if(trioPart1index > trioPart2index){
      int tmp = trioPart2index;
      trioPart2index = trioPart1index;
      trioPart1index = tmp;
    }
    
    if(trioPart1index == 0 && trioPart2index == 1) pairMinv = trio->MInv12();
    if(trioPart1index == 0 && trioPart2index == 2) pairMinv = trio->MInv31();
    if(trioPart1index == 1 && trioPart2index == 2) pairMinv = trio->MInv23();
    
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

void AliFemtoTrioCut::SetIncludeTrioOnly(double mass, double delta)
{
  fIncludeTrioMass = mass;
  fIncludeTrioDelta = delta;
}
