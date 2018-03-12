///
/// \file    AliFemtoTrioCut.h
/// \author  Jeremi Niedziela


#ifndef AliFemtoTrioCut_H
#define AliFemtoTrioCut_H

#include "AliFemtoTrio.h"
#include "AliFemtoString.h"

#include <TList.h>
#include <vector>

//
// AliFemtoTrioCut - cut on a set of three tracks
//
class AliFemtoTrioCut
{
public:
  AliFemtoTrioCut();
  ~AliFemtoTrioCut();

  bool Pass(AliFemtoTrio *trio);
  
  void SetExcludePair(double mass, double delta, AliFemtoTrio::EPart type1, AliFemtoTrio::EPart type2);
  
  void SetIncludeTrioOnly(double mass, double delta);
  
  AliFemtoString Report(){return "";}
private:
  int fNfailed;
  int fNpassed;
  
  std::vector<double> fExcludedPairsMasses; // pair masses that should be excluded
  std::vector<double> fExcludedPairsDeltas; // delta for pair masses to exclude
  std::vector<AliFemtoTrio::EPart> fExcludedPairsType1; // first paticle from pair to exclude
  std::vector<AliFemtoTrio::EPart> fExcludedPairsType2; // second particle from pair to exlude
  
  double fIncludeTrioMass;
  double fIncludeTrioDelta;
  
#ifdef __ROOT__
  ClassDef(AliFemtoTrioCut, 0)
#endif
};

#endif
