///
/// \file    AliFemtoTrioCut.h
/// \author  Jeremi Niedziela


#ifndef AliFemtoTrioCut_H
#define AliFemtoTrioCut_H

#include "AliFemtoTrio.h"
#include "AliFemtoString.h"

#include <TList.h>

//
// AliFemtoTrioCut - cut on a set of three tracks
//
class AliFemtoTrioCut
{
public:
  AliFemtoTrioCut();
  ~AliFemtoTrioCut();

  bool Pass(AliFemtoTrio *trio);
  
  AliFemtoString Report(){return "";}
private:
  int fNfailed;
  int fNpassed;
  
#ifdef __ROOT__
  ClassDef(AliFemtoTrioCut, 0)
#endif
};

#endif
