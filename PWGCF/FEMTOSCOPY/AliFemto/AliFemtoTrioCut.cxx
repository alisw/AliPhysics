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
  // add here checks for all cuts
  
  fNpassed++;
  return true;
}
